from __future__ import annotations

import csv
import math
import random
from collections import Counter
from copy import deepcopy
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import yaml
from Bio import Phylo
from Bio.Phylo.BaseTree import Clade, Tree

from .analysis import analyze_locus


DEFAULT_SCENARIOS = [
    "concordant_singlecopy",
    "branch_length_only",
    "singlecopy_discordance",
    "multicopy_duplication",
    "pseudoorthology",
    "missing_taxa",
    "gene_tree_error",
]


@dataclass(frozen=True, slots=True)
class ScenarioTruth:
    name: str
    description: str
    truth_class: str
    expected_statuses: tuple[str, ...]
    truth_has_duplication: bool
    truth_single_copy: bool


SCENARIO_TRUTH: dict[str, ScenarioTruth] = {
    "concordant_singlecopy": ScenarioTruth(
        name="concordant_singlecopy",
        description="Single-copy gene tree with the same topology as the species tree.",
        truth_class="clean_singlecopy",
        expected_statuses=("PASS",),
        truth_has_duplication=False,
        truth_single_copy=True,
    ),
    "branch_length_only": ScenarioTruth(
        name="branch_length_only",
        description=(
            "Single-copy gene tree with species-tree topology but independently perturbed branch lengths."
        ),
        truth_class="clean_singlecopy",
        expected_statuses=("PASS",),
        truth_has_duplication=False,
        truth_single_copy=True,
    ),
    "singlecopy_discordance": ScenarioTruth(
        name="singlecopy_discordance",
        description=(
            "Single-copy topology conflict with no simulated duplication; a structural stress case for discordance."
        ),
        truth_class="singlecopy_discordance",
        expected_statuses=("REVIEW",),
        truth_has_duplication=False,
        truth_single_copy=True,
    ),
    "multicopy_duplication": ScenarioTruth(
        name="multicopy_duplication",
        description=(
            "One sampled species contains two descendant copies created by an explicit terminal duplication."
        ),
        truth_class="multicopy_duplication",
        expected_statuses=("RESOLVE",),
        truth_has_duplication=True,
        truth_single_copy=False,
    ),
    "pseudoorthology": ScenarioTruth(
        name="pseudoorthology",
        description=(
            "One sequence per species sampled from two deep hidden lineages, mimicking pseudo-orthology risk."
        ),
        truth_class="hidden_paralog_structure",
        expected_statuses=("REVIEW",),
        truth_has_duplication=True,
        truth_single_copy=True,
    ),
    "missing_taxa": ScenarioTruth(
        name="missing_taxa",
        description=(
            "Single-copy concordant topology after pruning enough taxa to fall below the coverage threshold."
        ),
        truth_class="low_coverage",
        expected_statuses=("REVIEW",),
        truth_has_duplication=False,
        truth_single_copy=True,
    ),
    "gene_tree_error": ScenarioTruth(
        name="gene_tree_error",
        description=(
            "Single-copy tree with seeded label permutations representing gene-tree estimation/topology error."
        ),
        truth_class="singlecopy_tree_error",
        expected_statuses=("REVIEW",),
        truth_has_duplication=False,
        truth_single_copy=True,
    ),
}


@dataclass(slots=True)
class BenchmarkRecord:
    benchmark_id: str
    scenario: str
    replicate: int
    seed: int
    taxa: int
    truth_class: str
    truth_has_duplication: bool
    truth_single_copy: bool
    expected_statuses: str
    observed_status: str
    acceptable: bool
    false_resolve: bool
    primary_concern: str
    coverage: float
    max_copies: int
    overlap_duplications: int
    lca_duplications: int
    inferred_losses: int
    normalized_rf: float | None
    species_tree: str
    gene_tree: str
    mapping: str = ""
    delimiter: str = "|"
    species_field: int = 0


@dataclass(slots=True)
class BenchmarkConfig:
    seed: int = 20260926
    replicates: int = 5
    taxa: tuple[int, ...] = (8, 16)
    scenarios: tuple[str, ...] = tuple(DEFAULT_SCENARIOS)
    min_coverage: float = 0.95
    discordance_review_threshold: float | None = None
    missing_fraction: float = 0.25
    branch_length_sigma: float = 0.35
    gene_tree_error_swaps: int = 2


def load_benchmark_config(path: str | Path | None = None) -> BenchmarkConfig:
    if path is None:
        return BenchmarkConfig()
    raw = yaml.safe_load(Path(path).read_text(encoding="utf-8")) or {}
    config = BenchmarkConfig(
        seed=int(raw.get("seed", 20260926)),
        replicates=int(raw.get("replicates", 5)),
        taxa=tuple(int(value) for value in raw.get("taxa", [8, 16])),
        scenarios=tuple(raw.get("scenarios", DEFAULT_SCENARIOS)),
        min_coverage=float(raw.get("min_coverage", 0.95)),
        discordance_review_threshold=raw.get("discordance_review_threshold"),
        missing_fraction=float(raw.get("missing_fraction", 0.25)),
        branch_length_sigma=float(raw.get("branch_length_sigma", 0.35)),
        gene_tree_error_swaps=int(raw.get("gene_tree_error_swaps", 2)),
    )
    _validate_config(config)
    return config


def _validate_config(config: BenchmarkConfig) -> None:
    if config.replicates < 1:
        raise ValueError("replicates must be >= 1")
    if not config.taxa or any(value < 4 for value in config.taxa):
        raise ValueError("taxa must contain integers >= 4")
    unknown = sorted(set(config.scenarios) - set(SCENARIO_TRUTH))
    if unknown:
        raise ValueError("Unknown benchmark scenarios: " + ", ".join(unknown))
    if not 0 < config.min_coverage <= 1:
        raise ValueError("min_coverage must be in (0, 1]")
    if not 0 < config.missing_fraction < 1:
        raise ValueError("missing_fraction must be in (0, 1)")
    if config.branch_length_sigma < 0:
        raise ValueError("branch_length_sigma must be >= 0")
    if config.gene_tree_error_swaps < 1:
        raise ValueError("gene_tree_error_swaps must be >= 1")
    if config.discordance_review_threshold is not None:
        config.discordance_review_threshold = float(config.discordance_review_threshold)
        if not 0 <= config.discordance_review_threshold <= 1:
            raise ValueError("discordance_review_threshold must be between 0 and 1")


def _random_binary_tree(taxa: list[str], rng: random.Random) -> Tree:
    nodes = [Clade(name=name, branch_length=rng.uniform(0.05, 0.25)) for name in taxa]
    rng.shuffle(nodes)
    while len(nodes) > 1:
        left = nodes.pop(rng.randrange(len(nodes)))
        right = nodes.pop(rng.randrange(len(nodes)))
        parent = Clade(branch_length=rng.uniform(0.05, 0.25), clades=[left, right])
        nodes.append(parent)
    return Tree(root=nodes[0], rooted=True)


def _rename_gene_tips(tree: Tree) -> None:
    for terminal in tree.get_terminals():
        terminal.name = f"{terminal.name}|g1"


def _cross_root_pair(tree: Tree, rng: random.Random) -> tuple[Clade, Clade]:
    if len(tree.root.clades) < 2:
        terminals = tree.get_terminals()
        return terminals[0], terminals[-1]
    left = tree.root.clades[0].get_terminals()
    right = tree.root.clades[1].get_terminals()
    return rng.choice(left), rng.choice(right)


def _swap_cross_root_labels(tree: Tree, rng: random.Random, swaps: int = 1) -> None:
    for _ in range(swaps):
        left, right = _cross_root_pair(tree, rng)
        left.name, right.name = right.name, left.name


def _jitter_branch_lengths(tree: Tree, rng: random.Random, sigma: float) -> None:
    for clade in tree.find_clades(order="preorder"):
        if clade.branch_length is None:
            continue
        multiplier = math.exp(rng.gauss(0.0, sigma))
        clade.branch_length = max(1e-6, clade.branch_length * multiplier)


def _find_parent(root: Clade, target: Clade) -> Clade | None:
    for child in root.clades:
        if child is target:
            return root
        found = _find_parent(child, target)
        if found is not None:
            return found
    return None


def _add_terminal_duplication(tree: Tree, rng: random.Random) -> None:
    target = rng.choice(tree.get_terminals())
    species = str(target.name)
    parent = _find_parent(tree.root, target)
    duplicate = Clade(
        branch_length=target.branch_length,
        clades=[
            Clade(name=f"{species}|copy1", branch_length=0.03),
            Clade(name=f"{species}|copy2", branch_length=0.03),
        ],
    )
    if parent is None:
        tree.root = duplicate
    else:
        index = parent.clades.index(target)
        parent.clades[index] = duplicate
    for terminal in tree.get_terminals():
        if "|" not in str(terminal.name):
            terminal.name = f"{terminal.name}|g1"


def _make_pseudoorthology_tree(taxa: list[str], rng: random.Random) -> Tree:
    shuffled = taxa[:]
    rng.shuffle(shuffled)
    midpoint = max(2, len(shuffled) // 2)
    group_a = shuffled[:midpoint]
    group_b = shuffled[midpoint:]
    if len(group_b) < 2:
        group_b = group_a[-2:]
        group_a = group_a[:-2]
    tree_a = _random_binary_tree(group_a, rng).root
    tree_b = _random_binary_tree(group_b, rng).root
    root = Clade(branch_length=0.02, clades=[tree_a, tree_b])
    tree = Tree(root=root, rooted=True)
    _rename_gene_tips(tree)
    return tree


def _prune_missing_taxa(tree: Tree, rng: random.Random, fraction: float) -> None:
    terminals = tree.get_terminals()
    remove_n = max(1, math.ceil(len(terminals) * fraction))
    remove_n = min(remove_n, len(terminals) - 4)
    for terminal in rng.sample(terminals, remove_n):
        tree.prune(terminal)


def _scenario_tree(
    species_tree: Tree,
    scenario: str,
    rng: random.Random,
    *,
    missing_fraction: float,
    branch_length_sigma: float,
    gene_tree_error_swaps: int,
) -> Tree:
    taxa = [str(node.name) for node in species_tree.get_terminals()]
    if scenario == "pseudoorthology":
        return _make_pseudoorthology_tree(taxa, rng)

    gene_tree = deepcopy(species_tree)
    if scenario == "branch_length_only":
        _jitter_branch_lengths(gene_tree, rng, branch_length_sigma)
    elif scenario == "singlecopy_discordance":
        _swap_cross_root_labels(gene_tree, rng, swaps=1)
    elif scenario == "multicopy_duplication":
        _add_terminal_duplication(gene_tree, rng)
        return gene_tree
    elif scenario == "missing_taxa":
        _prune_missing_taxa(gene_tree, rng, missing_fraction)
    elif scenario == "gene_tree_error":
        _swap_cross_root_labels(gene_tree, rng, swaps=gene_tree_error_swaps)
    elif scenario != "concordant_singlecopy":
        raise ValueError(f"Unsupported scenario: {scenario}")
    _rename_gene_tips(gene_tree)
    return gene_tree


def _write_tree(tree: Tree, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(path), "newick")


def _scenario_seed(base_seed: int, taxon_count: int, scenario_index: int, replicate: int) -> int:
    return base_seed + taxon_count * 100_000 + scenario_index * 10_000 + replicate


def run_native_benchmark(
    config: BenchmarkConfig,
    outdir: str | Path,
) -> list[BenchmarkRecord]:
    outdir = Path(outdir)
    simulation_dir = outdir / "simulations"
    simulation_dir.mkdir(parents=True, exist_ok=True)
    records: list[BenchmarkRecord] = []

    for taxon_count in config.taxa:
        species_names = [f"S{index:03d}" for index in range(1, taxon_count + 1)]
        for scenario_index, scenario in enumerate(config.scenarios):
            truth = SCENARIO_TRUTH[scenario]
            for replicate in range(1, config.replicates + 1):
                seed = _scenario_seed(config.seed, taxon_count, scenario_index, replicate)
                rng = random.Random(seed)
                species_tree = _random_binary_tree(species_names, rng)
                gene_tree = _scenario_tree(
                    species_tree,
                    scenario,
                    rng,
                    missing_fraction=config.missing_fraction,
                    branch_length_sigma=config.branch_length_sigma,
                    gene_tree_error_swaps=config.gene_tree_error_swaps,
                )
                benchmark_id = f"{scenario}_n{taxon_count}_r{replicate:04d}"
                replicate_dir = simulation_dir / scenario / f"n{taxon_count}" / f"r{replicate:04d}"
                species_path = replicate_dir / "species_tree.nwk"
                gene_path = replicate_dir / "gene_tree.nwk"
                _write_tree(species_tree, species_path)
                _write_tree(gene_tree, gene_path)

                result = analyze_locus(
                    species_path,
                    gene_path,
                    family=benchmark_id,
                    delimiter="|",
                    species_field=0,
                    min_coverage=config.min_coverage,
                    discordance_review_threshold=config.discordance_review_threshold,
                )
                acceptable = result.verdict.status in truth.expected_statuses
                false_resolve = (
                    result.verdict.status == "RESOLVE" and not truth.truth_has_duplication
                )
                records.append(
                    BenchmarkRecord(
                        benchmark_id=benchmark_id,
                        scenario=scenario,
                        replicate=replicate,
                        seed=seed,
                        taxa=taxon_count,
                        truth_class=truth.truth_class,
                        truth_has_duplication=truth.truth_has_duplication,
                        truth_single_copy=truth.truth_single_copy,
                        expected_statuses="|".join(truth.expected_statuses),
                        observed_status=result.verdict.status,
                        acceptable=acceptable,
                        false_resolve=false_resolve,
                        primary_concern=result.verdict.primary_concern,
                        coverage=result.coverage,
                        max_copies=result.max_copies,
                        overlap_duplications=result.species_overlap_duplications,
                        lca_duplications=result.reconciliation.lca_duplications,
                        inferred_losses=result.reconciliation.inferred_losses,
                        normalized_rf=result.topology.normalized_rf,
                        species_tree=str(species_path.relative_to(outdir)),
                        gene_tree=str(gene_path.relative_to(outdir)),
                        mapping="",
                        delimiter="|",
                        species_field=0,
                    )
                )
    write_benchmark_outputs(records, outdir, config=config, source="native")
    return records


def _read_bool(value: str) -> bool:
    return value.strip().lower() in {"1", "true", "yes", "y"}


def run_manifest_benchmark(
    manifest_path: str | Path,
    outdir: str | Path,
    *,
    min_coverage: float = 0.95,
    discordance_review_threshold: float | None = None,
) -> list[BenchmarkRecord]:
    manifest_path = Path(manifest_path)
    outdir = Path(outdir)
    records: list[BenchmarkRecord] = []
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "benchmark_id",
            "scenario",
            "species_tree",
            "gene_tree",
            "expected_statuses",
            "truth_has_duplication",
            "truth_single_copy",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError("Truth manifest is missing columns: " + ", ".join(sorted(missing)))
        for row_index, row in enumerate(reader, start=1):
            species_path = (manifest_path.parent / row["species_tree"]).resolve()
            gene_path = (manifest_path.parent / row["gene_tree"]).resolve()
            mapping = row.get("mapping") or ""
            mapping_path = (manifest_path.parent / mapping).resolve() if mapping else None
            delimiter = row.get("delimiter") or None
            species_field = int(row.get("species_field") or 0)
            expected = tuple(value for value in row["expected_statuses"].split("|") if value)
            result = analyze_locus(
                species_path,
                gene_path,
                mapping_path=mapping_path,
                family=row["benchmark_id"],
                delimiter=delimiter,
                species_field=species_field,
                min_coverage=min_coverage,
                discordance_review_threshold=discordance_review_threshold,
            )
            truth_has_duplication = _read_bool(row["truth_has_duplication"])
            truth_single_copy = _read_bool(row["truth_single_copy"])
            records.append(
                BenchmarkRecord(
                    benchmark_id=row["benchmark_id"],
                    scenario=row["scenario"],
                    replicate=int(row.get("replicate") or row_index),
                    seed=int(row.get("seed") or 0),
                    taxa=int(row.get("taxa") or result.species_total),
                    truth_class=row.get("truth_class") or row["scenario"],
                    truth_has_duplication=truth_has_duplication,
                    truth_single_copy=truth_single_copy,
                    expected_statuses="|".join(expected),
                    observed_status=result.verdict.status,
                    acceptable=result.verdict.status in expected,
                    false_resolve=result.verdict.status == "RESOLVE" and not truth_has_duplication,
                    primary_concern=result.verdict.primary_concern,
                    coverage=result.coverage,
                    max_copies=result.max_copies,
                    overlap_duplications=result.species_overlap_duplications,
                    lca_duplications=result.reconciliation.lca_duplications,
                    inferred_losses=result.reconciliation.inferred_losses,
                    normalized_rf=result.topology.normalized_rf,
                    species_tree=str(species_path),
                    gene_tree=str(gene_path),
                    mapping=str(mapping_path) if mapping_path else "",
                    delimiter=delimiter or "",
                    species_field=species_field,
                )
            )
    write_benchmark_outputs(records, outdir, config=None, source="manifest")
    return records


def _safe_divide(numerator: int, denominator: int) -> float:
    return numerator / denominator if denominator else 0.0


def _class_metrics(records: list[BenchmarkRecord]) -> list[dict[str, str | int | float]]:
    metrics: list[dict[str, str | int | float]] = []
    for status in ("PASS", "REVIEW", "RESOLVE"):
        tp = sum(
            record.observed_status == status and status in record.expected_statuses.split("|")
            for record in records
        )
        fp = sum(
            record.observed_status == status and status not in record.expected_statuses.split("|")
            for record in records
        )
        fn = sum(
            record.observed_status != status and status in record.expected_statuses.split("|")
            for record in records
        )
        precision = _safe_divide(tp, tp + fp)
        recall = _safe_divide(tp, tp + fn)
        f1 = _safe_divide(2 * precision * recall, precision + recall) if precision + recall else 0.0
        metrics.append(
            {
                "status": status,
                "tp": tp,
                "fp": fp,
                "fn": fn,
                "precision": precision,
                "recall": recall,
                "f1": f1,
            }
        )
    return metrics


def _scenario_metrics(records: list[BenchmarkRecord]) -> list[dict[str, str | int | float]]:
    grouped: dict[str, list[BenchmarkRecord]] = {}
    for record in records:
        grouped.setdefault(record.scenario, []).append(record)
    rows: list[dict[str, str | int | float]] = []
    for scenario, values in sorted(grouped.items()):
        rows.append(
            {
                "scenario": scenario,
                "n": len(values),
                "acceptable": sum(value.acceptable for value in values),
                "accuracy": _safe_divide(sum(value.acceptable for value in values), len(values)),
                "false_resolve": sum(value.false_resolve for value in values),
            }
        )
    return rows


def benchmark_summary(records: list[BenchmarkRecord]) -> dict[str, float | int]:
    total = len(records)
    acceptable = sum(record.acceptable for record in records)
    true_resolve = [
        record
        for record in records
        if record.truth_has_duplication and not record.truth_single_copy
    ]
    no_dup = [record for record in records if not record.truth_has_duplication]
    single_copy_no_dup = [
        record
        for record in records
        if record.truth_single_copy and not record.truth_has_duplication
    ]
    resolve_tp = sum(record.observed_status == "RESOLVE" for record in true_resolve)
    resolve_fp = sum(record.observed_status == "RESOLVE" for record in no_dup)
    return {
        "total": total,
        "acceptable": acceptable,
        "accuracy": _safe_divide(acceptable, total),
        "resolve_sensitivity_multicopy_truth": _safe_divide(resolve_tp, len(true_resolve)),
        "false_resolve_rate_no_duplication_truth": _safe_divide(resolve_fp, len(no_dup)),
        "false_resolve_rate_singlecopy_no_dup_truth": _safe_divide(
            sum(record.observed_status == "RESOLVE" for record in single_copy_no_dup),
            len(single_copy_no_dup),
        ),
    }


def _write_tsv(path: Path, rows: Iterable[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def _plot_scenario_accuracy(rows: list[dict[str, str | int | float]], path: Path) -> None:
    if not rows:
        return
    import matplotlib.pyplot as plt

    labels = [str(row["scenario"]) for row in rows]
    values = [float(row["accuracy"]) for row in rows]
    fig, ax = plt.subplots(figsize=(max(7, len(labels) * 1.2), 4.5))
    ax.bar(labels, values)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Acceptable verdict fraction")
    ax.set_xlabel("Known-truth scenario")
    ax.set_title("GeneSpeciesTreeVerdict benchmark performance")
    ax.tick_params(axis="x", rotation=35)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def write_benchmark_outputs(
    records: list[BenchmarkRecord],
    outdir: str | Path,
    *,
    config: BenchmarkConfig | None,
    source: str,
) -> Path:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rows = [asdict(record) for record in records]
    fields = list(BenchmarkRecord.__dataclass_fields__)
    _write_tsv(outdir / "locus_results.tsv", rows, fields)
    _write_tsv(
        outdir / "truth_manifest.tsv",
        rows,
        [
            "benchmark_id",
            "scenario",
            "replicate",
            "seed",
            "taxa",
            "truth_class",
            "truth_has_duplication",
            "truth_single_copy",
            "expected_statuses",
            "species_tree",
            "gene_tree",
            "mapping",
            "delimiter",
            "species_field",
        ],
    )
    failures = [row for row in rows if not row["acceptable"]]
    _write_tsv(outdir / "failure_cases.tsv", failures, fields)

    scenario_rows = _scenario_metrics(records)
    _write_tsv(
        outdir / "scenario_metrics.tsv",
        scenario_rows,
        ["scenario", "n", "acceptable", "accuracy", "false_resolve"],
    )
    class_rows = _class_metrics(records)
    _write_tsv(
        outdir / "class_metrics.tsv",
        class_rows,
        ["status", "tp", "fp", "fn", "precision", "recall", "f1"],
    )
    confusion_rows: list[dict[str, str | int]] = []
    confusion = Counter(
        (record.expected_statuses.split("|")[0], record.observed_status) for record in records
    )
    for expected in ("PASS", "REVIEW", "RESOLVE"):
        for observed in ("PASS", "REVIEW", "RESOLVE"):
            confusion_rows.append(
                {
                    "expected": expected,
                    "observed": observed,
                    "count": confusion.get((expected, observed), 0),
                }
            )
    _write_tsv(
        outdir / "confusion_matrix.tsv",
        confusion_rows,
        ["expected", "observed", "count"],
    )

    summary = benchmark_summary(records)
    _write_tsv(
        outdir / "performance_summary.tsv",
        [{"metric": key, "value": value} for key, value in summary.items()],
        ["metric", "value"],
    )

    report = [
        "# GeneSpeciesTreeVerdict benchmark report",
        "",
        f"**Source:** `{source}`",
        "",
        "## Headline metrics",
        "",
        "| Metric | Value |",
        "|---|---:|",
    ]
    for key, value in summary.items():
        rendered = f"{value:.3f}" if isinstance(value, float) else str(value)
        report.append(f"| {key} | {rendered} |")
    report.extend(
        [
            "",
            "## Scenario performance",
            "",
            "| Scenario | n | Accuracy | False RESOLVE |",
            "|---|---:|---:|---:|",
        ]
    )
    for row in scenario_rows:
        report.append(
            f"| {row['scenario']} | {row['n']} | {row['accuracy']:.3f} | {row['false_resolve']} |"
        )
    report.extend(
        [
            "",
            "## Interpretation boundary",
            "",
            "The native benchmark is a structural stress test with known copy number, coverage, and topology manipulations. "
            "It is not a full generative model of ILS, HGT, recombination, sequence evolution, or gene-tree estimation. "
            "Process-level validation should use the external truth-manifest mode with simulators such as SimPhy or Zombi.",
        ]
    )
    if config is not None:
        report.extend(
            [
                "",
                "## Reproducibility",
                "",
                f"Base seed: `{config.seed}`  ",
                f"Replicates per taxon/scenario cell: `{config.replicates}`  ",
                f"Taxon counts: `{', '.join(map(str, config.taxa))}`  ",
                f"Scenarios: `{', '.join(config.scenarios)}`",
            ]
        )
    (outdir / "benchmark_report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    _plot_scenario_accuracy(scenario_rows, outdir / "scenario_accuracy.png")
    return outdir
