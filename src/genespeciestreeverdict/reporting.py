from __future__ import annotations

import csv
import json
from collections import Counter
from pathlib import Path
from typing import Iterable

from .models import LocusResult


def write_locus_outputs(result: LocusResult, outdir: str | Path) -> Path:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    json_path = outdir / "result.json"
    json_path.write_text(json.dumps(result.to_dict(), indent=2) + "\n", encoding="utf-8")

    with (outdir / "summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "family",
                "verdict",
                "topology_signal",
                "coverage",
                "max_copies",
                "overlap_duplications",
                "lca_duplications",
                "inferred_losses",
                "rf_distance",
                "normalized_rf",
            ]
        )
        writer.writerow(
            [
                result.family,
                result.verdict.status,
                result.verdict.topology_signal,
                f"{result.coverage:.6f}",
                result.max_copies,
                result.species_overlap_duplications,
                result.reconciliation.lca_duplications,
                result.reconciliation.inferred_losses,
                result.topology.rf_distance if result.topology.rf_distance is not None else "NA",
                (
                    f"{result.topology.normalized_rf:.6f}"
                    if result.topology.normalized_rf is not None
                    else "NA"
                ),
            ]
        )

    report = [
        f"# GeneSpeciesTreeVerdict report: {result.family}",
        "",
        f"**Verdict:** `{result.verdict.status}`  ",
        f"**Topology signal:** `{result.verdict.topology_signal}`",
        "",
        "## Core metrics",
        "",
        "| Metric | Value |",
        "|---|---:|",
        f"| Species in reference tree | {result.species_total} |",
        f"| Gene-tree leaves | {result.gene_leaves} |",
        f"| Unique mapped species | {result.unique_species} |",
        f"| Coverage | {result.coverage:.1%} |",
        f"| Maximum copies/species | {result.max_copies} |",
        f"| Species-overlap duplication candidates | {result.species_overlap_duplications} |",
        f"| LCA-DL duplications | {result.reconciliation.lca_duplications} |",
        f"| LCA-DL inferred losses | {result.reconciliation.inferred_losses} |",
        f"| RF distance | {result.topology.rf_distance if result.topology.rf_distance is not None else 'NA'} |",
        (
            f"| Normalized RF | {result.topology.normalized_rf:.3f} |"
            if result.topology.normalized_rf is not None
            else "| Normalized RF | NA |"
        ),
        "",
        "## Why this verdict?",
        "",
    ]
    report.extend(f"- {reason}" for reason in result.verdict.reasons)
    report.extend(
        [
            "",
            "## Recommended action",
            "",
            result.verdict.recommended_action,
            "",
            "> " + result.verdict.caution,
        ]
    )
    if result.warnings:
        report.extend(["", "## Interpretation warnings", ""])
        report.extend(f"- {warning}" for warning in result.warnings)
    (outdir / "report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    return outdir


def write_batch_outputs(
    results: Iterable[LocusResult],
    errors: list[dict[str, str]],
    outdir: str | Path,
) -> Path:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    results = list(results)

    with (outdir / "locus_summary.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            [
                "family",
                "verdict",
                "topology_signal",
                "coverage",
                "gene_leaves",
                "unique_species",
                "max_copies",
                "multicopy_species",
                "overlap_duplications",
                "lca_duplications",
                "inferred_losses",
                "rf_distance",
                "normalized_rf",
                "species_split_recovery",
            ]
        )
        for result in results:
            writer.writerow(
                [
                    result.family,
                    result.verdict.status,
                    result.verdict.topology_signal,
                    f"{result.coverage:.6f}",
                    result.gene_leaves,
                    result.unique_species,
                    result.max_copies,
                    ",".join(result.multicopy_species),
                    result.species_overlap_duplications,
                    result.reconciliation.lca_duplications,
                    result.reconciliation.inferred_losses,
                    result.topology.rf_distance if result.topology.rf_distance is not None else "NA",
                    (
                        f"{result.topology.normalized_rf:.6f}"
                        if result.topology.normalized_rf is not None
                        else "NA"
                    ),
                    (
                        f"{result.topology.species_split_recovery:.6f}"
                        if result.topology.species_split_recovery is not None
                        else "NA"
                    ),
                ]
            )

    (outdir / "results.json").write_text(
        json.dumps([result.to_dict() for result in results], indent=2) + "\n",
        encoding="utf-8",
    )

    if errors:
        with (outdir / "errors.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["family", "path", "error"], delimiter="\t")
            writer.writeheader()
            writer.writerows(errors)

    # Aggregate support for full-taxon, single-copy trees only. This avoids pretending that
    # partially sampled gene trees test a reference split that they cannot actually resolve.
    full_single = [
        result
        for result in results
        if result.coverage == 1.0
        and result.max_copies == 1
        and result.topology.comparable
        and result.topology.species_splits
    ]
    if full_single:
        reference_splits = sorted(set(full_single[0].topology.species_splits))
        with (outdir / "reference_split_support.tsv").open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["reference_split", "genes_tested", "genes_supporting", "support_fraction"])
            for split in reference_splits:
                supporting = sum(split in result.topology.gene_splits for result in full_single)
                writer.writerow([split, len(full_single), supporting, f"{supporting / len(full_single):.6f}"])

        alternatives: Counter[str] = Counter()
        for result in full_single:
            alternatives.update(result.topology.gene_only_splits)
        with (outdir / "alternative_splits.tsv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["alternative_split", "gene_count", "fraction_of_tested_genes"])
            for split, count in alternatives.most_common():
                writer.writerow([split, count, f"{count / len(full_single):.6f}"])

    verdict_counts = Counter(result.verdict.status for result in results)
    report = [
        "# GeneSpeciesTreeVerdict batch report",
        "",
        f"Analyzed **{len(results)}** gene trees; **{len(errors)}** failed validation/analysis.",
        "",
        "## Verdict summary",
        "",
        "| Verdict | Loci |",
        "|---|---:|",
    ]
    report.extend(f"| {status} | {count} |" for status, count in sorted(verdict_counts.items()))
    report.extend(
        [
            "",
            "## How to read this report",
            "",
            "- `PASS`: structurally suitable candidate for conventional single-copy concatenation.",
            "- `REVIEW`: requires interpretation before inclusion; discordance is not automatically paralogy.",
            "- `RESOLVE`: multi-copy/duplication evidence means orthology should be resolved before conventional concatenation.",
            "",
            "For full-taxon single-copy loci, `reference_split_support.tsv` asks the reverse question too: "
            "how consistently do independent gene trees support each branch of the supplied species tree?",
        ]
    )
    (outdir / "report.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    return outdir
