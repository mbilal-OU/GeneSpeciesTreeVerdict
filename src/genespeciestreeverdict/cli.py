from __future__ import annotations

from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.console import Console
from rich.table import Table

from . import __version__
from .analysis import analyze_locus
from .batch import analyze_directory
from .plotting import plot_tree_comparison
from .reporting import write_locus_outputs
from .tutorial import write_tutorial_dataset

app = typer.Typer(
    add_completion=False,
    no_args_is_help=True,
    help="Diagnose gene-tree/species-tree discordance and evaluate locus handling.",
)
console = Console()


def _show_result(result) -> None:
    table = Table(title=f"GeneSpeciesTreeVerdict — {result.family}")
    table.add_column("Metric")
    table.add_column("Value")
    table.add_row("Verdict", result.verdict.status)
    table.add_row("Topology", result.verdict.topology_signal)
    table.add_row("Coverage", f"{result.coverage:.1%}")
    table.add_row("Maximum copies/species", str(result.max_copies))
    table.add_row("Species-overlap duplications", str(result.species_overlap_duplications))
    table.add_row("LCA-DL duplications", str(result.reconciliation.lca_duplications))
    table.add_row("Inferred DL losses", str(result.reconciliation.inferred_losses))
    table.add_row(
        "Normalized RF",
        f"{result.topology.normalized_rf:.3f}"
        if result.topology.normalized_rf is not None
        else "NA",
    )
    console.print(table)
    console.print("\n[bold]Recommended action[/bold]")
    console.print(result.verdict.recommended_action)
    console.print("\n[dim]" + result.verdict.caution + "[/dim]")


@app.command()
def analyze(
    species_tree: Annotated[Path, typer.Option("--species-tree", exists=True, dir_okay=False)],
    gene_tree: Annotated[Path, typer.Option("--gene-tree", exists=True, dir_okay=False)],
    mapping: Annotated[
        Optional[Path], typer.Option("--mapping", exists=True, dir_okay=False)
    ] = None,
    family: Annotated[Optional[str], typer.Option("--family")] = None,
    delimiter: Annotated[Optional[str], typer.Option("--delimiter")] = None,
    species_field: Annotated[int, typer.Option("--species-field", min=0)] = 0,
    min_coverage: Annotated[float, typer.Option("--min-coverage", min=0.01, max=1.0)] = 0.95,
    rooted_comparison: Annotated[
        bool,
        typer.Option(
            "--rooted-comparison", help="Compare rooted clades instead of unrooted splits."
        ),
    ] = False,
    no_reconcile: Annotated[
        bool,
        typer.Option("--no-reconcile", help="Skip rooted duplication-loss LCA reconciliation."),
    ] = False,
    discordance_review_threshold: Annotated[
        Optional[float],
        typer.Option(
            "--discordance-review-threshold",
            min=0.0,
            max=1.0,
            help="Optional nRF threshold that upgrades PASS to REVIEW. No threshold is applied by default.",
        ),
    ] = None,
    outdir: Annotated[Path, typer.Option("--outdir")] = Path("gstv_result"),
    plot: Annotated[bool, typer.Option("--plot/--no-plot")] = True,
) -> None:
    """Analyze one gene tree against a reference species tree."""

    result = analyze_locus(
        species_tree,
        gene_tree,
        mapping_path=mapping,
        family=family,
        delimiter=delimiter,
        species_field=species_field,
        min_coverage=min_coverage,
        rooted_comparison=rooted_comparison,
        reconcile=not no_reconcile,
        discordance_review_threshold=discordance_review_threshold,
    )
    write_locus_outputs(result, outdir)
    if plot:
        plot_tree_comparison(
            species_tree, gene_tree, outdir / "tree_comparison.png", title=result.family
        )
    _show_result(result)
    console.print(f"\nOutputs: [bold]{outdir}[/bold]")


@app.command(name="batch")
def batch_command(
    species_tree: Annotated[Path, typer.Option("--species-tree", exists=True, dir_okay=False)],
    gene_tree_dir: Annotated[Path, typer.Option("--gene-tree-dir", exists=True, file_okay=False)],
    mapping: Annotated[
        Optional[Path], typer.Option("--mapping", exists=True, dir_okay=False)
    ] = None,
    delimiter: Annotated[Optional[str], typer.Option("--delimiter")] = None,
    species_field: Annotated[int, typer.Option("--species-field", min=0)] = 0,
    min_coverage: Annotated[float, typer.Option("--min-coverage", min=0.01, max=1.0)] = 0.95,
    rooted_comparison: Annotated[bool, typer.Option("--rooted-comparison")] = False,
    no_reconcile: Annotated[bool, typer.Option("--no-reconcile")] = False,
    discordance_review_threshold: Annotated[
        Optional[float], typer.Option("--discordance-review-threshold", min=0.0, max=1.0)
    ] = None,
    recursive: Annotated[bool, typer.Option("--recursive")] = False,
    outdir: Annotated[Path, typer.Option("--outdir")] = Path("gstv_batch"),
) -> None:
    """Analyze a directory of gene trees and summarize reference-tree support."""

    results, errors = analyze_directory(
        species_tree,
        gene_tree_dir,
        mapping_path=mapping,
        delimiter=delimiter,
        species_field=species_field,
        min_coverage=min_coverage,
        rooted_comparison=rooted_comparison,
        reconcile=not no_reconcile,
        discordance_review_threshold=discordance_review_threshold,
        recursive=recursive,
        outdir=outdir,
    )
    counts: dict[str, int] = {}
    for result in results:
        counts[result.verdict.status] = counts.get(result.verdict.status, 0) + 1
    table = Table(title="Batch summary")
    table.add_column("Verdict")
    table.add_column("Loci", justify="right")
    for status, count in sorted(counts.items()):
        table.add_row(status, str(count))
    table.add_row("ERROR", str(len(errors)))
    console.print(table)
    console.print(f"Outputs: [bold]{outdir}[/bold]")


@app.command()
def tutorial(
    outdir: Annotated[Path, typer.Option("--outdir")] = Path("gstv_tutorial"),
) -> None:
    """Write the deterministic 10-genome teaching dataset."""

    write_tutorial_dataset(outdir)
    console.print(f"Tutorial dataset written to [bold]{outdir}[/bold]")
    console.print(
        "Try: gstv batch --species-tree "
        f"{outdir}/species_tree.nwk --gene-tree-dir {outdir}/gene_trees "
        f"--mapping {outdir}/mapping.tsv --outdir {outdir}/results"
    )


@app.command()
def version() -> None:
    """Print the installed version."""

    console.print(__version__)


if __name__ == "__main__":
    app()
