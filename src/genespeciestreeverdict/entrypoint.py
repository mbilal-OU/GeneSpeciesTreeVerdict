from __future__ import annotations

from pathlib import Path
from typing import Annotated, Optional

import typer
from rich.table import Table

from .benchmark import (
    benchmark_summary,
    load_benchmark_config,
    run_manifest_benchmark,
    run_native_benchmark,
)
from .cli import app, console


@app.command(name="benchmark")
def benchmark_command(
    config: Annotated[Optional[Path], typer.Option("--config", exists=True, dir_okay=False)] = None,
    truth_manifest: Annotated[
        Optional[Path], typer.Option("--truth-manifest", exists=True, dir_okay=False)
    ] = None,
    outdir: Annotated[Path, typer.Option("--outdir")] = Path("gstv_benchmark"),
    min_coverage: Annotated[float, typer.Option("--min-coverage", min=0.01, max=1.0)] = 0.95,
    discordance_review_threshold: Annotated[
        Optional[float],
        typer.Option("--discordance-review-threshold", min=0.0, max=1.0),
    ] = None,
    fail_on_mismatch: Annotated[
        bool,
        typer.Option(
            "--fail-on-mismatch",
            help="Exit non-zero if any benchmark case falls outside its expected workflow state(s).",
        ),
    ] = False,
) -> None:
    """Benchmark GSTV decisions against known truth or an external simulation manifest."""

    if config is not None and truth_manifest is not None:
        console.print("[red]Use either --config or --truth-manifest, not both.[/red]")
        raise typer.Exit(code=2)

    if truth_manifest is not None:
        records = run_manifest_benchmark(
            truth_manifest,
            outdir,
            min_coverage=min_coverage,
            discordance_review_threshold=discordance_review_threshold,
        )
        mode = "external truth manifest"
    else:
        benchmark_config = load_benchmark_config(config)
        if min_coverage != 0.95:
            benchmark_config.min_coverage = min_coverage
        if discordance_review_threshold is not None:
            benchmark_config.discordance_review_threshold = discordance_review_threshold
        records = run_native_benchmark(benchmark_config, outdir)
        mode = "native known-truth structural stress test"

    summary = benchmark_summary(records)
    table = Table(title="GeneSpeciesTreeVerdict benchmark")
    table.add_column("Metric")
    table.add_column("Value", justify="right")
    table.add_row("Mode", mode)
    table.add_row("Cases", str(summary["total"]))
    table.add_row("Expected-state matches", str(summary["acceptable"]))
    table.add_row("Overall match fraction", f"{summary['accuracy']:.3f}")
    table.add_row(
        "RESOLVE sensitivity (sampled multicopy truth)",
        f"{summary['resolve_sensitivity_multicopy_truth']:.3f}",
    )
    table.add_row(
        "False RESOLVE (no-dup truth)",
        f"{summary['false_resolve_rate_no_duplication_truth']:.3f}",
    )
    table.add_row(
        "False RESOLVE (single-copy/no-dup truth)",
        f"{summary['false_resolve_rate_singlecopy_no_dup_truth']:.3f}",
    )
    console.print(table)
    console.print(f"Outputs: [bold]{outdir}[/bold]")
    console.print(
        "[dim]Native simulations validate structural decision logic, not mechanistic ILS/HGT inference. "
        "Use --truth-manifest with an independent simulator for process-level validation.[/dim]"
    )

    if fail_on_mismatch and summary["acceptable"] != summary["total"]:
        raise typer.Exit(code=1)


if __name__ == "__main__":
    app()
