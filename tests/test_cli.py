from typer.testing import CliRunner

from genespeciestreeverdict.entrypoint import app


runner = CliRunner()


def test_version_command():
    result = runner.invoke(app, ["version"])
    assert result.exit_code == 0
    assert "0.3.0" in result.stdout


def test_tutorial_command(tmp_path):
    out = tmp_path / "demo"
    result = runner.invoke(app, ["tutorial", "--outdir", str(out)])
    assert result.exit_code == 0
    assert (out / "species_tree.nwk").exists()
    assert (out / "mapping.tsv").exists()


def test_explain_command():
    result = runner.invoke(app, ["explain", "paralogy"])
    assert result.exit_code == 0
    assert "Orthology and paralogy" in result.stdout
    assert "discordant single-copy" in result.stdout


def test_explain_list_command():
    result = runner.invoke(app, ["explain", "list"])
    assert result.exit_code == 0
    assert "discordance" in result.stdout
    assert "marker-selection" in result.stdout


def test_analyze_cli(tutorial_dir, tmp_path):
    out = tmp_path / "single"
    result = runner.invoke(
        app,
        [
            "analyze",
            "--species-tree",
            str(tutorial_dir / "species_tree.nwk"),
            "--gene-tree",
            str(tutorial_dir / "gene_trees" / "geneA_ortholog_concordant.nwk"),
            "--mapping",
            str(tutorial_dir / "mapping.tsv"),
            "--outdir",
            str(out),
            "--no-plot",
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert "PASS" in result.stdout
    assert "Top next question" in result.stdout
    assert (out / "result.json").exists()
    assert (out / "evidence_matrix.tsv").exists()
    assert (out / "next_analyses.tsv").exists()


def test_batch_cli(tutorial_dir, tmp_path):
    out = tmp_path / "batch"
    result = runner.invoke(
        app,
        [
            "batch",
            "--species-tree",
            str(tutorial_dir / "species_tree.nwk"),
            "--gene-tree-dir",
            str(tutorial_dir / "gene_trees"),
            "--mapping",
            str(tutorial_dir / "mapping.tsv"),
            "--outdir",
            str(out),
        ],
    )
    assert result.exit_code == 0, result.stdout
    assert "RESOLVE" in result.stdout
    assert (out / "locus_summary.tsv").exists()
    assert (out / "evidence_matrix.tsv").exists()
    assert (out / "recommendations.tsv").exists()


def test_benchmark_cli(tmp_path):
    config = tmp_path / "benchmark.yaml"
    config.write_text(
        "seed: 44\n"
        "replicates: 1\n"
        "taxa: [8]\n"
        "scenarios: [concordant_singlecopy, multicopy_duplication]\n",
        encoding="utf-8",
    )
    out = tmp_path / "benchmark"
    result = runner.invoke(
        app,
        ["benchmark", "--config", str(config), "--outdir", str(out)],
    )

    assert result.exit_code == 0, result.stdout
    assert "GeneSpeciesTreeVerdict benchmark" in result.stdout
    assert "False RESOLVE" in result.stdout
    assert (out / "performance_summary.tsv").exists()
    assert (out / "truth_manifest.tsv").exists()
    assert (out / "scenario_accuracy.png").exists()
