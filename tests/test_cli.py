from typer.testing import CliRunner

from genespeciestreeverdict.cli import app


runner = CliRunner()


def test_version_command():
    result = runner.invoke(app, ["version"])
    assert result.exit_code == 0
    assert "0.1.0" in result.stdout


def test_tutorial_command(tmp_path):
    out = tmp_path / "demo"
    result = runner.invoke(app, ["tutorial", "--outdir", str(out)])
    assert result.exit_code == 0
    assert (out / "species_tree.nwk").exists()
    assert (out / "mapping.tsv").exists()


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
    assert (out / "result.json").exists()


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
