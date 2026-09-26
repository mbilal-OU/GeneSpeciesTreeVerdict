import json

from genespeciestreeverdict.analysis import analyze_locus
from genespeciestreeverdict.plotting import plot_tree_comparison
from genespeciestreeverdict.reporting import write_locus_outputs


def test_single_locus_reporting_and_plot(tutorial_dir, tmp_path):
    result = analyze_locus(
        tutorial_dir / "species_tree.nwk",
        tutorial_dir / "gene_trees" / "geneA_ortholog_concordant.nwk",
        mapping_path=tutorial_dir / "mapping.tsv",
    )
    out = tmp_path / "out"
    write_locus_outputs(result, out)
    plot_tree_comparison(
        tutorial_dir / "species_tree.nwk",
        tutorial_dir / "gene_trees" / "geneA_ortholog_concordant.nwk",
        out / "comparison.png",
        title="demo",
    )
    assert (out / "summary.tsv").exists()
    assert (out / "evidence_matrix.tsv").exists()
    assert (out / "next_analyses.tsv").exists()
    assert (out / "report.md").exists()
    assert (out / "comparison.png").exists()
    payload = json.loads((out / "result.json").read_text())
    assert payload["verdict"]["status"] == "PASS"
    assert payload["evidence"]
    assert payload["next_analyses"]
    report = (out / "report.md").read_text()
    assert "What the current evidence supports" in report
    assert "does **not** establish" in report
    assert "Next analyses" in report
