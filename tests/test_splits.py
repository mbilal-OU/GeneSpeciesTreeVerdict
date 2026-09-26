from pathlib import Path

from genespeciestreeverdict.io import read_tree
from genespeciestreeverdict.splits import compare_topologies


def test_identical_four_taxon_tree_has_zero_rf(tmp_path: Path):
    a = tmp_path / "a.nwk"
    b = tmp_path / "b.nwk"
    a.write_text("((A,B),(C,D));\n")
    b.write_text("((A,B),(C,D));\n")
    species = read_tree(a)
    gene = read_tree(b)
    mapping = {x: x for x in "ABCD"}
    metrics = compare_topologies(species, gene, mapping)
    assert metrics.rf_distance == 0
    assert metrics.normalized_rf == 0.0


def test_conflicting_four_taxon_tree_has_maximum_unrooted_nrf(tmp_path: Path):
    a = tmp_path / "a.nwk"
    b = tmp_path / "b.nwk"
    a.write_text("((A,B),(C,D));\n")
    b.write_text("((A,C),(B,D));\n")
    species = read_tree(a)
    gene = read_tree(b)
    mapping = {x: x for x in "ABCD"}
    metrics = compare_topologies(species, gene, mapping)
    assert metrics.rf_distance == 2
    assert metrics.normalized_rf == 1.0
