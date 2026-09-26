from pathlib import Path

import pytest

from genespeciestreeverdict.io import MappingRecord, mapping_for_family, read_mapping_table, read_tree


def test_three_column_mapping_with_header(tutorial_dir):
    records = read_mapping_table(tutorial_dir / "mapping.tsv")
    assert records
    assert all(isinstance(record, MappingRecord) for record in records)
    assert {record.family for record in records if record.family} >= {"geneA_ortholog_concordant"}


def test_identity_mapping_when_gene_and_species_labels_match(tmp_path: Path):
    species = tmp_path / "species.nwk"
    gene = tmp_path / "gene.nwk"
    species.write_text("((A,B),(C,D));\n")
    gene.write_text("((A,B),(C,D));\n")
    st = read_tree(species)
    gt = read_tree(gene)
    assert mapping_for_family(gt, st) == {"A": "A", "B": "B", "C": "C", "D": "D"}


def test_missing_mapping_is_actionable(tutorial_dir, tmp_path: Path):
    bad = tmp_path / "bad.tsv"
    bad.write_text("gene_id\tspecies_id\nG1|A\tG1\n")
    st = read_tree(tutorial_dir / "species_tree.nwk")
    gt = read_tree(tutorial_dir / "gene_trees" / "geneA_ortholog_concordant.nwk")
    records = read_mapping_table(bad)
    with pytest.raises(ValueError, match="missing from the mapping table"):
        mapping_for_family(gt, st, records=records, family="geneA_ortholog_concordant")


def test_delimiter_mapping(tutorial_dir):
    st = read_tree(tutorial_dir / "species_tree.nwk")
    gt = read_tree(tutorial_dir / "gene_trees" / "geneA_ortholog_concordant.nwk")
    mapping = mapping_for_family(gt, st, delimiter="|", species_field=0)
    assert mapping["G1|A"] == "G1"
    assert len(mapping) == 10


def test_bad_mapping_width_raises(tmp_path: Path):
    bad = tmp_path / "mapping.tsv"
    bad.write_text("a\tb\tc\td\n")
    with pytest.raises(ValueError, match="2 columns"):
        read_mapping_table(bad)


def test_unknown_species_mapping_raises(tutorial_dir, tmp_path: Path):
    bad = tmp_path / "mapping.tsv"
    lines = ["gene_id\tspecies_id"]
    for i in range(1, 11):
        species = "NOPE" if i == 10 else f"G{i}"
        lines.append(f"G{i}|A\t{species}")
    bad.write_text("\n".join(lines) + "\n")
    st = read_tree(tutorial_dir / "species_tree.nwk")
    gt = read_tree(tutorial_dir / "gene_trees" / "geneA_ortholog_concordant.nwk")
    records = read_mapping_table(bad)
    with pytest.raises(ValueError, match="absent from the species tree"):
        mapping_for_family(gt, st, records=records)
