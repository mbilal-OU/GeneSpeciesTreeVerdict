from genespeciestreeverdict.batch import analyze_directory


def test_batch_outputs_reference_branch_evidence(tutorial_dir, tmp_path):
    out = tmp_path / "results"
    results, errors = analyze_directory(
        tutorial_dir / "species_tree.nwk",
        tutorial_dir / "gene_trees",
        mapping_path=tutorial_dir / "mapping.tsv",
        outdir=out,
    )
    assert len(results) == 5
    assert errors == []
    assert (out / "locus_summary.tsv").exists()
    assert (out / "results.json").exists()
    assert (out / "reference_split_support.tsv").exists()
    assert (out / "alternative_splits.tsv").exists()


def test_batch_records_bad_tree_as_error(tutorial_dir, tmp_path):
    gene_dir = tmp_path / "genes"
    gene_dir.mkdir()
    # a valid identity-mapped tree and an invalid file
    (gene_dir / "good.nwk").write_text("(((G1,G2),(G3,G4)),((G5,G6),((G7,G8),(G9,G10))));\n")
    (gene_dir / "bad.nwk").write_text("not a newick tree\n")
    out = tmp_path / "out"
    results, errors = analyze_directory(
        tutorial_dir / "species_tree.nwk",
        gene_dir,
        outdir=out,
    )
    assert len(results) == 1
    assert len(errors) == 1
    assert (out / "errors.tsv").exists()
