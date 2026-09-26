from genespeciestreeverdict.analysis import analyze_locus


def _run(tutorial_dir, family):
    return analyze_locus(
        tutorial_dir / "species_tree.nwk",
        tutorial_dir / "gene_trees" / f"{family}.nwk",
        mapping_path=tutorial_dir / "mapping.tsv",
    )


def test_concordant_ortholog_passes(tutorial_dir):
    result = _run(tutorial_dir, "geneA_ortholog_concordant")
    assert result.verdict.status == "PASS"
    assert result.verdict.topology_signal == "CONCORDANT"
    assert result.max_copies == 1
    assert result.topology.normalized_rf == 0.0
    assert result.reconciliation.lca_duplications == 0


def test_singlecopy_discordance_is_not_called_paralogy(tutorial_dir):
    result = _run(tutorial_dir, "geneD_singlecopy_discordant")
    assert result.max_copies == 1
    assert result.species_overlap_duplications == 0
    assert result.verdict.topology_signal == "DISCORDANT"
    assert result.verdict.status == "REVIEW"
    assert any("model-dependent" in reason for reason in result.verdict.reasons)


def test_recent_duplication_requires_resolution(tutorial_dir):
    result = _run(tutorial_dir, "geneB_recent_duplication")
    assert result.max_copies == 2
    assert result.multicopy_species == ["G3"]
    assert result.species_overlap_duplications >= 1
    assert result.verdict.status == "RESOLVE"
    assert not result.topology.comparable


def test_ancient_duplication_loss_requires_resolution(tutorial_dir):
    result = _run(tutorial_dir, "geneC_ancient_duplication_loss")
    assert result.max_copies == 2
    assert result.reconciliation.lca_duplications >= 1
    assert result.reconciliation.inferred_losses >= 1
    assert result.verdict.status == "RESOLVE"


def test_mixed_pseudoortholog_is_reviewed_not_declared_paralog(tutorial_dir):
    result = _run(tutorial_dir, "geneP_mixed_pseudoortholog")
    assert result.max_copies == 1
    assert result.species_overlap_duplications == 0
    assert result.topology.normalized_rf == 1.0
    assert result.verdict.status == "REVIEW"


def test_user_discordance_threshold_can_upgrade_pass_to_review(tutorial_dir):
    result = analyze_locus(
        tutorial_dir / "species_tree.nwk",
        tutorial_dir / "gene_trees" / "geneA_ortholog_concordant.nwk",
        mapping_path=tutorial_dir / "mapping.tsv",
        discordance_review_threshold=0.0,
    )
    assert result.verdict.status == "REVIEW"


def test_reconciliation_can_be_disabled(tutorial_dir):
    result = analyze_locus(
        tutorial_dir / "species_tree.nwk",
        tutorial_dir / "gene_trees" / "geneD_singlecopy_discordant.nwk",
        mapping_path=tutorial_dir / "mapping.tsv",
        reconcile=False,
    )
    assert not result.reconciliation.attempted
    assert result.verdict.status == "PASS"
