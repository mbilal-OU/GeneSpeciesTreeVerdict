from genespeciestreeverdict.analysis import analyze_locus
from genespeciestreeverdict.knowledge import explain_topic


def _run(tutorial_dir, family):
    return analyze_locus(
        tutorial_dir / "species_tree.nwk",
        tutorial_dir / "gene_trees" / f"{family}.nwk",
        mapping_path=tutorial_dir / "mapping.tsv",
    )


def test_concordant_locus_has_transparent_support_and_routine_qc(tutorial_dir):
    result = _run(tutorial_dir, "geneA_ortholog_concordant")
    evidence = {item.key: item for item in result.evidence}
    assert evidence["copy_number"].signal == "SUPPORT"
    assert "does not guarantee orthology" in evidence["copy_number"].limitation
    assert evidence["topology"].signal == "SUPPORT"
    assert result.next_analyses[0].priority == "ROUTINE"
    assert "otherwise suitable" in result.next_analyses[0].question


def test_multicopy_locus_routes_to_orthology_resolution(tutorial_dir):
    result = _run(tutorial_dir, "geneB_recent_duplication")
    evidence = {item.key: item for item in result.evidence}
    assert evidence["copy_number"].signal == "STRONG_FLAG"
    assert result.verdict.primary_concern == "MULTICOPY_OR_DUPLICATION_SIGNAL"
    assert result.next_analyses[0].priority == "HIGH"
    assert "orthologous" in result.next_analyses[0].question
    assert any("GeneRax" in tool for tool in result.next_analyses[0].tools)


def test_singlecopy_discordance_does_not_claim_paralogy(tutorial_dir):
    result = _run(tutorial_dir, "geneD_singlecopy_discordant")
    evidence = {item.key: item for item in result.evidence}
    assert evidence["topology"].signal == "FLAG"
    assert "not cause-specific" in evidence["topology"].limitation
    assert any(
        "uniquely demonstrated" in statement for statement in result.verdict.unsupported_conclusions
    )
    assert any("well supported" in item.question for item in result.next_analyses)


def test_pseudoorthology_example_preserves_singlecopy_limitation(tutorial_dir):
    result = _run(tutorial_dir, "geneP_mixed_pseudoortholog")
    copy_evidence = next(item for item in result.evidence if item.key == "copy_number")
    assert "pseudo-orthologous" in copy_evidence.limitation
    assert result.verdict.status == "REVIEW"


def test_knowledge_base_explains_rf_without_causal_overclaim():
    entry = explain_topic("rf")
    assert "topological difference" in entry["decision"]
    assert "automatic paralogy filter" in entry["decision"]
