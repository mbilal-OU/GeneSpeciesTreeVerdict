from __future__ import annotations

TOPICS: dict[str, dict[str, str]] = {
    "discordance": {
        "title": "Gene-tree/species-tree discordance",
        "answer": (
            "A gene tree can differ from a species tree for biological and technical reasons. Common "
            "possibilities include duplication/loss, incomplete lineage sorting, horizontal transfer, "
            "recombination, gene-tree error, model misspecification, taxon sampling, or an uncertain "
            "reference species tree. Discordance is an observation, not a diagnosis."
        ),
        "decision": (
            "Inspect copy number, branch support, the pattern of conflict across loci, and model-specific "
            "evidence before assigning a cause or excluding the locus."
        ),
    },
    "paralogy": {
        "title": "Orthology and paralogy",
        "answer": (
            "Orthologs are gene copies related by speciation; paralogs are related by duplication. These "
            "terms describe relationships between copies, not whole gene families. A duplicated family can "
            "contain orthologous subgroups across species."
        ),
        "decision": (
            "Direct multi-copy structure is a strong reason to resolve the family before conventional "
            "single-copy concatenation. A discordant single-copy gene tree is not enough to call paralogy."
        ),
    },
    "single-copy": {
        "title": "Why one sequence per species does not prove orthology",
        "answer": (
            "An ancient duplication followed by differential loss can leave exactly one surviving copy in "
            "each sampled species while those surviving copies belong to different paralog lineages. This is "
            "often called pseudo-orthologous sampling."
        ),
        "decision": (
            "Treat single-copy status as a useful structural filter, not as proof of orthology. Reconciliation "
            "and broader gene-family context can still matter."
        ),
    },
    "rf": {
        "title": "Robinson-Foulds distance",
        "answer": (
            "RF distance counts bipartitions present in one tree but not the other after restricting the trees "
            "to comparable taxa. Normalized RF scales this disagreement by the number of splits being compared."
        ),
        "decision": (
            "RF measures topological difference. It does not identify the biological process that caused the "
            "difference and should not be used as an automatic paralogy filter."
        ),
    },
    "reconciliation": {
        "title": "Gene-tree/species-tree reconciliation",
        "answer": (
            "Reconciliation maps a gene-family history into a species-tree history and asks what events are "
            "needed under a chosen model, such as duplication-loss or duplication-transfer-loss."
        ),
        "decision": (
            "Interpret inferred events conditionally on the model, rooting, tree uncertainty, and sampling. "
            "A reconciliation is an evolutionary explanation under assumptions, not direct observation."
        ),
    },
    "ils": {
        "title": "Incomplete lineage sorting (ILS)",
        "answer": (
            "ILS occurs when ancestral allelic lineages fail to coalesce within the immediately preceding "
            "species branches, allowing a gene genealogy to differ from the species history. It is especially "
            "relevant around short internal species-tree branches and requires multi-locus context to evaluate."
        ),
        "decision": (
            "A single discordant gene tree cannot establish ILS. Look for repeated genome-wide patterns and "
            "use an appropriate multispecies-coalescent framework when ILS is biologically plausible."
        ),
    },
    "hgt": {
        "title": "Horizontal gene transfer (HGT)",
        "answer": (
            "HGT can place a gene history in conflict with vertical species relationships, particularly in "
            "microbial datasets. However, similar topological conflict can arise from other processes."
        ),
        "decision": (
            "Do not infer HGT from RF distance alone. Evaluate support, gene-family context, genomic context "
            "where available, and a transfer-aware reconciliation model when appropriate."
        ),
    },
    "marker-selection": {
        "title": "Using a locus in a conventional concatenated marker set",
        "answer": (
            "Core prevalence, copy number, orthology, alignment quality, substitutional behavior, and gene-tree "
            "history are separate properties. A core gene is not automatically a good phylogenomic marker."
        ),
        "decision": (
            "For a conventional single-copy concatenation workflow, favor loci with adequate coverage and "
            "defensible orthology, then perform sequence/alignment/model QC separately."
        ),
    },
}


def topic_names() -> list[str]:
    return sorted(TOPICS)


def explain_topic(topic: str) -> dict[str, str]:
    key = topic.strip().lower()
    if key not in TOPICS:
        options = ", ".join(topic_names())
        raise ValueError(f"Unknown topic '{topic}'. Available topics: {options}")
    return TOPICS[key]
