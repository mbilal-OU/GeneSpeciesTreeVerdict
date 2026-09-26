from __future__ import annotations

from pathlib import Path


SPECIES_TREE = "(((G1:0.10,G2:0.10)98:0.12,(G3:0.09,G4:0.09)97:0.13)99:0.20,((G5:0.08,G6:0.08)96:0.14,((G7:0.07,G8:0.07)95:0.10,(G9:0.06,G10:0.06)94:0.11)98:0.12)99:0.18);\n"

GENE_TREES = {
    "geneA_ortholog_concordant": "(((G1|A:0.14,G2|A:0.11)99:0.10,(G3|A:0.08,G4|A:0.12)96:0.16)98:0.19,((G5|A:0.10,G6|A:0.07)94:0.12,((G7|A:0.05,G8|A:0.09)97:0.11,(G9|A:0.08,G10|A:0.05)95:0.09)96:0.14)99:0.17);\n",
    "geneD_singlecopy_discordant": "(((G1|D:0.10,G3|D:0.10)98:0.12,(G2|D:0.09,G4|D:0.09)96:0.13)97:0.20,((G5|D:0.08,G6|D:0.08)95:0.14,((G7|D:0.07,G8|D:0.07)94:0.10,(G9|D:0.06,G10|D:0.06)93:0.11)97:0.12)99:0.18);\n",
    "geneB_recent_duplication": "(((G1|B1:0.10,G2|B1:0.10)98:0.12,((G3|B1:0.04,G3|B2:0.04)99:0.04,G4|B1:0.09)97:0.09)99:0.20,((G5|B1:0.08,G6|B1:0.08)96:0.14,((G7|B1:0.07,G8|B1:0.07)95:0.10,(G9|B1:0.06,G10|B1:0.06)94:0.11)98:0.12)99:0.18);\n",
    "geneC_ancient_duplication_loss": "(((G1|C2:0.08,G2|C2:0.08)96:0.09,(G3|C2:0.07,G4|C2:0.07)95:0.10)99:0.22,(((G1|C1:0.09,G2|C1:0.09)98:0.10,(G3|C1:0.08,G4|C1:0.08)96:0.11)97:0.16,((G5|C1:0.08,G6|C1:0.08)95:0.12,((G7|C1:0.07,G8|C1:0.07)94:0.09,(G9|C1:0.06,G10|C1:0.06)93:0.10)96:0.11)98:0.15)99:0.20);\n",
    "geneP_mixed_pseudoortholog": "(((G1|P2:0.06,G3|P2:0.07)97:0.08,(G5|P2:0.06,(G7|P2:0.05,G9|P2:0.05)95:0.05)96:0.08)99:0.25,((G2|P1:0.06,G4|P1:0.07)97:0.08,(G6|P1:0.06,(G8|P1:0.05,G10|P1:0.05)95:0.05)96:0.08)99:0.25);\n",
}


def write_tutorial_dataset(outdir: str | Path) -> Path:
    outdir = Path(outdir)
    gene_dir = outdir / "gene_trees"
    gene_dir.mkdir(parents=True, exist_ok=True)
    (outdir / "species_tree.nwk").write_text(SPECIES_TREE, encoding="utf-8")

    rows = ["family\tgene_id\tspecies_id"]
    for family, newick in GENE_TREES.items():
        (gene_dir / f"{family}.nwk").write_text(newick, encoding="utf-8")
        import re

        leaves = re.findall(r"([A-Za-z0-9]+\|[A-Za-z0-9]+):", newick)
        for gene in leaves:
            species = gene.split("|", 1)[0]
            rows.append(f"{family}\t{gene}\t{species}")
    (outdir / "mapping.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")

    readme = """# Simulated 10-genome tutorial

This deterministic tree-level tutorial is designed to teach what the software measures before you analyze real data.

- `geneA_ortholog_concordant`: single-copy ortholog; topology matches the reference.
- `geneD_singlecopy_discordant`: one copy/species, but the gene tree conflicts with the species tree. Discordance alone is **not** paralogy.
- `geneB_recent_duplication`: a recent extra copy in G3.
- `geneC_ancient_duplication_loss`: an ancient duplicated family with differential copy retention.
- `geneP_mixed_pseudoortholog`: one sequence/species but inconsistent ancient copies; illustrates why one tip/species is not sufficient proof of orthology.

These are topology-level simulations, not sequence simulations.
"""
    (outdir / "README.md").write_text(readme, encoding="utf-8")
    return outdir
