from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from Bio import Phylo

from .io import read_tree


def plot_tree_comparison(
    species_tree_path: str | Path,
    gene_tree_path: str | Path,
    output_path: str | Path,
    *,
    title: str | None = None,
) -> Path:
    species_tree = read_tree(species_tree_path)
    gene_tree = read_tree(gene_tree_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(1, 2, figsize=(14, 7))
    Phylo.draw(species_tree, axes=axes[0], do_show=False, show_confidence=True)
    axes[0].set_title("Reference species tree")
    axes[0].set_ylabel("")
    Phylo.draw(gene_tree, axes=axes[1], do_show=False, show_confidence=True)
    axes[1].set_title("Gene tree")
    axes[1].set_ylabel("")
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return output_path
