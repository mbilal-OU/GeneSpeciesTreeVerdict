# Real-data workflow

## Recommended workflow

```text
reference species tree
        +
per-family gene trees
        +
gene→species mapping
        ↓
GeneSpeciesTreeVerdict
        ↓
coverage / copy number / topology / DL diagnostic
        ↓
PASS / REVIEW / RESOLVE + reasons
        ↓
manual or method-specific follow-up
```

## Step 1 — prepare the species tree
Use a Newick tree with unique terminal labels. If rooted duplication-loss diagnostics are interpreted, the root must be biologically meaningful.

## Step 2 — prepare gene trees
Each family should be in its own Newick-like file (`.nwk`, `.newick`, `.tree`, `.treefile`, `.tre`). Gene-copy labels must be unique within a tree.

## Step 3 — map gene copies to species
For real data, explicit mappings are safer than guessing from sequence names. Use either a two-column per-family table or a three-column combined table.

## Step 4 — run batch analysis

```bash
gstv batch \
  --species-tree species_tree.nwk \
  --gene-tree-dir gene_trees/ \
  --mapping mappings.tsv \
  --min-coverage 0.95 \
  --outdir verdict_results
```

## Step 5 — interpret by evidence class
- `PASS`: candidate conventional single-copy marker.
- `REVIEW`: inspect conflict, support, rooting and alternative processes.
- `RESOLVE`: multi-copy/duplication evidence means orthology resolution is needed before ordinary concatenation.

## Step 6 — do not stop at the verdict
For reviewed/resolved loci, use appropriate follow-up analyses such as GeneRax/AleRax/Treerecs/Notung, ASTRAL-Pro/ASTER, synteny/context checks, sequence re-alignment, alternative rooting or additional taxon sampling.

## Step 7 — compare marker-set sensitivity
A strong study should compare downstream species-tree results before and after marker QC rather than assuming filtering is harmless. Record topology changes, branch support, alignment length and loci removed.
