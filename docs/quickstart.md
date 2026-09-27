# Quick start

## Install

GeneSpeciesTreeVerdict requires Python 3.10 or newer and is currently tested on Python 3.10–3.13. Check your interpreter before creating the environment:

```bash
python --version
```

Then create an isolated environment and install the development dependencies:

```bash
git clone https://github.com/mbilal-OU/GeneSpeciesTreeVerdict.git
cd GeneSpeciesTreeVerdict
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e ".[dev]"
```

Windows PowerShell activation:

```powershell
.venv\Scripts\Activate.ps1
```

## Learn with the bundled simulation

```bash
gstv tutorial --outdir tutorial_demo

gstv batch \
  --species-tree tutorial_demo/species_tree.nwk \
  --gene-tree-dir tutorial_demo/gene_trees \
  --mapping tutorial_demo/mapping.tsv \
  --outdir tutorial_demo/results
```

## Analyze one real locus

```bash
gstv analyze \
  --species-tree species_tree.nwk \
  --gene-tree gene001.treefile \
  --mapping mapping.tsv \
  --outdir gene001_verdict
```

## Analyze a directory

```bash
gstv batch \
  --species-tree species_tree.nwk \
  --gene-tree-dir gene_trees/ \
  --mapping mappings.tsv \
  --min-coverage 0.95 \
  --outdir verdict_results/
```

Always inspect the reasons and warnings in addition to the headline verdict.
