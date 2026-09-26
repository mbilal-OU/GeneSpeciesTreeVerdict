.PHONY: install test lint format docs tutorial

install:
	python -m pip install -e ".[dev]"

test:
	pytest --cov=genespeciestreeverdict --cov-report=term-missing

lint:
	ruff check .
	ruff format --check .

format:
	ruff check --fix .
	ruff format .

docs:
	mkdocs build --strict

tutorial:
	gstv tutorial --outdir gstv_tutorial
	gstv batch --species-tree gstv_tutorial/species_tree.nwk --gene-tree-dir gstv_tutorial/gene_trees --mapping gstv_tutorial/mapping.tsv --outdir gstv_tutorial/results
