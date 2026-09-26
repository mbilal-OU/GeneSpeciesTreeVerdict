# Contributing

Thank you for helping improve GeneSpeciesTreeVerdict.

## Before coding
For nontrivial changes, open an issue first so the scientific assumption and intended behavior can be discussed before implementation.

## Development setup
```bash
git clone https://github.com/mbilal-OU/GeneSpeciesTreeVerdict.git
cd GeneSpeciesTreeVerdict
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e ".[dev]"
```

## Quality checks
```bash
make test
make lint
make docs
```

## Scientific contribution standard
Changes to reconciliation, topology metrics or verdict rules must include:
- the biological/computational assumption;
- a citation when the method is established in the literature;
- automated tests;
- a minimal or simulated example when behavior is non-obvious;
- documented limitations and alternative explanations.

## Pull requests
Keep each PR focused. Describe the problem, method, tests, user-visible changes and any change in scientific interpretation.

## Bug reports
Provide the smallest species tree, gene tree, mapping and exact command that reproduce the behavior whenever sharing those data is permissible.
