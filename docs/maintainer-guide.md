# Maintainer guide

This repository is designed to model good research-software maintenance, not only produce phylogenetic output.

## Branch model

- `main` should remain releasable.
- Develop nontrivial changes on feature branches.
- Merge through pull requests.
- Prefer squash merge for small focused PRs or regular merge when preserving a meaningful multi-commit history matters.

Recommended branch protection for `main`:
- require a pull request before merging;
- require CI checks;
- require branches to be up to date when practical;
- block force pushes and deletion;
- require conversation resolution.

## Issue triage

Use issues to distinguish:
- software bugs;
- scientific-correctness concerns;
- feature requests;
- documentation problems;
- benchmark/validation tasks.

A scientific-correctness issue should receive higher priority than cosmetic work because incorrect biological interpretation is a functional defect.

## Pull-request review checklist

Before merging:
1. Does the change alter a scientific assumption?
2. Is that assumption documented and cited?
3. Are edge cases tested?
4. Does the simulated tutorial still behave as expected?
5. Are CLI/output changes documented?
6. Does `make test` pass?
7. Does `make lint` pass?
8. Does `mkdocs build --strict` pass?
9. Does the package build successfully?
10. Is the changelog updated when user-visible behavior changes?

## Versioning

Use semantic versioning once releases stabilize:
- patch: bug/documentation fixes without intended API break;
- minor: backward-compatible features;
- major: incompatible CLI/output/API changes.

During alpha, document breaking changes explicitly even if the major version remains zero.

## Release checklist

1. confirm CI on supported Python versions;
2. run the tutorial end-to-end;
3. run package build checks;
4. update `CHANGELOG.md`;
5. update version in `pyproject.toml`, package `__version__`, and `CITATION.cff`;
6. create a signed/annotated tag when possible;
7. create a GitHub release with release notes;
8. archive the stable release to Zenodo once configured;
9. publish to PyPI/Bioconda only after packaging is stable;
10. record benchmark data and exact tool versions for manuscript releases.

## Dependencies

Dependabot is configured for Python and GitHub Actions. Review dependency PRs rather than blindly auto-merging scientific-stack changes. Run the full tutorial and tests after Biopython changes because tree parsing/behavior is scientifically important.

## Documentation discipline

Any new verdict rule must answer:
- What evidence triggers it?
- What does it mean?
- What does it **not** mean?
- What alternative evolutionary explanations remain?
- Which downstream methods are appropriate?

## Reproducibility

For publications, record:
- GeneSpeciesTreeVerdict version and commit SHA;
- Python version;
- dependency versions;
- species-tree file checksum;
- gene-tree/mapping provenance;
- CLI command and non-default thresholds;
- downstream phylogenetic software versions and models.

## Repository health

Regularly review:
- failing CI;
- stale issues/PRs;
- dependency/security alerts;
- documentation drift;
- unsupported Python versions;
- benchmark regressions;
- whether roadmap claims still match implemented functionality.
