# Release checklist

Use this checklist before every tagged release.

- [ ] All supported Python versions pass CI.
- [ ] `ruff check`, formatting, and `pytest` are clean.
- [ ] Documentation builds with `mkdocs build --strict`.
- [ ] Tutorial scenarios are rerun and expected verdicts verified.
- [ ] CLI help and output examples match the current version.
- [ ] `CHANGELOG.md` is updated.
- [ ] `CITATION.cff` and package version are updated.
- [ ] Scientific assumptions and changed decision rules are documented.
- [ ] Breaking changes are explicitly identified.
- [ ] Release notes include limitations and known issues.
- [ ] Tag uses semantic versioning.
- [ ] After a stable public release, archive the tagged version with a DOI service such as Zenodo and update citation metadata.

A release is not just a packaging event: for scientific software, it is a frozen statement of algorithms, defaults, assumptions, and expected behavior.
