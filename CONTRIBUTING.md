# Contributing

Contributions to this adaptation are welcome through issues and pull requests on [Eamonmca/GuideHarbour](https://github.com/Eamonmca/GuideHarbour). Upstream ProbeHarbour provenance is recorded in `UPSTREAM.md`.

1. Create a focused branch from `main`.
2. Install the development environment with `python -m pip install -e '.[test,dev]'`.
3. Add or update tests for behavioural changes.
4. If an nbdev source notebook changes, export it and include the generated Python change.
5. Run `pytest`, `python -m build`, and `twine check dist/*` before opening a pull request.

Do not commit patient data, vendor-confidential exports, internal calibration panels, genome references, BLAST databases, generated ProbeStores or machine-specific paths. Use small synthetic or redistributable fixtures in public tests.

Please report security-sensitive issues privately to the maintainers rather than in a public issue.
