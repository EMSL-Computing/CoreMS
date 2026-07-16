# Contributing to CoreMS

Thanks for contributing. This guide covers expectations for external and internal contributors, how changes land, and how to test locally.

For cutting a release, see [RELEASE.md](./RELEASE.md).

## Table of Contents

- [Contributor paths](#contributor-paths)
- [Branch model](#branch-model)
- [Getting started](#getting-started)
- [Local testing](#local-testing)
- [PR / MR checklist](#pr--mr-checklist)
- [Code style](#code-style)
- [Issue reporting](#issue-reporting)
- [License](#license)

## Contributor paths

CoreMS is developed primarily on **GitLab** and mirrored to **GitHub**.

| | External (public) | Internal (EMSL / PNNL) |
|---|---|---|
| Where you work | [GitHub](https://github.com/EMSL-Computing/CoreMS) | GitLab (`code.emsl.pnl.gov`, mass-spectrometry/corems) |
| How you contribute | Fork, then open a **pull request** | Branch on GitLab, then open a **merge request** |
| Target branch | `dev` only | `dev` only |
| Review | Maintainer review on GitHub | Maintainer review on GitLab |
| What not to do | Do not open PRs against `master` | Do not open MRs against `master` |

GitHub is the public mirror of the GitLab project. Prefer GitLab for internal work so CI and review stay on the source of truth.

## Branch model

```
feature / fix branch  -->  dev  -->  master (releases only)
```

- **`dev`**: integration branch. All PRs and MRs target `dev`.
- **`master`**: release branch. Only maintainers merge `dev` into `master` when cutting a release. See [RELEASE.md](./RELEASE.md).

Do not open feature work against `master`.

## Getting started

1. Open an issue describing the bug or feature (unless one already exists).
2. Install dependencies. See [README.md](./README.md) and [Installing CoreMS.md](./Installing%20CoreMS.md).
3. Create a branch (or fork, for external contributors).
4. Make your changes. Add or update tests and docs as needed.
5. Run local tests (below).
6. Open a PR (GitHub) or MR (GitLab) **into `dev`**. Reference the issue (e.g. `closes #23`).
7. Address review feedback until CI is green and a maintainer approves.

Version bumps and packaging are handled at release time by maintainers, not on every feature PR/MR.

## Local testing

Activate your existing virtualenv first (with CoreMS, pytest, and Thermo/.NET support already set up). Then run tests from the repo root **without** reinstalling the package.

Prefer these targets for day-to-day work. They use your current env and do not run `pip install`:

| Command | What it runs |
|---|---|
| `make test-pytest-xdist` | pytest with xdist (recommended default) |
| `make test-notebooks` | Example notebook tests |
| `make ci-test` | Both of the above |
| `make download-lipidomics-db` | LC-MS lipidomics SQLite used by some tests (skipped if already present) |

Before opening a PR or MR:

```bash
make ci-test
# or source tests only:
make test-pytest-xdist
```

## PR / MR checklist

Reviewers check these before merging into **`dev`**:

1. CI is green (tests pass).
2. Unit tests cover new or changed behavior.
3. Docs and docstrings are updated when the public API or user-facing behavior changes.
4. Related issues/PRs/MRs are referenced.
5. Merge request title and/or description states if the PR/MR are a new feature, bug fix, or other change.
6. Target branch is `dev` (not `master`).

## Code style

- Docstrings follow the [NumPy style](https://numpydoc.readthedocs.io/en/latest/format.html).
- API docs are built with [pdoc](https://github.com/mitmproxy/pdoc) (`make docu`).

## Issue reporting

Report bugs and feature requests in the issue tracker on the platform you use (GitHub or GitLab). Include steps to reproduce, expected vs actual behavior, and relevant logs or versions.

## License

By contributing, you agree that your contributions are licensed as described in [LICENSE](./LICENSE).
