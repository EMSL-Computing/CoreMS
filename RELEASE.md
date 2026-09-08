# Releasing CoreMS

Maintainers only. Contributors merge work into `dev`; releases are cut by merging `dev` into `master` on **GitLab** (source of truth; mirrored to GitHub).

## Branch model

| Branch | Role |
|---|---|
| `dev` | Integration. All feature PRs/MRs land here. |
| `master` | Stable / released code. Updated only by a release MR. |

```
PRs/MRs  -->  dev  --(release MR)-->  master  -->  tag  -->  CI publishes to PyPI
```

## When to release

Release when `dev` has a coherent set of changes ready for users (features, fixes, or both). Coordinate version bump size with the change set:

| Bump | Use when |
|---|---|
| `make patch` | Backwards-compatible bug fixes |
| `make minor` | Backwards-compatible new features |
| `make major` | Incompatible API changes |

Each of those updates version metadata (see `.bumpversion.cfg`) and regenerates docs via `make docu` (UML class diagrams via `make uml` / pyreverse + Graphviz, then pdoc). The install how-to source is `docs/user/installation.md`, included on the package landing page through `corems.__doc__`. Maintainers need the `dev` extra (`ruff` for `make lint`, `pylint` for `pyreverse` / `make uml`, `pdoc`, …) and system Graphviz (`dot`) so format and docs steps succeed during a version bump.

## Release steps (GitLab)

All release steps should be done on the `dev` branch and then merged into `master` via a merge request (MR). Only bump the version on `dev` and merge into `master`. 

1. **Ensure `dev` is ready**
   - Gitlab CI green on `dev`.
   - Changelog or release notes drafted (these will be copied into the MR description and later into the release on GitHub).
   - No open blockers for the intended version.
   - **Format and auto-fix the package** (release prep; not a hard CI gate). This **rewrites files**.
     1. Activate the same environment you use for CoreMS development and install/update dev extras so `ruff` is available, e.g.:
        ```bash
        pip install -e ".[dev]"
        ```
     2. From the repo root run:
        ```bash
        make lint
        ```
        This runs `ruff format` then `ruff check --fix` on `corems` (config in `pyproject.toml`). Use the venv’s Python (activate the venv first, or `make lint PYTHON=.venv/bin/python`). The first run after ruff is introduced will produce a large diff; later releases should be small.
     3. Optional broader pass (tests and support scripts; also rewrites):
        ```bash
        make lint-all
        ```
     4. Review `git diff` and commit the formatted/fixed files (a dedicated commit is fine, especially the first time). Findings that `--fix` cannot apply may remain and may make the command exit non-zero; that does not block the release. Do not use `--unsafe-fixes`.

2. **Bump version on `dev` (or a short-lived release branch from `dev`)**
   ```bash
   git checkout dev
   git pull
   make patch   # or: make minor / make major
   git add -u   # or: git add -A if you want to include new files AND your repo is clean
   git commit -m "Bump version for release x.y.z"
   git push origin dev
   ```

3. **Open a release MR on GitLab**
   - Source: `dev`
   - Target: `master`
   - Title: "Bump version for release x.y.z" 
   - Description: Release notes / changelog (this will be copied into the GitHub release later!)
   - Require CI green and maintainer approval

4. **Merge the MR into `master`**

5. **Tag the release**
   ```bash
   git checkout master
   git pull
   make tag
   ```
   `make tag` creates an annotated tag from `.bumpversion.cfg` and pushes it.

6. **Publish via CI**
   - After the tag is on `master`, CI/CD publishes the package to PyPI. Do not run `make pypi` by hand.
   - Wait and verify that the release is visible on [PyPI](https://pypi.org/project/corems/) and CI is green on the tag pipeline.

7. **Sync `dev`**
   If `master` has commits that are not yet on `dev` (for example a version bump or a fix that landed only on `master`), update `dev` locally and push:
   ```bash
   git checkout dev
   git pull
   git merge master
   git push origin dev
   ```
8. **Create a GitHub release**
   - Wait until the GitLab → GitHub mirror has the new tag (and `master`).
   - On GitHub: **Releases** → **Draft a new release**.
   - Choose the existing tag (do not create a new one).
   - Paste the release notes from the GitLab release MR description into the release body.
   - Publish the release.

