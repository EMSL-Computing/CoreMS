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

Each of those updates version metadata (see `.bumpversion.cfg`) and regenerates docs via `make docu`.

## Release steps (GitLab)

Do this on GitLab, not as an ad hoc push to `master` on GitHub.

1. **Ensure `dev` is ready**
   - CI green on `dev`.
   - Changelog or release notes drafted (these will be copied into the MR description and later into the release on GitHub).
   - No open blockers for the intended version.

2. **Bump version on `dev` (or a short-lived release branch from `dev`)**
   ```bash
   git checkout dev
   git pull
   make patch   # or: make minor / make major
   git add -u   # or: git add -A if you want to include new files and your repo is clean
   git commit -m "Bump version for release"
   git push origin dev
   ```

3. **Open a release MR on GitLab**
   - Source: `dev`
   - Target: `master`
   - Title: version 
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

