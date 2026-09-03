# Installing CoreMS

How to install **CoreMS 4.0+** for normal use, development, optional databases, and Thermo `.raw` support.

This is the **canonical install guide**. The [README](../../README.md) has a short quickstart only. This file is loaded into `corems.__doc__` (together with the README) so `make docu` / pdoc shows it on the package landing page with the same styling as the API docs.

## Requirements

| Requirement | Notes |
|---|---|
| **Python 3.10–3.13** | From `requires-python` in `pyproject.toml` |
| **pip** | Prefer upgrading: `python -m pip install -U pip` |
| **Optional: Docker** | PostgreSQL via `docker-compose`, or the full CoreMS image |
| **Optional: .NET 8 runtime** | Host Thermo `.raw` via Python.NET (see [system dependencies](#system-dependencies-thermo-raw)) |

Declared Python dependencies (including **pythonnet**) are installed with the package. You do **not** need a separate `pip install pythonnet`.

## Quick install (PyPI)

Use an isolated environment (venv or conda). Example with venv:

```bash
python -m venv venv
# Windows: venv\Scripts\activate
source venv/bin/activate

python -m pip install -U pip
pip install corems
```

Verify:

```bash
python -c "import corems; print(corems.__version__)"
```

By default the molecular formula database uses **SQLite** (local file; no extra services).

## Optional extras

Install extras from `pyproject.toml` as needed:

```bash
pip install "corems[dev]"          # tests, pylint, docs tooling (incl. pdoc)
```

With `corems[dev]`, maintainers can run package lint via `make lint` (see [RELEASE.md](../../RELEASE.md)).

Check the installed package metadata or `pyproject.toml` for the current list of extras.

## Install from source

```bash
git clone https://github.com/EMSL-Computing/CoreMS.git
cd CoreMS

python -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate

python -m pip install -U pip
pip install .
# editable install for development:
# pip install -e ".[dev]"
```

Internal GitLab (same install after clone):

```text
https://code.emsl.pnl.gov/mass-spectrometry/corems.git
```

## Molecular formula database

### SQLite (default)

No extra services. Suitable for most users and light workloads.

### PostgreSQL (optional)

For multi-user or higher-load molecular formula search:

1. Install and start [Docker](https://www.docker.com/).
2. From the CoreMS repository root:

```bash
docker-compose up -d
```

3. Point CoreMS at the database (either approach):

```bash
export COREMS_DATABASE_URL="postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp"
```

Or set `MSParameters.molecular_search.url_database` to the same URL.

Keep Docker running while you use this database.

## System dependencies (Thermo `.raw`)

Reading Thermo `.raw` files uses **Python.NET** (a package dependency) plus a **.NET runtime on the host**. This is separate from `pip install`.

### Recommended: .NET 8 (CoreCLR)

Matches CoreMS CI and the project Dockerfile (**.NET 8**, not .NET 9).

1. Install the **.NET 8 runtime** for your OS and CPU (x64 or arm64): [Download .NET 8](https://dotnet.microsoft.com/download/dotnet/8.0).
2. Prefer CoreCLR:

```bash
export PYTHONNET_RUNTIME=coreclr
# If the runtime is not on the default path:
# export DOTNET_ROOT=/usr/local/dotnet
# export PATH="$PATH:$DOTNET_ROOT"
```

3. Confirm:

```bash
dotnet --list-runtimes
# should list a Microsoft.NETCore.App 8.x runtime
```

On **Apple Silicon**, install the **arm64** .NET 8 runtime when running Python natively on arm64. The CoreMS **Docker image** already includes .NET 8, so host .NET is optional if you only run inside the image.

### Optional: Mono

Some users still use [Mono](https://www.mono-project.com/) with Python.NET. That path is **not** required for CoreMS 4.0+. If you use Mono:

```bash
export PYTHONNET_RUNTIME=mono
```

You may also need Mono’s libraries discoverable on your system (distribution-specific). Prefer CoreCLR unless you already depend on Mono.

### Windows

- `pip install corems` installs pythonnet with the package.
- Install **.NET 8** if you need Thermo `.raw` support; use `PYTHONNET_RUNTIME=coreclr` when required.

## Docker image (full application)

To run CoreMS in a container that already includes **.NET 8** and dependencies, see [Building and Running the CoreMS Docker Image](../../README.md#docker-image) in the README (`make build-image-local` / `make build-image-mac-local` on Apple Silicon).

## Example notebooks

Workflow examples live under `examples/notebooks/`. These are intended as demonstrations and starting points for your own analysis. 

## Troubleshooting

| Symptom | What to check |
|---|---|
| `import corems` fails / wrong Python | Active venv; `python --version` is 3.10+ |
| Thermo `.raw` / pythonnet errors | .NET **8** installed; `PYTHONNET_RUNTIME=coreclr`; `dotnet --list-runtimes` |
| Using Mono deliberately | `PYTHONNET_RUNTIME=mono`; Mono installed |
| PostgreSQL connection errors | `docker-compose up -d`; `COREMS_DATABASE_URL` / URL in parameters |
| Apple Silicon Docker builds | macOS image build targets in the README (`linux/amd64`) |

## Building the documentation site

```bash
make docu
```

Opens as `docs/corems.html`. The landing module page includes this install guide via `corems/__init__.py`.

## See also

- [README.md](../../README.md) — overview, short install, Docker image
- [CONTRIBUTING.md](../../CONTRIBUTING.md) — development workflow and local tests
- [RELEASE.md](../../RELEASE.md) — versioning and releases (maintainers)
- [Package / API docs](https://emsl-computing.github.io/CoreMS/corems.html) (generated)
