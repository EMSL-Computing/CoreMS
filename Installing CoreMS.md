# Installing CoreMS

These instructions are for **CoreMS 4.0 and later**. They cover a normal pip install, optional PostgreSQL, and Thermo `.raw` support.

For a short overview, see the [Installation](./README.md#installation) section in the README. For contributing and local tests, see [CONTRIBUTING.md](./CONTRIBUTING.md).

## Requirements

| Requirement | Notes |
|---|---|
| **Python 3.10–3.13** | Specified as `requires-python >=3.10` in `pyproject.toml` |
| **pip** | Upgrade recommended: `python -m pip install -U pip` |
| **Optional: Docker** | PostgreSQL for molecular formula search (`docker-compose`), or the full CoreMS image |
| **Optional: .NET 8 runtime** | Needed on the host for Thermo `.raw` files via Python.NET (see [Thermo](#thermo-raw-file-support)) |

CoreMS depends on **pythonnet** (`pythonnet>=3.0.3`) and other packages declared in `pyproject.toml`. A normal `pip install` pulls them in; you do **not** need a separate `pip install pythonnet` step.

## Quick start (PyPI)

```bash
python -m venv venv
# Windows: venv\Scripts\activate
source venv/bin/activate

python -m pip install -U pip
pip install corems
```

Development and test extras:

```bash
pip install "corems[dev]"
```

Verify:

```bash
python -c "import corems; print(corems.__version__)"
```

By default the molecular formula database uses **SQLite** (local file; simplest setup).

## Install from source

```bash
git clone https://github.com/EMSL-Computing/CoreMS.git
cd CoreMS

python -m venv venv
source venv/bin/activate   # Windows: venv\Scripts\activate

python -m pip install -U pip
pip install .
# or, for development:
# pip install -e ".[dev]"
```

Internal GitLab clone (same install commands after clone):

```text
https://code.emsl.pnl.gov/mass-spectrometry/corems.git
```

## Molecular formula database

### SQLite (default)

No extra services. Suitable for most users and light workloads.

### PostgreSQL (optional)

For better multi-user / higher-load molecular formula search performance:

1. Install and start [Docker](https://www.docker.com/).
2. From the CoreMS repository root:

```bash
docker-compose up -d
```

3. Point CoreMS at the database (either approach works):

- Set the environment variable:

```bash
export COREMS_DATABASE_URL="postgresql+psycopg2://coremsappdb:coremsapppnnl@localhost:5432/coremsapp"
```

- Or set `MSParameters.molecular_search.url_database` to the same URL in code / config.

Keep Docker running while you use this database.

## Thermo `.raw` file support

Reading Thermo `.raw` files uses **Python.NET** (`pythonnet`), which is installed with CoreMS. The .NET **runtime** on the machine must match how Python.NET is configured.

### Recommended: .NET 8 (CoreCLR)

This matches CoreMS CI and the project Dockerfile (**.NET 8**, not .NET 9).

1. Install the **.NET 8 runtime** for your OS and CPU (x64 or arm64), for example from Microsoft’s [.NET 8 download page](https://dotnet.microsoft.com/download/dotnet/8.0).
2. Ensure `dotnet` is on your `PATH`, and set:

```bash
export PYTHONNET_RUNTIME=coreclr
# If the runtime is not on the default path, also set e.g.:
# export DOTNET_ROOT=/usr/local/dotnet
# export PATH="$PATH:$DOTNET_ROOT"
```

3. Confirm:

```bash
dotnet --list-runtimes
# should list a Microsoft.NETCore.App 8.x runtime
```

On **Apple Silicon** (M1/M2/M3/…), install the **arm64** .NET 8 runtime when running Python natively on arm64. If you only use the CoreMS Docker image, host .NET is not required (the image already includes .NET 8).

### Optional: Mono (macOS / some Linux setups)

The preferred path is CoreCLR + .NET 8. Some users still use **Mono** with Python.NET instead.

- Mono project / install options: [https://www.mono-project.com/](https://www.mono-project.com/) (e.g. Homebrew `mono` on macOS is a common source, but is **not** required for CoreMS 4.0+).
- If you use Mono, point Python.NET at it explicitly:

```bash
export PYTHONNET_RUNTIME=mono
```

You may also need Mono’s libraries discoverable on your system (distribution-specific). Prefer CoreCLR unless you already rely on Mono.

### Windows notes

- `pip install corems` installs pythonnet with the package.
- Install **.NET 8** if you need Thermo `.raw` support, and use `PYTHONNET_RUNTIME=coreclr` when required by your environment.

## Docker image (full application)

To run CoreMS in a container that already includes **.NET 8** and dependencies, see [Building and Running the CoreMS Docker Image](./README.md#docker-image) in the README (`make build-image-local` / `make build-image-mac-local` on Apple Silicon).

## Example notebooks

Workflow examples live under `examples/notebooks/`. Use any environment that can run Jupyter if you want to execute them interactively; installing an IDE is optional and not covered here.

## Troubleshooting

| Symptom | What to check |
|---|---|
| `import corems` fails / wrong Python | Active venv; `python --version` is 3.10+ |
| Thermo `.raw` import / pythonnet errors | .NET **8** runtime installed; `PYTHONNET_RUNTIME=coreclr`; `dotnet --list-runtimes` |
| Using Mono deliberately | `PYTHONNET_RUNTIME=mono`; Mono actually installed |
| PostgreSQL connection errors | `docker-compose up -d`; `COREMS_DATABASE_URL` / URL in parameters |
| Apple Silicon Docker builds | Use the macOS image build targets in the README (`linux/amd64` platform) |

## See also

- [README.md](./README.md) — overview, short install, Docker image
- [CONTRIBUTING.md](./CONTRIBUTING.md) — development workflow and local tests
- [RELEASE.md](./RELEASE.md) — versioning and releases (maintainers)
