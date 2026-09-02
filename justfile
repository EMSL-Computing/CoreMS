# CoreMS justfile
#
# Run `just` (or `just --list`) to see all recipes.
# Migrated from Makefile; behaviour preserved. CI-tunable knobs come
# from environment variables so `FOO=bar just recipe` mirrors the old
# `make recipe FOO=bar` UX.

lipidomics_sqlite_url := env("LIPIDOMICS_SQLITE_URL", "https://nmdcdemo.emsl.pnnl.gov/lipidomics/parameter_files/202412_lipid_ref.sqlite")
lipidomics_sqlite_path := env("LIPIDOMICS_SQLITE_PATH", "tests/tests_data/lcms/202412_lipid_ref.sqlite")
skip_lipidomics_db := env("SKIP_LIPIDOMICS_DB", "0")
skip_molecular_db := env("SKIP_MOLECULAR_DB", "0")

python := env("PYTHON", if os_family() == "windows" { "python" } else { "python3" })

# ----------------------------------------------------------------------
# Version metadata (read from .bumpversion.cfg via Python for portability).
# Evaluated lazily — only recipes that reference {{version}}/{{stage}} pay
# the cost, so recipes that don't need them still work when Python is not
# on PATH.
# ----------------------------------------------------------------------
version := if os_family() == "windows" { `python -c "import configparser; c=configparser.ConfigParser(); c.read('.bumpversion.cfg'); print(c.get('bumpversion', 'current_version', fallback='').strip())"` } else { `python3 -c "import configparser; c=configparser.ConfigParser(); c.read('.bumpversion.cfg'); print(c.get('bumpversion', 'current_version', fallback='').strip())"` }
stage := if os_family() == "windows" { `python -c "import configparser; c=configparser.ConfigParser(); c.read('.bumpversion.cfg'); print(c.get('bumpversion:part:release', 'optional_value', fallback='').strip())"` } else { `python3 -c "import configparser; c=configparser.ConfigParser(); c.read('.bumpversion.cfg'); print(c.get('bumpversion:part:release', 'optional_value', fallback='').strip())"` }

# Default: list available recipes
default:
    @just --list

# ----------------------------------------------------------------------
# Lipidomics reference DB
# ----------------------------------------------------------------------

[unix]
download-lipidomics-db:
    @if [ -f "{{lipidomics_sqlite_path}}" ]; then \
        echo "LC-MS lipidomics database already exists at {{lipidomics_sqlite_path}}"; \
    else \
        echo "Downloading LC-MS lipidomics database"; \
        mkdir -p "$(dirname "{{lipidomics_sqlite_path}}")"; \
        curl --retry 3 --retry-delay 5 --connect-timeout 30 --max-time 300 -L -o "{{lipidomics_sqlite_path}}" "{{lipidomics_sqlite_url}}"; \
        echo "LC-MS lipidomics database downloaded"; \
    fi

[windows]
download-lipidomics-db:
    @{{python}} -c "from pathlib import Path; from urllib.request import urlretrieve; p = Path(r'{{lipidomics_sqlite_path}}'); p.parent.mkdir(parents=True, exist_ok=True); url = r'{{lipidomics_sqlite_url}}'; print(f'LC-MS lipidomics database already exists at {p}') if p.exists() else (print('Downloading LC-MS lipidomics database'), urlretrieve(url, p), print('LC-MS lipidomics database downloaded'))"

# ----------------------------------------------------------------------
# Profiling helpers
# Usage: `just cpu path/to/output.prof`, `just mem path/to/script.py`
# (was `make cpu file=...` / `make mem script=...` under Make)
# ----------------------------------------------------------------------

cpu file:
    pyprof2calltree -k -i {{file}}

mem script:
    mprof run --multiprocess {{script}}
    mprof plot

# ----------------------------------------------------------------------
# Version bumps — regenerate docs after (`&& docu` runs post-body)
# ----------------------------------------------------------------------

major: && docu
    @bumpversion major --allow-dirty

minor: && docu
    @bumpversion minor --allow-dirty

patch: && docu
    @bumpversion patch --allow-dirty

# ----------------------------------------------------------------------
# PyPI publish (manual — CI publishes on tag)
# ----------------------------------------------------------------------

pypi_test:
    @rm -rf build dist *.egg-info
    @python3 -m build
    @twine upload --repository testpypi dist/*

pypi:
    @rm -rf build dist *.egg-info
    @python3 -m build
    @twine upload dist/*

# ----------------------------------------------------------------------
# Git tag from .bumpversion.cfg
# ----------------------------------------------------------------------

tag:
    @git tag -a {{version}}.{{stage}} -m "version {{version}}.{{stage}}"
    @git push origin {{version}}.{{stage}}
    @echo "tagged {{version}}.{{stage}} and pushed"

# ----------------------------------------------------------------------
# Docker image
# ----------------------------------------------------------------------

build-image-local:
    @echo corems:{{version}}
    @docker build -t corems:{{version}} .

build-image:
    @echo corilo/corems:{{version}}
    @docker build -t corilo/corems:{{version}} .

build-image-mac:
    @echo corilo/corems:{{version}}
    @docker build --platform linux/amd64 -t corilo/corems:{{version}} .

build-image-mac-local:
    @echo corems:{{version}}
    @docker build --platform linux/amd64 -t corems:{{version}} .

push-image:
    @docker push corilo/corems:{{version}}
    @docker image tag corilo/corems:{{version}} corilo/corems:latest
    @docker push corilo/corems:latest

image-run-mac:
    @docker run -it --platform linux/amd64 corilo/corems:{{version}}

image-run-mac-local:
    @docker run -it --platform linux/amd64 corems:{{version}}

image-run:
    @docker run -it corilo/corems:{{version}}

image-run-local:
    @docker run -it corems:{{version}}

# ----------------------------------------------------------------------
# Local dev database (docker-compose)
# ----------------------------------------------------------------------

db-up:
    @docker-compose up -d

db-down:
    @docker-compose down

db-logs:
    @docker-compose logs -f

db-connect:
    @docker exec -it molformdb psql -U postgres

# ----------------------------------------------------------------------
# Docs / UML
# ----------------------------------------------------------------------

uml:
    @{{python}} docs/generate_uml.py

docu: uml
    pdoc --output-dir docs --docformat numpy corems

# ----------------------------------------------------------------------
# Lint (release prep). Advisory — see RELEASE.md.
# Requires corems[dev] (pylint). Config: pyproject.toml.
# ----------------------------------------------------------------------

lint:
    {{python}} -m pylint corems

# Broader first-party Python (package + tests + support scripts). Still advisory.
lint-all:
    {{python}} -m pylint corems tests support_code

# ----------------------------------------------------------------------
# CI test recipes (install + test in one go)
# ----------------------------------------------------------------------

ci-test-source:
    #!/usr/bin/env bash
    set -euo pipefail
    {{python}} -V
    {{python}} -m pip install --upgrade pip
    # networking: networkx for molecular network plot/cluster unit tests
    {{python}} -m pip install -e ".[dev,networking]"
    skip_flags=""
    if [ "{{skip_lipidomics_db}}" = "1" ]; then
        echo "Skipping lipidomics DB download (SKIP_LIPIDOMICS_DB=1)"
        skip_flags="$skip_flags --skip-lipidomics-db"
    else
        just download-lipidomics-db
    fi
    if [ "{{skip_molecular_db}}" = "1" ]; then
        echo "Skipping molecular formula database tests (SKIP_MOLECULAR_DB=1)"
        skip_flags="$skip_flags --skip-molecular-db"
    fi
    PYTHONNET_RUNTIME=coreclr COREMS_LIPIDOMICS_SQLITE_PATH="{{lipidomics_sqlite_path}}" \
        {{python}} -m pytest --cache-clear -p no:warnings -n 4 --dist=loadfile --no-cov $skip_flags

ci-test-notebooks:
    #!/usr/bin/env bash
    set -euo pipefail
    {{python}} -V
    {{python}} -m pip install --upgrade pip
    # dev: pytest tooling; networking: networkx + ipysigma for plot/cluster/HTML
    {{python}} -m pip install -e ".[dev,networking]"
    {{python}} -m pip install --no-cache-dir jupyter nbconvert
    if [ "{{skip_lipidomics_db}}" = "1" ]; then
        echo "Skipping lipidomics DB download (SKIP_LIPIDOMICS_DB=1)"
        PYTHONNET_RUNTIME=coreclr {{python}} examples/test_notebooks.py
    else
        just download-lipidomics-db
        PYTHONNET_RUNTIME=coreclr COREMS_LIPIDOMICS_SQLITE_PATH="{{lipidomics_sqlite_path}}" \
            {{python}} examples/test_notebooks.py
    fi

ci-test-all: ci-test-source ci-test-notebooks

# ----------------------------------------------------------------------
# Local test recipes (use the current venv, do not reinstall)
# ----------------------------------------------------------------------

test-pytest-xdist: download-lipidomics-db
    pytest -n auto --no-cov --cache-clear -p no:warnings

test-notebooks: download-lipidomics-db
    python3 -m pip install --no-cache-dir jupyter nbconvert ipywidgets
    python3 -m pip install -e ".[networking]"
    cd examples && python3 test_notebooks.py

ci-test: test-pytest-xdist test-notebooks
