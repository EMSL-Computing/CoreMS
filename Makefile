app_name = CoreMS
parameters_path = parameter.json

LIPIDOMICS_SQLITE_URL ?= https://nmdcdemo.emsl.pnnl.gov/lipidomics/parameter_files/202412_lipid_ref.sqlite
LIPIDOMICS_SQLITE_PATH ?= tests/tests_data/lcms/202412_lipid_ref.sqlite

.PHONY: download-lipidomics-db ci-test-source ci-test-notebooks ci-test-all test-pytest-xdist test-notebooks ci-test

# ----------------------------------------------------------------------
# Platform-specific logic
#
# Default path: macOS/Linux/CI use POSIX shell tools.
# Windows path: avoid assuming Bash, WSL, PowerShell, curl, grep, awk,
# dirname, or mkdir -p. Use Python instead.
# ----------------------------------------------------------------------

ifeq ($(OS),Windows_NT)

# Windows: Python fallback
PYTHON ?= python

version := $(shell $(PYTHON) -c "import configparser; c=configparser.ConfigParser(); c.read('.bumpversion.cfg'); print(c['bumpversion']['current_version'].strip())")
stage := $(shell $(PYTHON) -c "import configparser; c=configparser.ConfigParser(); c.read('.bumpversion.cfg'); print(c['bumpversion:part:release']['optional_value'].strip())")

download-lipidomics-db:
	@$(PYTHON) -c "from pathlib import Path; from urllib.request import urlretrieve; p = Path(r'$(LIPIDOMICS_SQLITE_PATH)'); p.parent.mkdir(parents=True, exist_ok=True); url = r'$(LIPIDOMICS_SQLITE_URL)'; print(f'LC-MS lipidomics database already exists at {p}') if p.exists() else (print('Downloading LC-MS lipidomics database'), urlretrieve(url, p), print('LC-MS lipidomics database downloaded'))"

else

# macOS/Linux/CI: POSIX shell path
PYTHON ?= python3

version := $(shell grep '^[[:space:]]*current_version[[:space:]]*=' .bumpversion.cfg | head -n 1 | cut -d= -f2 | tr -d ' ')
stage := $(shell awk 'BEGIN{section=0} /^\[bumpversion:part:release\]/{section=1; next} /^\[/{section=0} section && /^[[:space:]]*optional_value[[:space:]]*=/{split($$0,a,"="); gsub(/[[:space:]]/,"",a[2]); print a[2]; exit}' .bumpversion.cfg)

download-lipidomics-db:
	@if [ -f "$(LIPIDOMICS_SQLITE_PATH)" ]; then \
		echo "LC-MS lipidomics database already exists at $(LIPIDOMICS_SQLITE_PATH)"; \
	else \
		echo "Downloading LC-MS lipidomics database"; \
		mkdir -p "$$(dirname "$(LIPIDOMICS_SQLITE_PATH)")"; \
		curl --retry 3 --retry-delay 5 --connect-timeout 30 --max-time 300 -L -o "$(LIPIDOMICS_SQLITE_PATH)" "$(LIPIDOMICS_SQLITE_URL)"; \
		echo "LC-MS lipidomics database downloaded"; \
	fi

endif


cpu:
	pyprof2calltree -k -i $(file)

mem: 

	mprof run --multiprocess $(script)
	mprof plot

major:
	
	@bumpversion major --allow-dirty
	@$(MAKE) docu

minor:
	
	@bumpversion minor --allow-dirty
	@$(MAKE) docu

patch:
	
	@bumpversion patch --allow-dirty
	@$(MAKE) docu

pypi_test:
	@rm -rf build dist *.egg-info
	@python3 -m build
	@twine upload --repository testpypi dist/*

pypi:
	@rm -rf build dist *.egg-info
	@python3 -m build
	@twine upload dist/*

tag:

	@git tag -a $(version).$(stage) -m "version $(version).$(stage)"
	@git push origin $(version).$(stage)
	@echo tagged $(version).$(stage) and pushed

build-image-local:

	@echo corems:$(version)
	@docker build -t corems:$(version) .

build-image:

	@echo corilo/corems:$(version)
	@docker build -t corilo/corems:$(version) .

build-image-mac:

	@echo corilo/corems:$(version)
	@docker build --platform linux/amd64 -t corilo/corems:$(version) .

build-image-mac-local:

	@echo corems:$(version)
	@docker build --platform linux/amd64 -t corems:$(version) .

push-image:

	@docker push corilo/corems:$(version)
	@docker image tag corilo/corems:$(version) corilo/corems:latest
	@docker push corilo/corems:latest

image-run-mac:

	@docker run -it --platform linux/amd64 corilo/corems:$(version)

image-run-mac-local:

	@docker run -it --platform linux/amd64 corems:$(version)

image-run:

	@docker run -it corilo/corems:$(version)

image-run-local:

	@docker run -it corems:$(version)

db-up:

	@docker-compose up -d 

db-down:

	@docker-compose down

db-logs:

	@docker-compose logs -f

db-connect:

	@docker exec -it molformdb psql -U postgres

docu:
	
	pdoc --output-dir docs --docformat numpy corems

SKIP_LIPIDOMICS_DB ?= 0
SKIP_MOLECULAR_DB ?= 0

ci-test-source:
	@$(PYTHON) -V
	@$(PYTHON) -m pip install --upgrade pip
	# networking: networkx for molecular network plot/cluster unit tests
	@$(PYTHON) -m pip install -e ".[dev,networking]"
	@skip_flags=""; \
	if [ "$(SKIP_LIPIDOMICS_DB)" = "1" ]; then \
		echo "Skipping lipidomics DB download (SKIP_LIPIDOMICS_DB=1)"; \
		skip_flags="$$skip_flags --skip-lipidomics-db"; \
	else \
		$(MAKE) download-lipidomics-db LIPIDOMICS_SQLITE_PATH="$(LIPIDOMICS_SQLITE_PATH)"; \
	fi; \
	if [ "$(SKIP_MOLECULAR_DB)" = "1" ]; then \
		echo "Skipping molecular formula database tests (SKIP_MOLECULAR_DB=1)"; \
		skip_flags="$$skip_flags --skip-molecular-db"; \
	fi; \
	PYTHONNET_RUNTIME=coreclr COREMS_LIPIDOMICS_SQLITE_PATH="$(LIPIDOMICS_SQLITE_PATH)" \
		$(PYTHON) -m pytest --cache-clear -p no:warnings -n 4 --dist=loadfile --no-cov $$skip_flags

ci-test-notebooks:
	@$(PYTHON) -V
	@$(PYTHON) -m pip install --upgrade pip
	# dev: pytest tooling; networking: networkx + ipysigma for plot/cluster/HTML
	@$(PYTHON) -m pip install -e ".[dev,networking]"
	@$(PYTHON) -m pip install --no-cache-dir jupyter nbconvert
	@if [ "$(SKIP_LIPIDOMICS_DB)" = "1" ]; then \
		echo "Skipping lipidomics DB download (SKIP_LIPIDOMICS_DB=1)"; \
		PYTHONNET_RUNTIME=coreclr $(PYTHON) examples/test_notebooks.py; \
	else \
		$(MAKE) download-lipidomics-db LIPIDOMICS_SQLITE_PATH="$(LIPIDOMICS_SQLITE_PATH)"; \
		PYTHONNET_RUNTIME=coreclr COREMS_LIPIDOMICS_SQLITE_PATH="$(LIPIDOMICS_SQLITE_PATH)" \
			$(PYTHON) examples/test_notebooks.py; \
	fi

ci-test-all: ci-test-source ci-test-notebooks

test-pytest-xdist: download-lipidomics-db
	pytest -n auto --no-cov --cache-clear -p no:warnings

test-notebooks: download-lipidomics-db
	python3 -m pip install --no-cache-dir jupyter nbconvert ipywidgets
	python3 -m pip install -e ".[networking]"
	cd examples && python3 test_notebooks.py

ci-test: test-pytest-xdist test-notebooks
