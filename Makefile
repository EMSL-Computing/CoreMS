app_name = CoreMS
parameters_path = parameter.json 
version := $(shell cat .bumpversion.cfg | grep current_version | cut -d= -f2 | tr -d ' ')
stage := $(shell cat .bumpversion.cfg | grep optional_value | cut -d= -f2 | tr -d ' ') 
LIPIDOMICS_SQLITE_URL ?= https://nmdcdemo.emsl.pnnl.gov/lipidomics/parameter_files/202412_lipid_ref.sqlite
LIPIDOMICS_SQLITE_PATH ?= tests/tests_data/lcms/202412_lipid_ref.sqlite
PYTHON ?= python3

.PHONY: download-lipidomics-db ci-test-source ci-test-notebooks ci-test-all

download-lipidomics-db:
	# Check if the file already exists before downloading
	@if [ -f "$(LIPIDOMICS_SQLITE_PATH)" ]; then \
		echo "LC-MS lipidomics database already exists at $(LIPIDOMICS_SQLITE_PATH)"; \
	else \
		echo "Downloading LC-MS lipidomics database"; \
		mkdir -p $$(dirname $(LIPIDOMICS_SQLITE_PATH)); \
		curl --retry 3 --retry-delay 5 --connect-timeout 30 --max-time 300 -L -o $(LIPIDOMICS_SQLITE_PATH) $(LIPIDOMICS_SQLITE_URL); \
		echo "LC-MS lipidomics database downloaded"; \
	fi
	
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
	@python3 setup.py sdist

pypi:	
	@rm -rf build dist *.egg-info
	@python3 setup.py sdist
	@twine upload dist/*

tag:

	@git tag -a $(version).$(stage) -m "version $(version).$(stage)"
	@git push origin $(version).$(stage)
	@echo tagged $(version).$(stage) and pushed

build-image:

	@echo corilo/corems:$(version).$(stage)
	@docker build -t corilo/corems:$(version) .
	@docker push corilo/corems:$(version)
	@docker image tag corilo/corems:$(version) corilo/corems:latest
	@docker push corilo/corems:latest

db-up:

	@docker-compose up -d 

db-down:

	@docker-compose down

db-logs:

	@docker-compose logs -f

db-connect:

	@docker exec -it molformdb psql -U postgres

all-up:
	
	@docker-compose -f docker-compose-jupyter.yml up	

fresh-stack-up:

	@docker build -t corems:local .
	@docker-compose up -d   
	@docker run --rm -v ./data:/home/CoreMS/data corems:local

docu:
	
	pdoc --output-dir docs --docformat numpy corems

ci-test-source:
	@$(PYTHON) -V
	@$(PYTHON) -m pip install --upgrade pip
	@$(PYTHON) -m pip install -r requirements.txt
	@$(PYTHON) -m pip install pytest pytest-cov psycopg2
	@$(MAKE) download-lipidomics-db LIPIDOMICS_SQLITE_PATH="$(LIPIDOMICS_SQLITE_PATH)"
	@$(PYTHON) -c "import pathlib; [p.unlink() for p in pathlib.Path('.').rglob('tests/win_only/__init__.py')]"
	@PYTHONNET_RUNTIME=coreclr COREMS_LIPIDOMICS_SQLITE_PATH="$(LIPIDOMICS_SQLITE_PATH)" $(PYTHON) -m pytest --cache-clear

ci-test-notebooks:
	@$(PYTHON) -V
	@$(PYTHON) -m pip install --upgrade pip
	@$(PYTHON) -m pip install -r requirements.txt
	@$(PYTHON) -m pip install jupyter nbconvert psycopg2
	@$(PYTHON) -m pip install -e .
	@PYTHONNET_RUNTIME=coreclr $(PYTHON) examples/test_notebooks.py

ci-test-all: ci-test-source ci-test-notebooks