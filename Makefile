app_name = CoreMS
parameters_path = parameter.json 
version := $(shell cat .bumpversion.cfg | grep current_version | cut -d= -f2 | tr -d ' ')
stage := $(shell cat .bumpversion.cfg | grep optional_value | cut -d= -f2 | tr -d ' ') 
LIPIDOMICS_SQLITE_URL ?= https://nmdcdemo.emsl.pnnl.gov/lipidomics/parameter_files/202412_lipid_ref.sqlite
LIPIDOMICS_SQLITE_PATH ?= tests/tests_data/lcms/202412_lipid_ref.sqlite

.PHONY: download-lipidomics-db

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

	@echo corems:$(version).$(stage)
	@docker build -t corems:$(version) .

build-image:

	@echo corilo/corems:$(version).$(stage)
	@docker build -t corilo/corems:$(version) .

build-image-mac:

	@echo corilo/corems:$(version).$(stage)
	@docker build --platform linux/amd64 -t corilo/corems:$(version) .

build-image-mac-local:

	@echo corems:$(version).$(stage)
	@docker build --platform linux/amd64 -t corems:$(version) .

push-image:

	@docker push corilo/corems:$(version)
	@docker image tag corilo/corems:$(version) corilo/corems:latest
	@docker push corilo/corems:latest

image-run-mac:

	@docker run -it --platform linux/amd64 corilo/corems:$(version)

image-run:

	@docker run -it corilo/corems:$(version)

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