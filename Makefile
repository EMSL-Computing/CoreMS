app_name = CoreMS
parameters_path = parameter.json 
version := $(shell cat .bumpversion.cfg | grep current_version | cut -d= -f2 | tr -d ' ')
stage := $(shell cat .bumpversion.cfg | grep optional_value | cut -d= -f2 | tr -d ' ') 

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