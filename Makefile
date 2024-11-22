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

minor:
	
	@bumpversion minor --allow-dirty

patch:
	
	@bumpversion patch --allow-dirty

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