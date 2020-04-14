app_name = CoreMS
parameters_path = parameter.json 
version := $(shell cat .bumpversion.cfg | grep current_version | cut -d= -f2 | tr -d ' ')
stage := $(shell cat .bumpversion.cfg | grep optional_value | cut -d= -f2 | tr -d ' ') 

	
major:
	
	@bumpversion major --allow-dirty

minor:
	
	@bumpversion minor --allow-dirty

patch:
	
	@bumpversion patch --allow-dirty

pypi:	
	
	@python3 setup.py sdist	
	@twine upload dist/*

tag:

	@git tag -a $(version).$(stage) -m "version $(version).$(stage)"
	@git push origin $(version).$(stage)
	@echo tagged $(version).$(stage) and pushed

docker:

	@echo corilo/corems:$(version).$(stage)
	@docker build -t corilo/corems:$(version) .
	@docker push corilo/corems:$(version)
	@docker image tag corilo/corems:$(version) corilo/corems:latest
	@docker push corilo/corems:latest


	
