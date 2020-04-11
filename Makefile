app_name = CoreMS
parameters_path = parameter.json 

build:
	@python3 setup.py sdist

bump-release:
	@bumpversion minor --allow-dirty
	@python3 setup.py sdist
	@twine upload dist/*

release:
	@python3 setup.py sdist
	@twine upload dist/*	