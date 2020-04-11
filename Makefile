app_name = CoreMS
parameters_path = parameter.json 

build:
	@python3 setup.py sdist

bump-release:
	@bumpversion patch --allow-dirty
	@python3 setup.py sdist
	@twine --repository corems upload dist/*
	@twine upload dist/*

release:
	@python3 setup.py sdist
	@twine --repository corems
	@twine upload dist/*	
