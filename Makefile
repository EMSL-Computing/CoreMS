init:
    pip install -r requirements.txt

test:
    python setup.py test

build:
    - python setup.py sdist bdist_wheel

.PHONY: init test build
