python = python

env:
	env create -f ./environment.yml -p ./env --quiet

test:
	${python} -m pytest -vrs tests

test-cov-html:
	${python} -m pytest -vrs --cov=rmsd --cov-report html tests
