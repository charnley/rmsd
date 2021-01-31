python = python

test:
	${python} -m pytest -vrs tests

test-cov-html:
	${python} -m pytest -vrs --cov=rmsd --cov-report html tests
