
env:
	conda env create -f ./environment.yml -p ./env --quiet
	./env/bin/python -m pip install -e .

test:
	python -m pytest -vrs tests

test-cov-html:
	python -m pytest -vrs --cov=rmsd --cov-report html tests
