python=./env/bin/python
pytest=./env/bin/pytest
package=rmsd

## Setup

env:
	conda env create -f ./environment.yml -p ./env --quiet
	${python} -m pre_commit install
	${python} -m pip install -e .

## Development

update-format:
	${python} -m pre_commit autoupdate

format:
	${python} -m pre_commit run --all-files

test:
	${python} -m pytest ./tests

test-dist:
	${python} -m twine check dist/*

types:
	${python} -m monkeytype run $$(which ${pytest}) ./tests
	${python} -m monkeytype list-modules | grep ${package} | parallel -j1 "${python} -m monkeytype apply {} > /dev/null && echo {}"

cov:
	${python} -m pytest --cov=${package} --cov-config .coveragerc --cov-report html tests

build:
	${python} -m build --sdist --skip-dependency-check  .

upload:
	${python} -m twine upload ./dist/*.tar.gz

## Version

## Github

## Clean

clean:
	find ./rmsd/ -type f \
		-name "*.so" \
		-name "*.pyc" \
		-name ".pyo" \
		-delete
	rm -rf ./rmsd/*.egg-info/
	rm -rf *.whl
	rm -rf ./build/ ./__pycache__/
	rm -rf ./dist/

clean-env:
	rm -rf ./env/
	rm ./.git/hooks/pre-commit
