
PROJECT_NAME := "pyddscat"
PROJECT_DIR := `pwd`
VERSION := ```python3 -c "import toml;print(toml.load('pyproject.toml')['project']['version'])"```
URL := ```python3 -c "import toml;print(toml.load('pyproject.toml')['project']['urls']['Repository'])"```
BRANCH := `git branch --show-current`
GITLAB_PROJECT_ID := "44436558"
GITLAB_GROUP_ID := "64380746"
LINT_FLAGS := "E501,F401,F403,F405,W503,E402,E203"


# List all recipes
list:
  just --list


# Install locally
install-loc:
	pip install -e .

# Install locally
install:
	pip install .


# Install doc requirements
doc-req:
	pip install -r doc/requirements.txt

# Build ddscat
ddscat:
	make -s ddscat



# Clean generated files
clean:
	cd doc && just clean
	cd ddscat && make clean
	rm -rf *.whl
	rm -rf wheelhouse build


# Build documentation
doc:
	cd doc && just build

# Build documentation and watch
autodoc:
	cd doc && just autobuild

# Push to gitlab
gl:
	@echo "Pushing to gitlab..."
	git add -A
	@read -p "Enter commit message: " MSG; \
	git commit -a -m "$MSG"
	git push origin main

# Clean, reformat and push to gitlab
save: clean format gl



# Install development dependencies
dev:
    pip install -r dev/requirements.txt

# Format with black
format:
    black .

# Lint with flake8
lint:
	@flake8 --exit-zero --ignore={{LINT_FLAGS}} {{PROJECT_NAME}} 
	
# Lint using flake8
lint-extra:
	@flake8 --exit-zero --ignore={{LINT_FLAGS}} {{PROJECT_NAME}}  test/ examples/ --exclude "dev*"

# Check for duplicated code
dup:
	@pylint --exit-zero -f colorized --disable=all --enable=similarities {{PROJECT_NAME}}


# Clean test coverage reports
cleantest:
	@rm -rf .coverage* htmlcov coverage.xml


# Install requirements for testing
test-req:
	@cd test && pip install -r requirements.txt


# Run the test suite
test:
	pytest test \
	--cov={{PROJECT_NAME}}/ --cov-report term --cov-report html --cov-report xml  \
	--durations=0

#  Update header
header:
	cd dev && python update_header.py

# Show html documentation in the default browser
show:
    cd doc && make -s show