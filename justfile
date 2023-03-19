
PROJECT_NAME := "pyddscat"

# List all recipes
list:
  just --list


# Install locally
install:
	pip install -e .


# Install doc requirements
doc-req:
	pip install -r doc/requirements.txt

# Build ddscat
ddscat:
	make -s fortran


# Run the test suite
test:
	pytest test \
	--cov={{PROJECT_NAME}} --cov-append --cov-report term \
	--durations=0

# Clean generated files
clean:
	cd doc && just clean
	cd ddscat && make clean
	rm -rf bin

# Format with black
format:
    black .

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
	git push origin master

# Clean, reformat and push to gitlab
save: clean format gl

# replace:
# 	find . -type f -exec sed -i 's/ScatPy/pyddscat/g' {} +