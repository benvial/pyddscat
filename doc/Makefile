# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build
HTML_FILES = `find $(BUILDDIR)/html/ -type f -name "*.html" -print`

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)



	
show:
	@python -c "import webbrowser; \
	webbrowser.get('firefox').open_new_tab('$(PWD)/$(BUILDDIR)/html/index.html')"

show-pdf:
	@xdg-open $(PWD)/$(BUILDDIR)/latex/klove.pdf

postpro:
	@python postprocess.py html $(HTML_FILES)
