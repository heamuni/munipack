# Build the documentation: make -C docs html   (or: make -C docs linkcheck)
# Needs the packages of requirements.txt (see README.md).

SPHINXBUILD ?= sphinx-build
SPHINXOPTS  ?= -W --keep-going
SOURCEDIR    = .
BUILDDIR     = _build

.PHONY: help html linkcheck clean

help:
	@echo "make html       build the HTML site in $(BUILDDIR)/html"
	@echo "make linkcheck  check the links"
	@echo "make clean      remove $(BUILDDIR)"

html:
	$(SPHINXBUILD) $(SPHINXOPTS) -b html "$(SOURCEDIR)" "$(BUILDDIR)/html"

linkcheck:
	$(SPHINXBUILD) $(SPHINXOPTS) -b linkcheck "$(SOURCEDIR)" "$(BUILDDIR)/linkcheck"

clean:
	rm -rf "$(BUILDDIR)"
