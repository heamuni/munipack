# Building the documentation

This site is built with [Sphinx](https://www.sphinx-doc.org/) from the
Markdown ([MyST](https://myst-parser.readthedocs.io/)) files in `docs/`.
To build it, run from the top of the repository:

```bash
python3 -m venv .venv-docs
.venv-docs/bin/pip install -r docs/requirements.txt
.venv-docs/bin/sphinx-build -W --keep-going -b html docs docs/_build/html
```

Then open `docs/_build/html/index.html` in a browser. With the virtual
environment activated, `make -C docs html` does the same, and
`make -C docs linkcheck` checks the links (this needs network access).
Any warning, such as a broken reference, makes the build fail, so problems
are caught early.

## Writing pages

Pages are ordinary Markdown with a few MyST extras, such as ` ```{note} `
boxes. Link to another page with `` {doc}`routines/aphot` `` or a relative
path, and to a section by putting a label like `(name)=` above its heading
and writing `` {ref}`name` ``. Add new pages to the table of contents in
`index.md`; a new routine page goes into the Routines toctree there. The
Tests page is pulled in from `tests/README.md`, so edit that file to change
it.

When you document an option, check it against `munipack <action> --help`
so the names and defaults stay correct.
