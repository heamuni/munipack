# Building the documentation

The documentation is a [Sphinx](https://www.sphinx-doc.org/) site written in
[MyST Markdown](https://myst-parser.readthedocs.io/), with the
[Furo](https://pradyunsg.me/furo/) theme (search box, light/dark toggle).
The sources are the `.md` files in `docs/`; `docs/tests.md` includes
`tests/README.md`.

From the top of the repository:

```bash
python3 -m venv .venv-docs
.venv-docs/bin/pip install -r docs/requirements.txt
.venv-docs/bin/sphinx-build -W --keep-going -b html docs docs/_build/html
```

Then open `docs/_build/html/index.html`. With the virtual environment
activated, `make -C docs html` does the same, and `make -C docs linkcheck`
checks the links (it needs network access for the external ones).

`-W --keep-going` turns every warning into an error and lists them all, so a
broken reference or toctree entry fails the build. The output in
`docs/_build/` is ignored by git.

## Writing pages

- Pages are MyST Markdown: fenced code blocks with a language (` ```bash `),
  admonitions as ` ```{note} `, tables in Markdown.
- The extensions `colon_fence`, `deflist`, `fieldlist` and `html_admonition`
  are enabled, and `sphinx_design` for grids and cards.
- Link other pages by a relative path (`[aphot](reference/aphot.md)`) or by
  `` {doc}`reference/aphot` ``; link a section by a label (`(name)=` above
  the heading, then `` {ref}`name` ``). Headings down to level 4 also have
  anchors (`page.md#heading-text`).
- Options, defaults and units must agree with `munipack <action> --help`.
- New pages go into a toctree of `index.md` or `reference/index.md`.
