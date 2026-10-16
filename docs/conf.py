# Sphinx configuration of the Munipack documentation.
# Build: sphinx-build -b html docs docs/_build/html  (see docs/README.md)

project = "Munipack"
author = "Filip Hroch and contributors"
copyright = "Filip Hroch and contributors; GNU GPL v3 or later"
version = "0.6.2"
release = "0.6.2+hg1653 (fork with fixes)"

extensions = [
    "myst_parser",
    "sphinx_copybutton",
    "sphinx_design",
]

source_suffix = {".md": "markdown"}
root_doc = "index"
language = "en"

# docs/README.md is the note on building the site, not a page of it.
exclude_patterns = ["_build", "README.md", "Thumbs.db", ".DS_Store", "requirements.txt"]

# -- MyST --------------------------------------------------------------------
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "fieldlist",
    "html_admonition",
]
# Anchors for headings, so that [text](page.md#heading) links resolve.
myst_heading_anchors = 4

# -- Copy button: strip the prompts of shell examples ------------------------
copybutton_prompt_text = r"\$ "
copybutton_prompt_is_regexp = True
copybutton_only_copy_prompt_lines = False

# -- HTML (furo: search box and light/dark toggle are built in) --------------
html_theme = "furo"
html_title = "Munipack"
html_short_title = "Munipack"
html_static_path = []
html_show_sourcelink = False
html_copy_source = False
html_last_updated_fmt = "%Y-%m-%d"

pygments_style = "friendly"
pygments_dark_style = "monokai"

html_theme_options = {
    "sidebar_hide_name": False,
    "navigation_with_keys": True,
    "top_of_page_buttons": ["view"],
    "source_repository": "https://github.com/heamuni/munipack/",
    "source_branch": "develop",
    "source_directory": "docs/",
    "light_css_variables": {
        "color-brand-primary": "#1f5fa8",
        "color-brand-content": "#1f5fa8",
    },
    "dark_css_variables": {
        "color-brand-primary": "#7fb2ff",
        "color-brand-content": "#7fb2ff",
    },
    "footer_icons": [
        {
            "name": "GitHub",
            "url": "https://github.com/heamuni/munipack",
            "html": (
                '<svg stroke="currentColor" fill="currentColor" stroke-width="0" '
                'viewBox="0 0 16 16"><path fill-rule="evenodd" d="M8 0C3.58 0 0 3.58 0 '
                "8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01"
                ".37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 "
                "1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89"
                "-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64"
                "-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92"
                ".08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 "
                "0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"
                '"></path></svg>'
            ),
            "class": "",
        },
    ],
}
