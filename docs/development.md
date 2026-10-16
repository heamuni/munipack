# Development

This page is for anyone who wants to change Munipack itself: where things
live, how to see what this fork changed, and how to test a change.

## Repository layout

| Path | Contents |
|---|---|
| `src/` | Munipack source, upstream plus the fixes |
| `src/doc/` | Upstream documentation; where it disagrees with these pages, these pages win |
| `deps/oakleaf/` | OakLeaf robust-statistics library, unmodified |
| `Dockerfile` | Builds the command-line tools (see {doc}`install`) |
| `tests/` | Tests on simulated frames (see {doc}`tests`) |
| `docs/` | This documentation (see {doc}`building`) |
| `.github/` | Continuous integration and the runner set-up |
| `BUGS.md` | Known bugs, compatibility notes and ideas |

## What this fork changed

The first commit is the untouched upstream Munipack and OakLeaf, and every
fix is a separate commit on top of it. To see all changes, or the reason for
the changes in one area:

```bash
git diff $(git rev-list --max-parents=0 HEAD) -- src
git log -- src/photometry
```

[BUGS.md](https://github.com/heamuni/munipack/blob/develop/BUGS.md) explains
each bug and points to the commit that fixes it. A summary for users is in
{doc}`changes`.

Several actions share common code (the sky estimate, the PSF models, the
photometric conversions), mostly in `src/lib/`, so a fix there applies to all
of them at once.

## Testing a change

```bash
docker build --target test .        # build and run all tests
tests/run.sh [t_name ...]           # run tests against the munipack on PATH
KEEP=1 tests/run.sh t_pipeline      # keep the working files for inspection
```

The Docker build first runs the small upstream test suite and stops if it
fails. The tests in `tests/` then check each area on frames simulated by
`artificial`, so the true answer is known and no network is needed. The whole
set takes a few minutes; {doc}`tests` describes what each test checks.

## Continuous integration

Every push to `main` or `develop`, and every pull request, is built and
tested automatically on the group's self-hosted runner; a run can also be
started by hand from the Actions tab. The result appears as a pass or fail
mark on the commit or pull request. How the runner is set up is in
[.github/RUNNER.md](https://github.com/heamuni/munipack/blob/develop/.github/RUNNER.md).

## Updating the documentation

The routine pages in `docs/routines/` follow `munipack <action> --help`. When
you change an option, update the help text in the source and the routine page
together; if the behaviour changes, also add a note to BUGS.md and
{doc}`changes`.

## Licence

Munipack is under the GNU General Public License, version 3 or later, and
OakLeaf under the GNU Lesser General Public License. The changes in this
repository use the same licences.
