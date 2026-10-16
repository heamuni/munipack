# Development

## Repository layout

| Path | Contents |
|---|---|
| `src/` | Munipack source: upstream hg1653 plus the fixes |
| `src/doc/` | Upstream documentation (HTML man pages `man_<action>.shtml`, data formats `dataform_*.shtml`); where it disagrees with these pages, these pages describe this fork |
| `deps/oakleaf/` | OakLeaf robust-statistics library, upstream hg snapshot (unmodified) |
| `Dockerfile` | Builds the command-line tools (no GUI) on Ubuntu 24.04; stages `build`, `runtime`, `test` (see {doc}`install`) |
| `tests/` | Tests of the installed tools on simulated frames (`tests/run.sh`, see {doc}`tests`) |
| `docs/` | This documentation (see {doc}`building`) |
| `.github/workflows/ci.yml` | Continuous integration |
| `.github/RUNNER.md` | How the self-hosted CI runner is set up |
| `BUGS.md` | Bug catalogue, compatibility notes, suggested improvements |

## History

The history starts with the untouched upstream tree (commit `003cb0d`:
Munipack hg1653 `ae465635f925` and the OakLeaf hg tip). Each fix is a
separate commit on top of it, so

```bash
git diff 003cb0d -- src       # every change to Munipack
git log -- src/photometry     # the reason for each change, by area
```

[BUGS.md](https://github.com/heamuni/munipack/blob/develop/BUGS.md) lists
every bug with its measured impact and the commit that fixes it; the
identifiers (A…, P…, S…) are used in the commit messages, the test comments
and {doc}`changes`.

## Shared routines

Several actions share one implementation, so a change applies to all of
them consistently:

| Routine | Used by |
|---|---|
| `src/lib/clipmean.f08` | The mean clipped by the noise model: the sky of `aphot`, `find`, `psfphot`, the default mean of `kombine`, the master of `dark` |
| `src/lib/psfmodel.f08` | Loading, interpolation and pixelisation of PSFs for `artificial`, `psfphot` and the enclosed energy of `aphot`, so a frame simulated with a PSF is measured with exactly that PSF |
| `src/lib/noise.f08` | Random deviates (Poisson, Normal) of `artificial` |
| `src/photometry/apee.f08` | Enclosed energy of the apertures (`aphot --apcorr`, the growth curve of `gphot`) |
| `src/photometry/exscatter.f08` | The zero point with the systematic error (`phcal --syserr auto`) and the colour terms (`phcal --color-term`, `phfotran --model mag`) |
| `src/lib/photoconv.f08` | Conversions between photons and the calibrated quantities (`phcal -q`) |

## Tests

```bash
docker build --target test .        # build, upstream make check, then tests/run.sh
tests/run.sh [t_name ...]           # the same against the munipack on PATH
KEEP=1 tests/run.sh t_pipeline      # keep the working directories
```

- The build runs the upstream suite (`make check` in `src/test/check.sh`:
  `find`, `aphot`, `astrometry`, `phcal` on the bundled Barnard frame); a
  failure stops the build. Its pass condition is CTPH within twice its error
  of 11.5 (it was 11–12, narrower than the error of about 1.3).
- `tests/` has one test per area, each on frames that `artificial` simulates,
  so the truth is known and no network is needed. The whole set takes 2–3
  minutes. Each test and what it checks is described in {doc}`tests`; the
  comments of each test name the BUGS.md entries it covers.

## Continuous integration

`.github/workflows/ci.yml` runs on every push to `main` or `develop`, on pull
requests and on demand. A newer push to the same branch cancels the older
run. The job:

1. builds the `runtime` stage, which includes the upstream `make check`;
2. builds the `test` stage, which runs `tests/run.sh`;
3. builds the final image and runs a smoke test (`--version`, one
   `artificial` frame);
4. removes the image; the build cache stays, so later runs reuse the
   wxWidgets and OakLeaf layers.

The job runs on the group's self-hosted runner `varda`
(`runs-on: [self-hosted, varda]`); a full run takes about 5 minutes. Every
push and pull request gets a pass or fail mark in the Actions tab. How the
runner is installed, checked and re-registered is in
[.github/RUNNER.md](https://github.com/heamuni/munipack/blob/develop/.github/RUNNER.md).

## Documentation

The pages are in `docs/` and are built by Sphinx (see {doc}`building`).
Options, defaults and units in the reference pages follow
`munipack <action> --help` of this fork. When an option changes, update the
help text in the source, the reference page and, for a change of behaviour,
the compatibility notes of BUGS.md and {doc}`changes`.

## Licence

Munipack is free software under the GNU General Public License, version 3 or
later (`src/COPYING`). OakLeaf is under the GNU Lesser General Public
License (`deps/oakleaf/COPYING.LESSER`). The changes in this repository are
distributed under the same licences.
