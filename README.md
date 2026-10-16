# Munipack (hg1653 with fixes)

[![CI](https://github.com/heamuni/munipack/actions/workflows/ci.yml/badge.svg?branch=develop)](https://github.com/heamuni/munipack/actions/workflows/ci.yml)

A fixed and extended copy of [Munipack](https://munipack.physics.muni.cz/),
the astronomical image-processing package by Filip Hroch (Masaryk
University), for optical and UV imaging: calibration frames, detection,
aperture, PSF and growth-curve photometry, astrometry, photometric
calibration, stacking and frame simulation, all from the command line.

It is upstream Mercurial revision **1653** (`ae465635f925`, reported as
`0.6.2+hg1653`) plus:

- about 60 bug fixes, each a separate commit and each listed with its
  measured impact in [BUGS.md](BUGS.md);
- new features: a detector-noise model in `artificial`, PSF photometry
  (`psfphot`), enclosed-energy aperture corrections, a rebuilt `gphot`, a
  systematic error and colour terms in `phcal`, astrometry by the header WCS
  or a pointing with an unknown roll, OpenMP and coarray parallelism;
- a Docker build of the command-line tools (no GUI).

## Status

Every push is tested by the CI (upstream `make check` and [tests/](tests/README.md)).
Some behaviour differs from upstream (e.g. `--scale` is in pixels per degree;
re-run `aphot` before the new `phcal`): see the
[compatibility notes](BUGS.md#5-compatibility-notes-behaviour-that-changed)
and the [open issues](BUGS.md#6-suggested-further-improvements-not-done).

## Build and run

```bash
docker build -t munipack:0.6.2-hg1653-fixes .
alias munipack='docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes'
munipack --help                    # the actions; munipack <action> --help for the options
docker build --target test .       # optionally, build and run the tests
```

The first build takes about 10 minutes. Native builds are described in
[docs/install.md](docs/install.md).

## Example

Detection, aperture photometry, astrometry and calibration of a V frame of
M67 against UCAC4 (APASS magnitudes) from VizieR:

```bash
munipack cone -r 0.2 -o cone.fits -- 132.85 11.81
munipack find frame.fits
munipack aphot frame.fits
munipack astrometry -c cone.fits frame.fits
munipack phcal -c cone.fits --photsys-ref Johnson -f V \
               --col-mag Vmag --col-magerr e_Vmag frame.fits   # writes frame_cal.fits
```

The [quick start](docs/quickstart.md) runs this chain on a simulated frame,
where the true magnitudes are known.

## Documentation

- **[docs/](docs/)**: the user documentation (Sphinx): installation, quick
  start, concepts, a reference page for every action, measured accuracy.
  Build it with `pip install -r docs/requirements.txt` and
  `sphinx-build -b html docs docs/_build/html` ([docs/README.md](docs/README.md)).
- **[tests/README.md](tests/README.md)**: what each test checks.
- **[BUGS.md](BUGS.md)**: every bug, its impact and fix; compatibility notes.
- **[docs/development.md](docs/development.md)**: repository layout,
  history, CI ([.github/RUNNER.md](.github/RUNNER.md)).

## Licence

Munipack is free software under the GNU General Public License, version 3 or
later (`src/COPYING`). OakLeaf is under the GNU Lesser General Public License
(`deps/oakleaf/COPYING.LESSER`). The changes in this repository are
distributed under the same licences.
