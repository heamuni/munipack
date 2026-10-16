# Munipack (hg1653 with fixes)

[![CI](https://github.com/heamuni/munipack/actions/workflows/ci.yml/badge.svg?branch=develop)](https://github.com/heamuni/munipack/actions/workflows/ci.yml)

A fixed and extended copy of [Munipack](https://munipack.physics.muni.cz/), Filip Hroch's command-line package for processing optical and UV images: calibration, photometry, astrometry, stacking and frame simulation.

## What is different from upstream

- Many bug fixes, mostly in photometry, astrometry and the frame simulator. Each is its own commit and is described in [BUGS.md](BUGS.md).
- A realistic detector-noise model in `artificial`, so simulated frames look like real ones.
- New PSF photometry (`psfphot`) and a rebuilt growth-curve photometry (`gphot`).
- Photometric calibration (`phcal`) with colour terms and a systematic error.
- Astrometry that can start from the header WCS or from a pointing with an unknown roll, and multi-threaded processing.
- A Docker image of the command-line tools (without the GUI).

Some results differ from upstream because of the fixes. For example, frames measured by an older `aphot` should be measured again before you run the new `phcal`. See the [compatibility notes](BUGS.md#5-compatibility-notes-behaviour-that-changed).

## Build and run

```bash
docker build -t munipack:0.6.2-hg1653-fixes .
alias munipack='docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes'
munipack --help          # list the actions; munipack <action> --help shows the options
```

The first build takes about 10 minutes. To build without Docker, see [docs/install.md](docs/install.md).

## Example

Find stars, measure them, then calibrate the frame's astrometry and photometry against a catalogue downloaded from VizieR:

```bash
munipack cone -r 0.2 -o cone.fits -- 132.85 11.81
munipack find frame.fits
munipack aphot frame.fits
munipack astrometry -c cone.fits frame.fits
munipack phcal -c cone.fits --photsys-ref Johnson -f V \
               --col-mag Vmag --col-magerr e_Vmag frame.fits   # writes frame_cal.fits
```

The [quick start](docs/quickstart.md) walks through the same steps on a simulated frame.

## More information

- [Documentation](https://heamuni.github.io/munipack/): the user guide, tutorials and a page for every action (sources in [docs/](docs/); to build it locally see [docs/README.md](docs/README.md)).
- [tests/README.md](tests/README.md): what the tests check.
- [BUGS.md](BUGS.md): the bugs that were fixed, and the open issues.
- [docs/development.md](docs/development.md): the repository layout and the CI, for developers.

## Licence

Munipack is under the GNU GPL, version 3 or later (`src/COPYING`), and OakLeaf is under the GNU LGPL (`deps/oakleaf/COPYING.LESSER`). The changes in this repository use the same licences.
