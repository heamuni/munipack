# Munipack

Munipack is a command-line package for reducing astronomical images, from raw
frames to calibrated magnitudes. Each step stores its results as FITS
extensions in the frame itself, so the next step simply reads what the
previous one wrote.

This is a fork of upstream Munipack 0.6.2 with many bug fixes and a few new
features; see {doc}`changes`. The easiest way to run it is the Docker image:

```bash
docker build -t munipack:0.6.2-hg1653-fixes .
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes find frame.fits
```

New here? Start with {doc}`install` and then the {doc}`quickstart`. The
tutorials then walk through whole tasks step by step:
{doc}`a single image <tutorials/single-image>` from raw frame to calibrated
magnitudes, {doc}`a time series <tutorials/time-series>` ending in a light
curve, and {doc}`stacking and colour images <tutorials/colour-image>`.

## What you can do

- **Calibrate frames** – make master frames with {doc}`routines/bias`,
  {doc}`routines/dark` and {doc}`routines/flat`, and apply them with
  {doc}`routines/phcorr`.
- **Find stars** with {doc}`routines/find`.
- **Measure brightness** – aperture ({doc}`routines/aphot`), PSF
  ({doc}`routines/psfphot`) or growth-curve ({doc}`routines/gphot`)
  photometry.
- **Calibrate positions** – match stars to a catalogue and fit the WCS with
  {doc}`routines/astrometry`.
- **Calibrate magnitudes** – zero point, colour terms and physical fluxes with
  {doc}`routines/phcal`; instrument transformations with
  {doc}`routines/phfotran`.
- **Stack and colour** – stacks and mosaics with {doc}`routines/kombine`,
  colour images with {doc}`routines/colouring`.
- **Simulate frames** – realistic test frames from a catalogue with
  {doc}`routines/artificial`.
- **Use catalogues** – download ({doc}`routines/cone`), resolve names
  ({doc}`routines/sesame`), convert ({doc}`routines/votable`) and
  cross-match ({doc}`routines/cross`).
- **Handle files** – inspect and edit FITS files ({doc}`routines/fits`) and
  make light curves ({doc}`routines/timeseries`).

## Routines

Every routine is run as `munipack <routine> [options] file(s)`, and each has
its own page in the sidebar. The ones used in everyday work are under
*Routines*; the less common ones (instrument transformations, name lookup,
catalogue conversion and cross-matching, FITS utilities) are under *Other
routines*. For the option list on your own machine, run
`munipack <routine> --help`.

`munipack` on its own takes three options: `--version` prints the version and
licence, `-h`/`--help` lists the routines, and `--verbose` makes the routine
report its progress. A routine exits with status 0 on success and 1 on an
error, with the reason printed to stderr, so it is easy to use in scripts.

```{toctree}
:caption: Getting started
:maxdepth: 2

install
quickstart
```

```{toctree}
:caption: User guide
:maxdepth: 2

Single image <tutorials/single-image>
Time series <tutorials/time-series>
Colour images <tutorials/colour-image>
concepts
```

```{toctree}
:caption: Routines
:maxdepth: 1

bias <routines/bias>
dark <routines/dark>
flat <routines/flat>
phcorr <routines/phcorr>
find <routines/find>
aphot <routines/aphot>
gphot <routines/gphot>
psfphot <routines/psfphot>
astrometry <routines/astrometry>
phcal <routines/phcal>
kombine <routines/kombine>
colouring <routines/colouring>
artificial <routines/artificial>
cone <routines/cone>
timeseries <routines/timeseries>
```

```{toctree}
:caption: Other routines
:maxdepth: 1

phfotran <routines/phfotran>
sesame <routines/sesame>
votable <routines/votable>
cross <routines/cross>
fits <routines/fits>
```

```{toctree}
:caption: Development
:maxdepth: 2

development
tests
changes
building
```
