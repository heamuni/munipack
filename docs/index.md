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

New here? Start with {doc}`install` and then the {doc}`quickstart`.

## What you can do

- **Calibrate frames** – make master bias, dark and flat frames and apply
  them (`bias`, `dark`, `flat`, `phcorr`): {doc}`reference/preprocessing`.
- **Find stars** (`find`): {doc}`reference/find`.
- **Measure brightness** – aperture (`aphot`, {doc}`reference/aphot`), PSF
  (`psfphot`, {doc}`reference/psfphot`) or growth-curve (`gphot`,
  {doc}`reference/gphot`) photometry.
- **Calibrate positions** – match stars to a catalogue and fit the WCS
  (`astrometry`): {doc}`reference/astrometry`.
- **Calibrate magnitudes** – zero point, colour terms, physical fluxes and
  instrument transformations (`phcal`, `phfotran`): {doc}`reference/phcal`.
- **Stack and colour** – combine frames into stacks, mosaics or colour images
  (`kombine`, `colouring`): {doc}`reference/kombine`.
- **Simulate frames** – make realistic test frames from a catalogue
  (`artificial`): {doc}`reference/artificial`.

The catalogue tools (`cone`, `sesame`, `votable`, `cross`) and the utilities
(`fits`, `timeseries`) are listed in the {doc}`reference/index`.

```{toctree}
:caption: Getting started
:maxdepth: 2

install
quickstart
```

```{toctree}
:caption: User guide
:maxdepth: 2

concepts
validation
```

```{toctree}
:caption: Reference
:maxdepth: 2

reference/index
```

```{toctree}
:caption: Development
:maxdepth: 2

development
tests
changes
building
```
