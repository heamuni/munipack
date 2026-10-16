# Munipack

Munipack is a command-line package for processing astronomical images in the
optical and ultraviolet: calibration frames, star detection, aperture, PSF and
growth-curve photometry, astrometric calibration, photometric calibration
against a catalogue, stacking, and a simulator of frames with a realistic
detector. Results are stored as FITS extensions of the frames themselves, so
each step reads what the previous one wrote.

This documentation covers a fork of upstream Munipack revision hg1653
(version `0.6.2+hg1653`) with about 60 bug fixes and several new features.
The changes are summarised in {doc}`changes`.

```bash
docker build -t munipack:0.6.2-hg1653-fixes .
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes find frame.fits
```

## Features

::::{grid} 1 2 2 3
:gutter: 2

:::{grid-item-card} Calibration frames
:link: reference/preprocessing
:link-type: doc
Master bias, dark and flat by means clipped by the detector noise model;
`phcorr` applies them with error propagation.
:::

:::{grid-item-card} Detection
:link: reference/find
:link-type: doc
`find`: a threshold set by a false-alarm rate, valid also on the skewed
Poisson noise of faint skies.
:::

:::{grid-item-card} Aperture photometry
:link: reference/aphot
:link-type: doc
`aphot`: unbiased sky at low counts, neighbour checks, and enclosed-energy
aperture corrections from a PSF or from the frame.
:::

:::{grid-item-card} PSF photometry
:link: reference/psfphot
:link-type: doc
`psfphot`: Gaussian, Moffat, tabulated or empirical PSFs, simultaneous fits
of groups, detection in the residuals.
:::

:::{grid-item-card} Growth-curve photometry
:link: reference/gphot
:link-type: doc
`gphot`: total fluxes by the minimum-variance combination of apertures.
:::

:::{grid-item-card} Astrometry
:link: reference/astrometry
:link-type: doc
`astrometry`: matching by the header WCS, by a pointing with an unknown
roll, or blind; TPV distortion; significance test.
:::

:::{grid-item-card} Photometric calibration
:link: reference/phcal
:link-type: doc
`phcal`, `phfotran`: zero point with a systematic error, colour terms,
physical flux units, instrument transformations.
:::

:::{grid-item-card} Stacking and colour
:link: reference/kombine
:link-type: doc
`kombine`: resampled, flux-conserving stacks and mosaics; `colouring`:
colour frames.
:::

:::{grid-item-card} Simulation
:link: reference/artificial
:link-type: doc
`artificial`: frames from a catalogue with Poisson, read and dark noise,
gain, bias, saturation and field-dependent PSFs.
:::

::::

The catalogue tools (`cone`, `sesame`, `votable`, `cross`) and the utilities
(`fits`, `timeseries`) are described in the {doc}`reference/index`.

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
