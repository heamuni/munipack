# artificial

`artificial` makes realistic synthetic images of the sky.

Use it to test a reduction pipeline, or to see how an observing effect behaves
when you know the true answer. Stars (from a catalogue or a random field),
galaxies, a variable star and the sky are seen through a model telescope, PSF
and detector, with noise, cosmic rays and, optionally, the atmosphere.

Every frame carries a WCS, the PSF it was rendered with, and a table of the
true star positions and fluxes, so you can compare your results with the truth.

## Usage

```bash
munipack artificial [options]
```

There is no input file. The output is `artificial.fits`, or
`artificial_001.fits`, `artificial_002.fits`, … for a series (`--count`).
Existing files with the same name are overwritten.

## Options

### PSF

| Option | Default | Description |
|---|---|---|
| `--psf` | `SEEING` | PSF shape: `SEEING` (Airy disk blurred by seeing, most natural), `MOFFAT` or `GAUSS` (fast). |
| `--psf-file` | — | Use a tabulated PSF from a FITS file instead; see {ref}`psf-files`. |
| `--hwhm` | 1 | PSF width in pixels: σ for `GAUSS` and `SEEING`, α for `MOFFAT` (not the true HWHM). |
| `--beta` | 2 | Moffat exponent β. |
| `--airy` | from aperture | Airy disk radius in pixels (`SEEING` only). |
| `--spread` | `AUTO` | How seeing is applied for `SEEING`: `FFT` (smooth), `RANDOM` (short-exposure look) or `AUTO` (`RANDOM` below 3 s). |
| `--eccentricity` | 0 | Elongate the PSF, e.g. to mimic tracking errors. |
| `--inclination` | 0 | Direction of the elongation [deg from the x axis]. |

### Telescope and detector

| Option | Default | Description |
|---|---|---|
| `--exptime` | 1 | Exposure time [s]. |
| `--diameter` | 1.12 | Aperture diameter [m]. |
| `--area` | from diameter | Collecting area [m²]. |
| `--qeff` | 1 | Total quantum efficiency (optics, filter, detector). |
| `--dark` | 0 | Dark current [e⁻/s/pix]. |
| `--readnoise` | 0 | Read-out noise [e⁻ rms]. |
| `--gain` | 1 | Gain [e⁻/ADU]; frames are written in ADU. |
| `--bias` | 0 | Bias level added after the gain [ADU]. |
| `--saturate` | max of `-B` | Saturation level [ADU]. |
| `--disable-noise` | off | Write the noise-free mean frame. |
| `--gnoise` | 0 | Deprecated extra background rate; use `--dark` and `--readnoise`. |
| `--flat-model-file` | — | FITS image of pixel sensitivity (flat field), same size as the frame. |
| `-B`, `--bitpix` | `-32` | Output type: `-32` (float), `32`, `16` or `8` (unsigned integers). |

### Cosmic rays

| Option | Default | Description |
|---|---|---|
| `--cosmics-freq` | 0 | Mean rate of cosmic-ray hits [1/s]. |
| `--cosmics-cmax` | 1000 | Peak electrons along a track. |
| `--cosmics-len` | 10 | Maximum track length [pix]. |

### Sky

| Option | Default | Description |
|---|---|---|
| `--sky-mag` | 21 | Sky brightness at the frame centre [mag/arcsec²]; `99` for no sky. |
| `--sky-grad-x`, `--sky-grad-y` | 0 | Sky gradient [mag/arcsec²/pix], positive is fainter; 0.0001–0.001 is typical. |

### Stars and galaxies

| Option | Default | Description |
|---|---|---|
| `-c`, `--cat` | — | Star catalogue (e.g. from `munipack cone`); without it a random field is made. |
| `--maglim` | 13 | Faintest magnitude of the random field. |
| `--col-ra`, `--col-dec` | from `--vocat` | Coordinate columns [deg]. |
| `--col-pm-ra`, `--col-pm-dec` | from `--vocat` | Proper-motion columns (must exist, but stars are not moved). |
| `--col-epoch-ra`, `--col-epoch-dec` | from `--vocat` | Epoch columns (accepted, not used). |
| `--col-mag` | `<filter>mag`, e.g. `Vmag` | Magnitude column in the filter `-f`. |
| `--vocat` | installed | Catalogue configuration with default column names. |
| `--galcat` | — | Catalogue of elliptical (Sérsic) galaxies. |

### Photometric system

| Option | Default | Description |
|---|---|---|
| `-f`, `--filter` | `V` | Filter of the frame (for stars and sky). |
| `--photsys` | `Johnson` | Photometric system of the filter, e.g. `Gunn` for ugriz. |
| `--phsystab` | installed | Table of photometric systems. |

### Geometry

| Option | Default | Description |
|---|---|---|
| `--width`, `--height` | 666, 666 | Frame size [pix]. |
| `--rcen`, `--dcen` | 0, 0 | RA and Dec of the frame centre [deg]. |
| `--scale` | 3600 | Plate scale in **pixels per degree** (3600 = 1″/px). |
| `--fov` | — | Field of view along the height [deg]; overrides `--scale`. |
| `--angle` | 0 | Rotation of the frame [deg], clockwise positive. |
| `-p`, `--projection` | `TAN` | WCS projection. |

### Variable star

A variable star is added when you give one of `--lc-model`, `--lc-table` or
`--lc-fourier`.

| Option | Default | Description |
|---|---|---|
| `--lc-model` | — | Built-in light curve: `wave`, `Cepheid` or `Algol`. |
| `--lc-table` | — | Light curve as a table of phase, relative flux and error. |
| `--lc-fourier` | — | Light curve as a table of Fourier coefficients. |
| `--lc-mag` | 15 | Magnitude of the variable. |
| `--lc-amp` | 1 | Amplitude of the curve. |
| `--lc-jd0` | 2451544.5 | Reference epoch (full JD, not MJD). |
| `--lc-per` | 1 | Period [days]. |
| `--lc-ra`, `--lc-dec` | frame centre | Position of the variable [deg]. |

### Atmosphere and time series

| Option | Default | Description |
|---|---|---|
| `--atmosphere` | off | Add extinction and seeing that grow with airmass. |
| `--extk` | 0.1 | Extinction coefficient. |
| `--long`, `--lat` | 0, 0 | Site longitude (+east) and latitude (+north) [deg]. |
| `--date` | 2000-01-01 | Start date (UT), `YYYY-MM-DD`. |
| `--time` | 00:00:00 | Start time (UT), `HH:MM:SS`. |
| `--count` | 1 | Number of frames. |
| `--timestep` | 60 | Time between frame starts [s]. |

### Output and runs

| Option | Default | Description |
|---|---|---|
| `--mask` | `artificial_???.fits` | Output name; a run of `?` is replaced by the frame number. |
| `--seed` | random | Random seed, for reproducible frames. |
| `--no-watermark` | off | Leave out the small mark in the top-right corner. |
| `--fits-keys-file` | — | Extra header cards, one `KEY = value / comment` per line; they win over the built-in ones. |
| `-mproc`, `--mproc` | 3 (1 in Docker) | Number of parallel processes for a series. |
| `--Xproc` | — | Extra options for the parallel launcher. |
| `--verbose` | off | Print parameters and progress. |
| `-h`, `--help` | — | Show help. |

## Output

Each file holds the image (in ADU) and two extensions: `PSF`, the PSF used
(not written with `--psf-file`), and `CATALOGUE`, the true `X`, `Y` and `COUNT`
(electrons, before the gain) of every star, with the variable first. The header records the
exposure, the detector settings, the WCS and, with `--atmosphere`, the airmass.

## Examples

A reproducible random field on a 16-bit CCD-like detector:

```bash
munipack artificial --width 1024 --height 1024 --maglim 16 \
    --psf GAUSS --hwhm 1.5 --exptime 60 \
    --gain 2 --readnoise 5 --bias 100 -B 16 \
    --seed 42 --no-watermark --mask field.fits
```

A real field from a catalogue, for testing `find`, `aphot` and `astrometry`
against the truth:

```bash
munipack cone -r 0.2 -o cat.fits -- 180.0 30.0
munipack artificial -c cat.fits --rcen 180 --dcen 30 \
    --col-ra RAJ2000 --col-dec DEJ2000 --col-pm-ra pmRA --col-pm-dec pmDE \
    --col-mag Vmag --width 640 --height 640 --exptime 60 --mask frame.fits
```

A night-long series of an eclipsing binary seen from a real site, writing
`algol_01.fits` … `algol_36.fits`:

```bash
munipack artificial --lc-model Algol --lc-mag 12 --lc-per 2.867 \
    --lc-jd0 2460000.42 --rcen 47.04 --dcen 40.96 \
    --atmosphere --long 16.6 --lat 49.2 --date 2023-02-24 --time 20:00:00 \
    --count 36 --timestep 300 --mask 'algol_??.fits'
```

(psf-files)=
## PSF files

`--psf-file` reads the PSF file format described in {ref}`psf-format`: one
PSF for the whole field, or several PSFs at different places that
`artificial` interpolates between. `psfphot --psf EMPIRICAL --psf-save`
writes this format, so you can re-render a field with the PSF measured on a
real frame. Use the same `--scale`, `--width` and `--height` as that frame,
because the PSF positions are tied to the field geometry. A PSF saved from a
frame without a WCS needs `--scale 1`.

## Notes

- Paired options (`--width`/`--height`, `--rcen`/`--dcen`, `--long`/`--lat`,
  `--lc-ra`/`--lc-dec`) must be given together or not at all.
- Give `?` in `--mask` for a series; without it every frame overwrites the
  same file.
- Give only one light-curve source. If several are given, `--lc-fourier`
  wins over `--lc-table`, which wins over `--lc-model`.
- Templates for the input tables are in the source tree: `src/doc/artgalaxy.lst`
  for `--galcat` and `src/doc/fourier.lst` for `--lc-fourier`. Turn them into
  FITS with `munipack fits --restore <file>.lst`.
- Use `--bias` with integer output (`-B 16`) and read noise, so that
  negative noise values are not clipped at zero.
- With `--atmosphere`, frames in which the field is below the horizon are
  skipped with a warning.
- With `--seed`, the frames are the same whatever `--mproc` you use.
- `find` detects the watermark as a few stars; use `--no-watermark` when you
  compare detections with the catalogue.
