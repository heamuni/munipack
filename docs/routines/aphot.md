# aphot

`aphot` measures the brightness of stars by summing their light in circular
apertures and subtracting the sky, which it measures in a ring around each
star.

It works on the stars that {doc}`find` detected, so run `find` first. The
results go into a new `APERPHOT` table in the same FITS file. This table is
what {doc}`phcal` calibrates, and it is also used by {doc}`gphot` and
{doc}`psfphot`. Optionally, `aphot` can also work out the aperture
corrections, so that `phcal` calibrates the total flux of each star rather
than just the part inside the aperture.

## Usage

```bash
munipack aphot [options] file(s)
```

The file is updated in place, and any earlier `APERPHOT` table is replaced.
Use `-t dir` to write the result to another directory instead.

## Options

### Apertures and sky

| Option | Default | Description |
|---|---|---|
| `--apertures r1,r2,...` | 12 radii, 1 to 20 px | Aperture radii in pixels, in increasing order. |
| `--ring ri,ro` | `20,30` | Inner and outer radius of the sky ring in pixels; keep the inner one beyond the largest aperture. |
| `--enable-ellipticity` | off | Use elliptic apertures, shaped by the star shapes from `find` (experimental). |
| `--eccentricity e` | from `find` | Eccentricity of the elliptic apertures (needs `--enable-ellipticity`). |
| `--inclination deg` | from `find` | Angle of the major axis from the x axis (needs `--enable-ellipticity`). |
| `--coo "x1,y1 x2,y2 ..."` | — | Measure these pixel positions instead; prints the results and leaves the file alone. |

### Aperture corrections from a known PSF

| Option | Default | Description |
|---|---|---|
| `--psf-file file` | — | Use a field-dependent PSF file ({ref}`format <psf-format>` as for `artificial` and `psfphot`); implies `--psf FILE`. |
| `--psf type` | — | Use a `GAUSS` or `MOFFAT` PSF, or `FILE` together with `--psf-file`. |
| `--fwhm pix` | from `find` | FWHM of the Gaussian or Moffat PSF. |
| `--beta b` | 2.5 | Moffat exponent. |
| `--scale pix/deg` | from the WCS | Plate scale for the PSF file, needed when the frame has no WCS. |
| `--psf-centre x,y` | frame centre | Field centre of the PSF file, also used by `--apcorr-radial`. |

### Aperture corrections from the frame's own stars

| Option | Default | Description |
|---|---|---|
| `--apcorr auto` | — | Derive the corrections from the growth curves of bright, isolated stars. |
| `--apcorr-var n` | 2 | Highest polynomial order (0, 1 or 2) of the variation over the field. |
| `--apcorr-radial` | off | Let the correction vary with distance from the field centre rather than with x,y. |

### Output and general

| Option | Description |
|---|---|
| `-t`, `--target-directory dir` | Write the output into this directory and leave the inputs unchanged. |
| `-O` | Turn on advanced output names (`--pattern`, `--mask`, `--format`). |
| `--pattern regex` | Pattern matched against input names (with `-O`). |
| `--mask string` | Output name built from the pattern's groups `\1`, `\2`, ... (with `-O`). |
| `--format fmt` | Numbered output names such as `out_%02d.fits` (with `-O`). |
| `--verbose` | Print progress messages. |
| `-h`, `--help` | Show the option summary. |

## How it works

**Sky.** The sky is the mean of the pixels in the ring, with the other stars
masked out and outliers (stars, cosmic rays, hot pixels) clipped. A mean is
used rather than a median because the aperture contains the sky's mean
level, and on faint skies the two differ noticeably.

**Counts and errors.** Pixels on the edge of an aperture count by the
fraction that lies inside it. The error combines the photon noise of the
star, the noise of the sky, and the uncertainty of the sky level. It uses
the GAIN from the header; if GAIN is missing, 1 e⁻/ADU is assumed and a
warning is printed.

**Invalid apertures** are marked with −1. This happens when an aperture runs
off the frame, contains saturated or bad pixels, has negative counts, or is
noticeably contaminated by a brighter neighbour. That way faint stars next
to bright ones do not spoil the calibration.

**Aperture corrections** are written only when you ask for them. For each
aperture, they give the fraction of the star's light that falls inside it
(the "enclosed energy", APEE). With a known PSF the corrections are exact.
With `--apcorr auto` they come from the growth curves of the frame's bright,
isolated stars and are allowed to vary smoothly over the field. If there are
too few suitable stars, `aphot` warns you and writes no corrections.

```{tip}
For `--apcorr auto`, let the apertures reach out to the sky ring, so that
the growth curves level off before the largest aperture, for example
`--apertures 1,1.5,2,3,4,5,7,9,12,15,20,23,26.4,30.4,35,40 --ring 45,60`.
It works best with about 30 or more bright, isolated stars and a correct
FWHM from `find`. On a nearly empty frame give a PSF with `--psf-file`
instead.
```

## Output

The `APERPHOT` table has one row per star, in the same order as in `find`:

| Column | Description |
|---|---|
| `X`, `Y` | Position of the star [pix]. |
| `SKY`, `SKYERR` | Sky level and its error [ADU/pix]. |
| `APCOUNT1` ... `APCOUNTn` | Sky-subtracted counts in each aperture [ADU]; −1 if invalid. |
| `APCOUNTERR1` ... `APCOUNTERRn` | Errors of the counts [ADU]. |
| `APEE1` ... `APEEn`, `APEEERR1` ... | Enclosed energy and its error (only with the correction options). |

The header records the aperture and ring radii (`APER1`..., `ANNULUS1`,
`ANNULUS2`), the star width `HWHM`, the gain used, and, when corrections are
computed, how they were obtained (keywords starting with `EE`).

With `--coo`, nothing is written to the file; one line per position is
printed to the terminal, starting with `=BYHAND>`.

## Examples

Detect the stars and measure them with the default apertures:

```bash
munipack find frame.fits
munipack aphot frame.fits
```

Choose your own apertures and sky ring, and keep the original file untouched:

```bash
munipack aphot --apertures 2,3,4,6,8,12 --ring 15,25 -t phot/ frame.fits
```

Add aperture corrections from a Gaussian PSF, or from the stars in the frame:

```bash
munipack aphot --psf GAUSS --fwhm 4.9 frame.fits
munipack aphot --apertures 1,1.5,2,3,4,5,7,9,12,15,20,23,26.4,30.4,35,40 \
               --ring 45,60 --apcorr auto frame.fits
```

## Notes

- `phcal` uses the aperture corrections automatically when they are present
  (`phcal --apcorr none` ignores them). If {doc}`gphot` has added a
  `GROWPHOT` table, `phcal` calibrates that by default instead
  (`phcal --photometry aper` selects the apertures).
- `--coo` cannot be combined with the aperture-correction options, and
  `--apcorr` cannot be combined with the PSF options.
- If your frames were measured with an older version of `aphot`, run it
  again before calibrating: older versions stored a different star width in
  `HWHM`, which changes the aperture that `phcal` picks.
