# gphot

`gphot` measures the total flux of every star in a frame from a growth curve.

The growth curve (the fraction of a star's light inside each aperture) is
built from the bright, isolated stars and then used to estimate the total
flux of all stars. Faint stars get small, low-noise apertures and bright
stars get large ones, all corrected to the same total, so the fluxes stay
accurate across a wide range of brightness. {doc}`phcal` then calibrates the
result; it picks these fluxes by default when the frame has them
(`--photometry grow`).

## Usage

```bash
munipack gphot [options] file(s)
```

Run {doc}`find` and {doc}`aphot` first. `aphot` should use many apertures,
from about 1 FWHM out to several FWHM, all inside the sky ring, so the curve
can reach its plateau. `gphot` needs at least 3 bright isolated stars and
stops with an error if it finds fewer.

By default the input file is updated in place. Use `-t` to write into another
directory, or give input and output as one quoted argument separated by a
space: `'frame.fits result.fits'`. The comma form shown by `--help` does not
work.

## Options

### Growth curve

| Option | Argument | Default | Description |
|---|---|---|---|
| `--var` | 0, 1, 2 | 2 | Highest polynomial order for how the curve varies over the field; raised only while significant. Not with `--legacy`. |
| `--radial` | — | off | Let the curve vary with distance from the field centre instead of x,y. Not with `--legacy`. |
| `--centre` | `x,y` [pix] | NAXIS/2 | Field centre for `--radial` (requires `--radial`). |
| `--legacy` | — | off | Use the growth-curve fit of earlier Munipack versions (before 0.6.2+hg1653). |
| `-th`, `--threshold` | real | 0.1 | With `--legacy` only: curve stars must be brighter than this times the sky in the aperture. |

### Output and general

| Option | Argument | Default | Description |
|---|---|---|---|
| `-t`, `--target-directory` | directory | — | Write results into this (existing) directory instead of updating the input. |
| `-O` | — | off | Turn on advanced output filenames. |
| `--pattern` | regex | `(.+)\.(.+)` | Pattern matched against input filenames (with `-O`). |
| `--mask` | string | — | Output filename built from the pattern's groups (with `-O`). |
| `--format` | printf format | — | Format for numbered output files, e.g. `out_%02d.fits` (with `-O`). |
| `--verbose` | — | off | Print the fitted growth curve and how it was used. |
| `-h`, `--help` | — | — | Show the help and exit. |

## Output

`gphot` adds two tables to the frame, replacing them if they already exist:

- **`GROWPHOT`**, one row per star: position `X`, `Y`, `SKY`, `SKYERR`, the
  total flux `GCOUNT` [ADU] with its error `GCOUNTERR`, and `GROWFLAG`
  (1 for stars used to build the curve, 0 otherwise). Stars without a valid
  flux have `GCOUNT = GCOUNTERR = -1`. The header records the curve's `HWHM`,
  the 90%-flux radius `RADFLX90`, and `EE*` keys describing the fit.
- **`GROWCURVE`**, one row per aperture: radius `R`, the growth curve
  `GROWCURVE` with its error, and the radial profile `RADIALPROFILE`.

## Examples

The full chain from detection to calibration:

```bash
munipack find frame.fits
munipack aphot --apertures 1,1.313,1.7241,2.2637,2.9724,3.9028,5.1245,6.7287,8.8349,11.601,15.232,20,23,26.4,30.4,35,40 \
               --ring 45,60 frame.fits
munipack gphot frame.fits
munipack phcal --photometry grow --photsys-ref Johnson -f V \
               --col-ra RAJ2000 --col-dec DEJ2000 --col-mag Vmag --col-magerr e_Vmag \
               -c catalogue.fits frame.fits
```

The calibrated stars end up in `frame_cal.fits`.

Wide-field optics, where the PSF changes with distance from the optical axis:

```bash
munipack gphot --radial --centre 1024,1024 frame.fits
```

A constant growth curve, written to a new file, with the fit printed:

```bash
munipack gphot --var 0 --verbose 'frame.fits frame_grow.fits'
```

## Notes

- `GCOUNT` is the flux inside the largest (plateau) apertures, not an
  extrapolation to infinity. Light in faint, extended PSF wings is missed,
  but it is the same fraction for every star, so the photometric calibration
  absorbs it.
- The growth curves come from bright, isolated stars, and isolation is
  judged from the FWHM written by {doc}`find`. If that FWHM is badly off, for
  example because stars were split into several detections, check the `find`
  result first. On a frame with very few stars, use `aphot --psf-file` or
  {doc}`psfphot` instead.
- {doc}`psfphot` is the other way to get total fluxes.
