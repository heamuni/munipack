# psfphot

`psfphot` measures stars by fitting a model of the point spread function (PSF)
to each one.

Use it instead of aperture photometry when the field is crowded or when you
want the best precision on faint stars. Close neighbours are fitted together,
so blended stars get sensible fluxes. The stars come from {doc}`find`, and the
results go into a new `PSFPHOT` table in the same file. The measured flux,
`PSFCOUNT`, is the total flux of the star, so {doc}`phcal` can calibrate it
with `--photometry psf`.

## Usage

```bash
munipack psfphot [options] file(s)
```

Run {doc}`find` first. Running {doc}`aphot` as well is recommended,
because `psfphot` takes the sky annulus and the starting sky from it; without
it, a default annulus just outside the fitting radius is used. The gain, read
noise and saturation level are read from the FITS header (`GAIN`, `READNS` or
`RDNOISE`, `SATURATE`); you can override them with the options below.

## Choosing a PSF

- **`GAUSS`** (default) or **`MOFFAT`**: an analytic profile. Its width (and
  the Moffat β) is fitted on bright, isolated stars unless you give it. Add
  `--elliptic` for stars that are not round.
- **`FILE`**: a PSF file that can change over the field, in the same format
  that {doc}`artificial` `--psf-file` reads (see {ref}`psf-files`). This is
  the natural choice when you have a model of your optics.
- **`EMPIRICAL`**: the PSF is built from the stars of the frame itself, as an
  analytic core plus a table of corrections. It can vary over the field
  (`--psf-var`, `--psf-var-radial`). If you also give `--psf-file`, the file is
  used as a starting point and the frame only corrects it.

An empirical PSF needs enough bright, isolated, unsaturated stars, and they
should be mostly single stars, since any close companion goes into the PSF. If there are
too few for the requested variation, `psfphot` lowers the order and warns you.

## Options

### PSF model

| Option | Default | Description |
|---|---|---|
| `--psf GAUSS\|MOFFAT\|FILE\|EMPIRICAL` | `GAUSS` (`FILE` with `--psf-file`) | The PSF type. |
| `--psf-file file` | — | A PSF file; with `EMPIRICAL`, the base PSF that the frame corrects. |
| `--fwhm pix` | fitted | FWHM of the analytic PSF. |
| `--beta n` | fitted (2.5 with `--fwhm`) | The Moffat exponent β. |
| `--elliptic` | off | Also fit the axis ratio and position angle. |
| `--scale pix/deg` | header (WCS) | Plate scale used to place stars in a PSF file. |
| `--psf-centre x,y` | `NAXIS/2` | Field centre of the PSF file [pix]. |

### Empirical PSF

| Option | Default | Description |
|---|---|---|
| `--psf-core MOFFAT\|GAUSS` | `MOFFAT` | Analytic core of the empirical PSF. |
| `--psf-zoom n` | 4 for FWHM < 5 px, else 2 | Samples per pixel of the correction table. |
| `--psf-radius pix` | 99.5 % of the flux, at least 4 × FWHM | Radius of the correction table. |
| `--psf-var 0\|1\|2` | 0 | Order of the variation in x, y over the field. |
| `--psf-var-radial 0-3` | 0 | Order of the variation with distance from the field centre. |
| `--psf-stars n` | 100 | Maximum number of PSF stars. |
| `--psf-list file` | automatic | Your own PSF stars: `x y` per line, or a FITS table with `X`, `Y`. |
| `--psf-save file` | — | Save the PSF as a PSF file (single input file only). |
| `--psf-smooth n` | by cross-validation | Smoothing of the correction table. |

### Fit

| Option | Default | Description |
|---|---|---|
| `--fitrad pix` | 1.5 × FWHM, min. 3 | Fitting radius. |
| `--sky fixed\|local` | `fixed` | `fixed` takes the sky from the `aphot` annulus; `local` fits it with each star. |
| `--group FWHM` | 1.5 | Stars closer than this are fitted together; `0` turns it off. |
| `--maxiter n` | 50 | Maximum iterations of one fit. |
| `--npass n` | 3 | Number of passes over all stars. |
| `--saturate ADU` | header | Saturation level; pixels at or above it are ignored. |
| `--gain e⁻/ADU` | header | Detector gain. |
| `--read-noise ADU` | header | Read noise. |

### Residuals

| Option | Default | Description |
|---|---|---|
| `--residual-detect n` | 0 | Passes that look for missed stars in the residual image. |
| `--residual-threshold σ` | from the false-alarm rate (about 5 for 2048 × 2048) | Detection threshold in the residuals. |
| `--residual file` | — | Save the image with all fitted stars subtracted (single file only). |

### Output and general

| Option | Description |
|---|---|
| `-t`, `--target-directory dir` | Write output files to this directory. |
| `-O` | Turn on advanced output filename processing. |
| `--pattern str` | Matching pattern (for `-O`). |
| `--mask str` | Mask for output files (for `-O`). |
| `--format str` | Format of output files (for `-O`). |
| `--verbose` | Print the fitted PSF, the PSF stars and the groups. |
| `-h`, `--help` | Show the help. |

## Output

The `PSFPHOT` table has one row per `FIND` star, in the same order. Stars that
`psfphot` adds itself (saturated stars that `find` skipped, and stars found in
the residuals) follow at the end.

| Column | Description |
|---|---|
| `X`, `Y`, `XERR`, `YERR` | Fitted position and its errors [pix]. |
| `PSFCOUNT`, `PSFCOUNTERR` | Total flux and its error [ADU]. |
| `SKY`, `SKYERR` | Sky and its error [ADU]. |
| `CHI2` | Reduced χ² of the fit, about 1 for a good fit. |
| `SHARP` | Positive for objects sharper than the PSF (cosmic rays), negative for extended ones (galaxies, blends). |
| `ROUND` | Ellipticity of the star minus that of the PSF. |
| `NPIX` | Number of fitted pixels. |
| `FLAG` | Bit flags, see below. |
| `GROUP`, `GSIZE` | The group the star was fitted in (0 if alone) and its size. |
| `PSFMODERR` | `EMPIRICAL` only: flux error due to the PSF model, already included in `PSFCOUNTERR`. |

| `FLAG` bit | Meaning |
|---|---|
| 1 | Saturated core; flux from the wings only. |
| 2 | Some pixels in the fitting radius were excluded. |
| 4 | Position not fitted (too faint, or the fit wandered off). |
| 8 | The fit did not converge. |
| 16 | The fit failed; the flux is not valid. |
| 32 | Another star lies within the fitting radius. |
| 64 | Saturated star added by `psfphot`. |
| 128 | Star found in the residuals. |
| 256 | S/N below 2. |

```{note}
`FLAG` is a sum of bits, so test the bits you care about rather than
`FLAG != 0`. Bits 32 and 256 are common in crowded fields and do not make a
measurement bad. `phcal` picks its stars by their relative error, not by `FLAG`.
```

The table header records the PSF that was used (type, FWHM, shape) and the
detector values.

## Examples

Fit a Gaussian PSF with the width found automatically, then calibrate:

```bash
munipack find frame.fits
munipack aphot frame.fits
munipack psfphot --verbose frame.fits
munipack phcal --photometry psf -c cat.fits --col-ra RAJ2000 --col-dec DEJ2000 \
    --col-mag Vmag --col-magerr e_Vmag --photsys-ref Johnson -f V frame.fits
```

Use an elliptical Moffat PSF with a locally fitted sky, and keep the residual
image to check the fit by eye:

```bash
munipack psfphot --psf MOFFAT --elliptic --sky local --residual resid.fits frame.fits
```

Build an empirical PSF that varies over the field, save it, and reuse it on
other frames:

```bash
munipack psfphot --psf EMPIRICAL --psf-var 2 --psf-save epsf.fits frame.fits
munipack psfphot --psf-file epsf.fits other.fits
```

## Notes

- A low-order variation in x, y cannot follow a PSF that changes shape
  strongly across the field. If the change follows the distance from the
  centre (typical of wide-field optics), use `--psf-var-radial`; otherwise
  start from a model PSF file with `--psf EMPIRICAL --psf-file`.
- A file saved with `--psf-save` can also be read by `artificial --psf-file`.
  If the frame had no astrometry, use the file with `--scale 1`.
- Stars closer than about half the FWHM cannot be separated.
