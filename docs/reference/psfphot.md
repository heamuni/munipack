# psfphot

`psfphot` fits a point spread function (PSF) to every star detected by
[`find`](find.md) and writes the fluxes to a new `PSFPHOT` extension. The PSF can be
analytic (Gaussian or Moffat, with the width, β and ellipticity fitted on bright
isolated stars or given), a field-dependent PSF file (the same format that
`artificial --psf-file` uses to simulate frames), or an empirical PSF built from the
stars of the frame: an analytic core or a PSF file plus residual tables, which
can vary over the field. Each star gets its flux and position (and, if you ask for it, a
local sky) from a weighted least-squares fit. The weights come from the noise
model of the detector. Fitted neighbours are subtracted, close stars are fitted
simultaneously in groups, and stars that `find` missed can be detected in
the residuals. `PSFCOUNT` is the total flux of the star (the integral of the PSF), the
same quantity as `GCOUNT` of `gphot` or an aperture-corrected `APCOUNT`.
[`phcal`](phcal.md) `--photometry psf` calibrates it.

## Synopsis

```bash
munipack psfphot [options] file(s)
```

The frames must have the `FIND` extension (`find`). The `APERPHOT` extension
([`aphot`](aphot.md)) provides the sky annulus and the starting sky of the stars. Without it,
`psfphot` uses an annulus that starts 2 px beyond the fitting radius and is at
least 5 px wide. The detector parameters come from the header: `SATURATE`
(default 2<sup>BITPIX</sup>−1 for integer frames), `GAIN` [e⁻/ADU] (1 when it is missing or
`BUNIT = 'COUNTS'`), `READNS` [ADU] or `RDNOISE` [e⁻], and `BIASLEV` [ADU]
(0 when missing). The keyword names can be redefined by `FITS_KEY_SATURATE`,
`FITS_KEY_GAIN` and `FITS_KEY_READNOISE`.

In a container:

```bash
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes \
    psfphot --verbose frame.fits
```

## Options

### PSF model

| Option | Argument | Default | Description |
|---|---|---|---|
| `--psf` | `GAUSS`, `MOFFAT`, `FILE`, `EMPIRICAL` | `GAUSS` (`FILE` when `--psf-file` is given) | The PSF type (see PSF types below). |
| `--psf-file` | file | — | The PSF file (the format of `artificial --psf-file`). With `EMPIRICAL`, this file is the base that the frame corrects (the correction mode). |
| `--fwhm` | pix | fitted | FWHM of the analytic PSF. If you omit it, the shape is fitted on up to 30 bright (S/N > 50), isolated, unsaturated stars (medians of their fits). |
| `--beta` | number | fitted | The Moffat exponent β. It is 2.5 when `--fwhm` is given and `--beta` is not. A fitted β above 1000 is treated as the Gaussian limit. |
| `--elliptic` | — | off | Also fits the axis ratio and the position angle of the analytic PSF. |
| `--scale` | pix/deg | header (WCS) | The plate scale used to place a star in a PSF file. It is in pixels per degree. |
| `--psf-centre` | `x,y` [pix] | `NAXIS/2` | The field centre of the PSF file (integer division, the same convention as `artificial`). |

### Empirical PSF (`--psf EMPIRICAL`)

| Option | Argument | Default | Description |
|---|---|---|---|
| `--psf-core` | `MOFFAT`, `GAUSS` | `MOFFAT` | The analytic core. Its shape and ellipticity are fitted, and β is capped at 10. The Moffat wings carry the light beyond the tables, so a Gaussian core on stars with Moffat wings loses that light. Ignored when `--psf-file` is given. |
| `--psf-zoom` | n | 4 for FWHM < 5 px, else 2 | The number of samples per pixel of the tables. In the correction mode the `ZOOM` of the file is used. |
| `--psf-radius` | pix | 99.5% of the flux, ≥ 4 FWHM | The radius of the residual tables. In the correction mode, the whole stamps of the file. |
| `--psf-var` | 0, 1, 2 | 0 | The order of the polynomial variation in x, y over the field (1, 3 or 6 terms; DAOPHOT `VARIABLE`). |
| `--psf-var-radial` | 0–3 | 0 (1 in the correction mode of a field-dependent file with ≥ 20 PSF stars) | Adds the powers ρ, ρ², ρ³ of the field radius, measured from the PSF-file centre and normalised to the farthest corner. |
| `--psf-stars` | n | 100 | The maximum number of PSF stars. |
| `--psf-list` | file | automatic | The PSF stars: `x y` per line (`#` starts a comment), or a FITS table with columns `X`, `Y`. The nearest star within 2 px is used. |
| `--psf-save` | file | — | Writes the empirical PSF as a PSF file (it is overwritten). Single input file only. |
| `--psf-smooth` | number | cross-validation | The regularisation of the tables, relative to the local information of each sample. By default it is chosen by 4-fold cross-validation (0.01 with fewer than 16 PSF stars). |

### Fit

| Option | Argument | Default | Description |
|---|---|---|---|
| `--fitrad` | pix | 1.5 × FWHM, min. 3 | The fitting radius (of the local PSF). |
| `--sky` | `fixed`, `local` | `fixed` | `fixed`: the clipped mean in the `aphot` annulus, re-estimated before each later pass on the frame with the fitted stars subtracted. `local`: a constant sky fitted for each star (or for each group). |
| `--maxiter` | n | 50 | The maximum number of iterations of one fit. |
| `--npass` | n | 3 | The number of passes over all stars. |
| `--group` | FWHM | 1.5 | Stars closer than this distance are fitted simultaneously. `0` turns grouping off (iterative subtraction only). |
| `--saturate` | ADU | header | The saturation level. Pixels at or above it are excluded. |
| `--gain` | e⁻/ADU | header | The gain. |
| `--read-noise` | ADU | header | The read noise. |

### Residual detection

| Option | Argument | Default | Description |
|---|---|---|---|
| `--residual-detect` | n | 0 | The number of passes that detect stars in the residual frame. Each pass adds the stars it finds (FLAG 128) and refits all stars (two more passes). |
| `--residual-threshold` | σ | by false-alarm rate | The detection threshold. The default is the `find` threshold (`THRESH`), but at least the Normal deviate for a false-alarm probability of 1/N<sub>pix</sub>: 4.75 for 1024², 5.0 for 2048², 5.3 for 4096². It is recorded as `RESTHRES`. |
| `--residual` | file | — | Saves the frame with all fitted stars subtracted (the sky is kept). Single input file only. |

### Common options

| Option | Argument | Default | Description |
|---|---|---|---|
| `-t`, `--target-directory` | dir | — | The directory for the output files. |
| `-O` | — | off | Turns on the advanced processing of output filenames. |
| `--pattern` | str | — | The matching pattern (for `-O`). |
| `--mask` | str | — | The mask for the output file(s) (for `-O`). |
| `--format` | str | — | The format of the output files (for `-O`). |
| `--verbose` | — | off | Verbose log: the fitted PSF shape, the PSF stars and the groups. |
| `-h`, `--help` | — | — | Shows the help. |

## PSF types

All PSFs are sampled with `ZOOM` samples per pixel (16 for analytic PSFs, the
file's `ZOOM` for a PSF file). A pixel value is the integral of the
interpolated samples over the pixel at the star's sub-pixel position, so the
flux is conserved exactly. `artificial` uses the same code, so a frame
simulated with a PSF is fitted by exactly that PSF.

| Type | Profile / source |
|---|---|
| `GAUSS` | exp(−x²/2σ² − y²/2σ²q²), FWHM = 2.3548 σ √q. |
| `MOFFAT` | (1 + x²/α² + y²/α²q²)<sup>−β</sup>. It is fitted as FWHM and u = 1/β, so Gaussian stars converge (u → 0). When the median β exceeds 1000, the Gaussian of the fitted FWHM is used (`PSFLIMIT = 'GAUSS'`). |
| `FILE` | A field-dependent PSF file: extensions with `RADCOO`, `ANGCOO`, interpolated linearly. The field position is the pixel offset from the centre (`--psf-centre`) divided by the scale (`--scale`). |
| `EMPIRICAL` | A base plus residual tables built from the frame's stars, see below. |

Here q is the axis ratio and x is along the major axis. The `artificial` options
`--hwhm` (σ for GAUSS, α for MOFFAT), `--eccentricity` (q = √(1 − e²)) and
`--inclination` describe the same functions.

### Empirical PSF

The model is PSF(x, y) = B + Σ<sub>t</sub> P<sub>t</sub>(x, y) R<sub>t</sub>, in the manner of
DAOPHOT. B is the base: either the analytic core (`--psf-core`), or, in the
**correction mode** (`--psf EMPIRICAL --psf-file lab.fits`), a PSF file. The
R<sub>t</sub> are residual tables, and the P<sub>t</sub> are the terms of the
variation: the polynomial in normalised x, y of order `--psf-var` (1, x, y, x², xy, y²),
followed by the radial powers ρ…ρ<sup>n</sup> up to `--psf-var-radial`.

- **PSF stars.** They come from `--psf-list`, or are the brightest stars (S/N ≥ 30, up to
  `--psf-stars`) that are unsaturated, free of bad pixels and fully inside the frame,
  with no brighter neighbour within the table radius and no neighbour above 10% within 3 FWHM. The variation needs three stars per term
  (order 1: 9, order 2: 18, plus three for each radial power). With fewer stars, the
  order is lowered and a warning is printed.
- **Fit of the tables.** A weighted linear least-squares fit with the exact
  pixel-integration model. Neighbours are subtracted and the sky of each PSF star is fitted on a
  ring outside the table. Pixels beyond 5σ are clipped. There are three iterations of
  fitting all stars and rebuilding the tables. The tables are constrained so that
  they do not move the centroid.
- **Flux normalisation.** The tables carry the light beyond their radius,
  measured in a large aperture (15–40 px, `PSFNAPER`) on isolated PSF stars, and
  the result is applied only when it is significant (3σ). As a result, the light beyond the
  tables is not taken from the wings of the analytic core.
- **Correction mode.** The correction shrinks towards zero in radial zones where
  the subsets of PSF stars do not show a significant change, so a
  correct base PSF is left as it is.
- **`--psf-save`** writes a PSF file: one stamp for a constant PSF with an
  analytic core, or rings × angles over the frame for a varying PSF or a
  correction. `psfphot --psf-file` and `artificial --psf-file` read the file back.

Limitations: a low-order polynomial in x, y cannot follow a PSF whose structure
changes over the field, such as a sharp core mixed with a wide halo. When the change follows the field
radius, use `--psf-var-radial`. Otherwise, use the correction mode with a model
PSF. The normalisation needs isolated bright stars, and the PSF stars must
be mostly single stars.

## Grouping and residual detection

Stars are fitted from the brightest down. The models of stars that are already fitted are subtracted,
and later passes refit each star with all the others subtracted. Stars that are linked by
distances below `--group` × FWHM, including chains of close stars, form a group. A group's members are fitted
together: all their fluxes and positions, plus a common sky with `--sky local`, over the union of their
fitting circles. The groups are formed again before each pass. A group with more than 15 stars is split
by shorter links. Iterative subtraction alone biases pairs at 1 FWHM by tens
of per cent, and the simultaneous fit removes that bias.

`--residual-detect n` searches the residual frame with the `find` detection, run
against the model noise plus the measured PSF mismatch. A candidate must also be
on a local maximum, at least 0.75 FWHM from every star, and off saturated areas. The stars it finds are added
and all stars are refitted, n times. `find` and the residual detection cannot separate stars
closer than about 0.5 FWHM.

Saturated stars are not in `FIND`. `psfphot` adds each connected saturated area that
has no detected star as a star (FLAG 64). That star is fitted from its wings, so its light is
subtracted from its neighbours.

## Outputs

`psfphot` adds a `PSFPHOT` binary table to each input frame. It has one row for each `FIND` row, in the
same order, followed by the added saturated stars and then the stars found in the residuals.

| Column | Unit | Description |
|---|---|---|
| `X`, `Y`, `XERR`, `YERR` | pix | The fitted position and its errors. The errors are −1 when the position is not fitted. |
| `PSFCOUNT`, `PSFCOUNTERR` | ADU | The total flux and its error. The error comes from the covariance of the fit, with the error of the fixed sky propagated, and it is −1 when the fit failed. With `EMPIRICAL`, it includes `PSFMODERR`. |
| `SKY`, `SKYERR` | ADU | The sky and its error. |
| `CHI2` | — | The reduced χ² of the fit (≈ 1 for a good fit). |
| `SHARP` | — | The relative residual inside 0.5 FWHM: positive for sources sharper than the PSF (cosmic rays, hot pixels), negative for extended ones (galaxies, blends). |
| `ROUND` | — | The ellipticity of the data minus that of the model. |
| `NPIX` | — | The number of fitted pixels. |
| `FLAG` | — | The bit flags, see below. |
| `GROUP`, `GSIZE` | — | The group of the simultaneous fit in the last pass (0 when the star was fitted alone) and its size. |
| `PSFMODERR` | ADU | **`EMPIRICAL` only**, after `GSIZE`. The flux error from the uncertainty of the PSF model: the split-sample variance over 8 subsets of the PSF stars, plus the error and the systematic scatter of the normalisation. It is included in `PSFCOUNTERR` in quadrature. |

`FLAG` bits:

| Bit | Meaning |
|---|---|
| 1 | Saturated pixel(s) in the core. The flux comes from the wings only. |
| 2 | Pixels within the fitting radius were excluded (bad, saturated, outside the frame). |
| 4 | The position was not fitted: S/N < 3 at the detected position, or the fit moved by more than half the fitting radius (at least 1 px). The flux is fitted at the detected position. |
| 8 | The fit did not converge. |
| 16 | The fit failed. The flux is not valid. |
| 32 | Another star lies within the fitting radius. |
| 64 | A saturated star that `find` did not detect, added by `psfphot`. |
| 128 | A star found in the residuals (`--residual-detect`). |
| 256 | S/N = `PSFCOUNT`/`PSFCOUNTERR` < 2. Treat the flux as an upper-limit-like estimate. |

```{note}
`FLAG` is a sum of bits. Test the individual bits rather than `FLAG != 0`: bit 32
(crowding) and bit 256 are common and do not invalidate a measurement.
`phcal` selects stars by their relative error (`-e`), not by `FLAG`.
```

Header keywords of `PSFPHOT`:

- PSF: `PSFTYPE`, `PSFFILE`, `FWHM`, `HWHM` [pix], `PSFSIGMA` or `PSFALPHA` and
  `PSFBETA` (`PSFLIMIT = 'GAUSS'` for the Gaussian limit of a Moffat PSF), `PSFAXRAT`, `PSFPA` [deg],
  `PSFNSTAR`, `PSFCEN1`, `PSFCEN2`, `PSFSCALE`, `PSFZOOM`.
- Fit: `FITRAD`, `SKYMODE`, `NPASS`, `GROUPDST` [FWHM], `NGROUPS`, `MAXGROUP`,
  `RESDETEC`, `NRESID`, `RESTHRES` [σ]; the detector values used are `GAIN`, `READNS`,
  `BIASLEV`, `SATURATE`.
- Empirical PSF: `PSFCORE` (`GAUSS`, `MOFFAT`, `FILE`), `PSFVAR`, `PSFVARR`,
  `PSFNSHAP`, `PSFRAD`, `PSFR995`, `PSFRESRM`, `PSFCHI2`, `PSFITER`, `PSFSMOOT`,
  `PSFNAPER`, `PSFNORMN`, `PSFNORM`, `PSFNORME`, `PSFSYSER`, `PSFSAVED`. Here `FWHM`
  is the empirical PSF at the field centre, and the shape keywords describe the core.

`--residual` writes the star-subtracted image. `--psf-save` writes a PSF file
whose primary header has `ZOOM`, `PSFCEN1`, `PSFCEN2`, `PSFSCALE`, `FWHM` and `PSFRAD`. Without a scale, `RADCOO` is in
pixels, and you should use the file with `--scale 1`.

## Examples

A Gaussian PSF with the width fitted, then calibration:

```bash
munipack find frame.fits
munipack aphot frame.fits
munipack psfphot --verbose frame.fits
munipack phcal --photometry psf -c cat.fits --col-ra RAJ2000 --col-dec DEJ2000 \
    --col-mag Vmag --col-magerr e_Vmag --photsys-ref Johnson -f V frame.fits
```

An elliptical Moffat PSF with a local sky and wider groups, keeping the residual frame:

```bash
munipack psfphot --psf MOFFAT --elliptic --sky local --group 2 \
    --residual resid.fits frame.fits
```

An empirical PSF that varies to second order over the field, saved and then reused on
other frames (or by `artificial`):

```bash
munipack psfphot --psf EMPIRICAL --psf-var 2 --psf-save epsf.fits --verbose frame.fits
munipack psfphot --psf-file epsf.fits other.fits
```

Wide-field optics: a PSF that changes with the field radius plus a linear tilt. A second command
corrects a laboratory PSF by the frame and searches for the stars that `find` missed:

```bash
munipack psfphot --psf EMPIRICAL --psf-var 1 --psf-var-radial 1 frame.fits
munipack psfphot --psf EMPIRICAL --psf-file lab.fits --residual-detect 1 frame.fits
```
