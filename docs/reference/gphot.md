# gphot — growth-curve photometry

`gphot` measures the total flux of every star in a frame from the frame's
growth curve. The growth curve is the enclosed energy EE<sub>k</sub> of
aperture *k*: the fraction of the total flux that falls inside it. It is
built the same way as the automatic aperture correction of
[`aphot --apcorr auto`](aphot.md). The curve uses bright (S/N > 30, or at
least the 30 brightest), unsaturated stars that have no neighbour brighter
than 1% within 5.5 FWHM.
Their total flux is the mean of three plateau apertures, and EE can vary over
the field as a polynomial in x,y or in the distance from the field centre.
For each star the counts are modelled as C<sub>k</sub> = GCOUNT · EE<sub>k</sub> + noise,
with the covariance of nested apertures. GCOUNT is the estimate with the
least variance among the single apertures C<sub>k</sub>/EE<sub>k</sub> and the best
linear combinations of apertures 1…m. Faint stars therefore get
profile-weighted apertures, with the scatter of PSF photometry, and bright
stars get the large apertures. GCOUNTERR includes the error of the curve.
`phcal` calibrates GCOUNT by default when the frame has a `GROWPHOT`
extension (`--photometry grow`).

## Synopsis

```bash
munipack gphot [options] file(s)
```

The input frames need:

- the `APERPHOT` extension written by `aphot`, with many apertures from about
  1 FWHM to several FWHM, all inside the sky ring;
- the FWHM of the `FIND` extension written by `find`, which sets the isolation
  radius of the curve stars (without it, twice the HWHM of `APERPHOT` is used).

For example:

```bash
munipack aphot --apertures 1,1.313,1.7241,2.2637,2.9724,3.9028,5.1245,6.7287,8.8349,11.601,15.232,20,23,26.4,30.4,35,40 \
               --ring 45,60 frame.fits
```

By default the file is updated in place. To write the result to a new file
instead, give the input and output names as one argument separated by a
space, `'frame.fits result.fits'` (quoted), or use `-t`. The comma form
`file,result` shown by `--help` is not split and fails.

## Options

### Growth curve

| Option | Argument | Default | Description |
|---|---|---|---|
| `--var` | 0, 1, 2 | 2 | Highest order of the polynomial for how the growth curve varies over the field. The fit starts at order 0 and moves up only while the higher terms are significant. With few stars the order is lowered. Not with `--legacy`. |
| `--radial` | — | off (x,y) | Model the variation as a function of distance from the field centre instead of x,y. Use it for PSFs that change with field radius. Not with `--legacy`. |
| `--centre` | `x,y` [pix] | NAXIS/2 | Field centre used by `--radial`; requires `--radial`. |
| `--legacy` | — | off | Use the growth-curve fit of versions before 0.6.2+hg1653: a profile model fitted to stars selected by HWHM. Gives the same results as those versions, bit for bit. |
| `-th`, `--threshold` | real | 0.1 | Only with `--legacy` (an error otherwise): curve stars must be brighter than this many times the sky in their aperture. |

### Output files

| Option | Argument | Default | Description |
|---|---|---|---|
| `-t`, `--target-directory` | directory | — | Write the output files into this directory (it must exist) instead of updating the input files. |
| `-O` | — | off | Turn on advanced output filenames (`--pattern`, `--mask`, `--format`). |
| `--pattern` | regex | `(.+)\.(.+)` | Regular expression matched against input filenames (with `-O`). |
| `--mask` | string | — | Output filename built from back-references `\0`, `\1`, … (with `-O`). |
| `--format` | printf format | — | Format for numbered output files, e.g. `out_%02d.fits` (with `-O`). |

### General

| Option | Argument | Default | Description |
|---|---|---|---|
| `--verbose` | — | off | Print the fit of the growth curve (EE of every aperture, the plateau, HWHM, the 90% radius) and, for each aperture, how many stars took GCOUNT up to it. |
| `-h`, `--help` | — | — | Show the help and exit. |

## Outputs

`gphot` adds two binary-table extensions to the frame. If either already
exists, it is replaced.

**`GROWPHOT`**, one row per star of `APERPHOT`:

| Column | Description |
|---|---|
| `X`, `Y` | Position [pix] |
| `SKY` | Sky from `APERPHOT` plus the residual sky *b* that the star's apertures estimate (prior 0 ± SKYERR) |
| `SKYERR` | Posterior error of SKY |
| `GCOUNT` | Total flux [ADU] within the plateau apertures |
| `GCOUNTERR` | Standard deviation of GCOUNT, including the error of the curve |
| `GROWFLAG` | 0: star not used for the curve; 1: curve star; 2: non-stellar object (`--legacy` only) |

Stars without a valid GCOUNT have GCOUNT = GCOUNTERR = −1 (their SKY and
SKYERR are those of `APERPHOT`).

Header keys of `GROWPHOT`:

- `HWHM`, `RADFLX90`: half width at half maximum of the curve, and the radius
  that contains 90% of the flux [pix];
- `GROWMETH`: `EE` (enclosed energy, the default) or `LEGACY`;
- the `EE*` keys, as in `APERPHOT` of `aphot --apcorr auto` (not written with `--legacy`):
  - `EESOURCE`, `EENSTAR`: source of the curve and the number of curve stars;
  - `EEVARREQ`, `EEVAR`, `EEVARCOO`, `EEFCX`, `EEFCY`: requested and fitted
    order of the variation, x,y or radial, and the field centre (`EEFCX`,
    `EEFCY` with `--radial` only);
  - `EERTOT1`, `EERTOT2`, `EEPLAT`, `EEPLATZ`, `EEPLATB`, `EEPLATBE`, `EER99`,
    `EERCAP`, `EESKYBIA`: the plateau apertures that define the total flux,
    how the plateau was found, the residual sky of the plateau, and the radius
    holding 99% of the stacked curve with its search cap;
  - `EEORDk`, `EESIGk`, `EELAMk`: for each aperture *k*, the order used, the
    scatter of the stars about the model, and the inflation of the model
    variance from the leave-one-out residuals.

**`GROWCURVE`**, one row per aperture:

| Column | Description |
|---|---|
| `R` | Aperture radius [pix] |
| `GROWCURVE` | Median EE of the stars; −1 for apertures without EE |
| `GROWCURVEERR` | Error of GROWCURVE |
| `RADIALPROFILE` | Mean surface brightness of the annuli, relative to the total flux |

```{note}
GCOUNT is the flux within the plateau apertures (EERTOT1…EERTOT2). This is
the same total that `aphot` gives as APCOUNT/APEE. It is not a profile
extrapolated to infinity, so on frames with extended PSF wings it is lower
by the light beyond the plateau. That fraction is constant and the
photometric calibration absorbs it.
```

```{note}
`gphot` needs at least 3 bright isolated stars and fails with an error if it
has fewer. Shallow frames with few stars often do not show a significant
variation of the PSF; the errors then account for it. In crowded fields SKY
includes the light of faint neighbours (about +0.3 SKYERR).
```

## Examples

The full chain from detection to calibration, with the curve varying in x,y:

```bash
munipack find frame.fits
munipack aphot --apertures 1,1.313,1.7241,2.2637,2.9724,3.9028,5.1245,6.7287,8.8349,11.601,15.232,20,23,26.4,30.4,35,40 \
               --ring 45,60 frame.fits
munipack gphot frame.fits
munipack phcal --photometry grow --photsys-ref Johnson -f V \
               --col-ra RAJ2000 --col-dec DEJ2000 --col-mag Vmag --col-magerr e_Vmag \
               -c catalogue.fits frame.fits
```

`phcal` writes the calibrated stars to `frame_cal.fits`.

A PSF that changes with field radius (for example, wide-field optics), with
the centre of the optical axis given explicitly:

```bash
munipack gphot --radial --centre 1024,1024 frame.fits
```

A constant growth curve, the result written to a new file, and the curve
printed:

```bash
munipack gphot --var 0 --verbose 'frame.fits frame_grow.fits'
```

In the container, many frames, results written to `grow/`:

```bash
mkdir -p grow
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes \
       gphot -t grow frame_*.fits
```
