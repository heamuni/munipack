# aphot — aperture photometry

`aphot` sums the counts of every star in a set of circular (optionally
elliptic) apertures and subtracts the sky. The sky is measured in an
annulus around each star. The stars are those that [`find`](find.md) stored
in the `FIND` extension of the frame. The results go into a new `APERPHOT`
table extension of the same file: one row per star, and for every aperture
the sky-subtracted counts and their error, which follows the noise model of
the detector. Optionally `aphot` also records the enclosed energy (APEE) of
every aperture: the fraction of the star's total flux that falls inside it.
The APEE comes either from a known PSF or from the growth curves of the
frame's own stars. [`phcal`](phcal.md) then calibrates the total flux
APCOUNT/APEE, so a PSF that changes over the field does not change the
zero point. The `APERPHOT` table is also the input of [`gphot`](gphot.md),
and its sky ring is used by [`psfphot`](psfphot.md).

## Synopsis

```bash
munipack aphot [options] file(s)
```

In the container:

```bash
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes aphot frame.fits
```

The frame must already have a `FIND` extension (run `munipack find frame.fits`
first). By default the file is updated in place. Use `-t dir` (or `-O` with
`--mask` or `--format`) to write the result to another file instead; the
`file,result` form of some other actions is not accepted by `aphot`. A
previous `APERPHOT` extension is replaced.

## Options

### Apertures and sky ring

| Option | Argument | Default | Description |
|---|---|---|---|
| `--apertures` | `r1,r2,...` [pix] | 12 radii from 1 to 20 px | Aperture radii in pixels, in increasing order, separated by commas. There can be any number of them. The default is a logarithmic sequence 20<sup>(i−1)/11</sup>: 1, 1.313, 1.724, 2.264, 2.972, 3.903, 5.125, 6.729, 8.835, 11.60, 15.23, 20. |
| `--ring` | `ri,ro` [pix] | `20,30` | Inner and outer radius of the sky annulus in pixels. The inner radius should be at least as large as the largest aperture. |
| `--enable-ellipticity` | — | off (circles) | Use elliptic apertures. Their shape is ECCENTR and INCL of the `FIND` extension, unless `--eccentricity`/`--inclination` override them. Experimental. In practice it gains little even for trailed stars. |
| `--eccentricity` | real, 0 ≤ e < 1 | ECCENTR of `FIND` | Eccentricity of the elliptic apertures. Only takes effect with `--enable-ellipticity` (without it the apertures are circles). |
| `--inclination` | real [deg] | INCL of `FIND` | Position angle of the major semi-axis, counter-clockwise from the x axis. Only takes effect with `--enable-ellipticity`. |
| `--coo` | `"x1,y1 x2,y2 ..."` [pix] | — | Measure the given pixel positions instead of the `FIND` stars. Positions are separated by spaces (or `\|`). The frame needs no `FIND` extension. The results are printed to standard output and the file is left unchanged (see "`--coo` output" below). Cannot be combined with the APEE options. |

### Enclosed energy from a PSF

Give one of these options to write APEE columns computed from a known PSF.
They cannot be combined with `--apcorr`.

| Option | Argument | Default | Description |
|---|---|---|---|
| `--psf-file` | file | — | Compute APEE exactly from a field-dependent PSF file, in the format used by `artificial --psf-file` and `psfphot --psf-file`. Implies `--psf FILE`; combined with `--psf GAUSS` or `MOFFAT` it is an error. |
| `--psf` | `GAUSS`, `MOFFAT`, `FILE` | — | Compute APEE from an analytic Gaussian or Moffat PSF, or from the PSF file (`FILE` requires `--psf-file`). |
| `--fwhm` | real [pix] | FWHM of `FIND` | FWHM of the GAUSS or MOFFAT PSF. If `FIND` has no FWHM, twice the HWHM measured by `aphot` is used. |
| `--beta` | real > 1 | 2.5 | Moffat exponent β. |
| `--scale` | real [pix/deg] | WCS of the header | Plate scale that converts pixel offsets into the field positions (RADCOO, degrees) of the PSF file. Required for a PSF file with more than one ring when the frame has no WCS. |
| `--psf-centre` | `x,y` [pix] | NAXIS/2 | Field centre of the PSF file, and also the centre used by `--apcorr-radial`. |

### Enclosed energy from the frame's stars

| Option | Argument | Default | Description |
|---|---|---|---|
| `--apcorr` | `auto` | — | Compute APEE from the growth curves of the bright, isolated stars of the frame. The only accepted value is `auto`. |
| `--apcorr-var` | 0, 1, 2 | 2 | Highest order of the polynomial that describes how APEE varies over the field. The order starts at 0 (constant) and goes up only while the higher terms are significant. Requires `--apcorr auto`. |
| `--apcorr-radial` | — | off (x,y) | Model the variation as a function of the distance from the field centre (`--psf-centre`, default NAXIS/2) instead of x,y. Use it for optics whose PSF changes with field radius. Requires `--apcorr auto`. |

### Output files and logging

| Option | Argument | Default | Description |
|---|---|---|---|
| `-t`, `--target-directory` | directory | — | Write the output files into this directory and leave the inputs unchanged. |
| `-O` | — | off | Turn on advanced output filenames (`--pattern`, `--mask`, `--format`). |
| `--pattern` | regex | `(.+)\.(.+)` | Regular expression matched against input filenames (with `-O`). |
| `--mask` | string | — | Output filename built from back-references `\0`, `\1`, … (with `-O`). |
| `--format` | printf format | — | Numbered output filenames, e.g. `out_%02d.fits` (with `-O`). |
| `--verbose` | — | off | Print progress: the number of stars, the number of apertures invalidated by neighbours, and the aperture-correction fit. |
| `-h`, `--help` | — | — | Show the option summary. |

## Method

### Sky

The aperture contains the sky pixels as well as the star, so the value
subtracted must be the **mean** sky level. A mode, median or robust mean
would bias the counts by the difference times the aperture area. The bias
is large on a sky of a few photo-electrons per pixel, whose distribution is
skewed and discrete. For each star the sky is computed as follows:

- **Pixels used:** the pixels of the ring `ri…ro`, leaving out the masks of
  the other `FIND` stars. A mask follows the star's light outwards, up to
  half the outer ring radius. If the masks would leave less than a quarter
  of the ring (crowded fields), the whole ring is used.
- **Clipped mean:** the mean of these pixels, iterated from the median and
  clipped to *t* − 5σ … *t* + *q* + √(*q*² + (5σ)²). Here σ is the larger
  of the pixel scatter and the noise model √((*t* − BIASLEV)/GAIN + READNS²)
  [ADU], and *q* = 25/6 photo-electrons. Read noise comes from READNS
  [ADU] or RDNOISE [e⁻]. With this asymmetric upper limit no sky photon is
  clipped, while stars, cosmic rays and hot pixels are. Without GAIN the
  pixel scatter alone is used.
- **SKYERR:** the scatter of the used pixels divided by the square root of
  their number.

The same estimator is used for the sky of `find` and `psfphot`, the default
stack of `kombine` and the master of `dark`.

### Counts, errors and invalid apertures

APCOUNT is the sum over the aperture minus SKY times the area
*A* = Σ*f*. Here *f* is the fraction of each pixel inside the aperture, so
partial pixels on the edge are weighted. The error is

```text
APCOUNTERR² = APCOUNT/GAIN + Σf² σ² + A² SKYERR²
```

The three terms are the Poisson noise of the star, the noise of the sky
pixels in the aperture (σ is the scatter of the ring, which includes the
read noise), and the error of the mean sky. GAIN is taken from the frame. It
is 1 for frames in electrons (`phcorr` output). If GAIN is missing, 1 e⁻/ADU
is assumed and a warning is printed. The gain actually used is recorded as
`APERRGN`.

An aperture that cannot be measured gets the value −1. This happens when the
aperture runs off the frame, contains saturated or bad pixels
(outside LOBAD…HIBAD), or has negative sky-subtracted counts. An aperture
is also set to −1 when it is contaminated by a brighter neighbour, that is,
when the neighbour's light inside the aperture exceeds three times the
aperture's error and is itself significant (three times its own error). The
neighbour's light is its flux times the mean profile of the frame (the
median growth curve of the isolated bright stars); with fewer than three
such stars the neighbour's own growth curve is used. A neighbour without
valid apertures (saturated, at the border) invalidates the apertures that
reach its mask. This drops faint stars next to bright ones, and duplicate
detections of the same star, from the calibration.

### Enclosed energy (aperture corrections)

APEE is written only when `--psf-file`, `--psf` or `--apcorr auto` is
given. Without these options the `APERPHOT` table is the same as before.

- **From a PSF** (`--psf-file`, `--psf GAUSS|MOFFAT`): the PSF at the
  star's field position is pixelised at its sub-pixel centre and summed with
  exactly the partial-pixel weights of the apertures. On noiseless frames
  simulated with the same PSF, APCOUNT/APEE equals the true flux to 10⁻⁶ in
  every aperture. APEEERR is zero. For a PSF file, the field position is
  the pixel offset from NAXIS/2 (or `--psf-centre`) divided by the scale
  (WCS or `--scale`). This is the same convention as `artificial` and
  `psfphot`.
- **From the frame** (`--apcorr auto`):
  - **Stars:** bright stars (S/N > 30 in the aperture of 2 FWHM) that are
    unsaturated, have all apertures valid, and have no neighbour brighter
    than 1% within 5.5 FWHM.
  - **Total flux T:** the mean of three consecutive plateau apertures.
    The plateau is the first such window, starting above 1.3 FWHM, whose
    growth-curve slope is consistent with zero within 2σ. The search stops at twice
    the radius that encloses 99% of the stacked growth curve
    (EEPLAT = SLOPE). If no window passes, the last apertures within that
    cap are used (EEPLAT = CAP).
  - **Field model:** APEE<sub>k</sub> = C<sub>k</sub>/T of each aperture is
    fitted by a weighted, robust polynomial over the field, of order up to
    `--apcorr-var`, in x,y or in field radius. The errors take the covariance
    of the nested apertures into account. The order is raised from the
    constant only when the higher terms are significant (99.9%), and it is
    lowered for few stars with a warning. The model is not extrapolated
    beyond the stars.
  - **APEEERR:** the model error combined in quadrature with the excess
    scatter of the stars.
  - **Failure:** if the fit fails (too few suitable stars), a warning is
    printed and no APEE is saved.

For `--apcorr auto`, include apertures out to the sky ring, so that the
plateau can be found, for example `--apertures ...,20,23,26.4,30.4,35,40 --ring 45,60`.
The method needs bright, isolated stars (about 30), and the FWHM in `FIND`
must be right.

## Outputs

### The `APERPHOT` extension

One row per `FIND` star, in the same order. All columns are double
precision. *n* is the number of apertures.

| Column | Unit | Description |
|---|---|---|
| `X`, `Y` | pix | Position of the star (from `FIND`). |
| `SKY` | ADU/pix | Mean sky in the ring. |
| `SKYERR` | ADU/pix | Error of SKY. |
| `APCOUNT1` … `APCOUNTn` | ADU | Sky-subtracted counts in each aperture; −1 if invalid. |
| `APCOUNTERR1` … `APCOUNTERRn` | ADU | Error of APCOUNT; −1 if invalid. |
| `APEE1` … `APEEn` | — | Enclosed energy of each aperture (only with the APEE options). These columns follow all the APCOUNT columns. |
| `APEEERR1` … `APEEERRn` | — | Error of APEE (zero for a PSF). |

Header keywords:

| Keyword | Description |
|---|---|
| `HWHM` | Half width at half maximum [pix], estimated from the second moments of the brightest stars (a true HWHM; versions before 0.6.2+hg1653 wrote σ). `phcal` uses it to select the aperture. |
| `ECCENTR`, `INCL` | Eccentricity and inclination [deg] of the apertures (0 for circles). |
| `NAPER`, `APER1` … `APERn` | Number of apertures and their radii [pix]. |
| `ANNULUS1`, `ANNULUS2` | Inner and outer sky ring radius [pix]. |
| `APERRGN` | Gain [e⁻/ADU] used for the star term of APCOUNTERR. |
| `EESOURCE` | Where APEE came from: `PSF` or `AUTO`. |
| `EEPSF`, `EEPSFFIL`, `EEPSFCX`, `EEPSFCY`, `EEPSFSCL` | From a PSF: the PSF type and file, the field centre [pix] and the scale [pix/deg]. |
| `EEFWHM`, `EEBETA` | From an analytic PSF: FWHM [pix] and Moffat β. |
| `EEVARREQ`, `EEVAR`, `EEVARCOO` | AUTO: requested order of the variation, the order allowed by the number of stars (at most EEVARREQ; the order of each aperture is EEORD*k*), and `XY` or `RADIAL`. |
| `EEFCX`, `EEFCY` | AUTO, radial: field centre [pix]. |
| `EENSTAR` | AUTO: number of stars used. |
| `EERTOT1`, `EERTOT2` | AUTO: smallest and largest plateau aperture [pix]. |
| `EEPLAT`, `EEPLATZ`, `EEPLATB`, `EEPLATBE` | AUTO: how the plateau was found (`SLOPE` or `CAP`), its slope in units of its error, and the residual sky of the plateau and its error [ADU/pix]. The residual sky is not removed. |
| `EER99`, `EERCAP` | AUTO: radius enclosing 99% of the stacked growth curve, and the search limit 2·EER99 [pix]. |
| `EESKYBIA` | AUTO: residual sky slope beyond 4 FWHM [ADU/pix]. A diagnostic only; it is not removed. |
| `EEORDk`, `EESIGk` | AUTO: significant order of the variation, and scatter of the stars over the model, for aperture *k*. |

### `--coo` output

With `--coo`, nothing is written to the file. One line per position is
printed to standard output:

```text
 =BYHAND>  i  x  y  SKY  SKYERR  APCOUNT1 ... APCOUNTn  APCOUNTERR1 ... APCOUNTERRn
```

## Examples

Detect stars and measure them in the default apertures:

```bash
munipack find frame.fits
munipack aphot frame.fits
```

Choose the apertures and the sky ring yourself, and write the result to a
separate directory:

```bash
munipack aphot --apertures 2,3,4,6,8,12 --ring 15,25 -t phot/ frame.fits
```

Aperture corrections from a field-dependent PSF file, from an analytic
Gaussian PSF, or from the frame's own stars with a radial variation. `phcal`
then calibrates the total flux (the column names are those of the
catalogue):

```bash
munipack aphot --psf-file psf.fits frame.fits
munipack aphot --psf GAUSS --fwhm 4.9 frame.fits
munipack aphot --apertures 1,1.5,2,3,4,5,7,9,12,15,20,23,26.4,30.4,35,40 --ring 45,60 \
               --apcorr auto --apcorr-radial frame.fits
munipack phcal --saper 7 -c cone.fits --photsys-ref Johnson -f V \
               --col-mag Vmag --col-magerr e_Vmag frame.fits
```

Counts at given positions, without modifying the file:

```bash
munipack aphot --coo "666,666 120.5,340.2" frame.fits
```

```{note}
`phcal` uses APEE automatically for aperture photometry whenever the
columns are present (`phcal --apcorr none` ignores them). If `gphot` has
added a `GROWPHOT` extension, `phcal` calibrates that by default instead
(`--photometry aper` selects the apertures). If a frame was processed by an older
`aphot`, run `aphot` again before calibrating with the current `phcal`. The
old HWHM keyword held σ, so `phcal` would select a different aperture.
```
