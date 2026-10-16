# Concepts

This page collects the conventions that every action shares: where the
results are stored, the units, how errors are propagated, how photometric
systems are described, and the PSF file format.

## Files and extensions

Munipack keeps the results of a frame in the frame's own FITS file. The
primary array stays as it is, and every step adds (or replaces) a binary-table
extension. Always address extensions by name (`frame.fits[APERPHOT]`), not by
index, because their order can vary.

| Step | Writes | Contents |
|---|---|---|
| `artificial` | new frame, `PSF`, `CATALOGUE` | the simulated image; the PSF used (with `ZOOM`; not written with `--psf-file`); the true `X`, `Y` and `COUNT` (mean electrons) of every star |
| `find` | `FIND` | detected stars: `X`, `Y`, `PEAKRATIO`, `SHARP`, `ROUND`, sorted by decreasing `PEAKRATIO` (the brightest first) |
| `aphot` | `APERPHOT` | per star: `X`, `Y`, `SKY`, `SKYERR`, `APCOUNTi`, `APCOUNTERRi` for each aperture *i*; optionally `APEEi`, `APEEERRi` (enclosed energy) |
| `gphot` | `GROWPHOT`, `GROWCURVE` | `X`, `Y`, `SKY`, `SKYERR`, `GCOUNT`, `GCOUNTERR`, `GROWFLAG`; the growth curve itself (`R`, `GROWCURVE`, `GROWCURVEERR`, `RADIALPROFILE`) |
| `psfphot` | `PSFPHOT` | `X`, `Y`, `XERR`, `YERR`, `PSFCOUNT`, `PSFCOUNTERR`, `SKY`, `SKYERR`, `CHI2`, `SHARP`, `ROUND`, `NPIX`, `FLAG`, `GROUP`, `GSIZE`; `PSFMODERR` for an empirical PSF |
| `astrometry` | primary header | the WCS (`CTYPEi`, `CRVALi`, `CRPIXi`, `CDi_j`, …), `RADESYS`, `JEPOCH`, `ASTRMETH`; `ASTROREPORT` with `--report` |
| `phcal` | new file `*_cal.fits` with `PHOTOMETRY` | the frame converted to photons, and a table `RAJ2000`, `DEJ2000`, `SKY`, `SKYERR`, `PHOTON`, `PHOTONERR`, then `Q`, `QERR` for every quantity of `-q` (default `FLUX`, `MAG`); `COLOR`, `COLORERR`, `COLORSRC` with a colour term; `PHRES` with `--advanced` |
| `phcorr`, `kombine` | primary array, `STDERR` | the corrected or combined frame and the standard error of every pixel |

`find`, `aphot`, `gphot` and `psfphot` update the input file in place;
`-t dir` writes the updated copy to `dir` instead (`-O` with `--pattern`,
`--mask`, `--format` sets other output names). `phcal` writes a new file, `frame.fits` →
`frame_cal.fits` by default.

A negative value (−1) in a count or error column means the value could not be
measured: the aperture crosses the frame edge, contains saturated or bad
pixels, the star is fainter than the sky, or the aperture is contaminated by a
brighter neighbour. Undefined magnitudes are 99.99999, with the error 9.99999.

### Header keywords

Keywords of the primary header that the actions read:

| Keyword | Units | Meaning | Used by |
|---|---|---|---|
| `GAIN` | e⁻/ADU | electrons per unit of the pixel values (1 for frames in electrons, also implied by `BUNIT = 'COUNTS'`) | find, aphot, psfphot, dark, kombine, phcorr |
| `READNS` | ADU | read noise in the units of the pixel values | find, aphot, psfphot, dark, kombine, phcorr |
| `RDNOISE` | e⁻ | read noise in electrons (used when `READNS` is missing) | as `READNS` |
| `BIASLEV` | ADU | level of zero electrons (bias pedestal) still in the values | find, aphot, psfphot, dark, kombine, phcorr |
| `SATURATE` | ADU | saturation level | find, aphot, psfphot, dark, flat, phcorr |
| `EXPTIME` | s | exposure time | phcal, phfotran, dark, phcorr, kombine |
| `AREA` | m² | collecting area (1 when missing) | phcal, phfotran |
| `FILTER` | — | filter name, matched to a row of the photometric system | phcal, phfotran |
| `PHOTSYS` | — | photometric system of the frame (instrumental) | phcal, phfotran |

`artificial` writes `GAIN`, `RDNOISE` (not `READNS`), `BIASLEV`, `SATURATE`,
`EXPTIME`, `AREA` and `FILTER`, plus `DARKCUR` [e⁻/s/pix], `ARTSEED`,
`WMARKX1/X2/Y1/Y2` (the region of the watermark) and `HWHM`, which there is
the PSF width parameter (σ for GAUSS and SEEING, α for MOFFAT), not a true
half width.

Keywords of the extension headers (selection):

| Extension | Keyword | Meaning |
|---|---|---|
| `FIND` | `FWHM` | typical FWHM of the stars [pix] |
| | `THRESH` | requested threshold [σ of Normal noise, as a false-alarm rate] |
| | `HTHRESH`, `SKEWFILT`, `KURTFILT` | threshold height under the noise model, skewness and kurtosis of the filtered sky |
| | `FALSERAT` | expected number of noise detections per frame |
| | `LOWBAD`, `HIGHBAD` | range of good data [ADU] |
| | `SHAPESIG`, `SHARPREF`, `ROUNDREF`, `SHRPLO`, `SHRPHI`, `RNDLO`, `RNDHI` | shape acceptance |
| `APERPHOT` | `NAPER`, `APERi`, `ANNULUS1`, `ANNULUS2` | aperture and sky-ring radii [pix] |
| | `HWHM` | true half width at half maximum [pix] (1.1774 σ of a Gaussian) |
| | `APERRGN` | gain [e⁻/ADU] applied to the Poisson term of `APCOUNTERR` |
| | `EESOURCE`, `EEPSF`, `EEVAR`, `EERTOT1/2`, `EER99`, … | source and model of the enclosed energy (`APEEi`) |
| `GROWPHOT` | `HWHM`, `RADFLX90`, `GROWMETH` | half width, radius of 90 % of the flux [pix], `EE` or `LEGACY` |
| `PSFPHOT` | `PSFTYPE`, `FWHM`, `PSFSIGMA` / `PSFALPHA`, `PSFBETA` | PSF model and its shape |
| | `GAIN`, `READNS`, `BIASLEV`, `SATURATE` | the noise model used for the weights |
| | `PSFNSTAR`, `PSFRAD`, `PSFZOOM`, `PSFVAR`, `PSFVARR`, `PSFSYSER` | empirical PSF: stars, table radius, sampling, variation, relative systematic error |
| `PHOTOMETRY` | `PHOTOTYP` | photometry calibrated: `APHOT`, `GPHOT` or `PSFPHOT` |
| | `CTPH`, `CTPHERR` | the calibration ratio (photons per count, see below) and its error |
| | `PHSYSERR`, `PHSYSERE`, `PHSYSNST` | systematic error [mag], its upper error, calibration stars used |
| | `APCORTYP`, `APCORR` | aperture correction: `EE`, `CONST` (the factor in `APCORR`) or `NONE` |
| | `APER`, `ANNULUS1`, `ANNULUS2` | calibration aperture and sky ring [deg] |
| | `COLTERM`, `COLTREF`, `COLTK1`, `COLTK1E`, … | colour term, if any |

The primary header of a calibrated frame also gets `BUNIT`, `PHOTFLAM`,
`PHOTZPT`, `PHOTPLAM` and `PHOTBW` (HST-style keywords), and `PHOTSYS`.

## Units

### Counts, ADU and electrons

*Counts* in Munipack are values in the units of the frame. `APCOUNT`, `GCOUNT`,
`PSFCOUNT` and `SKY` are therefore in ADU for a raw frame, and in electrons
for a frame that `phcorr` has converted (it converts when the gain is known,
from `GAIN` or `-gain`). No other action multiplies by the gain: the gain only
enters the **noise model**, through `GAIN`.

- A raw frame: `GAIN` is the detector gain, `READNS` [ADU] or `RDNOISE` [e⁻],
  and `BIASLEV` if the pedestal has not been subtracted.
- A frame in electrons: `GAIN = 1`, the read noise in electrons, and no
  `BIASLEV`. `phcorr` writes exactly this header, with `BUNIT = 'COUNTS'`, and
  keeps the conversion gain in `PHCGAIN`.

```{note}
A wrong `GAIN` gives wrong errors and a wrong detection threshold, but the
fluxes are not scaled by it. If the header has no `GAIN`, the sky clipping
falls back to the observed scatter of the sky, and `aphot` takes 1 e⁻/ADU for
the Poisson term of the errors (with a warning).
```

`CTPH` converts the counts to photons, *n* = `CTPH` · *c*. Despite its name
and its header comment ("counts per photons"), it is photons per count, the
inverse of the efficiency. It absorbs everything between the photons arriving
at the telescope and the values in the frame: throughput, quantum efficiency,
the gain (for ADU frames), the light outside an uncorrected aperture, and an
`AREA` that is not set. For a frame in electrons with the right `AREA` it is
therefore at least 1; for an ADU frame it is about `GAIN` / efficiency.

### Angular scale

`--scale` is in **pixels per degree** in `artificial`, `astrometry`,
`kombine`, `aphot` and `psfphot`. Older Munipack versions used degrees per
pixel in `artificial`, `astrometry` and `kombine`.

| Plate scale | `--scale` [pix/deg] |
|---|---|
| 2″/pix | 1800 |
| 1″/pix | 3600 |
| 0.756″/pix | 4761.9 |
| 0.396″/pix | 9090.9 |

In general `--scale` = 3600 / (arcsec per pixel).

Other units:

| Quantity | Units |
|---|---|
| aperture and ring radii, FWHM, `--fitrad`, `--psf-radius` | pixels |
| `APER`, `ANNULUS1/2` in `PHOTOMETRY` | degrees (with an astrometric solution) |
| `SKY`, `SKYERR` in `PHOTOMETRY` | photons per pixel |
| `artificial --hwhm` | pixels: Gaussian σ (GAUSS, SEEING) or Moffat α |
| `artificial --sky-mag` | mag/arcsec² |
| `astrometry --wcs-offset` | pixels, or with units: `30arcsec`, `2arcmin`, `0.1deg`, `20%` (of the frame) |
| `phcal --tol` | degrees |
| `phcal --syserr` | magnitudes |

## Noise model and errors

### Pixels

The variance of a pixel of value *t* [ADU] is taken as

```{math}
\sigma^2 = \frac{\max(t - \mathrm{BIASLEV},\,0)}{\mathrm{GAIN}} + R^2 ,
```

where *R* is the read noise in ADU (`READNS`, or `RDNOISE`/`GAIN`). The first
term is the Poisson noise of the electrons, the second is Gaussian. At a few
electrons per pixel the Poisson distribution is skewed, and Munipack accounts
for this where it matters:

- **Sky and averages** (the sky of `aphot`, `find`, `psfphot`; the master of
  `dark` and the default mean of `kombine`): the *mean* of the pixels, clipped
  by limits that follow the noise model (asymmetric, so that the Poisson tail
  is kept). It is not a median or a robust mean, because those are biased low
  on a faint sky.
- **Detection** (`find -th t`): the threshold height is chosen so that the
  skewed noise of the filtered sky exceeds it as rarely as Normal noise exceeds
  *t*σ. `--false-rate N` sets it by the expected number of noise detections per
  frame instead.

### Aperture photometry

For aperture *k* with area *A* = Σ*f* (where *f* is the fraction of each pixel
inside the aperture), sky-subtracted counts *C*, sky scatter σ_sky and error of
the mean sky σ_S (`SKYERR`):

```{math}
\sigma_C^2 = \sigma_\mathrm{sky}^2 \sum f^2 \;+\; \frac{C}{\mathrm{APERRGN}} \;+\; A^2\,\sigma_S^2 .
```

The terms are, in order: the noise of the sky pixels in the aperture (read
noise included), the Poisson noise of the star, and the sky-subtraction term.
The last term is the same for all stars that share a sky level, and it grows
with the area of the aperture. `APERRGN` is the frame's `GAIN` (1 when it is
unknown).

### PSF photometry

`psfphot` fits by weighted least squares. The variance of each pixel comes
from the model *M* (star + sky), not from the data:
max(*M* − `BIASLEV`, 0)/`GAIN` + *R*². There is no one-electron floor, so
χ² ≈ 1 on a sparse sky. For an empirical PSF, `PSFMODERR` (the uncertainty of
the PSF model) is included in `PSFCOUNTERR`.

### Calibration

`phcal` converts the counts *c* to photons, *n* = `CTPH` · *c*. Here *c* is
`APCOUNT`/`APEE` (or `APCOUNT` without EE), `GCOUNT` or `PSFCOUNT`, as chosen
by `--photometry`. The relative errors add in quadrature:

```{math}
\left(\frac{\sigma_n}{n}\right)^2 =
\left(\frac{\sigma_c}{c}\right)^2 +
\left(\frac{\sigma_\mathrm{EE}}{\mathrm{EE}}\right)^2 +
\left(\frac{\sigma_\mathrm{CTPH}}{\mathrm{CTPH}}\right)^2 +
s^2,
\qquad s = \frac{\mathrm{PHSYSERR}}{1.0857},
```

and `MAGERR` = 1.0857 σ_n/n. The colour terms, when present, add the errors of
their coefficients and of the star's colour.

- **Systematic error** `PHSYSERR` (`--syserr auto`, the default) is the scatter
  of the calibration stars beyond their own errors. It is fitted together with
  `CTPH` by a mixture likelihood, so outliers (blends, variables,
  misidentifications) carry no weight. It is ≈ 0 on well-behaved frames. With
  fewer than 10 calibration stars `phcal` warns; in that case pass a fixed
  `--syserr` value measured on many frames. `--syserr 0` turns it off.
- `CTPHERR` is common to all stars of a frame. It is a zero-point error,
  correlated between stars, and it cancels in differential photometry on one
  frame. `PHSYSERR` acts per star.

### Pre-processing and stacking

`phcorr` and `kombine` write the error of every output pixel to the `STDERR`
extension. For `phcorr` it includes the Poisson noise of the electrons (dark
current too), the read noise, the errors of the bias, dark and flat masters,
and the quantisation of integer frames. For `kombine` it is the standard error
of the combined value, on the scale of the output image.

## Photometric systems

A photometric system is a FITS table extension, one row per filter. The
`FILTER` column is matched against the `FILTER` keyword of the frames. The
values are in SI units, except `FLAM_REF` (per nm). The extension name is the
name of the system; its header also needs `FILTREF` (a reference filter of the
system) and `REFID` (a reference string).

| Column | Meaning | Units |
|---|---|---|
| `FILTER` | filter name | — |
| `LAM_EFF` | effective wavelength λ_eff | m |
| `LAM_FWHM` | passband width Δλ, taken as **2σ** of a Gaussian | m |
| `NU_EFF`, `NU_FWHM` | the same in frequency (informative only) | Hz |
| `FNU_REF` | f_ν of a zero-magnitude star (informative only) | W/m²/Hz |
| `FLAM_REF` | f_λ of a zero-magnitude star | W/m²/nm |

Only `LAM_EFF`, `LAM_FWHM` and `FLAM_REF` enter the calibration. The passband
is modelled as a Gaussian with σ = `LAM_FWHM`/2, so its equivalent (effective)
width is

```{math}
W = \sqrt{2\pi}\,\sigma = \sqrt{\pi/2}\;\mathrm{LAM\_FWHM} = 1.2533\,\mathrm{LAM\_FWHM}.
```

A Gaussian whose true FWHM was `LAM_FWHM` would give 1.0645 `LAM_FWHM`
instead. The expected photon flux of a star of magnitude *m* is

```{math}
\phi = \mathrm{FLAM\_REF}\; 10^{-0.4m}\; W\; \frac{\mathrm{LAM\_EFF}}{hc}
\quad [\mathrm{ph/s/m^2}], \qquad n = \phi\,A\,T .
```

Systems bundled with Munipack (`phcal --list`):

| System | Filters | Reference |
|---|---|---|
| `Johnson` | U, B, V, R, I | Bessell (1992) |
| `Gunn` | g, r, i, z | Schneider, Gunn & Hoessel (1983) |
| `Stromgren` | u, v, b, y, Hbetan, Hbetaw | Fiorucci & Munari (2003) |

For example, Johnson V has `LAM_EFF` = 553.8 nm and `LAM_FWHM` = 75.0 nm, so
*W* = 94.0 nm. To add a system, such as a UV band of your own instrument,
copy `photosystems.lst` (in the source tree `src/photometry/`, installed as
`share/munipack/photosystems.lst`, `/opt/munipack/share/munipack/` in the
Docker image). Append an extension (`# BEGIN HDU n` … `# END HDU n`) with the
same columns, your `EXTNAME`, `FILTREF`, `REFID` and one row per filter, then
convert the file to FITS (`file.lst` gives `file.fits`) and pass it with
`--phsystab`. The file replaces the bundled table, so keep the systems you
still need in it:

```bash
munipack fits --restore my_photosystems.lst      # writes my_photosystems.fits
munipack phcal --phsystab my_photosystems.fits --list
munipack phcal --phsystab my_photosystems.fits --photsys-ref MySystem -f NUV \
    -c cat.fits --col-ra RAJ2000 --col-dec DEJ2000 \
    --col-mag NUVmag --col-magerr e_NUVmag frame.fits
```

`artificial --phsystab my_photosystems.fits --photsys MySystem -f NUV`
simulates frames in the same band.

### Calibrated quantities

`phcal -q` selects the output columns. All of them are derived from the
photons with the same *W*. The frequency width is *W*_ν = *c W* / λ_eff².

| Quantity | Definition | Units |
|---|---|---|
| `PHOTON` | *n* | photons |
| `PHOTRATE` | ϕ = *n* / (*A T*) | ph/s/m² |
| `PHOTNU` | ϕ / *W*_ν | ph/s/m²/Hz |
| `PHOTLAM` | ϕ / *W* | ph/s/m²/nm |
| `FLUX` | *F* = ϕ *hc* / λ_eff | eV/s/m² |
| `FNU` | *F* / *W*_ν | W/m²/Hz |
| `FLAM` | *F* / *W* | W/m²/nm |
| `MAG` | −2.5 log₁₀(*F* / (`FLAM_REF` · *W*)) | mag |
| `ABMAG` | −2.5 log₁₀(f_ν / 3.631×10⁻²³ W/m²/Hz) | mag |
| `STMAG` | −2.5 log₁₀(f_λ / 3.6335×10⁻¹¹ W/m²/nm) | mag |

## PSF file format

One format describes a PSF, constant or varying over the field, for
`artificial --psf-file`, `psfphot --psf-file`, `aphot --psf-file` (the
enclosed energy), and the output of `psfphot --psf-save`. All of these use the
same code to interpolate and pixelise it, so a frame simulated with a PSF file
is fitted with exactly that PSF.

| Where | Keyword | Meaning |
|---|---|---|
| primary header | `ZOOM` | integer, PSF samples per pixel along each axis (the same for all PSFs) |
| each image extension | `RADCOO` | angular distance of this PSF from the field centre [deg] |
| | `ANGCOO` | position angle atan2(*y*, *x*) of the offset from the field centre [deg], counter-clockwise from the *x* axis |
| | `CRPIX1`, `CRPIX2` | optional: the sample (1-based) at the centre of the star; default (`NAXIS1`/2, `NAXIS2`/2), integer division |

- **Grid.** Every extension that has both `RADCOO` and `ANGCOO` is one PSF.
  `RADCOO = 0` is the centre. PSFs with equal `RADCOO` form a ring. The number
  of rings, the number of angles per ring, the sizes of the stamps and the
  order of the extensions are all free. Extensions without both keywords are
  ignored. A file with only a 2-D primary image and `ZOOM` is one PSF for the
  whole field.
- **Normalisation.** The samples should **sum to ZOOM²** (for example 4 for
  `ZOOM = 2`); a star then gets its whole flux. A sum *s* gives the fraction
  *s*/ZOOM², which can model light outside the stamp. A sum that differs from
  ZOOM² by more than 0.1 % triggers a warning.
- **Field position.** A star at pixel offsets (*x*, *y*) from the frame centre
  (pixel `NAXIS`/2) is at radius √(*x*² + *y*²) / *scale* [deg], with *scale*
  in pix/deg (the WCS, or `--scale`), and at angle atan2(*y*, *x*).
  `--psf-centre x,y` of `aphot` and `psfphot` moves the centre.
- **Interpolation.** On a ring, linear in angle (periodic). Between the centre
  and the first ring, and between rings, linear in radius. Beyond the last
  ring, the last ring is used.
- **Pixelisation.** The image of a star is the exact integral over each
  detector pixel of the interpolated samples: linear between samples for even
  `ZOOM`, constant per sample for odd `ZOOM`. Flux is conserved at any
  sub-pixel position. The analytic PSFs of `artificial` and `psfphot` are
  sampled with `ZOOM = 16` and integrated in the same way.

`psfphot --psf-save` adds `PSFCEN1`, `PSFCEN2` (the field centre [pix]),
`PSFSCALE` [pix/deg], `FWHM` and `PSFRAD` [pix] to the primary header. With
no WCS on the frame, `PSFSCALE = 1`: `RADCOO` is then in pixels, and the file
is used with `--scale 1`.

A minimal file with one Gaussian PSF (σ = 1.5 px), written with Astropy:

```python
import numpy as np
from astropy.io import fits

zoom, sigma = 4, 1.5                 # samples per pixel, sigma [pix]
k = 12 * zoom                        # half size: 12 pixels
y, x = (np.mgrid[-k:k + 1, -k:k + 1]) / zoom
psf = np.exp(-(x**2 + y**2) / (2 * sigma**2))
psf *= zoom**2 / psf.sum()           # normalise to ZOOM^2

pri = fits.PrimaryHDU()
pri.header['ZOOM'] = zoom
ext = fits.ImageHDU(psf.astype(np.float32))
ext.header.update(RADCOO=0.0, ANGCOO=0.0, CRPIX1=k + 1, CRPIX2=k + 1)
fits.HDUList([pri, ext]).writeto('psf.fits', overwrite=True)
```

To make it vary over the field, add more extensions with the same `ZOOM`.
For example, four PSFs at `RADCOO = 0.3` and `ANGCOO = 0, 90, 180, 270`
make one ring at 0.3° from the centre.
