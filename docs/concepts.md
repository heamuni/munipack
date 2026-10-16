# Concepts

The conventions shared by all Munipack actions: where results go, units,
errors, photometric systems and the PSF file format.

## Files and extensions

Munipack keeps everything about a frame in the frame's own FITS file. The
image itself is left untouched, and each processing step adds its results as
a new table extension. Refer to extensions by name (`frame.fits[APERPHOT]`),
not by number, because their order can change.

| Step | Adds | What is in it |
|---|---|---|
| `artificial` | new frame, `PSF`, `CATALOGUE` | the simulated image, the PSF used, and the true positions and fluxes of the stars |
| `find` | `FIND` | detected stars with positions and shape parameters, brightest first |
| `aphot` | `APERPHOT` | sky and counts with errors for every aperture |
| `gphot` | `GROWPHOT`, `GROWCURVE` | total counts from the growth curve, and the curve itself |
| `psfphot` | `PSFPHOT` | positions, counts, sky and fit quality from the PSF fit |
| `astrometry` | primary header | the WCS solution |
| `phcal` | new file `*_cal.fits` with `PHOTOMETRY` | the frame in photons and a table of calibrated fluxes and magnitudes |
| `phcorr`, `kombine` | image, `STDERR` | the corrected or combined frame and the error of each pixel |

`find`, `aphot`, `gphot` and `psfphot` update the input file in place; use
`-t dir` to write the updated copy elsewhere. `phcal` writes a new file
(by default `frame.fits` gives `frame_cal.fits`).

A value of −1 in a count or error column means the star could not be measured
(it is at the frame edge, saturated, fainter than the sky, or blended with a
brighter neighbour). Undefined magnitudes are written as 99.99999.

### Header keywords

These keywords of the primary header describe the detector. Set them
correctly, because they drive the errors and the detection threshold.

| Keyword | Units | Meaning |
|---|---|---|
| `GAIN` | e⁻/ADU | electrons per unit of the pixel values (1 for a frame in electrons) |
| `READNS` | ADU | read noise in the units of the frame |
| `RDNOISE` | e⁻ | read noise in electrons (used when `READNS` is missing) |
| `BIASLEV` | ADU | bias level still present in the values |
| `SATURATE` | ADU | saturation level |
| `EXPTIME` | s | exposure time |
| `AREA` | m² | collecting area of the telescope (1 when missing) |
| `FILTER` | — | filter name, matched to the photometric system |
| `PHOTSYS` | — | instrumental photometric system of the frame |

## Units

### ADU and electrons

Counts (`APCOUNT`, `GCOUNT`, `PSFCOUNT`, `SKY`) are always in the units of
the frame: ADU for a raw frame, electrons after `phcorr` has converted it
(it does so when the gain is known). No other action rescales the data; the
gain is used only to estimate the noise.

- **Raw frame:** `GAIN` is the detector gain, the read noise is in `READNS` or
  `RDNOISE`, and `BIASLEV` is set if the bias has not been subtracted.
- **Frame in electrons:** `GAIN = 1`, the read noise is in electrons and there
  is no `BIASLEV`. `phcorr` writes exactly this header for you.

```{note}
A wrong `GAIN` does not change the fluxes, but it gives wrong errors and a
wrong detection threshold. When `GAIN` is missing, `aphot` assumes
1 e⁻/ADU and warns.
```

After calibration, `CTPH` in the `PHOTOMETRY` header converts counts to
photons (photons = `CTPH` × counts). Despite its header comment ("counts
per photons"), it is photons per count. It absorbs throughput, quantum
efficiency, gain, light lost outside the aperture and an unset `AREA`.

### Angular scale

`--scale` is in **pixels per degree** in every action that takes it
(`artificial`, `astrometry`, `kombine`, `aphot`, `psfphot`). Compute it as
3600 / (arcsec per pixel). Older Munipack versions used degrees per pixel in
`artificial`, `astrometry` and `kombine`, so update old scripts.

| Plate scale | `--scale` |
|---|---|
| 2″/pix | 1800 |
| 1″/pix | 3600 |
| 0.5″/pix | 7200 |

Aperture radii, FWHM and fitting radii are in pixels. In the `PHOTOMETRY`
extension of a calibrated frame, the aperture radii are in degrees (when the
frame has an astrometric solution) and the sky is in photons per pixel.

## How errors are estimated

Every pixel's noise is the Poisson noise of its electrons plus the read
noise, both taken from the header keywords above. Sky levels are clipped
means rather than medians, because a median is biased low on a faint sky.

- **Aperture photometry:** the error combines the sky noise inside the
  aperture, the Poisson noise of the star, and the uncertainty of the sky
  level. The last part grows with the aperture size.
- **PSF photometry:** a weighted least-squares fit, with the noise of each
  pixel predicted from the model. For an empirical PSF, the uncertainty of the
  PSF itself is included.
- **Calibration:** `phcal` adds the error of the star's counts, of the
  calibration ratio and a systematic error. By default (`--syserr auto`) the
  systematic error is the extra scatter of the calibration stars, fitted in a
  way that ignores outliers such as variables or blends. With few calibration
  stars (`phcal` warns below 10), give a fixed `--syserr` measured on many
  frames; `--syserr 0` turns it off.
- **Pre-processing and stacking:** `phcorr` and `kombine` write the error of
  each output pixel to the `STDERR` extension.

The error of the calibration ratio is common to all stars on a frame, so it
cancels in differential photometry.

## Photometric systems

A photometric system is a FITS table with one row per filter. The `FILTER`
column is matched against the `FILTER` keyword of your frames.

| Column | Meaning | Units |
|---|---|---|
| `FILTER` | filter name | — |
| `LAM_EFF` | effective wavelength | m |
| `LAM_FWHM` | passband width | m |
| `NU_EFF`, `NU_FWHM` | the same in frequency (informative only) | Hz |
| `FNU_REF` | flux of a zero-magnitude star per Hz (informative only) | W/m²/Hz |
| `FLAM_REF` | flux of a zero-magnitude star per nm | W/m²/nm |

Only `LAM_EFF`, `LAM_FWHM` and `FLAM_REF` are used in the calibration.
The passband is treated as a Gaussian with σ = `LAM_FWHM`/2, so for your own
filters `LAM_FWHM` means 2σ rather than the true FWHM. The
extension name is the name of the system, and its header also needs `FILTREF`
(a reference filter) and `REFID` (a literature reference).

Systems bundled with Munipack (see `phcal --list`):

| System | Filters | Reference |
|---|---|---|
| `Johnson` | U, B, V, R, I | Bessell (1992) |
| `Gunn` | g, r, i, z | Schneider, Gunn & Hoessel (1983) |
| `Stromgren` | u, v, b, y, Hbetan, Hbetaw | Fiorucci & Munari (2003) |

To add your own system, copy `photosystems.lst` (installed in
`share/munipack/`; `/opt/munipack/share/munipack/` in the Docker image),
append a block `# BEGIN HDU n` … `# END HDU n` for your system with the same
columns, convert it to FITS and pass it with `--phsystab`. Your file replaces the
bundled one, so keep in it the systems you still need:

```bash
munipack fits --restore my_photosystems.lst      # writes my_photosystems.fits
munipack phcal --phsystab my_photosystems.fits --photsys-ref MySystem -f NUV \
    -c cat.fits --col-ra RAJ2000 --col-dec DEJ2000 \
    --col-mag NUVmag --col-magerr e_NUVmag frame.fits
```

The same file works with `artificial --phsystab` to simulate frames in that
band.

### Calibrated quantities

`phcal -q` chooses which quantities appear in the output table (by default
`FLUX` and `MAG`):

| Quantity | Meaning | Units |
|---|---|---|
| `PHOTON` | number of photons | photons |
| `PHOTRATE` | photon rate per area | ph/s/m² |
| `PHOTNU` | photon rate per Hz | ph/s/m²/Hz |
| `PHOTLAM` | photon rate per nm | ph/s/m²/nm |
| `FLUX` | energy flux in the band | eV/s/m² |
| `FNU` | flux density per Hz | W/m²/Hz |
| `FLAM` | flux density per nm | W/m²/nm |
| `MAG` | magnitude in the photometric system | mag |
| `ABMAG` | AB magnitude | mag |
| `STMAG` | ST magnitude | mag |

(psf-format)=
## PSF file format

One file format describes a PSF, either constant or varying over the field.
It is read by `artificial`, `aphot` and `psfphot` (`--psf-file`) and written
by `psfphot --psf-save`, so a frame simulated with a PSF file is fitted with
exactly the same PSF.

| Where | Keyword | Meaning |
|---|---|---|
| primary header | `ZOOM` | PSF samples per pixel along each axis |
| each image extension | `RADCOO` | distance of this PSF from the field centre [deg] |
| | `ANGCOO` | position angle of this PSF, counter-clockwise from the *x* axis [deg] |
| | `CRPIX1`, `CRPIX2` | optional: the sample at the star's centre (default: the middle of the stamp) |

Each image extension with `RADCOO` and `ANGCOO` holds one PSF; those with
the same `RADCOO` form a ring, and Munipack interpolates between them.
The field centre is the centre of the frame; `--psf-centre x,y` of `aphot`
and `psfphot` moves it. A file
with just a primary image and `ZOOM` is one PSF for the whole field. The
samples should sum to ZOOM² (4 for `ZOOM = 2`) so that a star gets all of its
flux.

A minimal file with one Gaussian PSF (σ = 1.5 pixels), written with Astropy:

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

To make the PSF vary over the field, add more extensions with the same
`ZOOM`. For example, four PSFs at `RADCOO = 0.3` and `ANGCOO = 0, 90, 180,
270` form a ring 0.3° from the centre.

```{note}
When the frame has no WCS, `psfphot --psf-save` writes `RADCOO` in pixels
instead of degrees; use such a file with `--scale 1`.
```
