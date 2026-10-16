# Accuracy and validation

How accurate is Munipack, and how do we know? Each method was tested on
frames made by `artificial`, where the true positions and fluxes of all stars
are known, and checked against independent calculations (astropy, numpy and
analytic formulas). Two kinds of frames were used: ground-based frames with a
bright sky, and faint low-count frames typical of a space ultraviolet
telescope such as QUVIK.

The short version: the results agree with the truth within their reported
errors, and the reported errors match the actual scatter. The table below
gives the key result for each area.

| Area | What to expect |
|---|---|
| Whole chain (`find` → `aphot` → `astrometry` → `phcal`) | Calibrated magnitudes agree with the catalogue to a few millimagnitudes on simulated frames, saturated 16-bit frames included. |
| Detection (`find`) | The threshold behaves as a false-alarm rate, so faint low-count frames do not fill up with noise detections. |
| Aperture photometry (`aphot`) | The sky is right even at a few electrons per pixel, and the flux errors agree with the actual scatter. |
| Total fluxes (`aphot --psf-file`, `--apcorr auto`) | Aperture fluxes corrected to the total light, also when the PSF changes over the field (with `--apcorr auto`, add `--apcorr-radial` when it changes with the distance from the centre). |
| PSF photometry (`psfphot`) | Unbiased fluxes and positions for well-measured stars; close pairs and crowded fields are handled far better than with apertures. |
| Growth curves (`gphot`) | Total fluxes within a fraction of a percent for bright stars, consistent with `psfphot`. |
| Calibration (`phcal`, `phfotran`) | Realistic zero-point errors even from a handful of stars; colour terms and instrument transformations are recovered accurately. |
| Astrometry | Typically better than 0.1 px, also with distortions; false matches are rejected. |
| Calibration frames (`kombine`, dark) | Unbiased master frames at low counts, with cosmic rays rejected. |
| Simulation (`artificial`) | Frames made with `--seed` are reproducible, whatever the number of processes. |

## Where the methods need care

A few methods depend on what is in the frame. Keep these in mind when you
work with sparse or crowded fields.

- **Automatic aperture correction** (`aphot --apcorr auto`) and **growth
  curves** (`gphot`) need enough bright, isolated stars. On a nearly empty
  frame give `aphot` a PSF file (`--psf-file`) instead.
- Both rely on the FWHM written by `find`. If the FWHM is badly off, for
  example because stars were split into several detections, isolation and
  the radial model go wrong.
- **PSF photometry** with a polynomial variation (`--psf-var`) cannot follow a
  PSF whose shape changes from the centre to the edge of the field. Use
  `--psf-var-radial` or a model PSF as the base in that case.
- The **empirical PSF** is built from the PSF stars, so these should be mostly
  single stars.

```{tip}
If the calibration stars scatter more than their photometric errors (for
example because the PSF varies over the field), `phcal` adds a systematic
error estimated from the residuals, so that the zero-point error stays
realistic. This is the default (`--syserr auto`); `--syserr 0` turns it off.
```

## Details

The individual tests, the frames they use and how to run them are described
in {doc}`tests`. The measured numbers for each fix, and the open issues, are
in [BUGS.md](https://github.com/heamuni/munipack/blob/develop/BUGS.md).
