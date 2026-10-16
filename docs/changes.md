# Differences from upstream

This fork is based on upstream Munipack Mercurial revision 1653 (`ae465635f925`, reported as `0.6.2+hg1653`). The untouched upstream tree is commit `003cb0d`, and each fix sits on top of it as its own commit, so `git diff 003cb0d -- src` shows every change to the code.

This page sums up the changes. The full catalogue, with the measured effect of each bug, how it was checked and the commit that fixes it, is in
[BUGS.md](https://github.com/heamuni/munipack/blob/develop/BUGS.md). The identifiers in brackets below (A…, P…, S…) are its entries. The accuracy measured after the fixes is collected in {doc}`validation`, and the open issues and suggested improvements are in [BUGS.md §6](https://github.com/heamuni/munipack/blob/develop/BUGS.md#6-suggested-further-improvements-not-done).

```{note}
Many of the photometric fixes matter most on **low-count frames**, where the sky is a few photo-electrons per pixel (for example space-UV imaging). On ground-based frames with a sky of hundreds of electrons, most photometric results change by less than their noise.
```

## Bug fixes

### Frame simulation (`artificial`)

- **Detector noise model** (§2.1). The old `--gnoise` raised the sky level and added no noise of the given σ. There was no read noise, and the sky noise was a clipped Gaussian. Frames now follow this chain: mean in electrons → Poisson → read noise → ÷ gain → + bias → saturation.
- **Unusable frames after the `--scale` unit change** (A1). In hg1653 the scale became pixels per degree, but the sky, galaxy and PSF-file formulas still used degrees per pixel.
- **Reproducibility** (A2). The new `--seed` gives the same frames whatever the `--mproc`.
- **Geometry and PSF** (A5, A6, A27, A28). Fixes for:
  - FFT PSF off-centre (a −0.095 px astrometric bias);
  - pixel integration counting edge subpixels twice;
  - `--psf-file` stars offset by −0.5 px;
  - stars rounded to the PSF sample grid (up to ±0.25 px).
- **Sky and time** (A3, A9, A11, A33):
  - the sky gradient had the wrong sign;
  - the sky darkened with airmass;
  - the light curve and airmass were taken at the start of the exposure, not mid-exposure;
  - `--sky-mag 99` still added about 1 e⁻/s/pixel.
- **Output** (A7, A8, A22, A34). Fixes for:
  - broken `BZERO` for `-B 8` and `-B 32`;
  - `SATURATE`/`GAIN` written with 3 significant digits;
  - crashes and truncated files caused by the watermark;
  - integer frames truncated instead of rounded.
- **Random star fields** (A23). The density per square degree was applied to the whole frame whatever its size.

### Detection and photometry (`find`, `aphot`, `psfphot`, `gphot`)

- **Biased sky at low counts** (P23). The sky of `aphot`, `find` and `psfphot` was a robust mean or a median. On a Poisson sky of a few electrons per pixel these lie below the true mean, so faint stars came out up to 23% too bright. The sky is now a mean clipped by the frame's noise model, and `aphot` masks the neighbouring stars in its sky ring.
- **`find` threshold and shape tests** (P13, P25, P32, P33, P37):
  - saturated pixels counted as valid;
  - noise peaks in the wings of bright stars detected as stars;
  - faint stars lost to the sharpness test;
  - far more noise detections than the threshold implied on skewed (Poisson) noise;
  - read noise counted twice.
- **Aperture errors** (P34, P36). `APCOUNTERR` ignored the frame's gain. The enclosed-energy errors treated nested apertures as independent.
- **Neighbour contamination** (P29, P46). Apertures that hold light from a brighter neighbour are now invalidated. A first version of this check also hit stars that had no real contamination, which is fixed too.
- **Header values** (P7, P8, P9). The growth-curve profile, the `HWHM` keyword (it stored σ) and the choice of aperture were wrong.
- **`gphot`** (P5, P6, P24, P38, P49):
  - errors too small;
  - sky correction with the wrong sign;
  - no valid fluxes on frames with extended PSFs.

  `gphot` is rebuilt on the enclosed energy of `aphot`; `--legacy` keeps the old method.
- **`psfphot`** (P17–P22, P35, P42, P50–P53). This is a new tool (see below), and its review fixed:
  - saturated neighbours that were not subtracted;
  - Moffat fits on Gaussian stars that did not converge;
  - close pairs;
  - the variance floor at low sky;
  - a high fixed sky;
  - flux lost beyond the empirical PSF table.

### Calibration (`phcal`, `phfotran`, `phcorr`, `dark`, `kombine`)

- **Physical quantities** (P1, §1.2). These `phcal -q` outputs were wrong: FLAM by ×0.0013, STMAG by +9.7 mag, PHOTLAM by ×7000. FNU, PHOTNU, ABMAG and `PHOTFLAM` were wrong too. MAG was always correct.
- **Several frames and filters** (P3, P14, P43). Fixes for:
  - reference photons scaled by the product of all exposures, followed by a crash;
  - frames paired with filters by position instead of by `FILTER`;
  - a spurious colour correction of the outputs (0.1–0.2 mag per mag of B − V).
- **Manual calibration and options** (P2, P4). `-C` left the transformation uninitialised. `-q PHOTRATE` was rejected.
- **`phfotran`** (P45). It was listed but never dispatched, and it crashed when run. `phcal --tratab` also used an inconsistent backward matrix.
- **Systematic error** (P44, P54, P55). The estimate of `--syserr auto` was biased by outliers and by a small number of stars, and its zero-point errors did not include the uncertainty of σ.
- **Calibration frames** (P11, P31, P39, P41):
  - `phcorr` ignored `-box` and propagated errors wrongly;
  - its output header kept the raw gain and the read noise in ADU;
  - master darks of a few electrons came out 15–23% low;
  - a single dark frame with `-bias` was not bias-subtracted.
- **`kombine`** (S9–S14, S20, P30):
  - the output grid was sampled at pixel corners (+0.5 px);
  - the bicubic stencil was wrong;
  - edge pixels were never resampled;
  - `--pm-ra`/`--pm-dec` had no effect;
  - the output was shifted when the first frame's CRPIX was off-centre;
  - `STDERR` was on the wrong scale;
  - the default mean was biased at low counts.

### Astrometry (`astrometry`)

- **Mirrored and rotated frames** (S1–S3). Mirrored frames never solved. `pix2uv` was not the inverse of `uv2pix`. CDELT/PC headers were read wrongly.
- **Distortion** (S4, S19). The TPV coefficients were written in the wrong frame. Large distortions are now fitted in stages: errors of 0.03–0.11 px where there were several pixels.
- **False solutions** (S16, S23). A match is now tested for significance against chance coincidences, counting the trials of the search. Before, saturated or dense frames could get a wrong WCS.
- **Robustness** (S5, S6, S8, S24). Fixes for:
  - `--mode=cross` never matching;
  - crashes after a failed RMS check;
  - the last star of each chunk dropped;
  - crashes on frames with fewer than five stars, which now fail on their own while the rest of the run continues.
- **Accuracy** (S7, S25, S26, S27):
  - RA residuals without cos δ;
  - pmRA not divided by cos δ;
  - catalogue and frame not sorted by brightness;
  - an unweighted final fit (rms error 0.099 → 0.055 px).

### Command line

- **Exit status** (A24): 0 on success and 1 on every error, with the reason on stderr. Before, errors gave 139 (it looked like a crash), 13 or 255, and `--help` gave 255.
- **Options** (A14–A17, P26, S17):
  - an unknown `--psf` or `--spread` name now gives an error (it silently fell back to the default);
  - a paired option given without its partner gives an error (it was silently dropped);
  - floating-point values reach the utilities at full precision;
  - `--mproc` below 1 is rejected;
  - `-m cross` is no longer parsed as the `cross` action.
- **`fits --update` or `--remove-keys` on an extension deleted the file** (S28), for example `'frame.fits[FIND]'`. It is now updated in place.

## New features

| Feature | Tool and options | Entries |
|---|---|---|
| Detector noise model | `artificial --dark --readnoise --bias --seed`, `--no-watermark` | §2.1, A31 |
| PSF photometry, with groups, an empirical PSF and detection in the residuals | `psfphot` (replaces the unused `pphot` stub) | P15–P22, P50–P53 |
| Enclosed-energy aperture correction | `aphot --psf-file`, `--psf`, `--apcorr auto`, `--apcorr-radial`; `phcal --apcorr` | P36, P40 |
| Systematic error of the calibration | `phcal --syserr auto` (default) | P44, P54, P55 |
| Colour terms, up to three colours | `phcal --color-term`, `--color-order`, `--color-ref` | P47, P56 |
| Magnitude transformation model | `phfotran --model mag` (default) | P48 |
| False-alarm-rate threshold | `find -th`, `--false-rate`, `--shape-sigma` | P32, P33 |
| Matching by the header WCS or by a pointing, roll known or unknown | `astrometry --wcs-offset`, `--wcs-rot`, `--wcs-scale`, `--roll-unknown`, `--reflex-unknown`, `--blind` | S21, S29, S31 |
| Scale prior for the blind match | `astrometry --scale-tol`, `--scale-free` | S22 |
| Parallel processing | `--mproc` (coarrays) for `artificial` and `astrometry`; OpenMP threads for the blind match | A25, S30 |
| Flux-conserving resampling | `kombine --flux-conserve` | S14 |
| Clipped mean for stacks and masters | `kombine`, `dark` (default); `--hampel` for the old mean | P30, P31 |

## Compatibility notes

Most of these behaviour changes affect results computed by older versions. The full list is in
[BUGS.md §5](https://github.com/heamuni/munipack/blob/develop/BUGS.md#5-compatibility-notes-behaviour-that-changed).

```{warning}
**Re-run `aphot` on frames processed by an older version** before using the new `phcal`. The `HWHM` keyword of APERPHOT is now a true HWHM (it was σ), and `APCOUNTERR` now uses the gain. Old tables make `phcal` pick a smaller aperture and weight the stars wrongly.
```

- **`--scale` is in pixels per degree** in `artificial`, `astrometry` and `kombine`. This is an upstream change in hg1652/1653; older versions used degrees per pixel.
- **`artificial --hwhm`** is σ for the GAUSS and SEEING PSFs and α for MOFFAT, not a half width. This upstream behaviour is kept and now documented.
- **`artificial` frames differ** from older versions:
  - Poisson sky, with galaxies noised;
  - the flat applied before the noise;
  - a correct sky gradient, and a sky that brightens with airmass;
  - exact sub-pixel star positions;
  - random fields with the documented density, which follows the frame area.

  Noiseless frames made with `--gnoise` and no gradient, airmass or flat are bit-identical to hg1653.
- **Sky and errors.** SKY of `aphot`, `find` and `psfphot` rises by up to 0.05–0.07 e⁻/px on faint skies, and faint fluxes drop. `APCOUNTERR` of raw frames changes by up to √(1/GAIN). The gain applied is recorded as `APERRGN`.
- **`find` detection counts change.** The threshold is a false-alarm rate, so on faint skies there are 2–5× fewer noise detections, and the shape tests accept more faint stars. Re-tune `-th` if you depend on a given depth.
- **`phcal`:**
  - adds a systematic error by default (`--syserr auto`; `--syserr 0` gives the old errors bit for bit);
  - divides aperture counts by `APEE` when present (`--apcorr none` ignores it);
  - with several filters, calibrates each filter as if it were alone.
- **`phfotran`** fits the magnitude model by default. `--model flux` gives the old matrix tables. Tables without `TRAMODEL` are read as flux tables.
- **`astrometry`:**
  - needs at least five stars;
  - rejects matches that are not significant (`--disable-rms-check` turns the test off);
  - matches by the header WCS, or by `--rcen`/`--dcen` with a scale (`--scale` or the WCS), before the blind match (`--blind` restores the old blind match, `--roll-known` the old order);
  - takes the scale from the header as a ±5% prior (`--scale-free` turns it off);
  - reads the whole catalogue and, with a known pointing, uses the brightest `--maxstars` stars near the frame, so cones larger than the frame give other stars than before.
- **`kombine` and `dark`** default to the clipped mean. `--hampel` restores the old robust mean, and `--arith` gives the plain mean. `kombine` `STDERR` is now on the scale of the image.
- **`phcorr` output** in electrons has `GAIN = 1` (the conversion gain is in `PHCGAIN`), `READNS`/`RDNOISE` in electrons and no `BIASLEV`.
- **`psfphot`** adds FLAG 256 (S/N < 2) and a `PSFMODERR` column, which is included in `PSFCOUNTERR`.
- **Exit status** is 0 or 1, as described above. A crashed utility still gives 256 minus the signal number.
