# Accuracy and validation

The fixes and new methods were validated on frames simulated by `artificial`,
where the true fluxes and positions are known, and against independent
calculations (astropy, numpy, analytic formulas). This page collects the
measured accuracy, the performance and the known limits. The bug identifiers
(A…, P…, S…) refer to
[BUGS.md](https://github.com/heamuni/munipack/blob/develop/BUGS.md); the
regression tests that keep these properties are described in {doc}`tests`.

Two kinds of simulated frames are used throughout:

- **Ground-based**: Gaussian or Moffat PSFs and a sky of hundreds of
  electrons per pixel (for example gain 2, read noise 5 e⁻, sky 560 e⁻/px).
- **Low-count space UV**: u-band frames of 2048² or 4096² pixels, typically
  16 bit, gain about 0.1 e⁻/ADU, read noise about 0.4 e⁻, 20–60 s exposures
  and a sky of a few tenths of an electron to a few electrons per pixel,
  modelled on the QUVIK mission. Some use
  `dpsf.fits`, a field-dependent PSF: a 1.7 px core mixed linearly with the
  field radius with an 8.4 px halo.

A pull is (measured − true)/error; its standard deviation (pull SD) is 1 when
the errors are right. Q is the scatter of the calibration stars about the
zero point in units of their errors, also 1 when the errors are right.

## End to end

`artificial` → `find` → `aphot` → `astrometry` → `phcal` on unsaturated
frames gives the same results before and after the fixes: astrometry RMS
0.08″. Saturated 16-bit frames, where the unfixed build failed the
astrometry, now solve: 0.12″ RMS, 0.05 px from the truth, MAG − V = +0.005 on
60 stars.

## Detection (`find`)

- **Noise detections** at `-th 4` on a 2048² frame: about 130 (1 at `-th 5`),
  where a sky of 0.4 e⁻/px gave 505 (21) before the threshold became a
  false-alarm rate (P33).
- **Shape tests** accept 100% of the stars at u = 15–17 on low-count frames,
  instead of 76–84% (P32).

## Sky and aperture photometry (`aphot`)

- **Sky**: on a sky of a few photo-electrons per pixel the old robust mean or
  median was 13–30% low, and faint stars came out up to 23% too bright
  (P23). On ground-based frames SKY changes within its noise.
- **Errors**: APCOUNTERR of raw frames with GAIN 0.1 grows up to 3× (P34);
  pulls ≈ 1 on simulated frames.
- **Enclosed energy from a PSF** (`--psf-file`, `--psf`): on noiseless frames
  APCOUNT/APEE equals the true flux to 10⁻⁶ in every aperture, also for
  Gaussian and Moffat PSFs.
- **Automatic enclosed energy** (`--apcorr auto`): the normalisation is within
  0.016–0.052% rms of the exact EE over sets of 8 frames; the old 0.5% plateau
  criterion let the plateau wander to 20–40 px on dense fields, with totals
  0.1–0.4% too large (P40).
- **Field-dependent PSF** (4096² space-UV frames, 4 seeds, aperture 5.1 px,
  stars < 17 mag): with one zero point the residuals against the field radius
  go from −0.12…−0.07 mag at the centre to +0.07…+0.08 mag at the edge. With
  `aphot --psf-file` every radial bin is within ±0.015 mag and the robust
  scatter drops from 0.073–0.077 to 0.041–0.047 mag; `--apcorr auto
  --apcorr-radial` gives the same (±0.016 mag, 0.040–0.046 mag). On Gaussian
  frames the automatic EE is constant and the results are unchanged.

```{note}
**Limits.** The automatic EE needs bright isolated stars (about 30). With so
few, a 3σ fluctuation in the annulus of the brightest one can send the plateau
to the cap (one of 8 sparse Gaussian frames: −0.14%, within its APEEERR). It
also needs a correct FWHM in `FIND`: on a frame with a 9.9 px Gaussian PSF
found with FWHM 5.2 (split detections), `--apcorr-radial` fitted a spurious
variation (scatter 0.051 → 0.070 mag), while the x,y model stayed constant.
```

## PSF photometry (`psfphot`)

Analytic and tabulated PSFs (8 seeds, realistic noise):

- flux bias < 0.02% and position bias < 0.001 px for S/N > 50; reported
  errors 0.90–1.10 × the scatter, for Gaussian and Moffat (shape fitted) and
  for PSF files (`dpsf.fits`, a two-ring file);
- stars with a neighbour of equal brightness within 2 FWHM: bias +0.07%,
  scatter 2.5% (apertures: +20%, 37%);
- χ² ≈ 1 on a sky of 0.4 e⁻/px, without the old one-electron variance floor
  (P35);
- end-to-end `phcal`: MAG − V = 0.000 ± 0.0003, scatter 0.0017 mag (apertures
  0.0020);
- Moffat fits of Gaussian stars converge: 454 → 10 s per 1024² frame (P18).

Empirical PSF, groups and residual detection (8 seeds):

- constant PSFs (Gaussian σ 1.5 and 3, Moffat β 2.5, elliptical Moffat, 16-bit
  with saturated stars): the empirical PSF agrees with the true one within
  0.05% in flux for 10 to 100 PSF stars; pull SD 0.91–1.08, position bias
  < 0.002 px;
- a smooth field variation (focus and radial astigmatism): `--psf-var 2`
  +0.005% (true −0.005%), < 0.07% in every radius bin;
- correction mode on `dpsf.fits` frames: the correct base gives the fluxes of
  `--psf-file` (|Δ| < 0.001% on 16 crowded 16-bit seeds); a 10% broader or an
  asymmetric base −0.03% / +0.13% (uncorrected: +5.1% / +0.9%);
- `--psf-save` round trip (99% of the stars): ≤ 0.03% for a constant PSF
  (Gaussian ≤ 0.001%), ≤ 0.04% for order 2;
- pairs at 1 FWHM: groups give +0.0% bias and pull SD 1.0, iterative
  subtraction +26% / −19%;
- residual detection finds 88–100% of the companions at 0.5–1 FWHM missed by
  `find`, with no spurious detections (7 runs × 3200 stars: float, 16-bit with
  saturated stars, Moffat; true, mismatched and empirical PSFs);
- crowded field (2500 stars on 1024²): stars with a neighbour within 2 FWHM
  scatter 0.87% instead of 1.32%; those with an undetected companion within
  1 FWHM +0.24% instead of +17.5% with one residual pass;
- end-to-end `phcal` with `EMPIRICAL`: MAG − V within ±0.0004 mag, scatter
  0.0016–0.0019 mag.

Normalisation and radial variation (P50–P53; 2048² space-UV frames with
`dpsf.fits`, 8 seeds, grids of 1764 stars u = 10–16 and an SDSS cone; bias of
S/N > 50 stars in field-radius bins < 500 / 500–900 / > 900 px):

| Model | Bias | By radius | Pull SD |
|---|---|---|---|
| `--psf-var-radial 1` | −0.022% | +0.002 / −0.041 / −0.007% | 0.99 (PSFMODERR 0.02%) |
| `--psf-file` (the true PSF) | −0.009% | | |
| `--psf-var 2` | −0.14% | ±0.9% | 1.0 (PSFMODERR 1%) |
| constant PSF | −1.6% | ±6% | 1.04 (PSFMODERR 4.6%) |
| before: analytic core, 4 FWHM table | −4.4% | +3.1 / −3.1 / −8.1% | 9 |
| before: `--psf-var 2` | −3.4% | | 3.1 |

Gaussian and Moffat grids are unchanged (−0.009%, −0.007%). In the correction
mode the correct base gives the fluxes of the file (0.000% on 16 seeds,
≤ 0.017% on 4 crowded 16-bit seeds); a 10% broader base −0.01%
(+0.008 / −0.030 / +0.007%; before +0.05%: −0.30 / −0.02 / +0.33%).
End-to-end `phcal` on the grid frames: MAG − u = −0.0000 [−0.0006, +0.0002],
scatter 0.0119 mag, the same as with the file (before: scatter 0.059 mag,
−0.076 / +0.049 mag at the centre / edge).

**Run time** on a 2048² frame: 5–7 s with 1700 stars (analytic PSF); with
about 2000 stars 6 s (`EMPIRICAL`), 9 s (order 2 and groups), 13.5 s (and one
residual pass). The normalisation and the radial variation make 2048² frames
1.3–2× slower (3.7–8.7 s).

```{note}
**Limits.** A polynomial in x, y cannot follow a PSF that changes structure
over the field (`dpsf.fits`: order 2 leaves ±0.8% in radial zones); use
`--psf-var-radial` when it changes with the field radius, or the correction
mode with a model PSF as the base. The normalisation needs isolated bright
stars (in a field of 16 px spacing: −0.27% with `--psf-var-radial 1`, pulls
1.0). PSF stars must be mostly single stars.
```

## Growth-curve photometry (`gphot`)

2048² frames, 8 seeds; space-UV frames with a Gaussian PSF and with
`dpsf.fits`, sparse and dense fields; ground-based Gaussian and Moffat frames:

- GCOUNT against the truth within ±0.3% at S/N > 30, and as `psfphot` below
  it (dense `dpsf.fits`: −0.8…+0.1% by field radius, `--radial` −0.5%; the old
  gphot: +5.4…−3.5%);
- pull SD 0.9–1.15 by S/N; 1.24–1.33 for the brightest ground-based stars, as
  `psfphot` (robust 0.85–1.12; the old gphot: robust 1.40 for bright
  ground-based stars, 3.05 with Moffat wings);
- 4096² `dpsf.fits` frames with the `find` defaults: GCOUNT of all 900 stars
  (none before), `phcal --photometry grow` Q 1.01–1.14.

```{note}
**Limits.** At least 3 bright isolated stars are needed (about 30 on shallow
2048² space-UV frames, where a PSF variation over the field is often not
significant: the errors then describe it, and the pull SD at S/N 10–30 is
1.14). The FWHM of `FIND` sets the isolation. SKY of crowded fields holds the
light of faint neighbours in the apertures (+0.3 SKYERR).
```

## Photometric calibration (`phcal`, `phfotran`)

**Systematic error** (`--syserr auto`, P44, P54, P55):

- Monte Carlo of 30–100 stars within 1–6% with up to 10% outliers at
  0.3–3 mag; with 10 stars and 10% outliers the zero point scatters by
  6 mmag (179 mmag with the previous clipped estimator);
- with few stars σ² is unbiased: σ = 0.06 mag from 5 stars gives rms 0.059
  (it was 0.054); σ = 0.03 mag from 10 stars gives 0.028 on average (0.026);
- CTPHERR covers the zero point in 65–70% of sets of 3–100 stars for
  σ = 0.01–0.06 mag (54–63% before for 3–10 stars), conservatively 70–87% for
  σ = 0; with 30 and more stars the errors are 1–6% larger than without the
  uncertainty of σ;
- on the field-dependent-PSF frames above, `--syserr auto` alone finds
  σ<sub>sys</sub> = 0.058–0.062 mag and brings Q from 1.55–1.63 to
  1.05–1.15; on Gaussian frames σ<sub>sys</sub> = 0–0.003 mag.

**Colour terms** (`--color-term`, P47, P55, P56):

- a simulated V frame calibrated against X = V + k (B − V): k = 0.0005,
  0.1007, 0.2996 ± 0.0009 for k = 0, 0.1, 0.3, with no residual colour slope;
- real SDSS colours (cones around RA 180°, Dec 0° and M67): a g frame
  calibrated against r (k = −1 exactly) gives −1.0003 ± 0.0023; a band
  between u and g (k = 0.4293) 0.4295 ± 0.0016;
- a 260–360 nm band against SDSS u: one colour u − g leaves 0.024–0.034 mag
  with the blue end 0.03–0.07 mag off, since the band depends on both u − g
  and g − r. Two terms (`--color-term upmag,gpmag,... --color-term
  gpmag,rpmag,...`) with an exact catalogue give k = −0.862 ± 0.020,
  +0.387 ± 0.041 for the true −0.856, +0.371, scatter 0.013 / 0.009 mag;
  against a catalogue with SDSS errors the coefficients are unbiased, Q
  0.98–0.99. Without the correction for the errors of the colours, two
  correlated coefficients are biased by 4% and 10%;
- full chains (`find` → `aphot` → `astrometry` → `phcal`) on 20–60 s
  space-UV frames: k within 2σ of the truth in 12 of 12 frames, Q 0.84–1.13;
  with two terms on 2048² frames of 15–67 stars, rms pull 0.78 / 0.72,
  Q 0.75–1.22.

Without `--color-term` the outputs are bit-identical to those without the
feature; one `--color-term` gives outputs bit-identical to the one-term code.

**Instrument transformation** (`phfotran`, P45, P48), on simulated B, V, R
frames:

- an instrument whose filters mix 4–13% of the neighbouring bands: the colour
  term −0.164 mag/mag in B drops to −0.0004 ± 0.0002 with `--model flux`;
- an instrument with colour terms linear in magnitudes, rms in B / V / R:

  | Calibration | B | V | R |
  |---|---|---|---|
  | no table | 0.058 | 0.023 | 0.019 |
  | `--model flux` | 0.0136 | 0.0011 | 0.0026 |
  | `--model mag` | 0.0012 | 0.0010 | 0.0008 |

  `--model flux` is still better for a strong mixing of the bands (a leak).

**Physical quantities** (`phcal -q`, P1): FLAM was wrong by ×0.0013, STMAG by
+9.7 mag; they are now correct. MAG was always correct.

## Astrometry

- **Distortion** (`-p TPV`, S19): distortions of 15 px are recovered to
  0.03–0.11 px (before: up to 6 px, 60 px with `--maxstars 1000`);
  undistorted frames give the same solution.
- **Weighted fit** (S27): rms error against the truth 0.099 → 0.055 px.
- **False solutions** (S23): on space-UV frames matched against other parts
  of the sky, 0 of 88 false solutions are accepted (4 before); with the
  largest pointing search, 0 of 192.
- **Footprint catalogues**, 512²–4096² frames: 66 of 81 solved (the rest have
  fewer than five stars; 63 and 5 crashes before), median 0.02 s, at most
  0.6 s; full 0.61° cones: 25 of 30 (14 before).
- **Pointing with an unknown roll** (S29, S31): frames of random roll
  (768²–4096²), the header or `--rcen`, `--dcen`, `--scale` off by 10″–5′,
  rotated at random and scaled by 1–2%: 135 of 140 solved by the guided
  match, 0.03–0.26 s; the true field moved by 15–50% with `--wcs-offset`: 80
  of 84; full 0.61° cones by the pointing: 54 of 56 (38 before).
- **Threads** (S30): a wrong sky of a dense 2048² frame without a scale
  prior takes 54 s with one thread, 11 s with eight; the solution does not
  depend on the number of threads.

## Calibration frames and stacking

- **Master dark** (P31): the old robust mean put a master dark of
  0.13 e⁻/pixel 18% low; the clipped mean is within ±0.0003 e⁻/pixel from
  0.05 to 20 e⁻, with cosmic rays rejected.
- **`kombine`**: the default clipped mean is unbiased at low counts and
  rejects cosmic rays; STDERR is on the scale of the image.

## Simulation (`artificial`)

- Frames made with `--seed` are bit-identical for 1, 2 and 4 processes of
  `--mproc`. Eight 2048² frames with the `SEEING` PSF: 69 s with one process,
  19.5 s with four.
- Noiseless frames made with `--gnoise` and no gradient, airmass or flat are
  bit-identical to upstream hg1653.
- Random star fields: 8 stars on the default 666² frame to the default
  `--maglim 13` (it was 239 on any frame, A23).
