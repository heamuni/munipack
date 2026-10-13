# Munipack hg1653: bugs found and fixed

This file covers upstream Munipack Mercurial revision 1653 (`ae465635f925`, 2026-07-15), reported as version `0.6.2+hg1653`. The untouched upstream tree is commit `003cb0d`, and every fix is a separate commit on top of it. Each bug below gives its measured impact and the commit that fixes it.

How the fixes were checked:
- Each bug was reproduced with the command-line tools in the Docker image, mostly using noiseless or seeded `artificial` frames.
- Results were checked against independent calculations (astropy, numpy and analytic formulas).
- The fix was then shown to remove the discrepancy.
- An end-to-end run (`artificial` → `find` → `aphot` → `astrometry` → `phcal`) gives the same results before and after the fixes on unsaturated frames: astrometry RMS 0.08″, `phcal` MAG − V = −0.016 mag on 156 stars.
- Saturated 16-bit frames, where the pre-fix build failed astrometry, now solve: 0.12″ RMS, 0.05 px from truth, MAG − V = +0.005 on 60 stars.

Status values:
- **fixed**: verified by experiment.
- **fixed (reading)**: the bug was confirmed by reading the code; the fix compiles and passes the regression runs.
- **documented**: behaviour kept, documentation corrected.
- **open**: not fixed.

Build: `docker build -t munipack:0.6.2-hg1653-fixes .` (see [Dockerfile](Dockerfile)).

---

## 1. The two issues from the QUVIK simulation project

### 1.1 `artificial --fov` "6× multiplier"
- **Cause:** a bug in Munipack **0.5.14**, the Ubuntu 22.04 package that the old scripts ran.
  - 0.5.14 converts `--fov` to a pixel scale using the *default* frame height of 666 px and ignores `--width/--height`.
  - A 4000 px frame therefore covered 4000/666 = 6.0× the requested field.
- **Status:** already fixed upstream (in 0.6.x). In hg1651 and hg1653 the header scale equals `--fov`/height, and stars land within 0.000 px of their WCS positions.
- **Remaining details:**
  - `--fov` sets the frame **height**, edge to edge. The width does not enter.
  - `CRPIX = naxes/2` uses integer division, so `--rcen/--dcen` falls on pixel 500 of a 1000 px frame, half a pixel from the geometric centre (500.5). The frame is self-consistent, so this is cosmetic. **open**
  - From hg1652/1653 on, `--scale` is in **pixels per degree** (it was degrees per pixel). The help text was wrong (fixed in `de1236c`). Scripts passing `FOV/(nsize-1)` must pass `(nsize-1)/FOV`.

### 1.2 Passband width ("sigma of the equivalent width")
- **The convention:** the calibration (`mag2ph0`/`ph2mag0` in `lib/photoconv.f08`) treats `LAM_FWHM` as **2σ** of a Gaussian passband. The effective (equivalent) width is therefore √(π/2)·LAM_FWHM = 1.2533·LAM_FWHM.
  - The photosystem table documents the column as an FWHM. A true Gaussian FWHM gives 1.0645·FWHM, so the two readings differ by 17.7% (0.18 mag) in absolute photon rates.
  - Only the wavelength columns (`LAM_EFF`, `LAM_FWHM`, `FLAM_REF`) enter the calibration; the frequency columns are ignored.
  - The convention cancels in every Munipack round trip (`artificial` → `phcal` MAG), so MAG results are correct.
  - **documented**, deliberately not changed, because existing simulated datasets depend on it: `87ce878` adds the description to photoconv.f08, dataform_photosys and dataform_photometry.
  - Suggestion: tabulate a true equivalent width (∫T dλ / T_max), or integrate real transmission curves. SDSS u is far from Gaussian.
- **The real bug was in the other calibrated quantities.** These used *different* widths and wrong formulas (fixed in `87ce878`). Median over 156 stars of an M67 test frame:

  | Quantity | Before | After |
  |---|---|---|
  | MAG − V | −0.013 | −0.013 (unchanged) |
  | FLAM / (FLAM_REF·10^(−0.4V)) | 0.0013 (dimensionally wrong factor h·λ/Δλ) | 1.012 |
  | STMAG − V | +9.745 (wrong factor, and ST zero point 10× too large) | +0.015 |
  | ABMAG − V | −0.062 (no √(π/2), frequency columns) | −0.010 |
  | PHOTLAM | multiplied by the bandwidth (×7000) | divided, consistent with FLAM |
  | FNU, PHOTNU | inconsistent width and columns | consistent with FLAM to 1e-6 |
  | `hstparam` PHOTFLAM | meaningless (2.4e-16) | 3.82e-15 |

---

## 2. `artificial` (frame simulator)

### 2.1 Detector noise model and `--gnoise`
The manual describes `--gnoise` as *"dispersion of gaussian noise per pixel … to simulate dark current, readout noise, noise of background … merged with skymag"*. The implementation does not do that:
- **`--gnoise` raises the sky level instead of adding noise.** It computes `back = sqrt(skyint² + gnoise²)`, adding a dispersion *in quadrature to a mean rate*. The result is then scaled by exposure × area × q.e. × pixel area. So `--gnoise` acts as an extra **background rate** whose noise is √level. It is never the σ given, and it grows with exposure time. Measured: `--gnoise 10` gave mean 9.85 and σ 3.1.
- **No read noise exists anywhere in the code.**
- **The sky noise is a clipped Gaussian, not Poisson** (`max(N(b,√b), ε)`). At 1 e⁻/px the mean is +8% too high, the variance is −25%, 16% of pixels sit at about 0, and values are non-integer.

**Fix (`7172d04`, `46f6935`).** The frame is now built as a mean in electrons, then noised:

*mean = flat · (sky + galaxies + stars) + dark·t* → **Poisson** → **+ N(0, readnoise)** → **÷ gain** → **+ bias** → clip at 0 (integer BITPIX only) → saturation.

| Option | Meaning |
|---|---|
| `--dark` [e⁻/s/pix] | dark current; raises the mean by dark·t, Poisson noise; not affected by gradient, airmass or flat |
| `--readnoise` [e⁻ rms] | Gaussian read noise; does not change the mean or scale with exposure |
| `--bias` [ADU] | pedestal, so that read noise on a low background is not clipped at 0 in unsigned integer frames |
| `--gnoise` | **kept for reproducibility**: noiseless frames are bit-identical to hg1653. Documented as what it really is, and deprecated |
| `--seed` | reproducible frames (see 2.2) |

New header keywords: `DARKCUR`, `RDNOISE`, `BIASLEV`, `ARTSEED`. The verification runs:

| Test | Measured | Expected |
|---|---|---|
| `--dark 0.5 --exptime 20` | mean 9.98, var 9.99 | 10, 10 |
| `--readnoise 5`, exptime 1 and 100 | var 25.06 and 25.10, mean 0.00 | 25, 0 |
| dark 5 e⁻ + RN 7, `-B 16`, no bias | mean 6.10 | biased by the clip at 0 |
| same with `--bias 100` | mean 104.998, var 53.97 | 105, 54 |
| sky ≈ 1 e⁻/px | mean 1.003, var 1.004, integer | 1, 1 |

### 2.2 Other `artificial` bugs

| # | Bug | Impact (measured) | Status | Commit |
|---|---|---|---|---|
| A1 | `scale` changed from deg/pix to pix/deg in hg1652/1653, but sky, galaxy and PSF-file formulas kept the old convention | every pixel ≈ 6×10¹⁵: unusable frames | fixed | `9a0d511` |
| A2 | No `--seed`; noise, random star fields, RANDOM spread and cosmics differ every run | not reproducible | fixed: reseeded per frame from (seed, frame index), so results do not depend on `--mproc` | `7172d04` |
| A3 | Sky gradient: `1.086·grad·back`, wrong sign and linear | +0.002 mag/″²/px requested gave −0.0024 | fixed: exponential, sky only | `7172d04` |
| A4 | RANDOM seeing spread √2 too wide | σ 2.72 vs FFT 2.03 px; PSF jumps 41% at the 3 s AUTO switch | fixed | `8ce40c7` |
| A5 | FFT seeing PSF off-centre by −1 subpixel | star centroids −0.095 px (systematic astrometric bias) | fixed: ≤0.004 px | `24fae93` |
| A6 | Pixel integration double-counts edge subpixels | flux 0.962–1.039 × true for σ = 0.3 px | fixed: trapezoid weights, 1.00000 | `3e178c8` |
| A7 | BITPIX 8 wrote BZERO = 128 (broken); BITPIX 32 BZERO = huge() | `-B 8` unusable; `-B 32` read as float | fixed: uint8/uint16/uint32 correct | `ba7177b` |
| A8 | SATURATE/GAIN written with 3 significant digits | header 65500 for 65535, so clipped pixels not flagged | fixed (15 digits) | `7172d04` |
| A9 | Light curve and airmass evaluated at exposure start; JD is mid-exposure | phase error of exptime/2 | fixed | `7172d04` |
| A10 | Cosmic-ray count drawn once per series | identical count in every frame | fixed: per frame | `7172d04` |
| A11 | Sky attenuated by extinction (darkens with airmass) | 914 → 737 e⁻ at X = 2.08 | fixed: b·X·e^(−k(X−1)), 914 → 1531 | `7172d04` |
| A12 | `--area` ignored for the default Airy radius | Airy radius from the default diameter | fixed | `df71158` |
| A13 | Flat applied after noise | variance ∝ flat²·λ (400 instead of 200 for flat 2) | fixed | `7172d04` |
| A14 | Unknown `--psf`/`--spread` names silently became SEEING/FFT | silent wrong PSF | fixed: error | `de1236c` |
| A15 | Front end passes doubles with `%lf` (6 decimals) | `--lc-per 0.0123456789` → 0.012346; also hit `--gnoise` | fixed (`%.17g`) | `de1236c` |
| A16 | Paired options (`--rcen/--dcen`, `--width/--height`, `--long/--lat`, `--lc-ra/--lc-dec`) silently dropped if only one given | silent | fixed: error | `de1236c` |
| A17 | `--mproc ≤ 0` accepted | `-n 0` | fixed | `de1236c` |
| A18 | `--scale` help said [deg/pix] (artificial, astrometry, kombine) | users following the help got CD ≈ 3600 | fixed | `de1236c` |
| A19 | CATALOGUE `COUNT` is in electrons (before gain/flat) while the image is in ADU | confusing | documented (TUNIT, comments) | `a6db7ca` |
| A20 | Galaxies had no photon noise | — | fixed (part of the new noise model) | `7172d04` |
| A21 | `--hwhm` is σ for GAUSS/SEEING and α for MOFFAT, not HWHM | GAUSS HWHM is 1.177× the input; MOFFAT β=2 is 0.64× | documented (behaviour kept; existing pipelines pass σ) | `7172d04` |
| A22 | Watermark added after the saturation clip, in electrons on an ADU frame. With integer `-B` and gain < 1, a bright star under the watermark (top-right corner) exceeded the integer range: FITSIO "datatype conversion overflow", then **SIGFPE and a truncated file**. Random: `--gain 0.5 -B 16` crashed in about half of the unseeded runs (also in the unpatched build) | crash, truncated output | fixed: watermark added in electrons before the conversion, so saturation applies to it; 0 failures in 20 runs that crashed before. The watermark can be switched off by `--no-watermark` (A31) | `83096f0` |
| A23 | `genstars` (no catalogue) applies a per-square-degree star count to the whole frame regardless of its size | wrong density for random fields: 2191 stars to `--maglim 15` on every frame, i.e. 64000 per deg² on the default 666² frame and 6770 per deg² on 2048² at 1″/px, instead of 2191 | fixed: the count of every magnitude bin is Σ(m) × width·height/scale² [deg²], rounded up with the probability of its fractional part; 8 seeds at 1″/px: 512² 43.9 stars (2169 ± 17 per deg²), 1024² 176.9 (2186 ± 5), 2048² 709.4 (2192 ± 2), 1024² at 2″/px 709; Σ(m) per bin recovered (m = 12–15: 54.8, 161.8, 484.3, 1458.8 against 54, 162, 487, 1465) | this branch |
| A24 | `wxLogFatalError` (the project's usual error call) exits with code 139 | looks like a crash: every error found by the front end (paired options, `--mproc 0`, unknown `--psf` of aphot/psfphot, `@missing.lst`, `-t /nonexistent`, `--vocat missing.xml`, …) gave 139; `munipack phfotran` (listed, but never dispatched) crashed on a null pointer; failures of the utilities gave 13, an invalid command line 255, and also `--help` and `--version` 255 | fixed: one helper (`MuniFatalError`) prints `Error: …` and exits with 1; the catalogue configuration (`VOCatConf`) reports an error instead of aborting; utility failures (`STOP` with a message, e.g. a missing FITS file, `artificial --psf unknown`) and an invalid command line or unknown action exit with 1, `--help`/`--version` with 0; a write error on the input pipe of a utility (not reproduced: a utility that does not read its input exits 1) stops it instead of aborting; 25 error cases checked | this branch |
| A25 | `--mproc` could not be tested: the Docker build has no coarray runtime (`-fcoarray=single`) | — | **open** (build) | — |
| A26 | `--gnoise` was merged with the sky into one level, so the sky gradient and the airmass scaled the `--gnoise` part too | `--sky-mag 21 --gnoise 1000 --sky-grad-x 0.001`: whole 986 e⁻ level scaled by the gradient | fixed: the merged level is split into the sky part (gradient, airmass) and a uniform `--gnoise` part; same total without gradient or airmass (bit-identical); left−right difference now −20.56 e⁻ with or without `--gnoise` | `932ac64` |
| A27 | `--psf-file`: the loader put the PSF centre at the sample NAXIS/2+1, while the PSF generators (`doc/artpsf.f08`, the field-PSF generator of `dpsf.fits`) and `CRPIX` of frames use NAXIS/2 | every star rendered with `dpsf.fits` offset by −0.50 px in x and y from its catalogue (WCS) position | fixed: the centre is NAXIS/2 or `CRPIX1/2` of the extension; centroid − catalogue position 0.0002 px | `c72d8a3` |
| A28 | Stars were pixelised at positions rounded to the PSF sample grid (`nint(zoom*dx)`) | position errors up to 1/(2·ZOOM) px: ±0.25 px with `dpsf.fits` (ZOOM=2), ±0.031 px (rms 0.018) for GAUSS/MOFFAT/SEEING (ZOOM=16) | fixed: exact pixel integral of the interpolated samples at any position (trapezoid/midpoint rule on the grid, so grid-aligned stars are unchanged: 1 of 262144 pixels differs by one float32 ulp); rms 0.0004 px | `c72d8a3` |
| A29 | `psfinterpol` supports one ring of exactly eight angles: `l = min(int(r/rgrid),0)` is always 0, the indices 2..9 are hard-coded, the loader reads 360/AGRID+1 extensions in file order | other grids (more rings, other angles, extensions in another order) silently wrong | fixed: N rings × M angles from RADCOO/ANGCOO (linear in angle with wrapping and in radius, centre PSF inside the first ring, outermost ring beyond), stamps of different sizes, other extensions ignored; agrees with an independent implementation to 2×10⁻⁵ of the peak (2 rings, 8 and 6 angles, ZOOM=3) | `c72d8a3` |
| A30 | PSF file documentation: described a single-image format the loader rejected (it required RGRID, AGRID) and recommended a sum of 1, while the loader takes a sum of ZOOM² as the whole flux; the template `doc/artpsf.f08` normalised to 1 | a template PSF (ZOOM=4) gives stars 16× too faint | fixed: format and normalisation documented (sum = ZOOM², a warning otherwise), the single-image format is accepted, the template normalised to ZOOM² | `c72d8a3`, this branch |
| A31 | The watermark of `artificial` (A22, top-right corner) is a real image feature that `find` detects as about two "stars" per frame; the first version of the `psfphot --residual-detect` added two more per frame there | spurious objects in the corner of every simulated frame | fixed: `--no-watermark` leaves it out; when drawn, its region is recorded in WMARKX1/X2/Y1/Y2 (first and last pixel); the default is unchanged (watermark on, image data bit-identical for a fixed `--seed`, 16-bit and float). QUVIK-like 1024² frames, `find -th 4 -f 4.9`, 4 seeds: 40 detections in the box with the watermark, 0 without (frames smaller than 44×40 px, which crashed with an empty file writing the mark outside the array, get no mark). `--residual-detect` no longer adds them (the local-scatter test, P20): 16 → 0 in 8 frames | this branch |
| A32 | `artificial` added only the positive pixels of the pixelised star (`if( f > 0 )` in `spreadstars`), so a PSF with negative values (measured or empirical PSFs: noise and background subtraction in the wings) gained flux | a frame rendered with an empirical PSF of `psfphot --psf-save` (wing samples of ±1e-5): a halo of +0.14% of the flux at 6–14 px; fitted by the same PSF −0.13% low against the frame, the empirical PSF of the frame recovered the halo | fixed: all pixels are added (the mean image is clipped at zero only for the noise, as before); the halo is gone (FILE −0.017%, EMPIRICAL −0.009%); frames with non-negative PSFs (all analytic, `dpsf.fits`) are bit-identical | this branch |

---

## 3. Photometry and calibration

| # | Bug | Impact (measured) | Status | Commit |
|---|---|---|---|---|
| P1 | FLAM, STMAG, PHOTLAM, FNU, PHOTNU, ABMAG formulas (see 1.2) | FLAM ×0.0013, STMAG +9.75 mag, PHOTLAM ×7000 | fixed | `87ce878` |
| P2 | `phcal -q PHOTRATE` (advertised) rejected: code expects `PHRATE` | crash, leaves `fitsio_scratch_*.fits` | fixed: both accepted, validated before writing | `ede0662` |
| P3 | `phcal` with several frames scales every frame's reference photons by the *product* of all exposure × area; catalogue errors never scaled | joint B,V gave ctph 6.346/3.174 instead of 1.058/1.058, then SIGSEGV | fixed: per frame; the crash came from 13-element system arrays overrunning frame arrays | `369b83b` |
| P4 | `phcal` manual calibration (`-C`): transformation uninitialised | inf/0 images, all MAG 99.99 | fixed | `36de6e6` |
| P5 | `gphot` GCOUNTERR overwritten by √counts | errors 1.5–3.7× too small (reported/observed 0.43) | fixed: from the aperture covariance, reported/observed 0.92–1.13 | `4458a13` |
| P6 | `gphot` sky correction subtracted instead of added; error divided by √2 | a sky error of +5 moved SKY by +9.9 instead of removing it | fixed: +0.09 | `bfa820a` |
| P7 | Growth-curve profile divided by 2π instead of π; area evaluated at the outer radius | GROWPHOT HWHM 6.01 for a true 2.355 | fixed: 2.32 | `b1c9f60` |
| P8 | APERPHOT `HWHM` keyword stores σ; brightest-star sample written to the wrong slot | HWHM 2.02 for a true 2.35 | fixed: true HWHM (see compatibility notes) | `90d8ff9` |
| P9 | `grow_phot` picks the aperture by radius, not by distance from HWHM; passes the uncompacted curve | GCOUNT bias +3.5% when inner apertures are invalid | fixed: +0.2% | `c29a59e` |
| P10 | `aphot --coo` low-bad threshold `max(0, sky−3σ)` | large apertures nulled in 6/8 bright stars; sky +24 on sky-subtracted frames | fixed: sky − 5σ, no clamp | `9cf6933` |
| P11 | `phcorr` median-mask box `min(box/2,1)` ignores `--box`; error propagation wrong (no ÷flat², Poisson from flat-fielded values, no read noise, DOUT = 0 for OUT ≤ 0) | 3×3 box always; error ratio 0.0–1.24 | fixed: box honoured, ratio 0.99–1.03, READNS read from header | `65bf89c` |
| P12 | `sfits`: every frame got the sky of the last file; RF90 ×scale instead of ÷scale (pix/deg); CRDER never read | SKYMEAN 278.6 for both B and V (202/277); RADFLX90 12630 "deg"; 18/79 stars matched | fixed: 79/79 | `34076f1` |
| P13 | `find`/`aphot` accept pixels *equal* to the saturation level as valid (`datum <= hibad`), while `daosky` uses `<` | saturated stars detected and measured on real 16-bit data (pixels read exactly 65535); 311 instead of 136 detections on a test frame | fixed | `5d96c1d` |
| P14 | `readframes` sorts frames by the photometric system's filter order, while catalogue columns follow the user's `-f` order | possible mismatch when `-f` is not in system order | **open** | — |
| P15 | `pphot` (2013 stub, built but no front-end action) did nothing | no PSF photometry | replaced by `psfphot` (see README, PSF photometry) | `f6afeda` |
| P16 | The `aphot` sky (robust mean in the 20–30 px annulus) includes the light of neighbours | +0.4 ADU on a 500-star Gaussian test field, +2.7 ADU with Moffat β=2.5 wings: −0.3% and −1.6% flux at S/N 50 when used as a fixed sky | `psfphot` re-estimates the sky on the frame with fitted stars subtracted (≤0.1%); `aphot` now masks the detected neighbours in the ring (P23) and invalidates apertures contaminated by brighter neighbours (P29); the light of undetected stars stays in the (mean) sky, as it is in the apertures | `53828cb`, `d805f61`, `a0f570d` |
| P17 | (psfphot, found in review) Saturated stars are not in FIND (find excludes saturated pixels, P13), so their light was not subtracted from neighbours | crowded 16-bit field with 36% saturated stars: flux errors underestimated, pull robust std 1.37 (fixed sky) / 1.16 (local); 5.2 for stars next to an undetected saturated one | fixed: saturated areas without a FIND star are added as stars (FLAG 64), fitted from the wings and subtracted; pull 0.92–1.06 | this branch |
| P18 | (psfphot) The Moffat shape fit on stars with Gaussian profiles runs β towards infinity (median β = 48 on a σ = 1.5 frame): every fit exhausts its iterations at ZOOM 16 | `--psf MOFFAT` on a 512² frame with 244 Gaussian stars: 254 s (GAUSS: 0.7 s); flux bias +0.45%, pull std 1.5 (a Moffat profile cannot become Gaussian) | **open** for `--psf MOFFAT` (behaviour kept: capping β at 10 is 19× faster but the bias grows to +0.99%; use `--psf GAUSS` or `--beta`); the core of the empirical PSF caps β at 10 and is fitted at the table zoom (the tables correct the profile): 1.5 s | this branch |
| P19 | (psfphot) Iterative neighbour subtraction (one star at a time, three passes) does not separate stars closer than about 1.5 FWHM: the passes do not converge to the joint solution | equal pair at 1 FWHM: +26% / −19%, pull std 25; at 1.5 FWHM: +0.40%, pull std 2.8; crowded field (2500 stars on 1024²), stars with a neighbour ≥ 0.1× within 2 FWHM: scatter 1.32%, pull std 2.44 | fixed: stars closer than 1.5 FWHM (`--group`) are fitted simultaneously: 1 FWHM −0.01% / +0.03%, pull std 1.0; 1.5 FWHM +0.01%, 1.04; crowded 0.87%, 1.76; `--group 0` reproduces the old results bit for bit | this branch |
| P20 | (psfphot) Only the stars of `find` are fitted: `find` merges a star with a neighbour closer than about 1 FWHM (and misses faint companions), whose light then goes to the fitted star | crowded field: stars with a companion ≥ 0.1× within 1 FWHM +17.5% (scatter 16%); pairs at 0.5–0.75 FWHM: 0% of the companions detected, 1:3 at 1 FWHM 6% | fixed (opt-in `--residual-detect N`): detection in the residuals against the model noise with the measured PSF mismatch, scaled by the observed local scatter; +0.24% (scatter 1.2%), 88–100% of the companions at 0.5–1 FWHM found (19% for 1:10 at 0.5 FWHM), 0 spurious detections in 7 × 3200 stars | this branch |
| P21 | (psfphot) A PSF file name longer than a FITS card (e.g. an absolute path) was copied into an 81-character buffer for PSFFILE | SIGSEGV after the whole fit, no PSFPHOT written (`--psf-file` with a 110-character path) | fixed: the keyword values keep the last 64 characters with a leading "..."; long comments of saved PSF files are split | this branch |
| P22 | (psfphot) The correction of a correct PSF file (`--psf EMPIRICAL --psf-file`) followed noise and frame systematics: a wing pedestal of ~0.05% of the flux per zone, amplified by the fixed-sky feedback | crowded 16-bit `dpsf.fits` frame (seed 21): −0.27% (pull median −1.27) against −0.007% by `--psf-file`; other seeds up to ±0.25% | fixed: significance shrinkage of the correction by 8 subsets of the PSF stars with a 0.05% floor per radial zone; 16 seeds: |correction − FILE| < 0.001%, pull median within ±0.09 | this branch |
| P23 | (aphot, find, psfphot) The sky of a star (aphot), of the frame (find) and of psfphot's annuli was a robust mean (rmean, Hampel), a median or a mean clipped at 3σ by the MAD; a sky of a few photo-electrons per pixel has a skewed, discrete (Poisson) distribution whose robust means lie below its mean, so every aperture kept bias × area of sky | QUVIK-like 2048² 16-bit frames (gain 0.101 e⁻/ADU, read noise 0.37 e⁻, 4–8 seeds): SKY − true background −0.050 e⁻/px (−13%) at 0.39 e⁻/px, −0.029 (−29%) at 0.1, −0.066 (−3.3%) at 2, −0.065 (−0.3%) at 20 e⁻/px; photometry at the catalogue positions (u < 18.5): aperture 7 offset +0.050 ± 0.009 e⁻/px, flux +5.3% at u = 17–18, +11.7% at 18–19, +23% at 19–20; SKYERR 0.79× the scatter between seeds; psfphot +41 ± 9 ADU; find's sky σ 4% low (7.05 instead of 7.30 ADU; 28% low, 5.26, on frames flattened by the user's pipeline): threshold 4 acted as 2.9σ, 8600 detections per frame | fixed: `skymean.f08`, the mean clipped at *t* − 5σ … *t* + *q* + √(*q*² + (5σ)²), *q* = 25/6 e⁻ (the Bernstein bound of the Poisson tail: no photon of the sky clipped even when σ < 1 e⁻), σ the larger of the observed scatter and the noise model √((*t* − BIASLEV)/GAIN + RN²); aphot masks the detected neighbours (a radius per star following its light). SKY − truth within ±0.0013 e⁻/px from 0.1 to 20 e⁻/px (−0.3%…+0.07%; 0.07 ± 0.04 e⁻ of 2000, gain 1.5: ±0.05 e⁻ of 200–2000); aperture offset −0.003 ± 0.009 e⁻/px (Gauss), +0.024 ± 0.011 (dpsf), −0.007 ± 0.010 (0.1 e⁻/px), −0.026 ± 0.028 (100 s); flux +0.7 ± 0.6% at 17–18, +1.6 ± 1.1% at 18–19, −0.9 ± 2.5% at 19–20; SKYERR/scatter 1.01–1.10 (8 seeds); psfphot +10 ± 9 ADU; end-to-end (the user's aperture route + phcal) MAG − u at 18–19 mag −0.165 → −0.057 (dpsf), −0.133 → −0.061 (Gauss, the rest is the detection bias at the limit; PSF photometry −0.030). find: the threshold is now in true σ (fewer noise peaks: 852 → 578 per Gaussian frame, real stars −5% at the limit). `--apcorr auto` no longer removes EESKYBIA (now 0.05–0.10 ADU/px, a diagnostic); APEE unchanged | `d805f61`, `e25280f` |
| P24 | (gphot) On frames with the field-dependent `dpsf.fits` PSF the growth curve is built from 4 stars, HWHM 0.80 px, and no star gets a valid GCOUNT | `phcal --photometry grow`: "Catalogue to frame joint stars not found" on all four test frames (the Gaussian frames work) | **open** | — |
| P25 | (find, seen in the calibration) Bright stars of the `dpsf.fits` frames (a 1.5 px core with a halo) get secondary detections 4–15 px from the centre; their apertures hold 70–97% of the star | 2–3 detections per frame matched to faint catalogue stars, 2.6–4.7 mag too bright; they dominate the unclipped error-calibration statistic Q (32.9 → 1.6 without them, see README). Cause: the halo carries the photon noise of its light (10× the sky noise at gain 0.1), its curvature adds to the convolution, and the threshold is set by the sky noise while the local-maximum test reaches only 0.637 FWHM: noise peaks in the halo pass. Some of the outliers were real faint stars in the wings (a u = 17.07 star 6.2 px from a u = 11.4 star, 4.4 mag too bright): see P29 | fixed: a candidate within 8 FWHM (20 px) of a candidate brighter twice, or of a saturated area, must exceed the convolution of the wing at its distance (annulus median) by the threshold times the noise of the wing (the sky noise raised by the photon noise of the wing light with GAIN, or the annulus scatter); 2048² frames, the user's find (auto FWHM, threshold 4), 8 seeds: detections near bright stars 9.25 → 1.88 per dpsf frame (half of the rest at the frame border), 1.00 → 0.62 Gaussian, 4.75 → 2.00 with saturated stars; companions at 1, 1.5, 2, 3 FWHM, 1:1–1:30, u = 12 and 14, Gauss and dpsf, 4 seeds: completeness identical in every cell, spurious near the pairs 879 → 190 (dpsf); end-to-end with P29: unclipped Q of the aperture route on 4 dpsf frames 56/60/0.99/47 → 1.11/1.02/0.99/1.08 (aperture + APEE by the file 0.90–0.93, `--apcorr auto --apcorr-radial` 0.90–0.94, psfphot 0.87–0.92), no outlier | `790c6a7`, `a0f570d` |
| P26 | Front end: `aphot --eccentricity/--inclination` passed with `%lf` (6 decimals), `phcal --tol/--area/-th/-e` with `%e` (7 digits) | `--tol` 1e-7 deg granularity | fixed (`%.17g`) | this branch |
| P27 | (aphot `--apcorr auto`) The warning about a lowered order of the EE variation wrote a character item with an `i0` descriptor | "Fortran runtime error: Expected INTEGER for item 4": aphot failed on every frame with too few bright stars for the requested order (all QUVIK-like 2048² frames with `--apcorr-radial`) | fixed | `4bb0f1c` |
| P28 | (find, aphot `--coo`) `daosky` packed every nstep-th pixel with a full-size mask (`pack(d(::n,::n), d < hibad)`) | for nstep > 1 (frames of 1.6·10⁷ pixels and more in find, e.g. 4096², and of 1.6·10⁵ in aphot `--coo`) the saturated pixels were selected by a mask of a different shape | fixed (P23 rewrite) | `d805f61` |
| P29 | (aphot) A faint star a few pixels from a bright one (or a second detection of the same star) kept the light of the bright one in its apertures | dpsf frames: a u = 17.07 star 6.2 px from a u = 11.4 star 4.4 mag too bright (2 of 4 frames, Q 56–60); Gaussian frames with the user's FWHM (2.7 for 4.9 px): twin detections 1.3–2.7 px apart | fixed: an aperture whose light from a brighter neighbour (by its growth curve) exceeds 3× its error (and 3× the error of the estimate) is −1; a saturated neighbour invalidates the apertures reaching its mask; no outlier left (P25); 0–1 fewer calibrated stars per frame | `a0f570d` |
| P30 | (kombine) The default robust per-pixel mean of the stack is biased on frames of a few electrons per pixel, and differently for sky and stars (non-linear) | 8 QUVIK-like frames (0.39 e⁻/px): sky −0.033 e⁻/px (−8.5%); `--arith` exact. The background subtraction by a robust mean of the frame leaves a constant (removed by aphot's sky) | **documented** (man page: use `--arith`); not changed: the stack has a few values of reprojected, rescaled pixels whose noise model is not carried through | `b72f985` |
| P31 | (dark) The robust per-pixel mean of dark frames is biased on a dark of a few electrons per pixel | 10 QUVIK-like dark frames (0.13 e⁻/px): master dark −0.024 e⁻/px (−15% of the dark signal) against the arithmetic mean | **documented**, open (no arithmetic option) | `b72f985` |
| P32 | (find) The sharpness test (0.2–1.0) uses the central pixel, whose Poisson noise at low counts moves the sharpness of real stars out of the limits | QUVIK-like Gaussian frame, the true FWHM 4.93: u = 15–16 (500–1000 e⁻) 10/15 found, 16–17 3/7; with `--shrplo 0 --shrphi 2` 15/15, 6/7 (+40% noise detections) | **open** (a noise-aware limit would need the gain in the shape tests) | — |
| P33 | (find) The threshold assumes Normal noise; the convolution of a sky of a few photons per pixel is skewed, so far more noise peaks pass | QUVIK-like 2048² frames, threshold 4: ~520 noise detections per frame with FWHM 2.7, ~200 with 4.9 (≈10 expected for Normal noise); the user's measured FWHM (2.5–2.8 px for a 4.93 px Gaussian PSF) doubles them | **open** (a Cornish–Fisher threshold by the skewness of the Poisson model is a candidate) | — |

---

## 4. Astrometry and `kombine`

| # | Bug | Impact (measured) | Status | Commit |
|---|---|---|---|---|
| S1 | Rotation gradient and Hessian of the affine fit ignore the reflection flag; Huber Hessian weights swapped | **every mirrored frame failed to solve** | fixed: mirrored frames at 0–179° solve with 0.12–0.15″ RMS | `08d2756` |
| S2 | `pix2uv` not the inverse of `uv2pix` for reflected, rotated frames | wrong cross-matches; written CD with the rotation sign reversed | fixed | `429a4af` |
| S3 | `get_affines` ignores the CDELT sign; PC matrix read transposed (`wcs_get_double2d`); reflection lost at 90°; rotation averaging wrong near ±180° | CDELT/PC headers treated as mirrored (kombine output mirrored E–W), 1.7–2″/px errors | fixed: exact | `5bd3980` |
| S4 | TPV coefficients fitted in the flipped (u = −ξ) frame but written unconverted | median residual 2.26 px with `-p TPV`, worse than plain TAN | fixed: 0.056 px | `20fe190` |
| S5 | `--mode=cross`: stale deg/pix scale and tolerance, angle and reflection never set | never matched | fixed: 0.17″ RMS | `e381f4d` |
| S6 | Failed RMS check left the solution marked reliable | SIGSEGV | fixed | `18cd173` |
| S7 | RA residuals without cos δ | RA RMS 0.348″ vs true 0.119″ at Dec +70° | fixed: equal to astropy | `d1045da` |
| S8 | Matching chunk off by one | last star of each chunk dropped | fixed | `db2271c` |
| S9 | `kombine` samples the output grid at pixel corners (`i − 0.5`) | stars offset +0.5 px from the output WCS | fixed: ≤0.02 px | `65387a9` |
| S10 | Bicubic: wrong stencil in `bi3conv`; wrong cross derivative | −0.29/−1.02 px shift; rms error 0.034 on a smooth image | fixed: −0.015/0.000 px, rms 0.018 | `a63187f`, `ededddb` |
| S11 | `kombine` never resamples edge pixels; counts out-of-frame pixels as valid; fill value averaged over empty pixels | border rms error 27; garbage values 0–128 | fixed: 0.45 | `ededddb`, `a84b702` |
| S12 | `kombine --pm-ra/--pm-dec` has no effect | proper motion ignored | fixed: shift 3.518/−1.782 px vs 3.524/−1.800 expected | `ededddb` |
| S13 | `kombine` output centre uses integer `naxes/2` with the first frame's CRVAL | output shifted, 60% empty when CRPIX is off-centre | fixed | `9a5306d` |
| S14 | `kombine` output is surface brightness: star flux scales by (s_in/s_out)² when rescaling | 0.24× flux at 2″/px from 1″/px | documented; new `--flux-conserve` switch (0.998) | `0325a92` |
| S15 | Two modules named `astrotrafo` (old `lib/astrotrafo.f08`, unbuilt `lsfit.f08`, `absfit.f08`) | accidental-linking hazard | marked with comments | `37e5787` |
| S16 | Saturated frames: astrometry accepted a **false** 6-object solution (scale 0.54×, rotated 50°, mirrored, 82 px off). The old build rejected the frame only because mirrored fits could not converge; no test checked whether a match was real. The saturated brightest stars are missing from the frame, so the first catalogue chunk can only match by chance | wrong WCS, `phcal` then matched 3 stars | fixed: Poisson significance test of the match count against chance coincidences (disabled by `--disable-rms-check`), plus a retry with the next star chunks. Now 120 objects, 0.12″ RMS, 0.05 px from truth; `phcal` 60 stars, MAG − V = +0.005 | `3341812`, `8871072` |
| S17 | `munipack astrometry -m cross` is parsed as the `cross` action (`--mode=cross` works): the front end took the last argument that is an action name, and removed every argument equal to it | `-m cross`, `--mode cross`: "Unknown option" (exit 255); a file named like an action was dropped | fixed: the action is the first argument that is an action name, and only it is removed; `-m cross`, `--mode cross`, `--mode=cross`, `-m CROSS` run the cross mode (0.107″ RMS on an M67 frame, as `--mode=cross` before), the match mode is unchanged; `munipack cross t1.fits find` passes the file `find`. | this branch |
| S18 | `wcsreport`/`--report` print debug output; index can overrun | — | **open** | — |
| S19 | Large TPV distortions (≈8 px) on rotated frames fit poorly at field edges | max 10 px | **open** | — |

---

## 5. Compatibility notes (behaviour that changed)

- **APERPHOT `HWHM` is now a true HWHM** (1.1774σ), previously σ. New `aphot` + new `phcal` pick the same calibration aperture as before. But new `phcal` on frames processed by the **old** `aphot` reads σ as HWHM and picks a smaller aperture: on one test, ctph went from 1.057 to 1.175. **Re-run `aphot`** on old frames, or set `--saper`.
- **`phcal` MAG shifts by 0.2–0.4 mmag**, because catalogue errors are now scaled correctly and the weights change (P3).
- **`artificial`:**
  - noise statistics changed: Poisson sky, galaxies noised, flat applied before noise;
  - the sky gradient is a different (correct) function;
  - the sky now rises with airmass;
  - RANDOM spread is narrower;
  - FFT PSF centring changed by up to 0.1 px.
  - stars are rendered at their exact sub-pixel positions (A28): changes of up to 1/(2·ZOOM) px, i.e. 0.031 px for the analytic PSFs; stars on the sample grid are unchanged;
  - `--psf-file`: with files centred at NAXIS/2 (such as `dpsf.fits`) the stars move by +0.5 px in x and y to their catalogue positions (A27). A file centred at NAXIS/2+1 keeps the old placement with `CRPIX1/2` = NAXIS/2+1 in its extensions.
  - Noiseless frames with `--gnoise` and no gradient, airmass or flat are bit-identical to hg1653.
  - Random star fields (no `-c` catalogue) have the density of the documented Σ(m) per square degree (A23), so their number of stars follows the frame area: to the default `--maglim 13`, 8 stars on the default 666² frame at 1″/px instead of 239, 77 on 2048² (raise `--maglim` for more). Random fields generated with the same `--seed` differ from older versions.
- **`phcal` adds a systematic error by default** (`--syserr auto`, this branch): the excess scatter of the calibration stars over their errors is added in quadrature to every MAGERR (PHSYSERR in the header) and CTPH is refitted with it. Frames whose stars scatter as their errors say get 0 (or a few mmag) and unchanged results; frames with unmodelled systematics (a field-dependent aperture loss) get larger, honest errors and CTPH changes within its error. `--syserr 0` restores the old errors exactly (verified bit-identical).
- **`phcal` divides aperture counts by APEE when APERPHOT has it** (`--apcorr auto`, default; `aphot --psf-file/--psf/--apcorr` write it). Old files have no APEE and give the same results; `--apcorr none` ignores it; a numeric `--apcorr` keeps its old meaning (a constant factor) and disables APEE. CTPH then refers to total counts.
- **The sky of `aphot`, `find`, `psfphot` is the clipped mean** (P23): on frames of a few photo-electrons per pixel SKY rises by up to 0.05–0.07 e⁻/px and APCOUNT of faint stars drops (−5…−20% at the detection limit); SKYERR is larger by up to 25%; on ground-based frames changes are within the noise. `find`'s threshold is now in the true scatter of the sky (it was a robust scale, too low on discrete data): fewer noise detections, a few % fewer stars at the limit; lower `-th` to keep the old depth. `aphot --apcorr auto` no longer subtracts EESKYBIA.
- **`find` rejects candidates in the wings of brighter stars** (P25): fewer detections around bright, haloed and saturated stars.
- **`aphot` invalidates apertures contaminated by a brighter neighbour** (P29, −1 as other invalid apertures): faint stars next to bright ones and twin detections drop out of phcal.
- **Exit status** (A24): 0 on success, 1 on every error (it was 139 for the errors found by the front end, 13 for the failures of the utilities, 255 for an invalid command line); `--help` and `--version` return 0 (were 255).
- **`astrometry` rejects statistically insignificant matches.** Solutions with few matched stars that the old code would have accepted can now fail; `--disable-rms-check` switches the test off.
- **`--scale` is pix/deg** in `artificial`, `astrometry` and `kombine` (an upstream change in hg1653; only the help text was fixed here).

---

## 6. Suggested further improvements (not done)

**Speed**
- `artificial`:
  - Recomputes the PSF for every frame, even without `--atmosphere`.
  - The FFT box is far too large (513² × 289 operations per star at hwhm 2). Caching the PSF and precomputing a 16×16 grid of sub-pixel-phase kernels would give orders of magnitude.
  - Poisson sampling uses Knuth's method (O(λ) draws per pixel for λ < 500); PTRS (Hörmann) is O(1) and exact.
- `jam.f95 cooid` and `growcurve select_valid_isolated` are O(N²) with an `acos` per pair; a k-d tree or declination-sorted search would remove that.
- `sfits.f08` transformation loops run in cache-hostile order.

**Precision**
- `kombine` point-samples, so downsampling aliases; area-weighted resampling (drizzle) would also make the output flux-conserving by construction.
- Use one magnitude-error constant 2.5/ln 10 everywhere (photoconv is done; phsysfits still uses 1.0857/1.086).
- Store a true equivalent width in the photosystem tables (see 1.2).
- The Airy-pattern centre value is 1, but the limit of (J1(r)/r)² is ¼ (negligible after convolution).

**Refactoring and dead code**
- `photoconv.f08`, about a third unused: `squantity`, `ctrate`, `ct2ph`, `cts2flux0`, `ph2fllux`, `ph2mag`, `x2mag`, `fnu2abmag`, `flam2stmag`, `mag2ph` (which also has a hard-coded 5-element array), `mag2ph3`, `ph2mag3`, and `flux2flam`, `flam2flux`, `ph2flam`, `flam2plux`, `tridig`, which only those call.
- `growcurve.f95`/`growmodel.f95` (5,200 lines): `goto 55` past a hard-coded curve, `do i=1,nid*0` loops.
- `phcal.f08`: a `do n=1,0*nfiles` loop.
- `gphot`: leftover debug output.
- `fits_aphot_image` duplicates `fits_aphot_read`; `sfits tra_frame` and `growfit snoise_graph` are unused.
- The old `astrotrafo` module and the unbuilt `lsfit.f08` and `absfit.f08`.
- The CATALOGUE/FITS column read loops use `l = min(i+srows, nrows)`, which reads overlapping chunks (harmless; should be `i+srows-1`).
