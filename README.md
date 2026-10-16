# Munipack (patched hg1653) — command-line build

A fixed and extended copy of [Munipack](https://munipack.physics.muni.cz/), the astronomical image-processing package by Filip Hroch (Masaryk University). It contains:
- upstream Mercurial revision **1653** (`ae465635f925`, 2026-07-15, reported as `0.6.2+hg1653`);
- about 60 bug fixes, a physically correct detector-noise model for `artificial`, and a Docker build of the command-line tools.

Every bug, its measured impact and the commit that fixes it are listed in **[BUGS.md](BUGS.md)**.

## Repository layout

| Path | Contents |
|---|---|
| `src/` | Munipack source: upstream hg1653 plus the fixes |
| `deps/oakleaf/` | OakLeaf robust-statistics library, upstream hg snapshot (unmodified) |
| `Dockerfile` | Builds the command-line tools (no GUI) on Ubuntu 24.04 |
| `tests/` | Tests of the installed tools on simulated frames (`tests/run.sh`) |
| `.github/workflows/ci.yml` | CI: builds the image and runs the tests on every push to `main`/`develop` and on pull requests |
| `BUGS.md` | Bug catalogue, compatibility notes, suggested improvements |

The history starts with the untouched upstream tree (commit `003cb0d`), and each fix is a separate commit on top of it. `git diff 003cb0d -- src` therefore shows every change to Munipack, and `git log` gives the reason for each one.

## Build

```bash
docker build -t munipack:0.6.2-hg1653-fixes .
```

What the build does:
- Downloads and checks wxWidgets 3.3.3 and builds only its non-GUI part. This revision needs the wxWidgets 3.3 API; Ubuntu ships 3.2.
- Builds OakLeaf from `deps/`, then Munipack, with OpenCoarrays (Ubuntu's `libcoarrays-openmpi-dev`) for `--mproc`.
- Produces an image of about 450 MB, most of it Open MPI for the parallel processes (`--mproc`) and the libraries its packages depend on (UCX with the AMD GPU runtime, hwloc plugins), which Munipack does not use.

The first build takes roughly 10 minutes; later builds reuse the cached dependency layers.

## Tests

```bash
docker build --target test .          # build, upstream make check, then tests/run.sh
tests/run.sh [t_name ...]             # the same against the munipack on PATH (KEEP=1 keeps the work files)
```

- The build runs the upstream test suite (`make check`: find, aphot, astrometry, phcal on the bundled Barnard frame); a failure stops the build. Its pass condition is CTPH within twice its error of 11.5 (it was 11–12, narrower than the error of 1.3).
- `tests/` has a test per area, each on frames that `artificial` simulates from a synthetic catalogue (one from 93 real SDSS stars embedded in the test), so the truth is known and no network is needed: `t_cli` (exit status, rejected options), `t_artificial` (seeds, the noise model, `--mproc`), `t_pipeline` (find → aphot → astrometry from scratch → phcal, and psfphot → phcal against the catalogue, with the pulls of the errors), `t_astrometry` (rotated frames at Dec +60), `t_kombine` (flux and positions of a stack), `t_crowding` (apertures beside brighter neighbours, P29/P46), `t_calib` (bias, dark at 0.13 e⁻/px, phcorr), `t_phcal_filters` (two filters in any order, no colour term), `t_colterm` (the colour term of one filter, phfotran's magnitude model and phcal `--tratab` with it), `t_colsdss` (the colour term with real SDSS colours, k = −1 exactly), `t_gphot` (GCOUNT and its errors against the truth, phcal of GCOUNT), `t_syserr` (the systematic error of few calibration stars, PHSYSERE), `t_psfempir` (the normalised empirical PSF with the radial variation, PSFMODERR). The comments of each test name the bugs of BUGS.md it covers. The whole set takes 2–3 minutes.
- CI runs on the group's self-hosted runner (`varda`): every push to `main` or `develop` and every pull request gets a ✓ or ✗ in the Actions tab.

## Usage

Mount the data directory at `/data` and run as your own user:

```bash
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes --help
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes aphot frame.fits
```

### Things that differ from upstream or older versions

- **`--scale` is in pixels per degree** in `artificial`, `astrometry` and `kombine`. This is an upstream change in hg1652/1653; older versions used degrees per pixel.
- **`artificial` has a detector-noise model.** The mean frame in electrons is
  flat · (sky + galaxies + stars) + dark · t,
  then Poisson noise, then Gaussian read noise, then division by gain, then the bias is added, then saturation is applied. The options:

  | Option | Units | Effect |
  |---|---|---|
  | `--dark` | e⁻/s/pixel | dark current; raises the level by dark·t, Poisson noise |
  | `--readnoise` | e⁻ rms | zero-mean Gaussian per pixel; independent of exposure time |
  | `--bias` | ADU | constant pedestal; use it with `--readnoise` and integer `-B`, so values are not clipped at 0 |
  | `--seed` | integer | reproducible frames, independent of `--mproc` |
  | `--gnoise` | ph/s/m²/arcsec² | **deprecated** extra background rate, kept for reproducibility of older datasets |

  The sky gradient (`--sky-grad-x/y`, mag/arcsec²/pixel) and the airmass act on the sky only.
  Random star fields (no `-c` catalogue) have a fixed density per square degree, so the number of stars follows the frame area (BUGS A23): 8 stars on the default 666² frame to the default `--maglim 13` (it was 239 on any frame).
  `--no-watermark` leaves out the mark of the artificial origin in the top-right corner, which `find` would detect as several stars; when it is drawn (the default), keywords WMARKX1/X2/Y1/Y2 give its region.

- **`artificial --hwhm` is σ for the GAUSS and SEEING PSFs**, and the scale parameter α for MOFFAT. This is the upstream meaning and is kept; the help text now says so.
- **The `HWHM` keyword written by `aphot` is now a true HWHM.** It used to be σ. Re-run `aphot` on frames processed by an older version before using the new `phcal` (see BUGS.md §5).
- **`phcal -q` quantities FLAM, FNU, PHOTLAM, PHOTNU, ABMAG and STMAG** are now physically correct. Previously FLAM was wrong by ×0.0013 and STMAG by +9.7 mag. MAG is unchanged.
- **The sky of `aphot`, `find` and `psfphot` is the mean of the sky**, clipped by the noise model of the frame (GAIN, READNS or RDNOISE, BIASLEV), not a robust mean or a median. On a sky of a few photo-electrons per pixel the old estimates were 13–30% low and faint stars too bright by up to +23% (BUGS P23). `aphot` leaves the other detected stars out of the sky ring and invalidates (−1) the apertures contaminated by a brighter neighbour (P29). On ground-based frames (hundreds of e⁻/pixel and more) SKY changes within its noise.
- **Exit status:** 0 on success, 1 on every error, with the reason on stderr (front-end errors gave 139, which looked like a crash, utility failures 13; BUGS A24).
- **`find` rejects detections explained by the wing of a brighter star** (the halo of an extended PSF, a saturated core; BUGS P25); its threshold now uses the clipped scatter of the sky, so on photon-starved frames it finds fewer noise peaks.
- **`find -th` is a false-alarm rate**: the threshold height is exceeded by the skewed (Poisson) noise of a faint sky as rarely as by Normal noise at *t*σ (by GAIN and BIASLEV; about 130 noise detections per 2048² at `-th 4`, 1 at `-th 5`, where a sky of 0.4 e⁻/px gave 505 and 21), and `--false-rate N` sets it by the expected number of noise detections per frame (BUGS P33). The sharpness and roundness tests accept the photon noise of faint stars (`--shape-sigma`, default 3): 100% instead of 76–84% of the stars at u = 15–17 on QUVIK-like frames (BUGS P32). FIND reports HTHRESH, SKEWFILT, FALSERAT, SHARPREF.
- **The errors follow the noise model of the detector** (BUGS P34–P39, S20): `aphot` applies the frame's GAIN to the Poisson noise of stars (APCOUNTERR of raw frames with GAIN 0.1 grows up to 3×; the gain used is APERRGN of APERPHOT, so pipelines that patched the errors and mark them by APERRGN do not patch them again) and the sky-pixel term uses Σf² of partial pixels; `psfphot` has no one-electron floor of the variance (χ² ≈ 1 on a sky of 0.4 e⁻/px) and `--residual-threshold` (default by the false-alarm rate, 5.0 on 2048²); the APEE errors of `--apcorr auto` account for the nesting of the apertures; `find` no longer adds READNS to the observed sky scatter; gphot's SKYERR no longer counts the ring error twice; `phcorr` STDERR includes the dark electrons, RDNOISE, no 1 e⁻ floor and writes GAIN = 1, READNS/RDNOISE in electrons and no BIASLEV for values in electrons; `kombine` STDERR is on the scale of the image (NaN for one value) and its default mean is clipped by the noise model (unbiased at low counts, rejects cosmic rays; the old one is `--hampel`). Validated by pulls ≈ 1 on simulated frames; see BUGS.md §5 for what changes in the outputs.
- **`phcal` pairs frames with filters by name** (BUGS P14): the catalogue columns (`--col-mag`, `--col-magerr`) and the `-C` values follow the order of `-f`, the frames are matched to them by their FILTER, so `-f V,B --col-mag Vmag,Bmag` equals `-f B,V --col-mag Bmag,Vmag` in any order of the files (before, an order other than that of the photometric system failed or paired a frame with another filter).
- **`phcal` with several filters has no colour term** (BUGS P43): it corrected the outputs by an "effective extinction" estimated from CTPH (0.1–0.2 mag per magnitude of B − V in B on simulated frames); every filter is now calibrated as alone (the stars must be on all frames), colour terms of an instrument are the transformation of `phfotran`.
- **`munipack phfotran` works** (BUGS P45; it was listed but not dispatched, and the utility crashed): it fits the linear transformation of the instrumental filters to the standard ones (photons to counts, a nearly tridiagonal matrix) on a calibration field, and `phcal --tratab` uses it for the calibration and its inverse for the photons: on simulated B,V,R frames of an instrument with filters shifted by mixing 4–13% of the neighbouring bands, the colour term −0.164 mag/mag in B drops to −0.0004 ± 0.0002 and the scatter to that of standard filters.
- **`phcal --color-term A,B[,errA,errB]`** (BUGS P47): the calibration of a filter whose band differs from that of the catalogue (QUVIK 260–360 nm against SDSS u or GALEX NUV) fits a colour term, m_cat − m_inst = ZP + k1 (c − c_ref) [+ k2 (c − c_ref)², `--color-order 2`], c = A − B of the catalogue, by the mixture likelihood of `--syserr auto` (outliers have no weight in k either; the errors of the colours are propagated). The term multiplies the photons, so MAG, FLUX and the other quantities stay consistent; CTPH is that of the reference colour (`--color-ref`, by default the mean colour of the calibration stars); every output star gets the correction of its catalogue colour, a star without one the mean correction of the calibration stars with its scatter as an error (COLOR, COLORERR, COLORSRC; COLTK1, COLTK1E, COLTREF, … in the header). On a simulated V frame calibrated against X = V + k (B − V): k = 0.0005, 0.1007, 0.2996 ± 0.0009 for k = 0, 0.1, 0.3, no residual colour slope, the scatter of the frame without a mismatch. With real SDSS colours (cones about RA 180, Dec 0 and M67; BUGS P55): a g frame calibrated against r (k = −1 exactly) gives −1.0003 ± 0.0023, a band between u and g (k = 0.4293) 0.4295 ± 0.0016, a QUVIK-like 310 nm band against u with u − g −0.686 ± 0.007 for the true −0.683, and QUVIK-like 20–60 s frames (find → aphot → astrometry → phcal) k within 2σ of the true value in 12 of 12 frames, Q 0.84–1.13. Without the option the outputs are bit-identical.
- **`phfotran --model mag`** (the default; BUGS P48): the standard transformation of every filter, m_std − m_inst = ZP + k1 (c − c_ref) [+ k2 (c − c_ref)²] with the colour of the filter and its neighbour (B − V for B and V, V − R for R), fitted as the colour term of phcal; `phcal --tratab` with such a table fits the zero points and solves the standard colours of every star from all its filters. An instrument with colour terms linear in magnitudes: rms 0.058/0.023/0.019 mag (B, V, R) without a table, 0.0136/0.0011/0.0026 with the flux model, 0.0012/0.0010/0.0008 with the magnitude model (as standard filters). `--model flux` is the photon-mixing matrix of P45, still better for a strong mixing of the bands (a leak); tables without TRAMODEL are flux tables.
- **`astrometry -p TPV` is fitted in stages** (BUGS S19): the gnomonic solution, a second-order deformation, all stars cross-matched again by it until the set is stable, then the third order iterated the same way; distortions of 15 px are recovered to 0.03–0.11 px (before: up to 6 px, 60 px with `--maxstars 1000`); undistorted frames give the same solution.
- **`artificial --sky-mag 99` gives no sky** (it added ≈1 e⁻/s/pixel) and integer frames are rounded (`-B 8` truncated; BUGS A33, A34).
- **`dark` averages by the mean clipped by the noise model of the frames** (GAIN, READNS or RDNOISE, BIASLEV or zero with `-bias`): the robust mean used before put a master dark of 0.13 e⁻/pixel 18% low (BUGS P31); now within ±0.0003 e⁻/pixel from 0.05 to 20 e⁻ with cosmic rays rejected, STDERR the error of the mean (n − 1). `--hampel` gives the old mean, `--arith` the plain one; `bias` keeps the robust mean (read noise is symmetric). One routine (`src/lib/clipmean.f08`) is the sky of aphot, find, psfphot, the default mean of kombine and the master of dark.

Example simulation of a 16-bit frame with a realistic detector:

```bash
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes \
  artificial -c cone.fits --rcen 180 --dcen 0 --width 4096 --height 4096 --scale 4762 \
  --exptime 20 --dark 6.7e-3 --readnoise 0.37 --gain 0.101 --bias 100 -B 16 --seed 1
```

## PSF photometry

`munipack psfphot` (new; it replaces the unused `pphot` stub) fits every star found by `find` with a PSF:

```bash
munipack find frame.fits
munipack aphot frame.fits                      # the sky annulus
munipack psfphot frame.fits                    # Gaussian, width fitted on bright isolated stars
munipack psfphot --psf MOFFAT frame.fits       # Moffat, FWHM and 1/beta fitted (beta > 1000: the Gaussian)
munipack psfphot --psf-file psf.fits frame.fits   # the field-dependent PSF of artificial --psf-file
munipack psfphot --psf EMPIRICAL --psf-var 2 --psf-save epsf.fits frame.fits   # empirical PSF, saved
munipack psfphot --psf EMPIRICAL --psf-var-radial 1 frame.fits   # a PSF changing with the field radius
munipack psfphot --psf EMPIRICAL --psf-file lab.fits frame.fits   # a lab PSF corrected by the frame
munipack psfphot --residual-detect 1 frame.fits   # add the stars missed by find
munipack phcal --photometry psf -c cone.fits ... frame.fits
```

- **One PSF code for simulation and analysis.** `src/lib/psfmodel.f08` loads, interpolates and pixelises PSFs for both `artificial` and `psfphot`, so frames simulated with a PSF file are fitted with exactly the same PSF.
- **PSF file:** primary `ZOOM`; every extension with `RADCOO` [deg] and `ANGCOO` [deg] is one PSF of the grid (any number of rings and angles). The PSF centre is the pixel NAXIS/2 or `CRPIX1/2` of the extension; the stamp should sum to ZOOM² (the whole flux). See `src/doc/dataform_artpsf.shtml`.
- **Field position** for a PSF file: pixel offsets from NAXIS/2 of the frame (as `artificial`), divided by the WCS scale [pix/deg]; `--psf-centre`, `--scale` override them.
- **Fit:** flux, position (and a constant sky with `--sky local`) by weighted least squares; the weights come from the model, max(model − BIASLEV, 0)/GAIN + read noise² (READNS [ADU] or RDNOISE [e⁻]), without a floor of one electron (BUGS P35). Neighbours are subtracted; three passes. The fixed sky is re-estimated in the `aphot` annulus on the frame with the fitted stars subtracted, leaving out the pixels where their light is significant (BUGS P42). The Moffat shape is fitted by the FWHM and u = 1/β, u = 0 being the Gaussian of the same FWHM, so Gaussian stars converge (BUGS P18; 454 → 10 s per 1024² frame); a median β above 1000 uses the Gaussian (PSFLIMIT = 'GAUSS').
- **Empirical PSF** (`--psf EMPIRICAL`, DAOPHOT-like): an analytic core (`--psf-core`, Moffat by default, shape and ellipticity fitted) or a PSF file (`--psf-file`, correction mode) plus residual tables sampled with ZOOM samples per pixel (`--psf-zoom`, 4 for FWHM < 5 px), of the radius of 99.5% of the flux by the stacked growth curve of the PSF stars (the largest over the field, at least 4 FWHM; `--psf-radius`), varying as a polynomial of order 0/1/2 over the field (`--psf-var`) and with the powers of the field radius from the PSF-file centre (`--psf-var-radial 0..3`; 1 by default for the correction of a field-dependent file with ≥ 20 PSF stars). The tables are solved by weighted linear least squares with the exact pixel-integration model on up to 100 automatically selected (or `--psf-list`) PSF stars, neighbours subtracted, their sky fitted on a ring outside the table, three iterations; no deposition onto the grid, hence no sub-pixel bias; the regularisation is relative to the local information of every sample and chosen by 4-fold cross-validation (`--psf-smooth` fixes it). The tables do not move the centroid; their flux is set by the **normalisation to a large aperture** (15–40 px) of the isolated PSF stars, applied when significant (3σ), so the light beyond the table is not taken from the wings of the analytic core (BUGS P50). **PSFMODERR** (a new column, included in PSFCOUNTERR) is the uncertainty of the model: the split-sample variance of the fluxes by the tables of 8 subsets of the PSF stars, the error of the normalisation and its systematic error (BUGS P52). In the correction mode the tables shrink towards zero unless significant (by the 8 subsets, per radial zone, thresholds by Student's t, with a 0.05%-of-flux floor), so a correct lab PSF is left unchanged. `--psf-save` writes the model in the PSF file format (one stamp, or rings × angles over the frame) for `psfphot --psf-file` and `artificial --psf-file`.
- **Groups:** stars closer than 1.5 FWHM (`--group`, 0 = off) are fitted simultaneously (fluxes and positions of all members, common local sky); groups above 15 stars are split by shorter links.
- **Residual detection** (`--residual-detect N`, default 0): find-like detection in the star-subtracted frame against the model noise plus the measured PSF mismatch, scaled by the observed local scatter, above `--residual-threshold` (by default find's threshold, at least the Normal deviate of the false-alarm probability 1/N_pix: 5.0 for 2048², RESTHRES; at 4σ the skewed noise of a 0.4 e⁻/px sky gives ~110 false detections per 2048² frame); added stars get FLAG 128 and all stars are refitted.
- **Output:** a `PSFPHOT` table: X, Y, XERR, YERR, PSFCOUNT, PSFCOUNTERR (total flux in ADU, like GCOUNT), SKY, SKYERR, CHI2, SHARP, ROUND, NPIX, FLAG (1 saturated core, 2 excluded pixels, 4 position fixed, 8 not converged, 16 failed, 32 crowded, 64 a saturated star not found by `find`, added so that it is subtracted from its neighbours, 128 found in the residuals, 256 S/N < 2: an upper-limit-like flux, BUGS P42), GROUP, GSIZE. `--residual` saves the star-subtracted frame.
- **Validation** on simulated frames (8 seeds, realistic noise): flux bias < 0.02% and position bias < 0.001 px for S/N > 50, reported errors 0.90–1.10 × the scatter, for Gauss and Moffat (shape fitted) and for PSF files (`dpsf.fits`, a 2-ring file). Stars with a neighbour of equal brightness within 2 FWHM: bias +0.07%, scatter 2.5% (aperture: +20%, 37%). End-to-end `phcal` MAG − V = 0.000 ± 0.0003, scatter 0.0017 (aperture 0.0020). A 2048² frame with 1700 stars takes 5–7 s.
- **Validation of the empirical PSF, groups and residual detection** (8 seeds): constant PSFs (Gauss σ 1.5 and 3, Moffat β 2.5, elliptical Moffat, 16-bit with saturated stars) — the empirical PSF agrees with the true one within 0.05% in flux (10 to 100 PSF stars), pull std 0.91–1.08, position bias < 0.002 px. A smooth field variation (focus + radial astigmatism): `--psf-var 2` +0.005% (true −0.005%), < 0.07% in every radius bin. Correction mode on `dpsf.fits` frames: the correct base gives the same fluxes as `--psf-file` (|Δ| < 0.001% on 16 crowded 16-bit seeds); a 10% broader / asymmetric base −0.03% / +0.13% (uncorrected: +5.1% / +0.9%). `--psf-save` round trip (99% of stars): ≤ 0.03% (constant PSF; Gauss ≤ 0.001%), ≤ 0.04% (order 2). Pairs at 1 FWHM: groups give +0.0% bias and pull std 1.0, iterative subtraction +26% / −19%; residual detection finds 88–100% of companions at 0.5–1 FWHM missed by `find`, with 0 spurious detections (7 runs × 3200 stars: float, 16-bit with saturated stars, Moffat; true, mismatched and empirical PSFs). Crowded field (2500 stars on 1024²): stars with a neighbour within 2 FWHM scatter 0.87% instead of 1.32%, those with an undetected companion within 1 FWHM +0.24% instead of +17.5% with one residual pass. End-to-end `phcal` with EMPIRICAL: MAG − V within ±0.0004, scatter 0.0016–0.0019. A 2048² frame with ~2000 stars: 6 s (EMPIRICAL), 9 s (order 2 + groups), 13.5 s (+ one residual pass).
- **Validation of the normalisation and the radial variation** (BUGS P50–P53; QUVIK-like 2048² u-band frames, 8 seeds, `dpsf.fits` — a 1.7 px core mixed linearly with the field radius with an 8.4 px halo; grids of 1764 stars u = 10–16 and the SDSS cone, bias of S/N > 50 stars in field-radius bins < 500 / 500–900 / > 900 px): `--psf-var-radial 1` −0.022% (+0.002 / −0.041 / −0.007%; `--psf-file` −0.009%), pull SD 0.99, PSFMODERR 0.02%; before, the analytic core with a 4 FWHM table gave −4.4% (+3.1 / −3.1 / −8.1%, pulls 9) and `--psf-var 2` −3.4% (pulls 3.1); now `--psf-var 2` −0.14% (±0.9% by radius, pulls 1.0 with PSFMODERR 1%) and the constant PSF −1.6% (±6%, pulls 1.04 with PSFMODERR 4.6%). Gauss and Moffat grids unchanged (−0.009%, −0.007%). Correction mode: the correct base gives the fluxes of the file (0.000% on 16 seeds, ≤ 0.017% on 4 crowded 16-bit seeds); a 10% broader base −0.01% (+0.008 / −0.030 / +0.007%; before +0.05%: −0.30 / −0.02 / +0.33%). End-to-end `phcal` (grid frames): MAG − u −0.0000 [−0.0006, +0.0002], scatter 0.0119 mag, the same as with the file (before: scatter 0.059 mag, −0.076 / +0.049 mag at the centre / edge). 2048² frames take 1.3–2× longer (3.7–8.7 s).
- **Limits:** the polynomial in x, y cannot follow a PSF that changes structure over the field (`dpsf.fits`: order 2 leaves ±0.8% in radial zones) — use `--psf-var-radial` when it changes with the field radius, or the correction mode with a model PSF as the base. The normalisation needs isolated bright stars (in a field of 16 px spacing: −0.27% with `--psf-var-radial 1`, pulls 1.0). PSF stars must be mostly single stars.

## Aperture corrections and calibration errors

A fixed aperture holds a fraction of a star's light, its enclosed energy (EE); with a PSF that changes over the field the aperture loss changes too, and one zero point cannot follow it. `aphot` optionally records EE of every star in every aperture (columns `APEEi`, `APEEERRi` of APERPHOT, appended after the old ones; keywords `EE*`), and `phcal` calibrates APCOUNT/APEE, the total flux, adding the relative EE error in quadrature:

```bash
munipack aphot --psf-file psf.fits frame.fits            # EE exactly from the PSF file (as psfphot/artificial)
munipack aphot --psf GAUSS --fwhm 4.9 frame.fits         # EE of an analytic PSF
munipack aphot --apertures 1,...,20,23,26.4,30.4,35,40 --ring 45,60 \
               --apcorr auto --apcorr-radial frame.fits   # EE from the frame's own growth curves
munipack phcal --saper 7 -c cone.fits ... frame.fits     # uses APEE when present (--apcorr none: ignore)
```

- **From a PSF:** the PSF at the star's field position is pixelised at its sub-pixel position (`psfmodel`) and summed with exactly the partial-pixel weights of the apertures (`daoapee`): on noiseless frames APCOUNT/APEE equals the true flux to 1e-6 in every aperture, also for Gauss and Moffat.
- **Automatic (`--apcorr auto`):** bright (S/N > 30), unsaturated stars without a neighbour brighter than 1% within 5.5 FWHM; total flux T, the mean of three plateau apertures: with the growth curve C = T·EE + b·A (b a residual sky), the first window of three apertures from 1.3 FWHM whose slope relative to T (per star, the covariance of the nested apertures, averaged with the excess scatter) is consistent with zero within 2σ, searched up to twice the radius enclosing 99% of the stacked growth curve (EEPLAT, EEPLATZ, EEPLATB, EER99, EERCAP; BUGS P40: the old 0.5% criterion let the plateau wander to 20–40 px on dense fields, T 0.1–0.4% too large; now the normalisation is within 0.016–0.052% rms of the exact EE over 8-frame sets); the slope of the growth curves beyond 4 FWHM recorded as a diagnostic of the sky (EESKYBIA; it was removed while the sky was biased, BUGS P23); EE of each aperture fitted by a polynomial of order `--apcorr-var` 0/1/2 (default 2) in x,y or, with `--apcorr-radial`, in the distance from the field centre; robust weighted fit (from the median) with the excess scatter added, the errors of C_k/T by the covariance of the nested apertures (BUGS P36), the order raised from the constant while the higher terms are significant (99.9%, by the weights of the lower order with its excess scatter) and limited for few stars (warning), no extrapolation. APEEERR = model error ⊕ star scatter.
- **Systematic error (`phcal --syserr auto|<mag>|0`, default auto):** the zero point (CTPH) and the excess scatter of the calibration stars over their errors are fitted together by the likelihood of a mixture of calibration stars and outliers (Hogg, Bovy & Lang 2010; EM from the median, the scatter by the restricted likelihood, unbiased for few stars; BUGS P44); outliers (misidentifications, blends, variables) have no weight in either. The scatter is added in quadrature to every output error (PHSYSERR, a COMMENT line); it is ≈ 0 on well-behaved frames. Monte Carlo: 30–100 stars within 1–6% with up to 10% outliers at 0.3–3 mag; with 10 stars and 10% outliers the zero point scatters by 6 mmag (179 mmag with the previous clipped estimator). With fewer than 16 stars (BUGS P54, P55; it was 8) the outlier model is fixed and the mode of the likelihood is chosen with a weak prior of σ, so σ² is unbiased (σ = 0.06 mag from 5 stars: rms 0.059, was 0.054) while gross outliers are still rejected. With 10 stars σ = 0.03 mag gives 0.028 on average (was 0.026). The 68% interval of σ (profile likelihood; an upper limit for σ = 0) is written as PHSYSERE (the upper error) with the number of stars PHSYSNST; fewer than 10 stars give a warning: for sparse fields (QUVIK UV frames with 2–8 stars) pass a fixed `--syserr` determined from many frames. With σ estimated, CTPHERR and the errors of the colour terms (and MAGERR by them) include its uncertainty (BUGS P55): the covariance of Knapp & Hartung (2003) with the floor of Röver et al. (2015), t² max(1, Q) times that of the weighted fit (Q the scatter of the stars in units of their errors, t the 68% quantile of Student's t of the degrees of freedom), so CTPHERR covers the zero point in 65–70% of sets of 3–100 stars for σ = 0.01–0.06 mag (54–63% before for 3–10 stars), conservatively 70–87% for σ = 0; 30 and more stars: errors 1–6% larger. `--syserr 0` or a value: the previous robust ratio, bit-identical. The zero point and the colour terms (`--color-term`, phfotran) share one implementation (`exscatter.f08`).
- **Validation** (QUVIK-like 4096² 16-bit u-band frames, 4 seeds each, aperture 7 = 5.1 px, stars < 17 mag): with the field-dependent PSF (`dpsf.fits`) the residuals versus field radius go from −0.12…−0.07 (centre) to +0.07…+0.08 mag (edge); with `--psf-file` every radial bin is within ±0.015 mag and the robust scatter drops from 0.073–0.077 to 0.041–0.047 mag (psfphot with the file: 0.042–0.053); `--apcorr auto --apcorr-radial` gives the same (±0.016 mag, 0.040–0.046). `--syserr auto` alone finds σ_sys = 0.058–0.062 mag and brings the error-calibration statistic Q (without the 2–3 misidentified detections per frame, BUGS P25) from 1.55–1.63 to 1.05–1.15. On Gaussian frames σ_sys = 0–0.003 mag and the automatic EE is constant, so results are unchanged.
- **Limits:** the automatic EE needs bright isolated stars (≈30 on these shallow 20 s frames; with so few, a 3σ fluctuation in the annulus of the brightest one can send the plateau to the cap: one of 8 sparse Gaussian frames, −0.14%, within its APEEERR) and a correct FIND FWHM: on a frame with a 9.9 px Gaussian PSF found with FWHM 5.2 (split detections), `--apcorr-radial` fitted a spurious variation (scatter 0.051 → 0.070 mag; the x,y model stayed constant). The additive sky bias of APCOUNT (BUGS P23) and the outliers that kept the unclipped Q of aperture routes high on the `dpsf.fits` frames (BUGS P25: split detections, faint stars in the wings of bright ones) are fixed.

## Growth-curve photometry

`munipack gphot` (rebuilt, BUGS P24, P49) measures the total flux GCOUNT of every star by the growth curve of the frame:

```bash
munipack aphot --apertures 1,...,20,23,26.4,30.4,35,40 --ring 45,60 frame.fits
munipack gphot frame.fits                  # the curve varying in x,y (--var 0/1/2, default 2)
munipack gphot --radial frame.fits         # varying with the distance from the field centre
munipack phcal --photometry grow -c cone.fits ... frame.fits
```

- **The growth curve** is the enclosed energy of the apertures by `aphot --apcorr auto` (bright isolated stars, the plateau by its slope, the variation over the field, the nested-aperture errors), its errors those of a prediction for the other stars: the model variance is inflated by the leave-one-out residuals of the curve stars (EELAMk) and the extrapolation outside their range adds its uncertainty.
- **GCOUNT** is the estimate of the least variance among the single apertures C_k/EE_k and the best linear estimates of the apertures 1..m by the covariance of nested apertures (sky, read noise, the photon noise by the model flux, the error of the mean sky), up to the plateau: faint stars by the profile-weighted apertures (the scatter of psfphot), bright stars by the large ones. GCOUNTERR includes the error of the curve (bounded as fully correlated between apertures). SKY is the sky of the ring corrected by the apertures (the prior 0 ± SKYERR), SKYERR its posterior error.
- **Validation** (2048² frames, 8 seeds; QUVIK-like u frames with a Gaussian PSF and `dpsf.fits`, sparse and dense fields; ground-like Gaussian and Moffat frames, gain 2, RN 5 e⁻, sky 560 e⁻/px): GCOUNT against the truth within ±0.3% at S/N > 30 and as psfphot's below (dense dpsf: −0.8…+0.1% by field radius, `--radial` −0.5%; old: +5.4…−3.5%), pull SD 0.9–1.15 by S/N, 1.24–1.33 for the brightest ground-based stars as psfphot's (robust 0.85–1.12; the old gphot: robust 1.40 for bright ground-based stars, 3.05 with Moffat wings); 4096² dpsf frames with `find` defaults: GCOUNT of all 900 stars (none before), `phcal --photometry grow` Q 1.01–1.14. `--legacy` runs the old growth-curve fit.
- **Limits:** at least 3 bright isolated stars (≈30 on the shallow QUVIK 2048² frames, where a PSF variation over the field is often not significant: the errors then describe it, the pull SD at S/N 10–30 is 1.14); FIND's FWHM sets the isolation. SKY of crowded fields holds the light of faint neighbours in the apertures (+0.3 SKYERR).

## Parallel processing

`artificial` and `astrometry` distribute the frames over `--mproc N` processes (Fortran coarrays on Open MPI by OpenCoarrays):

```bash
docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes \
  artificial -c cone.fits --count 8 --seed 3 --mproc 4
```

- One process (the default of the image, `MUNIPACK_CAF_MPROC=1`) runs the programs built for one process, exactly as before, without MPI. `--mproc N` > 1 runs their parallel variants (`artificial-caf`, `astrometry-caf`) by `mpirun`.
- The results do not depend on N: `artificial --seed` frames are bit-identical for 1, 2 and 4 processes (8 frames of 2048² with the SEEING PSF: 69 s with one process, 19.5 s with four), `astrometry` solutions are identical.
- The image sets `OMPI_*` variables for Open MPI in a container (see the Dockerfile): no UCX one-sided communication, `mpirun` as root, more processes than cores, no help banners.

## Known limitations

- The open issues are listed in BUGS.md.

## Licence

Munipack is free software under the GNU General Public License, version 3 or later (`src/COPYING`). OakLeaf is under the GNU Lesser General Public License (`deps/oakleaf/COPYING.LESSER`). The changes in this repository are distributed under the same licences.
