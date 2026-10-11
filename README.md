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
| `BUGS.md` | Bug catalogue, compatibility notes, suggested improvements |

The history starts with the untouched upstream tree (commit `003cb0d`), and each fix is a separate commit on top of it. `git diff 003cb0d -- src` therefore shows every change to Munipack, and `git log` gives the reason for each one.

## Build

```bash
docker build -t munipack:0.6.2-hg1653-fixes .
```

What the build does:
- Downloads and checks wxWidgets 3.3.3 and builds only its non-GUI part. This revision needs the wxWidgets 3.3 API; Ubuntu ships 3.2.
- Builds OakLeaf from `deps/`, then Munipack.
- Produces an image of about 170 MB.

The first build takes roughly 10 minutes; later builds reuse the cached dependency layers.

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

- **`artificial --hwhm` is σ for the GAUSS and SEEING PSFs**, and the scale parameter α for MOFFAT. This is the upstream meaning and is kept; the help text now says so.
- **The `HWHM` keyword written by `aphot` is now a true HWHM.** It used to be σ. Re-run `aphot` on frames processed by an older version before using the new `phcal` (see BUGS.md §5).
- **`phcal -q` quantities FLAM, FNU, PHOTLAM, PHOTNU, ABMAG and STMAG** are now physically correct. Previously FLAM was wrong by ×0.0013 and STMAG by +9.7 mag. MAG is unchanged.

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
munipack psfphot --psf MOFFAT frame.fits       # Moffat, width and beta fitted
munipack psfphot --psf-file psf.fits frame.fits   # the field-dependent PSF of artificial --psf-file
munipack phcal --photometry psf -c cone.fits ... frame.fits
```

- **One PSF code for simulation and analysis.** `src/lib/psfmodel.f08` loads, interpolates and pixelises PSFs for both `artificial` and `psfphot`, so frames simulated with a PSF file are fitted with exactly the same PSF.
- **PSF file:** primary `ZOOM`; every extension with `RADCOO` [deg] and `ANGCOO` [deg] is one PSF of the grid (any number of rings and angles). The PSF centre is the pixel NAXIS/2 or `CRPIX1/2` of the extension; the stamp should sum to ZOOM² (the whole flux). See `src/doc/dataform_artpsf.shtml`.
- **Field position** for a PSF file: pixel offsets from NAXIS/2 of the frame (as `artificial`), divided by the WCS scale [pix/deg]; `--psf-centre`, `--scale` override them.
- **Fit:** flux, position (and a constant sky with `--sky local`) by weighted least squares; the weights come from the model, (model − BIASLEV)/GAIN + read noise² (READNS [ADU] or RDNOISE [e⁻]). Neighbours are subtracted; three passes. The fixed sky is re-estimated in the `aphot` annulus on the frame with the fitted stars subtracted.
- **Output:** a `PSFPHOT` table: X, Y, XERR, YERR, PSFCOUNT, PSFCOUNTERR (total flux in ADU, like GCOUNT), SKY, SKYERR, CHI2, SHARP, ROUND, NPIX, FLAG (1 saturated core, 2 excluded pixels, 4 position fixed, 8 not converged, 16 failed, 32 crowded, 64 a saturated star not found by `find`, added so that it is subtracted from its neighbours). `--residual` saves the star-subtracted frame.
- **Validation** on simulated frames (8 seeds, realistic noise): flux bias < 0.02% and position bias < 0.001 px for S/N > 50, reported errors 0.90–1.10 × the scatter, for Gauss and Moffat (shape fitted) and for PSF files (`dpsf.fits`, a 2-ring file). Stars with a neighbour of equal brightness within 2 FWHM: bias +0.07%, scatter 2.5% (aperture: +20%, 37%). End-to-end `phcal` MAG − V = 0.000 ± 0.0003, scatter 0.0017 (aperture 0.0020). A 2048² frame with 1700 stars takes 5–7 s.

## Known limitations

- `--mproc` (coarray multiprocessing) is not available: the image is built with `-fcoarray=single`.
- The remaining open issues are listed in BUGS.md.

## Licence

Munipack is free software under the GNU General Public License, version 3 or later (`src/COPYING`). OakLeaf is under the GNU Lesser General Public License (`deps/oakleaf/COPYING.LESSER`). The changes in this repository are distributed under the same licences.
