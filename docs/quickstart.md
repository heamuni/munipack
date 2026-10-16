# Quick start

This page goes from nothing to calibrated magnitudes on a single frame. Because the frame is simulated, the true magnitudes are known and you can check the result. The steps are:

1. get a reference catalogue (`cone`, or your own FITS table);
2. simulate a frame with realistic detector noise (`artificial`);
3. detect stars (`find`);
4. measure aperture photometry with an aperture correction (`aphot --apcorr auto`);
5. calibrate the astrometry (`astrometry`);
6. calibrate the photometry against the catalogue (`phcal`);
7. optionally, use PSF photometry instead of apertures (`psfphot`).

The commands of steps 1–7 were run as shown. The numbers quoted come from one such run, and other builds can differ in the last digits. The field is the open cluster M67 (RA 132.85°, Dec +11.81°), seen in V on a 1024 × 1024 detector at 1″ per pixel.

## Setup

Run Munipack from the container, with the current directory mounted at `/data` and running as your own user. A shell function lets you write `munipack ...` exactly as you would with a native install:

```bash
munipack() {
  docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/data \
    munipack:0.6.2-hg1653-fixes "$@"
}
munipack --help
```

Every action works on files in the current directory. Most of them **add extensions to the input frame** instead of writing new files, so one FITS file collects the whole history of its reduction.

```{note}
Some runs print `Note: The following floating-point exceptions are signalling: IEEE_DENORMAL`. This is a notice from the Fortran runtime, not an error. The exit status is 0 on success and 1 on any error.
```

## 1. Reference catalogue

### From a Virtual Observatory cone search

```bash
munipack cone -r 0.2 -o cone.fits -- 132.85 11.81
```

**Output:** `cone.fits`, a FITS table of UCAC4 (the default catalogue) within 0.2° of the centre. The extension is named after the VizieR identifier, `I/322A`. It has the columns `RAJ2000`, `DEJ2000`, `pmRA`, `pmDE` and the APASS magnitudes `Bmag`, `Vmag`, `gmag`, `rmag`, `imag`, each with its error `e_*`.

- The radius `-r` is in degrees (default 0.1°). Choose it to cover the frame: this frame spans 0.28° and its half-diagonal is 0.2°.
- `-c` selects another catalogue (`munipack cone --list-catalogues` lists them, e.g. `UCAC5`, `GaiaDR3`, `APASS`), and `--id` takes any VizieR catalogue, e.g. `--id J/A+A/544/A31`.
- Put `--` before the coordinates. Otherwise a negative declination is read as an option: `munipack cone -r 0.2 -- 245.9 -26.5`.

`artificial`, `astrometry` and `phcal` recognise a catalogue that `cone` downloaded, so its coordinate and proper-motion columns need no `--col-*` options.

### Or a local FITS table

You can use any FITS binary table with the coordinates in degrees and a magnitude column, for example one exported from astroquery or TOPCAT. Name its columns explicitly:

```bash
C="--col-ra RA --col-dec DEC --col-pm-ra pmRA --col-pm-dec pmDE --col-mag V"
munipack artificial -c stars.fits $C -f V --rcen 132.85 --dcen 11.81 ...   # other options as in step 2
munipack astrometry -c stars.fits $C frame.fits
munipack phcal -c stars.fits --col-ra RA --col-dec DEC --col-mag V --col-magerr e_V \
  --photsys-ref Johnson -f V frame.fits
```

```{note}
For a table it does not recognise, `artificial` needs proper-motion columns, given by `--col-pm-ra` and `--col-pm-dec`. If the motions are unknown, add zero columns. `astrometry` and `phcal` work without them. `astrometry` applies proper motions only when the epoch of the catalogue is known, from the `JEPOCH` header keyword or the `--col-epoch-ra`, `--col-epoch-dec` columns. Without it, the motions are ignored and a warning is printed.
```

The rest of this page uses `cone.fits`.

## 2. Simulate a frame

```bash
munipack artificial -c cone.fits -f V --rcen 132.85 --dcen 11.81 \
  --scale 3600 --width 1024 --height 1024 --exptime 10 \
  --psf GAUSS --hwhm 1.5 \
  --gain 2 --readnoise 5 --bias 100 --dark 0.01 -B 16 \
  --seed 1 --no-watermark --mask frame.fits
```

**Output:** `frame.fits`. The primary HDU is a 16-bit image in ADU with a WCS header (TAN projection) and the keywords of the noise model (`GAIN`, `RDNOISE`, `BIASLEV`, `SATURATE`, `FILTER`, `EXPTIME`). Two extensions follow: `PSF`, the PSF that was used, and `CATALOGUE`, which holds the true position (`X`, `Y`) and the mean number of electrons (`COUNT`) of every star rendered.

| Option | Meaning |
|---|---|
| `-f V` | filter. The magnitudes are read from the column `Vmag` (the filter name followed by `mag`) unless `--col-mag` says otherwise |
| `--scale 3600` | pixels per degree, i.e. 1″ per pixel |
| `--hwhm 1.5` | for `GAUSS` (and `SEEING`) this is the Gaussian **σ** in pixels, not the half width. Here FWHM = 2.355 σ = 3.5 px |
| `--gain`, `--readnoise`, `--bias`, `--dark` | e⁻/ADU, e⁻ rms, ADU, e⁻/s/pixel |
| `-B 16` | 16-bit integer output. The bias keeps the read noise from being clipped at zero |
| `--seed 1` | reproducible noise |
| `--no-watermark` | leaves out the mark in the top-right corner, which `find` would detect as stars |

The default telescope is a 1.12 m aperture (`--area 0.985` m²) under a V sky of 21 mag/arcsec² (`--sky-mag`).

## 3. Detect stars

```bash
munipack find -f 3.5 frame.fits
```

**Output:** a `FIND` table extension in `frame.fits` (about 300 stars here), with the positions `X`, `Y`, the peak height relative to the threshold `PEAKRATIO`, and the shape indices `SHARP` and `ROUND`. The header records `FWHM` and the thresholds (`THRESH`, `HTHRESH`).

- `-f` is the FWHM in pixels. Give a value close to the real one.
- `-th` (default 7) is the threshold as a false-alarm rate: a noise peak exceeds it as rarely as Normal noise exceeds *t*σ, even on a faint Poisson sky. `--false-rate N` instead sets the expected number of noise detections per frame.

## 4. Aperture photometry

```bash
munipack aphot --apcorr auto frame.fits
```

**Output:** an `APERPHOT` table extension in `frame.fits`, with `APCOUNTi` and `APCOUNTERRi` (counts in ADU and their errors) for every aperture *i*, the sky (`SKY`, `SKYERR`) and, because of `--apcorr auto`, the enclosed energy of each aperture, `APEEi` and `APEEERRi`. The header holds the aperture radii (`APER1`…`APER12`), a true `HWHM`, and the `EE*` keywords of the growth-curve model (`EENSTAR` stars used, the plateau radii `EERTOT1`, `EERTOT2`).

- By default there are 12 apertures from 1 to 20 px (DAOPHOT-like spacing). Set them with `--apertures 1,1.5,2,...` and the sky ring with `--ring ri,ro`.
- `--apcorr auto` measures the growth curves of bright isolated stars and fits their enclosed energy over the field (`--apcorr-var 0|1|2`, or `--apcorr-radial` for a variation with the distance from the field centre). If the PSF is known, use `--psf-file psf.fits` or `--psf GAUSS|MOFFAT --fwhm ...` instead.
- The errors follow the noise model in the header (gain, read noise, bias).

## 5. Astrometry

### With the WCS in the header

`artificial` writes a WCS, and real frames often carry an approximate one from the telescope:

```bash
munipack astrometry -c cone.fits frame.fits
```

**Output:** the fitted WCS (`CRVAL*`, `CRPIX*`, `CD*_*`) is written into the primary header of `frame.fits`. `ASTRMETH = 'WCS'` records how the match was found, and `COMMENT` lines report the number of stars and the RMS (0.12″ on 141 stars here). No new extension is added.

The header WCS guides the match within a tolerance: `--wcs-offset` (default 10% of the frame, at least 100 px; units are accepted, e.g. `30arcsec`, `2arcmin`, `0.1deg`, `20%`), `--wcs-rot` (default 3°) and `--wcs-scale` (default 3%). `--blind` ignores the header.

### From the pointing, without a WCS

If the frame has no WCS, give the pointing and the scale. The roll of the field may be unknown:

```bash
munipack astrometry --remove frame.fits       # only to demonstrate: drop the WCS
munipack astrometry -c cone.fits --rcen 132.86 --dcen 11.80 --scale 3600 frame.fits
```

The match tries any rotation (`ASTRMETH = 'POINTING'`) and gives the same solution as above. Useful options:

- `--angle` for an approximately known rotation. It is searched within `--wcs-rot` first, then any rotation is tried. Add `--roll-known` to skip the search over any rotation, or use `--roll-unknown` to search any rotation from the start;
- `--reflex-unknown` for a possibly mirrored frame;
- `--scale-tol` for an uncertain scale (default 0.05, i.e. 5%; `--scale-free` accepts any scale).

The catalogue may be larger than the frame: the brightest `--maxstars` stars near the pointing are used.

## 6. Photometric calibration

```bash
munipack phcal -c cone.fits --photsys-ref Johnson -f V \
  --col-mag Vmag --col-magerr e_Vmag --syserr auto frame.fits
```

**Output:** a **new file** `frame_cal.fits`. Its primary HDU is the image converted to photons (`BUNIT = 'photon'`). Its `PHOTOMETRY` table has, for every star, `RAJ2000`, `DEJ2000`, `SKY`, `SKYERR`, `PHOTON`, `FLUX`, `MAG` and the error of each. Stars that cannot be measured (saturated, contaminated) have `MAG = 99.99999`. The header holds the zero point `CTPH ± CTPHERR` (photons per count), the systematic error `PHSYSERR` with its upper error `PHSYSERE` and the number of stars `PHSYSNST`, and `APCORTYP = 'EE'`, which says the counts were divided by `APEE`.

- `-f` and `--photsys-ref` name the filter and the photometric system. The frame's `FILTER` keyword must be one of the filters in `-f`.
- `--apcorr auto` is the default: when APEE columns are present, the total flux APCOUNT/APEE is calibrated. Use `--apcorr none` to calibrate the plain aperture counts, or give a constant factor.
- `--syserr auto` (the default) fits the zero point together with the scatter of the calibration stars in excess of their errors, and rejects outliers such as blends or variables. The scatter is added in quadrature to every error. It is 0 here because the frame was rendered from `Vmag` exactly. For sparse fields (fewer than about 10 stars), pass a fixed value determined from many frames, e.g. `--syserr 0.02`.
- `--saper N` chooses the aperture. `-q` chooses the output quantities (default `FLUX,MAG`; also `PHOTON`, `FLAM`, `FNU`, `ABMAG`, `STMAG`, ...).

Here MAG − V has a median of +0.001 mag and a robust scatter of 0.003 mag on 285 stars.

### Optional: a colour term

If the instrument's band differs from the catalogue band, fit a colour term from the catalogue colour A − B (with the errors of both columns):

```bash
munipack phcal -c cone.fits --photsys-ref Johnson -f V \
  --col-mag Vmag --col-magerr e_Vmag \
  --color-term Bmag,Vmag,e_Bmag,e_Vmag -O --mask '\1_ct.\2' frame.fits
```

**Output:** `frame_ct.fits`. `-O --mask '\1_ct.\2'` sets the output name: `\1` is the base name and `\2` the extension. Its `PHOTOMETRY` table adds the columns `COLOR`, `COLORERR` and `COLORSRC`, and the header records the term: `COLTERM = 'Bmag-Vmag'`, `COLTK1 ± COLTK1E` (mag/mag) and the reference colour `COLTREF`. The fit is m_cat − m_inst = ZP + k₁ (c − c_ref). This simulated band is exactly V, so k₁ comes out small, a few hundredths of a mag/mag at most. A real colour term of an instrument is usually larger.

- `--color-order 2` adds a quadratic term.
- `--color-ref` fixes the reference colour (default: the mean of the calibration stars).
- The option can be repeated, up to three colours (e.g. u − g and g − r).

## 7. PSF photometry (alternative to apertures)

`psfphot` fits every star detected by `find`. It uses the sky measured by `aphot`, so run it after steps 3–4:

```bash
munipack psfphot frame.fits
munipack phcal -c cone.fits --photsys-ref Johnson -f V \
  --col-mag Vmag --col-magerr e_Vmag --photometry psf -O --mask '\1_psf.\2' frame.fits
```

**Output:** `psfphot` adds a `PSFPHOT` table extension to `frame.fits`. It has `X`, `Y`, `PSFCOUNT`, `PSFCOUNTERR` (total flux in ADU), `SKY`, `CHI2`, `FLAG`, `GROUP`, and more. `phcal --photometry psf` then writes `frame_psf.fits` with `PHOTOTYP = 'PSFPHOT'` in its `PHOTOMETRY` table. Here MAG − V has a median of +0.001 mag and a robust scatter of 0.0025 mag on 320 stars.

The default PSF is a Gaussian whose width is fitted. Other choices:

```bash
munipack psfphot --psf MOFFAT frame.fits                     # Moffat, FWHM and beta fitted
munipack psfphot --psf-file psf.fits frame.fits              # a known PSF file (as for artificial --psf-file)
munipack psfphot --psf EMPIRICAL --psf-var 2 frame.fits      # empirical PSF varying over the field
munipack psfphot --residual-detect 1 frame.fits              # also fit stars missed by find
```

## Inspecting the results

```bash
munipack fits -K CTPH,CTPHERR,PHSYSERR 'frame_cal.fits[PHOTOMETRY]'   # header keywords
munipack fits -lh 'frame.fits[APERPHOT]'                              # a whole header
munipack fits -lt 'frame_cal.fits[PHOTOMETRY]'                        # table rows as text
```

Any FITS reader works as well, e.g. `astropy.io.fits` or TOPCAT.

## Files and extensions

| Step | Command | Writes |
|---|---|---|
| catalogue | `cone` | `cone.fits` (table) |
| simulation | `artificial` | `frame.fits`: image + `PSF`, `CATALOGUE` |
| detection | `find` | `FIND` extension in the frame |
| aperture photometry | `aphot` | `APERPHOT` extension (with `APEEi` for `--apcorr auto`) |
| astrometry | `astrometry` | WCS in the primary header of the frame |
| calibration | `phcal` | new file `frame_cal.fits`: image in photons + `PHOTOMETRY` |
| PSF photometry | `psfphot` | `PSFPHOT` extension in the frame |

Re-running an action replaces its extension. `phcal` overwrites its output file.

## Next steps

- [Concepts](concepts.md): files and extensions, units, the noise model and the errors, photometric systems, the PSF file format.
- [Accuracy and validation](validation.md): what accuracy to expect from each method, measured on simulated frames.
- The reference pages list every option: [artificial](reference/artificial.md), [find](reference/find.md), [aphot](reference/aphot.md), [psfphot](reference/psfphot.md), [gphot](reference/gphot.md), [astrometry](reference/astrometry.md), [phcal](reference/phcal.md), [catalogues](reference/catalogues.md), [preprocessing](reference/preprocessing.md) (bias, dark, flat), [kombine](reference/kombine.md), [utilities](reference/utilities.md).
