# phcal

`phcal` converts the instrumental counts of a frame into photons, fluxes and
magnitudes, calibrated against a reference catalogue or an already calibrated
frame.

It matches the stars of the frame with the reference, fits the calibration
ratio **CTPH** (photons per count) and rewrites the image and the photometry
table in calibrated units. To correct for the difference between your filters
and a standard system, measure a transformation table with {doc}`phfotran`
and pass it with `--tratab`.

## Usage

```bash
munipack phcal [options] file(s)[,result(s)]
```

Each frame needs an astrometric solution ({doc}`astrometry`), instrumental
photometry ({doc}`aphot`, {doc}`gphot` or {doc}`psfphot`), the date of
observation and a `FILTER` keyword. The exposure time and telescope area are
read from `EXPTIME` and `AREA`; if they are missing, `phcal` warns and uses
1 s and 1 m², so give `--area` when your headers lack it.

By default `frame.fits` is written as `frame_cal.fits`. Use `-t` to write
into another directory, or give input and output as one quoted argument
separated by a space: `'frame.fits result.fits'`. The comma form shown by
`--help` does not work.

## Options

### Calibration source

Give exactly one of `-c`, `-r` or `-C`, always together with `--photsys-ref`.

| Option | Argument | Default | Description |
|---|---|---|---|
| `-c`, `--cat` | file | – | Reference catalogue (FITS table, e.g. from `munipack cone`). Needs `-f`. |
| `-r`, `--ref` | file | – | A frame already calibrated by `phcal`, used as the reference (relative photometry, no colour terms). |
| `-C`, `--cal` | value(s) | – | Set CTPH by hand, one value per filter in the order of `-f`. |

### Catalogue and filters

| Option | Argument | Default | Description |
|---|---|---|---|
| `-f`, `--filters` | list | – | Filters to calibrate, comma separated; frames are paired by their `FILTER` keyword. |
| `--col-ra` | name | `RAJ2000` | Right Ascension column. |
| `--col-dec` | name | `DEJ2000` | Declination column. |
| `--col-mag` | list | – | Magnitude column(s), in the order of `-f`. |
| `--col-magerr` | list | none | Magnitude error column(s); without them errors come from photon counts. |
| `--vocat` | file | bundled | VO catalogue configuration (column IDs of known catalogues). |

### Star selection and matching

| Option | Argument | Default | Description |
|---|---|---|---|
| `-th`, `--threshold` | number | 0.1 | Use only stars brighter than this multiple of the sky in the aperture. |
| `-e`, `--maxerr` | number | 0.1 | Use only stars with a relative flux error below this. |
| `--tol` | deg | HWHM of the frame | Cross-match radius. |

### Which photometry

| Option | Argument | Default | Description |
|---|---|---|---|
| `--photometry` | `aper`, `grow`, `psf` | `grow` if present, else `aper` | Photometry to calibrate: `aphot`, `gphot` or `psfphot`. |
| `--saper` | integer | 0 (automatic) | Aperture number to use from the aperture photometry. |
| `--apcorr` | `auto`, `none`, factor | `auto` | Aperture correction: enclosed energy when available, none, or a constant factor. |

### Errors and colour terms

| Option | Argument | Default | Description |
|---|---|---|---|
| `--syserr` | `auto` or mag | `auto` | Systematic error added to all errors: fitted, or a fixed value (`0` turns it off). |
| `--color-term` | `A,B[,errA,errB]` | none | Fit a colour term in the catalogue colour A − B; repeat for up to 3 colours. |
| `--color-order` | `1`, `2`, or list | `1` | Linear or quadratic colour term (a list gives one per `--color-term`). |
| `--color-ref` | `mean`, mag, or list | `mean` | Colour at which CTPH is defined; fix it (e.g. `0`) to compare CTPH between fields. |
| `--color-err-correct` | `yes`, `no` | `yes` | Correct the colour term for the errors of the catalogue colours. |

### Photometric systems and output units

| Option | Argument | Default | Description |
|---|---|---|---|
| `--photsys-ref` | name | – | Photometric system of the catalogue (required); bundled are `Johnson`, `Gunn` and `Stromgren`. |
| `--photsys-instr` | name | `PHOTSYS` of the frames | Instrumental photometric system. |
| `--area` | m² | `AREA` of the frames | Area of the telescope's entrance aperture. |
| `-q`, `--quantity` | list | `FLUX,MAG` | Output quantities (see {ref}`output-quantities`). |
| `--tratab` | file | none | Transformation table made by `phfotran`. |
| `--phsystab` | file | bundled | Table of photometric system definitions. |
| `--list` | – | – | List the available photometric systems and exit. |

### Output and general

| Option | Argument | Default | Description |
|---|---|---|---|
| `--advanced` | – | off | Also keep the input extensions and add `PHRES`, the residuals of the calibration stars. |
| `-t`, `--target-directory` | dir | – | Write outputs into this directory. |
| `-O` | – | off | Turn on advanced output names (`--pattern`, `--mask`, `--format`). |
| `--pattern` | regex | – | Pattern for input names (with `-O`). |
| `--mask` | mask | – | Output name mask, e.g. `'\1_MAG.\2'` (with `-O`). |
| `--format` | format | – | Output file format (with `-O`). |
| `--verbose` | – | off | Print the fit results. |
| `-h`, `--help` | – | – | Show the help. |

## How it works

`phcal` converts each catalogue magnitude into an expected number of
photons, using the passband of the photometric system, the exposure time and
the telescope area. CTPH is the ratio of these photons to the measured
counts; for counts in ADU it also contains the gain. The fit is robust: misidentified, blended or variable stars are
recognised as outliers and get no weight, so you do not need to clean the
catalogue by hand.

**Calibrate the total flux.** An aperture holds only part of a star's light,
and that part changes across the field when the PSF does. To calibrate the
total flux, use `gphot` or `psfphot` photometry, or run `aphot` with an
enclosed-energy model and keep the default `--apcorr auto`.

**Systematic error.** By default `phcal` measures the extra scatter of the
calibration stars beyond their quoted errors and adds it to every output
error. With very few stars (it warns below 10) this estimate is little more
than a guess; then measure it from many frames of the same setup and pass a
fixed value, e.g. `--syserr 0.02`.

**Colour terms.** Your filter is never exactly the catalogue's, so the
calibration depends slightly on the colour of the star. `--color-term`
fits this dependence in a catalogue colour. Stars without a catalogue colour
get the average correction, with a correspondingly larger error. If several
colours are given they must be independent (not u − g, g − r and u − r).
Colour terms need a catalogue (`-c`) and cannot be combined with `-r` or
`--tratab`.

**Noisy catalogue colours.** Errors in the catalogue colours flatten the
fitted colour term. When the colour columns come with error columns
(`--color-term u,g,e_u,e_g`), `phcal` corrects for this. With fewer than
about 17 stars the correction is unreliable, so it is skipped with a
warning. If your catalogue colours are more precise than their quoted
errors (e.g. frames simulated from the catalogue), use
`--color-err-correct no`, otherwise the colour term is over-corrected.

**Several filters.** Each filter is calibrated on its own, but only stars
present on all frames are used and written out. Colour corrections between
filters come only from `--color-term` or `--tratab`.

(output-quantities)=

## Output quantities

All quantities follow from the photon count n = CTPH × counts. Undefined
magnitudes are written as 99.99999.

| Quantity | Meaning | Unit |
|---|---|---|
| `PHOTON` | photon count | photon |
| `PHOTRATE` | photon flux | photon/s/m² |
| `PHOTNU` | photon flux per frequency | photon/s/m²/Hz |
| `PHOTLAM` | photon flux per wavelength | photon/s/m²/nm |
| `FLUX` | energy flux in the band | eV/s/m² |
| `FNU` | flux density per frequency | W/m²/Hz |
| `FLAM` | flux density per wavelength | W/m²/nm |
| `MAG` | magnitude in the photometric system | mag |
| `ABMAG` | AB magnitude | mag |
| `STMAG` | ST magnitude | mag |

## Output

One calibrated file per frame; existing output files are overwritten. The
image is converted to photons. The `PHOTOMETRY` table gets the photon counts
and a value/error column pair for each `-q` quantity (e.g. `MAG`, `MAGERR`).
Its header holds `CTPH`, `CTPHERR`, the systematic error `PHSYSERR`, any
colour-term coefficients, and a list of the calibration stars.

## Examples

Calibrate a V frame against a catalogue, writing magnitudes, flux densities
and AB magnitudes:

```bash
munipack cone -c UCAC4 -o cat.fits -r 0.2 -- 110.47 71.34
munipack phcal -c cat.fits --photsys-ref Johnson -f V --col-mag Vmag --col-magerr e_Vmag \
    -q MAG,FLAM,ABMAG frame_V.fits
```

Calibrate PSF photometry of B and V frames, fitting a colour term in B − V:

```bash
munipack phcal --photometry psf -c cat.fits --photsys-ref Johnson --area 0.3 \
    -f V,B --col-mag Vmag,Bmag --col-magerr e_Vmag,e_Bmag \
    --color-term Bmag,Vmag,e_Bmag,e_Vmag --verbose -t calibrated frame_B.fits frame_V.fits
```

Calibrate a later frame relative to one already calibrated:

```bash
munipack phcal --photsys-ref Johnson -r frame_001_cal.fits frame_002.fits
```
