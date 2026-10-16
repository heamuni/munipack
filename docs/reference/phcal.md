# phcal, phfotran

Photometric calibration of frames against a reference catalogue (`phcal`), and the transformation of an instrumental set of filters to a standard photometric system (`phfotran`).

```{note}
The examples call `munipack` directly. With the container, put
`docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes`
in place of `munipack`, and keep the files in the current directory.
```

## phcal

`phcal` converts the instrumental counts of a frame into physical quantities. It cross-matches the stars of the frame with a reference catalogue, or with a previously calibrated frame, by their sky coordinates. It converts the catalogue magnitudes into expected photon rates, using the passband of the photometric system and the telescope area, and fits the calibration ratio **CTPH**, the number of photons per count (despite the header comment "counts per photons"). The fit uses a mixture likelihood, so outliers carry no weight. Optionally it also fits a systematic error, one or several colour terms, or applies a transformation table from `phfotran`. The image and the photometry table are then rewritten in photons and in the quantities chosen by `-q` (magnitudes, energy fluxes, flux densities, …). Every frame needs an astrometric solution (`astrometry`), instrumental photometry (`aphot`, `gphot` or `psfphot`), the date of observation and a `FILTER` keyword. The exposure time and the area are read from the header (`EXPTIME`, `AREA`); when they are missing, `phcal` warns and uses 1 s and 1 m².

### Synopsis

```bash
munipack phcal [options] file(s)[,result(s)]
```

Without `-O` or `-t`, `frame.fits` is written as `frame_cal.fits`. To name the output, give the input and output names as one argument separated by a space, `'frame.fits result.fits'` (quoted). The comma form `file,result` shown by `--help` is not split and fails.

### Options

#### Calibration source

You need exactly one of these three sources. Every mode also needs `--photsys-ref`, and `-c` and `-C` need `-f`.

| Option | Argument | Default | Description |
|---|---|---|---|
| `-c`, `--cat` | file | – | Reference catalogue: a FITS table with RA, Dec, magnitudes (and optionally their errors). It can come from `munipack cone` or be built by hand. |
| `-r`, `--ref` | file | – | Reference frame: a frame already calibrated by `phcal`. Its stars become the reference (relative photometry; no colour transformation). Cannot be combined with `--color-term`. |
| `-C`, `--cal` | value(s) | – | Manual calibration: the CTPH value(s), one per filter, in the order of `-f`. |

#### Catalogue columns and filters

| Option | Argument | Default | Description |
|---|---|---|---|
| `-f`, `--filters` | list | – | Filters to calibrate, comma separated, one per frame (required with `-c` and `-C`). Frames are paired with filters by their `FILTER` keyword, so neither the order of `-f` nor the order of the files matters. |
| `--col-ra` | name | `RAJ2000` | Right Ascension column of the catalogue (for a catalogue known to `--vocat`, its RA column). |
| `--col-dec` | name | `DEJ2000` | Declination column of the catalogue (for a catalogue known to `--vocat`, its Dec column). |
| `--col-mag` | list | – | Magnitude column(s), in the order of `-f` (`-f V,B --col-mag Vmag,Bmag`). |
| `--col-magerr` | list | none | Magnitude error column(s), in the order of `-f`. Without them, the errors are estimated from the photon counts. |
| `--vocat` | file | bundled `VOcat_conf.xml` | VO catalogue configuration (column IDs of known catalogues). |

#### Star selection and matching

| Option | Argument | Default | Description |
|---|---|---|---|
| `-th`, `--threshold` | number | 0.1 | Keep stars whose counts exceed this multiple of the sky counts in the same aperture. Applied to the stars of the calibrated frames, not to the catalogue or to a reference frame (`-r`). |
| `-e`, `--maxerr` | number | 0.1 | Keep stars whose relative flux error is below this value. 0.01 selects only precise stars. Applied as `-th`. |
| `--tol` | deg | HWHM of the frame | Cross-match radius in degrees. By default it is the `HWHM` keyword of the photometry converted to degrees by the WCS scale. |

#### Photometry used and aperture correction

| Option | Argument | Default | Description |
|---|---|---|---|
| `--photometry` | `aper`, `grow`, `psf` | `grow` if GROWPHOT exists, else `aper` | `aper`: APERPHOT (`aphot`). `grow`: GROWPHOT, the total flux GCOUNT (`gphot`). `psf`: PSFPHOT, PSFCOUNT (`psfphot`); stars with a saturated core or a failed fit are not used for the calibration. |
| `--saper` | integer | 0 (automatic) | Aperture number in APERPHOT. 0 picks the first aperture larger than 2σ of a Gaussian of the frame's HWHM. A positive value selects APERPHOT even when GROWPHOT exists, unless `--photometry grow` is given. |
| `--apcorr` | `auto`, `none`, factor | `auto` | `auto`: divide the aperture counts by their enclosed energy `APEEi` when APERPHOT has it (`aphot --psf-file`, `--psf`, `--apcorr`), so the total flux is calibrated, and add the relative APEE error in quadrature. Frames without APEE are calibrated by the aperture flux. `none`: ignore APEE. A number: a constant total/aperture factor applied to CTPH. |

#### Systematic error

| Option | Argument | Default | Description |
|---|---|---|---|
| `--syserr` | `auto` or mag | `auto` | `auto`: fit the excess scatter of the calibration stars together with CTPH and add it in quadrature to every output error. A value ≥ 0 [mag] is applied as given, with CTPH from the robust ratio. `0` turns it off and gives the errors of older versions. |

#### Colour terms

| Option | Argument | Default | Description |
|---|---|---|---|
| `--color-term` | `A,B[,errA,errB]` | none | Fit a colour term in the catalogue colour c = A − B, with optional error columns. Repeat the option for up to three colours. Needs `-c`; cannot be combined with `-r` or `--tratab`. |
| `--color-order` | `1`, `2`, or a list | `1` | Order of each colour term. A single value applies to all terms; a list gives one value per `--color-term`, in the same order (`2,1`). |
| `--color-ref` | `mean`, mag, or a list | `mean` | Reference colour c_ref at which CTPH is defined. `mean` is the weighted mean colour of the calibration stars. Takes one value or a list (`mean,0.5`). |
| `--color-err-correct` | `yes`, `no` | `yes` | Correct the colour fit for the errors of the catalogue colours (regression dilution). Only applies when the error columns are given. |

#### Photometric systems and output quantities

| Option | Argument | Default | Description |
|---|---|---|---|
| `--photsys-ref` | name | – | Reference (standard) photometric system of the catalogue. Required. The bundled table has `Johnson`, `Gunn` and `Stromgren` (`--list`); other systems need their own `--phsystab`. |
| `--photsys-instr` | name | `PHOTSYS` of the frames | Instrumental photometric system of the frames. Checked against `PHOTSYS1` of a `--tratab` table. |
| `--area` | m² | `AREA` of the frames | Area of the telescope's entrance aperture. Also written to the header. |
| `-q`, `--quantity` | list | `FLUX,MAG` | Output quantities: `PHOTON`, `PHOTRATE`, `PHOTNU`, `PHOTLAM`, `FLUX`, `FNU`, `FLAM`, `MAG`, `ABMAG`, `STMAG` (see below). |
| `--tratab` | file | none | Instrumental-to-standard transformation table made by `phfotran`. |
| `--phsystab` | file | bundled `photosystems.fits` | Table of photometric system definitions (filters, effective wavelengths, widths, zero-magnitude fluxes). |
| `--list` | – | – | List the available photometric systems and exit. |

#### Output

| Option | Argument | Default | Description |
|---|---|---|---|
| `--advanced` | – | off | Also write the extensions of the input frame and `PHRES`, the residuals of the calibration stars. The file is roughly twice as large. |
| `-t`, `--target-directory` | dir | – | Write the outputs into this directory, keeping the input names. |
| `-O` | – | off | Turn on advanced output file name processing (`--pattern`, `--mask`, `--format`). |
| `--pattern` | regex | – | Matching pattern for the input names (with `-O`). |
| `--mask` | mask | – | Output name mask, e.g. `'\1_MAG.\2'` (with `-O`). |
| `--format` | format | – | Output file format (with `-O`). |
| `--verbose` | – | off | Print the fit: CTPH, the systematic error with its interval, the outlier fraction and the colour coefficients. |
| `-h`, `--help` | – | – | Show the help. |

### The calibration

**Zero point.** For a star of catalogue magnitude m in a band of effective wavelength λ, zero-magnitude flux density f_λ and effective width W = √(π/2) × FWHM (from the photometric system table), the expected photon rate is ϕ = f_λ W (λ/hc) 10^(−0.4 m). The observed counts are converted to rates by the exposure time and the area. CTPH is the ratio of photons to counts. For counts in electrons it is the inverse of the efficiency of the telescope, filter and detector (CTPH ≥ 1 for a correct area). The gain is not applied, so for counts in ADU CTPH is about the gain divided by that efficiency. Every pixel and every star is then converted as n = CTPH · counts. The comment block of the PHOTOMETRY header lists every calibration star with its photons, its rate and its relative residual.

**Which photometry.** The aperture flux is only a fraction of the total flux, and that fraction changes over the field when the PSF does. Use one of these ways to calibrate the total flux:

- `aphot` with an enclosed-energy model, then the default `--apcorr auto`
- `gphot`, then `--photometry grow`
- `psfphot`, then `--photometry psf`

A constant `--apcorr` factor corrects only the mean aperture loss.

**Systematic error (`--syserr auto`).** Each star's log ratio ln(photons/counts) is modelled either as a calibration star, N(ln CTPH, e² + σ²), or as an outlier from a broad distribution (a misidentification, blend, variable or saturated star). The model is a mixture likelihood (Hogg, Bovy & Lang 2010). EM starts from the median, and σ comes from the restricted likelihood, so it is not biased low for few stars. No star is clipped, and outliers carry no weight in CTPH or σ. Below 16 stars (plus the number of colour coefficients), the outlier model is fixed and a weak prior on σ selects the mode of the likelihood. This keeps σ² unbiased while gross outliers are still rejected. The header records the following:

- `PHSYSERR`: σ.
- `PHSYSERE`: the 68% upper error of σ, or its upper limit when σ = 0.
- `PHSYSNST`: the number of stars used.

CTPHERR and the colour-term errors include the uncertainty of σ (the Knapp & Hartung covariance with the floor of Röver et al.).

```{note}
With fewer than 10 calibration stars `phcal` warns that σ is little more than a guess. For sparse fields, for example space-UV frames with only a handful of catalogue stars, estimate σ from many frames of the same setup (e.g. the mean of PHSYSERR² over the frames) and pass it as a fixed value (`--syserr 0.02`).
```

**Colour terms (`--color-term`).** The passband of an instrument is never exactly that of the catalogue, so the counts-to-photons ratio depends on the colour of the star. With `--color-term A,B` the calibration is

    m_cat − m_inst = ZP + Σ_j [ k1_j (c_j − c_ref,j) + k2_j (c_j − c_ref,j)² ]

where c_j is a catalogue colour, the quadratic terms are present only with `--color-order 2`, and no cross terms are fitted. The fit uses the same mixture likelihood as the systematic error. Colour errors are propagated, including their covariance with the catalogue magnitude of the calibration band and between colours that share a band (u − g and g − r). The term multiplies the photons, so all quantities stay consistent. CTPH is the calibration at the reference colour. With `--color-ref mean` it is nearly uncorrelated with k1; use a fixed value (`--color-ref 0`) when CTPH must be compared between fields.

How each output star is corrected:

- A star with a catalogue colour gets its own correction and the error of the calibration at its colour.
- A star without a colour gets the mean correction of the calibration stars, with their scatter as its error.
- With several colours, a missing colour is predicted from the known ones by the colour–colour relation of the calibration stars.

Only stars with all their colours calibrate. The colours must be independent: u − g, g − r, u − r is refused. Use several colours when the band depends on more than one, for example a near-UV band against SDSS u. Two correlated colours need many calibration stars (tens): their coefficients are strongly anti-correlated (`COLTC12` ≈ −0.95). With several filters, every filter gets its own colour term in the same colour.

**Errors of the catalogue colours (`--color-err-correct`).** A catalogue colour scatters about the true colour by its error. This biases k towards zero by the reliability ratio λ = var(true)/var(observed) (regression dilution). It also adds a bias when the colour contains the calibration band (u − g against u). With the error columns and the default `yes`, the moments of the weighted fit are corrected (Fuller 1987), for one or several colours and for order 2. The coefficient errors grow by about 1/λ. The correction is capped when the colour range is hardly wider than its errors; `phcal` then warns, and `COLTCAP` < 1. The correction assumes the catalogue was measured with the quoted errors. Use `--color-err-correct no`, or omit the error columns, when:

- the catalogue colours are exact, or more precise than their errors say, e.g. frames simulated directly from the catalogue magnitudes;
- the simulation draws the true magnitudes about the catalogue; it is over-corrected, strongly for two correlated colours;
- the errors are not measurement errors of the true colour, e.g. the scatter of a colour transformed from another system.

**Several filters.** Each filter is calibrated by its own catalogue magnitude, independently, as if the frames were processed one at a time. The only difference is that the calibration stars, and the stars of the output tables, must be present on all frames. Without `--tratab` or `--color-term` there is no colour correction between filters.

**Transformation table (`--tratab`).** The table comes from `phfotran`, and its filters and instrumental system (`PHOTSYS1`) must match the calibration.

- **Magnitude table** (`TRAMODEL = 'MAG'`): CTPH of every filter is fitted with the table's colour terms fixed. Every star's standard magnitudes are then solved from all its filters together. MAGERR includes the coefficient errors.
- **Flux table** (`TRAMODEL = 'FLUX'` or no keyword): the forward matrix converts the expected photons into instrumental counts for the calibration, so CTPH is a correction near 1. Its inverse converts the counts of all stars into standard photons.

Stars missing on a neighbouring filter get no calibrated values.

### Output quantities

All quantities derive from the photon count n = CTPH · counts. W is the effective width of the passband and W_ν = c W / λ². Undefined magnitudes are 99.99999, with the error 9.99999.

| Quantity | Meaning | Unit |
|---|---|---|
| `PHOTON` | photon count n | photon |
| `PHOTRATE` | photon flux ϕ = n / (A T) | photon/s/m² |
| `PHOTNU` | photon rate per frequency ϕ / W_ν | photon/s/m²/Hz |
| `PHOTLAM` | photon rate per wavelength ϕ / W | photon/s/m²/nm |
| `FLUX` | energy flux in the band F = ϕ hc / λ | eV/s/m² |
| `FNU` | flux density per frequency F / W_ν | W/m²/Hz |
| `FLAM` | flux density per wavelength F / W | W/m²/nm |
| `MAG` | magnitude in the system, −2.5 log₁₀ F / (f_λ W) | mag |
| `ABMAG` | AB magnitude, −2.5 log₁₀ f_ν / 3.631×10⁻²³ W/m²/Hz | mag |
| `STMAG` | ST magnitude, −2.5 log₁₀ f_λ / 3.6335×10⁻¹¹ W/m²/nm | mag |

### Outputs

One calibrated FITS file per input frame (`frame_cal.fits` by default). The input files are not changed, but existing output files are overwritten.

- **Primary image**: converted to photons. Header `BUNIT`, `PHOTSYS` (the reference system), and the HST-like `PHOTFLAM`, `PHOTZPT`, `PHOTPLAM`, `PHOTBW`; `AREA` updated when `--area` is given.
- **`PHOTOMETRY` table**: `RAJ2000`, `DEJ2000`, `SKY`, `SKYERR` (photons per pixel), `PHOTON`, `PHOTONERR`, then a pair `Q`, `QERR` for each `-q` quantity (e.g. `MAG`, `MAGERR`, `FLUX`, `FLUXERR`). Errors include CTPHERR, PHSYSERR and the colour errors. With colour terms the columns `COLOR`, `COLORERR`, `COLORSRC` and, with several colours, `COLOR2`, `COLOR2ERR` (`COLOR3`, …) follow. `COLORSRC` is 1 for a catalogue colour and 0 for none; with several colours it is a bit mask (bit j − 1 set for catalogue colour j); with a magnitude `--tratab` it is 2 (solved from the filters).
- **Header of `PHOTOMETRY`**:
  - calibration: `CTPH`, `CTPHERR`
  - systematic error: `PHSYSERR`, `PHSYSERE`, `PHSYSNST`
  - aperture photometry only: the aperture correction `APCORTYP` (`EE`, `CONST`, `NONE`), `APCORR`, `EESOURCE`, and the aperture geometry `APER`, `ANNULUS1`, `ANNULUS2` [deg]
  - a comment block listing the calibration stars
- **Colour-term keywords**:
  - mode and colour: `COLTMODE` (`CATALOGUE` or `TRANSFORM`), `COLTERM`, `COLTREF`
  - coefficients: `COLTK1`, `COLTK1E`, `COLTK2`, `COLTK2E`, `COLTCOV`
  - calibration-star colours: `COLTMEAN`, `COLTSTD`, `COLTNSTA`
  - error correction: `COLTDILU`, `COLTCAP`, `COLTLAM` (with the error columns only)
  - several colours: `COLTNTRM`, `COLTERMj`, `COLTREFj`, `COLTjK1`, `COLTjK1E`, …, `COLTMEAj`, `COLTSTDj`, `COLTLAMj` and the correlations `COLTC12`, `COLTC13`, `COLTC23`
- **`PHRES`** (with `--advanced`): RA, Dec, photons, counts and the relative residual of every calibration star.

### Examples

Calibrate one V frame against a catalogue cone and write magnitudes and energy fluxes:

```bash
munipack cone -c UCAC4 -o cat.fits -r 0.2 -- 110.47 71.34
munipack phcal -c cat.fits --col-ra RAJ2000 --col-dec DEJ2000 \
    --photsys-ref Johnson -f V --col-mag Vmag --col-magerr e_Vmag \
    -q MAG,FLAM,ABMAG frame_V.fits
```

Calibrate the total flux of PSF photometry for B and V frames, in any order. Missing header values are given on the command line:

```bash
munipack phcal --photometry psf -c cat.fits --photsys-ref Johnson \
    --photsys-instr MYCAM --area 0.3 \
    -f V,B --col-mag Vmag,Bmag --col-magerr e_Vmag,e_Bmag \
    -t calibrated frame_B.fits frame_V.fits
```

Fit a linear colour term in B − V, and a near-UV band against SDSS u with two colours. The near-UV band is not in the bundled table, so the second command gives its own system table (`mysystems.fits`, with a system `MYUV` and a filter `nuv`). It keeps the error correction off because the catalogue colours are exact (a simulation):

```bash
munipack phcal -c cat.fits --photsys-ref Johnson -f V --col-mag Vmag --col-magerr e_Vmag \
    --color-term Bmag,Vmag,e_Bmag,e_Vmag --verbose frame.fits

munipack phcal -c sdss.fits --phsystab mysystems.fits --photsys-ref MYUV -f nuv \
    --col-mag upmag --col-magerr e_upmag \
    --color-term upmag,gpmag,e_upmag,e_gpmag --color-term gpmag,rpmag,e_gpmag,e_rpmag \
    --color-err-correct no frame.fits
```

Use a fixed systematic error on a sparse field, and calibrate a later frame relative to a calibrated one:

```bash
munipack phcal -c cat.fits --photsys-ref Johnson -f V --col-mag Vmag --syserr 0.02 frame_001.fits
munipack phcal --photsys-ref Johnson -r frame_001_cal.fits frame_002.fits
```

## phfotran

`phfotran` finds the transformation between the instrumental filters of a telescope and the standard filters of a photometric system. It uses frames of a calibration field (e.g. a Landolt or Stetson standard field), one frame per filter. It writes the table that `phcal --tratab` applies. Two models are available:

- **`--model mag`** (the default) gives standard colour terms of the magnitudes, filter by filter, fitted as `phcal` colour terms (mixture likelihood, systematic error, propagated colour errors).
- **`--model flux`** gives a linear mixing of photon rates, fitted as a nearly tridiagonal matrix.

The inputs are the same as for `phcal`: frames with astrometry and instrumental photometry, paired with `-f` and the catalogue columns by their `FILTER` keyword.

### Synopsis

```bash
munipack phfotran [options] file(s)
```

### Options

| Option | Argument | Default | Description |
|---|---|---|---|
| `-c`, `--cat` | file | – | Reference catalogue of the calibration field (as for `phcal`). |
| `-f`, `--filters` | list | – | Filters, one per frame (the frames are paired by `FILTER`). |
| `--col-ra` | name | `RAJ2000` | Right Ascension column. |
| `--col-dec` | name | `DEJ2000` | Declination column. |
| `--col-mag` | list | – | Magnitude columns, in the order of `-f`. |
| `--col-magerr` | list | none | Magnitude error columns, in the order of `-f`. |
| `--model` | `mag`, `flux` | `mag` | `mag`: colour terms of the magnitudes. `flux`: linear mixing of photon rates. |
| `--color-order` | `1`, `2` | `1` | Order of the colour terms of the `mag` model. |
| `-C`, `--cal` | list | 1 | Factors multiplying the counts of the filters, in the order of `-f`. |
| `-E`, `--extin` | list | 0 | Extinction coefficients k of the filters: counts × exp(k × airmass), giving the transformation outside the atmosphere. The airmass needs the geographic coordinates of the observatory in the headers; without them it is undefined (`AIRMASS = -1`). |
| `--photsys-ref` | name | – | Standard photometric system of the catalogue (`PHOTSYS2`). |
| `--photsys-instr` | name | `PHOTSYS` of the frames | Name of the instrumental system (`PHOTSYS1`). `phcal` checks it against its own `--photsys-instr` or against the frames. |
| `--area` | m² | `AREA` of the frames | Area of the telescope's entrance aperture. |
| `--tol` | deg | HWHM of the frames | Cross-match radius. |
| `--phsystab` | file | bundled | Photometric system definitions. |
| `--list` | – | – | List the available photometric systems and exit. |
| `--vocat` | file | bundled | VO catalogue configuration (column IDs). |
| `-o`, `--output` | file | `phfotran.fits` | Output table (overwritten). |
| `--verbose` | – | off | Print the coefficients, the systematic error, the outlier fraction and the scatter with and without the terms. |
| `-h`, `--help` | – | – | Show the help. |

### Models

**Magnitude model (`--model mag`).** For every filter n:

    m_std − m_inst = ZP + k1 (c − c_ref) [+ k2 (c − c_ref)²]

where m_inst = −2.5 log₁₀ of the count rate [counts/s/m²], c is the standard colour of the filter and its neighbour in the order of the system (B − V for B and V, V − R for R, R − I for I), and c_ref is the mean colour of the stars. This is the usual transformation of standard photometry, and its coefficients are comparable to published instrument colour terms. `phcal --tratab` fits its own zero points with the coefficients fixed.

**Flux model (`--model flux`).** The count rates c in the instrumental filters are a linear combination of the photon rates p in the standard filters, c_i = Σ_j T_ij p_j. The forward matrix T holds the efficiency on its diagonal and the band mixing off it. The backward matrix is fitted independently but `phcal` does not use it. The flux model is exact for a passband that is a mix of standard ones (a leak of a neighbouring band). The magnitude model is better for shifted or reshaped passbands and needs `--color-order 2` for a strong mixing of 10% or more.

```{note}
The numbers of instrumental and standard filters must be equal. Stars are selected by fixed limits (counts above 5 × sky in the aperture, relative error below 0.1), and the colours of the magnitude model are fixed (a filter and its neighbour). Tables written before the magnitude model existed have no `TRAMODEL` and are read as flux tables.
```

### Outputs

A FITS file (`phfotran.fits` by default). The primary header holds `PHOTSYS1` (instrumental), `PHOTSYS2` (standard), `TRAMODEL` (`MAG` or `FLUX`), `AIRMASS` (the mean of the frames), `REFCI` (flux model) and the list of input frames in comments.

- `MAG` model: extension `FOTRAN_MAG`, one row per filter with columns `FILTER`, `COLOR_A`, `COLOR_B`, `CREF`, `ZP`, `ZPERR`, `K1`, `K1ERR`, `K2`, `K2ERR`, `SYSERR`, `RMS`, `NSTARS`.
- `FLUX` model: extensions `FOTRAN_FORWARD` (photons → counts), `FOTRAN_BACKWARD`, `FOTRAN_FORWARD_ERR`, `FOTRAN_BACKWARD_ERR`; one column and one row per filter.
- Both: `FOTRANdata`, the data of the fit (RA, Dec, catalogue photons `PH_f`, `ePH_f` and count rates `CTS_f`, `eCTS_f` per filter). A table given as the only input, with no catalogue, is fitted again from these data.

### Examples

Fit the magnitude transformation of an instrument on the standard field T Phe, then calibrate frames with it:

```bash
munipack phfotran --area 1.86 --photsys-instr DK154 --photsys-ref Johnson \
    -c T_Phe.fits --col-ra RA --col-dec DEC -f B,V,R,I --col-mag B,V,R,I \
    T_Phe_B.fits T_Phe_V.fits T_Phe_R.fits T_Phe_I.fits

munipack phcal --tratab phfotran.fits --photsys-instr DK154 --photsys-ref Johnson \
    -c field.fits --col-ra RA --col-dec DEC -f B,V,R,I --col-mag B,V,R,I \
    field_B.fits field_V.fits field_R.fits field_I.fits
```

Fit the photon-mixing matrix outside the atmosphere for an instrument with a band leak:

```bash
munipack phfotran --model flux -E 0.25,0.15,0.10 -c cat.fits --photsys-ref Johnson \
    -f B,V,R --col-mag Bmag,Vmag,Rmag --col-magerr e_Bmag,e_Vmag,e_Rmag \
    -o tra_flux.fits B.fits V.fits R.fits
```

Fit quadratic colour terms and inspect the table:

```bash
munipack phfotran --color-order 2 --verbose -c cat.fits --photsys-ref Johnson \
    -f B,V --col-mag Bmag,Vmag --col-magerr e_Bmag,e_Vmag -o tra.fits B.fits V.fits
munipack fits -lt "tra.fits[FOTRAN_MAG]"
```
