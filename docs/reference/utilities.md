# Utilities — fits, timeseries

Two helper actions that work on frames already processed by the other
actions: `fits` reads and edits FITS files (headers, tables, images,
extensions), and `timeseries` collects quantities from a set of frames into a
time-series table, typically a light curve.

The examples call `munipack` directly. With the container image, define once:

```bash
alias munipack='docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes'
```

(ref-fits)=
## fits — FITS file utility

`fits` lists the structure of a FITS file, prints header records, table rows
and image pixels, adds, updates and removes header records, removes
extensions, copies or filters files through the CFITSIO extended file-name
syntax, and converts a whole file to and from a plain-text form. It is mostly
a shell wrapper of CFITSIO calls and is meant for inspection and for shell
scripts. One run does one operation; with no operation option it prints the
list of HDUs.

### Synopsis

```bash
munipack fits [options] file(s)
```

Every file argument may carry the CFITSIO extended syntax, such as
`'frame.fits[1]'`, `'frame_cal.fits[PHOTOMETRY]'`, an image section
`'frame.fits[100:199,100:199]'` or a row filter
`'frame.fits[FIND][PEAKRATIO>0.5]'`. Quote such names, since `[` and `]` are
special to the shell, and write them without spaces (`[MAG<14]`, not
`[MAG < 14]`): a space separates the source and target names (below).

An argument of the form `'source target'` (both names in one quoted argument,
separated by a space) writes the result to *target* instead of the input file;
this is how `--cat` makes a new file.

### Options

When several operation options are given, only one is done, in this order of
precedence: `-K`, `--remove-keys`, `--update`, `-lh`, `-lt`, `-li`, `--dump`,
`--restore`, `--cat`, `--remove-extensions`.

#### Listing

| Option | Argument | Default | Description |
|---|---|---|---|
| (none) | — | — | Prints the HDUs: index, `EXTNAME`, type, `BITPIX` and size. |
| `-lh`, `--header` | — | — | Prints the whole header of the selected HDU (the primary HDU unless an extension is given), as stored. |
| `-lt`, `--table` | — | — | Prints all rows of a table extension, one row per line, without column names (use `-lh` and the `TTYPEn` records to see the column order). |
| `-li`, `--image` | — | — | Prints the image pixels, one per line, as `i1 i2 … value` (1-based indices). Use an image section to limit the output. |
| `-K`, `--keys` | `key[,key,…]` | — | Prints the records of the given keywords. Records are printed in header order, not in the order given; a missing keyword prints nothing. |
| `--value` | — | off | With `-K`: prints only the values (strings keep their apostrophes). |
| `--shell` | — | off | With `-K`: prints `KEY=VALUE` lines for `eval` in a shell. |

#### Editing the header

| Option | Argument | Default | Description |
|---|---|---|---|
| `--update` | — | — | Adds a record or updates an existing one, given by `--key`, `--val`, `--com`, or several records given by `--templ`. Works on the primary header or on an extension (`'file.fits[1]'`, `'file.fits[FIND]'`); the file keeps all its extensions. |
| `--key` | keyword | — | Keyword of the updated record (with `--update`). Use upper-case letters, digits, `-` and `_`. |
| `--val` | value | — | Value of the record (with `--update`). The type follows FITS rules: `T`/`F` logical, digits only integer, a number with `.` or an exponent real, a quoted string string (see the note below). |
| `--com` | text | — | Comment of the record (with `--update`). A unit in brackets, `[s] exposure time`, is the FITS convention. |
| `--templ` | file | — | CFITSIO template file with one `KEY = value / comment` record per line (with `--update`); updates several records in one run. |
| `--remove-keys` | `key[,key,…]` | — | Removes the records of the given keywords. A missing keyword prints a FITSIO `202` message; the exit status stays 0. |

```{note}
The shell removes one level of quotes. To store a string, quote it twice:
`--val "'M 67'"` stores `'M 67'`, while `--val 30` stores the integer 30,
`--val 30.0` the real 30. and `--val "'30'"` the string `'30'`. An apostrophe
inside a string is doubled: `"'Barnard''s star'"`.
```

#### Extensions, copies and text form

| Option | Argument | Default | Description |
|---|---|---|---|
| `--remove-extensions` | `extname[,extname,…]` | — | Removes the extension(s) with the given `EXTNAME`. See the caveat below. |
| `--cat` | — | — | Copies the input (with any extended syntax applied) to the output, like CFITSIO `fitscopy`. Used to cut image sections, select HDUs, filter rows or columns of tables. |
| `--dump` | — | — | Writes the whole file (all HDUs, headers and data) as text to `name.lst` (`*.fits` → `*.lst`). Never to standard output. |
| `--restore` | — | — | Creates a FITS file from a text file in the `--dump` format (`*.lst` → `*.fits`). |

```{warning}
`--remove-extensions` removes the first listed extension, then stops with
`FITSIO status = 301: illegal HDU number` and exit status 1, even when the
removal succeeded; further names in the list are not removed. Remove one
extension per run and check the result with `munipack fits file.fits`.
```

#### Output files and general

| Option | Argument | Default | Description |
|---|---|---|---|
| `-t`, `--target-directory` | directory | — | Writes the output to this directory and leaves the input untouched (for `--update`, `--remove-keys`, `--cat`). Ignored by `--dump` and `--restore`; use `-O --mask` there. |
| `-O` | — | off | Turns on advanced output file names (`--pattern`, `--mask`, `--format`). |
| `--pattern` | regex | `(.+)\.(.+)` | Regular expression matched against the input names; its groups are `\1`, `\2`, … (for `-O`). Without `-O`, `--dump` uses `(.*).fits` and `--restore` `(.*).lst`. |
| `--mask` | string | — | Output name built from the groups of `--pattern`, e.g. `'\1_copy.\2'` (for `-O`; giving it implies `-O`). Without `-O`, `--dump` uses `\1.lst` and `--restore` `\1.fits`. |
| `--format` | printf format | — | Numbered output names, e.g. `'frame%03d.fits'`, counted from 0 in the input order (for `-O`). |
| `-h`, `--help` | — | — | Prints the option summary. |
| `--verbose` | — | off | Prints diagnostic messages. |

### Text format of `--dump` / `--restore`

Each HDU is enclosed in `# BEGIN HDU n` and `# END HDU n` lines. The header
follows the first line and ends with `END`; the data, if any, follow. Table
rows are in Fortran free format: strings in apostrophes, numbers separated by
spaces or commas. When rows or columns are added by hand, update `NAXIS2`
(rows), `TFIELDS` (columns) and `NAXIS1` (bytes per row: 8 for `D`, 4 for `E`
or `J`, *n* for `nA`).

```text
# BEGIN HDU 0
SIMPLE  =                    T
BITPIX  =                    8
NAXIS   =                    0
EXTEND  =                    T
END
# END HDU 0
# BEGIN HDU 1
XTENSION= 'BINTABLE'
BITPIX  =                    8
NAXIS   =                    2
NAXIS1  =                   16
NAXIS2  =                    2
PCOUNT  =                    0
GCOUNT  =                    1
TFIELDS =                    2
TTYPE1  = 'RAJ2000 '
TFORM1  = '1D      '
TUNIT1  = 'deg     '
TTYPE2  = 'DEJ2000 '
TFORM2  = '1D      '
TUNIT2  = 'deg     '
END
180.0090 29.9575
179.9824 30.0462
# END HDU 1
```

This is the simplest way to make a small catalogue, a coordinate list for
`timeseries --cat`, or a test frame.

### Outputs

- Listings (`-lh`, `-lt`, `-li`, `-K`, no option) go to standard output.
- `--update`, `--remove-keys` and `--remove-extensions` modify the file in
  place, unless `-t` or a `'source target'` argument names another file; then
  the whole file is copied and the copy is modified. A copy to a new file may
  print a harmless `failed to find or open` message for the target.
- `--cat` writes the target file; `--dump` writes `*.lst`; `--restore` writes
  `*.fits` (an existing file is replaced).

### Examples

Structure, header and selected keywords of a processed frame:

```bash
munipack fits frame.fits
munipack fits -lh 'frame_cal.fits[PHOTOMETRY]'
munipack fits -K DATE-OBS,EXPTIME,FILTER frame.fits

# one value per keyword, for scripts
exptime=$(munipack fits -K EXPTIME --value frame.fits)
for f in *.fits; do
    echo "$f $(munipack fits -K DATE-OBS,FILTER --value "$f" | tr '\n' ' ')"
done
```

Add and fix header records, in the primary header and in an extension:

```bash
munipack fits --update --key OBJECT --val "'NGC 188'" --com "object name" frame.fits
munipack fits --update --key EXPTIME --val 60.0 --com "[s] exposure time" frame.fits
munipack fits --update --key SATURATE --val 60000 'frame.fits[APERPHOT]'
munipack fits --remove-keys OBSERVER,TELESCOP frame.fits

# several records from a template
cat > site.tpl <<'EOF'
LONGITUD = 16.6 / [deg] geographic longitude
LATITUDE = 49.2 / [deg] geographic latitude
EOF
munipack fits --update --templ site.tpl frame.fits
```

Cut, select and filter with the extended syntax:

```bash
munipack fits -li 'frame.fits[200:210,300:305]'           # a few pixels
munipack fits -lt 'frame_cal.fits[PHOTOMETRY][MAG<14]'     # rows by a condition
munipack fits --cat 'frame.fits[100:499,100:499] cut.fits' # a sub-image to a new file
munipack fits --cat -t backup frame.fits                   # a copy in backup/
```

Edit a file as text:

```bash
munipack fits --dump catalogue.fits       # writes catalogue.lst
$EDITOR catalogue.lst
munipack fits --restore catalogue.lst     # writes catalogue.fits
munipack fits -O --pattern '(.+)\.fits' --mask '\1.txt' --dump frame.fits
```

(ref-timeseries)=
## timeseries — time series

`timeseries` lists, for each input frame, the time of the exposure, selected
header values and selected table columns of chosen objects, and stores them
as rows of a FITS table. The objects are given by coordinates (on the command
line or in a small catalogue) and are cross-matched in the table extension of
every frame. With calibrated frames (`phcal` output, extension `PHOTOMETRY`)
it gives light curves in magnitudes or fluxes; with instrumental tables
(`APERPHOT`, `PSFPHOT`) and pixel coordinates it gives instrumental series.
It can also list header quantities alone (e.g. the airmass or the sky level
during a night). It is meant for tens to thousands of frames, not for
archives.

### Synopsis

```bash
munipack timeseries [options] file(s)
```

The time is taken from `DATE-OBS` (start of the exposure, also the legacy
`DATE-OBS` + `TIME-OBS` form) and `EXPTIME`. Horizontal coordinates and
airmass also need `LONGITUD` and `LATITUDE` [deg] and a WCS. The environment
variables `FITS_KEY_DATEOBS`, `FITS_KEY_EXPTIME`, `FITS_KEY_LONGITUDE` and
`FITS_KEY_LATITUDE` rename these keywords (in the container, pass them with
`docker run -e FITS_KEY_DATEOBS=DATE-BEG …`).

### Options

#### Objects

| Option | Argument | Default | Description |
|---|---|---|---|
| `-c`, `--coo` | `"x,y x,y …"` | none | Coordinates of the objects, one `x,y` pair per object, pairs separated by spaces (or `\|`), the whole list in one quoted argument. Degrees (RA, Dec) by default, pixels with `--coo-type RECT`. Use `;` between *x* and *y* in a locale with a decimal comma. |
| `--cat` | file | — | Objects from a FITS table with proper motions (see below). Used only when `-c` is not given. |
| `--coo-type` | `DEG`, `RECT` | `DEG` (implicit) | Type of the coordinates: spherical (RA, Dec) or rectangular (pixels). For spherical coordinates leave the option out: an explicit `DEG` does not work (see the caveat below). |
| `--coo-col` | `col,col` | `RAJ2000,DEJ2000` | Columns of the frame table compared with the object coordinates, e.g. `X,Y` with `--coo-type RECT`. |
| `--tol` | value | see description | Search radius of the cross-match, in degrees (spherical) or pixels (`RECT`). Spherical default: 5 × the larger of the astrometric error and one pixel, or 5 × HWHM, or 5 pixels, by what the frame provides (1″ as the last resort). `RECT` default: 5 pixels. The nearest table object within the radius is taken. |
| `--extname` | name | `PHOTOMETRY` | Table extension read from every frame, e.g. `APERPHOT`, `GROWPHOT`, `PSFPHOT`. |

#### Quantities

| Option | Argument | Default | Description |
|---|---|---|---|
| `-l`, `--col` | `col[,col,…]` | none | Table columns to list for every object; names must match the columns of the frames exactly, e.g. `MAG,MAGERR` or `FLUX,FLUXERR`. |
| `-K`, `--keys` | `key[,key,…]` | none | Header keywords to list for every frame, e.g. `FILTER,AIRMASS`. The first HDU that has the keyword gives the value. |
| `--enable-filename` | — | off | Adds the input file name (column `FILENAME`). |
| `--enable-horizon` | — | off | Adds the azimuth and the zenith distance of the field centre (columns `AZIMUTH`, `ZENITD`, deg). |
| `--enable-airmass` | — | off | Adds the airmass of the field centre (column `AIRMASS`). |
| `--disable-timetype` | — | off | Omits the time column. |

#### Time

| Option | Argument | Default | Description |
|---|---|---|---|
| `-T`, `--time-type` | `JD`, `MJD`, `HJD`, `PHASE` | `JD` | Time quantity: Julian date (UTC), modified JD (JD − 2400000.5), heliocentric JD (corrected for the field centre `CRVAL`), or phase. Case-insensitive. |
| `--time-stamp` | `MID`, `BEGIN`, `END` | `MID` | Point of the exposure the time refers to: start + `EXPTIME`/2, start, or start + `EXPTIME`. |
| `--lc-epoch` | JD | 0 | Reference epoch *e*₀ of the light-curve elements. Giving it (or `--lc-period`) without `-T` selects `PHASE`. See the caveat below. |
| `--lc-period` | days | 1 | Period *P*. The phase is the fractional part of (*t* − *e*₀)/*P*, negative for *t* < *e*₀. |

```{warning}
`--lc-epoch` and `--lc-period` are passed to the computing module with seven
significant digits, so a Julian date is rounded to a whole day:
`--lc-epoch 2460311.4` is used as 2460311.0 and shifts every phase by
0.4/*P*. The period keeps seven digits (0.35214567 → 0.3521457). For accurate
phases, list `-T JD` (or `MJD`, `HJD`) and compute the phase outside.
```

#### Output and general

| Option | Argument | Default | Description |
|---|---|---|---|
| `-o`, `--output` | file | `timeseries.fits` | Output FITS file; an existing file is replaced. |
| `--stdout` | — | off | Also prints one line per frame to standard output (file name, time, keyword values, horizon, airmass, then the columns of each object). |
| `-h`, `--help` | — | — | Prints the option summary. |
| `--verbose` | — | off | Reports objects not found (with the radius used) and missing keywords per frame. |

```{warning}
`--coo-type DEG` given explicitly does not work in this build: the
cross-match finds no object (`Error: tolerance is unset`). Spherical
coordinates are the default, so omit the option; use `--coo-type RECT` only
for pixel coordinates. HJD, horizontal coordinates and airmass need the WCS
and are meaningless with `RECT`.
```

### Catalogue for `--cat`

A FITS file whose first table extension (any `EXTNAME`) has exactly four
columns:

| Column | Description | Unit (`TUNITn`) |
|---|---|---|
| `RAJ2000` | right ascension α₀ | deg |
| `DEJ2000` | declination δ₀ | deg |
| `pmRA` | proper motion in α (of the coordinate itself, not μ<sub>α</sub> cos δ) | `arcsec/year` or `deg/day` |
| `pmDE` | proper motion in δ | the same as `pmRA` |

Both proper-motion columns must use the same unit, `arcsec/year` (stars) or
`deg/day` (Solar-system bodies). The header keyword `EPOCH` [JD] is the
reference epoch of α₀, δ₀ (J2000.0 when missing); the position at time *t* is
α₀ + μ<sub>α</sub>(*t* − `EPOCH`), δ₀ + μ<sub>δ</sub>(*t* − `EPOCH`). Zero
proper motions are fine for most stars. The file is easily made with
`fits --restore` from a text file in the `--dump` format (see the `fits`
section); `src/doc/timeserie_cat.lst` is a template, but as distributed it
mixes the two units and its `NAXIS2` disagrees with the rows, so set both
`TUNIT3` and `TUNIT4` to the same unit and `NAXIS2` to the number of rows
(and `NAXIS1 = 32`, `TFIELDS = 4`).

### Outputs

The output file has three HDUs:

| HDU | Content |
|---|---|
| primary | empty; `HDUNAME = 'TIMESERIES'`, `CREATOR` |
| `TIMESERIES` (binary table) | one row per valid input frame, in the input order |
| `CATALOGUE` (binary table) | the object coordinates `RAJ2000`, `DEJ2000` [deg] at `EPOCH` [JD], `REFRAME = 'ICRS'`; with `--coo-type RECT` the same columns hold the pixel coordinates |

Columns of `TIMESERIES`, in this order: `FILENAME` (`--enable-filename`),
`TIME` (unit `day`, none for phase), one column per `-K` keyword (logical,
integer, real or string by the type of the header value; the unit from the
`[unit]` of its comment), `AZIMUTH`, `ZENITD`, `AIRMASS` (when enabled), and
`<COL>_<n>` for each `-l` column and object *n* = 1, 2, … in the order of `-c`
or `--cat` (e.g. `MAG_1`, `MAGERR_1`, `MAG_2`, …), with the units of the
frame table.

Header of `TIMESERIES`: `TIMETYPE` (`JD`, `MJD`, `HJD`, `PHASE`),
`TIMESYS = 'UTC'`, `TIMEREF = 'GEOCENTER'`, `TIMESTMP`, and the input file
names as `COMMENT` records. `TIMESTMP` is always written as `'MIDPOINT'`,
whatever `--time-stamp` was; the `TIME` values do follow `--time-stamp`.

An object not found in a frame gets 0 in the table, and 99.999 (`MAG`),
9.999 (`MAGERR`) or 0 in the `--stdout` listing; a warning is printed. Frames
that cannot be read are skipped.

Since the primary HDU is also named `TIMESERIES`, address the table by its
index: `'timeseries.fits[1]'`.

### Examples

Light curves of two stars from calibrated frames, magnitudes with errors:

```bash
munipack timeseries -c "180.00896,29.95754 179.98240,30.04617" \
    -l MAG,MAGERR -o lc.fits frame_*_cal.fits
munipack fits -lt 'lc.fits[1]'      # TIME MAG_1 MAGERR_1 MAG_2 MAGERR_2
```

Phased light curve of a variable, heliocentric times on the terminal:

```bash
munipack timeseries -c "180.00896,29.95754" -l MAG,MAGERR \
    --lc-epoch 2460311.4 --lc-period 0.3521 -o phase.fits frame_*_cal.fits
# (the epoch is used as 2460311.0, see the caveat under Time)
munipack timeseries -c "180.00896,29.95754" -l MAG,MAGERR -T HJD --stdout \
    -o hjd.fits frame_*_cal.fits
```

Observing conditions of a night, no objects:

```bash
munipack timeseries -K FILTER,EXPTIME -T MJD --time-stamp BEGIN \
    --enable-filename --enable-horizon --enable-airmass --stdout \
    -o night.fits frame_*.fits
```

Instrumental counts by pixel positions (frames aligned, no astrometry), and
objects from a catalogue:

```bash
munipack timeseries --extname APERPHOT --coo-type RECT --coo-col X,Y --tol 3 \
    -c "172.07,47.15 254.85,366.20" -l APCOUNT1,APCOUNTERR1 \
    -o counts.fits frame_*.fits

munipack fits --restore targets.lst          # RAJ2000, DEJ2000, pmRA, pmDE
munipack timeseries --cat targets.fits -l FLUX,FLUXERR -o flux.fits frame_*_cal.fits
```
