# Utilities — fits, timeseries

Two helpers for frames you have already processed. `fits` looks inside FITS
files and edits them; `timeseries` gathers values from many frames into one
table, usually a light curve.

(ref-fits)=
## fits — inspect and edit FITS files

Use `fits` to see what is in a file, read or change header keywords, cut out
part of an image or table, or turn a file into editable text and back. It
does one operation per run; with no option it lists the HDUs of the file.

```bash
munipack fits [options] file(s)
```

File names accept the CFITSIO extended syntax: an extension
(`'frame_cal.fits[PHOTOMETRY]'`), an image section
(`'frame.fits[100:199,100:199]'`) or a row filter (`'frame.fits[FIND][PEAKRATIO>0.5]'`).
Quote these names and write them without spaces, because a space separates a
source from a target: `'source.fits target.fits'` writes the result to a new
file instead of changing the input.

### Options

#### Listing

| Option | Description |
|---|---|
| (none) | List the HDUs: index, name, type and size. |
| `-lh`, `--header` | Print the header (primary, or the extension you select). |
| `-lt`, `--table` | Print the rows of a table extension. |
| `-li`, `--image` | Print image pixels, one per line, as `x y value`. |
| `-K`, `--keys` `key[,key,…]` | Print the records of these keywords. |
| `--value` | With `-K`: print only the values. |
| `--shell` | With `-K`: print `KEY=VALUE` lines ready for `eval`. |

#### Editing the header

| Option | Description |
|---|---|
| `--update` | Add or update a record given by `--key`/`--val`/`--com`, or by `--templ`. |
| `--key` keyword | Keyword of the record to update. |
| `--val` value | Value of the record (type follows FITS rules, see below). |
| `--com` text | Comment of the record, e.g. `"[s] exposure time"`. |
| `--templ` file | Template file with `KEY = value / comment` lines, to update many records at once. |
| `--remove-keys` `key[,key,…]` | Remove these keywords. |

```{note}
The shell eats one level of quotes, so a string value needs two:
`--val "'M 67'"` stores the string `'M 67'`. Without inner quotes,
`--val 30` stores an integer and `--val 30.0` a real number.
```

#### Extensions, copies and text form

| Option | Description |
|---|---|
| `--remove-extensions` `extname[,…]` | Remove the extension with this name. |
| `--cat` | Copy the input to the output, applying any extended syntax (sections, filters). |
| `--dump` | Write the whole file as text to `name.lst`. |
| `--restore` | Build a FITS file from a `.lst` text file. |

#### Output names and general

| Option | Description |
|---|---|
| `-t`, `--target-directory` dir | Write results to this directory and leave the input untouched. |
| `-O` | Turn on custom output names (`--pattern`, `--mask`, `--format`). |
| `--pattern` regex | Regular expression matched against input names (default `(.+)\.(.+)`). |
| `--mask` string | Output name built from the pattern groups, e.g. `'\1_copy.\2'`. |
| `--format` format | Numbered output names, e.g. `'frame%03d.fits'`. |
| `-h`, `--help` | Print the option summary. |
| `--verbose` | Print diagnostic messages. |

Header edits change the file in place unless you give `-t` or a target name.

### Examples

Check the header of a frame, or pull one value into a shell variable:

```bash
munipack fits -K DATE-OBS,EXPTIME,FILTER frame.fits
exptime=$(munipack fits -K EXPTIME --value frame.fits)
```

Set the object name and fix an exposure time:

```bash
munipack fits --update --key OBJECT --val "'NGC 188'" --com "object name" frame.fits
munipack fits --update --key EXPTIME --val 60.0 --com "[s] exposure time" frame.fits
```

Cut a sub-image into a new file, or list only the bright stars of a table:

```bash
munipack fits --cat 'frame.fits[100:499,100:499] cut.fits'
munipack fits -lt 'frame_cal.fits[PHOTOMETRY][MAG<14]'
```

### Notes

- `--dump` and `--restore` are the easiest way to make a small catalogue or
  test file: dump an existing table, edit the text, restore it. If you add
  rows by hand, update `NAXIS2` to the new number of rows.
- `--remove-extensions` removes only the first name in a list and then exits
  with an error, even though that removal worked. Remove one extension per
  run and check the result with `munipack fits file.fits`.

(ref-timeseries)=
## timeseries — light curves and other time series

`timeseries` goes through a set of frames and writes one table row per frame:
the time of the exposure, any header values you ask for, and table columns
(such as magnitudes) of the objects you choose. Run it on calibrated frames
from `phcal` to get light curves, on `APERPHOT` or `PSFPHOT` tables for
instrumental counts, or with no objects at all to follow airmass or other
conditions through the night.

```bash
munipack timeseries [options] file(s)
```

Objects are given by coordinates and matched to the nearest star in each
frame. The time comes from `DATE-OBS` and `EXPTIME`; airmass and horizontal
coordinates also need `LONGITUD`, `LATITUDE` and a WCS. If your camera uses
other keyword names, set the environment variables `FITS_KEY_DATEOBS`,
`FITS_KEY_EXPTIME`, `FITS_KEY_LONGITUDE` or `FITS_KEY_LATITUDE`.

### Options

#### Objects

| Option | Description |
|---|---|
| `-c`, `--coo` `"x,y x,y …"` | Object coordinates, one pair per object, in one quoted argument (RA,Dec in degrees by default). |
| `--cat` file | Read the objects from a FITS catalogue instead (see Notes). |
| `--coo-type` `DEG`, `RECT` | Sky coordinates (default) or pixel coordinates. |
| `--coo-col` `col,col` | Table columns to match against (default `RAJ2000,DEJ2000`; use `X,Y` for pixels). |
| `--tol` value | Match radius, in degrees or pixels; by default set from the frame. |
| `--extname` name | Table extension to read (default `PHOTOMETRY`). |

```{warning}
Leave out `--coo-type` for sky coordinates: an explicit `--coo-type DEG`
finds no objects in this build. Use `--coo-type RECT` only for pixels.
```

#### What to list

| Option | Description |
|---|---|
| `-l`, `--col` `col[,col,…]` | Table columns to list for each object, e.g. `MAG,MAGERR`. |
| `-K`, `--keys` `key[,key,…]` | Header keywords to list for each frame, e.g. `FILTER,AIRMASS`. |
| `--enable-filename` | Add the file name. |
| `--enable-horizon` | Add azimuth and zenith distance of the field centre. |
| `--enable-airmass` | Add the airmass of the field centre. |
| `--disable-timetype` | Leave out the time column. |

#### Time

| Option | Description |
|---|---|
| `-T`, `--time-type` | `JD` (default), `MJD`, `HJD` or `PHASE`. |
| `--time-stamp` | Which moment of the exposure: `MID` (default), `BEGIN` or `END`. |
| `--lc-epoch` JD | Reference epoch for the phase; without `-T` it switches to `PHASE`. |
| `--lc-period` days | Period for the phase. |

#### Output and general

| Option | Description |
|---|---|
| `-o`, `--output` file | Output FITS table (default `timeseries.fits`). |
| `--stdout` | Also print one line per frame on the terminal. |
| `-h`, `--help` | Print the option summary. |
| `--verbose` | Report objects or keywords missing in a frame. |

### Output

The table is in the first extension, so read it as `'timeseries.fits[1]'`.
Columns are named `TIME`, then the keywords, then `<COL>_<n>` for each listed
column and object, e.g. `MAG_1`, `MAGERR_1`, `MAG_2`. An object not found in a
frame gets 0 and a warning.

### Examples

Light curves of two stars from calibrated frames:

```bash
munipack timeseries -c "180.00896,29.95754 179.98240,30.04617" \
    -l MAG,MAGERR -o lc.fits frame_*_cal.fits
munipack fits -lt 'lc.fits[1]'
```

Observing conditions through the night, no objects needed:

```bash
munipack timeseries -K FILTER,EXPTIME -T MJD --enable-filename \
    --enable-airmass --stdout -o night.fits frame_*.fits
```

Instrumental counts by pixel position, for aligned frames without astrometry:

```bash
munipack timeseries --extname APERPHOT --coo-type RECT --coo-col X,Y --tol 3 \
    -c "172.07,47.15 254.85,366.20" -l APCOUNT1,APCOUNTERR1 -o counts.fits frame_*.fits
```

### Notes

- For accurate phases, list `-T JD` (or `HJD`) and fold the light curve
  yourself: `--lc-epoch` is stored with limited precision and gets rounded
  to a whole day.
- A catalogue for `--cat` is a FITS table with four columns, `RAJ2000`,
  `DEJ2000` (deg), `pmRA` and `pmDE` (both in `arcsec/year`, or both in
  `deg/day` for Solar-system bodies). Zero proper motions are fine for most
  stars. Make one with `fits --restore` from a text file.
