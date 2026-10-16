(fits)=
# fits

`fits` shows what is in a FITS file, reads or changes header keywords, cuts
out part of an image or table, and turns a file into editable text and back.

## Usage

```bash
munipack fits [options] file(s)
```

It does one operation per run. With no option it lists the HDUs of the file.

File names accept the CFITSIO extended syntax: an extension
(`'frame_cal.fits[PHOTOMETRY]'`), an image section
(`'frame.fits[100:199,100:199]'`) or a row filter (`'frame.fits[FIND][PEAKRATIO>0.5]'`).
Quote these names and write them without spaces, because a space separates a
source from a target: `'source.fits target.fits'` writes the result to a new
file instead of changing the input.

## Options

### Listing

| Option | Default | Description |
|---|---|---|
| (none) | — | List the HDUs: index, name, type and size. |
| `-lh`, `--header` | off | Print the header (primary, or the extension you select). |
| `-lt`, `--table` | off | Print the rows of a table extension. |
| `-li`, `--image` | off | Print image pixels, one per line, as `x y value`. |
| `-K`, `--keys` *key[,key,…]* | — | Print the records of these keywords. |
| `--value` | off | With `-K`: print only the values. |
| `--shell` | off | With `-K`: print `KEY=VALUE` lines ready for `eval`. |

### Editing the header

| Option | Default | Description |
|---|---|---|
| `--update` | off | Add or update a record given by `--key`/`--val`/`--com`, or by `--templ`. |
| `--key` *keyword* | — | Keyword of the record to update. |
| `--val` *value* | — | Value of the record (type follows FITS rules, see below). |
| `--com` *text* | — | Comment of the record, e.g. `"[s] exposure time"`. |
| `--templ` *file* | — | Template file with `KEY = value / comment` lines, to update many records at once. |
| `--remove-keys` *key[,key,…]* | — | Remove these keywords. |

```{note}
The shell eats one level of quotes, so a string value needs two:
`--val "'M 67'"` stores the string `'M 67'`. Without inner quotes,
`--val 30` stores an integer and `--val 30.0` a real number.
```

### Extensions, copies and text form

| Option | Default | Description |
|---|---|---|
| `--remove-extensions` *extname[,…]* | — | Remove the extension with this name (see Notes). |
| `--cat` | off | Copy the input to the output, applying any extended syntax (sections, filters). |
| `--dump` | off | Write the whole file as text to `name.lst`. |
| `--restore` | off | Build a FITS file from a `.lst` text file. |

### Output names and general

| Option | Default | Description |
|---|---|---|
| `-t`, `--target-directory` *dir* | — | Write results to this directory and leave the input untouched. |
| `-O` | off | Turn on custom output names (`--pattern`, `--mask`, `--format`). |
| `--pattern` *regex* | `(.+)\.(.+)` | Regular expression matched against input names. |
| `--mask` *string* | — | Output name built from the pattern groups, e.g. `'\1_copy.\2'`. |
| `--format` *format* | — | Numbered output names, e.g. `'frame%03d.fits'`. |
| `--verbose` | off | Print diagnostic messages. |
| `-h`, `--help` | — | Show the help. |

## Output

Listings go to the terminal. Header edits and removals change the file in
place, unless you give `-t` or a target name. `--cat` writes the target
file, `--dump` writes `name.lst`, and `--restore` writes the FITS file back.

## Examples

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

## Notes

- `--dump` and `--restore` are the easiest way to make a small catalogue or
  test file: dump an existing table, edit the text, restore it. If you add
  rows by hand, update `NAXIS2` to the new number of rows.
- `--remove-extensions` removes only the first name in a list and then exits
  with an error, even though that removal worked. Remove one extension per
  run and check the result with `munipack fits file.fits`.
- To collect values from many frames into one table, use {doc}`timeseries`.
