# Catalogues and Virtual Observatory

These four actions get reference catalogues and handle them:

| Action | Purpose |
|---|---|
| [`cone`](#cone) | cone search in a VizieR catalogue; writes a FITS table or a VOTable |
| [`sesame`](#sesame) | turns an object name into RA, Dec |
| [`votable`](#votable) | converts a VOTable to FITS, CSV, plain text or an SVG chart |
| [`cross`](#cross) | cross-matches two FITS tables by position |

`cone` and `sesame` need network access. `votable` and `cross` work offline.

```{note}
The examples call `munipack`. To use the container instead, replace
`munipack` with
`docker run --rm -u $(id -u):$(id -g) -v "$PWD":/data munipack:0.6.2-hg1653-fixes`.
Files are read from and written to the current directory, which the
container sees as `/data`.
```

(cone)=
## cone

`cone` sends a cone search to a VizieR server (CDS by default, or a mirror)
and saves the stars within a radius of a position. A predefined catalogue
(UCAC4 by default) can be chosen by its alias, or any VizieR catalogue by its
identifier. The rows are sorted by the catalogue's magnitude column, and
magnitude limits are sent to the server as a constraint on that column. The
result is the reference catalogue that `astrometry`, `phcal` and `artificial`
read: give those actions the catalogue's column names with `--col-ra`,
`--col-dec`, `--col-mag`, … (for example `RAJ2000`, `DEJ2000`, `pmRA`,
`pmDE`, `Vmag` for UCAC4).

### Synopsis

```bash
munipack cone [options] [--] RA DEC
munipack cone --list-catalogues | --list-servers
```

RA and Dec are the centre in decimal degrees (ICRS); sexagesimal values are
not accepted. Put a negative declination after `--` (`-- 245.9 -26.5`), or
the parser reads it as an option and stops with `Unknown option '26'`.

### Options

**Catalogue and server**

| Option | Argument | Default | Description |
|---|---|---|---|
| `-c`, `--cat` | alias | `UCAC4` | Predefined catalogue (see `--list-catalogues`). A prefix is enough (case-sensitive): `-c Gaia` selects `GaiaDR3`, and `-c UCAC` selects the first match, `UCAC5`. |
| `--id` | VizieR ID | — | Any VizieR catalogue, e.g. `V/154` (SDSS DR16) or `J/A+A/544/A31`. Do not combine it with `-c`: when both are given, `-c` wins and its catalogue is queried. Munipack does not know the catalogue's magnitude column, so give `--sort` as well when you use `--magmin`/`--magmax`; without it, the limits are ignored. |
| `--server` | alias | `CDS, France` | VizieR mirror (see `--list-servers`), selected by a case-sensitive prefix of its alias, e.g. `--server Adac`. Some mirrors are often unreachable. |
| `--url` | URL | — | Full query URL in place of the catalogue and server. `${RA}`, `${DEC}`, `${SR}` (and `${ID}` with `--id`) are replaced by the values from the command line. `-c` is ignored, and so is `--server` (an unknown alias is still an error). No magnitude column is known, so `--magmin`/`--magmax` need `--sort`. |
| `--vocat` | file | installed `VOcat_conf.xml` | Catalogue definition file: the aliases, VizieR IDs, column names and mirrors. In the container it is `/opt/munipack/share/munipack/VOcat_conf.xml`. |
| `--list-catalogues` | — | — | List the predefined catalogues and their RA, Dec, proper-motion and magnitude columns. The current catalogue is marked `*`. |
| `--list-servers` | — | — | List the VizieR mirrors. The current one is marked `*`. |

**Selection**

| Option | Argument | Default | Description |
|---|---|---|---|
| `-r`, `--radius` | deg | `0.1` (6 arcmin) | Search radius. |
| `--magmin` | mag | — | Bright limit: only stars fainter than this, measured in the sort column. |
| `--magmax` | mag | — | Faint limit: only stars brighter than this, measured in the sort column. With both limits the server receives `col=min..max`. |
| `-s`, `--sort` | column | catalogue magnitude | Sort the rows by this column, in ascending order. It is also the column that `--magmin` and `--magmax` constrain. |
| `--par` | `key=value` | — | Extra VizieR constraint appended to the query, e.g. `'Vmag=11..13'` or `'-out.max=500'`. Only the last `--par` takes effect. To pass several constraints, join them with `&`: `'Vmag=11..13&-out.max=500'`. All columns are always requested (`VERB=3`). |
| `--Johnson-patch` | — | off | UCAC4 with FITS output only: replaces the APASS Gunn `rmag`, `imag` columns by Johnson `Rmag`, `Imag` (Lupton 2005). For other catalogues it only prints a warning. |

**Output**

| Option | Argument | Default | Description |
|---|---|---|---|
| `-o`, `--output` | file | `cone.fits` | Output file. Any suffix other than `.fits` gives a VOTable. |
| `--type` | `fits`, `votable` | by the suffix of `-o` | Output format. Any value other than `fits` writes VOTable XML. |
| `--verbose` | — | off | Print the query and the HTTP response. |
| `-h`, `--help` | — | — | Show the help. |

### Output

- **FITS** (`.fits`):
  - The primary HDU has no data. Its `COMMENT` cards record the search
    centre and radius, the constraints, the full query URL and the server's
    description.
  - One binary table extension holds the result. Its `EXTNAME` is the
    VizieR ID (e.g. `I/322A`, `I/355`). `REFURL`, `RADESYS` and `JEPOCH`
    (the epoch of the positions: 2000.0 for UCAC4, 2016.0 for Gaia DR3) are
    in its header when the VOTable gives them; some catalogues have no
    epoch (GSC: no `RADESYS`, no `JEPOCH`; SDSS DR16: `JEPOCH = 0`).
  - The table has every catalogue column, plus the columns that VizieR
    computes: `_r` (distance from the centre, deg), `_RAJ2000` and
    `_DEJ2000`.
- **VOTable** (`.xml` or any other suffix): the server's VOTable, with the
  rows sorted the same way. `votable` can convert it to other formats.
- When no star matches, `cone` prints `No objects found by given constrains`,
  writes no file and exits with status 0.

### Examples

UCAC4 (the default catalogue) around M67, within 0.3°:

```bash
munipack cone -r 0.3 -o m67.fits -- 132.85 11.81
```

Gaia DR3 stars of G < 17 at a southern position, saved as a VOTable:

```bash
munipack cone -c Gaia --magmax 17 -r 0.2 -o field.xml -- 245.9 -26.5
```

Any VizieR catalogue by its ID, here SDSS DR16 with u < 20 (`--sort` gives
the column of the limit):

```bash
munipack cone --id V/154 -s umag --magmax 20 -r 0.2 -o sdss.fits -- 180.0 0.0
```

Resolve a name and search around it in one command:

```bash
munipack cone -c GSC -r 0.25 -o ngc5139.fits -- $(munipack sesame "NGC 5139")
```

To give a catalogue an alias of its own, copy `VOcat_conf.xml`, add a row to
the `CATALOGUE` table (alias, VizieR ID, RA, Dec, pmRA, pmDE, EpRA, EpDE and
magnitude columns; leave unknown columns empty), and pass the file with
`--vocat`:

```bash
munipack cone --vocat myvocat.xml -c SDSS16 --magmax 20 -r 0.2 -o sdss.fits -- 180.0 0.0
```

(sesame)=
## sesame

`sesame` resolves an object name with the CDS Sesame service (which searches
SIMBAD, NED and VizieR). It prints the ICRS position as two numbers, RA and
Dec in decimal degrees. The output can go straight into `cone`, or into
`astrometry --rcen/--dcen`.

### Synopsis

```bash
munipack sesame [options] "object name"
```

Quote names that contain spaces (`"NGC 7000"`). Otherwise the second word is
an unexpected parameter.

### Options

| Option | Argument | Default | Description |
|---|---|---|---|
| `object` | name | — | Object name (required). |
| `--verbose` | — | off | Print the request and the server's full answer. |
| `-h`, `--help` | — | — | Show the help. |

### Output

`sesame` writes `RA Dec` in degrees, seven decimals, on one line of standard
output (e.g. `132.8458333 11.8138889` for M67). It writes no file. An
unknown name prints an error (`*** Nothing found ***`) and exits with
status 1.

### Examples

```bash
munipack sesame "beta per"
# 47.0422186 40.9556467
```

Store the position in shell variables:

```bash
read -r ra dec <<< "$(munipack sesame M67)"
munipack cone -r 0.3 -o m67.fits -- "$ra" "$dec"
```

(votable)=
## votable

`votable` converts a VOTable file, such as the output of `cone -o file.xml`
or a table saved from another VO tool, to a FITS table, CSV, tab-separated
text, an SVG star chart, or VOTable XML again. The input must be a VOTable:
FITS input stops with an XML parsing error.

### Synopsis

```bash
munipack votable [options] -o output file.xml
```

### Options

**Output**

| Option | Argument | Default | Description |
|---|---|---|---|
| `-o`, `--output` | file | `votable.fits` | Output file. Always give it: without it, the output goes to `votable.fits` whatever its format (XML unless `--type` says otherwise). |
| `--type` | `fits`, `csv`, `txt`, `svg`, `xml` | by the suffix of `-o` | Output format. An unknown suffix (e.g. `.dat`) without `--type` is an error. |
| `--verbose` | — | off | Verbose messages. |
| `-h`, `--help` | — | — | Show the help. |

**SVG chart** (used with `--type svg` or `-o *.svg`)

| Option | Argument | Default | Description |
|---|---|---|---|
| `-pt` | projection | `TAN` | Projection. Only the gnomonic projection is implemented. |
| `-pa` | deg | — (give it) | RA of the chart centre. |
| `-pd` | deg | — (give it) | Dec of the chart centre. |
| `-ps` | px/deg | — (give it) | Chart scale. |
| `-pw` | px | `500` | Canvas width. It takes effect only together with `-ph`. |
| `-ph` | px | `500` | Canvas height. It takes effect only together with `-pw`. |
| `-ml` | mag | — (give it) | Reference magnitude for the symbol size. A star of magnitude *m* is drawn as a circle of radius 3·10^(0.11(ml − m)) px, so stars fainter than `-ml` are smaller, not omitted. |
| `--col-ra` | field | — (required) | RA field name in the VOTable, e.g. `RAJ2000`. |
| `--col-dec` | field | — (required) | Dec field name, e.g. `DEJ2000`. |
| `--col-mag` | field | — (required) | Magnitude field name, e.g. `f.mag`, `Vmag`. |

```{warning}
`-pa`, `-pd`, `-ps` and `-ml` have no usable defaults. If you leave them out,
the chart is drawn with undefined values (typically every star at the centre
of the canvas). The three `--col-*` names must match FIELD names exactly,
otherwise the conversion fails with `RA, Dec of magnitude column has no
data`.
```

### Output

| Type | Content |
|---|---|
| `fits` | The same layout as the FITS output of `cone`: a primary HDU with the VOTable description as `COMMENT` cards, and one binary table whose `EXTNAME` is the catalogue ID, with `REFURL`, `RADESYS`, `JEPOCH`. |
| `csv` | A header line of the field names, then one comma-separated line per row. |
| `txt` | A header line of the field names, then one tab-separated line per row. Suitable for `awk`. |
| `svg` | An `<svg>` element of `-pw` × `-ph` px with one `<circle>` per star inside the canvas. North is up and east is left: x = W/2 − ps·(α − pa)·cos δ₁, y = H/2 − ps·(δ − pd), where δ₁ is the Dec of the first row. Only the first table of the VOTable is used. |
| `xml` | The VOTable, rewritten without indentation. |

### Examples

```bash
munipack cone -c Hipparcos -r 3 -o pleiades.xml -- 56.75 24.12
munipack votable -o pleiades.fits pleiades.xml     # FITS table
munipack votable -o pleiades.csv pleiades.xml      # spreadsheet
munipack votable -o pleiades.txt pleiades.xml      # text for awk
```

A finding chart of a UCAC4 cone, 4000 px/deg, 800 × 600 px:

```bash
munipack cone -r 0.1 -o m67.xml -- 132.85 11.81
munipack votable -o m67.svg -pa 132.85 -pd 11.81 -ps 4000 -ml 16 -pw 800 -ph 600 \
    --col-ra RAJ2000 --col-dec DEJ2000 --col-mag f.mag m67.xml
```

(cross)=
## cross

`cross` matches the rows of two FITS tables by angular distance. Each row of
the second table is paired with the nearest row of the first table, provided
the distance is under `--tol`. The tables can be two catalogues, or a
catalogue and the `PHOTOMETRY` table of a calibrated frame. The result is a
table of row indices, not a merged table. Positions are compared on the
sphere (haversine formula), so RA 0°/360° and high declinations are handled
correctly.

### Synopsis

```bash
munipack cross --tol DEG [options] file1 file2
```

Each file is read from its first table extension. To read another extension,
use the CFITSIO syntax, e.g. `'frame_cal.fits[PHOTOMETRY]'` (quoted for the
shell).

### Options

| Option | Argument | Default | Description |
|---|---|---|---|
| `--tol` | deg | machine epsilon | Match radius. The default matches nothing, so the option is required in practice: e.g. `3e-4` (≈ 1.1″) for catalogues and calibrated frames. |
| `--col-ra` | `col1,col2` | `RA,RA` | RA columns of file1 and file2 (in degrees). Give both names, even when they are equal: a single name stops with a Fortran runtime error. Names are not case-sensitive. |
| `--col-dec` | `col1,col2` | `DEC,DEC` | Dec columns of file1 and file2 (in degrees). Give both names. |
| `--col-pm-ra` | `col1,col2` | — | Proper-motion columns in RA (mas/yr). Do not use them in this version; see the warning below. |
| `--col-pm-dec` | `col1,col2` | — | Proper-motion columns in Dec (mas/yr). Do not use them in this version; see the warning below. |
| `--col-epoch-ra` | column | — | Accepted but ignored. The epoch comes from the `JEPOCH` header keyword. |
| `--col-epoch-dec` | column | — | Accepted but ignored. |
| `--col-mag` | `col1,col2` | — | Magnitude columns. They must exist when given, but the current match does not use them. |
| `--ftol` | ratio | `1` | Relative flux tolerance. Accepted, but the current match does not use it. |
| `--vocat` | file | — | Accepted but ignored. |
| `-o`, `--output` | file | `cross.fits` | Output file. An existing file is overwritten. |
| `--verbose` | — | off | Verbose messages. |
| `-h`, `--help` | — | — | Show the help. |

```{warning}
The proper-motion correction of `cross` is broken in this version. It
compares the `JEPOCH` keyword (a year) with J2000 expressed as a Julian date.
With `--col-pm-ra`/`--col-pm-dec`, positions in a table that has `JEPOCH`
(every `cone` output has one) are moved by about −6700 years of proper
motion, which gives few or no matches. In a table without `JEPOCH`, the
proper motions have no effect. Rows with a null proper motion are never
matched. Leave these options out. `astrometry` and `phcal` apply proper
motions themselves.
```

A table without the `JEPOCH` keyword (e.g. a frame's `PHOTOMETRY`
extension) prints `Epoch by 'JEPOCH' keyword not found (assumed J2000.0)` on
standard error. The match is unaffected.

### Output

The output FITS file has an empty primary HDU and one binary table with
`EXTNAME = CROSSMATCH`:

| Item | Content |
|---|---|
| `FILE1`, `FILE2` | The input file names as given on the command line. |
| column 1 (`TTYPE1` = file1 name, `1J`) | Row *i*: the 1-based row of file1 matched to row *i* of file2, or 0 when no row of file1 is within `--tol`. |
| column 2 (`TTYPE2` = file2 name, `1J`) | Row *j*: the 1-based row of file2 matched to row *j* of file1, or 0. |
| rows | max(N1, N2). The rows beyond the length of the other table are 0. |

The match is not one-to-one. Two rows of file2 can pick the same row of
file1, for example close pairs that the second table resolves. Column 1
lists that row of file1 for both of them, and column 2 keeps the last one.
Hence the number of non-zero entries can differ between the two columns.

### Examples

UCAC4 against Gaia DR3 over the same field:

```bash
munipack cone -r 0.1 -o ucac4.fits -- 132.85 11.81
munipack cone -c Gaia --magmax 17 -r 0.1 -o gaia.fits -- 132.85 11.81
munipack cross --tol 3e-4 --col-ra RAJ2000,RA_ICRS --col-dec DEJ2000,DE_ICRS \
    -o ucac4_gaia.fits ucac4.fits gaia.fits
```

The stars of a calibrated frame (output of `phcal`) against the reference
catalogue:

```bash
munipack cross --tol 3e-4 --col-ra RAJ2000,RAJ2000 --col-dec DEJ2000,DEJ2000 \
    -o matched.fits 'frame_cal.fits[PHOTOMETRY]' ucac4.fits
```

To build the matched pairs from the index table, for example with Astropy:

```python
from astropy.io import fits
idx = fits.getdata("ucac4_gaia.fits", "CROSSMATCH")
u = fits.getdata("ucac4.fits", 1)
g = fits.getdata("gaia.fits", 1)
j2 = idx.field(1)[: len(u)]           # row of gaia.fits for each UCAC4 row
ok = j2 > 0
pairs = list(zip(u["f.mag"][ok], g["Gmag"][j2[ok] - 1]))
```
