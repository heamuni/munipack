(timeseries)=
# timeseries

`timeseries` goes through a set of frames and writes one table row per frame,
usually to make a light curve.

## Usage

```bash
munipack timeseries [options] file(s)
```

Each row holds the time of the exposure, any header values you ask for, and
table columns (such as magnitudes) of the objects you choose. Run it on
calibrated frames from {doc}`phcal` to get light curves, on `APERPHOT` or
`PSFPHOT` tables for instrumental counts, or with no objects at all to follow
airmass or other conditions through the night.

Objects are given by coordinates and matched to the nearest star in each
frame. The time comes from `DATE-OBS` and `EXPTIME`; airmass and horizontal
coordinates also need `LONGITUD`, `LATITUDE` and a WCS. If your camera uses
other keyword names, set the environment variables `FITS_KEY_DATEOBS`,
`FITS_KEY_EXPTIME`, `FITS_KEY_LONGITUDE` or `FITS_KEY_LATITUDE`.

## Options

### Objects

| Option | Default | Description |
|---|---|---|
| `-c`, `--coo` *"x,y x,y …"* | — | Object coordinates, one pair per object, in one quoted argument (RA,Dec in degrees by default). |
| `--cat` *file* | — | Read the objects from a FITS catalogue instead (see Notes). |
| `--coo-type` *DEG\|RECT* | `DEG` | Sky coordinates or pixel coordinates (see the warning below). |
| `--coo-col` *col,col* | `RAJ2000,DEJ2000` | Table columns to match against; use `X,Y` for pixels. |
| `--tol` *value* | from the frame | Match radius, in degrees or pixels (5 px for pixel coordinates). |
| `--extname` *name* | `PHOTOMETRY` | Table extension to read. |

```{warning}
Leave out `--coo-type` for sky coordinates: an explicit `--coo-type DEG`
finds no objects in this build. Use `--coo-type RECT` only for pixels.
```

### What to list

| Option | Default | Description |
|---|---|---|
| `-l`, `--col` *col[,col,…]* | — | Table columns to list for each object, e.g. `MAG,MAGERR`. |
| `-K`, `--keys` *key[,key,…]* | — | Header keywords to list for each frame, e.g. `FILTER,AIRMASS`. |
| `--enable-filename` | off | Add the file name. |
| `--enable-horizon` | off | Add azimuth and zenith distance of the field centre. |
| `--enable-airmass` | off | Add the airmass of the field centre. |
| `--disable-timetype` | off | Leave out the time column. |

### Time

| Option | Default | Description |
|---|---|---|
| `-T`, `--time-type` *JD\|MJD\|HJD\|PHASE* | `JD` | Kind of time in the time column. |
| `--time-stamp` *MID\|BEGIN\|END* | `MID` | Which moment of the exposure. |
| `--lc-epoch` *JD* | — | Reference epoch for the phase; without `-T` it switches to `PHASE`. |
| `--lc-period` *days* | 1 | Period for the phase. |

### Output and general

| Option | Default | Description |
|---|---|---|
| `-o`, `--output` *file* | `timeseries.fits` | Output FITS table. |
| `--stdout` | off | Also print one line per frame on the terminal. |
| `--verbose` | off | Report objects or keywords missing in a frame. |
| `-h`, `--help` | — | Show the help. |

## Output

The table is in the first extension, so read it as `'timeseries.fits[1]'`.
Columns are named `TIME`, then the keywords, then `<COL>_<n>` for each listed
column and object, e.g. `MAG_1`, `MAGERR_1`, `MAG_2`. An object not found in a
frame gets 0 and a warning.

## Examples

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

## Notes

- For accurate phases, list `-T JD` (or `HJD`) and fold the light curve
  yourself: `--lc-epoch` is stored with limited precision and gets rounded
  to a whole day.
- A catalogue for `--cat` is a FITS table with four columns, `RAJ2000`,
  `DEJ2000` (deg), `pmRA` and `pmDE` (both in `arcsec/year`, or both in
  `deg/day` for Solar-system bodies). Zero proper motions are fine for most
  stars. Make one with `fits --restore` from a text file (see {doc}`fits`).
