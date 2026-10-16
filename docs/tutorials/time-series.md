(tutorial-time-series)=
# Time series photometry

This tutorial turns a night of frames of one field into the light curve of a
variable star. You will find and measure the stars on every frame, solve the
astrometry and calibrate the photometry of the whole series, and then pick
out the variable and two comparison stars with {doc}`../routines/timeseries`.

You need Munipack (see {doc}`../install`; the commands below assume the
`munipack` alias from that page), network access for the catalogue, and
Python with astropy and matplotlib or gnuplot for the plot. Work in an empty
directory, because the commands read and write files there.

## Get example data

To try the tutorial without your own data, simulate the series with
{doc}`../routines/artificial`. The field is the open cluster M67 in V. First
download a catalogue of the field; it serves as the star list for the
simulation and later as the reference for astrometry and photometry:

```bash
munipack cone -r 0.2 -o cone.fits -- 132.85 11.81
```

This is UCAC4, with B and V magnitudes. Without network access, use any FITS
table of the field instead and name its columns with `--col-ra`, `--col-dec`,
`--col-mag` and so on (see {doc}`../routines/cone`); for `artificial` it also
needs proper-motion columns (`--col-pm-ra`, `--col-pm-dec`), and zeros will
do.

Then render 40 frames, one every 5 minutes for a bit over three hours, as
seen from a site in central Europe on a February evening:

```bash
munipack artificial -c cone.fits -f V --rcen 132.85 --dcen 11.81 \
  --scale 3600 --width 1024 --height 1024 --exptime 10 \
  --psf GAUSS --hwhm 1.5 --gain 2 --readnoise 5 --bias 100 --dark 0.01 -B 16 \
  --lc-model wave --lc-ra 132.83288 --lc-dec 11.81048 \
  --lc-mag 13 --lc-amp 0.3 --lc-per 0.09 --lc-jd0 2460351.35 \
  --atmosphere --long 16.6 --lat 49.2 --date 2024-02-10 --time 20:00:00 \
  --count 40 --timestep 300 --seed 1 --no-watermark --mask 'frame_??.fits'
```

This writes `frame_01.fits` … `frame_40.fits` and takes a couple of minutes.
The `--lc-*` options add a variable star at the position of a faint catalogue
star. Its brightness follows a sine wave between V = 12.7 and 13.3 with a
period of 0.09 days (about two hours), so the series covers one and a half
cycles, with maxima at the epoch `--lc-jd0` and one period later.
`--atmosphere` adds extinction and seeing that change with the altitude of
the field, and records the site and the airmass in the header. The small
bias level only keeps the noise from being clipped at zero in the 16-bit
files; the photometry subtracts it with the sky, so these frames need no
bias, dark or flat correction.

With your own data, put your calibrated frames (see
{ref}`tutorial-single-image` for bias, dark and flat) in the directory instead,
and download a catalogue that covers your field.

## 1. Find the stars

Every action takes many files at once, so a shell glob processes the whole
series in one command:

```bash
munipack find -f 4 frame_??.fits
```

Each frame gets a `FIND` extension with its stars. Give `-f` close to the
real FWHM; here the seeing makes it about 4 pixels.

## 2. Aperture photometry

For long series, a list file is handier than a glob: one name per line,
passed as `@list`. Make it once and reuse it:

```bash
ls frame_??.fits > frames.lst
munipack aphot --apcorr auto @frames.lst
```

Each frame gets an `APERPHOT` extension. `--apcorr auto` measures the growth
curve of bright stars on every frame, so the calibration later works with the
total flux even though the seeing changes during the night.

## 3. Astrometry

The astrometry gives every star its sky coordinates, so that the same star
can be found on every frame by its RA and Dec, even if the telescope drifts.

```bash
munipack astrometry -c cone.fits @frames.lst
munipack fits -K ASTRMETH --value frame_??.fits | wc -l
```

The simulated frames already carry an approximate WCS, which `astrometry`
uses as its starting guess. The fitted WCS replaces it in each frame header,
followed by a short report in `COMMENT` lines. The keyword `ASTRMETH` is
written only by a successful fit, so the second command counts the solved
frames and should print 40. If your own frames have no WCS, give the
pointing and plate scale with `--rcen`, `--dcen` and `--scale`, as in
{ref}`tutorial-single-image`.

## 4. Photometric calibration

`phcal` calibrates each frame against the catalogue stars in the field.
Several files given to one `phcal` call are treated as one set of frames in
*different* filters, so run it once per frame in a shell loop:

```bash
for f in frame_??.fits; do
  munipack phcal -c cone.fits --photsys-ref Johnson -f V \
    --col-mag Vmag --col-magerr e_Vmag "$f"
done
```

This writes `frame_01_cal.fits` … `frame_40_cal.fits`, each with a
`PHOTOMETRY` table of calibrated magnitudes. Because every frame is fitted on
its own against hundreds of stars, the changing extinction and transparency
are taken out frame by frame. The variable itself does not spoil the fit: it
does not match its catalogue magnitude, so the robust fit gives it no
weight. `phcal` may print a note about `IEEE_DENORMAL` floating-point
exceptions for each frame; it comes from the Fortran runtime and is harmless.

## 5. Extract the light curves

Give `timeseries` the coordinates of the stars you want, here the variable
first, then a comparison and a check star of similar brightness from the
catalogue:

```bash
munipack timeseries -c "132.83288,11.81048 132.85036,11.80608 132.85571,11.79285" \
  -l MAG,MAGERR -o lc.fits frame_??_cal.fits
munipack fits -lt 'lc.fits[1]' > lc.txt
```

`lc.fits` holds one row per frame: the mid-exposure Julian date `TIME`, then
`MAG_1`, `MAGERR_1` for the variable, `MAG_2`, `MAGERR_2` for the comparison
and `MAG_3`, `MAGERR_3` for the check star. `lc.txt` is the same table as
plain text, in the same column order.

```{tip}
For a real variable, look up its coordinates with {doc}`../routines/sesame`
or in the catalogue. A star that is not found on a frame gets 0 and a
warning; run with `--verbose` to see which frames.
```

## 6. Plot the light curve

A few lines of Python; save them as `plot_lc.py` and run
`python3 plot_lc.py` to get `lc.png`:

```python
from astropy.io import fits
import matplotlib.pyplot as plt

t = fits.getdata('lc.fits', 1)
jd = t['TIME'] - 2460351
plt.errorbar(jd, t['MAG_1'], t['MAGERR_1'], fmt='o', label='variable')
plt.errorbar(jd, t['MAG_2'], t['MAGERR_2'], fmt='s', label='comparison')
plt.errorbar(jd, t['MAG_3'], t['MAGERR_3'], fmt='^', label='check')
plt.gca().invert_yaxis()
plt.xlabel('JD − 2460351')
plt.ylabel('V [mag]')
plt.legend()
plt.savefig('lc.png')
```

Or gnuplot, from the text table:

```bash
gnuplot -e "set terminal png; set output 'lc_gnuplot.png'; set yrange [*:*] reverse; \
  plot 'lc.txt' using (\$1-2460351):2:3 with yerrorbars title 'variable'"
```

The variable traces the sine wave with its maxima where the model puts them,
and its mean is a few hundredths of a magnitude brighter than 13 because the
faint catalogue star underneath adds its light. The comparison and check
stars stay flat at their catalogue magnitudes, V = 12.85 and 12.57, to within
a few millimagnitudes.

## 7. Differential photometry

`timeseries` lists magnitudes; it does not subtract stars. Because `phcal`
already calibrates each frame against all catalogue stars, `MAG_1` is in
effect differential photometry against the whole field. For the classical
form against your chosen comparison star, subtract the columns yourself, for
example with this script (save it as `diff_lc.py` and run
`python3 diff_lc.py`):

```python
from astropy.io import fits

t = fits.getdata('lc.fits', 1)
dv = t['MAG_1'] - t['MAG_2']      # variable - comparison
dk = t['MAG_3'] - t['MAG_2']      # check - comparison
print(f"variable - comparison: {dv.min():.3f} .. {dv.max():.3f}")
print(f"check - comparison: {dk.mean():.3f} +- {dk.std():.3f}")
```

The check star minus the comparison should be constant; its scatter tells
you how precise the light curve is. Without a catalogue you can do the same
with instrumental counts from the `APERPHOT` tables, which give pixel
positions only; the {doc}`../routines/timeseries` page has an example.

````{tip}
`timeseries` can also follow the observing conditions. Without any objects
it lists header keywords per frame:

```bash
munipack timeseries -K CTPH --enable-airmass --stdout -o night.fits frame_??_cal.fits
```

Each line shows the time, the calibration `CTPH` and the airmass, and the
same table goes to `night.fits`. Here `CTPH` falls as the field rises and
the airmass drops.
````

## Next steps

- {doc}`../routines/timeseries`: phases, other time types and pixel
  coordinates.
- {doc}`../routines/phcal`: colour terms and other calibrated quantities.
- {doc}`../routines/aphot` and {doc}`../routines/psfphot`: aperture and PSF
  photometry.
- {doc}`../routines/astrometry`: frames without a WCS in the header.
- {doc}`../routines/find`, {doc}`../routines/artificial`,
  {doc}`../routines/cone` and {doc}`../routines/fits`.
