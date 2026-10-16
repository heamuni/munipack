(tutorial-single-image)=
# Processing a single image

This tutorial takes one raw CCD frame all the way to calibrated magnitudes.
You will make master bias, dark and flat frames, calibrate the science frame,
find and measure the stars, fit the astrometry against a catalogue, and
calibrate the photometry.

You need Munipack (see {doc}`../install`; the commands below assume the
`munipack` alias from that page) and network access for the catalogue. Work
in an empty directory, because the commands read and write files there.

## Get example data

To try the tutorial without your own data, simulate everything with
{doc}`../routines/artificial`. The field is the open cluster M67 seen in V.
First download a catalogue of the field (UCAC4, with B and V magnitudes); it
serves as the star list for the simulation and later as the reference for
astrometry and photometry:

```bash
munipack cone -r 0.2 -o cone.fits -- 132.85 11.81
```

Without network access, use any FITS table of the field instead and name its
columns with `--col-ra`, `--col-dec`, `--col-mag` and so on (see
{doc}`../routines/cone`). For `artificial`, such a table also needs
proper-motion columns (`--col-pm-ra`, `--col-pm-dec`); zeros will do.

All frames come from the same camera: 1024 × 1024 pixels on a 20 cm
telescope, with gain, read noise, bias and dark current, and a Gaussian PSF
(which also keeps the files small). A shell variable saves typing:

```bash
CAM="--width 1024 --height 1024 --diameter 0.2 --gain 2 --readnoise 5 \
  --bias 300 --dark 0.2 -B 16 --psf GAUSS --hwhm 1.5 --no-watermark"
```

So that the flat has something to correct, make a mild vignetting of the
optics: a smooth frame that is brighter in one corner, which `flat`
normalises to a mean of 1.

```bash
munipack artificial --width 1024 --height 1024 --maglim 0 --sky-mag 15 \
  --sky-grad-x 0.0003 --sky-grad-y 0.0002 --disable-noise --psf GAUSS \
  --no-watermark --mask vignetting.fits
munipack flat -o optics.fits vignetting.fits
```

Now the calibration frames and the science frame. `--maglim 0` leaves out
the stars, `--sky-mag 99` the sky. The bias frames have zero exposure, the
darks match the 30 s science exposure, and the flats are short exposures of a
bright uniform sky seen through the vignetting:

```bash
munipack artificial $CAM --maglim 0 --sky-mag 99 --exptime 0 \
  --count 7 --seed 11 --mask 'bias_?.fits'
munipack artificial $CAM --maglim 0 --sky-mag 99 --exptime 30 \
  --count 7 --seed 21 --mask 'dark_?.fits'
munipack artificial $CAM --maglim 0 --sky-mag 10 --exptime 1 -f V \
  --flat-model-file optics.fits --count 7 --seed 31 --mask 'flat_V_?.fits'
munipack artificial $CAM -c cone.fits -f V --rcen 132.85 --dcen 11.81 \
  --scale 3600 --angle 1 --exptime 30 --flat-model-file optics.fits \
  --seed 1 --mask m67_V.fits
```

A raw frame from a telescope usually knows only where the telescope pointed,
not an exact WCS. Remove the WCS that `artificial` wrote, so that the
astrometry has real work to do:

```bash
munipack fits --remove-keys CTYPE1,CTYPE2,CRPIX1,CRPIX2,CRVAL1,CRVAL2,CD1_1,CD1_2,CD2_1,CD2_2 m67_V.fits
```

With your own data, skip this section: put your bias, dark, flat and science
frames in the directory and use their names in the steps below. The `FILTER`
keyword of your flats and science frame should be set.

## 1. Make the master frames

Averaging many calibration frames gives masters with little noise. The dark
is made with the bias removed, so it can be scaled to any exposure time.

```bash
munipack bias -o bias.fits bias_?.fits
munipack dark -bias bias.fits -o dark.fits dark_?.fits
munipack flat -bias bias.fits -dark dark.fits -o flat_V.fits flat_V_?.fits
```

You get `bias.fits`, `dark.fits` and `flat_V.fits`; the flat has a mean of 1
and should look like the vignetting you put in.

```{tip}
Expose your flats to well below saturation, about half of the range at
most, and make one flat per filter.
```

## 2. Calibrate the science frame

`phcorr` subtracts the bias and the scaled dark, divides by the flat and
converts the frame to electrons.

```bash
munipack phcorr -bias bias.fits -dark dark.fits -flat flat_V.fits \
  -O --mask '\1_c.\2' m67_V.fits
```

The raw frame stays untouched; the result is `m67_V_c.fits`, with a `STDERR`
extension holding the error of every pixel. The next steps add their
results to this file as new extensions.

## 3. Find the stars

```bash
munipack find -f 3.5 m67_V_c.fits
```

`-f` is the FWHM of the stars in pixels (about 3.5 for the simulated frame);
a value close to the real one works best. The detected stars go into a `FIND`
extension.

## 4. Measure the stars

Aperture photometry sums the light of each star in a series of circular
apertures. `--apcorr auto` also measures how the light grows with the
aperture, so that the total flux of every star is known.

```bash
munipack aphot --apcorr auto m67_V_c.fits
```

The result is an `APERPHOT` extension with the counts in every aperture.

## 5. Solve the astrometry

Give the pointing of the telescope (a little off, as in real life) and the
plate scale in pixels per degree; 3600 means 1″ per pixel. `astrometry`
matches the stars to the catalogue and fits the WCS.

```bash
munipack astrometry -c cone.fits --rcen 132.86 --dcen 11.80 --scale 3600 m67_V_c.fits
```

The WCS is written into the header of `m67_V_c.fits`, followed by a short
report in `COMMENT` lines: the rotation, the scale and the RMS of the fit.
Here the rotation should come out close to 1° and the scale close to 1″ per
pixel, as simulated. If your frame
already has an approximate WCS, leave out `--rcen`, `--dcen` and `--scale`;
it is then used as the starting guess.

## 6. Calibrate the magnitudes

`phcal` compares the measured fluxes with the catalogue magnitudes and finds
the zero point. `--syserr auto` (the default) estimates any scatter beyond
the quoted errors, rejects outliers such as variable stars, and adds that
scatter to every error.

```bash
munipack phcal -c cone.fits --photsys-ref Johnson -f V \
  --col-mag Vmag --col-magerr e_Vmag --syserr auto m67_V_c.fits
```

This writes a new file, `m67_V_c_cal.fits`: the image in photons and a
`PHOTOMETRY` table with the coordinates, fluxes and magnitudes of all stars.
`phcal` may print a note about `IEEE_DENORMAL` floating-point exceptions; it
comes from the Fortran runtime and is harmless.

Real filters never match the catalogue exactly, so the zero point depends a
little on the colour of the star. A colour term in B − V corrects this:

```bash
munipack phcal -c cone.fits --photsys-ref Johnson -f V \
  --col-mag Vmag --col-magerr e_Vmag --syserr auto \
  --color-term Bmag,Vmag,e_Bmag,e_Vmag --color-err-correct no \
  -O --mask '\1_ct.\2' m67_V_c.fits
munipack fits -K COLTERM,COLTK1,COLTK1E 'm67_V_c_ct.fits[PHOTOMETRY]'
```

Here the slope `COLTK1` comes out close to zero, because the frame was
simulated from `Vmag` itself. `--color-err-correct no` suits simulated
frames; leave it out for real data.

## 7. Look at the results

```bash
munipack fits -lh m67_V_c.fits                                     # header with the WCS
munipack fits -K CTPH,CTPHERR,PHSYSERR 'm67_V_c_cal.fits[PHOTOMETRY]'  # calibration
munipack fits -lt 'm67_V_c_cal.fits[PHOTOMETRY]'                   # the table as text
```

`CTPH` is the calibration (counts per photon), and `PHSYSERR` the systematic
error that was added; it is zero here, because the simulated stars scatter
no more than their errors. The table has one row per star, and its last two
columns are `MAG` and `MAGERR`. A few stars show 99.99999, which means that
they could not be measured, usually because they sit at the edge of the
frame or next to a brighter star.

To check the result, match the table with `cone.fits` by coordinates, for
example in TOPCAT. For the brighter stars `MAG` agrees with `Vmag` to better
than a hundredth of a magnitude, all over the frame, which shows that the
flat removed the vignetting. TOPCAT, or `astropy.io.fits` in Python, also
reads every other table in these files.

## 8. PSF photometry (optional)

Instead of apertures, you can fit a model of the star image to every star,
which copes better with crowded fields. `psfphot` uses the stars from `find`
and the sky from `aphot`.

```bash
munipack psfphot m67_V_c.fits
munipack phcal -c cone.fits --photsys-ref Johnson -f V \
  --col-mag Vmag --col-magerr e_Vmag --photometry psf \
  -O --mask '\1_psf.\2' m67_V_c.fits
```

The fits go into a `PSFPHOT` extension, and the calibrated table into
`m67_V_c_psf.fits`, so the aperture result is kept.

## Next steps

- Each routine has its own page with every option: {doc}`../routines/bias`,
  {doc}`../routines/dark`, {doc}`../routines/flat`, {doc}`../routines/phcorr`,
  {doc}`../routines/find`, {doc}`../routines/aphot`,
  {doc}`../routines/astrometry`, {doc}`../routines/phcal`,
  {doc}`../routines/psfphot` and {doc}`../routines/fits`.
- {ref}`calibration-sequence` explains the two ways of using bias and dark
  frames.
- {doc}`../concepts` describes the files, extensions, units and errors.
