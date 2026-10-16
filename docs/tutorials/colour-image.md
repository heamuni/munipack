# Stacking and colour images

This tutorial takes several short exposures of one field in each of the
Johnson B, V and R filters and turns them into a deep frame per filter and
a colour picture. Along the way you calibrate the deep frames
photometrically, which also gives the colour image a correct white balance.

The steps are: calibrate the positions of every frame, stack each filter,
measure the stacks, compose the colour image, and convert it to PNG.

You need Munipack (the `munipack` alias from {doc}`../install`) and, for the
last step, the small `fitspng` converter. If you have not used Munipack
before, the {doc}`../quickstart` explains the single-frame steps in more
detail.

## Get example data

To try the tutorial without your own data, simulate the frames. First
download a catalogue of the field, the open cluster M67:

```bash
munipack cone --Johnson-patch -r 0.25 -o cone.fits -- 132.85 11.81
```

This needs a network connection. The default catalogue, UCAC4, carries B
and V magnitudes and Sloan r and i; `--Johnson-patch` converts the latter to
Johnson R and I, so `cone.fits` has `Bmag`, `Vmag` and `Rmag` columns.

Then render three 20-second exposures per filter, each pointed a few
arcseconds away from the others, as when you dither the telescope:

```bash
seed=0
for f in B V R; do
  n=0
  for pos in "132.844 11.806" "132.850 11.810" "132.856 11.814"; do
    set -- $pos                       # $1 = RA, $2 = Dec
    n=$((n + 1)); seed=$((seed + 1))
    munipack artificial -c 'cone.fits[1][Bmag<50 && Rmag<50]' -f $f --col-mag ${f}mag \
      --rcen $1 --dcen $2 --width 800 --height 800 --scale 3600 \
      --exptime 20 --psf GAUSS --hwhm 1.5 --cosmics-freq 0.1 \
      --gain 2 --readnoise 5 --bias 100 -B 16 \
      --seed $seed --no-watermark --mask m67_${f}_$n.fits
  done
done
```

This writes nine frames, `m67_B_1.fits` … `m67_R_3.fits`, of 800 × 800
pixels of one arcsecond each (`--scale` is in pixels per degree). Each frame
has its own noise (a different `--seed`) and a few cosmic rays, so the
stacking has something to reject. The filter in brackets keeps only stars
that have all three magnitudes; a star missing from one band would come out
with a false colour.

**With your own data**, put your frames in place of these, one file per
exposure with a `FILTER` keyword of `B`, `V` or `R`, and calibrate them
first with {doc}`../routines/phcorr` (bias, dark and flat). Download the
catalogue with `cone` around your own field.

## 1. Calibrate the positions

Stacking aligns the frames by their world coordinates (WCS), so every frame
needs an accurate one. The simulated frames, like those from most
telescopes, already carry a rough WCS; `astrometry` refines it. Find the
stars, measure them, and match them to the catalogue:

```bash
munipack find -f 3.5 m67_?_?.fits
munipack aphot m67_?_?.fits
munipack astrometry -c cone.fits m67_?_?.fits
```

`astrometry` writes the fitted WCS into the header of each frame. Check the
error of the fit for one of them:

```bash
munipack fits -K CRDER1,CRDER2 m67_V_1.fits
```

The values are in degrees; here they are around 0.1″, a tenth of a pixel,
which is plenty for stacking.

```{tip}
If your frames have no WCS at all, tell `astrometry` where the telescope
pointed and the plate scale, for example
`--rcen 132.85 --dcen 11.81 --scale 3600`; see {doc}`../routines/astrometry`.
```

## 2. Stack each filter

{doc}`../routines/kombine` reprojects the frames of one filter onto a common
grid and averages them. Give the same grid to every filter, so that the
three stacks line up pixel for pixel, as colouring needs. Here it is centred
on the cluster, as large as the frames, with north up (`--angle 0`) and not
mirrored (`--reflex no`):

```bash
for f in B V R; do
  munipack kombine --rcen 132.85 --dcen 11.81 --width 800 --height 800 \
    --scale 3600 --angle 0 --reflex no -o m67_$f.fits m67_${f}_?.fits
done
```

You get `m67_B.fits`, `m67_V.fits` and `m67_R.fits`. The image is the sum of
the exposures with the background removed, and the total exposure time is in
`EXPTIME`. The `STDERR` and `MAXDEV` extensions hold the error of each pixel
and its largest deviation:

```bash
munipack fits m67_V.fits
munipack fits -K FILTER,EXPTIME m67_V.fits
```

The default average is a mean clipped by the noise model of the frames,
which rejects cosmic rays and hot pixels. `--arith` gives a plain mean that
needs far less memory, for long series already free of cosmic rays. Here the
output keeps the scale of the inputs; if you resample to larger pixels and
want total star fluxes for photometry, add `--flux-conserve`.

## 3. Measure the stacks

The stacks are ordinary frames, so you can measure them like any other.
`kombine` does not copy the `GAIN` and `AREA` keywords, so look them up in a
single frame:

```bash
munipack fits -K GAIN,AREA m67_V_1.fits
```

The stack is a sum in ADU, so its gain is that of the single frames. Put it
back, so that `aphot` gets the errors right:

```bash
munipack fits --update --key GAIN --val 2.0 --com "[e-/ADU] gain" m67_B.fits m67_V.fits m67_R.fits
munipack find -f 3.5 m67_B.fits m67_V.fits m67_R.fits
munipack aphot m67_B.fits m67_V.fits m67_R.fits
```

Then calibrate each band against its catalogue column:

```bash
for f in B V R; do
  munipack phcal -c cone.fits --photsys-ref Johnson --area 0.99 -f $f \
    --col-mag ${f}mag --col-magerr e_${f}mag --syserr auto m67_$f.fits
done
munipack fits -K CTPH,CTPHERR 'm67_V_cal.fits[PHOTOMETRY]'
```

The area of the telescope, 0.99 m² from the single frames, is given with
`--area`. `phcal` may print a note about floating-point exceptions; it is
harmless.

This writes `m67_B_cal.fits`, `m67_V_cal.fits` and `m67_R_cal.fits`, with the
image in photons and a `PHOTOMETRY` table of calibrated magnitudes. The last
command shows the calibration factor `CTPH` of the V stack and its error.
The magnitudes agree with the catalogue to about a hundredth of a
magnitude.

## 4. Compose the colour image

{doc}`../routines/colouring` transforms the bands into the CIE 1931 XYZ
colour space. List them from short to long wavelength, B, V, R:

```bash
munipack colouring -o m67_colour.fits m67_B_cal.fits m67_V_cal.fits m67_R_cal.fits
munipack fits -K NAXIS3,CSPACE m67_colour.fits
```

`m67_colour.fits` is a cube with one plane per colour component; the
`CSPACE` keyword names the colour space. Because the inputs are calibrated,
the bands are already balanced and the colours are natural: the giants of
M67 come out yellow-orange, the hot stars bluish.

```{tip}
You can also colour the uncalibrated stacks `m67_B.fits` … directly, but
then balance them yourself, with `--white-star x,y` on a star of B − V near
zero or with weights `-w`. Without either, `colouring` warns that the
weights are undefined and the colours false.
```

## 5. View and convert

An ordinary FITS viewer shows the cube as three grey planes. To see the
colours, convert it to PNG with `fitspng` (packaged as `fitspng` in Debian
and Ubuntu; it is not in the Docker image, so run it on the host):

```bash
fitspng -o m67_colour.png m67_colour.fits
```

`m67_colour.png` is an sRGB image ready for any viewer or web page.
`fitspng --help` lists options for the tone curve and saturation, useful when
the faint parts of a nebula need lifting.

## Next steps

- {doc}`../routines/kombine`: output grids, mosaics, moving targets and
  resampling.
- {doc}`../routines/colouring`: colour spaces, weights and white balance.
- {doc}`../routines/astrometry` and {doc}`../routines/phcal`: matching and
  calibration options.
- {doc}`../routines/phcorr`: bias, dark and flat corrections for real frames.
- {doc}`../routines/artificial` and {doc}`../routines/cone`: simulated frames
  and catalogues.
