# Quick start

This page takes you from an empty directory to calibrated magnitudes on one frame. The frame is simulated, so you know the true magnitudes and can check the result yourself. The field is the open cluster M67, seen in V.

The steps are: get a catalogue, simulate a frame, find the stars, measure them, solve the astrometry, and calibrate the photometry.

## Setup

Munipack runs from a container. This shell function lets you type `munipack ...` as if it were installed natively; it mounts the current directory and runs as your own user:

```bash
munipack() {
  docker run --rm -u "$(id -u):$(id -g)" -v "$PWD":/data \
    munipack:0.6.2-hg1653-fixes "$@"
}
munipack --help
```

Most actions add an extension to the input frame instead of writing a new file. One FITS file therefore collects the whole history of its reduction.

Some runs print a note about `IEEE_DENORMAL` floating-point exceptions. It comes from the Fortran runtime and is harmless; an action that fails exits with status 1.

## 1. Get a reference catalogue

```bash
munipack cone -r 0.2 -o cone.fits -- 132.85 11.81
```

This downloads UCAC4 stars (with APASS magnitudes) within 0.2° of the field centre into `cone.fits`. Choose the radius `-r` so that it covers your frame. Always put `--` before the coordinates, otherwise a negative declination is read as an option.

You can also use your own FITS table, for example one exported from TOPCAT or astroquery. Then tell each action which columns to use with `--col-ra`, `--col-dec`, `--col-mag` and so on. A catalogue downloaded by `cone` is recognised automatically. For your own table, `artificial` also needs proper-motion columns (`--col-pm-ra`, `--col-pm-dec`); add zero columns if the motions are unknown.

## 2. Simulate a frame

```bash
munipack artificial -c cone.fits -f V --rcen 132.85 --dcen 11.81 \
  --scale 3600 --width 1024 --height 1024 --exptime 10 \
  --psf GAUSS --hwhm 1.5 \
  --gain 2 --readnoise 5 --bias 100 --dark 0.01 -B 16 \
  --seed 1 --no-watermark --mask frame.fits
```

This renders `frame.fits`: a 1024 × 1024 image at 1″ per pixel (`--scale` is in pixels per degree), with a WCS and realistic detector noise. The bias keeps the read noise from being clipped at zero in the 16-bit output. `--no-watermark` leaves out the corner mark, which would otherwise be detected as stars.

```{note}
For `--psf GAUSS`, `--hwhm` is the Gaussian sigma in pixels, not the half width. The FWHM here is about 3.5 px.
```

## 3. Find the stars

```bash
munipack find -f 3.5 frame.fits
```

`-f` is the FWHM in pixels; give a value close to the real one. The detected stars go into a `FIND` extension of the frame.

## 4. Aperture photometry

```bash
munipack aphot --apcorr auto frame.fits
```

This measures every star in a series of apertures and writes an `APERPHOT` extension. With `--apcorr auto`, it also measures the growth curve of bright isolated stars, so that the total flux of each star can be recovered later.

## 5. Astrometry

```bash
munipack astrometry -c cone.fits frame.fits
```

This matches the stars to the catalogue and writes the fitted WCS into the frame header. The WCS already in the header is used as a starting guess.

If your frame has no WCS, give the pointing and the scale instead. The rotation of the field may be unknown:

```bash
munipack astrometry -c cone.fits --rcen 132.86 --dcen 11.80 --scale 3600 frame.fits
```

## 6. Photometric calibration

```bash
munipack phcal -c cone.fits --photsys-ref Johnson -f V \
  --col-mag Vmag --col-magerr e_Vmag --syserr auto frame.fits
```

This writes a new file, `frame_cal.fits`, with the image in photons and a `PHOTOMETRY` table of calibrated magnitudes and fluxes. `--syserr auto` also estimates how much the calibration stars scatter beyond their errors, rejects outliers such as variables, and adds that scatter to every error. In sparse fields with only a handful of stars, give a fixed value instead, for example `--syserr 0.02`.

The `FILTER` keyword of the frame must be one of the filters given by `-f`.

Compare `MAG` with the catalogue `Vmag`: since the frame was simulated from that column, they should agree closely.

If your filter differs from the catalogue band, add a colour term, for example `--color-term Bmag,Vmag,e_Bmag,e_Vmag`.

## 7. PSF photometry (optional)

Instead of apertures, you can fit a PSF to every star. `psfphot` uses the stars from `find` and the sky from `aphot`, so run it after steps 3 and 4:

```bash
munipack psfphot frame.fits
munipack phcal -c cone.fits --photsys-ref Johnson -f V \
  --col-mag Vmag --col-magerr e_Vmag --photometry psf -O --mask '\1_psf.\2' frame.fits
```

The default PSF is a Gaussian with a fitted width; `--psf MOFFAT` and `--psf EMPIRICAL` are alternatives. The `-O --mask '\1_psf.\2'` part names the output `frame_psf.fits`, so it does not overwrite the aperture result.

## Look at the results

```bash
munipack fits -K CTPH,CTPHERR 'frame_cal.fits[PHOTOMETRY]'   # selected header keywords
munipack fits -lh 'frame.fits[APERPHOT]'                      # a whole header
munipack fits -lt 'frame_cal.fits[PHOTOMETRY]'                # a table as text
```

Any FITS reader works too, such as TOPCAT or `astropy.io.fits`. Re-running an action replaces its extension, and `phcal` overwrites its output file.

## Next steps

- [Concepts](concepts.md): files and extensions, units, errors and the PSF file format.
- [Accuracy and validation](validation.md): what accuracy to expect from each method.
- The reference pages list every option: [artificial](reference/artificial.md), [find](reference/find.md), [aphot](reference/aphot.md), [psfphot](reference/psfphot.md), [gphot](reference/gphot.md), [astrometry](reference/astrometry.md), [phcal](reference/phcal.md), [catalogues](reference/catalogues.md), [preprocessing](reference/preprocessing.md), [kombine](reference/kombine.md), [utilities](reference/utilities.md).
