#!/usr/bin/env python

import sys

# https://gaia.esac.esa.int/documentation/GDR1/Data_processing/chap_cu5phot/sec_phot_calibr.html


def ucac_jmuc(ucac,jmuc):

    # converts UCAC4 r,i in Gunn to R,I of Johnson's

    fits = fitsio.open(ucac)
    rmag = fits[1].data['rmag']
    imag = fits[1].data['imag']
    e_rmag = fits[1].data['e_rmag']
    e_imag = fits[1].data['e_imag']

    # http://www.sdss.org/dr4/algorithms/sdssUBVRITransform.html
    # Lupton (2005)
    rjmag = rmag - 0.2936*(rmag - imag) - 0.1439
    ijmag = rmag - 1.2444*(rmag - imag) - 0.3820

    cr = fitsio.Column('Rmag','1E','mag',array=rjmag)
    ci = fitsio.Column('Imag','1E','mag',array=ijmag)
    cre = fitsio.Column('e_Rmag','1E','cmag',array=e_rmag)
    cie = fitsio.Column('e_Imag','1E','cmag',array=e_imag)

    # removing old ones to prevent duplicity (case sensitivity)
    d = fitsio.ColDefs(fits[1].data)
    d.del_col('rmag')
    d.del_col('imag')
    d.del_col('e_rmag')
    d.del_col('e_imag')

    # and replace its with Johnson's approximations
    d.add_col(cr)
    d.add_col(cre)
    d.add_col(ci)
    d.add_col(cie)

    # write to the new table
    f = fitsio.BinTableHDU.from_columns(d)
    f.writeto(jmuc,overwrite=True)



if __name__ == "__main__":

    try:
        import astropy.io.fits as fitsio
    except:
        print("{0}".format("Required python module `astropy' is missing. Please, install this module."))
        sys.exit(0)

    try:
        ucac_jmuc(sys.argv[1],sys.argv[2])
    except:
        print("{0}".format("Usage: python ucac_jmuc input.fits output.fits"))
        sys.exit(0)
