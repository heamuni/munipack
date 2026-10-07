#!/usr/bin/env python

import sys
import numpy

def display(image,output):

    # file
    f = astropy.io.fits.open(image)

    # works on non-compressed files only
    w = f[0].header['NAXIS1']
    h = f[0].header['NAXIS2']
    img = f[0].data

    # growth photometry extension
    data = f['GROWPHOT'].data
    f.close()

    nrows = f['GROWPHOT'].header['NAXIS2']
    
    # select objects
    growflag = data['GROWFLAG']
    xcol = data['X']
    ycol = data['Y']
    
    # converts intensity to interval 0 .. 1
    med = numpy.median(img)
    mad = numpy.median(abs(img - med))
    img = (img - med - 0.5*mad)/(30*mad)
    numpy.clip(img,0.0,1.0,img)

    # convert the image to sRGB (gamma transfer)
    img = 1.055*img**(1.0/2.4) - 0.055 
    numpy.clip(img,0.0,1.0,img)

    # plot
    plot.axis([0,w,0,h])
    plot.set_cmap('Greys')    
    plot.imshow(img,origin='lower')
    x = []
    y = []
    for i in range(nrows):
        if growflag[i] == 1:
            x.append(xcol[i])
            y.append(ycol[i])
    x = numpy.array(x)
    y = numpy.array(y)
    n = len(x)
    # these offsets centers coordinates and the image each other
    l1, = plot.plot(x-1,y+1,'o',mfc='none',mec='b',label='accepted stars ('+str(n)+')')
    x = []
    y = []
    for i in range(nrows):
        if growflag[i] == 2:
            x.append(xcol[i])
            y.append(ycol[i])
    x = numpy.array(x)
    y = numpy.array(y)

    l2, = plot.plot(x,y,'D',mfc='none',mec='r',label='non-stars')
    plot.legend(handles=[l1,l2],numpoints=1,ncol=2,bbox_to_anchor=(0.,-0.15,1.0,0.1),
                mode="expand",loc=3,frameon=False)
    plot.title(image)
    plot.savefig(output)

    # http://matplotlib.org/users/image_tutorial.html



def curve(filename):

    # file
    f = astropy.io.fits.open(filename)

    # flags
    growflag = f['GROWPHOT'].data['GROWFLAG']

    # data and residuals
    data = f['GROWDATA'].data
    nrows = f['GROWDATA'].header['NAXIS2']
    naper = f['GROWDATA'].header['NAPER']
    raper = []
    for i in range(naper):
        key = 'APER' + str(i+1)
        raper.append(float(f['GROWDATA'].header[key]))

    grow = numpy.zeros((nrows,naper))
    res = numpy.zeros((nrows,naper))
    for i in range(naper):
        col = 'GROWCURVE' + str(i+1)
        grow[:,i] = f['GROWDATA'].data[col]
        col = 'RESGROW' + str(i+1)
        res[:,i] = f['GROWDATA'].data[col]

    f.close()

    raper = numpy.array(raper)
    zero = numpy.zeros(naper)
    x = numpy.zeros(nrows*naper)
    y = numpy.zeros(nrows*naper)
    z = numpy.zeros(nrows*naper)
    f = numpy.zeros(nrows*naper,numpy.ubyte)
    n = 0
    for i in range(nrows):
        q = extinvdist(zero,0.03*raper)
        for j in range(naper):
            if grow[i,j] > 0:
                x[n] = raper[j] + q[j]
                y[n] = grow[i,j]
                z[n] = res[i,j]
                f[n] = growflag[i]
                n = n + 1

    for i in range(n):
        print("{0:5.2f} {1:7.4f} {2:7.4f} {3}".format(x[i],y[i],z[i],f[i]))


        
def extinvdist(mean,sig):

    # returns random variable from exponential distribution
    # https://en.wikipedia.org/wiki/Laplace_distribution
    
    #real, intent(in) :: mean,sig
    #real :: x

    x = numpy.random.uniform(size=len(sig))
    x = x - 0.5
    x = mean - sig*numpy.sign(x)*numpy.log10(1-2*abs(x))
    return x

            
if __name__ == "__main__":

    try:
        import astropy.io.fits
    except:
        print("{0}".format("Required python module `astropy' is missing. Please install the module."))
        sys.exit(0)


    try:
        import matplotlib.pyplot as plot
    except:
        print("{0}".format("Required python module `matplotlib' is missing. Please install the module."))
        sys.exit(0)

        
    try:
        if sys.argv[1] == "curve":
            curve(sys.argv[2])

        if sys.argv[1] == "display":
            display(sys.argv[2],sys.argv[3])
            
    except:
        print("{0}".format("Usage: python grow_report.py [curve|display] image.fits [output.png]"))
        sys.exit(0)
