
# the astrometry report visualiser

import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"] = (20,15)
import numpy as np
import scipy.special
import astropy.visualization
import astropy.io.fits
import sys

def ecdf(x):
    n = np.size(x)
    xcdf = np.sort(x)
    ycdf = np.linspace(0,1,n)
    return xcdf,ycdf

fitsname = sys.argv[1]
if len(sys.argv) > 2:
    plotname = sys.argv[2]
else:
    plotname = ''

def agamma(img):
    return np.clip(1.055*np.clip(img,0,1e30)**(1/2.4) - 0.055,0,1)
perc = astropy.visualization.PercentileInterval(99.5)

fits = astropy.io.fits.open(fitsname)
obr = fits[0].data
tab = fits['ASTROREPORT'].data

# image scale
cd = np.zeros((2,2))
cd[0,0] = fits[0].header['CD1_1']
cd[0,1] = fits[0].header['CD1_2']
cd[1,0] = fits[0].header['CD2_1']
cd[1,1] = fits[0].header['CD2_2']
c = (abs(cd[0,0]) + abs(cd[1,1])) / 2
s = (abs(cd[1,0]) + abs(cd[0,1])) / 2
sc = np.sqrt(c**2 + s**2) # deg / pix

# stdsig
#print(fits[0].header["COMMENT"])
#rms = None
#for card in fits[0].header["COMMENT"]:
#    if "RMS" in card:
##        print(card.split('='))
#        line = card.split('=')[1]
##        print(line.split('['))
#        rms = float(line.split('[')[0]) / 3600
#        break
#print('rms=',rms)
rms = None
for card in fits[0].header["COMMENT"]:
    if "RMS" in card:
        sl = card.split(' ')
        for a in sl:
            if a:
                try:
                    rms = float(a) / 3600
                    break
                except ValueError:
                    continue
        break
print('rms=',rms)



# unfortunatelly, authors of astropy didn't understand the FITS structure
nrows = int(fits['ASTROREPORT'].header['NAXIS2'])
print(np.shape(tab),nrows)
mask = np.zeros((nrows))
n = 0
irow = 0
for row in tab:
    irow = irow + 1
    if row[6] > 0:
        mask[row[6]] = irow
        n = n + 1

idx = mask[1:n+1]
print(mask[1:n+1])
x = np.zeros((n))
y = np.zeros((n))
for i in range(0,n):
    irow = int(mask[i+1])
    row = tab[irow-1]
    x[i] = row[2]
    y[i] = row[3]
    print(i+1,x[i],y[i],irow)

# the plot of match
img = agamma(perc(obr))
plt.imshow(img,cmap='gray',origin='lower')
plt.scatter(x,y,s=50,facecolors='none', edgecolors='yellow')
plt.plot(x,y,color='gray')
#for i in range(0,len(mask)):
for i in range(0,len(idx)):
    if i == 0:
        color='red'
    else:
        color = 'blue'
    plt.text(x[i],y[i],str(int(idx[i])),color=color,fontsize=20)
if plotname != '':
    plt.savefig(plotname)
plt.show()

# get stars
rec = np.zeros((nrows))
dec = np.zeros((nrows))
x = np.zeros((nrows))
y = np.zeros((nrows))
dx = np.zeros((nrows))
dy = np.zeros((nrows))
q = 3600 / sc
for i in range(0,nrows):
    row = tab[i]
    rec[i] = row[0]
    dec[i] = row[1]
    x[i] = row[2]
    y[i] = row[3]
    dx[i] = q*row[4]
    dy[i] = q*row[5]

# get residuals in pixels
ares = np.zeros((nrows))
dres = np.zeros((nrows))
for i in range(0,nrows):
    row = tab[i]
#    dec = row[1]
    ares[i] = row[4] / sc #* np.sin(dec/57.3)**2
    dres[i] = row[5] / sc

# plot of astrometry stars
plt.imshow(img,cmap='gray',origin='lower')
plt.scatter(x,y,s=50,facecolors='none', edgecolors='yellow')
plt.show()

# plot of residuals
plt.scatter(ares,dres,s=10,facecolors='blue', edgecolors='blue')
plt.xlabel("residuals [pix]")
plt.ylabel("residuals [pix]")
plt.xlim([-5,5])
plt.ylim([-5,5])
plt.gca().set_aspect('equal')
plt.show()

# residuals
plt.imshow(img,cmap='gray',origin='lower')
plt.quiver(x,y,dx,dy,angles='xy',width=0.001,color='yellow')
plt.show()

# vector field of deformations


# eCDF
xres = np.linspace(np.min(ares),np.max(ares),666)
sig = rms / sc
cdf = (1 + scipy.special.erf((xres/sig)/np.sqrt(2))) / 2
cauchy = 0.5 + np.arctan2(3.14*xres,sig) / 3.14
#rdist = np.sqrt((rec-299.9329739)**2*np.cos(dec/57.3)**2 + (dec - 22.7729120)**2)
#ures = ares/np.cos(rdist/57.3)**2
#plt.plot(rdist,ares,'.')
#plt.plot(rdist,ures,'.')
#plt.show()
print(sig,sig / 3.14,rms,sc)

acdf_x, acdf_y = ecdf(ares)
dcdf_x, dcdf_y = ecdf(dres)
#fig, ax = plt.subplots()
plt.grid(True)
#plt.legend("",loc=2)
plt.xlabel("residuals [pix]")
plt.ylabel("probability")
#plt.ecdf(ares,label='R.A.')
#plt.ecdf(dres,label='Dec.')
plt.plot(acdf_x,acdf_y,label='R.A.')
plt.plot(dcdf_x,dcdf_y,label='Dec.')
plt.plot(xres,cdf,label='Normal')
plt.plot(xres,cauchy,label='Cauchy')
#plt.ecdf(ures,label='x.')
plt.legend(loc=2)
plt.xlim([-5,5])
plt.show()


# CDF for distribution of Rayleigh's and the Uniform
rres = np.hypot(ares,dres)
rcdf_x, rcdf_y = ecdf(rres)
plt.plot(rcdf_x,rcdf_y,label='Rayleigh')
#plt.ecdf(rres,label='Rayleigh')
xres = np.linspace(0,np.max(rres),666)
cdf = 1 - np.exp(-(xres/sig)**2/np.sqrt(2))
plt.plot(xres,cdf,label='CDF')
plt.legend(loc=2)
plt.xlim([-0.1,5])
plt.show()

#ures = (180/np.pi)*np.atan2(dres,ares)
ures = (180/np.pi)*np.arctan2(dres,ares)
ucdf_x, ucdf_y = ecdf(ures)
xres = np.linspace(-180,180,666)
cdf = 0.5 + xres / 360
#plt.ecdf(ures,label='Uniform')
plt.plot(ucdf_x, ucdf_y,label='Uniform')
plt.plot(xres,cdf,label='CDF')
plt.legend(loc=2)
plt.show()

#r = np.hypot(x - 2600, y - 1730)
#plt.plot(r,rres,'*')
#plt.show()
#
#plt.plot(rres,ures,'*')
#plt.show()
#
