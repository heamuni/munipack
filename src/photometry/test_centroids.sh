
# Test determination of centroids

set -e

ECC=0
INCL=0

munipack artificial --psf seeing --verbose --hwhm 1 --exptime 10 --maglim 14 \
	 --eccentricity $ECC --inclination $INCL
munipack fits -lt artificial.fits\[CATALOGUE\] > art.txt

munipack find -f 4 artificial.fits
munipack fits -lt artificial.fits\[FIND\] > find.txt

R --vanilla <<EOF
art <- read.table("art.txt",header=FALSE)
find <- read.table("find.txt",header=FALSE)
x0 = art[,1]
y0 = art[,2]
x = find[,1]
y = find[,2]

# cross-match
cross <- matrix(0,nrow=length(x),ncol=2)
n = 0
for (i in 1:length(x)) {
  rmin = 10 # max. uncertainity
  m = 0
  for (j in 1:length(x0)) {
    r = sqrt((x[i]-x0[j])**2 + (y[i]-y0[j])**2)
    if( r < rmin ) {
      rmin = r
      m = j
    }
  }
  if( m > 0 ) {
    n = n + 1
    cross[n,1] = i
    cross[n,2] = m
  }
}
#print(cross)


dx <- vector("numeric",n)
dy <- vector("numeric",n)
for (i in 1:n) {
  if( cross[i,1] > 0 & cross[i,2] > 0 ) {
     k = cross[i,1]
     l = cross[i,2]
     dx[i] = x[k] - x0[l]
     dy[i] = y[k] - y0[l]
  }
}

png('dx.png')
hist(dx,60,main="dx: Ecc $ECC, Incl $INCL")
dev.off()
png('dy.png')
hist(dy,60,main="dy: Ecc $ECC, Incl $INCL")
#dev.off()


EOF
