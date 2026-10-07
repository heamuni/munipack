
export VOLTAGE=high

hist()
{
    N=$1
    NAME=$2
    R --vanilla <<EOF
data <- read.table('${NAME}',header=FALSE)
png('${NAME}_mean.png')
hist(data[,1],300,freq=FALSE)
x <- seq(-3, 3, length.out=100)
sig = sqrt($N)
y <- dnorm(x,0,1/sig)
#*sig*sqrt(2*pi)
lines(x,y,col='green')
dev.off()
png('${NAME}_stderr.png')
hist(data[,2],300)
dev.off()
png('${NAME}_stdsig.png')
hist(data[,3],300,freq=FALSE)
x <- seq(min(data[,3]), max(data[,3]), length.out=100)
dch <- function(x,df) { return(2*x*dchisq(x**2,df)) }
y = sqrt($N-1)*dch(x*sqrt($N-1),$N-1)
lines(x,y,col='green')
dev.off()
png('${NAME}_scale.png')
hist(data[,4],300)
# half-Normal distribution, mean is sqrt(2/pi) = 0.8
# for N=2
sig = 1/sqrt(2)
y <- sqrt(2/pi)/sig*exp(-x**2/2/sig**2)
lines(x,y,col='blue')
dev.off()
EOF
}

echo "smoke:" > /tmp/smoke.log
echo "N  T[s] dT[s]" >> /tmp/smoke.log
for N in 2 3 5 10 20 100; do
#for N in 2 3 4 5 6 7 8 9 10 20 50 100; do
#for N in 10; do
    NAME="/tmp/smoke$N"
    T1=$(date +'%s')
    smoke $N $NAME > /tmp/smoke${N}.log
    T=$(tail /tmp/smoke${N}.log | awk -F'=' '{if($0~/^Total/) { print $2; }}')
    T2=$(date +'%s')
    hist $N $NAME
    DT=$(echo "$T2 - $T1" | bc)
    echo $N $T $DT >> /tmp/smoke.log
    rm -f $NAME
done
