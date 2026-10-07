
testphrate > testphrate.out
XERR=$?

PASS=$(awk 'BEGIN { n=0; pass=0;} {if( $0 ~/^TEST/ ) { n++; pass += $4 == "T" ? 1 : 0;}} END {print !(n == pass);}' < testphrate.out)

exit $PASS && $XERR



gnuplot <<EOF
t1 = 1e4
t2 = 1e4
sig1 = 2e3
sig2 = 2e3
Phi(x)=(1+erf(x/sqrt(2)))/2
a(x)=sqrt(x**2/sig1**2 + 1/sig2**2)
b(x)=(t1*x/sig1**2 + t2/sig2**2)
c=t1**2/sig1**2 + t2**2/sig2**2
d(x)=exp((b(x)**2 - c*a(x)**2)/a(x)**2/2)
f(x)=b(x)*d(x)/sqrt(2*pi)/sig1/sig2/a(x)**3 * \
     (Phi(b(x)/a(x)) - Phi(-b(x)/a(x))) + exp(-c/2)/pi/sig1/sig2/a(x)**2
g(x)=b(x)*d(x)/sqrt(2*pi)/sig1/sig2/a(x)**3
plot f(x), g(x), d(x)
h(x)=(t1*x*sig1**2+t2*sig2**2)/(sig1**2*x**2+sig2**2)**1.5*exp(-((t2*x-t1)/sqrt(sig1**2*x**2+sig2**2))**2/2)
plot g(x), h(x)
EOF
