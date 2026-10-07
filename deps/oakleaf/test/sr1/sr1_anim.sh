
set -xe

SRC=../../src/
FLAGS="-std=f2008 -Wall -g -p -no-pie -fcheck=all -fimplicit-none -fbacktrace -ffpe-trap=invalid,zero,overflow"

gfortran $SRC/rfun_hampel.f08 $SRC/rfun.f08 $SRC/rlike.f08 ../noise.f08 \
	 -o sr1_anim sr1_anim.f08


#if [ ! -f /tmp/tsr1mean ]; then
#    ./testsr1mean.sh
#fi

DIR=/tmp/sr1mean_anim
mkdir -p $DIR
rm -f $DIR/a?????.png

./sr1_anim $DIR


# animate
#mkdir -p $DIR
#for A in `seq 0 359`; do
#    F=$(printf "$DIR/a%03d.png" $A)
#gnuplot <<EOF
#set term pngcairo size 512,512
#set output '$F'
#set view 60,$A,,
#unset colorbox
#unset key
#splot '$FCN' w pm3d
#EOF
#done

cd $DIR
png2yuv -I p -f 30 -b 0 -j a%05d.png > a.yuv
vpxenc --good --cpu-used=0 --auto-alt-ref=1 --lag-in-frames=16 --end-usage=vbr --passes=2 --threads=4 --target-bitrate=3000 --fps=71700/1700 -o a.webm a.yuv
ffmpeg -y -i a.webm -crf 2 a.mp4
cd -
