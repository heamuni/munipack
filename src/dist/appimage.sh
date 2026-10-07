#!/bin/sh
#
# AppImage packager, https://appimage.org/
#
# Usage:
#
# $ appimage.sh srcdir
#
#  *  the output package will be placed in the curent dir as Munipack-*.AppImage
#
# Description in detail:
#
#     https://munipack.physics.muni.cz/appimage.html
#
# These software tools are reqiured for the build, at least:
#  * g++ and gfortran compilers
#  * gtk-dev package
#  * make, cmake
#  * fuse, https://github.com/AppImage/AppImageKit/wiki/FUSE
#
# Mandatory distdir content (use the latest versions):
#  munipack-*.tar.gz
#  cfitsio*.tar.gz
#  wxWidgets-*.tar.bz2
#  plplot-*.tar.gz
#  oakleaf-*.tar.gz
#  minpack-*.tar.gz
#
# AppImage utilities:
#
# https://github.com/AppImage/AppImageKit/releases/download/continuous/appimagetool-x86_64.AppImage
# https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/linuxdeploy-x86_64.AppImage
# https://raw.githubusercontent.com/linuxdeploy/linuxdeploy-plugin-gtk/master/linuxdeploy-plugin-gtk.sh

set -xe

DISTDIR=$1

if [ -z "$DISTDIR" ]; then
    echo "$0 distdir"
    exit 0
fi

# environment setup
TOPDIR=$PWD
APPDIR=AppDir
PREFIX=$TOPDIR/$APPDIR
NCPU=$(nproc)

LDFLAGS="-L$PREFIX/lib"
PKG_CONFIG_PATH=$PREFIX/lib/pkgconfig
LD_LIBRARY_PATH=$PREFIX/lib
export LDFLAGS PKG_CONFIG_PATH LD_LIBRARY_PATH

# remove older builds
rm -rf $PREFIX/

# cFITSIO, http://heasarc.gsfc.nasa.gov/fitsio/
tar zxf $DISTDIR/cfitsio-*.tar.gz
cd cfitsio-*/
./configure --prefix=$PREFIX --enable-reentrant
make
make install
cd $TOPDIR

# plplot, https://plplot.sourceforge.net/
tar zxf $DISTDIR/plplot-*.tar.gz
cd plplot-*/
mkdir build_dir
cd build_dir
cmake -DCMAKE_INSTALL_PREFIX=$PREFIX \
  -DENABLE_DYNDRIVERS:BOOL=OFF -DHAVE_QHULL:BOOL=OFF \
  -DDEFAULT_NO_BINDINGS:BOOL=ON -DENABLE_cxx:BOOL=ON \
  -DENABLE_wxwidgets:BOOL=OFF $TOPDIR/plplot-*/
make
make install
cd $TOPDIR
PLDIR=$(ls $PREFIX/share/ | grep plplot)

# minpack, http://deb.debian.org/debian/pool/main/m/minpack/minpack_19961126+dfsg1.orig.tar.gz
tar zxf $DISTDIR/minpack_*.tar.gz
cd minpack-*/
chmod +x configure
./configure --prefix=$PREFIX
make install
cd $TOPDIR


# OakLeaf, https://integral.physics.muni.cz/oakleaf
tar zxf $DISTDIR/oakleaf-*.tar.gz
cd oakleaf-*/
./configure FCFLAGS="-O2 -ffpe-summary=invalid,zero,overflow" \
	    --prefix=$PREFIX
make install
cd $TOPDIR


# wxWidgets, http://www.wxwidgets.org/downloads/
tar jxf $DISTDIR/wxWidgets-*.tar.bz2
cd wxWidgets-*/
./configure --with-zlib=builtin \
    --with-libpng=builtin --with-libjpeg=builtin --with-libtiff=builtin \
    --with-regex=builtin --with-expat=builtin  \
    --enable-graphics_ctx  --enable-stl --enable-baseevtloop \
    --disable-compat28 \
    --without-opengl --without-libmspack  \
    --disable-xrc --disable-aui --disable-richtext  \
    --disable-ribbon --disable-stc --disable-propgrid  \
    --prefix=$PREFIX
make -j $NCPU
make install
cd $TOPDIR

# set PATH to use wx-config created by us in precedence to the system wide one
PATH=$PREFIX/bin:$PATH

# Munipack
STATIC_LDFLAGS="-static-libgfortran -static-libgcc -static-libstdc++"
tar zxf $DISTDIR/munipack-*.tar.gz
cd munipack-*/
./configure CPPFLAGS="-I$PREFIX/include" \
	    FCFLAGS="-O2 -ffpe-summary=invalid,zero,overflow -I$PREFIX/include" \
	    CFLAGS="-O2 -DNDEBUG" CXXFLAGS="-O2 -DNDEBUG" \
	    LDFLAGS="-L$PREFIX/lib $STATIC_LDFLAGS" --prefix=$PREFIX
make -j $NCPU
make install
cd $TOPDIR


cd $PREFIX
cat >AppRun  <<EOF
#!/bin/sh

#echo \$APPIMAGE
#echo \$APPDIR
#echo \$OWD
#echo \$ARGV0

LD_LIBRARY_PATH=\$APPDIR/lib
PLPLOT_LIB=\$APPDIR/share/$PLDIR
MUNIPACK_LIBEXEC_PATH=\$APPDIR/libexec/munipack
MUNIPACK_DATA_DIR=\$APPDIR/share/munipack/
MUNIPACK_HTML_DIR=\$APPDIR/share/doc/munipack
MUNIPACK_ICON_DIR=\$APPDIR/share/pixmaps/munipack
export LD_LIBRARY_PATH PLPLOT_LIB MUNIPACK_LIBEXEC_PATH
export MUNIPACK_DATA_DIR MUNIPACK_HTML_DIR MUNIPACK_ICON_DIR

ARGV0NAME=\$(basename \$ARGV0)
if [ \$ARGV0NAME = "munipack" ]; then
   exec \$APPDIR/bin/munipack "\$@"
else
   exec \$APPDIR/bin/xmunipack "\$@"
fi
EOF
chmod +x AppRun
ln -s share/applications/xmunipack.desktop
ln -s share/icons/hicolor/128x128/apps/munipack_icon.png
ln -s share/metainfo/cz.muni.physics.munipack.metainfo.xml share/metainfo/xmunipack.appdata.xml
cd $TOPDIR

# appimagetool-x86_64.AppImage $APPDIR
DEPLOY_GTK_VERSION=3 linuxdeploy-x86_64.AppImage --appdir $APPDIR --plugin gtk --output appimage
