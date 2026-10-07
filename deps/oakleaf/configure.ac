dnl
dnl  Oakleaf
dnl
dnl  Copyright © 2018 - 2026 F.Hroch (hroch@physics.muni.cz)
dnl
dnl  This file is part of Oakleaf.
dnl
dnl  Oakleaf is free software: you can redistribute it and/or modify
dnl  it under the terms of the GNU Lesser General Public License as published by
dnl  the Free Software Foundation, either version 3 of the License, or
dnl  (at your option) any later version.
dnl
dnl  Oakleaf is distributed in the hope that it will be useful,
dnl  but WITHOUT ANY WARRANTY; without even the implied warranty of
dnl  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
dnl  GNU Lesser General Public License for more details.
dnl
dnl  You should have received a copy of the GNU Lesser General Public License
dnl  along with Oakleaf.  If not, see <http://www.gnu.org/licenses/>.
dnl

#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([Oakleaf], [1.0.1], [hroch@physics.muni.cz],
        [oakleaf],[https://integral.physics.muni.cz/oakleaf])
AM_INIT_AUTOMAKE
LT_INIT

# The library version ***should*** be updated in src/Makefile.am!

# Checks for programs.
AC_PROG_FC
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_AWK
AC_PROG_SED

# Checks for libraries.
AC_CHECK_LIB([minpack], [lmder_])

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.
AC_FC_MODULE_FLAG
AC_FC_MODULE_OUTPUT_FLAG
AC_FC_MODULE_EXTENSION

# Checks for support of REAL128 for Fortran 2008 or later
AC_FC_SRCEXT([f08])
AC_FC_PP_SRCEXT([f08])
AC_LANG([Fortran])
AC_CACHE_CHECK(
  [for Fortran REAL128 support],
  [ac_cv_type_real128],
[AC_COMPILE_IFELSE([
program kind
 use iso_fortran_env
 real(REAL128) :: x
end program kind
],[ac_cv_type_real128=yes],[ac_cv_type_real128=no])
])
if test $ac_cv_type_real128 = yes; then
   AC_DEFINE([HAVE_REAL128])
fi

AC_ARG_ENABLE([real128],
   AS_HELP_STRING([--disable-real128],[disable support of REAL128]),
   disable_real128=disable,disable_real128=enable)
AM_CONDITIONAL([REAL128],[test x$disable_real128 = xenable && test $ac_cv_type_real128 = yes])

# clean precompiled modules of fortran
AC_SUBST(MOSTLYCLEANFILES,[*.${FC_MODEXT}])

AC_CONFIG_MACRO_DIRS([m4])

dnl There is an alternative way to generate Fortran's REAL32, REAL64
dnl and REAL128 files. It can utilise shell's loop methods, but any
dnl change in the procedures requires running autoreconf && ./configure!
dnl
dnl AC_CONFIG_COMMANDS([reals],[
dnl for T in real32 real64 real128
dnl do
dnl U=`echo $T | tr a-z A-Z`
dnl for A in rfun_hampel.f08 rfun_tukey.f08 rfun_huber.f08 \
dnl         rnfun_square.f08 median.f08 \
dnl 	rlike.f08 qsort.f08 qmed.f08 qmean.f08 fmin.f08 q2eig.f08 \
dnl 	qtrstep.f08 sqpfmean.f08 rmean.f08
dnl do
dnl     ${SED} s/REAL64/${U}/g <src/${A} >src/${A%.f08}_${T}.f08
dnl done
dnl done
dnl ])

AC_CONFIG_FILES([src/minpack/Makefile src/fmm/Makefile src/simplex/Makefile
                 src/Makefile test/Makefile doc/Makefile Makefile oakleaf.pc])
AC_OUTPUT
