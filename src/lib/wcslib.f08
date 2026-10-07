
! The Fortan 2008+ interface for WCSLIB by Mark Calabretta
! https://www.atnf.csiro.au/computing/software/wcs/
! The interface is incomplete covering only wcs.h and wcsheader.h.
!
! The reference version of WCSLIB is 8.7.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   This file must be updated if some component of WCSLIB has been changed.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Copyright © 2026 F.Hroch (hroch@physics.muni.cz)
!
! This file is part of Munipack.
!
! Munipack is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! Munipack is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with Munipack.  If not, see <http://www.gnu.org/licenses/>.
!

module wcslib

  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  implicit none

  integer, parameter, private :: PVN = 30

  interface wcs_put
     module procedure wcs_put_double, wcs_put_char, wcs_put_int
  end interface wcs_put

  interface wcs_get
     module procedure wcs_get_double1, wcs_get_double, wcs_get_double2d, &
          wcs_get_char, wcs_get_char1, wcs_get_int
  end interface wcs_get

  interface dis_get
     module procedure dis_get_double1, dis_get_int
  end interface dis_get

  interface dis_put
     module procedure dis_put_int, dis_put_char
  end interface dis_put

  interface prj_put
     module procedure prj_put_double, prj_put_char, prj_put_int
  end interface prj_put

  interface prj_get
     module procedure prj_get_double, prj_get_double2d, prj_get_char, prj_get_char1
  end interface prj_get

  interface cel_put
     module procedure cel_put_double, cel_put_char, cel_put_int
  end interface cel_put

  interface cel_get
     module procedure cel_get_double, cel_get_double2d, cel_get_char, cel_get_char1
  end interface cel_get


! extracted from wcslib-8.7/Fortran/wcs_f.c:
  ! awk '{if( $0 ~/^#define/) print "integer, parameter :: " $2  "= " $3;}' < wcs_f.c
  integer, parameter :: WCS_FLAG = 100
  integer, parameter :: WCS_NAXIS = 101
  integer, parameter :: WCS_CRPIX = 102
  integer, parameter :: WCS_PC = 103
  integer, parameter :: WCS_CDELT = 104
  integer, parameter :: WCS_CRVAL = 105
  integer, parameter :: WCS_CUNIT = 106
  integer, parameter :: WCS_CTYPE = 107
  integer, parameter :: WCS_LONPOLE = 108
  integer, parameter :: WCS_LATPOLE = 109
  integer, parameter :: WCS_RESTFRQ = 110
  integer, parameter :: WCS_RESTWAV = 111
  integer, parameter :: WCS_NPV = 112
  integer, parameter :: WCS_NPVMAX = 113
  integer, parameter :: WCS_PV = 114
  integer, parameter :: WCS_NPS = 115
  integer, parameter :: WCS_NPSMAX = 116
  integer, parameter :: WCS_PS = 117
  integer, parameter :: WCS_CD = 118
  integer, parameter :: WCS_CROTA = 119
  integer, parameter :: WCS_ALTLIN = 120
  integer, parameter :: WCS_VELREF = 121
  integer, parameter :: WCS_ALT = 122
  integer, parameter :: WCS_COLNUM = 123
  integer, parameter :: WCS_COLAX = 124
  integer, parameter :: WCS_CNAME = 125
  integer, parameter :: WCS_CRDER = 126
  integer, parameter :: WCS_CSYER = 127
  integer, parameter :: WCS_CZPHS = 128
  integer, parameter :: WCS_CPERI = 129
  integer, parameter :: WCS_WCSNAME = 130
  integer, parameter :: WCS_TIMESYS = 131
  integer, parameter :: WCS_TREFPOS = 132
  integer, parameter :: WCS_TREFDIR = 133
  integer, parameter :: WCS_PLEPHEM = 134
  integer, parameter :: WCS_TIMEUNIT = 135
  integer, parameter :: WCS_DATEREF = 136
  integer, parameter :: WCS_MJDREF = 137
  integer, parameter :: WCS_TIMEOFFS = 138
  integer, parameter :: WCS_DATEOBS = 139
  integer, parameter :: WCS_DATEBEG = 140
  integer, parameter :: WCS_DATEAVG = 141
  integer, parameter :: WCS_DATEEND = 142
  integer, parameter :: WCS_MJDOBS = 143
  integer, parameter :: WCS_MJDBEG = 144
  integer, parameter :: WCS_MJDAVG = 145
  integer, parameter :: WCS_MJDEND = 146
  integer, parameter :: WCS_JEPOCH = 147
  integer, parameter :: WCS_BEPOCH = 148
  integer, parameter :: WCS_TSTART = 149
  integer, parameter :: WCS_TSTOP = 150
  integer, parameter :: WCS_XPOSURE = 151
  integer, parameter :: WCS_TELAPSE = 152
  integer, parameter :: WCS_TIMSYER = 153
  integer, parameter :: WCS_TIMRDER = 154
  integer, parameter :: WCS_TIMEDEL = 155
  integer, parameter :: WCS_TIMEPIXR = 156
  integer, parameter :: WCS_OBSGEO = 157
  integer, parameter :: WCS_OBSORBIT = 158
  integer, parameter :: WCS_RADESYS = 159
  integer, parameter :: WCS_EQUINOX = 160
  integer, parameter :: WCS_SPECSYS = 161
  integer, parameter :: WCS_SSYSOBS = 162
  integer, parameter :: WCS_VELOSYS = 163
  integer, parameter :: WCS_ZSOURCE = 164
  integer, parameter :: WCS_SSYSSRC = 165
  integer, parameter :: WCS_VELANGL = 166
  integer, parameter :: WCS_RSUN_REF = 167
  integer, parameter :: WCS_DSUN_OBS = 168
  integer, parameter :: WCS_CRLN_OBS = 169
  integer, parameter :: WCS_HGLN_OBS = 170
  integer, parameter :: WCS_HGLT_OBS = 171
  integer, parameter :: WCS_A_RADIUS = 172
  integer, parameter :: WCS_B_RADIUS = 173
  integer, parameter :: WCS_C_RADIUS = 174
  integer, parameter :: WCS_BLON_OBS = 175
  integer, parameter :: WCS_BLAT_OBS = 176
  integer, parameter :: WCS_BDIS_OBS = 177
  integer, parameter :: WCS_NTAB = 200
  integer, parameter :: WCS_NWTB = 201
  integer, parameter :: WCS_TAB = 202
  integer, parameter :: WCS_WTB = 203
  integer, parameter :: WCS_LNGTYP = 204
  integer, parameter :: WCS_LATTYP = 205
  integer, parameter :: WCS_LNG = 206
  integer, parameter :: WCS_LAT = 207
  integer, parameter :: WCS_SPEC = 208
  integer, parameter :: WCS_TIME = 209
  integer, parameter :: WCS_CUBEFACE = 210
  integer, parameter :: WCS_TYPES = 211
  integer, parameter :: WCS_LIN = 212
  integer, parameter :: WCS_CEL = 213
  integer, parameter :: WCS_SPC = 214
  integer, parameter :: WCS_ERR = 215

  integer, parameter :: WCSHDR_all = 2**20 - 1

  integer, parameter :: PRJ_FLAG = 100
  integer, parameter :: PRJ_CODE = 101
  integer, parameter :: PRJ_R0 = 102
  integer, parameter :: PRJ_PV = 103
  integer, parameter :: PRJ_PHI0 = 104
  integer, parameter :: PRJ_THETA0 = 105
  integer, parameter :: PRJ_BOUNDS = 106
  integer, parameter :: PRJ_NAME = 200
  integer, parameter :: PRJ_CATEGORY = 201
  integer, parameter :: PRJ_PVRANGE = 202
  integer, parameter :: PRJ_SIMPLEZEN = 203
  integer, parameter :: PRJ_EQUIAREAL = 204
  integer, parameter :: PRJ_CONFORMAL = 205
  integer, parameter :: PRJ_GLOBAL = 206
  integer, parameter :: PRJ_DIVERGENT = 207
  integer, parameter :: PRJ_X0 = 208
  integer, parameter :: PRJ_Y0 = 209
  integer, parameter :: PRJ_ERR = 210
  integer, parameter :: PRJ_W = 211
  integer, parameter :: PRJ_N = 212

  integer, parameter :: CEL_FLAG = 100
  integer, parameter :: CEL_OFFSET = 101
  integer, parameter :: CEL_PHI0 = 102
  integer, parameter :: CEL_THETA0 = 103
  integer, parameter :: CEL_REF = 104
  integer, parameter :: CEL_PRJ = 105
  integer, parameter :: CEL_EULER = 200
  integer, parameter :: CEL_LATPRQ = 201
  integer, parameter :: CEL_ISOLAT = 202
  integer, parameter :: CEL_ERR = 203

  integer, parameter :: DIS_FLAG   = 100
  integer, parameter :: DIS_NAXIS  = 101
  integer, parameter :: DIS_DTYPE  = 102
  integer, parameter :: DIS_NDP    = 103
  integer, parameter :: DIS_NDPMAX = 104
  integer, parameter :: DIS_DP     = 105
  integer, parameter :: DIS_TOTDIS = 106
  integer, parameter :: DIS_MAXDIS = 107

  integer, parameter :: LIN_DISSEQ = 106


!!$  type, bind(c) :: wcserr
!!$     integer(kind=C_INT) :: status
!!$     integer(kind=C_INT) :: line_no
!!$     type(c_ptr) :: function  ! character(kind=C_CHAR)
!!$     type(c_ptr) :: file      ! character(kind=C_CHAR)
!!$     type(c_ptr) :: msg       ! character(kind=C_CHAR)
!!$  end type wcserr
!!$
!!$  type, bind(c) :: pvcard
!!$     integer(kind=C_INT) :: i, m
!!$     real(kind=C_DOUBLE) :: value
!!$  end type pvcard
!!$
!!$  type, bind(c) :: pscard
!!$     integer(kind=C_INT) :: i, m
!!$     character(kind=C_CHAR), dimension(72) :: value
!!$  end type pscard
!!$
!!$  type, bind(c) :: auxprm
!!$     real(kind=C_DOUBLE) :: rsun_ref
!!$     real(kind=C_DOUBLE) :: dsun_obs
!!$     real(kind=C_DOUBLE) :: crln_obs
!!$     real(kind=C_DOUBLE) :: hgln_obs
!!$     real(kind=C_DOUBLE) :: hglt_obs
!!$     real(kind=C_DOUBLE) :: a_radius
!!$     real(kind=C_DOUBLE) :: b_radius
!!$     real(kind=C_DOUBLE) :: c_radius
!!$     real(kind=C_DOUBLE) :: blon_obs
!!$     real(kind=C_DOUBLE) :: blat_obs
!!$     real(kind=C_DOUBLE) :: bdis_obs
!!$     real(kind=C_DOUBLE), dimension(2) :: dummy
!!$  end type auxprm
!!$
!!$  type, bind(c) :: tabprm
!!$     integer(kind=C_INT) :: flag
!!$     integer(kind=C_INT) :: M
!!$     type(c_ptr) :: K     ! integer(kind=C_INT), dimension(:) :: K
!!$     type(c_ptr) :: map   ! integer(kind=C_INT), dimension(:) :: map
!!$     type(c_ptr) :: crval ! real(kind=C_DOUBLE), dimension(*) :: crval
!!$     type(c_ptr) :: index ! real(kind=C_DOUBLE), dimension(**) :: index
!!$     type(c_ptr) :: coord !r eal(kind=C_DOUBLE), dimension(*) :: coord
!!$     integer(kind=C_INT) :: nc
!!$     integer(kind=C_INT) :: padding
!!$     type(c_ptr) :: sense   ! integer(kind=C_INT), dimension(:) :: sense
!!$     type(c_ptr) :: p0      ! integer(kind=C_INT), dimension(:) :: p0
!!$     type(c_ptr) :: delta   ! real(kind=C_DOUBLE), dimension(:) :: delta
!!$     type(c_ptr) :: extrema ! real(kind=C_DOUBLE), dimension(:) :: extrema
!!$     type(c_ptr) :: err     !type(wcserr), dimension(:) :: err
!!$     integer(kind=C_INT) :: m_flag
!!$     integer(kind=C_INT) :: m_M
!!$     integer(kind=C_INT) :: m_N
!!$     integer(kind=C_INT) :: set_M
!!$     type(c_ptr) :: m_K     ! integer(kind=C_INT), dimension(:) :: m_K
!!$     type(c_ptr) :: m_map   ! integer(kind=C_INT), dimension(:) :: m_map
!!$     type(c_ptr) :: m_crval ! real(kind=C_DOUBLE), dimension(:) :: m_crval
!!$     type(c_ptr) :: m_index, m_indxs !real(kind=C_DOUBLE), dimension(**)
!!$     type(c_ptr) :: m_coord ! real(kind=C_DOUBLE), dimension(:) :: m_coord
!!$  end type tabprm
!!$
!!$  type, bind(c) :: wtbarr
!!$     integer(kind=C_INT) :: i,m
!!$     integer(kind=C_INT) :: kind
!!$     character(kind=C_CHAR), dimension(72) :: extnam
!!$     integer(kind=C_INT) :: extver
!!$     integer(kind=C_INT) :: extlev
!!$     character(kind=C_CHAR), dimension(72) :: ttype
!!$     integer(kind=C_LONG) :: row
!!$     integer(kind=C_INT) :: ndim
!!$     type(c_ptr) :: dimlen  ! integer(kind=C_INT), dimension(:) :: dimlen
!!$     type(c_ptr) :: arrayp  ! real(kind=C_DOUBLE), dimension(**) :: arrayp
!!$  end type wtbarr
!!$
  type, bind(c) :: dpkey
     character(kind=C_CHAR), dimension(72) :: field
     integer(kind=C_INT) :: i
     integer(kind=C_INT) :: type
     real(kind=C_DOUBLE) :: value ! union with i
!     integer(kind=C_INT) :: ii ! union with f
  end type dpkey

  type, bind(c) :: disprm
     integer(kind=C_INT) :: flag
     integer(kind=C_INT) :: naxis
     character(kind=C_CHAR), dimension(72) :: dtype
     integer(kind=C_INT) :: ndp
     integer(kind=C_INT) :: ndpmax
     type(c_ptr) :: dp   !type(dpkey) :: dp
     real(kind=C_DOUBLE) :: todis
     type(c_ptr) :: maxdis  ! real(kind=C_DOUBLE), dimension(:) :: maxdis
     type(c_ptr) :: docorr  ! integer(kind=C_INT), dimension(:) :: docorr
     type(c_ptr) :: Nhat    ! integer(kind=C_INT), dimension(:) :: Nhat
     type(c_ptr) :: axmap   ! integer(kind=C_INT), dimension(**) :: axmap
     type(c_ptr) :: offset  ! real(kind=C_DOUBLE), dimension(**) :: offset
     type(c_ptr) :: scale   ! real(kind=C_DOUBLE), dimension(**) :: scale
     type(c_ptr) :: iparm   ! integer(kind=C_INT), dimension(**) :: iparm
     type(c_ptr) :: dparm   ! real(kind=C_DOUBLE), dimension(**) :: dparm
     integer(kind=C_INT) :: i_naxis
     integer(kind=C_INT) :: ndis
     type(c_ptr) :: err ! type(wcserr), dimension(*) :: err
     type(c_funptr) :: disp2x   ! integer(kind=C_INT), dimension(**) :: disp2x
     type(c_funptr) :: disx2p   ! integer(kind=C_INT), dimension(**) :: disx2p
     integer(kind=C_INT) :: m_flag
     integer(kind=C_INT) :: m_naxis
     character(kind=C_CHAR), dimension(72) :: m_dtype
     type(dpkey) :: m_dp
     type(c_ptr) :: m_maxdis ! integer(kind=C_DOUBLE), dimension(*) :: m_maxdis
  end type disprm

  type, bind(c) :: linprm
     integer(kind=C_INT) :: flag
     integer(kind=C_INT) :: naxis
     type(c_ptr) :: crpix   ! real(kind=C_DOUBLE), dimension(:) :: crpix
     type(c_ptr) :: pc      ! real(kind=C_DOUBLE), dimension(:) :: pc
     type(c_ptr) :: cdelt   ! real(kind=C_DOUBLE), dimension(:) :: cdelt
     type(c_ptr) :: dispre  ! type(disprm) :: dispre
     type(c_ptr) :: disseq  ! type(disprm) :: disseq
     type(c_ptr) :: piximg  ! real(kind=C_DOUBLE), dimension(:) :: piximg
     type(c_ptr) :: imgpix  ! real(kind=C_DOUBLE), dimension(:) :: imgpix
     integer(kind=C_INT) :: i_naxis
     integer(kind=C_INT) :: unity
     integer(kind=C_INT) :: affine
     integer(kind=C_INT) :: simple
     type(c_ptr) :: err    ! type(wcserr), dimension(*) :: err
     type(c_ptr) :: dummy  ! real(kind=C_DOUBLE), dimension(:) :: dummy
     integer(kind=C_INT) :: m_flag
     integer(kind=C_INT) :: m_naxis
     type(c_ptr) :: m_crpix  ! real(kind=C_DOUBLE), dimension(:) :: m_crpix
     type(c_ptr) :: m_pc     ! real(kind=C_DOUBLE), dimension(:) :: m_pc
     type(c_ptr) :: m_cdelt  ! real(kind=C_DOUBLE), dimension(:) :: m_cdelt
     type(c_ptr) :: m_dispre ! type(disprm) :: m_dispre
     type(c_ptr) :: m_disseq ! type(disprm) :: m_disseq
  end type linprm

  type, bind(c) :: prjprm
     integer(kind=C_INT) :: flag
     character(kind=C_CHAR), dimension(4) :: code
     real(kind=C_DOUBLE) :: r0
     real(kind=C_DOUBLE), dimension(PVN) :: pv
     real(kind=C_DOUBLE) :: phi0
     real(kind=C_DOUBLE) :: theta0
     integer(kind=C_INT) :: bounds
     character(kind=C_CHAR), dimension(40) :: name
     integer(kind=C_INT) :: category
     integer(kind=C_INT) :: pvrange
     integer(kind=C_INT) :: simplezen
     integer(kind=C_INT) :: equiareal
     integer(kind=C_INT) :: conformal
     integer(kind=C_INT) :: global
     integer(kind=C_INT) :: divergent
     real(kind=C_DOUBLE) :: x0
     real(kind=C_DOUBLE) :: y0
     type(c_ptr) :: err  !type(wcserr) :: err
     type(c_ptr) :: padding
     real(kind=C_DOUBLE), dimension(10) :: w
     integer(kind=C_INT) :: m
     integer(kind=C_INT) :: n
     type(c_funptr) :: prjx2s
     type(c_funptr) :: prjs2x
  end type prjprm

  type, bind(c) :: celprm
     integer(kind=C_INT) :: flag
     integer(kind=C_INT) :: offset
     real(kind=C_DOUBLE) :: phi0
     real(kind=C_DOUBLE) :: theta0
     real(kind=C_DOUBLE), dimension(4) :: ref
     type(prjprm) :: prj
     real(kind=C_DOUBLE), dimension(5) :: euler
     integer(kind=C_INT) :: latpreq
     integer(kind=C_INT) :: isolat
     type(c_ptr) :: err    ! type(wcserr), dimension(*) :: err
     type(c_ptr) :: padding
  end type celprm

  type, bind(c) :: spcprm
     integer(kind=C_INT) :: flag
     character(kind=C_CHAR), dimension(8) :: type
     character(kind=C_CHAR), dimension(4) :: code
     real(kind=C_DOUBLE) :: crval
     real(kind=C_DOUBLE) :: restfrq
     real(kind=C_DOUBLE) :: restwav
     real(kind=C_DOUBLE), dimension(7) :: pv
     real(kind=C_DOUBLE), dimension(6) :: w
     integer(kind=C_INT) :: isGrism
     integer(kind=C_INT) :: padding1
     type(c_ptr) :: err    ! type(wcserr), dimension(*) :: err
     type(c_ptr) :: padding2
     type(c_funptr) :: spxX2P
     type(c_funptr) :: spxP2S
     type(c_funptr) :: spxS2P
     type(c_funptr) :: spxP2X
  end type spcprm

  type, bind(c) :: wcsprm
     integer(kind=C_INT) :: flag
     integer(kind=C_INT) :: naxis
     type(c_ptr) :: crpix !  real(kind=C_DOUBLE), dimension(*) :: crpix
     type(c_ptr) :: pc    ! real(kind=C_DOUBLE), dimension(*) :: pc
     type(c_ptr) :: cdelt ! real(kind=C_DOUBLE), dimension(*) :: cdelt
     type(c_ptr) :: crval ! real(kind=C_DOUBLE), dimension(*) :: crval
     type(c_ptr) :: cunit ! character(kind=C_CHAR), dimension(*) :: cunit
     type(c_ptr) :: ctype ! character(kind=C_CHAR), dimension(*) :: ctype
     real(kind=C_DOUBLE) :: longpole
     real(kind=C_DOUBLE) :: latpole
     real(kind=C_DOUBLE) :: restfrq
     real(kind=C_DOUBLE) :: restwav
     integer(kind=C_INT) :: npv
     integer(kind=C_INT) :: npvmax
     type(c_ptr) :: pv ! type(pvcard), dimension(:) :: pv
     integer(kind=C_INT) :: nps
     integer(kind=C_INT) :: npsmax
     type(c_ptr) :: ps ! type(pscard), dimension(:) :: pv
     type(c_ptr) :: cd ! real(kind=C_DOUBLE), dimension(*) :: cd
     type(c_ptr) :: crota ! real(kind=C_DOUBLE), dimension(*) :: crota
     integer(kind=C_INT) :: altlin
     integer(kind=C_INT) :: velref
     character(kind=C_CHAR), dimension(4) :: alt
     integer(kind=C_INT) :: colnum
     type(c_ptr) :: colax ! integer(kind=C_INT), dimension(:) :: colax
     type(c_ptr) :: cname ! character(kind=C_CHAR), dimension(*) :: cname
     type(c_ptr) :: crder ! real(kind=C_DOUBLE), dimension(*) :: crder
     type(c_ptr) :: csyer ! real(kind=C_DOUBLE), dimension(*) :: csyer
     type(c_ptr) :: czphs ! real(kind=C_DOUBLE), dimension(*) :: czphs
     type(c_ptr) :: cperi ! real(kind=C_DOUBLE), dimension(*) :: cperi
     character(kind=C_CHAR), dimension(72) :: wcsname
     character(kind=C_CHAR), dimension(72) :: timesys
     character(kind=C_CHAR), dimension(72) :: trefpos
     character(kind=C_CHAR), dimension(72) :: trefdir
     character(kind=C_CHAR), dimension(72) :: plephem
     character(kind=C_CHAR), dimension(72) :: timeunit
     character(kind=C_CHAR), dimension(72) :: dateref
     real(kind=C_DOUBLE), dimension(2) :: mjdref
     real(kind=C_DOUBLE) :: timeoffs
     character(kind=C_CHAR), dimension(72) :: dateobs
     character(kind=C_CHAR), dimension(72) :: datebeg
     character(kind=C_CHAR), dimension(72) :: dateavg
     character(kind=C_CHAR), dimension(72) :: dateend
     real(kind=C_DOUBLE) :: mjdobs
     real(kind=C_DOUBLE) :: mjdbeg
     real(kind=C_DOUBLE) :: mjdavg
     real(kind=C_DOUBLE) :: mjdend
     real(kind=C_DOUBLE) :: jepoch
     real(kind=C_DOUBLE) :: bepoch
     real(kind=C_DOUBLE) :: tstart
     real(kind=C_DOUBLE) :: tstop
     real(kind=C_DOUBLE) :: xposure
     real(kind=C_DOUBLE) :: telapse
     real(kind=C_DOUBLE) :: timsyer
     real(kind=C_DOUBLE) :: timrder
     real(kind=C_DOUBLE) :: timedel
     real(kind=C_DOUBLE) :: timepixr
     real(kind=C_DOUBLE), dimension(6) :: obsgeo
     character(kind=C_CHAR), dimension(72) :: obsorbit
     character(kind=C_CHAR), dimension(72) :: radesys
     real(kind=C_DOUBLE) :: equinox
     character(kind=C_CHAR), dimension(72) :: specsys
     character(kind=C_CHAR), dimension(72) :: ssysobs
     real(kind=C_DOUBLE) :: velosys
     real(kind=C_DOUBLE) :: zsource
     character(kind=C_CHAR), dimension(72) :: ssysrc
     real(kind=C_DOUBLE) :: velangl
     type(c_ptr) :: aux  ! type(auxprm), dimension(*) :: aux
     integer(kind=C_INT) :: ntab
     integer(kind=C_INT) :: nwtb
     type(c_ptr) :: tab  ! type(tabprm), dimension(*) :: tab
     type(c_ptr) :: wtb  ! type(wtbarr), dimension(*) :: wtb
     character(kind=C_CHAR), dimension(8) :: lngtyp
     character(kind=C_CHAR), dimension(8) :: lattyp
     integer(kind=C_INT) :: lng
     integer(kind=C_INT) :: lat
     integer(kind=C_INT) :: spec
     integer(kind=C_INT) :: time
     integer(kind=C_INT) :: cubeface
     integer(kind=C_INT) :: chksum
     type(c_ptr) :: types   ! integer(kind=C_INT), dimension(:) :: types
     type(linprm) :: lin
     type(celprm) :: cel
     type(spcprm) :: spc
     type(c_ptr) :: err   ! type(wcserr), dimension(*) :: err
     integer(kind=C_INT) :: m_flag
     integer(kind=C_INT) :: m_naxis
     type(c_ptr) :: m_crpix  ! integer(kind=C_DOUBLE), dimension(*) :: m_crpix
     type(c_ptr) :: m_pc     ! real(kind=C_DOUBLE), dimension(*) :: m_pc
     type(c_ptr) :: m_cdelt  ! real(kind=C_DOUBLE), dimension(*) :: m_cdelt
     type(c_ptr) :: m_crval  ! real(kind=C_DOUBLE), dimension(*) :: m_crval
     character(kind=C_CHAR), dimension(72) :: m_cunit
     character(kind=C_CHAR), dimension(72) :: m_ctype
     type(c_ptr) :: m_pv  ! type(pvcard), dimension(:) :: m_pv
     type(c_ptr) :: m_ps  ! type(pvcard), dimension(:) :: m_ps
     type(c_ptr) :: m_cd  ! real(kind=C_DOUBLE), dimension(:) :: m_cd
     type(c_ptr) :: m_croat ! real(kind=C_DOUBLE), dimension(:) :: m_croat
     type(c_ptr) :: m_colax ! integer(kind=C_INT), dimension(*) :: m_colax
     character(kind=C_CHAR), dimension(72) :: m_cname
     type(c_ptr) :: m_crder ! real(kind=C_DOUBLE), dimension(:) :: m_crder
     type(c_ptr) :: m_csyer ! real(kind=C_DOUBLE), dimension(:) :: m_csyer
     type(c_ptr) :: m_czphs ! real(kind=C_DOUBLE), dimension(:) :: m_czphs
     type(c_ptr) :: m_cperi ! real(kind=C_DOUBLE), dimension(:) :: m_cperi
     type(c_ptr) :: m_aux   ! type(auxprm), dimension(*) :: m_aux
     type(c_ptr) :: m_tab   ! type(tabprm), dimension(*) :: m_tab
     type(c_ptr) :: m_wtb   ! type(wtbarr), dimension(*) :: m_wtb
  end type wcsprm

  type :: prjpar
     type(c_ptr) :: ptr = C_NULL_PTR
     type(prjprm), pointer :: prj => null()
  end type prjpar

  type :: celpar
     type(c_ptr) :: ptr = C_NULL_PTR
     type(celprm), pointer :: cel => null()
  end type celpar

  type :: wcspar
     type(c_ptr) :: ptr = C_NULL_PTR
     type(wcsprm), pointer :: wcs => null()
  end type wcspar

  logical, private :: debug = .false.


  interface

     integer(kind=C_INT) function wcsprt(wcs) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: wcs
     end function wcsprt

     integer(kind=C_INT) function wcsprintf_set(ptr) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: ptr
     end function wcsprintf_set

     integer(kind=C_INT) function wcsperr(wcs,prefix) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: wcs
       character(kind=C_CHAR), dimension(*), intent(in) :: prefix
     end function wcsperr

     subroutine wcsdealloc(wcs) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: wcs
     end subroutine wcsdealloc

     integer(kind=C_INT) function wcspih(header, nkeyrec, relax, crtl, nreject, &
          nwcs, wcs) bind(c)
       use, intrinsic :: iso_c_binding
       character(kind=C_CHAR), dimension(*), intent(in) :: header
       integer(kind=C_INT), value, intent(in) :: nkeyrec, relax, crtl
       integer(kind=C_INT), intent(out) :: nreject, nwcs
       type(c_ptr), intent(out) :: wcs
     end function wcspih

     integer(kind=C_INT) function wcsidx(nwcs, wcs, alts) bind(c)
       use, intrinsic :: iso_c_binding
       integer(kind=C_INT), value, intent(in) :: nwcs
!       type(c_ptr), value, intent(in) :: wcs
       type(c_ptr), intent(in out) :: wcs
       integer(kind=C_INT), dimension(*), intent(out) :: alts
     end function wcsidx

     integer(kind=C_INT) function wcsset(wcs) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: wcs
     end function wcsset

     integer(kind=C_INT) function wcsfix(crtl, naxis, wcs, stat) bind(c)
       use, intrinsic :: iso_c_binding
       integer(kind=C_INT), value, intent(in) :: crtl
       integer(kind=C_LONG), dimension(*), intent(in) :: naxis
       type(c_ptr), value :: wcs
       integer(kind=C_INT), dimension(*), intent(out) :: stat
     end function wcsfix

     integer(kind=C_INT) function wcsini(alloc,naxis,wcs) bind(c)
       use, intrinsic :: iso_c_binding
       integer(kind=C_INT), value, intent(in) :: alloc,naxis
       type(c_ptr), value :: wcs
     end function wcsini

     integer(kind=C_INT) function wcsfree(wcs) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: wcs
     end function wcsfree

     integer(kind=C_INT) function wcshdo(ctrl,wcs,nkeyrec,header) bind(c)
       use, intrinsic :: iso_c_binding
       integer(kind=C_INT), value, intent(in) :: ctrl
       type(c_ptr), value, intent(in) :: wcs
       integer(kind=C_INT), intent(out) :: nkeyrec
       type(c_ptr), intent(out) :: header
     end function wcshdo

     integer(kind=C_INT) function wcsp2s(wcs,ncoord,nelem,pixcrd, &
          imgcrd, phi, theta, world, stat) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: wcs
       integer(kind=C_INT), value, intent(in) :: ncoord, nelem
       real(kind=C_DOUBLE), dimension(*), intent(in) :: pixcrd
       real(kind=C_DOUBLE), dimension(*), intent(out) :: imgcrd, phi, theta, world
       integer(kind=C_INT), dimension(*), intent(out) :: stat
     end function wcsp2s

     integer(kind=C_INT) function wcss2p(wcs,ncoord,nelem,world,phi,theta, &
          imgcrd,pixcrd,stat) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: wcs
       integer(kind=C_INT), value, intent(in) :: ncoord,nelem
       real(kind=C_DOUBLE), dimension(*), intent(in) :: world
       real(kind=C_DOUBLE), dimension(*), intent(out) :: imgcrd, phi, theta, pixcrd
       integer(kind=C_INT), dimension(*), intent(out) :: stat
     end function wcss2p

     integer(kind=C_INT) function wcsput(wcs,what,val,i,j) bind(c,name='wcsput_')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: wcs, val
       integer(kind=C_INT) :: what,i,j
     end function wcsput

     integer(kind=C_INT) function wcsget(wcs,what,val) bind(c,name='wcsget_')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: wcs
       integer(kind=C_INT) :: what
       type(c_ptr), value :: val
     end function wcsget

     integer(kind=C_INT) function disini(deref,naxis,dis) bind(c,name='disini_')
       use, intrinsic :: iso_c_binding
       integer(kind=C_INT), intent(in) :: deref,naxis
       type(c_ptr), value :: dis
     end function disini

     integer(kind=C_INT) function disinit(deref,naxis,dis,ndp) bind(c,name='disinit_')
       use, intrinsic :: iso_c_binding
       integer(kind=C_INT), intent(in) :: deref,naxis,ndp
       type(c_ptr), value :: dis
     end function disinit

     integer(kind=C_INT) function disset(deref,wcs) bind(c,name='disset_')
       use, intrinsic :: iso_c_binding
       integer(kind=C_INT), intent(in) :: deref
       type(c_ptr), value :: wcs
     end function disset

     integer(kind=C_INT) function disput(deref,dis,what,val,i,j) bind(c,name='disput_')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: dis, val
       integer(kind=C_INT) :: deref,what,i,j
     end function disput

     integer(kind=C_INT) function disget(deref,dis,what,val) bind(c,name='disget_')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: dis
       integer(kind=C_INT) :: deref,what
       type(c_ptr), value :: val
     end function disget

     integer(kind=C_INT) function dpfill(dp,keyword,field,j,type,ival,fval) bind(c,name='dpfill_')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: dp, keyword,field
       integer(kind=C_INT) :: j,type,ival
       real(kind=C_DOUBLE) :: fval
     end function dpfill

     integer(kind=C_INT) function linget(lin,what,val) bind(c,name='linget_')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: lin
       integer(kind=C_INT) :: what
       type(c_ptr), value :: val
     end function linget

     integer(kind=C_INT) function prjput(prj,what,val,k) bind(c,name='prjput_')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: prj, val
       integer(kind=C_INT) :: what,k
     end function prjput

     integer(kind=C_INT) function prjget(prj,what,val) bind(c,name='prjget_')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: prj
       integer(kind=C_INT) :: what
       type(c_ptr), value :: val
     end function prjget

     integer(kind=C_INT) function prjini(prj) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: prj
     end function prjini

     integer(kind=C_INT) function prjfree(prj) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: prj
     end function prjfree

     integer(kind=C_INT) function prjperr(prj,prefix) bind(c,name='prjperr_')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: prj
       character(kind=C_CHAR), dimension(*), intent(in) :: prefix
     end function prjperr

     integer(kind=C_INT) function prjset(prj) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: prj
     end function prjset

     integer(kind=C_INT) function prjs2x(prj,npi,ntheta,spt,sxy,phi,theta,x,y,stat) &
          bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: prj
       integer(kind=C_INT), value, intent(in) :: npi,ntheta,spt,sxy
       real(kind=C_DOUBLE), dimension(*), intent(in) :: phi,theta
       real(kind=C_DOUBLE), dimension(*), intent(out) :: x,y
       integer(kind=C_INT), dimension(*), intent(out) :: stat
     end function prjs2x

     integer(kind=C_INT) function wcserr_enable(enable) bind(c)
       use, intrinsic :: iso_c_binding
       integer(kind=C_INT), value, intent(in) :: enable
     end function wcserr_enable

     integer(kind=C_INT) function celini(cel) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: cel
     end function celini

     integer(kind=C_INT) function celfree(cel) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: cel
     end function celfree

     integer(kind=C_INT) function celperr(cel,prefix) bind(c,name='celperr_')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: cel
       character(kind=C_CHAR), dimension(*), intent(in) :: prefix
     end function celperr

     integer(kind=C_INT) function celset(cel) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: cel
     end function celset

     integer(kind=C_INT) function celput(cel,what,val,k) bind(c,name='celput_')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: cel, val
       integer(kind=C_INT) :: what,k
     end function celput

     integer(kind=C_INT) function celget(cel,what,val) bind(c,name='celget_')
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: cel
       integer(kind=C_INT) :: what
       type(c_ptr), value :: val
     end function celget

     integer(kind=C_INT) function cels2x(cel,npi,ntheta,spt,sxy,phi,theta,x,y,stat) &
          bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value, intent(in) :: cel
       integer(kind=C_INT), value, intent(in) :: npi,ntheta,spt,sxy
       real(kind=C_DOUBLE), dimension(*), intent(in) :: phi,theta
       real(kind=C_DOUBLE), dimension(*), intent(out) :: x,y
       integer(kind=C_INT), dimension(*), intent(out) :: stat
     end function cels2x

     ! Both functions are standard libc's ones. They are required for wcsprm.
     ! Ignore gfortran's warning for shadowing of their intrinsics counterparts.
     type(c_ptr) function malloc(size) bind(c)
       use, intrinsic :: iso_c_binding
       integer(c_size_t), value :: size
     end function malloc
     subroutine free(ptr) bind(c)
       use, intrinsic :: iso_c_binding
       type(c_ptr), value :: ptr
     end subroutine free

  end interface

  private :: fchar, cchar

contains

  subroutine wcs_prt(wcs)
    type(wcspar), intent(in) :: wcs
    integer(kind=C_INT) :: status
    status = wcsprt(wcs % ptr)  ! yes, the return value is ignored.
  end subroutine wcs_prt

  subroutine wcs_perr(wcs)
    type(wcspar), intent(in) :: wcs
    character(kind=C_CHAR), dimension(88) :: prf
    integer(kind=C_INT) :: status

    call cchar("Error: ",prf)
    status = wcsperr(wcs % ptr,prf)

  end subroutine wcs_perr

  integer function wcs_image_header(header, nkeyrec, relax, ctrl, nreject, nwcs, wcs)

    character(len=*), intent(in) :: header
    integer, intent(in) :: nkeyrec, relax, ctrl
    integer, intent(out) :: nreject, nwcs
    type(wcspar), intent(out) :: wcs

    character(kind=C_CHAR), dimension(len(header)+1) :: hdr
    integer(C_INT) :: rlx, ctr, nr, nw, nh
    integer :: status

    call cchar(header,hdr)
    rlx = relax
    ctr = ctrl
    nh = nkeyrec
    status = wcspih(hdr,nh,rlx,ctr,nr,nw,wcs % ptr)
    if( status == 0 ) then
       if( nw /= 1 ) write(error_unit,*) 'Warning in wcs_image_header(): ', &
            'The number of coordinate representations found is ',nw, &
            ', while processing only of a single one is implemented.'
       nreject = nr
       nwcs = nw
    end if
    wcs_image_header = status

  end function wcs_image_header

  integer function wcs_index(nwcs, wcs, alts)

    integer, intent(in) :: nwcs
    type(wcspar), intent(in out), target :: wcs
    integer, dimension(:), intent(out) :: alts

    integer(kind=C_INT) :: nw
    integer(kind=C_INT), dimension(size(alts)) :: als

    nw = nwcs
    wcs_index = wcsidx(nw, wcs % ptr, als)
    alts = als

  end function wcs_index

  integer function wcs_set(wcs)

    type(wcspar), intent(in out) :: wcs

    wcs_set = wcsset(wcs % ptr)

  end function wcs_set


  function wcs_header_output(ctrl,wcs,nkeyrec,header) result(status)

    integer, parameter :: LEN_CARD = 80

    integer, intent(in) :: ctrl
    type(wcspar), intent(in) :: wcs
    integer, intent(out) :: nkeyrec
    character(len=:), allocatable, intent(out) :: header
    type(c_ptr) :: phdr = C_NULL_PTR
    character, dimension(:), pointer :: hdr
    integer(kind=C_INT) :: ctr, nk
    integer :: status

    ctr = ctrl
    status = wcshdo(ctr,wcs%ptr,nk,phdr)
    if( status == 0 ) then
       nkeyrec = nk
       call c_f_pointer(phdr,hdr,[nkeyrec*LEN_CARD+1])
       allocate(character(len=nkeyrec*LEN_CARD+1) :: header)
       call fchar(hdr,header)
    end if
    call wcsdealloc(phdr)

  end function wcs_header_output

  integer function wcs_pixel2world(wcs,pixcrd,imgcrd,phi,theta,world,stat)

    type(wcspar), intent(in) :: wcs
    real(REAL64), dimension(:,:), intent(in) :: pixcrd
    real(REAL64), dimension(:,:), intent(out) :: imgcrd, world
    real(REAL64), dimension(:), intent(out) :: phi, theta
    integer, dimension(:), intent(out) :: stat

    real(kind=C_DOUBLE), dimension(size(pixcrd)) :: pix
    real(kind=C_DOUBLE), dimension(size(imgcrd)) :: img
    real(kind=C_DOUBLE), dimension(size(phi)) :: fi
    real(kind=C_DOUBLE), dimension(size(theta)) :: the
    real(kind=C_DOUBLE), dimension(size(world)) :: word
    integer(kind=C_INT), dimension(size(stat)) :: st
    integer(kind=C_INT) :: ncoord, nelem
    integer :: i,j, j1, j2

    ncoord = size(pixcrd,1)
    nelem = size(pixcrd,2)

    do i = 1, ncoord
       j = (i-1)*nelem
       j1 = j+1
       j2 = j+nelem
       pix(j1:j2) = pixcrd(i,:)
    end do

    wcs_pixel2world = wcsp2s(wcs % ptr,ncoord,nelem,pix,img,fi,the,word,st)

    if( wcs_pixel2world == 0 ) then
       do i = 1, ncoord
          j = (i-1)*nelem
          j1 = j+1
          j2 = j+nelem
          world(i,:) = word(j1:j2)
          imgcrd(i,:) = img(j1:j2)
       end do
       phi = fi
       theta = the
       stat = st
    end if

  end function wcs_pixel2world

  integer function wcs_world2pixel(wcs,world,phi,theta,imgcrd,pixcrd,stat)

    type(wcspar), intent(in) :: wcs
    real(REAL64), dimension(:,:), intent(in) :: world
    real(REAL64), dimension(:,:), intent(out) :: imgcrd, pixcrd
    real(REAL64), dimension(:), intent(out) :: phi, theta
    integer, dimension(:), intent(out) :: stat

    real(kind=C_DOUBLE), dimension(size(pixcrd)) :: pix
    real(kind=C_DOUBLE), dimension(size(imgcrd)) :: img
    real(kind=C_DOUBLE), dimension(size(phi)) :: fi
    real(kind=C_DOUBLE), dimension(size(theta)) :: the
    real(kind=C_DOUBLE), dimension(size(world)) :: word
    integer(kind=C_INT), dimension(size(stat)) :: st
    integer(kind=C_INT) :: ncoord, nelem
    integer :: i,j, j1, j2

    ncoord = size(world,1)
    nelem = size(world,2)

    do i = 1, ncoord
       j = (i-1)*nelem
       j1 = j+1
       j2 = j+nelem
       word(j1:j2) = world(i,:)
    end do

    wcs_world2pixel = wcss2p(wcs % ptr,ncoord,nelem,word,fi,the,img,pix,st)

    do i = 1, ncoord
       j = (i-1)*nelem
       j1 = j+1
       j2 = j+nelem
       pixcrd(i,:) = pix(j1:j2)
       imgcrd(i,:) = img(j1:j2)
    end do
    phi = fi
    theta = the
    stat = st

  end function wcs_world2pixel



  pure subroutine fchar(cstring,fstring)

    character(kind=C_CHAR), dimension(*), intent(in) :: cstring
    character(len=*), intent(out) :: fstring

    integer :: i,n

    n = 1
    do i = 1, len(fstring)
       n = i
       if( cstring(n) == C_NULL_CHAR ) exit
       fstring(n:n) = cstring(n)
    end do
    fstring(n:) = ' '

  end subroutine fchar

  pure subroutine cchar(fstring,cstring)

    character(len=*), intent(in) :: fstring
    character(kind=C_CHAR), dimension(*), intent(out) :: cstring
    integer :: i, n

    n = len_trim(fstring)
    do i = 1, n
       cstring(i) = fstring(i:i)
    end do
    cstring(n+1) = C_NULL_CHAR

  end subroutine cchar

  subroutine wcs_free(wcs)
    type(wcspar) :: wcs
    integer :: status

    status = wcsfree(wcs % ptr)
    ! the status can be safety ignored as it is non-zero only for null pointer

    ! This function is the libc's free().
    call free(wcs % ptr)

  end subroutine wcs_free

  function wcs_alloc()  result(wcs)
    type(wcspar) :: wcs

    ! The function is stolen from the standard libc.
    wcs % ptr = malloc(c_sizeof(wcs % wcs))
    call c_f_pointer(wcs % ptr, wcs % wcs)
    if( .not. c_associated(wcs % ptr) ) stop 'A memory allocation failed.'

  end function wcs_alloc

  function wcs_ini(naxis,wcs) result(status)
    type(wcspar), intent(in out) :: wcs
    integer, intent(in) :: naxis
    integer :: status
    integer(kind=C_INT) :: n
    integer(kind=C_INT), target :: m

    n = naxis
    m = -1
    status = wcsput(wcs % ptr, WCS_FLAG, C_LOC(m), 0, 0)
    if( status == 0 ) then
       status = wcsini(1,n, wcs % ptr)
    else
       stop 'An attempt to initialise of WCSLIB failed.'
    end if

  end function wcs_ini

  integer function wcs_fix(ctrl,naxis,wcs,stat)
    integer, intent(in) :: ctrl
    integer, dimension(:), intent(in) :: naxis
    type(wcspar), intent(in out) :: wcs
    integer, dimension(:), intent(out) :: stat
    integer(kind=C_INT) :: c
    integer(kind=C_INT), dimension(size(stat)) :: st
    integer(kind=C_LONG), dimension(size(naxis)) :: na

    c = ctrl
    na = naxis
    wcs_fix = wcsfix(c,na,wcs % ptr, st)
    stat = st

  end function wcs_fix

  function wcs_put_int(wcs,what,val,i,j) result(status)
    type(wcspar), intent(in out) :: wcs
    integer, intent(in) :: val
    integer(kind=C_INT) :: what,i,j
    integer :: status
    integer(kind=C_INT), target :: v
    integer(kind=C_INT) :: w,k,l

    v = val
    w = what
    k = i
    l = j
    status = wcsput(wcs % ptr,w,C_LOC(v),k,l)

  end function wcs_put_int

  function wcs_put_double(wcs,what,val,i,j) result(status)
    type(wcspar), intent(in out) :: wcs
    real(REAL64), intent(in) :: val
    integer(kind=C_INT) :: what,i,j
    integer :: status
    real(kind=C_DOUBLE), target :: v
    integer(kind=C_INT) :: w,k,l

    v = val
    w = what
    k = i
    l = j
    status = wcsput(wcs % ptr,w,C_LOC(v),k,l)

  end function wcs_put_double

  function wcs_put_char(wcs,what,val,i,j) result(status)
    type(wcspar), intent(in out) :: wcs
    character(len=*), intent(in) :: val
    integer(kind=C_INT) :: what,i,j
    integer :: status
    character(kind=C_CHAR), dimension(len(val)+1), target :: a
    integer(kind=C_INT) :: w,k,l

    w = what
    call cchar(val,a)
    k = i
    l = j
    status = wcsput(wcs % ptr,w,C_LOC(a),k,l)

  end function wcs_put_char

  function wcs_get_char(wcs,what,val) result(status)
    type(wcspar), intent(in out) :: wcs
    character(len=*), dimension(:), intent(out) :: val
    integer, intent(in) :: what

    integer :: LEN_VALUE = 72
    integer :: status
    character(len=len(val)) :: b
    character(kind=C_CHAR), dimension(:), pointer :: a
    integer(kind=C_INT) :: w
    type(c_ptr) :: ptr = C_NULL_PTR
    integer :: i,n,j,l1,l2,k

    select case(what)
    case(WCS_ALT, WCS_LNGTYP, WCS_LATTYP)
       LEN_VALUE = 4
    case default
       LEN_VALUE = 72
    end select

    n = size(val)
    allocate(a(n*LEN_VALUE))
    ptr = C_LOC(a)
    w = what
    status = wcsget(wcs % ptr,w,ptr)
    if( status == 0 ) then
       call c_f_pointer(ptr,a,[n*LEN_VALUE])
       do i = 1, n
          j = (i-1)*LEN_VALUE
          l1 = j+1
          l2 = j + min(LEN_VALUE,len(b))
          do j = l1, l2
             if( a(j) == C_NULL_CHAR ) exit
             k = j-l1+1
             b(k:k) = a(j)
          end do
          b(k+1:) = ' '
          val(i) = b
       end do
    else
       val = ''
    end if
    deallocate(a)

  end function wcs_get_char

  function wcs_get_char1(wcs,what,val) result(status)
    type(wcspar), intent(in out) :: wcs
    character(len=*), intent(out) :: val
    integer, intent(in) :: what

    integer :: LEN_VALUE = 72
    integer :: status
    character(len=len(val)) :: b
    character(kind=C_CHAR), dimension(:), pointer :: a
    integer(kind=C_INT) :: w
    type(c_ptr) :: ptr = C_NULL_PTR
    integer :: j

    select case(what)
    case(WCS_ALT, WCS_LNGTYP, WCS_LATTYP)
       LEN_VALUE = 4
    case default
       LEN_VALUE = 72
    end select

    allocate(a(LEN_VALUE))
    ptr = C_LOC(a)
    w = what
    status = wcsget(wcs % ptr,w,ptr)
    if( status == 0 ) then
       call c_f_pointer(ptr,a,[LEN_VALUE])
       do j = 1, min(LEN_VALUE,len(b))
          if( a(j) == C_NULL_CHAR ) exit
          b(j:j) = a(j)
       end do
       b(j+1:) = ' '
       val = b
    else
       val = ''
    end if
    deallocate(a)

  end function wcs_get_char1

  function wcs_get_int(wcs,what,val) result(status)
    type(wcspar), intent(in out) :: wcs
    integer, intent(out) :: val
    integer, intent(in) :: what
    integer :: status
    integer(kind=C_INT), dimension(:), pointer :: d
    integer(kind=C_INT) :: w
    type(c_ptr) :: ptr = C_NULL_PTR

    allocate(d(1))
    ptr = C_LOC(d)
    w = what
    status = wcsget(wcs % ptr,w,ptr)
    if( status == 0 ) then
       call c_f_pointer(ptr,d,[1])
       val = d(1)
    else
       val = 0
    end if
    deallocate(d)

  end function wcs_get_int

  function wcs_get_double1(wcs,what,val) result(status)
    type(wcspar), intent(in out) :: wcs
    real(REAL64), intent(out) :: val
    integer, intent(in) :: what
    integer :: status
    real(kind=C_DOUBLE), dimension(:), pointer :: d
    integer(kind=C_INT) :: w
    type(c_ptr) :: ptr = C_NULL_PTR

    allocate(d(1))
    ptr = C_LOC(d)
    w = what
    status = wcsget(wcs % ptr,w,ptr)
    if( status == 0 ) then
       call c_f_pointer(ptr,d,[1])
       val = d(1)
    else
       val = 0
    end if
    deallocate(d)

  end function wcs_get_double1

  function wcs_get_double(wcs,what,val) result(status)
    type(wcspar), intent(in out) :: wcs
    real(REAL64), dimension(:), intent(out) :: val
    integer, intent(in) :: what
    integer :: status
    real(kind=C_DOUBLE), dimension(:), pointer :: d
    integer(kind=C_INT) :: w
    type(c_ptr) :: ptr = C_NULL_PTR
    integer :: i,n

    n = size(val)
    allocate(d(n))
    ptr = C_LOC(d)
    w = what
    status = wcsget(wcs % ptr,w,ptr)
    if( status == 0 ) then
       call c_f_pointer(ptr,d,[n])
       do i = 1, n
          val(i) = d(i)
       end do
    else
       val = 0
    end if
    deallocate(d)

  end function wcs_get_double

  function wcs_get_double2d(wcs,what,val) result(status)
    type(wcspar), intent(in out) :: wcs
    real(REAL64), dimension(:,:), intent(out) :: val
    integer, intent(in) :: what
    integer :: status
    real(kind=C_DOUBLE), dimension(:), pointer :: d
    integer(kind=C_INT) :: w
    type(c_ptr) :: ptr = C_NULL_PTR
    integer :: i,j,n,m

    n = size(val,1)
    m = size(val,2)
    allocate(d(n*m))
    ptr = C_LOC(d)
    w = what
    status = wcsget(wcs % ptr,w,ptr)
    if( status == 0 ) then
       call c_f_pointer(ptr,d,[n*m])
       do i = 1, n
          j = (i-1)*m
          val(i,:) = d(j+1:j+m)
       end do
    else
       val = 0
    end if
    deallocate(d)

  end function wcs_get_double2d

  !- prj

  integer function prj_init(prj)
    type(prjpar), intent(out) :: prj
!    type(prjprm), pointer :: p

    allocate(prj % prj)
    prj % ptr = C_LOC(prj % prj)
    prj_init = prjini(prj % ptr)

!    prj_init = palloc(prj % ptr)
!    prj % prj =
!    prj % ptr = C_LOC(prj % prj)

  end function prj_init

  integer function prj_free(prj)
    type(prjpar), intent(in out) :: prj

    if( debug ) &
         write(error_unit,*) c_associated(prj % ptr),associated(prj % prj)
    prj_free = prjfree(prj % ptr)
    deallocate(prj % prj)
!    prj % prj => null()
!    prj % ptr = C_NULL_PTR

  end function prj_free

  integer function prj_perr(prj,prefix)
    type(prjpar), intent(in) :: prj
    character(len=*), intent(in) :: prefix
    character(kind=C_CHAR), dimension(len(prefix)+1) :: prf

    call cchar(prefix,prf)
    prj_perr = prjperr(prj % ptr,prf)

  end function prj_perr

  integer function prj_set(prj)
    type(prjpar), intent(in out) :: prj

    prj_set = prjset(prj % ptr)

  end function prj_set

  integer function prj_cel2cart(prj,alpha,delta,u,v,stat)

    type(prjpar), intent(in out) :: prj
    real(REAL64), dimension(:), intent(in) :: alpha,delta
    real(REAL64), dimension(:), intent(out) :: u,v
    integer, dimension(:), intent(out) :: stat

    real(kind=C_DOUBLE), dimension(size(alpha)) :: a, d, x, y
    integer(kind=C_INT), dimension(size(stat)) :: st
    integer(kind=C_INT) :: n
    integer(kind=C_INT), parameter :: is = 1
    integer(kind=C_INT), parameter :: ntheta = 1
!    type(c_ptr) :: xptr,yptr
!    integer :: i

    n = size(alpha)
!    is = 1

!    prj % prj % phi0 = -33
!    prj % prj % code(1:1) = 'q'
!!$    do i = 1, n
!!$       a(i) = alpha(i)
!!$       d(i) = delta(i)
!!$       x(i) = 0
!!$       y(i) = 0
!!$       st(i) = 0
!!$    end do
    a = alpha
    d = delta
!    write(*,*) n,a
!    return
    prj_cel2cart = prjs2x(prj % ptr,n,ntheta,is,is,a,d,x,y,st)
    !    ptr = C_NULL_PTR
!    xptr = c_LOC(x)
!    yptr = c_LOC(y)
    !prj_cel2cart = tans2xx(prj%prj,n,ntheta,is,is,a,d,x,y,st)
!    write(*,*) 'prj_cel2cart',prj_cel2cart
    if( prj_cel2cart == 0 ) then
       u = x
       v = y
       stat = st
!!$       do i = 1, n
!!$          u(i) = x(i)
!!$          v(i) = y(i)
!!$          stat(i) = st(i)
!!$          write(*,*) u(i),v(i),stat(i)
!!$       end do
    end if

  end function prj_cel2cart

  function prj_put_double(prj,what,val,m) result(status)
    type(prjpar), intent(in out) :: prj
    real(REAL64), intent(in) :: val
    integer(kind=C_INT) :: what,m
    integer :: status
    real(kind=C_DOUBLE), target :: v
    integer(kind=C_INT) :: w,k

    v = val
    w = what
    k = m
    status = prjput(prj % ptr,w,C_LOC(v),k)

  end function prj_put_double

  function prj_put_int(prj,what,val,m) result(status)
    type(prjpar), intent(in out) :: prj
    integer, intent(in) :: val
    integer(kind=C_INT) :: what,m
    integer :: status
    integer(kind=C_INT), target :: v
    integer(kind=C_INT) :: w,k
    type(c_ptr) :: ptr

    v = val
    ptr = C_LOC(v)
    w = what
    k = m
    status = prjput(prj % ptr,w,ptr,k)

  end function prj_put_int

  function prj_put_char(prj,what,val,m) result(status)
    type(prjpar), intent(in out) :: prj
    character(len=*), intent(in) :: val
    integer(kind=C_INT) :: what,m
    integer :: status
    character(kind=C_CHAR), dimension(len(val)+1), target :: a
    integer(kind=C_INT) :: w,k

    w = what
    call cchar(val,a)
    k = m
    status = prjput(prj % ptr,w,C_LOC(a),k)

  end function prj_put_char

  function prj_get_char(prj,what,val) result(status)
    type(prjpar), intent(in out) :: prj
    character(len=*), dimension(:), intent(out) :: val
    integer, intent(in) :: what

    integer :: LEN_VALUE = 72
    integer :: status
    character(len=len(val)) :: b
    character(kind=C_CHAR), dimension(:), pointer :: a
    integer(kind=C_INT) :: w
    type(c_ptr) :: ptr = C_NULL_PTR
    integer :: i,n,j,l1,l2,l

    select case(what)
    case(WCS_ALT, WCS_LNGTYP, WCS_LATTYP)
       LEN_VALUE = 4
    case default
       LEN_VALUE = 72
    end select

    n = size(val)
    allocate(a(n*LEN_VALUE))
    ptr = C_LOC(a)
    w = what
!    write(*,*) c_associated(ptr),c_associated(prj%ptr)
    status = prjget(prj % ptr,w,ptr)
    if( status == 0 ) then
       call c_f_pointer(ptr,a,[n*LEN_VALUE])
       do i = 1, n
          j = (i-1)*LEN_VALUE
          l1 = j+1
          l2 = j + min(LEN_VALUE,len(b))
          do j = l1, l2
             if( a(j) == C_NULL_CHAR ) exit
             l = j-l1+1
!             write(*,*) i,l,l1,l2
             b(l:l) = a(j)
          end do
          b(l+1:) = ' '
          val(i) = b
       end do
    else
       val = ''
    end if
    deallocate(a)

  end function prj_get_char

  integer function prj_get_char1(prj,what,val)
    type(prjpar), intent(in out) :: prj
    character(len=*), intent(out) :: val
    integer, intent(in) :: what

    character(len=len(val)), dimension(1) :: a
    prj_get_char1 = prj_get(prj,what,a)
    val = a(1)

  end function prj_get_char1


  function prj_get_double(prj,what,val) result(status)
    type(prjpar), intent(in out) :: prj
    real(REAL64), dimension(:), intent(out) :: val
    integer, intent(in) :: what
    integer :: status
    real(kind=C_DOUBLE), dimension(:), pointer :: d
    integer(kind=C_INT) :: w
    type(c_ptr) :: ptr = C_NULL_PTR
    integer :: i,n

    n = size(val)
    allocate(d(n))
    ptr = C_LOC(d)
    w = what
    status = prjget(prj % ptr,w,ptr)
    if( status == 0 ) then
       call c_f_pointer(ptr,d,[n])
       do i = 1, n
          val(i) = d(i)
       end do
    else
       val = 0
    end if
    deallocate(d)

  end function prj_get_double

  function prj_get_double2d(prj,what,val) result(status)
    type(prjpar), intent(in out) :: prj
    real(REAL64), dimension(:,:), intent(out) :: val
    integer, intent(in) :: what
    integer :: status
    real(kind=C_DOUBLE), dimension(:), pointer :: d
    integer(kind=C_INT) :: w
    type(c_ptr) :: ptr = C_NULL_PTR
    integer :: i,j,n,m

    n = size(val,1)
    m = size(val,2)
    allocate(d(n*m))
    ptr = C_LOC(d)
    w = what
    status = prjget(prj % ptr,w,ptr)
    if( status == 0 ) then
       call c_f_pointer(ptr,d,[n*m])
       do i = 1, n
          j = (i-1)*m
          val(i,:) = d(j+1:j+m)
       end do
    else
       val = 0
    end if
    deallocate(d)

  end function prj_get_double2d

  ! -- cel

  integer function cel_init(cel)
    type(celpar), intent(out) :: cel

    allocate(cel % cel)
    cel % ptr = C_LOC(cel % cel)
    cel_init = celini(cel % ptr)

  end function cel_init

  integer function cel_free(cel)
    type(celpar), intent(in out) :: cel

!    write(*,*) c_associated(cel % ptr),associated(cel % cel)
    cel_free = celfree(cel % ptr)
    deallocate(cel % cel)

  end function cel_free

  integer function cel_perr(cel,prefix)
    type(celpar), intent(in) :: cel
    character(len=*), intent(in) :: prefix
    character(kind=C_CHAR), dimension(len(prefix)+1) :: prf

    call cchar(prefix,prf)
    cel_perr = celperr(cel % ptr,prf)

  end function cel_perr

  integer function cel_set(cel)
    type(celpar), intent(in out) :: cel

    cel_set = celset(cel % ptr)

  end function cel_set

  function cel_put_double(cel,what,val,m) result(status)
    type(celpar), intent(in out) :: cel
    real(REAL64), intent(in) :: val
    integer(kind=C_INT) :: what,m
    integer :: status
    real(kind=C_DOUBLE), target :: v
    integer(kind=C_INT) :: w,k

    v = val
    w = what
    k = m
    status = celput(cel % ptr,w,C_LOC(v),k)

  end function cel_put_double

  function cel_put_int(cel,what,val,m) result(status)
    type(celpar), intent(in out) :: cel
    integer, intent(in) :: val
    integer(kind=C_INT) :: what,m
    integer :: status
    integer(kind=C_INT), target :: v
    integer(kind=C_INT) :: w,k
    type(c_ptr) :: ptr

    v = val
    ptr = C_LOC(v)
    w = what
    k = m
    status = celput(cel % ptr,w,ptr,k)

  end function cel_put_int

  function cel_put_char(cel,what,val,m) result(status)
    type(celpar), intent(in out) :: cel
    character(len=*), intent(in) :: val
    integer(kind=C_INT) :: what,m
    integer :: status
    character(kind=C_CHAR), dimension(len(val)+1), target :: a
    integer(kind=C_INT) :: w,k

    w = what
    call cchar(val,a)
    k = m
    status = celput(cel % ptr,w,C_LOC(a),k)

  end function cel_put_char

  function cel_get_char(cel,what,val) result(status)
    type(celpar), intent(in out) :: cel
    character(len=*), dimension(:), intent(out) :: val
    integer, intent(in) :: what

    integer :: LEN_VALUE = 72
    integer :: status
    character(len=len(val)) :: b
    character(kind=C_CHAR), dimension(:), pointer :: a
    integer(kind=C_INT) :: w
    type(c_ptr) :: ptr = C_NULL_PTR
    integer :: i,n,j,l1,l2,l

    select case(what)
    case(WCS_ALT, WCS_LNGTYP, WCS_LATTYP)
       LEN_VALUE = 4
    case default
       LEN_VALUE = 72
    end select

    n = size(val)
    allocate(a(n*LEN_VALUE))
    ptr = C_LOC(a)
    w = what
!    write(*,*) c_associated(ptr),c_associated(cel%ptr)
    status = celget(cel % ptr,w,ptr)
    if( status == 0 ) then
       call c_f_pointer(ptr,a,[n*LEN_VALUE])
       do i = 1, n
          j = (i-1)*LEN_VALUE
          l1 = j+1
          l2 = j + min(LEN_VALUE,len(b))
          do j = l1, l2
             if( a(j) == C_NULL_CHAR ) exit
             l = j-l1+1
!             write(*,*) i,l,l1,l2
             b(l:l) = a(j)
          end do
          b(l+1:) = ' '
          val(i) = b
       end do
    else
       val = ''
    end if
    deallocate(a)

  end function cel_get_char

  integer function cel_get_char1(cel,what,val)
    type(celpar), intent(in out) :: cel
    character(len=*), intent(out) :: val
    integer, intent(in) :: what

    character(len=len(val)), dimension(1) :: a
    cel_get_char1 = cel_get(cel,what,a)
    val = a(1)

  end function cel_get_char1


  function cel_get_double(cel,what,val) result(status)
    type(celpar), intent(in out) :: cel
    real(REAL64), dimension(:), intent(out) :: val
    integer, intent(in) :: what
    integer :: status
    real(kind=C_DOUBLE), dimension(:), pointer :: d
    integer(kind=C_INT) :: w
    type(c_ptr) :: ptr = C_NULL_PTR
    integer :: i,n

    n = size(val)
    allocate(d(n))
    ptr = C_LOC(d)
    w = what
    status = celget(cel % ptr,w,ptr)
    if( status == 0 ) then
       call c_f_pointer(ptr,d,[n])
       do i = 1, n
          val(i) = d(i)
       end do
    else
       val = 0
    end if
    deallocate(d)

  end function cel_get_double

  function cel_get_double2d(cel,what,val) result(status)
    type(celpar), intent(in out) :: cel
    real(REAL64), dimension(:,:), intent(out) :: val
    integer, intent(in) :: what
    integer :: status
    real(kind=C_DOUBLE), dimension(:), pointer :: d
    integer(kind=C_INT) :: w
    type(c_ptr) :: ptr = C_NULL_PTR
    integer :: i,j,n,m

    n = size(val,1)
    m = size(val,2)
    allocate(d(n*m))
    ptr = C_LOC(d)
    w = what
    status = celget(cel % ptr,w,ptr)
    if( status == 0 ) then
       call c_f_pointer(ptr,d,[n*m])
       do i = 1, n
          j = (i-1)*m
          val(i,:) = d(j+1:j+m)
       end do
    else
       val = 0
    end if
    deallocate(d)

  end function cel_get_double2d

  integer function cel_cel2cart(cel,alpha,delta,u,v,stat)

    type(celpar), intent(in out) :: cel
    real(REAL64), dimension(:), intent(in) :: alpha,delta
    real(REAL64), dimension(:), intent(out) :: u,v
    integer, dimension(:), intent(out) :: stat

    real(kind=C_DOUBLE), dimension(size(alpha)) :: a, d, x, y
    integer(kind=C_INT), dimension(size(stat)) :: st
    integer(kind=C_INT) :: n
    integer(kind=C_INT), parameter :: is = 1
    integer(kind=C_INT), parameter :: ntheta = 1

    n = size(alpha)
    a = alpha
    d = delta
    cel_cel2cart = cels2x(cel % ptr,n,ntheta,is,is,a,d,x,y,st)
    if( cel_cel2cart == 0 ) then
       u = x
       v = y
       stat = st
    end if

  end function cel_cel2cart

  ! -- distortions

  function dis_put_int(dis,what,val,i,j) result(status)
    type(c_ptr), intent(in out) :: dis
    integer, intent(in) :: val,what,i,j
    integer :: status
    integer(kind=C_INT), target :: v
    integer(kind=C_INT) :: w,k,l

    v = val
    w = what
    k = i
    l = j
    status = disput(0,dis,w,C_LOC(v),k,l)

  end function dis_put_int

  function dis_put_char(dis,what,val,i,j) result(status)
    type(c_ptr), intent(in out) :: dis
    character(len=*), intent(in) :: val
    integer, intent(in) :: what,i,j
    integer :: status
    character(kind=C_CHAR), dimension(len(val)+1), target :: a
    integer(kind=C_INT) :: w,k,l

    call cchar(val,a)
    w = what
    k = i
    l = j
    status = disput(0,dis,w,C_LOC(a),k,l)

  end function dis_put_char

  function dis_get_int(dis,what,val) result(status)
    type(c_ptr), intent(in out) :: dis
    integer, intent(out) :: val
    integer, intent(in) :: what
    integer :: status
    integer(kind=C_INT), dimension(:), pointer :: d
    integer(kind=C_INT) :: w
    type(c_ptr) :: ptr = C_NULL_PTR

    allocate(d(1))
    ptr = C_LOC(d)
    w = what
    status = disget(1,dis,w,ptr)
    if( status == 0 ) then
       call c_f_pointer(ptr,d,[1])
       val = d(1)
    else
       val = 0
    end if
    deallocate(d)
  end function dis_get_int

  function dis_get_double1(dis,what,val) result(status)
    type(c_ptr), intent(in out) :: dis
    real(REAL64), intent(out) :: val
    integer, intent(in) :: what
    integer :: status
    real(kind=C_DOUBLE), dimension(:), pointer :: d
    integer(kind=C_INT) :: w
    type(c_ptr) :: ptr = C_NULL_PTR

    allocate(d(1))
    ptr = C_LOC(d)
    w = what
    status = disget(1,dis,w,ptr)
    if( status == 0 ) then
       call c_f_pointer(ptr,d,[1])
       val = d(1)
    else
       val = 0
    end if
    deallocate(d)
  end function dis_get_double1

  function dis_get_pvkey(dis,what,val) result(status)
    type(c_ptr), intent(in out) :: dis
    type(dpkey), dimension(:), intent(out) :: val
    integer, intent(in) :: what
    integer :: status, n, i
!    real(kind=C_DOUBLE), dimension(:), pointer :: d
    type(dpkey), dimension(:), pointer :: d
    integer(kind=C_INT) :: w
    type(c_ptr) :: ptr = C_NULL_PTR
!    integer, parameter :: DPLEN = 22

!    n = DPLEN
    n = size(val)!*DPLEN
    allocate(d(n))
    ptr = C_LOC(d)
    w = what
    status = disget(1,dis,w,ptr)
    if( status == 0 ) then
       call c_f_pointer(ptr,d,[n])
       do i = 1, n
          val(i) = d(i) !% value
       end do
!    else
!       val = 0
    end if
    deallocate(d)
  end function dis_get_pvkey

  function dis_put_pvkey(dis,what,dp) result(status)
    type(c_ptr), intent(in out) :: dis
    type(dpkey), dimension(:), intent(in), target :: dp
    integer, intent(in) :: what
    integer :: status
!    real(kind=C_DOUBLE), dimension(:), pointer :: d
!    type(dpkey), dimension(:), pointer :: d
    integer(kind=C_INT) :: w,j,dummy
!    type(c_ptr) :: ptr = C_NULL_PTR
!    integer, parameter :: DPLEN = 22

    w = what
    j = 0
    dummy = 0
    status = disput(0,dis,w,C_LOC(dp),j,dummy)

  end function dis_put_pvkey

  subroutine wcs_get_dis(wcs,pv)
    type(wcspar), intent(in out) :: wcs
    real(REAL64), dimension(:,:), allocatable, intent(out) :: pv

    type(c_ptr) :: pdis = C_NULL_PTR
    type(c_ptr) :: plin = C_NULL_PTR
    type(linprm), pointer :: lin
    type(disprm), pointer :: dis
!    type(wcsprm), pointer
!    real(REAL64), dimension(:), allocatable :: dp
    type(dpkey), dimension(:), allocatable :: dpk
    integer :: n, status, ndp, naxis, i,j

    write(error_unit,*) c_associated(wcs % ptr),'aaaaaaaaaaaaa'
    allocate(lin)
    plin = C_LOC(lin)
    status = wcsget(wcs % ptr,WCS_LIN,plin)
!    lin => wcs % wcs % lin
    write(error_unit,*) status,'1'
!    if( status /= 0 ) return
!    lin = C_LOC(wcs % wcs % lin)

    allocate(dis)
    pdis = C_LOC(dis)
    status = linget(plin,LIN_DISSEQ,pdis)
!!    status = linget(wcs % wcs % lin,LIN_DISSEQ,pdis)
    write(error_unit,*) status,'2'
    !    if( status /= 0 ) return
!    dis = lin %

!    lin = wcs % wcs % lin
!    pdis = wcs % wcs % lin % disseq
!    pdis = lin % disseq
!    call c_f_pointer(wcs % wcs % lin % disseq,dis)

    status = dis_get(pdis,DIS_NAXIS,naxis)
    write(error_unit,*) status,'3',naxis
    status = dis_get(pdis,DIS_NDP,ndp)
    write(error_unit,*) status,'4',ndp

    allocate(dpk(ndp))
    status = dis_get_pvkey(pdis,DIS_DP,dpk)
    write(error_unit,*) status,'5'

    if( status == 0 ) then
       n = ndp / naxis
       allocate(pv(naxis,n-3))
       do i = 1, naxis
          do j = 1,n - 3
             pv(i,j) = dpk(j+3) % value
          end do
!          j = (i-1)*n
!          write(error_unit,*) j+4,j+n
!          pv(i,:) = dp(j+4:j+n)
       end do
!       pv(2,:) = dp(n+1:)
    end if
!    deallocate(dp)

    write(error_unit,*) 'aaa'
!    plin = C_NULL_PTR
!    pdis = C_NULL_PTR
    deallocate(lin,dis)

  end subroutine wcs_get_dis

  subroutine wcs_put_dis(wcs,dtype,pv)
    type(wcspar), intent(in out) :: wcs
    character(len=*), dimension(:), intent(in) :: dtype
    real(REAL64), dimension(:,:), intent(in) :: pv

    type(c_ptr) :: pdis = C_NULL_PTR
    type(c_ptr) :: plin = C_NULL_PTR
    type(linprm), pointer :: lin
    type(disprm), pointer :: dis
!    real(REAL64), dimension(:), pointer :: dp
    type(dpkey), dimension(:), pointer :: dpk
    character(kind=C_CHAR), dimension(88), target :: ckeyword, cfield
    character(len=88) :: keyword,field
    real(kind=C_DOUBLE) :: fval
    integer(kind=C_INT) :: idx, type, ival, w!, ii, jj
    !   real(kind=C_DOUBLE), pointer :: pval
    integer :: i,j,ndp,naxis,status

    ndp = size(pv,2)
    naxis = size(pv,1)

    w = WCS_PV

    do i = 1, naxis
!       ii = i
       do j = 1, ndp
!          jj = j
          status = wcs_put(wcs,w, pv(i,j), i, j)
       end do
    end do

    return

    allocate(lin)
    plin = C_LOC(lin)
    status = wcsget(wcs % ptr,WCS_LIN,plin)

    allocate(dis)
    pdis = C_LOC(dis)
    status = linget(plin,LIN_DISSEQ,pdis)

    ndp = size(pv,2)
    naxis = size(pv,1)

    status = dis_put(pdis,DIS_FLAG, -1, 0, 0)
    status = disinit(0,naxis,pdis,ndp)

    do j = 1, naxis
       status = dis_put(pdis,DIS_DTYPE,dtype(j), i, 0)
    end do

!    allocate(dpk(ndp))
!    status = dis_get_pvkey(pdis,DIS_DP,dpk)
!    write(error_unit,*) status,'5'

    type =  1
    ival =  0
    allocate(dpk(ndp))
    do i = 1, naxis
       write(keyword,'(a,i0,a)') 'DQ',i,C_NULL_CHAR
       call cchar(keyword,ckeyword)
       do j = 1, ndp / naxis
          write(field,'(a,i0,a)') 'TPV',j,C_NULL_CHAR
          call cchar(field,cfield)
          idx = j + 3
          fval = pv(i,j)
  !        pval = C_LOC(val)
          status = dpfill(C_LOC(dpk(j)),C_LOC(ckeyword),C_LOC(cfield),&
               idx,type,ival,fval)
       end do
    end do
    status = dis_put_pvkey(pdis,DIS_DP,dpk)

    status = disset(0,pdis)

    deallocate(dpk,dis,lin)

  end subroutine wcs_put_dis


!!$  subroutine wcs_ini_dis(wcs) !result(status)
!!$    type(wcspar), intent(in out) :: wcs
!!$    type(c_ptr) :: pdis = C_NULL_PTR
!!$    type(c_ptr) :: plin = C_NULL_PTR
!!$    type(linprm), pointer :: lin
!!$    type(disprm), pointer :: dis
!!$    integer :: status
!!$
!!$!    status = wcsget(wcs % wcs,WCS_LIN,lin)
!!$    !    if( status /= 0 ) return
!!$!    lin => wcs % wcs % lin
!!$!    lin = C_LOC(wcs % wcs % lin)
!!$
!!$    allocate(lin)
!!$    plin = C_LOC(lin)
!!$    status = wcsget(wcs % ptr,WCS_LIN,plin)
!!$    write(error_unit,*) status,'1'
!!$
!!$    allocate(dis)
!!$    pdis = C_LOC(dis)
!!$    status = linget(plin,LIN_DISSEQ,pdis)
!!$!    status = linget(lin,LIN_DISSEQ,dis)
!!$    write(error_unit,*) status,'2'
!!$    if( status /= 0 ) return
!!$
!!$!    status = disget(1,dist,DIS_NAXIS,naxis)
!!$!    status = dis_put(pdis,DIS_NDP,0,0,0)
!!$!    write(error_unit,*) status,'3'
!!$
!!$!    status = distget(dis,DIS_DP,dp)
!!$!    do i = 0, ndp-1
!!$!    end do
!!$
!!$    status = disini(0,2,pdis)
!!$    write(error_unit,*) status,'4'
!!$
!!$    deallocate(lin,dis)
!!$
!!$  end subroutine wcs_ini_dis

end module wcslib
