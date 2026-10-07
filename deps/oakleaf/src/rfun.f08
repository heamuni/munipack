!
!  Robust minimization functions, the generic interfaces
!
!  Copyright © 2011 - 2025 F.Hroch (hroch@physics.muni.cz)
!
!  This file is part of Oakleaf.
!
!  Oakleaf is free software: you can redistribute it and/or modify
!  it under the terms of the GNU Lesser General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  Oakleaf is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU Lesser General Public License for more details.
!
!  You should have received a copy of the GNU Lesser General Public License
!  along with Oakleaf.  If not, see <http://www.gnu.org/licenses/>.
!

module rfun

!  use rfun_REAL32
!  use rfun_REAL64
  use rfun_hampel_REAL32
  use rfun_hampel_REAL64
  use rfun_tukey_REAL32
  use rfun_tukey_REAL64
  use rfun_huber_REAL32
  use rfun_huber_REAL64
  use rnfun_square_REAL32
  use rnfun_square_REAL64

#ifdef HAVE_REAL128
!  use rfun_REAL128
  use rfun_hampel_REAL128
  use rfun_tukey_REAL128
  use rfun_huber_REAL128
  use rnfun_square_REAL128
#endif

  real,parameter,public :: etukey = 2.1227  ! a temporary

  interface huber
     procedure huber_REAL32
     procedure huber_REAL64
#ifdef HAVE_REAL128
     procedure huber_REAL128
#endif
  end interface huber

  interface dhuber
     procedure dhuber_REAL32
     procedure dhuber_REAL64
#ifdef HAVE_REAL128
     procedure dhuber_REAL128
#endif
  end interface dhuber

  interface ihuber
     procedure ihuber_REAL32
     procedure ihuber_REAL64
#ifdef HAVE_REAL128
     procedure ihuber_REAL128
#endif
  end interface ihuber

  interface pshuber
     procedure pshuber_REAL32
     procedure pshuber_REAL64
#ifdef HAVE_REAL128
     procedure pshuber_REAL128
#endif
  end interface pshuber

  interface dpshuber
     procedure dpshuber_REAL32
     procedure dpshuber_REAL64
#ifdef HAVE_REAL128
     procedure dpshuber_REAL128
#endif
  end interface dpshuber

  interface ipshuber
     procedure ipshuber_REAL32
     procedure ipshuber_REAL64
#ifdef HAVE_REAL128
     procedure ipshuber_REAL128
#endif
  end interface ipshuber

  interface hubersqpf
     procedure hubersqpf_REAL32
     procedure hubersqpf_REAL64
#ifdef HAVE_REAL128
     procedure hubersqpf_REAL128
#endif
  end interface hubersqpf


  interface tukey
     procedure tukey_REAL32
     procedure tukey_REAL64
#ifdef HAVE_REAL128
     procedure tukey_REAL128
#endif
  end interface tukey

  interface dtukey
     procedure dtukey_REAL32
     procedure dtukey_REAL64
#ifdef HAVE_REAL128
     procedure dtukey_REAL128
#endif
  end interface dtukey

  interface itukey
     procedure itukey_REAL32
     procedure itukey_REAL64
#ifdef HAVE_REAL128
     procedure itukey_REAL128
#endif
  end interface itukey

  interface tukeys
     procedure tukeys_REAL32
     procedure tukeys_REAL64
#ifdef HAVE_REAL128
     procedure tukeys_REAL128
#endif
  end interface tukeys

  interface tukeysqpf
     procedure tukeysqpf_REAL32
     procedure tukeysqpf_REAL64
#ifdef HAVE_REAL128
     procedure tukeysqpf_REAL128
#endif
  end interface tukeysqpf


  interface square
     procedure square_REAL32
     procedure square_REAL64
#ifdef HAVE_REAL128
     procedure square_REAL128
#endif
  end interface square

  interface dsquare
     procedure dsquare_REAL32
     procedure dsquare_REAL64
#ifdef HAVE_REAL128
     procedure dsquare_REAL128
#endif
  end interface dsquare

  interface isquare
     procedure isquare_REAL32
     procedure isquare_REAL64
#ifdef HAVE_REAL128
     procedure isquare_REAL128
#endif
  end interface isquare

  interface pshampel
     procedure pshampel_REAL32
     procedure pshampel_REAL64
#ifdef HAVE_REAL128
     procedure pshampel_REAL128
#endif
  end interface pshampel

  interface dpshampel
     procedure dpshampel_REAL32
     procedure dpshampel_REAL64
#ifdef HAVE_REAL128
     procedure dpshampel_REAL128
#endif
  end interface dpshampel

  interface ipshampel
     procedure ipshampel_REAL32
     procedure ipshampel_REAL64
#ifdef HAVE_REAL128
     procedure ipshampel_REAL128
#endif
  end interface ipshampel

  interface shampel
     procedure shampel_REAL32
     procedure shampel_REAL64
#ifdef HAVE_REAL128
     procedure shampel_REAL128
#endif
  end interface shampel

  interface dshampel
     procedure dshampel_REAL32
     procedure dshampel_REAL64
#ifdef HAVE_REAL128
     procedure dshampel_REAL128
#endif
  end interface dshampel

  interface ishampel
     procedure ishampel_REAL32
     procedure ishampel_REAL64
#ifdef HAVE_REAL128
     procedure ishampel_REAL128
#endif
  end interface ishampel

  interface hampel
     procedure hampel_REAL32
     procedure hampel_REAL64
#ifdef HAVE_REAL128
     procedure hampel_REAL128
#endif
  end interface hampel

  interface dhampel
     procedure dhampel_REAL32
     procedure dhampel_REAL64
#ifdef HAVE_REAL128
     procedure dhampel_REAL128
#endif
  end interface dhampel

  interface ihampel
     procedure ihampel_REAL32
     procedure ihampel_REAL64
#ifdef HAVE_REAL128
     procedure ihampel_REAL128
#endif
  end interface ihampel

  interface hampelsqpf
     procedure hampelsqpf_REAL32
     procedure hampelsqpf_REAL64
#ifdef HAVE_REAL128
     procedure hampelsqpf_REAL128
#endif
  end interface hampelsqpf

end module rfun
