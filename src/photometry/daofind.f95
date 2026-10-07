!
!     Find  - Search stars on image
!     Copyright (C) 1999, 2010, 2012, 2016, 2018-2025 Filip Hroch,
!                   Masaryk University
!     Copyright (C) 1991 P.B. Stetson, Dominon Astrophysical Observatory
!
!
!  This file is part of Munipack.
!
!  Munipack is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  Munipack is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with Munipack.  If not, see <http://www.gnu.org/licenses/>.
!
!
! 1999, september ..  to F90, no a disk scratch file+other many changes
! 1999, november  ..  rewriting completed
!
!====================================================================
!
!  The original code of DAOPHOT is absolutelly impressive.
!  I just try to improve it.
!

! TODO:
!  * to double precision


module mdaofind

  logical, parameter, private :: debug = .false.

  private :: centres, moments

contains

  subroutine daofind (d, verbose, plog, hibad, skymod, skysig,  &
       fwhm, lothresh, threshold, shrplo, shrphi, rndlo, rndhi, readns, &
       nstar, lobad, hmin )

  implicit none

  real, dimension(:,:), intent(in) :: d
  logical, intent(in) :: verbose, plog
  real, intent(in) :: hibad, fwhm, shrplo, shrphi, rndlo, rndhi, &
       readns, lothresh, threshold, skymod, skysig
  integer, intent(out) :: nstar
  real, intent(out) :: lobad, hmin

  ! Gain: DAOPHOT scales data on-the-fly by the parameter,
  ! while Munipack assumes that the data has already gain correction done.
  real, parameter :: phpadu = 1.0

!
!=======================================================================
!
! This subroutine is supposed to find small, positive brightness
! perturbations in a two-dimensional image.
!
!                OFFICIAL DAO VERSION:  1991 April 18
!
! First, FIND reads in several rows' worth of image data.  For each
! pixel it computes a least-squares fit of an analytic Gaussian function
! to a roughly circular array of pixels surrounding the pixel in
! question.  The overall bias level (sky brightness in that vicinity)
! is removed by the calculation and, since the function is
! symmetric about the central pixel, a smooth gradient in the sky
! brightness cancels out exactly.  This means that the user does
! not have to specify an absolute brightness threshold for star
! detection, and if the mean background brightness varies over the
! frame, to the extent that the variations are smooth and large-scale,
! to first order they will have no effect on the detection limit.
!    The derived peak heights of the Gaussian functions are stored in a
! scratch disk image file.  Later they will be read back in, and local
! maxima in the peak values will be sought.  After undergoing a few
! tests designed to select against bad pixels and bad columns, these
! local maxima will be considered to be astronomical objects, better
! image centroids will be computed, and the objects will be assigned
! sequential ID numbers and will be written to a disk data file.
!    The user is asked to specify a "lowest good data-value"-- any pixel
! whose value is found to fall below this level or above the HIBAD
! value which is passed as an argument is presumed bad, and
! is ignored during all computations in this routine.  The numerical
! value of this bad pixel ceiling will be written out in the header
! of the output data file, and will be used in other DAOPHOT routines
! as well.
!
! Arguments
!
!     FWHM (INPUT) is the estimated full width at half-maximum of the
!          objects for which the algorithm is to be optimized.  It will
!          be used (a) to determine the size of the roughly circular
!          array which will be used to compute the brightness
!          enhancements and to define local maxima, and (b) to define
!          the coefficient assigned to each pixel in the computation
!          of the brightness enhancements.
!
!    WATCH (INPUT) governs whether information relating to the progress
!          of the star-finding is to be typed on the terminal screen
!          during execution.
!
! SHRPLO, SHRPHI (INPUT) are numerical cutoffs on the image-sharpness
!          statistic, designed to eliminate brightness maxima which
!          appear to be due to bad pixels, rather than to astronomical
!          objects.
!
! RNDLO, RNDHI (INPUT) are numerical cutoffs on the image-roundness
!          statistic, designed to eliminate brightness maxima which
!          appear to be due to bad rows or columns, rather than to
!          astronomical objects.
!
! HIBAD  is the highest valid data-value-- the level above which the
!          CCD chip is presumed to be non-linear.
!
! All of the above arguments are user-definable optional parameters,
! whose numerical values may be changed by a DEFAULT.OPT file, or by
! the OPTION command.  (WATCH may also the set by the MONITOR and
! NOMONITOR commands.)
!
!=======================================================================
!
! Parameters:
!
!!!!!  integer :: maxsky, nopt
!
! MAXBOX is the length of the side of the largest subarray that you plan
!        to need for computing the brightness enhancement in each pixel.
! Warning: maxbox is no more need, replaced by dynamic allocation of arrays
!
! MAX/MAXBOX is the length in the x-direction of the largest picture
!        you can try to reduce.
!
!-----------------------------------------------------------------------
!
! DIMENSIONS
!
! Arrays
!

!  real, dimension(size(d,1),size(d,2)) :: h
  real, dimension(:,:), allocatable :: h, g
!  real, dimension(:), allocatable :: sh
! data(2), opt(nopt)
!  real, dimension(:,:), allocatable :: g
  logical, dimension(:,:), allocatable :: skip
!
! Variables
!
!!!!  character(len=*) :: coofil
!!!!!  character(len=5) :: line
!!!!  real pixels, radius, fwhm, sigsq, rsq, relerr, skylvl, temp
!!!!  real hmin, lobad, hibad, watch, p, datum, height, denom, sgop
  real :: pixels, radius, sigsq, rsq, relerr, skylvl, temp, ssig
  real :: p, datum, height, denom, sgop
!!!  real :: lobad, hmin, p, datum, height, denom, sgop
  real :: sharp, round!, shrplo, shrphi, rndlo, rndhi
  real :: sumg, sumgsq, sumgd, sumd, sg, sgsq, sgd, sd, wt, hx, hy
  real :: dgdx, sdgdx, sdgdxs, sddgdx, sgdgdx
  real :: xcen, ycen, dx, dy, shape, ecc, incl
!!!!  real xcen, ycen, dx, dy, phpadu, readns, skymod, skymn, skymed
  integer :: nhalf, nbox, middle, lastcl, lastro, ncol, nrow, jsq
!!!  integer :: istat, nstar
  integer :: i, j, n, ix, iy, jx, jy, kx, ky
  logical :: inside, inrange

!!!!!!  hibad = opt(4)
!!!!!!  fwhm = opt(5)
!!!!!!  shrplo = opt(7)
!!!!!!  shrphi = opt(8)
!!!!!!  rndlo = opt(9)
!!!!!!  rndhi = opt(10)
!!!!!!  watch = opt(11)
!!!!!!

!write(*,*) hibad,fwhm,shrplo,shrphi,rndlo,rndhi

!-----------------------------------------------------------------------
!
! section 1
!
! Setup the necessary variables and arrays, particularly the constants
! to be used in the convolutions.
!
! The brightness enhancement will be computed on the basis only of those
! pixels within 1.5 sigma = 0.637*FWHM of the central pixel.  However,
! in the limit of infinitely small FWHM the brightness enhancement will
! be based on no fewer than the following subarray of pixels:
!
!                                .
!                                .
!                                .
!
!                          -  -  +  -  -
!                          -  +  +  +  -
!                 .  .  .  +  +  X  +  +  .  .  .
!                          -  +  +  +  -
!                          -  -  +  -  -
!                                .
!                                .
!                                .
!
! This represents a 5 x 5 subarray taken out of the original picture.
! The X represents the pixel for which the brightness enhancement is
! currently being computed and the +'s represent other pixels included
! in the calculation; the -'s and all pixels lying outside this 5 x 5
! subarray will not be used in computing the brightness enhancement in
! the central pixel.  In the limit of infinitely large FWHM, only those
! pixels lying within a MAXBOX x MAXBOX square subarray centered on the
! pixel in question will be used in computing its brightness
! enhancement.
!
! Compute the size of the subarray needed.  The radius of the circular
! area desired is MAX (2.0, 0.637*FWHM), so the distance from the
! central pixel to the center of an edge pixel is the integer smaller
! than this.
!
  radius = max(2.001, 0.637*fwhm)
  nhalf = int(radius)
  nbox = 2*nhalf + 1                ! length of the side of the subarray
  middle = nhalf + 1

  allocate(g(nbox,nbox),skip(nbox,nbox))

  ncol = size(d,1)
  nrow = size(d,2)

!
! Just for future reference--
!
! MIDDLE is the index of the central pixel of the box in both x and y,
!        where the corner of the box is considered to be at (1,1).
!
!  NHALF is the number of pixels between the central pixel (exclusive)
!        and the edge of the box (inclusive).  For example, if NBOX = 7,
!        MIDDLE = 4 and NHALF = 3.  Note that all the way around the
!        picture being reduced there will be a border NHALF pixels wide
!        where define brightness enhancements can't be defined, because
!        the box would extend beyond the boundaries of the frame.  We
!        will thus be able to compute brightness enhancements only for
!        MIDDLE <= x <= LASTCL,   MIDDLE <= y <= LASTRO, where...
!
  lastro = nrow - nhalf
  lastcl = ncol - nhalf
!
!-----------------------------------------------------------------------
!
! Compute the values of a bivariate circular Gaussian function with
! unit height and the specified value of the FWHM.
!
  sigsq = (fwhm/2.35482)**2
  radius = radius**2
!
! RADIUS is now the square of the radius of the circle to be used.
!
!-----------------------------------------------------------------------
!
! EXPLANATION:
!
! The approach taken by this star-finding algorithm is defined by this
! question:  "Assuming for the moment that there is a star with a
! Gaussian light distribution centered in the central pixel of this
! subarray, then how bright is it?"  Having answered that question for
! every pixel MIDDLE <= x <= LASTCL, MIDDLE <= y <= LASTRO, we will
! then go through the picture looking for places where the numerical
! answer to the question achieves local maxima.  For the region around
! each pixel, then, we want to solve this equation via least squares:
!
!                     D(i,j) = h * G(i,j) + s
!
! where D is the observed brightness in some pixel of the subarray, G
! is the value of the Gaussian function of unit central height in the
! in that pixel
!
! G(i,j) = exp{[(i-MIDDLE)**2 + (j-MIDDLE)**2]/(2 * sigma**2)}, for
!
!                      (i-MIDDLE)**2 + (j-MIDDLE)**2 < (1.5 * sigma)**2
!
! (the center of the subarray has relative coordinates i = j = MIDDLE).
!
!      The parameters  h  (= central brightness of the hypothetical
! star centered in the central pixel of the subarray), and s (= the
! local sky background) are unknowns.  The least-squares solution
! for this system of equations is given by
!
!         [G*D] - [G] [D]/n
!    h =  ----------------- ,         s = {[D] - h [G]}/n
!         [G**2] - [G]**2/n
!
! where the square brackets denote summation (Gauss's notation).
!
! For use in solving for the many values of  h, we will save the
! array G(i,j) (= G(I,J)) and the constants [G] (= SUMG, meaning
! "sum of the Gaussian"), [G**2] (= SUMGSQ), n (= PIXELS); also the
! denominator of the fraction for  h (= DENOM), and [G]/n (= SGOP).
! [G*D] and [D] will have to be computed each time.
!
! It is possible to show that each of these least-squares problems can
! be reduced to a linear function of the image data D(i,j), and that the
! entire ensemble of least-squares problems is arithmetically identical
! with a convolution of the original image data with a truncated,
! lowered Gaussian function.  Hence, I will occasionally refer to the
! generation of the array of values h(i,j) as a "convolution."
!
!-----------------------------------------------------------------------
!
! Loop over the pixels in the subarray, computing the value of the
! Gaussian function G(i,j) at each point.  Also, accumulate the sum of
! the values of the Gaussian and the sum of the squares of the values
! of the Gaussian.  These will be held for later use in the convolution.
!
  sumg = 0.0
  sumgsq = 0.0
  pixels = 0.0
  do j = 1, nbox
     jsq = (j - middle)**2
     do i = 1, nbox
        rsq = (i - middle)**2 + jsq
        g(i,j) = exp(-0.5*rsq/sigsq)
        if (rsq <= radius) then
           skip(i,j) = .false.
           sumg = sumg + g(i,j)
           sumgsq = sumgsq + g(i,j)**2
           pixels = pixels + 1.0
        else
           skip(i,j)= .true.
        end if
     end do
  end do
  denom = sumgsq - (sumg**2)/pixels
  sgop = sumg/pixels
!
! At this point the two-dimensional array G(I,J) contains the values of
! a unit Gaussian function, with the input value of FWHM, at each point
! in the SQUARE subarray.
!
! SUMG   contains the sum of the values of the Gaussian function over
!        the CIRCULAR area which will be used in the convolution.
!
! SUMGSQ contains the sum of the squares of the values of the Gaussian
!        function over the CIRCULAR area which will be used in the
!        convolution.
!
! PIXELS contains the number of pixels in the CIRCULAR area which will
!        be used in the convolution.
!
! DENOM  contains the denominator of the fraction defining  h.
!
! SGOP   contains [G]/n
!
! Using our knowledge of least squares, we can compute the standard
! error of the coefficient  h  in terms of the standard error of the
! brightness in a single pixel:
!
!      sigma**2(h) = sigma**2(1 pixel) / ([G**2] - [G]**2/n)
!
  relerr = 1.0/denom
  relerr = sqrt(relerr)
  !  call daosky(ncol, nrow, d, min(maxsky, (ncol*nrow)/3), opt(4), skymn, skymed, skymod, ix)
  !  call daosky(d,min(maxsky,(ncol*nrow)/3),verbose,hibad,skymod,skyerr,skysig)

!  if( ncol*nrow > maxsky ) then
!     n = 2
!  else
!     n = 1
!  end if
!  call daosky(d,n,verbose,hibad,skymod,skyerr,skysig)
!  write(*,*) skymod,skyerr,skysig,relerr

  if( verbose ) write (*,"(23X, 'Relative error = ', F5.2/)") relerr

!  if( abs(skymod) < epsilon(skymod) ) then
!     write(*,*) 'Zero sky standard deviations: defective frame(?).'
!     nstar = 0
!     return
!  end if

!
! Now ask the user for a star-detection threshold and a bad pixel
! ceiling.
!
!      call getdat ('Number of frames averaged, summed:', DATA, 2)
!!!!!  if (data(1) < 0.5 .or. data(2) < 0.5) return
!!!!!  readns = opt(1)**2*data(2)/data(1)
!!!!!  phpadu = opt(2)*data(1)

  if( .not. threshold > 0 ) &
       stop 'Error in FIND: threshold > 0 assumption failed.'

  ssig = (phpadu*readns)**2 + skysig**2
  if( ssig > 0 ) then
     ssig = sqrt(ssig)
  else if( skymod > 0 ) then
     ! suitable for a test data without included Poisson noise
     ssig = sqrt(readns**2 + max(0.0,skymod)/phpadu) ! assuming: var(d)==mean(d)
  else
     ssig = 1 !sqrt(epsilon(hmin))
  end if
!  lobad = 0.1*nint(10.*(skymod - lothresh*hmin))
! lobad = max(0.0,skymod - lothresh*hmin)  ! bad idea !!!
  lobad = skymod - lothresh*ssig
!  write(*,*) 0.1*nint(10.*(skymod - lothresh*hmin)),max(0.0,skymod - lothresh*hmin)
!  hmin = 0.01*nint(100.*threshold*relerr*hmin)
!  write(*,*) skymod,lobad,hmin
  hmin = threshold*relerr*ssig
!  write(*,*) lobad,hmin,relerr,ssig

  ! simplified code:
!  hmin = sqrt(readns**2 + skysig**2) ! replacing skymod->skysig**2 is realiable
!  lobad = max(0.0,skymod - lothresh*hmin)
!  hmin = threshold*relerr*hmin
!  write(*,*) lobad,hmin,relerr
!!!!!  lobad = 0.1*nint(10.*(skymod - opt(3)*hmin))
!!!!!  hmin = 0.01*nint(100.*opt(6)*relerr*hmin)
!!!!!  readns = sqrt(readns)
!
! Later on, the threshold HMIN will be the minimum value of the local
! brightness enhancement that will be considered when searching for
! local maxima, and any pixel whose brightness value is less than LOBAD
! or greater than HIBAD will be ignored in the computations.
!
! Open the input and scratch disk files.
!
! Open output data file for newly-discovered stars.
!
!!!!!!  call outfil (3, coofil, istat)
!!!!!!  if (istat /= 0) then
!!!!!!     call stupid ('Error opening output file '//trim(coofil))
!!!!!!     return
!!!!!!  end if
!!!!!!  if (watch > 0.5) then
!!!!!!     call tblank
!!!!!!!     call ovrwrt('  Row', 1)
!!!!!!  end if
!!!!!!!
!-----------------------------------------------------------------------
!
! SECTION 2
!
! Read the raw image data in, holding only a few rows' worth of data in
! memory at any one time.  Convolve the data with the appropriate
! Gaussian function, and write the resulting numbers into the scratch
! disk picture.
!
! .... (censored :-)
!
! The cylinder buffer D now contains the actual image data for the
! first NHALF rows of the picture.  We will soon create the file
! containing the derived values of  h  (see above) one row at a time.
!
! Now we will step through the picture row by row. JY remembers which
! row in the big picture we are working on.  For each row JY, the
! convolved data will be accumulated in the vector H(i,2), and then
! written into the JY-th row of the scratch picture.
!
!      jy = 0
! 2020 jy = jy + 1                              ! Increment image-row pointer
!      if (jy > nrow) go to 2100         ! Have we reached the bottom?
!
! ... (censored :-)
! Note that at any given time we have only NBOX rows of the original
! image in memory, contained in the cylinder buffer D(i,j),
! j = 1, ..., NBOX, but not necessarily in that order.  For instance,
! if NBOX = 5, when JY = 1,
!
!      row:     1    2    3    4    5     of G is to be fitted to
!
!      row:     *    *    1    2    3     of the original picture which
!                                         is contained in
!      row:     1    2    3    4    5     of D.
!
! When row 1 of the picture is done, JY is set to 2, and row 4
! of the original picture is read into row 1 of the cylinder buffer,
! D, overwriting the null values which we put there before.
! Hence, when JY = 2,
!
!      row:     1    2    3    4    5     of G is to be fitted to
!
!      row:     *    1    2    3    4     of the original picture which
!                                         is contained in
!      row:     2    3    4    5    1     of D.
!
! As a final example, consider the situation for JY = 7:
!
!      row:     1    2    3    4    5     of G is to be fitted to
!
!      row:     5    6    7    8    9     of the original picture which
!                                         is contained in
!      row:     2    3    4    5    1     of D.
!
! In other words:
!
!      row:     1    2    3    4    5     variable J
!
!      row:     5    6    7    8    9     variable JY
!
!      row:     2    3    4    5    1     vector JCYLN(J)
!
! The cylinder buffer, D, just rolls down through the picture like a
! caterpillar tread, dropping off rows of data when they are no longer
! necessary and picking up new ones in their place.  The data are
! handled in this way (a) to minimize the amount of memory required,
! by storing only those rows that are immediately wanted, consistent
! with (b) minimizing the number of data transfers.  Now, for the
! CURRENT value of JY, which row of the cylinder buffer is to be fitted
! to each row of G?  The answers will be contained in the vector
! JCYLN.
!
! JCYLN(MIDDLE) is the row in the cylinder buffer where we will find
! the data for row JY of the big picture, which is to be fitted to row
! MIDDLE of G.  Similarly, JCYLN returns the position in the cylinder
! buffer of the row to be fitted to the J-th row of G
! (J = 1, ..., NBOX).
!
! Now that this is all straight, read in the data for row JY+NHALF
! (overwriting the data for row JY-NHALF-1, which is no longer needed).
!
!      do j = 1, nbox
!         iy = jy + (j - middle)
!
! iy is that row of the big picture which is to be matched up against
! row J of the Gaussian function, during the convolution of this
! row JY of the big picture.
!
! Which row of the cylinder buffer contains row IY of the big picture?
!
!         i = iy + nhalf
!
! i now represents the position that row IY of the big picture would
! have had in the cylinder buffer if the cylinder buffer were
! arbitrarily long, i.e. row 1 of the image in row 3 of D, row 2
! in row 4, row 3 in row 5, row 4 in row 6, ... in the examples
! above.  Now we wrap this around.
!
!         jcyln(j) = mod(i-1,nbox) + 1
!      end do
!
!      ly = jy+nhalf
!      if (ly >= nrow) then
!         call rdaray ('DATA', lx, ly, ncol, nrows, maxcol, d(1,jcyln(nbox)), istat)
!      else
!         k = jcyln(nbox)
!         do ix=1,ncol
!            d(ix,k) = -1.1e38
!         end do
!      end if
!
!      if (istat /= 0) then
!         call stupid ('Error reading image data from disk file.')
!         return
!      end if
!
! Compute the local brightness enhancement for each pixel in the row,
! The enhancement is computed from a circular region contained
! within an NBOX x NBOX array centered on the current pixel, using the
! array, G(I,J), and the constants SUMG, SUMGSQ, and PIXELS computed
! above.  (These constants will need to be modified if the circular
! region used in the calculation contains any bad pixels; we will use
! the variables SG, SGSQ, and P for temporary storage of these
! constants, and SGD and SD for the accumulation of [G*D] and [D] which
! are also needed.)
!

!  if( plog ) write(*,'(a)') "=FIND> Convolution ..."
  allocate(h(ncol,nrow))

  do jy = 1, nrow

     do jx = 1, ncol

        sgd = 0.0
        sd = 0.0
        sgsq = sumgsq
        sg = sumg
        p = pixels

        do  ix = jx - nhalf, jx + nhalf
           i = middle + (ix - jx)
           do  iy = jy - nhalf, jy + nhalf
              j = middle + (iy - jy)
              if ( .not. skip(i,j)) then
                 inside = (1 <= ix  .and. ix <= ncol ) .and. &
                      ( 1 <= iy  .and. iy <= nrow )
                 inrange = .false.
                 if( inside ) &
                      inrange = (lobad <= d(ix,iy) .and. d(ix,iy) <= hibad )
!                 if( ( 1 <= ix  .and. ix <= ncol ) .and. &
!                     ( 1 <= iy  .and. iy <= nrow ) ) then
                 !                    if( (lobad <= d(ix,iy) .and. d(ix,iy) <= hibad ) )then
!                 inrange = .true.
                 if( inside .and. inrange ) then
                       datum = d(ix,iy)
                       sgd = sgd + g(i,j)*datum
                       sd = sd + datum
!                    end if
                 else
                    sgsq = sgsq - g(i,j)**2
                    sg = sg - g(i,j)
                    p = p - 1.0
                 endif
              endif
           enddo
        enddo
!
! compute the central height of the best fitting Gaussian function,
! temporarily storing it in the variable, then putting it into array
! element H(JX, 2).
!
        if (p > 1.5) then
           if (p < pixels) then
              sgsq = sgsq - (sg**2)/p
              if (sgsq /= 0.0) then
                 sgd = (sgd - sg*sd/p)/sgsq
              else
                 sgd = 0.0
              end if
           else
              sgd = (sgd - sgop*sd)/denom
           end if
        else
           sgd = 0.0
        end if
        h(jx,jy) = sgd
     enddo !jx
!
! Write this newly-computed row of brightness enhancements to the
! scratch output picture.
!
!!!!!!      if (watch > 0.5) then
!!!!!!         write (line,"(I5)") jy
!!!!!!!         call ovrwrt (line(1:5), 2)
!!!!!!      end if
  enddo ! jy

!  write(*,*) h(50,136),lobad,hibad

  ! save convoluted image
!  if( debug ) then
!     i = 0
!     Call ftinit(26,'!h.fits',1,i)
!     Call ftphps(26,-32,2,(/ncol,nrow/),i)
!     Call ftp2de(26,1,ncol,ncol,nrow,h,i)
!     call ftclos(26,i)
!  end if


!!!!!!  call ovrwrt (' ', 4)
!!!!!!!
! Later on, when we try to decide whether a local maximum represents
! a stellar profile or a delta function ( = bright bad pixel), we will
! compare the brightness of the central pixel to the average of the
! surrounding pixels.  To be ready for that, we here modify SKIP to
! skip over the central pixel, and set PIXELS equal to the number of
! pixels in the circular area not counting the central pixel.
!
  skip(middle,middle) = .true.
  pixels = pixels - 1.0

!  do i = 1,nbox
!     write(*,*) skip(i,:)
!  end do
!  stop 0
!
!-----------------------------------------------------------------------
!
! SECTION 3
!
! Read in both the convolved data from the scratch disk file and the raw
! data from the original picture.  Search for local maxima in the
! convolved brightness data.  When these are found, compute image-shape
! statistics from the raw data to eliminate non-stellar brightness
! enhancements (as well as possible) and estimate the position of the
! centroid of the brightness enhancement.
!
! 3000 continue
!
! Now the star search may begin.  The original image data will be read
! into the cylinder buffer D again, just as before.  At the same time,
! the brightness enhancements will be read from the scratch disk file
! into another cylinder buffer, H.  The brightness enhancements will
! then be searched for local maxima.  When these are found, functions
! of the original image data will be used to derive shape parameters
! designed to identify bad pixels and bad columns or rows.
!
!!!!!!  call wrhead (3, 1, ncol, nrow, 6, lobad, hibad, hmin, 0., phpadu, readns, 0.)
!!!!!!  if (watch > 0.5) then
!     call ovrwrt (' ', 4)
  if( verbose ) then
!     write (*,"(6X, '                              MAGS')" )
!     write (*,"(6X, '                              FROM ')")
     write (*,"(6X, '                             PEAK')")
     write (*,"(6X, '       STAR     X      Y     RATIO      SHARP    ROUND')")
  end if
!!!!!!!
! .... (censored :-)
!
! Now step through the picture row by row.  Again JY is the image-row
! counter.
!
  nstar = 0

  do jy = 1, nrow

!
! .... (censored :-)
!
!
! Now step across the row, pixel by pixel.
!
     jx = 1
     do
        height = h(jx,jy)
!
! sieve to locate a local maximum in the brightness enhancement.  To
! be a local maximum, the brightness enhancement in a given pixel must
! be above the threshold, and it must also be greater than the
! brightness enhancement of any pixel within a radius equal to
! 1.5 sigma.
!
        if (height < hmin) go to 3200
        do ix = jx - nhalf, jx + nhalf
           if ((1 <= ix) .and. (ix <= ncol)) then
              i = middle + (ix - jx)
              do iy = jy - nhalf, jy + nhalf
                 j = middle + (iy - jy)
!                 inside = 1 <= ix .and. ix <= ncol .and. &
!                      1 <= iy .and. iy <= nrow
                 if( 1 <= iy .and. iy <= nrow ) then
!                    if ( .not. skip(i,j) .and. ) then
                    if ( .not. skip(i,j) ) then
                       if (height < h(ix,iy)) go to 3200
                    endif
                 end if
              enddo
           endif
        enddo
!
! The brightness enhancement of this pixel is now confirmed to be above
! the threshold, and to be larger than in any other pixel within a
! radius of 1.5 sigma.
!
! Now we derive the shape indices.  First, is the object much more
! sharply peaked than the input FWHM?  Compare the central pixel to
! the mean of the surrounding (non-bad) pixels.  If this difference is
! greater than the originally estimated height of the Gaussian or less
! than two-tenths the height of the Gaussian, reject the star
! (assuming SHRPLO and SHRPHI have the default values of 0.2 and
! 1.0; otherwise, muta mutandis.)
!
!
! ********** IF THE CENTRAL PIXEL IS BAD SKIP THIS TEST. **********
!
!
!D     TYPE *, JX, JY
!D     DO 1666 J=1,NBOX
!D1666 TYPE 6661, (JNINT(D(I,JCYLN(J))),
!D    .     I=MAX0(1,JX-NHALF),MIN0(NCOL,IX+NHALF)),
!D    .     (JNINT(H(I,JCYLN(J))), I=IX-NHALF,IX+NHALF)
!D6661 FORMAT(1X, <NBOX>I6, 1X, <NBOX>I6)
!
! As one final nuance, for this and subsequent calculations I propose
! to subtract off the modal sky level.  Otherwise, for faint stars on
! bright backgrounds in large boxes, it is barely possible that
! truncation error could affect the numerical results of the analysis.
!
!        allocate(sh(nbox*nbox))
!        n = 0
        sharp = 0.0
        datum = d(jx,jy)
        if( lobad <= datum .and. datum <= hibad ) then
           p = 0.0
           do ix = jx - nhalf, jx + nhalf
              if( 1 <= ix .and. ix <= ncol ) then
                 i = middle + (ix - jx)
                 temp = 0.0
                 do iy = jy - nhalf, jy + nhalf
                    j = middle + (iy - jy)
                    if( 1 <= iy .and. iy <= nrow ) then
                       !                    if( 1 <= i .and. i <= ncol .and. 1 <= j .and. j <= nrow)then
                       if ( .not. skip(i,j) )then
                          datum = d(ix,iy)
                          if ((datum >= lobad) .and. (datum <= hibad)) then
                             temp = temp + (datum - skymod)
                             p = p + 1.0
!                             n = n + 1
!                             sh(n) = datum - skymod
                          end if
                       endif
                    end if
                 enddo
                 sharp = sharp + temp
              endif
           enddo

           ! the original estimator
           if( height > 0 .and. p > 0 ) then
              sharp = (d(jx,jy) - skymod - sharp/p)/height
           else
              sharp = 0
           end if

           ! new experimental method for estimate of local background
           ! with outliers, looks more reliable
!           call rmean(sh(1:n),temp)
!           write(*,*) jx,jy,skymod,(temp/(sharp/p))
!           sharp = (d(jx,jy) - skymod - temp) / height

        endif
!        deallocate(sh)
        if ( sharp < shrplo .or. sharp > shrphi ) go to 3200
!           write(*,*) sharp,skysig,datum

!
! Now check to see whether the object is strongly elongated either
! along the row or along the column.  Compute the height of a Gaussian
! function of x and a Gaussian function of y by least-squares fits to
! the marginal distributions of the image data.  That is, fit the
! sum over y of the actual brightness values to the sum over y of the
! values of the array G, as functions of x.  If a bad pixel is found
! omit both the picture datum and the value of G for that pixel from
! their respective sums.  If the computed height of either the
! x-marginal or the y-marginal is non-positive, or if the central
! heights of the two marginals differ by more than their average
! (assuming that RNDLO and RNDHI have their default values
! of -1.0 and 1.0; otherwise, etc.), reject the star.
!
! We will now compute the height of the one-dimensional Gaussian
! distribution which best fits the x-marginal distribution of the
! brightness.  The equation which will be used will be the same as
! in the comments above ( h = ...) except that the symbol D in the
! equation now represents stands for the brightness data in the NBOX by
! NBOX square array summed over the y spatial direction, and the
! symbol G now stands for a one-dimensional Gaussian function (= the
! two-dimensional function G(i,j) also summed over the y spatial
! direction.  This sum is actually carried out numerically, rather
! than being done analytically, in order to permit the omission of
! "bad" pixels.)  At the same time, we will set up the necessary sums
! to permit the computation of a first-order correction to the centroid
! of the Gaussian profile in x:
!
!                -[G'*(D-G)]          [G*G']-[D*G']
! Delta x = -------------------- = -------------------,
!            [G'**2] - [G']**2/n   [G'**2] - [G']**2/n
!
! where G is the one-dimensional Gaussian profile, G' = (dG/dx), and
! D = the summed actual image data.  (There would normally be a
! [G']*[(D-G)]/n term in the numerator, but because G is already the
! "best fitting" Gaussian, [(D-G)] = 0.)  We will use
!
! SD      for the marginal sum of the actual image data
!                    (mnemonic:  "temporary sum of the data")
! SG      for the marginal sum of the 2-D Gaussian function
!                               ("temporary sum of the Gaussian")
! SUMGD   for [G*D]             ("sum of the Gaussian times the data")
! SUMG    for [G]               ("sum of the Gaussian")
! SUMD    for [D]               ("sum of the data")
! SUMGSQ  for [G**2]            ("sum of the Gaussian squared")
! SDGDX   for [G']              ("sum of d(Gaussian)/dx")
! SDGDXS  for [G'**2]           ("sum of {d(Gaussian)/dx}**2")
! SDDGDX  for [D*G']            ("sum of data times d(Gaussian)/dx")
! SGDGDX  for [G*G']            ("sum of Gaussian times d(Gaussian)/dx")
!
! In addition, for these calculations, pixels will arbitrarily be
! assigned weights ranging from unity at the corners of the box to
! MIDDLE**2 at the center (e.g. if NBOX = 5 or 7, the weights will be
!
!                                 1   2   3   4   3   2   1
!      1   2   3   2   1          2   4   6   8   6   4   2
!      2   4   6   4   2          3   6   9  12   9   6   3
!      3   6   9   6   3          4   8  12  16  12   8   4
!      2   4   6   4   2          3   6   9  12   9   6   3
!      1   2   3   2   1          2   4   6   8   6   4   2
!                                 1   2   3   4   3   2   1
!
! respectively).  This is done to desensitize the derived parameters to
! possible neighboring, brighter stars.
!
! The temporary variable P will be used to accumulate the sum of the
! weights, and N will count the number of points in the marginal
! distribution that actually get used.
!
! SKIP ALL THIS IF THE STAR IS TOO NEAR THE EDGE OF THE FRAME!
!
        if ((jx < middle) .or. (jx > lastcl) .or.    &
             (jy < middle) .or. (jy > lastro)) then
           xcen = jx
           ycen = jy
           round = 0
           ecc = 0
           incl = 0

           go to 3190
!           goto 3200
        end if

        ix = jx - middle
        iy = jy - middle

        sumgd = 0.0
        sumgsq = 0.0
        sumg = 0.0
        sumd = 0.0
        sdgdx = 0.0
        sdgdxs = 0.0
        sddgdx = 0.0
        sgdgdx = 0.0
        p = 0.0
        n = 0
        do i = 1,nbox
           sg = 0.0
           sd = 0.0
           kx = ix + i
           do j = 1,nbox
              wt = middle - abs(j - middle)
              ky = iy + j
              if( 1 <= kx .and. kx <= ncol .and. 1 <= ky .and. ky <= nrow ) then
                 datum = d(kx,ky)
                 if ((datum >= lobad) .and. (datum <= hibad)) then
                    sd = sd + (datum - skymod)*wt
                    sg = sg + g(i,j)*wt
                 end if
              end if
           enddo
           if (sg > 0.0) then
              wt = middle - abs(i - middle)
              sumgd = sumgd + wt*sg*sd
              sumgsq = sumgsq + wt*sg**2
              sumg = sumg + wt*sg
              sumd = sumd + wt*sd
              p = p + wt
              n = n + 1
              dgdx = sg*(middle - i)
              sdgdxs = sdgdxs + wt*dgdx**2
              sdgdx = sdgdx + wt*dgdx
              sddgdx = sddgdx + wt*sd*dgdx
              sgdgdx = sgdgdx + wt*sg*dgdx
           end if
        enddo
!
! we need at least three points to estimate the height and position
! of the star, and the local sky brightness.
!
        if (n <= 2) go to 3200
        hx = (sumgd - sumg*sumd/p)/(sumgsq - (sumg**2)/p)
!
! dx is the height of the best-fitting marginal Gaussian.  If this is
! non-positive, this is not an acceptable star.
!
        if (hx <= 0.0) go to 3200
!
! compute the first-order correction to the x-centroid of the star.
! Note that a factor of HX/SIGSQ is missing from SDGDX, SDDGDX, and
! SGDGDX, and a factor of (HX/SIGSQ)**2 is missing from SDGDXS.
!
        skylvl = (sumd - hx*sumg)/p
        dx = (sgdgdx - (sddgdx - sdgdx*(hx*sumg + skylvl*p)))/(hx*sdgdxs/sigsq)
        xcen = jx + dx/(1.0 + abs(dx))
!
! if the best estimate of the star's center falls outside the image,
! reject it.
!
        if( xcen < 0.5 .or. xcen > ncol - 0.5 ) go to 3200
!
! Compute the height of the y-marginal Gaussian distribution.
!
        sumgd = 0.0
        sumgsq = 0.0
        sumg = 0.0
        sumd = 0.0
        sdgdx = 0.0
        sdgdxs = 0.0
        sddgdx = 0.0
        sgdgdx = 0.0
        p = 0.0
        n = 0
        do j = 1,nbox
           ky = iy + j
           sg = 0.0
           sd = 0.0
           do i = 1, nbox
              wt = middle - abs(i - middle)
              kx = ix + i
              if( 1 <= kx .and. kx <= ncol .and. 1 <= ky .and. ky <= nrow ) then
                 datum = d(kx,ky)
                 if ( lobad <= datum .and. datum <= hibad ) then
                    sd = sd + (datum - skymod)*wt
                    sg = sg + g(i,j)*wt
                 end if
              end if
           enddo

           if (sg > 0.0) then
              wt = middle - abs(j - middle)
              sumgd = sumgd + wt*sg*sd
              sumgsq = sumgsq + wt*sg**2
              sumg = sumg + wt*sg
              sumd = sumd + wt*sd
              p = p + wt
              dgdx = sg*(middle - j)
              sdgdx = sdgdx + wt*dgdx
              sdgdxs = sdgdxs + wt*dgdx**2
              sddgdx = sddgdx + wt*sd*dgdx
              sgdgdx = sgdgdx + wt*sg*dgdx
              n = n + 1
           end if
        enddo

        if (n <= 2) go to 3200
        hy = (sumgd - sumg*sumd/p)/(sumgsq - (sumg**2)/p)
        if (hy <= 0.0) go to 3200
        skylvl = (sumd - hy*sumg)/p
        dy = (sgdgdx - (sddgdx - sdgdx*(hy*sumg + skylvl*p)))/(hy*sdgdxs/sigsq)
        ycen = jy + dy/(1.0 + abs(dy))
        if ((ycen < 0.5) .or. (ycen > nrow-0.5)) go to 3200

        round = 2.0*(hx - hy)/(hx + hy)
        if( round < rndlo .or. round > rndhi ) go to 3200

        ! estimate of stars centroids
        call centres(d,jx,jy,nbox,middle,skymod,lobad,hibad,xcen,ycen)

        ! estimate of the shape parameter
        call moments(d,jx,jy,nbox,middle,skymod,lobad,hibad,xcen,ycen,shape,ecc,incl)

!
! the fully verified and located star may now be dignified with its own
! ID number.
!
3190    nstar = nstar + 1

!        height = -2.5*log10(height/hmin)
        height = height / max(hmin,1.0)

!!!!!!        if (watch > 0.5) then
        if( verbose ) then
           write (*,"(12X, I5, 2F7.1, 2x,1pg9.3, 0p2F9.2)") &
                nstar, xcen, ycen, height, sharp, round
        endif
!!!!!!
!        write (3,"(I6, 14F9.3)") nstar, xcen, ycen, height, sharp, round
!!!!!!

        ! the data are saved to the scratch file and re-read later
        write(3) xcen, ycen, height, sharp, round, ecc, incl

        if( plog ) write(*,'(a,i0,1x,2(f0.3,1x))') '=FIND> ',nstar,xcen,ycen
!        if( verbose ) write(*,*) nstar,xcen,ycen

3200    continue
!
! If the sieve above (between statements 3040 and 3050) has detected a
! local maximum in the brightness enhancement, whether this enhancement
! was subsequently confirmed to be a star or not, then there is no need
! to check the other pixels in this row between JX+1 and JX+NHALF,
! inclusive, since we know there can't be a local maximum there.
!
!        jx = jx + nhalf
       jx = jx + 1
!
! Have we passed the last pixel in the row?  If not, work on this
! pixel.  If so, go to next row.
!
        if( jx > ncol ) exit

        enddo !jx

      enddo !jy
!
!-----------------------------------------------------------------------
!
! SECTION 4
!
! Find out whether the user is happy.  If so, delete the scratch picture
! and close up shop.  If not, return to the beginning of Section 3.
!
!!!!!!      if (watch <= 0.5) write (*,"(//1X, I5, ' stars.')") nstar
      if ( verbose ) write (*,"(//1X, I5, ' stars.')") nstar
!      if( plog ) write(*,'(a,i0,a)') "=FIND> Found ",nstar," stars."
!!!!!!
!!!!!!      call clfile (3)
!!!!!!      call tblank                                    ! type a blank line
!!!!!!      call tblank                                    ! type a blank line
!!!!!!
!-----------------------------------------------------------------------
!
! Normal return.
!

      deallocate(h,g,skip)

    end subroutine daofind


    subroutine centres(d,jx,jy,nbox,middle,skymod,lobad,hibad,xcen,ycen)

      real, dimension(:,:), intent(in) :: d
      integer, intent(in) :: jx,jy,nbox,middle
      real, intent(in) :: skymod,lobad,hibad
      real, intent(inout) :: xcen,ycen

      integer :: i,j,ix,iy,ky,kx,ncol,nrow
      real :: wt,w,s

      ix = jx - middle
      iy = jy - middle
      ncol = size(d,1)
      nrow = size(d,2)

      w = 0
      s = 0
      do j = 1,nbox
         ky = iy + j
         if( 1 <= ky .and. ky <= nrow ) then
            do i = 1, nbox
               kx = ix + i
               if( 1 <= kx .and. kx <= ncol ) then
                  wt = d(kx,ky) - skymod
                  if ( lobad < d(kx,ky) .and. d(kx,ky) < hibad .and. wt > 0 ) then
!                     wt = wt*(middle - abs(i - middle))
                     w = w + wt
                     s = s + wt*kx
                  end if
               end if
            end do
         end if
      end do
      if( w > 0 ) then
         xcen = s / w
      else
         xcen = jx
      end if

      w = 0
      s = 0
      do j = 1,nbox
         ky = iy + j
         if( 1 <= ky .and. ky <= nrow ) then
            do i = 1, nbox
               kx = ix + i
               if( 1 <= kx .and. kx <= ncol ) then
                  wt = d(kx,ky) - skymod
                  if ( lobad < d(kx,ky) .and. d(kx,ky) < hibad .and. wt > 0 ) then
!                     wt = wt*(middle - abs(j - middle))
                     w = w + wt
                     s = s + wt*ky
                  end if
               end if
            enddo
         end if
      end do
      if( w > 0 ) then
         ycen = s / w
      else
         ycen = jy
      end if

      if( abs(xcen - jx) > middle .or. abs(ycen - jy) > middle ) then
         write(*,*) 'Warning in FIND: ',jx,jy,xcen,ycen
      end if


!!$      goto 667
!!$      p = 0
!!$        xcen = 0
!!$        do j = 1,nbox
!!$           ky = iy + j
!!$           do i = 1, nbox
!!$              wt = middle - abs(i - middle)
!!$              wt = 1
!!$              kx = ix + i
!!$              datum = d(kx,ky)
!!$              if ( lobad <= datum .and. datum <= hibad ) then
!!$                 xcen = xcen + (datum - skymod)*wt*kx
!!$                 p = p + (datum - skymod)*wt
!!$              end if
!!$           enddo
!!$        end do
!!$        xcen = xcen / p
!!$
!!$        p = 0
!!$        ycen = 0
!!$        do j = 1,nbox
!!$           ky = iy + j
!!$           do i = 1, nbox
!!$              wt = middle - abs(j - middle)
!!$              wt = 1
!!$              kx = ix + i
!!$              datum = d(kx,ky)
!!$              if ( lobad <= datum .and. datum <= hibad  ) then
!!$                 ycen = ycen + (datum - skymod)*wt*ky
!!$                 p = p + (datum - skymod)*wt
!!$              end if
!!$           enddo
!!$        end do
!!$        ycen = ycen / p
!!$        667 continue


    end subroutine centres

    subroutine moments(d,jx,jy,nbox,middle,skymod,lobad,hibad,xcen,ycen,shape,ecc,incl)

      real, dimension(:,:), intent(in) :: d
      integer, intent(in) :: jx,jy,nbox,middle
      real, intent(in) :: skymod,lobad,hibad,xcen,ycen
      real, intent(out) :: shape,ecc,incl

      real, parameter :: rad = 57.29577951308232

      integer :: i,j,ix,iy,ky,kx,ncol,nrow
      real :: wt,a,b,h,hmean,s,w,hx,hy,hxy

      ix = jx - middle
      iy = jy - middle
      ncol = size(d,1)
      nrow = size(d,2)

      s = 0
      w = 0
      do j = 1,nbox
         ky = iy + j
         do i = 1, nbox
            kx = ix + i
            if( (1 <= kx .and. kx <= ncol) .and. (1 <= ky .and. ky <= nrow))then
               wt = d(kx,ky) - skymod
               if ( lobad < d(kx,ky) .and. d(kx,ky) < hibad .and. wt > 0 ) then
                  w = w + wt
                  s = s + wt*(kx - xcen)**2
               end if
            end if
         enddo
      end do
      if( w > 0 ) then
         if( s > 0 ) then
            hx = sqrt(s / w)
         else
            hx = 0
         end if
      else
         hx = -ncol
      end if

      w = 0
      s = 0
      do j = 1,nbox
         ky = iy + j
         do i = 1, nbox
            kx = ix + i
            if( (1 <= kx .and. kx <= ncol) .and. (1 <= ky .and. ky <= nrow))then
               wt = d(kx,ky) - skymod
               if ( lobad < d(kx,ky) .and. d(kx,ky) < hibad .and. wt > 0 ) then
                  w = w + wt
                  s = s + wt*(ky - ycen)**2
               end if
            end if
         enddo
      end do
      if( w > 0 ) then
         if( s > 0 ) then
            hy = sqrt(s / w)
         else
            hy = 0
         end if
      else
         hy = -nrow
      end if

      w = 0
      s = 0
      do j = 1,nbox
         ky = iy + j
         do i = 1, nbox
            kx = ix + i
            if( (1 <= kx .and. kx <= ncol) .and. (1 <= ky .and. ky <= nrow))then
               wt = d(kx,ky) - skymod
               if ( lobad < d(kx,ky) .and. d(kx,ky) < hibad .and. wt > 0 ) then
                  w = w + wt
                  s = s + wt*(kx - xcen)*(ky - ycen)
               end if
            end if
         enddo
      end do
      if( w > 0 ) then
         hxy = s / w
      else
         hxy = 0 !?
      end if


      hmean = (hx + hy) / 2
      if( hmean > 0 ) then
         if( hx > 0 .or. hy > 0 ) then
            h = sqrt((hx - hy)**2 + 4*hxy**2)
         else
            h = 0
         end if
         shape = h / hmean
         a = hmean + h/2
         b = hmean - h/2
         if( a > b .and. a > 0 ) then
            ecc = sqrt(a**2 - b**2) / a
         else
            ecc = 0
         end if
         if( h > 0 ) then
            if( hmean > abs(hxy) .and. 2*abs(hxy) < h ) then
               incl = rad*asin(2*hxy / h) / 2
            else if( hmean <= abs(hxy) .and. abs(hx - hy) < h ) then
               incl = sign(rad,hxy)*acos((hx - hy) / h) / 2
            else
               incl = 0
            end if
         else
            incl = 0
         end if
      else
         shape = 0
         ecc = 0
         incl = 0
      end if
!     write(*,*) 'shape:',shape,hx,hy,hxy,ecc,incl

    end subroutine moments


end module mdaofind
