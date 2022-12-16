!
! Copyright (C) 2002 FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


!----------------------------------------------------------------------!
! FFT scalar drivers Module.
! Written by Carlo Cavazzoni 
! Last update April 2003
!----------------------------------------------------------------------!


#include "machine.h"


!=----------------------------------------------------------------------=!
   MODULE fft_scalar
!=----------------------------------------------------------------------=!

        USE kinds

        IMPLICIT NONE
        SAVE

        PRIVATE
        PUBLIC :: cft_1z, cft_2xy, cft_b, cfft3d, cfft3ds
        PUBLIC :: good_fft_dimension, allowed, good_fft_order

! ...   Local Parameter

        INTEGER, PARAMETER :: ndims = 3          

        !   ndims   Number of different FFT tables that the module 
        !           could keep into memory without reinitialization

#if defined __AIX

        INTEGER, PARAMETER :: nfftx = 2049
        INTEGER, PARAMETER :: lwork = 20000 + ( 2 * nfftx + 256 ) * 64 + 3 * nfftx
        INTEGER, PARAMETER :: ltabl = 20000 + 3 * nfftx

        !   see the ESSL manual ( DCFT ) for the workspace and table lenght formulas

#else

        INTEGER, PARAMETER :: nfftx = 1025
        INTEGER, PARAMETER :: lwork = 20 * nfftx
        INTEGER, PARAMETER :: ltabl = 4 * nfftx

#endif

#if defined __SGI64 || defined __ALPHA

        INTEGER, PARAMETER :: ipt = 8

#else

        INTEGER, PARAMETER :: ipt = 4

#endif

        !   ipt   Size of integer that store "C" pointers
        !         ipt = 4 for 32bit executables
        !         ipt = 8 for 64bit executables

#if defined __FFTW

        INTEGER ( kind=ipt ) :: fw_plan( 3, ndims ) = 0
        INTEGER ( kind=ipt ) :: bw_plan( 3, ndims ) = 0

        !   Pointers to the "C" structures containing FFT factors ( PLAN )

#elif defined __AIX

        REAL (dbl) :: fw_table( ltabl, 3, ndims )
        REAL (dbl) :: bw_table( ltabl, 3, ndims )
        REAL (dbl) :: work( lwork ) 

        !   Array containing FFT factors + workspace

#elif defined __SGI 

        REAL (dbl) :: work(lwork) 
        REAL (dbl) :: tablez(ltabl,ndims)
        REAL (dbl) :: tablex(ltabl,ndims)
        REAL (dbl) :: tabley(ltabl,ndims)

        !   Array containing FFT factors + workspace

#endif

        REAL (dbl) :: scale



!=----------------------------------------------------------------------=!
   CONTAINS
!=----------------------------------------------------------------------=!

!
!=----------------------------------------------------------------------=!
!
!
!
!         FFT along "z" 
!
!
!
!=----------------------------------------------------------------------=!
!

   SUBROUTINE cft_1z(c, nsl, nz, ldc, sgn, cout)

!     driver routine for m 1d complex fft's 
!     nx=n+1 is allowed (in order to avoid memory conflicts)
!     A separate initialization is stored each combination of input sizes
!     NOTA BENE: the output in fout !

     INTEGER, INTENT(IN) :: sgn
     INTEGER, INTENT(IN) :: nsl, nz, ldc
     COMPLEX (dbl) :: c(:), cout(:) 
     REAL(dbl)  :: tscale
     INTEGER    :: i, j, idir, ip, isign
     INTEGER, SAVE :: zdims( 3, ndims ) = -1
     INTEGER, SAVE :: icurrent = 1


     IF( nsl < 0 ) THEN
       CALL errore(" fft_scalar: cft_1 ", " nsl out of range ", nsl)
     END IF

     isign = -sgn

     !
     !   Here initialize table only if necessary
     !

     ip = -1
     DO i = 1, ndims

       !   first check if there is already a table initialized
       !   for this combination of parameters

       IF( ( nz == zdims(1,i) ) .and. ( nsl == zdims(2,i) ) .and. ( ldc == zdims(3,i) ) ) THEN
         ip = i
         EXIT
       END IF

     END DO

     IF( ip == -1 ) THEN

       !   no table exist for these parameters
       !   initialize a new one

       ! WRITE(6, fmt="('DEBUG cft_1z, reinitializing tables ', I3)" ) icurrent

#if defined __FFTW

       IF( fw_plan( 3, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( fw_plan( 3, icurrent) )
       IF( bw_plan( 3, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( bw_plan( 3, icurrent) )
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 3, icurrent), nz, idir) 
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 3, icurrent), nz, idir) 

#elif defined __SGI

       CALL ZFFT1DI( nz, tablez(1,icurrent) )

#elif defined __AIX

       tscale = 1.0d0 / nz
       CALL DCFT ( 1, c(1), 1, ldc, cout(1), 1, ldc, nz, nsl,  1, &
          tscale, fw_table(1, 3, icurrent), ltabl, work(1), lwork)
       CALL DCFT ( 1, c(1), 1, ldc, cout(1), 1, ldc, nz, nsl, -1, &
          1.0d0, bw_table(1, 3, icurrent), ltabl, work(1), lwork)

#else 

       CALL errore(' cft_1 ',' no scalar fft driver specified ', 1)

#endif

       zdims(1,icurrent) = nz; zdims(2,icurrent) = nsl; zdims(3,icurrent) = ldc;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !

#if defined __FFTW

     IF (isign > 0) THEN
       tscale = 1.0d0 / nz
       CALL FFT_Z_STICK(fw_plan( 3, ip), c(1), ldc, nsl)
       CALL ZDSCAL( ldc * nsl, tscale, c(1), 1)
     ELSE IF (isign < 0) THEN
       CALL FFT_Z_STICK(bw_plan( 3, ip), c(1), ldc, nsl)
     END IF

#elif defined __SGI

     IF (isign /= 0) THEN
       IF( isign < 0 ) idir = +1
       IF( isign > 0 ) idir = -1
       CALL zfftm1d( idir, nz, nsl, c(1), 1, ldc, tablez(1,ip) )
       IF (isign > 0) THEN
         tscale = 1.0d0 / nz
         CALL ZDSCAL( ldc * nsl, tscale, c(1), 1)
       END IF
     END IF

#elif defined __AIX

     IF( isign > 0 ) THEN
       tscale = 1.0d0 / nz
       idir   = 1
       CALL DCFT (0, c(1), 1, ldc, cout(1), 1, ldc, nz, nsl, idir, &
          tscale, fw_table(1, 3, ip), ltabl, work, lwork)
     ELSE IF( isign < 0 ) THEN
       idir   = -1
       tscale = 1.0d0
       CALL DCFT (0, c(1), 1, ldc, cout(1), 1, ldc, nz, nsl, idir, &
          tscale, bw_table(1, 3, ip), ltabl, work, lwork)
     END IF

#else 
                                                                                                      
     CALL errore(' cft_1 ',' no scalar fft driver specified ', 1)

#endif

#if ! defined __AIX
     cout( 1 : ldc * nsl ) = c( 1 : ldc * nsl )
#endif


     RETURN
   END SUBROUTINE cft_1z

!
!
!=----------------------------------------------------------------------=!
!
!
!
!         FFT along "x" and "y" direction
!
!
!
!=----------------------------------------------------------------------=!
!
!


   SUBROUTINE cft_2xy(r, nzl, nx, ny, ldx, ldy, sgn, pl2ix)

!     driver routine for nzl 2d complex fft's of lengths nx and ny
!     (sparse grid, both charge and wavefunctions) 
!     on input, sgn=+/-1 for charge density, sgn=+/-2 for wavefunctions
!     ldx is the actual dimension of f (may differ from n)
!     for compatibility: ldy is not used
!     A separate initialization is stored for each combination of input parameters

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: sgn, ldx, ldy, nx, ny, nzl
     INTEGER, OPTIONAL, INTENT(IN) :: pl2ix(:)
     COMPLEX (dbl) :: r( : )
     INTEGER :: i, k, j, idir, ip, isign, kk
     REAL(dbl) :: tscale
     INTEGER, SAVE :: icurrent = 1
     INTEGER, SAVE :: dims( 4, ndims) = -1
     LOGICAL :: dofft( nfftx )


     isign = - sgn

     dofft( 1 : nx ) = .TRUE.
     IF( PRESENT( pl2ix ) ) THEN
       IF( SIZE( pl2ix ) < nx ) &
         CALL errore( ' cft_2xy ', ' wrong dimension for arg no. 8 ', 1 )
       DO i = 1, nx
         IF( pl2ix(i) < 1 ) dofft( i ) = .FALSE.
       END DO
     END IF

     ! WRITE(6,*) 'DEBUG: ', COUNT( dofft )

     !
     !   Here initialize table only if necessary
     !

     ip = -1
     DO i = 1, ndims
            
       !   first check if there is already a table initialized
       !   for this combination of parameters

       IF( ( ny == dims(1,i) ) .and. ( ldx == dims(2,i) ) .and. &
           ( nx == dims(3,i) ) .and. ( nzl == dims(4,i) ) ) THEN
         ip = i
         EXIT
       END IF

     END DO

     IF( ip == -1 ) THEN

       !   no table exist for these parameters
       !   initialize a new one 

       ! WRITE(6, fmt="('DEBUG cft_2xy, reinitializing tables ', I3)" ) icurrent

#if defined __FFTW

       IF( fw_plan( 2, icurrent) /= 0 )   CALL DESTROY_PLAN_1D( fw_plan( 2, icurrent) )
       IF( bw_plan( 2, icurrent) /= 0 )   CALL DESTROY_PLAN_1D( bw_plan( 2, icurrent) )
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 2, icurrent), ny, idir)
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 2, icurrent), ny, idir)

       IF( fw_plan( 1, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( fw_plan( 1, icurrent) )
       IF( bw_plan( 1, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( bw_plan( 1, icurrent) )
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 1, icurrent), nx, idir) 
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 1, icurrent), nx, idir) 

#elif defined __AIX

       tscale = 1.0d0 / ( nx * ny )
       CALL DCFT ( 1, r(1), ldx, 1, r(1), ldx, 1, ny, 1,  1, 1.0d0, &
          fw_table( 1, 2, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, r(1), ldx, 1, r(1), ldx, 1, ny, 1, -1, 1.0d0, &
          bw_table(1, 2, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, r(1), 1, ldx, r(1), 1, ldx, nx, ny,  1, &
          tscale, fw_table( 1, 1, icurrent), ltabl, work(1), lwork)
       CALL DCFT ( 1, r(1), 1, ldx, r(1), 1, ldx, nx, ny, -1, &
          1.0d0, bw_table(1, 1, icurrent), ltabl, work(1), lwork)

#elif defined __SGI

       CALL ZFFT1DI( ny, tabley(1, icurrent) )
       CALL ZFFT1DI( nx, tablex(1, icurrent) )

#else

       CALL errore(' fft_y ',' no scalar fft driver specified ', 1)

#endif

       dims(1,icurrent) = ny; dims(2,icurrent) = ldx; 
       dims(3,icurrent) = nx; dims(4,icurrent) = nzl;
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1

     END IF

     !
     !   Now perform the FFTs using machine specific drivers
     !

#if defined __FFTW

     IF( isign > 0 ) THEN

       CALL FFT_X_STICK( fw_plan( 1, ip), r(1), nx, ny, nzl, ldx, ldy ) 
       do i = 1, nx
         do k = 1, nzl
           IF( dofft( i ) ) THEN
             j = i + ldx*ldy * ( k - 1 )
             call FFT_Y_STICK(fw_plan( 2, ip), r(j), ny, ldx) 
           END IF
         end do
       end do
       tscale = 1.0d0 / ( nx * ny )
       CALL ZDSCAL( ldx * ldy * nzl, tscale, r(1), 1)

     ELSE IF( isign < 0 ) THEN

       do i = 1, nx
         do k = 1, nzl
           IF( dofft( i ) ) THEN
             j = i + ldx*ldy * ( k - 1 )
             call FFT_Y_STICK( bw_plan( 2, ip), r(j), ny, ldx) 
           END IF
         end do
       end do
       CALL FFT_X_STICK( bw_plan( 1, ip), r(1), nx, ny, nzl, ldx, ldy ) 

     END IF

#elif defined __AIX

     IF( isign > 0 ) THEN

       idir = 1
       tscale = 1.0d0 / ( nx * ny )
       do k = 1, nzl
         kk = 1 + ( k - 1 ) * ldx * ldy
         CALL DCFT ( 0, r(kk), 1, ldx, r(kk), 1, ldx, nx, ny, idir, &
           tscale, fw_table( 1, 1, ip ), ltabl, work( 1 ), lwork)
         do i = 1, nx
           IF( dofft( i ) ) THEN
             kk = i + ( k - 1 ) * ldx * ldy
             call DCFT ( 0, r( kk ), ldx, 1, r( kk ), ldx, 1, ny, 1, &
               idir, 1.0d0, fw_table(1, 2, ip), ltabl, work( 1 ), lwork)
           END IF
         end do
       end do

     ELSE IF( isign < 0 ) THEN

       idir = -1
       DO k = 1, nzl
         do i = 1, nx
           IF( dofft( i ) ) THEN
             kk = i + ( k - 1 ) * ldx * ldy
             call DCFT ( 0, r( kk ), ldx, 1, r( kk ), ldx, 1, ny, 1, &
               idir, 1.0d0, bw_table(1, 2, ip), ltabl, work( 1 ), lwork)
           END IF
         end do
         kk = 1 + ( k - 1 ) * ldx * ldy
         CALL DCFT ( 0, r( kk ), 1, ldx, r( kk ), 1, ldx, nx, ny, idir, &
           1.0d0, bw_table(1, 1,  ip), ltabl, work( 1 ), lwork)
       END DO
         
     END IF

#elif defined __SGI

     IF( isign > 0 ) THEN
       idir =  -1
       DO i = 1, nzl
         k = 1 + ( i - 1 ) * ldx * ldy
         call zfftm1d( idir, nx, ny, r(k), 1, ldx, tablex(1,ip) )
       END DO
       do i = 1, nx
         IF( dofft( i ) ) THEN
           call zfftm1d( idir, ny, nzl, r(i), ldx, ldx*ldy, tabley(1, ip) )
         END IF
       end do
       tscale = 1.0d0 / ( nx * ny )
       CALL ZDSCAL( ldx * ldy * nzl, tscale, r(1), 1)
     ELSE IF( isign < 0 ) THEN
       idir = 1
       do i = 1, nx
         IF( dofft( i ) ) THEN
           call zfftm1d( idir, ny, nzl, r(i), ldx, ldx*ldy, tabley(1, ip) )
         END IF
       end do
       DO i = 1, nzl
         k = 1 + ( i - 1 ) * ldx * ldy
         call zfftm1d( idir, nx, ny, r(k), 1, ldx, tablex(1,ip) )
       END DO
     END IF

#else

     CALL errore(' cft_2xy ',' no scalar fft driver specified ', 1)

#endif


     return
   end subroutine cft_2xy


!
!=----------------------------------------------------------------------=!
!
!
!
!         3D scalar FFTs 
!
!
!
!=----------------------------------------------------------------------=!
!

   SUBROUTINE cfft3d( f, nr1, nr2, nr3, nr1x, nr2x, nr3x, sgn )

     IMPLICIT NONE

     INTEGER, INTENT(IN) :: nr1, nr2, nr3, nr1x, nr2x, nr3x, sgn 
     COMPLEX (dbl) :: f(:)
     INTEGER :: i, k, j, idir, ip, isign
     REAL(dbl) :: tscale
     INTEGER, SAVE :: icurrent = 1
     INTEGER, SAVE :: dims(3,ndims) = -1

#if defined __FFTW

     integer(kind=ipt), save :: fw_plan(ndims) = 0
     integer(kind=ipt), save :: bw_plan(ndims) = 0

#elif defined __AIX

#elif defined __SGI

      real(kind=8), save :: table( 3 * ltabl,  ndims )

#endif


     isign = -sgn

     !
     !   Here initialize table only if necessary
     !

     ip = -1
     DO i = 1, ndims

       !   first check if there is already a table initialized
       !   for this combination of parameters

       IF ( ( nr1 == dims(1,i) ) .and. ( nr2 == dims(2,i) ) .and. ( nr3 == dims(3,i) ) ) THEN
         ip = i
         EXIT
       END IF
     END DO

     IF( ip == -1 ) THEN

       !   no table exist for these parameters
       !   initialize a new one

#if defined __FFTW

       IF ( nr1 /= nr1x .or. nr2 /= nr2x .or. nr3 /= nr3x ) &
         call errore('cfft3','not implemented',1)

       IF( fw_plan(icurrent) /= 0 ) CALL DESTROY_PLAN_3D( fw_plan(icurrent) )
       IF( bw_plan(icurrent) /= 0 ) CALL DESTROY_PLAN_3D( bw_plan(icurrent) )
       idir = -1; CALL CREATE_PLAN_3D( fw_plan(icurrent), nr1, nr2, nr3, idir) 
       idir =  1; CALL CREATE_PLAN_3D( bw_plan(icurrent), nr1, nr2, nr3, idir) 

#elif defined __AIX

#elif defined __SGI

       CALL zfft3di( nr1, nr2, nr3, table(1,icurrent) )

#else

       CALL errore(' cfft3d ',' no scalar fft driver specified ', 1)

#endif

       dims(1,icurrent) = nr1; dims(2,icurrent) = nr2; dims(3,icurrent) = nr3
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1

     END IF

     !
     !   Now perform the 3D FFT using the machine specific driver
     !

#if defined __FFTW

     IF( isign > 0 ) THEN

       call FFTW_INPLACE_DRV_3D( fw_plan(ip), 1, f(1), 1, 1 )

       tscale = 1.0d0 / DBLE( nr1 * nr2 * nr3 )
       call ZDSCAL( nr1 * nr2 * nr3, tscale, f(1), 1)

     ELSE IF( isign < 0 ) THEN

       call FFTW_INPLACE_DRV_3D( bw_plan(ip), 1, f(1), 1, 1 )

     END IF

#elif defined __AIX

     if ( isign > 0 ) then
       tscale = 1.0d0 / ( nr1 * nr2 * nr3 )
     else
       tscale = 1.0d0
     end if
 
     call dcft3( f(1), nr1x, nr1x*nr2x, f(1), nr1x, nr1x*nr2x, nr1, nr2, nr3,  &
       isign, tscale, work(1), lwork)

#elif defined __SGI

     IF( isign > 0 ) idir = -1
     IF( isign < 0 ) idir = +1
     IF( isign /= 0 ) &
       CALL zfft3d( idir, nr1, nr2, nr3, f(1), nr1x, nr2x, table(1,ip) )
     IF( isign > 0 ) THEN
       tscale = 1.0d0 / DBLE( nr1 * nr2 * nr3 )
       call ZDSCAL( nr1x * nr2x * nr3x, tscale, f(1), 1)
     END IF

#endif

      
     RETURN
   END SUBROUTINE

!
!=----------------------------------------------------------------------=!
!
!
!
!         3D scalar FFTs,  but using sticks!
!
!
!
!=----------------------------------------------------------------------=!
!


!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------

subroutine cfft3ds (f, nr1, nr2, nr3, nrx1, nrx2, nrx3, sign, do_fft_x, do_fft_y)
  !
  !     driver routine for 3d complex "reduced" fft
  !     sign > 0 : f(G) => f(R)   ; sign < 0 : f(R) => f(G)
  !
  !     The 3D fft are computed only on lines and planes which have
  !     non zero elements. These lines and planes are defined by
  !     the two vectors do_fft_x and do_fft_y 
  !
  !     The routine is implemented for:
  !
  !       IBM      : essl library
  !
  !----------------------------------------------------------------------
  !
  implicit none

  integer :: nr1, nr2, nr3, nrx1, nrx2, nrx3, sign
  !
  !   logical dimensions of the fft
  !   physical dimensions of the f array
  !   sign of the transformation

  complex(dbl) :: f ( nrx1 * nrx2 * nrx3 )
  integer :: do_fft_x(:), do_fft_y(:)
  !
  ! the fft array
  !
  ! ESSL fft's require a different initialization for sign=-1 and sign=1
  ! aux1 contains the initialized quantities
  ! aux2 is work space
  !
  integer :: m, incx1, incx2
  INTEGER :: i, k, j, idir, ip, isign, ii, jj
  REAL(dbl) :: tscale
  INTEGER, SAVE :: icurrent = 1
  INTEGER, SAVE :: dims(3,ndims) = -1


  tscale = 1.d0
  isign = - sign   !  here we follow ESSL convention

  !
  ! ESSL sign convention for fft's is the opposite of the "usual" one
  !

  ! write(6, fmt="('DEBUG cfft3ds :',6I6)") nr1, nr2, nr3, nrx1, nrx2, nrx3
  ! write(6, fmt="('DEBUG cfft3ds :',24I2)") do_fft_x
  ! write(6, fmt="('DEBUG cfft3ds :',24I2)") do_fft_y

  IF( nr2 /= nrx2 ) &
    CALL errore(' cfft3ds ', ' wrong dimensions: nr2 /= nrx2 ', 1 )

     ip = -1
     DO i = 1, ndims

       !   first check if there is already a table initialized
       !   for this combination of parameters

       IF( ( nr1 == dims(1,i) ) .and. ( nr2 == dims(2,i) ) .and. &
           ( nr3 == dims(3,i) ) ) THEN
         ip = i
         EXIT
       END IF

     END DO

     IF( ip == -1 ) THEN

       !   no table exist for these parameters
       !   initialize a new one

#if defined __FFTW

       IF( fw_plan( 1, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( fw_plan( 1, icurrent) )
       IF( bw_plan( 1, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( bw_plan( 1, icurrent) )
       IF( fw_plan( 2, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( fw_plan( 2, icurrent) )
       IF( bw_plan( 2, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( bw_plan( 2, icurrent) )
       IF( fw_plan( 3, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( fw_plan( 3, icurrent) )
       IF( bw_plan( 3, icurrent) /= 0 ) CALL DESTROY_PLAN_1D( bw_plan( 3, icurrent) )
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 1, icurrent), nr1, idir) 
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 1, icurrent), nr1, idir) 
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 2, icurrent), nr2, idir) 
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 2, icurrent), nr2, idir) 
       idir = -1; CALL CREATE_PLAN_1D( fw_plan( 3, icurrent), nr3, idir) 
       idir =  1; CALL CREATE_PLAN_1D( bw_plan( 3, icurrent), nr3, idir) 

#elif defined __AIX

       tscale = 1.0d0 
       !  x - direction
       incx1 = 1; incx2 = nrx1; m = 1
       CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nr1, m,  1, 1.0d0, &
          fw_table( 1, 1, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nr1, m, -1, 1.0d0, &
          bw_table(1, 1, icurrent), ltabl, work(1), lwork )
       !  y - direction
       incx1 = nrx1; incx2 = 1; m = nr1;
       CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nr2, m,  1, 1.0d0, &
          fw_table( 1, 2, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nr2, m, -1, 1.0d0, &
          bw_table(1, 2, icurrent), ltabl, work(1), lwork )
       !  z - direction
       incx1 = nrx1 * nrx2; incx2 = 1; m = nrx1 * nr2
       CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nr3, m,  1, 1.0d0, &
          fw_table(1, 3, icurrent), ltabl, work(1), lwork )
       CALL DCFT ( 1, f(1), incx1, incx2, f(1), incx1, incx2, nr3, m, -1, 1.0d0, &
          bw_table(1, 3, icurrent), ltabl, work(1), lwork )

#else

       CALL errore(' cfft3ds ',' no scalar fft driver specified ', 1)

#endif

       dims(1,icurrent) = nr1; dims(2,icurrent) = nr2; dims(3,icurrent) = nr3
       ip = icurrent
       icurrent = MOD( icurrent, ndims ) + 1

     END IF


     IF ( isign < 0 ) THEN
   
        !
        !  i - direction ...
        !

        incx1 = 1;  incx2 = nrx1;  m = 1

        do k = 1, nr3
           do j = 1, nr2
              jj = j + ( k - 1 ) * nrx2 
              ii = 1 + nrx1 * ( jj - 1 ) 
              if ( do_fft_x( jj ) == 1 ) THEN
#if defined __FFTW
                call FFTW_INPLACE_DRV_1D( bw_plan( 1, ip), m, f( ii ), incx1, incx2 )
#elif defined __AIX
                call dcft (0, f ( ii ), incx1, incx2, f ( ii ), incx1, incx2, nr1, m, &
                  isign, 1.0d0, bw_table ( 1, 1,  ip ), ltabl, work( 1 ), lwork)
#else
                call errore(' cfft3ds ',' no scalar fft driver specified ', 1)
#endif
              endif
           enddo
        enddo

        !
        !  ... j-direction ...
        !

        incx1 = nrx1;  incx2 = 1;  m = nr1

        do k = 1, nr3
           ii = 1 + nrx1 * nrx2 * ( k - 1 ) 
           if ( do_fft_y( k ) == 1 ) then
#if defined __FFTW
             call FFTW_INPLACE_DRV_1D( bw_plan( 2, ip), m, f( ii ), incx1, incx2 )
#elif defined __AIX
             call dcft (0, f ( ii ), incx1, incx2, f ( ii ), incx1, incx2, nr2, m, &
               isign, 1.0d0, bw_table ( 1, 2,  ip ), ltabl, work( 1 ), lwork)
#else
             call errore(' cfft3ds ',' no scalar fft driver specified ', 1)
#endif
           endif
        enddo

        !
        !     ... k-direction
        !

        incx1 = nrx1 * nrx2;  incx2 = 1;  m = nrx1 * nr2

#if defined __FFTW
        call FFTW_INPLACE_DRV_1D( bw_plan( 3, ip), m, f( 1 ), incx1, incx2 )
#elif defined __AIX
        call dcft (0, f( 1 ), incx1, incx2, f( 1 ), incx1, incx2, nr3, m, &
          isign, 1.0d0, bw_table ( 1, 3, ip ), ltabl, work( 1 ), lwork)
#endif

     ELSE

        !
        !     ... k-direction
        !

        incx1 = nrx1 * nr2;  incx2 = 1;  m = nrx1 * nr2

#if defined __FFTW
        call FFTW_INPLACE_DRV_1D( fw_plan( 3, ip), m, f( 1 ), incx1, incx2 )
#elif defined __AIX
        call dcft (0, f( 1 ), incx1, incx2, f( 1 ), incx1, incx2, nr3, m, &
          isign, 1.0d0, fw_table ( 1, 3, ip ), ltabl, work( 1 ), lwork)
#endif

        !
        !     ... j-direction ...
        !

        incx1 = nrx1;  incx2 = 1;  m = nr1

        do k = 1, nr3
           ii = 1 + nrx1 * nrx2 * ( k - 1 ) 
           if ( do_fft_y ( k ) == 1 ) then
#if defined __FFTW
             call FFTW_INPLACE_DRV_1D( fw_plan( 2, ip), m, f( ii ), incx1, incx2 )
#elif defined __AIX
             call dcft (0, f ( ii ), incx1, incx2, f ( ii ), incx1, incx2, nr2, m, &
               isign, 1.0d0, fw_table ( 1, 2, ip ), ltabl, work( 1 ), lwork)
#else
             call errore(' cfft3ds ',' no scalar fft driver specified ', 1)
#endif
           endif
        enddo

        !
        !     i - direction ...
        !

        incx1 = 1;  incx2 = nrx1;  m = 1

        do k = 1, nr3
           do j = 1, nr2
              jj = j + ( k - 1 ) * nrx2 
              ii = 1 + nrx1 * ( jj - 1 ) 
              if ( do_fft_x( jj ) == 1 ) then
#if defined __FFTW
                call FFTW_INPLACE_DRV_1D( fw_plan( 1, ip), m, f( ii ), incx1, incx2 )
#elif defined __AIX
                call dcft (0, f ( ii ), incx1, incx2, f ( ii ), incx1, incx2, nr1, m, &
                   isign, 1.0d0, fw_table ( 1, 1, ip ), ltabl, work( 1 ), lwork)
#else
                call errore(' cfft3ds ',' no scalar fft driver specified ', 1)
#endif
              endif
           enddo
        enddo

#if defined __AIX || defined __FFTW
        call DSCAL (2 * nrx1 * nrx2 * nr3, 1.0d0 / (nr1 * nr2 * nr3), f( 1 ), 1)
#endif

     END IF

     RETURN
   END SUBROUTINE 

!
!=----------------------------------------------------------------------=!
!
!
!
!         3D parallel FFT on sub-grids
!
!
!
!=----------------------------------------------------------------------=!
!
   SUBROUTINE cft_b ( f, n1, n2, n3, n1x, n2x, n3x, imin3, imax3, sgn )

!     driver routine for 3d complex fft's on box grid - ibm essl
!     fft along xy is done only on planes that correspond to
!     dense grid planes on the current processor, i.e. planes
!     with imin3 .le. n3 .le. imax3
!
      implicit none
      integer :: k,n1,n2,n3,n1x,n2x,n3x,imin3,imax3,sgn
      complex(kind=8) :: f(:)

      integer isign, nplanes, nstart
      real(dbl) :: tscale

      integer :: ip, i
      integer, save :: icurrent = 1
      integer, save :: dims( 4, ndims ) = -1

#if defined __FFTW

      integer(kind=ipt), save :: bw_planz(  ndims ) = 0
      integer(kind=ipt), save :: bw_planxy( ndims ) = 0

#elif defined __AIX

      real(kind=8), save :: aux3( ltabl, ndims )
      real(kind=8), save :: aux2( ltabl, ndims )
      real(kind=8), save :: aux1( ltabl, ndims )

#elif defined __SGI

      real(kind=8), save :: bw_coeffz( ltabl,  ndims )
      real(kind=8), save :: bw_coeffy( ltabl,  ndims )
      real(kind=8), save :: bw_coeffx( ltabl,  ndims )

#endif


      isign = -sgn
      tscale = 1.d0

      if ( isign > 0 ) then
         call errore('cft_b','not implemented',isign)
      end if
!
! 2d fft on xy planes - only needed planes are transformed
! note that all others are left in an unusable state
!
      nplanes = imax3 - imin3 + 1
      nstart  = ( imin3 - 1 ) * n1x * n2x + 1

      !
      !   Here initialize table only if necessary
      !

      ip = -1
      DO i = 1, ndims

        !   first check if there is already a table initialized
        !   for this combination of parameters

        IF ( ( n1 == dims(1,i) ) .and. ( n2 == dims(2,i) ) .and. &
             ( n3 == dims(3,i) ) .and. ( nplanes == dims(4,i) ) ) THEN
           ip = i
           EXIT
        END IF

      END DO

      IF( ip == -1 ) THEN

        !   no table exist for these parameters
        !   initialize a new one

#if defined __FFTW

        if ( bw_planz(icurrent) /= 0 ) call DESTROY_PLAN_1D( bw_planz(icurrent) )
        call CREATE_PLAN_1D( bw_planz(icurrent), n3, 1 )

        if ( bw_planxy(icurrent) /= 0 ) call DESTROY_PLAN_2D( bw_planxy(icurrent) )
        call CREATE_PLAN_2D( bw_planxy(icurrent), n1, n2, 1 )
!
#elif defined __AIX

         if( n3 /= dims(3,icurrent) ) then
           call dcft( 1, f(1), n1x*n2x, 1, f(1), n1x*n2x, 1, n3, n1x*n2x, isign,          &
     &        tscale, aux3(1,icurrent), ltabl, work(1), lwork)
         end if
         call dcft( 1, f(1), 1, n1x, f(1), 1, n1x, n1, n2x*nplanes, isign,              &
     &        tscale, aux1(1,icurrent), ltabl, work(1), lwork)
         if( n2 /= dims(2,icurrent) ) then
           call dcft( 1, f(1), n1x, 1, f(1), n1x, 1, n2, n1x, isign,                      &
     &        tscale, aux2(1,icurrent), ltabl, work(1), lwork)
         end if

#elif defined __SGI

         call zfft1di( n3, bw_coeffz( 1, icurrent ) )
         call zfft1di( n2, bw_coeffy( 1, icurrent ) )
         call zfft1di( n1, bw_coeffx( 1, icurrent ) )

#else

        CALL errore(' cft_b ',' no scalar fft driver specified ', 1)
 

#endif

        dims(1,icurrent) = n1; dims(2,icurrent) = n2
        dims(3,icurrent) = n3; dims(4,icurrent) = nplanes
        ip = icurrent
        icurrent = MOD( icurrent, ndims ) + 1

      END IF


#if defined __FFTW

      call FFTW_INPLACE_DRV_1D( bw_planz(ip), n1x*n2x, f(1), n1x*n2x, 1 )
      call FFTW_INPLACE_DRV_2D( bw_planxy(ip), nplanes, f(nstart), 1, n1x*n2x )

#elif defined __AIX


      !   fft in the z-direction...

      call dcft( 0, f(1), n1x*n2x, 1, f(1), n1x*n2x, 1, n3, n1x*n2x, isign,             &
     &        tscale, aux3(1,ip), ltabl, work(1), lwork)

      !   x-direction

      call dcft( 0, f(nstart), 1, n1x, f(nstart), 1, n1x, n1, n2x*nplanes, isign,  &
     &        tscale, aux1(1,ip), ltabl, work(1), lwork)
     
      !   y-direction
     
      DO k = imin3, imax3
        nstart = ( k - 1 ) * n1x * n2x + 1
        call dcft( 0, f(nstart), n1x, 1, f(nstart), n1x, 1, n2, n1x, isign,        &
     &        tscale, aux2(1,ip), ltabl, work(1), lwork)
      END DO

#elif defined __SGI

      call zfftm1d( 1, n3, n1x*n2x, f(1), n1x*n2x, 1, bw_coeffz(1, ip) )
      call zfftm1d( 1, n1, n2x*nplanes, f(nstart), 1, n1x, bw_coeffx(1, ip) )
      DO K = imin3, imax3
        nstart = ( k - 1 ) * n1x * n2x + 1
        call zfftm1d( 1, n2, n1x, f(nstart), n1x, 1, bw_coeffy(1, ip) )
      END DO

#endif

     RETURN
   END SUBROUTINE

!
!=----------------------------------------------------------------------=!
!
!
!
!         FFT support Functions/Subroutines
!
!
!
!=----------------------------------------------------------------------=!
!

!
integer function good_fft_dimension (n)
  !
  ! Determines the optimal maximum dimensions of fft arrays
  ! Useful on some machines to avoid memory conflicts
  !
  use parameters
  implicit none
  integer :: n, nx
  ! this is the default: max dimension = fft dimension
  nx = n
#if defined(__AIX) || defined(DXML)
  if ( n==8 .or. n==16 .or. n==32 .or. n==64 .or. n==128 .or. n==256) &
       nx = n + 1
#endif
#if defined(CRAYY) || defined(__SX4)
  if (mod (nr1, 2) ==0) nx = n + 1
#endif
  good_fft_dimension = nx
  return
end function good_fft_dimension

!=----------------------------------------------------------------------=!

function allowed (nr)


  ! find if the fft dimension is a good one
  ! a "bad one" is either not implemented (as on IBM with ESSL)
  ! or implemented but with awful performances (most other cases)

  use parameters

  implicit none
  integer :: nr

  logical :: allowed
  integer :: pwr (5)
  integer :: mr, i, fac, p, maxpwr
  integer :: factors( 5 ) = (/ 2, 3, 5, 7, 11 /)

  ! find the factors of the fft dimension

  mr  = nr
  pwr = 0
  factors_loop: do i = 1, 5
     fac = factors (i)
     maxpwr = NINT ( LOG( REAL (mr) ) / LOG( REAL (fac) ) ) + 1
     do p = 1, maxpwr
        if ( mr == 1 ) EXIT factors_loop
        if ( MOD (mr, fac) == 0 ) then
           mr = mr / fac
           pwr (i) = pwr (i) + 1
        endif
     enddo
  end do factors_loop

  IF ( nr /= ( mr * 2**pwr (1) * 3**pwr (2) * 5**pwr (3) * 7**pwr (4) * 11**pwr (5) ) ) &
     CALL errore (' allowed ', ' what ?!? ', 1 )

  if ( mr /= 1 ) then

     ! fft dimension contains factors > 11 : no good in any case

     allowed = .false.

  else

#ifdef __AIX

     ! IBM machines with essl libraries

     allowed =  ( pwr(1) >= 1 ) .and. ( pwr(2) <= 2 ) .and. ( pwr(3) <= 1 ) .and. &
                ( pwr(4) <= 1 ) .and. ( pwr(5) <= 1 ) .and. &
                ( ( (pwr(2) == 0 ) .and. ( pwr(3) + pwr(4) + pwr(5) ) <= 2 ) .or. &
                  ( (pwr(2) /= 0 ) .and. ( pwr(3) + pwr(4) + pwr(5) ) <= 1 ) )
#else

     ! fftw and all other cases: no factors 7 and 11

     allowed = ( ( pwr(4) == 0 ) .and. ( pwr(5) == 0 ) )

#endif

  endif

  return
end function allowed

!=----------------------------------------------------------------------=!

   INTEGER FUNCTION good_fft_order( nr, np )

!    
!    This function find a "good" fft order value grather or equal to "nr"
!
!    nr  (input) tentative order n of a fft
!            
!    np  (optional input) if present restrict the search of the order
!        in the ensamble of multiples of np
!            
!    Output: the same if n is a good number
!         the closest higher number that is good
!         an fft order is not good if not implemented (as on IBM with ESSL)
!         or implemented but with awful performances (most other cases)
!

     IMPLICIT NONE
     INTEGER, INTENT(IN) :: nr
     INTEGER, OPTIONAL, INTENT(IN) :: np
     INTEGER :: new

     new = nr
     IF( PRESENT( np ) ) THEN
       DO WHILE( ( ( .NOT. allowed( new ) ) .OR. ( MOD( new, np ) /= 0 ) ) .AND. ( new <= nfftx ) )
         new = new + 1
       END DO
     ELSE
       DO WHILE( ( .NOT. allowed( new ) ) .AND. ( new <= nfftx ) )
         new = new + 1
       END DO
     END IF

     IF( new > nfftx ) &
       CALL errore( ' good_fft_order ', ' fft order too large ', new )

     good_fft_order = new
  
     RETURN
   END FUNCTION


!=----------------------------------------------------------------------=!
   END MODULE fft_scalar
!=----------------------------------------------------------------------=!
