!
! Copyright (C) 2004 PWSCF-CP-FPMD group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
#include "machine.h"
!
subroutine invmat (n, a, a_inv, da)
  !-----------------------------------------------------------------------
  ! computes the inverse "a_inv" of matrix "a", both dimensioned (n,n)
  ! if the matrix is dimensioned 3x3, it also computes determinant "da"
  ! matrix "a" is unchanged on output - LAPACK
  !
  USE kinds, ONLY : dbl
  implicit none
  integer :: n
  real(kind=dbl), DIMENSION (n,n) :: a, a_inv
  real(kind=dbl) :: da
  !
  integer :: info, lda, lwork, ipiv (n)
  ! info=0: inversion was successful
  ! lda   : leading dimension (the same as n)
  ! ipiv  : work space for pivoting (assumed of length lwork=n)
  real(kind=dbl) :: work (n) 
  ! more work space
  !
  lda = n
  lwork=n
  !
  a_inv(:,:) = a(:,:)
  !
  CALL DGETRF (n, n, a_inv, lda, ipiv, info)
     IF (info/=0) CALL errore ('invmat', 'error in DGETRF', abs (info) )
  CALL DGETRI (n, a_inv, lda, ipiv, work, lwork, info)
     IF (info/=0) CALL errore ('invmat', 'error in DGETRI', abs (info) )
  !
  if (n == 3) then
     da = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) + &
          a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) + &
          a(1,3)*(a(2,1)*a(3,2)-a(3,1)*a(2,2))
     IF (ABS(da) < 1.d-10) CALL errore(' invmat ',' singular matrix ', 1)
  else
     da = 0.d0
  end if

  return
end subroutine invmat
