!
! Copyright (C) 2006 LEPES-CNRS/2007 Institut Néel Grenoble
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!***********************************************************************************
   SUBROUTINE green_brutal( dim_mat, cu, eig, ene, g)
   !***********************************************************************************
   !
   !  Construct green's functions
   !     
   !
   USE kinds
   USE constants, ONLY : CZERO, CONE
   USE timing_module, ONLY : timing
   USE util_module,      ONLY : mat_mul, mat_sv
   IMPLICIT NONE

   !
   ! I/O variables
   !
   INTEGER,      INTENT(in)    :: dim_mat
   COMPLEX(dbl), INTENT(in)    :: cu(dim_mat,dim_mat)
   REAL(dbl), INTENT(in)       :: eig(dim_mat)
   COMPLEX(dbl), INTENT(in)    :: ene
   COMPLEX(dbl), INTENT(out)   :: g(dim_mat,dim_mat)

   !
   ! local variables
   !
   INTEGER                   :: i, ierr
   COMPLEX(dbl)              :: work(dim_mat,dim_mat)
   COMPLEX(dbl)              :: g_work(dim_mat,dim_mat)
   COMPLEX(dbl)              :: cu_dagger(dim_mat,dim_mat)

!
!----------------------------------------
! main Body
!----------------------------------------
!
   CALL timing('green_brutal',OPR='start')
!  local variables
   cu_dagger(:,:) = CONJG( TRANSPOSE( cu(:,:)  ))
   g_work(:,:) = CZERO
   work(:,:) = CZERO
   g(:,:)=CZERO
!  output variable


!  initialisation ene-h(i,i)

   DO i=1,dim_mat
      !
         g_work(i,i)= CONE / (ene - eig(i))
      !
   ENDDO
!



! Rotation

    work(:,:) = CZERO

    work(:,:) = matmul(cu, g_work)
    g(:,:) = matmul(work, cu_dagger )

!

   CALL timing('green_brutal',OPR='stop')

END SUBROUTINE green_brutal

