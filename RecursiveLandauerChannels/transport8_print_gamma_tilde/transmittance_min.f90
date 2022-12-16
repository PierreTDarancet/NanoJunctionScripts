!
!      Copyright (C) 2005 WanT Group
!
!      This file is distributed under the terms of the
!      GNU General Public License. See the file `License'
!      in the root directory of the present distribution,
!      or http://www.gnu.org/copyleft/gpl.txt .
!
!***********************************************
   SUBROUTINE transmittance_min(dimC, gL, gR, gintr1N, conduct)
   !***********************************************
   !dim_subspace, gamma_tilde_L, gamma_tilde_R, P1_g_CC_PN(:,:), PN_g_CC_P1(:,:), conduct(:,ie)
   ! Calculates the matrix involved in the quantum transmittance, 
   ! returned in CONDUCT variable, following
   ! LANDAUER formula in the MEAN FILED case
   ! otherwise uses another formula derived in the
   ! paper PRL 94, 
   !
   USE kinds,          ONLY : dbl
   USE constants,      ONLY : CZERO, CONE, CI, ZERO , EPS_m5
   USE util_module,    ONLY : mat_mul
   USE timing_module,  ONLY : timing
   IMPLICIT NONE

   !
   ! input/output variables
   !
   INTEGER,      INTENT(in) ::  dimC
   COMPLEX(dbl), INTENT(in) ::  gL(dimC,dimC), gR(dimC,dimC)
   COMPLEX(dbl), INTENT(in) ::  gintr1N(dimC,dimC)
   REAL(dbl),    INTENT(out)::  conduct(dimC)

   !
   ! local variables
   !
   COMPLEX(dbl), ALLOCATABLE :: tmp(:,:), tmp1(:,:), ga_aux(:,:)
   INTEGER :: i, j, ierr
   !
   ! end of declarations
   !

!
!------------------------------
! main body
!------------------------------
!
   CALL timing('transmittance_min', OPR='start')

   ALLOCATE( tmp(dimC,dimC), tmp1(dimC,dimC), STAT=ierr )
      IF (ierr/=0) CALL errore('transmittance','allocating tmp,tm1',ABS(ierr))
   ALLOCATE( ga_aux(dimC,dimC), STAT=ierr )
      IF (ierr/=0) CALL errore('transmittance','allocating ga_aux',ABS(ierr))

  ! ordinary formula
  !
    ! 
   ! adding the identity matrix
   ! 


   !
   ! gL * gr -> tmp1
   !
   CALL mat_mul(tmp1, gL, 'N', gintr1N, 'N', dimC, dimC, dimC)
   !
   ! gL * gr * gR -> tmp
   !
   CALL mat_mul(tmp, tmp1, 'N', gR, 'N', dimC, dimC, dimC)
   !
   ! gL * gr * gR * ga -> tmp
   !
   ga_aux(:,:) =  CONJG( TRANSPOSE( gintr1N(:,:) ))
   CALL mat_mul(tmp1, tmp, 'N', ga_aux, 'N', dimC, dimC, dimC)
       
      
   DO i=1,dimC
      conduct(i) = REAL( tmp1(i,i) )
   ENDDO

!
! local memopry clean
!
   DEALLOCATE( tmp, tmp1, STAT=ierr )
      IF (ierr/=0) CALL errore('transmittance','deallocating tmp,tm1',ABS(ierr))
   DEALLOCATE( ga_aux, STAT=ierr )
      IF (ierr/=0) CALL errore('transmittance','deallocating ga_aux',ABS(ierr))
      
   CALL timing('transmittance_min', OPR='stop')
END SUBROUTINE transmittance_min

