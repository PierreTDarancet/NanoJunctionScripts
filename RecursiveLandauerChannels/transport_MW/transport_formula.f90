!
! Copyright (C) 2009 Molecular Foundry Berkeley
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!***********************************************
   MODULE  transport_formula_module 
   !***********************************************
  USE kinds, ONLY : dbl
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI, EPS_m3, EPS_m5
  USE util_module,          ONLY : mat_mul, mat_sv
   IMPLICIT NONE
   PRIVATE 
   SAVE
!


   ! Public
   ! Private
   ! Public variables

   ! Public routines:

   PUBLIC                 :: transmittance_calculation 
   PUBLIC                 :: transmittance_mw_calculation		 
   CONTAINS 
!***********************************************
   SUBROUTINE  transmittance_calculation(conductance, gC, gamR, gamL, dimx)
   !***********************************************
   !                       CALL  transmittance_calculation( cond_aux,  g_C(grid_condtogreen(ie),:,:),   gamma_R(grid_condtogreen(ie),:,:), gamma_L(grid_condtogreen(ie),:,:), dimC)
  IMPLICIT NONE
      INTEGER, INTENT(in)   :: dimx
      REAL(dbl), INTENT(out) ::  conductance

      COMPLEX(dbl), INTENT(in) ::  gC(dimx,dimx) 
      COMPLEX(dbl), INTENT(in) ::  gamR(dimx,dimx) 
      COMPLEX(dbl), INTENT(in) ::  gamL(dimx,dimx) 

      ! Local
       CHARACTER(25)      :: subname="transmittance_calculation"
       INTEGER :: ierr, idiag
       COMPLEX(dbl)           :: tmp(dimx,dimx), tmp1(dimx,dimx)
       REAL(dbl) ::  conduct(dimx)

        IF( dimx < 1 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in Conductor dimension", " dimC  ", dimx
            STOP
         ENDIF 

   !
   ! gL * gintr -> tmp
   ! 
   CALL mat_mul(tmp, gamL, 'N', gC, 'N', dimx, dimx, dimx)
   !
   ! gL * gintr * gR -> tmp1
   !
   CALL mat_mul(tmp1, tmp, 'N', gamR, 'N', dimx, dimx, dimx)
   !
   ! gL * gintr * gR * ginta  -->  tmp
   !
   CALL mat_mul(tmp, tmp1, 'N', gC, 'C', dimx, dimx, dimx)
   !
       
   conductance = ZERO
   DO idiag=1,dimx
      conduct(idiag) = REAL( tmp(idiag,idiag) )
   ENDDO
   conductance = SUM(conduct(:))
       !
    END SUBROUTINE transmittance_calculation

!***********************************************
   SUBROUTINE  transmittance_mw_calculation(conductance, conduct,  f_L, f_R, A_C, gC, gamR, gamL, dimx)
   !***********************************************
   !  transmittance_mw_calculation( cond_aux, transmittance_perchannel(ik,ie,:) ,  f_L(grid_condtogreen(ie)), f_R(grid_condtogreen(ie)), A_C(grid_condtogreen(ie),:,:) , g_lesser(grid_condtogreen(ie),:.:), gamma_R(grid_condtogreen(ie),:,:), gamma_L(grid_condtogreen(ie),:,:), dimC)
  IMPLICIT NONE
      INTEGER, INTENT(in)   :: dimx
      REAL(dbl), INTENT(out) ::  conductance
      REAL(dbl), INTENT(out) ::  conduct(dimx) !cond per channel 
      REAL(dbl), INTENT(in)  ::  f_L
      REAL(dbl), INTENT(in)  ::  f_R
      COMPLEX(dbl), INTENT(in) ::  gC(dimx,dimx) !Glesser
      REAL(dbl), INTENT(in)    ::  A_C(dimx,dimx)  !spectral function
      COMPLEX(dbl), INTENT(in) ::  gamR(dimx,dimx) 
      COMPLEX(dbl), INTENT(in) ::  gamL(dimx,dimx) 

      ! Local
       CHARACTER(25)      :: subname="transmittance_calculation"
       INTEGER :: ierr, idiag, irow, icol
       COMPLEX(dbl)           :: tmp(dimx,dimx), tmp1(dimx,dimx),  tmp2(dimx,dimx), tmp3(dimx,dimx), work_AC(dimx,dimx)

        IF( dimx < 1 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in Conductor dimension", " dimC  ", dimx
            STOP
         ENDIF 

    work_AC(:,:) = A_C(:,:)

   DO icol=1, dimx
      DO irow=1, dimx
         tmp2(irow,icol)= f_L * gamL(irow,icol) - f_R * gamR(irow,icol)
         tmp3(irow,icol)= gamL(irow,icol) - gamR(irow,icol)
      ENDDO
   ENDDO

   !
   ! (f_L \gamma_L - f_R \gamma_R ) * A_C -> tmp
   ! 
   CALL mat_mul(tmp, tmp2, 'N', work_AC, 'N', dimx, dimx, dimx)
   !
   !  ( \gamma_L -  \gamma_R ) * gC -> tmp1
   !
   CALL mat_mul(tmp1, tmp3, 'N', gC, 'N', dimx, dimx, dimx)

   DO icol=1, dimx
      DO irow=1, dimx
         tmp2(irow,icol)= tmp(irow,icol) - CI * tmp1(irow,icol)
      ENDDO
   ENDDO


       
   conductance = ZERO
   DO idiag=1,dimx
      conduct(idiag) = REAL( tmp2(idiag,idiag) )
   ENDDO
   conductance = SUM(conduct(:))
       !
    END SUBROUTINE transmittance_mw_calculation
  END  MODULE transport_formula_module 

