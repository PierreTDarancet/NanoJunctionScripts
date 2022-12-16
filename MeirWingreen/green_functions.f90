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
   MODULE T_green_operation_module
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
   REAL(dbl), PARAMETER :: delta_lead =  EPS_m3
   ! Public variables

   ! Public routines:
   PUBLIC                 :: calcul_g_C
   PUBLIC                 :: calcul_g_lesser
   PUBLIC                 :: calcul_spectralfunction
   CONTAINS 

!***********************************************
   SUBROUTINE calcul_g_C( gC , h00, sigma_l , sigma_r, ene, dimC)
   !***********************************************
  IMPLICIT NONE
      INTEGER, INTENT(in)   :: dimC
      COMPLEX(dbl), INTENT(out) ::  gC(dimC,dimC) 
      COMPLEX(dbl), INTENT(in) ::  h00(dimC,dimC) 
      COMPLEX(dbl), INTENT(in) ::  sigma_l(dimC,dimC) 
      COMPLEX(dbl), INTENT(in) ::  sigma_r(dimC,dimC) 
      COMPLEX(dbl), INTENT(in) ::  ene

      ! Local
       CHARACTER(10)      :: subname="calcul_g_C"
       INTEGER :: ierr, idiag
       COMPLEX(dbl) ::  work1(dimC,dimC)
       COMPLEX(dbl) ::  work2(dimC,dimC)
         !
         IF( dimC < 1 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in Conductor dimension", " dim C = ", dimC
            STOP
         ENDIF 
         !
         work1(:,:)= CZERO
         work2(:,:)= CZERO
         gC(:,:)= CZERO
         !
         DO idiag = 1, dimC
             work1(idiag,idiag)= CONE * ene
             gC(idiag,idiag)= CONE 
         ENDDO
         !
         work2(:,:) = work1(:,:) - h00(:,:) -  sigma_l(:,:) - sigma_r(:,:)
         !
         CALL mat_sv(dimC,dimC, work2, gC)
         !
    END SUBROUTINE calcul_g_C

	!***********************************************
   SUBROUTINE calcul_spectralfunction( A_C, gC , h00, sigma_l , sigma_r, sigma_c, ene, dimC)
   !***********************************************
	!  calcul_spectralfunction( A_C(ie,:,:) , g_C(ie,:,:), h_00(:,:), sigma_L(ie,:,:), sigma_R(ie,:,:) , sigma_C_r(ie,:,:) , ene, dimC)
  IMPLICIT NONE
      INTEGER, INTENT(in)   :: dimC
      REAL(dbl), INTENT(out) ::  A_C(dimC,dimC) 
      COMPLEX(dbl), INTENT(out) ::  gC(dimC,dimC) 
      COMPLEX(dbl), INTENT(in) ::  h00(dimC,dimC) 
      COMPLEX(dbl), INTENT(in) ::  sigma_l(dimC,dimC) 
      COMPLEX(dbl), INTENT(in) ::  sigma_r(dimC,dimC) 
      COMPLEX(dbl), INTENT(in) ::  sigma_c(dimC,dimC) 
      COMPLEX(dbl), INTENT(in) ::  ene

      ! Local
       CHARACTER(10)      :: subname="calcul_g_C"
       INTEGER :: ierr, idiag
       COMPLEX(dbl) ::  work1(dimC,dimC)
       COMPLEX(dbl) ::  work2(dimC,dimC)
         !
         IF( dimC < 1 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in Conductor dimension", " dim C = ", dimC
            STOP
         ENDIF 
         !
         work1(:,:)= CZERO
         work2(:,:)= CZERO
         gC(:,:)= CZERO
         A_C(:,:)= ZERO
         !
         DO idiag = 1, dimC
             work1(idiag,idiag)= CONE * ene
             gC(idiag,idiag)= CONE 
         ENDDO
         !
         work2(:,:) = work1(:,:) - h00(:,:) -  sigma_l(:,:) - sigma_r(:,:) - sigma_c(:,:)
         !
         CALL mat_sv(dimC,dimC, work2, gC)
         !
         A_C(:,:) = CI * (  gC(:,:) - CONJG( TRANSPOSE(gC(:,:)) )   )

    END SUBROUTINE calcul_spectralfunction
	              


!***********************************************
   SUBROUTINE calcul_g_lesser( glesser , sigma_l , sigma_r , sigma_c , gC , dimC)
   !***********************************************
   !calcul_g_lesser( glesser(ie,:,:) , sigma_L_lesser(ie,:,:), sigma_R_lesser(ie,:,:) , sigma_C_lesser(ie,:,:) , g_C(ie,:,:), dimC)
  IMPLICIT NONE
      INTEGER, INTENT(in)   :: dimC
      COMPLEX(dbl), INTENT(in) ::  gC(dimC,dimC) 
      COMPLEX(dbl), INTENT(in) ::  sigma_l(dimC,dimC) 
      COMPLEX(dbl), INTENT(in) ::  sigma_r(dimC,dimC) 
      COMPLEX(dbl), INTENT(in) ::  sigma_c(dimC,dimC) 
      COMPLEX(dbl), INTENT(out) ::  glesser(dimC,dimC) 

      ! Local
       CHARACTER(15)      :: subname="calcul_g_lesser"
       INTEGER :: ierr, irow, icol
       COMPLEX(dbl) ::  work(dimC,dimC)
       COMPLEX(dbl) ::  work2(dimC,dimC)
         !
         IF( dimC < 1 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in Conductor dimension", " dim C = ", dimC
            STOP
         ENDIF 
         !
         work(:,:)= CZERO
         DO icol=1, dimC
            DO irow=1, dimC
               work2( irow,icol)= sigma_l(irow,icol) +  sigma_r(irow,icol) + sigma_c(irow,icol) 
            ENDDO
         ENDDO
         !  gintr * sigma_lesser -> work
         ! 
         CALL mat_mul(work, gC, 'N', work2, 'N', dimC, dimC, dimC)
         !
         ! gintr * sigma_lesser * ginta  -->  glesser
         !
         CALL mat_mul( glesser, work, 'N', gC, 'C', dimC, dimC, dimC)
         !
         !
    END SUBROUTINE calcul_g_lesser
!***********************************************
   SUBROUTINE calcul_spectralfunction( A_C , g_C, h_00, sigma_l , sigma_r , sigma_c ,  ene, dimC)
   !***********************************************
  IMPLICIT NONE
      INTEGER, INTENT(in)   :: dimC
      COMPLEX(dbl), INTENT(in) ::  h_00(dimC,dimC) 
      COMPLEX(dbl), INTENT(in) ::  sigma_l(dimC,dimC) 
      COMPLEX(dbl), INTENT(in) ::  sigma_r(dimC,dimC) 
      COMPLEX(dbl), INTENT(in) ::  sigma_c(dimC,dimC) 
      COMPLEX(dbl), INTENT(in) ::  ene
      REAL(dbl), INTENT(out) ::  A_C(dimC,dimC) 
      COMPLEX(dbl), INTENT(out) ::  g_C(dimC,dimC) 
      ! Local
       CHARACTER(23)      :: subname="calcul_spectralfunction"
       INTEGER :: ierr, idiag
       COMPLEX(dbl) ::  work1(dimC,dimC)
       COMPLEX(dbl) ::  work2(dimC,dimC)

         !
         IF( dimC < 1 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in Conductor dimension", " dim C = ", dimC
            STOP
         ENDIF 
         !

         work1(:,:)= CZERO
         work2(:,:)= CZERO
         g_C(:,:)= CZERO
         A_C(:,:)= CZERO
         !
         DO idiag = 1, dimC
             work1(idiag,idiag)= CONE * ene
             gC(idiag,idiag)= CONE 
         ENDDO
         !
         work2(:,:) = work1(:,:) - h00(:,:) -  sigma_l(:,:) - sigma_r(:,:) - sigma_c(:,:)
         !
         ! g_C^r
         CALL mat_sv(dimC, dimC, work2, g_C)
    

          A_C(:,:) = CI * (  g_C(:,:) - CONJG( TRANSPOSE( g_C(:,:)) )   )

    END SUBROUTINE calcul_spectralfunction

  END  MODULE T_green_operation_module

 

 


