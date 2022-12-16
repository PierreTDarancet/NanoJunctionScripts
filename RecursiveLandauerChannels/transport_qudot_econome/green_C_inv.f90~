!
! Copyright (C) 2007 Institut Néel Grenoble
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!*********************************************
   SUBROUTINE  green_CC_inv( g_CC, dos, ene, H_matrix_in, sgm_L, sgm_R, dimC) 
!*********************************************
   USE parameters,      ONLY : nstrx
   USE kinds
   USE constants,            ONLY : ZERO, CZERO, CONE, ONE, EPS_m6, EPS_m4, EPS_m2, CI, PI
   USE util_module,          ONLY : mat_sv
   USE timing_module, ONLY : timing
   IMPLICIT NONE

  
! Le but est d'obtenir l'element de matrice P1 G_CC (z) PN
! les 1 et N sont censés correspondre aux indices des matrices A et B
! Sgm_L/R est la self energy des contacts L et R et est donc projetee sur P1/PN
! ATTENTION les g_diag et g_off_diag NE SONT dans l'ensemble PAS les elements de matrice
! de la fonction Green TOTALE MAIS des elements de matrices de la fonction de Green 
! d'une restriction du Hamiltonien


   !
   ! Input variables
   !
   INTEGER,      INTENT(in)    :: dimC
   !
   !
   COMPLEX(dbl), INTENT(in)    :: H_matrix_in(dimC,dimC)
   COMPLEX(dbl), INTENT(in)    :: sgm_L(dimC,dimC)
   COMPLEX(dbl), INTENT(in)    :: sgm_R(dimC,dimC)
   COMPLEX(dbl), INTENT(in)    :: ene
   !
   ! Output variables
   !
   COMPLEX(dbl), INTENT(out)   :: g_CC(dimC,dimC)
   REAL(dbl),    INTENT(out)   :: dos(dimC)
   !
   ! local variables
   !
   INTEGER                   :: i_diag
   COMPLEX(dbl)              :: i_work(dimC,dimC)
   COMPLEX(dbl)              :: h_work(dimC,dimC)
   COMPLEX(dbl)              :: g_aux(dimC,dimC)

   CHARACTER(12)             :: subname='green_CC_inv'
  
!
!----------------------------------------
! main Body
!----------------------------------------
!
   CALL timing('green_inv',OPR='start')
   ! Init

      !
      dos(:) = ZERO
      g_aux(:,:) = CZERO
      i_work(:,:) =CZERO
      h_work(:,:) =CZERO
      !
      DO i_diag = 1, dimC
         !
         g_aux(i_diag,i_diag)= CONE
         i_work(i_diag,i_diag) = ene
         !
      ENDDO

      h_work(:,:) = i_work(:,:) - H_matrix_in(:,:) -  sgm_L(:,:) -  sgm_R(:,:)

      CALL mat_sv(dimC, dimC, h_work, g_aux(:,:))
      !

      DO i_diag = 1, dimC
         !
         dos(i_diag) =  - AIMAG( g_aux(i_diag,i_diag)) / PI
         !
      ENDDO

     g_CC(:,:) = g_aux(:,:)

   ! Final step
   !


   CALL timing('green_inv',OPR='stop')
   !

END SUBROUTINE green_CC_inv
