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
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI
   PRIVATE 
   SAVE
!
   IMPLICIT NONE

   ! Public
   ! Private

   ! Public variables

   ! Public routines:
   PUBLIC                 :: lead_green_function
   PUBLIC                 :: 
   PUBLIC                 :: 
			 
   CONTAINS 
!***********************************************
   SUBROUTINE lead_green_function(  gL(:,:), h00_L, h01_L, ene, dimL)
   !***********************************************
  IMPLICIT NONE
    CHARACTER(15)      :: subname="photongrid_init"
       INTEGER :: ierr, ie

        IF( .NOT. alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Not allocated"
            STOP
         ENDIF 
       !
    END SUBROUTINE


  END  MODULE T_green_operation_module

 
			 CALL calcul_gamma( gamma_L(ie,:,:), sigma_L(:,:) ,  gL(:,:), h_LC(:,:), h_CL(:,:), dimL, dimC )
 






   CALL gamma_allocate()



 gamma_L , gamma_R, sigma_L, sigma_R, lead_green_function, calcul_gamma





         aux00_L(:,:)  = h00_L(:,:) - (biasgrid(ibias)/2.00) !-ene * s00_L(:,:)
          aux01_L(:,:)  = h01_L(:,:) 
          !
          aux00_R(:,:)  = h00_R(:,:) + (biasgrid(ibias)/2.00) ! -ene * s00_R(:,:)
          aux01_R(:,:)  = h01_R(:,:) 
          !
          aux00_C(:,:)  = h00_C(:,:) -ene * s00_C(:,:) 
          aux_LC(:,:) = h_LC(:,:)
          aux_CR(:,:) = h_CR(:,:)
          !
          aux_CL(:,:) = CONJG( TRANSPOSE( h_LC(:,:) ))
          aux_RC(:,:) = CONJG( TRANSPOSE( h_CR(:,:) ))




         CALL mat_mul(work, aux_CR, 'N', gR, 'N', Nb_states, 1,1)
          CALL mat_mul(sigma_aux_R, work, 'N', aux_RC, 'N', Nb_states,Nb_states, 1)
           !
          CALL mat_mul(work, aux_CL, 'N', gL, 'N', Nb_states, 1, 1)
          CALL mat_mul(sigma_aux_L, work, 'N', aux_LC, 'N',Nb_states,Nb_states, 1) 
          !
          ! gamma_L and gamma_R
          !
          gamma_aux_L(:,:) = CI * (  sigma_aux_L(:,:) - CONJG( TRANSPOSE(sigma_aux_L(:,:)) )   )
          gamma_aux_R(:,:) = CI * (  sigma_aux_R(:,:) - CONJG( TRANSPOSE(sigma_aux_R(:,:)) )   )


    CALL gamma_deallocate() 
