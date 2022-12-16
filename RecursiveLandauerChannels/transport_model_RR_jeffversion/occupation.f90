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
   MODULE occupation_module
   !***********************************************
  USE kinds, ONLY : dbl
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI
  USE T_egrid_module, ONLY : nemax => ne

   IMPLICIT NONE
   PRIVATE 
   SAVE
!


   ! Public
   REAL(dbl), ALLOCATABLE ::  f_L(:) ! Occupation function fL (\omega - bias/2) where f is a fermi dirac distribution
   REAL(dbl), ALLOCATABLE ::  f_R(:) ! Occupation function fR (\omega + bias/2) where f is a fermi dirac distribution
  
   ! Private
   LOGICAL :: alloc = .FALSE.
   !REAL(dbl), PARAMETER :: kBT= 0.0001*ONE

   ! Public variables
   PUBLIC                 :: f_L
   PUBLIC                 :: f_R
 
   ! Public routines:
   PUBLIC                 :: occupation_allocate
   PUBLIC                 :: occupation_deallocate
   PUBLIC                 :: occupation_function_init
   PUBLIC                 :: occupation_function_onoff_init

   CONTAINS 
!***********************************************
   SUBROUTINE occupation_allocate
   !***********************************************
   IMPLICIT NONE
       CHARACTER(19)      :: subname="occupation_allocate"

       INTEGER  :: ierr
       !
        IF( alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Already allocated"
            STOP
         ENDIF 

        IF( nemax <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in nemax definition", " nemax ", nemax
            STOP
         ENDIF 

       ALLOCATE( f_L(nemax), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating fL"
            STOP
       ENDIF 
       ALLOCATE( f_R(nemax), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating fR"
            STOP
       ENDIF 
   

       f_R(:) = ZERO
       f_L(:) = ZERO


       alloc = .TRUE.


  END SUBROUTINE occupation_allocate
!***********************************************
   SUBROUTINE occupation_deallocate 
   !***********************************************
  IMPLICIT NONE
       CHARACTER(21)      :: subname="occupation_deallocate"
       INTEGER :: ierr


       IF ( ALLOCATED(f_R) ) THEN
            DEALLOCATE(f_R, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating f_R"
               STOP
           ENDIF 
       ENDIF

       IF ( ALLOCATED(f_L) ) THEN
            DEALLOCATE(f_L, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating f_L"
               STOP
           ENDIF 
       ENDIF
  

       alloc = .FALSE.

   END SUBROUTINE occupation_deallocate

!***********************************************
   SUBROUTINE occupation_function_init(f_L_aux, f_R_aux, ene, bias, shift, kBT)
   !***********************************************
  IMPLICIT NONE
   REAL(dbl), INTENT(out) ::  f_L_aux ! Occupation function fL (\omega - bias/2) where f is a fermi dirac distribution
   REAL(dbl), INTENT(out) ::  f_R_aux ! Occupation function fR (\omega + bias/2) where f is a fermi dirac distribution
   REAL(dbl), INTENT(in)  ::  ene
   REAL(dbl), INTENT(in)  ::  bias
   REAL(dbl), INTENT(in)  ::   shift
   REAL(dbl), INTENT(in)  ::  kBT


    CHARACTER(24)      :: subname="occupation_function_init"
    INTEGER :: ierr

        IF( .NOT. alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Not allocated"
            STOP
         ENDIF 
        !
      

        !(exp{(omega - \mu)/k_b T} + 1 )^{-1}
        ! Only Fermi Dirac distributions for the moment
	f_L_aux = ONE / ( EXP( ( ene - shift - ( bias/(2.0*ONE) ) )/ (kBT) )  + ONE) ! Left Electrode at the potential E_Fermi (+0eV) + shift + bias/2
	f_R_aux = ONE / ( EXP( ( ene - shift + ( bias/(2.0*ONE) ) )/ (kBT) )  + ONE) ! Right Electrode at the potential E_Fermi (+0eV) + shift - bias/2

   END SUBROUTINE occupation_function_init
!***********************************************
   SUBROUTINE occupation_function_onoff_init(f_L_aux, f_R_aux, ene, bias, shift, kBT)
   !***********************************************
  IMPLICIT NONE
   REAL(dbl), INTENT(out) ::  f_L_aux ! Occupation function fL (\omega - bias/2) where f is a fermi dirac distribution
   REAL(dbl), INTENT(out) ::  f_R_aux ! Occupation function fR (\omega + bias/2) where f is a fermi dirac distribution
   REAL(dbl), INTENT(in)  ::  ene
   REAL(dbl), INTENT(in)  ::  bias
   REAL(dbl), INTENT(in)  ::   shift
   REAL(dbl), INTENT(in)  ::  kBT


    CHARACTER(24)      :: subname="occupation_function_init"
    INTEGER :: ierr

        IF( .NOT. alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Not allocated"
            STOP
         ENDIF 
        !
      

        !(exp{(omega - \mu)/k_b T} + 1 )^{-1}
        ! Only Fermi Dirac distributions for the moment
        IF ( ene >  shift + (bias/(2.0*ONE)) ) THEN
             f_L_aux = ZERO
        ELSE  
             f_L_aux = ONE
        ENDIF 

       IF ( ene >  shift - (bias/(2.0*ONE)) ) THEN
             f_R_aux = ZERO
        ELSE  
             f_R_aux = ONE
        ENDIF 



   END SUBROUTINE occupation_function_onoff_init


  END  MODULE occupation_module

 
