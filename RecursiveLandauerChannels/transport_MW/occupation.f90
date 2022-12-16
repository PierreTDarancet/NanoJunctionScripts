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
  USE T_egrid_module, ONLY : ne_green, fermi_energy, electronic_temperature
   IMPLICIT NONE
   PRIVATE 
   SAVE
!


   ! Public
   REAL(dbl), ALLOCATABLE ::  f_L(:) ! Occupation function fL (\omega - bias/2) where f is a fermi dirac distribution
   REAL(dbl), ALLOCATABLE ::  f_R(:) ! Occupation function fR (\omega + bias/2) where f is a fermi dirac distribution
 
   ! Private
   LOGICAL :: alloc = .FALSE.


   ! Public variables
   PUBLIC                 :: f_L
   PUBLIC                 :: f_R

   ! Public routines:
   PUBLIC                 :: occupation_allocate
   PUBLIC                 :: occupation_deallocate
   PUBLIC                 :: occupation_function_init

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

        IF( ne_green <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in nemax definition", " ne_green", ne_green
            STOP
         ENDIF 


         ALLOCATE( f_L(ne_green), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating fL"
            STOP
       ENDIF 
       ALLOCATE( f_R(ne_green), STAT=ierr )
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
   SUBROUTINE occupation_function_init(egrid)
   !***********************************************
  IMPLICIT NONE
    REAL(dbl), INTENT(in) :: egrid(ne_green)
    CHARACTER(24)      :: subname="occupation_function_init"
    INTEGER :: ierr, ie

        IF( .NOT. alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Not allocated"
            STOP
         ENDIF 
        !
        DO ie=1, ne_green
            f_L(ie) = ONE / ( EXP(  (egrid(ie)-fermi_energy)/electronic_temperature )  + ONE )
        ENDDO

        f_R(:) = f_L(:)  
        !(exp{(omega - \mu)/k_b T} + 1 )^{-1}
        ! Only Fermi Dirac distributions for the moment

   END SUBROUTINE occupation_function_init


  END  MODULE occupation_module

 
