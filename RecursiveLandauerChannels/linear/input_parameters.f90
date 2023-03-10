!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!********************************************
   MODULE input_parameters_module
!********************************************
   USE kinds,         ONLY : dbl
   USE constants,     ONLY : ZERO, EPS_m2, EPS_m5 
   USE parameters,    ONLY : nstrx
   USE parser_module, ONLY : change_case
   IMPLICIT NONE
   PRIVATE
   SAVE
!
! This module contains the definitions of the parameters in the
! input file and of thier default values (when any). 
! These data are then exported to the input module where the
! main routine controls the IO and after that exports to the
! final modules where internal data are stored.
!
! Here are also the routine reading and checking the NAMELIST
!
! routines in this module:
! SUBROUTINE  read_namelist_Control(unit)
!

!
! ... declarations

!
!======================================== 
! INPUT_CONDUCTOR Namelist parameters
!======================================== 
!
   CHARACTER(nstrx) :: title = "Wannier Transport Calculation"
       ! the title of the calculation

   CHARACTER(nstrx) :: prefix = "WanT"
       ! specifies the prefix for the names of all the output and input (data) files
       ! INPUT  files:   "prefix".aaa
       ! OUTPUT files:   "prefix""postfix".bbb

   CHARACTER(nstrx) :: postfix = " "
       ! specifies the second part of the names of the output files
   
   CHARACTER(nstrx) :: work_dir = "./"
       ! the directory in which produced data files are written 

 
   CHARACTER(nstrx) :: datafile_H = ' '
!  datafile for hamiltonian

!
   INTEGER :: max_iter=200
! max number of iteration

!
   REAL(dbl) :: conv_criteria=1.0
! convergency criteria
!
   INTEGER :: dim_subspace=1
! dimension of the initial state for the recursion
!
   INTEGER :: dim_recursion
! dimension of the hamiltonian used for the recursion
!
   INTEGER :: dimwan
! dimension of the total hamiltonian, given in inputs
!
   REAL(dbl) :: x_limit
! begin of the recursion chain
!
   LOGICAL :: perform_recursion_x_greater =.FALSE.
   LOGICAL :: perform_recursion_x_lesser =.FALSE.
! take hamiltonian data with x >/< x0 to generate the recursion hamiltonian 
!




   NAMELIST / CONTROL /    title, prefix, postfix, work_dir, datafile_H, max_iter, conv_criteria, dim_subspace, dim_recursion, &
                           dimwan, x_limit, perform_recursion_x_greater, perform_recursion_x_lesser

   PUBLIC :: title, prefix, postfix, work_dir
   PUBLIC ::  datafile_H
   PUBLIC ::  max_iter
   PUBLIC :: conv_criteria
   PUBLIC ::  dim_subspace
   PUBLIC ::  dim_recursion
   PUBLIC ::  dimwan
   PUBLIC ::  x_limit
   PUBLIC :: perform_recursion_x_greater
   PUBLIC :: perform_recursion_x_lesser

   PUBLIC :: CONTROL

   PUBLIC ::  read_namelist_control


CONTAINS

!**********************************************************
   SUBROUTINE read_namelist_control(unit)
   !**********************************************************
   !
   ! reads CONTROL namelist
   !
   IMPLICIT NONE
      INTEGER, INTENT(in)   :: unit

      CHARACTER(21) :: subname='read_namelist_control'
      INTEGER :: ios

      READ(unit, CONTROL, IOSTAT=ios )
         IF (ios/=0) CALL errore(subname,'reading CONTROL namelist',ABS(ios))

      !
      ! ... checking parameters
      !

        ! Check the completeness of the input's data with respect to problem's dimensionality
        IF ( LEN_TRIM(datafile_H) == 0 ) &
              CALL errore(subname,'datafile_H unspecified',2)
        IF ( dimwan == 0 ) &
              CALL errore(subname,'incorrect or unspecified dimwan',2)
        IF ( dim_recursion == 0 ) &
              CALL errore(subname,'incorrect or unspecified dim_recursion',3)
        IF ( dim_subspace == 0 ) &
              CALL errore(subname,'incorrect or unspecified dim_subspace',4)
        IF ( dim_recursion >  dimwan ) &
              CALL errore(subname,'dim_recursion > dimwan',5)
        IF ( dim_subspace > dim_recursion ) &
              CALL errore(subname,'dim_subspace > dim_recursion',6)
        IF ( conv_criteria <= 0.0 ) &
              CALL errore(subname,'bad criteria for convergency',7)
        IF ( max_iter == 0 ) &
              CALL errore(subname,'bad maximum of iteration',8)
        IF ( x_limit < 0.0 ) &
              CALL errore(subname,'bad x for recursion hamiltonian',9)
        IF ( (.NOT. perform_recursion_x_greater).AND.(.NOT. perform_recursion_x_lesser) ) &
              CALL errore(subname,'perform_recursion_x_greater and recursion_x_lesser are false',11)
        IF ( (perform_recursion_x_greater).AND.(perform_recursion_x_lesser) ) &
              CALL errore(subname,'both perform_recursion_x_greater and recursion_x_lesser are true',12)
        IF ( max_iter*dim_subspace > dim_recursion  ) &
              CALL errore(subname,'too big number of iterations',13)


   END SUBROUTINE read_namelist_control


END MODULE input_parameters_module

