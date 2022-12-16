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
   MODULE convert4_input_parameters_module
!********************************************
   USE kinds,         ONLY : dbl
   USE constants,     ONLY : ZERO, EPS_m2, EPS_m5 
   USE parameters,    ONLY : nstrx
   !USE parser_module, ONLY : change_case
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
   CHARACTER(nstrx) :: title = "Convert to Wannier Transport Calculation"
       ! the title of the calculation

   CHARACTER(nstrx) :: prefix = "WanT"
       ! specifies the prefix for the names of all the output and input (data) files
       ! INPUT  files:   "prefix".aaa
       ! OUTPUT files:   "prefix""postfix".bbb

   CHARACTER(nstrx) :: postfix = " "
       ! specifies the second part of the names of the output files
   
   CHARACTER(nstrx) :: work_dir = "./"
       ! the directory in which produced data files are written 

!
   CHARACTER(nstrx) :: out_datafile_C = ' '
!  datafile for hamiltonian for central region
   CHARACTER(nstrx) :: out_datafile_L = ' '
!  datafile for hamiltonian left lead
   CHARACTER(nstrx) :: out_datafile_R = ' '
!  datafile for hamiltonian right lead
!
   CHARACTER(nstrx) :: in_datafile_L = ' '
!  datafile for hamiltonian left lead
   CHARACTER(nstrx) :: in_datafile_R = ' '
!  datafile for hamiltonian right lead

!
  INTEGER :: out_iter_C=0
!
!
  INTEGER :: out_iter_R=0
!
!
  INTEGER :: out_iter_L=0
!
!
   INTEGER :: n_iter_L=0
! max number of iteration
!
   INTEGER :: n_iter_R=0
! max number of iteration

   INTEGER :: dim_subspace=1
! dimension of the initial state for the recursion
!


   NAMELIST / CONTROL /    title, prefix, postfix, work_dir, out_datafile_C, &
                           out_datafile_L, out_datafile_R, in_datafile_L, in_datafile_R, &
                           n_iter_L, n_iter_R, dim_subspace, out_iter_C, out_iter_L, out_iter_R

   PUBLIC :: title, prefix, postfix, work_dir
   PUBLIC ::  out_datafile_L
   PUBLIC ::  out_datafile_R
   PUBLIC ::  out_datafile_C
   PUBLIC ::  in_datafile_R
   PUBLIC ::  in_datafile_L
   PUBLIC ::  dim_subspace
   PUBLIC :: n_iter_L
   PUBLIC :: n_iter_R
   PUBLIC :: out_iter_C
   PUBLIC :: out_iter_L
   PUBLIC :: out_iter_R



   PUBLIC :: CONTROL

   PUBLIC ::  convert4_read_namelist_control


CONTAINS

!**********************************************************
   SUBROUTINE convert4_read_namelist_control(unit)
   !**********************************************************
   !
   ! reads CONTROL namelist
   !
   IMPLICIT NONE
      INTEGER, INTENT(in)   :: unit

      CHARACTER(30) :: subname='convert4_read_namelist_control'
      INTEGER :: ios

      READ(unit, CONTROL, IOSTAT=ios )
         IF (ios/=0) CALL errore(subname,'reading CONTROL namelist',ABS(ios))

      !
      ! ... checking parameters
      !

        ! Check the completeness of the input's data with respect to problem's dimensionality
        IF ( LEN_TRIM(out_datafile_C) == 0 ) &
              CALL errore(subname,'out_datafile_C unspecified',1)
        IF ( LEN_TRIM(out_datafile_R) == 0 ) &
              CALL errore(subname,'out_datafile_R unspecified',2)
        IF ( LEN_TRIM(out_datafile_L) == 0 ) &
              CALL errore(subname,'out_datafile_L unspecified',3)
        IF ( LEN_TRIM(in_datafile_R) == 0 ) &
              CALL errore(subname,'in_datafile_R unspecified',4)
        IF ( LEN_TRIM(in_datafile_L) == 0 ) &
              CALL errore(subname,'in_datafile_L unspecified',5)

        IF ( dim_subspace == 0 ) &
              CALL errore(subname,'incorrect or unspecified dim_subspace',6)

        IF ( out_iter_L == 0 ) &
              CALL errore(subname,'bad number of iteration in L part',7)
        IF ( out_iter_R == 0 ) &
              CALL errore(subname,'bad number of iteration in R part',7)
        IF ( out_iter_C == 0 ) &
              CALL errore(subname,'bad number of iteration in C part',7)

        IF ( n_iter_L == 0 ) &
              CALL errore(subname,'bad number of iteration in L part',7)
        IF ( n_iter_R == 0 ) &
              CALL errore(subname,'bad number of iteration in R part',8)

   END SUBROUTINE convert4_read_namelist_control


END MODULE convert4_input_parameters_module

