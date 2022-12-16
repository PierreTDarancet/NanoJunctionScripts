!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!               2007 Institut Neel CNRS/UJF Grenoble
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
!  datafile for hamiltonian for central region

   INTEGER :: n_orb=0
! dimension of the initial state for the recursion
!
   REAL(dbl) :: cell_size=ZERO
! useful to shift the cell

   REAL(dbl) :: limit_0=ZERO
! lower limit
   REAL(dbl) :: limit_1=ZERO
! limit for first neightbor
   REAL(dbl) :: limit_2=ZERO
! limit for second neightbor
   REAL(dbl) :: limit_3=ZERO
! limit for third neightbor

   LOGICAL :: use_second =.FALSE.
! allows to deal with a larger coupling (second neightbor)
!
   LOGICAL :: use_third =.FALSE.
! allows to deal with a larger coupling
!
   LOGICAL :: print_distance =.FALSE.
!
!
   LOGICAL :: print_orbital =.FALSE.
!
!
   LOGICAL :: print_hamiltonian =.TRUE.
!
!



   NAMELIST / CONTROL /    title, prefix, postfix, work_dir, datafile_H, &
                           cell_size, n_orb, limit_0, limit_1, limit_2, limit_3, use_second, use_third, &
                           print_orbital, print_distance,  print_hamiltonian


   PUBLIC :: title, prefix, postfix, work_dir
   PUBLIC ::  datafile_H
   PUBLIC :: cell_size
   PUBLIC :: n_orb
   PUBLIC :: limit_0
   PUBLIC :: limit_1
   PUBLIC :: limit_2
   PUBLIC :: limit_3
   PUBLIC :: use_second
   PUBLIC :: use_third
   PUBLIC :: print_orbital
   PUBLIC :: print_distance
   PUBLIC :: print_hamiltonian


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
              CALL errore(subname,'datafile_H unspecified',1)
!    
!          IF ( LEN_TRIM(datafile_L) == 0 ) &
!                CALL errore(subname,'datafile_L unspecified',2)
!    
!          IF ( LEN_TRIM(datafile_R) == 0 ) &
!                CALL errore(subname,'datafile_R unspecified',3)

        IF ( cell_size < 0.0 ) &
              CALL errore(subname,'bad cell_size',4)

        IF ( limit_0 == 0.0 ) &
              CALL errore(subname,'limit 0 unspecified ',5)

        IF ( limit_1 == 0.0 ) &
              CALL errore(subname,'limit 1 unspecified ',6)

        IF ( (use_second) .AND. (limit_2 == 0.0 ) ) &
              CALL errore(subname,'limit 2 unspecified',7)

        IF ( (use_third) .AND. (limit_3 == 0.0 ) ) &
              CALL errore(subname,'limit 3 unspecified',8)

        IF (  n_orb == 0 ) &
              CALL errore(subname,'n_orb unspecified',9)

        IF ( (.NOT. print_distance) .AND. (.NOT. print_orbital) .AND. (.NOT. print_hamiltonian)) &
              CALL errore(subname,'no output files',10)

   END SUBROUTINE read_namelist_control


END MODULE input_parameters_module

