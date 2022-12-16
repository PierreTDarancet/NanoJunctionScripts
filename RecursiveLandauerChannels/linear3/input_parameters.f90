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
   USE constants,     ONLY : ZERO, EPS_m2, EPS_m5, EPS_m10 
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
 !  CHARACTER(nstrx) :: datafile_H = ' '
!  datafile for hamiltonian for central region
   CHARACTER(nstrx) :: calculation_mode = 'change_V'
!  method to delete the phase factor in the states 
!   no_condition | change_V | diagonal A | ... | ...
!  no condition should be a debug option
!  change_V enforces the calculated recursion states to have positive real part
!  diagonal A calculates recursion states which diagonalize the A matrix
   CHARACTER(nstrx) :: variation_mode = 'diagonal'
!  method to determine the convergency diagonal | full
! 
   INTEGER :: max_iter=200
! max number of iteration
   INTEGER :: min_iter=0
! min number of iteration
   INTEGER :: dim_subspace=1
! dimension of the initial state for the recursion
   INTEGER :: dim_recursion=0
! dimension of the hamiltonian used for the recursion
   INTEGER :: n_orb=0
! dimension of the total hamiltonian, given in inputs
   INTEGER :: nb_max_first=1
! max number of first neighbor
   INTEGER :: nb_max_second=1
!
   INTEGER :: nb_max_third=1
!

!
   REAL(dbl) :: x_limit
! begin of the recursion chain
   REAL(dbl)  :: numerical_cut_off = EPS_m10
!
   REAL(dbl) :: limit_0=ZERO
! lower limit
   REAL(dbl) :: limit_1=ZERO
! limit for first neightbor
   REAL(dbl) :: limit_2=ZERO
! limit for second neightbor
   REAL(dbl) :: limit_3=ZERO
! limit for third neightbor
   REAL(dbl) :: conv_criteria=1.0
! convergency criteria

!
   LOGICAL :: use_second =.FALSE.
! allows to deal with a larger coupling (second neightbor)
   LOGICAL :: use_third =.FALSE.
! allows to deal with a larger coupling
   LOGICAL :: perform_recursion_x_lesser =.FALSE.
! take hamiltonian data with x >/< x0 to generate the recursion hamiltonian 

! Printed output files
   LOGICAL :: print_hamiltonian=.FALSE.
!
   LOGICAL :: print_center=.FALSE.
!
   LOGICAL :: print_state=.FALSE.
!
   LOGICAL :: print_variation=.TRUE.
!
   LOGICAL :: print_A_sum=.FALSE.
!
   LOGICAL :: print_B_sum=.FALSE.
!
   LOGICAL :: print_matrix=.TRUE.
!
   LOGICAL :: print_overlap=.FALSE.
!
   LOGICAL :: debug_mode=.FALSE.
!
   LOGICAL :: print_distance =.FALSE.
!
   LOGICAL :: print_orbital =.FALSE.
!
   LOGICAL :: perform_recursion_x_greater =.FALSE.
!
! Debug mode. no ending before max_iter error report

   NAMELIST / CONTROL /    title, prefix, postfix, work_dir, &
                           max_iter, conv_criteria, dim_subspace, &
                           dim_recursion, &
                           n_orb, x_limit, perform_recursion_x_greater, &
                           perform_recursion_x_lesser, &
                           variation_mode, calculation_mode, &
                           print_hamiltonian, &
                           print_center, print_state, &
                           print_variation, print_A_sum, &
                           print_B_sum, print_matrix, &
                           print_overlap, debug_mode, &
                           numerical_cut_off, min_iter, &
                           limit_0, limit_1, limit_2, &
                           limit_3, use_second, use_third, &
                           print_orbital, print_distance, &
                           nb_max_first, nb_max_second, nb_max_third


   PUBLIC :: title, prefix, postfix, work_dir
!   PUBLIC ::  datafile_H
   PUBLIC :: max_iter
   PUBLIC :: min_iter
   PUBLIC :: conv_criteria
   PUBLIC :: dim_subspace
   PUBLIC :: dim_recursion
   PUBLIC :: n_orb
   PUBLIC :: x_limit
   PUBLIC :: perform_recursion_x_greater
   PUBLIC :: perform_recursion_x_lesser
   PUBLIC :: variation_mode
   PUBLIC :: calculation_mode
   PUBLIC :: limit_0
   PUBLIC :: limit_1
   PUBLIC :: limit_2
   PUBLIC :: limit_3
   PUBLIC :: use_second
   PUBLIC :: use_third
   PUBLIC :: print_orbital
   PUBLIC :: print_distance
   PUBLIC :: print_hamiltonian
   PUBLIC :: print_center
   PUBLIC :: print_state
   PUBLIC :: print_variation
   PUBLIC :: print_A_sum
   PUBLIC :: print_B_sum
   PUBLIC :: print_matrix
   PUBLIC :: print_overlap
   PUBLIC :: debug_mode
   PUBLIC :: numerical_cut_off
   PUBLIC :: nb_max_first
   PUBLIC :: nb_max_second
   PUBLIC :: nb_max_third

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
!                IF ( LEN_TRIM(datafile_H) == 0 ) &
!                      CALL errore(subname,'datafile_H unspecified',2)
        IF ( n_orb == 0 ) &
              CALL errore(subname,'incorrect or unspecified n_orb',2)
        IF ( dim_recursion == 0 ) &
              CALL errore(subname,'incorrect or unspecified dim_recursion',3)
        IF ( dim_subspace == 0 ) &
              CALL errore(subname,'incorrect or unspecified dim_subspace',4)
!         IF ( dim_recursion >  dimwan ) &
!               CALL errore(subname,'dim_recursion > dimwan',5)
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

        IF ( limit_0 == 0.0 ) &
              CALL errore(subname,'limit 0 unspecified ',15)

        IF ( limit_1 == 0.0 ) &
              CALL errore(subname,'limit 1 unspecified ',16)

        IF ( (use_second) .AND. (limit_2 == 0.0 ) ) &
              CALL errore(subname,'limit 2 unspecified',17)

        IF ( (use_third) .AND. (limit_3 == 0.0 ) ) &
              CALL errore(subname,'limit 3 unspecified',18)


        IF (LEN_TRIM(variation_mode) == 0 ) &
              CALL errore(subname,'variation_mode unspecified',22)
        IF ((TRIM(variation_mode) /= 'full' )  .AND. (TRIM(variation_mode) /= 'diagonal' ) ) &
              CALL errore(subname,'incorrect variation_mode',23)

        IF (LEN_TRIM(calculation_mode) == 0 ) &
              CALL errore(subname,'calculation_mode unspecified',24)
        IF ((TRIM(calculation_mode) /= 'diagonal_A' )  .AND. (TRIM(calculation_mode) /= 'no_condition' ) &
                   .AND. (TRIM(calculation_mode) /= 'change_V' ) .AND. (TRIM(calculation_mode) /= 'diag_A_cond_V' ) &
                   .AND. (TRIM(calculation_mode) /= 'diag_A_sum_B' )  .AND. (TRIM(calculation_mode) /= 'diag_A_sgn_B')  &
                   .AND. (TRIM(calculation_mode) /= 'diag_B_B_dagger' )  ) &
              CALL errore(subname,'incorrect calculation_mode',25)


   END SUBROUTINE read_namelist_control


END MODULE input_parameters_module





 


