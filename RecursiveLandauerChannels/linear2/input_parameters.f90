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

 
   CHARACTER(nstrx) :: datafile_H = ' '
!  datafile for hamiltonian for central region
   CHARACTER(nstrx) :: datafile_L = ' '
!  datafile for hamiltonian left lead
   CHARACTER(nstrx) :: datafile_R = ' '
!  datafile for hamiltonian right lead

   CHARACTER(nstrx) :: calculation_mode = 'change_V'
!  method to delete the phase factor in the states 
!   no_condition | change_V | diagonal A
!  no condition should be a debug option
!  change_V enforces the calculated recursion states to have positive real part
!  diagonal A calculates recursion states which diagonalize the A matrix

   CHARACTER(nstrx) :: variation_mode = 'diagonal'
!  method to determine the convergency diagonal | full
! 

!
   INTEGER :: max_iter=200
! max number of iteration
   INTEGER :: min_iter=0
! min number of iteration

!
   REAL(dbl) :: conv_criteria=1.0
! convergency criteria
!
   INTEGER :: dim_subspace=1
! dimension of the initial state for the recursion
!
   INTEGER :: dim_recursion=0
! dimension of the hamiltonian used for the recursion
!
   INTEGER :: dimwan=0
! dimension of the total hamiltonian, given in inputs
!
   REAL(dbl) :: x_limit
! begin of the recursion chain
!
!
   REAL(dbl) :: cell_size_C=0
! useful to shift the wannier center when use_extension_L/R = .TRUE.
   REAL(dbl) :: cell_size_R=0
! useful to shift the wannier center when use_extension_L/R = .TRUE. and nb_replica /=0
   REAL(dbl) :: cell_size_L=0
! useful to shift the wannier center when use_extension_L/R = .TRUE. and nb_replica /=0
!
   REAL(dbl)  :: numerical_cut_off = EPS_m10
!
!

   INTEGER :: dim_L=0
! dimension of the left lead hamiltonian (1 replica)
!
   INTEGER :: dim_R=0
! dimension of the right lead hamiltonian (1 replica)
!
   INTEGER :: nb_replica_L=0
! left  lead replica's number used in the recursion calculation
!
   INTEGER :: nb_replica_R=0
! right lead replica's number used in the recursion calculation
!
   LOGICAL :: use_extension_R =.FALSE.
! allows to deal with a bigger hamiltonian in the right part
!
   LOGICAL :: use_extension_L =.FALSE.
! allows to deal with a bigger hamiltonian in the left part
!
   LOGICAL :: perform_recursion_x_greater =.FALSE.
   LOGICAL :: perform_recursion_x_lesser =.FALSE.
! take hamiltonian data with x >/< x0 to generate the recursion hamiltonian 
!
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
! Debug mode. no ending before max_iter error report

   NAMELIST / CONTROL /    title, prefix, postfix, work_dir, datafile_H, max_iter, conv_criteria, dim_subspace, dim_recursion, &
                           dimwan, x_limit, perform_recursion_x_greater, perform_recursion_x_lesser, datafile_L, datafile_R, &
                           nb_replica_L, nb_replica_R, dim_L, dim_R, use_extension_L, use_extension_R, cell_size_C, &
                           cell_size_R, cell_size_L, variation_mode, calculation_mode, print_hamiltonian, &
                           print_center, print_state, print_variation, print_A_sum, print_B_sum, print_matrix, &
                           print_overlap, debug_mode, numerical_cut_off, min_iter

   PUBLIC :: title, prefix, postfix, work_dir
   PUBLIC ::  datafile_H
   PUBLIC ::  max_iter
   PUBLIC ::  min_iter
   PUBLIC :: conv_criteria
   PUBLIC ::  dim_subspace
   PUBLIC ::  dim_recursion
   PUBLIC ::  dimwan
   PUBLIC ::  x_limit
   PUBLIC :: perform_recursion_x_greater
   PUBLIC :: perform_recursion_x_lesser
   PUBLIC :: datafile_L
   PUBLIC :: datafile_R
   PUBLIC :: nb_replica_L
   PUBLIC :: nb_replica_R
   PUBLIC :: dim_L
   PUBLIC :: dim_R
   PUBLIC :: use_extension_L
   PUBLIC :: use_extension_R
   PUBLIC :: cell_size_C
   PUBLIC :: cell_size_R
   PUBLIC :: cell_size_L
   PUBLIC :: variation_mode
   PUBLIC :: calculation_mode
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
        IF ( (perform_recursion_x_greater).AND.(use_extension_L) ) &
              CALL errore(subname,'both perform_recursion_x_greater and use_extension_L are true',14)
        IF ( (perform_recursion_x_lesser).AND.(use_extension_R) ) &
              CALL errore(subname,'both perform_recursion_x_lesser and use_extension_R are true',14)

        IF ( use_extension_R ) THEN
            !
            IF ( ( dim_R == 0 ) .AND. perform_recursion_x_greater) &
                  CALL errore(subname,'incorrect or unspecified dim_R ',15)
            IF ( perform_recursion_x_greater .AND. ( LEN_TRIM(datafile_R) == 0 ) ) &
                  CALL errore(subname,'incorrect or unspecified datafile_R ',16)
          !  IF ( perform_recursion_x_greater .AND. ( nb_replica_R == 0 ) ) &
          !        CALL errore(subname,'incorrect or unspecified nb_replica_R ',17)
           IF ( perform_recursion_x_greater .AND. ( cell_size_C == 0 ) ) &
                  CALL errore(subname,'incorrect or unspecified cell_size_C ',17)
           IF ( ( nb_replica_R /= 0 ) .AND. ( cell_size_R == 0 ) ) &
                  CALL errore(subname,'incorrect or unspecified cell_size_R ',17)

        ENDIF

        IF ( use_extension_L ) THEN
            !
            IF ( ( dim_L == 0 ) .AND. perform_recursion_x_lesser) &
                  CALL errore(subname,'incorrect or unspecified dim_L ',18)
            IF ( perform_recursion_x_lesser .AND. ( LEN_TRIM(datafile_L) == 0 ) ) &
                  CALL errore(subname,'incorrect or unspecified datafile_L ',19)
          !  IF ( perform_recursion_x_lesser .AND. ( nb_replica_L == 0 ) ) &
          !        CALL errore(subname,'incorrect or unspecified nb_replica_L ',20)
           IF ( perform_recursion_x_lesser .AND. ( cell_size_C == 0 ) ) &
                  CALL errore(subname,'incorrect or unspecified cell_size_C ',20)
           IF ( ( nb_replica_L /= 0 ) .AND. ( cell_size_L == 0 ) ) &
                  CALL errore(subname,'incorrect or unspecified cell_size_L ',20)

        ENDIF

        IF (( use_extension_L ) .AND. ( use_extension_R )) &
           CALL errore(subname,'both use_extension_L and use_extension_R are true',21)

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

