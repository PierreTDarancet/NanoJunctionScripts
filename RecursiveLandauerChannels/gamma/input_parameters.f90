! 
! Copyright (C) 2005 WanT Group
! 
! This file is distributed under the terms of the 
! GNU General Public License. See the file `License' 
! in the root directory of the present distribution, 
! or http://www.gnu.org/copyleft/gpl.txt . 
! 
!********************************************
   MODULE T_input_parameters_module
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
! The T character in front of the module name stands there
! in order to distinguish from wannier modules.
!
! Here are also the routine reading and checking the NAMELIST
!
! routines in this module:
! SUBROUTINE  read_namelist_input_conductor(unit)
!

!
! ... declarations

!
!======================================== 
! INPUT_CONDUCTOR Namelist parameters
!======================================== 

! General information



! Input files



   CHARACTER(nstrx) :: in_datafile = ' '
!  datafile for hamiltonian lead


   CHARACTER(nstrx) :: method_sigma = ' '
!  datafile for hamiltonian lead

! 

!
   INTEGER :: in_max_iter=0
! max number of iteration
!
   INTEGER :: dim_subspace = 0
       ! WF number in the grenn region

! Leads SE part


  INTEGER :: max_iter_term = 0
       ! recursion iter  


   !REAL(dbl) :: conv_criterion = 0.000001
       ! lower bound for recursion
!   REAL(dbl) :: numerical_cut_off = 0.00000001
!

  INTEGER :: max_iter_renorm = 0
       !  recursion iter

   INTEGER :: nprint = 20 
       ! every nprint energy step write to stdout


   !LOGICAL :: cut_chain = .FALSE.

   LOGICAL :: print_gamma = .TRUE.


   LOGICAL :: print_gamma_tilde = .TRUE.

! Transmittance part


   INTEGER :: ne = 1000  
       ! the dimension of the energy grid

   REAL(dbl) :: emin = -10.0
       ! lower bound of the energy range

   REAL(dbl) :: emax =  10.0
       ! upper bound of the energy range

   REAL(dbl) :: delta_lead =  EPS_m5
       ! i\delta broadening of green functions

   REAL(dbl) :: a_analytique=0.00
       ! 

   REAL(dbl) :: b_analytique = 4.45
       !



   NAMELIST / INPUT_CONDUCTOR /  ne, emin, emax, nprint, delta_lead, &
                 in_datafile, in_max_iter, dim_subspace,  &
                 max_iter_term, max_iter_renorm, &
                 print_gamma_tilde, print_gamma, method_sigma, &
                 a_analytique, b_analytique
    !

   PUBLIC :: dim_subspace, max_iter_term, max_iter_renorm
   PUBLIC :: ne, emin, emax, nprint, delta_lead
   PUBLIC :: in_datafile, in_max_iter
   PUBLIC :: INPUT_CONDUCTOR
   PUBLIC :: a_analytique, b_analytique
!   PUBLIC :: dimR, dimL
   PUBLIC :: print_gamma
   PUBLIC :: print_gamma_tilde
   PUBLIC :: method_sigma
!
!
! ... subroutines
   PUBLIC :: read_namelist_input_conductor

CONTAINS

!**********************************************************
   SUBROUTINE read_namelist_input_conductor(unit)
   !**********************************************************
   !
   ! reads INPUT_CONDUCTOR namelist
   !
   IMPLICIT NONE
      INTEGER, INTENT(in)   :: unit

      CHARACTER(29) :: subname='read_namelist_input_conductor'
      LOGICAL :: allowed
      INTEGER :: i, ios, itest

      READ(unit, INPUT_CONDUCTOR, IOSTAT=ios )
         IF (ios/=0) CALL errore(subname,'reading INPUT_CONDUCTOR namelist',ABS(ios))

      !
      ! ... checking parameters
      !

      !IF ( dimC <= 0) CALL errore(subname,'Invalid dimC',1)

     ! IF ( LEN_TRIM(datafile_C) == 0 ) &
     !      CALL errore(subname,'datafile_C unspecified',1)

      IF ( emax <= emin ) CALL errore(subname,'Invalid EMIN EMAX',1)
      IF ( ne <= 1 ) CALL errore(subname,'Invalid NE',2)
      IF ( nprint <= 0) CALL errore(subname, ' nprint must be > 0 ', -nprint+1 )
      !IF ( delta < ZERO ) CALL errore(subname,'Invalid DELTA',3)
      !IF ( delta > EPS_m2 ) CALL errore(subname,'DELTA too large',4)


           !
      IF ( LEN_TRIM(in_datafile) == 0 ) &
            CALL errore(subname,'datafile unspecified',5)
      IF (max_iter_renorm == 0) &
             CALL errore(subname,'max_iter_renorm unspecified',11)

      IF (max_iter_term == 0) &
             CALL errore(subname,'max_iter_terminaison unspecified',12)

      IF (in_max_iter == 0) &
             CALL errore(subname,'in_max_iter unspecified',13)

      IF ( (max_iter_renorm  + max_iter_term) /=  in_max_iter ) &
             CALL errore(subname,'term + renorm /= in max_iter',14)


           !
      IF ( (LEN_TRIM(method_sigma) == 0) .OR. ((TRIM(method_sigma) /= 'WanT' ) .AND. (TRIM(method_sigma) /= 'analytique' ))) &
            CALL errore(subname,'method_sigma',15)

           !
   END SUBROUTINE read_namelist_input_conductor


END MODULE T_input_parameters_module

