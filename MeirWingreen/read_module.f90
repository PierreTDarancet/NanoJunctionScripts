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
   MODULE read_module
   !***********************************************
  USE kinds, ONLY : dbl
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI

   IMPLICIT NONE
   PRIVATE 
   SAVE
   !
   ! Public
   !
   ! Public routines:

   PUBLIC                 :: read_mol_selfenergy
   PUBLIC                 :: read_imagecharge_selfenergy
   PUBLIC                 :: read_lead_selfenergy
   PUBLIC                 :: read_lead_gamma
   PUBLIC                 :: read_electron_photon


   CONTAINS 

!***********************************************
   SUBROUTINE read_mol_selfenergy(sigma_ee_mol_r,dimC, method)
   !***********************************************
  IMPLICIT NONE

    COMPLEX(dbl), INTENT(out) :: sigma_ee_mol_r(dimC,dimC)
    INTEGER, INTENT(in)       :: dimC
    CHARACTER(3), INTENT(in)         ::  method
    CHARACTER(15)         ::     fname_in="SIGMA_MOL_GW.in"
    CHARACTER(19)      :: subname="read_mol_selfenergy"
    INTEGER :: ierr, i_dim

        IF ((TRIM(method)=="EXC" ).OR.((TRIM(method)=="SIG")) THEN
		!
	 	PRINT*, "     ... OPENING File:", fname_in
		OPEN(100, FILE=TRIM(fname_in),FORM='formatted')
		!
		sigma_ee_mol_r(:,:) = CZERO
	  	PRINT*, "     ...... READING GW CORRECTIONS"
	  	PRINT*, "     ...... WARNING: Only diagonal corrections are accepted at this point"
		DO i_dim=1,dimC
		   !
		   READ(100,"(a11,I10)")  sigma_ee_mol_r(i_dim,i_dim) 
		   PRINT*, "          n=", i_dim, " GW correction= ",  sigma_ee_mol_r(i_dim,i_dim) 
		   !
		ENDDO
		!
	  	PRINT*, "     ... CLOSING File:", fname_in
		CLOSE(100)
        ELSE
	  	PRINT*, "     ...... WARNING: DFT Calculation, no Self energy read, values set to ZERO"
		sigma_ee_mol_r(:,:) = CZERO

        ENDIF

   END SUBROUTINE read_mol_selfenergy


!***********************************************
   SUBROUTINE read_imagecharge_selfenergy(sigma_ee_contact_r(:,:),dimC, method)
   !***********************************************
  IMPLICIT NONE
    COMPLEX(dbl), INTENT(out) :: sigma_ee_contact_r(dimC,dimC)
    INTEGER, INTENT(in)       :: dimC
    CHARACTER(3), INTENT(in)         ::  method
    CHARACTER(19)         ::     fname_in="SIGMA_CONTACT_GW.in"
    CHARACTER(27)      :: subname="read_imagecharge_selfenergy"
    INTEGER :: ierr, i_dim

        IF ((TRIM(method)=="EXC" ).OR.((TRIM(method)=="SIG")) THEN
		!
	 	PRINT*, "     ... OPENING File:", fname_in
		OPEN(100, FILE=TRIM(fname_in),FORM='formatted')
		!
		sigma_ee_contact_r(:,:)= CZERO
	  	PRINT*, "     ...... READING IMAGE CHARGE CORRECTIONS"
	  	PRINT*, "     ...... WARNING: Only diagonal corrections are accepted at this point"
		DO i_dim=1,dimC
		   !
		   READ(100,"(a11,I10)")  sigma_ee_contact_r(i_dim,i_dim)
		   PRINT*, "          n=", i_dim, " GW correction= ",  sigma_ee_contact_r(i_dim,i_dim)
		   !
		ENDDO
		!
	  	PRINT*, "     ... CLOSING File:", fname_in
		CLOSE(100)

        ELSE
	  	PRINT*, "     ...... WARNING: DFT Calculation, no Self energy read, values set to ZERO"
		sigma_ee_contact_r(:,:)= CZERO

        ENDIF

   END SUBROUTINE read_imagecharge_selfenergy


!***********************************************
   SUBROUTINE read_lead_selfenergy( sigma_L, sigma_R, ik, ne_green, dimC)
   !***********************************************
  IMPLICIT NONE

    COMPLEX(dbl), INTENT(out) :: sigma_L(ne_green,dimC,dimC)
    COMPLEX(dbl), INTENT(out) :: sigma_R(ne_green,dimC,dimC)
    INTEGER, INTENT(in)       :: dimC
    INTEGER, INTENT(in)       :: ne_green
    INTEGER, INTENT(in)       :: ik
    CHARACTER(100)         ::   fname_in 
    CHARACTER(20)      :: subname="read_lead_selfenergy"
    INTEGER :: ierr, i_dim, ie, irow, icol, i_test


        WRITE (fname_in, *) "SIGMA_L", ik, ".in"
        !
 	PRINT*, "     ... OPENING File:", TRIM(fname_in)
	OPEN(100, FILE=TRIM(fname_in),FORM='formatted')
	!
        sigma_L(:,:,:)= CZERO
  	PRINT*, "     ...... READING SIGMA_L"
  	PRINT*, "     ...... WARNING: The grid is assumed to be equal to the one used in ne_green"
        i_test =0
        DO ie=1,ne_green
        PRINT*, "Energy Point", ie
        DO i_col=1,dimC
           DO i_row=1,dimC
           !
           READ(100,"(a11,I10)")  sigma_L(ie,i_row,i_col)
           i_test = i_test + 1
           !
           ENDDO
        ENDDO
        ENDDO
        !

  	PRINT*, "     ... CLOSING File:", fname_in
        CLOSE(100)

      WRITE (fname_in, *) "SIGMA_R", ik, ".in"
        !
 	PRINT*, "     ... OPENING File:", TRIM(fname_in)
	OPEN(100, FILE=TRIM(fname_in),FORM='formatted')
	!
        sigma_R(:,:,:)= CZERO
  	PRINT*, "     ...... READING SIGMA_R"
  	PRINT*, "     ...... WARNING: The grid is assumed to be equal to the one used in ne_green"
        i_test =0
        DO ie=1,ne_green
        PRINT*, "Energy Point", ie
        DO i_col=1,dimC
           DO i_row=1,dimC
           !
           READ(100,"(a11,I10)")  sigma_R(ie,i_row,i_col)
           i_test = i_test + 1
           !
           ENDDO
        ENDDO
        ENDDO
        !

  	PRINT*, "     ... CLOSING File:", fname_in
        CLOSE(100)


   END SUBROUTINE read_lead_selfenergy

!***********************************************
   SUBROUTINE read_lead_gamma( gamma_L(:,:,:), gamma_R(:,:,:), ik,  ne_green,dimC)
   !***********************************************
  IMPLICIT NONE
  
    COMPLEX(dbl), INTENT(out) :: gamma_L(ne_green,dimC,dimC)
    COMPLEX(dbl), INTENT(out) :: gamma_R(ne_green,dimC,dimC)
    INTEGER, INTENT(in)       :: dimC
    INTEGER, INTENT(in)       :: ne_green
    INTEGER, INTENT(in)       :: ik
    CHARACTER(100)         ::   fname_in 
    CHARACTER(15)      :: subname="read_lead_gamma"
    INTEGER :: ierr, i_dim, ie, irow, icol, i_test


        WRITE (fname_in, *) "GAMMA_L", ik, ".in"
        !
 	PRINT*, "     ... OPENING File:", TRIM(fname_in)
	OPEN(100, FILE=TRIM(fname_in),FORM='formatted')
	!
        gamma_L(:,:,:)= CZERO
  	PRINT*, "     ...... READING GAMMA_L"
  	PRINT*, "     ...... WARNING: The grid is assumed to be equal to the one used in ne_green"
        i_test =0
        DO ie=1,ne_green
        PRINT*, "Energy Point", ie
        DO i_col=1,dimC
           DO i_row=1,dimC
           !
           READ(100,"(a11,I10)")  gamma_L(ie,i_row,i_col)
           i_test = i_test + 1
           !
           ENDDO
        ENDDO
        ENDDO
        !

  	PRINT*, "     ... CLOSING File:", fname_in
        CLOSE(100)

      WRITE (fname_in, *) "GAMMA_R", ik, ".in"
        !
 	PRINT*, "     ... OPENING File:", TRIM(fname_in)
	OPEN(100, FILE=TRIM(fname_in),FORM='formatted')
	!
        gamma_R(:,:,:)= CZERO
  	PRINT*, "     ...... READING GAMMA_R"
  	PRINT*, "     ...... WARNING: The grid is assumed to be equal to the one used in ne_green"
        i_test =0
        DO ie=1,ne_green
        PRINT*, "Energy Point", ie
        DO i_col=1,dimC
           DO i_row=1,dimC
           !
           READ(100,"(a11,I10)")  gamma_R(ie,i_row,i_col)
           i_test = i_test + 1
           !
           ENDDO
        ENDDO
        ENDDO
        !

  	PRINT*, "     ... CLOSING File:", fname_in
        CLOSE(100)



  END SUBROUTINE read_lead_gamma

!***********************************************
   SUBROUTINE read_electron_photon( dipole_matrix, dimC, ep_method)
   !***********************************************
  IMPLICIT NONE
    REAL(dbl), INTENT(out) :: dipole_matrix(dimC,dimC)
    INTEGER, INTENT(in)       :: dimC
    CHARACTER(2), INTENT(in)  :: ep_method
    CHARACTER(9)         ::     fname_in="DIPOLE.in"
    CHARACTER(20)        :: subname="read_electron_photon"
    INTEGER :: ierr, i_row, i_col



        dipole_matrix(:,:) = ZERO

        IF ( (TRIM(method)=="AM") .OR. (TRIM(method)=="GN") ) THEN ! Aeberhard Morf PRB 77 (2008) - SCBA - Galperin Nitzan PRL 2005 with dipole transition/
		!
		PRINT*, "     ... OPENING File:", fname_in
		OPEN(100, FILE=TRIM(fname_in),FORM='formatted')
		PRINT*, "     ...... READING DIPOLE"
		!
		DO i_col=1,dimC
		   DO i_row=1,dimC
		   !
		   READ(100,"(a11,I10)")  dipole_matrix(i_row,i_col)
		   !
		   ENDDO
		ENDDO
		!
		PRINT*, "     ... CLOSING File:", fname_in
		!
                CLOSE(100)

        ELSE IF (TRIM(method)=="G0") THEN !  Galperin Nitzan PRL 2005 without dipole transition -set to one
		PRINT*, "     ...... WARNING ---"
		PRINT*, "     ...... METHOD is G0 ---"
		PRINT*, "     ...... dipole transition is set to one off diagonal ---"
		DO i_col=1,dimC
		   DO i_row=1,dimC
		   !
                   IF ( ( irow <= INT(dimC/2.00) ) .AND. ( icol > INT(dimC/2.00)) ) THEN
             		   dipole_matrix(i_row,i_col) = ONE
                   ELSE IF  ( ( icol <= INT(dimC/2.00) ) .AND. ( irow > INT(dimC/2.00)) ) THEN
             		   dipole_matrix(i_row,i_col) = ONE
                   ENDIF
		   !
                   ENDDO
		ENDDO
        ELSE
		PRINT*, "     ...... WARNING ---"
		PRINT*, "     ...... No dipole transition are considered ---"
		PRINT*, "     ...... Electron Photon SE set to 0 ---"
        ENDIF     
   END SUBROUTINE read_electron_photon

 

  END  MODULE read_module

 
