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
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI, EPS_m3

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
    INTEGER, INTENT(in)       :: dimC
    COMPLEX(dbl), INTENT(out) :: sigma_ee_mol_r(dimC,dimC)

    CHARACTER(3), INTENT(in)         ::  method
    CHARACTER(15)         ::     fname_in="SIGMA_MOL_GW.in"
    CHARACTER(19)      :: subname="read_mol_selfenergy"
    INTEGER :: ierr, i_dim
    REAL(dbl) :: eigenvalue

        IF ( (TRIM(method)=="EXC") .OR. (TRIM(method)=="SIG")) THEN
		!
	 	PRINT*, "     ... OPENING File:", fname_in
		OPEN(100, FILE=TRIM(fname_in),FORM='formatted')
		!
		sigma_ee_mol_r(:,:) = CZERO
	  	PRINT*, "     ...... READING GW CORRECTIONS"
	  	PRINT*, "     ...... WARNING: Only diagonal corrections are accepted at this point"
		DO i_dim=1,dimC
		   !
		   READ(100,"(f11.6)")  eigenvalue 
                   sigma_ee_mol_r(i_dim,i_dim) =eigenvalue
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
   SUBROUTINE read_imagecharge_selfenergy(sigma_ee_contact_r,dimC, method)
   !***********************************************
  IMPLICIT NONE
    INTEGER, INTENT(in)       :: dimC
    COMPLEX(dbl), INTENT(out) :: sigma_ee_contact_r(dimC,dimC)
    CHARACTER(3), INTENT(in)         ::  method

    CHARACTER(19)         ::     fname_in="SIGMA_CONTACT_GW.in"
    CHARACTER(27)      :: subname="read_imagecharge_selfenergy"
    INTEGER :: ierr, i_dim
    REAL(dbl) :: eigenvalue

        IF ((TRIM(method)=="EXC" ).OR.(TRIM(method)=="SIG")) THEN
		!
	 	PRINT*, "     ... OPENING File:", fname_in
		OPEN(100, FILE=TRIM(fname_in),FORM='formatted')
		!
		sigma_ee_contact_r(:,:)= CZERO
	  	PRINT*, "     ...... READING IMAGE CHARGE CORRECTIONS"
	  	PRINT*, "     ...... WARNING: Only diagonal corrections are accepted at this point"
		DO i_dim=1,dimC
		   !
		   READ(100,"(f11.6)")  eigenvalue
                   sigma_ee_contact_r(i_dim,i_dim) = eigenvalue
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
    INTEGER, INTENT(in)       :: dimC
    COMPLEX(dbl), INTENT(out) :: sigma_L(ne_green,dimC,dimC)
    COMPLEX(dbl), INTENT(out) :: sigma_R(ne_green,dimC,dimC)
    REAL(dbl) :: sigr, sigi
    INTEGER, INTENT(in)       :: ne_green
    INTEGER, INTENT(in)       :: ik
    CHARACTER(100)         ::   fname_in 
    CHARACTER(20)      :: subname="read_lead_selfenergy"
    INTEGER :: ierr, i_dim, ie, i_row, i_col, i_test


        WRITE (fname_in, *) "SIGMA_L", ik, ".in"
        !
 	PRINT*, "     ... OPENING File:", TRIM(fname_in)
	OPEN(100, FILE=TRIM(fname_in),FORM='formatted')
	!
        sigma_L(:,:,:)= CZERO
  	PRINT*, "     ...... READING SIGMA_L"
  	PRINT*, "     ...... WARNING: The grid is assumed to be equal to the one used in ne_green"
  	PRINT*, "     ...... WARNING: SIGMA IS SUPPOSED TO BE DIAGONAL"
        i_test =0
        PRINT*, "Number of points", ne_green
        DO ie=1,ne_green
        !PRINT*, "Energy Point", ie
!        DO i_col=1,dimC
!           DO i_row=1,dimC
!           !
!           READ(100,"(f11.6)")  sigma_L(ie,i_row,i_col)
!           i_test = i_test + 1
!           !
!           ENDDO
!        ENDDO
        sigma_L(ie,:,:)= CZERO
        DO i_col=1,dimC

           !
           READ(100,*)  sigr, sigi
            sigma_L(ie,i_col,i_col) = CMPLX( sigr, sigi)
           i_test = i_test + 1
           !

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
  	PRINT*, "     ...... WARNING: SIGMA IS SUPPOSED TO BE DIAGONAL"
        i_test =0
        DO ie=1,ne_green
        !PRINT*, "Energy Point", ie
        sigma_R(ie,:,:)=CZERO
        DO i_col=1,dimC
   !        DO i_row=1,dimC
           !
           READ(100,*) sigr, sigi
           sigma_R(ie,i_col,i_col) = CMPLX( sigr, sigi)
           i_test = i_test + 1
           !
  !         ENDDO
        ENDDO
        ENDDO
        !

  	PRINT*, "     ... CLOSING File:", fname_in
        CLOSE(100)


   END SUBROUTINE read_lead_selfenergy

!***********************************************
   SUBROUTINE read_lead_gamma( gamma_L, gamma_R, ik,  ne_green,dimC)
   !***********************************************
  IMPLICIT NONE
      INTEGER, INTENT(in)       :: dimC
    COMPLEX(dbl), INTENT(out) :: gamma_L(ne_green,dimC,dimC)
    COMPLEX(dbl), INTENT(out) :: gamma_R(ne_green,dimC,dimC)

    INTEGER, INTENT(in)       :: ne_green
    INTEGER, INTENT(in)       :: ik
    CHARACTER(100)         ::   fname_in 
    CHARACTER(15)      :: subname="read_lead_gamma"
    INTEGER :: ierr, i_dim, ie, i_row, i_col, i_test
    REAL(dbl) :: gamma_aux(dimC,dimC)


        WRITE (fname_in, *) "GAMMA_L", ik, ".in"
        !
 	PRINT*, "     ... OPENING File:", TRIM(fname_in)
	OPEN(100, FILE=TRIM(fname_in),FORM='formatted')
	!
        gamma_L(:,:,:)= CZERO
  	PRINT*, "     ...... READING GAMMA_L"
  	PRINT*, "     ...... WARNING: The grid is assumed to be equal to the one used in ne_green"
  	PRINT*, "     ...... WARNING: GAMMA IS SUPPOSED TO BE DIAGONAL"
        i_test =0
        DO ie=1,ne_green
		!PRINT*, "Energy Point", ie
!		DO i_col=1,dimC
!		   DO i_row=1,dimC
			   !
!			   READ(100,"(f11.6)")  gamma_aux(i_row,i_col)
!			    i_test = i_test + 1
!			   !
!		   ENDDO
!		ENDDO

                gamma_aux(:,:)=ZERO
		DO i_col=1,dimC

			   !
			   READ(100,*)  gamma_aux(i_col,i_col)
			    i_test = i_test + 1
			   !

		ENDDO

		gamma_L(ie,:,:) =  gamma_aux(:,:)
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
  	PRINT*, "     ...... WARNING: GAMMA IS SUPPOSED TO BE DIAGONAL"
        i_test =0
        DO ie=1,ne_green
		!PRINT*, "Energy Point", ie
                gamma_aux(:,:)=ZERO
!		DO i_col=1,dimC
!		   DO i_row=1,dimC
!			   !
!			   READ(100,"(f11.6)")  gamma_aux(i_row,i_col)
!			   i_test = i_test + 1
!			   !
!	           ENDDO
!		ENDDO
		DO i_col=1,dimC

			   !
			   READ(100,*)  gamma_aux(i_col,i_col)
			    i_test = i_test + 1
			   !

		ENDDO




		gamma_R(ie,:,:)=  gamma_aux(:,:)
        ENDDO
        !

  	PRINT*, "     ... CLOSING File:", fname_in
        CLOSE(100)



  END SUBROUTINE read_lead_gamma

!***********************************************
   SUBROUTINE read_electron_photon( dipole_matrix, dimC, method, ihomo, h_00)
   !***********************************************
  IMPLICIT NONE

    INTEGER, INTENT(in)       :: dimC
    INTEGER, INTENT(in)       :: ihomo
    REAL(dbl), INTENT(out) :: dipole_matrix(dimC,dimC)
    COMPLEX(dbl), INTENT(in) :: h_00(dimC,dimC)
    CHARACTER(2), INTENT(in)  :: method

    CHARACTER(20)         ::     fname_in="PositionsAO.in"
    CHARACTER(20)         ::     fname="STATES_on_AO.in"
    CHARACTER(50)         ::     fname_out="Calculated_dipole_matrix.out"
    CHARACTER(50)         ::    chr
    CHARACTER(20)        :: subname="read_electron_photon"
    INTEGER :: ierr, i_row, i_col, i_st

    REAL(dbl), ALLOCATABLE :: position_ao(:,:), states(:,:)
    INTEGER :: iorb, norb, icheck, iorb2
    REAL(dbl) :: dipole_matrix_z(dimC,dimC), deltapos, posmin, posmax!,  dipole_matrix_x(dimC,dimC),  dipole_matrix_y(dimC,dimC) NOT YET implemented
    REAL(dbl) :: eigenvaldiff
        dipole_matrix(:,:) = ZERO

        IF ( (TRIM(method)=="AM") .OR. (TRIM(method)=="GN") ) THEN ! Aeberhard Morf PRB 77 (2008) - SCBA - Galperin Nitzan PRL 2005 with dipole transition/
		!
		PRINT*, "     ... OPENING File:", fname_in
		OPEN(100, FILE=TRIM(fname_in),FORM='formatted')
		PRINT*, "     ...... READING POSITIONS AO"
		!
                READ(100,*) norb
		PRINT*, "     ......... Number of AO:", norb

  	     	       ALLOCATE( position_ao(3,norb)  , STAT=ierr )
	       IF( ierr /=0 ) THEN 
		    PRINT*, "Error in Routine ",  TRIM(subname) 
		    PRINT*, "Error allocating  position_ao"
		    STOP
	       ENDIF 
               posmin=ZERO
               posmax=ZERO
               DO iorb=1, norb
                   READ(100,*)  position_ao(1,iorb), position_ao(2,iorb), position_ao(3,iorb)
                   IF ( position_ao(3,iorb) < posmin) posmin=position_ao(3,iorb)
                   IF ( position_ao(3,iorb) > posmax) posmax=position_ao(3,iorb)
                   !PRINT*, position_ao(1,iorb), position_ao(2,iorb), position_ao(3,iorb)
               ENDDO
               CLOSE(100)
        	PRINT*, "     ... CLOSING File:", fname_in
                ALLOCATE( states(dimC,norb)  , STAT=ierr )
       	        IF( ierr /=0 ) THEN 
 		    PRINT*, "Error in Routine ",  TRIM(subname) 
		    PRINT*, "Error allocating states"
		    STOP
	        ENDIF 
		PRINT*, "     ... OPENING File:", fname

		OPEN(100, FILE=TRIM(fname),FORM='formatted')
     		PRINT*, "     ...... READING STATES ON AO"
                DO i_st=1, dimC
		        READ(100,*)  chr
	     		PRINT*, "     .......... READING:", TRIM(chr)
		        DO iorb = 1, norb
		           READ(100,*) icheck, states(i_st,iorb)
		           IF (icheck /= iorb) THEN
		     		    PRINT*, "Error in Routine ",  TRIM(subname) 
		     		    PRINT*, "State,", i_st, "on AO ", iorb, "indice=", icheck
		                    STOP
                           ENDIF
		        ENDDO
                ENDDO
                CLOSE(100)
        	PRINT*, "     ... CLOSING File:", fname

		PRINT*, "     ... Calculating Dipole transition matrix:"
                deltapos = (posmax-posmin)/(REAL(norb) )
     		PRINT*, "     ...... dr=", deltapos

		DO i_col=1,dimC
		   DO i_row=i_col,dimC
		   !
                   dipole_matrix_z(i_row,i_col)= ZERO
                   IF (i_col == i_row) THEN
              		   dipole_matrix_z(i_row,i_col)= ZERO
                   ELSE IF ((i_col <= ihomo) .AND. (i_row <=ihomo) ) THEN
              		   dipole_matrix_z(i_row,i_col)= ZERO
                   ELSE IF ((i_col > ihomo) .AND. (i_row > ihomo) ) THEN
              		   dipole_matrix_z(i_row,i_col)= ZERO
                   ELSE
                       !PRINT*, "     ... Calculating Dipole transition for states:", i_row, i_col
                       DO iorb = 1,norb
                           !
                           DO iorb2 = 1,norb

                               !
                               IF (( position_ao(3,iorb) >=  position_ao(3,iorb2) - EPS_m3 ) .AND. (position_ao(3,iorb) <=  position_ao(3,iorb2) + EPS_m3 )) THEN
                                   dipole_matrix_z(i_row,i_col) = dipole_matrix_z(i_row,i_col) + deltapos*position_ao(3,iorb)*states(i_row,iorb)*states(i_col,iorb2)                                 
                                   !PRINT*, "Conv ", iorb, iorb2
                              ENDIF

                           ENDDO
                       ENDDO
                   ENDIF
		   !
		   ENDDO
		ENDDO
          	PRINT*, "     .......  Dipole transition matrix done"
		DO i_col=1,dimC
		   DO i_row=1,i_col
                   !
                   dipole_matrix_z(i_row,i_col) = dipole_matrix_z(i_col,i_row)
                   !
                   ENDDO
                ENDDO

                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SIMPLE CASE, for more general, implement dipole_matrix_y and others
		DO i_col=1,dimC
		   DO i_row=1,dimC
			   !
		           eigenvaldiff =ABS( REAL( h_00(i_row, i_row)) - REAL(h_00(i_col, i_col) ) )
			   dipole_matrix(i_row,i_col)= ABS (dipole_matrix_z(i_row,i_col) ) *  eigenvaldiff  /27.21138386
			   !
		   ENDDO
		ENDDO
		!
		PRINT*, "     ... CLOSING File:", fname_in
		!
                CLOSE(100)


                PRINT*, "     ... PRINTING Dipole matrix in file ", TRIM(fname_out)
		PRINT*, "     ... OPENING File:", fname_out

		OPEN(100, FILE=TRIM(fname_out),FORM='formatted')
                
                WRITE(100,*)  (dimC*dimC)
	        DO i_col=1,dimC
		   !
                   DO i_row=1, dimC
   		   WRITE(100,*) dipole_matrix(i_row,i_col)
                   ENDDO
		   !
		ENDDO


  		PRINT*, "     ... CLOSING File:", fname_out
		!
                CLOSE(100)
              

        ELSE IF (TRIM(method)=="G0") THEN !  Galperin Nitzan PRL 2005 without dipole transition -set to one
		PRINT*, "     ...... WARNING ---"
		PRINT*, "     ...... METHOD is G0 ---"
		PRINT*, "     ...... dipole transition is set to one off diagonal ---"
		DO i_col=1,dimC
		   DO i_row=1,dimC
		   !
                   IF ( ( i_row <= INT(dimC/2.00) ) .AND. ( i_col > INT(dimC/2.00)) ) THEN
             		   dipole_matrix(i_row,i_col) = ONE
                   ELSE IF  ( ( i_col <= INT(dimC/2.00) ) .AND. ( i_row > INT(dimC/2.00)) ) THEN
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


       IF ( ALLOCATED( position_ao  ) ) THEN
            DEALLOCATE( position_ao, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating  position_ao"
               STOP
           ENDIF 
       ENDIF
   



       IF ( ALLOCATED( states  ) ) THEN
            DEALLOCATE( states, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating  states"
               STOP
           ENDIF 
       ENDIF
   
   END SUBROUTINE read_electron_photon

 

  END  MODULE read_module

 
