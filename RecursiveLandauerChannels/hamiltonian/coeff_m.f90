
! Copyright (C) 2006 LEPES-CNRS Grenoble
!               2007 Institut Neel CNRS/UJF Grenoble
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!*********************************************
   MODULE coeff_module
!*********************************************
   USE parameters,      ONLY : nstrx
   USE kinds
   USE constants,            ONLY : CZERO, CONE
   USE identity_module,  ONLY :  orbital_def
   USE data_module, ONLY : orb_id, n_id

   IMPLICIT NONE
   PRIVATE 
   SAVE


   PUBLIC :: find_coeff
   PUBLIC :: find_coeff_onsite

   CONTAINS



!*********************************************
 SUBROUTINE find_coeff(name1, name2, rank, ene)
!*********************************************
  ! INPUTS
  CHARACTER(2), INTENT(in) :: name1
  CHARACTER(2), INTENT(in) :: name2
  INTEGER, INTENT(in) :: rank
  ! outputs
  COMPLEX(dbl),INTENT(out) :: ene
  ! local variables
  CHARACTER(10)  :: subname='find_coeff'
  LOGICAL :: find
  INTEGER :: i_name1, index_orb_id
   !
   ! end of declarations
   !

!
!----------------------------------------
! main Body
!----------------------------------------
!


    ! test input
    IF ((rank>3) .OR. (rank<1)) CALL errore (subname, 'invalid rank' , 1 )


  find=.FALSE.

   !find index1
    DO i_name1=1, n_id

        IF (TRIM(orb_id(i_name1)%name) == TRIM(name1) ) THEN
                IF (find) CALL errore (subname, 'problem in type def'//TRIM(name1) , 2 )
                index_orb_id=i_name1
                find=.TRUE.
        ENDIF
   ENDDO
   IF (.NOT.find) CALL errore (subname, 'invalid type 1= '//TRIM(name1) , 3 )

   !find hopping

            SELECT CASE( TRIM(name2) )
            CASE( "pz" )
                ene = orb_id(index_orb_id)%pz(rank)
            CASE( "px")
                ene = orb_id(index_orb_id)%px(rank)

            CASE DEFAULT
                CALL errore(subname, 'invalid type 2= '//TRIM(name2), 5 )
            END SELECT
            !
END SUBROUTINE find_coeff

!*********************************************
 SUBROUTINE find_coeff_onsite(name1, ene)
!*********************************************
  ! INPUTS
  CHARACTER(2), INTENT(in) :: name1
  ! outputs
  COMPLEX(dbl),INTENT(out) :: ene
  ! local variables
  CHARACTER(17)  :: subname='find_coeff_onsite'
  LOGICAL :: find
  INTEGER :: i_name1, index_orb_id
   !
   ! end of declarations
   !

!
!----------------------------------------
! main Body
!----------------------------------------
!


  
  find=.FALSE.

   !find index1
    DO i_name1=1, n_id

        IF (TRIM(orb_id(i_name1)%name) == TRIM(name1) ) THEN
                IF (find) CALL errore (subname, 'problem in type def'//TRIM(name1) , 2 )
                index_orb_id=i_name1
                find=.TRUE.
        ENDIF
   ENDDO
   IF (.NOT.find) CALL errore (subname, 'invalid type 1= '//TRIM(name1) , 3 )

   !find onsite

   ene = orb_id(index_orb_id)%onsite

 !
            !
END SUBROUTINE find_coeff_onsite

END MODULE coeff_module