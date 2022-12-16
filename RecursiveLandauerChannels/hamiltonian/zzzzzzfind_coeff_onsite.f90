
! Copyright (C) 2006 LEPES-CNRS Grenoble
!               2007 Institut Neel CNRS/UJF Grenoble
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!*********************************************
 SUBROUTINE find_coeff_onsite(name1, ene)
!*********************************************
   USE parameters,      ONLY : nstrx
   USE kinds
   USE constants,            ONLY : CZERO, CONE
   USE identity_module,  ONLY orb_id, orbital_def, n_id
   IMPLICIT NONE



  ! INPUTS

  CHARACTER(2), INTENT(in) :: name1


  ! outputs
  COMPLEX(dbl),INTENT(out) :: ene


  ! local variables
  CHARACTER(17)  :: subname='find_coeff_onsite'
  LOGICAL :: find
  INTEGER :: i_name1
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