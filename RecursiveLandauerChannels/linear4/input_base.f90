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
   MODULE input_base_module
!********************************************
   USE kinds, ONLY : dbl
   USE constants, ONLY : ZERO
   USE io_module, ONLY : stdout
   IMPLICIT NONE
   PRIVATE
!
! This module contains basic routines to read
! cards from input
!
! routines in this module:
! SUBROUTINE read_cards(unit)
! SUBROUTINE card_wannier_centers(unit,line)
! 

   PUBLIC :: read_cards
   PUBLIC :: card_orb


CONTAINS

!
! subroutines
!

!**********************************************************
   SUBROUTINE read_cards(unit)
   !**********************************************************
   !
   ! read the following cards:
   ! - WANNIER_CENTERS
   !
   USE parser_module, ONLY: read_line, capital
   IMPLICIT NONE
      INTEGER, INTENT(in):: unit

      CHARACTER(LEN=256) :: input_line
      CHARACTER(LEN=80)  :: card
      LOGICAL            :: lend, lstop
      LOGICAL            :: wannier_centers_found
      INTEGER :: i
      !

      !
      lstop = .FALSE.
      wannier_centers_found = .FALSE.
 100  CALL read_line(unit, input_line, END_OF_FILE=lend )
      !
      IF( lend .OR. lstop ) GO TO 120
      IF( input_line == ' ' .OR. input_line(1:1) == '#' ) GO TO 100
      !
      READ (input_line, *) card
       
      DO i = 1, LEN_TRIM( input_line )
         input_line( i : i ) = capital( input_line( i : i ) )
      ENDDO
       
      !
      IF ( TRIM(card) == 'ORB_INFO' ) THEN
         !
         CALL card_orb(unit,input_line)
         lstop = .TRUE.
         wannier_centers_found = .TRUE.
         !
      ELSE
         !
         WRITE( stdout,'(A)') 'Warning: card '//TRIM(input_line)//' ignored'
         !
      ENDIF
      !
      ! ...     END OF LOOP ... !
      !
      GOTO 100
      !
120   CONTINUE
      !

      IF ( .NOT. wannier_centers_found ) &
            CALL errore('read_cards','Card ORB_INFO not found',1)
    RETURN
  END SUBROUTINE read_cards


!**********************************************************
   SUBROUTINE card_orb(unit, input_line )
   !**********************************************************
      USE parser_module,  ONLY : read_line, matches, change_case
      USE orbital_module, ONLY : list => orb
      USE identity_module, ONLY : orbitale
      IMPLICIT NONE

      !
      INTEGER,            INTENT(in) :: unit
      CHARACTER(LEN=256), INTENT(in) :: input_line

      LOGICAL, SAVE      :: tread = .FALSE.
      CHARACTER(LEN=256) :: tmp_line
      CHARACTER(LEN=8) :: subname='card_orb'
      INTEGER            :: dim
      INTEGER            :: iwann, ierr
      CHARACTER(LEN=10)  :: adum, units
      !
      !
      IF ( tread ) CALL errore( 'card_orb', ' two occurrences ', 2 )

      dim = SIZE( list )


      !
      ! through the trial centers
      !
      DO iwann = 1, dim

           !
           ! ... init center
           !CALL trial_center_init(list(iwann))
           list(iwann)%coord(:) = ZERO
           IF (list(iwann)%init)       CALL errore(subname,'orb already init' , iwann )

           CALL read_line(unit, tmp_line )
           READ(tmp_line,*, IOSTAT=ierr) list(iwann)%nature
           IF (ierr/=0) CALL errore(subname,'reading line I',ABS(ierr))
           !
           CALL change_case(list(iwann)%nature,'lower')
 
           !
           READ(tmp_line,*, IOSTAT=ierr) adum, list(iwann)%coord(1:3)

              IF (ierr/=0) CALL errore(subname,'reading line II',ABS(ierr))

           !
           !
           list(iwann)%init = .TRUE.

       ENDDO

  END SUBROUTINE card_orb

END MODULE input_base_module

