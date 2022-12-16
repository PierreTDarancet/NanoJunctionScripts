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
   MODULE  T_output_module
   !***********************************************
  USE kinds, ONLY : dbl
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI, EPS_m3, EPS_m5

   PRIVATE 
   SAVE
!
   IMPLICIT NONE

   ! Public
   ! Private
   ! Default Values

   ! Other Module Variables
   CHARACTER(15) :: modulename="T_output_module"
   ! Public variables

   ! Public routines:
   PUBLIC                 :: write_output
   ! private routines
 
  CONTAINS 
!***********************************************
   SUBROUTINE write_output(Current(nalpha,nbias), Photocurrent(nphoton,nalpha,nbias), Photocurrentmax(nalpha,nbias), RR(nalpha, (INT(nbias/2))),  nalpha, alphagrid(nalpha),  nbias, biasgrid(bbias), nphotonfrequency, photonfrequency(nphotonfrequency))
   !***********************************************
   ! read the inp file located in the local directory
  IMPLICIT NONE
      ! Local
       CHARACTER(12)      :: subname="write_output"
 
       INTEGER :: ierr
 

        IF( dimx < 1 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in Lead dimension", " dim Lead = ", dimx
            STOP
         ENDIF 
  END MODULE  T_output_module
   USE constants,                          ONLY : PI, ZERO, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   USE T_control_module,                   ONLY : nprint, calculate_RR, debug_level, print_RR, print_photocurrent, print !ok
   USE T_output_module,                    ONLY : write_output

    END SUBROUTINE write_output

  PRINT*, "Writing output files"

  CALL write_output(Current(:,:), Photocurrent(:,:,:), Photocurrentmax(:,:), RR(:,:),  nalpha, alphagrid(:),  nbias, biasgrid(:), nphotonfrequency, photonfrequency(:))

 PRINT*, "Writing results"


   OPEN ( 12, FILE='Current.dat', FORM='formatted' )
   DO ialpha = 1, nalpha
	DO ibias = 1, nbias
		Current(ialpha,ibias) = Current(ialpha,ibias) * 6623617.82 / (27.2113845)   !!! nA
		WRITE (12, '(3(f15.9))' ) alphagrid(ialpha) , biasgrid(ibias),   Current(ialpha,ibias)
	ENDDO
        WRITE (12, '(2(e15.9))' )
   ENDDO
   CLOSE( 12 )

!   IF ( photoncouplingconstant  =/ ZERO)
!   DO ialpha = 1, nalpha
!	DO ibias = 1, nbias
!		Current(ialpha,ibias) = Current(ialpha,ibias) * 6623617.82 / (27.2113845)   !!! nA
!		WRITE (12, '(3(f15.9))' ) alphagrid(ialpha) , biasgrid(ibias),   Current(ialpha,ibias)
!	ENDDO
!        WRITE (12, '(2(e15.9))' )
!   ENDDO


   ENDIF     

   IF (calculate_RR )  THEN 
   OPEN ( 13, FILE='RectificationRatio.dat', FORM='formatted' )
   DO ivoc = 1, nvoc
      DO ibias = 1, INT(nbias/2)
       WRITE (13, '(3(f15.9))' ) work_scal2, ABS(biasgrid(ibias)),   RR(ialpha,ibias)
      ENDDO
       WRITE (13, '(2(e15.9))' )

   ENDDO
   CLOSE( 13 )
   ENDIF

   IF ( write_biasgrid) THEN
   ENDIF

   IF ( alphagrid) THEN
   ENDIF


!
!...  free memory
!

  PRINT*, "Deallocate Hamiltonian"
  CALL hamiltonian_deallocate()
  CALL gamma_deallocate()
  CALL occupation_deallocate()

  PRINT*, "Deallocate grids"
  CALL alphagrid_deallocate()
  CALL egrid_deallocate()
  CALL biasgrid_deallocate()
  CALL photongrid_deallocate()
  !

  PRINT*, "Deallocate transport variables"

    DEALLOCATE ( Photocurrentmax, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating Photocurrentmax"
            STOP
         ENDIF 

    DEALLOCATE ( Photocurrent, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating Photocurrent"
            STOP
         ENDIF 
     DEALLOCATE ( Current, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating current"
            STOP
         ENDIF 
     DEALLOCATE ( RR, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating rectification ratio"
            STOP
         ENDIF 

     DEALLOCATE ( omegamax, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating omegamax"
            STOP
         ENDIF 
    DEALLOCATE ( Photocurrent_aux, STAT=ierr )
         IF( ierr /=0 ) THEN 
            PRINT*, "Error deallocating Photocurrent_aux"
            STOP
         ENDIF 

   

 !    DEALLOCATE ( Emittedlight, STAT=ierr )
 !        IF( ierr /=0 ) THEN 
 !           PRINT*, "Error deallocating emitted light"
 !           STOP
 !        ENDIF 
 PRINT*, "End"


END PROGRAM conductor
