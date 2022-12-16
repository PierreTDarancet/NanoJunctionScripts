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
   MODULE T_egrid_module
   !***********************************************
  USE kinds, ONLY : dbl
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI
   IMPLICIT NONE
   PRIVATE 
   SAVE
!


   ! Public
   INTEGER                :: nemax     ! dimension of the maximum energy grid (input)
   INTEGER                :: ne        ! dimension of the energy grid for each bias
   REAL(dbl)              :: de        
   REAL(dbl)              :: deref        ! for debug purpose only   
   REAL(dbl), ALLOCATABLE :: egrid(:)  ! grid values
   REAL(dbl)              :: delta     ! i\delta for GFs
   REAL(dbl)              :: emin      !  Minimum energy (input) 
                                       ! (Should be automatized) like H00_C(i,i) min - biasmin - photon_frequency - a bit more
   REAL(dbl)              :: emax      !  Maximum energy (input) 
                                       ! (Should be automatized) like H00_C(i,i) max - biasmax + photon_frequency + a bit more

   ! Private
   LOGICAL :: alloc = .FALSE.

   ! Public variables
   PUBLIC                 :: nemax 
   PUBLIC                 :: ne 
   PUBLIC                 :: de
   PUBLIC                 :: deref ! for debug purpose only 
   PUBLIC                 :: egrid
   PUBLIC                 :: delta
   PUBLIC                 :: emin
   PUBLIC                 :: emax

   ! Public routines:
   PUBLIC                 :: egrid_allocate
   PUBLIC                 :: egrid_deallocate
   PUBLIC                 :: egrid_init

   CONTAINS 
!***********************************************
   SUBROUTINE egrid_allocate
   !***********************************************
   IMPLICIT NONE

       CHARACTER(14)      :: subname="egrid_allocate"
       INTEGER  :: ie, ierr
       !
        IF( alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Already allocated"
            STOP
         ENDIF 

        IF( nemax <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in nemax definition", " nemax = ", nemax
            STOP
         ENDIF 

       ALLOCATE( egrid(nemax), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating egrid"
            STOP
       ENDIF 
       egrid(:) = ZERO
       deref = (emax - emin) / REAL(nemax-1, dbl)
       alloc = .TRUE.


  END SUBROUTINE egrid_allocate
!***********************************************
   SUBROUTINE egrid_deallocate 
   !***********************************************
  IMPLICIT NONE
       CHARACTER(16)      :: subname="egrid_deallocate"
       INTEGER :: ierr


       IF ( ALLOCATED(egrid) ) THEN
            DEALLOCATE(egrid, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating egrid"
               STOP
           ENDIF 
       ENDIF
       alloc = .FALSE.

   END SUBROUTINE egrid_deallocate

!***********************************************
   SUBROUTINE egrid_init()
   !***********************************************
  IMPLICIT NONE

   CHARACTER(10)      :: subname="egrid_init"
   INTEGER      :: ie, ierr
   REAL(dbl)    :: eminlocal
   REAL(dbl)    :: emaxlocal
   REAL(dbl)    :: biasref

 
       egrid(:) = ZERO
       emaxlocal = emax! This formulation should allow to deal with asymmetric bias cases: 
       eminlocal = emin 

       !
       IF( emaxlocal <  eminlocal ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error emaxlocal is smaller than eminlocal"
            PRINT*, "========================================="
            PRINT*, "=============STOPPING===================="
            PRINT*, "========================================="
            PRINT*, "Entered    emax     = ", emax, "  eV"
            PRINT*, "Entered    emin     = ", emin, "  eV"
            PRINT*, "Calculated emax     = ", emaxlocal, "  eV"
            PRINT*, "Calculated emin     = ", eminlocal, "  eV"
            PRINT*, "========================================="
            STOP
       ENDIF 
   
       ne = INT( (  (emaxlocal - eminlocal) / (emax - emin)  )  * nemax ) 
       IF ( ne <= 1 ) THEN 
            PRINT*, "========================================="
            PRINT*, "==============WARNING===================="
            PRINT*, "========================================="
            PRINT*, "         ne will be set to 1             "
            PRINT*, "========================================="
            PRINT*, "============Continuing==================="
            PRINT*, "The results are probably going to be affected..."

            PRINT*, "Entered    emax     = ", emax, "  eV"
            PRINT*, "Entered    emin     = ", emin, "  eV"
            PRINT*, "Calculated emax     = ", emaxlocal, "  eV"
            PRINT*, "Calculated emin     = ", eminlocal, "  eV"
            PRINT*, "Nemax         = ", nemax
            PRINT*, "Ne         = ", ne
            ne = 1
            PRINT*, "========================================="
            PRINT*, "============Continuing==================="
       ENDIF
       IF( ne > nemax ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error ne of iteration is greater than nemax"
            PRINT*, "========================================="
            PRINT*, "=============STOPPING===================="
            PRINT*, "========================================="

            PRINT*, "Entered    emax     = ", emax, "  eV"
            PRINT*, "Entered    emin     = ", emin, "  eV"
            PRINT*, "Calculated emax     = ", emaxlocal, "  eV"
            PRINT*, "Calculated emin     = ", eminlocal, "  eV"
            PRINT*, "Nemax         = ", nemax
            PRINT*, "Calculated ne = ", ne
            PRINT*, "========================================="
            STOP
       ENDIF 


          de = (emaxlocal - eminlocal) / REAL(ne-1, dbl)
          DO ie = 1, ne
              egrid(ie) = eminlocal + REAL(ie -1, dbl) * de
          ENDDO
 
          ! check this part another time
   END SUBROUTINE egrid_init

  END  MODULE T_egrid_module

 
