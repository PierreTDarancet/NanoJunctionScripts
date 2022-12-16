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
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI, EPS_m5
   IMPLICIT NONE
   PRIVATE 
   SAVE
!


   ! Public
   INTEGER                :: ne_green    ! dimension of the maximum energy grid for green functions
   INTEGER                :: ne        ! dimension of the energy grid for transmission (input)
   REAL(dbl)              :: de        
   REAL(dbl)              :: fermi_energy
   REAL(dbl)              :: electronic_temperature
   REAL(dbl), ALLOCATABLE :: egrid(:)  ! grid values
   INTEGER,   ALLOCATABLE :: egridplusomega(:)  ! grid values
   INTEGER,   ALLOCATABLE :: egridminusomega(:)  ! grid values
   REAL(dbl), ALLOCATABLE :: green_egrid(:)  ! grid values
   INTEGER,   ALLOCATABLE :: grid_condtogreen(:)  ! grid values
   REAL(dbl)              :: delta     ! i\delta for GFs
   REAL(dbl)              :: emin      !  Minimum energy (input) 

   REAL(dbl)              :: emax      !  Maximum energy (input) 

   REAL(dbl)              :: green_emin      !  Minimum energy (input) 

   REAL(dbl)              :: green_emax      !  Maximum energy (input) 

   ! Private
   LOGICAL :: alloc = .FALSE.

   ! Public variables
   PUBLIC                 :: ne_green
   PUBLIC                 :: ne 
   PUBLIC                 :: de

   PUBLIC                 :: egrid
   PUBLIC                 :: delta
   PUBLIC                 :: emin
   PUBLIC                 :: emax
   PUBLIC                 :: green_emin
   PUBLIC                 :: green_emax

   PUBLIC                 :: fermi_energy
   PUBLIC                 :: electronic_temperature

   PUBLIC                 :: green_egrid
   PUBLIC                 :: egridplusomega
   PUBLIC                 :: egridminusomega
   PUBLIC                 :: grid_condtogreen
   ! Public routines:
   PUBLIC                 :: egrid_allocate
   PUBLIC                 :: egrid_deallocate
   PUBLIC                 :: egrid_init
   PUBLIC                 :: green_egrid_allocate
   PUBLIC                 :: green_egrid_deallocate
   PUBLIC                 :: green_egrid_init
   PUBLIC                 :: egridshifted_init 
   PUBLIC                 :: evaluate_negreen

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

        IF( ne <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in ne definition", " ne = ", ne
            STOP
         ENDIF 

       ALLOCATE( egrid(ne), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating egrid"
            STOP
       ENDIF 
       egrid(:) = ZERO
       alloc = .TRUE.


  END SUBROUTINE egrid_allocate
!***********************************************
   SUBROUTINE green_egrid_allocate
   !***********************************************
   IMPLICIT NONE

       CHARACTER(20)      :: subname="green_egrid_allocate"
       INTEGER  :: ie, ierr
       !
   
        IF( ne_green <= 0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in ne definition", " ne_green = ", ne_green
            STOP
         ENDIF 

       ALLOCATE( green_egrid(ne_green), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating green_egrid"
            STOP
       ENDIF 
       ALLOCATE( grid_condtogreen(ne), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating grid_condtogreen"
            STOP
       ENDIF 
       ALLOCATE(  egridplusomega(ne_green), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating egridplusomega"
            STOP
       ENDIF 
       ALLOCATE( egridminusomega(ne_green), STAT=ierr )
       IF( ierr /=0 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error allocating egridminusomega"
            STOP
       ENDIF 

       grid_condtogreen(:) = 0
       green_egrid(:) = ZERO
       egridminusomega(:) = 0
       egridplusomega(:) = 0
       delta=EPS_m5
  END SUBROUTINE green_egrid_allocate
!***********************************************
   SUBROUTINE evaluate_negreen(omega_photon)
   !***********************************************
   IMPLICIT NONE
       REAL(dbl), INTENT(in)    :: omega_photon
       CHARACTER(16)      :: subname="evaluate_negreen"
       INTEGER  :: ie, ierr
       !

       de = (emax - emin) / REAL(ne-1, dbl)

       green_emin = emin - ABS(omega_photon)
       green_emax = emax + ABS(omega_photon)

       ne_green =  1+ INT( (green_emax - green_emin) /  de  )

  END SUBROUTINE evaluate_negreen
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
   SUBROUTINE green_egrid_deallocate
   !***********************************************
   IMPLICIT NONE

       CHARACTER(22)      :: subname="green_egrid_deallocate"
       INTEGER  :: ie, ierr
       !
   
 
       IF ( ALLOCATED(green_egrid) ) THEN
            DEALLOCATE(green_egrid, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating green_egrid"
               STOP
           ENDIF 
       ENDIF
      IF ( ALLOCATED(grid_condtogreen) ) THEN
            DEALLOCATE(grid_condtogreen, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating grid_condtogreen"
               STOP
           ENDIF 
       ENDIF
      IF ( ALLOCATED(egridplusomega) ) THEN
            DEALLOCATE(egridplusomega, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating egridplusomega"
               STOP
           ENDIF 
       ENDIF
      IF ( ALLOCATED(egridminusomega) ) THEN
            DEALLOCATE(egridminusomega, STAT=ierr)
            IF( ierr /=0 ) THEN 
               PRINT*, "Error in Routine ",  TRIM(subname) 
               PRINT*, "Error deallocating egridminusomega"
               STOP
           ENDIF 
       ENDIF

  END SUBROUTINE green_egrid_deallocate

!***********************************************
   SUBROUTINE egrid_init
   !***********************************************
  IMPLICIT NONE

   CHARACTER(10)      :: subname="egrid_init"
   INTEGER      :: ie, ierr
 
       IF( .NOT. alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Not allocated"
            STOP
         ENDIF 

          de = (emax - emin) / REAL(ne-1, dbl)
          DO ie = 1, ne
              egrid(ie) = emin + REAL(ie -1, dbl) * de
          ENDDO
 
   END SUBROUTINE egrid_init

!***********************************************
   SUBROUTINE green_egrid_init
   !***********************************************
  IMPLICIT NONE

   CHARACTER(16)      :: subname="green_egrid_init"
   INTEGER      :: ie, ierr, ie_test
 
       IF( .NOT. alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Not allocated"
            STOP
         ENDIF 
          ie_test= 1
          de = (emax - emin) / REAL(ne-1, dbl)
          DO ie = 1, ne_green
              green_egrid(ie) = green_emin + REAL(ie -1, dbl) * de
              IF ((green_egrid(ie) >= emin).AND.(green_egrid(ie) <= emax).AND. (ie_test<=(ne+1))) THEN
                        grid_condtogreen(ie_test)=ie
                        ie_test=ie_test+1
              ENDIF
          ENDDO


   END SUBROUTINE green_egrid_init
!***********************************************
   SUBROUTINE  egridshifted_init(omega_photon)
   !***********************************************
  IMPLICIT NONE
       REAL(dbl), INTENT(in)    :: omega_photon
   REAL(dbl)    ::  ene_min, ene_max
   CHARACTER(17)      :: subname="egridshifted_init"
   INTEGER      :: ie, ierr, ie_test


       IF( .NOT. alloc ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Not allocated"
            STOP
         ENDIF 


          ie_test = INT(omega_photon/de)

          DO ie = 1, ne_green
                  IF ((ie +  ie_test) <ne_green) THEN
                     egridplusomega(ie)=ie + ie_test 
                  ELSE
                     egridplusomega(ie)=ne_green
                  ENDIF
          ENDDO


          DO ie = 1, ne_green
                  IF ((ie -  ie_test) > 0) THEN
                   egridminusomega(ie)=ie - ie_test 
                  ELSE
                   egridminusomega(ie)=1
                  ENDIF
          ENDDO




   END SUBROUTINE  egridshifted_init


  END  MODULE T_egrid_module

 
