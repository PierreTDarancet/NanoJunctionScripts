 
!
!      Copyright (C) 2005 WanT Group
!
!      This file is distributed under the terms of the
!      GNU General Public License. See the file `License'
!      in the root directory of the present distribution,
!      or http://www.gnu.org/copyleft/gpl.txt .
!
!*********************************************
   MODULE hamiltonian_module
!*********************************************
   USE parameters,      ONLY : nstrx
    USE dim_module,        ONLY : dimwann
   USE kinds
   USE constants,            ONLY : CZERO, CONE
   USE io_global_module,     ONLY : stdin
 !  USE kpoints_module,     ONLY :  nrtot_par
   USE iotk_module



   IMPLICIT NONE
   PRIVATE 
   SAVE

   COMPLEX(dbl), ALLOCATABLE ::  H_total
   PUBLIC :: H_total
   PUBLIC :: hamiltonian_init

CONTAINS

!*******************************************************************
   SUBROUTINE hamiltonian_init
   !*******************************************************************
   !
   ! Initialize hamiltonian

   !
   IMPLICIT NONE

   ! 
   ! input variables
   !
   !
   ! local variables
   !
   CHARACTER(16) :: subname="hamiltonian_init"
   COMPLEX(dbl), ALLOCATABLE :: aux(:,:)
   INTEGER       :: i, ierr

   !
   ! end of declarations
   !

!
!----------------------------------------
! main Body
!----------------------------------------
!

   !
   ! allocations
   !
   !
   CALL hamiltonian_allocate()
   !

   ALLOCATE( aux(dimwann,dimwann), STAT=ierr)
      IF ( ierr/=0 ) CALL errore(subname,'allocating aux',ABS(ierr))

   !
   ! open the IOTK tag
   !
   CALL iotk_scan_begin( stdin, 'HAMILTONIAN_DATA', IERR=ierr )
      IF (ierr/=0) CALL errore(subname,'searching HAMILTONIAN_DATA',ABS(ierr))


   !
   ! read basic quantities
   !
   CALL read_matrix( datafile_H, 'H_total', dimwan, dimwan, aux, dimx, dimx)
   !


   CALL iotk_scan_end( stdin, 'HAMILTONIAN_DATA', IERR=ierr )
      IF (ierr/=0) CALL errore(subname,'searching end for HAMILTONIAN_DATA',ABS(ierr))


   test_hamiltonian()

   !
   ! local cleaning
   !
   DEALLOCATE( aux, STAT=ierr)
      IF ( ierr/=0 ) CALL errore(subname,'deallocating aux',ABS(ierr))

END SUBROUTINE hamiltonian_init






!*******************************************************************
   SUBROUTINE hamiltonian_allocate
   !*******************************************************************
      CHARACTER(20)      :: subname="hamiltonian_allocate"
      INTEGER :: ierr


END SUBROUTINE hamiltonian_allocate
!**********************************************************
   SUBROUTINE hamiltonian_deallocate()
   !**********************************************************
   IMPLICIT NONE
      CHARACTER(22)      :: subname="hamiltonian_deallocate"
      INTEGER :: ierr

!       DEALLOCATE ( h00_L, h01_L, STAT=ierr )
!            IF( ierr /=0 ) CALL errore(subname, 'deallocating h00_L, h01_L', 1 )
   !
   ! read basic quantities
   !

  DEALLOCATE(,  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'deallocating H24',ABS(ierr))
   !
   END SUBROUTINE hamiltonian_deallocate


END MODULE hamiltonian_module


