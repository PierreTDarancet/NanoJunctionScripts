!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!               2007 Institut Neel CNRS/UJF Grenoble
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!*******************************************************************
   SUBROUTINE construct_hamiltonian
   !*******************************************************************
   USE kinds
   USE hamiltonian_module,   ONLY : hamiltonian_allocate, hamiltonian_alloc, &
                                    init_diagonal, init_off_diagonal, &
                                    init_hopping_hamiltonian
   USE distance_module,      ONLY : distance_allocate, distance_alloc, &
                                    init_metric, &
                                    init_metric_LC, init_metric_CR
   IMPLICIT NONE

   CHARACTER(21)      :: subname="construct_hamiltonian"
   INTEGER :: ierr


     hamiltonian_alloc=.FALSE.

     distance_alloc=.FALSE.

     CALL hamiltonian_allocate()

     CALL distance_allocate()

     CALL init_metric() 

     CALL init_metric_CR()

     CALL init_metric_LC()

     CALL init_diagonal()

     CALL init_off_diagonal()

     CALL init_hopping_hamiltonian()


END SUBROUTINE construct_hamiltonian




