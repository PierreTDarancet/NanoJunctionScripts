!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!               2007 Institut Neel CNRS/UJF Grenoble
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!*********************************************
   MODULE hamiltonian_module
!*********************************************
   USE parameters,      ONLY : nstrx
   USE kinds
   USE constants,            ONLY : CZERO, CONE, ZERO, EPS_m11
   USE io_module,     ONLY :  ham_unit, rec_ham_unit => aux3_unit, ioname
   USE iotk_module
   USE control_variable_module,    ONLY : use_second, use_third, datafile_H
   USE dim_variable_module,        ONLY : n_orb, limit_0,  limit_1, limit_2, limit_3
   USE identity_module,      ONLY : orbitale
   USE orbital_module,       ONLY : orb
   USE distance_module,      ONLY : distance, distance_LC, distance_CR
   USE coeff_module,         ONLY : find_coeff, find_coeff_onsite
   USE files_module, ONLY : file_open, file_close
   USE io_global_module,     ONLY : stdin, stdout

   IMPLICIT NONE
   PRIVATE 
   SAVE

! Variables
    COMPLEX(dbl), ALLOCATABLE :: hamiltonian(:,:)
!   
    COMPLEX(dbl), ALLOCATABLE :: hamiltonian_CR(:,:)
!
    COMPLEX(dbl), ALLOCATABLE :: hamiltonian_LC(:,:)
!

   LOGICAL :: hamiltonian_alloc
!

! Status

   PUBLIC :: hamiltonian
!
   PUBLIC :: hamiltonian_CR
!
   PUBLIC :: hamiltonian_LC
!
!
!

   PUBLIC:: hamiltonian_alloc
!



! functions


  PUBLIC :: hamiltonian_allocate
  PUBLIC :: hamiltonian_deallocate
  PUBLIC :: init_diagonal
  PUBLIC :: init_off_diagonal
  PUBLIC :: init_hopping_hamiltonian
  PUBLIC :: print_hamiltonian




CONTAINS


!*******************************************************************
   SUBROUTINE hamiltonian_allocate
   !*******************************************************************
      CHARACTER(20)      :: subname="hamiltonian_allocate"
      INTEGER :: ierr

   !
   ! allocate basic quantities
   !
   IF  (hamiltonian_alloc) CALL errore(subname,'hamiltonian already allocated',1)
!


   ALLOCATE( hamiltonian(n_orb, n_orb),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating hamiltonian',ABS(ierr))
   !
   ALLOCATE( hamiltonian_CR(n_orb, n_orb),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating hamiltonian_CR',ABS(ierr))
   !
   ALLOCATE( hamiltonian_LC(n_orb, n_orb),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating hamiltonian_LC',ABS(ierr))

   hamiltonian(:,:)=CZERO

   hamiltonian_CR(:,:)=CZERO

   hamiltonian_LC(:,:)=CZERO

    hamiltonian_alloc=.TRUE.

END SUBROUTINE hamiltonian_allocate

!**********************************************************
   SUBROUTINE hamiltonian_deallocate
   !**********************************************************
   IMPLICIT NONE
      CHARACTER(22)      :: subname="hamiltonian_deallocate"
      INTEGER :: ierr

    IF (.NOT. hamiltonian_alloc) CALL errore(subname,'hamiltonian not allocated',1)

     IF ( ALLOCATED( hamiltonian  ) ) THEN
           DEALLOCATE ( hamiltonian, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating hamiltonian', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( hamiltonian_CR  ) ) THEN
           DEALLOCATE ( hamiltonian_CR, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating hamiltonian_CR', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( hamiltonian_LC  ) ) THEN
           DEALLOCATE ( hamiltonian_LC, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating hamiltonian_LC', ABS(ierr) )
      ENDIF
 
     hamiltonian_alloc=.FALSE.

   END SUBROUTINE hamiltonian_deallocate

!*******************************************************************
   SUBROUTINE init_diagonal
   !*******************************************************************
      CHARACTER(13)      :: subname="init_diagonal"
      INTEGER :: ierr
      INTEGER :: i_orb, j_orb
      REAL(dbl) :: test
      COMPLEX(dbl) :: ene
      !
      !

      DO i_orb=1, n_orb
            test=ZERO
            !
            DO j_orb=1, n_orb
                   test = test + ABS(hamiltonian(i_orb,j_orb))
            ENDDO
            !

            IF (( test <= - EPS_m11 ).OR. ( test >= EPS_m11 ) ) &
                        CALL errore(subname,'off diagonal element are present ',i_orb)
           !
           !
      ENDDO


      DO i_orb=1, n_orb
            hamiltonian(i_orb,i_orb) = CZERO

            CALL find_coeff_onsite((TRIM(orb(i_orb)%nature)), ene)
            hamiltonian(i_orb,i_orb) = ene
      ENDDO

      ! test



END SUBROUTINE init_diagonal

!*******************************************************************
   SUBROUTINE init_off_diagonal
   !*******************************************************************
      CHARACTER(17)      :: subname="init_off_diagonal"
      INTEGER :: ierr
      INTEGER :: i_orb, j_orb
      REAL(dbl) :: test
      COMPLEX(dbl)  :: ene
      !
      DO i_orb=1, n_orb-1
            !
            DO j_orb=i_orb+1,n_orb
               !
               IF (i_orb == j_orb) CALL errore(subname,'i_orb and j_orb are supposed to be different',i_orb)
               !

               IF ((distance(i_orb,j_orb) <= limit_1) .AND. (distance(i_orb,j_orb) >= limit_0) ) THEN
                     CALL find_coeff((TRIM(orb(i_orb)%nature)), (TRIM(orb(j_orb)%nature)), 1, ene)
                     hamiltonian(i_orb,j_orb) = ene
                     hamiltonian(j_orb,i_orb) = CONJG (hamiltonian(i_orb,j_orb))
               !
               ELSE IF ((distance(i_orb,j_orb) <= limit_2) .AND. (distance(i_orb,j_orb) >= limit_1).AND. use_second ) THEN
                     CALL find_coeff((TRIM(orb(i_orb)%nature)), (TRIM(orb(j_orb)%nature)), 2, ene)

                     hamiltonian(i_orb,j_orb) = ene
                     hamiltonian(j_orb,i_orb) = CONJG (hamiltonian(i_orb,j_orb))
               !
               ELSE IF ((distance(i_orb,j_orb) <= limit_3) .AND. (distance(i_orb,j_orb) >= limit_2).AND. use_third )  THEN
                     CALL find_coeff((TRIM(orb(i_orb)%nature)), (TRIM(orb(j_orb)%nature)), 3, ene)

                     hamiltonian(i_orb,j_orb) = ene
                     hamiltonian(j_orb,i_orb) = CONJG (hamiltonian(i_orb,j_orb))
               !
               ELSE
                     hamiltonian(i_orb,j_orb) = CZERO
                     hamiltonian(j_orb,i_orb) = CONJG (hamiltonian(i_orb,j_orb))
               !
               ENDIF
               !

            ENDDO
            !
      ENDDO

      test=ZERO
      !
      DO i_orb=1, n_orb
            !
            DO j_orb=1, n_orb
                !
                test = test + ABS(hamiltonian(i_orb,j_orb) - hamiltonian(j_orb,i_orb))
                !
            ENDDO
      ENDDO

      IF (( test <= - EPS_m11 ).OR. ( test >=  EPS_m11 ) ) &
                        CALL errore(subname,'hamiltonian not hermitian',1)
      !


END SUBROUTINE init_off_diagonal
!*******************************************************************
   SUBROUTINE init_hopping_hamiltonian
   !*******************************************************************
      CHARACTER(24)      :: subname="init_hopping_hamiltonian"
      INTEGER :: ierr
      COMPLEX(dbl) :: ene
      INTEGER :: i_orb,j_orb
       ! LC
      !
      DO i_orb=1, n_orb
            !
            DO j_orb=1,n_orb
               !
               IF ((distance_LC(i_orb,j_orb) <= limit_1) .AND. (distance_LC(i_orb,j_orb) >= limit_0) ) THEN
                     CALL find_coeff((TRIM(orb(i_orb)%nature)), (TRIM(orb(j_orb)%nature)), 1, ene)
                     hamiltonian_LC(i_orb,j_orb) = ene
               ELSE IF ((distance_LC(i_orb,j_orb) <= limit_2) .AND. (distance_LC(i_orb,j_orb) >= limit_1).AND. use_second ) THEN
                     CALL find_coeff((TRIM(orb(i_orb)%nature)), (TRIM(orb(j_orb)%nature)), 2, ene)
                     hamiltonian_LC(i_orb,j_orb) = ene
               ELSE IF ((distance_LC(i_orb,j_orb) <= limit_3) .AND. (distance_LC(i_orb,j_orb) >= limit_2).AND. use_third )  THEN
                     CALL find_coeff((TRIM(orb(i_orb)%nature)), (TRIM(orb(j_orb)%nature)), 3, ene)
                     hamiltonian_LC(i_orb,j_orb) = ene
               ELSE
                     hamiltonian_LC(i_orb,j_orb) = CZERO
               ENDIF
               !
            ENDDO
            !
      ENDDO


     !CR
      DO i_orb=1, n_orb
            !
            DO j_orb=1,n_orb
               !
               IF ((distance_CR(i_orb,j_orb) <= limit_1) .AND. (distance_CR(i_orb,j_orb) >= limit_0) ) THEN
                     CALL find_coeff((TRIM(orb(i_orb)%nature)), (TRIM(orb(j_orb)%nature)), 1, ene)
                     hamiltonian_CR(i_orb,j_orb) = CONJG(ene)
               ELSE IF ((distance_CR(i_orb,j_orb) <= limit_2) .AND. (distance_CR(i_orb,j_orb) >= limit_1).AND. use_second ) THEN
                     CALL find_coeff((TRIM(orb(i_orb)%nature)), (TRIM(orb(j_orb)%nature)), 2, ene)
                     hamiltonian_CR(i_orb,j_orb) = CONJG(ene)
               ELSE IF ((distance_CR(i_orb,j_orb) <= limit_3) .AND. (distance_CR(i_orb,j_orb) >= limit_2).AND. use_third )  THEN
                     CALL find_coeff((TRIM(orb(i_orb)%nature)), (TRIM(orb(j_orb)%nature)), 3, ene)
                     hamiltonian_CR(i_orb,j_orb) = CONJG(ene)
               ELSE
                     hamiltonian_CR(i_orb,j_orb) = CZERO
               ENDIF
               !
            ENDDO
            !
      ENDDO



END SUBROUTINE init_hopping_hamiltonian



!*******************************************************************
   SUBROUTINE print_hamiltonian
   !*******************************************************************
      CHARACTER(17)      :: subname="print_hamiltonian"
      INTEGER :: ierr
      REAL(dbl) :: test
      CHARACTER(nstrx)   :: attr
      CHARACTER(nstrx)  :: name
      CHARACTER( LEN=nstrx )  :: filename
      REAL(dbl), ALLOCATABLE :: position(:,:)
      COMPLEX(dbl), ALLOCATABLE :: temp(:,:)
      INTEGER            :: nrtot_C
      INTEGER            :: nr_C(3)
      INTEGER            :: ivr_C(3,3)
      INTEGER            :: ik, ir, i_orb

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! initialization PART
         name='HAMILTONIAN'
   !
         nrtot_C = 3
         nr_C(1) = 3
         nr_C(2) = 1
         nr_C(3) = 1 
   !
         !ivr_C(3,nrtot)
         ivr_C(:,:)=0
         ivr_C(1,1)=-1
         ivr_C(1,2)=0
         ivr_C(1,3)=1

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ALLOCATE( position(3, n_orb),  STAT=ierr)
                IF ( ierr /=0 ) CALL errore(subname,'allocating position',ABS(ierr))
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         !
         DO i_orb=1, n_orb
            !
            position(:,i_orb) = orb(i_orb)%coord(:)
            !
         ENDDO
         !
         CALL file_open(ham_unit,TRIM(datafile_H),PATH="/",ACTION="write", FORM='formatted')
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CONDUCTOR PART
            !
            CALL iotk_write_begin(ham_unit,TRIM(name))
            CALL iotk_write_attr(attr,"dimwann",n_orb,FIRST=.TRUE.) 
            CALL iotk_write_attr(attr,"nkpts",nrtot_C) 
            CALL iotk_write_attr(attr,"nk",nr_C) 
            CALL iotk_write_attr(attr,"nrtot",nrtot_C) 
            CALL iotk_write_attr(attr,"nr",nr_C) 
            CALL iotk_write_empty(ham_unit,"DATA",ATTR=attr)
                  !
            CALL iotk_write_dat(ham_unit,"IVR", ivr_C, ATTR=attr, COLUMNS=3, IERR=ierr) 
                  IF (ierr/=0) CALL errore(subname,'writing ivr',ABS(ierr))

      !
            CALL iotk_write_begin(ham_unit,"RHAM")
            CALL iotk_write_dat(ham_unit,"VR"//TRIM(iotk_index(1)), hamiltonian_LC(:,:))
            CALL iotk_write_dat(ham_unit,"VR"//TRIM(iotk_index(2)), hamiltonian(:,:))
            CALL iotk_write_dat(ham_unit,"VR"//TRIM(iotk_index(3)), hamiltonian_CR(:,:))
            CALL iotk_write_end(ham_unit,"RHAM")
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            CALL iotk_write_dat(ham_unit,"WANCENTER", position, ATTR=attr, COLUMNS=3, IERR=ierr) 
                 IF (ierr/=0) CALL errore(subname,'writing position',ABS(ierr))
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            CALL iotk_write_end(ham_unit,TRIM(name))
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         CALL file_close(ham_unit,PATH="/",ACTION="write")
         WRITE( stdout,"(/,'  Hamiltonian of the Central region written on file : ',a)") TRIM(datafile_H)




         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         ALLOCATE( temp(n_orb , n_orb),  STAT=ierr)
                IF ( ierr /=0 ) CALL errore(subname,'allocating temp',ABS(ierr))
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         temp(:,:)= CONJG (TRANSPOSE (hamiltonian_LC(:,:)))

         CALL ioname('rechamiltonian',filename)


         CALL file_open(rec_ham_unit,TRIM(filename),PATH="/",ACTION="write", FORM='formatted')
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CONDUCTOR PART
            !
            CALL iotk_write_begin(rec_ham_unit,TRIM(name))
            CALL iotk_write_attr(attr,"dimwann",n_orb,FIRST=.TRUE.) 
            CALL iotk_write_attr(attr,"nkpts",nrtot_C) 
            CALL iotk_write_attr(attr,"nk",nr_C) 
            CALL iotk_write_attr(attr,"nrtot",nrtot_C) 
            CALL iotk_write_attr(attr,"nr",nr_C) 
            CALL iotk_write_empty(rec_ham_unit,"DATA",ATTR=attr)
                  !
            CALL iotk_write_dat(rec_ham_unit,"IVR", ivr_C, ATTR=attr, COLUMNS=3, IERR=ierr) 
                  IF (ierr/=0) CALL errore(subname,'writing ivr',ABS(ierr))

      !
            CALL iotk_write_begin(rec_ham_unit,"RHAM")
            CALL iotk_write_dat(rec_ham_unit,"VR"//TRIM(iotk_index(1)), temp(:,:))
            CALL iotk_write_dat(rec_ham_unit,"VR"//TRIM(iotk_index(2)), hamiltonian(:,:))
            CALL iotk_write_dat(rec_ham_unit,"VR"//TRIM(iotk_index(3)), hamiltonian_CR(:,:))
            CALL iotk_write_end(rec_ham_unit,"RHAM")
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
            CALL iotk_write_dat(rec_ham_unit,"WANCENTER", position, ATTR=attr, COLUMNS=3, IERR=ierr) 
                 IF (ierr/=0) CALL errore(subname,'writing position',ABS(ierr))

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            CALL iotk_write_end(rec_ham_unit,TRIM(name))

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

         CALL file_close(rec_ham_unit,PATH="/",ACTION="write")
         CALL ioname('rechamiltonian',filename,LPATH=.FALSE.)
         WRITE( stdout,"(/,'  Rec hamiltonian written on file : ',5x,a)") TRIM(filename)

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IF ( ALLOCATED( temp  ) ) THEN
           DEALLOCATE ( temp, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating temp', ABS(ierr) )
      ENDIF

         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    IF ( ALLOCATED( position  ) ) THEN
           DEALLOCATE ( position, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating position', ABS(ierr) )
      ENDIF
 
         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
END SUBROUTINE print_hamiltonian


END MODULE hamiltonian_module


