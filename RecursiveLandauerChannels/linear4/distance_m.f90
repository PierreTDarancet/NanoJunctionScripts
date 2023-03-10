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
   MODULE distance_module
!*********************************************
   USE parameters,      ONLY : nstrx
   USE dim_variable_module,   ONLY : dim_recursion, nb_max_first,  nb_max_second, nb_max_third, limit_0, &
                                     limit_1, limit_2, limit_3, n_orb
   USE control_variable_module,    ONLY : print_d => print_distance,  use_second,   use_third
   USE kinds
   USE constants,            ONLY : CZERO, CONE, ZERO, EPS_m11
   USE identity_module,      ONLY : orbitale
   USE orbital_module,       ONLY : orb_recursion
   USE io_module,            ONLY : met_unit => aux_unit, &
                                    ioname
   USE iotk_module
   USE files_module,         ONLY : file_open, file_close
   USE io_global_module,     ONLY : stdin, stdout
   USE coeff_module,         ONLY : find_coeff, find_coeff_onsite

   IMPLICIT NONE
   PRIVATE 
   SAVE

! Variables
    REAL(dbl), ALLOCATABLE :: distance(:,:)
!
    INTEGER, ALLOCATABLE :: nb_first(:)
!
    INTEGER, ALLOCATABLE :: nb_second(:)
!
    INTEGER, ALLOCATABLE :: nb_third(:)
!
    INTEGER, ALLOCATABLE :: id_first(:,:)
!
    INTEGER, ALLOCATABLE :: id_second(:,:)
!
    INTEGER, ALLOCATABLE :: id_third(:,:)
!
    COMPLEX(dbl), ALLOCATABLE :: ene_first(:,:)
!
    COMPLEX(dbl), ALLOCATABLE :: ene_second(:,:)
!
    COMPLEX(dbl), ALLOCATABLE :: ene_third(:,:)
!
    COMPLEX(dbl), ALLOCATABLE :: ene_onsite(:)
!
    LOGICAL :: distance_alloc

! Status

!
   PUBLIC :: distance
!
   PUBLIC :: nb_first
!
   PUBLIC :: nb_second
!
   PUBLIC :: nb_third
!
   PUBLIC :: id_first
!
   PUBLIC :: id_second
!
   PUBLIC :: id_third
!
   PUBLIC :: ene_first
!
   PUBLIC :: ene_second
!
   PUBLIC :: ene_third
!
   PUBLIC :: ene_onsite
!

!
   PUBLIC :: distance_alloc



! functions


  PUBLIC :: distance_allocate
  PUBLIC :: distance_deallocate
  PUBLIC :: init_metric
  PUBLIC :: print_distance



CONTAINS


!*******************************************************************
   SUBROUTINE distance_allocate
   !*******************************************************************
      CHARACTER(17)      :: subname="distance_allocate"
      INTEGER :: ierr

   !
   ! allocate basic quantities
   !

   IF  (distance_alloc) CALL errore(subname,'distance already allocated',1)


   IF (print_d) THEN
      ALLOCATE( distance(dim_recursion, dim_recursion),  STAT=ierr)
         IF ( ierr /=0 ) CALL errore(subname,'allocating distance',ABS(ierr))
   ENDIF
   !

   ALLOCATE(  nb_first(dim_recursion),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating nb_first',ABS(ierr))
   !
   ALLOCATE( nb_second(dim_recursion),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating nb_second',ABS(ierr))
   !
   ALLOCATE( nb_third(dim_recursion),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating nb_third',ABS(ierr))
   !
   ALLOCATE( id_first(nb_max_first, dim_recursion),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating id_first',ABS(ierr))
   !
   ALLOCATE( id_second(nb_max_second, dim_recursion),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating id_second',ABS(ierr))
   !
   ALLOCATE( id_third(nb_max_third, dim_recursion),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating id_third',ABS(ierr))
   !
   ALLOCATE( ene_first(nb_max_first, dim_recursion),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating ene_first',ABS(ierr))
   !
   ALLOCATE( ene_second(nb_max_second, dim_recursion),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating ene_second',ABS(ierr))
   !
   ALLOCATE( ene_third(nb_max_third, dim_recursion),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating ene_third',ABS(ierr))
   !
   ALLOCATE( ene_onsite(dim_recursion),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating ene_onsite',ABS(ierr))
   !

    distance_alloc=.TRUE.



END SUBROUTINE distance_allocate



!**********************************************************
   SUBROUTINE distance_deallocate()
   !**********************************************************
   IMPLICIT NONE
      CHARACTER(19)      :: subname="distance_deallocate"
      INTEGER :: ierr

    IF (.NOT. distance_alloc) CALL errore(subname,'distance not allocated',1)

     IF ( ALLOCATED( distance  ) ) THEN
           DEALLOCATE ( distance, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating distance', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( nb_first ) ) THEN
           DEALLOCATE ( nb_first, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating nb_first', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( nb_second ) ) THEN
           DEALLOCATE ( nb_second, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating nb_second', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( nb_third ) ) THEN
           DEALLOCATE ( nb_third , STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating nb_third', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( id_first ) ) THEN
           DEALLOCATE ( id_first , STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating id_first', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( id_second ) ) THEN
           DEALLOCATE ( id_second , STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating id_second', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( id_third ) ) THEN
           DEALLOCATE ( id_third, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating id_third', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( ene_first ) ) THEN
           DEALLOCATE ( ene_first , STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating ene_first', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( ene_second ) ) THEN
           DEALLOCATE ( ene_second , STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating ene_second', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( ene_third ) ) THEN
           DEALLOCATE ( ene_third, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating ene_third', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( ene_onsite ) ) THEN
           DEALLOCATE ( ene_onsite, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating ene_onsite', ABS(ierr) )
      ENDIF
     !
     distance_alloc=.FALSE.

   END SUBROUTINE distance_deallocate


!*******************************************************************
   SUBROUTINE init_metric()
   !*******************************************************************
      CHARACTER(11)      :: subname="init_metric"
      INTEGER :: ierr
      INTEGER :: i_rec, j_rec
      INTEGER :: i_coord
      REAL(dbl) :: test
      REAL(dbl) :: distance_aux
      COMPLEX(dbl) :: ene


    IF (.NOT. distance_alloc) CALL errore(subname,'distance not allocated',1)
     !init
     nb_first(:)     =  0
     nb_second(:)    =  0
     nb_third(:)     =  0
     !
     id_first(:,:)   =  0
     id_second(:,:)  =  0
     id_third(:,:)   =  0
     !
     ene_first(:,:)  =  CZERO
     ene_second(:,:) =  CZERO
     ene_third(:,:)  =  CZERO
     !

     IF (print_d) THEN
      distance(:,:) = ZERO
     ENDIF

      !
      ene = CZERO
      !
      DO i_rec=1, dim_recursion
            !
            CALL find_coeff_onsite((TRIM(orb_recursion(i_rec)%nature)), ene)
            !
            ene_onsite(i_rec) = ene
            !
      ENDDO
      !

      !
      DO i_rec=1, dim_recursion-1
            !
            DO j_rec=i_rec+1, dim_recursion
               !
               distance_aux = ZERO
               DO i_coord = 1 , 3
                     !
                     distance_aux =  distance_aux + (orb_recursion(i_rec)%coord(i_coord) - orb_recursion(j_rec)%coord(i_coord))**2
                     !
               ENDDO
               !
               distance_aux = SQRT(distance_aux )
               !
               IF ( print_d ) THEN
                  distance(i_rec,j_rec) = distance_aux
                  distance(j_rec,i_rec) = distance(i_rec,j_rec)
               ENDIF
               !
               IF ((distance_aux <= limit_1) .AND. (distance_aux > limit_0) ) THEN
                   !
                   !
                   nb_first(i_rec) = nb_first(i_rec) + 1
                   nb_first(j_rec) = nb_first(j_rec) + 1
                   !
                   IF (nb_first(i_rec) > nb_max_first)  CALL errore(subname,' nb first > nb max ',i_rec)
                   IF (nb_first(j_rec) > nb_max_first)  CALL errore(subname,' nb first > nb max ',j_rec)
                   !
                   id_first(nb_first(i_rec),i_rec) = j_rec 
                   id_first(nb_first(j_rec),j_rec) = i_rec 
                   !
                   !
                   CALL find_coeff((TRIM(orb_recursion(i_rec)%nature)), (TRIM(orb_recursion(j_rec)%nature)), 1, ene)
                   !
                   !PRINT*, ene
                   ene_first(nb_first(i_rec),i_rec) = ene
                   ene_first(nb_first(j_rec),j_rec) = CONJG(ene)
                   !
               ELSE IF ((distance_aux <= limit_2) .AND. (distance_aux > limit_1) .AND. use_second ) THEN
                   !
                   !
                   nb_second(i_rec) = nb_second(i_rec) + 1
                   nb_second(j_rec) = nb_second(j_rec) + 1
                   !
                   IF (nb_second(i_rec) > nb_max_second)  CALL errore(subname,' nb second > nb max ',i_rec)
                   IF (nb_second(j_rec) > nb_max_second)  CALL errore(subname,' nb second > nb max ',j_rec)
                   !
                   id_second(nb_second(i_rec),i_rec) = j_rec 
                   id_second(nb_second(j_rec),j_rec) = i_rec 
                   !
                   !
                   CALL find_coeff((TRIM(orb_recursion(i_rec)%nature)), (TRIM(orb_recursion(j_rec)%nature)), 2, ene)
                   !
                   ene_second(nb_second(i_rec),i_rec) = ene
                   ene_second(nb_second(j_rec),j_rec) = CONJG(ene)
                   !
               ELSE IF ((distance_aux <= limit_3) .AND. (distance_aux > limit_2).AND. use_third )  THEN
                   !
                   !
                   nb_third(i_rec) = nb_third(i_rec) + 1
                   nb_third(j_rec) = nb_third(j_rec) + 1
                   !
                   IF (nb_third(i_rec) > nb_max_third)  CALL errore(subname,' nb third > nb max ',i_rec)
                   IF (nb_third(j_rec) > nb_max_third)  CALL errore(subname,' nb third > nb max ',j_rec)
                   !
                   id_third(nb_third(i_rec),i_rec) = j_rec 
                   id_third(nb_third(j_rec),j_rec) = i_rec 
                   !
                   !
                   CALL find_coeff((TRIM(orb_recursion(i_rec)%nature)), (TRIM(orb_recursion(j_rec)%nature)), 3, ene)
                   !a
                   ene_third(nb_third(i_rec),i_rec) = ene
                   ene_third(nb_third(j_rec),j_rec) = CONJG(ene)
                   !
               ELSE
                   !
                   !
                   IF ( i_rec == j_rec )  CALL errore(subname,' i_rec and j_rec are supposed to be different ',i_rec)
                   !
                   !
               ENDIF
               !

               !
            ENDDO
      ENDDO
    ! TESTS
!       PRINT*, '-----------------------------------------------------------------'
!       PRINT*, '                        nb_first'
!       PRINT*, '-----------------------------------------------------------------'
!       PRINT*, nb_first(:)
!       PRINT*, '-----------------------------------------------------------------'
!       PRINT*, '-----------------------------------------------------------------'
!       PRINT*, '-----------------------------------------------------------------'
!       DO i_rec=1, dim_recursion
!          PRINT*, '                        ID_first'
!          PRINT*, ' i_rec'
!          PRINT*, i_rec
!          PRINT*, ' id'
!          PRINT*, id_first(:,i_rec)
!          PRINT*, '-----------------------------------------------------------------'
!       ENDDO


END SUBROUTINE init_metric


!*******************************************************************
   SUBROUTINE print_distance
   !*******************************************************************
       CHARACTER(14)      :: subname="print_distance"
       CHARACTER(8) :: name="DISTANCE"
       CHARACTER(nstrx)   :: attr
       INTEGER            :: iter
       INTEGER            :: ierr, i_sub
      CHARACTER( LEN=nstrx )  :: filename
      INTEGER            :: nrtot_C
      INTEGER            :: nr_C(3)
      INTEGER            :: ivr_C(3,3)



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
 
       CALL ioname('distance',filename)

       CALL file_open(met_unit,TRIM(filename),PATH="/",ACTION="write", &
                              FORM='formatted')

         !
       CALL iotk_write_begin(met_unit,TRIM(name))
       CALL iotk_write_attr(attr,"dimwann",n_orb,FIRST=.TRUE.) 
       CALL iotk_write_attr(attr,"nkpts",nrtot_C) 
       CALL iotk_write_attr(attr,"nk",nr_C) 
       CALL iotk_write_attr(attr,"nrtot",nrtot_C) 
       CALL iotk_write_attr(attr,"nr",nr_C) 
       CALL iotk_write_empty(met_unit,"DATA",ATTR=attr)
        !
       CALL iotk_write_dat(met_unit,"IVR", ivr_C, ATTR=attr, COLUMNS=3, IERR=ierr) 
               IF (ierr/=0) CALL errore(subname,'writing ivr',ABS(ierr))
       !
       CALL iotk_write_begin(met_unit,"DIST_MAT")
       CALL iotk_write_dat(met_unit,"DISTANCE", distance(:,:))
       CALL iotk_write_end(met_unit,"DIST_MAT")
         !

       CALL iotk_write_end(met_unit,TRIM(name))

      CALL file_close(met_unit,PATH="/",ACTION="write")

      CALL ioname('distance',filename,LPATH=.FALSE.)
      WRITE( stdout,"(/,'  Distance written on file : ',5x,a)") TRIM(filename)


END SUBROUTINE print_distance



END MODULE distance_module


