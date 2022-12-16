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
   USE dim_variable_module,        ONLY : n_orb, cell_size
   USE kinds
   USE constants,            ONLY : CZERO, CONE, ZERO, EPS_m11
   USE identity_module,      ONLY : orbitale
   USE orbital_module,       ONLY : orb
   USE io_module,            ONLY : met_unit => aux_unit, &
                                    ioname
   USE iotk_module
   USE files_module, ONLY : file_open, file_close
   USE io_global_module,     ONLY : stdin, stdout

   IMPLICIT NONE
   PRIVATE 
   SAVE

! Variables
    REAL(dbl), ALLOCATABLE :: distance(:,:)
!
    REAL(dbl), ALLOCATABLE :: distance_LC(:,:)
!
    REAL(dbl), ALLOCATABLE :: distance_CR(:,:)
!
    LOGICAL :: distance_alloc

! Status

!
   PUBLIC :: distance
!
   PUBLIC :: distance_CR
!
   PUBLIC :: distance_LC
!
   PUBLIC :: distance_alloc



! functions


  PUBLIC :: distance_allocate
  PUBLIC :: distance_deallocate
  PUBLIC :: init_metric
  PUBLIC :: init_metric_CR
  PUBLIC :: init_metric_LC
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

   ALLOCATE( distance(n_orb, n_orb),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating distance',ABS(ierr))
   !
   ALLOCATE( distance_CR(n_orb, n_orb),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating distance_CR',ABS(ierr))
   !
   ALLOCATE( distance_LC(n_orb, n_orb),  STAT=ierr)
       IF ( ierr /=0 ) CALL errore(subname,'allocating distance_LC',ABS(ierr))

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

     IF ( ALLOCATED( distance_CR  ) ) THEN
           DEALLOCATE ( distance_CR, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating distance_CR', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( distance_LC  ) ) THEN
           DEALLOCATE ( distance_LC, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating distance_LC', ABS(ierr) )
      ENDIF
 
     distance_alloc=.FALSE.

   END SUBROUTINE distance_deallocate


!*******************************************************************
   SUBROUTINE init_metric()
   !*******************************************************************
      CHARACTER(11)      :: subname="init_metric"
      INTEGER :: ierr
      INTEGER :: i_orb, j_orb
      INTEGER :: i_coord
      REAL(dbl) :: test


      DO i_orb=1, n_orb
            !
            DO j_orb=i_orb, n_orb
               distance(i_orb, j_orb) = ZERO
               !
               DO i_coord = 1 , 3
                     distance(i_orb, j_orb) =  distance(i_orb, j_orb) + (orb(i_orb)%coord(i_coord)- orb(j_orb)%coord(i_coord))**2
               ENDDO
               !
               distance(i_orb,j_orb) = SQRT(distance(i_orb,j_orb) )
               distance(j_orb,i_orb) = distance(i_orb,j_orb)
               !
            ENDDO
      ENDDO
    ! TESTS


      ! test diagonal
      test=ZERO
      DO i_orb=1, n_orb
          test= test + distance(i_orb,i_orb)
      ENDDO

      IF (( test <= - EPS_m11 ).OR. ( test >= EPS_m11 ) ) CALL errore(subname,'distance between i and i /= 0 ',1)

      ! test (i,j)=(j,i)
      test = ZERO
      DO i_orb=1, n_orb
            !
            DO j_orb=i_orb, n_orb
               !
               test =  test + ABS( distance(i_orb, j_orb) -  distance(j_orb, i_orb))
               !
            ENDDO
            !
      ENDDO

      IF (( test <= - EPS_m11 ).OR. ( test >= EPS_m11 ) ) CALL errore(subname,'distance btwn i,j /= j,i ',2)

END SUBROUTINE init_metric


!*******************************************************************
   SUBROUTINE init_metric_CR()
   !*******************************************************************
      CHARACTER(14)      :: subname="init_metric_CR"
      INTEGER :: ierr
      INTEGER :: i_orb, j_orb
      INTEGER :: i_coord
      REAL(dbl) :: test
      REAL(dbl) :: coord_temp(3)

      DO i_orb=1, n_orb
            !

            ! ATTENTION DIFFERENT DU CAS C-C j_orb commence à 1
            DO j_orb=1, n_orb
               distance_CR(i_orb, j_orb) = ZERO
               !
               coord_temp(1) = cell_size + orb(j_orb)%coord(1)
               coord_temp(2) = orb(j_orb)%coord(2)
               coord_temp(3) = orb(j_orb)%coord(3)

               ! ATTENTION I_ORB EST DANS C ET J_ORB EST DANS R
               DO i_coord = 1 , 3

                     distance_CR(i_orb, j_orb) =  distance_CR(i_orb, j_orb) + (orb(i_orb)%coord(i_coord)- coord_temp(i_coord))**2
               ENDDO
               !
               distance_CR(i_orb,j_orb) = SQRT(distance_CR(i_orb,j_orb) )

               !
            ENDDO
      ENDDO
    ! TESTS

END SUBROUTINE init_metric_CR

!*******************************************************************
   SUBROUTINE init_metric_LC()
   !*******************************************************************
      CHARACTER(14)      :: subname="init_metric_LC"
      INTEGER :: ierr
      INTEGER :: i_orb, j_orb
      INTEGER :: i_coord
      REAL(dbl) :: test
      REAL(dbl) :: coord_temp(3)


      DO i_orb=1, n_orb
            !
            coord_temp(1) = orb(i_orb)%coord(1) - cell_size
            coord_temp(2) = orb(i_orb)%coord(2)
            coord_temp(3) = orb(i_orb)%coord(3)

            ! ATTENTION DIFFERENT DU CAS C-C j_orb commence à 1
            DO j_orb=1, n_orb
               distance_LC(i_orb, j_orb) = ZERO
               !

               ! ATTENTION I_ORB EST DANS L ET J_ORB EST DANS C
               DO i_coord = 1 , 3
                     distance_LC(i_orb, j_orb) =  distance_LC(i_orb, j_orb) + ( coord_temp(i_coord) - orb(j_orb)%coord(i_coord))**2
               ENDDO
               !
               distance_LC(i_orb,j_orb) = SQRT(distance_LC(i_orb,j_orb) )

               !
            ENDDO
      ENDDO
    ! TESTS

END SUBROUTINE init_metric_LC

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
       CALL iotk_write_dat(met_unit,"DISTANCE_LC", distance_LC(:,:))
       CALL iotk_write_dat(met_unit,"DISTANCE_CR", distance_CR(:,:))
       CALL iotk_write_end(met_unit,"DIST_MAT")
         !

       CALL iotk_write_end(met_unit,TRIM(name))

      CALL file_close(met_unit,PATH="/",ACTION="write")

      CALL ioname('distance',filename,LPATH=.FALSE.)
      WRITE( stdout,"(/,'  Distance written on file : ',5x,a)") TRIM(filename)


END SUBROUTINE print_distance



END MODULE distance_module


