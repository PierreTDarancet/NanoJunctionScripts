!
! Copyright (C) 2006 LEPES-CNRS Grenoble
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!*********************************************
   MODULE recursion_module
!*********************************************
   USE parameters,      ONLY : nstrx
   USE dim_variable_module,        ONLY : dimwan, dim_subspace, dim_recursion
   USE kinds
   USE constants,            ONLY : ZERO, CZERO, CONE, ONE
   USE io_global_module,     ONLY : stdin, stdout
   USE io_module,     ONLY :  ham_unit, mat_unit => aux_unit, &
                                    sta_unit => aux2_unit, cen_unit => aux3_unit, ioname
   USE files_module, ONLY : file_open, file_close
   USE parser_module
   USE iotk_module
   USE control_variable_module,      ONLY : max_iter, conv_criteria
   USE recursion_function_module, ONLY :  calculate_A_matrix, calculate_recursion_state, &
                                          normalize_recursion_state, calculate_variation


   IMPLICIT NONE
   PRIVATE 
   SAVE

   ! variables for recursion
            ! For first iteration
   COMPLEX(dbl), ALLOCATABLE ::  A_matrix_0(:,:), B_matrix_0(:,:), recursion_state_0(:,:)

   COMPLEX(dbl), ALLOCATABLE ::  A_matrix(:,:,:), B_matrix(:,:,:), recursion_state(:,:,:)

   COMPLEX(dbl), ALLOCATABLE ::  final_hamiltonian(:,:)

   COMPLEX(dbl), ALLOCATABLE ::  recursion_hamiltonian(:,:)

   COMPLEX(dbl), ALLOCATABLE ::  total_hamiltonian(:,:)

   COMPLEX(dbl), ALLOCATABLE ::  linear_hamiltonian(:,:)

   REAL(dbl), ALLOCATABLE :: wannier_coordinates(:,:)
   REAL(dbl), ALLOCATABLE :: wannier_coordinates_sub_recursion(:,:)
   REAL(dbl), ALLOCATABLE :: norm(:)


   INTEGER, ALLOCATABLE   :: wan_num(:)
   INTEGER, ALLOCATABLE   :: subrecur2tot(:)
   INTEGER  :: final_iter

   LOGICAL :: alloc


   PUBLIC :: A_matrix_0, B_matrix_0, recursion_state_0
   PUBLIC :: A_matrix, B_matrix, recursion_state
   PUBLIC :: final_hamiltonian
   PUBLIC :: recursion_hamiltonian, total_hamiltonian
   PUBLIC :: wannier_coordinates
   PUBLIC :: wan_num
   PUBLIC :: wannier_coordinates_sub_recursion
   PUBLIC :: subrecur2tot
   PUBLIC :: alloc
   PUBLIC :: final_iter

   ! general functions
   PUBLIC :: recursion_allocate
   PUBLIC :: recursion_deallocate
   PUBLIC :: initial_value
   PUBLIC :: first_iteration
   PUBLIC :: recursion_loop
   PUBLIC :: init_wan_num
   PUBLIC :: build_ham_and_coor
   PUBLIC :: select_wf
   PUBLIC :: reduce_hamiltonian
   PUBLIC :: init_first_recursion_state
   PUBLIC :: print_hamiltonian
   PUBLIC :: print_matrix
   PUBLIC :: print_states
   PUBLIC :: print_wan_center
CONTAINS



!*******************************************************************
   SUBROUTINE recursion_allocate
   !*******************************************************************
      CHARACTER(18)      :: subname="recursion_allocate"
      INTEGER :: ierr


    ! Allocate hamiltonian
   ALLOCATE ( total_hamiltonian(dimwan,dimwan), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating total_hamiltonian', ABS(ierr) )
   ALLOCATE ( recursion_hamiltonian(dim_recursion,dim_recursion), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating recursion_hamiltonian', ABS(ierr) )
   ALLOCATE ( final_hamiltonian(((max_iter+2)*dim_subspace),((max_iter+2)*dim_subspace)), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating final_hamiltonian', ABS(ierr) )


    ! Allocate A-type matrices
   ALLOCATE ( A_matrix_0(dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating A_matrix_0', ABS(ierr) )
   ALLOCATE ( A_matrix(max_iter,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating A_matrix', ABS(ierr) )
    ! Allocate B-type matrices
   ALLOCATE ( B_matrix_0(dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix_0', ABS(ierr) )
   ALLOCATE ( B_matrix(max_iter,dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating B_matrix', ABS(ierr) )

    ! Allocate Wave functions
   ALLOCATE ( recursion_state_0(dim_subspace, dim_recursion), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating recursion_state_0', ABS(ierr) )
   ALLOCATE ( recursion_state(max_iter+2,dim_subspace, dim_recursion), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating recursion_state', ABS(ierr) )

    ! Allocate Wannier Centers coordinates
   ALLOCATE ( wannier_coordinates(3,dimwan), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating wannier_coordinates', ABS(ierr) )

    ! Allocate Wan_num i.e. WF's num of the first  recursion states or subspace if dim_subspace > 1
   ALLOCATE ( wan_num(dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating wan_num', ABS(ierr) )


    ! libs to allow subspace <=> total
   ALLOCATE ( subrecur2tot(dim_recursion), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating subrecur2tot', ABS(ierr) )
   ALLOCATE ( wannier_coordinates_sub_recursion(3,dim_recursion), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating wannier_coordinates_sub_recursion', ABS(ierr) )


    ! allocate norm table
   ALLOCATE ( norm(max_iter+2), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating norm', ABS(ierr) )


   ! set alloc to true
    alloc=.TRUE.

END SUBROUTINE recursion_allocate

!**********************************************************
   SUBROUTINE recursion_deallocate()
   !**********************************************************
   IMPLICIT NONE
      CHARACTER(20)      :: subname="recursion_deallocate"
      INTEGER :: ierr

     IF ( ALLOCATED( total_hamiltonian  ) ) THEN
           DEALLOCATE ( total_hamiltonian, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating total_hamiltonian', ABS(ierr) )
      ENDIF

!   PRINT*, 'deallocating total_hamiltonian'

     IF ( ALLOCATED( final_hamiltonian  ) ) THEN
           DEALLOCATE (final_hamiltonian, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating final_hamiltonian', ABS(ierr) )
      ENDIF
!   PRINT*, 'deallocating fin_hamiltonian'

     IF ( ALLOCATED( A_matrix_0  ) ) THEN
           DEALLOCATE ( A_matrix_0, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating A_matrix_0', ABS(ierr) )
      ENDIF
!   PRINT*, 'deallocating A 0'

     IF ( ALLOCATED( A_matrix  ) ) THEN
           DEALLOCATE ( A_matrix , STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating A_matrix', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix_0  ) ) THEN
           DEALLOCATE ( B_matrix_0, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix_0', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( B_matrix  ) ) THEN
           DEALLOCATE ( B_matrix , STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating B_matrix', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( recursion_state_0  ) ) THEN
           DEALLOCATE ( recursion_state_0 , STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating recursion_state_0', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( recursion_state  ) ) THEN
           DEALLOCATE ( recursion_state, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating recursion_state', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( wannier_coordinates  ) ) THEN
           DEALLOCATE ( wannier_coordinates, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating wannier_coordinates', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( wan_num  ) ) THEN
           DEALLOCATE ( wan_num, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating wan_num', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( subrecur2tot  ) ) THEN
           DEALLOCATE ( subrecur2tot, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating subrecur2tot', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( wannier_coordinates_sub_recursion ) ) THEN
           DEALLOCATE ( wannier_coordinates_sub_recursion, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating wannier_coordinates_sub_recursion', ABS(ierr) )
      ENDIF

     IF ( ALLOCATED( norm ) ) THEN
           DEALLOCATE ( norm, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating norm', ABS(ierr) )
      ENDIF


     IF ( ALLOCATED( recursion_hamiltonian  ) ) THEN
!            PRINT*, 'dans rec hamiltonian'
           DEALLOCATE (recursion_hamiltonian, STAT=ierr )
!            PRINT*, 'dans rec hamiltonian 2'

           IF( ierr /=0 ) CALL errore(subname, 'deallocating recursion_hamiltonian', ABS(ierr) )
      ENDIF

!    PRINT*, 'deallocating rec_hamiltonian'


   END SUBROUTINE recursion_deallocate

!**********************************************************
   SUBROUTINE initial_value()
   !**********************************************************
   IMPLICIT NONE
      CHARACTER(13)      :: subname="initial_value"
      INTEGER :: ierr

    ! matrix
     A_matrix_0(:,:)=CZERO
     A_matrix(:,:,:)=CZERO
     B_matrix_0(:,:)=CZERO
     B_matrix(:,:,:)=CZERO


    ! Wave function
    recursion_state_0(:,:)=CZERO
    recursion_state(:,:,:)=CZERO

    ! Hamiltonian
    total_hamiltonian(:,:)=CZERO
    recursion_hamiltonian(:,:)=CZERO
    final_hamiltonian(:,:)=CZERO

    ! Wannier coordinates
    wannier_coordinates(:,:)=ZERO

    ! WF indice for the first subspace
    wan_num(:)=0
    final_iter=0
    subrecur2tot(:)=0
    wannier_coordinates_sub_recursion(:,:)=ZERO
    !
    ! set norm
    norm(:)=ZERO

   END SUBROUTINE initial_value

!**********************************************************
   SUBROUTINE init_wan_num()
   !**********************************************************
   IMPLICIT NONE
      CHARACTER(12)      :: subname="init_wan_num"
      CHARACTER(nstrx)    :: wf_ind
       CHARACTER(nstrx)   :: attr
      INTEGER :: ierr, n_wf_ind, i, irecur, num
      INTEGER, ALLOCATABLE   :: wan_num_aux(:)
      LOGICAL :: found


      CALL iotk_scan_begin( stdin, 'SUBSPACE_WF', IERR=ierr )
         IF (ierr/=0) CALL errore(subname,'SUBSPACE_WF',ABS(ierr))
      CALL iotk_scan_empty(stdin, 'WF', ATTR=attr, IERR=ierr)
         IF (ierr/=0) CALL errore(subname, 'searching for WF', ABS(ierr) )
      !
      CALL iotk_scan_attr(attr, 'nb', wf_ind, FOUND=found, IERR=ierr)
         IF (ierr/=0) CALL errore(subname, 'searching for nb', ABS(ierr) )
         IF( .NOT. found ) wf_ind = '1'

      !
      CALL parser_replica( wf_ind, n_wf_ind, IERR=ierr)
      IF ( ierr/=0 ) CALL errore(subname,'wrong FMT in wf_ind string I',ABS(ierr))

      IF ( n_wf_ind /= dim_subspace ) CALL errore(subname,'invalid number of WF for subspace',3)

      ALLOCATE ( wan_num_aux(dim_subspace), STAT=ierr )
         IF( ierr /=0 ) CALL errore(subname, 'allocating wan_num_aux', ABS(ierr) )

      CALL parser_replica( wf_ind, n_wf_ind, wan_num_aux, IERR=ierr)

      IF ( ierr/=0 ) CALL errore(subname,'wrong FMT in nb string II',ABS(ierr))


      CALL iotk_scan_end( stdin, 'SUBSPACE_WF', IERR=ierr )
         IF (ierr/=0) CALL errore(subname,'searching end for SUBSPACE_WF',ABS(ierr))

      DO i=1,dim_subspace
          wan_num(i)=wan_num_aux(i)
      ENDDO
     !       PRINT*, 'wan_num'
     !       PRINT*, wan_num(:)
     !       PRINT*, 'subrecur2tot(:)'
     !
     !       PRINT*, subrecur2tot(:)
        ! check if initial state is included in the recursion subspace
      DO i=1,dim_subspace
          num = wan_num(i)
          found = .FALSE.
              DO irecur=1,dim_recursion
                  IF (num == subrecur2tot(irecur)) THEN
                    found=.TRUE.
                  ENDIF
             ENDDO
         IF (.NOT. found) CALL errore(subname,'initial subspace is not included in the recursion subspace',ABS(num))
      ENDDO



END SUBROUTINE init_wan_num

!*******************************************************************
   SUBROUTINE build_ham_and_coor(aux, wanaux)
   !*******************************************************************
     !input variables
      REAL(dbl), INTENT(in) :: wanaux(3,dimwan)
      COMPLEX(dbl), INTENT(in) :: aux(dimwan,dimwan)

     ! local variables
      CHARACTER(18)      :: subname="build_ham_and_coor"
      INTEGER :: ierr, iwan,jwan

   DO iwan=1,dimwan
       DO jwan=1,dimwan
          total_hamiltonian(jwan, iwan)= aux(jwan, iwan)

           !rajouter  test hermitainite
           !pour hamiltonien

       ENDDO
   wannier_coordinates(:,iwan)=wanaux(:,iwan)
   ENDDO

   !    DO iwan=1,dimwan
   !        DO jwan=1,dimwan
   !           IF (total_hamiltonian(jwan, iwan) /= CONJG (aux(iwan, jwan))) &
   !                  CALL errore(subname,'problem with hermitianity of the total hamiltonian',1)
   ! 
   !        ENDDO
   !    ENDDO



END SUBROUTINE build_ham_and_coor
!*******************************************************************
   SUBROUTINE select_wf(limit,lesser,greater)
   !*******************************************************************

      REAL(dbl), INTENT(in) :: limit
      LOGICAL, INTENT(in) ::lesser, greater
      CHARACTER(9)      :: subname="select_wf"
      INTEGER :: ierr, icount, iwan
       !should be generalized to non x-only transport, with transport_dir input
      INTEGER :: transport_dir=1

     ! initialize subrecur2tot
   IF (lesser) THEN
      icount=0

      DO  iwan=1,dimwan
         IF (wannier_coordinates(transport_dir,iwan) <= limit) THEN
            subrecur2tot(icount+1)=iwan
            icount=icount+1
         ENDIF
      ENDDO
      IF (icount /= dim_recursion) CALL errore(subname,'problem of dimensionnality in select_wf',ABS(icount))

   ELSEIF (greater) THEN
      icount=0

      DO  iwan=1,dimwan
         IF (wannier_coordinates(transport_dir,iwan) >= limit) THEN
            subrecur2tot(icount+1)=iwan
            icount=icount+1
         ENDIF
      ENDDO
      IF (icount /= dim_recursion) CALL errore(subname,'problem of dimensionnality in select_wf',ABS(icount))

   ELSE
     CALL errore(subname,'wrong boolean value of greater and lesser in select_wf',10)
   ENDIF


     ! initialize the coordinates of the WF center kept in the recursion
   DO  icount=1,dim_recursion
      wannier_coordinates_sub_recursion(:,icount) = wannier_coordinates(:,subrecur2tot(icount))
   ENDDO



END SUBROUTINE select_wf
!*******************************************************************
   SUBROUTINE reduce_hamiltonian()
   !*******************************************************************

      CHARACTER(18)      :: subname="reduce_hamiltonian"
      INTEGER :: ierr, irecur,jrecur


   DO  jrecur=1,dim_recursion
      DO  irecur=1,dim_recursion
         recursion_hamiltonian(irecur,jrecur) = total_hamiltonian(subrecur2tot(irecur),subrecur2tot(jrecur))
      ENDDO
   ENDDO

   ! test hermitianity to be added


END SUBROUTINE reduce_hamiltonian

!*******************************************************************
   SUBROUTINE init_first_recursion_state
   !*******************************************************************
      CHARACTER(26)      :: subname="init_first_recursion_state"
      INTEGER :: ierr, isub, irecur
 
     !  prepare subspace

    DO isub = 1, dim_subspace
        DO irecur =1, dim_recursion
           IF  ( wan_num(isub) == subrecur2tot(irecur) ) THEN
               recursion_state_0(isub,irecur) = ONE
           ELSE 
               recursion_state_0(isub,irecur) = ZERO
           ENDIF
        ENDDO
    ENDDO

END SUBROUTINE init_first_recursion_state

 !*******************************************************************
    SUBROUTINE first_iteration
    !*******************************************************************
    USE summary_module, ONLY : summary_recursion_input
    CHARACTER(15) :: subname="first_iteration"
    COMPLEX(dbl), ALLOCATABLE :: aux1(:,:), aux2(:,:)
    INTEGER       :: ierr, i_sub, j_sub
    INTEGER       :: i_sub_shift0, j_sub_shift0, i_sub_shift1, j_sub_shift1


!
!--------------------------------------------
! ... Startup
!--------------------------------------------
!
   ALLOCATE ( aux1(dim_subspace,dim_subspace), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating aux1', ABS(ierr) )
   ALLOCATE ( aux2(dim_subspace,dim_recursion), STAT=ierr )
        IF( ierr /=0 ) CALL errore(subname, 'allocating aux2', ABS(ierr) )


  aux1(:,:)=CZERO
  aux2(:,:)=CZERO

    !PRINT*, 'Debut, premi?re iteration'

     !
     CALL calculate_A_matrix( recursion_hamiltonian(:,:), recursion_state_0(:,:), A_matrix_0(:,:) ) 
     !

     !
     CALL calculate_recursion_state( recursion_hamiltonian(:,:), recursion_state_0(:,:), A_matrix_0(:,:), &
                                aux1(:,:), aux2(:,:), recursion_state(1,:,:) )
     !
    !PRINT*, 'Debut, normalize_recursion state 1'
     !
     CALL normalize_recursion_state( recursion_state(1,:,:), B_matrix_0(:,:), norm(1) )
     !
    !PRINT*, 'Fin, normalize_recursion state 1'
     DO j_sub=1, dim_subspace
         j_sub_shift0 = j_sub + dim_subspace
      DO i_sub=1, dim_subspace
         i_sub_shift0 = i_sub + dim_subspace
         final_hamiltonian(i_sub, j_sub)= A_matrix_0(i_sub, j_sub)
         final_hamiltonian(i_sub, j_sub_shift0)= B_matrix_0(i_sub, j_sub)
         final_hamiltonian(i_sub_shift0, j_sub)= CONJG (B_matrix_0(j_sub, i_sub))
      ENDDO
     ENDDO
     !
     CALL calculate_A_matrix( recursion_hamiltonian(:,:) , recursion_state(1,:,:) ,  A_matrix(1,:,:) ) 
     !

     !
     CALL calculate_recursion_state( recursion_hamiltonian(:,:), recursion_state(1,:,:), A_matrix(1,:,:), &
                         B_matrix_0(:,:), recursion_state_0(:,:), recursion_state(2,:,:) )
     !

     !
     CALL normalize_recursion_state(recursion_state(2,:,:),  B_matrix(1,:,:), norm(2) )
     !
     DO j_sub=1, dim_subspace
         j_sub_shift0 = j_sub + dim_subspace
         j_sub_shift1 = j_sub + 2*dim_subspace
      DO i_sub=1, dim_subspace
         i_sub_shift0 = i_sub + dim_subspace
         i_sub_shift1 = i_sub + 2*dim_subspace
         final_hamiltonian(i_sub_shift0, j_sub_shift0)= A_matrix(1, i_sub, j_sub)
         final_hamiltonian(i_sub_shift0, j_sub_shift1)= B_matrix(1,i_sub, j_sub)
         final_hamiltonian(i_sub_shift1, j_sub_shift0)= CONJG (B_matrix(1,j_sub, i_sub))
      ENDDO
     ENDDO
     CALL summary_recursion_input(stdout, A_matrix_0, B_matrix_0, recursion_state_0,recursion_hamiltonian, total_hamiltonian, wannier_coordinates, &
                                    wan_num, wannier_coordinates_sub_recursion, subrecur2tot, norm(1),  alloc)

     IF ( ALLOCATED( aux1  ) ) THEN
           DEALLOCATE ( aux1, STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating aux1', ABS(ierr) )
      ENDIF
     IF ( ALLOCATED( aux2  ) ) THEN
           DEALLOCATE ( aux2 , STAT=ierr )
           IF( ierr /=0 ) CALL errore(subname, 'deallocating aux2', ABS(ierr) )
      ENDIF




 END SUBROUTINE first_iteration
 !*******************************************************************
    SUBROUTINE recursion_loop
    !*******************************************************************
   !
    USE summary_module, ONLY : summary_recursion

    CHARACTER(14) :: subname="recursion_loop"
    REAL(dbl) :: variation(dim_subspace), var
    INTEGER       :: iter, ierr, i_sub, j_sub
    INTEGER       :: i_sub_shift0, j_sub_shift0, i_sub_shift1, j_sub_shift1
    INTEGER       :: iteration

!
!--------------------------------------------
! ... Startup
!--------------------------------------------
!

  iteration=1
  CALL summary_recursion(stdout, iteration, A_matrix(iteration,:,:), B_matrix(iteration,:,:), variation, var, norm(iteration+1))


  iteration_loop : &
  DO iter = 1, max_iter-1
     ! real iteration indice
     iteration = iter+1
     !
     CALL calculate_A_matrix(  recursion_hamiltonian(:,:), recursion_state(iteration,:,:) , A_matrix(iteration,:,:) ) 
     !
     !
     CALL calculate_recursion_state( recursion_hamiltonian(:,:), recursion_state(iteration,:,:), A_matrix(iteration,:,:), &
                         B_matrix(iteration-1,:,:), recursion_state(iteration-1,:,:), recursion_state(iteration+1,:,:) )
     !

     !
     CALL normalize_recursion_state(  recursion_state(iteration+1,:,:), B_matrix(iteration,:,:), norm(iteration+1) )
     !

      !ajouter fonction add data to recursion_hamiltonian
      ! Call add data to recursion_hamiltonian(recursion_hamiltonian, A_matrix(iter+1,:,:) , B_matrix(iter+1), iter)
     DO j_sub=1, dim_subspace
         j_sub_shift0 = j_sub + (iteration)*dim_subspace
         j_sub_shift1 = j_sub + (iteration+1)*dim_subspace
      DO i_sub=1, dim_subspace
         i_sub_shift0 = i_sub + (iteration)*dim_subspace
         i_sub_shift1 = i_sub + (iteration+1)*dim_subspace
         final_hamiltonian(i_sub_shift0, j_sub_shift0)= A_matrix(iteration, i_sub, j_sub)
         final_hamiltonian(i_sub_shift0, j_sub_shift1)= B_matrix(iteration,i_sub, j_sub)
         final_hamiltonian(i_sub_shift1, j_sub_shift0)= CONJG (B_matrix(iteration,j_sub, i_sub))
      ENDDO
     ENDDO



     CALL calculate_variation( A_matrix(iteration,:,:), A_matrix(iteration-1,:,:), variation) 
     !
     !
      ! out process if it is the case
     final_iter = iteration
     var=ZERO
     DO i_sub = 1, dim_subspace
       var = ABS(variation(i_sub)) + var
     ENDDO
     CALL summary_recursion(stdout, iteration, A_matrix(iteration,:,:), B_matrix(iteration,:,:), variation, var, norm(iteration+1))

     IF ( ABS( var ) <  conv_criteria ) EXIT iteration_loop


  ENDDO iteration_loop

    ! init recursion_hamiltonian


 END SUBROUTINE recursion_loop



!*******************************************************************
   SUBROUTINE print_hamiltonian
   !*******************************************************************
       CHARACTER(17)      :: subname="print_hamiltonian"
       CHARACTER(11) :: name="HAMILTONIAN"
       CHARACTER(nstrx)   :: attr
       INTEGER            :: ik=1, ir=1, iopt,jopt
       INTEGER            :: ierr
      CHARACTER( LEN=nstrx )  :: filename


       CALL ioname('hamiltonian',filename)

       CALL file_open(ham_unit,TRIM(filename),PATH="/",ACTION="write", &
                              FORM='formatted')


       CALL iotk_write_begin(ham_unit,TRIM(name))

       CALL iotk_write_attr(attr,"dimwann",dimwan,FIRST=.TRUE.) 
       CALL iotk_write_attr(attr,"dim_rec",dim_recursion) 
       CALL iotk_write_attr(attr,"dim_sub",dim_subspace) 
       CALL iotk_write_empty(ham_unit,"DATADIM",ATTR=attr)

       CALL iotk_write_attr(attr,"max_iter",max_iter,FIRST=.TRUE.) 
       CALL iotk_write_attr(attr,"final_iter",final_iter) 
       CALL iotk_write_attr(attr,"dim_final",((max_iter+1)*dim_subspace)) 
       CALL iotk_write_empty(ham_unit,"DATAFINAL",ATTR=attr)

       !

       CALL iotk_write_begin(ham_unit,"TOTHAM")
       CALL iotk_write_dat(ham_unit,"TOT"//TRIM(iotk_index(ik)), total_hamiltonian(:,:))
       CALL iotk_write_end(ham_unit,"TOTHAM")
       !
       CALL iotk_write_begin(ham_unit,"RECHAM")
       CALL iotk_write_dat(ham_unit,"REC"//TRIM(iotk_index(ir)), recursion_hamiltonian(:,:))
       CALL iotk_write_end(ham_unit,"RECHAM")

       CALL iotk_write_begin(ham_unit,"FINHAM")
       CALL iotk_write_dat(ham_unit,"FIN"//TRIM(iotk_index(ik)), final_hamiltonian(:,:))
       CALL iotk_write_end(ham_unit,"FINHAM")

     ! to be added : optimized ham
    !       CALL iotk_write_begin(ham_unit,"OPTFINHAM")
    !
    !   CALL iotk_write_dat(ham_unit,"OFIN"//TRIM(iotk_index(ik)), (final_hamiltonian(iopt,jopt),&
    !                                      iopt=1,(final_iter*dim_subspace)
    !   ENDDO
    !   CALL iotk_write_end(ham_unit,"OPTFINHAM")


       CALL iotk_write_end(ham_unit,TRIM(name))

      CALL file_close(ham_unit,PATH="/",ACTION="write")

      CALL ioname('hamiltonian',filename,LPATH=.FALSE.)
      WRITE( stdout,"(/,'  Hamiltonian on WF and recursion bases written on file : '5x,a)") TRIM(filename)


END SUBROUTINE print_hamiltonian

!*******************************************************************
   SUBROUTINE print_matrix
   !*******************************************************************
       CHARACTER(12)      :: subname="print_matrix"
       CHARACTER(6) :: name="MATRIX"
       CHARACTER(nstrx)   :: attr
       INTEGER            :: iter
       INTEGER            :: ierr
      CHARACTER( LEN=nstrx )  :: filename


       CALL ioname('matrix',filename)

       CALL file_open(mat_unit,TRIM(filename),PATH="/",ACTION="write", &
                              FORM='formatted')


       CALL iotk_write_begin(mat_unit,TRIM(name))

       CALL iotk_write_attr(attr,"dimwann",dimwan,FIRST=.TRUE.) 
       CALL iotk_write_attr(attr,"dim_rec",dim_recursion) 
       CALL iotk_write_attr(attr,"dim_sub",dim_subspace) 
       CALL iotk_write_empty(mat_unit,"DATADIM",ATTR=attr)

       CALL iotk_write_attr(attr,"max_iter",max_iter,FIRST=.TRUE.) 
       CALL iotk_write_attr(attr,"final_iter",final_iter) 
       CALL iotk_write_attr(attr,"dim_final",((max_iter+1)*dim_subspace)) 
       CALL iotk_write_empty(mat_unit,"DATAFINAL",ATTR=attr)

       !
       CALL iotk_write_begin(mat_unit,"A_MAT")
       CALL iotk_write_dat(mat_unit,"ITER"//TRIM(iotk_index(0)), A_matrix_0(:,:))
       DO iter=1, max_iter
        CALL iotk_write_dat(mat_unit,"ITER"//TRIM(iotk_index(iter)), A_matrix(iter,:,:))
       ENDDO
       CALL iotk_write_end(mat_unit,"A_MAT")
       !
       !
       CALL iotk_write_begin(mat_unit,"B_MAT")
       CALL iotk_write_dat(mat_unit,"ITER"//TRIM(iotk_index(0)), B_matrix_0(:,:))

       DO iter=1, max_iter
        CALL iotk_write_dat(mat_unit,"ITER"//TRIM(iotk_index(iter)), B_matrix(iter,:,:))
       ENDDO
       CALL iotk_write_end(mat_unit,"B_MAT")
    ! 
    ! to be added : optimized ham
    !       CALL iotk_write_begin(mat_unit,"OPTFINHAM")
    !
    !   CALL iotk_write_dat(mat_unit,"OFIN"//TRIM(iotk_index(ik)), (final_hamiltonian(iopt,jopt),&
    !                                      iopt=1,(final_iter*dim_subspace)
    !   ENDDO
    !   CALL iotk_write_end(mat_unit,"OPTFINHAM")


       CALL iotk_write_end(mat_unit,TRIM(name))

      CALL file_close(mat_unit,PATH="/",ACTION="write")

      CALL ioname('matrix',filename,LPATH=.FALSE.)
      WRITE( stdout,"(/,'  Matrices on recursion basis written on file : ',5x,a)") TRIM(filename)


END SUBROUTINE print_matrix

!*******************************************************************
   SUBROUTINE print_states
   !*******************************************************************
       CHARACTER(12)      :: subname="print_states"
       CHARACTER(6) :: name="STATES"
       CHARACTER(nstrx)   :: attr
       INTEGER            :: iter, i_sub
       INTEGER            :: ierr
      CHARACTER( LEN=nstrx )  :: filename


       CALL ioname('states',filename)

       CALL file_open(sta_unit,TRIM(filename),PATH="/",ACTION="write", &
                              FORM='formatted')


       CALL iotk_write_begin(sta_unit,TRIM(name))

       CALL iotk_write_attr(attr,"dimwann",dimwan,FIRST=.TRUE.) 
       CALL iotk_write_attr(attr,"dim_rec",dim_recursion) 
       CALL iotk_write_attr(attr,"dim_sub",dim_subspace) 
       CALL iotk_write_empty(sta_unit,"DATADIM",ATTR=attr)

       CALL iotk_write_attr(attr,"max_iter",max_iter,FIRST=.TRUE.) 
       CALL iotk_write_attr(attr,"final_iter",final_iter) 
       CALL iotk_write_attr(attr,"dim_final",((max_iter+1)*dim_subspace)) 
       CALL iotk_write_empty(sta_unit,"DATAFINAL",ATTR=attr)

       !
       CALL iotk_write_begin(sta_unit,"STA")
       DO i_sub=1,dim_subspace
         CALL iotk_write_dat(sta_unit,"ITER"//TRIM(iotk_index(0))//TRIM(iotk_index(i_sub)), recursion_state_0(i_sub,:))
       ENDDO 
       DO iter=1, max_iter
         DO i_sub=1,dim_subspace
            CALL iotk_write_dat(sta_unit,"ITER"//TRIM(iotk_index(iter))//TRIM(iotk_index(i_sub)), recursion_state(iter,i_sub,:))
         ENDDO
       ENDDO
       !
       CALL iotk_write_end(sta_unit,"STA")

       CALL iotk_write_end(sta_unit,TRIM(name))

      CALL file_close(sta_unit,PATH="/",ACTION="write")

      CALL ioname('states',filename,LPATH=.FALSE.)
      WRITE( stdout,"(/,' Recursion states on WF basis written on file : ',5x,a)") TRIM(filename)


END SUBROUTINE print_states

!*******************************************************************
   SUBROUTINE print_wan_center
   !*******************************************************************
       CHARACTER(16)      :: subname="print_wan_center"
       CHARACTER(7) :: name="CENTERS"
       CHARACTER(nstrx)   :: attr
       INTEGER            :: iter, i_sub
       INTEGER            :: ierr
      CHARACTER( LEN=nstrx )  :: filename


       CALL ioname('wannier_center',filename)

       CALL file_open(cen_unit,TRIM(filename),PATH="/",ACTION="write", &
                              FORM='formatted')


       CALL iotk_write_begin(cen_unit,TRIM(name))

       CALL iotk_write_attr(attr,"dimwann",dimwan,FIRST=.TRUE.) 
       CALL iotk_write_attr(attr,"dim_rec",dim_recursion) 
       CALL iotk_write_attr(attr,"dim_sub",dim_subspace) 
       CALL iotk_write_empty(cen_unit,"DATADIM",ATTR=attr)

       CALL iotk_write_attr(attr,"max_iter",max_iter,FIRST=.TRUE.) 
       CALL iotk_write_attr(attr,"final_iter",final_iter) 
       CALL iotk_write_attr(attr,"dim_final",((max_iter+1)*dim_subspace)) 
       CALL iotk_write_empty(cen_unit,"DATAFINAL",ATTR=attr)

       DO i_sub=1, dim_subspace
         CALL iotk_write_attr(attr,"wan_num",wan_num(i_sub),FIRST=.TRUE.) 
         CALL iotk_write_attr(attr,"COOR",wannier_coordinates(:,wan_num(i_sub)), IERR=ierr) 
            IF (ierr/=0) CALL errore(subname,'writing wannier_coordinates of wan_num centers',ABS(ierr))
         CALL iotk_write_empty(cen_unit,"SUBCENTER"//TRIM(iotk_index(i_sub)),ATTR=attr)
       ENDDO

       !

       CALL iotk_write_dat(cen_unit,"REC2TOT", subrecur2tot(:))

       CALL iotk_write_dat(cen_unit,"INITCENTER", wannier_coordinates, COLUMNS=3, IERR=ierr) 
            IF (ierr/=0) CALL errore(subname,'writing wannier_coordinates',ABS(ierr))
      !

       CALL iotk_write_dat(cen_unit,"RECCENTER", wannier_coordinates_sub_recursion,  COLUMNS=3, IERR=ierr) 
            IF (ierr/=0) CALL errore(subname,'writing wannier_coordinates_sub_recursion',ABS(ierr))

       CALL iotk_write_end(cen_unit,TRIM(name))

       CALL file_close(cen_unit,PATH="/",ACTION="write")

       CALL ioname('wannier_center',filename,LPATH=.FALSE.)
       WRITE( stdout,"(/,' WF centers written on file : ',a)") TRIM(filename)


 END SUBROUTINE print_wan_center

 END MODULE recursion_module


