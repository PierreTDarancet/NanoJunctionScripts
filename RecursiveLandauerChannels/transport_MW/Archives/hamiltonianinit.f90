







!***********************************************
   SUBROUTINE hamiltonian_init()
   !***********************************************
  INTEGER, PARAMETER    ::   Nb_states=4
  INTEGER, PARAMETER    ::   Nb_states_lead=1

   COMPLEX(dbl), INTENT(in)       ::   HOMO_left, HOMO_right
   COMPLEX(dbl)       ::   LUMO_left, LUMO_right
   COMPLEX(dbl)    ::   H_lL
   COMPLEX(dbl)    ::     H_lR
   COMPLEX(dbl)    ::    H_rL
   COMPLEX(dbl)    ::    H_rR


!
!------------------------------
! main body
!------------------------------
!
!
! init
!

   nbias=100
   biasmin=-4.00 !!!!!!!!!!!!!!! biasmin = - biasmax
   biasmax=4.00 !!!!!!!!!!!!!!!!! otherwise RR is false

   delta=0.00001 !idelta for gf


   PRINT*, 'End of grid parameter initialization'

   PRINT*, 'Beginning grid allocations'

   !
   ! bias grid
   !


      IF ( (nbias<0)) THEN
          STOP 
       END IF 
       ALLOCATE( biasgrid(nbias), STAT=ierr )

       !
       ! setting the energy grid
       !
       dbias = (biasmax - biasmin) / REAL(nbias -1, dbl)
       !
       DO ibias = 1, nbias
          biasgrid(ibias) = biasmin + REAL(ibias -1, dbl) * dbias
       ENDDO


   !
   ! energy grid
   !
       ne= 1001
!       !ne= nbias
       IF ((ne<0)) THEN
          STOP
       END IF 
       !
!       de = (biasmax) / REAL(nbias -1, dbl)

!       de = (biasmax - biasmin) / REAL(nbias -1, dbl)

       de = (biasmax) / REAL(ne-1, dbl)
       !
       ALLOCATE( egrid(ne), STAT=ierr )!!!!!!!taille nb bias max 

!!1       egridalloc = .TRUE.

       !
       ! setting the energy grid


   !
   ! Voc grid
   !

       IF ( (nvoc<0)) THEN
          STOP
       END IF 
       ALLOCATE( vocgrid(nvoc), STAT=ierr )

       !
       ! setting the energy grid
       !
       dvoc = (vocmax - vocmin) / REAL(nvoc -1, dbl)
       !
       DO ivoc = 1, nvoc
          vocgrid(ivoc) = vocmin + REAL(ivoc -1, dbl) * dvoc
       ENDDO

   !!    vocgridalloc = .TRUE.


   PRINT*, 'End of grid allocations'

   ! Set up the hamiltonian
   !

   PRINT*, 'Allocate Hamiltonian'
 !  ALLOCATE ( gamma_L_safe(nvoc,nbias,ne,Nb_states,Nb_states), STAT=ierr )
 !  ALLOCATE ( gamma_R_safe(nvoc,nbias,ne,Nb_states,Nb_states), STAT=ierr )

  ! ALLOCATE ( conduct(nvoc,nbias,Nb_states,ne), STAT=ierr )
  ! ALLOCATE ( dos(nvoc,nbias, Nb_states,ne), STAT=ierr )


   ALLOCATE( h00_C(Nb_states,Nb_states), STAT=ierr )
   ALLOCATE( h_CR(Nb_states,1), STAT=ierr )
   ALLOCATE( h_LC(1,Nb_states), STAT=ierr )
   ALLOCATE( h00_L(1,1), STAT=ierr )
   ALLOCATE( h00_R(1,1), STAT=ierr )
   ALLOCATE( h01_L(1,1), STAT=ierr )
   ALLOCATE( h01_R(1,1), STAT=ierr )
   ALLOCATE( gamma_R(Nb_states,Nb_states), STAT=ierr )
   ALLOCATE( sigma_R(Nb_states,Nb_states), STAT=ierr )
   ALLOCATE( gamma_L(Nb_states,Nb_states), STAT=ierr )
   ALLOCATE( sigma_L(Nb_states,Nb_states), STAT=ierr )

   ALLOCATE ( total_conduct(ne), STAT=ierr )

   ALLOCATE ( cond_aux(Nb_states), STAT=ierr )
   ALLOCATE ( work(Nb_states,Nb_states), STAT=ierr )
   ALLOCATE ( work2(Nb_states,Nb_states), STAT=ierr )
   ALLOCATE (  gC(Nb_states,Nb_states), STAT=ierr )
   ALLOCATE (  aux00_C(Nb_states,Nb_states), STAT=ierr )
   ALLOCATE ( aux00_R(1,1), STAT=ierr )
   ALLOCATE (  aux00_L(1,1), STAT=ierr )
   ALLOCATE (aux_LC(1,Nb_states), STAT=ierr )
   ALLOCATE ( aux_CL(Nb_states,1), STAT=ierr )
   ALLOCATE ( aux_RC(1,Nb_states), STAT=ierr )
   ALLOCATE ( aux_CR(Nb_states,1), STAT=ierr )
   ALLOCATE (aux01_R(1,1), STAT=ierr )
   ALLOCATE (aux01_L(1,1), STAT=ierr )
   ALLOCATE (  Current(nvoc,nbias), STAT=ierr )
   ALLOCATE (  RR(nvoc,(INT(nbias/2))), STAT=ierr )
   ALLOCATE ( s00_C(Nb_states,Nb_states), STAT=ierr )


  ! dos(:,:,:,:) = ZERO
  ! conduct(:,:,:,:) = ZERO
  ! gamma_L_safe(:,:,:,:,:) = CZERO
  ! gamma_R_safe(:,:,:,:,:) = CZERO
   PRINT*, 'End of allocations'

 
   ! define unity matrix
   s00_C(:,:) = CZERO
 
   DO i=1, Nb_States
      s00_C(i,i) = CONE
   ENDDO 

    ! define leads

   h00_L(1,1)=CZERO
   h00_R(1,1)=CZERO


   h01_L(1,1)=4.000*CONE
   h01_R(1,1)=4.000*CONE
   !!!!!!!!VOC LOOP

   Do ivoc = 1, nvoc    

     HOMO_right = HOMO_left+vocgrid(ivoc)             ! indicated wrt to \mu_R
     LUMO_right = LUMO_left+vocgrid(ivoc)             ! Let's take Right as the Donor

     PRINT*, "Voc", vocgrid(ivoc)            

   !!!!!!BIAS LOOP


   Do ibias = 1, nbias    

        PRINT*, "Beginning of the bias loop"
        PRINT*, "ibias=", ibias, "bias=", biasgrid(ibias)


   h00_C(:,:) = CZERO
   h00_C(1,1) = HOMO_left-biasgrid(ibias)/2
   h00_C(2,2) = LUMO_left-biasgrid(ibias)/2
   h00_C(3,3) = HOMO_right+biasgrid(ibias)/2
   h00_C(4,4) = LUMO_right+biasgrid(ibias)/2

 
   h_CR(1,1) = H_lR
   h_CR(2,1) = H_lR
   h_CR(3,1) = H_rR 
   h_CR(4,1) = H_rR 
   h_LC(1,1) = H_lL
   h_LC(1,2) = H_lL
   h_LC(1,3) = H_rL
   h_LC(1,4) = H_rL

!PRINT*, h00_C(:,:)

   Current(ivoc,ibias)=ZERO

 !!!!!!!!!!DEFINE NUMBER OF TRANSMISSION POINT NEEDED IN HTE BIAS WINDOW
  !ne_bias=ABS(INT(biasgrid(ibias)/de) + 1 )

!  ne_bias_sum(ibias) = ne_bias

!   PRINT*, "ne_bias = ", ne_bias
 
      de = ABS(biasgrid(ibias))/ REAL(ne-1, dbl)
   PRINT*, "de= ", de
       egrid(:) = ZERO
       DO ie = 1, ne
!       DO ie = 1, ne_bias
          egrid(ie) = -ABS(biasgrid(ibias))/2.00 + REAL(ie -1, dbl) * de
       ENDDO


  !     IF (biasgrid(ibias) == )
!
! main loop over frequency
! 

 !  PRINT*, "beginning energy loop"
 
!   energy_loop: &
!   DO ie = 1, ne_bias


   energy_loop: &
   DO ie = 1, ne


      ncount = ie
      !
      ! grids and misc
      !
      ene =  egrid(ie)  + delta * CI

   !   dos(ivoc,ibias,:,ie) = ZERO
   !   conduct(ivoc,ibias,:,ie) = ZERO
          ! init
          !
          aux00_L(:,:)  = h00_L(:,:) - (biasgrid(ibias)/2.00) !-ene * s00_L(:,:)
          aux01_L(:,:)  = h01_L(:,:) 
          !
          aux00_R(:,:)  = h00_R(:,:) + (biasgrid(ibias)/2.00) ! -ene * s00_R(:,:)
          aux01_R(:,:)  = h01_R(:,:) 
          !
          aux00_C(:,:)  = h00_C(:,:) -ene * s00_C(:,:) 
          aux_LC(:,:) = h_LC(:,:)
          aux_CR(:,:) = h_CR(:,:)
          !
          aux_CL(:,:) = CONJG( TRANSPOSE( h_LC(:,:) ))
          aux_RC(:,:) = CONJG( TRANSPOSE( h_CR(:,:) ))
          ! 
 
          !IF (ibias>=500) THEN
          !     PRINT*, "Energy", ie,  egrid(ie), REAL(ene)
          !     PRINT*, "biasgrid(ibias)",   (biasgrid(ibias)/2.0)
          !ENDIF


 
          ! construct leads self-energies 
          ! 
          ! ene +bias/2
          ! Right
           gR(:,:) = CZERO
           spectra_min= aux00_R(1,1)  - 2 * ABS (aux01_R(1,1))
           spectra_max= aux00_R(1,1)  + 2 * ABS (aux01_R(1,1))
          !
         ! IF (ibias>=500) THEN
         !      PRINT*, "Spectra_min",           spectra_min
         !      PRINT*, "Spectra_max",          spectra_max
         ! ENDIF

           IF ( ( REAL(ene) <= spectra_max ) .AND. ( REAL (ene) >=  spectra_min )) THEN
                !
                ! b² - 4ac
                !
                delta_spectra = (aux00_R(1,1) - ene)**2 - (4 * (aux01_R(1,1))**2)
                IF (REAL(delta_spectra) > ZERO )  PRINT*, "POSITIVE DELTA!!!!!!!!!!!" 
                !
                work_scal = (  -(aux00_R(1,1)- ene) - CI* SQRT ( -delta_spectra ) ) / (2* (aux01_R(1,1)**2))
                !
                gR(1,1) = CONE / (ene - aux00_R(1,1) - ((aux01_R(1,1) **2) * work_scal))
                !
            ELSE 
                !
                PRINT*, "OUt of spectrum R"
                gR(1,1) = CZERO
                PRINT*, "spectra_min"
                PRINT*, spectra_min
                PRINT*, "spectra_max"
                PRINT*, spectra_max
                PRINT*, REAL( ene )
                PRINT*, ie 
                PRINT*, egrid(ie) 
                !
            ENDIF

  
          ! ene - bias/2
          ! Left

           gL(:,:) = CZERO
           spectra_min= aux00_L(1,1)  - 2 * ABS (aux01_L(1,1))
           spectra_max= aux00_L(1,1)  + 2 * ABS (aux01_L(1,1))
      !
  
        !IF (ibias>=500) THEN
        !      PRINT*, "Left"
        !       PRINT*, "Spectra_min",           spectra_min
        !       PRINT*, "Spectra_max",          spectra_max
        !  ENDIF

         IF ( ( REAL(ene) <= spectra_max ) .AND. ( REAL (ene) >=  spectra_min )) THEN
                !
                ! b² - 4ac
                !
                delta_spectra = (aux00_L(1,1) - ene)**2 - (4 * (aux01_L(1,1))**2)
                IF (REAL(delta_spectra) > ZERO )  PRINT*, "POSITIVE DELTA!!!!!!!!!!!" 
                !
                work_scal = (  -(aux00_L(1,1)- ene) - CI* SQRT ( -delta_spectra ) ) / (2* (aux01_L(1,1)**2))
                !
                gL(1,1) = CONE / (ene - aux00_L(1,1) - ((aux01_L(1,1) **2) * work_scal))
                !
            ELSE 
                !
                PRINT*, "OUt of spectrum L"
              PRINT*, "spectra_min"
                PRINT*, spectra_min
                PRINT*, "spectra_max"
                PRINT*, spectra_max
  

                gL(1,1) = CZERO
                !
            ENDIF

            !
            !

 
        !IF (ibias>=500) THEN
        !      PRINT*, "End of g_leads"
        !  ENDIF

             !
          CALL mat_mul(work, aux_CR, 'N', gR, 'N', Nb_states, 1,1)
          CALL mat_mul(sigma_R, work, 'N', aux_RC, 'N', Nb_states,Nb_states, 1)
           !
          CALL mat_mul(work, aux_CL, 'N', gL, 'N', Nb_states, 1, 1)
          CALL mat_mul(sigma_L, work, 'N', aux_LC, 'N',Nb_states,Nb_states, 1) 
 
          !
          ! gamma_L and gamma_R
          !
          gamma_L(:,:) = CI * (  sigma_L(:,:) - CONJG( TRANSPOSE(sigma_L(:,:)) )   )
          gamma_R(:,:) = CI * (  sigma_R(:,:) - CONJG( TRANSPOSE(sigma_R(:,:)) )   )

    !      gamma_L_safe(ivoc,ibias,ie,:,:) = gamma_L(:,:)
    !      gamma_R_safe(ivoc,ibias,ie,:,:) = gamma_R(:,:)
          !
          ! Construct the conductor green's function
          ! gC = work^-1  (retarded)
          !


        !IF (ibias>=500) THEN
        !      PRINT*, "End of leads"
        !  ENDIF

          work(:,:) = -aux00_C(:,:) -sigma_L(:,:) -sigma_R(:,:)

          gC(:,:) = CZERO
          DO i = 1, Nb_states
             gC(i,i)= CONE
          ENDDO
 


          CALL mat_sv(Nb_states,Nb_states, work, gC)
          !




        !IF (ibias>=500) THEN
        !      PRINT*, "End of diag"
        !  ENDIF

          ! Compute density of states for the conductor layer
          !
     !     DO i = 1, Nb_states
     !        dos(ivoc,ibias,i,ie) = dos(ivoc,ibias,i,ie) -(AIMAG( gC(i,i) ) / PI)
     !     ENDDO
          cond_aux(:)=ZERO
          !
          ! evaluate the transmittance according to the Fisher-Lee formula


        !IF (ibias>=500) THEN
        !      PRINT*, "End of DOS"
        !  ENDIF

          

          !
          ! gL * gintr -> tmp
          !
          CALL mat_mul(work, gamma_L, 'N', gC, 'N', Nb_states,  Nb_states,  Nb_states)
          !
           ! gL * gintr * gR -> tmp1
          !
          CALL mat_mul(work2, work, 'N', gamma_R, 'N', Nb_states,  Nb_states,  Nb_states)
          !
          ! gL * gintr * gR * lambda * ginta -> work
          !
          CALL mat_mul(work, work2, 'N', gC, 'C', Nb_states,  Nb_states,  Nb_states)
          !
          DO i=1,Nb_states
             cond_aux(i) = REAL( work(i,i) )
          ENDDO
          !



        !IF (ibias>=500) THEN
       !       PRINT*, "End of mat mul"
       !   ENDIF

      !    conduct(ivoc,ibias,:,ie) =  cond_aux(:)

          total_conduct(ie) = SUM(cond_aux(:))

        !IF (ibias>=500) THEN
         !     PRINT*, "End of conduct"
         ! ENDIF

           ! Integrate the transmission
!          Current(ivoc,ibias)=Current(ivoc,ibias) + (total_conduct(ie)/REAL(ne_bias,dbl)) *de

          Current(ivoc,ibias)=Current(ivoc,ibias) + (total_conduct(ie)/REAL(ne,dbl)) *de
          !

        !IF (ibias>=500) THEN
         !     PRINT*, "End of loop"
         ! ENDIF

   ENDDO energy_loop

          Current(ivoc,ibias)=Current(ivoc,ibias) * ( biasgrid(ibias) / (ABS(biasgrid(ibias))) )

           PRINT*, "End of the bias loop"
           PRINT*, "ibias=", ibias, "bias=", biasgrid(ibias)


   ENDDO

           PRINT*, "Calculating RR"
!!!!!ENDBIAS LOOP
    DO ibias=1, INT(nbias/2)  !!!!!!!1bias is assumed to be symmetric !!!!!!
       RR(ivoc, ibias)=ABS( Current(ivoc,(nbias - ibias+1)) / Current(ivoc,ibias) )
    ENDDO 

!!!!!!!!!!!!!!!!
           PRINT*, "Calculating RR"
           PRINT*, "End of the Voc loop"
   ENDDO
!!!!!ENDVOC LOOP







!
! ... write DOS and CONDUCT data on files
!
! AND OTHER STUFF

  PRINT*, "Writing results"


   OPEN ( 12, FILE='Current.dat', FORM='formatted' )
   DO ivoc = 1, nvoc
       work_scal2 = LUMO_left - HOMO_left - vocgrid(ivoc) !!!!!!
   DO ibias = 1, nbias
       Current(ivoc,ibias) = Current(ivoc,ibias) * 6623617.82 / (27.2113845)   !!! nA
       WRITE (12, '(3(f15.9))' ) work_scal2, biasgrid(ibias),   Current(ivoc,ibias)
   ENDDO
       WRITE (12, '(2(e15.9))' )

   ENDDO
   CLOSE( 12 )
   OPEN ( 13, FILE='RectificationRatio.dat', FORM='formatted' )
   DO ivoc = 1, nvoc
      work_scal2 = LUMO_left - HOMO_left - vocgrid(ivoc) !!!!!!
      DO ibias = 1, INT(nbias/2)
       WRITE (13, '(3(f15.9))' ) work_scal2, ABS(biasgrid(ibias)),   RR(ivoc,ibias)
      ENDDO
       WRITE (13, '(2(e15.9))' )

   ENDDO
   CLOSE( 13 )
!   OPEN ( 12, FILE='dos.dat', FORM='formatted' )
!   DO ivoc = 1, nvoc
!   DO ibias = 1, nbias
!       DO ie = 1, ne
!           WRITE (12, '(10(f15.9))' ) vocgrid(ivoc), biasgrid(ibias), egrid(ie), SUM(dos(ivoc,ibias,:,ie)), (dos(ivoc,ibias,i,ie), i=1,Nb_states)
!       ENDDO
!       WRITE (12, '(2(e15.9))' )
!   ENDDO
!       WRITE (12, '(2(e15.9))' )

!   ENDDO
!   CLOSE( 12 )

!   OPEN ( 12, FILE='conduct.dat', FORM='formatted' )
!   DO ivoc = 1, nvoc
!   DO ibias = 1, nbias
!       DO ie = 1, ne
!           WRITE (12, '(10(f15.9))' ) vocgrid(ivoc), biasgrid(ibias), egrid(ie), SUM(conduct(ivoc,ibias,:,ie)), (conduct(ivoc,ibias,i,ie), i=1,Nb_states)
!       ENDDO
!       WRITE (12, '(2(e15.9))' )
!   ENDDO
!       WRITE (12, '(2(e15.9))' )

!   ENDDO
!   CLOSE( 12 )

 ! OPEN ( 12, FILE='gamma_L.dat', FORM='formatted' )
 !  DO ivoc = 1, nvoc
 !  DO ibias = 1, nbias
 !      DO ie = 1, ne
 !          WRITE (12, '(10(f15.9))' ) vocgrid(ivoc), biasgrid(ibias), egrid(ie), (gamma_L_safe(ivoc,ibias,ie,i,i), i=1,Nb_states)
 !      ENDDO
 !      WRITE (12, '(2(e15.9))' )
 !  ENDDO
 !      WRITE (12, '(2(e15.9))' )

!   ENDDO
!   CLOSE( 12 )


!  OPEN ( 12, FILE='gamma_L_nondiag.dat', FORM='formatted' )
!   DO ivoc = 1, nvoc
!   DO ibias = 1, nbias
!       DO ie = 1, ne
!           WRITE (12, '(10(f15.9))' ) vocgrid(ivoc), biasgrid(ibias), egrid(ie), (gamma_L_safe(ivoc,ibias,ie,i,1), i=1,Nb_states)
!           WRITE (12, '(10(f15.9))' ) vocgrid(ivoc), biasgrid(ibias), egrid(ie), (gamma_L_safe(ivoc,ibias,ie,i,2), i=1,Nb_states)
!           WRITE (12, '(10(f15.9))' ) vocgrid(ivoc), biasgrid(ibias), egrid(ie), (gamma_L_safe(ivoc,ibias,ie,i,3), i=1,Nb_states)
!       ENDDO
!       WRITE (12, '(2(e15.9))' )
!   ENDDO
!       WRITE (12, '(2(e15.9))' )

!   ENDDO
!   CLOSE( 12 )


 ! OPEN ( 12, FILE='gamma_R.dat', FORM='formatted' )
 !  DO ivoc = 1, nvoc
 !  DO ibias = 1, nbias
 !      DO ie = 1, ne
  !         WRITE (12, '(10(f15.9))' ) vocgrid(ivoc), biasgrid(ibias), egrid(ie), (gamma_R_safe(ivoc,ibias,ie,i,i), i=1,Nb_states)
!       ENDDO
!       WRITE (12, '(2(e15.9))' )
!   ENDDO
!       WRITE (12, '(2(e15.9))' )

!   ENDDO
!   CLOSE( 12 )


!   OPEN ( 12, FILE='conduct.dat', FORM='formatted' )
!   DO ivoc = 1, nvoc
!   DO ibias = 1, nbias
!       DO ie = 1, ne
!           WRITE (12, '(10(f15.9))' ) vocgrid(ivoc), biasgrid(ibias), egrid(ie), SUM(conduct(ivoc,ibias,:,ie)), (conduct(ivoc,ibias,i,ie), i=1,Nb_states)
!       ENDDO
!       WRITE (12, '(2(e15.9))' )
!   ENDDO
!       WRITE (12, '(2(e15.9))' )

!   ENDDO
!   CLOSE( 12 )


!
!...  free memory
!

  PRINT*, "Deallocate Hamiltonian"
  DEALLOCATE ( h00_C, STAT=ierr )
  DEALLOCATE ( h_CR , STAT=ierr )
  DEALLOCATE ( h_LC , STAT=ierr )
  DEALLOCATE ( h00_L , STAT=ierr )
  DEALLOCATE ( h00_R , STAT=ierr )
  DEALLOCATE ( h01_L , STAT=ierr )
  DEALLOCATE ( h01_R , STAT=ierr )
  PRINT*, "Deallocate grids"
  DEALLOCATE ( vocgrid , STAT=ierr )
  DEALLOCATE ( biasgrid , STAT=ierr ) 
  DEALLOCATE ( egrid , STAT=ierr )

  PRINT*, "Deallocate transport variables"
  DEALLOCATE ( gamma_R, STAT=ierr )
  DEALLOCATE ( sigma_R , STAT=ierr )
  DEALLOCATE ( gamma_L , STAT=ierr )
  DEALLOCATE ( sigma_L , STAT=ierr )
!  DEALLOCATE ( dos , STAT=ierr )
!  DEALLOCATE ( conduct , STAT=ierr )

 ! DEALLOCATE ( gamma_L_safe, STAT=ierr )
 ! DEALLOCATE ( gamma_R_safe , STAT=ierr )
  DEALLOCATE ( total_conduct , STAT=ierr )
  DEALLOCATE ( cond_aux , STAT=ierr )
  DEALLOCATE ( work , STAT=ierr )
  DEALLOCATE ( work2 , STAT=ierr )
  DEALLOCATE ( gC , STAT=ierr )
  DEALLOCATE ( aux00_C, STAT=ierr )
  DEALLOCATE ( aux00_R, STAT=ierr )
  DEALLOCATE ( aux00_L, STAT=ierr )
  DEALLOCATE ( aux_LC, STAT=ierr )
  DEALLOCATE ( aux_CL, STAT=ierr )
  DEALLOCATE ( aux_RC, STAT=ierr ) 
  DEALLOCATE ( aux_CR, STAT=ierr )
  DEALLOCATE ( aux01_R, STAT=ierr )
  DEALLOCATE ( aux01_L, STAT=ierr )

