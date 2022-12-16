!      Copyright (C) 2007 WanT Group
!
!      This file is distributed under the terms of the
!      GNU General Public License. See the file `License\'
!      in the root directory of the present distribution,
!      or http://www.gnu.org/copyleft/gpl.txt .
!
!***********************************************
   PROGRAM current
   !***********************************************
   USE kinds,                ONLY : dbl
   USE parameters,           ONLY : nstrx 

!   USE T_smearing_module,    ONLY : smearing_init

   IMPLICIT NONE
   INTEGER, PARAMETER   ::  curr_unit =  11
   INTEGER, PARAMETER   ::  cond_unit = 12

   !
   REAL(dbl), ALLOCATABLE   :: curr(:)                 ! current
   REAL(dbl), ALLOCATABLE   :: curr2(:)                 ! current
   !REAL(dbl), ALLOCATABLE   :: diff_cond(:)                 ! current
   REAL(dbl)                :: mu_L, mu_R              ! chemical potentials
   REAL(dbl)                :: mu_L_aux, mu_R_aux      ! 
   REAL(dbl)                :: sigma                   ! broadening
   REAL(dbl), ALLOCATABLE   :: ftemp_L(:), ftemp_R(:)  ! temperature smearing functions
   REAL(dbl), ALLOCATABLE   :: transm(:)               ! transmittance from data file
   !
   CHARACTER(99)         :: filein                  ! input  filename (transmittance)
   CHARACTER(99)         :: fileout                 ! output filename (current)
   CHARACTER(99)         :: chr                 !
   CHARACTER(99)         :: chr2                 !


   !
   ! energy grid
   !
   INTEGER                  :: ne                      ! dimension of the energy grid
   REAL(dbl)                :: de         
   REAL(dbl)                :: de_old         
   REAL(dbl), ALLOCATABLE   :: egrid(:)                ! energy grid
   REAL(dbl)                :: E_fermi        
   REAL(dbl)                :: norm       
   REAL(dbl)                :: norm2
          !
   ! bias grid
   !
   INTEGER                  :: nV                      ! dimension of the bias grid
   REAL(dbl)                :: Vmin                    ! Vgrid extrema
   REAL(dbl)                :: Vmax                    !
   REAL(dbl)                :: dV         
   REAL(dbl), ALLOCATABLE   :: Vgrid(:)                ! bias grid

   !
   ! interpolation variables
   !
   INTEGER                  :: ndiv, ndim_new          ! interpolation grid dimension
   REAL(dbl)                :: de_new                  ! interpolation grid step 
   REAL(dbl), ALLOCATABLE   :: egrid_new(:)            ! interpolation energy grid
   REAL(dbl), ALLOCATABLE   :: transm_new(:)           ! interpolated transmittance
   INTEGER                  :: i_min, i_max, inew      ! 

   !
   ! local variables
   !
   INTEGER                  :: ie, iv, ierr, ios
   INTEGER                  :: i_start, i_end, ndim    ! integration extrema
   REAL(dbl), ALLOCATABLE   :: funct(:)                ! auxiliary vectors for integration

   !
   ! input 
!
!------------------------------
! main body
!------------------------------
!
!   CALL startup(version_number,'current')
      filein ='INPUTCURRENT'
!
! ... Read INPUT namelist from stdin
!
   mu_L                        =  -0.5
   mu_R                        =   0.5
   sigma                       =   0.0000001    ! eV
   Vmin                        =  -0.80
   Vmax                        =   0.80
   nV                          =  200
   E_fermi  = -4.554538523489600
!
! init
!
   !
   ! get energy grid and transmittance from data file
   !

    PRINT*, 'Beginning'
    PRINT*, "Opening File:", TRIM(filein)
   !
   OPEN ( cond_unit, FILE=TRIM(filein), FORM='formatted', IOSTAT=ierr )
   IF ( ierr/=0 ) THEN
        PRINT*, 'current','opening file = '//TRIM(filein), ABS(ierr)
        STOP
   ENDIF
    PRINT*, "Skipping first line"
    READ ( cond_unit, *, IOSTAT=ierr ) chr2
   !
   ie = 0
   !
   DO WHILE ( .TRUE. ) 
      !
      READ ( cond_unit, *, IOSTAT=ierr )
      !
      IF ( ierr /= 0 ) EXIT
      ie = ie + 1
      !
   ENDDO
   !
   ne = ie 
   !
    PRINT*, "Total number of lines=", ne
   ALLOCATE ( egrid(ne), STAT=ierr )
   IF( ierr /=0 ) THEN 
          PRINT*, 'current','allocating egrid', ABS(ierr)
          STOP 
   ENDIF !!!!!!'current','allocating egrid', ABS(ierr) )
   !
   ALLOCATE ( transm(ne), STAT=ierr )
   IF( ierr /=0 ) THEN 
          PRINT*, 'current','allocating transmittance', ABS(ierr)
          STOP 
   ENDIF !!!!!!'current','allocating transmittance', ABS(ierr) )
   !
   REWIND ( cond_unit )
   !
    PRINT*, "Skipping first line"
    READ ( cond_unit, *, IOSTAT=ierr ) chr2
    PRINT*, "Reading Transmission: ..."
   DO ie = 1, ne
       READ ( cond_unit, *, IOSTAT=ios ) egrid(ie), transm(ie)
       IF ( ios/=0 ) THEN
          PRINT*, 'current','reading T',ie, egrid(ie),  transm(ie)
          STOP 
       ENDIF !!!!!!'current','reading T',ie)
   ENDDO
   !
   CLOSE( cond_unit )
    PRINT*, "Reading Transmission: ... Done"


   DO ie = 1, ne
      egrid(ie) = egrid(ie) - E_fermi
   ENDDO
   !
   ! allocate
   !
   ALLOCATE ( Vgrid(nV), STAT=ierr )
   IF( ierr /=0 ) THEN 
          PRINT*, 'current','allocating Vgrid', ABS(ierr)
          STOP 
   ENDIF !!!!!!'current','allocating Vgrid', ABS(ierr) )
   !
   ALLOCATE ( curr(nV), STAT=ierr )
   IF( ierr /=0 ) THEN 
          PRINT*, 'current','allocating current', ABS(ierr)
          STOP 
   ENDIF !!!!!!'current','allocating current', ABS(ierr) )
  ALLOCATE ( curr2(nV), STAT=ierr )
   IF( ierr /=0 ) THEN 
          PRINT*, 'current','allocating current2', ABS(ierr)
          STOP 
   ENDIF !!!!!!'current','allocating current', ABS(ierr) )

   !ALLOCATE ( diff_cond(nV), STAT=ierr )
   !IF( ierr /=0 ) THEN 
   !       PRINT*, 'current','allocating diff_cond', ABS(ierr)
   !       STOP 
   !ENDIF !!!!!!'current','allocating diff_cond', ABS(ierr) )
    PRINT*, "Constructing Bias Grid: ..."
   !
   ! bias grid
   !
   dV = (Vmax - Vmin)/REAL(nV-1, dbl)
   !
   DO iv = 1, nV
      Vgrid(iv) = Vmin + REAL(iv-1, dbl) * dV
   ENDDO 
   !
    PRINT*, "Constructing Bias Grid: ... Done"
!
! current calculation
!
   !
   curr(:) = 0.0
   curr2(:) = 0.0
   de_old = (egrid(ne) - egrid(1))/REAL(ne-1, dbl)
   !
   PRINT*, "Integration: ... "
   DO iv = 1, nV
      !
      mu_L_aux = mu_L * Vgrid(iv)
      mu_R_aux = mu_R * Vgrid(iv)
      !
      ! integration extrema
      i_start=0
      i_end=0
      DO ie = 1, ne-1
      !CALL locate( egrid, ne, MIN( mu_L_aux, mu_R_aux ) -sigma -3.0_dbl*de_old, i_start)
           IF ((egrid(ie) < (MIN( mu_L_aux, mu_R_aux ) -sigma -3.0_dbl*de_old) ).AND.(egrid(ie+1) > (MIN( mu_L_aux, mu_R_aux ) -sigma -3.0_dbl*de_old) )) THEN
              i_start =ie
           ELSE IF ((egrid(ie) < (MAX( mu_L_aux, mu_R_aux ) +sigma +3.0_dbl*de_old) ).AND.(egrid(ie+1) > (MAX( mu_L_aux, mu_R_aux ) +sigma +3.0_dbl*de_old) )) THEN
              i_end = ie+1
           ENDIF

      !
      !CALL locate( egrid, ne, MAX( mu_R_aux, mu_L_aux ) +sigma +3.0_dbl*de_old, i_end)

      ENDDO
      IF ( i_start == 0 .OR. i_start == ne ) THEN
 
           PRINT*, 'current','invalid i_start',(MIN( mu_L_aux, mu_R_aux ) -sigma -3.0_dbl*de_old), i_start
           STOP
      ENDIF
      IF ( i_end == 0 .OR. i_end == ne ) THEN
          PRINT*, 'current','invalid i_end',5
          STOP
      ENDIF
      !
      !
      ! simpson routine requires that ndim is an odd number
      !IF ( MOD(ndim, 2) == 0 ) i_end = i_end - 1
      !
      ndim = i_end - i_start + 1
      de = (egrid(i_end) - egrid(i_start))/REAL(ndim-1, dbl)

      !
      ! redefinition of the integration mesh for a better interpolation
      !
      ndiv = NINT( de / (2.0_dbl*sigma) )
      IF (ndiv == 0) ndiv = 1      
      !
      de_new   = de / REAL(ndiv, dbl)
      ndim_new = (ndim - 1) * ndiv + 1 
      !
      ALLOCATE ( transm_new(ndim_new), STAT=ierr )
      IF( ierr /=0 ) THEN 
          PRINT*, 'current','allocating transm_new', ABS(ierr) 
          STOP
      ENDIF
      !
      ALLOCATE ( egrid_new(ndim_new), STAT=ierr )
      IF( ierr /=0 ) THEN 
          PRINT*, 'current','allocating egrid_new', ABS(ierr) 
          STOP
      ENDIF

      !
      ! new integration mesh
      egrid_new(1) = egrid(i_start)
      !
      DO inew = 2, ndim_new   
         egrid_new(inew) = egrid_new(1) + de_new * REAL( inew -1, dbl)
      ENDDO
      !
      ! Transmittance interpolated on the new grid
      !
      IF (ndiv /= 1) THEN
         !
         DO ie = i_start, i_end - 1
            !
            i_min = (ie-i_start) * ndiv + 1 
            i_max = (ie-i_start+1) * ndiv + 1 
            !
            transm_new(i_min) = transm(ie)
            transm_new(i_max) = transm(ie+1)
            !
            DO inew = i_min+1, i_max-1
               !
               transm_new(inew) = transm_new(i_max)*(egrid_new(inew)-egrid_new(i_min))/  &
                                  (egrid_new(i_max)-egrid_new(i_min)) -                  &
                                  transm_new(i_min)*(egrid_new(inew)-egrid_new(i_max))/  &
                                  (egrid_new(i_max)-egrid_new(i_min))
            ENDDO
         ENDDO
         !
      ELSE 
         !
         ! ndiv == 1
         !
         transm_new(:) = transm(i_start:i_end)
      ENDIF
      !
      ! auxiliary vectors for integral calculation
      !
      ALLOCATE ( funct(ndim_new), STAT=ierr )
      IF( ierr /=0 ) THEN 
          STOP 
          ENDIF !!!!!!'current','allocating funct', ABS(ierr) )
      !
      ALLOCATE ( ftemp_L(ndim_new), ftemp_R(ndim_new), STAT=ierr )
      IF( ierr /=0 ) THEN 
          STOP 
          ENDIF !!!!!!'current','allocating ftemp', ABS(ierr) )
      !
      !
      DO ie = 1, ndim_new
          !
          ftemp_L(ie) = 1.0 / ( EXP( -(egrid_new( ie ) -mu_L_aux ) / sigma) + 1.0 )
          ftemp_R(ie) = 1.0 / ( EXP( -(egrid_new( ie ) -mu_R_aux ) / sigma) + 1.0 )
          !
      ENDDO
      !
      ! perform the integration
      !
      ! if you want to use the simpson routine for integration 
      ! uncomment the following line and comment the next one   
      !
!      funct(:) = ( ftemp_L(:) - ftemp_R(:) ) * transm(i_start:i_end)
      funct(:) = transm_new(1:ndim_new)
      !
      ! if you want to use the simpson routine for integration 
      ! uncomment the following line and comment the next ones   
      !
!      CALL simpson (ndim, funct, rab, curr(iv) )
      DO ie = 1, ndim_new-1
         curr(iv) = curr(iv) + ( ftemp_L(ie) - ftemp_R(ie) )*funct(ie)*de_new/3.0 + &
                    ( ftemp_L(ie+1) - ftemp_R(ie+1) )*funct(ie+1)*de_new/3.0 +      &
                    ( ftemp_L(ie) - ftemp_R(ie) )*funct(ie+1)*de_new/6.0 +          &
                    ( ftemp_L(ie+1) - ftemp_R(ie+1) )*funct(ie)*de_new/6.0
         curr2(iv) = curr2(iv) + ( ftemp_L(ie) - ftemp_R(ie) )*funct(ie)*de_new
      ENDDO

      !
      DEALLOCATE ( transm_new, STAT=ierr )
      IF( ierr /=0 ) THEN 
          PRINT*, 'current','deallocating transm_new', ABS(ierr)
          STOP 
      ENDIF !!!!!!'current','deallocating transm_new', ABS(ierr) )
      !
      DEALLOCATE ( egrid_new, STAT=ierr )
      IF( ierr /=0 ) THEN 
          PRINT*, 'current','deallocating transm_new', ABS(ierr)
          STOP 
      ENDIF !!!!!!'current','deallocating transm_new', ABS(ierr) )
      DEALLOCATE ( funct, STAT=ierr )
      IF( ierr /=0 ) THEN 
          PRINT*, 'current','deallocating funct', ABS(ierr)
          STOP 
      ENDIF !!!!!!'current','deallocating funct', ABS(ierr) )
      DEALLOCATE ( ftemp_L, ftemp_R, STAT=ierr )
      IF( ierr /=0 ) THEN 
          PRINT*, 'current','deallocating ftemp', ABS(ierr)
          STOP 
      ENDIF !!!!!!'current','deallocating ftemp', ABS(ierr) )
      !
   ENDDO
   PRINT*, "Integration: ... Done"
   !
   !     fileout = 'diff_cond.dat'
  !OPEN ( curr_unit, FILE=TRIM(fileout), FORM='formatted' )
   !
   !norm=(curr(INT(nV/2.0) +1)- curr(INT(nV/2.0)-1))/(2*dV)
   !norm2=(curr2(INT(nV/2.0) +1)- curr2(INT(nV/2.0)-1))/(2*dV)
   !DO iv = 2, nV-1
      !
    !  diff_cond(iv)=  (curr(iv+1)- curr(iv-1))/(2*dV)

 
  !    WRITE ( curr_unit, '(4(f15.9))' ) Vgrid(iv)*1000,  ((curr(iv+1)- curr(iv-1))/(2*dV*norm)), diff_cond(iv),    ((curr2(iv+1)- curr2(iv-1))/(2*dV*norm2))
  ! ENDDO
   !
  ! CLOSE( curr_unit )

   !
   ! write input data on the output file
   !

        fileout = 'current.dat'
   !
   OPEN ( curr_unit, FILE=TRIM(fileout), FORM='formatted' )
   !
   DO iv = 1, nV
       WRITE ( curr_unit,*) Vgrid(iv)*1000, " ", (curr(iv) * (1000.0 / (3.1415*27.2113845)) *1.60217653*1000.0*4.35974417 / 1.05457168 ) !nA
   ENDDO
   !
   CLOSE( curr_unit )
   !

   
   !
   ! deallocate
   !
   DEALLOCATE ( egrid, STAT=ierr )
   IF( ierr /=0 ) THEN 
          PRINT*, 'current','deallocating egrid', ABS(ierr) 
          STOP 
   ENDIF !!!!!!'current','deallocating egrid', ABS(ierr) )
   !
   DEALLOCATE ( Vgrid, STAT=ierr )
   IF( ierr /=0 ) THEN 
          PRINT*, 'current','deallocating Vgrid', ABS(ierr) 
          STOP 
   ENDIF !!!!!!'current','deallocating Vgrid', ABS(ierr) )
   !
   DEALLOCATE ( curr, STAT=ierr )
   IF( ierr /=0 ) THEN 
          PRINT*, 'current','deallocating current', ABS(ierr)
          STOP 
   ENDIF !!!!!!'current','deallocating current', ABS(ierr) )
   !
   DEALLOCATE ( curr2, STAT=ierr )
   IF( ierr /=0 ) THEN 
          PRINT*, 'current','deallocating current', ABS(ierr)
          STOP 
   ENDIF !!!!!!'current','deallocating current', ABS(ierr) )
   !
  ! DEALLOCATE ( diff_cond, STAT=ierr )
  ! IF( ierr /=0 ) THEN 
  !        PRINT*, 'current','deallocating diff_cond', ABS(ierr) 
  !        STOP 
  ! ENDIF !!!!!!'current','deallocating diff_cond', ABS(ierr) )
   !
   DEALLOCATE ( transm, STAT=ierr )
   IF( ierr /=0 ) THEN 
          PRINT*, 'current','deallocating transmittance', ABS(ierr) 
          STOP 
   ENDIF !!!!!!'current','deallocating transmittance', ABS(ierr) )

   !

END PROGRAM current
  
