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
   MODULE  transport_formula_module 
   !***********************************************
  USE kinds, ONLY : dbl
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI, EPS_m3, EPS_m5
  USE util_module,          ONLY : mat_mul, mat_sv
   IMPLICIT NONE
   PRIVATE 
   SAVE
!


   ! Public
   ! Private
   ! Public variables

   ! Public routines:
   PUBLIC                 :: photocurrent_calculation
   PUBLIC                 :: transmittance_calculation 
   PUBLIC                 :: current_calculation 
			 
   CONTAINS 
!***********************************************
   SUBROUTINE  transmittance_calculation(conductance, gC, gamR, gamL, dimx)
   !***********************************************
   ! Calculate the surface green function associated with a 1D semi infinite and periodic chain
  IMPLICIT NONE
      INTEGER, INTENT(in)   :: dimx
      REAL(dbl), INTENT(out) ::  conductance

      COMPLEX(dbl), INTENT(in) ::  gC(dimx,dimx) 
      COMPLEX(dbl), INTENT(in) ::  gamR(dimx,dimx) 
      COMPLEX(dbl), INTENT(in) ::  gamL(dimx,dimx) 

      ! Local
       CHARACTER(25)      :: subname="transmittance_calculation"
       INTEGER :: ierr, idiag
       COMPLEX(dbl)           :: tmp(dimx,dimx), tmp1(dimx,dimx)
       REAL(dbl) ::  conduct(dimx)

        IF( dimx < 1 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in Conductor dimension", " dimC  ", dimx
            STOP
         ENDIF 

   !
   ! gL * gintr -> tmp
   ! 
   CALL mat_mul(tmp, gamL, 'N', gC, 'N', dimx, dimx, dimx)
   !
   ! gL * gintr * gR -> tmp1
   !
   CALL mat_mul(tmp1, tmp, 'N', gamR, 'N', dimx, dimx, dimx)
   !
   ! gL * gintr * gR * ginta  -->  tmp
   !
   CALL mat_mul(tmp, tmp1, 'N', gC, 'C', dimx, dimx, dimx)
   !
       
   conductance = ZERO
   DO idiag=1,dimx
      conduct(idiag) = REAL( tmp(idiag,idiag) )
   ENDDO
   conductance = SUM(conduct(:))
       !
    END SUBROUTINE transmittance_calculation
!***********************************************
    SUBROUTINE photocurrent_calculation( pcurrent, pcurrentmax, ene, gamma_L, gamma_R, f_L_e, f_R_e,  f_L_e_minus_omega,  f_R_e_minus_omega, f_L_e_minus_omega0, f_R_e_minus_omega0, photongrid, v0ph, dimx, nphotonfrequency, omegamax, ihomo, ilumo, ehomo, elumo,  gamma_broadening) ! Add HOMO projector and LUMO projector in the future
   !***********************************************
  IMPLICIT NONE
      INTEGER, INTENT(in)   :: dimx
      INTEGER, INTENT(in)   :: nphotonfrequency

      REAL(dbl), INTENT(out) ::  pcurrent(nphotonfrequency)
      REAL(dbl), INTENT(out) ::  pcurrentmax
      !
      COMPLEX(dbl), INTENT(in) :: gamma_L(dimx,dimx)
      COMPLEX(dbl), INTENT(in) :: gamma_R(dimx,dimx)
      REAL(dbl), INTENT(in) ::  ene
      REAL(dbl), INTENT(in) ::  f_L_e
      REAL(dbl), INTENT(in) ::  f_R_e
      REAL(dbl), INTENT(in) ::  f_L_e_minus_omega(nphotonfrequency)
      REAL(dbl), INTENT(in) ::  f_R_e_minus_omega(nphotonfrequency)
      REAL(dbl), INTENT(in) ::  f_L_e_minus_omega0
      REAL(dbl), INTENT(in) ::  f_R_e_minus_omega0
      REAL(dbl), INTENT(in) ::  photongrid(nphotonfrequency)
      REAL(dbl), INTENT(in) ::  v0ph
      REAL(dbl), INTENT(in) ::  omegamax
      REAL(dbl), INTENT(in) ::  gamma_broadening
      INTEGER, INTENT(in)   :: ihomo(dimx)
      INTEGER, INTENT(in)   :: ilumo(dimx)
      REAL(dbl), INTENT(in) ::  ehomo(dimx)
      REAL(dbl), INTENT(in) ::  elumo(dimx)

!      INTEGER, INTENT(in)   :: ihomo(INT(dimx/2))
!      INTEGER, INTENT(in)   :: ilumo(INT(dimx/2))
!      REAL(dbl), INTENT(in) ::  ehomo(INT(dimx/2))
!      REAL(dbl), INTENT(in) ::  elumo(INT(dimx/2))
      ! Local
       CHARACTER(24)      :: subname="photocurrent_calculation"
       INTEGER :: ierr,  iph, ilev1, ilev2
      REAL(dbl) ::  num
      REAL(dbl) ::  denom
      REAL(dbl) ::  golden_rule_broad


        IF( dimx < 1 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in Conductor dimension", " dimC  ", dimx
            STOP
         ENDIF 

         IF( nphotonfrequency < 1 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in Number of photon frequencies definition", " nphotonfrequency  ", nphotonfrequency
            STOP
         ENDIF 

         IF(  v0ph < ZERO ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in Coupling constant definition", " photoncoupling  ", v0ph
            STOP
         ENDIF 


         IF(  v0ph == ZERO ) THEN 
            PRINT*, "-----WARNING----- Routine: ",  TRIM(subname) 
            PRINT*, "Coupling constant == 0"
            PRINT*, "The photocurrent will not be calculated"
            RETURN
         ENDIF 

          pcurrent(:) = ZERO
          pcurrentmax = ZERO


	DO iph = 1, nphotonfrequency
            ! SUm over all the transitions 
            DO ilev1 = 1, dimx
               DO ilev2 =1, dimx  
     !            DO ilev1 = 1, (INT(dimx/2))  
     !               DO ilev2 =1, (INT(dimx/2))  

                 num =ZERO
                 denom =ONE
                 !IF ( (elumo(ilev2) - ehomo(ilev1))  < ( photongrid(iph) + 10.00 )  ) THEN 
		         golden_rule_broad =  gamma_broadening * gamma_broadening / ( (photongrid(iph) - (elumo(ilev2) - ehomo(ilev1)))**(2)  + gamma_broadening**2 )
		          ! put BROADENING OF THE TRANSITION HERE + THRESHOLD
		         num = ( gamma_L(ihomo(ilev1),ihomo(ilev1) ) * gamma_R(ilumo(ilev2),ilumo(ilev2)) * f_L_e_minus_omega(iph) * (ONE -  f_R_e ) ) -   gamma_L(ilumo(ilev2),ilumo(ilev2)) *  gamma_R(ihomo(ilev1),ihomo(ilev1))  * f_R_e_minus_omega(iph) * ( ONE -  f_L_e ) 
		         denom =  ( ( ene - elumo(ilev2))**2 + ( (gamma_R(ilumo(ilev2),ilumo(ilev2)) + gamma_L(ilumo(ilev2),ilumo(ilev2)))/2.00 )**2  ) * ( (ene - photongrid(iph) - ehomo(ilev1))**2 + ( ( gamma_R(ihomo(ilev1),ihomo(ilev1)) + gamma_L(ihomo(ilev1),ihomo(ilev1)) ) /2.00 )**2 )
                 !ENDIF
                 pcurrent(iph) =  pcurrent(iph) + ( v0ph**2 ) * ( num / denom ) * golden_rule_broad
              ENDDO 
           ENDDO 

        ENDDO

        num =ZERO
        denom =ZERO

        num = ( gamma_L(ihomo(1),ihomo(1)) * gamma_R(ilumo(1),ilumo(1)) * f_L_e_minus_omega0 * (ONE -  f_R_e ) ) -   gamma_L(ilumo(1),ilumo(1)) *  gamma_R(ihomo(1),ihomo(1))  * f_R_e_minus_omega0 * (ONE -  f_L_e ) 
        denom =  ( ( ene - elumo(1))**2 + ( (gamma_R(ilumo(1),ilumo(1)) + gamma_L(ilumo(1),ilumo(1)))/2.00 )**2  ) * ( ( ene - omegamax -  ehomo(1))**2 + ( ( gamma_R(ihomo(1),ihomo(1)) + gamma_L(ihomo(1),ihomo(1)) ) /2.00 )**2 )

        pcurrentmax = ( v0ph**2 ) * ( num / denom )

          !

    END SUBROUTINE photocurrent_calculation

!***********************************************
   SUBROUTINE current_calculation(current, transmittance, f_L_e, f_R_e)
   !***********************************************
   ! Calculate the surface green function associated with a 1D semi infinite and periodic chain
  IMPLICIT NONE
      REAL(dbl), INTENT(out) ::  current
      REAL(dbl), INTENT(in) :: transmittance
      REAL(dbl), INTENT(in) :: f_L_e
      REAL(dbl), INTENT(in) :: f_R_e
      
      ! Local
       CHARACTER(19)      :: subname="current_calculation"
       INTEGER :: ierr
       !

       current = ZERO
       current = transmittance * (f_L_e - f_R_e)
       
       
    END SUBROUTINE current_calculation
  END  MODULE transport_formula_module 

 !***********************************************
  ! SUBROUTINE photocurrent_calculation( pcurrent, pcurrentmax, egrid, gamma_L, gamma_R, f_L_e, f_R_e,  f_L_e_minus_omega,  f_R_e_minus_omega, f_L_e_minus_omega0, f_R_e_minus_omega0, photongrid, v0ph, deltaenergy, nenergy, dimx, nphotonfrequency, omegamax, ihomo, ilumo, ehomo, elumo) ! Add HOMO projector and LUMO projector in the future
   !***********************************************
  !IMPLICIT NONE
   !   REAL(dbl), INTENT(inout) ::  pcurrent(nphotonfrequency)
   !   REAL(dbl), INTENT(inout) ::  pcurrentmax
      !
   !   COMPLEX(dbl), INTENT(in) :: gamma_L(nenergy,dimx,dimx)
   !   COMPLEX(dbl), INTENT(in) :: gamma_R(nenergy,dimx,dimx)
   !   REAL(dbl), INTENT(in) ::  egrid(nenergy)
   !   REAL(dbl), INTENT(in) ::  f_L_e(nenergy)
   !   REAL(dbl), INTENT(in) ::  f_R_e(nenergy)
   !   REAL(dbl), INTENT(in) ::  f_L_e_minus_omega(nphotonfrequency,nenergy)
   !   REAL(dbl), INTENT(in) ::  f_R_e_minus_omega(nphotonfrequency,nenergy)
  !    REAL(dbl), INTENT(in) ::  f_L_e_minus_omega0(nenergy)
  !   REAL(dbl), INTENT(in) ::  f_R_e_minus_omega0(nenergy)
  !    REAL(dbl), INTENT(in) ::  photongrid(nphotonfrequency)
  !    REAL(dbl), INTENT(in) ::  v0ph
  !    REAL(dbl), INTENT(in) ::  deltaenergy
  !    INTEGER, INTENT(in)   :: nenergy
  !    INTEGER, INTENT(in)   :: dimx
 !     INTEGER, INTENT(in)   :: nphotonfrequency
      !REAL(dbl), INTENT(in) ::  omegamax
     !INTEGER, INTENT(in)   :: ihomo
     ! INTEGER, INTENT(in)   :: ilumo
     ! REAL(dbl), INTENT(in) ::  ehomo
     ! REAL(dbl), INTENT(in) ::  elumo
      ! Local
     !  CHARACTER(24)      :: subname="photocurrent_calculation"
     !  INTEGER :: ierr, ie, iph
     ! REAL(dbl) ::  num
     ! REAL(dbl) ::  denom

       ! IF( dimx < 1 ) THEN 
         !   PRINT*, "Error in Routine ",  TRIM(subname) 
         !   PRINT*, "Error in Conductor dimension", " dimC  ", dimx
         !   STOP
         !ENDIF 

        !IF( nenergy < 1 ) THEN 
        !    PRINT*, "Error in Routine ",  TRIM(subname) 
        !    PRINT*, "Error in Number of energies definition", " Ne  ", nenergy
        !    STOP
        ! ENDIF 

       ! IF( nphotonfrequency < 1 ) THEN 
       !     PRINT*, "Error in Routine ",  TRIM(subname) 
       !     PRINT*, "Error in Number of photon frequencies definition", " nphotonfrequency  ", nphotonfrequency
       !     STOP
       !  ENDIF 

      !   IF(  v0ph < ZERO ) THEN 
      !      PRINT*, "Error in Routine ",  TRIM(subname) 
      !      PRINT*, "Error in Coupling constant definition", " photoncoupling  ", v0ph
     !       STOP
        ! ENDIF 


       !  IF(  v0ph == ZERO ) THEN 
       !     PRINT*, "-----WARNING----- Routine: ",  TRIM(subname) 
     !       PRINT*, "Coupling constant == 0"
      !      PRINT*, "The photocurrent will not be calculated"
    !        RETURN
   !      ENDIF 

  !        pcurrent(:) = ZERO
 !         pcurrentmax = ZERO


!		DO iph = 1, nphotonfrequency
   !                 num =ZERO
  !                  denom =ZERO
  !                  num = ( gamma_L(ie,ihomo,ihomo) * gamma_R(ie,ilumo,ilumo) * f_L_e_minus_omega(iph,ie) * (ONE -  f_R_e(ie)) ) -   gamma_L(ie,ilumo,ilumo) *  gamma_R(ie,ihomo,ihomo)  * f_R_e_minus_omega(iph,ie) * (ONE -  f_L_e(ie)) 
  !                  denom =  ( (egrid(ie) - elumo)**2 + ( (gamma_R(ie,ilumo,ilumo) + gamma_L(ie,ilumo,ilumo))/2.00 )**2  ) * ( (egrid(ie) - photongrid(iph) - ehomo)**2 + ( ( gamma_R(ie,ihomo,ihomo) + gamma_L(ie,ihomo,ihomo) ) /2.00 )**2 )

 !                   pcurrent(iph) =  pcurrent(iph) + ( ( v0ph**2 ) * ( num / denom ) * (ONE/REAL(ne,dbl)) *de  )
 !               ENDDO
!                num =ZERO
!                denom =ZERO

 !               num = ( gamma_L(ie,ihomo,ihomo) * gamma_R(ie,ilumo,ilumo) * f_L_e_minus_omega0(ie) * (ONE -  f_R_e(ie)) ) -   gamma_L(ie,ilumo,ilumo) *  gamma_R(ie,ihomo,ihomo)  * f_R_e_minus_omega0(ie) * (ONE -  f_L_e(ie)) 
 !               denom =  ( (egrid(ie) - elumo)**2 + ( (gamma_R(ie,ilumo,ilumo) + gamma_L(ie,ilumo,ilumo))/2.00 )**2  ) * ( (egrid(ie) - omegamax -  ehomo)**2 + ( ( gamma_R(ie,ihomo,ihomo) + gamma_L(ie,ihomo,ihomo) ) /2.00 )**2 )


!                pcurrentmax = pcurrentmax + ( ( v0ph**2 ) *  ( num / denom ) * (ONE/REAL(ne,dbl)) *de  )

          !

!    END SUBROUTINE photocurrent_calculation

                !CALL photocurrent_calculation( Photocurrent_aux(:), Photocurrentmax_aux, egrid(1:ne), gamma_L(1:ne,:,:), gamma_R(1:ne,:,:), f_L(1:ne), f_R(1:ne),  f_shifted_L(:,1:ne), f_shifted_R(:,1:ne), f_shifted_max_L(1:ne), f_shifted_max_R(1:ne), photonfrequency(:), photoncouplingconstant, de, ne, dimC, nphotonfrequency, omegamax(ialpha, ibias), ihomo, ilumo, ehomo, elumo)
                    !CALL calcul_emittedlight(Emittedlight(ialpha,ibias), gamma_L(:,:,:), gamma_R(:,:,:), fL(:), f_R(:), nomega, ne, dimC)         
 


