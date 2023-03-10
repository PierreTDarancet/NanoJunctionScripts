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
   MODULE self_energy_formula_module
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
   REAL(dbl), PARAMETER :: delta_lead =  EPS_m3
   ! Public variables

   ! Public routines:
   PUBLIC                 ::  lead_green_function ! Not used
   PUBLIC                 ::  calcul_gamma ! Not Used
   PUBLIC                 ::  calcul_electron_photon_sigma
   PUBLIC                 ::  calcul_sigma_ee_retarded

   CONTAINS 
!***********************************************
   SUBROUTINE lead_green_function(  g_lead, h00, h01, ene, dimx)
   !***********************************************
   ! Calculate the surface green function associated with a 1D semi infinite and periodic chain
  IMPLICIT NONE
      INTEGER, INTENT(in)   :: dimx

      COMPLEX(dbl), INTENT(out) ::  g_lead(dimx,dimx) 
      COMPLEX(dbl), INTENT(in) ::  h00(dimx,dimx) 
      COMPLEX(dbl), INTENT(in) ::  h01(dimx,dimx) 
      REAL(dbl), INTENT(in) ::  ene 
      ! Local
       CHARACTER(19)      :: subname="lead_green_function"
       INTEGER :: ierr, ie
       REAL(dbl)              :: spectra_min, spectra_max
       COMPLEX(dbl)           :: delta_spectra     ! discriminant
       COMPLEX(dbl)           :: enework, work_scal

        IF( dimx < 1 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in Lead dimension", " dim Lead = ", dimx
            STOP
         ENDIF 


        g_lead(:,:) = CZERO

        IF (dimx == 1) THEN 
           enework = ene + (delta_lead*CI)
           spectra_min= h00(1,1) - 2 * ABS (h01(1,1))
           spectra_max= h00(1,1) + 2 * ABS (h01(1,1))
           !
           
           IF ( ( ene <= spectra_max ) .AND. ( ene >=  spectra_min )) THEN
                ! b?? - 4ac
                delta_spectra = (h00(1,1) - enework )**2 - (4 * (h01(1,1))**2)
                IF ( REAL(delta_spectra) > ZERO )  THEN 
		    PRINT*, "Error in Routine ",  TRIM(subname) 
		    PRINT*, "Positive discriminant in the spectra", " Energy= ", ene, " Discriminant", REAL(delta_spectra)
		    STOP
                ENDIF
                work_scal = (  -(h00(1,1)- ene) - CI* SQRT ( -delta_spectra ) ) / (2* (h01(1,1)**2))
                g_lead(1,1) = CONE / (ene - h00(1,1) - ((h01(1,1) **2) * work_scal))
            ELSE 
                ! OUtside spectra
                g_lead(1,1) = CZERO
                !
            ENDIF
        !
        ELSE IF( dimx > 1 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Method for dealing with dim > 1 Leads not yet implemented", " dim_Lead = ", dimx
            STOP
       !
        ELSE
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Problem with Lead dimension", " dim_Lead = ", dimx
            STOP
         ENDIF 
       !
    END SUBROUTINE lead_green_function
!***********************************************
   SUBROUTINE calcul_gamma( gamma_lead, sigma_lead ,  g_lead, h_lead_C, h_C_lead, dimx, dimC )
   !***********************************************
  IMPLICIT NONE
      INTEGER, INTENT(in)   :: dimx
      INTEGER, INTENT(in)   :: dimC
      COMPLEX(dbl), INTENT(out) ::  gamma_lead(dimC,dimC) 
      COMPLEX(dbl), INTENT(out) ::  sigma_lead(dimC,dimC) 
      COMPLEX(dbl), INTENT(in) ::  g_lead(dimx,dimx) 
      COMPLEX(dbl), INTENT(in) ::  h_lead_C(dimx,dimC) 
      COMPLEX(dbl), INTENT(in) ::  h_C_lead(dimC,dimx) 
      ! Local
       CHARACTER(12)      :: subname="calcul_gamma"
       INTEGER :: ierr
     COMPLEX(dbl)  ::  work(dimC,dimx) 
 
      IF( dimx < 1 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in Lead dimension", " dim Lead = ", dimx
            STOP
         ENDIF 
      IF( dimC < 1 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in Conductor dimension", " dim C = ", dimC
            STOP
         ENDIF 

          !
          sigma_lead(:,:) = CZERO
          ! 
          gamma_lead(:,:) = CZERO
          !
          CALL mat_mul(work, h_C_lead , 'N', g_lead, 'N', dimC, dimx, dimx)
          CALL mat_mul(sigma_lead, work, 'N', h_lead_C, 'N', dimC, dimC, dimx) 
          !
          ! gamma_L and gamma_R
          !
          gamma_lead(:,:) = CI * (  sigma_lead(:,:) - CONJG( TRANSPOSE(sigma_lead(:,:)) )   )

    END SUBROUTINE calcul_gamma

!***********************************************
   SUBROUTINE calcul_electron_photon_sigma(sigma_r, sigma_lesser, sigma_greater, g_r_eneplusomega, glesser_eneplusomega, ggreater_eneplusomega, g_r_eneminusomega, glesser_eneminusomega, ggreater_eneminusomega, amplitude, omegaphoton, dipole_matrix, dimC, method)
   !***********************************************
   !   CALL calcul_electron_photon_sigma(sigma_ep_r(ie,:,:), sigma_ep_lesser(ie,:,:), sigma_ep_greater(ie,:,:),  g_C(egridplusomega(ie),:,:), glesser(egridplusomega(ie),:,:), ggreater(egridplusomega(ie),:,:), g_C(egridminusomega(ie),:,:), glesser(egridminusomega(ie),:,:), ggreater(egridminusomega(ie),:,:), amplitude, photonfrequency, dipole_matrix(:,:), dimC, ep_method)   

  IMPLICIT NONE
      INTEGER, INTENT(in)   :: dimC
      COMPLEX(dbl), INTENT(out) ::  sigma_r(dimC,dimC) 
      COMPLEX(dbl), INTENT(out) ::  sigma_lesser(dimC,dimC) 
      COMPLEX(dbl), INTENT(out) ::  sigma_greater(dimC,dimC) 
      COMPLEX(dbl), INTENT(in) ::  g_r_eneplusomega(dimC,dimC)
      COMPLEX(dbl), INTENT(in) ::  glesser_eneplusomega(dimC,dimC)
      COMPLEX(dbl), INTENT(in) ::  ggreater_eneplusomega(dimC,dimC)
      COMPLEX(dbl), INTENT(in) ::  g_r_eneminusomega(dimC,dimC)
      COMPLEX(dbl), INTENT(in) ::  glesser_eneminusomega(dimC,dimC)
      COMPLEX(dbl), INTENT(in) ::  ggreater_eneminusomega(dimC,dimC)
      REAL(dbl), INTENT(in) ::  amplitude
      REAL(dbl), INTENT(in) ::  omegaphoton
      REAL(dbl), INTENT(in) ::  dipole_matrix(dimC,dimC)
      CHARACTER(2), INTENT(in)   :: method
      ! Local
       CHARACTER(28)      :: subname="calcul_electron_photon_sigma"
       INTEGER :: ierr
     COMPLEX(dbl)  ::  work1(dimC,dimC) 
     COMPLEX(dbl)  ::  work2(dimC,dimC) 
     COMPLEX(dbl)  ::  work3(dimC,dimC) 
     COMPLEX(dbl)  ::  work4(dimC,dimC) 
     COMPLEX(dbl)  ::  work_dm(dimC,dimC) 


 
         IF( dimC < 1 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) 
            PRINT*, "Error in Conductor dimension", " dim C = ", dimC
            STOP
         ENDIF 

             work_dm(:,:) = dipole_matrix(:,:)
          ! sigma_lesser

          IF (( TRIM(method)=="AM") .OR. ((TRIM(method)=="GN"))) THEN ! Aeberhard Morf PRB 77 (2008) - SCBA
              work4(:,:) =  glesser_eneplusomega(:,:) +  glesser_eneminusomega(:,:)  
              CALL mat_mul(work3, work4 , 'N', work_dm, 'N', dimC, dimC, dimC)
              CALL mat_mul(work1, work_dm, 'N', work3, 'N', dimC, dimC, dimC)
              !sigma_lesser(:,:) = CI*work1(:,:)* (amplitude**(2))
              sigma_lesser(:,:) = work1(:,:)* (amplitude**(2))
              work4(:,:) =  ggreater_eneplusomega(:,:) +  ggreater_eneminusomega(:,:)  
              CALL mat_mul(work3, work4 , 'N', work_dm, 'N', dimC, dimC, dimC)
              CALL mat_mul(work1, work_dm, 'N', work3, 'N', dimC, dimC, dimC)
              ! sigma_greater(:,:) = CI*work1(:,:)* (amplitude**(2))
                sigma_greater(:,:) = work1(:,:)* (amplitude**(2))
              !sigma_lesser(:,:) = work1(:,:)* (amplitude**(2))
          ELSE IF (TRIM(method)=="G0" ) THEN ! Galperin Nitzan PRL 2005 with dipole transition/ Galperin Nitzan PRL 2005 without dipole transition -set to one
              work4(:,:) =  glesser_eneplusomega(:,:) +  glesser_eneminusomega(:,:)  
              CALL mat_mul(work3, work4 , 'N', work_dm, 'N', dimC, dimC, dimC)
              CALL mat_mul(work1, work_dm, 'N', work3, 'N', dimC, dimC, dimC)
              !sigma_lesser(:,:) = CI*work1(:,:)* (amplitude**(2))
              sigma_lesser(:,:) = work1(:,:)* (amplitude**(2))
              work4(:,:) =  ggreater_eneplusomega(:,:) +  ggreater_eneminusomega(:,:)  
              CALL mat_mul(work3, work4 , 'N', work_dm, 'N', dimC, dimC, dimC)
              CALL mat_mul(work1, work_dm, 'N', work3, 'N', dimC, dimC, dimC)
               !sigma_greater(:,:) = CI*work1(:,:)* (amplitude**(2))
               sigma_greater(:,:) = work1(:,:)* (amplitude**(2))

          ELSE
              sigma_lesser(:,:)=CZERO
          ENDIF 



          ! sigma_retarded

          IF (TRIM(method)=="AM") THEN ! Aeberhard Morf PRB 77 (2008) - SCBA
              work4(:,:) =  g_r_eneplusomega(:,:) +  g_r_eneminusomega(:,:)  
              work3(:,:) = (ONE/(ONE*2.0)) * (  glesser_eneminusomega(:,:) - glesser_eneplusomega(:,:)  )
              work2(:,:) = work4(:,:) + work3(:,:)
              CALL mat_mul(work3, work2 , 'N', work_dm, 'N', dimC, dimC, dimC)
              CALL mat_mul(work1, work_dm, 'N', work3, 'N', dimC, dimC, dimC)
              sigma_r(:,:) = CI*work1(:,:)* (amplitude**(2))
          ELSE IF ((TRIM(method)=="GN").OR.(TRIM(method)=="G0")) THEN ! Galperin Nitzan PRL 2005
              sigma_r(:,:)=CZERO
          ELSE
              sigma_r(:,:)=CZERO
          ENDIF 

    END SUBROUTINE calcul_electron_photon_sigma

!***********************************************
   SUBROUTINE calcul_sigma_ee_retarded(  glesser, A_C, sigma_ee_mol_corr, fermi_energy, electronic_temperature, egrid, dimC, ne,  method )
   !***********************************************
   ! CALL calcul_sigma_ee_retarded( glesser(:,:,:),  A_C(:,:,:), sigma_ee_mol_corr(:), fermi_energy,  electronic_temperature, green_egrid(:) ,  dimC, ne_green, ee_method)
  IMPLICIT NONE
      INTEGER, INTENT(in)   :: dimC
      COMPLEX(dbl), INTENT(out) :: sigma_ee_mol_corr(dimC)
      COMPLEX(dbl), INTENT(in)  :: glesser(ne,dimC,dimC)
      REAL(dbl), INTENT(in)     :: A_C(ne,dimC,dimC)
      REAL(dbl), INTENT(in)     :: egrid(ne)
      REAL(dbl), INTENT(in)     :: fermi_energy
      REAL(dbl), INTENT(in)     :: electronic_temperature
      CHARACTER(3), INTENT(in)  :: method

      INTEGER, INTENT(in)   :: ne
      INTEGER   :: i_dim, ie

      REAL, PARAMETER:: Exciton_binding=2.2*ONE

          sigma_ee_mol_corr(:)= CZERO
          IF (TRIM(method)=="EXC") THEN ! INject a part of excitonic effects accounting for the coulomb hole
            DO i_dim = 1, dimC
            sigma_ee_mol_corr(i_dim)=CZERO
                DO ie=1, ne
                    sigma_ee_mol_corr(i_dim)= sigma_ee_mol_corr(i_dim) + SIGN(fermi_energy-egrid(ie),ONE) *  Exciton_binding*( ABS(  AIMAG( glesser(ie,i_dim,i_dim) ) ) - ONE/(EXP((egrid(ie)-fermi_energy)/electronic_temperature)+ ONE) * A_C(ie,i_dim,i_dim) )
                ENDDO
            ENDDO
          ELSE
            RETURN
          ENDIF 


    END SUBROUTINE calcul_sigma_ee_retarded
  END  MODULE self_energy_formula_module

 

 


