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
   MODULE T_input_module
   !***********************************************
  USE kinds, ONLY : dbl
  USE constants, ONLY : ZERO, ONE, CZERO, CONE, CI,PI, EPS_m3, EPS_m5
   IMPLICIT NONE
   PRIVATE 
   SAVE
!
   ! Public
   ! Private
   ! Default Values
   INTEGER :: debug_level=0
   INTEGER :: nmaxiter_gf=0
   REAL(dbl) :: amplitude=EPS_m5
   REAL(dbl) :: photonfrequency
   CHARACTER(2) :: ep_method
   CHARACTER(3) :: ee_method
   INTEGER ::  ne=101
   REAL(dbl) :: fermi_energy=ZERO
   REAL(dbl) :: emin
   REAL(dbl) :: emax
   REAL(dbl) :: electronic_temperature=EPS_m5
   INTEGER :: nk=0
   INTEGER :: dimC=0


   CHARACTER(14) :: modulename="T_input_module"

   ! Public variables

   ! Public routines:
   PUBLIC                 :: input_manager
   ! private routines
   ! read_input
   ! set_control  USE T_control_module,                  ONLY : nprint, nprint_PDOS, nprint_transmission, calculate_RR, debug_level
   ! set_kgrid
   ! set_egrid   USE T_egrid_module,                     ONLY : emin, emax, nemax
   ! set_hamiltonian_workspace     USE hamiltonian_workspace_module,       ONLY : dimL, dimC, dimR
   ! input_summarize

  CONTAINS 
!***********************************************
   SUBROUTINE read_input
   !***********************************************
   ! read the inp file located in the local directory
  IMPLICIT NONE
      ! Local
       CHARACTER(10)      :: subname="read_input"
       INTEGER :: ierr
        CHARACTER(100)  :: chr
       CHARACTER(3)      :: fname_in="inp"


   PRINT*, "... READING Input File"
   OPEN(100, FILE=TRIM(fname_in),FORM='formatted')
   !
   PRINT*, "     ... READING debug_level"
         READ(100,"(a11,I10)") chr, debug_level
           IF ( TRIM(chr) /= "debug_level" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for debug_level", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "          debug_level=", debug_level
   !
   PRINT*, "     ... READING nmaxiter_gf"
         READ(100,"(a11,I10)") chr, nmaxiter_gf
           IF ( TRIM(chr) /= "nmaxiter_gf" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for nmaxiter_gf", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "          nmaxiter_gf=", nmaxiter_gf
   !

    PRINT*, "     ... READING  amplitude"
         READ(100,"(a9,F10.10)") chr, amplitude
           IF ( TRIM(chr) /= "amplitude" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for amplitude", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "          amplitude=",  amplitude
   !
    PRINT*, "     ... READING  photonfrequency"
         READ(100,"(a15,F10.10)") chr, photonfrequency
           IF ( TRIM(chr) /= "photonfrequency" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for photonfrequency", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "          photonfrequency=",  photonfrequency


    PRINT*, "     ... READING   ep_method"
         READ(100,"(a10,a2)") chr, ep_method
           IF ( TRIM(chr) /= "ep_method" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for ep_method", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "          ep_method=", ep_method

    PRINT*, "     ... READING   ee_method"
         READ(100,"(a10,a3)") chr, ee_method
           IF ( TRIM(chr) /= "ee_method" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for ee_method", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "          ee_method=", ee_method

   !
   PRINT*, "     ... READING ne"
         READ(100,"(a2,I10)") chr, ne
           IF ( TRIM(chr) /= "ne" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for ne", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "          ne=", ne
   !
   PRINT*, "     ... READING fermi_energy"
         READ(100,"(a12,F10.10)") chr,  fermi_energy
           IF ( TRIM(chr) /= "fermi_energy" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for fermi_energy", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "          fermi_energy=",  fermi_energy
   
 
   PRINT*, "     ... READING emin"
         READ(100,"(a4,F10.10)") chr, emin
           IF ( TRIM(chr) /= "emin" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for emin", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "         emin=", emin
   !
   PRINT*, "     ... READING emax"
         READ(100,"(a4,F10.10)") chr, emax
           IF ( TRIM(chr) /= "emax" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for emax", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "         emax=", emax
!
 !
   PRINT*, "     ... READING  electronic_temperature"
         READ(100,"(a22,F10.10)") chr,   electronic_temperature
           IF ( TRIM(chr) /= "electronic_temperature" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for electronic_temperature", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "          electronic_temperature=",  electronic_temperature
 
   PRINT*, "     ... READING nk"
         READ(100,"(a2,I10)") chr, nk
           IF ( TRIM(chr) /= "nk" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for nk", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "          nk=", nk



   PRINT*, "     ... READING dimC"
         READ(100,"(a4,I10)") chr, dimC
           IF ( TRIM(chr) /= "dimC" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for dimC", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "          dimC=", dimC

 
   PRINT*, "... END of READING Input File"
   CLOSE(100)

    END SUBROUTINE read_input
!***********************************************
   SUBROUTINE input_manager
   !***********************************************
   ! 
  IMPLICIT NONE
      ! Local
       CHARACTER(13)      :: subname="input_manager"
       INTEGER :: ierr
 
        CALL read_input()


        CALL set_control()
        CALL set_kgrid()

        CALL set_egrid()


        CALL set_hamiltonian_workspace()


        CALL input_summarize()

    END SUBROUTINE input_manager
!***********************************************
   SUBROUTINE  set_control  
   !***********************************************
   ! 
    USE T_control_module,       ONLY : debug_level_ => debug_level, &
                                     nmaxiter_gf_ => nmaxiter_gf, &
                                     amplitude_ => amplitude, &
                                     photonfrequency_ => photonfrequency, &
                                     ee_method_ => ee_method, &
                                     ep_method_ => ep_method 
  IMPLICIT NONE
      ! Local
       CHARACTER(11)      :: subname="set_control"
       INTEGER :: ierr
         
      
  
        IF ( ( debug_level < 0 ) .OR. ( debug_level > 5 )) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) , " in module ", TRIM(modulename)
            PRINT*, "Error in debug level definition", " debug level is intended to be 6 >  dl > -1 ", debug_level
            STOP
        ENDIF 

        debug_level_ = debug_level
        nmaxiter_gf_ = nmaxiter_gf
        amplitude_ = amplitude
        photonfrequency_ = photonfrequency
        ee_method_ = ee_method
        ep_method_ = ep_method    

    END SUBROUTINE set_control
!***********************************************
   SUBROUTINE   set_kgrid
   !***********************************************
   ! 
  USE kpoint_module,          ONLY : nk_  =>  nk
   IMPLICIT NONE
      ! Local
       CHARACTER(10)      :: subname="set_kgrid"
       INTEGER :: ierr
         
      
        IF( nk < 1 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) , " in module ", TRIM(modulename)
            PRINT*, "Error in nk:", "nk is intended to be > = 1", nk
            STOP
        ENDIF 
     
        nk_ = nk


    END SUBROUTINE set_kgrid




!***********************************************
   SUBROUTINE   set_egrid
   !***********************************************
   !
  USE T_egrid_module,                     ONLY : ne_           => ne, &
                                                 fermi_energy_ => fermi_energy, &
                                                 emin_         =>  emin, &
                                                 emax_         =>  emax, &
                                                 electronic_temperature_ => electronic_temperature
   IMPLICIT NONE
      ! Local
       CHARACTER(9)      :: subname="set_egrid"
       INTEGER :: ierr
         
        !      
        IF( ne < 1 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) , " in module ", TRIM(modulename)
            PRINT*, "Error in nemax definition:", "nemax is intended to be > = 1", ne
            STOP
        ENDIF 
        !
        IF (   emax < emin  ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) , " in module ", TRIM(modulename)
            PRINT*, "Error in energy definition:", " emax is supposed > emin "
            PRINT*, " emax=", emax
            PRINT*, " emin=", emin
            STOP
        ENDIF 
        !
        IF (  ( emax == emin ) .AND. ( ne > 1 )  ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) , " in module ", TRIM(modulename)
            PRINT*, "Error in energy definition:", " emax == emin while nemax > 1"
            PRINT*, " emax=", emax
            PRINT*, " emin=", emin
            PRINT*, " nemax", ne
            STOP
        ENDIF 

        emin_    =  emin
        emax_    =  emax
        ne_   =  ne
        fermi_energy_ = fermi_energy
        electronic_temperature_ = electronic_temperature
    END SUBROUTINE set_egrid

!***********************************************
   SUBROUTINE   set_hamiltonian_workspace
   !***********************************************
   ! set_hamiltonian_workspace     USE hamiltonian_workspace_module,       ONLY : dimL, dimC, dimR
  USE hamiltonian_workspace_module,         ONLY : dimC_          =>  dimC

   IMPLICIT NONE
      ! Local
       CHARACTER(25)      :: subname="set_hamiltonian_workspace"
       INTEGER :: ierr
         
        !      
        IF( dimC < 1 ) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) , " in module ", TRIM(modulename)
            PRINT*, "Error in dimC definition:", " dimC is intended to be > = 1", dimC
            STOP
        ENDIF 
        !

        dimC_   =  dimC
        !
    END SUBROUTINE set_hamiltonian_workspace


!***********************************************
   SUBROUTINE  input_summarize  
   !***********************************************
   ! 
  IMPLICIT NONE
      ! Local
       CHARACTER(15)      :: subname="input_summarize"
       INTEGER :: ierr

       PRINT*, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
       PRINT*, "%                      General Control Variables                              %"
       PRINT*, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
       PRINT*, "%%% debug_level  ---  debug level set to ",  debug_level
   IF (debug_level == 0) THEN
       PRINT*, "%%% Debug Mode: OFF"
   ELSE
       PRINT*, "%%% Debug Mode:  ON - Not Yet fully Implemented"
       PRINT*, "%%% Debug Mode>= 1  - All the standard output files are printed"
       PRINT*, "%%% Debug Mode>= 2  - All output transmission files are printed"
       PRINT*, "%%% Debug Mode>= 3  - All the Grids                 are printed"
   ENDIF
       PRINT*, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
       PRINT*, "%                            Method Variables                                 %"
       PRINT*, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
       PRINT*, "%%% nmaxiter_gf  ---  Iterative calculation of SE ",  nmaxiter_gf
       PRINT*, "%%% ep_method    ---  Electron Photon SE formula  ",  ep_method
   IF (TRIM(ep_method) == "G0") THEN
       PRINT*, "%%% ep_method = Galperin Nitzan with no dipole transition matrix"
   ELSE IF (TRIM(ep_method) == "GN") THEN
       PRINT*, "%%% ep_method = Galperin Nitzan with dipole transition matrix"
   ELSE IF (TRIM(ep_method) == "AM") THEN
       PRINT*, "%%% ep_method = Aeberhard Morf - SCBA SE for Electron-Photon SE"
   ELSE
       PRINT*, "%%% ep_method = None - Electron -photon SE will be set to 0"
   ENDIF
       PRINT*, "%%% ee_method    ---  Electron Electron SE formula",  ee_method
   IF (TRIM(ee_method) == "DFT") THEN
       PRINT*, "%%% ee_method = No Electron-electron correlations: DFT calculation"
   ELSE IF (TRIM(ee_method) == "SIG") THEN
       PRINT*, "%%% ee_method = Neaton method DFT+Sigma"
   ELSE IF (TRIM(ee_method) == "EXC") THEN
       PRINT*, "%%% ee_method = DFT+SIGMA+EXC Approximative inclusion of"
       PRINT*, "%%%             Excitonic effects on the spectral function"
       PRINT*, "%%% WARNING     Not fully implemented yet-need glesser"
   ELSE
       PRINT*, "%%% ee_method = Bad Value of the Electron-Electron SE"
       STOP
   ENDIF

       PRINT*, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
       PRINT*, "%                     Electronic Structure variables                          %"
       PRINT*, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
       PRINT*, "%%% fermi_energy           --- Fermi Energy           ",   fermi_energy
       PRINT*, "%%% electronic_temperature --- Temperature used in FD ", electronic_temperature
       PRINT*, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
       PRINT*, "%                          Photon coupling variables                          %"
       PRINT*, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
       PRINT*, "%%% amplitude       --- Coupling Constant   ",   amplitude
       PRINT*, "%%% photonfrequency --- Energy of the photon",  photonfrequency

       PRINT*, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
       PRINT*, "%                              Dimensions                                     %"
       PRINT*, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
       PRINT*, "%%% dimC   --- Number of molecular levels  --- ",   dimC

   IF ( dimC == 1) THEN
       PRINT*, "%%% WARNING %%%%%%               --- Only one Molecular level used"
   ENDIF
       PRINT*, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
       PRINT*, "%                        Grid Control Variables                               %"
       PRINT*, "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
       PRINT*, "%                 Energy grid Control Variables "
       PRINT*, "%%% ne       --- Maximum Number of points used for the Transmission energy grid ",   ne
       PRINT*, "%%% emax     --- Maximum value of Transmission the energy grid",   emax
       PRINT*, "%%% emin     --- Minimum value of Transmission the energy grid",   emin
        PRINT*, "%%% nk       --- Number of k points used for the Transmission",   nk
   
    END SUBROUTINE input_summarize  
  END  MODULE T_input_module

 

 


