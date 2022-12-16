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
   REAL(dbl) :: alpha=ZERO
   REAL(dbl) :: gamma_big=ZERO

   REAL(dbl) :: biasmin=ZERO
   REAL(dbl) :: biasmax=ZERO
   INTEGER ::   nbias=0
   REAL(dbl) :: gapmin=ZERO
   REAL(dbl) :: gapmax=ZERO
   INTEGER ::   ngap=1
   REAL(dbl) :: gammamin=ZERO
   REAL(dbl) :: gammamax=ZERO
   INTEGER :: ngamma=1

   ! Other Module Variables
   CHARACTER(14) :: modulename="T_input_module"
   ! Public variables

   ! Public routines:
   PUBLIC                 :: input_manager
   ! private routines
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
   PRINT*, "     ... READING  alpha"
         READ(100,*) chr,  alpha
           IF ( TRIM(chr) /= "alpha" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for  alpha", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "          alpha=",  alpha

   PRINT*, "     ... READING  gamma_big"
         READ(100,*) chr,  gamma_big
           IF ( TRIM(chr) /= "gamma_big" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for  gamma_big", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "          gamma_big=",  gamma_big

   !
   PRINT*, "     ... READING  biasmin"
         READ(100,*) chr,  biasmin
           IF ( TRIM(chr) /= "biasmin" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for  biasmin", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "          biasmin=",  biasmin
   !
   PRINT*, "     ... READING biasmax"
         READ(100,*) chr, biasmax
           IF ( TRIM(chr) /= "biasmax" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for biasmax", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "         biasmax=", biasmax
   !
   PRINT*, "     ... READING nbias"
         READ(100,*) chr, nbias
           IF ( TRIM(chr) /= "nbias" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for nbias", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "          nemax=", nbias
   !
 !

   PRINT*, "     ... READING gapmin"
         READ(100,*) chr, gapmin
           IF ( TRIM(chr) /= "gapmin" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for gapmin", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "         gapmin=", gapmin
   !
   PRINT*, "     ... READING gapmax"
         READ(100,*) chr, gapmax
           IF ( TRIM(chr) /= "gapmax" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for gapmax", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "         gapmax=", gapmax
   !
   PRINT*, "     ... READING ngap"
         READ(100,"(a4,I10)") chr, ngap
           IF ( TRIM(chr) /= "ngap" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for ngap", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "          ngap=", ngap
   !
   PRINT*, "     ... READING gammamin"
         READ(100,*) chr, gammamin
           IF ( TRIM(chr) /= "gammamin" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for gammamin", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "         gammamin=", gammamin
   !
   PRINT*, "     ... READING gammamax"
         READ(100,*) chr, gammamax
           IF ( TRIM(chr) /= "gammamax" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for gammamax", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "         gammamax=", gammamax
   !
   PRINT*, "     ... READING ngamma"
         READ(100,*) chr, ngamma
           IF ( TRIM(chr) /= "ngamma" ) THEN
            PRINT*, "Error in Routine ",  TRIM(subname),  " of ", TRIM(modulename)
            PRINT*, "Waiting for ngamma", " Input = ", TRIM(chr)
            STOP
         ENDIF 
   PRINT*, "          ngamma=", ngamma
 
   CLOSE(100)
   PRINT*, "... END of READING Input File"

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


        CALL set_biasgrid()
        CALL set_gap()
        CALL set_gamma()
  
        CALL input_summarize()

    END SUBROUTINE input_manager
!***********************************************
   SUBROUTINE  set_control  
   !***********************************************
   ! 
  USE T_control_module,       ONLY :    debug_level_ => debug_level, &
					gamma_big_ => gamma_big, &
					alpha_ => alpha
					!gamma_small_ => gamma_small


  IMPLICIT NONE
      ! Local
       CHARACTER(11)      :: subname="set_control"
       INTEGER :: ierr
         
      
   
        IF ( ( debug_level < 0 ) .OR. ( debug_level > 5 )) THEN 
            PRINT*, "Error in Routine ",  TRIM(subname) , " in module ", TRIM(modulename)
            PRINT*, "Error in debug level definition", " debug level is intended to be 6 >  dl > -1 ", debug_level
            STOP
        ENDIF 

        gamma_big_ = gamma_big
        alpha_ = alpha
        debug_level_ = debug_level

    END SUBROUTINE set_control



!***********************************************
   SUBROUTINE   set_biasgrid
   !***********************************************
   !
  USE T_biasgrid_module,         ONLY : biasmin_         =>  biasmin, &
                                     biasmax_         =>  biasmax, &
                                     nbias_           =>  nbias
   IMPLICIT NONE
      ! Local
       CHARACTER(9)      :: subname="set_egrid"
       INTEGER :: ierr
         
        !      
 
        biasmin_    =  biasmin
        biasmax_    =  biasmax
        nbias_   =  nbias

    END SUBROUTINE set_biasgrid
!***********************************************
   SUBROUTINE   set_gap
   !***********************************************
   ! set_gap    USE T_gapgrid_module,                  ONLY :  ngap, gapmax, gapmin
  USE T_gapgrid_module,         ONLY : gapmin_      =>  gapmin, &
                                        gapmax_      =>   gapmax, &
                                        ngap_        =>  ngap
   IMPLICIT NONE
      ! Local
       CHARACTER(8)      :: subname="set_gap"
       INTEGER :: ierr
         
        !      
         !
        gapmin_   =  gapmin
        gapmax_   =  gapmax
        ngap_     =  ngap

    END SUBROUTINE set_gap
!***********************************************
   SUBROUTINE   set_gamma
   !***********************************************
    ! set_gamma   USE T_gammagrid_module,                 ONLY : ngamma, gammamax, gammamin  
  USE T_gammagrid_module,         ONLY : ngamma_          =>  ngamma, &
                                        gammamax_        =>  gammamax, &
                                        gammamin_        =>  gammamin
   IMPLICIT NONE
      ! Local
       CHARACTER(9)      :: subname="set_gamma"
       INTEGER :: ierr
         
         !
        ngamma_   =  ngamma
        gammamax_   =  gammamax
        gammamin_   =  gammamin
        !
    END SUBROUTINE set_gamma

!
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
       PRINT*, "%%% alpha        ",   alpha
       PRINT*, "%%% gamma_big        ",   gamma_big
       PRINT*, "%%% nbias                              --- Number of points used for the bias grid ",   nbias
       PRINT*, "%%% biasmax                            --- Maximum value of the bias grid",   biasmax
       PRINT*, "%%% biasmin                            --- Minimum value of the bias grid",   biasmin
       PRINT*, "%%% ngap                               --- Number of points used for the gap grid ",   ngap
       PRINT*, "%%% gapmax                             --- Maximum value of the gap grid",   gapmax
       PRINT*, "%%% gapmin                             --- Minimum value of the gap grid",   gapmin
       PRINT*, "%%% ngamma                              --- Number of points used for the gamma grid ",   ngamma
       PRINT*, "%%% gammamax                            --- Maximum value of the gamma grid",   gammamax
       PRINT*, "%%% gammamin                            --- Minimum value of the gamma grid",   gammamin
 
  
    END SUBROUTINE input_summarize  
  END  MODULE T_input_module

 

 


