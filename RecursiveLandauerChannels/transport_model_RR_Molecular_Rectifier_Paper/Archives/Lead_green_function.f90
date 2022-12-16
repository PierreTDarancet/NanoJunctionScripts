
   REAL(dbl)              :: spectra_min, spectra_max
   REAL(dbl)              :: delta_spectra     ! discriminant

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




           gL(:,:) = CZERO
           spectra_min= aux00_L(1,1)  - 2 * ABS (aux01_L(1,1))
           spectra_max= aux00_L(1,1)  + 2 * ABS (aux01_L(1,1))



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


