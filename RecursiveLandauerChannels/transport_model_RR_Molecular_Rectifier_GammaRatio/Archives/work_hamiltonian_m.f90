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


    ! define leads

   h00_L(1,1)=CZERO
   h00_R(1,1)=CZERO


   h01_L(1,1)=4.000*CONE
   h01_R(1,1)=4.000*CONE
     HOMO_right = HOMO_left+vocgrid(ivoc)             ! indicated wrt to \mu_R
     LUMO_right = LUMO_left+vocgrid(ivoc)             ! Let's take Right as the Donor

     PRINT*, "Voc", vocgrid(ivoc)            


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
  ! define unity matrix
   s00_C(:,:) = CZERO
 
   DO i=1, Nb_States
      s00_C(i,i) = CONE
   ENDDO 


  DEALLOCATE ( gamma_R, STAT=ierr )
  DEALLOCATE ( sigma_R , STAT=ierr )
  DEALLOCATE ( gamma_L , STAT=ierr )
  DEALLOCATE ( sigma_L , STAT=ierr )

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


