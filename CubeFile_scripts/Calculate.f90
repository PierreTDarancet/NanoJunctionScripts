!
! Copyright (C) 2011 The Molecular Foundry Berkeley
!
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors  : Pierre DARANCET
!***********************************************
   MODULE calculatequantities_Module
   !***********************************************
   IMPLICIT NONE
   PRIVATE 
   SAVE


   ! Public routines:
   PUBLIC                 ::  CalculateVolume
   PUBLIC                 ::  CalculateCenterOfMass
   PUBLIC                 ::  CalculateCenterOfCharge
   PUBLIC                 ::  CalculateDipoleOperator
   PUBLIC                 ::  CalculateDipoleOperatorWF
   PUBLIC                 ::  CalculateFirstOrderEnergy

 
CONTAINS

!***********************************************
      SUBROUTINE CalculateFirstOrderEnergy(ExpectationValue, Integrated_volume, Integrated_density, dimensions,  X0Cell, dV, cell, mesh, Wavefunction, field)
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE
   ! Calculte < \Psi |V | \Psi >
   REAL(dbl),   INTENT(out)        ::  Integrated_density
   REAL(dbl),   INTENT(out)        ::  Integrated_volume
   REAL(dbl),   INTENT(out)        ::  ExpectationValue
   INTEGER,     INTENT(in)         ::  dimensions
   INTEGER,     INTENT(in)         ::  mesh(dimensions)
   REAL(dbl),   INTENT(in)         ::  X0Cell(dimensions)
   REAL(dbl),   INTENT(in)         ::  cell(dimensions,dimensions)
   REAL(dbl),   INTENT(in)         ::  dV
   REAL(dbl),   INTENT(in)         ::   Wavefunction( mesh(1), mesh(2), mesh(3) )
   REAL(dbl),   INTENT(in)         ::   field( mesh(1), mesh(2), mesh(3) )

   INTEGER :: ia, ib, ic, idimension    

   PRINT*, "    ...Calculating First Order energy correction"
   ExpectationValue = ZERO
   Integrated_density= ZERO
   Integrated_volume   = ZERO


   DO ic=1,mesh(3)
	   DO ib=1,mesh(2)
		   DO  ia=1,mesh(1) 
                          Integrated_density = Integrated_density + (dV * Wavefunction(ia , ib, ic) * Wavefunction(ia , ib, ic))
                          Integrated_volume  = Integrated_volume + dV
                          ExpectationValue = ExpectationValue +  ( dV  *  Wavefunction(ia , ib, ic) * Wavefunction(ia , ib, ic) *field(ia , ib, ic) )
                   ENDDO     
             ENDDO
   ENDDO
   PRINT*, "    ......done "
END SUBROUTINE CalculateFirstOrderEnergy

!***********************************************
      SUBROUTINE CalculateDipoleOperatorWF(DipoleOperator, Integrated_volume, Integrated_density, dimensions,  X0Cell, dV, cell, mesh, Wavefunction)
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE

   REAL(dbl),   INTENT(out)        ::  Integrated_density
   REAL(dbl),   INTENT(out)        ::  Integrated_volume
   REAL(dbl),   INTENT(out)        ::  DipoleOperator(3)
   INTEGER,     INTENT(in)         ::  dimensions
   INTEGER,     INTENT(in)         ::  mesh(dimensions)
   REAL(dbl),   INTENT(in)         ::  X0Cell(dimensions)
   REAL(dbl),   INTENT(in)         ::  cell(dimensions,dimensions)
   REAL(dbl),   INTENT(in)         ::  dV
   REAL(dbl),   INTENT(in)         ::  Wavefunction( mesh(1), mesh(2), mesh(3) )
   REAL(dbl)   :: PositionXYZ(dimensions)
   INTEGER :: ia, ib, ic, idimension    

   PRINT*, "    ...Calculating Dipole operator associated with the pseudo-charge density"
   DipoleOperator(:) = ZERO
   Integrated_density= ZERO
   Integrated_volume   = ZERO
    PositionXYZ(:)  = ZERO

   DO ic=1,mesh(3)
	   DO ib=1,mesh(2)
		   DO  ia=1,mesh(1) 
                          Integrated_density = Integrated_density + (dV * Wavefunction(ia , ib, ic) * Wavefunction(ia , ib, ic))
                          Integrated_volume  = Integrated_volume + dV

                          DO idimension=1, dimensions
                                 PositionXYZ(idimension) = ((ia-1.5)*cell(1,idimension))  + ((ib-1.5)*cell(2,idimension))   + ((ic-1.5)*cell(3,idimension))  +X0Cell(idimension)   
                                 !PositionXYZ(idimension) = ((ia-1)*cell(1,idimension)/mesh(1))  + ((ib-1)*cell(2,idimension)/mesh(2))   + ((ic-1)*cell(3,idimension)/mesh(3))  +X0Cell(idimension)   
                                 DipoleOperator(idimension) = DipoleOperator(idimension) + ( PositionXYZ(idimension)* dV  * Wavefunction(ia , ib, ic) * Wavefunction(ia , ib, ic) )
                          ENDDO
                   ENDDO     
             ENDDO
   ENDDO
   PRINT*, "    ......done "
END SUBROUTINE CalculateDipoleOperatorWF


!***********************************************
      SUBROUTINE CalculateDipoleOperator(DipoleOperator, Integrated_volume, Integrated_density, dimensions,  X0Cell, dV, cell, mesh, density)
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE

   REAL(dbl),   INTENT(out)        ::  Integrated_density
   REAL(dbl),   INTENT(out)        ::  Integrated_volume
   REAL(dbl),   INTENT(out)        ::  DipoleOperator(3)
   INTEGER,     INTENT(in)         ::  dimensions
   INTEGER,     INTENT(in)         ::  mesh(dimensions)
   REAL(dbl),   INTENT(in)         ::  X0Cell(dimensions)
   REAL(dbl),   INTENT(in)         ::  cell(dimensions,dimensions)
   REAL(dbl),   INTENT(in)         ::  dV
   REAL(dbl),   INTENT(in)         ::  density( mesh(1), mesh(2), mesh(3) )
   REAL(dbl)   :: PositionXYZ(dimensions)
   INTEGER :: ia, ib, ic, idimension    

   PRINT*, "    ...Calculating Dipole operator "
   DipoleOperator(:) = ZERO
   Integrated_density= ZERO
   Integrated_volume   = ZERO


   DO ic=1,mesh(3)
	   DO ib=1,mesh(2)
		   DO  ia=1,mesh(1) 
                          Integrated_density = Integrated_density + dV * density(ia , ib, ic) 
                          Integrated_volume  = Integrated_volume + dV

                          DO idimension=1, dimensions
                                 PositionXYZ(idimension) = ((ia-1.5)*cell(1,idimension))  + ((ib-1.5)*cell(2,idimension))   + ((ic-1.5)*cell(3,idimension))  +X0Cell(idimension)   
                                 DipoleOperator(idimension) = DipoleOperator(idimension) + ( PositionXYZ(idimension)* dV  * density(ia , ib, ic) )
                          ENDDO
                   ENDDO     
             ENDDO
   ENDDO
   PRINT*, "    ......done "
END SUBROUTINE CalculateDipoleOperator



!***********************************************
      SUBROUTINE CalculateDipoleComponentDENSITYDIFF(DipoleComponent, Integrated_volume, Integrated_density, dimensions,  X0Cell, dV, cell, mesh, CenterOfCharges, density)
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE

   REAL(dbl),   INTENT(out)        ::  Integrated_density
   REAL(dbl),   INTENT(out)        ::  Integrated_volume
   REAL(dbl),   INTENT(out)        ::  DipoleComponent(mesh(3))
   REAL(dbl),   INTENT(out)        ::  RenormalizedPosition(mesh(3))
   INTEGER,     INTENT(in)         ::  dimensions
   INTEGER,     INTENT(in)         ::  mesh(dimensions)
   REAL(dbl),   INTENT(in)         ::  X0Cell(dimensions)
   REAL(dbl),   INTENT(in)         ::  cell(dimensions,dimensions)
   REAL(dbl),   INTENT(in)         ::  dV
   REAL(dbl),   INTENT(in)         ::  CenterOfCharges(dimensions)
   REAL(dbl),   INTENT(in)         ::  density( mesh(1), mesh(2), mesh(3) )
   REAL(dbl)   :: PositionXYZ(dimensions)
   INTEGER :: ia, ib, ic, idimension    

   PRINT*, "    ...Calculating Dipole operator "
   DipoleComponent(:) = ZERO
   Integrated_density= ZERO
   Integrated_volume   = ZERO
   idimension=3
   DO ic=1,mesh(3)
	   DipoleComponent(ic) = ZERO
	   DO ib=1,mesh(2)
		   DO  ia=1,mesh(1) 
                          Integrated_density = Integrated_density + dV * density(ia , ib, ic) 
                          Integrated_volume  = Integrated_volume + dV
	                  PositionXYZ(idimension) = ((ia-1.5)*cell(1,idimension))  + ((ib-1.5)*cell(2,idimension))   + ((ic-1.5)*cell(3,idimension)) + X0Cell(idimension)   - CenterOfCharges(idimension)
			  DipoleComponent(ic)=DipoleComponent(ic)+(dV*PositionXYZ(idimension)*density(ia,ib,ic))

                   ENDDO     
             ENDDO

             RenormalizedPosition(ic)= ((ic-1.5)*cell(3,idimension)) + X0Cell(idimension)   - CenterOfCharges(idimension)
   ENDDO
   PRINT*, "    ......done "
   ! 

END SUBROUTINE CalculateDipoleComponentDENSITYDIFF



!***********************************************
      SUBROUTINE  CalculateVolume(Volume, dV, dimensions, cell, grid)
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE

   INTEGER,     INTENT(in)         ::  dimensions
   INTEGER,     INTENT(in)         ::  grid(dimensions)
   REAL(dbl),   INTENT(in)         ::  cell(dimensions,dimensions)
   REAL(dbl),   INTENT(out)        ::  Volume
   REAL(dbl),   INTENT(out)        ::  dV
   REAL(dbl)         :: Crossaux(dimensions), aux(dimensions)

        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...A vector: ",  (grid(1))*cell(1,1), " | " ,  (grid(1))*cell(1,2), " | " ,  (grid(1))*cell(1,3) 
        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...B vector: ",  (grid(2))*cell(2,1), " | " ,  (grid(2))*cell(2,2), " | " ,  (grid(2))*cell(2,3) 
        PRINT'(A18,F12.6,A3,F12.6,A3,F12.6)', "    ...C vector: ",  (grid(3))*cell(3,1), " | " ,  (grid(3))*cell(3,2), " | " ,  (grid(3))*cell(3,3) 
        Crossaux(1) = ((grid(1))*cell(1,2)*(grid(2))*cell(2,3))  -((grid(2))*cell(2,2)*(grid(1))*cell(1,3) )
        Crossaux(2) = ((grid(1))*cell(1,3)*(grid(2))*cell(2,1))  -((grid(2))*cell(2,3)*(grid(1))*cell(1,1))
        Crossaux(3) = ((grid(1))*cell(1,1)*(grid(2))*cell(2,2))  -((grid(2))*cell(2,1)*(grid(1))*cell(1,2))
        aux(1) = (grid(3))*cell(3,1)
        aux(2) = (grid(3))*cell(3,2)
        aux(3) = (grid(3))*cell(3,3)
        Volume = dot_product(Crossaux(:), aux(:))
        dV = Volume/  (grid(1)*grid(2)*grid(3))
        PRINT'(A28,F12.4)', "    ...Volume of the cell [Ang^3]: ", Volume
        PRINT'(A28,F12.6)', "    ...Voxel [Ang^3]:              ", dV

END  SUBROUTINE CalculateVolume



!***********************************************
      SUBROUTINE  CalculateCenterOfMass(CenterOfMass, dimensions, Nb_atoms, Znumber, atoms_position)
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE

   !
   INTEGER,     INTENT(in)         ::  dimensions
   REAL(dbl),   INTENT(out)         ::  CenterOfMass(3)

   INTEGER,     INTENT(in)         ::  Nb_atoms
   INTEGER,   INTENT(in)           ::  Znumber(Nb_atoms)     !Znumber( Nb_atoms )
   REAL(dbl),   INTENT(in)         ::  atoms_position( Nb_atoms, dimensions)

   REAL(dbl) ::  Mass_atoms( Nb_atoms )
   REAL(dbl) ::  Total_Mass
   INTEGER :: iatom, idimension

   DO iatom=1, Nb_atoms 
           IF (  Znumber(iatom) == 1 ) THEN 
                Mass_atoms(iatom) =  1.00794
           ELSE IF (  Znumber(iatom) == 2 ) THEN 
                Mass_atoms(iatom) =  4.002602
           ELSE IF ( Znumber(iatom) == 3 ) THEN 
                Mass_atoms(iatom) =   6.941
           ELSE IF (  Znumber(iatom) == 4  ) THEN 
                Mass_atoms(iatom) =  9.012182
           ELSE IF (    Znumber(iatom) == 5) THEN 
                Mass_atoms(iatom) =10.811
           ELSE IF (   Znumber(iatom) == 6  ) THEN 
                Mass_atoms(iatom) = 12.0107
           ELSE IF ( Znumber(iatom) == 7 ) THEN 
                Mass_atoms(iatom) =   14.0067
           ELSE IF (Znumber(iatom) == 8  ) THEN 
                Mass_atoms(iatom) =    15.9994
           ELSE IF (   Znumber(iatom) == 9  ) THEN 
                Mass_atoms(iatom) = 18.9984032
           ELSE IF ( Znumber(iatom) ==10  ) THEN 
                Mass_atoms(iatom) =   20.1797
           ELSE IF ( Znumber(iatom) == 11  ) THEN 
                Mass_atoms(iatom) =   22.98976928
           ELSE IF ( Znumber(iatom) == 12  ) THEN 
                Mass_atoms(iatom) =   24.3050
           ELSE IF (Znumber(iatom) == 13 ) THEN 
                Mass_atoms(iatom) =    26.9815386
           ELSE IF (Znumber(iatom) == 14  ) THEN 
                Mass_atoms(iatom) =    28.0855
           ELSE IF (Znumber(iatom) == 15  ) THEN 
                Mass_atoms(iatom) =    30.973762
           ELSE IF ( Znumber(iatom) == 16  ) THEN 
                Mass_atoms(iatom) =    32.065
           ELSE IF (  Znumber(iatom) == 17 ) THEN 
                Mass_atoms(iatom) =   35.453
           ELSE IF (  Znumber(iatom) == 18 ) THEN 
                Mass_atoms(iatom) = 39.948    
           ELSE IF ( Znumber(iatom) == 79 ) THEN 
                Mass_atoms(iatom) =196.966569    
           ELSE 
               PRINT*, "Warning, atoms not identified",iatom, Znumber(iatom)
               STOP
           ENDIF
   ENDDO
   CenterOfMass(:) = ZERO
   Total_Mass = ZERO

   DO iatom=1, Nb_atoms 
       DO idimension=1, dimensions
             CenterOfMass(idimension) = CenterOfMass(idimension)+(atoms_position(iatom,idimension)  * Mass_atoms(iatom) )
       ENDDO 
        Total_Mass =  Total_Mass + Mass_atoms(iatom) 
   ENDDO
   CenterOfMass(:) =   CenterOfMass(:) / Total_Mass

        PRINT'(A28,3F12.4)', "    ...Center of Mass [Ang]:           ",CenterOfMass(1) ,CenterOfMass(2) , CenterOfMass(3)
        PRINT'(A28,F12.6)', "    ...Total Mass of the system [AMU]: ", Total_Mass


END  SUBROUTINE CalculateCenterOfMass

!***********************************************
      SUBROUTINE  CalculateCenterOfCharge(CenterOfCharge, Total_Charge,  dimensions, Nb_atoms, Znumber, atoms_position)
   !***********************************************
   USE constants,                          ONLY : PI, ZERO, ONE, CZERO, CONE, CI, EPS_m5, EPS_m4, EPS_m2, EPS_m1 !ok
   USE kinds,                              ONLY : dbl !ok
   IMPLICIT NONE

   !
   REAL(dbl),   INTENT(out)         ::  CenterOfCharge(3)
   REAL(dbl),   INTENT(out)         ::  Total_Charge
   INTEGER,     INTENT(in)         ::  dimensions
   INTEGER,     INTENT(in)         ::  Nb_atoms
   INTEGER,   INTENT(in)           ::  Znumber(Nb_atoms)     !Znumber( Nb_atoms )
   REAL(dbl),   INTENT(in)         ::  atoms_position( Nb_atoms, dimensions)

   REAL(dbl) ::  Charge_atoms( Nb_atoms )

   INTEGER :: iatom, idimension

   DO iatom=1, Nb_atoms 
           IF (  Znumber(iatom) == 1 ) THEN 
                Charge_atoms(iatom) =  1*ONE
           ELSE IF (  Znumber(iatom) == 2 ) THEN 
                Charge_atoms(iatom) =  2*ONE
           ELSE IF ( Znumber(iatom) == 3 ) THEN 
                Charge_atoms(iatom) =  1*ONE
           ELSE IF (  Znumber(iatom) == 4  ) THEN 
                Charge_atoms(iatom) = 2*ONE
           ELSE IF (    Znumber(iatom) == 5) THEN 
                Charge_atoms(iatom) =3*ONE
           ELSE IF (   Znumber(iatom) == 6  ) THEN 
                Charge_atoms(iatom) =4*ONE
           ELSE IF ( Znumber(iatom) == 7 ) THEN 
                Charge_atoms(iatom) =   5*ONE
           ELSE IF (Znumber(iatom) == 8  ) THEN 
                Charge_atoms(iatom) =   6*ONE
           ELSE IF (   Znumber(iatom) == 9  ) THEN 
                Charge_atoms(iatom) = 7*ONE
           ELSE IF ( Znumber(iatom) ==10  ) THEN 
                Charge_atoms(iatom) =   8*ONE
           ELSE IF ( Znumber(iatom) == 11  ) THEN 
                Charge_atoms(iatom) =   1*ONE
           ELSE IF ( Znumber(iatom) == 12  ) THEN 
                Charge_atoms(iatom) =   2*ONE
           ELSE IF (Znumber(iatom) == 13 ) THEN 
                Charge_atoms(iatom) =    3*ONE
           ELSE IF (Znumber(iatom) == 14  ) THEN 
                Charge_atoms(iatom) =    4*ONE
           ELSE IF (Znumber(iatom) == 15  ) THEN 
                Charge_atoms(iatom) =    5*ONE
           ELSE IF ( Znumber(iatom) == 16  ) THEN 
                Charge_atoms(iatom) =    6*ONE
           ELSE IF (  Znumber(iatom) == 17 ) THEN 
                Charge_atoms(iatom) =   7*ONE
           ELSE IF (  Znumber(iatom) == 18 ) THEN 
                Charge_atoms(iatom) = 8 *ONE  
           ELSE IF ( Znumber(iatom) == 79 ) THEN 
                Charge_atoms(iatom) =11*ONE   
           ELSE 
               PRINT*, "Warning, atoms not identified", iatom, Znumber(iatom)
               STOP
           ENDIF
   ENDDO
   CenterOfCharge(:) = ZERO
   Total_Charge = ZERO

   DO iatom=1, Nb_atoms 
       DO idimension=1, dimensions
             CenterOfCharge(idimension) = CenterOfCharge(idimension) +(atoms_position(iatom,idimension)  * Charge_atoms(iatom) )
       ENDDO 
        Total_Charge =  Total_Charge + Charge_atoms(iatom) 
   ENDDO
   CenterOfCharge(:) =   CenterOfCharge(:) / Total_Charge

        PRINT'(A28,3F12.4)', "    ...Center of Charge [Ang]:           ",CenterOfCharge(1) ,CenterOfCharge(2) , CenterOfCharge(3)
        PRINT'(A28,F12.6)', "    ...Total Charge of the system [AMU]: ", Total_Charge


END  SUBROUTINE CalculateCenterOfCharge

END MODULE calculatequantities_Module
