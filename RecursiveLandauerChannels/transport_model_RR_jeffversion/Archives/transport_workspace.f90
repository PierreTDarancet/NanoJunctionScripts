   COMPLEX(dbl), ALLOCATABLE :: work(:,:), work2(:,:)


!
! Contains Green function 
 COMPLEX(dbl), ALLOCATABLE :: gC(:,:)
 COMPLEX(dbl), ALLOCATABLE :: gR(:,:)
 COMPLEX(dbl), ALLOCATABLE :: gL(:,:)

!
! Contains local variables used in the calculation

   COMPLEX(dbl), ALLOCATABLE :: gamma_L(:,:)
   COMPLEX(dbl), ALLOCATABLE :: gamma_R(:,:)

   COMPLEX(dbl), ALLOCATABLE :: sigma_L(:,:)
   COMPLEX(dbl), ALLOCATABLE :: sigma_R(:,:)

   COMPLEX(dbl), ALLOCATABLE :: aux00_C(:,:)
   COMPLEX(dbl), ALLOCATABLE :: aux00_R(:,:)
   COMPLEX(dbl), ALLOCATABLE :: aux00_L(:,:)
   COMPLEX(dbl), ALLOCATABLE :: s00_C(:,:) 
   COMPLEX(dbl), ALLOCATABLE :: aux_LC(:,:)
   COMPLEX(dbl), ALLOCATABLE :: aux_CL(:,:)

   COMPLEX(dbl), ALLOCATABLE :: aux_RC(:,:)
   COMPLEX(dbl), ALLOCATABLE :: aux_CR(:,:)

   COMPLEX(dbl), ALLOCATABLE :: aux01_R(:,:)
   COMPLEX(dbl), ALLOCATABLE :: aux01_L(:,:)
   COMPLEX(dbl) :: work_scal
   REAL(dbl)    :: work_scal2


   ALLOCATE( gamma_R(Nb_states,Nb_states), STAT=ierr )
   ALLOCATE( sigma_R(Nb_states,Nb_states), STAT=ierr )
   ALLOCATE( gamma_L(Nb_states,Nb_states), STAT=ierr )
   ALLOCATE( sigma_L(Nb_states,Nb_states), STAT=ierr )
   ALLOCATE ( s00_C(Nb_states,Nb_states), STAT=ierr )

   ALLOCATE ( cond_aux(Nb_states), STAT=ierr )
   ALLOCATE ( work(Nb_states,Nb_states), STAT=ierr )
   ALLOCATE ( work2(Nb_states,Nb_states), STAT=ierr )
   ALLOCATE (  gC(Nb_states,Nb_states), STAT=ierr )
   ALLOCATE (  aux00_C(Nb_states,Nb_states), STAT=ierr )
   ALLOCATE ( aux00_R( Nb_states_lead, Nb_states_lead), STAT=ierr )
   ALLOCATE (  aux00_L( Nb_states_lead, Nb_states_lead), STAT=ierr )
   ALLOCATE (aux_LC( Nb_states_lead,Nb_states), STAT=ierr )
   ALLOCATE ( aux_CL(Nb_states, Nb_states_lead), STAT=ierr )
   ALLOCATE ( aux_RC( Nb_states_lead,Nb_states), STAT=ierr )
   ALLOCATE ( aux_CR(Nb_states, Nb_states_lead), STAT=ierr )
   ALLOCATE (aux01_R( Nb_states_lead, Nb_states_lead), STAT=ierr )
   ALLOCATE (aux01_L( Nb_states_lead, Nb_states_lead), STAT=ierr )
