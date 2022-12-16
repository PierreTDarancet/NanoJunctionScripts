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


