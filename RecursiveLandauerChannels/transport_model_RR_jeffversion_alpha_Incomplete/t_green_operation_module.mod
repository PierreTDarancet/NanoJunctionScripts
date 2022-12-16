  *  ^   k820309              10.1        ÐĄĖK                                                                                                           
       green_functions.f90 T_GREEN_OPERATION_MODULE              LEAD_GREEN_FUNCTION CALCUL_GAMMA CALCUL_G_C                   @ @                    
       DBL                                             
       ZERO ONE CZERO CONE CI PI EPS_M3 EPS_M5                      @                      
       MAT_MUL MAT_SV                                              u #ZMAT_MUL    #DMAT_MUL    #         @     @                                    #ZMAT_MUL%SIZE    #ZMAT_MUL%CONJG    #C    #A    #OPA 	   #B 
   #OPB    #M    #N    #K                  @                         SIZE               @                         CONJG            @                                         "              &                   &                                                     
  @                                                       &                   &                                                     
   @                      	                                     
  @                     
                    !             &                   &                                                     
   @                                                           
   @                                           
   @                                           
   @                                 #         @     @                                     #DMAT_MUL%SIZE    #C    #A    #OPA    #B    #OPB    #M    #N    #K                  @                         SIZE            @                                        
 %              &                   &                                                     
  @                                        
 #             &                   &                                                     
   @                                                           
  @                                        
 $             &                   &                                                     
   @                                                           
   @                                           
   @                                           
   @                                                                              u #ZMAT_SV    #ZMAT_SV_1 "   #DMAT_SV +   #DMAT_SV_1 4   #         @     @                                    #ZMAT_SV%SIZE    #ZMAT_SV%ABS    #ZMAT_SV%PRESENT    #N    #NRHS    #A    #B     #IERR !                 @                         SIZE               @                         ABS               @                         PRESENT           
   @                                           
   @                                           
  @                                                      &                   &                                                     
 @                                                        &                   &                                                       @                      !            #         @     @                   "                  #ZMAT_SV_1%SIZE #   #ZMAT_SV_1%ABS $   #ZMAT_SV_1%PRESENT %   #N &   #NRHS '   #A (   #B )   #IERR *                 @                    #     SIZE               @                    $     ABS               @                    %     PRESENT           
   @                      &                     
   @                      '                     
   @                     (                                 &                   &                                                     
 @                     )                                  &                                                       @                      *            #         @     @                   +                  #DMAT_SV%SIZE ,   #DMAT_SV%ABS -   #DMAT_SV%PRESENT .   #N /   #NRHS 0   #A 1   #B 2   #IERR 3                 @                    ,     SIZE               @                    -     ABS               @                    .     PRESENT           
   @                      /                     
   @                      0                     
  @                     1                   
              &                   &                                                     
 @                     2                   
               &                   &                                                       @                      3            #         @     @                   4                  #DMAT_SV_1%SIZE 5   #DMAT_SV_1%ABS 6   #DMAT_SV_1%PRESENT 7   #N 8   #NRHS 9   #A :   #B ;   #IERR <                 @                    5     SIZE               @                    6     ABS               @                    7     PRESENT           
   @                      8                     
   @                      9                     
   @                     :                   
              &                   &                                                     
 @                     ;                   
               &                                                       @                      <                                                 =                                                         #         @                          >                  #LEAD_GREEN_FUNCTION%SQRT ?   #LEAD_GREEN_FUNCTION%ABS @   #LEAD_GREEN_FUNCTION%TRIM A   #LEAD_GREEN_FUNCTION%REAL B   #G_LEAD C   #H00 E   #H01 F   #ENE G   #DIMX D                 @                    ?     SQRT               @                    @     ABS               @                    A     TRIM               @                    B     REAL          D  @                     C                           p        5  p        r D   p          5  p        r D     5  p        r D       5  p        r D     5  p        r D                              
   @                     E                          p        5  p        r D   p          5  p        r D     5  p        r D       5  p        r D     5  p        r D                              
   @                     F                          p        5  p        r D   p          5  p        r D     5  p        r D       5  p        r D     5  p        r D                               
   @                     G     
                
   @                      D           #         @                          H                  #CALCUL_GAMMA%TRANSPOSE I   #CALCUL_GAMMA%CONJG J   #CALCUL_GAMMA%TRIM K   #GAMMA_LEAD L   #SIGMA_LEAD N   #G_LEAD O   #H_LEAD_C Q   #H_C_LEAD R   #DIMX P   #DIMC M                 @                    I     TRANSPOSE               @                    J     CONJG               @                    K     TRIM          D  @                     L                           p        5  p        r M   p          5  p        r M     5  p        r M       5  p        r M     5  p        r M                              D @@                     N                           p        5  p        r M   p          5  p        r M     5  p        r M       5  p        r M     5  p        r M                              
@ @@                     O                          p        5  p        r P   p          5  p        r P     5  p        r P       5  p        r P     5  p        r P                              
@ @@                     Q                          p        5  p        r P   p          5  p        r P     5  p        r M       5  p        r P     5  p        r M                              
@ @@                     R                          p        5  p        r M   p          5  p        r M     5  p        r P       5  p        r M     5  p        r P                               
@ @@                      P                     
@ @@                      M           #         @                          S                  #CALCUL_G_C%TRIM T   #GC U   #H00 W   #SIGMA_L X   #SIGMA_R Y   #ENE Z   #DIMC V                 @                    T     TRIM          D @@                     U                     
      p        5  p        r V   p          5  p        r V     5  p        r V       5  p        r V     5  p        r V                              
   @                     W                          p        5  p        r V   p          5  p        r V     5  p        r V       5  p        r V     5  p        r V                              
   @                     X                          p        5  p        r V   p          5  p        r V     5  p        r V       5  p        r V     5  p        r V                              
   @                     Y                          p        5  p        r V   p          5  p        r V     5  p        r V       5  p        r V     5  p        r V                               
   @                     Z                     
@ @@                      V                  5      fn#fn .   Õ   <   b   uapp(T_GREEN_OPERATION_MODULE      <   J  KINDS    M  `   J  CONSTANTS    ­  G   J  UTIL_MODULE (   ô  T       gen@MAT_MUL+UTIL_MODULE %   H  Ģ      ZMAT_MUL+UTIL_MODULE /   ë  5      ZMAT_MUL%SIZE+UTIL_MODULE=SIZE 1      6      ZMAT_MUL%CONJG+UTIL_MODULE=CONJG '   V     e   ZMAT_MUL%C+UTIL_MODULE '   ō     e   ZMAT_MUL%A+UTIL_MODULE )     H   e   ZMAT_MUL%OPA+UTIL_MODULE '   Ö     e   ZMAT_MUL%B+UTIL_MODULE )   r  H   e   ZMAT_MUL%OPB+UTIL_MODULE '   š  8   e   ZMAT_MUL%M+UTIL_MODULE '   ō  8   e   ZMAT_MUL%N+UTIL_MODULE '   *  8   e   ZMAT_MUL%K+UTIL_MODULE %   b        DMAT_MUL+UTIL_MODULE /   ņ  5      DMAT_MUL%SIZE+UTIL_MODULE=SIZE '   &     e   DMAT_MUL%C+UTIL_MODULE '   Â     e   DMAT_MUL%A+UTIL_MODULE )   ^  H   e   DMAT_MUL%OPA+UTIL_MODULE '   Ķ     e   DMAT_MUL%B+UTIL_MODULE )   B	  H   e   DMAT_MUL%OPB+UTIL_MODULE '   	  8   e   DMAT_MUL%M+UTIL_MODULE '   Â	  8   e   DMAT_MUL%N+UTIL_MODULE '   ú	  8   e   DMAT_MUL%K+UTIL_MODULE '   2
  p       gen@MAT_SV+UTIL_MODULE $   Ē
  Ą      ZMAT_SV+UTIL_MODULE .   C  5      ZMAT_SV%SIZE+UTIL_MODULE=SIZE ,   x  4      ZMAT_SV%ABS+UTIL_MODULE=ABS 4   Ž  8      ZMAT_SV%PRESENT+UTIL_MODULE=PRESENT &   ä  8   e   ZMAT_SV%N+UTIL_MODULE )     8   e   ZMAT_SV%NRHS+UTIL_MODULE &   T     e   ZMAT_SV%A+UTIL_MODULE &   ð     e   ZMAT_SV%B+UTIL_MODULE )     8   e   ZMAT_SV%IERR+UTIL_MODULE &   Ä  §      ZMAT_SV_1+UTIL_MODULE 0   k  5      ZMAT_SV_1%SIZE+UTIL_MODULE=SIZE .      4      ZMAT_SV_1%ABS+UTIL_MODULE=ABS 6   Ô  8      ZMAT_SV_1%PRESENT+UTIL_MODULE=PRESENT (     8   e   ZMAT_SV_1%N+UTIL_MODULE +   D  8   e   ZMAT_SV_1%NRHS+UTIL_MODULE (   |     e   ZMAT_SV_1%A+UTIL_MODULE (        e   ZMAT_SV_1%B+UTIL_MODULE +     8   e   ZMAT_SV_1%IERR+UTIL_MODULE $   Ô  Ą      DMAT_SV+UTIL_MODULE .   u  5      DMAT_SV%SIZE+UTIL_MODULE=SIZE ,   Š  4      DMAT_SV%ABS+UTIL_MODULE=ABS 4   Þ  8      DMAT_SV%PRESENT+UTIL_MODULE=PRESENT &     8   e   DMAT_SV%N+UTIL_MODULE )   N  8   e   DMAT_SV%NRHS+UTIL_MODULE &        e   DMAT_SV%A+UTIL_MODULE &   "     e   DMAT_SV%B+UTIL_MODULE )   ū  8   e   DMAT_SV%IERR+UTIL_MODULE &   ö  §      DMAT_SV_1+UTIL_MODULE 0     5      DMAT_SV_1%SIZE+UTIL_MODULE=SIZE .   Ō  4      DMAT_SV_1%ABS+UTIL_MODULE=ABS 6     8      DMAT_SV_1%PRESENT+UTIL_MODULE=PRESENT (   >  8   e   DMAT_SV_1%N+UTIL_MODULE +   v  8   e   DMAT_SV_1%NRHS+UTIL_MODULE (   Ū     e   DMAT_SV_1%A+UTIL_MODULE (   J     e   DMAT_SV_1%B+UTIL_MODULE +   Î  8   e   DMAT_SV_1%IERR+UTIL_MODULE      h       DBL+KINDS $   n  č       LEAD_GREEN_FUNCTION )   V  5      LEAD_GREEN_FUNCTION%SQRT (     4      LEAD_GREEN_FUNCTION%ABS )   ŋ  5      LEAD_GREEN_FUNCTION%TRIM )   ô  5      LEAD_GREEN_FUNCTION%REAL +   )    a   LEAD_GREEN_FUNCTION%G_LEAD (   E    a   LEAD_GREEN_FUNCTION%H00 (   a    a   LEAD_GREEN_FUNCTION%H01 (   }  8   a   LEAD_GREEN_FUNCTION%ENE )   ĩ  8   a   LEAD_GREEN_FUNCTION%DIMX    í  į       CALCUL_GAMMA '   Ô  :      CALCUL_GAMMA%TRANSPOSE #     6      CALCUL_GAMMA%CONJG "   D  5      CALCUL_GAMMA%TRIM (   y    a   CALCUL_GAMMA%GAMMA_LEAD (       a   CALCUL_GAMMA%SIGMA_LEAD $   ą     a   CALCUL_GAMMA%G_LEAD &   Í!    a   CALCUL_GAMMA%H_LEAD_C &   é"    a   CALCUL_GAMMA%H_C_LEAD "   $  8   a   CALCUL_GAMMA%DIMX "   =$  8   a   CALCUL_GAMMA%DIMC    u$         CALCUL_G_C     %  5      CALCUL_G_C%TRIM    =%    a   CALCUL_G_C%GC    Y&    a   CALCUL_G_C%H00 #   u'    a   CALCUL_G_C%SIGMA_L #   (    a   CALCUL_G_C%SIGMA_R    ­)  8   a   CALCUL_G_C%ENE     å)  8   a   CALCUL_G_C%DIMC 