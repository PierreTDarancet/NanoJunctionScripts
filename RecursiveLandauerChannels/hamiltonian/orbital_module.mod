  �  <   k820309    v          9.0         yF                                                                                                       
       orbitale_m.f90 ORBITAL_MODULE          ORB ORBITAL_ALLOC ORBITAL_ALLOCATE ORBITAL_DEALLOCATE ORBITAL_INIT PRINT_ORBITAL                                          
       NSTRX              � @ @                     
                                                
       CZERO CONE ZERO                  @                       
       ORBITALE                  @                       
       N_ORB              � @ @                      
                        @                       
       FILE_OPEN FILE_CLOSE                  @                       
      ORB_UNIT AUX2_UNIT IONAME                  @                  	     
       STDIN STDOUT               !@                   
     '(            #NATURE    #COORD    #INIT             � $                                             � $                                      
  p      p        p                       � $                                                  @                        '            #DUMMY             � $                                            @@   !                            #     @                                        
   #FILE_OPEN%TRIM    #FILE_OPEN%PRESENT    #FILE_OPEN%ASSOCIATED    #UNIT    #FILENAME    #PATH    #ROOT    #STATUS    #ACCESS    #RECL    #FORM    #POSITION    #ACTION              @                         TRIM           @                         PRESENT           @                         ASSOCIATED       
   @                                   
   @                                 1       
  @                                 1       
  @                                 1       
  @                                 1       
  @                                 1       
  @                                   
  @                                 1       
  @                                 1       
  @                                 1 #     @                                           #FILE_CLOSE%TRIM     #FILE_CLOSE%PRESENT !   #FILE_CLOSE%ASSOCIATED "   #UNIT #   #PATH $   #ACTION %             @                          TRIM           @                    !     PRESENT           @                    "     ASSOCIATED       
   @                      #             
  @                     $            1       
  @                     %            1 #     @                         &                  #IONAME%PRESENT '   #IONAME%TRIM (   #IONAME%LEN_TRIM )   #DATA *   #FILENAME +   #LPOSTFIX ,   #LPATH -             @                    '     PRESENT           @                    (     TRIM           @                    )     LEN_TRIM       
   @                     *            1         @                     +             1       
  @                      ,             
  @                      -                                        .                                      5                                 /                                      6     @@                       0        (            &                       #ORBITALE 
          @                        1        #     @                          2                   #ORBITAL_ALLOCATE%ABS 3             @                    3     ABS #     @                          4                   #ORBITAL_DEALLOCATE%ALLOCATED 5   #ORBITAL_DEALLOCATE%ABS 6             @                    5     ALLOCATED           @                    6     ABS #     @                          7                    #     @                          8                   #PRINT_ORBITAL%TRIM 9   #PRINT_ORBITAL%ABS :             @                    9     TRIM           @                    :     ABS    �   &      fn#fn $   �   ]   b   uapp(ORBITAL_MODULE      :   J  PARAMETERS    Y  4   J  KINDS    �  D   J  CONSTANTS     �  =   J  IDENTITY_MODULE $     :   J  DIM_VARIABLE_MODULE    H  4   J  IOTK_MODULE    |  I   J  FILES_MODULE    �  N   J  IO_MODULE !     A   J  IO_GLOBAL_MODULE )   T  ]       ORBITALE+IDENTITY_MODULE 0   �  8   a   ORBITALE%NATURE+IDENTITY_MODULE /   �  p   a   ORBITALE%COORD+IDENTITY_MODULE .   Y  8   a   ORBITALE%INIT+IDENTITY_MODULE )   �  G       IOTK_DUMMYTYPE+IOTK_BASE /   �  8   a   IOTK_DUMMYTYPE%DUMMY+IOTK_BASE *     0       N_ORB+DIM_VARIABLE_MODULE '   @  �       FILE_OPEN+FILES_MODULE 1   3  1      FILE_OPEN%TRIM+FILES_MODULE=TRIM 7   d  4      FILE_OPEN%PRESENT+FILES_MODULE=PRESENT =   �  7      FILE_OPEN%ASSOCIATED+FILES_MODULE=ASSOCIATED ,   �  0   e   FILE_OPEN%UNIT+FILES_MODULE 0   �  8   e   FILE_OPEN%FILENAME+FILES_MODULE ,   7  8   e   FILE_OPEN%PATH+FILES_MODULE ,   o  8   e   FILE_OPEN%ROOT+FILES_MODULE .   �  8   e   FILE_OPEN%STATUS+FILES_MODULE .   �  8   e   FILE_OPEN%ACCESS+FILES_MODULE ,     0   e   FILE_OPEN%RECL+FILES_MODULE ,   G  8   e   FILE_OPEN%FORM+FILES_MODULE 0     8   e   FILE_OPEN%POSITION+FILES_MODULE .   �  8   e   FILE_OPEN%ACTION+FILES_MODULE (   �  �       FILE_CLOSE+FILES_MODULE 2   �	  1      FILE_CLOSE%TRIM+FILES_MODULE=TRIM 8   �	  4      FILE_CLOSE%PRESENT+FILES_MODULE=PRESENT >   �	  7      FILE_CLOSE%ASSOCIATED+FILES_MODULE=ASSOCIATED -   /
  0   e   FILE_CLOSE%UNIT+FILES_MODULE -   _
  8   e   FILE_CLOSE%PATH+FILES_MODULE /   �
  8   e   FILE_CLOSE%ACTION+FILES_MODULE !   �
  �       IONAME+IO_MODULE 1   v  4      IONAME%PRESENT+IO_MODULE=PRESENT +   �  1      IONAME%TRIM+IO_MODULE=TRIM 3   �  5      IONAME%LEN_TRIM+IO_MODULE=LEN_TRIM &     8   e   IONAME%DATA+IO_MODULE *   H  8   e   IONAME%FILENAME+IO_MODULE *   �  0   e   IONAME%LPOSTFIX+IO_MODULE '   �  0   e   IONAME%LPATH+IO_MODULE '   �  U       STDIN+IO_GLOBAL_MODULE (   5  U       STDOUT+IO_GLOBAL_MODULE    �  j       ORB    �  0       ORBITAL_ALLOC !   $  V       ORBITAL_ALLOCATE %   z  0      ORBITAL_ALLOCATE%ABS #   �  z       ORBITAL_DEALLOCATE -   $  6      ORBITAL_DEALLOCATE%ALLOCATED '   Z  0      ORBITAL_DEALLOCATE%ABS    �  <       ORBITAL_INIT    �  k       PRINT_ORBITAL #   1  1      PRINT_ORBITAL%TRIM "   b  0      PRINT_ORBITAL%ABS 