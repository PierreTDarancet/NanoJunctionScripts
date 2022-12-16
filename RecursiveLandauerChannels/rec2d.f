       program rec2d
c       recursion sur un réseau 2D
!nb iter
        parameter(NIter=40)
!Energie on site
        parameter(E0=3.1)
!hopping
        parameter(beta=0.9)

!gc
        implicit complex (R)


!declaration
       REAL     ::  WF_n(-NIter-2:NIter+2,-NIter-2:NIter+2)
       REAL     ::  WF_n_m(-NIter-2:NIter+2,-NIter-2:NIter+2)
       REAL     ::  WF_n_p(-NIter-2:NIter+2,-NIter-2:NIter+2)
       REAL     ::  A(0:NIter)
       REAL     ::  B(0:NIter+1)


       !Initialisation des variables
       WF_n(:,:)=0
       WF_n_m(:,:)=0
       WF_n_p(:,:)=0
       WF_n(0,0)=1.0
       B(0)=0.
       do n=0,NIter
          do i=-n-1,n+1,1
             do j=-n-1,n+1,1
                WF_n_m(i,j)=E0*WF_n(i,j)+beta*(WF_n(i-1,j))
                WF_n_m(i,j)=WF_n_m(i,j)+beta*(WF_n(i+1,j))
                WF_n_m(i,j)=WF_n_m(i,j)+beta*(WF_n(i,j+1))
                WF_n_m(i,j)=WF_n_m(i,j)+beta*(WF_n(i,j-1))
                WF_n_m(i,j)=WF_n_m(i,j)-b(n)*WF_n_p(i,j)
                A(n)=A(n)+WF_n_m(i,j)*WF_n(i,j)
             enddo
          enddo
          do i=-n-1,n+1
             do j=-n-1,n+1
c      WF_n_m(i,j)=WF_n_m(i,j)-a(n)*WF_n(i,j)-b(n)*WF_n_p(i,j)
                WF_n_m(i,j)=WF_n_m(i,j)-a(n)*WF_n(i,j)
                b(n+1)=b(n+1)+WF_n_m(i,j)*WF_n_m(i,j)
             enddo
          enddo
       b(n+1)=sqrt(b(n+1))
          do i=-n-1,n+1
             do j=-n-1,n+1
                WF_n_m(i,j)=WF_n_m(i,j)/b(n+1)
                WF_n_p(i,j)=WF_n(i,j)
                WF_n(i,j)=WF_n_m(i,j)
             enddo
         enddo
c      print*,n,a(n),b(n+1)
       enddo
       do m=0,200
           E=-0.5+m/25.
           RR=E-a(NIter)
           RR = REAL(1.0/(2*b(NIter+1)*b(NIter+1)))*(RR-sqrt(RR**2-4*b(NIter+1)**2))
           do n=NIter,0,-1
              RR=1/(E-a(n)-b(n+1)**2*RR)
           enddo
       DOS=-AIMAG(RR)/3.141592
       print*,E,DOS 
       enddo
       end
