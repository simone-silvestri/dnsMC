
      subroutine inter_c_s(n,f,df,dx)
      real f(0:N),df(N),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      call interpol_c_s6('n',n,f,df,dx,t1,a1,r)
      end 

      subroutine inter_s_c(n,f,df,dx)
      real F(N),dF(0:N),dx,t1(n+1,n+1),a1(n+1,n),r(n+1,n),td(n+1,n+1)
      call interpol_s_c6('n',n,f,df,dx,t1,a1,r)
      end

      subroutine deriv_c_s(n,f,df,dx)
      real f(0:N),df(N),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      call der1_c_s6('n',n,f,df,dx,t1,a1,r)
      end

      subroutine deriv_s_c(n,f,df,dx)
      real F(N),dF(0:N),dx ,t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      call der1_s_c6('n',n,f,df,dx,t1,a1,r)
      end

      subroutine inter_s_c_m(n,nrhs,f,df,dx)
      real F(N,nrhs),dF(0:N,nrhs),dx,t1(n+1,n+1),a1(n+1,n),r(n+1,n),td(n+1,n+1)
      call interpol_s_c6_m('n',n,nrhs,f,df,dx,t1,a1,r)
      end

      subroutine inter_c_s_m(n,nrhs,f,df,dx)
      real f(0:N,nrhs),df(N,nrhs),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      call interpol_c_s6_m('n',n,nrhs,f,df,dx,t1,a1,r)
      end 

      subroutine deriv_c_s_m(n,nrhs,f,df,dx)
      real f(0:N,nrhs),df(N,nrhs),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      call der1_c_s6_m('n',n,nrhs,f,df,dx,t1,a1,r)
      end

      subroutine deriv_s_c_m(n,nrhs,f,df,dx)
      real F(N),dF(0:N),dx ,t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      call der1_s_c6_m('n',n,nrhs,f,df,dx,t1,a1,r)
      end

      subroutine der2_s_v(n,f,df,dx)
      integer n
      real f(n),df(n),a(n),b(n),c(n),rhs(n),d(n)
      real ft(n)
      val = 0
      do i=2,n-1
	    b(i)=1.
        a(i)=1./10
        c(i)=1./10
        rhs(i) =6*(f(i+1)-2*f(i)+f(i-1))/(5*dx*dx)
       enddo
      do i=3,n-2
	    b(i)=1.
        a(i)=2./11
        c(i)=2./11
        rhs(i) =12*(f(i+1)-2*f(i)+f(i-1))/(11*dx*dx)+3*(f(i+2)-2*f(i)+f(i-2))/(44*dx*dx)
       enddo
       b(1)=1.
       a(1)=0.
       c(1)=0.5
       val = 0.0
       rhs(1)=(16.*val/5. -9./2.*f(1)+f(2)+(3./10.)*f(3))/(dx*dx)
       b(n)=1.
       c(n)=0.
       a(n)=0.5
       val = 1.0
       rhs(n)=(16.*val/5. -9./2.*f(n)+f(n-1)+(3./10.)*f(n-2))/(dx*dx)
       call trid(n,a,b,c,rhs,d)
       do i=1,n	
         df(i)=rhs(i)
       enddo
      end

      subroutine der2_s_n(n,f,df,dx)
      integer n
      real f(n),df(n),a(n),b(n),c(n),rhs(n),d(n)
      real ft(n)
      val = 0
      do i=2,n-1
	    b(i)=1.
        a(i)=1./10
        c(i)=1./10
        rhs(i) =6*(f(i+1)-2*f(i)+f(i-1))/(5*dx*dx)
       enddo
      do i=3,n-2
	    b(i)=1.
        a(i)=2./11
        c(i)=2./11
        rhs(i) =12*(f(i+1)-2*f(i)+f(i-1))/(11*dx*dx)+3*(f(i+2)-2*f(i)+f(i-2))/(44*dx*dx)
       enddo
       b(1)=1.
       a(1)=0.
       c(1)=-11./23
       val = 0.0
       rhs(1)=-24.*val/(23.*dx) +( -36./23.*f(1)+48*f(2)/23.-12.*f(3)/23.)/(dx*dx)

       b(n)=1.
       c(n)=0.
       a(n)=-11./23.
       val = 0.
       rhs(n)=24.*val/(23.*dx) +( -36./23.*f(n)+48*f(n-1)/23.-12.*f(n-2)/23.)/(dx*dx)
       call trid(n,a,b,c,rhs,d)
       do i=1,n	
         df(i)=rhs(i)
       enddo
      end

      subroutine der2_s(n,f,df,dx)
      integer n
      real f(n),df(n),a(n),b(n),c(n),rhs(n),d(n)
      real ft(n)
      val = 0
      do i=2,n-1
	    b(i)=1.
        a(i)=1./10
        c(i)=1./10
        rhs(i) =6*(f(i+1)-2*f(i)+f(i-1))/(5*dx*dx)
       enddo
      do i=3,n-2
	    b(i)=1.
        a(i)=2./11
        c(i)=2./11
        rhs(i) =12*(f(i+1)-2*f(i)+f(i-1))/(11*dx*dx)+3*(f(i+2)-2*f(i)+f(i-2))/(44*dx*dx)
       enddo
       b(1)=1.
       a(1)=0.
       c(1)=0.5
       rhs(1)=(16.*val/5. -9./2*f(1)+f(2)+3./10*f(3))/(dx*dx)
       b(n)=1.
       c(n)=0.
       a(n)=0.5
       rhs(n)=(16.*val/5. -9./2*f(n)+f(n-1)+3./10*f(n-2))/(dx*dx)
       call trid(n,a,b,c,rhs,d)
       do i=1,n	
         df(i)=rhs(i)
       enddo
c         do i=2,n-1
c         df(i)=(f(i+1)-2*f(i)+f(i-1))/(dx*dx)
c         enddo
c         df(1)=(-3*f(1)+f(2))/(dx*dx)
c         df(n)=(-3*f(n)+f(n-1))/(dx*dx)
      end

      subroutine der2_c(n,f,df,dx)
      integer n
      real f(n+1),df(n+1),a(n+1),b(n+1),c(n+1),d(n+1),rhs(n+1)

      do i=3,n-1
	    a(i)=2./11
	    b(i)=1.
         c(i)=2./11
         rhs(i) =12*(f(i+1)-2*f(i)+f(i-1))/(11*dx*dx)+3*(f(i+2)-2*f(i)+f(i-2))/(44*dx*dx)
      enddo
	    b(2)=1.
        a(2)=1./10
        c(2)=1./10
        rhs(2) =6*(f(3)-2*f(2)+f(1))/(5*dx*dx)
	    b(n)=1.
        a(n)=1./10
        c(n)=1./10
        rhs(n) =6*(f(n-1)-2*f(n)+f(n+1))/(5*dx*dx)
        a(1)=0
        b(1)=1
        c(1)=11. 
        rhs(1)=(13*f(1)-27*f(2)+15*f(3)-f(4))/(dx*dx)
        a(n+1)=11
        b(n+1)=1
        c(n+1)=0. 
        rhs(n+1)=(13*f(n+1)-27*f(n)+15*f(n-1)-f(n-2))/(dx*dx)
       call trid(n+1,a,b,c,rhs,d)
       do i=1,n+1	
         df(i)=rhs(i)
       enddo
      end



      subroutine interpol_c_s6(flag,N,F,dF,dx,t1,a1,r)
      character*1 flag
      real f(0:N),df(N),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      integer N
      real A(N),B(N),C(N),D(N),RHS(N)
      dxi=1./dx
      do i=2,n-1 
      A  (i  )=1./6
      B  (i  )=1.
      C  (i  )=1./6
      RHS(i  )=2./3*(F(i)+F(i-1  ))


      enddo
      alpha =.3
      fac1=(1./8)*(9+10*alpha)
      fac2=(1./8)*(6*alpha-1 )


      do i=3,n-2
        A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
        RHS (i) = fac1*(F(i)+F(i-1))*.5 + fac2*(F(i+1)+F(i-2))*.5
      enddo 
      if (flag.eq.'m') then
      do i=2,n-1
      a1(i,i )=2./3
      a1(i,i+1)=2./3
      enddo
      do i=3,n-2
	    a1(i,i-1)=fac2*.5
	    a1(i,i  )=fac1*.5
	    a1(i,i+1)=fac1*.5
	    a1(i,i+2)=fac2*.5
      enddo
      endif 
      A(1)=0
      B(1)=1.
      C(1)=1.
      RHS(1) =(1./4)*f(0)+(3./2)*f(1)+f(2)/4
      a1(1,1)=1./4
      a1(1,2)=3./2
      a1(1,3)=1./4
      a1(n,n+1)=1./4
      a1(n,n  )=3./2
      a1(n,n-1)=1./4
      C(N)=0
      B(N)=1
      A(N)=1
      RHS(n) =(1./4)*f(n)+(3./2)*f(n-1)+f(n-2)/4
      if (flag.eq.'m') then 
       do i=2,n-1
	  t1(i,i)=b(i)
          t1(i,i-1)=a(i)
          t1(i,i+1)=c(i)
       enddo
	  
      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n,n)=b(n)
      t1(n,n-1)=a(n)
      td =0
      call gaussj(t1,n,n,td,n,n)
      r=0
      do j=1,n
	     do k=1,n+1
         do i=1,n
	  r(j,k)=r(j,k)+t1(j,i)*a1(i,k)
         enddo
       enddo
      enddo
      endif
      call trid(N,A,B,C,RHS,D)
      do i=1,n
       df(i)=rhs(i)
      enddo
      end 
      
      subroutine interpol_c_s6_m(flag,N,NRHS,F,dF,dx,t1,a1,r)
      character*1 flag
      real f(0:N,nrhs),df(N,nrhs),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      integer N
      real A(N),B(N),C(N),D(N),RHS(N,NRHS)
      dxi=1./dx
      do i=2,n-1 
      A  (i  )=1./6
      B  (i  )=1.
      C  (i  )=1./6
      do j=1,nrhs
      RHS(i,j  )=2./3*(F(i,j)+F(i-1,j  ))
      enddo

      enddo
      alpha =.3
      fac1=(1./8)*(9+10*alpha)
      fac2=(1./8)*(6*alpha-1 )


      do i=3,n-2
        A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
        do j=1,nrhs
        RHS (i,j) = fac1*(F(i,j)+F(i-1,j))*.5 + fac2*(F(i+1,j)+F(i-2,j))*.5
        enddo
      enddo 
      if (flag.eq.'m') then
      do i=2,n-1
      a1(i,i )=2./3
      a1(i,i+1)=2./3
      enddo
      do i=3,n-2
	    a1(i,i-1)=fac2*.5
	    a1(i,i  )=fac1*.5
	    a1(i,i+1)=fac1*.5
	    a1(i,i+2)=fac2*.5
      enddo
      endif 
      A(1)=0
      B(1)=1.
      C(1)=1.
      do j=1,nrhs
      RHS(1,j) =(1./4)*f(0,j)+(3./2)*f(1,j)+f(2,j)/4
      enddo
      a1(1,1)=1./4
      a1(1,2)=3./2
      a1(1,3)=1./4
      a1(n,n+1)=1./4
      a1(n,n  )=3./2
      a1(n,n-1)=1./4
      C(N)=0
      B(N)=1
      A(N)=1
      do j=1,nrhs
      RHS(n,j) =(1./4)*f(n,j)+(3./2)*f(n-1,j)+f(n-2,j)/4
      enddo
      if (flag.eq.'m') then 
       do i=2,n-1
	  t1(i,i)=b(i)
          t1(i,i-1)=a(i)
          t1(i,i+1)=c(i)
       enddo
	  
      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n,n)=b(n)
      t1(n,n-1)=a(n)
      td =0
      call gaussj(t1,n,n,td,n,n)
      r=0
      do j=1,n
	     do k=1,n+1
         do i=1,n
	  r(j,k)=r(j,k)+t1(j,i)*a1(i,k)
         enddo
       enddo
      enddo
      endif
      call dgtsv(N,NRHS,A(2),B,C,RHS,N,info)
      do i=1,n
       do j=1,nrhs
       df(i,j)=rhs(i,j)
      enddo
      enddo
      end 
      
      subroutine der1_c_s6_m(flag,N,NRHS,F,dF,dx,t1,a1,r)
      real f(0:N,NRHS),df(N,NRHS),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      integer N
      character*1 flag
      real A(N),B(N),C(N),D(N),RHS(N,NRHS)
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      dxi=1./dx
      do i=2,n-1 
      A  (i  )=1./22
      B  (i  )=1.
      C  (i  )=1./22 
      do j=1,nrhs
      RHS(i ,j )=3*(3-1./11)/8.*(F(i,j)-F(i-1,j  ))*dxi
      enddo
      enddo
      do i=3,n-2
        A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
      do j=1,nrhs
        RHS (i,j) = fac1*(F(i,j)-F(i-1,j))*dxi + fac2*(F(i+1,j)-F(i-2,j))*dxi /3.
      enddo
      enddo

      A(1)=0
      B(1)=1.
      C(1)=-1.
      do j=1,nrhs
      RHS(1,j) = dxi*(-f(0,j)+2*f(1,j)-f(2,j))
      enddo
      C(N)=0
      B(N)=1
      A(N)=-1.
      do j=1,nrhs
      RHS(N,j)=dxi*((f(n,j)-2*f(n-1,j)+f(n-2,j)))
      enddo
      call dgtsv(N,NRHS,A(2),b,c,rhs,n,info)
      do i=1,n
      do j=1,nrhs
       df(i,j)=rhs(i,j)
      enddo
      enddo
      end 
      


      subroutine der1_s_c6(flag,N,F,dF,dx,t1,a1,r)
c
c
c    |       |       |        |         |
c    1       2       3        4         5
c      df(1) using values in between
c 
      implicit none
      integer i,n,j,k
      character flag
      real F(N),dF(0:N),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1),D(N+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      real t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      a1=0
      t1=0
      r=0
      dxi = 1./dx
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      do i=2,n
      A  (i  )=1./22
      B  (i  )=1
      C  (i  )=1./22
      RHS(i  )=3*(3-1./11)/8.*(F(i)-F(i-1  ))*dxi
      enddo
CCCCCCC   n-2 replaced by n-1
      do i=3,n-1
        A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
        RHS (i) = fac1*(F(i)-F(i-1))*dxi + fac2*(F(i+1)-F(i-2))*dxi /3.
      enddo
     
      if (flag .eq. 'm') then
      do i=2,n
      a1(i,i-1)=-3*(3-1./11)/8.*dxi
      a1(i,i  )= 3*(3-1./11)/8.*dxi
      enddo
      do i=3,n-1
      a1(i,i-2)=-fac2*dxi/3
      a1(i,i-1)=-fac1*dxi
      a1(i,i  )= fac1*dxi
      a1(i,i+1)= fac2*dxi/3.
      enddo
      endif

      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=23.
      ac =-(23**2+71)/24.
      bc = (7*23+47)/8.
      cc = (23-31)/8.
      RHS(1  )=(ac*F(1)+bc*F(2  )+cc*F(3  ))*dxi
      if (flag.eq. 'm') then
      a1(1,1)=ac*dxi
      a1(1,2)=bc*dxi
      a1(1,3)=cc*dxi
      endif
      A  (N+1  )=23
      B  (N+1  )=1.
      C  (N+1  )=0.
      RHS(N+1  )=(-ac*F(N)-bc*F(N-1)-cc*F(N-2))*dxi
      if (flag.eq. 'm') then
      a1(n+1,n    )=-ac*dxi
      a1(n+1,n-1  )=-bc*dxi
      a1(n+1,n-2  )=-cc*dxi

      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n+1,n+1)=b(n+1)
      t1(n+1,n  )=a(n+1)
      do i=2,n
    	t1(i,i-1)=a(i)
        t1(i,i  )=b(i)
        t1(i,i+1)=c(i)
      enddo
      td =0
      call gaussj(t1,n+1,n+1,td,n+1,n+1) 
      endif
      call trid(N+1,A,B,C,RHS,D)
      do i=1,N+1
       dF(i-1)=RHS(i)
      enddo
      if (flag.eq. 'm') then
      do i=1,n
       do j=1,n+1
	  r(j,i)=0
          do k=1,n+1
	    r(j,i)=r(j,i)+t1(j,k)*a1(k,i)
          enddo
       enddo
      enddo
      endif
      end
  
      subroutine der1_s_c6_m(flag,N,nrhs,F,dF,dx,t1,a1,r)
c
c
c    |       |       |        |         |
c    1       2       3        4         5
c      df(1) using values in between
c 
      implicit none
      integer i,n,j,k,info,nrhs
      character flag
      real F(N,nrhs),dF(0:N,nrhs),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1,nrhs),D(N+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      real t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      a1=0
      t1=0
      r=0
      dxi = 1./dx
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      do i=2,n
      A  (i  )=1./22
      B  (i  )=1
      C  (i  )=1./22
      do j=1,nrhs
      RHS(i ,j )=3*(3-1./11)/8.*(F(i,j)-F(i-1 ,j ))*dxi
      enddo
      enddo
CCCCCCC   n-2 replaced by n-1
      do i=3,n-1
        A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
        do j=1,nrhs
        RHS (i,j) = fac1*(F(i,j)-F(i-1,j))*dxi + fac2*(F(i+1,j)-F(i-2,j))*dxi /3.
        enddo
      enddo
     
      if (flag .eq. 'm') then
      do i=2,n
      a1(i,i-1)=-3*(3-1./11)/8.*dxi
      a1(i,i  )= 3*(3-1./11)/8.*dxi
      enddo
      do i=3,n-1
      a1(i,i-2)=-fac2*dxi/3
      a1(i,i-1)=-fac1*dxi
      a1(i,i  )= fac1*dxi
      a1(i,i+1)= fac2*dxi/3.
      enddo
      endif

      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=23.
      ac =-(23**2+71)/24.
      bc = (7*23+47)/8.
      cc = (23-31)/8.
      do j=1,nrhs
      RHS(1,j  )=(ac*F(1,j)+bc*F(2,j  )+cc*F(3,j  ))*dxi
      enddo
      if (flag.eq. 'm') then
      a1(1,1)=ac*dxi
      a1(1,2)=bc*dxi
      a1(1,3)=cc*dxi
      endif
      A  (N+1  )=23
      B  (N+1  )=1.
      C  (N+1  )=0.   
      do j=1,nrhs
      RHS(N+1 ,j )=(-ac*F(N,j)-bc*F(N-1,j)-cc*F(N-2,j))*dxi
      enddo
      if (flag.eq. 'm') then
      a1(n+1,n    )=-ac*dxi
      a1(n+1,n-1  )=-bc*dxi
      a1(n+1,n-2  )=-cc*dxi

      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n+1,n+1)=b(n+1)
      t1(n+1,n  )=a(n+1)
      do i=2,n
	    t1(i,i-1)=a(i)
        t1(i,i  )=b(i)
        t1(i,i+1)=c(i)
      enddo
      td =0
      call gaussj(t1,n+1,n+1,td,n+1,n+1) 
      endif
      call dgtsv(N+1,nrhs,A(2),B,C,RHS,N+1,info)
      do i=1,N+1
       do j=1,nrhs
       dF(i-1,j)=RHS(i,j)
      enddo
      enddo
      end
  
      subroutine derw_s_c(flag,N,F,dF,dx,t1,a1,r)
c
c
c    |       |       |        |         |
c    1       2       3        4         5
c      df(1) using values in between
c 
      implicit none
      integer i,n,j,k
      character flag
      real F(N),dF(0:N),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1),D(N+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      real t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      a1=0
      t1=0
      r=0
      dxi = 1./dx
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      do i=2,n
      A  (i  )=1./22
      B  (i  )=1
      C  (i  )=1./22
      RHS(i  )=3*(3-1./11)/8.*(F(i)-F(i-1  ))*dxi
      enddo
CCCCCCC   n-2 replaced by n-1
      do i=3,n-1
        A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
        RHS (i) = fac1*(F(i)-F(i-1))*dxi + fac2*(F(i+1)-F(i-2))*dxi /3.
      enddo
     

      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=23.
      ac =-(23**2+71)/24.
      bc = (7*23+47)/8.
      cc = (23-31)/8.
      RHS(1  )=(ac*F(1)+bc*F(2  )+cc*F(3  ))*dxi
      if (flag.eq. 'm') then
      a1(1,1)=ac*dxi
      a1(1,2)=bc*dxi
      a1(1,3)=cc*dxi
      endif
      A  (N+1  )=23
      B  (N+1  )=1.
      C  (N+1  )=0.
      RHS(N+1  )=(-ac*F(N)-bc*F(N-1)-cc*F(N-2))*dxi
      if (flag.eq. 'm') then
      a1(n+1,n    )=-ac*dxi
      a1(n+1,n-1  )=-bc*dxi
      a1(n+1,n-2  )=-cc*dxi

      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n+1,n+1)=b(n+1)
      t1(n+1,n  )=a(n+1)
      do i=2,n
	    t1(i,i-1)=a(i)
        t1(i,i  )=b(i)
        t1(i,i+1)=c(i)
      enddo
      td =0
      call gaussj(t1,n+1,n+1,td,n+1,n+1) 
      endif
      call trid(N+1,A,B,C,RHS,D)
      do i=1,N+1
       dF(i-1)=RHS(i)
      enddo
      if (flag.eq. 'm') then
      do i=1,n
       do j=1,n+1
	  r(j,i)=0
          do k=1,n+1
	    r(j,i)=r(j,i)+t1(j,k)*a1(k,i)
          enddo
       enddo
      enddo
      endif
      end
  
      subroutine interpol_s_c6(flag,N,f,df,dx,t1,a1,r)
      implicit none
      integer i,n,j,k
      character*1 flag
      real F(N),dF(0:N),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1),D(N+1),t1(n+1,n+1),a1(n+1,n),r(n+1,n),td(n+1,n+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      t1 =0.
      a1 =0
      r  =0
      do i=2,N
	    A   (i) = 1./6
        B   (i) = 1.
        C   (i) = 1./6
        RHS (i) = 2./3*(F(i)+F(i-1))
      enddo
      alpha =.3
      fac1=(1./8)*(9+10*alpha)
      fac2=(1./8)*(6*alpha-1 )
      do i=3,n-1
        A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
        RHS (i) = fac1*(F(i)+F(i-1))*.5 + fac2*(F(i+1)+F(i-2))*.5
      enddo

	 
      if (flag.eq.'m') then
      do i=2,n
      a1(i,i-1)=2./3
      a1(i,i  )=2./3
      enddo
      do i=3,n-1
      a1(i,i-2)=fac2*.5
      a1(i,i-1)=fac1*.5
      a1(i,i  )=fac1*.5
      a1(i,i+1)=fac2*.5
      enddo
      endif
      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=0
      ac =15./8
      bc =-5./4
      cc = 3./8
      if (flag.eq.'m') then
      a1(1,1)=ac
      a1(1,2)=bc
      a1(1,3)=cc
      endif
      RHS(1  )=(ac*F(1)+bc*F(2  )+cc*F(3  ))
      A  (N+1  )=0
      B  (N+1  )=1.
      C  (N+1  )=0.
      RHS(N+1  )=(ac*F(N)+bc*F(N-1)+cc*F(N-2))
      if (flag.eq.'m') then
      a1(n+1,n  )=ac
      a1(n+1,n-1)=bc
      a1(n+1,n-2)=cc

      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n+1,n+1)=b(n+1)
      t1(n+1,n  )=a(n+1)
      do i=2,n
	    t1(i,i-1)=a(i)
        t1(i,i  )=b(i)
        t1(i,i+1)=c(i)
      enddo
      td = 0
      call gaussj(t1,n+1,n+1,td,n+1,n+1)
      do i=1,n
       do j=1,n+1
	  r(j,i)=0
          do k=1,n+1
	    r(j,i)=r(j,i)+t1(j,k)*a1(k,i)
          enddo
       enddo
      enddo
      endif
      call trid(N+1,A,B,C,RHS,D)
      do i=1,N+1
       dF(i-1)=RHS(i)
      enddo
      end

      subroutine interpol_s_c6_m(flag,N,nrhs,f,df,dx,t1,a1,r)
      implicit none
      integer i,n,j,k,info,nrhs
      character*1 flag
      real F(N,nrhs),dF(0:N,nrhs),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1,nrhs),D(N+1),t1(n+1,n+1),a1(n+1,n),r(n+1,n),td(n+1,n+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      t1 =0.
      a1 =0
      r  =0
      do i=2,N
	    A   (i) = 1./6
        B   (i) = 1.
        C   (i) = 1./6
        do j=1,nrhs
        RHS (i,j) = 2./3*(F(i,j)+F(i-1,j))
        enddo
      enddo
      alpha =.3
      fac1=(1./8)*(9+10*alpha)
      fac2=(1./8)*(6*alpha-1 )
      do i=3,n-1
        A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
        do j=1,nrhs
        RHS (i,j) = fac1*(F(i,j)+F(i-1,j))*.5 + fac2*(F(i+1,j)+F(i-2,j))*.5
        enddo
      enddo

	 
      if (flag.eq.'m') then
      do i=2,n
      a1(i,i-1)=2./3
      a1(i,i  )=2./3
      enddo
      do i=3,n-1
      a1(i,i-2)=fac2*.5
      a1(i,i-1)=fac1*.5
      a1(i,i  )=fac1*.5
      a1(i,i+1)=fac2*.5
      enddo
      endif
      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=0
      ac =15./8
      bc =-5./4
      cc = 3./8
      if (flag.eq.'m') then
      a1(1,1)=ac
      a1(1,2)=bc
      a1(1,3)=cc
      endif
      do j=1,nrhs
      RHS(1,j  )=(ac*F(1,j)+bc*F(2,j  )+cc*F(3 ,j ))
      enddo
      A  (N+1  )=0
      B  (N+1  )=1.
      C  (N+1  )=0.
      do j=1,nrhs
      RHS(N+1 ,j )=(ac*F(N,j)+bc*F(N-1,j)+cc*F(N-2,j))
      enddo
      if (flag.eq.'m') then
      a1(n+1,n  )=ac
      a1(n+1,n-1)=bc
      a1(n+1,n-2)=cc

      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n+1,n+1)=b(n+1)
      t1(n+1,n  )=a(n+1)
      do i=2,n
	    t1(i,i-1)=a(i)
        t1(i,i  )=b(i)
        t1(i,i+1)=c(i)
      enddo
      td = 0
      call gaussj(t1,n+1,n+1,td,n+1,n+1)
      do i=1,n
       do j=1,n+1
	  r(j,i)=0
          do k=1,n+1
	    r(j,i)=r(j,i)+t1(j,k)*a1(k,i)
          enddo
       enddo
      enddo
      endif
      call dgtsv(N+1,nrhs,A(2),B,C,RHS,n+1,info)
      do i=1,N+1
       do j=1,nrhs
       dF(i-1,j)=RHS(i,j)
       enddo
      enddo
      end

      subroutine interpol_c_s4(flag,N,F,dF,dx,t1,a1,r)
      character*1 flag
      real f(0:N),df(N),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      integer N
      real A(N),B(N),C(N),D(N),RHS(N)
      dxi=1./dx
      do i=2,n-1 
      A  (i  )=1./6
      B  (i  )=1.
      C  (i  )=1./6
      RHS(i  )=2./3*(F(i)+F(i-1  ))
      enddo
      if (flag.eq.'m') then
      do i=2,n-1
      a1(i,i )=2./3
      a1(i,i+1)=2./3
      enddo
      endif 
      A(1)=0
      B(1)=1.
      C(1)=1.
      RHS(1) =(1./4)*f(0)+(3./2)*f(1)+f(2)/4
      a1(1,1)=1./4
      a1(1,2)=3./2
      a1(1,3)=1./4
      a1(n,n+1)=1./4
      a1(n,n  )=3./2
      a1(n,n-1)=1./4
      C(N)=0
      B(N)=1
      A(N)=1
      RHS(n) =(1./4)*f(n)+(3./2)*f(n-1)+f(n-2)/4
      if (flag.eq.'m') then 
       do i=2,n-1
	  t1(i,i)=b(i)
          t1(i,i-1)=a(i)
          t1(i,i+1)=c(i)
       enddo
      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n,n)=b(n)
      t1(n,n-1)=a(n)
      td =0
      call gaussj(t1,n,n,td,n,n)
      r=0
      do j=1,n
	     do k=1,n+1
         do i=1,n
	  r(j,k)=r(j,k)+t1(j,i)*a1(i,k)
         enddo
       enddo
      enddo
      endif
      call trid(N,A,B,C,RHS,D)
      do i=1,n
       df(i)=rhs(i)
      enddo
      end 
      
      subroutine der1_c_s4(flag,N,F,dF,dx,t1,a1,r)
      real f(0:N),df(N),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      integer N
      character*1 flag
      real A(N),B(N),C(N),D(N),RHS(N)
      dxi=1./dx
      do i=2,n-1 
      A  (i  )=1./22
      B  (i  )=1.
      C  (i  )=1./22 
      RHS(i  )=3*(3-1./11)/8.*(F(i)-F(i-1  ))*dxi
      enddo
      if (flag.eq.'m') then
      do i=2,n-1
      a1(i,i )=-3*(3-1./11)/8.*dxi
      a1(i,i+1)=3*(3-1./11)/8.*dxi
      enddo
      endif 
      A(1)=0
      B(1)=1.
      C(1)=-1.
      RHS(1) = dxi*(-f(0)+2*f(1)-f(2))
      C(N)=0
      B(N)=1
      A(N)=-1.
      RHS(N)=dxi*((f(n)-2*f(n-1)+f(n-2)))
      if (flag.eq.'m') then
      a1(1,1)=-dxi
      a1(1,2)=2*dxi
      a1(1,3)=-dxi
      a1(n,n+1)=dxi
      a1(n,n  )=-2*dxi
      a1(n,n-1)=dxi
       do i=2,n-1
	  t1(i,i)  =b(i)
          t1(i,i-1)=a(i)
          t1(i,i+1)=c(i)
       enddo
      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n,n)=b(n)
      t1(n,n-1)=a(n)
      td=0
      call gaussj(t1,n,n,td,n,n)
      r=0
      do j=1,n
    	 do k=1,n+1
         do i=1,n
	  r(j,k)=r(j,k)+t1(j,i)*a1(i,k)
         enddo
       enddo
      enddo
      endif
      call trid(N,A,B,C,RHS,D)
      do i=1,n
       df(i)=rhs(i)
      enddo
      end 
      


      subroutine der1_s_c4(flag,N,F,dF,dx,t1,a1,r)
c
c
c    |       |       |        |         |
c    1       2       3        4         5
c      df(1) using values in between
c 
      implicit none
      integer i,n,j,k
      character flag
      real F(N),dF(0:N),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1),D(N+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      real t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      a1=0
      t1=0
      r=0
      dxi = 1./dx
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      do i=2,n
      A  (i  )=1./22
      B  (i  )=1
      C  (i  )=1./22
      RHS(i  )=3*(3-1./11)/8.*(F(i)-F(i-1  ))*dxi
      enddo
     
      if (flag .eq. 'm') then
      do i=2,n
      a1(i,i-1)=-3*(3-1./11)/8.*dxi
      a1(i,i  )= 3*(3-1./11)/8.*dxi
      enddo
      endif

      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=23.
      ac =-(23**2+71)/24.
      bc = (7*23+47)/8.
      cc = (23-31)/8.
      RHS(1  )=(ac*F(1)+bc*F(2  )+cc*F(3  ))*dxi
      if (flag.eq. 'm') then
      a1(1,1)=ac*dxi
      a1(1,2)=bc*dxi
      a1(1,3)=cc*dxi
      endif
      A  (N+1  )=23
      B  (N+1  )=1.
      C  (N+1  )=0.
      RHS(N+1  )=(-ac*F(N)-bc*F(N-1)-cc*F(N-2))*dxi
      if (flag.eq. 'm') then
      a1(n+1,n    )=-ac*dxi
      a1(n+1,n-1  )=-bc*dxi
      a1(n+1,n-2  )=-cc*dxi

      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n+1,n+1)=b(n+1)
      t1(n+1,n  )=a(n+1)
      do i=2,n
	t1(i,i-1)=a(i)
        t1(i,i  )=b(i)
        t1(i,i+1)=c(i)
      enddo
      td =0
      call gaussj(t1,n+1,n+1,td,n+1,n+1) 
      endif
      call trid(N+1,A,B,C,RHS,D)
      do i=1,N+1
       dF(i-1)=RHS(i)
      enddo
      if (flag.eq. 'm') then
      do i=1,n
       do j=1,n+1
	  r(j,i)=0
          do k=1,n+1
	    r(j,i)=r(j,i)+t1(j,k)*a1(k,i)
          enddo
       enddo
      enddo
      endif
      end
  

      subroutine interpol_s_c4(flag,N,f,df,dx,t1,a1,r)
      implicit none
      integer i,n,j,k
      character*1 flag
      real F(N),dF(0:N),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1),D(N+1),t1(n+1,n+1),a1(n+1,n),r(n+1,n),td(n+1,n+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      alpha =1./6
      fac1=2./3.
      t1 =0.
      a1 =0
      r  =0
      do i=2,N
	    A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
        RHS (i) = fac1*(F(i)+F(i-1))
      enddo
      if (flag.eq.'m') then
      do i=2,n
      a1(i,i-1)=fac1
      a1(i,i  )=fac1
      enddo
      endif
      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=0
      ac =15./8
      bc =-5./4
      cc = 3./8
      if (flag.eq.'m') then
      a1(1,1)=ac
      a1(1,2)=bc
      a1(1,3)=cc
      endif
      RHS(1  )=(ac*F(1)+bc*F(2  )+cc*F(3  ))
      A  (N+1  )=0
      B  (N+1  )=1.
      C  (N+1  )=0.
      RHS(N+1  )=(ac*F(N)+bc*F(N-1)+cc*F(N-2))
      if (flag.eq.'m') then
      a1(n+1,n  )=ac
      a1(n+1,n-1)=bc
      a1(n+1,n-2)=cc

      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n+1,n+1)=b(n+1)
      t1(n+1,n  )=a(n+1)
      do i=2,n
	    t1(i,i-1)=a(i)
        t1(i,i  )=b(i)
        t1(i,i+1)=c(i)
      enddo
      td = 0
      call gaussj(t1,n+1,n+1,td,n+1,n+1)
      do i=1,n
       do j=1,n+1
	  r(j,i)=0
          do k=1,n+1
	    r(j,i)=r(j,i)+t1(j,k)*a1(k,i)
          enddo
       enddo
      enddo
      endif
      call trid(N+1,A,B,C,RHS,D)
      do i=1,N+1
       dF(i-1)=RHS(i)
      enddo
      end

      SUBROUTINE gaussj(a,n,np,b,m,mp)
      INTEGER m,mp,n,np,NMAX
      REAL a(np,np),b(np,mp)
      PARAMETER (NMAX=900)
      INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      REAL big,dum,pivinv
      do 11 j=1,n
        ipiv(j)=0
11    continue
      do 22 i=1,n
        big=0.
        do 13 j=1,n
          if(ipiv(j).ne.1)then
            do 12 k=1,n
              if (ipiv(k).eq.0) then
                if (abs(a(j,k)).ge.big)then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if (ipiv(k).gt.1) then
                pause 'singular matrix in gaussj'
              endif
12          continue
          endif
13      continue
        ipiv(icol)=ipiv(icol)+1
        if (irow.ne.icol) then
          do 14 l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
14        continue
          do 15 l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
15        continue
        endif
        indxr(i)=irow
        indxc(i)=icol
        if (a(icol,icol).eq.0.) pause 'singular matrix in gaussj'
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do 16 l=1,n
          a(icol,l)=a(icol,l)*pivinv
16      continue
        do 17 l=1,m
          b(icol,l)=b(icol,l)*pivinv
17      continue
        do 21 ll=1,n
          if(ll.ne.icol)then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do 18 l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
18          continue
            do 19 l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
19          continue
          endif
21      continue
22    continue
      do 24 l=n,1,-1
        if(indxr(l).ne.indxc(l))then
          do 23 k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
23        continue
        endif
24    continue
      return
      END


      subroutine der1w_s_c6(N,F,dF,dx)
c
c
c    |       |       |        |         |
c    1       2       3        4         5
c      df(1) using values in between
c 
      implicit none
      integer i,n,j,k
      character flag
      real F(N),dF(0:N),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1),D(N+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      real t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      a1=0
      t1=0
      r=0
      dxi = 1./dx
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      do i=2,n
      A  (i  )=1./22
      B  (i  )=1
      C  (i  )=1./22
      RHS(i  )=3*(3-1./11)/8.*(F(i)-F(i-1  ))*dxi
      enddo
      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=15!23!15.
      ac =-(23**2+71)/24.
      bc = (7*23+47)/8.
      cc = (23-31)/8.
c      RHS(1  )=(ac*F(1)+bc*F(2  )+cc*F(3  ))*dxi
      RHS(1  )=(-15*F(1)+50*F(2  )/3-0.6*F(3 ))*dxi
      A  (N+1  )=15
      B  (N+1  )=1.
      C  (N+1  )=0.
      RHS(N+1  )=(15*f(n)-50*f(n-1)/3+0.6*f(n-2))*dxi
      call trid(N+1,A,B,C,RHS,D)
      do i=1,N+1
       dF(i-1)=RHS(i)
      enddo
      end


      subroutine der1c_s_c6(N,F,dF,dx,Twh,Twc)
c
c
c    |       |       |        |         |
c    1       2       3        4         5
c      df(1) using values in between
c 
      implicit none
      integer i,n,j,k
      character flag
      real F(N),dF(0:N),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1),D(N+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      real t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1),Twh,Twc
      a1=0
      t1=0
      r=0
      dxi = 1./dx
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      do i=2,n
      A  (i  )=1./22
      B  (i  )=1
      C  (i  )=1./22

      RHS(i  )=3*(3-1./11)/8.*(F(i)-F(i-1  ))*dxi
      enddo
      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=15.
      ac =-(23**2+71)/24.
      bc = (7*23+47)/8.
      cc = (23-31)/8.
c     RHS(1  )=(ac*F(1)+bc*F(2  )+cc*F(3  ))*dxi
      RHS(1  )=(-15*F(1)+50*F(2  )/3-0.6*F(3 )-16./15.*Twh)*dxi
      A  (N+1  )=15
      B  (N+1  )=1.
      C  (N+1  )=0.
      RHS(N+1  )=(15*f(n)-50*f(n-1)/3+0.6*f(n-2)+16./15.*Twc)*dxi
      call trid(N+1,A,B,C,RHS,D)
      do i=1,N+1
       dF(i-1)=RHS(i)
      enddo
      end
  

      subroutine der1c_s_c6_m(N,NRHS,F,dF,dx,Twh,Twc)
c
c
c    |       |       |        |         |
c    1       2       3        4         5
c      df(1) using values in between
c 
      implicit none
      integer i,n,j,k,info,nrhs
      character flag
      real F(N,NRHS),dF(0:N,NRHS),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1,NRHS),D(N+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      real t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      real Qval,mr_c(0:N+1),Twh,Twc
      a1=0
      t1=0
      r=0
      dxi = 1./dx
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      do i=2,n
      A  (i  )=1./22
      B  (i  )=1
      C  (i  )=1./22
      do j=1,nrhs
      RHS(i,j  )=3*(3-1./11)/8.*(F(i,j)-F(i-1,j  ))*dxi
      enddo
      enddo
      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=15!23.
      ac =-(23**2+71)/24.
      bc = (7*23+47)/8.
      cc = (23-31)/8.
      do j=1,nrhs
c      RHS(1,j  )=(ac*F(1,j)+bc*F(2,j)+cc*F(3,j))*dxi
      RHS(1,j  )=(-15*F(1,j)+50*F(2,j)/3-0.6*F(3,j)-16./15.*Twh)*dxi
      enddo
      A  (N+1  )=15.
      B  (N+1  )=1.
      C  (N+1  )=0.
      do j=1,nrhs
       RHS(N+1,j  )=(15*f(n,j)-50*f(n-1,j)/3+0.6*f(n-2,j)+16./15.*Twc)*dxi
      enddo

c       A  (1)=0
c       B  (1)=1
c       C  (1)=0
c       A  (N+1)=0
c       B  (N+1)=1
c       C  (N+1)=0
c       do j=1,nrhs
c       RHS(1,j)= -Qval/mr_c(0)*1
c       RHS(N+1,j)= Qval/mr_c(N)*1
c       enddo

c      call trid(N+1,A,B,C,RHS,D)
      call dgtsv(N+1,nrhs,A(2),B,C,RHS,n+1,info)
      do i=1,N+1
       do j=1,nrhs
       dF(i-1,j)=RHS(i,j)
      enddo
      enddo
      end
     
      subroutine der1mu_s_c6_m(N,NRHS,F,dF,dx)
c
c
c    |       |       |        |         |
c    1       2       3        4         5
c      df(1) using values in between
c 
      implicit none
      integer i,n,j,k,info,nrhs
      character flag
      real F(N,NRHS),dF(0:N,NRHS),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1,NRHS),D(N+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      real t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      real Qval,mr_c(0:N+1)
      a1=0
      t1=0
      r=0
      dxi = 1./dx
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      do i=2,n
      A  (i  )=1./22
      B  (i  )=1
      C  (i  )=1./22
      do j=1,nrhs
      RHS(i,j  )=3*(3-1./11)/8.*(F(i,j)-F(i-1,j  ))*dxi
      enddo
      enddo
      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=15.
c      ac =-(23**2+71)/24.
c      bc = (7*23+47)/8.
c      cc = (23-31)/8.
      do j=1,nrhs
      RHS(1,j  )=(-15*F(1,j)+50*F(2,j)/3-0.6*F(3,j)-16./15.0*0.1/40)*dxi
      enddo
      A  (N+1  )=15.
      B  (N+1  )=1.
      C  (N+1  )=0.
      do j=1,nrhs
       RHS(N+1,j  )=(15*f(n,j)-50*f(n-1,j)/3+0.6*f(n-2,j)+16.0/15.0*1./40.)*dxi
      enddo

c      call trid(N+1,A,B,C,RHS,D)
      call dgtsv(N+1,nrhs,A(2),B,C,RHS,n+1,info)
      do i=1,N+1
       do j=1,nrhs
       dF(i-1,j)=RHS(i,j)
      enddo
      enddo
      
      end


      subroutine der1w_s_c6_m(N,NRHS,F,dF,dx)
c
c
c    |       |       |        |         |
c    1       2       3        4         5
c      df(1) using values in between
c 
      implicit none
      integer i,n,j,k,info,nrhs
      character flag
      real F(N,NRHS),dF(0:N,NRHS),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1,NRHS),D(N+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      real t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      a1=0
      t1=0
      r=0
      dxi = 1./dx
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      do i=2,n
      A  (i  )=1./22
      B  (i  )=1
      C  (i  )=1./22
      do j=1,nrhs
      RHS(i,j  )=3*(3-1./11)/8.*(F(i,j)-F(i-1,j  ))*dxi
      enddo
      enddo
      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=15!23.!15.
      ac =-(23**2+71)/24.
      bc = (7*23+47)/8.
      cc = (23-31)/8.
      do j=1,nrhs
c      RHS(1 ,j )=(ac*F(1,j)+bc*F(2,j  )+cc*F(3,j  ))*dxi
      RHS(1,j  )=(-15*F(1,j)+50*F(2,j)/3-0.6*F(3,j))*dxi
      enddo
      A  (N+1  )=15
      B  (N+1  )=1.
      C  (N+1  )=0.
      do j=1,nrhs
      RHS(N+1,j  )=(15*f(n,j)-50*f(n-1,j)/3+0.6*f(n-2,j))*dxi
      enddo
c      call trid(N+1,A,B,C,RHS,D)
      call dgtsv(N+1,nrhs,A(2),B,C,RHS,n+1,info)
      do i=1,N+1
       do j=1,nrhs
       dF(i-1,j)=RHS(i,j)
      enddo
      enddo
      end
  
      subroutine der1l_s_c6_m(N,NRHS,F,dF,dx)
c
c
c    |       |       |        |         |
c    1       2       3        4         5
c      df(1) using values in between
c 
      implicit none
      integer i,n,j,k,info,nrhs
      character flag
      real F(N,NRHS),dF(0:N,NRHS),dx 
      real A(N+1),B(N+1),C(N+1),RHS(N+1,NRHS),D(N+1)
      real dxi,ac,bc,cc,alpha,fac1,fac2
      real t1(n+1,n+1) ,a1(n+1,n),r(n+1,n),td(n+1,n+1)
      real Qval,mr_c(0:N+1)
      a1=0
      t1=0
      r=0
      dxi = 1./dx
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      do i=2,n
      A  (i  )=1./22
      B  (i  )=1
      C  (i  )=1./22
      do j=1,nrhs
      RHS(i,j  )=3*(3-1./11)/8.*(F(i,j)-F(i-1,j  ))*dxi
      enddo
      enddo
      A  (1  )=0.
      B  (1  )=1.
      C  (1  )=15.
c      ac =-(23**2+71)/24.
c      bc = (7*23+47)/8.
c      cc = (23-31)/8.
      do j=1,nrhs
      RHS(1,j  )=(-15*F(1,j)+50*F(2,j)/3-0.6*F(3,j)-16./15.0*1 )*dxi
      enddo
      A  (N+1  )=15.
      B  (N+1  )=1.
      C  (N+1  )=0.
      do j=1,nrhs
       RHS(N+1,j  )=(15*f(n,j)-50*f(n-1,j)/3+0.6*f(n-2,j)+16.0/15.0*1)*dxi
      enddo

c      call trid(N+1,A,B,C,RHS,D)
      call dgtsv(N+1,nrhs,A(2),B,C,RHS,n+1,info)
      do i=1,N+1
       do j=1,nrhs
       dF(i-1,j)=RHS(i,j)
      enddo
      enddo
      end 


      subroutine der1_c_s6(flag,N,F,dF,dx,t1,a1,r)
      real f(0:N),df(N),dx,t1(n,n),a1(n,n+1),r(n,n+1),td(n,n)
      integer N
      character*1 flag
      real A(N),B(N),C(N),D(N),RHS(N)
      alpha= 9./62
      fac1 = (3./8)*(3-2*alpha)
      fac2 = (1./8)*(-1+22*alpha)
      dxi=1./dx
      i=2 
      A  (i  )=1./22
      B  (i  )=1.
      C  (i  )=1./22 
      RHS(i  )=3*(3-1./11)/8.*(F(i)-F(i-1  ))*dxi
      i=n-1
      A  (i  )=1./22
      B  (i  )=1.
      C  (i  )=1./22 
      RHS(i  )=3*(3-1./11)/8.*(F(i)-F(i-1  ))*dxi
      do i=3,n-2
        A   (i) = alpha
        B   (i) = 1.
        C   (i) = alpha
        RHS (i) = fac1*(F(i)-F(i-1))*dxi + fac2*(F(i+1)-F(i-2))*dxi /3.
      enddo

      if (flag.eq.'m') then
      do i=2,n-1
      a1(i,i )=-3*(3-1./11)/8.*dxi
      a1(i,i+1)=3*(3-1./11)/8.*dxi
      enddo

      do i=3,n-2
      a1(i,i-1)=-fac2*dxi/3.
      a1(i,i  )=-fac1*dxi
      a1(i,i+1)=fac1*dxi
      a1(i,i+2)=fac2*dxi/3.
      enddo
      endif 
      A(1)=0
      B(1)=1.
      C(1)=-1.
      RHS(1) = dxi*(-f(0)+2*f(1)-f(2))
      C(N)=0
      B(N)=1
      A(N)=-1.
      RHS(N)=dxi*((f(n)-2*f(n-1)+f(n-2)))
      if (flag.eq.'m') then
      a1(1,1)=-dxi
      a1(1,2)=2*dxi
      a1(1,3)=-dxi
      a1(n,n+1)=dxi
      a1(n,n  )=-2*dxi
      a1(n,n-1)=dxi
       do i=2,n-1
	  t1(i,i)  =b(i)
          t1(i,i-1)=a(i)
          t1(i,i+1)=c(i)
       enddo
      t1(1,1)=b(1)
      t1(1,2)=c(1)
      t1(n,n)=b(n)
      t1(n,n-1)=a(n)
      td=0
      call gaussj(t1,n,n,td,n,n)
      r=0
      do j=1,n
	 do k=1,n+1
         do i=1,n
	  r(j,k)=r(j,k)+t1(j,i)*a1(i,k)
         enddo
       enddo
      enddo
      endif
      call trid(N,A,B,C,RHS,D)
      do i=1,n
       df(i)=rhs(i)
      enddo
      end 
      
