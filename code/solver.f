
      subroutine solver_com(ini,p,ru,rp,dr,dtheta,dz,mr_s,mr_c,rank)
      use decomp_2d


      implicit none
      include 'param.txt'
      include 'mpif.h' 
      character*1 even,odd
      real     p(imax,jmax/p_row,kmax/p_col),ru(0:i1),rp(0:i1),dtheta,dz
      real     pp (0:i1,jmax,kmax)
      real    wj(2*jmax+15),trt(jmax),pi,pdi(imax)
      real    wk(2*kmax+15),zrt(kmax),t_som
      real    a(imax),b(imax),c(imax),bb(imax),d(imax),dd(imax),rhs(imax),x(imax),rhs2(imax,2)
      real    ft(jmax),dft(jmax),fz(kmax),dfz(kmax),som,stime,dr
      real    hi_c(0:imax,jmax/p_row,kmax/p_col),hi_s(imax,jmax/p_row,kmax/p_col),resid,adiag(imax)
c     NOTE: MATSTOR DIM REDUCED ONLY IF 2ND ORDER PRECISION IS USED
c      real    matstor(imax,imax,(jmax/p_row)/2+1,(kmax/p_col)/2+1),mr_s(imax),mr_c(0:imax)
      real    matstor(1,1,1,1),mr_s(imax),mr_c(0:imax)
      integer istor(imax,jmax/p_row,kmax/p_col),jdex,kdex,j_dex,k_dex,i_dex
      integer ier,mxmv,indx(imax),kk,info,ini,ii,indj(jmax),indk(kmax),ktel,jtel
   

c
c   
      real matx(imax,imax),mat(imax,imax)
      real t1x(imax+1,imax+1),a1x(imax+1,imax),r1x(imax+1,imax)
      real t2x(imax  ,imax  ),a2x(imax,imax+1),r2x(imax,imax+1)

      real p1x(0:imax+1,jmax/p_row,kmax/p_col)
      real p1y(0:imx,   jmax,kmax/p_col)
      real p1z(0:imx,   jmax/p_col,kmax)


      common /mstor/matstor
      common /mstoi/istor
c solves poisson equation in cylindrical coordinates, using FFT's' in theta and z
c   and staggered compact finite difference in the radial direction.
c   b.j.boersma, 2010
c
c   modified 17-3-2012  with odd/even storage of L and U 
c                       proper values of jtel and ktel  
c
      common /mats/matx,adiag
      stime = MPI_WTIME()
      pi =    4.*atan(1.)

cSet  Setup lookup tables.

      call    vrffti(jmax,wj)
      call    vrffti(kmax,wk)

      do  j=1,jmax
         trt(j) =-( (j/2)*1.)**2   * (8*atan(1.)/(jmax*dtheta   ))**2 
      enddo
      trt(1) = 0!-1e-9
      do k=1,kmax
        zrt(k)= - ( (k/2)*1.)**2   * (8*atan(1.)/(kmax*dz   ))**2
      enddo
      zrt(1) = 0!-1e-6
      if (ini.eq.0) then
      t1x = 0
      a1x = 0
      r1x = 0
      t2x = 0
      a2x = 0
      r2x = 0 
      matx  =0 
      call der1_s_c6('m',imax,hi_s,hi_c,dr,t1x,a1x,r1x)
      call der1_c_s6('m',imax,hi_c,hi_s,dr,t2x,a2x,r2x)
      do j=1,imax
        do i=1,imax+1
	 r1x(i,j)=r1x(i,j)*ru(i-1)
        enddo
      enddo
      do i=1,imax
	r1x(1,i)=0
        r1x(imax+1,i)=0
      enddo
      do j=1,imax
           do k=1,imax
             do i=1,imax+1
            matx(j,k)=matx(j,k)+r2x(j,i)*r1x(i,k)*mr_s(j)*mr_c(i-1)
            enddo
          enddo
        enddo
      do j=1,imax
	do i=1,imax
	  matx(i,j)=matx(i,j)/rp(i)
        enddo
      enddo
      endif


      do k=1,kmax/p_col
	do j=1,jmax/p_row
	  do i=1,imax
	   p1x(i,j,k)=p(i,j,k)
          enddo
         enddo
      enddo
      call transpose_x_to_y(p1x,p1y)
      do k=1,kmax/p_col
        do i=0,imx
          do j=1,jmax
           ft(j)=p1y(i,j,k)
           enddo
           call vrfftf(1,jmax,ft,dft,1,wj)
           do j=1,jmax
         p1y(i,j,k)=ft(j)
           enddo
         enddo
      enddo
      call transpose_y_to_z(p1y,p1z)
      do i=0,imx
        do j=1,jmax/p_col
          do k=1,kmax
          fz(k)=p1z(i,j,k)
           enddo
           call vrfftf(1,kmax,fz,dfz,1,wk)
          do k=1,kmax
          p1z(i,j,k)=fz(k)
           enddo
         enddo
      enddo
      call transpose_z_to_y(p1z,p1y)
      call transpose_y_to_x(p1y,p1x)
c********************************************************************
       
      ktel  = ( nrank - (nrank/p_col)*p_col) * (kmax/p_col)
      jtel  = (  (nrank/p_col)) * (jmax/p_row)
      if (mod(1+jtel,2).eq.0) even  = 'n'
      if (mod(1+jtel,2).eq.1) even  = 'y'
      if (ini.eq.0) then

       
        do k=1,kmax/p_col
        do j=1,jmax/p_row
         mat = matx
          do i=1,imax
          mat(i,i) = mat(i,i)+zrt(k+ ktel)+trt(j+jtel)/rp(i)**2
          enddo
c          if ((k+ktel).eq.1.and.(j+jtel).eq.1) mat(1, 1) = mat(1,1)*(1.0+1E-2)
             call dgetrf   (      imax,imax,mat,imax,indx,info)
             if (even.eq. 'y') then
             do ii=1,imax
              do i=1,imax
	        matstor(i,ii,(j)/2+1,(k)/2+1)=mat(i,ii)
              enddo
                  istor(  ii,(j),(k)/2+1)=indx(ii)
             enddo
             endif
             if (even.eq. 'n') then
             do ii=1,imax
              do i=1,imax
	        matstor(i,ii,(j+1)/2,(k)/2+1)=mat(i,ii)
              enddo
                  istor(  ii,(j),(k)/2+1)=indx(ii)
             enddo
            endif
        enddo
        enddo
      endif

C

        do k=1,1
         do  j =1,jmax/p_row
         do i=1,imax
	 rhs(i)   = p1x(i,j,k)
         enddo

             if (even.eq.'y') then
             do ii=1,imax
              do i=1,imax
	        mat(i,ii)=matstor(i,ii,(j)/2+1,(k)/2+1)
              enddo
                indx(ii)=istor(ii,(j),(k)/2+1)
             enddo
          endif
             if (even.eq.'n') then
             do ii=1,imax
              do i=1,imax
	        mat(i,ii)=matstor(i,ii,(j+1)/2,(k)/2+1)
              enddo
                indx(ii)=istor(ii,(j),(k)/2+1)
             enddo
          endif

          call dgetrs   (  'N',imax,1,mat,imax,indx,rhs,imax,info)
          do i=1,imax
	  p1x(i,j  ,k)=rhs(i)
          enddo
        enddo
       enddo

        do k=2,kmax/p_col-2,2
         do  j =1,jmax/p_row
         do i=1,imax
	 rhs2(i,1)   = p1x(i,j,k)
	 rhs2(i,2)   = p1x(i,j,k+1)
         enddo
             if (even.eq.'y') then
             do ii=1,imax
              do i=1,imax
	        mat(i,ii)=matstor(i,ii,(j)/2+1,(k)/2+1)
              enddo
                indx(ii)=istor(ii,(j),(k)/2+1)
             enddo
          endif
             if (even.eq.'n') then
             do ii=1,imax
              do i=1,imax
	        mat(i,ii)=matstor(i,ii,(j+1)/2,(k)/2+1)
              enddo
                indx(ii)=istor(ii,(j),(k)/2+1)
             enddo
          endif
          call dgetrs   (  'N',imax,2,mat,imax,indx,rhs2,imax,info)
          do i=1,imax
	  p1x(i,j  ,k)=rhs2(i,1)
	  p1x(i,j  ,k+1)=rhs2(i,2)
          enddo
        enddo
      enddo
        do k=kmax/p_col,kmax/p_col
         do  j =1,jmax/p_row
         do i=1,imax
	 rhs(i)   = p1x(i,j,k)
         enddo
             if (even.eq.'y') then
             do ii=1,imax
              do i=1,imax
	        mat(i,ii)=matstor(i,ii,(j)/2+1,(k)/2+1)
              enddo
                indx(ii)=istor(ii,(j),(k)/2+1)
             enddo
          endif
             if (even.eq.'n') then
             do ii=1,imax
              do i=1,imax
	        mat(i,ii)=matstor(i,ii,(j+1)/2,(k)/2+1)
              enddo
                indx(ii)=istor(ii,(j),(k)/2+1)
             enddo
           endif
          call dgetrs   (  'N',imax,1,mat,imax,indx,rhs,imax,info)
          do i=1,imax
	  p1x(i,j  ,k)=rhs(i)
          enddo
        enddo
       enddo


c********************************************************************
      call transpose_x_to_y(p1x,p1y)
      do k=1,kmax/p_col
        do i=0,imx
          do j=1,jmax
           ft(j)=p1y(i,j,k)
           enddo
           call vrfftb(1,jmax,ft,dft,1,wj)
           do j=1,jmax
         p1y(i,j,k)=ft(j)
           enddo
         enddo
      enddo
      call transpose_y_to_z(p1y,p1z)
      do i=0,imx
        do j=1,jmax/p_col
          do k=1,kmax
          fz(k)=p1z(i,j,k)
           enddo
           call vrfftb(1,kmax,fz,dfz,1,wk)
          do k=1,kmax
          p1z(i,j,k)=fz(k)
           enddo
         enddo
      enddo
      call transpose_z_to_y(p1z,p1y)
      call transpose_y_to_x(p1y,p1x)
      do k=1,kmax/p_col
        do j=1,jmax/p_row
          do i=1,imax
            p(i,j,k)=p1x(i,j,k)
          enddo
        enddo
      enddo
 
      end



      subroutine solver_2(ini,p,ru,rp,dr,dtheta,dz,mr_s,mr_c,rank)
      use decomp_2d


      implicit none
      include 'param.txt'
      include 'mpif.h' 
      character*1 even,odd
      real     p(imax,jmax/p_row,kmax/p_col),ru(0:i1),rp(0:i1),dtheta,dz
      real     pp (0:i1,jmax,kmax)
      real    wj(2*jmax+15),trt(jmax),pi,pdi(imax)
      real    wk(2*kmax+15),zrt(kmax),t_som
      real    a(imax),b(imax),c(imax),bb(imax),d(imax),dd(imax),rhs(imax),x(imax),rhs2(imax,2)
      real    ft(jmax),dft(jmax),fz(kmax),dfz(kmax),som,stime,dr
      real    hi_c(0:imax,jmax/p_row,kmax/p_col),hi_s(imax,jmax/p_row,kmax/p_col),resid,adiag(imax)
      real    matstor(imax,imax,(jmax/p_row)/2+1,(kmax/p_col)/2+1),mr_s(imax),mr_c(0:imax)
      integer istor(imax,jmax/p_row,kmax/p_col),jdex,kdex,j_dex,k_dex,i_dex
      integer ier,mxmv,indx(imax),kk,info,ini,ii,indj(jmax),indk(kmax),ktel,jtel
   

c
c   
      real matx(imax,imax),mat(imax,imax)
      real t1x(imax+1,imax+1),a1x(imax+1,imax),r1x(imax+1,imax)
      real t2x(imax  ,imax  ),a2x(imax,imax+1),r2x(imax,imax+1)

      real p1x(0:imax+1,jmax/p_row,kmax/p_col)
      real p1y(0:imx,   jmax,kmax/p_col)
      real p1z(0:imx,   jmax/p_col,kmax)


      common /mstor/matstor
      common /mstoi/istor
c solves poisson equation in cartesian coordinates, using FFT's' in theta and z
c   and staggered compact finite difference in the radial direction.
c   b.j.boersma, 2010
c
c   modified 17-3-2012  with odd/even storage of L and U 
c                       proper values of jtel and ktel  
c
      common /mats/matx,adiag
      stime = MPI_WTIME()
      pi =    4.*atan(1.)

cSet  Setup lookup tables.

      call    vrffti(jmax,wj)
      call    vrffti(kmax,wk)

      do  j=1,jmax
         trt(j) =-( (j/2)*1.)**2   * (8*atan(1.)/(jmax*dtheta   ))**2 
      enddo
      trt(1) = 0!-1e-5
      do k=1,kmax
        zrt(k)= - ( (k/2)*1.)**2   * (8*atan(1.)/(kmax*dz   ))**2
      enddo
      zrt(1) = 0!-1e-6

      do i=1,imax
       c(i)=1./(Rp(i+1)-Rp(i)) * (1. /((ru(i)-ru(i-1))))
       a(i)=1./(Rp(i)-Rp(i-1)) * (1. /((ru(i)-ru(i-1))))
       b(i)=-(1./(Rp(i+1)-Rp(i))+1./(Rp(i)-Rp(i-1)))* (1. /((ru(i)-ru(i-1))))

!       pdi(i) = 1./(Rp(i)**2)
       a(i)= 1./((Rp(I)-Rp(I-1))*(Ru(I)-Ru(I-1)))
       b(i)=-(1./(Rp(I+1)-Rp(I))+1./(Rp(I)-Rp(I-1)))/
     $      ((Ru(I)-Ru(I-1)))
       c(i)= 1. /((Rp(I+1)-Rp(I))*(Ru(I)-Ru(I-1)))
       end do
       b(1)=   -(1./(Rp(2)-Rp(1)))/(Ru(1)-Ru(0))
       b(imax)=-(1./(Rp(IMAX)-Rp(IMAX-1))) /
     $          ((Ru(IMAX)-Ru(IMAX-1)))
       c(imax)=0.
       a(1)=0.

       do i=1,imax
c       write(*,*) a(i), b(i), c(i)
       enddo

      do k=1,kmax/p_col
	    do j=1,jmax/p_row
	  do i=1,imax
	   p1x(i,j,k)=p(i,j,k)
          enddo
         enddo
      enddo
      call transpose_x_to_y(p1x,p1y)
      do k=1,kmax/p_col
        do i=0,imx
          do j=1,jmax
           ft(j)=p1y(i,j,k)
           enddo
           call vrfftf(1,jmax,ft,dft,1,wj)
           do j=1,jmax
         p1y(i,j,k)=ft(j)
           enddo
         enddo
      enddo
      call transpose_y_to_z(p1y,p1z)
      do i=0,imx
        do j=1,jmax/p_col
          do k=1,kmax
          fz(k)=p1z(i,j,k)
           enddo
           call vrfftf(1,kmax,fz,dfz,1,wk)
          do k=1,kmax
          p1z(i,j,k)=fz(k)
           enddo
         enddo
      enddo
      call transpose_z_to_y(p1z,p1y)
      call transpose_y_to_x(p1y,p1x)
c********************************************************************
       
      ktel  = ( nrank - (nrank/p_col)*p_col) * (kmax/p_col)
      jtel  = (  (nrank/p_col)) * (jmax/p_row)
c      if (mod(1+jtel,2).eq.0) even  = 'n'
c      if (mod(1+jtel,2).eq.1) even  = 'y'
c      if (ini.eq.0) then

       
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1, imax
!         bb(i)=b(i)+(trt(j+jtel)/rp(i)**2+zrt(k+ktel))
          bb(i)=b(i)+(trt(j+jtel)+zrt(k+ktel))
         d(i)= p1x(i,j,k)
c         write(*,*) bb(i), d(i)
        enddo
        call trid_2(imax,a,bb,c,d,dd)
        do i=1, imax
        p1x(i,j,k) = d(i)
        enddo
       enddo
      enddo
c********************************************************************
      call transpose_x_to_y(p1x,p1y)
      do k=1,kmax/p_col
        do i=0,imx
          do j=1,jmax
           ft(j)=p1y(i,j,k)
           enddo
           call vrfftb(1,jmax,ft,dft,1,wj)
           do j=1,jmax
         p1y(i,j,k)=ft(j)
           enddo
         enddo
      enddo
      call transpose_y_to_z(p1y,p1z)
      do i=0,imx
        do j=1,jmax/p_col
          do k=1,kmax
          fz(k)=p1z(i,j,k)
           enddo
           call vrfftb(1,kmax,fz,dfz,1,wk)
          do k=1,kmax
          p1z(i,j,k)=fz(k)
           enddo
         enddo
      enddo
      call transpose_z_to_y(p1z,p1y)
      call transpose_y_to_x(p1y,p1x)
      do k=1,kmax/p_col
        do j=1,jmax/p_row
          do i=1,imax
            p(i,j,k)=p1x(i,j,k)
          enddo
        enddo
      enddo
 
      end



      SUBROUTINE TRID_2 (MR,A,B,C,Y,D)
      DIMENSION       A(1)       ,B(1)       ,C(1)       ,Y(1)       ,
     1                D(1)
      M = MR
      MM1 = M-1
      Z = 1./B(1)
      D(1) = C(1)*Z
      Y(1) = Y(1)*Z
      DO 101 I=2,MM1
         Z = 1./(B(I)-A(I)*D(I-1))
         D(I) = C(I)*Z
         Y(I) = (Y(I)-A(I)*Y(I-1))*Z
  101 CONTINUE
      Z = B(M)-A(M)*D(MM1)
      IF (Z .NE. 0.) GO TO 102
      Y(M) = 0.
      GO TO 103
  102 Y(M) = (Y(M)-A(M)*Y(MM1))/Z
  103 CONTINUE
      DO 104 IP=1,MM1
         I = M-IP
         Y(I) = Y(I)-D(I)*Y(I+1)
  104 CONTINUE
      RETURN
      END
