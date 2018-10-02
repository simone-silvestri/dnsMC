c     ***************      
c     VERSION 2.1
c     ***************

      subroutine momx(dfr,dft,dfz,dfc,u,v,w,c,imax,jmax,kmax,imx,p_row,p_col,ru,rp,dr,dtheta,dz,mr_s,mr_c,rank,px,mu,
     ^                 rho,lambda,e,Twh,Twc)
      use decomp_2d
      integer imax,jmax,kmax,ip,im,rank,px,p_row,p_col,imx
      real dfr(0:imax+1,jmax/p_row,kmax/p_col),dtheta,dz
      real dft(0:imax+1,jmax/p_row,kmax/p_col)
      real dfz(0:imax+1,jmax/p_row,kmax/p_col)
      real dfc(0:imax+1,jmax/p_row,kmax/p_col)
      real u(0:imax+1,jmax/p_row,kmax/p_col)
      real v(0:imax+1,jmax/p_row,kmax/p_col)
      real w(0:imax+1,jmax/p_row,kmax/p_col)
      real c(0:imax+1,jmax/p_row,kmax/p_col)
      real e(0:imax+1,jmax/p_row,kmax/p_col)
      real mu(0:imax+1,jmax/p_row,kmax/p_col)
      real muc(0:imax+1,jmax/p_row,kmax/p_col)
      real rho(0:imax+1,jmax/p_row,kmax/p_col)
      real rhoc(0:imax+1,jmax/p_row,kmax/p_col)
      real lambda(0:imax+1,jmax/p_row,kmax/p_col)
      real rp(0:imax+1),ru(0:imax+1)
      real ft1(jmax),dft1(jmax),ddft1(jmax),visc
      real ft2(jmax),dft2(jmax),ddft2(jmax)
      real ft3(jmax),dft3(jmax),ddft3(jmax),dr
      real ftc(jmax),dfrc(jmax),ddftc(jmax)
      real hi_s0(imax,jmax/p_row,kmax/p_col),hi_c0(0:imax,jmax/p_row,kmax/p_col)
      real hi_s1(imax,jmax/p_row,kmax/p_col),hi_c1(0:imax,jmax/p_row,kmax/p_col)
      real hi_s2(imax,jmax/p_row,kmax/p_col),hi_c2(0:imax,jmax/p_row,kmax/p_col)
      real hic_s2(imax,jmax/p_row,kmax/p_col),hic_c2(0:imax,jmax/p_row,kmax/p_col)
      real hi_s3(imax,jmax/p_row,kmax/p_col),hi_c3(0:imax,jmax/p_row,kmax/p_col)
      real hi_s4(imax,jmax/p_row,kmax/p_col),hi_c4(0:imax,jmax/p_row,kmax/p_col)
      real hi_s5(imax,jmax/p_row,kmax/p_col),hi_c5(0:imax,jmax/p_row,kmax/p_col)
      real hi_s6(imax,jmax/p_row,kmax/p_col),hi_c6(0:imax,jmax/p_row,kmax/p_col)
      real hi_s7(imax,jmax/p_row,kmax/p_col),hi_c7(0:imax,jmax/p_row,kmax/p_col)
      real hi_s8(imax,jmax/p_row,kmax/p_col),hi_c8(0:imax,jmax/p_row,kmax/p_col)
      real hi_s9(imax,jmax/p_row,kmax/p_col),hi_c9(0:imax,jmax/p_row,kmax/p_col)
      real hi_s10(imax,jmax/p_row,kmax/p_col),hi_c10(0:imax,jmax/p_row,kmax/p_col)
      real vdmuridr_s(imax,jmax/p_row,kmax/p_col)
      real hic_s6(imax,jmax/p_row,kmax/p_col),hic_c6(0:imax,jmax/p_row,kmax/p_col)
      real hil_s6(imax,jmax/p_row,kmax/p_col),hil_c6(0:imax,jmax/p_row,kmax/p_col)
      real himu_s6(imax,jmax/p_row,kmax/p_col),himu_c6(0:imax,jmax/p_row,kmax/p_col)
      real himut_s6(imax,jmax/p_row,kmax/p_col),himut_c6(0:imax,jmax/p_row,kmax/p_col)
      real himur_c6(0:imax,jmax/p_row,kmax/p_col)
      real hi_tmp(0:imax,jmax/p_row,kmax/p_col)
      real hir_s(imax,jmax/p_row,kmax/p_col),hir_c(0:imax,jmax/p_row,kmax/p_col)
      real rui(0:imax+1),rpi(0:imax+1),mr_s(imax),mr_c(0:imax)
      real u_s(0:imax+1,jmax/p_row,kmax/p_col)
      real u_st (0:imx,jmax      ,kmax/p_col)
      real u_sf (0:imx,jmax/p_col,kmax)
      real mu_t (0:imx,jmax      ,kmax/p_col)
      real mu_f (0:imx,jmax/p_col,kmax)
      real mdudt_t (0:imx,jmax,kmax/p_col)
      real mdudz_f (0:imx,jmax/p_col,kmax),fz1(kmax),dfz1(kmax)
      real mdudz_t(0:imx,jmax,kmax/p_col)
      real mdudz(0:imax+1,jmax/p_row,kmax/p_col),mdudt(0:imax+1,jmax/p_row,kmax/p_col) 
      real Twh,Twc
        
!      rui = 1./(ru)
!      rpi = 1./rp
      dri =1./dr
       do k=1,kmax/p_col
          do j=1,jmax/p_row
            do i=1,imax
            hi_s1(i,j,k)=v(i,j,k)
            hi_s5(i,j,k)=v(i,j,k)
            hi_s2(i,j,k)=w(i,j,k)
            hic_s2(i,j,k)=e(i,j,k)
            hi_s6(i,j,k)=w(i,j,k)
            hic_s6(i,j,k)=c(i,j,k)
            himu_s6(i,j,k)=mu(i,j,k)
            hir_s(i,j,k)=rho(i,j,k)
            hil_s6(i,j,k)=lambda(i,j,k)
            hi_s7(i,j,k) = mu(i,j,k)
            enddo
            do i=0,imax
            hi_c3(i,j,k)=u(i,j,k)!*ru(i)
            hi_c0(i,j,k)=u(i,j,k)
            enddo
         enddo
       enddo
       call inter_c_s_m (imax,jmax*kmax/px,hi_c0,hi_s0,dr)
       do k=1,kmax/p_col
        do j=1,jmax/p_row
         do i=1,imax
          u_s(i,j,k) = hi_s0(i,j,k)
         enddo
        enddo
       enddo 
       call transpose_x_to_y(u_s,u_st)      
       call transpose_y_to_z(u_st,u_sf)      
       call transpose_x_to_y(mu,mu_t)      
       call transpose_y_to_z(mu_t,mu_f)      

       call inter_s_c_m (imax,jmax*kmax/px,hi_s1,hi_c1,dr) !v on staggered
       call inter_s_c_m (imax,jmax*kmax/px,hi_s2,hi_c2,dr) !w on staggered
       call inter_s_c_m (imax,jmax*kmax/px,hic_s2,hic_c2,dr)
       call der1w_s_c6_m(imax,jmax*kmax/px,hi_s5,hi_c5,dr)   !dv/dx
       call der1w_s_c6_m(imax,jmax*kmax/px,hi_s6,hi_c6,dr)   !dw/dx
       call der1c_s_c6_m(imax,jmax*kmax/px,hic_s6,hic_c6,dr,Twh,Twc) !dT/dx
       call inter_s_c_m (imax,jmax*kmax/px,himu_s6,hi_tmp,dr)
       call deriv_s_c_m (imax,jmax*kmax/px,himu_s6,himu_c6,dr) !dmu/dx
       call deriv_c_s_m (imax,jmax*kmax/px,hi_c0,hi_s0,dr)   !du/dx
       call inter_s_c_m (imax,jmax*kmax/px,hir_s,hir_c,dr)
       call deriv_s_c_m(imax,jmax*kmax/px,hil_s6,hil_c6,dr)  !dlam/dx
       call deriv_s_c_m(imax,jmax*kmax/px,hi_s7,hi_c7,dr)   !dmu/dx repeat

       
      do k=1,kmax/p_col
	     do j=1,jmax/p_row
           do i=0,imax
       hi_c5(i,j,k) = hi_c5(i,j,k)*mr_c(i)  !dv/dr
       hi_c6(i,j,k) = hi_c6(i,j,k)*mr_c(i)  !dw/dr
       hic_c6(i,j,k) = hic_c6(i,j,k)*mr_c(i) !dT/dr
       muc(i,j,k) = hi_tmp(i,j,k)
       hil_c6(i,j,k) = hil_c6(i,j,k)*mr_c(i) !dlam/dr
       himu_c6(i,j,k) = himu_c6(i,j,k)*mr_c(i) !dmu/dr
       hi_c7(i,j,k) = hi_c7(i,j,k)*mr_c(i)      !dmu/dr
           enddo
           do i=1,imax        
       hi_s0(i,j,k) = hi_s0(i,j,k)*mr_s(i)   !du/dr on central 
           enddo
        enddo
       enddo
      call inter_c_s_m(imax,jmax*kmax/px,hi_c7,hi_s7,dr)  !dmu/dr central
      call inter_s_c_m(imax,jmax*kmax/px,hi_s0,hi_c0,dr)  !du/dr staggered

       do k=1,kmax/p_col
        do j=1,jmax/p_row
         do i=0,imax
       hil_c6(i,j,k) = hil_c6(i,j,k)*hic_c6(i,j,k)    !dlam/dr * dT/dr
       himur_c6(i,j,k) = himu_c6(i,j,k)*hi_c0(i,j,k)  !dmu/dr * dU/dr
       himut_c6(i,j,k) = himu_c6(i,j,k)*hi_c5(i,j,k)  !dmu/dr * dv/dr
       himu_c6(i,j,k) = himu_c6(i,j,k)*hi_c6(i,j,k)   !dmu/dr * dw/dr
         enddo
         do i=1,imax
        vdmuridr_s(i,j,k) = hi_s7(i,j,k)*v(i,j,k)*rpi(i) !not needed for channel
         enddo
        enddo
       enddo

       do k=1,kmax/p_col
	  do j=1,jmax/p_row
            do i=0,imax
	      hi_c1(i,j,k)=hi_c1(i,j,k)*hi_c3(i,j,k)*hir_c(i,j,k)  !vu(r)rho r removed see 78
	      hi_c2(i,j,k)=hi_c2(i,j,k)*hi_c3(i,j,k)*hir_c(i,j,k)  !wu(r)rho r removed see 78
	      hic_c2(i,j,k)=hic_c2(i,j,k)*hi_c3(i,j,k)*hir_c(i,j,k)!hu(r)rho r removed see 78
              hi_c3(i,j,k)=hi_c3(i,j,k)*u(i,j,k)*hir_c(i,j,k)  !uu(r)rho r removed see 78
              hi_c5(i,j,k)=hi_c5(i,j,k)!*ru(i)     !dv/dr
              hi_c6(i,j,k)=hi_c6(i,j,k)!*ru(i)     !dw/dr
              hic_c6(i,j,k)=hic_c6(i,j,k)!*ru(i)    !dT/dr
            enddo
          enddo
        enddo
            call inter_c_s_m(imax,jmax*kmax/px,hi_c3,hi_s3,dr)
            call deriv_s_c_m(imax,jmax*kmax/px,hi_s3,hi_c3,dr)  !drho uu/dx
            call deriv_c_s_m(imax,jmax*kmax/px,hi_c1,hi_s1,dr)  !drho uv/dx
            call deriv_c_s_m(imax,jmax*kmax/px,hi_c2,hi_s2,dr)  !drho uw/dx
            call deriv_c_s_m(imax,jmax*kmax/px,hic_c2,hic_s2,dr) !drho uh/dx
            call deriv_c_s_m(imax,jmax*kmax/px,hi_c5,hi_s5,dr)    !d2v/dr2
            call deriv_c_s_m(imax,jmax*kmax/px,hi_c6,hi_s6,dr)    !d2w/dr2
            call deriv_c_s_m(imax,jmax*kmax/px,hic_c6,hic_s6,dr)   !d2T/dr2
            call inter_c_s_m(imax,jmax*kmax/px,himu_c6,himu_s6,dr)
            call inter_c_s_m(imax,jmax*kmax/px,himut_c6,himut_s6,dr)
            call inter_c_s_m(imax,jmax*kmax/px,hil_c6,hil_s6,dr)
 
          do k=1,kmax/p_col
	        do j=1,jmax/p_row
	          do i=0,imax
            hi_c3(i,j,k)=hi_c3(i,j,k)*mr_c(i)           !drho uu/dx
             enddo
             do i=1,imax
            hi_s1(i,j,k)=hi_s1(i,j,k)*mr_s(i)          !drho uv/dx
            hi_s2(i,j,k)=hi_s2(i,j,k)*mr_s(i)          !drho uw/dx
            hic_s2(i,j,k)=hic_s2(i,j,k)*mr_s(i)         !drho uh/dx
            hi_s5(i,j,k)=hi_s5(i,j,k)*mr_s(i)            !d2v/dr2
            hi_s6(i,j,k)=hi_s6(i,j,k)*mr_s(i)            !d2w/dr2
            hic_s6(i,j,k)=hic_s6(i,j,k)*mr_s(i)           !d2T/dr2
             enddo
           enddo
         enddo
        do k=1,kmax/p_col
	  do j=1,jmax/p_row
            do i=1,imax
            im = i -1
            ip = i +1
!            dfr(i,j,k)=dfr(i,j,k)-hi_c3(i,j,k)*Rui(i)+2*muc(i,j,k)*Rui(i)*( rp(ip)*(u(ip,j,k)-u(i ,j,k))/(Ru(ip)-Ru(i)) -
!     ^                                                              rp( i)*(u(i ,j,k)-u(im,j,k))/(Ru(i)-Ru(im)))/(Rp(ip)-Rp(i))
!     ^                                               +2*himur_c6(i,j,k)
!            dft(i,j,k)=dft(i,j,k)-hi_s1(i,j,k)*Rpi(i)+hi_s5(i,j,k)*rpi(i)*mu(i,j,k)
!     ^                                                +himut_s6(i,j,k) - vdmuridr_s(i,j,k)
!            dfz(i,j,k)=dfz(i,j,k)-hi_s2(i,j,k)*Rpi(i)+hi_s6(i,j,k)*rpi(i)*mu(i,j,k)
!     ^                                               +himu_s6(i,j,k)
!            dfc(i,j,k)=dfc(i,j,k)-hic_s2(i,j,k)*Rpi(i)+hic_s6(i,j,k)*rpi(i)*lambda(i,j,k)+hil_s6(i,j,k)
            dfr(i,j,k)=dfr(i,j,k)-hi_c3(i,j,k)+2*muc(i,j,k)*( (u(ip,j,k)-u(i ,j,k))/(Ru(ip)-Ru(i)) -
     ^                                                              (u(i ,j,k)-u(im,j,k))/(Ru(i)-Ru(im)))/(Rp(ip)-Rp(i))    !checkk
     ^                                               +2*himur_c6(i,j,k)
            dft(i,j,k)=dft(i,j,k)-hi_s1(i,j,k)+hi_s5(i,j,k)*mu(i,j,k)
     ^                                                +himut_s6(i,j,k)
            dfz(i,j,k)=dfz(i,j,k)-hi_s2(i,j,k)+hi_s6(i,j,k)*mu(i,j,k)
     ^                                               +himu_s6(i,j,k)
            dfc(i,j,k)=dfc(i,j,k)-hic_s2(i,j,k)+hic_s6(i,j,k)*lambda(i,j,k)+hil_s6(i,j,k)
           enddo
          enddo
         enddo 
   
      do i=0,imx
       do k=1,kmax/p_col  
        do j=1,jmax
         ft1(j)=u_st(i,j,k)
        enddo
        call four1(jmax,ft1,dft1,dtheta)  !du/dy
        do j=1,jmax
         mdudt_t(i,j,k) = mu_t(i,j,k)*dft1(j)  !mudu/dy
        enddo  
       enddo
      enddo
      call transpose_y_to_x(mdudt_t,mdudt)


      do j=1,jmax/p_col
       do i=0,imx
        do k=1,kmax
         fz1(k)=u_sf(i,j,k)
        enddo
        call four1(kmax,fz1,dfz1,dz)  !du/dz
        do k=1,kmax
        mdudz_f(i,j,k) = mu_f(i,j,k)*dfz1(k)  !mudu/dz
        enddo
       enddo 
      enddo
      call transpose_z_to_y(mdudz_f,mdudz_t)
      call transpose_y_to_x(mdudz_t,mdudz)
      
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
         hi_s8(i,j,k) = mdudz(i,j,k)
         hi_s9(i,j,k) = mdudt(i,j,k)!*rpi(i)
         hi_s10(i,j,k) = mdudz(i,j,k)!*rpi(i)
        enddo
c         hi_s9(1,j,k) = 0.0 
       enddo
      enddo
      call deriv_s_c_m(imax,jmax*kmax/px,hi_s8,hi_c8,dr) !d(mudu/dz)/dx
      call deriv_s_c_m(imax,jmax*kmax/px,hi_s9,hi_c9,dr) !d(mudu/dy)/dx
        
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=0,imax
         hi_c8(i,j,k) = hi_c8(i,j,k)*mr_c(i)     !d(mudu/dz)/dr
         hi_c9(i,j,k) = hi_c9(i,j,k)*mr_c(i)     !d(mudu/dy)/dr
        enddo
       enddo
      enddo
      call inter_c_s_m(imax,kmax*jmax/px,hi_c8,hi_s8,dr)
      call inter_c_s_m(imax,kmax*jmax/px,hi_c9,hi_s9,dr)
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax 
          dft(i,j,k) = dft(i,j,k) + hi_s9(i,j,k)
          dfz(i,j,k) = dfz(i,j,k) + hi_s8(i,j,k) !not needed for channel --- + hi_s10(i,j,k)
        enddo
       enddo
      enddo

      end


      subroutine momy(dfr,dft,dfz,dfc,dfr_f,dft_f,dfz_f,dfc_f,u,v,w,c,u_f,v_f,w_f,c_f,
     ^ imax,jmax,kmax,imx,p_row,p_col,ru,rp,dr,dtheta,dz,rank,px,rho,rho_f,mu,mu_f,e,e_f,lambda,lambda_f,mr_s,mr_c)
      use decomp_2d
      integer imax,jmax,kmax,itr,ifil(0:imax+1),rank,px,p_row,p_col,imx,idex,i_dex
      real dfr(0:imax+1,jmax/p_row,kmax/p_col)
      real dft(0:imax+1,jmax/p_row,kmax/p_col)
      real dfz(0:imax+1,jmax/p_row,kmax/p_col)
      real dfc(0:imax+1,jmax/p_row,kmax/p_col)
      real u  (0:imax+1,jmax/p_row,kmax/p_col)
      real v  (0:imax+1,jmax/p_row,kmax/p_col)
      real w  (0:imax+1,jmax/p_row,kmax/p_col)
      real c  (0:imax+1,jmax/p_row,kmax/p_col)
      real rho  (0:imax+1,jmax/p_row,kmax/p_col)
      real mu  (0:imax+1,jmax/p_row,kmax/p_col)
      real muc  (0:imax+1,jmax/p_row,kmax/p_col)
      real e  (0:imax+1,jmax/p_row,kmax/p_col)
      real e_f (0:imx,jmax      ,kmax/p_col)
      real lambda  (0:imax+1,jmax/p_row,kmax/p_col)
      real lambda_f (0:imx,jmax      ,kmax/p_col)
      
      real v_c(0:imax+1,jmax/p_row,kmax/p_col)
      real rhoc(0:imax+1,jmax/p_row,kmax/p_col)


      real u_s(0:imax+1,jmax/p_row,kmax/p_col)

      real u_f (0:imx,jmax      ,kmax/p_col)
      real v_f (0:imx,jmax      ,kmax/p_col)
      real v_ff (0:imx,jmax/p_col,kmax)
      real w_f (0:imx,jmax      ,kmax/p_col)
      real c_f (0:imx,jmax      ,kmax/p_col)
      real rho_f (0:imx,jmax      ,kmax/p_col)
      real v_cf(0:imx,jmax      ,kmax/p_col)
      real rhoc_f(0:imx,jmax      ,kmax/p_col)
      real u_sf(0:imx,jmax      ,kmax/p_col)
      real mu_f (0:imx,jmax      ,kmax/p_col)
      real mu_ff (0:imx,jmax/p_col,kmax)
      real muc_f(0:imx,jmax      ,kmax/p_col)

      real dfr_f(0:imx,jmax,kmax/p_col)
      real dft_f(0:imx,jmax,kmax/p_col)
      real dfz_f(0:imx,jmax,kmax/p_col)
      real dfc_f(0:imx,jmax,kmax/p_col)
  
      real rp(0:imax+1),ru(0:imax+1)

      real ft1(jmax),dft1(jmax),ddft1(jmax),visc
      real ft2(jmax),dft2(jmax),ddft2(jmax)
      real ft3(jmax),dft3(jmax),ddft3(jmax),dtheta,dz
      real ftc(jmax),dftc(jmax),ddftc(jmax)
      real fte(jmax),dfte(jmax)
      real ftl(jmax),dftl(jmax)
      real hi_s(imax,jmax/p_row,kmax/p_col),hi_c(0:imax,jmax/p_row,kmax/p_col)
      real hi2_s(imax,jmax/p_row,kmax/p_col),hi2_c(0:imax,jmax/p_row,kmax/p_col)
      real hi3_s(imax,jmax/p_row,kmax/p_col),hi3_c(0:imax,jmax/p_row,kmax/p_col)
      real hir_s(imax,jmax/p_row,kmax/p_col),hir_c(0:imax,jmax/p_row,kmax/p_col)
      real him_s(imax,jmax/p_row,kmax/p_col),him_c(0:imax,jmax/p_row,kmax/p_col)

      real rui(0:imax+1),rpi(0:imax+1)
      real mudvdr_c(0:imax+1,jmax/p_row,kmax/p_col)      
      real mudvdr_cf(0:imx,jmax,kmax/p_col)      
      real mr_s(imax),mr_c(0:imax)

      real ft1r(jmax),dft1r(jmax)
      real ft2r(jmax),dft2r(jmax)
      real ft3r(jmax),dft3r(jmax)
      real ftm(kmax),dftm(kmax)
      real ftmc(kmax),dftmc(kmax)
      real mridwdth_t(0:imx,jmax,kmax/p_col)
      real mridwdth_f(0:imx,jmax/p_col,kmax)
      real dmridwdth_t(0:imx,jmax,kmax/p_col)
      real dmridwdth_f(0:imx,jmax/p_col,kmax)
      real mudvdz_t(0:imx,jmax,kmax/p_col)
      real mudvdz_f(0:imx,jmax/p_col,kmax)

      real fz1(kmax),dfz1(kmax)
      real fz2(kmax),dfz2(kmax)


!      idex = i_dex(nrank)
      
!      RUI = 1./(RU)
!      RPI = 1./RP
      do k=1,kmax/p_col
	  do j=1,jmax/p_row
	  do i=1,imax
	    hi_s(i,j,k)=v(i,j,k)
	    hi3_s(i,j,k)=v(i,j,k)!*rpi(i)
	    hir_s(i,j,k)=rho(i,j,k)
	    him_s(i,j,k)=mu(i,j,k)
          enddo
        enddo
       enddo
 
           call der1w_s_c6_m(imax,jmax*kmax/px,hi3_s,hi2_c,dr)
!          call deriv_s_c_m(imax,jmax*kmax/(p_col*p_row),hi3_s,hi2_c,dr)  !dv/dx
       do k=1,kmax/p_col
        do j=1,jmax/p_row 
         do i=0,imax
          hi2_c(i,j,k) = hi2_c(i,j,k)*mr_c(i)  !dv/dr
         enddo 
        enddo
       enddo
          call inter_s_c_m(imax,jmax*kmax/(p_col*p_row),hi_s,hi_c,dr)    !v_c
          call inter_s_c_m(imax,jmax*kmax/(p_col*p_row),hir_s,hir_c,dr)  !rho_c
          call inter_s_c_m(imax,jmax*kmax/(p_col*p_row),him_s,him_c,dr)  !mu_c
      do k=1,kmax/p_col
    	do j=1,jmax/p_row
         do i=0,imax
	  v_c(i,j,k)=hi_c(i,j,k)
	  rhoc(i,j,k)=hir_c(i,j,k)
	  muc(i,j,k)=him_c(i,j,k)
          mudvdr_c(i,j,k) = muc(i,j,k)*hi2_c(i,j,k)!*ru(i)   !mudv/dr
         enddo
        enddo
      enddo
        
      do k=1,kmax/p_col
    	do j=1,jmax/p_row
	      do i=0,imax
	    hi_c(i,j,k)=u(i,j,k)
          enddo
        enddo
       enddo
      call inter_c_s_m(imax,jmax*kmax/px,hi_c,hi_s,dr)
       do k=1,kmax/p_col
    	do j=1,jmax/p_row
          do i=1,imax
	  u_s(i,j,k)=hi_s(i,j,k)
          enddo
        enddo
      enddo
      call transpose_x_to_y(mudvdr_c,mudvdr_cf)
      call transpose_x_to_y(v_c,v_cf)
      call transpose_x_to_y(u_s,u_sf)
      call transpose_x_to_y(rhoc,rhoc_f)
      call transpose_x_to_y(muc,muc_f)
      call transpose_y_to_z(v_f,v_ff)
      call transpose_y_to_z(mu_f,mu_ff)

      do i=1,imax
	  ifil(i)=16+12*(i-1)!7+8*(i-1)
      enddo
      do i=1,imax
    	 if (ifil(i).gt. jmax) ifil(i)=jmax
      enddo
      ifil(0)=ifil(1)
      ifil(imax+1)=ifil(imax)
      do i=0,imx
        do k=1,kmax/p_col
          do j=1,jmax
            ft1(j)=U_f(i,j,k)!*(rui(i+idex)**2)
            ft2(j)=V_f(i,j,k)!*(rpi(i+idex)**2)
            ft3(j)=W_f(i,j,k)!*(rpi(i+idex)**2)
            ftc(j)=c_f(i,j,k)!*(rpi(i+idex)**2)

          enddo
           itr = jmax!ifil (i+idex)
          call four2r(jmax,itr,ft1,ddft1,dtheta)  !d2u/dy2
          call four2r(jmax,itr,ft2,ddft2,dtheta)  !d2v/dy2
          call four2r(jmax,itr,ft3,ddft3,dtheta)  !d2w/dy2
          call four2r(jmax,itr,ftc,ddftc,dtheta)  !d2T/dy2
          do j=1,jmax
	     dfr_f(i,j,k)=dfr_f(i,j,k)+muc_f(i,j,k)*ddft1(j)
	     dft_f(i,j,k)=dft_f(i,j,k)+2*mu_f(i,j,k)*ddft2(j)
	     dfz_f(i,j,k)=dfz_f(i,j,k)+mu_f(i,j,k)*ddft3(j)
             dfc_f(i,j,k)=dfc_f(i,j,k)+lambda_f(i,j,k)*ddftc(j)
          enddo
        enddo
       enddo
       do i=0,imx
	    do k=1,kmax/p_col
	     do j=1,jmax
	       ft1(j)=    rhoc_f(i,j,k)*U_f(i,j,k)*v_cf (i,j,k)  !rho uv
           ft2(j)=    rho_f(i,j,k)*v_f(i,j,k)*v_f  (i,j,k)   !rho vv
           ft3(j)=    rho_f(i,j,k)*w_f(i,j,k)*v_f  (i,j,k)   !rho wv
           fte(j)=    rho_f(i,j,k)*e_f(i,j,k)*v_f  (i,j,k)   !rho hv
         enddo
          call four1(jmax,ft1,dft1,dtheta)   !d(rho uv)/dy
          call four1(jmax,ft2,dft2,dtheta)   !d(rho vv)/dy
          call four1(jmax,ft3,dft3,dtheta)   !d(rho wv)/dy
          call four1(jmax,fte,dfte,dtheta)   !d(rho hv)/dy
          do j=1,jmax
	     dfr_f(i,j,k)=dfr_f(i,j,k)-0.5*dft1(j)!*rui(i+idex) !not needed ---- +           rhoc_f(i,j,k)*v_cf(i,j,k)**2*rui(i+idex)
	     dft_f(i,j,k)=dft_f(i,j,k)-0.5*dft2(j)!*rpi(i+idex) !not needed ---- -rho_f(i,j,k)*v_f(i,j,k)*u_sf(i,j,k)   *rpi(i+idex)
	     dfz_f(i,j,k)=dfz_f(i,j,k)-0.5*dft3(j)!*rpi(i+idex)
         dfc_f(i,j,k)=dfc_f(i,j,k)-0.5*dfte(j)!*rpi(i+idex)
          enddo
          enddo
          enddo
          do i=0,imx
            do k=1,kmax/p_col
	            do j=1,jmax
	      ft1(j) = u_f(i,j,k)
	      ft2(j) = v_f(i,j,k)
	      ft3(j) = w_f(i,j,k)
          ftc(j) = c_f(i,j,k)
          fte(j) = e_f(i,j,k)
	      ft1r(j) = rhoc_f(i,j,k)*u_f(i,j,k)
	      ft2r(j) = rho_f(i,j,k)*v_f(i,j,k)
	      ft3r(j) = rho_f(i,j,k)*w_f(i,j,k)
              ftm(j) = mu_f(i,j,k)
              ftmc(j) = muc_f(i,j,k)
              ftl(j)=lambda_f(i,j,k)
         enddo
          call four1(jmax,ft1,dft1,dtheta) !du/dy
          call four1(jmax,ft2,dft2,dtheta) !dv/dy
          call four1(jmax,ft3,dft3,dtheta) !dw/dy
          call four1(jmax,ftc,dftc,dtheta) !dT/dy
          call four1(jmax,ft1r,dft1r,dtheta) !drho u/dy
          call four1(jmax,ft2r,dft2r,dtheta) !drho v/dy
          call four1(jmax,ft3r,dft3r,dtheta) !drho w/dy
          call four1(jmax,ftm,dftm,dtheta)   !dmu/dy
          call four1(jmax,ftmc,dftmc,dtheta) !dmu/dy !staggerd
          call four1(jmax,fte,dfte,dtheta)   !dh/dy
          call four1(jmax,ftl,dftl,dtheta)   !dlam/dy
            do j=1,jmax
!	      dfr_f(i,j,k)=dfr_f(i,j,k)-0.5*    (v_cf(i,j,k)    )*dft1r(j)*rui(i+idex)   !vdrho u/dy
!     ^                                 +dftm(j)*dft1(j)*(rpi(i+idex))**2  !dmu/dy du/dy
!	      dft_f(i,j,k)=dft_f(i,j,k)-0.5*    (v_f(i,j,k)     )*dft2r(j)*rpi(i+idex)-0.5*rho_f(i,j,k)*v_f(i,j,k)*rpi(i+idex)*dft2(j)
!     ^                                 +2*dftm(j)*dft2(j)*(rpi(i+idex))**2
!	      dfz_f(i,j,k)=dfz_f(i,j,k)-0.5*    (v_f(i,j,k)     )*dft3r(j)*rpi(i+idex)-0.5*rho_f(i,j,k)*w_f(i,j,k)*rpi(i+idex)*dft2(j)
!     ^                                 +dftm(j)*dft3(j)*(rpi(i+idex))**2
!              dfc_f(i,j,k)=dfc_f(i,j,k)-0.5*(rho_f(i,j,k)*v_f(i,j,k))*dfte(j)*rpi(i+idex)-0.5*e_f(i,j,k)*rpi(i+idex)*dft2r(j)
!     ^                                 +dftl(j)*dftc(j)*(rpi(i+idex))**2
          dfr_f(i,j,k)=dfr_f(i,j,k)-0.5*    (v_cf(i,j,k)    )*dft1r(j)   !vdrho u/dy
     ^                                 +dftm(j)*dft1(j)  !dmu/dy du/dy
          dft_f(i,j,k)=dft_f(i,j,k)-0.5*    (v_f(i,j,k)     )*dft2r(j)-0.5*rho_f(i,j,k)*v_f(i,j,k)*dft2(j)  !-0.5vdrho v/dy-0.5rho vdv/dy
     ^                                 +2*dftm(j)*dft2(j)  !2dmu/dydv/dy
          dfz_f(i,j,k)=dfz_f(i,j,k)-0.5*    (v_f(i,j,k)     )*dft3r(j)-0.5*rho_f(i,j,k)*w_f(i,j,k)*dft2(j)  !-0.5vdrho w/dy-0.5rho wdv/dy
     ^                                 +dftm(j)*dft3(j)                  !dmu/dydw/dy
          dfc_f(i,j,k)=dfc_f(i,j,k)-0.5*(rho_f(i,j,k)*v_f(i,j,k))*dfte(j)-0.5*e_f(i,j,k)*dft2r(j)   !-0.5rho vdh/dy-0.5hdrho v/dy
     ^                                 +dftl(j)*dftc(j)   !dlam/dydT/dy

            enddo
           enddo
          enddo
          do i=0,imx
            do k=1,kmax/p_col
	    do j=1,jmax
	      ft1(j) = v_cf (i,j,k)
	      ft2(j) = muc_f(i,j,k)*v_cf (i,j,k)  !muv
         enddo
          call four1(jmax,ft1,dft1,dtheta)
          call four1(jmax,ft2,dft2,dtheta)
            do j=1,jmax
!	      dfr_f(i,j,k)=dfr_f(i,j,k)-0.5*rhoc_f(i,j,k)*u_f(i,j,k)*Rui(i+idex)*dft1(j)                   !Skew completion
!     ^                                 -(2*u_f(i,j,k)*Rui(i+idex)**2+2*Rui(i+idex)**2*dft1(j))*muc_f(i,j,k)  ! White pagina 793
!     ^                                 -rui(i+idex)**2*dft2(j)
          dfr_f(i,j,k)=dfr_f(i,j,k)-0.5*rhoc_f(i,j,k)*u_f(i,j,k)*dft1(j)       !Skew completion -0.5rho u dv/dy checkkkkkkkkk
            enddo
           enddo
         enddo
               

         do k=1,kmax  /  p_col  
            do i=0,imx
              do j=1,jmax
                 ft1(j)=u_sf(i,j,k)
                 ft2(j)=mu_f(i,j,k)*u_sf(i,j,k)
                 ft3(j)=w_f(i,j,k)
              enddo
              call four1(jmax,ft1,dft1,dtheta)
              call four1(jmax,ft2,dft2,dtheta)
              call four1(jmax,ft3,dft3,dtheta)
              do j=1,jmax
!              mridwdth_t(i,j,k) = dft3(j)*mu_f(i,j,k)*rpi(i+idex)
!                dft_f(i,j,k)=dft_f(i,j,k)+(-1*v_f(i,j,k)*Rpi(i+idex)**2 + 2*Rpi(i+idex)**2 *  dft1(j))*mu_f(i,j,k)  ! White pagina 793
!     ^                                   +2*rpi(i+idex)**2*dft2(j)
              mridwdth_t(i,j,k) = dft3(j)*mu_f(i,j,k)  !mudw/dy
                dft_f(i,j,k)=dft_f(i,j,k)
             enddo
           enddo
         enddo

        call transpose_y_to_z(mridwdth_t,mridwdth_f)
        do j=1,jmax/p_col
         do i=0,imx
          do k=1,kmax
           fz1(k) = mridwdth_f(i,j,k)
           fz2(k) = v_ff(i,j,k)
          enddo
          call four1(kmax,fz1,dfz1,dz)    !d(mudw/dy)/dz
          call four1(kmax,fz2,dfz2,dz)    !dv/dz
          do k=1,kmax
          dmridwdth_f(i,j,k) = dfz1(k)
          mudvdz_f(i,j,k) = mu_ff(i,j,k)*dfz2(j)  !mudv/dz
          enddo
         enddo 
        enddo
        call transpose_z_to_y(dmridwdth_f,dmridwdth_t) !d(mudw/dy)/dz
        call transpose_z_to_y(mudvdz_f,mudvdz_t)
 
        do k=1,kmax/p_col
         do i=0,imx
          do j=1,jmax
           ft1(j)=mudvdr_cf(i,j,k)
           ft2(j)=mudvdz_t(i,j,k)
          enddo
          call four1(jmax,ft1,dft1,dtheta)  !d(mudv/dr)/dy
          call four1(jmax,ft2,dft2,dtheta)  !d(mudv/dz)/dy
          do j=1,jmax
!           dfr_f(i,j,k) = dfr_f(i,j,k)+rui(i+idex)*dft1(j)               !From transpose of shear
!           dft_f(i,j,k) = dft_f(i,j,k)+dmridwdth_t(i,j,k)
!           dfz_f(i,j,k) = dfz_f(i,j,k)+rpi(i+idex)*dft2(j)
           dfr_f(i,j,k) = dfr_f(i,j,k)+dft1(j)               !From transpose of shear
           dft_f(i,j,k) = dft_f(i,j,k)+dmridwdth_t(i,j,k)
           dfz_f(i,j,k) = dfz_f(i,j,k)+dft2(j)
          enddo
         enddo
        enddo

 
        call transpose_y_to_x(dfr_f,dfr)
        call transpose_y_to_x(dft_f,dft)
        call transpose_y_to_x(dfz_f,dfz)
        call transpose_y_to_x(dfc_f,dfc)
      end


      subroutine momz(dfr,dft,dfz,dfc,drr,dtt,dzz,dzc,u,v,w,c,u_t,v_t,w_t,c_t,imax,jmax,kmax,imx,p_row,
     ^                                           p_col,ru,rp,dr,dtheta,dz,rank,px,rho,rho_t,mu,mu_t,e,e_t,lambda,lambda_t,mr_s,mr_c)
      use decomp_2d
      integer imax,jmax,kmax,rank,px,p_col,p_row
      real dfr(0:imax+1,jmax/p_row,kmax/p_col),dtheta,dz
      real dft(0:imax+1,jmax/p_row,kmax/p_col)
      real dfz(0:imax+1,jmax/p_row,kmax/p_col)
      real dfc(0:imax+1,jmax/p_row,kmax/p_col)
      real u_s  (imax+1,jmax/p_row,kmax/p_col)
      real u  (0:imax+1,jmax/p_row,kmax/p_col)
      real v  (0:imax+1,jmax/p_row,kmax/p_col)
      real w  (0:imax+1,jmax/p_row,kmax/p_col)
      real c  (0:imax+1,jmax/p_row,kmax/p_col)                  !Added scalar transport variable c
      real e  (0:imax+1,jmax/p_row,kmax/p_col)                  !Added scalar transport variable c
      real rho (0:imax+1,jmax/p_row,kmax/p_col)                 !Added scalar transport variable c
      real mu  (0:imax+1,jmax/p_row,kmax/p_col)                 !Added scalar transport variable c
      real lambda  (0:imax+1,jmax/p_row,kmax/p_col)             !Added scalar transport variable c
      real wc (0:imax+1,jmax/p_row,kmax/p_col)
      real rhoc (0:imax+1,jmax/p_row,kmax/p_col)                !Added scalar transport variable c
      real muc (0:imax+1,jmax/p_row,kmax/p_col)                 !Added scalar transport variable c
      real u_f(0:imx,jmax/p_col,kmax)
      real u_sf(0:imx,jmax/p_col,kmax)
      real v_f(0:imx,jmax/p_col,kmax)
      real w_f(0:imx,jmax/p_col,kmax)
      real c_f(0:imx,jmax/p_col,kmax)
      real e_f(0:imx,jmax/p_col,kmax)
      real rho_f(0:imx,jmax/p_col,kmax)
      real mu_f(0:imx,jmax/p_col,kmax)
      real lambda_f(0:imx,jmax/p_col,kmax)

      real dfr_f(0:imx,jmax/p_col,kmax)
      real dft_f(0:imx,jmax/p_col,kmax)
      real dfz_f(0:imx,jmax/p_col,kmax)
      real dfc_f(0:imx,jmax/p_col,kmax)

      real w_cc(0:imx,jmax/p_col,kmax)
      real rhoc_f(0:imx,jmax/p_col,kmax)
      real muc_f(0:imx,jmax/p_col,kmax)

      real drr(0:imx,jmax,kmax/p_col)
      real dtt(0:imx,jmax,kmax/p_col)
      real dzz(0:imx,jmax,kmax/p_col)
      real dzc(0:imx,jmax,kmax/p_col)

      real u_t(0:imx,jmax,kmax/p_col)
      real u_st(0:imx,jmax,kmax/p_col)
      real v_t(0:imx,jmax,kmax/p_col)
      real w_t(0:imx,jmax,kmax/p_col)
      real c_t(0:imx,jmax,kmax/p_col)
      real e_t(0:imx,jmax,kmax/p_col)
      real rho_t(0:imx,jmax,kmax/p_col)
      real mu_t(0:imx,jmax,kmax/p_col)
      real lambda_t(0:imx,jmax,kmax/p_col)
 
      real tmp (0:imx,jmax,kmax/p_col)
      real hi_s(imax,jmax/p_row,kmax/p_col),hi_c(0:imax,jmax/p_row,kmax/p_col)
      real hi_s3(imax,jmax/p_row,kmax/p_col),hi_c3(0:imax,jmax/p_row,kmax/p_col)
      real hi2_s(imax,jmax/p_row,kmax/p_col),hi2_c(0:imax,jmax/p_row,kmax/p_col)   !Note: the dimension of hi* important for deriv and inter
      real mudwdr_c(0:imax+1,jmax/p_row,kmax/p_col)
      real mudwdr_cf(0:imx,jmax/p_col,kmax)
      real hir_s(imax,jmax/p_row,kmax/p_col),hir_c(0:imax,jmax/p_row,kmax/p_col)
      real him_s(imax,jmax/p_row,kmax/p_col),him_c(0:imax,jmax/p_row,kmax/p_col)
      real mududz_f(0:imx,jmax/p_col,kmax) 
      real mududz_t(0:imx,jmax,kmax/p_col)

      real rp(0:imax+1),ru(0:imax+1)

      real fz1(kmax),dfz1(kmax),ddfz1(kmax)
      real fz2(kmax),dfz2(kmax),ddfz2(kmax)
      real fz3(kmax),dfz3(kmax),ddfz3(kmax)
      real fzc(kmax),dfzc(kmax),ddfzc(kmax)
      real fze(kmax),dfze(kmax)
 
      real fz1r(kmax),dfz1r(kmax)
      real fz2r(kmax),dfz2r(kmax)
      real fz3r(kmax),dfz3r(kmax)
      real fzm(kmax),dfzm(kmax)
      real fzmc(kmax),dfzmc(kmax)
      real fzl(kmax),dfzl(kmax)
      real mr_s(imax),mr_c(0:imax)
      real rui(0:imax+1),rpi(0:imax+1)

      idex = i_dex(nrank)
      
!      RUI = 1./(RU)
!      RPI = 1./RP
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=0,imax
         hi_c3(i,j,k) = u(i,j,k)
        enddo
       enddo
      enddo
      call inter_c_s_m(imax,jmax*kmax/px,hi_c3,hi_s3,dr)
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
         u_s(i,j,k) = hi_s3(i,j,k)
        enddo
       enddo 
      enddo
      
      call transpose_x_to_y(u,u_t)
      call transpose_y_to_z(u_t,u_f)
      call transpose_x_to_y(v,v_t)
      call transpose_y_to_z(v_t,v_f)

      call transpose_x_to_y(c,c_t)
      call transpose_y_to_z(c_t,c_f)

      call transpose_x_to_y(w,w_t)
      call transpose_y_to_z(w_t,w_f)
      

      call transpose_x_to_y(rho,rho_t)
      call transpose_y_to_z(rho_t,rho_f)

      call transpose_x_to_y(mu,mu_t)
      call transpose_y_to_z(mu_t,mu_f)

      call transpose_x_to_y(e,e_t)
      call transpose_y_to_z(e_t,e_f)
      call transpose_x_to_y(lambda,lambda_t)
      call transpose_y_to_z(lambda_t,lambda_f)

      call transpose_x_to_y(u_s,u_st)
      call transpose_y_to_z(u_st,u_sf)



      
      do k=1,kmax/p_col
	  do j=1,jmax/p_row
	  do i=1,imax
	    hi_s(i,j,k)=w(i,j,k)
	    hir_s(i,j,k)=rho(i,j,k)
	    him_s(i,j,k)=mu(i,j,k)
          enddo
        enddo
       enddo
       call der1w_s_c6_m(imax,jmax*kmax/px,hi_s,hi2_c,dr)
!      call deriv_s_c_m(imax,jmax*kmax/px,hi_s,hi2_c,dr)   !dw/dx
       do k=1,kmax/p_col
        do j=1,jmax/p_row
         do i=0,imax
          hi2_c(i,j,k) = hi2_c(i,j,k)*mr_c(i)       !dw/dr
         enddo
        enddo
       enddo
          call inter_s_c_m(imax,jmax*kmax/px,hi_s,hi_c,dr)
          call inter_s_c_m(imax,jmax*kmax/px,hir_s,hir_c,dr)
          call inter_s_c_m(imax,jmax*kmax/px,him_s,him_c,dr)
       do k=1,kmax/p_col
         do j=1,jmax/p_row
          do i=0,imax
	  wc(i,j,k)=hi_c(i,j,k)
	  rhoc(i,j,k)=hir_c(i,j,k)
	  muc(i,j,k)=him_c(i,j,k)
          mudwdr_c(i,j,k) = muc(i,j,k)*hi2_c(i,j,k)    !mudw/dr
          enddo
        enddo
      enddo

      call transpose_x_to_y(mudwdr_c,tmp)
      call transpose_y_to_z(tmp,mudwdr_cf)

      call transpose_x_to_y(rhoc,tmp)
      call transpose_y_to_z(tmp,rhoc_f)
      call transpose_x_to_y(wc,tmp)
      call transpose_y_to_z(tmp,w_cc)
      call transpose_x_to_y(muc,tmp)
      call transpose_y_to_z(tmp,muc_f)

       do j=1,jmax/p_col
        do i=0,imx
	  do k=1,kmax
	   fz1(k)=    rhoc_f(i,j,k)*u_f(i,j,k)*w_cc(i,j,k)   !rho u w
           fz2(k)=    rho_f(i,j,k)*v_f(i,j,k)*w_f(i,j,k) !rho vw
           fz3(k)=    rho_f(i,j,k)*w_f(i,j,k)**2.         !rho ww
           fze(k)=    rho_f(i,j,k)*w_f(i,j,k)*e_f(i,j,k) !rho wh                        !Added scalar transport variable c
         enddo
          call four1(kmax,fz1,dfz1,dz)    !drho u w/dz
          call four1(kmax,fz2,dfz2,dz)    !drho vw /dz
          call four1(kmax,fz3,dfz3,dz)    !drho ww /dz
          call four1(kmax,fze,dfze,dz)    !drho wh /dz                       !Added
          do k=1,kmax
	     dfr_f(i,j,k)=-0.5*dfz1(k)
	     dft_f(i,j,k)=-0.5*dfz2(k)
	     dfz_f(i,j,k)=-0.5*dfz3(k)
             dfc_f(i,j,k)=-0.5*dfze(k)                            !Added
          enddo
          enddo
          enddo
       do i=0,imx
    	do j=1,jmax/p_col
	     do k=1,kmax

	   fz1(k)=u_f(i,j,k)
           fz2(k)=v_f(i,j,k)
           fz3(k)=w_f(i,j,k)
           fzc(k)=c_f(i,j,k)                                      !Added
           fze(k)=e_f(i,j,k)                                      !Added
           fzl(k)=lambda_f(i,j,k)                                      !Added
           fzm(k)=mu_f(i,j,k)                                      !Added
           fzmc(k)=muc_f(i,j,k)                                      !Added
	   fz1r(k)=rhoc_f(i,j,k)*u_f(i,j,k)
           fz2r(k)=rho_f(i,j,k)*v_f(i,j,k)
           fz3r(k)=rho_f(i,j,k)*w_f(i,j,k)
         enddo
          call four12(kmax,fz1,dfz1,ddfz1,dz)   !du/dz d2u/dz2
          call four12(kmax,fz2,dfz2,ddfz2,dz)   !dv/dz d2v/dz2
          call four12(kmax,fz3,dfz3,ddfz3,dz)   !dw/dz d2w/dz2
          call four12(kmax,fzc,dfzc,ddfzc,dz)   !dT/dz d2T/dz2            !Added
          call four1(kmax,fz1r,dfz1r,dz)        !drho u/dz
          call four1(kmax,fz2r,dfz2r,dz)        !drho v/dz
          call four1(kmax,fz3r,dfz3r,dz)        !drho w/dz
          call four1(kmax,fzm,dfzm,dz)          !dmu/dz
          call four1(kmax,fzmc,dfzmc,dz)        !dmu/dz
          call four1(kmax,fzl,dfzl,dz)          !dlam/dz
          call four1(kmax,fze,dfze,dz)          !dh/dz
         do k=1,kmax
!	   dfr_f(i,j,k)=dfr_f(i,j,k)+muc_f(i,j,k)*ddfz1(k)-0.5*    (w_cc(i,j,k)     )*dfz1r(k)
!     ^                              +dfzmc(k)*dfz1(k)
!	   dft_f(i,j,k)=dft_f(i,j,k)+mu_f(i,j,k)*ddfz2(k)-0.5*    (w_f (i,j,k)     )*dfz2r(k)-0.5* rho_f(i,j,k)*v_f(i,j,k)*dfz3(k)
!     ^                              +dfzm(k)*dfz2(k)
!	   dfz_f(i,j,k)=dfz_f(i,j,k)+2*mu_f(i,j,k)*ddfz3(k)-0.5*    (w_f (i,j,k)     )*dfz3r(k)-0.5* rho_f(i,j,k)*w_f(i,j,k)*dfz3(k)
!     ^                              +2*dfzm(k)*dfz3(k)
!           dfc_f(i,j,k)=dfc_f(i,j,k)+lambda_f(i,j,k)*ddfzc(k)
!     ^                              -0.5*(rho_f(i,j,k)*w_f (i,j,k))*dfze(k)-0.5*(e_f(i,j,k))*dfz3r(k)          !Added
!     ^                              +dfzl(k)*dfzc(k)
       dfr_f(i,j,k)=dfr_f(i,j,k)+muc_f(i,j,k)*ddfz1(k)-0.5*    (w_cc(i,j,k)     )*dfz1r(k)
     ^                              +dfzmc(k)*dfz1(k)
       dft_f(i,j,k)=dft_f(i,j,k)+mu_f(i,j,k)*ddfz2(k)-0.5*    (w_f (i,j,k)     )*dfz2r(k)-0.5* rho_f(i,j,k)*v_f(i,j,k)*dfz3(k)
     ^                              +dfzm(k)*dfz2(k)
       dfz_f(i,j,k)=dfz_f(i,j,k)+2*mu_f(i,j,k)*ddfz3(k)-0.5*    (w_f (i,j,k)     )*dfz3r(k)-0.5* rho_f(i,j,k)*w_f(i,j,k)*dfz3(k)
     ^                              +2*dfzm(k)*dfz3(k)
       dfc_f(i,j,k)=dfc_f(i,j,k)+lambda_f(i,j,k)*ddfzc(k)
     ^                              -0.5*(rho_f(i,j,k)*w_f (i,j,k))*dfze(k)-0.5*(e_f(i,j,k))*dfz3r(k)          !Added
     ^                              +dfzl(k)*dfzc(k)
         enddo
         enddo
         enddo
       do i=0,imx
	    do j=1,jmax/p_col
   	     do k=1,kmax
	   fz1(k)=     w_cc(i,j,k)
	   fz2(k)=     mudwdr_cf(i,j,k)
         enddo
          call four1(kmax,fz1,dfz1,dz)  !dw/dz
          call four1(kmax,fz2,dfz2,dz)  !d(mudwdr)/dz
         do k=1,kmax
	  dfr_f(i,j,k)=dfr_f(i,j,k)-0.5* rhoc_f(i,j,k)*u_f(i,j,k)*dfz1(k)
         dfr_f(i,j,k)=dfr_f(i,j,k) + dfz2(k)
         enddo
        enddo
       enddo

         call transpose_z_to_y(dfr_f,drr)
         call transpose_z_to_y(dft_f,dtt)
         call transpose_z_to_y(dfz_f,dzz)
         call transpose_z_to_y(dfc_f,dzc)

       end
