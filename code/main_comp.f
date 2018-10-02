      use decomp_2d
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'

      integer ini,iload,ploc,nstap,ierr,istap,ini_sol,icount,i_dex,j_dex
      real ut(jmax),dut(jmax),ener,t_som,bulk_tot, stress_tot,f_colebrook,f_shokling
      real Wbulk, ReB,bulk,stress,instress,time,dpdx,Lz,ti_s(imax),ti_c(0:imax),vm(imax)
      real Cbulk,fluxPipe,fluxAnnulus,Qval,cpipe,dp_ref,meanVisc,source,cannulus
      real hir_s(  imax,jmax/p_row,kmax/p_col)
      real hir_c(0:imax,jmax/p_row,kmax/p_col)

      real int4,error2,int5,errormax,int6,int7
      real start,finish,start1,finish1,opthck
      real num,deno,temp,Voldom,small,toll,Ncr,coeff,coeff2
      real(kind=4) :: qr (0:imax+1,0:jmax/p_row+1,0:kmax+1)
      real(kind=4) :: vr (0:imax+1,0:jmax/p_row+1,0:kmax+1)
      real(kind=4) :: TMC(0:imax+1,0:jmax+1,0:kmax+1)
      real kP(nTemp),Tnb(nTemp)


      integer result_proc, color, MPI_COMM_NODE
      integer node_rank, node_size

      logical periodic_bc(3),period(2)
      character*(MPI_MAX_PROCESSOR_NAME) node_name

      call mpi_init(ierr)

      periodic_bc(1)=.false.
      periodic_bc(2)=.true.
      periodic_bc(3)=.true.

      call decomp_2d_init(imax+2,jmax,kmax,p_row,p_col,periodic_bc)
      rank = nrank
      if (nrank.eq.0) then
      write(6,*) xsize(1),xsize(2),xsize(3),xsize(1)*xsize(2)*xsize(3)
      write(6,*) ysize(1),ysize(2),ysize(3),ysize(1)*ysize(2)*ysize(3)
      write(6,*) zsize(1),zsize(2),zsize(3),zsize(1)*zsize(2)*zsize(3)
      endif

c********************************************************************************
c     IDENTIFYING PROCESSING NODES, DEFINING NEW COMMUNICATOR 
c********************************************************************************

      call mpi_get_processor_name(node_name, result_proc, ierr)

      color = xstart(2)
      call mpi_comm_split(MPI_COMM_WORLD, color, nrank, MPI_COMM_NODE,ierr)
      call mpi_comm_rank(MPI_COMM_NODE, node_rank, ierr)
      call mpi_comm_size(MPI_COMM_NODE, node_size, ierr)
      ! now all the nodes are devided into MPI_COMM_NODE with rank
      ! node_rank
      call mpi_barrier(MPI_COMM_WORLD, ierr)

c********************************************************************************
c     STARTING CODE
c********************************************************************************

      !unew,vnew,wnew are velocity in wall-normal,spanwise and streamwise
      !direction
      !phirnew,phitnew,phiznew correspond to rho*u,rho*v,rho*w 
      !enew is enthalpy and cnew is temperature which are equal since specific
      !heat is constant in my simulations.

      hwallh =1.0 !4./3.!1.6   !wall enthalpy
      hwallc =0.   !2./3.!1.!0.4
      Twallh =1.0 !4./3.!1.6   !wall temperature 
      Twallc =0.   !2./3.!1.!0.4
      Tc = 573.
      Re = 3750  !10000! = uplus*Dh/visc
      Lz = 16.*atan(1.)

      Tplus = 1.5
      dt =1e-6
      nstap = 999999 

      Rin = 0.0
      Rout =2.0 !half channel height is unity.
    
      iload = 1               !iload=0 is for starting from scratch

      visc = 1./Re
      Pr = 1.0
      Pl = 0.03
      call read_kP(kP,Tnb,nrank)     
 
      time = 0.
      ini = 0
      ini_sol=0
      Ich=0.
      divQ=0.
      qflux=0.
      qfluy=0.
      qfluz=0.
      G=0.      

      dpdx = 0.005

      open(11,file='icount')
      read(11,*)icount,time
      close(11)
   
      call mkgrid(Lz,nrank)


      if (iload.eq. 0)   call init (nrank,dpdx)
      if (iload.ne. 0)   call loadd(0,nrank,icount)

      call stateIG()
      do k=1,kmax/p_col
         do j=1,jmax/p_row
          do i=1,imax
            hir_s(i,j,k)=rhonew(i,j,k)
          enddo
        enddo
       enddo

      call inter_s_c_m(imax,jmax*kmax/(p_col*p_row),hir_s,hir_c,dr)

      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
            phirnew(i,j,k) = unew(i,j,k)*hir_c(i,j,k)
        enddo
       enddo
      enddo

!      do i=1,imax
!       cnew(i,:,:)=0.1 !(1-(rp(i)/2.)**()
!      enddo
!      enew=cnew

      renew = enew*rhonew
      phitnew = vnew*rhonew
      phiznew = wnew*rhonew
      rho_o=rhonew 
      rho_oo=rhonew 
      rho_ooo=rhonew 

      Twb = Twallh
      Twt = Twallc
      call bound(unew,vnew,wnew,nrank)
      call bound(phirnew,phitnew,phiznew,nrank)
      call bounds(enew,cnew,hwallh,hwallc,Twallh,Twallc,nrank)
      call output(0,rank,cpipe,dpdx)
      call chkdt(nrank)
      if (nrank.eq.0) call report

c***********************************************************************
c     Kickoff Radiative calculation
c***********************************************************************
      
      call fix_temp(cnew,TMC,Tc,Tplus,xstart,xend)
      if(node_rank.eq.0) then
        call MC_GPU(TMC, xstart(2))
      endif

c***********************************************************************
c     Time cycle
c***********************************************************************

      do istap=1,nstap

      stime1 = MPI_WTIME()
 
       if (mod(istap,10).eq. 0) then
        call cmpbs(bulk,stress,instress,Cbulk,fluxPipe,fluxAnnulus,Qval,cpipe,meanVisc,cannulus,opthck)
        dpdx = 0.98*dpdx - (bulk-1) * 50  
        Ck = Ck - (opthck-10)*0.2  
       endif

       if (mod(istap,10).eq. 0) then
c
c     Old levels of Adamsb-Bashforth or lost after restart. First and second step with '  ini = 0 '
c
      if (nrank.eq.0) write(6,111) time, Re*bulk , stress,instress, istap!dpdx! ,dr,0.5/imax
111   format (' tijd = ', f16.6,' RE_bulk  =  ',F16.6, ' Stress = ',1E15.5, ' instress = ', 1E15.5, ' dpdx = ', I6.6)!,2E15.5)
      if (nrank.eq.0) write(6,200) dt, Cbulk , fluxPipe, fluxAnnulus,opthck,Ck
200   format (' tijd = ' f16.6,' Tbulk    =  ',F16.6, ' FluxPi = ',1E15.5, ' FluxAnnu = ', 1E15.5,' opthck = '1E15.5,
     +        ' Ck = ',1E15.5)
       endif

      stime = MPI_WTIME()
     
c***********************************************************************
c     Radiation loop
c***********************************************************************

      stime2 = MPI_WTIME()
      if(mod(istap,5) == 0) then
        if(node_rank.eq.0) then
           call GET_RESULTS(qr, vr)       
        endif
        call mpi_barrier(MPI_COMM_NODE,ierr)       
        call fix_rad(divQ,vcore,qr,vr,Tc,MPI_COMM_NODE,xstart,xend)
        call fix_kappa(kP,Tnb)
        G = 4*kappa*(cnew/Tplus+1)**4. - divQ
        call fix_temp(cnew,TMC,Tc,Tplus,xstart,xend)       
        call mpi_barrier(MPI_COMM_WORLD,ierr)       
        if(node_rank.eq.0) then
           call MC_GPU(TMC, xstart(2)) 
        endif
      endif
      stime3 = MPI_WTIME()

c***********************************************************************
c     End of radiation loop
c***********************************************************************
 
      call adamsb(ini,dpdx,rank)
      ini = 1
      if (mod(istap,1000).eq. 0) then
      if (rank.eq.0) write(6,*) ' Adamsb tijd = ', MPI_WTIME()-stime
      stime = MPI_WTIME()
      endif
      
      call bound(dudt,dvdt,dwdt,rank)
      call bound(phirnew,phitnew,phiznew,nrank)
!      call bounds(enew,cnew,hwallh,hwallc,Twallh,Twallc,nrank)
      call trunc(dudt,rank)
      call trunc(dvdt,rank)
      call trunc(dwdt,rank)
!      call trunc(renew,rank)
             if (mod(istap,1000).eq. 0) then
             if (rank.eq.0) write(6,*) ' trunc  tijd = ', MPI_WTIME()-stime
             stime = MPI_WTIME()
             endif
!      rhonew=rho_ooo     
      call fillps(ini,dt,1,istap)
     
             if (mod(istap,1000).eq. 0) then
             if (rank.eq.0) write(6,*) ' fillps tijd = ', MPI_WTIME()-stime
             stime = MPI_WTIME()
             endif 

       if (ipois.eq.2) then 
         call solver_2(ini_sol,p,ru,rp,dr,dtheta,dz,mr_s,mr_c,nrank)
       endif
       if (ipois.ne.2) then 
         call solver_com(ini_sol,p,ru,rp,dr,dtheta,dz,mr_s,mr_c,nrank)
         ini_sol =10
       endif
      
            if (mod(istap,1000).eq. 0) then
            if (rank.eq.0) write(6,*) ' Solver tijd = ', MPI_WTIME()-stime
            stime = MPI_WTIME()
            endif
 
       call correc()
!       enew = renew/rho_star
!       call stateIG()    
         
             if (mod(istap,1000).eq. 0) then
             if (rank.eq.0) write(6,*) ' correc tijd = ', MPI_WTIME()-stime
             endif
       
      call bound(unew,vnew,wnew,rank)
      call bound(phirnew,phitnew,phiznew,rank)
      call bounds(enew,cnew,hwallh,hwallc,Twallh,Twallc,nrank)
      
      if (mod(istap,100).eq. 0) call chkdiv(istap,rank)
      if (mod(istap, 100).eq. 0) call chkdt (rank)

      if (rank.eq.0     .and.  mod(istap,1000).eq.0  ) write(6,*) ' tijd     =    ',time,dt
456   format( ' The turbulent energy | stress | bulk | dp/dx  = ',F13.6 ,' | ', F13.6, ' | ', F13.6, '|', F13.6)
      call mpi_barrier(MPI_COMM_WORLD,ierr)

      if (mod(istap,5).eq.0) call output(istap/500,rank,cpipe,dpdx)
      if (mod(istap,1500).eq.0) then
      icount = icount + 1
      call loadd(1,nrank,icount)
      if (nrank.eq.0) then
      open(11,file='icount')
      write(11,*) icount,time
      close(11)
      endif
      endif
      

      time = time + dt
      stime4 = MPI_WTIME(); 
      if(nrank.eq.0 .and.mod(istap,1).eq.0) write(6,117)istap,stime4-stime1,stime3-stime2
117   format('Iteration:',i5,' Cyle time: ',f10.5,' Radiation time: ',f10.5)

      enddo
      call decomp_2d_finalize
      call mpi_finalize(ierr)
      stop
      end
 

      subroutine loadd(ini,nnrank,istap)
      use decomp_2d
      use decomp_2d_io
      include 'param.txt'
      include 'common.txt'
      real utmp(0:i1,jmax/p_row,kmax/p_col)
      real vtmp(0:i1,jmax/p_row,kmax/p_col)
      real wtmp(0:i1,jmax/p_row,kmax/p_col)
      real ptmp(0:i1,jmax/p_row,kmax/p_col)
      real ctmp(0:i1,jmax/p_row,kmax/p_col),yplus
      
      integer istap
      character*5 cha
      character*5 cha2
      call cnvstr(nnrank,cha)
      call cnvstr(istap,cha2)

       do i=1,imax
       ptmp(i,:,:)=p(i,:,:)
       enddo

 
      if (ini.eq.0) then
      call decomp_2d_read_one(1,unew,'DNS/u.'//cha2//'.dns')
      call decomp_2d_read_one(1,vnew,'DNS/v.'//cha2//'.dns')
      call decomp_2d_read_one(1,wnew,'DNS/w.'//cha2//'.dns')
      call decomp_2d_read_one(1,cnew,'DNS/c.'//cha2//'.dns')
!      call decomp_2d_read_one(1,enew,'DNS/e.'//cha2//'.dns')
!      call decomp_2d_read_one(1,rhonew,'DNS/rho.'//cha2//'.dns')
!      call decomp_2d_read_one(1,munew,'DNS/mu.'//cha2//'.dns')
!      call decomp_2d_read_one(1,lambdanew,'DNS/lambda.'//cha2//'.dns')
!      call decomp_2d_read_one(1,p,'DNS/p.'//cha2//'.dns')
!       wnew=1.002*wnew 
       enew=cnew      
!       do i=1,imax    
!        yplus=min((2.-rp(i))*Re,(rp(i))*Re)
!           if  (yplus .lt. 11.6)  enew(i,:,:)= (yplus)*4/Re+1.
!           if  (yplus .gt. 11.6)  enew(i,:,:)= (2.5*log(yplus)+6.)*4/Re+1
!   
!           cnew(i,:,:)=enew(i,:,:)
!      enddo
      endif

      if (ini.eq.1) then
      call decomp_2d_write_one(1,unew,'DNS/u.'//cha2//'.dns')
      call decomp_2d_write_one(1,vnew,'DNS/v.'//cha2//'.dns')
      call decomp_2d_write_one(1,wnew,'DNS/w.'//cha2//'.dns')
      call decomp_2d_write_one(1,cnew,'DNS/c.'//cha2//'.dns')
      call decomp_2d_write_one(1,divQ,'DNS/q.'//cha2//'.dns')
      call decomp_2d_write_one(1,kappa,'DNS/k.'//cha2//'.dns')
      call decomp_2d_write_one(1,G,'DNS/G.'//cha2//'.dns')
      call decomp_2d_write_one(1,qflux,'DNS/fx.'//cha2//'.dns')
      call decomp_2d_write_one(1,qfluy,'DNS/fy.'//cha2//'.dns')
      call decomp_2d_write_one(1,qfluz,'DNS/fz.'//cha2//'.dns')
!      call decomp_2d_write_one(1,enew,'DNS/e.'//cha2//'.dns')
!      call decomp_2d_write_one(1,rhonew,'DNS/rho.'//cha2//'.dns')
!      call decomp_2d_write_one(1,munew,'DNS/mu.'//cha2//'.dns')
!      call decomp_2d_write_one(1,lambdanew,'DNS/lambda.'//cha2//'.dns')
      call decomp_2d_write_one(1,ptmp,'DNS/p.'//cha2//'.dns')
      endif

      if (ini.eq.1) write(6,*) nrank
      end



      subroutine output(istap,rank,cpipe,dpdx)
      use decomp_2d
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer istap
      character*5 cha
      character*5 cha2
      real um(imax),vm(imax),wm(imax),cm(imax),em(imax),phiwm(imax),rm(imax),dpdzm(imax),mum(imax),lm(imax)
      real umm(imax),vmm(imax),wmm(imax),cmm(imax),emm(imax),phiwmm(imax),rmm(imax),dpdzmm(imax),mumm(imax),lmm(imax)
      real ur(imax),vr(imax),wr(imax),cr(imax),uw(imax),hi_s(imax),hi_c(0:imax),totstressm(imax),uphiwr(imax),ue(imax)
      real wmdiff(imax),wmmdiff(imax),hi_s2(imax,jmax/p_row,kmax/p_col),hi_c2(0:imax,jmax/p_row,kmax/p_col)
      real hi_s3(imax),hi_c3(0:imax),cpm(imax),cpmm(imax)
      real vcm(imax)
      real hi_s4(imax,jmax/p_row,kmax/p_col),hi_c4(0:imax,jmax/p_row,kmax/p_col)
      real mudwm(imax),mudwmm(imax),w_fluc(imax,jmax/p_row,kmax/p_col) 
      real lambdadTm(imax),lambdadTmm(imax),Gm(imax),Gr(imax) 
      real integral,dpdx,totflux(imax),integral2
      real fz (kmax), dfz(kmax)
      real p1x(0:imax+1,jmax/p_row,kmax/p_col)
      real p1y(0:imx   ,jmax      ,kmax/p_col)
      real p1z(0:imx   ,jmax/p_col,kmax      )
      real qm(imax),km(imax),fxm(imax),fym(imax),fzm(imax),dfxm(imax),dfxm_c(0:imax)
      real qr(imax),kr(imax),Rdif

      Rdif = 1./(Re*Pr*Pl)
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
         p1x(i,j,k)=p(i,j,k)
        enddo 
       enddo
      enddo

      call transpose_x_to_y(p1x,p1y)
      call transpose_y_to_z(p1y,p1z)
      do j=1,jmax/p_col
       do i=0,imx
        do k=1,kmax
         fz(k)=p1z(i,j,k)
        enddo
        call four1(kmax,fz,dfz,dz)
        do k=1,kmax
         p1z(i,j,k)=dfz(k)
        enddo
       enddo
      enddo 
      call transpose_z_to_y(p1z,p1y)
      call transpose_y_to_x(p1y,p1x)

      fxm=0
      fym=0
      fzm=0
      qm=0
      Gm=0
      qr=0
      km=0
      kr=0
      um=0
      vm=0
      wm=0
      vcm=0
      wmdiff=0
      wmmdiff=0
      cm=0
      em=0
      rm=0
      mum=0
      lm=0
      phiwm=0
      dpdzm=0
      ur=0
      Gr=0
      vr=0
      wr=0
      cr=0
      uw=0
      ue=0
      uphiwr=0
      mudwm = 0
      lambdadTm = 0
      cpm = 0
      do k=1,kmax/p_col
       do j=1,jmax/p_row
         do i=1,imax
              um(i)=um(i)+phirnew(i,j,k)/(jmax*kmax)
              vm(i)=vm(i)+phitnew(i,j,k)/(jmax*kmax)
              vcm(i)=vcm(i)+vcore(i,j,k)/(jmax*kmax)
              wm(i)=wm(i)+phiznew(i,j,k)/(jmax*kmax)  !Favre averaging
              cm(i)=cm(i)+cnew(i,j,k)/(jmax*kmax)
              em(i)=em(i)+renew(i,j,k)/(jmax*kmax)
              rm(i)=rm(i)+rhonew(i,j,k)/(jmax*kmax)
              mum(i)=mum(i)+munew(i,j,k)/(jmax*kmax)
              lm(i)=lm(i)+lambdanew(i,j,k)/(jmax*kmax)
              phiwm(i)=phiwm(i)+rhonew(i,j,k)*wnew(i,j,k)/(jmax*kmax)
              dpdzm(i)=dpdzm(i)+p1x(i,j,k)/(jmax*kmax)
              cpm(i)=cpm(i)+cpnew(i,j,k)/(jmax*kmax)
              qm(i)=qm(i)+divQ(i,j,k)/(jmax*kmax)
              Gm(i)=Gm(i)+G(i,j,k)/(jmax*kmax)
              km(i)=km(i)+kappa(i,j,k)/(jmax*kmax)
              fxm(i)=fxm(i)+qflux(i,j,k)/(jmax*kmax)
              fym(i)=fym(i)+qfluy(i,j,k)/(jmax*kmax)
              fzm(i)=fzm(i)+qfluz(i,j,k)/(jmax*kmax)
           enddo
         enddo
      enddo
      call mpi_allreduce(um,umm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(vm,vmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(wm,wmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(cm,cmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(em,emm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(rm,rmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(mum,mumm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(lm,lmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(phiwm,phiwmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(dpdzm,dpdzmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(cpm,cpmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)

      um = umm/rmm
      vm = vmm/rmm
      wm = wmm/rmm
      cm = cmm
      em = emm/rmm
      rm = rmm
      mum = mumm
      lm = lmm
      cpm = cpmm
      phiwm = phiwmm
      dpdzm = dpdzmm
      call mpi_allreduce(qm,cmm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      qm = cmm
      call mpi_allreduce(km,cmm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      km = cmm
      call mpi_allreduce(fzm,cmm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      fzm = cmm
      call mpi_allreduce(fym,cmm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      fym = cmm
      call mpi_allreduce(vcm,cmm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      vcm = cmm
      call mpi_allreduce(Gm,cmm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      Gm = cmm

      do k=1,kmax/p_col                         !Calculate derivative in r of velocity fluctuations
       do j=1,jmax/p_row
        do i=1,imax
         w_fluc(i,j,k) = wnew(i,j,k)-wm(i)!wmdiff(i)
         hi_s2(i,j,k) = w_fluc(i,j,k)
         hi_s4(i,j,k) = cnew(i,j,k) - cm(i)
        enddo
       enddo
      enddo
      call der1w_s_c6_m(imax,jmax*kmax/px,hi_s2,hi_c2,dr)
      call deriv_s_c_m(imax,jmax*kmax/px,hi_s4,hi_c4,dr)
      do k=1,kmax/p_col
       do j=1,jmax/p_row 
        do i=1,imax
         hi_c2(i,j,k) = hi_c2(i,j,k) * mr_c(i)
         hi_c4(i,j,k) = hi_c4(i,j,k) * mr_c(i)
        enddo
       enddo
      enddo
      call inter_c_s_m(imax,jmax*kmax/px,hi_c2,hi_s2,dr)
      call inter_c_s_m(imax,jmax*kmax/px,hi_c4,hi_s4,dr)

      do k=1,kmax/p_col
        do j=1,jmax/p_row
          do i=1,imax
            ur(i)=ur(i)+(unew(i,j,k)-um(i))**2
            vr(i)=vr(i)+(vnew(i,j,k)-vm(i))**2
            wr(i)=wr(i)+(wnew(i,j,k)-wm(i))**2
            cr(i)=cr(i)+(cnew(i,j,k)-cm(i))**2
            qr(i)=qr(i)+(divQ(i,j,k)-qm(i))**2
            Gr(i)=Gr(i)+(G(i,j,k)-Gm(i))**2
            kr(i)=kr(i)+(kappa(i,j,k)-km(i))**2
            uw(i)=uw(i)+rhonew(i,j,k)*(wnew(i,j,k)-wm(i))*(unew(i,j,k)+Unew(i-1,j,k)-um(i)-um(i-1))*.5
            ue(i)=ue(i)+rhonew(i,j,k)*(enew(i,j,k)-em(i))*(unew(i,j,k)+Unew(i-1,j,k)-um(i)-um(i-1))*.5
            mudwm(i) = mudwm(i) + (munew(i,j,k)-mum(i))*hi_s2(i,j,k)
            lambdadTm(i) = lambdadTm(i) + (lambdanew(i,j,k)-lm(i))*hi_s4(i,j,k)
           enddo
        enddo
      enddo
      call mpi_allreduce(ur,umm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(vr,vmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(wr,wmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(cr,cmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      ur = umm
      vr = vmm
      wr = wmm
      cr = cmm
      call mpi_allreduce(qr,cmm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      qr = cmm
      call mpi_allreduce(kr,cmm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      kr = cmm
      call mpi_allreduce(Gr,cmm,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      Gr = cmm
      call mpi_allreduce(uw,wmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(ue,cmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(uphiwr,phiwmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(mudwm,mudwmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(lambdadTm,lambdadTmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      uw = wmm 
      ue = cmm 
      uphiwr = phiwmm 
      mudwm = mudwmm
      lambdadTm = lambdadTmm

      hi_s = wm
      call der1w_s_c6(imax,hi_s,hi_c,dr)
      hi_c = hi_c * mr_c
      call inter_c_s(imax,hi_c,hi_s,dr)
      hi_s3 = cm
      call der1c_s_c6_m(imax,1,hi_s3,hi_c3,dr,Twallh,Twallc)
      hi_c3 = hi_c3 * mr_c
      call inter_c_s_m(imax,1,hi_c3,hi_s3,dr)
      call deriv_s_c_m(imax,1,fxm,dfxm_c,dr)
      dfxm_c=dfxm_c*mr_c
      call inter_c_s_m(imax,1,dfxm_c,dfxm,dr)
  
      do i=1,imax
       integral = integral - (-(dpdx+dpdzm(i))*rp(i)+1./(1*0.01)*4*0.224*rm(i)*rp(i))*(ru(i)-ru(i-1))
       integral2 = integral2 - 4*((wm(i)*rp(i))*(ru(i)-ru(i-1)))
       totstressm(i) = integral/rp(i)
       totflux(i) = integral2/rp(i)
      enddo  

      if (rank.eq.0) then
      open(17,file = 'prof')
      open(18,file = 'prof_log')
      open(19,file = 'profTemp')
      open(20,file = 'profProp')
      do i=1,imax
       write(17,112) rp(i),um(i),vm(i),wm(i),
     1 sqrt(ur(i)/(jmax*kmax)),sqrt(vr(i)/(jmax*kmax)),sqrt(wr(i)/(jmax*kmax)),(uw(i)/(jmax*kmax)),
     2 -mum(i)*hi_s(i),(uw(i)/(jmax*kmax))-mum(i)*hi_s(i),mudwm(i)/(jmax*kmax)!,uphiwr(i)/(jmax*kmax),uphiwr(i)/(jmax*kmax)-visc*hi_c(i),
c     3 totstressm(i)
      enddo
      do i=imax,1,-1
       write(18,112) (0.5-rp(i))/visc,um(i),vm(i),wm(i),
     1 sqrt(ur(i)/(jmax*kmax)),sqrt(vr(i)/(jmax*kmax)),sqrt(wr(i)/(jmax*kmax)),(uw(i)/(jmax*kmax))
      enddo
      do i=1,imax
        write(19,112) rp(i),em(i),cm(i),sqrt(cr(i)/(jmax*kmax)),ue(i)/jmax/kmax,-lm(i)*hi_s3(i)
     ^               ,ue(i)/(jmax*kmax)-lm(i)*hi_s3(i),lambdadTm(i)/(jmax*kmax),rm(i)*cpm(i)*wm(i)*(ru(i)-ru(i-1))/(lm(i))
      enddo
      do i=1,imax
        write(20,112) rp(i),em(i),cm(i),rm(i),mum(i)*Re,lm(i)*Re*Pr,cpm(i),
     +                qm(i),sqrt(qr(i)/(jmax*kmax)),vcm(i),Gm(i),sqrt(Gr(i)/(jmax*kmax)),
     +                km(i),sqrt(kr(i)/(jmax*kmax)),
     +                fxm(i),fym(i),fzm(i),dfxm(i)
      enddo
      close(17)
      close(18)
      close(19)
      close(20)

112   FORMAT(19E16.5) 
      endif
      end

      subroutine fillps(ini,dtt,switch,istap)
c
c  right hand side of the poisson equation
c
c       1 dru*    1 dv*    dw* 
c       - --- +   - ---  + --- = 0
c       r d r     r d t    d z
c
c
c         (1)      (2)      (3)
c
      use decomp_2d
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer im,ier,idex,i_dex,ini,switch,istap
      real ft(jmax),dft(jmax),fz(kmax),dfz(kmax),som,stime
      real hi_s(imax,jmax/p_row,kmax/p_col),hi_c(0:imax,jmax/p_row,kmax/p_col),rpi(0:imax+1),rui(0:imax+1)
      real tmp (0:imx,jmax,kmax/p_col)
      real tmp2(0:imx,jmax,kmax/p_col)
      real tmp3(0:imx,jmax/p_col,kmax)
      real p1(0:imax+1,jmax/p_row,kmax/p_col)
      real vol1, pcorr1, vol2, pcorr, pcorr_o, pcorr_oo, pcorr_ooo
      real dtt,mindrt 
!      idex = i_dex(nrank)
!      rpi = 1./rp
!      rui = 1./ru
!      p = 0
      pcorr1 = 0 
      vol1 = 0
      mindrt=10.

      if (switch.eq.0) then 
      drhodt =1./(2*dt)*(3*rhonew-4*rho_o+rho_oo)!1./(2*dt)*(3*rho_star-4*rhonew+rho_o)!1./(dt)*(rho_star-rhonew)
      endif
      if (switch.eq.1) then
      drhodt =1./(2*dt)*(3*rhonew-4*rho_o+rho_oo)!1./(2*dt)*(3*rho_star-4*rhonew+rho_o)!1./(dt)*(rho_star-rhonew)
      endif
      call trunc(drhodt,nrank)

      call transpose_x_to_y(dvdt,tmp) 

      do i=0,imx
	    do k=1,kmax/p_col
	  do j=1,jmax
	    ft(j)=tmp(i,j,k)
	  enddo
	  call four1(jmax,ft,dft,dtheta)     !drhov/dy
	  do j=1,jmax
	   tmp(i,j,k)=dft(j)!*rpi(i+idex)
	  enddo
	    enddo
      enddo
      call transpose_x_to_y(dwdt,tmp2)
      call transpose_y_to_z(tmp2,tmp3)
c term (1)
      do i=0,imx
	    do j=1,jmax/p_col
	   do k=1,kmax
	     fz(k)=tmp3(i,j,k)
	   enddo
	   call four1(kmax,fz,dfz,dz)       !drho w/dz
	   do k=1,kmax
	    tmp3(i,j,k)=dfz(k)
	   enddo
	     enddo
	    enddo
        call transpose_z_to_y(tmp3,tmp2)
        tmp = tmp + tmp2 
        call transpose_y_to_x(tmp,p1)

         if (ipois.ne.2) then        
         do k=1,kmax/p_col
	  do j=1,jmax/p_row
	   do i=0,imax
!            hi_c(i,j,k)=ru(i)*dudt(i,j,k)
             hi_c(i,j,k)=dudt(i,j,k)
           enddo
          enddo
         enddo
           call deriv_c_s_m(imax,jmax*kmax/px,hi_c,hi_s,dr)  !drho u/dx
         do k=1,kmax/p_col
	    do j=1,jmax/p_row
           do i=1,imax
	     p(i,j,k)=p1(i,j,k)+hi_s(i,j,k)*mr_s(i)!*rpi(i)     !drho u/dr
             pcorr1 = pcorr1-(p(i,j,k))*(ru(i)-ru(i-1))*dtheta*dz!*rp(i)   !wronggg
             vol1 = vol1+(ru(i)-ru(i-1))*dtheta*dz!*rp(i)
           enddo
         enddo
         enddo
      endif

         if (ipois.eq.2) then        
         
         do k=1,kmax/p_col
	    do j=1,jmax/p_row
           do i=1,imax

!	     p(i,j,k)=p1(i,j,k)+(ru(i)*dudt(i,j,k)-ru(i-1)*dudt(i-1,j,k))/(rp(i)*(ru(i)-ru(i-1)))
         p(i,j,k)=p1(i,j,k)+(dudt(i,j,k)-dudt(i-1,j,k))/((ru(i)-ru(i-1)))
         pcorr1 = pcorr1+(p(i,j,k)+drhodt(i,j,k))*(ru(i)-ru(i-1))*dtheta*dz
             vol1 = vol1+(ru(i)-ru(i-1))*dtheta*dz!*rp(i)
         mindrt=min(mindrt,abs(drhodt(i,j,k))) 
          enddo
         enddo
         enddo
      endif

!      do k=1,kmax/p_col
!       do j=1,jmax/p_row
!        do i=1,imax
!         pcorr1 = pcorr1+(p(i,j,k)+drhodt(i,j,k))*(ru(i)-ru(i-1))*dtheta*dz     !*rp(i)
!        enddo
!       enddo
!      enddo

      call mpi_allreduce(pcorr1  ,pcorr      ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(vol1  ,vol2      ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      if (mod(istap,100).eq. 0) then
      if (nrank.eq.0) write(*,*) ' press  corr   =',(pcorr)/vol2,mindrt
      endif      

      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
         drhodt(i,j,k) = drhodt(i,j,k) - pcorr/vol2
c         if (switch.eq.0) drhodt(i,j,k) = drhodt(i,j,k) - pcorr/vol
c         if (switch.eq.1) then
c         rhonew(i,j,k) = rhonew(i,j,k) -2./3.*pcorr*dt/vol
c         drhodt = 1./(2*dt)*(3*rhonew-4*rho_o+rho_oo)
c         endif
        enddo 
       enddo
      enddo

       do k=1,kmax/p_col
        do j=1,jmax/p_row
         do i=1,imax
          p(i,j,k) = p(i,j,k)+drhodt(i,j,k)
         enddo
        enddo
       enddo

       p = p /dt

       if (switch.eq.1) then
       rho_ooo = rho_oo
       rho_oo = rho_o
       rho_o = rhonew

       endif
 
       end
     
     
      subroutine chkdiv(istap,rank)
      use decomp_2d
c    calculates divergence
c
c
c       1 dru     1 d v    d w 
c       - --- +   - ---  + --- = 0
c       r  dr     r dt     dz
c
c term    (1)      (2)     (3)
c

      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'

      integer im,ier,idex,i_dex,ierr,istap
      real ft(jmax),dft(jmax),fz(kmax),dfz(kmax),som,dmaxx!,div(0:i1,jmax/p_row,kmax/p_col)
      real t_dmaxx
      real t_som,hi_c(0:imax),hi_s(imax)
      real tmpx(0:i1 ,jmax/p_row,kmax/p_col)         
      real tmpy(0:imx,jmax      ,kmax/p_col)         
      real tmpz(0:imx,jmax/p_col,kmax      )    
!      idex = i_dex(nrank)
      div = 0.
      dmaxx =0.

      call transpose_x_to_y(phiznew,tmpy)
      call transpose_y_to_z(tmpy,tmpz)
      do i=0,imx
	    do j=1,jmax/p_col
	   do k=1,kmax
	    fz(k) =tmpz(i,j,k)
	   enddo
	   call four1(kmax,fz,dfz,dz)   !d(rho w)/dz
	   do k=1,kmax
	    tmpz(i,j,k)=dfz(k)
	   enddo
	    enddo
      enddo
      call transpose_z_to_y(tmpz,tmpy)
      call transpose_y_to_x(tmpy,div)
     
      call transpose_x_to_y(phitnew,tmpy)
      do i=0,imx
	    do k=1,kmax/p_col
	  do j=1,jmax
	    ft(j)=tmpy(i,j,k)
	  enddo
	  call four1(jmax,ft,dft,dtheta)  !d(rho v)/dy
	  do j=1,jmax
	   tmpy(i,j,k)=dft(j)!/rp(i+idex)
	  enddo
	    enddo
      enddo
      call transpose_y_to_x(tmpy,tmpx)
      div = div  + tmpx
       som = 0 
c term (1)
       dmaxx = 0.
	   do j=1,jmax/p_row
	     do k=1,kmax/p_col
	    do i=0,imax
            hi_c(i)=phirnew(i,j,k)!*ru(i)
         enddo
         call deriv_c_s(imax,hi_c,hi_s,dr)  !d(rho u)/dx
         hi_s = hi_s*mr_s                   !d(rho u)/dr
          do i=1,imax
	  if (ipois.ne.2)    div(i,j,k)=div(i,j,k)+hi_s(i)! /rp(i)
!	  if (ipois.eq.2)   div(i,j,k)=div(i,j,k)+(ru(i)*phirnew(i,j,k)-ru(i-1)*phirnew(i-1,j,k))/(rp(i)*(ru(i)-ru(i-1)))
      if (ipois.eq.2)   div(i,j,k)=div(i,j,k)+(phirnew(i,j,k)-phirnew(i-1,j,k))/((ru(i)-ru(i-1)))
	    enddo
	     enddo
       enddo

       div_star_o = div_star       
       div_star = div       
       dmaxx = 0
       som =0
       do k=1,kmax/p_col
         do j=1,jmax/p_row
	   do i=1,imax
             div(i,j,k)=(div(i,j,k)+drhodt(i,j,k))*(ru(i)-ru(i-1))*dtheta*dz!*rp(i)
             dmaxx = max (dmaxx,abs(div(i,j,k)))
             som = som + div(i,j,k)
c             if (abs(div(i,j,k)).gt.1e-6) 
c              write(*,*) i,j,k,abs(div(i,j,k))
           enddo
         enddo
       enddo
       call mpi_allreduce(som,t_som    ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
       call mpi_allreduce(dmaxx,t_dmaxx,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,ier)
       if (mod(istap,100).eq. 0) then
       if (nrank.eq.0) write(*,*) 'divtot = ',t_som, 'divmax = ',t_dmaxx
       endif
       if (abs(t_dmaxx).gt.1e-0) then
       call mpi_finalize(ierr)
       stop 'chkdiv'
       endif 
       end
     
     



      subroutine bound(u,v,w,rank)
      include 'param.txt'
      include 'common.txt'
      real u(0:i1,jmax/p_row,kmax/p_col)
      real v(0:i1,jmax/p_row,kmax/p_col)
      real w(0:i1,jmax/p_row,kmax/p_col)
       do k=1,kmax/p_col
	    do j=1,jmax/p_row
	  u(imax,j,k)=0.
	  v(i1  ,j,k)=-v(imax,j,k)	
	  w(i1  ,j,k)=-w(imax,j,k)	

	  u(0   ,j,k)= 0!-u(1,j,k)
	  v(0   ,j,k)= -v(1,j,k)
	  w(0   ,j,k)= -w(1,j,k)

	    enddo
       enddo
       end


      subroutine bounds(e,c,hwh,hwc,Twh,Twc,rank)
      include 'param.txt'
      include 'common.txt'
      real c(0:i1,jmax/p_row,kmax/p_col)
      real e(0:i1,jmax/p_row,kmax/p_col)
      real Twh,Twc,hwh,hwc
      do k=1,kmax/p_col
        do j=1,jmax/p_row

!      e(i1,  j,k)=-e(imax,j,k)
!      c(i1,  j,k)=-c(imax,j,k)
      e(i1,  j,k)=(hwc-e(imax,j,k))/(ru(imax)-rp(imax))*(rp(i1) -rp(imax)) + e(imax,j,k)
      c(i1,  j,k)=(Twc-c(imax,j,k))/(ru(imax)-rp(imax))*(rp(i1) -rp(imax)) + c(imax,j,k)
      c(0   ,j,k)= (c(1,j,k)-Twh)/(rp(1))*(rp(0) - rp(1)) + c(1,j,k)
      e(0   ,j,k)= (e(1,j,k)-hwh)/(rp(1))*(rp(0) - rp(1)) + e(1,j,k)
        enddo
       enddo
       end


      subroutine correc()
      use decomp_2d
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      real som
      integer idex,i_dex,ier
      real fz(kmax),dfz(kmax)
      real ft(jmax),dft(jmax)
      real hi_s(  imax,jmax/p_row,kmax/p_col)
      real hi_c(0:imax,jmax/p_row,kmax/p_col)
      real rpi(0:imax+1),rui(0:imax+1)

      real p1x(0:imax+1,jmax/p_row,kmax/p_col)
      real p1y(0:imx   ,jmax      ,kmax/p_col)
      real p1z(0:imx   ,jmax/p_col,kmax      )

      real hir_s(  imax,jmax/p_row,kmax/p_col)
      real hir_c(0:imax,jmax/p_row,kmax/p_col)
      real bulk,dpdx
      
!      rpi = 1./rp
!      rui = 1./ru 
      do k=1,kmax/p_col
	    do j=1,jmax/p_row
	  do i=1,imax
            hi_s(i,j,k)=p(i,j,k)
            hir_s(i,j,k)=rhonew(i,j,k)!rho_star(i,j,k)
          enddo
        enddo
      enddo

      call inter_s_c_m(imax,jmax*kmax/(p_col*p_row),hir_s,hir_c,dr)
      call deriv_s_c_m(imax,jmax*kmax/px,hi_s,hi_c,dr)                  !dp/dx

      do k=1,kmax/p_col
        do j=1,jmax/p_row
          do i=1,imax-1 
	    if (ipois.ne.2) phirnew(i,j,k)=dudt(i,j,k)-hi_c(i,j,k)*dt*mr_c(i)       !dp/dr
	    if (ipois.eq.2) phirnew(i,j,k)=dudt(i,j,k)-dt*(p(i+1,j,k)-p(i,j,k))/(rp(i+1)-rp(i))
            if (ipois.ne.2) unew(i,j,k) = phirnew(i,j,k)/hir_c(i,j,k)
            if (ipois.eq.2) unew(i,j,k) = phirnew(i,j,k)/hir_c(i,j,k)!((rhonew(i+1,j,k)+rhonew(i,j,k))/2)
	  enddo
	    enddo
      enddo 

      do k=1,kmax/p_col
	    do j=1,jmax/p_row
	  do i=1,imax
	   p1x(i,j,k)=p(i,j,k)
          enddo
        enddo
      enddo

      call transpose_x_to_y(p1x,p1y)
      call transpose_y_to_z(p1y,p1z)
      do k=1,kmax/p_col
	    do i=0,imx
	   do j=1,jmax
	     ft(j)=p1y(i,j,k)
	   enddo
	   call four1(jmax,ft,dft,dtheta)                   !dp/dy
           do j=1,jmax
	     p1y(i,j,k)=dft(j)
           enddo
        enddo
      enddo
      call transpose_y_to_x(p1y,p1x)

      do k=1,kmax/p_col
	     do j=1,jmax/p_row
           do i=1,imax 
	   phitnew(i,j,k)=(dvdt(i,j,k)-p1x(i,j,k)*dt)  !*rpi(i)
	   enddo
	    enddo
      enddo 
      do k=1,kmax/p_col
         do j=1,jmax/p_row
	   do i=1,imax
	    p1x(i,j,k)=p(i,j,k)
           enddo
         enddo
      enddo

      do j=1,jmax/p_col
    	do i=0,imx
	  do k=1,kmax
	    fz(k)=p1z(i,j,k)
	  enddo
	  call four1(kmax,fz,dfz,dz)     !dp/dz
	  do k=1,kmax
	   p1z(i,j,k)=dfz(k)
	  enddo
    	enddo
      enddo
      call transpose_z_to_y(p1z,p1y)
      call transpose_y_to_x(p1y,p1x)

      do k=1,kmax/p_col
    	do j=1,jmax/p_row
         do i=1,imax
	  phiznew(i,j,k)=(dwdt(i,j,k)-p1x(i,j,k)*dt)!*1.0/bulk
    	 enddo
	    enddo
      enddo

      vnew = phitnew/rhonew!rho_star
      wnew = phiznew/rhonew!rho_star
      end
      
      subroutine init(rank,dpdx)
      include 'param.txt'
      include 'common.txt'
      real yplus 
      
      unew= 0
      vnew= 0
      wnew= 0.0
      cnew= 0.0
      enew= 0.0
      do k=1,kmax/p_col
        do j=1,jmax/p_row
          do i=1,imax
              yplus=min((Rout-rp(i))/visc,(rp(i)-Rin)/visc) 
              enew(i,j,k)=hwallh*(1.- rp(i)/2.)
            if  (yplus .lt. 11.6)  wnew(i,j,k)= (yplus)
            if  (yplus .gt. 11.6)  wnew(i,j,k)= (2.5*log(yplus)+5.5)
            vnew(i,j,k)=3*sin(1.*j*j*k)
            unew(i,j,k)=0 
            enddo
          enddo
       enddo
c      wnew = 0.04917634*wnew*1.024
c      vnew = vnew/20
      call stateIG()
      renew = rhonew*enew
      phirnew=unew*rhonew
      phitnew=vnew*rhonew
      phiznew=wnew*rhonew
      end  

      subroutine mkgrid(Lz,rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      real const,Lz,rp_c(0:imax),rp_s(imax),x,dx,rnorm,rmax
      dz =Lz/(kmax)
      dr = (Rout-Rin)/imax
      dtheta =8.*atan(1.)/jmax
       rmax = Rout-Rin 
       ru(0)=0. 
!       do i=1,imax+1
!        x = i*1./imax
!    	dx= 2*x-1.8*x**2!1-3*(x-0.55)**2
!       ru(i)=ru(i-1)+dx
!       enddo
!       rnorm = Ru(0)
!       do i=0,imax+1
!       ru(i)=ru(i)-rnorm
!       enddo
!       ru(0)=0.!1e-49
!       rnorm = Ru(imax)
!       do i=0,imax+1
!	  ru(i)=ru(i)/(2*rnorm)
!       enddo
!       do i=0,imax+1
!         ru(i)=ru(i)+Rin
!       enddo

!      ru(0) = Rin
!        do i=1,imax/2
!          x  = (i*2./imax)
!          dx = dr!tanh(x*5.) + 1.
!          ru(i)=ru(i-1)+dx
!        enddo
!        rnorm = ru(imax/2)
!        do i=1,imax/2
!          ru(i)=Rout/2.*ru(i)/rnorm
!        enddo
!        do i=imax,imax/2+1,-1
!          ru(i)=Rout-ru(imax-i)
!        enddo

!       do i=1,imax/2
!         x  = i
!         dx = 0.00145*x**3. -5.966*x**2. +734.37*x +1493.81
!         ru(i)=ru(i-1)+dx
!       enddo
!       rnorm = ru(imax/2)
!       do i=1,imax/2
!         ru(i)=Rout/2.*ru(i)/rnorm
!       enddo
!       do i=imax,imax/2+1,-1
!         ru(i)=Rout-ru(imax-i)
!       enddo
      
      do i=1,imax/2
          x  = 1.*i/imax
          dx = 0.5-1.45*(x-0.5)**2.!0.00145*x**3. -5.966*x**2. +734.37*x+1493.81
          ru(i)=ru(i-1)+dx
        enddo
        rnorm = ru(imax/2)
        do i=1,imax/2
          ru(i)=Rout/2.*ru(i)/rnorm
        enddo
       do i=imax,imax/2+1,-1
         ru(i)=Rout-ru(imax-i)
       enddo


      do i=1,imax
        rp(i)=0.5*(ru(i)+ru(i-1))
        delta(i)=ru(i)-ru(i-1)
      enddo


!      do i=1,imax
!        rp(i)=0.5*(ru(i)+ru(i-1))
!        delta(i)=ru(i)-ru(i-1)
!      enddo
c      rp(0 )=ru(0)-rp(1)
      rp(0) = Rin - (rp(1)-Rin)
      rp(i1)=ru(imax)+(Ru(imax)-rp(imax))
c      write(*,*) 'hello  ',rp(0), ru(0), rp(1)
      do i=1,imax
	  rp_s(i)=rp(i)
      enddo
      call deriv_s_c(imax,rp_s,mr_c,dr)
      call inter_c_s(imax,mr_c,mr_s,dr)
      mr_c = 1./mr_c
      mr_s = 1./mr_s
!      if (rank.eq.0) then
!      do i=1,imax
!      if (rank.eq.0)                write(6,111) i,(rp(i)-Rin),(Rout-rp(i)), (Ru(i)-Ru(i-1))   !outer wall y+, dy+
!      enddo
!      endif
!111   format ('Grid node =  ',i5, ' y+an =', F15.5, ' y+Pipe  =  ', F15.5, ' d y+  =  ', F15.5)
      if (rank.eq.0) then
      open(11,file = 'grid.txt')
      write(11,*) Re,Ru(imax)
      do i=1,imax
         write(11,'(i5,4F12.6)') i,Ru(i),Rp(i),delta(i),delta(i)*395
      enddo
      endif
      close(11)

      end 



      subroutine adamsb(ini,dpdx,nrank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'

      real dfr_n1(0:imax+1,jmax/p_row,kmax/p_col)
      real dft_n1(0:imax+1,jmax/p_row,kmax/p_col)
      real dfz_n1(0:imax+1,jmax/p_row,kmax/p_col)
      real dfc_n1(0:imax+1,jmax/p_row,kmax/p_col)
      real dfr_n(0:imax+1,jmax/p_row,kmax/p_col)
      real dft_n(0:imax+1,jmax/p_row,kmax/p_col)
      real dfz_n(0:imax+1,jmax/p_row,kmax/p_col)
      real dfc_n(0:imax+1,jmax/p_row,kmax/p_col)
      real dfr_o(0:imax+1,jmax/p_row,kmax/p_col)
      real dft_o(0:imax+1,jmax/p_row,kmax/p_col)
      real dfz_o(0:imax+1,jmax/p_row,kmax/p_col)
      real dfc_o(0:imax+1,jmax/p_row,kmax/p_col)
      real dfr_oo(0:imax+1,jmax/p_row,kmax/p_col)
      real dft_oo(0:imax+1,jmax/p_row,kmax/p_col)
      real dfz_oo(0:imax+1,jmax/p_row,kmax/p_col)
      real dfc_oo(0:imax+1,jmax/p_row,kmax/p_col)
 
      real uint(0:imax+1,jmax/p_row,kmax/p_col)
      real vint(0:imax+1,jmax/p_row,kmax/p_col)
      real wint(0:imax+1,jmax/p_row,kmax/p_col)
      real eint(0:imax+1,jmax/p_row,kmax/p_col)
      real cint(0:imax+1,jmax/p_row,kmax/p_col)
      real phirint(0:imax+1,jmax/p_row,kmax/p_col)
      real phitint(0:imax+1,jmax/p_row,kmax/p_col)
      real phizint(0:imax+1,jmax/p_row,kmax/p_col)
      real reint(0:imax+1,jmax/p_row,kmax/p_col)
      real rhoint(0:imax+1,jmax/p_row,kmax/p_col)

      real u_t(0:imx,jmax,kmax/p_col) 
      real v_t(0:imx,jmax,kmax/p_col) 
      real w_t(0:imx,jmax,kmax/p_col) 
      real c_t(0:imx,jmax,kmax/p_col) 
      real e_t(0:imx,jmax,kmax/p_col) 
      real rho_t(0:imx,jmax,kmax/p_col) 
      real lambda_t(0:imx,jmax,kmax/p_col) 
      real mu_t(0:imx,jmax,kmax/p_col) 

      real ww1(0:imx,jmax,kmax/p_col) 
      real ww2(0:imx,jmax,kmax/p_col) 
      real ww3(0:imx,jmax,kmax/p_col) 
      real cc(0:imx,jmax,kmax/p_col) 
c      real iets(0:imax+1,jmax/p_row,kmax/p_col)
      real Fri,cbulk,bulk,source
      real meanVisc,dpdx 
 
      common /ab3/dfr_n1,dft_n1,dfz_n1,dfc_n1,dfr_n,dft_n,dfz_n,dfc_n,dfr_o,dft_o,dfz_o,dfc_o,dfr_oo,dft_oo,dfz_oo,dfc_oo
      integer order,ini,ini_sol,nrank
!      stime = MPI_WTIME()
      
      Fri = 1.0!0.0288/(65*0.098704*0.00306)*0.672
      
!     uint = unew
!     vint = vnew
!     cint = cnew
!     eint = enew
!     wint = wnew
      phirint = phirnew
      phitint = phitnew
      phizint = phiznew
      reint = renew
      rhoint=rhonew
!      source = 1.0/(Re*Pr)/bulk
!     fac1 = 23./12.
!     fac2 =-4./3.
!     fac3 = 5./12.

!     if (ini.eq.0) then
!     dfr_oo = dfr_n
!     dft_oo = dft_n
!     dfz_oo = dfz_n
!     dfc_oo = dfc_n
!     dfr_o = dfr_n
!     dft_o = dft_n
!     dfz_o = dfz_n
!     dfc_o = dfc_n
!      drhodt_o = drhodt
!     endif
      if (ini.eq.0) then
      dfr_oo = 0.!dfr_n
      dft_oo = 0.!dft_n
      dfz_oo = 0.!dfz_n
      dfc_oo = 0.!dfc_n
      dfr_o =dfr_n
      dft_o =dft_n
      dfz_o =dfz_n
      dfc_o =dfc_n
      rho_oo = rhonew
      rho_o = rhonew
!       drhodt_o = drhodt
      endif


      call momz(dfr_n,dft_n,dfz_n,dfc_n,ww1,ww2,ww3,cc,unew,vnew,wnew,cnew,
     ^ u_t,v_t,w_t,c_t,imax,jmax,kmax,imx,p_row,p_col,ru,rp,dr,dtheta,dz,nrank,px,rhonew,rho_t,munew,mu_t,enew,e_t,
     ^ lambdanew,lambda_t,mr_s,mr_c)
      call momy(dfr_n,dft_n,dfz_n,dfc_n,ww1,ww2,ww3,cc,unew,vnew,wnew,cnew,
     ^ u_t,v_t,w_t,c_t,imax,jmax,kmax,imx,p_row,p_col,ru,rp,dr,dtheta,dz,nrank,px,rhonew,rho_t,munew,mu_t,
     ^ enew,e_t,lambdanew,lambda_t,mr_s,mr_c)
      call momx(dfr_n,dft_n,dfz_n,dfc_n,unew,vnew,wnew,cnew,
     ^ imax,jmax,kmax,imx,p_row,p_col,ru,rp,dr,dtheta,dz,mr_s,mr_c,nrank,px,munew,rhonew,lambdanew,enew,Twallh,Twallc)

      call divTrace(rank)
      dfz_n = dfz_n - divdivW + dpdx!- divdivW! - Fri*rhonew 
      dft_n = dft_n - divdivV 
      dfr_n = dfr_n - divdivU 
      dfc_n = dfc_n - divQ/(Re*Pr*Pl) 

      renew = reint  + dt * (3./2.*dfc_n-1./2.*dfc_o+0*dfc_oo)
      call trunc(renew,rank)
      rho_star = rhonew + dt*( drhodt)
      enew = renew/rho_star

      dudt = phirint + dt * (3./2.*dfr_n-1./2.*dfr_o+0*dfr_oo)
      dvdt = phitint + dt * (3./2.*dft_n-1./2.*dft_o+0*dft_oo)
      dwdt = phizint + dt * (3./2.*dfz_n-1./2.*dfz_o+0*dfz_oo)
      
      call bound(dudt,dvdt,dwdt,rank)
      call bounds(enew,cnew,hwallh,hwallc,Twallh,Twallc,nrank)
      call stateIG()

      dfr_oo = dfr_o
      dft_oo = dft_o
      dfz_oo = dfz_o
      dfc_oo = dfc_o
      dfr_o = dfr_n!0.5*(dfr_n+dfr_n1)
      dft_o = dft_n!0.5*(dft_n+dft_n1)
      dfz_o = dfz_n!0.5*(dfz_n+dfz_n1)
      dfc_o = dfc_n!0.5*(dfc_n+dfc_n1)


      end




      subroutine cmpbs(bulk,stress,instress,Cbulk,fluxPipe,fluxAnnulus,Qval,cpipe,meanVisc,cannulus,opthck)
c      implicit none

      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      real um(imax),rm(imax),vm(imax),wm(imax),wm2(imax),wmmm2(imax),bulk1,hi_s(imax),hi_c(0:imax),wmmm(imax)
      real mum(imax),mummm(imax),mum_c(0:imax),lambdam(imax),lambdammm(imax),lambdam_c(0:imax)
      real ur(imax),vr(imax),wr(imax),uw(imax),ener1,instress,meanVisc,meanVisc1,cannulus1
      real CTbulk,cm(imax),cmmm(imax),wcm(imax),wcmmm(imax),hic_s(imax),hic_c(0:imax),fluxAnnulus,fluxPipe
      real Qval,km(imax)
      real cpipe,cannulus,opthck1 

      wm=0
      wm2=0
      uw=0
      wcm=0
      cm=0
      rm=0
      km=0
      mum=0
      lambdam=0
      do k=1,kmax/p_col
       do j=1,jmax/p_row
        do i=1,imax
         rm(i)=rm(i)+rhonew(i,j,k)/(jmax*kmax/(p_col*p_row))
         wm(i)=wm(i)+phiznew(i,j,k)/(jmax*kmax/(p_col*p_row))
         km(i)=km(i)+kappa(i,j,k)/(jmax*kmax/(p_col*p_row))
         wm2(i)=wm2(i)+wnew(i,j,k)/(jmax*kmax/(p_col*p_row))
         mum(i)=mum(i)+munew(i,j,k)/(jmax*kmax/(p_col*p_row))
         lambdam(i)=lambdam(i)+lambdanew(i,j,k)/(jmax*kmax/(p_col*p_row))
         wcm(i)=wcm(i)+phiznew(i,j,k)*cnew(i,j,k)/(jmax*kmax/(p_col*p_row))
         cm(i)=cm(i)+cnew(i,j,k)/(jmax*kmax/(p_col*p_row))
        enddo
       enddo
      enddo
      opthck = 0
      opthck1 = 0
      bulk1 = 0
      meanVisc1 = 0 
      CTbulk = 0
      do i=1,imax
       opthck1 = opthck1 + (Ru(i)-Ru(i-1))/2.*km(i) !Replaced 8 with 2/(Rout**2-Rin**2)
       bulk1 = bulk1 + (Ru(i)-Ru(i-1))/2.*wm(i)/rm(i) !Replaced 8 with 2/(Rout**2-Rin**2)
       meanVisc1 = meanVisc1 + (Ru(i)-Ru(i-1))/2.*mum(i) !Replaced 8 with 2/(Rout**2-Rin**2)
       CTbulk = CTbulk + (Ru(i)-Ru(i-1))/2.*wcm(i)/rm(i)
      enddo
      CTbulk = CTbulk/(bulk1)
      bulk1=bulk1!*CTbulk
      call mpi_allreduce(wm,wmmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(wm2,wmmm2  ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(wcm,wcmmm  ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(cm,cmmm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(mum,mummm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(lambdam,lambdammm    ,imax,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      wm = wmmm/px
      wm2 = wmmm2/px
      mum  = mummm/px
      lambdam  = lambdammm/px
      wcm = wcmmm/px
      cm = cmmm/px 
      hi_s = wm2
      hic_s = cm
      call inter_s_c_m(imax,1,mum,mum_c,dr)
      call inter_s_c_m(imax,1,lambdam,lambdam_c,dr)
      call der1w_s_c6_m(imax,1,hi_s,hi_c,dr)
      hi_c = hi_c*mr_c
      stress1 = -hi_c(imax)*mum_c(imax) !c is staggered grid
      stress2 = -hi_c(0)*mum_c(0)       !s is central grid
      call der1c_s_c6_m(imax,1,hic_s,hic_c,dr,Twallh,Twallc) 
      hic_c = hic_c*mr_c
      fluxPipe1 = hic_c(imax)*lambdam_c(imax) 
      fluxAnnulus1 = hic_c(0)*lambdam_c(0)
       
      call inter_s_c(imax,hic_s,hic_c,dr)
      
      cpipe1=hic_c(imax)
      cannulus1 = hic_c(0) 
      call mpi_allreduce(opthck1  ,opthck  ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(bulk1  ,bulk      ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(meanVisc1  ,meanVisc      ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(stress1,stress    ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(stress2,instress  ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(CTbulk ,Cbulk     ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(fluxPipe1,fluxPipe  ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(fluxAnnulus1,fluxAnnulus  ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(cpipe1,cpipe  ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      call mpi_allreduce(cannulus1,cannulus  ,1,MPI_REAL8,MPI_SUM,MPI_COMM_WORLD,ier)
      opthck = opthck/ px
      bulk = bulk/ px
      stress = stress /px
      instress = instress /px 
      Cbulk = Cbulk/px
      fluxPipe = fluxPipe /px
      fluxAnnulus = fluxAnnulus /px
      cpipe = cpipe /px
      cannulus = cannulus /px
      meanVisc = meanVisc /px
      end

      subroutine chkdt()
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      real tmp1,tmp2,dt1,dr1,umax,dtmax
      integer ier
      dtmax = 1.5e-3
      umax =0
      do k=1,kmax/p_col
        do j=1,jmax/p_row
	  do i=1,imax 
      dr1 = Ru(i)-Ru(i-1)
      dt1 = dtheta!*rp(i)
      
!      tmp1 = (abs(unew(i,j,k))/dr1 + 3*abs(vnew(i,j,k))/dt1 + 3*abs(wnew(i,j,k))/dz)
!      tmp2 =10*(munew(i,j,k)/rhonew(i,j,k))*(1./dr1**2+1./dz**2)
      tmp1 = (abs(unew(i,j,k))/dr1 + abs(vnew(i,j,k))/dt1 + abs(wnew(i,j,k))/dz)
      dt1 =0.2/(tmp1)
      umax = max (unew(i,j,k),vnew(i,j,k),wnew(i,j,k),umax)
      enddo
      enddo
      enddo
      call mpi_allreduce(dt1,dt    ,1,MPI_REAL8,MPI_MIN,MPI_COMM_WORLD,ier)
      dt  = min ( dt , dtmax)
c      dt = dtmax
      end
 


      subroutine report
      include 'param.txt'
      include 'common.txt'
      write(6,*) '**************************************************************'
      write(6,*) '**************************************************************'
      write(6,*) '***  DNS of turbulent pipe flow                            ***'
      write(6,111) imax,jmax,kmax
 111  format(     '***  resolution (Nr,Nt,Nz) = ' , 3i7    , '            ***')
      write(6,112) px
 112  format(     '***  Number of processors  = ' , i7              , '         ***')
      write(6,*)
      write(6,*) '**************************************************************'
      write(6,*) '***** number of points  = ', imax*jmax*kmax
      end

      function f_colebrook(Re,bulk)
      real cole
      f=0.0001
      do i=1,20
         Reb = Re*bulk
	f = f - cole(f,Reb)/((cole(f+1e-5,Reb)-cole(f-1e-5,Reb))/2e-5)
      enddo
      f_colebrook=f
      end
      function cole(f,Reb)
       cole = ( 1/sqrt(f) +2.*log10(2.51/(Reb*sqrt(f))))
      end
      function f_shokling(Re,bulk)
      real shok 
      f=0.0001
      do i=1,20
         Reb = Re*bulk
	f = f - shok(f,Reb)/((cole(f+1e-5,Reb)-cole(f-1e-5,Reb))/2e-5)
      enddo
      f_shokling=f
      end
      function shok(f,Reb)
       shok= ( 1/sqrt(f) -1.93*log10(Reb*sqrt(f))+0.537)
      end


      subroutine trunc(f,rank)
      use decomp_2d
      include 'param.txt'
      real wt(2*jmax+15),wk(2*kmax+15)
      real f   (0:i1          ,jmax/p_row,kmax/p_col),fz(kmax),ft(jmax),dk(kmax),dj(jmax)
      real t_f (0:imx,jmax    ,kmax/p_col)
      real t_fz(0:imx,jmax/p_col,kmax    )
      integer mask(imax),ifil(0:imax+1),i_dex
      do i=1,imax
         ifil(i)= jmax!16+12*(i-1)!7+12(i-1)
      enddo
      idex = i_dex(nrank)
      ifil(0)=ifil(1)
      ifil(imax+1)=ifil(imax)
      do i=1,imax
         if (ifil(i).gt. jmax) ifil(i)=jmax
      enddo
      call vrffti(jmax,wt)
      call vrffti(kmax,wk)
      call transpose_x_to_y(f,t_f)
      do i=0,imx
        do k=1,kmax/p_col
         do j=1,jmax
           ft(j)=t_f(i,j,k)
         enddo
         call vrfftf(1,jmax,ft,dj,1,wt)
         do j=ifil(i+idex),jmax
         ft(j)=0
         enddo
         call vrfftb(1,jmax,ft,dj,1,wt)
          do j=1,jmax
          t_f(i,j,k)=ft(j)
          enddo
        enddo
       enddo
      call transpose_y_to_z(t_f,t_fz)
      do i=0,imx
	  do j=1,jmax/p_col
	  do k=1,kmax
	    fz(k)=t_fz(i,j,k)
          enddo
          call vrfftf(1,kmax,fz,dk,1,wk)
          fz(kmax)=0
          call vrfftb(1,kmax,fz,dk,1,wk)
          do k=1,kmax
             t_fz(i,j,k)=fz(k)
          enddo
        enddo
       enddo
       call transpose_z_to_y(t_fz,t_f)
       call transpose_y_to_x(t_f ,  f)
       end

	   

       function i_dex (nrank)
       include 'param.txt'
       i_dex = (nrank/p_col) *(imx+1)
       end  

       function j_dex (nrank)
       include 'param.txt'
       stop
       j_dex = (nrank/p_col) * (jmax/p_row)
       end  

       function k_dex (nrank)
       include 'param.txt'
       stop
       k_dex = nrank/p_row * (kmax/p_col)
       end  
       
      subroutine wait(n,x)
      x =1000.
      do i=1,n*1000
         x = sin(x)+cos(x)
      enddo
      end

      subroutine readTable(rank)
      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'
      integer ierr

c        if (rank.eq.0) then
          open(31,file='SCO2_8MPa_300_309b.txt')
            do i=1,nTab
              READ (31,*) enthTab(i),tempTab(i),rhoTab(i),muTab(i),lamTab(i),cpTab(i)
            enddo
          close(31)
c        endif

c      call MPI_BCAST(enthTab,   nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
c      call MPI_BCAST(tempTab,   nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
c      call MPI_BCAST(rhoTab,    nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
c      call MPI_BCAST(muTab,     nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)
c      call MPI_BCAST(lamTab,    nTab,mpi_real8,0,MPI_COMM_WORLD,ierr)

      call spline(enthTab, tempTab, nTab, temp2Tab)
      call spline(enthTab, rhoTab,  nTab, rho2Tab)
      call spline(enthTab, muTab,   nTab, mu2Tab)
      call spline(enthTab, lamTab,  nTab, lam2Tab)
      call spline(enthTab, cpTab,   nTab, cp2Tab)

      end


      subroutine stateIG()
      implicit none
      include 'param.txt'
      include 'common.txt'

       do k=1,kmax/p_col
          do j=1,jmax/p_row
            do i=0,imax+1
              cnew(i,j,k)=enew(i,j,k)!+1.
              rhonew(i,j,k)=Tplus/(Tplus+cnew(i,j,k))
              cpnew(i,j,k)=1.
              munew(i,j,k)  = 1./Re
              lambdanew(i,j,k) = 1./(Re*Pr)
            enddo
          enddo
        enddo
      return
      end

      subroutine divTrace(rank)
      use decomp_2d
c    calculates divergence
c
c
c       1 dru     1 d v    d w 
c       - --- +   - ---  + --- = 0
c       r  dr     r dt     dz
c
c term    (1)      (2)     (3)
c

      implicit none
      include 'param.txt'
      include 'common.txt'
      include 'mpif.h'

      integer im,ier,idex,i_dex,ierr
      real ft(jmax),dft(jmax),fz(kmax),dfz(kmax),som,dmaxx!,div(0:i1,jmax/p_row,kmax/p_col)
      real t_dmaxx
      real t_som,hi_c(0:imax),hi_s(imax)

      real tmpx(0:i1 ,jmax/p_row,kmax/p_col)         
      real tmpy(0:imx,jmax      ,kmax/p_col)         
      real tmpz(0:imx,jmax/p_col,kmax      )    
!      idex = i_dex(nrank)

      trace = 0

      call transpose_x_to_y(wnew,tmpy)
      call transpose_y_to_z(tmpy,tmpz)
      do i=0,imx
	    do j=1,jmax/p_col
	   do k=1,kmax
	    fz(k) =tmpz(i,j,k)
	   enddo
	   call four1(kmax,fz,dfz,dz)   !dw/dz
	   do k=1,kmax
	    tmpz(i,j,k)=dfz(k)
	   enddo
	  enddo
      enddo
      call transpose_z_to_y(tmpz,tmpy)
      call transpose_y_to_x(tmpy,trace)
     
      call transpose_x_to_y(vnew,tmpy)
      do i=0,imx
	    do k=1,kmax/p_col
	  do j=1,jmax
	    ft(j)=tmpy(i,j,k)
	  enddo
	  call four1(jmax,ft,dft,dtheta)   !dv/dy
	  do j=1,jmax
	   tmpy(i,j,k)=dft(j)!/rp(i+idex)
	  enddo
	    enddo
      enddo
      call transpose_y_to_x(tmpy,tmpx)
      trace = trace  + tmpx
	   
         do j=1,jmax/p_row
	     do k=1,kmax/p_col
	    do i=0,imax
            hi_c(i)=unew(i,j,k)!*ru(i)
         enddo
         call deriv_c_s(imax,hi_c,hi_s,dr)  !du/dx
         hi_s = hi_s*mr_s                   !du/dr
          do i=1,imax
!	  if (ipois.ne.2) trace(i,j,k)=trace(i,j,k)+hi_s(i) !/rp(i)
!	  if (ipois.eq.2) trace(i,j,k)=trace(i,j,k)+(unew(i,j,k)-unew(i-1,j,k))/((ru(i)-ru(i-1)))!+(ru(i)*unew(i,j,k)-ru(i-1)*unew(i-1,j,k))/(rp(i)*(ru(i)-ru(i-1)))
	  trace(i,j,k)=trace(i,j,k)+hi_s(i)
	    enddo
	    enddo
       enddo
       
      trace = 2.0/3.0*munew*trace


      call transpose_x_to_y(trace,tmpy)
      call transpose_y_to_z(tmpy,tmpz)
      do i=0,imx
	     do j=1,jmax/p_col
	   do k=1,kmax
	    fz(k) =tmpz(i,j,k)
	   enddo
	   call four1(kmax,fz,dfz,dz)       !d(2/3mu.grad)/dz
	   do k=1,kmax
	    tmpz(i,j,k)=dfz(k)
	   enddo
	     enddo
      enddo
      call transpose_z_to_y(tmpz,tmpy)
      call transpose_y_to_x(tmpy,tmpx)
      divdivW = tmpx
     
      call transpose_x_to_y(trace,tmpy)
      do i=0,imx
	    do k=1,kmax/p_col
	  do j=1,jmax
	    ft(j)=tmpy(i,j,k)
	  enddo
	  call four1(jmax,ft,dft,dtheta)   !d(2/3mu.grad)/dy
	  do j=1,jmax
	   tmpy(i,j,k)=dft(j)!/rp(i+idex)
	  enddo
	    enddo
      enddo
      call transpose_y_to_x(tmpy,tmpx)
      divdivV = tmpx      

	   do j=1,jmax/p_row
	     do k=1,kmax/p_col
	     do i=1,imax
            hi_s(i)=trace(i,j,k)
         enddo
         call deriv_s_c(imax,hi_s,hi_c,dr)  !d(2/3mu.grad)/dx
         hi_c = hi_c*mr_c                   !d(2/3mu.grad)/dr
          do i=0,imax
!	  if (ipois.eq.2)   divdivU(i,j,k)=(trace(i,j,k)-trace(i-1,j,k))/((ru(i)-ru(i-1)))   !not in U position
         divdivU(i,j,k)=hi_c(i)
	    enddo
	    enddo
       enddo
       
       end
     

c********************************************************************
c     spline (numerical recipes)
c********************************************************************
      subroutine spline(x, y, n, y2)
!   use nrtype
!
! Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e.
! y(i)=f(x(i)), with x(1)<x(2)<...<x(n), and given values yp1 and ypn for
! the first derivative of the interpolating function at points 1 and n,
! respectively, this routine returns an array y2(1:n) of length n which
! contains the second derivatives of the interpolating function at the
! tabulated points x(i).  If yp1 and/or ypn are equal to 1.e30 or larger,
! the routine is signaled to set the corresponding boundary condition for a
! natural spline with zero second derivative on that boundary.
! Parameter: nmax is the largest anticipiated value of n
! (adopted from Numerical Recipes in FORTRAN 77)
!
      INTEGER, PARAMETER :: DP=KIND(1.0D0)
      INTEGER:: n
      INTEGER, PARAMETER:: nmax=5000
      REAL(DP):: yp1, ypn, x(n), y(n), y2(n)
      INTEGER:: i, k
      REAL(DP):: p, qn, sig, un, u(nmax)

        y2(1) = 0.
        u(1)  = 0.
        do i=2, n-1
           sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
           p=sig*y2(i-1)+2.
           y2(i)=(sig-1.)/p
           u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/
     &          (x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
        enddo

        qn=0.
        un=0.
        y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)

        do k=n-1, 1, -1
           y2(k)=y2(k)*y2(k+1)+u(k)
        enddo

        return
      end


      subroutine splint(xa,ya,y2a,n,x,y,khi,klo)
      INTEGER n
      REAL x,y,xa(n),y2a(n),ya(n)
      INTEGER k,khi,klo
      REAL a,b,h

      if ((khi.eq.0) .and. (klo.eq.0)) then
        klo=1
        khi=n
1       if (khi-klo.gt.1) then
          k=(khi+klo)/2
          if(xa(k).gt.x)then
            khi=k
          else
            klo=k
          endif
        goto 1
        endif

      endif

      h=xa(khi)-xa(klo)
      if (h.eq.0.) stop 'bad xa input in splint'
      a=(xa(khi)-x)/h
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
      return
      END

