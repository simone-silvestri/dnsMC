
      module gpu_module
        INTERFACE
          subroutine MC_GPU() BIND (C, NAME='mc_gpu')
            USE ISO_C_BINDING
            implicit none
          end subroutine MC_GPU
        END INTERFACE
      end module gpu_module

      program main

      use decomp_2d
      include 'mpif.h'
      integer imax,jmax,kmax,p_col,p_row 
      integer grid_num
      real Lx,Ly
      parameter(imax = 64, jmax = 64, kmax = 64)
      parameter(p_col= 1 , p_row= 1 ) 
      parameter(Lx= 1.0, Ly= 1.0)
      real, parameter :: pi = 3.1415927
      real  Tcore(0:imax+1,0:jmax/p_row+1,0:kmax/p_col+1)
      real  Tpart(0:imax+1,0:jmax/p_row+1,0:kmax/p_col+1)
      real  Ccore(0:imax+1,0:jmax/p_row+1,0:kmax/p_col+1)
      real  Qcore(0:imax+1,0:jmax/p_row+1,0:kmax/p_col+1)
      real  Vcore(0:imax+1,0:jmax/p_row+1,0:kmax/p_col+1)
      real  T   (0:imax+1,0:jmax+1,0:kmax+1)
      real  Tbuf (0:imax+1,0:jmax+1,0:kmax+1)
      real  Qbuf(0:imax+1,0:jmax/p_row+1,0:kmax+1)
      real  Qr(0:imax+1,0:jmax/p_row+1,0:kmax+1)
      real  Vr(0:imax+1,0:jmax/p_row+1,0:kmax+1)
      real  x(0:imax+1),dx,y(0:jmax+1),dy
      real  Qrm(0:imax+1),Vrm(0:imax+1),Vrm2(0:imax+1)
      integer i,j,k,ierr 
      integer istat(mpi_status_size)
      integer result_proc, color, MPI_COMM_NODE
      integer node_rank, node_size
 
      logical periodic_bc(3)
      character*200 filename
      character*(MPI_MAX_PROCESSOR_NAME) node_name  

      call mpi_init(ierr)

      periodic_bc(1) = .false.
      periodic_bc(2) = .true.
      periodic_bc(3) = .true.

      call decomp_2d_init(imax,jmax,kmax,p_row,p_col,periodic_bc)

      call mpi_get_processor_name(node_name, result_proc, ierr)

      color = xstart(2)
      call mpi_comm_split(MPI_COMM_WORLD, color, nrank, MPI_COMM_NODE, ierr)
      call mpi_comm_rank(MPI_COMM_NODE, node_rank, ierr)
      call mpi_comm_size(MPI_COMM_NODE, node_size, ierr)
 
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      if (nrank.eq.0) then
       write(6,*) xstart(1),xstart(2),xstart(3),xend(1),xend(2),xend(3)
      endif
 
      call mkgrid(x,Lx,imax)
      do i=0,imax+1
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         Tcore(i,j,k) =500 - 2000*(x(i)/Lx)**2 + 2000*(x(i)/Lx) !400 + 1400 * (sin(2 * pi * x(i) /Lx))**2. !1000 + 1504*x(i)**2 - 1704*x(i) !
        enddo
       enddo
      enddo 


      stime3 = MPI_WTIME()
      Tbuf = 0
      do i=xstart(1),xend(1)
       do j=xstart(2),xend(2)
        do k=xstart(3),xend(3)
         Tbuf(i,j,k) = Tcore(i-xstart(1)+1,j-xstart(2)+1,k-xstart(3)+1)
        enddo
       enddo
      enddo

      call mpi_allreduce(Tbuf , T , (imax+2)*(jmax+2)*(kmax+2), MPI_REAL,
     1                    MPI_SUM, MPI_COMM_WORLD, ierr)

      if(node_rank.eq.0) then

      open(unit=1,file='temp.txt')
      do i=0,imax+1
       do j=0,jmax+1
        do k=0,kmax+1
         write(1,*) i,j,k,T(i,j,k)
        enddo
       enddo
      enddo
      close(1)
 
      stime1 = MPI_WTIME()
      call MC_GPU(T, xstart(2)) 

      call GET_RESULTS(Qr, Vr) 
      stime2 = MPI_WTIME()

      endif


       
      call mpi_barrier(MPI_COMM_NODE,ierr)
      Qbuf = 0
      Qcore= 0
      call mpi_allreduce(Qr, Qbuf, (imax+2)*(jmax/p_row+2)*(kmax+2), MPI_REAL,
     1                    MPI_SUM, MPI_COMM_NODE, ierr)
    
      do i=1,imax
       do j=1,jmax/p_row
        do k=xstart(3),xend(3)
         Qcore(i,j,k-xstart(3)+1) = Qbuf(i,j,k)
        enddo
       enddo
      enddo


      stime4 = MPI_WTIME()

      if(node_rank.eq.0) then
      write(filename,'(A7,A6,A4)') 'resfort',node_name,'.txt'
      open(unit=1,file=filename)
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax
         write(1,*) i,j,k,Qr(i,j,k),Vr(k,j,i)
        enddo
       enddo
      enddo
      close(1)
      endif
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if(nrank.eq.0) then
      Qrm = 0
      Vrm = 0
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax
         Qrm(i) = Qrm(i) + Qr(i,j,k) / (kmax*jmax/p_row)
         Vrm(i) = Vrm(i) + Vr(i,j,k) / (kmax*jmax/p_row)
        enddo
       enddo
      enddo
      Vrm2 = 0
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax
         Vrm2(i) = Vrm2(i) + (Qr(i,j,k) - Qrm(i))**2. / (kmax*jmax/p_row)
        enddo
       enddo
      enddo
      open(unit=1,file='resfortm.txt')
      do i=1,imax
        write(1,*) x(i),Qrm(i),Vrm(i),Vrm2(i)**(0.5)
      enddo
      close(1)
      endif
 
      if(node_rank.eq.0) then     
      write(6,*) 'Radiation time step : ', stime2 - stime1 
      write(6,*) 'Total time step : ', stime4 - stime3 
      endif

      if(nrank.lt.100) then
       write(filename,'(A7,I2.2,A4)') 'rescore',nrank,'.txt'
      else
       write(filename,'(A7,I3.3,A4)') 'rescore',nrank,'.txt'
      endif
      open(unit = 1, file = filename)
      do i=1,imax
       do j=xstart(2),xend(2)
        do k=xstart(3),xend(3)
         write(1,*) i,j,k,Qcore(i,j-xstart(2)+1,k-xstart(3)+1)
        enddo
       enddo
      enddo
      close(1)


      call decomp_2d_finalize
      call mpi_finalize(ierr)

      end

      subroutine mkgrid(rp,Lx,imax)
      implicit none
      real rmax,Lx,rp(0:imax+1),ru(0:imax),delta(0:imax),dr,x,dx,rnorm
      integer imax,i
      dr = Lx/imax
      rmax = Lx
      ru(0)=0. 
      do i=1,imax/2
          x  = 1.*i/imax
          dx = 0.5-1.45*(x-0.5)**2.
          ru(i)=ru(i-1)+dx
        enddo
        rnorm = ru(imax/2)
        do i=1,imax/2
          ru(i)=Lx/2.*ru(i)/rnorm
        enddo
       do i=imax,imax/2+1,-1
         ru(i)=Lx-ru(imax-i)
       enddo

      do i=1,imax
        rp(i)=0.5*(ru(i)+ru(i-1))
        delta(i)=ru(i)-ru(i-1)
      enddo

      rp(0) =-rp(1)
      rp(imax+1)=ru(imax)+(ru(imax)-rp(imax))

      end 


