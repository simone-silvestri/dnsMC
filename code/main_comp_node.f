      program main

      use decomp_2d
      implicit none
      include 'mpif.h'
      integer imax,jmax,kmax,p_col,p_row 
      integer grid_num
      parameter(imax = 166, jmax = 168, kmax = 192)
      parameter(p_col= 24 , p_row=  8 )
      parameter(grid_num = 4)
      real  Tcore(0:imax+1,0:jmax/p_row+1,0:kmax/p_col+1)
      real  Qcore(0:imax+1,0:jmax/p_row+1,0:kmax/p_col+1)
      real  Vcore(0:imax+1,0:jmax/p_row+1,0:kmax/p_col+1)
      real  T   (0:imax+1,0:jmax+1,0:kmax+1)
      real  Tbuf(0:imax+1,0:jmax+1,0:kmax+1)
      real  Qbuf(0:imax+1,0:jmax+1,0:kmax+1)
      real  Qr(0:imax+1,0:jmax/p_row+1,0:kmax+1)
      real  Vr(0:imax+1,0:jmax/p_row+1,0:kmax+1)
      real  x(0:imax+1),dx
      real  Qrm(0:imax+1),Vrm(0:imax+1)
      integer i,j,k,ierr,istat(mpi_status_size),result_proc,color
      integer rank1, size1
      integer node_comm

      logical periodic_bc(3)
      character*200 filename
      character*(MPI_MAX_PROCESSOR_NAME) name_proc
  
      call mpi_init(ierr)

      periodic_bc(1) = .false.
      periodic_bc(2) = .true.
      periodic_bc(3) = .true.

      call decomp_2d_init(imax,jmax,kmax,p_row,p_col,periodic_bc)

      call mpi_get_processor_name(name_proc, result_proc,ierr)

      color = xstart(2)
      call mpi_comm_split(MPI_COMM_WORLD, color, nrank, node_comm, ierr)
      call mpi_comm_rank(node_comm, rank1, ierr)
      call mpi_comm_size(node_comm, size1, ierr)

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      do i=0,p_row*p_col-1
       if(nrank.eq.i) then
        write(*,'(A20,4I5)') name_proc, nrank, rank1, xstart(2), xstart(3)
       endif
       call mpi_barrier(MPI_COMM_WORLD,ierr)
      enddo 

      stop

      dx = 1.0/imax
      x(0) = -dx / 2.0;
      do i=1,imax+1
       x(i) = x(i-1) + dx
      enddo
      do i=0,imax+1
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         Tcore(i,j,k) = 500 - 2000*x(i)**2 + 2000*x(i)
        enddo
       enddo
      enddo 
      Tbuf = 0
      do i=xstart(1),xend(1)
       do j=xstart(2),xend(2)
        do k=xstart(3),xend(3)
         Tbuf(i,j,k) = Tcore(i-xstart(1)+1,j-xstart(2)+1,k-xstart(3)+1)
        enddo
       enddo
      enddo

      call mpi_allreduce(Tbuf, T, (imax+2)*(jmax+2)*(kmax+2), MPI_REAL,
     1                    MPI_SUM, MPI_COMM_WORLD, ierr)

      if(rank1.eq.0) then
 
      open(unit=1,file='tempfort.txt')
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax
         write(1,*) x(i),j,k,T(i,j,k)
        enddo
       enddo
      enddo
      close(1)

      open(unit=1,file='resfort.txt')
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax
         write(1,*) x(i),j,k,Qr(i,j,k),Vr(k,j,i)
        enddo
       enddo
      enddo
      close(1)
      Qrm = 0
      Vrm = 0
      do i=1,imax
       do i=j,jmax/p_row
        do k=1,kmax
         Qrm(i) = Qrm(i) + Qr(i,j,k) / (kmax*jmax)
         Vrm(i) = Vrm(i) + Vr(i,j,k) / (kmax*jmax)
        enddo
       enddo
      enddo
      open(unit=1,file='resfortm.txt')
      do i=1,imax
       write(1,*) x(i),Qrm(i),Vrm(i)
      enddo
      close(1)
      endif
        
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      Qbuf = 0
      Qcore= 0
      call mpi_allreduce(Qr, Qbuf, (imax+2)*(jmax/p_row+2)*(kmax+2), MPI_REAL,
     1                    MPI_SUM, node_comm, ierr)
     
      do i=1,imax
       do j=1,jmax/p_row
        do k=xstart(3),xend(3)
         Qcore(i,j,k-xstart(3)+1) =Qbuf(i,j,k)
        enddo
       enddo
      enddo
      write(filename,'(A7,I2.2,A4)') 'rescore',nrank,'.txt'
      open(unit = 1, file = filename)
      do i=1,imax
       do j=1,jmax/p_row
        do k=xstart(3),xend(3)
       write(1,*) i,j,k,Qcore(i,j,k-xstart(3)+1)
        enddo
       enddo
      enddo
      close(1)


      call decomp_2d_finalize
      call mpi_finalize(ierr)

      end
      
      subroutine mkgrid(x,y,z,im,jm,km)
      implicit none
      integer im,jm,km
      real x(0:im+1)
      real y(0:jm+1)
      real z(0:km+1)
      real dx,dy,dz
      integer i,j,k
      
      x(0) = -dx/2.0
      y(0) = -dy/2.0
      z(0) = -dz/2.0
      do i=1,im+1
       x(i) = x(i-1) + dx
      enddo
      do i=1,jm+1
       y(i) = y(i-1) + dy
      enddo
      do i=1,km+1
       z(i) = z(i-1) + dz
      enddo

      end subroutine

