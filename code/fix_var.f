
      subroutine read_kP(kP,Tnb,nrank)
      implicit none
      include "param.txt"
      include "common.txt"
  
      real kP(nTemp),Tnb(nTemp)
      integer nrank

      open(unit=1,file="../distr/tables/planck-mean.txt")
      do i=1,nTemp       
        read(1,*) Tnb(i),kP(i)
      enddo
      close(1)
      
      end

      subroutine fix_kappa(kP,Tnb)
      implicit none
      include "param.txt"
      include "common.txt"
  
      real kP(nTemp),Tnb(nTemp)
      real ctemp(0:imax+1,jmax/p_row,kmax/p_col)
      
      ctemp = cnew * Tc/Tplus + Tc
      do i=1,imax
       do j=1,jmax/p_row
        do k=1,kmax/p_col
         call linear_int(kP,Tnb,ctemp(i,j,k),kappa(i,j,k),nTemp)
        enddo
       enddo
      enddo

      end

      subroutine linear_int(y,x,x0,y0,n)
      implicit none
      include "param.txt"
      real y(n),x(n),x0,y0
      integer t,n
 
      t = int(x0 - x(1)) / int(x(2) - x(1)) + 1

      y0 = (y(t+1) - y(t)) / (x(t+1) - x(t)) * (x0 - x(t)) + y(t)
 
      end

      subroutine fix_temp(cnew, TMC, Tc, Tplus, startc, endc)
      implicit none
      include "mpif.h"
      include "param.txt"

      real cnew(0:imax+1,jmax/p_row,kmax/p_col)
      real(kind=4) :: TMC(0:imax+1,0:jmax+1,0:kmax+1)
      real Tbuf (0:imax+1,0:jmax+1,0:kmax+1)
      real Tbuf2(0:imax+1,0:jmax+1,0:kmax+1)
      real Tc, Tplus

      integer istat(mpi_status_size), startc(3), endc(3),ierr

      Tbuf = 0
      do i=startc(1),endc(1)
       do j=startc(2),endc(2)
        do k=startc(3),endc(3)
         Tbuf(i,j,k) = cnew(i-startc(1)+1,j-startc(2)+1,k-startc(3)+1)*Tc/Tplus+ Tc
        enddo
       enddo
      enddo

      call mpi_allreduce(Tbuf, Tbuf2, (imax+2)*(jmax+2)*(kmax+2),MPI_REAL8,
     1                    MPI_SUM, MPI_COMM_WORLD, ierr)

      TMC = Tbuf2      

      end

      subroutine fix_rad(divQ, vcore, qr, vr, Tc, MPI_COMM_NODE, startc, endc)
      implicit none
      include "mpif.h"
      include "param.txt"

      real divQ (0:imax+1,jmax/p_row,kmax/p_col)
      real vcore(0:imax+1,jmax/p_row,kmax/p_col)
      real(kind=4) :: qr  (0:imax+1,0:jmax/p_row+1,0:kmax+1)
      real(kind=4) :: vr  (0:imax+1,0:jmax/p_row+1,0:kmax+1)
      real(kind=4) :: Qbuf(0:imax+1,0:jmax/p_row+1,0:kmax+1)
      real Tc
 
      integer istat(mpi_status_size), startc(3), endc(3),ierr
      integer MPI_COMM_NODE 
   
      Qbuf = 0
      divQ= 0
      call mpi_allreduce(qr, Qbuf, (imax+2)*(jmax/p_row+2)*(kmax+2), MPI_REAL,
     1                    MPI_SUM, MPI_COMM_NODE, ierr)

      do i=1,imax
       do j=1,jmax/p_row
        do k=startc(3),endc(3)
         divQ(i,j,k-startc(3)+1) = Qbuf(i,j,k)/stefan/(Tc**4.)
        enddo
       enddo
      enddo
   
      call mpi_barrier(MPI_COMM_NODE,ierr)

      Qbuf = 0
      vcore= 0
      call mpi_allreduce(vr, Qbuf, (imax+2)*(jmax/p_row+2)*(kmax+2),MPI_REAL,
     1                    MPI_SUM, MPI_COMM_NODE, ierr)

      do i=1,imax
       do j=1,jmax/p_row
        do k=startc(3),endc(3)
         vcore(i,j,k-startc(3)+1) = Qbuf(i,j,k)/stefan/(Tc**4.)
        enddo
       enddo
      enddo

      end
