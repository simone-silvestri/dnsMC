      subroutine four1(n,a,da,dx)
      implicit none
      integer n,i, j,jwav
      real a(n),da(n),dx,h,rn
      real wk(2*n+15)
      call vrffti(n,wk)
c
c    Goto Fourier space
c
      do i=1,n
       da(i)=a(i)
      enddo
      call vrfftf(1,n,da,a,1,wk)
      da(1)=0.
      da(n)=0.
      do j=2,n-1,2
          jwav=j/2
          h =da(j+1)
          da(j+1)=da(j)*jwav
          da(j  )= - h*jwav
      enddo
      call vrfftb(1,n,da,a,1,wk)
      rn = 8*atan(1.)/(n*dx)
      do i=1,n
	da(i)=da(i)*rn
      enddo
      end



      subroutine four12(n,a,da,dda,dx)
      implicit none
      integer n,i, j,jwav
      real a(n),da(n),dx,h,rn,dda(n)
      real wk(2*n+15)
      call vrffti(n,wk)
c
c    Goto Fourier space
c
      do i=1,n
       da(i)=a(i)
      enddo
      call vrfftf(1,n,da,a,1,wk)
      dda = da
      da(1)=0.
      da(n)=0.
      do j=2,n-1,2
          jwav=j/2
          h =da(j+1)
          da(j+1)=da(j)*jwav
          da(j  )= - h*jwav
      enddo
      do j=1,n
         jwav = j/2
	 dda(j)=-dda(j)*(jwav)**2
      enddo
      call vrfftb(1,n,da,a,1,wk)
      call vrfftb(1,n,dda,a,1,wk)
      rn = 8*atan(1.)/(n*dx)
      do i=1,n
	da(i)=da(i)*rn
	dda(i)=dda(i)*rn*rn
      enddo
      end
      subroutine four1r(n,itr,a,da,dx)
      implicit none
      integer n,i, j,jwav,itr
      real a(n),da(n),dx,h,rn
      real wk(2*n+15)
      call vrffti(n,wk)

c
c    Goto Fourier space
c
      do i=1,n
       da(i)=a(i)
      enddo
      call vrfftf(1,n,da,a,1,wk)
      da(1)=0.
      da(n)=0.
      do j=2,n-1,2
          jwav=j/2
          h =da(j+1)
          da(j+1)=da(j)*jwav
          da(j  )= - h*jwav
      enddo

      do j=itr,n
      da(j)=0
      enddo

      call vrfftb(1,n,da,a,1,wk)
      rn = 8*atan(1.)/(n*dx)
      do i=1,n
	da(i)=da(i)*rn
      enddo
      end

      subroutine four2(n,a,dda,dx)
      integer n
      real a(n),dda(n)
      real wk(2*n+15)
      call vrffti(n,wk)

c
c    Goto Fourier space
c
      do i=1,n
	 dda(i)=a(i)
      enddo
      call vrfftf(1,n,dda,a,1,wk)

      do j=1,n
         jwav = j/2  
         dda(j) =- dda(j)*jwav**2
      enddo

c
c    Back to physical space
c
      call vrfftb(1,n,dda,a,1,wk)
      rn = (8*atan(1.)/(n*dx))**2

      do j=1,n
	   dda(j)=dda(j)*rn
      enddo
      end

      subroutine four2r(n,itr,a,dda,dx)
      integer n,itr
      real a(n),dda(n)
      real wk(2*n+15)
      call vrffti(n,wk)

c
c    Goto Fourier space
c
      dda =  0
      do i=1,n
	 dda(i)=a(i)
      enddo
      call vrfftf(1,n,dda,a,1,wk)

      do j=1,n
         jwav = j/2  
         dda(j) =- dda(j)*jwav**2
      enddo
      do j=itr,n
       dda(j)=0
      enddo

c
c    Back to physical space
c
      call vrfftb(1,n,dda,a,1,wk)
c    
      rn = (8*atan(1.)/(n*dx))**2

      do j=1,n
	   dda(j)=dda(j)*rn
      enddo
c
c   all done
c 
      end

       function pimach()
       pimach = 4.*atan(1.)
       end


      SUBROUTINE VRADB2 (MP,IDO,L1,CC,CH,MDIMC,WA1)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION  CC(MDIMC,IDO,2,L1)    ,CH(MDIMC,IDO,L1,2),
     1                WA1(IDO)
      DO 101 K=1,L1
          DO 1001 M=1,MP
         CH(M,1,K,1) = CC(M,1,1,K)+CC(M,IDO,2,K)
         CH(M,1,K,2) = CC(M,1,1,K)-CC(M,IDO,2,K)
 1001     CONTINUE
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
               DO 1002 M=1,MP
            CH(M,I-1,K,1) = CC(M,I-1,1,K)+CC(M,IC-1,2,K)
            CH(M,I,K,1) = CC(M,I,1,K)-CC(M,IC,2,K)
            CH(M,I-1,K,2) = WA1(I-2)*(CC(M,I-1,1,K)-CC(M,IC-1,2,K))
     1      -WA1(I-1)*(CC(M,I,1,K)+CC(M,IC,2,K))
            CH(M,I,K,2) = WA1(I-2)*(CC(M,I,1,K)+CC(M,IC,2,K))+WA1(I-1)
     1      *(CC(M,I-1,1,K)-CC(M,IC-1,2,K))
 1002          CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
          DO 1003 M=1,MP
         CH(M,IDO,K,1) = CC(M,IDO,1,K)+CC(M,IDO,1,K)
         CH(M,IDO,K,2) = -(CC(M,1,2,K)+CC(M,1,2,K))
 1003     CONTINUE
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE VRADB3 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION  CC(MDIMC,IDO,3,L1)    ,CH(MDIMC,IDO,L1,3),
     1                WA1(IDO)   ,WA2(IDO)
      ARG=2.*PIMACH(1.0)/3.
      TAUR=COS(ARG)
      TAUI=SIN(ARG)
      DO 101 K=1,L1
          DO 1001 M=1,MP
         CH(M,1,K,1) = CC(M,1,1,K)+2.*CC(M,IDO,2,K)
         CH(M,1,K,2) = CC(M,1,1,K)+(2.*TAUR)*CC(M,IDO,2,K)
     1   -(2.*TAUI)*CC(M,1,3,K)
         CH(M,1,K,3) = CC(M,1,1,K)+(2.*TAUR)*CC(M,IDO,2,K)
     1   +2.*TAUI*CC(M,1,3,K)
 1001     CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
               DO 1002 M=1,MP
            CH(M,I-1,K,1) = CC(M,I-1,1,K)+(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
            CH(M,I,K,1) = CC(M,I,1,K)+(CC(M,I,3,K)-CC(M,IC,2,K))
            CH(M,I-1,K,2) = WA1(I-2)*
     1 ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))-
     * (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K))))
     2                   -WA1(I-1)*
     3 ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))+
     * (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))))
            CH(M,I,K,2) = WA1(I-2)*
     4 ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))+
     8 (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))))
     5                  +WA1(I-1)*
     6 ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))-
     8 (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K))))
              CH(M,I-1,K,3) = WA2(I-2)*
     7 ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))+
     8 (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K))))
     8                      -WA2(I-1)*
     9 ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))-
     8 (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))))
            CH(M,I,K,3) = WA2(I-2)*
     1 ((CC(M,I,1,K)+TAUR*(CC(M,I,3,K)-CC(M,IC,2,K)))-
     8 (TAUI*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))))
     2                 +WA2(I-1)*
     3 ((CC(M,I-1,1,K)+TAUR*(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))+
     8 (TAUI*(CC(M,I,3,K)+CC(M,IC,2,K))))
 1002          CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
      SUBROUTINE VRADB4 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2,WA3)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION  CC(MDIMC,IDO,4,L1)  ,CH(MDIMC,IDO,L1,4)    ,
     1                WA1(IDO)  ,WA2(IDO)  ,WA3(IDO)
      SQRT2=SQRT(2.)
      DO 101 K=1,L1
          DO 1001 M=1,MP
         CH(M,1,K,3) = (CC(M,1,1,K)+CC(M,IDO,4,K))
     1   -(CC(M,IDO,2,K)+CC(M,IDO,2,K))
         CH(M,1,K,1) = (CC(M,1,1,K)+CC(M,IDO,4,K))
     1   +(CC(M,IDO,2,K)+CC(M,IDO,2,K))
         CH(M,1,K,4) = (CC(M,1,1,K)-CC(M,IDO,4,K))
     1   +(CC(M,1,3,K)+CC(M,1,3,K))
         CH(M,1,K,2) = (CC(M,1,1,K)-CC(M,IDO,4,K))
     1   -(CC(M,1,3,K)+CC(M,1,3,K))
 1001     CONTINUE
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
               DO 1002 M=1,MP
            CH(M,I-1,K,1) = (CC(M,I-1,1,K)+CC(M,IC-1,4,K))
     1      +(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
            CH(M,I,K,1) = (CC(M,I,1,K)-CC(M,IC,4,K))
     1      +(CC(M,I,3,K)-CC(M,IC,2,K))
            CH(M,I-1,K,2)=WA1(I-2)*((CC(M,I-1,1,K)-CC(M,IC-1,4,K))
     1      -(CC(M,I,3,K)+CC(M,IC,2,K)))-WA1(I-1)
     1      *((CC(M,I,1,K)+CC(M,IC,4,K))+(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))
            CH(M,I,K,2)=WA1(I-2)*((CC(M,I,1,K)+CC(M,IC,4,K))
     1      +(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))+WA1(I-1)
     1      *((CC(M,I-1,1,K)-CC(M,IC-1,4,K))-(CC(M,I,3,K)+CC(M,IC,2,K)))
            CH(M,I-1,K,3)=WA2(I-2)*((CC(M,I-1,1,K)+CC(M,IC-1,4,K))
     1      -(CC(M,I-1,3,K)+CC(M,IC-1,2,K)))-WA2(I-1)
     1      *((CC(M,I,1,K)-CC(M,IC,4,K))-(CC(M,I,3,K)-CC(M,IC,2,K)))
            CH(M,I,K,3)=WA2(I-2)*((CC(M,I,1,K)-CC(M,IC,4,K))
     1      -(CC(M,I,3,K)-CC(M,IC,2,K)))+WA2(I-1)
     1      *((CC(M,I-1,1,K)+CC(M,IC-1,4,K))-(CC(M,I-1,3,K)
     1      +CC(M,IC-1,2,K)))
            CH(M,I-1,K,4)=WA3(I-2)*((CC(M,I-1,1,K)-CC(M,IC-1,4,K))
     1      +(CC(M,I,3,K)+CC(M,IC,2,K)))-WA3(I-1)
     1     *((CC(M,I,1,K)+CC(M,IC,4,K))-(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))
            CH(M,I,K,4)=WA3(I-2)*((CC(M,I,1,K)+CC(M,IC,4,K))
     1      -(CC(M,I-1,3,K)-CC(M,IC-1,2,K)))+WA3(I-1)
     1      *((CC(M,I-1,1,K)-CC(M,IC-1,4,K))+(CC(M,I,3,K)+CC(M,IC,2,K)))
 1002          CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
               DO 1003 M=1,MP
         CH(M,IDO,K,1) = (CC(M,IDO,1,K)+CC(M,IDO,3,K))
     1   +(CC(M,IDO,1,K)+CC(M,IDO,3,K))
         CH(M,IDO,K,2) = SQRT2*((CC(M,IDO,1,K)-CC(M,IDO,3,K))
     1   -(CC(M,1,2,K)+CC(M,1,4,K)))
         CH(M,IDO,K,3) = (CC(M,1,4,K)-CC(M,1,2,K))
     1   +(CC(M,1,4,K)-CC(M,1,2,K))
         CH(M,IDO,K,4) = -SQRT2*((CC(M,IDO,1,K)-CC(M,IDO,3,K))
     1   +(CC(M,1,2,K)+CC(M,1,4,K)))
 1003          CONTINUE
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE VRADB5 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2,WA3,WA4)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION  CC(MDIMC,IDO,5,L1)    ,CH(MDIMC,IDO,L1,5),
     1             WA1(IDO)     ,WA2(IDO)     ,WA3(IDO)     ,WA4(IDO)
      ARG=2.*PIMACH(1.0)/5.
      TR11=COS(ARG)
      TI11=SIN(ARG)
      TR12=COS(2.*ARG)
      TI12=SIN(2.*ARG)
      DO 101 K=1,L1
      DO 1001 M=1,MP
         CH(M,1,K,1) = CC(M,1,1,K)+2.*CC(M,IDO,2,K)+2.*CC(M,IDO,4,K)
         CH(M,1,K,2) = (CC(M,1,1,K)+TR11*2.*CC(M,IDO,2,K)
     1   +TR12*2.*CC(M,IDO,4,K))-(TI11*2.*CC(M,1,3,K)
     1   +TI12*2.*CC(M,1,5,K))
         CH(M,1,K,3) = (CC(M,1,1,K)+TR12*2.*CC(M,IDO,2,K)
     1   +TR11*2.*CC(M,IDO,4,K))-(TI12*2.*CC(M,1,3,K)
     1   -TI11*2.*CC(M,1,5,K))
         CH(M,1,K,4) = (CC(M,1,1,K)+TR12*2.*CC(M,IDO,2,K)
     1   +TR11*2.*CC(M,IDO,4,K))+(TI12*2.*CC(M,1,3,K)
     1   -TI11*2.*CC(M,1,5,K))
         CH(M,1,K,5) = (CC(M,1,1,K)+TR11*2.*CC(M,IDO,2,K)
     1   +TR12*2.*CC(M,IDO,4,K))+(TI11*2.*CC(M,1,3,K)
     1   +TI12*2.*CC(M,1,5,K))
 1001          CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
      DO 1002 M=1,MP
            CH(M,I-1,K,1) = CC(M,I-1,1,K)+(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +(CC(M,I-1,5,K)+CC(M,IC-1,4,K))
            CH(M,I,K,1) = CC(M,I,1,K)+(CC(M,I,3,K)-CC(M,IC,2,K))
     1      +(CC(M,I,5,K)-CC(M,IC,4,K))
            CH(M,I-1,K,2) = WA1(I-2)*((CC(M,I-1,1,K)+TR11*
     1      (CC(M,I-1,3,K)+CC(M,IC-1,2,K))+TR12
     1      *(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))-(TI11*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))+TI12*(CC(M,I,5,K)+CC(M,IC,4,K))))
     1      -WA1(I-1)*((CC(M,I,1,K)+TR11*(CC(M,I,3,K)-CC(M,IC,2,K))
     1      +TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))+(TI11*(CC(M,I-1,3,K)
     1      -CC(M,IC-1,2,K))+TI12*(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
            CH(M,I,K,2) = WA1(I-2)*((CC(M,I,1,K)+TR11*(CC(M,I,3,K)
     1      -CC(M,IC,2,K))+TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      +(TI11*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))+TI12
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))+WA1(I-1)
     1      *((CC(M,I-1,1,K)+TR11*(CC(M,I-1,3,K)
     1      +CC(M,IC-1,2,K))+TR12*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))
     1      -(TI11*(CC(M,I,3,K)+CC(M,IC,2,K))+TI12
     1      *(CC(M,I,5,K)+CC(M,IC,4,K))))
            CH(M,I-1,K,3) = WA2(I-2)
     1      *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))-(TI12*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K))))
     1     -WA2(I-1)
     1     *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-
     1      CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      +(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
            CH(M,I,K,3) = WA2(I-2)
     1     *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-
     1      CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      +(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
     1      +WA2(I-1)
     1      *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))-(TI12*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K))))
            CH(M,I-1,K,4) = WA3(I-2)
     1      *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI12*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K))))
     1      -WA3(I-1)
     1     *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-
     1      CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      -(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
            CH(M,I,K,4) = WA3(I-2)
     1     *((CC(M,I,1,K)+TR12*(CC(M,I,3,K)-
     1      CC(M,IC,2,K))+TR11*(CC(M,I,5,K)-CC(M,IC,4,K)))
     1      -(TI12*(CC(M,I-1,3,K)-CC(M,IC-1,2,K))-TI11
     1      *(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
     1      +WA3(I-1)
     1      *((CC(M,I-1,1,K)+TR12*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR11*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI12*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))-TI11*(CC(M,I,5,K)+CC(M,IC,4,K))))
            CH(M,I-1,K,5) = WA4(I-2)
     1      *((CC(M,I-1,1,K)+TR11*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR12*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI11*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))+TI12*(CC(M,I,5,K)+CC(M,IC,4,K))))
     1      -WA4(I-1)
     1      *((CC(M,I,1,K)+TR11*(CC(M,I,3,K)-CC(M,IC,2,K))
     1      +TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))-(TI11*(CC(M,I-1,3,K)
     1      -CC(M,IC-1,2,K))+TI12*(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
            CH(M,I,K,5) = WA4(I-2)
     1      *((CC(M,I,1,K)+TR11*(CC(M,I,3,K)-CC(M,IC,2,K))
     1      +TR12*(CC(M,I,5,K)-CC(M,IC,4,K)))-(TI11*(CC(M,I-1,3,K)
     1      -CC(M,IC-1,2,K))+TI12*(CC(M,I-1,5,K)-CC(M,IC-1,4,K))))
     1      +WA4(I-1)
     1      *((CC(M,I-1,1,K)+TR11*(CC(M,I-1,3,K)+CC(M,IC-1,2,K))
     1      +TR12*(CC(M,I-1,5,K)+CC(M,IC-1,4,K)))+(TI11*(CC(M,I,3,K)
     1      +CC(M,IC,2,K))+TI12*(CC(M,I,5,K)+CC(M,IC,4,K))))
 1002          CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
      SUBROUTINE VRADBG (MP,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
     *                 MDIMC,WA)
      DIMENSION    CH(MDIMC,IDO,L1,IP)    ,CC(MDIMC,IDO,IP,L1) ,
     1           C1(MDIMC,IDO,L1,IP)     ,C2(MDIMC,IDL1,IP),
     2                CH2(MDIMC,IDL1,IP)       ,WA(IDO)
      TPI=2.*PIMACH(1.0)
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IPP2 = IP+2
      IPPH = (IP+1)/2
      IF (IDO .LT. L1) GO TO 103
      DO 102 K=1,L1
         DO 101 I=1,IDO
            DO 1001 M=1,MP
            CH(M,I,K,1) = CC(M,I,1,K)
 1001       CONTINUE
  101    CONTINUE
  102 CONTINUE
      GO TO 106
  103 DO 105 I=1,IDO
         DO 104 K=1,L1
            DO 1004 M=1,MP
            CH(M,I,K,1) = CC(M,I,1,K)
 1004       CONTINUE
  104    CONTINUE
  105 CONTINUE
  106 DO 108 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 107 K=1,L1
            DO 1007 M=1,MP
            CH(M,1,K,J) = CC(M,IDO,J2-2,K)+CC(M,IDO,J2-2,K)
            CH(M,1,K,JC) = CC(M,1,J2-1,K)+CC(M,1,J2-1,K)
 1007       CONTINUE
  107    CONTINUE
  108 CONTINUE
      IF (IDO .EQ. 1) GO TO 116
      IF (NBD .LT. L1) GO TO 112
      DO 111 J=2,IPPH
         JC = IPP2-J
         DO 110 K=1,L1
            DO 109 I=3,IDO,2
               IC = IDP2-I
               DO 1009 M=1,MP
               CH(M,I-1,K,J) = CC(M,I-1,2*J-1,K)+CC(M,IC-1,2*J-2,K)
               CH(M,I-1,K,JC) = CC(M,I-1,2*J-1,K)-CC(M,IC-1,2*J-2,K)
               CH(M,I,K,J) = CC(M,I,2*J-1,K)-CC(M,IC,2*J-2,K)
               CH(M,I,K,JC) = CC(M,I,2*J-1,K)+CC(M,IC,2*J-2,K)
 1009          CONTINUE
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
      GO TO 116
  112 DO 115 J=2,IPPH
         JC = IPP2-J
         DO 114 I=3,IDO,2
            IC = IDP2-I
            DO 113 K=1,L1
               DO 1013 M=1,MP
               CH(M,I-1,K,J) = CC(M,I-1,2*J-1,K)+CC(M,IC-1,2*J-2,K)
               CH(M,I-1,K,JC) = CC(M,I-1,2*J-1,K)-CC(M,IC-1,2*J-2,K)
               CH(M,I,K,J) = CC(M,I,2*J-1,K)-CC(M,IC,2*J-2,K)
               CH(M,I,K,JC) = CC(M,I,2*J-1,K)+CC(M,IC,2*J-2,K)
 1013          CONTINUE
  113       CONTINUE
  114    CONTINUE
  115 CONTINUE
  116 AR1 = 1.
      AI1 = 0.
      DO 120 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 117 IK=1,IDL1
            DO 1017 M=1,MP
            C2(M,IK,L) = CH2(M,IK,1)+AR1*CH2(M,IK,2)
            C2(M,IK,LC) = AI1*CH2(M,IK,IP)
 1017       CONTINUE
  117    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 119 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 118 IK=1,IDL1
               DO 1018 M=1,MP
               C2(M,IK,L) = C2(M,IK,L)+AR2*CH2(M,IK,J)
               C2(M,IK,LC) = C2(M,IK,LC)+AI2*CH2(M,IK,JC)
 1018          CONTINUE
  118       CONTINUE
  119    CONTINUE
  120 CONTINUE
      DO 122 J=2,IPPH
         DO 121 IK=1,IDL1
            DO 1021 M=1,MP
            CH2(M,IK,1) = CH2(M,IK,1)+CH2(M,IK,J)
 1021       CONTINUE
  121    CONTINUE
  122 CONTINUE
      DO 124 J=2,IPPH
         JC = IPP2-J
         DO 123 K=1,L1
            DO 1023 M=1,MP
            CH(M,1,K,J) = C1(M,1,K,J)-C1(M,1,K,JC)
            CH(M,1,K,JC) = C1(M,1,K,J)+C1(M,1,K,JC)
 1023       CONTINUE
  123    CONTINUE
  124 CONTINUE
      IF (IDO .EQ. 1) GO TO 132
      IF (NBD .LT. L1) GO TO 128
      DO 127 J=2,IPPH
         JC = IPP2-J
         DO 126 K=1,L1
            DO 125 I=3,IDO,2
               DO 1025 M=1,MP
               CH(M,I-1,K,J) = C1(M,I-1,K,J)-C1(M,I,K,JC)
               CH(M,I-1,K,JC) = C1(M,I-1,K,J)+C1(M,I,K,JC)
               CH(M,I,K,J) = C1(M,I,K,J)+C1(M,I-1,K,JC)
               CH(M,I,K,JC) = C1(M,I,K,J)-C1(M,I-1,K,JC)
 1025          CONTINUE
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      GO TO 132
  128 DO 131 J=2,IPPH
         JC = IPP2-J
         DO 130 I=3,IDO,2
            DO 129 K=1,L1
               DO 1029 M=1,MP
               CH(M,I-1,K,J) = C1(M,I-1,K,J)-C1(M,I,K,JC)
               CH(M,I-1,K,JC) = C1(M,I-1,K,J)+C1(M,I,K,JC)
               CH(M,I,K,J) = C1(M,I,K,J)+C1(M,I-1,K,JC)
               CH(M,I,K,JC) = C1(M,I,K,J)-C1(M,I-1,K,JC)
 1029          CONTINUE
  129       CONTINUE
  130    CONTINUE
  131 CONTINUE
  132 CONTINUE
      IF (IDO .EQ. 1) RETURN
      DO 133 IK=1,IDL1
         DO 1033 M=1,MP
         C2(M,IK,1) = CH2(M,IK,1)
 1033    CONTINUE
  133 CONTINUE
      DO 135 J=2,IP
         DO 134 K=1,L1
            DO 1034 M=1,MP
            C1(M,1,K,J) = CH(M,1,K,J)
 1034       CONTINUE
  134    CONTINUE
  135 CONTINUE
      IF (NBD .GT. L1) GO TO 139
      IS = -IDO
      DO 138 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 137 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 136 K=1,L1
               DO 1036 M=1,MP
               C1(M,I-1,K,J) = WA(IDIJ-1)*CH(M,I-1,K,J)-WA(IDIJ)*
     1          CH(M,I,K,J)
               C1(M,I,K,J) = WA(IDIJ-1)*CH(M,I,K,J)+WA(IDIJ)*
     1          CH(M,I-1,K,J)
 1036          CONTINUE
  136       CONTINUE
  137    CONTINUE
  138 CONTINUE
      GO TO 143
  139 IS = -IDO
      DO 142 J=2,IP
         IS = IS+IDO
         DO 141 K=1,L1
            IDIJ = IS
            DO 140 I=3,IDO,2
               IDIJ = IDIJ+2
               DO 1040 M=1,MP
               C1(M,I-1,K,J) = WA(IDIJ-1)*CH(M,I-1,K,J)-WA(IDIJ)*
     1          CH(M,I,K,J)
               C1(M,I,K,J) = WA(IDIJ-1)*CH(M,I,K,J)+WA(IDIJ)*
     1          CH(M,I-1,K,J)
 1040          CONTINUE
  140       CONTINUE
  141    CONTINUE
  142 CONTINUE
  143 RETURN
      END
      SUBROUTINE VRADF2 (MP,IDO,L1,CC,CH,MDIMC,WA1)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION   CH(MDIMC,IDO,2,L1)  ,CC(MDIMC,IDO,L1,2)     ,
     1                WA1(IDO)
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = CC(M,1,K,1)+CC(M,1,K,2)
         CH(M,IDO,2,K) = CC(M,1,K,1)-CC(M,1,K,2)
 1001    CONTINUE
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            DO 1003 M=1,MP
            CH(M,I,1,K) = CC(M,I,K,1)+(WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))
            CH(M,IC,2,K) = (WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))-CC(M,I,K,1)
            CH(M,I-1,1,K) = CC(M,I-1,K,1)+(WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))
            CH(M,IC-1,2,K) = CC(M,I-1,K,1)-(WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))
 1003       CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 DO 106 K=1,L1
         DO 1006 M=1,MP
         CH(M,1,2,K) = -CC(M,IDO,K,2)
         CH(M,IDO,1,K) = CC(M,IDO,K,1)
 1006    CONTINUE
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE VRADF3 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION   CH(MDIMC,IDO,3,L1)  ,CC(MDIMC,IDO,L1,3)     ,
     1                WA1(IDO)     ,WA2(IDO)
      ARG=2.*PIMACH(1.0)/3.
      TAUR=COS(ARG)
      TAUI=SIN(ARG)
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = CC(M,1,K,1)+(CC(M,1,K,2)+CC(M,1,K,3))
         CH(M,1,3,K) = TAUI*(CC(M,1,K,3)-CC(M,1,K,2))
         CH(M,IDO,2,K) = CC(M,1,K,1)+TAUR*
     1      (CC(M,1,K,2)+CC(M,1,K,3))
 1001    CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DO 1002 M=1,MP
            CH(M,I-1,1,K) = CC(M,I-1,K,1)+((WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))+(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3)))
            CH(M,I,1,K) = CC(M,I,K,1)+((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3)))
            CH(M,I-1,3,K) = (CC(M,I-1,K,1)+TAUR*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))))+(TAUI*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA2(I-2)*
     1       CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3))))
            CH(M,IC-1,2,K) = (CC(M,I-1,K,1)+TAUR*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))))-(TAUI*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA2(I-2)*
     1       CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3))))
            CH(M,I,3,K) = (CC(M,I,K,1)+TAUR*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))))+(TAUI*((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2))))
            CH(M,IC,2,K) = (TAUI*((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2))))-(CC(M,I,K,1)+TAUR*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))))
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
      SUBROUTINE VRADF4 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2,WA3)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION    CC(MDIMC,IDO,L1,4)   ,CH(MDIMC,IDO,4,L1)     ,
     1                WA1(IDO)     ,WA2(IDO)     ,WA3(IDO)
      HSQT2=SQRT(2.)/2.
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = (CC(M,1,K,2)+CC(M,1,K,4))
     1      +(CC(M,1,K,1)+CC(M,1,K,3))
         CH(M,IDO,4,K) = (CC(M,1,K,1)+CC(M,1,K,3))
     1      -(CC(M,1,K,2)+CC(M,1,K,4))
         CH(M,IDO,2,K) = CC(M,1,K,1)-CC(M,1,K,3)
         CH(M,1,3,K) = CC(M,1,K,4)-CC(M,1,K,2)
 1001    CONTINUE
  101 CONTINUE
      IF (IDO-2) 107,105,102
  102 IDP2 = IDO+2
      DO 104 K=1,L1
         DO 103 I=3,IDO,2
            IC = IDP2-I
            DO 1003 M=1,MP
            CH(M,I-1,1,K) = ((WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2))+(WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4)))+(CC(M,I-1,K,1)+(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))
            CH(M,IC-1,4,K) = (CC(M,I-1,K,1)+(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))-((WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))+(WA3(I-2)*CC(M,I-1,K,4)+
     1       WA3(I-1)*CC(M,I,K,4)))
            CH(M,I,1,K) = ((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))+(CC(M,I,K,1)+(WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3)))
            CH(M,IC,4,K) = ((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))-(CC(M,I,K,1)+(WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3)))
            CH(M,I-1,3,K) = ((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))+(CC(M,I-1,K,1)-(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))
            CH(M,IC-1,2,K) = (CC(M,I-1,K,1)-(WA2(I-2)*CC(M,I-1,K,3)+
     1       WA2(I-1)*CC(M,I,K,3)))-((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))
            CH(M,I,3,K) = ((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))+(CC(M,I,K,1)-(WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3)))
            CH(M,IC,2,K) = ((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))-(CC(M,I,K,1)-(WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3)))
 1003       CONTINUE
  103    CONTINUE
  104 CONTINUE
      IF (MOD(IDO,2) .EQ. 1) RETURN
  105 CONTINUE
      DO 106 K=1,L1
         DO 1006 M=1,MP
            CH(M,IDO,1,K) = (HSQT2*(CC(M,IDO,K,2)-CC(M,IDO,K,4)))+
     1       CC(M,IDO,K,1)
            CH(M,IDO,3,K) = CC(M,IDO,K,1)-(HSQT2*(CC(M,IDO,K,2)-
     1       CC(M,IDO,K,4)))
            CH(M,1,2,K) = (-HSQT2*(CC(M,IDO,K,2)+CC(M,IDO,K,4)))-
     1       CC(M,IDO,K,3)
            CH(M,1,4,K) = (-HSQT2*(CC(M,IDO,K,2)+CC(M,IDO,K,4)))+
     1       CC(M,IDO,K,3)
 1006    CONTINUE
  106 CONTINUE
  107 RETURN
      END
      SUBROUTINE VRADF5 (MP,IDO,L1,CC,CH,MDIMC,WA1,WA2,WA3,WA4)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION  CC(MDIMC,IDO,L1,5)    ,CH(MDIMC,IDO,5,L1)     ,
     1           WA1(IDO)     ,WA2(IDO)     ,WA3(IDO)     ,WA4(IDO)
      ARG=2.*PIMACH(1.0)/5.
      TR11=COS(ARG)
      TI11=SIN(ARG)
      TR12=COS(2.*ARG)
      TI12=SIN(2.*ARG)
      DO 101 K=1,L1
         DO 1001 M=1,MP
         CH(M,1,1,K) = CC(M,1,K,1)+(CC(M,1,K,5)+CC(M,1,K,2))+
     1    (CC(M,1,K,4)+CC(M,1,K,3))
         CH(M,IDO,2,K) = CC(M,1,K,1)+TR11*(CC(M,1,K,5)+CC(M,1,K,2))+
     1    TR12*(CC(M,1,K,4)+CC(M,1,K,3))
         CH(M,1,3,K) = TI11*(CC(M,1,K,5)-CC(M,1,K,2))+TI12*
     1    (CC(M,1,K,4)-CC(M,1,K,3))
         CH(M,IDO,4,K) = CC(M,1,K,1)+TR12*(CC(M,1,K,5)+CC(M,1,K,2))+
     1    TR11*(CC(M,1,K,4)+CC(M,1,K,3))
         CH(M,1,5,K) = TI12*(CC(M,1,K,5)-CC(M,1,K,2))-TI11*
     1    (CC(M,1,K,4)-CC(M,1,K,3))
 1001    CONTINUE
  101 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IDP2 = IDO+2
      DO 103 K=1,L1
         DO 102 I=3,IDO,2
            IC = IDP2-I
            DO 1002 M=1,MP
            CH(M,I-1,1,K) = CC(M,I-1,K,1)+((WA1(I-2)*CC(M,I-1,K,2)+
     1       WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*
     1       CC(M,I,K,5)))+((WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))+(WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4)))
            CH(M,I,1,K) = CC(M,I,K,1)+((WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*
     1       CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4)))
            CH(M,I-1,3,K) = CC(M,I-1,K,1)+TR11*
     1      ( WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2)
     1       +WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5))+TR12*
     1      ( WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3)
     1       +WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))+TI11*
     1      ( WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2)
     1       -(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*CC(M,I-1,K,5)))+TI12*
     1      ( WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3)
     1       -(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*CC(M,I-1,K,4)))
            CH(M,IC-1,2,K) = CC(M,I-1,K,1)+TR11*
     1      ( WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2)
     1       +WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5))+TR12*
     1     ( WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3)
     1      +WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))-(TI11*
     1      ( WA1(I-2)*CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2)
     1       -(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*CC(M,I-1,K,5)))+TI12*
     1      ( WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*CC(M,I-1,K,3)
     1       -(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*CC(M,I-1,K,4))))
            CH(M,I,3,K) = (CC(M,I,K,1)+TR11*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR12*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))+(TI11*((WA4(I-2)*CC(M,I-1,K,5)+
     1       WA4(I-1)*CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))+TI12*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))
            CH(M,IC,2,K) = (TI11*((WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*
     1       CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))+TI12*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))-(CC(M,I,K,1)+TR11*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR12*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
            CH(M,I-1,5,K) = (CC(M,I-1,K,1)+TR12*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*
     1       CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5)))+TR11*((WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))+(WA3(I-2)*
     1       CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))))+(TI12*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA4(I-2)*CC(M,I,K,5)-
     1       WA4(I-1)*CC(M,I-1,K,5)))-TI11*((WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
            CH(M,IC-1,4,K) = (CC(M,I-1,K,1)+TR12*((WA1(I-2)*
     1       CC(M,I-1,K,2)+WA1(I-1)*CC(M,I,K,2))+(WA4(I-2)*
     1       CC(M,I-1,K,5)+WA4(I-1)*CC(M,I,K,5)))+TR11*((WA2(I-2)*
     1       CC(M,I-1,K,3)+WA2(I-1)*CC(M,I,K,3))+(WA3(I-2)*
     1       CC(M,I-1,K,4)+WA3(I-1)*CC(M,I,K,4))))-(TI12*((WA1(I-2)*
     1       CC(M,I,K,2)-WA1(I-1)*CC(M,I-1,K,2))-(WA4(I-2)*CC(M,I,K,5)-
     1       WA4(I-1)*CC(M,I-1,K,5)))-TI11*((WA2(I-2)*CC(M,I,K,3)-
     1       WA2(I-1)*CC(M,I-1,K,3))-(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
            CH(M,I,5,K) = (CC(M,I,K,1)+TR12*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR11*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))+(TI12*((WA4(I-2)*CC(M,I-1,K,5)+
     1       WA4(I-1)*CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))-TI11*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))
            CH(M,IC,4,K) = (TI12*((WA4(I-2)*CC(M,I-1,K,5)+WA4(I-1)*
     1       CC(M,I,K,5))-(WA1(I-2)*CC(M,I-1,K,2)+WA1(I-1)*
     1       CC(M,I,K,2)))-TI11*((WA3(I-2)*CC(M,I-1,K,4)+WA3(I-1)*
     1       CC(M,I,K,4))-(WA2(I-2)*CC(M,I-1,K,3)+WA2(I-1)*
     1       CC(M,I,K,3))))-(CC(M,I,K,1)+TR12*((WA1(I-2)*CC(M,I,K,2)-
     1       WA1(I-1)*CC(M,I-1,K,2))+(WA4(I-2)*CC(M,I,K,5)-WA4(I-1)*
     1       CC(M,I-1,K,5)))+TR11*((WA2(I-2)*CC(M,I,K,3)-WA2(I-1)*
     1       CC(M,I-1,K,3))+(WA3(I-2)*CC(M,I,K,4)-WA3(I-1)*
     1       CC(M,I-1,K,4))))
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      RETURN
      END
      SUBROUTINE VRADFG (MP,IDO,IP,L1,IDL1,CC,C1,C2,CH,CH2,MDIMC,WA)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION     CH(MDIMC,IDO,L1,IP)   ,CC(MDIMC,IDO,IP,L1)  ,
     1            C1(MDIMC,IDO,L1,IP)    ,C2(MDIMC,IDL1,IP),
     2                CH2(MDIMC,IDL1,IP)           ,WA(IDO)
      TPI=2.*PIMACH(1.0)
      ARG = TPI/FLOAT(IP)
      DCP = COS(ARG)
      DSP = SIN(ARG)
      IPPH = (IP+1)/2
      IPP2 = IP+2
      IDP2 = IDO+2
      NBD = (IDO-1)/2
      IF (IDO .EQ. 1) GO TO 119
      DO 101 IK=1,IDL1
         DO 1001 M=1,MP
         CH2(M,IK,1) = C2(M,IK,1)
 1001    CONTINUE
  101 CONTINUE
      DO 103 J=2,IP
         DO 102 K=1,L1
            DO 1002 M=1,MP
            CH(M,1,K,J) = C1(M,1,K,J)
 1002       CONTINUE
  102    CONTINUE
  103 CONTINUE
      IF (NBD .GT. L1) GO TO 107
      IS = -IDO
      DO 106 J=2,IP
         IS = IS+IDO
         IDIJ = IS
         DO 105 I=3,IDO,2
            IDIJ = IDIJ+2
            DO 104 K=1,L1
               DO 1004 M=1,MP
               CH(M,I-1,K,J) = WA(IDIJ-1)*C1(M,I-1,K,J)+WA(IDIJ)
     1           *C1(M,I,K,J)
               CH(M,I,K,J) = WA(IDIJ-1)*C1(M,I,K,J)-WA(IDIJ)
     1           *C1(M,I-1,K,J)
 1004          CONTINUE
  104       CONTINUE
  105    CONTINUE
  106 CONTINUE
      GO TO 111
  107 IS = -IDO
      DO 110 J=2,IP
         IS = IS+IDO
         DO 109 K=1,L1
            IDIJ = IS
            DO 108 I=3,IDO,2
               IDIJ = IDIJ+2
               DO 1008 M=1,MP
               CH(M,I-1,K,J) = WA(IDIJ-1)*C1(M,I-1,K,J)+WA(IDIJ)
     1           *C1(M,I,K,J)
               CH(M,I,K,J) = WA(IDIJ-1)*C1(M,I,K,J)-WA(IDIJ)
     1           *C1(M,I-1,K,J)
 1008          CONTINUE
  108       CONTINUE
  109    CONTINUE
  110 CONTINUE
  111 IF (NBD .LT. L1) GO TO 115
      DO 114 J=2,IPPH
         JC = IPP2-J
         DO 113 K=1,L1
            DO 112 I=3,IDO,2
               DO 1012 M=1,MP
               C1(M,I-1,K,J) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               C1(M,I-1,K,JC) = CH(M,I,K,J)-CH(M,I,K,JC)
               C1(M,I,K,J) = CH(M,I,K,J)+CH(M,I,K,JC)
               C1(M,I,K,JC) = CH(M,I-1,K,JC)-CH(M,I-1,K,J)
 1012          CONTINUE
  112       CONTINUE
  113    CONTINUE
  114 CONTINUE
      GO TO 121
  115 DO 118 J=2,IPPH
         JC = IPP2-J
         DO 117 I=3,IDO,2
            DO 116 K=1,L1
               DO 1016 M=1,MP
               C1(M,I-1,K,J) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               C1(M,I-1,K,JC) = CH(M,I,K,J)-CH(M,I,K,JC)
               C1(M,I,K,J) = CH(M,I,K,J)+CH(M,I,K,JC)
               C1(M,I,K,JC) = CH(M,I-1,K,JC)-CH(M,I-1,K,J)
 1016          CONTINUE
  116       CONTINUE
  117    CONTINUE
  118 CONTINUE
      GO TO 121
  119 DO 120 IK=1,IDL1
         DO 1020 M=1,MP
         C2(M,IK,1) = CH2(M,IK,1)
 1020    CONTINUE
  120 CONTINUE
  121 DO 123 J=2,IPPH
         JC = IPP2-J
         DO 122 K=1,L1
            DO 1022 M=1,MP
            C1(M,1,K,J) = CH(M,1,K,J)+CH(M,1,K,JC)
            C1(M,1,K,JC) = CH(M,1,K,JC)-CH(M,1,K,J)
 1022       CONTINUE
  122    CONTINUE
  123 CONTINUE
C
      AR1 = 1.
      AI1 = 0.
      DO 127 L=2,IPPH
         LC = IPP2-L
         AR1H = DCP*AR1-DSP*AI1
         AI1 = DCP*AI1+DSP*AR1
         AR1 = AR1H
         DO 124 IK=1,IDL1
            DO 1024 M=1,MP
            CH2(M,IK,L) = C2(M,IK,1)+AR1*C2(M,IK,2)
            CH2(M,IK,LC) = AI1*C2(M,IK,IP)
 1024       CONTINUE
  124    CONTINUE
         DC2 = AR1
         DS2 = AI1
         AR2 = AR1
         AI2 = AI1
         DO 126 J=3,IPPH
            JC = IPP2-J
            AR2H = DC2*AR2-DS2*AI2
            AI2 = DC2*AI2+DS2*AR2
            AR2 = AR2H
            DO 125 IK=1,IDL1
               DO 1025 M=1,MP
               CH2(M,IK,L) = CH2(M,IK,L)+AR2*C2(M,IK,J)
               CH2(M,IK,LC) = CH2(M,IK,LC)+AI2*C2(M,IK,JC)
 1025          CONTINUE
  125       CONTINUE
  126    CONTINUE
  127 CONTINUE
      DO 129 J=2,IPPH
         DO 128 IK=1,IDL1
            DO 1028 M=1,MP
            CH2(M,IK,1) = CH2(M,IK,1)+C2(M,IK,J)
 1028       CONTINUE
  128    CONTINUE
  129 CONTINUE
C
      IF (IDO .LT. L1) GO TO 132
      DO 131 K=1,L1
         DO 130 I=1,IDO
            DO 1030 M=1,MP
            CC(M,I,1,K) = CH(M,I,K,1)
 1030       CONTINUE
  130    CONTINUE
  131 CONTINUE
      GO TO 135
  132 DO 134 I=1,IDO
         DO 133 K=1,L1
            DO 1033 M=1,MP
            CC(M,I,1,K) = CH(M,I,K,1)
 1033       CONTINUE
  133    CONTINUE
  134 CONTINUE
  135 DO 137 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 136 K=1,L1
            DO 1036 M=1,MP
            CC(M,IDO,J2-2,K) = CH(M,1,K,J)
            CC(M,1,J2-1,K) = CH(M,1,K,JC)
 1036       CONTINUE
  136    CONTINUE
  137 CONTINUE
      IF (IDO .EQ. 1) RETURN
      IF (NBD .LT. L1) GO TO 141
      DO 140 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 139 K=1,L1
            DO 138 I=3,IDO,2
               IC = IDP2-I
               DO 1038 M=1,MP
               CC(M,I-1,J2-1,K) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               CC(M,IC-1,J2-2,K) = CH(M,I-1,K,J)-CH(M,I-1,K,JC)
               CC(M,I,J2-1,K) = CH(M,I,K,J)+CH(M,I,K,JC)
               CC(M,IC,J2-2,K) = CH(M,I,K,JC)-CH(M,I,K,J)
 1038          CONTINUE
  138       CONTINUE
  139    CONTINUE
  140 CONTINUE
      RETURN
  141 DO 144 J=2,IPPH
         JC = IPP2-J
         J2 = J+J
         DO 143 I=3,IDO,2
            IC = IDP2-I
            DO 142 K=1,L1
               DO 1042 M=1,MP
               CC(M,I-1,J2-1,K) = CH(M,I-1,K,J)+CH(M,I-1,K,JC)
               CC(M,IC-1,J2-2,K) = CH(M,I-1,K,J)-CH(M,I-1,K,JC)
               CC(M,I,J2-1,K) = CH(M,I,K,J)+CH(M,I,K,JC)
               CC(M,IC,J2-2,K) = CH(M,I,K,JC)-CH(M,I,K,J)
 1042          CONTINUE
  142       CONTINUE
  143    CONTINUE
  144 CONTINUE
      RETURN
      END
      SUBROUTINE VRFFTB(M,N,R,RT,MDIMR,WSAVE)
C***BEGIN PROLOGUE  VRFFTB
C***DATE WRITTEN   850801   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM, 
C             FOURIER SYNTHESIS, BACKWARD TRANSFORM, MULTIPLE SEQUENCES
C***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
C***PURPOSE  Backward real periodic transform, M sequences.
C***DESCRIPTION
C
C  Subroutine VRFFTB computes the synthesis (backward transform) of a
C  number of real periodic sequences from their Fourier coefficients. 
C  Specifically, for each set of independent Fourier coefficients
C  F(K), the corresponding real periodic sequence is computed. 
C
C  The array WSAVE which is used by subroutine VRFFTB must be
C  initialized by calling subroutine VRFFTI(N,WSAVE).
C
C
C  Input Parameters
C
C  M       the number of sets of coefficients.
C
C  N       the length of the sequences of coefficients to be 
C          transformed.  The method is most efficient when N is a
C          product of small primes, however n may be any positive 
C          integer.
C
C  R       areal two-dimensional array of size MDIMX x N containing the
C          coefficients to be transformed.  Each set of coefficients
C          F(K), K\0,1,..,N-1, is stored as a ROW of R.  Specifically,
C          the I-th set of independent Fourier coefficients is stored
C
C                R(I,1) = REAL( F(I,0) ),
C
C                R(I,2*K) = REAL( F(I,K) )
C
C                R(I,2*K+1) = IMAG( F(I,K) )
C
C                   for K = 1, 2, . . . , M-1,
C
C                and, when N is even,
C
C                R(I,N) = REAL( F(I,N/2) ).
C
C  RT      a real two-dimensional work array of size MDIMX x N.
C
C  MDIMR   the row (or first) dimension of the arrays R and RT exactly 
C          as they appear in the calling program.  This parameter is 
C          used to specify the variable dimension of these arrays.
C
C  WSAVE   a real one-dimensional work array which must be dimensioned
C          at least N+15.  The WSAVE array must be initialized by 
C          calling subroutine VRFFTI.  A different WSAVE array must be
C          used for each different value of N.  This initialization does
C          not have to be repeated so long as N remains unchanged.  The
C          same WSAVE array may be used by VRFFTB and VRFFTB.
C
C  Output Parameters
C
C  R       contains M real periodic sequences corresponding to the given
C          coefficients.  Specifically, the I-th row of R contains the 
C          real periodic sequence corresponding to the I-th set of
C          independent Fourier coefficients F(I,K) stored as
C
C               R(I,J) = X(I,J-1) ,   J = 1, 2, . . . , N, where
C
C               X(I,J) = SQRT(1/N)* F(I,0) + (-1)**J*F(I,N/2)
C                        + 2*SUM(K=1,M)[ REAL(F(I,2K))*COS(2K*J*PI/N)
C                        - IMAG(F(I,2K+1))*SIN(2K*J*PI/N) ]  ,
C
C                 when N is even, and
C
C               X(I,J) = SQRT(1/N)* F(I,0) +
C                        2*SUM(K=1,M)[ REAL(F(I,2K))*COS(2K*J*PI/N)
C                        - IMAG(F(I,2K+1))*SIN(2K*J*PI/N) ]  ,
C
C                 when N is odd.
C
C  WSAVE   contains results which must not be destroyed between calls
C          to VRFFTF or VRFFTB.
C
C  -----------------------------------------------------------------
C
C  NOTE  -  A call of VRFFTF followed immediately by a call of
C           of VRFFTB will return the original sequences R.  Thus,
C           VRFFTB is the correctly normalized inverse of VRFFTF.
C
C  -----------------------------------------------------------------
C
C  VRFFTB is a straightforward extension of the subprogram RFFTB to
C  handle M simultaneous sequences.  RFFTB was originally developed
C  by P. N. Swarztrauber of NCAR.
C
C
C              * * * * * * * * * * * * * * * * * * * * *
C              *                                       *
C              *         PROGRAM SPECIFICATIONS        *
C              *                                       *
C              * * * * * * * * * * * * * * * * * * * * *
C
C
C     DIMENSION OF    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
C     ARGUMENTS
C
C     LATEST          AUGUST 1, 1985
C     REVISION
C
C     SUBPROGRAMS     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
C     REQUIRED        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
C                     VRADB3, VRADB4, VRADB5, VRADBG, PIMACH
C
C     SPECIAL         NONE
C     CONDITIONS
C
C     COMMON          NONE
C     BLOCKS
C
C     I/O             NONE
C
C     PRECISION       SINGLE
C
C     SPECIALIST      ROLAND SWEET
C
C     LANGUAGE        FORTRAN
C
C     HISTORY         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
C                     NATIONAL BUREAU OF STANDARDS (BOULDER).
C
C     ALGORITHM       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
C                     OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.
C
C     PORTABILITY     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
C                     THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
C                     THE FUNCTION PIMACH.
C
C     REQUIRED        COS,SIN
C     RESIDENT
C     ROUTINES
C
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFTB1
C***END PROLOGUE  VRFFTB
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION     R(MDIMR,N),RT(MDIMR,N),WSAVE(N+15)
      IF (N .EQ. 1) RETURN
      CALL VRFTB1 (M,N,R,RT,MDIMR,WSAVE(1),WSAVE(N+1))
      RETURN
      END
      SUBROUTINE VRFFTF (M,N,R,RT,MDIMR,WSAVE)
C***BEGIN PROLOGUE  VRFFTF
C***DATE WRITTEN   850801   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM, 
C             FOURIER ANALYSIS, FORWARD TRANSFORM, MULTIPLE SEQUENCES
C***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
C***PURPOSE  Forward real periodic transform, M sequences.
C***DESCRIPTION
C
C  Subroutine VRFFTF computes the Fourier coefficients (forward 
C  transform) of a number of real periodic sequences.  Specifically,
C  for each sequence the subroutine claculates the independent
C  Fourier coefficients described below at output parameter R.
C
C  The array WSAVE which is used by subroutine VRFFTF must be
C  initialized by calling subroutine VRFFTI(N,WSAVE).
C
C
C  Input Parameters
C
C  M       the number of sequences to be transformed.
C
C  N       the length of the sequences to be transformed.  The method
C          is most efficient when N is a product of small primes,
C          however n may be any positive integer.
C
C  R       areal two-dimensional array of size MDIMX x N containing the
C          the sequences to be transformed.  The sequences are stored
C          in the ROWS of R.  Thus, the I-th sequence to be transformed,
C          X(I,J), J=0,1,...,N-1, is stored as
C
C               R(I,J) = X(I,J-1) , J=1, 2, . . . , N.
C
C  RT      a real two-dimensional work array of size MDIMX x N.
C
C  MDIMR   the row (or first) dimension of the arrays R and RT exactly 
C          as they appear in the calling program.  This parameter is 
C          used to specify the variable dimension of these arrays.
C
C  WSAVE   a real one-dimensional work array which must be dimensioned
C          at least N+15.  The WSAVE array must be initialized by 
C          calling subroutine VRFFTI.  A different WSAVE array must be
C          used for each different value of N.  This initialization does
C          not have to be repeated so long as N remains unchanged.  The
C          same WSAVE array may be used by VRFFTF and VRFFTB.
C
C  Output Parameters
C
C  R       contains the Fourier coefficients F(K) for each of the M 
C          input sequences.  Specifically, row I of R, R(I,J), 
C          J=1,2,..,N, contains the independent Fourier coefficients
C          F(I,K), for the I-th input sequence stored as
C
C             R(I,1) = REAL( F(I,0) ),
C                    = SQRT(1/N)*SUM(J=0,N-1)[ X(I,J) ],
C
C             R(I,2*K) = REAL( F(I,K) )
C                      = SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*COS(2J*K*PI/N)]
C
C             R(I,2*K+1) = IMAG( F(I,K) )
C                        =-SQRT(1/N)*SUM(J=0,N-1)[X(I,J)*SIN(2J*K*PI/N)]
C
C                   for K = 1, 2, . . . , M-1,
C
C              and, when N is even,
C
C              R(I,N) = REAL( F(I,N/2) ).
C                     = SQRT(1/N)*SUM(J=0,N-1)[ (-1)**J*X(I,J) ].
C
C  WSAVE   contains results which must not be destroyed between calls
C          to VRFFTF or VRFFTB.
C
C  -----------------------------------------------------------------
C
C  NOTE  -  A call of VRFFTF followed immediately by a call of
C           of VRFFTB will return the original sequences R.  Thus,
C           VRFFTB is the correctly normalized inverse of VRFFTF.
C
C  -----------------------------------------------------------------
C
C  VRFFTF is a straightforward extension of the subprogram RFFTF to
C  handle M simultaneous sequences.  RFFTF was originally developed
C  by P. N. Swarztrauber of NCAR.
C
C
C              * * * * * * * * * * * * * * * * * * * * *
C              *                                       *
C              *         PROGRAM SPECIFICATIONS        *
C              *                                       *
C              * * * * * * * * * * * * * * * * * * * * *
C
C
C     DIMENSION OF    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
C     ARGUMENTS
C
C     LATEST          AUGUST 1, 1985
C     REVISION
C
C     SUBPROGRAMS     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
C     REQUIRED        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
C                     VRADB3, VRADB4, VRADB5, VRADBG, PIMACH
C
C     SPECIAL         NONE
C     CONDITIONS
C
C     COMMON          NONE
C     BLOCKS
C
C     I/O             NONE
C
C     PRECISION       SINGLE
C
C     SPECIALIST      ROLAND SWEET
C
C     LANGUAGE        FORTRAN
C
C     HISTORY         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
C                     NATIONAL BUREAU OF STANDARDS (BOULDER).
C
C     ALGORITHM       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
C                     OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.
C
C     PORTABILITY     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
C                     THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
C                     THE FUNCTION PIMACH.
C
C     REQUIRED        COS,SIN
C     RESIDENT
C     ROUTINES
C
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFTF1
C***END PROLOGUE  VRFFTF
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION       R(MDIMR,N)  ,RT(MDIMR,N)    ,WSAVE(N+15)
C***FIRST EXECUTABLE STATEMENT  VRFFTF
      IF (N .EQ. 1) RETURN
      CALL VRFTF1 (M,N,R,RT,MDIMR,WSAVE(1),WSAVE(N+1))
      RETURN
      END
      SUBROUTINE VRFFTI (N,WSAVE)
C***BEGIN PROLOGUE  VRFFTI
C***DATE WRITTEN   860701   (YYMMDD)
C***REVISION DATE  900509   (YYMMDD)
C***CATEGORY NO.  J1A1
C***KEYWORDS  FAST FOURIER TRANSFORM, REAL PERIODIC TRANSFORM,
C             MULTIPLE SEQUENCES
C***AUTHOR  SWEET, R.A. (NIST) AND LINDGREN, L.L. (NIST)
C***PURPOSE  Initialization for VRFFTF and VRFFTB.
C***DESCRIPTION
C
C  Subroutine VRFFTI initializes the array WSAVE which is used in
C  both VRFFTF and VRFFTB.  The prime factorization of N together with
C  a tabulation of certain trigonometric functions are computed and
C  stored in the array WSAVE.
C
C  Input Parameter
C
C  N       the length of the sequence to be transformed.  There is no
C          restriction on N.
C
C  Output Parameter
C
C  WSAVE   a work array which must be dimensioned at least N+15.
C          The same work array can be used for both VRFFTF and VRFFTB
C          as long as N remains unchanged.  Different WSAVE arrays
C          are required for different values of N.  The contents of
C          WSAVE must not be changed between calls of VRFFTF or VRFFTB.
C
C
C              * * * * * * * * * * * * * * * * * * * * *
C              *                                       *
C              *         PROGRAM SPECIFICATIONS        *
C              *                                       *
C              * * * * * * * * * * * * * * * * * * * * *
C
C
C     DIMENSION OF    R(MDIMR,N), RT(MDIMR,N), WSAVE(N+15)
C     ARGUMENTS
C
C     LATEST          AUGUST 1, 1985
C     REVISION
C
C     SUBPROGRAMS     VRFFTI, VRFTI1, VRFFTF, VRFTF1, VRADF2, VRADF3,
C     REQUIRED        VRADF4, VRADF5, VRADFG, VRFFTB, VRFTB1, VRADB2,
C                     VRADB3, VRADB4, VRADB5, VRADBG, PIMACH
C
C     SPECIAL         NONE
C     CONDITIONS
C
C     COMMON          NONE
C     BLOCKS
C
C     I/O             NONE
C
C     PRECISION       SINGLE
C
C     SPECIALIST      ROLAND SWEET
C
C     LANGUAGE        FORTRAN
C
C     HISTORY         WRITTEN BY LINDA LINDGREN AND ROLAND SWEET AT THE
C                     NATIONAL BUREAU OF STANDARDS (BOULDER).
C
C     ALGORITHM       A REAL VARIANT OF THE STOCKHAM AUTOSORT VERSION
C                     OF THE COOLEY-TUKEY FAST FOURIER TRANSFORM.
C
C     PORTABILITY     AMERICAN NATIONAL STANDARDS INSTITUTE FORTRAN 77.
C                     THE ONLY MACHINE DEPENDENT CONSTANT IS LOCATED IN
C                     THE FUNCTION PIMACH.
C
C     REQUIRED        COS,SIN
C     RESIDENT
C     ROUTINES
C
C
C***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
C               Computations, (G. Rodrigue, ed.), Academic Press, 1982,
C               pp. 51-83.
C***ROUTINES CALLED  VRFTI1
C***END PROLOGUE  VRFFTI
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION       WSAVE(N+15)
C***FIRST EXECUTABLE STATEMENT  VRFFTI
      IF (N .EQ. 1) RETURN
      CALL VRFTI1 (N,WSAVE(1),WSAVE(N+1))
      RETURN
      END
      SUBROUTINE VRFTB1 (M,N,C,CH,MDIMC,WA,FAC)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION       CH(MDIMC,N), C(MDIMC,N), WA(N) ,FAC(15)
      NF = FAC(2)
      NA = 0
      L1 = 1
      IW = 1
      DO 116 K1=1,NF
         IP = FAC(K1+2)
         L2 = IP*L1
         IDO = N/L2
         IDL1 = IDO*L1
         IF (IP .NE. 4) GO TO 103
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL VRADB4 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2),WA(IX3))
         GO TO 102
  101    CALL VRADB4 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
         GO TO 115
  103    IF (IP .NE. 2) GO TO 106
         IF (NA .NE. 0) GO TO 104
         CALL VRADB2 (M,IDO,L1,C,CH,MDIMC,WA(IW))
         GO TO 105
  104    CALL VRADB2 (M,IDO,L1,CH,C,MDIMC,WA(IW))
  105    NA = 1-NA
         GO TO 115
  106    IF (IP .NE. 3) GO TO 109
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 107
         CALL VRADB3 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2))
         GO TO 108
  107    CALL VRADB3 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2))
  108    NA = 1-NA
         GO TO 115
  109    IF (IP .NE. 5) GO TO 112
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 110
      CALL VRADB5 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 111
  110 CALL VRADB5 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
         GO TO 115
  112    IF (NA .NE. 0) GO TO 113
         CALL VRADBG (M,IDO,IP,L1,IDL1,C,C,C,CH,CH,MDIMC,WA(IW))
         GO TO 114
  113    CALL VRADBG (M,IDO,IP,L1,IDL1,CH,CH,CH,C,C,MDIMC,WA(IW))
  114    IF (IDO .EQ. 1) NA = 1-NA
  115    L1 = L2
         IW = IW+(IP-1)*IDO
  116 CONTINUE
      SCALE=SQRT(1./N)
      IF (NA .EQ. 0) GO TO 118
      DO 117 J=1,N
      DO 117 I=1,M
         C(I,J) = SCALE*CH(I,J)
  117 CONTINUE
      RETURN
  118 DO 119 J=1,N
      DO 119 I=1,M
         C(I,J)=SCALE*C(I,J)
  119 CONTINUE
      RETURN
      END
      SUBROUTINE VRFTF1 (M,N,C,CH,MDIMC,WA,FAC)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION       CH(MDIMC,N) ,C(MDIMC,N)  ,WA(N)   ,FAC(15)
      NF = FAC(2)
      NA = 1
      L2 = N
      IW = N
      DO 111 K1=1,NF
         KH = NF-K1
         IP = FAC(KH+3)
         L1 = L2/IP
         IDO = N/L2
         IDL1 = IDO*L1
         IW = IW-(IP-1)*IDO
         NA = 1-NA
         IF (IP .NE. 4) GO TO 102
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IF (NA .NE. 0) GO TO 101
         CALL VRADF4 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  101    CALL VRADF4 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2),WA(IX3))
         GO TO 110
  102    IF (IP .NE. 2) GO TO 104
         IF (NA .NE. 0) GO TO 103
         CALL VRADF2 (M,IDO,L1,C,CH,MDIMC,WA(IW))
         GO TO 110
  103    CALL VRADF2 (M,IDO,L1,CH,C,MDIMC,WA(IW))
         GO TO 110
  104    IF (IP .NE. 3) GO TO 106
         IX2 = IW+IDO
         IF (NA .NE. 0) GO TO 105
         CALL VRADF3 (M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2))
         GO TO 110
  105    CALL VRADF3 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2))
         GO TO 110
  106    IF (IP .NE. 5) GO TO 108
         IX2 = IW+IDO
         IX3 = IX2+IDO
         IX4 = IX3+IDO
         IF (NA .NE. 0) GO TO 107
      CALL VRADF5(M,IDO,L1,C,CH,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  107 CALL VRADF5 (M,IDO,L1,CH,C,MDIMC,WA(IW),WA(IX2),WA(IX3),WA(IX4))
         GO TO 110
  108    IF (IDO .EQ. 1) NA = 1-NA
         IF (NA .NE. 0) GO TO 109
         CALL VRADFG (M,IDO,IP,L1,IDL1,C,C,C,CH,CH,MDIMC,WA(IW))
         NA = 1
         GO TO 110
  109    CALL VRADFG (M,IDO,IP,L1,IDL1,CH,CH,CH,C,C,MDIMC,WA(IW))
         NA = 0
  110    L2 = L1
  111 CONTINUE
      SCALE=SQRT(1./N)
      IF (NA .EQ. 1) GO TO 113
      DO 112 J=1,N
      DO 112 I=1,M
         C(I,J) = SCALE*CH(I,J)
  112 CONTINUE
      RETURN
  113 DO 114 J=1,N
      DO 114 I=1,M
         C(I,J)=SCALE*C(I,J)
  114 CONTINUE
      RETURN
      END
      SUBROUTINE VRFTI1 (N,WA,FAC)
C
C     VRFFTPK, VERSION 1, AUGUST 1985
C
      DIMENSION       WA(N)      ,FAC(15)    ,NTRYH(4)
      DATA NTRYH(1),NTRYH(2),NTRYH(3),NTRYH(4)/4,2,3,5/
      NL = N
      NF = 0
      J = 0
  101 J = J+1
      IF (J-4) 102,102,103
  102 NTRY = NTRYH(J)
      GO TO 104
  103 NTRY = NTRY+2
  104 NQ = NL/NTRY
      NR = NL-NTRY*NQ
      IF (NR) 101,105,101
  105 NF = NF+1
      FAC(NF+2) = NTRY
      NL = NQ
      IF (NTRY .NE. 2) GO TO 107
      IF (NF .EQ. 1) GO TO 107
      DO 106 I=2,NF
         IB = NF-I+2
         FAC(IB+2) = FAC(IB+1)
  106 CONTINUE
      FAC(3) = 2
  107 IF (NL .NE. 1) GO TO 104
      FAC(1) = N
      FAC(2) = NF
      TPI = 2.*PIMACH(1.0)
      ARGH = TPI/FLOAT(N)
      IS = 0
      NFM1 = NF-1
      L1 = 1
      IF (NFM1 .EQ. 0) RETURN
      DO 110 K1=1,NFM1
         IP = FAC(K1+2)
         LD = 0
         L2 = L1*IP
         IDO = N/L2
         IPM = IP-1
         DO 109 J=1,IPM
            LD = LD+L1
            I = IS
            ARGLD = FLOAT(LD)*ARGH
            FI = 0.
            DO 108 II=3,IDO,2
               I = I+2
               FI = FI+1.
               ARG = FI*ARGLD
               WA(I-1) = COS(ARG)
               WA(I) = SIN(ARG)
  108       CONTINUE
            IS = IS+IDO
  109    CONTINUE
         L1 = L2
  110 CONTINUE
      RETURN
      END

      subroutine cnvstr(ing,cha)
      integer intg,i1,i2,i3,ing
      character*5 cha
      character c1,c2,c3,c4,c5
      intg=ing
      i1=intg/10000
      call cn(i1,c1)
      intg=intg-i1*10000
      i1=intg/1000
      call cn(i1,c2)
      intg=intg-i1*1000
      i1=intg/100
      call cn(i1,c3)
      intg=intg-i1*100
      i1=intg/10
      call cn(i1,c4)
      intg=intg-i1*10
      i1=intg
      call cn(i1,c5)
      cha=c1//c2//c3//c4//c5
      end

      subroutine cn(inte,ch)
      integer inte
      character ch
      if (inte .eq. 0) ch='0'
      if (inte .eq. 1) ch='1'
      if (inte .eq. 2) ch='2'
      if (inte .eq. 3) ch='3'
      if (inte. eq. 4) ch='4'
      if (inte. eq. 5) ch='5'
      if (inte. eq. 6) ch='6'
      if (inte. eq. 7) ch='7'
      if (inte. eq. 8) ch='8'
      if (inte. eq. 9) ch='9'
      end

      SUBROUTINE TRID (MR,A,B,C,Y,D)
      DIMENSION       A(1)       ,B(1)       ,C(1)       ,Y(1)       ,
     1                D(1)
      call dgtsv(mr,1,a(2),b,c,y,mr,info)
      if (info.ne.0) stop 'error in tridagional systems ' 
      END  

