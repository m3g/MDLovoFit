!
! The following subroutines were not written by me 
!

      subroutine jacobi(a,v,n) 
!*********************************************************************** 
! 
!     Diagonalisation of real symmetric matices by Jacobi method 
! 
!*********************************************************************** 
      implicit double precision(a-h,o-z) 
      integer i, j, k, n, m 
      double precision v1, v2, v3, u, omg, c, rho, tes, scl, tem, s
      double precision a(4,4),v(4,4)
      rho=1.0d-12 
      tes=0. 
      scl=0. 
      do 10 i=1,n 
   10 scl=scl+a(i,i)**2 
      scl=dsqrt(scl)/dfloat(n) 
      do 20 i=1,n 
      do 20 j=1,n 
   20 a(i,j)=a(i,j)/scl 
      do 30 i=1,n 
      do 30 j=1,n 
      v(i,j)=0. 
      if(i.eq.j)v(i,j)=1. 
   30 continue 
      do 100 i=2,n 
      do 100 j=1,i-1 
  100 tes=tes+2.*a(i,j)*a(i,j) 
      tes=dsqrt(tes) 
      m=0 
  105 tes=tes/dfloat(n) 
      if(tes.lt.rho)tes=rho 
  110 do 165 i=2,n 
      do 165 j=1,i-1 
      if(abs(a(i,j))-tes)165,115,115 
  115 m=1 
      v1=a(j,j) 
      v2=a(i,j) 
      v3=a(i,i) 
      u=0.5*(v1-v3) 
      if(abs(u)-rho)120,125,125 
  120 omg=-1. 
      go to 130 
  125 omg=-v2/dsqrt(v2*v2+u*u) 
      if(u.lt.0.)omg=-omg 
  130 s=omg/dsqrt(2.*(1.+dsqrt(1.-omg*omg))) 
      c=dsqrt(1.-s*s) 
      do 160 k=1,n 
      if(k-i)140,135,135 
  135 tem=a(k,j)*c-a(k,i)*s 
      a(k,i)=a(k,j)*s+a(k,i)*c 
      a(k,j)=tem 
      go to 155 
  140 if(k-j)145,150,150 
  145 tem=a(j,k)*c-a(i,k)*s 
      a(i,k)=a(j,k)*s+a(i,k)*c 
      a(j,k)=tem 
      go to 155 
  150 tem=a(k,j)*c-a(i,k)*s 
      a(i,k)=a(k,j)*s+a(i,k)*c 
      a(k,j)=tem 
  155 tem=v(k,j)*c-v(k,i)*s 
      v(k,i)=v(k,j)*s+v(k,i)*c 
      v(k,j)=tem 
  160 continue 
      a(j,j)=v1*c*c+v3*s*s-2.*v2*s*c 
      a(i,i)=v1*s*s+v3*c*c+2.*v2*s*c 
      a(i,j)=(v1-v3)*s*c+v2*(c*c-s*s) 
  165 continue 
      if(m-1)175,170,170 
  170 m=0 
      go to 110 
  175 if(tes-rho)180,180,105 
  180 do 190 i=1,n 
      do 190 j=1,n 
  190 a(i,j)=scl*a(i,j) 
      return 
      end 

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                             c
!     Subroutine Flashsort
!     SORTS ARRAY A WITH N ELEMENTS BY USE OF INDEX VECTOR L  c
!     OF DIMENSION M WITH M ABOUT 0.1 N.                      c
!     Karl-Dietrich Neubert, FlashSort1 Algorithm             c
!     in  Dr. Dobb's Journal Feb.1998,p.123                   c
!                                                             c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine flashsort(A, N, L, M, ind)

      double precision a(*), anmin, c1, hold, flash
      integer L(*), ind(*), i, n, nmax, m, k, ihold, nmove, j, iflash
!     ============================ CLASS FORMATION ===== 


      do i = 1, n
      ind(i) = i
      end do

      ANMIN=A(1)
      NMAX=1 
      DO I=1,N
         IF( A(I).LT.ANMIN) ANMIN=A(I)
         IF( A(I).GT.A(NMAX)) NMAX=I
      END DO

      IF (ANMIN.EQ.A(NMAX)) RETURN
      C1=(M - 1) / (A(NMAX) - ANMIN)
      DO K=1,M  
         L(K)=0
      END DO 
      DO I=1,N
         K=1 + INT(C1 * (A(I) - ANMIN))
         L(K)=L(K) + 1
      END DO
      DO K=2,M
         L(K)=L(K) + L(K - 1)
      END DO
      HOLD=A(NMAX)
      A(NMAX)=A(1) 
      A(1)=HOLD

      ihold = ind(nmax)
      ind(nmax) = ind(1)
      ind(1) = ihold


!     =============================== PERMUTATION ===== 
      NMOVE=0 
      J=1
      K=M
      DO WHILE (NMOVE.LT.N - 1)
         DO WHILE (J.GT.L(K)) 
            J=J + 1 
            K=1 + INT(C1 * (A(J) - ANMIN)) 
         END DO  
         FLASH=A(J)
         iflash=ind(j)

         DO WHILE (.NOT.(J.EQ.L(K) + 1)) 
            K=1 + INT(C1 * (FLASH - ANMIN))
            HOLD=A(L(K)) 
            ihold = ind(L(k))
            A(L(K))=FLASH
            ind(L(k)) = iflash
            iflash = ihold
            FLASH=HOLD
            L(K)=L(K) - 1
            NMOVE=NMOVE + 1 
         END DO
      END DO

!     ========================= STRAIGHT INSERTION =====
      DO I=N-2,1,-1
         IF  (A(I + 1).LT.A(I)) THEN
            HOLD=A(I)
            ihold = ind(i)
            J=I
            DO WHILE  (A(J + 1).LT.HOLD)
               A(J)=A(J + 1)
               ind(j) = ind(j+1) 
               J=J + 1 
            END DO
            A(J)=HOLD 
            ind(j) = ihold
         ENDIF
      END DO

!     =========================== RETURN,END FLASH1 =====
      RETURN
      END                               

