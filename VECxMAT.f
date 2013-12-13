      SUBROUTINE VECxMAT(A,B,M,N,IB,C)
C
C     core routine on file VECxMAT.f
C
c     called by INFBOTM, RADW, and RADZETA
C
C     Multiply row vector A times matrix B and return the result
c     in row vector C.
C     A is 1 by M, B is M by N, and C is 1 by N.
C     IB gives the row dimension of B in the calling program.
C
      DIMENSION A(M),B(IB,N),C(N)
C               
         DO j=1,N
           SUM = 0.0
           DO k=1,M
                SUM = SUM + A(k)*B(k,j)
           ENDDO
           C(j)=SUM
         ENDDO
C
      RETURN
      END
