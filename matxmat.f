      SUBROUTINE MATxMAT(A,B,L,M,N,IA,IB,C,IC)

C     core routine on file matxmat.f (canned routine)
C
c     called by INFBOTM, RADW, and RADZETA
C
C     Multiply matrix A times matrix B and return the result in matrix C.
C     A is L by M, B is M by N, and C is L by N.
C     IA, IB, and IC give the row dimensions of A, B, and C in
C     the calling program.
C
      DIMENSION A(IA,M),B(IB,N),C(IC,N)
C               
      DO i=1,L
         DO j=1,N
           SUM = 0.0
           DO k=1,M
                SUM = SUM + A(i,k)*B(k,j)
           ENDDO
           C(i,j)=SUM
         ENDDO
      ENDDO
C
      RETURN
      END
