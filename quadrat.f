c     This file contains public-code routines for numerical integration
c     of a given function
c
      SUBROUTINE QAG (F, A, B, EPSABS, EPSREL, KEY, RESULT, ABSERR,
     +   NEVAL, IER, LIMIT, LENW, LAST, IWORK, WORK)
C***BEGIN PROLOGUE  QAG
C***PURPOSE  The routine calculates an approximation result to a given
C            definite integral I = integral of F over (A,B),
C            hopefully satisfying following claim for accuracy
C            ABS(I-RESULT)LE.MAX(EPSABS,EPSREL*ABS(I)).
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A1
C***TYPE      SINGLE PRECISION (QAG-S, DQAG-D)
C***KEYWORDS  AUTOMATIC INTEGRATOR, GAUSS-KRONROD RULES,
C             GENERAL-PURPOSE, GLOBALLY ADAPTIVE, INTEGRAND EXAMINATOR,
C             QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Computation of a definite integral
C        Standard fortran subroutine
C        Real version
C
C            F      - Real
C                     Function subprogram defining the integrand
C                     Function F(X). The actual name for F needs to be
C                     Declared E X T E R N A L in the driver program.
C
C            A      - Real
C                     Lower limit of integration
C
C            B      - Real
C                     Upper limit of integration
C
C            EPSABS - Real
C                     Absolute accuracy requested
C            EPSREL - Real
C                     Relative accuracy requested
C                     If  EPSABS.LE.0
C                     And EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                     The routine will end with IER = 6.
C
C            KEY    - Integer
C                     Key for choice of local integration rule
C                     A GAUSS-KRONROD PAIR is used with
C                       7 - 15 POINTS If KEY.LT.2,
C                      10 - 21 POINTS If KEY = 2,
C                      15 - 31 POINTS If KEY = 3,
C                      20 - 41 POINTS If KEY = 4,
C                      25 - 51 POINTS If KEY = 5,
C                      30 - 61 POINTS If KEY.GT.5.
C
C         ON RETURN
C            RESULT - Real
C                     Approximation to the integral
C
C            ABSERR - Real
C                     Estimate of the modulus of the absolute error,
C                     Which should EQUAL or EXCEED ABS(I-RESULT)
C
C            NEVAL  - Integer
C                     Number of integrand evaluations
C
C            IER    - Integer
C                     IER = 0 Normal and reliable termination of the
C                             routine. It is assumed that the requested
C                             accuracy has been achieved.
C                     IER.GT.0 Abnormal termination of the routine
C                             The estimates for RESULT and ERROR are
C                             Less reliable. It is assumed that the
C                             requested accuracy has not been achieved.
C                      ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved. One can allow more
C                             subdivisions by increasing the value of
C                             LIMIT (and taking the according dimension
C                             adjustments into account). HOWEVER, If
C                             this yield no improvement it is advised
C                             to analyze the integrand in order to
C                             determine the integration difficulties.
C                             If the position of a local difficulty can
C                             be determined (I.E. SINGULARITY,
C                             DISCONTINUITY WITHIN THE INTERVAL) One
C                             will probably gain from splitting up the
C                             interval at this point and calling the
C                             INTEGRATOR on the SUBRANGES. If possible,
C                             AN APPROPRIATE SPECIAL-PURPOSE INTEGRATOR
C                             should be used which is designed for
C                             handling the type of difficulty involved.
C                         = 2 The occurrence of roundoff error is
C                             detected, which prevents the requested
C                             tolerance from being achieved.
C                         = 3 Extremely bad integrand behaviour occurs
C                             at some points of the integration
C                             interval.
C                         = 6 The input is invalid, because
C                             (EPSABS.LE.0 AND
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28))
C                             OR LIMIT.LT.1 OR LENW.LT.LIMIT*4.
C                             RESULT, ABSERR, NEVAL, LAST are set
C                             to zero.
C                             EXCEPT when LENW is invalid, IWORK(1),
C                             WORK(LIMIT*2+1) and WORK(LIMIT*3+1) are
C                             set to zero, WORK(1) is set to A and
C                             WORK(LIMIT+1) to B.
C
C         DIMENSIONING PARAMETERS
C            LIMIT - Integer
C                    Dimensioning parameter for IWORK
C                    Limit determines the maximum number of subintervals
C                    in the partition of the given integration interval
C                    (A,B), LIMIT.GE.1.
C                    If LIMIT.LT.1, the routine will end with IER = 6.
C
C            LENW  - Integer
C                    Dimensioning parameter for work
C                    LENW must be at least LIMIT*4.
C                    IF LENW.LT.LIMIT*4, the routine will end with
C                    IER = 6.
C
C            LAST  - Integer
C                    On return, LAST equals the number of subintervals
C                    produced in the subdivision process, which
C                    determines the number of significant elements
C                    actually in the WORK ARRAYS.
C
C         WORK ARRAYS
C            IWORK - Integer
C                    Vector of dimension at least limit, the first K
C                    elements of which contain pointers to the error
C                    estimates over the subintervals, such that
C                    WORK(LIMIT*3+IWORK(1)),... , WORK(LIMIT*3+IWORK(K))
C                    form a decreasing sequence with K = LAST If
C                    LAST.LE.(LIMIT/2+2), and K = LIMIT+1-LAST otherwise
C
C            WORK  - Real
C                    Vector of dimension at least LENW
C                    on return
C                    WORK(1), ..., WORK(LAST) contain the left end
C                    points of the subintervals in the partition of
C                     (A,B),
C                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain the
C                     right end points,
C                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
C                     the integral approximations over the subintervals,
C                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST) contain
C                     the error estimates.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  QAGE, XERMSG
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C***END PROLOGUE  QAG
      REAL A,ABSERR,B,EPSABS,EPSREL,F,RESULT,WORK
      INTEGER IER,IWORK,KEY,LENW,LIMIT,LVL,L1,L2,L3,NEVAL
C
      DIMENSION IWORK(*),WORK(*)
C
      EXTERNAL F
C***FIRST EXECUTABLE STATEMENT  QAG
      IER = 6
      NEVAL = 0
      LAST = 0
      RESULT = 0.0E+00
      ABSERR = 0.0E+00
      IF (LIMIT.GE.1 .AND. LENW.GE.LIMIT*4) THEN
C
C        PREPARE CALL FOR QAGE.
C
         L1 = LIMIT+1
         L2 = LIMIT+L1
         L3 = LIMIT+L2
C
         CALL QAGE(F,A,B,EPSABS,EPSREL,KEY,LIMIT,RESULT,ABSERR,NEVAL,
     1     IER,WORK(1),WORK(L1),WORK(L2),WORK(L3),IWORK,LAST)
C
C        CALL ERROR HANDLER IF NECESSARY.
C
         LVL = 0
      ENDIF
C
      IF (IER.EQ.6) LVL = 1
clks102799      IF (IER .NE. 0) CALL XERMSG ('SLATEC', 'QAG',
      IF (IER .gt. 1) CALL XERMSG ('SLATEC', 'QAG',
     +   'ABNORMAL RETURN', IER, LVL)
      RETURN
      END
*DECK QAGE
      SUBROUTINE QAGE (F, A, B, EPSABS, EPSREL, KEY, LIMIT, RESULT,
     +   ABSERR, NEVAL, IER, ALIST, BLIST, RLIST, ELIST, IORD, LAST)
C***BEGIN PROLOGUE  QAGE
C***PURPOSE  The routine calculates an approximation result to a given
C            definite integral   I = Integral of F over (A,B),
C            hopefully satisfying following claim for accuracy
C            ABS(I-RESLT).LE.MAX(EPSABS,EPSREL*ABS(I)).
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A1
C***TYPE      SINGLE PRECISION (QAGE-S, DQAGE-D)
C***KEYWORDS  AUTOMATIC INTEGRATOR, GAUSS-KRONROD RULES,
C             GENERAL-PURPOSE, GLOBALLY ADAPTIVE, INTEGRAND EXAMINATOR,
C             QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Computation of a definite integral
C        Standard fortran subroutine
C        Real version
C
C        PARAMETERS
C         ON ENTRY
C            F      - Real
C                     Function subprogram defining the integrand
C                     function F(X). The actual name for F needs to be
C                     declared E X T E R N A L in the driver program.
C
C            A      - Real
C                     Lower limit of integration
C
C            B      - Real
C                     Upper limit of integration
C
C            EPSABS - Real
C                     Absolute accuracy requested
C            EPSREL - Real
C                     Relative accuracy requested
C                     If  EPSABS.LE.0
C                     and EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                     the routine will end with IER = 6.
C
C            KEY    - Integer
C                     Key for choice of local integration rule
C                     A Gauss-Kronrod pair is used with
C                          7 - 15 points if KEY.LT.2,
C                         10 - 21 points if KEY = 2,
C                         15 - 31 points if KEY = 3,
C                         20 - 41 points if KEY = 4,
C                         25 - 51 points if KEY = 5,
C                         30 - 61 points if KEY.GT.5.
C
C            LIMIT  - Integer
C                     Gives an upper bound on the number of subintervals
C                     in the partition of (A,B), LIMIT.GE.1.
C
C         ON RETURN
C            RESULT - Real
C                     Approximation to the integral
C
C            ABSERR - Real
C                     Estimate of the modulus of the absolute error,
C                     which should equal or exceed ABS(I-RESULT)
C
C            NEVAL  - Integer
C                     Number of integrand evaluations
C
C            IER    - Integer
C                     IER = 0 Normal and reliable termination of the
C                             routine. It is assumed that the requested
C                             accuracy has been achieved.
C                     IER.GT.0 Abnormal termination of the routine
C                             The estimates for result and error are
C                             less reliable. It is assumed that the
C                             requested accuracy has not been achieved.
C            ERROR MESSAGES
C                     IER = 1 Maximum number of subdivisions allowed
C                             has been achieved. One can allow more
C                             subdivisions by increasing the value
C                             of LIMIT.
C                             However, if this yields no improvement it
C                             is rather advised to analyze the integrand
C                             in order to determine the integration
C                             difficulties. If the position of a local
C                             difficulty can be determined(e.g.
C                             SINGULARITY, DISCONTINUITY within the
C                             interval) one will probably gain from
C                             splitting up the interval at this point
C                             and calling the integrator on the
C                             subranges. If possible, an appropriate
C                             special-purpose integrator should be used
C                             which is designed for handling the type of
C                             difficulty involved.
C                         = 2 The occurrence of roundoff error is
C                             detected, which prevents the requested
C                             tolerance from being achieved.
C                         = 3 Extremely bad integrand behaviour occurs
C                             at some points of the integration
C                             interval.
C                         = 6 The input is invalid, because
C                             (EPSABS.LE.0 and
C                              EPSREL.LT.MAX(50*REL.MACH.ACC.,0.5D-28),
C                             RESULT, ABSERR, NEVAL, LAST, RLIST(1) ,
C                             ELIST(1) and IORD(1) are set to zero.
C                             ALIST(1) and BLIST(1) are set to A and B
C                             respectively.
C
C            ALIST   - Real
C                      Vector of dimension at least LIMIT, the first
C                       LAST  elements of which are the left
C                      end points of the subintervals in the partition
C                      of the given integration range (A,B)
C
C            BLIST   - Real
C                      Vector of dimension at least LIMIT, the first
C                       LAST  elements of which are the right
C                      end points of the subintervals in the partition
C                      of the given integration range (A,B)
C
C            RLIST   - Real
C                      Vector of dimension at least LIMIT, the first
C                       LAST  elements of which are the
C                      integral approximations on the subintervals
C
C            ELIST   - Real
C                      Vector of dimension at least LIMIT, the first
C                       LAST  elements of which are the moduli of the
C                      absolute error estimates on the subintervals
C
C            IORD    - Integer
C                      Vector of dimension at least LIMIT, the first K
C                      elements of which are pointers to the
C                      error estimates over the subintervals,
C                      such that ELIST(IORD(1)), ...,
C                      ELIST(IORD(K)) form a decreasing sequence,
C                      with K = LAST if LAST.LE.(LIMIT/2+2), and
C                      K = LIMIT+1-LAST otherwise
C
C            LAST    - Integer
C                      Number of subintervals actually produced in the
C                      subdivision process
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  QK15, QK21, QK31, QK41, QK51, QK61, QPSRT, R1MACH
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890831  Modified array declarations.  (WRB)
C   890831  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  QAGE
C
      REAL A,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,A2,B,BLIST,
     1  B1,B2,DEFABS,DEFAB1,DEFAB2,R1MACH,ELIST,EPMACH,
     2  EPSABS,EPSREL,ERRBND,ERRMAX,ERROR1,ERROR2,ERRO12,ERRSUM,F,
     3  RESABS,RESULT,RLIST,UFLOW
      INTEGER IER,IORD,IROFF1,IROFF2,K,KEY,KEYF,LAST,
     1  LIMIT,MAXERR,NEVAL,NRMAX
C
      DIMENSION ALIST(*),BLIST(*),ELIST(*),IORD(*),
     1  RLIST(*)
C
      EXTERNAL F
C
C            LIST OF MAJOR VARIABLES
C            -----------------------
C
C           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
C                       CONSIDERED UP TO NOW
C           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
C                      (ALIST(I),BLIST(I))
C           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
C           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST
C                       ERROR ESTIMATE
C           ERRMAX    - ELIST(MAXERR)
C           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
C           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
C           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
C                       ABS(RESULT))
C           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
C           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
C           LAST      - INDEX FOR SUBDIVISION
C
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH  IS THE LARGEST RELATIVE SPACING.
C           UFLOW  IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  QAGE
      EPMACH = R1MACH(4)
      UFLOW = R1MACH(1)
C
C           TEST ON VALIDITY OF PARAMETERS
C           ------------------------------
C
      IER = 0
      NEVAL = 0
      LAST = 0
      RESULT = 0.0E+00
      ABSERR = 0.0E+00
      ALIST(1) = A
      BLIST(1) = B
      RLIST(1) = 0.0E+00
      ELIST(1) = 0.0E+00
      IORD(1) = 0
      IF(EPSABS.LE.0.0E+00.AND.
     1  EPSREL.LT.MAX(0.5E+02*EPMACH,0.5E-14)) IER = 6
      IF(IER.EQ.6) GO TO 999
C
C           FIRST APPROXIMATION TO THE INTEGRAL
C           -----------------------------------
C
      KEYF = KEY
      IF(KEY.LE.0) KEYF = 1
      IF(KEY.GE.7) KEYF = 6
      NEVAL = 0
      IF(KEYF.EQ.1) CALL QK15(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF(KEYF.EQ.2) CALL QK21(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF(KEYF.EQ.3) CALL QK31(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF(KEYF.EQ.4) CALL QK41(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF(KEYF.EQ.5) CALL QK51(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      IF(KEYF.EQ.6) CALL QK61(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
      LAST = 1
      RLIST(1) = RESULT
      ELIST(1) = ABSERR
      IORD(1) = 1
C
C           TEST ON ACCURACY.
C
      ERRBND = MAX(EPSABS,EPSREL*ABS(RESULT))
      IF(ABSERR.LE.0.5E+02*EPMACH*DEFABS.AND.ABSERR.GT.
     1  ERRBND) IER = 2
      IF(LIMIT.EQ.1) IER = 1
      IF(IER.NE.0.OR.(ABSERR.LE.ERRBND.AND.ABSERR.NE.RESABS)
     1  .OR.ABSERR.EQ.0.0E+00) GO TO 60
C
C           INITIALIZATION
C           --------------
C
C
      ERRMAX = ABSERR
      MAXERR = 1
      AREA = RESULT
      ERRSUM = ABSERR
      NRMAX = 1
      IROFF1 = 0
      IROFF2 = 0
C
C           MAIN DO-LOOP
C           ------------
C
      DO 30 LAST = 2,LIMIT
C
C           BISECT THE SUBINTERVAL WITH THE LARGEST ERROR ESTIMATE.
C
        A1 = ALIST(MAXERR)
        B1 = 0.5E+00*(ALIST(MAXERR)+BLIST(MAXERR))
        A2 = B1
        B2 = BLIST(MAXERR)
        IF(KEYF.EQ.1) CALL QK15(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF(KEYF.EQ.2) CALL QK21(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF(KEYF.EQ.3) CALL QK31(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF(KEYF.EQ.4) CALL QK41(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF(KEYF.EQ.5) CALL QK51(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF(KEYF.EQ.6) CALL QK61(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
        IF(KEYF.EQ.1) CALL QK15(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
        IF(KEYF.EQ.2) CALL QK21(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
        IF(KEYF.EQ.3) CALL QK31(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
        IF(KEYF.EQ.4) CALL QK41(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
        IF(KEYF.EQ.5) CALL QK51(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
        IF(KEYF.EQ.6) CALL QK61(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
C
C           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
C           AND ERROR AND TEST FOR ACCURACY.
C
        NEVAL = NEVAL+1
        AREA12 = AREA1+AREA2
        ERRO12 = ERROR1+ERROR2
        ERRSUM = ERRSUM+ERRO12-ERRMAX
        AREA = AREA+AREA12-RLIST(MAXERR)
        IF(DEFAB1.EQ.ERROR1.OR.DEFAB2.EQ.ERROR2) GO TO 5
        IF(ABS(RLIST(MAXERR)-AREA12).LE.0.1E-04*ABS(AREA12)
     1  .AND.ERRO12.GE.0.99E+00*ERRMAX) IROFF1 = IROFF1+1
        IF(LAST.GT.10.AND.ERRO12.GT.ERRMAX) IROFF2 = IROFF2+1
    5   RLIST(MAXERR) = AREA1
        RLIST(LAST) = AREA2
        ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
        IF(ERRSUM.LE.ERRBND) GO TO 8
C
C           TEST FOR ROUNDOFF ERROR AND EVENTUALLY
C           SET ERROR FLAG.
C
        IF(IROFF1.GE.6.OR.IROFF2.GE.20) IER = 2
C
C           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF
C           SUBINTERVALS EQUALS LIMIT.
C
        IF(LAST.EQ.LIMIT) IER = 1
C
C           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
C           AT A POINT OF THE INTEGRATION RANGE.
C
        IF(MAX(ABS(A1),ABS(B2)).LE.(0.1E+01+0.1E+03*
     1  EPMACH)*(ABS(A2)+0.1E+04*UFLOW)) IER = 3
C
C           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
C
    8   IF(ERROR2.GT.ERROR1) GO TO 10
        ALIST(LAST) = A2
        BLIST(MAXERR) = B1
        BLIST(LAST) = B2
        ELIST(MAXERR) = ERROR1
        ELIST(LAST) = ERROR2
        GO TO 20
   10   ALIST(MAXERR) = A2
        ALIST(LAST) = A1
        BLIST(LAST) = B1
        RLIST(MAXERR) = AREA2
        RLIST(LAST) = AREA1
        ELIST(MAXERR) = ERROR2
        ELIST(LAST) = ERROR1
C
C           CALL SUBROUTINE QPSRT TO MAINTAIN THE DESCENDING ORDERING
C           IN THE LIST OF ERROR ESTIMATES AND SELECT THE
C           SUBINTERVAL WITH THE LARGEST ERROR ESTIMATE (TO BE
C           BISECTED NEXT).
C
   20   CALL QPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
C ***JUMP OUT OF DO-LOOP
        IF(IER.NE.0.OR.ERRSUM.LE.ERRBND) GO TO 40
   30 CONTINUE
C
C           COMPUTE FINAL RESULT.
C           ---------------------
C
   40 RESULT = 0.0E+00
      DO 50 K=1,LAST
        RESULT = RESULT+RLIST(K)
   50 CONTINUE
      ABSERR = ERRSUM
   60 IF(KEYF.NE.1) NEVAL = (10*KEYF+1)*(2*NEVAL+1)
      IF(KEYF.EQ.1) NEVAL = 30*NEVAL+15
  999 RETURN
      END
*DECK QK15
      SUBROUTINE QK15 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
C***BEGIN PROLOGUE  QK15
C***PURPOSE  To compute I = Integral of F over (A,B), with error
C                           estimate
C                       J = integral of ABS(F) over (A,B)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A2
C***TYPE      SINGLE PRECISION (QK15-S, DQK15-D)
C***KEYWORDS  15-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C           Integration rules
C           Standard fortran subroutine
C           Real version
C
C           PARAMETERS
C            ON ENTRY
C              F      - Real
C                       Function subprogram defining the integrand
C                       FUNCTION F(X). The actual name for F needs to be
C                       Declared E X T E R N A L in the calling program.
C
C              A      - Real
C                       Lower limit of integration
C
C              B      - Real
C                       Upper limit of integration
C
C            ON RETURN
C              RESULT - Real
C                       Approximation to the integral I
C                       Result is computed by applying the 15-POINT
C                       KRONROD RULE (RESK) obtained by optimal addition
C                       of abscissae to the 7-POINT GAUSS RULE(RESG).
C
C              ABSERR - Real
C                       Estimate of the modulus of the absolute error,
C                       which should not exceed ABS(I-RESULT)
C
C              RESABS - Real
C                       Approximation to the integral J
C
C              RESASC - Real
C                       Approximation to the integral of ABS(F-I/(B-A))
C                       over (A,B)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  QK15
C
      REAL A,ABSC,ABSERR,B,CENTR,DHLGTH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,
     1  FV1,FV2,HLGTH,RESABS,RESASC,RESG,RESK,RESKH,RESULT,R1MACH,UFLOW,
     2  WG,WGK,XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(7),FV2(7),WG(4),WGK(8),XGK(8)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 7-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 7-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE
C
      SAVE XGK, WGK, WG
      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),XGK(8)/
     1     0.9914553711208126E+00,   0.9491079123427585E+00,
     2     0.8648644233597691E+00,   0.7415311855993944E+00,
     3     0.5860872354676911E+00,   0.4058451513773972E+00,
     4     0.2077849550078985E+00,   0.0E+00              /
      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),WGK(8)/
     1     0.2293532201052922E-01,   0.6309209262997855E-01,
     2     0.1047900103222502E+00,   0.1406532597155259E+00,
     3     0.1690047266392679E+00,   0.1903505780647854E+00,
     4     0.2044329400752989E+00,   0.2094821410847278E+00/
      DATA WG(1),WG(2),WG(3),WG(4)/
     1     0.1294849661688697E+00,   0.2797053914892767E+00,
     2     0.3818300505051189E+00,   0.4179591836734694E+00/
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  QK15
      EPMACH = R1MACH(4)
      UFLOW = R1MACH(1)
C
      CENTR = 0.5E+00*(A+B)
      HLGTH = 0.5E+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      FC = F(CENTR)
      RESG = FC*WG(4)
      RESK = FC*WGK(8)
      RESABS = ABS(RESK)
      DO 10 J=1,3
        JTW = J*2
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,4
        JTWM1 = J*2-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5E+00
      RESASC = WGK(8)*ABS(FC-RESKH)
      DO 20 J=1,7
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0E+00.AND.ABSERR.NE.0.0E+00)
     1  ABSERR = RESASC*MIN(0.1E+01,
     2  (0.2E+03*ABSERR/RESASC)**1.5E+00)
      IF(RESABS.GT.UFLOW/(0.5E+02*EPMACH)) ABSERR = MAX
     1  ((EPMACH*0.5E+02)*RESABS,ABSERR)
      RETURN
      END
*DECK QK21
      SUBROUTINE QK21 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
C***BEGIN PROLOGUE  QK21
C***PURPOSE  To compute I = Integral of F over (A,B), with error
C                           estimate
C                       J = Integral of ABS(F) over (A,B)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A2
C***TYPE      SINGLE PRECISION (QK21-S, DQK21-D)
C***KEYWORDS  21-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C           Integration rules
C           Standard fortran subroutine
C           Real version
C
C           PARAMETERS
C            ON ENTRY
C              F      - Real
C                       Function subprogram defining the integrand
C                       FUNCTION F(X). The actual name for F needs to be
C                       Declared E X T E R N A L in the driver program.
C
C              A      - Real
C                       Lower limit of integration
C
C              B      - Real
C                       Upper limit of integration
C
C            ON RETURN
C              RESULT - Real
C                       Approximation to the integral I
C                       RESULT is computed by applying the 21-POINT
C                       KRONROD RULE (RESK) obtained by optimal addition
C                       of abscissae to the 10-POINT GAUSS RULE (RESG).
C
C              ABSERR - Real
C                       Estimate of the modulus of the absolute error,
C                       which should not exceed ABS(I-RESULT)
C
C              RESABS - Real
C                       Approximation to the integral J
C
C              RESASC - Real
C                       Approximation to the integral of ABS(F-I/(B-A))
C                       over (A,B)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  QK21
C
      REAL A,ABSC,ABSERR,B,CENTR,DHLGTH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,
     1  FV1,FV2,HLGTH,RESABS,RESG,RESK,RESKH,RESULT,R1MACH,UFLOW,WG,WGK,
     2  XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(10),FV2(10),WG(5),WGK(11),XGK(11)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 21-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 10-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 10-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 21-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 10-POINT GAUSS RULE
C
      SAVE XGK, WGK, WG
      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),
     1  XGK(8),XGK(9),XGK(10),XGK(11)/
     2         0.9956571630258081E+00,     0.9739065285171717E+00,
     3     0.9301574913557082E+00,     0.8650633666889845E+00,
     4     0.7808177265864169E+00,     0.6794095682990244E+00,
     5     0.5627571346686047E+00,     0.4333953941292472E+00,
     6     0.2943928627014602E+00,     0.1488743389816312E+00,
     7     0.0000000000000000E+00/
C
      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),
     1  WGK(8),WGK(9),WGK(10),WGK(11)/
     2     0.1169463886737187E-01,     0.3255816230796473E-01,
     3     0.5475589657435200E-01,     0.7503967481091995E-01,
     4     0.9312545458369761E-01,     0.1093871588022976E+00,
     5     0.1234919762620659E+00,     0.1347092173114733E+00,
     6     0.1427759385770601E+00,     0.1477391049013385E+00,
     7     0.1494455540029169E+00/
C
      DATA WG(1),WG(2),WG(3),WG(4),WG(5)/
     1     0.6667134430868814E-01,     0.1494513491505806E+00,
     2     0.2190863625159820E+00,     0.2692667193099964E+00,
     3     0.2955242247147529E+00/
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 10-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 21-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  QK21
      EPMACH = R1MACH(4)
      UFLOW = R1MACH(1)
C
      CENTR = 0.5E+00*(A+B)
      HLGTH = 0.5E+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 21-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      RESG = 0.0E+00
      FC = F(CENTR)
      RESK = WGK(11)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,5
        JTW = 2*J
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,5
        JTWM1 = 2*J-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5E+00
      RESASC = WGK(11)*ABS(FC-RESKH)
      DO 20 J=1,10
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0E+00.AND.ABSERR.NE.0.0E+00)
     1  ABSERR = RESASC*MIN(0.1E+01,
     2  (0.2E+03*ABSERR/RESASC)**1.5E+00)
      IF(RESABS.GT.UFLOW/(0.5E+02*EPMACH)) ABSERR = MAX
     1  ((EPMACH*0.5E+02)*RESABS,ABSERR)
      RETURN
      END
*DECK QK31
      SUBROUTINE QK31 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
C***BEGIN PROLOGUE  QK31
C***PURPOSE  To compute I = Integral of F over (A,B) with error
C                           estimate
C                       J = Integral of ABS(F) over (A,B)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A2
C***TYPE      SINGLE PRECISION (QK31-S, DQK31-D)
C***KEYWORDS  31-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C           Integration rules
C           Standard fortran subroutine
C           Real version
C
C           PARAMETERS
C            ON ENTRY
C              F      - Real
C                       Function subprogram defining the integrand
C                       FUNCTION F(X). The actual name for F needs to be
C                       Declared E X T E R N A L in the calling program.
C
C              A      - Real
C                       Lower limit of integration
C
C              B      - Real
C                       Upper limit of integration
C
C            ON RETURN
C              RESULT - Real
C                       Approximation to the integral I
C                       RESULT is computed by applying the 31-POINT
C                       GAUSS-KRONROD RULE (RESK), obtained by optimal
C                       addition of abscissae to the 15-POINT GAUSS
C                       RULE (RESG).
C
C              ABSERR - Real
C                       Estimate of the modulus of the modulus,
C                       which should not exceed ABS(I-RESULT)
C
C              RESABS - Real
C                       Approximation to the integral J
C
C              RESASC - Real
C                       Approximation to the integral of ABS(F-I/(B-A))
C                       over (A,B)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  QK31
      REAL A,ABSC,ABSERR,B,CENTR,DHLGTH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,
     1  FV1,FV2,HLGTH,RESABS,RESASC,RESG,RESK,RESKH,RESULT,R1MACH,UFLOW,
     2  WG,WGK,XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(15),FV2(15),XGK(16),WGK(16),WG(8)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 31-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 15-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 15-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 31-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 15-POINT GAUSS RULE
C
      SAVE XGK, WGK, WG
      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),XGK(8),
     1  XGK(9),XGK(10),XGK(11),XGK(12),XGK(13),XGK(14),XGK(15),
     2  XGK(16)/
     3     0.9980022986933971E+00,   0.9879925180204854E+00,
     4     0.9677390756791391E+00,   0.9372733924007059E+00,
     5     0.8972645323440819E+00,   0.8482065834104272E+00,
     6     0.7904185014424659E+00,   0.7244177313601700E+00,
     7     0.6509967412974170E+00,   0.5709721726085388E+00,
     8     0.4850818636402397E+00,   0.3941513470775634E+00,
     9     0.2991800071531688E+00,   0.2011940939974345E+00,
     1     0.1011420669187175E+00,   0.0E+00               /
      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),WGK(8),
     1  WGK(9),WGK(10),WGK(11),WGK(12),WGK(13),WGK(14),WGK(15),
     2  WGK(16)/
     3     0.5377479872923349E-02,   0.1500794732931612E-01,
     4     0.2546084732671532E-01,   0.3534636079137585E-01,
     5     0.4458975132476488E-01,   0.5348152469092809E-01,
     6     0.6200956780067064E-01,   0.6985412131872826E-01,
     7     0.7684968075772038E-01,   0.8308050282313302E-01,
     8     0.8856444305621177E-01,   0.9312659817082532E-01,
     9     0.9664272698362368E-01,   0.9917359872179196E-01,
     1     0.1007698455238756E+00,   0.1013300070147915E+00/
      DATA WG(1),WG(2),WG(3),WG(4),WG(5),WG(6),WG(7),WG(8)/
     1     0.3075324199611727E-01,   0.7036604748810812E-01,
     2     0.1071592204671719E+00,   0.1395706779261543E+00,
     3     0.1662692058169939E+00,   0.1861610000155622E+00,
     4     0.1984314853271116E+00,   0.2025782419255613E+00/
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 15-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 31-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  QK31
      EPMACH = R1MACH(4)
      UFLOW = R1MACH(1)
C
      CENTR = 0.5E+00*(A+B)
      HLGTH = 0.5E+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 31-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      FC = F(CENTR)
      RESG = WG(8)*FC
      RESK = WGK(16)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,7
        JTW = J*2
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,8
        JTWM1 = J*2-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5E+00
      RESASC = WGK(16)*ABS(FC-RESKH)
      DO 20 J=1,15
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0E+00.AND.ABSERR.NE.0.0E+00)
     1  ABSERR = RESASC*MIN(0.1E+01,
     2  (0.2E+03*ABSERR/RESASC)**1.5E+00)
      IF(RESABS.GT.UFLOW/(0.5E+02*EPMACH)) ABSERR = MAX
     1  ((EPMACH*0.5E+02)*RESABS,ABSERR)
      RETURN
      END
*DECK QK41
      SUBROUTINE QK41 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
C***BEGIN PROLOGUE  QK41
C***PURPOSE  To compute I = Integral of F over (A,B), with error
C                           estimate
C                       J = Integral of ABS(F) over (A,B)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A2
C***TYPE      SINGLE PRECISION (QK41-S, DQK41-D)
C***KEYWORDS  41-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C           Integration rules
C           Standard fortran subroutine
C           Real version
C
C           PARAMETERS
C            ON ENTRY
C              F      - Real
C                       Function subprogram defining the integrand
C                       FUNCTION F(X). The actual name for F needs to be
C                       declared E X T E R N A L in the calling program.
C
C              A      - Real
C                       Lower limit of integration
C
C              B      - Real
C                       Upper limit of integration
C
C            ON RETURN
C              RESULT - Real
C                       Approximation to the integral I
C                       RESULT is computed by applying the 41-POINT
C                       GAUSS-KRONROD RULE (RESK) obtained by optimal
C                       addition of abscissae to the 20-POINT GAUSS
C                       RULE (RESG).
C
C              ABSERR - Real
C                       Estimate of the modulus of the absolute error,
C                       which should not exceed ABS(I-RESULT)
C
C              RESABS - Real
C                       Approximation to the integral J
C
C              RESASC - Real
C                       Approximation to the integral of ABS(F-I/(B-A))
C                       over (A,B)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  QK41
C
      REAL A,ABSC,ABSERR,B,CENTR,DHLGTH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,
     1  FV1,FV2,HLGTH,RESABS,
     2  RESASC,RESG,RESK,RESKH,RESULT,R1MACH,UFLOW,
     3  WG,WGK,XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(20),FV2(20),XGK(21),WGK(21),WG(10)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 41-POINT GAUSS-KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 20-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 20-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 41-POINT GAUSS-KRONROD RULE
C
C           WG     - WEIGHTS OF THE 20-POINT GAUSS RULE
C
      SAVE XGK, WGK, WG
      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),XGK(8),
     1  XGK(9),XGK(10),XGK(11),XGK(12),XGK(13),XGK(14),XGK(15),
     2  XGK(16),XGK(17),XGK(18),XGK(19),XGK(20),XGK(21)/
     3     0.9988590315882777E+00,   0.9931285991850949E+00,
     4     0.9815078774502503E+00,   0.9639719272779138E+00,
     5     0.9408226338317548E+00,   0.9122344282513259E+00,
     6     0.8782768112522820E+00,   0.8391169718222188E+00,
     7     0.7950414288375512E+00,   0.7463319064601508E+00,
     8     0.6932376563347514E+00,   0.6360536807265150E+00,
     9     0.5751404468197103E+00,   0.5108670019508271E+00,
     1     0.4435931752387251E+00,   0.3737060887154196E+00,
     2     0.3016278681149130E+00,   0.2277858511416451E+00,
     3     0.1526054652409227E+00,   0.7652652113349733E-01,
     4     0.0E+00               /
      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),WGK(8),
     1  WGK(9),WGK(10),WGK(11),WGK(12),WGK(13),WGK(14),WGK(15),WGK(16),
     2  WGK(17),WGK(18),WGK(19),WGK(20),WGK(21)/
     3     0.3073583718520532E-02,   0.8600269855642942E-02,
     4     0.1462616925697125E-01,   0.2038837346126652E-01,
     5     0.2588213360495116E-01,   0.3128730677703280E-01,
     6     0.3660016975820080E-01,   0.4166887332797369E-01,
     7     0.4643482186749767E-01,   0.5094457392372869E-01,
     8     0.5519510534828599E-01,   0.5911140088063957E-01,
     9     0.6265323755478117E-01,   0.6583459713361842E-01,
     1     0.6864867292852162E-01,   0.7105442355344407E-01,
     2     0.7303069033278667E-01,   0.7458287540049919E-01,
     3     0.7570449768455667E-01,   0.7637786767208074E-01,
     4     0.7660071191799966E-01/
      DATA WG(1),WG(2),WG(3),WG(4),WG(5),WG(6),WG(7),WG(8),WG(9),WG(10)/
     1     0.1761400713915212E-01,    0.4060142980038694E-01,
     2     0.6267204833410906E-01,    0.8327674157670475E-01,
     3     0.1019301198172404E+00,    0.1181945319615184E+00,
     4     0.1316886384491766E+00,    0.1420961093183821E+00,
     5     0.1491729864726037E+00,    0.1527533871307259E+00/
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 20-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 41-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO MEAN VALUE OF F OVER (A,B), I.E.
C                    TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  QK41
      EPMACH = R1MACH(4)
      UFLOW = R1MACH(1)
C
      CENTR = 0.5E+00*(A+B)
      HLGTH = 0.5E+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 41-POINT GAUSS-KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      RESG = 0.0E+00
      FC = F(CENTR)
      RESK = WGK(21)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,10
        JTW = J*2
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,10
        JTWM1 = J*2-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5E+00
      RESASC = WGK(21)*ABS(FC-RESKH)
      DO 20 J=1,20
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0E+00.AND.ABSERR.NE.0.E+00)
     1  ABSERR = RESASC*MIN(0.1E+01,
     2  (0.2E+03*ABSERR/RESASC)**1.5E+00)
      IF(RESABS.GT.UFLOW/(0.5E+02*EPMACH)) ABSERR = MAX
     1  ((EPMACH*0.5E+02)*RESABS,ABSERR)
      RETURN
      END
*DECK QK51
      SUBROUTINE QK51 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
C***BEGIN PROLOGUE  QK51
C***PURPOSE  To compute I = Integral of F over (A,B) with error
C                           estimate
C                       J = Integral of ABS(F) over (A,B)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A2
C***TYPE      SINGLE PRECISION (QK51-S, DQK51-D)
C***KEYWORDS  51-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C           Integration rules
C           Standard fortran subroutine
C           Real version
C
C           PARAMETERS
C            ON ENTRY
C              F      - Real
C                       Function subroutine defining the integrand
C                       function F(X). The actual name for F needs to be
C                       declared E X T E R N A L in the calling program.
C
C              A      - Real
C                       Lower limit of integration
C
C              B      - Real
C                       Upper limit of integration
C
C            ON RETURN
C              RESULT - Real
C                       Approximation to the integral I
C                       RESULT is computed by applying the 51-point
C                       Kronrod rule (RESK) obtained by optimal addition
C                       of abscissae to the 25-point Gauss rule (RESG).
C
C              ABSERR - Real
C                       Estimate of the modulus of the absolute error,
C                       which should not exceed ABS(I-RESULT)
C
C              RESABS - Real
C                       Approximation to the integral J
C
C              RESASC - Real
C                       Approximation to the integral of ABS(F-I/(B-A))
C                       over (A,B)
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  QK51
C
      REAL A,ABSC,ABSERR,B,CENTR,DHLGTH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,
     1  FV1,FV2,HLGTH,RESABS,RESASC,RESG,RESK,RESKH,RESULT,R1MACH,UFLOW,
     2  WG,WGK,XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(25),FV2(25),XGK(26),WGK(26),WG(13)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
C           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
C           CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK    - ABSCISSAE OF THE 51-POINT KRONROD RULE
C                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 25-POINT
C                    GAUSS RULE
C                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
C                    ADDED TO THE 25-POINT GAUSS RULE
C
C           WGK    - WEIGHTS OF THE 51-POINT KRONROD RULE
C
C           WG     - WEIGHTS OF THE 25-POINT GAUSS RULE
C
      SAVE XGK, WGK, WG
      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),XGK(8),
     1  XGK(9),XGK(10),XGK(11),XGK(12),XGK(13),XGK(14)/
     2     0.9992621049926098E+00,   0.9955569697904981E+00,
     3     0.9880357945340772E+00,   0.9766639214595175E+00,
     4     0.9616149864258425E+00,   0.9429745712289743E+00,
     5     0.9207471152817016E+00,   0.8949919978782754E+00,
     6     0.8658470652932756E+00,   0.8334426287608340E+00,
     7     0.7978737979985001E+00,   0.7592592630373576E+00,
     8     0.7177664068130844E+00,   0.6735663684734684E+00/
       DATA XGK(15),XGK(16),XGK(17),XGK(18),XGK(19),XGK(20),XGK(21),
     1  XGK(22),XGK(23),XGK(24),XGK(25),XGK(26)/
     2     0.6268100990103174E+00,   0.5776629302412230E+00,
     3     0.5263252843347192E+00,   0.4730027314457150E+00,
     4     0.4178853821930377E+00,   0.3611723058093878E+00,
     5     0.3030895389311078E+00,   0.2438668837209884E+00,
     6     0.1837189394210489E+00,   0.1228646926107104E+00,
     7     0.6154448300568508E-01,   0.0E+00               /
      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),WGK(8),
     1  WGK(9),WGK(10),WGK(11),WGK(12),WGK(13),WGK(14)/
     2     0.1987383892330316E-02,   0.5561932135356714E-02,
     3     0.9473973386174152E-02,   0.1323622919557167E-01,
     4     0.1684781770912830E-01,   0.2043537114588284E-01,
     5     0.2400994560695322E-01,   0.2747531758785174E-01,
     6     0.3079230016738749E-01,   0.3400213027432934E-01,
     7     0.3711627148341554E-01,   0.4008382550403238E-01,
     8     0.4287284502017005E-01,   0.4550291304992179E-01/
       DATA WGK(15),WGK(16),WGK(17),WGK(18),WGK(19),WGK(20),WGK(21)
     1  ,WGK(22),WGK(23),WGK(24),WGK(25),WGK(26)/
     2     0.4798253713883671E-01,   0.5027767908071567E-01,
     3     0.5236288580640748E-01,   0.5425112988854549E-01,
     4     0.5595081122041232E-01,   0.5743711636156783E-01,
     5     0.5868968002239421E-01,   0.5972034032417406E-01,
     6     0.6053945537604586E-01,   0.6112850971705305E-01,
     7     0.6147118987142532E-01,   0.6158081806783294E-01/
      DATA WG(1),WG(2),WG(3),WG(4),WG(5),WG(6),WG(7),WG(8),WG(9),
     1  WG(10),WG(11),WG(12),WG(13)/
     2     0.1139379850102629E-01,    0.2635498661503214E-01,
     3     0.4093915670130631E-01,    0.5490469597583519E-01,
     4     0.6803833381235692E-01,    0.8014070033500102E-01,
     5     0.9102826198296365E-01,    0.1005359490670506E+00,
     6     0.1085196244742637E+00,    0.1148582591457116E+00,
     7     0.1194557635357848E+00,    0.1222424429903100E+00,
     8     0.1231760537267155E+00/
C
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 25-POINT GAUSS FORMULA
C           RESK   - RESULT OF THE 51-POINT KRONROD FORMULA
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
C                    I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  QK51
      EPMACH = R1MACH(4)
      UFLOW = R1MACH(1)
C
      CENTR = 0.5E+00*(A+B)
      HLGTH = 0.5E+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 51-POINT KRONROD APPROXIMATION TO
C           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      FC = F(CENTR)
      RESG = WG(13)*FC
      RESK = WGK(26)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,12
        JTW = J*2
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J = 1,13
        JTWM1 = J*2-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
   15 CONTINUE
      RESKH = RESK*0.5E+00
      RESASC = WGK(26)*ABS(FC-RESKH)
      DO 20 J=1,25
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0E+00.AND.ABSERR.NE.0.0E+00)
     1  ABSERR = RESASC*MIN(0.1E+01,
     2  (0.2E+03*ABSERR/RESASC)**1.5E+00)
      IF(RESABS.GT.UFLOW/(0.5E+02*EPMACH)) ABSERR = MAX
     1  ((EPMACH*0.5E+02)*RESABS,ABSERR)
      RETURN
      END
*DECK QK61
      SUBROUTINE QK61 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
C***BEGIN PROLOGUE  QK61
C***PURPOSE  To compute I = Integral of F over (A,B) with error
C                           estimate
C                       J = Integral of ABS(F) over (A,B)
C***LIBRARY   SLATEC (QUADPACK)
C***CATEGORY  H2A1A2
C***TYPE      SINGLE PRECISION (QK61-S, DQK61-D)
C***KEYWORDS  61-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
C***AUTHOR  Piessens, Robert
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C           de Doncker, Elise
C             Applied Mathematics and Programming Division
C             K. U. Leuven
C***DESCRIPTION
C
C        Integration rule
C        Standard fortran subroutine
C        Real version
C
C
C        PARAMETERS
C         ON ENTRY
C           F      - Real
C                    Function subprogram defining the integrand
C                    function F(X). The actual name for F needs to be
C                    declared E X T E R N A L in the calling program.
C
C           A      - Real
C                    Lower limit of integration
C
C           B      - Real
C                    Upper limit of integration
C
C         ON RETURN
C           RESULT - Real
C                    Approximation to the integral I
C                    RESULT is computed by applying the 61-point
C                    Kronrod rule (RESK) obtained by optimal addition of
C                    abscissae to the 30-point Gauss rule (RESG).
C
C           ABSERR - Real
C                    Estimate of the modulus of the absolute error,
C                    which should equal or exceed ABS(I-RESULT)
C
C           RESABS - Real
C                    Approximation to the integral J
C
C           RESASC - Real
C                    Approximation to the integral of ABS(F-I/(B-A))
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  R1MACH
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890531  Changed all specific intrinsics to generic.  (WRB)
C   890531  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C***END PROLOGUE  QK61
C
      REAL A,ABSC,ABSERR,B,CENTR,DHLGTH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,
     1  FV1,FV2,HLGTH,RESABS,RESASC,RESG,RESK,RESKH,RESULT,R1MACH,UFLOW,
     2  WG,WGK,XGK
      INTEGER J,JTW,JTWM1
      EXTERNAL F
C
      DIMENSION FV1(30),FV2(30),XGK(31),WGK(31),WG(15)
C
C           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE
C           INTERVAL (-1,1). BECAUSE OF SYMMETRY ONLY THE POSITIVE
C           ABSCISSAE AND THEIR CORRESPONDING WEIGHTS ARE GIVEN.
C
C           XGK   - ABSCISSAE OF THE 61-POINT KRONROD RULE
C                   XGK(2), XGK(4)  ... ABSCISSAE OF THE 30-POINT
C                   GAUSS RULE
C                   XGK(1), XGK(3)  ... OPTIMALLY ADDED ABSCISSAE
C                   TO THE 30-POINT GAUSS RULE
C
C           WGK   - WEIGHTS OF THE 61-POINT KRONROD RULE
C
C           WG    - WEIGHTS OF THE 30-POINT GAUSS RULE
C
      SAVE XGK, WGK, WG
      DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),XGK(8),
     1   XGK(9),XGK(10)/
     2     0.9994844100504906E+00,     0.9968934840746495E+00,
     3     0.9916309968704046E+00,     0.9836681232797472E+00,
     4     0.9731163225011263E+00,     0.9600218649683075E+00,
     5     0.9443744447485600E+00,     0.9262000474292743E+00,
     6     0.9055733076999078E+00,     0.8825605357920527E+00/
      DATA XGK(11),XGK(12),XGK(13),XGK(14),XGK(15),XGK(16),
     1  XGK(17),XGK(18),XGK(19),XGK(20)/
     2     0.8572052335460611E+00,     0.8295657623827684E+00,
     3     0.7997278358218391E+00,     0.7677774321048262E+00,
     4     0.7337900624532268E+00,     0.6978504947933158E+00,
     5     0.6600610641266270E+00,     0.6205261829892429E+00,
     6     0.5793452358263617E+00,     0.5366241481420199E+00/
      DATA XGK(21),XGK(22),XGK(23),XGK(24),
     1  XGK(25),XGK(26),XGK(27),XGK(28),XGK(29),XGK(30),XGK(31)/
     2     0.4924804678617786E+00,     0.4470337695380892E+00,
     3     0.4004012548303944E+00,     0.3527047255308781E+00,
     4     0.3040732022736251E+00,     0.2546369261678898E+00,
     5     0.2045251166823099E+00,     0.1538699136085835E+00,
     6     0.1028069379667370E+00,     0.5147184255531770E-01,
     7     0.0E+00                   /
      DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),WGK(8),
     1  WGK(9),WGK(10)/
     2     0.1389013698677008E-02,     0.3890461127099884E-02,
     3     0.6630703915931292E-02,     0.9273279659517763E-02,
     4     0.1182301525349634E-01,     0.1436972950704580E-01,
     5     0.1692088918905327E-01,     0.1941414119394238E-01,
     6     0.2182803582160919E-01,     0.2419116207808060E-01/
      DATA WGK(11),WGK(12),WGK(13),WGK(14),WGK(15),WGK(16),
     1  WGK(17),WGK(18),WGK(19),WGK(20)/
     2     0.2650995488233310E-01,     0.2875404876504129E-01,
     3     0.3090725756238776E-01,     0.3298144705748373E-01,
     4     0.3497933802806002E-01,     0.3688236465182123E-01,
     5     0.3867894562472759E-01,     0.4037453895153596E-01,
     6     0.4196981021516425E-01,     0.4345253970135607E-01/
      DATA WGK(21),WGK(22),WGK(23),WGK(24),
     1  WGK(25),WGK(26),WGK(27),WGK(28),WGK(29),WGK(30),WGK(31)/
     2     0.4481480013316266E-01,     0.4605923827100699E-01,
     3     0.4718554656929915E-01,     0.4818586175708713E-01,
     4     0.4905543455502978E-01,     0.4979568342707421E-01,
     5     0.5040592140278235E-01,     0.5088179589874961E-01,
     6     0.5122154784925877E-01,     0.5142612853745903E-01,
     7     0.5149472942945157E-01/
      DATA WG(1),WG(2),WG(3),WG(4),WG(5),WG(6),WG(7),WG(8)/
     1     0.7968192496166606E-02,     0.1846646831109096E-01,
     2     0.2878470788332337E-01,     0.3879919256962705E-01,
     3     0.4840267283059405E-01,     0.5749315621761907E-01,
     4     0.6597422988218050E-01,     0.7375597473770521E-01/
      DATA WG(9),WG(10),WG(11),WG(12),WG(13),WG(14),WG(15)/
     1     0.8075589522942022E-01,     0.8689978720108298E-01,
     2     0.9212252223778613E-01,     0.9636873717464426E-01,
     3     0.9959342058679527E-01,     0.1017623897484055E+00,
     4     0.1028526528935588E+00/
C
C           LIST OF MAJOR VARIABLES
C           -----------------------
C
C           CENTR  - MID POINT OF THE INTERVAL
C           HLGTH  - HALF-LENGTH OF THE INTERVAL
C           ABSC   - ABSCISSA
C           FVAL*  - FUNCTION VALUE
C           RESG   - RESULT OF THE 30-POINT GAUSS RULE
C           RESK   - RESULT OF THE 61-POINT KRONROD RULE
C           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F
C                    OVER (A,B), I.E. TO I/(B-A)
C
C           MACHINE DEPENDENT CONSTANTS
C           ---------------------------
C
C           EPMACH IS THE LARGEST RELATIVE SPACING.
C           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
C
C***FIRST EXECUTABLE STATEMENT  QK61
      EPMACH = R1MACH(4)
      UFLOW = R1MACH(1)
C
      CENTR = 0.5E+00*(B+A)
      HLGTH = 0.5E+00*(B-A)
      DHLGTH = ABS(HLGTH)
C
C           COMPUTE THE 61-POINT KRONROD APPROXIMATION TO THE
C           INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
C
      RESG = 0.0E+00
      FC = F(CENTR)
      RESK = WGK(31)*FC
      RESABS = ABS(RESK)
      DO 10 J=1,15
        JTW = J*2
        ABSC = HLGTH*XGK(JTW)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTW) = FVAL1
        FV2(JTW) = FVAL2
        FSUM = FVAL1+FVAL2
        RESG = RESG+WG(J)*FSUM
        RESK = RESK+WGK(JTW)*FSUM
        RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
      DO 15 J=1,15
        JTWM1 = J*2-1
        ABSC = HLGTH*XGK(JTWM1)
        FVAL1 = F(CENTR-ABSC)
        FVAL2 = F(CENTR+ABSC)
        FV1(JTWM1) = FVAL1
        FV2(JTWM1) = FVAL2
        FSUM = FVAL1+FVAL2
        RESK = RESK+WGK(JTWM1)*FSUM
        RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
  15    CONTINUE
      RESKH = RESK*0.5E+00
      RESASC = WGK(31)*ABS(FC-RESKH)
      DO 20 J=1,30
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
      RESULT = RESK*HLGTH
      RESABS = RESABS*DHLGTH
      RESASC = RESASC*DHLGTH
      ABSERR = ABS((RESK-RESG)*HLGTH)
      IF(RESASC.NE.0.0E+00.AND.ABSERR.NE.0.0E+00)
     1  ABSERR = RESASC*MIN(0.1E+01,
     2  (0.2E+03*ABSERR/RESASC)**1.5E+00)
      IF(RESABS.GT.UFLOW/(0.5E+02*EPMACH)) ABSERR = MAX
     1  ((EPMACH*0.5E+02)*RESABS,ABSERR)
      RETURN
      END
*DECK QPSRT
      SUBROUTINE QPSRT (LIMIT, LAST, MAXERR, ERMAX, ELIST, IORD, NRMAX)
C***BEGIN PROLOGUE  QPSRT
C***SUBSIDIARY
C***PURPOSE  Subsidiary to QAGE, QAGIE, QAGPE, QAGSE, QAWCE, QAWOE and
C            QAWSE
C***LIBRARY   SLATEC
C***TYPE      SINGLE PRECISION (QPSRT-S, DQPSRT-D)
C***KEYWORDS  SEQUENTIAL SORTING
C***AUTHOR  (UNKNOWN)
C***DESCRIPTION
C
C 1.        QPSRT
C           Ordering Routine
C              Standard FORTRAN Subroutine
C              REAL Version
C
C 2.        PURPOSE
C              This routine maintains the descending ordering
C              in the list of the local error estimates resulting from
C              the interval subdivision process. At each call two error
C              estimates are inserted using the sequential search
C              method, top-down for the largest error estimate
C              and bottom-up for the smallest error estimate.
C
C 3.        CALLING SEQUENCE
C              CALL QPSRT(LIMIT,LAST,MAXERR,ERMAX,ELIST,IORD,NRMAX)
C
C           PARAMETERS (MEANING AT OUTPUT)
C              LIMIT  - INTEGER
C                       Maximum number of error estimates the list
C                       can contain
C
C              LAST   - INTEGER
C                       Number of error estimates currently
C                       in the list
C
C              MAXERR - INTEGER
C                       MAXERR points to the NRMAX-th largest error
C                       estimate currently in the list
C
C              ERMAX  - REAL
C                       NRMAX-th largest error estimate
C                       ERMAX = ELIST(MAXERR)
C
C              ELIST  - REAL
C                       Vector of dimension LAST containing
C                       the error estimates
C
C              IORD   - INTEGER
C                       Vector of dimension LAST, the first K
C                       elements of which contain pointers
C                       to the error estimates, such that
C                       ELIST(IORD(1)),... , ELIST(IORD(K))
C                       form a decreasing sequence, with
C                       K = LAST if LAST.LE.(LIMIT/2+2), and
C                       K = LIMIT+1-LAST otherwise
C
C              NRMAX  - INTEGER
C                       MAXERR = IORD(NRMAX)
C
C***SEE ALSO  QAGE, QAGIE, QAGPE, QAGSE, QAWCE, QAWOE, QAWSE
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   800101  DATE WRITTEN
C   890831  Modified array declarations.  (WRB)
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900328  Added TYPE section.  (WRB)
C***END PROLOGUE  QPSRT
C
      REAL ELIST,ERMAX,ERRMAX,ERRMIN
      INTEGER I,IBEG,IDO,IORD,ISUCC,J,JBND,JUPBN,K,LAST,LIMIT,MAXERR,
     1  NRMAX
      DIMENSION ELIST(*),IORD(*)
C
C           CHECK WHETHER THE LIST CONTAINS MORE THAN
C           TWO ERROR ESTIMATES.
C
C***FIRST EXECUTABLE STATEMENT  QPSRT
      IF(LAST.GT.2) GO TO 10
      IORD(1) = 1
      IORD(2) = 2
      GO TO 90
C
C           THIS PART OF THE ROUTINE IS ONLY EXECUTED
C           IF, DUE TO A DIFFICULT INTEGRAND, SUBDIVISION
C           INCREASED THE ERROR ESTIMATE. IN THE NORMAL CASE
C           THE INSERT PROCEDURE SHOULD START AFTER THE
C           NRMAX-TH LARGEST ERROR ESTIMATE.
C
   10 ERRMAX = ELIST(MAXERR)
      IF(NRMAX.EQ.1) GO TO 30
      IDO = NRMAX-1
      DO 20 I = 1,IDO
        ISUCC = IORD(NRMAX-1)
C ***JUMP OUT OF DO-LOOP
        IF(ERRMAX.LE.ELIST(ISUCC)) GO TO 30
        IORD(NRMAX) = ISUCC
        NRMAX = NRMAX-1
   20    CONTINUE
C
C           COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO
C           BE MAINTAINED IN DESCENDING ORDER. THIS NUMBER
C           DEPENDS ON THE NUMBER OF SUBDIVISIONS STILL
C           ALLOWED.
C
   30 JUPBN = LAST
      IF(LAST.GT.(LIMIT/2+2)) JUPBN = LIMIT+3-LAST
      ERRMIN = ELIST(LAST)
C
C           INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN,
C           STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)).
C
      JBND = JUPBN-1
      IBEG = NRMAX+1
      IF(IBEG.GT.JBND) GO TO 50
      DO 40 I=IBEG,JBND
        ISUCC = IORD(I)
C ***JUMP OUT OF DO-LOOP
        IF(ERRMAX.GE.ELIST(ISUCC)) GO TO 60
        IORD(I-1) = ISUCC
   40 CONTINUE
   50 IORD(JBND) = MAXERR
      IORD(JUPBN) = LAST
      GO TO 90
C
C           INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP.
C
   60 IORD(I-1) = MAXERR
      K = JBND
      DO 70 J=I,JBND
        ISUCC = IORD(K)
C ***JUMP OUT OF DO-LOOP
        IF(ERRMIN.LT.ELIST(ISUCC)) GO TO 80
        IORD(K+1) = ISUCC
        K = K-1
   70 CONTINUE
      IORD(I) = LAST
      GO TO 90
   80 IORD(K+1) = LAST
C
C           SET MAXERR AND ERMAX.
C
   90 MAXERR = IORD(NRMAX)
      ERMAX = ELIST(MAXERR)
      RETURN
      END
