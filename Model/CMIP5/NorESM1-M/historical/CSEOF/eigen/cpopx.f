C *****************************************************************************
C
C       PROGRAM NAME : cpopx.f
C       PROGRAMMER : Dr. Kwang Y. Kim
C       CODE IDENTIFICATION : CPOPX/VERSION 1.0
C       CODE CLASSIFICATION : scientific computer code
C       CREATION DATE : June 19, 1998
C       REVISION DATE : not revised
C       REVISION INFORMATION : not applicable
C
C *****************************************************************************


C     This program computes principal oscillation patterns of a given
C     data satisfying 
C
C          X(t) = C X(t-1) + N(t-1)
C
C          C = A + iB       X(t) = T(t) + H(T(t))
C
C     where C is called the (complex) system matrix and H(*) denotes the 
C     Hermitian transformation.  CPOPs (complex principal oscillation
C     patterns), W = P + iQ, are the eigenvectors of the system matrix C.
C     The problem can be redefined in terms of real quantities:
C
C          D U = lambda U
C
C     where the new, real and symmetric covariance matrix D and the real
C     eigenvectors are written as
C
C          D = ( A  -B )         U = ( P )
C              ( B   A )             ( Q )
C
C     Thus the dimension of the problem is (2N x 2N).  Note that (-Q, P)
C     is also an eigenvector for an eigenvalue lambda.  They are identical 
C     up to the essential phase.
C
C     The CPOP coefficients Z are given by
C
C          Z = (W')* X
C
C     where W' is the adjoint of W, i.e.,
C
C          W* W' = (W')* W = I
C
C     Alternatively, CPOP coefficients are given by
C
C          Z(t) = lambda Z(t-1) + R(t-1)
C
C     where
C
C          R = (W')* N
C
C     Note that the system matrix is in general not symmetric.  Thus, one 
C     may expect complex eigenvalues and complex eigenfunctions.


      PARAMETER (NTMAX=600, NMAX=100, NN=2*NMAX, MST=615)
      COMMON /WORK/ Y(NTMAX), A(0:NTMAX), B(0:NTMAX)
      DIMENSION TSER(NTMAX,NMAX), THLB(NTMAX,NMAX), ARR(NTMAX),
     &          AVG(NMAX), DMTRX(NN,NN), SIG0(NN,NN), SIG1(NN,NN)
      DIMENSION COVR(NMAX,NMAX), COVI(NMAX,NMAX), WR(NMAX), WI(NMAX),
     &          ZR(NMAX,NMAX), ZI(NMAX,NMAX), ORTR(NMAX), ORTI(NMAX)
      DIMENSION SCALE(NMAX), INDX(NN)
      DIMENSION W(NMAX), V(NMAX,NMAX)
      DIMENSION EFLD(NMAX), PERD(NMAX)
      DIMENSION EOF(MST,NMAX)
      COMPLEX EGF(NMAX,NMAX), EGV(NMAX), PCTS(NTMAX), POP(MST), JIMAG
      COMPLEX TCMP(NTMAX,NMAX), EGFA(NMAX,NMAX)
      CHARACTER*50 FILNM1, FORMT1, FILNM2, FORMT2

      DATA JIMAG / (0.,1.) /

      PI = 4.0*ATAN(1.0)
      TPI = 2.0*PI
      RTD = 180./PI


5     FORMAT(A50)
      PRINT *, '  Type name and format of the input PC file.'
      READ 5, FILNM1, FORMT1
      PRINT *, '  Type name and format of the eigenfunction file.'
      READ 5, FILNM2, FORMT2
      PRINT *, '  Type the number of modes to be retained.'
      READ *, NMD
      PRINT *, '  Type the dimension (NX,NY) of sampling stations.'
      READ *, NX, NY
      NST = NX*NY
      PRINT *, '  Type the number of samples at each station.'
      READ *, NPTS
      PRINT *, '  Type the lag of the covariance matrix.'
      READ *, LTAU
      PRINT *, '  Type smoothing option (0:No, M:Moving average lag).'
      READ *, LAG
      PRINT *, '  Type POP scaling factor.'
      READ *, SCL
      PRINT *, '  Type the eigenvalue sorting option (0:No, 1:Yes).'
      READ *, ISRT
      PRINT *, '  Type the number of eigenfunctions to be printed.'
      READ *, NPRT
      PRINT *, '  Type output option.'
      PRINT *
      PRINT *, '    1: POP real and imaginary parts'
      PRINT *, '    2: POP amplitude and phase'
      PRINT *, '    4: PC time series real and imaginary parts'
      PRINT *, '    8: PC time series amplitude and phase'
      PRINT *, '   16: Adjoint POP real and imaginary parts'
      PRINT *, '   32: Adjoint POP amplitude and phase'
      PRINT *
      READ *, IOUT
      DT = 1.0*FLOAT(LTAU)

C ------- Read Input File
      IF (FORMT1.EQ.'UNF') THEN
        OPEN(UNIT=4, FILE=FILNM1, STATUS='OLD', FORM='UNFORMATTED')
        DO K=1,NMD
          READ(4)  (TSER(I,K), I=1,NPTS)
        END DO
      ELSE
        OPEN(UNIT=4, FILE=FILNM1, STATUS='OLD')
        DO K=1,NMD
          READ(4,FORMT1)  (TSER(I,K), I=1,NPTS)
        END DO
      END IF

C ------- Open Output Files
      OPEN(UNIT=7, FILE='inform.d', STATUS='UNKNOWN')
      OPEN(UNIT=8, FILE='eigen.d', STATUS='UNKNOWN')
      OPEN(UNIT=9, FILE='aigen.d', STATUS='UNKNOWN')
      OPEN(UNIT=10, FILE='pcts.d', STATUS='UNKNOWN')
      OPEN(UNIT=11, FILE='emode.d', STATUS='UNKNOWN')
      OPEN(UNIT=12, FILE='amode.d', STATUS='UNKNOWN')

C ------- Smoothing
      DO 10 J=1,NMD
        DO I=1,NPTS
          KS = MAX(I-LAG,1)
          KE = MIN(I+LAG,NPTS)
          KN = KE-KS+1
          SUM = 0.0
          DO K=KS,KE
            SUM = SUM + TSER(K,J)
          END DO
          ARR(I) = SUM/FLOAT(KN)
        END DO
        DO I=1,NPTS
          TSER(I,J) = ARR(I)
        END DO
10    CONTINUE

C ------- Remove Mean
      DO J=1,NMD
        SUM = 0.0
        DO I=1,NPTS
          SUM = SUM + TSER(I,J)
        END DO
        AVG(J) = SUM/FLOAT(NPTS)
        DO I=1,NPTS
          TSER(I,J) = TSER(I,J) - AVG(J)
        END DO
      END DO
C      WRITE(9,'(5E15.7)')  (AVG(J), J=1,NMD)

      NH = (NPTS+1)/2
      DO J=1,NMD
C ------- Fourier Analysis
        DO I=1,NPTS
          Y(I) = TSER(I,J)
        END DO
        CALL FOURIER(NPTS,0,NH)

C ------- Hilbert Transformation
        DO I=1,NPTS
          THLB(I,J) = A(0)
        END DO
        DO K=1,NH
          FRQ = TPI*FLOAT(K)/FLOAT(NPTS)
          DO I=1,NPTS
            T = FLOAT(I-1)
            THLB(I,J) = THLB(I,J) - A(K)*SIN(FRQ*T) + B(K)*COS(FRQ*T)
          END DO
        END DO

C ------- New Complex Variable
        DO J=1,NPTS
          TCMP(I,J) = TSER(I,J) + JIMAG*THLB(I,J)
        END DO
      END DO

C ------- Lagged Covariance Matrices (Complex)
      DO 20 J=1,NMD
        JJ = J+NMD
      DO 20 I=1,NMD
        II = I+NMD
        SUM1 = 0.0
        SUM2 = 0.0
        DO K=1,NPTS
          SUM1 = SUM1 + TSER(K,I)*TSER(K,J) + THLB(K,I)*THLB(K,J)
          SUM2 = SUM2 - TSER(K,I)*THLB(K,J) + THLB(K,I)*TSER(K,J)
        END DO
        SIG0(I,J) = SUM1/FLOAT(NPTS)
        SIG0(II,J) = SUM2/FLOAT(NPTS)
        SIG0(I,JJ) = -SUM2/FLOAT(NPTS)
        SIG0(II,JJ) = SUM1/FLOAT(NPTS)
        SUM1 = 0.0
        SUM2 = 0.0
        DO K=1,NPTS-LTAU
          KK = K+LTAU
          SUM1 = SUM1 + TSER(K,I)*TSER(KK,J) + THLB(K,I)*THLB(KK,J)
          SUM2 = SUM2 - TSER(K,I)*THLB(KK,J) + THLB(K,I)*TSER(KK,J)
        END DO
        SIG1(I,J) = SUM1/FLOAT(NPTS)
        SIG1(II,J) = SUM2/FLOAT(NPTS)
20    CONTINUE

C ------- System Matrix
      CALL LUDCMP(SIG0,2*NMD,NN,INDX,SIGN)
      DO I=1,NMD
        II = I+NMD
        DO J=1,2*NMD
          DMTRX(I,J) = 0.0
          DMTRX(II,J) = 0.0
        END DO
        DMTRX(I,I) = 1.0
      END DO
      DO J=1,NMD
        CALL LUBKSB(SIG0,2*NMD,NN,INDX,DMTRX(1,J))
      END DO

      DO 30 J=1,NMD
      DO 30 I=1,NMD
        II = I+NMD
        SUM1 = 0.0
        SUM2 = 0.0
        DO K=1,NMD
          KK = K+NMD
          SUM1 = SUM1 + SIG1(I,K)*DMTRX(K,J) - SIG1(II,K)*DMTRX(KK,J)
          SUM2 = SUM2 + SIG1(II,K)*DMTRX(K,J) + SIG1(I,K)*DMTRX(KK,J)
        END DO
        COVR(I,J) = SUM1
        COVI(I,J) = SUM2
30    CONTINUE

C ------- Total Variance
      TVAR = 0.0
      DO 35 I=1,NMD
        TVAR = TVAR + COVR(I,I)
35    CONTINUE
      WRITE(7,40) TVAR
40    FORMAT(5X,'TOTAL VARIANCE = ',E15.7,///)

C ------- Call Eigenvalue Routines
      CALL CBAL(NMAX,NMD,COVR,COVI,LOW,IGH,SCALE)
      CALL CORTH(NMAX,NMD,LOW,IGH,COVR,COVI,ORTR,ORTI)
      CALL COMQR2(NMAX,NMD,LOW,IGH,ORTR,ORTI,COVR,COVI,WR,WI,ZR,ZI,IERR)
      CALL CBABK2(NMAX,NMD,LOW,IGH,SCALE,NMD,ZR,ZI)

C ------- Adjoint Patterns
      DO 50 J=1,NMD
        JJ = J+NMD
      DO 50 I=1,NMD
        II = I+NMD
        SIG0(I,J) = ZR(J,I)
	SIG0(II,J) = ZI(J,I)
        SIG0(I,JJ) = -ZI(J,I)
        SIG0(II,JJ) = ZR(J,I)
50    CONTINUE

      CALL LUDCMP(SIG0,2*NMD,NN,INDX,SIGN)
      DO I=1,NMD
        II = I+NMD
        DO J=1,2*NMD
          DMTRX(I,J) = 0.0
          DMTRX(II,J) = 0.0
        END DO
        DMTRX(I,I) = 1.0
      END DO
      DO J=1,NMD
        CALL LUBKSB(SIG0,2*NMD,NN,INDX,DMTRX(1,J))
      END DO

C ------- Rearrange Eigenmodes
      NPOP = 0
      DO 55 J=1,NMD
        NPOP = NPOP + 1
        EGV(NPOP) = WR(J) + JIMAG*WI(J)
        EFLD(NPOP) = -DT/ALOG(SQRT(WR(J)**2+WI(J)**2))
        PERD(NPOP) = TPI*DT/ATAN2(ABS(WI(J)), ABS(WR(J)))
        DO I=1,NMD
          II = I+NMD
          EGF(I,NPOP) = ZR(I,J) + JIMAG*ZI(I,J)
          EGFA(I,NPOP) = DMTRX(I,J) + JIMAG*DMTRX(II,J)
        END DO
55    CONTINUE
      print *, '  3: passed'

C ------- Rotation and Normalization
      PRINT *, NPOP
      DO 60 J=1,NPOP
        SUM = 0.0
        DO I=1,NMD
          SUM = SUM + EGF(I,J)*CONJG(EGF(I,J))
        END DO
        ANORM = SQRT(SUM)
        DO I=1,NMD
          EGF(I,J) = EGF(I,J)/ANORM
          EGFA(I,J) = EGFA(I,J)*ANORM
        END DO
60    CONTINUE
      print *, '  4: passed'

C ------- Sorting by Eigenvalues
      IF (ISRT.EQ.1)  CALL EIGSRT(EGV,EGF,EGFA,EFLD,PERD,NPOP,NMD,NMAX)
      print *, '  5: passed'

C ------- Write Eigenmodes and Modal Contributions
      DO 75 I=1,MIN(NPOP,NPRT)
        WRITE(7,65) EGV(I), CABS(EGV(I)), EFLD(I), PERD(I)
65      FORMAT(5X,5E13.5,/)
        WRITE(8,70) (EGF(K,I)*SCL, K=1,NMD)
        WRITE(9,70) (EGFA(K,I)*SCL, K=1,NMD)
70      FORMAT(6E13.5)
        NMODE = I
75    CONTINUE
80    PRINT *, NMODE
      print *, '  6: passed'

C ------- Write POP Patterns
      IF ((MOD(IOUT,2).EQ.0) .AND. (MOD(IOUT/2,2).EQ.0) .AND.
     &    (MOD(IOUT/16,2).EQ.0) .AND. (MOD(IOUT/32,2).EQ.0))  GO TO 90
      IF (FORMT2.EQ.'UNF') THEN
        OPEN(UNIT=3, FILE=FILNM2, STATUS='OLD', FORM='UNFORMATTED')
        DO J=1,NMD
          READ(3)  (EOF(K,J), K=1,NST)
        END DO
      ELSE
        OPEN(UNIT=3, FILE=FILNM2, STATUS='OLD')
        DO J=1,NMD
        DO L=1,NY
          LL = (L-1)*NX
          READ(3,FORMT2)  (EOF(K+LL,J), K=1,NX)
        END DO
        END DO
      END IF
      print *, '  7: passed'

      DO I=1,MIN(NPOP,NPRT)
        DO K=1,NST
          POP(K) = (0.,0.)
        END DO
        DO 85 J=1,NMD
        DO 85 K=1,NST
          POP(K) = POP(K) + EGF(J,I)*EOF(K,J)
85      CONTINUE
        DO L=1,NY
          LL = (L-1)*NX
          WRITE(11,70)  (POP(K+LL)*SCL, K=1,NX)
        END DO
      END DO
      print *, '  8: passed'

90    CONTINUE
C ------- Write Adjoint POP Patterns
      IF ((MOD(IOUT/16,2).EQ.0) .AND. (MOD(IOUT/32,2).EQ.0))  GO TO 100
      DO I=1,MIN(NPOP,NPRT)
        DO K=1,NST
          POP(K) = (0.,0.)
        END DO
        DO 95 J=1,NMD
        DO 95 K=1,NST
          POP(K) = POP(K) + EGFA(J,I)*EOF(K,J)
95      CONTINUE
        DO L=1,NY
          LL = (L-1)*NX
          WRITE(12,70)  (POP(K+LL)*SCL, K=1,NX)
        END DO
      END DO
      print *, '  9: passed'

100   CONTINUE
C ------- PC Time Series
      IF ((MOD(IOUT/4,2).EQ.0) .AND. (MOD(IOUT/8,2).EQ.0))  GO TO 200
      DO 105 J=1,MIN(NPOP,NPRT)
        DO I=1,NPTS
          SUM1 = 0.0
          SUM2 = 0.0
          DO K=1,NMD
            SUM1 = SUM1 + TSER(I,K)*REAL(EGFA(K,J))
            SUM2 = SUM2 + TSER(I,K)*AIMAG(EGFA(K,J))
          END DO
          PCTS(I) = SUM1 + JIMAG*SUM2
        END DO

        IF (MOD(IOUT/4,2).EQ.1)  WRITE(10,70)  (PCTS(I), I=1,NPTS)
        IF (MOD(IOUT/8,2).EQ.0)  GO TO 105
        DO I=1,NPTS
          AMPL = CABS(PCTS(I))
          PHSE = ATAN2(AIMAG(PCTS(I)), REAL(PCTS(I)))*RTD
          PCTS(I) = AMPL + JIMAG*PHSE
        END DO
        WRITE(10,70)  (REAL(PCTS(I)), I=1,NPTS)
        WRITE(10,70)  (AIMAG(PCTS(I)), I=1,NPTS)
105   CONTINUE
      print *, '  10: passed'

200   STOP
      END

      SUBROUTINE CBAL(NM,N,AR,AI,LOW,IGH,SCALE)
C
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      REAL AR(NM,N),AI(NM,N),SCALE(N)
      REAL C,F,G,R,S,B2,RADIX
      LOGICAL NOCONV
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE
C     CBALANCE, WHICH IS A COMPLEX VERSION OF BALANCE,
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
C     THIS SUBROUTINE BALANCES A COMPLEX MATRIX AND ISOLATES
C     EIGENVALUES WHENEVER POSSIBLE.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX MATRIX TO BE BALANCED.
C
C     ON OUTPUT
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE BALANCED MATRIX.
C
C        LOW AND IGH ARE TWO INTEGERS SUCH THAT AR(I,J) AND AI(I,J)
C          ARE EQUAL TO ZERO IF
C           (1) I IS GREATER THAN J AND
C           (2) J=1,...,LOW-1 OR I=IGH+1,...,N.
C
C        SCALE CONTAINS INFORMATION DETERMINING THE
C           PERMUTATIONS AND SCALING FACTORS USED.
C
C     SUPPOSE THAT THE PRINCIPAL SUBMATRIX IN ROWS LOW THROUGH IGH
C     HAS BEEN BALANCED, THAT P(J) DENOTES THE INDEX INTERCHANGED
C     WITH J DURING THE PERMUTATION STEP, AND THAT THE ELEMENTS
C     OF THE DIAGONAL MATRIX USED ARE DENOTED BY D(I,J).  THEN
C        SCALE(J) = P(J),    FOR J = 1,...,LOW-1
C                 = D(J,J)       J = LOW,...,IGH
C                 = P(J)         J = IGH+1,...,N.
C     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,
C     THEN 1 TO LOW-1.
C
C     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.
C
C     THE ALGOL PROCEDURE EXC CONTAINED IN CBALANCE APPEARS IN
C     CBAL  IN LINE.  (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS
C     K,L HAVE BEEN REVERSED.)
C
C     ARITHMETIC IS REAL THROUGHOUT.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      RADIX = 16.0E0
C
      B2 = RADIX * RADIX
      K = 1
      L = N
      GO TO 100
C     .......... IN-LINE PROCEDURE FOR ROW AND
C                COLUMN EXCHANGE ..........
   20 SCALE(M) = J
      IF (J .EQ. M) GO TO 50
C
      DO 30 I = 1, L
         F = AR(I,J)
         AR(I,J) = AR(I,M)
         AR(I,M) = F
         F = AI(I,J)
         AI(I,J) = AI(I,M)
         AI(I,M) = F
   30 CONTINUE
C
      DO 40 I = K, N
         F = AR(J,I)
         AR(J,I) = AR(M,I)
         AR(M,I) = F
         F = AI(J,I)
         AI(J,I) = AI(M,I)
         AI(M,I) = F
   40 CONTINUE
C
   50 GO TO (80,130), IEXC
C     .......... SEARCH FOR ROWS ISOLATING AN EIGENVALUE
C                AND PUSH THEM DOWN ..........
   80 IF (L .EQ. 1) GO TO 280
      L = L - 1
C     .......... FOR J=L STEP -1 UNTIL 1 DO -- ..........
  100 DO 120 JJ = 1, L
         J = L + 1 - JJ
C
         DO 110 I = 1, L
            IF (I .EQ. J) GO TO 110
            IF (AR(J,I) .NE. 0.0E0 .OR. AI(J,I) .NE. 0.0E0) GO TO 120
  110    CONTINUE
C
         M = L
         IEXC = 1
         GO TO 20
  120 CONTINUE
C
      GO TO 140
C     .......... SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
C                AND PUSH THEM LEFT ..........
  130 K = K + 1
C
  140 DO 170 J = K, L
C
         DO 150 I = K, L
            IF (I .EQ. J) GO TO 150
            IF (AR(I,J) .NE. 0.0E0 .OR. AI(I,J) .NE. 0.0E0) GO TO 170
  150    CONTINUE
C
         M = K
         IEXC = 2
         GO TO 20
  170 CONTINUE
C     .......... NOW BALANCE THE SUBMATRIX IN ROWS K TO L ..........
      DO 180 I = K, L
  180 SCALE(I) = 1.0E0
C     .......... ITERATIVE LOOP FOR NORM REDUCTION ..........
  190 NOCONV = .FALSE.
C
      DO 270 I = K, L
         C = 0.0E0
         R = 0.0E0
C
         DO 200 J = K, L
            IF (J .EQ. I) GO TO 200
            C = C + ABS(AR(J,I)) + ABS(AI(J,I))
            R = R + ABS(AR(I,J)) + ABS(AI(I,J))
  200    CONTINUE
C     .......... GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW ..........
         IF (C .EQ. 0.0E0 .OR. R .EQ. 0.0E0) GO TO 270
         G = R / RADIX
         F = 1.0E0
         S = C + R
  210    IF (C .GE. G) GO TO 220
         F = F * RADIX
         C = C * B2
         GO TO 210
  220    G = R * RADIX
  230    IF (C .LT. G) GO TO 240
         F = F / RADIX
         C = C / B2
         GO TO 230
C     .......... NOW BALANCE ..........
  240    IF ((C + R) / F .GE. 0.95E0 * S) GO TO 270
         G = 1.0E0 / F
         SCALE(I) = SCALE(I) * F
         NOCONV = .TRUE.
C
         DO 250 J = K, N
            AR(I,J) = AR(I,J) * G
            AI(I,J) = AI(I,J) * G
  250    CONTINUE
C
         DO 260 J = 1, L
            AR(J,I) = AR(J,I) * F
            AI(J,I) = AI(J,I) * F
  260    CONTINUE
C
  270 CONTINUE
C
      IF (NOCONV) GO TO 190
C
  280 LOW = K
      IGH = L
      RETURN
      END

      SUBROUTINE CORTH(NM,N,LOW,IGH,AR,AI,ORTR,ORTI)
C
      INTEGER I,J,M,N,II,JJ,LA,MP,NM,IGH,KP1,LOW
      REAL AR(NM,N),AI(NM,N),ORTR(IGH),ORTI(IGH)
      REAL F,G,H,FI,FR,SCALE,PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF A COMPLEX ANALOGUE OF
C     THE ALGOL PROCEDURE ORTHES, NUM. MATH. 12, 349-368(1968)
C     BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     GIVEN A COMPLEX GENERAL MATRIX, THIS SUBROUTINE
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
C     UNITARY SIMILARITY TRANSFORMATIONS.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX INPUT MATRIX.
C
C     ON OUTPUT
C
C        AR AND AI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE HESSENBERG MATRIX.  INFORMATION
C          ABOUT THE UNITARY TRANSFORMATIONS USED IN THE REDUCTION
C          IS STORED IN THE REMAINING TRIANGLES UNDER THE
C          HESSENBERG MATRIX.
C
C        ORTR AND ORTI CONTAIN FURTHER INFORMATION ABOUT THE
C          TRANSFORMATIONS.  ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      LA = IGH - 1
      KP1 = LOW + 1
      IF (LA .LT. KP1) GO TO 200
C
      DO 180 M = KP1, LA
         H = 0.0E0
         ORTR(M) = 0.0E0
         ORTI(M) = 0.0E0
         SCALE = 0.0E0
C     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) ..........
         DO 90 I = M, IGH
   90    SCALE = SCALE + ABS(AR(I,M-1)) + ABS(AI(I,M-1))
C
         IF (SCALE .EQ. 0.0E0) GO TO 180
         MP = M + IGH
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
         DO 100 II = M, IGH
            I = MP - II
            ORTR(I) = AR(I,M-1) / SCALE
            ORTI(I) = AI(I,M-1) / SCALE
            H = H + ORTR(I) * ORTR(I) + ORTI(I) * ORTI(I)
  100    CONTINUE
C
         G = SQRT(H)
         F = PYTHAG(ORTR(M),ORTI(M))
         IF (F .EQ. 0.0E0) GO TO 103
         H = H + F * G
         G = G / F
         ORTR(M) = (1.0E0 + G) * ORTR(M)
         ORTI(M) = (1.0E0 + G) * ORTI(M)
         GO TO 105
C
  103    ORTR(M) = G
         AR(M,M-1) = SCALE
C     .......... FORM (I-(U*UT)/H) * A ..........
  105    DO 130 J = M, N
            FR = 0.0E0
            FI = 0.0E0
C     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
            DO 110 II = M, IGH
               I = MP - II
               FR = FR + ORTR(I) * AR(I,J) + ORTI(I) * AI(I,J)
               FI = FI + ORTR(I) * AI(I,J) - ORTI(I) * AR(I,J)
  110       CONTINUE
C
            FR = FR / H
            FI = FI / H
C
            DO 120 I = M, IGH
               AR(I,J) = AR(I,J) - FR * ORTR(I) + FI * ORTI(I)
               AI(I,J) = AI(I,J) - FR * ORTI(I) - FI * ORTR(I)
  120       CONTINUE
C
  130    CONTINUE
C     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) ..........
         DO 160 I = 1, IGH
            FR = 0.0E0
            FI = 0.0E0
C     .......... FOR J=IGH STEP -1 UNTIL M DO -- ..........
            DO 140 JJ = M, IGH
               J = MP - JJ
               FR = FR + ORTR(J) * AR(I,J) - ORTI(J) * AI(I,J)
               FI = FI + ORTR(J) * AI(I,J) + ORTI(J) * AR(I,J)
  140       CONTINUE
C
            FR = FR / H
            FI = FI / H
C
            DO 150 J = M, IGH
               AR(I,J) = AR(I,J) - FR * ORTR(J) - FI * ORTI(J)
               AI(I,J) = AI(I,J) + FR * ORTI(J) - FI * ORTR(J)
  150       CONTINUE
C
  160    CONTINUE
C
         ORTR(M) = SCALE * ORTR(M)
         ORTI(M) = SCALE * ORTI(M)
         AR(M,M-1) = -G * AR(M,M-1)
         AI(M,M-1) = -G * AI(M,M-1)
  180 CONTINUE
C
  200 RETURN
      END

      SUBROUTINE COMQR2(NM,N,LOW,IGH,ORTR,ORTI,HR,HI,WR,WI,ZR,ZI,IERR)
C
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,NM,NN,IGH,IP1,
     X        ITN,ITS,LOW,LP1,ENM1,IEND,IERR
      REAL HR(NM,N),HI(NM,N),WR(N),WI(N),ZR(NM,N),ZI(NM,N),
     X       ORTR(IGH),ORTI(IGH)
      REAL SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,TST1,TST2,
     X       PYTHAG
C
C     THIS SUBROUTINE IS A TRANSLATION OF A UNITARY ANALOGUE OF THE
C     ALGOL PROCEDURE  COMLR2, NUM. MATH. 16, 181-204(1970) BY PETERS
C     AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C     THE UNITARY ANALOGUE SUBSTITUTES THE QR ALGORITHM OF FRANCIS
C     (COMP. JOUR. 4, 332-345(1962)) FOR THE LR ALGORITHM.
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A COMPLEX UPPER HESSENBERG MATRIX BY THE QR
C     METHOD.  THE EIGENVECTORS OF A COMPLEX GENERAL MATRIX
C     CAN ALSO BE FOUND IF  CORTH  HAS BEEN USED TO REDUCE
C     THIS GENERAL MATRIX TO HESSENBERG FORM.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY THE BALANCING
C          SUBROUTINE  CBAL.  IF  CBAL  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        ORTR AND ORTI CONTAIN INFORMATION ABOUT THE UNITARY TRANS-
C          FORMATIONS USED IN THE REDUCTION BY  CORTH, IF PERFORMED.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.  IF THE EIGENVECTORS
C          OF THE HESSENBERG MATRIX ARE DESIRED, SET ORTR(J) AND
C          ORTI(J) TO 0.0E0 FOR THESE ELEMENTS.
C
C        HR AND HI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE COMPLEX UPPER HESSENBERG MATRIX.
C          THEIR LOWER TRIANGLES BELOW THE SUBDIAGONAL CONTAIN FURTHER
C          INFORMATION ABOUT THE TRANSFORMATIONS WHICH WERE USED IN THE
C          REDUCTION BY  CORTH, IF PERFORMED.  IF THE EIGENVECTORS OF
C          THE HESSENBERG MATRIX ARE DESIRED, THESE ELEMENTS MAY BE
C          ARBITRARY.
C
C     ON OUTPUT
C
C        ORTR, ORTI, AND THE UPPER HESSENBERG PORTIONS OF HR AND HI
C          HAVE BEEN DESTROYED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  IF AN ERROR
C          EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVECTORS.  THE EIGENVECTORS
C          ARE UNNORMALIZED.  IF AN ERROR EXIT IS MADE, NONE OF
C          THE EIGENVECTORS HAS BEEN FOUND.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C     CALLS CSROOT FOR COMPLEX SQUARE ROOT.
C     CALLS PYTHAG FOR  SQRT(A*A + B*B) .
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
C     .......... INITIALIZE EIGENVECTOR MATRIX ..........
      DO 101 J = 1, N
C
         DO 100 I = 1, N
            ZR(I,J) = 0.0E0
            ZI(I,J) = 0.0E0
  100    CONTINUE
         ZR(J,J) = 1.0E0
  101 CONTINUE
C     .......... FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS
C                FROM THE INFORMATION LEFT BY CORTH ..........
      IEND = IGH - LOW - 1
      IF (IEND) 180, 150, 105
C     .......... FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
  105 DO 140 II = 1, IEND
         I = IGH - II
         IF (ORTR(I) .EQ. 0.0E0 .AND. ORTI(I) .EQ. 0.0E0) GO TO 140
         IF (HR(I,I-1) .EQ. 0.0E0 .AND. HI(I,I-1) .EQ. 0.0E0) GO TO 140
C     .......... NORM BELOW IS NEGATIVE OF H FORMED IN CORTH ..........
         NORM = HR(I,I-1) * ORTR(I) + HI(I,I-1) * ORTI(I)
         IP1 = I + 1
C
         DO 110 K = IP1, IGH
            ORTR(K) = HR(K,I-1)
            ORTI(K) = HI(K,I-1)
  110    CONTINUE
C
         DO 130 J = I, IGH
            SR = 0.0E0
            SI = 0.0E0
C
            DO 115 K = I, IGH
               SR = SR + ORTR(K) * ZR(K,J) + ORTI(K) * ZI(K,J)
               SI = SI + ORTR(K) * ZI(K,J) - ORTI(K) * ZR(K,J)
  115       CONTINUE
C
            SR = SR / NORM
            SI = SI / NORM
C
            DO 120 K = I, IGH
               ZR(K,J) = ZR(K,J) + SR * ORTR(K) - SI * ORTI(K)
               ZI(K,J) = ZI(K,J) + SR * ORTI(K) + SI * ORTR(K)
  120       CONTINUE
C
  130    CONTINUE
C
  140 CONTINUE
C     .......... CREATE REAL SUBDIAGONAL ELEMENTS ..........
  150 L = LOW + 1
C
      DO 170 I = L, IGH
         LL = MIN0(I+1,IGH)
         IF (HI(I,I-1) .EQ. 0.0E0) GO TO 170
         NORM = PYTHAG(HR(I,I-1),HI(I,I-1))
         YR = HR(I,I-1) / NORM
         YI = HI(I,I-1) / NORM
         HR(I,I-1) = NORM
         HI(I,I-1) = 0.0E0
C
         DO 155 J = I, N
            SI = YR * HI(I,J) - YI * HR(I,J)
            HR(I,J) = YR * HR(I,J) + YI * HI(I,J)
            HI(I,J) = SI
  155    CONTINUE
C
         DO 160 J = 1, LL
            SI = YR * HI(J,I) + YI * HR(J,I)
            HR(J,I) = YR * HR(J,I) - YI * HI(J,I)
            HI(J,I) = SI
  160    CONTINUE
C
         DO 165 J = LOW, IGH
            SI = YR * ZI(J,I) + YI * ZR(J,I)
            ZR(J,I) = YR * ZR(J,I) - YI * ZI(J,I)
            ZI(J,I) = SI
  165    CONTINUE
C
  170 CONTINUE
C     .......... STORE ROOTS ISOLATED BY CBAL ..........
  180 DO 200 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 200
         WR(I) = HR(I,I)
         WI(I) = HI(I,I)
  200 CONTINUE
C
      EN = IGH
      TR = 0.0E0
      TI = 0.0E0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUE ..........
  220 IF (EN .LT. LOW) GO TO 680
      ITS = 0
      ENM1 = EN - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
  240 DO 260 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 300
         TST1 = ABS(HR(L-1,L-1)) + ABS(HI(L-1,L-1))
     X            + ABS(HR(L,L)) + ABS(HI(L,L))
         TST2 = TST1 + ABS(HR(L,L-1))
         IF (TST2 .EQ. TST1) GO TO 300
  260 CONTINUE
C     .......... FORM SHIFT ..........
  300 IF (L .EQ. EN) GO TO 660
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .EQ. 10 .OR. ITS .EQ. 20) GO TO 320
      SR = HR(EN,EN)
      SI = HI(EN,EN)
      XR = HR(ENM1,EN) * HR(EN,ENM1)
      XI = HI(ENM1,EN) * HR(EN,ENM1)
      IF (XR .EQ. 0.0E0 .AND. XI .EQ. 0.0E0) GO TO 340
      YR = (HR(ENM1,ENM1) - SR) / 2.0E0
      YI = (HI(ENM1,ENM1) - SI) / 2.0E0
      CALL CSROOT(YR**2-YI**2+XR,2.0E0*YR*YI+XI,ZZR,ZZI)
      IF (YR * ZZR + YI * ZZI .GE. 0.0E0) GO TO 310
      ZZR = -ZZR
      ZZI = -ZZI
  310 CALL CDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI)
      SR = SR - XR
      SI = SI - XI
      GO TO 340
C     .......... FORM EXCEPTIONAL SHIFT ..........
  320 SR = ABS(HR(EN,ENM1)) + ABS(HR(ENM1,EN-2))
      SI = 0.0E0
C
  340 DO 360 I = LOW, EN
         HR(I,I) = HR(I,I) - SR
         HI(I,I) = HI(I,I) - SI
  360 CONTINUE
C
      TR = TR + SR
      TI = TI + SI
      ITS = ITS + 1
      ITN = ITN - 1
C     .......... REDUCE TO TRIANGLE (ROWS) ..........
      LP1 = L + 1
C
      DO 500 I = LP1, EN
         SR = HR(I,I-1)
         HR(I,I-1) = 0.0E0
         NORM = PYTHAG(PYTHAG(HR(I-1,I-1),HI(I-1,I-1)),SR)
         XR = HR(I-1,I-1) / NORM
         WR(I-1) = XR
         XI = HI(I-1,I-1) / NORM
         WI(I-1) = XI
         HR(I-1,I-1) = NORM
         HI(I-1,I-1) = 0.0E0
         HI(I,I-1) = SR / NORM
C
         DO 490 J = I, N
            YR = HR(I-1,J)
            YI = HI(I-1,J)
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            HR(I-1,J) = XR * YR + XI * YI + HI(I,I-1) * ZZR
            HI(I-1,J) = XR * YI - XI * YR + HI(I,I-1) * ZZI
            HR(I,J) = XR * ZZR - XI * ZZI - HI(I,I-1) * YR
            HI(I,J) = XR * ZZI + XI * ZZR - HI(I,I-1) * YI
  490    CONTINUE
C
  500 CONTINUE
C
      SI = HI(EN,EN)
      IF (SI .EQ. 0.0E0) GO TO 540
      NORM = PYTHAG(HR(EN,EN),SI)
      SR = HR(EN,EN) / NORM
      SI = SI / NORM
      HR(EN,EN) = NORM
      HI(EN,EN) = 0.0E0
      IF (EN .EQ. N) GO TO 540
      IP1 = EN + 1
C
      DO 520 J = IP1, N
         YR = HR(EN,J)
         YI = HI(EN,J)
         HR(EN,J) = SR * YR + SI * YI
         HI(EN,J) = SR * YI - SI * YR
  520 CONTINUE
C     .......... INVERSE OPERATION (COLUMNS) ..........
  540 DO 600 J = LP1, EN
         XR = WR(J-1)
         XI = WI(J-1)
C
         DO 580 I = 1, J
            YR = HR(I,J-1)
            YI = 0.0E0
            ZZR = HR(I,J)
            ZZI = HI(I,J)
            IF (I .EQ. J) GO TO 560
            YI = HI(I,J-1)
            HI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
  560       HR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            HR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            HI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  580    CONTINUE
C
         DO 590 I = LOW, IGH
            YR = ZR(I,J-1)
            YI = ZI(I,J-1)
            ZZR = ZR(I,J)
            ZZI = ZI(I,J)
            ZR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
            ZI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
            ZR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
            ZI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  590    CONTINUE
C
  600 CONTINUE
C
      IF (SI .EQ. 0.0E0) GO TO 240
C
      DO 630 I = 1, EN
         YR = HR(I,EN)
         YI = HI(I,EN)
         HR(I,EN) = SR * YR - SI * YI
         HI(I,EN) = SR * YI + SI * YR
  630 CONTINUE
C
      DO 640 I = LOW, IGH
         YR = ZR(I,EN)
         YI = ZI(I,EN)
         ZR(I,EN) = SR * YR - SI * YI
         ZI(I,EN) = SR * YI + SI * YR
  640 CONTINUE
C
      GO TO 240
C     .......... A ROOT FOUND ..........
  660 HR(EN,EN) = HR(EN,EN) + TR
      WR(EN) = HR(EN,EN)
      HI(EN,EN) = HI(EN,EN) + TI
      WI(EN) = HI(EN,EN)
      EN = ENM1
      GO TO 220
C     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
C                VECTORS OF UPPER TRIANGULAR FORM ..........
  680 NORM = 0.0E0
C
      DO 720 I = 1, N
C
         DO 720 J = I, N
            TR = ABS(HR(I,J)) + ABS(HI(I,J))
            IF (TR .GT. NORM) NORM = TR
  720 CONTINUE
C
      IF (N .EQ. 1 .OR. NORM .EQ. 0.0E0) GO TO 1001
C     .......... FOR EN=N STEP -1 UNTIL 2 DO -- ..........
      DO 800 NN = 2, N
         EN = N + 2 - NN
         XR = WR(EN)
         XI = WI(EN)
         HR(EN,EN) = 1.0E0
         HI(EN,EN) = 0.0E0
         ENM1 = EN - 1
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
         DO 780 II = 1, ENM1
            I = EN - II
            ZZR = 0.0E0
            ZZI = 0.0E0
            IP1 = I + 1
C
            DO 740 J = IP1, EN
               ZZR = ZZR + HR(I,J) * HR(J,EN) - HI(I,J) * HI(J,EN)
               ZZI = ZZI + HR(I,J) * HI(J,EN) + HI(I,J) * HR(J,EN)
  740       CONTINUE
C
            YR = XR - WR(I)
            YI = XI - WI(I)
            IF (YR .NE. 0.0E0 .OR. YI .NE. 0.0E0) GO TO 765
               TST1 = NORM
               YR = TST1
  760          YR = 0.01E0 * YR
               TST2 = NORM + YR
               IF (TST2 .GT. TST1) GO TO 760
  765       CONTINUE
            CALL CDIV(ZZR,ZZI,YR,YI,HR(I,EN),HI(I,EN))
C     .......... OVERFLOW CONTROL ..........
            TR = ABS(HR(I,EN)) + ABS(HI(I,EN))
            IF (TR .EQ. 0.0E0) GO TO 780
            TST1 = TR
            TST2 = TST1 + 1.0E0/TST1
            IF (TST2 .GT. TST1) GO TO 780
            DO 770 J = I, EN
               HR(J,EN) = HR(J,EN)/TR
               HI(J,EN) = HI(J,EN)/TR
  770       CONTINUE
C
  780    CONTINUE
C
  800 CONTINUE
C     .......... END BACKSUBSTITUTION ..........
      ENM1 = N - 1
C     .......... VECTORS OF ISOLATED ROOTS ..........
      DO  840 I = 1, ENM1
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840
         IP1 = I + 1
C
         DO 820 J = IP1, N
            ZR(I,J) = HR(I,J)
            ZI(I,J) = HI(I,J)
  820    CONTINUE
C
  840 CONTINUE
C     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                VECTORS OF ORIGINAL FULL MATRIX.
C                FOR J=N STEP -1 UNTIL LOW+1 DO -- ..........
      DO 880 JJ = LOW, ENM1
         J = N + LOW - JJ
         M = MIN0(J,IGH)
C
         DO 880 I = LOW, IGH
            ZZR = 0.0E0
            ZZI = 0.0E0
C
            DO 860 K = LOW, M
               ZZR = ZZR + ZR(I,K) * HR(K,J) - ZI(I,K) * HI(K,J)
               ZZI = ZZI + ZR(I,K) * HI(K,J) + ZI(I,K) * HR(K,J)
  860       CONTINUE
C
            ZR(I,J) = ZZR
            ZI(I,J) = ZZI
  880 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
      END

      SUBROUTINE CBABK2(NM,N,LOW,IGH,SCALE,M,ZR,ZI)
C
      INTEGER I,J,K,M,N,II,NM,IGH,LOW
      REAL SCALE(N),ZR(NM,M),ZI(NM,M)
      REAL S
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE
C     CBABK2, WHICH IS A COMPLEX VERSION OF BALBAK,
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A COMPLEX GENERAL
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     BALANCED MATRIX DETERMINED BY  CBAL.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY  CBAL.
C
C        SCALE CONTAINS INFORMATION DETERMINING THE PERMUTATIONS
C          AND SCALING FACTORS USED BY  CBAL.
C
C        M IS THE NUMBER OF EIGENVECTORS TO BE BACK TRANSFORMED.
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVECTORS TO BE
C          BACK TRANSFORMED IN THEIR FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        ZR AND ZI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE TRANSFORMED EIGENVECTORS
C          IN THEIR FIRST M COLUMNS.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IF (M .EQ. 0) GO TO 200
      IF (IGH .EQ. LOW) GO TO 120
C
      DO 110 I = LOW, IGH
         S = SCALE(I)
C     .......... LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
C                IF THE FOREGOING STATEMENT IS REPLACED BY
C                S=1.0E0/SCALE(I). ..........
         DO 100 J = 1, M
            ZR(I,J) = ZR(I,J) * S
            ZI(I,J) = ZI(I,J) * S
  100    CONTINUE
C
  110 CONTINUE
C     .......... FOR I=LOW-1 STEP -1 UNTIL 1,
C                IGH+1 STEP 1 UNTIL N DO -- ..........
  120 DO 140 II = 1, N
         I = II
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 140
         IF (I .LT. LOW) I = LOW - II
         K = SCALE(I)
         IF (K .EQ. I) GO TO 140
C
         DO 130 J = 1, M
            S = ZR(I,J)
            ZR(I,J) = ZR(K,J)
            ZR(K,J) = S
            S = ZI(I,J)
            ZI(I,J) = ZI(K,J)
            ZI(K,J) = S
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END

      SUBROUTINE EIGSRT(D,V,VA,ARR,BRR,N,M,NP)

      COMPLEX D(NP),V(NP,NP),VA(NP,NP),P
      DIMENSION ARR(NP),BRR(NP)

      DO 13 I=1,N-1
        K=I
        P=D(I)
        DO 11 J=I+1,N
          IF(CABS(D(J)).GE.CABS(P))THEN
            K=J
            P=D(J)
          ENDIF
11      CONTINUE
        IF(K.NE.I)THEN
          D(K)=D(I)
          D(I)=P
          DO 12 J=1,M
            P=V(J,I)
            V(J,I)=V(J,K)
            V(J,K)=P
            P=VA(J,I)
            VA(J,I)=VA(J,K)
            VA(J,K)=P
12        CONTINUE
          TMP=ARR(I)
          ARR(I)=ARR(K)
          ARR(K)=TMP
          TMP=BRR(I)
          BRR(I)=BRR(K)
          BRR(K)=TMP
        ENDIF
13    CONTINUE

      RETURN
      END

      SUBROUTINE LUDCMP(A,N,NP,INDX,D)

C       Given an N x N matrix A, with physical dimension NP, this routine
C       replaces it by the LU decomposition of a rowwise permutation of
C       itself.  A and N are input.  A is output; INDX is an output vector 
C       which records the row permutation effected by the partial pivoting;
C       D is output as +/-1 depending on whether the number of row 
C       interchanges was even or odd, respectively.  This routine is used
C       in combination with LUBKSBR to solve linear equations or invert a 
C       matrix.

      PARAMETER (NMAX=100, TINY=1.0E-20)
      DIMENSION INDX(N), VV(NMAX), A(NP,NP)


        D = 1.
C ---------- Find maximum entry in each row -- scaling factor
      DO 10 I=1,N
          AAMAX = 0.
        DO 5 J=1,N
          IF (ABS(A(I,J)) .GT. AAMAX)  AAMAX = ABS(A(I,J))
5       CONTINUE
        IF (AAMAX .EQ. 0.)  STOP
        VV(I) = 1./AAMAX
10    CONTINUE

C ---------- LU decomposition
      DO 50 J=1,N
        IF (J .GT. 1)  THEN
          DO 20 I=1,J-1
            SUM = A(I,J)
            IF (I .GT. 1)  THEN
              DO 15 K=1,I-1
                SUM = SUM - A(I,K)*A(K,J)
15            CONTINUE
              A(I,J) = SUM
            END IF
20        CONTINUE
        END IF

          AAMAX = 0.
        DO 30 I=J,N
          SUM = A(I,J)
          IF (J .GT. 1)  THEN
            DO 25 K=1,J-1
              SUM = SUM - A(I,K)*A(K,J)
25          CONTINUE
            A(I,J) = SUM
          END IF
          DUM = VV(I)*ABS(SUM)
          IF (DUM .GE. AAMAX)  THEN
            IMAX = I
            AAMAX = DUM
          END IF
30      CONTINUE

        IF (J .NE. IMAX)  THEN
          DO 35 K=1,N
            CUM = A(IMAX,K)
            A(IMAX,K) = A(J,K)
            A(J,K) = CUM
35        CONTINUE
          D = -D
          VV(IMAX) = VV(J)
        END IF
        INDX(J) = IMAX
        IF (J .NE. N)  THEN
          IF (A(J,J) .EQ. 0.)  A(J,J) = TINY
          CUM = 1./A(J,J)
          DO 40 I=J+1,N
            A(I,J) = A(I,J)*CUM
40        CONTINUE
        END IF
50    CONTINUE
      IF (A(N,N) .EQ. 0.)  A(N,N) = TINY

      RETURN
      END

      SUBROUTINE LUBKSB(A,N,NP,INDX,B)

C       Solves the set of N linear equations AX = B.  Here A is input, not
C       as the matrix A but as its LU decomposition, determined by the routine 
C       LUDCMP.  INDX is input as the permutation vector returned by LUDCMP.
C       B is input as the right-hand side vector B, and returns with the 
C       solution vector X.  A, N, NP and INDX are not modified by this routine
C       and can be left in place for successive calls with different right-
C       hand sides B.  This routine takes into account the possibility that 
C       B will begin with many zero elements, so it is efficient for use in
C       matrix inversion.

      DIMENSION INDX(N), A(NP,NP), B(N)


C ---------- Ly = b
        II = 0
      DO 10 I=1,N
        LL = INDX(I)
        SUM = B(LL)
        B(LL) = B(I)
        IF (II .NE. 0)  THEN
          DO 5 J=II,I-1
            SUM = SUM - A(I,J)*B(J)
5         CONTINUE
        ELSE IF (SUM .NE. 0.)  THEN
          II = I
        END IF
        B(I) = SUM
10    CONTINUE

C ---------- Ux = y
      DO 20 I=N,1,-1
        SUM = B(I)
        IF (I .LT. N)  THEN
          DO 15 J=I+1,N
            SUM = SUM - A(I,J)*B(J)
15        CONTINUE
        END IF
        B(I) = SUM/A(I,I)
20    CONTINUE

      RETURN
      END

      SUBROUTINE FOURIER(NFPTS,NHS,NHE)

      PARAMETER (NTMAX=600)
      COMMON /WORK/ Y(NTMAX), A(0:NTMAX), B(0:NTMAX)

      DO 15 NFP=NHS,NHE
        COEF = 2.0/FLOAT(NFPTS)
        CONST = 3.1415926536*COEF*NFP
        U1 = 0.0
        U2 = 0.0
        C = COS(CONST)
        S = SIN(CONST)
        I = NFPTS

5       U0 = Y(I) + 2.0*C*U1 - U2
        U2 = U1
        U1 = U0
        I = I - 1
        IF (I-1) 10, 10, 5

10      A(NFP) = COEF*(Y(1) + C*U1 - U2)
        B(NFP) = COEF*S*U1
        IF (MOD(2*NFP,NFPTS).EQ.0)  A(NFP) = 0.5*A(NFP)
15    CONTINUE

      RETURN
      END

      SUBROUTINE CDIV(AR,AI,BR,BI,CR,CI)
      REAL AR,AI,BR,BI,CR,CI
C
C     COMPLEX DIVISION, (CR,CI) = (AR,AI)/(BR,BI)
C
      REAL S,ARS,AIS,BRS,BIS
      S = ABS(BR) + ABS(BI)
      ARS = AR/S
      AIS = AI/S
      BRS = BR/S
      BIS = BI/S
      S = BRS**2 + BIS**2
      CR = (ARS*BRS + AIS*BIS)/S
      CI = (AIS*BRS - ARS*BIS)/S
      RETURN
      END

      SUBROUTINE CSROOT(XR,XI,YR,YI)
      REAL XR,XI,YR,YI
C
C     (YR,YI) = COMPLEX SQRT(XR,XI)
C     BRANCH CHOSEN SO THAT YR .GE. 0.0 AND SIGN(YI) .EQ. SIGN(XI)
C
      REAL S,TR,TI,PYTHAG
      TR = XR
      TI = XI
      S = SQRT(0.5E0*(PYTHAG(TR,TI) + ABS(TR)))
      IF (TR .GE. 0.0E0) YR = S
      IF (TI .LT. 0.0E0) S = -S
      IF (TR .LE. 0.0E0) YI = S
      IF (TR .LT. 0.0E0) YR = 0.5E0*(TI/YI)
      IF (TR .GT. 0.0E0) YI = 0.5E0*(TI/YR)
      RETURN
      END

      REAL FUNCTION PYTHAG(A,B)
      REAL A,B
C
C     FINDS SQRT(A**2+B**2) WITHOUT OVERFLOW OR DESTRUCTIVE UNDERFLOW
C
      REAL P,R,S,T,U
      P = AMAX1(ABS(A),ABS(B))
      IF (P .EQ. 0.0E0) GO TO 20
      R = (AMIN1(ABS(A),ABS(B))/P)**2
   10 CONTINUE
         T = 4.0E0 + R
         IF (T .EQ. 4.0E0) GO TO 20
         S = R/T
         U = 1.0E0 + 2.0E0*S
         P = U*P
         R = (S/U)**2 * R
      GO TO 10
   20 PYTHAG = P
      RETURN
      END
