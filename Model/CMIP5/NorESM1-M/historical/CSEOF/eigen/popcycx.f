C *****************************************************************************
C
C       PROGRAM NAME : popcycx.f
C       PROGRAMMER : Dr. Kwang Y. Kim
C       CODE IDENTIFICATION : POPCYCX/VERSION 1.0
C       CODE CLASSIFICATION : scientific computer code
C       CREATION DATE : August 22, 1997
C       REVISION DATE : not revised
C       REVISION INFORMATION : not applicable
C
C *****************************************************************************


C     This program computes principal oscillation patterns of a given
C     data satisfying 
C
C          T(t,t'+1) = A(t') T(t,t') + N(t,t')
C
C     where A is the system matrix.  The time series T(t,t') denotes a
C     set of principal component (PC) time series of a regular EOF
C     analysis.  The dimension of the problem is determined by the
C     truncation order.  Because of the imposed periodicity there are 
C     n models of this form, namely for different t'.  Then
C
C          B(t) = A(t+n-1) A(t+n-2) ... A(t+1) A(t)
C
C     so that
C
C          T(t+1,t') = B(t') T(t,t') + N(t,t')
C
C     Cyclostationary POPs (principal oscillation patterns), P, are the 
C     eigenvectors of the matrix B(t).  Note that there are n different 
C     models for which eigenvalues are the same.  If p is the eigenvector
C     of B(t), then, A(t) p is the eigenvector of B(t+1).
C
C     The POP coefficients Z are given by
C
C          Z = (P')* T
C
C     where P' is the adjoint of P, i.e.,
C
C          P* P' = (P')* P = I
C
C     Alternatively, POP coefficients are given by
C
C          Z(t) = lambda Z(t-1) + R(t-1)
C
C     where
C
C          R = (P')* N
C
C
C     Note that the system matrix is in general not symmetric.  Thus, one 
C     may expect complex eigenvalues and complex eigenfunctions.


      PARAMETER (NTMAX=600, NMAX=100, MCYC=12, MST=615)
      DIMENSION TSER(NTMAX,NMAX), ARR(NTMAX), AVG(NMAX,MCYC),
     &          DMTRX(NMAX,NMAX), SIG0(NMAX,NMAX), SIG1(NMAX,NMAX)
      DIMENSION COV(NMAX,NMAX,MCYC), WR(NMAX), WI(NMAX), Z(NMAX,NMAX)
      DIMENSION SCALE(NMAX), INDX(NMAX)
      DIMENSION W(NMAX), V(NMAX,NMAX)
      DIMENSION EFLD(NMAX), PERD(NMAX)
      DIMENSION EOF(MST,NMAX)
      COMPLEX EGF(NMAX,NMAX), EGV(NMAX), FCT(NMAX), PCTS(NTMAX,NMAX),
     &        EGFA(NMAX,NMAX), FCTA(NMAX), POP(MST), JIMAG, CSUM
      CHARACTER*50 FILNM1, FORMT1, FILNM2, FORMT2

      DATA JIMAG / (0.,1.) /

      PI = 4.0*ATAN(1.0)
      TPI = 2.0*PI
      RTD = 180./PI
      DT = 1.0


5     FORMAT(A50)
      PRINT *, '  Type name and format of the PC input file.'
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
      PRINT *, '  Type the nested periodicity of the time series.'
      READ *, NCYC
      PRINT *, '  Type smoothing option (0:No, M:Moving average lag).'
      READ *, LAG
      PRINT *, '  Type detrending option (0:Mean, 1:Ann Cyc).'
      READ *, IDTR
      PRINT *, '  Type adjoint computation option (0:No, 1:Yes).'
      READ *, IADJ
      PRINT *, '  Type EOF scaling factor.'
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
      PRINT *
      READ *, IOUT

      DT = DT*FLOAT(NCYC)

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
      OPEN(UNIT=9, FILE='avg.d', STATUS='UNKNOWN')
      OPEN(UNIT=10, FILE='pcts.d', STATUS='UNKNOWN')
      OPEN(UNIT=11, FILE='emode.d', STATUS='UNKNOWN')

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
      NPTS = INT(NPTS/NCYC)*NCYC
      PRINT *, '  Number of points = ', NPTS
      LCYC = NCYC
      IF (IDTR.EQ.0)  LCYC=1
      DO 15 J=1,NMD
      DO 15 IM=1,LCYC
        SUM = 0.0
        DO I=IM,NPTS,LCYC
          SUM = SUM + TSER(I,J)
        END DO
        AVG(J,IM) = SUM/FLOAT(NPTS/LCYC)
        DO I=IM,NPTS,LCYC
          TSER(I,J) = TSER(I,J) - AVG(J,IM)
        END DO
15    CONTINUE
      DO IM=1,LCYC
        WRITE(9,'(5E15.7)')  (AVG(J,IM), J=1,NMD)
      END DO

      DO IM=1,NCYC
        PRINT *, 'Cycle = ', IM
C ------- Lagged Covariance Matrices
        DO 20 J=1,NMD
        DO 20 I=1,NMD
          SUM = 0.0
          DO K=IM,NPTS,NCYC
            SUM = SUM + TSER(K,I)*TSER(K,J)
          END DO
          SIG0(I,J) = SUM/FLOAT(NPTS/NCYC)
          SUM = 0.0
          DO K=IM,NPTS-NCYC,NCYC
            SUM = SUM + TSER(K+1,I)*TSER(K,J)
          END DO
          SIG1(I,J) = SUM/FLOAT(NPTS/NCYC)
20      CONTINUE

C ------- System Matrix
        CALL LUDCMP(SIG0,NMD,NMAX,INDX,SIGN)
        DO I=1,NMD
          DO J=1,NMD
            DMTRX(I,J) = 0.0
          END DO
          DMTRX(I,I) = 1.0
        END DO
        DO J=1,NMD
          CALL LUBKSB(SIG0,NMD,NMAX,INDX,DMTRX(1,J))
        END DO

        DO 30 J=1,NMD
        DO 30 I=1,NMD
          SUM = 0.0
          DO K=1,NMD
            SUM = SUM + SIG1(I,K)*DMTRX(K,J)
          END DO
          COV(I,J,IM) = SUM
30      CONTINUE
      END DO

C ------- Cyclostationary Covariance Matrix
      DO 35 J=1,NMD
      DO 35 I=1,NMD
        DMTRX(I,J) = COV(I,J,1)
35    CONTINUE

      DO IM=2,NCYC
        DO 40 J=1,NMD
        DO 40 I=1,NMD
          SUM = 0.0
          DO K=1,NMD
            SUM = SUM + COV(I,K,IM)*DMTRX(K,J)
          END DO
          SIG0(I,J) = SUM
40      CONTINUE
        DO 45 J=1,NMD
        DO 45 I=1,NMD
          DMTRX(I,J) = SIG0(I,J)
45      CONTINUE
      END DO

C ------- Total Variance
      TVAR = 0.0
      DO I=1,NMD
        TVAR = TVAR + DMTRX(I,I)
      END DO
      WRITE(7,50) TVAR
50    FORMAT(5X,'TOTAL VARIANCE = ',E15.7,///)

C ------- Call Eigenvalue Routines
      CALL BALANC(NMAX,NMD,DMTRX,LOW,IGH,SCALE)
      CALL ELMHES(NMAX,NMD,LOW,IGH,DMTRX,INDX)
      CALL ELTRAN(NMAX,NMD,LOW,IGH,DMTRX,INDX,Z)
      CALL HQR2(NMAX,NMD,LOW,IGH,DMTRX,WR,WI,Z,IERR)
      CALL BALBAK(NMAX,NMD,LOW,IGH,SCALE,NMD,Z)

C ------- Adjoint Patterns
      IF (IADJ.EQ.1) THEN
        DO 55 J=1,NMD
        DO 55 I=1,NMD
          SIG0(I,J) = Z(J,I)
55      CONTINUE

        CALL LUDCMP(SIG0,NMD,NMAX,INDX,SIGN)
        DO I=1,NMD
          DO J=1,NMD
            DMTRX(I,J) = 0.0
          END DO
          DMTRX(I,I) = 1.0
        END DO
        DO J=1,NMD
          CALL LUBKSB(SIG0,NMD,NMAX,INDX,DMTRX(1,J))
        END DO
      END IF

C ------- Rearrange Eigenmodes
      JCON = 0
      JPOP = 0
      DO 60 J=1,NMD
        IF (JCON.EQ.1) THEN
          JCON = 0
          GO TO 60
        END IF

        JPOP = JPOP + 1
        EGV(JPOP) = WR(J) + JIMAG*WI(J)
        EFLD(JPOP) = -DT/ALOG(SQRT(WR(J)**2+WI(J)**2))
        IF (WI(J).EQ.0.) THEN
          PERD(JPOP) = 0.0
          DO I=1,NMD
            EGF(I,JPOP) = Z(I,J)
            EGFA(I,JPOP) = DMTRX(I,J)
          END DO
        ELSE
          PERD(JPOP) = TPI*DT/ATAN2(ABS(WI(J)), ABS(WR(J)))
          DO I=1,NMD
            EGF(I,JPOP) = Z(I,J) + JIMAG*Z(I,J+1)
            EGFA(I,JPOP) = 0.5*(DMTRX(I,J) + JIMAG*DMTRX(I,J+1))
          END DO
          JCON = 1
        END IF
C        IF (CABS(EGV(JPOP)).GE.1.0)  JPOP = JPOP-1
60    CONTINUE

C ------- Rotation and Normalization
      NPOP = JPOP
      PRINT *, NPOP
      DO 70 J=1,NPOP
C ----------> rotate such that the first mode is real
        CSUM = EGF(1,J)
        DO I=1,NMD
          EGF(I,J) = EGF(I,J)/CSUM
          EGFA(I,J) = EGFA(I,J)/CONJG(CSUM)
        END DO
C        SUM1 = 0.0
C        SUM2 = 0.0
C        SUMC = 0.0
C        DO I=1,NMD
C          SUM1 = SUM1 + REAL(EGF(I,J))**2
C          SUM2 = SUM2 + AIMAG(EGF(I,J))**2
C          SUMC = SUMC + REAL(EGF(I,J))*AIMAG(EGF(I,J))
C        END DO
C        IF (SUMC.EQ.0.) THEN
C          COSH = 1.0
C          SINH = 0.0
C        ELSE
C          COT2H = (SUM1 - SUM2)/(2.*SUMC)
C          TANH = COT2H + SQRT(COT2H*COT2H + 1.)
C          COSH = SQRT(1./(1.+TANH*TANH))
C          SINH = TANH*COSH
C        END IF
C
C        SUM1 = 0.0
C        SUM2 = 0.0
C        DO I=1,NMD
C          TMPR = REAL(EGF(I,J))
C          TMPI = AIMAG(EGF(I,J))
C          VALR = COSH*TMPR - SINH*TMPI
C          VALI = SINH*TMPR + COSH*TMPI
C          EGF(I,J) = CMPLX(VALR,VALI)
C          SUM1 = SUM1 + VALR**2
C          SUM2 = SUM2 + VALI**2
C        END DO
C        ANORM = SQRT(SUM1+SUM2)
        SUM = 0.0
        DO I=1,NMD
          SUM = SUM + EGF(I,J)*CONJG(EGF(I,J))
        END DO
        ANORM = SQRT(SUM)
        DO I=1,NMD
          EGF(I,J) = EGF(I,J)/ANORM
          EGFA(I,J) = EGFA(I,J)*ANORM
        END DO
70    CONTINUE
      IF (ISRT.EQ.1)  CALL EIGSRT(EGV,EGF,EGFA,EFLD,PERD,NPOP,NMD,NMAX)

C ------- Write Eigenmodes and Modal Contributions
      DO 85 I=1,MIN(NPOP,NPRT)
        WRITE(7,75) EGV(I), CABS(EGV(I)), EFLD(I), PERD(I)
75      FORMAT(5X,5E13.5,/)
        WRITE(8,80) (EGF(K,I), K=1,NMD)
80      FORMAT(6E13.5)
        NMODE = I
85    CONTINUE
      PRINT *, NMODE

C ------- Read EOF Patterns if POPs are Requested
      IF ((MOD(IOUT,2).EQ.1) .OR. (MOD(IOUT/2,2).EQ.1)) THEN
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
      END IF

      DO 130 IM=1,NCYC
        IU = IM+10
C ------- Eigenfunctions of Different Models
        IF (IM.EQ.1)  GO TO 100
        DO J=1,NPOP
          IF (PERD(J).EQ.0.) THEN
            ETAN = 0.0
          ELSE
            ETAN = TPI*DT/PERD(J)/FLOAT(NCYC)
          END IF
          DO I=1,NMD
            CSUM = (0.,0.)
            DO K=1,NMD
              CSUM = CSUM + COV(I,K,IM-1)*EGF(K,J)
            END DO
            FCT(I) = CSUM*CEXP(-JIMAG*ETAN)
          END DO
          DO I=1,NMD
            EGF(I,J) = FCT(I)
          END DO
        END DO

C ------- Adjoint Patterns
        JPOP = 1
        DO 90 J=1,NPOP
          DO I=1,NMD
            SIG0(JPOP,I) = REAL(EGF(I,J))
          END DO
          JPOP = JPOP + 1
          IF (AIMAG(EGV(J)).EQ.0.)  GO TO 90
          DO I=1,NMD
            SIG0(JPOP,I) = AIMAG(EGF(I,J))
          END DO
          JPOP = JPOP + 1
90      CONTINUE
        IF (JPOP-1.NE.NMD)  PRINT *, '  ERROR in computing adjoint.'
C ---------- Adjoint as an inverse of eigenmatrix
        CALL LUDCMP(SIG0,NMD,NMAX,INDX,SIGN)
        DO I=1,NMD
          DO J=1,NMD
            DMTRX(I,J) = 0.0
          END DO
          DMTRX(I,I) = 1.0
        END DO
        DO J=1,NMD
          CALL LUBKSB(SIG0,NMD,NMAX,INDX,DMTRX(1,J))
        END DO

        JCON = 0
        JPOP = 0
        DO 95 J=1,NMD
C ---------- Rearrange adjoint vectors
          IF (JCON.EQ.1) THEN
            JCON = 0
            GO TO 95
          END IF

          JPOP = JPOP + 1
          IF (AIMAG(EGV(JPOP)).EQ.0.) THEN
            DO I=1,NMD
              EGFA(I,JPOP) = DMTRX(I,J)
            END DO
          ELSE
            DO I=1,NMD
              EGFA(I,JPOP) = 0.5*(DMTRX(I,J) + JIMAG*DMTRX(I,J+1))
            END DO
            JCON = 1
          END IF

C ---------- Normalization
          SUM = 0.0
          DO I=1,NMD
            SUM = SUM + EGF(I,JPOP)*CONJG(EGF(I,JPOP))
          END DO
          ANORM = SQRT(SUM)
          DO I=1,NMD
            EGF(I,JPOP) = EGF(I,JPOP)/ANORM
            EGFA(I,JPOP) = EGFA(I,JPOP)*ANORM
          END DO
95      CONTINUE

C ------- Write POP Patterns
100     CONTINUE
        IF ((MOD(IOUT,2).EQ.0) .AND. (MOD(IOUT/2,2).EQ.0))  GO TO 110
        DO J=1,MIN(NPOP,NPRT)
          DO K=1,NST
            POP(K) = (0.,0.)
          END DO
          DO 105 I=1,NMD
          DO 105 K=1,NST
            POP(K) = POP(K) + EGF(I,J)*EOF(K,I)
105       CONTINUE
          DO L=1,NY
            LL = (L-1)*NX
            WRITE(IU,80) (POP(K+LL)*SCL, K=1,NX)
          END DO
        END DO

110     CONTINUE
C ------- PC Time Series
        IF ((MOD(IOUT/4,2).EQ.0) .AND. (MOD(IOUT/8,2).EQ.0))  GO TO 130
        DO 120 J=1,MIN(NPOP,NPRT)
          DO I=IM,NPTS,NCYC
            IF (IADJ.EQ.0) THEN
              SUMT1 = 0.0
              SUMT2 = 0.0
              SUM1 = 0.0
              SUM2 = 0.0
              SUMC = 0.0
              DO K=1,NMD
                SUMT1 = SUMT1 + TSER(I,K)*REAL(EGF(K,J))
                SUMT2 = SUMT2 + TSER(I,K)*AIMAG(EGF(K,J))
                SUM1 = SUM1 + REAL(EGF(K,J))*REAL(EGF(K,J))
                SUM2 = SUM2 + AIMAG(EGF(K,J))*AIMAG(EGF(K,J))
                SUMC = SUMC + REAL(EGF(K,J))*AIMAG(EGF(K,J))
              END DO
              IF (SUM2.EQ.0.) THEN
                PCTS(I,J) = SUMT1 / SUM1
              ELSE IF (SUM1.EQ.0.) THEN
                PCTS(I,J) = JIMAG*SUMT2 / SUM2
              ELSE
                PCTS(I,J) = ((SUMT1*SUM2 - SUMT2*SUMC)
     &                      + JIMAG*(SUMT2*SUM1 - SUMT1*SUMC))
     &                    / (SUM1*SUM2 - SUMC**2)
              END IF
            ELSE
              SUM1 = 0.0
              SUM2 = 0.0
              DO K=1,NMD
                SUM1 = SUM1 + TSER(I,K)*REAL(EGFA(K,J))
                SUM2 = SUM2 + TSER(I,K)*AIMAG(EGFA(K,J))
              END DO
              PCTS(I,J) = SUM1 + JIMAG*SUM2
            END IF
          END DO
120     CONTINUE
130   CONTINUE
      
      IF ((MOD(IOUT/4,2).EQ.0) .AND. (MOD(IOUT/8,2).EQ.0))  GO TO 150
      DO 140 J=1,MIN(NPOP,NPRT)
        IF (MOD(IOUT/4,2).EQ.1)  WRITE(10,80)  (PCTS(I,J), I=1,NPTS)
        IF (MOD(IOUT/8,2).EQ.0)  GO TO 140
        DO I=1,NPTS
          AMPL = CABS(PCTS(I,J))
          PHSE = ATAN2(AIMAG(PCTS(I,J)), REAL(PCTS(I,J)))*RTD
          PCTS(I,J) = AMPL + JIMAG*PHSE
        END DO
        WRITE(10,80)  (REAL(PCTS(I,J)), I=1,NPTS)
        WRITE(10,80)  (AIMAG(PCTS(I,J)), I=1,NPTS)
140   CONTINUE

150   STOP
      END

      SUBROUTINE BALANC(NM,N,A,LOW,IGH,SCALE)
C
      INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
      REAL A(NM,N),SCALE(N)
      REAL C,F,G,R,S,B2,RADIX
      LOGICAL NOCONV
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BALANCE,
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
C     THIS SUBROUTINE BALANCES A REAL MATRIX AND ISOLATES
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
C        A CONTAINS THE INPUT MATRIX TO BE BALANCED.
C
C     ON OUTPUT
C
C        A CONTAINS THE BALANCED MATRIX.
C
C        LOW AND IGH ARE TWO INTEGERS SUCH THAT A(I,J)
C          IS EQUAL TO ZERO IF
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
C                 = D(J,J),      J = LOW,...,IGH
C                 = P(J)         J = IGH+1,...,N.
C     THE ORDER IN WHICH THE INTERCHANGES ARE MADE IS N TO IGH+1,
C     THEN 1 TO LOW-1.
C
C     NOTE THAT 1 IS RETURNED FOR IGH IF IGH IS ZERO FORMALLY.
C
C     THE ALGOL PROCEDURE EXC CONTAINED IN BALANCE APPEARS IN
C     BALANC  IN LINE.  (NOTE THAT THE ALGOL ROLES OF IDENTIFIERS
C     K,L HAVE BEEN REVERSED.)
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
         F = A(I,J)
         A(I,J) = A(I,M)
         A(I,M) = F
   30 CONTINUE
C
      DO 40 I = K, N
         F = A(J,I)
         A(J,I) = A(M,I)
         A(M,I) = F
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
            IF (A(J,I) .NE. 0.0E0) GO TO 120
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
            IF (A(I,J) .NE. 0.0E0) GO TO 170
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
            C = C + ABS(A(J,I))
            R = R + ABS(A(I,J))
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
  250    A(I,J) = A(I,J) * G
C
         DO 260 J = 1, L
  260    A(J,I) = A(J,I) * F
C
  270 CONTINUE
C
      IF (NOCONV) GO TO 190
C
  280 LOW = K
      IGH = L
      RETURN
      END

      SUBROUTINE ELMHES(NM,N,LOW,IGH,A,INT)
C
      INTEGER I,J,M,N,LA,NM,IGH,KP1,LOW,MM1,MP1
      REAL A(NM,N)
      REAL X,Y
      INTEGER INT(IGH)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ELMHES,
C     NUM. MATH. 12, 349-368(1968) BY MARTIN AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
C
C     GIVEN A REAL GENERAL MATRIX, THIS SUBROUTINE
C     REDUCES A SUBMATRIX SITUATED IN ROWS AND COLUMNS
C     LOW THROUGH IGH TO UPPER HESSENBERG FORM BY
C     STABILIZED ELEMENTARY SIMILARITY TRANSFORMATIONS.
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
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        A CONTAINS THE INPUT MATRIX.
C
C     ON OUTPUT
C
C        A CONTAINS THE HESSENBERG MATRIX.  THE MULTIPLIERS
C          WHICH WERE USED IN THE REDUCTION ARE STORED IN THE
C          REMAINING TRIANGLE UNDER THE HESSENBERG MATRIX.
C
C        INT CONTAINS INFORMATION ON THE ROWS AND COLUMNS
C          INTERCHANGED IN THE REDUCTION.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
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
         MM1 = M - 1
         X = 0.0E0
         I = M
C
         DO 100 J = M, IGH
            IF (ABS(A(J,MM1)) .LE. ABS(X)) GO TO 100
            X = A(J,MM1)
            I = J
  100    CONTINUE
C
         INT(M) = I
         IF (I .EQ. M) GO TO 130
C     .......... INTERCHANGE ROWS AND COLUMNS OF A ..........
         DO 110 J = MM1, N
            Y = A(I,J)
            A(I,J) = A(M,J)
            A(M,J) = Y
  110    CONTINUE
C
         DO 120 J = 1, IGH
            Y = A(J,I)
            A(J,I) = A(J,M)
            A(J,M) = Y
  120    CONTINUE
C     .......... END INTERCHANGE ..........
  130    IF (X .EQ. 0.0E0) GO TO 180
         MP1 = M + 1
C
         DO 160 I = MP1, IGH
            Y = A(I,MM1)
            IF (Y .EQ. 0.0E0) GO TO 160
            Y = Y / X
            A(I,MM1) = Y
C
            DO 140 J = M, N
  140       A(I,J) = A(I,J) - Y * A(M,J)
C
            DO 150 J = 1, IGH
  150       A(J,M) = A(J,M) + Y * A(J,I)
C
  160    CONTINUE
C
  180 CONTINUE
C
  200 RETURN
      END
 
      SUBROUTINE ELTRAN(NM,N,LOW,IGH,A,INT,Z)
C
      INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
      REAL A(NM,IGH),Z(NM,N)
      INTEGER INT(IGH)
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE ELMTRANS,
C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     THIS SUBROUTINE ACCUMULATES THE STABILIZED ELEMENTARY
C     SIMILARITY TRANSFORMATIONS USED IN THE REDUCTION OF A
C     REAL GENERAL MATRIX TO UPPER HESSENBERG FORM BY  ELMHES.
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
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        A CONTAINS THE MULTIPLIERS WHICH WERE USED IN THE
C          REDUCTION BY  ELMHES  IN ITS LOWER TRIANGLE
C          BELOW THE SUBDIAGONAL.
C
C        INT CONTAINS INFORMATION ON THE ROWS AND COLUMNS
C          INTERCHANGED IN THE REDUCTION BY  ELMHES.
C          ONLY ELEMENTS LOW THROUGH IGH ARE USED.
C
C     ON OUTPUT
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE
C          REDUCTION BY  ELMHES.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
C     .......... INITIALIZE Z TO IDENTITY MATRIX ..........
      DO 80 J = 1, N
C
         DO 60 I = 1, N
   60    Z(I,J) = 0.0E0
C
         Z(J,J) = 1.0E0
   80 CONTINUE
C
      KL = IGH - LOW - 1
      IF (KL .LT. 1) GO TO 200
C     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
      DO 140 MM = 1, KL
         MP = IGH - MM
         MP1 = MP + 1
C
         DO 100 I = MP1, IGH
  100    Z(I,MP) = A(I,MP-1)
C
         I = INT(MP)
         IF (I .EQ. MP) GO TO 140
C
         DO 130 J = MP, IGH
            Z(MP,J) = Z(I,J)
            Z(I,J) = 0.0E0
  130    CONTINUE
C
         Z(I,MP) = 1.0E0
  140 CONTINUE
C
  200 RETURN
      END
      SUBROUTINE HQR2(NM,N,LOW,IGH,H,WR,WI,Z,IERR)
C
      INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NA,NM,NN,
     X        IGH,ITN,ITS,LOW,MP2,ENM2,IERR
      REAL H(NM,N),WR(N),WI(N),Z(NM,N)
      REAL P,Q,R,S,T,W,X,Y,RA,SA,VI,VR,ZZ,NORM,TST1,TST2
      LOGICAL NOTLAS
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE HQR2,
C     NUM. MATH. 16, 181-204(1970) BY PETERS AND WILKINSON.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
C
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS
C     OF A REAL UPPER HESSENBERG MATRIX BY THE QR METHOD.  THE
C     EIGENVECTORS OF A REAL GENERAL MATRIX CAN ALSO BE FOUND
C     IF  ELMHES  AND  ELTRAN  OR  ORTHES  AND  ORTRAN  HAVE
C     BEEN USED TO REDUCE THIS GENERAL MATRIX TO HESSENBERG FORM
C     AND TO ACCUMULATE THE SIMILARITY TRANSFORMATIONS.
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
C          SUBROUTINE  BALANC.  IF  BALANC  HAS NOT BEEN USED,
C          SET LOW=1, IGH=N.
C
C        H CONTAINS THE UPPER HESSENBERG MATRIX.
C
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED BY  ELTRAN
C          AFTER THE REDUCTION BY  ELMHES, OR BY  ORTRAN  AFTER THE
C          REDUCTION BY  ORTHES, IF PERFORMED.  IF THE EIGENVECTORS
C          OF THE HESSENBERG MATRIX ARE DESIRED, Z MUST CONTAIN THE
C          IDENTITY MATRIX.
C
C     ON OUTPUT
C
C        H HAS BEEN DESTROYED.
C
C        WR AND WI CONTAIN THE REAL AND IMAGINARY PARTS,
C          RESPECTIVELY, OF THE EIGENVALUES.  THE EIGENVALUES
C          ARE UNORDERED EXCEPT THAT COMPLEX CONJUGATE PAIRS
C          OF VALUES APPEAR CONSECUTIVELY WITH THE EIGENVALUE
C          HAVING THE POSITIVE IMAGINARY PART FIRST.  IF AN
C          ERROR EXIT IS MADE, THE EIGENVALUES SHOULD BE CORRECT
C          FOR INDICES IERR+1,...,N.
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGENVECTORS.
C          IF THE I-TH EIGENVALUE IS REAL, THE I-TH COLUMN OF Z
C          CONTAINS ITS EIGENVECTOR.  IF THE I-TH EIGENVALUE IS COMPLEX
C          WITH POSITIVE IMAGINARY PART, THE I-TH AND (I+1)-TH
C          COLUMNS OF Z CONTAIN THE REAL AND IMAGINARY PARTS OF ITS
C          EIGENVECTOR.  THE EIGENVECTORS ARE UNNORMALIZED.  IF AN
C          ERROR EXIT IS MADE, NONE OF THE EIGENVECTORS HAS BEEN FOUND.
C
C        IERR IS SET TO
C          ZERO       FOR NORMAL RETURN,
C          J          IF THE LIMIT OF 30*N ITERATIONS IS EXHAUSTED
C                     WHILE THE J-TH EIGENVALUE IS BEING SOUGHT.
C
C     CALLS CDIV FOR COMPLEX DIVISION.
C
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO BURTON S. GARBOW,
C     MATHEMATICS AND COMPUTER SCIENCE DIV, ARGONNE NATIONAL LABORATORY
C
C     THIS VERSION DATED AUGUST 1983.
C
C     ------------------------------------------------------------------
C
      IERR = 0
      NORM = 0.0E0
      K = 1
C     .......... STORE ROOTS ISOLATED BY BALANC
C                AND COMPUTE MATRIX NORM ..........
      DO 50 I = 1, N
C
         DO 40 J = K, N
   40    NORM = NORM + ABS(H(I,J))
C
         K = I
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 50
         WR(I) = H(I,I)
         WI(I) = 0.0E0
   50 CONTINUE
C
      EN = IGH
      T = 0.0E0
      ITN = 30*N
C     .......... SEARCH FOR NEXT EIGENVALUES ..........
   60 IF (EN .LT. LOW) GO TO 340
      ITS = 0
      NA = EN - 1
      ENM2 = NA - 1
C     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
C                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
   70 DO 80 LL = LOW, EN
         L = EN + LOW - LL
         IF (L .EQ. LOW) GO TO 100
         S = ABS(H(L-1,L-1)) + ABS(H(L,L))
         IF (S .EQ. 0.0E0) S = NORM
         TST1 = S
         TST2 = TST1 + ABS(H(L,L-1))
         IF (TST2 .EQ. TST1) GO TO 100
   80 CONTINUE
C     .......... FORM SHIFT ..........
  100 X = H(EN,EN)
      IF (L .EQ. EN) GO TO 270
      Y = H(NA,NA)
      W = H(EN,NA) * H(NA,EN)
      IF (L .EQ. NA) GO TO 280
      IF (ITN .EQ. 0) GO TO 1000
      IF (ITS .NE. 10 .AND. ITS .NE. 20) GO TO 130
C     .......... FORM EXCEPTIONAL SHIFT ..........
      T = T + X
C
      DO 120 I = LOW, EN
  120 H(I,I) = H(I,I) - X
C
      S = ABS(H(EN,NA)) + ABS(H(NA,ENM2))
      X = 0.75E0 * S
      Y = X
      W = -0.4375E0 * S * S
  130 ITS = ITS + 1
      ITN = ITN - 1
C     .......... LOOK FOR TWO CONSECUTIVE SMALL
C                SUB-DIAGONAL ELEMENTS.
C                FOR M=EN-2 STEP -1 UNTIL L DO -- ..........
      DO 140 MM = L, ENM2
         M = ENM2 + L - MM
         ZZ = H(M,M)
         R = X - ZZ
         S = Y - ZZ
         P = (R * S - W) / H(M+1,M) + H(M,M+1)
         Q = H(M+1,M+1) - ZZ - R - S
         R = H(M+2,M+1)
         S = ABS(P) + ABS(Q) + ABS(R)
         P = P / S
         Q = Q / S
         R = R / S
         IF (M .EQ. L) GO TO 150
         TST1 = ABS(P)*(ABS(H(M-1,M-1)) + ABS(ZZ) + ABS(H(M+1,M+1)))
         TST2 = TST1 + ABS(H(M,M-1))*(ABS(Q) + ABS(R))
         IF (TST2 .EQ. TST1) GO TO 150
  140 CONTINUE
C
  150 MP2 = M + 2
C
      DO 160 I = MP2, EN
         H(I,I-2) = 0.0E0
         IF (I .EQ. MP2) GO TO 160
         H(I,I-3) = 0.0E0
  160 CONTINUE
C     .......... DOUBLE QR STEP INVOLVING ROWS L TO EN AND
C                COLUMNS M TO EN ..........
      DO 260 K = M, NA
         NOTLAS = K .NE. NA
         IF (K .EQ. M) GO TO 170
         P = H(K,K-1)
         Q = H(K+1,K-1)
         R = 0.0E0
         IF (NOTLAS) R = H(K+2,K-1)
         X = ABS(P) + ABS(Q) + ABS(R)
         IF (X .EQ. 0.0E0) GO TO 260
         P = P / X
         Q = Q / X
         R = R / X
  170    S = SIGN(SQRT(P*P+Q*Q+R*R),P)
         IF (K .EQ. M) GO TO 180
         H(K,K-1) = -S * X
         GO TO 190
  180    IF (L .NE. M) H(K,K-1) = -H(K,K-1)
  190    P = P + S
         X = P / S
         Y = Q / S
         ZZ = R / S
         Q = Q / P
         R = R / P
         IF (NOTLAS) GO TO 225
C     .......... ROW MODIFICATION ..........
         DO 200 J = K, N
            P = H(K,J) + Q * H(K+1,J)
            H(K,J) = H(K,J) - P * X
            H(K+1,J) = H(K+1,J) - P * Y
  200    CONTINUE
C
         J = MIN0(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
         DO 210 I = 1, J
            P = X * H(I,K) + Y * H(I,K+1)
            H(I,K) = H(I,K) - P
            H(I,K+1) = H(I,K+1) - P * Q
  210    CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
         DO 220 I = LOW, IGH
            P = X * Z(I,K) + Y * Z(I,K+1)
            Z(I,K) = Z(I,K) - P
            Z(I,K+1) = Z(I,K+1) - P * Q
  220    CONTINUE
         GO TO 255
  225    CONTINUE
C     .......... ROW MODIFICATION ..........
         DO 230 J = K, N
            P = H(K,J) + Q * H(K+1,J) + R * H(K+2,J)
            H(K,J) = H(K,J) - P * X
            H(K+1,J) = H(K+1,J) - P * Y
            H(K+2,J) = H(K+2,J) - P * ZZ
  230    CONTINUE
C
         J = MIN0(EN,K+3)
C     .......... COLUMN MODIFICATION ..........
         DO 240 I = 1, J
            P = X * H(I,K) + Y * H(I,K+1) + ZZ * H(I,K+2)
            H(I,K) = H(I,K) - P
            H(I,K+1) = H(I,K+1) - P * Q
            H(I,K+2) = H(I,K+2) - P * R
  240    CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
         DO 250 I = LOW, IGH
            P = X * Z(I,K) + Y * Z(I,K+1) + ZZ * Z(I,K+2)
            Z(I,K) = Z(I,K) - P
            Z(I,K+1) = Z(I,K+1) - P * Q
            Z(I,K+2) = Z(I,K+2) - P * R
  250    CONTINUE
  255    CONTINUE
C
  260 CONTINUE
C
      GO TO 70
C     .......... ONE ROOT FOUND ..........
  270 H(EN,EN) = X + T
      WR(EN) = H(EN,EN)
      WI(EN) = 0.0E0
      EN = NA
      GO TO 60
C     .......... TWO ROOTS FOUND ..........
  280 P = (Y - X) / 2.0E0
      Q = P * P + W
      ZZ = SQRT(ABS(Q))
      H(EN,EN) = X + T
      X = H(EN,EN)
      H(NA,NA) = Y + T
      IF (Q .LT. 0.0E0) GO TO 320
C     .......... REAL PAIR ..........
      ZZ = P + SIGN(ZZ,P)
      WR(NA) = X + ZZ
      WR(EN) = WR(NA)
      IF (ZZ .NE. 0.0E0) WR(EN) = X - W / ZZ
      WI(NA) = 0.0E0
      WI(EN) = 0.0E0
      X = H(EN,NA)
      S = ABS(X) + ABS(ZZ)
      P = X / S
      Q = ZZ / S
      R = SQRT(P*P+Q*Q)
      P = P / R
      Q = Q / R
C     .......... ROW MODIFICATION ..........
      DO 290 J = NA, N
         ZZ = H(NA,J)
         H(NA,J) = Q * ZZ + P * H(EN,J)
         H(EN,J) = Q * H(EN,J) - P * ZZ
  290 CONTINUE
C     .......... COLUMN MODIFICATION ..........
      DO 300 I = 1, EN
         ZZ = H(I,NA)
         H(I,NA) = Q * ZZ + P * H(I,EN)
         H(I,EN) = Q * H(I,EN) - P * ZZ
  300 CONTINUE
C     .......... ACCUMULATE TRANSFORMATIONS ..........
      DO 310 I = LOW, IGH
         ZZ = Z(I,NA)
         Z(I,NA) = Q * ZZ + P * Z(I,EN)
         Z(I,EN) = Q * Z(I,EN) - P * ZZ
  310 CONTINUE
C
      GO TO 330
C     .......... COMPLEX PAIR ..........
  320 WR(NA) = X + P
      WR(EN) = X + P
      WI(NA) = ZZ
      WI(EN) = -ZZ
  330 EN = ENM2
      GO TO 60
C     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
C                VECTORS OF UPPER TRIANGULAR FORM ..........
  340 IF (NORM .EQ. 0.0E0) GO TO 1001
C     .......... FOR EN=N STEP -1 UNTIL 1 DO -- ..........
      DO 800 NN = 1, N
         EN = N + 1 - NN
         P = WR(EN)
         Q = WI(EN)
         NA = EN - 1
         IF (Q) 710, 600, 800
C     .......... REAL VECTOR ..........
  600    M = EN
         H(EN,EN) = 1.0E0
         IF (NA .EQ. 0) GO TO 800
C     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
         DO 700 II = 1, NA
            I = EN - II
            W = H(I,I) - P
            R = 0.0E0
C
            DO 610 J = M, EN
  610       R = R + H(I,J) * H(J,EN)
C
            IF (WI(I) .GE. 0.0E0) GO TO 630
            ZZ = W
            S = R
            GO TO 700
  630       M = I
            IF (WI(I) .NE. 0.0E0) GO TO 640
            T = W
            IF (T .NE. 0.0E0) GO TO 635
               TST1 = NORM
               T = TST1
  632          T = 0.01E0 * T
               TST2 = NORM + T
               IF (TST2 .GT. TST1) GO TO 632
  635       H(I,EN) = -R / T
            GO TO 680
C     .......... SOLVE REAL EQUATIONS ..........
  640       X = H(I,I+1)
            Y = H(I+1,I)
            Q = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I)
            T = (X * S - ZZ * R) / Q
            H(I,EN) = T
            IF (ABS(X) .LE. ABS(ZZ)) GO TO 650
            H(I+1,EN) = (-R - W * T) / X
            GO TO 680
  650       H(I+1,EN) = (-S - Y * T) / ZZ
C
C     .......... OVERFLOW CONTROL ..........
  680       T = ABS(H(I,EN))
            IF (T .EQ. 0.0E0) GO TO 700
            TST1 = T
            TST2 = TST1 + 1.0E0/TST1
            IF (TST2 .GT. TST1) GO TO 700
            DO 690 J = I, EN
               H(J,EN) = H(J,EN)/T
  690       CONTINUE
C
  700    CONTINUE
C     .......... END REAL VECTOR ..........
         GO TO 800
C     .......... COMPLEX VECTOR ..........
  710    M = NA
C     .......... LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
C                EIGENVECTOR MATRIX IS TRIANGULAR ..........
         IF (ABS(H(EN,NA)) .LE. ABS(H(NA,EN))) GO TO 720
         H(NA,NA) = Q / H(EN,NA)
         H(NA,EN) = -(H(EN,EN) - P) / H(EN,NA)
         GO TO 730
  720    CALL CDIV(0.0E0,-H(NA,EN),H(NA,NA)-P,Q,H(NA,NA),H(NA,EN))
  730    H(EN,NA) = 0.0E0
         H(EN,EN) = 1.0E0
         ENM2 = NA - 1
         IF (ENM2 .EQ. 0) GO TO 800
C     .......... FOR I=EN-2 STEP -1 UNTIL 1 DO -- ..........
         DO 795 II = 1, ENM2
            I = NA - II
            W = H(I,I) - P
            RA = 0.0E0
            SA = 0.0E0
C
            DO 760 J = M, EN
               RA = RA + H(I,J) * H(J,NA)
               SA = SA + H(I,J) * H(J,EN)
  760       CONTINUE
C
            IF (WI(I) .GE. 0.0E0) GO TO 770
            ZZ = W
            R = RA
            S = SA
            GO TO 795
  770       M = I
            IF (WI(I) .NE. 0.0E0) GO TO 780
            CALL CDIV(-RA,-SA,W,Q,H(I,NA),H(I,EN))
            GO TO 790
C     .......... SOLVE COMPLEX EQUATIONS ..........
  780       X = H(I,I+1)
            Y = H(I+1,I)
            VR = (WR(I) - P) * (WR(I) - P) + WI(I) * WI(I) - Q * Q
            VI = (WR(I) - P) * 2.0E0 * Q
            IF (VR .NE. 0.0E0 .OR. VI .NE. 0.0E0) GO TO 784
               TST1 = NORM * (ABS(W) + ABS(Q) + ABS(X)
     X                      + ABS(Y) + ABS(ZZ))
               VR = TST1
  783          VR = 0.01E0 * VR
               TST2 = TST1 + VR
               IF (TST2 .GT. TST1) GO TO 783
  784       CALL CDIV(X*R-ZZ*RA+Q*SA,X*S-ZZ*SA-Q*RA,VR,VI,
     X                H(I,NA),H(I,EN))
            IF (ABS(X) .LE. ABS(ZZ) + ABS(Q)) GO TO 785
            H(I+1,NA) = (-RA - W * H(I,NA) + Q * H(I,EN)) / X
            H(I+1,EN) = (-SA - W * H(I,EN) - Q * H(I,NA)) / X
            GO TO 790
  785       CALL CDIV(-R-Y*H(I,NA),-S-Y*H(I,EN),ZZ,Q,
     X                H(I+1,NA),H(I+1,EN))
C
C     .......... OVERFLOW CONTROL ..........
  790       T = AMAX1(ABS(H(I,NA)), ABS(H(I,EN)))
            IF (T .EQ. 0.0E0) GO TO 795
            TST1 = T
            TST2 = TST1 + 1.0E0/TST1
            IF (TST2 .GT. TST1) GO TO 795
            DO 792 J = I, EN
               H(J,NA) = H(J,NA)/T
               H(J,EN) = H(J,EN)/T
  792       CONTINUE
C
  795    CONTINUE
C     .......... END COMPLEX VECTOR ..........
  800 CONTINUE
C     .......... END BACK SUBSTITUTION.
C                VECTORS OF ISOLATED ROOTS ..........
      DO 840 I = 1, N
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 840
C
         DO 820 J = I, N
  820    Z(I,J) = H(I,J)
C
  840 CONTINUE
C     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
C                VECTORS OF ORIGINAL FULL MATRIX.
C                FOR J=N STEP -1 UNTIL LOW DO -- ..........
      DO 880 JJ = LOW, N
         J = N + LOW - JJ
         M = MIN0(J,IGH)
C
         DO 880 I = LOW, IGH
            ZZ = 0.0E0
C
            DO 860 K = LOW, M
  860       ZZ = ZZ + Z(I,K) * H(K,J)
C
            Z(I,J) = ZZ
  880 CONTINUE
C
      GO TO 1001
C     .......... SET ERROR -- ALL EIGENVALUES HAVE NOT
C                CONVERGED AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
      END

      SUBROUTINE BALBAK(NM,N,LOW,IGH,SCALE,M,Z)
C
      INTEGER I,J,K,M,N,II,NM,IGH,LOW
      REAL SCALE(N),Z(NM,M)
      REAL S
C
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE BALBAK,
C     NUM. MATH. 13, 293-304(1969) BY PARLETT AND REINSCH.
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 315-326(1971).
C
C     THIS SUBROUTINE FORMS THE EIGENVECTORS OF A REAL GENERAL
C     MATRIX BY BACK TRANSFORMING THOSE OF THE CORRESPONDING
C     BALANCED MATRIX DETERMINED BY  BALANC.
C
C     ON INPUT
C
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM
C          DIMENSION STATEMENT.
C
C        N IS THE ORDER OF THE MATRIX.
C
C        LOW AND IGH ARE INTEGERS DETERMINED BY  BALANC.
C
C        SCALE CONTAINS INFORMATION DETERMINING THE PERMUTATIONS
C          AND SCALING FACTORS USED BY  BALANC.
C
C        M IS THE NUMBER OF COLUMNS OF Z TO BE BACK TRANSFORMED.
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE EIGEN-
C          VECTORS TO BE BACK TRANSFORMED IN ITS FIRST M COLUMNS.
C
C     ON OUTPUT
C
C        Z CONTAINS THE REAL AND IMAGINARY PARTS OF THE
C          TRANSFORMED EIGENVECTORS IN ITS FIRST M COLUMNS.
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
  100    Z(I,J) = Z(I,J) * S
C
  110 CONTINUE
C     ......... FOR I=LOW-1 STEP -1 UNTIL 1,
C               IGH+1 STEP 1 UNTIL N DO -- ..........
  120 DO 140 II = 1, N
         I = II
         IF (I .GE. LOW .AND. I .LE. IGH) GO TO 140
         IF (I .LT. LOW) I = LOW - II
         K = SCALE(I)
         IF (K .EQ. I) GO TO 140
C
         DO 130 J = 1, M
            S = Z(I,J)
            Z(I,J) = Z(K,J)
            Z(K,J) = S
  130    CONTINUE
C
  140 CONTINUE
C
  200 RETURN
      END

      SUBROUTINE EIGSRT(D,V,VV,ARR,BRR,N,M,NP)

      COMPLEX D(NP),V(NP,NP),VV(NP,NP),P
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
            P=VV(J,I)
            VV(J,I)=VV(J,K)
            VV(J,K)=P
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
