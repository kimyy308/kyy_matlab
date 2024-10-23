C *****************************************************************************
C
C       PROGRAM NAME : reof.f
C       PROGRAMMER : Dr. Kwang Y. Kim
C       CODE IDENTIFICATION NUMBER = REOF/VERSION 1.0
C       CODE CLASSIFICATION = Scientific Computer Code
C       CREATION DATE = February 10, 1998
C       REVISION DATE = not revised
C       REVISION INFORMATION = not applicable
C
C *****************************************************************************


C     This program computes rotated EOFs based on the varimax rotation.
C     The program consists of calling the varimax rotation subroutine for
C     a given set of guess patterns (EOFs).  The subroutine call is:
C
C          CALL TVARMX(NST,NMD,AMAT,NMAX,D,H,MRAW)
C
C     where
C
C          AMAT:  original loading matrix (eigenvectors)
C          RMAT:  rotated loading matrix (eigenvectors)
C          NST:   number of stations
C          NMD:   number of modes used in rotation
C          NMAX:  first physical dimension of AMAT
C          H:     vector containing communalities of attributes
C                    output must have at least NST elements
C          MRAW:  an integer concerning case of varimax
C                    = 0  for unweighted (raw) varimax solution
C                    = 1  for standardized varimax solution


      PARAMETER (NMAX=615, MMD=50, NTMAX=600)
      DIMENSION AMAT(NMAX,MMD), RMAT(NMAX,MMD), H(NMAX), D(51)
      DIMENSION ARR(NMAX,NTMAX), PCT(NTMAX,MMD), EGV(MMD)
      DIMENSION SIG(MMD,MMD), SOL(MMD,MMD), INDX(MMD)
      CHARACTER*50 FILNM1, FILNM2, FILNM3, FILNM4,
     &             FORMT1, FORMT2, FORMT3, FORMT4


5     FORMAT(A50)
      PRINT *, '  Type name and format of the input eigenfunction file.'
      READ 5, FILNM1, FORMT1
      PRINT *, '  Type name and format of the output REOF file.'
      READ 5, FILNM3, FORMT3
      PRINT *, '  Type the spatial dimensions (NX x NY) of the data.'
      READ *, NX, NY
      NST = NX*NY
      PRINT *, '  Type the number of retained modes.'
      READ *, NMD
      PRINT *, '  Type variance loading option.'
      PRINT *, '     0: unloaded normalized eigenvector.'
      PRINT *, '     1: eigenvector multiplied by sqrt of eigenvalue.'
      READ *, IVAR
      IF (IVAR.EQ.1) THEN
        PRINT *, '  Type name and format of the input eigenvalue file.'
        READ 5, FILNM2, FORMT2
      END IF
      PRINT *, '  Type the varimax standardization option.'
      PRINT *, '     0: for unweighted (raw) varimax solution'
      PRINT *, '     1: for standardized varimax solution'
      READ *, MRAW
      PRINT *, '  Type the REOF scaling factor.'
      READ *, SCL
      PRINT *, '  Type the REOF orthogonality test option.'
      READ *, ITST
      PRINT *, '  Do you want PC time series?  0: No,  1: Yes'
      READ *, IPC
      IF (IPC.EQ.1) THEN
        PRINT *, '  Type name and format of the input dataset.'
        READ 5, FILNM4, FORMT4
        PRINT *, '  Type smoothing option (0:No, M:Moving average lag).'
        READ *, LAG
        PRINT *, '  Type mean removing option (0:Annual, 1:Ann Cycle).'
        READ *, IDTR
      END IF


C ------- Read eigenfunctions and eigenvectors
      OPEN(UNIT=1, FILE=FILNM1, STATUS='OLD')
      DO 10 I=1,NMD
      DO L=1,NY
        LL = (L-1)*NX
        READ(1,FORMT1)  (AMAT(K+LL,I), K=1,NX)
        END DO
10    CONTINUE
      IF (IVAR.EQ.1) THEN
        OPEN(UNIT=2, FILE=FILNM2, STATUS='OLD')
        READ(2,'(22X,E15.7,///)')  TVAR
        DO I=1,NMD
          READ(2,FORMT2)  EGV(I)
        END DO
      END IF

C ------- Copy AMAT into RMAT
      DO 15 J=1,NMD
      DO 15 I=1,NST
        IF (IVAR.EQ.0) THEN
          RMAT(I,J) = AMAT(I,J)
        ELSE
          RMAT(I,J) = AMAT(I,J)*SQRT(EGV(J))
        END IF
15    CONTINUE

C ------- Call TVARMX
      CALL TVARMX(NST,NMD,RMAT,NMAX,D,H,MRAW)

C ------- Normalization of rotated EOFs
      DO I=1,NMD
        SUM = 0.0
        DO 20 L=1,NY
          LL = (L-1)*NX
        DO 20 K=1,NX
          SUM = SUM + RMAT(K+LL,I)**2
20      CONTINUE
C        SUM = SUM/FLOAT(NST)

        DO 25 L=1,NY
          LL = (L-1)*NX
        DO 25 K=1,NX
          RMAT(K+LL,I) = RMAT(K+LL,I)/SQRT(SUM)
25      CONTINUE
      END DO

      IF (IPC.EQ.0 .AND. ITST.EQ.0)  GO TO 60
C ------- Adjoint is defined by Inv(P'P) P'
      IF (IVAR.EQ.0) THEN
        DO 30 IM=1,NMD
        DO 30 K=1,NST
          AMAT(K,IM) = RMAT(K,IM)
30      CONTINUE
      ELSE
        DO 35 IM=1,NMD
        DO 35 JM=1,NMD
          SUM = 0.0
          DO K=1,NST
            SUM = SUM + RMAT(K,IM)*RMAT(K,JM)
          END DO
          SIG(IM,JM) = SUM
35      CONTINUE
        CALL LUDCMP(SIG,NMD,MMD,INDX,SIGN)
        DO I=1,NMD
          DO J=1,NMD
            SOL(I,J) = 0.0
          END DO
          SOL(I,I) = 1.0
        END DO
        DO J=1,NMD
          CALL LUBKSB(SIG,NMD,MMD,INDX,SOL(1,J))
        END DO

        DO 40 IM=1,NMD
        DO 40 K=1,NST
          SUM = 0.0
          DO JM=1,NMD
            SUM = SUM + SOL(IM,JM)*RMAT(K,JM)
          END DO
          AMAT(K,IM) = SUM
40      CONTINUE
      END IF

C ------- test of orthogonality
      IF (ITST.EQ.1) THEN
      DO 50 IM=1,NMD
      DO 50 JM=1,NMD
        SUM = 0.0
        DO 45 L=1,NY
          LL = (L-1)*NX
        DO 45 K=1,NX
          SUM = SUM + RMAT(K+LL,IM)*AMAT(K+LL,JM)
45      CONTINUE
C        SUM = SUM/FLOAT(NST)
        PRINT *, IM, JM, SUM
50    CONTINUE
      END IF

60    IF (IPC.EQ.0)  GO TO 100
C ------- Read input data
      IF (FORMT4.EQ.'UNF') THEN
        ICNT = 1
        OPEN(UNIT=3, FILE=FILNM4, STATUS='OLD', FORM='UNFORMATTED')
65      READ(3,END=80)  ((ARR(K+(L-1)*NX,ICNT), K=1,NX), L=1,NY)
        ICNT = ICNT+1
        GO TO 65
      ELSE
        ICNT = 1
        OPEN(UNIT=3, FILE=FILNM4, STATUS='OLD')
70      DO L=1,NY
          LL = (L-1)*NX
          READ(3,FORMT4,END=80)  (ARR(K+LL,ICNT), K=1,NX)
        END DO
        ICNT = ICNT+1
        GO TO 70
      END IF

C ------- Smoothing
80    CONTINUE
      NCNT = ICNT-1
      DO L=1,NMAX
        DO I=1,NCNT
          KS = MAX(I-LAG,1)
          KE = MIN(I+LAG,NCNT)
          KN = KE-KS+1
          SUM = 0.0
          DO K=KS,KE
            SUM = SUM + ARR(L,K)
          END DO
          PCT(I,1) = SUM/FLOAT(KN)
        END DO
        DO I=1,NCNT
          ARR(L,I) = PCT(I,1)
        END DO
      END DO

C ------- Subtract mean
      NCYC = 1
      IF (IDTR.NE.0)  NCYC=12
      DO 85 L=1,NMAX
      DO 85 IM=1,NCYC
        SUM = 0.0
        DO K=IM,NCNT,NCYC
          SUM = SUM + ARR(L,K)
        END DO
        AVG = SUM/FLOAT(NCNT/NCYC)
        DO K=IM,NCNT,NCYC
          ARR(L,K) = ARR(L,K) - AVG
        END DO
85    CONTINUE

C ------- PC time series
      DO 95 I=1,NCNT
      DO 95 IM=1,NMD
        SUM = 0.0
        DO 90 L=1,NY
          LL = (L-1)*NX
        DO 90 K=1,NX
          SUM = SUM + ARR(K+LL,I)*AMAT(K+LL,IM)
90      CONTINUE
        PCT(I,IM) = SUM/FLOAT(NST)
95    CONTINUE

C ------- Sorting by eigenvalues
      TVAR = 0.0
      DO IM=1,NMD
        AVG = 0.0
        VAR = 0.0
        DO I=1,NCNT
          AVG = AVG + PCT(I,IM)
          VAR = VAR + PCT(I,IM)**2
        END DO
        AVG = AVG/FLOAT(NCNT)
        VAR = VAR/FLOAT(NCNT) - AVG**2
        EGV(IM) = VAR
        TVAR = TVAR + VAR
      END DO
      CALL EIGSRT(EGV,RMAT,PCT,NMAX,NMD,NTMAX)

100   CONTINUE
C ------- Write the rotated EOFs
      OPEN(UNIT=7, FILE=FILNM3, STATUS='UNKNOWN')
      DO 105 I=1,NMD
      DO 105 L=1,NY
        LL = (L-1)*NX
        WRITE(7,FORMT3)  (RMAT(K+LL,I)*SCL, K=1,NX)
105   CONTINUE

      IF (IPC.EQ.0)  GO TO 200
C ------- Write the eigenvalues
      OPEN(UNIT=8, FILE='inform.d', STATUS='UNKNOWN')
      WRITE(8,110)  TVAR*FLOAT(NST)
110   FORMAT(5X,'TOTAL VARIANCE = ',E15.7,///)
      SUM = 0.0
      DO IM=1,NMD
        VAR = EGV(IM)/TVAR
        SUM = SUM + VAR
        WRITE(8,115)  VAR, SUM
115     FORMAT(5X,'VARIANCE AND CUMULATIVE VARIANCE = ',2E16.7,/)
      END DO

C ------- Write the PC time series
      OPEN(UNIT=9, FILE='pcts.d', STATUS='UNKNOWN')
      DO IM=1,NMD
        WRITE(9,120)  (PCT(I,IM), I=1,NCNT)
120     FORMAT(6E13.5)
      END DO

200   STOP
      END

      SUBROUTINE TVARMX(M,K,A,NC,TV,H,KRAW)

C        Performs the varimax rotation (Revision of VARMX from IBM SSP).
C
C          M:     number of attributes; input/output
C          K:     number of factors; input/output
C          A:     factor matrix to be rotated on input
C                    on output contains the rotated factor matrix
C          NC:    the number of rows in main of array A; input
C                    on output it containsthe number of iterations
C          TV:    a scratch vector with at least 51 elements
C                    on output contains variances in varimax solution
C          H:     vector containing communalities on output
C                    must have at least M elements
C          KRAW:  an integer concerning case of varimax solution
C                    = 0  for unweighted (raw) solution
C                    = 1  for normalized solution


      DIMENSION A(NC*K), H(NC), TV(51)


      MAX = NC
      RAW = KRAW
      MC = NC
C ------- Initialization
      EPS = 0.00116
      TVLT = 0.0
      LL = K-1
      NV = 1
      MC = 0
      FN = M
      FFN = FN*FN
      CONS = 0.7071066

C ------- Calculate original communalities
      DO 10 I=1,M
        H(I) = 0.0
      DO 10 J=1,K
        L = MAX*(J-1)+I
        H(I) = H(I) + A(L)*A(L)
10    CONTINUE

C ------- Calculate normalized factor matrix
      DO 20 I=1,M
C ------ Kim 02/24/99 (extraordinary case encountered)
C        H(I) = SQRT(H(I))
        H(I) = SQRT(AMAX1(H(I),1.E-8))
      DO 20 J=1,K
        L = MAX*(J-1)+I
        A(L) = A(L)/(RAW*H(I)+1.-RAW)
20    CONTINUE
      GO TO 35

C ------- Calculate variance for factor matrix
30    NV = NV+1
      TVLT = TV(NV-1)
35    TV(NV) = 0.0
      DO J=1,K
        AA = 0.0
        BB = 0.0
        LB = MAX*(J-1)
        DO 40 I=1,M
          L = LB+I
          CC = A(L)*A(L)
          AA = AA + CC
          BB = BB + CC*CC
40      CONTINUE
        TV(NV) = TV(NV) + (FN*BB-AA*AA)/FFN
      END DO

      IF (NV.GE.51)  GO TO 300
C ------- Perform convergence test
      IF ( (TV(NV)-TVLT)-(1.E-7) ) 110,110,120
110   MC = MC+1
      IF (MC-3) 120,120,300

C -------> rotation of two factors continues up to the statement 200
120   DO 200 J=1,LL
        L1 = MAX*(J-1)
        II = J+1
C ------- Calculate NUM and DEN
      DO 200 K1=II,K
        L2 = MAX*(K1-1)
        AA = 0.0
        BB = 0.0
        CC = 0.0
        DD = 0.0
        DO I=1,M
          L3 = L1+I
          L4 = L2+I
          U = (A(L3)+A(L4))*(A(L3)-A(L4))
          T = A(L3)*A(L4)
          T = T+T
          CC = CC+(U+T)*(U-T)
          DD = DD+2.0*U*T
          AA = AA+U
          BB = BB+T
        END DO
        T = DD - 2.0*AA*BB/FN
        B = CC - (AA*AA-BB*BB)/FN

C ------- Comparison of NUM and DEN
        IF (T-B) 135,130,150
130     IF ((T+B) .LT. EPS)  GO TO 200
C ---------- NUM + DEN is greater than or equal to the tolerance factor
        COS4T = CONS
        SIN4T = CONS
        GO TO 160
C ---------- NUM is less than DEN
135     TAN4T = ABS(T)/ABS(B)
        IF (TAN4T .LT. EPS)  GO TO 140
        COS4T = 1.0 / SQRT(1.0+TAN4T*TAN4T)
        SIN4T = TAN4T*COS4T
        GO TO 160
140     IF (B.GE.0.)  GO TO 200
        SINP = CONS
        COSP = CONS
        GO TO 180
C ---------- NUM is greater than DEN
150     CTN4T = ABS(T/B)
        IF (CTN4T .LT. EPS)  GO TO 155
        SIN4T = 1.0 / SQRT(1.0+CTN4T*CTN4T)
        COS4T = CTN4T*SIN4T
        GO TO 160
155     COS4T = 0.0
        SIN4T = 1.0

C ------- Determine COS(THETA) and SIN(THETA)
160     COS2T = SQRT((1.0+COS4T)/2.0)
        SIN2T = SIN4T/(2.0*COS2T)
        COST = SQRT((1.0+COS2T)/2.0)
        SINT = SIN2T/(2.0*COST)

C ------- Determine COS(PHI) and SIN(PHI)
        IF (B.LE.0.)  GO TO 170
        COSP = COST
        SINP = SINT
        GO TO 175
170     COSP = CONS*COST + CONS*SINT
        SINP = ABS(CONS*COST-CONS*SINT)
175     IF (T.GT.0.)  GO TO 180
        SINP = -SINP

C ---------- perform rotation
180     DO I=1,M
          L3 = L1+I
          L4 = L2+I
          AA = A(L3)*COSP+A(L4)*SINP
          A(L4) = -A(L3)*SINP+A(L4)*COSP
          A(L3) = AA
        END DO
200   CONTINUE
      GO TO 30

C ------- Denormalize varimax loadings
300   DO 320 J=1,K
        KEEP = MAX*(J-1)
        CSUM = 0.
        DO 305 I=1,M
          L = KEEP + I
          A(L) = A(L)*(RAW*H(I)+1.-RAW)
          CSUM = CSUM + A(L)
305     CONTINUE
C -------> Check for negative column sums -- If found, perform 180 deg rotation
        IF (CSUM .GE. 0.)  GO TO 320
        DO 310 I=1,M
          L = KEEP + I
          A(L) = -A(L)
310     CONTINUE
320   CONTINUE

C ------- Reconstitute communalities
      MC = NV-1
      DO 330 I=1,M
        H(I) = H(I)*H(I)
330   CONTINUE

      RETURN
      END

      SUBROUTINE EIGSRT(D,V,PC,NP,NM,NT)

      DIMENSION D(NM), V(NP,NM), PC(NT,NM)

      DO 20 I=1,NM-1
        K = I
        P = D(I)
        DO 5 J=I+1,NM
          IF (D(J).GE.P) THEN
            K = J
            P = D(J)
          END IF
5       CONTINUE
        IF (K.NE.I) THEN
          D(K) = D(I)
          D(I) = P
          DO 10 J=1,NP
            P = V(J,I)
            V(J,I) = V(J,K)
            V(J,K) = P
10        CONTINUE
          DO 15 J=1,NT
            P = PC(J,I)
            PC(J,I) = PC(J,K)
            PC(J,K) = P
15        CONTINUE
        END IF
20    CONTINUE

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

      PARAMETER (NMAX=615, TINY=1.0E-20)
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
