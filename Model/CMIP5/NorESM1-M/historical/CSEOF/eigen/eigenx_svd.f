C *****************************************************************************
C
C        PROGRAM NAME : eigenx.f
C        PROGRAMMER : Dr. Kwang Y. Kim
C        CODE IDENTIFICATION NUMBER = EIGENX/VERSION 1.0
C        CODE CLASSIFICATION = Scientific Computer Code
C        CREATION DATE = January 07, 1990
C        REVISION DATE = not revised
C        REVISION INFORMATION = not applicable
C
C *****************************************************************************


      ALLOCATABLE DMTRX(:,:), ARR(:), AVG(:,:)
      ALLOCATABLE COV(:,:), EIGV(:), EIGF(:,:), PCTS(:)
      ALLOCATABLE TSER(:,:), EOF(:)
      CHARACTER*50 FILNM1, FILNM2, FILNM3, FORMT1, FORMT2, FORMT3
      CHARACTER*3 ANS, SEQ(10)

      DATA SEQ / '0TH', '1ST', '2ND', '3RD', '4TH',
     &           '5TH', '6TH', '7TH', '8TH', '9TH' /


      PI = 4.0*ATAN(1.0)
      TPI = 2.0*PI

      PRINT *, '  Type the name and format of the time series data.'
      READ 5, FILNM1, FORMT1
      PRINT *, '  Type the size (NX x NY) of station array.'
      READ *, NX, NY
      NST = NX*NY
      PRINT *, '  Type the number of samples at each station.'
      READ *, NPTS
      FNPTS = FLOAT(NPTS)
      SNPTS = SQRT(FNPTS)
      PRINT *, '  First index of the input array: (1: time, 2: space).'
      READ *, IRD
      PRINT *, '  Type smoothing option (0:No, M:Moving average lag).'
      READ *, LAG
      PRINT *, '  Type cycle period for detrending (0: No).'
      READ *, IDTR
      PRINT *, '  Type the area adjustment option (0: No, 1: Yes).'
      READ *, IARA
      PRINT *, '  Type the starting latitude and increment.'
      READ *, SLAT, DLAT
      PRINT *, '  Type percent variance to be achieved.'
      READ *, PVAR
      PRINT *, '  Type pattern scaling factor.'
      READ *, SCL
      PRINT *, '  type the number of eigenfunctions to be printed.'
      READ *, NMD
      PRINT *, '  Type the PC normalization option (0: no, 1: yes)'
      READ *, INORM
      PRINT *, '  Type the name and format of the EOF output file.'
      READ 5, FILNM2, FORMT2
      PRINT *, '  Type the name and format of the PC time series file.'
      READ 5, FILNM3, FORMT3
5     FORMAT(A50)


C ------- Allocate space for dynamic arrays
      ALLOCATE(DMTRX(NPTS,NPTS))
      ALLOCATE(ARR(NPTS))
      ALLOCATE(AVG(IDTR,NST))
      ALLOCATE(COV(NPTS,NPTS))
      ALLOCATE(EIGV(NPTS))
      ALLOCATE(EIGF(NPTS,NPTS))
      ALLOCATE(PCTS(NPTS))
      ALLOCATE(TSER(NPTS,NST))
      ALLOCATE(EOF(NST))


C ------- Read Dataset
      SELECT CASE (FORMT1)
      CASE ('DIR')
        IF (IRD.EQ.1) THEN
          OPEN(UNIT=11, FILE=FILNM1, STATUS='OLD',
     &         ACCESS='DIRECT', RECL=NPTS*4)
          DO L=1,NST
            READ(11,REC=L)  (TSER(I,L), I=1,NPTS)
          END DO
        ELSE
          OPEN(UNIT=11, FILE=FILNM1, STATUS='OLD',
     &         ACCESS='DIRECT', RECL=NST*4)
          DO I=1,NPTS
            READ(11,REC=I)  (TSER(I,L), L=1,NST)
          END DO
        END IF
      CASE ('SEQ')
        OPEN(UNIT=11, FILE=FILNM1, STATUS='OLD', FORM='UNFORMATTED')
        IF (IRD.EQ.1) THEN
          DO L=1,NST
            READ(11)  (TSER(I,L), I=1,NPTS)
          END DO
        ELSE
          DO I=1,NPTS
            READ(11)  (TSER(I,L), L=1,NST)
          END DO
        END IF
      CASE DEFAULT
        OPEN(UNIT=11, FILE=FILNM1, STATUS='OLD')
        IF (IRD.EQ.1) THEN
          DO L=1,NST
            READ(11,FORMT1)  (TSER(I,L), I=1,NPTS)
          END DO
        ELSE
          DO I=1,NPTS
          DO L=1,NY
            READ(11,FORMT1)  (TSER(I,K+(L-1)*NX), K=1,NX)
          END DO
          END DO
        END IF
      END SELECT

C ------- Smoothing
      DO 10 J=1,NST
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

C ------- Subtract (Monthly) Mean
      IF (IDTR.NE.0) THEN
        OPEN(UNIT=9, FILE='avg.d', STATUS='UNKNOWN')
        DO J=1,NST
          DO 15 IM=1,IDTR
            SUM = 0.0
            DO K=IM,NPTS,IDTR
              SUM = SUM + TSER(K,J)/FLOAT((NPTS-IM)/IDTR+1)
            END DO
            AVG(IM,J) = SUM
            DO K=IM,NPTS,IDTR
              TSER(K,J) = TSER(K,J) - AVG(IM,J)
            END DO
15        CONTINUE
        END DO
        DO IM=1,IDTR
          WRITE(9,'(6E13.5)')  (AVG(IM,J), J=1,NST)
        END DO
      END IF

C ------- Area adjustment
      IF (IARA.EQ.1) THEN
        TAREA = 0.0
        DO J=1,NY
          JJ = J-1
          ALAT = SLAT + FLOAT(JJ)*DLAT
          TAREA = TAREA + COS(ALAT*PI/180.)
        END DO
        TAREA = TAREA*FLOAT(NX)

        DO K=1,NPTS
        DO J=1,NY
          JJ = J-1
          ALAT = SLAT + FLOAT(JJ)*DLAT
          SCLF = SQRT(COS(ALAT*PI/180.)/TAREA)*SQRT(FLOAT(NST))
          DO I=1,NX
            IJ = (J-1)*NX+I
            TSER(K,IJ) = TSER(K,IJ)*SCLF
          END DO
        END DO
        END DO
      END IF

C ------- Covariance Matrix
      DO 20 J=1,NPTS
      DO 20 I=1,NPTS
        SUMC = 0.0
        DO K=1,NST
          SUMC = SUMC + TSER(I,K)*TSER(J,K)
        END DO
        SUMC = SUMC/FLOAT(NST)
        COV(I,J) = SUMC
20    CONTINUE

C ------- Total Variance
      OPEN(UNIT=7, FILE='inform.d', STATUS='UNKNOWN')
      TVAR = 0.0
      DO 25 I=1,NPTS
        TVAR = TVAR + COV(I,I)
25    CONTINUE
      TVAR = TVAR/FNPTS
      WRITE(7,30) TVAR
30    FORMAT(5X,'TOTAL VARIANCE = ',E15.7,///)

C ------- Call Eigenvalue Routines
      CALL SVDCMP(COV,NPTS,NPTS,NPTS,NPTS,EIGV,EIGF)
      CALL EIGSRT(EIGV,EIGF,NPTS,NPTS)

C ------- Write Eigenmodes (PCs) and Modal Contributions
      SELECT CASE (FORMT3)
      CASE ('DIR')
        OPEN(UNIT=22, FILE=FILNM3, STATUS='UNKNOWN',
     &       ACCESS='DIRECT', RECL=NPTS*4)
      CASE ('SEQ')
        OPEN(UNIT=22, FILE=FILNM3, STATUS='UNKNOWN',
     &       FORM='UNFORMATTED')
      CASE DEFAULT
        OPEN(UNIT=22, FILE=FILNM3, STATUS='UNKNOWN')
      END SELECT

      SUM = 0.0
      DO 45 I=1,NPTS
        EIGV(I) = EIGV(I)/FNPTS
        STD = SQRT(EIGV(I))
        VAR = EIGV(I)/TVAR
        SUM = SUM + VAR
        WRITE(7,35) VAR, SUM
35      FORMAT(5X,'VARIANCE AND CUMULATIVE VARIANCE = ',2E16.7,/)
        DO J=1,NPTS
          IF (INORM.EQ.1) THEN
            PCTS(J) = EIGF(J,I)*SNPTS
          ELSE
            PCTS(J) = EIGF(J,I)*SNPTS*STD
          END IF
        END DO
        NMODE = I

        SELECT CASE (FORMT3)
        CASE ('DIR')
          WRITE(22,REC=I)  (PCTS(J), J=1,NPTS)
        CASE ('SEQ')
          WRITE(22)  (PCTS(J), J=1,NPTS)
        CASE DEFAULT
          WRITE(22,FORMT3)  (PCTS(J), J=1,NPTS)
        END SELECT

        IF (SUM*100. .GE. PVAR)  GO TO 50
45    CONTINUE
50    PRINT *, NMODE

C ------- Write EOF patterns
      SELECT CASE (FORMT2)
      CASE ('DIR')
        OPEN(UNIT=21, FILE=FILNM2, STATUS='UNKNOWN',
     &       ACCESS='DIRECT', RECL=NST*4)
      CASE ('SEQ')
        OPEN(UNIT=21, FILE=FILNM2, STATUS='UNKNOWN',
     &       FORM='UNFORMATTED')
      CASE DEFAULT
        OPEN(UNIT=21, FILE=FILNM2, STATUS='UNKNOWN')
      END SELECT

      DO IM=1,MIN(NMODE,NMD)
C ------- EOF Pattern
        DO J=1,NST
          SUM = 0.0
          DO K=1,NPTS
            SUM = SUM + TSER(K,J)*EIGF(K,IM)
          END DO
          EOF(J) = SUM/SNPTS
        END DO

C ------- Normalization
C        SUM = 0.0
C        DO J=1,NST
C          SUM = SUM + EOF(J)*EOF(J)
C        END DO
C        STD = SQRT(SUM/FLOAT(NST))
C        STD = SQRT(SUM)
C        STD = SQRT(EIGV(IM))
        DO L=1,NY
          IF (IARA.EQ.1) THEN
            LL = L-1
            ALAT = SLAT + FLOAT(LL)*DLAT
            SCLF = SQRT(COS(ALAT*PI/180.)/TAREA)*SQRT(FLOAT(NST))
          ELSE
            SCLF = 1.0
          END IF
        DO K=1,NX
          J = (L-1)*NX+K
          IF (INORM.EQ.1) THEN
            EOF(J) = EOF(J)/SCLF
          ELSE
            EOF(J) = EOF(J)/SCLF/SQRT(EIGV(IM))
          END IF
        END DO
        END DO

        SELECT CASE (FORMT2)
        CASE ('DIR')
          WRITE(21,REC=IM)  (EOF(J)*SCL, J=1,NST)
        CASE ('SEQ')
          WRITE(21)  (EOF(J)*SCL, J=1,NST)
        CASE DEFAULT
          DO L=1,NY
            LL = (L-1)*NX
            WRITE(21,FORMT2) (EOF(K+LL)*SCL, K=1,NX)
          END DO
        END SELECT
      END DO

      STOP
      END

      SUBROUTINE EIGSRT(D,V,N,NP)

      DIMENSION D(NP),V(NP,NP)


      DO 13 I=1,N-1
        K=I
        P=D(I)
        DO 11 J=I+1,N
          IF(D(J).GE.P)THEN
            K=J
            P=D(J)
          ENDIF
11      CONTINUE
        IF(K.NE.I)THEN
          D(K)=D(I)
          D(I)=P
          DO 12 J=1,N
            P=V(J,I)
            V(J,I)=V(J,K)
            V(J,K)=P
12        CONTINUE
        ENDIF
13    CONTINUE
      RETURN
      END

      SUBROUTINE SVDCMP(A,MM,NN,MP,NP,W,V)

C       Given a matrix A, with logical dimensions M by N and physical 
C       dimensions MP by NP, this routine computes its singular value 
C       decomposition, A = U*W*V'.  The matrix U replaces A on output.
C       This diagonal matrix of singular values W is output as a vector 
C       W.  The matrix V (not the transpose V') is output as V.  M must 
C       be greater or equal to N; if it is smaller, then A should be 
C       filled up to square with zero rows.

      PARAMETER (NMAX=5000)
      DIMENSION A(MP,NP), W(NP), V(NP,NP), RV1(NMAX)

      G = 0.0
      SCALE = 0.0
      ANORM = 0.0
      DO 25 I=1,NN
        L = I+1
        RV1(I) = SCALE*G
        G = 0.0
        S = 0.0
        SCALE = 0.0
        IF (I.LE.MM) THEN
          DO 11 K=I,MM
            SCALE = SCALE + ABS(A(K,I))
11        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 12 K=I,MM
              A(K,I) = A(K,I)/SCALE
              S = S + A(K,I)*A(K,I)
12          CONTINUE
            F = A(I,I)
            G = -SIGN(SQRT(S),F)
            H = F*G - S
            A(I,I) = F-G
            IF (I.NE.NN) THEN
              DO 15 J=L,NN
                S = 0.0
                DO 13 K=I,MM
                  S = S + A(K,I)*A(K,J)
13              CONTINUE
                F = S/H
                DO 14 K=I,MM
                  A(K,J) = A(K,J) + F*A(K,I)
14              CONTINUE
15            CONTINUE
            ENDIF
            DO 16 K=I,MM
              A(K,I) = SCALE*A(K,I)
16          CONTINUE
          ENDIF
        ENDIF
        W(I) = SCALE*G
        G = 0.0
        S = 0.0
        SCALE = 0.0
        IF ((I.LE.MM).AND.(I.NE.NN)) THEN
          DO 17 K=L,NN
            SCALE = SCALE + ABS(A(I,K))
17        CONTINUE
          IF (SCALE.NE.0.0) THEN
            DO 18 K=L,NN
              A(I,K) = A(I,K)/SCALE
              S = S + A(I,K)*A(I,K)
18          CONTINUE
            F = A(I,L)
            G = -SIGN(SQRT(S),F)
            H = F*G - S
            A(I,L) = F-G
            DO 19 K=L,NN
              RV1(K) = A(I,K)/H
19          CONTINUE
            IF (I.NE.MM) THEN
              DO 23 J=L,MM
                S = 0.0
                DO 21 K=L,NN
                  S = S + A(J,K)*A(I,K)
21              CONTINUE
                DO 22 K=L,NN
                  A(J,K) = A(J,K) + S*RV1(K)
22              CONTINUE
23            CONTINUE
            ENDIF
            DO 24 K=L,NN
              A(I,K) = SCALE*A(I,K)
24          CONTINUE
          ENDIF
        ENDIF
        ANORM = MAX(ANORM,(ABS(W(I))+ABS(RV1(I))))
25    CONTINUE

      DO 32 I=NN,1,-1
        IF (I.LT.NN) THEN
          IF (G.NE.0.0) THEN
            DO 26 J=L,NN
              V(J,I) = (A(I,J)/A(I,L))/G
26          CONTINUE
            DO 29 J=L,NN
              S = 0.0
              DO 27 K=L,NN
                S = S + A(I,K)*V(K,J)
27            CONTINUE
              DO 28 K=L,NN
                V(K,J) = V(K,J) + S*V(K,I)
28            CONTINUE
29          CONTINUE
          ENDIF
          DO 31 J=L,NN
            V(I,J) = 0.0
            V(J,I) = 0.0
31        CONTINUE
        ENDIF
        V(I,I) = 1.0
        G = RV1(I)
        L = I
32    CONTINUE

      DO 39 I=NN,1,-1
        L = I+1
        G = W(I)
        IF (I.LT.NN) THEN
          DO 33 J=L,NN
            A(I,J) = 0.0
33        CONTINUE
        ENDIF
        IF (G.NE.0.0) THEN
          G = 1.0/G
          IF (I.NE.NN) THEN
            DO 36 J=L,NN
              S = 0.0
              DO 34 K=L,MM
                S = S + A(K,I)*A(K,J)
34            CONTINUE
              F = (S/A(I,I))*G
              DO 35 K=I,MM
                A(K,J) = A(K,J) + F*A(K,I)
35            CONTINUE
36          CONTINUE
          ENDIF
          DO 37 J=I,MM
            A(J,I) = A(J,I)*G
37        CONTINUE
        ELSE
          DO 38 J=I,MM
            A(J,I) = 0.0
38        CONTINUE
        ENDIF
        A(I,I) = A(I,I) + 1.0
39    CONTINUE

      DO 49 K=NN,1,-1
        DO 48 ITS=1,30
          DO 41 L=K,1,-1
            NM = L-1
            IF ((ABS(RV1(L))+ANORM).EQ.ANORM)  GO TO 2
            IF ((ABS(W(NM))+ANORM).EQ.ANORM)  GO TO 1
41        CONTINUE
1         C = 0.0
          S = 1.0
          DO 43 I=L,K
            F = S*RV1(I)
            IF ((ABS(F)+ANORM).NE.ANORM) THEN
              G = W(I)
              H = SQRT(F*F+G*G)
              W(I) = H
              H = 1.0/H
              C = (G*H)
              S = -(F*H)
              DO 42 J=1,MM
                Y = A(J,NM)
                Z = A(J,I)
                A(J,NM) = (Y*C)+(Z*S)
                A(J,I) = -(Y*S)+(Z*C)
42            CONTINUE
            ENDIF
43        CONTINUE
2         Z = W(K)
          IF (L.EQ.K) THEN
            IF (Z.LT.0.0) THEN
              W(K) = -Z
              DO 44 J=1,NN
                V(J,K) = -V(J,K)
44            CONTINUE
            END IF
            GO TO 3
          END IF
          IF (ITS.EQ.100) STOP 'No convergence in 100 iterations'
          X = W(L)
          NM = K-1
          Y = W(NM)
          G = RV1(NM)
          H = RV1(K)
          F = ((Y-Z)*(Y+Z)+(G-H)*(G+H))/(2.0*H*Y)
          G = SQRT(F*F+1.0)
          F = ((X-Z)*(X+Z)+H*((Y/(F+SIGN(G,F)))-H))/X
          C = 1.0
          S = 1.0
          DO 47 J=L,NM
            I = J+1
            G = RV1(I)
            Y = W(I)
            H = S*G
            G = C*G
            Z = SQRT(F*F + H*H)
            RV1(J) = Z
            C = F/Z
            S = H/Z
            F = (X*C) + (G*S)
            G = -(X*S) + (G*C)
            H = Y*S
            Y = Y*C
            DO 45 NM=1,NN
              X = V(NM,J)
              Z = V(NM,I)
              V(NM,J) = (X*C) + (Z*S)
              V(NM,I) = -(X*S) + (Z*C)
45          CONTINUE
            Z = SQRT(F*F + H*H)
            W(J) = Z
            IF (Z.NE.0.0) THEN
              Z = 1.0/Z
              C = F*Z
              S = H*Z
            END IF
            F = (C*G) + (S*Y)
            X = -(S*G) + (C*Y)
            DO 46 NM=1,MM
              Y = A(NM,J)
              Z = A(NM,I)
              A(NM,J) = (Y*C) + (Z*S)
              A(NM,I) = -(Y*S) + (Z*C)
46          CONTINUE
47        CONTINUE
          RV1(L) = 0.0
          RV1(K) = F
          W(K) = X
48      CONTINUE
3       CONTINUE
49    CONTINUE

      RETURN
      END
