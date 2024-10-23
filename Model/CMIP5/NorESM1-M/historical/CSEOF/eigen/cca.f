C *****************************************************************************
C
C       PROGRAM NAME : cca.f
C       PROGRAMMER : Dr. Kwang Y. Kim
C       CODE IDENTIFICATION : CCA/VERSION 1.0
C       CODE CLASSIFICATION : scientific computer code
C       CREATION DATE : October 31, 1997
C       REVISION DATE : October 31, 2002
C       REVISION INFORMATION : f90 upgrade (October 31, 2002)
C
C *****************************************************************************


C     This program computes canonical correlation patterns of two vector
C     time series:
C
C          X = sum(k) a_k(t) p_k        Y = sum(k) b_k(t) q_k
C
C     where p_k and q_k are eigenvectors (canonical patterns) and a_k(t)
C     and b_k(t) are the corresponding PC time series.  The eigenvectors
C     are obtained from their adjoint patterns.  The adjoint patterns are
C     the eigenvectors of
C
C          A = (S_x)^I S_xy (S_y)^I (S_xy)^T
C
C          B = (S_y)^I (S_xy)^T (S_x)^I S_xy
C
C     where S_x, S_y are the covariance matrices of X and Y and S_xy are
C     is the cross-covariance matrix.  The superscripts I and T represent 
C     inverse and transpose, respectively.
C
C     Once adjoint patterns are computed the eigenvectors are given by
C
C          p_k = S_x (p_k)^A            q_k = S_y (q_k)^A
C
C     where the superscript A indicates adjoint pattern.  Finally, the PC
C     time series are computed as
C
C          a_k(t) = < (p_k)^A X >       b_k(t) = < (q_k)^A Y >
C
C     where < > represents inner product of two vectors.


      ALLOCATABLE TSER1(:,:), TSER2(:,:), SER1(:,:), SER2(:,:)
      ALLOCATABLE ARR(:), AVG(:), AVG1(:), AVG2(:), VAR1(:), VAR2(:)
      ALLOCATABLE SIGX(:,:), SIGY(:,:), SIGXY(:,:),
     &            SIGXI(:,:), SIGYI(:,:)
      ALLOCATABLE COVA(:,:), COVB(:,:), COV1(:,:), COV2(:,:)
      ALLOCATABLE INDX(:)
      ALLOCATABLE EGF1(:,:), EGV1(:), EGF2(:,:), EGV2(:)
      ALLOCATABLE PCT1(:,:), PCT2(:,:)
      CHARACTER*50 FILNM1, FORMT1, FILNM2, FORMT2

      PI = 4.0*ATAN(1.0)
      TPI = 2.0*PI
      RTD = 180./PI


5     FORMAT(A50)
      PRINT *, '  Type name and format of the first input file.'
      READ 5, FILNM1, FORMT1
      PRINT *, '  First index of the input array (1: time,  2: space).'
      READ *, IRD
      PRINT *, '  Type the dimension (NX,NY) of sampling stations.'
      READ *, NX, NY
      NST = NX*NY
      PRINT *, '  Type name and format of the second input file.'
      READ 5, FILNM2, FORMT2
      PRINT *, '  First index of the input array (1: time,  2: space).'
      READ *, JRD
      PRINT *, '  Type the dimension (MX,MY) of sampling stations.'
      READ *, MX, MY
      MST = MX*MY
      PRINT *, '  Type the length of time series.'
      READ *, NTOT
      PRINT *, '  Type the number of samples for CCA calculation.'
      READ *, NPTS
      PRINT *, '  Type smoothing option (0:No, M:Moving average lag).'
      READ *, LAG
      PRINT *, '  Type the cycle period for detrending  (0: No).'
      READ *, IDTR
      PRINT *, '  Type percent variance to be retained.'
      READ *, PVAR
      PRINT *, '  Type maximum lead time for the second time series.'
      READ *, MLEAD
      NPTS = MIN(NPTS,NTOT-MLEAD)

C ------- Allocate Dynamic Arrays
      ALLOCATE(TSER1(NTOT,NST))
      ALLOCATE(TSER2(NTOT,MST))
      ALLOCATE(SER1(NTOT,NST))
      ALLOCATE(SER2(NTOT,MST))
      ALLOCATE(ARR(NTOT))
      ALLOCATE(AVG(IDTR))
      ALLOCATE(AVG1(NST))
      ALLOCATE(AVG2(MST))
      ALLOCATE(VAR1(NST))
      ALLOCATE(VAR2(MST))
      ALLOCATE(SIGX(NST,NST))
      ALLOCATE(SIGY(MST,MST))
      ALLOCATE(SIGXY(NST,MST))
      ALLOCATE(SIGXI(NST,NST))
      ALLOCATE(SIGYI(MST,MST))
      ALLOCATE(COVA(NST,MST))
      ALLOCATE(COVB(MST,NST))
      ALLOCATE(COV1(NST,NST))
      ALLOCATE(COV2(MST,MST))
      ALLOCATE(INDX(MAX(NST,MST)))
      ALLOCATE(EGF1(NST,NST))
      ALLOCATE(EGV1(NST))
      ALLOCATE(EGF2(MST,MST))
      ALLOCATE(EGV2(MST))
      ALLOCATE(PCT1(NTOT,NST))
      ALLOCATE(PCT2(NTOT,MST))


C ------- Read Input File
      SELECT CASE (FORMT1)
      CASE('DIR')
        IF (IRD.EQ.1) THEN
          OPEN(UNIT=4, FILE=FILNM1, STATUS='OLD',
     &         ACCESS='DIRECT', RECL=NTOT)
          DO K=1,NST
            READ(4,REC=K)  (TSER1(I,K), I=1,NTOT)
          END DO
        ELSE
          OPEN(UNIT=4, FILE=FILNM1, STATUS='OLD',
     &         ACCESS='DIRECT', RECL=NST)
          DO I=1,NTOT
            READ(4,REC=L)  (TSER1(I,K), K=1,NST)
          END DO
        END IF
      CASE ('SEQ')
        IF (IRD.EQ.1) THEN
          OPEN(UNIT=4, FILE=FILNM1, STATUS='OLD', FORM='UNFORMATTED')
          DO K=1,NST
            READ(4)  (TSER1(I,K), I=1,NTOT)
          END DO
        ELSE
          OPEN(UNIT=4, FILE=FILNM1, STATUS='OLD', FORM='UNFORMATTED')
          DO I=1,NTOT
            READ(4)  (TSER1(I,K), K=1,NST)
          END DO
        END IF
      CASE DEFAULT
        OPEN(UNIT=4, FILE=FILNM1, STATUS='OLD')
        IF (IRD.EQ.1) THEN
          DO K=1,NST
            READ(4,FORMT1)  (TSER1(I,K), I=1,NTOT)
          END DO
        ELSE
          DO I=1,NTOT
          DO L=1,NY
            READ(4,FORMT1)  (TSER1(I,K+(L-1)*NX), K=1,NX)
          END DO
          END DO
        END IF
      END SELECT
      CLOSE(UNIT=4)

      SELECT CASE (FORMT2)
      CASE('DIR')
        IF (JRD.EQ.1) THEN
          OPEN(UNIT=4, FILE=FILNM2, STATUS='OLD',
     &         ACCESS='DIRECT', RECL=NTOT)
          DO K=1,MST
            READ(4,REC=K)  (TSER2(I,K), I=1,NTOT)
          END DO
        ELSE
          OPEN(UNIT=4, FILE=FILNM2, STATUS='OLD',
     &         ACCESS='DIRECT', RECL=MST)
          DO I=1,NTOT
            READ(4,REC=L)  (TSER2(I,K), K=1,MST)
          END DO
        END IF
      CASE ('SEQ')
        IF (JRD.EQ.1) THEN
          OPEN(UNIT=4, FILE=FILNM2, STATUS='OLD', FORM='UNFORMATTED')
          DO K=1,MST
            READ(4)  (TSER2(I,K), I=1,NTOT)
          END DO
        ELSE
          OPEN(UNIT=4, FILE=FILNM2, STATUS='OLD', FORM='UNFORMATTED')
          DO I=1,NTOT
            READ(4)  (TSER2(I,K), K=1,MST)
          END DO
        END IF
      CASE DEFAULT
        OPEN(UNIT=4, FILE=FILNM2, STATUS='OLD')
        IF (JRD.EQ.1) THEN
          DO K=1,MST
            READ(4,FORMT2)  (TSER2(I,K), I=1,NTOT)
          END DO
        ELSE
          DO I=1,NTOT
          DO L=1,MY
            READ(4,FORMT2)  (TSER2(I,K+(L-1)*MX), K=1,MX)
          END DO
          END DO
        END IF
      END SELECT
      CLOSE(UNIT=4)

C ------- Open Output Files
      OPEN(UNIT=7, FILE='inform.d', STATUS='UNKNOWN')
      OPEN(UNIT=8, FILE='amode.d', STATUS='UNKNOWN')
      OPEN(UNIT=9, FILE='pct1.d', STATUS='UNKNOWN')
      OPEN(UNIT=10, FILE='pct2.d', STATUS='UNKNOWN')
      OPEN(UNIT=11, FILE='emode.d', STATUS='UNKNOWN')

C ------- Smoothing
      DO J=1,NST
        DO I=1,NTOT
          KS = MAX(I-LAG,1)
          KE = MIN(I+LAG,NTOT)
          KN = KE-KS+1
          SUM = 0.0
          DO K=KS,KE
            SUM = SUM + TSER1(K,J)
          END DO
          ARR(I) = SUM/FLOAT(KN)
        END DO
        DO I=1,NTOT
          TSER1(I,J) = ARR(I)
        END DO
      END DO

      DO J=1,MST
        DO I=1,NTOT
          KS = MAX(I-LAG,1)
          KE = MIN(I+LAG,NTOT)
          KN = KE-KS+1
          SUM = 0.0
          DO K=KS,KE
            SUM = SUM + TSER2(K,J)
          END DO
          ARR(I) = SUM/FLOAT(KN)
        END DO
        DO I=1,NTOT
          TSER2(I,J) = ARR(I)
        END DO
      END DO


C       { INTRODUCING DIFFERENT LAGS }

      DO 200 LEAD=0,MLEAD

C ------- Remove Mean or Seasonal Cycle
        IF (IDTR.NE.0) THEN
          DO J=1,NST
            DO IM=1,IDTR
              SUM = 0.0
              DO I=IM,NPTS,IDTR
                SUM = SUM + TSER1(I,J)
              END DO
              SUM = SUM/FLOAT((NTOT-IM)/IDTR+1)
              AVG(IM) = SUM
              DO I=IM,NTOT,IDTR
                TSER1(I,J) = TSER1(I,J) - AVG(IM)
              END DO
            END DO
C            WRITE(29,'(6E13.5)')  (AVG(IM), IM=1,IDTR)
          END DO

          DO J=1,MST
            DO IM=1,IDTR
              SUM = 0.0
              DO I=IM,NPTS,IDTR
                SUM = SUM + TSER2(I,J)
              END DO
              SUM = SUM/FLOAT((NTOT-IM)/IDTR+1)
              AVG(IM) = SUM
              DO I=IM,NTOT,IDTR
                TSER2(I,J) = TSER2(I,J) - AVG(IM)
              END DO
            END DO
C            WRITE(29,'(6E13.5)')  (AVG(IM), IM=1,IDTR)
          END DO
        END IF

C ------- Calculate Mean and Variance
        DO J=1,NST
          SUM1 = 0.0
          SUM2 = 0.0
          DO I=1,NPTS
            SUM1 = SUM1 + TSER1(I,J)
            SUM2 = SUM2 + TSER1(I,J)*TSER1(I,J)
          END DO
          AVG1(J) = SUM1/FLOAT(NPTS)
          VAR1(J) = SUM2/FLOAT(NPTS) - AVG1(J)**2
          DO I=1,NTOT
            SER1(I,J) = (TSER1(I,J) - AVG1(J)) / SQRT(VAR1(J))
          END DO
        END DO
C        WRITE(30,'(5E15.7)')  (AVG1(J), J=1,NST)
C        WRITE(30,'(5E15.7)')  (VAR1(J), J=1,NST)

        DO J=1,MST
          SUM1 = 0.0
          SUM2 = 0.0
          DO I=1,NPTS
            SUM1 = SUM1 + TSER2(I+LEAD,J)
            SUM2 = SUM2 + TSER2(I+LEAD,J)*TSER2(I+LEAD,J)
          END DO
          AVG2(J) = SUM1/FLOAT(NPTS)
          VAR2(J) = SUM2/FLOAT(NPTS) - AVG2(J)**2
          DO I=1,NTOT
            SER2(I,J) = (TSER2(I,J) - AVG2(J)) / SQRT(VAR2(J))
          END DO
        END DO
C        WRITE(31,'(5E15.7)')  (AVG2(J), J=1,MST)
C        WRITE(31,'(5E15.7)')  (VAR2(J), J=1,MST)

C ------- Covariance Matrices
        DO 10 J=1,NST
        DO 10 I=J,NST
          SUM = 0.0
          DO K=1,NPTS
            SUM = SUM + SER1(K,I)*SER1(K,J)
          END DO
          SIGX(I,J) = SUM/FLOAT(NPTS)
          SIGX(J,I) = SIGX(I,J)
10      CONTINUE

        DO 15 J=1,MST
        DO 15 I=J,MST
          SUM = 0.0
          DO K=1,NPTS
            SUM = SUM + SER2(K+LEAD,I)*SER2(K+LEAD,J)
          END DO
          SIGY(I,J) = SUM/FLOAT(NPTS)
          SIGY(J,I) = SIGY(I,J)
15      CONTINUE

        DO 20 J=1,MST
        DO 20 I=1,NST
          SUM = 0.0
          DO K=1,NPTS
            SUM = SUM + SER1(K,I)*SER2(K+LEAD,J)
          END DO
          SIGXY(I,J) = SUM/FLOAT(NPTS)
20      CONTINUE

C ------- Eigen Matrix
        CALL LUDCMP(SIGX,NST,NST,INDX,SIGN)
        DO I=1,NST
          DO J=1,NST
            SIGXI(I,J) = 0.0
          END DO
          SIGXI(I,I) = 1.0
        END DO
        DO J=1,NST
          CALL LUBKSB(SIGX,NST,NST,INDX,SIGXI(1,J))
        END DO

        CALL LUDCMP(SIGY,MST,MST,INDX,SIGN)
        DO I=1,MST
          DO J=1,MST
            SIGYI(I,J) = 0.0
          END DO
          SIGYI(I,I) = 1.0
        END DO
        DO J=1,MST
          CALL LUBKSB(SIGY,MST,MST,INDX,SIGYI(1,J))
        END DO

        DO 25 J=1,MST
        DO 25 I=1,NST
          SUM = 0.0
          DO K=1,NST
            SUM = SUM + SIGXI(I,K)*SIGXY(K,J)
          END DO
          COVA(I,J) = SUM
25      CONTINUE

        DO 30 J=1,NST
        DO 30 I=1,MST
          SUM = 0.0
          DO K=1,MST
            SUM = SUM + SIGYI(I,K)*SIGXY(J,K)
          END DO
          COVB(I,J) = SUM
30      CONTINUE

C ------- Adjoint Patterns of X
        DO 35 J=1,NST
        DO 35 I=1,NST
          SUM = 0.0
          DO K=1,MST
            SUM = SUM + COVA(I,K)*COVB(K,J)
          END DO
          COV1(I,J) = SUM
35      CONTINUE

C ------- Total Variance
        TVAR = 0.0
        DO 40 I=1,NST
          TVAR = TVAR + COV1(I,I)
40      CONTINUE
        WRITE(7,45) TVAR
45      FORMAT(///,5X,'TOTAL VARIANCE = ',E15.7,/)

C ------- Call Eigenvalue Routines
        CALL JACOBI(COV1,NST,NST,EGV1,EGF1,NROT)
        CALL EIGSRT(EGV1,EGF1,NST,NST)

C ------- Write Eigenmodes and Modal Contributions
        SUM = 0.0
        DO 60 I=1,NST
          VAR = EGV1(I)/TVAR
          SUM = SUM + VAR
          WRITE(7,50) VAR, SUM
50        FORMAT(5X,'VARIANCE AND CUMULATIVE VARIANCE = ',2E16.7)
          WRITE(8,55) (EGF1(K,I)/SQRT(VAR1(K)), K=1,NST)
55        FORMAT(6E13.5)
          NMD = I
          IF (SUM*100. .GE. PVAR)  GO TO 65
60      CONTINUE
65      PRINT *, '  Number of retained modes = ', NMD

C ------- Adjoint Patterns of Y
        DO 75 J=1,MST
        DO 75 I=1,MST
          SUM = 0.0
          DO K=1,NST
            SUM = SUM + COVB(I,K)*COVA(K,J)
          END DO
          COV2(I,J) = SUM
75      CONTINUE

C ------- Total Variance
        TVAR = 0.0
        DO 80 I=1,MST
          TVAR = TVAR + COV2(I,I)
80      CONTINUE
        WRITE(7,85) TVAR
85      FORMAT(//,5X,'TOTAL VARIANCE = ',E15.7,/)
        WRITE(8,55)

C ------- Call Eigenvalue Routines
        CALL JACOBI(COV2,MST,MST,EGV2,EGF2,NROT)
        CALL EIGSRT(EGV2,EGF2,MST,MST)

C ------- Write Eigenmodes and Modal Contributions
        SUM = 0.0
        DO 90 I=1,MST
          VAR = EGV2(I)/TVAR
          SUM = SUM + VAR
          WRITE(7,50) VAR, SUM
          WRITE(8,55) (EGF2(K,I)/SQRT(VAR2(K)), K=1,MST)
          MMD = I
          IF (SUM*100. .GE. PVAR)  GO TO 95
90      CONTINUE
95      PRINT *, '  Number of retained modes = ', MMD

C ------- PC Time Series
        DO 100 J=1,NMD
          DO I=1,NPTS
            SUM = 0.0
            DO K=1,NST
              SUM = SUM + SER1(I,K)*EGF1(K,J)
            END DO
            PCT1(I,J) = SUM
          END DO
          WRITE(9,55)  (PCT1(I,J), I=1,NPTS)
100     CONTINUE

        DO 105 J=1,MMD
          DO I=1,NPTS
            SUM = 0.0
            DO K=1,MST
              SUM = SUM + SER2(I,K)*EGF2(K,J)
            END DO
            PCT2(I,J) = SUM
          END DO
          WRITE(10,55)  (PCT2(I,J), I=1,NPTS)
105     CONTINUE

C ------- Eigenfunctions
        DO 110 J=1,NST
        DO 110 I=J,NST
          SUM = 0.0
          DO K=1,NPTS
            SUM = SUM + SER1(K,I)*SER1(K,J)
          END DO
          SIGX(I,J) = SUM/FLOAT(NPTS)
          SIGX(J,I) = SIGX(I,J)
110     CONTINUE

        DO 120 J=1,NMD
        DO 120 I=1,NST
          SUM = 0.0
          DO K=1,NST
            SUM = SUM + SIGX(I,K)*EGF1(K,J)
          END DO
          COVA(I,J) = SUM
120     CONTINUE
        DO J=1,NMD
          WRITE(11,55)  (COVA(I,J)*SQRT(VAR1(I)), I=1,NST)
        END DO 

        DO 130 J=1,MST
        DO 130 I=J,MST
          SUM = 0.0
          DO K=1,NPTS
            SUM = SUM + SER2(K+LEAD,I)*SER2(K+LEAD,J)
          END DO
          SIGY(I,J) = SUM/FLOAT(NPTS)
          SIGY(J,I) = SIGY(I,J)
130     CONTINUE

        DO 140 J=1,MMD
        DO 140 I=1,MST
          SUM = 0.0
          DO K=1,MST
            SUM = SUM + SIGY(I,K)*EGF2(K,J)
          END DO
          COVB(I,J) = SUM
140     CONTINUE
        DO J=1,MMD
          WRITE(11,55)  (COVB(I,J)*SQRT(VAR2(I)), I=1,MST)
        END DO

200   CONTINUE

      STOP
      END

      SUBROUTINE JACOBI(A,N,NP,D,V,NROT)

*******************************************************************************
*
*     This subroutine computes all eigenvalues and eigenvectors of a real
*     symmetric matrix A, which is of size N by N, stored in a physical
*     NP by NP array.  On output, elements of A above the diagonal are
*     destroyed.  D returns eigenvalues of A in its first N elements.
*     V is a matrix with the same logical and physical dimensions as A
*     whose columns contain, on output, the normalized eigenvectors of A.
*     NROT returns the number of Jacobi rotation which were required.
*
*******************************************************************************

      PARAMETER (NMX=3000)
      DIMENSION A(NP,NP), D(NP), V(NP,NP), B(NMX), Z(NMX)


      DO 12 IP=1,N
        DO 11 IQ=1,N
          V(IP,IQ)=0.
11      CONTINUE
        V(IP,IP)=1.
12    CONTINUE
      DO 13 IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=0.
13    CONTINUE
      NROT=0
      DO 24 I=1,50
        SM=0.
        DO 15 IP=1,N-1
          DO 14 IQ=IP+1,N
            SM=SM+ABS(A(IP,IQ))
14        CONTINUE
15      CONTINUE
        IF(SM.EQ.0.)RETURN
        IF(I.LT.4)THEN
          TRESH=0.2*SM/N**2
        ELSE
          TRESH=0.
        ENDIF
        DO 22 IP=1,N-1
          DO 21 IQ=IP+1,N
            G=100.*ABS(A(IP,IQ))
            IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP)))
     *         .AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ))))THEN
              A(IP,IQ)=0.
            ELSE IF(ABS(A(IP,IQ)).GT.TRESH)THEN
              H=D(IQ)-D(IP)
              IF(ABS(H)+G.EQ.ABS(H))THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=0.5*H/A(IP,IQ)
                T=1./(ABS(THETA)+SQRT(1.+THETA**2))
                IF(THETA.LT.0.)T=-T
              ENDIF
              C=1./SQRT(1+T**2)
              S=T*C
              TAU=S/(1.+C)
              H=T*A(IP,IQ)
              Z(IP)=Z(IP)-H
              Z(IQ)=Z(IQ)+H
              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H
              A(IP,IQ)=0.
              DO 16 J=1,IP-1
                G=A(J,IP)
                H=A(J,IQ)
                A(J,IP)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
16            CONTINUE
              DO 17 J=IP+1,IQ-1
                G=A(IP,J)
                H=A(J,IQ)
                A(IP,J)=G-S*(H+G*TAU)
                A(J,IQ)=H+S*(G-H*TAU)
17            CONTINUE
              DO 18 J=IQ+1,N
                G=A(IP,J)
                H=A(IQ,J)
                A(IP,J)=G-S*(H+G*TAU)
                A(IQ,J)=H+S*(G-H*TAU)
18            CONTINUE
              DO 19 J=1,N
                G=V(J,IP)
                H=V(J,IQ)
                V(J,IP)=G-S*(H+G*TAU)
                V(J,IQ)=H+S*(G-H*TAU)
19            CONTINUE
              NROT=NROT+1
            ENDIF
21        CONTINUE
22      CONTINUE
        DO 23 IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=0.
23      CONTINUE
24    CONTINUE
      PAUSE '50 iterations should never happen'
      RETURN
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
