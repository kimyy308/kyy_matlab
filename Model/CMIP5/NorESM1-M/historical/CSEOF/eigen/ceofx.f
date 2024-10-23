C *****************************************************************************
C
C       PROGRAM NAME : ceofx.f
C       PROGRAMMER : Dr. Kwang Y. Kim
C       CODE IDENTIFICATION NUMBER = CEOFX/VERSION 1.0
C       CODE CLASSIFICATION = Scientific Computer Code
C       CREATION DATE = July 31, 1997
C       REVISION DATE = October 31, 2002
C       REVISION INFORMATION = f90 version (October 31, 2002)
C
C *****************************************************************************


C     This program computes complex EOFs of a given data:
C
C          C W = lambda W
C
C          C = A + iB          W = U + iV
C
C     where C is an (N x N) Hermitian matrix, W is an eigenvector and lambda 
C     is the corresponding real eigenvalue.  A Hermitian covariance matrix
C     is modified such that the new matrix is symmetric and real.  This is
C     achieved by redefining the problem as
C
C          D X = lambda X
C
C     where the new, real, and symmetric covariance matrix D and the real
C     eigenvectors are written as
C
C          D = ( A -B )        X = ( U )
C              ( B  A )            ( V )
C 
C     Thus the dimension of the problem is (2N x 2N).  Note that (-V, U)
C     is also an eigenvector for an eigenvalue lambda.  These are identical
C     up to the essential phase.


      ALLOCATABLE TSER(:,:), AMP(:,:), PHS(:,:), AVG(:,:)
      ALLOCATABLE THLB(:), ARR(:)
      ALLOCATABLE COV(:,:), V(:,:), D(:)
      ALLOCATABLE INDX(:)
      ALLOCATABLE Y(:), A(:), B(:)
      COMPLEX, DIMENSION(:,:), ALLOCATABLE :: TCMP, DMTRX, W
      COMPLEX, DIMENSION(:), ALLOCATABLE :: PC
      CHARACTER*50 FILNM1, FORMT1, FILNM2, FORMT2
      COMPLEX JIMAG, SUMC
      DATA JIMAG / (0.,1.) /

      PI = 4.0*ATAN(1.0)
      TPI = 2.0*PI
      RTD = 180./PI


5     FORMAT(A50)
      PRINT *, '  Type name and format of the input file.'
      READ 5, FILNM1, FORMT1
      PRINT *, '  Type the dimension (NX,NY) of sampling stations.'
      READ *, NX, NY
      NST = NX*NY
      PRINT *, '  Type the number of samples at each station.'
      READ *, NPTS
      NNPT = 2*NPTS
      PRINT *, '  First index of the input array  (1: time, 2: space).'
      READ *, IRD
      PRINT *, '  Type smoothing option (0:No, M:Moving average lag).'
      READ *, LAG
      PRINT *, '  Type the cycle period for detrending  (0: No).'
      READ *, IDTR
      PRINT *, '  Type percent variance to be achieved.'
      READ *, PVAR
      PRINT *, '  Type CEOF scaling factor.'
      READ *, SCL
      PRINT *, '  Type the number of eigenfunctions to be printed.'
      READ *, NPRT
      PRINT *, '  Type output option.'
      PRINT *
      PRINT *, '    1: PC time series real and imaginary parts'
      PRINT *, '    2: PC time series amplitude and phase'
      PRINT *, '    4: EOF real and imaginary parts'
      PRINT *, '    8: EOF amplitude and phase'
      PRINT *
      READ *, IOUT


C ------- Allocate Dynamic Arrays
      ALLOCATE(TSER(NST,NPTS))
      ALLOCATE(AMP(NPTS,NPTS))
      ALLOCATE(PHS(NPTS,NPTS))
      ALLOCATE(THLB(NPTS))
      ALLOCATE(ARR(NPTS))
      ALLOCATE(AVG(IDTR,NST))
      ALLOCATE(COV(NNPT,NNPT))
      ALLOCATE(V(NNPT,NNPT))
      ALLOCATE(D(NNPT))
      ALLOCATE(INDX(NNPT))
      ALLOCATE(Y(NPTS))
      ALLOCATE(A(0:NPTS))
      ALLOCATE(B(0:NPTS))
      ALLOCATE(TCMP(NST,NPTS))
      ALLOCATE(DMTRX(NPTS,NPTS))
      ALLOCATE(W(NPTS,NPTS))
      ALLOCATE(PC(NST))


C ------- Read Input File
      SELECT CASE (FORMT1)
      CASE('DIR')
        IF (IRD.EQ.1) THEN
          OPEN(UNIT=4, FILE=FILNM1, STATUS='OLD',
     &         ACCESS='DIRECT', RECL=NPTS)
          DO K=1,NST
            READ(4,REC=K)  (TSER(K,I), I=1,NPTS)
          END DO
        ELSE
          OPEN(UNIT=4, FILE=FILNM1, STATUS='OLD',
     &         ACCESS='DIRECT', RECL=NST)
          DO I=1,NPTS
            READ(4,REC=L)  (TSER(K,I), K=1,NST)
          END DO
        END IF
      CASE ('SEQ')
        IF (IRD.EQ.1) THEN
          OPEN(UNIT=4, FILE=FILNM1, STATUS='OLD', FORM='UNFORMATTED')
          DO K=1,NST
            READ(4)  (TSER(K,I), I=1,NPTS)
          END DO
        ELSE
          OPEN(UNIT=4, FILE=FILNM1, STATUS='OLD', FORM='UNFORMATTED')
          DO I=1,NPTS
            READ(4)  (TSER(K,I), K=1,NST)
          END DO
        END IF
      CASE DEFAULT
        OPEN(UNIT=4, FILE=FILNM1, STATUS='OLD')
        IF (IRD.EQ.1) THEN
          DO K=1,NST
            READ(4,FORMT1)  (TSER(K,I), I=1,NPTS)
          END DO
        ELSE
          DO I=1,NPTS
          DO L=1,NY
            READ(4,FORMT1)  (TSER(K+(L-1)*NX,I), K=1,NX)
          END DO
          END DO
        END IF
      END SELECT

C ------- Open Output Files
      OPEN(UNIT=7, FILE='inform.d', STATUS='UNKNOWN')
      OPEN(UNIT=8, FILE='pcts.d', STATUS='UNKNOWN')
      OPEN(UNIT=9, FILE='avg.d', STATUS='UNKNOWN')
      OPEN(UNIT=10, FILE='emode.d', STATUS='UNKNOWN')

C ------- Smoothing
      DO 10 J=1,NST
        DO I=1,NPTS
          KS = MAX(I-LAG,1)
          KE = MIN(I+LAG,NPTS)
          KN = KE-KS+1
          SUM = 0.0
          DO K=KS,KE
            SUM = SUM + TSER(J,K)
          END DO
          ARR(I) = SUM/FLOAT(KN)
        END DO
        DO I=1,NPTS
          TSER(J,I) = ARR(I)
        END DO
10    CONTINUE

C ------- Remove Mean
      IF (IDTR.NE.0) THEN
        DO J=1,NST
        DO IM=1,IDTR
          SUM = 0.0
          DO I=IM,NPTS,IDTR
            SUM = SUM + TSER(J,I)
          END DO
          SUM = SUM/FLOAT((NPTS-IM)/IDTR+1)
          AVG(IM,J) = SUM
          DO I=IM,NPTS,IDTR
            TSER(J,I) = TSER(J,I) - AVG(IM,J)
          END DO
        END DO
        END DO
        DO IM=1,IDTR
          WRITE(9,'(6E13.5)')  (AVG(IM,J), J=1,NST)
        END DO
      END IF

      NH = (NPTS+1)/2
      DO J=1,NST
C ------- Fourier Analysis
        DO I=1,NPTS
          Y(I) = TSER(J,I)
        END DO
        CALL FOURIER(Y,A,B,NPTS,0,NH)

C ------- Hilbert Transformation
        DO I=1,NPTS
          THLB(I) = A(0)
        END DO
        DO K=1,NH
          FRQ = TPI*FLOAT(K)/FLOAT(NPTS)
          DO I=1,NPTS
            T = FLOAT(I-1)
            THLB(I) = THLB(I) - A(K)*SIN(FRQ*T) + B(K)*COS(FRQ*T)
          END DO
        END DO

C ------- New Complex Variable
        DO I=1,NPTS
          TCMP(J,I) = TSER(J,I) + JIMAG*THLB(I)
        END DO
      END DO

C ------- Covariance Matrix (Hermitian)
      DO 20 I=1,NPTS
      DO 20 J=I,NPTS
        SUMC = (0.,0.)
        DO K=1,NST
          SUMC = SUMC + TCMP(K,I)*CONJG(TCMP(K,J))
        END DO
        DMTRX(I,J) = SUMC/FLOAT(NST)
        DMTRX(J,I) = CONJG(DMTRX(I,J))
20    CONTINUE

C ------- Conversion to a Real Matrix
      DO 30 J=1,NPTS
        JJ = J+NPTS
      DO 30 I=1,NPTS
        II = I+NPTS
        COV(I,J) = REAL(DMTRX(I,J))
        COV(II,JJ) = COV(I,J)
        COV(I,JJ) = -AIMAG(DMTRX(I,J))
        COV(II,J) = -COV(I,JJ)
30    CONTINUE

C ------- Total Variance
      TVAR = 0.0
      DO 35 I=1,NPTS
        TVAR = TVAR + COV(I,I)
35    CONTINUE
      WRITE(7,40) TVAR
40    FORMAT(5X,'TOTAL VARIANCE = ',E15.7,///)

C ------- Call Eigenvalue Routines
      CALL JACOBI(COV,NNPT,NNPT,D,V,NROT)
      CALL EIGSRT(D,V,NNPT,NNPT)
      PRINT *, '# OF JOCOBI ROTATION :', NROT
      PRINT *

C ------- Complex Eigenvectors
      DO 50 I=1,NPTS
        II = 2*I-1
        D(I) = D(II)
      DO 50 J=1,NPTS
        JJ = J + NPTS
        W(J,I) = V(J,II) + JIMAG*V(JJ,II)
50    CONTINUE

C ------- Amplitude and Phase
      IF (MOD(IOUT/2,2).EQ.1) THEN
        DO 60 I=1,NPTS
        DO 60 J=1,NPTS
          AMP(J,I) = CABS(W(J,I))
          PHS(J,I) = ATAN2(AIMAG(W(J,I)), REAL(W(J,I)))*RTD
60      CONTINUE
      END IF

C ------- Write Eigenmodes and Modal Contributions
      SUM = 0.0
      DO 75 I=1,NPTS
        STD = SQRT(D(I)*FLOAT(NST))
        VAR = D(I)/TVAR
        SUM = SUM + VAR
        WRITE(7,65) VAR, SUM
65      FORMAT(5X,'VARIANCE AND CUMULATIVE VARIANCE = ',2E16.7,/)
        IF (MOD(IOUT,2).EQ.1) THEN
          WRITE(8,70) (W(J,I)*STD, J=1,NPTS)
70        FORMAT(6E13.5)
        END IF
        IF (MOD(IOUT/2,2).EQ.1) THEN
          WRITE(8,70) (AMP(J,I)*STD, J=1,NPTS)
          WRITE(8,70) (PHS(J,I), J=1,NPTS)
        END IF
        NMODE = I
        IF (SUM*100. .GE. PVAR)  GO TO 80
75    CONTINUE
80    PRINT *, NMODE

      IF ((MOD(IOUT/4,2).EQ.0) .AND. (MOD(IOUT/8,2).EQ.0))  GO TO 100
      DO 90 J=1,MIN(NMODE,NPRT)
C ------- PC Patterns
        DO K=1,NST
          SUMC = (0.,0.)
          DO I=1,NPTS
            SUMC = SUMC + TCMP(K,I)*CONJG(W(I,J))
          END DO
          PC(K) = SUMC
        END DO

C ------- Normalization
        SUM = 0.0
        DO K=1,NST
          SUM = SUM + PC(K)*CONJG(PC(K))
        END DO
        STD = SQRT(SUM)
        DO K=1,NST
          PC(K) = PC(K)/STD
        END DO

C ------- Outputs
        IF (MOD(IOUT/4,2).EQ.0)  GO TO 85
        DO L=1,NY
          LL = (L-1)*NX
          WRITE(10,70)  (PC(K+LL)*SCL, K=1,NX)
        END DO
85      CONTINUE
        IF (MOD(IOUT/8,2).EQ.0)  GO TO 90
        DO I=1,NST
          AMPL = CABS(PC(I))
          PHSE = ATAN2(AIMAG(PC(I)), REAL(PC(I)))*RTD
          PC(I) = AMPL + JIMAG*PHSE
        END DO
        DO L=1,NY
          LL = (L-1)*NX
          WRITE(10,70)  (REAL(PC(K+LL))*SCL, K=1,NX)
        END DO
        DO L=1,NY
          LL = (L-1)*NX
          WRITE(10,70)  (AIMAG(PC(K+LL)), K=1,NX)
        END DO
90    CONTINUE

100   STOP
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

      PARAMETER (NMX=2000)
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

      SUBROUTINE FOURIER(Y,A,B,NFPTS,NHS,NHE)

      DIMENSION Y(NFPTS), A(0:NFPTS), B(0:NFPTS)


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
