C *****************************************************************************
C
C        PROGRAM NAME : exteof.f
C        PROGRAMMER : Dr. Kwang Y. Kim
C        CODE IDENTIFICATION NUMBER = EXTEOF/VERSION 1.0
C        CODE CLASSIFICATION = Scientific Computer Code
C        CREATION DATE = August 4, 1997
C        REVISION DATE = November 5, 2002
C        REVISION INFORMATION = f90 upgrade (November 5, 2002)
C
C *****************************************************************************


C     This program computes extended EOFs of a given data
C
C          X(i,j), i=1,NX, j=1,NT
C
C     Then, the covariance matrix is given by
C
C          C = < X(i,t-k) X(j,t-l) >,  i,j=1,NX  and  k,l=0,nl
C
C     where nl being maximum lag.  Thus, the rank of the matrix is NX*(nl+1).
C     The resulting eigen problem is
C
C          C U = lambda U
C
C     where U is the gienvector.  The PC time series, then, are given by
C
C          Z(i,j) = U(i,0)*X(i,j) + U(i,1)*X(i,j-1) + ... + U(i,nl)*X(i,j-nl)


      ALLOCATABLE TSER(:,:), AVG(:,:), ARR(:)
      ALLOCATABLE DMTRX(:,:), PCTS(:)
      ALLOCATABLE COV(:,:), D(:), V(:,:), INDX(:)
      CHARACTER*50 FILNMI, FORMTI, FILNMO, FORMTO
      CHARACTER*3 ANS, SEQ(10)
      DATA SEQ / '0TH', '1ST', '2ND', '3RD', '4TH',
     &           '5TH', '6TH', '7TH', '8TH', '9TH' /


      PI = 4.0*ATAN(1.0)
      TPI = 2.0*PI


5     FORMAT(A50)
      PRINT *, '  Type name and format of the input file.'
      READ 5, FILNMI, FORMTI
      PRINT *, '  First index of the input array (1: time,  2: space).'
      READ *, IRD
      PRINT *, '  Type the dimension (NX,NY) of sampling stations.'
      READ *, NX, NY
      NST = NX*NY
      PRINT *, '  Type the maximum lag and lag interval of analysis.'
      READ *, NL, LD
      NNST = NST*(NL/LD+1)
      PRINT *, '  Type the length of time series at each station.'
      READ *, NTOT
      NPTS = NTOT-NL
      PRINT *, '  Type the cycle period for detrending  (0: No).'
      READ *, IDTR
      PRINT *, '  Type smoothing option (0:No, M:Moving average lag).'
      READ *, LAG
      PRINT *, '  Type percent variance to be achieved.'
      READ *, PVAR
      PRINT *, '  Type extended EOF scale factor.'
      READ *, SCL


C ------- Allocate Dynamic Memory
      ALLOCATE(TSER(NTOT,NST))
      ALLOCATE(AVG(NST,IDTR))
      ALLOCATE(ARR(NTOT))
      ALLOCATE(DMTRX(NST,NST))
      ALLOCATE(PCTS(NPTS))
      ALLOCATE(COV(NNST,NNST))
      ALLOCATE(D(NNST))
      ALLOCATE(V(NNST,NNST))
      ALLOCATE(INDX(NNST))

C ------- Read Input File
      SELECT CASE (FORMTI)
      CASE('DIR')
        IF (IRD.EQ.1) THEN
          OPEN(UNIT=4, FILE=FILNMI, STATUS='OLD',
     &         ACCESS='DIRECT', RECL=NTOT)
          DO K=1,NST
            READ(4,REC=K)  (TSER(I,K), I=1,NTOT)
          END DO
        ELSE
          OPEN(UNIT=4, FILE=FILNMI, STATUS='OLD',
     &         ACCESS='DIRECT', RECL=NST)
          DO I=1,NTOT
            READ(4,REC=I)  (TSER(I,K), K=1,NST)
          END DO
        END IF
      CASE ('SEQ')
        IF (IRD.EQ.1) THEN
          OPEN(UNIT=4, FILE=FILNMI, STATUS='OLD', FORM='UNFORMATTED')
          DO K=1,NST
            READ(4)  (TSER(I,K), I=1,NTOT)
          END DO
        ELSE
          OPEN(UNIT=4, FILE=FILNMI, STATUS='OLD', FORM='UNFORMATTED')
          DO I=1,NTOT
            READ(4)  (TSER(I,K), K=1,NST)
          END DO
        END IF
      CASE DEFAULT
        OPEN(UNIT=4, FILE=FILNMI, STATUS='OLD')
        IF (IRD.EQ.1) THEN
          DO K=1,NST
            READ(4,FORMTI)  (TSER(I,K), I=1,NTOT)
          END DO
        ELSE
          DO I=1,NTOT
          DO L=1,NY
            READ(4,FORMTI)  (TSER(I,K+(L-1)*NX), K=1,NX)
          END DO
          END DO
        END IF
      END SELECT
      CLOSE(UNIT=4)

C ------- Open Output Files
      OPEN(UNIT=7, FILE='inform.d', STATUS='UNKNOWN')
      OPEN(UNIT=8, FILE='emode.d', STATUS='UNKNOWN')
      OPEN(UNIT=9, FILE='avg.d', STATUS='UNKNOWN')
      OPEN(UNIT=10, FILE='pcts.d', STATUS='UNKNOWN')

C ------- Remove Mean or Seasonal Cycle
      IF (IDTR.NE.0) THEN
        DO J=1,NST
          DO IM=1,IDTR
            SUM = 0.0
            DO I=IM,NTOT,IDTR
              SUM = SUM + TSER(I,J)
            END DO
            AVG(J,IM) = SUM/FLOAT((NTOT-IM)/IDTR+1)
            DO I=IM,NTOT,IDTR
              TSER(I,J) = TSER(I,J) - AVG(J,IM)
            END DO
          END DO
        END DO
        DO IM=1,IDTR
          WRITE(9,'(6E13.5)')  (AVG(J,IM), J=1,NST)
        END DO
      END IF

C ------- Smoothing
      DO J=1,NST
        DO I=1,NTOT
          KS = MAX(I-LAG,1)
          KE = MIN(I+LAG,NTOT)
          KN = KE-KS+1
          SUM = 0.0
          DO K=KS,KE
            SUM = SUM + TSER(K,J)
          END DO
          ARR(I) = SUM/FLOAT(KN)
        END DO
        DO I=1,NTOT
          TSER(I,J) = ARR(I)
        END DO
      END DO

C ------- Covariance Matrix (Lag Extended)
      DO 30 JLG=0,NL,LD
      DO 30 ILG=JLG,NL,LD
        LAG = IABS(ILG-JLG)
        DO 10 J=1,NST
        DO 10 I=1,NST
          SUM = 0.0
          DO K=1,NPTS
            SUM = SUM + TSER(K,I)*TSER(K+LAG,J)
          END DO
          DMTRX(I,J) = SUM/FLOAT(NPTS)
10      CONTINUE

        DO 20 J=1,NST
          JJ = J + (JLG/LD)*NST
        DO 20 I=1,NST
          II = I + (ILG/LD)*NST
          COV(II,JJ) = DMTRX(I,J)
          COV(JJ,II) = DMTRX(I,J)
20      CONTINUE
30    CONTINUE

C ------- Total Variance
      TVAR = 0.0
      DO 35 I=1,NNST
        TVAR = TVAR + COV(I,I)
35    CONTINUE
      WRITE(7,40) TVAR
40    FORMAT(5X,'TOTAL VARIANCE = ',E15.7,///)

C ------- Call Eigenvalue Routines
      CALL JACOBI(COV,NNST,NNST,D,V,NROT)
      CALL EIGSRT(D,V,NNST,NNST)
      PRINT *, '# OF JOCOBI ROTATION :', NROT
      PRINT *

C ------- Write Eigenmodes and Modal Contributions
      SUM = 0.0
      DO I=1,NNST
        VAR = D(I)/TVAR
        SUM = SUM + VAR
        WRITE(7,45) VAR, SUM
45      FORMAT(5X,'VARIANCE AND CUMULATIVE VARIANCE = ',2E16.7,/)
        DO 55 LAG=0,NL,LD
        DO 55 L=1,NY
          LL = (L-1)*NX + LAG/LD*NST
          WRITE(8,50) (SCL*V(K+LL,I), K=1,NX)
50        FORMAT(6E13.5)
55      CONTINUE
        NMODE = I
        IF (SUM*100. .GE. PVAR)  GO TO 60
      END DO
60    PRINT *, NMODE

C ------- PC Time Series
      DO 70 J=1,MIN(NST,10)
      DO I=1,NPTS
        SUM = 0.0
        DO 65 LAG=0,NL,LD
        DO 65 K=1,NST
          KK = K + LAG/LD*NST
          SUM = SUM + TSER(I+LAG,K)*V(KK,J)
65      CONTINUE
        PCTS(I) = SUM
      END DO
        WRITE(10,'(6E13.5)')  (PCTS(I), I=1,NPTS)
70    CONTINUE

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

      PARAMETER (NMX=4000)
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
