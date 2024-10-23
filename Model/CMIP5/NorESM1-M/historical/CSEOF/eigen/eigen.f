C ***************************************************************************** C
C        PROGRAM NAME : eigen.f
C        PROGRAMMER : Dr. Kwang Y. Kim
C        CODE IDENTIFICATION NUMBER = EIGEN/VERSION 1.0
C        CODE CLASSIFICATION = Scientific Computer Code
C        CREATION DATE = August 1, 1997
C        REVISION DATE = not revised
C        REVISION INFORMATION = not applicable
C
C *****************************************************************************


C     This program computes EOFs of a given data:
C
C          C x = lambda x
C
C     where C is an (N x N) real symmetric matrix, x is an eigenvector, and
C     lambda is the corresponding eigenvalue.  The eigen equation is solved
C     via Jacobi rotation method.


      ALLOCATABLE TSER(:,:), ARR(:), AVG(:,:)
      ALLOCATABLE COV(:,:), EIGV(:), EIGF(:,:), EIGC(:), INDX(:)
      ALLOCATABLE PCTS(:)
      CHARACTER*50 FILNM1, FILNM2, FILNM3, FORMT1, FORMT2, FORMT3

      PI = 4.0*ATAN(1.0)
      TPI = 2.0*PI


5     FORMAT(A50)
      PRINT *, '  Type the name and format of the time series data.'
      READ 5, FILNM1, FORMT1
      PRINT *, '  Type the dimension (NX,NY) of sampling stations.'
      READ *, NX, NY
      NST = NX*NY
      FNST = FLOAT(NST)
      SNST = SQRT(FNST)
      PRINT *, '  Type the number of samples at each station.'
      READ *, NPTS
      PRINT *, '  First index of the input array: (1: time, 2: space).'
      READ *, IRD
      PRINT *, '  Type smoothing option (0:No, M:Moving average lag).'
      READ *, LAG
      PRINT *, '  Type the cycle period for detrending (0: No).'
      READ *, IDTR
      PRINT *, '  Type the area adjustment option (0: No, 1: Yes).'
      READ *, IARA
      PRINT *, '  Type the starting latitude and increment.'
      READ *,  SLAT, DLAT
      PRINT *, '  Type percent variance to be achieved.'
      READ *, PVAR
      PRINT *, '  Type the eigenfunction scale factor.'
      READ *, SCL
      PRINT *, '  Type the number of eigenfunctions to be printed.'
      READ *, NPRT
      PRINT *, '  Type the PC normalization option (0: no, 1: yes)'
      READ *, INORM
      PRINT *, '  Type the name and format of the EOF output file.'
      READ 5, FILNM2, FORMT2
      PRINT *, '  Type the name and format of the PC time series file.'
      READ 5, FILNM3, FORMT3

C ------- Allocate space for dynamic arrays
      ALLOCATE(TSER(NPTS,NST))
      ALLOCATE(ARR(NPTS))
      ALLOCATE(AVG(NST,IDTR))
      ALLOCATE(COV(NST,NST))
      ALLOCATE(EIGV(NST))
      ALLOCATE(EIGF(NST,NST))
      ALLOCATE(EIGC(NST))
      ALLOCATE(INDX(NST))
      ALLOCATE(PCTS(NPTS))


C ------- Read Input Data
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

C ------- Subtract Mean
      IF (IDTR.NE.0) THEN
        OPEN(UNIT=9, FILE='avg.d', STATUS='UNKNOWN')
        DO 15 J=1,NST
        DO 15 IM=1,IDTR
          SUM = 0.0
          DO K=IM,NPTS,IDTR
            SUM = SUM + TSER(K,J)
          END DO
          AVG(J,IM) = SUM/FLOAT((NPTS-IM)/IDTR+1)
          DO K=IM,NPTS,IDTR
            TSER(K,J) = TSER(K,J) - AVG(J,IM)
          END DO
15      CONTINUE
        DO IM=1,IDTR
          WRITE(9,'(5E15.7)')  (AVG(J,IM), J=1,NST)
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
          SCLF = SQRT(COS(ALAT*PI/180.)/TAREA)*SNST
          DO I=1,NX
            IJ = (J-1)*NX+I
            TSER(K,IJ) = TSER(K,IJ)*SCLF
          END DO
        END DO
        END DO
      END IF

C ------- Covariance Matrix
      DO 20 J=1,NST
      DO 20 I=J,NST
        SUMC = 0.0
        DO K=1,NPTS
          SUMC = SUMC + TSER(K,I)*TSER(K,J)
        END DO
        SUMC = SUMC/FLOAT(NPTS)
        COV(I,J) = SUMC
        COV(J,I) = SUMC
20    CONTINUE

C ------- Total Variance
      OPEN(UNIT=7, FILE='inform.d', STATUS='UNKNOWN')
      TVAR = 0.0
      DO 25 I=1,NST
        TVAR = TVAR + COV(I,I)
25    CONTINUE
      TVAR = TVAR/FNST
      WRITE(7,30) TVAR
30    FORMAT(5X,'TOTAL VARIANCE = ',E15.7,///)

C ------- Call Eigenvalue Routines
      CALL JACOBI(COV,NST,NST,EIGV,EIGF,NROT)
      CALL EIGSRT(EIGV,EIGF,NST,NST)
      PRINT *, '# OF JOCOBI ROTATION :', NROT
      PRINT *

C ------- Write Eigenmodes and Modal Contributions
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

      SUM = 0.0
      DO 45 I=1,NST
        EIGV(I) = EIGV(I)/FNST
        VAR = EIGV(I)/TVAR
        SUM = SUM + VAR
        WRITE(7,35) VAR, SUM
35      FORMAT(5X,'VARIANCE AND CUMULATIVE VARIANCE = ',2E16.7,/)
        DO L=1,NY
          LL = (L-1)*NX
          IF (IARA.EQ.1) THEN
            ALAT = SLAT + FLOAT(L-1)*DLAT
            SCLF = SQRT(COS(ALAT*PI/180.)/TAREA)*SNST
          ELSE
            SCLF = 1.0
          END IF
          DO K=1,NX
            IF (INORM.EQ.1)  THEN
              EIGC(K+LL) = EIGF(K+LL,I)*SNST/SCLF*SQRT(EIGV(I))*SCL
            ELSE
              EIGC(K+LL) = EIGF(K+LL,I)*SNST/SCLF*SCL
            END IF
          END DO
        END DO
        NMODE = I

        SELECT CASE (FORMT2)
        CASE ('DIR')
          WRITE(21,REC=I)  (EIGC(K), K=1,NST)
        CASE ('SEQ')
          WRITE(21)  (EIGC(K), K=1,NST)
        CASE DEFAULT
          DO L=1,NY
            LL = (L-1)*NX
            WRITE(21,FORMT2) (EIGC(K+LL), K=1,NX)
          END DO 
        END SELECT

        IF (SUM*100. .GE. PVAR)  GO TO 50
45    CONTINUE
50    PRINT *, NMODE

C ------- PC time series
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

      DO 60 J=1,MIN(NMODE,NPRT)
        DO I=1,NPTS
          SUM = 0.0
          DO K=1,NST
            SUM = SUM + TSER(I,K)*EIGF(K,J)
          END DO
          SUM = SUM/SNST
          IF (INORM.EQ.1) THEN
            PCTS(I) = SUM/SQRT(EIGV(J))
          ELSE
            PCTS(I) = SUM
          END IF
        END DO
        SELECT CASE (FORMT3) 
        CASE ('DIR') 
          WRITE(22,REC=J)  (PCTS(I), I=1,NPTS)
        CASE ('SEQ') 
          WRITE(22)  (PCTS(I), I=1,NPTS)
        CASE DEFAULT 
          WRITE(22,FORMT3)  (PCTS(I), I=1,NPTS)
        END SELECT
60    CONTINUE

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

      PARAMETER (NMX=6000)
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
