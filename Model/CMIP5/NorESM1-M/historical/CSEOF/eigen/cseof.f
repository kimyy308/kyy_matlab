C ***************************************************************************** 
C 
C        PROGRAM NAME : cseof.f 
C        PROGRAMMER : Dr. Kwang Y. Kim 
C        CODE IDENTIFICATION = CSEOF/VERSION 1.0 
C        CODE CLASSIFICATION = Scientific Computer Code 
C        CREATION DATE = May 27, 1997
C        REVISION DATE = March 28, 2001
C        REVISION INFORMATION = f90 upgrade
C 
C *****************************************************************************

C     This program computes the eigenfunctions of a harmonizable time series
C
C          X(r,t) = Sum_k(0,T-1) [ a_k(r,t) exp(tpi ikt/d) ],
C
C     where
C
C          a_k(r,t) = Int [ w(t-s) X(r,s) exp(-tpi iks/d) ds ]
C
C     and
C
C          w(t) = sin(pi t/d) / (pi t).
C
C     It then follows that the eigenfunctions are 
C
C          f_nm(r,t) = exp(2pi*i*n*t/N) * g_m(r,t),
C
C     where g_m(r,t) is an eigenfunction of the covariance matrix
C
C          C(k,l) = < a_k(r,t) a_l(r,t) >.
C
C     In this version of the code, covariance matrix is computed based on
C     time series at each sampling station.
C
C     REF: Cyclostationarity by Gardner (1994, IEEE Press)
C          EOFs of Harmonizable Cyclostationary Processes
C                by Kim and North (1997, JAS)


      ALLOCATABLE TSER(:,:), TAVG(:,:), TCPY(:), TCOF(:,:,:),
     &            CF(:,:), SF(:,:)
      ALLOCATABLE COV(:,:), D(:), V(:,:), INDX(:)
      ALLOCATABLE W(:,:,:)
      ALLOCATABLE PCTS(:), EOF(:,:)
      ALLOCATABLE WGTS(:)
      CHARACTER*50 FILNM1, FILNM2, FILNM3, FORMT1, FORMT2, FORMT3


      PI = 4.0*ATAN(1.0)
      TPI = 2.0*PI
      DTR = PI/180.0

1     FORMAT(A50)
      PRINT *, '  Do you want to'
      PRINT *, '   0: compute,  1: read Coef TS, or  2: read Cov Mat?'
      READ *, IOPT
      IF (IOPT.EQ.0) THEN
        PRINT *, '  Type name and format of input expansion file.'
        READ 1, FILNM1, FORMT1
      ELSE IF (IOPT.EQ.1) THEN
        PRINT *, '  Type name of coefficient time series file.'
        READ 1, FILNM1
      ELSE
        PRINT *, '  Type name of coefficient time series file.'
        READ 1, FILNM1
        PRINT *, '  Type name of covariance matrix file.'
        READ 1, FILNM2
      END IF
      PRINT *, '  Type the size (NX x NY) of station array.'
      READ *, NX, NY
      NST = NX*NY
      PRINT *, '  Type the total length of time series.'
      READ *, NTOT
      PRINT *, '  First index of the input array: (1: time, 2: space).'
      READ *, IRD
      PRINT *, '  Type the period of the nested cycle.'
      READ *, ICYC
      PRINT *, '  Type the number of spectral points.'
      READ *, NPTS
      PRINT *, '  Type the number of interval subdivisions.'
      READ *, NSUB
      PRINT *, '  Type the cycle period for detrending  (0: No)'
      READ *, IDTR
      PRINT *, '  Type the number of points for covariance calculation.'
      READ *, NCOR
      PRINT *, '  Type the percent variance to be achieved.'
      READ *, PVAR
      PRINT *, '  Type the number of eigenfunctions to be printed.'
      READ *, NPRT
      PRINT *, '  Type the eof scaling factor.'
      READ *, SCL
      PRINT *, '  Type the output option.'
      READ *, IOUT
      IF (2*NPTS .GT. ICYC) THEN
        PRINT *, '  # of spectral points cannot be greater than half',
     &           ' the period of the nested cycle.'
        PRINT *, '  NPTS is now adjusted.'
        NPTS = ICYC/2
      END IF
      MPTS = 2*NPTS+1
      IF (2*NPTS.EQ.ICYC)  MPTS = 2*NPTS
      MATS = NST*MPTS
      NYR = NTOT/ICYC


C ------- Allocate space for each array
      ALLOCATE(TSER(NTOT,NST))
      ALLOCATE(TAVG(IDTR,NST))
      ALLOCATE(TCPY(NTOT))
      ALLOCATE(TCOF(NTOT,ICYC,NST))
      ALLOCATE(CF(NTOT,0:ICYC))
      ALLOCATE(SF(NTOT,0:ICYC))
      ALLOCATE(COV(MATS,MATS))
      ALLOCATE(D(MATS))
      ALLOCATE(V(MATS,MATS))
      ALLOCATE(INDX(MATS))
      ALLOCATE(W(ICYC,NST,MATS))
      ALLOCATE(PCTS(NTOT))
      ALLOCATE(EOF(NST,ICYC))
      ALLOCATE(WGTS(NTOT))


      IF (IOPT.NE.0)  GO TO 20
C ------- Read input data
      SELECT CASE (FORMT1)
      CASE ('DIR')
        IF (IRD.EQ.1) THEN
          OPEN(UNIT=11, FILE=FILNM1, STATUS='OLD',
     &         ACCESS='DIRECT', RECL=NTOT*4)
          DO L=1,NST
            READ(11,REC=L)  (TSER(I,L), I=1,NTOT)
          END DO
        ELSE
          OPEN(UNIT=11, FILE=FILNM1, STATUS='OLD',
     &         ACCESS='DIRECT', RECL=NST*4)
          DO I=1,NTOT
            READ(11,REC=I)  (TSER(I,L), L=1,NST)
          END DO
        END IF
      CASE ('SEQ')
        OPEN(UNIT=11, FILE=FILNM1, STATUS='OLD', FORM='UNFORMATTED')
        IF (IRD.EQ.1) THEN
          DO L=1,NST
            READ(11)  (TSER(I,L), I=1,NTOT)
          END DO
        ELSE
          DO I=1,NTOT
            READ(11)  (TSER(I,L), L=1,NST)
          END DO
        END IF
      CASE DEFAULT
        OPEN(UNIT=11, FILE=FILNM1, STATUS='OLD')
        IF (IRD.EQ.1) THEN
          DO L=1,NST
            READ(11,FORMT1)  (TSER(I,L), I=1,NTOT)
          END DO
        ELSE
          DO I=1,NTOT
          DO L=1,NY
            READ(11,FORMT1)  (TSER(I,K+(L-1)*NX), K=1,NX)
          END DO
          END DO
        END IF
      END SELECT

C ------- Remove mean
      IF (IDTR.NE.0) THEN
        OPEN(UNIT=9, FILE='avg.d', STATUS='UNKNOWN')
        DO L=1,NST
        DO IM=1,IDTR
          SUM = 0.0
          DO I=IM,NTOT,IDTR
            SUM = SUM + TSER(I,L)
          END DO
          SUM = SUM/FLOAT((NTOT-IM)/IDTR+1)
          TAVG(IM,L) = SUM
          DO I=IM,NTOT,IDTR
            TSER(I,L) = TSER(I,L) - TAVG(IM,L)
          END DO
        END DO
          WRITE(9,'(6E13.5)')  (TAVG(IM,L), IM=1,IDTR)
        END DO
      END IF

C ------- Calculate the coefficients
C ----------> modification begins
        NWGT = 2*ICYC
        NSUB = 1
C ----------> modification ends
        CALL WEIGHT(WGTS,ICYC,NPTS,NWGT,NSUB)
      DO L=1,NST
        print *, '  Station =', l
        CALL INTEGR(TSER(1,L),CF,SF,WGTS,NTOT,ICYC,NTOT,NPTS,NWGT,NSUB)
        DO I=1,NTOT
          TCOF(I,1,L) = CF(I,0)
        END DO
        DO K=2,MPTS
          KK = K/2
          IF (MOD(K,2).EQ.0) THEN
            DO I=1,NTOT
              TCOF(I,K,L) = SQRT(2.0)*CF(I,KK)
            END DO
          ELSE
            DO I=1,NTOT
              TCOF(I,K,L) = -SQRT(2.0)*SF(I,KK)
            END DO
          END IF
        END DO
        IF (MPTS.EQ.ICYC) THEN
          DO I=1,NTOT
            TCOF(I,MPTS,L) = CF(I,NPTS)
          END DO
        END IF
      END DO

C ------- Print
      OPEN(UNIT=10, FILE='hcoef.d', STATUS='UNKNOWN',
     &     FORM='UNFORMATTED')
      DO 10 L=1,NST
      DO 10 K=1,MPTS
        WRITE(10)  (TCOF(I,K,L), I=1,NTOT)
10    CONTINUE

C ------- Reconstructed time series
      IF (MOD(IOUT,2).EQ.0)  GO TO 30
      OPEN(UNIT=13, FILE='rc_ts.d', STATUS='UNKNOWN')
      DO L=1,NST
        DO I=1,NTOT
          TCPY(I) = TCOF(I,1,L)
        END DO
        DO 15 K=2,MPTS
          FACT = SQRT(2.0)
          IF (K.EQ.ICYC)  FACT = 1.0
          KK = K/2
          FRQ = TPI*FLOAT(KK)/FLOAT(ICYC)
          IF (MOD(K,2).EQ.0) THEN
            DO I=1,NTOT
              T = FLOAT(I-1)
              TCPY(I) = TCPY(I) + FACT*TCOF(I,K,L)*COS(FRQ*T)
            END DO
          ELSE
            DO I=1,NTOT
              T = FLOAT(I-1)
              TCPY(I) = TCPY(I) + FACT*TCOF(I,K,L)*SIN(FRQ*T)
            END DO
          END IF
15      CONTINUE
        WRITE(13,'(6E13.5)')  (TCPY(I), I=1,NTOT)
      END DO
      GO TO 30

C ------- Read coefficient time series
20    CONTINUE
      IF (IOPT.NE.1)  GO TO 50
      OPEN(UNIT=10, FILE=FILNM1, STATUS='OLD', FORM='UNFORMATTED')
      DO 25 L=1,NST
      DO 25 K=1,MPTS
        READ(10)  (TCOF(I,K,L), I=1,NTOT)
25    CONTINUE

C ------- Covariance matrix
30    CONTINUE
      DO 40 L=1,NST
        LL = (L-1)*MPTS
      DO 40 J=1,MPTS
        JJ = LL + J
      DO 35 M=L,NST
        MM = (M-1)*MPTS
      DO 35 K=1,MPTS
        KK = MM + K

        SUM = 0.0
        DO I=1,NCOR
          SUM = SUM + TCOF(I,J,L)*TCOF(I,K,M)
        END DO
        COV(JJ,KK) = SUM/FLOAT(NCOR)
        COV(KK,JJ) = COV(JJ,KK)
35    CONTINUE
40    CONTINUE

      OPEN(UNIT=31, FILE='covm.d', STATUS='UNKNOWN', FORM='UNFORMATTED')
C      OPEN(UNIT=31, FILE='covm.d', STATUS='UNKNOWN')
      DO J=1,MATS
        WRITE(31)  (COV(I,J), I=1,MATS)
C        WRITE(31,'(6E13.5)')  (COV(I,J), I=1,MATS)
      END DO
      GO TO 60

C ------- Read covariance matrix
50    CONTINUE
      OPEN(UNIT=10, FILE=FILNM1, STATUS='OLD', FORM='UNFORMATTED')
      DO 55 L=1,NST
      DO 55 K=1,MPTS
        READ(10)  (TCOF(I,K,L), I=1,NTOT)
55    CONTINUE
      OPEN(UNIT=31, FILE=FILNM2, STATUS='OLD', FORM='UNFORMATTED')
      DO J=1,MATS
        READ(31)  (COV(I,J), I=1,MATS)
      END DO

60    CONTINUE
C ------- Open output files
      OPEN(UNIT=7, FILE='inform.d', STATUS='UNKNOWN')
      OPEN(UNIT=8, FILE='eigen.d', STATUS='UNKNOWN')
      OPEN(UNIT=32, FILE='Bloch.d', STATUS='UNKNOWN', 
     &     FORM='UNFORMATTED')

C ------- Total variance
      TVAR = 0.0
      DO I=1,MATS
        TVAR = TVAR + COV(I,I)
      END DO
      WRITE(7,65) TVAR
65    FORMAT(5X,'TOTAL VARIANCE = ',E15.7,///)

C ------- Call eigenfunction routines
      CALL JACOBI(COV,MATS,MATS,D,V,NROT)
      CALL EIGSRT(D,V,MATS,MATS)
      PRINT *, '# OF JACOBI ROTATION :', NROT
      PRINT *

C ------- Write eigenmodes and modal contributions
      SUM = 0.0
      DO IM=1,MATS
        VAR = D(IM)/TVAR
        SUM = SUM + VAR
        WRITE(7,70)  VAR, SUM
70      FORMAT(5X,'VARIANCE AND CUMULATIVE VARIANCE = ',2E16.7)
        DO 71 J=1,NST
        DO 71 I=1,MPTS
          IJ = (J-1)*MPTS+I
          W(I,J,IM) = V(IJ,IM)
71      CONTINUE
        DO J=1,NST
          WRITE(8,75)  (W(I,J,IM), I=1,MPTS)
        END DO
75      FORMAT(6E13.5)
        NMODE = IM
        IF (SUM*100. GE. PVAR)  GO TO 80
      END DO
80    PRINT *, NMODE

      IF (MOD(IOUT/2,2).NE.0)
     &    OPEN(UNIT=12, FILE='pcts.d', STATUS='UNKNOWN')
      DO 90 IM=1,MIN(NMODE,NPRT)
C ------- Compute Bloch functions
        DO L=1,NST
          DO I=1,ICYC
            EOF(L,I) = W(1,L,IM)
          END DO
          DO K=2,MPTS
            FACT = SQRT(2.0)
            IF (K.EQ.ICYC)  FACT = 1.0
            KK = K/2
            FRQ = FLOAT(KK)/FLOAT(ICYC)*TPI
            IF (MOD(K,2).EQ.0) THEN
              DO I=1,ICYC
                ANG = FRQ*FLOAT(I-1)
                EOF(L,I) = EOF(L,I) + FACT*W(K,L,IM)*COS(ANG)
              END DO
            ELSE
              DO I=1,ICYC
                ANG = FRQ*FLOAT(I-1)
                EOF(L,I) = EOF(L,I) + FACT*W(K,L,IM)*SIN(ANG)
              END DO
            END IF
          END DO
        END DO
        DO I=1,ICYC
          WRITE(32)  (EOF(L,I), L=1,NST)
        END DO

C ------- PC time series
        IF (MOD(IOUT/2,2).EQ.0)  GO TO 90
        DO I=1,NTOT
          SUM = 0.0
          DO 85 L=1,NST
          DO 85 K=1,MPTS
            SUM = SUM + W(K,L,IM)*TCOF(I,K,L)
85        CONTINUE
          PCTS(I) = SUM
        END DO
        WRITE(12,'(6E13.5)')  (PCTS(I), I=1,NTOT)
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

      PARAMETER (NMAX=5000)
      DIMENSION A(NP,NP), D(NP), V(NP,NP), B(NMAX), Z(NMAX)


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

      SUBROUTINE INTEGR(TSER,CF,SF,WGTS,NTMAX,ICYC,NTOT,NPTS,NWGT,NSUB)

      DIMENSION TSER(NTMAX), CF(NTMAX,0:ICYC), SF(NTMAX,0:ICYC)
      DIMENSION WGTS(0:NWGT*NSUB)

      PI = 4.0*ATAN(1.0)
      TPI = 2.0*PI
      MWGT = NWGT*NSUB

C ------- Modification to include seasonal cycle (May 16, 2000)
      AVG = 0.0
      DO I=1,NTOT
        AVG = AVG + TSER(I)
      END DO
      AVG = AVG/FLOAT(NTOT)

      DO 10 K=0,NPTS
        FRQ = TPI*FLOAT(K)/FLOAT(ICYC)
      DO 10 I=1,NTOT
        T = FLOAT(I-1)
        SUM1 = 0.0
        SUM2 = 0.0
        DO 5 J=-MWGT,MWGT
          TP = T + FLOAT(J)/FLOAT(NSUB)
          IJ = IABS(J)
          J1 = I + J/NSUB
C ------- Turn on the IF STATEMENT below and comment out up to ELSE STATEMENT
C          IF (J1.LT.1 .OR. J1.GT.NTOT)  GO TO 5
          IF (J1.LT.1) THEN
            JJ = MOD(J1+2*ICYC-1,ICYC)+1
            VAL = TSER(JJ)
          ELSE IF (J1.GE.NTOT) THEN
            JJ = NTOT-MOD(NTOT-J1+2*ICYC,ICYC)
            VAL = TSER(JJ)
          ELSE
            J2 = J1 + 1
            W2 = MOD(J+MWGT,NSUB)/FLOAT(NSUB)
            W1 = 1. - W2
            VAL = TSER(J1)*W1 + TSER(J2)*W2
          END IF
          SUM1 = SUM1 + WGTS(IJ)*VAL*COS(FRQ*TP)
          SUM2 = SUM2 - WGTS(IJ)*VAL*SIN(FRQ*TP)
5       CONTINUE
        CF(I,K) = SUM1
        SF(I,K) = SUM2
10    CONTINUE

      RETURN
      END

      SUBROUTINE WEIGHT(WGTS,ICYC,NPTS,NWGT,NSUB)

      DIMENSION WGTS(0:NWGT*NSUB)

      PI = 4.0*ATAN(1.0)
      WGTS(0) = 1./FLOAT(ICYC)
      DO 10 I=1,NWGT*NSUB
        T = FLOAT(I)/FLOAT(NSUB)
        WGT = SIN(PI*T/FLOAT(ICYC))/(PI*T)
        WGTS(I) = WGT
10    CONTINUE

C ------- Normalization
      SWGT = FLOAT(NSUB)
      DO I=0,NWGT*NSUB
        WGTS(I) = WGTS(I)/SWGT
      END DO

      RETURN
      END
