C ***************************************************************************** 
C 
C        PROGRAM NAME : cseof_h.f 
C        PROGRAMMER : Dr. Kwang Y. Kim 
C        CODE IDENTIFICATION = CSEOF_H/VERSION 1.0 
C        CODE CLASSIFICATION = Scientific Computer Code 
C        CREATION DATE = April 08, 1996
C        REVISION DATE = not revised 
C        REVISION INFORMATION = not applicable 
C 
C *****************************************************************************

C     This program computes the eigenfunctions of a harmonizable time series
C
C          X(t) = Sum_k(0,T-1) [ a_k(t) exp(tpi ikt/d) ],
C
C     where
C
C          a_k(t) = Int [ w(t-s) X(s) exp(-tpi iks/d) ds ]
C
C     and
C
C          w(t) = sin(pi t/d) / (pi t).
C
C     It then follows that the eigenfunctions are 
C
C          f_nm(t) = exp(2pi*i*n*t/N) * g_m(t),
C
C     where g_m(t) is an eigenfunction of the covariance matrix
C
C          C(k,l) = < a_k(t) a_l(t) >.
C
C     REF: Cyclostationarity by Gardner (1994, IEEE Press)


      ALLOCATABLE TSER(:), TCPY(:), TCOF(:,:), CF(:,:), SF(:,:),
     &            EIGB(:,:), EIGR(:), EIGI(:)
      ALLOCATABLE COV(:,:), D(:), V(:,:), INDX(:)
      ALLOCATABLE EGV(:), NMD(:,:), PCT(:,:), COVR(:,:), CVAR(:,:)
      ALLOCATABLE EGVS(:), NMDS(:,:), IORD(:), JORD(:)
      ALLOCATABLE WNDOW(:), CN(:), WGTS(:)
      ALLOCATABLE Y(:), AN(:), BN(:)
      COMPLEX, DIMENSION(:,:), ALLOCATABLE :: COF
      COMPLEX JIMAG
      CHARACTER*50 FILENM, FORMAT
      DATA JIMAG / (0.,1.) /


      PI = 4.0*ATAN(1.0)
      TPI = 2.0*PI
      DTR = PI/180.0

5     FORMAT(A50)
      PRINT *, '  Type name and format of the input file.'
      READ 5, FILENM, FORMAT
      PRINT *, '  Type the percent variance to be achieved.'
      READ *, PVAR
      PRINT *, '  Type the total number of time series.'
      READ *, NTOT
      PRINT *, '  Type the period of the nested cycle.'
      READ *, ICYC
      PRINT *, '  Type the number of spectral points.'
      READ *, NPTS
      PRINT *, '  Type the number of interval subdivisions.'
      READ *, NSUB
      PRINT *, '  Type the number of eigenfunctions to be printed.'
      READ *, NPRT
      PRINT *, '  Type the length of the output.'
      READ *, LPRT
      PRINT *, '  Type the resolution of the output.'
      READ *, NRES
      PRINT *, '  Type the cycle period for detrending  (0: No).'
      READ *, IDTR
      PRINT *, '  Type the type of smoothing window.'
      READ *, IWNDW
      PRINT *, '  Type the lag of the smoothing window.'
      READ *, NLAG
      PRINT *, '  Type the output option.'
      READ *, IOUT
      IF (2*NPTS .GT. ICYC) THEN
        PRINT *, '  Number of spectral points cannot be greater than',
     &           '     half the period of the nested cycle.'
        PRINT *, '  NPTS is now adjusted.'
        NPTS = ICYC/2
      END IF
      MTOT = NTOT/ICYC

C ------- Allocate space for dynamic variables
      ALLOCATE(TSER(NTOT))
      ALLOCATE(TCPY(NTOT))
      ALLOCATE(TCOF(NTOT,ICYC))
      ALLOCATE(CF(NTOT,0:ICYC))
      ALLOCATE(SF(NTOT,0:ICYC))
      ALLOCATE(EIGB(NTOT,ICYC))
      ALLOCATE(EIGR(NTOT))
      ALLOCATE(EIGI(NTOT))
      ALLOCATE(COV(ICYC,ICYC))
      ALLOCATE(D(ICYC))
      ALLOCATE(V(ICYC,ICYC))
      ALLOCATE(INDX(ICYC))
      ALLOCATE(EGV(NTOT))
      ALLOCATE(NMD(2,NTOT))
      ALLOCATE(PCT(NTOT,ICYC))
      ALLOCATE(COVR(ICYC,0:NTOT))
      ALLOCATE(CVAR(MTOT,MTOT))
      ALLOCATE(EGVS(NTOT))
      ALLOCATE(NMDS(2,NTOT))
      ALLOCATE(IORD(NTOT))
      ALLOCATE(JORD(NTOT))
      ALLOCATE(WNDOW(0:NTOT))
      ALLOCATE(CN(0:NTOT))
      ALLOCATE(WGTS(NTOT))
      ALLOCATE(COF(0:ICYC,ICYC))
      ALLOCATE(Y(NTOT))
      ALLOCATE(AN(0:NTOT))
      ALLOCATE(BN(0:NTOT))


C ------- Read data
      SELECT CASE (FORMAT)
      CASE ('DIR')
        OPEN(UNIT=11, FILE=FILENM, STATUS='OLD',
     &       ACCESS='DIRECT', RECL=NTOT*4)
        READ(11,REC=1)  (TSER(I), I=1,NTOT)
      CASE ('SEQ')
        OPEN(UNIT=11, FILE=FILENM, STATUS='OLD', FORM='UNFORMATTED')
        READ(11)  (TSER(I), I=1,NTOT)
      CASE DEFAULT
        OPEN(UNIT=11, FILE=FILENM, STATUS='OLD')
        READ(11,FORMAT)  (TSER(I), I=1,NTOT)
      END SELECT

C ------- Remove mean
      IF (IDTR.NE.0) THEN
        DO IM=1,IDTR
          SUM = 0.0
          DO I=IM,NTOT,IDTR
            SUM = SUM + TSER(I)
          END DO
          SUM = SUM/FLOAT((NTOT-IM)/IDTR+1)
          DO I=IM,NTOT,IDTR
            TSER(I) = TSER(I) - SUM
          END DO
          write(17,'(6e13.5)')  sum
        END DO
      END IF

C ------- Calculate the coefficients
      NWGT = 2*ICYC
      NSUB = 1
      CALL WEIGHT(WGTS,ICYC,NPTS,NWGT,NSUB)
      CALL INTEGR(TSER,CF,SF,WGTS,NTOT,ICYC,NTOT,NPTS,NWGT,NSUB)

C ------- Coefficient time series
      MPTS = 2*NPTS+1
      IF (2*NPTS.EQ.ICYC)  MPTS = 2*NPTS
      DO I=1,NTOT
        TCOF(I,1) = CF(I,0)
      END DO
      DO K=2,MPTS
        KK = K/2
        IF (MOD(K,2).EQ.0) THEN
          DO I=1,NTOT
            TCOF(I,K) = SQRT(2.0)*CF(I,KK)
          END DO
        ELSE
          DO I=1,NTOT
            TCOF(I,K) = -SQRT(2.0)*SF(I,KK)
          END DO
        END IF
      END DO
      IF (MPTS.EQ.ICYC) THEN
        DO I=1,NTOT
          TCOF(I,MPTS) = CF(I,NPTS)
        END DO
      END IF

C ------- Test
      IF (MOD(IOUT,2).EQ.0)  GO TO 20
      DO I=1,NTOT
        TCPY(I) = TCOF(I,1)
      END DO
      DO 10 K=2,MPTS
        FACT = SQRT(2.0)
        IF (K.EQ.ICYC)  FACT = 1.0
        KK = K/2
        FRQ = TPI*FLOAT(KK)/FLOAT(ICYC)
        IF (MOD(K,2).EQ.0) THEN
          DO I=1,NTOT
            T = FLOAT(I-1)
            TCPY(I) = TCPY(I) + FACT*TCOF(I,K)*COS(FRQ*T)
          END DO
        ELSE 
          DO I=1,NTOT
            T = FLOAT(I-1)
            TCPY(I) = TCPY(I) + FACT*TCOF(I,K)*SIN(FRQ*T)
          END DO
        END IF
10    CONTINUE
      OPEN(UNIT=7, FILE='rec_ts.d', STATUS='UNKNOWN')
      WRITE(7,'(6E13.5)')  (TCPY(I), I=1,NTOT)

C ------- Print
20    CONTINUE
      OPEN(UNIT=8, FILE='hcoef.d', STATUS='UNKNOWN')
      DO K=1,MPTS
        WRITE(8,'(6E13.5)')  (TCOF(I,K), I=1,NTOT)
      END DO

C ------- Covariance function
      DO 30 K=1,MPTS
      DO 30 L=1,MPTS
        SUM = 0.0
        DO I=1,NTOT
          SUM = SUM + TCOF(I,K)*TCOF(I,L)
        END DO
        COV(K,L) = SUM/FLOAT(NTOT)
30    CONTINUE

C ------- Open output files
      OPEN(UNIT=9, FILE='inform.d', STATUS='UNKNOWN')
      OPEN(UNIT=10, FILE='eigen.d', STATUS='UNKNOWN')

C ------- Total variance
      TVAR = 0.0
      DO I=1,MPTS
        TVAR = TVAR + COV(I,I)
      END DO
      WRITE(9,35) TVAR
35    FORMAT(5X,'TOTAL VARIANCE = ',E15.7,///)

C ------- Call eigenfunction routines
      CALL JACOBI(COV,MPTS,ICYC,D,V,NROT)
      CALL EIGSRT(D,V,MPTS,ICYC)
      PRINT *, '# OF JACOBI ROTATION :', NROT
      PRINT *

C ------- Write eigenmodes and modal contributions
      SUM = 0.0
      DO 50 I=1,MPTS
        VAR = D(I)/TVAR
        SUM = SUM + VAR
        WRITE(9,40)  VAR, SUM
40      FORMAT(5X,'VARIANCE AND CUMULATIVE VARIANCE = ',2E16.7)
        WRITE(10,45)  (V(J,I), J=1,MPTS)
45      FORMAT(5X,4E15.7)
        NMODE = I
        IF (SUM*100. GE. PVAR)  GO TO 60
50    CONTINUE
60    PRINT *, NMODE

C ------- Write the Bloch function
      IF (MOD(IOUT/2,2).EQ.0 .AND. IOUT/8.EQ.0)  GO TO 70
      OPEN(UNIT=11, FILE='bloch.d', STATUS='UNKNOWN')
      DO IM=1,NMODE
        DO I=1,NRES*ICYC+1
          EIGB(I,IM) = V(1,IM)
        END DO
        DO K=2,MPTS
          FACT = SQRT(2.0)
          IF (K.EQ.ICYC)  FACT = 1.0
          KK = K/2
          FRQ = FLOAT(KK)/FLOAT(ICYC)*TPI
          IF (MOD(K,2).EQ.0) THEN
            DO I=1,NRES*ICYC+1
              ANG = FRQ*FLOAT(I-1)/FLOAT(NRES)
              EIGB(I,IM) = EIGB(I,IM) + FACT*V(K,IM)*COS(ANG)
            END DO
          ELSE
            DO I=1,NRES*ICYC+1
              ANG = FRQ*FLOAT(I-1)/FLOAT(NRES)
              EIGB(I,IM) = EIGB(I,IM) + FACT*V(K,IM)*SIN(ANG)
            END DO
          END IF
        END DO
        IF (MOD(IOUT/2,2).EQ.1)
C     &      WRITE(11,'(6E13.5)')  (EIGB(I,IM), I=1,NRES*ICYC+1)
     &      WRITE(11,'(6E13.5)')  (EIGB(I,IM), I=1,NRES*ICYC)
      END DO

C ------- Write the PC time series
70    CONTINUE
      IF (MOD(IOUT/4,2).EQ.0 .AND. IOUT/8.EQ.0)  GO TO 80
      OPEN(UNIT=12, FILE='pc_ts.d', STATUS='UNKNOWN')
      DO IM=1,NMODE
        DO I=1,NTOT
          PCT(I,IM) = 0.0
        END DO
        DO 75 K=1,MPTS
        DO 75 I=1,NTOT
          PCT(I,IM) = PCT(I,IM) + V(K,IM)*TCOF(I,K)
75      CONTINUE
        IF (MOD(IOUT/4,2).EQ.1)
     &      WRITE(12,'(6E13.5)')  (PCT(I,IM), I=1,NTOT)
      END DO

C ------- Covariance functions
80    CONTINUE
      IF (MOD(IOUT/8,2).EQ.0)  GO TO 90
      OPEN(UNIT=13, FILE='cov.d', STATUS='UNKNOWN')
      DO 85 IM=1,NMODE
      DO 85 K=1,ICYC
        K1 = (K-1)*NRES + 1
      DO 85 LAG=0,NTOT-1
        K2 = MOD(K+LAG-1,ICYC)*NRES + 1
        COVR(K,LAG) = 0.0
        SUM = 0.0
        DO L=K,NTOT-LAG,ICYC
          SUM = SUM + PCT(L,IM)*PCT(L+LAG,IM)
        END DO
        COVR(K,LAG) = COVR(K,LAG)
     &              + SUM/FLOAT(NTOT/ICYC)*EIGB(K1,IM)*EIGB(K2,IM)
85    CONTINUE

      DO I=1,MIN(30,NTOT)
      DO 86 LAG=0,MIN(NTOT-1,30)
        II = MOD(I-1,ICYC) + 1
        J = I+LAG
        IF (J.GT.30)  GO TO 86
        CVAR(I,J) = COVR(II,LAG)
        CVAR(J,I) = CVAR(I,J)
86    CONTINUE
      END DO
      DO J=1,30
C       WRITE(13,'(I5)')  J
        WRITE(13,'(6E13.5)')  (CVAR(I,J), I=1,30)
      END DO

C ------- Eigenvalues
90    CONTINUE
      IF (MOD(IOUT/16,2).EQ.0 .AND. IOUT/32.EQ.0)  GO TO 100
      NHARM = NTOT/ICYC
      MHARM = NHARM/2
      IF (ICYC.EQ.1) THEN
        NHARM = NTOT/2
        MHARM = NHARM
      END IF
      MMODE = (MHARM+1)*NMODE
      DF = 0.5/MHARM
C ------- Spectral window generator
      AMP = TPI*FLOAT(NLAG)
      IF (IWNDW.EQ.0) THEN
C ---------- Truncated Periodogram Window
        WNDOW(0) = AMP/PI
        DO J=1,MHARM
          FRQ = AMP*FLOAT(J)*DF
          WNDOW(J) = AMP*SIN(FRQ)/(PI*FRQ)
        END DO
      ELSE IF (IWNDW.EQ.1) THEN
C ---------- Bartlett Window
        WNDOW(0) = AMP/TPI
        DO J=1,MHARM
          FRQ = AMP*FLOAT(J)*DF
          WNDOW(J) = AMP*(SIN(FRQ/2.)/(FRQ/2.))**2/TPI
        END DO
      ELSE IF (IWNDW.EQ.2) THEN
C ---------- Parzen Window
        WNDOW(0) = AMP*3./(8.*PI)
        DO J=1,MHARM
          FRQ = AMP*FLOAT(J)*DF
          WNDOW(J) = AMP*3.*(SIN(FRQ/4.)/(FRQ/4.))**4/(8.*PI)
        END DO
      ELSE IF (IWNDW.EQ.3) THEN
C ---------- Hanning Window
        A = 0.5
        B = 0.5
        WNDOW(0) = AMP*A/PI
        DO J=1,MHARM
          FRQ = AMP*FLOAT(J)*DF
          WNDOW(J) = AMP*SIN(FRQ)/(PI*FRQ)*(A + B*FRQ**2/(PI**2-FRQ**2))
          IF (ABS(PI-FRQ) .LT. 1.E-5) WNDOW(J) = AMP*B/(TPI)
        END DO
      ELSE IF (IWNDW.EQ.4) THEN
C ---------- Tukey Window
        A = 0.54
        B = 0.46
        WNDOW(0) = AMP*A/PI
        DO J=1,MHARM
          FRQ = AMP*FLOAT(J)*DF
          WNDOW(J) = AMP*SIN(FRQ)/(PI*FRQ)*(A + B*FRQ**2/(PI**2-FRQ**2))
          IF (ABS(PI-FRQ) .LT. 1.E-5) WNDOW(J) = AMP*B/(TPI)
        END DO
      END IF

      OPEN(UNIT=14, FILE='eigv.d', STATUS='UNKNOWN')
      DO IM=1,NMODE
        DO I=1,NTOT
          Y(I) = PCT(I,IM)
        END DO
        CALL FOURIER(NTOT,0,MHARM,Y,AN,BN)
        DO I=0,MHARM
          II = 2*MHARM-I
          CN(I) = 0.5*(AN(I)**2 + BN(I)**2)
          CN(II) = CN(I)
        END DO
        DO I=0,MHARM
          II = I*NMODE + IM
          SUM = 0.0
          DO 91 J=0,2*MHARM
            JJ = IABS(I-J)
            IF (JJ.GT.MHARM)  GO TO 91
            SUM = SUM + CN(J)*WNDOW(JJ)
91        CONTINUE
          EGV(II) = SUM*DF
          NMD(1,II) = I
          NMD(2,II) = IM
        END DO
      END DO

      DO I=1,NTOT
        Y(I) = TSER(I)
      END DO
      CALL FOURIER(NTOT,0,NTOT,Y,AN,BN)
      SUM1 = 0.0
      SUM2 = 0.0
      DO IM=1,NMODE
        NM = (IM-1)*NHARM + MHARM
        II = MHARM*NMODE + IM
        SUM1 = SUM1 + EGV(II)
        SUM2 = SUM2 + 0.5*(AN(NM)**2 + BN(NM)**2)
      END DO
      SCL = SUM1/SUM2
      DO IM=1,NMODE
        NM = (IM-1)*NHARM + MHARM
        II = MHARM*NMODE + IM
        EGV(II) = 0.5*(AN(NM)**2 + BN(NM)**2)*SCL
        NMD(1,II) = MHARM
        NMD(2,II) = IM
      END DO

C ------- Ordering and rank of values
      CALL SORT(MMODE,EGV,EGVS,IORD,JORD)
      CUMV = 0.0
      DO I=MMODE,1,-1
        NMDS(1,I) = NMD(1,IORD(I))
        NMDS(2,I) = NMD(2,IORD(I))
        EGVS(I) = EGVS(I)/TVAR
        CUMV = CUMV + EGVS(I)
        IF (MOD(IOUT/16,2).EQ.1)
     &      WRITE(14,95)  NMDS(1,I), NMDS(2,I), EGVS(I), CUMV
95      FORMAT('VAR AND CUM VAR OF MODE (',I4,',',I2,') = ',2E15.7)
      END DO

C ------- Eigenfunctions
      IF (MOD(IOUT/32,2).EQ.0)  GO TO 100
      OPEN(UNIT=15, FILE='emode.d', STATUS='UNKNOWN')
      DO IM=1,NMODE
        COF(0,IM) = V(1,IM)
        IF (2*NPTS.EQ.ICYC)  COF(NPTS,IM) = V(ICYC,IM)
        DO K=1,(MPTS-1)/2
          KK = MPTS - K
          K1 = 2*K
          K2 = 2*K+1
          COF(K,IM) = (V(K1,IM) - JIMAG*V(K2,IM))/SQRT(2.0)
          COF(KK,IM) = CONJG(COF(K,IM))
        END DO
      END DO

      DO NM=1,NPRT
        DO I=1,NRES*LPRT+1
          EIGR(I) = 0.0
          EIGI(I) = 0.0
        END DO
        MN = MMODE-NM+1
        IN = NMDS(1,MN)
        IM = NMDS(2,MN)
        FACT = SQRT(2.0)
        IF (IN.EQ.0)  FACT = 1.0
        DO K=0,MPTS-1
          FRQ = (FLOAT(K)/FLOAT(ICYC) + FLOAT(IN)/FLOAT(NTOT))*TPI
          DO I=1,NRES*LPRT+1
            ANG = FRQ*FLOAT(I-1)/FLOAT(NRES)
            EIGR(I) = EIGR(I) + FACT*
     &        (REAL(COF(K,IM))*COS(ANG) - AIMAG(COF(K,IM))*SIN(ANG))
            EIGI(I) = EIGI(I) + FACT*
     &        (REAL(COF(K,IM))*SIN(ANG) + AIMAG(COF(K,IM))*COS(ANG))
          END DO
        END DO
        IF (IN.EQ.MHARM) THEN
          FRQ = FLOAT(IM*NHARM - MHARM)/FLOAT(NTOT)*TPI
          DO I=1,NRES*LPRT+1
            ANG = FRQ*FLOAT(I-1)/FLOAT(NRES)
            EIGR(I) = COS(ANG) + SIN(ANG)
          END DO
        END IF

        WRITE(15,'(6E13.5)')  (EIGR(I), I=1,NRES*LPRT+1)
C        IF (IN.NE.0 .AND. IN.NE.MHARM)
        WRITE(15,'(6E13.5)')  (EIGI(I), I=1,NRES*LPRT+1)
      END DO

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

      PARAMETER (NMAX=3000)
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
C          IF (J1.LT.1 .OR. J1.GT.NTOT)  GO TO 5
          IF (J1.LT.1) THEN
            JJ = MOD(J1+2*ICYC-1,ICYC)+1
            VAL = TSER(JJ)
          ELSE IF (J1.GT.NTOT) THEN
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

      SUBROUTINE FOURIER(NFPTS,NHS,NHE,Y,A,B)

      DIMENSION Y(NFPTS), A(0:NFPTS), B(0:NFPTS)

      DO 15 NFP=NHS,NHE
        COEF = 2.D0/FLOAT(NFPTS)
        CONST = 3.1415926536D0*COEF*NFP
        U1 = 0.D0
        U2 = 0.D0
        C = COS(CONST)
        S = SIN(CONST)
        I = NFPTS

5       U0 = Y(I) + 2.D0*C*U1 - U2
        U2 = U1
        U1 = U0
        I = I - 1
        IF (I-1) 10, 10, 5

10      A(NFP) = COEF*(Y(1) + C*U1 - U2)
        B(NFP) = COEF*S*U1
        IF (MOD(2*NFP,NFPTS).EQ.0)  A(NFP) = 0.5D0*A(NFP)
15    CONTINUE

      RETURN
      END

      SUBROUTINE SORT(N,RA,RB,INDX,JNDX)
      DIMENSION RA(N), RB(N), INDX(N), JNDX(N)

C       This program sorts N elements of an array RA in an ascending order 
C       using the heap sorting algorithm.

C ------- initialization
        L = N/2 + 1
        IR = N
      DO 10 I=1,N
        INDX(I) = I
10    CONTINUE

20    CONTINUE
C ------- hiring phase
      IF (L.GT.1)  THEN
        L = L-1
        KNDX = INDX(L)
        RRA = RA(KNDX)
C ------- retirement and promotion phase
      ELSE
        KNDX = INDX(IR)
        RRA = RA(KNDX)
        INDX(IR) = INDX(1)
        IR = IR-1
        IF (IR.EQ.1)  THEN
          INDX(1) = KNDX
          GO TO 40
        END IF
      END IF

        I = L
        J = L+L
30    CONTINUE
C ------- formation of heaps
      IF (J.LE.IR)  THEN
        J1 = MIN0(J+1,IR)
        IF (RA(INDX(J)).LT.RA(INDX(J1)))  J = J1
        IF (RRA.LT.RA(INDX(J)))  THEN
          INDX(I) = INDX(J)
          I = J
          J = I+I
        ELSE
          J = IR+1
        END IF
        GO TO 30
      END IF
        INDX(I) = KNDX
        GO TO 20

C ------- rank
40    DO 50 I=1,N
        JNDX(INDX(I)) = I
        RB(I) = RA(INDX(I))
50    CONTINUE

      RETURN
      END
