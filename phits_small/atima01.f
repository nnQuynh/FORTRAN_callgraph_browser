      SUBROUTINE APPDIF (ALDIF, ALPHA, XI, N, K, NCOMP, M, MSTAR)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ALDIF(*), ALPHA(*), XI(*), M(*)
      KD = K * NCOMP
      INCOMP = 0
      K3 = 0
      K4 = 0
      DO 130 J=1,NCOMP
           MJ = M(J)
           KMJ = K - MJ
           MJM1 = MJ - 1
           KMR = K + MJ
           NALPHJ = N * K + MJ
           INN = INCOMP
           K1 = MSTAR
           K2 = KD
           K5 = INN + 1
           DO 10 L=1,MJ
             ALDIF(K5) = ALPHA(K3+L)
   10      K5 = K5 + 1
           DO 50 I = 1, N
             IF (KMJ .EQ. 0)                        GO TO 30
             DO 20 L = 1, KMJ
               ALDIF(K5) = ALPHA(K1+K4+L)
   20        K5 = K5 + 1
   30        DO 40 L = 1, MJ
               ALDIF(K5) = ALPHA(K2+K3+L)
   40        K5 = K5 + 1
             K1 = K1 + KD
             K2 = K2 + KD
   50      CONTINUE
           IF (MJM1 .EQ. 0)                         GO TO 120
           DO 110 NR = 1, MJM1
             INN1 = INN + NALPHJ
             KMR = KMR - 1
             MJR = MJ - NR
             KMJR = K - MJR
             XIP1 = XI(1)
             DNK2 = DFLOAT(KMR) / (XI(2) - XIP1)
             DO 60 L=1,NR
   60        ALDIF(INN1+L) = 0.D0
             DO 70 L = NR, MJM1
               L1 = L + 1
   70        ALDIF(INN1+L) = 0.D0
   60        ALDIF(INN1+L) = 0.D0
             IBEG1 = MJ
             IBEG2 = K + NR
             DO 100 I = 1, N
               XII = XIP1
               XIP1 = XI(I+1)
               DNK1 = DFLOAT(KMR) / (XIP1 - XII)
               IF (I .LT. N) DNK2 = DFLOAT(KMR) / (XI(I+2) - XII)
               IF (I .EQ. N) DNK2 = DNK1
               DO 80 L = 1, KMJR
                 L1 = IBEG1 + L
   80          ALDIF(INN1+L1) = (ALDIF(INN+L1) - ALDIF(INN+L1-1)) 
               DO 90 L = 1, MJR
                 L1 = IBEG2 + L
   90          ALDIF(INN1+L1) = (ALDIF(INN+L1) - ALDIF(INN+L1-1)) 
               IBEG1 = IBEG1 + K
               IBEG2 = IBEG2 + K
  100        CONTINUE
             INN = INN1
  110      CONTINUE
  120      CONTINUE
           K3 = K3 + MJ
           K4 = K4 + KMJ
           INCOMP = INCOMP + NALPHJ * MJ
  130 CONTINUE
      RETURN
      END
      SUBROUTINE APPROX (I, X, Z, VN, XI, N, ALDIF, K, NCOMP,&
                 M, MSTAR, MODE, DMVAL, MODHI)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /NONLN/ PRECIS,NONLIN,ITER,LIMIT,ICARE,IPRINT,IGUESS,IFREEZ
!$OMP THREADPRIVATE(/NONLN/)
      COMMON /SIDE/  ZETA(40), ALEFT, ARIGHT, IZETA, IWR
!$OMP THREADPRIVATE(/SIDE/)
      DIMENSION Z(*), VN(*), XI(*), ALDIF(*), M(*), DMVAL(*)
      GO TO (10, 60, 70, 10, 80), MODE
   10 CONTINUE
      IF (X .GE. XI(1)-PRECIS .AND. X .LE. XI(N+1)+PRECIS)&
                                                   GO TO 20
      IF (IPRINT .LT. 1) WRITE(IWR,160) X, XI(1), XI(N+1)
      IF (X .LT. XI(1)) X = XI(1)
      IF (X .GT. XI(N+1)) X = XI(N+1)
   20 IF (I .GT. N .OR. I .LT. 1) I = (N+1) / 2
      ILEFT = I
      IF (X .LT. XI(ILEFT))                         GO TO 40
      DO 30 L=ILEFT,N
           I = L
           IF (X .LT. XI(L+1))                      GO TO 60
   30 CONTINUE
      GO TO 60
   40 IRIGHT = ILEFT - 1
      DO 50 L=1,IRIGHT
           I = IRIGHT + 1 - L
           IF (X .GE. XI(I))                        GO TO 60
   50 CONTINUE
   60 IF (MODE .EQ. 4)                              RETURN
      RHOX = (XI(I+1) - X) / (XI(I+1) - XI(I))
      CALL BSPFIX (RHOX, VN, K, NCOMP, M)
   70 CALL BSPVAR (I, X, VN, XI, N, K, NCOMP, M)
   80 DO 90 L=1,MSTAR
   90 Z(L) = 0.D0
      INDIF = 0
      K5 = 1
      IF (MODHI .EQ. 0)                             GO TO 110
      IVNHI = K * (K-1) / 2
      DNK2 = DFLOAT(K) / (XI(I+1) - XI(I))
      INCOMP = 0
      DO 100 J=1,NCOMP
  100 DMVAL(J) = 0.D0
  110 DO 150 J = 1, NCOMP
           MJ = M(J)
           NALPHJ = N * K + MJ
           KMR = K + MJ
           IVN = KMR * (KMR - 1) / 2
           DO 130 NR = 1, MJ
             LEFT = I * K + MJ - KMR
             DO 120 L = 1, KMR
               LEFTPL = LEFT + L
  120        Z(K5) = Z(K5) + ALDIF(INDIF+LEFTPL) * VN(IVN+L)
             KMR = KMR - 1
             IVN = IVN - KMR
             K5 = K5 + 1
  130      INDIF = INDIF + NALPHJ
           IF (MODHI .EQ. 0)                        GO TO 150
           INCOMP = INCOMP + (MJ-1) * NALPHJ
           LEFT = (I-1) * K + MJ - 1
           DO 140 L = 1, K
  140      DMVAL(J) = DMVAL(J) + DNK2 * (ALDIF(INCOMP+LEFT+L+1) -&
           ALDIF(INCOMP+LEFT+L)) * VN(IVNHI+L)
           INCOMP = INCOMP + NALPHJ
  150 CONTINUE
      RETURN
  160 FORMAT(37H ****** DOMAIN ERROR IN APPROX ******&
            /4H X =,D20.10, 10H   ALEFT =,D20.10,&
            11H   ARIGHT =,D20.10)
      END
      double precision function  bethek(zp,zt,mp,mt,energy,rho&
           ,gas,fntp,pot)
      implicit real*8 (a-h,o-z)
      include 'atimacnt.inc'
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      INTEGER*4 ZP,ZT,I,GAS
      INTEGER*4 ZPTMP
      REAL*8    MP,MT,ENERGY,RHO,FNTP,POT,SEOUT
      REAL*8    E,E1,GAMMA,BSQ,BETA,ZPEFF,ETA,ZETA
      REAL*8    F1,F2,F3,F4,F5,F6,DELTA
      REAL*8    BF4,bf1
      real*8    dl_lisoe
      external  dl_lisoe
      e = energy/mp
      e1 = energy/(1000.d0 * mp)
      zptmp = zp              ! S.Abe 2017/12/08
      gamma = 1.d0 + e1/AMU
      bsq   = 1.d0 -1.d0/(gamma**2)
      beta  = dsqrt(bsq)
      zeta = 1.d0 - exp( -130.d0 * beta / (iabs(zp) ** (2.d0/3.d0)) )   ! S.Abe 2017/12/08, zp is changed to iabs(zp)
      zpeff = zeta * dble(zp)
      eta   = beta * gamma
      f1 = 4.d-5*PI*ECHARGE2**2*N_AVO/EMASS * zpeff**2 * dble(zt)&
          / (bsq * mt)
      f2 =  dlog(2.d0 * emass * 1.d6 * bsq / pot)
      if  (eta .ge. 0.13d0) then
         c =   (.422377d0*eta**(-2) + .0304043d0*eta**(-4)&
            - .00038106d0*(eta**(-6)))*1.d-6*(pot**2)&
            + (3.858019d0*(eta**(-2)) - 0.1667989d0*(eta**(-4))&
            + 0.00157955d0*(eta**(-6))) * 1.0d-9*(pot**3)
         f2 = f2 - c / dble(zt)
      endif
      f4    = bf4(zt,eta,zpeff)
      f6    = 2.d0*dlog(gamma)-bsq
      call dichteeff(beta, zt, mt, pot, rho, delta)
      bethek  = f1 * ((f2+f6)*f4&
            + dl_lisoe(zp, mp, beta)-delta/2.d0)/1000.d0
      zp = zptmp    ! S.Abe 2017/12/08
      return
      end
      DOUBLE PRECISION FUNCTION BF4(ZT,ETA,ZPEFF)
      implicit real*8 (a-h,o-z)
      include 'atimacnt.inc'
      INTEGER*4 ZT
      INTEGER*4 I
      REAL*8 ETA
      REAL*8 ZPEFF
      REAL*8 V1
      REAL*8 V2FV
      REAL*8 VA(4)
      REAL*8 V2FVA(4)
      DATA V2FVA /0.33D0,0.30D0,0.26D0,0.23D0/
      DATA VA  /1.D0,2.D0,3.D0,4.D0/
      V1 = ETA/(ALPHA*DSQRT(DBLE(ZT)))
      IF (V1 .GE. 4.D0) THEN
        V2FV = 0.45D0/DSQRT(V1)
      ELSE IF ((V1 .GE. 1.D0).AND.(V1 .LT. 4.D0)) THEN
        DO 2,I = 1,4
          IF (VA(I) .GE. V1) GOTO 100
2       CONTINUE
100     V2FV = V2FVA(I-1)+(V1-VA(I-1))*(V2FVA(I)-V2FVA(I-1))&
                /(VA(I)-VA(I-1))
      ELSE
         V2FV = 0.0D0
         WRITE (*,*)'BARKAS REDUCED VELOCITY SMALLER THAN ONE'&
         ,V1, ETA, ALPHA,ZT
         STOP
      END IF
      BF4 = 1.D0+2.0D0 * ZPEFF * V2FV/((V1**2)*DSQRT(DBLE(ZT)))
      RETURN
      END
      SUBROUTINE BLDBLK (I, X, LL, Q, NROW, NC, Z, DF, NCOMP,&
                XI, ALPHO, IALPHO, MODE, DFSUB, DGSUB)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /APPR/ N,NOLD,NMAX,NALPHA,MSHFLG,MSHNUM,MSHLMT,MSHALT
!$OMP THREADPRIVATE(/APPR/)
      COMMON /ORDER/  K, ND, MSTAR, KD, KDM, MNSUM, M(20)
!$OMP THREADPRIVATE(/ORDER/)
      COMMON /SIDE/   ZETA(40), ALEFT, ARIGHT, IZETA, IWR
!$OMP THREADPRIVATE(/SIDE/)
      COMMON /NONLN/ PRECIS,NONLIN,ITER,LIMIT,ICARE,IPRINT,IGUESS,IFREEZ
!$OMP THREADPRIVATE(/NONLN/)
      COMMON /BSPLIN/  VNCOL(66,7), VNSAVE(66,5), VN(66)
!$OMP THREADPRIVATE(/BSPLIN/)
      DIMENSION Q(NROW,*), Z(*), DF(NCOMP,1)
      DIMENSION XI(*), BASEF(620), ALPHO(*), DG(40)
      NK = N      
      IF (MODE .EQ. 2)  NK = NC + NCOMP - 1
      DO 10 J=NC,NK
           DO 10 L=1,KDM
   10 Q(J,L)=0.D0
      GO TO (20, 130), MODE
   20 CALL BSPDER (VN, XI, N, X, I, BASEF, 2)
      CALL DGSUB (IZETA, Z, DG)
      IF (ITER .GE. 1 .OR. NONLIN .EQ. 0)           GO TO 40
      VALUE = 0.D0
      DO 30 J=1,MSTAR
   30 VALUE = VALUE + DG(J) * Z(J)
      ALPHO(IALPHO) =  VALUE
   40 IQ = 0
      IQM = MSTAR
      IDG = 0
      IBASEF = 0
      ID = N      
      DO 120 JCOMP=1,NCOMP
           MJ = M(JCOMP)
           MJ1 = MJ + 1
           KMJ = K - MJ
           DO 60 L=1,MJ
             DO 50 J=1,MJ
   50        Q(ID, IQ+L) = Q(ID, IQ+L) + DG(IDG+J) *&
            BASEF(IBASEF + J)
   60      IBASEF = IBASEF + MJ1
           IF (KMJ .LE. 0)                          GO TO 90
           DO 80 L=1,KMJ
             DO 70 J=1,MJ
   70        Q(ID, IQM+L) = Q(ID, IQM+L) + DG(IDG+J) *&
            BASEF(IBASEF+J)
   80      IBASEF = IBASEF + MJ1
   90      DO 110 L=1,MJ
             DO 100 J=1,MJ
  100        Q(ID, IQ+KD+L) = Q(ID, IQ+KD+L) + DG(IDG+J) *&
            BASEF(IBASEF+J)
  110      IBASEF = IBASEF + MJ1
           IDG =IDG + MJ
           IQ = IQ + MJ
           IQM = IQM + KMJ
  120 CONTINUE
      RETURN
  130 CALL BSPDER (VNCOL(1,LL), XI, N, X , I, BASEF, 3)
      CALL DFSUB (X, Z, DF)
      DO 240 JJ=1,NCOMP
           IF (ITER .GE. 1 .OR. NONLIN .EQ. 0)      GO TO 150
           IALPHO = IALPHO + 1
           VALUE = 0.D0
           DO 140 J=1,MSTAR
  140      VALUE = VALUE + DF(JJ,J) * Z(J)
           ALPHO(IALPHO) = ALPHO(IALPHO) - VALUE
  150      ID = JJ + NC - 1
           IQ=0
           IQM=MSTAR
           IDF=0
           IBASEF=0
           DO 230 JCOMP=1,NCOMP
             MJ = M(JCOMP)
             MJ1 = MJ + 1
             KMJ = K - MJ
             DO 170 L=1,MJ
               IF (JCOMP . EQ. JJ) Q(ID, IQ+L) = BASEF(IBASEF+MJ1)
               DO 160 J=1,MJ
  160          Q(ID, IQ+L) = Q(ID, IQ+L) - DF(JJ, IDF+J) * &
               BASEF(IBASEF+J)
  170        IBASEF = IBASEF + MJ1
             IF (KMJ .LE. 0)                        GO TO 200
             DO 190 L=1,KMJ
               IF (JCOMP .EQ. JJ) Q(ID, IQM+L) = BASEF(IBASEF+MJ1)
               DO 180 J=1,MJ
  180          Q(ID, IQM+L) = Q(ID, IQM+L) -&
              DF(JJ, IDF+J) * BASEF(IBASEF+J)
  190        IBASEF = IBASEF + MJ1
  200        DO 220 L=1,MJ
               IF (JCOMP .EQ. JJ) Q(ID, IQ+KD+L) = BASEF(IBASEF+MJ1)
               DO 210 J=1,MJ
  210          Q(ID, IQ+KD+L) = Q(ID, IQ+KD+L) - DF(JJ, IDF+J) *&
              BASEF(IBASEF+J)
  220        IBASEF = IBASEF + MJ1
             IDF = IDF + MJ
             IQ = IQ + MJ
             IQM = IQM + KMJ
  230      CONTINUE
  240 CONTINUE
      RETURN
      END
      SUBROUTINE BSCOEF(I,FSPACE,ISPACE,T,BCOEF,NB,KPM)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION FSPACE(*),ISPACE(*),T(*),BCOEF(*)
      N=ISPACE(1)
      K=ISPACE(2)
      M=ISPACE(7+I)
      KPM=K+M
      CALL KNOTS(FSPACE(1),N,M,KPM,T,NB)
      NCOMP=ISPACE(3)
      INCOMP=N+1
      DO 10 J=1,NCOMP
      MJ=ISPACE(7+J)
      NALPHJ=N*K+MJ
      IF(J.EQ.I) GO TO 20
   10 INCOMP=INCOMP+NALPHJ*MJ
      WRITE(6,15)
   15 FORMAT('  FEHLER IN BCOEF: I=',I2,' NCOMP=',I2)
      STOP
   20 DO 30 J=1,NB
   30 BCOEF(J)=FSPACE(INCOMP+J)
      RETURN
      END
      SUBROUTINE BSPDER (VN, XMESH, N, X, I, BASEF, MODE)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /ORDER/ K, NCOMP, MSTAR, KD, KDM, MNSUM, M(20)
!$OMP THREADPRIVATE(/ORDER/)
      COMMON /HI/    DN1, DN2, DN3
!$OMP THREADPRIVATE(/HI/)
      COMMON /EQORD/ IND(5), INEQ(20), MND(5), ND, NEQ
!$OMP THREADPRIVATE(/EQORD/)
      DIMENSION BASEF(*), VN(*), XMESH(*)
      DIMENSION ALPHD(80), ALPHDO(80), ALPHN(280) , ALPHNO(280)
      GO TO (10, 20, 30, 40), MODE
   10 XIL = XMESH(1)
      IF (I .GT. 1)  XIL = XMESH(I-1)
      XIR = XMESH(N+1)
      IF (I .LT. N)  XIR = XMESH(I+2)
      DN1 = 1.D0 / (XMESH(I+1) - XIL)
      DN2 = 1.D0 / (XMESH(I+1) - XMESH(I))
      DN3 = 1.D0 / (XIR - XMESH(I))
   20 RHOX = (XMESH(I+1) - X) * DN2
      CALL BSPFIX (RHOX, VN, K, NCOMP, M)
   30 CALL BSPVAR (I, X, VN, XMESH, N, K, NCOMP, M)
   40 MD = MND(ND)
      KMD = K + MD
      KMD1 = KMD + 1
      MD1 = MD + 1
      MD2M2 = MD * 2 - 2
      MD2M1 = MD2M2 + 1
      INL = KMD * 2
      DO 50 J=1,KMD
           ALPHDO (J)  = 0.D0
   50 ALPHDO(J+KMD) = 1.D0
      KUP = KMD * MD
      DO 60 J=1,KUP
   60 ALPHDO(J+INL) = 0.D0
      NDM1 = ND - 1
      NREST = MD2M2  - K
      INN = 0
      IF (NREST .LE. 0)                             GO TO 100
      IF (ND .EQ. 1)                                GO TO 100
      INL = 2 * MD2M2
      DO 90 NN = 1,NDM1
           MN2 = MND(NN) + 2
           DO 70 J = 1,MD2M2
             ALPHNO(J+INN) = 0.D0
   70      ALPHNO(J+INN+MD2M2) = 1.D0
           KUP = MD2M2 * MND(NN)
           DO 80 J=1,KUP
   80      ALPHNO(J+INN+INL) = 0.D0
   90 INN = INN + MN2 * MD2M2
  100 INNS = INN
      DO 120 J=1,ND
           K1 = IND(J)
           MJ = MND(J)
           KMJ = K + MJ
           MJ1 = MJ + 1
           IVN = KMJ * (KMJ-1) / 2
           DO 120 L=1,KMJ
             BASEF(K1) = VN(IVN+L)
             DO 110 JJ=1,MJ
  110        BASEF(K1+JJ) = 0.D0
  120 K1 = K1 + MJ1
      DO 310 NR=1,MD
           NR1 = NR + 1
           MDR = MD - NR
           K1 = IND(ND) + NR
           KMDR = K + MDR
           IVN = KMDR * (KMDR-1) /2
           IF (MDR .EQ. 0)                          GO TO 150
           DO 140 J=1,MDR
             JR = J + NR
             JIN = JR + NR1 * KMD
             JINK = JIN + K
             DO 130 L=J,JR
               JIN1 = JIN - KMD1
               JINK1 = JINK - KMD1
               ALPHD(JIN) = DN1 * (ALPHDO(JIN) - ALPHDO(JIN1))
               ALPHD(JINK) = DN3 * (ALPHDO(JINK) - ALPHDO(JINK1))
               IN = K1 + (L-1) * MD1
               BASEF(IN) = BASEF(IN) + ALPHD(JIN) * VN(IVN+J)
               IN = IN + K * MD1
               BASEF(IN) = BASEF(IN) + ALPHD(JINK) * VN(IVN+J+K)
               JIN = JIN - KMD
  130        JINK = JINK - KMD
  140      CONTINUE
  150      MDR1 = MDR + 1
           IF ( MDR1 .GT. K)                        GO TO 180
           DO 170 J = MDR1,K
             JR = J + NR
             JIN = JR + NR1 * KMD
             DO 160 L = J,JR
               JIN1 = JIN - KMD1
               ALPHD(JIN) = DN2 * (ALPHDO(JIN) - ALPHDO(JIN1))
               IN = K1 + (L-1) * MD1
               BASEF(IN) = BASEF(IN) + ALPHD(JIN) * VN(IVN+J)
  160        JIN = JIN - KMD
  170      CONTINUE
  180      CONTINUE
           IF (ND .EQ. 1)                           GO TO 230
           INN = INNS
           DO 220 NN=1,NDM1
             NJ = ND - NN
             MJ = MND(NJ)
             INN = INN - (MJ+2) * MD2M2
             IF (NR .GT. MJ)                        GO TO 230
             KMJR = K + MJ - NR
             K1 = IND(NJ)+ NR
             IVN = KMJR * (KMJR-1) / 2
             MJ1 = MJ + 1
             JR1 = KMJR - MD + 1
             JR1 = MIN0 (JR1, MD-1)
             DO 190 J=1,JR1
               JR = J + NR
               JIN = JR + NR1 * KMD + MD - MJ
               DO 190 L=J,JR
                 IN = K1 + (L-1) * MJ1
                 BASEF(IN) = BASEF(IN) + ALPHD(JIN) * VN(IVN+J)
                 JIN = JIN - KMD
  190        CONTINUE
             DO 200 J=MD,KMJR
               JR = J + NR
               JIN = JR + NR1 * KMD
               DO 200 L=J,JR
                 IN = K1 + (L-1) * MJ1
                 BASEF(IN) = BASEF(IN) + ALPHD(JIN) * VN(IVN+J)
                 JIN = JIN - KMD
  200        CONTINUE
             JR2 = MD2M2 - KMJR
             IF (JR2 .LE. 0)                        GO TO 220
             DO 210 JJ=1,JR2
               J = JJ + JR1
               JR = J + NR
               JIN = JR + NR1 * MD2M2 + INN
               DO 210 L=J,JR
                 JIN1 = JIN - MD2M1
                 ALPHN(JIN) = DN2 * (ALPHNO(JIN) - ALPHNO(JIN1))
                 IN = K1 + (L-1) * MJ1
                 BASEF(IN) = BASEF(IN) + ALPHN(JIN) * VN(IVN+J)
                 JIN = JIN - MD2M2
  210        CONTINUE
  220      CONTINUE
  230      CONTINUE
           IF (NR .EQ. MD)                          GO TO 300
           NR2 = NR + 2
           INJ = NR
           DO 240 L=2,NR2
             INJ = INJ + KMD
             DO 240 J=1,KMDR
  240      ALPHDO(J+INJ) = ALPHD(J+INJ)
           IF (ND .EQ. 1)                           GO TO 300
           IF (NREST .LE. 0)                        GO TO 300
           INN = 0
           DO 290 NN = 1,NDM1
             MN = MND(NN)
             IF (MN .LE. NR)                        GO TO 280
             KMNR = K + MN - NR
             JR1 = MIN0 (KMNR-MD+1, MD-1)
             INJ = NR + INN
             INL = NR + MD - MN
             DO 250 L=2,NR2
               INJ = INJ + MD2M2
               INL = INL + KMD
               DO 250 J=1,JR1
  250        ALPHNO(INJ+J) = ALPHD(INL+J)
             MUP = MIN0 (KMNR, MD2M2)
             INJ = NR + INN
             INL = NR
             DO 260 L=2,NR2
               INJ = INJ + MD2M2
               INL = INL + KMD
               DO 260 J =MD,MUP
  260        ALPHNO(INJ+J) = ALPHD(INL+J)
             JR2 = MD2M2 - KMNR
             IF (JR2 .LE. 0)                        GO TO 280
             INJ = NR + INN
             DO 270 L=2,NR2
               INJ = INJ + MD2M2
               DO 270 JJ = 1,JR2
                 JIN = INJ + JJ + JR1
  270        ALPHNO(JIN) = ALPHN(JIN)
  280        INN = INN + (MN+2) * MD2M2
  290      CONTINUE
  300      CONTINUE
  310 CONTINUE
      DO 320 J=1,ND
           IN = IND(J)
           ICONS = 1
           MJ = MND(J)
           KMJ = K + MJ
           MJ1 = MJ + 1
           DO 320 NR = 1,MJ
             ICONS = ICONS * (KMJ-NR)
             IN = IN + 1
             DO 320 L=1,KMJ
               LBASEF = IN + (L-1) * MJ1
               BASEF(LBASEF) = BASEF(LBASEF) * DFLOAT(ICONS)
  320 CONTINUE
      IF (NEQ .EQ. 0)                               RETURN
      JD = 1
      DO 360 J=1,NEQ
           IN1 = INEQ(J)
  330      IF (IN1 .LT. IND(JD+1))                  GO TO 340
           JD = JD + 1
           GO TO 330
  340      MJ = MND(JD)
           NTOT = (K+MJ)*(1+MJ)
           IN2 = IND(JD)
           DO 350 L=1,NTOT
  350      BASEF(IN1-1+L) = BASEF(IN2-1+L)
  360 CONTINUE
      RETURN
      END
      SUBROUTINE BSPFIX (RHOX, VN, K, NCOMP, M)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VN(*), M(*)
      XRHO = 1.D0 - RHOX
      IVN = 0
      VN(1) = 1.D0
      DO 20 L=1,K
           IVN = IVN + L
           VNP = 0.D0
           DO 10 J=1,L
             REP =  VN(IVN-L+J)
             VN(IVN+J) = VNP + REP * RHOX
   10      VNP = REP * XRHO
   20 VN(IVN+L+1) = VNP
      MD1 = M(NCOMP) - 1
      IF (MD1 .LE. 0)                               RETURN
      DO 40 L=1,MD1
           IVN = IVN + K + L
           INC = L + 2
           VNP = VN(IVN+1-K) * XRHO
           IF (K .LT. INC)                          RETURN
           DO 30 J=INC,K
             REP = VN(IVN-L-K+J)
             VN(IVN+J) = VNP + REP  * RHOX
   30      VNP = REP * XRHO
   40 VN(IVN+K+1) = VNP
      RETURN
      END
      SUBROUTINE BSPVAR (I, X, VN, XI, N, K, NCOMP, M)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION VN(*), XI(*), M(*)
      MD1 = M(NCOMP) -1
      IF(MD1 .LE. 0)                                RETURN
      XIL = XI(1)
      IF (I .GT. 1) XIL = XI(I-1)
      XIR = XI(N+1)
      IF (I .LT. N) XIR = XI(I+2)
      RHO1 = (XI(I+1) - X) / (XI(I+1) - XI(I))
      RHO2 = (XI(I+1) - X) / (XI(I+1) - XIL)
      RHO3 = (XIR - X) / (XIR - XI(I))
      XRHO1 = 1.D0 - RHO1
      XRHO2 = 1.D0 - RHO2
      XRHO3 = 1.D0 - RHO3
      IVN = K * (K+1) / 2
      DO 30 L=1,MD1
           IVN = IVN + K + L
           VNP = 0.D0
           DO 10 J=1,L
             REP = VN(IVN-L-K+J)
             VN(IVN+J) = VNP + REP  * RHO2
   10      VNP = REP * XRHO2
           VN(IVN+L+1) = VNP + RHO1 * VN(IVN-K+1)
           VNP = VN(IVN-L) * XRHO1
           DO 20 J=1,L
             REP = VN(IVN+J-L)
             VN(IVN+K+J) = VNP + REP * RHO3
   20      VNP = REP * XRHO3
   30 VN(IVN+K+L+1) = VNP
      RETURN
      END
      real*8 function bvalue ( t, bcoef, n, k, x, jderiv )
      implicit real*8 (a-h,o-z)
      integer kmax
      parameter (kmax = 20)
      integer*4 jderiv,k,n,i,ilo,imk,j,jc,jcmin,jcmax,jj,kmj,km1,mflag&
                          ,nmi,jdrvp1
      real*8 bcoef(n),t,x,   aj(kmax),dl(kmax),dr(kmax),fkmj
      dimension t(n+k)
      bvalue = 0.
      if (jderiv .ge. k)                go to 99
      call interv ( t, n+k, x, i, mflag )
      if (mflag .ne. 0)                 go to 99
      km1 = k - 1
      if (km1 .gt. 0)                   go to 1
      bvalue = bcoef(i)
                                        go to 99
    1 jcmin = 1
      imk = i - k
      if (imk .ge. 0)                   go to 8
      jcmin = 1 - imk
      do 5 j=1,i
    5    dl(j) = x - t(i+1-j)
      do 6 j=i,km1
         aj(k-j) = 0.
    6    dl(j) = dl(i)
                                        go to 10
    8 do 9 j=1,km1
    9    dl(j) = x - t(i+1-j)
   10 jcmax = k
      nmi = n - i
      if (nmi .ge. 0)                   go to 18
      jcmax = k + nmi
      do 15 j=1,jcmax
   15    dr(j) = t(i+j) - x
      do 16 j=jcmax,km1
         aj(j+1) = 0.
   16    dr(j) = dr(jcmax)
                                        go to 20
   18 do 19 j=1,km1
   19    dr(j) = t(i+j) - x
   20 do 21 jc=jcmin,jcmax
   21    aj(jc) = bcoef(imk + jc)
      if (jderiv .eq. 0)                go to 30
      do 23 j=1,jderiv
         kmj = k-j
         fkmj = float(kmj)
         ilo = kmj
         do 23 jj=1,kmj
            aj(jj) = ((aj(jj+1) - aj(jj))/(dl(ilo) + dr(jj)))*fkmj
   23       ilo = ilo - 1
   30 if (jderiv .eq. km1)              go to 39
      jdrvp1 = jderiv + 1
      do 33 j=jdrvp1,km1
         kmj = k-j
         ilo = kmj
         do 33 jj=1,kmj
            aj(jj) = (aj(jj+1)*dl(ilo) + aj(jj)*dr(jj))/(dl(ilo)+dr(jj))
   33       ilo = ilo - 1
   39 bvalue = aj(1)
   99                                   return
      end
      subroutine calcint
      implicit real*8 (a-h,o-z)
      parameter ( kspc=10000, ispc=1000 )
      include 'atimasys.inc'
      include 'atimadim.inc'
      dimension fspace(kspc), ispace(ispc)
      external ffrange, ffrstr, ffastr, fftof
        iprint = 1
        ndimf  = ksp        
        ndimi  = isp
            a = 1.d-3
            b = 4.5d+5
            c = a
            tol = 1.d-6
            call integr(ffrange,a,b,c,tol,iprint,fspace,ispace,&
                       ndimf,ndimi,iflag)
            call bscoef(1,fspace,ispace,t1,bcoef1,nb1,kpm1)
            a = 1.d-3
            b = 4.5d+5
            c = a
            tol = 1.d-6
            call integr(ffrstr,a,b,c,tol,iprint,fspace,ispace,&
                      ndimf,ndimi,iflag)
            call bscoef(1,fspace,ispace,t2,bcoef2,nb2,kpm2)
            a = 1.d-3
            b = 4.5d+5
            c = b
            tol = 1.d-8
            call integr(ffastr,a,b,c,tol,iprint,fspace,ispace,&
                        ndimf,ndimi,iflag)
            call bscoef(1,fspace,ispace,t3,bcoef3,nb3,kpm3)
            a = 1.d-3
            b = 4.5d+5
            c = a
            tol = 1.d-5
            call integr(fftof,a,b,c,tol,iprint,fspace,ispace,&
                      ndimf,ndimi,iflag)
            call bscoef(1,fspace,ispace,t4,bcoef4,nb4,kpm4)
      return
      end
      double precision function frange(x)
      implicit real*8 (a-h,o-z)
      include 'atimasys.inc'
      double precision x
      double precision dedxver
      frange = mp / dedxver(x)
      return
      end
      subroutine  ffrange(x,z,f)
      implicit real*8 (a-h,o-z)
      double precision x, z(1), f(1)
      double precision frange
      f(1) = frange(x)
      return
      end
      double precision function frstragg(x)
      implicit real*8 (a-h,o-z)
      include 'atimasys.inc'
      double precision x
      double precision dedxver,enlostv
      frstragg = mp * enlostv(x) / dedxver(x)**3
      return
      end
      subroutine  ffrstr(x,z,f)
      implicit real*8 (a-h,o-z)
      double precision x, z(1), f(1)
      double precision frstragg
      f(1) = frstragg(x)
      return
      end
      double precision function fastragg(x)
      implicit real*8 (a-h,o-z)
      include 'atimasys.inc'
      double precision x
      double precision dedxver, omthver
      fastragg = mp * omthver(x) / dedxver(x)
      return
      end
      subroutine  ffastr(x,z,f)
      implicit real*8 (a-h,o-z)
      double precision x, z(1), f(1)
      double precision fastragg
      f(1) = fastragg(x)
      return
      end
      double precision function ftof(x)
      implicit real*8 (a-h,o-z)
      include 'atimasys.inc'
      include 'atimacnt.inc'
      double precision x
      double precision dedxver, frange
      ftof =  frange(x)/(CLIGHT*rho*1000*sqrt(1-1/(x/AMU + 1)**2))
      return
      end
      subroutine  fftof(x,z,f)
      implicit real*8 (a-h,o-z)
      double precision x, z(1), f(1)
      double precision ftof
      f(1) = ftof(x)
      return
      end
      SUBROUTINE COLSYS (NCOMP, M, ALEFT, ARIGHT, ZETA, IPAR, LTOL,&
                TOL, FIXPNT, ISPACE, FSPACE, IFLAG, FSUB,&
                DFSUB, GSUB, DGSUB, SOLUTN)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /ORDER/ K,NC,MSTAR,KD,KDM,MNSUM,MT(20)
!$OMP THREADPRIVATE(/ORDER/)
      COMMON /APPR/ N,NOLD,NMAX,NALPHA,MSHFLG,MSHNUM,MSHLMT,MSHALT
!$OMP THREADPRIVATE(/APPR/)
      COMMON /SIDE/ TZETA(40),TLEFT,TRIGHT,IZETA,IWR
!$OMP THREADPRIVATE(/SIDE/)
      COMMON /NONLN/ PRECIS,NONLIN,ITER,LIMIT,ICARE,IPRINT,IGUESS,IFREEZ
!$OMP THREADPRIVATE(/NONLN/)
      COMMON /EQORD/  IND(5), INEQ(20), MND(5), ND, NEQ
!$OMP THREADPRIVATE(/EQORD/)
      COMMON /ERRORS/ TTL(40),WGTMSH(40),TOLIN(40),ROOT(40),&
            JTOL(40),LTTOL(40),NTOL
!$OMP THREADPRIVATE(/ERRORS/)
      EXTERNAL FSUB, DFSUB, GSUB, DGSUB, SOLUTN
      DIMENSION M(*), ZETA(*), IPAR(*), LTOL(*), TOL(*),&
          FIXPNT(*), ISPACE(*), FSPACE(*)
      IWR = 6
      PRECIS = 1.D0
   10 PRECIS = PRECIS / 2.D0
      PRECP1 = PRECIS + 1.D0
      IF (PRECP1 .GT. 1.D0)                         GO TO 10
      PRECIS = PRECIS * 100.D0
      IFLAG = -3
      IF (NCOMP .LT. 1 .OR. NCOMP .GT. 20)          RETURN
      IF (M(1) .LT. 1 .OR. M(NCOMP) .GT. 4)         RETURN
      IF (NCOMP .EQ. 1)                             GO TO 30
      DO 20 I=2,NCOMP
           IF (M(I-1) .GT. M(I))                    RETURN
   20 CONTINUE
   30 CONTINUE
      NONLIN = IPAR(1)
      K = IPAR(2)
      IF (K .EQ. 0)   K = MAX0( M(NCOMP)+1, 5-M(NCOMP) )
      N = IPAR(3)
      IF (N .EQ. 0)  N = 5
      IREAD = IPAR(8)
      IGUESS = IPAR(9)
      IF (NONLIN .EQ. 0 .AND. IGUESS .EQ. 1) IGUESS = 0
      IF (IGUESS .GE. 2 .AND. IREAD .EQ. 0)  IREAD = 1
      ICARE = IPAR(10)
      NTOL = IPAR(4)
      NDIMF = IPAR(5)
      NDIMI = IPAR(6)
      NFXPNT = IPAR(11)
      IPRINT = IPAR(7)
      MSTAR = 0
      MNSUM = 0
      DO  40 I=1,NCOMP
           MNSUM = MNSUM + M(I)**2
   40 MSTAR = MSTAR + M(I)
      DO 50 I=1,NCOMP
   50 MT(I) = M(I)
      DO 60 I=1,MSTAR
   60 TZETA(I) = ZETA(I)
      DO 70 I=1,NTOL
           LTTOL(I) = LTOL(I)
   70 TOLIN(I) = TOL(I)
      TLEFT = ALEFT
      TRIGHT = ARIGHT
      NC = NCOMP
      KD = K * NCOMP
      KDM = KD + MSTAR
      IF (IPRINT .GT. (-1))                         GO TO 100
      IF (NONLIN .GT. 0)                            GO TO 80
      WRITE (IWR,260) NCOMP, (M(IP), IP=1,NCOMP)
      GO TO 90
   80 WRITE(IWR,270) NCOMP, (M(IP), IP=1,NCOMP)
   90 WRITE (IWR,280) (ZETA(IP), IP=1,MSTAR)
      WRITE (IWR,290) K
      WRITE (IWR,300) (LTOL(IP), IP=1,NTOL)
      WRITE (IWR,310) (TOL(IP), IP=1,NTOL)
      IF (IGUESS .GE. 2) WRITE (IWR,320)
      IF (IREAD .EQ. 2) WRITE (IWR,330)
      IF (NFXPNT .GT. 0) WRITE (IWR,340) NFXPNT,&
           (FIXPNT(IP), IP=1,NFXPNT)
  100 CONTINUE
      IF (K .LT. 0 .OR. K .GT. 7)                   RETURN
      IF (N .LT. 0)                                 RETURN
      IF (IREAD .LT. 0 .OR. IREAD .GT. 2)           RETURN
      IF (IGUESS .LT. 0 .OR. IGUESS .GT. 4)         RETURN
      IF (ICARE .LT. 0 .OR. ICARE .GT. 2)           RETURN
      IF (NTOL .LT. 0 .OR. NTOL .GT. MSTAR)         RETURN
      IF (NFXPNT .LT. 0)                            RETURN
      IF (IPRINT .LT. (-1) .OR. IPRINT .GT. 1)      RETURN
      IF (MSTAR .LT. 0 .OR. MSTAR .GT. 40)          RETURN
      MSHLMT = 3
      MSHFLG = 0
      MSHNUM = 1
      MSHALT = 1
      LIMIT = 40
      NREC = 0
      DO 110 II=1,MSTAR
           I = MSTAR + 1 - II
           IF (ZETA(I) .LT. ARIGHT)                 GO TO 110
           NREC = II
  110 CONTINUE
      NFIXI = NRE      
      NSIZEI = 3 + KDM - NRE      
      NFIXF = NREC * (KDM+1) + 2 * MNSUM + 2 * MSTAR + 3
      NSIZEF = 4 + K + 2 * KD + (4+2*K) * MSTAR +&
     (KDM-NREC) * (KDM+1)
      NMAXF = (NDIMF - NFIXF) / NSIZEF
      NMAXI = (NDIMI - NFIXI) / NSIZEI
      IF (IPRINT .LT. 1) WRITE(IWR,350) NMAXF, NMAXI
      NMAX = MIN0(NMAXF,NMAXI)
      IF (NMAX .LT. N)                              RETURN
      IF (NMAX .LT. NFXPNT+1)                       RETURN
      IF (NMAX .LT. 2*NFXPNT+2  .AND.  IPRINT .LT. 1) WRITE(IWR,360)
      LXI = 1
      LA = LXI + NMAX + 1
      LXIOLD = LA + KDM * (NMAX * (KDM-NREC) + NREC)
      LXIJ = LXIOLD + NMAX + 1
      LALPHA = LXIJ + K * NMAX
      LDLPHA = LALPHA + NMAX * KD + MSTAR
      LELPHA = LDLPHA + NMAX * KD + MSTAR
      LALDIF = LELPHA + NMAX * K * MSTAR + MNSUM
      LRHS = LALDIF + NMAX * K * MSTAR + MNSUM
      LVALST = LRHS + NMAX * (KDM - NREC) + NRE      
      LSLOPE = LVALST + 4 * MSTAR * NMAX
      LACCUM = LSLOPE + NMAX
      LIPIV = 1
      LINTEG = LIPIV + (LVALST - LRHS)
      IF (IGUESS .LT. 2)                            GO TO 160
      NOLD = N
      IF (IGUESS .EQ. 4)  NOLD = ISPACE(1)
      NALDIF =  NOLD * K * MSTAR + MNSUM
      NP1 = N + 1
      IF (IGUESS .EQ. 4)  NP1 = NP1 + NOLD + 1
      DO 120 I=1,NALDIF
  120 FSPACE( LALDIF+I-1 )  =  FSPACE( NP1+I )
      NP1 = NOLD + 1
      IF (IGUESS .EQ. 4)                            GO TO 140
      DO 130 I=1,NP1
  130 FSPACE( LXIOLD+I-1 )  =  FSPACE( LXI+I-1 )
      GO TO 160
  140 DO 150 I=1,NP1
  150 FSPACE( LXIOLD+I-1 )  =  FSPACE( N+1+I )
  160 CONTINUE
      CALL CONSTS
      CALL NEWMSH (3+IREAD, FSPACE(LXI), FSPACE(LXIOLD),&
          FSPACE(LXIJ), DUM1, DUM2, DUM3, DUM4,&
          NFXPNT, FIXPNT)
      IND(1) = 1
      MND(1) = M(1)
      ND = 1
      NEQ = 0
      IG = (M(1)+1) * (M(1)+K) + 1
      IF (NCOMP .LE. 1)                             GO TO 200
      DO 190 J=2,NCOMP
           MJ = M(J)
           IF (MJ .EQ. M(J-1))                      GO TO 170
           ND = ND + 1
           IND(ND) = IG
           MND(ND) = MJ
           GO TO 180
  170      NEQ = NEQ + 1
           INEQ(NEQ) = IG
  180      IG = IG + (MJ+1) * (MJ+K)
  190 CONTINUE
      IND(ND+1) =IND(ND) + IG
  200 CONTINUE
      IF (IGUESS .GE. 2)                            GO TO 230
      NP1 = N + 1
      DO 210 I = 1,NP1
  210 FSPACE(I + LXIOLD - 1) = FSPACE(I + LXI - 1)
      NOLD = N
      IF (NONLIN .EQ. 0 .OR. IGUESS .EQ. 1)         GO TO 230
      DO 220 I = 1,NALPHA
  220 FSPACE(I + LALPHA - 1) = 0.D0
      CALL APPDIF (FSPACE(LALDIF), FSPACE(LALPHA), FSPACE(LXI),&
          N, K, NC, MT, MSTAR)
  230 CONTINUE
      IF (IGUESS .GE. 2)  IGUESS = 0
      CALL CONTRL (FSPACE(LXI),FSPACE(LXIOLD),FSPACE(LXIJ),&
          FSPACE(LALPHA),FSPACE(LALDIF),FSPACE(LRHS),&
          FSPACE(LDLPHA), FSPACE(LELPHA),&
          FSPACE(LA),FSPACE(LVALST),FSPACE(LSLOPE),&
          FSPACE(LACCUM),ISPACE(LIPIV),ISPACE(LINTEG),&
          NFXPNT,FIXPNT,IFLAG,FSUB,DFSUB,GSUB,DGSUB,&
          SOLUTN)
      ISPACE(1) = N
      ISPACE(2) = K
      ISPACE(3) = NCOMP
      ISPACE(4) = MSTAR
      NALDIF = N * K * MSTAR + MNSUM
      ISPACE(5) = NALDIF
      ISPACE(6) = NALDIF + N + 2
      ISPACE(7) = ISPACE(6) + 65
      DO 240 I=1,NCOMP
  240 ISPACE(7+I) = M(I)
      DO 250 I=1,NALDIF
  250 FSPACE(N+1+I) = FSPACE(LALDIF-1+I)
      RETURN
  260 FORMAT(/// 37H THE NUMBER OF (LINEAR) DIFF EQNS IS , I3/ 1X,&
            16HTHEIR ORDERS ARE, 20I3)
  270 FORMAT(/// 40H THE NUMBER OF (NONLINEAR) DIFF EQNS IS , I3/ 1X,&
            16HTHEIR ORDERS ARE, 20I3)
  280 FORMAT(27H SIDE CONDITION POINTS ZETA, 8F10.6, 4( / 27X, 8F10.6))
  290 FORMAT(37H NUMBER OF COLLOC PTS PER INTERVAL IS, I3)
  300 FORMAT(39H COMPONENTS OF Z REQUIRING TOLERANCES -,8(7X,I2,1X),&
            4(/38X,8I10))
  310 FORMAT(33H CORRESPONDING ERROR TOLERANCES -,6X,8D10.2,&
            4(/39X,8D10.2))
  320 FORMAT(44H INITIAL MESH(ES) AND ALPHA PROVIDED BY USER)
  330 FORMAT(27H NO ADAPTIVE MESH SELECTION)
  340 FORMAT(10H THERE ARE ,I5,27H FIXED POINTS IN THE MESH - ,&
            10(6D12.4/))
  350 FORMAT(44H THE MAXIMUM NUMBER OF SUBINTERVALS IS MIN (, I4,&
            23H (ALLOWED FROM FSPACE),,I4, 24H (ALLOWED FROM ISPACE) ))
  360 FORMAT(/53H INSUFFICIENT SPACE TO DOUBLE MESH FOR ERROR ESTIMATE)
      END
      SUBROUTINE CONSTS
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /COLLOC/ RHO(7),WGTERR(40)
!$OMP THREADPRIVATE(/COLLOC/)
      COMMON /ORDER/ K,NCOMP,MSTAR,KD,KDM,MNSUM,M(20)
!$OMP THREADPRIVATE(/ORDER/)
      COMMON /BSPLIN/ VNCOL(66,7), VNSAVE(66,5), VN(66)
!$OMP THREADPRIVATE(/BSPLIN/)
      COMMON /ERRORS/ TOL(40),WGTMSH(40),TOLIN(40),ROOT(40),&
            JTOL(40),LTOL(40),NTOL
!$OMP THREADPRIVATE(/ERRORS/)
      COMMON /NONLN/ PRECIS,NONLIN,ITER,LIMIT,ICARE,IPRINT,IGUESS,IFREEZ
!$OMP THREADPRIVATE(/NONLN/)
      DIMENSION CNSTS1(28), CNSTS2(28)
      DATA CNSTS1 / .25D0,   .625D-1, .72169D-1,    1.8342D-2,&
          1.9065D-2, 5.8190D-2,   5.4658D-3, 5.3370D-3, 1.8890D-2,&
          2.7792D-2,   1.6095D-3, 1.4964D-3, 7.5938D-3, 5.7573D-3,&
          1.8342D-2,  4.673D-3, 4.150D-4, 1.919D-3, 1.468D-3,&
          6.371D-3, 4.610D-3,  1.342D-4, 1.138D-4, 4.889D-4,&
          4.177D-4, 1.374D-3, 1.654D-3, 2.863D-3  /
      DATA CNSTS2 / 1.25D-1,    2.604D-3,  8.019D-3,  2.170D-5,&
          7.453D-5, 5.208D-4, 9.689D-8, 3.689D-7, 3.100D-6, 2.451D-5,&
          2.691D-10, 1.120D-9,   1.076D-8,  9.405D-8,  1.033D-6,&
          5.097D-13, 2.290D-12, 2.446D-11, 2.331D-10, 2.936D-9,&
          3.593D-8,  7.001D-16, 3.363D-15, 3.921D-14, 4.028D-13,&
          5.646D-12, 7.531D-11, 1.129D-9  /
      KOFF = K * ( K + 1 ) / 2
      IEXTRA = 1
      DO 10 J = 1,NCOMP
           LIM = M(J)
           DO 10 L = 1,LIM
             WGTERR(IEXTRA) = CNSTS1(KOFF - LIM + L)
             IEXTRA = IEXTRA + 1
   10 CONTINUE
      JCOMP = 1
      MTOT = M(1)
      DO 40 I=1,NTOL
           LTOLI = LTOL(I)
   20      CONTINUE
           IF (LTOLI .LE. MTOT)                     GO TO 30
           JCOMP = JCOMP + 1
           MTOT = MTOT + M(JCOMP)
           GO TO 20
   30      CONTINUE
           JTOL(I) = JCOMP
           WGTMSH(I) = 1.D1 * CNSTS2(KOFF+LTOLI-MTOT) / TOLIN(I)
           ROOT(I) = 1.D0 / DFLOAT(K+MTOT-LTOLI+1)
   40 CONTINUE
      GO TO (50,60,70,80,90,100,110),K
   50 RHO(1) = 0.D0
      GO TO 120
   60 RHO(2) = .57735026918962576451D0
      RHO(1) = - RHO(2)
      GO TO 120
   70 RHO(3) = .77459666924148337704D0
      RHO(2) = .0D0
      RHO(1) = - RHO(3)
      GO TO 120
   80 RHO(1) = -.86113631159405257523D0
      RHO(2) = -.33998104358485626480D0
      RHO(3) = - RHO(2)
      RHO(4) = - RHO(1)
      GO TO 120
   90 RHO(5) = .90617984593866399280D0
      RHO(4) = .53846931010568309104D0
      RHO(3) = .0D0
      RHO(2) = - RHO(4)
      RHO(1) = - RHO(5)
      GO TO 120
  100 RHO(6) = .93246951420315202781D0
      RHO(5) = .66120938646626451366D0
      RHO(4) = .23861918608319690863D0
      RHO(3) = -RHO(4)
      RHO(2) = -RHO(5)
      RHO(1) = -RHO(6)
      GO TO 120
  110 RHO(7) = .949107991234275852452D0
      RHO(6) = .74153118559939443986D0
      RHO(5) = .40584515137739716690D0
      RHO(4) = 0.D0
      RHO(3) = -RHO(5)
      RHO(2) = -RHO(6)
      RHO(1) = -RHO(7)
  120 CONTINUE
      DO 130 J=1,K
           ARG = .5D0 * (1.D0 - RHO(J))
           CALL BSPFIX (ARG, VNCOL(1,J), K, NCOMP, M)
  130 CONTINUE
      CALL BSPFIX (1.D0, VNSAVE(1,1), K, NCOMP, M)
      CALL BSPFIX (5.D0/6.D0, VNSAVE(1,2), K, NCOMP, M)
      CALL BSPFIX (2.D0/3.D0, VNSAVE(1,3), K, NCOMP, M)
      CALL BSPFIX (1.D0/3.D0, VNSAVE(1,4), K, NCOMP, M)
      CALL BSPFIX (1.D0/6.D0, VNSAVE(1,5), K, NCOMP, M)
      RETURN
      END
      SUBROUTINE CONTRL(XI, XIOLD, XIJ, ALPHA, ALDIF, RHS,&
                DALPHA, EALPHA, A, VALSTR, SLOPE,&
                ACCUM, IPIV, INTEGS, NFXPNT, FIXPNT, IFLAG,&
                FSUB, DFSUB, GSUB, DGSUB, SOLUTN)
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL FSUB, DFSUB, GSUB, DGSUB, SOLUTN
      DIMENSION XI(*), XIOLD(*), XIJ(*), ALPHA(*), ALDIF(*), RHS(*)
      DIMENSION A(*), VALSTR(*), SLOPE(*), ACCUM(*), IPIV(*), INTEGS(*)
      DIMENSION DALPHA(*), EALPHA(*) , FIXPNT(*)
      COMMON /ORDER/ K,NCOMP,MSTAR,KD,KDM,MNSUM,M(20)
!$OMP THREADPRIVATE(/ORDER/)
      COMMON /APPR/ N,NOLD,NMAX,NALPHA,MSHFLG,MSHNUM,MSHLMT,MSHALT
!$OMP THREADPRIVATE(/APPR/)
      COMMON /SIDE/ ZETA(40),ALEFT,ARIGHT,IZETA,IWR
!$OMP THREADPRIVATE(/SIDE/)
      COMMON /NONLN/ PRECIS,NONLIN,ITER,LIMIT,ICARE,IPRINT,IGUESS,IFREEZ
!$OMP THREADPRIVATE(/NONLN/)
      COMMON /EQORD/  IND(5), INEQ(20), MND(5), ND, NEQ
!$OMP THREADPRIVATE(/EQORD/)
      COMMON /ERRORS/ TOL(40),WGTMSH(40),TOLIN(40),ROOT(40),&
            JTOL(40),LTOL(40),NTOL
!$OMP THREADPRIVATE(/ERRORS/)
      RELMIN = 1.D-3
      RSTART = 1.D-2
      LMTFRZ = 4
      CHECK = 0.D0
      DO 10 I=1,NTOL
   10 CHECK = DMAX1 (TOLIN(I), CHECK )
      FALPHA = DFLOAT(NALPHA)
      IMESH = 1
      ICON = 0
      IF (NONLIN .EQ. 0) ICON=1
      ICOR = 0
      LCONV = 0
   20      CONTINUE
           ITER = 0
           NALDIF = N * K * MSTAR + MNSUM
           IF (NONLIN .GT. 0)                       GO TO 60
           CALL LSYSLV (IFLAG, XI, XIOLD, XIJ, ALPHA, ALDIF, RHS,&
               EALPHA, A, IPIV, INTEGS, RNORM, 0, FSUB,&
               DFSUB, GSUB, DGSUB, SOLUTN)
           IF (IFLAG .NE. 0)                        GO TO 40
   30      IF (IPRINT .LT. 1) WRITE (IWR,490 )
           RETURN
   40      NP1 = N + 1
           DO 50 I=1,NP1
   50      XIOLD(I) = XI(I)
           NOLD = N
           CALL APPDIF (ALDIF, ALPHA, XI, N, K, NCOMP, M, MSTAR)
           GO TO 450
   60      RELAX = 1.D0
           IF (ICARE .EQ. 1 .OR. ICARE .EQ. (-1))  RELAX = RSTART
           IF (ICON .EQ. 0)                         GO TO 140
           IFREEZ = 0
           CALL LSYSLV (IFLAG, XI, XIOLD, XIJ, DALPHA, ALDIF, RHS,&
               ALPHA, A, IPIV, INTEGS, RNORM, 1, FSUB,&
               DFSUB, GSUB, DGSUB, SOLUTN)
           IF (IPRINT .LT. 0  .AND.  ITER .EQ. 0)  WRITE(IWR,530)
   70      IF (IPRINT .LT. 0)  WRITE (IWR,510) ITER, RNORM
           RNOLD = RNORM
           CALL LSYSLV (IFLAG, XI, XIOLD, XIJ, DALPHA, ALDIF, RHS,&
               ALPHA, A, IPIV, INTEGS, RNORM, 2+IFREEZ ,&
               FSUB, DFSUB, GSUB, DGSUB, SOLUTN)
           IF (IFLAG .EQ. 0)                        GO TO 30
           IF (IFREEZ .EQ. 1)                       GO TO 90
           ITER = ITER + 1
           IFRZ = 0
           NP1 = N + 1
           DO 80 I=1,NP1
   80      XIOLD(I) = XI(I)
           NOLD = N
   90      CONTINUE
           DO 100 I=1,NALPHA
  100      ALPHA(I) = ALPHA(I) + DALPHA(I)
           CALL APPDIF (ALDIF, ALPHA, XI, N, K, NCOMP, M, MSTAR)
           CALL LSYSLV (IFLAG, XI, XIOLD, XIJ, DALPHA, ALDIF, RHS,&
               ALPHA, A, IPIV, INTEGS, RNORM, 1, FSUB,&
               DFSUB, GSUB, DGSUB, SOLUTN)
           IF (RNORM .LT. PRECIS)                   GO TO 405
           IF (RNORM .LE. RNOLD)                    GO TO 120
           IF (IPRINT .LT. 0)  WRITE (IWR,510) ITER, RNORM
           IF (IPRINT .LT. 0)  WRITE (IWR,540)
           ICON = 0
           RELAX = RSTART
           DO 110 I=1,NALPHA
  110      ALPHA(I) = ALPHA(I) - DALPHA(I)
           CALL APPDIF (ALDIF, ALPHA, XI, N, K, NCOMP, M, MSTAR)
           ITER = 0
           GO TO 140
  120      IF (IFREEZ .EQ. 1)                       GO TO 130
           IFREEZ = 1
           GO TO 70
  130      IFRZ = IFRZ + 1
           IF (IFRZ .GE. LMTFRZ) IFREEZ = 0
           IF (RNOLD .LT. 4.D0*RNORM) IFREEZ = 0
           GO TO 300
  140      IF(IPRINT .LT. 0)  WRITE (IWR,500)
           CALL LSYSLV (IFLAG, XI, XIOLD, XIJ, DALPHA, ALDIF, RHS,&
               ALPHA, A, IPIV, INTEGS, RNORM, 1, FSUB,&
               DFSUB, GSUB, DGSUB, SOLUTN)
  150      RNOLD = RNORM
           IF (ITER .GE. LIMIT)                     GO TO 420
           CALL LSYSLV (IFLAG, XI, XIOLD, XIJ, DALPHA, ALDIF, RHS,&
               ALPHA, A, IPIV, INTEGS, RNORM, 2, FSUB,&
               DFSUB, GSUB, DGSUB, SOLUTN)
           IF (IFLAG .EQ. 0)                        GO TO 30
           IF (ITER .GT. 0)                         GO TO 170
           IF ( IGUESS .EQ. 1) IGUESS = 0
           NP1 = N + 1
           DO 160 I=1,NP1
  160      XIOLD(I) = XI(I)
           NOLD = N
           GO TO 190
  170      CONTINUE
           ANDIF = 0.D0
           DO 180 I=1,NALPHA
  180      ANDIF = ANDIF + (EALPHA(I) - DALPHA(I))**2&
          / (ALPHA(I)*ALPHA(I) + PRECIS)
           RELAX = RELAX * ANSCL / DMAX1( DSQRT(ANDIF/FALPHA),&
          PRECIS)
           IF (RELAX .GT. 1.D0)  RELAX = 1.D0
  190      RLXOLD = RELAX
           IPRED = 1
           ITER = ITER + 1
           DO 200 I=1,NALPHA
  200      ALPHA(I) = ALPHA(I) + RELAX * DALPHA(I)
  210      CALL APPDIF (ALDIF, ALPHA, XI, N, K, NCOMP, M, MSTAR)
           CALL LSYSLV (IFLAG, XI, XIOLD, XIJ, DALPHA, ALDIF, RHS,&
               ALPHA, A, IPIV, INTEGS, RNORM, 1, FSUB,&
               DFSUB, GSUB, DGSUB, SOLUTN)
           CALL LSYSLV (IFLAG, XI, XIOLD, XIJ, EALPHA, ALDIF, RHS,&
               ALPHA, A, IPIV, INTEGS, RNORM, 3, FSUB,&
               DFSUB, GSUB, DGSUB, SOLUTN)
           ANORM = 0.D0
           ANFIX = 0.D0
           ANSCL = 0.D0
           DO 220 I=1,NALPHA
             ANSCL = ANSCL + DALPHA(I) * DALPHA(I) /&
            (ALPHA(I)*ALPHA(I) + PRECIS)
             SCALE = ALPHA(I) - RELAX*DALPHA(I)
             SCALE = 1.D0 / (SCALE*SCALE + PRECIS)
             ANORM = ANORM + DALPHA(I) * DALPHA(I) * SCALE
  220      ANFIX = ANFIX + EALPHA(I) * EALPHA(I) * SCALE
           ANORM = DSQRT(ANORM / FALPHA)
           ANFIX = DSQRT(ANFIX / FALPHA)
           ANSCL = DSQRT(ANSCL / FALPHA)
           IF (ICOR .EQ. 1)                         GO TO 230
           IF (IPRINT .LT. 0)  WRITE (IWR,520) ITER, RELAX, ANORM,&
                ANFIX, RNOLD, RNORM
           GO TO 240
  230      IF (IPRINT .LT. 0) WRITE (IWR,550) RELAX, ANORM, ANFIX,&
                RNOLD, RNORM
  240      ICOR = 0
           IF (ANFIX.LT.PRECIS .OR. RNORM.LT.PRECIS)GO TO 405
           IF (ANFIX .GT. ANORM)                    GO TO 250
           IF (ANFIX .LE. CHECK)                    GO TO 290
           IF (IPRED .NE. 1)                        GO TO 150
  250      IF (ITER .GE. LIMIT)                     GO TO 420
           IPRED = 0
           ARG = (ANFIX/ANORM - 1.D0) / RELAX + 1.D0
           IF (ARG .LT. 0.D0)                         GO TO 150
           IF (ARG .LE. .25D0*RELAX+.125D0*RELAX**2 ) GO TO 260
           FACTOR = -1.D0 + DSQRT (1.D0+8.D0 * ARG)
           IF ( DABS(FACTOR-1.D0) .LT. .1D0*FACTOR )  GO TO 150
           RELAX = RELAX / FACTOR
           GO TO 270
  260      IF (RELAX .GE. .9D0)                       GO TO 150
           RELAX = 1.D0
  270      ICOR = 1
           IF (RELAX .LT. RELMIN)                     GO TO 430
           DO 280 I=1,NALPHA
  280      ALPHA(I) = ALPHA(I) + (RELAX-RLXOLD) * DALPHA(I)
           RLXOLD = RELAX
           GO TO 210
  290      CALL APPDIF (A, EALPHA, XI, N, K, NCOMP, M, MSTAR)
           GO TO 310
  300      CALL APPDIF (EALPHA, DALPHA, XI, N, K, NCOMP, M, MSTAR)
  310      CONTINUE
           INN = 0
           JCOL = 0
           JINIT = 1
           DO 380 I = 1, NTOL
             JEND = JTOL(I) - 1
             IF (JEND .LT. JINIT)                   GO TO 330
             DO 320 J = JINIT, JEND
               MJ = M(J)
               NALPHJ = N * K + MJ
               JCOL = JCOL + MJ
               INN = INN + MJ * NALPHJ
  320        CONTINUE
  330        JINIT = JEND + 1
             NALPHJ = N * K + M(JINIT)
             INN1 = INN
             JCOL1 = JCOL + 1
  340        IF (JCOL1 .EQ. LTOL(I))                GO TO 350
             INN1 = INN1 + NALPHJ
             JCOL1 = JCOL1 + 1
             GO TO 340
  350        IINIT = JCOL1 - JCOL
             DO 370 II = IINIT, NALPHJ
               IN = INN1 + II
               IF (ICON .EQ. 1)                     GO TO 360
               IF (DABS(A(IN)) .GT. TOLIN(I) *&
                   (DABS(ALDIF(IN)) +1.D0))        GO TO 410
               GO TO 370
  360          IF (DABS(EALPHA(IN)) .GT. TOLIN(I) *&
                   (DABS(ALDIF(IN)) + 1.D0))       GO TO 410
  370        CONTINUE
  380      CONTINUE
           IF (IPRINT .LT. 1) WRITE (IWR,560) ITER
           IF (ICON .EQ. 1)                         GO TO 450
           DO 390 I=1,NALDIF
  390      ALDIF(I) = ALDIF(I) + A(I)
           DO 400 I=1,NALPHA
  400      ALPHA(I) = ALPHA(I) + EALPHA(I)
  405      IF ((ANFIX.LT.PRECIS.OR.RNORM.LT.PRECIS).AND.IPRINT.LT.1)&
          WRITE (IWR,560) ITER
           ICON = 1
           IF (ICARE .EQ. (-1))  ICARE = 0
           GO TO 450
  410      IF ( ICON .EQ. 0)                        GO TO 150
           GO TO 70
  420      IF(IPRINT .LT. 1) WRITE (IWR,570) ITER
           GO TO 440
  430      IF(IPRINT .LT. 1)  WRITE(IWR,580) RELAX, RELMIN
  440      IFLAG = -2
           LCONV = LCONV + 1
           IF (ICARE .EQ. 2 .AND. LCONV .GT. 1)     RETURN
           IF (ICARE .EQ. 0)  ICARE = -1
           GO TO 460
  450      CALL ERRCHK(IMESH,XIOLD,ALDIF,VALSTR,A,MSTAR,IFIN)
           IF (IMESH .EQ. 1 .OR. IFIN .EQ. 0 .AND.&
               ICARE .NE. 2)                       GO TO 460
           IFLAG = 1
           RETURN
  460      IMESH = 1
           IF (ICON .EQ. 0 .OR. MSHNUM .GE. MSHLMT&
          .OR. MSHALT .GE. MSHLMT)  IMESH = 2
           IF (MSHALT .GE. MSHLMT .AND. MSHNUM .LT. MSHLMT)&
          MSHALT = 1
           CALL NEWMSH(IMESH, XI, XIOLD, XIJ, ALDIF, VALSTR,&
               SLOPE, ACCUM, NFXPNT, FIXPNT)
           IF (N .LE. NMAX)                         GO TO 470
           N = N / 2
           IFLAG = -1
           IF (ICON .EQ. 0 .AND. IPRINT .LT. 1) WRITE (IWR,590)
           IF (ICON .EQ. 1 .AND. IPRINT .LT. 1) WRITE (IWR,600)
           RETURN
  470      IF (ICON .EQ. 0)  IMESH = 1
           IF (ICARE .EQ. 1)  ICON = 0
      GO TO 20
  490 FORMAT(//24H THE MATRIX IS SINGULAR )
  500 FORMAT(/30H FULL DAMPED NEWTON ITERATION,)
  510 FORMAT(13H ITERATION = , I3, 15H  NORM (RHS) = , D10.2)
  520 FORMAT(13H ITERATION = ,I3,22H  RELAXATION FACTOR = ,D10.2&
            /33H NORM OF SCALED RHS CHANGES FROM ,D10.2,3H TO,D10.2&
            /33H NORM   OF   RHS  CHANGES  FROM  ,D10.2,3H TO,D10.2)
  530 FORMAT(/27H FIXED JACOBIAN ITERATIONS,)
  540 FORMAT(/35H SWITCH TO DAMPED NEWTON ITERATION,)
  550 FORMAT(40H RELAXATION FACTOR CORRECTED TO RELAX = , D10.2&
            /33H NORM OF SCALED RHS CHANGES FROM ,D10.2,3H TO,D10.2&
            /33H NORM   OF   RHS  CHANGES  FROM  ,D10.2,3H TO,D10.2)
  560 FORMAT(/18H CONVERGENCE AFTER , I3,11H ITERATIONS /)
  570 FORMAT(/22H NO CONVERGENCE AFTER , I3, 11H ITERATIONS/)
  580 FORMAT(/37H NO CONVERGENCE.  RELAXATION FACTOR =,D10.3&
            ,24H IS TOO SMALL (LESS THAN, D10.3, 1H)/)
  590 FORMAT(18H  (NO CONVERGENCE) )
  600 FORMAT(50H  (PROBABLY TOLERANCES TOO STRINGENT, OR NMAX TOO&
            ,6HSMALL) )
      END
      DOUBLE PRECISION FUNCTION DB2VAL(XVAL,YVAL,IDX,IDY,TX,TY,NX,NY,&
       KX,KY,BCOEF,WORK)
      INTEGER  IDX, IDY, NX, NY, KX, KY
      DOUBLE PRECISION XVAL,YVAL,TX(*), TY(*), BCOEF(NX,NY), WORK(*)
      INTEGER  ILOY, INBVX, INBV, K, LEFTY, MFLAG, KCOL, IW
      DATA ILOY /1/,  INBVX /1/
      SAVE ILOY    ,  INBVX
!$OMP THREADPRIVATE(ILOY,INBVX)
      DB2VAL = 0.0D0
      CALL DINTRV(TY,NY+KY,YVAL,ILOY,LEFTY,MFLAG)
      IF (MFLAG .NE. 0)  GO TO 100
         IW = KY + 1
         KCOL = LEFTY - KY
         DO 50 K=1,KY
            KCOL = KCOL + 1
            WORK(K) = DBVALU(TX,BCOEF(1,KCOL),NX,KX,IDX,XVAL,INBVX,&
                           WORK(IW))
   50    CONTINUE
         INBV = 1
         KCOL = LEFTY - KY + 1
         DB2VAL = DBVALU(TY(KCOL),WORK,KY,KY,IDY,YVAL,INBV,WORK(IW))
  100 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DBVALU(T,A,N,K,IDERIV,X,INBV,WORK)
      INTEGER I,IDERIV,IDERP1,IHI,IHMKMJ,ILO,IMK,IMKPJ, INBV, IPJ,&
      IP1, IP1MJ, J, JJ, J1, J2, K, KMIDER, KMJ, KM1, KPK, MFLAG, N
      DOUBLE PRECISION A, FKMJ, T, WORK, X
      DIMENSION T(*), A(N), WORK(*)
      DBVALU = 0.0D0
      IF(K.LT.1) GO TO 102
      IF(N.LT.K) GO TO 101
      IF(IDERIV.LT.0 .OR. IDERIV.GE.K) GO TO 110
      KMIDER = K - IDERIV
      KM1 = K - 1
      CALL DINTRV(T, N+1, X, INBV, I, MFLAG)
      IF (X.LT.T(K)) GO TO 120
      IF (MFLAG.EQ.0) GO TO 20
      IF (X.GT.T(I)) GO TO 130
   10 IF (I.EQ.K) GO TO 140
      I = I - 1
      IF (X.EQ.T(I)) GO TO 10
   20 IMK = I - K
      DO 30 J=1,K
        IMKPJ = IMK + J
        WORK(J) = A(IMKPJ)
   30 CONTINUE
      IF (IDERIV.EQ.0) GO TO 60
      DO 50 J=1,IDERIV
        KMJ = K - J
        FKMJ = DBLE(FLOAT(KMJ))
        DO 40 JJ=1,KMJ
          IHI = I + JJ
          IHMKMJ = IHI - KMJ
          WORK(JJ) = (WORK(JJ+1)-WORK(JJ))/(T(IHI)-T(IHMKMJ))*FKMJ
   40   CONTINUE
   50 CONTINUE
   60 IF (IDERIV.EQ.KM1) GO TO 100
      IP1 = I + 1
      KPK = K + K
      J1 = K + 1
      J2 = KPK + 1
      DO 70 J=1,KMIDER
        IPJ = I + J
        WORK(J1) = T(IPJ) - X
        IP1MJ = IP1 - J
        WORK(J2) = X - T(IP1MJ)
        J1 = J1 + 1
        J2 = J2 + 1
   70 CONTINUE
      IDERP1 = IDERIV + 1
      DO 90 J=IDERP1,KM1
        KMJ = K - J
        ILO = KMJ
        DO 80 JJ=1,KMJ
          WORK(JJ) = (WORK(JJ+1)*WORK(KPK+ILO)+WORK(JJ)&
                   *WORK(K+JJ))/(WORK(KPK+ILO)+WORK(K+JJ))
          ILO = ILO - 1
   80   CONTINUE
   90 CONTINUE
  100 DBVALU = WORK(1)
      RETURN
  101 CONTINUE
      CALL XERROR( ' DBVALU,  N DOES NOT SATISFY N.GE.K',35,2,1)
      RETURN
  102 CONTINUE
      CALL XERROR( ' DBVALU,  K DOES NOT SATISFY K.GE.1',35,2,1)
      RETURN
  110 CONTINUE
      CALL XERROR( ' DBVALU,  IDERIV DOES NOT SATISFY 0.LE.IDERIV.LT.K',&
      50, 2, 1)
      RETURN
  120 CONTINUE
      CALL XERROR( ' DBVALU,  X IS N0T GREATER THAN OR EQUAL TO T(K)',&
      48, 2, 1)
      RETURN
  130 CONTINUE
      CALL XERROR( ' DBVALU,  X IS NOT LESS THAN OR EQUAL TO T(N+1)',&
      47, 2, 1)
      RETURN
  140 CONTINUE
      CALL XERROR( ' DBVALU,  A LEFT LIMITING VALUE CANN0T BE OBTAINED A&
      T T(K)',    58, 2, 1)
      RETURN
      END
      DOUBLE PRECISION FUNCTION DEBOHR(ZP,MP,ZT,MT,ENERGY)
      implicit real*8 (a-h,o-z)
      include 'atimacnt.inc'
      INTEGER*4  ZP
      INTEGER*4  ZT
      REAL*8     MT, MP
      REAL*8     ENERGY
      REAL*8     GAMMA
      REAL*8     BSQ , GSQ
      REAL*8     BETA
      REAL*8     BOHR
      REAL*8     ZETA, ZPEFF
      real*8     titeica
      real*8    x_lisoe
      external  x_lisoe
      GAMMA = 1.D0 + ENERGY/AMU
      BSQ = 1.D0 - 1.D0/GAMMA**2
      GSQ = GAMMA*GAMMA
      BETA = DSQRT(BSQ)
      ZETA = 1.d0 - exp( -130.d0 * BETA / (ZP ** (2.d0/3.d0)) )
      ZPEFF = ZETA * DBLE(ZP)
      titeica = 24.89 * ZT**1.2324 / (EMASS*1.d6 * bsq) *&
               DLOG( 2.d0 * EMASS*1.d6 * bsq/(33.05*ZT**1.6364) )
      titeica = MAX( titeica, 0.d0 )
      BOHR =  4.d-8*PI*ECHARGE2**2 * N_AVO  * ZPEFF**2 * DBLE(ZT)/MT
      DEBOHR = BOHR * ( x_lisoe(ZP, MP, BETA) * GSQ + titeica )
      RETURN
      END
      DOUBLE PRECISION FUNCTION DEDX&
       (ZP,ZT,MP,MT,POT,RHO,FNTP,ENERGY,GAS)
      implicit real*8 (a-h,o-z)
      INTEGER*4 ZP,ZT,GAS
      REAL*8 MP,MT,POT,RHO,FNTP,ENERGY
      REAL*8 DDEDX,SEOUT,ENERGS,ZETA
      REAL*8 BETHE,DEDXELA,BETHEK
      real*8 factor
      common /paraspar/ mat_spar    !S.Abe 2017/09/13
      common /fixcharge/ifixchg     !T.Sato 2019/02/17
      ENERGS = ENERGY*1000.D0*MP
      IF (ENERGY .LE. 10.d0) THEN
         if( ZT .le. 92 ) then
            CALL SEZI(ZP,ZT,MP,MT,ENERGS,RHO,SEOUT,POT,ZETA)
         else
            if( ZP .eq. 1 ) then
                if( MP .ge. 2.0d0 ) then
                  ipty = 1
                elseif( MP .ge. 0.20d0 .and. MP .lt. 2.0d0 ) then
                  ipty = 2
                elseif( MP .ge. 0.12d0 .and. MP .lt. 0.2d0 ) then
                  ipty = 3
                elseif( MP .lt. 0.12d0 ) then
                  ipty = 4
                endif
            else
               ipty = 1
            endif
            EN0 = ENERGY*MP
            CALL DEDXspar_elec(ipty,MP,dble(ZP),EN0,mat_spar,SEOUT)
            SEOUT = SEOUT / 1.d3 / RHO
         endif
         DDEDX = SEOUT    + DEDXELA(ZP,ZT,MP,MT,ENERGS)
      ELSE IF (ENERGY .LE. 30.d0 .or. ifixchg.ge.1) THEN ! T.Sato, for fixed charge mode, always use SEZI
         factor = DMIN1(1.0d0,0.05d0 * ( energy - 10.d0 )) ! maximum factor is 1
         if( ZT .le. 92 ) then
            CALL SEZI(ZP,ZT,MP,MT,ENERGS,RHO,SEOUT,POT,ZETA)
         else
            if( ZP .eq. 1 ) then
                if( MP .ge. 2.0d0 ) then
                  ipty = 1
                elseif( MP .ge. 0.20d0 .and. MP .lt. 2.0d0 ) then
                  ipty = 2
                elseif( MP .ge. 0.12d0 .and. MP .lt. 0.2d0 ) then
                  ipty = 3
                elseif( MP .lt. 0.12d0 ) then
                  ipty = 4
                endif
            else
               ipty = 1
            endif
            EN0 = ENERGY*MP
            CALL DEDXspar_elec(ipty,MP,dble(ZP),EN0,mat_spar,SEOUT)
            SEOUT = SEOUT / 1.d3 / RHO
         endif
         ddedx = (1.d0 - factor) * seout +&
            factor * BETHEK(ZP,ZT,MP,MT,ENERGS,RHO,GAS,FNTP,POT)&
             + DEDXELA(ZP,ZT,MP,MT,ENERGS)
      ELSE
         DDEDX = BETHEK(ZP,ZT,MP,MT,ENERGS,RHO,GAS,FNTP,POT)&
        + DEDXELA(ZP,ZT,MP,MT,ENERGS)
      ENDIF
      DEDX = DDEDX
      RETURN
      END
      SUBROUTINE DEDXspar_elec(ITYP,XM,Z,E,MED,STP)
      use MEMBANKMOD !FURUTA
      use moddas_material
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'
      common /kmat1a/ mxmat, mxmat0, mxnel
      common /sparcmn/ khb, kh1, kh2, kh3, khe(kvlmax), khp,&
                      kdf, kdg, kro
      common /kmat1g/ kmat(kvlmax)
      real(8),parameter:: Q2D3=2.0d0/3.0d0
      real(8),parameter:: Q1D3=1.0d0/3.0d0
      DATA MOLD/0/,ZOLD/0./,XMOLD/0./,IOLD/0/
      DIMENSION XMF(4)
      DATA XMF/0.,1.,0.1488,0.1129/
      save MOLD,ZOLD,XMOLD,IOLD !FURUTA
!$OMP THREADPRIVATE(MOLD,ZOLD,XMOLD,IOLD)
       save XKI             !FURUTA
!$OMP THREADPRIVATE(XKI)
       DATA XKI/0./         !FURUTA
      DATA NFIRST/0/
      save NFIRST !FURUTA
!$OMP THREADPRIVATE(NFIRST)
! T.Sato 2016/2/17 Shielding distribution calculation mode
      common /paraj/  mstz(300), parz(300)
      ismode=mstz(116)
      ROsum=0.0d0
       MOLD = 0
       ZOLD = 0.
       XMOLD = 0.
       IOLD = 0
      if( nfirst .eq. 0 ) then
         ErrCha = ''
         ErrID = 'L:3280/R:DEDXspar_elec/F:atima01.f' !W04_001_001
         call ErrWrite(ErrID,ErrCha)
         write(*,'(2a)') "*** Warning: No original data found in ATIMA",&
                        " for heavy target nuclei (Z>92)."
         write(*,'(2a)') "    Extended parameters were used for those ",&
                        "calculations (See PHITS manual)."
         write(*,'(a)') " "
         call prepspar
         nfirst = nfirst + 1
      end if
     
      IF(E.EQ.0.) GO TO 1001
      IF(ITYP.EQ.1) GO TO 120
      ZI = 1.0
      XMI = XMF(ITYP)
      XMAS = 938.232 * XMI
      BCUT1 = .07
      BCUT2 = .0046
      GO TO 130
  120 ZI = Z
      XMI = XM
      XMAS = 931.141 * XMI
      BCUT1 = .07 * Z**Q2D3
      BCUT2 = .0046 * Z**Q1D3
  130 ESTAR = E / XMAS
      TEMP1 = (ESTAR + 1.D0)**2
      ESQ = TEMP1 - 1.D0
      BSQ = ESQ / TEMP1
      B = DSQRT(BSQ)
      IF(B.LT.BCUT1) GO TO 135
      ZSQ = ZI * ZI
      SNUC = 0.
      IFROM = 1
      if(ismode.ne.1) GO TO 150 ! T.Sato 2016/2/17, should not goto 150 for shielding calculation mode
  135 CONTINUE
      ISAME = 1
      IF(ITYP.EQ.IOLD .AND.&
        XM.EQ.XMOLD.AND.Z.EQ.ZOLD.AND.MED.EQ.MOLD) GO TO 144
  137 XMOLD = XM
      ZOLD = Z
  138 IOLD = ITYP
      MOLD = MED
      ISAME = 0
      NELM = nint( dnel_das(kmat0+med) )
      ZIQ2D3 = ZI**Q2D3
      DO 140 J = 1,NELM
      Q4 = a_das(kmat(med)+j)
      Q1 = SQRT(ZIQ2D3 + zz_das(kmat(med)+j)**Q2D3)
      Q2 = Q4 + XMI
      Q3 = ZI*zz_das(kmat(med)+j)
      F(kdf+J) = 3.255D4 * Q4 / (Q2 * Q3 * Q1)
      G(kdg+J) = 1.96D-4 * Q4 * Q2 * Q1 / (Q3 * XMI)
      RO(kro+J) = Q4*den_das(kmat(med)+j)/0.602
      ROsum=ROsum+RO(kro+J) ! T.Sato 2016/2/17
  140 CONTINUE
      IF( denh_das(kmat0+med) .EQ.0. ) GO TO 142
      NELM = NELM + 1
      J = NELM
      Q2 = 1.0 + XMI
      Q1 = SQRT(ZIQ2D3 + 1.0)
      F(kdf+J) = 3.255D4 / (Q2 * ZI * Q1)
      G(kdg+J) = 1.96D-4 * Q2 * Q1 / (ZI * XMI)
      RO(kro+J) = denh_das(kmat0+med)/0.602
      ROsum=ROsum+RO(kro+J) ! T.Sato 2016/2/17
  142 CONTINUE
      if(ismode.eq.1) then ! density calculation mode, T.Sato 2016/2/17
       STP=ROsum
       return
      endif
  144 CONTINUE
      IF(B.LT.BCUT2) GO TO 200
      SNUC = 0.
      IFROM = 1
  145 ZSTAR = ZI * (1.0 - EXP(-125.0*B/ZIQ2D3))
      ZSQ = ZSTAR * ZSTAR
  150 A6 = DLOG(ESQ) - BSQ + 0.0217615D0
      DEL = DLOG(1.378D-9 * S1(kh1+MED) * ESQ)&
         - 2. * S2(kh2+MED)/S1(kh1+MED) - 1.
      IF(DEL.LT.0.) DEL = 0.
      CALL SHELLspar(BSQ,S1(kh1+MED),S2(kh2+MED),S3(kh3+MED),COZ)
      A7 = 5.0985D-1 * ((A6-DEL/2.-COZ)*S1(kh1+MED)- S2(kh2+MED))/ BSQ
      IF(IFROM.EQ.2) GO TO 215
      STP = ZSQ * A7 + SNU      GO TO 999
  200 CONTINUE
      IF(ISAME.EQ.1) GO TO 220
      IF(ITYP.EQ.1) GO TO 205
      B = 0.0046
      GO TO 210
  205 B = 0.0046 * Z ** Q1D3
  210 BSQ = B * B
      ESQ = BSQ / (1. - BSQ)
      EI = (DSQRT(1./(1.-BSQ)) -1.) * XMAS
      IFROM = 2
      GO TO 145
  215 SI = A7 * ZSQ
      XKI = SI / SQRT(EI)
  220 CONTINUE
      SE = XKI * SQRT(E)
      SNUC = 0.
      STP = SE + SNU
  999 RETURN
 1001 STP = 0.
      RETURN
      END
      DOUBLE PRECISION  FUNCTION DEDXELA(ZP,ZT,MP,MT,ENERGY)
      implicit real*8 (a-h,o-z)
      include 'atimacnt.inc'
      INTEGER*4    ZT,ZP
      REAL*8       MP,MT,ENERGY,EPSIL,A,SN
      EPSIL = 32.53D0*MT*ENERGY&
        /(DBLE(ZP)*DBLE(ZT)*(MP+MT)*(DBLE(ZP)**.23D0+DBLE(ZT)**.23D0))
      IF (EPSIL .GE. 30.0D0) THEN
        SN = DLOG(EPSIL)/(2.0D0*EPSIL)
      ELSE
        A = (.01321D0*(EPSIL**.21226D0))+(.19593D0*(EPSIL**.5D0))
        SN = .5D0*DLOG(1.0D0+1.1383D0*EPSIL)/(EPSIL+A)
      END IF
      SN = SN*DBLE(ZP)*DBLE(ZT)*MP*8.462D0/&
           ((MP+MT)*(DBLE(ZP)**.23D0+DBLE(ZT)**.23D0))
      SN = SN * .1*N_AVO / MT
      DEDXELA = SN
      RETURN
      END
      double precision function dedxver(energy)
      use dedx_file
      implicit real*8 (a-h,o-z)
      double precision energy
      include 'atimasys.inc'
      double precision m, sum
      integer i
      double precision dedx
      integer kfno
      m = 0.d0
      do i = 1, nnu        m = m + anuc(i) * mt(i)
      enddo
      kfno = 0
      if( dedx_filename(dedxmat) /= ' ' ) then
       call get_dedx_kfcode(dedxmat,kfno,zp,mp)
      endif
      if( kfno /= 0 ) then
       dedxver = 0.d0
       call dedx_file_select(m,kfno,energy,dedxver)
      else
       sum = 0.d0
       dedxver = 0.d0
       do i = 1, nnu       sum = sum+anuc(i) * mt(i) *&
         dedx(zp, zt(i), mp, mt(i), pot(i), rho, fntp, energy, gas) / m
       enddo
       dedxver = sum
      endif
      return
      end
      DOUBLE PRECISION FUNCTION  DEFIRSO (Z1,Z2,M2,ENERGY)
      implicit real*8 (a-h,o-z)
      include 'atimacnt.inc'
      INTEGER*4 Z1
      INTEGER*4 Z2
      REAL*8    M2
      REAL*8    ENERGY
      REAL*8    GAMMA
      REAL*8    BSQ
      REAL*8    BETA
      REAL*8    FIRSOV
      GAMMA = 1.D0+ENERGY/AMU
      BSQ = 1.D0-1.D0/GAMMA**2
      BETA = SQRT(BSQ)
      FIRSOV = 4.8184D-6*DBLE(Z1+Z2)**(8.D0/3.D0)/M2
      DEFIRSO = FIRSOV * ( BETA**2  / ALPHA**2 )
      RETURN
      END
      SUBROUTINE DFF(X,Z,DF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Z(1),DF(1,1)
      DF(1,1)=0.D0
      RETURN
      END
      SUBROUTINE DGF(I,Z,DG)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Z(1),DG(1,1)
      DG(1,1)=1.D0
      RETURN
      END
      SUBROUTINE DICHTEEFF(d_beta,d_z2,d_a2,d_pot,d_rho,
     $                     d_del_84)
      implicit real*8 (a-h,o-z)
      double precision d_del_71,d_del_84
      double precision d_beta
      integer*4 d_z2
      double precision d_a2
      double precision d_pot
      double precision d_rho
      double precision d_gamma
      double precision d_x
      double precision  d_del_0 (97)
      double precision  d_m (97)
      double precision  d_x1 (97)
      double precision  d_x0 (97)
      double precision  d_a (97)
      double precision  d_c (97)
      d_gamma = 1. / sqrt(1. - d_beta**2)
      d_x = log(d_beta * d_gamma) / 2.302585
  	d_x0(1)     = 1.8639
 	d_x0(2)     = 2.2017
      	d_x0(3)     = 0.1304
      	d_x0(4)     = 0.0592
        d_x0(5)     = 0.0305
        d_x0(6)     = -.0178
        d_x0(7)     = 1.7378
        d_x0(8)     = 1.7541
        d_x0(9)     = 1.8433
        d_x0(10)    = 2.0735
        d_x0(11)    = 0.2880
        d_x0(12)    = 0.1499
        d_x0(13)    = 0.1708
        d_x0(14)    = 0.2014
        d_x0(15)    = 0.1696
        d_x0(16)    = 0.1580
        d_x0(17)    = 1.5555
        d_x0(18)    = 1.7635
        d_x0(19)    = 0.3851
        d_x0(20)    = 0.3228
        d_x0(21)    = 0.1640
        d_x0(22)    = 0.0957
        d_x0(23)    = 0.0691
        d_x0(24)    = 0.0340
        d_x0(25)    = 0.0447
        d_x0(26)    = -.0012
        d_x0(27)    = -.0187
        d_x0(28)    = -.0566
        d_x0(29)    = -.0254
        d_x0(30)    = 0.0049
        d_x0(31)    = 0.2267
        d_x0(32)    = 0.3376
        d_x0(33)    = 0.1767
        d_x0(34)    = 0.2258
        d_x0(35)    = 1.5262
        d_x0(36)    = 1.7158
        d_x0(37)    = 0.5737
        d_x0(38)    = 0.4585
        d_x0(39)    = 0.3608
        d_x0(40)    = 0.2957
        d_x0(41)    = 0.1785
        d_x0(42)    = 0.2267
        d_x0(43)    = 0.0949
        d_x0(44)    = 0.0599
        d_x0(45)    = 0.0576
        d_x0(46)    = 0.0563
        d_x0(47)    = 0.0657
        d_x0(48)    = 0.1281
        d_x0(49)    = 0.2406
        d_x0(50)    = 0.2879
        d_x0(51)    = 0.3189
        d_x0(52)    = 0.3296
        d_x0(53)    = 0.0549
        d_x0(54)    = 1.5630
        d_x0(55)    = 0.5473
        d_x0(56)    = 0.4190
        d_x0(57)    = 0.3161
        d_x0(58)    = 0.2713
        d_x0(59)    = 0.2333
        d_x0(60)    = 0.1984
        d_x0(61)    = 0.1627
        d_x0(62)    = 0.1520
        d_x0(63)    = 0.1888
        d_x0(64)    = 0.1058
        d_x0(65)    = 0.0947
        d_x0(66)    = 0.0822
        d_x0(67)    = 0.0761
        d_x0(68)    = 0.0648
        d_x0(69)    = 0.0812
        d_x0(70)    = 0.1199
        d_x0(71)    = 0.1560
        d_x0(72)    = 0.1965
        d_x0(73)    = 0.2117
        d_x0(74)    = 0.2167
        d_x0(75)    = 0.0559
        d_x0(76)    = 0.0891
        d_x0(77)    = 0.0819
        d_x0(78)    = 0.1484
        d_x0(79)    = 0.2021
        d_x0(80)    = 0.2756
        d_x0(81)    = 0.3491
        d_x0(82)    = 0.3776
        d_x0(83)    = 0.4152
        d_x0(84)    = 0.4267
        d_x0(86)    = 1.5368
        d_x0(88)    = 0.5991
        d_x0(89)    = 0.4559
        d_x0(90)    = 0.4202
        d_x0(91)    = 0.3144
        d_x0(92)    = 0.2260
        d_x0(93)    = 0.1869
        d_x0(94)    = 0.1557
        d_x0(95)    = 0.2274
        d_x0(96)    = 0.2484
        d_x0(97)    = 0.2378
        d_x1(2)     = 3.6122
        d_x1(3)     = 1.6397
        d_x1(4)     = 1.6922
        d_x1(5)     = 1.9688
        d_x1(6)     = 2.3415
        d_x1(7)     = 4.1323
        d_x1(8)     = 4.3213
        d_x1(9)     = 4.4096
        d_x1(10)    = 4.6421
        d_x1(11)    = 3.1962
        d_x1(12)    = 3.0668
        d_x1(13)    = 3.0127
        d_x1(14)    = 2.8715
        d_x1(15)    = 2.7815
        d_x1(16)    = 2.7159
        d_x1(17)    = 4.2994
        d_x1(18)    = 4.4855
        d_x1(19)    = 3.1724
        d_x1(20)    = 3.1191
        d_x1(21)    = 3.0593
        d_x1(22)    = 3.0386
        d_x1(23)    = 3.0322
        d_x1(24)    = 3.0451
        d_x1(25)    = 3.1074
        d_x1(26)    = 3.1531
        d_x1(27)    = 3.1790
        d_x1(28)    = 3.1851
        d_x1(29)    = 3.2792
        d_x1(30)    = 3.3668
        d_x1(31)    = 3.5434
        d_x1(32)    = 3.6096
        d_x1(33)    = 3.5702
        d_x1(34)    = 3.6264
        d_x1(35)    = 4.9899
        d_x1(36)    = 5.0748
        d_x1(37)    = 3.7995
        d_x1(38)    = 3.6778
        d_x1(39)    = 3.5542
        d_x1(40)    = 3.4890
        d_x1(41)    = 3.2201
        d_x1(42)    = 3.2784
        d_x1(43)    = 3.1253
        d_x1(44)    = 3.0834
        d_x1(45)    = 3.1069
        d_x1(46)    = 3.0555
        d_x1(47)    = 3.1074
        d_x1(48)    = 3.1667
        d_x1(49)    = 3.2032
        d_x1(50)    = 3.2959
        d_x1(51)    = 3.3489
        d_x1(52)    = 3.4418
        d_x1(53)    = 3.2596
        d_x1(54)    = 4.7371
        d_x1(55)    = 3.5914
        d_x1(56)    = 3.4547
        d_x1(57)    = 3.3293
        d_x1(58)    = 3.3432
        d_x1(59)    = 3.2773
        d_x1(60)    = 3.3063
        d_x1(61)    = 3.3199
        d_x1(62)    = 3.3460
        d_x1(63)    = 3.4633
        d_x1(64)    = 3.3932
        d_x1(65)    = 3.4224
        d_x1(66)    = 3.4474
        d_x1(67)    = 3.4782
        d_x1(68)    = 3.4922
        d_x1(69)    = 3.5085
        d_x1(70)    = 3.6246
        d_x1(71)    = 3.5218
        d_x1(72)    = 3.4337
        d_x1(73)    = 3.4805
        d_x1(74)    = 3.4960
        d_x1(75)    = 3.4845
        d_x1(76)    = 3.5414
        d_x1(77)    = 3.5480
        d_x1(78)    = 3.6212
        d_x1(79)    = 3.6979
        d_x1(80)    = 3.7275
        d_x1(81)    = 3.8044
        d_x1(82)    = 3.8073
        d_x1(83)    = 3.8248
        d_x1(84)    = 3.8293
        d_x1(86)    = 4.9889
        d_x1(88)    = 3.9428
        d_x1(89)    = 3.7966
        d_x1(90)    = 3.7681
        d_x1(91)    = 3.5079
        d_x1(92)    = 3.3721
        d_x1(93)    = 3.3690
        d_x1(94)    = 3.3981
        d_x1(95)    = 3.5021
        d_x1(96)    = 3.5160
        d_x1(97)    = 3.5186
        d_a(2)      = 0.13443
        d_a(3)      = 0.95136
        d_a(4)      = 0.80392
        d_a(5)      = 0.56224
        d_a(6)      = 0.26142
        d_a(7)      = 0.15349
        d_a(8)      = 0.11778
        d_a(9)      = 0.11083
        d_a(10)     = 0.08064
        d_a(11)     = 0.07772
        d_a(12)     = 0.08163
        d_a(13)     = 0.08024
        d_a(14)     = 0.14921
        d_a(15)     = 0.23610
        d_a(16)     = 0.33992
        d_a(17)     = 0.19849
        d_a(18)     = 0.19714
        d_a(19)     = 0.19827
        d_a(20)     = 0.15643
        d_a(21)     = 0.15754
        d_a(22)     = 0.15662
        d_a(23)     = 0.15436
        d_a(24)     = 0.15419
        d_a(25)     = 0.14973
        d_a(26)     = 0.14680
        d_a(27)     = 0.14474
        d_a(28)     = 0.16496
        d_a(29)     = 0.14339
        d_a(30)     = 0.14714
        d_a(31)     = 0.09440
        d_a(32)     = 0.07188
        d_a(33)     = 0.06633
        d_a(34)     = 0.06568
        d_a(35)     = 0.06335
        d_a(36)     = 0.07446
        d_a(37)     = 0.07261
        d_a(38)     = 0.07165
        d_a(39)     = 0.07138
        d_a(40)     = 0.07177
        d_a(41)     = 0.13883
        d_a(42)     = 0.10525
        d_a(43)     = 0.16572
        d_a(44)     = 0.19342
        d_a(45)     = 0.19205
        d_a(46)     = 0.24178
        d_a(47)     = 0.24585
        d_a(48)     = 0.24609
        d_a(49)     = 0.23879
        d_a(50)     = 0.18689
        d_a(51)     = 0.16652
        d_a(52)     = 0.13815
        d_a(53)     = 0.23766
        d_a(54)     = 0.23314
        d_a(55)     = 0.18233
        d_a(56)     = 0.18268
        d_a(57)     = 0.18591
        d_a(58)     = 0.18885
        d_a(59)     = 0.23265
        d_a(60)     = 0.23530
        d_a(61)     = 0.24280
        d_a(62)     = 0.24698
        d_a(63)     = 0.24448
        d_a(64)     = 0.25109
        d_a(65)     = 0.24453
        d_a(66)     = 0.24665
        d_a(67)     = 0.24638
        d_a(68)     = 0.24823
        d_a(69)     = 0.24889
        d_a(70)     = 0.25295
        d_a(71)     = 0.24033
        d_a(72)     = 0.22918
        d_a(73)     = 0.17798
        d_a(74)     = 0.15509
        d_a(75)     = 0.15184
        d_a(76)     = 0.12751
        d_a(77)     = 0.12690
        d_a(78)     = 0.11128
        d_a(79)     = 0.09756
        d_a(80)     = 0.11014
        d_a(81)     = 0.09455
        d_a(82)     = 0.09359
        d_a(83)     = 0.09410
        d_a(84)     = 0.09282
        d_a(86)     = 0.20798
        d_a(88)     = 0.08804
        d_a(89)     = 0.08567
        d_a(90)     = 0.08655
        d_a(91)     = 0.14770
        d_a(92)     = 0.19677
        d_a(93)     = 0.19741
        d_a(94)     = 0.20419
        d_a(95)     = 0.20308
        d_a(96)     = 0.20257
        d_a(97)     = 0.20192
        d_c(2)      = 11.1393
        d_c(3)      = 3.1221
        d_c(4)      = 2.7847
        d_c(5)      = 2.8477
        d_c(6)      = 2.8680
        d_c(7)      = 10.5400
        d_c(8)      = 10.7004
        d_c(9)      = 10.9653
        d_c(10)     = 11.9041
        d_c(11)     = 5.0526
        d_c(12)     = 4.5297
        d_c(13)     = 4.2395
        d_c(14)     = 4.4351
        d_c(15)     = 4.5214
        d_c(16)     = 4.6659
        d_c(17)     = 11.1421
        d_c(18)     = 11.9480
        d_c(19)     = 5.6423
        d_c(20)     = 5.0396
        d_c(21)     = 4.6949
        d_c(22)     = 4.4450
        d_c(23)     = 4.2659
        d_c(24)     = 4.1781
        d_c(25)     = 4.2702
        d_c(26)     = 4.2911
        d_c(27)     = 4.2601
        d_c(28)     = 4.3115
        d_c(29)     = 4.4190
        d_c(30)     = 4.6906
        d_c(31)     = 4.9353
        d_c(32)     = 5.1411
        d_c(33)     = 5.0510
        d_c(34)     = 5.3210
        d_c(35)     = 11.7307
        d_c(36)     = 12.5115
        d_c(37)     = 6.4776
        d_c(38)     = 5.9867
        d_c(39)     = 5.4801
        d_c(40)     = 5.1774
        d_c(41)     = 5.0141
        d_c(42)     = 4.8793
        d_c(43)     = 4.7769
        d_c(44)     = 4.7694
        d_c(45)     = 4.8008
        d_c(46)     = 4.9358
        d_c(47)     = 5.0630
        d_c(48)     = 5.2727
        d_c(49)     = 5.5211
        d_c(50)     = 5.5340
        d_c(51)     = 5.6241
        d_c(52)     = 5.7131
        d_c(53)     = 5.9488
        d_c(54)     = 12.7281
        d_c(55)     = 6.9135
        d_c(56)     = 6.3153
        d_c(57)     = 5.7850
        d_c(58)     = 5.7837
        d_c(59)     = 5.8096
        d_c(60)     = 5.8290
        d_c(61)     = 5.8224
        d_c(62)     = 5.8597
        d_c(63)     = 6.2278
        d_c(64)     = 5.8738
        d_c(65)     = 5.9045
        d_c(66)     = 5.9183
        d_c(67)     = 5.9587
        d_c(68)     = 5.9521
        d_c(69)     = 5.9677
        d_c(70)     = 6.3325
        d_c(71)     = 5.9785
        d_c(72)     = 5.7139
        d_c(73)     = 5.5262
        d_c(74)     = 5.4059
        d_c(75)     = 5.3445
        d_c(76)     = 5.3083
        d_c(77)     = 5.3418
        d_c(78)     = 5.4732
        d_c(79)     = 5.5747
        d_c(80)     = 5.9605
        d_c(81)     = 6.1365
        d_c(82)     = 6.2018
        d_c(83)     = 6.3505
        d_c(84)     = 6.4003
        d_c(86)     = 13.2839
        d_c(88)     = 7.0452
        d_c(89)     = 6.3742
        d_c(90)     = 6.2473
        d_c(91)     = 6.0327
        d_c(92)     = 5.8694
        d_c(93)     = 5.8149
        d_c(94)     = 5.8748
        d_c(95)     = 6.2813
        d_c(96)     = 6.3097
        d_c(97)     = 6.2912
        d_m(2)      = 5.8347
        d_m(3)      = 2.4993
        d_m(4)      = 2.4339
        d_m(5)      = 2.4512
        d_m(6)      = 2.8697
        d_m(7)      = 3.2125
        d_m(8)      = 3.2913
        d_m(9)      = 3.2962
        d_m(10)     = 3.5771
        d_m(11)     = 3.6452
        d_m(12)     = 3.6166
        d_m(13)     = 3.6345
        d_m(14)     = 3.2546
        d_m(15)     = 2.9158
        d_m(16)     = 2.6456
        d_m(17)     = 2.9702
        d_m(18)     = 2.9618
        d_m(19)     = 2.9233
        d_m(20)     = 3.0745
        d_m(21)     = 3.0517
        d_m(22)     = 3.0302
        d_m(23)     = 3.0163
        d_m(24)     = 2.9896
        d_m(25)     = 2.9796
        d_m(26)     = 2.9632
        d_m(27)     = 2.9502
        d_m(28)     = 2.8430
        d_m(29)     = 2.9044
        d_m(30)     = 2.8652
        d_m(31)     = 3.1314
        d_m(32)     = 3.3306
        d_m(33)     = 3.4176
        d_m(34)     = 3.4317
        d_m(35)     = 3.4670
        d_m(36)     = 3.4051
        d_m(37)     = 3.4177
        d_m(38)     = 3.4435
        d_m(39)     = 3.4585
        d_m(40)     = 3.4533
        d_m(41)     = 3.0930
        d_m(42)     = 3.2549
        d_m(43)     = 2.9738
        d_m(44)     = 2.8707
        d_m(45)     = 2.8633
        d_m(46)     = 2.7239
        d_m(47)     = 2.6899
        d_m(48)     = 2.6772
        d_m(49)     = 2.7144
        d_m(50)     = 2.8576
        d_m(51)     = 2.9319
        d_m(52)     = 3.0354
        d_m(53)     = 2.7276
        d_m(54)     = 2.7414
        d_m(55)     = 2.8866
        d_m(56)     = 2.8906
        d_m(57)     = 2.8828
        d_m(58)     = 2.8592
        d_m(59)     = 2.7331
        d_m(60)     = 2.7050
        d_m(61)     = 2.6674
        d_m(62)     = 2.6403
        d_m(63)     = 2.6245
        d_m(64)     = 2.5977
        d_m(65)     = 2.6056
        d_m(66)     = 2.5849
        d_m(67)     = 2.5726
        d_m(68)     = 2.5573
        d_m(69)     = 2.5469
        d_m(70)     = 2.5141
        d_m(71)     = 2.5643
        d_m(72)     = 2.6155
        d_m(73)     = 2.7623
        d_m(74)     = 2.8447
        d_m(75)     = 2.8627
        d_m(76)     = 2.9608
        d_m(77)     = 2.9658
        d_m(78)     = 3.0417
        d_m(79)     = 3.1101
        d_m(80)     = 3.0519
        d_m(81)     = 3.1450
        d_m(82)     = 3.1608
        d_m(83)     = 3.1671
        d_m(84)     = 3.1830
        d_m(86)     = 2.7409
        d_m(88)     = 3.2454
        d_m(89)     = 3.2683
        d_m(90)     = 3.2610
        d_m(91)     = 2.9845
        d_m(92)     = 2.8171
        d_m(93)     = 2.8082
        d_m(94)     = 2.7679
        d_m(95)     = 2.7615
        d_m(96)     = 2.7579
        d_m(97)     = 2.7560
        d_del_0(2)  = 0.0
        d_del_0(3)  = 0.14
        d_del_0(4)  = 0.14
        d_del_0(5)  = 0.14
        d_del_0(6)  = 0.12
        d_del_0(7)  = 0.0
        d_del_0(8)  = 0.0
        d_del_0(9)  = 0.0
        d_del_0(10) = 0.0
        d_del_0(11) = 0.08
        d_del_0(12) = 0.08
        d_del_0(13) = 0.12
        d_del_0(14) = 0.14
        d_del_0(15) = 0.14
        d_del_0(16) = 0.14
        d_del_0(17) = 0.0
        d_del_0(18) = 0.0
        d_del_0(19) = 0.10
        d_del_0(20) = 0.14
        d_del_0(21) = 0.10
        d_del_0(22) = 0.12
        d_del_0(23) = 0.14
        d_del_0(24) = 0.14
        d_del_0(25) = 0.14
        d_del_0(26) = 0.12
        d_del_0(27) = 0.12
        d_del_0(28) = 0.10
        d_del_0(29) = 0.08
        d_del_0(30) = 0.08
        d_del_0(31) = 0.14
        d_del_0(32) = 0.14
        d_del_0(33) = 0.08
        d_del_0(34) = 0.10
        d_del_0(35) = 0.0
        d_del_0(36) = 0.0
        d_del_0(37) = 0.14
        d_del_0(38) = 0.14
        d_del_0(39) = 0.14
        d_del_0(40) = 0.14
        d_del_0(41) = 0.14
        d_del_0(42) = 0.14
        d_del_0(43) = 0.14
        d_del_0(44) = 0.14
        d_del_0(45) = 0.14
        d_del_0(46) = 0.14
        d_del_0(47) = 0.14
        d_del_0(48) = 0.14
        d_del_0(49) = 0.14
        d_del_0(50) = 0.14
        d_del_0(51) = 0.14
        d_del_0(52) = 0.14
        d_del_0(53) = 0.0
        d_del_0(54) = 0.0
        d_del_0(55) = 0.14
        d_del_0(56) = 0.14
        d_del_0(57) = 0.14
        d_del_0(58) = 0.14
        d_del_0(59) = 0.14
        d_del_0(60) = 0.14
        d_del_0(61) = 0.14
        d_del_0(62) = 0.14
        d_del_0(63) = 0.14
        d_del_0(64) = 0.14
        d_del_0(65) = 0.14
        d_del_0(66) = 0.14
        d_del_0(67) = 0.14
        d_del_0(68) = 0.14
        d_del_0(69) = 0.14
        d_del_0(70) = 0.14
        d_del_0(71) = 0.14
        d_del_0(72) = 0.14
        d_del_0(73) = 0.14
        d_del_0(74) = 0.14
        d_del_0(75) = 0.08
        d_del_0(76) = 0.10
        d_del_0(77) = 0.10
        d_del_0(78) = 0.12
        d_del_0(79) = 0.14
        d_del_0(80) = 0.14
        d_del_0(81) = 0.14
        d_del_0(82) = 0.14
        d_del_0(83) = 0.14
        d_del_0(84) = 0.14
        d_del_0(86) = 0.0
        d_del_0(88) = 0.14
        d_del_0(89) = 0.14
        d_del_0(90) = 0.14
        d_del_0(91) = 0.14
        d_del_0(92) = 0.14
        d_del_0(93) = 0.14
        d_del_0(94) = 0.14
        d_del_0(95) = 0.14
        d_del_0(96) = 0.14
        d_del_0(97) = 0.14
      d_c(d_z2) = - d_c(d_z2)
      if (d_x .lt. d_x0(d_z2) ) then
        if (d_del_0(d_z2) .gt. 0.) then
          d_del_84 = d_del_0(d_z2) * 1.d1**( 2.d0*(d_x-d_x0(d_z2)) )
        elseif (d_del_0(d_z2) .eq. 0.) then
          d_del_84 = 0.
        endif
      elseif ( (d_x0(d_z2) .le. d_x) .and. (d_x .le. d_x1(d_z2)) )then
        d_del_84 = 4.6052 * d_x +
     &             d_a(d_z2) * (d_x1(d_z2) - d_x)**d_m(d_z2) + d_c(d_z2)
      elseif (d_x .gt. d_x1(d_z2) ) then
        d_del_84 = 4.6052 * d_x + d_c(d_z2)
      endif
      return
      end
      SUBROUTINE DINTRV(XT,LXT,X,ILO,ILEFT,MFLAG)
      INTEGER IHI, ILEFT, ILO, ISTEP, LXT, MFLAG, MIDDLE
      DOUBLE PRECISION X, XT
      DIMENSION XT(LXT)
      IHI = ILO + 1
      IF (IHI.LT.LXT) GO TO 10
      IF (X.GE.XT(LXT)) GO TO 110
      IF (LXT.LE.1) GO TO 90
      ILO = LXT - 1
      IHI = LXT
   10 IF (X.GE.XT(IHI)) GO TO 40
      IF (X.GE.XT(ILO)) GO TO 100
      ISTEP = 1
   20 IHI = ILO
      ILO = IHI - ISTEP
      IF (ILO.LE.1) GO TO 30
      IF (X.GE.XT(ILO)) GO TO 70
      ISTEP = ISTEP*2
      GO TO 20
   30 ILO = 1
      IF (X.LT.XT(1)) GO TO 90
      GO TO 70
   40 ISTEP = 1
   50 ILO = IHI
      IHI = ILO + ISTEP
      IF (IHI.GE.LXT) GO TO 60
      IF (X.LT.XT(IHI)) GO TO 70
      ISTEP = ISTEP*2
      GO TO 50
   60 IF (X.GE.XT(LXT)) GO TO 110
      IHI = LXT
   70 MIDDLE = (ILO+IHI)/2
      IF (MIDDLE.EQ.ILO) GO TO 100
      IF (X.LT.XT(MIDDLE)) GO TO 80
      ILO = MIDDLE
      GO TO 70
   80 IHI = MIDDLE
      GO TO 70
   90 MFLAG = -1
      ILEFT = 1
      RETURN
  100 MFLAG = 0
      ILEFT = ILO
      RETURN
  110 MFLAG = 1
      ILEFT = LXT
      RETURN
      END
      DOUBLE PRECISION FUNCTION ENLOST (ZP,MP,ZT,MT,ENERGY)
      implicit real*8 (a-h,o-z)
      INTEGER*4  ZP
      INTEGER*4  ZT
      REAL*8     MT, MP
      REAL*8     ENERGY
      REAL*8     BOHR
      REAL*8     FIRSOV
      REAL*8     DEBOHR,DEFIRSO,DETIT
      BOHR = DEBOHR(ZP,MP,ZT,MT,ENERGY)
      FIRSOV = DEFIRSO(ZP,ZT,MT,ENERGY)
      ENLOST = DMIN1(BOHR,FIRSOV)
      RETURN
      END
      double precision function enlostv (energy)
      implicit real*8 (a-h,o-z)
      include 'atimasys.inc'
      double precision energy
      integer i
      double precision  m
      double precision  enlost
      m = 0.d0
      do i = 1, nnu        m = m + anuc(i) * mt(i)
      enddo
      enlostv = 0.d0
      do i = 1, nnu        enlostv = enlostv + (anuc(i)*mt(i)/m) *
     &           enlost(zp,mp,zt(i),mt(i),energy)
      enddo
      return
      end
      SUBROUTINE ERRCHK(IMESH,XIOLD,ALDIF,VALSTR,WORK,MSTAR,IFIN)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ERR(40),Z(40),ERREST(40) ,DMVAL(20)
      COMMON /ORDER/K,NCOMP,MSTR,KD,KDM,MNSUM,M(20)
!$OMP THREADPRIVATE(/ORDER/)
      COMMON /APPR/ N,NOLD,NMAX,NALPHA,MSHFLG,MSHNUM,MSHLMT,MSHALT
!$OMP THREADPRIVATE(/APPR/)
      COMMON /SIDE/  ZETA(40), ALEFT, ARIGHT, IZETA, IWR
!$OMP THREADPRIVATE(/SIDE/)
      COMMON /ERRORS/ TOL(40),WGTMSH(40),TOLIN(40),ROOT(40),&
            JTOL(40),LTOL(40),NTOL
!$OMP THREADPRIVATE(/ERRORS/)
      COMMON /COLLOC/ RHO(7),WGTERR(40)
!$OMP THREADPRIVATE(/COLLOC/)
      COMMON /NONLN/ PRECIS,NONLIN,ITER,LIMIT,ICARE,IPRINT,IGUESS,IFREEZ
!$OMP THREADPRIVATE(/NONLN/)
      COMMON /BSPLIN/ VNCOL(66,7), VNSAVE(66,5), VN(66)
!$OMP THREADPRIVATE(/BSPLIN/)
      DIMENSION XIOLD(*), ALDIF(*), VALSTR(*), WORK(MSTAR,1)
      IFIN = 1
      NOLDP1 = NOLD + 1
      IF ( IPRINT .GE. 0 )                          GO TO 30
      DO 10 I = 1, NOLD
   10 CALL APPROX(I, XIOLD(I), WORK(1,I), VNSAVE(1,1), XIOLD,&
          NOLD, ALDIF, K, NCOMP, M, MSTAR, 3, DUMM, 0)
      CALL APPROX(NOLD, XIOLD(NOLDP1), WORK(1,NOLDP1), VN, XIOLD,&
          NOLD, ALDIF, K, NCOMP, M, MSTAR, 2, DUMM, 0)
      DO 20 I = 1, MSTAR
           WRITE(IWR,140) I
   20 WRITE(IWR,150) (WORK(I,J), J=1, NOLDP1)
   30 CONTINUE
      IF (IMESH.EQ.1)                               RETURN
      DO 40 J = 1,MSTAR
   40 ERREST(J) = 0.D0
      DO 100 IBACK = 1,NOLD
           I = NOLD +1 -IBACK
           MSHFLG = 1
           DO 50 J = 1,MSTAR
             Z(J) = 0.D0
   50      ERR(J) = 0.D0
           DO 60 J = 1,2
             JJ = 5 - J
             KNEW = ( 4 * (I-1) + 3 - J ) * MSTAR + 1
             KSTORE = ( 2 * (I-1) + 2 - J ) * MSTAR + 1
             X = XIOLD(I) + DFLOAT(3-J)/3.D0*(XIOLD(I+1)-XIOLD(I))
             CALL APPROX (I, X, VALSTR(KNEW), VNSAVE(1,JJ), XIOLD,&
                 NOLD, ALDIF, K, NCOMP, M, MSTAR, 3,DUMM,0)
             DO 60 L = 1,MSTAR
               ERR(L) = ERR(L) + WGTERR(L) * DABS(VALSTR(KNEW) -&
              VALSTR(KSTORE))
               Z(L) = Z(L) + .5D0 * DABS(VALSTR(KNEW))
               KNEW = KNEW + 1
               KSTORE = KSTORE + 1
   60      CONTINUE
           IF (IFIN .EQ. 0)                         GO TO 80
           DO 70 J = 1, NTOL
             LTOLJ = LTOL(J)
   70      IF ( ERR(LTOLJ) .GT. TOLIN(J) * (Z(LTOLJ)+1.D0) )  IFIN = 0
   80      DO 90 L = 1,MSTAR
   90      ERREST(L) = DMAX1(ERREST(L),ERR(L))
  100 CONTINUE
      IF (IPRINT .LT. 1) WRITE(IWR,130)
      LJ = 1
      DO 110 J = 1,NCOMP
           MJ = LJ - 1 + M(J)
           IF (IPRINT .LT. 1) WRITE(IWR,120) J, (ERREST(L), L= LJ, MJ)
           LJ = MJ + 1
  110 CONTINUE
      RETURN
  120 FORMAT (3H U(, I2, 3H) -,4D12.4)
  130 FORMAT (/26H THE ESTIMATED ERRORS ARE,)
  140 FORMAT( 19H MESH VALUES FOR Z(, I2, 2H), )
  150 FORMAT(1H , 8D15.7)
      END
      real*8 FUNCTION e_out(t,bcoef,n,k,range,ei,thick)
      use ion_track_structure, only : lflgTS
      implicit real*8 (a-h, o-z)
       include 'atimasys.inc'
      real*8 t(*), bcoef(*)
      PARAMETER (ITMAX=100)!,TOL=1.e-6)!,EPS=3.E-8)
      FUNC(X) = range - bvalue(t,bcoef,n,k,x,0) - thick
      EPS=3.E-8
      TOL=1.e-6
      if(lflgTS .eq. 1) then
          EPS=EPS*1.d-4
          TOL=TOL*1.e-4
      endif
      e_out = 0.0d0
      if( thick .ge. range ) return
      A=ei
      B=0.d0
      FA=FUNC(A)
      FB=FUNC(B)
      FC=FB
      DO 11 ITER=1,ITMAX
        IF(FB*FC.GT.0.) THEN
          C=A
          FC=FA
          D=B-A
          E=D
        ENDIF
        IF(ABS(FC).LT.ABS(FB)) THEN
          A=B
          B=          C=A
          FA=FB
          FB=F          FC=FA
        ENDIF
        TOL1=2.*EPS*ABS(B)+0.5*TOL
        XM=.5*(C-B)
        IF(ABS(XM).LE.TOL1 .OR. FB.EQ.0.)THEN
          if( ITER .eq. 2 ) then
            S=FB/FA
            IF(A.EQ.C) THEN
              P=2.*XM*S
              Q=1.-S
            ELSE
              Q=FA/F              R=FB/F              P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
              Q=(Q-1.)*(R-1.)*(S-1.)
            ENDIF
            IF(P.GT.0.) Q=-Q
            P=ABS(P)
            IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
              E=D
              D=P/Q
            ELSE
              D=XM
              E=D
            ENDIF
            B=B+D
          endif
          e_out = b
          RETURN
        ENDIF
        IF(ABS(E).GE.TOL1 .AND. ABS(FA).GT.ABS(FB)) THEN
          S=FB/FA
          IF(A.EQ.C) THEN
            P=2.*XM*S
            Q=1.-S
          ELSE
            Q=FA/F            R=FB/F            P=S*(2.*XM*Q*(Q-R)-(B-A)*(R-1.))
            Q=(Q-1.)*(R-1.)*(S-1.)
          ENDIF
          IF(P.GT.0.) Q=-Q
          P=ABS(P)
          IF(2.*P .LT. MIN(3.*XM*Q-ABS(TOL1*Q),ABS(E*Q))) THEN
            E=D
            D=P/Q
          ELSE
            D=XM
            E=D
          ENDIF
        ELSE
          D=XM
          E=D
        ENDIF
        A=B
        FA=FB
        IF(ABS(D) .GT. TOL1) THEN
          B=B+D
        ELSE
          B=B+SIGN(TOL1,XM)
        ENDIF
        FB=FUNC(B)
11    CONTINUE
          e_out = b
      RETURN
      END
      SUBROUTINE FACTRB ( W, IPIVOT, D, NROW, NCOL, LAST, IFLAG )
      INTEGER IPIVOT(NROW),NCOL,LAST,IFLAG, I,IPIVI,IPIVK,J,K,KP1
      DOUBLE PRECISION W(NROW,NCOL),D(NROW), AWIKDI,COLMAX,RATIO,ROWMAX
      DOUBLE PRECISION DABS,DMAX1
      DO 20 I=1,NROW
           IPIVOT(I) = I
           ROWMAX = 0.D0
           DO 10 J=1,NCOL
   10      ROWMAX = DMAX1(ROWMAX, DABS(W(I,J)))
           IF (ROWMAX .EQ. 0.D0)                    GO TO 90
   20 D(I) = ROWMAX
      K = 1
   30      IPIVK = IPIVOT(K)
           IF (K .EQ. NROW)                         GO TO 80
           J = K
           KP1 = K+1
           COLMAX = DABS(W(IPIVK,K))/D(IPIVK)
           DO 40 I=KP1,NROW
             IPIVI = IPIVOT(I)
             AWIKDI = DABS(W(IPIVI,K))/D(IPIVI)
             IF (AWIKDI .LE. COLMAX)                GO TO 40
             COLMAX = AWIKDI
             J = I
   40      CONTINUE
           IF (J .EQ. K)                            GO TO 50
           IPIVK = IPIVOT(J)
           IPIVOT(J) = IPIVOT(K)
           IPIVOT(K) = IPIVK
           IFLAG = -IFLAG
   50      CONTINUE
           IF (DABS(W(IPIVK,K))+D(IPIVK) .LE. D(IPIVK))&
                                                   GO TO 90
           DO 60 I=KP1,NROW
             IPIVI = IPIVOT(I)
             W(IPIVI,K) = W(IPIVI,K)/W(IPIVK,K)
             RATIO = -W(IPIVI,K)
             DO 60 J=KP1,NCOL
   60      W(IPIVI,J) = RATIO*W(IPIVK,J) + W(IPIVI,J)
           K = KP1
           IF (K .LE. LAST)                         GO TO 30
      RETURN
   80 IF(DABS(W(IPIVK,NROW))+D(IPIVK) .GT. D(IPIVK)) RETURN
   90 IFLAG = 0
      RETURN
      END
      SUBROUTINE FCBLOK ( BLOKS, INTEGS, NBLOKS, IPIVOT, SCRTCH, IFLAG )
      INTEGER INTEGS(3,NBLOKS),IPIVOT(*),IFLAG, I,INDEX,INDEXB,INDEXN,&
             LAST,NCOL,NROW
      DOUBLE PRECISION BLOKS(*),SCRTCH(*)
      IFLAG = 1
      INDEXB = 1
      INDEXN = 1
      I = 1
   10      INDEX = INDEXN
           NROW = INTEGS(1,I)
           NCOL = INTEGS(2,I)
           LAST = INTEGS(3,I)
           CALL FACTRB ( BLOKS(INDEX), IPIVOT(INDEXB), SCRTCH, NROW,&
               NCOL, LAST, IFLAG)
           IF (IFLAG .EQ. 0 .OR. I .EQ. NBLOKS)     RETURN
           I = I+1
           INDEXN = NROW*NCOL + INDEX
           CALL SHIFTB ( BLOKS(INDEX), IPIVOT(INDEXB), NROW, NCOL,&
               LAST, BLOKS(INDEXN), INTEGS(1,I), INTEGS(2,I) )
           INDEXB = INDEXB + NROW
      GO TO 10
      END
      SUBROUTINE FDUMP
      RETURN
      END
      SUBROUTINE GF(I,Z,G)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Z(1)
      G=Z(1)
      RETURN
      END
      SUBROUTINE HORDER (I, UHIGH, XIOLD, ALDIF)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /APPR/ N,NOLD,NMAX,NALPHA,MSHFLG,MSHNUM,MSHLMT,MSHALT
!$OMP THREADPRIVATE(/APPR/)
      COMMON /ORDER/ K,NCOMP,MSTAR,KD,KDM,MNSUM,M(20)
!$OMP THREADPRIVATE(/ORDER/)
      DIMENSION UHIGH(*) , AR(20), ARM1(20)
      DIMENSION ALDIF(*), XIOLD(*)
      DN2 = 1.D0 / (XIOLD(I+1) - XIOLD(I))
      INCOMP = 0
      DO 50 J = 1,NCOMP
           MJ = M(J)
           NALPHJ = K * NOLD + MJ
           KPMJ = K + MJ
           KMR = K + 1
           MJM1 = MJ - 1
           INCOMP = INCOMP + MJM1 * NALPHJ
           LEFT = I * K + MJ - KMR
           DO 10 L=1,KMR
             LEFTPL = LEFT + L
   10      ARM1(L+MJ-1) = ALDIF(INCOMP+LEFTPL)
           INCOMP = INCOMP + NALPHJ
           KPMJ1 = KPMJ - 1
           DO 40 NR = MJ,KPMJ1
             KMR = KMR - 1
             DNK2 = DN2 * DFLOAT(KMR)
             DO 20 L = 1,KMR
   20        AR(L+NR) = DNK2 * (ARM1(L+NR) - ARM1(L+NR-1))
             DO 30 L=NR,KPMJ1
   30        ARM1(L+1) = AR(L+1)
   40      CONTINUE
           UHIGH(J) = AR(KPMJ)
   50 CONTINUE
      RETURN
      END
      INTEGER FUNCTION I1MACH(I)
      INTEGER IMACH(16),OUTPUT
      EQUIVALENCE (IMACH(4),OUTPUT)
       DATA IMACH( 1) /    5 /
       DATA IMACH( 2) /    6 /
       DATA IMACH( 3) /    7 /
       DATA IMACH( 4) /    6 /
       DATA IMACH( 5) /   32 /
       DATA IMACH( 6) /    4 /
       DATA IMACH( 7) /    2 /
       DATA IMACH( 8) /   31 /
       DATA IMACH( 9) / 2147483647 /
       DATA IMACH(10) /    2 /
       DATA IMACH(11) /   24 /
       DATA IMACH(12) / -125 /
       DATA IMACH(13) /  128 /
       DATA IMACH(14) /   53 /
       DATA IMACH(15) / -1021 /
       DATA IMACH(16) /  1024 /
      IF (I .LT. 1  .OR.  I .GT. 16)&
        CALL XERROR ( 'I1MACH -- I OUT OF BOUNDS',25,1,2)
      I1MACH=IMACH(I)
      RETURN
      END
      SUBROUTINE INTEGR(FF,A,B,C,TOL,IPRINT,FSPACE,ISPACE,NDIMF,NDIMI,&
                       IFLAG)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION ZETA(1),IPAR(11),FSPACE(NDIMF),ISPACE(NDIMI),Z(1)
      EXTERNAL FF,DFF,GF,DGF,SF
*********** 1. PARAMETER FUER COLSYS
      M=1
      ZETA(1)=
      DO 10 I=1,11
   10 IPAR(I)=0
      IPAR(4)=1
      IPAR(5)=NDIMF
      IPAR(6)=NDIMI
      IPAR(7)=IPRINT
      LTOL=1
*********** 2. COLSYS AUFRUF
      CALL COLSYS(1,M,A,B,ZETA,IPAR,LTOL,TOL,FIXPNT,ISPACE,FSPACE,
     &            IFLAG,FF,DFF,GF,DGF,SF)
       RETURN
       END
      subroutine interv ( xt, lxt, x, left, mflag )
      integer*4 left,lxt,mflag,   ihi,ilo,istep,middle
      real*8 x,xt(lxt)
      data ilo /1/
      save ilo
!$OMP THREADPRIVATE(ilo)
      ihi = ilo + 1
      if (ihi .lt. lxt)                 go to 20
         if (x .ge. xt(lxt))            go to 110
         if (lxt .le. 1)                go to 90
         ilo = lxt - 1
         ihi = lxt
   20 if (x .ge. xt(ihi))               go to 40
      if (x .ge. xt(ilo))               go to 100
      istep = 1
   31    ihi = ilo
         ilo = ihi - istep
         if (ilo .le. 1)                go to 35
         if (x .ge. xt(ilo))            go to 50
         istep = istep*2
                                        go to 31
   35 ilo = 1
      if (x .lt. xt(1))                 go to 90
                                        go to 50
   40 istep = 1
   41    ilo = ihi
         ihi = ilo + istep
         if (ihi .ge. lxt)              go to 45
         if (x .lt. xt(ihi))            go to 50
         istep = istep*2
                                        go to 41
   45 if (x .ge. xt(lxt))               go to 110
      ihi = lxt
   50 middle = (ilo + ihi)/2
      if (middle .eq. ilo)              go to 100
      if (x .lt. xt(middle))            go to 53
         ilo = middle
                                        go to 50
   53    ihi = middle
                                        go to 50
   90 mflag = -1
      left = 1
                                        return
  100 mflag = 0
      left = ilo
                                        return
  110 mflag = 1
	  if (x .eq. xt(lxt)) mflag = 0
      left = lxt
  111 if (left .eq. 1)                  return
	  left = left - 1
	  if (xt(left) .lt. xt(lxt))        return
										go to 111
      end
      FUNCTION J4SAVE(IWHICH,IVALUE,ISET)
      LOGICAL ISET
      INTEGER IPARAM(9)
      SAVE IPARAM
!$OMP THREADPRIVATE(IPARAM)
      DATA IPARAM(1),IPARAM(2),IPARAM(3),IPARAM(4)/0,2,0,10/
      DATA IPARAM(5)/1/
      DATA IPARAM(6),IPARAM(7),IPARAM(8),IPARAM(9)/0,0,0,0/
      J4SAVE = IPARAM(IWHICH)
      IF (ISET) IPARAM(IWHICH) = IVALUE
      RETURN
      END
      subroutine knots ( break, l, m, kpm, t, n )
      integer*4 l,kpm,n,   iside,j,jj,jjj,k,ll,m
      real*8 break(*),t(*)
      k = kpm - m
      n = l*k + m
      jj = n + kpm
      jjj = l + 1
      do 11 ll=1,kpm
         t(jj) = break(jjj)
   11    jj = jj - 1
      do 12 j=1,l
      jjj = jjj - 1
         do 12 ll=1,k
            t(jj) = break(jjj)
   12       jj = jj - 1
      do 13 ll=1,kpm
   13    t(ll) = break(1)
                                        return
      end
      SUBROUTINE LSYSLV (IFLAG, XI, XIOLD, XIJ, ALPHA, ALDIF,&
                RHS, ALPHO, A, IPIV, INTEGS, RNORM,&
                MODE, FSUB, DFSUB, GSUB, DGSUB, SOLUTN)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /ORDER/ K, NCOMP, MSTAR, KD, KDM, MNSUM, M(20)
!$OMP THREADPRIVATE(/ORDER/)
      COMMON /SIDE/  ZETA(40), ALEFT, ARIGHT, IZETA, IWR
!$OMP THREADPRIVATE(/SIDE/)
      COMMON /BSPLIN/  VNCOL(66,7), VNSAVE(66,5), VN(66)
!$OMP THREADPRIVATE(/BSPLIN/)
      COMMON /APPR/ N,NOLD,NMAX,NALPHA,MSHFLG,MSHNUM,MSHLMT,MSHALT
!$OMP THREADPRIVATE(/APPR/)
      COMMON /NONLN/ PRECIS,NONLIN,ITER,LIMIT,ICARE,IPRINT,IGUESS,IFREEZ
!$OMP THREADPRIVATE(/NONLN/)
      COMMON /HI/   DN1, DN2, DN3
!$OMP THREADPRIVATE(/HI/)
      EXTERNAL DFSUB, DGSUB
      DIMENSION  ALPHO(*), XI(*), XIOLD(*), XIJ(*), ALPHA(*)
      DIMENSION ALDIF(*), RHS(*), A(*), IPIV(*), INTEGS(3,*)
      DIMENSION Z(40), F(40), DF(800), DMVAL(20)
      M1 = MODE + 1
      GO TO (10, 30, 30, 310), M1
   10 DO 20 I=1,MSTAR
   20 Z(I) = 0.D0
   30 IRHS = 0
      IA = 1
      IZETA = 1
      LSIDE = 0
      RNORM = 0.D0
      IALPHO = 0
      IF (ITER .GE. 1 .OR. MODE .EQ. 2)             GO TO 80
      DO 70 I = 1,N
           INTEGS(2,I) = KDM
           IF (I .LT. N)                            GO TO 40
           INTEGS(3,I) = KDM
           LSIDE = MSTAR
           GO TO 60
   40      INTEGS(3,I) = KD
   50      IF( LSIDE .EQ. MSTAR )                   GO TO 60
           IF ( ZETA(LSIDE+1) .GE. XI(I+1) )        GO TO 60
           LSIDE = LSIDE + 1
           GO TO 50
   60      NROW = KD + LSIDE
   70 INTEGS(1,I) = NROW
   80 DO 290 I=1,N
           XIL = XI(1)
           IF (I .GT. 1) XIL = XI(I-1)
           XIR = XI(N+1)
           IF (I .LT. N) XIR = XI(I+2)
           DN1 = 1.D0 / (XI(I+1) - XIL)
           DN2 = 1.D0 / (XI(I+1) - XI(I))
           DN3 = 1.D0 / (XIR - XI(I))
           NROW = INTEGS(1,I)
           II = I
           ICOLC = (I-1) * K
           ID = IRHS + IZETA - 1
           DO 270  LL=1,K
             XX = XIJ (ICOLC+LL)
  100        IF ( IZETA .GT. MSTAR )                GO TO 160
             IF ( ZETA(IZETA) .GE. XX )             GO TO 160
  110        ID = ID + 1
             IALPHO = IALPHO + 1
             IF (MODE .EQ. 0)                       GO TO 130
             IF (IGUESS .NE. 1)                     GO TO 120
             CALL SOLUTN (ZETA(IZETA), Z, DMVAL)
             GO TO (130, 140), MODE
  120        CALL APPROX (II, ZETA(IZETA), Z, VN, XIOLD, NOLD,&
                 ALDIF, K, NCOMP, M, MSTAR, 1, DUMMY, 0)
             IF (MODE .EQ. 2)                       GO TO 140
  130        CALL GSUB (IZETA, Z, G)
             RHS(ID) = -G
             RNORM = RNORM + G**2
             IF (MODE .EQ. 1)                       GO TO 150
  140        CALL BLDBLK (I, ZETA(IZETA), LL, A(IA), NROW, ID-IRHS, Z,&
                 DF, NCOMP, XI, ALPHO, IALPHO, 1, DFSUB, DGSUB)
  150        IZETA = IZETA + 1
             IF (IZETA .GT. MSTAR .AND.&
                ZETA(MSTAR) .GE. DMIN1(XX,ARIGHT)) GO TO 280
             IF (XX .GT. XI(N+1))                   GO TO 260
             GO TO 100
  160        IF (IGUESS .NE. 1)  GO TO (210, 170, 230), M1
             CALL SOLUTN(XX, Z, DMVAL)
             GO TO (190, 250), MODE
  170        IF (ITER .GE. 1 )                      GO TO 180
             CALL APPROX (II, XX, Z, VN, XIOLD, NOLD, ALDIF, K,&
                 NCOMP, M, MSTAR, 1, DMVAL, 1)
             GO TO 190
  180        CALL APPROX (I, XX, Z, VNCOL(1,LL), XIOLD, NOLD,&
                 ALDIF, K, NCOMP, M, MSTAR, 3, DMVAL, 1)
  190        CALL FSUB (XX, Z, F)
             DO 200 J=1,NCOMP
               ID = ID + 1
               VALUE = DMVAL(J) - F(J)
               RHS(ID) = -VALUE
               RNORM = RNORM + VALUE**2
               IF (ITER .GE. 1)                     GO TO 200
               IALPHO = IALPHO + 1
               ALPHO(IALPHO) = DMVAL(J)
  200        CONTINUE
             GO TO 260
  210        CALL FSUB (XX, Z, F)
             DO 220 J=1,NCOMP
               ID = ID + 1
  220        RHS(ID) = F(J)
             ID = ID - NCOMP
             GO TO 250
  230        IF (ITER .GE. 1 )                      GO TO 240
             CALL APPROX (II, XX, Z, VN, XIOLD, NOLD, ALDIF, K,&
                 NCOMP, M, MSTAR, 1, DUMMY, 0)
             GO TO 250
  240        CALL APPROX (I, XX, Z, VNCOL(1,LL), XIOLD, NOLD,&
                 ALDIF, K, NCOMP, M, MSTAR, 3, DUMMY, 0)
  250        CALL BLDBLK (I, XX, LL, A(IA), NROW, ID-IRHS+1, Z,&
                 DF, NCOMP, XI, ALPHO, IALPHO, 2, DFSUB, DGSUB)
             ID = ID + NCOMP
  260        IF (LL .LT. K)                         GO TO 270
             IF (I .LT. N .OR. IZETA .GT. MSTAR)    GO TO 280
             XX = XI(N+1) + 1.D0
             GO TO 110
  270      CONTINUE
  280      IRHS = IRHS + NROW
           IA = IA + NROW * KDM
  290 CONTINUE
      IF (MODE .NE. 1)                              GO TO 300
      RNORM = DSQRT(RNORM / DFLOAT(NALPHA))
      RETURN
  300 CALL FCBLOK (A, INTEGS, N, IPIV, ALPHA, IFLAG)
      IF(IFLAG .EQ. 0)                              RETURN
  310 CALL SBBLOK (A, INTEGS, N, IPIV, RHS, ALPHA)
      IF (ITER .GE. 1 .OR. MODE .NE. 2)             RETURN
      IALPHO = 0
      IRHS = 0
      ISTO = 0
      DO 320 I=1,N
           NROW = INTEGS(1,I)
           IRHS = IRHS + ISTO
           ISTART = ISTO + 1
           ISTO = NROW - KD
           DO 320 J=ISTART,NROW
             IRHS = IRHS + 1
             IALPHO = IALPHO + 1
             RHS(IRHS) = RHS(IRHS) + ALPHO(IALPHO)
  320 CONTINUE
      CALL SBBLOK (A, INTEGS, N, IPIV, RHS, ALPHO)
      DO 330 I=1,NALPHA
  330 ALPHO(I) = ALPHO(I) - ALPHA(I)
      RETURN
      END
      SUBROUTINE NEWMSH (MODE, XI, XIOLD, XIJ, ALDIF, VALSTR,&
                SLOPE, ACCUM, NFXPNT, FIXPNT)
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /ORDER/K,NCOMP,MSTAR,KD,KDM,MNSUM,M(20)
!$OMP THREADPRIVATE(/ORDER/)
      COMMON /APPR/ N,NOLD,NMAX,NALPHA,MSHFLG,MSHNUM,MSHLMT,MSHALT
!$OMP THREADPRIVATE(/APPR/)
      COMMON /ERRORS/ TOL(40),WGTMSH(40),TOLIN(40),ROOT(40),&
            JTOL(40),LTOL(40),NTOL
!$OMP THREADPRIVATE(/ERRORS/)
      COMMON /COLLOC/ RHO(7),WGTERR(40)
!$OMP THREADPRIVATE(/COLLOC/)
      COMMON /SIDE/  ZETA(40), ALEFT, ARIGHT, IZETA, IWR
!$OMP THREADPRIVATE(/SIDE/)
      COMMON /NONLN/ PRECIS,NONLIN,ITER,LIMIT,ICARE,IPRINT,IGUESS,IFREEZ
!$OMP THREADPRIVATE(/NONLN/)
      COMMON /BSPLIN/ VNCOL(66,7), VNSAVE(66,5), VN(66)
!$OMP THREADPRIVATE(/BSPLIN/)
      DIMENSION D1(40), D2(40), ZV(40), SLOPE(*), ACCUM(*), VALSTR(*)
      DIMENSION XI(*), XIOLD(*), XIJ(*), ALDIF(*), FIXPNT(*)
      NFXP1 = NFXPNT +1
      GO TO (180, 100, 50, 20, 10), MODE
   10 MSHLMT = 1
   20 IF (IGUESS .LT. 2)                            GO TO 40
      NOLDP1 = NOLD + 1
      IF (IPRINT .LT. 1)  WRITE(IWR,360) NOLD, (XIOLD(I), I=1,NOLDP1)
      IF (IGUESS .NE. 3)                            GO TO 40
      N = NOLD /2
      I = 0
      DO 30 J = 1, NOLD, 2
           I = I + 1
   30 XI(I) = XIOLD(J)
   40 CONTINUE
      NP1 = N + 1
      XI(1) = ALEFT
      XI(NP1) = ARIGHT
      GO TO 320
   50 IF ( N .LT. NFXP1 ) N = NFXP1
      NP1 = N + 1
      XI(1) = ALEFT
      ILEFT = 1
      XLEFT = ALEFT
      DO 90 J = 1,NFXP1
           XRIGHT = ARIGHT
           IRIGHT = NP1
           IF ( J .EQ. NFXP1 )                      GO TO 60
           XRIGHT = FIXPNT(J)
           NMIN = (XRIGHT-ALEFT) / (ARIGHT-ALEFT) * DFLOAT(N) + 1.5D0
           IF (NMIN .GT. N-NFXPNT+J)  NMIN = N - NFXPNT + J
           IRIGHT = MAX0 (ILEFT+1, NMIN)
   60      XI(IRIGHT) = XRIGHT
           NREGN = IRIGHT - ILEFT - 1
           IF ( NREGN .EQ. 0 )                      GO TO 80
           DX = (XRIGHT - XLEFT) / DFLOAT(NREGN+1)
           DO 70 I = 1,NREGN
   70      XI(ILEFT+I) = XLEFT + DFLOAT(I) * DX
   80      ILEFT = IRIGHT
           XLEFT = XRIGHT
   90 CONTINUE
      GO TO 320
  100 N2 = 2 * N
      IF (N2 .LE. NMAX)                             GO TO 120
      IF (MODE .EQ. 2)                              GO TO 110
      N = NMAX / 2
      GO TO 220
  110 IF (IPRINT .LT. 1)  WRITE(IWR,370)
      N = N2
      RETURN
  120 IF (MSHFLG .EQ. 0)                            GO TO 140
      KSTORE = 1
      DO 130 I = 1,NOLD
           HD6 = (XIOLD(I+1) - XIOLD(I)) / 6.D0
           X = XIOLD(I) + HD6
           CALL APPROX (I, X, VALSTR(KSTORE), VNSAVE(1,2), XIOLD,&
               NOLD, ALDIF, K, NCOMP, M, MSTAR, 3,DUMM,0)
           X = X + 4.D0 * HD6
           KSTORE = KSTORE + 3 * MSTAR
           CALL APPROX (I, X, VALSTR(KSTORE), VNSAVE(1,5),&
               XIOLD, NOLD, ALDIF, K, NCOMP, M, MSTAR, 3,DUMM,0)
           KSTORE = KSTORE + MSTAR
  130 CONTINUE
      GO TO 160
  140 KSTORE = 1
      DO 150 I = 1,N
           X = XI(I)
           HD6 = (XI(I+1) - XI(I)) / 6.D0
           DO 150 J = 1,4
             X = X + HD6
             IF ( J.EQ.3 ) X = X + HD6
             CALL APPROX (I, X, VALSTR(KSTORE), VNSAVE(1,J+1),&
                 XIOLD, NOLD, ALDIF, K, NCOMP, M, MSTAR, 3,&
                 DUMM,0)
             KSTORE = KSTORE + MSTAR
  150 CONTINUE
  160 MSHFLG = 0
      MSHNUM = 1
      MODE = 2
      J = 2
      DO 170 I = 1,N
           XI(J) = (XIOLD(I) + XIOLD(I+1)) / 2.D0
           XI(J+1) = XIOLD(I+1)
  170 J = J + 2
      N = N2
      GO TO 320
  180 IF ( NOLD .EQ. 1 )                            GO TO 100
      IF (NOLD .LE. 2*NFXPNT)                       GO TO 100
      I = 1
      CALL HORDER (1, D1, XIOLD, ALDIF)
      CALL HORDER (2, D2, XIOLD, ALDIF)
      CALL APPROX (I, XIOLD(I), ZV, VNSAVE(1,1), XIOLD, NOLD,&
          ALDIF, K, NCOMP, M, MSTAR, 3, DUMM, 0)
      ACCUM(1) = 0.D0
      SLOPE(1) = 0.D0
      ONEOVH = 2.D0 / ( XIOLD(3) - XIOLD(1) )
      DO 190 J = 1,NTOL
           JJ = JTOL(J)
           JV = LTOL(J)
  190 SLOPE(1) = DMAX1(SLOPE(1),(DABS(D2(JJ)-D1(JJ))*WGTMSH(J)*&
      ONEOVH / (1.D0 + DABS(ZV(JV)))) **ROOT(J))
      SLPHMX = SLOPE(1) * (XIOLD(2) - XIOLD(1))
      ACCUM(2) = SLPHMX
      IFLIP = 1
      DO 210 I = 2,NOLD
           IF ( IFLIP .EQ. (-1) ) CALL HORDER ( I, D1, XIOLD, ALDIF)
           IF ( IFLIP .EQ. 1 ) CALL HORDER ( I, D2, XIOLD, ALDIF)
           CALL APPROX (I, XIOLD(I), ZV, VNSAVE(1,1), XIOLD, NOLD,&
               ALDIF, K, NCOMP, M, MSTAR, 3, DUMM, 0)
           ONEOVH = 2.D0 / ( XIOLD(I+1) - XIOLD(I-1) )
           SLOPE(I) = 0.D0
           DO 200 J = 1,NTOL
             JJ = JTOL(J)
             JV = LTOL(J)
  200      SLOPE(I) = DMAX1(SLOPE(I),(DABS(D2(JJ)-D1(JJ))*WGTMSH(J)*&
          ONEOVH / (1.D0 + DABS(ZV(JV)))) **ROOT(J))
           TEMP = SLOPE(I) * (XIOLD(I+1)-XIOLD(I))
           SLPHMX = DMAX1(SLPHMX,TEMP)
           ACCUM(I+1) = ACCUM(I) + TEMP
  210 IFLIP = - IFLIP
      AVRG = ACCUM(NOLD+1) / DFLOAT(NOLD)
      DEGEQU = AVRG / DMAX1(SLPHMX,PRECIS)
      NACCUM = ACCUM(NOLD+1) + 1.D0
      IF (IPRINT .LT. 0)  WRITE(IWR,350) DEGEQU, NACCUM
      IF (AVRG .LT. PRECIS)                         GO TO 100
      IF (DEGEQU .GE. .5D0)                         GO TO 100
      NMX = MAX0 (NOLD+1, NACCUM) / 2
      NMAX2 = NMAX / 2
      N = MIN0 (NMAX2, NOLD, NMX)
  220 NOLDP1 = NOLD + 1
      IF (N .LT. NFXP1) N=NFXP1
      MSHNUM = MSHNUM + 1
      IF (N .LT. NOLD) MSHNUM = MSHLMT
      IF (N .GT. NOLD/2)  MSHALT = 1
      IF (N .EQ. NOLD/2)  MSHALT = MSHALT + 1
      MSHFLG = 0
      IN = 1
      ACCL = 0.D0
      LOLD =2
      XI(1) = ALEFT
      XI(N+1) = ARIGHT
      DO 310 I = 1, NFXP1
           IF (I .EQ. NFXP1)                        GO TO 250
           DO 230 J = LOLD, NOLDP1
             LNEW = J
             IF (FIXPNT(I) .LE. XIOLD(J))           GO TO 240
  230      CONTINUE
  240      CONTINUE
           ACCR = ACCUM(LNEW) + (FIXPNT(I)-XIOLD(LNEW))*SLOPE(LNEW-1)
           NREGN = (ACCR-ACCL) / ACCUM(NOLDP1) * DFLOAT(N) - .5D0
           NREGN = MIN0(NREGN, N - IN - NFXP1 + I)
           XI(IN + NREGN + 1) = FIXPNT(I)
           GO TO 260
  250      ACCR = ACCUM(NOLDP1)
           LNEW = NOLDP1
           NREGN = N - IN
  260      IF (NREGN .EQ. 0)                        GO TO 300
           TEMP = ACCL
           TSUM = (ACCR - ACCL) / DFLOAT(NREGN+1)
           DO 290 J = 1, NREGN
             IN = IN + 1
             TEMP = TEMP + TSUM
             DO 270 L = LOLD, LNEW
               LCARRY = L
               IF (TEMP .LE. ACCUM(L))              GO TO 280
  270        CONTINUE
  280        CONTINUE
             LOLD = LCARRY
  290      XI(IN) = XIOLD(LOLD-1) + (TEMP - ACCUM(LOLD-1)) /&
          SLOPE(LOLD-1)
  300      IN = IN + 1
           ACCL = ACCR
           LOLD = LNEW
  310 CONTINUE
      MODE = 1
  320 CONTINUE
      K2 = 1
      DO 330 I = 1,N
           H = (XI(I+1) - XI(I)) / 2.D0
           XM = (XI(I+1) + XI(I)) / 2.D0
           DO 330 J = 1,K
             XIJ(K2) = RHO(J) * H + XM
             K2 = K2 + 1
  330 CONTINUE
      NP1 = N + 1
      IF (IPRINT .LT. 1)  WRITE(IWR,340) N, (XI(I),I=1,NP1)
      NALPHA = N * K * NCOMP + MSTAR
      RETURN
  340 FORMAT(/17H THE NEW MESH (OF,I5,16H SUBINTERVALS), ,100(/8F12.6))
  350 FORMAT(/21H MESH SELECTION INFO,/30H DEGREE OF EQUIDISTRIBUTION =&
            , F8.5, 28H PREDICTION FOR REQUIRED N = , I8)
  360 FORMAT(/20H THE FORMER MESH (OF,I5,15H SUBINTERVALS),,&
            100(/8F12.6))
  370 FORMAT (/23H  EXPECTED N TOO LARGE  )
      END
       double precision function omthver(energy)
       implicit real*8 (a-h,o-z)
       double precision energy
       include 'atimasys.inc'
       double precision m
       integer i
       double precision omth
       m = 0.d0
       omthver = 0.d0
       do i = 1, nnu          m = m + anuc(i) * mt(i)
       enddo
       do i = 1, nnu         omthver = omthver + anuc(i) * omth(zp, mp, zt(i), mt(i),
     &       energy) * mt(i) / m
       enddo
       return
       end
       DOUBLE PRECISION FUNCTION OMTH(ZP,MP,ZT,MT,ENERGY)
       implicit real*8 (a-h,o-z)
       include 'atimacnt.inc'
       integer*4 zp
       integer*4 zt
       double precision mp
       double precision mt
       double precision energy
       double precision beta
       double precision p
       double precision lr
       double precision sirale
       p        = dsqrt( energy * (energy + 2.d0 * amu)) * mp
       beta     = p / ( (energy + amu) * mp )
       lr       = sirale(zt,mt)
       omth     = 198.81d0 * dble(zp)**2 / (lr * (p * beta)**2)
       return
       end
      DOUBLE PRECISION FUNCTION RPSTOP(Z,E)
      implicit real*8 (a-h,o-z)
      INTEGER*4  Z
      INTEGER*4  I,J
      REAL*8     VELPWR,PCOEF(8 ,97),PE0,PE,SH,SL,SP,E,B2
      DATA   ( (PCOEF(I,J), I= 1 ,8), J=1 ,9 )&
           / .0091827d0, .0053496d0, .69741d0, .48493d0,&
               316.07d0,1.0143d0, 9329.3d0, .0539890d0,&
             .11393d0, .0051984d0, 1.0822d0, .39252d0,&
             1081.0d0, 1.0645d0, 4068.5d0, .0176990d0,&
             .85837d0, .0050147d0, 1.6044d0, .38844d0,&
             1337.3d0, 1.047d0, 2659.2d0, .01898d0,&
             .8781d0,  .0051049d0, 5.4232d0, .2032d0,&
            1200.6d0, 1.0211d0, 1401.8d0, .0385290d0,&
            1.4608d0,  .0048836d0, 2.338d0,  .44249d0,&
            1801.3d0, 1.0352d0, 1784.1d0, .02024d0,&
            3.2579d0,  .0049148d0, 2.7156d0, .36473d0,&
            2092.2d0, 1.0291d0, 2643.6d0, .0182370d0,&
             1.59674d0, .0050837d0, 4.2073d0, .30612d0,&
             2394.2d0, 1.0255d0, 4892.1d0, .0160060d0,&
             1.75253d0, .0050314d0, 4.0824d0, .30067d0,&
            2455.8d0, 1.0181d0, 5069.7d0, .0174260d0,&
            1.226d0,   .0051385d0, 3.2246d0, .32703d0,&
            2525.0d0, 1.0142d0, 7563.6d0, .0194690d0/
      DATA   ( (PCOEF(I,J), I= 1 ,8), J=10 ,18)&
          / 1.0332d0,  .0051645d0, 3.004d0,  .33889d0,&
            2338.6d0, .99997d0, 6991.2d0, .0217990d0,&
            6.0972d0,  .0044292d0, 3.1929d0, .45763d0,&
           1363.3d0, .95182d0, 2380.6d0, .0818350d0,&
           14.013d0,   .0043646d0, 2.2641d0, .36326d0,&
           2187.4d0, .99098d0, 6264.8d0, .0462d0,&
             .039001d0, .0045415d0, 5.5463d0, .39562d0,&
            1589.2d0, .95316d0, 816.16d0, .0474840d0,&
            2.072d0,   .0044516d0, 3.5585d0, .53933d0,&
           1515.2d0, .93161d0, 1790.3d0, .0351980d0,&
           17.575d0,   .0038346d0, .078694d0, 1.2388d0,&
           2806.0d0, .97284d0, 1037.6d0, .0128790d0,&
           16.126d0,   .0038315d0, .054164d0, 1.3104d0,&
           2813.3d0, .96587d0, 1251.4d0, .0118470d0,&
            3.217d0,   .0044579d0, 3.6696d0, .5091d0,&
            2734.6d0, .96253d0, 2187.5d0, .0169070d0,&
            2.0379d0,  .0044775d0, 3.0743d0, .54773d0,&
            3505.0d0, .97575d0, 1714.00d0, .0117010d0/
      DATA   ( (PCOEF(I,J), I= 1 ,8), J=19 ,27)
     &    /   .74171d0, .0043051d0, 1.1515d0, .95083d0,
     &       917.21d0, .8782d0, 389.93d0, .18926d0,
     &       9.1316d0,  .0043809d0, 5.4611d0, .31327d0,
     &       3891.8d0, .97933d0, 6267.9d0, .0151960d0,
     &       7.2247d0,  .0043718d0, 6.1017d0, .37511d0,
     &       2829.2d0, .95218d0, 6376.1d0, .0203980d0,
     &        .147d0,   .0048456d0, 6.3485d0, .41057d0,
     &       2164.1d0, .94028d0, 5292.6d0, .0502630d0,
     &       5.0611d0,  .0039867d0, 2.6174d0, .57957d0,
     &       2218.9d0, .92361d0, 6323.00d0, .0256690d0,
     &        .53267d0, .0042968d0, .39005d0, 1.2725d0,
     &        1872.7d0, .90776d0, 64.166d0, .0301070d0,
     &        .47697d0, .0043038d0, .31452d0, 1.3289d0,
     &        1920.5d0, .90649d0, 45.576d0, .0274690d0,
     &        .027426d0, .0035443d0, .031563d0, 2.1755d0,
     &        1919.5d0, .90099d0, 23.902d0, .0253630d0,
     &        .16383d0, .0043042d0, .073454d0, 1.8592d0,
     &       1918.4d0, .89678d0, 27.61d0, .0231840d0/
      DATA   ( (PCOEF(I,J), I= 1 ,8), J=28 ,36)
     &    /  4.2562d0,  .0043737d0, 1.5606d0, .72067d0,
     &       1546.8d0, .87958d0, 302.02d0, .0409440d0,
     &       2.3508d0,  .0043237d0, 2.882d0,  .50113d0,
     &       1837.7d0, .89992d0, 2377.00d0, .04965d0,
     &       3.1095d0,  .0038455d0, .11477d0, 1.5037d0,
     &      2184.7d0, .89309d0, 67.306d0, .0165880d0,
     &      15.322d0,   .0040306d0, .65391d0, .67668d0,
     &      3001.7d0, .92484d0, 3344.2d0, .0163660d0,
     &       3.6932d0,  .0044813d0, 8.608d0,  .27638d0,
     &       2982.7d0, .9276d0, 3166.6d0, .0308740d0,
     &       7.1373d0,  .0043134d0, 9.4247d0, .27937d0,
     &       2725.8d0, .91597d0, 3166.1d0, .0250080d0,
     &       4.8979d0,  .0042937d0, 3.7793d0, .50004d0,
     &       2824.5d0, .91028d0, 1282.4d0, .0170610d0,
     &       1.3683d0,  .0043024d0, 2.5679d0, .60822d0,
     &       6907.8d0, .9817d0, 628.01d0, .0068055d0,
     &       1.8301d0,  .0042983d0, 2.9057d0, .6038d0,
     &       4744.6d0, .94722d0, 936.64d0, .0092242d0/
      DATA   ( (PCOEF(I,J), I= 1 ,8), J=37 ,45)
     &    /   .42056d0, .0041169d0, .01695d0, 2.3616d0,
     &      2252.7d0, .89192d0, 39.752d0, .0277570d0,
     &      30.78d0,    .0037736d0, .55813d0, .76816d0,
     &      7113.2d0, .97697d0, 1604.4d0, .0065268d0,
     &      11.576d0,   .0042119d0, 7.0244d0, .37764d0,
     &      4713.5d0, .94264d0, 2493.2d0, .01127d0,
     &       6.2406d0,  .0041916d0, 5.2701d0, .49453d0,
     &       4234.6d0, .93232d0, 2063.9d0, .0118440d0,
     &        .33073d0, .0041243d0, 1.7246d0, 1.1062d0,
     &        1930.2d0, .86907d0, 27.416d0, .0382080d0,
     &        .017747d0, .0041715d0, .14586d0, 1.7305d0,
     &       1803.6d0, .86315d0, 29.669d0, .0321230d0,
     &       3.7229d0,  .0041768d0, 4.6286d0, .56769d0,
     &       1678.0d0, .86202d0, 3094.00d0, .06244d0,
     &        .13998d0, .0041329d0, .25573d0, 1.4241d0,
     &        1919.3d0, .86326d0, 72.797d0, .0322350d0,
     &        .2859d0,  .0041386d0, .31301d0, 1.3424d0,
     &        1954.8d0, .86175d0, 115.18d0, .0293420d0/
      DATA   ( (PCOEF(I,J), I= 1 ,8), J=46 ,54)
     &     /  .76d0,    .0042179d0, 3.386d0,  .76285d0,
     &       1867.4d0, .85805d0, 69.994d0, .0364480d0,
     &       6.3957d0,  .0041935d0, 5.4689d0, .41378d0,
     &       1712.6d0, .85397d0,18493.00d0, .0564710d0,
     &       3.4717d0,  .0041344d0, 3.2337d0, .63788d0,
     &       1116.4d0, .81959d0, 4766.0d0, .1179d0,
     &       2.5265d0,  .0042282d0, 4.532d0,  .53562d0,
     &       1030.8d0, .81652d0,16252.0d0, .19722d0,
     &       7.3683d0,  .0041007d0, 4.6791d0, .51428d0,
     &       1160.0d0, .82454d0,17965.0d0, .13316d0,
     &       7.7197d0,  .004388d0, 3.242d0,   .68434d0,
     &      1428.1d0, .83398d0, 1786.7d0, .0665120d0,
     &      16.78d0,    .0041918d0, 9.3198d0, .29568d0,
     &      3370.9d0, .90289d0, 7431.7d0, .02616d0,
     &       4.2132d0,  .0042098d0, 4.6753d0, .57945d0,
     &       3503.9d0, .89261d0, 1468.9d0, .0143590d0,
     &       4.0818d0,  .004214d0, 4.4425d0,  .58393d0,
     &       3945.3d0, .90281d0, 1340.5d0, .0134140d0/
      DATA   ( (PCOEF(I,J), I= 1 ,8), J=55 ,63)
     &     /  .18517d0, .0036215d0, .00058788d0,3.5315d0,
     &       2931.3d0, .88936d0, 26.18d0, .0263930d0,
     &       4.8248d0,  .0041458d0, 6.0934d0, .57026d0,
     &       2300.1d0, .86359d0, 2980.7d0, .0386790d0,
     &        .49857d0, .0041054d0, 1.9775d0, .95877d0,
     &       786.55d0, .78509d0, 806.6d0, .40882d0,
     &       3.2754d0,  .0042177d0, 5.768d0,  .54054d0,
     &       6631.3d0, .94282d0, 744.07d0, .0083026d0,
     &       2.9978d0,  .0040901d0, 4.5299d0, .62025d0,
     &       2161.2d0, .85669d0, 1268.6d0, .0430310d0,
     &       2.8701d0,  .004096d0, 4.2568d0,  .6138d0,
     &      2130.4d0, .85235d0, 1704.1d0, .0393850d0,
     &      10.853d0,   .0041149d0, 5.8907d0, .46834d0,
     &      2857.2d0, .8755d0, 3654.2d0, .0299550d0,
     &       3.6407d0,  .0041782d0, 4.8742d0, .57861d0,
     &      1267.7d0, .82211d0, 3508.2d0, .24174d0,
     &      17.645d0,   .0040992d0, 6.5855d0, .32734d0,
     &      3931.3d0, .90754d0, 5156.7d0, .0362780d0/
      DATA   ( (PCOEF(I,J), I= 1 ,8), J=64 ,72)
     &    /  7.5309d0,  .0040814d0, 4.9389d0, .50679d0,
     &       2519.7d0, .85819d0, 3314.6d0, .0305140d0,
     &       5.4742d0,  .0040829d0, 4.897d0,  .51113d0,
     &       2340.1d0, .85296d0, 2342.7d0, .0356620d0,
     &       4.2661d0,  .0040667d0, 4.5032d0, .55257d0,
     &       2076.4d0, .84151d0, 1666.6d0, .0408010d0,
     &       6.8313d0,  .0040486d0, 4.3987d0, .51675d0,
     &       2003.0d0, .83437d0, 1410.4d0, .03478d0,
     &       1.2707d0,  .0040553d0, 4.6295d0, .57428d0,
     &       1626.3d0, .81858d0, 995.68d0, .0553190d0,
     &       5.7561d0,  .0040491d0, 4.357d0,  .52496d0,
     &      2207.3d0, .83796d0, 1579.5d0, .0271650d0,
     &      14.127d0,   .0040596d0, 5.8304d0, .37755d0,
     &      3645.9d0, .87823d0, 3411.8d0, .0163920d0,
     &       6.6948d0,  .0040603d0, 4.9361d0, .47961d0,
     &       2719.0d0, .85249d0, 1885.8d0, .0197130d0,
     &       3.0619d0,  .0040511d0, 3.5803d0, .59082d0,
     &      2346.1d0, .83713d0, 1222.0d0, .0200720d0/
      DATA   ( (PCOEF(I,J), I= 1 ,8), J=73 ,81)
     &    / 10.811d0,   .0033008d0, 1.3767d0, .76512d0,
     &      2003.7d0, .82269d0, 1110.6d0, .0249580d0,
     &       2.7101d0,  .0040961d0, 1.2289d0, .98598d0,
     &       1232.4d0, .79066d0, 155.42d0, .0472940d0,
     &        .52345d0, .0040244d0, 1.4038d0, .8551d0,
     &        1461.4d0, .79677d0, 503.34d0, .0367890d0,
     &        .4616d0,  .0040203d0, 1.3014d0, .87043d0,
     &        1473.5d0, .79687d0, 443.09d0, .0363010d0,
     &        .97814d0, .0040374d0, 2.0127d0, .7225d0,
     &       1890.8d0, .81747d0, 930.7d0, .02769d0,
     &       3.2086d0,  .004051d0, 3.6658d0,  .53618d0,
     &       3091.2d0, .85602d0, 1508.1d0, .0154010d0,
     &       2.0035d0,  .0040431d0, 7.4882d0, .3561d0,
     &      4464.3d0, .88836d0, 3966.5d0, .0128390d0,
     &      15.43d0,    .0039432d0, 1.1237d0, .70703d0,
     &      4595.7d0, .88437d0, 1576.5d0, .0088534d0,
     &       3.1512d0,  .0040524d0, 4.0996d0, .5425d0,
     &       3246.3d0, .85772d0, 1691.8d0, .0150580d0/
      DATA   ( (PCOEF(I,J), I= 1 ,8), J=82 ,90)
     &     / 7.1896d0,  .0040588d0, 8.6927d0, .35842d0,
     &       4760.6d0, .88833d0, 2888.3d0, .0110290d0,
     &       9.3209d0,  .004054d0, 11.543d0,  .32027d0,
     &      4866.2d0, .89124d0, 3213.4d0, .0119350d0,
     &      29.242d0,   .0036195d0, .16864d0, 1.1226d0,
     &      5688.0d0, .89812d0, 1033.3d0, .0071303d0,
     &       1.8522d0,  .0039973d0, 3.1556d0, .65096d0,
     &       3755.0d0, .86383d0, 1602.0d0, .0120420d0,
     &       3.222d0,   .0040041d0, 5.9024d0, .52678d0,
     &       4040.2d0, .86804d0, 1658.4d0, .0117470d0,
     &       9.3412d0,  .0039661d0, 7.921d0,  .42977d0,
     &      5180.9d0, .88773d0, 2173.2d0, .0092007d0,
     &      36.183d0,   .0036003d0, .58341d0, .86747d0,
     &      6990.2d0, .91082d0, 1417.1d0, .0062187d0,
     &       5.9284d0,  .0039695d0, 6.4082d0, .52122d0,
     &       4619.5d0, .88083d0, 2323.5d0, .0116270d0,
     &       5.2454d0,  .0039744d0, 6.7969d0, .48542d0,
     &      4586.3d0, .87794d0, 2481.5d0, .0112820d0/
      DATA   ( (PCOEF(I,J), I= 1 ,8), J=91 ,97)
     &   /  33.702d0,   .0036901d0, .47257d0, .89235d0,
     &      5295.7d0, .8893d0, 2053.3d0, .0091908d0,
     &       2.7589d0,  .0039806d0, 3.2092d0, .66122d0,
     &     2505.4d0, .82863d0, 2065.1d0, .0228160d0,
     &       2.7589d0,  .0039806d0, 3.2092d0, .66122d0,
     &     2505.4d0, .82863d0, 2065.1d0, .0228160d0,
     &       2.7589d0,  .0039806d0, 3.2092d0, .66122d0,
     &     2505.4d0, .82863d0, 2065.1d0, .0228160d0,
     &       2.7589d0,  .0039806d0, 3.2092d0, .66122d0,
     &     2505.4d0, .82863d0, 2065.1d0, .0228160d0,
     &       2.7589d0,  .0039806d0, 3.2092d0, .66122d0,
     &     2505.4d0, .82863d0, 2065.1d0, .0228160d0,
     &       2.7589d0,  .0039806d0, 3.2092d0, .66122d0,
     &     2505.4d0, .82863d0, 2065.1d0, .0228160d0/
      PE0 = 25.D0
      PE    = DMAX1(PE0,E)
      SL = (PCOEF(1 ,Z)*
     &     (PE**PCOEF(2 ,Z))) + PCOEF(3 ,Z)*(PE**PCOEF(4 ,Z))
      SH = PCOEF(5,Z)/
     &     (PE**PCOEF(6,Z))*DLOG((PCOEF(7,Z)/PE) + PCOEF(8,Z)*PE)
      SP = SL * SH/(SL+SH)
      IF (E .LE. PE0) THEN
         IF (Z .GT. 6) THEN
            VELPWR = 0.45D0
            SP = SP * ((E/PE0)**VELPWR)
         ELSE IF (Z .LE. 6) THEN
            VELPWR = 0.25D0
            SP = SP*((E/PE0)**VELPWR)
         END IF
      END IF
      RPSTOP = SP
      RETURN
      END
      SUBROUTINE SBBLOK ( BLOKS, INTEGS, NBLOKS, IPIVOT, B, X )
      INTEGER INTEGS(3,NBLOKS),IPIVOT(*), I,INDEX,INDEXB,INDEXX,J,LAST,
     1        NBP1,NCOL,NROW
      DOUBLE PRECISION BLOKS(*),B(*),X(*)
      INDEX = 1
      INDEXB = 1
      INDEXX = 1
      DO 10 I=1,NBLOKS
           NROW = INTEGS(1,I)
           LAST = INTEGS(3,I)
           CALL SUBFOR ( BLOKS(INDEX), IPIVOT(INDEXB), NROW, LAST,
     1          B(INDEXB), X(INDEXX) )
           INDEX = NROW*INTEGS(2,I) + INDEX
           INDEXB = INDEXB + NROW
   10 INDEXX = INDEXX + LAST
      NBP1 = NBLOKS + 1
      DO 20 J=1,NBLOKS
           I = NBP1 - J
           NROW = INTEGS(1,I)
           NCOL = INTEGS(2,I)
           LAST = INTEGS(3,I)
           INDEX = INDEX - NROW*NCOL
           INDEXB = INDEXB - NROW
           INDEXX = INDEXX - LAST
   20 CALL SUBBAK ( BLOKS(INDEX), IPIVOT(INDEXB), NROW, NCOL,
     1     LAST, X(INDEXX) )
      RETURN
      END
      SUBROUTINE SCOEF(Z,ATRHO,VFERMI,LFCTR,POT)
      implicit real*8 (a-h,o-z)
     
      common /ionpot/ rh2o, ih2o
!$OMP THREADPRIVATE(/ionpot/)
     
      INTEGER*4  Z
      INTEGER*4  I,J
      REAL*8     ATRHO,VFERMI,LFCTR,POT
      REAL*8     TABLE(4,97)
       DATA ( (TABLE(I,J), I=1,4), J=1,18)
     &         /4.271d0, 1.0309d0,  1.00d0, 19.2d0,
     &          1.894d0, 0.15976d0, 1.00d0, 41.8d0,
     &          4.597d0, 0.59782d0, 1.10d0, 40.0d0,
     &         12.046d0, 1.0781d0,  1.06d0, 63.7d0,
     &         13.093d0, 1.0486d0,  1.01d0, 76.d0,
     &         11.364d0, 1.00d0,    1.03d0, 78.d0,
     &          3.481d0, 1.058d0,   1.04d0, 82.d0,
     &          4.302d0, 0.93942d0, 0.99d0, 95.0d0,
     &          3.522d0, 0.74562d0, 0.95d0, 115.d0,
     &          3.585d0, 0.3424d0,  0.90d0, 137.d0,
     &          2.541d0, 0.45259d0, 0.82d0, 149.d0,
     &          4.302d0, 0.71074d0, 0.81d0, 156.d0,
     &          6.023d0, 0.90519d0, 0.83d0, 166.d0,
     &          4.977d0, 0.97411d0, 0.88d0, 173.d0,
     &          3.542d0, 0.97184d0, 1.00d0, 173.d0,
     &          3.885d0, 0.89852d0, 0.95d0, 180.d0,
     &          3.22d0,  0.70827d0, 0.97d0, 174.d0,
     &          2.488d0, 0.39816d0, 0.99d0, 188.d0/
       DATA( (TABLE(I,J), I=1,4), J=19,36)
     &         /1.329d0, 0.36552d0, 0.98d0, 190.d0,
     &          2.014d0, 0.62712d0, 0.97d0, 191.d0,
     &          4.015d0, 0.81707d0, 0.98d0, 191.d0,
     &          5.682d0, 0.9943d0,  0.97d0, 233.d0,
     &          7.213d0, 1.1423d0,  0.96d0, 245.d0,
     &          8.33d0,  1.2381d0,  0.93d0, 257.d0,
     &          8.15d0,  1.1222d0,  0.91d0, 272.d0,
     &          8.483d0, 0.92705d0, 0.9d0,  286.d0,
     &          8.989d0, 1.0047d0,  0.88d0, 297.d0,
     &          9.125d0, 1.2d0,     0.9d0,  311.d0,
     &          8.483d0, 1.0661d0,  0.9d0,  322.d0,
     &          6.546d0, 0.97411d0, 0.9d0,  330.d0,
     &          5.104d0, 0.84912d0, 0.9d0,  334.d0,
     &          4.428d0, 0.95d0,    0.85d0, 350.d0,
     &          4.597d0, 1.0903d0,  0.9d0,  347.d0,
     &          3.65d0,  1.0429d0,  0.9d0,  348.d0,
     &          2.562d0, 0.49715d0, 0.91d0, 343.d0,
     &          1.87d0,  0.37755d0, 0.92d0, 352.d0/
       DATA ( (TABLE(I,J), I=1,4), J=37,54)
     &         /1.077d0, 0.35211d0, 0.9d0,  363.d0,
     &          1.787d0, 0.57801d0, 0.9d0,  366.d0,
     &          3.041d0, 0.77773d0, 0.9d0,  379.d0,
     &          4.271d0, 1.0207d0,  0.9d0,  393.d0,
     &          5.576d0, 1.029d0,   0.9d0,  417.d0,
     &          6.407d0, 1.2542d0,  0.88d0, 424.d0,
     &          0.0d0,   1.122d0,   0.9d0,  428.d0,
     &          7.256d0, 1.1241d0,  0.88d0, 441.d0,
     &          7.256d0, 1.0882d0,  0.88d0, 449.d0,
     &          6.767d0, 1.2709d0,  0.9d0,  470.d0,
     &          5.847d0, 1.2542d0,  0.9d0,  470.d0,
     &          4.597d0, 0.90094d0, 0.88d0, 469.d0,
     &          3.836d0, 0.74093d0, 0.9d0,  488.d0,
     &          3.695d0, 0.86054d0, 0.9d0,  488.d0,
     &          3.273d0, 0.93155d0, 0.9d0,  487.d0,
     &          2.938d0, 1.0047d0,  0.9d0,  485.d0,
     &          2.343d0, 0.55379d0, 0.96d0, 474.d0,
     &          1.403d0, 0.43289d0, 1.2d0,  482.d0/
       DATA ( (TABLE(I,J), I=1,4), J=55,72)
     &          /.860d0, 0.32636d0, 0.9d0,  488.d0,
     &          1.544d0, 0.5131d0,  0.88d0, 491.d0,
     &          2.676d0, 0.6950d0,  0.88d0, 501.d0,
     &          2.868d0, 0.72591d0, 0.85d0, 523.d0,
     &          2.895d0, 0.71202d0, 0.90d0, 535.d0,
     &          2.923d0, 0.67413d0, 0.90d0, 546.d0,
     &          0.00d0,  0.71418d0, 0.92d0, 560.d0,
     &          3.026d0, 0.71453d0, 0.95d0, 574.d0,
     &          2.084d0, 0.5911d0,  0.99d0, 580.d0,
     &          3.026d0, 0.70263d0, 1.03d0, 591.d0,
     &          3.136d0, 0.68049d0, 1.05d0, 614.d0,
     &          3.170d0, 0.68203d0, 1.07d0, 628.d0,
     &          3.220d0, 0.68121d0, 1.08d0, 650.d0,
     &          3.273d0, 0.68532d0, 1.10d0, 658.d0,
     &          3.327d0, 0.68715d0, 1.08d0, 674.d0,
     &          2.428d0, 0.61884d0, 1.08d0, 684.d0,
     &          3.383d0, 0.71801d0, 1.08d0, 694.d0,
     &          4.428d0, 0.83048d0, 1.08d0, 705.d0/
       DATA ( (TABLE(I,J), I=1,4), J=73,90)
     &         /5.525d0, 1.1222d0,  1.09d0, 718.d0,
     &          6.320d0, 1.2381d0,  1.09d0, 727.d0,
     &          6.805d0, 1.045d0,   1.10d0, 736.d0,
     &          7.144d0, 1.0733d0,  1.11d0, 746.d0,
     &          7.052d0, 1.0953d0,  1.12d0, 757.d0,
     &          6.618d0, 1.2381d0,  1.13d0, 790.d0,
     &          5.904d0, 1.2879d0,  1.14d0, 790.d0,
     &          4.069d0, 0.78654d0, 1.15d0, 800.d0,
     &          3.501d0, 0.66401d0, 1.17d0, 810.d0,
     &          3.291d0, 0.84912d0, 1.20d0, 823.d0,
     &          2.827d0, 0.88433d0, 1.18d0, 823.d0,
     &          2.653d0, 0.80746d0, 1.17d0, 830.d0,
     &          0.00d0,  0.43357d0, 1.17d0, 825.d0,
     &          0.00d0,  0.41923d0, 1.16d0, 794.d0,
     &          0.00d0,  0.43638d0, 1.16d0, 827.d0,
     &          1.338d0, 0.51464d0, 1.16d0, 826.d0,
     &          0.00d0,  0.73087d0, 1.16d0, 841.d0,
     &          3.026d0, 0.81065d0, 1.16d0, 847.d0/
       DATA ( (TABLE(I,J), I=1,4), J=91,92)
     &         /4.015d0, 1.9578d0, 1.16d0, 878.d0,
     &          4.818d0, 1.0257d0, 1.16d0, 890.d0/
       DATA ( (TABLE(I,J), I=1,4), J=93,97)
     &         /0.000d0, 1.0257d0, 1.16d0, 902.d0,
     &          4.818d0, 0.41923d0, 1.16d0, 890.d0,
     &          0.000d0, 1.0257d0, 1.16d0, 934.d0,
     &          0.000d0, 1.0257d0, 1.16d0, 939.d0,
     &          0.000d0, 1.0257d0, 1.16d0, 952.d0/
      ATRHO   = TABLE(1,min(97,Z))*1.D22
      VFERMI  = TABLE(2,min(97,Z))
      LFCTR   = TABLE(3,min(97,Z))
      POT     = TABLE(4,min(97,Z))
     
     
         if( ih2o .gt. 0 ) pot = pot * rh2o
     
      RETURN
      END
      SUBROUTINE SEZI (ZP,ZT,MP,MT,ENERGY,RHO,SEOUT,POT,ZETA)
      implicit real*8 (a-h,o-z)
      INTEGER*4 ZP,ZT
      common /fixcharge/ifixchg     !T.Sato 2019/02/17
      REAL*8 MP,MT,ENERGY,RHO,SEOUT,POT,E,DUMMY1,DUMMY2,DUMMY3
      REAL*8 LFCTR,ATRHO,VFERMI,HE0,SP,HELP1,HELP2,HELP3,HELP4
      REAL*8 HE,B,A,HEH,SE,YRMIN,V,VR,HELP,YR,Q,L,L0,Q1,L1,ZETA
      REAL*8 VRMIN,VMIN,EEE,POWER,DABS,DEXP,DMAX1,DMIN1,HELP5
      REAL*8 RPSTOP
      SAVE LFCTR,ATRHO,VFERMI
!$OMP THREADPRIVATE(LFCTR,ATRHO,VFERMI)
      CALL SCOEF(ZP,DUMMY1,DUMMY2,lfctr,DUMMY3)
      CALL SCOEF(ZT,ATRHO,VFERMI,dummy1,pot)
      E = ENERGY / MP
      IF (ZP .EQ. 1) THEN
        ZETA = DBLE(ZP)
        SEOUT = RPSTOP(ZT,E) * 0.60222D0 / MT
      ELSE IF (ZP .EQ. 2) THEN
        HE0 = 1.0D0
        HE = DMAX1(HE0, E)
        B = DLOG (HE)
        A = .2865D0 + .1266D0*B - .001429D0*(B**2) +&
          0.02402D0*B**3 - 0.01135D0*(B**4) + 0.001475D0*B**5
        HEH =1.0D0 - DEXP (- DMIN1 (30.0D0, A))
        A = (1.0D0 + (0.007D0 + 0.00005D0 * DBLE(ZT)) *&
             DEXP(-(7.6D0 - DMAX1(0.0D0 , DLOG(HE)))**2))
        HEH = HEH *(A**2)
        ZETA = DSQRT(HEH)
        SE = RPSTOP(ZT,E)*HEH*4.0D0
        IF (E.LE.HE0) THEN
          SE = SE * DSQRT(E/HE0)
        END IF
        SEOUT = SE * 0.60222D0/MT
      ELSE IF (ZP .GT. 2) THEN
        YRMIN = 0.130D0
        VRMIN = 1.0D0
        V = DSQRT (E/25.0D0)/VFERMI
        IF (V .GE. 1.0D0) THEN
          VR = V * VFERMI * (1.D0 + 1.D0 / (5.D0 * (V**2)))
        ELSE
          VR = (3.D0*VFERMI/4.D0)*&
         (1.D0 + (2.D0/3.D0*V**2) - (V**4)/15.D0)
        END IF
        HELP = VR / (DBLE(ZP)**0.6667D0)
        YR = DMAX1 (YRMIN,HELP)
        HELP = VRMIN/(DBLE(ZP)**0.6667D0)
        YR = DMAX1(YR,HELP)
        A = -0.803D0*(YR**0.3D0) + 1.3167D0*(YR**0.6D0) +&
              0.38157D0*YR + 0.008983D0*(YR**2)
        Q = DMIN1(1.0D0, DMAX1(0.0D0 , (1.0D0 -&
             DEXP (-DMIN1(A, 50.0D0)))))
        if(ifixchg.ge.1) Q=1.0d0 ! fixed charge mode, T.Sato 2019/02/17
        B = (DMIN1 (0.43D0,DMAX1 (0.32D0,0.12D0 + 0.025D0*DBLE(ZP))))&
              / (DBLE(ZP)**0.3333D0)
        L0 = (.8D0 - Q *&
       DMIN1 (1.2D0,0.6D0 + DBLE(ZP)/30.0D0))/(DBLE(ZP)**0.3333D0)
        IF (Q .LT. 0.2D0) THEN
          L1 = 0.0D0
        ELSE IF (Q .LT. (DMAX1(0.0D0,0.9D0 - 0.025D0*DBLE(ZP)))) THEN
          Q1 = 0.2D0
          L1 = B*(Q - 0.2D0)/DABS(DMAX1(0.D0,.9D0 - .025*DBLE(ZP))&
               - 0.2000001D0)
         ELSE IF (Q .LT.&
          (DMAX1(0.0D0,1.0D0 - 0.025D0*DMIN1(16.0D0,DBLE(ZP))))) THEN
           L1 = B
         ELSE
           L1 = B*(1.0D0 - Q)/(0.025D0*DMIN1(16.0D0 , DBLE(ZP)))
         END IF
         L = DMAX1(L1,L0*LFCTR)
         ZETA = Q + (1.0D0/(2.0D0*(VFERMI**2)))*(1.0D0 - Q)&
             * DLOG(1.0D0 + (4.0D0*L*VFERMI/1.919D0)**2)
         A = -((7.6D0 - DMAX1(0.0D0, DLOG(E)))**2)
         ZETA = ZETA*(1.D0 + (1.D0/(DBLE(ZP)**2))*&
            (0.18D0 + .0015D0*DBLE(ZT))*DEXP(A))
         IF (YR .LE. (DMAX1(YRMIN, VRMIN/(DBLE(ZP)**0.6667D0)))) THEN
           VRMIN=DMAX1(VRMIN, YRMIN*(DBLE(ZP)**0.6667D0))
           VMIN=.5D0*&
             (VRMIN + DSQRT(DMAX1(0.0D0,VRMIN**2 - .8D0*(VFERMI**2))))
           EEE = 25.0D0*VMIN **2
           POWER = 0.5D0
           IF ((ZT .EQ. 6) .OR. (((ZT .EQ. 14) .OR. (ZT .EQ. 32)) .AND.&
            (ZP .LE. 19))) THEN
             POWER = 0.35D0
           END IF
           SP = RPSTOP(ZT,EEE)
           HELP1 = ZETA * DBLE(ZP)
           HELP2 = HELP1 ** 2
           HELP3 = E / EEE
           HELP4 = HELP3 ** POWER
           HELP5 = HELP2 * HELP4
           SE = SP * HELP5
           SEOUT=SE*0.60222D0/MT
         ELSE
           SEOUT = RPSTOP(ZT,E)*((ZETA*DBLE(ZP))**2)*0.60222D0/MT
         END IF
      ELSE
        STOP
      END IF
      RETURN
      END
      SUBROUTINE SF(X,Z,DMVL)
      RETURN
      END
      SUBROUTINE SHIFTB ( AI, IPIVOT, NROWI, NCOLI, LAST,&
                AI1, NROWI1, NCOLI1 )
      INTEGER IPIVOT(NROWI),LAST, IP,J,JMAX,JMAXP1,M,MMAX
      DOUBLE PRECISION AI(NROWI,NCOLI),AI1(NROWI1,NCOLI1)
      MMAX = NROWI - LAST
      JMAX = NCOLI - LAST
      IF (MMAX .LT. 1 .OR. JMAX .LT. 1)             RETURN
      DO 10 M=1,MMAX
           IP = IPIVOT(LAST+M)
           DO 10 J=1,JMAX
   10 AI1(M,J) = AI(IP,LAST+J)
      IF (JMAX .EQ. NCOLI1)                         RETURN
      JMAXP1 = JMAX + 1
      DO 20 J=JMAXP1,NCOLI1
           DO 20 M=1,MMAX
   20 AI1(M,J) = 0.D0
      RETURN
      END
       DOUBLE PRECISION FUNCTION SIRALE (Z,M)
       implicit real*8 (a-h,o-z)
       include 'atimacnt.inc'
      INTEGER*4 Z
      REAL*8    M
      REAL*8    LR
      IF (Z .EQ. 1) THEN
        LR=61280.D0
      ELSE IF (Z .EQ. 2) THEN
        LR=94000.D0
      ELSE IF (Z .EQ. 3) THEN
        LR=82760.D0
      ELSE IF (Z .EQ. 4) THEN
        LR=65190.D0
      ELSE IF (Z .EQ. 6) THEN
        LR=42700.D0
      ELSE IF (Z .EQ. 7) THEN
        LR=37990.D0
      ELSE IF (Z .EQ. 8) THEN
        LR=34240.D0
      ELSE IF (Z .EQ. 10) THEN
        LR=28940.D0
      ELSE IF (Z .EQ. 13) THEN
        LR=24010.D0
      ELSE IF (Z .EQ. 14) THEN
        LR=21820.D0
      ELSE IF (Z .EQ. 18) THEN
        LR=19550.D0
      ELSE IF (Z .EQ. 26) THEN
        LR=13840.D0
      ELSE IF (Z .EQ. 29) THEN
        LR=12860.D0
      ELSE IF (Z .EQ. 50) THEN
        LR=8820.D0
      ELSE IF (Z .EQ. 54) THEN
        LR=8480.D0
      ELSE IF (Z .EQ. 74) THEN
        LR=6760.D0
      ELSE IF (Z .EQ. 82) THEN
        LR=6370.D0
      ELSE IF (Z .EQ. 92) THEN
        LR=6000.D0
      ELSE
        LR = 1000.D0/((DLOG(184.15D0/DBLE(Z)**(1.D0/3.D0))&
          +DLOG(1194.D0/DBLE(Z)**(2.D0/3.D0))/DBLE(Z)&
          -1.202D0*ALPHA**2*DBLE(Z)**2+1.0369D0*ALPHA**4*DBLE(Z)**47
          -1.008D0*ALPHA**6*DBLE(Z)**6/(1.D0+ALPHA**2*DBLE(Z)**2)))&
          *716.405D0*M/DBLE(Z)**2
      END IF
      SIRALE = LR
      RETURN
      END
      SUBROUTINE SUBBAK ( W, IPIVOT, NROW, NCOL, LAST, X )
      INTEGER IPIVOT(NROW),LAST,  IP,J,K,KP1
      DOUBLE PRECISION W(NROW,NCOL),X(NCOL), SUM
      K = LAST
      IP = IPIVOT(K)
      SUM = 0.D0
      IF (K .EQ. NCOL)                              GO TO 30
      KP1 = K+1
   10 DO 20 J=KP1,NCOL
   20 SUM = W(IP,J)*X(J) + SUM
   30 X(K) = (X(K) - SUM)/W(IP,K)
      IF (K .EQ. 1)                                 RETURN
      KP1 = K
      K = K-1
      IP = IPIVOT(K)
      SUM = 0.D0
      GO TO 10
      END
      SUBROUTINE SUBFOR ( W, IPIVOT, NROW, LAST, B, X )
      INTEGER IPIVOT(NROW), IP,JMAX,K
      DOUBLE PRECISION W(NROW,LAST),B(*),X(NROW),SUM
      IP = IPIVOT(1)
      X(1) = B(IP)
      IF (NROW .EQ. 1)                              GO TO 40
      DO 20 K=2,NROW
           IP = IPIVOT(K)
           JMAX = MIN0(K-1,LAST)
           SUM = 0.D0
           DO 10 J=1,JMAX
   10      SUM = W(IP,J)*X(J) + SUM
   20 X(K) = B(IP) - SUM
      NROWML = NROW - LAST
      IF (NROWML .EQ. 0)                            GO TO 40
      LASTP1 = LAST+1
      DO 30 K=LASTP1,NROW
   30 B(NROWML+K) = X(K)
   40 RETURN
      END
      SUBROUTINE XERABT(MESSG,NMESSG)
      CHARACTER*(*) MESSG
      STOP
      END
      SUBROUTINE XERCTL(MESSG1,NMESSG,NERR,LEVEL,KONTRL)
      CHARACTER*20 MESSG1
      RETURN
      END
      SUBROUTINE XERPRT(MESSG,NMESSG)
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
      CALL XGETUA(LUN,NUNIT)
      LENMES = LEN(MESSG)
      DO 20 KUNIT=1,NUNIT
         IUNIT = LUN(KUNIT)
         IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
         DO 10 ICHAR=1,LENMES,72
            LAST = MIN0(ICHAR+71 , LENMES)
            WRITE (IUNIT,'(1X,A)') MESSG(ICHAR:LAST)
   10    CONTINUE
   20 CONTINUE
      RETURN
      END
      SUBROUTINE XERROR(MESSG,NMESSG,NERR,LEVEL)
      CHARACTER*(*) MESSG
      CALL XERRWV(MESSG,NMESSG,NERR,LEVEL,0,0,0,0,0.,0.)
      RETURN
      END
      SUBROUTINE XERRWV(MESSG,NMESSG,NERR,LEVEL,NI,I1,I2,NR,R1,R2)
      CHARACTER*(*) MESSG
      CHARACTER*20 LFIRST
      CHARACTER*37 FORM
      DIMENSION LUN(5)
      LKNTRL = J4SAVE(2,0,.FALSE.)
      MAXMES = J4SAVE(4,0,.FALSE.)
      IF ((NMESSG.GT.0).AND.(NERR.NE.0).AND.
     1    (LEVEL.GE.(-1)).AND.(LEVEL.LE.2)) GO TO 10
         IF (LKNTRL.GT.0) CALL XERPRT('FATAL ERROR IN...',17)
         CALL XERPRT('XERROR -- INVALID INPUT',23)
         IF (LKNTRL.GT.0) CALL FDUMP
         IF (LKNTRL.GT.0) CALL XERPRT('JOB ABORT DUE TO FATAL ERROR.',
     1  29)
         IF (LKNTRL.GT.0) CALL XERSAV(' ',0,0,0,KDUMMY)
         CALL XERABT('XERROR -- INVALID INPUT',23)
         RETURN
   10 CONTINUE
      JUNK = J4SAVE(1,NERR,.TRUE.)
      CALL XERSAV(MESSG,NMESSG,NERR,LEVEL,KOUNT)
      LFIRST = MESSG
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      CALL XERCTL(LFIRST,LMESSG,LERR,LLEVEL,LKNTRL)
      LMESSG = NMESSG
      LERR = NERR
      LLEVEL = LEVEL
      LKNTRL = MAX0(-2,MIN0(2,LKNTRL))
      MKNTRL = IABS(LKNTRL)
      IF ((LLEVEL.LT.2).AND.(LKNTRL.EQ.0)) GO TO 100
      IF (((LLEVEL.EQ.(-1)).AND.(KOUNT.GT.MIN0(1,MAXMES)))
     1.OR.((LLEVEL.EQ.0)   .AND.(KOUNT.GT.MAXMES))
     2.OR.((LLEVEL.EQ.1)   .AND.(KOUNT.GT.MAXMES).AND.(MKNTRL.EQ.1))
     3.OR.((LLEVEL.EQ.2)   .AND.(KOUNT.GT.MAX0(1,MAXMES)))) GO TO 100
         IF (LKNTRL.LE.0) GO TO 20
            CALL XERPRT(' ',1)
            IF (LLEVEL.EQ.(-1)) CALL XERPRT
     1('WARNING MESSAGE...THIS MESSAGE WILL ONLY BE PRINTED ONCE.',57)
            IF (LLEVEL.EQ.0) CALL XERPRT('WARNING IN...',13)
            IF (LLEVEL.EQ.1) CALL XERPRT
     1      ('RECOVERABLE ERROR IN...',23)
            IF (LLEVEL.EQ.2) CALL XERPRT('FATAL ERROR IN...',17)
   20    CONTINUE
         CALL XERPRT(MESSG,LMESSG)
         CALL XGETUA(LUN,NUNIT)
         ISIZEI = LOG10(FLOAT(I1MACH(9))) + 1.0
         ISIZEF = LOG10(FLOAT(I1MACH(10))**I1MACH(11)) + 1.0
         DO 50 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
            DO 22 I=1,MIN(NI,2)
               WRITE (FORM,21) I,ISIZEI
   21          FORMAT ('(11X,21HIN ABOVE MESSAGE, I',I1,'=,I',I2,')   ')
               IF (I.EQ.1) WRITE (IUNIT,FORM) I1
               IF (I.EQ.2) WRITE (IUNIT,FORM) I2
   22       CONTINUE
            DO 24 I=1,MIN(NR,2)
               WRITE (FORM,23) I,ISIZEF+10,ISIZEF
   23          FORMAT ('(11X,21HIN ABOVE MESSAGE, R',I1,'=,E',
     1         I2,'.',I2,')')
               IF (I.EQ.1) WRITE (IUNIT,FORM) R1
               IF (I.EQ.2) WRITE (IUNIT,FORM) R2
   24       CONTINUE
            IF (LKNTRL.LE.0) GO TO 40
               WRITE (IUNIT,30) LERR
   30          FORMAT (15H ERROR NUMBER =,I10)
   40       CONTINUE
   50    CONTINUE
         IF (LKNTRL.GT.0) CALL FDUMP
  100 CONTINUE
      IFATAL = 0
      IF ((LLEVEL.EQ.2).OR.((LLEVEL.EQ.1).AND.(MKNTRL.EQ.2)))
     1IFATAL = 1
      IF (IFATAL.LE.0) RETURN
      IF ((LKNTRL.LE.0).OR.(KOUNT.GT.MAX0(1,MAXMES))) GO TO 120
         IF (LLEVEL.EQ.1) CALL XERPRT
     1   ('JOB ABORT DUE TO UNRECOVERED ERROR.',35)
         IF (LLEVEL.EQ.2) CALL XERPRT
     1   ('JOB ABORT DUE TO FATAL ERROR.',29)
         CALL XERSAV(' ',-1,0,0,KDUMMY)
  120 CONTINUE
      IF ((LLEVEL.EQ.2).AND.(KOUNT.GT.MAX0(1,MAXMES))) LMESSG = 0
      CALL XERABT(MESSG,LMESSG)
      RETURN
      END
      SUBROUTINE XERSAV(MESSG,NMESSG,NERR,LEVEL,ICOUNT)
      INTEGER LUN(5)
      CHARACTER*(*) MESSG
      CHARACTER*20 MESTAB(10),MES
      DIMENSION NERTAB(10),LEVTAB(10),KOUNT(10)
      SAVE MESTAB,NERTAB,LEVTAB,KOUNT,KOUNTX
!$OMP THREADPRIVATE(MESTAB,NERTAB,LEVTAB,KOUNT,KOUNTX)
      DATA KOUNT(1),KOUNT(2),KOUNT(3),KOUNT(4),KOUNT(5),
     1     KOUNT(6),KOUNT(7),KOUNT(8),KOUNT(9),KOUNT(10)
     2     /0,0,0,0,0,0,0,0,0,0/
      DATA KOUNTX/0/
      IF (NMESSG.GT.0) GO TO 80
         IF (KOUNT(1).EQ.0) RETURN
         CALL XGETUA(LUN,NUNIT)
         DO 60 KUNIT=1,NUNIT
            IUNIT = LUN(KUNIT)
            IF (IUNIT.EQ.0) IUNIT = I1MACH(4)
            WRITE (IUNIT,10)
   10       FORMAT (32H0          ERROR MESSAGE SUMMARY/
     1      51H MESSAGE START             NERR     LEVEL     COUNT)
            DO 20 I=1,10
               IF (KOUNT(I).EQ.0) GO TO 30
               WRITE (IUNIT,15) MESTAB(I),NERTAB(I),LEVTAB(I),KOUNT(I)
   15          FORMAT (1X,A20,3I10)
   20       CONTINUE
   30       CONTINUE
            IF (KOUNTX.NE.0) WRITE (IUNIT,40) KOUNTX
   40       FORMAT (41H0OTHER ERRORS NOT INDIVIDUALLY TABULATED=,I10)
            WRITE (IUNIT,50)
   50       FORMAT (1X)
   60    CONTINUE
         IF (NMESSG.LT.0) RETURN
         DO 70 I=1,10
   70       KOUNT(I) = 0
         KOUNTX = 0
         RETURN
   80 CONTINUE
      MES = MESSG
      DO 90 I=1,10
         II = I
         IF (KOUNT(I).EQ.0) GO TO 110
         IF (MES.NE.MESTAB(I)) GO TO 90
         IF (NERR.NE.NERTAB(I)) GO TO 90
         IF (LEVEL.NE.LEVTAB(I)) GO TO 90
         GO TO 100
   90 CONTINUE
         KOUNTX = KOUNTX+1
         ICOUNT = 1
         RETURN
  100    KOUNT(II) = KOUNT(II) + 1
         ICOUNT = KOUNT(II)
         RETURN
  110    MESTAB(II) = MES
         NERTAB(II) = NERR
         LEVTAB(II) = LEVEL
         KOUNT(II)  = 1
         ICOUNT = 1
         RETURN
      END
      SUBROUTINE XGETUA(IUNITA,N)
      DIMENSION IUNITA(5)
      N = J4SAVE(5,0,.FALSE.)
      DO 30 I=1,N
         INDEX = I+4
         IF (I.EQ.1) INDEX = 3
         IUNITA(I) = J4SAVE(INDEX,0,.FALSE.)
   30 CONTINUE
      RETURN
      END
