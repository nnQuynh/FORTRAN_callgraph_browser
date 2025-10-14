!  SCINFUL related subroutines used in PHITS
!  Last update on February 25, 2025.
************************************************************************
*                                                                      *
      subroutine scinmode(Eneut) ! called from nreac.f
*                                                                      *
*     Eneut: neutron energy in MeV                                     *
*     determine secondary particle information using SCINFUL           *
************************************************************************
      implicit real*8 (a-h,o-z)

! define direction, should be always forward
      U = 0.0d0
      V = 0.0d0
      W = 1.0d0

      call CX_prepr(Eneut, U, V, W)
      return
      end

C =====================================================================
C     LIST OF ROUTIENS IN THIS FILE
C ---------------------------------------------------------------------
C     S: CX_initialize
C     S: CX_prepr
C     S: CX_score
C     S: CX_anal
C =====================================================================
      subroutine CX_init
C ---------------------------------------------------------------------
      common /LGNDRE/ IC, E(220), F(6,220)
      common /DATA_PSCAT/ A1(220), B1(220), B2(220), B3(220),
     +                    B4(220), B5(220), B6(220)
C ---------------------------------------------------------------------
C     SUBROUTINE PSCAT
C     Set-up needed only once -- the first time through.
      Do 2 I=1,Ic
      A1(I)=E(I)
      B1(I)=F(1,I)
      B2(I)=F(2,I)
      B3(I)=F(3,I)
      B4(I)=F(4,I)
      B5(I)=F(5,I)
    2 B6(I)=F(6,I)
C ---------------------------------------------------------------------
      return
      end subroutine CX_init
C =====================================================================
C =====================================================================
      subroutine CX_prepr(Eneut, U, V, W)
C ---------------------------------------------------------------------
      real*8 totalx, Eneut ! T.Sato 2020/02/16
      real*8 U, V, W
C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------
      common /NEUTRN/ En, Un, Vn, Wn
!$OMP THREADPRIVATE(/NEUTRN/)
      Common /NEUTR2/ Eneut2, U2,V2,W2
!$OMP THREADPRIVATE(/NEUTR2/)
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
!$OMP THREADPRIVATE(/VECTOR/)
      Common /CSCATT/ Xp,Yp,Zp, Vnca,Vcca,Vcom
!$OMP THREADPRIVATE(/CSCATT/)
      Common /P3ALF/ Pbar(13)
!$OMP THREADPRIVATE(/P3ALF/)
      common /PROB/ P(12), indx
!$OMP THREADPRIVATE(/PROB/)
C ---------------------------------------------------------------------
C     Copy the data to /NEUTRN/
      En = sngl(Eneut)
      Un = sngl(U)
      Vn = sngl(V)
      Wn = sngl(W)
C ---------------------------------------------------------------------
C     Initialization for /VECTOR/
      Xpn = 0.d0
      Ypn = 0.d0
      Zpn = 0.d0
      Xn  = Un
      Yn  = Vn
      Zn  = Wn
C ---------------------------------------------------------------------
C     Initialization for /NEUTR2/
      Eneut2=0
      U2=0
      V2=0
      W2=0
C ---------------------------------------------------------------------
C     Initialization for /CSCATT/
      Xpn = 0.d0
      Ypn = 0.d0
      Zpn = 0.d0
      Xn  = 0.d0
      Yn  = 0.d0
      Zn  = 0.d0
C ---------------------------------------------------------------------
C     Initialization for /P3ALF/
      Pbar(:) = 0.d0
C ---------------------------------------------------------------------
C     Initialization for /PROB/
      indx = 0
      P(:) = 0.d0
C ---------------------------------------------------------------------
      KF(:)   = 0
      DENG(:) = 0.d0
      DVX(:)  = 0.d0
      DVY(:)  = 0.d0
      DVZ(:)  = 0.d0
      Npart   = 0
C ---------------------------------------------------------------------
      totx = TOTALX(Eneut)
      itype= 0
      call IBOX(itype)
C ---------------------------------------------------------------------
      k=itype

      select case(k)
      case(2)
        call pscat
      case(3)
        call inelas
      case(4)
        call nalpha
      case(5)
        call nn3alf
      case(6)
        call n3he
      case(7)
        call npx
      case(8)
        call n2n
      case(9)
        call nd
      case(10)
        call nt
      case default
      end select
C ---------------------------------------------------------------------
      return
      end subroutine CX_prepr
C =====================================================================
C =====================================================================
      subroutine CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------
      Npart = Npart + 1
C ---------------------------------------------------------------------
      KF(Npart)   = IKF
      DENG(Npart) = DDENG
      DVX(Npart)  = DDVX
      DVY(Npart)  = DDVY
      DVZ(Npart)  = DDVZ
C ---------------------------------------------------------------------
      return
      end subroutine CX_score
C =====================================================================
      subroutine CX_anal
!     output outgoing particle information into iclust, jclust, qclust
*---- in common -------------------------------------------------------*
*                                                                      *
*        nclst   : total number of out going particles and nuclei      *
*        iclust(nclst)                                                 *
*        jclust(i,nclst)                                               *
*        qclust(i,nclst)                                               *
*        numpat(i) : total number of out going particles or nuclei     *
*        mathz, mathn : z and n of mather nucleus                      *
C =====================================================================
      use MMBANKMOD ! need to use wt
      use NGSDATAMOD, only : bindeg

      implicit real*8 (a-h,o-z)

      include 'err.inc'
      include 'param00.inc'  ! need to use nnn
      include 'param-physcnst.inc'

C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /clustp/ rumpat(0:20), numpat(0:20)
!$OMP THREADPRIVATE(/clustp/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

      if( Npart <= 0 ) then
       write(ErrCha,'(''ERROR in CX_anal: Npart = ,'',i10)') Npart
       ErrID = 'L:213/R:CX_anal/F:scinucl.f' !E00_000_000
       call ErrWrite(ErrID,ErrCha)
       call parastop(121)
      end if
C ---------------------------------------------------------------------
!  following part are taken from subroutine nfrgdat in ncasc.f, T.Sato 2020/02/16
      mathz=6  ! mother is always 12C
      mathn=6  ! mother is always 12C
      wgt=1.0d0 ! weight ratio between primary and secondary particles
      nclst=0  ! reset nclst

      do i = 1, npart
       eout=deng(i) ! energy of outgoing particle in MeV

       if( kf(i) .eq. 2212 ) then ! proton
        ibary = 1
        ipid  = 1
        ippad = 1
        ipprt = 1
        ipneu = 0
        ipchg = 1
        bene = 0.0d0
        rms = rstms(1)
       else if( kf(i) .eq. 2112 ) then ! neutron
        ibary = 1
        ipid  = 2
        ippad = 2
        ipprt = 0
        ipneu = 1
        ipchg = 0
        bene = 0.0d0
        rms = rstms(2)
       else if( kf(i) .eq. 22 ) then ! photon
        ibary = 0
        ipid  = 4
        ippad = 14
        ipprt = 0
        ipneu = 0
        ipchg = 0
        bene = 0d0
        rms = rstms(14)
       else if( kf(i) .eq. 1000002 ) then ! deuteron
        ibary = 2
        ipid  = 0
        ippad = 15
        ipprt = 1
        ipneu = 1
        ipchg = 1
        bene = bindeg(ipprt,ipneu)
        rms = rstms(15)
       else if( kf(i) .eq. 1000003 ) then ! triton
        ibary = 3
        ipid  = 0
        ippad = 16
        ipprt = 1
        ipneu = 2
        ipchg = 1
        bene = bindeg(ipprt,ipneu)
        rms = rstms(16)
       else if( kf(i) .eq. 2000003 ) then ! 3He
        ibary = 3
        ipid  = 0
        ippad = 17
        ipprt = 2
        ipneu = 1
        ipchg = 2
        bene = bindeg(ipprt,ipneu)
        rms = rstms(17)
       else if( kf(i) .eq. 2000004 ) then ! alpha
        ibary = 4
        ipid  = 0
        ippad = 18
        ipprt = 2
        ipneu = 2
        ipchg = 2
        bene = bindeg(ipprt,ipneu)
        rms = rstms(18)
       else
        ipid  = 0
        ibary = mod(kf(i),1000000)
        ipprt = ( kf(i) - ibary ) / 1000000 + 0.1
        ipneu = ibary - ipprt
        ippad = 19
        ipchg = ipprt
        bene = bindeg(ipprt,ipneu)
        rms = 0.93827d0 * ipprt + 0.93958d0 * ipneu
     &                - bene / 1000.0d0
       end if

       tmas = rms / 1000.d0
       etot = tmas + eout / 1000.d0
       psqr = dsqrt( etot**2 - tmas**2 )
       pclx = psqr * dvx(i)
       pcly = psqr * dvy(i)
       pclz = psqr * dvz(i)

       nclst = nclst + 1
       ii    = nclst

       iclust(ii)    = ipid

       jclust(0,ii)  = 0
       jclust(1,ii)  = ipprt
       jclust(2,ii)  = ipneu
       jclust(3,ii)  = ippad
       jclust(4,ii)  = 0
       jclust(5,ii)  = ipchg
       jclust(6,ii)  = ibary
       jclust(7,ii)  = kf(i)
       jclust(8,ii)  = 0

       qclust(0,ii)  = 1.0
       qclust(1,ii)  = pclx
       qclust(2,ii)  = pcly
       qclust(3,ii)  = pclz
       qclust(4,ii)  = etot
       qclust(5,ii)  = tmas
       qclust(6,ii)  = 0.0d0
       qclust(7,ii)  = eout
       qclust(8,ii)  = wgt
       qclust(9,ii)  = 0.0d0
       qclust(10,ii) = 0.0d0
       qclust(11,ii) = 0.0d0
       qclust(12,ii) = 0.0d0

       numpat(ippad) = numpat(ippad) + 1
       rumpat(ippad) = rumpat(ippad) + wgt

      end do

C ---------------------------------------------------------------------
      return
      end subroutine CX_anal
C =====================================================================
C =====================================================================
C     LIST OF ROUTIENS IN THIS FILE
C ---------------------------------------------------------------------
C     F: exterp
C     F: velocity
C     F: efromv
C     S: dircos
C     S: cmtran
C     S: labtran
C     S: rvect
C     S: transvec
C     F: rcke
C =====================================================================
C =====================================================================
C    This is file EXTERP.FOR
C
C      Purpose is to provide interpolation between two points in a table
C      of independent data, A, to get a value from a table of dependent
C      data, B, for an input value C(A).
C          N=dimension of tables A and B
C          Ntype variable:                 FOR:     A-array   B-array
C                                       Ntype =1    Linear    Linear
C                                             =2    Linear    Log
C                                             =3    Log       Linear
C                                             =4    Log       Log
C ---------------------------------------------------------------------
      FUNCTION EXTERP(A,B,C,N,Ntype)
      Real*8 Ac,Ae,Af,Amin,R,Xb,Xc,Y,Z
      Dimension A(1),B(1)
      Data Amin/3.72008D-24/

      Exterp=B(1)
      If (C.LE.A(1)) Return
      Exterp=B(N)
      If (C.GE.A(N)) Return
C
C     First 4 statements take care of case where input value of C was
C     outside of the limits of the tabular values of the A array.
C
      Do 2 I=2,N
      IF (C - A(I)) 3,9,2
    9 Exterp=B(I)
      Return
    2 Continue
    3 J=I
      I=J-1

      L=Ntype
      Ac=C
      Ae=A(I)
      Af=A(J)
      Xb=B(I)
      Xc=B(J)
      If (L.LE.2) goto 4
      IF (Ac .LT. Amin) Ac=Amin
      Ac=Dlog(Ac)
      IF (Ae .LT. Amin) Ae=Amin
      Ae=Dlog(Ae)
      IF (Af .LT. Amin) Af=Amin
      Af=Dlog(Af)
    4 IF (L.LE.1 .OR. L.EQ.3) goto 7
      IF (Xb .LT. Amin) Xb=Amin
      Xb=Dlog(Xb)
      IF (Xc .LT. Amin) Xc=Amin
      Xc=Dlog(Xc)
    7 R=(Ac-Ae)/(Af-Ae)
      Y=Xb+R*(Xc-Xb)
      Z=Y
      IF (L.EQ.2 .OR. L.EQ.4) Z=Dexp(Y)
      Exterp=Z
      Return
      END
C =====================================================================
C =====================================================================
      Function VELOCITY(N,E)
C
C       N=1  Particle=neutron
C         2     "    =proton
C         3     "    =alpha
C         4     "    =9-Be
C         5     "    =11-B
C         6     "    =12-B
C         7     "    =12-C
C         8     "    =8-Be
C         9     "    =11-C
C        10     "    =deuteron
C        11     "    =5-He
C        12     "    =10-B
C        13     "    =triton
C        14     "    =11-Be
C        15     "    =10-Be
C        16     "    =8-Li
C        17     "    =7-Li
C        18     "    =6-Li
C        19     "    =3-He
C        20     "    =pion+ (Added by satoh)
C        21     "    =pion- (Added by satoh)
C
      Real*8 G,Gsq,B,D1
      Common /MASSES/ Emass(21)
      Data D1/1.0D+0/,C/2.997925E+10/
      Data Emass /939.583,938.272,3728.43,8394.86,10255.2,11191.4,
     & 11178.0,7456.95,10257.17,1876.14,4668.9,9327.05,2809.45,
     & 10266.69,9327.62,7472.96,6535.42,5603.097,2809.44,139.570,
     & 139.570/

      V=0.0
      IF (E .LE. 0.0) goto 2
      IF (N.GE.1 .AND. N.LE.20) goto 1
    2 continue
      Goto 7
    1 M=N
      Gamma=E/Emass(M)
      G=Gamma
      G=G+D1
      Gsq=G*G
      B=DSQRT(D1-D1/Gsq)
      Beta=B
      V=C*Beta
    7 Velocity=V
      Return
      END
C =====================================================================
C =====================================================================
      Function EFROMV(N,V)
C     Purpose is to compute particle energy (in MeV) from input velocity V
C     (in cm/sec)
C     N = same as in function VELOCITY

      Real*8 B,G,D1
      Common /MASSES/ Emass(21)
      Data D1/1.0D+0/,C/2.997925E+10/

      E=0.0

      IF (N.GE.1 .AND. N.LE.20) goto 1

      write(*,4)N,V
    4 Format(/'  *** Error in Function EFROMV; N,V = 'I12,1PE10.3/)
      Goto 7
    1 M=N
      Beta=V/C
      B=Beta
      G=D1/DSQRT(D1 - B*B)
      G=G-D1
      Gc=G
      E=Emass(M)*Gc
    7 Efromv=E
      Return
      End
C =====================================================================
C =====================================================================
      SUBROUTINE DIRCOS(V1,V2,V3, C1,C2,C3)
C     Purpose is to return with direction cosines, C1,C2,C3, for input
C     velocities V1,V2,V3.

      Real*8 C,Csq,X,Y,Z
      Vsq=V1*V1 + V2*V2 + V3*V3
      V=SQRT(Vsq)
    1 C1=V1/V
      C2=V2/V
      C3=V3/V
      X=C1
      Y=C2
      Z=C3
      Csq=X*X + Y*Y + Z*Z
      IF (Csq .LE. 1.0D+0) Return
      C=DSQRT(Csq)
      CC=C
      IF (CC .LT. 1.00001) CC=1.00001
      C1=C1/CC
      C2=C2/CC
      C3=C3/CC
      Return
      END
C =====================================================================
C =====================================================================
      Subroutine CMTRAN(N1,N2,V, Vcom,E1,E2)
C
      Real*8 B,Bcom,E,Em1,Em2,D1,Etot,G,Gem1,P1
      Common /MASSES/ Emass(21)
      Data D1/1.0D+0/, C/2.997925E+10/

      Em1=Emass(N1)
      Em2=Emass(N2)
C ---------------------------------------------------------------------
C     Initialize
      Vcom = 0.
      E1   = 0.
      E2   = 0.

C     Let particle of mass Em1 approach particle of mass Em2 with
C     velocity V.  Particle of mass Em2 is at rest.  Then compute
C     velocity of the center of mass of the Em1-Em2 system with respect
C     to the mass-Em2-at-rest system, and compute the relativistic
C     kinetic energies, E1 and E2, of the two particles (respectively)
C     in the center-of-mass system.

      Beta=V/C
      B=Beta
      G=D1/DSQRT(D1-B*B)
      Gem1=G*Em1
      Etot=Gem1+Em2
      P1=Gem1*B
      Bcom=P1/Etot
      Beta=Bcom
      Vcom=Beta*C
      G=D1/DSQRT(D1-Bcom*Bcom)
      E=G*Em2 - Em2
      E2=E
      E=G*(Gem1 - Bcom*P1) - Em1
      E1=E

      Return
      END
C =====================================================================
C =====================================================================
      Subroutine LABTRAN(Vxp,Vyp,Vzp, V, Vx,Vy,Vz)
C     Given velocity components Vxp, Vyp, and Vzp in the "primed"
C     coordinate system which is moving with velocity V along the
C     Z-axis of the laboratory system, compute the relativistically
C     correct velocity components Vx, Vy, and Vz of the velocity
C     in the laboratory system.

      Real*8 Vxx,Vyy,Vzz,Vv,Beta,C,Denom,D1,Sqr
      Data D1/1.0D+0/, C/2.99792458D+10/

      Vxx=Vxp
      Vyy=Vyp
      Vzz=Vzp
      Vv=V
      Beta=Vv/C
      Denom=D1+Beta*Vzz/C
      Sqr=DSQRT(D1-Beta*Beta)
      Vxx=Vxx*Sqr/Denom
      Vyy=Vyy*Sqr/Denom
      Vzz=(Vzz+Vv)/Denom
      Vx=Vxx
      Vy=Vyy
      Vz=Vzz

      Return
      END
C =====================================================================
C =====================================================================
      Subroutine RVECT(Xn,Yn,Zn)
C     Purpose is to obtain a random unit vector with components
C     (i.e. direction cosines) Xn, Yn, Zn

      Real*8 Zet,S,D1
      real*8 UNIRN

      Data D1/1.0D+0/

      Z=1.0-2.0*sngl(UNIRN(dummy))
      ZN=Z
      Zet=Z

      THETA=6.283185*sngl(UNIRN(dummy))

      S=DSQRT(1.0-Zet*Zet)
      SINZ=S

      XN=COS(THETA)*SINZ
      YN=SIN(THETA)*SINZ

      Return
      END
C =====================================================================
C =====================================================================
      Subroutine TRANSVEC(Zx,Zy,Zz)
C
C     Given a vector with direction cosines Vxp,Vyp, and Vzp in a
C     coordinate system having its Z-axis direction cosines as
C     Zx, Zy, and Zz in a second coordinate system, determine the
C     directions cosines of the vector in the second coordinate
C     system. Since the X and Y axes corresponding to the given
C     Z axis are not given, it is assumed that the Y axis lies
C     in the x-y plane of the second coordinate system. The
C     resulting direction cosines of the vector in the second
C     coordinate system are Vx, Vy, and Vz.

      Real*8 D1,Dum1,Dum2,dcut
      Common /VECTOR/ Vxp,Vyp,Vzp, Vx,Vy,Vz
!$OMP THREADPRIVATE(/VECTOR/)
      Data D1/1.0D+0/
      data dcut/-1.0d-2/

      IF (Zz.GE.1.0) goto 2
      Dum1=Zz
      Dum2=DSQRT(D1-Dum1*Dum1)
      Xz=-Dum2

      IF (Xz .GT. dcut ) goto 2
      Yx=Zy/Xz
      Yy=-Zx/Xz
      Xx=Yy*Zz
      Xy=-Yx*Zz

      Vx=Xx*Vxp + Yx*Vyp + Zx*Vzp
      Vy=Xy*Vxp + Yy*Vyp + Zy*Vzp
      Vz=Xz*Vxp          + Zz*Vzp
      Return

    2 Vx=Vxp
      Vy=Vyp
      Vz=Vzp
      Return
      END
C =====================================================================
C =====================================================================
      Function RCKE(N1,N2,TKE)
C     RCKE = Relativistically Correct Kinetic Energy
C     TKE is the total kinetic energy available in the center of mass
C     of the system of particles having masses EM1 and EM2.  RCKE
C     computes the kinetic energy for the particle having mass EM1.

      Real*8 Em1,Em2,E,Etot,Dhalf
      Common /MASSES/ Emass(21)
      Data Dhalf/0.5E+0/

      Em1=Emass(N1)
      Em2=Emass(N2)
      E=TKE
      Etot=E+Em1+Em2
      E=Dhalf*(Etot*Etot + Em1*Em1 - Em2*Em2)/Etot
      E=E-Em1
      RCKE=E
      Return
      END
C =====================================================================
C =====================================================================
C     LIST OF ROUTIENS IN THIS FILE
C ---------------------------------------------------------------------
C     F: totalx
C     F: ibox
C =====================================================================
      Function TOTALX(Eneut)
C     Purpose is to determine total cross section in unit of barns and
C     also the fractional parts of this total due to the separate partial
C     corss sections for incident neutron of energy Eneut.
      real*8 totalx,Eneut ! T.Sato 2020/02/16
      common /PROB/ P(12), indx
!$OMP THREADPRIVATE(/PROB/)
      dimension S(10)
      En=sngl(Eneut)
C ---------------------------------------------------------------------
C     Cross sections
      S(1) = 0.0 ! not use
      S(2) = SIGCELAS(En)
      S(3) = SIGCINEL(En)
      S(4) = SIGCNAL(En)
      S(5) = SIGCNN3A(En)
      S(6) = SIGCN3HE(En)
      S(7) = (SIGCNPN(En)+SIGCNP(En))
      S(8) = SIGCN2N(En)
      S(9) = SIGCND(En)
      S(10)= SIGCNT(En)
C ---------------------------------------------------------------------
      Sigt = S(1)
      do J = 2, 10
        P(J+2) = 1.0
        Sigt = Sigt+S(J)
      end do
C     Sigt = total cross section of n + 12-C in barns

      P(1) = 2.0001
      P(2) = Sigt
      P(3) = 0.0 ! not use

      do J = 2, 10
        P(J+2)=P(J+1)+S(J)/Sigt
        if(P(J+2) .LE. 0.999999) P(1)=P(1)+1.0
      end do

      TOTALX = dble(Sigt)
      Return
      END
C =====================================================================
C =====================================================================
      Subroutine IBOX(itype)
C     Purpose is to determine the event type by random number.

      common /PROB/ P(12), indx
!$OMP THREADPRIVATE(/PROB/)
      real*8 UNIRN

      V=sngl(UNIRN(dummy))

      Ip1=IFIX(P(1)) -1
      K=2

      do I = 1, Ip1
        if( V <= P(I+2) ) exit
        K=K+1
      end do

      itype = K-1

      Return
      END
C =====================================================================
C =====================================================================
C     LIST OF ROUTIENS IN THIS FILE
C ---------------------------------------------------------------------
C     F: sigcelas
C     F: sigcinel
C     F: sigcn2n
C     F: sigcn3he
C     F: sigcnal
C     F: sigcnd
C     F: sigcnn3a
C     F: sigcnp
C     F: sigcnpn
C     F: sigcnt
C =====================================================================
C =====================================================================
      FUNCTION SIGCELAS(En)
C     Purpose is to return a cross section (in barns) for elastic scattering
C     of Carbon by an incident neutron of energy En (in MeV).  These data are
C     taken from the ENDF/B-V Evaluation for En up to about 16 MeV.
C     Other data for En > 16 MeV were obtained by subtracting total
C     non-elastic cross sections from total cross sections.

      Dimension E(541),Sig(541)
      Data NN/541/
      Data (E(I),I=1,173)/
     W  2.530E-08, 1.0E-07, 1.0E-05, 0.001, 0.005, 0.01,
     X  0.015,  0.02,   0.04,   0.05,   0.075,  0.10,   0.125,  0.225,
     X  0.3250, 0.4250, 0.5250, 0.6250, 0.7250, 0.8500, 1.0000, 1.1000,
     Y  1.2000, 1.3000, 1.4000, 1.5000, 1.6000, 1.7000, 1.8000, 1.9000,
     Z  2.0000, 2.0250, 2.0500, 2.0520, 2.0540, 2.0560, 2.0580, 2.0600,
     X  2.0610, 2.0620, 2.0630, 2.0640, 2.0650, 2.0660, 2.0670, 2.0680,
     X  2.0690, 2.0700, 2.0710, 2.0720, 2.0730, 2.0740, 2.0750, 2.0760,
     X  2.0770, 2.0780, 2.0790, 2.0800, 2.0810, 2.0820, 2.0830, 2.0840,
     X  2.0850, 2.0860, 2.0870, 2.0880, 2.0890, 2.0900, 2.0910, 2.0920,
     X  2.0930, 2.0940, 2.0950, 2.0960, 2.0980, 2.1000, 2.1020, 2.1040,
     A  2.1060, 2.1080, 2.1100, 2.1140, 2.1180, 2.1200, 2.1400, 2.1600,
     B  2.1800, 2.2000, 2.2400, 2.2800, 2.3200, 2.36,   2.4000, 2.4400,
     C  2.4800, 2.5200, 2.5600, 2.6000, 2.6400, 2.6600, 2.6800, 2.7000,
     X  2.7200, 2.7400, 2.7600, 2.7800, 2.7900, 2.8000, 2.8040, 2.8060,
     P  2.8080, 2.8100, 2.8110, 2.8120, 2.8130, 2.8140, 2.8150, 2.8160,
     Q  2.8170, 2.8180, 2.8190, 2.8200, 2.8220, 2.8240, 2.8260, 2.8300,
     R  2.835,  2.8400, 2.8500, 2.8600, 2.8650, 2.8700, 2.8750, 2.8800,
     S  2.8850, 2.8900, 2.8950, 2.9000, 2.9050, 2.9100, 2.9150, 2.9200,
     T  2.9280, 2.9320, 2.9360, 2.9400, 2.9440, 2.9480, 2.9520, 2.9560,
     U  2.9600, 2.9640, 2.9680, 2.9720, 2.9760, 2.9800, 2.9840, 2.9880,
     V  2.9920, 2.9960, 3.0000, 3.0100, 3.0200, 3.0300, 3.0400, 3.0600,
     W  3.0800, 3.1000, 3.1200, 3.14,   3.1600, 3.1800, 3.20/
      Data (E(I), I=174,377)/
     X  3.2200, 3.2600, 3.3000, 3.3200, 3.3400, 3.36,   3.3800, 3.4000,
     X  3.4200, 3.4400, 3.4600, 3.4800, 3.5000, 3.52,   3.5400, 3.56,
     Y  3.5800, 3.6000, 3.6400, 3.66,   3.6800, 3.7000, 3.7400, 3.7800,
     Z  3.8000, 3.8400, 3.9200, 3.9400, 3.9800, 4.0000, 4.02,
     a  4.0400, 4.0600, 4.0800, 4.1000, 4.1100,
     b  4.1200, 4.1300, 4.1400, 4.1500, 4.1600, 4.1700, 4.1800, 4.1900,
     c  4.2000, 4.2200, 4.2300, 4.2400, 4.2500, 4.2600, 4.2700, 4.2800,
     d  4.2900, 4.3000, 4.3100, 4.3300, 4.3500, 4.3700, 4.3900, 4.4000,
     e  4.4200, 4.4400, 4.4600, 4.4800, 4.5000, 4.5400, 4.5600, 4.6000,
     f  4.6400, 4.6800, 4.7000, 4.7400, 4.8000, 4.8200, 4.8500, 4.9000,
     K  4.9194, 4.9300, 4.9345, 4.9355, 4.9367, 4.9385, 4.9395, 4.9432,
     L  4.9481, 4.9546, 4.9584, 4.9698, 4.9800, 5.0000, 5.0300, 5.1000,
     M  5.1200, 5.1800, 5.2000, 5.2300, 5.2800, 5.3000, 5.3300, 5.3350,
     N  5.3400, 5.3600, 5.3620, 5.3710, 5.3780, 5.3900, 5.4000, 5.4100,
     O  5.4200, 5.4400, 5.4600, 5.5000, 5.5500, 5.5530, 5.6000, 5.6500,
     P  5.7,    5.8000, 5.9000, 6.0000, 6.0500, 6.1250, 6.1600, 6.1800,
     Q  6.2000, 6.2100, 6.2200, 6.2300, 6.2400, 6.2500, 6.2850, 6.2950,
     R  6.3030, 6.3100, 6.3200, 6.3300, 6.3400, 6.3500, 6.3600, 6.3700,
     S  6.3900, 6.4000, 6.4100, 6.4200, 6.4300, 6.4400, 6.4500, 6.4700,
     T  6.4900, 6.5100, 6.5400, 6.5600, 6.5700, 6.5800, 6.5900, 6.6000,
     U  6.6200, 6.6400, 6.6575, 6.6650, 6.6800, 6.7000, 6.7500, 6.8100,
     V  6.9200, 7.0000, 7.1000, 7.1400, 7.1800, 7.2000, 7.2200, 7.2250,
     W  7.2500, 7.2700, 7.3400, 7.3500, 7.3700, 7.4000, 7.4200, 7.4700,
     X  7.5417, 7.5937, 7.6200, 7.6500, 7.6674, 7.6800, 7.7000, 7.7250,
     Y  7.7450, 7.7500, 7.7700, 7.7887, 7.8100, 7.8191, 7.8600, 7.8884,
     Z  7.8971, 7.9300, 8.0000, 8.0440, 8.0800, 8.1000, 8.1050, 8.120/
      Data (E(I), I=378,521)/
     A  8.1380, 8.1660, 8.2000, 8.2100, 8.2400, 8.2800, 8.2960, 8.3200,
     B  8.3300, 8.4000, 8.4260, 8.4500, 8.5000, 8.5200, 8.6000, 8.6641,
     C  8.7000, 8.7681, 8.8000, 8.8500, 8.9198, 8.9400, 8.9800, 9.0000,
     D  9.005,  9.0200, 9.0300, 9.0450, 9.0800, 9.1490, 9.1625, 9.1800,
     E  9.2189, 9.2500, 9.2535, 9.3000, 9.3600, 9.4500, 9.5000, 9.5220,
     F  9.5600, 9.5900, 9.63,   9.6400, 9.6800, 9.6920, 9.7000, 9.7259,
     G  9.7400, 9.7500, 9.8000, 9.8299, 9.9000, 9.9209, 10.000, 10.050,
     H  10.170, 10.250, 10.300, 10.372, 10.400, 10.500, 10.550, 10.620,
     I  10.690, 10.830, 10.940, 11.000, 11.053, 11.100, 11.170, 11.250,
     J  11.400, 11.500, 11.700, 11.800, 11.900, 11.917, 12.000, 12.050,
     K  12.100, 12.224, 12.250, 12.300, 12.400, 12.500, 12.599, 12.70,
     L  12.990, 13.000, 13.100, 13.12,  13.250, 13.300, 13.540, 13.587,
     M  13.700, 13.822, 13.830, 13.965, 14.000, 14.182, 14.250, 14.419,
     N  14.500, 14.566, 14.694, 14.750, 14.767, 14.812, 14.837, 14.863,
     O  14.888, 14.927, 14.962, 15.000, 15.045, 15.093, 15.250, 15.477,
     P  15.731, 15.970, 16.000, 16.068, 16.256, 16.440, 16.695, 16.820,
     Q  16.974, 17.138, 17.300, 17.467, 17.687, 17.900, 18.087, 18.273,
     R  18.632, 18.833, 19.030, 19.185, 19.346, 19.511, 19.660, 20.00/
      Data (Sig(I),I=1,173)/4.7392,4.7392,4.7391,4.7346,4.7161, 4.6991,
     X  4.6821, 4.6653, 4.5989, 4.5662, 4.4862, 4.4084, 4.3326, 4.0491,
     A  3.7937, 3.5626, 3.3527, 3.1615, 2.9868, 2.7888, 2.5774, 2.4503,
     B  2.3331, 2.2249, 2.1250, 2.0328, 1.9479, 1.8698, 1.7981, 1.7321,
     C  1.6704, 1.6591, 1.6849, 1.6951, 1.7089, 1.7278, 1.7542, 1.7918,
     X  1.8165, 1.8467, 1.8836, 1.9294, 1.9868, 2.0595, 2.1525, 2.2736,
     X  2.4328, 2.6449, 2.9297, 3.3122, 3.8176, 4.4511, 5.1565, 5.7664,
     X  6.0435, 5.8719, 5.3728, 4.7654, 4.1989, 3.7282, 3.3559, 3.0658,
     X  2.8399, 2.6626, 2.5218, 2.4085, 2.3166, 2.2407, 2.1776, 2.1245,
     X  2.0794, 2.0407, 2.0073, 1.9783, 1.9303, 1.8926, 1.8624, 1.8376,
     A  1.8170, 1.7996, 1.7848, 1.7609, 1.7424, 1.7347, 1.6863, 1.6611,
     B  1.6445, 1.6320, 1.6136, 1.6004, 1.5911, 1.5853, 1.5829, 1.5843,
     C  1.5899, 1.6008, 1.6165, 1.6390, 1.6709, 1.6905, 1.7132, 1.7396,
     X  1.7702, 1.8059, 1.8478, 1.8979, 1.9270, 1.9607, 1.9786, 1.9903,
     P  2.0070, 2.0389, 2.0718, 2.1400, 2.3280, 3.1605, 5.0709, 2.7762,
     Q  2.2895, 2.1589, 2.1090, 2.0863, 2.0697, 2.0673, 2.0702, 2.0825,
     R  2.1029, 2.1265, 2.1795, 2.2425, 2.2774, 2.3149, 2.3555, 2.3994,
     S  2.4470, 2.4987, 2.5548, 2.6157, 2.6818, 2.7532, 2.8296, 2.9103,
     T  3.0432, 3.1067, 3.1633, 3.2070, 3.2297, 3.2212, 3.1707, 3.0686,
     U  2.9106, 2.7013, 2.4557, 2.1960, 1.9458, 1.7235, 1.5395, 1.3962,
     V  1.2908, 1.2176, 1.1704, 1.1295, 1.1503, 1.1976, 1.2547, 1.3725,
     W  1.4825, 1.5826, 1.6744, 1.76,   1.8404, 1.9171, 1.9904/
      Data (Sig(I),I=174,377)/
     X  2.0609, 2.1936, 2.3132, 2.3678, 2.4183, 2.4644, 2.5058, 2.5420,
     X  2.5729, 2.5984, 2.6182, 2.6324, 2.6410, 2.6442, 2.6422, 2.6353,
     Y  2.6238, 2.6081, 2.5654, 2.5393, 2.5105, 2.4795, 2.4120, 2.3394,
     Z  2.3019, 2.2257, 2.0736, 2.0367, 1.9657, 1.9319, 1.9004,
     a  1.8711, 1.8456, 1.8255, 1.8126, 1.8098,
     b  1.8103, 1.8145, 1.8233, 1.8372, 1.8572, 1.8839, 1.9175, 1.9579,
     c  2.0043, 2.1066, 2.1559, 2.1988, 2.2319, 2.2528, 2.2609, 2.2568,
     d  2.2422, 2.2195, 2.1909, 2.1240, 2.0534, 1.9854, 1.9226, 1.8933,
     e  1.8390, 1.7899, 1.7454, 1.7049, 1.6677, 1.6015, 1.5717, 1.5175,
     f  1.4688, 1.4245, 1.4036, 1.3639, 1.3078, 1.2877, 1.2600, 1.2230,
     g  1.2387, 1.3107, 1.3415, 1.7193, 1.9671, 1.9407, 1.7372, 1.4172,
     L  1.2776, 1.2301, 1.2072, 1.1818, 1.1720, 1.1583, 1.1448, 1.1201,
     M  1.1125, 1.0746, 1.0580, 1.0350,0.99397,0.97397,0.96297,0.97447,
     N  1.0360, 1.5038, 1.5516, 1.7101, 1.5508, 1.3542, 1.2424, 1.1506,
     O  1.1088, 1.0483, 1.0158,0.99097,0.98096,0.98012, 0.971, 0.94796,
     P 0.92596,0.90196,0.89546,0.88695,0.88595,0.89595,0.91222,0.94508,
     Q  1.0129, 1.0430, 1.0830, 1.1529, 1.2430, 1.4230, 2.1414, 2.2325,
     R  2.1251, 1.8340, 1.5599, 1.3796, 1.2490, 1.1593, 1.0766, 1.0023,
     S 0.91364, 0.8635, 0.8280, 0.8120,0.82164, 0.8280, 0.8390, 0.8535,
     T 0.83905, 0.7746,0.63544, 0.5551, 0.5216,0.50227,0.49794, 0.5236,
     U 0.60744,0.68294,0.71994,0.71294,0.67155, 0.6448, 0.6387, 0.6270,
     V  0.5963, 0.5696,0.56836, 0.5822, 0.6180,0.64442, 0.6829, 0.6913,
     W  0.7734, 0.8646, 1.3667, 1.4166, 1.4483, 1.4353, 1.4247, 1.4085,
     X  1.4056, 1.3920, 1.3775, 1.3966, 1.4427, 1.4692, 1.5413, 1.7184,
     Y  1.9003, 1.9170, 1.7906, 1.6925, 1.5909, 1.5687, 1.4917, 1.4642,
     Z  1.4517, 1.4123, 1.3597, 1.2955, 1.27765,1.2459, 1.2464, 1.223/
      Data (Sig(I),I=378,541)/
     A 1.16900,1.09276,0.99394,0.95854,0.88064,0.82764,0.81564,0.79411,
     B 0.78703,0.77581,0.7821, 0.77684,0.77778,0.7679, 0.75438,0.74355,
     C 0.7327, 0.72074,0.70652,0.68353,0.65027,0.63938,0.63623,0.64866,
     D 0.64878,0.68113,0.68937,0.67623,0.66755,0.70213,0.71266,0.72048,
     E 0.69418,0.69261,0.69549,0.75252,0.7096, 0.67788,0.6738, 0.66661,
     F 0.6542, 0.66894,0.63361,0.63353,0.68819,0.69409,0.69591,0.66616,
     G 0.64781,0.64182,0.64082,0.63852,0.62361,0.61916,0.62124,0.61361,
     H 0.56131,0.55811,0.57448,0.61246,0.62573,0.67312,0.68875,0.75263,
     I 0.78372,0.82552,0.88923,0.85707,0.88929,0.91773,0.85653,0.86087,
     J 0.893,  0.89375,0.85125,0.86114,0.90718,0.91882,0.9785, 1.0188,
     K 1.01411,0.91976,0.89953,0.87388,0.8526, 0.84632,0.85985,0.87126,
     L 0.91038,0.9069, 0.873,  0.86616,0.8512, 0.8448, 0.86513,0.84551,
     M 0.79844,0.80251,0.80274,0.80207,0.79686,0.78104,0.78273,0.78597,
     N 0.7946, 0.80265,0.83275,0.8548, 0.85801,0.88437,0.91573,0.92163,
     O 0.90352,0.87374,0.88093,0.89,   0.89651,0.90612,0.9129, 0.934,
     P 0.94242,0.9473, 0.9483, 0.9514, 0.9448, 0.93634,0.92169,0.9023,
     Q 0.8771, 0.8562, 0.8483, 0.8527, 0.8711, 0.8894, 0.9005, 0.9073,
     R 0.9144, 0.921,  0.9339, 0.952,  0.9823, 1.0076, 1.0103, 1.0055,
C     (Next 3 rows for E(neutron) > 20 MeV, as given in next data
C     statement.)
     S 0.9475, 0.930,  0.9205, 0.913,  0.885,  0.91,   0.938,  0.964,
     T 0.884,  0.830,  0.780,  0.730,  0.65,   0.605,  0.545,  0.5,
     U 0.46,   0.389,  0.325,  0.214/
      Data (E(I), I=522,541)/
     S 20.8,   22.0,   24.0,   26.0,   28.6,   29.0,   29.25,  29.59,
     T 30.0,   35.0,   40.0,   45.0,   50.0,   55.0,   60.0,   65.0,
     U 70.0,   80.0,   90.0,  110.0/

      Enn=En
      Sigma=EXTERP(E,Sig,Enn,NN,1)
      IF (Enn .LE. E(10)) Sigma=EXTERP(E,Sig,Enn,NN,4)
      Sigcelas=Sigma
      Return
      END
C =====================================================================
C =====================================================================
      FUNCTION SIGCINEL(En)
C     Purpose is to return a cross section (in barns) for inelastic
C     scattering to the 4.4-MeV level in 12-C by an incident neutron
C     of energy En (in MeV).  Data for En < 15 MeV from ENDF/B-V
C     evaluation.
C
C     For E(neutron) = 20.8, 22.0, 24.0 and 26.0 MeV see Phys. Med.
C     Biol. 29 (1984) 643 for measurements; for En > 26 MeV see same
C     reference for nuclear model predictions.

      Dimension E(122),Sig(122)
      Data NN/122/, Nterp/1/
      Data E/ 4.812, 4.850, 4.900, 4.920, 4.930, 4.940, 4.9500, 4.9800,
     X  5.0000, 5.0300, 5.1000, 5.1200, 5.1500, 5.1800, 5.2000, 5.2300,
     X  5.2800, 5.3600, 5.3700, 5.3800, 5.4300, 5.5000, 5.5500, 5.6000,
     X  5.6500, 5.9000, 6.0500, 6.2000, 6.2500, 6.3200, 6.3400, 6.3500,
     X  6.3600, 6.3900, 6.4100, 6.4300, 6.4500, 6.5400, 6.5600, 6.6200,
     X  6.6400, 6.6700, 6.7500, 6.8100, 6.9200, 7.1400, 7.1800, 7.2200,
     X  7.2500, 7.3600, 7.4200, 7.4700, 7.5417, 7.5937, 7.6674, 7.7887,
     X  7.8191, 7.8884, 7.9361, 8.0000, 8.0140, 8.0440, 8.1000, 8.1380,
     X  8.1660, 8.2000, 8.2400, 8.3200, 8.4260, 8.5000, 8.7500, 8.8330,
     Y  9.0000, 9.0450, 9.1490, 9.2500, 9.5000, 9.6920, 9.7500, 10.000,
     Z  10.250, 10.500, 10.690, 10.750, 10.830, 11.000, 11.250, 11.500,
     A  11.750, 11.909, 12.000, 12.224, 12.599, 13.000, 13.250, 13.500,
     B  13.748, 14.000, 14.500, 14.750, 14.807, 14.863, 14.909, 14.954,
     C  16.443, 18.6,   19.5,   20.000, 20.8,   22.0,   24.0,   26.0,
     D  30.0,   35.0,   40.0,   45.0,   50.0,   60.0,   65.0,   70.000,
     E  80.0,   110.0/
      Data Sig/ 0.0, 0.008, 0.022, 0.028, 0.032, 0.035, 0.037, 0.047,
     X 0.048,  0.047,  0.038,  0.036,  0.040,  0.045,  0.052,  0.066,
     X 0.092,  0.148,  0.150,  0.149,  0.133,  0.124,  0.124,  0.124,
     X 0.1370, 0.1970, 0.2360, 0.2520, 0.2770, 0.340,  0.3510, 0.3490,
     X 0.340,  0.2880, 0.2650, 0.2550, 0.2520, 0.2720, 0.2640, 0.200,
     X 0.1870, 0.1750, 0.1620, 0.1560, 0.1560, 0.1670, 0.1750, 0.1920,
     X 0.2120, 0.310,  0.3510, 0.3480, 0.32156,0.30943,0.31376,0.3857,
     X 0.38917,0.35883,0.3597, 0.3770, 0.400,  0.450,  0.490,  0.460,
     X 0.430,  0.400,  0.390,  0.3450, 0.2660, 0.2450, 0.260,  0.2650,
     A 0.280,  0.310,  0.3140, 0.290,  0.280,  0.3220, 0.350,  0.3250,
     B 0.3150, 0.310,  0.330,  0.3450, 0.3650, 0.360,  0.3350, 0.270,
     C 0.260,  0.250,  0.2390, 0.2270, 0.2140, 0.210,  0.2050, 0.2030,
     D 0.1950 ,0.190,  0.170,  0.1650, 0.190,  0.2150, 0.1983, 0.1817,
     E 0.1630, 0.126,  0.1192, 0.1305, 0.102,  0.093,  0.079,  0.07,
     F 0.0577, 0.0475, 0.0397, 0.0331, 0.0275, 0.0195, 0.0166, 0.015,
     G 0.0126, 0.01/
c
      Enn=En
      Sigma=EXTERP(E,Sig,Enn,NN,Nterp)
      Sigcinel=Sigma
      Return
      END
C =====================================================================
C =====================================================================
      FUNCTION SIGCN2N(En)
C     Purpose is to return a cross-section value (in barns) for the reaction
C     n + 12-C --> n + n + 11-C  for an incident neutron having energy
C     En (MeV).  Data for En to 34 MeV from Zeit. fur Physik A301 (1981)
C     353. See also Phys. Rev. 73 (1947) 262; ibid 265.
C
C     4/87. Add in "cross sections" to account for added capability of
C     the N2N routine to compute some (n,2np) reactions.  Threshold for
C     the 12-C(n,2np)10-B reaction is about 34.5 MeV.

      Dimension E(21),Sig(21)
      Data NN/21/, Nterp/1/
      Data E/  20.3,   22.0,  22.8,  23.9,  25.0,  26.0,  26.7,  28.0,
     U  30.0,  32.0,   36.3,  37.5,  40.0,  42.0,  45.0,  50.0,  56.0,
     W  60.0,  75.0,   90.0, 110.0/
      Data Sig/ 0.0,  0.002, 0.003,0.0063,0.0114,0.0139,0.0163,0.0207,
     W0.0235,0.0256,  0.028, 0.029,0.0313,0.0323, 0.033,0.0333,0.0325,
     X0.0316,0.0283,  0.022, 0.0145/

      Enn=En
      Sigma=EXTERP(E,Sig,Enn,NN,Nterp)
      Sigcn2n=Sigma
      Return
      END
C =====================================================================
C =====================================================================
      FUNCTION SIGCN3HE(En)
C     Purpose to return a cross section (in barns) for  n + 12-C -->
C     3-He + 10-Be  reactions.  Data for En = 39.7 and 60.7 MeV
C     extracted from 3-He spectra published in Phys. Rev. C28, 521
C     (1983).  Data at En = 90 MeV for all 3-He reactions yield a
C     cross section of about 6 mb.

      Dimension E(22), Sig(22)
      Data E/22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 33.5, 39.7, 42.5, 45.0,
     x 47.0, 49.0, 51.0, 53.0,  55.0,  57.5,  60.7, 65.0, 70.0, 76.0,
     y 90.0, 110.0/
C     Next array -- cross sections in mb, not barns.
      Data Sig/0.0, 0.5, 1.2,  2.4,  3.3,  4.6,  5.8,  10.9, 12.7, 14.0,
     x 14.8, 15.5, 15.9, 16.15, 16.25, 16.25, 16.0, 15.3, 14.0, 11.7,
     y 6.0, 2.5/
      Data N/22/,Nterp/1/

      Enn=En
      S=EXTERP(E,Sig, Enn, N,Nterp)
      Sigcn3He=0.001*S
      Return
      END
C =====================================================================
C =====================================================================
      FUNCTION SIGCNAL(En)
C     Purpose is to return a cross-section value (in barns) for the reaction
C     n + 12-C --> alpha + 9-Be (ground state)  by an incident neutron
C     of energy En (MeV).

      Dimension E(45),Sig(45)
      Data NN/ 45/, Nterp/1/
      Data E/ 6.186, 6.340, 7.180, 7.280, 7.340, 7.400, 7.5287, 7.6977,
     X  7.8971, 8.0054, 8.1008, 8.2958, 8.4995, 8.6641, 8.7681, 8.9198,
     X  9.1625, 9.2189, 9.2535, 9.3099, 9.7259, 9.8212, 9.9209, 10.372,
     X  11.004, 11.4,   12.5,   13.5,   14.5,   15.0,   17.0,   19.0,
     Y  22.0,   23.0,   24.0,   27.0,   30.0,   32.5,   40.0,   50.,
     Z  55.0,   60.0,   70.0,   90.0,  110.0/
      Data Sig/ 0.0, 1.0E-05, 0.001,0.006,0.011, 0.018,0.04391,0.09476,
     X  0.1589, 0.171,  0.1589, 0.0965, 0.05952,0.05952,0.07338,0.135,
     X  0.248,  0.293,  0.299,  0.286,  0.186,  0.178,  0.189,  0.1306,
     X  0.08263,0.0785, 0.074,  0.0716, 0.066,  0.062,  0.0392, 0.026,
     Y  0.0175, 0.0155, 0.0142, 0.0112, 0.0095, 0.0084, 0.00605,0.004,
     Z  0.0034, 0.0029, 0.00215,0.0019, 0.0017/

      Enn=En
      IF (Enn .LT. E(3)) goto 2
      Sigma = EXTERP(E,Sig,Enn,NN,Nterp)
      Sigcnal=Sigma
      Return
    2 Sigma = EXTERP(E,Sig,Enn,NN,2)
      Sigcnal=Sigma
      Return
      END
C =====================================================================
C =====================================================================
      FUNCTION SIGCND(En)
C     Purpose is to return a cross-section value (in barns) for the reaction
C     n + 12-C --> d + 11-B for an incident neutron energy En (in MeV).
C     This reaction can occur for neutrons having energies too small for the
C     n + 12-C --> n + p + 11-B reaction; the ENDF/B-V evaluation gives
C     cross sections up to 20 MeV based on indirect evidence (it appears).
C
C     Also tried to fit energy distributions at 27.4, 39.7, and 60.7 MeV
C     of Subramanian et al. Phys. Rev. C28, 521 (1983).  The fits are
C     quite reasonable.
C
C     (2006.01.26 by d.satoh)
C     Cross sections between 110 and 150 MeV were added based on the calculation
C     results of TALYS code. (n, d) + (n, nd)

      Dimension E(25),Sig(25)
      Data NN/25/
      Data E/15.25, 15.48, 15.97, 16.44, 16.97, 17.9,  18.46, 18.75,
     P       19.0,  19.51, 20.0,  21.0,  23.0,  26.5,  33.5,
     Q       42.7,  53.0,  61.0,  70.0,  90.0,
     +       110.0, 120.0, 130.0, 140.0, 150.0/
      Data Sig/0.0, 0.002, 0.02,  0.03,  0.04,  0.06,   0.07, 0.071,
     P       0.068, 0.06,  0.053, 0.0505,0.049, 0.0485, 0.047,
     Q       0.0449,0.042, 0.0395,0.0356,0.025,
     +       0.010501, 0.0085761, 0.0075167, 0.0062766, 0.0054025/

      Enn=En
      Sigma=EXTERP(E,Sig,Enn,NN,1)
      Sigcnd=Sigma
      Return
      End
C =====================================================================
C =====================================================================
      FUNCTION SIGCNN3A(En)
C     Purpose is to return a cross section (in barns) for the reaction
C     n + 12-C --> n + 3 alphas for incident neutron of energy
C     En (in MeV).  For En < 13 MeV the sum of partial (n,n') cross
C     sections are being used.
C
C     It appears that values for 11 to 35 MeV gotton from Nuclear
C     Physics A394 (1983) 87 are too large, especially for En > 17 MeV.

      Dimension E(31),Sig(31)
      Data NN/31/, Nterp/1/
      Data E/  8.4,  8.7,   9.0,   9.5,   10.0,  10.5,  11.0,  12.0,
     W 13.0,  13.5,  14.0,  14.5,  15.0,  16.0,  17.0,  18.0,  20.0,
     Y 23.0,  25.0,  27.0,  29.0,  31.0,  35.0,  40.0,  45.0,  50.0,
     Z 55.0,  60.0,  70.0,  90.0, 110.0/
      Data Sig/0.0,  0.006, 0.013, 0.035, 0.05,  0.061, 0.09,  0.15,
     W 0.205, 0.235, 0.265, 0.292, 0.306, 0.314, 0.308, 0.292, 0.258,
     X 0.21,  0.179, 0.152, 0.134, 0.117, 0.094, 0.073, 0.06,  0.048,
     Z 0.04,  0.034, 0.0267,0.019, 0.012/
      Enn=En
      IF (Enn .LT. E(4)) goto 2
      Sigma=EXTERP(E,Sig,Enn,NN,Nterp)
      Sigcnn3a=Sigma
      Return
    2 Sigcnn3a=EXTERP(E,Sig,Enn,NN,2)
      Return
      END
C =====================================================================
C =====================================================================
      FUNCTION SIGCNP(En)
C     Purpose is to return a cross-section value (in barns) for the reaction
C     n + 12-C --> p + 12-B (the particle stable levels only) for incident
C     neutron of energy En (MeV).
C
C     Cross sections between 15 and 20 MeV are essentially ENDF/B-V values
C     which, in turn, were taken from Rimmer and Fisher, Nuclear Physics
C     A108 (1968) 567.  Values at 27.4, 39.7, and 60.7 MeV adjusted
C     to satisfy spectral data of Subramanian et al. Physical Review C28,
C     521 (1983).  These adjustments took into consideration protons
C     generated in the other subroutines.

      Dimension E(27),Sig(27)
      Data NN/27/
      Data E/13.665, 14.0,   14.5,   15.0, 15.477, 15.966, 16.443,
     Y  16.974,17.467,17.901,18.46,  19.034, 19.522, 20.0, 21.0,
     Z  22.0,  23.0,  25.0,  27.0,   28.4,   30.0,  34.8,  40.0,
     A  45.0,  55.0,  74.0, 110.0/
      Data Sig/0.0, 0.00014, 0.0004, 0.001,  0.004, 0.008, 0.011,
     Y  0.013, 0.016, 0.019, 0.019,  0.018,  0.015, 0.013, 0.0088,
     Z  0.0088,0.01,  0.0121,0.0133, 0.0135, 0.0133,0.012, 0.0108,
     A  0.0098,0.0085,0.0068,0.005/

      Enn=En
      Sigma=EXTERP(E,Sig,Enn,NN,3)
      Sigcnp=Sigma
      Return
      END
C =====================================================================
C =====================================================================
      FUNCTION SIGCNPN(En)
C     Purpose is to return a cross-section value (in barns) for the generic
C     reaction n + 12-C --> p + n + 11-B for incident neutrons having
C     energy En (in MeV).  It includes reactions in the program
C     following breakup of the 11-B when energetically available.
C     It also includes other reactions following breakup of 12-B when
C     energetically available, e.g.  n + 12-C --> p + p + 11-Be.
C     Datum at 90 MeV derived from Kellogg, Phys. Rev. 90, 224 (1953).

      Dimension E(28),Sig(28)
      Data NN/28/, Nterp/1/

      Data E/ 17.35, 18.0,  18.7,  19.3,  20.0,  21.0,  22.0, 23.0,
     W 24.5,  26.0,  28.0,  30.0,  33.0,  36.0,  40.0,  42.0, 45.0,
     X 47.0,  50.0,  53.0,  58.0,  65.0,  70.0,  80.0,  85.0, 90.0,
     Y 100.0, 110.0/
      Data Sig/ 0.0, 0.005, 0.009, 0.012, 0.016, 0.02,  0.0235,0.0275,
     W 0.0355, 0.04, 0.0465,0.05,  0.0545,0.0575,0.06,  0.0615,0.0625,
     X 0.0627,0.0625,0.062, 0.062, 0.0635,0.065, 0.074, 0.083, 0.095,
     Y 0.117, 0.141/

      Enn=En
      Sigma=EXTERP(E,Sig,Enn,NN,Nterp)
      Sigcnpn=Sigma
      Return
      END
C =====================================================================
C =====================================================================
      FUNCTION SIGCNT(En)
C     Purpose to return cross section for (N,T) reaction.  Data values
C     for En = 27.4, 39.7, and 60.7 MeV deduced from graphs shown for
C     double differential cross sections in paper of Subramanian, et al,
C     Physical Review C28, 521 (1983).  Datum at 90 MeV from Kellogg's
C     paper, Phys. Rev. 90, 224 (l953) which includes all reactions
C     resulting in a triton.  And there are quite a few of them.

      Dimension  E(12),Sig(12)
      Data E/ 21.5,22.5,25.0,27.4,30.0, 39.7, 45., 50., 54., 60.7, 90.0,
     k  110.0/
      Data Sig/0.0, 1.0, 3.8, 7.1,10.5, 23.3, 28.5,30.7,31.3,30.5, 22.0,
     k  13.0/

C     Sig values are in millibarns, not barns as in the other functions.
      Data Ndata/12/, Nterp/1/

      Enn=En
      Sigma=EXTERP(E, Sig, Enn, Ndata, Nterp)
      Sigcnt=0.001*Sigma
      Return
      END
C =====================================================================
C =====================================================================
C     LIST OF ROUTIENS IN THIS FILE
C ---------------------------------------------------------------------
C     S: pscat
C     S: cscat
C     S: inelas
C     S: photon
C     S: nalpha
C     S: nn3alf
C     S: np
C     S: npx
C     S: n2n
C     S: n3he
C     B: legendre
C =====================================================================
C =====================================================================
      SUBROUTINE PSCAT
C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------
      Common /NEUTRN/ Eneut, U, V, W
      Common /CSCATT/ Xp,Yp,Zp, Vnca,Vcca,Vcm
!$OMP THREADPRIVATE(/NEUTRN/)
!$OMP THREADPRIVATE(/CSCATT/)

      Common /LGNDRE/ IC, E(220), F(6,220)
C     See file LEGENDRE.FOR for the data in common /LGNDRE/ arrays.
      common /DATA_PSCAT/ A1(220), B1(220), B2(220), B3(220),
     +                    B4(220), B5(220), B6(220)
C     See subroutine CX_init for the data in common /DATA_PSCAT/ arrays. (daiki, 2019.12.02)
      real*8 UNIRN

      Dimension Fi(6)
      Data Nterp/3/, Kn/1/, KC/7/, Nfi/6/

C     For the present programming to obtain the neutron's polar scattering
C     angle the Legendre polynomial "technique" is used.  The BLOCK DATA have
C     coefficients for E(neut) up to 20 MeV from the ENDF/B evaluation;
C     whether the values are correct or not I can't say.  For larger
C     E(neut) values are included which have been deduced from comparisons
C     with experimental angular distributions as 20.8, 26 and 40 MeV.
C
C     It appears, from the O5S coding, that the values of the Legendre
C     coefficients in the Block Data are tabulated for specific incident
C     neutron energies in MeV in the laboratory frame of reference, and the
C     values themselves are for angular distributions in the center of mass.
C
C     Set-up needed only once -- the first time through.
C     CHOOSE THE COSINE OF THE POLAR ANGLE OF THE OUTGOING NEUTRON.
C
   10 En=Eneut
      Vn=VELOCITY(Kn,En)
      CALL CMTRAN(Kn,KC,Vn, Vcom,Enc,Ecc)
      Vcm=Vcom
      Vnca=VELOCITY(Kn,Enc)
      Vcca=VELOCITY(KC,Ecc)
C     Vnca, Vcca = velocities of neutron and Carbon ion in center of mass

      Fi(1)=EXTERP(A1,B1,En,Ic,Nterp)
      Fi(2)=EXTERP(A1,B2,En,Ic,Nterp)
      Fi(3)=EXTERP(A1,B3,En,Ic,Nterp)
      Fi(4)=EXTERP(A1,B4,En,Ic,Nterp)
      Fi(5)=EXTERP(A1,B5,En,Ic,Nterp)
      Fi(6)=EXTERP(A1,B6,En,Ic,Nterp)
C
C     Get polar scattering angle from function CHOOSL
C
      FMU = CHOOSL(Fi,NFi)

   24 SINPSI = SQRT ( 1. - FMU * FMU )

      RANN = 6.28318530*sngl(UNIRN(dummy))
      SINETA = SIN(RANN)
      COSETA = COS(RANN)

      Xp=Sinpsi*Coseta
      Yp=Sinpsi*Sineta
      Zp=Fmu

      CALL CSCAT
C     CSCAT finishes up -- gets new neutron energy, dir. cosines.
      Nelm=2

C ---------------------------------------------------------------------
C     NEUTRON
C ---------------------------------------------------------------------
      IKF   = 2112
      DDENG = Eneut
      DDVX  = U
      DDVY  = V
      DDVZ  = W
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

      IDRTYP = Nelm
      CALL CX_anal

      Return
      END
C =====================================================================
C =====================================================================
      Subroutine CSCAT
C
C     Finish up  n + 12-C  elastic and inelastic scattering.
C
C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------

      Common /CSCATT/ Xp,Yp,Zp, Vnca,Vcca,Vcom
!$OMP THREADPRIVATE(/CSCATT/)
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
!$OMP THREADPRIVATE(/VECTOR/)
      Common /NEUTRN/ Eneut,U,V,W
!$OMP THREADPRIVATE(/NEUTRN/)

      Data Kn/1/, K12C/7/

C     Get neutron velocity components in "neutron" center-of-mass
C     coordinate system:
      Vnx=Xp*Vnca
      Vny=Yp*Vnca
      Vnz=Zp*Vnca

C     Get, similarly, carbon-ion velocity components:
      Vcx=-Xp*Vcca
      Vcy=-Yp*Vcca
      Vcz=-Zp*Vcca

C     Transform neutron velocity components from "neutron" center-of-mass
C     coordinate system to "neutron" laboratory coordinate system:
      CALL LABTRAN(Vnx,Vny,Vnz, Vcom, Vx,Vy,Vz)

C     Now get V(neutron) and E(neutron). V(neutron) has the same value
C     in "neutron" lab. coordinates as in "detector" lab. coordinates.
      Vneut=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Eneut=EFROMV(Kn,Vneut)
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
C     Xpn,Ypn and Zpn are neutron's dir. cosines in "neutron" lab. coords.

C     Rotate neutron's dir. cosines into "detector" lab. coordinates
C     and the information is saved in common block labelled NEUTRN.
      Zx=U
      Zy=V
      Zz=W
      CALL TRANSVEC(U,V,W)
      U=Xn
      V=Yn
      W=Zn

C     Now get information on the Carbon ion:
      CALL LABTRAN(Vcx,Vcy,Vcz, Vcom, Vx,Vy,Vz)
      Vcarb=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Ecarb=EFROMV(K12C,Vcarb)

C ---------------------------------------------------------------------
C     RECOIL CARBON
C ---------------------------------------------------------------------
      if( Ecarb .le. 0.0) return

      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
      IKF   = 6*1000000+12
      DDENG = Ecarb
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C
C     All Done!
C
      Return
      END
C =====================================================================
C =====================================================================
      SUBROUTINE INELAS
C     INELAS DOES THE CALCULATIONS FOR THE 12C INELASTIC REACTION.

C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------
      Common /NEUTRN/ Eneut, Ax, Ay, Az
!$OMP THREADPRIVATE(/NEUTRN/)
      Common /CSCATT/ Xpn,Ypn,Zpn,Vnca,Vcca,Vcm
!$OMP THREADPRIVATE(/CSCATT/)
      real*8 UNIRN

      Dimension A1(10),B1(10),B2(10),B3(10),B4(10),Fi(4)

      Data Ipoly/10/, Nterp/1/, Nfi/4/
      Data Q/4.433/, Rmass/0.922461/, Kn/1/, K12C/7/, TenMeV/10.0/

C     Next arrays to get anisotropic scattering distribution.
C     (See Glasgow, et al, Nuclear Science and Engineering, 61
C     (1976) 521 for data between 9.19 and 13 MeV.)
      Data A1/ 6.0, 8.56, 9.19, 10.69, 10.96, 11.73,
     m  12.95, 14.6, 20.8, 26.0/
      Data B1/0.07784, 0.03034, -0.0396, 0.166, 0.2027, 0.258,
     m 0.1733, 0.21885, 0.39757, 0.55621/
      Data B2/0.04819, 0.17527, 0.1836, 0.203, 0.2504, 0.2406,
     m 0.2843, 0.209, 0.20374, 0.24515/
      Data B3/0.004335, -0.03078, -0.00117, 0.0587, 0.0376, 0.0448,
     m 0.0489, 0.03678, 0.07678, 0.09518/
      Data B4/0.0, 0.01355, 0.0182, 0.0186, 0.0276, 0.0128,
     m 0.0233, -0.001315, 0.01595, 0.01973/

C     Enter with incident neutron energy, direction cosines.
C     Exit with final neutron energy, carbon-ion energy, new direction
C     cosines of scattered neutron.  Spot of interaction is unchanged.

      En=Eneut
      Vn=VELOCITY(Kn,En)
      CALL CMTRAN(Kn,K12C,Vn, Vcom,Enc,Ecc)
      Vcm=Vcom
      Tec=Enc+Ecc
      Ta=Tec-Q
      IF (Ta .GT. 0.0) goto 10

      write(*,5) en
    5 Format(/'  *** Error Subroutine INELAS; E(neut) at entry = '
     V    1PE11.4/)
      Eneut=0.0
      Return

   10 Ena=Rmass*Ta
      IF (Ta .GT. TenMeV) Ena=RCKE(Kn,K12C,Ta)
      Eca=Ta-Ena
      Vnca=VELOCITY(Kn,Ena)
      Vcca=VELOCITY(K12C,Eca)

C     O5S treats inelastic scattering as isotropic for all incident
C     neutron energies.  We will not do the same but will instead
C     use a relatively crude grid of Legendre coefficients to get
C     a handle on the angular distribution of the scattered neutrons.

      Fi(1)=EXTERP(A1,B1,En,Ipoly,Nterp)
      Fi(2)=EXTERP(A1,B2,En,Ipoly,Nterp)
      Fi(3)=EXTERP(A1,B3,En,Ipoly,Nterp)
      Fi(4)=EXTERP(A1,B4,En,Ipoly,Nterp)
      Fmu=CHOOSL(Fi,Nfi)
      Sinpsi=SQRT(1.0 - Fmu*Fmu)
      Phi=6.283185*sngl(UNIRN(dummy))
      Xpn=Sinpsi*COS(Phi)
      Ypn=Sinpsi*SIN(Phi)
      Zpn=Fmu
      CALL CSCAT

      Nelm=3

C ---------------------------------------------------------------------
C     NEUTRON
C ---------------------------------------------------------------------
      IKF   = 2112
      DDENG = Eneut
      DDVX  = Ax
      DDVY  = Ay
      DDVZ  = Az
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

C     Now to check for possible Compton scattering of photon in detector

      Egamma=Q
      CALL PHOTON(Egamma)

      IDRTYP = Nelm
      CALL CX_anal

      Return
      END
C =====================================================================
C =====================================================================
      Subroutine PHOTON(Egamma)

C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------

      E=Egamma
      CALL RVECT(Cx,Cy,Cz)
C ---------------------------------------------------------------------
C     PHOTON
C ---------------------------------------------------------------------
      IKF   = 22
      DDENG = E
      DDVX  = Cx
      DDVY  = Cy
      DDVZ  = Cz
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

      Return
      END
C =====================================================================
C =====================================================================
      SUBROUTINE NALPHA
C     NALPHA DOES THE CALCULATIONS FOR THE REACTION N + 12C -> ALPHA + 9BE.
C     12-C(N,ALPHA)9-BE COLLISION

C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------
      COMMON /VECTOR/ Cax,Cay,Caz, Cx,Cy,Cz
!$OMP THREADPRIVATE(/VECTOR/)
      Common /NEUTRN/ Eneut,U,V,W
!$OMP THREADPRIVATE(/NEUTRN/)
      real*8 UNIRN

      Dimension A1(10),B1(10),B2(10),B3(10),B4(10),B5(10),B6(10),Fi(6)

      Data Q/5.71/, Kn/1/, K12C/7/, K9Be/4/, Ka/3/
      Data Ipoly/10/, Nfi/6/, Nterp/1/

C     Next are normalized Legendre polynomial coefficients, B1 thru B6,
C     for En = A1.  Data for En .LE. 9.83 from G. Dietze et al. in Nuclear
C     Data for Science and Technology, 6-10 Sept 1982, Antwerp.
C     Data for En = 13.9 and 15.6 deduced from angular distributions
C     measured by
C     Data for En = 14.1 MeV deduced from angular distribution measured
C     by Haight et al, Nuclear Science & Eng. 87 (1984) 41.
C     Data for En = 11.5 MeV included to bridge the gap and has no
C     experimental basis.
      Data A1/8.0, 8.64, 8.99, 9.22, 9.41, 9.83, 11.5, 13.9,14.1,15.6/
      Data B1/0.2447, 0.226, -0.1538, -0.1006, -0.1708, -0.04401,
     X 0.06217, 0.1525, 0.24192, 0.06217/
      Data B2/0.07407, -0.1602, 0.06007, 0.1265, 0.1833, 0.1678,
     X -0.058, 0.01694, 0.10631, -0.058/
      Data B3/0.04101, -0.0905, 0.00231, 0.0048, 0.02966, 0.01752,
     X -0.0324, -0.00457, 0.05494, -0.03239/
      Data B4/0.0236, 0.02001, 0.08777, 0.138, 0.1263, 0.04874,
     X 0.0497, 0.04144, 0.09454, 0.04966/
      Data B5/ 0.0, 0.02022, 0.0, 0.00375, -0.01591, -0.105,
     X 0.0259, 0.05025, 0.07588, 0.02585/
      Data B6/0.0, 0.0, 0.0, 0.01643, 0.03383, -0.0293, -0.02105,
     X -0.0102, 0.015443, -0.02105/

      En=Eneut
      Eneut=0.0
      Vn=VELOCITY(Kn,En)
      CALL CMTRAN(Kn,K12C,Vn, Vcom,Enc,Ecc)
      Tec=Enc+Ecc
      Ta=Tec-Q
      IF (Ta .GT. 0.0) goto 10

      write(*,5)En
    5 Format(/'   *** Error in Subroutine NALPHA; E(neut) at entry ='
     Y 1PE11.3/10x,'Set E(Alpha) = 0 and exit'/)
      Return

   10 Ealpha=RCKE(Ka,K9Be,Ta)
      Vac=VELOCITY(Ka,Ealpha)

C     Now get cosine of polar scattering angle.
      Fi(1)=EXTERP(A1,B1,En,Ipoly,Nterp)
      Fi(2)=EXTERP(A1,B2,En,Ipoly,Nterp)
      Fi(3)=EXTERP(A1,B3,En,Ipoly,Nterp)
      Fi(4)=EXTERP(A1,B4,En,Ipoly,Nterp)
      Fi(5)=EXTERP(A1,B5,En,Ipoly,Nterp)
      Fi(6)=EXTERP(A1,B6,En,Ipoly,Nterp)
      Fmu=CHOOSL(Fi,Nfi)
      Sinpsi=SQRT(1.0 - Fmu*Fmu)
      Phi=6.283185*sngl(UNIRN(dummy))
      Xa=Sinpsi*COS(Phi)
      Ya=Sinpsi*SIN(Phi)
      Za=Fmu
      Vax=Xa*Vac
      Vay=Ya*Vac
      Vaz=Za*Vac
      CALL LABTRAN(Vax,Vay,Vaz, Vcom, Vx,Vy,Vz)
      Va=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Ealfa=EFROMV(Ka,Va)
C ---------------------------------------------------------------------
C     Alpha
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Cax,Cay,Caz)
      CALL TRANSVEC(U,V,W)
      IKF   = 2000004
      DDENG = Ealfa
      DDVX  = Cx
      DDVY  = Cy
      DDVZ  = Cz
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     9-Be
      E9Bec=Ta-Ealpha
      V9Bec=VELOCITY(K9Be,E9Bec)
      V9x=-Xa*V9Bec
      V9y=-Ya*V9Bec
      V9z=-Za*V9Bec
      CALL LABTRAN(V9x,V9y,V9z, Vcom, Vx,Vy,Vz)
      V9Be=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      E9Be=EFROMV(K9Be,V9Be)
C ---------------------------------------------------------------------
C     9-Be
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Cax,Cay,Caz)
      CALL TRANSVEC(U,V,W)
      IKF   = 4*1000000+9
      DDENG = E9Be
      DDVX  = Cx
      DDVY  = Cy
      DDVZ  = Cz
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

      Nelm=4
      IDRTYP = Nelm
      CALL CX_anal

      Return
      END
C =====================================================================
C =====================================================================
      Subroutine NN3ALF
C     Purpose is to compute the reaction  n + 12-C --> n' + 3 alphas
C
C ---------------------------------------------------------------------
C     The calculation follows one of several reaction schemes:
C ---------------------------------------------------------------------
C    1. (a)        n + 12-C  -->  n' + 12-C (excited)
C       (b)  12-C (excited)  -->  alpha + 8-Be (ground state)
C       (c)  8-Be (grnd st)  -->  2 alphas.
C
C    2. (a)  same as 1. (a)
C       (b)  12-C (excited)  -->  alpha + 8-Be (excited state at 3.0 MeV)
C       (c)  8-Be (excited)  -->  2 alphas.
C
C    3. (a)  same as 1. (a)
C       (b)  12-C (excited)  -->  3 alphas via 3-body breakup
C
C    4. (a)        n + 12-C  -->  alpha + 9-Be (excited)
C       (b)  9-Be (excited)  -->  n + 8-Be (ground state)
C       (c)  same as 1. (c)
C
C    5. (a)  same as 4. (a)
C       (b)  9-Be (excited)  -->  n + 8-Be (excited state at 3.0 MeV)
C       (c)  same as 2. (c)
C
C    6. (a)  same as 4. (a)
C       (b)  9-Be (excited)  -->  alpha + 5-He
C       (c)  5-He            -->  n + alpha
C
C    7. (a)  same as 4. (a)
C       (b)  9-Be (excited!) --> p + 8-Li
C       (c)  8-Li --> several modes.
C ---------------------------------------------------------------------

C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------

      Common /P3ALF/ Pbar(13)
!$OMP THREADPRIVATE(/P3ALF/)
      Common /NEUTRN/ Eneut, Vnprx, Vnpry, Vnprz
!$OMP THREADPRIVATE(/NEUTRN/)
      Common /VECTOR/ Xb,Yb,Zb, Xn,Yn,Zn
!$OMP THREADPRIVATE(/VECTOR/)
      Common /NEUTR2/ Eneut2, U2,V2,W2
!$OMP THREADPRIVATE(/NEUTR2/)
      Common /MASSES/ Emass(21)
      Common /PPLI8/  Nx9,Ex9(14),Pprot(14)
      real*8 UNIRN

      Dimension Q(14)

      Data Rmass/0.922461/, Rmass2/0.666670/
      Data K5He/11/, K9Be/4/, Kn/1/, KC/7/, K8Be/8/, Ka/3/
      Data TenMeV/10.0/, Q5He/0.88/, Qcont/17.9/
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
      Data Qnn/7.3665/, E8Bex/3.0/, QBa/2.46/, QBn/1.666/
      Data E8Begs/0.092/, Qaa/22.5/, Qna/5.71/

C     Next arrays for determining proton branching following
C     n + 12-C --> alpha + 9-Be.
      Data Ex9/ 16.89, 17.11, 17.33, 18.1, 18.7, 18.94, 19.7,
     a  19.94, 21.1, 26.0, 32.0, 38.6, 44.0, 51.3/
      Data Pprot/ 0.0, .000036, .0134, .0484, .1226, .0965, .156,
     a  .1237, .2,   .22,  .175, .18,   .2,  .25/
      Data Nx9/14/, Nterp/1/, K8Li/16/, Li8pp/16.888/, Li7pn/2.033/
      Data Kp/2/

C ---------------------------------------------------------------------
C     Initial set up.
C ---------------------------------------------------------------------
      Cx=Vnprx
      Cy=Vnpry
      Cz=Vnprz
C ---------------------------------------------------------------------
      Egamma=0.0
      Npgo=0
      L3body=0
      N8Bex=0
C ---------------------------------------------------------------------
      Q(1) = 7.656
      Q(2) = 9.641
      Q(3) = 8.13
      Q(4) = 8.54
      Q(5) = 10.84
      Q(6) = 10.4
      Q(7) = 11.84
      Q(8) = 12.46
      Q(9) = 14.08
      Q(10)= 16.1
      Q(11)= 0.0
      Q(12)= 16.99
      Q(13)= 19.7
      Q(14)= 0.0
C ---------------------------------------------------------------------
C     First step: transform to center-of-mass coordinates
      En=Eneut
      Vn=VELOCITY(Kn,En)
      CALL CMTRAN(Kn,KC,Vn, Vcom,Enc,Ecc)
      Tec=Enc+Ecc
      Ta=Tec-Q(1)
      IF (Ta .GT. 0.0) goto 10
C ---------------------------------------------------------------------
      write(*,5)En
    5 Format(/'   *** Error in Subroutine NN3ALF; E(neut) at entry ='
     & 1PE11.3/10x,'Set E(neut) = E(all alphas) = 0 and exit'/)
      Eneut=0.0
      Return
C     Last is an "error" return.
C ---------------------------------------------------------------------
C     Next is to determine which of 14 reactions to use in rest of
C     computation.
   10 CALL P3ALPH(En)

C     That sets up the -Pbar- array of the P3ALF Common area
C     Choose branching by random number from this array.
      Pran=sngl(UNIRN(dummy))
      Nbrnch=1
      Do 12 J=1,13
      IF (Pran .LE. Pbar(J)) goto 14
   12 Nbrnch=Nbrnch+1
   14 Continue

C     Now have Nbrnch between 1 and 14.  The chosen reaction is
C     as follows:
C
C       Nbrnch    From (n,n')   or   From (n,alpha)
C       ------    -----------        --------------
C          1      Ex = 7.65 MeV
C          2          9.64
C          3                         Ex = 2.43 MeV
C          4                           2.8 + 3.05
C          5         10.84
C          6                              4.70
C          7      11.8+12.7+13.3
C          8                              6.76
C          9         14.08
C         10      16.1 - 18.0
C         11       Continuum
C         12                           11.5 group
C         13                           14.0 group
C         14                           Continuum

      Goto (20,20,40,40,20,40,20,40,20,20,15,40,40,38), Nbrnch

C     Next step is for the 12-C continuum, i.e. excitation of a 12-C
C     "excited state" having Ex > 18 MeV.  We get an approximate
C     n' energy from Function CHOOSN, and then get from that an
C     effective excitation energy for a "level" in the continuum
C     presumably excited in the reaction.

   15 Ta=Tec-Qcont
      Encom=RCKE(Kn,KC,Ta)
C     That's the maximum energy for the "continuum" neutron.

      F=0.065+0.001*En
      Temp=F*En
C     -Temp- determined empirically so that the CHOOSN function
C     deals with a neutron continuum spectrum similar to that
C     deduced by the nuclear model code TNG.
      Elow=0.0
      Entry=CHOOSN(Elow,Encom,Temp)
      Try=Entry*(Emass(Kn) + Emass(KC))/Emass(KC)
      Q(11)=Tec-Try

C     Next step is to get the outgoing neutron energy and direction
C     cosines in detector (laboratory) coordinates.

   20 Ta=Tec-Q(Nbrnch)
      Ena=Rmass*Ta
      IF (Ta .GT. TenMeV) Ena=RCKE(Kn,KC,Ta)
      Eca=Ta-Ena
      Vnca=VELOCITY(Kn,Ena)
      Vcca=VELOCITY(KC,Eca)

C     Choose neutron velocity direction in the center-of-mass
C     coordinates by random number.  An isotropic
C     distribution is assumed.
      CALL RVECT(Xb,Yb,Zb)
      Zx=Vnprx
      Zy=Vnpry
      Zz=Vnprz

C     Get neutron velocity components in "neutron" center-of-mass
C     coordinates:
      Vnxp=Xb*Vnca
      Vnyp=Yb*Vnca
      Vnzp=Zb*Vnca

C     Get carbon velocity components in "neutron" c.o.m. coordinates.
      Vcxp=-Xb*Vcca
      Vcyp=-Yb*Vcca
      Vczp=-Zb*Vcca

C     Transform neutron velocity components into "neutron" laboratory
C     coordinates:
      CALL LABTRAN(Vnxp,Vnyp,Vnzp, Vcom, Vnwx,Vnwy,Vnwz)
      Vneut=SQRT(Vnwx*Vnwx + Vnwy*Vnwy +Vnwz*Vnwz)

C     Now can get E(Neut) and dir. cosines in lab coordinates
C     for the outgoing neutron.
      Eneut=EFROMV(Kn,Vneut)
      CALL DIRCOS(Vnwx,Vnwy,Vnwz, Xb,Yb,Zb)
C     Xb,Yb,Zb are dir. cosines in "neutron" laboratory coordinates.

C     Now get dir. cosines in "detector" laboratory coordinates:
      CALL TRANSVEC(Zx,Zy,Zz)
      Vnprx=Xn
      Vnpry=Yn
      Vnprz=Zn
C ---------------------------------------------------------------------
C     NEUTRON
C ---------------------------------------------------------------------
      IKF   = 2112
      DDENG = Eneut
      DDVX  = Vnprx
      DDVY  = Vnpry
      DDVZ  = Vnprz
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     That takes care of the outgoing neutron, n', of reaction 1. (a)

C     Next step: the "fission" of the excited carbon ion.
C
C     The center-of-mass for this step is
C     the moving excited carbon ion.
      CALL LABTRAN(Vcxp,Vcyp,Vczp, Vcom, Vcx,Vcy,Vcz)
      Vcarb=SQRT(Vcx*Vcx + Vcy*Vcy + Vcz*Vcz)

C     Now test to see if ground state or excited state of 8-Be is
C     involved.
      Ta=Q(Nbrnch)-Qnn-E8Begs
      N8Bex=0
      IF (Nbrnch-5) 25,17,18
   17 IF (sngl(UNIRN(dummy)) - 0.6) 25,25,19
   18 IF (Nbrnch .EQ. 11) goto 25
   19 Ta=Ta-E8Bex
      N8Bex=1
C     N8Bex=0 means ground state involved; N8Bex=1 means excited state.

   25 Eaa=Rmass2*Ta
      IF (Ta .GT. TenMeV) Eaa=RCKE(Ka,K8Be,Ta)

C     Check to see if a contiuum reaction -- if so test for the
C     3-body breakup mode -- if so get an alpha energy (by
C     random number in the Function CHOOSA).
      L3body=0
      IF (Nbrnch .LT. 11) goto 22
      Pcomp=(En - 20.0)/50.
C     The assumption of variable Pcomp: a simple increase in the
C     probability of 3-body breakup with increasing E(neutron)

      IF (sngl(UNIRN(dummy)) .GE. Pcomp) goto 22

C     If it gets to here it's a 3-body breakup reaction.
      L3body=1
      Emax=Eaa
      Eaa=CHOOSA(Emax)
   22 Eba=Ta-Eaa
      Vaa=VELOCITY(Ka,Eaa)

C     The next portion gets the energy of the alpha from the fission
C     reaction 1. (b) or 2. (b) or of the the first alpha from the 3-body
C     breakup reaction 3. (b).
C
C     The angular distribution of the "fission" products is, perforce,
C     randomly chosen.

      CALL RVECT(Xc,Yc,Zc)

C     Velocities are slow enough for these heavier ions to do the
C     transformations non-relativistically.
      Va1x=Xc*Vaa
      Va1y=Yc*Vaa
      Va1z=Zc*Vaa + Vcarb
C     Va1x, Va1y, Va1z are components of alpha velocity in the "carbon-ion"
C     laboratory coordinate system.  Nonrelativistic transformation.

      Valph1=SQRT(Va1x*Va1x + Va1y*Va1y + Va1z*Va1z)
C     Valph1 is the same in the "detector" laboratory coordinate system as
C     it is in the "carbon-ion" laboratory coordinate system.

      Ealph1=EFROMV(Ka,Valph1)

C ---------------------------------------------------------------------
C     1ST ALPHA
C ---------------------------------------------------------------------
      CALL DIRCOS(Vcx,Vcy,Vcz,Xb,Yb,Zb)
      CALL TRANSVEC(Cx,Cy,Cz)
C     (Xn,Yn,Zn) are components of carbon dir cos in the detector lab sys.
      Zx = Xn
      Zy = Yn
      Zz = Zn

      CALL DIRCOS(Va1x,Va1y,Va1z,Xb,Yb,Zb)
      CALL TRANSVEC(Zx,Zy,Zz)
      IKF   = 2000004
      DDENG = Ealph1
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     and that's it for the "first" alpha particle.

C     Now recheck for the 3-body breakup reaction again.
      IF (L3body .EQ. 0) goto 26

C     Next is for 3-body breakup reaction 3. (b).  Get 2d alpha information.
      Ea2=0.5*Eba
      Va2=VELOCITY(Ka,Ea2)
      Vopp=0.5*Vaa
C     (That's the component of velocity of the 2d and 3rd
C     alphas along a z-axis defined by the velocity of the
C     first alpha in the 3-body breakup.)

      Dirzet=Vopp/Va2
      Sintheta=SQRT(1.0-Dirzet*Dirzet)
      Phi=3.1415926*sngl(UNIRN(dummy))

C     For the 2d alpha; the 3rd will scatter at Phi + Pi
      Sinphi=SIN(Phi)
      Cosphi=COS(Phi)

C     Need to transform coordinates from those defined by the z-axis
C     determined by the velocity vector of the first alpha to the
C     "carbon-ion" coordinates.
      Xb=Sintheta*Cosphi
      Yb=Sintheta*Sinphi
      Zb=Dirzet
      Zx=-Xc
      Zy=-Yc
      Zz=-Zc
      CALL TRANSVEC(Zx,Zy,Zz)
C     Xn,Yn,Zn are now direction cosines of 2d alpha in "carbon-ion"
C     center-of-mass coordinates.
      Va2x=Va2*Xn
      Va2y=Va2*Yn
      Va2z=Va2*Zn + Vcarb
C     (Non-relativistic coord. transformation, again.)

      Valph2=SQRT(Va2x*Va2x + Va2y*Va2y + Va2z*Va2z)
      Ealph2=EFROMV(Ka,Valph2)
C ---------------------------------------------------------------------
C     2ND ALPHA
C ---------------------------------------------------------------------
      CALL DIRCOS(Vcx,Vcy,Vcz,Xb,Yb,Zb)
      CALL TRANSVEC(Cx,Cy,Cz)

C     (Xn,Yn,Zn) are components of carbon dir cos in the detector lab sys.
      Zx = Xn
      Zy = Yn
      Zz = Zn

      CALL DIRCOS(Va2x,Va2y,Va2z,Xb,Yb,Zb)
      CALL TRANSVEC(Zx,Zy,Zz)

      IKF   = 2000004
      DDENG = EFROMV(Ka,Valph2)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn

      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     And that's it for the 2d alpha.  Now for the third alpha --
C     its velocity, Va3, equals Va2 -- but in a
C     different direction

      Phi=Phi+3.1415926
      Sinphi=SIN(Phi)
      Cosphi=COS(Phi)

C     Get dir. cosines of 3d alpha in "carbon-ion" c.o.m. coordinates.
      Xb=Sintheta*Cosphi
      Yb=Sintheta*Sinphi
C     Zb is the same as for 2d alpha; so are Zx, Zy, and Zz.
C
      Zb=Dirzet
      Zx=-Xc
      Zy=-Yc
      Zz=-Zc

      CALL TRANSVEC(Zx,Zy,Zz)
C     Now get velocity components of 3d alpha in "carbon-ion" lab. coords.
      Va3x=Va2*Xn
      Va3y=Va2*Yn
      Va3z=Va2*Zn + Vcarb
      Valph3=SQRT(Va3x*Va3x + Va3y*Va3y + Va3z*Va3z)

C ---------------------------------------------------------------------
C     3RD ALPHA
C ---------------------------------------------------------------------
      CALL DIRCOS(Vcx,Vcy,Vcz,Xb,Yb,Zb)
      CALL TRANSVEC(Cx,Cy,Cz)
C     (Xn,Yn,Zn) are components of carbon dir cos in the detector lab sys.
      Zx = Xn
      Zy = Yn
      Zz = Zn

      CALL DIRCOS(Va3x,Va3y,Va3z,Xb,Yb,Zb)
      CALL TRANSVEC(Zx,Zy,Zz)
      IKF   = 2000004
      DDENG = EFROMV(Ka,Valph3)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn

      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      Goto 60
C     End of 3-body breakup energetics and kinematics.

C     Next branch is "fission" of the moving 8-Be ion, 1. (b) or 2. (b).
C     The center of mass is the moving 8-Be ion.

   26 Vba=VELOCITY(K8Be,Eba)
      Vbwx=-Xc*Vba
      Vbwy=-Yc*Vba
      Vbwz=-Zc*Vba + Vcarb
      VBe = SQRT(Vbwx*Vbwx + Vbwy*Vbwy + Vbwz*Vbwz)
C     Now recall whether Ground state or excited state of 8-Be
C     is involved.
   21 Ta=E8Begs
      IF (N8Bex .EQ. 1) Ta=Ta+E8Bex
   28 Ea2=0.5*Ta
   29 Va2=VELOCITY(Ka,Ea2)
C     Again, perforce, isotropy in angular distribution of the two
C     alphas in the center-of-mass coordinates.
      CALL RVECT(Xc,Yc,Zc)
      Va2x=Xc*Va2
      Va2y=Yc*Va2
      Va2z=Zc*Va2+Vbe
      Valph2=SQRT(Va2x*Va2x + Va2y*Va2y + Va2z*Va2z)
      Ealph2=EFROMV(Ka,Valph2)

C ---------------------------------------------------------------------
C     2ND ALPHA
C ---------------------------------------------------------------------
      CALL DIRCOS(Vbwx,Vbwy,Vbwz,Xb,Yb,Zb)
      CALL TRANSVEC(Zx,Zy,Zz)
C     (Xn,Yn,Zn) are components of 8-Be dir cos in the detector lab sys.
      Zx = Xn
      Zy = Yn
      Zz = Zn

      CALL DIRCOS(Va2x,Va2y,Va2z,Xb,Yb,Zb)
      CALL TRANSVEC(Zx,Zy,Zz)
      IKF   = 2000004
      DDENG = Ealph2
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn

      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     and that's it for the second alpha.
C
C     Va3x=-Va2x and Va3y=-Va2y

      Va3x=-Xc*Va2
      Va3y=-Yc*Va2
      Va3z=-Zc*Va2+Vbe
      Valph3=SQRT(Va3x*Va3x + Va3y*Va3y + Va3z*Va3z)
      Ealph3= EFROMV(Ka,Valph3)

C ---------------------------------------------------------------------
C     3RD ALPHA
C ---------------------------------------------------------------------
      CALL DIRCOS(Va3x,Va3y,Va3z,Xb,Yb,Zb)
      CALL TRANSVEC(Zx,Zy,Zz)
      IKF   = 2000004
      DDENG = Ealph3
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      Goto 60
C     All done there, except to tidy up at very end.

C     Next sections are for the reactions  n + 12-C --> alpha + 9-Be.
C
C     Next 5 statements obtain an effective 9-Be "continuum" energy
C     level very approximately.  Two rationalizations: (1) probability
C     for this branch is small, a few percent of the  n + 3 alpha
C     reaction; and (2) the "continuum" is not known experimentally,
C     and the calculated "continuum" using TNG is very similar to that
C     calculated for the proton "continuum".

   38 Ta=Tec-Qaa
      Eamax=RCKE(Ka,K9Be,Ta)
      Temp=4.0
      Eatry=CHOOSP(Eamax,Temp)
      Try=Eatry*(Emass(Ka) + Emass(K9Be))/Emass(K9Be)
      Etry=RCKE(Ka,K9Be,Try)
      Try=Try*Eatry/Etry
      Q(14)=Tec-Try

C     Before continuing on alpha + 9-Be, determine if subsequent 9-Be
C     decay will be by proton emission.
      Npgo=0
      Excit9=Q(14)-Qna
      IF (Excit9 .LE. Ex9(1)) goto 40
      Prob=EXTERP(Ex9,Pprot,Excit9,Nx9,Nterp)
      IF (sngl(UNIRN(dummy)) .LE. Prob) Npgo=1
C     OK.  Npgo=1 means proton decay of 9-Be; Npgo = 0 means alpha decay
C     of 9-Be.  In either case first work on  alpha + 9-Be energetics.

C     Get information on the alpha.
   40 Ta=Tec-Q(Nbrnch)
      Ealpha=RCKE(Ka,K9Be,Ta)
      Va1=VELOCITY(Ka,Ealpha)
      CALL RVECT(Xa,Ya,Za)
      Va1x=Xa*Va1
      Va1y=Ya*Va1
      Va1z=Za*Va1
C     Now transform to neutron laboratory coordinates.
      CALL LABTRAN(Va1x,Va1y,Va1z, Vcom, Vx,Vy,Vz)
C     ... done.  Vx, Vy, Vz are velocity components of alpha in the
C     neutron lab. coordinates.

      Va=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Ea1=EFROMV(Ka,Va)

C ---------------------------------------------------------------------
C     1ST ALPHA
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xb,Yb,Zb)
      CALL TRANSVEC(Cx,Cy,Cz)
      IKF   = 2000004
      DDENG = Ea1
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     That's it for the (first) alpha.

C     Now get 9-Be motion into detector coordinates.  We have to get back
C     to detector coordinates because we still have to get information
C     on the neutron in those coordinates.
      E9Be=Ta-Ealpha
      V9Bec=VELOCITY(K9Be,E9Be)
      V9xc=-Xa*V9Bec
      V9yc=-Ya*V9Bec
      V9zc=-Za*V9Bec
      CALL LABTRAN(V9xc,V9yc,V9zc, Vcom, V9x,V9y,V9z)
      V9Be=SQRT(V9x*V9x + V9y*V9y + V9z*V9z)
      CALL DIRCOS(V9x,V9y,V9z, Xb,Yb,Zb)
      Zx=Vnprx
      Zy=Vnpry
      Zz=Vnprz
      CALL TRANSVEC(Zx,Zy,Zz)
C     Now have 9-Be in detector laboratory coordinates; its direction
C     cosines are  Xn, Yn, and Zn.

C     First decide if next reaction is reaction 7. (b).
      IF (Npgo .EQ. 1) goto 75

C     Program counter to here, then not 7. (b), so try for 6. (b).
      IF (Nbrnch .EQ. 3) goto 48
      IF (Nbrnch .EQ. 6  .AND.  sngl(UNIRN(dummy)) .LT. 0.87) goto 50
      IF (Nbrnch .EQ. 8  .AND.  sngl(UNIRN(dummy)) .LT. 0.45) goto 50
      IF (Nbrnch .EQ. 12 .AND.  sngl(UNIRN(dummy)) .LT. 0.88) goto 50
      IF (Nbrnch .EQ. 13 .AND.  sngl(UNIRN(dummy)) .LT. 0.60) goto 50

C     If program gets to here it wasn't  6. (b).
C     Choice of  4. (b)  or  5. (b)  depends only on the 9-Be level.
      Qq=Qnn+E8Begs
C     So, next is "fission" of 9-Be --> n + 8-Be.  Determine which
C     8-Be state is involved.
      N8Bex=0
      IF (Nbrnch .LT. 8) goto42
      N8Bex=1
      Qq=Qq+E8Bex
   42 Ta=Q(Nbrnch)-Qq

C     Get information on the neutron first.
      Enn=RCKE(Kn,K8Be,Ta)
      Vnn=VELOCITY(Kn,Enn)
      CALL RVECT(Xa,Ya,Za)
      Vnx=Xa*Vnn
      Vny=Ya*Vnn
      Vnz=Za*Vnn

C     Transform velocity components to  9-Be  laboratory coordinates.
      CALL LABTRAN(Vnx,Vny,Vnz, V9Be, Vx,Vy,Vz)
C     ... done.  Vx, Vy, Vz are the desired velocity components.
      Vneutr=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Eneut=EFROMV(Kn,Vneutr)
      CALL DIRCOS(Vx,Vy,Vz, Xb,Yb,Zb)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
C     and now have neutron's direction cosines in detector coordinates

      Vnprx=Xn
      Vnpry=Yn
      Vnprz=Zn

C ---------------------------------------------------------------------
C     NEUTRON
C ---------------------------------------------------------------------
      IKF   = 2112
      DDENG = Eneut
      DDVX  = Vnprx
      DDVY  = Vnpry
      DDVZ  = Vnprz
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     so that takes care of the neutron for reactions  4. (b) and 5. (b)

C     Now need to get motion of 8-Be in the 9-Be lab. coordinates.
      E8Be=Ta-Enn
      V8Be=VELOCITY(K8Be,E8Be)
      V8x=-Xa*V8Be
      V8y=-Ya*V8Be
      V8z=-Za*V8Be + V9Be
      VBe=SQRT(V8x*V8x + V8y*V8y + V8z*V8z)
C     Rest of computation same as for 1. (c) or 2. (c)

C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
      Vbwx=V8x
      Vbwy=V8y
      Vbwz=V8z
C ---------------------------------------------------------------------
      Goto 21

C     Coming down the home stretch!!  Next is for  9-Be --> alpha + 5-He
C     reaction  6. (b).

C     First look at case where the Ex=2.43 MeV state of 9-Be decays into
C     alpha + 5-He.  This is clearly an unusual decay mode, since the Q-
C     value for decay of 9-Be into alpha + 5-He is 2.46 MeV, or some 30
C     keV larger than Ex.  However both the level in 9-Be and the ground
C     state in 5-He are quite broad, and the decay mechanism must take
C     these widths into account.  Programming in such widths seems a bit
C     much for the present purpose.  Instead we will simply note that on
C     the average there is no residual kinetic energy in the 9-Be center
C     of mass to be shared by the alpha and the 5-He ion after the decay.
C     Hence V(alpha2) = V(5-He) = V(9-Be) in the lab. coordinates

   48 Ea2=EFROMV(Ka,V9Be)
C ---------------------------------------------------------------------
C     2ND ALPHA
C ---------------------------------------------------------------------
C     (Xn,Yn,Zn) are components of 9-Be dir cos in the detector lab sys.

      CALL DIRCOS(V9x,V9y,V9z,Xb,Yb,Zb)
      Zx = Xn
      Zy = Yn
      Zz = Zn
      CALL TRANSVEC(Zx,Zy,Zz)
      IKF   = 2000004
      DDENG = Ea2
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

      V5He=V9Be
C     Also the dir. cosines of the 5-He ion in the detector lab. coords.
C     are the same as the dir. cosines of the 9-Be (i.e., Xn,Yn,Zn).

      Goto 55

C     Decay of other 9-Be levels will involve some c.o.m. kinetic energy.
   50 Qq=Qnn+QBa-QBn
      Ta=Q(Nbrnch)-Qq
      Ea2c=RCKE(Ka,K5He,Ta)
C     That's the energy of the 2d alpha in the 9-Be center of mass.
      Va2c=VELOCITY(Ka,Ea2c)
      CALL RVECT(Xa,Ya,Za)
C     Again, isotropy in the decay process in the center of mass.
      Va2x=Xa*Va2c
      Va2y=Ya*Va2c
      Va2z=Za*Va2c + V9Be
C     and again, a non-relativistic transformation to lab. coordinates.

      Va2=SQRT(Va2x*Va2x + Va2y*Va2y + Va2z*Va2z)
      Ea2=EFROMV(Ka,Va2)

C ---------------------------------------------------------------------
C     2ND ALPHA
C ---------------------------------------------------------------------
C     (Xn,Yn,Zn) are components of 9-Be dir cos in the detector lab sys.

      CALL DIRCOS(Va2x,Va2y,Va2z,Xb,Yb,Zb)
      Zx = Xn
      Zy = Yn
      Zz = Zn
      CALL TRANSVEC(Zx,Zy,Zz)
      IKF   = 2000004
      DDENG = Ea2
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C     Recall 9-Be ion has dir. cosines Xn,Yn,Zn in detector lab. coords
      Xn = Zx
      Yn = Zy
      Zn = Zy
C ---------------------------------------------------------------------
C     So -- finished with alpha from 9-Be decay.  Now set up the motion
C     of the 5-He ion in detector lab. coordinates (remember, we still
C     have the neutron's parameters to get into detector lab coords.).

      E5Hec=Ta-Ea2c
      V5Hec=VELOCITY(K5He,E5Hec)
      V5Hex=-Xa*V5Hec
      V5Hey=-Ya*V5Hec
      V5Hez=-Za*V5Hec + V9Be
C     and once more a non-relativistic transformation, this time to get
C     the motion of the 5-He ion in the 9-Be lab. coordinates.

      V5He=SQRT(V5Hex*V5Hex + V5Hey*V5Hey + V5Hez*V5Hez)
      CALL DIRCOS(V5Hex,V5Hey,V5Hez, Xb,Yb,Zb)
C     Recall 9-Be ion has dir. cosines Xn,Yn,Zn in detector lab. coords.

      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)

C     Now we have the 5-He ion dir. cosines in detector lab. coordinates;
C     the direction cosines are Xn, Yn, Zn.  So now to the decay of the
C     5-He into a neutron and an alpha.  Get neutron information first.
   55 Ta=Q5He
      Ena=RCKE(Kn,Ka,Ta)
      Vnc=VELOCITY(Kn,Ena)
      CALL RVECT(Xa,Ya,Za)
C     once more, isotropy in the decay.
      Vnxc=Xa*Vnc
      Vnyc=Ya*Vnc
      Vnzc=Za*Vnc
C     Transform neutron's velocity components from 5-He c.o.m. coords.
C     to 5-He lab. coords.
      CALL LABTRAN(Vnxc,Vnyc,Vnzc, V5He, Vx,Vy,Vz)
C     Got them.  Vx,Vy,Vz are neutron's velocity components in 5-He c.o.m.
      Vneut=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Eneut=EFROMV(Kn,Vneut)
C     Eneut=E(neutron) in lab. coords.  Now get dir. cosines in detector
C     lab. coordinates.

      CALL DIRCOS(Vx,Vy,Vz, Xb,Yb,Zb)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
      Vnprx=Xn
      Vnpry=Yn
      Vnprz=Zn

C ---------------------------------------------------------------------
C     NEUTRON
C ---------------------------------------------------------------------
      IKF   = 2112
      DDENG = Eneut
      DDVX  = Vnprx
      DDVY  = Vnpry
      DDVZ  = Vnprz
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     Whew!!  That's it for the neutron.  Now for the final alpha.
      Ea3c=Ta-Ena
      Va3c=VELOCITY(Ka,Ea3c)
      Va3x=-Xa*Va3c
      Va3y=-Ya*Va3c
      Va3z=-Za*Va3c
      CALL LABTRAN(Va3x,Va3y,Va3z, V5He, Vx,Vy,Vz)
      Valph3=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C ---------------------------------------------------------------------
C     3RD ALPHA
C ---------------------------------------------------------------------
C     (Zx,Zy,Zz) are components of 5He dir cos in the detector lab sys.

      CALL DIRCOS(Vx,Vy,Vz,Xb,Yb,Zb)
      CALL TRANSVEC(Zx,Zy,Zz)
      IKF   = 2000004
      DDENG = EFROMV(Ka,Valph3)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

C     All that's left is to finish up on the third alpha.
   60 Ealph3=EFROMV(Ka,Valph3)
      Goto 90
C     END of (n,n 3-alpha) reaction computation.

C ---------------------------------------------------------------------
C     Next portion of programming for 9-Be --> p + 8-Li decay.
C     Set outgoing neutron energy to zero to begin with.
   75 Eneut=0.0
      Egamma=0.0
      Eleftov=0.0
      CALL PP8LI(Excit9,V9Be,Egamma)
C     That's it.

C     Final bit of business.

   90 Nelm=5
      Eneut2=ABS(Eneut2)
      IF (Egamma .GT. 0.0) CALL PHOTON(Egamma)

      IDRTYP = Nelm
      CALL CX_anal

      Return
      END
C =====================================================================
C =====================================================================
      Subroutine NP
C     Purpose is to compute the energy of the outgoing proton from the
C     reaction  n + 12-C --> p + 12B.
C     The ground state and first four excited states of 12-B are stable
C     against particle emission; the programming chooses which of these
C     five states is the state excited by the incident neutron from a
C     table of probabilities generated from cross sections derived from
C     the Hauser-Feshback code TNG developed by C. Y. Fu.

C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
!$OMP THREADPRIVATE(/VECTOR/)
      Common /NEUTRN/ SPDSQ, U, V, W
!$OMP THREADPRIVATE(/NEUTRN/)
      Common /MASSES/ Emass(21)
      Common /PPOLAR/ X(37),Ptheta(37),Npolar
      real*8 UNIRN

      Dimension Ei(12),Pgs(12),P1ex(12),P2ex(12),P3ex(12),Qex(5)
      Dimension Px(4)
      Equivalence (Qex(1),Q)

C     Next are Q values for the 12-B ground and first 4 excited states.
      Data Qex/12.613, 13.566, 14.287, 15.233, 15.333/

C     Next arrays are used to determine which bound state of 12-B was
C     populated.  Pgs, P1ex, P2ex and P3ex are cumulative probabilities.
      Data Ei/14.75, 15.5, 16.5, 17.0, 18.0, 18.25, 19.0, 20.0, 21.0,
     Z  23.0, 25.0, 40.0/
      Data Pgs/1.0, 0.691, 0.468, 0.368, 0.287, 0.279, 0.246, 0.224,
     Z  0.215, 0.203, 0.198, 0.180/
      Data P1ex/1.0, 1.0, 0.85, 0.779, 0.658, 0.646, 0.592, 0.552,
     Z  0.53, 0.511, 0.502, 0.475/
      Data P2ex/1.0, 1.0, 1.0, 0.985, 0.904, 0.894, 0.857, 0.827,
     Z  0.808, 0.79, 0.782, 0.762/
      Data P3ex/1.0, 1.0, 1.0, 0.998, 0.976, 0.974, 0.966, 0.959,
     Z  0.955, 0.95, 0.947, 0.939/

      Data Npolar/37/
C     Next array is used also in Subroutines NPN and N2N.
C     This table is COS(Theta) for theta every 5 degrees.
      Data X /1.0, .99619, .98481, .965926, .939693,
     &  0.90631, .866025, .81915, .76604, .707107, .64279, .57358,
     &  0.5, .42262, .34202, .25882, .17305, .081557, 0.0, -.081557,
     &  -.17305, -.25882, -.34202, -.42262, -0.5, -.57358, -.64279,
     &  -.707107, -.76604, -.81915, -.866025, -.90631, -.939693,
     &  -.965926, -.98481, -.99619, -1.0/

C     The next array is used to determine the polar scattering angle
C     of the proton, and is taken from measurements reported in
C     Nuclear Instruments and Methods 129 (1975) 241.  These data are
C     for En(lab)=56 MeV corresponding to En(c.o.m.)=51.7 MeV; the
C     measured angular distribution is strongly forward peaked for
C     all measured proton energies.
      Data Ptheta/0.0, 0.0091, 0.0385, 0.091, 0.1601, 0.2333, 0.3069,
     S 0.3782, 0.4456, 0.5078, 0.5644, 0.6151, 0.660, 0.6994, 0.7335,
     T 0.7629, 0.7879, 0.8091, 0.8269, 0.842, 0.857, 0.8717, 0.8861,
     U 0.9001, 0.9134, 0.9262, 0.9382, 0.9495, 0.9595, 0.9687, 0.9768,
     V 0.9838, 0.9896, 0.9941, 0.9974, 0.9993, 1.0/

      Data Npx/12/, Nterp/1/
      Data Kn/1/, K12C/7/, Kp/2/, K12B/6/

      En=SPDSQ
      Vn = VELOCITY(Kn,En)
      CALL CMTRAN(Kn,K12C,Vn, Vcom,Eneut,Ecc)
      Tec=Ecc+Eneut
      Ta=Tec-Q
      IF (Ta .GT. 0.0) goto 2

      write(*,100)En
  100 Format(/'   ***  Error in Subroutine NP; E(neut) at entry ='
     b 1PE11.3/10x,'Set E(Neut) = E(Proton) = 0 and exit'/)
      SPDSQ=0.0
      Return

C     CHOOSE C.M. PROTON ENERGY
    2 Px(1)=EXTERP(Ei,Pgs,En,Npx,Nterp)
      Px(2)=EXTERP(Ei,P1ex,En,Npx,Nterp)
      Px(3)=EXTERP(Ei,P2ex,En,Npx,Nterp)
      Px(4)=EXTERP(Ei,P3ex,En,Npx,Nterp)
      Eran=sngl(UNIRN(dummy))
      Eleftov=0.0
      K=1
      Do 4 J=1,4
      IF (Eran .LE. Px(J)) goto 6
    4 K=K+1
    6 Ta=Tec-Qex(K)
      Egamma=Qex(K)-Qex(1)
      IF (K .NE. 4) goto 7
C     Third excited state of 12-B decays primarily to the first excited
C     state.  Other excited states decay primarily to the ground state.
      Egamma=Qex(2)-Qex(1)
      Eleftov=Qex(4)-Qex(2)
    7 Eprot=RCKE(Kp,K12B,Ta)

C     FIND THE COSINE  OF THE POLAR DIRECTION OF THE PROTON(C.M.)
   15 Rand=sngl(UNIRN(dummy))
   18 Ctheta = EXTERP(Ptheta,X,Rand,Npolar,Nterp)

C     Assumed angular distribution at En(c.o.m.) .GE. 51.7 MeV, but
C     modified toward isotropy for En(c.o.m.) < 51.7 MeV.
      IF (Tec .GE. 51.7) goto 20
      Ciso=1.0-2.0*Rand
      Slope=(Ctheta-Ciso)/(51.7-Qex(K))
      Ctheta=Slope*Ta
C     FIND THE SINE AND COSINE OF THE AZIMUTHAL DIRECTION OF THE
C     PROTON (C.M.) -- By Random-number Choice.
   20 Rann = 6.2831853*sngl(UNIRN(dummy))
      Sinphi = SIN(Rann)
      Cosphi = COS(Rann)
C     FIND THE C. M. DIRECTION COSINES OF THE PROTON.
      Sintheta=SQRT(1.-Ctheta*Ctheta)
      Xp=Sintheta*Cosphi
      Yp=Sintheta*Sinphi
      Zp=Ctheta
      Vpr=VELOCITY(Kp,Eprot)
C     NOW GET C. M. VELOCITY COMPONENTS OF THE PROTON.
      Vpx=Xp*Vpr
      Vpy=Yp*Vpr
      Vpz=Zp*Vpr
      CALL LABTRAN(Vpx,Vpy,Vpz, Vcom, Vx,Vy,Vz)
C     THAT HAS GOT THE VELOCITY COMPONENTS OF THE PROTON IN THE
C     "NEUTRON" LABORATORY COORDINATE SYSTEM.
      Vprot=SQRT(Vx*Vx +Vy*Vy +Vz*Vz)
C     FIND THE LAB ENERGY OF THE PROTON.
      Epl=EFROMV(Kp,Vprot)

C ---------------------------------------------------------------------
C     PROTON
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)
      IKF   = 2212
      DDENG = Epl
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

C     FIND THE COMPONENTS OF THE 12B VELOCITY IN THE LAB.
   70 E12B=RCKE(K12B,Kp,Ta)
      V12B=VELOCITY(K12B,E12B)
      Vbcx=-Xp*V12B
      Vbcy=-Yp*V12B
      Vbcz=-Zp*V12B
C     (Recall that the direction cosines of the 12-B ion in the center
C     of mass are negatives of the direction cosines of the proton.)

      CALL LABTRAN(Vbcx,Vbcy,Vbcz, Vcom, Vbx,Vby,Vbz)
      V12B=SQRT(Vbx*Vbx + Vby*Vby + Vbz*Vbz)
C     FIND THE LAB ENERGY OF THE 12B.
      E12B = EFROMV(K12B,V12B)
C ---------------------------------------------------------------------
C     12B
C ---------------------------------------------------------------------
      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)
      IKF   = 5*1000000 + 12
      DDENG = E12B
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      Nelm=6
      SPDSQ=0.0

C     Check for possible gamma-ray interactions:
      IF (Egamma .GT. 0.0) then
        CALL PHOTON(Egamma)
      ELSE IF (Eleftov .GT. 0.0) then
        CALL PHOTON(Eleftov)
      END IF

      IDRTYP = Nelm
      CALL CX_anal

      Return
      END
C =====================================================================
C =====================================================================
      Subroutine NPX
C     Purpose to compute energetics of the reactions starting with
C     n + 12-C --> p + 12-B.  Program first checks for populating the
C     residual 12-B ion in a low-lying (particle-bound) energy level,
C     and if so program reverts to Subroutine NP.  If not, then
C     emission, and if neither of those then consider neutron emission.
C     The several possible "next" reactions are:
C     12-B --> p + 11-Be,
C     12-B --> alpha + 8-Li, or
C     12-B --> n + 11-B.
C     Following each of these reactions the program will consider
C     possible further particle emissions (if the excitation energy
C     is large enough) of the residual nuclei 11-Be, 8-Li or 11-B.

C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
!$OMP THREADPRIVATE(/VECTOR/)
      Common /NEUTRN/ SPDSQ, U, V, W
!$OMP THREADPRIVATE(/NEUTRN/)
      Common /MASSES/ Emass(21)
      Common /PPOLAR/ X(37),Ptheta(37),N37
C     The arrays X and Ptheta are given in Subroutine NP
      Common /EXC11B/ Exc11(14),TenBpn,Li7pa
      Common /NEUTR2/ Eneut2, U2,V2,W2
!$OMP THREADPRIVATE(/NEUTR2/)
      Common /EXC10B/ Ex10B(9)
C     Array Ex10B is in subroutine NT.
      Common /EXC8LI/ E8Li(10)
C     Array E8Li is in Function EX8LI.
      Common /EXCTBE/ Ex10Be(4),Q9Bepn
      real*8 UNIRN

      Dimension Ex12B(18), PP2P(18)
      Dimension Ex12B1(12),Ppalph(12)

C     Next array has a representation of excited states of 11-B.
      Data Exc11/0.0,2.125,4.445,5.02,6.743,6.792,7.286,7.978,
     a 8.56,8.92,9.2,9.88,10.33,11.25/
C     The last 5 decay by both alpha and gamma.
      Data Qnp/12.613/, Be11pp/14.096/, Be10pn/0.503/, Ka/3/
      Data Q9Bepn/6.812/, Nterp/1/, Npp2p/18/, B11pn/3.37/
      Data Kn/1/, K12C/7/, Kp/2/, K12B/6/, K11B/5/, K11Be/14/, K10Be/15/
      Data K10B/12/, TenBpn/11.46/, Li7pa/8.665/, Li6pa/4.46/, K9Be/4/
      Data Li7pn/2.033/, Li8pa/10.002/, K8Li/16/, Nppalf/12/

C     Next data array represent particle-bound states of 10-Be.
      Data Ex10Be/0.0, 3.368, 5.959, 6.263/

C     Next arrays to determine probability of second proton emission.
C     Calculations using TNG show that  (n,2p)/(n,p)  is a function
C     of the 12-B excitation, but essentially independent of incident
C     neutron energy.
      Data Ex12B/14.1, 14.45, 15.06,  16.0,  16.8,  18.0,  19.2, 20.0,
     a  22.0, 24.0, 25.6, 27.6, 29.4, 32.4, 34.6, 36.4, 39.0, 45.0/
      Data PP2P/0.0,0.0012, 0.0039, 0.011, 0.019, 0.032, 0.042,0.05,
     a 0.068, 0.085, 0.1, 0.13, 0.17, 0.24, 0.30, 0.35, 0.40, 0.42/

C     Next arrays to determine probability of alpha emission following
C     first proton emission.  Again, TNG calculations show dependence
C     on 12-B excitation, but essentially not on incident neutron energy.
      Data Ex12B1/10.003, 10.14, 11.4, 12.5, 13.6, 19.6, 22.0, 26.0,
     a   30.0, 38.0, 50.0, 100.0/
      Data Ppalph/0.0, 0.00036, 0.047, 0.1, 0.15, 0.15, 0.1, 0.075,
     a   0.048, 0.03, 0.024, 0.018/

C     The reaction scheme is:
C
C     (a)  n + 12-C --> p + 12-B;      then:
C     Decide if residual 12-B is in a bound state, and if so
C     go to Subroutine NP.  If not, then:
C     (b)      12-B --> n + 11-B;
C     Check for possible further neutron emission; if so, then:
C     (c)      11-B --> n + 10-B;
C
C     New additions (5/87):
C     (b')     12-B --> p + 11-Be, with possible further emission:
C     (c')     11-Be --> n + 10-Be, with possible further emission:
C     (d')     10-Be --> n + 9-Be (ground state, only)
C
C     Further additions (5/87):
C     (b")     12-B --> alpha + 8-Li, with possible further emission:
C     (c")     8-Li --> gamma or neutron decay.
C
C     Added 9/87:
C     (c3)    11-B --> alpha + 7-Li, with possible further 7-Li breakup.

C     Initial set-up:
      En2=0.0
      K  = 0
      En=SPDSQ
      Vn = VELOCITY(Kn,En)
      CALL CMTRAN(Kn,K12C,Vn, Vcom,Eneut,Ecc)
      Tec=Ecc+Eneut
      Ta=Tec-Qnp
      IF (Ta .GT. 0.0) goto 2

      write(*,222)En
  222 Format(/'   ***  Error in Subroutine NPX; E(neut) at entry ='
     b 1PE11.3/10x,'Set E(Neut) = E(Prot) = 0 and exit.'/)
      SPDSQ=0.0
      Return

C     First step -- decide if particle-stable state of 12-B is involved.
    2 Ratio=SIGCNP(En)/(SIGCNP(En) + SIGCNPN(En))
      IF (sngl(UNIRN(dummy)) .GT. Ratio) goto 7
    3 CALL NP
      Return
C     (The point to this last bit is that originally the (n,p) and
C     (n,pn) reactions were considered entirely separately.  The
C     present amalgamation is, thus, a hybrid...)

C     Set up for first proton emission:
C     CHOOSE C.M. PROTON ENERGY using Function CHOOSP
    7 Tapn=Ta-B11pn
      IF (Tapn .LE. 0.0) goto 3
C     (Last statement just a check; probably shouldn't be needed.)

C     At this point we consider the proton spectrum for En = 60 MeV neutron
C     interactions with 12-C as reported by Brady et al, Journal of Physics
C     G 10 (1984) 363, which shows sufficient enhancement of the population
C     of (unbound) levels at about 4.5 and 8.3 MeV in 12-B to be worth
C     programming.  Thus ...
      IF (En .LE. 34.0) goto 8
      Enhanc=Ratio*(En-34.0)/50.0
      Excit=4.5
      Try=Ta-Excit
      Rani=sngl(UNIRN(dummy))
      IF (Rani .LE. Enhanc) goto 12
      IF (En .LE. 39.0) goto 8
      Excit=8.3
      Try=Ta-Excit
      Enhanc2=Enhanc+Ratio*(En-39.0)/30.0
      IF (Rani .LE. Enhanc2) goto 12

C     Program counter to here -- treat proton output as "continuum"
    8 Epcom=RCKE(Kp,K12B,Tapn)
C     Epcom is the maximum energy the proton can have in the center of mass
      F=0.1245+0.001*ABS(En-45.0)
      Temp=F*En
C     Temp value chosen empirically so that CHOOSP gives a reasonable
C     representation of the proton "continuum" as calculated by the
C     nuclear model code TNG of C. Y. Fu.
      Eprot = CHOOSP(Epcom,Temp)
C     The difference in energy between Epcom and Eprot could be
C     interpreted as the energy of a "virtual" energy level in the 12-B
C     nucleus.  The effect of this "level" should be to increase
C     the value of Q by the "virtual" excitation energy.
      Try=Eprot*(Emass(Kp)+Emass(K12B))/Emass(K12B)
C     That's approximately the "total available energy" that
C     would yield Eprot from the function RCKE.
      Eptry=RCKE(Kp,K12B,Try)
      Try=Try*Eprot/Eptry
C     (Iterate once on TRY for slightly more accuracy.)
      Excit=Ta-Try
C     EXCIT  should be close to the "virtual" excitation energy
C     of the residual 12-B ion.
      IF (Excit .LE. B11pn) goto 3
C     Last is a check, mainly on possible single-precision round-off error.

   12 Qex=Qnp + Excit
      Eprot=RCKE(Kp,K12B,Try)

C     FIND THE COSINE  OF THE POLAR DIRECTION OF THE PROTON(C.M.)
   15 Rand=sngl(UNIRN(dummy))
   18 Ctheta = EXTERP(Ptheta,X,Rand,N37,Nterp)
      IF (Tec .GE. 51.7) goto 20
      Ciso=1.0-2.0*Rand
      Slope=(Ctheta-Ciso)/(51.7-Qex)
      Ctheta=Slope*(Tec-Qex)
C     FIND THE SINE AND COSINE OF THE AZIMUTHAL DIRECTION OF THE
C     PROTON (C.M.) -- By Random-number Choice.
   20 Rann = 6.2831853*sngl(UNIRN(dummy))
      Sinphi = SIN(Rann)
      Cosphi = COS(Rann)
C     FIND THE C. M. DIRECTION COSINES OF THE PROTON.
      Sintheta=SQRT(1.-Ctheta*Ctheta)
      Xp=Sintheta*Cosphi
      Yp=Sintheta*Sinphi
      Zp=Ctheta
      Vpr=VELOCITY(Kp,Eprot)
C     NOW GET C. M. VELOCITY COMPONENTS OF THE PROTON.
      Vpx=Xp*Vpr
      Vpy=Yp*Vpr
      Vpz=Zp*Vpr
      CALL LABTRAN(Vpx,Vpy,Vpz, Vcom, Vx,Vy,Vz)
C     THAT HAS GOT THE VELOCITY COMPONENTS OF THE PROTON IN THE
C     "NEUTRON" LABORATORY COORDINATE SYSTEM.
      Vprot=SQRT(Vx*Vx +Vy*Vy +Vz*Vz)
C     FIND THE LAB ENERGY OF THE PROTON.
      Epl=EFROMV(Kp,Vprot)

c     Get proton dir. cosines in "detector" laboratory coord. system.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)

C ---------------------------------------------------------------------
C     PROTON
C ---------------------------------------------------------------------
      IKF   = 2212
      DDENG = Epl
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     Done.

C     FIND THE COMPONENTS OF THE 12B VELOCITY IN THE LAB.
C     First adjust total available energy for the energy of the
C     "virtual" excited state.
   70 Taex=Tec-Qex
      E12B=RCKE(K12B,Kp,Taex)
      V12B=VELOCITY(K12B,E12B)
      Vbcx=-Xp*V12B
      Vbcy=-Yp*V12B
      Vbcz=-Zp*V12B
C     (Recall that the direction cosines of the 12-B ion in the center
C     of mass are negatives of the direction cosines of the proton.)

      CALL LABTRAN(Vbcx,Vbcy,Vbcz, Vcom, Vbx,Vby,Vbz)
      V12B=SQRT(Vbx*Vbx + Vby*Vby + Vbz*Vbz)

      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)
C     Now variables Xn,Yn,Zn = direction cosines of the motion of the
C     12-B in "detector" laboratory coordinate system.

C     New addition 5/87: check for second proton or alpha emission.
      IF (Excit .LE. Ex12B(1)) goto 72

C     SECOND PROTON OR ALPHA EMISSION
C     If program gets to here, enough excitation in 12-B for second
C     proton or possibly alpha emission.
      Fr=EXTERP(Ex12B,PP2P,Excit,Npp2p,Nterp)
      Rani=sngl(UNIRN(dummy))
      IF (Rani .LE. Fr) goto 71
      Fr=Fr+EXTERP(Ex12B1,Ppalph,Excit,Nppalf,Nterp)
      IF (Rani .GT. Fr) goto 72

C ---------------------------------------------------------------------
C     ALPHA PLUS 8-LI
C     Program counter to here means alpha emission:  12-B --> alpha + 8-Li
C ---------------------------------------------------------------------
      Ta=Excit-Li8pa
      Excita=EX8LI(Excit)
      Egamma=0.0
      Eleftov=0.0
      IF (Excita .GT. Li7pn) goto 50
      Egamma=Excita
C     Either ground state or 1st excited state of 8-Li, and no neutron.
      SPDSQ=0.0
   50 Taex=Ta-Excita
      Epa=RCKE(Ka,K8Li,Taex)
      CALL RVECT(Xpa,Ypa,Zpa)
      Vpa=VELOCITY(Ka,Epa)
      Vacx=Xpa*Vpa
      Vacy=Ypa*Vpa
      Vacz=Zpa*Vpa
      CALL LABTRAN(Vacx,Vacy,Vacz, V12B, Vx,Vy,Vz)
      Valf=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Ealf=EFROMV(Ka,Valf)
C ---------------------------------------------------------------------
C     ALPHA
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      Cx = Xn
      Cy = Yn
      Cz = Zn
      CALL TRANSVEC(Cx,Cy,Cz)
      IKF   = 2000004
      DDENG = Ealf
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
      Xn = Cx
      Yn = Cy
      Zn = Cz
C ---------------------------------------------------------------------
C     That's it for the alpha; now get information on 8-Li ion.
      E8Lic=Taex-Epa
      V8Lic=VELOCITY(K8Li,E8Lic)
      Vax=-Xpa*V8Lic
      Vay=-Ypa*V8Lic
      Vaz=-Zpa*V8Lic + V12B
      V8Li=SQRT(Vax*Vax + Vay*Vay + Vaz*Vaz)
C     Check for further neutron emission by 8-Li ion.
      IF (SPDSQ .GT. 0.0) goto 52
C     Program counter to here = no neutron emission, so tidy up and exit.
C ---------------------------------------------------------------------
C     8Li
C ---------------------------------------------------------------------
      CALL DIRCOS(Vax,Vay,Vaz,Xpn,Ypn,Zpn)
      Cx = Xn
      Cy = Yn
      Cz = Zn
      CALL TRANSVEC(Cx,Cy,Cz)
      IKF   = 3*1000000 + 8
      DDENG = EFROMV(K8Li,V8Li)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
      Xn = Cx
      Yn = Cy
      Zn = Cz
C ---------------------------------------------------------------------
      Goto 90
C ---------------------------------------------------------------------

C     Next, set up for decay of 8-Li.
   52 CALL DIRCOS(Vax,Vay,Vaz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
C     Ready to go to 8-Li decay subroutine.
      CALL ELI8DK(Excita,V8Li,Egamma)
C     Finished breaking up 8-Li, so finish off.
C ---------------------------------------------------------------------
      Goto 90
C ---------------------------------------------------------------------

C ---------------------------------------------------------------------
C     PROTON PLUS 11-BE
C     Got to here: 2d PROTON emission to consider.
C ---------------------------------------------------------------------
   71 Ta=Excit-Be11pp
      Ep2com=RCKE(Kp,K11Be,Ta)
      Ehat=En-Qnp
      F=0.1245+0.001*ABS(Ehat-45.0)
      Temp=F*Ehat
      Eprot2=CHOOSP(Ep2com,Temp)
      Try2=Eprot2*(Emass(Kp) + Emass(K11Be))/Emass(K11Be)
      Ep2try=RCKE(Kp,K11Be,Try2)
C     and as done before, iterate once on Try2 for better precision
      Try2=Try2*Eprot2/Ep2try
      Excitp=Ta-Try2
      Egamma=0.0
      Eleftov=0.0
      IF (Excitp .GT. 1.6) goto 30
C     Low-energy excitation of 11-Be.  Check for E(gamma).
      IF (Excitp .LT. 0.3918) goto 29
      Excitp = 0.3198
      Egamma = 0.3198
C     Only one excited bound state of 11-Be.  However, no subsequent
C       neutron emission -- just (n,2p).
   29 SPDSQ=0.0
   30 Taex=Ta-Excitp
C     Get information on outgoing 2nd proton.
      Epp=RCKE(Kp,K11Be,Taex)
C     Isotropy for 2nd proton emission in center-of-mass coordinates (of
C     moving 12-B ion)
      CALL RVECT(Xp,Yp,Zp)
      Vpp=VELOCITY(Kp,Epp)
      Vpcx=Xp*Vpp
      Vpcy=Yp*Vpp
      Vpcz=Zp*Vpp
      CALL LABTRAN(Vpcx,Vpcy,Vpcz, V12B, Vpx,Vpy,Vpz)
      Vprot=SQRT(Vpx*Vpx + Vpy*Vpy + Vpz*Vpz)
      Eprot=EFROMV(Kp,Vprot)
      CALL DIRCOS(Vpx,Vpy,Vpz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)

C ---------------------------------------------------------------------
C     2ND PROTON
C ---------------------------------------------------------------------
      IKF   = 2212
      DDENG = Eprot
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     Now have 2nd proton direction cosines in "detector" lab. coord. system.
C     Its direction cosines are Xn, Yn, Zn.

C     Now get 11-Be velocity and energy in lab. coord. system
   32 E11Be=Taex-Epp
      V11Be=VELOCITY(K11Be,E11Be)
      Vbx=-Xp*V11Be
      Vby=-Yp*V11Be
      Vbz=-Zp*V11Be + V12B
C     Again, non-relativistic transformation from c.o.m. to lab. coordinates
C     of the moving 12-B ion.

      V11Be=SQRT(Vbx*Vbx + Vby*Vby + Vbz*Vbz)
      IF (SPDSQ .GT. 0.0) goto 35

C     Program counter to here means no further neutron emission.  We
C     excited either the ground state or the 320-keV state in 11-Be.
   33 E11Be=EFROMV(K11Be,V11Be)

C ---------------------------------------------------------------------
C     11Be
C ---------------------------------------------------------------------
      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
      IKF   = 4*1000000 + 11
      DDENG = E11Be
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      Goto 90
C ---------------------------------------------------------------------

C ---------------------------------------------------------------------
C     NEUTRON PLUS 10-BE
C     If program goes to next step, then have (n,2pn) reaction to finish up.
C     Schematic will look something like the drawing in the N2N subroutine.
C ---------------------------------------------------------------------

   35 CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
C     Now have 11-Be ion dir. cosines (i.e. Xn,Yn,Zn) in "detector" lab.
C     coordinate system.

      Ta=Excitp-Be10pn
      Encom=RCKE(Kn,K10Be,Ta)
      Ehat=Ehat-Be11pp
      Fn=0.065+0.001*Ehat
      Temp=Fn*Ehat
      Eloww=0.0
C     Test if Eneut2 has a value left over from a previous "2n" reaction --
C     if so then force present neutron emission to populate a state in
C     10-Be having Ex < binding energy of 9-Be + n.
      IF (Eneut2 .LE. 0.0) goto 36
      Tryy=Encom-Q9Bepn
      IF (Tryy .GT. 0.0) Eloww=Tryy
   36 Enn=CHOOSN(Eloww,Encom,Temp)
      IF (Enn.GT.0.0) goto 37
      SPDSQ=0.0
      Goto 33
   37 Try=Enn*(Emass(Kn) + Emass(K10Be))/Emass(K10Be)
C     Once again, an iteration for precision.
      Entry=RCKE(Kn,K10Be,Try)
      Try=Try*Enn/Entry
      Excitn=Ta-Try
      Egamma=0.0
      Eleftov=0.0
      IF (Excitn .GT. Q9Bepn) goto 44

      K=1
C     Now pair -Excitn- with a level energy for 10-Be.
      Do 40 I=2,4
      IF (Excitn .LT. Ex10Be(I)) goto 42
   40 K=K+1
   42 Excitn=Ex10Be(K)
C     Set up gamma decay characteristics:
      IF (K .GE. 2) Egamma=Ex10Be(2)
      IF (K.EQ.3 .AND. sngl(UNIRN(dummy)).GE.0.5) Egamma=Excitn
      Eleftov=Excitn-Egamma

C     Now get lab. coord. energy of emitted neutron.
C     Again, isotropy of emission in center of mass (of the 11-Be ion).
   44 Taex=Ta-Excitn
      Enn=RCKE(Kn,K10Be,Taex)
      CALL RVECT(Xp,Yp,Zp)
      Vnn=VELOCITY(Kn,Enn)
      Vncx=Xp*Vnn
      Vncy=Yp*Vnn
      Vncz=Zp*Vnn
      CALL LABTRAN(Vncx,Vncy,Vncz, V11Be, Vx,Vy,Vz)
      Vneut=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Eneut=EFROMV(Kn,Vneut)
      SPDSQ=Eneut
C     Done.  Now get neutron dir. cosines into "Detector" lab. coordinates.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
C     Done.  Xn,Yn,Zn are dir. cosines of neutron in "detector" lab. coords.
      U=Xn
      V=Yn
      W=Zn
C ---------------------------------------------------------------------
C     NEUTRON
C ---------------------------------------------------------------------
      IKF   = 2112
      DDENG = SPDSQ
      DDVX  = U
      DDVY  = V
      DDVZ  = W
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

C     Now get 10-Be energy in laboratory coordinates.
      E10Be=Taex-Enn
      V10Be=VELOCITY(K10Be,E10Be)
      Vbx=-Xp*V10Be
      Vby=-Yp*V10Be
      Vbz=-Zp*V10Be + V11Be
C     Non-relativistic transformation from c.o.m. to lab. coords.
      V10Bel=SQRT(Vbx*Vbx + Vby*Vby + Vbz*Vbz)

C ---------------------------------------------------------------------
C     NEUTRON PLUS 9-BE
C     Test for possible SECOND NEUTRON emission.
C ---------------------------------------------------------------------
      IF (Excitn .LE. Q9Bepn) goto 60
C ---------------------------------------------------------------------
C     SECOND NEUTRON:  10-Be --> n + 9-Be
      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
C     Now Xn,Yn,Zn are dir. cosines of 10-Be ion in "detector" lab. coords.

      Ta=Excitn-Q9Bepn
      En2c=RCKE(Kn,K9Be,Ta)
      CALL RVECT(Xn2,Yn2,Zn2)
      Vn2=VELOCITY(Kn,En2c)
      Vn2x=Xn2*Vn2
      Vn2y=Yn2*Vn2
      Vn2z=Zn2*Vn2
      CALL LABTRAN(Vn2x,Vn2y,Vn2z, V10Bel, Vx,Vy,Vz)
      Vneut=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      En2=EFROMV(Kn,Vneut)
      Eneut2=-En2
C     This negative is temporary, just to alert BANKER routine.

C     Now get direction cosines of 2nd neutron into "detector" coords.
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
C     Done.  Transfer to 2nd neutron variables.
      U2=Xn
      V2=Yn
      W2=Zn
C ---------------------------------------------------------------------
C     2ND NEUTRON
C ---------------------------------------------------------------------
      IKF   = 2112
      DDENG = En2
      DDVX  = U2
      DDVY  = V2
      DDVZ  = W2
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     and that's done.  Get energy of 9-Be.
      E9be=Ta-En2c
      V9Be=VELOCITY(K9Be,E9Be)
      Vbx=-Xn2*V9Be
      Vby=-Yn2*V9Be
      Vbz=-Zn2*V9Be
      CALL LABTRAN(Vbx,Vby,Vbz, V10Bel, Vx,Vy,Vz)
      V9Belab=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C ---------------------------------------------------------------------
C     9Be
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
      IKF   = 4*1000000 + 9
      DDENG = EFROMV(K9Be,V9Belab)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      Goto 90
C ---------------------------------------------------------------------
C     END SECTION ON 10-Be --> n + 9-Be

C     Finish case where 10-Be is final heavy ion.
   60 E10Be=EFROMV(K10Be,V10Bel)
C ---------------------------------------------------------------------
C     10Be
C ---------------------------------------------------------------------
      CALL DIRCOS(Vbx,Vby,Vbz,Xpn,Ypn,Zpn)
      Cx = Zx
      Cy = Zy
      Cz = Zz
      CALL TRANSVEC(Cx,Cy,Cz)
      IKF   = 4*1000000 + 10
      DDENG = E10Be
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
      Zx = Cx
      Zy = Cy
      Zz = Cz
C ---------------------------------------------------------------------
      Goto 90
C ---------------------------------------------------------------------
C     END OF THIS SECTION

C ---------------------------------------------------------------------
C     NEUTRON PLUS 11-B
C     Next section for 12-B --> n + 11-B
C ---------------------------------------------------------------------
   72 Ta=Excit-B11pn
      Encom=RCKE(Kn,K11B,Ta)
      Ehat=En-Qnp
      F=0.065 + 0.001*Ehat
      Temp=F*Ehat
C     Temp value sort of empirically determined so that CHOOSN gives a
C     reasonable representation of the neutron "continuum" as calculated
C     by the nuclear model code TNG.
      Eloww=0.0
C     Test if Eneut2 has a value left over from a previous "2n" reaction --
C     if so then force the present neutron emission to populate a state
C     in 11-B having Ex < binding energy of 10-B + n.
      IF (Eneut2 .LE. 0.0) goto 73
      Tryy = Encom-TenBpn
      IF (Tryy .GT. 0.0) Eloww=Tryy
   73 Enn=CHOOSN(Eloww,Encom,Temp)

C     AND As observed when choosing E(proton) from a distribution (as done above)
C     a "Virtual" energy level in 11-B.

      Try=Enn*(Emass(Kn)+Emass(K11B))/Emass(K11B)
      Entry=RCKE(Kn,K11B,Try)
      Try=Try*Enn/Entry
      Excitn=Ta-Try
C     If EXCITN is relatively small, pair it up with a known excited
C     state in 11-B.
      Egamma=0.0
      Eleftov=0.0
      IF (Excitn .LT. TenBpn) goto 75
      IF (Excitn .GT. TenBpn+0.5) goto 80
      Excitn=TenBpn+0.25
      K=15
      Goto 79
C     (That covers alpha emission in competition with neutron emission...)
   75 K=1
      Do 74 I=2,14
      IF (Excitn .LT. Exc11(I)) goto 76
   74 K=K+1
   76 Excitn=Exc11(K)
      Egamma=Exc11(K)
C     Most 11-B levels decay predominantly by ground-state transitions;
C     however, for three levels the dominant decay transition is through
C     an excited state.  In such cases, there are at least two gamma rays
C     to be followed, including the ground-state transition of
C     the intermediate excited state.
      IF (K.EQ.8 .OR. K.EQ.9) Egamma=Exc11(2)
      IF (K .EQ. 11) Egamma=Exc11(3)
      Eleftov=Exc11(K)-Egamma
C     That's the "other" gamma radiation.

   79 Taex=Ta-Excitn
      Enn=RCKE(Kn,K11B,Taex)

   80 CALL RVECT(Xp,Yp,Zp)
      Vnn=VELOCITY(Kn,Enn)
      Vncx=Xp*Vnn
      Vncy=Yp*Vnn
      Vncz=Zp*Vnn
      CALL LABTRAN(Vncx,Vncy,Vncz, V12B, Vnx,Vny,Vnz)
      Vneut=SQRT(Vnx*Vnx + Vny*Vny + Vnz*Vnz)
      Eneut=EFROMV(Kn,Vneut)
      SPDSQ=Eneut
C     Now get "new" neutron direction cosines in the
C     "detector" laboratory coordinates.
      CALL DIRCOS(Vnx,Vny,Vnz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
      U=Xn
      V=Yn
      W=Zn

C ---------------------------------------------------------------------
C     NEUTRON
C ---------------------------------------------------------------------
      IKF   = 2112
      DDENG = SPDSQ
      DDVX  = U
      DDVY  = V
      DDVZ  = W
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     and, finally, get the 11-B energy.

      Taex=Ta-Excitn
      E11B=RCKE(K11B,Kn,Taex)
      V11B=VELOCITY(K11B,E11B)
      Vbx=-Xp*V11B
      Vby=-Yp*V11B
      Vbz=-Zp*V11B
      CALL LABTRAN(Vbx,Vby,Vbz, V12B, Vx,Vy,Vz)
      V11B=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)

C     At this point want to look for second neutron emission,
C     i.e. (n,p2n) reaction.
C     First check -EXCITN- of 11-B for possible alpha emission.
      IF (Excitn .LE. Li7pa) goto 88
      IF (K .EQ. 10) goto 88
      Tryn=sngl(UNIRN(dummy))
      IF (K.EQ.11 .AND. Tryn.LT.0.9) goto 88
      IF (K.EQ.13 .AND. Tryn.LT.0.0001) goto 88
C ---------------------------------------------------------------------
C     If program counter gets to here, could be alpha or second
C     neutron emission.
C     First get moving 11-B ion into "detector" coordinates.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
C     Done.  Xn,Yn,Zn are dir. cosines of 11-B motion in "detector" coords.

      IF (K.EQ.15 .AND. Tryn.LT.0.56) goto 81
      IF (Excitn .GT. TenBpn) goto 82

C ---------------------------------------------------------------------
C     ALPHA PLUS 7-LI
C     Program counter to here, it's  11-B --> 7-Li + alpha
C ---------------------------------------------------------------------
   81 Egamma=0.0
      Eleftov=0.0
      Ta=Excitn-Li7pa
      CALL AP7LI(Ta,V11B,Egamma)
      Goto 90
C     END SECTION ON ALPHA PRODUCTION
C ---------------------------------------------------------------------

C ---------------------------------------------------------------------
C     NEUTRON PLUS 10-B
C     Next is case of second neutron emission, 11-B --> n + 10-B.
C ---------------------------------------------------------------------
   82 Ta=Excitn-TenBpn
      En2com=RCKE(Kn,K10B,Ta)
      Ehat=Ehat-B11pn
      Fn=0.065 + 0.001*Ehat
      Temp=Fn*Ehat
      Eloww=0.0
      En2=CHOOSN(Eloww,En2com,Temp)
      Try=En2*(Emass(Kn) + Emass(K10B))/Emass(K10B)
C     Iterate once for greater precision (since it costs very little)
      Entry=RCKE(Kn,K10B,Try)
      Try=Try*En2/Entry
      Excitn2=Ta-Try
      Egamma=0.0
      Eleftov=0.0
      Level=10
      IF (Excitn2 .GT. 6.11) goto 86
      K=1
C     Pair -EXCITN2- with 10-B energy level.
      Do 84 I=2,9
      IF (Excitn2 .LT. Ex10B(I)) goto 85
   84 K=K+1
   85 Excitn2=Ex10B(K)
      Level=K
C     See discussion following statement No. 18 of Subroutine NT for
C     gamma decay of levels in 10-B.
      IF (K .GT. 1) Egamma=Ex10B(2)
      Eleftov=Ex10B(K)-Egamma
C     Now get information on 2nd neutron into "detector" lab. coords.
   86 Taex=Ta-Excitn2
      Enn=RCKE(Kn,K10B,Taex)
      CALL RVECT(Xn2,Yn2,Zn2)
      Vn2=VELOCITY(Kn,Enn)
      V2x=Xn2*Vn2
      V2y=Yn2*Vn2
      V2z=Zn2*Vn2
      CALL LABTRAN(V2x,V2y,V2z, V11B, Vx,Vy,Vz)
      Vneut2=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      En2=EFROMV(Kn,Vneut2)
      Eneut2=-En2
c     The negative is temporary, just to alert BANKER routine.
c     Now get direction cosines of 2d neutron
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
C ---------------------------------------------------------------------
C     2ND NEUTRON
C ---------------------------------------------------------------------
      IKF   = 2112
      DDENG = En2
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
c     Done.  Xn,Yn,Zn are the needed direction cosines

c     Fill in rest of data for 2nd neutron.
      U2=Xn
      V2=Yn
      W2=Zn
C     Done with 2nd neutron; now get 10-B ion information.
      E10B=Taex-Enn
      V10B=VELOCITY(K10B,E10B)

      Vbx=-Xn2*V10B
      Vby=-Yn2*V10B
      Vbz=-Zn2*V10B + V11B
C     Again, non-relativistic transformation of c.o.m. --> lab. coords.

      V10Blab=SQRT(Vbx*Vbx + Vby*Vby + Vbz*Vbz)

C     Test for 10-B particle-emission decay:
      IF (Level .LE. 6) goto 87
      IF (Level.EQ.8 .AND. sngl(UNIRN(dummy)).LT. 0.5) goto 87

C     Program counter to here, prepare for 10-B particle decay.
      Eleftov=0.0
      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
C     All prepared.
C     (Xn,Yn,Zn) are the direction cosine of the 10B in the detector lab system.
      CALL TENBDK(Excitn2,V10Blab,Egamma)
C ---------------------------------------------------------------------
      Goto 90
C ---------------------------------------------------------------------

C ---------------------------------------------------------------------
C     Tidy up (n,p2n) results:
   87 E10B=EFROMV(K10B,V10Blab)
C ---------------------------------------------------------------------
C     10B
C ---------------------------------------------------------------------
      CALL DIRCOS(Vbx,Vby,Vbz,Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
      IKF   = 5*1000000 + 10
      DDENG = E10B
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      Goto 90
C     END THIS SECTION OF SECOND NEUTRON EMISSION
C ---------------------------------------------------------------------

C ---------------------------------------------------------------------
C     Tidy up (n,pn) results:
   88 E11B=EFROMV(K11B,V11B)
C ---------------------------------------------------------------------
C     11B
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
      IKF   = 5*1000000 + 11
      DDENG = E11B
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
   90 Nelm=7
      IDRTYP = Nelm

C     Now correct the Eneut2 variable.
      Eneut2=ABS(Eneut2)
C     Now check for possible gamma-ray contributions
      IF (Egamma .GT. 0.0) then
        CALL PHOTON(Egamma)
      ELSE IF (Eleftov .GT. 0.0) then
        CALL PHOTON(Eleftov)
      END IF

      CALL CX_anal

      Return
      END
C =====================================================================
C =====================================================================
      Subroutine N2N
C     Purpose is to compute the energies and directions of scatter
C     of the two neutrons and the energy of the 11-C ion from the
C     reaction  n + 12-C --> n' + n' + 11-C.
C
C     New features added 4/87: (1) fix excitation of residual 11-C to
C     be less than (a third) neutron separation energy; (2) allow
C     residual 11-C to decay by proton emission, if there is a
C     sufficient amount of excitation energy; and (3) in any case,
C     check for possible gamma-ray decay of particle-bound excited
C     states in either 11-C or 10-B, the latter if proton emission
C     is allowed.

C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------
      Common /NEUTR2/ Eneut2, U2,V2,W2
!$OMP THREADPRIVATE(/NEUTR2/)
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
!$OMP THREADPRIVATE(/VECTOR/)
      Common /NEUTRN/ SPDSQ, U, V, W
!$OMP THREADPRIVATE(/NEUTRN/)
C
      Common /MASSES/ Emass(21)
      Common /EXC10B/ Exc10(9)
C     Exc10 array in Subroutine NT
      Common /PPOLAR/ X(37),Pt(37),N37
      real*8 UNIRN

      Dimension Ptheta(37), Exc11C(11)

C     Next array is used to determine the polar scattering angle for
C     the first neutron.  It was taken from the formalism of the O5S
C     code, and represents a mildly forward-peaked scattering.  The
C     X-array which is its companion is given in Subroutine NP.
      Data Ptheta/0.0, 0.0332, 0.06605, 0.0986, 0.1309, 0.1628,
     A 0.1944, 0.2258, 0.2568, 0.2875, 0.3179, 0.348, 0.3778, 0.4073,
     B 0.4364, 0.4653, 0.4938, 0.5221, 0.550, 0.5776, 0.6049, 0.6319,
     C 0.6586, 0.685, 0.7111, 0.7369, 0.7624, 0.7875, 0.8124, 0.8369,
     D 0.8611, 0.885, 0.9086, 0.9319, 0.9549, 0.9776, 1.0/
C     Next array gives energies of particle-stable excited states in 11-C
      Data Exc11C/0.0, 2.0, 4.32, 4.8, 6.34, 6.48, 6.9, 7.5, 8.1, 8.42,
     a   8.66/
      Data Q/18.72/, Nterp/1/, TenCpn/13.124/, TenBpp/8.691/
      Data Kn/1/, Kp/2/, K12C/7/, K11C/9/, K10B/12/

C ---------------------------------------------------------------------
C     The reaction scheme for purposes of computation is:
C
C     n + 12-C --> n + 12-C;      then:
C     12-C --> n + 11-C;
C     We take E(neutron) from CHOOSN.
C
C     Addition 4/87.  Consider possibility of 11-C --> p + 10-B
C ---------------------------------------------------------------------

C     Start computation: go to Center-of-mass coordinates
      En=SPDSQ
      Vn = VELOCITY(Kn,En)
      CALL CMTRAN(Kn,K12C,Vn, Vcom,Eneut,Ecc)
      Tec=Ecc+Eneut
      Ta=Tec-Q
      Epl=0.0

      Egamma  = 0.0
      Eleftov = 0.0

      IF (Ta .GT. 0.0) goto 2

      write(*,555)En
  555 Format(/'   ***  Error in Subroutine N2N; E(neut) at entry ='
     b 1PE11.3/10x,'Set E(neut1) = E(neut2) = 0 and exit.'/)
      SPDSQ=0.0
      Eneut2=0.0
      Return

C     CHOOSE C.M. NEUTRON ENERGY
    2 Encom=RCKE(Kn,K12C,Ta)
C     Encom is the maximum energy the neutron can have in the center of mass
C     for the (n,2n) reaction.
      Eloww=0.0
      F=0.065+0.001*En
      Temp=F*En
C     Temp value chosen so that the "continuum" neutron spectrum calcu-
C     lated in CHOOSN is reasonably representative of that computed by
C     the nuclear model code TNG.
      Enc1 = CHOOSN(Eloww,Encom,Temp)

C     The difference in energy between Encom and Enc1 could be
C     interpreted as the added energy (larger than Q) of a "virtual"
C     energy level in the 12-C nucleus decaying by neutron emission.
      Try=Enc1*(Emass(Kn)+Emass(K12C))/Emass(K12C)
C     (That's approximately the "total available energy" that
C     would yield Enc1 from the function RCKE.)
      Entry=RCKE(Kn,K12C,Try)
      Try=Try*Enc1/Entry
C     (Iterate once on TRY for slightly more accuracy.)
      Excit=Ta-Try
C     (EXCIT  should be close to the added "virtual" excitation energy.)
      Qex=Q + Excit

C     FIND THE COSINE OF THE POLAR DIRECTION OF THE First Neutron (C.M.)
   15 RAND=sngl(UNIRN(dummy))
   18 CTHETA = EXTERP(Ptheta,X,Rand,N37,Nterp)
C     FIND THE SINE AND COSINE OF THE AZIMUTHAL DIRECTION OF THE
C     Neutron (C.M.) -- By Random-number Choice.
      RANN = 6.2831853*sngl(UNIRN(dummy))
      SINPHI = SIN(RANN)
      COSPHI = COS(RANN)

C     FIND THE C. M. VELOCITY COMPONENTS OF THE Neutron.
      Sintheta=SQRT(1.-Ctheta*Ctheta)
      Xp=Sintheta*Cosphi
      Yp=Sintheta*Sinphi
      Zp=Ctheta
      Vnr=VELOCITY(Kn,Enc1)
      Vnx=Xp*Vnr
      Vny=Yp*Vnr
      Vnz=Zp*Vnr
      CALL LABTRAN(Vnx,Vny,Vnz, Vcom, Vx,Vy,Vz)
      Vneut=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)

C     FIND THE LAB ENERGY OF THE Neutron.
      Enl=EFROMV(Kn,Vneut)
C     Store information
      SPDSQ=Enl
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=U
      Zy=V
      Zz=W
      CALL TRANSVEC(Zx,Zy,Zz)
      U=Xn
      V=Yn
      W=Zn
C ---------------------------------------------------------------------
C     NEUTRON
C ---------------------------------------------------------------------
      IKF   = 2112
      DDENG = SPDSQ
      DDVX  = U
      DDVY  = V
      DDVZ  = W
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

C     That takes care of neutron number one.
C
C     FIND THE COMPONENTS OF THE 12-C VELOCITY IN THE LAB.
C     First adjust total available energy for the energy of the
C     "virtual" excited state.
   70 Taex=Tec-Qex
      E12C=RCKE(K12C,Kn,Taex)
      V12C=VELOCITY(K12C,E12C)
      Vbcx=-Xp*V12C
      Vbcy=-Yp*V12C
      Vbcz=-Zp*V12C
C     (Recall that the direction cosines of the 12-C ion in the center
C     of mass are negatives of the direction cosines of the neutron.)
      CALL LABTRAN(Vbcx,Vbcy,Vbcz, Vcom, Vbx,Vby,Vbz)
      V12C=SQRT(Vbx*Vbx + Vby*Vby + Vbz*Vbz)

C     Now worry about the second neutron energy.  The maximum center-of-mass
C     kinetic energy available is given by the variable EXCIT where now the
C     center of mass is the 12-C ion, and its motion is the center-of-mass
C     motion.  Need to get motion of 12-C ion in detector coordinates.

      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
C     Now variables Xn,Yn,Zn = direction cosines of the motion of the
C     12-C in "detector" laboratory coordinate system.

      Ta=Excit
      Encom=RCKE(Kn,K11C,Ta)
      Ehat=En-Q
      F=0.065+0.001*Ehat
      Temp=F*Ehat

C     Now the added assumption -- following the neutron decay of the highly
C     excited 12-C nucleus, the residual 11-C ion is stable against further
C     neutron emission.  However, the possibility of proton decay of the
C     11-C to (as it turns out) a particle-bound state in 10-B will
C     be included.
      Tryy = Encom-TenCpn
      IF (Tryy .GT. 0.0) Eloww=Tryy
      Enn=CHOOSN(Eloww,Encom,Temp)

C     AND, once again, when choosing from a distribution the difference
C     in energy between Encom and Enn could be interpreted as a "Virtual"
C     energy level in 11-C.
      Try=Enn*(Emass(Kn)+Emass(K11C))/Emass(K11C)
      Entry=RCKE(Kn,K11C,Try)
      Try=Try*Enn/Entry
      Excitn=Ta-Try

C     If variable -Excitn- is small enough, pair it with a known particle-
C     stable state in 11-C
      Egamma=0.0
      IF (Excitn .GT. TenBpp) goto 9
      K=1
      Do 7 I=2,11
      IF (Excitn .LT. Exc11C(I)) goto 8
    7 K=K+1
    8 Excitn=Exc11C(K)
      Egamma=Excitn
      Taex=Ta-Excitn
      Enn=RCKE(Kn,K11C,Taex)
    9 Continue

C     SECOND NEUTRON EMISSION
C     For second neutron emission, the center-of-mass scattering is
C     isotropic.
      CALL RVECT(Xp,Yp,Zp)
      Vnn=VELOCITY(Kn,Enn)
      Vncx=Xp*Vnn
      Vncy=Yp*Vnn
      Vncz=Zp*Vnn
      CALL LABTRAN(Vncx,Vncy,Vncz, V12C, Vnx,Vny,Vnz)
      Vneut=SQRT(Vnx*Vnx + Vny*Vny + Vnz*Vnz)
      Eneut=EFROMV(Kn,Vneut)
      Eneut2=Eneut

C     Now get "new" neutron direction cosines in the
C     "detector" laboratory coordinates.
      CALL DIRCOS(Vnx,Vny,Vnz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
      U2=Xn
      V2=Yn
      W2=Zn
C ---------------------------------------------------------------------
C     2ND NEUTRON
C ---------------------------------------------------------------------
      IKF   = 2112
      DDENG = Eneut2
      DDVX  = U2
      DDVY  = V2
      DDVZ  = W2
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

C     now get the 11-C energy.
      Taex=Ta-Excitn
      E11C=RCKE(K11C,Kn,Taex)
      V11C=VELOCITY(K11C,E11C)
      Vbx=-Xp*V11C
      Vby=-Yp*V11C
      Vbz=-Zp*V11C+V12C
C     For Vbz including "LABTRAN" nonrelativistically.

C     Now have components of 11-C velocity in 12-C ion laboratory
C     coordinate system.
      V11C=SQRT(Vbx*Vbx + Vby*Vby + Vbz*Vbz)
      E11C=EFROMV(K11C,V11C)
      IF (Excitn .LE. TenBpp) goto 30

C ---------------------------------------------------------------------
C     PROTON EMISSION
C     If program gets to here, 11-C is energetic enough to decay by
C     proton emission (actually, also by alpha emission, but we don't
C     attempt that case) to particle-stable levels in 10-B.

C     First get 11-C motion into "detector" laboratory coordinate system.
      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
C     Done!  Xn,Yn,Zn = direction cosines of 11-C in "detector"
C     laboratory coordinate system.
      Ta=Excitn-TenBpp
      Epcom=RCKE(Kp,K10B,Ta)
      F=0.1245+0.001*ABS(En-45.0)
      Temp=F*En
      Eprot=CHOOSP(Epcom,Temp)
C     and as done before, one iteration for slightly improved accuracy..
      Try=Eprot*(Emass(Kp) + Emass(K10B))/Emass(K10B)
      Eptry=RCKE(Kp,K10B,Try)
      Try=Try*Eprot/Eptry
      Excit=Ta-Try
C     Now pair with an excited state in 10-B (see Function NT).
      Eleftov=0.0
      Egamma=0.0
      K=1
C     By present construction only 5 levels in 10-B are energetically
C     available, and they are all particle stable.
      Do 22 I=2,5
      IF (Excit .LT. Exc10(I)) goto 23
   22 K=K+1
   23 Excit=Exc10(K)

C     See comment following statement no. 18 in Subroutine NT for
C     gamma decay of levels in 10-B.
      IF (K .GT. 1) Egamma=Exc10(2)
      Eleftov=Exc10(K)-Egamma
      Taex=Ta-Excit
      Eprot=RCKE(Kp,K10B,Taex)

C     Get proton velocity in 11-C center-of-mass coordinates.
   25 Vprot=VELOCITY(Kp,Eprot)
C     Done.  "Fission" of 11-C --> p + 10-B is isotropic, perforce.
      CALL RVECT(Xp, Yp, Zp)
      Vpx=Xp*Vprot
      Vpy=Yp*Vprot
      Vpz=Zp*Vprot
      CALL LABTRAN(Vpx,Vpy,Vpz, V11C, Vx,Vy,Vz)
      Vprotl=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C     Vprotl=velocity of proton in 11-C laboratory coordinates.

      Epl=EFROMV(Kp,Vprotl)
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
C     Now Xn,Yn,Zn = direction cosines of proton in detector coordinates.
C ---------------------------------------------------------------------
C     PROTON
C ---------------------------------------------------------------------
      IKF   = 2212
      DDENG = Epl
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

   34 E10B=Taex-Eprot
      V10B=VELOCITY(K10B,E10B)
      Vbcx=-Xp*V10B
      Vbcy=-Yp*V10B
      Vbcz=-Zp*V10B +V11C
C     coordinate system.

      V10Blab=SQRT(Vbcx*Vbcx + Vbcy*Vbcy + Vbcz*Vbcz)
      E10B=EFROMV(K10B,V10Blab)
C ---------------------------------------------------------------------
C     10B
C ---------------------------------------------------------------------
      CALL DIRCOS(Vbcx,Vbcy,Vbcz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
      IKF   = 5*1000000 + 10
      DDENG = E10B
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      Goto 31
C ---------------------------------------------------------------------

C ---------------------------------------------------------------------
   30 continue
C ---------------------------------------------------------------------
C     11C
C ---------------------------------------------------------------------
      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
      IKF   = 6*1000000 + 11
      DDENG = E11C
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

   31 Nelm=8
      IDRTYP = Nelm

      IF (Egamma .GT. 0.0) then
         CALL PHOTON(Egamma)
      ELSE IF (Eleftov .GT. 0.0) then
        CALL PHOTON(Eleftov)
      END IF

      CALL CX_anal

      Return
      END
C =====================================================================
C =====================================================================
      Subroutine ND
C     Purpose to compute energies of the deuteron and 11-B coming from
C     the collision  n + 12-C --> d + 11-B.
C
C     One assumption is that the deuteron "continuum" has the
C     same relative energy distribution as computed for protons by
C     the nuclear model code TNG.  In addition, the programming
C     specifies some amount of ground-state (n,d) collisions for
C     all En, more or less consistent with measured spectra.
C
C     Added in 5/87: Consider possibility of 11-B --> n + 10-B
C     following the (n,d) reaction to a highly-excited "state" in 10-B.

C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------
      Common /NEUTRN/ Eneut, U,V,W
!$OMP THREADPRIVATE(/NEUTRN/)
      Common /MASSES/ Emass(21)
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
!$OMP THREADPRIVATE(/VECTOR/)
      Common /PPOLAR/ Xmu(37), Pth(37), Npolar
      Common /EXC11B/ Exc11(14),TenBpn,Li7pa
C     Exc11 array tabulated in subroutine NPX
      Common /EXC10B/ Ex10B(9)
C     Ex10B array tabulated in Subroutine NT
      real*8 UNIRN

      Dimension Ptheta(37)
C     The next array is used to determine the polar scattering
C     angle of the deuteron and is taken from the measurements
C     in Nuclear Instruments & Methods 129 (1975) 241 for
C     En(lab)=56 MeV.
      Data Ptheta/0.0, 0.013, 0.0596, 0.1318, 0.2128, 0.2959, 0.3768,
     Q  0.4527, 0.522, 0.584, 0.6387, 0.6861, 0.7267, 0.7611, 0.7901,
     P  0.8141, 0.8341, 0.8504, 0.8636, 0.8755, 0.8873, 0.8989,
     N  0.9103, 0.9212, 0.9318, 0.9418, 0.9513, 0.960, 0.9681,
     M  0.9753, 0.9817, 0.9872, 0.9918, 0.9954, 0.9979,0.9995, 1.0/
      Data Q/13.732/, Kn/1/, K12C/7/, Kd/10/, K11B/5/, Nterp/1/
      Data K10B/12/

C     Start computation.
      K  = 0
C     Go to center-of-mass coordinates.  Check to be sure incident
C     neutron energy is large enough.
      En=Eneut
      Eneut=0.0
      Vn=VELOCITY(Kn,En)
      CALL CMTRAN(Kn,K12C,Vn, Vcom,Enc,Ecc)
      Tec=Enc+Ecc
      Ta=Tec-Q

      IF (Ta .GT. 0.0) goto 2
      write(*,333)En
  333 Format(/'   ***  Error in Subroutine ND; E(neut) at entry = '
     a  1PE11.3/10x,'Set E(neut) = E(deuteron) = 0 and exit.'/)
      Return

C     Get center-of-mass energies, then transform to lab coord. energies
    2 Edcom=RCKE(Kd,K11B,Ta)
      Edeut=Edcom
      Excit=0.0
      Egamma=0.0
      Eleftov=0.0

      IF (Edcom .LT. 2.22) goto 9

      GSenhanc=0.00543*En
      IF (En .GT. 30.0) GSenhanc=0.08*SQRT(En-25.85)
      IF (sngl(UNIRN(dummy)) .LT. GSenhanc) goto 9
C     Last is to enhance ground-state transition, and is ad hoc.

      Temp=0.13*En
      IF (En.LT.30.0) Temp=12.6-1.07*En+0.026*En*En
C     Ad hoc formula for "nuclear temperature" to be used in CHOOSP
      Edeut=CHOOSP(Edcom,Temp)
C     Assume same "continuum" distribution of deuterons as we
C     have for protons.  It's an assumption.
      Try=Edeut*(Emass(Kd) + Emass(K11B))/Emass(K11B)
      Edtry=RCKE(Kd,K11B,Try)
C     Same "continuum" level excitation energy determination procedure
C     as used for protons in subroutine NPX.
      Try=Try*Edeut/Edtry
      Excit=Ta-Try
C     If variable Excit is relatively small, pair it with a known
C     excited state in 11-B.
      IF (Excit .LT. TenBpn)  goto 4
      IF (Excit .GT. TenBpn+0.5) goto 9
      Excit=TenBpn+0.25
      K=15
      Goto 8
    4 K=1
      Do 17 I=2,14
      IF (Excit .LT. Exc11(I)) goto 18
   17 K=K+1
   18 Excit=Exc11(K)
C     For gamma-ray energy assignments see comment following statement
C     No. 76 in Subroutine NPN
      Egamma=Exc11(K)
      IF (K.EQ.8 .OR. K.EQ.9) Egamma=Exc11(2)
      IF (K.EQ. 11) Egamma=Exc11(3)
      Eleftov=Exc11(K)-Egamma

    8 Taex=Ta-Excit
      Edeut=RCKE(Kd,K11B,Taex)

    9 Qex=Q + Excit

C     Next get polar scattering angle for the deuteron from array
C     Ptheta, above, and array Xmu (=array X in subroutine NP).
      Rand=sngl(UNIRN(dummy))
      Ctheta=EXTERP(Ptheta,Xmu,Rand,Npolar,Nterp)
      IF (Tec .GE. 51.7) goto 20
      Ciso=1.0-2.0*Rand
      Slope=(Ctheta-Ciso)/(51.7-Qex)
C     For En(c.o.m.) < 51.7 MeV interpolate between Ptheta results
C     and isotropy to get Ctheta.
      Ctheta=Slope*(Tec-Qex)
      Ctheta= EXTERP(Ptheta,Xmu,Rand,Npolar,Nterp)
C     Now get azimuthal scattering and direction cosines
   20 Rann = 6.283185*sngl(UNIRN(dummy))
      Sinphi=SIN(Rann)
      Cosphi=COS(Rann)
      Sintheta=SQRT(1.0 - Ctheta*Ctheta)
      Cx=Sintheta*Cosphi
      Cy=Sintheta*Sinphi
      Cz=Ctheta
      Vdeut=VELOCITY(Kd,Edeut)
      Vdx=Cx*Vdeut
      Vdy=Cy*Vdeut
      Vdz=Cz*Vdeut
      CALL LABTRAN(Vdx,Vdy,Vdz, Vcom, Vx,Vy,Vz)
      Vd=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Edt=EFROMV(Kd,Vd)
C     That's the energy for the deuteron.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)
C     Xn,Yn,Zn are deuteron's direction cosines in detector coords.
C ---------------------------------------------------------------------
C     DEUTERON
C ---------------------------------------------------------------------
      IKF   = 1000002
      DDENG = Edt
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

C     Next: get energy of recoil 11-B ion.
   10 Taex=Tec-Qex
      E11B=RCKE(K11B,Kd,Taex)
      V11B=VELOCITY(K11B,E11B)
      Vbx=-Cx*V11B
      Vby=-Cy*V11B
      Vbz=-Cz*V11B
      CALL LABTRAN(Vbx,Vby,Vbz, Vcom, Vx,Vy,Vz)
      V11B=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)

C     Okey to here for straight (n,d).  Now test for subsequent neutron decay,
C     i.e., (n,dn) reaction leading to bound states of 7-Li
      IF (Excit .LE. Li7pa) goto 40
      IF (K .EQ. 10) goto 40
      Tryn=sngl(UNIRN(dummy))
      IF (K.EQ.11 .AND. Tryn.LT.0.9) goto 40
      IF (K.EQ.13 .AND. Tryn.LT.0.0001) goto 40
C ---------------------------------------------------------------------

C ---------------------------------------------------------------------
C     Program counter to here, have alpha or neutron emission.
C     First get 11-B ion into "detector" coordinates.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)
C     Done.  Xn,Yn,Zn = dir. cosines of moving 11-B ion in "detector"
C     laboratory coordinates.
      IF (K.EQ.15 .AND. Tryn.LT.0.56) goto 14
      IF (Excit .GT. TenBpn) goto 15

C ---------------------------------------------------------------------
C     Program to here, have 11-B --> 7-Li + alpha decay
   14 Egamma=0.0
      Eleftov=0.0
      Ta=Excit-Li7pa
      CALL AP7LI(Ta,V11B,Egamma)
      Goto 45
C     END ALPHA DECAY SECTION
C ---------------------------------------------------------------------

C ---------------------------------------------------------------------
C     Program to here have 11-B --> n + 10-B decay to consider.
   15 Ta=Excit-TenBpn
      Encom=RCKE(Kn,K10B,Ta)
      Ehat=En-Q
      Fn=0.065 + 0.001*Ehat
      Temp=Fn*Ehat
C     Assumption of last three statements is that the "continuum"
C     neutron distribution from decay of the excited 11-B ion is the
C     same as used for decay of an excited 12-C ion after adjusting
C     the "incident" neutron energy.
      Eloww=0.0
      Enn=CHOOSN(Eloww,Encom,Temp)
      Try=Enn*(Emass(Kn) + Emass(K10B))/Emass(K10B)
C     Iterate once for improved accuracy on -Try-.
      Entry=RCKE(Kn,K10B,Try)
      Excitn=Ta-Try
      Egamma=0.0
      Eleftov=0.0
      Level=10
      IF (Excitn .GE. 6.11) goto 26
      K=1
      Do 24 I=2,9
      IF (Excitn .LT. Ex10B(I)) goto 25
   24 K=K+1
   25 Excitn=Ex10B(K)
      Level=K
C     For gamma decay of excited states of 10-B see comments following
C     Statement No. 18 in Subroutine NT
      IF (K.GT.1) Egamma=Ex10B(2)
      Eleftov=Ex10B(K)-Egamma

C     Now get information on neutron into "detector" laboratory coords.
   26 Taex=Ta-Excitn
      Enn=RCKE(Kn,K10B,Taex)
      CALL RVECT(Xnn,Ynn,Znn)
      Vnn=VELOCITY(Kn,Enn)
      Vxc=Xnn*Vnn
      Vyc=Ynn*Vnn
      Vzc=Znn*Vnn
      CALL LABTRAN(Vxc,Vyc,Vzc, V11B, Vx,Vy,Vz)
      Vneutl=Sqrt(Vx*Vx + Vy*Vy + Vz*Vz)
      Eneut=EFROMV(Kn,Vneutl)
C     That gets outgoing neutron energy.  Now get dir. cosines.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
C     Done.  Save in /NEUTRN/ Common variable locations.
      U=Xn
      V=Yn
      W=Zn
C ---------------------------------------------------------------------
C     NEUTRON
C ---------------------------------------------------------------------
      IKF   = 2112
      DDENG = Eneut
      DDVX  = U
      DDVY  = V
      DDVZ  = W
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

C     Now get information on 10-B ion.
      E10B=Taex-Enn
      V10B=VELOCITY(K10B,E10B)
      Vbx=-Xnn*V10B
      Vby=-Ynn*V10B
      Vbz=-Znn*V10B + V11B
C     Non-relativistic c.o.m. --> lab. coord. transformation for
C     heavy ion motion.
      V10Blab=SQRT(Vbx*Vbx + Vby*Vby + Vbz*Vbz)
C     Now test for possible (highly-excited) 10-B ion particle decay.
      IF (Level .LE. 6) goto 37
      IF (Level.EQ.8 .AND. sngl(UNIRN(dummy)).LT.0.5) goto 37

      Eleftov=0.0
C     Get 10-B ion dir. cosines in "detector" lab. coordinates
      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
C     Done.  Now do 10-B decay.
      CALL TENBDK(Excitn,V10Blab,Egamma)
      Goto 45

   37 E10B=EFROMV(K10B,V10Blab)
C ---------------------------------------------------------------------
C     10B
C ---------------------------------------------------------------------
      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
      IKF   = 5*1000000 + 10
      DDENG = E10B
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      Goto 45
C     END SECTION for neutron emission.
C ---------------------------------------------------------------------

C ---------------------------------------------------------------------
C     Next statement is the tag end of the (n,d) computation.
   40 EB=EFROMV(K11B,V11B)
C ---------------------------------------------------------------------
C     11B
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)
      IKF   = 5*1000000 + 11
      DDENG = EB
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     and that's it for the the 11-B ion

   45 Nelm=9
      IDRTYP = Nelm
C     Now check for possible gamma-ray contributions:
      IF (Egamma .GT. 0.0) then
        CALL PHOTON(Egamma)
      ELSE IF (Eleftov .GT. 0.0) then
        CALL PHOTON(Eleftov)
      END IF

      CALL CX_anal

      Return
      END
C =====================================================================
C =====================================================================
      Subroutine NT
C     Purpose to compute energies of the triton and 10-B ion coming from
C     the collision reaction  n + 12-C --> t + 10-B which (reaction)
C     has cross sections in the 30 mb region around En = 40 MeV.
C
C     If the c.o.m. total kinetic energy is not too large, the
C     program will fix the 10-B excitation at one of the known low-
C     lying energy levels.  The triton scattering angle (or more
C     precisely, the cosine of the scattering angle) is chosen by
C     random number from Legendre polynomial fits to data for En
C     between 27 and 61 MeV from Subramanian, et al, Physical
C     Review C28, 521.

C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------

      Common /NEUTRN/ Eneut, U,V,W
!$OMP THREADPRIVATE(/NEUTRN/)
      Common /MASSES/ Emass(21)
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
!$OMP THREADPRIVATE(/VECTOR/)
      Common /EXC10B/ Exc10(9)
      real*8 UNIRN

      Dimension Epoly(3),B1(3),B2(3),B3(3),B4(3),B5(3),B6(3),Fi(6)

C     Next array represents low-lying excited states of 10-B
      Data Exc10/0.0, 0.7183, 1.7402, 2.154, 3.587, 4.774,
     a   5.11, 5.17, 5.92/

C     Next arrays are for Legendre polynomials representing angular
C     distributions of tritons.  These distributions are forward
C     peaked at all En.
      Data Ipoly/3/, Epoly/27.4, 39.7, 60.7/, B1/0.42035, 0.3074,
     x  0.3681/, B2/0.209, 0.1459, 0.1619/, B3/0.0925, 0.1106, 0.1018/,
     y  B4/0.03296, 0.0228, 0.03733/, B5/0.01292, -0.01512, -0.00044/,
     z  B6/0.00382, -0.00692, -0.00318/
      Data Nfi/6/, Q/18.93/, Kn/1/, K12C/7/, Kt/13/, K10B/12/,Nterp/1/

C     START COMPUTATION.
C     Go to center of mass coordinates.
      En=Eneut
      Eneut=0.0
      Vn=VELOCITY(Kn,En)
      CALL CMTRAN(Kn,K12C,Vn, Vcom,Enc,Ecc)
      Tec=Enc+Ecc
      Ta=Tec-Q
      IF (Ta .GT. 0.0) goto 2

      write(*,444)En
  444 Format(/'   ***  Error in Subroutine NT; E(neut) at entry = '
     a  1PE11.3/10x,'Set E(neut) = E(Triton) = 0 and exit.'/)
      Return
C     Last is an error return.

C     Get center-of-mass energies, then transform to lab coordinate
C     energies.  First set up Legendre coeff. for En.
    2 Fi(1)=EXTERP(Epoly,B1,En,Ipoly,Nterp)
      Fi(2)=EXTERP(Epoly,B2,En,Ipoly,Nterp)
      Fi(3)=EXTERP(Epoly,B3,En,Ipoly,Nterp)
      Fi(4)=EXTERP(Epoly,B4,En,Ipoly,Nterp)
      Fi(5)=EXTERP(Epoly,B5,En,Ipoly,Nterp)
      Fi(6)=EXTERP(Epoly,B6,En,Ipoly,Nterp)
      Etcom=RCKE(Kt,K10B,Ta)
      Etrit=Etcom
      Excit=0.0
      Egamma=0.0
      Eleftov=0.0
      Level=1
      IF (Etcom .LT. 0.72) goto 9
      IF (sngl(UNIRN(dummy)) .LT. 0.05) goto 9
C     Last is to enhance slightly the ground-state transition,
C     and is ad hoc.

      Level=10
      Temp=0.5+0.04*En
C     Assume same "continuum" distribution of tritons as is used for
C     protons.  It's an assumption!
      Etrit=CHOOSP(Etcom,Temp)
      Try=Etrit*(Emass(Kt) + Emass(K10B))/Emass(K10B)
      Etry=RCKE(Kt,K10B,Try)

C     Same level excitation energy determination procedure as used
C     for deuterons in subroutine ND.
      Try=Try*Etrit/Etry
      Excit=Ta-Try
C     If variable Excit is relatively small, pair it with a known
C     excited state in 10-B
      IF (Excit .GT. 6.11) goto 9
      K=1
      Do 17 I=2,9
      IF (Excit .LT. Exc10(I)) goto 18
   17 K=K+1

   18 Excit=Exc10(K)
      Level=K

C     Gamma decay of levels in 10-B strongly proceed via the first-excited
C     state of 10-B.  This means at least two photons per decay of an
C     excited state, one of which corresponds to the ground-state
C     transition of the first excited state.
      IF (K .GT. 1) Egamma=Exc10(2)
      Eleftov=Exc10(K)-Egamma
      Taex=Ta-Excit
      Etrit=RCKE(Kt,K10B,Taex)

    9 Qex=Q + Excit

C     Next get polar scattering angle for the triton from array Fi
      Ctheta=CHOOSL(Fi,Nfi)
C     Now get azimuthal scattering angle.
   20 Rann = 6.283185*sngl(UNIRN(dummy))
      Sinphi=SIN(Rann)
      Cosphi=COS(Rann)
      Sintheta=SQRT(1.0 - Ctheta*Ctheta)
      Cx=Sintheta*Cosphi
      Cy=Sintheta*Sinphi
      Cz=Ctheta

      Vtri=VELOCITY(Kt,Etrit)
      Vtx=Cx*Vtri
      Vty=Cy*Vtri
      Vtz=Cz*Vtri
      CALL LABTRAN(Vtx,Vty,Vtz, Vcom, Vx,Vy,Vz)
      Vt=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Ett=EFROMV(Kt,Vt)
C ---------------------------------------------------------------------
C     TRITON
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)

      Call TRANSVEC(U,V,W)

      IKF   = 1000003
      DDENG = Ett
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn

      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     That's the energy for the triton.

C     So next get the energy of the recoil 10-B ion.
   10 Taex=Tec-Qex
      E10B=RCKE(K10B,Kt,Taex)
      V10B=VELOCITY(K10B,E10B)
      Vbx=-Cx*V10B
      Vby=-Cy*V10B
      Vbz=-Cz*V10B + Vcom
C     Last transforms to "neutron" Lab. coords. non-relativistically.

      V10B=SQRT(Vbx*Vbx + Vby*Vby + Vbz*Vbz)

C     Test for possible (highly-excited) 10-B ion particle decay.
      IF (Level .LE. 6) goto 12
      IF (Level.EQ.8 .AND. sngl(UNIRN(dummy)).LT.0.5) goto 12

      Eleftov=0.0
C     Get 10-B ion direction cosines in "detector" lab. coordinates.
      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)
C     Done.  Now do decay of highly excited 10-B ion.

      CALL TENBDK(Excit,V10B,Egamma)
      Goto 16
C ---------------------------------------------------------------------
   12 E10B=EFROMV(K10B,V10B)
C ---------------------------------------------------------------------
C     10B
C ---------------------------------------------------------------------
      CALL DIRCOS(Vbx,Vby,Vbz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(U,V,W)
      IKF   = 5*1000000 + 10
      DDENG = E10B
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     and that's it for the 10-B ion.
C ---------------------------------------------------------------------

   16 Nelm=11
      IDRTYP = Nelm

C     Now check for possible gamma-ray contribution.

      IF (Egamma .GT. 0.0) then
        CALL PHOTON(Egamma)
      ELSE IF (Eleftov .GT. 0.0) then
        CALL PHOTON(Eleftov)
      END IF

      CALL CX_anal

      Return
      End
C =====================================================================
C =====================================================================
      Subroutine N3HE
C     Purpose is to compute energetics for reaction
C     n + 12-C --> 3-He + 10-Be, with possible subsequent decay
C     10-Be --> n + 9-Be.            As of 6/87.
C     Added in 9-Be decay modes -- 8/87.

C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------
      Common /NEUTRN/ Eneut, U,V,W
!$OMP THREADPRIVATE(/NEUTRN/)
      Common /MASSES/ Emass(21)
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
!$OMP THREADPRIVATE(/VECTOR/)

      Common /EXCTBE/ Ex10Be(4),Be9pn
C     Array Ex10Be is in subroutine NPX.

      Common /NEUTR2/ Eneut2, U2,V2,W2
!$OMP THREADPRIVATE(/NEUTR2/)
      Common /PPLI8/  Nx9,Ex9(14),Pp(14)
      Common /BENINE/ Exc9Be(6),Be8pn
C     See subroutine NN3ALF for values of variables in -PPLI8- common.
      real*8 UNIRN

      Dimension Epoly(3),B1(3),B2(3),B3(3),B4(3),B5(3),B6(3),Fi(6)

      Data Kn/1/, K12C/7/, K3He/19/, K10Be/15/, K9Be/4/, K8Be/8/
      Data Q/19.47/, Nfi/6/, Be8pn/1.6652/, Be8gs/0.092/, Ka/3/

C     Next array represents 9-Be level structure well enough for present use.
      Data Exc9Be/0.0, 2.9, 4.7, 6.8, 11.5, 14.0/

C     Next arrays will be used to determine angular distributions of
C     outgoing 3-He ions using Legendre coefficients determined from
C     study of the data provided by Subramanian, et al, Physical
C     Review C28, 521 (1983).
      Data Ipoly/3/, Epoly/27.4, 39.7, 60.7/
      Data B1/0.4204, 0.4207, 0.3993/,  B2/0.209, 0.1761, 0.1569/
      Data B3/0.0925, 0.046,  0.09874/, B4/0.033, 0.0234, 0.04945/
      Data B5/0.001,  0.00243,0.01223/, B6/0.004, 0.01271,0.00558/

C     Start computation.  Initialize variables.
      En2=0.0

C     Go to center-of-mass and check E(neutron)
      En=Eneut
      Eneut=0.0
      Vn=VELOCITY(Kn,En)
      CALL CMTRAN(Kn,K12C,Vn, Vcom,Enc,Ecc)
      Tec=Enc+Ecc
      Ta=Tec-Q
      IF (Ta .GT. 0.0) goto 2

      write(*,100)En
  100 Format(/'   ***  Error in Subroutine N3HE, En = '1PE11.4/
     x  10x,'Set E(neutron) = 0.0 and Exit... ')
      Return
C     Last is an "error" return.

C     Set up Legendre coeff. computation.
    2 Fi(1)=EXTERP(Epoly,B1,En,Ipoly,1)
      Fi(2)=EXTERP(Epoly,B2,En,Ipoly,1)
      Fi(3)=EXTERP(Epoly,B3,En,Ipoly,1)
      Fi(4)=EXTERP(Epoly,B4,En,Ipoly,1)
      Fi(5)=EXTERP(Epoly,B5,En,Ipoly,1)
      Fi(6)=EXTERP(Epoly,B6,En,Ipoly,1)

C     Get center-of-mass energies:
      Ehcom=RCKE(K3He,K10Be,Ta)
      Ehe=Ehcom
      Excit=0.0
      Egamma=0.0
      Eleftov=0.0
      Level=1
      IF (Ehcom .LT. Ex10Be(2)) goto 8
      IF (sngl(UNIRN(dummy)) .LT. 0.04) goto 8
C     Last is to slightly enhance ground-state population of 10-Be.

      Level=5
      Temp=3.0
C     Last is an assumption, and probably not very critical since the
C     3-He production is small.

      Ehe=CHOOSP(Ehcom,Temp)
      Try=Ehe*(Emass(K3He) + Emass(K10Be))/Emass(K10Be)
      Etry=RCKE(K3He,K10Be,Try)
      Try=Try*Ehe/Etry

C     As usual, one iteration on -Try- for improved precision.
      Excit=Ta-Try
      IF (Excit .GT. Be9pn) goto 6

C     If -Excit- is < 9-Be + n, then pair -Excit- with "level" in 10-Be.
      K=1
      Do 3 I=2,4
      IF (Excit .LT. Ex10Be(I)) goto 4
    3 K=K+1
    4 Excit=Ex10Be(K)
      Level=K

C     See subroutine NPX, comments following statement no. 42, for decay
C     characteristics of 10-Be levels.

      IF (K.GE.2) Egamma=Ex10Be(2)
      IF (K.EQ.3 .AND. sngl(UNIRN(dummy)).GE.0.5) Egamma=Excit
      Eleftov=Excit-Egamma

    6 Taex=Ta-Excit
      Ehe=RCKE(K3He,K10Be,Taex)

    8 Qex=Q+Excit

C     Now get 3-He scattering angle by random choice from Legendre
C     coefficients tabulated above.
      Ctheta=CHOOSL(Fi,Nfi)
C     Done.  Now choose azimuthal angle and get direction cosines of 3-He.

      Rann=6.283185*sngl(UNIRN(dummy))
      Sinphi=SIN(Rann)
      Cosphi=COS(Rann)
      Sintheta=SQRT(1.0-Ctheta*Ctheta)
      Cx=Sintheta*Cosphi
      Cy=Sintheta*Sinphi
      Cz=Ctheta
C     Done.  Get velocity components of 3-He ion in center-of-mass coords.

      Vhe=VELOCITY(K3He,Ehe)
      Vhx=Cx*Vhe
      Vhy=Cy*Vhe
      Vhz=Cz*Vhe
C     Done.  Transform to velocity components in laboratory coords.

      CALL LABTRAN(Vhx,Vhy,Vhz, Vcom, Vx,Vy,Vz)
C     Done.  Get Velocity and Energy of 3-He ion in laboratory coords.

      V3He=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      E3He=EFROMV(K3He,V3He)
C ---------------------------------------------------------------------
C     3He
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      Zx = U
      Zy = V
      Zz = W
      CALL TRANSVEC(Zx,Zy,Zz)
      IKF   = 2000003
      DDENG = E3He
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     That's the energy of the 3-He ion.  Now study 10-Be ion motion.
      Taex=Tec-Qex
      EBe=RCKE(K10Be,K3He,Taex)
      VBe=VELOCITY(K10Be,EBe)
      Vbx=-Cx*VBe
      Vby=-Cy*VBe
      Vbz=-Cz*VBe
C     (Recall that the direction cosines of the 10-Be ion in the center
C     of mass are negatives of the direction cosines of the 3-He ion.)
      CALL LABTRAN(Vbx,Vby,Vbz, Vcom, Vx,Vy,Vz)
      V10Be=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C     Check for 10-Be decay.  If not, save energy of 10-Be ion and exit.

      IF (Level .LE. 4) goto 20
C ---------------------------------------------------------------------

C ---------------------------------------------------------------------
C     Program counter to here, study 10-Be --> n + 9-Be.
C
C     First rotate 10-Be motion into "detector" lab. coordinates.
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
C     Done.  Dir. cosines of 10-Be motion are Xn,Yn,Zn in "detector"
C     laboratory coordinate system.

      Ta=Excit-Be9pn
      Enmax=RCKE(Kn,K9Be,Ta)
      Enn=Enmax
      Ehat=En-Q
      Fn=0.065+0.001*Ehat
      Temp=Fn*Ehat
      Elow=0.0
      Entry=CHOOSN(Elow,Enmax,Temp)
      Try=Entry*(Emass(Kn) + Emass(K9Be))/Emass(K9Be)
      Etry=RCKE(Kn,K9Be,Try)
      Try=Try*Entry/Etry
      Excit=Ta-Try

C     At this point check for possible proton decay.
      Npgo=0
      IF (Excit .LE. Ex9(1)) goto 10
      Prob=EXTERP(Ex9,Pp,Excit,Nx9,1)
      IF (sngl(UNIRN(dummy)) .GT. Prob) goto 10

      Npgo=1
C     Program to here, have  9-Be --> p + 8-Li
      Goto 14

C     Next check for 2d neutron already waiting to be processed.  If
C     go directly to ground state of 9-Be.

   10 IF (Eneut2 .GT. 0.0) goto 15

C ---------------------------------------------------------------------
C     Program counter to here, then have  9-Be --> n + 8-Be
C     followed by                         8-Be --> 2 alphas.
C
C     Restrict 8-Be to ground state for simplicity -- justification is
C     small cross sections just to get this far, so won't be very
C     many of these events.
      Npgo=-1
      IF (Excit .GE. 16.0) goto 14
      K=1
      Do 12 I=2,6
      IF (Excit .LT. Exc9Be(I)) goto 13
   12 K=K+1

   13 Excit=Exc9Be(K)

C     Check for  -EXCIT-  =0, i.e., ground state of 9-Be.
      IF (K .EQ. 1) Npgo=0
      Try=Ta-Excit
   14 Enn=RCKE(Kn,K9Be,Try)
   15 Vn=VELOCITY(Kn,Enn)

C     Isotropy of  n + 9-Be  in 10-Be center-of-mass coordinates.
      CALL RVECT(Xp,Yp,Zp)
      Vnx=Xp*Vn
      Vny=Yp*Vn
      Vnz=Zp*Vn

C     Transform velocity components to 10-Be lab. coords.
      CALL LABTRAN(Vnx,Vny,Vnz, V10Be, Vx,Vy,Vz)
C     Done.  Get velocity, energy of neutron in 10-Be lab. coords.

      Vnlab=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Enlab=EFROMV(Kn,Vnlab)
      Eneut=Enlab
C     Done.  Rotate neutron velocity into "detector" lab. coords.

      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Z1=Xn
      ZB=Yn
      Z3=Zn
      CALL TRANSVEC(Z1,ZB,Z3)
C ---------------------------------------------------------------------
C     NEUTRON
C ---------------------------------------------------------------------
      IKF   = 2112
      DDENG = Eneut
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     Done.  Put new dir. cosines into /NEUTRN/ Common variables.
      U=Xn
      V=Yn
      W=Zn
C     (Same position coordinates as for  n + 12-C assumed.)

C     Now get information on energy of 9-Be to save.
      IF (Npgo .EQ. 0) Try=Ta
      E9Bec=Try-Enn
      V9Bec=VELOCITY(K9Be,E9Bec)
      V9x=-Xp*V9Bec
      V9y=-Yp*V9Bec
      V9z=-Zp*V9Bec

C     Transform those c.o.m. components to 10-Be laboratory coordinates.
      CALL LABTRAN(V9x,V9y,V9z, V10Be, Vx,Vy,Vz)
C     Done.  Get energy of 9-Be.

      V9Be=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)

C     Check for further neutron or proton decay of 9-Be
      IF (Npgo .EQ. 0) goto 19

C ---------------------------------------------------------------------
C     Program to here, further particle decay.  First rotate 9-Be
C     velocity vector to "detector" coordinates.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Z1,ZB,Z3)
C     Done.  Set up and test for proton decay of 9-Be.
C     (Xn,Yn,Zn) are dir cos of 9-Be in detector lab sys

      Egamma=0.0
      Eleftov=0.0
      IF (Npgo .EQ. 1) goto 16

C     Program to here, it's  9-Be --> n + 8-Be.
      Taex=Excit-Be8gs-Be8pn
      Enc=RCKE(Kn,K8Be,Taex)
      Vnc=VELOCITY(Kn,Enc)
      CALL RVECT(Xp,Yp,Zp)
      Vnx=Xp*Vnc
      Vny=Yp*Vnc
      Vnz=Zp*Vnc
      CALL LABTRAN(Vnx,Vny,Vnz, V9Be, Vx,Vy,Vz)
      Vnlab=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Eneut2=EFROMV(Kn,Vnlab)

C     Rotate neutron velocity into "detector" lab. coordinates.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Z1=Xn
      ZB=Yn
      Z3=Zn
      CALL TRANSVEC(Z1,ZB,Z3)
C ---------------------------------------------------------------------
C     2ND NEUTRON
C ---------------------------------------------------------------------
      IKF   = 2112
      DDENG = Eneut2
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     Done.  Put values into  /NEUTR2/  common variables.
      U2=Xn
      V2=Yn
      W2=Zn

C     Now get 8-Be ion motion.
      E8Bec=Taex-Enc
      V8Bec=VELOCITY(K8Be,E8Bec)
      V8x=-Xp*V8Bec
      V8y=-Yp*V8Bec
      V8z=-Zp*V8Bec
      CALL LABTRAN(V8x,V8y,V8z, V9Be, Vx,Vy,Vz)
      V8Be=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C     Done.
C ---------------------------------------------------------------------
C     Rotate 8-Be ion velocity into "detector" lab. coordinates.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Z1,ZB,Z3)
      C1 = Xn
      CB = Yn
      C3 = Zn
C     (C1,CB,C3) are the direction cosine of the 8-Be ion in the detector lab. sys.
C ---------------------------------------------------------------------

C     Now  8-Be --> 2 alphas.
      Ta=Be8gs
      Ea=0.5*Ta
      Vac=VELOCITY(Ka,Ea)
      CALL RVECT(Xp,Yp,Zp)
      Vax=Xp*Vac
      Vay=Yp*Vac
      Vaz=Zp*Vac
      CALL LABTRAN(Vax,Vay,Vaz, V8Be, Vx,Vy,Vz)
      Valph=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Ealph=EFROMV(Ka,Valph)
C ---------------------------------------------------------------------
C     1ST ALPHA
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      CALL TRANSVEC(C1,CB,C3)
      IKF   = 2000004
      DDENG = Ealph
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C    That's one alpha.  Now get the other.
      Vax2=-Vax
C    (By symmetry.)
      Vay2=-Vay
      Vaz2=-Vaz
      CALL LABTRAN(Vax2,Vay2,Vaz2, V8Be, Vx,Vy,Vz)
      Valph2=SQRT(Vx*Vx + Vy*Vy +Vz*Vz)
      Ealph2=EFROMV(Ka,Valph2)
C ---------------------------------------------------------------------
C     2ND ALPHA
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      CALL TRANSVEC(C1,CB,C3)
      IKF   = 2000004
      DDENG = Ealph2
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     Done!!
      Goto 25
C ---------------------------------------------------------------------

C ---------------------------------------------------------------------
C     Next is proton decay, i.e.,  9-Be --> p + 8-Li
C     A complexity at this point.  If there is a second neutron
C     waiting to be processed, then Subroutine PP8LI must be
C     alerted to inhibit further neutron decay of 8-Li.  In
C     addition, Subroutine ELI8DK (which is called by PP8LI)
C     puts the information on the neutron from  8-Li --> n + 7-Li
C     into the /NEUTRN/ common; however, this common area already
C     has information for the neutron from  10-Be --> n + 9-Be
C     reaction just completed.  So need to do some preparation.

   16 V9=V9Be
      IF (Eneut2 .GT. 0.0) goto 17

C     If program counter to here, no second neutron waiting to be
C     processed, so shift first neutron information into /NEUTR2/ common.
      Eneut2 = Eneut
      U2 = U
      V2 = V
      W2 = W
      Eneut = 0.0
      I2=0
      Goto 18

C     But if there is a second neutron, alert -PP8LI- through variable V9Be.
   17 V9=-V9Be
      I2=1
   18 CALL PP8LI(Excit,V9,Egamma)

C     If both neutrons determined by present routine, set Eneut2 to
C     negative to alert -BANKR2- ; will return to >0 before exit from
C     routine.
      IF (I2 .EQ. 0) Eneut2=-Eneut2
      Goto 25
C ---------------------------------------------------------------------

C ---------------------------------------------------------------------
   19 E9Be=EFROMV(K9Be,V9Be)
C ---------------------------------------------------------------------
C     9Be
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Z1,ZB,Z3)
      IKF   = 4*1000000 + 9
      DDENG = E9Be
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     That's it for 9-Be ground-state results!
      Goto 25
C     END  n + 9-Be SECTION.
C ---------------------------------------------------------------------

C ---------------------------------------------------------------------
C     Tidy up 10-Be ion.
   20 E10Be=EFROMV(K10Be,V10Be)
C ---------------------------------------------------------------------
C     10Be
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
C     Now variables Xn,Yn,Zn = direction cosines of the motion of the
C     10-Be in "detector" laboratory coordinate system.
      IKF   = 4*1000000 + 10
      DDENG = E10Be
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

   25 Nelm=12
      IDRTYP = Nelm

      IF (Eneut2 .LT. 0.0) Eneut2=-Eneut2

C     Check for gamma ray interaction.
      IF (Egamma .GT. 0.0) then
        CALL PHOTON(Egamma)
      ELSE IF (Eleftov.GT. 0.0) then
        CALL PHOTON(Eleftov)
      END IF

      CALL CX_anal

      Return
      END
C =====================================================================
C =====================================================================
      BLOCK DATA LEGNDR
C     These are block data for  n + 12-C  elastic scattering angular
C     distributions.  For E(neut) (lab) up to 20 MeV these are data
C     in ENDF/B-V.  For the three higher energies, the values are
C     taken from fits to data exhibited in Phys. Med. Bio. 29 (1984) 643.

      COMMON /LGNDRE/ IC, E(220), F(6,220)

      DATA IC /220/

      DATA E /
     a.10000E-10,.10000E-02,.50000E-02,.10000E-01,.50000E-01,.10000E+00,
     a.20000E+00,.30000E+00,.40000E+00,.50000E+00,.60000E+00,.70000E+00,
     a.80000E+00,.90000E+00,.10000E+01,.11000E+01,.12000E+01,.13000E+01,
     a.14000E+01,.15000E+01,.16000E+01,.17000E+01,.18000E+01,.18500E+01,
     a.19000E+01,.19250E+01,.19500E+01,.19600E+01,.19700E+01,.19800E+01,
     a.19900E+01,.20000E+01,.20200E+01,.20400E+01,.20600E+01,.20700E+01,
     a.20800E+01,.20900E+01,.21000E+01,.21100E+01,.21200E+01,.21300E+01,
     a.21400E+01,.21500E+01,.23000E+01,.24500E+01,.26000E+01,.27500E+01,
     a.27800E+01,.28000E+01,.28050E+01,.28100E+01,.28150E+01,.28200E+01,
     a.28300E+01,.28400E+01,.28800E+01,.29000E+01,.29200E+01,.29400E+01,
     a.29600E+01,.29800E+01,.30000E+01,.30200E+01,.30400E+01,.30800E+01,
     a.31200E+01,.32000E+01,.33000E+01,.34000E+01,.35000E+01,.36000E+01,
     a.37000E+01,.38000E+01,.39000E+01,.40000E+01,.40500E+01,.41000E+01,
     a.41500E+01,.42000E+01,.42500E+01,.43000E+01,.43500E+01,.44000E+01,
     a.44500E+01,.45000E+01,.46000E+01,.47000E+01,.48000E+01,.48800E+01,
     a.49000E+01,.49200E+01,.49250E+01,.49300E+01,.49350E+01,.49400E+01,
     a.49600E+01,.49800E+01,.50000E+01,.52000E+01,.52500E+01,.52850E+01,
     a.53000E+01,.53200E+01,.53250E+01,.53500E+01,.53600E+01,.53650E+01,
     a.53800E+01,.54000E+01,.54200E+01,.54700E+01,.55000E+01,.55500E+01,
     a.57000E+01,.59500E+01,.60500E+01,.61000E+01,.61500E+01,.61850E+01,
     a.62000E+01,.62100E+01,.62250E+01,.62300E+01,.62500E+01,.62600E+01,
     a.62750E+01,.62850E+01,.63000E+01,.63100E+01,.63200E+01,.63250E+01,
     a.63400E+01,.63500E+01,.63600E+01,.63900E+01,.64000E+01,.64300E+01,
     a.64500E+01,.64600E+01,.64700E+01,.64800E+01,.64900E+01,.65000E+01,
     a.65300E+01,.65400E+01,.65500E+01,.65700E+01,.66000E+01,.66100E+01,
     a.66200E+01,.66500E+01,.66800E+01,.67000E+01,.67300E+01,.68000E+01,
     a.68500E+01,.69500E+01,.69700E+01,.70000E+01,.70500E+01,.71100E+01,
     a.71300E+01,.71500E+01,.71600E+01,.71700E+01,.72000E+01,.72500E+01,
     a.73000E+01,.73250E+01,.73500E+01,.73900E+01,.74000E+01,.74400E+01,
     a.74500E+01,.75000E+01,.75500E+01,.75800E+01,.76000E+01,.76350E+01,
     a.76700E+01,.77000E+01,.77400E+01,.77500E+01,.77700E+01,.77800E+01,
     a.77900E+01,.78000E+01,.78100E+01,.78300E+01,.78500E+01,.80000E+01,
     a.81500E+01,.82000E+01,.82500E+01,.83500E+01,.84500E+01,.85000E+01,
     a.87050E+01,.90000E+01,.95000E+01,.10000E+02,.10500E+02,.11000E+02,
     a.11500E+02,.12000E+02,.12500E+02,.13000E+02,.13500E+02,.14000E+02,
     a.14500E+02,.15000E+02,.16000E+02,.17000E+02,.18000E+02,.19000E+02,
     a.20000E+02, 20.8, 26.0, 40.0/

      DATA ((F(I,J),I=1,6),J=1,25)/
     a 0.00000E+00, 0.00000E+00, 0.00000E+00,
     a 0.00000E+00, 0.00000E+00, 0.00000E+00,
     a 0.14011E-03, 0.00000E+00, 0.00000E+00,
     a 0.00000E+00, 0.00000E+00, 0.00000E+00,
     a 0.69820E-03, 0.00000E+00, 0.00000E+00,
     a 0.00000E+00, 0.00000E+00, 0.00000E+00,
     a 0.13906E-02, 0.00000E+00, 0.00000E+00,
     a 0.00000E+00, 0.00000E+00, 0.00000E+00,
     a 0.67282E-02, 0.74988E-04, 0.00000E+00,
     a 0.00000E+00, 0.00000E+00, 0.00000E+00,
     a 0.12923E-01, 0.27931E-03, 0.00000E+00,
     a 0.00000E+00, 0.00000E+00, 0.00000E+00,
     a 0.23883E-01, 0.97357E-03, 0.00000E+00,
     a 0.00000E+00, 0.00000E+00, 0.00000E+00,
     a 0.33165E-01, 0.19172E-02, 0.63380E-04,
     a 0.00000E+00, 0.00000E+00, 0.00000E+00,
     a 0.40990E-01, 0.29952E-02, 0.12848E-03,
     a 0.00000E+00, 0.00000E+00, 0.00000E+00,
     a 0.47529E-01, 0.41298E-02, 0.22415E-03,
     a 0.00000E+00, 0.00000E+00, 0.00000E+00,
     a 0.52916E-01, 0.52714E-02, 0.34924E-03,
     a 0.00000E+00, 0.00000E+00, 0.00000E+00,
     a 0.57252E-01, 0.63923E-02, 0.49974E-03,
     a-0.57108E-04, 0.00000E+00, 0.00000E+00,
     a 0.60616E-01, 0.74832E-02, 0.66841E-03,
     a-0.10257E-03, 0.00000E+00, 0.00000E+00,
     a 0.63065E-01, 0.85502E-02, 0.84457E-03,
     a-0.17208E-03, 0.00000E+00, 0.00000E+00,
     a 0.64638E-01, 0.96140E-02, 0.10135E-02,
     a-0.27366E-03, 0.00000E+00, 0.00000E+00,
     a 0.65355E-01, 0.10710E-01, 0.11562E-02,
     a-0.41675E-03, 0.00000E+00, 0.00000E+00,
     a 0.65220E-01, 0.11887E-01, 0.12476E-02,
     a-0.61223E-03, 0.00000E+00, 0.00000E+00,
     a 0.64218E-01, 0.13213E-01, 0.12560E-02,
     a-0.87217E-03, 0.00000E+00, 0.00000E+00,
     a 0.62315E-01, 0.14775E-01, 0.11394E-02,
     a-0.12092E-02, 0.61831E-04, 0.00000E+00,
     a 0.59450E-01, 0.16680E-01, 0.84180E-03,
     a-0.16345E-02, 0.89279E-04, 0.00000E+00,
     a 0.55528E-01, 0.19063E-01, 0.28256E-03,
     a-0.21536E-02, 0.12619E-03, 0.00000E+00,
     a 0.50385E-01, 0.22081E-01,-0.66742E-03,
     a-0.27532E-02, 0.17476E-03, 0.00000E+00,
     a 0.43701E-01, 0.25897E-01,-0.22609E-02,
     a-0.33576E-02, 0.23674E-03, 0.00000E+00,
     a 0.39559E-01, 0.28143E-01,-0.34747E-02,
     a-0.35833E-02, 0.27253E-03, 0.00000E+00,
     a 0.34589E-01, 0.30612E-01,-0.52021E-02,
     a-0.36296E-02, 0.31000E-03, 0.00000E+00/

        DATA ((F(I,J),I=1,6),J=26,50) /
     a 0.31624E-01, 0.31924E-01,-0.63914E-02,
     a-0.35053E-02, 0.32820E-03, 0.00000E+00,
     a 0.28147E-01, 0.33294E-01,-0.79505E-02,
     a-0.31843E-02, 0.34451E-03, 0.00000E+00,
     a 0.26550E-01, 0.33863E-01,-0.87269E-02,
     a-0.29673E-02, 0.34999E-03, 0.00000E+00,
     a 0.24790E-01, 0.34453E-01,-0.96255E-02,
     a-0.26751E-02, 0.35449E-03, 0.00000E+00,
     a 0.22820E-01, 0.35077E-01,-0.10682E-01,
     a-0.22814E-02, 0.35761E-03, 0.00000E+00,
     a 0.20574E-01, 0.35757E-01,-0.11949E-01,
     a-0.17464E-02, 0.35876E-03, 0.00000E+00,
     a 0.17948E-01, 0.36542E-01,-0.13506E-01,
     a-0.10079E-02, 0.35702E-03, 0.00000E+00,
     a 0.10797E-01, 0.38940E-01,-0.18061E-01,
     a 0.15711E-02, 0.33794E-03, 0.00000E+00,
     a-0.20627E-02, 0.46005E-01,-0.26880E-01,
     a 0.79342E-02, 0.26785E-03, 0.00000E+00,
     a-0.35183E-01, 0.96915E-01,-0.49691E-01,
     a 0.32880E-01, 0.00000E+00, 0.00000E+00,
     a-0.60305E-01, 0.27559E+00,-0.61130E-01,
     a 0.77912E-01,-0.35190E-03, 0.00000E+00,
     a 0.54072E-01, 0.42562E+00, 0.35823E-01,
     a 0.44666E-01, 0.52985E-03, 0.00000E+00,
     a 0.77445E-01, 0.19784E+00, 0.44762E-01,
     a-0.10256E-01, 0.94775E-03, 0.00000E+00,
     a 0.63154E-01, 0.12317E+00, 0.30788E-01,
     a-0.17406E-01, 0.91420E-03, 0.00000E+00,
     a 0.52841E-01, 0.96240E-01, 0.22028E-01,
     a-0.17897E-01, 0.87527E-03, 0.00000E+00,
     a 0.45834E-01, 0.84015E-01, 0.16513E-01,
     a-0.17438E-01, 0.85294E-03, 0.00000E+00,
     a 0.40733E-01, 0.77678E-01, 0.12757E-01,
     a-0.16914E-01, 0.84260E-03, 0.00000E+00,
     a 0.36770E-01, 0.74179E-01, 0.10022E-01,
     a-0.16482E-01, 0.84020E-03, 0.00000E+00,
     a 0.33525E-01, 0.72239E-01, 0.79238E-02,
     a-0.16156E-01, 0.84324E-03, 0.00000E+00,
     a 0.62067E-02, 0.84697E-01,-0.46031E-02,
     a-0.16971E-01, 0.11097E-02, 0.00000E+00,
     a-0.15727E-01, 0.11803E+00,-0.11828E-01,
     a-0.20894E-01, 0.15418E-02, 0.00000E+00,
     a-0.36559E-01, 0.16847E+00,-0.18632E-01,
     a-0.25506E-01, 0.20407E-02, 0.00000E+00,
     a-0.47913E-01, 0.24159E+00,-0.22075E-01,
     a-0.29241E-01, 0.21029E-02, 0.00000E+00,
     a-0.43743E-01, 0.25907E+00,-0.20197E-01,
     a-0.30103E-01, 0.15532E-02, 0.61947E-04,
     a-0.28166E-01, 0.26849E+00,-0.13528E-01,
     a-0.32241E-01,-0.21437E-03, 0.20818E-03/

      DATA ((F(I,J),I=1,6),J=51,75)/
     a-0.13929E-01, 0.26828E+00,-0.71639E-02,
     a-0.34018E-01,-0.17762E-02, 0.33610E-03,
     a 0.30595E-01, 0.26054E+00, 0.14740E-01,
     a-0.38154E-01,-0.63259E-02, 0.71018E-03,
     a 0.10102E+00, 0.28222E+00, 0.18667E+00,
     a 0.64051E-01, 0.10858E-01,-0.54770E-03,
     a-0.11610E+00, 0.31003E+00,-0.41429E-01,
     a-0.11141E-01, 0.11611E-01,-0.74585E-03,
     a-0.76780E-01, 0.30548E+00,-0.31554E-01,
     a-0.21945E-01, 0.59704E-02,-0.29350E-03,
     a-0.66288E-01, 0.31086E+00,-0.27855E-01,
     a-0.23491E-01, 0.47364E-02,-0.19503E-03,
     a-0.46647E-01, 0.34458E+00,-0.19399E-01,
     a-0.21409E-01, 0.33409E-02,-0.99015E-04,
     a-0.33850E-01, 0.36445E+00,-0.13283E-01,
     a-0.17719E-01, 0.28564E-02,-0.81989E-04,
     a-0.12984E-01, 0.38360E+00,-0.30164E-02,
     a-0.11111E-01, 0.21751E-02,-0.69374E-04,
     a 0.24035E-01, 0.39194E+00, 0.15780E-01,
     a 0.94189E-03, 0.10389E-02,-0.60773E-04,
     a 0.87275E-01, 0.35021E+00, 0.49476E-01,
     a 0.21660E-01,-0.81547E-03,-0.65698E-04,
     a 0.14055E+00, 0.18448E+00, 0.82260E-01,
     a 0.39004E-01,-0.21678E-02,-0.11007E-03,
     a 0.78723E-01, 0.40297E-01, 0.56224E-01,
     a 0.18446E-01, 0.00000E+00, 0.00000E+00,
     a-0.48864E-02, 0.66214E-01, 0.13522E-01,
     a-0.90949E-02, 0.26058E-02,-0.16588E-03,
     a-0.45902E-01, 0.12959E+00,-0.89423E-02,
     a-0.22539E-01, 0.38383E-02,-0.15366E-03,
     a-0.70098E-01, 0.21934E+00,-0.23488E-01,
     a-0.30369E-01, 0.45766E-02,-0.13402E-03,
     a-0.72297E-01, 0.27160E+00,-0.25723E-01,
     a-0.30953E-01, 0.46885E-02,-0.12334E-03,
     a-0.62053E-01, 0.33097E+00,-0.21337E-01,
     a-0.27232E-01, 0.44871E-02,-0.11382E-03,
     a-0.40777E-01, 0.36825E+00,-0.10147E-01,
     a-0.19502E-01, 0.38918E-02,-0.11154E-03,
     a-0.16572E-01, 0.38230E+00, 0.39313E-02,
     a-0.10273E-01, 0.30287E-02,-0.11634E-03,
     a 0.67391E-02, 0.38034E+00, 0.19113E-01,
     a-0.58746E-03, 0.19750E-02,-0.12861E-03,
     a 0.26326E-01, 0.36741E+00, 0.34108E-01,
     a 0.88341E-02, 0.80500E-03,-0.14952E-03,
     a 0.40168E-01, 0.34751E+00, 0.48019E-01,
     a 0.17593E-01,-0.42459E-03,-0.18031E-03,
     a 0.46705E-01, 0.32332E+00, 0.60218E-01,
     a 0.25589E-01,-0.16814E-02,-0.22222E-03,
     a 0.44481E-01, 0.29578E+00, 0.70007E-01,
     a 0.32971E-01,-0.29549E-02,-0.27657E-03/

         DATA ((F(I,J),I=1,6),J=76,100)/
     a 0.33608E-01, 0.26270E+00, 0.75760E-01,
     a 0.40105E-01,-0.42390E-02,-0.34410E-03,
     a 0.28649E-01, 0.24144E+00, 0.75771E-01,
     a 0.43674E-01,-0.48635E-02,-0.38144E-03,
     a 0.35068E-01, 0.21374E+00, 0.72274E-01,
     a 0.47017E-01,-0.54139E-02,-0.41658E-03,
     a 0.81827E-01, 0.17730E+00, 0.64112E-01,
     a 0.48987E-01,-0.57314E-02,-0.43758E-03,
     a 0.21365E+00, 0.14195E+00, 0.55898E-01,
     a 0.46744E-01,-0.55866E-02,-0.42680E-03,
     a 0.37952E+00, 0.13808E+00, 0.62233E-01,
     a 0.40288E-01,-0.52980E-02,-0.40759E-03,
     a 0.44461E+00, 0.16064E+00, 0.80988E-01,
     a 0.36307E-01,-0.55954E-02,-0.43560E-03,
     a 0.43050E+00, 0.17880E+00, 0.97444E-01,
     a 0.36235E-01,-0.63249E-02,-0.50010E-03,
     a 0.39809E+00, 0.18667E+00, 0.10869E+00,
     a 0.37708E-01,-0.71495E-02,-0.57577E-03,
     a 0.36838E+00, 0.18805E+00, 0.11650E+00,
     a 0.39460E-01,-0.79622E-02,-0.65449E-03,
     a 0.34466E+00, 0.18601E+00, 0.12236E+00,
     a 0.41100E-01,-0.87514E-02,-0.73543E-03,
     a 0.31151E+00, 0.17745E+00, 0.13099E+00,
     a 0.43785E-01,-0.10298E-01,-0.90741E-03,
     a 0.28997E+00, 0.16662E+00, 0.13728E+00,
     a 0.45753E-01,-0.11871E-01,-0.10993E-02,
     a 0.27433E+00, 0.15454E+00, 0.14127E+00,
     a 0.47200E-01,-0.13551E-01,-0.13199E-02,
     a 0.26290E+00, 0.14104E+00, 0.13845E+00,
     a 0.48329E-01,-0.15130E-01,-0.15339E-02,
     a 0.25907E+00, 0.13397E+00, 0.13224E+00,
     a 0.48853E-01,-0.15677E-01,-0.16052E-02,
     a 0.25282E+00, 0.11273E+00, 0.10271E+00,
     a 0.50113E-01,-0.16651E-01,-0.17219E-02,
     a 0.25252E+00, 0.95467E-01, 0.71494E-01,
     a 0.50491E-01,-0.17088E-01,-0.17714E-02,
     a 0.28510E+00, 0.69485E-01,-0.16901E-01,
     a 0.45508E-01,-0.16106E-01,-0.16737E-02,
     a 0.36293E+00, 0.23091E+00, 0.18762E+00,
     a 0.28516E-01,-0.89499E-02,-0.93233E-03,
     a 0.29646E+00, 0.19665E+00, 0.20709E+00,
     a 0.40482E-01,-0.12977E-01,-0.13551E-02,
     a 0.26574E+00, 0.15642E+00, 0.17270E+00,
     a 0.46306E-01,-0.15567E-01,-0.16416E-02,
     a 0.26018E+00, 0.14753E+00, 0.16541E+00,
     a 0.47138E-01,-0.16277E-01,-0.17332E-02,
     a 0.24600E+00, 0.12440E+00, 0.10800E+00,
     a 0.40000E-01,-0.30000E-02, 0.00000E+00,
     a 0.23500E+00, 0.10630E+00, 0.92000E-01,
     a 0.37160E-01,-0.56700E-02, 0.00000E+00/

      DATA ((F(I,J),I=1,6),J=101,125)/
     a 0.23400E+00, 0.10295E+00, 0.77000E-01,
     a 0.36450E-01,-0.65000E-02, 0.00000E+00,
     a 0.23960E+00, 0.10060E+00, 0.50000E-01,
     a 0.35950E-01,-0.70800E-02, 0.00000E+00,
     a 0.24200E+00, 0.10000E+00, 0.47000E-01,
     a 0.35740E-01,-0.73300E-02, 0.00000E+00,
     a 0.29320E+00, 0.10400E+00, 0.46330E-01,
     a 0.35450E-01,-0.76700E-02, 0.00000E+00,
     a 0.30600E+00, 0.11700E+00, 0.47000E-01,
     a 0.35380E-01,-0.77500E-02, 0.00000E+00,
     a 0.37000E+00, 0.18200E+00, 0.80130E-01,
     a 0.34720E-01,-0.81700E-02, 0.00000E+00,
     a 0.38100E+00, 0.18467E+00, 0.93380E-01,
     a 0.34450E-01,-0.83300E-02, 0.00000E+00,
     a 0.38050E+00, 0.18600E+00, 0.10000E+00,
     a 0.34320E-01,-0.84200E-02, 0.00000E+00,
     a 0.37900E+00, 0.19000E+00, 0.10943E+00,
     a 0.33920E-01,-0.86700E-02, 0.00000E+00,
     a 0.31950E+00, 0.18831E+00, 0.12200E+00,
     a 0.33390E-01,-0.90000E-02, 0.00000E+00,
     a 0.26000E+00, 0.18663E+00, 0.12680E+00,
     a 0.32860E-01,-0.96000E-02, 0.00000E+00,
     a 0.22500E+00, 0.18241E+00, 0.13640E+00,
     a 0.31530E-01,-0.11100E-01, 0.00000E+00,
     a 0.21862E+00, 0.17987E+00, 0.14000E+00,
     a 0.30730E-01,-0.12000E-01, 0.00000E+00,
     a 0.20800E+00, 0.17566E+00, 0.14025E+00,
     a 0.29400E-01,-0.13500E-01, 0.00000E+00,
     a 0.19433E+00, 0.16300E+00, 0.14100E+00,
     a 0.24600E-01,-0.18000E-01, 0.00000E+00,
     a 0.17156E+00, 0.16150E+00, 0.14400E+00,
     a 0.16600E-01,-0.28000E-01, 0.00000E+00,
     a 0.16300E+00, 0.16300E+00, 0.16000E+00,
     a 0.14000E-01,-0.32940E-01, 0.00000E+00,
     a 0.15500E+00, 0.16700E+00, 0.17900E+00,
     a 0.13000E-01,-0.35880E-01, 0.20000E-02,
     a 0.13692E+00, 0.17100E+00, 0.19800E+00,
     a 0.14000E-01,-0.38820E-01, 0.70000E-02,
     a 0.12427E+00, 0.17931E+00, 0.21641E+00,
     a 0.16000E-01,-0.40000E-01, 0.13130E-01,
     a 0.11885E+00, 0.18287E+00, 0.22430E+00,
     a 0.20800E-01,-0.40000E-01, 0.17000E-01,
     a 0.11523E+00, 0.18525E+00, 0.22956E+00,
     a 0.24000E-01,-0.38400E-01, 0.20000E-01,
     a 0.10981E+00, 0.18881E+00, 0.23744E+00,
     a 0.46800E-01,-0.36000E-01, 0.24500E-01,
     a 0.10800E+00, 0.19000E+00, 0.24007E+00,
     a 0.54400E-01,-0.29700E-01, 0.26000E-01,
     a 0.10600E+00, 0.20500E+00, 0.25059E+00,
     a 0.84800E-01,-0.45000E-02, 0.32000E-01/

        DATA ((F(I,J),I=1,6),J=126,150)/
     a 0.11200E+00, 0.21250E+00, 0.25585E+00,
     a 0.10000E+00, 0.81000E-02, 0.32800E-01,
     a 0.14920E+00, 0.22650E+00, 0.26374E+00,
     a 0.12280E+00, 0.27000E-01, 0.34000E-01,
     a 0.17400E+00, 0.23950E+00, 0.26900E+00,
     a 0.13800E+00, 0.39600E-01, 0.33200E-01,
     a 0.21120E+00, 0.25900E+00, 0.26800E+00,
     a 0.16080E+00, 0.58500E-01, 0.32000E-01,
     a 0.23600E+00, 0.27200E+00, 0.26734E+00,
     a 0.17600E+00, 0.71100E-01, 0.29800E-01,
     a 0.26080E+00, 0.28500E+00, 0.26667E+00,
     a 0.18200E+00, 0.83700E-01, 0.27600E-01,
     a 0.27320E+00, 0.28700E+00, 0.26500E+00,
     a 0.18200E+00, 0.90000E-01, 0.26500E-01,
     a 0.31040E+00, 0.29300E+00, 0.23500E+00,
     a 0.18200E+00, 0.92000E-01, 0.23200E-01,
     a 0.33520E+00, 0.29100E+00, 0.21500E+00,
     a 0.18200E+00, 0.90000E-01, 0.21000E-01,
     a 0.36000E+00, 0.28900E+00, 0.20720E+00,
     a 0.17800E+00, 0.83500E-01, 0.18800E-01,
     a 0.36261E+00, 0.28300E+00, 0.18380E+00,
     a 0.15571E+00, 0.64000E-01, 0.12200E-01,
     a 0.36348E+00, 0.27675E+00, 0.17600E+00,
     a 0.14759E+00, 0.61360E-01, 0.10000E-01,
     a 0.36609E+00, 0.25800E+00, 0.17100E+00,
     a 0.12324E+00, 0.53450E-01, 0.34000E-02,
     a 0.36783E+00, 0.25667E+00, 0.17150E+00,
     a 0.10700E+00, 0.48180E-01,-0.10000E-02,
     a 0.36800E+00, 0.25600E+00, 0.17175E+00,
     a 0.10650E+00, 0.45550E-01,-0.20000E-02,
     a 0.38850E+00, 0.26850E+00, 0.17200E+00,
     a 0.10600E+00, 0.42910E-01,-0.30000E-02,
     a 0.40900E+00, 0.28100E+00, 0.17806E+00,
     a 0.10550E+00, 0.40270E-01,-0.40000E-02,
     a 0.42950E+00, 0.29350E+00, 0.18413E+00,
     a 0.10600E+00, 0.37640E-01,-0.43300E-02,
     a 0.45000E+00, 0.30600E+00, 0.19019E+00,
     a 0.11150E+00, 0.35000E-01,-0.46700E-02,
     a 0.45150E+00, 0.31000E+00, 0.21075E+00,
     a 0.12800E+00, 0.29600E-01,-0.43300E-02,
     a 0.45200E+00, 0.30900E+00, 0.21800E+00,
     a 0.12900E+00, 0.27800E-01,-0.40000E-02,
     a 0.45000E+00, 0.28583E+00, 0.21867E+00,
     a 0.13000E+00, 0.26000E-01,-0.25000E-02,
     a 0.39000E+00, 0.23950E+00, 0.22000E+00,
     a 0.12800E+00, 0.22400E-01, 0.50000E-03,
     a 0.30000E+00, 0.17000E+00, 0.19500E+00,
     a 0.11675E+00, 0.17000E-01, 0.50000E-02,
     a 0.27000E+00, 0.16000E+00, 0.18980E+00,
     a 0.11300E+00, 0.15600E-01, 0.60000E-02/

      DATA ((F(I,J),I=1,6),J=151,175)/
     a 0.26257E+00, 0.15000E+00, 0.18460E+00,
     a 0.10925E+00, 0.14200E-01, 0.70000E-02,
     a 0.24029E+00, 0.12000E+00, 0.16900E+00,
     a 0.98000E-01, 0.10000E-01, 0.50000E-02,
     a 0.21800E+00, 0.10500E+00, 0.15340E+00,
     a 0.89000E-01, 0.58000E-02, 0.32000E-02,
     a 0.21380E+00, 0.95000E-01, 0.14300E+00,
     a 0.83000E-01, 0.30000E-02, 0.20000E-02,
     a 0.20751E+00, 0.80000E-01, 0.13580E+00,
     a 0.74000E-01,-0.30000E-03, 0.11400E-02,
     a 0.19987E+00, 0.64830E-01, 0.11900E+00,
     a 0.57000E-01,-0.80000E-02, 0.00000E+00,
     a 0.19644E+00, 0.54000E-01, 0.11000E+00,
     a 0.46000E-01,-0.12000E-01, 0.00000E+00,
     a 0.19586E+00, 0.45000E-01, 0.87000E-01,
     a 0.31000E-01,-0.14000E-01, 0.00000E+00,
     a 0.19700E+00, 0.44200E-01, 0.81400E-01,
     a 0.28600E-01,-0.14200E-01, 0.00000E+00,
     a 0.20675E+00, 0.43000E-01, 0.73000E-01,
     a 0.25000E-01,-0.13000E-01, 0.00000E+00,
     a 0.22300E+00, 0.45000E-01, 0.57500E-01,
     a 0.21000E-01,-0.10000E-01, 0.00000E+00,
     a 0.26400E+00, 0.56800E-01, 0.40330E-01,
     a 0.19000E-01,-0.52000E-02, 0.00000E+00,
     a 0.29000E+00, 0.64400E-01, 0.37000E-01,
     a 0.19000E-01,-0.36000E-02, 0.00000E+00,
     a 0.31600E+00, 0.72000E-01, 0.37000E-01,
     a 0.19000E-01,-0.20000E-02, 0.00000E+00,
     a 0.32900E+00, 0.85600E-01, 0.37000E-01,
     a 0.19400E-01,-0.10500E-02, 0.00000E+00,
     a 0.34200E+00, 0.99200E-01, 0.39500E-01,
     a 0.19800E-01,-0.10000E-03, 0.00000E+00,
     a 0.40985E+00, 0.14000E+00, 0.47000E-01,
     a 0.21000E-01, 0.27500E-02, 0.00000E+00,
     a 0.52292E+00, 0.22667E+00, 0.70000E-01,
     a 0.23000E-01, 0.75000E-02, 0.00000E+00,
     a 0.63600E+00, 0.31333E+00, 0.10150E+00,
     a 0.26000E-01, 0.12250E-01, 0.00000E+00,
     a 0.64500E+00, 0.35667E+00, 0.11725E+00,
     a 0.27500E-01, 0.14630E-01, 0.00000E+00,
     a 0.63700E+00, 0.39500E+00, 0.13300E+00,
     a 0.29000E-01, 0.17000E-01, 0.00000E+00,
     a 0.51433E+00, 0.43000E+00, 0.15380E+00,
     a 0.33000E-01, 0.20200E-01, 0.00000E+00,
     a 0.48367E+00, 0.43000E+00, 0.15900E+00,
     a 0.34000E-01, 0.21000E-01, 0.00000E+00,
     a 0.36100E+00, 0.44000E+00, 0.17420E+00,
     a 0.38000E-01, 0.21800E-01, 0.00000E+00,
     a 0.35367E+00, 0.44200E+00, 0.17800E+00,
     a 0.39000E-01, 0.22000E-01, 0.00000E+00/

        DATA ((F(I,J),I=1,6),J=176,200)/
     a 0.31700E+00, 0.43800E+00, 0.19150E+00,
     a 0.48000E-01, 0.22000E-01, 0.57000E-03,
     a 0.29800E+00, 0.41800E+00, 0.20500E+00,
     a 0.69180E-01, 0.22000E-01, 0.20000E-02,
     a 0.28660E+00, 0.40500E+00, 0.21014E+00,
     a 0.81880E-01, 0.22600E-01, 0.53000E-02,
     a 0.27900E+00, 0.40586E+00, 0.21357E+00,
     a 0.90350E-01, 0.23000E-01, 0.75000E-02,
     a 0.27220E+00, 0.40737E+00, 0.21957E+00,
     a 0.10518E+00, 0.30000E-01, 0.11350E-01,
     a 0.26800E+00, 0.41956E+00, 0.22557E+00,
     a 0.12000E+00, 0.45830E-01, 0.13670E-01,
     a 0.27010E+00, 0.43000E+00, 0.23071E+00,
     a 0.12480E+00, 0.59390E-01, 0.13500E-01,
     a 0.27500E+00, 0.47520E+00, 0.23757E+00,
     a 0.12925E+00, 0.77480E-01, 0.10690E-01,
     a 0.28125E+00, 0.48650E+00, 0.23929E+00,
     a 0.12988E+00, 0.82000E-01, 0.95400E-02,
     a 0.29375E+00, 0.50910E+00, 0.24271E+00,
     a 0.13112E+00, 0.86000E-01, 0.72300E-02,
     a 0.30000E+00, 0.52040E+00, 0.24443E+00,
     a 0.13175E+00, 0.86000E-01, 0.60800E-02,
     a 0.32091E+00, 0.53170E+00, 0.24614E+00,
     a 0.13238E+00, 0.86000E-01, 0.49200E-02,
     a 0.34182E+00, 0.54300E+00, 0.24786E+00,
     a 0.13300E+00, 0.84000E-01, 0.37700E-02,
     a 0.36273E+00, 0.54767E+00, 0.24957E+00,
     a 0.13307E+00, 0.82000E-01, 0.26200E-02,
     a 0.40455E+00, 0.55700E+00, 0.25300E+00,
     a 0.13320E+00, 0.74330E-01, 0.31000E-03,
     a 0.44636E+00, 0.55433E+00, 0.25643E+00,
     a 0.13333E+00, 0.66670E-01,-0.20000E-02,
     a 0.50000E+00, 0.50000E+00, 0.24000E+00,
     a 0.10000E+00, 0.24000E-01,-0.12000E-01,
     a 0.63370E+00, 0.52990E+00, 0.35260E+00,
     a 0.21080E+00, 0.80000E-01, 0.80000E-02,
     a 0.65480E+00, 0.53600E+00, 0.36120E+00,
     a 0.21800E+00, 0.80000E-01, 0.80000E-02,
     a 0.65100E+00, 0.50000E+00, 0.35230E+00,
     a 0.20670E+00, 0.66250E-01, 0.73300E-02,
     a 0.56740E+00, 0.44000E+00, 0.33210E+00,
     a 0.16890E+00, 0.40000E-01, 0.51100E-02,
     a 0.50330E+00, 0.40700E+00, 0.31870E+00,
     a 0.14100E+00, 0.30000E-01, 0.40000E-02,
     a 0.46000E+00, 0.36000E+00, 0.30000E+00,
     a 0.11000E+00, 0.40000E-02, 0.10000E-02,
     a 0.37450E+00, 0.33380E+00, 0.30190E+00,
     a 0.10150E+00, 0.13710E-01, 0.16480E-01,
     a 0.28000E+00, 0.23000E+00, 0.29000E+00,
     a 0.90000E-01, 0.10000E-02, 0.50000E-02/

      DATA ((F(I,J),I=1,6),J=201,220)/
     a 0.46000E+00, 0.38000E+00, 0.31000E+00,
     a 0.14000E+00, 0.40000E-01, 0.50000E-02,
     a 0.55000E+00, 0.42000E+00, 0.33000E+00,
     a 0.15000E+00, 0.30000E-01, 0.60000E-02,
     a 0.53000E+00, 0.43000E+00, 0.33000E+00,
     a 0.13000E+00, 0.40000E-01, 0.35000E-01,
     a 0.60000E+00, 0.46000E+00, 0.32000E+00,
     a 0.13000E+00, 0.40000E-01, 0.10000E-01,
     a 0.56000E+00, 0.44000E+00, 0.31000E+00,
     a 0.14000E+00, 0.30000E-01, 0.20000E-01,
     a 0.54000E+00, 0.41000E+00, 0.32000E+00,
     a 0.14000E+00, 0.50000E-01, 0.30000E-01,
     a 0.51000E+00, 0.40000E+00, 0.32000E+00,
     a 0.16000E+00, 0.50000E-01, 0.20000E-01,
     a 0.58000E+00, 0.46000E+00, 0.35000E+00,
     a 0.18000E+00, 0.70000E-01, 0.30000E-01,
     a 0.57000E+00, 0.44000E+00, 0.34000E+00,
     a 0.18000E+00, 0.60000E-01, 0.20000E-01,
     a 0.60000E+00, 0.45000E+00, 0.32000E+00,
     a 0.20000E+00, 0.70000E-01, 0.40000E-01,
     a 0.59000E+00, 0.47000E+00, 0.32000E+00,
     a 0.20000E+00, 0.80000E-01, 0.40000E-01,
     a 0.63000E+00, 0.46000E+00, 0.31000E+00,
     a 0.20000E+00, 0.90000E-01, 0.47000E-01,
     a 0.66000E+00, 0.48000E+00, 0.32000E+00,
     a 0.21000E+00, 0.10000E+00, 0.20000E-01,
     a 0.68000E+00, 0.49000E+00, 0.33000E+00,
     a 0.20000E+00, 0.10700E+00, 0.50000E-01,
     a 0.70000E+00, 0.49000E+00, 0.33000E+00,
     a 0.20000E+00, 0.11100E+00, 0.50000E-01,
     a 0.72000E+00, 0.51000E+00, 0.35000E+00,
     a 0.20000E+00, 0.11200E+00, 0.43000E-01,
     a 0.73000E+00, 0.52000E+00, 0.35000E+00,
     a 0.20000E+00, 0.11300E+00, 0.40000E-01,
     a 0.77009,     0.53175,     0.36398,
     a 0.22328,     0.098031,    0.022452,
     b 0.78288,     0.55059,     0.36917,
     b 0.22776,     0.10595,     0.027254,
     c 0.78302,     0.51968,     0.29398,
     c 0.1419,      0.064051,    0.020662/
       END
C =====================================================================
C =====================================================================
C     LIST OF ROUTIENS IN THIS FILE
C ---------------------------------------------------------------------
C     F: choosa
C     F: choosl
C     F: choosm
C     F: choosn
C     F: choosp
C     S: p3alph
C =====================================================================
C =====================================================================
      Function CHOOSA(Emax)
C     Purpose is to choose an alpha energy from a simple formalism
C     modestly representing a distribution of the energy of one
C     of the alphas (the first in our case) from the 12-C --> 3 alphas
C     3-body breakup reaction.
C
C     The formalism used is:   Phi = E**2 * (Emax - E)**2
      real*8 UNIRN

      Dimension E(101), S(101)
      Data N/100/, M/101/, Nterp/1/, Eold/0.0/

      Em=Emax
C     IF (Em .EQ. Eold)  goto 8
      Eold=Em
      D=0.01*Em
      Halfd=0.5*D
      E(1)=0.0
      S(1)=0.0
      Do 5 I=1,N
      Ek=E(I)+Halfd
      E(I+1)=Ek+Halfd
      Sqrtphi=(Em-Ek)*Ek
    5 S(I+1)=S(I) + Sqrtphi*Sqrtphi

      Sx=S(M)

    8 R=Sx*sngl(UNIRN(dummy))
      Chosen=EXTERP(S,E,R,M,Nterp)
      CHOOSA=Chosen
      Return
      END
C =====================================================================
C =====================================================================
      Function CHOOSL(Fi,N)
C     Purpose to choose an azimuthal angle, or actually the cosine of the
C     chosen angle, from a distribution described by Legendre polynomial
C     expansion.  The technique is from an old "theorem" in which one
C     computes the maximum value the function can have within the
C     independent variable's domain of interest (in this case from -1
C     to +1) and then compares that to a value of the function for a
C     randomly chosen value of the independent variable.

      real*8 UNIRN

      Dimension P(10),Fi(1)
C     There are N (where N .LE. 10, so be careful) coefficients Fi at
C     entry.

C     Next portion of programming is from O5S code.
      A = 1.
      DO 19 J = 1, N
   19 A = A + FLOAT(J+J+1) * ABS (FI(J))
C     In this case -A- is .GE. the maximum value that can be
C     obtained from this polynomial expansion, and not just
C     equal to the maximum value.


   20 P(1) = 2. * sngl(UNIRN(dummy)) - 1.
      P(2) = ( 3. * P(1)**2 - 1. ) / 2.
      NM = N - 1
      FL = 1.
      DO 21 L = 2, NM
      FL = FL + 1.
      P(L+1) = ( P(1) * ( FL + FL + 1. ) * P(L) - FL * P(L-1) ) /
     1  ( FL + 1. )
   21 CONTINUE
      X = 1.
      FL = 0.
      DO 22 L = 1, N
      FL = FL + 1.
      X = X + ( FL + FL +1. ) * FI(L) * P(L)
   22 CONTINUE
      R=sngl(UNIRN(dummy))
      IF ( X .LT. A * R )  GOTO 20
      CHOOSL = P(1)
      Return
      END
C =====================================================================
C =====================================================================
      Function CHOOSM(Enmin,Enmax,Temp)
C     Purpose is to choose a neutron energy by random number generation from
C     a Maxwellian distribution of energies given by:
C
C     Phi(En) = SQRT(En) * EXP(-En/T),
C
C     where T is a "temperature".

      real*8 UNIRN
      Dimension En(401),Sn(401)
      Data Eold/0.0/, Elow/0.0/, Told/0.0/

      N=100
      Elo=Enmin
      T=Temp
      E=Enmax
      IF (E .EQ. Eold  .AND.  T .EQ. Told) goto 8

C     Set up new En, Sn arrays only for new value of Enmax and/or Temp
      Told=T
      Eold=E
      Edel=E
      IF (Edel .LE. 10.0) goto 2
      IF (Edel .GT. 20.0*T) Edel=20.0*T
      N=10*IFIX(Edel)
      IF (N .GT. 400) N=400
    2 D=Edel/FLOAT(N)
      En(1)=0.0
      Sn(1)=0.0
      Halfd=0.5*D
      Do 5 I=1,N
      Ek=En(I)+Halfd
      En(I+1)=Ek+Halfd
      Phi=SQRT(Ek)*EXP(-Ek/T)
    5 Sn(I+1)=Sn(I) + Phi
      M=N+1
      Snx=Sn(M)

    8 IF (Elo .EQ. Elow) goto 10
      Elow=Elo
      Sny=EXTERP(En,Sn,Elo,M,1)
   10 R=Sny + (Snx-Sny)*sngl(UNIRN(dummy))
      Chosen=EXTERP(Sn,En,R,M,1)
   15 CHOOSM=Chosen
      Return
      End
C =====================================================================
C =====================================================================
      Function CHOOSN(Enmin,Enmax,Temp)
C     Purpose is to choose a neutron energy by random number generation from
C     a distribution of energies given by:
C
C     Phi(En) = EXP(-En/T) * En**1.5,
C
C     where T is a "temperature".  This functional dependence approximates
C     the relative "continuum" (of) 12-C excitations by neutrons having
C     20 MeV < En < 70 MeV as computed by the Hauser-Feshbach compound-
C     nucleus plus pre-compound formalism in the code TNG.

      real*8 UNIRN
      Dimension En(401),Sn(401)
      Data Eold/0.0/, Elow/0.0/, Told/0.0/

      N=100
      Elo=Enmin
      T=Temp
      E=Enmax

      SNX=0.0
      SNY=0.0

      IF (E .EQ. Eold  .AND.  T .EQ. Told) then
        continue
      end if

C     Set up new En, Sn arrays only for new value of Enmax and/or Temp
      Told=T
      Eold=E
      Edel=E
      IF (Edel .LE. 10.0) goto 2
      IF (Edel .GT. 20.0*T) Edel=20.0*T
      N=10*IFIX(Edel)
      IF (N .GT. 400) N=400
    2 D=Edel/FLOAT(N)
      En(1)=0.0
      Sn(1)=0.0
      Halfd=0.5*D
      Do 5 I=1,N
      Ek=En(I)+Halfd
      En(I+1)=Ek+Halfd
      Phi=Ek**1.5 * EXP(-Ek/T)
    5 Sn(I+1)=Sn(I) + Phi
      M=N+1
      Snx=Sn(M)

    8 IF (Elo .EQ. Elow) goto 10

      Elow=Elo
      Sny=EXTERP(En,Sn,Elo,M,1)
   10 R=Sny + (Snx-Sny)*sngl(UNIRN(dummy))
      Chosen=EXTERP(Sn,En,R,M,1)
   15 CHOOSN=Chosen
      Return
      End
C =====================================================================
C =====================================================================
      Function CHOOSP(Epmax,Temp)
C     Purpose is to choose a proton energy from a simple formalism for the
C     (n,p) and (n,np) reactions.
C
C     The formulation used is:
C
C     Phi(Eproton) = Eprot*EXP(-Eprot/Temp)*(1.0-EXP(-K*Eprot/B))
C
C     where Temp = "Nuclear Temperature", B=barrier penetration
C     factor , and K=an empirical parameter selected to give the
C     best fit of (1.0-EXP(-K*Eprot/B)) to the barrier penetration
C     probability of 1.5 MeV.  The value for B was deduced to be
C     2.9 MeV for the O5S program, and it appears that the value
C     of K was taken to be 1.5 MeV.  That's why the "B" in the
C     Data statement below was set = (2.9/1.5).  The given
C     functional dependence reproduces adequately the shape of the
C     n + 12-C --> p + 12-B proton continuum (in 12-B) excitation
C     as computed by TNG.

      real*8 UNIRN
      Dimension Epr(401),Sp(401)
      Data B/1.93333/, Eold/0.0/, Told/0.0/, Nterp/1/

      N=100
      E=Epmax
      IF (E .EQ. Eold .AND. Temp .EQ. Told) then
        CONTINUE
      end if

C     Set up Epr and Sp arrays for new energy Epmax or Temp
      Told= Temp
      Eold=E
      IF (E .LE. 10.0) goto 2
      N=10*IFIX(E)
      IF (N .GT. 400) N=400
    2 D=E/FLOAT(N)
      Epr(1)=0.0
      Sp(1)=0.0
      Halfd=0.5*D

      Do 5 J=1,N
      Eprj=Epr(J)+Halfd
      Epr(J+1)=Eprj+Halfd
      Phi=Eprj*EXP(-Eprj/Temp) * (1.0 - EXP(-Eprj/B))
      Sp(J+1)=Sp(J)+Phi
    5 Continue
      M=N+1
      Spx=Sp(M)

   10 Rspx=Spx*sngl(UNIRN(dummy))
      Chosen=EXTERP(Sp,Epr,Rspx,M,Nterp)
   15 CHOOSP=Chosen
      Return
      END
C =====================================================================
C =====================================================================
      SUBROUTINE P3ALPH(E)
C     Purpose to return array PBAR, the summed probabilities for the 14 (in
C     effect) sub-reactions leading to the (n,n' 3alpha) break-up reaction.
C     See Subroutine NN3ALF for the specific branching designations.

      Common /P3ALF/ Pbar(13)
!$OMP THREADPRIVATE(/P3ALF/)

      Dimension Eofp(26),P1(26),P2(26),P3(26),P4(26),P5(26),P6(26),
     X P7(26),P8(26),P9(26),P10(26),P11(26),P12(26),P13(26)

      Data Eofp/9.7, 10.2,  10.5,  11.0,  11.85, 12.0,  12.5, 13.0,
     V 13.2,  13.7,  14.3,  15.3,  16.2,  17.5,  18.5,  19.6,  20.,
     W 22.,   23.2,  25.,   30.,   35.,   40.,   50.,   60.,   70./
      Data P1/1.0,   0.86,  0.672, 0.244, 0.105, 0.095, 0.074, 0.059,
     V 0.054, 0.045, 0.04,  0.032, 0.029, 0.031, 0.31,  0.034, 0.034,
     W 0.022, 0.017, 0.02,  0.026, 0.023, 0.016, 0.01,  0.0,   0.0/
      Data P2/1.0,   0.86,  0.672, 0.478, 0.476, 0.473, 0.493, 0.469,
     V 0.455, 0.43,  0.413, 0.369, 0.337, 0.322, 0.271, 0.238, 0.222,
     W 0.174, 0.180, 0.181, 0.183, 0.201, 0.196, 0.195, 0.195, 0.195/
      Data P3/1.0,   1.0,   0.943, 0.822, 0.748, 0.738, 0.729, 0.666,
     V 0.636, 0.577, 0.531, 0.449, 0.397, 0.369, 0.309, 0.272, 0.253,
     W 0.196, 0.196, 0.194, 0.191, 0.206, 0.196, 0.195, 0.195, 0.195/
      Data P4/1.0,   1.0,   1.0,   1.0,   1.0,   0.99,  0.955, 0.863,
     V 0.819, 0.73,  0.659, 0.536, 0.463, 0.422, 0.351, 0.309, 0.288,
     W 0.220, 0.214, 0.209, 0.199, 0.21,  0.196, 0.195, 0.195, 0.195/
      Data P5/1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   0.991, 0.921,
     V 0.88,  0.79,  0.729, 0.601, 0.52,  0.472, 0.391, 0.346, 0.322,
     W 0.242, 0.229, 0.218, 0.203, 0.212, 0.196, 0.195, 0.195, 0.195/
      Data P6/1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     W 0.994, 0.875, 0.938, 0.802, 0.696, 0.601, 0.508, 0.447, 0.416,
     W 0.305, 0.278, 0.254, 0.222, 0.222, 0.201, 0.195, 0.195, 0.195/
      Data P7/1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     W 1.0,   1.0,   0.997, 0.907, 0.827, 0.766, 0.662, 0.598, 0.563,
     W 0.407, 0.349, 0.301, 0.276, 0.256, 0.217, 0.197, 0.195, 0.195/
      Data P8/1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     W 1.0,   1.0,   1.0,   1.0,   0.986, 0.96,  0.833, 0.755, 0.709,
     W 0.505, 0.425, 0.358, 0.301, 0.269, 0.224, 0.197, 0.195, 0.195/
      Data P9/1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     W 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   0.94,  0.792, 0.744,
     W 0.535, 0.452, 0.384, 0.317, 0.278, 0.229, 0.197, 0.195, 0.195/
      Data P10/1.0,  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     W 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   0.985, 0.93,
     W 0.653, 0.533, 0.441, 0.373, 0.307, 0.247, 0.203, 0.195, 0.195/
      Data P11/1.0,  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     W 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   0.985, 0.977,
     W 0.92,  0.894, 0.86,  0.81,  0.764, 0.745, 0.767, 0.765, 0.75/
      Data P12/1.0,  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     W 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     X 1.0,   0.98,  0.94,  0.845, 0.780, 0.753, 0.767, 0.765, 0.75/
      Data P13/1.0,  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     W 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     X 1.0,   1.0,   1.0,   0.96,  0.865, 0.809, 0.79,  0.765, 0.75/

      Data N/26/, Nterp/1/

      En=E
      Pbar(1) = EXTERP(Eofp, P1,En,N,Nterp)
      Pbar(2) = EXTERP(Eofp, P2,En,N,Nterp)
      Pbar(3) = EXTERP(Eofp, P3,En,N,Nterp)
      Pbar(4) = EXTERP(Eofp, P4,En,N,Nterp)
      Pbar(5) = EXTERP(Eofp, P5,En,N,Nterp)
      Pbar(6) = EXTERP(Eofp, P6,En,N,Nterp)
      Pbar(7) = EXTERP(Eofp, P7,En,N,Nterp)
      Pbar(8) = EXTERP(Eofp, P8,En,N,Nterp)
      Pbar(9) = EXTERP(Eofp, P9,En,N,Nterp)
      Pbar(10)= EXTERP(Eofp,P10,En,N,Nterp)
      Pbar(11)= EXTERP(Eofp,P11,En,N,Nterp)
      Pbar(12)= EXTERP(Eofp,P12,En,N,Nterp)
      Pbar(13)= EXTERP(Eofp,P13,En,N,Nterp)

      Return
      END
C =====================================================================
C =====================================================================
C     LIST OF ROUTIENS IN THIS FILE
C ---------------------------------------------------------------------
C     S: alfapd
C     S: alfapt
C     S: eli8dk
C     S: dpli7
C =====================================================================
      Subroutine ALFAPD(Excita,V6Li)
C     6-Li --> d + alpha.

C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
!$OMP THREADPRIVATE(/VECTOR/)
C     (Xn,Yn,Zn) are dir cos components of 6Li in detector lab system
C ---------------------------------------------------------------------

      Data Kd/10/, Ka/3/, He4pd/1.474/

      Ta=Excita-He4pd
      Ed=RCKE(Kd,Ka,Ta)
      Vd=VELOCITY(Kd,Ed)
      CALL RVECT(Xd,Yd,Zd)
      Vdx=Xd*Vd
      Vdy=Yd*Vd
      Vdz=Zd*Vd
      CALL LABTRAN(Vdx,Vdy,Vdz, V6Li, Vx,Vy,Vz)
      Vdl=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C ---------------------------------------------------------------------
C     DEUTERON
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      Cx = Xn
      Cy = Yn
      Cz = Zn
      CALL TRANSVEC(Cx,Cy,Cz)
      IKF   = 1000002
      DDENG = EFROMV(Kd,Vdl)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
C     That's it for the deuteron.  Now get alpha energy.
      Eac=Ta-Ed
      Vac=VELOCITY(Ka,Eac)
      Vax=-Xd*Vac
      Vay=-Yd*Vac
      Vaz=-Zd*Vac
      CALL LABTRAN(Vax,Vay,Vaz, V6Li, Vx,Vy,Vz)
      Va=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C ---------------------------------------------------------------------
C     ALPHA
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      CALL TRANSVEC(Cx,Cy,Cz)
      IKF   = 2000004
      DDENG = EFROMV(Ka,Va)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      Return
      END
C =====================================================================
C =====================================================================
      Subroutine ALFAPT(Excitn,V7Li)
C     Does energetics for  7-Li --> alpha + triton.

C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
!$OMP THREADPRIVATE(/VECTOR/)
C     (Xn,Yn,Zn) are dir cos components of 7Li in detector lab system
C ---------------------------------------------------------------------
c
      Data Ka/3/, Kt/13/, He4pt/2.467/
c
      Ta=Excitn-He4pt
      Et=RCKE(Kt,Ka,Ta)
      CALL RVECT(Xt,Yt,Zt)
      Vtc=VELOCITY(Kt,Et)
      Vtx=Xt*Vtc
      Vty=Yt*Vtc
      Vtz=Zt*Vtc
      CALL LABTRAN(Vtx,Vty,Vtz, V7Li, Vx,Vy,Vz)
      Vt=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C ---------------------------------------------------------------------
C     TRITON
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      Cx = Xn
      Cy = Yn
      Cz = Zn
      CALL TRANSVEC(Cx,Cy,Cz)
      IKF   = 1000003
      DDENG = EFROMV(Kt,Vt)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
      Xn = Cx
      Yn = Cy
      Zn = Cz
C ---------------------------------------------------------------------
C     Triton done.  Get alpha energy.
      Ealf=Ta-Et
      Valf=VELOCITY(Ka,Ealf)
      Vax=-Xt*Valf
      Vay=-Yt*Valf
      Vaz=-Zt*Valf
      CALL LABTRAN(Vax,Vay,Vaz, V7Li, Vx,Vy,Vz)
      Va=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C ---------------------------------------------------------------------
C     ALPHA
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      Cx = Xn
      Cy = Yn
      Cz = Zn
      CALL TRANSVEC(Cx,Cy,Cz)
      IKF   = 1000003
      DDENG = EFROMV(Ka,Va)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
      Xn = Cx
      Yn = Cy
      Zn = Cz
C ---------------------------------------------------------------------
C     and that gets the alpha energy.

      Return
      END
C =====================================================================
C =====================================================================
      Subroutine ELI8DK(Ex,V8Li,Eg)
C     Particle decay of highly excited 8-Li ion.

C     Decay modes considered in this subroutine are:
C     (a)   8-Li --> n + 7-Li; then,
C     (b)   7-Li --> gamma + 7-Li (ground state);
C     or (b')  7-Li --> triton + alpha;
C     or (b")  7-Li --> n + 6-Li; then,
C        (c")  6-Li --> gamma + 6-Li (ground state);
C     or (c'") 6-Li --> d + alpha.

      Common /MASSES/ Emass(21)
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
!$OMP THREADPRIVATE(/VECTOR/)
C     At entry, Xn,Yn,Zn are direction cosines of the moving
C     8-Li ion in "detector" laboratory coordinates.

      Common /NEUTRN/ Eneut, Dc1,Dc2,Dc3
!$OMP THREADPRIVATE(/NEUTRN/)
      Common /EXC7LI/ E7Li(6),Li6pn
      Common /EXC6LI/ ExLi6(6),He4pd
C     Array ExLi6 and variable He4pd given in Subroutine TENBDK

C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------
      Common /NEUTR2/ Eneut2, U2,V2,W2
!$OMP THREADPRIVATE(/NEUTR2/)
      real*8 UNIRN

      Data Kn/1/, K7Li/17/, K6Li/18/
      Data Li7pn/2.033/, Li6pn/7.251/
C     Next array represents level structure of 7-Li:
      Data E7Li/0.0, 0.478, 4.63, 6.68, 7.456, 9.7/

C     For 7-Li ground state is stable, 0.478-MeV level decays by
C     gamma emission, 4.63- and 6.68-MeV levels decay by triton +
C     alpha, and higher-lying levels decay (primarily) by neutron
C     emission to levels in 6-Li.

      Eg=0.0
      Ta=Ex
C     First decay (a) 8-Li --> n + 7-Li:
      Taex=Ta-Li7pn
      IF (Taex .LE. 0.0) Return
      Enmax=RCKE(Kn,K7Li,Taex)
      Ehat=Enmax+5.0
C     (No magic to 5.0; it's just a value chosen to somewhat lessen
C     the "chosen" 8-Li excitation energies, on the average.)

      Fn=0.065 + 0.001*Ehat
      Temp=Fn*Ehat
      Eloww=0.0

C     Now check for existing 2d neutron waiting to be processed.  If
C     there is one, restrict present decay to neutron-bound levels
C     of 7-Li.

      IF (Eneut2 .LE. 0.0) goto 5

      Tryy=Enmax-Li6pn
      IF (Tryy .GT. 0.0) Eloww=Tryy
    5 Enn=CHOOSN(Eloww,Enmax,Temp)
      Try=Enn*(Emass(Kn) + Emass(K7Li))/Emass(K7Li)
      Entry=RCKE(Kn,K7Li,Try)
      Try=Try*Enn/Entry
      Excitn=Taex-Try
      Level=7
C     Possibly pair -EXCITN- with level in 7-Li.
      IF (Excitn .GT. 11.0) goto 10
      K=1
      Do 8 I=2,6
      IF (Excitn .LT. E7Li(I)) goto 9
    8 K=K+1
    9 Excitn=E7Li(K)
      Level=K
   10 Taex=Taex-Excitn
      Enn=RCKE(Kn,K7Li,Taex)
      CALL RVECT(Xp,Yp,Zp)
      Vnn=VELOCITY(Kn,Enn)
      Vnx=Xp*Vnn
      Vny=Yp*Vnn
      Vnz=Zp*Vnn
      CALL LABTRAN(Vnx,Vny,Vnz, V8Li, Vx,Vy,Vz)
      Vneut=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Eneut=EFROMV(Kn,Vneut)
C     Rotate neutron coordinates to "detector" system.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
C     Done. Store information in NEUTRN common area.
      Dc1=Xn
      Dc2=Yn
      Dc3=Zn
C ---------------------------------------------------------------------
C     NEUTRON
C ---------------------------------------------------------------------
      IKF   = 2112
      DDENG = Eneut
      DDVX  = Dc1
      DDVY  = Dc2
      DDVZ  = Dc3
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

C     Now get 7-Li ion information.
      E7Lic=Taex-Enn
      V7Lic=VELOCITY(K7Li,E7Lic)
      V7x=-Xp*V7Lic
      V7y=-Yp*V7Lic
      V7z=-Zp*V7Lic
      CALL LABTRAN(V7x,V7y,V7z, V8Li, Vx,Vy,Vz)
      V7Li=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C     Test for 7-Li decay mode (b') or (b"):
      IF (Level .GE. 3) goto 15

C     Program counter to here means ground state or 0.478-MeV level
C     so tidy up 7-Li.

      Eg=Excitn
C ---------------------------------------------------------------------
C     7Li
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
      IKF   = 3*1000000 + 7
      DDENG = EFROMV(K7Li,V7Li)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      Return
C ---------------------------------------------------------------------
C     Next, decide between decay modes (b') and (b").
   15 IF (Level .GE. 5) goto 20

C     Next section for  7-Li --> alpha + triton.
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
C     (Xn,Yn,Zn) are dir. cos. of 7Li in detector system.
C ---------------------------------------------------------------------
      CALL ALFAPT(Excitn,V7Li)
      Return
C     Ends section on  7-Li --> alpha + triton.

C     Next is for second neutron emission to levels in 6-Li.
C     First transform motion of 7-Li ion to detector coordinates.
   20 CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
C     Done -- now have 7-Li ion motion in "detector" coordinates.
      Ta=Excitn-Li6pn
      Enmax=RCKE(Kn,K6Li,Ta)
      Ehat=Enmax+5.0
      Fn=0.065 + 0.001*Ehat
      Temp=Fn*Ehat
      Eloww=0.0
      Enn=CHOOSN(Eloww,Enmax,Temp)
      Try=Enn*(Emass(Kn) + Emass(K6Li))/Emass(K6Li)
      Entry=RCKE(Kn,K6Li,Try)
      Try=Try*Enn/Entry
      Excitn2=Ta-Try
C     Maybe pair -EXCITN2- with energy of a level in 6-Li.
      Level=7
      IF (Excitn2 .GT. 20.0) goto 25
      K=1
      Do 22 I=2,6
      IF (Excitn2 .LT. ExLi6(I)) goto 23
   22 K=K+1
   23 Excitn2=ExLi6(K)
      Level=K
      IF (K. EQ. 3) Eg=Excitn2

   25 Taex=Ta-Excitn2
      En2c=RCKE(Kn,K6Li,Taex)
      Vn2c=VELOCITY(Kn,En2c)
      CALL RVECT(Xp,Yp,Zp)
      Vn2x=Xp*Vn2c
      Vn2y=Yp*Vn2c
      Vn2z=Zp*Vn2c
      CALL LABTRAN(Vn2x,Vn2y,Vn2z, V7Li, Vx,Vy,Vz)
      Vn2=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      En2=EFROMV(Kn,Vn2)
      Eneut2=-En2
C     Negative value to alert Subroutine BANKER; will be corrected
C     at exit of calling routine.

C     Now rotate neutron coordinates into "detector" system.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
C     Done.

C ---------------------------------------------------------------------
C     2ND NEUTRON
C ---------------------------------------------------------------------
      IKF   = 2112
      DDENG = En2
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

C     Okay.  Now get information on the 6-Li ion.
      E6Li=Taex-En2c
      V6Lic=VELOCITY(K6Li,E6Li)
      V6x=-Xp*V6Lic
      V6y=-Yp*V6Lic
      V6z=-Zp*V6Lic
      CALL LABTRAN(V6x,V6y,V6z, V7Li, Vx,Vy,Vz)
      V6Li=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C     See subroutine TENBDK for decay characteristics of levels of 6-Li.
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
C     (Zx,Zy,Zz) are dir. cos. of 7Li in detecotr lab system.
      CALL TRANSVEC(Zx,Zy,Zz)
C     (Xn,Yn,Zn) are dir. cos. of 6Li in detecotr lab system.
C ---------------------------------------------------------------------

      IF (Level .EQ. 2) goto 30
      IF (Level .GE. 4) goto 30

C     Program counter to here for Level=1 or 3.  Tidy up 6-Li.
C ---------------------------------------------------------------------
C     6Li
C ---------------------------------------------------------------------
      IKF   = 3*1000000 + 6
      DDENG = EFROMV(K6Li,V6Li)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      Return

C     Last section: 6-Li --> d + alpha.
   30 CALL ALFAPD(Excitn2,V6Li)
      Return

      END
C =====================================================================
C =====================================================================
      Subroutine DPLI7(Excit,V9Be)
C     Does energetics for  9-Be --> d + 7-Li
C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
!$OMP THREADPRIVATE(/VECTOR/)
C     (Xn,Yn,Zn) are dir cos components of 9Be in detector lab system
C ---------------------------------------------------------------------
      Data  Kd/10/, K7Li/17/, Dp7Li/16.6965/

      Ta=Excit- Dp7Li
      Ed=RCKE(Kd,K7Li,Ta)
      CALL RVECT(Xd,Yd,Zd)
      Vdc=VELOCITY(Kd,Ed)
      Vdx=Xd*Vdc
      Vdy=Yd*Vdc
      Vdz=Zd*Vdc
      CALL LABTRAN(Vdx,Vdy,Vdz, V9Be, Vx,Vy,Vz)
      Vd=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C       Negative value of E(deuteron) to alert Subroutine BANKER.
C ---------------------------------------------------------------------
C     DEUTERON
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      Cx = Xn
      Cy = Yn
      Cz = Zn
      CALL TRANSVEC(Cx,Cy,Cz)
      IKF   = 1000002
      DDENG = EFROMV(Kd,Vd)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      E7Li=Ta-Ed
      V7c=VELOCITY(K7Li,E7Li)
      V7x=-Xd*V7c
      V7y=-Yd*V7c
      V7z=-Zd*V7c
      CALL LABTRAN(V7x,V7y,V7z, V9Be, Vx,Vy,Vz)
      V7=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C     That gets the 7-Li ion energy.
C ---------------------------------------------------------------------
C     7Li
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      CALL TRANSVEC(Cx,Cy,Cz)
      IKF   = 3*1000000 + 7
      DDENG = EFROMV(K7Li,V7)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      Return
      END
C =====================================================================
C =====================================================================
C     LIST OF ROUTIENS IN THIS FILE
C ---------------------------------------------------------------------
C     S: ap7li
C     S: ex8li
C     S: exli8
C     S: pp8li
C     S: tenbdk
C =====================================================================
C =====================================================================
      Subroutine AP7LI(Ta,V11B,Eg)
C     Purpose is to follow  11-B --> 7-Li + alpha  decay, including
C     any subsequent 7-Li excited ion decay.

      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
!$OMP THREADPRIVATE(/VECTOR/)
C     Xn,Yn,Zn at entry are direction cosines of 11-B ion in
C     "detector" coordinates.
C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------
      Common /MASSES/ Emass(21)
      Common /NEUTRN/ Eneut, U,V,W
!$OMP THREADPRIVATE(/NEUTRN/)
      Common /NEUTR2/ Eneut2, U2,V2,W2
!$OMP THREADPRIVATE(/NEUTR2/)
      Common /EXC7LI/ E7Li(6),Li6pn
C     See Subroutine ELI8DK for values in EXC7LI common.

      Common /EXC6LI/ E6Li(6),He4pd
C     See Subroutine TENBDK for values in EXC6LI common.
      real*8 UNIRN

      Data Ka/3/, K7Li/17/, Kn/1/, K6Li/18/

      Eax=RCKE(Ka,K7Li,Ta)
      Eamin=0.0
      Ndx=0
      IF (Eneut .LE. 0.0) goto 5

C     Program counter to here means call to AP7LI came from NPX subroutine.
C     Otherwise call came from subroutine ND.
      Ndx=1
C     Check to see if there is a 2nd neutron waiting to be processed.  If so
C     then fix 11-B decay only to n-stable states of 7-Li.

      IF (Eneut2 .LE. 0.0) goto 5
C     Fix 11-B ion decay by limiting minimum value of energy of out-
C     going alpha particle.

      Eamin=Ta-Li6pn
    5 Temp=4.0

C     See comments just prior to statement 38 in subroutine NN3ALF for
C     estimating "alpha" continuum.
    6 Eatry=CHOOSP(Eax,Temp)
      IF (Eatry .LT. Eamin) goto 6
      Try=Eatry*(Emass(Ka) + Emass(K7Li))/Emass(K7Li)
      Etry=RCKE(Ka,K7Li,Try)
      Try=Try*Eatry/Etry
      Excit=Ta-Try

C     Now, maybe, set  -EXCIT-  to specific energy level in 7-Li.
      Level=7
      IF (Excit .GT. 11.0) goto 10
      K=1
      Do 8 I=2,6
      IF (Excit .LT. E7Li(I)) goto 9
    8 K=K+1
    9 Excit=E7Li(K)
      Level=K
   10 Taex=Ta-Excit
      Eaa=RCKE(Ka,K7Li,Taex)
      CALL RVECT(Xp,Yp,Zp)
      Vaa=VELOCITY(Ka,Eaa)
      Vax=Xp*Vaa
      Vay=Yp*Vaa
      Vaz=Zp*Vaa
      CALL LABTRAN(Vax,Vay,Vaz, V11B, Vx,Vy,Vz)
      Valf=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Ealf=EFROMV(Ka,Valf)
C ---------------------------------------------------------------------
C     ALPHA
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      Cx = Xn
      Cy = Yn
      Cz = Zn
      CALL TRANSVEC(Cx,Cy,Cz)
      IKF   = 2000004
      DDENG = Ealf
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
      Xn = Cx
      Yn = Cy
      Zn = Cz
C ---------------------------------------------------------------------
C     That gets the alpha information.

C     Now get 7-Li ion information
      E7Lic=Taex-Eaa
      V7Lic=VELOCITY(K7Li,E7Lic)
      V7x=-Xp*V7Lic
      V7y=-Yp*V7Lic
      V7z=-Zp*V7Lic
      CALL LABTRAN(V7x,V7y,V7z, V11B, Vx,Vy,Vz)
      V7Li=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      IF (Level .GE. 3) goto 15

C ---------------------------------------------------------------------
C     Program to here have bound state of 7-Li to tidy up.
      Eg=Excit
C ---------------------------------------------------------------------
C     7Li
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      Cx = Xn
      Cy = Yn
      Cz = Zn
      CALL TRANSVEC(Cx,Cy,Cz)
      IKF   = 3*1000000 + 7
      DDENG = EFROMV(K7Li,V7Li)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      Return
C ---------------------------------------------------------------------

C     Next: decide 7-Li particle decay mode.
   15 IF (Level .GE. 5) goto 20
C     If program to here have  7-Li --> alpha + triton  decay mode.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Cx=Xn
      Cy=Yn
      Cz=Zn
      CALL TRANSVEC(Cx,Cy,Cz)

      CALL ALFAPT(Excit,V7Li)

      Return
C ---------------------------------------------------------------------

C     Next is for  7-Li --> n + 6-Li  decay mode.  First get 7-Li motion
C     into detector coordinates.
   20 CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
C     Done. (Xn,Yn,Zn) are the direction cosine of the 7Li in the detector lab system.

C     Now get neutron information.
      Ta=Excit-Li6pn
      Enmax=RCKE(Kn,K6Li,Ta)
      Ehat=Enmax+5.0
      Fn=0.065 + 0.001*Ehat
      Temp=Fn*Ehat
      Eloww=0.0
      Enn=CHOOSN(Eloww,Enmax,Temp)
      Try=Enn*(Emass(Kn) + Emass(K6Li))/Emass(K6Li)
      Entry=RCKE(Kn,K6Li,Try)
      Try=Try*Enn/Entry
      Excitn2=Ta-Try
C     Maybe pair -EXCITN2-  with energy of a level in 6-Li
      Level=7
      IF (Excitn2 .GT. 20.0) goto 25
      K=1
      Do 22 I=2,6
      IF (Excitn2 .LT. E6Li(I)) goto 23
   22 K=K+1
   23 Excitn2=E6Li(K)
      Level=K
      IF (K .EQ. 3) Eg=Excitn2

   25 Taex=Ta-Excitn2
      En2c=RCKE(Kn,K6Li,Taex)
      Vn2c=VELOCITY(Kn,En2c)
      CALL RVECT(Xp,Yp,Zp)
      Vn2x=Xp*Vn2c
      Vn2y=Yp*Vn2c
      Vn2z=Zp*Vn2c
      CALL LABTRAN(Vn2x,Vn2y,Vn2z, V7Li, Vx,Vy,Vz)
      Vn2=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      En2=EFROMV(Kn,Vn2)

C     Transform neutron dir. cos. into "detector" coordinates.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
C ---------------------------------------------------------------------
C     2ND NEUTRON
C ---------------------------------------------------------------------
      IKF   = 2112
      DDENG = En2
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

C     Now get 6-Li ion information.
   28 E6Lix=Taex-En2c
      V6Lic=VELOCITY(K6Li,E6Lix)
      V6x=-Xp*V6Lic
      V6y=-Yp*V6Lic
      V6z=-Zp*V6Lic
      CALL LABTRAN(V6x,V6y,V6z, V7Li, Vx,Vy,Vz)
      V6Li=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C     See Subroutine TENBDK for decay characteristics of levels of 6-Li

C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      CALL TRANSVEC(Zx,Zy,Zz)
C     (Xn,Yn,Zn) are the direction cosine of the 6Li in the detector lab system.
C ---------------------------------------------------------------------

      IF (Level .EQ. 2) goto 30
      IF (Level .GE. 4) goto 30

C     Program to here, bound state of 6-Li.  Tidy up.
C ---------------------------------------------------------------------
C     6Li
C ---------------------------------------------------------------------
      IKF   = 3*1000000 + 6
      DDENG = EFROMV(K6Li,V6Li)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      Return
C ---------------------------------------------------------------------

C     Last section -- 6-Li --> d + alpha
   30 CALL ALFAPD(Excitn2,V6Li)
      Return
      END
C =====================================================================
C =====================================================================
      Function EX8LI(Exc12B)
C     Purpose: for alpha decay of highly-excited 12-B, determine 8-Li
C     excitation energy.  Present programming mocks up TNG results.
C
C     Contains routines involving "decay" of Li ions.

      Common /EXC8LI/ E8Li(10)
      real*8 UNIRN

      Dimension Ex(28),P(10),P1(28),P2(28),P3(28),P4(28),P5(28),
     x    P6(28),P7(28),P8(28),P9(28),P10(28)

C     Next array represents level structure of 8-Li:
      Data E8Li/0.0,0.981,2.255, 3.21, 5.4, 6.1, 6.53, 7.1, 9.0, 10.82/
      Data N/28/, Nterp/1/, Temp/4.0/, Li8pa/10.002/

C     Next arrays relate probabilities for population of individual levels
C     vs excitation energy (Ex) in 12-B.
      Data Ex/ 11.0,  11.4,  12.3,  13.3,  14.2,  15.0,  15.5,  16.5,
     a  17.5,  18.5,  19.5,  20.0,  21.0,  22.0,  23.0,  24.0,  25.0,
     b  26.0,  27.0,  28.0,  30.0,  32.0,  34.0,  36.0,  38.0,  40.0,
     c  44.0,  48.0/
      Data P1/  1.0,  .9815, .985,  .8606, .6478, .574,  .513,  .41,
     a  .341,  .288,  .3448, .338,  .315,  .28,   .25,   .226,  .215,
     b  .2066, .2004, .186,  .1582, .141,  .1235, .1087, .0933, .089,
     c  .079,  .0665/
      Data P2/  1.0,   1.0,   1.0,  .9996, .9808, .899,  .813,  .66,
     a  .553,  .483,  .5268, .518,  .49,   .45,   .415,  .388,  .373,
     b  .3616, .3514, .333,  .3002, .275,  .2505, .2257, .1983, .179,
     c  .149,  .1115/
      Data P3/  1.0,   1.0,   1.0,   1.0,  .9998, .989,  .983,  .95,
     a  .893,  .783,  .7798, .758,  .71,   .66,   .616,  .583,  .56,
     b  .5416, .5264, .502,  .4572, .42,   .3845, .3507, .3133, .284,
     c  .235,  .1815/
      Data P4/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  .984,
     a  .968,  .903,  .9048, .888,  .838,  .7843, .733,  .6943, .665,
     b  .6406, .6204, .591,  .5382, .494,  .4535, .4147, .3733, .338,
     c  .276,  .2115/
      Data P5/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a  .999,  .988,  .9798, .958,  .918,  .8773, .838,  .8093, .786,
     b  .7656, .7454, .715,  .6532, .595,  .5405, .4977, .4433, .4,
     c  .322,  .2445/
      Data P6/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,  .994,  .9898, .973,  .947,  .9233, .903,  .8833, .866,
     b  .8486, .8284, .796,  .7294, .666,  .6075, .5557, .5033, .455,
     c  .366,  .2785/
      Data P7/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,   1.0,  .9998, .988,  .97,   .9553, .943,  .9313, .921,
     b  .9086, .8924, .863,  .7954, .73,   .6675, .6127, .5583, .507,
     c  .414,  .3235/
      Data P8/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,   1.0,   1.0,  .998,  .99,   .9843, .978,  .9723, .9665,
     b  .9576, .9424, .913,  .8444, .777,  .7115, .6527, .5953, .542,
     c  .4433, .3475/
      Data P9/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,   1.0,   1.0,   1.0,   1.0,  .9998, .9985, .9973, .9947,
     b  .9886, .9754, .9472, .8804, .8152, .7495, .6897, .6303, .5762,
     c  .4733, .3755/
      Data P10/ 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  .9999, .9995, .998,
     b  .993,  .981,  .954,  .889,  .825,  .76,   .7,    .64,   .585,
     c  .48,   .38/

      E=Exc12B
      P(1) =EXTERP(Ex, P1,E,N,Nterp)
      P(2) =EXTERP(Ex, P2,E,N,Nterp)
      P(3) =EXTERP(Ex, P3,E,N,Nterp)
      P(4) =EXTERP(Ex, P4,E,N,Nterp)
      P(5) =EXTERP(Ex, P5,E,N,Nterp)
      P(6) =EXTERP(Ex, P6,E,N,Nterp)
      P(7) =EXTERP(Ex, P7,E,N,Nterp)
      P(8) =EXTERP(Ex, P8,E,N,Nterp)
      P(9) =EXTERP(Ex, P9,E,N,Nterp)
      P(10)=EXTERP(Ex,P10,E,N,Nterp)

      Pran=sngl(UNIRN(dummy))
      Nb=1
      Do 1 J=1,10
      IF (Pran .LE. P(J)) goto 2
    1 Nb=Nb+1
    2 IF (Nb .EQ. 11) goto 3
      Ex8Li=E8Li(Nb)
      Return

C     Program to here--highly excited 8-Li, E(level) from "continuum"
    3 Ex8Li=E8Li(10)
      Emax=E-22.5
      IF (Emax .LE. 0.0) Return
      Eatry=CHOOSP(Emax,Temp)
      Ex8Li=E-Eatry-Li8pa
      Return
      END
C =====================================================================
C =====================================================================
      Function EXLI8(Exc9Be)
C     Purpose: for proton decay of highly-excited 9-Be, determine 8-Li
C     excitation energy...program mocks up TNG results.

      real*8 UNIRN
      Common /EXC8LI/ E8Li(10)
C     E8Li array in Function EX8LI
      Dimension Ex(28),P(10),P1(28),P2(28),P3(28),P4(28),P5(28),
     x     P6(28),P7(28),P8(28),P9(28),P10(28)
      Data N/28/, Nterp/1/, Li8pp/16.888/

C     Next arrays relate probabilities for population of individual
C     levels vs excitation energy (Ex) in 9-Be.
      Data Ex/ 17.9,  18.25, 19.2,  19.6,  20.15, 20.6,  22.35, 22.7,
     a  23.05, 23.45, 24.15, 25.0,  25.9,  26.8,  27.7,  28.7,  30.0,
     b  32.0,  34.0,  36.0,  28.0,  40.0,  42.0,  44.0,  46.0,  48.0,
     c  50.0,  55.0/
      Data P1/  1.0,  .99,   .67,   .629,  .585,  .551,  .514,  .497,
     a  .417,  .388,  .3367, .323,  .31,   .3042, .3042, .2822, .28,
     b  .236,  .195,  .151,  .113,  .08,   .056,  .038,  .025,  .0165,
     c  .0106, .003/
      Data P2/  1.0,   1.0,   1.0,  .989,  .97,   .936,  .834,  .795,
     a  .697,  .655,  .5667, .525,  .483,  .4472, .4282, .3842, .358,
     b  .288,  .2295, .273,  .127,  .089,  .0616, .0415, .0271, .0177,
     c  .0113, .003/
      Data P3/  1.0,   1.0,   1.0,   1.0,   1.0,  .991,  .93,   .908,
     a  .885,  .857,  .7847, .745,  .697,  .6472, .6102, .5442, .493,
     b  .386,  .2985, .22,   .158,  .109,  .0743, .0492, .0318, .0204,
     c  .0128, .003/
      Data P4/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  .998,
     b  .995,  .989,  .9687, .942,  .894,  .8352, .7702, .6962, .623,
     c  .484,  .3695, .27,   .193,  .1325, .0909, .0597, .0388, .0248,
     d  .0156, .0035/
      Data P5/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,  .9976, .9917, .98,   .952,  .9142, .8672, .8062, .741,
     b  .594,  .4564, .335,  .241,  .1665, .1149, .0762, .0503, .0328,
     c  .021,  .0054/
      Data P6/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,   1.0,  .998,  .992,  .977,  .9562, .9292, .8922, .841,
     b  .697,  .5485, .407,  .294,  .2045, .1419, .0947, .0628, .0416,
     c  .0268, .0074/
      Data P7/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,   1.0,   1.0,  .998,  .992,  .9822, .9687, .9462, .908,
     b  .769,  .6135, .46,   .335,  .2363, .1659, .1127, .0758, .0511,
     c  .0336, .0104/
      Data P8/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,   1.0,   1.0,   1.0,   1.0,  .9997, .9967, .9882, .966,
     b  .8375, .6815, .518,  .381,  .2713, .1916, .1307, .0878, .0591,
     c  .0386, .012/
      Data P9/  1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  .9997, .991,
     b  .8855, .7395, .576,  .431,  .3108, .2216, .1532, .1044, .0711,
     c  .0471, .0153/
      Data P10/ 1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,
     a   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,   1.0,  .992,
     b  .8917, .7532, .5925, .447,  .3252, .2336, .1627, .1114, .0759,
     c  .0501, .0163/

      E=Exc9Be
      P(1) =EXTERP(Ex, P1,E,N,Nterp)
      P(2) =EXTERP(Ex, P2,E,N,Nterp)
      P(3) =EXTERP(Ex, P3,E,N,Nterp)
      P(4) =EXTERP(Ex, P4,E,N,Nterp)
      P(5) =EXTERP(Ex, P5,E,N,Nterp)
      P(6) =EXTERP(Ex, P6,E,N,Nterp)
      P(7) =EXTERP(Ex, P7,E,N,Nterp)
      P(8) =EXTERP(Ex, P8,E,N,Nterp)
      P(9) =EXTERP(Ex, P9,E,N,Nterp)
      P(10)=EXTERP(Ex,P10,E,N,Nterp)

      Pran=sngl(UNIRN(dummy))
      M=1
      Do 1 J=1,10
      IF (Pran .LE. P(J)) goto 2
    1 M=M+1
    2 IF (M .EQ. 11) goto 3
      ExLi8=E8Li(M)
      Return

C     Or get highly-excited 8-Li, E(Level) from "continuum"
    3 ExLi8=E8Li(10)
      Emax=E-29.4
      IF (Emax .LE. 0.0) Return
      Ehat=Emax+5.0
      F=0.1245 + 0.001*ABS(Ehat-45.0)
      Temp=F*Ehat
      Ep=CHOOSP(Emax,Temp)
      ExLi8=E-Ep-Li8pp
      Return
      END
C =====================================================================
C =====================================================================
      Subroutine PP8LI(Excit,V9Be,Eg)
C     Next portion of programming for 9-Be --> p + 8-Li decay.

      Common /MASSES/ Emass(21)
      Common /NEUTRN/ Eneut, U,V,W
!$OMP THREADPRIVATE(/NEUTRN/)
      Common /VECTOR/ Xb,Yb,Zb, Xn,Yn,Zn
!$OMP THREADPRIVATE(/VECTOR/)
C     At entry, Xn,Yn,Zn are dir. cos. of 9-Be ion motion in "detector"
C     coordinates.
C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------
      real*8 UNIRN

      Data Li8pp/16.888/, Li7pn/2.033/, Kp/2/, K8Li/16/

      Ta=Excit-Li8pp
C     Test for possible alert to this routine to do only proton decay,
C     i.e., no further decay of 8-Li to n + 7-Li.

      IF (V9Be .GT. 0.0) goto 4
      V9Be=-V9Be
      Excitp=0.0
C     Choose between ground state and first-excited state of 8-Li.
      IF (sngl(UNIRN(dummy)) .LT. 0.3) Excitp=0.9808
      Goto 6

C     Next is for no limits on 8-Li excitation energy.
    4 Excitp=EXLI8(Excit)
C     That gets 8-Li excitation energy.

      Ix=1
      IF (Excitp .GT. Li7pn) goto 7
    6 Eg=Excitp
      Ix=0
    7 Taex=Ta-Excitp
      Epp=RCKE(Kp,K8Li,Taex)
      Vpp=VELOCITY(Kp,Epp)
      CALL RVECT(Xp,Yp,Zp)
      Vpx=Xp*Vpp
      Vpy=Yp*Vpp
      Vpz=Zp*Vpp
      CALL LABTRAN(Vpx,Vpy,Vpz, V9Be, Vx,Vy,Vz)
      Vprot=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Eprot=EFROMV(Kp,Vprot)

      CALL DIRCOS(Vx,Vy,Vz, Xb,Yb,Zb)
      Zx=Xn
      Zy=Yn
      Zz=Zn
      CALL TRANSVEC(Zx,Zy,Zz)
C ---------------------------------------------------------------------
C     PROTON
C ---------------------------------------------------------------------
      IKF   = 2212
      DDENG = Eprot
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

C     Now get 8-Li ion:
      E8Li=Taex-Epp
      V8Lic=VELOCITY(K8Li,E8Li)
      V8x=-Xp*V8Lic
      V8y=-Yp*V8Lic
      V8z=-Zp*V8Lic
      CALL LABTRAN(V8x,V8y,V8z, V9Be, Vx,Vy,Vz)
      V8Li=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      IF (Ix .EQ. 1) goto 2
C     If Ix = 0, then no additional neutron decay. So tidy up and exit.

C ---------------------------------------------------------------------
C     8LI
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz, Xb,Yb,Zb)
      CALL TRANSVEC(Zx,Zy,Zz)
      IKF   = 3*1000000+8
      DDENG = EFROMV(K8Li,V8Li)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      Return
C ---------------------------------------------------------------------

C     Set up for 8-Li decay.
    2 CALL DIRCOS(Vx,Vy,Vz, Xb,Yb,Zb)
      CALL TRANSVEC(Zx,Zy,Zz)
C     Set up done.  Get 8-Li decay.
      CALL ELI8DK(Excitp,V8Li,Egamma)
      Eg=Egamma
C     That's that!

      Return
      END
C =====================================================================
C =====================================================================
      Subroutine TENBDK(Excit,V10B,Eg)
C     Purpose is to follow the decay of particle-unstable excited
C     states of 10-B.  This routine is used in Subroutines ND, NT,
C     and NPX.
C
C     Choices of decay modes are:
C     (a)   10-B --> p + 9-Be (ground state)
C     or,
C     (a')  10-B --> d + 8-Be
C     (b')   8-Be--> alpha + alpha
C     or,
C     (a")  10-B --> alpha + 6-Li
C     (b")   6-Li--> gamma + 6-Li (ground state)
C     or,
C     (b'")  6-Li--> d + alpha
C
C     or, under certain circumstances,
C     (a3)  10-B --> p + 9-Be (excited)
C     (b3)   9-Be--> n + 8-Be
C     (c3)   8-Be--> alpha + alpha
C
C     or, fractionally,
C     (b3')  9-Be--> d + 7-Li
C ---------------------------------------------------------------------
      Common /MASSES/ Emass(21)
      real*8 UNIRN
      Common /VECTOR/ Xpn,Ypn,Zpn, Xn,Yn,Zn
!$OMP THREADPRIVATE(/VECTOR/)
C     At entry to this routine Xn,Yn,Zn are direction cosines of
C     moving 10-B ion in "detector" lab. coordinates.

      Common /NEUTRN/ Eneut, Dc1,Dc2,Dc3
!$OMP THREADPRIVATE(/NEUTRN/)
      Common /EXC6LI/ ExLi6(6),He4pd
      Common /BENINE/ Ex9Be(6),Be8pn
C     See Subroutine N3He for values in BENINE common.
C ---------------------------------------------------------------------
C     Modified on February 25, 2025
C ---------------------------------------------------------------------
      common /PARTICLE/ IDRTYP, KF(1:20), Npart,
     &        DENG(1:20), DVX(1:20), DVY(1:20), DVZ(1:20)
      real(8)   :: DENG, DVX, DVY, DVZ
      integer(4):: IDRTYP, KF, Npart
!$OMP THREADPRIVATE(/PARTICLE/)
C ---------------------------------------------------------------------

C     Next array represents level structure of 6-Li.
      Data ExLi6/0.0,2.185, 3.563, 4.31, 5.37, 5.65/

      Data Kp/2/, K6Li/18/, Kd/10/, Ka/3/, K9Be/4/, Kn/1/, K8Be/8/
      Data Qnt/18.93/, Be8pd/6.026/, Be9pp/6.585/, Be8aa/0.092/
      Data He4pd/1.474/, Li6pa/4.46/, Ntrp/1/

      Eg=0.0
      Ta=Excit
C     Choose which decay mode to follow.

      IF (Excit .LE. Be8pd) goto 10
      IF (Excit .LE. Be9pp) goto 5

C     Such experimental information as exists suggests that about half
C     of the highly-excited levels in 10-B decay predominantly by
C         proton emission.
      IF (sngl(UNIRN(dummy)) .GE. 0.5) goto 5

C     Study p + 9-Be mode first.
      Taex=Ta-Be9pp
      Try=Taex
      Ep=RCKE(Kp,K9Be,Taex)
      Npgo=-1
      Excit=0.0
      IF (sngl(UNIRN(dummy)) .LE. 0.04) goto 4
C     Last is to slightly enhance 9-Be ground state -- 4% is ad hoc.

C     Next test variable -ENEUT-; if it's zero then can set up for
C     possible 9-Be --> n + 8-Be calculation.

      IF (Eneut .GT. 0.0) goto 4
C     -ENEUT- is zero, so look for decay mode (a3), above.

      Ehat=Ep+5.0
      F=0.1245 + 0.001*ABS(Ehat-45.0)
      T=F*Ehat
      Eprot=CHOOSP(Ep,T)
      Try=Eprot*(Emass(Kp) + Emass(K9Be))/Emass(K9Be)
      Eptry=RCKE(Kp,K9Be,Try)
      Try=Try*Eprot/Eptry
      Excit=Taex-Try
C     (No test for 9-Be --> p + 8-Li as already have one proton's info.)

C     Next check possible  9-Be --> d + 7-Li decay mode.
      Rani=sngl(UNIRN(dummy))
      IF (Excit .LT. 17.4 .OR. Rani .GT. 0.09) goto 2
      Npgo=1
      Goto 9
    2 Npgo=0
      IF (Excit .GE. 16.0) goto 9
      K=1
      Do 6 I=2,6
      IF (Excit .LT. Ex9Be(I)) goto 7
    6 K=K+1
    7 Excit=Ex9Be(K)
C     Recheck for 9-Be ground state at this point.
      IF (K .EQ. 1) Npgo=-1
      Try=Taex-Excit
    9 Ep=RCKE(Kp,K9Be,Try)
    4 Vp=VELOCITY(Kp,Ep)
C     "Fission" of 10-B is isotropic in center of mass coords.
      CALL RVECT(Xp,Yp,Zp)
      Vpx=Xp*Vp
      Vpy=Yp*Vp
      Vpz=Zp*Vp
      CALL LABTRAN(Vpx,Vpy,Vpz, V10B, Vx,Vy,Vz)
      Vpl=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Epl=EFROMV(Kp,Vpl)

      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Za=Xn
      Zb=Yn
      Zc=Zn
      CALL TRANSVEC(Za,Zb,Zc)
C ---------------------------------------------------------------------
C     PROTON
C ---------------------------------------------------------------------
      IKF   = 2212
      DDENG = Epl
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

C     Next: Get information for 9-Be ion.
    3 E9Bec=Try-Ep
      V9Bec=VELOCITY(K9Be,E9Bec)
      V9x=-Xp*V9Bec
      V9y=-Yp*V9Bec
      V9z=-Zp*V9Bec
      CALL LABTRAN(V9x,V9y,V9z, V10B, Vx,Vy,Vz)
      V9Be=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C
      IF (Npgo .EQ. -1) goto 13
C ---------------------------------------------------------------------
      IF (Npgo .EQ. 1) THEN
C       Program counter to here have 9-Be --> d + 7-Li decay mode
C       Rotate 9-Be ion motion into "detector" corrdinates:
        CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
        CALL TRANSVEC(Za,Zb,Zc)
C       Now (Xn, Yn, Zn) are dir cos components of 9Be in detector lab sys
        goto 25
      END IF
C ---------------------------------------------------------------------

C     Program counter to here have 9-Be --> n + 8-Be decay mode
C     First rotate 9-Be ion motion into "detector" coordinates:
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Za,Zb,Zc)
C     Done.

      Tax=Excit-Be8aa-Be8pn
      Enc=RCKE(Kn,K8Be,Tax)
      Vnc=VELOCITY(Kn,Enc)
      CALL RVECT(Xp,Yp,Zp)
      Vnx=Xp*Vnc
      Vny=Yp*Vnc
      Vnz=Zp*Vnc
      CALL LABTRAN(Vnx,Vny,Vnz, V9Be, Vx,Vy,Vz)
      Vnlab=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
      Enn=EFROMV(Kn,Vnlab)
C     That's the neutron's energy.  Store it in /NEUTRN/ common area.
      Eneut=Enn
C     Rotate neutron's motion into "detector" coordinates and store them.
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      Za=Xn
      Zb=Yn
      Zc=Zn
      CALL TRANSVEC(Za,Zb,Zc)
      Dc1=Xn
      Dc2=Yn
      Dc3=Zn
C ---------------------------------------------------------------------
C     NEUTRON
C ---------------------------------------------------------------------
      IKF   = 2112
      DDENG = Eneut
      DDVX  = Dc1
      DDVY  = Dc2
      DDVZ  = Dc3
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
      Cx = Za
      Cy = Zb
      Cz = Zc
C     (Cx, Cy, Cz) are dir cos components of 9Be in detector lab sys
C ---------------------------------------------------------------------

C     Now get 8-Be ion motion.
      E8Bec=Tax-Enc
      V8Bec=VELOCITY(K8Be,E8Bec)
      V8x=-Xp*V8Bec
      V8y=-Yp*V8Bec
      V8z=-Zp*V8Bec
      CALL LABTRAN(V8x,V8y,V8z, V9Be, Vx,Vy,Vz)
      V8Be=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C     Done.  Now get  8-Be --> 2 alphas.
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz, Xpn,Ypn,Zpn)
      CALL TRANSVEC(Cx,Cy,Cz)
C     (Xn, Yn, Zn) are dir cos components of 8Be in detector lab sys
C ---------------------------------------------------------------------
      Iechrg=3
      Goto 16
C ---------------------------------------------------------------------

C     Next it to tidy up  10-B --> p + 9-Be (ground state) decay mode.
   13 E9Be=EFROMV(K9Be,V9Be)
C ---------------------------------------------------------------------
C     9Be
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      CALL TRANSVEC(Za,Zb,Zc)
      IKF   = 4*1000000 + 9
      DDENG = E9Be
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      Return
C     END SECTION on p + 9-Be
C ---------------------------------------------------------------------

C     Choosing between the other two modes of 10-B decay appears a toss up.
    5 IF (sngl(UNIRN(dummy)) .GE. 0.5) goto 10

C ---------------------------------------------------------------------
C     Do 10-B --> d + 8-Be decay mode next.  Ignore 8-Be excited states
C     as we are already high in (effective) 12-C excitation, and these
C     branches of decay are relatively minor.

      Taex=Ta-Be8pd
      Ed=RCKE(Kd,K8Be,Taex)
      Vd=VELOCITY(Kd,Ed)
      CALL RVECT(Xd,Yd,Zd)
      Vdx=Xd*Vd
      Vdy=Yd*Vd
      Vdz=Zd*Vd
      CALL LABTRAN(Vdx,Vdy,Vdz, V10B, Vx,Vy,Vz)
      Vdl=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C ---------------------------------------------------------------------
C     DEUTERON
C ---------------------------------------------------------------------
C     (Xn, Yn, Zn) are dir cos components of 10B in detector lab sys
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      Cx = Xn
      Cy = Yn
      Cz = Zn
      CALL TRANSVEC(Cx,Cy,Cz)
      IKF   = 1000002
      DDENG = EFROMV(Kd,Vdl)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
      Xn = Cx
      Yn = Cy
      Zn = Cz
C     (Cx, Cy, Cz) are dir cos components of 10B in detector lab sys
C ---------------------------------------------------------------------

      E8Bec=Taex-Ed
      V8Bec=VELOCITY(K8Be,E8Bec)
      V8x=-Xd*V8Bec
      V8y=-Yd*V8Bec
      V8z=-Zd*V8Bec
      CALL LABTRAN(V8x,V8y,V8z, V10B, Vx,Vy,Vz)
      V8Be=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)

C ---------------------------------------------------------------------
C     Now get 8-Be --> alpha + alpha(2)
      Iechrg=4
   16 Ta=Be8aa
      Ea=0.5*Ta
C     Alphas share the available energy.
      Va=VELOCITY(Ka,Ea)
      CALL RVECT(Xc,Yc,Zc)
      Vax=Xc*Va
      Vay=Yc*Va
      Vaz=Zc*Va
      CALL LABTRAN(Vax,Vay,Vaz, V8Be, Vx,Vy,Vz)
      Valph=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C ---------------------------------------------------------------------
C     ALPHA
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      Cx = Xn
      Cy = Yn
      Cz = Zn
      CALL TRANSVEC(Cx,Cy,Cz)
      IKF   = 2000004
      DDENG = EFROMV(Ka,Valph)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------

C     Recall symmetry so that Va2x=-Vax, etc, for the other
C     alpha.
      Va2x=-Vax
      Va2y=-Vay
      Va2z=-Vaz
      CALL LABTRAN(Va2x,Va2y,Va2z, V8Be, Vx,Vy,Vz)
      Valph2=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C ---------------------------------------------------------------------
C     ALPHA
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      CALL TRANSVEC(Cx,Cy,Cz)
      IKF   = 2000004
      DDENG = EFROMV(Ka,Valph2)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      Return
C     END SECTION for 10-B --> d + 8-Be
C             and for 10-B --> p + n + 8-Be.
C ---------------------------------------------------------------------

C ---------------------------------------------------------------------
C     That leaves 10-B --> alpha + 6-Li decay modes to do.
   10 Taex=Ta-Li6pa
      Eax=RCKE(Ka,K6Li,Taex)
      Eg=0.0
      Excita=0.0
      K=1
      IF (sngl(UNIRN(dummy)) .LE. 0.03) goto 18
C     As done above, a slight (3%, ad hoc) enhancement of ground-state.

      Temp=4.0
C     See comments for alpha "continuum" computation just prior to statement
C     no. 38 in subroutine NN3ALF.

      Eatry=CHOOSP(Eax,Temp)
      Try=Eatry*(Emass(Ka) + Emass(K6Li))/Emass(K6Li)
      Etry=RCKE(Ka,K6Li,Try)
      Try=Try*Eatry/Etry
      Excita=Taex-Try
C     Now pair -Excita- with a level in 6-Li
      K=7
      IF (Excita .GT. 20.0) goto 18

      K=1
      Do 14 I=2,6
      IF (Excita .LT. ExLi6(I)) goto 15
   14 K=K+1
   15 Excita=ExLi6(K)
C ---------------------------------------------------------------------
C     6-Li level decay scheme is as follows:
C
C               K=1     Ex=0.0          Stable
C               K=2     Ex=2.185        d + alpha (very weak gamma)
C               K=3     Ex=3.563        ground-state gamma ray
C               K=4     Ex=4.31         d + alpha
C               K=5     Ex=5.366        weak g.s. gamma reported, but
C                                       p or n decay most likely.
C               K=6     Ex=5.65         d + alpha
C               K>6     Ex>20           triton + 3-He
C
C     For the present programming, take K=3 gamma decay and the remaining
C     K>1 excited state decay to be d + alpha.
C ---------------------------------------------------------------------
      IF (K .EQ. 3) Eg=Excita

C     Now get information on alpha.
   18 Taex=Taex-Excita
      Eaa=RCKE(Ka,K6Li,Taex)
      Vaa=VELOCITY(Ka,Eaa)
      CALL RVECT(Xa,Ya,Za)
      Vax=Xa*Vaa
      Vay=Ya*Vaa
      Vaz=Za*Vaa
      CALL LABTRAN(Vax,Vay,Vaz, V10B, Vx,Vy,Vz)
      Va=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)
C ---------------------------------------------------------------------
C     ALPHA
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      Cx = Xn
      Cy = Yn
      Cz = Zn
      CALL TRANSVEC(Cx,Cy,Cz)
      IKF   = 2000004
      DDENG = EFROMV(Ka,Va)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
      Xn = Cx
      Yn = Cy
      Zn = Cz
C     (Xn, Yn, Zn) are dir cos components of 10B in detector lab sys
C ---------------------------------------------------------------------

C     Now set up 6-Li ion.
      E6Lic=RCKE(K6Li,Ka,Taex)
      V6Lic=VELOCITY(K6Li,E6Lic)
      V6x=-Xa*V6Lic
      V6y=-Ya*V6Lic
      V6z=-Za*V6Lic
      CALL LABTRAN(V6x,V6y,V6z, V10B, Vx,Vy,Vz)
      V6Li=SQRT(Vx*Vx + Vy*Vy + Vz*Vz)

      IF(K .EQ. 2) goto 20
      IF (K. GE. 4) goto 20

C     For K=1 or 3 tidy up 6-Li energy.
C ---------------------------------------------------------------------
C     6Li
C ---------------------------------------------------------------------
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      Cx = Xn
      Cy = Yn
      Cz = Zn
      CALL TRANSVEC(Cx,Cy,Cz)
      IKF   = 3*1000000 + 6
      DDENG = EFROMV(K6Li,V6Li)
      DDVX  = Xn
      DDVY  = Yn
      DDVZ  = Zn
      CALL CX_score(IKF, DDENG, DDVX, DDVY, DDVZ)
C ---------------------------------------------------------------------
      Return

C ---------------------------------------------------------------------
C     Penultimate section: 6-Li --> d + alpha
C ---------------------------------------------------------------------
   20 CONTINUE
      CALL DIRCOS(Vx,Vy,Vz,Xpn,Ypn,Zpn)
      Cx = Xn
      Cy = Yn
      Cz = Zn
      CALL TRANSVEC(Cx,Cy,Cz)
C     Now (Xn, Yn, Zn) are dir cos components of 6Li in detector lab sys
C ---------------------------------------------------------------------
      CALL ALFAPD(Excita,V6Li)
      Return

C ---------------------------------------------------------------------
C     Last section:  9-Be --> d + 7-Li  (no deuteron escape considered)
   25 CALL DPLI7(Excit,V9Be)
      Return

      END
C =====================================================================
C =====================================================================
C     LIST OF ROUTIENS IN THIS FILE
C ---------------------------------------------------------------------
C     S: hydrog (2021.05.13)
C     B: hyddat (2021.05.13)
C     B: cardat (2021.05.24)
C =====================================================================
      subroutine hydrog(plab,angcm)
C ---------------------------------------------------------------------
C     This subroutine is called by s:jamangel to return the direction
C     cosine of elastic scatterng nucleons in CM system.
C ---------------------------------------------------------------------
C     (IN) plab:  momentum of incidnet nucleon in lab system (GeV)
C     (OUT)angcm: cosine of angle relative to incident angle (-)
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
      implicit double precision(a-h, o-z)
      parameter(ihemax=70)
      parameter(ihamax=21)
      common/hydcmm/hmu(ihamax), hydrene(ihemax), hydrang(ihemax,ihamax)
      dimension prob(1:ihamax)
      dimension accp(1:ihamax)
C ---------------------------------------------------------------------
C     mass of nucleon (GeV)
      data dm/0.938272046d0/ ! T.Sato 2022/11/03
C ---------------------------------------------------------------------
      angcm = -1.d0
      indx = -1
C ---------------------------------------------------------------------
C     kinetic energy of incident nucleon in lab sys (MeV)
      tlab = ( sqrt(plab**2 + dm**2) - dm ) * 1.0d+3
C ---------------------------------------------------------------------
      do i = 1, ihemax-1
        if( tlab >= hydrene(i) .and. tlab < hydrene(i+1) ) then
          indx=i
          exit
        end if
      end do
      if(indx.eq.-1) then ! energy is too high, T.Sato 2022/11/03
       angcm=-101.0d0 ! index that angcm cannot be determined in this subroutine
       return
      endif
C ---------------------------------------------------------------------
C     probability densities of angular distribution at the energy of
C     tlab (MeV). P(tlab,mu), -1 =< mu =< 1
      prob(:) = (hydrang(indx+1,:) - hydrang(indx,:)) /
     & (hydrene(indx+1) - hydrene(indx)) * (tlab - hydrene(indx))
     & + hydrang(indx,:)
C ---------------------------------------------------------------------
C     integrating the probability densities
      accp(1) = 0.d0
      do j = 2, ihamax
        a = ( prob(j-1) - prob(j) ) / ( hmu(j-1) - hmu(j) )
        b = prob(j) - a*hmu(j)
        pup = 1.d0/2.d0*a*hmu(j)**2 + b*hmu(j)
        plw = 1.d0/2.d0*a*hmu(j-1)**2 + b*hmu(j-1)
        accp(j) = accp(j-1) + (pup - plw)
      end do
C ---------------------------------------------------------------------
      rannum = UNIRN(dummy)
      xprob = rannum*accp(ihamax)
C ---------------------------------------------------------------------
       do j = 1, ihamax-1
         if( xprob >= accp(j) .and. xprob < accp(j+1) ) then
           indx = j
         else if( xprob >= accp(ihamax) ) then
           xprob = accp(ihamax)
           indx = ihamax-1
         end if
       end do
C ---------------------------------------------------------------------
      angcm = (hmu(indx+1)-hmu(indx))/(accp(indx+1)-accp(indx))
     & *(xprob-accp(indx)) + hmu(indx)
C ---------------------------------------------------------------------
      return
      end subroutine hydrog
C ---------------------------------------------------------------------
C =====================================================================
      block data hyddat
C ---------------------------------------------------------------------
C     Data for angular distribution of low-energy np elastic scattering.
C     Numerical values were taken from JENDL/HE-2007.
C ---------------------------------------------------------------------
C     Refference:
C     CHIBA, S. ET AL.: NUCL. SCI. TECHNOL. 33 (1996)654.
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
      implicit double precision(a-h, o-z)
      parameter(ihemax=70)
      parameter(ihamax=21)
      common/hydcmm/hmu(ihamax), hydrene(ihemax), hydrang(ihemax,ihamax)
C ---------------------------------------------------------------------
      data (hmu(j),j=1,ihamax)/
     & -1.00,-0.90,-0.80,-0.70,-0.60,-0.50,-0.40,-0.30,-0.20,-0.10,0.00,
     &  0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00
     & /
C ---------------------------------------------------------------------
      data (hydrene(i),i=1,ihemax)/
     & 16.0,18.0,20.0,25.0,30.0,35.0,40.0,45.0,50.0,55.0,
     & 60.0,65.0,70.0,75.0,80.0,85.0,90.0,95.0,100.0,110.0,
     & 120.0,130.0,140.0,150.0,160.0,170.0,180.0,190.0,200.0,220.0,
     & 240.0,260.0,280.0,300.0,320.0,340.0,360.0,380.0,400.0,420.0,
     & 440.0,460.0,480.0,500.0,520.0,540.0,560.0,580.0,600.0,620.0,
     & 640.0,660.0,680.0,700.0,720.0,740.0,760.0,780.0,800.0,820.0,
     & 840.0,860.0,880.0,900.0,920.0,940.0,960.0,980.0,1000.0,3000.0
     & /
C ---------------------------------------------------------------------
C     [  1]     16.0  MeV
      data (hydrang(  1,j),j=1,ihamax)/
     & 0.567004,0.552002,0.539494,0.529109,0.520509,0.513390,0.507478,
     & 0.502535,0.498353,0.494761,0.491617,0.488814,0.486278,0.483966,
     & 0.481871,0.480016,0.478459,0.477291,0.476633,0.476644,0.477511
     & /
C ---------------------------------------------------------------------
C     [  2]     18.0  MeV
      data (hydrang(  2,j),j=1,ihamax)/
     & 0.579361,0.560043,0.544187,0.531263,0.520786,0.512324,0.505491,
     & 0.499953,0.495423,0.491667,0.488495,0.485772,0.483408,0.481364,
     & 0.479650,0.478327,0.477503,0.477336,0.478034,0.479853,0.483101
     & /
C ---------------------------------------------------------------------
C     [  3]     20.0  MeV
      data (hydrang(  3,j),j=1,ihamax)/
     & 0.559543,0.539790,0.527195,0.518510,0.512071,0.506927,0.502511,
     & 0.498516,0.494827,0.491476,0.488603,0.486413,0.485125,0.484888,
     & 0.485691,0.487274,0.489133,0.490689,0.491709,0.492986,0.497098
     & /
C ---------------------------------------------------------------------
C     [  4]     25.0  MeV
      data (hydrang(  4,j),j=1,ihamax)/
     & 0.576326,0.548767,0.531991,0.520501,0.511761,0.504519,0.498153,
     & 0.492407,0.487280,0.482948,0.479704,0.477890,0.477810,0.479604,
     & 0.483104,0.487728,0.492536,0.496610,0.499915,0.504678,0.517061
     & /
C ---------------------------------------------------------------------
C     [  5]     30.0  MeV
      data (hydrang(  5,j),j=1,ihamax)/
     & 0.592495,0.556629,0.535854,0.521742,0.510723,0.501269,0.492764,
     & 0.485067,0.478338,0.472942,0.469362,0.468113,0.469606,0.473972,
     & 0.480846,0.489234,0.497641,0.504723,0.510716,0.519760,0.542848
     & /
C ---------------------------------------------------------------------
C     [  6]     35.0  MeV
      data (hydrang(  6,j),j=1,ihamax)/
     & 0.609196,0.564185,0.539450,0.522820,0.509486,0.497656,0.486782,
     & 0.476896,0.468366,0.461774,0.457825,0.457233,0.460538,0.467860,
     & 0.478605,0.491292,0.503766,0.514170,0.523076,0.536975,0.572885
     & /
C ---------------------------------------------------------------------
C     [  7]     40.0  MeV
      data (hydrang(  7,j),j=1,ihamax)/
     & 0.627561,0.572229,0.543437,0.524317,0.508577,0.494158,0.480645,
     & 0.468299,0.457730,0.449764,0.445347,0.445408,0.450634,0.461138,
     & 0.476071,0.493404,0.510227,0.524093,0.535957,0.555061,0.605588
     & /
C ---------------------------------------------------------------------
C     [  8]     45.0  MeV
      data (hydrang(  8,j),j=1,ihamax)/
     & 0.648736,0.581568,0.548481,0.526822,0.508523,0.491255,0.474791,
     & 0.459676,0.446794,0.437229,0.432175,0.432789,0.439918,0.453670,
     & 0.472930,0.495069,0.516340,0.533634,0.548324,0.572765,0.639383
     & /
C ---------------------------------------------------------------------
C     [  9]     50.0  MeV
      data (hydrang(  9,j),j=1,ihamax)/
     & 0.673861,0.593002,0.555243,0.530919,0.509853,0.489424,0.469657,
     & 0.451429,0.435924,0.424487,0.418560,0.419532,0.428415,0.445326,
     & 0.468872,0.495788,0.521423,0.541935,0.559140,0.588825,0.672691
     & /
C ---------------------------------------------------------------------
C     [ 10]     55.0  MeV
      data (hydrang( 10,j),j=1,ihamax)/
     & 0.703760,0.607116,0.564207,0.537038,0.512958,0.489031,0.465588,
     & 0.443888,0.425430,0.411819,0.404734,0.405790,0.416175,0.436022,
     & 0.463664,0.495172,0.524937,0.548319,0.567589,0.602251,0.704243
     & /
C ---------------------------------------------------------------------
C     [ 11]     60.0  MeV
      data (hydrang( 11,j),j=1,ihamax)/
     & 0.738004,0.623620,0.575139,0.544982,0.517685,0.489980,0.462554,
     & 0.437086,0.415403,0.399357,0.390850,0.391720,0.403338,0.425871,
     & 0.457384,0.493274,0.526923,0.552836,0.573746,0.613117,0.734016
     & /
C ---------------------------------------------------------------------
C     [ 12]     65.0  MeV
      data (hydrang( 12,j),j=1,ihamax)/
     & 0.775849,0.642004,0.587626,0.554396,0.523747,0.492062,0.460433,
     & 0.430987,0.405880,0.387195,0.377042,0.377478,0.390067,0.415034,
     & 0.450193,0.490257,0.527566,0.555720,0.577908,0.621763,0.762291
     & /
C ---------------------------------------------------------------------
C     [ 13]     70.0  MeV
      data (hydrang( 13,j),j=1,ihamax)/
     & 0.816547,0.661762,0.601254,0.564928,0.530855,0.495068,0.459101,
     & 0.425553,0.396897,0.375429,0.363446,0.363219,0.376525,0.403672,
     & 0.442248,0.486284,0.527052,0.557203,0.580369,0.628528,0.789355
     & /
C ---------------------------------------------------------------------
C     [ 14]     75.0  MeV
      data (hydrang( 14,j),j=1,ihamax)/
     & 0.859347,0.682378,0.615608,0.576222,0.538721,0.498787,0.458437,
     & 0.420746,0.388494,0.364157,0.350199,0.349104,0.362880,0.391952,
     & 0.433709,0.481519,0.525566,0.557514,0.581424,0.633749,0.815483
     & /
C ---------------------------------------------------------------------
C     [ 15]     80.0  MeV
      data (hydrang( 15,j),j=1,ihamax)/
     & 0.903512,0.703347,0.630274,0.587926,0.547056,0.503012,0.458317,
     & 0.416529,0.380705,0.353471,0.337433,0.335286,0.349293,0.380031,
     & 0.424735,0.476124,0.523295,0.556888,0.581369,0.637768,0.840968
     & /
C ---------------------------------------------------------------------
C     [ 16]     85.0  MeV
      data (hydrang( 16,j),j=1,ihamax)/
     & 0.948303,0.724177,0.644852,0.599690,0.555571,0.507528,0.458612,
     & 0.412856,0.373559,0.343460,0.325277,0.321915,0.335919,0.368065,
     & 0.415477,0.470259,0.520422,0.555561,0.580511,0.640939,0.866102
     & /
C ---------------------------------------------------------------------
C     [ 17]     90.0  MeV
      data (hydrang( 17,j),j=1,ihamax)/
     & 0.992965,0.744339,0.658914,0.611156,0.563979,0.512130,0.459208,
     & 0.409697,0.367103,0.334226,0.313872,0.309155,0.322931,0.356225,
     & 0.406101,0.464093,0.517135,0.553760,0.579133,0.643584,0.891159
     & /
C ---------------------------------------------------------------------
C     [ 18]     95.0  MeV
      data (hydrang( 18,j),j=1,ihamax)/
     & 1.036726,0.763309,0.672036,0.621964,0.571989,0.516611,0.459985,
     & 0.407023,0.361382,0.325877,0.303366,0.297176,0.310503,0.344680,
     & 0.396773,0.457790,0.513618,0.551713,0.577523,0.646030,0.916397
     & /
C ---------------------------------------------------------------------
C     [ 19]    100.0  MeV
      data (hydrang( 19,j),j=1,ihamax)/
     & 1.078887,0.780618,0.683833,0.631774,0.579316,0.520753,0.460805,
     & 0.404775,0.356411,0.318482,0.293868,0.286107,0.298776,0.333572,
     & 0.387634,0.451505,0.510057,0.549663,0.576002,0.648654,0.942147
     & /
C ---------------------------------------------------------------------
C     [ 20]    110.0  MeV
      data (hydrang( 20,j),j=1,ihamax)/
     & 1.156540,0.809052,0.702535,0.647638,0.591288,0.527537,0.462250,
     & 0.401391,0.348707,0.306695,0.278167,0.267079,0.277872,0.313159,
     & 0.370437,0.439513,0.503358,0.546178,0.573914,0.655057,0.995249
     & /
C ---------------------------------------------------------------------
C     [ 21]    120.0  MeV
      data (hydrang( 21,j),j=1,ihamax)/
     & 1.225802,0.830024,0.715596,0.659288,0.600228,0.532512,0.463261,
     & 0.399019,0.343348,0.298259,0.266336,0.251903,0.260341,0.295357,
     & 0.355043,0.428707,0.497582,0.543699,0.572947,0.662090,1.047671
     & /
C ---------------------------------------------------------------------
C     [ 22]    130.0  MeV
      data (hydrang( 22,j),j=1,ihamax)/
     & 1.287767,0.844855,0.724409,0.667950,0.606990,0.536035,0.463673,
     & 0.397062,0.339473,0.292262,0.257610,0.240101,0.246042,0.280329,
     & 0.341818,0.419519,0.493084,0.542351,0.572761,0.668420,1.095658
     & /
C ---------------------------------------------------------------------
C     [ 23]    140.0  MeV
      data (hydrang( 23,j),j=1,ihamax)/
     & 1.343568,0.854883,0.730372,0.674859,0.612436,0.538466,0.463320,
     & 0.394916,0.336214,0.287781,0.251212,0.231178,0.234821,0.268227,
     & 0.331124,0.412384,0.490223,0.542258,0.573018,0.672726,1.135488
     & /
C ---------------------------------------------------------------------
C     [ 24]    150.0  MeV
      data (hydrang( 24,j),j=1,ihamax)/
     & 1.394228,0.861377,0.734838,0.681220,0.617420,0.540174,0.462064,
     & 0.392013,0.332742,0.283938,0.246410,0.224689,0.226568,0.259246,
     & 0.323352,0.407748,0.489352,0.543520,0.573337,0.673619,1.163335
     & /
C ---------------------------------------------------------------------
C     [ 25]    160.0  MeV
      data (hydrang( 25,j),j=1,ihamax)/
     & 1.440719,0.865523,0.738927,0.687948,0.622553,0.541406,0.459774,
     & 0.387902,0.328397,0.280010,0.242566,0.220195,0.221096,0.253430,
     & 0.318708,0.405880,0.490692,0.546208,0.573489,0.670215,1.176622
     & /
C ---------------------------------------------------------------------
C     [ 26]    170.0  MeV
      data (hydrang( 26,j),j=1,ihamax)/
     & 1.483296,0.867906,0.742648,0.694665,0.627436,0.541962,0.456484,
     & 0.382758,0.323363,0.276084,0.239621,0.217490,0.218092,0.250406,
     & 0.316796,0.406389,0.493905,0.550145,0.573659,0.663375,1.177268
     & /
C ---------------------------------------------------------------------
C     [ 27]    180.0  MeV
      data (hydrang( 27,j),j=1,ihamax)/
     & 1.522068,0.868955,0.745729,0.700678,0.631424,0.541531,0.452263,
     & 0.376907,0.318034,0.272450,0.237655,0.216423,0.217207,0.249690,
     & 0.317063,0.408718,0.498522,0.555099,0.574137,0.654392,1.168348
     & /
C ---------------------------------------------------------------------
C     [ 28]    190.0  MeV
      data (hydrang( 28,j),j=1,ihamax)/
     & 1.557152,0.869121,0.747908,0.705292,0.633868,0.539799,0.447178,
     & 0.370669,0.312794,0.269382,0.236733,0.216829,0.218082,0.250795,
     & 0.318957,0.412313,0.504071,0.560841,0.575226,0.644585,1.152949
     & /
C ---------------------------------------------------------------------
C     [ 29]    200.0  MeV
      data (hydrang( 29,j),j=1,ihamax)/
     & 1.588621,0.868821,0.748905,0.707804,0.634120,0.536459,0.441304,
     & 0.364382,0.308045,0.267179,0.236946,0.218567,0.220382,0.253249,
     & 0.321934,0.416621,0.510077,0.567130,0.577202,0.635231,1.134114
     & /
C ---------------------------------------------------------------------
C     [ 30]    220.0  MeV
      data (hydrang( 30,j),j=1,ihamax)/
     & 1.641054,0.867937,0.747025,0.705296,0.626898,0.524497,0.427518,
     & 0.352472,0.300888,0.265890,0.240729,0.225376,0.228082,0.260705,
     & 0.329501,0.425777,0.522100,0.580595,0.584428,0.621307,1.094756
     & /
C ---------------------------------------------------------------------
C     [ 31]    240.0  MeV
      data (hydrang( 31,j),j=1,ihamax)/
     & 1.679389,0.866367,0.741503,0.695240,0.611758,0.507145,0.411784,
     & 0.341467,0.296349,0.267947,0.248014,0.235573,0.238833,0.270522,
     & 0.338315,0.434912,0.533452,0.594420,0.595157,0.613850,1.056111
     & /
C ---------------------------------------------------------------------
C     [ 32]    260.0  MeV
      data (hydrang( 32,j),j=1,ihamax)/
     & 1.703416,0.863458,0.734466,0.681263,0.592141,0.486702,0.395047,
     & 0.331204,0.293565,0.272164,0.257528,0.247906,0.251451,0.281641,
     & 0.347519,0.443378,0.543510,0.607682,0.608110,0.612217,1.020226
     & /
C ---------------------------------------------------------------------
C     [ 33]    280.0  MeV
      data (hydrang( 33,j),j=1,ihamax)/
     & 1.712959,0.858573,0.728031,0.666979,0.571485,0.465475,0.378262,
     & 0.321521,0.291665,0.277330,0.267969,0.261097,0.264736,0.293001,
     & 0.356263,0.450534,0.551657,0.619454,0.622017,0.615790,0.989181
     & /
C ---------------------------------------------------------------------
C     [ 34]    300.0  MeV
      data (hydrang( 34,j),j=1,ihamax)/
     & 1.707881,0.851087,0.724340,0.656023,0.553232,0.445753,0.362359,
     & 0.312235,0.289765,0.282235,0.278039,0.273871,0.277480,0.303522,
     & 0.363677,0.455731,0.557279,0.628830,0.635626,0.623959,0.965081
     & /
C ---------------------------------------------------------------------
C     [ 35]    320.0  MeV
      data (hydrang( 35,j),j=1,ihamax)/
     & 1.688703,0.840574,0.724642,0.650765,0.539752,0.429163,0.347997,
     & 0.303178,0.287154,0.285901,0.286685,0.285208,0.288773,0.312482,
     & 0.369296,0.458754,0.560246,0.635445,0.648147,0.636066,0.949372
     & /
C ---------------------------------------------------------------------
C     [ 36]    340.0  MeV
      data (hydrang( 36,j),j=1,ihamax)/
     & 1.658974,0.827568,0.726781,0.648637,0.529153,0.414594,0.334615,
     & 0.294092,0.283692,0.288188,0.293742,0.295015,0.298763,0.320425,
     & 0.374138,0.461079,0.562406,0.641244,0.660819,0.651440,0.941223
     & /
C ---------------------------------------------------------------------
C     [ 37]    360.0  MeV
      data (hydrang( 37,j),j=1,ihamax)/
     & 1.622960,0.812801,0.727721,0.645832,0.518486,0.400266,0.321361,
     & 0.284709,0.279395,0.289178,0.299288,0.303451,0.307885,0.328243,
     & 0.379614,0.464616,0.566100,0.648728,0.675347,0.669360,0.939229
     & /
C ---------------------------------------------------------------------
C     [ 38]    380.0  MeV
      data (hydrang( 38,j),j=1,ihamax)/
     & 1.584894,0.797019,0.724432,0.638530,0.504788,0.384392,0.307386,
     & 0.274775,0.274287,0.288961,0.303400,0.310676,0.316575,0.336819,
     & 0.387129,0.471268,0.573661,0.660399,0.693454,0.689126,0.941916
     & /
C ---------------------------------------------------------------------
C     [ 39]    400.0  MeV
      data (hydrang( 39,j),j=1,ihamax)/
     & 1.549050,0.780952,0.713877,0.622924,0.485109,0.365189,0.291834,
     & 0.264020,0.268387,0.287623,0.306162,0.316853,0.325267,0.347038,
     & 0.398089,0.482939,0.587430,0.678764,0.716851,0.710017,0.947875
     & /
C ---------------------------------------------------------------------
C     [ 40]    420.0  MeV
      data (hydrang( 40,j),j=1,ihamax)/
     & 1.519935,0.765193,0.694969,0.597713,0.458550,0.342266,0.274694,
     & 0.252609,0.261849,0.285191,0.307438,0.321722,0.333633,0.358541,
     & 0.412138,0.499365,0.607412,0.704302,0.746457,0.732919,0.957860
     & /
C ---------------------------------------------------------------------
C     [ 41]    440.0  MeV
      data (hydrang( 41,j),j=1,ihamax)/
     & 1.503023,0.749671,0.674349,0.571640,0.432468,0.320829,0.259333,
     & 0.242440,0.255396,0.281487,0.306285,0.323384,0.338315,0.366003,
     & 0.421871,0.511585,0.624280,0.729386,0.779970,0.765052,0.981401
     & /
C ---------------------------------------------------------------------
C     [ 42]    460.0  MeV
      data (hydrang( 42,j),j=1,ihamax)/
     & 1.503409,0.734069,0.660221,0.555492,0.415869,0.307170,0.249735,
     & 0.235704,0.249838,0.276295,0.301620,0.319650,0.335403,0.363178,
     & 0.418583,0.509049,0.627009,0.744909,0.814450,0.814698,1.029647
     & /
C ---------------------------------------------------------------------
C     [ 43]    480.0  MeV
      data (hydrang( 43,j),j=1,ihamax)/
     & 1.510740,0.716149,0.652177,0.549054,0.408257,0.300450,0.244932,
     & 0.231612,0.244769,0.269618,0.293774,0.311124,0.325844,0.351533,
     & 0.404442,0.494696,0.619070,0.753976,0.850702,0.878372,1.100272
     & /
C ---------------------------------------------------------------------
C     [ 44]    500.0  MeV
      data (hydrang( 44,j),j=1,ihamax)/
     & 1.499344,0.691730,0.641182,0.541134,0.399653,0.292698,0.238981,
     & 0.226371,0.238565,0.261694,0.284505,0.301204,0.315424,0.340204,
     & 0.392450,0.484932,0.618446,0.771960,0.893294,0.940809,1.177564
     & /
C ---------------------------------------------------------------------
C     [ 45]    520.0  MeV
      data (hydrang( 45,j),j=1,ihamax)/
     & 1.451080,0.658590,0.620849,0.523476,0.382668,0.278025,0.227463,
     & 0.217105,0.229948,0.252691,0.275234,0.292665,0.308778,0.336325,
     & 0.392476,0.491796,0.637750,0.808960,0.944484,0.990916,1.249173
     & /
C ---------------------------------------------------------------------
C     [ 46]    540.0  MeV
      data (hydrang( 46,j),j=1,ihamax)/
     & 1.380849,0.622599,0.596840,0.501504,0.361947,0.260081,0.212853,
     & 0.205061,0.219220,0.242512,0.265872,0.285360,0.305131,0.337651,
     & 0.399927,0.507505,0.665627,0.851931,0.996626,1.036235,1.318739
     & /
C ---------------------------------------------------------------------
C     [ 47]    560.0  MeV
      data (hydrang( 47,j),j=1,ihamax)/
     & 1.311435,0.591547,0.577630,0.483817,0.344916,0.244729,0.199241,
     & 0.192465,0.207067,0.230973,0.255909,0.278359,0.302319,0.339626,
     & 0.406741,0.519631,0.685240,0.882670,1.039984,1.088592,1.393777
     & /
C ---------------------------------------------------------------------
C     [ 48]    580.0  MeV
      data (hydrang( 48,j),j=1,ihamax)/
     & 1.255852,0.569105,0.566078,0.472849,0.333505,0.233599,0.188062,
     & 0.180467,0.194223,0.218316,0.245025,0.270644,0.298425,0.339255,
     & 0.408864,0.523390,0.691782,0.897417,1.073127,1.149920,1.479188
     & /
C ---------------------------------------------------------------------
C     [ 49]    600.0  MeV
      data (hydrang( 49,j),j=1,ihamax)/
     & 1.217304,0.554823,0.559417,0.464839,0.324123,0.224081,0.178089,
     & 0.169139,0.181471,0.205191,0.233066,0.261092,0.291766,0.335113,
     & 0.406259,0.521651,0.692490,0.906871,1.102973,1.212330,1.577263
     & /
C ---------------------------------------------------------------------
C     [ 50]    620.0  MeV
      data (hydrang( 50,j),j=1,ihamax)/
     & 1.196107,0.547102,0.554559,0.456305,0.313596,0.213808,0.168097,
     & 0.158375,0.169381,0.192189,0.220143,0.249242,0.281608,0.326748,
     & 0.399520,0.517114,0.693196,0.919332,1.135153,1.270960,1.688315
     & /
C ---------------------------------------------------------------------
C     [ 51]    640.0  MeV
      data (hydrang( 51,j),j=1,ihamax)/
     & 1.182502,0.540334,0.548000,0.445822,0.301235,0.202006,0.157220,
     & 0.147460,0.157586,0.179533,0.207374,0.237286,0.271004,0.317418,
     & 0.391187,0.510580,0.692127,0.931096,1.168889,1.335043,1.805463
     & /
C ---------------------------------------------------------------------
C     [ 52]    660.0  MeV
      data (hydrang( 52,j),j=1,ihamax)/
     & 1.164840,0.528262,0.536435,0.432658,0.287104,0.188434,0.144813,
     & 0.135639,0.145577,0.167396,0.196105,0.227972,0.263775,0.311108,
     & 0.384123,0.502312,0.685685,0.935622,1.201712,1.416459,1.920325
     & /
C ---------------------------------------------------------------------
C     [ 53]    680.0  MeV
      data (hydrang( 53,j),j=1,ihamax)/
     & 1.145835,0.512510,0.523027,0.419589,0.273530,0.175125,0.132476,
     & 0.124050,0.134251,0.156604,0.186983,0.221536,0.259719,0.307411,
     & 0.377896,0.491665,0.672642,0.930880,1.230067,1.508278,2.032605
     & /
C ---------------------------------------------------------------------
C     [ 54]    700.0  MeV
      data (hydrang( 54,j),j=1,ihamax)/
     & 1.142606,0.502646,0.517450,0.412925,0.265122,0.166400,0.124045,
     & 0.115699,0.125876,0.148586,0.179896,0.215642,0.254549,0.301478,
     & 0.368741,0.477074,0.654140,0.919368,1.249362,1.584846,2.150206
     & /
C ---------------------------------------------------------------------
C     [ 55]    720.0  MeV
      data (hydrang( 55,j),j=1,ihamax)/
     & 1.163666,0.505735,0.525523,0.416013,0.264557,0.165140,0.122205,
     & 0.112661,0.121876,0.143993,0.174245,0.208039,0.244685,0.289507,
     & 0.353950,0.457881,0.632031,0.904404,1.258178,1.629625,2.276119
     & /
C ---------------------------------------------------------------------
C     [ 56]    740.0  MeV
      data (hydrang( 56,j),j=1,ihamax)/
     & 1.180581,0.517503,0.536503,0.419698,0.266375,0.168045,0.124616,
     & 0.112908,0.120010,0.140207,0.167548,0.197207,0.229997,0.272602,
     & 0.335621,0.437259,0.610635,0.891857,1.268083,1.665963,2.392462
     & /
C ---------------------------------------------------------------------
C     [ 57]    760.0  MeV
      data (hydrang( 57,j),j=1,ihamax)/
     & 1.157803,0.530896,0.536073,0.412270,0.263355,0.170337,0.127679,
     & 0.113409,0.117185,0.133906,0.157001,0.181968,0.211367,0.253194,
     & 0.317080,0.419291,0.594761,0.887921,1.293063,1.725676,2.477320
     & /
C ---------------------------------------------------------------------
C     [ 58]    780.0  MeV
      data (hydrang( 58,j),j=1,ihamax)/
     & 1.107315,0.540021,0.522684,0.394542,0.254306,0.168768,0.128090,
     & 0.111935,0.112125,0.124576,0.143345,0.164978,0.193217,0.236311,
     & 0.302673,0.407110,0.586003,0.891256,1.328084,1.806346,2.535667
     & /
C ---------------------------------------------------------------------
C     [ 59]    800.0  MeV
      data (hydrang( 59,j),j=1,ihamax)/
     & 1.088855,0.540253,0.507654,0.379935,0.244112,0.161620,0.122808,
     & 0.107012,0.105313,0.114459,0.130810,0.152672,0.183451,0.229509,
     & 0.297712,0.402862,0.582740,0.893009,1.349117,1.871340,2.599496
     & /
C ---------------------------------------------------------------------
C     [ 60]    820.0  MeV
      data (hydrang( 60,j),j=1,ihamax)/
     & 1.143453,0.528774,0.500166,0.378360,0.236793,0.148763,0.110838,
     & 0.098500,0.097835,0.105831,0.122818,0.149387,0.186688,0.236669,
     & 0.304304,0.406353,0.582375,0.886850,1.341023,1.894998,2.690551
     & /
C ---------------------------------------------------------------------
C     [ 61]    840.0  MeV
      data (hydrang( 61,j),j=1,ihamax)/
     & 1.229384,0.510025,0.498818,0.383632,0.231920,0.136141,0.099285,
     & 0.091382,0.092762,0.100529,0.118829,0.150296,0.193661,0.246373,
     & 0.311480,0.408125,0.578944,0.877852,1.327607,1.901641,2.785135
     & /
C ---------------------------------------------------------------------
C     [ 62]    860.0  MeV
      data (hydrang( 62,j),j=1,ihamax)/
     & 1.285610,0.489962,0.499357,0.385638,0.227920,0.130942,0.096943,
     & 0.091662,0.093491,0.100180,0.117358,0.148560,0.192166,0.244061,
     & 0.305687,0.397032,0.566065,0.873738,1.341175,1.926863,2.849752
     & /
C ---------------------------------------------------------------------
C     [ 63]    880.0  MeV
      data (hydrang( 63,j),j=1,ihamax)/
     & 1.284413,0.468184,0.493283,0.376839,0.222684,0.134310,0.105144,
     & 0.099657,0.099583,0.104136,0.117791,0.143693,0.182119,0.230524,
     & 0.288880,0.375883,0.546619,0.877152,1.385555,1.977697,2.886191
     & /
C ---------------------------------------------------------------------
C     [ 64]    900.0  MeV
      data (hydrang( 64,j),j=1,ihamax)/
     & 1.231119,0.437829,0.467762,0.352219,0.213553,0.141340,0.117783,
     & 0.110032,0.106822,0.109541,0.120450,0.141601,0.175579,0.221922,
     & 0.278540,0.361480,0.532784,0.885669,1.436061,2.032597,2.931518
     & /
C ---------------------------------------------------------------------
C     [ 65]    920.0  MeV
      data (hydrang( 65,j),j=1,ihamax)/
     & 1.138570,0.394578,0.415447,0.310389,0.199229,0.147799,0.129649,
     & 0.118573,0.111996,0.114045,0.125258,0.146550,0.181618,0.230336,
     & 0.287444,0.365825,0.532753,0.895516,1.471168,2.074654,3.013222
     & /
C ---------------------------------------------------------------------
C     [ 66]    940.0  MeV
      data (hydrang( 66,j),j=1,ihamax)/
     & 1.043278,0.345147,0.351507,0.263970,0.184025,0.153289,0.140505,
     & 0.126657,0.116716,0.117846,0.130498,0.155150,0.195208,0.248784,
     & 0.306771,0.379259,0.537149,0.898392,1.486915,2.110567,3.114464
     & /
C ---------------------------------------------------------------------
C     [ 67]    960.0  MeV
      data (hydrang( 67,j),j=1,ihamax)/
     & 0.986028,0.298764,0.296241,0.228616,0.173467,0.158371,0.151364,
     & 0.136985,0.123699,0.121733,0.134093,0.162254,0.208054,0.265879,
     & 0.322745,0.387188,0.532717,0.885390,1.483962,2.152224,3.207945
     & /
C ---------------------------------------------------------------------
C     [ 68]    980.0  MeV
      data (hydrang( 68,j),j=1,ihamax)/
     & 0.971908,0.260567,0.258981,0.209740,0.169539,0.164097,0.163133,
     & 0.150031,0.132822,0.125093,0.134927,0.165997,0.217340,0.277907,
     & 0.331050,0.385474,0.517417,0.859151,1.468001,2.195769,3.279608
     & /
C ---------------------------------------------------------------------
C     [ 69]   1000.0  MeV
      data (hydrang( 69,j),j=1,ihamax)/
     & 0.969884,0.231350,0.237988,0.202458,0.170656,0.172023,0.176686,
     & 0.164142,0.141216,0.125987,0.132894,0.167833,0.225761,0.288860,
     & 0.336870,0.380456,0.500439,0.833878,1.449768,2.221603,3.328851
     & /
C ---------------------------------------------------------------------
C     [ 70]   3000.0  MeV
      data (hydrang( 70,j),j=1,ihamax)/
     & 0.969884,0.231350,0.237988,0.202458,0.170656,0.172023,0.176686,
     & 0.164142,0.141216,0.125987,0.132894,0.167833,0.225761,0.288860,
     & 0.336870,0.380456,0.500439,0.833878,1.449768,2.221603,3.328851
     & /
C ---------------------------------------------------------------------
      end
C =====================================================================
C =====================================================================
C ---------------------------------------------------------------------
      function carbonsigma(iflag,emev)
C ---------------------------------------------------------------------
C     Purpose is to retrun the total, elastic, inelastic cross sections
C     of carbon for neutrons from 80 to 3000 MeV.
C ---------------------------------------------------------------------
C     This function is called from s:sigrc when the iscinful = 1 and
C     high-energy models are used for nuclear reaction simulations.
C ---------------------------------------------------------------------
C     (IN) emev:  incidnet neutron energy in lab system (MeV)
C          iflag: flag for the mode
C                 =1 total cross section
C                 =2 elastic cross section
C                 =3 inelastic cross section
C     (OUT)carbonsigma: cross section of neutron-carbon collsiont (barn)
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
      parameter(icemax=314)
      common/carcmm/csigenerg(icemax), csigtotal(icemax),
     & csiginela(icemax)
C ---------------------------------------------------------------------
      real*8 carbonsigma
      real*8 emev
      dimension sig(icemax), e(icemax)
C ---------------------------------------------------------------------
      nterp = 1
      nn = icemax
      enn=sngl(emev)

      e(:) = csigenerg(:)

      if( iflag == 1 ) then
        sig(:) = csigtotal(:)
      else if( iflag == 2 ) then
        sig(:) = csigtotal(:) - csiginela(:)
      else if( iflag == 3 ) then
        sig(:) = csiginela(:)
      else
        call parastop(122)
      end if

      sigma=EXTERP(e,sig,enn,nn,nterp)
      carbonsigma=dble(sigma) / 1.0e+03  ! (barn)

      return
      end
C =====================================================================
C ---------------------------------------------------------------------
      block data cardat
C ---------------------------------------------------------------------
C     Data for total and inelastic cross sections from 80-3000 MeV.
C     Numerical values were taken from JENDL/HE-2007.
C ---------------------------------------------------------------------
C     2021.06.01
C     The data from 150-200 MeV in JENDL/HE-2007 were omitted to
C     avoid strange structure in inelastic cross-section curve.
C ---------------------------------------------------------------------
C ---------------------------------------------------------------------
      parameter(icemax=314)
      common/carcmm/csigenerg(icemax), csigtotal(icemax),
     & csiginela(icemax)
C ---------------------------------------------------------------------
C    Neutron energy (MeV)
      data (csigenerg(i),i=1,icemax)/
     & 80.000,81.129,82.000,82.273,83.157,83.434,84.000,
     & 84.611,85.000,86.000,86.199,87.213,87.415,88.000,
     & 88.649,89.242,89.899,90.000,91.270,92.557,93.863,
     & 95.000,96.340,97.699,99.078,100.000,101.411,102.841,
     & 104.292,105.764,107.256,108.769,110.000,111.552,113.126,
     & 114.722,116.340,117.982,119.646,120.000,121.693,123.410,
     & 125.151,126.917,128.707,130.000,131.834,133.694,135.580,
     & 137.493,139.433,140.000,141.975,143.978,146.009,148.069,
     & 150.000,200.000,202.822,205.683,208.585,211.528,214.512,
     & 217.538,220.000,223.104,226.251,229.443,232.680,235.963,
     & 239.292,240.000,243.386,246.820,250.000,250.302,253.527,
     & 253.833,257.104,257.414,260.000,263.668,267.388,271.160,
     & 274.986,278.865,280.000,283.950,287.956,292.019,296.139,
     & 300.000,304.232,308.525,312.877,317.291,320.000,324.515,
     & 329.093,333.736,338.444,340.000,344.797,349.661,350.000,
     & 354.594,354.938,359.597,359.945,360.000,365.079,370.230,
     & 375.453,380.000,385.361,390.798,396.311,400.000,405.643,
     & 411.366,417.170,420.000,423.055,425.925,429.024,431.934,
     & 435.077,438.028,440.000,441.215,446.208,447.439,450.000,
     & 452.503,456.349,458.887,460.000,462.787,466.490,469.316,
     & 473.071,475.937,479.745,480.000,482.652,486.772,489.461,
     & 493.639,496.367,500.000,507.054,514.208,520.000,521.462,
     & 527.336,528.819,534.776,536.280,540.000,543.846,547.618,
     & 550.000,555.344,557.760,560.000,565.629,567.901,573.608,
     & 575.913,580.000,581.701,588.183,589.908,596.481,598.230,
     & 600.000,608.465,617.049,625.755,634.583,643.536,650.000,
     & 652.615,659.170,661.822,668.470,671.159,677.901,680.628,
     & 687.465,690.231,697.164,699.968,700.000,709.876,719.891,
     & 730.047,740.347,750.000,750.792,760.581,761.384,771.312,
     & 772.126,782.193,783.019,793.229,794.066,800.000,811.287,
     & 822.732,834.340,846.111,850.000,858.048,861.992,870.153,
     & 874.153,882.430,886.486,894.879,898.993,900.000,912.697,
     & 925.574,938.632,950.000,951.875,963.403,965.304,976.995,
     & 978.922,990.778,992.733,1000.000,1014.108,1028.415,1042.925,
     & 1057.638,1072.560,1087.692,1103.037,1118.599,1134.380,1150.384,
     & 1166.614,1183.073,1199.764,1216.691,1233.856,1251.263,1268.916,
     & 1286.819,1304.973,1323.384,1342.055,1360.989,1380.190,1399.662,
     & 1419.408,1439.434,1459.741,1480.336,1500.000,1521.162,1542.623,
     & 1564.387,1586.457,1608.840,1631.537,1654.555,1677.898,1701.570,
     & 1725.576,1749.921,1774.609,1799.646,1825.036,1850.784,1876.895,
     & 1903.375,1930.228,1957.460,1985.076,2000.000,2028.216,2056.831,
     & 2085.849,2115.277,2145.119,2175.383,2206.074,2237.198,2268.760,
     & 2300.769,2333.228,2366.146,2399.528,2433.381,2467.712,2502.527,
     & 2537.833,2573.637,2609.946,2646.768,2684.109,2721.977,2760.379,
     & 2799.323,2838.817,2878.867,2919.483,2960.672,3000.000
     & /
C ---------------------------------------------------------------------
C    Total cross sections of carbon nucleus (mb)
      data (csigtotal(i),i=1,icemax)/
     & 618.310,611.538,606.310,603.832,595.817,593.304,588.170,
     & 583.941,581.250,574.330,572.923,565.758,564.331,560.200,
     & 555.443,551.093,546.269,545.530,539.095,532.569,525.951,
     & 520.190,513.397,506.509,499.524,494.850,488.053,481.159,
     & 474.169,467.080,459.891,452.600,446.670,441.386,436.027,
     & 430.593,425.082,419.493,413.825,412.620,407.475,402.257,
     & 396.966,391.600,386.159,382.230,378.395,374.506,370.562,
     & 366.562,362.506,361.320,357.532,353.690,349.794,345.843,
     & 342.140,293.760,292.578,291.379,290.163,288.930,287.679,
     & 286.411,285.380,284.551,283.711,282.859,281.994,281.118,
     & 280.229,280.040,279.869,279.696,279.535,279.520,279.357,
     & 279.341,279.176,279.161,279.030,278.957,278.882,278.807,
     & 278.730,278.653,278.630,278.095,277.552,277.001,276.443,
     & 275.920,276.760,277.612,278.476,279.352,279.890,280.258,
     & 280.631,281.009,281.393,281.520,282.530,283.554,283.625,
     & 284.592,284.664,285.645,285.719,285.730,286.578,287.438,
     & 288.311,289.070,288.126,287.170,286.199,285.550,288.908,
     & 292.313,295.766,297.450,297.532,297.610,297.694,297.772,
     & 297.857,297.937,297.990,298.282,299.483,299.779,300.395,
     & 300.997,301.922,302.532,302.800,302.521,302.151,301.868,
     & 301.493,301.206,300.825,300.800,302.098,304.115,305.431,
     & 307.476,308.811,310.590,315.217,319.910,323.710,323.720,
     & 323.761,323.772,323.813,323.824,323.850,323.913,323.976,
     & 324.015,324.103,324.143,324.180,325.041,325.389,326.262,
     & 326.615,327.240,328.518,333.389,334.686,339.625,340.940,
     & 342.270,342.626,342.986,343.352,343.723,344.099,344.371,
     & 344.459,344.680,344.770,344.994,345.085,345.312,345.404,
     & 345.635,345.728,345.962,346.057,346.058,346.323,346.592,
     & 346.864,347.141,347.400,347.417,347.625,347.642,347.853,
     & 347.870,348.084,348.102,348.319,348.337,348.463,348.652,
     & 348.844,349.039,349.237,349.302,349.408,349.460,349.568,
     & 349.620,349.729,349.783,349.894,349.948,349.961,350.092,
     & 350.225,350.360,350.478,350.493,350.586,350.602,350.696,
     & 350.712,350.807,350.823,350.882,350.919,350.957,350.995,
     & 351.034,351.073,351.113,351.153,351.194,351.236,351.278,
     & 351.321,351.364,351.408,351.452,351.498,351.543,351.590,
     & 351.637,351.685,351.733,351.782,351.832,351.883,351.934,
     & 351.986,352.039,352.092,352.146,352.198,352.202,352.207,
     & 352.212,352.216,352.221,352.226,352.231,352.236,352.241,
     & 352.246,352.251,352.256,352.262,352.267,352.272,352.278,
     & 352.284,352.289,352.295,352.301,352.304,352.304,352.304,
     & 352.305,352.305,352.305,352.305,352.306,352.306,352.306,
     & 352.306,352.307,352.307,352.307,352.307,352.308,352.308,
     & 352.308,352.309,352.309,352.309,352.309,352.310,352.310,
     & 352.310,352.311,352.311,352.311,352.312,352.312
     & /
C ---------------------------------------------------------------------
C    Inelastic cross sections of carbon nucleus (mb)
      data (csiginela(i),i=1,icemax)/
     & 250.799,251.252,251.601,250.872,248.513,247.774,246.263,
     & 245.946,245.744,245.144,244.995,244.240,244.089,243.653,
     & 242.995,242.393,241.725,241.623,240.733,239.831,238.916,
     & 238.120,237.343,236.555,235.756,235.221,234.617,234.005,
     & 233.384,232.755,232.116,231.469,230.942,230.573,230.199,
     & 229.820,229.436,229.046,228.650,228.566,228.351,228.134,
     & 227.913,227.690,227.463,227.299,227.250,227.200,227.149,
     & 227.098,227.046,227.031,227.784,228.547,229.321,230.105,
     & 230.841,191.987,191.861,191.733,191.603,191.471,191.337,
     & 191.202,191.092,191.051,191.010,190.968,190.925,190.882,
     & 190.838,190.829,191.298,191.774,192.215,192.247,192.586,
     & 192.619,192.963,192.996,193.268,193.622,193.982,194.346,
     & 194.716,195.090,195.200,194.972,194.742,194.508,194.270,
     & 194.048,195.103,196.173,197.259,198.359,199.035,199.546,
     & 200.065,200.591,201.125,201.301,202.413,203.540,203.619,
     & 204.658,204.735,205.789,205.868,205.880,206.775,207.681,
     & 208.601,209.402,208.472,207.528,206.572,205.932,209.266,
     & 212.647,216.075,217.748,217.817,217.883,217.953,218.020,
     & 218.091,218.159,218.204,218.491,219.670,219.961,220.566,
     & 221.147,222.040,222.628,222.887,222.584,222.183,221.876,
     & 221.469,221.158,220.745,220.717,221.993,223.975,225.268,
     & 227.278,228.590,230.338,234.899,239.525,243.271,243.267,
     & 243.253,243.250,243.236,243.232,243.223,243.251,243.278,
     & 243.295,243.336,243.354,243.372,244.183,244.511,245.334,
     & 245.666,246.256,247.519,252.333,253.614,258.496,259.795,
     & 261.109,261.405,261.705,262.009,262.317,262.630,262.855,
     & 262.925,263.100,263.170,263.347,263.419,263.598,263.671,
     & 263.853,263.927,264.112,264.186,264.187,264.404,264.625,
     & 264.848,265.075,265.287,265.300,265.461,265.474,265.637,
     & 265.650,265.816,265.829,265.997,266.010,266.108,266.262,
     & 266.418,266.577,266.738,266.791,266.872,266.911,266.993,
     & 267.034,267.117,267.158,267.242,267.283,267.294,267.400,
     & 267.508,267.617,267.713,267.724,267.795,267.806,267.878,
     & 267.890,267.962,267.974,268.019,268.048,268.077,268.107,
     & 268.138,268.168,268.200,268.231,268.263,268.296,268.329,
     & 268.362,268.396,268.431,268.466,268.501,268.537,268.574,
     & 268.611,268.648,268.686,268.725,268.764,268.803,268.844,
     & 268.884,268.926,268.968,269.010,269.051,269.054,269.058,
     & 269.061,269.065,269.069,269.073,269.076,269.080,269.084,
     & 269.088,269.092,269.097,269.101,269.105,269.109,269.114,
     & 269.118,269.123,269.127,269.132,269.134,269.134,269.135,
     & 269.135,269.135,269.135,269.135,269.136,269.136,269.136,
     & 269.136,269.136,269.137,269.137,269.137,269.137,269.137,
     & 269.138,269.138,269.138,269.138,269.139,269.139,269.139,
     & 269.139,269.140,269.140,269.140,269.140,269.141
     & /
C ---------------------------------------------------------------------
      end
C =====================================================================
