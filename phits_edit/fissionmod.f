************************************************************************
*                                                                      *
      module fission_mod
*                                                                      *
*                                                                      *
*        Last Revised:     2020 05 14  by T.Ogawa                      *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              fission physics module                                  *
*                                                                      *
*                                                                      *
*                                                                      *
*       subroutines and functions                                      *
*                                                                      *
*        subroutine fis(a,z,u)                                  *
*          Old fission algorithm. Determine scission product           *
*          configuration                                               *
*                                                                      *
*        function fprob (z,a,e)
*          Old fission probability calculation algorithm               *
*                                                                      *
*                                                                      *
*        subroutine fis2(a,z,u)                                  *
*          New fission algorithm. Determine scission product           *
*          configuration                                               *
*                                                                      *
*        function fprob2 (z,a,e)
*          New fission probability calculation algorithm               *
*                                                                      *
*                                                                      *
************************************************************************

      use NGSDATAMOD, only : energm
      implicit doubleprecision(a-h, o-z)
      include 'err.inc'
      include 'param-physcnst.inc'

      common /options/alev, rcal, ifis
!$OMP THREADPRIVATE(/options/)
      common /mom/ erec,bet0(3),ekin,bet1(3)
!$OMP THREADPRIVATE(/mom/)

      double precision, private :: ef
      logical fisinh
      integer ifiss

      double precision, dimension(2) :: zfis,afis,ufis,er
      double precision, dimension(2,3) :: betf
!$OMP THREADPRIVATE(fisinh,ef,zfis,afis,ufis,er,betf)

      contains

************************************************************************
*                                                                      *
      subroutine fis(a,z,u)
C/////////////////////////////////////////////////////////////////////
C****************This routine is originally in LAHET code**************
C  FIS
c      Pick post fission parameters such as mass, charge, kinetic energy
c      and excitation energy
C=====================================================================
C <variables>
C     a    :   mass of nucleus before fission                     (IN)
C     z    :   charge  of nucleus before fission                  (IN)
C     u    :   excitation energy of nucleus before fission        (IN)
C   ifiss  :   fission model version. 1:LAHET, 2:Iwamoto          (IN)
C  zfis    :  charge of the fission fragment                     (OUT)
C  afis    :  mass of the fission fragment                       (OUT)
C  ufis    :  excitation energy of the fission fragment          (OUT)
C   er     :  recoil energy of the fission fragment              (OUT)
C  betf    :  recoil direction of the fission fragment           (OUT)
C/////////////////////////////////////////////////////////////////////

      implicit doubleprecision(a-h,o-z)

      parameter (pi=3.1415926535898d0)

      dimension evodba(4)
c     following data for assymetric vs symetric picking
      data evodba/18.8d0,18.1d0,18.1d0,18.5d0/
      data sigmaa, aamean /6.5d0,140.d0/
c     asym gauss width,1.0/level density parameter and high mass mean.
      data afact, zfact, enmass /8.071323d0,0.782354d0,939.56563d0/
c     mass diff from 1 rstms(0) for neutron,diff p & n masses and n mass

      if (.not.fisinh) then
c***********************************************************************
c
c  logic error fissed called with fisinh unset.
c
c***********************************************************************
       write (*,60) a,z,u,erec
       return
      endif
c***********************************************************************
c
c       pick the masses
c
c***********************************************************************
      jz=nint(z)
      ja=nint(a)
      nck2=0
   10 continue
      nck2=nck2+1
      if (nck2.gt.10) then
c***********************************************************************
c
c     fission failure
c
c***********************************************************************
        write (*,50) ja,jz,u,erec,zfis(1),afis(1),zfis(2),afis(2)
        fisinh=.false.
        return
      endif
      temp=z*z/a
      if (temp.le.34.d0) then
        if (jz.le.88) go to 20
      elseif (u.le.62.d0) then
c***********************************************************************
c
c   high z fission mass distribution. competition for sym vs assym
c  simple symmetric to assymetric data fit
c
c***********************************************************************
        arg=-0.36d0*u
        arg=4.87d+03*exp(arg)
        proba=arg/(1.d0+arg)
        if (rn(0).le.proba) then
c***********************************************************************
c
c      assymetric fission
c
c***********************************************************************
          a1=gaussn(aamean,sigmaa)
          go to 30
        endif
      endif
c***********************************************************************
c
c     find assymetric barrier for width computation
c     assymetric barrier from seaborg and vandenbosch
c     phys.rev. 88,507 (1952) phys.rev. 110,507 (1958)
c     find if ee eo oe or oo nucleus
c     na=1 odd-odd,2 even-odd,3 odd-even,4 even-even
c
c***********************************************************************
      if(a.gt.240) then        ! furi May 4, 1999
       in=ja-jz
       na=1
       if (jz.eq.2*(jz/2)) na=na+1
       if (in.eq.2*(in/2)) na=na+2
       temp=z*z/a
       ef=evodba(na)-0.36d0*temp
      endif                    ! furi May 4, 1999
   20 continue
C  furi ... May 4, 1999
      if(ifis.ne.0) then
       upr=u-ef
       if (upr.gt.1.d+02) upr=1.d+02
       sigmas=0.425d0*(upr*(1.d0-0.005d0*upr)+9.35d0)

      else

       xx=z*z/a
       bf=efms(z,a)
       if(bf.lt.0)bf=ef
       upr=min(u-bf,450.0d0)

       sigmas = 0.122d0*xx**2.d0 - 7.77 * xx + 134.0d0 + 3.32d-2 * upr

      endif
c***********************************************************************
c
c     sigmas is symmetrin fission mass width.taken from systematic of
c     neuzil & fairhall phys.rev. 129,2705,(1963)
c
c      low z fission is always symmetric
c
c      high z fission is sometimes.loop back here from 1 loop if
c      symmetric fission predicted.
c
c***********************************************************************
C  furi ... May 4, 1999
      amean=0.5d0*a
      a1=gaussn(amean,sigmas)
   30 continue
c***********************************************************************
c
c      1 loop for assymmetrin fission returns to here.
c
c***********************************************************************
      afis(1)=aint(a1)
c***********************************************************************
c
c    check for low final a
c
c***********************************************************************
      if (afis(1).lt.5.d0) afis(1)=5.d0
      if (a-afis(1).lt.5.d0) afis(1)=a-5.d0
      afis(2)=a-afis(1)
c***********************************************************************
c
c        pick the charge
c
c***********************************************************************
      z1=65.5d0*afis(1)/(131.d0+afis(1)**0.666667d0)
      z2=65.5d0*afis(2)/(131.d0+afis(2)**0.666667d0)
      z1=z1+.5d0*(z-z1-z2)
c***********************************************************************
c
c       we use constant charge density with a 2 unit gaussian smearing
c
c***********************************************************************
c   sigma = 0.75;   z1=gaussn(z1,dph)
c                                by furi 18/DEC/1997
c***********************************************************************
      if(ifis.ne.0) then
       sigz=2.0d0
      else
       sigz=0.75d0
      endif
      z1=gaussn(z1,sigz)
      zfis(1)=aint(z1)
      zfis(2)=z-zfis(1)
c***********************************************************************
c
c   check for reasonable z a combinations.
c
c***********************************************************************
      if (zfis(1).ge.afis(1)) go to 10
      if (zfis(2).ge.afis(2)) go to 10
      if (zfis(1).lt.1.d0) go to 10
      if (zfis(2).lt.1.d0) go to 10
c***********************************************************************
c
c      compute binding energy and actual masses of fragments
c
c***********************************************************************
      be0=afact*a-zfact*z-energm(jz,ja,2)
      rm0=enmass*a-zfact*z-be0
      iaf=nint(afis(1))
      izf=nint(zfis(1))
      be1=afact*afis(1)-zfact*zfis(1)-energm(izf,iaf,2)
      rm1=enmass*afis(1)-zfact*zfis(1)-be1
      iaf=nint(afis(2))
      izf=nint(zfis(2))
      be2=afact*afis(2)-zfact*zfis(2)-energm(izf,iaf,2)
      rm2=enmass*afis(2)-zfact*zfis(2)-be2
c***********************************************************************
c
c      pick recoil kinetic energy.use systematic of ...........
c      unik et.al. proc.3rd iaea symp.on phy.& chem. fision,rochester.vo
c
c***********************************************************************
      if(ifis.ne.0)then
       totkm=0.13323d0*z*z/a**0.33333333d0-11.4d0
      else
C... 4 May 1999 by furi
       x=z*z/a**0.33333d0
       if(x.le.900.d0) then
        totkm=0.131d0*x
       else if(x.le.1800.d0) then
        totkm=0.104d0*x+24.3d0
       else
        write(ErrCha,*)'Error in totk in subroutine fiss'
        ErrID = 'L:267/R:fis/F:fissionmod.f' !E00_002_001
        call ErrWrite(ErrID,ErrCha)
       endif
      endif
c***********************************************************************
c
c      use a width of 15% value at half height.
c
c***********************************************************************
      if(ifis.ne.0)then
       sigmak=0.084d0*totkm
C... 4 May 1999 by furi
      else
       if(x.lt.1000)then
        sigmak=86.5d0
       else if (x.lt.1800.d0)then
        sigmak=5.70d-4*(x-1000.d0)**2.d0+86.5d0
       else
        write(ErrCha,*)'Error in sigmak in subroutine fiss'
        ErrID = 'L:286/R:fis/F:fissionmod.f' !E00_003_001
        call ErrWrite(ErrID,ErrCha)
       endif
      sigmak=sqrt(sigmak)
      endif
c***********************************************************************
c
c   check event is energetically possible
c
c***********************************************************************
      temp2=u+be1+be2-be0
      nck=0
   40 continue
      totke=gaussn(totkm,sigmak)
      if (nck.gt.10) go to 10
      nck=nck+1
      if (totke.gt.temp2) go to 40
c***********************************************************************
c
c      pick excitation from equidistribution of original plus energy bal
c
c***********************************************************************
      temp=(temp2-totke)/a
      ufis(1)=afis(1)*temp
      ufis(2)=afis(2)*temp
c     find total masses, including excitation energies, at evap time
      amcf=rm0+u
      amc1=rm1+ufis(1)
      amc2=rm2+ufis(2)
      amdiff=amcf-amc1-amc2
      if (amdiff.lt.0.d0) go to 40
c***********************************************************************
c
c     amdiff= ekin should be satisfied
c
c***********************************************************************
C

C... Velocity of pre-fission nucleus : vres*bet0
      pres=sqrt(erec*2+2.d0*a*rstms(0)*erec)
      vres=pres/(Erec+a*rstms(0))

C... Momentum of CM system
      redm=rstms(0)*(afis(1)*afis(2))/a
      pcm=sqrt(totke**2+2.d0*redm*totke)
      th  = acos( 2.0d0 * rn(0) - 1.0d0 )
      ph  = 2.0d0 * pi  * rn(0)

      pcmx = pcm* dsin(th) * dcos(ph)
      pcmy = pcm* dsin(th) * dsin(ph)
      pcmz = pcm* dcos(th)

C... Velocity of fission fragment 1 in CM system
      v1x=pcmx/sqrt(pcm**2 + (rstms(0)*afis(1))**2)
      v1y=pcmy/sqrt(pcm**2 + (rstms(0)*afis(1))**2)
      v1z=pcmz/sqrt(pcm**2 + (rstms(0)*afis(1))**2)

C... Velocity of fission fragment 2 in CM system
      v2x=-pcmx/sqrt(pcm**2 + (rstms(0)*afis(2))**2)
      v2y=-pcmy/sqrt(pcm**2 + (rstms(0)*afis(2))**2)
      v2z=-pcmz/sqrt(pcm**2 + (rstms(0)*afis(2))**2)

C... Boost to Lab system
      v1 = sqrt(v1x**2 + v1y**2 + v1z**2)
      v2 = sqrt(v2x**2 + v2y**2 + v2z**2)

      v1x= ( v1x * sqrt(1.d0 - vres**2 * ( bet0(2)**2 + bet0(3)**2)) +
     & vres * bet0(1) ) / (1.d0 + v1 * vres )
      v1y= ( v1y * sqrt(1.d0 - vres**2 * ( bet0(1)**2 + bet0(3)**2)) +
     & vres * bet0(2) ) / (1.d0 + v1 * vres )
      v1z= ( v1z * sqrt(1.d0 - vres**2 * ( bet0(2)**2 + bet0(1)**2)) +
     & vres * bet0(3) ) / (1.d0 + v1 * vres )

      v2x= ( v2x * sqrt(1.d0 - vres**2 * ( bet0(2)**2 + bet0(3)**2)) +
     & vres * bet0(1) ) / (1.d0 + v2 * vres )
      v2y= ( v2y * sqrt(1.d0 - vres**2 * ( bet0(1)**2 + bet0(3)**2)) +
     & vres * bet0(2) ) / (1.d0 + v2 * vres )
      v2z= ( v2z * sqrt(1.d0 - vres**2 * ( bet0(2)**2 + bet0(1)**2)) +
     & vres * bet0(3) ) / (1.d0 + v2 * vres )

C... Kinetic Energy and velocity of fission fragment 1 in Lab system
      gam=1.d0/sqrt(1.d0-v1x**2-v1y**2-v1z**2)
      er(1)= afis(1)*rstms(0) *(gam -1.d0)

      betf(1,1)=v1x/sqrt(v1x**2+v1y**2+v1z**2)
      betf(1,2)=v1y/sqrt(v1x**2+v1y**2+v1z**2)
      betf(1,3)=v1z/sqrt(v1x**2+v1y**2+v1z**2)

C... Kinetic Energy and velocity of fission fragment 2 in Lab system
      gam=1.d0/sqrt(1.d0-v2x**2-v2y**2-v2z**2)
      er(2)= afis(2)*rstms(0) *(gam -1.d0)

      betf(2,1)=v2x/sqrt(v2x**2+v2y**2+v2z**2)
      betf(2,2)=v2y/sqrt(v2x**2+v2y**2+v2z**2)
      betf(2,3)=v2z/sqrt(v2x**2+v2y**2+v2z**2)

      return
c
   50 format ('---> fission failed: ja=',i5,'  jz=',i5/'-       u=',
     1 1pe10.3,'      erec=',e10.3/'    zfis1=',e10.3,'     afis1=',
     2 e10.3/'     zfis2=',e10.3,'     afis2=',e10.3)
   60 format (//'  logic error in fiss.called with fisinh flag',' unse
     1t.',2i10,2f10.5)
      end subroutine



************************************************************************
*                                                                      *
      function fprob (z,a,e)
C/////////////////////////////////////////////////////////////////////
C****************This routine is originally in LAHET code**************
C  FPROB
C    Calculate fission probability
C=====================================================================
C <variables>
C     a    :   mass of nucleus before fission                     (IN)
C     z    :   charge  of nucleus before fission                  (IN)
C     e    :   excitation energy of nucleus before fission        (IN)
C    ifiss :   Fission version number (1 : LAHET original version,
C                                      2 : Iwamoto-Nakano version)(IN)
C  frob    :   fission probability                               (OUT)
C   ef     :   fission barrier                                   (OUT)
C/////////////////////////////////////////////////////////////////////

      implicit doubleprecision(a-h,o-z)
      parameter (inn=150, iiz=98)
      logical isz, isn
      common /cook/ sz(iiz), sn(inn), con(2), amean(240), pz(iiz),
     1 pn(inn), isz(iiz), isn(inn)
!$OMP THREADPRIVATE(/cook/)

      common /std1/ r(70),s(70),sigma,rr(70)
!$OMP THREADPRIVATE(/std1/)
      common /emitr/ gj(70),q(70),V(70),delta(70),smalla(70)
!$OMP THREADPRIVATE(/emitr/)
      common /fiss/ beta,alp
!$OMP THREADPRIVATE(/fiss/)

      real*8 i0,i1
c
c   function to compute the fission probability.
c     for z<90 uses ..........
c    uses statistical model fits.
      parameter (lnfis=18)
      dimension slope(lnfis), anort(lnfis)
      data slope /0.23d0,0.233d0,0.12225d0,0.14727d0,0.13559d0,0.15735d0
     1 ,0.16597d0,0.17589d0,0.18018d0,0.19568d0,0.16313d0,0.17123d0,
     2 6*0.17d0/
      data anort /219.4d0,226.9d0,229.75d0,234.04d0,238.88d0,241.34d0,
     1 243.04d0,245.52d0,246.84d0,250.18d0,254.0d0,257.8d0,261.3d0,
     2 264.8d0,268.3d0,271.8d0,275.3d0,278.8d0/

      data a1, a2, a3 /0.2185024d0,16.70314d0,321.175d0/
      data const /0.3518099d0/
      data c1, c2, c3 /1.089257d0,0.01097896d0,31.08551d0/

      real*8 alpha(2), xi(2,0:3), c_prok(4,3), prok(4)
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)

C   Initialization
      fprob=0.d0
      iz=nint(z)
      ia=nint(a)
      in=nint(a-z)
      u=e

      if (iz.gt.88) then
c
c     high z fission probability
c     for z> 89 and <100 ........
c     use the systematics of vandenbosch & huizenga with a ball park
c     observation that fission probability drops for most nuclei at 6 m

       iz=iz-88
       if (iz.gt.lnfis.or.u.lt.6.0) go to 20
       gamnf=slope(iz)*(a-anort(iz))
       gamnf=1.d1**gamnf
       agoes=1.d0
       go to 10
      endif

      izk = min(iz,iiz)
      ink = min(in,inn)

Calculate -Qn + paring energy + shell correction
      se=energm(iz,ia-1,2)+energm(0,1,2)-energm(iz,ia,2)-sz(izk)-
     & sn(ink-1)
      if ((2*(iz/2).ne.iz).or.(2*(in/2).ne.in)) se=se-pz(izk)-pn(ink-1)
      if ((2*(iz/2).ne.iz).and.(2*(in/2).ne.in)) se=se-pz(izk)-pn(ink-1)
Calculate fission barrier
      x=z*z/a
      ef=x*(a1*x-a2)+a3+se
      if (ef.gt.u) go to 20    !Excited energy is below fission barrier
Calculate level density paameter for neutron emission
      an=0.125d0*(a-1.d0)
      an1=0.5d0/an
      an2=0.25d0*an1/an
Calculate level density parameter for fission
      af=x-c3
      af=an*(c1+c2*af*af)
      a1thrd=a**0.33333333d0
      ss=2.d0*sqrt(an*(u-se))
Calculate Gamma_n/Gamma_f
      if (ss.gt.10.d0) then
C I0=J0
        i0=an1*(ss-1.d0)
C I1=J1
        i1=ss*an2*((ss+ss-6.d0)+6.d0)
        gamnf=const*(((0.76d0*i1-5.d-02*i0)*a1thrd+1.93d0*i1)*a1thrd+
     1       1.66d0*i0)
        s2=2.d0*sqrt(af*(u-ef))
        exps=0.d0
        eag = ss-s2
        if (eag.gt.-150.d0) exps=exp(eag)

       gamnf=gamnf*exps*af/(s2-1.d0)
      else
        e1=exp(ss)
        i0=((ss-1.d0)*e1+1.d0)*an1
        i1=((6.d0+ss*(ss+ss-6.d0))*e1+ss*ss-6.d0)*an2
        gamnf=const*(((0.76d0*i1-5.d-02*i0)*a1thrd+1.93d0*i1)*a1thrd+
     1  1.66d0*i0)
        ss=2.d0*sqrt(af*(u-ef))
        e2=((ss-1.d0)*exp(ss)+1.d0)/af
        gamnf=gamnf/e2
      endif
C....Since we could not find any reference of the following suppression
C     factor for high excited energy fission, i.e. agoes, we did not
C     mention agoes in the paper.
C     But we include agoes because it is used in LAHET code.
      call drein1 (1,s(1),smalla(1),eye1,eye0)
      call drein2 (s(1),smalla(1),eye2)
      if( eye1*eye0 .ne. 0.0d0 ) then
      epsav=(eye2+beta*eye1)/(eye1+beta*eye0)
      agoes=(u-7.d0)/(epsav+7.d0)
      else
      agoes = 1.0
      end if
      if (agoes.lt.1.d0) agoes=1.d0
      goto 10

 10   continue
      if(ifiss .eq. 1) then
       fprob=1.d0/(1.d0+gamnf)/agoes
      elseif(ifiss .eq. 2) then
       xi(1,0) = -42.181d+00
       xi(1,1) =   7.0882d-01
       xi(1,2) =  -3.9383d-03
       xi(1,3) =   7.2523155d-06
       xi(2,0) = 119.13d+00
       xi(2,1) =  -1.9639d+00
       xi(2,2) =   1.0909d-02
       xi(2,3) =  -2.0086d-05

       do k = 1, 2
          alpha(k) = 0.0d0
          do l = 0, 3
              alpha(k) = alpha(k) + xi(k,l) * a**l
          end do
       end do

       ep = alpha(1) * e**(alpha(2))
       if( ep .lt. 0 ) ep = 1.0e-6

       delta_mass    = 0.0443 * e - 0.982
       delta_charge  = 0.0221 * e - 0.977
       delta_neutron = 0.0222 * e - 0.052

       a_fic = a + delta_mass
       z_fic = z + delta_charge

       eta = z_fic * z_fic / a_fic

      ! Table 2, Prokofiev, NIMA 463 (2001), 557
       c_prok(1,1) =  486.9316d0
       c_prok(1,2) =   16.7668d0
       c_prok(1,3) = -180.1085d0
       c_prok(2,1) =   11.0081d0
       c_prok(2,2) =   -0.1455d0
       c_prok(2,3) =   -0.4247d0
       c_prok(3,1) =  -17.9890d0
       c_prok(3,2) =    0.9642d0
       c_prok(3,3) =   -3.1994d0
       c_prok(4,1) =   -1.1188d0
       c_prok(4,2) =    0.0352d0
       c_prok(4,3) =    0.0000d0

       do k=1,3
         prok(k) =
     &        + exp( c_prok(k,1)
     &        + c_prok(k,2)*eta
     &        + c_prok(k,3)*sqrt(eta) )
       end do
       prok(4) = c_prok(4,1)+c_prok(4,2)*eta

       sigmaf = prok(1) * ( 1.0 - exp( -prok(3)*(ep-prok(2) ) ) )
     &  * (1.0-prok(4)*log(ep))
       if( sigmaf .lt. 0.0 ) sigmaf = 0.0

       CALL sigrc(itype,ep,int(a_fic),int(z_fic),sigt,sigr,sigs) ! itype=1: proton
       sigmar = sigr*1000.0
       fission_prob = sigmaf/sigmar
       if( fission_prob .gt. 1. ) fission_prob = 1.

      ! supression factor by hiroki iwamoto
       gamma  = 43.0
       gamma0 = (32.32692/eta)**3.0/2
       gamma1 = (32.32692/eta)**3.5*gamma
       gamma2 = 0.08
       sfactor = ( 1.0/(1.0+exp(-gamma2*(e-gamma1)))-gamma0 )*0.34

       if( sfactor .lt. 0.0 ) sfactor =0.0

       fprob = fission_prob*sfactor
      endif
 20   continue
      return
      end function

************************************************************************
*                                                                      *
      function efms(z,a)
C/////////////////////////////////////////////////////////////////////
C  EFMS
C  Fission barrier given by Myer & Swaiteski (PRC60,014606,1999)
C=====================================================================
C <variables>
C     a   :   the mass of a fissioning nucleus      (IN)
C     z   :   the charge of a fissioning nucleus    (IN)
C   efms  :   fission barrier  [MeV]                (OUT)
C/////////////////////////////////////////////////////////////////////

      implicit doubleprecision(a-h,o-z)

      common /fiss4/ shell(150,250)
!$OMP THREADPRIVATE(/fiss4/)

C... 8/15/1999
      parameter(x0=48.5428d0, x1=34.15d0)
      f1(t)=1.99749d-4*(x0-t)**3
      f2(t)=5.95553d-1-0.124136*(t-x1)

      efms=0.d0
      iz=nint(z)
      ia=nint(a)

      c=1.9+(z-80)/75
      ai=1.-2*(z/a)
      xx=1-c*ai**2
      ss=a**.66667*xx
      x=z**2/a/xx
      if(ia-iz.le.0.or.ia-iz.gt.250.or.iz.gt.150.or.iz.lt.1) then
       sh=0.d0
      else
       sh=shell(iz,ia-iz)
      endif
      if(x.ge.x1.and.x.le.x0) then
       efms=ss*f1(x)-sh
      else if(x.ge.20.and.x.lt.x0)then
       efms=ss*f2(x)-sh
      else
       efms=-1.0
      endif

      return
      end function


************************************************************************
*                                                                      *
      function fprob2 (z,a,e,sigma)
C/////////////////////////////////////////////////////////////////////
C/////////////////////////////////////////////////////////////////////
C  FPROB2
C    Calculate fission decay width r(nimax-1)
C    based on Bohr-Wheeler formula modified by K.-H. Shmidt
C=====================================================================
C <variables>
C     a    :   mass of nucleus before fission                     (IN)
C     z    :   charge  of nucleus before fission                  (IN)
C     e    :   excitation energy of nucleus before fission        (IN)
C     sigma:   hadronic evaporation decay width                   (IN)
C  frob2   :   fission probability                               (OUT)
C/////////////////////////////////////////////////////////////////////
      use levdenmod, only : constT, rho_levden
      use NGSDATAMOD
      implicit doubleprecision(a-h,o-z)
      parameter( pi = 3.141592653589793d0 )

      sp = 0.d0 ! consider spin of nucleus someday

      Ba   = fisBarr(int(z),int(a),1) ! Inner fission barrier
      Bb   = fisBarr(int(z),int(a),2) ! Outer fission barrier. Currently the same as Ba
      Bmax = max(Ba,Bb)

      if(u .le. Bmax) then ! Do not consider sub barrier fission
       fprob2 = 0.d0
       return
      endif

      Tct  = constT(int(z),int(a-z)) ! constant temperature
      Tf   = 0.01d0/( Log(rho_levden(int(a),int(z),u-Bmax+0.01d0,sp,2))
     &              - Log(rho_levden(int(a),int(z),u-Bmax       ,sp,2))) ! fermi gas temperature (at barrier?).
      Tequi= 1.d0/(2.d0*pi)         ! temperature related to barrier frequency. Barrier curverture / 2pi
      x    = z**2 / a               ! fissility parameter
      splim = 7.93d0 * sqrt(Tf / Tct / (1.d0 - x) )
      Frot = exp((sp/splim)**2)

      Fa   = 1.d0 / ( 1.d0 + exp( -(e-Ba) / Tequi ) )
      Fb   = 1.d0 / ( 1.d0 + exp( -(e-Bb) / Tequi ) )
      Ga   = Fa * 0.14d0 / sqrt( pi / 2.d0 )
      Gb   = Fb / 2.d0

      integ = 0.d0
      de = (u - Bmax)/1000.d0
      do i = 1, 1000
       integ = integ + rho_levden(int(a),int(z),u-Bmax-de*dble(i),sp,2) ! this level must be replaced by saddle level density
      enddo

      fiswid = Frot /
     & (Ga * exp(( Ba - Bmax )/ Tf) + Gb * exp(( Bb - Bmax )/ Tf))
     & / rho_levden(int(a),int(z),u,sp,2)
     & * integ

      fprob2 = fiswid / amu * hbarc**2 * pi
     & * rho_levden(ia,iz,u,sp,2)

      end function

      end module
