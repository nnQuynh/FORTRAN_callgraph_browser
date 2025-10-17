************************************************************************
*                                                                      *
      subroutine xstgmm(mm,sigt,ein,mk)
*                                                                      *
*     calculate cross section for the muon pair production.            *
*     Last modified by Y.Sakaki on 2019/07/02                          *
*                                                                      *
************************************************************************
      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      use moddas_material

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      common /kmat1g/ kmat(kvlmax)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)

*-----------------------------------------------------------------------
*     initialization of this material and particle
*-----------------------------------------------------------------------

      erg  = ein
      imat = mk
      averaged_xsec = 0.0d0

*-----------------------------------------------------------------------
*     expand natural nucleus and calculate average x-section
*-----------------------------------------------------------------------
      ! "jmd(1+mk)" is basically 1
      ! "jmd(1+mk+1) - 1" is the total number of RIs in the considering material.
      ! For example, considering a material "H 2 O 1", it's water,
      ! "jmd(1+mk+1) - 1" is 3, because registered RIs in PHITS are H-1, H-2 and O-16.
      ! "fme(m)" is somehow a contribution ratio for each RI;
      !  fme(m) =
      !     0.66656 (for H-1),
      !     0.00010 (for H-2),
      !     0.33333 (for O-16).
      ! Therefore this number is the product of [Atomic nuber ratio in the material]
      ! and [Abundance ratio for each atom].
      ! Note that the contribution from "H" is double compared to "O",
      ! since the number of "H" in water is double compared to "O".
      ! We can get the atomic number (Z) and the mass number (A) as below.
      ! The point to pay attention is that some elements have zero mass-number (A==0).
      ! For example it's carbon (C).
      ! How to get the mass number and the abundance ratio for such "zero-mass" elements like C is different from "non-zero-mass" elements.
      ! See the detail of how to do it below "if( ims .eq. 0 )".
      do m = jmd(1+mk), jmd(1+mk+1) - 1
       ie = lme(2,m)
       izm = iza(m) / 1000       ! Z: proton number
       ims = iza(m) - 1000 * izm ! A: mass number
       icntm = m - jmd(1+mk) + 1
       ! write(*,"(I5,I5,f12.6)") izm, ims, fme(m)

       ! -------------------------
       ! For "zero-mass" elements.
       if( ims .eq. 0 ) then

        lemm = nint( dnel_das(kmat0+imat) )
        sekt = 0.0d0 ! total of Unnormalized Abundance Ratio
        sekc = 0.0d0 ! averaged cross section with the Unnormalized Abundance Ratio
        do lem = 1, lemm
         ! This loop runs not only for elements having "ims==0" but also "ims!=0".
         ! In order to count up xsc only for elements having "ims==0",
         ! we need the following condition "izt .eq. izm".
         izt = zz_das(kmat(imat)+lem)
         if( izt .eq. izm ) then
           unnorm_AR = den_das(kmat(imat)+lem) ! Unnormalized Abundance Ratio
           sekt = sekt + unnorm_AR
           ims = nint( a_das(kmat(imat)+lem) ) ! mass number
           xsc = getxstgmm_sum(erg,izm,ims)
           sekc = sekc + xsc * unnorm_AR
           ! write(*,"(I10,I5,f10.6)") izt,ims,unnorm_AR
         end if
        end do
        xsc = sekc / sekt ! correct to the averaged cross section with "normalized" abundance ratio

       !-------------------------
       ! For "non-zero-mass" elements.
       else ! here
        xsc = getxstgmm_sum(erg,izm,ims)
       end if

       averaged_xsec = averaged_xsec + xsc * fme(m)

      end do

*-----------------------------------------------------------------------

      sigt  = averaged_xsec

*-----------------------------------------------------------------------
      return
      end









************************************************************************
*                                                                      *
      function getxstgmm_sum(erg,iZ,iA)
*                                                                      *
*     calculate cross section for the muon pair production.            *
*     The contribution from all scattering modes are summed up         *
*     The unit is barn.                                                *
*     erg: incident photon energy [MeV]                                *
*     iZ: atomic number of a target particle                           *
*     iA: mass   number of a target particle                           *
*     Exact Born Formula is used for the estimation.                   *
*     Last modified by Y.Sakaki on 2019/08/14                          *
*                                                                      *
************************************************************************
      implicit real*8 (a-h,o-z)
      igmm_coh   = 1 ! sakaki parameter
      igmm_quasi = 1
      igmm_inel  = 1
      tmp=0.0d0
      if(igmm_coh   .ne. 0) tmp = tmp + getxstgmm(1,erg,iZ,iA)
      if(igmm_quasi .ne. 0) tmp = tmp + getxstgmm(2,erg,iZ,iA)
      if(igmm_inel  .ne. 0) tmp = tmp + getxstgmm(3,erg,iZ,iA)
      getxstgmm_sum=tmp
      return
      end
************************************************************************
*                                                                      *
      function getxstgmm(isw,erg,iZ,iA)
*                                                                      *
*     calculate cross section for the muon pair production.            *
*     The unit is barn.                                                *
*     isw=1: coherent scattering                                       *
*        =2: quasi-elastic scattering                                  *
*        =3: inelastic scattering                                      *
*     erg: incident photon energy [MeV]                                *
*     iZ: atomic number of a target particle                           *
*     iA: mass   number of a target particle                           *
*     Exact Born Formula is used for the estimation.                   *
*     Last modified by Y.Sakaki on 2019/08/14                          *
*                                                                      *
************************************************************************
      implicit real*8 (a-h,o-z)
      real*8 lnk,lnk0
      parameter ( pi = 3.141592653589793d0 )
      parameter ( rm = 105.6583715d0 )    ! muon mass [MeV/c**2]
      parameter ( rpmass = 938.272046d0 ) ! proton mass [MeV/c**2]
      parameter ( rnmass = 939.565379d0 ) ! neutron mass [MeV/c**2]
      parameter ( hbarc = 197.32696d-13 ) ! MeV -> cm


      ! conversion constant from "micro barn" to "barn"
      conv = 1.0d-6


      ! set the mass of target particle
      if(isw.eq.1) then ! coherent scattering
        rmi=rpmass*iA   ! or rpmass*iZ + rnmass*(iA-iZ)
      else if(isw.eq.2) then ! quasi-elastic scattering
        rmi=rpmass
      else if(isw.eq.3) then ! inelastic scattering
        rmi=rpmass
      else
        print*,"set correct isw"
        stop
      endif


      ! require threshold energy
      rk   = erg
      rk0= 2.*rm*(rm+rmi)/rmi ! threshold energy of muon pair production
      if(rk.le.rk0) then
        getxstgmm=0.0d0
        return
      endif


      ! etc
      if(iZ.gt.iA) then
        print*,"Error: iZ > iA"
       !stop
      endif


      ! ----------------------------------------------------------------
      ! calculate the cross section
      ! ----------------------------------------------------------------
      ! useful variables
      lnk =log(rk /1000.) ! lnk = ln(k/GeV)
      lnk0=log(rk0/1000.)

      ! ----------------------------------------------------------------
      ! coherent scattering
      if(isw.eq.1) then
        ! (low energy fit)  ----------------
        if(lnk.lt.1) then
          if(iA.lt.20) then
            c =fit_gmm1(iA, -4.50253,5.18323,-0.0294001)
            p1=fit_gmm1(iA, 1.8709,0.191658,0.435094)
          else
            c =fit_gmm1(iA, -0.0911805,1.06669,-0.383089)
            p1=fit_gmm1(iA, 1.45742,0.51752,0.256067)
          endif
          getxstgmm=conv*(iZ**2)*0.01*c*(lnk-lnk0)**p1
          return
        ! (high energy fit) ----------------
        else
          if(iA.lt.20) then
            c =fit_gmm1(iA, -18.6336,22.1782,-0.0119774)
            p1=fit_gmm1(iA, 0.846844,0.218372,0.1153)
            p2=fit_gmm1(iA, 0.233853,-0.0770553,0.52948)
          else
            c =fit_gmm1(iA, -13.7519,17.5143,-0.0195081)
            p1=fit_gmm1(iA, 2.89131,-1.89271,-0.0282275)
            p2=fit_gmm1(iA, -6.87297,7.33152,-0.0278267)
          endif
          getxstgmm=conv*(iZ**2)*0.01*c*(lnk**p1+p2)
          return
        ! (end) ----------------------------
        endif

      ! ----------------------------------------------------------------
      ! quasi-elastic scattering
      else if(isw.eq.2) then
        ! elastic cross section for proton
        tmp = 0.00614518+0.00144848*exp(-1.06378*lnk)
        tmp = tmp +0.0246819*tanh(0.400696*lnk)
        xst_p_el=tmp
        ! elastic cross section for neutron
        tmp = 0.000205765+0.000128124*exp(-0.788059*lnk)
        tmp = tmp +0.00103824*tanh(0.462178*lnk)
        xst_n_el=tmp
        getxstgmm=conv*(iZ*xst_p_el + (iA-iZ)*xst_n_el)
        return

      ! ----------------------------------------------------------------
      ! inelastic scattering
      else if(isw.eq.3) then
        ! inelastic cross section for proton
        tmp = -0.000837713
        tmp = tmp +0.00320051 *exp((-0.910448-0.102039*lnk)*lnk)
        tmp = tmp +0.031662 *asinh(0.194865*lnk)
        xst_p_inel = tmp
        getxstgmm=conv*iA*xst_p_inel
        return

      ! ----------------------------------------------------------------
      else
        print*,"set correct isw"
        stop
      endif

      return
      end
************************************************************************
*                                                                      *
      function getxstgmm_max(isw,erg,iZ,iA)
*                                                                      *
*     calculate maximum value of the cross section for the muon pair   *
*     production, dsigma/(dp*dphi*deta*du*dmf2).                       *
*     The unit is (1/MeV2)/(MeV*1*1*1*MeV2)=1/MeV5                     *
*     isw=1: coherent scattering                                       *
*        =2: quasi-elastic scattering                                  *
*        =3: inelastic scattering                                      *
*     erg: incident photon energy [MeV]                                *
*     iZ: atomic number of a target particle                           *
*     iA: mass   number of a target particle                           *
*     Exact Born Formula is used for the estimation                    *
*     Last modified by Y.Sakaki on 2019/08/14                          *
*                                                                      *
************************************************************************
      implicit real*8 (a-h,o-z)
      real*8 lnk,lnk0
      parameter ( pi = 3.141592653589793d0 )
      parameter ( rm = 105.6583715d0 )    ! muon mass [MeV/c**2]
      parameter ( rpmass = 938.272046d0 ) ! proton mass [MeV/c**2]
      parameter ( rnmass = 939.565379d0 ) ! neutron mass [MeV/c**2]
      parameter ( hbarc = 197.32696d-13 ) ! MeV -> cm


      ! set the mass of target particle
      if(isw.eq.1) then ! coherent scattering
        rmi=rpmass*iA   ! or rpmass*iZ + rnmass*(iA-iZ)
      else if(isw.eq.2) then ! quasi-elastic scattering
        rmi=rpmass
      else if(isw.eq.3) then ! inelastic scattering
        rmi=rpmass
      else
        print*,"set correct isw"
        stop
      endif


      ! require threshold energy
      rk   = erg
      rk0= 2.*rm*(rm+rmi)/rmi ! threshold energy of muon pair production
      if(rk.le.rk0) then
        getxstgmm=0.0d0
        return
      endif


      ! etc
      if(iZ.gt.iA) then
        print*,"Error: iZ > iA"
      endif


      ! ----------------------------------------------------------------
      ! calculate the cross section
      ! ----------------------------------------------------------------
      ! useful variables
      lnk =log(rk /1000.) ! lnk = ln(k/GeV)
      lnk0=log(rk0/1000.)

      ! ----------------------------------------------------------------
      ! coherent scattering
      if(isw.eq.1) then
        ! (low energy fit)  ----------------
        if(lnk.lt.1) then
          if(iA.lt.20) then ! (A<20)
            p1=fit_gmm1(iA, -9.18295,11.9109,-0.0619237)
            p2=fit_gmm1(iA, -7.25862,7.8054,-0.013548)
          else !              (A>=20)
            p1=fit_gmm1(iA, -0.0400572,2.88225,-0.426769)
            p2=fit_gmm1(iA, -0.223057,0.758513,-0.16994)
          endif
          tmp=(iZ**2)*(1.0d-14)*p1*(lnk-lnk0)*exp(-p2*(lnk-lnk0)**2)
          getxstgmm_max=tmp
          return
        ! (high energy fit) ----------------
        else
          if(iA.lt.20) then ! (A<20)
            p1=fit_gmm1(iA, -5.46223,6.45955,-0.0159198)
            p2=fit_gmm1(iA, 1.16374,-0.0070101,0.324688)
          else !              (A>=20)
            p1=fit_gmm1(iA, -5.70449,6.84888,-0.0225444)
            p2=fit_gmm1(iA, 0.211317,0.953672,-0.00665326)
          endif
          tmp=(iZ**2)*(1.0d-14)*(p1+0.7*lnk)*exp(-p2*lnk)
          getxstgmm_max=tmp
          return
        ! (end) ----------------------------
        endif

      ! ----------------------------------------------------------------
      ! quasi-elastic scattering
      else if(isw.eq.2) then
        ! (low energy fit)  ----------------
        if(lnk.lt.1) then
          p1=1.20236
          p2=1.161
          p3=0.40472
          xst_p_el=(1.0d-14)*p1*(lnk-lnk0)*exp(-p2*lnk-p3*lnk**2)
          p1=0.168811
          p2=0.809264
          p3=0.571756
          xst_n_el=(1.0d-14)*p1*(lnk-lnk0)*exp(-p2*lnk-p3*lnk**2)
          getxstgmm_max = iZ*xst_p_el + (iA-iZ)*xst_n_el
          return
        ! (high energy fit) ----------------
        else
          p1=1.90531
          p2=0.945597
          xst_p_el=(1.0d-14)*p1*exp(-p2*lnk)
          p1=0.30888
          p2=0.936533
          xst_n_el=(1.0d-14)*p1*exp(-p2*lnk)
          getxstgmm_max = iZ*xst_p_el + (iA-iZ)*xst_n_el
          return
        ! (end) ----------------------------
        endif

      ! ----------------------------------------------------------------
      ! inelastic scattering
      else if(isw.eq.3) then
        ! (low energy fit)  ----------------
        if(lnk.lt.1) then
          p1=1.82643
          p2=1.48225
          p3=0.113904
          xst_p_inel=(1.0d-20)*p1*(lnk-lnk0)*exp(-p2*lnk-p3*lnk**2)
          getxstgmm_max = iA*xst_p_inel
          return
        ! (high energy fit) ----------------
        else
          p1=3.64122
          p2=1.14994
          xst_p_inel=(1.0d-20)*p1*exp(-p2*lnk)
          getxstgmm_max = iA*xst_p_inel
          return
        ! (end) ----------------------------
        endif

      ! ----------------------------------------------------------------
      else
        print*,"set correct isw"
        stop
      endif

      return
      end
************************************************************************
      function fit_gmm1(iA,b1,b2,b3)
************************************************************************
      implicit real*8 (a-h,o-z)
      real b1,b2,b3
      fit_gmm1=b1+b2*iA**b3
      return
      end
************************************************************************
      subroutine gen_gmm(mtd,isw,ein,iZ,iA,            event)
************************************************************************
      implicit real*8 (a-h,o-z) ! i,j,k,l,m,n
      parameter ( pi = 3.141592653589793d0 )
      real*8 :: event(20)
      character(len=*) :: mtd
      c1=0.
      c2=ein
      c3=0.
      c4=pi
      call gen_gmm_cut(mtd,isw,ein,iZ,iA,c1,c2,c3,c4,event)
      return
      end
************************************************************************
      subroutine gen_gmm_cut(mtd,isw,ein,iZ,iA,c1,c2,c3,c4,event)
      ! generate 4-momenta of final states in gamma-->mu+mu-
      ! p = absolute value of 3 momentum for muon-
      ! theta = emission angle for muon-
      ! isw = 1:coherent, 2:quasi-elastic, 3: inelastic, 4: deep-inelastic
      ! ein = incident photon energy [MeV]
      ! iZ = Z
      ! iA = A
      ! c1 < p     < c2
      ! c3 < theta < c4
************************************************************************
      implicit real*8 (a-h,o-z) ! i,j,k,l,m,n
      parameter ( pi = 3.141592653589793d0 )
      parameter ( rm = 105.6583715d0 )    ! muon mass [MeV/c**2]
      parameter ( rpmass = 938.272046d0 ) ! proton mass [MeV/c**2]
      parameter ( rnmass = 939.565379d0 ) ! neutron mass [MeV/c**2]
      parameter ( alpha = 7.297353d-3 )   ! fine-structure constant
      parameter ( hbarc = 197.32696d-13 ) ! MeV -> cm
      real*8 :: event(20)
      real*8 :: dsigma_max
      integer :: nPass
      logical is_elastic, is_inelastic
      logical is_Born, is_WW
      character(len=*) :: mtd


      ! ----------------------------------------
      ! setup
      iverbose=0

      if     (mtd .eq. "Born") then
        is_Born=.true.
      else if(mtd .eq. "WW"  ) then
        is_Born=.false.
      else
        print*,"set correct mtd [gen_gmm]"
        call abort ()
      endif
      is_WW=.not.is_Born

      do i = 1, 20
      event(i)=-1.0d+10
      end do


      ! ----------------------------------------
      ! phase space cut
      rpmin=c1       ! 0.
      rpmax=c2       ! ein
      rthetamin=c3       ! 0.
      rthetamax=c4       ! pi

      rdp=rpmax-rpmin
      rp_mid=(rpmax+rpmin)/2.
      rdtheta=rthetamax-rthetamin
      rtheta_mid=(rthetamax+rthetamin)/2.


      ! ----------------------------------------
      ! initial parameters

      rk=ein ! k: incident photon energy [MeV]
      rmi=0.0d0 ! mi
      if ( isw.eq.1 ) rmi=rpmass*iZ+rnmass*(iA-iZ) !rmi=rpmass*iA
      if ( isw.eq.2 ) rmi=rpmass
      if ( isw.eq.3 ) rmi=rpmass
      if ( isw.eq.4 ) rmi=rpmass

      safety_factor = 3.0d0
      dsigma_max = safety_factor * getxstgmm_max(isw,rk,iZ,iA)



      ! ----------------------------------------
      ! useful variables
      rm2=rm*rm
      rmi2=rmi*rmi
      rk2=rk*rk
      rsqrts=sqrt(rmi2 + 2.*rmi*rk)

      if (isw.eq.1 .or. isw.eq.2 ) then
        is_elastic  =.true.
      else
        is_elastic  =.false.
      endif
      is_inelastic=.not.is_elastic



      ! ----------------------------------------
      ! etc
      if ( iZ.gt.iA ) then
        print*,"[error]: Z>A"
      end if




      ! ----------------------------------------
      ! MC parameters range and Area
      rphimin=0.
      rphimax=2.*pi
      if(is_Born) then
      tmp1 = 2.*rk2*rmi -4.*rm2*(rk+rmi)
      tmp2 = 2.*rk*sqrt(4.*rm**4+rk2*rmi2-4.*rm2*rmi*(rk+rmi))
      epstmp=4.*rm2/(rk2*rmi2)*(rk*rmi+rmi2-rm2)
      if(epstmp<0.05) then
      tmp=    -2.*(rm**8)
      tmp=tmp +6.*(rm**6)* rmi    *(rk+rmi)
      tmp=tmp -1.*(rm**4)*(rmi**2)*(5.*rk2+12.*rk*rmi+6.*rmi2)
      tmp=tmp +2.*(rm**2)*(rmi**4)*(rk+rmi)*(2.*rk+rmi)
      tmp=tmp +rk2*(rmi**5)*(2.*rk+rmi)
      tmp=tmp*4.*(rm**4)/( (rk**4)*(rmi**5)*(2.*rk+rmi) )
      rtminMC= tmp
      else
      rtminMC= ( tmp1 - tmp2 )/( 2.*rk+rmi )
      endif
      rtmaxMC= ( tmp1 + tmp2 )/( 2.*rk+rmi )
      else ! WW method
      rtminMC=4.*rm**4/ein**2
      rtmaxMC=rm2*(1.+ein**2*pi**2/rm2)**2
      endif
      rumin=log(rtminMC)
      rumax=log(rtmaxMC)
      rmf2min=rmi2
      rmf2max=(rsqrts-2.*rm)**2



      ! ----------------------------------------------------------------
      ! generate
      ! ----------------------------------------------------------------
      if      (isw.eq.1) then; ntrial=20000;
      else if (isw.eq.2) then; ntrial=20000;
      else if (isw.eq.3) then; ntrial=1000000;
      else; print*,"set correct isw [gen_gmm]"; stop;
      endif
      ! ----------------------------------------------------------------

    ! ntrial=100 ! for test
      do i = 1, iabs(ntrial)


        ! ----------------------
        ! generate kinetic variables

        rp=randmy(rpmin,rpmax)
        retamin=-0.5/( 1.+(rp**2+rm2)*(rthetamin**2)/rm2 )
        retamax=-0.5/( 1.+(rp**2+rm2)*(rthetamax**2)/rm2 )
        reta=randmy(retamin,retamax)
        rphi=randmy(rphimin,rphimax)
        ru  =randmy(rumin  ,rumax  )


        if ( is_inelastic ) then
          rmf2=randmy(rmf2min,rmf2max)
        else
          rmf2=rmi2
        endif


        if ( ntrial .lt. 0 ) then ! fix kinetic variables for test
          rp=1000.0d0 ! p=1GeV
          reta=-0.5/( 1.+(rp**2+rm2)*(0.01)**2/rm2 ) ! theta=0.01
          ru  =log(0.01*1000**2) ! t=0.01GeV^2
          rmf2=rmi2*1.1 ! mf2
        endif




        ! useful variables
        rE=sqrt(rp**2+rm2) ! muon- energy
        rtheta=rm/rE*sqrt( 1./(-2.*reta)-1. ) ! theta
        cost=cos(rtheta)
        sint=sin(rtheta)
        rmf=sqrt(rmf2)
        rt=exp(ru)
        rl=(rE**2)*(rtheta**2)/rm2
        rx=rE/rk


        ! ----------------------
        ! phase space conditions
        ! @@@ p
        tmp=rmi*(2.*rk+rmi)-rmf*(2.*rm+rmf)
        tmp1= -8.*(rm**3)*rmf + tmp**2 -4.*(rm2)*(rmf2+tmp)
        tmp1=tmp1 -4.*(rk2)*(rm2)*sint**2
        if ( tmp1 .le. 0 ) cycle
        tmp1=rk*tmp*cost + (rk+rmi)*sqrt(tmp1)
        tmp1=tmp1/2./(2.*rm*rmf+(rmf2)+tmp+(rk2)*sint**2)
        rpmax_=tmp1
        if ( rp .gt. rpmax_ ) cycle

        ! @@@ mf2
        rkp=rE*rk -cost*rk*rp
        U2 = -2.*rkp +rm2 +2.*(-rE+rk)*rmi +rmi2
        if ( U2 .lt. 0 ) cycle
        sqrtU2 = sqrt(U2)
        if ( rmi .gt. sqrtU2-rm ) cycle
        if ( rmf .gt. sqrtU2-rm ) cycle

        ! @@@ t
        rEps=(rm2 -rmf2 +U2) / sqrtU2 / 2.
        rEs =( rkp -rm2 +rE*rmi) / sqrtU2
        rks =(-rkp      +rk*rmi) / sqrtU2
        rpis=(rmi*sqrt(rk2 -2.0*cost*rk*rp +rp**2)) / sqrtU2
        rpps=sqrt(rEps**2 - rm2)
        tmin =2.*rkp +2.*rEps*(-rEs+rks) -2.*rm2 -2.*rpis*rpps
        tmax =2.*rkp +2.*rEps*(-rEs+rks) -2.*rm2 +2.*rpis*rpps
        tminp=((1.+rl)**2)*rm**4/(4.*(rk2)*((1.-rx)**2)*(rx**2))
        tup  =((1.+rl)**2)*rm**2
        if(is_WW) then
          if(rt .gt. tup ) cycle
          if(is_elastic   .and. tminp.gt. rt) cycle
          if(is_inelastic .and. tmin .gt. rt) cycle
        else ! if Born
          if(tmin .gt. rt .or. rt .gt. tmax) cycle
        endif


        ! --------------------------------------------------------------
        ! calculating differential cross section


        ! if Born --------
        if(is_Born) then
        q2=-rt
        Delta=(rmf2-rmi2)/(2.*rmi)
        costp=( -rt-(2.*rm2-2.*rkp-2.*rEps*(rks-rEs)) )/(2.*rpps*rpis)

        costk=(rks-rEs)/rpis +(rkp)/(rks*rpis)
        rH0=-rm2*(0.5*q2*(1.-2.*rE/rmi) +2.*rE**2 +2.*rE*Delta)
        rH1= rm2*(2.*rm2+q2)
        W=rEps-rpps*costp*costk
        Y=sqrt( rm2*(1.-costk**2)+(rpps*costp-rEps*costk)**2 )

        B0_=(2.*rE*(rE-rk) +0.5*q2*((rk-2.*rE)/rmi+1.)+(2.*rE-rk)*Delta)
        B0_=-(2./rkp)*( (rm2-q2/2.)*B0_ -0.5*q2*rk2 )
        B0_=B0_ +(q2/rmi)*(rmi+rE-rk-0.5*q2/rmi)
        B0_=B0_ -2.*Delta*(Delta-rk+rE-q2/rmi) +rkp
        B1=-( (q2**2-4.*rm**4)/rkp +2.*q2 +2.*rkp +4.*rm2 )
        C0_=-(rm2/(rkp**2))*(2.*(rk-rE-Delta+q2/(2.*rmi))*(rk-rE)+q2/2.)
        C0_=C0_+(1./rkp)*(q2*(1.-rE/rmi)+2.*rE*Delta)
        C1_=rm2*(2.*rm2+q2)/(rkp**2) -2.*(2.*rm2+q2)/rkp
        D0_= 1./rkp
        D1_=-2./rkp

        XX2=(rH0*W)/( (Y**3)*(rks**2) ) +B0_/(Y*rks) +C0_ +D0_*rks*W
        XX1=(rH1*W)/( (Y**3)*(rks**2) ) +B1/(Y*rks) +C1_ +D1_*rks*W
        XX2=XX2 * W2_em(isw,rt,rmf2,iZ,iA)
        XX1=XX1 * W1_em(isw,rt,rmf2,iZ,iA)
        XXX=XX2+XX1
        XXX=-XXX*(alpha**3)/(2.*pi)/(2.*rpis*rpps)
        XXX=XXX*(rp**2)/(sqrtU2*rk*rE)*rpps/(rt**2)
        XXX=XXX*(sint/rtheta)*(rm2*(1+rl)**2/rE**2)*rt ! jacobian
        dsigma = XXX
        if(dsigma > dsigma_max) then
          print*,"[caution] --------------------------------------"
          print*,"[caution] dsigma > dsigma_max."
          print*,"[caution] Larger 'safety_factor' maybe suitable."
          print*,"[caution] dsigma     = ",dsigma
          print*,"[caution] dsigma_max = ",dsigma_max
          print*,"[caution] ",rp
          print*,"[caution] ",reta
          print*,"[caution] ",ru
          print*,"[caution] ",rmf2
          print*,"[caution] --------------------------------------"
        endif


        ! if WW method --------
        else
        if (is_elastic  ) tmin0=tminp
        if (is_inelastic) tmin0=tmin
        YYY = 2.*(alpha**3)/pi/rk/rm2 *sint/rtheta
        YYY = YYY*( 1.-2.*rx*(1.-rx)*( (1.+rl**2)/(1.+rl)**2 ) )
        XXX =      1./rt*(rt-tmin0)*W2_em(isw,rt,rmf2,iZ,iA)
        XXX = XXX +1./rt*(2.*tminp)*W1_em(isw,rt,rmf2,iZ,iA)
        XXX = XXX/(2.*rmi)
        dsigma = YYY*XXX
        endif


        if(dsigma/dsigma_max .gt. randmy(0.0d0,1.0d0)) then
          ! rotation will be applied in outside of here
          ! p1: muon-
          p1E=rE
          p1x=rp*sint*cos(rphi)
          p1y=rp*sint*sin(rphi)
          p1z=rp*cost
          ! p2: muon+
          rpf2=( rt+(rmf-rmi)**2 )*( rt+(rmf+rmi)**2 )/(4.*rmi2)
          if(rpf2 .lt. 0.0d0 .and.iverbose.gt.0) then
            print*,"[CAUTION] rpf2 .lt. 0.0d0: ",rpf2
          endif
          if(rpf2 .lt. 0.0d0) cycle
          rpf=sqrt(rpf2)
          rEf=sqrt(rpf2+rmf2)
          rWWW=(rmf2+rmi2)/2.+rEf*(rE-rk-rmi)+rk*rmi-rE*(rk+rmi)+rk*p1z
          pL=rk-p1z
          rootin= -rk**3*p1z -rk*rpf2*p1z +rpf2*(p1x**2-pL*p1z)
          rootin= rootin+rk2*(rpf2+p1z*(pL+p1z)) -rWWW**2
          if(rootin .lt. 0.0d0) cycle
          costf1=(-pL*rWWW+dabs(p1x)*sqrt(rootin))/(rpf*(pL**2+p1x**2))
          costf2=(-pL*rWWW-dabs(p1x)*sqrt(rootin))/(rpf*(pL**2+p1x**2))
          sintf1=sqrt(1.0d0-costf1**2)
          sintf2=sqrt(1.0d0-costf2**2)
          p2x1=   -p1x -rpf*sintf1
          p2x2=   -p1x -rpf*sintf2
          p2y1=   -p1y
          p2y2=   -p1y
          p2z1=rk -p1z -rpf*costf1
          p2z2=rk -p1z -rpf*costf2
          p2E1=sqrt(p2x1**2 +p2y1**2 +p2z1**2 +rm2)
          p2E2=sqrt(p2x2**2 +p2y2**2 +p2z2**2 +rm2)
          tmp1=dabs(rk+rmi-rE-rEf-p2E1)
          tmp2=dabs(rk+rmi-rE-rEf-p2E2)
          select=0; small=1.0d+10;
          if(tmp1 .lt. tmp2) then; select=1; small=tmp1;
          else                   ; select=2; small=tmp2;
          endif
          if(small.gt.1.0d-6) cycle
          if(select.eq.1) then; p2E=p2E1; p2x=p2x1; p2y=p2y1; p2z=p2z1;
          else                ; p2E=p2E2; p2x=p2x2; p2y=p2y2; p2z=p2z2;
          endif

          ! rotate azimuthal angles along the incident photon direction
          rphi_ = randmy(0.0d0,2.0d0*pi)
          ! p1 (muon-)
          tmp_x = cos(rphi_)*p1x - sin(rphi_)*p1y
          tmp_y = sin(rphi_)*p1x + cos(rphi_)*p1y
          p1x = tmp_x
          p1y = tmp_y
          ! p2 (muon+)
          tmp_x = cos(rphi_)*p2x - sin(rphi_)*p2y
          tmp_y = sin(rphi_)*p2x + cos(rphi_)*p2y
          p2x = tmp_x
          p2y = tmp_y

          ! fill p1 and p2 in
          event(1)=p1E; event(5)=p2E;
          event(2)=p1x; event(6)=p2x;
          event(3)=p1y; event(7)=p2y;
          event(4)=p1z; event(8)=p2z;
          return
        endif

      enddo

      print*,"[CAUTION] An event was not sampled with the current setup"
      print*,"[CAUTION] of 'ntiral'. So, an event is choosed ramdomly"
      print*,"[CAUTION] on the flat phase space. If you don't like it,"
      print*,"[CAUTION] increase 'ntiral' in this subroutine. "
      print*
      event(1)=0.0d0; event(5)=0.0d0
      event(2)=0.0d0; event(6)=0.0d0
      event(3)=0.0d0; event(7)=0.0d0
      event(4)=0.0d0; event(8)=0.0d0

      ! ----------------------------------------------------------------
      return
      end
************************************************************************
      subroutine rotate2gamma_dir(event,k,uuu,vvv,www)
************************************************************************
      implicit real*8 (a-z)
      real*8 :: event(20)
      kx=k*uuu
      ky=k*vvv
      kz=k*www
      kL=sqrt(kx*kx+ky*ky)

      if(kL .gt. 0.0d0) then
        R11=kx*kz/(k*kL); R12=-ky/kL; R13=kx/k;
        R21=ky*kz/(k*kL); R22= kx/kL; R23=ky/k;
        R31=-kL/k       ; R32= 0.0d0; R33=kz/k;
      else
        R11=1.0d0; R12=0.0d0; R13=0.0d0;
        R21=0.0d0; R22=1.0d0; R23=0.0d0;
        R31=0.0d0; R32=0.0d0; R33=1.0d0;
      endif
      px=event(2)
      py=event(3)
      pz=event(4)
      event(2)= R11*px + R12*py + R13*pz
      event(3)= R21*px + R22*py + R23*pz
      event(4)= R31*px + R32*py + R33*pz
      px=event(6)
      py=event(7)
      pz=event(8)
      event(6)= R11*px + R12*py + R13*pz
      event(7)= R21*px + R22*py + R23*pz
      event(8)= R31*px + R32*py + R33*pz
      return
      end
************************************************************************
      function G2p_el(rt)
************************************************************************
      implicit real*8 (a-h,o-z)
      tau = rt/(4.*938.272046d0**2)
      rFactor = 1.0/(1.0+rt/(0.71*1.0d+6))**4
      G2p_el = rFactor*(1+7.7841*tau)/(1.+tau)
      return
      end
************************************************************************
      function G1p_el(rt)
************************************************************************
      implicit real*8 (a-h,o-z)
      tau = rt/(4.*938.272046d0**2)
      rFactor = 1.0/(1.0+rt/(0.71*1.0d+6))**4
      G1p_el = rFactor*(7.7841*tau)
      return
      end
************************************************************************
      function G2n_el(rt)
************************************************************************
      implicit real*8 (a-h,o-z)
      tau = rt/(4.*938.272046d0**2)
      rFactor = 1.0/(1.0+rt/(0.71*1.0d+6))**4
      G2n_el = rFactor*(3.6481*tau)/(1.+tau)
      return
      end
************************************************************************
      function G1n_el(rt)
************************************************************************
      implicit real*8 (a-h,o-z)
      tau = rt/(4.*938.272046d0**2)
      rFactor = 1.0/(1.0+rt/(0.71*1.0d+6))**4
      G1n_el = rFactor*(3.6481*tau)
      return
      end


************************************************************************
      function W1p_inel(rt,rmf2)
************************************************************************
      implicit real*8 (a-h,o-z)
      parameter ( pi = 3.141592653589793d0 )
      parameter ( alpha = 7.297353d-3 )   ! fine-structure constant
      rmrho  = 775.26d0
      rmrho2 = rmrho**2
      rmp    = 938.272046d0
      rmp2   = rmp**2
      rnu = (rmf2-rmp2+rt)/(2.*rmp)
      rx = rt/(2.*rmp*rnu+rmp2)
      rC = (1.0d-6)*(1.0d-4)/( (0.197**2)*alpha*(pi**2)*8*rmp )

      tmp = (rmrho**4)*(rmf2-rmp2)/((rmrho2+rt)**2)*97.5
      tmp = tmp +250.6*rmp2*((1.-rx)**4)/(1. -1.26*rx +0.96*(rx**2))
      tmp = tmp*rC
      W1p_inel = tmp
      return
      end
************************************************************************
      function W2p_inel(rt,rmf2)
************************************************************************
      implicit real*8 (a-h,o-z)
      parameter ( pi = 3.141592653589793d0 )
      parameter ( alpha = 7.297353d-3 )   ! fine-structure constant
      rmrho  = 775.26d0
      rmrho2 = rmrho**2
      rmp    = 938.272046d0
      rmp2   = rmp**2
      rnu = (rmf2-rmp2+rt)/(2.*rmp)
      rx = rt/(2.*rmp*rnu+rmp2)
      rC = (1.0d-6)*(1.0d-4)/( (0.197**2)*alpha*(pi**2)*8*rmp )

      tmp = rC*( 56.3*(rmf2-rmp2)*rt*rmrho2 )/((rmrho2+rt)**2)
      tmp = tmp * (1.-rt/( 2.*rmp*rnu ))**2
      tmp = tmp + W1p_inel(rt,rmf2)
      tmp = tmp / ( 1.+(rnu**2)/rt )
      W2p_inel = tmp
      return
      end
************************************************************************
      function W2_em(isw,rt,rmf2,iZ,iA)
************************************************************************
      implicit real*8 (a-h,o-z)
      parameter ( rpmass = 938.272046d0 ) ! proton mass [MeV/c**2]
      parameter ( rnmass = 939.565379d0 ) ! neutron mass [MeV/c**2]

      ! --------- coherent
      if ( isw .eq. 1 ) then
        rmi=rpmass*iZ+rnmass*(iA-iZ) !rmi=rpmass*iA
        rd=(1.0d+6)*0.164*iA**(-2./3.)
        W2_em=2.*rmi*(iZ**2)/( 1.+rt/rd )**2

      ! --------- quasi-elastic
      else if ( isw .eq. 2 ) then
        rmi=rpmass
        rPauli=1.0d0
        rQ = sqrt( rt+(rt**2)/(4.*rpmass**2) )
        rPF = 250.0d0 ! [MeV]
        if( rQ .lt. 2.*rPF ) then
          rPauli=(3.*rQ)/(4.*rPF)*( 1.-(1./12.)*(rQ/rPF)**2 )
        end if
        W2_em=2.*rmi*rPauli*( (iA-iZ)*G2n_el(rt) + iZ*G2p_el(rt) )

      ! --------- inelastic
      else if ( isw .eq. 3 ) then ! quasi-elastic
        W2_em=iA*W2p_inel(rt,rmf2)

      else
        print*,"set correct isw"
        call abort ()
      end if
      return
      end
************************************************************************
      function W1_em(isw,rt,rmf2,iZ,iA)
************************************************************************
      implicit real*8 (a-h,o-z)
      parameter ( rpmass = 938.272046d0 ) ! proton mass [MeV/c**2]

      ! --------- coherent
      if ( isw .eq. 1 ) then
        W1_em=0.0d0

      ! --------- quasi-elastic
      else if ( isw .eq. 2 ) then
        rmi=rpmass   ! S.Abe 2019/11/07
        rPauli=1.0d0
        rQ = sqrt( rt+(rt**2)/(4.*rpmass**2) )
        rPF = 250.0d0 ! [MeV]
        if( rQ .lt. 2.*rPF ) then
          rPauli = (3.*rQ)/(4.*rPF)*( 1.-(1./12.)*(rQ/rPF)**2 )
        end if
        W1_em=2.*rmi*rPauli*( (iA-iZ)*G1n_el(rt) + iZ*G1p_el(rt) )

      ! --------- inelastic
      else if ( isw .eq. 3 ) then ! quasi-elastic
        W1_em=iA*W1p_inel(rt,rmf2)

      else
        print*,"set correct isw"
        call abort ()
      end if
      return
      end
************************************************************************
      function randmy(rmin,rmax)
************************************************************************
      implicit real*8 (a-h,o-z)

    ! call random_number(rnd)
    ! randmy = rmin+(rmax-rmin)*rnd

    ! randmy = rmin+(rmax-rmin)*grnd()

      randmy = rmin+(rmax-rmin)*unirn(dummy)

      return
      end

************************************************************************
*                                                                      *
      subroutine sctgmm(eein,wgti,ireg,imat)
*                                                                      *
*        calculate a collision of a photon with an atom.               *
*        last modified by Y.Sakaki on 2019/08/21                       *
*                                                                      *
************************************************************************
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      parameter ( rpmass = 938.27d0, rnmass = 939.58d0 )
      real*8 event(20)

*-----------------------------------------------------------------------

      include 'param00.inc'

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn),  qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)
      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)


      common /kmat1g/ kmat(kvlmax)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)

*-----------------------------------------------------------------------

      common /wparm/  swtm(20), wc1(20), wc2(20) !FURUTA
      common /wparm0/ wc01(20), wc02(20)         !FURUTA
!$OMP THREADPRIVATE(/wparm0/)

      common /kcount/ rncnt(40), rnint(200), rnintr(200),
     &                rnpnt(40), rnpntr(40)

      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)

*-----------------------------------------------------------------------
      common /adjoint/ iadjnt

*-----------------------------------------------------------------------

      parameter ( rmumas = 105.6583715d0 )      ! muon mass [MeV/c**2]

*-----------------------------------------------------------------------

      dimension usave(3), suu(3,2)
      dimension xv(3)

*-----------------------------------------------------------------------

            nclst = -1
            nclsts = 0

*-----------------------------------------------------------------------

            if( mcal .ne. 0 ) return

            nter = 0
            ipt  = 2
            idx  = 0

            erg  = eein
            eg0  = eein
            icl  = ireg
            wgt  = wgti
            wg0  = wgti
            wga  = 0.0d0
            tme  = 0.0d0
            ecmp = 0.0d0

            usave(1) = uuu
            usave(2) = vvv
            usave(3) = www

            do i = 0, 20
               numpal(i) = 0
               rumpal(i) = 0.0d0
            enddo

*-----------------------------------------------------------------------
*     save the incoming direction.
*-----------------------------------------------------------------------

            uold(1) = 0.0
            uold(2) = 0.0
            uold(3) = 1.0

            ntyn = 0
            jsu  = 0

            mk   = imat
            m    = jmd(1+mk)
            m1   = m

*-----------------------------------------------------------------------
*     pair production
*-----------------------------------------------------------------------

            ! --------------------------------------------------------
            ! choose a scattering mode
            izm = iza(m) / 1000       ! Z: proton number
            ims = iza(m) - 1000 * izm ! A: mass number
            igmm_coh   = 1 ! sakaki parameter
            igmm_quasi = 1
            igmm_inel  = 1
            xsc1 = 0.0d0
            xsc2 = 0.0d0
            xsc3 = 0.0d0
            if(igmm_coh   .ne. 0) xsc1 = getxstgmm(1,erg,izm,ims)
            if(igmm_quasi .ne. 0) xsc2 = getxstgmm(2,erg,izm,ims)
            if(igmm_inel  .ne. 0) xsc3 = getxstgmm(3,erg,izm,ims)
            xsc_tot = xsc1 + xsc2 + xsc3

            rnd=randmy(0.0d0,1.0d0)
            if     (rnd .lt. (xsc1     )/xsc_tot) then; isw=1;
            else if(rnd .lt. (xsc1+xsc2)/xsc_tot) then; isw=2;
            else                                      ; isw=3;
            endif

            ! --------------------------------------------------------
            ! generate 4-momenta of the final states
            do i = 1,20; event(i)=0.0d0; enddo;
            call gen_gmm("Born",isw,erg,izm,ims,event)
          ! PHITS rotate the frame automatically... We just need to generate events on z-axis frame.
          ! call rotate2gamma_dir(event,erg,uuu,vvv,www)


            ! set cutoff weight very small to zero
            wc2(6) = 0.0 ! for muon+
            wc2(7) = 0.0 ! for muon-



            ! --------------------------------------------------------
            ! Fill information of the final states

            rmg = rmumas / 1000.0 ! muon mass [GeV]

            iflip=0
            if(unirn(dummy) .lt. 0.5) iflip=1

            do k = 1, 2

              if(k.eq.1) then
                if(iflip .eq. 0) then
                  ityp_ = 7
                  ikfcode_ = 13
                  icharge = -1
                else
                  ityp_ = 6
                  ikfcode_ = -13
                  icharge = +1
                endif
                etot = event(1)/1000. ! GeV
                pxl  = event(2)
                pyl  = event(3)
                pzl  = event(4)
                ekin = (etot - rmg) * 1000.0 ! (MeV)
                xxx1=ekin
              else
                if(iflip .eq. 0) then
                  ityp_ = 6
                  ikfcode_ = -13
                  icharge = +1
                else
                  ityp_ = 7
                  ikfcode_ = 13
                  icharge = -1
                endif
                etot = event(5)/1000. ! GeV
                pxl  = event(6)
                pyl  = event(7)
                pzl  = event(8)
                ekin = (etot - rmg) * 1000.0 ! (MeV)
                xxx2=ekin
              end if

              nclsts = nclsts + 1
              iclusts(nclsts) = 6   ! 4:photon, 6: muon, 7:others(e-,e+)

              jclusts(0,nclsts) = 0 !ipsc
              jclusts(1,nclsts) = 0
              jclusts(2,nclsts) = 0
              jclusts(3,nclsts) = ityp_ ! ityp: 12:e-, 13:e+, 14:photon, 6:muon+, 7:muon-
              jclusts(4,nclsts) = 0
              jclusts(5,nclsts) = icharge
              jclusts(6,nclsts) = 0
              jclusts(7,nclsts) = ikfcode_ ! kf: 22:photon, 11:e-, 13:muon-
              jclusts(8,nclsts) = 0

              qclusts(0,nclsts)  = 0.0
              qclusts(1,nclsts)  = pxl
              qclusts(2,nclsts)  = pyl
              qclusts(3,nclsts)  = pzl
              qclusts(4,nclsts)  = etot  ! total energy (GeV)
              qclusts(5,nclsts)  = rmg   ! rest mass (GeV)
              qclusts(6,nclsts)  = 0.0   ! exitation energy
              qclusts(7,nclsts)  = ekin  ! kinetic energy (MeV)
              qclusts(8,nclsts)  = wgt / wg0
              qclusts(9,nclsts)  = 0.0
              qclusts(10,nclsts) = 0d0 ! x
              qclusts(11,nclsts) = 0d0 ! y
              qclusts(12,nclsts) = 0d0 ! z
            end do

            ! --------------------------------------------------------
      return
      end


















































