************************************************************************
*                                                                      *
      module ion_track_structure
*                                                                      *
*   Track Stracture Code for low energy ions based on                  *
*   M. E. Rudd., Physical Review A, Vol. 38, 6129-6137 (1988)          *
*                                          by T.Ogawa                  *
*                                                                      *
*   ITSART                                                             *
*  (Ion Track Structure model for Arbitrary Radiation and Targets)     *
*                                                                      *
*                                                                      *
************************************************************************

      use ELEDATAMOD
      use GGMARRAYMOD, only : jemi !  gas/condenced switch
      use moddas_material
      implicit double precision(a-h, o-z)

      include 'param-physcnst.inc'
      include 'param.inc'

      parameter (pi = 3.141592653589793d0 )

      parameter (nch_exc = 3 ) ! number of excitation channels

      private t, vf, wc, smj, f1, f2, sdX, 
     & itsreac1, density_calc, iruthreac, nuc_recoil_dataup,
     & its_dataup, itsreac2, Ion_paramset, Exc_paramset, Exc_Xsection,
     & Elastic_Xsec_Collparam
      double precision, private :: rat_reac(2) ! ionization/excitation/elastic reaction branching ratio

      double precision, private :: eIhigh, eps
      double precision, public  :: eIlow
      double precision, private :: xsecscale ! T.Sato 2022/11/30

      double precision, private :: S                ! 4pi a0^2 N (R/eI0)^2
      double precision, private :: eI0              ! ionization potential I0 in the literature
      double precision, private, parameter :: gamma = 2.d0     ! gamma parameter (=2)
      double precision, private, parameter :: R     = 13.6d0   ! binding energy
      double precision, private, parameter :: a0    = 5.29d-11 ! first Bohr orbit radius (m)
      double precision, private, parameter :: w_lbnd= 1.d0     ! produced electron energy lower minimum (eV)
      double precision, private, parameter :: slbnd = 1.d-20   ! electron production cross section lower bound
      double precision, private :: alp,A1,B1,C1,D1,E1,A2,B2,C2,D2      ! parameter in Table 1 in the literature
      double precision, private :: sigma(1:28)
!$OMP THREADPRIVATE(eIlow,eIhigh,eps,S,eI0,xsecscale) ! T.Sato 2022/11/30, add xsecscale
!$OMP THREADPRIVATE(alp,A1,B1,C1,D1,E1,A2,B2,C2,D2)
!$OMP THREADPRIVATE(sigma, rat_reac)
      common / tsminmax  / tsmax

c Parameters of Mueller expansion of Moliere elastic scattering cross section
      double precision, private, parameter :: aMoli0 = -20.45d0
      double precision, private, parameter :: aMoli1 = -71.d0
      double precision, private, parameter :: aMoli2 = 422.097d0
      double precision, private, parameter :: aMoli3 = -1429.7d0
      double precision, private, parameter :: bMoli1 = 7.d-3
      double precision, private, parameter :: bMoli2 = 3.87d-2
      double precision, private, parameter :: bMoli3 = 8.26d-1

      double precision, private, parameter :: emscale =3.d6 ! scaling energy in eV/n, T.Sato 2022/11/22

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /kmat1g/ kmat(kvlmax)

      double precision, allocatable, private :: xsec_brkdwn(:,:)
!$OMP THREADPRIVATE(xsec_brkdwn)
      
      common / ruthang/ rathe1, rathe2 

      integer iflnel ! no electron emission flag. 1:skip electron emission, 0:normal
!$OMP THREADPRIVATE(iflnel)      

! flag to distinguish KURBUC (lflgTS is not used so lflgTS = 0) and this TS mode (lflgTS = 1). 
! Energy recovery is conducted in ecol if lflgTS == 2. lflgTS is doubled when ecol is called from partrs (not doubled when called from talls**). 
      integer lflgTS
!$OMP THREADPRIVATE(lflgTS)

      data lflgTS /0/

      data iflnel /0/

      contains

************************************************************************
*                                                                      *
*                                                                      *
*     |                                                                *
*     |      Track structure common tools                              *
*   \   /                                                              *
*    \ /                                                               *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

c=======================================================================
      subroutine pts_flt2(eng,sig_macro)
*                                                                      *
*     determine track structure mean free path                         *
*                                                                      *
*     input                                                            *
*           eng  : incoming particle energy (eV)                       *
*                                                                      *
*     output                                                           *
*                                                                      *
*     sig_macro : macroscopic cross section (/cm)                      *
*     rat_reac(2) : ionization, excitation, elastic reaction ratio     *
*                   passed via module variable                         *
c=================================================================

      implicit double precision (a-h,o-z)
      common / etsexe  / dexc_ene
      common /eparm/  esmax, esmin, emin(20)
!$OMP THREADPRIVATE(/etsexe/)
      double precision, save :: ensav, sig_macrosav
      integer,          save :: matsav, itypsav
!$OMP THREADPRIVATE(ensav, sig_macrosav, matsav)

      data matsav /0/
      data ensav  /0.d0/
      data sig_macrosav /0.d0/

      double precision Exc_Xsec(nch_exc) ! excitation state energy

      if(matsav   .eq.    mat .and. itypsav  .eq. ityp .and.
     &   ensav/eng .gt. 0.95d0 .and. ensav/eng .lt. 1.05d0) then ! reuse result
          sig_macro = sig_macrosav
          return
      endif

      xlm       = 0.d0
      sig_macro = 0.d0
      eIlow     = 0.d0
      eIhigh    = 0.d0
      rat_reac  = 0.d0
      e_n        = min(eng/dble(ibryf(ityp,ktyp)), emscale) ! X-section above 3 MeV/n follows dedx
      en = e_n*dble(ibryf(ityp,ktyp))

***** Ionization cross section by Rudd formula *********************

      call getsdXsum(e_n,sdXsum,sdXsumdelta,eweighsum,ewid) ! T.Sato 2022/11/26
      sig_macro1 = sdXsum*getdenst(dummy)*1.d28

***** Excitation cross section by Miller-Green formula ************

      call getsig_EX(en, e_n, sig_macro2, eng_sig_EX) ! get energy weighted excitation cross section

***** Elastic cross section by Moliere formula (advanced version of Rutherford) ************
      sig_macro3 = 0.d0
      hydro = denh_das(kmat0+mat)
      if(hydro .gt. 0.d0) then
          call Elastic_Xsec_Collparam(1, 1, eng*1.d-6, xsec, 1.d10,1,1
     &     ,1.d0) ! b is dimensionless collision parameter t
          sig_macro3 = xsec * hydro ! Contribution of 1H
      endif

      lem   = nint( dnel_das(kmat0+mat) )
      do i = 1, lem
        itz  = nint( zz_das(kmat(mat)+i) )
        ita  = nint( a_das(kmat(mat)+i) )
        dens = den_das(kmat(mat)+i)
        call Elastic_Xsec_Collparam(itz, ita, eng*1.d-6, xsec,1.d10,1,1
     &   ,1.d0)
        sig_macro3 = sig_macro3 + xsec * dens   ! rutherford scattering cross section
      enddo
**************************************************************

** T.Sato 2022/12/02, apply correction above 3 MeV/n and/or A>1 (all ions except for proton) ******
      if(eng/dble(mtyp).gt.emscale .or. mtyp.gt.1) then
        eneref=min(emscale*1.d-6,eng/dble(mtyp)*1.d-6) ! reference energy (below 3 MeV/n)
        call atima(1.d0,1.d0,eneref,rmtyp(1,2212),mat,
     &  rng,delt,eccref,5) ! dedx at reference radiation (proton below 3 MeV).   iways = 6 to detour cutoff by estmin
        call atima(dble(mtyp),ctyp,eng*1.d-6,rmtyp(ityp,ktyp),mat,
     &  rng,delt,ecc ,5) ! dedx at current energy
        if(eng/dble(mtyp).gt.emscale) then
          e2nd_emscale=eweighsum/sdXsum  ! mean 2nd electron energy produced by ionization at 3 MeV/n
          call getsdXsum(eng/dble(ibryf(ityp,ktyp)),sdXsum,sdXsumdelta,
     &    eweighsum,ewid) ! get sdXsum and weighsum for acutal projectile energy
          e2nd_mean=eweighsum/sdXsum  ! mean 2nd electron energy produced by ionization at actual projectile energy
          e2ndscale=(sig_macro1*e2nd_emscale + eng_sig_EX) /
     &              (sig_macro1*e2nd_mean    + eng_sig_EX)
         else
          e2ndscale=1.0d0
         endif
        xsecscale=ecc/eccref*e2ndscale ! T.Sato 2022/11/30, scaling factor for ion and exc cross section, will be used in range.f
        sig_macro1 = sig_macro1 * xsecscale !  2022/7/28 Correct cross section to agree with dE/dx
        sig_macro2 = sig_macro2 * xsecscale !  2022/7/28 Correct cross section to agree with dE/dx
      else
        xsecscale=1.0d0
      endif
***************************************************************

      sig_macro = sig_macro1 + sig_macro2 + sig_macro3 ! total cross section of ionization and elastic
      rat_reac(1) = sig_macro2/sig_macro              ! excitation probability
      rat_reac(2) = (sig_macro1+sig_macro2)/sig_macro ! ionization probability

      ensav = eng
      matsav = mat
      itypsav = ityp
      sig_macrosav =  sig_macro
      if(sig_macro .le. 0.d0 .or. eng .le. eIlow) return ! This particle is stopped. To be scored as deposit is "en"

      return
      end subroutine


c=================================================================
      subroutine itsreac(mat1, ierr)
*     mat1 is ntscell ID
*                                                                      *
*     switch electron track structure and Rutherford scattering        *
c=================================================================
      implicit double precision (a-h,o-z)
      rr = unirn(dummy)
      if(rr .lt. rat_reac(1)) then ! excitation
          call itsreac2(mat1, ierr)
          if(ierr .eq. 0) return
      endif
          
      if(rr .ge. rat_reac(2) .or. ierr .eq. 2) then ! elastic. ierr = 2 when reaction rejected owing to energy below threshold
          call iruthreac(ierr)
          if(ierr .eq. 0) return
      endif

      if(rr .lt. rat_reac(2) .or. ierr .eq. 3) then  ! ionization. ierr = 3 when dexe_ene is negative
          call itsreac1(mat1, ierr)
          if(ierr .eq. 0) return
      endif
      
      end subroutine


************************************************************************
*                                                                      *
      function e_mean(en)
*                                                                      *
*       Average energy loss per reaction (eV)                          *
*                                                                      *
*   input  :                                                           *
*                                                                      *
*        en  : projectile energy (eV)                                  *
*                                                                      *
*   output :                                                           *
*                                                                      *
*      e_mean : Mean of particle energy loss (eV)                      *
*                = outgoing electron energy mean + shell energy mean   *
************************************************************************

      implicit double precision (a-h,o-z)
      double precision, save :: ensav, e_meansav,xsecscalesav
      integer,          save :: matsav, itypsav
!$OMP THREADPRIVATE(ensav, e_meansav, xsecscalesav,matsav,itypsav)

      data matsav /0/
      data itypsav /0/
      data ensav  /0.d0/
      data xsecscalesav /0.d0/

      if(matsav .eq. mat .and. e_meansav .ne. 0.d0 .and.
     &  itypsav .eq. ityp    .and. xsecscalesav.eq.xsecscale. and.
     &  ensav/en .gt. 0.95d0 .and. ensav/en .lt. 1.05d0) then ! reuse result
          e_mean = e_meansav
          return
      endif

******* energy loss by ionization from here************************

      e_mean    = 0.d0
      eweighsum = 0.d0
      sdXsum    = 0.d0
      eubnd     = 0.d0
      sigdX0    = 0.d0

************************************************
* ideally eI0 must be switched depending on shell inside sdX calculation. This is approximative
* hope this does not make big difference of mean delta-ray energy
      itz  = nint( zz_das(kmat(mat)+1) ) ! first material
      call shell_eng_highlow(ichem(mat,1),itz)
      eI0     = eIlow  ! use lowest potential
************************************************

      if(en .le. eIlow) return

      e_n = en / dble(ibryf(ityp,ktyp)) ! T.Sato 2022/11/26, en is changed to e_n in this routine

      e_nscale = min(e_n, emscale) ! X-section above 3 MeV/n follows dedx
      enscale = e_n*dble(ibryf(ityp,ktyp))

****** Ioniztaion part ***********
      call getsdXsum(e_nscale,sdXsum,sdXsumdelta,eweighsum,ewid) ! T.Sato 2022/11/26, get cross section at scaling energy (3 MeV/n)
      sig_TS=sdXsum*getdenst(dummy)*1.d28*xsecscale ! sig_TS is macroscopic cross section

      call getsdXsum(e_n,sdXsum,sdXsumdelta,eweighsum,ewid) ! T.Sato 2022/11/26, get sdXsum & ewieghsum at actual energy

****** Excitation part ***********
      call getsig_EX(enscale, e_nscale, sig_EX, eng_sig_EX) ! T.Sato 2022/12/01, get cross section at scaling energy (3MeV/n)
      sig_EX=sig_EX*xsecscale
      eng_sig_EX=eng_sig_EX*xsecscale

******* energy loss by atomic elastic collision by Moliere from here ****************************

       sig_el = 0.d0 ! cross section (b)
       eng_sig_el = 0.d0

       hydro = denh_das(kmat0+mat)
       if(hydro .gt. 0.d0) then
         call Elastic_Xsec_Collparam(1, 1, en*1.d-6, xsec, 1.d10,1,1
     &    ,1.d0)
         sig_el     = xsec * hydro
         call moliere_stopXsec(1, 1, en*1.d-6, e_xsec)
         eng_sig_el = e_xsec * hydro
       endif

       lem   = nint( dnel_das(kmat0+mat) )

       do i = 1, lem
        itchar = nint( zz_das(kmat(mat)+i) )
        itmass = nint( a_das(kmat(mat)+i) )
        call Elastic_Xsec_Collparam(itchar,itmass,en*1.d-6, xsec, 1.d10
     &     ,1,1,1.d0)
        sig_el     = sig_el + xsec * den_das(kmat(mat)+i) ! cross section (b) * density (10^24 atom/cm^3)
        call moliere_stopXsec(itchar, itmass, en*1.d-6, e_xsec)
        eng_sig_el = eng_sig_el + e_xsec * den_das(kmat(mat)+i)
       enddo

       sig_el=sig_el
       eng_sig_el=eng_sig_el

******* energy loss by rutherford to here ******************************

      e_mean   = (eng_sig_el + eng_sig_Ex + eweighsum /sdXsum *sig_TS)
     &          / (sig_el + sig_Ex + sig_TS)

      itypsav   = ityp
      e_meansav = e_mean
      ensav     = en
      matsav    = mat
      xsecscalesav = xsecscale

      return
      end function

c=====================================================================
      subroutine nuc_recoil_dataup(utarg, vtarg, wtarg, ptarg, itchar,
     & itmass, uscat, vscat, wscat, pscat, icycle, weifc)
*                                                                      *
*       record projectile and nuclear recoil final state kinematics    *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       *targ : quantities of recoiled target nucleus                  *
*       *scat : quantities of scattered projectile                     *
*       icycle: cycle id (2 for forced Rutherford angular bias)        *
*       weifc : weight of forced Rutherford scattered particle         *
*                                                                      *
*     output :                                                         *
*                                                                      *
*    jclusts, qclusts                                                  *
*                                                                      *
c====================================================================
      use MMBANKMOD

      implicit double precision (a-h,o-z)

      include 'param00.inc'

      common / clustt / nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common / clustw / jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /pnsave/ egs, uus, vvs, wws, wts, tms, nms, nct(3)
!$OMP THREADPRIVATE(/pnsave/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)


      numpal = 0
      rumpal = 0.d0

      escat  = sqrt(pscat**2 + rmtyp(ityp,ktyp)**2) - rmtyp(ityp,ktyp)

      nclsts = 2 * icycle ! if icycle = 2, angular biased. nclsts is 2 pairs of Recoil + proj (forced and rest)
      i      = 2 * icycle - 1

      iclusts( i ) =  ipatf(ityp,ktyp)
      jclusts(  0, i ) = 0
      jclusts(  1, i ) = ichgf(ityp,ktyp)
      jclusts(  2, i ) = ibryf(ityp,ktyp) - ichgf(ityp,ktyp)
      jclusts(  3, i ) = ityp
      jclusts(  4, i ) = 0
      jclusts(  5, i ) = nint(ctyp)  !  store charge state
      jclusts(  6, i ) = ibryf(ityp,ktyp)
      jclusts(  7, i ) = ktyp
      jclusts(  8, i ) = 0

      qclusts(  0, i ) = 0.d0                      ! impact parameter
      qclusts(  1, i ) = pscat * uscat * 1.d-3       ! (GeV/c)
      qclusts(  2, i ) = pscat * vscat * 1.d-3       ! (GeV/c)
      qclusts(  3, i ) = pscat * wscat * 1.d-3       ! (GeV/c)
      qclusts(  4, i ) = (escat + rmtyp(ityp,ktyp)) * 1.d-3 ! (GeV)
      qclusts(  5, i ) = rmtyp(ityp,ktyp) * 1.d-3  ! (GeV)
      qclusts(  6, i ) = 0.d0                      ! excitation energy (MeV)
      qclusts(  7, i ) = escat                     ! kinetic energy (MeV)
      qclusts(  8, i ) = weifc                     ! weight change
      qclusts(  9, i ) = 0.d0                      ! time delay (ns)
      qclusts( 10, i ) = 0.d0                      ! x displacement (cm)
      qclusts( 11, i ) = 0.d0                      ! y displacement (cm)
      qclusts( 12, i ) = 0.d0                      ! z displacement (cm)

      uus  =     uscat
      vvs  =     vscat
      wws  =     wscat
      egs  =    escat
      wts  =    wt(ibkwt+no,ipomp+1)
      tms  =    tc(ibktc+no,ipomp+1)
      nms  =  name(ibknam+no,ipomp+1)

      if(itchar .eq. 1 .and. itmass .eq. 1) then
        ktyptar = 2212
        ityptar = 1
      else
        ktyptar = itchar * 1000000 + itmass
        ityptar = kftp(ktyptar)
      endif
      weitar  = rmtyp(ityptar,ktyptar)
      etarg  = sqrt(ptarg**2 + weitar**2) - weitar

c // Recoiled target //
      iclusts(     i+1 ) = ipatf(ityptar,ktyptar)
      jclusts(  0, i+1 ) = 0
      jclusts(  1, i+1 ) = itchar
      jclusts(  2, i+1 ) = itmass - itchar
      jclusts(  3, i+1 ) = ityptar
      jclusts(  4, i+1 ) = 0
      jclusts(  5, i+1 ) = itchar
      jclusts(  6, i+1 ) = itmass
      jclusts(  7, i+1 ) = ktyptar

      qclusts(  0, i+1 ) = 0.d0                      ! impact parameter
      qclusts(  1, i+1 ) = ptarg * utarg * 1.d-3     ! (GeV/c)
      qclusts(  2, i+1 ) = ptarg * vtarg * 1.d-3     ! (GeV/c)
      qclusts(  3, i+1 ) = ptarg * wtarg * 1.d-3     ! (GeV/c)
      qclusts(  4, i+1 ) = (etarg + rmtyp(ityptar,ktyptar)) * 1.d-3 ! (GeV)
      qclusts(  5, i+1 ) = rmtyp(ityptar,ktyptar) * 1.d-3      ! (GeV)
      qclusts(  6, i+1 ) = 0.d0                      ! excitation energy (MeV)
      qclusts(  7, i+1 ) = etarg                     ! kinetic energy (MeV)
      qclusts(  8, i+1 ) = weifc                     ! weight change
      qclusts(  9, i+1 ) = 0.d0                      ! (ns)
      qclusts( 10, i+1 ) = 0.d0                      ! (cm)
      qclusts( 11, i+1 ) = 0.d0                      ! (cm)
      qclusts( 12, i+1 ) = 0.d0                      ! (cm)

      if(icycle .eq. 1) then ! Do it only once. 
      numpal(ityptar)       = numpal(ityptar) + 1
      rumpal(ityptar)       = rumpal(ityptar) + wt(ibkwt+no,ipomp+1)
      endif

      return
      end subroutine


************************************************************************
*                                                                      *
      function chrg_eff(e,ityp,ktyp)
*                                                                      *
*       effective charge calculated using ZBL formula                  *
*       Radiation Physics and Chemistry 182 (2021) 109352              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       e     : projectile energy in MeV/n                             *
*       kf    : projectile kf code                                     *
*                                                                      *
*     output :                                                         *
*                                                                      *
*    chrg_eff : effective charge                                       *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      include 'param-physcnst.inc'
      double precision, save :: ensav, chrgsav
      integer,          save :: ktypsav, matsav
      double precision, save :: ensav3, chrgsav3
      integer,          save :: ktypsav3, mat3
!$OMP THREADPRIVATE(ensav,  chrgsav,  ktypsav,  matsav)
!$OMP THREADPRIVATE(ensav3, chrgsav3, ktypsav3, mat3)

      data ensav    /0.d0/
      data chrgsav  /0.d0/
      data ktypsav  /0/
      data matsav   /0/
      data ensav3   /0.d0/
      data chrgsav3 /0.d0/
      data ktypsav3 /0/
      data matsav3  /0/
      
      if(ityp .ge. 1 .and. ityp .le. 14)  then
          chrg_eff = ichgf(ityp,ktyp) ! don't consider charge reduction
          return
      endif

      if(ktypsav .eq. ktyp .and. matsav .eq. mat .and. 
     &       ensav/e .gt. 0.95d0 .and. ensav/e .lt. 1.05d0) then ! reuse result 
          chrg_eff = chrgsav
          return
      elseif(ktypsav3 .eq. ktyp .and. matsav3 .eq. mat .and. 
     &       ensav3/e .gt. 0.95d0 .and. ensav3/e .lt. 1.05d0) then ! reuse result 
          chrg_eff = chrgsav3
          return
      endif

      v0 = 1.d0/physc(4)  ! Bohr velocity in unit of c

      v = Sqrt((e*dble(mtyp))**2 + 2.d0 * rmtyp(ityp,ktyp)*e*dble(mtyp))
     &                        /(e*dble(mtyp) + rmtyp(ityp,ktyp))   ! Speed in unit of c

      del = 0.d0 ! electron density (electron/cm**3)

      lem   = nint( dnel_das(kmat0+mat) )
      do i = 1, lem
        itz  = nint( zz_das(kmat(mat)+i) )
        dens = den_das(kmat(mat)+i)
        del = del + dens * itz
      enddo

      vfe = (3.d0*physc(1)**3*del*1.d24)**(1.d0/3.d0) * physc(3)
     &      * 1.d-13 / (rstms(12)* 1.d3 )        ! Fermi velocity in unit of c

      vr = v * (1.d0 + vfe**2 / (5.d0 * v**2) )

      yr = vr / (v0 * ichgf(ityp,ktyp)**(2.d0/3.d0) )

      q = 1 - Exp(0.803d0   * yr**0.3d0 - 1.3167   * yr**0.6d0 -
     &            0.38157d0 * yr        - 0.008983 * yr**2     )

      lambda= c2 * a0 * (1.d0 - q)**(2.d0/3.d0) /
     & (ichgf(ityp,ktyp)**(1.d0/3.d0) * (1.d0 - (1.d0 - q)/7.d0))

      gam = q + (1.d0 - q)/2.d0 * (v0/vfe)**2 * Log(1.d0 + (c1 *
     & lambda/a0 * vfe/v0 )**2 )

      chrg_eff = max(gam * ichgf(ityp,ktyp), 1.d0) ! effective charge minimum is uncertain.

      ktypsav  = ktyp
      matsav   = mat
      ensav    = e
      chrgsav  = chrg_eff
      
      if(e_n .gt. 2.99d6 .and. e_n .lt. 3.01d6) then ! store result at normalization energy
          ktypsav3  = ktyp
          matsav3   = mat
          ensav3    = e
          chrgsav3  = chrg_eff
      endif
      
      return
      end function

************************************************************************
*                                                                      *
      function density_calc(ichemID)
*                                                                      *
*       calculate molecular density                                    *
*                                                                      *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*      ichemID : Chemical compound ID                                  *
*                                                                      *
*     output :                                                         *
*                                                                      *
*   density_calc : density in 10^24 (molecule/cm^3)                    *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      denst = denh_das(kmat0+mat)
      lem   = nint( dnel_das(kmat0+mat) )
      do i = 1, lem
          denst = denst + den_das(kmat(mat)+i)
      enddo

      Select Case (ichemID)
        Case(1)       ! H2O
          density_calc = denst/3.d0
        Case(2)       ! CO2
          density_calc = denst/3.d0
        Case(3)       ! NH2
          density_calc = denst/3.d0
        Case(4)       ! NH3
          density_calc = denst/4.d0
        Case(5)       ! SF6
          density_calc = denst/7.d0
        Case(6)       ! TeF6
          density_calc = denst/7.d0

        Case(10001)       ! CH4
          density_calc = denst/5.d0
        Case(10002)       ! CH3
          density_calc = denst/4.d0
        Case(10003)       ! C2H2
          density_calc = denst/4.d0
        Case(10004)       ! C2H4
          density_calc = denst/6.d0
        Case(10005)       ! C2H6
          density_calc = denst/8.d0
        Case(10006)       ! C6H6
          density_calc = denst/12.d0
        Case(10007)       ! CH3-2NH
          density_calc = denst/8.d0

        Case(20010)       ! H2
          density_calc = denst/2.d0
        Case(20070)       ! N2
          density_calc = denst/2.d0
        Case(20080)       ! O2
          density_calc = denst/2.d0

        Case default
          density_calc = denst
      end select

      return
      end function


************************************************************************
*                                                                      *
      function getdenst(dummy)
*                                                                      *
*     get density for calculating macroscopic cross section            *
*     output : getdenst (10^24 atoms/cm^3 or molecules/cm^3)           *
*                                                                      *
************************************************************************
      common /kmat1a/ mxmat, mxmat0, mxnel

      double precision, allocatable,save ::  getdenstmem(:) ! remember getdenst calculated already
      logical, allocatable,save          :: lgetdenstmem(:) ! flag of getdenst calculated already
      integer ,save                      :: lfirst ! flag of getdenst calculated already
      data  lfirst /1/
      
!$OMP CRITICAL (getdenst0_crit)
      if(lfirst .eq. 1) then ! array allocation for the first time
        allocate( getdenstmem(mxmat))
        allocate(lgetdenstmem(mxmat))
        getdenstmem  = 0.d0
        lgetdenstmem = .false. 
        lfirst = 0
      endif
!$OMP END CRITICAL (getdenst0_crit)

      if(lgetdenstmem(mat)) then ! use calculated data if exists
        getdenst = getdenstmem(mat)
        return
      endif
      
      getdenstmem(mat) = 0.d0
      
***** Atomic density ************
      hydro = denh_das(kmat0+mat)
      denst = 0.d0
      if(hydro .gt. 0.d0) denst = hydro ! Contribution of 1H
      lem = nint( dnel_das(kmat0+mat) )
      do i = 1, lem
        dens = den_das(kmat(mat)+i)
        denst = denst + dens 
      enddo
**********************************      
      
      if(frac(mat,1) .ne. 0.d0) then ! chem specified
          i          = 1
          atomnumb   = 0.d0 ! (number of constituent atoms per molecule) averaged over chemicals
          do while(frac(mat,i) .gt. 0.d0)
            atomnumb = atomnumb + denst/density_calc(ichem(mat,i)) 
     &                  * frac(mat,i) 
            i = i + 1
          enddo
          getdenst = denst / atomnumb

      else ! chem unspecified. Calculated as mixture of pure elements

          getdenst = denst
          
      endif

!$OMP CRITICAL (getdenst1_crit)
      getdenstmem(mat) = getdenst
      lgetdenstmem(mat) = .true.
!$OMP END CRITICAL (getdenst1_crit)

      return
      
      end function

************************************************************************
*                                                                      *
*    / \                                                               *
*   /   \                                                              *
*     |      Common tools                                              *
*     |                                                                *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*     |                                                                *
*     |      Ionization part                                           *
*   \   /                                                              *
*    \ /                                                               *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

c=================================================================
      subroutine itsreac1(mat1, ierr)
*                                                                      *
*     determine kinematics of ionization reaction                      *
*                                                                      *
*     input                                                            *
*                                                                      *
*      mat1    : track structure mode switch                           *
*                                                                      *
*                                                                      *
*     output                                                           *
*                                                                      *
*      outgoing electron/proton energy/momentum sent through jclusts   *
*      ierr : error identifier                                         *
*                                                                      *
*                                                                      *
*    Note that charge state is not considered because it scales        *
*    cross sections uniformly.                                         *
*                                                                      *
c=================================================================
c=================================================================
      use MMBANKMOD ! , only : x,y,z,e,ec,nty,nkf ! commented out because there might be others

      implicit double precision (a-h,o-z)

      include 'param00.inc'
      include 'param.inc'
      include 'param-physcnst.inc'

      parameter (pi=3.141592654d0)

      common / ele2nd  / e2nd,u2nd,v2nd,w2nd
!$OMP THREADPRIVATE(/ele2nd/)
      common / etsexe  / dexc_ene
!$OMP THREADPRIVATE(/etsexe/)
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)
      common /kmat1g/ kmat(kvlmax)

      dimension vel(3)

      double precision, save::e_nsav, sdXsum, ewid,sdXsumdelta,eweighsum
      integer,          save::matsav, itypsav
!$OMP THREADPRIVATE(e_nsav, sdXsum, ewid , sdXsumdelta, eweighsum)
!$OMP THREADPRIVATE(matsav, itypsav)

      data itypsav /0/
      data matsav /0/
      data e_nsav  /0.d0/
      data sdXsum /0.d0/
      data sdXsumdelta /0.d0/
      data eweighsum /0.d0/
      data ewid   /0.d0/

      ierr = 0
      ec(ibkec+no,ipomp+1) = ec(ibkec+no,ipomp+1) * 1.d+6   ! from (MeV) to (eV)

      en  = ec(ibkec+no,ipomp+1) ! incident energy
      it = nty(ibknty+no,ipomp+1)
      kf = nkf(ibknkf+no,ipomp+1)

      e_n = en / dble(ibryf(ityp,ktyp)) ! T.Sato 2022/11/26, en is changed to e_n in this routine

      e2nd      = 0.d0
      eatom     = 0.d0
      
      if(iflnel .eq. 1 ) then ! no emission event. determined by restricted LET energy conservation
          pejc = sqrt( en**2 + 2.d0 * rmtyp(it,kf) * 1.d+6 * en)
          call itsdataup(0.d0,0.d0,0.d0,1.d0,pejc,en,0.d0,0.d0,0.d0,
     &                   0.d0,0.d0,0.d0,0,0)
          iflnel = 0
          ec(ibkec+no,ipomp+1) = en * 1.d-6 ! from (eV) to (MeV)
          dexc_ene     = e(ibke+no,ipomp+1) - ec(ibkec+no,ipomp+1)
          return
      endif

      if(matsav   .ne.    mat .or. itypsav  .ne. ityp .or.
     &   e_nsav/e_n .lt. 0.95d0 .or. e_nsav/e_n .gt. 1.05d0) then ! not reuse result

        matsav = mat
        e_nsav  = e_n
        itypsav = ityp


        call getsdXsum(e_n,sdXsum,sdXsumdelta,eweighsum,ewid)

      endif

	if( sdXsumdelta .lt. sdXsum * unirn(dummy) ) then ! TS model or Delta-ray model      sdXran  = sdXsum * unirn(dummy)

      eng   = w_lbnd      
      sdXran  = (sdXsum - sdXsumdelta) * unirn(dummy)

      do ii = 0, 999

        call sdX(eng,e_nsav,1,xsec,ebmean) ! 2024/4/25 use e_nsav because sdXsum was calculated for e_nsav. If en is used instead, sum(xsec) does not agree as en goes down and xsec goes up. 
        sdXran  = sdXran - xsec * eng * (ewid - 1.d0)
        if(sdXran .lt. 0.d0) goto 90

        eng     = eng * ewid
      enddo

  90  e2nd     = eng + eng * (ewid - 1.d0) * unirn(dummy)

      else
          
        emind = 1.d3
        rn = unirn(dummy)
        e2nd = 1.d0/(1.d0 / emind * (1.d0 - rn) + 1.d0 /emaxd(e_n) * rn) ! Delta-ray energy sampling taken from subroutine delprd
        call sdX(1.d3,e_nsav,1,xsec,ebmean) ! 2024/6/24 call sdX to calculate shell-wise partial cross section (xsec_brkdwn) used to determine shell. Its "dkin" is used in angular distribution

      endif    

!<-20211126murofushi update
!--      sdXtot   = sum(xsec_brkdwn(1:28,1:nint(dnel(kmat(mat)+1)))) 
!--     &           * unirn(dummy)
      sdXtot   = sum(xsec_brkdwn(1:28,0:nint(dnel_das(kmat0+mat)))) 
     &           * unirn(dummy)

      if(sdXtot .le. 0.d0) then
          ierr = 1
          ec(ibkec+no,ipomp+1) = en * 1.d-6 ! from (eV) to (MeV)
          return
      endif

      lng_rel = 0
      dkin    = 0.d0 ! Kinetic energy of recoiled electron in the orbit. This migh be zero if e- is from outermost shell
      j = 0
      if(frac(mat,1) .le. 0.d0) then ! unnecessary anymore? 
        do j = 1, 28
         do i = 0, nint( dnel_das(kmat0+mat) )
          sdXtot = sdXtot - xsec_brkdwn(j,i)
          if(sdXtot .le. 0.d0) goto 100
         enddo
        enddo
100     if(i .eq. 0) then
          eI0 = pot_elem(j,1)
        else
          eI0 = pot_elem(j,nint( zz_das(kmat(mat)+i) ))
          If( nint(zz_das(kmat(mat)+i)) .le. 5 ) goto 110 ! No relaxation data below Z = 5
          do jj = 28, 1, -1 ! atomic relaxation 
            if(pot_elem(jj,nint(zz_das(kmat(mat)+i))) .lt. eI0 .and.
     &         pot_elem(jj,nint(zz_das(kmat(mat)+i))) .ne. 0) then
              call do_atom_relax(nint(zz_das(kmat(mat)+i)),j)
              exit
            endif
          enddo
        endif
      else
        i = 1
        do while(frac(mat,i) .gt. 0.d0)
         if(frac(mat,i) .le. 0.d0) goto 110
         do j = 1, 28
           sdXtot = sdXtot - xsec_brkdwn(j,i)
           if(sdXtot .le. 0.d0) then
            eI0 = pot_mol(j,ichem(mat,i))

**** Find the ionized atom ******************************************
            edif = 1.d100
            do ii = 1, nint( dnel_das(kmat0+mat) ) 
             do jj = 1, 28
              If( abs(pot_elem(jj,nint(zz_das(kmat(mat)+ii)))-eI0) 
     &         .lt. edif) then
               If(pot_elem(jj,nint(zz_das(kmat(mat)+ii))) .eq. 0.d0) 
     &          cycle
               edif=abs(pot_elem(jj,nint(zz_das(kmat(mat)+ii)))-eI0)
               ir = ii
               jr = jj ! find closest level in molecule
              endif
             enddo
            enddo
            do jj = 1, 28
             If(pot_elem(jj,nint(zz_das(kmat(mat)+ir))).lt.eI0 .and.
     &          pot_elem(jj,nint(zz_das(kmat(mat)+ir))) .ne. 0) then
              call do_atom_relax(nint(zz_das(kmat(mat)+ir)),jr)
              exit
             endif
            enddo
*********************************************************************

            goto 110
           endif
         enddo
         i = i + 1
        enddo
      endif

 110  idelem = i
      idshel = j
      dkin = dkin_elem(j,nint(zz_das(max(kmat(mat)+i,1))))
      
      pel  = sqrt(e2nd**2 + 2.d0 * rmtyp(12,11) * 1.d+6 * e2nd) ! electron momentum
      
      theelc = elec_ioniz_angle(e_n*1.d-6, it, kf, e2nd*1.d-6, 
     & dkin*1.d-6) ! electron emission angle (radian)
      pel_l  = pel * cos(theelc) ! electron longitudinal momentum 
      pel_t  = pel * sin(theelc) ! electron transverse momentum 
      pinc   = sqrt(  en**2 + 2.d0 * rmtyp(it,kf) * 1.d+6 * en)   ! projectile momentum

      deatom = 0.d0
c      prec = (-(en - eI0 - e2nd) * rmtyp(it,kf) * 2.d+6 + pinc**2)
c     &/(2.d0 * pinc)
c      pejc_l = pinc - prec  ! projectile longitudinal momentum. 
c      pejc_t =        prec  ! projectile transverse momentum. dirfac is just a fitting parameter
      eejc = en - eI0 - e2nd ! target recoil is small. Disregarded 
      pejc = sqrt(eejc**2 + 2.d+6 * rmtyp(it,kf) * eejc)
      pejcsin = 1.d-1**(2.5d0 + 6.d0 * unirn(dummy)) ! angle factor
      pejc_l = pejc * (1.d0 - pejcsin**2)  ! projectile longitudinal momentum. pejcsin is fitting parameter
      pejc_t = pejc * pejcsin  ! projectile transverse momentum. 

      do i = 1, 100
       pejc= sqrt(pejc_t**2 + pejc_l**2 + rmtyp(it,kf) * 2.d+6 * deatom)
       enp = sqrt(pejc**2 + (rmtyp(it,kf) * 1.d+6)**2) - rmtyp(it,kf) * 
     & 1.d+6
       eatom = en - eI0 - e2nd - enp
       if(abs(eatom) .lt. en * 1.d-13 ) then ! when en and enp are equal, noise dominates eatom
         eatom = 0.d0  
         exit
       elseif(eatom .gt. 0.d0 .and. eatom .lt. unirn(dummy) * 5.d0) then
         exit
       endif    
       deatom = deatom + eatom
      enddo
      
      eatom = max(eatom,0.d0)
        
      theejc = atan(pejc_t/max(pejc_l,1.d-10))

      if(enp .le. 0.d0) then ! No reaction at the end of range
          ierr = 1
          ec(ibkec+no,ipomp+1) = en * 1.d-6 ! from (eV) to (MeV)
          return
      endif

      if(i .ge. 1000) then ! No reaction at the end of range
          write(*,*) 'itsmode warning', en, eI0, ichemID, sdXtot
          ierr = 1
          ec(ibkec+no,ipomp+1) = en * 1.d-6 ! from (eV) to (MeV)
          return
      endif

      phielc = 2.d0 * unirn(dummy) * pi

      call vel_vec(vel,theelc,phielc)

      u2nd = vel(1)
      v2nd = vel(2)
      w2nd = vel(3)

      phiejc = 2.d0 * unirn(dummy) * pi

      call vel_vec(vel,theejc,phiejc)

      uejc = vel(1)
      vejc = vel(2)
      wejc = vel(3)

! Reaction score
      atmrc(10,2) = atmrc(10,2) + 1.d0
*-----------------------------------------------------------------------
*        data up in bank
*-----------------------------------------------------------------------

      call itsdataup(xsec_brkdwn(idshel,idelem), uejc, vejc, wejc, pejc,
     &enp, u2nd, v2nd, w2nd, pel, e2nd, eatom, 
     & nint(zz_das(max(kmat(mat)+idelem,1))),
     & nint( a_das(max(kmat(mat)+idelem,1))))

      ec(ibkec+no,ipomp+1) = enp * 1.d-6 ! from (eV) to (MeV)

      dexc_ene =
     & e(ibke+no,ipomp+1) - ec(ibkec+no,ipomp+1) - e2nd * 1.d-6 - eatom ! energy deposit to be scored by T-deposit = Continuous energy loss + electron binding energy
     &  - sum(eng_rel(1:lng_rel))
      lflgTS = 1


      end subroutine

c=================================================================
      function elec_ioniz_angle(ep, it, kf, ee, ekin_elem)
*                                                                      *
*     determine angle of electron from ionization reactions            *
*                                                                      *
*     input                                                            *
*                                                                      *
*      ep : Projectile energy (MeV/n)                                  *
*      it : Projectile ityp                                            *
*      kf : Projectile kfcode                                          *
*      ee : Ejected electron energy (MeV)                              *
*      kin_elem : kinetic energy of electron in orbit (MeV)            * 
*                                                                      *
*                                                                      *
*     output                                                           *
*                                                                      *
*      elec_ioniz_angle : electron emission angle(radian)              *
*                                                                      *
*                                                                      *
*    Note that charge state is not considered because it scales        *
*    cross sections uniformly.                                         *
*                                                                      *
c=================================================================
c=================================================================
      implicit double precision (a-h,o-z)
      include 'param-physcnst.inc'
      
      double precision prob(180) 
      
! smoothly switch by outgoing electron energy. 
! e2nd < 10  keV -> 1000% Ogawa original formula, 
! e2nd > 100 keV  ->   0%, Not Ogawa formula. Use Niita's formula of delta-ray
      if(Log10(ee*1.d2) .lt. unirn(dummy)) then 
      
        rel = sqrt(ekin_elem**2 + 2.d0 * ekin_elem * rstms(12)*1.d3)
     & /sqrt(ee**2 + 2.d0 * ee * rstms(12)*1.d3)
        ael = (rel + rel**2/2.d0)/(2.d0+Log10(ep))
        bel =acos(min(1.d0,0.5*sqrt(rstms(1)*ee/(ep*rstms(12)))))
        cel = 0.0017/ep**1.27d0*(ee*1.d6)**(0.482 + 0.278*Log10(ep))
        del = 0.1d0 ! Fermi distribution switching angle (radian)

c Sampling from angular distribution      
        probsum = 0.d0
        do iel = 1, 180
          ang = (dble(iel)-0.5d0) * physc(1) / 180.d0 
          arg1 = min(cel*(ang-del)*180.d0/physc(1),708.d0)
          prob(iel) = ael/((ang - bel)**2 + ael**2)
     &       /(exp(arg1) + 1.d0) * sin(ang) ! *sin(ang) because function is /sr
          probsum = probsum + prob(iel)
        enddo
      
        probran = probsum * unirn(dummy)
        do iel = 1, 180
          probran = probran - prob(iel)
          if(probran .le. 0.d0) exit
        enddo
        elec_ioniz_angle = (dble(iel)- unirn(dummy)) * physc(1) / 180.d0 
      
      else
      
        rmase = rstms(12)*1.d3  
        rmass = rmtyp(it,kf)

        efin = ep * mtyp - ee 
        pprj = sqrt( efin**2 + 2.0 * rmass * efin )
        pdel = sqrt( ee**2 + 2.0 * rmase * ee )

        etot = efin + rmass - rmase
         

        elec_ioniz_angle = acos(min(ee * etot / pprj / pdel,1.d0))
          
      endif    
          
      end function
      
************************************************************************
*                                                                      *
      subroutine getsdXsum(e_n, sdXsum, sdXsumdelta, eweighsum,ewid)
*                                                                      *
*     get integral of single scattering cross section, Eq.(6) of Rudd  *
*     microscopic cross section is returned in different from getsig_EX *
*   input  :                                                           *
*        e_n : projectile energy (eV/n)                                *
*                                                                      *
*   output :                                                           *
*    ewid        : step size of energy (used for energy determination) *
*    sdXsum      : integral cross section (m^2)                        *
*                  integrated over full energy range                   *
*    sdXsumdelta : integral cross section (m^2)                        *
*                  integrated above 1 keV using delta-ray formula      *
*    eweighsum   : energy-weighted integral cross section (eV.m^2)     *
*                                                                      *
*         Attention! To get energy-weighted                            *
*         Energy = electron kinetic energy + biding energy             *
*                  Atom recoil is disregarded                          *
*                                                                      *
************************************************************************
      implicit real*8(a-h,o-z)

      double precision, save :: e_nsav, sdXsumsav, sdXsumdeltasav,
     & eweighsumsav, ewidsav
      double precision, save :: e_nsav3,sdXsumsav3,sdXsumdeltasav3,
     & eweighsumsav3,ewidsav3
      integer,          save :: matsav
!$OMP THREADPRIVATE(matsav,e_nsav,sdXsumsav,sdXsumdeltasav)
!$OMP THREADPRIVATE(eweighsumsav,ewidsav)
      
      data matsav /0/
      data e_nsav  /0.d0/

      data sdXsumsav /0.d0/
      data sdXsumdeltasav /0.d0/
      data eweighsumsav /0.d0/
      data ewidsav /0.d0/

      data sdXsumsav3 /0.d0/
      data sdXsumdeltasav3 /0.d0/
      data eweighsumsav3 /0.d0/
      data ewidsav3 /0.d0/

      if(matsav .eq. mat .and. 
     &  e_nsav/e_n .gt. 0.95d0 .and. e_nsav/e_n .lt. 1.05d0) then ! reuse previous result 
          sdXsum=sdXsumsav
          sdXsumdelta = sdXsumdeltasav
          eweighsum = eweighsumsav
          ewid=ewidsav
          return
      elseif(matsav .eq. mat .and. sdXsumsav3 .ne. 0.d0 .and. 
     &  e_n .gt. 2.99d6 .and. e_n .lt. 3.01d6) then ! reuse result at normalization energy
          sdXsum=sdXsumsav3
          sdXsumdelta = sdXsumdeltasav3
          eweighsum = eweighsumsav3
          ewid=ewidsav3
          return
      endif

! determine sdXref0, cross section reference point
      sdXref0   = 0.d0
      
      eIlowsav = 1.d10
      eIhighsav = 0.d0
      if(frac(mat,1) .gt. 0.d0) then
       i = 1
       do while ( frac(mat,i) .gt. 0.d0 )  
        call shell_eng_highlow(ichem(mat,i), 0)
        eIlowsav  = min(eIlow,eIlowsav)
        eIhighsav = max(eIlow,eIhighsav)
        i = i + 1
       enddo
      else
       do i = 1, nint( dnel_das(kmat0+mat) )  
        call shell_eng_highlow(-1, nint( zz_das(kmat(mat)+i) ))
        eIlowsav  = min(eIlow,eIlowsav)
        eIhighsav = max(eIlow,eIhighsav)
       enddo
      endif
      
      eI0 = eIlow  ! w_v_ubnd is en dependent (insensitive to ichemID )
      call sdX(1.d-2,e_n,2,xsec,ebmean)
      sdXref0 = xsec ! chem unspecified
      
! check energy upper bound.       
      w_v_ubnd  = 0.d0
      xsecsum   = 1.d100
      w_v_ubnd_max  = 1.d3/ vf(e_n)/ eI0 ! upper bound corresponding to 1 keV 
      do while (xsecsum .ge. slbnd * sdXref0 .and. 
     &          w_v_ubnd .lt. w_v_ubnd_max) 
         w_v_ubnd = w_v_ubnd + 1.d1
         xsecsum   = 0.d0
         ! determine w_v_ubnd
         call sdX(w_v_ubnd,e_n,2,xsec,ebmean)  
         xsecsum = xsecsum + xsec
      enddo
      
      eI0   = eIhigh
      eubnd = min(w_v_ubnd * vf(e_n) * eI0, 1.d3)
      
      ewid  = (eubnd/w_lbnd)**0.001d0
      eng   = w_lbnd
      sdXsum = 0.d0
      eweighsum = 0.d0

! calculate energy weighted sum below 1 keV
      do ii = 1, 1000
      
        sigdX   = 0.d0
        ebmean1 = 0.d0
        call sdX(eng,e_n,1,xsec,ebmean)  
        sigdX   = sigdX + xsec
        ebmean1 = ebmean1 + ebmean * xsec
        
        eweighsum = eweighsum + 
     &        (eng * sigdX + ebmean1) * eng * (ewid - 1.d0) ! "eng * (ewid - 1.d0)" is energy bin integration width
        sdXsum    = sdXsum    + sigdX * eng * (ewid - 1.d0) ! arbitrary unit. Normalized later
        
        if(ii .eq. 1) sigdX0 = sigdX  
        if(sigdX .lt. sigdX0 * 1.d-20) exit ! CPU time saving
      
        eng     = eng * ewid
      enddo
      
      
! calculate energy weighted sum above 1 keV or energy upper bound of Rudd model
      if(emaxd(e_n) .gt. eubnd) then 
          sdXsumdelta = sigdX * eubnd**2 * ewid * 
     &                                    (1.d0/eubnd -1.d0/emaxd(e_n)) ! smoothly joint X-section from Rudd model
          sdXsum = sdXsum + sdXsumdelta 
          eweighsum = eweighsum + sigdX * eubnd**2 * ewid * 
     &                                          log(emaxd(e_n) / eubnd)
      endif
      
      sdXsumsav = sdXsum
      eweighsumsav = eweighsum
      sdXsumdeltasav = sdXsumdelta
      ewidsav = ewid
      e_nsav = e_n
      matsav = mat

      if(e_n .gt. 2.99d6 .and. e_n .lt. 3.01d6) then ! store result at normalization energy
          sdXsum3 = sdXsumsav
          sdXsumdelta3 = sdXsumdeltasav
          eweighsum3 = eweighsumsav
          ewid3 = ewidsav
          matsav = mat
      endif

      return

      end subroutine
      
c================================================
      subroutine itsdataup(xsec, uejc, vejc, wejc, pejc,
     & enp, u2nd, v2nd, w2nd, pel, e2nd, eatom, iz, ia)
c===============================================
      use MMBANKMOD

      implicit double precision (a-h,o-z)

      include 'param00.inc'

      common / clustt / nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common / clustw / jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /pnsave/ egs, uus, vvs, wws, wts, tms, nms, nct(3)
!$OMP THREADPRIVATE(/pnsave/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)

! hirata 2024-02-08 Counting generated electron with eps value
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)
      common /etsart/ bgets(kvlmax), ebgets(kvlmax),
     &                wvets(kvlmax), ewvets(kvlmax)
      common /ets_wvalue/ wvlin
!$OMP THREADPRIVATE(/ets_wvalue/)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common / etsminmax / etsmin, etsmax
! hirata 2024-02-08
      ehnum  = 0.d0
      wvlin = ewvets( idgr(iblz(ibkblz+no,ipomp+1)) )

      numpal = 0
      rumpal = 0.d0
      nclsts = 1

      iclusts( 1 ) =  ipatf(ityp,ktyp)
      jclusts(  0, 1 ) = 0
      jclusts(  1, 1 ) = ichgf(ityp,ktyp)
      jclusts(  2, 1 ) = ibryf(ityp,ktyp) - ichgf(ityp,ktyp)
      jclusts(  3, 1 ) = ityp
      jclusts(  4, 1 ) = 0
      jclusts(  5, 1 ) = nint(ctyp)  !  store charge state
      jclusts(  6, 1 ) = ibryf(ityp,ktyp)
      jclusts(  7, 1 ) = ktyp
      jclusts(  8, 1 ) = 0


      qclusts(  0, 1 ) = 0.d0                      ! impact parameter
      qclusts(  1, 1 ) = pejc * uejc * 1.d-9       ! (GeV/c)
      qclusts(  2, 1 ) = pejc * vejc * 1.d-9       ! (GeV/c)
      qclusts(  3, 1 ) = pejc * wejc * 1.d-9       ! (GeV/c)
      qclusts(  4, 1 ) = enp * 1.d-9 + rmtyp(ityp,ktyp) * 1.d-3 ! (GeV)
      qclusts(  5, 1 ) = rmtyp(ityp,ktyp) * 1.d-3  ! (GeV)
      qclusts(  6, 1 ) = 0.d0                      ! from (eV) to (MeV)
      qclusts(  7, 1 ) = enp * 1.d-6               ! from (eV) to (MeV)
      qclusts(  8, 1 ) = 1.d0                      ! weight change
      qclusts(  9, 1 ) = 0.d0                      ! (ns)
      qclusts( 10, 1 ) = 0.d0                      ! x displacement (cm)
      qclusts( 11, 1 ) = 0.d0                      ! y displacement (cm)
      qclusts( 12, 1 ) = 0.d0                      ! z displacement (cm)

      uus  =     uejc
      vvs  =     vejc
      wws  =     wejc
      egs  =    enp * 1.d-6
      wts  =    wt(ibkwt+no,ipomp+1)
      tms  =    tc(ibktc+no,ipomp+1)
      nms  =  name(ibknam+no,ipomp+1)

      if(iflnel .eq. 1) then
          return
      endif


! 2024-02-08 epsilon value correction for calculating charged particles
      if(e2nd .le. etsmin * 1.d+6)then
         eleemit = 1.d0
         if(wvlin.gt.0.d0) then
           if(e2nd.gt.wvlin) then
            eleemit = dble(int(e2nd/wvlin))
            if(eleemit.lt.1.d0) eleemit = 1.d0
           endif
           ehnum = ehnum + eleemit -1.d0 ! -1.d0 : cancel atmrc(10,2) in itsreac1
           e2nd = 0.d0
         endif
      endif


      ! interaction radius
      nclsts = nclsts + 1
      radis = 0.d0 ! remote hit for ion-ionduced ion production rad * sqrt( unirn(dummy) )
      costh = 2.d0 * unirn(dummy) - 1.d0
      sinth = sqrt(1.d0 - costh**2)
      phi   = 2.d0 * pi * unirn(dummy)

c // Ionized electrons //
      iclusts(     nclsts ) = 7
      jclusts(  0, nclsts ) = 0
      jclusts(  1, nclsts ) = 0
      jclusts(  2, nclsts ) = 0
      jclusts(  3, nclsts ) = 12
      jclusts(  4, nclsts ) = 0
      jclusts(  5, nclsts ) = -1
      jclusts(  6, nclsts ) = 0
      jclusts(  7, nclsts ) = 11

      qclusts(  0, nclsts ) = 0.d0                      ! impact parameter
      qclusts(  1, nclsts ) = pel * u2nd * 1.d-9        ! (GeV/c)
      qclusts(  2, nclsts ) = pel * v2nd * 1.d-9        ! (GeV/c)
      qclusts(  3, nclsts ) = pel * w2nd * 1.d-9        ! (GeV/c)
      qclusts(  4, nclsts ) = e2nd * 1.d-9 + rmtyp(12,11) * 1.d-3 ! (GeV)
      qclusts(  5, nclsts ) = rmtyp(12,11) * 1.d-3      ! (GeV)
      qclusts(  6, nclsts ) = 0.d0                      ! from (eV) to (MeV)
      qclusts(  7, nclsts ) = e2nd * 1.d-6              ! from (eV) to (MeV)
      qclusts(  8, nclsts ) = 1.d0                      ! weight change
      qclusts(  9, nclsts ) = 0.d0                      ! (ns)
      qclusts( 10, nclsts ) = radis * sinth * cos(phi)  ! (cm)
      qclusts( 11, nclsts ) = radis * sinth * sin(phi)  ! (cm)
      qclusts( 12, nclsts ) = radis * costh             ! (cm)

      numpal(12) = numpal(12) + 1
      rumpal(12) = rumpal(12) + wt(ibkwt+no,ipomp+1)

      if(eatom .gt. 0.d0)then ! target atom recoil
       nclsts = nclsts + 1

        if( ia .eq. 1) then ! recoil proton
         iclusts( nclsts )    =  1
         jclusts( 3, nclsts ) =  1
         jclusts( 7, nclsts ) = 2212
         numpal(1)       = numpal(1) + 1
         rumpal(1)       = rumpal(1) + wt(ibkwt+no,ipomp+1)
        else                                    ! recoil heavier atoms
         iclusts( nclsts )    =  0
         jclusts( 3, nclsts ) = 19
!<-20211126murofushi update
!--         jclusts( 7, 3 ) = nint( zz(kmat(mat)+3) )*1000000 
!--     &                   + nint( a(kmat(mat)+4) )
         jclusts( 7, nclsts ) = iz*1000000 + ia
         numpal(19)      = numpal(19) + 1
         rumpal(19)      = rumpal(19) + wt(ibkwt+no,ipomp+1)
        endif

       jclusts(  0, nclsts ) = 0
       jclusts(  1, nclsts ) = iz
       jclusts(  2, nclsts ) = ia - iz
       jclusts(  4, nclsts ) = 0
       jclusts(  5, nclsts ) = iz
       jclusts(  6, nclsts ) = ia
       jclusts(  8, nclsts ) = 0

       rmass = rmtyp(ityp,jclusts(7,3)) * 1.d-3 ! GeV
c       eatom = (sqrt((patom*1.d-9)**2 + rmass**2)-rmass)*1.d3
       eatom = eatom * 1.d-6 ! (eV to MeV)
       patom = sqrt(eatom**2 + 2.d0 * eatom * rmass* 1.d3)* 1.d-3 ! MeV to GeV 

       qclusts(  0, nclsts ) = 0.d0                      ! impact parameter
       qclusts(  1, nclsts ) = 0.d0                      ! (GeV/c)
       qclusts(  2, nclsts ) = 0.d0                      ! (GeV/c)
       qclusts(  3, nclsts ) = patom                     ! (GeV/c)
       qclusts(  4, nclsts ) = sqrt( patom **2 + rmass**2 ) ! (GeV)
       qclusts(  5, nclsts ) = rmass                     ! (GeV)
       qclusts(  6, nclsts ) = 0.d0                      ! from (eV) to (MeV)
       qclusts(  7, nclsts ) = eatom                     ! from (eV) to (MeV)
       qclusts(  8, nclsts ) = 1.d0                      ! weight change
       qclusts(  9, nclsts ) = 0.d0                      ! (ns)
       qclusts( 10, nclsts ) = 0.d0                      ! x displacement (cm)
       qclusts( 11, nclsts ) = 0.d0                      ! y displacement (cm)
       qclusts( 12, nclsts ) = 0.d0                      ! z displacement (cm)

c       egs  =    enp * 1.d-6 - eatom
c       qclusts(  4, 1 ) = (enp * 1.d-6 - eatom + rmtyp(ityp,ktyp) ) 
c     &                     * 1.d-3 ! (GeV)
c       qclusts(  7, 1 ) = enp * 1.d-6 - eatom      ! from (eV) to (MeV)
      endif

      do i = 1, lng_rel ! atomic relaxation

! 2024-02-08 epsilon value correction for calculating charged particles
       if( ktp_rel(i) .ne. 22) then ! Auger electron
         if(eng_rel(i) .le. etsmin )then
           eleemit = 1.d0
           if(wvlin.gt.0.d0) then
             if(eng_rel(i)*1.d6 .gt. wvlin) then
              eleemit = dble(int(eng_rel(i)*1.d6/wvlin))
              if(eleemit.lt.1.d0) eleemit = 1.d0
             endif
             eng_rel(i) = 0.d0
             ehnum = ehnum + eleemit
           endif
          endif
        endif

       nclsts = nclsts + 1

        if( ktp_rel(i) .eq. 22) then ! emit X-ray
         iclusts( nclsts )    =  4
         jclusts( 3, nclsts ) = 14
         jclusts( 5, nclsts ) =  0
         numpal(14)           = numpal(14) + 1
         rumpal(14)           = rumpal(14) + wt(ibkwt+no,ipomp+1)
        else                         ! emit Auger electron
         iclusts( nclsts )    =  0
         jclusts( 3, nclsts ) = 12
         jclusts( 5, nclsts ) = -1
         numpal(12)      = numpal(12) + 1
         rumpal(12)      = rumpal(12) + wt(ibkwt+no,ipomp+1)
        endif

       jclusts(  0, nclsts ) = 0
       jclusts(  1, nclsts ) = 0
       jclusts(  2, nclsts ) = 0
       jclusts(  4, nclsts ) = 0
       jclusts(  6, nclsts ) = 0
       jclusts(  7, nclsts ) = ktp_rel(i)
       jclusts(  8, nclsts ) = 0

       rmass = rmtyp(ityp,jclusts(7,nclsts))
       p_ejc = sqrt(eng_rel(i)**2 + 2.d0 * rmass * eng_rel(i)) * 1.d-3
       rmass = rmass * 1.d-3

       phi_r = 2.d0 * physc(1) * unirn(dummy)
       cos1 = 1.d0 - 2.d0 * unirn(dummy)
       sin1 = sqrt( 1.d0 - cos1**2 )
       px = p_ejc * sin1 * sin(phi_r)
       py = p_ejc * sin1 * cos(phi_r)
       pz = p_ejc * cos1

       qclusts(  0, nclsts ) = 0.d0                      ! impact parameter
       qclusts(  1, nclsts ) = px                        ! (GeV/c)
       qclusts(  2, nclsts ) = py                        ! (GeV/c)
       qclusts(  3, nclsts ) = pz                        ! (GeV/c)
       qclusts(  4, nclsts ) = sqrt(p_ejc**2 + rmass**2) ! (GeV)
       qclusts(  5, nclsts ) = rmass                     ! (GeV)
       qclusts(  6, nclsts ) = 0.d0                      ! from (eV) to (MeV)
       qclusts(  7, nclsts ) = eng_rel(i)                ! from (eV) to (MeV)
       qclusts(  8, nclsts ) = 1.d0                      ! weight change
       qclusts(  9, nclsts ) = 0.d0                      ! (ns)
       qclusts( 10, nclsts ) = 0.d0                      ! x displacement (cm)
       qclusts( 11, nclsts ) = 0.d0                      ! y displacement (cm)
       qclusts( 12, nclsts ) = 0.d0                      ! z displacement (cm)

      enddo
! hirata electron generation count
      if(wvlin .gt. 0.d0 .and. ehnum.ne.0.d0) then
           atmrc(10,2) = atmrc(10,2) + ehnum
      endif

      return
      end subroutine
c================================================
      subroutine ele_dataup(ndele,eene)
c===============================================
      use MMBANKMOD

      implicit double precision (a-h,o-z)
      include 'param00.inc'
      common / clustt / nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common / clustw / jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /pnsave/ egs, uus, vvs, wws, wts, tms, nms, nct(3)
!$OMP THREADPRIVATE(/pnsave/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)
     
      numpal = 0
      rumpal = 0.d0
      
      do i = 1, ndele ! atomic relaxation
       nclsts = nclsts + 1
       iclusts( nclsts )    =  0
       jclusts( 3, nclsts ) = 12
       jclusts( 5, nclsts ) = -1
       numpal(12)      = numpal(12) + 1
       rumpal(12)      = rumpal(12) + wt(ibkwt+no,ipomp+1)

       jclusts(  0, nclsts ) = 0
       jclusts(  1, nclsts ) = 0
       jclusts(  2, nclsts ) = 0
       jclusts(  4, nclsts ) = 0
       jclusts(  6, nclsts ) = 0
       jclusts(  7, nclsts ) = 11
       jclusts(  8, nclsts ) = 0

       rmass = rmtyp(ityp,jclusts(7,nclsts))
       p_ejc = sqrt(eene**2 + 2.d0 * rmass * eene) * 1.d-3
       rmass = rmass * 1.d-3
       
       phi_r = 2.d0 * physc(1) * unirn(dummy)
       cos1 = 1.0 - 2.0 * unirn(dummy)
       sin1 = sqrt( 1.0 - cos1**2 )
       
       px = p_ejc * cos1 * sin(phi_r)
       py = p_ejc * cos1 * cos(phi_r)
       pz = p_ejc * sin1

       qclusts(  0, nclsts ) = 0.d0                      ! impact parameter
       qclusts(  1, nclsts ) = px                        ! (GeV/c)
       qclusts(  2, nclsts ) = py                        ! (GeV/c)
       qclusts(  3, nclsts ) = pz                        ! (GeV/c)
       qclusts(  4, nclsts ) = sqrt(p_ejc**2 + rmass**2) ! (GeV)
       qclusts(  5, nclsts ) = rmass                     ! (GeV)
       qclusts(  6, nclsts ) = 0.d0                      ! from (eV) to (MeV)
       qclusts(  7, nclsts ) = eene                ! from (eV) to (MeV)
       qclusts(  8, nclsts ) = 1.d0                      ! weight change
       qclusts(  9, nclsts ) = 0.d0                      ! (ns)
       qclusts( 10, nclsts ) = 0.d0                      ! x displacement (cm)
       qclusts( 11, nclsts ) = 0.d0                      ! y displacement (cm)
       qclusts( 12, nclsts ) = 0.d0                      ! z displacement (cm)
      
      enddo
      return
      end subroutine
************************************************************************
*                                                                      *
      subroutine Ion_paramset(j)
*                                                                      *
*       Set parameters a,A1,B1,C1,D1,E1,A2,B2,C2,D2                    *
*                                                                      *
*       Data source:                                                   *
*                                                                      *
*       Noble gas : M.E.Rudd, Review of Modern Physics, 64 (1992)      *
*                    441-491                                           *
*                                                                      *
*     input  :                                                         *
*          j : shell number                                            *
*                                                                      *
*     output :                                                         *
*                                                                      *
*          a,A1,B1,C1,D1,E1,A2,B2,C2,D2                                *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      Select Case (j)

        Case(0)       ! inner shell
          A1 = 1.25d0
          B1 = 0.50d0
          C1 = 1.00d0
          D1 = 1.00d0
          E1 = 3.00d0
          A2 = 1.10d0
          B2 = 1.30d0
          C2 = 1.00d0
          D2 = 0.00d0
         alp=  0.66d0

        Case(1,-20020)       ! He
          A1 = 1.02d0
          B1 = 2.40d0
          C1 = 0.70d0
          D1 = 1.15d0
          E1 = 0.70d0
          A2 = 0.84d0
          B2 = 6.00d0
          C2 = 0.70d0
          D2 = 0.50d0
         alp=  0.86d0

        Case(2:4,-20100)       ! Ne
          A1 = 0.58d0
          B1 = 65.0d0
          C1 = 0.23d0
          D1 = 0.55d0
          E1 = 0.16d0
          A2 = 1.40d0
          B2 = 0.00d0
          C2 = 0.72d0
          D2 = 1.35d0
         alp=  0.57d0

        Case(5:10,-20180)        ! Ar
          A1 = 1.20d0
          B1 = 8.00d0
          C1 = 0.86d0
          D1 = 0.00d0
          E1 = 0.80d0
          A2 = 0.90d0
          B2 = 2.70d0
          C2 = 0.75d0
          D2 = 0.80d0
         alp = 0.71d0

        Case(11:17,-20360)        ! Kr
          A1 = 1.46d0
          B1 = 5.70d0
          C1 = 0.65d0
          D1 =-0.55d0
          E1 = 1.00d0
          A2 = 1.30d0
          B2 = 22.0d0
          C2 = 0.95d0
          D2 =-1.00d0
         alp=  0.78d0

        Case(-1)       ! Water
          A1 = 0.97d0
          B1 = 82.0d0
          C1 = 0.40d0
          D1 =-0.30d0
          E1 = 0.38d0
          A2 = 1.04d0
          B2 = 17.3d0
          C2 = 0.76d0
          D2 = 0.04d0
         alp=  0.64d0

        Case(-2)       ! CO2
          A1 = 1.09d0
          B1 = 25.0d0
          C1 = 0.75d0
          D1 = 0.75d0
          E1 = 0.65d0
          A2 = 0.78d0
          B2 = 3.00d0
          C2 = 0.70d0
          D2 = 0.85d0
         alp=  0.53d0

        Case(-10001)       ! CH4
          A1 = 1.15d0
          B1 = 14.0d0
          C1 = 0.35d0
          D1 = 0.50d0
          E1 = 3.00d0
          A2 = 0.60d0
          B2 = 3.80d0
          C2 = 1.20d0
          D2 = 0.45d0
         alp=  0.61d0

        Case(-20010)       ! H2
          A1 = 0.96d0
          B1 = 2.60d0
          C1 = 0.38d0
          D1 = 0.23d0
          E1 = 2.20d0
          A2 = 1.04d0
          B2 = 5.90d0
          C2 = 1.15d0
          D2 = 0.20d0
         alp=  0.87d0

        Case(-20070)       ! N2
          A1 = 1.05d0
          B1 = 12.0d0
          C1 = 0.74d0
          D1 =-0.39d0
          E1 = 0.80d0
          A2 = 0.95d0
          B2 = 1.20d0
          C2 = 1.00d0
          D2 = 1.30d0
         alp=  0.70d0

        Case(-20080)       ! O2
          A1 = 1.02d0
          B1 = 50.0d0
          C1 = 0.40d0
          D1 = 0.12d0
          E1 = 0.30d0
          A2 = 1.00d0
          B2 = 5.00d0
          C2 = 0.55d0
          D2 = 0.00d0
         alp=  0.59d0

        Case default !(18:28)       ! Not given in the literature. Data of inner shell
          A1 = 1.25d0
          B1 = 0.50d0
          C1 = 1.00d0
          D1 = 1.00d0
          E1 = 3.00d0
          A2 = 1.10d0
          B2 = 1.30d0
          C2 = 1.00d0
          D2 = 0.00d0
         alp=  0.66d0

      end select

      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine shell_eng_highlow(ichemID,izt)
*                                                                      *
*       Select the highest and the lowest chemical potential           *
*                                                                      *
*     input                                                            *
*                                                                      *
*          ichemID : Chemical compound ID                              *
*          izt     : element atomic number                             *
*                                                                      *
*     output :                                                         *
*                                                                      *
*          eIlow  : lowest ionization energy                           *
*          eIhigh : highest ionization energy                          *
*                                                                      *
************************************************************************

      if(ichemID .le. 0) then
       eIlow =minval(pot_elem(1:28,izt),MASK=pot_elem(1:28,izt).gt.0.d0)
       eIhigh=maxval(pot_elem(1:28,izt))
      else
       eIlow  = 1.d30
       eIhigh = 0.d0
        do jj = 1, 28
          if(pot_mol(jj,ichemID) .eq. 0.d0) exit
          eIlow  = min(eIlow,  pot_mol(jj,ichemID))
          eIhigh = max(eIhigh, pot_mol(jj,ichemID))
        enddo
      endif

      return
      end subroutine


************************************************************************
*                                                                      *
      function emaxd(e_n)
*                                                                      *
*       delta-ray energy maximum                                       *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       e_n : projectile energy in eV/n                                *
*                                                                      *
*     output :                                                         *
*                                                                      *
*      emaxd   : delta-ray energy max (eV)                             *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      include 'param-physcnst.inc'

      rmase = rstms(12) * 1.d9
      rmass = rtyp * 1.d6
      eein  = e_n*ibryf(ityp,ktyp)
      
      emaxd = 4. * rmase / rmass * eein * ( 1. + eein / 2. / rmass )
     &     / ( ( 1. + rmase /rmass )**2 + 2. * rmase * eein / rmass**2 ) ! Delta-ray energy max taken from subroutine delprd

      return
      end function

************************************************************************
*                                                                      *
      function t(eng)
*                                                                      *
*       T function in the literature                                   *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       eng  : projectile energy in eV                                 *
*                                                                      *
*     output :                                                         *
*                                                                      *
*          t   : t function                                            *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      include 'param-physcnst.inc'

      t = eng /( rstms(1)*1.d3 / rmtyp(12,11) )

      return
      end function

************************************************************************
*                                                                      *
      function vf(eng)
*                                                                      *
*       v function in the literature                                   *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       eng  : projectile energy in eV                                 *
*                                                                      *
*     output :                                                         *
*                                                                      *
*        vf  : v function                                              *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      vf = sqrt( t(eng) / eI0 )

      return
      end function


************************************************************************
*                                                                      *
      function wc(e)
*                                                                      *
*       w_c function in the literature                                 *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       e     : projectile energy in eV                                *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       wc    : w_c function                                           *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      wc = 4.d0 * vf(e) **2 - 2.d0 * vf(e) - R / 4.d0 / eI0

      return
      end function

************************************************************************
*                                                                      *
      function smj(e)
*                                                                      *
*       J function in the literature                                   *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       e     : projectile energy in eV                                *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       j     : J function                                             *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      smj = Exp(alp * wc(e) / vf(e))

      return
      end function


************************************************************************
*                                                                      *
      function f1(e)
*                                                                      *
*       F1 function in the literature                                  *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       e     : projectile energy in eV                                *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       f1    : F1 function                                            *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      f1 = C1 * vf(e)**D1/(1.d0 + E1 * vf(e)**(D1 + 4.d0))
     &   + A1 * Log(1.d0 + vf(e)**2)/(vf(e)**2 + B1/vf(e)**2)

      return
      end function


************************************************************************
*                                                                      *
      function f2(e)
*                                                                      *
*       F2 function in the literature                                  *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       e     : projectile energy in eV                                *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       f2    : F2 function                                            *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      f2 = C2 * vf(e)**D2 * (A2/vf(e)**2 + B2/vf(e)**4)
     &   /(C2 * vf(e)**D2 + (A2/vf(e)**2 + B2/vf(e)**4))

      return
      end function

************************************************************************
*                                                                      *
      subroutine sdX(w_arg, e_n, IDcall, xsec, ebmean)
*                                                                      *
*       Single differential cross section  Eq.(6)                      *
*                                                                      *
*   input  :                                                           *
*                                                                      *
*      w_arg : reduced outgoing electron energy                        *
*              if IDcall = 1   (w_arg = w*eI0)  unreduced energy       *
*              if IDcall = 2   (w_arg = w/vf(e_n)) w/f reduced         *
*        e_n : projectile energy (eV/n)                                *
*    IDcall  : call position ID. Switches w_arg definition             *
*                                                                      *
*   output :                                                           *
*                                                                      *
*       xsec : calculated single differential cross section  (m^2)     *
*              To get absolute X-sec, multiply by (projectile charge)^2 *
*     ebmean : binding energy mean (eV)                                *
************************************************************************

      implicit double precision (a-h,o-z)

      xsec   = 0.d0
      ebmean = 0.d0
      
      if(frac(mat,1) .le. 0.d0) then ! unregistered material

        lem   = nint( dnel_das(kmat0+mat) )
        dens  = 0.d0
        denst = 0.d0
        itzt  = 0

        if( .not. allocated(xsec_brkdwn)) then
            allocate(xsec_brkdwn(28,0:lem) )
            xsec_brkdwn = 0.d0
        elseif( size(xsec_brkdwn,dim=2) .ne. lem+1) then
            deallocate(xsec_brkdwn)
            allocate(xsec_brkdwn(28,0:lem) )
        endif
        xsec_brkdwn = 0.d0 ! T.Sato 2022/11/30, should be initialize everytime, otherwise elastic cross section is contaminated

        do i = 0, lem

            if(i .eq. 0) then
              itz  = 1
              dens = denh_das(kmat0+mat)
            else
              itz  = nint( zz_das(kmat(mat)+i) )
              dens = den_das(kmat(mat)+i)
            endif

            if(dens .eq. 0.d0) cycle

            call shell_eng_highlow(-1,itz)

            do j = 1, 28

              elnum  = elec_elem(j,itz)

              if(elnum .eq. 0.d0) cycle

              eI0  = pot_elem(j,itz)   ! (sub)shell energy
              S    = 4.d0 * Pi * a0**2 * elnum * (R/eI0)**2

              if(eI0 .gt. e_n*dble(ibryf(ityp,ktyp))) cycle

              if(eI0/2.d0 .le. eIlow) then    ! outermost shell or shallow inner shell. Parameter is shell specific 
                  call Ion_paramset(j)
              else                               ! inner shell
                  call Ion_paramset(0)
              endif

              if( IDcall .eq. 1) then
                w = w_arg / eI0
              elseif( IDcall .eq. 2) then
                w = w_arg * vf(e_n)
              endif

              if(w .gt. 1.d3) then 
                xsec_brkdwn(j,i) = delmfp(rtyp,ctyp,dens*elnum*10**-24,
     &           e_n*dble(ibryf(ityp,ktyp))*1.d-6,1.d-3)
              else    
               arg1 = alp *(w - wc(e_n))/vf(e_n)
               if(arg1 .le. 708.d0) then ! avoid Exponetial overflow
!               xsec_brkdwn(j,i) = dens * elnum * S/eI0 * (f1(e_n) +   ! T.Sato 2022/11/28, elnum is already multiplied in S
                xsec_brkdwn(j,i) = dens * S/eI0 * (f1(e_n) +  
     &         f2(e_n) * w)/ (1.d0 + w)**3 / (1.d0 + Exp(arg1))
               else
                xsec_brkdwn(j,i) = 0.d0
               endif    
              endif                   
              
! 2021/08/30 Ogawa. This conversion is needed for OpenMP-Single agreement because rounding error depends on OpenMP flag.
              xsec_brkdwn(j,i) = real(xsec_brkdwn(j,i)) 
              ebmean = ebmean + eI0 * xsec_brkdwn(j,i)

              xsec = xsec + xsec_brkdwn(j,i)
              
            enddo

            denst = denst + dens
            itzt  = itzt + itz

        enddo

        ebmean = ebmean / max(xsec,1.d-78) ! 2024/3/18 avoid division by zero

        xsec = xsec / denst               ! T.Sato 2022/11/27, xsec should be microscopic

      else ! registered molecules

         if( .not. allocated(xsec_brkdwn)) then
            allocate(xsec_brkdwn(28,0:maxcomp) )
            xsec_brkdwn = 0.d0
         elseif( size(xsec_brkdwn,dim=2) .ne. maxcomp) then
            deallocate(xsec_brkdwn)
            allocate(xsec_brkdwn(28,0:maxcomp) )
         endif
         xsec_brkdwn = 0.d0 ! T.Sato 2022/11/30, should be initialize everytime, otherwise elastic cross section is contaminated

         i = 1

         do while(frac(mat,i) .gt. 0.d0) ! compounds

         call shell_eng_highlow(ichem(mat,i),0)

         do j = 1, 28 ! shell-wise cross sections

             elnum   = elec_mol(j,ichem(mat,i))
             if(elnum .eq. 0.d0) cycle

             eI0   = pot_mol(j,ichem(mat,i))   ! (sub)shell energy
             S     = 4.d0 * Pi * a0**2 * elnum * (R/eI0)**2 ! Absolute value is important T.Sato 2022/11/24

             if(eI0/2.d0 .le. eIlow) then    ! outermost shell or shallow inner shell. Parameter is shell specific
                call Ion_paramset(-ichem(mat,i))
             else                               ! inner shell
                call Ion_paramset(0)
             endif

             if( IDcall .eq. 1) then
               w = w_arg / eI0
             elseif( IDcall .eq. 2) then
               w = w_arg * vf(e_n)
             endif

              if(w .gt. 1.d3) then 
                xsec_brkdwn(j,i) = delmfp(rtyp,ctyp,frac(mat,i)*elnum
     &           *1.d-24,e_n*dble(ibryf(ityp,ktyp))*1.d-6,1.d-3)
               else    
              arg1 = alp *(w - wc(e_n))/vf(e_n)
              if(arg1 .le. 708.d0) then ! avoid Exponetial overflow
                xsec_brkdwn(j,i) = frac(mat,i) * S/eI0 * (f1(e_n) + 
     &                 f2(e_n) * w)/ (1.d0+w)**3 / (1.d0 + Exp(arg1))
                xsec = xsec + xsec_brkdwn(j,i)
              else
                xsec_brkdwn(j,i) = 0.d0
              endif                   
             endif

              ebmean = ebmean + eI0 * xsec_brkdwn(j,i)
         enddo

         i = i + 1

         enddo

         ebmean = ebmean / max(xsec,1.d-78) ! 2024/3/18 avoid division by zero

      endif

! T.Sato 2022/12/02, unknown correction factor to make dEdx(ITSART) < dEdx(ATIMA) in range.f
      xsec=xsec*0.9d0 ! this correction factor is determined for water, and should be investigated for other materials in the future

      return
      end subroutine


************************************************************************
*                                                                      *
*    / \                                                               *
*   /   \                                                              *
*     |                                                                *
*     |      Ionization part                                           *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*     |                                                                *
*     |      Excitation part                                           *
*   \   /                                                              *
*    \ /                                                               *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************



************************************************************************
*                                                                      *
      subroutine Exc_Xsection(izt, e_n, ichemID, Exc_Xsec)
*                                                                      *
*       Excitation cross section calculation                           *
*                                                                      *
*       Reference:                                                     *
*          J.H.Miller, A.E.S.Green, Radiation Research, 54 (1973)      *
*                    343-363                                           *
*                                                                      *
*   input  :                                                           *
*                                                                      *
*       izt  : target atomic number                                    *
*       e_n  : projectile energy (eV/n)                                *
*    ichemID : chemical compound ID (-1: no particular data)           *
*                                                                      *
*   output :                                                           *
*                                                                      *
*    Exc_Xsec(2) : calculated cross section (b)                        *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      double precision Exc_Xsec(nch_exc) ! excitation state energy

      Exc_Xsec = 0.d0

      if(ichemID .eq. -1) then ! T.Sato 2022/11/27
       elnum   = dble(izt)
      else
       elnum   = 0.d0
       do ielem = 1, 28
           elnum = elnum + dble(elec_mol(ielem,ichemID))
       enddo
      endif

      do j = 1, nch_exc ! Rydberg AB, CD, and Diffuse band

        call Exc_paramset(j, ichemID, izt, wMG, aMG, gJMG, omgMG,
     &                    gammMG)

        if(e_n - wMG .le. 0.d0) cycle
        Exc_Xsec(j) = (elnum*aMG)**omgMG*(e_n*1.d-3 - wMG*1.d-3)**gammMG ! Threshold of heavy ions are high (e.g. Hg-200 * 10 keV(wMG) = 2 MeV)
     &   / (gJMG**(omgMG+gammMG) + (e_n*1.d-3)**(omgMG+gammMG)) * 1.d8
     &   * 4.2d0 ! correction to consider unimplemented excitation channels. T.Ogawa 2024/01/22, change 20 to 10. 
                 ! change correction factor (2.0 -> 4.2) Y.Matsuya 2024/02/14

      enddo

      return
      end subroutine

************************************************************************
*                                                                      *
      subroutine Exc_paramset(imode, ichemID, izt, egap, par1, par2,
     & par3,par4)
*                                                                      *
*       Set parameters for excitation/charge-exchange cross sections   *
*                                                                      *
*       Data source:                                                   *
*                                                                      *
*          J.H.Miller, A.E.S.Green, Radiation Research, 54 (1973)      *
*                    343-363                                           *
*                                                                      *
*     input  :                                                         *
*      imode : mode number                                             *
*          1 -> Rydgerg AB                                             *
*          2 -> Rydgerg CD                                             *
*                                                                      *
*     output : cross section fitting parameters                        *
*          egap : energy gap (eV), wMG in excitation,                  *
*                 ionization energy in charge-changing reaction        *
*                                                                      *
*          par1 : fit param1, aMG in excitation (keV)                  *
*                                                                      *
*          par2 : fit param2, gJMG in excitation (keV)                 *
*                                                                      *
*          par3 : fit param3, omgMG in excitation                      *
*                                                                      *
*          par4 : fit param4, gammMG in excitation                     *
*                                                                      *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

c      Select Case (ichemID) To be extended later
c      Select Case (izt)     To be extended later

      Select Case (imode)

        Case(1)       ! Rydgerg AB
          If(jemi(mat) .eq. 1) then
            egap   = 10.0d0
          else
            egap   = 11.26d0
          endif
          par1   = 0.249d0
          par2   = 19.6d0
          par3   = 0.75d0
          par4   = 1.00d0

        Case(2)       ! Rydgerg CD
          If(jemi(mat) .eq. 1) then
            egap   = 11.1d0
          else
            egap   = 11.93d0
          endif
          par1   = 0.555d0
          par2   = 21.6d0
          par3   = 0.75d0
          par4   = 1.00d0

        Case(3)       ! Diffuse band
          egap   = 13.32d0
          par1   = 7.28d0
          par2   = 26.0d0
          par3   = 0.75d0
          par4   = 1.00d0

        Case default
          write(*,*) 'Error. Excitation of this material is undefined in
     & ITSART'

      end select

      return
      end subroutine

************************************************************************
*                                                                      *
      subroutine getsig_EX(en, e_n, sig_EX, eng_sig_EX)
*                                                                      *
*     get integral of excitation cross section                         *
*     macroscopic cross section is returned in different from getsdXsum *
*   input  :                                                           *
*        en  : projectile energy (eV)                                  *
*        e_n : projectile energy (eV/n)                                *
*                                                                      *
*   output :                                                           *
*    sig_EX  : macroscopic cross section of excitation (cm^-1)         *
*    eng_sig_EX : energy-weighted macroscopic cross section (eV.cm^-1)  *
************************************************************************

      double precision Exc_Xsec(nch_exc) ! excitation state energy
      double precision Xsec_Exc(nch_exc,maxcomp) 

      sig_EX = 0
      eng_sig_EX = 0

      if(frac(mat,1) .ne. 0.d0) then ! chem specified
          i          = 1
          do while(frac(mat,i) .gt. 0.d0)
            call Exc_Xsection(0,e_n,ichem(mat,i),Exc_Xsec)
            do j = 1, nch_exc
             xsec_exc(j,i) = Exc_Xsec(j)* getdenst(dummy) * frac(mat,i) ! Exc_Xsec (b) * getdenst(10^24 atom/cm^3) = (/cm)

             call Exc_paramset(j, ichemID, izt, wMG, aMG, gJMG, omgMG,
     &                    gammMG)
!<-20211126murofushi update
!--             tmass = denh(kmat(mat)+2) * rmtyp(1,2212) ! formula weight
!--             do k = 1, nint( dnel(kmat(mat)+1) ) ! number of elements
!--               dens = den(kmat(mat)+(k-1)*3+5)
!--               itz  = nint( zz(kmat(mat)+(k-1)*3+3) )
!--               ita  = nint( zz(kmat(mat)+(k-1)*3+4) )
             tmass = denh_das(kmat0+mat) * rmtyp(1,2212) ! formula weight
             do k = 1, nint( dnel_das(kmat0+mat) ) ! number of elements
               dens = den_das(kmat(mat)+k)
               itz  = nint( zz_das(kmat(mat)+k) )
               ita  = nint( a_das(kmat(mat)+k) )
!--->
               elmass =rmtyp(kftp(itz*1000000+ita),itz*1000000+ita)
               tmass = tmass + elmass * dens
             enddo
             tmass = tmass / density_calc(ichem(mat,i)) ! mass of molecule   MeV/c / Molecule 

             eprj  =  en - wMG                   ! Target recoil disregarded at first
             pprj  = sqrt(eprj**2 + 2.d0 * rtyp*1.d6 * eprj)
             ptarg = sqrt(en**2 + 2.d0 * rtyp*1.d6 * en) - pprj 
             etarg = sqrt( ptarg**2 + (tmass*1.d6)**2 ) - tmass*1.d6

             sig_Ex = sig_Ex + xsec_exc(j,i)
             eng_sig_Ex = eng_sig_Ex + (etarg + wMG) * xsec_exc(j,i)
            enddo  
            i = i + 1
          enddo

      else ! chem unspecified. Calculated as mixture of pure elements

!<-20211126murofushi update
!--          hydro = denh(kmat(mat)+2)
          hydro = denh_das(kmat0+mat)
!--->
          if(hydro .gt. 0.d0) then
            call Exc_Xsection(1,e_n,-1,Exc_Xsec)
            do j = 1, nch_exc
              xsec_exc(j,1) = Exc_Xsec(j) * hydro ! Contribution of 1H

              call Exc_paramset(j, -1, 1, wMG, aMG, gJMG, omgMG, gammMG) 
              tmass = rmtyp(1,2212)

              eprj  =  en - wMG                   ! Target recoil disregarded at first
              pprj  = sqrt(eprj**2 + 2.d0 * rtyp*1.d6 * eprj)
              ptarg = sqrt(en**2 + 2.d0 * rtyp*1.d6 * en) - pprj 
              etarg = sqrt( ptarg**2 + (tmass*1.d6)**2 ) - tmass*1.d6

              sig_Ex = sig_Ex + xsec_exc(j,1)
              eng_sig_Ex = eng_sig_Ex + (etarg + wMG) * xsec_exc(j,1)
            end do
          endif    

!<-20211126murofushi update
!--          lem   = nint( dnel(kmat(mat)+1) )
          lem   = nint( dnel_das(kmat0+mat) )
!--->
          if(lem .gt. 99) write(*,*) 'Error : maxcomp is insufficient.'
          do i = 1, lem

!<-20211126murofushi update
!--            itz  = nint( zz(kmat(mat)+(i-1)*3+3) )
!--            ita  = nint( zz(kmat(mat)+(i-1)*3+4) )
!--            dens = den(kmat(mat)+(i-1)*3+5)
            itz  = nint( zz_das(kmat(mat)+i) )
            ita  = nint( a_das(kmat(mat)+i) )
            dens = den_das(kmat(mat)+i)
!--->
            call Exc_Xsection(itz, e_n, -1,Exc_Xsec)

            do j = 1, nch_exc
              xsec_exc(j,1+i) = Exc_Xsec(j) * dens  

              call Exc_paramset(j, -1, 1, wMG, aMG, gJMG, omgMG, gammMG) 
              tmass = rmtyp(kftp(itz*1000000+ita),itz*1000000+ita)

              eprj  =  en - wMG                   ! Target recoil disregarded at first
              pprj  = sqrt(eprj**2 + 2.d0 * rtyp*1.d6 * eprj)
              ptarg = sqrt(en**2 + 2.d0 * rtyp*1.d6 * en) - pprj 
              etarg = sqrt( ptarg**2 + (tmass*1.d6)**2 ) - tmass*1.d6

              sig_Ex = sig_Ex + xsec_exc(j,1+i)
              eng_sig_Ex = eng_sig_Ex + (etarg + wMG) * xsec_exc(j,1+i)
            end do
          enddo
      endif

      end subroutine

c=================================================================
      subroutine itsreac2(mat1, ierr)
*                                                                      *
*     determine kinematics of excitation reaction                      *
*                                                                      *
*     input                                                            *
*                                                                      *
*      mat1 : track structure mode switch                           *
*                                                                      *
*                                                                      *
*     output                                                           *
*                                                                      *
*      outgoing proton energy/momentum sent through jclusts            *
*      ierr : error identifier                                         *
*                                                                      *
*                                                                      *
c=================================================================
c=================================================================
      use MMBANKMOD

      implicit double precision (a-h,o-z)

      include 'param00.inc'
      include 'param.inc'
      include 'param-physcnst.inc'

      parameter (pi=3.141592654d0)

      common / etsexe  / dexc_ene
!$OMP THREADPRIVATE(/etsexe/)
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)
      common /  tscreg   / ntscell(kvlmax)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      double precision xsec_exc(nch_exc,maxcomp) ! excitation cross section. AB and CD are considered
      double precision Exc_Xsec(nch_exc) ! excitation state energy

      dimension vel(3)


! hirata 2024-02-08 Counting generated electron with eps value
      common /etsart/ bgets(kvlmax), ebgets(kvlmax),
     &                wvets(kvlmax), ewvets(kvlmax)
      common /ets_wvalue/ wvlin
!$OMP THREADPRIVATE(/ets_wvalue/)
      common / etsminmax / etsmin, etsmax
! hirata 2024-02-08
      ehnum  = 0.d0
      wvlin = ewvets( idgr(iblz(ibkblz+no,ipomp+1)) )
      ii_chemID = 0  ! selected ichemID


      ierr = 0
      ec(ibkec+no,ipomp+1) = ec(ibkec+no,ipomp+1) * 1.d+6   ! from (MeV) to (eV)

      en  = ec(ibkec+no,ipomp+1) ! incident energy
      it = nty(ibknty+no,ipomp+1)
      kf = nkf(ibknkf+no,ipomp+1)

      e_n = en / dble(ibryf(ityp,ktyp)) ! T.Sato 2022/11/26, en is changed to e_n in this routine


      exmin = 1.d10
      do i = 1, nch_exc
       call Exc_paramset(i, ichemID, izt, egap, par1, par2, par3,par4)
       exmin = min(exmin, egap)
      enddo

      if(iflnel .eq. 1 .or. en .lt. exmin) then ! no reaction event. determined by restricted LET energy conservation
          pejc = sqrt( en**2 + 2.d0 * rmtyp(it,kf) * 1.d+6 * en)
          call itsdataup(0.d0,0.d0,0.d0,1.d0,pejc,en,0.d0,0.d0,0.d0,
     &                   0.d0,0.d0,0.d0,0,0)
          iflnel = 0
          ec(ibkec+no,ipomp+1) = en * 1.d-6 ! from (eV) to (MeV)
          dexc_ene     = e(ibke+no,ipomp+1) - ec(ibkec+no,ipomp+1)
          return
      endif


! determine cross section sum

      xsecSum = 0.d0
      if(frac(mat,1) .ne. 0.d0) then ! chem specified
        i         = 1
        do while(frac(mat,i) .gt. 0.d0)
          call Exc_Xsection(0,e_n,ichem(mat,i),Exc_Xsec)
          do j = 1, nch_exc
            xsec_exc(j,i) = Exc_Xsec(j)* density_calc(ichem(mat,i))
     &                * frac(mat,i)
            xsecSum = xsecSum + xsec_exc(j,i)
          enddo
          i = i + 1
        enddo

      else ! chem unspecified. Calculated as mixture of pure elements

        hydro = denh_das(kmat0+mat)
        if(hydro .gt. 0.d0) then
            call Exc_Xsection(1,e_n,-1,Exc_Xsec)
            do j = 1, nch_exc
              xsec_exc(j,1) = Exc_Xsec(j) * hydro ! Contribution of 1H
              xsecSum = xsecSum + xsec_exc(j,1)
            enddo
        endif

        lem   = nint( dnel_das(kmat0+mat) )
        if(lem .gt. 99) write(*,*) 'Error : maxcomp is insufficient.'
        do i = 1, lem

          itz  = nint( zz_das(kmat(mat)+i) )
          dens = den_das(kmat(mat)+i)
          call Exc_Xsection(itz, e_n, -1,Exc_Xsec)

          do j = 1, nch_exc
            xsec_exc(j,1+i) = Exc_Xsec(j) * dens
            xsecSum = xsecSum + xsec_exc(j,1+i)
        enddo

        enddo
      endif
        
      if(xsecSum .le. 0.d0) then ! this happens if energy is below threshold
          ec(ibkec+no,ipomp+1) = en * 1.d-6 ! from (eV) to (MeV)
          ierr = 2
          return
      endif    
      
      
! Randomly choose excitation channel and molecule/atom  START ****************************************

      xsecRan = xsecSum * unirn(dummy)
      if(frac(mat,1) .ne. 0.d0) then ! chem specified
          i          = 1
          do while(frac(mat,i) .gt. 0.d0)
          call Exc_Xsection(0,e_n,ichem(mat,i),Exc_Xsec)
          do j = 1, nch_exc
            xsecRan = xsecRan - Exc_Xsec(j)*
     &                density_calc(ichem(mat,i)) * frac(mat,i)
            if(xsecRan .lt. 0.d0) then
                tmass = denh_das(kmat0+mat) * rmtyp(1,2212) ! formula weight
                do k = 1, nint( dnel_das(kmat0+mat) ) ! number of elements
                  dens = den_das(kmat(mat)+k)
                  itz  = nint( zz_das(kmat(mat)+k) )
                  ita  = nint( a_das(kmat(mat)+k) )
                  elmass =rmtyp(kftp(itz*1000000+ita),itz*1000000+ita)
                  tmass = tmass + elmass * dens
                enddo
                tmass = tmass / density_calc(ichem(mat,i)) ! mass of molecule   MeV/c / Molecule
                goto 100
            endif
          enddo
            i = i + 1
          enddo

 100      call Exc_paramset(j, ichem(mat,i), 0, wMG, aMG, gJMG, omgMG,
     & gammMG)
        excE = wMG ! excitation energy of chosen channel
        ii_chemID = ichem(mat,i)
      else ! chem unspecified. Calculated as mixture of pure elements
        ii_chemID = 0
        itz  = 1
          hydro = denh_das(kmat0+mat)
        if(hydro .gt. 0.d0) then
          call Exc_Xsection(1,e_n,-1,Exc_Xsec)

          do j = 1, nch_exc
            xsecRan = xsecRan - Exc_Xsec(j)* hydro
            if(xsecRan .lt. 0.d0) then
                itz=1 ! T.Sato 2022/11/21
                ita=1 ! T.Sato 2022/11/21
                tmass = rmtyp(1,2212)
                goto 110
            endif
          enddo

        endif

        lem   = nint( dnel_das(kmat0+mat) )
        do i = 1, lem

            itz  = nint( zz_das(kmat(mat)+i) )
            dens = den_das(kmat(mat)+i)
            call Exc_Xsection(itz, e_n, -1,Exc_Xsec)

          do j = 1, nch_exc
            xsecRan = xsecRan - Exc_Xsec(j)* dens
            if(xsecRan .lt. 0.d0 .or. (j .eq. nch_exc .and. i .eq. lem)) ! condition in () because xsecRan did not go below 0. 
     &        then
                ita  = nint( a_das(kmat(mat)+i) )
                tmass  = rmtyp(kftp(itz*1000000+ita),itz*1000000+ita)
                goto 110
            endif
          enddo

        enddo

 110    call Exc_paramset(j, -1, itz, wMG, aMG, gJMG, omgMG, gammMG)
        excE = wMG ! excitation energy of chosen channel

      endif





      eprj  =  en - excE ! Target recoil disregarded at first
      pprj  = sqrt(eprj**2 + 2.d0 * rtyp*1.d6 * eprj)
      ptarg = sqrt(en**2 + 2.d0 * rtyp*1.d6 * en) - pprj
      etarg = sqrt( ptarg**2 + (tmass*1.d6)**2 ) - tmass*1.d6


! Reaction score
      atmrc(10,3) = atmrc(10,3) + 1.d0
*-----------------------------------------------------------------------
*        return data
*-----------------------------------------------------------------------

      ec(ibkec+no,ipomp+1) = (eprj - etarg)* 1.d-6  ! from (eV) to (MeV)

      dexc_ene = e(ibke+no,ipomp+1) - ec(ibkec+no,ipomp+1) ! energy deposit scored by T-deposit 
      lflgTS = 1

      if(dexc_ene .le. 0.d0) then ! this happens if energy loss is too small compared with energy recovery
          ec(ibkec+no,ipomp+1) = en * 1.d-6 ! from (eV) to (MeV)
          ierr = 3
          return
      endif
      
      call nuc_recoil_dataup(0.d0, 0.d0, 1.d0, ptarg*1.d-6, itz, ita,  
     &0.d0, 0.d0, 1.d0, pprj*1.d-6,1,1.d0)

! hirata 20240210 excitation electron treatment when mID = -1 and ichemID != 1. Excitement currently targets water only
      if(ii_chemID .ne. 1 .and. mat1 .eq. -1) then 
        eleemit = 1.d0
        if(wvlin.gt.0.d0) then
          eleemit = dble(int((dexc_ene*1.d+6)/wvlin))
          ehnum = ehnum + eleemit
          atmrc(10,2) = atmrc(10,2) + ehnum

        else ! currently 1 electron will emit when eps is not setted.
          bg  = ebgets( idgr(iblz(ibkblz+no,ipomp+1)) )
          if(bg.eq.0.d0 .or. (dexc_ene*1.d+6) .le. bg )  then
             atmrc(10,2) = atmrc(10,2) + 1
          else
            eemit = dexc_ene - bg*1.d-6  !eemit (MeV)
            if(eemit .gt. etsmin) then         !(MeV)
              call ele_dataup(1,eemit)
              dexc_ene = bg*1.d-6
            else
              atmrc(10,2) = atmrc(10,2) + 1
            endif
          endif
        endif
      endif
      

c      write(*,*) 'e',dexc_ene, e(ibke+no,ipomp+1), ec(ibkec+no,ipomp+1), 
c     & e2nd * 1.d-6, eatom, sum(eng_rel(1:lng_rel))
      
      end subroutine


************************************************************************
*                                                                      *
*    / \                                                               *
*   /   \                                                              *
*     |                                                                *
*     |      Excitation part                                           *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*     |                                                                *
*     |      Atom elastic collision part                               *
*   \   /                                                              *
*    \ /                                                               *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

************************************************************************
*                                                                      *
      subroutine Elastic_Xsec_Collparam(izt, ita, en, xsec, t1, lcycle, 
     & icycle, wei12)
*                                                                      *
*       Elastic scattering cross section calculated with Moliere       *
*       Moliere formula expanded by Mueller's method                   *
*                                                                      *
*       G.P.MUELLER, Radiation Effects, 50, 87-92 (1980)               *
*                                                                      *
*                                                                      *
*   input  :                                                           *
*                                                                      *
*       izt  : target atomic number                                    *
*       iza  : target mass number                                      *
*        en  : projectile energy (MeV)                                 *
*    lcycle  : If 1 no angle bias. If 2, angle bias.                   *
*    icycle  : If 1, forced particle. If 2 remainder particle          *
*                                                                      *
*                                                                      *
*   output :                                                           *
*                                                                      *
*    xsec    : cross section (b)                                       *
*    t1      : dimensionless t parameter (if 1.d10, t1 is not needed)  *
*    wei12   : Weight of forced particle                               *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      r2 = rmtyp(19, izt*1000000 + ita)

      as = (9.d0 * pi**2 / 128.d0)**(1.d0/3.d0) * a0 /
     &sqrt(dble(izt)**(2.d0/3.d0) + dble(ichgf(ityp,ktyp))**(2.d0/3.d0)) ! screening radius

      eps = en * as * r2 * 4.d0 * pi * physc(6) * 1.d6/
     &(dble(izt) * ichgf(ityp,ktyp) * physc(9) * (r2 + rtyp)) ! passed to ruth_momentum as a module variable 
      
      tmax = 4.d0 * r2 * rtyp * en * 1.d6 / (r2 + rtyp)**2 ! maximum transferred energy in eV

c Take threshold energy from T-DPA 
cc      if(frac(mat,1) .le. 0.d0 .or. frac(mat,1) .eq. 1.d0 ) then ! No chemical form specification or pure material
c        call dpath(dble(izt),ttd,1)
c        tmin = ttd
c       elseif(frac(mat,1) .eq. 1.d0 .and. ichem(mat,1) .gt. 20000) then ! Pure material. Use binding energy of the material
c        call dpath(dble((ichemID-20000)/10),ttd,1)
c        tmin = ttd
cc      elseif(frac(mat,1) .eq. 1.d0 ) then ! Pure material. Use binding energy of the material
cc       tmin = binding_energy_call_function_to_be_made_in_future(ichem(mat,1))
cc      else ! if material is complexed mixture, consider which molecule the atom belongs to.
c      endif
       tmin = 20.d0 ! minimum transferred energy (in eV) = crystal binding energy. Should be taken from DPA tally

      if(tmin .gt. tmax) then
         if(t1 .eq. 1.d10) then
           xsec = 0
         else
           t1 = 1.d-10
      endif
         return ! too low en. Return dummy output
      endif

      t0 = eps **2 * tmin/ tmax
      xsec = ruth_Xsec(as,sqrt(t0),eps) 
      
      if(t1 .eq. 1.d10) return ! only xsec is needed. Skip t1 sampling

      if(lcycle .eq. 2 .and. icycle .eq. 1) then ! bias run
c          thetmp = 
c          ttmp = eps**2 * sin(thetmp/2.d0)**2
            xsec1  = ruth_Xsec(as,eps*sin(rathe1/2.d0),eps) 

c          thetmp = 
c          ttmp = eps**2 * sin(thetmp/2.d0)**2
              xsec2  = ruth_Xsec(as,eps*sin(rathe2/2.d0),eps) 
          
              xsec12 = xsec1 - xsec2
              wei12  = xsec12/xsec
              xsec   = xsec12 
              ndiv = 49
              rat = (sin(rathe2/2.d0)/sin(rathe1/2.d0))
     &              **(1.d0/dble(ndiv+1))
              s   = eps*sin(rathe1/2.d0)
              t0  = s**2
      else
          if(icycle .eq. 2) wei12  = 1.d0 - wei12 ! bias run. Remainder part
          ndiv = 999
          rat = (eps/sqrt(t0))**(1.d0/dble(ndiv+1))
          s   = sqrt(t0)
      endif      
      
      s_save = s ! save s in case goto 10
 10   xsecr = xsec * unirn(dummy)
      s = s_save 

      do is = 0, ndiv
        s = s * rat
        xsecr = xsecr - ruth_Xsec(as,s,s*rat) 
        if(xsecr .le. 0 .or. is .eq. ndiv) exit
      enddo

      t1 = t0 * rat**(2*is) * rat**(2.d0*unirn(dummy))
      
      if(2.d0*asin(sqrt(t1)/eps) .gt. rathe1 .and. icycle .eq. 2 .and. 
     &   2.d0*asin(sqrt(t1)/eps) .lt. rathe2 ) goto 10 ! If within the biased angle, resample to aboid double counting. 

      return
      end subroutine

************************************************************************
*                                                                      *
      function ruth_Xsec(as,t1,t2)
*                                                                      *
*       Ruthrford particial cross section                              *
*                                                                      *
*   input  :                                                           *
*                                                                      *
*      t1  : dimensionless energy parameter                            *
*      t2  : dimensionless energy parameter                            *
*                                                                      *
*   output :                                                           *
*                                                                      *
*    XMoli   : X function of Mueller's expansion                       *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      
      ruth_Xsec = pi * (as * 1.d15)**2 * 1.d-2 * (XMoli(t1)-XMoli(t2))
      
      return
      end function
      

************************************************************************
*                                                                      *
      function XMoli(s)
*                                                                      *
*       Mueller's expansion of elastic cross section                   *
*                                                                      *
*   input  :                                                           *
*                                                                      *
*       s  : sqrt(t), magnitude of energy transfer                     *
*                                                                      *
*   output :                                                           *
*                                                                      *
*    XMoli   : X function of Mueller's expansion                       *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      if(s .le. 6.d-2) then
       XMoli = 148.298d0 - (0.5d0 * aMoli0 * (log(s))**2 + aMoli1
     & * log(s) + aMoli2 * s + 0.5d0 * aMoli3 * s**2 )
      else
       XMoli = ((bMoli2 * bMoli3 - 2.d0 * bMoli1 * (bMoli3**2 -2.d0 *
     &  bMoli2)) * pMoli2(s) + 4.d0 * bMoli1 * bMoli2 / s + (bMoli2 -
     &  2.d0 * bMoli1 * bMoli3) * (pMoli1(s) - 2.d0 * log(s) ))
     &  / (4.d0 * bMoli2**2)
      endif

      return
      end function

************************************************************************
*                                                                      *
      subroutine moliere_stopXsec(izt, ita, en, e_xsec)
*                                                                      *
*       Elastic scattering cross section weighted by energy loss       *
*       (stopping cross section)                                       *
*       calculated with Moliere formula expanded by Mueller's method   *
*                                                                      *
*       G.P.MUELLER, Radiation Effects, 50, 87-92 (1980)               *
*                                                                      *
*                                                                      *
*   input  :                                                           *
*                                                                      *
*       izt  : target atomic number                                    *
*       iza  : target mass number                                      *
*        en  : projectile energy (MeV)                                 *
*                                                                      *
*   output :                                                           *
*                                                                      *
*    xsec    : stopping cross section (eV * b)                         *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      r2 = rmtyp(19, izt*1000000 + ita)

      as = (9.d0 * pi**2 / 128.d0)**(1.d0/3.d0) * a0 /
     &sqrt(dble(izt)**(2.d0/3.d0) + dble(ichgf(ityp,ktyp))**(2.d0/3.d0)) ! screening radius

      eps = en * as * r2 * 4.d0 * pi * physc(6) * 1.d6/
     &(dble(izt) * ichgf(ityp,ktyp) * physc(9) * (r2 + rtyp)) ! passed to ruth_momentum as a module variable 
      
      tmax = 4.d0 * r2 * rtyp * en * 1.d6 / (r2 + rtyp)**2 ! maximum transferred energy in eV

      tmin = 20.d0 ! minimum transferred energy (in eV) = crystal binding energy. Should be taken from DPA tally

      t0 = eps **2 * tmin/ tmax

      e_xsec = pi * (as * 1.d15)**2 * 1.d-2 * (tmax / eps **2)
     & * ( YMoli(eps) - YMoli(sqrt(t0)) )
c     & * (en/1.d-3)**0.1d0 ! unknown correction factor to match stopping cross section of Miller and Green

      return
      end subroutine

************************************************************************
*                                                                      *
      function YMoli(s)
*                                                                      *
*       Mueller's expansion of elastic stopping cross section          *
*                                                                      *
*   input  :                                                           *
*                                                                      *
*       s  : sqrt(t), magnitude of energy transfer                     *
*                                                                      *
*   output :                                                           *
*                                                                      *
*    YMoli   : Y function of Mueller's expansion                       *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      if(s .le. 6.d-2) then
       YMoli = 0.5d0 * aMoli0 * s ** 2 * (log(s) - 0.5d0) + 0.5d0
     &  * aMoli1 * s**2 + aMoli2 * s**3 / 3.d0 + aMoli3 * s**4 / 4.d0
      else
       YMoli = 5.9298d-2 + 1.d0/4.d0 * ((4.d0 * bMoli1 - bMoli3)
     &  * pMoli2(s) + pMoli1(s))
      endif

      return
      end function

************************************************************************
*                                                                      *
      function pMoli1(s)
*                                                                      *
*       Mueller's expansion function                                   *
*                                                                      *
*   input  :                                                           *
*                                                                      *
*       s  : sqrt(t), magnitude of energy transfer                     *
*                                                                      *
*   output :                                                           *
*                                                                      *
*    pMoli1   : p1 function of Mueller's expansion                     *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      pMoli1 = log(bMoli2 + s * bMoli3 + s**2)

      return
      end function

************************************************************************
*                                                                      *
      function pMoli2(s)
*                                                                      *
*       Mueller's expansion function                                   *
*                                                                      *
*   input  :                                                           *
*                                                                      *
*       s  : sqrt(t), magnitude of energy transfer                     *
*                                                                      *
*   output :                                                           *
*                                                                      *
*    pMoli2   : p2 function of Mueller's expansion                     *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      pMoli0 = sqrt(bMoli3**2 - 4.d0 * bMoli2)

      pMoli2 = log((bMoli3 + 2.d0 * s - pMoli0)/
     &             (bMoli3 + 2.d0 * s + pMoli0)) / pMoli0

      return
      end function

************************************************************************
*                                                                      *
      subroutine ruth_momentum(b, aa, en, weitar, pscat, ptarg, thet)
*                                                                      *
*       Rutherford scattering cross section                            *
*                                                                      *
*   input  :                                                           *
*                                                                      *
*         b  : dimensionless collision parameter t                    *
*        aa  : electric potential per distance (fm*MeV)                *
*        en  : projectile energy (MeV)                                 *
*    weitar  : target mass (MeV/c**2)                                  *
*                                                                      *
*   output :                                                           *
*                                                                      *
*    pscat   : projectile momentum (MeV/c)                             *
*    ptarg   : target momentum (MeV/c)                                 *
*     thet   : projectile scattering angle (rad)                       *
************************************************************************

      implicit double precision (a-h,o-z)

      thet = 2.d0 * asin(sqrt(b/eps**2))
      thet = atan(sin(thet)/(rtyp/weitar + cos(thet))) ! CM to Lab conversion

      if(thet .lt. 0.d0 ) thet = thet + pi

      etarg = b/eps**2 * 4.d0 * rtyp * weitar / (rtyp + weitar)**2
     &   * en
      ptarg = sqrt(etarg**2 + 2.d0 * etarg * weitar)
      escat = en - etarg
      pscat = sqrt(escat**2 + 2.d0 * escat * rtyp)

      return
      end subroutine

c=================================================================
      subroutine iruthreac(ierr)
*                                                                      *
*     Determine Rutherford scattering kinematics                       *
*     momentum and energy of target nucleus and scattered projectile   *
*                                                                      *
*     input                                                            *
*      None                                                            *
*                                                                      *
*     output                                                           *
*      outgoing target/projectile energy/momentum sent through jclusts *
*      ierr : error identifier                                         *
*                                                                      *
c=================================================================
c=================================================================
      use MMBANKMOD

      implicit double precision (a-h,o-z)

      include 'param00.inc'
      include 'param.inc'
      include 'param-physcnst.inc'

      common / etsexe  / dexc_ene
!$OMP THREADPRIVATE(/etsexe/)
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)

      dimension vel(3)

      ierr = 0
      en = ec(ibkec+no,ipomp+1) ! incident energy

***** Cross section array allocation ************
      if( .not. allocated(xsec_brkdwn)) then
         allocate(xsec_brkdwn(28,0:maxcomp) )
      elseif( size(xsec_brkdwn,dim=2) .ne. maxcomp) then
         deallocate(xsec_brkdwn)
         allocate(xsec_brkdwn(28,0:maxcomp) )
      endif
      xsec_brkdwn = 0.d0

***** Rutherford scattering Cross section of elements ************
      hydro = denh_das(kmat0+mat)
      if(hydro .gt. 0.d0) then
        call Elastic_Xsec_Collparam(1, 1, en, xsec, 1.d10,1,1,1.d0)
        xsec_brkdwn(1,1) = xsec * hydro ! Contribution of 1H
      endif

      lem   = nint( dnel_das(kmat0+mat) )
      if(lem .gt. maxcomp) write(*,*) 'Increase maxcomp. Risk of access
     & violation'
      do i = 1, lem
        itz  = nint( zz_das(kmat(mat)+i) )
        ita  = nint( a_das(kmat(mat)+i) )
        dens = den_das(kmat(mat)+i)
        call Elastic_Xsec_Collparam(itz, ita, en, xsec, 1.d10,1,1,1.d0)
        xsec_brkdwn(1,1+i) = xsec * dens   ! rutherford scattering cross section
      enddo
**************************************************************

      sigtot = sum(xsec_brkdwn(1,1:1+lem)) * unirn(dummy)
      do i = 0, lem
        sigtot = sigtot - xsec_brkdwn(1,1+i)   ! rutherford scattering cross section
        if(sigtot .le. 0.d0 .or. i .eq. lem) then
         if(i .eq. 0) then
           itchar = 1
           itmass = 1
           ktyptar = 2212
           ityptar = 1
         else
           itchar = nint( zz_das(kmat(mat)+i) )
           itmass = nint( a_das(kmat(mat)+i) )
           ktyptar = itchar * 1000000 + itmass
           ityptar = kftp(ktyptar)
         endif
         exit
        endif
      enddo

      lcycle = 1 
      if(rathe1 .lt. rathe2) lcycle = 2 ! angular bias is active 
      weifc = 1.d0
      
      do ii = 1, lcycle 
      call Elastic_Xsec_Collparam(itchar, itmass, en, xsec, t1, lcycle, 
     & ii, weifc)

      weitar  = rmtyp(ityptar,ktyptar)

      call ruth_momentum(t1, aa, en, weitar, pscat, ptarg, thet)

      etarg = sqrt(ptarg**2 + weitar**2) - weitar
      thettarg = 0.d0


      phitarg = 2.d0 * unirn(dummy) * pi

      call vel_vec(vel,thettarg,phitarg)

      utarg = vel(1)
      vtarg = vel(2)
      wtarg = vel(3)

      phiscat = 2.d0 * unirn(dummy) * pi

      call vel_vec(vel,thet,phiscat)

      uscat = vel(1)
      vscat = vel(2)
      wscat = vel(3)

      atmrc(10,1) = atmrc(10,1) + 1.d0
*-----------------------------------------------------------------------
*        data up in bank
*-----------------------------------------------------------------------

      call nuc_recoil_dataup(utarg, vtarg, wtarg, ptarg, itchar, itmass,
     &uscat, vscat, wscat, pscat, ii, weifc)

      enddo ! end of angular bias
      
      ec(ibkec+no,ipomp+1) = sqrt(pscat**2 + rmtyp(ityp,ktyp)**2)
     & - rmtyp(ityp,ktyp)

      dexc_ene = e(ibke+no,ipomp+1) - ec(ibkec+no,ipomp+1) - etarg
      
      lflgTS = 1

      end subroutine


      end module       
