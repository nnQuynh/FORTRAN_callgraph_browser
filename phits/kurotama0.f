c=======================================================================
      subroutine KUROTAMA0 (ap0, zp0, Tp, at0, zt0, L, bs_xsec, bmax)
c              ----- Improved Hybrid-Kurotama Model -----
c                                            coded by K. Iida     (2007)
c                                              released  (Jul. 24, 2009)
c                 modified to function form by H. Iwase  (Aug. 25, 2011)
c                                  modified by A. Kohama (Sep. 20, 2011)
c                                  modified by A. Kohama (Nov. 04, 2011)
c                                  modified by A. Kohama (Dec. 20, 2011)
c                                  modified by A. Kohama (Jan. 20, 2012)
c                                  modified by L. Sihver (Feb. 24, 2014)
c                                  modified by A. Kohama (Mar. 01, 2014)
c     see eq. (2.7)
c     The formula is applicable for 100 MeV =< e =< 5 GeV
c     Ref.1: K. Iida, A. Kohama, and K. Oyamatsu,
c            J. Phys. Soc. Japan 76, 044201 (2007).
c
c     Make Kurotama hybrid with the Tripathi formula.
c     Ref.2: A. Kohama, K. Iida, K. Oyamatsu, H. Iwase,
c            S. Hashimoto, and K. Niita
c     Subroutine Kurotama in PHITS
c     RIKEN Accelerator Progress Report 46, p. 55 (2012).
c        http://www.nishina.riken.jp/researcher/APR/
c
c     The formulation in detail.
c     Ref.3: A. Kohama, K. Iida, and K. Oyamatsu,
c     Energy and mass-number dependence of hadron-nucleus
c     total reaction cross sections, arXiv:1411.7737
c
c     Improve hybrid-kurotama for the reactions
c     involving alpha particle ie., p+alpha, alpha+A
c     by putting ecut = 400 MeV/nucleon.
c     other corr.: mass number can be non-integer (natural abundance)
c                  energy >= 5 GeV/nucleon without limit
c     Ref.4: L. Sihver A. Kohama, K. Iida, K. Oyamatsu,
c            H. Iwase, S. Hashimoto, K. Niita,     RIKEN-NC-NP-125
c     Nucl. Instr. Meth., B 334 (2014), pp. 34-39
c
c     input : (Z,A,Tp,a0,L) Tp (MeV/nucleon); a0 (fm); L (MeV)
c              Ap0,At0: projectile/target mass number.
c              Zp0,Zt0: atomic number.
c              if Ap0 > At0
c              then make it Ap0 < At0 automatically in this code.
c     output: bs_xsec (fm^2), bmax (fm)
c     Note  : use subroutine nasa(...)
c             Currently, no contribution from the Coulomb breakup
c             is included,
c             which affects the reactions between heavy nuclei.
c=======================================================================
c       L: symmetry energy coefficient (0-100) MeV
c       stable nuclei L can be 0
c       (L should be a PHITS parameter)
c=======================================================================

      Implicit NONE
      real*8 tp,a0,L,bs1,bs2,bs_radius,a1,a   ! revised
      real*8 z,z1,ap0,zp0,at0,zt0             ! revised
      real*8 bs_xsec,bmax,tmp,ep
      Real*8 PI/3.1415926535897932384626434D0/
      Integer i, i1            ! j
      Integer a2,a3                           ! revised
      real*8 fn_a0
      real*8 a01
      real*8 signe, sigel
      real*8 cut1, cut2, cutoff_p, cutoff_m, ecut, width
      real*8 bs10, bs20, signe0, epcut

      a01= -50.0d0           ! fake parameter for bs_radius(...)
      a0 = -50.0d0           ! fake parameter for bs_radius(...)
c----------------------------------------------------------------------
c     this routine requires that at > ap
c----------------------------------------------------------------------
c  correction by L. Sihver  (Feb. 17, 2014)
      if(ap0.gt.at0)then
             a2 = nint(at0)
             a1 = real(a2)
             z1 = zt0
             a3 = nint(ap0)
             a = real(a3)
             z = zp0
      else
             a2 = nint(ap0)
             a1 = real(a2)
             z1 = zp0
             a3 = nint(at0)
             a =  real(a3)
             z = zt0
      end if

c set a1 < a :  target must be heavier than projectile

c recommended values of "a0", i.e., "a" at 800 MeV.
      if( a1.eq. 4  .and. z1.eq. 2 ) a01 = 1.43d0   ! tentative
      if( a .eq. 4  .and. z .eq. 2 ) a0  = 1.43d0   ! tentative
      if( a1.eq. 9  .and. z1.eq. 4 ) a01 = 2.55d0   !+-0.05d0
      if( a .eq. 9  .and. z .eq. 4 ) a0  = 2.55d0   !+-0.05d0
      if( a1.eq. 12 .and. z1.eq. 6 ) a01 = 2.70d0   !+-0.05d0
      if( a .eq. 12 .and. z .eq. 6 ) a0  = 2.70d0   !+-0.05d0
      if( a1.eq. 20 .and. z1.eq.10 ) a01 = 3.298d0  !+-0.07d0
      if( a .eq. 20 .and. z .eq.10 ) a0  = 3.298d0  !+-0.07d0
      if (a1.eq. 40 .and. z1.eq.20 ) a01 = 4.25d0
      if (a .eq. 40 .and. z .eq.20 ) a0  = 4.25d0
      if (a1.eq. 48 .and. z1.eq.20 ) a01 = 4.50d0
      if (a .eq. 48 .and. z .eq.20 ) a0  = 4.50d0
      if( a1.eq. 90 .and. z1.eq.40 ) a01 = 5.70d0
      if( a .eq. 90 .and. z .eq.40 ) a0  = 5.70d0
      if( a1.eq.208 .and. z1.eq.82 ) a01 = 7.45d0
      if( a .eq.208 .and. z .eq.82 ) a0  = 7.45d0

c cutoff parameters and cutoff functions
c introduced by A. Kohama, K. Iida, and K. Oyamatsu
c when make kurotama hybrid
      ecut  = 115.0d0   ! best choice
      if( a1 .eq.4 .or. a .eq.4 ) ecut = 400.0d0
c improved by L. Sihver  (Feb. 20, 2014)
      width = 1.0d0
      cut1  = cutoff_m(tp, ecut, width)
      cut2  = cutoff_p(tp, ecut, width)

c choice of A-dep of "a_0" for fn_a0(...)
        if ( a1 .ge. 50.0d0) then
           i1 = 3
        else   ! a1 < 50.0d0
           i1 = 2
        end if
        if ( a .ge. 50.0d0) then
           i = 3
        else   ! a < 50.0d0
           i = 2
        end if

c n + A case    tentative: use only Tripathi formula
      if( a1.eq.1  .and. z1.eq.0 )then
        ep      = a1*tp
        call nasa(a1,z1,ep,   a,z,signe, sigel,bmax)  ! return in barn
        bs_xsec = signe*100.0d0
        bmax    = dsqrt(bs_xsec/pi)
        return
      endif

c p + A case
      if( a1.eq.1  .and. z1.eq.1 )then
         if ( a .eq. 2 .and. z .eq. 1) a0 = -50.0d0
         if ( a .eq. 4 .and. z .eq. 2) a0 = -50.0d0
         tmp  = fn_a0(A,a0,i)
         bs1  = 0.0d0
         bs10 = 0.0d0
         bs2  = bs_radius(z,a,Tp,  tmp,L)
         bs20 = bs_radius(z,a,ecut,tmp,L)
c d + A case (in future: tentative prescription can be found above)
c 3He + A case (in future: tentatively Carlson-1 is adopted)
c a + A case (in future: tentative prescription can be found below)
c A + A case
      else
         tmp  = fn_a0(A1,a01,i1)
         bs1  = bs_radius(z1,a1,Tp,  tmp,L)
         bs10 = bs_radius(z1,a1,ecut,tmp,L)
         tmp  = fn_a0(A,a0,i)
         bs2  = bs_radius(z,a,Tp,  tmp,L)
         bs20 = bs_radius(z,a,ecut,tmp,L)
      endif

c Tripathi formula
      ep    = a1*tp
      epcut = a1*ecut
      call nasa(a1,z1,epcut,a,z,signe0,sigel,bmax)  ! return in barn
      call nasa(a1,z1,ep,   a,z,signe, sigel,bmax)  ! return in barn

c hybrid kurotama
      bs_xsec = cut1*pi*(bs1+bs2)**2.0
     +        + cut2*signe*pi*(bs10+bs20)**2.0/signe0

c d + A case (tentatively use the Tripathi formula) not good
c a + A case (tentatively use the Tripathi formula)
c following to the correction of ecut by L. Sihver  (Feb. 20, 2014)

      bmax = dsqrt(bs_xsec/pi)

      return
      end

c=======================================================================
        function bs_radius(Z,A,Tp,a0,L)
c                                            coded by K. Iida     (2007)
c                                              released  (Jul. 24, 2009)
c                 modified to function form by H. Iwase  (Aug. 25, 2011)
c                                  modified by A. Kohama (Sep. 20, 2011)
c                                  modified by A. Kohama (Nov. 04, 2011)
c                                  modified by A. Kohama (Dec. 20, 2011)
c                                  modified by A. Kohama (Feb. xx, 2014)
c     see eq. (2.7)
c     in K. Iida, A. Kohama, and K. Oyamatsu,
c         J. Phys. Soc. Japan 76, 044201 (2007).
c     applicable for 100 MeV =< e =< 5 GeV
c
c     improved:  energy >= 5 GeV/nucleon without limit  201402
c
c     input: (Z,A,Tp,a0,L,i) Tp (MeV); a0 (fm); L (MeV)
c     output: bs_radius (fm)
c=======================================================================

        Implicit NONE
        Real*8 Z,A,N,sigpn,sigpp,Tp
        Real*8 a0,sigpp0,sigpn0,tau,D,rho0,C,R0,R,L0,dLda0,dsig,sig0
        Real*8 da,f,sigmar,nc0
        Real*8 L,ff,rp,rn,rp0,rn0,rm,rm0,delta,rcfs,rmfs
        Real*8 PI/3.1415926535897932384626434D0/
        real*8 bs_radius
        real*8 sigpp_fit, sigpn_fit
        Integer i,j

        if(A .gt. 270.d0) then
            bs_radius = A**(1.d0/3.d0) ! avoid floating point error
            return
        endif

        sigpp = sigpp_fit(tp)
        sigpn = sigpn_fit(tp)

        N    =A-Z
        delta=1D0-2D0*Z/A

        rp  =0.915D0*A**(1D0/3D0)-0.102D0+0.389D0*(delta-0.880D0)**2
        rn  =0.880D0*A**(1D0/3D0)*(1D0+0.00635D0*L*delta**2
     &       -0.000172D0*L**2*delta**4)+0.302D0+0.193D0*delta
        rm  =sqrt(Z/A*rp**2+N/A*rn**2)
        rp0 =0.915D0*A**(1D0/3D0)-0.102D0+0.389D0*(-0.880D0)**2
        rn0 =0.880D0*A**(1D0/3D0)+0.302D0
        rm0 =sqrt(Z/A*rp0**2+N/A*rn0**2)
        ff  =(rm/rm0)**2
        rmfs=sqrt(rm**2+1.5D0*0.65D0**2)
        rcfs=sqrt(rp**2+1.5D0*0.65D0**2)



        sigpp0=4.7D0
        sigpn0=3.8D0

        tau  =0.9D0    ! 0.8D0
        D    =2.2D0
        rho0 =0.16D0
        C    =rho0/D
        R0   =(4D0*PI*rho0/3D0)**(-1D0/3D0)*A**(1D0/3D0)
        R    =R0+D/2-R0/(1D0+12D0*R0**2/D**2)
        L0   =2D0*sqrt(R**2-a0**2) ! path length
        dLda0=2D0/sqrt(R**2-a0**2)
     &        *((R0/a0)*R
     &        *(1D0-1D0/(1D0+12D0*R0**2/D**2)
     &        +24D0*R0**2/D**2/(1D0+12D0*R0**2/D**2)**2)-a0)

        dsig=Z/A*sigpp+(A-Z)/A*sigpn-Z/A*sigpp0-(A-Z)/A*sigpn0
        sig0=Z/A*sigpp0+(A-Z)/A*sigpn0
        nc0 =tau/L0/sig0
        da  =a0*dsig/sig0/(a0*C/nc0-a0*dLda0/L0)  ! /tau
        f   =(1D0+da/a0)**2



        bs_radius = a0 * dsqrt(f*ff)           ! fm

        End
c=======================================================================
        function fn_a0(A, a0, i)
c                                     coded by A. Kohama (Sep. 15, 2011)
c     see
c     K. Iida, A. Kohama, and K. Oyamatsu,
c         J. Phys. Soc. Japan 76, 044201 (2007).
c         Phys. Rev. C 72, 024602 (2005).
c
c     input: A    mass number
c            a0 [fm]    see below
c     output: fn_a0 (fm)   black-sphere radius at 800 MeV
c=======================================================================
c       a0 is automatically selected from table
c       if a0 not in table, choose scaling by i
c       (i should be a PHITS parameter)
c
c       i = 1 simple   scaling  suitable for general use
c             BS Scaling
c         = 2 modified scaling1 suitable for A =< 50
c             BS-fit1
c         = 3 modified scaling2 suitable for A >  50
c             BS-fit2
c=======================================================================

      implicit none
      include 'err.inc'

      real*8   A, a0, fn_a0
      integer  i

        if (i.ge.4 .or. i.le.0) then
         write(ErrCha,*)'error in i; i must be 1, 2, or 3'
         ErrID = 'L:294/R:fn_a0/F:kurotama0.f' !E00_007_001
         call ErrWrite(ErrID,ErrCha)
         stop
      endif

        if (a0.gt.0.0d0) then
          fn_a0 = a0
c if a0 not inputted, use scaling with 3 options
          else if (i.eq.1) then  ! simple scaling
            fn_a0 = 1.2135D0*A**(1D0/3D0)
          else if (i.eq.2) then ! modified scaling1
            fn_a0 = 1.2671D0*A**(1D0/3D0) - 0.152D0  !+-0.096d0
          else        ! i = 3   ! modified scaling2
            fn_a0 = 1.33D0*A**(1D0/3D0)   - 0.35D0
        end if
        return
        end

c=======================================================================
      function sigpp_fit(e)
c                                     coded by H. Iwase  (Aug. 25, 2011)
c                                  modified by A. Kohama (Sep. 12, 2011)
c                                  modified by A. Kohama (Dec. 08, 2011)
c                                  modified by A. Kohama (Dec. 20, 2011)
c                                  modified by A. Kohama (Feb. 04, 2014)
c     Returns the value of sigma_pp^total as a fn. of energy in lab.
c     see eq. (1)
c     in PRC.81.064603 C. A. Bertulani and C. De Conti (2010)
c     applicable for 10 MeV =< e =< 5 GeV
c     Note:  0 MeV < e < 10 MeV: Returns the value at 10 MeV
c            5 GeV < e         : Error   (up to v5.2)
c            5 GeV < e         : Returns the value at 5 GeV
c              Following to the suggestion by L. Sihver  (Feb. 04, 2014)
c     input: e (MeV)   output: sigpp_fit (fm^2)
c=======================================================================
      implicit none
      real*8   e,sigpp_fit

      if    ( e .lt. 10.0d0 )then
         sigpp_fit = 326.7006d0
      elseif( e .ge. 10.0d0 .and. e .lt. 280.0d0 )then
         sigpp_fit = 19.6d0 + 4253.0d0/e - 375.0d0/sqrt(e) + 3.86d-2*e

      elseif( e .ge. 280.0d0  .and.  e .lt. 840.0d0 )then
         sigpp_fit = 32.7d0 - 5.52d-2*e + 3.53d-7*e*e*e
     $               - 2.97d-10*e*e*e*e

      elseif( e .ge. 840.0d0  .and.  e .le. 5000.0d0 )then
         sigpp_fit = 50.9d0 - 3.8d-3*e + 2.78d-7*e*e
     $               + 1.92d-15*e*e*e*e
      else
      sigpp_fit = 40.05d0
      endif

      sigpp_fit = sigpp_fit/10.0d0

      end

c=======================================================================
      function sigpn_fit(e)
c                                     coded by H. Iwase  (Aug. 25, 2011)
c                                  modified by A. Kohama (Sep. 12, 2011)
c                                  modified by A. Kohama (Dec. 08, 2011)
c                                  modified by A. Kohama (Dec. 20, 2011)
c                                  modified by A. Kohama (Feb. 04, 2014)
c     Returns the value of sigma_pn^total as a fn. of energy in lab.
c     see eq. (2)
c     in PRC.81.064603 C. A. Bertulani and C. De Conti (2010)
c     applicable for 10 MeV =< e =< 5 GeV
c     Note:  0 MeV < e < 10 MeV: Returns the value at 10 MeV
c            5 GeV < e         : Error   (up to v5.2)
c            5 GeV < e         : Returns the value at 5 GeV
c             Following to the suggestion by L. Sihver  (Feb. 04, 2014)
c     input: e (MeV)   output: sigpn_fit (fm^2)
c=======================================================================
      implicit none
      real*8   e,sigpn_fit

      if    ( e .lt. 10.0d0 )then
         sigpn_fit = 924.4888d0
      elseif( e .ge. 10.0d0 .and. e .lt. 300.0d0 )then
         sigpn_fit = 89.4d0 -2025.0d0/sqrt(e)
     $               + 19108.0d0/e -43535.0d0/e/e

      elseif( e .ge. 300.0d0  .and.  e .lt. 700.0d0 )then
         sigpn_fit = 14.2d0 + 5436.0d0/e + 3.72d-5*e*e - 7.55d-9*e*e*e

      elseif( e .ge. 700.0d0  .and.  e .le. 5000.0d0 )then
         sigpn_fit = 33.9d0 + 6.1d-3*e - 1.55d-6*e*e
     $               + 1.3d-10*e*e*e
      else
      sigpn_fit = 41.90d0
      endif

      sigpn_fit = sigpn_fit/10.0d0

      end

c=======================================================================
      function cutoff_p(e, ecut, width)
c                                     coded by A. Kohama (Dec. 06, 2011)
c     smooth cut-off for e > ecut
c     input : e (MeV), ecut (MeV), width (MeV)
c     output: cutoff (no-dim)
c=======================================================================
      implicit none
      real*8   e, ecut, width, cutoff_p

      if ( (e - ecut)/width .gt. 700d0) then
         cutoff_p = 0.0d0
         return
      end if
      cutoff_p = 1.0d0/(1.0d0 + dexp((e - ecut)/width))

      end
c=======================================================================
      function cutoff_m(e, ecut, width)
c                                     coded by A. Kohama (Dec. 06, 2011)
c     smooth cut-off for e < ecut
c     input : e (MeV), ecut (MeV), width (MeV)
c     output: cutoff (no-dim)
c=======================================================================
      implicit none
      real*8   e, ecut, width, cutoff_m

      cutoff_m = 1.0d0/(1.0d0 + dexp((-e + ecut)/width))

      end
************************************************************************
