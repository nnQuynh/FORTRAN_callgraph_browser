!-----------------------------------------------------------------------
      double precision function sigma_d( x)
!     Deuteron disintegration cross section below photon energy of 250 MeV
!     x: incident photon energy (MeV)
!     returns cross section (b)
!     S.Noda 08/08/2013
!-----------------------------------------------------------------------
      implicit none
      double precision x

cABE change @2014/07/28
cABE      if ( x > 2.224 ) then
      if ( 250.0d0 .ge. x .and. x .gt. 2.224d0 ) then
cABE end
         sigma_d = 61.2*( x- 2.224)**1.5/ x**3/ 1e3

      else
         sigma_d = 0.

      end if

      return
      end function


!-----------------------------------------------------------------------
      double precision function pbf( x)
!     Pauli blocking function
!     x: incident photon energy (MeV)
!     S.Noda 08/08/2013
!-----------------------------------------------------------------------
      implicit none
      double precision x

      if ( x < 17 ) then
         pbf = 0.

      else if ( x > 150. ) then
         pbf = 0.88665462

      else
         pbf = 8.3714e-2- 9.8343e-3* x+ 4.1222e-4* x**2- 3.4762e-6* x**3
     $        +9.3537e-9* x**4

      end if

      return
      end function


!-----------------------------------------------------------------------
      double precision function sigma_qd( Z, N, x)
!     Quasideuteron disintegration cross section
!     x: incident photon energy (MeV)
!     returns cross section (b)
!     S.Noda 08/08/2013
!-----------------------------------------------------------------------
      implicit none
      double precision x
      double precision sigma_d, pbf         ! functions
      double precision L(100)               ! Levinger constant
      integer A, Z, N           ! Mass, proton & neutron number
      integer eval(100)         ! Nucleus with evaluated data
      integer i, k              ! counter

cABE add @2014/05/28, for bisection method
      integer ilo, ihi, middle

      ilo = 0
      ihi = 0
      middle = 0
cABE end

      k = 0
      A = Z + N

      ! Nucleus with evaluated data (without 25055)
      data ( eval(i), i=1, 67) /1002, 2003, 3006, 3007, 4009, 5010, 5011
     $     , 6012, 7014, 8016, 9019, 11023, 12024, 12025, 12026, 13027,
     $     14028, 14029, 14030, 15031, 20040, 20048, 22046, 23051 ,
     $     24052, 26054, 26056 , 27059, 28058, 28060, 29063, 29065,
     $     30064, 40090, 41093, 42092, 42094, 42096, 42098, 42100 ,
     $     55133, 64152, 64154, 64155 , 64156, 64157, 64158, 64160,
     $     73181, 74182, 74184, 74186, 79197, 80196, 80198, 80199, 80200
     $     , 80201, 80202, 80204, 82206 , 82207, 82208, 83209, 92235,
     $     92238, 93237/

      ! Levinger constant obtained from fitting JENDL/PD-2004
      data ( L(i), i=1, 67) /3.1984, 4.59619, 9.58239, 8.12356, 7.36906,
     $     7.89692, 9.10849, 10.2498, 10.1789, 10.1721, 9.46442, 7.46377
     $     , 13.9527, 8.51766, 5.21752, 10.1772, 6.86479, 2.8598,
     $     4.19423, 7.64929, 15.0024, 7.89163, 8.20111, 17.2354, 8.14885
     $     , 4.33232, 6.52553, 7.58507, 6.04506, 8.12711, 7.05181,
     $     7.30983, 6.47168, 7.38043, 7.90619, 8.46574, 8.6036, 7.89508,
     $     7.62862, 7.81672, 7.72712, 7.6744, 7.67794, 7.67611, 7.67451,
     $     7.67373, 7.67043, 7.6643, 7.36136, 6.50001, 6.49997, 6.56219,
     $     7.16243, 7.6513, 7.64892, 7.64598, 7.64664, 7.6432, 7.64451,
     $     7.64035, 6.91661, 6.87741, 6.89187, 7.66061, 9.15804, 9.21951
     $     , 8.86247/

      ! Levinger constant averaged L(i) for the nucleus without evaluation
      L(99) = 7.87779

      ! find nucleus
cABE change @2014/05/28, bisection method
      ilo = 1
      ihi = 67

      if (eval(ilo) .gt. Z*1000+A .or.
     &    eval(ihi) .lt. Z*1000+A) then
         k = 0
         goto 110
      endif

  100 middle = (ilo + ihi) / 2

      if (eval(middle) .eq. Z*1000+A)then
         k = middle
         goto 110
      elseif (eval(middle) .lt. Z*1000+A) then
         if (ilo .eq. middle) then
            k = 0
            goto 110
         endif
         ilo = middle
         goto 100
      else
         if (ihi .eq. middle) then
            k = 0
            goto 110
         endif
         ihi = middle
         goto 100
      endif

  110 continue
cABE end

      if ( k /= 0 ) then
         L(98) = L(k)
      else if ( k == 0 ) then
         L(98) = L(99)
      end if

      sigma_qd = L(98)/ A* N* Z* sigma_d(x)* pbf(x)

      return
      end function


************************************************************************
* Quasideuteron disintegration reaction
      subroutine qdeuteron
* Done by S. Noda
* Last revised: 02/18/2014
************************************************************************

*-----------------------------------------------------------------------
      implicit real*8(a-h,o-z)

      include 'param00.inc'
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      double precision r, p
      common /coodrp/ r(5,nnn), p(6,nnn)
!$OMP THREADPRIVATE(/coodrp/)

      common /gg005/ xxx, yyy, zzz, uuu, vvv, www, tme, erg,
     &               dls, wgt, vel, dtc,
     &               icl, iii, jjj, kkk, jsu, iap, jgp, ipt,
     &               mtp, iexp, iex, idx,
     &               npa, ncp
!$OMP THREADPRIVATE(/gg005/)

*-----------------------------------------------------------------------
cABE add @2014/06/13, flag for legacy JQMD or JQMD Ver.2~
      common /qmdflg/ irqmd

cABE add @2015/04/10
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)
!$OMP THREADPRIVATE(/coodid/)

*-----------------------------------------------------------------------
cABE end

      double precision u, v, w  ! unit vector of incident photon
      double precision e        ! energy of incident photon (MeV)
      double precision pc(5,2)  ! p in center-of-mass system
      double precision pp(5)  ! momentum and energy of incident photon
      double precision beta(3), gamma ! beta & gamma for Lorentz transform
      integer pid, nid        ! proton and neutron id to absorb photon
      double precision ecm(3) ! energy of CM system
      double precision pcm    ! momentum in CM
      double precision theta, phi ! for Lorentz transform
      integer k, l                ! counters
      double precision rand       ! random number
      double precision unirn      ! defined function
cABE add @2014/06/13, to calculate (mass^2 + 2*mass*epot)
      double precision epotp, epotn       ! potential energy for nucleon
cABE end

      u = uuu
      v = vvv
      w = www
      k = 0
      l = 0

*-----------------------------------------------------------------------
!     Obtain deuteron-like proton & neutron pair
      call qdpair( pid, nid, k)

      e = erg/ 1e3
      pp(1) = u* e
      pp(2) = v* e
      pp(3) = w* e
      pp(4) = e
      pp(5) = 0.

!     ---------- ---------- ---------- ---------- ----------
!     (1) Laboratory system (Lab) ==> Center of Mass system (CM)
!     pc(x,1): proton in CM
!     pc(x,2): neutron in CM
!     beta(1): parameter beta for x-axis
!     beta(2): parameter beta for y-axis
!     beta(3): parameter beta for z-axis
!     gamma: parameter gamma

cABE change @2014/06/13, p(4,*) is changed energy -> energy(pot included) in JQMD Ver.2~
      if(irqmd .eq. 0) then
         beta(1) = (p(1,pid) + p(1,nid) + pp(1))
     &             / (p(4,pid) + p(4,nid) + pp(4))
         beta(2) = (p(2,pid) + p(2,nid) + pp(2))
     &             / (p(4,pid) + p(4,nid) + pp(4))
         beta(3) = (p(3,pid) + p(3,nid)+ pp(3))
     &             / (p(4,pid) + p(4,nid) + pp(4))
      else
         beta(1) = (p(1,pid) + p(1,nid) + pp(1))
     &             / (p(6,pid) + p(6,nid) + pp(4))
         beta(2) = (p(2,pid) + p(2,nid) + pp(2))
     &             / (p(6,pid) + p(6,nid) + pp(4))
         beta(3) = (p(3,pid) + p(3,nid)+ pp(3))
     &             / (p(6,pid) + p(6,nid) + pp(4))
      endif
cABE end

      gamma = 1./ sqrt( 1.- beta(1)* beta(1)- beta(2)* beta(2)- beta(3)*
     $     beta(3))

!     (2) in CM, quasideuteron absorbs photon
!     ecm(1): energy of center-of-mass system
cABE change @2014/06/13, p(4,*) is changed energy -> energy(pot included) in JQMD Ver.2~
      if (irqmd .eq. 0) then
         ecm(1) = sqrt((p(4,pid) + p(4,nid) + pp(4))**2
     &                 - (p(1,pid) + p(1,nid) + pp(1))**2
     &                 - (p(2,pid) + p(2,nid) + pp(2))**2
     &                 - (p(3,pid) + p(3,nid) + pp(3))**2)
      else
         ecm(1) = sqrt((p(6,pid) + p(6,nid) + pp(4))**2
     &                 - (p(1,pid) + p(1,nid) + pp(1))**2
     &                 - (p(2,pid) + p(2,nid) + pp(2))**2
     &                 - (p(3,pid) + p(3,nid) + pp(3))**2)
      endif
cABE end

!     ecm(2), ecm(3): energies of proton & neutron in CM
      ecm(2) = 1./ 2.*( ecm(1)+( p(5,pid)**2- p(5,nid)**2)/ ecm(1))
      ecm(3) = 1./ 2.*( ecm(1)-( p(5,pid)**2- p(5,nid)**2)/ ecm(1))

!     absolute value of proton & neutron's momentum in CM
      pcm = sqrt( ecm(2)**2- p(5,pid)**2)

!     (3) CM ==> Lab
      rand = unirn()
      theta = pi* rand
      rand = unirn()
      phi = 2.* pi* rand

      pc(1,1) = pcm* sin( theta)* cos( phi)
      pc(2,1) = pcm* sin( theta)* sin( phi)
      pc(3,1) = pcm* cos( theta)
      pc(4,1) = ecm(2)

      pc(1,2) = -1.* pc(1,1)
      pc(2,2) = -1.* pc(2,1)
      pc(3,2) = -1.* pc(3,1)
      pc(4,2) = ecm(3)

!     Lorentz transform
!     The value below may be slightly different from the one next below.
c$$$      p(4,pid) = gamma* pc(4,1)+ beta(1)* gamma* pc(1,1)+ beta(2)* gamma
c$$$     $     * pc(2,1)+ beta(3)* gamma* pc(3,1)
      p(1,pid) = beta(1)* gamma* pc(4,1)+ gamma* pc(1,1)
      p(2,pid) = beta(2)* gamma* pc(4,1)+ gamma* pc(2,1)
      p(3,pid) = beta(3)* gamma* pc(4,1)+ gamma* pc(3,1)
cABE change @2014/06/13, p(4,*) is changed energy -> energy(pot included)
      if (irqmd .eq. 0) then
         p(4,pid) = sqrt(p(1,pid)**2 + p(2,pid)**2 + p(3,pid)**2
     &                   + p(5,pid)**2)
      else
         epotp = (p(4,pid)**2 - p(1,pid)**2 - p(2,pid)**2 - p(3,pid)**2
     &            - p(5,pid)**2) * 0.5d0 / p(5,pid)
         p(4,pid) = sqrt(p(1,pid)**2 + p(2,pid)**2 + p(3,pid)**2
     &                   + p(5,pid)**2 + 2.0d0 * p(5,pid) * epotp)
         p(6,pid) = sqrt(p(1,pid)**2 + p(2,pid)**2 + p(3,pid)**2
     &                   + p(5,pid)**2)
      endif
cABE end

c$$$      p(4,nid) = gamma* pc(4,2)+ beta(1)* gamma* pc(1,2)+ beta(2)* gamma
c$$$     $     * pc(2,2)+ beta(3)* gamma* pc(3,2)
      p(1,nid) = beta(1)* gamma* pc(4,2)+ gamma* pc(1,2)
      p(2,nid) = beta(2)* gamma* pc(4,2)+ gamma* pc(2,2)
      p(3,nid) = beta(3)* gamma* pc(4,2)+ gamma* pc(3,2)
cABE change @2014/06/13, p(4,*) is changed energy -> energy(pot included)
      if (irqmd .eq. 0) then
         p(4,nid) = sqrt(p(1,nid)**2 + p(2,nid)**2 + p(3,nid)**2
     &                   + p(5,nid)**2)
      else
         epotn = (p(4,nid)**2 - p(1,nid)**2 - p(2,nid)**2 - p(3,nid)**2
     &            - p(5,nid)**2) * 0.5d0 / p(5,nid)
         p(4,nid) = sqrt(p(1,nid)**2 + p(2,nid)**2 + p(3,nid)**2
     &                   + p(5,nid)**2 + 2.0d0 * p(5,nid) * epotn)
         p(6,nid) = sqrt(p(1,nid)**2 + p(2,nid)**2 + p(3,nid)**2
     &                   + p(5,nid)**2)
      endif
cABE end

cABE add @2015/04/10
      iavd(pid) = -1
      iavd(nid) = -1

!     Quasideuteron has been made.

      return
      end subroutine qdeuteron


!------------------------------------------------------------------------
      double precision function phxs( Z, N, E)
! returns cross section of photonuclear reaction
! As of Oct 2013, GDR, QD and pion production reactions.
! Z: atomic number
! N: neutron number
! E: incident photon energy (MeV)
! by S.Noda 10/25/2013
! Last modified by T.Ogawa 08/25/2014
!------------------------------------------------------------------------
      implicit none
      common /pnint/  ipnint
      integer Z, N, ipnint
      double precision E
      double precision getGDRxsec, sigma_qd, sigma_delta, sigma_nstar,
     & sigmaNRF
                                ! functions

      if ( E .le. 0.01d0 .or. (E .le. 0.5d0 .and. ipnint .le. 1) ) then        ! ~10 keV phxs = 0.  2014/8/25 ogawa.  For NRF
         phxs = 0.d0
      else
cABE 2016/03/09
         phxs = getGDRxsec(Z, N, E)+ sigma_qd(Z, N, E)
     $        + sigma_delta(dble(Z+N),dble(E),1,5)
     $        + sigmaNRF(Z, N, E)
cABE end
      endif



      return
      end function phxs


!------------------------------------------------------------------------
      integer function jgdrqd( Z, N, E)
! judges photonuclear reaction (giant-diple resonance, quasideuteron or
! pion production) according to the ratio of cross sections.
! jgdrqd =  2 ==> NRF
! jgdrqd =  1 ==> GDR
! jgdrqd =  0 ==> QD
! jgdrqd = -1 ==> PD
! Z: atomic number
! N: neutron number
! E: incident photon energy (MeV)
! by S.Noda 10/25/2013
!------------------------------------------------------------------------
      implicit none
      integer Z, N              ! atomic number, neutron number
      double precision E        ! photon energy (MeV)
      double precision xsnrf, xsgdr, xsqd, xspion ! cross sections of NRF, QDR, QD
                                           ! PD
      double precision sigmaNRF, getGDRxsec, sigma_qd, sigma_delta,
     & sigma_nstar
                                ! functions, cross sections of NRF, QDR, QD and PD
      double precision ratio0, ratio1, ratio2 ! ratios of cross sections
      double precision random   ! random number
      double precision unirn    ! function, random number
cABE add @2014/07/22
      double precision xsstring, xsptot, ratio3
cABE end

      xsnrf = sigmaNRF( Z, N, E)  ! 2014/8/25 ogawa
      xsgdr = getGDRxsec( Z, N, E)
      xsqd = sigma_qd( Z, N, E)
      xspion = sigma_delta( dble(Z+N), dble(E),1,4)
cABE add @2014/07/30, to use threshold energy of string production
      xsstring = sigma_delta( dble(Z+N), dble(E),5,5)
      xsptot = xsnrf + xsgdr + xsqd + xspion + xsstring

      if (xsptot .le. 0.0d0) then
         jgdrqd = 10
         return
      endif
cABE end

cABE change @2014/07/30
      ratio0 = xsnrf  / xsptot            ! 2014/8/25 ogawa
      ratio1 = xsgdr  / xsptot + ratio0
      ratio2 = xsqd   / xsptot + ratio1
      ratio3 = xspion / xsptot + ratio2
cABE end

      random = unirn()

      ! 0--<NRF>--ratio0--<GDR>--ratio1--<QD>--ratio2--<PD>--1
      if ( random .le. ratio0 .and. Z .ge. 3 .and. N .ge. 3) then ! NRF
         jgdrqd = 2
      else if ( ratio0 .lt. random .and. random .le. ratio1 ) then ! GDR
         jgdrqd = 1
cABE add @2014/05/09
         if (Z .eq. 1 .and. N .eq. 1) jgdrqd = 0
cABE end

      else if ( ratio1 .lt. random .and. random .le. ratio2 ) then ! QD
         jgdrqd = 0

cABE change @ 2014/07/22
      elseif (ratio2 .lt. random .and. random .le. ratio3) then         ! PD
         jgdrqd = -1
      else
         jgdrqd = -2
cABE end
      end if

      return
      end function jgdrqd

!-------------------------------------------------------------------------------
      subroutine qdpair( pid, nid, k)
!     deuteron-like proton neutron pairing
!     pid: proton id
!     nid: neutron id
!     k: counter, if k>20 then go back to re-make the ground state
!-------------------------------------------------------------------------------
c$$$      include 'param00.inc'
c$$$      include 'param01.inc'
      integer idnta, idnpr, massta, masspr, mstapr, msprpr
      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr
!$OMP THREADPRIVATE(/const0/)
c$$$      double precision r, p
c$$$      common /coodrp/ r(5,nnn), p(6,nnn)
c$$$!$OMP THREADPRIVATE(/coodrp/)

!     implicit none
c$$$      double precision integ    ! integral value of wave function
c$$$      double precision x, dx    ! integral inderval
c$$$      double precision rand     ! random number
c$$$      double precision deltap   ! difference of p & n's moemntums
c$$$      double precision proto    ! proportional part \sigma_\mathrm{qd}
c$$$      double precision, allocatable :: qdr(:,:) ! store p-n distance
      integer pid, nid            ! index of proton & neutron
      integer i, j, k             ! counters
c$$$      integer d(2)                ! receiver
c$$$      double precision u2         ! defined function
c$$$      double precision sigmaproto ! defined function
      double precision unirn      ! defined function

c$$$      allocate( qdr( mstapr, massta- mstapr))
      pid = int( unirn()* mstapr)+ 1
      nid = int( unirn()*( massta- mstapr))+ mstapr+ 1

c$$$C Choose deutron-like pair carefully
c$$$      pid = 0
c$$$      nid = 0
c$$$
c$$$!     -------- -------- -------- -------- --------
c$$$!     20 tries to determine the pair
c$$$!     (plus 10 tries of remaking ground state)
c$$$!     -------- -------- -------- -------- --------
c$$$      do k = 1, 20
c$$$
c$$$!     -------- -------- -------- -------- --------
c$$$!     (1) distance
c$$$!     -------- -------- -------- -------- --------
c$$$!     (1.1) Find x (fm), corresponding distance between p & n
c$$$         x = 0.
c$$$         dx = 0.01
c$$$         integ = 0.
c$$$         d(1:2) = 0
c$$$         rand = unirn()
c$$$
c$$$         do while (integ < rand)
c$$$            x = x+ dx
c$$$            integ = integ+ dx*( u2(x)+ u2(x-dx))/ 2.
c$$$         end do
c$$$
c$$$!     (1.2) Find good p-n pair of smallest abs( qdr(i,j)- x)
c$$$!     qdr(i,j): distances of all p-n pairs in nucleus
c$$$!     x: deuteron-like distance
c$$$         qdr(1:mstapr,1:massta-mstapr) = 0.
c$$$         do i = 1, mstapr           ! i: proton
c$$$            do j = mstapr+1, massta ! j: neutron
c$$$               qdr(i,j-mstapr) = sqrt(( r(1,i)- r(1,j))**2+( r(2,i)- r(2
c$$$     $              ,j))**2+( r(3,i)- r(3,j))**2)
c$$$               qdr(i,j-mstapr) = abs( qdr(i,j-mstapr)- x)
c$$$            end do
c$$$         end do
c$$$
c$$$!     -------- -------- -------- -------- --------
c$$$!     (2) momenta
c$$$!     -------- -------- -------- -------- --------
c$$$         proto = 0.
c$$$         rand = unirn()
c$$$         d(1:2) = minloc( qdr)
c$$$         pid = d(1)
c$$$         nid = d(2)+ mstapr
c$$$
c$$$         deltap = sqrt(( p(1,pid)- p(1,nid))**2+ ( p(2,pid)- p(2,nid))
c$$$     $        **2+( p(3,pid)- p(3,nid))**2+( p(5,pid)- p(5,nid))**2)
c$$$
c$$$         proto = sigmaproto( deltap, massta)
c$$$
c$$$         if ( rand <= proto) then
c$$$            exit
c$$$         else
c$$$            qdr(d(1),d(2)) = 999.+ dble(i) ! Move to next resonable
c$$$                                           ! pair of p-n
c$$$         end if
c$$$      end do

      end subroutine qdpair


!*******************************************************************************
      double precision function u2(x)
!     Hulthen wave function
!     x: distance between proton & neutron in fm
!*******************************************************************************
      implicit none
      double precision x
      double precision N2, gamma, beta
      N2 = 0.783
      gamma = 0.2316
      beta = 5.98* gamma

      u2 = ( sqrt( N2)*( exp( -1.* gamma* x)- exp( -1.* beta* x)))**2

      return
      end function


!*******************************************************************************
      double precision function sigmaproto(x, A)
!     Levinger, Phys Rev 84 (1951).
!     \sigma_\mathrm{qd}( \Delta{}p) \propto \frac{1}{\left(\alpha^2+
!     \frac{1}{4\hbar^2}\Delta{}p^2\right) V}
!     \Delta{}p = |p_{p}- p_{n}|
!     x: difference of p& n's momentums (GeV/c)
!     A: mass number
!*******************************************************************************
      implicit none
      double precision x
      integer A
      double precision alpha    ! scattering length (fm^{-1})
      double precision hbar
      double precision V        ! volume
      double precision c        ! speed of light
      double precision pi
      double precision integral_sigmaproto

      alpha = 0.232             ! fm^{-1}
      hbar = 6.58211928e-16     ! eV s
      c = 2.99792458e23         ! fm/ s
      pi = 3.14159265
      integral_sigmaproto = 64.1699

      V = 4.* pi/ 3.* A*( 1.2e-1)**3 ! fm^{3}
      sigmaproto = 1./(( alpha**2+ 1./( 4.* hbar**2)*( x* 1e9/ c)**2)*
     $     V)/ integral_sigmaproto

      return
      end function


! ----------------------------------------------------------------------
!     Modify the decay width of Li-6, B-10, C-12, N-14 and O-16 in GDR
      subroutine suppressalpha( iz, in, u)
!     iz: atomic number
!     in: neutron number
!     u: incident photon energy (MeV)
! ----------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

      integer iz, in
      double precision u

! Li-6
      if ( iz == 3 .and. in == 3 ) then
         call suppressalpha_li6( u)
      end if

! B-10
      if ( iz == 5 .and. in == 5 ) then
         call suppressalpha_b10( u)
      end if

! C-12
      if ( iz == 6 .and. in == 6 ) then
         call suppressalpha_c12( u)
      end if

! N-14
      if ( iz == 7 .and. in == 7 ) then
         call suppressalpha_n14( u)
      end if

! O-16
      if ( iz == 8 .and. in == 8 ) then
         call suppressalpha_o16( u)
      end if

      end subroutine suppressalpha


************************************************************************
*                                                                      *
      double precision function getGDRxsec(iZ,iN,Ene)
*                                                                      *
*     iZ: atomic number                                                *
*     iN: neutron number                                               *
*     Ene: incident photon energy (MeV)                                *
*     getGDRxsec returns giant-dipole resonace cross section (barn)    *
*     last modified by S. Noda, February 22 2013.                      *
*                                                                      *
************************************************************************

      use NGSDATAMOD, only: shellE
      implicit none
      integer iZ, iN, iA
      double precision Ene      ! in MeV
      double precision Enev     ! in eV
      double precision p0, p1, p2, p3   ! Loretzian parameters
      integer i, j, k           ! counters

!     Parameters in Lorentz distribution
      double precision center, width, height

      integer nbins(100)        ! number of bins
      double precision eg(100,1000)         ! photon energy (eV)
      double precision gdr(100,1000)        ! cross section (b)
      double precision sigma_qd                ! function, cross
                                               ! section of
                                               ! quasideuteron
      integer middle, ilo, ihi

C     Nucleus with evaluated data (without 25055)
      integer eval(100)
      data ( eval(i), i=1, 67) /1002, 2003, 3006, 3007, 4009, 5010, 5011
     $     , 6012, 7014, 8016, 9019, 11023, 12024, 12025, 12026, 13027,
     $     14028, 14029, 14030, 15031, 20040, 20048, 22046, 23051 ,
     $     24052, 26054, 26056 , 27059, 28058, 28060, 29063, 29065,
     $     30064, 40090, 41093, 42092, 42094, 42096, 42098, 42100 ,
     $     55133, 64152, 64154, 64155 , 64156, 64157, 64158, 64160,
     $     73181, 74182, 74184, 74186, 79197, 80196, 80198, 80199, 80200
     $     , 80201, 80202, 80204, 82206 , 82207, 82208, 83209, 92235,
     $     92238, 93237/

C     1002, JENDL evaluated data
      nbins(1) = 28
      data ( eg(1,i), i=1, 28) /2224000.0, 2300000.0, 2600000.0,
     $     3000000.0, 3500000.0, 4000000.0, 4500000.0, 5000000.0,
     $     5500000.0, 6000000.0, 8000000.0, 1.0D7, 1.5D7, 2.0D7, 2.5D7,
     $     3.0D7, 3.5D7, 4.0D7, 4.5D7, 5.0D7, 5.5D7, 6.0D7, 7.0D7, 8.0D7
     $     , 9.0D7, 1.0D8, 1.2D8, 1.4D8/

      data ( gdr(1,i), i=1, 28) /0.0, 6.581D-4, 0.001183, 0.001828,
     $     0.00228, 0.002454, 0.002474, 0.002416, 0.002318, 0.002205,
     $     0.001753, 0.001373, 8.528D-4, 5.882D-4, 4.35D-4, 3.376D-4,
     $     2.714D-4, 2.242D-4, 1.893D-4, 1.629D-4, 1.425D-4, 1.264D-4,
     $     1.031D-4, 8.744D-5, 7.622D-5, 6.764D-5, 5.483D-5, 4.558D-5/


C     2003, JENDL evaluated data
      nbins(2) = 43
      data ( eg(2,i), i=1, 43) /5493500.0, 6000000.0, 6500000.0,
     $     7000000.0, 7500000.0, 8000000.0, 8500000.0, 9000000.0,
     $     9500000.0, 1.0D7, 1.05D7, 1.1D7, 1.15D7, 1.2D7, 1.25D7, 1.3D7
     $     , 1.35D7, 1.4D7, 1.5D7, 1.6D7, 1.7D7, 1.8D7, 1.9D7, 2.0D7,
     $     2.2D7, 2.4D7, 2.6D7, 2.8D7, 3.0D7, 3.5D7, 4.0D7, 4.5D7, 5.0D7
     $     , 6.0D7, 7.0D7, 8.0D7, 9.0D7, 1.0D8, 1.1D8, 1.2D8, 1.3D8,
     $     1.4D8, 1.5D8/

      data ( gdr(2,i), i=1, 43) /0.0, 4.1254D-4, 5.6976D-4, 6.802D-4,
     $     7.6218D-4, 8.2666D-4, 8.7832D-4, 9.2248D-4, 0.0010895,
     $     0.0013577, 0.0015307, 0.0016697, 0.0017784, 0.0018637,
     $     0.0019314, 0.0019844, 0.0019902, 0.0019872, 0.0019676,
     $     0.001938, 0.001885, 0.0018331, 0.0017789, 0.0016873,
     $     0.0014814, 0.0012963, 0.0011391, 0.0010021, 8.824D-4, 6.6608D
     $     -4, 5.2704D-4, 4.3109D-4, 3.6159D-4, 2.696D-4, 2.1337D-4,
     $     1.7665D-4, 1.5114D-4, 1.3226D-4, 1.1739D-4, 1.0516D-4,
     $     9.4862D-5, 8.611D-5, 7.8662D-5/


C     3006, JENDL evaluated data
      nbins(3) = 215
      data ( eg(3,i), i=1, 215) /3697000.0, 3800000.0, 4000000.0,
     $     4200000.0, 4400000.0, 4600000.0, 4800000.0, 5000000.0,
     $     5200000.0, 5400000.0, 5600000.0, 5800000.0, 6000000.0,
     $     6200000.0, 6400000.0, 6600000.0, 6800000.0, 7000000.0,
     $     7200000.0, 7400000.0, 7600000.0, 7800000.0, 8000000.0,
     $     8200000.0, 8400000.0, 8600000.0, 8800000.0, 9000000.0,
     $     9200000.0, 9400000.0, 9600000.0, 9800000.0, 1.0D7, 1.02D7,
     $     1.04D7, 1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7, 1.18D7
     $     , 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7, 1.32D7,
     $     1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7, 1.48D7
     $     , 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7, 1.62D7,
     $     1.64D7, 1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7, 1.78D7
     $     , 1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7, 1.92D7,
     $     1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.05D7, 2.1D7, 2.15D7, 2.2D7,
     $     2.25D7, 2.3D7, 2.35D7, 2.4D7, 2.45D7, 2.5D7, 2.55D7, 2.6D7,
     $     2.65D7, 2.7D7, 2.75D7, 2.8D7, 2.85D7, 2.9D7, 2.95D7, 3.0D7,
     $     3.05D7, 3.1D7, 3.15D7, 3.2D7, 3.25D7, 3.3D7, 3.35D7, 3.4D7,
     $     3.45D7, 3.5D7, 3.55D7, 3.6D7, 3.65D7, 3.7D7, 3.75D7, 3.8D7,
     $     3.85D7, 3.9D7, 3.95D7, 4.0D7, 4.05D7, 4.1D7, 4.15D7, 4.2D7,
     $     4.25D7, 4.3D7, 4.35D7, 4.4D7, 4.45D7, 4.5D7, 4.55D7, 4.6D7,
     $     4.65D7, 4.7D7, 4.75D7, 4.8D7, 4.85D7, 4.9D7, 4.95D7, 5.0D7,
     $     5.05D7, 5.1D7, 5.15D7, 5.2D7, 5.25D7, 5.3D7, 5.35D7, 5.4D7,
     $     5.45D7, 5.5D7, 5.55D7, 5.6D7, 5.65D7, 5.7D7, 5.75D7, 5.8D7,
     $     5.85D7, 5.9D7, 5.95D7, 6.0D7, 6.05D7, 6.1D7, 6.15D7, 6.2D7,
     $     6.25D7, 6.3D7, 6.35D7, 6.4D7, 6.45D7, 6.5D7, 6.55D7, 6.6D7,
     $     6.65D7, 6.7D7, 6.75D7, 6.8D7, 6.85D7, 6.9D7, 6.95D7, 7.0D7,
     $     7.05D7, 7.1D7, 7.15D7, 7.2D7, 7.25D7, 7.3D7, 7.35D7, 7.4D7,
     $     7.45D7, 7.5D7, 7.55D7, 7.6D7, 7.65D7, 7.7D7, 7.75D7, 7.8D7,
     $     7.85D7, 7.9D7, 7.95D7, 8.0D7, 8.5D7, 9.0D7, 9.5D7, 1.0D8,
     $     1.05D8, 1.1D8, 1.15D8, 1.2D8, 1.25D8, 1.3D8, 1.35D8, 1.4D8/

      data ( gdr(3,i), i=1, 215) /0.0, 2.5859D-5, 1.504D-4, 2.6764D-4,
     $     3.9646D-4, 5.3899D-4, 6.8865D-4, 8.4888D-4, 0.0010303,
     $     0.001229, 0.0014373, 0.0016281, 0.0017944, 0.0019468,
     $     0.0020833, 0.0022059, 0.002319, 0.0024275, 0.0025347,
     $     0.0026413, 0.0027456, 0.0028361, 0.0029107, 0.0029752,
     $     0.0030275, 0.0030672, 0.0030952, 0.0031136, 0.0031255,
     $     0.0031344, 0.0031432, 0.0031506, 0.0031586, 0.0031715,
     $     0.0031895, 0.0032118, 0.003237, 0.0032631, 0.0032878,
     $     0.0033087, 0.0033232, 0.0033275, 0.0033206, 0.0033034,
     $     0.0032758, 0.0032385, 0.0031926, 0.0031398, 0.0030823,
     $     0.0030221, 0.0029616, 0.0029016, 0.0028441, 0.0027923,
     $     0.0027482, 0.002714, 0.0026919, 0.0026842, 0.0026929,
     $     0.0027182, 0.0027561, 0.0027947, 0.0028136, 0.0027938,
     $     0.0027294, 0.0026373, 0.0025382, 0.0024466, 0.0023685,
     $     0.0023047, 0.0022536, 0.002213, 0.0021808, 0.0021554,
     $     0.0021354, 0.0021198, 0.002107, 0.0020965, 0.0020884,
     $     0.0020823, 0.0020778, 0.0020746, 0.0020725, 0.0020708,
     $     0.002072, 0.0020733, 0.0020737, 0.002073, 0.0020703,
     $     0.0020652, 0.0020566, 0.0020447, 0.0020298, 0.0020118,
     $     0.0019908, 0.0019668, 0.0019402, 0.0019114, 0.0018807,
     $     0.0018483, 0.0018145, 0.0017796, 0.001744, 0.0017083,
     $     0.0016725, 0.0016368, 0.0016014, 0.0015666, 0.0015323,
     $     0.0014988, 0.001466, 0.0014342, 0.0014033, 0.0013734,
     $     0.0013445, 0.0013167, 0.0012898, 0.001264, 0.0012393,
     $     0.0012155, 0.0011928, 0.0011709, 0.00115, 0.0011301, 0.001111
     $     , 0.0010927, 0.0010752, 0.0010585, 0.0010426, 0.0010274,
     $     0.0010128, 9.9884D-4, 9.855D-4, 9.7276D-4, 9.6059D-4, 9.4895D
     $     -4, 9.3783D-4, 9.272D-4, 9.1703D-4, 9.0731D-4, 8.98D-4,
     $     8.8909D-4, 8.8055D-4, 8.7237D-4, 8.6453D-4, 8.57D-4, 8.4978D
     $     -4, 8.4285D-4, 8.3618D-4, 8.2977D-4, 8.2354D-4, 8.1751D-4,
     $     8.1169D-4, 8.0608D-4, 8.0066D-4, 7.9543D-4, 7.9038D-4,
     $     7.8549D-4, 7.8076D-4, 7.7618D-4, 7.7174D-4, 7.6743D-4,
     $     7.6324D-4, 7.5918D-4, 7.5522D-4, 7.5137D-4, 7.4762D-4,
     $     7.4396D-4, 7.4039D-4, 7.3691D-4, 7.3346D-4, 7.3007D-4,
     $     7.2675D-4, 7.235D-4, 7.2031D-4, 7.1719D-4, 7.1411D-4, 7.1109D
     $     -4, 7.0812D-4, 7.052D-4, 7.0232D-4, 6.9948D-4, 6.9668D-4,
     $     6.9391D-4, 6.9118D-4, 6.8848D-4, 6.8581D-4, 6.8317D-4,
     $     6.8055D-4, 6.7796D-4, 6.7537D-4, 6.7279D-4, 6.7023D-4, 6.677D
     $     -4, 6.6517D-4, 6.6267D-4, 6.6018D-4, 6.5771D-4, 6.5525D-4,
     $     6.528D-4, 6.5036D-4, 6.4793D-4, 6.4552D-4, 6.2171D-4, 5.9823D
     $     -4, 5.8498D-4, 5.7004D-4, 5.5394D-4, 5.3717D-4, 5.2019D-4,
     $     5.0342D-4, 4.8725D-4, 4.7207D-4, 4.5824D-4, 4.4608D-4/


C     3007, JENDL evaluated data
      nbins(4) = 221
      data ( eg(4,i), i=1, 221) /2467000.0, 2600000.0, 2800000.0,
     $     3000000.0, 3200000.0, 3400000.0, 3600000.0, 3800000.0,
     $     4000000.0, 4200000.0, 4400000.0, 4600000.0, 4800000.0,
     $     5000000.0, 5200000.0, 5400000.0, 5600000.0, 5800000.0,
     $     6000000.0, 6200000.0, 6400000.0, 6600000.0, 6800000.0,
     $     7000000.0, 7200000.0, 7400000.0, 7600000.0, 7800000.0,
     $     8000000.0, 8200000.0, 8400000.0, 8600000.0, 8800000.0,
     $     9000000.0, 9200000.0, 9400000.0, 9600000.0, 9800000.0, 1.0D7,
     $     1.02D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7
     $     , 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7,
     $     1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7
     $     , 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7,
     $     1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7
     $     , 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7,
     $     1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.05D7, 2.1D7, 2.15D7,
     $     2.2D7, 2.25D7, 2.3D7, 2.35D7, 2.4D7, 2.45D7, 2.5D7, 2.55D7,
     $     2.6D7, 2.65D7, 2.7D7, 2.75D7, 2.8D7, 2.85D7, 2.9D7, 2.95D7,
     $     3.0D7, 3.05D7, 3.1D7, 3.15D7, 3.2D7, 3.25D7, 3.3D7, 3.35D7,
     $     3.4D7, 3.45D7, 3.5D7, 3.55D7, 3.6D7, 3.65D7, 3.7D7, 3.75D7,
     $     3.8D7, 3.85D7, 3.9D7, 3.95D7, 4.0D7, 4.05D7, 4.1D7, 4.15D7,
     $     4.2D7, 4.25D7, 4.3D7, 4.35D7, 4.4D7, 4.45D7, 4.5D7, 4.55D7,
     $     4.6D7, 4.65D7, 4.7D7, 4.75D7, 4.8D7, 4.85D7, 4.9D7, 4.95D7,
     $     5.0D7, 5.05D7, 5.1D7, 5.15D7, 5.2D7, 5.25D7, 5.3D7, 5.35D7,
     $     5.4D7, 5.45D7, 5.5D7, 5.55D7, 5.6D7, 5.65D7, 5.7D7, 5.75D7,
     $     5.8D7, 5.85D7, 5.9D7, 5.95D7, 6.0D7, 6.05D7, 6.1D7, 6.15D7,
     $     6.2D7, 6.25D7, 6.3D7, 6.35D7, 6.4D7, 6.45D7, 6.5D7, 6.55D7,
     $     6.6D7, 6.65D7, 6.7D7, 6.75D7, 6.8D7, 6.85D7, 6.9D7, 6.95D7,
     $     7.0D7, 7.05D7, 7.1D7, 7.15D7, 7.2D7, 7.25D7, 7.3D7, 7.35D7,
     $     7.4D7, 7.45D7, 7.5D7, 7.55D7, 7.6D7, 7.65D7, 7.7D7, 7.75D7,
     $     7.8D7, 7.85D7, 7.9D7, 7.95D7, 8.0D7, 8.5D7, 9.0D7, 9.5D7,
     $     1.0D8, 1.05D8, 1.1D8, 1.15D8, 1.2D8, 1.25D8, 1.3D8, 1.35D8,
     $     1.4D8/

      data ( gdr(4,i), i=1, 221) /0.0, 4.2774D-6, 1.5847D-5, 2.5942D-5,
     $     3.6124D-5, 4.6475D-5, 5.6139D-5, 6.5979D-5, 7.678D-5, 8.8596D
     $     -5, 1.0149D-4, 1.13D-4, 1.2408D-4, 1.3596D-4, 1.4874D-4,
     $     1.625D-4, 1.7735D-4, 1.9344D-4, 2.1096D-4, 2.3019D-4, 2.5154D
     $     -4, 2.7322D-4, 2.971D-4, 3.2588D-4, 3.6122D-4, 4.0102D-4,
     $     4.3295D-4, 4.5003D-4, 4.6541D-4, 4.8831D-4, 5.1969D-4,
     $     5.5698D-4, 6.0099D-4, 6.5344D-4, 7.1573D-4, 7.897D-4, 8.775D
     $     -4, 9.8077D-4, 0.0010929, 0.0011841, 0.001237, 0.0012848,
     $     0.0013543, 0.0014475, 0.0015576, 0.0016757, 0.0017922,
     $     0.001899, 0.0019924, 0.0020749, 0.0021537, 0.0022352,
     $     0.0023262, 0.0024302, 0.0025457, 0.0026679, 0.0027892,
     $     0.0029005, 0.002993, 0.00306, 0.0030987, 0.0031109, 0.0031027
     $     , 0.0030819, 0.0030559, 0.0030346, 0.0030244, 0.00303,
     $     0.0030546, 0.0031001, 0.003167, 0.0032548, 0.0033606,
     $     0.0034778, 0.0035946, 0.0036958, 0.0037692, 0.0038142,
     $     0.0038373, 0.0038448, 0.0038394, 0.0038211, 0.0037896,
     $     0.003745, 0.0036889, 0.0036237, 0.0035527, 0.0034797,
     $     0.0034082, 0.0032558, 0.0031743, 0.0031794, 0.0032483,
     $     0.0033068, 0.0032609, 0.0030938, 0.0028803, 0.002704,
     $     0.0026064, 0.0025968, 0.0026768, 0.0028376, 0.0030467,
     $     0.0032286, 0.0032879, 0.0032, 0.003035, 0.0028758, 0.0027616,
     $     0.0026964, 0.0026724, 0.0026783, 0.0027041, 0.0027411,
     $     0.0027818, 0.0028209, 0.0028534, 0.0028751, 0.002883,
     $     0.0028751, 0.0028507, 0.0028104, 0.0027553, 0.0026877,
     $     0.0026099, 0.0025249, 0.0024354, 0.0023438, 0.0022524,
     $     0.0021628, 0.0020767, 0.001995, 0.0019186, 0.0018479,
     $     0.0017829, 0.0017242, 0.0016719, 0.0016261, 0.0015868,
     $     0.001554, 0.0015277, 0.0015079, 0.0014946, 0.0014878,
     $     0.0014873, 0.0014928, 0.0015037, 0.0015189, 0.0015367,
     $     0.0015549, 0.0015706, 0.0015806, 0.0015817, 0.0015717,
     $     0.0015494, 0.0015156, 0.0014723, 0.0014223, 0.0013686,
     $     0.0013139, 0.0012602, 0.0012092, 0.0011615, 0.0011176,
     $     0.0010776, 0.0010413, 0.0010085, 9.7897D-4, 9.5231D-4,
     $     9.2824D-4, 9.0647D-4, 8.8673D-4, 8.688D-4, 8.5244D-4, 8.3735D
     $     -4, 8.2351D-4, 8.1079D-4, 7.9907D-4, 7.8823D-4, 7.7818D-4,
     $     7.6884D-4, 7.6013D-4, 7.5198D-4, 7.4435D-4, 7.3718D-4,
     $     7.3043D-4, 7.2405D-4, 7.1802D-4, 7.1229D-4, 7.0685D-4,
     $     7.0167D-4, 6.9672D-4, 6.9199D-4, 6.8745D-4, 6.8303D-4,
     $     6.7879D-4, 6.747D-4, 6.7075D-4, 6.6693D-4, 6.6324D-4, 6.5965D
     $     -4, 6.5617D-4, 6.5279D-4, 6.4949D-4, 6.4627D-4, 6.4313D-4,
     $     6.4005D-4, 6.3704D-4, 6.3409D-4, 6.069D-4, 5.8238D-4, 5.6912D
     $     -4, 5.5419D-4, 5.3809D-4, 5.2132D-4, 5.0433D-4, 4.8756D-4,
     $     4.714D-4, 4.5622D-4, 4.4239D-4, 4.3023D-4/


C     4009, JENDL evaluated data
      nbins(5) = 387
      data ( eg(5,i), i=1, 387) /1573500.0, 1660000.0, 1670000.0,
     $     1680000.0, 1690000.0, 1700000.0, 1710000.0, 1720000.0,
     $     1730000.0, 1740000.0, 1750000.0, 1760000.0, 1770000.0,
     $     1780000.0, 1790000.0, 1800000.0, 1810000.0, 1820000.0,
     $     1830000.0, 1840000.0, 1850000.0, 1860000.0, 1870000.0,
     $     1880000.0, 1890000.0, 1900000.0, 1910000.0, 1920000.0,
     $     1930000.0, 1940000.0, 1950000.0, 1960000.0, 1970000.0,
     $     1980000.0, 1990000.0, 2000000.0, 2010000.0, 2020000.0,
     $     2030000.0, 2040000.0, 2050000.0, 2060000.0, 2070000.0,
     $     2080000.0, 2090000.0, 2100000.0, 2110000.0, 2120000.0,
     $     2130000.0, 2140000.0, 2150000.0, 2160000.0, 2170000.0,
     $     2180000.0, 2190000.0, 2200000.0, 2210000.0, 2220000.0,
     $     2230000.0, 2240000.0, 2250000.0, 2260000.0, 2270000.0,
     $     2280000.0, 2290000.0, 2300000.0, 2310000.0, 2320000.0,
     $     2330000.0, 2340000.0, 2350000.0, 2360000.0, 2370000.0,
     $     2380000.0, 2390000.0, 2400000.0, 2410000.0, 2420000.0,
     $     2430000.0, 2440000.0, 2450000.0, 2460000.0, 2470000.0,
     $     2480000.0, 2490000.0, 2500000.0, 2510000.0, 2520000.0,
     $     2530000.0, 2540000.0, 2550000.0, 2560000.0, 2570000.0,
     $     2580000.0, 2590000.0, 2600000.0, 2650000.0, 2700000.0,
     $     2750000.0, 2800000.0, 2850000.0, 2900000.0, 2950000.0,
     $     3000000.0, 3050000.0, 3100000.0, 3150000.0, 3200000.0,
     $     3250000.0, 3300000.0, 3350000.0, 3400000.0, 3450000.0,
     $     3500000.0, 4000000.0, 4500000.0, 5000000.0, 5500000.0,
     $     6000000.0, 6500000.0, 7000000.0, 7500000.0, 8000000.0,
     $     8500000.0, 9000000.0, 9500000.0, 1.0D7, 1.05D7, 1.1D7, 1.15D7
     $     , 1.2D7, 1.25D7, 1.3D7, 1.35D7, 1.4D7, 1.45D7, 1.5D7, 1.55D7,
     $     1.6D7, 1.65D7, 1.7D7, 1.75D7, 1.8D7, 1.85D7, 1.9D7, 1.95D7,
     $     2.0D7, 2.05D7, 2.1D7, 2.15D7, 2.2D7, 2.25D7, 2.3D7, 2.35D7,
     $     2.4D7, 2.45D7, 2.5D7, 2.55D7, 2.6D7, 2.65D7, 2.7D7, 2.75D7,
     $     2.8D7, 2.85D7, 2.9D7, 2.95D7, 3.0D7, 3.05D7, 3.1D7, 3.15D7,
     $     3.2D7, 3.25D7, 3.3D7, 3.35D7, 3.4D7, 3.45D7, 3.5D7, 3.55D7,
     $     3.6D7, 3.65D7, 3.7D7, 3.75D7, 3.8D7, 3.85D7, 3.9D7, 3.95D7,
     $     4.0D7, 4.05D7, 4.1D7, 4.15D7, 4.2D7, 4.25D7, 4.3D7, 4.35D7,
     $     4.4D7, 4.45D7, 4.5D7, 4.55D7, 4.6D7, 4.65D7, 4.7D7, 4.75D7,
     $     4.8D7, 4.85D7, 4.9D7, 4.95D7, 5.0D7, 5.05D7, 5.1D7, 5.15D7,
     $     5.2D7, 5.25D7, 5.3D7, 5.35D7, 5.4D7, 5.45D7, 5.5D7, 5.55D7,
     $     5.6D7, 5.65D7, 5.7D7, 5.75D7, 5.8D7, 5.85D7, 5.9D7, 5.95D7,
     $     6.0D7, 6.05D7, 6.1D7, 6.15D7, 6.2D7, 6.25D7, 6.3D7, 6.35D7,
     $     6.4D7, 6.45D7, 6.5D7, 6.55D7, 6.6D7, 6.65D7, 6.7D7, 6.75D7,
     $     6.8D7, 6.85D7, 6.9D7, 6.95D7, 7.0D7, 7.05D7, 7.1D7, 7.15D7,
     $     7.2D7, 7.25D7, 7.3D7, 7.35D7, 7.4D7, 7.45D7, 7.5D7, 7.55D7,
     $     7.6D7, 7.65D7, 7.7D7, 7.75D7, 7.8D7, 7.85D7, 7.9D7, 7.95D7,
     $     8.0D7, 8.05D7, 8.1D7, 8.15D7, 8.2D7, 8.25D7, 8.3D7, 8.35D7,
     $     8.4D7, 8.45D7, 8.5D7, 8.55D7, 8.6D7, 8.65D7, 8.7D7, 8.75D7,
     $     8.8D7, 8.85D7, 8.9D7, 8.95D7, 9.0D7, 9.05D7, 9.1D7, 9.15D7,
     $     9.2D7, 9.25D7, 9.3D7, 9.35D7, 9.4D7, 9.45D7, 9.5D7, 9.55D7,
     $     9.6D7, 9.65D7, 9.7D7, 9.75D7, 9.8D7, 9.85D7, 9.9D7, 9.95D7,
     $     1.0D8, 1.005D8, 1.01D8, 1.015D8, 1.02D8, 1.025D8, 1.03D8,
     $     1.035D8, 1.04D8, 1.045D8, 1.05D8, 1.055D8, 1.06D8, 1.065D8,
     $     1.07D8, 1.075D8, 1.08D8, 1.085D8, 1.09D8, 1.095D8, 1.1D8,
     $     1.105D8, 1.11D8, 1.115D8, 1.12D8, 1.125D8, 1.13D8, 1.135D8,
     $     1.14D8, 1.145D8, 1.15D8, 1.155D8, 1.16D8, 1.165D8, 1.17D8,
     $     1.175D8, 1.18D8, 1.185D8, 1.19D8, 1.195D8, 1.2D8, 1.205D8,
     $     1.21D8, 1.215D8, 1.22D8, 1.225D8, 1.23D8, 1.235D8, 1.24D8,
     $     1.245D8, 1.25D8, 1.255D8, 1.26D8, 1.265D8, 1.27D8, 1.275D8,
     $     1.28D8, 1.285D8, 1.29D8, 1.295D8, 1.3D8, 1.305D8, 1.31D8,
     $     1.315D8, 1.32D8, 1.325D8, 1.33D8, 1.335D8, 1.34D8, 1.345D8,
     $     1.35D8, 1.355D8, 1.36D8, 1.365D8, 1.37D8, 1.375D8, 1.38D8,
     $     1.385D8, 1.39D8, 1.395D8, 1.4D8/

      data ( gdr(5,i), i=1, 387) /0.0, 0.0, 0.0010032, 0.0013043,
     $     0.001257, 9.7665D-4, 7.1369D-4, 5.3049D-4, 4.0969D-4, 3.2887D
     $     -4, 2.7305D-4, 2.3318D-4, 2.0383D-4, 1.8165D-4, 1.6452D-4,
     $     1.5105D-4, 1.4028D-4, 1.3157D-4, 1.2446D-4, 1.186D-4, 1.1374D
     $     -4, 1.0971D-4, 1.0635D-4, 1.0355D-4, 1.0122D-4, 9.9305D-5,
     $     9.774D-5, 9.6482D-5, 9.5496D-5, 9.4754D-5, 9.4231D-5, 9.391D
     $     -5, 9.3773D-5, 9.381D-5, 9.401D-5, 9.4367D-5, 9.4876D-5,
     $     9.5533D-5, 9.6336D-5, 9.7287D-5, 9.8386D-5, 9.9637D-5,
     $     1.0105D-4, 1.0262D-4, 1.0437D-4, 1.063D-4, 1.0843D-4, 1.1077D
     $     -4, 1.1335D-4, 1.1619D-4, 1.1931D-4, 1.2275D-4, 1.2656D-4,
     $     1.3077D-4, 1.3544D-4, 1.4065D-4, 1.4648D-4, 1.5304D-4,
     $     1.6045D-4, 1.6886D-4, 1.7849D-4, 1.8957D-4, 2.0243D-4,
     $     2.1747D-4, 2.3522D-4, 2.5636D-4, 2.818D-4, 3.1274D-4, 3.5076D
     $     -4, 3.9799D-4, 4.572D-4, 5.3188D-4, 6.2615D-4, 7.4392D-4,
     $     8.866D-4, 0.0010483, 0.0012087, 0.00133, 0.001371, 0.0013167,
     $     0.0011912, 0.0010368, 8.8697D-4, 7.5756D-4, 6.5224D-4,
     $     5.6892D-4, 5.0377D-4, 4.5301D-4, 4.1345D-4, 3.8259D-4,
     $     3.5852D-4, 3.398D-4, 3.2537D-4, 3.144D-4, 3.0627D-4, 3.0053D
     $     -4, 2.977D-4, 3.2637D-4, 3.8429D-4, 4.7993D-4, 6.3163D-4,
     $     8.6571D-4, 0.0011915, 0.0015102, 0.001575, 0.0013307,
     $     0.0010106, 7.6058D-4, 5.9245D-4, 4.8278D-4, 4.1073D-4,
     $     3.6272D-4, 3.3045D-4, 3.0882D-4, 3.1635D-4, 5.3367D-4,
     $     0.0010842, 0.0011317, 0.001182, 0.001424, 0.0014591,
     $     0.0014434, 0.0014864, 0.0015183, 0.0014822, 0.001379,
     $     0.0012438, 0.00111, 9.9605D-4, 9.0806D-4, 8.4577D-4, 8.0688D
     $     -4, 7.8901D-4, 7.9103D-4, 8.1997D-4, 9.5956D-4, 9.0251D-4,
     $     9.7207D-4, 0.0010686, 0.0011942, 0.0013563, 0.0015677,
     $     0.0018595, 0.0024105, 0.0029585, 0.0031394, 0.0037143,
     $     0.0043916, 0.0050402, 0.0055, 0.0056594, 0.0055486, 0.0053071
     $     , 0.0050732, 0.0049238, 0.0048786, 0.0049217, 0.0050173,
     $     0.0051204, 0.0051873, 0.0051881, 0.0051158, 0.0049871,
     $     0.0048333, 0.0046891, 0.0045829, 0.0045306, 0.0045298,
     $     0.0045512, 0.0045435, 0.0044608, 0.0042987, 0.0040954,
     $     0.0038966, 0.003731, 0.003608, 0.0035258, 0.003478, 0.0034568
     $     , 0.0034548, 0.0034652, 0.0034817, 0.0034988, 0.0035117,
     $     0.0035163, 0.0035098, 0.0034899, 0.0034556, 0.003407,
     $     0.0033449, 0.0032708, 0.0031866, 0.0030946, 0.0029969,
     $     0.0028958, 0.0027932, 0.0026906, 0.0025895, 0.0024909,
     $     0.0023956, 0.0023041, 0.0022168, 0.0021339, 0.0020555,
     $     0.0019815, 0.0019118, 0.0018463, 0.0017848, 0.0017271,
     $     0.001673, 0.0016222, 0.0015747, 0.0015301, 0.0014883,
     $     0.001449, 0.0014121, 0.0013774, 0.0013448, 0.0013142,
     $     0.0012852, 0.001258, 0.0012322, 0.0012079, 0.0011849,
     $     0.0011631, 0.0011425, 0.001123, 0.0011044, 0.0010867,
     $     0.0010699, 0.0010539, 0.0010387, 0.0010241, 0.0010102,
     $     9.9695D-4, 9.8424D-4, 9.7206D-4, 9.6038D-4, 9.4918D-4,
     $     9.3842D-4, 9.2807D-4, 9.1811D-4, 9.0852D-4, 8.9928D-4,
     $     8.9036D-4, 8.8174D-4, 8.7341D-4, 8.6536D-4, 8.5756D-4, 8.5D-4
     $     , 8.4268D-4, 8.3556D-4, 8.2866D-4, 8.2194D-4, 8.1541D-4,
     $     8.0905D-4, 8.0286D-4, 7.9682D-4, 7.9093D-4, 7.8518D-4,
     $     7.7956D-4, 7.7407D-4, 7.6869D-4, 7.6343D-4, 7.5828D-4,
     $     7.5324D-4, 7.4829D-4, 7.4343D-4, 7.3867D-4, 7.3399D-4,
     $     7.2939D-4, 7.2486D-4, 7.2042D-4, 7.1604D-4, 7.1173D-4,
     $     7.0748D-4, 7.0329D-4, 6.9917D-4, 6.951D-4, 6.9108D-4, 6.8712D
     $     -4, 6.832D-4, 6.7934D-4, 6.7552D-4, 6.7174D-4, 6.6801D-4,
     $     6.6432D-4, 6.6314D-4, 6.6195D-4, 6.6074D-4, 6.5951D-4,
     $     6.5826D-4, 6.5699D-4, 6.557D-4, 6.5439D-4, 6.5307D-4, 6.5173D
     $     -4, 6.5037D-4, 6.49D-4, 6.4762D-4, 6.4622D-4, 6.448D-4,
     $     6.4338D-4, 6.4194D-4, 6.4048D-4, 6.3902D-4, 6.3754D-4,
     $     6.3606D-4, 6.3456D-4, 6.3305D-4, 6.3153D-4, 6.3001D-4,
     $     6.2847D-4, 6.2693D-4, 6.2538D-4, 6.2382D-4, 6.2225D-4,
     $     6.2068D-4, 6.1911D-4, 6.1752D-4, 6.1594D-4, 6.1434D-4,
     $     6.1275D-4, 6.1115D-4, 6.0954D-4, 6.0794D-4, 6.0633D-4,
     $     6.0472D-4, 6.0311D-4, 6.0149D-4, 5.9988D-4, 5.9826D-4,
     $     5.9665D-4, 5.9503D-4, 5.9342D-4, 5.9181D-4, 5.902D-4, 5.8859D
     $     -4, 5.8698D-4, 5.8538D-4, 5.8378D-4, 5.8218D-4, 5.8059D-4,
     $     5.79D-4, 5.7742D-4, 5.7584D-4, 5.7427D-4, 5.727D-4, 5.7114D-4
     $     , 5.6959D-4, 5.6804D-4, 5.665D-4, 5.6496D-4, 5.6344D-4,
     $     5.6192D-4, 5.6042D-4, 5.5892D-4, 5.5743D-4, 5.5595D-4,
     $     5.5448D-4, 5.5302D-4, 5.5157D-4, 5.5013D-4, 5.4871D-4,
     $     5.4729D-4, 5.4589D-4, 5.445D-4, 5.4313D-4, 5.4176D-4, 5.4041D
     $     -4, 5.3908D-4, 5.3775D-4, 5.3645D-4, 5.3515D-4, 5.3387D-4,
     $     5.3261D-4, 5.3136D-4, 5.3013D-4, 5.2892D-4, 5.2772D-4,
     $     5.2654D-4, 5.2537D-4, 5.2423D-4, 5.231D-4, 5.2199D-4, 5.2089D
     $     -4, 5.1982D-4/


C     5010, JENDL evaluated data
      nbins(6) = 241
      data ( eg(6,i), i=1, 241) /4461000.0, 4600000.0, 4800000.0,
     $     5000000.0, 5200000.0, 5400000.0, 5600000.0, 5800000.0,
     $     6000000.0, 6200000.0, 6400000.0, 6600000.0, 6800000.0,
     $     7000000.0, 7200000.0, 7400000.0, 7600000.0, 7800000.0,
     $     8000000.0, 8200000.0, 8400000.0, 8600000.0, 8800000.0,
     $     9000000.0, 9200000.0, 9400000.0, 9600000.0, 9800000.0, 1.0D7,
     $     1.02D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7
     $     , 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7,
     $     1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7
     $     , 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7,
     $     1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7
     $     , 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7,
     $     1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7
     $     , 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7,
     $     2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7
     $     , 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7,
     $     2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7
     $     , 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7,
     $     2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7
     $     , 2.98D7, 3.0D7, 3.05D7, 3.1D7, 3.15D7, 3.2D7, 3.25D7, 3.3D7,
     $     3.35D7, 3.4D7, 3.45D7, 3.5D7, 3.55D7, 3.6D7, 3.65D7, 3.7D7,
     $     3.75D7, 3.8D7, 3.85D7, 3.9D7, 3.95D7, 4.0D7, 4.05D7, 4.1D7,
     $     4.15D7, 4.2D7, 4.25D7, 4.3D7, 4.35D7, 4.4D7, 4.45D7, 4.5D7,
     $     4.55D7, 4.6D7, 4.65D7, 4.7D7, 4.75D7, 4.8D7, 4.85D7, 4.9D7,
     $     4.95D7, 5.0D7, 5.05D7, 5.1D7, 5.15D7, 5.2D7, 5.25D7, 5.3D7,
     $     5.35D7, 5.4D7, 5.45D7, 5.5D7, 5.55D7, 5.6D7, 5.65D7, 5.7D7,
     $     5.75D7, 5.8D7, 5.85D7, 5.9D7, 5.95D7, 6.0D7, 6.05D7, 6.1D7,
     $     6.15D7, 6.2D7, 6.25D7, 6.3D7, 6.35D7, 6.4D7, 6.45D7, 6.5D7,
     $     6.55D7, 6.6D7, 6.65D7, 6.7D7, 6.75D7, 6.8D7, 6.85D7, 6.9D7,
     $     6.95D7, 7.0D7, 7.05D7, 7.1D7, 7.15D7, 7.2D7, 7.25D7, 7.3D7,
     $     7.35D7, 7.4D7, 7.45D7, 7.5D7, 7.55D7, 7.6D7, 7.65D7, 7.7D7,
     $     7.75D7, 7.8D7, 7.85D7, 7.9D7, 7.95D7, 8.0D7, 8.5D7, 9.0D7,
     $     9.5D7, 1.0D8, 1.05D8, 1.1D8, 1.15D8, 1.2D8, 1.25D8, 1.3D8,
     $     1.35D8, 1.4D8/

      data ( gdr(6,i), i=1, 241) /0.0, 8.4372D-6, 2.9316D-5, 4.7866D-5,
     $     6.678D-5, 8.6213D-5, 1.0452D-4, 1.235D-4, 1.4452D-4, 1.6773D
     $     -4, 1.9333D-4, 2.1657D-4, 2.3983D-4, 2.6596D-4, 2.9842D-4,
     $     3.8375D-4, 4.3891D-4, 5.354D-4, 4.3296D-4, 4.6366D-4, 5.0551D
     $     -4, 5.5097D-4, 6.4057D-4, 6.7514D-4, 6.9684D-4, 7.5367D-4,
     $     8.2071D-4, 8.9789D-4, 9.8738D-4, 0.0010928, 0.0012199,
     $     0.0013725, 0.001567, 0.0018305, 0.0022063, 0.0027662, 0.00359
     $     , 0.0045778, 0.0050926, 0.0046853, 0.0039937, 0.003512,
     $     0.0032748, 0.003208, 0.0032514, 0.003369, 0.0035396,
     $     0.0037497, 0.0039895, 0.0042501, 0.0045224, 0.004793,
     $     0.0050529, 0.0052919, 0.0054987, 0.0056639, 0.0057817,
     $     0.0058511, 0.0058762, 0.0058654, 0.0058306, 0.0057855,
     $     0.0057441, 0.005719, 0.0057198, 0.0057601, 0.0058495,
     $     0.0059958, 0.0062043, 0.006476, 0.0068052, 0.0071752,
     $     0.0075559, 0.0079039, 0.008171, 0.0083201, 0.0083393,
     $     0.0082529, 0.008109, 0.0079618, 0.0078583, 0.0078322,
     $     0.0079036, 0.0080803, 0.008358, 0.0087188, 0.009127,
     $     0.0095273, 0.0098496, 0.010024, 0.010011, 0.0098108,
     $     0.0094656, 0.0090428, 0.0086101, 0.0082212, 0.007913,
     $     0.0077088, 0.0076243, 0.0076718, 0.0078619, 0.0081969,
     $     0.0086637, 0.0092045, 0.0096935, 0.0099575, 0.0098761,
     $     0.0094797, 0.0089185, 0.0083462, 0.0078519, 0.007465,
     $     0.0071818, 0.0069864, 0.0068602, 0.0067891, 0.0067602,
     $     0.006763, 0.0067889, 0.006831, 0.0068833, 0.0069408, 0.006999
     $     , 0.007054, 0.0071024, 0.0071413, 0.0071676, 0.0071797,
     $     0.007176, 0.0070923, 0.0069005, 0.006613, 0.0062535,
     $     0.0058495, 0.0054275, 0.005009, 0.004609, 0.0042364,
     $     0.0038952, 0.0035873, 0.0033119, 0.0030671, 0.0028501,
     $     0.0026574, 0.0024871, 0.0023363, 0.0022029, 0.0020844,
     $     0.0019788, 0.0018847, 0.0018008, 0.0017259, 0.0016587,
     $     0.0015981, 0.0015436, 0.0014944, 0.0014501, 0.0014099,
     $     0.001373, 0.0013395, 0.001309, 0.0012812, 0.0012559,
     $     0.0012327, 0.0012115, 0.001192, 0.0011741, 0.0011577,
     $     0.0011426, 0.0011286, 0.0011157, 0.0011038, 0.0010927,
     $     0.0010824, 0.0010729, 0.001064, 0.0010557, 0.001048,
     $     0.0010405, 0.0010336, 0.001027, 0.0010209, 0.0010151,
     $     0.0010097, 0.0010045, 9.9969D-4, 9.951D-4, 9.9074D-4, 9.8659D
     $     -4, 9.8264D-4, 9.7887D-4, 9.7526D-4, 9.7181D-4, 9.6848D-4,
     $     9.6529D-4, 9.622D-4, 9.5922D-4, 9.5633D-4, 9.5344D-4, 9.5063D
     $     -4, 9.4789D-4, 9.4523D-4, 9.4262D-4, 9.4006D-4, 9.3756D-4,
     $     9.3509D-4, 9.3267D-4, 9.3028D-4, 9.2791D-4, 9.2557D-4,
     $     9.2326D-4, 9.2096D-4, 9.1867D-4, 9.164D-4, 9.1414D-4, 9.1189D
     $     -4, 9.0964D-4, 9.0739D-4, 9.051D-4, 9.0282D-4, 9.0054D-4,
     $     8.9825D-4, 8.9597D-4, 8.9367D-4, 8.9138D-4, 8.8907D-4,
     $     8.8676D-4, 8.8444D-4, 8.821D-4, 8.5812D-4, 8.3264D-4, 8.1281D
     $     -4, 7.9045D-4, 7.6636D-4, 7.4126D-4, 7.1584D-4, 6.9073D-4,
     $     6.6654D-4, 6.4383D-4, 6.2312D-4, 6.0493D-4/


C     5011, JENDL evaluated data
      nbins(7) = 220
      data ( eg(7,i), i=1, 220) /8665000.0, 8800000.0, 9000000.0,
     $     9200000.0, 9400000.0, 9600000.0, 9800000.0, 1.0D7, 1.02D7,
     $     1.04D7, 1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7, 1.18D7
     $     , 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7, 1.32D7,
     $     1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7, 1.48D7
     $     , 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7, 1.62D7,
     $     1.64D7, 1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7, 1.78D7
     $     , 1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7, 1.92D7,
     $     1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7, 2.08D7
     $     , 2.1D7, 2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7, 2.22D7,
     $     2.24D7, 2.26D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7, 2.38D7
     $     , 2.4D7, 2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7, 2.52D7,
     $     2.54D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7, 2.68D7
     $     , 2.7D7, 2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7, 2.82D7,
     $     2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7, 2.98D7
     $     , 3.0D7, 3.05D7, 3.1D7, 3.15D7, 3.2D7, 3.25D7, 3.3D7, 3.35D7,
     $     3.4D7, 3.45D7, 3.5D7, 3.55D7, 3.6D7, 3.65D7, 3.7D7, 3.75D7,
     $     3.8D7, 3.85D7, 3.9D7, 3.95D7, 4.0D7, 4.05D7, 4.1D7, 4.15D7,
     $     4.2D7, 4.25D7, 4.3D7, 4.35D7, 4.4D7, 4.45D7, 4.5D7, 4.55D7,
     $     4.6D7, 4.65D7, 4.7D7, 4.75D7, 4.8D7, 4.85D7, 4.9D7, 4.95D7,
     $     5.0D7, 5.05D7, 5.1D7, 5.15D7, 5.2D7, 5.25D7, 5.3D7, 5.35D7,
     $     5.4D7, 5.45D7, 5.5D7, 5.55D7, 5.6D7, 5.65D7, 5.7D7, 5.75D7,
     $     5.8D7, 5.85D7, 5.9D7, 5.95D7, 6.0D7, 6.05D7, 6.1D7, 6.15D7,
     $     6.2D7, 6.25D7, 6.3D7, 6.35D7, 6.4D7, 6.45D7, 6.5D7, 6.55D7,
     $     6.6D7, 6.65D7, 6.7D7, 6.75D7, 6.8D7, 6.85D7, 6.9D7, 6.95D7,
     $     7.0D7, 7.05D7, 7.1D7, 7.15D7, 7.2D7, 7.25D7, 7.3D7, 7.35D7,
     $     7.4D7, 7.45D7, 7.5D7, 7.55D7, 7.6D7, 7.65D7, 7.7D7, 7.75D7,
     $     7.8D7, 7.85D7, 7.9D7, 7.95D7, 8.0D7, 8.5D7, 9.0D7, 9.5D7,
     $     1.0D8, 1.05D8, 1.1D8, 1.15D8, 1.2D8, 1.25D8, 1.3D8, 1.35D8,
     $     1.4D8/

      data ( gdr(7,i), i=1, 220) /0.0, 4.9764D-5, 1.7447D-4, 3.0614D-4,
     $     4.6463D-4, 6.5783D-4, 8.7954D-4, 0.0011321, 0.0014069,
     $     0.0016997, 0.0020275, 0.0023715, 0.0027426, 0.0031348,
     $     0.003512, 0.0038355, 0.0040803, 0.004247, 0.0043594,
     $     0.0044515, 0.0045536, 0.0046662, 0.0048046, 0.004975,
     $     0.0051654, 0.0053569, 0.0055272, 0.0056558, 0.0057295,
     $     0.0057458, 0.005713, 0.00564, 0.0055491, 0.0054629, 0.0053989
     $     , 0.0053721, 0.0053951, 0.005479, 0.0056345, 0.0058707,
     $     0.0061931, 0.0065872, 0.0070252, 0.0074427, 0.0077512,
     $     0.0078969, 0.0079031, 0.0078431, 0.0077823, 0.0077519,
     $     0.0077549, 0.0077759, 0.0078019, 0.0078198, 0.0078174,
     $     0.0077869, 0.0077249, 0.0076319, 0.0075121, 0.0073719,
     $     0.0072192, 0.0070621, 0.0069082, 0.0067634, 0.0066299,
     $     0.0065132, 0.0064137, 0.0063296, 0.0062569, 0.0061905,
     $     0.0061252, 0.0060573, 0.0059863, 0.0059148, 0.0058485,
     $     0.0057949, 0.0057597, 0.0057512, 0.005776, 0.0058382,
     $     0.0059395, 0.0060789, 0.0062523, 0.0064514, 0.0066634,
     $     0.0068701, 0.0070498, 0.0071792, 0.0072379, 0.0072124,
     $     0.0071019, 0.0069154, 0.0066706, 0.0063897, 0.0060948,
     $     0.0058052, 0.0055365, 0.0053, 0.0051045, 0.0049568, 0.0048631
     $     , 0.0048283, 0.004858, 0.004953, 0.0051012, 0.0052663,
     $     0.0053843, 0.0053876, 0.0049039, 0.0044063, 0.0042935,
     $     0.0045659, 0.0052258, 0.0062558, 0.0073563, 0.0077524,
     $     0.007035, 0.0058385, 0.0048032, 0.0040801, 0.0036083,
     $     0.0033059, 0.0031106, 0.0029795, 0.002883, 0.0028005,
     $     0.002719, 0.0026318, 0.0025364, 0.0024342, 0.0023281,
     $     0.002222, 0.0021194, 0.0020227, 0.0019333, 0.0018519,
     $     0.0017785, 0.0017129, 0.0016547, 0.0016031, 0.0015573,
     $     0.0015168, 0.0014809, 0.0014492, 0.001421, 0.0013958,
     $     0.0013733, 0.0013533, 0.0013354, 0.0013194, 0.001305,
     $     0.0012921, 0.0012805, 0.0012701, 0.0012606, 0.0012521,
     $     0.0012443, 0.0012372, 0.0012308, 0.0012249, 0.0012195,
     $     0.0012146, 0.00121, 0.0012058, 0.0012019, 0.0011982,
     $     0.0011947, 0.0011914, 0.0011883, 0.0011854, 0.0011826,
     $     0.00118, 0.0011775, 0.0011752, 0.0011729, 0.0011707,
     $     0.0011686, 0.0011665, 0.0011645, 0.0011625, 0.0011606,
     $     0.0011587, 0.0011568, 0.0011549, 0.0011531, 0.0011512,
     $     0.0011492, 0.0011473, 0.0011454, 0.0011435, 0.0011416,
     $     0.0011397, 0.0011377, 0.0011357, 0.0011337, 0.0011317,
     $     0.0011297, 0.0011276, 0.0011256, 0.0011235, 0.0011213,
     $     0.0011191, 0.001117, 0.0011147, 0.0011125, 0.0011102,
     $     0.0011078, 0.0011054, 0.00108, 0.0010515, 0.0010259, 9.9718D
     $     -4, 9.6618D-4, 9.3389D-4, 9.0119D-4, 8.6889D-4, 8.3777D-4,
     $     8.0855D-4, 7.8191D-4, 7.585D-4/


C     6012, JENDL evaluated data
      nbins(8) = 492
      data ( eg(8,i), i=1, 492) /7367000.0, 7400000.0, 7450000.0,
     $     7500000.0, 7550000.0, 7600000.0, 7650000.0, 7700000.0,
     $     7750000.0, 7800000.0, 7850000.0, 7900000.0, 7950000.0,
     $     8000000.0, 8050000.0, 8100000.0, 8150000.0, 8200000.0,
     $     8250000.0, 8300000.0, 8350000.0, 8400000.0, 8450000.0,
     $     8500000.0, 8550000.0, 8600000.0, 8650000.0, 8700000.0,
     $     8750000.0, 8800000.0, 8850000.0, 8900000.0, 8950000.0,
     $     9000000.0, 9050000.0, 9100000.0, 9150000.0, 9200000.0,
     $     9250000.0, 9300000.0, 9350000.0, 9400000.0, 9450000.0,
     $     9500000.0, 9550000.0, 9600000.0, 9650000.0, 9700000.0,
     $     9750000.0, 9800000.0, 9850000.0, 9900000.0, 9950000.0, 1.0D7,
     $     1.005D7, 1.01D7, 1.015D7, 1.02D7, 1.025D7, 1.03D7, 1.035D7,
     $     1.04D7, 1.045D7, 1.05D7, 1.055D7, 1.06D7, 1.065D7, 1.07D7,
     $     1.075D7, 1.08D7, 1.085D7, 1.09D7, 1.095D7, 1.1D7, 1.105D7,
     $     1.11D7, 1.115D7, 1.12D7, 1.125D7, 1.13D7, 1.135D7, 1.14D7,
     $     1.145D7, 1.15D7, 1.155D7, 1.16D7, 1.165D7, 1.17D7, 1.175D7,
     $     1.18D7, 1.185D7, 1.19D7, 1.195D7, 1.2D7, 1.205D7, 1.21D7,
     $     1.215D7, 1.22D7, 1.225D7, 1.23D7, 1.235D7, 1.24D7, 1.245D7,
     $     1.25D7, 1.255D7, 1.26D7, 1.265D7, 1.27D7, 1.275D7, 1.28D7,
     $     1.285D7, 1.29D7, 1.295D7, 1.3D7, 1.305D7, 1.31D7, 1.315D7,
     $     1.32D7, 1.325D7, 1.33D7, 1.335D7, 1.34D7, 1.345D7, 1.35D7,
     $     1.355D7, 1.36D7, 1.365D7, 1.37D7, 1.375D7, 1.38D7, 1.385D7,
     $     1.39D7, 1.395D7, 1.4D7, 1.405D7, 1.41D7, 1.415D7, 1.42D7,
     $     1.425D7, 1.43D7, 1.435D7, 1.44D7, 1.445D7, 1.45D7, 1.455D7,
     $     1.46D7, 1.465D7, 1.47D7, 1.475D7, 1.48D7, 1.485D7, 1.49D7,
     $     1.495D7, 1.5D7, 1.505D7, 1.51D7, 1.515D7, 1.52D7, 1.525D7,
     $     1.53D7, 1.535D7, 1.54D7, 1.545D7, 1.55D7, 1.555D7, 1.56D7,
     $     1.565D7, 1.57D7, 1.575D7, 1.58D7, 1.585D7, 1.59D7, 1.595D7,
     $     1.6D7, 1.605D7, 1.61D7, 1.615D7, 1.62D7, 1.625D7, 1.63D7,
     $     1.635D7, 1.64D7, 1.645D7, 1.65D7, 1.655D7, 1.66D7, 1.665D7,
     $     1.67D7, 1.675D7, 1.68D7, 1.685D7, 1.69D7, 1.695D7, 1.7D7,
     $     1.705D7, 1.71D7, 1.715D7, 1.72D7, 1.725D7, 1.73D7, 1.735D7,
     $     1.74D7, 1.745D7, 1.75D7, 1.755D7, 1.76D7, 1.765D7, 1.77D7,
     $     1.775D7, 1.78D7, 1.785D7, 1.79D7, 1.795D7, 1.8D7, 1.805D7,
     $     1.81D7, 1.815D7, 1.82D7, 1.825D7, 1.83D7, 1.835D7, 1.84D7,
     $     1.845D7, 1.85D7, 1.855D7, 1.86D7, 1.865D7, 1.87D7, 1.875D7,
     $     1.88D7, 1.885D7, 1.89D7, 1.895D7, 1.9D7, 1.905D7, 1.91D7,
     $     1.915D7, 1.92D7, 1.925D7, 1.93D7, 1.935D7, 1.94D7, 1.945D7,
     $     1.95D7, 1.955D7, 1.96D7, 1.965D7, 1.97D7, 1.975D7, 1.98D7,
     $     1.985D7, 1.99D7, 1.995D7, 2.0D7, 2.005D7, 2.01D7, 2.015D7,
     $     2.02D7, 2.025D7, 2.03D7, 2.035D7, 2.04D7, 2.045D7, 2.05D7,
     $     2.055D7, 2.06D7, 2.065D7, 2.07D7, 2.075D7, 2.08D7, 2.085D7,
     $     2.09D7, 2.095D7, 2.1D7, 2.105D7, 2.11D7, 2.115D7, 2.12D7,
     $     2.125D7, 2.13D7, 2.135D7, 2.14D7, 2.145D7, 2.15D7, 2.155D7,
     $     2.16D7, 2.165D7, 2.17D7, 2.175D7, 2.18D7, 2.185D7, 2.19D7,
     $     2.195D7, 2.2D7, 2.205D7, 2.21D7, 2.215D7, 2.22D7, 2.225D7,
     $     2.23D7, 2.235D7, 2.24D7, 2.245D7, 2.25D7, 2.255D7, 2.26D7,
     $     2.265D7, 2.27D7, 2.275D7, 2.28D7, 2.285D7, 2.29D7, 2.295D7,
     $     2.3D7, 2.305D7, 2.31D7, 2.315D7, 2.32D7, 2.325D7, 2.33D7,
     $     2.335D7, 2.34D7, 2.345D7, 2.35D7, 2.355D7, 2.36D7, 2.365D7,
     $     2.37D7, 2.375D7, 2.38D7, 2.385D7, 2.39D7, 2.395D7, 2.4D7,
     $     2.405D7, 2.41D7, 2.415D7, 2.42D7, 2.425D7, 2.43D7, 2.435D7,
     $     2.44D7, 2.445D7, 2.45D7, 2.455D7, 2.46D7, 2.465D7, 2.47D7,
     $     2.475D7, 2.48D7, 2.485D7, 2.49D7, 2.495D7, 2.5D7, 2.505D7,
     $     2.51D7, 2.515D7, 2.52D7, 2.525D7, 2.53D7, 2.535D7, 2.54D7,
     $     2.545D7, 2.55D7, 2.555D7, 2.56D7, 2.565D7, 2.57D7, 2.575D7,
     $     2.58D7, 2.585D7, 2.59D7, 2.595D7, 2.6D7, 2.605D7, 2.61D7,
     $     2.615D7, 2.62D7, 2.625D7, 2.63D7, 2.635D7, 2.64D7, 2.645D7,
     $     2.65D7, 2.655D7, 2.66D7, 2.665D7, 2.67D7, 2.675D7, 2.68D7,
     $     2.685D7, 2.69D7, 2.695D7, 2.7D7, 2.705D7, 2.71D7, 2.715D7,
     $     2.72D7, 2.725D7, 2.73D7, 2.735D7, 2.74D7, 2.745D7, 2.75D7,
     $     2.755D7, 2.76D7, 2.765D7, 2.77D7, 2.775D7, 2.78D7, 2.785D7,
     $     2.79D7, 2.795D7, 2.8D7, 2.805D7, 2.81D7, 2.815D7, 2.82D7,
     $     2.825D7, 2.83D7, 2.835D7, 2.84D7, 2.845D7, 2.85D7, 2.855D7,
     $     2.86D7, 2.865D7, 2.87D7, 2.875D7, 2.88D7, 2.885D7, 2.89D7,
     $     2.895D7, 2.9D7, 2.905D7, 2.91D7, 2.915D7, 2.92D7, 2.925D7,
     $     2.93D7, 2.935D7, 2.94D7, 2.945D7, 2.95D7, 2.955D7, 2.96D7,
     $     2.965D7, 2.97D7, 2.975D7, 2.98D7, 2.985D7, 2.99D7, 2.995D7,
     $     3.0D7, 3.1D7, 3.2D7, 3.3D7, 3.4D7, 3.5D7, 3.6D7, 3.7D7, 3.8D7
     $     , 3.9D7, 4.0D7, 4.1D7, 4.2D7, 4.3D7, 4.4D7, 4.5D7, 4.6D7,
     $     4.7D7, 4.8D7, 4.9D7, 5.0D7, 5.5D7, 6.0D7, 6.5D7, 7.0D7, 7.5D7
     $     , 8.0D7, 8.5D7, 9.0D7, 9.5D7, 1.0D8, 1.05D8, 1.1D8, 1.15D8,
     $     1.2D8, 1.25D8, 1.3D8, 1.35D8, 1.4D8/

      data ( gdr(8,i), i=1, 492) /6.655D-7, 6.6884D-7, 6.719D-7, 6.7467D
     $     -7, 6.7739D-7, 6.8016D-7, 6.8301D-7, 6.8593D-7, 6.8893D-7,
     $     6.9279D-7, 6.9726D-7, 7.0186D-7, 7.0659D-7, 7.1686D-7,
     $     7.3094D-7, 7.4523D-7, 7.5974D-7, 7.8692D-7, 8.2275D-7,
     $     8.5898D-7, 8.9562D-7, 9.6775D-7, 1.0641D-6, 1.1613D-6,
     $     1.2594D-6, 1.3584D-6, 1.4585D-6, 1.5595D-6, 1.6615D-6,
     $     1.7646D-6, 1.8688D-6, 1.974D-6, 2.0804D-6, 2.1879D-6, 2.2967D
     $     -6, 2.4066D-6, 2.5179D-6, 2.6304D-6, 2.7442D-6, 2.8594D-6,
     $     2.9761D-6, 3.0914D-6, 3.2063D-6, 3.3226D-6, 3.4404D-6,
     $     3.5595D-6, 3.6802D-6, 3.8025D-6, 3.9263D-6, 4.0518D-6, 4.179D
     $     -6, 4.3081D-6, 4.4389D-6, 4.5717D-6, 4.7064D-6, 4.8432D-6,
     $     4.9822D-6, 5.1234D-6, 5.267D-6, 5.413D-6, 5.5616D-6, 5.7129D
     $     -6, 5.867D-6, 6.0241D-6, 6.1843D-6, 6.3479D-6, 6.5151D-6,
     $     6.686D-6, 6.861D-6, 7.0402D-6, 7.2241D-6, 7.413D-6, 7.6073D-6
     $     , 7.8075D-6, 8.0141D-6, 8.2277D-6, 8.4491D-6, 8.679D-6,
     $     8.9185D-6, 9.1687D-6, 9.431D-6, 9.6825D-6, 9.9322D-6, 1.0199D
     $     -5, 1.0484D-5, 1.0792D-5, 1.1127D-5, 1.1493D-5, 1.1898D-5,
     $     1.235D-5, 1.286D-5, 1.3443D-5, 1.4118D-5, 1.4911D-5, 1.5858D
     $     -5, 1.7011D-5, 1.8439D-5, 2.0245D-5, 2.2572D-5, 2.5623D-5,
     $     2.9663D-5, 3.4988D-5, 4.1743D-5, 4.9448D-5, 5.6281D-5,
     $     5.9282D-5, 5.6849D-5, 5.0643D-5, 4.358D-5, 3.7431D-5, 3.2662D
     $     -5, 2.9139D-5, 2.6585D-5, 2.4745D-5, 2.3426D-5, 2.2488D-5,
     $     2.1834D-5, 2.1395D-5, 2.1123D-5, 2.0983D-5, 2.0952D-5,
     $     2.0992D-5, 2.1099D-5, 2.1274D-5, 2.1514D-5, 2.1814D-5,
     $     2.2173D-5, 2.2592D-5, 2.3074D-5, 2.3624D-5, 2.4251D-5,
     $     2.4964D-5, 2.578D-5, 2.6719D-5, 2.781D-5, 2.9093D-5, 3.0624D
     $     -5, 3.2487D-5, 3.4801D-5, 3.7751D-5, 4.1621D-5, 4.6859D-5,
     $     5.4185D-5, 6.4705D-5, 7.9877D-5, 1.0047D-4, 1.2267D-4,
     $     1.3383D-4, 1.2473D-4, 1.0495D-4, 8.696D-5, 7.4493D-5, 6.6824D
     $     -5, 6.2606D-5, 6.0881D-5, 6.1124D-5, 6.3171D-5, 6.7172D-5,
     $     7.3624D-5, 8.3445D-5, 9.7922D-5, 1.1772D-4, 1.3945D-4,
     $     1.5191D-4, 1.4701D-4, 1.3373D-4, 1.2382D-4, 1.2144D-4,
     $     1.2729D-4, 1.4226D-4, 1.6822D-4, 2.0604D-4, 2.4775D-4,
     $     2.8513D-4, 2.7442D-4, 2.4113D-4, 2.1069D-4, 1.9378D-4,
     $     1.8657D-4, 1.859D-4, 1.8974D-4, 2.0174D-4, 2.1739D-4, 2.3555D
     $     -4, 2.5609D-4, 2.8177D-4, 3.1066D-4, 3.4257D-4, 3.7782D-4,
     $     4.174D-4, 4.6142D-4, 5.1039D-4, 5.6499D-4, 6.2067D-4, 6.8186D
     $     -4, 7.4994D-4, 8.2531D-4, 9.0805D-4, 9.9771D-4, 0.001093,
     $     0.0011916, 0.0012896, 0.0013823, 0.0014639, 0.0015288,
     $     0.0015726, 0.0015929, 0.0015904, 0.0015678, 0.00153,
     $     0.0014821, 0.0014289, 0.0013745, 0.0013134, 0.0012554,
     $     0.0012028, 0.0011561, 0.0011154, 0.0010803, 0.0010506,
     $     0.0010259, 0.0010058, 9.8983D-4, 9.7777D-4, 9.6933D-4,
     $     9.6429D-4, 9.6252D-4, 9.6396D-4, 0.0010085, 0.0010638,
     $     0.00112, 0.0011805, 0.001236, 0.0012903, 0.0013547, 0.0014336
     $     , 0.0015321, 0.0016609, 0.0018363, 0.0020665, 0.0023154,
     $     0.0024628, 0.0024175, 0.0022677, 0.0021321, 0.0020481,
     $     0.0020104, 0.0020058, 0.0020205, 0.0020488, 0.0020898,
     $     0.0021411, 0.0022015, 0.0022669, 0.0023399, 0.0024213,
     $     0.0025116, 0.0026115, 0.0027221, 0.0028446, 0.0029805,
     $     0.0031318, 0.0033007, 0.0034899, 0.0037023, 0.0039413,
     $     0.0042103, 0.0045127, 0.0048492, 0.0052206, 0.0056248,
     $     0.0060532, 0.0064908, 0.006917, 0.00731, 0.0076547, 0.0079502
     $     , 0.0082138, 0.0084774, 0.0087799, 0.0091586, 0.0096402,
     $     0.01023, 0.010898, 0.011569, 0.012152, 0.0126, 0.012959,
     $     0.013314, 0.013732, 0.014234, 0.014797, 0.015379, 0.015923,
     $     0.016383, 0.016728, 0.016946, 0.017055, 0.017091, 0.017108,
     $     0.01716, 0.017296, 0.017552, 0.017946, 0.018473, 0.019098,
     $     0.019755, 0.02035, 0.020785, 0.020982, 0.020921, 0.020652,
     $     0.020282, 0.019941, 0.019754, 0.019806, 0.020112, 0.020576,
     $     0.020971, 0.021021, 0.020599, 0.019848, 0.019055, 0.018456,
     $     0.018157, 0.018141, 0.0183, 0.018471, 0.018495, 0.018299,
     $     0.017927, 0.017498, 0.017126, 0.016861, 0.01669, 0.016543,
     $     0.016325, 0.015963, 0.015436, 0.014781, 0.014068, 0.013368,
     $     0.012731, 0.012184, 0.011735, 0.011383, 0.01112, 0.010936,
     $     0.010818, 0.010758, 0.010752, 0.010791, 0.010871, 0.010987,
     $     0.011135, 0.011311, 0.011508, 0.011718, 0.011923, 0.012105,
     $     0.012253, 0.01237, 0.012465, 0.012546, 0.012614, 0.012669,
     $     0.012708, 0.012726, 0.012722, 0.012693, 0.012641, 0.012565,
     $     0.012465, 0.012343, 0.012202, 0.012045, 0.011875, 0.011696,
     $     0.011509, 0.01132, 0.011129, 0.010938, 0.010749, 0.010565,
     $     0.010386, 0.010214, 0.010049, 0.0098923, 0.0097437, 0.0096033
     $     , 0.0094717, 0.0093488, 0.0092342, 0.0091277, 0.0090289,
     $     0.0089375, 0.0088534, 0.0087753, 0.0087021, 0.0086363,
     $     0.0085785, 0.0085281, 0.0084546, 0.0083967, 0.0083454,
     $     0.0082988, 0.008256, 0.0082166, 0.0081804, 0.0081471,
     $     0.008116, 0.0080875, 0.0080616, 0.0080384, 0.0080179,
     $     0.0080001, 0.007985, 0.0079726, 0.0079614, 0.0079517,
     $     0.0079435, 0.0079378, 0.0079346, 0.0079326, 0.0079327,
     $     0.007935, 0.0079379, 0.0079421, 0.0079478, 0.0079547,
     $     0.0079628, 0.0079721, 0.0079822, 0.0079929, 0.0080031,
     $     0.0080132, 0.0080232, 0.0080328, 0.0080416, 0.0080495,
     $     0.0080565, 0.008062, 0.008067, 0.0080706, 0.0080724,
     $     0.0080723, 0.00807, 0.0080654, 0.0080585, 0.008049, 0.008037,
     $     0.0080223, 0.0080049, 0.0079848, 0.0071106, 0.0061401,
     $     0.0058063, 0.0060132, 0.0060795, 0.0055147, 0.0048223,
     $     0.0044294, 0.0043181, 0.0043527, 0.0043923, 0.0043207,
     $     0.0041104, 0.0038146, 0.0035173, 0.0032802, 0.0031248,
     $     0.0030395, 0.0029757, 0.0028678, 0.0019567, 0.0016333,
     $     0.001358, 0.0013698, 0.001369, 0.0013576, 0.0013374,
     $     0.0013102, 0.0012773, 0.0012402, 0.0012003, 0.0011587,
     $     0.0011166, 0.001075, 0.0010349, 9.9727D-4, 9.6296D-4, 9.3281D
     $     -4/


C     7014, JENDL evaluated data
      nbins(9) = 488
      data ( eg(9,i), i=1, 488) /7551000.0, 7600000.0, 7650000.0,
     $     7700000.0, 7750000.0, 7800000.0, 7850000.0, 7900000.0,
     $     7950000.0, 8000000.0, 8050000.0, 8100000.0, 8150000.0,
     $     8200000.0, 8250000.0, 8300000.0, 8350000.0, 8400000.0,
     $     8450000.0, 8500000.0, 8550000.0, 8600000.0, 8650000.0,
     $     8700000.0, 8750000.0, 8800000.0, 8850000.0, 8900000.0,
     $     8950000.0, 9000000.0, 9050000.0, 9100000.0, 9150000.0,
     $     9200000.0, 9250000.0, 9300000.0, 9350000.0, 9400000.0,
     $     9450000.0, 9500000.0, 9550000.0, 9600000.0, 9650000.0,
     $     9700000.0, 9750000.0, 9800000.0, 9850000.0, 9900000.0,
     $     9950000.0, 1.0D7, 1.005D7, 1.01D7, 1.015D7, 1.02D7, 1.025D7,
     $     1.03D7, 1.035D7, 1.04D7, 1.045D7, 1.05D7, 1.055D7, 1.06D7,
     $     1.065D7, 1.07D7, 1.075D7, 1.08D7, 1.085D7, 1.09D7, 1.095D7,
     $     1.1D7, 1.105D7, 1.11D7, 1.115D7, 1.12D7, 1.125D7, 1.13D7,
     $     1.135D7, 1.14D7, 1.145D7, 1.15D7, 1.155D7, 1.16D7, 1.165D7,
     $     1.17D7, 1.175D7, 1.18D7, 1.185D7, 1.19D7, 1.195D7, 1.2D7,
     $     1.205D7, 1.21D7, 1.215D7, 1.22D7, 1.225D7, 1.23D7, 1.235D7,
     $     1.24D7, 1.245D7, 1.25D7, 1.255D7, 1.26D7, 1.265D7, 1.27D7,
     $     1.275D7, 1.28D7, 1.285D7, 1.29D7, 1.295D7, 1.3D7, 1.305D7,
     $     1.31D7, 1.315D7, 1.32D7, 1.325D7, 1.33D7, 1.335D7, 1.34D7,
     $     1.345D7, 1.35D7, 1.355D7, 1.36D7, 1.365D7, 1.37D7, 1.375D7,
     $     1.38D7, 1.385D7, 1.39D7, 1.395D7, 1.4D7, 1.405D7, 1.41D7,
     $     1.415D7, 1.42D7, 1.425D7, 1.43D7, 1.435D7, 1.44D7, 1.445D7,
     $     1.45D7, 1.455D7, 1.46D7, 1.465D7, 1.47D7, 1.475D7, 1.48D7,
     $     1.485D7, 1.49D7, 1.495D7, 1.5D7, 1.505D7, 1.51D7, 1.515D7,
     $     1.52D7, 1.525D7, 1.53D7, 1.535D7, 1.54D7, 1.545D7, 1.55D7,
     $     1.555D7, 1.56D7, 1.565D7, 1.57D7, 1.575D7, 1.58D7, 1.585D7,
     $     1.59D7, 1.595D7, 1.6D7, 1.605D7, 1.61D7, 1.615D7, 1.62D7,
     $     1.625D7, 1.63D7, 1.635D7, 1.64D7, 1.645D7, 1.65D7, 1.655D7,
     $     1.66D7, 1.665D7, 1.67D7, 1.675D7, 1.68D7, 1.685D7, 1.69D7,
     $     1.695D7, 1.7D7, 1.705D7, 1.71D7, 1.715D7, 1.72D7, 1.725D7,
     $     1.73D7, 1.735D7, 1.74D7, 1.745D7, 1.75D7, 1.755D7, 1.76D7,
     $     1.765D7, 1.77D7, 1.775D7, 1.78D7, 1.785D7, 1.79D7, 1.795D7,
     $     1.8D7, 1.805D7, 1.81D7, 1.815D7, 1.82D7, 1.825D7, 1.83D7,
     $     1.835D7, 1.84D7, 1.845D7, 1.85D7, 1.855D7, 1.86D7, 1.865D7,
     $     1.87D7, 1.875D7, 1.88D7, 1.885D7, 1.89D7, 1.895D7, 1.9D7,
     $     1.905D7, 1.91D7, 1.915D7, 1.92D7, 1.925D7, 1.93D7, 1.935D7,
     $     1.94D7, 1.945D7, 1.95D7, 1.955D7, 1.96D7, 1.965D7, 1.97D7,
     $     1.975D7, 1.98D7, 1.985D7, 1.99D7, 1.995D7, 2.0D7, 2.005D7,
     $     2.01D7, 2.015D7, 2.02D7, 2.025D7, 2.03D7, 2.035D7, 2.04D7,
     $     2.045D7, 2.05D7, 2.055D7, 2.06D7, 2.065D7, 2.07D7, 2.075D7,
     $     2.08D7, 2.085D7, 2.09D7, 2.095D7, 2.1D7, 2.105D7, 2.11D7,
     $     2.115D7, 2.12D7, 2.125D7, 2.13D7, 2.135D7, 2.14D7, 2.145D7,
     $     2.15D7, 2.155D7, 2.16D7, 2.165D7, 2.17D7, 2.175D7, 2.18D7,
     $     2.185D7, 2.19D7, 2.195D7, 2.2D7, 2.205D7, 2.21D7, 2.215D7,
     $     2.22D7, 2.225D7, 2.23D7, 2.235D7, 2.24D7, 2.245D7, 2.25D7,
     $     2.255D7, 2.26D7, 2.265D7, 2.27D7, 2.275D7, 2.28D7, 2.285D7,
     $     2.29D7, 2.295D7, 2.3D7, 2.305D7, 2.31D7, 2.315D7, 2.32D7,
     $     2.325D7, 2.33D7, 2.335D7, 2.34D7, 2.345D7, 2.35D7, 2.355D7,
     $     2.36D7, 2.365D7, 2.37D7, 2.375D7, 2.38D7, 2.385D7, 2.39D7,
     $     2.395D7, 2.4D7, 2.405D7, 2.41D7, 2.415D7, 2.42D7, 2.425D7,
     $     2.43D7, 2.435D7, 2.44D7, 2.445D7, 2.45D7, 2.455D7, 2.46D7,
     $     2.465D7, 2.47D7, 2.475D7, 2.48D7, 2.485D7, 2.49D7, 2.495D7,
     $     2.5D7, 2.505D7, 2.51D7, 2.515D7, 2.52D7, 2.525D7, 2.53D7,
     $     2.535D7, 2.54D7, 2.545D7, 2.55D7, 2.555D7, 2.56D7, 2.565D7,
     $     2.57D7, 2.575D7, 2.58D7, 2.585D7, 2.59D7, 2.595D7, 2.6D7,
     $     2.605D7, 2.61D7, 2.615D7, 2.62D7, 2.625D7, 2.63D7, 2.635D7,
     $     2.64D7, 2.645D7, 2.65D7, 2.655D7, 2.66D7, 2.665D7, 2.67D7,
     $     2.675D7, 2.68D7, 2.685D7, 2.69D7, 2.695D7, 2.7D7, 2.705D7,
     $     2.71D7, 2.715D7, 2.72D7, 2.725D7, 2.73D7, 2.735D7, 2.74D7,
     $     2.745D7, 2.75D7, 2.755D7, 2.76D7, 2.765D7, 2.77D7, 2.775D7,
     $     2.78D7, 2.785D7, 2.79D7, 2.795D7, 2.8D7, 2.805D7, 2.81D7,
     $     2.815D7, 2.82D7, 2.825D7, 2.83D7, 2.835D7, 2.84D7, 2.845D7,
     $     2.85D7, 2.855D7, 2.86D7, 2.865D7, 2.87D7, 2.875D7, 2.88D7,
     $     2.885D7, 2.89D7, 2.895D7, 2.9D7, 2.905D7, 2.91D7, 2.915D7,
     $     2.92D7, 2.925D7, 2.93D7, 2.935D7, 2.94D7, 2.945D7, 2.95D7,
     $     2.955D7, 2.96D7, 2.965D7, 2.97D7, 2.975D7, 2.98D7, 2.985D7,
     $     2.99D7, 2.995D7, 3.0D7, 3.1D7, 3.2D7, 3.3D7, 3.4D7, 3.5D7,
     $     3.6D7, 3.7D7, 3.8D7, 3.9D7, 4.0D7, 4.1D7, 4.2D7, 4.3D7, 4.4D7
     $     , 4.5D7, 4.6D7, 4.7D7, 4.8D7, 4.9D7, 5.0D7, 5.5D7, 6.0D7,
     $     6.5D7, 7.0D7, 7.5D7, 8.0D7, 8.5D7, 9.0D7, 9.5D7, 1.0D8,
     $     1.05D8, 1.1D8, 1.15D8, 1.2D8, 1.25D8, 1.3D8, 1.35D8, 1.4D8/

      data ( gdr(9,i), i=1, 488) /1.846D-6, 2.1417D-6, 2.3573D-6,
     $     2.5948D-6, 2.8864D-6, 8.5224D-6, 1.833D-5, 3.8521D-5, 9.3247D
     $     -5, 4.9567D-4, 0.0081448, 0.0019207, 5.5736D-4, 3.4416D-4,
     $     2.7593D-4, 2.5693D-4, 2.6421D-4, 2.9663D-4, 3.5026D-4,
     $     4.3032D-4, 5.4515D-4, 6.9322D-4, 8.8373D-4, 0.0010953,
     $     0.001268, 0.0013286, 0.00126, 0.0011158, 9.5798D-4, 8.1966D-4
     $     , 7.0929D-4, 6.2458D-4, 5.6049D-4, 5.1215D-4, 4.7573D-4,
     $     4.4833D-4, 4.2785D-4, 4.128D-4, 4.0243D-4, 4.0136D-4, 3.9074D
     $     -4, 3.8328D-4, 3.7898D-4, 3.7661D-4, 3.758D-4, 3.7633D-4,
     $     3.7803D-4, 3.808D-4, 3.8454D-4, 3.8921D-4, 3.9477D-4, 4.0122D
     $     -4, 4.0857D-4, 4.1686D-4, 4.2615D-4, 4.3657D-4, 4.4823D-4,
     $     4.6131D-4, 4.7607D-4, 4.9287D-4, 5.1215D-4, 5.8488D-4,
     $     6.7158D-4, 7.7157D-4, 8.8933D-4, 9.9309D-4, 0.0011143,
     $     0.0012599, 0.0014279, 0.0015961, 0.0017425, 0.0018326,
     $     0.0018471, 0.0017935, 0.0017088, 0.0016224, 0.0015497,
     $     0.0014925, 0.0014537, 0.0014311, 0.0014219, 0.0014159,
     $     0.0014186, 0.0014289, 0.0014455, 0.0014677, 0.0014948,
     $     0.0015263, 0.0015619, 0.0016014, 0.0016446, 0.0016916,
     $     0.0017421, 0.0017963, 0.0018542, 0.0019157, 0.0019808,
     $     0.0020497, 0.0021222, 0.0021975, 0.0022744, 0.0023529,
     $     0.0024344, 0.0025185, 0.0026047, 0.0026923, 0.0027804,
     $     0.0028676, 0.0029434, 0.003016, 0.0030842, 0.0031468,
     $     0.0032024, 0.0032502, 0.003289, 0.0033157, 0.003294,
     $     0.0032624, 0.0032216, 0.0031726, 0.0031166, 0.0030536,
     $     0.0029866, 0.0029144, 0.0028025, 0.0026936, 0.0025889,
     $     0.0024894, 0.0023959, 0.0023089, 0.0022287, 0.0021542,
     $     0.0020685, 0.0019915, 0.001923, 0.0018625, 0.0018096,
     $     0.0017639, 0.001725, 0.0016922, 0.0016623, 0.0016353,
     $     0.0016135, 0.0015968, 0.001585, 0.0015777, 0.0015746,
     $     0.0015756, 0.0015804, 0.0015889, 0.001601, 0.0016167,
     $     0.0016358, 0.0016585, 0.0016848, 0.0017147, 0.0017485,
     $     0.0017862, 0.0018282, 0.0018747, 0.0019259, 0.0019814,
     $     0.0020421, 0.0021081, 0.0021794, 0.0022556, 0.0023356,
     $     0.0024181, 0.0025006, 0.0025799, 0.002652, 0.0027129,
     $     0.0027589, 0.0027879, 0.0027995, 0.0027957, 0.00278,
     $     0.0027567, 0.0027302, 0.0027041, 0.0026799, 0.0026615,
     $     0.0026495, 0.0026442, 0.0026457, 0.0026538, 0.0026681,
     $     0.0026884, 0.0027143, 0.0027452, 0.0027807, 0.0028204,
     $     0.0028638, 0.0029106, 0.0029605, 0.0030138, 0.0030707,
     $     0.0031317, 0.0031971, 0.0032674, 0.003343, 0.0034235,
     $     0.0035101, 0.0036034, 0.0037039, 0.0038123, 0.0039295,
     $     0.0040562, 0.0041936, 0.0043428, 0.0045052, 0.0046821,
     $     0.004875, 0.0050855, 0.005315, 0.0055645, 0.0058343,
     $     0.0061236, 0.00643, 0.0067485, 0.007071, 0.0073852, 0.0076778
     $     , 0.0079334, 0.0081379, 0.0082822, 0.0083639, 0.0083888,
     $     0.0083691, 0.0083205, 0.0082591, 0.008199, 0.0081509,
     $     0.0081224, 0.008118, 0.0081399, 0.0081887, 0.0082641,
     $     0.0083649, 0.0084897, 0.0086369, 0.0088048, 0.0089913,
     $     0.0091944, 0.009412, 0.0096418, 0.0098811, 0.010127, 0.010377
     $     , 0.010627, 0.010875, 0.011115, 0.011345, 0.011564, 0.011769,
     $     0.011958, 0.012131, 0.012288, 0.01243, 0.01256, 0.012684,
     $     0.012802, 0.012918, 0.013037, 0.013164, 0.013302, 0.013452,
     $     0.013613, 0.013783, 0.013954, 0.01412, 0.01427, 0.014394,
     $     0.014488, 0.01455, 0.014586, 0.014607, 0.014627, 0.014661,
     $     0.014721, 0.014818, 0.014961, 0.015156, 0.015408, 0.015722,
     $     0.016105, 0.016556, 0.01708, 0.017678, 0.018353, 0.019104,
     $     0.019927, 0.020813, 0.021743, 0.022692, 0.02362, 0.02448,
     $     0.025221, 0.025794, 0.026166, 0.026326, 0.026284, 0.02608,
     $     0.025761, 0.025377, 0.024971, 0.024571, 0.024198, 0.023863,
     $     0.02357, 0.023318, 0.023099, 0.022911, 0.022747, 0.022601,
     $     0.022465, 0.022339, 0.022221, 0.022099, 0.021971, 0.021835,
     $     0.021697, 0.021546, 0.021382, 0.021202, 0.020999, 0.020785,
     $     0.020562, 0.020331, 0.020095, 0.019855, 0.019613, 0.019371,
     $     0.019132, 0.018896, 0.018665, 0.01844, 0.018224, 0.018016,
     $     0.017818, 0.017628, 0.017449, 0.017281, 0.017124, 0.016977,
     $     0.016839, 0.016712, 0.016594, 0.016486, 0.016387, 0.016295,
     $     0.016209, 0.016128, 0.016052, 0.01598, 0.01591, 0.015841,
     $     0.015772, 0.015703, 0.015632, 0.015558, 0.015481, 0.0154,
     $     0.015316, 0.015228, 0.015136, 0.015041, 0.014943, 0.014842,
     $     0.014739, 0.014635, 0.014531, 0.014426, 0.014322, 0.01422,
     $     0.014121, 0.014024, 0.013931, 0.013842, 0.013757, 0.013677,
     $     0.013602, 0.013532, 0.013467, 0.013406, 0.013349, 0.013296,
     $     0.013247, 0.013202, 0.01316, 0.013121, 0.013083, 0.013048,
     $     0.013013, 0.01298, 0.012946, 0.012912, 0.012878, 0.012843,
     $     0.012807, 0.012771, 0.012732, 0.012693, 0.012654, 0.012614,
     $     0.012574, 0.012535, 0.012497, 0.01246, 0.012425, 0.012392,
     $     0.012361, 0.012333, 0.012309, 0.012289, 0.012272, 0.01226,
     $     0.012251, 0.012248, 0.012248, 0.012253, 0.012263, 0.012276,
     $     0.012294, 0.012315, 0.012339, 0.012366, 0.012396, 0.012429,
     $     0.012464, 0.0125, 0.012538, 0.012575, 0.012613, 0.012651,
     $     0.012687, 0.012721, 0.012753, 0.012782, 0.012808, 0.01283,
     $     0.012847, 0.01286, 0.012867, 0.012868, 0.012863, 0.012852,
     $     0.012835, 0.012811, 0.012781, 0.011128, 0.0093885, 0.0085326,
     $     0.0081057, 0.0076757, 0.0070688, 0.0063491, 0.0056535,
     $     0.0050768, 0.0046583, 0.0043897, 0.0042411, 0.0041599,
     $     0.0040745, 0.0039233, 0.0036855, 0.0033943, 0.0030999,
     $     0.0028364, 0.002617, 0.0020339, 0.0018386, 0.0015733,
     $     0.0015869, 0.001586, 0.0015729, 0.0015495, 0.0015179,
     $     0.0014798, 0.0014369, 0.0013906, 0.0013425, 0.0012937,
     $     0.0012455, 0.001199, 0.0011554, 0.0011157, 0.0010807/


C     8016, JENDL evaluated data
      nbins(10) = 496
      data ( eg(10,i), i=1, 496) /7162000.0, 7200000.0, 7250000.0,
     $     7300000.0, 7350000.0, 7400000.0, 7450000.0, 7500000.0,
     $     7550000.0, 7600000.0, 7650000.0, 7700000.0, 7750000.0,
     $     7800000.0, 7850000.0, 7900000.0, 7950000.0, 8000000.0,
     $     8050000.0, 8100000.0, 8150000.0, 8200000.0, 8250000.0,
     $     8300000.0, 8350000.0, 8400000.0, 8450000.0, 8500000.0,
     $     8550000.0, 8600000.0, 8650000.0, 8700000.0, 8750000.0,
     $     8800000.0, 8850000.0, 8900000.0, 8950000.0, 9000000.0,
     $     9050000.0, 9100000.0, 9150000.0, 9200000.0, 9250000.0,
     $     9300000.0, 9350000.0, 9400000.0, 9450000.0, 9500000.0,
     $     9550000.0, 9600000.0, 9650000.0, 9700000.0, 9750000.0,
     $     9800000.0, 9850000.0, 9900000.0, 9950000.0, 1.0D7, 1.005D7,
     $     1.01D7, 1.015D7, 1.02D7, 1.025D7, 1.03D7, 1.035D7, 1.04D7,
     $     1.045D7, 1.05D7, 1.055D7, 1.06D7, 1.065D7, 1.07D7, 1.075D7,
     $     1.08D7, 1.085D7, 1.09D7, 1.095D7, 1.1D7, 1.105D7, 1.11D7,
     $     1.115D7, 1.12D7, 1.125D7, 1.13D7, 1.135D7, 1.14D7, 1.145D7,
     $     1.15D7, 1.155D7, 1.16D7, 1.165D7, 1.17D7, 1.175D7, 1.18D7,
     $     1.185D7, 1.19D7, 1.195D7, 1.2D7, 1.205D7, 1.21D7, 1.215D7,
     $     1.22D7, 1.225D7, 1.23D7, 1.235D7, 1.24D7, 1.245D7, 1.25D7,
     $     1.255D7, 1.26D7, 1.265D7, 1.27D7, 1.275D7, 1.28D7, 1.285D7,
     $     1.29D7, 1.295D7, 1.3D7, 1.305D7, 1.31D7, 1.315D7, 1.32D7,
     $     1.325D7, 1.33D7, 1.335D7, 1.34D7, 1.345D7, 1.35D7, 1.355D7,
     $     1.36D7, 1.365D7, 1.37D7, 1.375D7, 1.38D7, 1.385D7, 1.39D7,
     $     1.395D7, 1.4D7, 1.405D7, 1.41D7, 1.415D7, 1.42D7, 1.425D7,
     $     1.43D7, 1.435D7, 1.44D7, 1.445D7, 1.45D7, 1.455D7, 1.46D7,
     $     1.465D7, 1.47D7, 1.475D7, 1.48D7, 1.485D7, 1.49D7, 1.495D7,
     $     1.5D7, 1.505D7, 1.51D7, 1.515D7, 1.52D7, 1.525D7, 1.53D7,
     $     1.535D7, 1.54D7, 1.545D7, 1.55D7, 1.555D7, 1.56D7, 1.565D7,
     $     1.57D7, 1.575D7, 1.58D7, 1.585D7, 1.59D7, 1.595D7, 1.6D7,
     $     1.605D7, 1.61D7, 1.615D7, 1.62D7, 1.625D7, 1.63D7, 1.635D7,
     $     1.64D7, 1.645D7, 1.65D7, 1.655D7, 1.66D7, 1.665D7, 1.67D7,
     $     1.675D7, 1.68D7, 1.685D7, 1.69D7, 1.695D7, 1.7D7, 1.705D7,
     $     1.71D7, 1.715D7, 1.72D7, 1.725D7, 1.73D7, 1.735D7, 1.74D7,
     $     1.745D7, 1.75D7, 1.755D7, 1.76D7, 1.765D7, 1.77D7, 1.775D7,
     $     1.78D7, 1.785D7, 1.79D7, 1.795D7, 1.8D7, 1.805D7, 1.81D7,
     $     1.815D7, 1.82D7, 1.825D7, 1.83D7, 1.835D7, 1.84D7, 1.845D7,
     $     1.85D7, 1.855D7, 1.86D7, 1.865D7, 1.87D7, 1.875D7, 1.88D7,
     $     1.885D7, 1.89D7, 1.895D7, 1.9D7, 1.905D7, 1.91D7, 1.915D7,
     $     1.92D7, 1.925D7, 1.93D7, 1.935D7, 1.94D7, 1.945D7, 1.95D7,
     $     1.955D7, 1.96D7, 1.965D7, 1.97D7, 1.975D7, 1.98D7, 1.985D7,
     $     1.99D7, 1.995D7, 2.0D7, 2.005D7, 2.01D7, 2.015D7, 2.02D7,
     $     2.025D7, 2.03D7, 2.035D7, 2.04D7, 2.045D7, 2.05D7, 2.055D7,
     $     2.06D7, 2.065D7, 2.07D7, 2.075D7, 2.08D7, 2.085D7, 2.09D7,
     $     2.095D7, 2.1D7, 2.105D7, 2.11D7, 2.115D7, 2.12D7, 2.125D7,
     $     2.13D7, 2.135D7, 2.14D7, 2.145D7, 2.15D7, 2.155D7, 2.16D7,
     $     2.165D7, 2.17D7, 2.175D7, 2.18D7, 2.185D7, 2.19D7, 2.195D7,
     $     2.2D7, 2.205D7, 2.21D7, 2.215D7, 2.22D7, 2.225D7, 2.23D7,
     $     2.235D7, 2.24D7, 2.245D7, 2.25D7, 2.255D7, 2.26D7, 2.265D7,
     $     2.27D7, 2.275D7, 2.28D7, 2.285D7, 2.29D7, 2.295D7, 2.3D7,
     $     2.305D7, 2.31D7, 2.315D7, 2.32D7, 2.325D7, 2.33D7, 2.335D7,
     $     2.34D7, 2.345D7, 2.35D7, 2.355D7, 2.36D7, 2.365D7, 2.37D7,
     $     2.375D7, 2.38D7, 2.385D7, 2.39D7, 2.395D7, 2.4D7, 2.405D7,
     $     2.41D7, 2.415D7, 2.42D7, 2.425D7, 2.43D7, 2.435D7, 2.44D7,
     $     2.445D7, 2.45D7, 2.455D7, 2.46D7, 2.465D7, 2.47D7, 2.475D7,
     $     2.48D7, 2.485D7, 2.49D7, 2.495D7, 2.5D7, 2.505D7, 2.51D7,
     $     2.515D7, 2.52D7, 2.525D7, 2.53D7, 2.535D7, 2.54D7, 2.545D7,
     $     2.55D7, 2.555D7, 2.56D7, 2.565D7, 2.57D7, 2.575D7, 2.58D7,
     $     2.585D7, 2.59D7, 2.595D7, 2.6D7, 2.605D7, 2.61D7, 2.615D7,
     $     2.62D7, 2.625D7, 2.63D7, 2.635D7, 2.64D7, 2.645D7, 2.65D7,
     $     2.655D7, 2.66D7, 2.665D7, 2.67D7, 2.675D7, 2.68D7, 2.685D7,
     $     2.69D7, 2.695D7, 2.7D7, 2.705D7, 2.71D7, 2.715D7, 2.72D7,
     $     2.725D7, 2.73D7, 2.735D7, 2.74D7, 2.745D7, 2.75D7, 2.755D7,
     $     2.76D7, 2.765D7, 2.77D7, 2.775D7, 2.78D7, 2.785D7, 2.79D7,
     $     2.795D7, 2.8D7, 2.805D7, 2.81D7, 2.815D7, 2.82D7, 2.825D7,
     $     2.83D7, 2.835D7, 2.84D7, 2.845D7, 2.85D7, 2.855D7, 2.86D7,
     $     2.865D7, 2.87D7, 2.875D7, 2.88D7, 2.885D7, 2.89D7, 2.895D7,
     $     2.9D7, 2.905D7, 2.91D7, 2.915D7, 2.92D7, 2.925D7, 2.93D7,
     $     2.935D7, 2.94D7, 2.945D7, 2.95D7, 2.955D7, 2.96D7, 2.965D7,
     $     2.97D7, 2.975D7, 2.98D7, 2.985D7, 2.99D7, 2.995D7, 3.0D7,
     $     3.1D7, 3.2D7, 3.3D7, 3.4D7, 3.5D7, 3.6D7, 3.7D7, 3.8D7, 3.9D7
     $     , 4.0D7, 4.1D7, 4.2D7, 4.3D7, 4.4D7, 4.5D7, 4.6D7, 4.7D7,
     $     4.8D7, 4.9D7, 5.0D7, 5.5D7, 6.0D7, 6.5D7, 7.0D7, 7.5D7, 8.0D7
     $     , 8.5D7, 9.0D7, 9.5D7, 1.0D8, 1.05D8, 1.1D8, 1.15D8, 1.2D8,
     $     1.25D8, 1.3D8, 1.35D8, 1.4D8/

      data ( gdr(10,i), i=1, 496) /1.692D-6, 1.6966D-6, 1.701D-6,
     $     1.7053D-6, 1.7096D-6, 1.7139D-6, 1.7182D-6, 1.7227D-6,
     $     1.7274D-6, 1.7321D-6, 1.737D-6, 1.742D-6, 1.7472D-6, 1.7527D
     $     -6, 1.7583D-6, 1.7641D-6, 1.7701D-6, 1.7769D-6, 1.784D-6,
     $     1.7914D-6, 1.799D-6, 1.8162D-6, 1.8366D-6, 1.8574D-6, 1.8787D
     $     -6, 1.9004D-6, 1.9225D-6, 1.9452D-6, 1.9684D-6, 1.9922D-6,
     $     2.0166D-6, 2.0419D-6, 2.0679D-6, 2.0949D-6, 2.1231D-6,
     $     2.1525D-6, 2.1836D-6, 2.2166D-6, 2.252D-6, 2.2906D-6, 2.3331D
     $     -6, 2.4344D-6, 2.5668D-6, 2.7191D-6, 2.8977D-6, 3.1081D-6,
     $     3.3499D-6, 3.6084D-6, 3.847D-6, 4.0192D-6, 4.1014D-6, 4.1102D
     $     -6, 4.0831D-6, 4.051D-6, 4.0302D-6, 4.0261D-6, 4.0386D-6,
     $     4.0656D-6, 4.1048D-6, 4.1538D-6, 4.2109D-6, 4.2747D-6, 4.344D
     $     -6, 4.418D-6, 4.496D-6, 4.5776D-6, 4.6623D-6, 4.7499D-6,
     $     4.8401D-6, 4.9329D-6, 5.028D-6, 5.1255D-6, 5.2252D-6, 5.3273D
     $     -6, 5.4317D-6, 5.5384D-6, 5.6477D-6, 5.7595D-6, 5.8741D-6,
     $     5.9917D-6, 6.1124D-6, 6.2218D-6, 6.3302D-6, 6.4425D-6,
     $     6.5593D-6, 6.6811D-6, 6.8085D-6, 6.9426D-6, 7.0844D-6,
     $     7.2353D-6, 7.397D-6, 7.5719D-6, 7.7628D-6, 7.9736D-6, 8.2092D
     $     -6, 8.4763D-6, 8.7835D-6, 9.1423D-6, 9.5673D-6, 1.0076D-5,
     $     1.548D-5, 2.05D-5, 2.4614D-5, 2.8468D-5, 4.5412D-5, 8.0676D-5
     $     , 1.1275D-4, 1.3734D-4, 1.9836D-4, 2.963D-4, 3.8157D-4,
     $     4.7139D-4, 6.3036D-4, 9.1646D-4, 0.0013644, 0.0021468,
     $     0.0036896, 0.0067724, 0.011407, 0.013383, 0.010457, 0.0071436
     $     , 0.0050116, 0.0037273, 0.0029248, 0.0023963, 0.002031,
     $     0.0017679, 0.001572, 0.001422, 0.0013046, 0.0012109,
     $     0.0011351, 0.0010729, 0.0010215, 9.7861D-4, 9.4265D-4,
     $     9.1239D-4, 8.8685D-4, 8.6531D-4, 8.4394D-4, 8.2189D-4,
     $     8.0294D-4, 7.8667D-4, 7.7276D-4, 7.609D-4, 7.5088D-4, 7.4249D
     $     -4, 7.3558D-4, 7.3D-4, 7.2564D-4, 7.2239D-4, 7.2019D-4,
     $     7.1896D-4, 7.1865D-4, 7.1921D-4, 7.2061D-4, 7.2283D-4,
     $     7.2586D-4, 7.2969D-4, 7.3434D-4, 7.3981D-4, 7.4615D-4,
     $     7.5342D-4, 7.6169D-4, 7.7107D-4, 7.8168D-4, 7.9371D-4,
     $     8.0738D-4, 8.2302D-4, 8.4103D-4, 8.981D-4, 9.5855D-4,
     $     0.0010228, 0.0010943, 0.0011638, 0.0012402, 0.0013297,
     $     0.0014328, 0.0015424, 0.0016454, 0.0017232, 0.0017638,
     $     0.0017678, 0.0017535, 0.0017413, 0.0017449, 0.0017701,
     $     0.0018209, 0.0018979, 0.0019967, 0.0021031, 0.0021951,
     $     0.0022466, 0.002246, 0.0022074, 0.0021602, 0.0021323,
     $     0.0021467, 0.0022321, 0.0024514, 0.0029994, 0.0045885,
     $     0.0072663, 0.0052079, 0.0044296, 0.0050276, 0.0055636,
     $     0.0047666, 0.0036513, 0.0029418, 0.0025456, 0.0023211,
     $     0.0021897, 0.0021111, 0.0020645, 0.0020384, 0.0020261,
     $     0.0020237, 0.0020287, 0.002039, 0.0020533, 0.0020716,
     $     0.0020936, 0.0021189, 0.0021477, 0.0021799, 0.002216,
     $     0.0022566, 0.0023025, 0.0023553, 0.0024171, 0.0024919,
     $     0.0025862, 0.0027117, 0.0028912, 0.0031695, 0.0036277,
     $     0.0043163, 0.004821, 0.0045463, 0.0041104, 0.0039568,
     $     0.0040919, 0.0045423, 0.0054547, 0.0069201, 0.0079203,
     $     0.0070036, 0.0056176, 0.0047589, 0.0043146, 0.004099,
     $     0.0040101, 0.0039971, 0.0040346, 0.0041101, 0.0042174,
     $     0.0043539, 0.0045195, 0.0047138, 0.0049369, 0.0051882,
     $     0.0054628, 0.00575, 0.0060333, 0.0062933, 0.0065153,
     $     0.0066986, 0.006862, 0.0070403, 0.0072797, 0.0076357,
     $     0.0081754, 0.0089776, 0.010085, 0.011314, 0.01205, 0.01184,
     $     0.011157, 0.010616, 0.010405, 0.010494, 0.01082, 0.011332,
     $     0.011998, 0.012784, 0.013652, 0.014553, 0.015428, 0.016215,
     $     0.016871, 0.017395, 0.017832, 0.018263, 0.018789, 0.019507,
     $     0.020498, 0.021808, 0.023432, 0.025273, 0.027119, 0.028663,
     $     0.029602, 0.029803, 0.029381, 0.0286, 0.027698, 0.026798,
     $     0.02593, 0.025102, 0.024352, 0.023752, 0.023391, 0.023356,
     $     0.023698, 0.024435, 0.025506, 0.026826, 0.028153, 0.029153,
     $     0.029507, 0.029079, 0.028, 0.026568, 0.025079, 0.023731,
     $     0.022618, 0.021767, 0.02118, 0.020854, 0.020781, 0.020946,
     $     0.021309, 0.021847, 0.022499, 0.023202, 0.023922, 0.024687,
     $     0.025589, 0.026731, 0.028044, 0.029168, 0.029439, 0.028493,
     $     0.026742, 0.024895, 0.023377, 0.022296, 0.021607, 0.021251,
     $     0.021169, 0.021314, 0.02164, 0.022104, 0.022582, 0.022946,
     $     0.023073, 0.022916, 0.022519, 0.021993, 0.021454, 0.020984,
     $     0.020626, 0.020401, 0.02031, 0.020349, 0.020499, 0.020723,
     $     0.020953, 0.021096, 0.021061, 0.020808, 0.020369, 0.019827,
     $     0.019263, 0.018735, 0.01827, 0.017873, 0.01754, 0.017263,
     $     0.017031, 0.016836, 0.016669, 0.016524, 0.016395, 0.016277,
     $     0.016167, 0.01606, 0.015954, 0.015847, 0.015737, 0.015621,
     $     0.0155, 0.01537, 0.015233, 0.015089, 0.014939, 0.014781,
     $     0.014618, 0.014448, 0.014273, 0.014095, 0.013913, 0.01373,
     $     0.013546, 0.013362, 0.013179, 0.012998, 0.01282, 0.012644,
     $     0.012473, 0.012306, 0.012144, 0.011988, 0.011837, 0.011691,
     $     0.011552, 0.011419, 0.011293, 0.011172, 0.011058, 0.01095,
     $     0.010848, 0.010753, 0.010663, 0.010579, 0.010501, 0.010428,
     $     0.010361, 0.0103, 0.010243, 0.010191, 0.010144, 0.010101,
     $     0.010062, 0.010027, 0.0099969, 0.00997, 0.0099467, 0.0099269,
     $     0.0099103, 0.0098969, 0.0098863, 0.0098785, 0.0098732,
     $     0.0098704, 0.0098698, 0.0098712, 0.0098745, 0.0098791,
     $     0.0098849, 0.009892, 0.0099003, 0.0099098, 0.0099201,
     $     0.0099312, 0.0099428, 0.010084, 0.0095804, 0.0084646,
     $     0.0073048, 0.0064492, 0.0059401, 0.0057255, 0.005716,
     $     0.0058207, 0.0059382, 0.0059677, 0.0058399, 0.0055433,
     $     0.0051286, 0.0046692, 0.0042245, 0.0038293, 0.0034941,
     $     0.0032175, 0.0029929, 0.00237, 0.0021361, 0.0017972,
     $     0.0018125, 0.0018114, 0.0017964, 0.0017697, 0.0017336,
     $     0.0016901, 0.0016411, 0.0015883, 0.0015332, 0.0014775,
     $     0.0014224, 0.0013694, 0.0013196, 0.0012742, 0.0012343/

c$$$C     8016, ENDF-B/VII.0 evaluated data
c$$$      nbins(10) = 96
c$$$      data ( eg(10,i), i=1, 96) / 7161900.0, 1.21275D7, 1.25D7, 1.3D7,
c$$$     $     1.325D7, 1.35D7, 1.4D7, 1.45D7, 1.5D7, 1.55D7, 1.56638D7,
c$$$     $     1.6D7, 1.65D7, 1.7D7, 1.75D7, 1.8D7, 1.85D7, 1.9D7, 1.95D7,
c$$$     $     2.0D7, 2.05D7, 2.1D7, 2.125D7, 2.15D7, 2.175D7, 2.2D7,
c$$$     $     2.225D7, 2.25D7, 2.275D7, 2.3D7, 2.325D7, 2.35D7, 2.375D7,
c$$$     $     2.4D7, 2.425D7, 2.45D7, 2.475D7, 2.5D7, 2.55D7, 2.6D7, 2.65D7
c$$$     $     , 2.7D7, 2.75D7, 2.8D7, 2.85D7, 2.9D7, 2.95D7, 3.0D7, 3.25D7,
c$$$     $     3.5D7, 3.75D7, 4.0D7, 4.25D7, 4.5D7, 4.75D7, 5.0D7, 5.25D7,
c$$$     $     5.5D7, 5.75D7, 6.0D7, 6.25D7, 6.5D7, 6.75D7, 7.0D7, 7.25D7,
c$$$     $     7.5D7, 7.75D7, 8.0D7, 8.25D7, 8.5D7, 8.75D7, 9.0D7, 9.25D7,
c$$$     $     9.5D7, 9.75D7, 1.0D8, 1.025D8, 1.05D8, 1.075D8, 1.1D8,
c$$$     $     1.125D8, 1.15D8, 1.175D8, 1.2D8, 1.225D8, 1.25D8, 1.275D8,
c$$$     $     1.3D8, 1.325D8, 1.35D8, 1.375D8, 1.4D8, 1.425D8, 1.45D8,
c$$$     $     1.475D8, 1.5D8/
c$$$
c$$$      data ( gdr(10,i), i=1, 96) / 0.0, 0.0, 0.0, 0.0, 0.0025, 0.00204,
c$$$     $     0.0, 1.0D-4, 1.0D-4, 1.0D-4, 1.0D-4, 1.0D-4, 1.0D-4, 9.99998D
c$$$     $     -5, 0.0, 1.0D-5, 0.0011, 8.999999D-5, 8.999999D-5, 0.00295,
c$$$     $     0.00275, 0.00457, 0.006905, 0.00448, 0.00312, 0.009039999,
c$$$     $     0.0097, 0.0118, 0.01789, 0.01755, 0.015715, 0.01525, 0.01367,
c$$$     $     0.01162, 0.012095, 0.01002, 0.013875, 0.01087, 0.009799999,
c$$$     $     0.00919, 0.0082, 0.00712, 0.0057, 0.0048, 0.0054, 0.00502,
c$$$     $     0.0053, 0.005419998, 0.005796018, 0.0060125, 0.00597394,
c$$$     $     0.005664772, 0.005105959, 0.004542546, 0.004167533, 0.003791,
c$$$     $     0.00307301, 0.002382815, 0.002091185, 0.001985684,
c$$$     $     0.001991926, 0.002054, 0.002338106, 0.002558663, 0.001901061,
c$$$     $     0.001082044, 7.333554D-4, 5.920832D-4, 7.012359D-4, 8.057648D
c$$$     $     -4, 7.979377D-4, 7.90439D-4, 9.532238D-4, 0.001138284,
c$$$     $     9.054484D-4, 6.659999D-4, 9.257685D-4, 0.001411012,
c$$$     $     0.001870532, 0.002053132, 0.002012404, 0.001970306,
c$$$     $     0.00192819, 0.00188741, 0.001847906, 0.001808842, 0.001770404
c$$$     $     , 0.001732777, 0.001818051, 0.002057859, 0.002316702,
c$$$     $     0.002459084, 0.002421418, 0.002245361, 0.001953279,
c$$$     $     0.001567535/


C     9019, JENDL evaluated data
      nbins(11) = 213
      data ( eg(11,i), i=1, 213) /4013000.0, 4200000.0, 4400000.0,
     $     4600000.0, 4800000.0, 5000000.0, 5200000.0, 5400000.0,
     $     5600000.0, 5800000.0, 6000000.0, 6200000.0, 6400000.0,
     $     6600000.0, 6800000.0, 7000000.0, 7200000.0, 7400000.0,
     $     7600000.0, 7800000.0, 8000000.0, 8200000.0, 8400000.0,
     $     8600000.0, 8800000.0, 9000000.0, 9200000.0, 9400000.0,
     $     9600000.0, 9800000.0, 1.0D7, 1.02D7, 1.04D7, 1.06D7, 1.08D7,
     $     1.1D7, 1.12D7, 1.14D7, 1.16D7, 1.18D7, 1.2D7, 1.22D7, 1.24D7,
     $     1.26D7, 1.28D7, 1.3D7, 1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7,
     $     1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7
     $     , 1.58D7, 1.6D7, 1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7,
     $     1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7
     $     , 1.88D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7,
     $     2.05D7, 2.1D7, 2.15D7, 2.2D7, 2.25D7, 2.3D7, 2.35D7, 2.4D7,
     $     2.45D7, 2.5D7, 2.55D7, 2.6D7, 2.65D7, 2.7D7, 2.75D7, 2.8D7,
     $     2.85D7, 2.9D7, 2.95D7, 3.0D7, 3.05D7, 3.1D7, 3.15D7, 3.2D7,
     $     3.25D7, 3.3D7, 3.35D7, 3.4D7, 3.45D7, 3.5D7, 3.55D7, 3.6D7,
     $     3.65D7, 3.7D7, 3.75D7, 3.8D7, 3.85D7, 3.9D7, 3.95D7, 4.0D7,
     $     4.05D7, 4.1D7, 4.15D7, 4.2D7, 4.25D7, 4.3D7, 4.35D7, 4.4D7,
     $     4.45D7, 4.5D7, 4.55D7, 4.6D7, 4.65D7, 4.7D7, 4.75D7, 4.8D7,
     $     4.85D7, 4.9D7, 4.95D7, 5.0D7, 5.05D7, 5.1D7, 5.15D7, 5.2D7,
     $     5.25D7, 5.3D7, 5.35D7, 5.4D7, 5.45D7, 5.5D7, 5.55D7, 5.6D7,
     $     5.65D7, 5.7D7, 5.75D7, 5.8D7, 5.85D7, 5.9D7, 5.95D7, 6.0D7,
     $     6.05D7, 6.1D7, 6.15D7, 6.2D7, 6.25D7, 6.3D7, 6.35D7, 6.4D7,
     $     6.45D7, 6.5D7, 6.55D7, 6.6D7, 6.65D7, 6.7D7, 6.75D7, 6.8D7,
     $     6.85D7, 6.9D7, 6.95D7, 7.0D7, 7.05D7, 7.1D7, 7.15D7, 7.2D7,
     $     7.25D7, 7.3D7, 7.35D7, 7.4D7, 7.45D7, 7.5D7, 7.55D7, 7.6D7,
     $     7.65D7, 7.7D7, 7.75D7, 7.8D7, 7.85D7, 7.9D7, 7.95D7, 8.0D7,
     $     8.5D7, 9.0D7, 9.5D7, 1.0D8, 1.05D8, 1.1D8, 1.15D8, 1.2D8,
     $     1.25D8, 1.3D8, 1.35D8, 1.4D8/

      data ( gdr(11,i), i=1, 213) /0.0, 1.8543D-5, 4.1119D-5, 6.3058D-5,
     $     8.5314D-5, 1.0796D-4, 1.2822D-4, 1.5012D-4, 1.7405D-4,
     $     2.0013D-4, 2.2845D-4, 2.5133D-4, 2.7519D-4, 3.0064D-4,
     $     3.2781D-4, 3.5682D-4, 3.8784D-4, 4.2104D-4, 4.5664D-4, 4.949D
     $     -4, 5.3614D-4, 5.7364D-4, 6.1384D-4, 6.5778D-4, 7.063D-4,
     $     7.6059D-4, 8.2247D-4, 8.9484D-4, 9.8272D-4, 0.0010954,
     $     0.0012519, 0.0014885, 0.0019253, 0.002899, 0.0048956,
     $     0.0049179, 0.003267, 0.0027273, 0.0030918, 0.0091917,
     $     0.0092999, 0.0090821, 0.0046047, 0.004304, 0.0040116,
     $     0.0053745, 0.0041918, 0.0044915, 0.0067281, 0.010494,
     $     0.007217, 0.0057254, 0.0063164, 0.0065598, 0.0063116,
     $     0.0069083, 0.0079403, 0.0091034, 0.0099502, 0.010088,
     $     0.0096303, 0.0090386, 0.0086482, 0.0085671, 0.0088002,
     $     0.0093324, 0.010129, 0.011116, 0.01213, 0.012916, 0.013226,
     $     0.013004, 0.012433, 0.011774, 0.01122, 0.010861, 0.010722,
     $     0.010822, 0.011193, 0.011902, 0.013067, 0.016724, 0.014907,
     $     0.014589, 0.015467, 0.016574, 0.01748, 0.017935, 0.017852,
     $     0.017313, 0.016524, 0.015729, 0.015165, 0.015058, 0.015583,
     $     0.016565, 0.01683, 0.015524, 0.013776, 0.012465, 0.011605,
     $     0.011009, 0.010527, 0.010069, 0.0095881, 0.0090732, 0.008531,
     $     0.0079778, 0.007432, 0.0069078, 0.0064176, 0.0059675, 0.00556
     $     , 0.0051945, 0.0048677, 0.0045777, 0.004321, 0.0040939,
     $     0.0038932, 0.0037148, 0.0035567, 0.0034165, 0.0032919,
     $     0.003181, 0.0030816, 0.0029928, 0.0029133, 0.002842,
     $     0.0027779, 0.0027195, 0.0026668, 0.0026192, 0.0025762,
     $     0.0025373, 0.002502, 0.0024699, 0.0024407, 0.0024142, 0.00239
     $     , 0.0023679, 0.0023477, 0.0023292, 0.0023123, 0.0022967,
     $     0.0022825, 0.0022693, 0.0022572, 0.002246, 0.0022356,
     $     0.0022256, 0.0022164, 0.0022077, 0.0021997, 0.0021922,
     $     0.0021852, 0.0021786, 0.0021724, 0.0021666, 0.0021611,
     $     0.0021559, 0.0021509, 0.0021462, 0.0021416, 0.0021373,
     $     0.0021331, 0.002129, 0.0021251, 0.0021213, 0.0021175,
     $     0.0021137, 0.00211, 0.0021063, 0.0021027, 0.0020991,
     $     0.0020956, 0.0020921, 0.0020886, 0.0020851, 0.0020816,
     $     0.0020781, 0.0020746, 0.002071, 0.0020675, 0.0020639,
     $     0.0020603, 0.0020567, 0.002053, 0.0020493, 0.0020455,
     $     0.0020417, 0.0020378, 0.0020338, 0.0020298, 0.0020258,
     $     0.0020217, 0.0020175, 0.0020134, 0.0020091, 0.0020048,
     $     0.0020005, 0.0019961, 0.0019489, 0.0018963, 0.0018505,
     $     0.0017988, 0.0017431, 0.001685, 0.0016262, 0.0015682,
     $     0.0015122, 0.0014597, 0.0014118, 0.0013698/


C     11023, JENDL evaluated data
      nbins(12) = 188
      data ( eg(12,i), i=1, 188) /1.26D7, 1.28D7, 1.3D7, 1.32D7,
     $     1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7, 1.48D7
     $     ,1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7, 1.62D7, 1.64D7
     $     ,1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7
     $     ,1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7, 1.92D7, 1.94D7,
     $     1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7, 2.08D7, 2.1D7
     $     ,2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7, 2.22D7, 2.24D7,
     $     2.26D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7, 2.38D7, 2.4D7
     $     ,2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7, 2.52D7, 2.54D7,
     $     2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7, 2.68D7, 2.7D7
     $     ,2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7, 2.82D7, 2.84D7,
     $     2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7, 2.98D7, 3.0D7
     $     ,3.02D7, 3.04D7, 3.06D7, 3.08D7, 3.1D7, 3.12D7, 3.14D7,
     $     3.16D7, 3.18D7, 3.2D7, 3.22D7, 3.24D7, 3.26D7, 3.28D7, 3.3D7
     $     ,3.32D7, 3.34D7, 3.36D7, 3.38D7, 3.4D7, 3.42D7, 3.44D7,
     $     3.46D7, 3.48D7, 3.5D7, 3.52D7, 3.54D7, 3.56D7, 3.58D7, 3.6D7
     $     ,3.62D7, 3.64D7, 3.66D7, 3.68D7, 3.7D7, 3.72D7, 3.74D7,
     $     3.76D7, 3.78D7, 3.8D7, 3.82D7, 3.84D7, 3.86D7, 3.88D7, 3.9D7
     $     ,3.92D7, 3.94D7, 3.96D7, 3.98D7, 4.0D7, 4.2D7, 4.4D7, 4.6D7
     $     ,4.8D7, 5.0D7, 5.2D7, 5.4D7, 5.6D7, 5.8D7, 6.0D7, 6.2D7,
     $     6.4D7, 6.6D7, 6.8D7, 7.0D7, 7.2D7, 7.4D7, 7.6D7, 7.8D7, 8.0D7
     $     ,8.2D7, 8.4D7, 8.6D7, 8.8D7, 9.0D7, 9.2D7, 9.4D7, 9.6D7,
     $     9.8D7, 1.0D8, 1.02D8, 1.04D8, 1.06D8, 1.08D8, 1.1D8, 1.12D8
     $     ,1.14D8, 1.16D8, 1.18D8, 1.2D8, 1.22D8, 1.24D8, 1.26D8,
     $     1.28D8, 1.3D8, 1.32D8, 1.34D8, 1.36D8, 1.38D8, 1.4D8/

      data ( gdr(12,i), i=1, 188) /7.137D-4, 0.0013, 0.001779, 0.002209,
     $     0.002638, 0.003112, 0.003662, 0.004343, 0.005182, 0.006099,
     $     0.006991, 0.007821, 0.008625, 0.009475, 0.01056, 0.01202,
     $     0.01379, 0.01554, 0.017, 0.01812, 0.01892, 0.01949, 0.01997,
     $     0.02049, 0.02097, 0.02128, 0.02145, 0.0216, 0.02186, 0.02215,
     $     0.02238, 0.02253, 0.02265, 0.02277, 0.02289, 0.02301, 0.02312
     $     , 0.02324, 0.02338, 0.02354, 0.02369, 0.02383, 0.02397,
     $     0.02378, 0.02378, 0.02378, 0.02377, 0.02375, 0.02373, 0.02368
     $     , 0.02357, 0.02344, 0.02337, 0.02341, 0.02346, 0.02339,
     $     0.02312, 0.02296, 0.02298, 0.02311, 0.02324, 0.02332, 0.02339
     $     , 0.02332, 0.0232, 0.02305, 0.02286, 0.02225, 0.02204,
     $     0.02183, 0.0216, 0.02134, 0.02094, 0.02039, 0.01976, 0.01908,
     $     0.01842, 0.01774, 0.017, 0.01626, 0.01556, 0.01489, 0.01423,
     $     0.01356, 0.01287, 0.01217, 0.01146, 0.01077, 0.009868,
     $     0.009265, 0.008705, 0.008194, 0.007738, 0.007342, 0.006998,
     $     0.006698, 0.006431, 0.00619, 0.005966, 0.005757, 0.005562,
     $     0.005379, 0.005206, 0.005042, 0.004887, 0.004741, 0.004604,
     $     0.004475, 0.004354, 0.004242, 0.004142, 0.004051, 0.003975,
     $     0.003909, 0.003852, 0.003798, 0.003746, 0.003691, 0.00363,
     $     0.003564, 0.003497, 0.003429, 0.003365, 0.003304, 0.00325,
     $     0.003204, 0.003161, 0.003122, 0.003085, 0.003049, 0.003016,
     $     0.002983, 0.00295, 0.002917, 0.002885, 0.002854, 0.002826,
     $     0.002805, 0.002638, 0.002535, 0.002436, 0.002371, 0.002319,
     $     0.002272, 0.002236, 0.002206, 0.002177, 0.002148, 0.002134,
     $     0.002117, 0.002102, 0.002094, 0.002088, 0.00207, 0.002044,
     $     0.002029, 0.002016, 0.001996, 0.001982, 0.001967, 0.001936,
     $     0.001909, 0.001901, 0.001874, 0.001839, 0.001812, 0.001791,
     $     0.001767, 0.001733, 0.001708, 0.001684, 0.001654, 0.001622,
     $     0.001591, 0.00157, 0.001548, 0.001521, 0.001495, 0.001468,
     $     0.001439, 0.001418, 0.001397, 0.001371, 0.001348, 0.001331,
     $     0.001312, 0.001291, 0.001272/


C     12024, JENDL evaluated data
      nbins(13) = 266
      data ( eg(13,i), i=1, 266) /1.1D7, 1.11D7, 1.12D7, 1.13D7,
     $     1.14D7, 1.15D7, 1.16D7, 1.17D7, 1.18D7, 1.19D7, 1.2D7, 1.21D7
     $     ,1.22D7, 1.23D7, 1.24D7, 1.25D7, 1.26D7, 1.27D7, 1.28D7
     $     ,1.29D7, 1.3D7, 1.31D7, 1.32D7, 1.33D7, 1.34D7, 1.35D7,
     $     1.36D7, 1.37D7, 1.38D7, 1.39D7, 1.4D7, 1.41D7, 1.42D7, 1.43D7
     $     ,1.44D7, 1.45D7, 1.46D7, 1.47D7, 1.48D7, 1.49D7, 1.5D7,
     $     1.51D7, 1.52D7, 1.53D7, 1.54D7, 1.55D7, 1.56D7, 1.57D7,
     $     1.58D7,1.59D7, 1.6D7, 1.61D7, 1.62D7, 1.63D7, 1.64D7, 1.65D7,
     $     1.66D7, 1.67D7, 1.68D7, 1.69D7, 1.7D7, 1.71D7, 1.72D7, 1.73D7
     $     ,1.74D7, 1.75D7, 1.76D7, 1.77D7, 1.78D7, 1.79D7, 1.8D7,
     $     1.81D7, 1.82D7, 1.83D7, 1.84D7, 1.85D7, 1.86D7, 1.87D7,
     $     1.88D7,1.89D7, 1.9D7, 1.91D7, 1.92D7, 1.93D7, 1.94D7, 1.95D7,
     $     1.96D7, 1.97D7, 1.98D7, 1.99D7, 2.0D7, 2.01D7, 2.02D7, 2.03D7
     $     ,2.04D7, 2.05D7, 2.06D7, 2.07D7, 2.08D7, 2.09D7, 2.1D7,
     $     2.11D7, 2.12D7, 2.13D7, 2.14D7, 2.15D7, 2.16D7, 2.17D7,
     $     2.18D7,2.19D7, 2.2D7, 2.21D7, 2.22D7, 2.23D7, 2.24D7, 2.25D7,
     $     2.26D7, 2.27D7, 2.28D7, 2.29D7, 2.3D7, 2.31D7, 2.32D7, 2.33D7
     $     ,2.34D7, 2.35D7, 2.36D7, 2.37D7, 2.38D7, 2.39D7, 2.4D7,
     $     2.41D7, 2.42D7, 2.43D7, 2.44D7, 2.45D7, 2.46D7, 2.47D7,
     $     2.48D7,2.49D7, 2.5D7, 2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7,
     $     2.62D7,2.64D7, 2.66D7, 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7,
     $     2.78D7, 2.8D7, 2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7
     $     ,2.94D7, 2.96D7, 2.98D7, 3.0D7, 3.02D7, 3.04D7, 3.06D7,
     $     3.08D7, 3.1D7, 3.12D7, 3.14D7, 3.16D7, 3.18D7, 3.2D7, 3.22D7
     $     ,3.24D7, 3.26D7, 3.28D7, 3.3D7, 3.32D7, 3.34D7, 3.36D7,
     $     3.38D7, 3.4D7, 3.42D7, 3.44D7, 3.46D7, 3.48D7, 3.5D7, 3.52D7
     $     ,3.54D7, 3.56D7, 3.58D7, 3.6D7, 3.62D7, 3.64D7, 3.66D7,
     $     3.68D7, 3.7D7, 3.72D7, 3.74D7, 3.76D7, 3.78D7, 3.8D7, 3.82D7
     $     ,3.84D7, 3.86D7, 3.88D7, 3.9D7, 3.92D7, 3.94D7, 3.96D7,
     $     3.98D7, 4.0D7, 4.2D7, 4.4D7, 4.6D7, 4.8D7, 5.0D7, 5.2D7,
     $     5.4D7,5.6D7, 5.8D7, 6.0D7, 6.2D7, 6.4D7, 6.6D7, 6.8D7, 7.0D7,
     $     7.2D7, 7.4D7, 7.6D7, 7.8D7, 8.0D7, 8.2D7, 8.4D7, 8.6D7, 8.8D7
     $     ,9.0D7, 9.2D7, 9.4D7, 9.6D7, 9.8D7, 1.0D8, 1.02D8, 1.04D8
     $     ,1.06D8, 1.08D8, 1.1D8, 1.12D8, 1.14D8, 1.16D8, 1.18D8, 1.2D8
     $     ,1.22D8, 1.24D8, 1.26D8, 1.28D8, 1.3D8, 1.32D8, 1.34D8,
     $     1.36D8, 1.38D8, 1.4D8/

      data ( gdr(13,i), i=1, 266) /0.003062, 0.003166, 0.003272, 0.00338
     $     , 0.00349, 0.003602, 0.003716, 0.003833, 0.003952, 0.004074,
     $     0.004198, 0.004325, 0.004455, 0.004588, 0.004724, 0.004863,
     $     0.005005, 0.005151, 0.0053, 0.005453, 0.00561, 0.00577,
     $     0.005935, 0.006104, 0.006277, 0.006454, 0.006637, 0.006824,
     $     0.007016, 0.007213, 0.007415, 0.007623, 0.007836, 0.008055,
     $     0.00828, 0.008512, 0.008749, 0.008993, 0.009243, 0.0095,
     $     0.009764, 0.01004, 0.01031, 0.0106, 0.01089, 0.01119, 0.0115,
     $     0.01182, 0.01215, 0.01248, 0.01282, 0.01317, 0.01353, 0.01389
     $     , 0.01427, 0.01465, 0.01504, 0.01544, 0.01585, 0.01627,
     $     0.01669, 0.01712, 0.01756, 0.018, 0.01845, 0.01891, 0.01937,
     $     0.01984, 0.02031, 0.02078, 0.02126, 0.02174, 0.02221, 0.02269
     $     , 0.02317, 0.02364, 0.02411, 0.02458, 0.02503, 0.02548,
     $     0.02592, 0.02635, 0.02677, 0.02717, 0.02756, 0.02793, 0.02828
     $     , 0.02862, 0.02893, 0.02922, 0.02949, 0.02974, 0.02996,
     $     0.03015, 0.03032, 0.03047, 0.03058, 0.03067, 0.03073, 0.03077
     $     , 0.03077, 0.03076, 0.03071, 0.03064, 0.03055, 0.03043,
     $     0.03028, 0.03012, 0.02994, 0.02973, 0.02951, 0.02927, 0.02902
     $     , 0.02875, 0.02846, 0.02817, 0.02786, 0.02755, 0.02722,
     $     0.02689, 0.02655, 0.02621, 0.02586, 0.02551, 0.02516, 0.0248,
     $     0.02444, 0.02409, 0.02373, 0.02338, 0.02302, 0.02267, 0.02233
     $     , 0.02198, 0.02164, 0.0213, 0.02097, 0.02064, 0.02031,
     $     0.01999, 0.01968, 0.01906, 0.01847, 0.0179, 0.01734, 0.01681,
     $     0.01631, 0.01582, 0.01535, 0.01491, 0.01448, 0.01407, 0.01368
     $     , 0.01331, 0.01295, 0.01261, 0.01229, 0.01198, 0.01169,
     $     0.0114, 0.01113, 0.01087, 0.01063, 0.01039, 0.01017, 0.009951
     $     , 0.009745, 0.009547, 0.009359, 0.009178, 0.009005, 0.008839,
     $     0.00868, 0.008528, 0.008382, 0.008242, 0.008107, 0.007978,
     $     0.007854, 0.007736, 0.007621, 0.007512, 0.007406, 0.007305,
     $     0.007207, 0.007113, 0.007023, 0.006936, 0.006852, 0.006772,
     $     0.006694, 0.006619, 0.006547, 0.006477, 0.00641, 0.006346,
     $     0.006283, 0.006223, 0.006165, 0.006109, 0.006054, 0.006002,
     $     0.005952, 0.005903, 0.005855, 0.00581, 0.005765, 0.005723,
     $     0.005681, 0.005641, 0.005603, 0.005565, 0.005529, 0.005494,
     $     0.00546, 0.005427, 0.005148, 0.004941, 0.004784, 0.004664,
     $     0.004569, 0.004493, 0.004431, 0.004378, 0.004331, 0.004289,
     $     0.00425, 0.004213, 0.004177, 0.00414, 0.004104, 0.004066,
     $     0.004028, 0.003989, 0.003948, 0.003907, 0.003863, 0.003819,
     $     0.003774, 0.003727, 0.003679, 0.003631, 0.003581, 0.003531,
     $     0.00348, 0.003429, 0.003377, 0.003325, 0.003273, 0.003221,
     $     0.00317, 0.003118, 0.003067, 0.003017, 0.002967, 0.002918,
     $     0.002871, 0.002824, 0.002778, 0.002734, 0.002692, 0.002651,
     $     0.002611, 0.002574, 0.002538, 0.002505/


C     12025, JENDL evaluated data
      nbins(14) = 214
      data ( eg(14,i), i=1, 214) /7500000.0, 7600000.0, 7800000.0,
     $     8000000.0, 8200000.0, 8400000.0, 8600000.0, 8800000.0,
     $     9000000.0, 9200000.0, 9400000.0, 9600000.0, 9800000.0, 1.0D7,
     $     1.02D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7
     $     , 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7,
     $     1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7
     $     , 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7,
     $     1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7
     $     , 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7,
     $     1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7
     $     , 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7,
     $     2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7
     $     , 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7,
     $     2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7
     $     , 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7,
     $     2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7
     $     , 2.98D7, 3.0D7, 3.02D7, 3.04D7, 3.06D7, 3.08D7, 3.1D7,
     $     3.12D7, 3.14D7, 3.16D7, 3.18D7, 3.2D7, 3.22D7, 3.24D7, 3.26D7
     $     , 3.28D7, 3.3D7, 3.32D7, 3.34D7, 3.36D7, 3.38D7, 3.4D7,
     $     3.42D7, 3.44D7, 3.46D7, 3.48D7, 3.5D7, 3.52D7, 3.54D7, 3.56D7
     $     , 3.58D7, 3.6D7, 3.62D7, 3.64D7, 3.66D7, 3.68D7, 3.7D7,
     $     3.72D7, 3.74D7, 3.76D7, 3.78D7, 3.8D7, 3.82D7, 3.84D7, 3.86D7
     $     , 3.88D7, 3.9D7, 3.92D7, 3.94D7, 3.96D7, 3.98D7, 4.0D7, 4.2D7
     $     , 4.4D7, 4.6D7, 4.8D7, 5.0D7, 5.2D7, 5.4D7, 5.6D7, 5.8D7,
     $     6.0D7, 6.2D7, 6.4D7, 6.6D7, 6.8D7, 7.0D7, 7.2D7, 7.4D7, 7.6D7
     $     , 7.8D7, 8.0D7, 8.2D7, 8.4D7, 8.6D7, 8.8D7, 9.0D7, 9.2D7,
     $     9.4D7, 9.6D7, 9.8D7, 1.0D8, 1.02D8, 1.04D8, 1.06D8, 1.08D8,
     $     1.1D8, 1.12D8, 1.14D8, 1.16D8, 1.18D8, 1.2D8, 1.22D8, 1.24D8,
     $     1.26D8, 1.28D8, 1.3D8, 1.32D8, 1.34D8, 1.36D8, 1.38D8, 1.4D8/

      data ( gdr(14,i), i=1, 214) /1.84D-4, 2.73D-4, 3.53D-4, 4.0D-4,
     $     4.66D-4, 5.53D-4, 6.48D-4, 7.48D-4, 8.5D-4, 9.54D-4, 0.001051
     $     , 0.00113, 0.001194, 0.001254, 0.001322, 0.001396, 0.001474,
     $     0.001554, 0.001638, 0.001726, 0.001818, 0.001915, 0.002017,
     $     0.002124, 0.002236, 0.002353, 0.002477, 0.002607, 0.002743,
     $     0.002886, 0.003036, 0.003194, 0.00336, 0.003535, 0.003727,
     $     0.003939, 0.004165, 0.004405, 0.004663, 0.004937, 0.005229,
     $     0.005541, 0.005875, 0.006232, 0.006603, 0.006995, 0.007433,
     $     0.007897, 0.008389, 0.00891, 0.00946, 0.01004, 0.01066,
     $     0.0113, 0.01199, 0.0127, 0.01345, 0.01424, 0.01507, 0.01594,
     $     0.01683, 0.01771, 0.01862, 0.01967, 0.02095, 0.02229, 0.02352
     $     , 0.02455, 0.02552, 0.02645, 0.02737, 0.02826, 0.02908,
     $     0.02986, 0.03057, 0.03121, 0.03177, 0.03225, 0.03263, 0.03291
     $     , 0.03302, 0.03316, 0.03321, 0.03312, 0.03293, 0.03266,
     $     0.03234, 0.03195, 0.03156, 0.03108, 0.03053, 0.0299, 0.02925,
     $     0.02857, 0.02786, 0.02714, 0.0264, 0.02565, 0.0249, 0.02417,
     $     0.02343, 0.02273, 0.02202, 0.02133, 0.02066, 0.02001, 0.01938
     $     , 0.01877, 0.01817, 0.01759, 0.01705, 0.01652, 0.01602,
     $     0.01552, 0.01506, 0.01458, 0.01413, 0.0137, 0.01329, 0.0129,
     $     0.01252, 0.01216, 0.01182, 0.0115, 0.01119, 0.01089, 0.0106,
     $     0.01033, 0.01008, 0.009828, 0.009591, 0.009363, 0.009144,
     $     0.008936, 0.008739, 0.008549, 0.008372, 0.008199, 0.008036,
     $     0.00788, 0.007732, 0.00759, 0.00745, 0.007315, 0.007187,
     $     0.007063, 0.006944, 0.006829, 0.006721, 0.006616, 0.006494,
     $     0.006397, 0.006301, 0.006208, 0.006119, 0.006041, 0.005968,
     $     0.005897, 0.00583, 0.005765, 0.005702, 0.005638, 0.005574,
     $     0.005511, 0.005057, 0.004692, 0.0044, 0.004153, 0.003934,
     $     0.003739, 0.003542, 0.003414, 0.003302, 0.003215, 0.003146,
     $     0.003087, 0.003034, 0.002994, 0.002957, 0.002914, 0.002884,
     $     0.002845, 0.002803, 0.002761, 0.002722, 0.002677, 0.00263,
     $     0.002584, 0.002535, 0.002483, 0.002433, 0.002386, 0.002335,
     $     0.002282, 0.002236, 0.002185, 0.002139, 0.002094, 0.002044,
     $     0.002001, 0.00196, 0.001917, 0.001874, 0.001832, 0.001796,
     $     0.001761, 0.001724, 0.00169, 0.00166, 0.001627, 0.001597,
     $     0.001569, 0.001542, 0.001518/


C     12026, JENDL evaluated data
      nbins(15) = 293
      data ( eg(15,i), i=1, 293) /1.12D7, 1.13D7, 1.14D7, 1.15D7,
     $     1.16D7, 1.17D7, 1.18D7, 1.19D7, 1.2D7, 1.21D7, 1.22D7, 1.23D7
     $     , 1.24D7, 1.25D7, 1.26D7, 1.27D7, 1.28D7, 1.29D7, 1.3D7,
     $     1.31D7, 1.32D7, 1.33D7, 1.34D7, 1.35D7, 1.36D7, 1.37D7,
     $     1.38D7, 1.39D7, 1.4D7, 1.41D7, 1.414D7, 1.42D7, 1.43D7,
     $     1.44D7, 1.45D7, 1.46D7, 1.47D7, 1.48D7, 1.49D7, 1.5D7, 1.51D7
     $     , 1.52D7, 1.53D7, 1.54D7, 1.55D7, 1.56D7, 1.57D7, 1.58D7,
     $     1.59D7, 1.6D7, 1.61D7, 1.62D7, 1.63D7, 1.64D7, 1.65D7, 1.66D7
     $     , 1.67D7, 1.68D7, 1.69D7, 1.7D7, 1.71D7, 1.72D7, 1.73D7,
     $     1.74D7, 1.75D7, 1.76D7, 1.77D7, 1.78D7, 1.79D7, 1.8D7, 1.81D7
     $     , 1.82D7, 1.83D7, 1.84D7, 1.842D7, 1.85D7, 1.86D7, 1.87D7,
     $     1.88D7, 1.89D7, 1.9D7, 1.91D7, 1.92D7, 1.93D7, 1.94D7, 1.95D7
     $     , 1.96D7, 1.97D7, 1.98D7, 1.99D7, 2.0D7, 2.01D7, 2.02D7,
     $     2.03D7, 2.04D7, 2.05D7, 2.06D7, 2.07D7, 2.08D7, 2.09D7,
     $     2.093D7, 2.1D7, 2.11D7, 2.12D7, 2.13D7, 2.14D7, 2.15D7,
     $     2.16D7, 2.17D7, 2.18D7, 2.19D7, 2.2D7, 2.21D7, 2.22D7, 2.23D7
     $     , 2.24D7, 2.25D7, 2.26D7, 2.27D7, 2.28D7, 2.29D7, 2.3D7,
     $     2.31D7, 2.32D7, 2.33D7, 2.34D7, 2.35D7, 2.36D7, 2.37D7,
     $     2.38D7, 2.39D7, 2.4D7, 2.41D7, 2.42D7, 2.43D7, 2.44D7, 2.45D7
     $     , 2.46D7, 2.47D7, 2.48D7, 2.49D7, 2.5D7, 2.51D7, 2.52D7,
     $     2.53D7, 2.54D7, 2.55D7, 2.56D7, 2.57D7, 2.58D7, 2.59D7, 2.6D7
     $     , 2.61D7, 2.62D7, 2.63D7, 2.64D7, 2.65D7, 2.66D7, 2.67D7,
     $     2.68D7, 2.69D7, 2.7D7, 2.71D7, 2.72D7, 2.73D7, 2.74D7, 2.75D7
     $     , 2.76D7, 2.77D7, 2.78D7, 2.79D7, 2.8D7, 2.81D7, 2.82D7,
     $     2.83D7, 2.84D7, 2.85D7, 2.86D7, 2.87D7, 2.88D7, 2.89D7, 2.9D7
     $     , 2.91D7, 2.92D7, 2.93D7, 2.94D7, 2.95D7, 2.96D7, 2.97D7,
     $     2.98D7, 2.99D7, 3.0D7, 3.02D7, 3.04D7, 3.06D7, 3.08D7, 3.1D7,
     $     3.12D7, 3.14D7, 3.16D7, 3.18D7, 3.2D7, 3.22D7, 3.24D7, 3.26D7
     $     , 3.28D7, 3.3D7, 3.32D7, 3.34D7, 3.36D7, 3.38D7, 3.4D7,
     $     3.42D7, 3.44D7, 3.46D7, 3.48D7, 3.495D7, 3.5D7, 3.52D7,
     $     3.54D7, 3.56D7, 3.58D7, 3.6D7, 3.62D7, 3.64D7, 3.66D7, 3.68D7
     $     , 3.7D7, 3.72D7, 3.74D7, 3.76D7, 3.78D7, 3.8D7, 3.82D7,
     $     3.84D7, 3.86D7, 3.88D7, 3.9D7, 3.92D7, 3.94D7, 3.96D7, 3.98D7
     $     , 4.0D7, 4.2D7, 4.4D7, 4.6D7, 4.8D7, 5.0D7, 5.2D7, 5.4D7,
     $     5.6D7, 5.8D7, 6.0D7, 6.2D7, 6.4D7, 6.6D7, 6.8D7, 7.0D7, 7.2D7
     $     , 7.4D7, 7.6D7, 7.8D7, 8.0D7, 8.2D7, 8.4D7, 8.6D7, 8.8D7,
     $     9.0D7, 9.2D7, 9.4D7, 9.6D7, 9.8D7, 1.0D8, 1.02D8, 1.04D8,
     $     1.06D8, 1.08D8, 1.1D8, 1.12D8, 1.14D8, 1.16D8, 1.18D8, 1.2D8,
     $     1.22D8, 1.24D8, 1.26D8, 1.28D8, 1.3D8, 1.32D8, 1.34D8, 1.36D8
     $     , 1.38D8, 1.4D8/

      data ( gdr(15,i), i=1, 293) /0.002272, 0.001307, 0.003674,
     $     0.002555, 0.002833, 0.007763, 0.008711, 0.003285, 0.002281,
     $     0.002399, 0.00321, 0.004632, 0.005201, 0.00413, 0.003405,
     $     0.003892, 0.006801, 0.008535, 0.004485, 0.002989, 0.003112,
     $     0.005564, 0.004765, 0.003384, 0.004503, 0.008121, 0.006766,
     $     0.008593, 0.004834, 0.003133, 0.003148, 0.003171, 0.003849,
     $     0.005247, 0.007335, 0.008464, 0.00731, 0.006101, 0.006176,
     $     0.008134, 0.01202, 0.01165, 0.007557, 0.005291, 0.004493,
     $     0.004956, 0.01254, 0.006179, 0.005702, 0.00736, 0.01137,
     $     0.01739, 0.01622, 0.01225, 0.01147, 0.01287, 0.01497, 0.01623
     $     , 0.017, 0.01844, 0.01948, 0.01873, 0.01773, 0.01853, 0.02183
     $     , 0.0261, 0.02665, 0.0226, 0.01927, 0.0256, 0.01958, 0.0233,
     $     0.02116, 0.02329, 0.0224, 0.01894, 0.01532, 0.01596, 0.01558,
     $     0.01351, 0.0114, 0.01028, 0.01057, 0.01234, 0.0119, 0.01152,
     $     0.01297, 0.01546, 0.01708, 0.0163, 0.01478, 0.0146, 0.01662,
     $     0.01848, 0.01663, 0.01552, 0.01648, 0.01879, 0.02117, 0.02183
     $     , 0.02145, 0.02064, 0.02017, 0.02217, 0.02114, 0.02771,
     $     0.03193, 0.0285, 0.02995, 0.03046, 0.02802, 0.03022, 0.03541,
     $     0.03173, 0.03543, 0.02717, 0.02388, 0.02299, 0.02278, 0.02311
     $     , 0.02387, 0.02454, 0.02558, 0.02658, 0.02741, 0.02798,
     $     0.02828, 0.0283, 0.02808, 0.02768, 0.02719, 0.02668, 0.02622,
     $     0.02584, 0.0256, 0.02549, 0.02554, 0.02575, 0.02606, 0.02648,
     $     0.02696, 0.02744, 0.02788, 0.02821, 0.02832, 0.02827, 0.02802
     $     , 0.02754, 0.02685, 0.026, 0.02501, 0.02394, 0.02285, 0.02179
     $     , 0.02079, 0.01992, 0.01921, 0.01877, 0.01874, 0.01949,
     $     0.0218, 0.02689, 0.03148, 0.02695, 0.02152, 0.019, 0.01829,
     $     0.01863, 0.01978, 0.02159, 0.02388, 0.02621, 0.02783, 0.02796
     $     , 0.02635, 0.02357, 0.02041, 0.01747, 0.01496, 0.01293,
     $     0.0113, 0.009992, 0.008937, 0.00808, 0.007373, 0.006783,
     $     0.006285, 0.00586, 0.005493, 0.005174, 0.004894, 0.004646,
     $     0.004224, 0.003884, 0.003603, 0.003366, 0.003164, 0.00299,
     $     0.002837, 0.002703, 0.002585, 0.002479, 0.002385, 0.002233,
     $     0.002159, 0.002091, 0.00203, 0.001975, 0.001923, 0.001877,
     $     0.001833, 0.001794, 0.001758, 0.001724, 0.001694, 0.001665,
     $     0.001643, 0.001637, 0.001613, 0.00159, 0.001569, 0.001548,
     $     0.00153, 0.001513, 0.001497, 0.001482, 0.001467, 0.001454,
     $     0.001443, 0.001431, 0.001421, 0.001411, 0.001402, 0.001393,
     $     0.001385, 0.001377, 0.00137, 0.001364, 0.001357, 0.001353,
     $     0.001347, 0.001342, 0.001338, 0.00131, 0.001301, 0.001305,
     $     0.001316, 0.00133, 0.001344, 0.001357, 0.001367, 0.001376,
     $     0.001383, 0.00139, 0.001398, 0.001408, 0.001419, 0.001429,
     $     0.001439, 0.001449, 0.001457, 0.001464, 0.001469, 0.001476,
     $     0.001478, 0.001475, 0.001472, 0.001466, 0.001456, 0.001444,
     $     0.00143, 0.001414, 0.001396, 0.001378, 0.001357, 0.001336,
     $     0.001314, 0.001291, 0.001268, 0.001246, 0.001223, 0.001199,
     $     0.001176, 0.001153, 0.001132, 0.001109, 0.001089, 0.001068,
     $     0.001049, 0.00103, 0.001012, 9.959D-4, 9.798D-4/


C     13027, JENDL evaluated data
      nbins(16) = 321
      data ( eg(16,i), i=1, 321) /8271040.0, 8700000.0, 8800000.0,
     $     8900000.0, 9000000.0, 9100000.0, 9200000.0, 9300000.0,
     $     9400000.0, 9500000.0, 9600000.0, 9700000.0, 9800000.0,
     $     9900000.0, 1.0D7, 1.01D7, 1.02D7, 1.03D7, 1.04D7, 1.05D7,
     $     1.06D7, 1.07D7, 1.08D7, 1.09D7, 1.1D7, 1.11D7, 1.12D7, 1.13D7
     $     , 1.14D7, 1.15D7, 1.16D7, 1.17D7, 1.18D7, 1.19D7, 1.2D7,
     $     1.21D7, 1.22D7, 1.23D7, 1.24D7, 1.25D7, 1.26D7, 1.27D7,
     $     1.28D7, 1.29D7, 1.3D7, 1.31D7, 1.32D7, 1.33D7, 1.34D7, 1.35D7
     $     , 1.36D7, 1.37D7, 1.38D7, 1.39D7, 1.4D7, 1.41D7, 1.42D7,
     $     1.43D7, 1.44D7, 1.45D7, 1.46D7, 1.47D7, 1.48D7, 1.49D7, 1.5D7
     $     , 1.51D7, 1.52D7, 1.53D7, 1.54D7, 1.55D7, 1.56D7, 1.57D7,
     $     1.58D7, 1.59D7, 1.6D7, 1.61D7, 1.62D7, 1.63D7, 1.64D7, 1.65D7
     $     , 1.66D7, 1.67D7, 1.68D7, 1.69D7, 1.7D7, 1.71D7, 1.72D7,
     $     1.73D7, 1.74D7, 1.75D7, 1.76D7, 1.77D7, 1.78D7, 1.79D7, 1.8D7
     $     , 1.81D7, 1.82D7, 1.83D7, 1.84D7, 1.85D7, 1.86D7, 1.87D7,
     $     1.88D7, 1.89D7, 1.9D7, 1.91D7, 1.92D7, 1.93D7, 1.94D7, 1.95D7
     $     , 1.96D7, 1.97D7, 1.98D7, 1.99D7, 2.0D7, 2.01D7, 2.02D7,
     $     2.03D7, 2.04D7, 2.05D7, 2.06D7, 2.07D7, 2.08D7, 2.09D7, 2.1D7
     $     , 2.11D7, 2.12D7, 2.13D7, 2.14D7, 2.15D7, 2.16D7, 2.17D7,
     $     2.18D7, 2.19D7, 2.2D7, 2.21D7, 2.22D7, 2.23D7, 2.24D7, 2.25D7
     $     , 2.26D7, 2.27D7, 2.28D7, 2.29D7, 2.3D7, 2.31D7, 2.32D7,
     $     2.33D7, 2.34D7, 2.35D7, 2.36D7, 2.37D7, 2.38D7, 2.39D7, 2.4D7
     $     , 2.41D7, 2.42D7, 2.43D7, 2.44D7, 2.45D7, 2.46D7, 2.47D7,
     $     2.48D7, 2.49D7, 2.5D7, 2.51D7, 2.52D7, 2.53D7, 2.54D7, 2.55D7
     $     , 2.56D7, 2.57D7, 2.58D7, 2.59D7, 2.6D7, 2.61D7, 2.62D7,
     $     2.63D7, 2.64D7, 2.65D7, 2.66D7, 2.67D7, 2.68D7, 2.69D7, 2.7D7
     $     , 2.71D7, 2.72D7, 2.73D7, 2.74D7, 2.75D7, 2.76D7, 2.77D7,
     $     2.78D7, 2.79D7, 2.8D7, 2.81D7, 2.82D7, 2.83D7, 2.84D7, 2.85D7
     $     , 2.86D7, 2.87D7, 2.88D7, 2.89D7, 2.9D7, 2.91D7, 2.92D7,
     $     2.93D7, 2.94D7, 2.95D7, 2.96D7, 2.97D7, 2.98D7, 2.99D7, 3.0D7
     $     , 3.02D7, 3.04D7, 3.06D7, 3.08D7, 3.1D7, 3.12D7, 3.14D7,
     $     3.16D7, 3.18D7, 3.2D7, 3.22D7, 3.24D7, 3.26D7, 3.28D7, 3.3D7,
     $     3.32D7, 3.34D7, 3.36D7, 3.38D7, 3.4D7, 3.42D7, 3.44D7, 3.46D7
     $     , 3.48D7, 3.5D7, 3.52D7, 3.54D7, 3.56D7, 3.58D7, 3.6D7,
     $     3.62D7, 3.64D7, 3.66D7, 3.68D7, 3.7D7, 3.72D7, 3.74D7, 3.76D7
     $     , 3.78D7, 3.8D7, 3.82D7, 3.84D7, 3.86D7, 3.88D7, 3.9D7,
     $     3.92D7, 3.94D7, 3.96D7, 3.98D7, 4.0D7, 4.05D7, 4.1D7, 4.2D7,
     $     4.3D7, 4.4D7, 4.5D7, 4.6D7, 4.7D7, 4.8D7, 4.9D7, 5.0D7, 5.2D7
     $     , 5.4D7, 5.6D7, 5.8D7, 6.0D7, 6.2D7, 6.4D7, 6.6D7, 6.8D7,
     $     7.0D7, 7.2D7, 7.4D7, 7.6D7, 7.8D7, 8.0D7, 8.2D7, 8.4D7, 8.6D7
     $     , 8.8D7, 9.0D7, 9.2D7, 9.4D7, 9.6D7, 9.8D7, 1.0D8, 1.02D8,
     $     1.04D8, 1.06D8, 1.08D8, 1.1D8, 1.12D8, 1.14D8, 1.16D8, 1.18D8
     $     , 1.2D8, 1.22D8, 1.24D8, 1.26D8, 1.28D8, 1.3D8, 1.32D8,
     $     1.34D8, 1.36D8, 1.38D8, 1.4D8/

      data ( gdr(16,i), i=1, 321) /0.0, 6.46D-5, 9.83D-5, 1.39D-4,
     $     1.867D-4, 2.414D-4, 3.031D-4, 3.719D-4, 4.477D-4, 5.305D-4,
     $     6.204D-4, 7.172D-4, 8.211D-4, 9.321D-4, 9.69D-4, 0.001166,
     $     0.001452, 0.001894, 0.002609, 0.00371, 0.004892, 0.004996,
     $     0.003958, 0.002936, 0.002287, 0.001914, 0.001702, 0.001585,
     $     0.001526, 0.001506, 0.001515, 0.001548, 0.001607, 0.001697,
     $     0.001835, 0.002057, 0.002451, 0.003243, 0.004911, 0.006702,
     $     0.005308, 0.003672, 0.002949, 0.002671, 0.002599, 0.002646,
     $     0.00279, 0.003045, 0.003467, 0.004166, 0.005296, 0.006719,
     $     0.007245, 0.006326, 0.005291, 0.004733, 0.00459, 0.004788,
     $     0.005352, 0.006358, 0.00751, 0.007702, 0.006953, 0.006363,
     $     0.006224, 0.006465, 0.007038, 0.007954, 0.009227, 0.01072,
     $     0.01198, 0.01251, 0.01252, 0.01261, 0.01312, 0.01403, 0.01509
     $     , 0.01592, 0.01624, 0.01609, 0.01589, 0.01613, 0.01726,
     $     0.01918, 0.02028, 0.01931, 0.01803, 0.01782, 0.01873, 0.02037
     $     , 0.02207, 0.02315, 0.02367, 0.02439, 0.02598, 0.02843,
     $     0.03073, 0.03148, 0.03059, 0.02939, 0.02893, 0.02934, 0.03022
     $     , 0.03097, 0.03121, 0.03102, 0.0308, 0.03094, 0.0316, 0.03277
     $     , 0.0343, 0.03595, 0.03745, 0.0385, 0.03895, 0.03879, 0.0382,
     $     0.03744, 0.03679, 0.03649, 0.03669, 0.03745, 0.03873, 0.04033
     $     , 0.04191, 0.04303, 0.04335, 0.04278, 0.04157, 0.04015,
     $     0.03893, 0.03815, 0.03785, 0.03789, 0.03797, 0.03773, 0.03693
     $     , 0.0356, 0.03398, 0.03237, 0.03099, 0.02995, 0.02924, 0.0288
     $     , 0.02849, 0.0282, 0.02781, 0.02727, 0.02661, 0.0259, 0.02522
     $     , 0.02466, 0.02425, 0.024, 0.02384, 0.02368, 0.02344, 0.02306
     $     , 0.02261, 0.02215, 0.02179, 0.02157, 0.02151, 0.02159,
     $     0.02175, 0.02191, 0.022, 0.02195, 0.02174, 0.02139, 0.02097,
     $     0.02054, 0.02014, 0.01979, 0.01951, 0.0193, 0.01915, 0.01905,
     $     0.019, 0.01897, 0.01897, 0.01899, 0.01901, 0.01904, 0.01906,
     $     0.01908, 0.01909, 0.01909, 0.01907, 0.01904, 0.01899, 0.01893
     $     , 0.01885, 0.01876, 0.01865, 0.01853, 0.01841, 0.01827,
     $     0.01812, 0.01798, 0.01782, 0.01767, 0.01751, 0.01735, 0.01719
     $     , 0.01704, 0.01688, 0.01673, 0.01659, 0.01644, 0.01631,
     $     0.01617, 0.01604, 0.01592, 0.0158, 0.01557, 0.01535, 0.01515,
     $     0.01497, 0.01479, 0.01461, 0.01445, 0.01429, 0.01412, 0.01396
     $     , 0.0138, 0.01364, 0.01347, 0.0133, 0.01313, 0.01295, 0.01277
     $     , 0.01259, 0.0124, 0.01222, 0.01203, 0.01183, 0.01164,
     $     0.01145, 0.01125, 0.01106, 0.01087, 0.01068, 0.01049, 0.0103,
     $     0.01011, 0.009932, 0.009753, 0.009578, 0.009405, 0.009236,
     $     0.009071, 0.00891, 0.008753, 0.0086, 0.00845, 0.008305,
     $     0.008164, 0.008027, 0.007894, 0.007765, 0.00764, 0.007518,
     $     0.0074, 0.007286, 0.007017, 0.006768, 0.006327, 0.005954,
     $     0.005636, 0.005366, 0.005135, 0.004937, 0.004766, 0.004619,
     $     0.00449, 0.004278, 0.004114, 0.003982, 0.003875, 0.003786,
     $     0.00371, 0.003644, 0.003585, 0.003531, 0.003481, 0.003434,
     $     0.003388, 0.003343, 0.0033, 0.003257, 0.003214, 0.003171,
     $     0.003128, 0.003085, 0.003041, 0.002997, 0.002953, 0.002909,
     $     0.002865, 0.002821, 0.002776, 0.002732, 0.002688, 0.002644,
     $     0.0026, 0.002557, 0.002514, 0.002472, 0.00243, 0.002389,
     $     0.00235, 0.002311, 0.002273, 0.002236, 0.002201, 0.002167,
     $     0.002134, 0.002103, 0.002074, 0.002046/


C     14028, JENDL evaluated data
      nbins(17) = 275
      data ( eg(17,i), i=1, 275) /1.01D7, 1.02D7, 1.03D7, 1.04D7,
     $     1.05D7, 1.06D7, 1.07D7, 1.08D7, 1.09D7, 1.1D7, 1.11D7, 1.12D7
     $     , 1.13D7, 1.14D7, 1.15D7, 1.16D7, 1.17D7, 1.18D7, 1.19D7,
     $     1.2D7, 1.21D7, 1.22D7, 1.23D7, 1.24D7, 1.25D7, 1.26D7, 1.27D7
     $     , 1.28D7, 1.29D7, 1.3D7, 1.31D7, 1.32D7, 1.33D7, 1.34D7,
     $     1.35D7, 1.36D7, 1.37D7, 1.38D7, 1.39D7, 1.4D7, 1.41D7, 1.42D7
     $     , 1.43D7, 1.44D7, 1.45D7, 1.46D7, 1.47D7, 1.48D7, 1.49D7,
     $     1.5D7, 1.51D7, 1.52D7, 1.53D7, 1.54D7, 1.55D7, 1.56D7, 1.57D7
     $     , 1.58D7, 1.59D7, 1.6D7, 1.61D7, 1.62D7, 1.63D7, 1.64D7,
     $     1.65D7, 1.66D7, 1.67D7, 1.68D7, 1.69D7, 1.7D7, 1.71D7, 1.72D7
     $     , 1.73D7, 1.74D7, 1.75D7, 1.76D7, 1.77D7, 1.78D7, 1.79D7,
     $     1.8D7, 1.81D7, 1.82D7, 1.83D7, 1.84D7, 1.85D7, 1.86D7, 1.87D7
     $     , 1.88D7, 1.89D7, 1.9D7, 1.91D7, 1.92D7, 1.93D7, 1.94D7,
     $     1.95D7, 1.96D7, 1.97D7, 1.98D7, 1.99D7, 2.0D7, 2.01D7, 2.02D7
     $     , 2.03D7, 2.04D7, 2.05D7, 2.06D7, 2.07D7, 2.08D7, 2.09D7,
     $     2.1D7, 2.11D7, 2.12D7, 2.13D7, 2.14D7, 2.15D7, 2.16D7, 2.17D7
     $     , 2.18D7, 2.19D7, 2.2D7, 2.21D7, 2.22D7, 2.23D7, 2.24D7,
     $     2.25D7, 2.26D7, 2.27D7, 2.28D7, 2.29D7, 2.3D7, 2.31D7, 2.32D7
     $     , 2.33D7, 2.34D7, 2.35D7, 2.36D7, 2.37D7, 2.38D7, 2.39D7,
     $     2.4D7, 2.41D7, 2.42D7, 2.43D7, 2.44D7, 2.45D7, 2.46D7, 2.47D7
     $     , 2.48D7, 2.49D7, 2.5D7, 2.52D7, 2.54D7, 2.56D7, 2.58D7,
     $     2.6D7, 2.62D7, 2.64D7, 2.66D7, 2.68D7, 2.7D7, 2.72D7, 2.74D7,
     $     2.76D7, 2.78D7, 2.8D7, 2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7,
     $     2.92D7, 2.94D7, 2.96D7, 2.98D7, 3.0D7, 3.02D7, 3.04D7, 3.06D7
     $     , 3.08D7, 3.1D7, 3.12D7, 3.14D7, 3.16D7, 3.18D7, 3.2D7,
     $     3.22D7, 3.24D7, 3.26D7, 3.28D7, 3.3D7, 3.32D7, 3.34D7, 3.36D7
     $     , 3.38D7, 3.4D7, 3.42D7, 3.44D7, 3.46D7, 3.48D7, 3.5D7,
     $     3.52D7, 3.54D7, 3.56D7, 3.58D7, 3.6D7, 3.62D7, 3.64D7, 3.66D7
     $     , 3.68D7, 3.7D7, 3.72D7, 3.74D7, 3.76D7, 3.78D7, 3.8D7,
     $     3.82D7, 3.84D7, 3.86D7, 3.88D7, 3.9D7, 3.92D7, 3.94D7, 3.96D7
     $     , 3.98D7, 4.0D7, 4.2D7, 4.4D7, 4.6D7, 4.8D7, 5.0D7, 5.2D7,
     $     5.4D7, 5.6D7, 5.8D7, 6.0D7, 6.2D7, 6.4D7, 6.6D7, 6.8D7, 7.0D7
     $     , 7.2D7, 7.4D7, 7.6D7, 7.8D7, 8.0D7, 8.2D7, 8.4D7, 8.6D7,
     $     8.8D7, 9.0D7, 9.2D7, 9.4D7, 9.6D7, 9.8D7, 1.0D8, 1.02D8,
     $     1.04D8, 1.06D8, 1.08D8, 1.1D8, 1.12D8, 1.14D8, 1.16D8, 1.18D8
     $     , 1.2D8, 1.22D8, 1.24D8, 1.26D8, 1.28D8, 1.3D8, 1.32D8,
     $     1.34D8, 1.36D8, 1.38D8, 1.4D8/

      data ( gdr(17,i), i=1, 275) /0.003639, 0.00392, 0.004152, 0.004311
     $     , 0.004385, 0.004371, 0.004277, 0.004123, 0.003931, 0.003722,
     $     0.003512, 0.003314, 0.003135, 0.002978, 0.002846, 0.002739,
     $     0.002656, 0.002599, 0.002569, 0.002569, 0.00261, 0.00271,
     $     0.002911, 0.003306, 0.004021, 0.004761, 0.00447, 0.003736,
     $     0.0033, 0.00311, 0.003054, 0.003072, 0.003138, 0.003239,
     $     0.00337, 0.003528, 0.003712, 0.003923, 0.004161, 0.004426,
     $     0.004718, 0.005034, 0.005371, 0.005721, 0.006074, 0.006415,
     $     0.006727, 0.006992, 0.007193, 0.007319, 0.007364, 0.007333,
     $     0.007237, 0.007093, 0.006918, 0.00673, 0.006545, 0.006374,
     $     0.006225, 0.006106, 0.00602, 0.00597, 0.005959, 0.00599,
     $     0.006067, 0.006193, 0.006375, 0.006623, 0.00695, 0.007374,
     $     0.007923, 0.008641, 0.009591, 0.01087, 0.01267, 0.01526,
     $     0.01913, 0.02498, 0.03323, 0.04172, 0.0442, 0.03956, 0.0342,
     $     0.03148, 0.03148, 0.03354, 0.03686, 0.04032, 0.04252, 0.04258
     $     , 0.04094, 0.03897, 0.038, 0.039, 0.04264, 0.04928, 0.05792,
     $     0.06437, 0.06327, 0.05551, 0.04674, 0.04041, 0.03722, 0.03713
     $     , 0.04007, 0.04538, 0.05025, 0.05172, 0.05184, 0.05345,
     $     0.05531, 0.05452, 0.05116, 0.04769, 0.04512, 0.04258, 0.0392,
     $     0.03537, 0.03199, 0.0296, 0.0283, 0.02792, 0.02814, 0.02855,
     $     0.02866, 0.02824, 0.02742, 0.0265, 0.02576, 0.02527, 0.02502,
     $     0.02492, 0.02486, 0.02475, 0.02455, 0.02423, 0.02381, 0.02334
     $     , 0.02287, 0.02243, 0.02206, 0.02176, 0.02155, 0.02141,
     $     0.02134, 0.02133, 0.02134, 0.02137, 0.02141, 0.02143, 0.02138
     $     , 0.02119, 0.02086, 0.02043, 0.01996, 0.0195, 0.01908,
     $     0.01873, 0.01845, 0.01822, 0.01804, 0.01789, 0.01776, 0.01763
     $     , 0.01748, 0.01732, 0.01712, 0.01689, 0.01662, 0.01632,
     $     0.01598, 0.01561, 0.01521, 0.01479, 0.01435, 0.0139, 0.01344,
     $     0.01298, 0.01253, 0.01208, 0.01164, 0.01122, 0.01081, 0.01041
     $     , 0.01003, 0.009665, 0.009316, 0.008983, 0.008667, 0.008367,
     $     0.008082, 0.007812, 0.007556, 0.007314, 0.007085, 0.006868,
     $     0.006663, 0.006469, 0.006285, 0.006112, 0.005948, 0.005792,
     $     0.005645, 0.005505, 0.005373, 0.005247, 0.005128, 0.005016,
     $     0.004908, 0.004807, 0.00471, 0.004618, 0.00453, 0.004447,
     $     0.004368, 0.004292, 0.00422, 0.004151, 0.004086, 0.004023,
     $     0.003964, 0.003907, 0.003852, 0.0038, 0.00375, 0.003354,
     $     0.003089, 0.002906, 0.002776, 0.002682, 0.002611, 0.002557,
     $     0.002514, 0.002479, 0.00245, 0.002424, 0.002401, 0.002379,
     $     0.002358, 0.002337, 0.002316, 0.002295, 0.002274, 0.002251,
     $     0.002229, 0.002205, 0.002181, 0.002155, 0.00213, 0.002103,
     $     0.002076, 0.002049, 0.002021, 0.001992, 0.001964, 0.001935,
     $     0.001906, 0.001876, 0.001847, 0.001818, 0.001789, 0.00176,
     $     0.001732, 0.001703, 0.001676, 0.001649, 0.001622, 0.001596,
     $     0.001571, 0.001547, 0.001524, 0.001501, 0.00148, 0.00146,
     $     0.001441/


C     14029, JENDL evaluated data
      nbins(18) = 412
      data ( eg(18,i), i=1, 412) /8473560.0, 8600000.0, 8700000.0,
     $     8800000.0, 8900000.0, 9000000.0, 9100000.0, 9200000.0,
     $     9300000.0, 9400000.0, 9500000.0, 9600000.0, 9700000.0,
     $     9800000.0, 9900000.0, 1.0D7, 1.01D7, 1.02D7, 1.03D7, 1.04D7,
     $     1.05D7, 1.06D7, 1.07D7, 1.08D7, 1.09D7, 1.1D7, 1.11D7,
     $     1.11272D7, 1.113D7, 1.12D7, 1.13D7, 1.14D7, 1.15D7, 1.16D7,
     $     1.17D7, 1.18D7, 1.19D7, 1.2D7, 1.21D7, 1.22D7, 1.23D7,
     $     1.233D7, 1.23334D7, 1.24D7, 1.25D7, 1.26D7, 1.27D7, 1.28D7,
     $     1.29D7, 1.3D7, 1.31D7, 1.32D7, 1.33D7, 1.34D7, 1.35D7, 1.36D7
     $     , 1.37D7, 1.38D7, 1.39D7, 1.4D7, 1.41D7, 1.42D7, 1.43D7,
     $     1.44D7, 1.45D7, 1.46D7, 1.47D7, 1.48D7, 1.49D7, 1.5D7, 1.51D7
     $     , 1.52D7, 1.53D7, 1.54D7, 1.55D7, 1.56D7, 1.57D7, 1.58D7,
     $     1.59D7, 1.6D7, 1.61D7, 1.62D7, 1.63D7, 1.64D7, 1.65D7, 1.66D7
     $     , 1.67D7, 1.68D7, 1.69D7, 1.7D7, 1.71D7, 1.72D7, 1.73D7,
     $     1.74D7, 1.75D7, 1.76D7, 1.77D7, 1.78D7, 1.783D7, 1.78339D7,
     $     1.79D7, 1.8D7, 1.81D7, 1.82D7, 1.83D7, 1.84D7, 1.85D7, 1.86D7
     $     , 1.87D7, 1.88D7, 1.89D7, 1.9D7, 1.91D7, 1.92D7, 1.93D7,
     $     1.93479D7, 1.94D7, 1.95D7, 1.96D7, 1.97D7, 1.98D7, 1.99D7,
     $     2.0D7, 2.01D7, 2.02D7, 2.0283D7, 2.03D7, 2.04D7, 2.05D7,
     $     2.06D7, 2.06118D7, 2.07D7, 2.08D7, 2.09D7, 2.1D7, 2.1086D7,
     $     2.11D7, 2.1105D7, 2.12D7, 2.13D7, 2.14D7, 2.15D7, 2.16D7,
     $     2.17D7, 2.18D7, 2.18865D7, 2.19D7, 2.2D7, 2.21D7, 2.22D7,
     $     2.23D7, 2.24D7, 2.25D7, 2.26D7, 2.27D7, 2.28D7, 2.29D7, 2.3D7
     $     , 2.31D7, 2.32D7, 2.33D7, 2.34D7, 2.35D7, 2.36D7, 2.37D7,
     $     2.38D7, 2.39D7, 2.4D7, 2.41D7, 2.42D7, 2.43D7, 2.44D7, 2.45D7
     $     , 2.46D7, 2.46345D7, 2.47D7, 2.48D7, 2.48599D7, 2.49D7, 2.5D7
     $     , 2.51D7, 2.51563D7, 2.52D7, 2.53D7, 2.54D7, 2.55D7, 2.56D7,
     $     2.565D7, 2.56532D7, 2.57D7, 2.58D7, 2.59D7, 2.6D7, 2.61D7,
     $     2.62D7, 2.62007D7, 2.63D7, 2.64D7, 2.64519D7, 2.65D7, 2.66D7,
     $     2.67D7, 2.68D7, 2.69D7, 2.7D7, 2.71D7, 2.72D7, 2.73D7, 2.74D7
     $     , 2.75D7, 2.76D7, 2.77D7, 2.77703D7, 2.78D7, 2.79D7, 2.8D7,
     $     2.81D7, 2.82D7, 2.83D7, 2.84D7, 2.85D7, 2.86D7, 2.87D7,
     $     2.88D7, 2.89D7, 2.9D7, 2.91D7, 2.92D7, 2.93D7, 2.94D7, 2.95D7
     $     , 2.96D7, 2.96402D7, 2.97D7, 2.98D7, 2.99D7, 3.0D7, 3.01D7,
     $     3.02D7, 3.03D7, 3.04D7, 3.05D7, 3.06D7, 3.07D7, 3.08D7,
     $     3.09D7, 3.1D7, 3.11D7, 3.12D7, 3.14D7, 3.14364D7, 3.16D7,
     $     3.16206D7, 3.17635D7, 3.18D7, 3.2D7, 3.22D7, 3.24D7,
     $     3.24584D7, 3.26D7, 3.28D7, 3.29239D7, 3.3D7, 3.32D7, 3.34D7,
     $     3.36D7, 3.38D7, 3.4D7, 3.40165D7, 3.42D7, 3.42933D7, 3.44D7,
     $     3.46D7, 3.48D7, 3.5D7, 3.51158D7, 3.52D7, 3.54D7, 3.56D7,
     $     3.58D7, 3.6D7, 3.62D7, 3.62527D7, 3.64D7, 3.66D7, 3.68D7,
     $     3.68595D7, 3.7D7, 3.72D7, 3.72156D7, 3.74D7, 3.76D7, 3.78D7,
     $     3.78577D7, 3.8D7, 3.82D7, 3.84D7, 3.86D7, 3.88D7, 3.89071D7,
     $     3.89644D7, 3.897D7, 3.9D7, 3.92D7, 3.94D7, 3.96D7, 3.97746D7,
     $     3.98D7, 4.0D7, 4.0657D7, 4.09272D7, 4.2D7, 4.4D7, 4.47684D7,
     $     4.6D7, 4.71039D7, 4.73758D7, 4.76085D7, 4.8D7, 4.89922D7,
     $     4.91375D7, 5.0D7, 5.2D7, 5.25567D7, 5.36938D7, 5.39435D7,
     $     5.4D7, 5.46056D7, 5.6D7, 5.8D7, 5.80056D7, 5.8674D7,
     $     5.86999D7, 5.99727D7, 6.0D7, 6.2D7, 6.23358D7, 6.34566D7,
     $     6.4D7, 6.44047D7, 6.50508D7, 6.58844D7, 6.6D7, 6.78739D7,
     $     6.8D7, 7.0D7, 7.0526D7, 7.09674D7, 7.2D7, 7.30064D7, 7.4D7,
     $     7.49795D7, 7.6D7, 7.74545D7, 7.8D7, 7.95459D7, 8.0D7,
     $     8.02452D7, 8.02839D7, 8.2D7, 8.4D7, 8.447D7, 8.6D7, 8.6019D7,
     $     8.8D7, 8.80089D7, 9.0D7, 9.01172D7, 9.2D7, 9.31259D7, 9.4D7,
     $     9.4095D7, 9.50365D7, 9.6D7, 9.8D7, 1.0D8, 1.02D8, 1.0212D8,
     $     1.033D8, 1.03728D8, 1.04D8, 1.04999D8, 1.06D8, 1.07592D8,
     $     1.08D8, 1.1D8, 1.10558D8, 1.10955D8, 1.12D8, 1.14D8, 1.147D8,
     $     1.16D8, 1.17583D8, 1.18D8, 1.2D8, 1.22D8, 1.24D8, 1.26D8,
     $     1.26167D8, 1.26582D8, 1.28D8, 1.28465D8, 1.3D8, 1.30283D8,
     $     1.32D8, 1.34D8, 1.35145D8, 1.36D8, 1.37007D8, 1.38D8,
     $     1.38627D8, 1.4D8/

      data ( gdr(18,i), i=1, 412) /2.84751D-4, 2.89D-4, 4.31D-4, 5.44D-4
     $     , 6.65D-4, 8.29D-4, 0.001075, 0.001398, 0.001505, 0.001374,
     $     0.001118, 8.51D-4, 6.87D-4, 7.4D-4, 0.001124, 0.001902,
     $     0.002884, 0.003807, 0.004405, 0.004416, 0.003573, 0.00192,
     $     0.001429, 0.001568, 0.001712, 0.001859, 0.002009, 0.00204978,
     $     0.002054, 0.002159, 0.00231, 0.00246, 0.002608, 0.002752,
     $     0.002893, 0.003028, 0.003156, 0.003278, 0.00339, 0.003493,
     $     0.003585, 0.003609, 0.00361179, 0.003666, 0.003733, 0.003786,
     $     0.003825, 0.003849, 0.00386, 0.003859, 0.003849, 0.00383,
     $     0.003804, 0.003772, 0.003737, 0.003698, 0.003659, 0.00362,
     $     0.003583, 0.00355, 0.003521, 0.003498, 0.003484, 0.003479,
     $     0.003484, 0.003502, 0.003534, 0.003581, 0.003644, 0.003726,
     $     0.003828, 0.003952, 0.004096, 0.004258, 0.004437, 0.004631,
     $     0.004835, 0.005049, 0.00527, 0.005496, 0.005723, 0.005949,
     $     0.006084, 0.006042, 0.005889, 0.005693, 0.00552, 0.005427,
     $     0.005426, 0.005508, 0.005666, 0.005889, 0.006172, 0.006503,
     $     0.006878, 0.007284, 0.007716, 0.008164, 0.008298, 0.00831584,
     $     0.008622, 0.009078, 0.009526, 0.009957, 0.01035, 0.01076,
     $     0.0112, 0.01167, 0.01218, 0.01271, 0.01325, 0.01381, 0.01438,
     $     0.01496, 0.01552, 0.01572, 0.01594, 0.01583, 0.01542, 0.01504
     $     , 0.01497, 0.01531, 0.0157, 0.0161, 0.01651, 0.0168499,
     $     0.01692, 0.01732, 0.01769, 0.01804, 0.0180728, 0.01832,
     $     0.01855, 0.01872, 0.0188, 0.0187914, 0.01879, 0.0187845,
     $     0.01868, 0.01846, 0.0181, 0.01761, 0.01696, 0.01615, 0.01546,
     $     0.0152522, 0.01522, 0.01521, 0.0152, 0.01498, 0.01438,
     $     0.01378, 0.01323, 0.01276, 0.01233, 0.01197, 0.01167, 0.0114,
     $     0.01117, 0.01098, 0.01083, 0.0107, 0.0106, 0.01051, 0.01043,
     $     0.01037, 0.0103, 0.01024, 0.01018, 0.01011, 0.01004, 0.00996,
     $     0.009875, 0.009782, 0.00972297, 0.009612, 0.009329,
     $     0.00931641, 0.009308, 0.009304, 0.009163, 0.00906169,
     $     0.008984, 0.00881, 0.008604, 0.008411, 0.008196, 0.008091,
     $     0.00808513, 0.008001, 0.007796, 0.007603, 0.007427, 0.007258,
     $     0.007151, 0.00715038, 0.00706, 0.007005, 0.00698941, 0.006975
     $     , 0.006989, 0.007038, 0.007151, 0.007318, 0.007576, 0.007916,
     $     0.008372, 0.008902, 0.009472, 0.009931, 0.01015, 0.01001,
     $     0.00973191, 0.009617, 0.009078, 0.008557, 0.008102, 0.007757,
     $     0.007496, 0.007326, 0.007212, 0.007163, 0.007148, 0.00718,
     $     0.007234, 0.007325, 0.007428, 0.007562, 0.007704, 0.007872,
     $     0.008042, 0.008235, 0.00831049, 0.008424, 0.008633, 0.008831,
     $     0.009043, 0.009239, 0.009441, 0.00962, 0.009799, 0.009945,
     $     0.01009, 0.01019, 0.01027, 0.01032, 0.01034, 0.01032, 0.01028
     $     , 0.0101, 0.0100443, 0.009799, 0.00975799, 0.00947883,
     $     0.009409, 0.008956, 0.008461, 0.007952, 0.00780021, 0.007445,
     $     0.006954, 0.00666115, 0.006488, 0.006052, 0.00565, 0.00528,
     $     0.004942, 0.004635, 0.00461128, 0.004356, 0.00423538,
     $     0.004102, 0.003872, 0.003664, 0.003474, 0.00337319, 0.003302,
     $     0.003145, 0.003001, 0.002871, 0.00275, 0.00264, 0.00261297,
     $     0.002539, 0.002447, 0.00236, 0.00233618, 0.002281, 0.002207,
     $     0.0022016, 0.002139, 0.002076, 0.002017, 0.00200123, 0.001963
     $     , 0.001913, 0.001865, 0.001821, 0.00178, 0.00175786,
     $     0.00174614, 0.001745, 0.001739, 0.001704, 0.001669, 0.001639,
     $     0.00161364, 0.00161, 0.001583, 0.00151356, 0.00148621,
     $     0.001384, 0.001269, 0.00124045, 0.001197, 0.00117171,
     $     0.00116566, 0.00116053, 0.001152, 0.00113585, 0.00113353,
     $     0.00112, 0.001098, 0.00109403, 0.0010861, 0.00108439,
     $     0.001084, 0.00107876, 0.001067, 0.001064, 0.001064,
     $     0.00106366, 0.00106365, 0.00106301, 0.001063, 0.001062,
     $     0.00106166, 0.00106054, 0.00106, 0.00105938, 0.00105841,
     $     0.00105717, 0.001057, 0.00105419, 0.001054, 0.001049,
     $     0.00104714, 0.00104558, 0.001042, 0.00104099, 0.00104,
     $     0.00103605, 0.001032, 0.00102396, 0.001021, 0.00101246,
     $     0.00101, 0.00100854, 0.00100831, 9.983D-4, 9.865D-4, 9.83114D
     $     -4, 9.723D-4, 9.72158D-4, 9.577D-4, 9.57637D-4, 9.437D-4,
     $     9.4283D-4, 9.291D-4, 9.20358D-4, 9.137D-4, 9.1297D-4,
     $     9.05803D-4, 8.986D-4, 8.828D-4, 8.67D-4, 8.513D-4, 8.50343D-4
     $     , 8.4103D-4, 8.37702D-4, 8.356D-4, 8.27787D-4, 8.201D-4,
     $     8.07874D-4, 8.048D-4, 7.897D-4, 7.85485D-4, 7.82515D-4,
     $     7.748D-4, 7.599D-4, 7.54731D-4, 7.453D-4, 7.337D-4, 7.307D-4,
     $     7.174D-4, 7.044D-4, 6.914D-4, 6.788D-4, 6.77818D-4, 6.75387D
     $     -4, 6.672D-4, 6.6428D-4, 6.548D-4, 6.53325D-4, 6.445D-4,
     $     6.333D-4, 6.27478D-4, 6.232D-4, 6.18686D-4, 6.143D-4,
     $     6.11163D-4, 6.044D-4/


C     14030, JENDL evaluated data
      nbins(19) = 346
      data ( eg(19,i), i=1, 346) /1.06092D7, 1.06433D7, 1.08D7, 1.09D7
     $     ,1.1D7, 1.11D7, 1.12D7, 1.13D7, 1.14D7, 1.15D7, 1.16D7,
     $     1.17D7, 1.18D7, 1.19D7, 1.2D7, 1.21D7, 1.22D7, 1.23D7, 1.24D7
     $     ,1.25D7, 1.26D7, 1.27D7, 1.28D7, 1.29D7, 1.3D7, 1.31D7,
     $     1.32D7, 1.33D7, 1.34D7, 1.35D7, 1.35063D7, 1.351D7, 1.36D7,
     $     1.37D7,1.38D7, 1.39D7, 1.4D7, 1.41D7, 1.42D7, 1.43D7, 1.44D7,
     $     1.45D7, 1.46D7, 1.47D7, 1.48D7, 1.49D7, 1.5D7, 1.51D7, 1.52D7
     $     ,1.53D7, 1.54D7, 1.55D7, 1.56D7, 1.57D7, 1.58D7, 1.59D7,
     $     1.6D7, 1.61D7, 1.62D7, 1.63D7, 1.64D7, 1.65D7, 1.66D7, 1.67D7
     $     ,1.68D7, 1.69D7, 1.7D7, 1.71D7, 1.72D7, 1.73D7, 1.74D7,
     $     1.75D7, 1.76D7, 1.77D7, 1.78D7, 1.79D7, 1.8D7, 1.81D7, 1.82D7
     $     ,1.83D7, 1.84D7, 1.85D7, 1.86D7, 1.87D7, 1.88D7, 1.89D7,
     $     1.9D7, 1.908D7, 1.90827D7, 1.91D7, 1.92D7, 1.93D7, 1.94D7,
     $     1.95D7,1.96D7, 1.97D7, 1.98D7, 1.99D7, 2.0D7, 2.01D7, 2.02D7,
     $     2.03D7, 2.04D7, 2.05D7, 2.06D7, 2.07D7, 2.0718D7, 2.072D7,
     $     2.08D7,2.09D7, 2.1D7, 2.11D7, 2.12D7, 2.13D7, 2.13502D7,
     $     2.14D7,2.15D7, 2.16D7, 2.17D7, 2.18D7, 2.19D7, 2.2D7,
     $     2.21858D7,2.22D7, 2.24D7, 2.26D7, 2.26264D7, 2.27158D7,
     $     2.28D7, 2.3D7,2.32D7, 2.34D7, 2.36D7, 2.36508D7, 2.38D7,
     $     2.39921D7, 2.4D7,2.42D7, 2.44D7, 2.46D7, 2.46358D7, 2.47776D7
     $     , 2.48D7, 2.5D7,2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7
     $     , 2.64D7, 2.66D7, 2.67489D7, 2.67543D7, 2.68D7, 2.7D7, 2.72D7
     $     , 2.74D7, 2.76D7, 2.78D7, 2.8D7, 2.80936D7, 2.82D7, 2.84D7,
     $     2.86D7, 2.88D7,2.9D7, 2.92D7, 2.94D7, 2.96D7, 2.98D7,
     $     2.98111D7, 2.99975D7,3.0D7, 3.00487D7, 3.01016D7, 3.02D7,
     $     3.04D7, 3.06D7, 3.08D7,3.1D7, 3.12D7, 3.14D7, 3.16D7,
     $     3.16135D7, 3.18D7, 3.2D7,3.22D7, 3.24D7, 3.26D7, 3.28D7,
     $     3.29797D7, 3.3D7, 3.32D7,3.34D7, 3.36D7, 3.38D7, 3.4D7,
     $     3.42D7, 3.44D7, 3.46D7, 3.48D7, 3.5D7, 3.50487D7, 3.52D7,
     $     3.54D7, 3.56D7, 3.58D7, 3.6D7,3.62D7, 3.62624D7, 3.64D7,
     $     3.64662D7, 3.66D7, 3.68D7,3.6801D7, 3.68598D7, 3.7D7,
     $     3.71614D7, 3.72D7, 3.74D7, 3.76D7, 3.76492D7, 3.78D7, 3.8D7,
     $     3.81503D7, 3.82D7, 3.84D7,3.8417D7, 3.86D7, 3.86752D7, 3.88D7
     $     , 3.9D7, 3.92D7, 3.94D7,3.96D7, 3.98D7, 4.0D7, 4.0197D7,
     $     4.07189D7, 4.2D7, 4.28506D7,4.29468D7, 4.39284D7, 4.4D7,
     $     4.50705D7, 4.6D7, 4.61335D7,4.63579D7, 4.72025D7, 4.8D7,
     $     4.95736D7, 5.0D7, 5.01477D7,5.2D7, 5.23511D7, 5.4D7,
     $     5.40186D7, 5.48161D7, 5.48881D7,5.56303D7, 5.56341D7, 5.6D7,
     $     5.8D7, 5.87769D7, 6.0D7,6.03702D7, 6.14083D7, 6.2D7, 6.4D7,
     $     6.4328D7, 6.52538D7,6.6D7, 6.62416D7, 6.8D7, 6.85802D7,
     $     6.86148D7, 6.88113D7,7.0D7, 7.01768D7, 7.2D7, 7.21441D7,
     $     7.4D7, 7.45765D7, 7.6D7,7.61629D7, 7.77202D7, 7.8D7,
     $     7.81574D7, 8.0D7, 8.12392D7,8.14048D7, 8.2D7, 8.36156D7,
     $     8.4D7, 8.58984D7, 8.6D7,8.71412D7, 8.8D7, 8.83077D7, 9.0D7,
     $     9.08134D7, 9.2D7,9.38889D7, 9.4D7, 9.51226D7, 9.6D7, 9.8D7,
     $     9.91292D7,9.98469D7, 1.0D8, 1.02D8, 1.04D8, 1.04704D8,
     $     1.05484D8,1.05557D8, 1.06D8, 1.08D8, 1.1D8, 1.12D8, 1.129D8,
     $     1.13191D8,1.14D8, 1.16D8, 1.1666D8, 1.18D8, 1.1988D8, 1.2D8,
     $     1.21167D8,1.21308D8, 1.21964D8, 1.22D8, 1.24D8, 1.25264D8,
     $     1.26D8,1.2697D8, 1.28D8, 1.3D8, 1.32D8, 1.34D8, 1.35328D8,
     $     1.35401D8, 1.36D8, 1.38D8, 1.39608D8, 1.4D8/

      data ( gdr(19,i), i=1, 346) /9.95102D-4, 9.98304D-4, 0.001013,
     $     5.91D-4, 2.32D-4, 2.38D-4, 5.23D-4, 9.69D-4, 0.001457,
     $     0.001869, 0.002087, 0.001992, 0.001492, 9.31D-4, 0.001017,
     $     0.002313, 0.004162, 0.005327, 0.00481, 0.00337, 0.002527,
     $     0.003146, 0.004772, 0.006789, 0.008584, 0.009541, 0.009046,
     $     0.006484, 0.002734, 0.007013, 0.00730957, 0.007486, 0.01179,
     $     0.01385, 0.0139, 0.01265, 0.01079, 0.009044, 0.008101,
     $     0.008668, 0.0108, 0.01339, 0.01519, 0.01493, 0.01154,
     $     0.006944, 0.005216, 0.006543, 0.009725, 0.01355, 0.01684,
     $     0.01837, 0.01818, 0.01774, 0.01711, 0.01635, 0.01549, 0.0146,
     $     0.0137, 0.01288, 0.01216, 0.01158, 0.01122, 0.0111, 0.01129,
     $     0.01185, 0.01279, 0.01421, 0.01612, 0.0186, 0.02058, 0.02035,
     $     0.01903, 0.01784, 0.01742, 0.01761, 0.01817, 0.01892, 0.01979
     $     , 0.02068, 0.02152, 0.02223, 0.0227, 0.02292, 0.02294, 0.0228
     $     , 0.02255, 0.02224, 0.02224, 0.02224, 0.02192, 0.02164,
     $     0.02146, 0.02139, 0.02139, 0.02145, 0.02156, 0.0217, 0.02186,
     $     0.02199, 0.0221, 0.02216, 0.02213, 0.02203, 0.02181, 0.02146,
     $     0.0214058, 0.0214, 0.02095, 0.02027, 0.0194, 0.0183, 0.01735,
     $     0.01699, 0.0168183, 0.01665, 0.01633, 0.01603, 0.01576,
     $     0.0155, 0.01525, 0.01503, 0.0146672, 0.01464, 0.01429,
     $     0.01401, 0.0139753, 0.0138586, 0.01375, 0.01354, 0.01333,
     $     0.01315, 0.01296, 0.0129165, 0.01279, 0.0126171, 0.01261,
     $     0.01242, 0.01225, 0.01207, 0.0120375, 0.01191, 0.01189,
     $     0.01171, 0.01153, 0.01134, 0.01116, 0.01098, 0.01078, 0.01059
     $     , 0.0104, 0.0102, 0.0100581, 0.010053, 0.01001, 0.009808,
     $     0.009608, 0.009408, 0.009209, 0.009012, 0.008818, 0.0087279,
     $     0.008627, 0.008438, 0.008253, 0.008072, 0.007892, 0.007698,
     $     0.007503, 0.007314, 0.007129, 0.00711894, 0.00695218, 0.00695
     $     , 0.00690708, 0.00686092, 0.006776, 0.006608, 0.006443,
     $     0.006283, 0.006128, 0.00596, 0.005813, 0.005672, 0.00566258,
     $     0.005534, 0.005401, 0.00527, 0.005143, 0.00502, 0.0049,
     $     0.0047956, 0.004784, 0.004671, 0.004562, 0.004456, 0.004352,
     $     0.004252, 0.004155, 0.00406, 0.003968, 0.003878, 0.003791,
     $     0.00377056, 0.003708, 0.003626, 0.003547, 0.003469, 0.003394,
     $     0.003322, 0.00329995, 0.003252, 0.00322961, 0.003185,
     $     0.003119, 0.00311868, 0.00309998, 0.003056, 0.00300665,
     $     0.002995, 0.002935, 0.002877, 0.00286311, 0.002821, 0.002767,
     $     0.00272704, 0.002714, 0.002664, 0.00265971, 0.002614,
     $     0.00259582, 0.002566, 0.002521, 0.002479, 0.00244, 0.002218,
     $     0.002183, 0.002148, 0.00211949, 0.00204639, 0.001881,
     $     0.00180893, 0.00180105, 0.00172346, 0.001718, 0.00166928,
     $     0.001629, 0.00162638, 0.00162201, 0.00160585, 0.001591,
     $     0.00158784, 0.001587, 0.00158775, 0.001597, 0.00159682,
     $     0.001596, 0.00159614, 0.00160217, 0.00160271, 0.00160826,
     $     0.00160829, 0.001611, 0.001626, 0.00163226, 0.001642,
     $     0.0016448, 0.0016526, 0.001657, 0.001668, 0.00166916,
     $     0.00167241, 0.001675, 0.00167537, 0.001678, 0.00168005,
     $     0.00168017, 0.00168086, 0.001685, 0.00168437, 0.001678,
     $     0.00167712, 0.001666, 0.00166104, 0.001649, 0.00164734,
     $     0.00163175, 0.001629, 0.00162691, 0.001603, 0.00158739,
     $     0.00158534, 0.001578, 0.00155529, 0.00155, 0.00152339,
     $     0.001522, 0.00150588, 0.001494, 0.0014893, 0.001464,
     $     0.00145206, 0.001435, 0.00140853, 0.001407, 0.00139057,
     $     0.001378, 0.001351, 0.00133619, 0.00132695, 0.001325,
     $     0.001298, 0.001272, 0.00126273, 0.00125261, 0.00125167,
     $     0.001246, 0.001222, 0.001197, 0.001174, 0.00116354, 0.0011602
     $     , 0.001151, 0.001129, 0.00112131, 0.001106, 0.00108718,
     $     0.001086, 0.00107365, 0.00107218, 0.00106537, 0.001065,
     $     0.001046, 0.00103328, 0.001026, 0.00101769, 0.001009, 9.91D-4
     $     , 9.742D-4, 9.586D-4, 9.48391D-4, 9.47836D-4, 9.433D-4,
     $     9.291D-4, 9.18778D-4, 9.163D-4/


C     15031, JENDL evaluated data
      nbins(20) = 227
      data ( eg(20,i), i=1, 227) /7296000.0, 7400000.0, 7600000.0,
     $     7800000.0, 8000000.0, 8200000.0, 8400000.0, 8600000.0,
     $     8800000.0, 9000000.0, 9200000.0, 9400000.0, 9600000.0,
     $     9800000.0, 1.0D7, 1.02D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7,
     $     1.12D7, 1.14D7, 1.16D7, 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7
     $     , 1.28D7, 1.3D7, 1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7,
     $     1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7
     $     , 1.58D7, 1.6D7, 1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7,
     $     1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7
     $     , 1.88D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7,
     $     2.02D7, 2.04D7, 2.06D7, 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7
     $     , 2.18D7, 2.2D7, 2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7,
     $     2.32D7, 2.34D7, 2.36D7, 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7
     $     , 2.48D7, 2.5D7, 2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7,
     $     2.62D7, 2.64D7, 2.66D7, 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7
     $     , 2.78D7, 2.8D7, 2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7,
     $     2.92D7, 2.94D7, 2.96D7, 2.98D7, 3.0D7, 3.05D7, 3.1D7, 3.15D7,
     $     3.2D7, 3.25D7, 3.3D7, 3.35D7, 3.4D7, 3.45D7, 3.5D7, 3.55D7,
     $     3.6D7, 3.65D7, 3.7D7, 3.75D7, 3.8D7, 3.85D7, 3.9D7, 3.95D7,
     $     4.0D7, 4.05D7, 4.1D7, 4.15D7, 4.2D7, 4.25D7, 4.3D7, 4.35D7,
     $     4.4D7, 4.45D7, 4.5D7, 4.55D7, 4.6D7, 4.65D7, 4.7D7, 4.75D7,
     $     4.8D7, 4.85D7, 4.9D7, 4.95D7, 5.0D7, 5.05D7, 5.1D7, 5.15D7,
     $     5.2D7, 5.25D7, 5.3D7, 5.35D7, 5.4D7, 5.45D7, 5.5D7, 5.55D7,
     $     5.6D7, 5.65D7, 5.7D7, 5.75D7, 5.8D7, 5.85D7, 5.9D7, 5.95D7,
     $     6.0D7, 6.05D7, 6.1D7, 6.15D7, 6.2D7, 6.25D7, 6.3D7, 6.35D7,
     $     6.4D7, 6.45D7, 6.5D7, 6.55D7, 6.6D7, 6.65D7, 6.7D7, 6.75D7,
     $     6.8D7, 6.85D7, 6.9D7, 6.95D7, 7.0D7, 7.05D7, 7.1D7, 7.15D7,
     $     7.2D7, 7.25D7, 7.3D7, 7.35D7, 7.4D7, 7.45D7, 7.5D7, 7.55D7,
     $     7.6D7, 7.65D7, 7.7D7, 7.75D7, 7.8D7, 7.85D7, 7.9D7, 7.95D7,
     $     8.0D7, 8.5D7, 9.0D7, 9.5D7, 1.0D8, 1.05D8, 1.1D8, 1.15D8,
     $     1.2D8, 1.25D8, 1.3D8, 1.35D8, 1.4D8/

      data ( gdr(20,i), i=1, 227) /0.0, 4.1424D-5, 2.2603D-4, 4.0042D-4,
     $     5.937D-4, 8.1154D-4, 0.0010472, 0.001313, 0.0016406,
     $     0.0020483, 0.0025621, 0.0031654, 0.0038899, 0.0048482,
     $     0.0061511, 0.0079764, 0.010616, 0.014542, 0.020451, 0.029013,
     $     0.039466, 0.047022, 0.04592, 0.038137, 0.029763, 0.023497,
     $     0.019399, 0.016941, 0.015702, 0.015492, 0.016368, 0.018692,
     $     0.02333, 0.030734, 0.035653, 0.031992, 0.027163, 0.025509,
     $     0.026843, 0.030956, 0.037681, 0.044381, 0.045631, 0.042297,
     $     0.039617, 0.039257, 0.040825, 0.043707, 0.047362, 0.05125,
     $     0.054789, 0.057405, 0.058737, 0.058784, 0.057961, 0.057057,
     $     0.057124, 0.059198, 0.062862, 0.063383, 0.056927, 0.048333,
     $     0.041772, 0.037581, 0.035147, 0.034123, 0.03443, 0.036159,
     $     0.038841, 0.040004, 0.037638, 0.034033, 0.031215, 0.02932,
     $     0.028014, 0.02703, 0.026205, 0.025451, 0.024725, 0.024003,
     $     0.02328, 0.022557, 0.021844, 0.021154, 0.020501, 0.0199,
     $     0.019368, 0.018917, 0.01856, 0.018292, 0.018108, 0.017973,
     $     0.017824, 0.017583, 0.017188, 0.016628, 0.015954, 0.015244,
     $     0.014571, 0.013983, 0.013501, 0.013128, 0.012852, 0.012658,
     $     0.012512, 0.012374, 0.012198, 0.011945, 0.011594, 0.011152,
     $     0.010647, 0.010115, 0.0095852, 0.0090812, 0.0086121,
     $     0.0076201, 0.0068532, 0.006252, 0.0057691, 0.0053713,
     $     0.0050375, 0.004755, 0.0045138, 0.0043063, 0.004126,
     $     0.0039684, 0.0038309, 0.0037104, 0.0036046, 0.0035109,
     $     0.0034277, 0.003354, 0.0032887, 0.0032307, 0.0031789,
     $     0.0031324, 0.0030909, 0.003054, 0.003021, 0.0029914,
     $     0.0029647, 0.0029408, 0.0029195, 0.0029005, 0.0028834,
     $     0.0028679, 0.0028541, 0.0028419, 0.0028309, 0.0028209,
     $     0.0028118, 0.0028037, 0.0027964, 0.00279, 0.0027843,
     $     0.0027792, 0.0027748, 0.0027708, 0.0027673, 0.0027643,
     $     0.0027615, 0.0027591, 0.002757, 0.0027551, 0.0027534,
     $     0.0027519, 0.0027505, 0.0027493, 0.0027482, 0.002747,
     $     0.0027457, 0.0027445, 0.0027433, 0.0027422, 0.0027411,
     $     0.0027399, 0.0027388, 0.0027376, 0.0027364, 0.0027352,
     $     0.0027339, 0.0027325, 0.0027311, 0.0027296, 0.002728,
     $     0.0027264, 0.0027246, 0.0027228, 0.0027208, 0.0027187,
     $     0.0027165, 0.0027141, 0.0027117, 0.0027091, 0.0027064,
     $     0.0027036, 0.0027007, 0.0026977, 0.0026946, 0.0026914,
     $     0.002688, 0.0026846, 0.002681, 0.0026773, 0.0026735,
     $     0.0026696, 0.0026655, 0.0026614, 0.0026571, 0.0026527,
     $     0.0026482, 0.0026435, 0.0026388, 0.0026339, 0.0026289,
     $     0.0025736, 0.002509, 0.0024479, 0.0023789, 0.0023045,
     $     0.0022271, 0.0021487, 0.0020712, 0.0019966, 0.0019265,
     $     0.0018627, 0.0018065/


C     20040, JENDL evaluated data
      nbins(21) = 230
      data ( eg(21,i), i=1, 230) /7039100.0, 7200000.0, 7400000.0,
     $     7600000.0, 7800000.0, 8000000.0, 8200000.0, 8400000.0,
     $     8600000.0, 8800000.0, 9000000.0, 9200000.0, 9400000.0,
     $     9600000.0, 9800000.0, 1.0D7, 1.02D7, 1.04D7, 1.06D7, 1.08D7,
     $     1.1D7, 1.12D7, 1.14D7, 1.16D7, 1.18D7, 1.2D7, 1.22D7, 1.24D7,
     $     1.26D7, 1.28D7, 1.3D7, 1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7,
     $     1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7
     $     , 1.58D7, 1.6D7, 1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7,
     $     1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7
     $     , 1.88D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7,
     $     2.02D7, 2.04D7, 2.06D7, 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7
     $     , 2.18D7, 2.2D7, 2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7,
     $     2.32D7, 2.34D7, 2.36D7, 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7
     $     , 2.48D7, 2.5D7, 2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7,
     $     2.62D7, 2.64D7, 2.66D7, 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7
     $     , 2.78D7, 2.8D7, 2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7,
     $     2.92D7, 2.94D7, 2.96D7, 2.98D7, 3.0D7, 3.05D7, 3.1D7, 3.15D7,
     $     3.2D7, 3.25D7, 3.3D7, 3.35D7, 3.4D7, 3.45D7, 3.5D7, 3.55D7,
     $     3.6D7, 3.65D7, 3.7D7, 3.75D7, 3.8D7, 3.85D7, 3.9D7, 3.95D7,
     $     4.0D7, 4.05D7, 4.1D7, 4.15D7, 4.2D7, 4.25D7, 4.3D7, 4.35D7,
     $     4.4D7, 4.45D7, 4.5D7, 4.55D7, 4.6D7, 4.65D7, 4.7D7, 4.75D7,
     $     4.8D7, 4.85D7, 4.9D7, 4.95D7, 5.0D7, 5.05D7, 5.1D7, 5.15D7,
     $     5.2D7, 5.25D7, 5.3D7, 5.35D7, 5.4D7, 5.45D7, 5.5D7, 5.55D7,
     $     5.6D7, 5.65D7, 5.7D7, 5.75D7, 5.8D7, 5.85D7, 5.9D7, 5.95D7,
     $     6.0D7, 6.05D7, 6.1D7, 6.15D7, 6.2D7, 6.25D7, 6.3D7, 6.35D7,
     $     6.4D7, 6.45D7, 6.5D7, 6.55D7, 6.6D7, 6.65D7, 6.7D7, 6.75D7,
     $     6.8D7, 6.85D7, 6.9D7, 6.95D7, 7.0D7, 7.05D7, 7.1D7, 7.15D7,
     $     7.2D7, 7.25D7, 7.3D7, 7.35D7, 7.4D7, 7.45D7, 7.5D7, 7.55D7,
     $     7.6D7, 7.65D7, 7.7D7, 7.75D7, 7.8D7, 7.85D7, 7.9D7, 7.95D7,
     $     8.0D7, 8.5D7, 9.0D7, 9.5D7, 1.0D8, 1.05D8, 1.1D8, 1.15D8,
     $     1.2D8, 1.25D8, 1.3D8, 1.35D8, 1.4D8, 1.45D8, 1.5D8/

      data ( gdr(21,i), i=1, 230) /0.0, 3.4291D-5, 8.8292D-5, 1.3886D-4,
     $     1.9127D-4, 2.459D-4, 2.9771D-4, 3.5338D-4, 4.1549D-4, 4.8453D
     $     -4, 5.6107D-4, 6.2995D-4, 7.0122D-4, 7.7902D-4, 8.6396D-4,
     $     9.5681D-4, 0.0010584, 0.0011698, 0.0012922, 0.001427,
     $     0.0015763, 0.0017252, 0.0018871, 0.0020691, 0.0022759,
     $     0.0025141, 0.0027934, 0.0031286, 0.0035433, 0.0040777,
     $     0.0048041, 0.0058376, 0.0074593, 0.010157, 0.014174, 0.01692,
     $     0.015389, 0.013344, 0.013054, 0.01453, 0.017501, 0.020493,
     $     0.020685, 0.018496, 0.016618, 0.015946, 0.016368, 0.017733,
     $     0.020017, 0.023244, 0.027273, 0.03146, 0.034893, 0.037494,
     $     0.040681, 0.046373, 0.055677, 0.066474, 0.073139, 0.07661,
     $     0.081633, 0.08706, 0.090951, 0.094794, 0.098701, 0.0993,
     $     0.094579, 0.087078, 0.080251, 0.07543, 0.072262, 0.069816,
     $     0.067234, 0.064023, 0.060119, 0.05577, 0.051329, 0.047141,
     $     0.043446, 0.040373, 0.037972, 0.036244, 0.03515, 0.034614,
     $     0.034495, 0.034571, 0.034546, 0.034126, 0.033146, 0.031651,
     $     0.02986, 0.02802, 0.026315, 0.024843, 0.023634, 0.022678,
     $     0.021946, 0.021405, 0.021022, 0.020763, 0.020601, 0.020504,
     $     0.020453, 0.020427, 0.020405, 0.02037, 0.020307, 0.020204,
     $     0.020053, 0.019851, 0.019596, 0.019291, 0.018942, 0.018556,
     $     0.018142, 0.017709, 0.016608, 0.015572, 0.014677, 0.013964,
     $     0.01344, 0.013113, 0.012978, 0.01303, 0.013259, 0.013645,
     $     0.014161, 0.01475, 0.01532, 0.015751, 0.015924, 0.015772,
     $     0.015306, 0.014608, 0.013794, 0.012963, 0.012185, 0.011496,
     $     0.01091, 0.010424, 0.010031, 0.0097212, 0.0094866, 0.009319,
     $     0.0092123, 0.0091603, 0.0091605, 0.0092093, 0.0093023,
     $     0.0094323, 0.0095861, 0.0097478, 0.0098933, 0.0099944,
     $     0.010026, 0.0099722, 0.009836, 0.0096338, 0.0093905,
     $     0.0091311, 0.0088752, 0.0086356, 0.0084188, 0.0082269,
     $     0.0080595, 0.0079145, 0.0077895, 0.0076819, 0.0075892,
     $     0.0075093, 0.0074392, 0.0073782, 0.0073251, 0.0072787,
     $     0.0072379, 0.0072019, 0.00717, 0.0071416, 0.0071162,
     $     0.0070934, 0.0070727, 0.0070539, 0.0070367, 0.0070209,
     $     0.0070062, 0.0069925, 0.0069797, 0.0069675, 0.0069559,
     $     0.0069449, 0.0069339, 0.0069233, 0.0069129, 0.0069028,
     $     0.0068928, 0.0068829, 0.006873, 0.0068632, 0.0068534,
     $     0.0068436, 0.0068336, 0.0068236, 0.0068135, 0.0068033,
     $     0.0067929, 0.0067823, 0.0067716, 0.0067607, 0.0067496,
     $     0.0067384, 0.0067268, 0.006715, 0.0067029, 0.0066907,
     $     0.0066783, 0.0066656, 0.0065266, 0.0063657, 0.0062092,
     $     0.0060329, 0.0058429, 0.0056449, 0.0054445, 0.0052465,
     $     0.0050557, 0.0048765, 0.0047132, 0.0045697, 0.0044498,
     $     0.0043571/


C     20048, JENDL evaluated data
      nbins(22) = 201
      data ( eg(22,i), i=1, 201) /1.0D7, 1.02D7, 1.04D7, 1.06D7,
     $     1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7, 1.18D7, 1.2D7, 1.22D7
     $     ,1.24D7, 1.26D7, 1.28D7, 1.3D7, 1.32D7, 1.34D7, 1.36D7,
     $     1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7, 1.52D7
     $     ,1.54D7, 1.56D7, 1.58D7, 1.6D7, 1.62D7, 1.64D7, 1.66D7,
     $     1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7, 1.82D7
     $     ,1.84D7, 1.86D7, 1.88D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7,
     $     1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7, 2.08D7, 2.1D7, 2.12D7
     $     ,2.14D7, 2.16D7, 2.18D7, 2.2D7, 2.22D7, 2.24D7, 2.26D7,
     $     2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7, 2.38D7, 2.4D7, 2.42D7
     $     ,2.44D7, 2.46D7, 2.48D7, 2.5D7, 2.52D7, 2.54D7, 2.56D7,
     $     2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7, 2.68D7, 2.7D7, 2.72D7
     $     ,2.74D7, 2.76D7, 2.78D7, 2.8D7, 2.82D7, 2.84D7, 2.86D7,
     $     2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7, 2.98D7, 3.0D7, 3.02D7
     $     ,3.04D7, 3.06D7, 3.08D7, 3.1D7, 3.12D7, 3.14D7, 3.16D7,
     $     3.18D7, 3.2D7, 3.22D7, 3.24D7, 3.26D7, 3.28D7, 3.3D7, 3.32D7
     $     ,3.34D7, 3.36D7, 3.38D7, 3.4D7, 3.42D7, 3.44D7, 3.46D7,
     $     3.48D7, 3.5D7, 3.52D7, 3.54D7, 3.56D7, 3.58D7, 3.6D7, 3.62D7
     $     ,3.64D7, 3.66D7, 3.68D7, 3.7D7, 3.72D7, 3.74D7, 3.76D7,
     $     3.78D7, 3.8D7, 3.82D7, 3.84D7, 3.86D7, 3.88D7, 3.9D7, 3.92D7
     $     ,3.94D7, 3.96D7, 3.98D7, 4.0D7, 4.2D7, 4.4D7, 4.6D7, 4.8D7
     $     ,5.0D7, 5.2D7, 5.4D7, 5.6D7, 5.8D7, 6.0D7, 6.2D7, 6.4D7,
     $     6.6D7, 6.8D7, 7.0D7, 7.2D7, 7.4D7, 7.6D7, 7.8D7, 8.0D7, 8.2D7
     $     ,8.4D7, 8.6D7, 8.8D7, 9.0D7, 9.2D7, 9.4D7, 9.6D7, 9.8D7,
     $     1.0D8, 1.02D8, 1.04D8, 1.06D8, 1.08D8, 1.1D8, 1.12D8, 1.14D8
     $     ,1.16D8, 1.18D8, 1.2D8, 1.22D8, 1.24D8, 1.26D8, 1.28D8, 1.3D8
     $     ,1.32D8, 1.34D8, 1.36D8, 1.38D8, 1.4D8/

      data ( gdr(22,i), i=1, 201) /2.34D-4, 0.001099, 0.001743, 0.002229
     $     , 0.002724, 0.0035, 0.004705, 0.005984, 0.006874, 0.007385,
     $     0.0080, 0.009076, 0.01038, 0.01155, 0.01248, 0.01334, 0.0143,
     $     0.01536, 0.0165, 0.01774, 0.01907, 0.02052, 0.0221, 0.02383,
     $     0.02569, 0.02773, 0.02994, 0.03235, 0.03496, 0.0378, 0.04088,
     $     0.0442, 0.04777, 0.0516, 0.05569, 0.06002, 0.06456, 0.06926,
     $     0.07408, 0.07894, 0.08373, 0.08834, 0.09263, 0.09646, 0.09971
     $     , 0.1023, 0.104, 0.105, 0.1049, 0.1041, 0.1026, 0.1008,
     $     0.09852, 0.09554, 0.09195, 0.08805, 0.08411, 0.08018, 0.07629
     $     , 0.07247, 0.06878, 0.06524, 0.06188, 0.0587, 0.0557, 0.05289
     $     , 0.05024, 0.04777, 0.04547, 0.04331, 0.04128, 0.03941,
     $     0.03765, 0.03601, 0.03449, 0.03306, 0.03172, 0.03047, 0.0293,
     $     0.0282, 0.02717, 0.0262, 0.02528, 0.02443, 0.02362, 0.02285,
     $     0.02215, 0.02147, 0.02082, 0.02021, 0.01961, 0.01905, 0.01851
     $     , 0.01801, 0.01753, 0.01708, 0.01663, 0.01621, 0.01581,
     $     0.01542, 0.01505, 0.01472, 0.01438, 0.01406, 0.01374, 0.01344
     $     , 0.01316, 0.01289, 0.01263, 0.01237, 0.01213, 0.0119,
     $     0.01168, 0.01147, 0.01126, 0.01107, 0.01088, 0.01069, 0.01052
     $     , 0.01036, 0.0102, 0.01004, 0.009889, 0.009742, 0.009603,
     $     0.00947, 0.009343, 0.009218, 0.009099, 0.008983, 0.008869,
     $     0.008757, 0.008651, 0.00855, 0.008455, 0.008362, 0.008271,
     $     0.008181, 0.008094, 0.00801, 0.00793, 0.007852, 0.007779,
     $     0.007708, 0.007636, 0.007566, 0.007497, 0.00743, 0.007366,
     $     0.007308, 0.007251, 0.006751, 0.006368, 0.006078, 0.005848,
     $     0.005649, 0.005484, 0.005346, 0.005225, 0.005119, 0.005034,
     $     0.00495, 0.004878, 0.004815, 0.004755, 0.0047, 0.004646,
     $     0.004591, 0.004542, 0.004492, 0.004437, 0.004379, 0.004326,
     $     0.004272, 0.00421, 0.004151, 0.00409, 0.004029, 0.003966,
     $     0.003901, 0.003836, 0.003772, 0.003706, 0.003644, 0.003579,
     $     0.003509, 0.003447, 0.003383, 0.003319, 0.003258, 0.003201,
     $     0.003143, 0.003085, 0.003029, 0.002977, 0.002924, 0.002876,
     $     0.00283, 0.002784, 0.00274, 0.002702/


C     22046, JENDL evaluated data
      nbins(23) = 185
      data ( eg(23,i), i=1, 185) /1.33D7, 1.34D7, 1.36D7, 1.38D7,
     $     1.4D7, 1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7, 1.52D7, 1.54D7
     $     ,1.56D7, 1.58D7, 1.6D7, 1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7
     $     ,1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7, 1.82D7, 1.84D7,
     $     1.86D7, 1.88D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7
     $     ,2.02D7, 2.04D7, 2.06D7, 2.08D7, 2.1D7, 2.12D7, 2.14D7,
     $     2.16D7, 2.18D7, 2.2D7, 2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7
     $     ,2.32D7, 2.34D7, 2.36D7, 2.38D7, 2.4D7, 2.42D7, 2.44D7,
     $     2.46D7, 2.48D7, 2.5D7, 2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7
     $     ,2.62D7, 2.64D7, 2.66D7, 2.68D7, 2.7D7, 2.72D7, 2.74D7,
     $     2.76D7, 2.78D7, 2.8D7, 2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7
     $     ,2.92D7, 2.94D7, 2.96D7, 2.98D7, 3.0D7, 3.02D7, 3.04D7,
     $     3.06D7, 3.08D7, 3.1D7, 3.12D7, 3.14D7, 3.16D7, 3.18D7, 3.2D7
     $     ,3.22D7, 3.24D7, 3.26D7, 3.28D7, 3.3D7, 3.32D7, 3.34D7,
     $     3.36D7, 3.38D7, 3.4D7, 3.42D7, 3.44D7, 3.46D7, 3.48D7, 3.5D7
     $     ,3.52D7, 3.54D7, 3.56D7, 3.58D7, 3.6D7, 3.62D7, 3.64D7,
     $     3.66D7, 3.68D7, 3.7D7, 3.72D7, 3.74D7, 3.76D7, 3.78D7, 3.8D7
     $     ,3.82D7, 3.84D7, 3.86D7, 3.88D7, 3.9D7, 3.92D7, 3.94D7,
     $     3.96D7, 3.98D7, 4.0D7, 4.2D7, 4.4D7, 4.6D7, 4.8D7, 5.0D7,
     $     5.2D7,5.4D7, 5.6D7, 5.8D7, 6.0D7, 6.2D7, 6.4D7, 6.6D7, 6.8D7,
     $     7.0D7, 7.2D7, 7.4D7, 7.6D7, 7.8D7, 8.0D7, 8.2D7, 8.4D7, 8.6D7
     $     ,8.8D7, 9.0D7, 9.2D7, 9.4D7, 9.6D7, 9.8D7, 1.0D8, 1.02D8
     $     ,1.04D8, 1.06D8, 1.08D8, 1.1D8, 1.12D8, 1.14D8, 1.16D8,
     $     1.18D8, 1.2D8, 1.22D8, 1.24D8, 1.26D8, 1.28D8, 1.3D8, 1.32D8
     $     ,1.34D8, 1.36D8, 1.38D8, 1.4D8/

      data ( gdr(23,i), i=1, 185) /0.002275, 0.0043, 0.007508, 0.01004,
     $     0.01214, 0.01358, 0.01468, 0.01595, 0.01764, 0.01979, 0.02239
     $     , 0.02523, 0.02809, 0.03037, 0.03111, 0.02969, 0.0274,
     $     0.02603, 0.02621, 0.02715, 0.02819, 0.02969, 0.03218, 0.03537
     $     , 0.03812, 0.03934, 0.03923, 0.03833, 0.03711, 0.03602,
     $     0.03536, 0.03494, 0.03447, 0.03385, 0.03323, 0.03272, 0.03225
     $     , 0.03172, 0.0311, 0.03045, 0.02982, 0.02919, 0.02858,
     $     0.02793, 0.02721, 0.02638, 0.02553, 0.0248, 0.02415, 0.02341,
     $     0.02246, 0.02148, 0.02067, 0.02009, 0.01966, 0.01928, 0.01878
     $     , 0.01802, 0.01704, 0.01609, 0.01527, 0.01469, 0.01427,
     $     0.01394, 0.01367, 0.01339, 0.01311, 0.01281, 0.01249, 0.01216
     $     , 0.01183, 0.0115, 0.01118, 0.01088, 0.01059, 0.01032,
     $     0.01007, 0.009835, 0.009618, 0.009413, 0.009218, 0.009014,
     $     0.008814, 0.008627, 0.008452, 0.008231, 0.008084, 0.00794,
     $     0.007808, 0.007683, 0.007567, 0.00746, 0.007357, 0.007261,
     $     0.00717, 0.007082, 0.006997, 0.006916, 0.006836, 0.006758,
     $     0.006681, 0.006608, 0.006536, 0.006466, 0.006397, 0.006331,
     $     0.006267, 0.006206, 0.006145, 0.006087, 0.006031, 0.005977,
     $     0.005925, 0.005878, 0.005834, 0.005792, 0.00575, 0.00571,
     $     0.005672, 0.005633, 0.005594, 0.005557, 0.005519, 0.005483,
     $     0.005445, 0.00541, 0.005373, 0.005344, 0.005313, 0.005284,
     $     0.005255, 0.005226, 0.005203, 0.005179, 0.005156, 0.004985,
     $     0.004826, 0.004706, 0.004623, 0.004577, 0.004553, 0.004551,
     $     0.004563, 0.004583, 0.004611, 0.004636, 0.004655, 0.004672,
     $     0.004679, 0.004676, 0.004668, 0.004652, 0.004631, 0.004609,
     $     0.004585, 0.004564, 0.004536, 0.004502, 0.004459, 0.004403,
     $     0.00433, 0.004249, 0.004158, 0.004067, 0.003979, 0.003898,
     $     0.003822, 0.00375, 0.003681, 0.003611, 0.003539, 0.003466,
     $     0.003394, 0.003323, 0.003257, 0.003194, 0.003133, 0.003076,
     $     0.003023, 0.002969, 0.002918, 0.002868, 0.002819, 0.002772,
     $     0.002725/


C     23051, JENDL evaluated data
      nbins(24) = 328
      data ( eg(24,i), i=1, 328) /8060610.0, 1.02908D7, 1.10513D7,
     $     1.13D7, 1.14D7, 1.16D7, 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7
     $     , 1.28D7, 1.3D7, 1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7,
     $     1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7
     $     , 1.58D7, 1.6D7, 1.62D7, 1.64D7, 1.66D7, 1.67752D7, 1.68D7,
     $     1.7D7, 1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7, 1.82D7, 1.84D7,
     $     1.86D7, 1.86603D7, 1.88D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7,
     $     1.98D7, 2.0D7, 2.02D7, 2.02232D7, 2.0384D7, 2.04D7, 2.05461D7
     $     , 2.06D7, 2.07415D7, 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7,
     $     2.1635D7, 2.18D7, 2.18251D7, 2.2D7, 2.22D7, 2.23968D7, 2.24D7
     $     , 2.26D7, 2.26029D7, 2.26359D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7
     $     , 2.34559D7, 2.36D7, 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7,
     $     2.48D7, 2.5D7, 2.50777D7, 2.52D7, 2.54D7, 2.54713D7,
     $     2.55245D7, 2.56D7, 2.56361D7, 2.58D7, 2.58256D7, 2.6D7,
     $     2.62D7, 2.64D7, 2.66D7, 2.68D7, 2.68309D7, 2.7D7, 2.72D7,
     $     2.74D7, 2.76D7, 2.78D7, 2.8D7, 2.81647D7, 2.82D7, 2.84D7,
     $     2.85238D7, 2.86D7, 2.87223D7, 2.88D7, 2.9D7, 2.90346D7,
     $     2.92D7, 2.92452D7, 2.94D7, 2.96D7, 2.98D7, 2.98497D7, 3.0D7,
     $     3.01365D7, 3.01911D7, 3.02D7, 3.03D7, 3.04D7, 3.06D7,
     $     3.06819D7, 3.08D7, 3.1D7, 3.1181D7, 3.12D7, 3.12257D7, 3.14D7
     $     , 3.15534D7, 3.16D7, 3.18D7, 3.19368D7, 3.2D7, 3.21568D7,
     $     3.22D7, 3.23246D7, 3.24D7, 3.25487D7, 3.26D7, 3.28D7, 3.3D7,
     $     3.31935D7, 3.32D7, 3.34D7, 3.36D7, 3.38D7, 3.39639D7, 3.4D7,
     $     3.40449D7, 3.42D7, 3.44D7, 3.46D7, 3.47352D7, 3.48D7,
     $     3.4802D7, 3.5D7, 3.5178D7, 3.52D7, 3.52903D7, 3.54D7, 3.56D7,
     $     3.58D7, 3.6D7, 3.62D7, 3.64D7, 3.66D7, 3.68D7, 3.68781D7,
     $     3.7D7, 3.72D7, 3.74D7, 3.76D7, 3.78D7, 3.8D7, 3.82D7, 3.84D7,
     $     3.86D7, 3.86046D7, 3.87284D7, 3.88D7, 3.88364D7, 3.9D7,
     $     3.92D7, 3.94D7, 3.96D7, 3.98D7, 4.0D7, 4.00587D7, 4.04919D7,
     $     4.05569D7, 4.12092D7, 4.2D7, 4.2042D7, 4.24789D7, 4.3141D7,
     $     4.34671D7, 4.4D7, 4.45299D7, 4.49061D7, 4.49132D7, 4.56565D7,
     $     4.6D7, 4.63742D7, 4.76079D7, 4.78033D7, 4.8D7, 4.87935D7,
     $     5.0D7, 5.01102D7, 5.0934D7, 5.2D7, 5.4D7, 5.43085D7,
     $     5.47746D7, 5.49191D7, 5.50543D7, 5.54802D7, 5.6D7, 5.6118D7,
     $     5.61238D7, 5.70542D7, 5.8D7, 5.81472D7, 5.89216D7, 5.96528D7,
     $     6.0D7, 6.08421D7, 6.08949D7, 6.2D7, 6.32622D7, 6.33666D7,
     $     6.4D7, 6.46883D7, 6.6D7, 6.6231D7, 6.72807D7, 6.8D7, 6.8221D7
     $     , 6.87518D7, 7.0D7, 7.0419D7, 7.15956D7, 7.17002D7, 7.2D7,
     $     7.23653D7, 7.32946D7, 7.4D7, 7.50172D7, 7.6D7, 7.7022D7,
     $     7.75355D7, 7.8D7, 8.0D7, 8.11195D7, 8.2D7, 8.26931D7, 8.4D7,
     $     8.44199D7, 8.48507D7, 8.54775D7, 8.58184D7, 8.6D7, 8.70461D7,
     $     8.8D7, 8.99844D7, 9.0D7, 9.03145D7, 9.2D7, 9.20781D7,
     $     9.23403D7, 9.4D7, 9.40056D7, 9.51382D7, 9.6D7, 9.71317D7,
     $     9.8D7, 9.80166D7, 9.87441D7, 1.0D8, 1.01789D8, 1.02D8,
     $     1.026D8, 1.04D8, 1.06D8, 1.06233D8, 1.07325D8, 1.08D8,
     $     1.09382D8, 1.1D8, 1.10172D8, 1.10273D8, 1.12D8, 1.13074D8,
     $     1.13592D8, 1.14D8, 1.16D8, 1.1667D8, 1.16744D8, 1.17923D8,
     $     1.18D8, 1.19018D8, 1.2D8, 1.22D8, 1.24D8, 1.2592D8, 1.26D8,
     $     1.27039D8, 1.28D8, 1.3D8, 1.31053D8, 1.3188D8, 1.32D8,
     $     1.32669D8, 1.34D8, 1.34096D8, 1.34977D8, 1.36D8, 1.36601D8,
     $     1.36964D8, 1.38D8, 1.4D8/

      data ( gdr(24,i), i=1, 328) /8.44581D-4, 0.00107825, 0.00115794,
     $     0.001184, 0.002167, 0.003603, 0.00463, 0.005572, 0.006451,
     $     0.007293, 0.008112, 0.008922, 0.009745, 0.0106, 0.0115,
     $     0.01248, 0.01355, 0.01473, 0.01605, 0.01754, 0.01921, 0.02107
     $     , 0.02319, 0.02559, 0.02831, 0.0314, 0.03491, 0.03888,
     $     0.04336, 0.04835, 0.05382, 0.0589318, 0.05969, 0.06574,
     $     0.0717, 0.07719, 0.08178, 0.08509, 0.08692, 0.08727, 0.08646,
     $     0.08469, 0.0839967, 0.08242, 0.07998, 0.07761, 0.07549,
     $     0.07367, 0.07217, 0.07097, 0.07, 0.0699147, 0.069328, 0.06927
     $     , 0.0688013, 0.06863, 0.0681832, 0.068, 0.06735, 0.06661,
     $     0.06552, 0.06423, 0.0639534, 0.06267, 0.0624417, 0.06088,
     $     0.0589, 0.0568521, 0.05682, 0.05466, 0.0546278, 0.0542609,
     $     0.05248, 0.05027, 0.04807, 0.04592, 0.0453188, 0.04381,
     $     0.0418, 0.03986, 0.03802, 0.03629, 0.03465, 0.03313, 0.03169,
     $     0.0311572, 0.03034, 0.02909, 0.0286665, 0.0283548, 0.02792,
     $     0.0277195, 0.02683, 0.0266969, 0.02581, 0.02487, 0.02396,
     $     0.02314, 0.02237, 0.0222583, 0.02166, 0.02099, 0.02036,
     $     0.01978, 0.01923, 0.01871, 0.0183136, 0.01823, 0.01779,
     $     0.0175223, 0.01736, 0.017114, 0.01696, 0.01658, 0.0165188,
     $     0.01623, 0.0161523, 0.01589, 0.01557, 0.01528, 0.0152098,
     $     0.015, 0.0148151, 0.0147418, 0.01473, 0.01462, 0.01449,
     $     0.01425, 0.0141552, 0.01402, 0.01379, 0.0135998, 0.01358,
     $     0.0135553, 0.01339, 0.0132439, 0.0132, 0.01302, 0.0129102,
     $     0.01286, 0.0127343, 0.0127, 0.0126063, 0.01255, 0.0124456,
     $     0.01241, 0.01228, 0.01215, 0.0120435, 0.01204, 0.01192,
     $     0.01181, 0.01171, 0.0116279, 0.01161, 0.0115897, 0.01152,
     $     0.01143, 0.01135, 0.0112958, 0.01127, 0.0112693, 0.0112,
     $     0.0111377, 0.01113, 0.0110983, 0.01106, 0.011, 0.01094,
     $     0.01088, 0.01082, 0.01077, 0.01072, 0.01068, 0.0106604,
     $     0.01063, 0.01059, 0.01055, 0.01051, 0.01048, 0.01044, 0.01041
     $     , 0.01038, 0.01035, 0.0103493, 0.0103307, 0.01032, 0.0103145,
     $     0.01029, 0.01027, 0.01025, 0.01023, 0.01021, 0.01019,
     $     0.010187, 0.0101649, 0.0101616, 0.0101288, 0.01009, 0.0100891
     $     , 0.0100802, 0.0100669, 0.0100605, 0.01005, 0.01005, 0.01005,
     $     0.01005, 0.01005, 0.01005, 0.0100538, 0.0100661, 0.0100681,
     $     0.01007, 0.0100981, 0.01014, 0.010145, 0.0101824, 0.01023,
     $     0.01033, 0.0103456, 0.010369, 0.0103763, 0.0103831, 0.0104043
     $     , 0.01043, 0.010436, 0.0104363, 0.010483, 0.01053, 0.0105367,
     $     0.0105718, 0.0106045, 0.01062, 0.0106539, 0.010656, 0.0107,
     $     0.0107317, 0.0107343, 0.01075, 0.0107569, 0.01077, 0.0107712,
     $     0.0107764, 0.01078, 0.0107766, 0.0107686, 0.01075, 0.0107415,
     $     0.010718, 0.0107159, 0.01071, 0.0106989, 0.0106709, 0.01065,
     $     0.010609, 0.01057, 0.0105185, 0.0104929, 0.01047, 0.01037,
     $     0.0103023, 0.01025, 0.0102044, 0.01012, 0.0100925, 0.0100646,
     $     0.0100243, 0.0100025, 0.009991, 0.00991871, 0.009854,
     $     0.00971408, 0.009713, 0.00969048, 0.009572, 0.00956624,
     $     0.00954695, 0.009427, 0.00942658, 0.00934326, 0.009281,
     $     0.00919773, 0.009135, 0.00913376, 0.00907968, 0.008988,
     $     0.00885806, 0.008843, 0.0087983, 0.008696, 0.008553,
     $     0.0085358, 0.00845635, 0.008408, 0.00830792, 0.008264,
     $     0.00825171, 0.00824446, 0.008123, 0.00804665, 0.00801035,
     $     0.007982, 0.007845, 0.00779855, 0.00779345, 0.00771315,
     $     0.007708, 0.00764075, 0.007577, 0.007446, 0.00732, 0.00720281
     $     , 0.007198, 0.00713516, 0.007078, 0.006964, 0.00690563,
     $     0.0068605, 0.006854, 0.00681854, 0.006749, 0.00674405,
     $     0.00669881, 0.006647, 0.00661847, 0.00660136, 0.006553,
     $     0.006463/


C     24052, JENDL evaluated data
      nbins(25) = 190
      data ( eg(25,i), i=1, 190) /1.22D7, 1.24D7, 1.26D7, 1.28D7,
     $     1.3D7, 1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7
     $     ,1.46D7, 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7
     $     ,1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7,
     $     1.76D7, 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7
     $     ,1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7,
     $     2.06D7, 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7
     $     ,2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7,
     $     2.36D7, 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7
     $     ,2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7,
     $     2.66D7, 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7
     $     ,2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7,
     $     2.96D7, 2.98D7, 3.0D7, 3.02D7, 3.04D7, 3.06D7, 3.08D7, 3.1D7
     $     ,3.12D7, 3.14D7, 3.16D7, 3.18D7, 3.2D7, 3.22D7, 3.24D7,
     $     3.26D7, 3.28D7, 3.3D7, 3.32D7, 3.34D7, 3.36D7, 3.38D7, 3.4D7
     $     ,3.42D7, 3.44D7, 3.46D7, 3.48D7, 3.5D7, 3.52D7, 3.54D7,
     $     3.56D7, 3.58D7, 3.6D7, 3.62D7, 3.64D7, 3.66D7, 3.68D7, 3.7D7
     $     ,3.72D7, 3.74D7, 3.76D7, 3.78D7, 3.8D7, 3.82D7, 3.84D7,
     $     3.86D7, 3.88D7, 3.9D7, 3.92D7, 3.94D7, 3.96D7, 3.98D7, 4.0D7,
     $     4.2D7, 4.4D7, 4.6D7, 4.8D7, 5.0D7, 5.2D7, 5.4D7, 5.6D7, 5.8D7
     $     ,6.0D7, 6.2D7, 6.4D7, 6.6D7, 6.8D7, 7.0D7, 7.2D7, 7.4D7,
     $     7.6D7, 7.8D7, 8.0D7, 8.2D7, 8.4D7, 8.6D7, 8.8D7, 9.0D7, 9.2D7
     $     ,9.4D7, 9.6D7, 9.8D7, 1.0D8, 1.02D8, 1.04D8, 1.06D8, 1.08D8
     $     ,1.1D8, 1.12D8, 1.14D8, 1.16D8, 1.18D8, 1.2D8, 1.22D8, 1.24D8
     $     ,1.26D8, 1.28D8, 1.3D8, 1.32D8, 1.34D8, 1.36D8, 1.38D8, 1.4D8
     $     /

      data ( gdr(25,i), i=1, 190) /0.001609, 0.003821, 0.006206, 0.00845
     $     , 0.01, 0.01057, 0.011, 0.01234, 0.01477, 0.0175, 0.01987,
     $     0.0218, 0.0233, 0.02457, 0.02601, 0.0281, 0.03158, 0.0372,
     $     0.04466, 0.05252, 0.05949, 0.06542, 0.07043, 0.07483, 0.07914
     $     , 0.08378, 0.08864, 0.0934, 0.09686, 0.09687, 0.09226,
     $     0.08614, 0.08272, 0.08298, 0.08454, 0.08534, 0.08533, 0.08499
     $     , 0.08464, 0.08407, 0.08323, 0.08215, 0.08088, 0.07938,
     $     0.07744, 0.07478, 0.07064, 0.0642, 0.05656, 0.05079, 0.04918,
     $     0.0503, 0.05173, 0.05215, 0.05131, 0.04919, 0.0464, 0.04368,
     $     0.04149, 0.0399, 0.03895, 0.03844, 0.03811, 0.03775, 0.03726,
     $     0.03648, 0.03527, 0.03371, 0.0319, 0.03007, 0.02848, 0.02701,
     $     0.02548, 0.0239, 0.02238, 0.02104, 0.01991, 0.01899, 0.01826,
     $     0.01765, 0.01706, 0.01652, 0.01599, 0.01547, 0.01496, 0.01448
     $     , 0.01402, 0.01358, 0.01317, 0.01278, 0.01251, 0.01217,
     $     0.01184, 0.01152, 0.01122, 0.01094, 0.01068, 0.01041, 0.01015
     $     , 0.009897, 0.009669, 0.009459, 0.009255, 0.009054, 0.008862,
     $     0.008683, 0.008516, 0.008363, 0.00822, 0.00808, 0.007939,
     $     0.0078, 0.007666, 0.00754, 0.00743, 0.00733, 0.007242,
     $     0.007147, 0.007044, 0.006944, 0.006855, 0.006774, 0.006692,
     $     0.006606, 0.006526, 0.006454, 0.006391, 0.006331, 0.006276,
     $     0.006222, 0.006169, 0.006116, 0.006062, 0.00601, 0.005961,
     $     0.005918, 0.005878, 0.00584, 0.005803, 0.005765, 0.005515,
     $     0.005429, 0.005247, 0.00517, 0.005117, 0.00512, 0.005173,
     $     0.005234, 0.005267, 0.00528, 0.005323, 0.005354, 0.005376,
     $     0.005403, 0.005415, 0.005406, 0.005387, 0.005359, 0.005326,
     $     0.005282, 0.005215, 0.005146, 0.005081, 0.004993, 0.004888,
     $     0.004817, 0.004746, 0.004653, 0.004549, 0.004456, 0.004371,
     $     0.004291, 0.004209, 0.004123, 0.004036, 0.003947, 0.003876,
     $     0.003803, 0.003722, 0.003649, 0.003575, 0.003511, 0.003441,
     $     0.003368, 0.00331, 0.003252, 0.0032, 0.003147, 0.003088,
     $     0.003047/


C     26054, JENDL evaluated data
      nbins(26) = 388
      data ( eg(26,i), i=1, 388) /8418860.0, 8853760.0, 1.31D7, 1.32D7
     $     ,1.33D7, 1.33785D7, 1.34D7, 1.35D7, 1.36D7, 1.37D7, 1.38D7
     $     ,1.39D7, 1.4D7, 1.41D7, 1.42D7, 1.43D7, 1.44D7, 1.45D7,
     $     1.46D7, 1.47D7, 1.48D7, 1.49D7, 1.5D7, 1.51D7, 1.52D7, 1.53D7
     $     ,1.54D7, 1.54136D7, 1.55D7, 1.56D7, 1.57D7, 1.58D7, 1.59D7
     $     ,1.6D7, 1.61D7, 1.62D7, 1.63D7, 1.64D7, 1.65D7, 1.66D7,
     $     1.67D7, 1.67967D7, 1.68D7, 1.69D7, 1.7D7, 1.70647D7, 1.71D7,
     $     1.72D7, 1.73D7, 1.74D7, 1.75D7, 1.76D7, 1.77D7, 1.77016D7,
     $     1.78D7,1.79D7, 1.8D7, 1.81D7, 1.82D7, 1.83D7, 1.84D7, 1.85D7,
     $     1.86D7, 1.8683D7, 1.87D7, 1.88D7, 1.89D7, 1.9D7, 1.91D7,
     $     1.92D7,1.93D7, 1.94D7, 1.95D7, 1.96D7, 1.97D7, 1.97349D7,
     $     1.98D7,1.99D7, 1.99711D7, 2.0D7, 2.01D7, 2.02D7, 2.03D7,
     $     2.04D7,2.05D7, 2.06D7, 2.07D7, 2.08D7, 2.09D7, 2.1D7, 2.11D7,
     $     2.12D7, 2.13D7, 2.14D7, 2.15D7, 2.16D7, 2.17D7, 2.18D7,
     $     2.19D7,2.2D7, 2.21D7, 2.22D7, 2.23D7, 2.23092D7, 2.24D7,
     $     2.25D7,2.2543D7, 2.26D7, 2.27D7, 2.28D7, 2.29D7, 2.29612D7,
     $     2.3D7,2.31D7, 2.32D7, 2.33D7, 2.34D7, 2.35D7, 2.36D7, 2.37D7
     $     ,2.38D7, 2.39D7, 2.4D7, 2.40619D7, 2.41D7, 2.42D7, 2.42359D7
     $     ,2.43D7, 2.44D7, 2.44222D7, 2.44914D7, 2.45D7, 2.46D7, 2.47D7
     $     ,2.48D7, 2.49D7, 2.5D7, 2.51D7, 2.52D7, 2.52359D7, 2.53D7
     $     ,2.54D7, 2.55D7, 2.55768D7, 2.56D7, 2.57D7, 2.58D7, 2.58601D7
     $     ,2.5862D7, 2.59D7, 2.59178D7, 2.6D7, 2.6035D7, 2.61D7
     $     ,2.61359D7, 2.62D7, 2.63D7, 2.64D7, 2.65D7, 2.66D7, 2.67D7
     $     ,2.68D7, 2.69D7, 2.7D7, 2.70862D7, 2.71D7, 2.72D7, 2.73D7
     $     ,2.74D7, 2.75D7, 2.75948D7, 2.76D7, 2.77D7, 2.78D7, 2.79D7
     $     ,2.8D7, 2.81D7, 2.82D7, 2.83D7, 2.84D7, 2.84905D7, 2.85D7
     $     ,2.85891D7, 2.86D7, 2.87D7, 2.87288D7, 2.88D7, 2.89D7, 2.9D7
     $     ,2.91D7, 2.91522D7, 2.92D7, 2.9215D7, 2.93D7, 2.94D7, 2.95D7
     $     ,2.96D7, 2.97D7, 2.98D7, 2.99D7, 3.0D7, 3.02D7, 3.04D7
     $     ,3.04779D7, 3.04834D7, 3.06D7, 3.08D7, 3.1D7, 3.10272D7
     $     ,3.12D7, 3.13067D7, 3.14D7, 3.16D7, 3.18D7, 3.2D7, 3.22D7
     $     ,3.23511D7, 3.24D7, 3.26D7, 3.28D7, 3.3D7, 3.31302D7, 3.32D7
     $     ,3.34D7, 3.36D7, 3.38D7, 3.39784D7, 3.4D7, 3.42D7, 3.44D7
     $     ,3.46D7, 3.48D7, 3.5D7, 3.52D7, 3.54D7, 3.56D7, 3.58D7, 3.6D7
     $     ,3.62D7, 3.64D7, 3.65338D7, 3.66D7, 3.68D7, 3.7D7, 3.72D7
     $     ,3.74D7, 3.74108D7, 3.74964D7, 3.76D7, 3.78D7, 3.8D7, 3.82D7
     $     ,3.84D7, 3.85652D7, 3.86D7, 3.88D7, 3.88453D7, 3.9D7, 3.92D7
     $     ,3.932D7, 3.94D7, 3.95547D7, 3.96D7, 3.97846D7, 3.98D7, 4.0D7
     $     ,4.01207D7, 4.02451D7, 4.16671D7, 4.2D7, 4.4D7, 4.41684D7
     $     ,4.43495D7, 4.46159D7, 4.56864D7, 4.6D7, 4.6141D7, 4.62919D7
     $     ,4.71021D7, 4.78063D7, 4.78835D7, 4.8D7, 4.82076D7, 4.93287D7
     $     ,5.0D7, 5.19405D7, 5.2D7, 5.38742D7, 5.4D7, 5.40622D7
     $     ,5.49288D7, 5.54719D7, 5.57675D7, 5.6D7, 5.62039D7, 5.74984D7
     $     ,5.8D7, 5.83755D7, 5.93858D7, 6.0D7, 6.01372D7, 6.09418D7
     $     ,6.2D7, 6.33537D7, 6.4D7, 6.41244D7, 6.54526D7, 6.6D7, 6.8D7
     $     ,6.87484D7, 6.91048D7, 6.91152D7, 7.0D7, 7.10473D7, 7.15743D7
     $     ,7.2D7, 7.2023D7, 7.24884D7, 7.4D7, 7.41939D7, 7.47332D7
     $     ,7.6D7, 7.65339D7, 7.76548D7, 7.8D7, 7.90203D7, 8.0D7, 8.2D7
     $     ,8.21608D7, 8.39577D7, 8.4D7, 8.40245D7, 8.48394D7, 8.6D7
     $     ,8.65684D7, 8.72806D7, 8.8D7, 8.82433D7, 8.83931D7, 9.0D7
     $     ,9.12856D7, 9.15234D7, 9.2D7, 9.35208D7, 9.4D7, 9.53059D7
     $     ,9.6D7, 9.8D7, 9.81754D7, 1.0D8, 1.02D8, 1.02806D8, 1.04D8
     $     ,1.06D8, 1.06125D8, 1.0682D8, 1.06961D8, 1.08D8, 1.08185D8
     $     ,1.08747D8, 1.1D8, 1.11211D8, 1.12D8, 1.14D8, 1.15852D8
     $     ,1.16D8, 1.18D8, 1.19573D8, 1.2D8, 1.20145D8, 1.21574D8
     $     ,1.22D8, 1.22336D8, 1.22468D8, 1.24D8, 1.25847D8, 1.26D8
     $     ,1.26529D8, 1.28D8, 1.3D8, 1.31182D8, 1.32D8, 1.34D8, 1.36D8
     $     ,1.37283D8, 1.37814D8, 1.38D8, 1.4D8/

      data ( gdr(26,i), i=1, 388) /0.0012699, 0.0013355, 0.001976,
     $     0.003952, 0.005929, 0.0074315, 0.007905, 0.008642, 0.009377,
     $     0.009736, 0.01088, 0.0124, 0.01417, 0.01613, 0.01822, 0.02042
     $     , 0.02269, 0.02502, 0.02737, 0.02972, 0.03232, 0.03537,
     $     0.03871, 0.04216, 0.04554, 0.04876, 0.05092, 0.0510468,
     $     0.05186, 0.05212, 0.05231, 0.05287, 0.05413, 0.05623, 0.05885
     $     , 0.06211, 0.0657, 0.06936, 0.07271, 0.07564, 0.07815,
     $     0.0804205, 0.0805, 0.083, 0.08582, 0.0885888, 0.09013,
     $     0.09468, 0.09936, 0.1045, 0.1105, 0.1182, 0.1282, 0.12837,
     $     0.1389, 0.1464, 0.1455, 0.1363, 0.1248, 0.1161, 0.1128,
     $     0.1151, 0.121, 0.124642, 0.1254, 0.1232, 0.1153, 0.107,
     $     0.1021, 0.1019, 0.1061, 0.1135, 0.1221, 0.1586, 0.1619,
     $     0.160492, 0.1579, 0.1486, 0.140977, 0.138, 0.129, 0.1229,
     $     0.12, 0.1199, 0.1222, 0.126, 0.13, 0.133, 0.1338, 0.1315,
     $     0.1264, 0.1191, 0.1105, 0.1015, 0.09299, 0.08533, 0.07872,
     $     0.07328, 0.06904, 0.06598, 0.06416, 0.06356, 0.06404,
     $     0.0641517, 0.06526, 0.06629, 0.0658053, 0.06517, 0.06214,
     $     0.05737, 0.05171, 0.0481693, 0.04606, 0.04053, 0.03607,
     $     0.03237, 0.02931, 0.02676, 0.02461, 0.02278, 0.02122, 0.01986
     $     , 0.01867, 0.0180122, 0.01762, 0.01668, 0.016369, 0.01583,
     $     0.01507, 0.0149184, 0.0144565, 0.0144, 0.01378, 0.01323,
     $     0.01272, 0.01225, 0.01182, 0.01143, 0.01106, 0.0109402,
     $     0.01073, 0.01042, 0.01014, 0.00992955, 0.009867, 0.009616,
     $     0.009383, 0.00924892, 0.00924471, 0.009161, 0.00912296,
     $     0.00895, 0.00888047, 0.008753, 0.00868573, 0.008567, 0.008388
     $     , 0.008218, 0.008044, 0.007877, 0.007718, 0.007567, 0.007419,
     $     0.007283, 0.00716822, 0.00715, 0.007022, 0.006901, 0.006787,
     $     0.006674, 0.00657349, 0.006568, 0.006466, 0.006367, 0.006273,
     $     0.006182, 0.006095, 0.006012, 0.00593, 0.005854, 0.00578607,
     $     0.005779, 0.0057148, 0.005707, 0.005636, 0.00561661, 0.005569
     $     , 0.005505, 0.005441, 0.005381, 0.00535012, 0.005322,
     $     0.00531356, 0.005266, 0.005213, 0.005161, 0.005108, 0.005058,
     $     0.005007, 0.004956, 0.004908, 0.004737, 0.004652, 0.00461943,
     $     0.00461716, 0.004569, 0.004493, 0.004423, 0.00441352,
     $     0.004354, 0.00432077, 0.004292, 0.004234, 0.004178, 0.004124,
     $     0.004072, 0.00403112, 0.004018, 0.003965, 0.003912, 0.003856,
     $     0.00382008, 0.003801, 0.003743, 0.003684, 0.003626,
     $     0.00357598, 0.00357, 0.003511, 0.003455, 0.003398, 0.003344,
     $     0.00329, 0.003238, 0.003187, 0.00314, 0.003095, 0.003051,
     $     0.003009, 0.002971, 0.00294618, 0.002934, 0.0029, 0.002867,
     $     0.002836, 0.002807, 0.00280553, 0.00279394, 0.00278, 0.002754
     $     , 0.002731, 0.00271, 0.00269, 0.00267429, 0.002671, 0.002653,
     $     0.00264936, 0.002637, 0.002621, 0.00261379, 0.002609,
     $     0.00259816, 0.002595, 0.00258392, 0.002583, 0.002573,
     $     0.00257176, 0.00257049, 0.00255625, 0.002553, 0.002524,
     $     0.00252885, 0.00253405, 0.00254169, 0.00257214, 0.002581,
     $     0.00259185, 0.00260347, 0.00266611, 0.00272088, 0.0027269,
     $     0.002736, 0.00275777, 0.0028767, 0.002949, 0.00317397,
     $     0.003181, 0.0033851, 0.003399, 0.00340537, 0.00349455,
     $     0.00355088, 0.00358169, 0.003606, 0.00362474, 0.00374435,
     $     0.003791, 0.00381763, 0.00388935, 0.003933, 0.00393825,
     $     0.00396892, 0.004009, 0.00401308, 0.004015, 0.00401208,
     $     0.0039814, 0.003969, 0.003891, 0.00386327, 0.00385024,
     $     0.00384987, 0.003818, 0.00377721, 0.00375708, 0.003741,
     $     0.00374004, 0.00372085, 0.00366, 0.00365177, 0.00362909,
     $     0.003577, 0.00355389, 0.00350637, 0.003492, 0.00344447,
     $     0.0034, 0.003307, 0.00329926, 0.00321494, 0.003213,
     $     0.00321182, 0.00317294, 0.003119, 0.00309353, 0.00306213,
     $     0.003031, 0.00301992, 0.00301314, 0.002942, 0.00288752,
     $     0.00287763, 0.002858, 0.0027945, 0.002775, 0.00272694,
     $     0.002702, 0.002631, 0.00262455, 0.002559, 0.002492,
     $     0.00246628, 0.002429, 0.00237, 0.00236625, 0.00234551,
     $     0.00234131, 0.002311, 0.0023057, 0.00228979, 0.002255,
     $     0.00222202, 0.002201, 0.00215, 0.00210269, 0.002099, 0.002052
     $     , 0.00201488, 0.002005, 0.00200183, 0.00197102, 0.001962,
     $     0.0019545, 0.00195154, 0.001918, 0.00188101, 0.001878,
     $     0.00186755, 0.001839, 0.001803, 0.00178217, 0.001768,
     $     0.001734, 0.001702, 0.00168264, 0.00167474, 0.001672,
     $     0.001645/


C     26056, JENDL evaluated data
      nbins(27) = 319
      data ( eg(27,i), i=1, 319) /7613120.0, 1.01836D7, 1.11973D7,
     $     1.14D7, 1.16D7, 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7
     $     , 1.3D7, 1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7,
     $     1.44D7, 1.46D7, 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7
     $     , 1.6D7, 1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7, 1.70557D7,
     $     1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7, 1.81855D7, 1.82D7,
     $     1.82506D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7, 1.91319D7, 1.92D7,
     $     1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.00824D7, 2.02D7, 2.02516D7,
     $     2.04D7, 2.04952D7, 2.05425D7, 2.06D7, 2.08D7, 2.08241D7,
     $     2.08672D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7,
     $     2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7
     $     , 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7,
     $     2.50741D7, 2.51042D7, 2.51436D7, 2.52D7, 2.53176D7, 2.54D7,
     $     2.54698D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.65502D7,
     $     2.66D7, 2.68D7, 2.7D7, 2.70169D7, 2.72D7, 2.74D7, 2.75496D7,
     $     2.75525D7, 2.76D7, 2.78D7, 2.78125D7, 2.79407D7, 2.8D7,
     $     2.81996D7, 2.82D7, 2.84D7, 2.84431D7, 2.86D7, 2.86851D7,
     $     2.88D7, 2.9D7, 2.90822D7, 2.92D7, 2.93534D7, 2.94D7, 2.96D7,
     $     2.98D7, 3.0D7, 3.02D7, 3.03891D7, 3.04D7, 3.0439D7, 3.06D7,
     $     3.06233D7, 3.08D7, 3.08931D7, 3.1D7, 3.12D7, 3.14D7, 3.16D7,
     $     3.18D7, 3.2D7, 3.22D7, 3.24D7, 3.26D7, 3.27494D7, 3.28D7,
     $     3.29865D7, 3.3D7, 3.31647D7, 3.32D7, 3.34D7, 3.36D7, 3.38D7,
     $     3.38737D7, 3.4D7, 3.41675D7, 3.41854D7, 3.42D7, 3.44D7,
     $     3.44838D7, 3.46D7, 3.4665D7, 3.48D7, 3.5D7, 3.52D7, 3.53367D7
     $     , 3.54D7, 3.56D7, 3.58D7, 3.58277D7, 3.6D7, 3.61541D7, 3.62D7
     $     , 3.64D7, 3.66D7, 3.68D7, 3.7D7, 3.70725D7, 3.72D7, 3.74D7,
     $     3.76D7, 3.78D7, 3.8D7, 3.82D7, 3.84D7, 3.86D7, 3.88D7, 3.9D7,
     $     3.90292D7, 3.92D7, 3.94D7, 3.9449D7, 3.96D7, 3.98D7, 4.0D7,
     $     4.02929D7, 4.09756D7, 4.13858D7, 4.17681D7, 4.19166D7, 4.2D7,
     $     4.21289D7, 4.3551D7, 4.4D7, 4.45572D7, 4.49164D7, 4.56145D7,
     $     4.6D7, 4.7175D7, 4.8D7, 4.81291D7, 4.88702D7, 4.89535D7,
     $     4.9198D7, 4.93837D7, 5.0D7, 5.14803D7, 5.15071D7, 5.16622D7,
     $     5.17342D7, 5.2D7, 5.21438D7, 5.3827D7, 5.4D7, 5.43613D7,
     $     5.6D7, 5.65773D7, 5.8D7, 5.98433D7, 6.0D7, 6.03868D7,
     $     6.07403D7, 6.09257D7, 6.18063D7, 6.2D7, 6.3466D7, 6.4D7,
     $     6.47634D7, 6.48811D7, 6.49391D7, 6.49424D7, 6.5836D7, 6.6D7,
     $     6.6867D7, 6.8D7, 6.95233D7, 7.0D7, 7.2D7, 7.30296D7, 7.4D7,
     $     7.41506D7, 7.43685D7, 7.45574D7, 7.46658D7, 7.53625D7, 7.6D7,
     $     7.71563D7, 7.73964D7, 7.8D7, 7.90991D7, 7.99401D7, 8.0D7,
     $     8.02738D7, 8.13037D7, 8.2D7, 8.4D7, 8.6D7, 8.6851D7, 8.8D7,
     $     9.0D7, 9.01493D7, 9.08516D7, 9.12009D7, 9.17219D7, 9.2D7,
     $     9.25183D7, 9.4D7, 9.45452D7, 9.55931D7, 9.56908D7, 9.6D7,
     $     9.8D7, 9.89326D7, 9.9041D7, 1.0D8, 1.02D8, 1.02423D8, 1.04D8,
     $     1.05099D8, 1.06D8, 1.0665D8, 1.07064D8, 1.07986D8, 1.08D8,
     $     1.08737D8, 1.1D8, 1.11D8, 1.12D8, 1.1305D8, 1.14D8, 1.15848D8
     $     , 1.16D8, 1.18D8, 1.19375D8, 1.2D8, 1.22D8, 1.22111D8,
     $     1.22806D8, 1.23668D8, 1.24D8, 1.26D8, 1.26337D8, 1.2662D8,
     $     1.28D8, 1.28059D8, 1.3D8, 1.30135D8, 1.32D8, 1.34D8, 1.36D8,
     $     1.38D8, 1.38786D8, 1.4D8/

      data ( gdr(27,i), i=1, 319) /8.54806D-4, 0.00114342, 0.00125724,
     $     0.00128, 0.00249, 0.0037, 0.005029, 0.006436, 0.007881,
     $     0.009343, 0.01081, 0.01228, 0.01361, 0.01502, 0.01664, 0.0187
     $     , 0.02117, 0.02404, 0.02732, 0.03101, 0.03511, 0.03962,
     $     0.04442, 0.04942, 0.05444, 0.05926, 0.0637, 0.06764, 0.07105,
     $     0.07394, 0.07633, 0.07816, 0.0784742, 0.07929, 0.07957,
     $     0.07897, 0.07766, 0.07601, 0.0745339, 0.07442, 0.0741106,
     $     0.07321, 0.07255, 0.07234, 0.07231, 0.0721513, 0.07207,
     $     0.07129, 0.06979, 0.06767, 0.06518, 0.0640929, 0.06258,
     $     0.0619338, 0.06012, 0.0590156, 0.0584773, 0.05783, 0.05569,
     $     0.0554264, 0.0549596, 0.05355, 0.05128, 0.04878, 0.04591,
     $     0.04291, 0.0399, 0.03696, 0.03414, 0.0315, 0.02905, 0.0268,
     $     0.02476, 0.02294, 0.02131, 0.01985, 0.01854, 0.01738, 0.01633
     $     , 0.0154, 0.01456, 0.01379, 0.0135333, 0.0134308, 0.0132978,
     $     0.01311, 0.0127353, 0.01248, 0.0122776, 0.01191, 0.01139,
     $     0.01091, 0.01047, 0.01007, 0.00979292, 0.009703, 0.009362,
     $     0.009047, 0.00902193, 0.008755, 0.008486, 0.0082981,
     $     0.00829454, 0.008236, 0.0080, 0.00798588, 0.00784233,
     $     0.007777, 0.00756636, 0.007566, 0.007371, 0.00733062,
     $     0.007186, 0.00711261, 0.007015, 0.006853, 0.00678996,
     $     0.006701, 0.00659249, 0.00656, 0.006426, 0.006301, 0.006184,
     $     0.006073, 0.00597459, 0.005969, 0.00594753, 0.00586,
     $     0.00584775, 0.005756, 0.00571055, 0.005659, 0.005566, 0.00548
     $     , 0.005398, 0.005322, 0.005248, 0.005179, 0.005112, 0.005048,
     $     0.00500157, 0.004986, 0.00493001, 0.004926, 0.00487735,
     $     0.004867, 0.004811, 0.004755, 0.004703, 0.00468373, 0.004651,
     $     0.00460908, 0.00460462, 0.004601, 0.004553, 0.00453322,
     $     0.004506, 0.00449163, 0.004462, 0.004419, 0.004378,
     $     0.00435198, 0.00434, 0.004304, 0.004269, 0.00426427, 0.004235
     $     , 0.00421108, 0.004204, 0.004174, 0.004144, 0.004116, 0.00409
     $     , 0.00408053, 0.004064, 0.00404, 0.004017, 0.003998, 0.003979
     $     , 0.003962, 0.003945, 0.00393, 0.003916, 0.003901, 0.00389924
     $     , 0.003889, 0.003878, 0.00387579, 0.003869, 0.003861,
     $     0.003851, 0.00384499, 0.00383119, 0.00382303, 0.00381552,
     $     0.00381262, 0.003811, 0.00381847, 0.00390035, 0.003926,
     $     0.00400547, 0.00405703, 0.0041579, 0.004214, 0.00444609,
     $     0.004613, 0.00464031, 0.00479879, 0.00481679, 0.00486983,
     $     0.00491031, 0.005046, 0.00534501, 0.0053505, 0.00538234,
     $     0.00539715, 0.005452, 0.00547543, 0.00575228, 0.005781,
     $     0.00582552, 0.006028, 0.00606801, 0.006166, 0.00618261,
     $     0.006184, 0.00617002, 0.00615735, 0.00615074, 0.00611974,
     $     0.006113, 0.0060, 0.00596, 0.00589156, 0.00588115, 0.00587602
     $     , 0.00587573, 0.00579803, 0.005784, 0.00570683, 0.005609,
     $     0.00548707, 0.00545, 0.005304, 0.00523675, 0.005175,
     $     0.00516521, 0.0051511, 0.00513894, 0.00513199, 0.00508777,
     $     0.005048, 0.00497906, 0.004965, 0.00493, 0.00486387,
     $     0.00481448, 0.004811, 0.00479563, 0.0047387, 0.004701,
     $     0.00459, 0.004487, 0.00444519, 0.00439, 0.004294, 0.00428684,
     $     0.00425347, 0.00423706, 0.00421283, 0.0042, 0.0041771,
     $     0.004113, 0.00408836, 0.0040418, 0.00403751, 0.004024,
     $     0.003936, 0.00389734, 0.0038929, 0.003854, 0.003774,
     $     0.00375746, 0.003697, 0.00365431, 0.00362, 0.00359565,
     $     0.00358029, 0.00354651, 0.003546, 0.00351914, 0.003474,
     $     0.00343816, 0.003403, 0.0033675, 0.003336, 0.00327215,
     $     0.003267, 0.003203, 0.00315944, 0.00314, 0.003081, 0.00307768
     $     , 0.00305697, 0.00303165, 0.003022, 0.002965, 0.00295594,
     $     0.00294838, 0.002912, 0.00291041, 0.002859, 0.00285551,
     $     0.002808, 0.002761, 0.002716, 0.002673, 0.00265714, 0.002633/


C     27059, JENDL evaluated data
      nbins(28) = 198
      data ( eg(28,i), i=1, 198) /1.06D7, 1.08D7, 1.1D7, 1.12D7,
     $     1.14D7, 1.16D7, 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7
     $     ,1.3D7, 1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7
     $     ,1.46D7, 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7
     $     ,1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7,
     $     1.76D7, 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7
     $     ,1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7,
     $     2.06D7, 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7
     $     ,2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7,
     $     2.36D7, 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7
     $     ,2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7,
     $     2.66D7, 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7
     $     ,2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7,
     $     2.96D7, 2.98D7, 3.0D7, 3.02D7, 3.04D7, 3.06D7, 3.08D7, 3.1D7
     $     ,3.12D7, 3.14D7, 3.16D7, 3.18D7, 3.2D7, 3.22D7, 3.24D7,
     $     3.26D7, 3.28D7, 3.3D7, 3.32D7, 3.34D7, 3.36D7, 3.38D7, 3.4D7
     $     ,3.42D7, 3.44D7, 3.46D7, 3.48D7, 3.5D7, 3.52D7, 3.54D7,
     $     3.56D7, 3.58D7, 3.6D7, 3.62D7, 3.64D7, 3.66D7, 3.68D7, 3.7D7
     $     ,3.72D7, 3.74D7, 3.76D7, 3.78D7, 3.8D7, 3.82D7, 3.84D7,
     $     3.86D7, 3.88D7, 3.9D7, 3.92D7, 3.94D7, 3.96D7, 3.98D7, 4.0D7,
     $     4.2D7, 4.4D7, 4.6D7, 4.8D7, 5.0D7, 5.2D7, 5.4D7, 5.6D7, 5.8D7
     $     ,6.0D7, 6.2D7, 6.4D7, 6.6D7, 6.8D7, 7.0D7, 7.2D7, 7.4D7,
     $     7.6D7, 7.8D7, 8.0D7, 8.2D7, 8.4D7, 8.6D7, 8.8D7, 9.0D7, 9.2D7
     $     ,9.4D7, 9.6D7, 9.8D7, 1.0D8, 1.02D8, 1.04D8, 1.06D8, 1.08D8
     $     ,1.1D8, 1.12D8, 1.14D8, 1.16D8, 1.18D8, 1.2D8, 1.22D8, 1.24D8
     $     ,1.26D8, 1.28D8, 1.3D8, 1.32D8, 1.34D8, 1.36D8, 1.38D8, 1.4D8
     $     /

      data ( gdr(28,i), i=1, 198) /0.003634, 0.004824, 0.0045, 0.005131,
     $     0.006639, 0.008327, 0.009813, 0.01104, 0.012, 0.01285,
     $     0.01379, 0.01489, 0.01612, 0.01749, 0.01898, 0.02064, 0.02252
     $     , 0.02463, 0.027, 0.02967, 0.03268, 0.03608, 0.04, 0.04447,
     $     0.04933, 0.05436, 0.05958, 0.06529, 0.07153, 0.07708, 0.08048
     $     , 0.08135, 0.08034, 0.07818, 0.07539, 0.07248, 0.07041,
     $     0.07067, 0.07404, 0.07838, 0.08086, 0.0807, 0.07915, 0.07741,
     $     0.0756, 0.07365, 0.07142, 0.06879, 0.06569, 0.06266, 0.06027,
     $     0.05863, 0.05733, 0.056, 0.05454, 0.05288, 0.05122, 0.04973,
     $     0.04864, 0.04778, 0.04689, 0.04584, 0.04475, 0.04362, 0.04255
     $     , 0.04174, 0.04118, 0.04072, 0.04029, 0.03969, 0.03877,
     $     0.03762, 0.03657, 0.03585, 0.03537, 0.03486, 0.03427, 0.03365
     $     , 0.03303, 0.03244, 0.03184, 0.03124, 0.03064, 0.03005,
     $     0.02946, 0.02887, 0.02829, 0.02772, 0.02715, 0.02656, 0.02597
     $     , 0.0254, 0.02483, 0.02427, 0.02371, 0.02315, 0.0226, 0.02208
     $     , 0.02169, 0.02129, 0.0209, 0.02052, 0.0201, 0.01961, 0.01911
     $     , 0.01859, 0.01807, 0.01759, 0.01716, 0.01677, 0.01642,
     $     0.0161, 0.01577, 0.01545, 0.01513, 0.0148, 0.01442, 0.01402,
     $     0.01357, 0.01311, 0.01263, 0.01217, 0.01173, 0.01134, 0.01098
     $     , 0.01065, 0.01034, 0.01003, 0.00972, 0.009414, 0.009121,
     $     0.008852, 0.008616, 0.008416, 0.008245, 0.008088, 0.007934,
     $     0.007767, 0.007579, 0.007385, 0.007199, 0.007043, 0.006928,
     $     0.006866, 0.006834, 0.006804, 0.006765, 0.006724, 0.006525,
     $     0.00631, 0.006101, 0.005979, 0.005878, 0.00582, 0.005755,
     $     0.005703, 0.005684, 0.005667, 0.00564, 0.005624, 0.005599,
     $     0.005573, 0.005567, 0.005533, 0.005495, 0.005457, 0.005407,
     $     0.005351, 0.005298, 0.005221, 0.005148, 0.005088, 0.00501,
     $     0.004929, 0.004847, 0.004768, 0.004694, 0.004621, 0.004541,
     $     0.004458, 0.004384, 0.004312, 0.004228, 0.00415, 0.004074,
     $     0.004003, 0.003932, 0.003856, 0.003789, 0.003726, 0.00366,
     $     0.003596, 0.003536, 0.00348, 0.003427, 0.003373, 0.003321,
     $     0.003279/


C     28058, JENDL evaluated data
      nbins(29) = 309
      data ( eg(29,i), i=1, 309) /6399520.0, 8172320.0, 1.22189D7,
     $     1.24D7, 1.26D7, 1.28D7, 1.3D7, 1.32D7, 1.34D7, 1.36D7, 1.38D7
     $     , 1.4D7, 1.41999D7, 1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.49102D7
     $     , 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7, 1.60977D7,
     $     1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.69393D7, 1.7D7, 1.72D7,
     $     1.72183D7, 1.73237D7, 1.74D7, 1.76D7, 1.76792D7, 1.78D7,
     $     1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7, 1.92D7, 1.94D7,
     $     1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7, 2.0668D7,
     $     2.08D7, 2.1D7, 2.11491D7, 2.12D7, 2.14D7, 2.16D7, 2.18D7,
     $     2.2D7, 2.22D7, 2.24D7, 2.2466D7, 2.26D7, 2.28D7, 2.3D7,
     $     2.31847D7, 2.32D7, 2.34D7, 2.36D7, 2.36082D7, 2.38D7, 2.4D7,
     $     2.41368D7, 2.42D7, 2.42765D7, 2.43411D7, 2.43835D7, 2.44D7,
     $     2.44991D7, 2.45477D7, 2.45948D7, 2.46D7, 2.467D7, 2.46825D7,
     $     2.48D7, 2.5D7, 2.52D7, 2.52706D7, 2.54D7, 2.56D7, 2.56946D7,
     $     2.58D7, 2.59418D7, 2.6D7, 2.62D7, 2.62452D7, 2.64D7, 2.66D7,
     $     2.67604D7, 2.68D7, 2.68871D7, 2.68937D7, 2.7D7, 2.72D7,
     $     2.74D7, 2.74346D7, 2.76D7, 2.77991D7, 2.78D7, 2.8D7, 2.82D7,
     $     2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.9488D7,
     $     2.96D7, 2.98D7, 2.99919D7, 3.0D7, 3.00143D7, 3.02D7, 3.04D7,
     $     3.06D7, 3.08D7, 3.1D7, 3.12D7, 3.14D7, 3.16D7, 3.18D7, 3.2D7,
     $     3.22D7, 3.24D7, 3.24506D7, 3.26D7, 3.28D7, 3.3D7, 3.32D7,
     $     3.34D7, 3.36D7, 3.37855D7, 3.38D7, 3.39265D7, 3.4D7, 3.42D7,
     $     3.44D7, 3.46D7, 3.48D7, 3.5D7, 3.52D7, 3.54D7, 3.56D7, 3.58D7
     $     , 3.58025D7, 3.58934D7, 3.6D7, 3.62D7, 3.62803D7, 3.64D7,
     $     3.6518D7, 3.66D7, 3.68D7, 3.7D7, 3.72D7, 3.74D7, 3.76D7,
     $     3.78D7, 3.8D7, 3.81455D7, 3.82D7, 3.84D7, 3.85477D7, 3.86D7,
     $     3.88D7, 3.9D7, 3.91071D7, 3.92D7, 3.94D7, 3.94516D7, 3.96D7,
     $     3.98D7, 4.0D7, 4.05786D7, 4.07758D7, 4.08519D7, 4.19609D7,
     $     4.2D7, 4.25212D7, 4.29384D7, 4.38912D7, 4.4D7, 4.48233D7,
     $     4.49116D7, 4.58446D7, 4.6D7, 4.6116D7, 4.75664D7, 4.8D7,
     $     4.87963D7, 4.91359D7, 5.0D7, 5.2D7, 5.28138D7, 5.33022D7,
     $     5.4D7, 5.44179D7, 5.44928D7, 5.45284D7, 5.51241D7, 5.54292D7,
     $     5.58285D7, 5.6D7, 5.73495D7, 5.8D7, 6.0D7, 6.00407D7,
     $     6.19155D7, 6.2D7, 6.30682D7, 6.31351D7, 6.4D7, 6.42762D7,
     $     6.53358D7, 6.55466D7, 6.6D7, 6.61008D7, 6.8D7, 6.8027D7,
     $     6.81708D7, 7.0D7, 7.1166D7, 7.12006D7, 7.2D7, 7.4D7,
     $     7.44695D7, 7.59796D7, 7.6D7, 7.62968D7, 7.64594D7, 7.70525D7,
     $     7.8D7, 7.87559D7, 7.89633D7, 7.92784D7, 8.0D7, 8.2D7,
     $     8.29881D7, 8.4D7, 8.5997D7, 8.6D7, 8.61416D7, 8.65704D7,
     $     8.8D7, 8.82959D7, 8.91357D7, 9.0D7, 9.0925D7, 9.2D7,
     $     9.23505D7, 9.4D7, 9.42884D7, 9.6D7, 9.7949D7, 9.8D7, 1.0D8,
     $     1.02D8, 1.03615D8, 1.03697D8, 1.04D8, 1.05283D8, 1.05969D8,
     $     1.06D8, 1.06202D8, 1.06884D8, 1.08D8, 1.08792D8, 1.1D8,
     $     1.12D8, 1.14D8, 1.14955D8, 1.16D8, 1.1634D8, 1.18D8,
     $     1.19061D8, 1.2D8, 1.21003D8, 1.22D8, 1.22875D8, 1.23199D8,
     $     1.23516D8, 1.24D8, 1.25359D8, 1.26D8, 1.28D8, 1.29257D8,
     $     1.3D8, 1.32D8, 1.34D8, 1.36D8, 1.36023D8, 1.3627D8, 1.37858D8
     $     , 1.38D8, 1.39676D8, 1.4D8/

      data ( gdr(29,i), i=1, 309) /2.62174D-4, 3.34802D-4, 5.00579D-4,
     $     5.08D-4, 0.001246, 0.002102, 0.003036, 0.004045, 0.005164,
     $     0.006445, 0.007867, 0.009331, 0.0107496, 0.01075, 0.01215,
     $     0.01353, 0.01491, 0.0156577, 0.01629, 0.01767, 0.01908,
     $     0.02057, 0.02217, 0.02382, 0.0246382, 0.02552, 0.02722,
     $     0.02894, 0.03062, 0.0317274, 0.03222, 0.0337, 0.03382,
     $     0.0345186, 0.03503, 0.03618, 0.0365389, 0.03709, 0.03778,
     $     0.0382, 0.03837, 0.03828, 0.03795, 0.03735, 0.03651, 0.03557,
     $     0.03472, 0.03403, 0.03344, 0.0329, 0.03239, 0.03185,
     $     0.0316579, 0.03129, 0.03072, 0.0303085, 0.03017, 0.02964,
     $     0.02913, 0.02861, 0.02809, 0.02758, 0.02706, 0.0268934,
     $     0.02656, 0.02608, 0.02568, 0.0253195, 0.02529, 0.02497,
     $     0.02468, 0.0246684, 0.0244, 0.02394, 0.0236652, 0.02354,
     $     0.0233819, 0.0232496, 0.0231633, 0.02313, 0.0229104,
     $     0.0228038, 0.0227012, 0.02269, 0.0225205, 0.0224902, 0.02221,
     $     0.02172, 0.02121, 0.021028, 0.0207, 0.0202, 0.0199614, 0.0197
     $     , 0.0193652, 0.01923, 0.01879, 0.0186917, 0.01836, 0.01796,
     $     0.0176222, 0.01754, 0.0173599, 0.0173463, 0.01713, 0.01671,
     $     0.01628, 0.0162028, 0.01584, 0.0154318, 0.01543, 0.01504,
     $     0.01467, 0.01432, 0.01397, 0.01363, 0.01329, 0.01295, 0.0126,
     $     0.0124534, 0.01227, 0.01195, 0.0116331, 0.01162, 0.0115997,
     $     0.01134, 0.01109, 0.01085, 0.01061, 0.01037, 0.01014,
     $     0.009912, 0.009675, 0.009453, 0.009234, 0.009015, 0.008798,
     $     0.00874249, 0.008581, 0.008364, 0.008151, 0.007936, 0.007724,
     $     0.007514, 0.00732454, 0.00731, 0.00718338, 0.007111, 0.006915
     $     , 0.006723, 0.006537, 0.006355, 0.00618, 0.006011, 0.005848,
     $     0.005689, 0.005537, 0.00553513, 0.00546682, 0.005388,
     $     0.005244, 0.00518641, 0.005102, 0.00501885, 0.004962,
     $     0.004824, 0.004687, 0.00455, 0.004413, 0.004276, 0.004139,
     $     0.004002, 0.00390401, 0.003868, 0.003737, 0.00364493,
     $     0.003613, 0.003498, 0.00339, 0.00333603, 0.00329, 0.003199,
     $     0.00317838, 0.00312, 0.003051, 0.002996, 0.0029397,
     $     0.00292093, 0.00291375, 0.00281246, 0.002809, 0.00279887,
     $     0.00279089, 0.00277301, 0.002771, 0.00278098, 0.00278205,
     $     0.00279317, 0.002795, 0.00279817, 0.0028374, 0.002849,
     $     0.00289244, 0.00291095, 0.002958, 0.003104, 0.00314516,
     $     0.00316983, 0.003205, 0.003222, 0.00322504, 0.00322649,
     $     0.00325064, 0.00326298, 0.00327909, 0.003286, 0.00337433,
     $     0.003417, 0.003582, 0.00358356, 0.00365481, 0.003658,
     $     0.00372261, 0.00372666, 0.003779, 0.00379243, 0.00384385,
     $     0.00385406, 0.003876, 0.00387925, 0.00394, 0.00394141,
     $     0.00394889, 0.004044, 0.00403931, 0.00403917, 0.004036,
     $     0.004056, 0.00406757, 0.0041045, 0.004105, 0.00410214,
     $     0.00410058, 0.00409493, 0.004086, 0.00405495, 0.00404653,
     $     0.0040338, 0.004005, 0.004004, 0.00401046, 0.004017,
     $     0.0039471, 0.003947, 0.00394115, 0.00392354, 0.003866,
     $     0.00386495, 0.003862, 0.003859, 0.00382956, 0.003796,
     $     0.00378119, 0.003713, 0.00370258, 0.003642, 0.0036069,
     $     0.003606, 0.003589, 0.003517, 0.00345213, 0.00344889,
     $     0.003437, 0.00340537, 0.00338874, 0.003388, 0.00338391,
     $     0.00337017, 0.003348, 0.003328, 0.003298, 0.003263, 0.00319,
     $     0.00315963, 0.003127, 0.00312201, 0.003098, 0.00308306,
     $     0.00307, 0.00305183, 0.003034, 0.00300481, 0.00299412,
     $     0.00298374, 0.002968, 0.00293452, 0.002919, 0.002886,
     $     0.00285249, 0.002833, 0.002799, 0.002749, 0.002713, 0.0027128
     $     , 0.00271068, 0.00269719, 0.002696, 0.00265484, 0.002647/


C     28060, JENDL evaluated data
      nbins(30) = 209
      data ( eg(30,i), i=1, 209) /1.1388D7, 1.14D7, 1.16D7, 1.18D7,
     $     1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7, 1.32D7, 1.34D7,
     $     1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7,
     $     1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7, 1.62D7, 1.64D7, 1.66D7
     $     , 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7,
     $     1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7
     $     , 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7, 2.08D7, 2.1D7,
     $     2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7, 2.22D7, 2.24D7, 2.26D7
     $     , 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7, 2.38D7, 2.4D7,
     $     2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7, 2.52D7, 2.54D7, 2.56D7
     $     , 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7, 2.68D7, 2.7D7,
     $     2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7, 2.82D7, 2.84D7, 2.86D7
     $     , 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7, 2.98D7, 3.0D7,
     $     3.05D7, 3.1D7, 3.15D7, 3.2D7, 3.25D7, 3.3D7, 3.35D7, 3.4D7,
     $     3.45D7, 3.5D7, 3.55D7, 3.6D7, 3.65D7, 3.7D7, 3.75D7, 3.8D7,
     $     3.85D7, 3.9D7, 3.95D7, 4.0D7, 4.05D7, 4.1D7, 4.15D7, 4.2D7,
     $     4.25D7, 4.3D7, 4.35D7, 4.4D7, 4.45D7, 4.5D7, 4.55D7, 4.6D7,
     $     4.65D7, 4.7D7, 4.75D7, 4.8D7, 4.85D7, 4.9D7, 4.95D7, 5.0D7,
     $     5.05D7, 5.1D7, 5.15D7, 5.2D7, 5.25D7, 5.3D7, 5.35D7, 5.4D7,
     $     5.45D7, 5.5D7, 5.55D7, 5.6D7, 5.65D7, 5.7D7, 5.75D7, 5.8D7,
     $     5.85D7, 5.9D7, 5.95D7, 6.0D7, 6.05D7, 6.1D7, 6.15D7, 6.2D7,
     $     6.25D7, 6.3D7, 6.35D7, 6.4D7, 6.45D7, 6.5D7, 6.55D7, 6.6D7,
     $     6.65D7, 6.7D7, 6.75D7, 6.8D7, 6.85D7, 6.9D7, 6.95D7, 7.0D7,
     $     7.05D7, 7.1D7, 7.15D7, 7.2D7, 7.25D7, 7.3D7, 7.35D7, 7.4D7,
     $     7.45D7, 7.5D7, 7.55D7, 7.6D7, 7.65D7, 7.7D7, 7.75D7, 7.8D7,
     $     7.85D7, 7.9D7, 7.95D7, 8.0D7, 8.5D7, 9.0D7, 9.5D7, 1.0D8,
     $     1.05D8, 1.1D8, 1.15D8, 1.2D8, 1.25D8, 1.3D8, 1.35D8, 1.4D8,
     $     1.45D8, 1.5D8/

      data ( gdr(30,i), i=1, 209) /0.0, 6.1441D-5, 3.2877D-4, 0.0010041,
     $     0.0025269, 0.0058965, 0.01337, 0.016077, 0.014742, 0.016729,
     $     0.020583, 0.024097, 0.025438, 0.025628, 0.025831, 0.026886,
     $     0.029144, 0.03275, 0.037803, 0.044384, 0.052515, 0.062025,
     $     0.071832, 0.081614, 0.09069, 0.098041, 0.10297, 0.1054,
     $     0.10585, 0.1052, 0.10431, 0.10383, 0.1037, 0.10419, 0.10546,
     $     0.1073, 0.10943, 0.11159, 0.1136, 0.1153, 0.11642, 0.11635,
     $     0.11408, 0.10922, 0.10247, 0.094978, 0.087753, 0.081419,
     $     0.076193, 0.072093, 0.069061, 0.067009, 0.065756, 0.065325,
     $     0.065604, 0.066353, 0.067404, 0.068445, 0.069186, 0.069321,
     $     0.068635, 0.067073, 0.064736, 0.061843, 0.058628, 0.055332,
     $     0.052154, 0.049192, 0.046524, 0.044206, 0.042247, 0.040639,
     $     0.039372, 0.038434, 0.037814, 0.037503, 0.037496, 0.037759,
     $     0.038325, 0.039189, 0.040323, 0.041714, 0.04332, 0.045068,
     $     0.046853, 0.048537, 0.049962, 0.05097, 0.051442, 0.051315,
     $     0.050612, 0.049433, 0.047913, 0.046206, 0.044473, 0.040749,
     $     0.038628, 0.038289, 0.039184, 0.039974, 0.038985, 0.035634,
     $     0.031002, 0.026424, 0.02256, 0.019516, 0.017165, 0.015341,
     $     0.013906, 0.012761, 0.011837, 0.011079, 0.01045, 0.0099168,
     $     0.0094631, 0.0090768, 0.0087434, 0.0084529, 0.0081965,
     $     0.0079701, 0.0077683, 0.0075905, 0.0074316, 0.0072879,
     $     0.0071583, 0.0070415, 0.006936, 0.0068402, 0.0067523,
     $     0.0066719, 0.0065986, 0.0065316, 0.0064701, 0.0064131,
     $     0.0063598, 0.0063119, 0.0062678, 0.0062268, 0.0061878,
     $     0.0061517, 0.0061184, 0.0060874, 0.0060587, 0.0060319,
     $     0.006007, 0.0059838, 0.0059621, 0.0059417, 0.0059226,
     $     0.0059047, 0.0058878, 0.0058718, 0.0058567, 0.0058421,
     $     0.0058277, 0.0058137, 0.0058003, 0.0057873, 0.005774,
     $     0.0057611, 0.0057487, 0.0057367, 0.005725, 0.0057136,
     $     0.0057024, 0.0056915, 0.0056807, 0.0056702, 0.0056597,
     $     0.0056494, 0.0056392, 0.0056291, 0.005619, 0.005609,
     $     0.0055984, 0.0055878, 0.0055771, 0.0055664, 0.0055553,
     $     0.0055442, 0.0055331, 0.005522, 0.0055108, 0.0054995,
     $     0.0054882, 0.0054768, 0.0054654, 0.0054538, 0.0054422,
     $     0.0054305, 0.0054187, 0.0054068, 0.0053948, 0.0053826,
     $     0.0053701, 0.0052373, 0.0050932, 0.0049777, 0.0048463,
     $     0.0047033, 0.0045537, 0.0044002, 0.0042485, 0.0041021,
     $     0.0039645, 0.0038384, 0.0037277, 0.0036353, 0.003564/


C     29063, JENDL evaluated data
      nbins(31) = 328
      data ( eg(31,i), i=1, 328) /5777470.0, 6122450.0, 1.0853D7,
     $     1.1D7, 1.12D7, 1.14D7, 1.16D7, 1.18D7, 1.2D7, 1.22D7, 1.24D7
     $     ,1.26D7, 1.28D7, 1.28114D7, 1.3D7, 1.32D7, 1.33787D7, 1.34D7
     $     ,1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.44951D7, 1.46D7
     $     ,1.48D7, 1.5D7, 1.51062D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7
     $     ,1.6D7, 1.60579D7, 1.62D7, 1.64D7, 1.65075D7, 1.66D7, 1.68D7
     $     ,1.7D7, 1.7121D7, 1.71906D7, 1.72D7, 1.72591D7, 1.74D7,
     $     1.76D7, 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7,
     $     1.88314D7,1.88631D7, 1.9D7, 1.92D7, 1.94D7, 1.94836D7,
     $     1.94868D7,1.94958D7, 1.96D7, 1.97392D7, 1.98D7, 1.98815D7,
     $     1.98828D7,2.0D7, 2.02D7, 2.04D7, 2.06D7, 2.06399D7, 2.08D7,
     $     2.1D7,2.10635D7, 2.12D7, 2.13725D7, 2.14D7, 2.16D7, 2.18D7,
     $     2.2D7,2.22D7, 2.23463D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7,
     $     2.30085D7,2.32D7, 2.34D7, 2.34991D7, 2.36D7, 2.37248D7,
     $     2.38D7,2.38306D7, 2.39387D7, 2.4D7, 2.42D7, 2.44D7, 2.44401D7
     $     ,2.46D7, 2.47893D7, 2.48D7, 2.5D7, 2.52D7, 2.54D7, 2.56D7
     $     ,2.58D7, 2.6D7, 2.60362D7, 2.62D7, 2.63279D7, 2.64D7
     $     ,2.64658D7, 2.66D7, 2.68D7, 2.7D7, 2.70303D7, 2.72D7, 2.74D7
     $     ,2.75003D7, 2.76D7, 2.78D7, 2.8D7, 2.82D7, 2.84D7, 2.86D7
     $     ,2.88D7, 2.9D7, 2.90474D7, 2.92D7, 2.93216D7, 2.94D7
     $     ,2.94301D7, 2.96D7, 2.98D7, 3.0D7, 3.02D7, 3.02379D7
     $     ,3.03528D7, 3.04D7, 3.06D7, 3.08D7, 3.09746D7, 3.1D7
     $     ,3.10939D7, 3.12D7, 3.14D7, 3.14489D7, 3.16D7, 3.18D7
     $     ,3.18839D7, 3.2D7, 3.22D7, 3.24D7, 3.26D7, 3.28D7, 3.3D7
     $     ,3.31855D7, 3.31874D7, 3.32D7, 3.33939D7, 3.34D7, 3.36D7
     $     ,3.38D7, 3.4D7, 3.42D7, 3.42817D7, 3.44D7, 3.46D7, 3.48D7
     $     ,3.5D7, 3.5049D7, 3.52D7, 3.54D7, 3.56D7, 3.56514D7, 3.58D7
     $     ,3.6D7, 3.62D7, 3.64D7, 3.66D7, 3.68D7, 3.7D7, 3.72D7
     $     ,3.72664D7, 3.74D7, 3.76D7, 3.78D7, 3.8D7, 3.81554D7, 3.82D7
     $     ,3.84D7, 3.86D7, 3.88D7, 3.9D7, 3.92D7, 3.92589D7, 3.94D7
     $     ,3.96D7, 3.96864D7, 3.98D7, 4.0D7, 4.10727D7, 4.15099D7
     $     ,4.1517D7, 4.2D7, 4.21871D7, 4.24156D7, 4.31504D7, 4.38883D7
     $     ,4.39145D7, 4.4D7, 4.41494D7, 4.46031D7, 4.6D7, 4.62183D7
     $     ,4.76699D7, 4.76851D7, 4.8D7, 4.95388D7, 5.0D7, 5.00902D7
     $     ,5.13644D7, 5.15229D7, 5.17278D7, 5.2D7, 5.4D7, 5.42728D7
     $     ,5.44052D7, 5.6D7, 5.66131D7, 5.8D7, 5.80223D7, 5.80429D7
     $     ,5.83892D7, 5.84506D7, 5.86266D7, 5.8926D7, 5.96662D7, 6.0D7
     $     ,6.17157D7, 6.2D7, 6.4D7, 6.49336D7, 6.56872D7, 6.6D7
     $     ,6.66988D7, 6.7848D7, 6.8D7, 6.81055D7, 6.92763D7, 7.0D7
     $     ,7.03826D7, 7.1747D7, 7.2D7, 7.34134D7, 7.4D7, 7.41765D7
     $     ,7.45635D7, 7.54549D7, 7.6D7, 7.74669D7, 7.8D7, 7.98215D7
     $     ,8.0D7, 8.19441D7, 8.2D7, 8.31875D7, 8.34744D7, 8.4D7
     $     ,8.51122D7, 8.6D7, 8.72842D7, 8.75795D7, 8.77832D7, 8.8D7
     $     ,9.0D7, 9.09506D7, 9.2D7, 9.35489D7, 9.4D7, 9.42268D7
     $     ,9.4461D7, 9.6D7, 9.71207D7, 9.8D7, 9.84862D7, 9.85227D7
     $     ,1.0D8, 1.00299D8, 1.01082D8, 1.02D8, 1.03367D8, 1.03581D8
     $     ,1.04D8, 1.06D8, 1.08D8, 1.1D8, 1.11405D8, 1.12D8, 1.1223D8
     $     ,1.13347D8, 1.13848D8, 1.14D8, 1.15793D8, 1.16D8, 1.16524D8
     $     ,1.17917D8, 1.18D8, 1.2D8, 1.21508D8, 1.22D8, 1.24D8, 1.26D8
     $     ,1.26333D8, 1.26752D8, 1.2763D8, 1.28D8, 1.3D8, 1.31486D8
     $     ,1.31827D8, 1.32D8, 1.32829D8, 1.34D8, 1.35597D8, 1.36D8
     $     ,1.37722D8, 1.38D8, 1.4D8/

      data ( gdr(31,i), i=1, 328) /0.00106831, 0.0011321, 0.00200681,
     $     0.002034, 0.004388, 0.007615, 0.009665, 0.01078, 0.01212,
     $     0.01362, 0.01523, 0.01691, 0.01864, 0.0187398, 0.02045,
     $     0.02236, 0.0241842, 0.02441, 0.02661, 0.02901, 0.03157,
     $     0.03521, 0.03931, 0.041424, 0.04387, 0.04891, 0.05442,
     $     0.0574717, 0.06029, 0.06641, 0.07258, 0.0785, 0.08389,
     $     0.0851902, 0.08845, 0.09199, 0.093271, 0.09438, 0.09572,
     $     0.09615, 0.0960105, 0.0959307, 0.09592, 0.0957238, 0.09526,
     $     0.09448, 0.09355, 0.09248, 0.0912, 0.08963, 0.08763, 0.08515,
     $     0.0846834, 0.0842136, 0.08223, 0.07893, 0.0754, 0.0738642,
     $     0.0738075, 0.0736446, 0.07179, 0.0692887, 0.06823, 0.0668481,
     $     0.0668265, 0.0649, 0.06189, 0.05901, 0.05661, 0.0562142,
     $     0.05466, 0.05331, 0.0530247, 0.05242, 0.0521005, 0.05205,
     $     0.05225, 0.05297, 0.05403, 0.05518, 0.055845, 0.05609,
     $     0.05637, 0.05569, 0.05395, 0.0538325, 0.05127, 0.04799,
     $     0.0462135, 0.04448, 0.0423113, 0.04106, 0.0405677, 0.0388811,
     $     0.03796, 0.03531, 0.03317, 0.0328501, 0.03161, 0.0306815,
     $     0.03063, 0.03025, 0.03035, 0.03058, 0.03038, 0.02918, 0.02701
     $     , 0.0265267, 0.02445, 0.0228452, 0.02199, 0.0212601, 0.01985,
     $     0.0181, 0.01666, 0.0164731, 0.01547, 0.0145, 0.0140769,
     $     0.01367, 0.01296, 0.01235, 0.0118, 0.01131, 0.01088, 0.01049,
     $     0.01015, 0.0100709, 0.009821, 0.00964396, 0.009532,
     $     0.00949179, 0.009269, 0.009025, 0.008801, 0.008619,
     $     0.00858275, 0.00847409, 0.00843, 0.008254, 0.00809,
     $     0.00795706, 0.007938, 0.00787089, 0.007796, 0.007662,
     $     0.00763116, 0.007537, 0.00742, 0.00737314, 0.007309, 0.007206
     $     , 0.007108, 0.007016, 0.006929, 0.006847, 0.00677463,
     $     0.00677388, 0.006769, 0.00669919, 0.006697, 0.006627,
     $     0.006563, 0.006501, 0.006442, 0.00641945, 0.006387, 0.006333,
     $     0.006284, 0.006235, 0.00622392, 0.00619, 0.006146, 0.006105,
     $     0.00609467, 0.006065, 0.006028, 0.005992, 0.005957, 0.005924,
     $     0.005893, 0.005862, 0.005833, 0.00582368, 0.005805, 0.005779,
     $     0.005753, 0.005729, 0.00571266, 0.005708, 0.005688, 0.005669,
     $     0.005652, 0.005634, 0.005618, 0.00561298, 0.005601, 0.005587,
     $     0.00558051, 0.005572, 0.005559, 0.0055018, 0.00547908,
     $     0.00547871, 0.005454, 0.00544698, 0.00543847, 0.00541148,
     $     0.00538496, 0.00538403, 0.005381, 0.00537825, 0.00536995,
     $     0.005345, 0.00534433, 0.00533997, 0.00533993, 0.005339,
     $     0.00535368, 0.005358, 0.00535942, 0.00537926, 0.0053817,
     $     0.00538484, 0.005389, 0.005434, 0.00544104, 0.00544444,
     $     0.005485, 0.00550047, 0.005535, 0.00553551, 0.00553598,
     $     0.00554385, 0.00554524, 0.00554922, 0.00555598, 0.00557257,
     $     0.00558, 0.00560664, 0.005611, 0.005627, 0.00562888,
     $     0.00563038, 0.005631, 0.00562676, 0.0056199, 0.005619,
     $     0.00561777, 0.00560424, 0.005596, 0.00558941, 0.00556624,
     $     0.005562, 0.00552861, 0.005515, 0.00551016, 0.00549958,
     $     0.00547552, 0.005461, 0.00541604, 0.0054, 0.0053407, 0.005335
     $     , 0.00526498, 0.005263, 0.00521932, 0.00520892, 0.00519,
     $     0.00514794, 0.005115, 0.00506587, 0.00505474, 0.0050471,
     $     0.005039, 0.004961, 0.00492452, 0.004885, 0.00482356,
     $     0.004806, 0.00479712, 0.004788, 0.004729, 0.00468437, 0.00465
     $     , 0.00463102, 0.0046296, 0.004573, 0.00456129, 0.00453098,
     $     0.004496, 0.00444306, 0.00443489, 0.004419, 0.004343,
     $     0.004268, 0.004194, 0.00414102, 0.004119, 0.0041107,
     $     0.00407091, 0.0040533, 0.004048, 0.00398243, 0.003975,
     $     0.00395642, 0.00390785, 0.003905, 0.003836, 0.00378603,
     $     0.00377, 0.003704, 0.003641, 0.0036307, 0.00361783,
     $     0.00359114, 0.00358, 0.003521, 0.00347996, 0.0034707,
     $     0.003466, 0.00344256, 0.00341, 0.00336995, 0.00336,
     $     0.00331687, 0.00331, 0.003265/


C     29065, JENDL evaluated data
      nbins(32) = 341
      data ( eg(32,i), i=1, 341) /6789590.0, 7452790.0, 9910240.0,
     $     1.002D7, 1.02D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7, 1.12D7,
     $     1.14D7, 1.16D7, 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7
     $     , 1.3D7, 1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7,
     $     1.44D7, 1.46D7, 1.47165D7, 1.48D7, 1.48862D7, 1.5D7, 1.52D7,
     $     1.54D7, 1.54151D7, 1.54668D7, 1.56D7, 1.58D7, 1.59353D7,
     $     1.59704D7, 1.6D7, 1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7,
     $     1.70017D7, 1.72D7, 1.74D7, 1.75317D7, 1.76D7, 1.78D7,
     $     1.78262D7, 1.8D7, 1.80821D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7,
     $     1.89473D7, 1.89994D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7, 1.98D7,
     $     2.0D7, 2.00006D7, 2.00185D7, 2.02D7, 2.03963D7, 2.04D7,
     $     2.06D7, 2.06876D7, 2.07628D7, 2.08D7, 2.1D7, 2.12D7, 2.14D7,
     $     2.16D7, 2.17018D7, 2.18D7, 2.19579D7, 2.2D7, 2.2146D7,
     $     2.21599D7, 2.22D7, 2.24D7, 2.26D7, 2.26302D7, 2.28D7, 2.3D7,
     $     2.31412D7, 2.32D7, 2.34D7, 2.36D7, 2.38D7, 2.4D7, 2.42D7,
     $     2.42072D7, 2.43851D7, 2.43868D7, 2.44D7, 2.46D7, 2.48D7,
     $     2.5D7, 2.52D7, 2.5348D7, 2.54D7, 2.55776D7, 2.56D7, 2.58D7,
     $     2.6D7, 2.62D7, 2.64D7, 2.66D7, 2.66435D7, 2.68D7, 2.7D7,
     $     2.7127D7, 2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.79153D7, 2.8D7,
     $     2.80235D7, 2.80302D7, 2.82D7, 2.84D7, 2.84091D7, 2.86D7,
     $     2.86792D7, 2.88D7, 2.88346D7, 2.89707D7, 2.9D7, 2.92D7,
     $     2.94D7, 2.96D7, 2.98D7, 2.99873D7, 3.0D7, 3.01611D7, 3.02D7,
     $     3.02287D7, 3.04D7, 3.06D7, 3.08D7, 3.1D7, 3.12D7, 3.14D7,
     $     3.16D7, 3.16126D7, 3.18D7, 3.18824D7, 3.2D7, 3.20474D7,
     $     3.22D7, 3.24D7, 3.26D7, 3.28D7, 3.29977D7, 3.3D7, 3.32D7,
     $     3.34D7, 3.35286D7, 3.36D7, 3.36624D7, 3.38D7, 3.3811D7, 3.4D7
     $     , 3.42D7, 3.44D7, 3.46D7, 3.48D7, 3.5D7, 3.52D7, 3.54D7,
     $     3.54451D7, 3.56D7, 3.58D7, 3.6D7, 3.62D7, 3.64D7, 3.66D7,
     $     3.68D7, 3.7D7, 3.72D7, 3.74D7, 3.75056D7, 3.75654D7, 3.76D7,
     $     3.76525D7, 3.77732D7, 3.77877D7, 3.78D7, 3.8D7, 3.82D7,
     $     3.84D7, 3.86D7, 3.88D7, 3.89948D7, 3.9D7, 3.92D7, 3.94D7,
     $     3.94787D7, 3.9518D7, 3.96D7, 3.96256D7, 3.98D7, 4.0D7,
     $     4.03565D7, 4.2D7, 4.28349D7, 4.31968D7, 4.38216D7, 4.4D7,
     $     4.46804D7, 4.59722D7, 4.6D7, 4.64559D7, 4.79724D7, 4.8D7,
     $     4.92751D7, 4.97552D7, 4.9934D7, 4.99902D7, 5.0D7, 5.00193D7,
     $     5.01798D7, 5.0217D7, 5.18072D7, 5.2D7, 5.20431D7, 5.23413D7,
     $     5.4D7, 5.49433D7, 5.6D7, 5.68784D7, 5.70758D7, 5.7592D7,
     $     5.79774D7, 5.8D7, 5.89052D7, 5.9016D7, 5.93361D7, 6.0D7,
     $     6.13396D7, 6.2D7, 6.2046D7, 6.30192D7, 6.4D7, 6.51517D7,
     $     6.51533D7, 6.6D7, 6.64144D7, 6.8D7, 6.81607D7, 7.0D7,
     $     7.00203D7, 7.2D7, 7.2099D7, 7.29462D7, 7.30817D7, 7.38773D7,
     $     7.4D7, 7.41254D7, 7.53433D7, 7.6D7, 7.75921D7, 7.8D7,
     $     7.87765D7, 7.98738D7, 8.0D7, 8.10119D7, 8.17191D7, 8.2D7,
     $     8.4D7, 8.4525D7, 8.47149D7, 8.50294D7, 8.6D7, 8.66354D7,
     $     8.67145D7, 8.8D7, 8.80723D7, 8.82557D7, 8.93737D7, 9.0D7,
     $     9.2D7, 9.4D7, 9.6D7, 9.68921D7, 9.71528D7, 9.8D7, 9.82785D7,
     $     9.85867D7, 9.96653D7, 1.0D8, 1.01301D8, 1.01819D8, 1.02D8,
     $     1.03158D8, 1.04D8, 1.04184D8, 1.06D8, 1.06564D8, 1.08D8,
     $     1.09182D8, 1.1D8, 1.12D8, 1.12485D8, 1.1295D8, 1.14D8,
     $     1.15219D8, 1.16D8, 1.16234D8, 1.16349D8, 1.18D8, 1.2D8,
     $     1.20132D8, 1.22D8, 1.22684D8, 1.2342D8, 1.24D8, 1.24591D8,
     $     1.26D8, 1.26475D8, 1.28D8, 1.3D8, 1.31422D8, 1.3149D8, 1.32D8
     $     , 1.34D8, 1.3435D8, 1.34684D8, 1.35497D8, 1.36D8, 1.38D8,
     $     1.38915D8, 1.4D8/

      data ( gdr(32,i), i=1, 341) /5.42083D-4, 5.95033D-4, 7.91237D-4,
     $     8.0D-4, 0.002864, 0.005023, 0.006855, 0.008368, 0.01, 0.01205
     $     , 0.01421, 0.0159, 0.01714, 0.01837, 0.02004, 0.022, 0.02398,
     $     0.02597, 0.02825, 0.03106, 0.03408, 0.03683, 0.03945, 0.04263
     $     , 0.0469, 0.05179, 0.05663, 0.0591834, 0.06107, 0.0628158,
     $     0.06518, 0.06903, 0.07263, 0.0728787, 0.0737341, 0.07597,
     $     0.07901, 0.0807701, 0.081231, 0.08162, 0.08373, 0.08528,
     $     0.08626, 0.08664, 0.08644, 0.0864335, 0.0857, 0.08447,
     $     0.0833712, 0.08281, 0.0808, 0.0804947, 0.07851, 0.0774703,
     $     0.07601, 0.07341, 0.07076, 0.06806, 0.0660504, 0.0653581,
     $     0.06535, 0.06269, 0.06004, 0.05735, 0.05474, 0.05255,
     $     0.0525448, 0.0524008, 0.05097, 0.0496545, 0.04963, 0.04803,
     $     0.047168, 0.0464431, 0.04609, 0.04389, 0.04182, 0.03988,
     $     0.038, 0.037066, 0.03619, 0.0348238, 0.03447, 0.0332929,
     $     0.0331835, 0.03287, 0.03139, 0.03001, 0.0298106, 0.02872,
     $     0.02753, 0.0267409, 0.02642, 0.0254, 0.02443, 0.02353,
     $     0.02271, 0.02194, 0.0219132, 0.0212633, 0.0212572, 0.02121,
     $     0.02053, 0.01988, 0.01928, 0.0187, 0.018306, 0.01817,
     $     0.0177162, 0.01766, 0.01718, 0.01673, 0.01629, 0.01589,
     $     0.0155, 0.0154207, 0.01514, 0.01478, 0.0145694, 0.01445,
     $     0.01414, 0.01384, 0.01354, 0.0133776, 0.01326, 0.0132279,
     $     0.0132187, 0.01299, 0.01272, 0.0127089, 0.01248, 0.0123882,
     $     0.01225, 0.0122097, 0.0120533, 0.01202, 0.0118, 0.01159,
     $     0.01139, 0.0112, 0.0110407, 0.01103, 0.0108928, 0.01086,
     $     0.0108353, 0.01069, 0.01052, 0.01037, 0.01023, 0.01009,
     $     0.009954, 0.009828, 0.00982062, 0.009712, 0.00966477,
     $     0.009598, 0.00957129, 0.009486, 0.009377, 0.009271, 0.009171,
     $     0.00907513, 0.009074, 0.008979, 0.008889, 0.00883228,
     $     0.008801, 0.00877528, 0.008719, 0.00871457, 0.008639,
     $     0.008563, 0.008488, 0.008413, 0.008339, 0.008271, 0.008208,
     $     0.00815, 0.00813617, 0.008089, 0.008027, 0.007967, 0.007912,
     $     0.00786, 0.00781, 0.007761, 0.007713, 0.007666, 0.007621,
     $     0.00759664, 0.00758291, 0.007575, 0.00756287, 0.00753513,
     $     0.00753181, 0.007529, 0.007486, 0.007447, 0.007411, 0.007376,
     $     0.007342, 0.00730888, 0.007308, 0.007271, 0.007238,
     $     0.00723011, 0.00722618, 0.007218, 0.0072171, 0.007211,
     $     0.007208, 0.00716823, 0.006992, 0.0068959, 0.00685523,
     $     0.00678635, 0.006767, 0.00671869, 0.00662987, 0.006628,
     $     0.00659741, 0.00649876, 0.006497, 0.00643066, 0.0064063,
     $     0.00639731, 0.00639449, 0.006394, 0.00639309, 0.00638552,
     $     0.00638377, 0.00631066, 0.006302, 0.00630023, 0.00628806,
     $     0.006222, 0.00619336, 0.006162, 0.00613933, 0.00613429,
     $     0.00612123, 0.00611156, 0.006111, 0.00609178, 0.00608945,
     $     0.00608275, 0.006069, 0.00604001, 0.006026, 0.00602504,
     $     0.00600491, 0.005985, 0.00596062, 0.00596059, 0.005943,
     $     0.00593354, 0.005898, 0.00589416, 0.005851, 0.00585048,
     $     0.005801, 0.00579823, 0.00577473, 0.005771, 0.00574932,
     $     0.005746, 0.00574204, 0.0057041, 0.005684, 0.00563446,
     $     0.005622, 0.00559766, 0.00556385, 0.00556, 0.00552477,
     $     0.00550054, 0.005491, 0.005423, 0.00540384, 0.00539696,
     $     0.00538561, 0.005351, 0.00532719, 0.00532425, 0.005277,
     $     0.00527442, 0.00526791, 0.00522865, 0.005207, 0.005134,
     $     0.005061, 0.004986, 0.00495176, 0.00494186, 0.00491,
     $     0.0048994, 0.00488772, 0.00484737, 0.004835, 0.00478461,
     $     0.00476487, 0.004758, 0.0047125, 0.00468, 0.00467297,
     $     0.004605, 0.00458328, 0.004529, 0.00448437, 0.004454,
     $     0.004379, 0.00436083, 0.00434354, 0.004305, 0.0042596,
     $     0.004231, 0.00422244, 0.00421826, 0.004159, 0.004088,
     $     0.00408338, 0.004019, 0.00399514, 0.00396976, 0.00395,
     $     0.00393088, 0.003886, 0.00387037, 0.003821, 0.003759,
     $     0.00371759, 0.00371562, 0.003701, 0.003644, 0.00363443,
     $     0.00362535, 0.00360343, 0.00359, 0.00354, 0.00351788,
     $     0.003492/


C     30064, JENDL evaluated data
      nbins(33) = 319
      data ( eg(33,i), i=1, 319) /3956340.0, 7712340.0, 1.03402D7,
     $     1.05867D7, 1.18616D7, 1.196D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7,
     $     1.27755D7, 1.28D7, 1.3D7, 1.32D7, 1.34D7, 1.36D7, 1.38D7,
     $     1.38348D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7,
     $     1.52D7, 1.54D7, 1.54548D7, 1.56D7, 1.58D7, 1.6D7, 1.62D7,
     $     1.62346D7, 1.63407D7, 1.64D7, 1.66D7, 1.6714D7, 1.68D7, 1.7D7
     $     , 1.72D7, 1.74D7, 1.76D7, 1.76798D7, 1.78D7, 1.8D7, 1.81448D7
     $     , 1.82D7, 1.84D7, 1.86D7, 1.8765D7, 1.88D7, 1.89698D7, 1.9D7,
     $     1.91475D7, 1.92D7, 1.94D7, 1.96D7, 1.96454D7, 1.98D7,
     $     1.98721D7, 1.99839D7, 2.0D7, 2.02D7, 2.02444D7, 2.02588D7,
     $     2.02969D7, 2.04D7, 2.06D7, 2.07088D7, 2.08D7, 2.09748D7,
     $     2.1D7, 2.10921D7, 2.12D7, 2.13694D7, 2.14D7, 2.15675D7,
     $     2.15933D7, 2.16D7, 2.18D7, 2.2D7, 2.22D7, 2.22252D7,
     $     2.23577D7, 2.24D7, 2.2499D7, 2.25801D7, 2.26D7, 2.28D7, 2.3D7
     $     , 2.32D7, 2.33706D7, 2.34D7, 2.36D7, 2.38D7, 2.4D7, 2.42D7,
     $     2.44D7, 2.46D7, 2.47698D7, 2.48D7, 2.49714D7, 2.5D7, 2.52D7,
     $     2.54D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7, 2.68D7
     $     , 2.7D7, 2.72D7, 2.72083D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7,
     $     2.82D7, 2.84D7, 2.86D7, 2.86673D7, 2.88D7, 2.9D7, 2.92D7,
     $     2.94D7, 2.96D7, 2.96755D7, 2.98D7, 3.0D7, 3.02D7, 3.03262D7,
     $     3.04D7, 3.05472D7, 3.06D7, 3.08D7, 3.09107D7, 3.1D7, 3.12D7,
     $     3.14D7, 3.16D7, 3.18D7, 3.18664D7, 3.19002D7, 3.2D7, 3.2015D7
     $     , 3.22D7, 3.22481D7, 3.24D7, 3.26D7, 3.28D7, 3.3D7, 3.32D7,
     $     3.34D7, 3.35861D7, 3.36D7, 3.37485D7, 3.38D7, 3.38711D7,
     $     3.4D7, 3.42D7, 3.43248D7, 3.44D7, 3.46D7, 3.48D7, 3.5D7,
     $     3.52D7, 3.52322D7, 3.54D7, 3.56D7, 3.58D7, 3.6D7, 3.60343D7,
     $     3.62D7, 3.63825D7, 3.64D7, 3.66D7, 3.68D7, 3.68783D7,
     $     3.69297D7, 3.7D7, 3.72D7, 3.74D7, 3.76D7, 3.78D7, 3.8D7,
     $     3.82D7, 3.84D7, 3.86D7, 3.88D7, 3.9D7, 3.91154D7, 3.92D7,
     $     3.93834D7, 3.94D7, 3.96D7, 3.98D7, 3.98099D7, 4.0D7, 4.2D7,
     $     4.20427D7, 4.27468D7, 4.33371D7, 4.36979D7, 4.4D7, 4.41017D7,
     $     4.4741D7, 4.5805D7, 4.6D7, 4.6482D7, 4.69713D7, 4.70418D7,
     $     4.8D7, 4.85833D7, 4.99864D7, 5.0D7, 5.0163D7, 5.10444D7,
     $     5.2D7, 5.24838D7, 5.4D7, 5.58118D7, 5.58401D7, 5.6D7,
     $     5.62033D7, 5.73559D7, 5.78026D7, 5.8D7, 5.90987D7, 5.94407D7,
     $     6.0D7, 6.14526D7, 6.1488D7, 6.2D7, 6.392D7, 6.4D7, 6.41562D7,
     $     6.49029D7, 6.6D7, 6.63166D7, 6.73924D7, 6.8D7, 6.84462D7,
     $     6.87334D7, 7.0D7, 7.14104D7, 7.2D7, 7.21344D7, 7.2646D7,
     $     7.32074D7, 7.4D7, 7.40842D7, 7.56034D7, 7.6D7, 7.61238D7,
     $     7.8D7, 7.94209D7, 8.0D7, 8.02066D7, 8.2D7, 8.21795D7, 8.4D7,
     $     8.51792D7, 8.6D7, 8.74528D7, 8.8D7, 8.98128D7, 8.99114D7,
     $     9.0D7, 9.00246D7, 9.02163D7, 9.03397D7, 9.2D7, 9.21226D7,
     $     9.4D7, 9.43755D7, 9.6D7, 9.8D7, 9.85715D7, 1.0D8, 1.02D8,
     $     1.02081D8, 1.03493D8, 1.04D8, 1.04842D8, 1.04878D8, 1.05952D8
     $     , 1.06D8, 1.08D8, 1.08876D8, 1.09088D8, 1.09824D8, 1.1D8,
     $     1.10265D8, 1.12D8, 1.14D8, 1.16D8, 1.18D8, 1.1839D8, 1.2D8,
     $     1.20041D8, 1.2025D8, 1.22D8, 1.2337D8, 1.23718D8, 1.24D8,
     $     1.24573D8, 1.26D8, 1.27678D8, 1.28D8, 1.3D8, 1.32D8, 1.34D8,
     $     1.3486D8, 1.35932D8, 1.36D8, 1.36585D8, 1.38D8, 1.4D8/

      data ( gdr(33,i), i=1, 319) /6.66889D-4, 0.00130001, 0.00174296,
     $     0.00178452, 0.00199941, 0.002016, 0.003115, 0.006895,
     $     0.009928, 0.0123, 0.0138766, 0.01411, 0.01544, 0.01657,
     $     0.01755, 0.01895, 0.02102, 0.021477, 0.02377, 0.02718,
     $     0.03125, 0.036, 0.04143, 0.04751, 0.05407, 0.06083, 0.0626002
     $     , 0.06751, 0.07373, 0.07921, 0.08359, 0.0840996, 0.0856729,
     $     0.08656, 0.08818, 0.0885109, 0.08876, 0.08871, 0.08849,
     $     0.08838, 0.08859, 0.0888377, 0.08921, 0.09018, 0.0910268,
     $     0.09135, 0.09253, 0.09346, 0.0938314, 0.09391, 0.0937486,
     $     0.09372, 0.0930545, 0.09282, 0.09123, 0.08906, 0.0884652,
     $     0.08648, 0.0854529, 0.0838919, 0.08367, 0.08048, 0.0797944,
     $     0.079574, 0.078994, 0.07745, 0.07472, 0.0734289, 0.07237,
     $     0.0706615, 0.07042, 0.0697419, 0.06896, 0.0682133, 0.06808,
     $     0.0677697, 0.0677223, 0.06771, 0.06777, 0.0682, 0.06885,
     $     0.0689455, 0.0694491, 0.06961, 0.0699812, 0.0702855, 0.07036,
     $     0.07088, 0.07097, 0.07059, 0.0698295, 0.0697, 0.06834,
     $     0.06659, 0.06452, 0.06228, 0.05995, 0.05765, 0.0556733,
     $     0.05533, 0.0535244, 0.05323, 0.05177, 0.05047, 0.04938,
     $     0.04853, 0.04792, 0.04754, 0.04732, 0.04719, 0.04697, 0.04663
     $     , 0.04598, 0.0459375, 0.04497, 0.04352, 0.04171, 0.0396,
     $     0.03732, 0.03499, 0.03265, 0.0318979, 0.03047, 0.02841,
     $     0.02653, 0.02483, 0.02331, 0.0227742, 0.02192, 0.02069,
     $     0.0193, 0.0186557, 0.01829, 0.017615, 0.01738, 0.01656,
     $     0.0161456, 0.01582, 0.01515, 0.01453, 0.01398, 0.01346,
     $     0.0133017, 0.0132221, 0.01299, 0.0129572, 0.01256, 0.0124623,
     $     0.01216, 0.0118, 0.01145, 0.01114, 0.01084, 0.01057,
     $     0.0103278, 0.01031, 0.0101311, 0.01007, 0.00999075, 0.009849,
     $     0.00964, 0.00951585, 0.009442, 0.009257, 0.009081, 0.008915,
     $     0.008759, 0.00873476, 0.00861, 0.00847, 0.008338, 0.008211,
     $     0.00819023, 0.008091, 0.00798686, 0.007977, 0.007868,
     $     0.007765, 0.0077264, 0.00770124, 0.007667, 0.007573, 0.007482
     $     , 0.007397, 0.007316, 0.007237, 0.007169, 0.007107, 0.007046,
     $     0.00699, 0.006935, 0.00690316, 0.00688, 0.00683229, 0.006828,
     $     0.006777, 0.006726, 0.00672356, 0.006677, 0.006258,
     $     0.00625147, 0.00614573, 0.00605976, 0.00600837, 0.005966,
     $     0.00595544, 0.00589006, 0.00578482, 0.005766, 0.00573289,
     $     0.00569982, 0.0056951, 0.005632, 0.00560761, 0.00555054,
     $     0.00555, 0.00554608, 0.00552518, 0.005503, 0.0054971,
     $     0.005479, 0.00547991, 0.00547992, 0.00548, 0.00548155,
     $     0.00549022, 0.00549354, 0.005495, 0.00550772, 0.00551164,
     $     0.005518, 0.00553841, 0.00553891, 0.005546, 0.00556713,
     $     0.005568, 0.00556919, 0.00557482, 0.005583, 0.00558284,
     $     0.0055823, 0.005582, 0.00557884, 0.00557682, 0.005568,
     $     0.00554887, 0.005541, 0.00553827, 0.00552793, 0.00551669,
     $     0.005501, 0.00549873, 0.00545836, 0.005448, 0.0054439,
     $     0.005383, 0.00533085, 0.00531, 0.00530148, 0.005229,
     $     0.00522114, 0.005143, 0.00509008, 0.005054, 0.00498606,
     $     0.004961, 0.00487472, 0.00487012, 0.004866, 0.00486478,
     $     0.00485531, 0.00484923, 0.004769, 0.00476306, 0.004674,
     $     0.00465586, 0.004579, 0.004484, 0.00445703, 0.004391,
     $     0.004299, 0.0042952, 0.00422999, 0.004207, 0.00416994,
     $     0.00416839, 0.00412207, 0.00412, 0.004033, 0.00399534,
     $     0.00398635, 0.00395537, 0.003948, 0.0039367, 0.003864,
     $     0.003784, 0.003705, 0.003629, 0.00361414, 0.003554,
     $     0.00355253, 0.003545, 0.003483, 0.00343548, 0.00342356,
     $     0.003414, 0.00339457, 0.003347, 0.00329315, 0.003283,
     $     0.003222, 0.003163, 0.003107, 0.00308356, 0.00305481,
     $     0.003053, 0.00303821, 0.003003, 0.002956/


C     40090, JENDL evaluated data
      nbins(34) = 336
      data ( eg(34,i), i=1, 336) /6671290.0, 8354800.0, 1.06089D7,
     $     1.17423D7, 1.19698D7, 1.21D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7,
     $     1.3D7, 1.3121D7, 1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7,
     $     1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7, 1.52D7, 1.54D7,
     $     1.54003D7, 1.54262D7, 1.56D7, 1.58D7, 1.58369D7, 1.6D7,
     $     1.62D7, 1.64D7, 1.65363D7, 1.66D7, 1.66618D7, 1.68D7, 1.7D7,
     $     1.72D7, 1.72939D7, 1.74D7, 1.76D7, 1.76066D7, 1.78D7, 1.8D7,
     $     1.82D7, 1.82791D7, 1.84D7, 1.86D7, 1.88D7, 1.88208D7, 1.9D7,
     $     1.90515D7, 1.92D7, 1.93418D7, 1.94D7, 1.96D7, 1.97674D7,
     $     1.98D7, 2.0D7, 2.01972D7, 2.02D7, 2.04D7, 2.06D7, 2.0701D7,
     $     2.08D7, 2.1D7, 2.12D7, 2.12868D7, 2.13679D7, 2.14D7, 2.16D7,
     $     2.18D7, 2.2D7, 2.20045D7, 2.21966D7, 2.22D7, 2.24D7,
     $     2.24219D7, 2.26D7, 2.28D7, 2.3D7, 2.31041D7, 2.32D7,
     $     2.32093D7, 2.34D7, 2.34812D7, 2.35015D7, 2.35615D7, 2.35918D7
     $     , 2.36D7, 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.45556D7, 2.46D7,
     $     2.46029D7, 2.47114D7, 2.48D7, 2.49298D7, 2.5D7, 2.52D7,
     $     2.54D7, 2.56D7, 2.56825D7, 2.57178D7, 2.58D7, 2.6D7,
     $     2.60398D7, 2.62D7, 2.64D7, 2.64723D7, 2.66D7, 2.67071D7,
     $     2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7, 2.82D7,
     $     2.84D7, 2.86D7, 2.86849D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7,
     $     2.96D7, 2.9736D7, 2.98D7, 3.0D7, 3.02D7, 3.02494D7, 3.04D7,
     $     3.06D7, 3.08D7, 3.1D7, 3.10753D7, 3.12D7, 3.13408D7, 3.14D7,
     $     3.16D7, 3.18D7, 3.2D7, 3.20962D7, 3.22D7, 3.24D7, 3.26D7,
     $     3.28D7, 3.3D7, 3.32D7, 3.32898D7, 3.34D7, 3.35251D7, 3.36D7,
     $     3.36341D7, 3.36556D7, 3.38D7, 3.4D7, 3.4193D7, 3.42D7, 3.44D7
     $     , 3.44997D7, 3.46D7, 3.46579D7, 3.48D7, 3.5D7, 3.52D7, 3.54D7
     $     , 3.54142D7, 3.56D7, 3.58D7, 3.6D7, 3.62D7, 3.64D7, 3.66D7,
     $     3.67362D7, 3.68D7, 3.7D7, 3.72D7, 3.74D7, 3.76D7, 3.76267D7,
     $     3.78D7, 3.8D7, 3.80812D7, 3.82D7, 3.82654D7, 3.84D7, 3.86D7,
     $     3.88D7, 3.9D7, 3.92D7, 3.94D7, 3.96D7, 3.96846D7, 3.98D7,
     $     4.0D7, 4.00925D7, 4.06645D7, 4.07306D7, 4.12925D7, 4.15883D7,
     $     4.2D7, 4.2143D7, 4.32482D7, 4.4D7, 4.43572D7, 4.52833D7,
     $     4.57975D7, 4.6D7, 4.61211D7, 4.61941D7, 4.62686D7, 4.66022D7,
     $     4.77543D7, 4.8D7, 4.8596D7, 4.89537D7, 4.96454D7, 5.0D7,
     $     5.09192D7, 5.2D7, 5.2255D7, 5.4D7, 5.59699D7, 5.6D7, 5.6349D7
     $     , 5.64734D7, 5.6552D7, 5.67345D7, 5.8D7, 5.80598D7, 5.85093D7
     $     , 5.87672D7, 5.95014D7, 5.96716D7, 6.0D7, 6.05437D7,
     $     6.19719D7, 6.2D7, 6.4D7, 6.4388D7, 6.56996D7, 6.57039D7,
     $     6.6D7, 6.6479D7, 6.72731D7, 6.74974D7, 6.8D7, 6.82675D7,
     $     6.89516D7, 6.95072D7, 6.98592D7, 6.99106D7, 7.0D7, 7.19282D7,
     $     7.2D7, 7.4D7, 7.52314D7, 7.6D7, 7.8D7, 7.84305D7, 7.88071D7,
     $     7.92915D7, 8.0D7, 8.14563D7, 8.19372D7, 8.2D7, 8.35164D7,
     $     8.35846D7, 8.37881D7, 8.39869D7, 8.4D7, 8.6D7, 8.72029D7,
     $     8.8D7, 8.87169D7, 8.87511D7, 8.89006D7, 8.91458D7, 9.0D7,
     $     9.15399D7, 9.2D7, 9.4D7, 9.6D7, 9.7464D7, 9.8D7, 9.85211D7,
     $     9.95234D7, 9.98445D7, 1.0D8, 1.00184D8, 1.00638D8, 1.02D8,
     $     1.02554D8, 1.04D8, 1.05565D8, 1.06D8, 1.07148D8, 1.07777D8,
     $     1.08D8, 1.0979D8, 1.1D8, 1.11707D8, 1.12D8, 1.14D8, 1.14104D8
     $     , 1.16D8, 1.16252D8, 1.18D8, 1.18158D8, 1.2D8, 1.20086D8,
     $     1.20124D8, 1.20134D8, 1.20873D8, 1.22D8, 1.24D8, 1.26D8,
     $     1.28D8, 1.29049D8, 1.3D8, 1.30195D8, 1.3026D8, 1.30469D8,
     $     1.32D8, 1.33307D8, 1.33779D8, 1.34D8, 1.36D8, 1.38D8,
     $     1.39448D8, 1.4D8/

      data ( gdr(34,i), i=1, 336) /0.00165404, 0.00207144, 0.00263031,
     $     0.0029113, 0.00296773, 0.0030, 0.005322, 0.009178, 0.01551,
     $     0.02359, 0.031, 0.0339798, 0.03606, 0.04037, 0.04627, 0.05426
     $     , 0.06307, 0.07187, 0.08122, 0.0921, 0.105, 0.1199, 0.1365,
     $     0.1542, 0.154227, 0.156412, 0.1718, 0.1885, 0.191301, 0.2041,
     $     0.2178, 0.2289, 0.233656, 0.2359, 0.236704, 0.2385, 0.2368,
     $     0.2311, 0.226865, 0.2222, 0.2109, 0.21047, 0.1983, 0.185,
     $     0.1718, 0.167015, 0.16, 0.1506, 0.1436, 0.142951, 0.1375,
     $     0.13587, 0.1313, 0.126868, 0.1251, 0.1193, 0.115012, 0.1142,
     $     0.1098, 0.105953, 0.1059, 0.1025, 0.09955, 0.0982102, 0.09692
     $     , 0.0945, 0.0922, 0.0912801, 0.090432, 0.0901, 0.08841,
     $     0.0872, 0.08632, 0.0863043, 0.0856416, 0.08563, 0.08499,
     $     0.0849118, 0.08428, 0.08325, 0.08152, 0.0801163, 0.07885,
     $     0.0787043, 0.07578, 0.074628, 0.0743424, 0.073509, 0.0730927,
     $     0.07298, 0.07036, 0.06694, 0.06208, 0.0566, 0.0527387,
     $     0.05169, 0.0516314, 0.0494999, 0.04783, 0.045904, 0.0449,
     $     0.043, 0.04105, 0.03897, 0.038028, 0.0376328, 0.03673,
     $     0.03448, 0.0340607, 0.03243, 0.03057, 0.0299431, 0.02887,
     $     0.0280004, 0.02727, 0.0257, 0.02409, 0.0225, 0.02098, 0.01961
     $     , 0.01842, 0.01744, 0.01665, 0.016, 0.0157679, 0.01546,
     $     0.01502, 0.01466, 0.01436, 0.01406, 0.0138619, 0.01377,
     $     0.01349, 0.01324, 0.0131776, 0.01299, 0.01274, 0.01251,
     $     0.01229, 0.0122103, 0.01208, 0.0119458, 0.01189, 0.0117,
     $     0.01152, 0.01134, 0.0112626, 0.01118, 0.01103, 0.01089,
     $     0.01075, 0.01061, 0.01046, 0.0104013, 0.01033, 0.010261,
     $     0.01022, 0.0102063, 0.0101976, 0.01014, 0.01005, 0.00993803,
     $     0.009934, 0.009808, 0.00975895, 0.00971, 0.00969111, 0.009645
     $     , 0.009587, 0.009514, 0.009433, 0.00942773, 0.009359,
     $     0.009295, 0.009235, 0.009177, 0.009121, 0.00907, 0.00903589,
     $     0.00902, 0.008969, 0.00891, 0.008849, 0.008788, 0.00878035,
     $     0.008731, 0.008678, 0.00865845, 0.00863, 0.00861655, 0.008589
     $     , 0.008553, 0.008524, 0.0085, 0.008481, 0.008468, 0.008458,
     $     0.00845546, 0.008452, 0.008445, 0.00843499, 0.00837383,
     $     0.00836686, 0.00830818, 0.00827778, 0.008236, 0.00822272,
     $     0.00812222, 0.008056, 0.0080333, 0.00797558, 0.00794422,
     $     0.007932, 0.0079246, 0.00792015, 0.00791562, 0.00789544,
     $     0.00782726, 0.007813, 0.00779329, 0.00778161, 0.00775929,
     $     0.007748, 0.00772474, 0.007698, 0.00769865, 0.007703,
     $     0.00773945, 0.00774, 0.00774759, 0.00775029, 0.00775199,
     $     0.00775594, 0.007783, 0.00778457, 0.00779638, 0.00780312,
     $     0.00782217, 0.00782656, 0.007835, 0.00785173, 0.00789515,
     $     0.007896, 0.007944, 0.00795145, 0.00797636, 0.00797644,
     $     0.007982, 0.00799312, 0.00801141, 0.00801654, 0.008028,
     $     0.00803421, 0.00805002, 0.00806276, 0.0080708, 0.00807197,
     $     0.008074, 0.00807496, 0.008075, 0.008014, 0.00797868,
     $     0.007957, 0.00792, 0.00790604, 0.00789391, 0.00787843,
     $     0.007856, 0.00779235, 0.00777169, 0.007769, 0.00770047,
     $     0.00769743, 0.00768839, 0.00767958, 0.007679, 0.007593,
     $     0.00753424, 0.007496, 0.00745391, 0.00745192, 0.00744322,
     $     0.00742901, 0.00738, 0.0072903, 0.007264, 0.007143, 0.007048,
     $     0.00698182, 0.006958, 0.00692257, 0.00685545, 0.00683423,
     $     0.006824, 0.00681104, 0.00677936, 0.006686, 0.00665346,
     $     0.00657, 0.00647575, 0.00645, 0.00638466, 0.0063494, 0.006337
     $     , 0.0062545, 0.006245, 0.00613715, 0.006119, 0.006001,
     $     0.00599555, 0.005898, 0.00588482, 0.005795, 0.00578642,
     $     0.005688, 0.00568337, 0.0056813, 0.00568081, 0.00564126,
     $     0.005582, 0.005494, 0.005399, 0.00531, 0.00527198, 0.005238,
     $     0.00522847, 0.00522532, 0.00521522, 0.005142, 0.00508681,
     $     0.00506717, 0.005058, 0.004987, 0.004916, 0.00486367,
     $     0.004844/


C     41093, JENDL evaluated data
      nbins(35) = 359
      data ( eg(35,i), i=1, 359) /1931550.0, 6043150.0, 8555570.0,
     $     8831110.0, 8950000.0, 9000000.0, 9200000.0, 9234380.0,
     $     9400000.0, 9600000.0, 9800000.0, 9982720.0, 1.0D7, 1.02D7,
     $     1.02845D7, 1.04D7, 1.06D7, 1.08D7, 1.0968D7, 1.1D7, 1.12D7,
     $     1.14D7, 1.16D7, 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.245D7,
     $     1.24533D7, 1.26D7, 1.28D7, 1.3D7, 1.32D7, 1.33906D7,
     $     1.33942D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.41102D7,
     $     1.41552D7, 1.41607D7, 1.42D7, 1.42979D7, 1.44D7, 1.44449D7,
     $     1.46D7, 1.47948D7, 1.48D7, 1.5D7, 1.52D7, 1.52882D7, 1.54D7,
     $     1.54404D7, 1.56D7, 1.56521D7, 1.58D7, 1.6D7, 1.62D7, 1.64D7,
     $     1.64167D7, 1.66D7, 1.671D7, 1.67136D7, 1.68D7, 1.7D7,
     $     1.71934D7, 1.72D7, 1.72479D7, 1.74D7, 1.75949D7, 1.76D7,
     $     1.78D7, 1.8D7, 1.82D7, 1.82931D7, 1.84D7, 1.86D7, 1.88D7,
     $     1.88062D7, 1.9D7, 1.92D7, 1.92777D7, 1.94D7, 1.96D7,
     $     1.96827D7, 1.98D7, 1.9828D7, 2.0D7, 2.00494D7, 2.02D7,
     $     2.02309D7, 2.03633D7, 2.04D7, 2.06D7, 2.07426D7, 2.08D7,
     $     2.09882D7, 2.1D7, 2.12D7, 2.14D7, 2.15792D7, 2.16D7,
     $     2.17861D7, 2.18D7, 2.2D7, 2.22D7, 2.24D7, 2.26D7, 2.28D7,
     $     2.2884D7, 2.3D7, 2.31338D7, 2.32D7, 2.34D7, 2.36D7, 2.363D7,
     $     2.38D7, 2.38982D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7, 2.48D7,
     $     2.5D7, 2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.58776D7, 2.6D7,
     $     2.62D7, 2.63219D7, 2.64D7, 2.66D7, 2.68D7, 2.7D7, 2.7041D7,
     $     2.70523D7, 2.72D7, 2.74D7, 2.76D7, 2.76689D7, 2.78D7, 2.8D7,
     $     2.81571D7, 2.81577D7, 2.82D7, 2.82671D7, 2.84D7, 2.86D7,
     $     2.86009D7, 2.86491D7, 2.876D7, 2.87658D7, 2.88D7, 2.9D7,
     $     2.92D7, 2.94D7, 2.96D7, 2.98D7, 3.0D7, 3.02D7, 3.04D7, 3.06D7
     $     , 3.08D7, 3.08566D7, 3.1D7, 3.11937D7, 3.12D7, 3.14D7,
     $     3.15184D7, 3.16D7, 3.18D7, 3.18107D7, 3.2D7, 3.22D7, 3.24D7,
     $     3.26D7, 3.28D7, 3.3D7, 3.32D7, 3.34D7, 3.36D7, 3.38D7, 3.4D7,
     $     3.40374D7, 3.42D7, 3.44D7, 3.45758D7, 3.46D7, 3.46539D7,
     $     3.48D7, 3.5D7, 3.52D7, 3.54D7, 3.56D7, 3.58D7, 3.6D7, 3.62D7,
     $     3.64D7, 3.64303D7, 3.66D7, 3.68D7, 3.68287D7, 3.7D7, 3.72D7,
     $     3.74D7, 3.74848D7, 3.76D7, 3.78D7, 3.79566D7, 3.8D7, 3.82D7,
     $     3.82015D7, 3.83461D7, 3.84D7, 3.86D7, 3.88D7, 3.89156D7,
     $     3.9D7, 3.92D7, 3.94D7, 3.96D7, 3.98D7, 4.0D7, 4.0419D7,
     $     4.0759D7, 4.09591D7, 4.2D7, 4.2016D7, 4.25985D7, 4.2994D7,
     $     4.31026D7, 4.3462D7, 4.36922D7, 4.39614D7, 4.4D7, 4.45462D7,
     $     4.6D7, 4.67238D7, 4.8D7, 4.86944D7, 4.97247D7, 5.0D7,
     $     5.11414D7, 5.15677D7, 5.19959D7, 5.2D7, 5.22416D7, 5.38688D7,
     $     5.4D7, 5.4119D7, 5.47362D7, 5.54758D7, 5.6D7, 5.61976D7,
     $     5.62946D7, 5.70851D7, 5.8D7, 5.94429D7, 6.0D7, 6.07233D7,
     $     6.12362D7, 6.14538D7, 6.19463D7, 6.2D7, 6.30037D7, 6.36907D7,
     $     6.4D7, 6.51028D7, 6.5729D7, 6.57782D7, 6.6D7, 6.64085D7,
     $     6.76906D7, 6.8D7, 7.0D7, 7.0978D7, 7.2D7, 7.3881D7, 7.4D7,
     $     7.4016D7, 7.52892D7, 7.6D7, 7.67498D7, 7.69953D7, 7.73D7,
     $     7.79513D7, 7.8D7, 7.93691D7, 8.0D7, 8.2D7, 8.20015D7,
     $     8.22467D7, 8.31997D7, 8.4D7, 8.44032D7, 8.44766D7, 8.46247D7,
     $     8.52362D7, 8.6D7, 8.8D7, 8.98699D7, 9.0D7, 9.0429D7, 9.2D7,
     $     9.29674D7, 9.34723D7, 9.4D7, 9.50299D7, 9.6D7, 9.62373D7,
     $     9.79716D7, 9.8D7, 9.871D7, 1.0D8, 1.01119D8, 1.02D8, 1.04D8,
     $     1.04196D8, 1.05121D8, 1.05811D8, 1.06D8, 1.07723D8, 1.08D8,
     $     1.08962D8, 1.1D8, 1.10861D8, 1.11303D8, 1.12D8, 1.12735D8,
     $     1.14D8, 1.14627D8, 1.15384D8, 1.16D8, 1.18D8, 1.2D8,
     $     1.20802D8, 1.2124D8, 1.22D8, 1.23019D8, 1.23284D8, 1.24D8,
     $     1.26D8, 1.26384D8, 1.2778D8, 1.28D8, 1.3D8, 1.30681D8, 1.32D8
     $     , 1.34D8, 1.34059D8, 1.34235D8, 1.36D8, 1.36605D8, 1.38D8,
     $     1.38595D8, 1.4D8/

      data ( gdr(35,i), i=1, 359) /4.41127D-4, 0.00138013, 0.00195392,
     $     0.00201685, 0.002044, 0.002852, 0.005592, 0.00591954,
     $     0.007764, 0.009471, 0.01081, 0.0118025, 0.0119, 0.01283,
     $     0.0131968, 0.01371, 0.01463, 0.01571, 0.0168296, 0.01705,
     $     0.01874, 0.0209, 0.02362, 0.02686, 0.0301, 0.0334, 0.03687,
     $     0.03777, 0.0378309, 0.0406, 0.04468, 0.04923, 0.05433,
     $     0.0598366, 0.0599457, 0.06012, 0.06671, 0.07423, 0.0828,
     $     0.0880533, 0.0902843, 0.0905567, 0.09255, 0.0978234, 0.1036,
     $     0.106319, 0.1162, 0.132245, 0.1327, 0.1504, 0.1667, 0.173051,
     $     0.1814, 0.183909, 0.1941, 0.196744, 0.2044, 0.2122, 0.2173,
     $     0.22, 0.220034, 0.2204, 0.2197, 0.219663, 0.2188, 0.2155,
     $     0.210759, 0.2106, 0.209115, 0.2045, 0.197481, 0.1973, 0.1891,
     $     0.1809, 0.1731, 0.169747, 0.166, 0.1594, 0.1532, 0.15302,
     $     0.1475, 0.142, 0.140071, 0.1371, 0.1323, 0.130457, 0.1279,
     $     0.127286, 0.1236, 0.122546, 0.1194, 0.11877, 0.11612, 0.1154,
     $     0.1112, 0.108184, 0.107, 0.103043, 0.1028, 0.09876, 0.09475,
     $     0.0912732, 0.09088, 0.0873468, 0.08709, 0.08342, 0.07986,
     $     0.07642, 0.0731, 0.0699, 0.0685949, 0.06684, 0.064849,
     $     0.06389, 0.06107, 0.05838, 0.0579842, 0.0558, 0.0545707,
     $     0.05333, 0.05099, 0.04874, 0.04662, 0.04458, 0.04266, 0.04082
     $     , 0.03906, 0.03737, 0.03575, 0.0351429, 0.03421, 0.03275,
     $     0.031907, 0.03138, 0.03007, 0.02885, 0.02769, 0.0274643,
     $     0.0274028, 0.02661, 0.0256, 0.02466, 0.0243488, 0.02377,
     $     0.02295, 0.0223348, 0.0223321, 0.02217, 0.0219286, 0.02146,
     $     0.02078, 0.0207775, 0.0206467, 0.02035, 0.0203254, 0.02018,
     $     0.01961, 0.01908, 0.01859, 0.01815, 0.01774, 0.01735, 0.017,
     $     0.01668, 0.01638, 0.01611, 0.0160358, 0.01585, 0.0156175,
     $     0.01561, 0.01539, 0.0152593, 0.01517, 0.01497, 0.0149592,
     $     0.01477, 0.01458, 0.0144, 0.01422, 0.01405, 0.01388, 0.01373,
     $     0.01357, 0.01343, 0.01329, 0.01316, 0.0131336, 0.01302,
     $     0.0129, 0.0127944, 0.01278, 0.0127475, 0.01266, 0.01255,
     $     0.01244, 0.01233, 0.01224, 0.01214, 0.01205, 0.01196, 0.01187
     $     , 0.0118578, 0.01179, 0.01171, 0.0116985, 0.01163, 0.01157,
     $     0.01149, 0.0114602, 0.01142, 0.01136, 0.0113051, 0.01129,
     $     0.01124, 0.0112395, 0.0111961, 0.01118, 0.01112, 0.01106,
     $     0.011031, 0.01101, 0.01095, 0.0109, 0.01086, 0.01081, 0.01076
     $     , 0.0105975, 0.0104687, 0.0103941, 0.01002, 0.0100176,
     $     0.00993275, 0.00987619, 0.0098608, 0.00981032, 0.00977833,
     $     0.00974128, 0.009736, 0.00967415, 0.009515, 0.00945487,
     $     0.009352, 0.00931031, 0.00924987, 0.009234, 0.00917979,
     $     0.00915993, 0.00914019, 0.00914, 0.0091315, 0.00907543,
     $     0.009071, 0.00906863, 0.00905646, 0.00904207, 0.009032,
     $     0.00902929, 0.00902796, 0.00901722, 0.009005, 0.00899195,
     $     0.008987, 0.00898152, 0.00897767, 0.00897604, 0.0089724,
     $     0.008972, 0.00896188, 0.00895505, 0.008952, 0.00893867,
     $     0.00893121, 0.00893062, 0.008928, 0.00892179, 0.00890258,
     $     0.008898, 0.008859, 0.0088373, 0.008815, 0.00877075, 0.008768
     $     , 0.00876752, 0.00872974, 0.008709, 0.00868476, 0.00867688,
     $     0.00866716, 0.00864653, 0.008645, 0.00859476, 0.008572,
     $     0.008492, 0.00849193, 0.00848073, 0.00843764, 0.008402,
     $     0.00838279, 0.0083793, 0.00837229, 0.00834351, 0.008308,
     $     0.008206, 0.00810584, 0.008099, 0.00807463, 0.007987,
     $     0.00793282, 0.00790491, 0.007876, 0.00781418, 0.007757,
     $     0.00774253, 0.00763867, 0.007637, 0.0075921, 0.007512,
     $     0.00744265, 0.007389, 0.007265, 0.00725266, 0.00719491,
     $     0.00715252, 0.007141, 0.00703304, 0.007016, 0.0069558,
     $     0.006892, 0.00683935, 0.00681265, 0.006771, 0.00672567,
     $     0.006649, 0.00661124, 0.0065662, 0.00653, 0.006414, 0.0063,
     $     0.00625624, 0.00623259, 0.006192, 0.0061365, 0.00612225,
     $     0.006084, 0.005977, 0.00595734, 0.00588693, 0.005876,
     $     0.005778, 0.00574533, 0.005683, 0.005594, 0.00559138,
     $     0.00558363, 0.005507, 0.00548195, 0.005425, 0.00540186,
     $     0.005348/


C     42092, JENDL evaluated data
      nbins(36) = 190
      data ( eg(36,i), i=1, 190) /1.28D7, 1.3D7, 1.32D7, 1.34D7,
     $     1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7
     $     ,1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7, 1.62D7, 1.64D7,
     $     1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7
     $     ,1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7, 1.92D7, 1.94D7,
     $     1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7, 2.08D7, 2.1D7
     $     ,2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7, 2.22D7, 2.24D7,
     $     2.26D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7, 2.38D7, 2.4D7
     $     ,2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7, 2.52D7, 2.54D7,
     $     2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7, 2.68D7, 2.7D7
     $     ,2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7, 2.82D7, 2.84D7,
     $     2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7, 2.98D7, 3.0D7
     $     ,3.02D7, 3.04D7, 3.06D7, 3.08D7, 3.1D7, 3.12D7, 3.14D7,
     $     3.16D7, 3.18D7, 3.2D7, 3.22D7, 3.24D7, 3.26D7, 3.28D7, 3.3D7
     $     ,3.32D7, 3.34D7, 3.36D7, 3.38D7, 3.4D7, 3.42D7, 3.44D7,
     $     3.46D7, 3.48D7, 3.5D7, 3.52D7, 3.54D7, 3.56D7, 3.58D7, 3.6D7
     $     ,3.62D7, 3.64D7, 3.66D7, 3.68D7, 3.7D7, 3.72D7, 3.74D7,
     $     3.76D7, 3.78D7, 3.8D7, 3.82D7, 3.84D7, 3.86D7, 3.88D7, 3.9D7
     $     ,3.92D7, 3.94D7, 3.96D7, 3.98D7, 4.0D7, 4.1D7, 4.2D7, 4.3D7
     $     ,4.4D7, 4.5D7, 4.6D7, 4.8D7, 5.0D7, 5.2D7, 5.4D7, 5.6D7,
     $     5.8D7, 6.0D7, 6.2D7, 6.4D7, 6.6D7, 6.8D7, 7.0D7, 7.2D7, 7.4D7
     $     ,7.6D7, 7.8D7, 8.0D7, 8.2D7, 8.4D7, 8.6D7, 8.8D7, 9.0D7,
     $     9.2D7, 9.4D7, 9.6D7, 9.8D7, 1.0D8, 1.02D8, 1.04D8, 1.06D8,
     $     1.08D8,1.1D8, 1.12D8, 1.14D8, 1.16D8, 1.18D8, 1.2D8, 1.22D8,
     $     1.24D8,1.26D8, 1.28D8, 1.3D8, 1.32D8, 1.34D8, 1.36D8, 1.38D8,
     $     1.4D8/

      data ( gdr(36,i), i=1, 190) /0.0043, 0.009969, 0.01435, 0.01859,
     $     0.02349, 0.02945, 0.03723, 0.04756, 0.06047, 0.07545, 0.09211
     $     , 0.1101, 0.1287, 0.1474, 0.1655, 0.1822, 0.1969, 0.209,
     $     0.2177, 0.2225, 0.2227, 0.2189, 0.2123, 0.2041, 0.1954,
     $     0.1875, 0.18, 0.1733, 0.1669, 0.1612, 0.1556, 0.1505, 0.1455,
     $     0.1407, 0.1359, 0.1311, 0.1261, 0.1209, 0.1154, 0.1095,
     $     0.1031, 0.09637, 0.08976, 0.08381, 0.07907, 0.07607, 0.07539,
     $     0.07758, 0.08075, 0.08271, 0.0836, 0.08359, 0.08276, 0.08124,
     $     0.07918, 0.0767, 0.07394, 0.07104, 0.06814, 0.06536, 0.06291,
     $     0.06095, 0.05965, 0.059, 0.05917, 0.05942, 0.05958, 0.05968,
     $     0.05972, 0.05966, 0.05953, 0.05932, 0.05905, 0.0587, 0.0583,
     $     0.05782, 0.0573, 0.05677, 0.0562, 0.05557, 0.05488, 0.05416,
     $     0.05341, 0.05261, 0.05179, 0.05091, 0.04967, 0.04873, 0.04772
     $     , 0.04669, 0.04562, 0.04452, 0.04338, 0.04225, 0.0411,
     $     0.03993, 0.03875, 0.03759, 0.03641, 0.03525, 0.03408, 0.03294
     $     , 0.03189, 0.03091, 0.03002, 0.02919, 0.02843, 0.02772,
     $     0.02708, 0.02649, 0.02595, 0.02532, 0.02502, 0.02471, 0.02443
     $     , 0.02416, 0.02389, 0.02363, 0.02337, 0.02313, 0.0229,
     $     0.02267, 0.02245, 0.02223, 0.02202, 0.02182, 0.02161, 0.02143
     $     , 0.02124, 0.02106, 0.02087, 0.02071, 0.02055, 0.02039,
     $     0.02023, 0.02008, 0.01993, 0.01917, 0.01857, 0.01803, 0.01754
     $     , 0.01709, 0.01667, 0.01593, 0.01532, 0.01483, 0.01441,
     $     0.01407, 0.01378, 0.01352, 0.01327, 0.01306, 0.01285, 0.01265
     $     , 0.01245, 0.01223, 0.01201, 0.01178, 0.01154, 0.01129,
     $     0.01104, 0.01078, 0.01052, 0.01025, 0.009988, 0.009725,
     $     0.009466, 0.009211, 0.008959, 0.008713, 0.008472, 0.008241,
     $     0.008012, 0.007792, 0.00758, 0.007372, 0.007174, 0.006981,
     $     0.006796, 0.006617, 0.006446, 0.006282, 0.006124, 0.005972,
     $     0.005827, 0.005687, 0.005553, 0.005426, 0.005307, 0.005192/


C     42094, JENDL evaluated data
      nbins(37) = 210
      data ( eg(37,i), i=1, 210) /9700000.0, 9800000.0, 9900000.0,
     $     1.0D7, 1.02D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7,
     $     1.16D7, 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7,
     $     1.32D7, 1.34D7, 1.35D7, 1.36D7, 1.37D7, 1.38D7, 1.39D7, 1.4D7
     $     , 1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7, 1.51D7, 1.52D7,
     $     1.54D7, 1.56D7, 1.58D7, 1.6D7, 1.62D7, 1.64D7, 1.66D7, 1.68D7
     $     , 1.7D7, 1.72D7, 1.74D7, 1.76D7, 1.775D7, 1.78D7, 1.8D7,
     $     1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7
     $     , 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7, 2.08D7, 2.1D7,
     $     2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7, 2.22D7, 2.24D7, 2.26D7
     $     , 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7, 2.38D7, 2.4D7,
     $     2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7, 2.52D7, 2.54D7, 2.56D7
     $     , 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7, 2.68D7, 2.7D7,
     $     2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7, 2.82D7, 2.84D7, 2.86D7
     $     , 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7, 2.98D7, 3.0D7,
     $     3.02D7, 3.04D7, 3.043D7, 3.06D7, 3.08D7, 3.1D7, 3.12D7,
     $     3.14D7, 3.16D7, 3.18D7, 3.2D7, 3.22D7, 3.24D7, 3.26D7, 3.28D7
     $     , 3.3D7, 3.32D7, 3.34D7, 3.36D7, 3.38D7, 3.4D7, 3.42D7,
     $     3.44D7, 3.46D7, 3.48D7, 3.5D7, 3.52D7, 3.54D7, 3.56D7, 3.58D7
     $     , 3.6D7, 3.62D7, 3.64D7, 3.66D7, 3.68D7, 3.7D7, 3.72D7,
     $     3.74D7, 3.76D7, 3.78D7, 3.8D7, 3.82D7, 3.84D7, 3.86D7, 3.88D7
     $     , 3.9D7, 3.92D7, 3.94D7, 3.96D7, 3.98D7, 4.0D7, 4.2D7, 4.4D7,
     $     4.6D7, 4.8D7, 5.0D7, 5.2D7, 5.4D7, 5.6D7, 5.8D7, 6.0D7, 6.2D7
     $     , 6.4D7, 6.6D7, 6.8D7, 7.0D7, 7.2D7, 7.4D7, 7.6D7, 7.8D7,
     $     8.0D7, 8.2D7, 8.4D7, 8.6D7, 8.8D7, 9.0D7, 9.2D7, 9.4D7, 9.6D7
     $     , 9.8D7, 1.0D8, 1.02D8, 1.04D8, 1.06D8, 1.08D8, 1.1D8, 1.12D8
     $     , 1.14D8, 1.16D8, 1.18D8, 1.2D8, 1.22D8, 1.24D8, 1.26D8,
     $     1.28D8, 1.3D8, 1.32D8, 1.34D8, 1.36D8, 1.38D8, 1.4D8/

      data ( gdr(37,i), i=1, 210) /0.0, 0.007499, 0.009201, 0.01021,
     $     0.01114, 0.01221, 0.01538, 0.01832, 0.01994, 0.02188, 0.02417
     $     , 0.02683, 0.02988, 0.03335, 0.03726, 0.04163, 0.04647,
     $     0.05183, 0.05771, 0.06415, 0.07116, 0.07457, 0.06606, 0.06484
     $     , 0.07157, 0.079, 0.08632, 0.1005, 0.1141, 0.1271, 0.1391,
     $     0.1501, 0.1552, 0.1599, 0.1684, 0.1755, 0.1811, 0.185, 0.1869
     $     , 0.1869, 0.1848, 0.1804, 0.1743, 0.1682, 0.1621, 0.1559,
     $     0.1514, 0.1497, 0.1435, 0.1374, 0.1314, 0.1254, 0.1195,
     $     0.1138, 0.1084, 0.1031, 0.09801, 0.09319, 0.08868, 0.08447,
     $     0.08063, 0.07718, 0.07414, 0.07154, 0.06942, 0.06759, 0.0658,
     $     0.06401, 0.06225, 0.06053, 0.05883, 0.05714, 0.05548, 0.05385
     $     , 0.05221, 0.05059, 0.04898, 0.04739, 0.04582, 0.04428,
     $     0.04276, 0.04127, 0.03982, 0.0384, 0.03701, 0.03568, 0.03437,
     $     0.0331, 0.03187, 0.03066, 0.02948, 0.02833, 0.02722, 0.02615,
     $     0.02513, 0.02416, 0.02323, 0.02236, 0.02153, 0.02075, 0.02004
     $     , 0.01936, 0.01874, 0.01818, 0.01765, 0.01716, 0.01672,
     $     0.01633, 0.01596, 0.01562, 0.01532, 0.01528, 0.01505, 0.01478
     $     , 0.01455, 0.01434, 0.01414, 0.01396, 0.01378, 0.01362,
     $     0.01345, 0.01329, 0.01315, 0.013, 0.01287, 0.01273, 0.0126,
     $     0.01248, 0.01235, 0.01223, 0.01212, 0.01201, 0.01191, 0.01181
     $     , 0.01172, 0.01162, 0.01153, 0.01144, 0.01136, 0.01128,
     $     0.0112, 0.01112, 0.01106, 0.01098, 0.01092, 0.01086, 0.0108,
     $     0.01074, 0.01068, 0.01062, 0.01056, 0.01051, 0.01045, 0.0104,
     $     0.01035, 0.0103, 0.01025, 0.0102, 0.01016, 0.01011, 0.009679,
     $     0.009437, 0.009311, 0.009307, 0.009375, 0.009465, 0.00953,
     $     0.009578, 0.009517, 0.009648, 0.009686, 0.009715, 0.009745,
     $     0.009762, 0.009776, 0.009773, 0.009764, 0.009755, 0.00973,
     $     0.009696, 0.009653, 0.009596, 0.009532, 0.009437, 0.009319,
     $     0.009187, 0.009042, 0.008885, 0.008717, 0.008547, 0.00838,
     $     0.008206, 0.008034, 0.007861, 0.00769, 0.007526, 0.007365,
     $     0.007207, 0.007055, 0.006908, 0.006768, 0.00663, 0.006498,
     $     0.00637, 0.006247, 0.006129, 0.006017, 0.005912, 0.005812,
     $     0.005718/


C     42096, JENDL evaluated data
      nbins(38) = 210
      data ( eg(38,i), i=1, 210) /9300000.0, 9400000.0, 9500000.0,
     $     9600000.0, 9700000.0, 9800000.0, 9900000.0, 1.0D7, 1.02D7,
     $     1.04D7, 1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7, 1.18D7
     $     , 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7, 1.32D7,
     $     1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7, 1.48D7
     $     , 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7, 1.62D7,
     $     1.64D7, 1.652D7, 1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7,
     $     1.76D7, 1.778D7, 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7,
     $     1.88D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.02D7,
     $     2.04D7, 2.06D7, 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7, 2.18D7
     $     , 2.2D7, 2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7, 2.32D7,
     $     2.34D7, 2.36D7, 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7, 2.48D7
     $     , 2.5D7, 2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7,
     $     2.64D7, 2.66D7, 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7, 2.78D7
     $     , 2.8D7, 2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7,
     $     2.94D7, 2.96D7, 2.98D7, 3.0D7, 3.02D7, 3.04D7, 3.06D7, 3.08D7
     $     , 3.1D7, 3.12D7, 3.14D7, 3.16D7, 3.18D7, 3.2D7, 3.22D7,
     $     3.24D7, 3.26D7, 3.28D7, 3.3D7, 3.32D7, 3.34D7, 3.36D7, 3.38D7
     $     , 3.4D7, 3.42D7, 3.44D7, 3.46D7, 3.48D7, 3.5D7, 3.52D7,
     $     3.54D7, 3.56D7, 3.58D7, 3.6D7, 3.62D7, 3.64D7, 3.66D7, 3.68D7
     $     , 3.7D7, 3.72D7, 3.74D7, 3.76D7, 3.78D7, 3.8D7, 3.82D7,
     $     3.84D7, 3.86D7, 3.88D7, 3.9D7, 3.92D7, 3.94D7, 3.96D7, 3.98D7
     $     , 4.0D7, 4.2D7, 4.4D7, 4.6D7, 4.8D7, 5.0D7, 5.2D7, 5.4D7,
     $     5.6D7, 5.8D7, 6.0D7, 6.2D7, 6.4D7, 6.6D7, 6.8D7, 7.0D7, 7.2D7
     $     , 7.4D7, 7.6D7, 7.8D7, 8.0D7, 8.2D7, 8.4D7, 8.6D7, 8.8D7,
     $     9.0D7, 9.2D7, 9.4D7, 9.6D7, 9.8D7, 1.0D8, 1.02D8, 1.04D8,
     $     1.06D8, 1.08D8, 1.1D8, 1.12D8, 1.14D8, 1.16D8, 1.18D8, 1.2D8,
     $     1.22D8, 1.24D8, 1.26D8, 1.28D8, 1.3D8, 1.32D8, 1.34D8, 1.36D8
     $     , 1.38D8, 1.4D8/

      data ( gdr(38,i), i=1, 210) /0.01636, 0.02115, 0.02225, 0.02147,
     $     0.02052, 0.01996, 0.01975, 0.01986, 0.02087, 0.02267, 0.02492
     $     , 0.02734, 0.0299, 0.03262, 0.03555, 0.03873, 0.0422, 0.046,
     $     0.05017, 0.05475, 0.05979, 0.06532, 0.07138, 0.07802, 0.08527
     $     , 0.09318, 0.1017, 0.1111, 0.1212, 0.1321, 0.1433, 0.1531,
     $     0.1617, 0.1689, 0.1748, 0.1794, 0.1827, 0.1845, 0.1851,
     $     0.1841, 0.1832, 0.1817, 0.178, 0.1727, 0.1659, 0.1588, 0.1522
     $     , 0.1466, 0.146, 0.1402, 0.1349, 0.1299, 0.1251, 0.1206,
     $     0.1161, 0.1118, 0.1076, 0.1033, 0.09913, 0.09474, 0.09019,
     $     0.08546, 0.08099, 0.07692, 0.0732, 0.06982, 0.06678, 0.06406,
     $     0.06164, 0.05948, 0.0576, 0.05594, 0.0545, 0.05326, 0.05218,
     $     0.05127, 0.05049, 0.04979, 0.04916, 0.0486, 0.04807, 0.04757,
     $     0.04706, 0.04654, 0.04598, 0.04538, 0.04469, 0.0439, 0.04301,
     $     0.04197, 0.04078, 0.03941, 0.03783, 0.03603, 0.03399, 0.03168
     $     , 0.02909, 0.0262, 0.02396, 0.02299, 0.02255, 0.02216,
     $     0.02179, 0.02141, 0.02107, 0.02073, 0.02039, 0.02006, 0.01976
     $     , 0.01945, 0.01917, 0.01889, 0.0186, 0.01834, 0.01808,
     $     0.01783, 0.01759, 0.01735, 0.01714, 0.01692, 0.0167, 0.0165,
     $     0.0163, 0.01612, 0.01593, 0.01575, 0.01558, 0.01541, 0.01525,
     $     0.0151, 0.01495, 0.01481, 0.01467, 0.01453, 0.0144, 0.01427,
     $     0.01415, 0.01403, 0.01391, 0.0138, 0.01368, 0.01359, 0.01348,
     $     0.01338, 0.01328, 0.01319, 0.0131, 0.01301, 0.01293, 0.01284,
     $     0.01276, 0.01268, 0.0126, 0.01253, 0.01246, 0.01239, 0.01232,
     $     0.01226, 0.01219, 0.01213, 0.01158, 0.01114, 0.01076, 0.01051
     $     , 0.01029, 0.01012, 0.009971, 0.009844, 0.009739, 0.009636,
     $     0.009542, 0.009455, 0.009375, 0.009299, 0.009226, 0.009155,
     $     0.009083, 0.009012, 0.008936, 0.008863, 0.008776, 0.008683,
     $     0.008585, 0.008477, 0.008364, 0.008246, 0.008124, 0.007997,
     $     0.007877, 0.007747, 0.007616, 0.007487, 0.007356, 0.007228,
     $     0.0071, 0.006972, 0.006847, 0.006724, 0.006603, 0.00649,
     $     0.006374, 0.006261, 0.006151, 0.006045, 0.005943, 0.005845,
     $     0.00575, 0.005662, 0.005576, 0.005496/


C     42098, JENDL evaluated data
      nbins(39) = 210
      data ( eg(39,i), i=1, 210) /8800000.0, 9000000.0, 9200000.0,
     $     9400000.0, 9600000.0, 9800000.0, 1.0D7, 1.02D7, 1.04D7,
     $     1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7, 1.18D7, 1.2D7,
     $     1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7, 1.32D7, 1.34D7, 1.36D7
     $     , 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7,
     $     1.52D7, 1.54D7, 1.546D7, 1.56D7, 1.564D7, 1.58D7, 1.6D7,
     $     1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7
     $     , 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7,
     $     1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7
     $     , 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7,
     $     2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7
     $     , 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7, 2.461D7, 2.48D7,
     $     2.5D7, 2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7,
     $     2.66D7, 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7,
     $     2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7
     $     , 2.98D7, 3.0D7, 3.02D7, 3.04D7, 3.06D7, 3.08D7, 3.1D7,
     $     3.12D7, 3.14D7, 3.16D7, 3.18D7, 3.2D7, 3.22D7, 3.24D7, 3.26D7
     $     , 3.28D7, 3.3D7, 3.32D7, 3.34D7, 3.36D7, 3.38D7, 3.4D7,
     $     3.42D7, 3.44D7, 3.46D7, 3.48D7, 3.5D7, 3.52D7, 3.54D7, 3.56D7
     $     , 3.58D7, 3.6D7, 3.62D7, 3.64D7, 3.66D7, 3.68D7, 3.7D7,
     $     3.72D7, 3.74D7, 3.76D7, 3.78D7, 3.8D7, 3.82D7, 3.84D7, 3.86D7
     $     , 3.88D7, 3.9D7, 3.92D7, 3.94D7, 3.96D7, 3.98D7, 4.0D7, 4.2D7
     $     , 4.4D7, 4.6D7, 4.8D7, 5.0D7, 5.2D7, 5.4D7, 5.6D7, 5.8D7,
     $     6.0D7, 6.2D7, 6.4D7, 6.6D7, 6.8D7, 7.0D7, 7.2D7, 7.4D7, 7.6D7
     $     , 7.8D7, 8.0D7, 8.2D7, 8.4D7, 8.6D7, 8.8D7, 9.0D7, 9.2D7,
     $     9.4D7, 9.6D7, 9.8D7, 1.0D8, 1.02D8, 1.04D8, 1.06D8, 1.08D8,
     $     1.1D8, 1.12D8, 1.14D8, 1.16D8, 1.18D8, 1.2D8, 1.22D8, 1.24D8,
     $     1.26D8, 1.28D8, 1.3D8, 1.32D8, 1.34D8, 1.36D8, 1.38D8, 1.4D8/

      data ( gdr(39,i), i=1, 210) /0.003711, 0.008138, 0.01204, 0.01551,
     $     0.01863, 0.02148, 0.02414, 0.02669, 0.02923, 0.03182, 0.03457
     $     , 0.03754, 0.04082, 0.04449, 0.04865, 0.05334, 0.05857,
     $     0.06429, 0.07049, 0.07713, 0.08417, 0.09159, 0.09936, 0.1074,
     $     0.1158, 0.1244, 0.1332, 0.1423, 0.1515, 0.1603, 0.1684,
     $     0.1756, 0.1816, 0.1862, 0.1871, 0.1892, 0.1895, 0.1906,
     $     0.1899, 0.1871, 0.1825, 0.1766, 0.1699, 0.1631, 0.1563,
     $     0.1496, 0.1429, 0.1364, 0.1299, 0.1236, 0.1175, 0.1115,
     $     0.1057, 0.1002, 0.09482, 0.08974, 0.08493, 0.08044, 0.07627,
     $     0.07245, 0.06899, 0.06592, 0.06328, 0.06108, 0.05935, 0.0581,
     $     0.05737, 0.05717, 0.05754, 0.0585, 0.06006, 0.06219, 0.06393,
     $     0.06501, 0.06548, 0.06541, 0.06482, 0.06378, 0.06235, 0.06059
     $     , 0.05857, 0.05634, 0.05622, 0.05396, 0.05149, 0.04901,
     $     0.04656, 0.04422, 0.04204, 0.04001, 0.03815, 0.03643, 0.03485
     $     , 0.0334, 0.03207, 0.03085, 0.02973, 0.0287, 0.02775, 0.02688
     $     , 0.02607, 0.02532, 0.02461, 0.02395, 0.02334, 0.02276,
     $     0.02223, 0.02171, 0.02124, 0.02078, 0.02036, 0.01997, 0.01958
     $     , 0.01922, 0.01886, 0.01854, 0.01822, 0.01793, 0.01765,
     $     0.01738, 0.01712, 0.01687, 0.01664, 0.0164, 0.01619, 0.01599,
     $     0.01578, 0.01559, 0.01541, 0.01522, 0.01505, 0.01489, 0.01473
     $     , 0.01458, 0.01443, 0.01429, 0.01415, 0.01402, 0.01389,
     $     0.01376, 0.01364, 0.01354, 0.01342, 0.01331, 0.0132, 0.01311,
     $     0.01301, 0.01291, 0.01283, 0.01274, 0.01265, 0.01257, 0.01248
     $     , 0.01241, 0.01233, 0.01226, 0.01219, 0.01212, 0.01205,
     $     0.01199, 0.01144, 0.01093, 0.01061, 0.01035, 0.01014,
     $     0.009966, 0.009815, 0.009683, 0.009566, 0.009461, 0.009368,
     $     0.009284, 0.009208, 0.009135, 0.009064, 0.008995, 0.008922,
     $     0.008856, 0.008779, 0.008698, 0.008608, 0.008509, 0.008405,
     $     0.008296, 0.008181, 0.008065, 0.007945, 0.00783, 0.007707,
     $     0.007581, 0.007457, 0.007333, 0.007208, 0.007083, 0.006958,
     $     0.006837, 0.006717, 0.006599, 0.006484, 0.006372, 0.006261,
     $     0.006151, 0.006047, 0.005945, 0.005846, 0.005752, 0.00566,
     $     0.005573, 0.005491, 0.005414/


C     42100, JENDL evaluated data
      nbins(40) = 212
      data ( eg(40,i), i=1, 212) /8400000.0, 8600000.0, 8800000.0,
     $     9000000.0, 9200000.0, 9400000.0, 9600000.0, 9800000.0, 1.0D7,
     $     1.02D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7
     $     , 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7,
     $     1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7
     $     , 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7,
     $     1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7
     $     , 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7,
     $     1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7
     $     , 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7,
     $     2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7
     $     , 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7,
     $     2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7
     $     , 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7,
     $     2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7
     $     , 2.98D7, 3.0D7, 3.02D7, 3.04D7, 3.06D7, 3.08D7, 3.1D7,
     $     3.12D7, 3.14D7, 3.16D7, 3.18D7, 3.2D7, 3.22D7, 3.24D7, 3.26D7
     $     , 3.28D7, 3.3D7, 3.32D7, 3.34D7, 3.36D7, 3.38D7, 3.4D7,
     $     3.42D7, 3.44D7, 3.46D7, 3.48D7, 3.5D7, 3.52D7, 3.54D7, 3.56D7
     $     , 3.58D7, 3.6D7, 3.62D7, 3.64D7, 3.66D7, 3.68D7, 3.7D7,
     $     3.72D7, 3.74D7, 3.76D7, 3.78D7, 3.8D7, 3.82D7, 3.84D7, 3.86D7
     $     , 3.88D7, 3.9D7, 3.92D7, 3.94D7, 3.96D7, 3.98D7, 4.0D7, 4.1D7
     $     , 4.2D7, 4.3D7, 4.4D7, 4.5D7, 4.6D7, 4.8D7, 5.0D7, 5.2D7,
     $     5.4D7, 5.6D7, 5.8D7, 6.0D7, 6.2D7, 6.4D7, 6.6D7, 6.8D7, 7.0D7
     $     , 7.2D7, 7.4D7, 7.6D7, 7.8D7, 8.0D7, 8.2D7, 8.4D7, 8.6D7,
     $     8.8D7, 9.0D7, 9.2D7, 9.4D7, 9.6D7, 9.8D7, 1.0D8, 1.02D8,
     $     1.04D8, 1.06D8, 1.08D8, 1.1D8, 1.12D8, 1.14D8, 1.16D8, 1.18D8
     $     , 1.2D8, 1.22D8, 1.24D8, 1.26D8, 1.28D8, 1.3D8, 1.32D8,
     $     1.34D8, 1.36D8, 1.38D8, 1.4D8/

      data ( gdr(40,i), i=1, 212) /0.001212, 0.002953, 0.005271, 0.0080,
     $     0.011, 0.01427, 0.01781, 0.02169, 0.026, 0.03072, 0.03538,
     $     0.03942, 0.043, 0.047, 0.05201, 0.05742, 0.06238, 0.06703,
     $     0.0725, 0.07949, 0.08678, 0.09273, 0.09749, 0.103, 0.1106,
     $     0.1193, 0.1271, 0.1337, 0.1397, 0.1457, 0.1514, 0.1566,
     $     0.1611, 0.1648, 0.1677, 0.1697, 0.1708, 0.171, 0.1703, 0.1688
     $     , 0.1665, 0.1636, 0.16, 0.1561, 0.1518, 0.1472, 0.1425,
     $     0.1376, 0.1327, 0.1278, 0.123, 0.1182, 0.1137, 0.1093, 0.105,
     $     0.1008, 0.09666, 0.0927, 0.08932, 0.08675, 0.08453, 0.082,
     $     0.07903, 0.07602, 0.07328, 0.07078, 0.06842, 0.06614, 0.06396
     $     , 0.06191, 0.05996, 0.0581, 0.05633, 0.05465, 0.05304, 0.0515
     $     , 0.05003, 0.04861, 0.04725, 0.04596, 0.04471, 0.04352,
     $     0.04237, 0.04128, 0.04023, 0.03924, 0.03829, 0.03738, 0.03651
     $     , 0.03567, 0.03487, 0.0341, 0.03335, 0.03263, 0.03196,
     $     0.03129, 0.03067, 0.03007, 0.02948, 0.02891, 0.02838, 0.02785
     $     , 0.02735, 0.02686, 0.02641, 0.02596, 0.02554, 0.02513,
     $     0.02474, 0.02435, 0.02398, 0.02362, 0.02327, 0.02294, 0.02261
     $     , 0.0223, 0.02199, 0.0217, 0.02141, 0.02114, 0.02088, 0.02062
     $     , 0.02037, 0.02013, 0.01989, 0.01967, 0.01944, 0.01923,
     $     0.01902, 0.01881, 0.01861, 0.01842, 0.01824, 0.01806, 0.01788
     $     , 0.01771, 0.01754, 0.01739, 0.01724, 0.01709, 0.01694,
     $     0.01679, 0.01666, 0.01651, 0.01638, 0.01625, 0.01614, 0.01601
     $     , 0.01588, 0.01577, 0.01565, 0.01554, 0.01544, 0.01533,
     $     0.01524, 0.01513, 0.01503, 0.01493, 0.01484, 0.01441, 0.01402
     $     , 0.01367, 0.01335, 0.01306, 0.01282, 0.01238, 0.01198,
     $     0.01167, 0.01139, 0.01115, 0.01095, 0.01076, 0.01059, 0.01043
     $     , 0.01028, 0.01014, 0.01001, 0.009872, 0.009739, 0.009603,
     $     0.009469, 0.009335, 0.009196, 0.009055, 0.008914, 0.008775,
     $     0.008637, 0.008495, 0.008355, 0.008216, 0.008074, 0.007932,
     $     0.007795, 0.007658, 0.00752, 0.007384, 0.00725, 0.00712,
     $     0.006991, 0.006864, 0.00674, 0.006616, 0.006499, 0.006384,
     $     0.006271, 0.006164, 0.00606, 0.005959, 0.005861, 0.005769,
     $     0.005683, 0.005602/


C     55133, JENDL evaluated data
      nbins(41) = 207
      data ( eg(41,i), i=1, 207) /9200000.0, 9400000.0, 9600000.0,
     $     9800000.0, 1.0D7, 1.02D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7,
     $     1.12D7, 1.14D7, 1.16D7, 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7
     $     , 1.28D7, 1.3D7, 1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7,
     $     1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7
     $     , 1.58D7, 1.6D7, 1.616D7, 1.62D7, 1.64D7, 1.66D7, 1.68D7,
     $     1.7D7, 1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7, 1.82D7, 1.84D7,
     $     1.86D7, 1.88D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7,
     $     2.02D7, 2.04D7, 2.06D7, 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7
     $     , 2.18D7, 2.2D7, 2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7,
     $     2.32D7, 2.34D7, 2.36D7, 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7
     $     , 2.48D7, 2.5D7, 2.52D7, 2.54D7, 2.545D7, 2.56D7, 2.58D7,
     $     2.6D7, 2.62D7, 2.64D7, 2.66D7, 2.68D7, 2.7D7, 2.72D7, 2.74D7,
     $     2.76D7, 2.78D7, 2.8D7, 2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7,
     $     2.92D7, 2.94D7, 2.96D7, 2.98D7, 3.0D7, 3.02D7, 3.04D7, 3.06D7
     $     , 3.08D7, 3.1D7, 3.12D7, 3.14D7, 3.16D7, 3.18D7, 3.2D7,
     $     3.22D7, 3.24D7, 3.26D7, 3.28D7, 3.3D7, 3.32D7, 3.34D7, 3.36D7
     $     , 3.38D7, 3.4D7, 3.42D7, 3.44D7, 3.46D7, 3.48D7, 3.5D7,
     $     3.52D7, 3.54D7, 3.56D7, 3.58D7, 3.6D7, 3.62D7, 3.64D7, 3.66D7
     $     , 3.68D7, 3.7D7, 3.72D7, 3.74D7, 3.76D7, 3.78D7, 3.8D7,
     $     3.82D7, 3.84D7, 3.86D7, 3.88D7, 3.9D7, 3.92D7, 3.94D7, 3.96D7
     $     , 3.98D7, 4.0D7, 4.2D7, 4.4D7, 4.6D7, 4.8D7, 5.0D7, 5.2D7,
     $     5.4D7, 5.6D7, 5.8D7, 6.0D7, 6.2D7, 6.4D7, 6.6D7, 6.8D7, 7.0D7
     $     , 7.2D7, 7.4D7, 7.6D7, 7.8D7, 8.0D7, 8.2D7, 8.4D7, 8.6D7,
     $     8.8D7, 9.0D7, 9.2D7, 9.4D7, 9.6D7, 9.8D7, 1.0D8, 1.02D8,
     $     1.04D8, 1.06D8, 1.08D8, 1.1D8, 1.12D8, 1.14D8, 1.16D8, 1.18D8
     $     , 1.2D8, 1.22D8, 1.24D8, 1.26D8, 1.28D8, 1.3D8, 1.32D8,
     $     1.34D8, 1.36D8, 1.38D8, 1.4D8/

      data ( gdr(41,i), i=1, 207) /0.0263, 0.03229, 0.04712, 0.04041,
     $     0.04841, 0.05247, 0.04885, 0.05044, 0.06052, 0.05964, 0.06147
     $     , 0.06888, 0.08744, 0.09651, 0.0921, 0.1089, 0.1113, 0.1239,
     $     0.1389, 0.1559, 0.175, 0.1959, 0.2174, 0.2381, 0.2566, 0.2731
     $     , 0.314, 0.2982, 0.3109, 0.3245, 0.3183, 0.331, 0.309, 0.3099
     $     , 0.3015, 0.2851, 0.2808, 0.2702, 0.2651, 0.2487, 0.2546,
     $     0.2367, 0.2305, 0.2026, 0.1916, 0.1839, 0.1786, 0.1731,
     $     0.1574, 0.1647, 0.1371, 0.1312, 0.1273, 0.124, 0.1207, 0.1172
     $     , 0.1132, 0.1089, 0.1043, 0.09966, 0.09513, 0.09098, 0.08744,
     $     0.08474, 0.08325, 0.08347, 0.08575, 0.08927, 0.09034, 0.08694
     $     , 0.08301, 0.08164, 0.08237, 0.08231, 0.07843, 0.07117,
     $     0.06336, 0.05668, 0.05143, 0.04735, 0.04413, 0.04153, 0.03936
     $     , 0.03893, 0.03752, 0.03592, 0.03449, 0.03323, 0.03208,
     $     0.03104, 0.03007, 0.02919, 0.02837, 0.02761, 0.0269, 0.02623,
     $     0.02561, 0.02502, 0.02448, 0.02396, 0.02347, 0.02301, 0.02257
     $     , 0.02215, 0.02176, 0.02139, 0.02103, 0.02071, 0.02038,
     $     0.02008, 0.01979, 0.01951, 0.01924, 0.01899, 0.01874, 0.01851
     $     , 0.01829, 0.01808, 0.01787, 0.01767, 0.01749, 0.0173,
     $     0.01713, 0.01697, 0.01681, 0.01666, 0.01651, 0.01637, 0.01624
     $     , 0.0161, 0.01597, 0.01586, 0.01574, 0.01562, 0.01552,
     $     0.01542, 0.01532, 0.01522, 0.01513, 0.01504, 0.01495, 0.01487
     $     , 0.01478, 0.01471, 0.01463, 0.01455, 0.01449, 0.01442,
     $     0.01435, 0.01429, 0.01423, 0.01416, 0.01411, 0.01405, 0.014,
     $     0.01395, 0.0139, 0.01338, 0.01309, 0.01289, 0.01274, 0.01259,
     $     0.01248, 0.01239, 0.01232, 0.01224, 0.01217, 0.01211, 0.01204
     $     , 0.01198, 0.01193, 0.01186, 0.01181, 0.01174, 0.01167,
     $     0.01159, 0.01151, 0.01143, 0.01134, 0.01123, 0.01112, 0.01099
     $     , 0.01086, 0.01072, 0.01058, 0.01043, 0.01028, 0.01012,
     $     0.009961, 0.009798, 0.009635, 0.009472, 0.009309, 0.009149,
     $     0.008996, 0.008839, 0.008685, 0.008533, 0.008385, 0.008242,
     $     0.008102, 0.007968, 0.007839, 0.007714, 0.007597, 0.007485,
     $     0.007379/


C     64152, JENDL evaluated data
      nbins(42) = 223
      data ( eg(42,i), i=1, 223) /8592400.0, 8600000.0, 8800000.0,
     $     9000000.0, 9200000.0, 9400000.0, 9600000.0, 9800000.0, 1.0D7,
     $     1.02D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7
     $     , 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7,
     $     1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7
     $     , 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7,
     $     1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7
     $     , 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7,
     $     1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7
     $     , 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7,
     $     2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7
     $     , 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7,
     $     2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7
     $     , 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7,
     $     2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7
     $     , 2.98D7, 3.0D7, 3.05D7, 3.1D7, 3.15D7, 3.2D7, 3.25D7, 3.3D7,
     $     3.35D7, 3.4D7, 3.45D7, 3.5D7, 3.55D7, 3.6D7, 3.65D7, 3.7D7,
     $     3.75D7, 3.8D7, 3.85D7, 3.9D7, 3.95D7, 4.0D7, 4.05D7, 4.1D7,
     $     4.15D7, 4.2D7, 4.25D7, 4.3D7, 4.35D7, 4.4D7, 4.45D7, 4.5D7,
     $     4.55D7, 4.6D7, 4.65D7, 4.7D7, 4.75D7, 4.8D7, 4.85D7, 4.9D7,
     $     4.95D7, 5.0D7, 5.05D7, 5.1D7, 5.15D7, 5.2D7, 5.25D7, 5.3D7,
     $     5.35D7, 5.4D7, 5.45D7, 5.5D7, 5.55D7, 5.6D7, 5.65D7, 5.7D7,
     $     5.75D7, 5.8D7, 5.85D7, 5.9D7, 5.95D7, 6.0D7, 6.05D7, 6.1D7,
     $     6.15D7, 6.2D7, 6.25D7, 6.3D7, 6.35D7, 6.4D7, 6.45D7, 6.5D7,
     $     6.55D7, 6.6D7, 6.65D7, 6.7D7, 6.75D7, 6.8D7, 6.85D7, 6.9D7,
     $     6.95D7, 7.0D7, 7.05D7, 7.1D7, 7.15D7, 7.2D7, 7.25D7, 7.3D7,
     $     7.35D7, 7.4D7, 7.45D7, 7.5D7, 7.55D7, 7.6D7, 7.65D7, 7.7D7,
     $     7.75D7, 7.8D7, 7.85D7, 7.9D7, 7.95D7, 8.0D7, 8.5D7, 9.0D7,
     $     9.5D7, 1.0D8, 1.05D8, 1.1D8, 1.15D8, 1.2D8, 1.25D8, 1.3D8,
     $     1.35D8, 1.4D8, 1.45D8, 1.5D8/

      data ( gdr(42,i), i=1, 223) /0.0, 5.4231D-5, 0.022668, 0.032179,
     $     0.039199, 0.044987, 0.050207, 0.054576, 0.060159, 0.067621,
     $     0.077697, 0.091045, 0.10505, 0.12293, 0.14489, 0.17061,
     $     0.19889, 0.2273, 0.25278, 0.27243, 0.28467, 0.28966, 0.28662,
     $     0.27948, 0.27039, 0.26155, 0.25486, 0.24955, 0.2449, 0.23989,
     $     0.23378, 0.22648, 0.21781, 0.20895, 0.20507, 0.22027, 0.22749
     $     , 0.23092, 0.23218, 0.23146, 0.22952, 0.22764, 0.22637,
     $     0.22636, 0.227, 0.2278, 0.22816, 0.22796, 0.22611, 0.22189,
     $     0.21504, 0.20575, 0.19441, 0.18178, 0.16847, 0.15522, 0.14267
     $     , 0.1311, 0.12029, 0.11061, 0.10211, 0.094692, 0.088272,
     $     0.082748, 0.078002, 0.073965, 0.070595, 0.067835, 0.065639,
     $     0.06397, 0.062761, 0.061924, 0.061532, 0.061559, 0.061978,
     $     0.062758, 0.063859, 0.065212, 0.066861, 0.068676, 0.070478,
     $     0.072103, 0.073395, 0.074111, 0.074215, 0.073695, 0.07255,
     $     0.07082, 0.068555, 0.065874, 0.062896, 0.059741, 0.056516,
     $     0.053312, 0.050199, 0.047215, 0.044376, 0.041728, 0.039275,
     $     0.037016, 0.034945, 0.033051, 0.03132, 0.029741, 0.028299,
     $     0.026982, 0.02578, 0.024679, 0.023664, 0.021514, 0.019775,
     $     0.018359, 0.017197, 0.016235, 0.015431, 0.01476, 0.014197,
     $     0.013724, 0.013323, 0.012968, 0.012665, 0.012408, 0.012189,
     $     0.012002, 0.01184, 0.011704, 0.011589, 0.011493, 0.011412,
     $     0.011365, 0.01133, 0.011305, 0.01129, 0.011282, 0.01128,
     $     0.011284, 0.011293, 0.011307, 0.011325, 0.011345, 0.011369,
     $     0.011395, 0.011423, 0.011453, 0.011481, 0.011511, 0.011543,
     $     0.011575, 0.011608, 0.011645, 0.011682, 0.01172, 0.011758,
     $     0.011796, 0.011833, 0.011871, 0.011908, 0.011944, 0.011981,
     $     0.012017, 0.012052, 0.012086, 0.01212, 0.012152, 0.012181,
     $     0.012208, 0.012235, 0.012262, 0.012287, 0.012311, 0.012335,
     $     0.012357, 0.012379, 0.012399, 0.012419, 0.012438, 0.012456,
     $     0.012472, 0.012488, 0.012503, 0.012517, 0.01253, 0.012542,
     $     0.012552, 0.01256, 0.012567, 0.012573, 0.012578, 0.012582,
     $     0.012585, 0.012587, 0.012589, 0.012589, 0.012588, 0.012587,
     $     0.012584, 0.012581, 0.012577, 0.012572, 0.012566, 0.01256,
     $     0.012552, 0.012544, 0.012534, 0.012523, 0.012512, 0.012499,
     $     0.012486, 0.012472, 0.012289, 0.012043, 0.011758, 0.011431,
     $     0.011075, 0.010702, 0.010323, 0.0099474, 0.0095849, 0.0092431
     $     , 0.0089304, 0.0086558, 0.0084268, 0.0082503/


C     64154, JENDL evaluated data
      nbins(43) = 221
      data ( eg(43,i), i=1, 221) /8894400.0, 9000000.0, 9200000.0,
     $     9400000.0, 9600000.0, 9800000.0, 1.0D7, 1.02D7, 1.04D7,
     $     1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7, 1.18D7, 1.2D7,
     $     1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7, 1.32D7, 1.34D7, 1.36D7
     $     , 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7,
     $     1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7, 1.62D7, 1.64D7, 1.66D7
     $     , 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7,
     $     1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7
     $     , 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7, 2.08D7, 2.1D7,
     $     2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7, 2.22D7, 2.24D7, 2.26D7
     $     , 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7, 2.38D7, 2.4D7,
     $     2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7, 2.52D7, 2.54D7, 2.56D7
     $     , 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7, 2.68D7, 2.7D7,
     $     2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7, 2.82D7, 2.84D7, 2.86D7
     $     , 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7, 2.98D7, 3.0D7,
     $     3.05D7, 3.1D7, 3.15D7, 3.2D7, 3.25D7, 3.3D7, 3.35D7, 3.4D7,
     $     3.45D7, 3.5D7, 3.55D7, 3.6D7, 3.65D7, 3.7D7, 3.75D7, 3.8D7,
     $     3.85D7, 3.9D7, 3.95D7, 4.0D7, 4.05D7, 4.1D7, 4.15D7, 4.2D7,
     $     4.25D7, 4.3D7, 4.35D7, 4.4D7, 4.45D7, 4.5D7, 4.55D7, 4.6D7,
     $     4.65D7, 4.7D7, 4.75D7, 4.8D7, 4.85D7, 4.9D7, 4.95D7, 5.0D7,
     $     5.05D7, 5.1D7, 5.15D7, 5.2D7, 5.25D7, 5.3D7, 5.35D7, 5.4D7,
     $     5.45D7, 5.5D7, 5.55D7, 5.6D7, 5.65D7, 5.7D7, 5.75D7, 5.8D7,
     $     5.85D7, 5.9D7, 5.95D7, 6.0D7, 6.05D7, 6.1D7, 6.15D7, 6.2D7,
     $     6.25D7, 6.3D7, 6.35D7, 6.4D7, 6.45D7, 6.5D7, 6.55D7, 6.6D7,
     $     6.65D7, 6.7D7, 6.75D7, 6.8D7, 6.85D7, 6.9D7, 6.95D7, 7.0D7,
     $     7.05D7, 7.1D7, 7.15D7, 7.2D7, 7.25D7, 7.3D7, 7.35D7, 7.4D7,
     $     7.45D7, 7.5D7, 7.55D7, 7.6D7, 7.65D7, 7.7D7, 7.75D7, 7.8D7,
     $     7.85D7, 7.9D7, 7.95D7, 8.0D7, 8.5D7, 9.0D7, 9.5D7, 1.0D8,
     $     1.05D8, 1.1D8, 1.15D8, 1.2D8, 1.25D8, 1.3D8, 1.35D8, 1.4D8,
     $     1.45D8, 1.5D8/

      data ( gdr(43,i), i=1, 221) /0.0, 2.1209D-4, 0.002331, 0.0055284,
     $     0.010073, 0.01622, 0.023947, 0.033611, 0.046522, 0.063503,
     $     0.085371, 0.11049, 0.13839, 0.17015, 0.20405, 0.23726, 0.2661
     $     , 0.28757, 0.30047, 0.30561, 0.3045, 0.29789, 0.28839,
     $     0.27901, 0.27057, 0.26467, 0.25929, 0.25309, 0.2455, 0.23664,
     $     0.22697, 0.21695, 0.20846, 0.21962, 0.22176, 0.22352, 0.22443
     $     , 0.22427, 0.22331, 0.22304, 0.224, 0.22642, 0.22978, 0.233,
     $     0.23611, 0.23672, 0.23486, 0.23013, 0.22252, 0.21231, 0.2001,
     $     0.18672, 0.17289, 0.15923, 0.1464, 0.13466, 0.12409, 0.11456,
     $     0.10592, 0.098421, 0.091947, 0.086398, 0.081675, 0.077687,
     $     0.074375, 0.071703, 0.069624, 0.068104, 0.067112, 0.066625,
     $     0.066476, 0.066751, 0.067442, 0.068504, 0.06989, 0.071529,
     $     0.073475, 0.075481, 0.077374, 0.078977, 0.080123, 0.080682,
     $     0.080507, 0.07955, 0.077927, 0.075693, 0.072937, 0.069797,
     $     0.066404, 0.062884, 0.059346, 0.055877, 0.052539, 0.049375,
     $     0.046409, 0.04362, 0.041038, 0.038668, 0.036497, 0.034513,
     $     0.032701, 0.031048, 0.029538, 0.028158, 0.026897, 0.025743,
     $     0.024684, 0.022421, 0.020575, 0.019068, 0.017828, 0.016799,
     $     0.015935, 0.015208, 0.014596, 0.014078, 0.01364, 0.013267,
     $     0.012948, 0.012678, 0.012448, 0.012253, 0.012086, 0.011943,
     $     0.011822, 0.011723, 0.011641, 0.011584, 0.01154, 0.011507,
     $     0.011485, 0.011471, 0.011463, 0.011461, 0.011465, 0.011474,
     $     0.011489, 0.011505, 0.011525, 0.011547, 0.011573, 0.0116,
     $     0.011627, 0.011654, 0.011683, 0.011713, 0.011744, 0.01178,
     $     0.011817, 0.011854, 0.011892, 0.011929, 0.011966, 0.012004,
     $     0.012041, 0.012077, 0.012113, 0.012149, 0.012185, 0.012219,
     $     0.012253, 0.012286, 0.012316, 0.012344, 0.012371, 0.012397,
     $     0.012423, 0.012447, 0.01247, 0.012492, 0.012513, 0.012533,
     $     0.012553, 0.012571, 0.012588, 0.012605, 0.01262, 0.012635,
     $     0.012648, 0.012661, 0.012672, 0.012683, 0.012691, 0.012698,
     $     0.012704, 0.012709, 0.012713, 0.012716, 0.012718, 0.012719,
     $     0.012719, 0.012719, 0.012717, 0.012715, 0.012712, 0.012708,
     $     0.012703, 0.012697, 0.01269, 0.012683, 0.012674, 0.012665,
     $     0.012654, 0.012643, 0.01263, 0.012617, 0.012603, 0.012414,
     $     0.012163, 0.011874, 0.011545, 0.011186, 0.010809, 0.010425,
     $     0.010046, 0.009679, 0.0093335, 0.0090181, 0.0087413,
     $     0.0085103, 0.0083324/


C     64155, JENDL evaluated data
      nbins(44) = 233
      data ( eg(44,i), i=1, 233) /6435400.0, 6600000.0, 6800000.0,
     $     7000000.0, 7200000.0, 7400000.0, 7600000.0, 7800000.0,
     $     8000000.0, 8200000.0, 8400000.0, 8600000.0, 8800000.0,
     $     9000000.0, 9200000.0, 9400000.0, 9600000.0, 9800000.0, 1.0D7,
     $     1.02D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7
     $     , 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7,
     $     1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7
     $     , 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7,
     $     1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7
     $     , 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7,
     $     1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7
     $     , 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7,
     $     2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7
     $     , 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7,
     $     2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7
     $     , 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7,
     $     2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7
     $     , 2.98D7, 3.0D7, 3.05D7, 3.1D7, 3.15D7, 3.2D7, 3.25D7, 3.3D7,
     $     3.35D7, 3.4D7, 3.45D7, 3.5D7, 3.55D7, 3.6D7, 3.65D7, 3.7D7,
     $     3.75D7, 3.8D7, 3.85D7, 3.9D7, 3.95D7, 4.0D7, 4.05D7, 4.1D7,
     $     4.15D7, 4.2D7, 4.25D7, 4.3D7, 4.35D7, 4.4D7, 4.45D7, 4.5D7,
     $     4.55D7, 4.6D7, 4.65D7, 4.7D7, 4.75D7, 4.8D7, 4.85D7, 4.9D7,
     $     4.95D7, 5.0D7, 5.05D7, 5.1D7, 5.15D7, 5.2D7, 5.25D7, 5.3D7,
     $     5.35D7, 5.4D7, 5.45D7, 5.5D7, 5.55D7, 5.6D7, 5.65D7, 5.7D7,
     $     5.75D7, 5.8D7, 5.85D7, 5.9D7, 5.95D7, 6.0D7, 6.05D7, 6.1D7,
     $     6.15D7, 6.2D7, 6.25D7, 6.3D7, 6.35D7, 6.4D7, 6.45D7, 6.5D7,
     $     6.55D7, 6.6D7, 6.65D7, 6.7D7, 6.75D7, 6.8D7, 6.85D7, 6.9D7,
     $     6.95D7, 7.0D7, 7.05D7, 7.1D7, 7.15D7, 7.2D7, 7.25D7, 7.3D7,
     $     7.35D7, 7.4D7, 7.45D7, 7.5D7, 7.55D7, 7.6D7, 7.65D7, 7.7D7,
     $     7.75D7, 7.8D7, 7.85D7, 7.9D7, 7.95D7, 8.0D7, 8.5D7, 9.0D7,
     $     9.5D7, 1.0D8, 1.05D8, 1.1D8, 1.15D8, 1.2D8, 1.25D8, 1.3D8,
     $     1.35D8, 1.4D8, 1.45D8, 1.5D8/

      data ( gdr(44,i), i=1, 233) /0.0, 7.71D-4, 0.0035748, 0.0074702,
     $     0.012455, 0.018231, 0.023898, 0.02953, 0.03484, 0.039575,
     $     0.043767, 0.046649, 0.04918, 0.051799, 0.054788, 0.058473,
     $     0.063211, 0.069382, 0.077388, 0.087658, 0.10064, 0.11534,
     $     0.13293, 0.15383, 0.17785, 0.20421, 0.23134, 0.257, 0.27833,
     $     0.29302, 0.30012, 0.29947, 0.29324, 0.28374, 0.27315, 0.26303
     $     , 0.25396, 0.24767, 0.24187, 0.23534, 0.22758, 0.21846,
     $     0.20874, 0.1992, 0.19041, 0.18579, 0.19373, 0.1961, 0.19744,
     $     0.19803, 0.19816, 0.19828, 0.19959, 0.20231, 0.20623, 0.21032
     $     , 0.21337, 0.21381, 0.21194, 0.20745, 0.20038, 0.19099,
     $     0.17991, 0.16788, 0.15563, 0.14366, 0.13232, 0.12197, 0.11269
     $     , 0.10429, 0.096664, 0.090055, 0.084365, 0.079506, 0.075396,
     $     0.071956, 0.069111, 0.066894, 0.065269, 0.064197, 0.063653,
     $     0.063613, 0.063888, 0.064619, 0.065797, 0.06738, 0.069288,
     $     0.071352, 0.073536, 0.075695, 0.077671, 0.079291, 0.080397,
     $     0.080862, 0.080573, 0.079499, 0.077739, 0.075339, 0.072443,
     $     0.069184, 0.065699, 0.062113, 0.058534, 0.055043, 0.0517,
     $     0.048543, 0.045592, 0.042828, 0.040281, 0.037946, 0.035811,
     $     0.033861, 0.032083, 0.030462, 0.028984, 0.027634, 0.026401,
     $     0.025272, 0.024239, 0.02203, 0.020241, 0.018785, 0.017589,
     $     0.016602, 0.015777, 0.015087, 0.014508, 0.014022, 0.013611,
     $     0.013247, 0.012936, 0.012672, 0.012447, 0.012256, 0.012092,
     $     0.011953, 0.011836, 0.011738, 0.011658, 0.011604, 0.011562,
     $     0.011532, 0.011511, 0.011499, 0.011492, 0.011491, 0.011497,
     $     0.011507, 0.011522, 0.011539, 0.011559, 0.011581, 0.011607,
     $     0.011635, 0.011661, 0.011689, 0.011718, 0.011747, 0.011778,
     $     0.011814, 0.011851, 0.011888, 0.011925, 0.011962, 0.011999,
     $     0.012036, 0.012073, 0.012109, 0.012145, 0.01218, 0.012215,
     $     0.012249, 0.012283, 0.012315, 0.012345, 0.012372, 0.012399,
     $     0.012425, 0.012449, 0.012475, 0.012499, 0.012522, 0.012544,
     $     0.012566, 0.012586, 0.012606, 0.012624, 0.012642, 0.012658,
     $     0.012673, 0.012688, 0.012701, 0.012714, 0.012725, 0.012734,
     $     0.012742, 0.012748, 0.012754, 0.012758, 0.012761, 0.012763,
     $     0.012764, 0.012764, 0.012764, 0.012762, 0.012759, 0.012756,
     $     0.012752, 0.012746, 0.01274, 0.012734, 0.012726, 0.012717,
     $     0.012708, 0.012697, 0.012685, 0.012672, 0.012659, 0.012644,
     $     0.012459, 0.01221, 0.011924, 0.011596, 0.011235, 0.010856,
     $     0.01047, 0.010089, 0.009721, 0.0093744, 0.0090573, 0.0087788,
     $     0.0085466, 0.0083676/


C     64156, JENDL evaluated data
      nbins(45) = 223
      data ( eg(45,i), i=1, 223) /8536400.0, 8600000.0, 8800000.0,
     $     9000000.0, 9200000.0, 9400000.0, 9600000.0, 9800000.0, 1.0D7,
     $     1.02D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7
     $     , 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7,
     $     1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7
     $     , 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7,
     $     1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7
     $     , 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7,
     $     1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7
     $     , 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7,
     $     2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7
     $     , 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7,
     $     2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7
     $     , 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7,
     $     2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7
     $     , 2.98D7, 3.0D7, 3.05D7, 3.1D7, 3.15D7, 3.2D7, 3.25D7, 3.3D7,
     $     3.35D7, 3.4D7, 3.45D7, 3.5D7, 3.55D7, 3.6D7, 3.65D7, 3.7D7,
     $     3.75D7, 3.8D7, 3.85D7, 3.9D7, 3.95D7, 4.0D7, 4.05D7, 4.1D7,
     $     4.15D7, 4.2D7, 4.25D7, 4.3D7, 4.35D7, 4.4D7, 4.45D7, 4.5D7,
     $     4.55D7, 4.6D7, 4.65D7, 4.7D7, 4.75D7, 4.8D7, 4.85D7, 4.9D7,
     $     4.95D7, 5.0D7, 5.05D7, 5.1D7, 5.15D7, 5.2D7, 5.25D7, 5.3D7,
     $     5.35D7, 5.4D7, 5.45D7, 5.5D7, 5.55D7, 5.6D7, 5.65D7, 5.7D7,
     $     5.75D7, 5.8D7, 5.85D7, 5.9D7, 5.95D7, 6.0D7, 6.05D7, 6.1D7,
     $     6.15D7, 6.2D7, 6.25D7, 6.3D7, 6.35D7, 6.4D7, 6.45D7, 6.5D7,
     $     6.55D7, 6.6D7, 6.65D7, 6.7D7, 6.75D7, 6.8D7, 6.85D7, 6.9D7,
     $     6.95D7, 7.0D7, 7.05D7, 7.1D7, 7.15D7, 7.2D7, 7.25D7, 7.3D7,
     $     7.35D7, 7.4D7, 7.45D7, 7.5D7, 7.55D7, 7.6D7, 7.65D7, 7.7D7,
     $     7.75D7, 7.8D7, 7.85D7, 7.9D7, 7.95D7, 8.0D7, 8.5D7, 9.0D7,
     $     9.5D7, 1.0D8, 1.05D8, 1.1D8, 1.15D8, 1.2D8, 1.25D8, 1.3D8,
     $     1.35D8, 1.4D8, 1.45D8, 1.5D8/

      data ( gdr(45,i), i=1, 223) /0.0, 0.005709, 0.027797, 0.037345,
     $     0.044615, 0.050775, 0.056265, 0.061507, 0.068363, 0.077596,
     $     0.090039, 0.10537, 0.12242, 0.14379, 0.16945, 0.19852,
     $     0.22903, 0.25796, 0.28198, 0.29875, 0.30753, 0.30861, 0.30269
     $     , 0.29362, 0.28358, 0.27403, 0.26543, 0.25726, 0.24843,
     $     0.24063, 0.23156, 0.22121, 0.21002, 0.1993, 0.21136, 0.21705,
     $     0.21871, 0.21839, 0.21687, 0.21603, 0.21541, 0.21545, 0.21659
     $     , 0.21886, 0.22092, 0.22352, 0.22564, 0.22615, 0.22422,
     $     0.21841, 0.20999, 0.19934, 0.18713, 0.17411, 0.1608, 0.14806,
     $     0.13622, 0.12544, 0.1158, 0.10728, 0.099802, 0.093147,
     $     0.087241, 0.08218, 0.077843, 0.074207, 0.071216, 0.06882,
     $     0.066979, 0.065661, 0.064839, 0.064493, 0.0646, 0.065135,
     $     0.065891, 0.066978, 0.068336, 0.069878, 0.071554, 0.07334,
     $     0.074923, 0.076153, 0.076893, 0.077036, 0.076512, 0.075348,
     $     0.073524, 0.071107, 0.068266, 0.065129, 0.061823, 0.058458,
     $     0.055125, 0.051894, 0.048812, 0.04591, 0.043203, 0.040696,
     $     0.038384, 0.036235, 0.034267, 0.03247, 0.030829, 0.02933,
     $     0.027962, 0.026711, 0.025566, 0.024517, 0.023555, 0.021507,
     $     0.019824, 0.018448, 0.017313, 0.016372, 0.015587, 0.014921,
     $     0.01436, 0.013886, 0.013485, 0.013141, 0.012844, 0.012591,
     $     0.012376, 0.012194, 0.01204, 0.011906, 0.011794, 0.011701,
     $     0.011624, 0.011578, 0.011541, 0.011514, 0.011497, 0.011489,
     $     0.011487, 0.01149, 0.011498, 0.011511, 0.011528, 0.011549,
     $     0.011571, 0.011596, 0.011624, 0.011654, 0.011685, 0.011714,
     $     0.011745, 0.011776, 0.011809, 0.011847, 0.011886, 0.011925,
     $     0.011964, 0.012004, 0.012043, 0.012082, 0.01212, 0.012158,
     $     0.012196, 0.012233, 0.01227, 0.012306, 0.012341, 0.012375,
     $     0.012409, 0.012438, 0.012466, 0.012493, 0.012519, 0.012544,
     $     0.012567, 0.01259, 0.012611, 0.012632, 0.012651, 0.01267,
     $     0.012687, 0.012704, 0.01272, 0.012734, 0.012748, 0.012761,
     $     0.012773, 0.012783, 0.012793, 0.0128, 0.012806, 0.012811,
     $     0.012815, 0.012817, 0.012819, 0.01282, 0.01282, 0.012819,
     $     0.012817, 0.012815, 0.012811, 0.012806, 0.012801, 0.012795,
     $     0.012788, 0.01278, 0.012771, 0.012762, 0.012751, 0.012739,
     $     0.012726, 0.012713, 0.012698, 0.012515, 0.012267, 0.011976,
     $     0.011644, 0.011282, 0.010902, 0.010515, 0.010133, 0.0097635,
     $     0.0094157, 0.0090973, 0.0088178, 0.0085846, 0.008405/


C     64157, JENDL evaluated data
      nbins(46) = 234
      data ( eg(46,i), i=1, 234) /6359400.0, 6400000.0, 6600000.0,
     $     6800000.0, 7000000.0, 7200000.0, 7400000.0, 7600000.0,
     $     7800000.0, 8000000.0, 8200000.0, 8400000.0, 8600000.0,
     $     8800000.0, 9000000.0, 9200000.0, 9400000.0, 9600000.0,
     $     9800000.0, 1.0D7, 1.02D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7,
     $     1.12D7, 1.14D7, 1.16D7, 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7
     $     , 1.28D7, 1.3D7, 1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7,
     $     1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7
     $     , 1.58D7, 1.6D7, 1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7,
     $     1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7
     $     , 1.88D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7,
     $     2.02D7, 2.04D7, 2.06D7, 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7
     $     , 2.18D7, 2.2D7, 2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7,
     $     2.32D7, 2.34D7, 2.36D7, 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7
     $     , 2.48D7, 2.5D7, 2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7,
     $     2.62D7, 2.64D7, 2.66D7, 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7
     $     , 2.78D7, 2.8D7, 2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7,
     $     2.92D7, 2.94D7, 2.96D7, 2.98D7, 3.0D7, 3.05D7, 3.1D7, 3.15D7,
     $     3.2D7, 3.25D7, 3.3D7, 3.35D7, 3.4D7, 3.45D7, 3.5D7, 3.55D7,
     $     3.6D7, 3.65D7, 3.7D7, 3.75D7, 3.8D7, 3.85D7, 3.9D7, 3.95D7,
     $     4.0D7, 4.05D7, 4.1D7, 4.15D7, 4.2D7, 4.25D7, 4.3D7, 4.35D7,
     $     4.4D7, 4.45D7, 4.5D7, 4.55D7, 4.6D7, 4.65D7, 4.7D7, 4.75D7,
     $     4.8D7, 4.85D7, 4.9D7, 4.95D7, 5.0D7, 5.05D7, 5.1D7, 5.15D7,
     $     5.2D7, 5.25D7, 5.3D7, 5.35D7, 5.4D7, 5.45D7, 5.5D7, 5.55D7,
     $     5.6D7, 5.65D7, 5.7D7, 5.75D7, 5.8D7, 5.85D7, 5.9D7, 5.95D7,
     $     6.0D7, 6.05D7, 6.1D7, 6.15D7, 6.2D7, 6.25D7, 6.3D7, 6.35D7,
     $     6.4D7, 6.45D7, 6.5D7, 6.55D7, 6.6D7, 6.65D7, 6.7D7, 6.75D7,
     $     6.8D7, 6.85D7, 6.9D7, 6.95D7, 7.0D7, 7.05D7, 7.1D7, 7.15D7,
     $     7.2D7, 7.25D7, 7.3D7, 7.35D7, 7.4D7, 7.45D7, 7.5D7, 7.55D7,
     $     7.6D7, 7.65D7, 7.7D7, 7.75D7, 7.8D7, 7.85D7, 7.9D7, 7.95D7,
     $     8.0D7, 8.5D7, 9.0D7, 9.5D7, 1.0D8, 1.05D8, 1.1D8, 1.15D8,
     $     1.2D8, 1.25D8, 1.3D8, 1.35D8, 1.4D8, 1.45D8, 1.5D8/

      data ( gdr(46,i), i=1, 234) /0.0, 1.3073D-5, 0.0017122, 0.0046946,
     $     0.0088451, 0.013997, 0.019728, 0.025298, 0.030896, 0.03611,
     $     0.040749, 0.044657, 0.04732, 0.049895, 0.05262, 0.05579,
     $     0.059739, 0.064832, 0.071455, 0.080016, 0.090952, 0.10438,
     $     0.11958, 0.13798, 0.15968, 0.18435, 0.211, 0.23786, 0.26238,
     $     0.28186, 0.29446, 0.29936, 0.29691, 0.28959, 0.27961, 0.26903
     $     , 0.25913, 0.25021, 0.24169, 0.2325, 0.22431, 0.21519,
     $     0.20485, 0.19424, 0.1886, 0.20313, 0.20851, 0.21032, 0.21031,
     $     0.20871, 0.20804, 0.20777, 0.20807, 0.20946, 0.21124, 0.21329
     $     , 0.21569, 0.21746, 0.21756, 0.21524, 0.2093, 0.20076,
     $     0.19019, 0.17828, 0.1656, 0.15289, 0.1408, 0.12961, 0.11948,
     $     0.11044, 0.10246, 0.095483, 0.08931, 0.083787, 0.079054,
     $     0.075039, 0.071688, 0.069041, 0.067025, 0.065605, 0.064751,
     $     0.064429, 0.064605, 0.065305, 0.066505, 0.068042, 0.069915,
     $     0.072004, 0.074289, 0.076622, 0.078839, 0.080761, 0.082208,
     $     0.083032, 0.083129, 0.082472, 0.081056, 0.078901, 0.076127,
     $     0.07291, 0.069393, 0.065713, 0.06199, 0.058323, 0.054782,
     $     0.051419, 0.048257, 0.045308, 0.042584, 0.040082, 0.03778,
     $     0.035674, 0.033754, 0.032003, 0.030407, 0.02895, 0.02762,
     $     0.026402, 0.025288, 0.024268, 0.022085, 0.020305, 0.01885,
     $     0.017652, 0.016661, 0.015836, 0.015141, 0.014557, 0.014064,
     $     0.013648, 0.013287, 0.012978, 0.012714, 0.01249, 0.0123,
     $     0.01214, 0.012002, 0.011886, 0.011789, 0.011709, 0.011659,
     $     0.011619, 0.011589, 0.01157, 0.011559, 0.011555, 0.011556,
     $     0.011563, 0.011574, 0.01159, 0.01161, 0.011632, 0.011656,
     $     0.011683, 0.011712, 0.011743, 0.011772, 0.011802, 0.011833,
     $     0.011865, 0.011902, 0.01194, 0.011978, 0.012016, 0.012054,
     $     0.012092, 0.012129, 0.012167, 0.012204, 0.012241, 0.012277,
     $     0.012312, 0.012347, 0.012381, 0.012415, 0.012448, 0.012476,
     $     0.012504, 0.01253, 0.012556, 0.012581, 0.012605, 0.012629,
     $     0.012652, 0.012673, 0.012694, 0.012713, 0.012732, 0.01275,
     $     0.012766, 0.012782, 0.012797, 0.01281, 0.012823, 0.012834,
     $     0.012845, 0.012852, 0.012859, 0.012865, 0.012869, 0.012872,
     $     0.012875, 0.012876, 0.012876, 0.012876, 0.012874, 0.012872,
     $     0.012868, 0.012864, 0.012859, 0.012853, 0.012846, 0.012839,
     $     0.01283, 0.012821, 0.012811, 0.012799, 0.012786, 0.012773,
     $     0.012758, 0.012567, 0.012314, 0.012023, 0.01169, 0.01133,
     $     0.010951, 0.010562, 0.010179, 0.0098068, 0.0094568, 0.0091364
     $     , 0.0088552, 0.0086205, 0.0084397/


C     64158, JENDL evaluated data
      nbins(47) = 226
      data ( eg(47,i), i=1, 226) /7938400.0, 8000000.0, 8200000.0,
     $     8400000.0, 8600000.0, 8800000.0, 9000000.0, 9200000.0,
     $     9400000.0, 9600000.0, 9800000.0, 1.0D7, 1.02D7, 1.04D7,
     $     1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7, 1.18D7, 1.2D7,
     $     1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7, 1.32D7, 1.34D7, 1.36D7
     $     , 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7,
     $     1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7, 1.62D7, 1.64D7, 1.66D7
     $     , 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7,
     $     1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7
     $     , 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7, 2.08D7, 2.1D7,
     $     2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7, 2.22D7, 2.24D7, 2.26D7
     $     , 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7, 2.38D7, 2.4D7,
     $     2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7, 2.52D7, 2.54D7, 2.56D7
     $     , 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7, 2.68D7, 2.7D7,
     $     2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7, 2.82D7, 2.84D7, 2.86D7
     $     , 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7, 2.98D7, 3.0D7,
     $     3.05D7, 3.1D7, 3.15D7, 3.2D7, 3.25D7, 3.3D7, 3.35D7, 3.4D7,
     $     3.45D7, 3.5D7, 3.55D7, 3.6D7, 3.65D7, 3.7D7, 3.75D7, 3.8D7,
     $     3.85D7, 3.9D7, 3.95D7, 4.0D7, 4.05D7, 4.1D7, 4.15D7, 4.2D7,
     $     4.25D7, 4.3D7, 4.35D7, 4.4D7, 4.45D7, 4.5D7, 4.55D7, 4.6D7,
     $     4.65D7, 4.7D7, 4.75D7, 4.8D7, 4.85D7, 4.9D7, 4.95D7, 5.0D7,
     $     5.05D7, 5.1D7, 5.15D7, 5.2D7, 5.25D7, 5.3D7, 5.35D7, 5.4D7,
     $     5.45D7, 5.5D7, 5.55D7, 5.6D7, 5.65D7, 5.7D7, 5.75D7, 5.8D7,
     $     5.85D7, 5.9D7, 5.95D7, 6.0D7, 6.05D7, 6.1D7, 6.15D7, 6.2D7,
     $     6.25D7, 6.3D7, 6.35D7, 6.4D7, 6.45D7, 6.5D7, 6.55D7, 6.6D7,
     $     6.65D7, 6.7D7, 6.75D7, 6.8D7, 6.85D7, 6.9D7, 6.95D7, 7.0D7,
     $     7.05D7, 7.1D7, 7.15D7, 7.2D7, 7.25D7, 7.3D7, 7.35D7, 7.4D7,
     $     7.45D7, 7.5D7, 7.55D7, 7.6D7, 7.65D7, 7.7D7, 7.75D7, 7.8D7,
     $     7.85D7, 7.9D7, 7.95D7, 8.0D7, 8.5D7, 9.0D7, 9.5D7, 1.0D8,
     $     1.05D8, 1.1D8, 1.15D8, 1.2D8, 1.25D8, 1.3D8, 1.35D8, 1.4D8,
     $     1.45D8, 1.5D8/

      data ( gdr(47,i), i=1, 226) /0.0, 0.001993, 0.022592, 0.031735,
     $     0.0381, 0.043171, 0.047237, 0.050492, 0.054354, 0.059329,
     $     0.065985, 0.074173, 0.083172, 0.095107, 0.11054, 0.12994,
     $     0.15351, 0.18092, 0.21095, 0.24136, 0.26907, 0.29029, 0.30232
     $     , 0.30663, 0.30444, 0.29755, 0.28812, 0.27815, 0.26885,
     $     0.26036, 0.25195, 0.24231, 0.2304, 0.21774, 0.21735, 0.22074,
     $     0.22421, 0.22632, 0.22618, 0.22439, 0.22215, 0.21958, 0.2169,
     $     0.21374, 0.21143, 0.21151, 0.2129, 0.21525, 0.21776, 0.21938,
     $     0.2191, 0.21626, 0.21068, 0.2023, 0.1913, 0.17861, 0.16554,
     $     0.15272, 0.1406, 0.12943, 0.11933, 0.11032, 0.10239, 0.095412
     $     , 0.089325, 0.08409, 0.079635, 0.075634, 0.07228, 0.069542,
     $     0.067376, 0.065744, 0.064617, 0.063964, 0.063765, 0.064005,
     $     0.064777, 0.06626, 0.068208, 0.070448, 0.072791, 0.075108,
     $     0.07729, 0.079153, 0.080514, 0.081236, 0.081151, 0.080326,
     $     0.07881, 0.076662, 0.073981, 0.070891, 0.067464, 0.063897,
     $     0.060308, 0.056786, 0.053372, 0.050131, 0.047096, 0.044276,
     $     0.041673, 0.03928, 0.037086, 0.03508, 0.033241, 0.031554,
     $     0.03001, 0.028599, 0.027307, 0.026127, 0.025047, 0.024055,
     $     0.02192, 0.020183, 0.018748, 0.017563, 0.016578, 0.015754,
     $     0.015061, 0.014471, 0.013971, 0.013547, 0.013201, 0.012908,
     $     0.012657, 0.012444, 0.012265, 0.012114, 0.011988, 0.011881,
     $     0.011794, 0.011722, 0.011673, 0.011636, 0.011608, 0.011589,
     $     0.011578, 0.011576, 0.011579, 0.011587, 0.011599, 0.011615,
     $     0.011636, 0.01166, 0.011685, 0.011712, 0.011742, 0.011773,
     $     0.011807, 0.011837, 0.011869, 0.011901, 0.01194, 0.011978,
     $     0.012017, 0.012056, 0.012095, 0.012134, 0.012173, 0.012211,
     $     0.012249, 0.012287, 0.012324, 0.01236, 0.012396, 0.012431,
     $     0.012465, 0.012499, 0.012532, 0.01256, 0.012587, 0.012613,
     $     0.012638, 0.012662, 0.012685, 0.012707, 0.012728, 0.012749,
     $     0.012768, 0.012786, 0.012803, 0.012819, 0.012834, 0.012848,
     $     0.012862, 0.012874, 0.012885, 0.012895, 0.012904, 0.012911,
     $     0.012916, 0.01292, 0.012923, 0.012924, 0.012925, 0.012925,
     $     0.012923, 0.012921, 0.012918, 0.012915, 0.01291, 0.012904,
     $     0.012898, 0.012891, 0.012883, 0.012874, 0.012864, 0.012853,
     $     0.012842, 0.012829, 0.012815, 0.0128, 0.012613, 0.012362,
     $     0.012071, 0.011737, 0.011373, 0.010991, 0.010602, 0.010217,
     $     0.0098455, 0.0094959, 0.0091748, 0.008893, 0.0086579,
     $     0.0084768/


C     64160, JENDL evaluated data
      nbins(48) = 228
      data ( eg(48,i), i=1, 228) /7452400.0, 7600000.0, 7800000.0,
     $     8000000.0, 8200000.0, 8400000.0, 8600000.0, 8800000.0,
     $     9000000.0, 9200000.0, 9400000.0, 9600000.0, 9800000.0, 1.0D7,
     $     1.02D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7
     $     , 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7,
     $     1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7
     $     , 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7,
     $     1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7
     $     , 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7,
     $     1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7
     $     , 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7,
     $     2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7
     $     , 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7,
     $     2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7
     $     , 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7,
     $     2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7
     $     , 2.98D7, 3.0D7, 3.05D7, 3.1D7, 3.15D7, 3.2D7, 3.25D7, 3.3D7,
     $     3.35D7, 3.4D7, 3.45D7, 3.5D7, 3.55D7, 3.6D7, 3.65D7, 3.7D7,
     $     3.75D7, 3.8D7, 3.85D7, 3.9D7, 3.95D7, 4.0D7, 4.05D7, 4.1D7,
     $     4.15D7, 4.2D7, 4.25D7, 4.3D7, 4.35D7, 4.4D7, 4.45D7, 4.5D7,
     $     4.55D7, 4.6D7, 4.65D7, 4.7D7, 4.75D7, 4.8D7, 4.85D7, 4.9D7,
     $     4.95D7, 5.0D7, 5.05D7, 5.1D7, 5.15D7, 5.2D7, 5.25D7, 5.3D7,
     $     5.35D7, 5.4D7, 5.45D7, 5.5D7, 5.55D7, 5.6D7, 5.65D7, 5.7D7,
     $     5.75D7, 5.8D7, 5.85D7, 5.9D7, 5.95D7, 6.0D7, 6.05D7, 6.1D7,
     $     6.15D7, 6.2D7, 6.25D7, 6.3D7, 6.35D7, 6.4D7, 6.45D7, 6.5D7,
     $     6.55D7, 6.6D7, 6.65D7, 6.7D7, 6.75D7, 6.8D7, 6.85D7, 6.9D7,
     $     6.95D7, 7.0D7, 7.05D7, 7.1D7, 7.15D7, 7.2D7, 7.25D7, 7.3D7,
     $     7.35D7, 7.4D7, 7.45D7, 7.5D7, 7.55D7, 7.6D7, 7.65D7, 7.7D7,
     $     7.75D7, 7.8D7, 7.85D7, 7.9D7, 7.95D7, 8.0D7, 8.5D7, 9.0D7,
     $     9.5D7, 1.0D8, 1.05D8, 1.1D8, 1.15D8, 1.2D8, 1.25D8, 1.3D8,
     $     1.35D8, 1.4D8, 1.45D8, 1.5D8/

      data ( gdr(48,i), i=1, 228) /0.0, 0.0040748, 0.015712, 0.024186,
     $     0.030456, 0.035216, 0.038666, 0.041743, 0.04498, 0.048716,
     $     0.053366, 0.057994, 0.063635, 0.071174, 0.081082, 0.09385,
     $     0.10994, 0.1297, 0.15312, 0.17966, 0.20792, 0.23391, 0.25605,
     $     0.27246, 0.28173, 0.28405, 0.28062, 0.27323, 0.26394, 0.25446
     $     , 0.24587, 0.24131, 0.24084, 0.24346, 0.24848, 0.25529,
     $     0.26141, 0.26828, 0.2744, 0.27848, 0.28183, 0.27927, 0.27445,
     $     0.26779, 0.26002, 0.25204, 0.24472, 0.23861, 0.23385, 0.23018
     $     , 0.22777, 0.22387, 0.21872, 0.21184, 0.20311, 0.19268, 0.181
     $     , 0.16868, 0.1563, 0.1443, 0.1327, 0.12199, 0.11239, 0.10386,
     $     0.096354, 0.089791, 0.084093, 0.079182, 0.074988, 0.07145,
     $     0.068512, 0.066193, 0.064337, 0.062943, 0.062073, 0.061694,
     $     0.061763, 0.062298, 0.063275, 0.064649, 0.066315, 0.068133,
     $     0.070061, 0.071957, 0.073648, 0.074836, 0.075528, 0.075625,
     $     0.075079, 0.073896, 0.072107, 0.069774, 0.06703, 0.063993,
     $     0.060784, 0.057513, 0.054267, 0.051089, 0.048039, 0.045169,
     $     0.042491, 0.040008, 0.037721, 0.035625, 0.033706, 0.031952,
     $     0.03035, 0.028887, 0.02755, 0.026327, 0.025199, 0.024164,
     $     0.023216, 0.022346, 0.020477, 0.018957, 0.017709, 0.016673,
     $     0.015809, 0.015087, 0.014483, 0.013973, 0.013541, 0.013173,
     $     0.012863, 0.012601, 0.01238, 0.012191, 0.012031, 0.011897,
     $     0.011786, 0.011695, 0.011619, 0.011557, 0.011523, 0.0115,
     $     0.011487, 0.011481, 0.011481, 0.011488, 0.011501, 0.011519,
     $     0.01154, 0.011563, 0.011591, 0.011621, 0.011654, 0.011688,
     $     0.011723, 0.011759, 0.011797, 0.011836, 0.011873, 0.011909,
     $     0.011951, 0.011993, 0.012035, 0.012077, 0.012119, 0.01216,
     $     0.012202, 0.012242, 0.012283, 0.012323, 0.012362, 0.0124,
     $     0.012438, 0.012475, 0.012511, 0.012547, 0.012581, 0.012615,
     $     0.012645, 0.012673, 0.012699, 0.012724, 0.012749, 0.012772,
     $     0.012795, 0.012816, 0.012837, 0.012856, 0.012875, 0.012892,
     $     0.012908, 0.012923, 0.012937, 0.012951, 0.012963, 0.012974,
     $     0.012984, 0.012993, 0.013, 0.013005, 0.013008, 0.01301,
     $     0.013011, 0.013012, 0.013011, 0.013009, 0.013007, 0.013004,
     $     0.012999, 0.012994, 0.012988, 0.012981, 0.012974, 0.012965,
     $     0.012956, 0.012945, 0.012934, 0.012922, 0.012909, 0.012895,
     $     0.012709, 0.012459, 0.012164, 0.011829, 0.011461, 0.011077,
     $     0.010684, 0.010296, 0.009921, 0.0095683, 0.009244, 0.0089593,
     $     0.0087219, 0.0085389/


C     73181, JENDL evaluated data
      nbins(49) = 81
      data ( eg(49,i), i=1, 81) /7800000.0, 8100000.0, 8400000.0,
     $     8600000.0, 8900000.0, 9200000.0, 9400000.0, 9700000.0, 1.0D7,
     $     1.04D7, 1.08D7, 1.12D7, 1.16D7, 1.2D7, 1.24D7, 1.28D7, 1.32D7
     $     , 1.36D7, 1.4D7, 1.44D7, 1.48D7, 1.52D7, 1.56D7, 1.6D7,
     $     1.64D7, 1.68D7, 1.72D7, 1.76D7, 1.8D7, 1.84D7, 1.88D7, 1.92D7
     $     , 1.96D7, 2.0D7, 2.04D7, 2.08D7, 2.12D7, 2.16D7, 2.2D7,
     $     2.24D7, 2.28D7, 2.32D7, 2.36D7, 2.4D7, 2.44D7, 2.48D7, 2.52D7
     $     , 2.56D7, 2.6D7, 2.64D7, 2.68D7, 2.72D7, 2.76D7, 2.8D7,
     $     2.84D7, 2.88D7, 2.92D7, 2.96D7, 3.0D7, 3.04D7, 3.08D7, 3.1D7,
     $     3.4D7, 3.7D7, 4.0D7, 4.3D7, 4.5D7, 5.0D7, 5.5D7, 6.0D7, 6.5D7
     $     , 7.0D7, 8.0D7, 9.0D7, 1.0D8, 1.1D8, 1.2D8, 1.3D8, 1.4D8,
     $     1.5D8, 1.6D8/

      data ( gdr(49,i), i=1, 81) /0.0169, 0.0278, 0.0286, 0.0346, 0.0426
     $     , 0.0553, 0.0602, 0.0689, 0.0837, 0.1085, 0.1443, 0.1965,
     $     0.2675, 0.3398, 0.3718, 0.3605, 0.3456, 0.345, 0.3546, 0.3693
     $     , 0.3712, 0.3603, 0.3374, 0.3068, 0.2732, 0.2404, 0.2105,
     $     0.1842, 0.1618, 0.1427, 0.1266, 0.113, 0.1014, 0.0916, 0.0833
     $     , 0.0761, 0.0699, 0.0645, 0.0598, 0.0556, 0.052, 0.0488,
     $     0.0459, 0.0434, 0.0411, 0.0391, 0.0372, 0.0356, 0.034, 0.0327
     $     , 0.0314, 0.0303, 0.0293, 0.0283, 0.0275, 0.0267, 0.0259,
     $     0.0253, 0.0246, 0.0241, 0.0235, 0.0233, 0.0205, 0.0188,
     $     0.0178, 0.0171, 0.0168, 0.0164, 0.0161, 0.0159, 0.0157,
     $     0.0154, 0.0148, 0.014, 0.0131, 0.0121, 0.0112, 0.0103, 0.0096
     $     , 0.0091, 0.0089/


C     74182, JENDL evaluated data
      nbins(50) = 280
      data ( eg(50,i), i=1, 280) /500000.0, 1000000.0, 1500000.0,
     $     2000000.0, 2500000.0, 3000000.0, 3500000.0, 4000000.0,
     $     4500000.0, 5000000.0, 5500000.0, 6000000.0, 6500000.0,
     $     7000000.0, 7500000.0, 8000000.0, 8500000.0, 9000000.0,
     $     9500000.0, 1.0D7, 1.05D7, 1.1D7, 1.15D7, 1.2D7, 1.25D7, 1.3D7
     $     , 1.35D7, 1.4D7, 1.45D7, 1.5D7, 1.55D7, 1.6D7, 1.65D7, 1.7D7,
     $     1.75D7, 1.8D7, 1.85D7, 1.9D7, 1.95D7, 2.0D7, 2.05D7, 2.1D7,
     $     2.15D7, 2.2D7, 2.25D7, 2.3D7, 2.35D7, 2.4D7, 2.45D7, 2.5D7,
     $     2.55D7, 2.6D7, 2.65D7, 2.7D7, 2.75D7, 2.8D7, 2.85D7, 2.9D7,
     $     2.95D7, 3.0D7, 3.05D7, 3.1D7, 3.15D7, 3.2D7, 3.25D7, 3.3D7,
     $     3.35D7, 3.4D7, 3.45D7, 3.5D7, 3.55D7, 3.6D7, 3.65D7, 3.7D7,
     $     3.75D7, 3.8D7, 3.85D7, 3.9D7, 3.95D7, 4.0D7, 4.05D7, 4.1D7,
     $     4.15D7, 4.2D7, 4.25D7, 4.3D7, 4.35D7, 4.4D7, 4.45D7, 4.5D7,
     $     4.55D7, 4.6D7, 4.65D7, 4.7D7, 4.75D7, 4.8D7, 4.85D7, 4.9D7,
     $     4.95D7, 5.0D7, 5.05D7, 5.1D7, 5.15D7, 5.2D7, 5.25D7, 5.3D7,
     $     5.35D7, 5.4D7, 5.45D7, 5.5D7, 5.55D7, 5.6D7, 5.65D7, 5.7D7,
     $     5.75D7, 5.8D7, 5.85D7, 5.9D7, 5.95D7, 6.0D7, 6.05D7, 6.1D7,
     $     6.15D7, 6.2D7, 6.25D7, 6.3D7, 6.35D7, 6.4D7, 6.45D7, 6.5D7,
     $     6.55D7, 6.6D7, 6.65D7, 6.7D7, 6.75D7, 6.8D7, 6.85D7, 6.9D7,
     $     6.95D7, 7.0D7, 7.05D7, 7.1D7, 7.15D7, 7.2D7, 7.25D7, 7.3D7,
     $     7.35D7, 7.4D7, 7.45D7, 7.5D7, 7.55D7, 7.6D7, 7.65D7, 7.7D7,
     $     7.75D7, 7.8D7, 7.85D7, 7.9D7, 7.95D7, 8.0D7, 8.05D7, 8.1D7,
     $     8.15D7, 8.2D7, 8.25D7, 8.3D7, 8.35D7, 8.4D7, 8.45D7, 8.5D7,
     $     8.55D7, 8.6D7, 8.65D7, 8.7D7, 8.75D7, 8.8D7, 8.85D7, 8.9D7,
     $     8.95D7, 9.0D7, 9.05D7, 9.1D7, 9.15D7, 9.2D7, 9.25D7, 9.3D7,
     $     9.35D7, 9.4D7, 9.45D7, 9.5D7, 9.55D7, 9.6D7, 9.65D7, 9.7D7,
     $     9.75D7, 9.8D7, 9.85D7, 9.9D7, 9.95D7, 1.0D8, 1.005D8, 1.01D8,
     $     1.015D8, 1.02D8, 1.025D8, 1.03D8, 1.035D8, 1.04D8, 1.045D8,
     $     1.05D8, 1.055D8, 1.06D8, 1.065D8, 1.07D8, 1.075D8, 1.08D8,
     $     1.085D8, 1.09D8, 1.095D8, 1.1D8, 1.105D8, 1.11D8, 1.115D8,
     $     1.12D8, 1.125D8, 1.13D8, 1.135D8, 1.14D8, 1.145D8, 1.15D8,
     $     1.155D8, 1.16D8, 1.165D8, 1.17D8, 1.175D8, 1.18D8, 1.185D8,
     $     1.19D8, 1.195D8, 1.2D8, 1.205D8, 1.21D8, 1.215D8, 1.22D8,
     $     1.225D8, 1.23D8, 1.235D8, 1.24D8, 1.245D8, 1.25D8, 1.255D8,
     $     1.26D8, 1.265D8, 1.27D8, 1.275D8, 1.28D8, 1.285D8, 1.29D8,
     $     1.295D8, 1.3D8, 1.305D8, 1.31D8, 1.315D8, 1.32D8, 1.325D8,
     $     1.33D8, 1.335D8, 1.34D8, 1.345D8, 1.35D8, 1.355D8, 1.36D8,
     $     1.365D8, 1.37D8, 1.375D8, 1.38D8, 1.385D8, 1.39D8, 1.395D8,
     $     1.4D8/

      data ( gdr(50,i), i=1, 280) /6.1364D-5, 2.4716D-4, 5.6261D-4,
     $     0.0010167, 0.0016227, 0.0023988, 0.0033695, 0.0045665,
     $     0.0060315, 0.0078188, 0.0099994, 0.012667, 0.015948, 0.020011
     $     , 0.025094, 0.03153, 0.039804, 0.050645, 0.065195, 0.085323,
     $     0.11423, 0.15747, 0.22361, 0.3166, 0.40204, 0.42257, 0.40988,
     $     0.40564, 0.40508, 0.39359, 0.36574, 0.32606, 0.28252, 0.24119
     $     , 0.20499, 0.17458, 0.14953, 0.12902, 0.11222, 0.098394,
     $     0.086936, 0.077367, 0.06931, 0.062473, 0.056627, 0.05159,
     $     0.047222, 0.043409, 0.04006, 0.037103, 0.034479, 0.032138,
     $     0.03004, 0.028153, 0.026448, 0.024903, 0.023497, 0.022214,
     $     0.021039, 0.019961, 0.018969, 0.018053, 0.017206, 0.016421,
     $     0.015691, 0.015012, 0.014379, 0.013787, 0.013233, 0.012714,
     $     0.012227, 0.011768, 0.011337, 0.01093, 0.010545, 0.010182,
     $     0.0098385, 0.0093651, 0.0094797, 0.009592, 0.0097021, 0.00981
     $     , 0.0099156, 0.010019, 0.01012, 0.010219, 0.010316, 0.01041,
     $     0.010502, 0.010592, 0.01068, 0.010766, 0.010849, 0.010931,
     $     0.01101, 0.011087, 0.011162, 0.011235, 0.011306, 0.011375,
     $     0.011442, 0.011507, 0.011569, 0.01163, 0.01169, 0.011747,
     $     0.011802, 0.011856, 0.011907, 0.011957, 0.012005, 0.012051,
     $     0.012096, 0.012138, 0.012179, 0.012219, 0.012256, 0.012292,
     $     0.012327, 0.01236, 0.012391, 0.012421, 0.012449, 0.012476,
     $     0.012501, 0.012525, 0.012547, 0.012568, 0.012587, 0.012605,
     $     0.012622, 0.012637, 0.012651, 0.012664, 0.012676, 0.012686,
     $     0.012695, 0.012703, 0.012709, 0.012714, 0.012719, 0.012722,
     $     0.012724, 0.012724, 0.012724, 0.012723, 0.01272, 0.012717,
     $     0.012713, 0.012707, 0.012701, 0.012693, 0.012685, 0.012676,
     $     0.012666, 0.012654, 0.012643, 0.01263, 0.012616, 0.012602,
     $     0.012586, 0.01257, 0.012553, 0.012536, 0.012517, 0.012498,
     $     0.012478, 0.012458, 0.012436, 0.012414, 0.012392, 0.012369,
     $     0.012345, 0.01232, 0.012295, 0.01227, 0.012243, 0.012216,
     $     0.012189, 0.012161, 0.012133, 0.012104, 0.012075, 0.012045,
     $     0.012014, 0.011984, 0.011952, 0.011921, 0.011889, 0.011856,
     $     0.011823, 0.01179, 0.011756, 0.011723, 0.011688, 0.011654,
     $     0.011619, 0.011584, 0.011548, 0.011512, 0.011476, 0.01144,
     $     0.011403, 0.011367, 0.01133, 0.011292, 0.011255, 0.011217,
     $     0.01118, 0.011142, 0.011104, 0.011065, 0.011027, 0.010988,
     $     0.01095, 0.010911, 0.010872, 0.010834, 0.010795, 0.010756,
     $     0.010717, 0.010677, 0.010638, 0.010599, 0.01056, 0.010521,
     $     0.010482, 0.010443, 0.010404, 0.010365, 0.010326, 0.010287,
     $     0.010248, 0.010209, 0.01017, 0.010132, 0.010093, 0.010055,
     $     0.010017, 0.0099784, 0.0099404, 0.0099026, 0.0098649,
     $     0.0098274, 0.0097901, 0.0097529, 0.009716, 0.0096792,
     $     0.0096426, 0.0096063, 0.0095702, 0.0095343, 0.0094987,
     $     0.0094633, 0.0094283, 0.0093934, 0.0093588, 0.0093246,
     $     0.0092906, 0.0092569, 0.0092235, 0.0091905, 0.0091577,
     $     0.0091253, 0.0090933, 0.0090616, 0.0090302, 0.0089992,
     $     0.0089686, 0.0089384, 0.0089086, 0.0088791, 0.00885,
     $     0.0088214, 0.0087932, 0.0087654, 0.008738, 0.0087111,
     $     0.0086846, 0.0086585/


C     74184, JENDL evaluated data
      nbins(51) = 280
      data ( eg(51,i), i=1, 280) /500000.0, 1000000.0, 1500000.0,
     $     2000000.0, 2500000.0, 3000000.0, 3500000.0, 4000000.0,
     $     4500000.0, 5000000.0, 5500000.0, 6000000.0, 6500000.0,
     $     7000000.0, 7500000.0, 8000000.0, 8500000.0, 9000000.0,
     $     9500000.0, 1.0D7, 1.05D7, 1.1D7, 1.15D7, 1.2D7, 1.25D7, 1.3D7
     $     , 1.35D7, 1.4D7, 1.45D7, 1.5D7, 1.55D7, 1.6D7, 1.65D7, 1.7D7,
     $     1.75D7, 1.8D7, 1.85D7, 1.9D7, 1.95D7, 2.0D7, 2.05D7, 2.1D7,
     $     2.15D7, 2.2D7, 2.25D7, 2.3D7, 2.35D7, 2.4D7, 2.45D7, 2.5D7,
     $     2.55D7, 2.6D7, 2.65D7, 2.7D7, 2.75D7, 2.8D7, 2.85D7, 2.9D7,
     $     2.95D7, 3.0D7, 3.05D7, 3.1D7, 3.15D7, 3.2D7, 3.25D7, 3.3D7,
     $     3.35D7, 3.4D7, 3.45D7, 3.5D7, 3.55D7, 3.6D7, 3.65D7, 3.7D7,
     $     3.75D7, 3.8D7, 3.85D7, 3.9D7, 3.95D7, 4.0D7, 4.05D7, 4.1D7,
     $     4.15D7, 4.2D7, 4.25D7, 4.3D7, 4.35D7, 4.4D7, 4.45D7, 4.5D7,
     $     4.55D7, 4.6D7, 4.65D7, 4.7D7, 4.75D7, 4.8D7, 4.85D7, 4.9D7,
     $     4.95D7, 5.0D7, 5.05D7, 5.1D7, 5.15D7, 5.2D7, 5.25D7, 5.3D7,
     $     5.35D7, 5.4D7, 5.45D7, 5.5D7, 5.55D7, 5.6D7, 5.65D7, 5.7D7,
     $     5.75D7, 5.8D7, 5.85D7, 5.9D7, 5.95D7, 6.0D7, 6.05D7, 6.1D7,
     $     6.15D7, 6.2D7, 6.25D7, 6.3D7, 6.35D7, 6.4D7, 6.45D7, 6.5D7,
     $     6.55D7, 6.6D7, 6.65D7, 6.7D7, 6.75D7, 6.8D7, 6.85D7, 6.9D7,
     $     6.95D7, 7.0D7, 7.05D7, 7.1D7, 7.15D7, 7.2D7, 7.25D7, 7.3D7,
     $     7.35D7, 7.4D7, 7.45D7, 7.5D7, 7.55D7, 7.6D7, 7.65D7, 7.7D7,
     $     7.75D7, 7.8D7, 7.85D7, 7.9D7, 7.95D7, 8.0D7, 8.05D7, 8.1D7,
     $     8.15D7, 8.2D7, 8.25D7, 8.3D7, 8.35D7, 8.4D7, 8.45D7, 8.5D7,
     $     8.55D7, 8.6D7, 8.65D7, 8.7D7, 8.75D7, 8.8D7, 8.85D7, 8.9D7,
     $     8.95D7, 9.0D7, 9.05D7, 9.1D7, 9.15D7, 9.2D7, 9.25D7, 9.3D7,
     $     9.35D7, 9.4D7, 9.45D7, 9.5D7, 9.55D7, 9.6D7, 9.65D7, 9.7D7,
     $     9.75D7, 9.8D7, 9.85D7, 9.9D7, 9.95D7, 1.0D8, 1.005D8, 1.01D8,
     $     1.015D8, 1.02D8, 1.025D8, 1.03D8, 1.035D8, 1.04D8, 1.045D8,
     $     1.05D8, 1.055D8, 1.06D8, 1.065D8, 1.07D8, 1.075D8, 1.08D8,
     $     1.085D8, 1.09D8, 1.095D8, 1.1D8, 1.105D8, 1.11D8, 1.115D8,
     $     1.12D8, 1.125D8, 1.13D8, 1.135D8, 1.14D8, 1.145D8, 1.15D8,
     $     1.155D8, 1.16D8, 1.165D8, 1.17D8, 1.175D8, 1.18D8, 1.185D8,
     $     1.19D8, 1.195D8, 1.2D8, 1.205D8, 1.21D8, 1.215D8, 1.22D8,
     $     1.225D8, 1.23D8, 1.235D8, 1.24D8, 1.245D8, 1.25D8, 1.255D8,
     $     1.26D8, 1.265D8, 1.27D8, 1.275D8, 1.28D8, 1.285D8, 1.29D8,
     $     1.295D8, 1.3D8, 1.305D8, 1.31D8, 1.315D8, 1.32D8, 1.325D8,
     $     1.33D8, 1.335D8, 1.34D8, 1.345D8, 1.35D8, 1.355D8, 1.36D8,
     $     1.365D8, 1.37D8, 1.375D8, 1.38D8, 1.385D8, 1.39D8, 1.395D8,
     $     1.4D8/

      data ( gdr(51,i), i=1, 280) /7.0091D-5, 2.8231D-4, 6.4258D-4,
     $     0.0011611, 0.0018531, 0.0027391, 0.0038469, 0.0052126,
     $     0.0068836, 0.008921, 0.011405, 0.014443, 0.018175, 0.022792,
     $     0.02856, 0.035851, 0.045202, 0.057417, 0.073736, 0.096154,
     $     0.12797, 0.17459, 0.24325, 0.33426, 0.41445, 0.43796, 0.42911
     $     , 0.42317, 0.41903, 0.40538, 0.37761, 0.33925, 0.29698,
     $     0.25622, 0.21988, 0.18883, 0.16286, 0.14135, 0.12354, 0.10875
     $     , 0.096405, 0.086033, 0.077256, 0.069775, 0.063354, 0.057805,
     $     0.052979, 0.048756, 0.04504, 0.041752, 0.038829, 0.036218,
     $     0.033875, 0.031765, 0.029857, 0.028125, 0.026549, 0.025109,
     $     0.02379, 0.022578, 0.021462, 0.020431, 0.019478, 0.018593,
     $     0.017771, 0.017006, 0.016291, 0.015624, 0.014999, 0.014413,
     $     0.013862, 0.013344, 0.012856, 0.012396, 0.011962, 0.011551,
     $     0.011163, 0.010794, 0.0107, 0.0106, 0.0105, 0.0104, 0.0102,
     $     0.010094, 0.010195, 0.010295, 0.010392, 0.010487, 0.01058,
     $     0.010671, 0.01076, 0.010846, 0.01093, 0.011012, 0.011092,
     $     0.011169, 0.011245, 0.011319, 0.01139, 0.011459, 0.011527,
     $     0.011592, 0.011656, 0.011717, 0.011777, 0.011834, 0.01189,
     $     0.011944, 0.011996, 0.012046, 0.012094, 0.012141, 0.012186,
     $     0.012229, 0.01227, 0.01231, 0.012348, 0.012384, 0.012419,
     $     0.012452, 0.012483, 0.012513, 0.012542, 0.012569, 0.012594,
     $     0.012618, 0.01264, 0.012661, 0.012681, 0.012699, 0.012716,
     $     0.012731, 0.012746, 0.012758, 0.01277, 0.01278, 0.012789,
     $     0.012797, 0.012804, 0.012809, 0.012813, 0.012816, 0.012818,
     $     0.012819, 0.012819, 0.012818, 0.012815, 0.012812, 0.012807,
     $     0.012802, 0.012795, 0.012788, 0.012779, 0.01277, 0.01276,
     $     0.012749, 0.012737, 0.012724, 0.01271, 0.012695, 0.01268,
     $     0.012664, 0.012647, 0.012629, 0.01261, 0.012591, 0.012571,
     $     0.01255, 0.012529, 0.012507, 0.012484, 0.012461, 0.012437,
     $     0.012412, 0.012387, 0.012361, 0.012334, 0.012307, 0.01228,
     $     0.012252, 0.012223, 0.012194, 0.012164, 0.012134, 0.012104,
     $     0.012073, 0.012041, 0.012009, 0.011977, 0.011944, 0.011911,
     $     0.011878, 0.011844, 0.01181, 0.011775, 0.01174, 0.011705,
     $     0.01167, 0.011634, 0.011598, 0.011562, 0.011525, 0.011488,
     $     0.011451, 0.011414, 0.011376, 0.011339, 0.011301, 0.011263,
     $     0.011225, 0.011186, 0.011148, 0.011109, 0.01107, 0.011031,
     $     0.010992, 0.010953, 0.010914, 0.010875, 0.010836, 0.010796,
     $     0.010757, 0.010718, 0.010678, 0.010639, 0.010599, 0.01056,
     $     0.01052, 0.010481, 0.010442, 0.010403, 0.010363, 0.010324,
     $     0.010285, 0.010246, 0.010207, 0.010168, 0.01013, 0.010091,
     $     0.010053, 0.010014, 0.0099763, 0.0099384, 0.0099006, 0.009863
     $     , 0.0098255, 0.0097883, 0.0097513, 0.0097145, 0.0096778,
     $     0.0096415, 0.0096053, 0.0095695, 0.0095338, 0.0094985,
     $     0.0094634, 0.0094285, 0.009394, 0.0093598, 0.0093258,
     $     0.0092922, 0.0092589, 0.0092259, 0.0091933, 0.009161,
     $     0.0091291, 0.0090975, 0.0090663, 0.0090354, 0.009005,
     $     0.0089749, 0.0089452, 0.0089159, 0.0088871, 0.0088587,
     $     0.0088306, 0.0088031, 0.0087759, 0.0087493, 0.008723/


C     74186, JENDL evaluated data
      nbins(52) = 104
      data ( eg(52,i), i=1, 104) /500000.0, 1000000.0, 1500000.0,
     $     2000000.0, 2500000.0, 3000000.0, 3500000.0, 4000000.0,
     $     4500000.0, 5000000.0, 5500000.0, 6000000.0, 6500000.0,
     $     7000000.0, 7500000.0, 8000000.0, 8500000.0, 9000000.0,
     $     9500000.0, 1.0D7, 1.05D7, 1.1D7, 1.15D7, 1.2D7, 1.25D7, 1.3D7
     $     , 1.35D7, 1.4D7, 1.45D7, 1.5D7, 1.55D7, 1.6D7, 1.65D7, 1.7D7,
     $     1.75D7, 1.8D7, 1.85D7, 1.9D7, 1.95D7, 2.0D7, 2.05D7, 2.1D7,
     $     2.15D7, 2.2D7, 2.25D7, 2.3D7, 2.35D7, 2.4D7, 2.45D7, 2.5D7,
     $     2.55D7, 2.6D7, 2.65D7, 2.7D7, 2.75D7, 2.8D7, 2.85D7, 2.9D7,
     $     2.95D7, 3.0D7, 3.1D7, 3.2D7, 3.3D7, 3.4D7, 3.5D7, 3.6D7,
     $     3.7D7, 3.8D7, 3.9D7, 4.0D7, 4.2D7, 4.4D7, 4.6D7, 4.8D7, 5.0D7
     $     , 5.2D7, 5.4D7, 5.6D7, 5.8D7, 6.0D7, 6.2D7, 6.4D7, 6.6D7,
     $     6.8D7, 7.0D7, 7.2D7, 7.4D7, 7.6D7, 7.8D7, 8.0D7, 8.2D7, 8.4D7
     $     , 8.6D7, 8.8D7, 9.0D7, 9.2D7, 9.4D7, 9.6D7, 9.8D7, 1.0D8,
     $     1.1D8, 1.2D8, 1.3D8, 1.4D8/

      data ( gdr(52,i), i=1, 104) /6.56D-5, 2.642D-4, 6.013D-4, 0.001086
     $     , 0.001734, 0.002562, 0.003597, 0.004873, 0.006433, 0.008334,
     $     0.01065, 0.01348, 0.01695, 0.02124, 0.02659, 0.03335, 0.04199
     $     , 0.05325, 0.06826, 0.08885, 0.1181, 0.1614, 0.2268, 0.3181,
     $     0.4028, 0.4284, 0.4252, 0.4293, 0.4339, 0.4239, 0.3948,
     $     0.3523, 0.3054, 0.2608, 0.2218, 0.189, 0.1619, 0.1398, 0.1216
     $     , 0.1067, 0.09436, 0.084, 0.07528, 0.06786, 0.06153, 0.05605,
     $     0.05131, 0.04716, 0.04352, 0.04029, 0.03745, 0.03489, 0.03261
     $     , 0.03055, 0.0287, 0.02702, 0.02549, 0.02409, 0.02281,
     $     0.02163, 0.01956, 0.0178, 0.01626, 0.01493, 0.01377, 0.01275,
     $     0.01183, 0.01102, 0.009522, 0.009752, 0.01019, 0.01058,
     $     0.01094, 0.01127, 0.01157, 0.01183, 0.01206, 0.01227, 0.01244
     $     , 0.01258, 0.01271, 0.0128, 0.01288, 0.01294, 0.01298, 0.013,
     $     0.013, 0.01299, 0.01295, 0.01291, 0.01284, 0.01277, 0.01268,
     $     0.01258, 0.01248, 0.01236, 0.01224, 0.01209, 0.01196, 0.0118,
     $     0.01102, 0.01022, 0.009482, 0.008864/


C     79197, JENDL evaluated data
      nbins(53) = 321
      data ( eg(53,i), i=1, 321) /5783030.0, 8150000.0, 8200000.0,
     $     8300000.0, 8400000.0, 8500000.0, 8600000.0, 8700000.0,
     $     8800000.0, 8900000.0, 9000000.0, 9100000.0, 9200000.0,
     $     9300000.0, 9400000.0, 9500000.0, 9600000.0, 9700000.0,
     $     9800000.0, 9900000.0, 1.0D7, 1.01D7, 1.02D7, 1.03D7, 1.04D7,
     $     1.05D7, 1.06D7, 1.07D7, 1.08D7, 1.09D7, 1.1D7, 1.11D7, 1.12D7
     $     , 1.13D7, 1.14D7, 1.15D7, 1.16D7, 1.17D7, 1.18D7, 1.19D7,
     $     1.2D7, 1.21D7, 1.22D7, 1.23D7, 1.24D7, 1.25D7, 1.26D7, 1.27D7
     $     , 1.28D7, 1.29D7, 1.3D7, 1.31D7, 1.32D7, 1.33D7, 1.34D7,
     $     1.35D7, 1.36D7, 1.37D7, 1.38D7, 1.39D7, 1.4D7, 1.41D7, 1.42D7
     $     , 1.43D7, 1.44D7, 1.45D7, 1.46D7, 1.47D7, 1.48D7, 1.49D7,
     $     1.5D7, 1.51D7, 1.52D7, 1.53D7, 1.54D7, 1.55D7, 1.56D7, 1.57D7
     $     , 1.58D7, 1.59D7, 1.6D7, 1.61D7, 1.62D7, 1.63D7, 1.64D7,
     $     1.65D7, 1.66D7, 1.67D7, 1.68D7, 1.69D7, 1.7D7, 1.71D7, 1.72D7
     $     , 1.73D7, 1.74D7, 1.75D7, 1.76D7, 1.77D7, 1.78D7, 1.79D7,
     $     1.8D7, 1.81D7, 1.82D7, 1.83D7, 1.84D7, 1.85D7, 1.86D7, 1.87D7
     $     , 1.88D7, 1.89D7, 1.9D7, 1.91D7, 1.92D7, 1.93D7, 1.94D7,
     $     1.95D7, 1.96D7, 1.97D7, 1.98D7, 1.99D7, 2.0D7, 2.01D7, 2.02D7
     $     , 2.03D7, 2.04D7, 2.05D7, 2.06D7, 2.07D7, 2.08D7, 2.09D7,
     $     2.1D7, 2.11D7, 2.12D7, 2.13D7, 2.14D7, 2.15D7, 2.16D7, 2.17D7
     $     , 2.18D7, 2.19D7, 2.2D7, 2.21D7, 2.22D7, 2.23D7, 2.24D7,
     $     2.25D7, 2.26D7, 2.27D7, 2.28D7, 2.29D7, 2.3D7, 2.31D7, 2.32D7
     $     , 2.33D7, 2.34D7, 2.35D7, 2.36D7, 2.37D7, 2.38D7, 2.39D7,
     $     2.4D7, 2.41D7, 2.42D7, 2.43D7, 2.44D7, 2.45D7, 2.46D7, 2.47D7
     $     , 2.48D7, 2.49D7, 2.5D7, 2.51D7, 2.52D7, 2.53D7, 2.54D7,
     $     2.55D7, 2.56D7, 2.57D7, 2.58D7, 2.59D7, 2.6D7, 2.61D7, 2.62D7
     $     , 2.63D7, 2.64D7, 2.65D7, 2.66D7, 2.67D7, 2.68D7, 2.69D7,
     $     2.7D7, 2.71D7, 2.72D7, 2.73D7, 2.74D7, 2.75D7, 2.76D7, 2.77D7
     $     , 2.78D7, 2.79D7, 2.8D7, 2.81D7, 2.82D7, 2.83D7, 2.84D7,
     $     2.85D7, 2.86D7, 2.87D7, 2.88D7, 2.89D7, 2.9D7, 2.91D7, 2.92D7
     $     , 2.93D7, 2.94D7, 2.95D7, 2.96D7, 2.97D7, 2.98D7, 2.99D7,
     $     3.0D7, 3.02D7, 3.04D7, 3.06D7, 3.08D7, 3.1D7, 3.12D7, 3.14D7,
     $     3.16D7, 3.18D7, 3.2D7, 3.22D7, 3.24D7, 3.26D7, 3.28D7, 3.3D7,
     $     3.32D7, 3.34D7, 3.36D7, 3.38D7, 3.4D7, 3.42D7, 3.44D7, 3.46D7
     $     , 3.48D7, 3.5D7, 3.52D7, 3.54D7, 3.56D7, 3.58D7, 3.6D7,
     $     3.62D7, 3.64D7, 3.66D7, 3.68D7, 3.7D7, 3.72D7, 3.74D7, 3.76D7
     $     , 3.78D7, 3.8D7, 3.82D7, 3.84D7, 3.86D7, 3.88D7, 3.9D7,
     $     3.92D7, 3.94D7, 3.96D7, 3.98D7, 4.0D7, 4.2D7, 4.4D7, 4.6D7,
     $     4.8D7, 5.0D7, 5.2D7, 5.4D7, 5.6D7, 5.8D7, 6.0D7, 6.2D7, 6.4D7
     $     , 6.6D7, 6.8D7, 7.0D7, 7.2D7, 7.4D7, 7.6D7, 7.8D7, 8.0D7,
     $     8.2D7, 8.4D7, 8.6D7, 8.8D7, 9.0D7, 9.2D7, 9.4D7, 9.6D7, 9.8D7
     $     , 1.0D8, 1.02D8, 1.04D8, 1.06D8, 1.08D8, 1.1D8, 1.12D8,
     $     1.14D8, 1.16D8, 1.18D8, 1.2D8, 1.22D8, 1.24D8, 1.26D8, 1.28D8
     $     , 1.3D8, 1.32D8, 1.34D8, 1.36D8, 1.38D8, 1.4D8/

      data ( gdr(53,i), i=1, 321) /0.0, 0.005002, 0.008048, 0.01353,
     $     0.01845, 0.02288, 0.02687, 0.03048, 0.03377, 0.0368, 0.03963,
     $     0.04231, 0.04491, 0.04748, 0.05009, 0.05278, 0.05562, 0.05868
     $     , 0.06199, 0.06564, 0.06967, 0.07414, 0.07911, 0.08465,
     $     0.0908, 0.09763, 0.1052, 0.1136, 0.1228, 0.1329, 0.144,
     $     0.1561, 0.1692, 0.1833, 0.1982, 0.2139, 0.2304, 0.2475,
     $     0.2654, 0.2838, 0.3027, 0.3221, 0.3419, 0.362, 0.3825, 0.4031
     $     , 0.424, 0.445, 0.4656, 0.4835, 0.4984, 0.5103, 0.5195,
     $     0.5261, 0.5304, 0.5324, 0.5325, 0.5306, 0.527, 0.522, 0.5156,
     $     0.508, 0.4994, 0.4899, 0.4798, 0.4692, 0.4583, 0.4473, 0.4362
     $     , 0.425, 0.4139, 0.4027, 0.3915, 0.3804, 0.3693, 0.3583,
     $     0.3473, 0.3364, 0.3256, 0.3149, 0.3044, 0.2939, 0.2837,
     $     0.2736, 0.2636, 0.2539, 0.2444, 0.2351, 0.226, 0.2173, 0.2087
     $     , 0.2006, 0.1926, 0.185, 0.1777, 0.1707, 0.1642, 0.158,
     $     0.1521, 0.1467, 0.1418, 0.1372, 0.1331, 0.1296, 0.1264,
     $     0.1237, 0.1217, 0.1201, 0.1191, 0.1184, 0.118, 0.1177, 0.117,
     $     0.116, 0.1143, 0.1117, 0.1081, 0.1031, 0.09669, 0.08848,
     $     0.07831, 0.07415, 0.07265, 0.07119, 0.06978, 0.06841, 0.06708
     $     , 0.0658, 0.06455, 0.06333, 0.06216, 0.06102, 0.05993,
     $     0.05887, 0.05783, 0.05684, 0.05587, 0.05493, 0.05403, 0.05315
     $     , 0.0523, 0.05148, 0.05069, 0.04992, 0.04917, 0.04844,
     $     0.04775, 0.04707, 0.04641, 0.04577, 0.04515, 0.04453, 0.04394
     $     , 0.04336, 0.0428, 0.04225, 0.04171, 0.04119, 0.04068,
     $     0.04018, 0.0397, 0.03922, 0.03875, 0.0383, 0.03786, 0.03743,
     $     0.03702, 0.03661, 0.03621, 0.03583, 0.035, 0.03464, 0.03429,
     $     0.03394, 0.0336, 0.03327, 0.03295, 0.03262, 0.03231, 0.03201,
     $     0.03172, 0.03143, 0.03115, 0.03088, 0.03061, 0.03035, 0.03009
     $     , 0.02984, 0.02959, 0.02935, 0.02911, 0.02888, 0.02865,
     $     0.02843, 0.02821, 0.028, 0.0278, 0.02759, 0.02739, 0.0272,
     $     0.02701, 0.02682, 0.02664, 0.02646, 0.02628, 0.02611, 0.02594
     $     , 0.02577, 0.0256, 0.02544, 0.02529, 0.02513, 0.02499,
     $     0.02484, 0.0247, 0.02456, 0.02442, 0.02429, 0.02416, 0.02402,
     $     0.02389, 0.0237, 0.02346, 0.02323, 0.02301, 0.0228, 0.02259,
     $     0.02239, 0.02221, 0.02202, 0.02184, 0.02168, 0.0215, 0.02135,
     $     0.0212, 0.02105, 0.02091, 0.02077, 0.02064, 0.02051, 0.02039,
     $     0.02027, 0.02016, 0.02006, 0.01995, 0.01984, 0.01974, 0.01965
     $     , 0.01956, 0.01947, 0.01939, 0.0193, 0.01922, 0.01915,
     $     0.01907, 0.019, 0.01893, 0.01887, 0.0188, 0.01874, 0.01867,
     $     0.01862, 0.01856, 0.01851, 0.01846, 0.0184, 0.01836, 0.0183,
     $     0.01826, 0.01822, 0.01817, 0.01782, 0.01757, 0.01739, 0.01726
     $     , 0.01715, 0.01706, 0.01698, 0.0169, 0.01683, 0.01676, 0.0167
     $     , 0.01663, 0.01655, 0.01647, 0.01638, 0.01627, 0.01616,
     $     0.01604, 0.01591, 0.01576, 0.0156, 0.01543, 0.01525, 0.01508,
     $     0.01489, 0.01469, 0.0145, 0.0143, 0.01409, 0.01389, 0.01367,
     $     0.01347, 0.01325, 0.01305, 0.01284, 0.01263, 0.01242, 0.01222
     $     , 0.01202, 0.01181, 0.01162, 0.01143, 0.01124, 0.01107,
     $     0.01089, 0.01072, 0.01056, 0.01042, 0.01027, 0.01013/


C     80196, JENDL evaluated data
      nbins(54) = 222
      data ( eg(54,i), i=1, 222) /8852400.0, 9000000.0, 9200000.0,
     $     9400000.0, 9600000.0, 9800000.0, 1.0D7, 1.02D7, 1.04D7,
     $     1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7, 1.18D7, 1.2D7,
     $     1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7, 1.32D7, 1.34D7, 1.36D7
     $     , 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7,
     $     1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7, 1.62D7, 1.64D7, 1.66D7
     $     , 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7,
     $     1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7
     $     , 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7, 2.08D7, 2.1D7,
     $     2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7, 2.22D7, 2.24D7, 2.26D7
     $     , 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7, 2.38D7, 2.4D7,
     $     2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7, 2.52D7, 2.54D7, 2.56D7
     $     , 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7, 2.68D7, 2.7D7,
     $     2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7, 2.82D7, 2.84D7, 2.86D7
     $     , 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7, 2.98D7, 3.0D7,
     $     3.05D7, 3.1D7, 3.15D7, 3.2D7, 3.25D7, 3.3D7, 3.35D7, 3.4D7,
     $     3.45D7, 3.5D7, 3.55D7, 3.6D7, 3.65D7, 3.7D7, 3.75D7, 3.8D7,
     $     3.85D7, 3.9D7, 3.95D7, 4.0D7, 4.05D7, 4.1D7, 4.15D7, 4.2D7,
     $     4.25D7, 4.3D7, 4.35D7, 4.4D7, 4.45D7, 4.5D7, 4.55D7, 4.6D7,
     $     4.65D7, 4.7D7, 4.75D7, 4.8D7, 4.85D7, 4.9D7, 4.95D7, 5.0D7,
     $     5.05D7, 5.1D7, 5.15D7, 5.2D7, 5.25D7, 5.3D7, 5.35D7, 5.4D7,
     $     5.45D7, 5.5D7, 5.55D7, 5.6D7, 5.65D7, 5.7D7, 5.75D7, 5.8D7,
     $     5.85D7, 5.9D7, 5.95D7, 6.0D7, 6.05D7, 6.1D7, 6.15D7, 6.2D7,
     $     6.25D7, 6.3D7, 6.35D7, 6.4D7, 6.45D7, 6.5D7, 6.55D7, 6.6D7,
     $     6.65D7, 6.7D7, 6.75D7, 6.772D7, 6.8D7, 6.85D7, 6.9D7, 6.95D7,
     $     7.0D7, 7.05D7, 7.1D7, 7.15D7, 7.2D7, 7.25D7, 7.3D7, 7.35D7,
     $     7.4D7, 7.45D7, 7.5D7, 7.55D7, 7.6D7, 7.65D7, 7.7D7, 7.75D7,
     $     7.8D7, 7.85D7, 7.9D7, 7.95D7, 8.0D7, 8.5D7, 9.0D7, 9.5D7,
     $     1.0D8, 1.05D8, 1.1D8, 1.15D8, 1.2D8, 1.25D8, 1.3D8, 1.35D8,
     $     1.4D8, 1.45D8, 1.5D8/

      data ( gdr(54,i), i=1, 222) /0.0, 0.0038812, 0.020273, 0.039021,
     $     0.058563, 0.077175, 0.092861, 0.10716, 0.12132, 0.13625,
     $     0.14873, 0.161, 0.17607, 0.19595, 0.22198, 0.2554, 0.29703,
     $     0.34677, 0.40295, 0.46201, 0.51898, 0.56484, 0.59869, 0.61609
     $     , 0.62068, 0.61713, 0.60611, 0.58723, 0.5601, 0.52295,
     $     0.47943, 0.43366, 0.38858, 0.34659, 0.30907, 0.27682, 0.26878
     $     , 0.26051, 0.25154, 0.24247, 0.23367, 0.22469, 0.21281,
     $     0.20134, 0.19071, 0.18082, 0.17156, 0.16303, 0.15526, 0.14831
     $     , 0.14222, 0.13701, 0.13238, 0.12861, 0.12577, 0.12374,
     $     0.12229, 0.12133, 0.12059, 0.11974, 0.11846, 0.11647, 0.11363
     $     , 0.10973, 0.10488, 0.099327, 0.093453, 0.087485, 0.081623,
     $     0.076016, 0.070759, 0.065902, 0.061459, 0.057425, 0.053775,
     $     0.05048, 0.047509, 0.044798, 0.042342, 0.040122, 0.038112,
     $     0.036288, 0.034592, 0.033044, 0.031627, 0.030328, 0.029135,
     $     0.028037, 0.027024, 0.026085, 0.025207, 0.024393, 0.023637,
     $     0.022933, 0.022279, 0.021668, 0.021099, 0.020566, 0.020068,
     $     0.019602, 0.019165, 0.018756, 0.018366, 0.017998, 0.017651,
     $     0.017325, 0.017018, 0.01637, 0.015819, 0.015351, 0.014942,
     $     0.014588, 0.014291, 0.01404, 0.013829, 0.013644, 0.013489,
     $     0.013394, 0.013322, 0.01327, 0.01323, 0.013204, 0.013192,
     $     0.013193, 0.013205, 0.013221, 0.013245, 0.013283, 0.013328,
     $     0.013379, 0.013428, 0.013479, 0.013534, 0.013593, 0.013655,
     $     0.013715, 0.013777, 0.013841, 0.013907, 0.013973, 0.014035,
     $     0.014097, 0.014159, 0.014221, 0.014283, 0.014346, 0.014408,
     $     0.014474, 0.014539, 0.014604, 0.014667, 0.014731, 0.014793,
     $     0.014854, 0.014914, 0.014974, 0.015032, 0.015089, 0.015144,
     $     0.015199, 0.015246, 0.015292, 0.015336, 0.015379, 0.015421,
     $     0.015461, 0.0155, 0.015538, 0.015573, 0.015608, 0.015641,
     $     0.015672, 0.015702, 0.015731, 0.015759, 0.015785, 0.015809,
     $     0.015832, 0.015854, 0.015875, 0.015891, 0.015906, 0.01591172,
     $     0.015919, 0.015931, 0.015942, 0.015951, 0.01596, 0.015966,
     $     0.015971, 0.015975, 0.015977, 0.015978, 0.015979, 0.015978,
     $     0.015975, 0.015972, 0.015968, 0.015962, 0.015955, 0.015948,
     $     0.015937, 0.015926, 0.015913, 0.015899, 0.015885, 0.015869,
     $     0.015852, 0.015633, 0.015329, 0.01497, 0.014557, 0.014107,
     $     0.013633, 0.013151, 0.012674, 0.012213, 0.011778, 0.01138,
     $     0.011031, 0.01074, 0.010516/


C     80198, JENDL evaluated data
      nbins(55) = 224
      data ( eg(55,i), i=1, 224) /8485400.0, 8600000.0, 8800000.0,
     $     9000000.0, 9200000.0, 9400000.0, 9600000.0, 9800000.0, 1.0D7,
     $     1.02D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7
     $     , 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7,
     $     1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7
     $     , 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7,
     $     1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7
     $     , 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7,
     $     1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7
     $     , 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7,
     $     2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7
     $     , 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7,
     $     2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7
     $     , 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7,
     $     2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7
     $     , 2.98D7, 3.0D7, 3.05D7, 3.1D7, 3.15D7, 3.2D7, 3.25D7, 3.3D7,
     $     3.35D7, 3.4D7, 3.45D7, 3.5D7, 3.55D7, 3.6D7, 3.65D7, 3.7D7,
     $     3.75D7, 3.8D7, 3.85D7, 3.9D7, 3.95D7, 4.0D7, 4.05D7, 4.1D7,
     $     4.15D7, 4.2D7, 4.25D7, 4.3D7, 4.35D7, 4.4D7, 4.45D7, 4.5D7,
     $     4.55D7, 4.6D7, 4.65D7, 4.7D7, 4.75D7, 4.8D7, 4.85D7, 4.9D7,
     $     4.95D7, 5.0D7, 5.05D7, 5.1D7, 5.15D7, 5.2D7, 5.25D7, 5.3D7,
     $     5.35D7, 5.4D7, 5.45D7, 5.5D7, 5.55D7, 5.6D7, 5.65D7, 5.7D7,
     $     5.75D7, 5.8D7, 5.85D7, 5.9D7, 5.95D7, 6.0D7, 6.05D7, 6.1D7,
     $     6.15D7, 6.2D7, 6.25D7, 6.3D7, 6.35D7, 6.4D7, 6.45D7, 6.5D7,
     $     6.55D7, 6.6D7, 6.65D7, 6.7D7, 6.75D7, 6.8D7, 6.85D7, 6.9D7,
     $     6.95D7, 7.0D7, 7.05D7, 7.1D7, 7.15D7, 7.2D7, 7.25D7, 7.3D7,
     $     7.35D7, 7.4D7, 7.45D7, 7.5D7, 7.55D7, 7.6D7, 7.65D7, 7.7D7,
     $     7.75D7, 7.8D7, 7.85D7, 7.9D7, 7.95D7, 8.0D7, 8.5D7, 8.995D7,
     $     9.0D7, 9.5D7, 1.0D8, 1.05D8, 1.1D8, 1.15D8, 1.2D8, 1.25D8,
     $     1.3D8, 1.35D8, 1.4D8, 1.45D8, 1.5D8/

      data ( gdr(55,i), i=1, 224) /0.0, 0.0010955, 0.010974, 0.024152,
     $     0.039957, 0.057084, 0.073384, 0.088665, 0.10359, 0.11824,
     $     0.13319, 0.14701, 0.1573, 0.17025, 0.18709, 0.20904, 0.23737,
     $     0.2732, 0.31709, 0.36856, 0.42545, 0.48106, 0.5309, 0.57229,
     $     0.60283, 0.61568, 0.61797, 0.61231, 0.59914, 0.578, 0.5487,
     $     0.50814, 0.46185, 0.4139, 0.36819, 0.3302, 0.3099, 0.29307,
     $     0.27885, 0.26621, 0.25412, 0.24295, 0.23293, 0.2231, 0.21059,
     $     0.19912, 0.1884, 0.17842, 0.16917, 0.16068, 0.15298, 0.14612,
     $     0.14015, 0.13509, 0.13097, 0.12768, 0.12519, 0.12313, 0.12171
     $     , 0.1207, 0.11981, 0.11873, 0.11713, 0.11479, 0.1116, 0.10726
     $     , 0.10222, 0.096714, 0.090952, 0.085089, 0.079319, 0.073845,
     $     0.068744, 0.064052, 0.059773, 0.055892, 0.052385, 0.049221,
     $     0.046367, 0.043791, 0.041463, 0.039357, 0.037417, 0.035656,
     $     0.033995, 0.032477, 0.031084, 0.029806, 0.028628, 0.027543,
     $     0.02654, 0.025611, 0.024751, 0.023952, 0.023204, 0.0225,
     $     0.021844, 0.021231, 0.020658, 0.020121, 0.019618, 0.019146,
     $     0.018703, 0.018287, 0.017896, 0.017528, 0.017182, 0.016847,
     $     0.016531, 0.015944, 0.015451, 0.015022, 0.014663, 0.014352,
     $     0.014089, 0.013871, 0.013692, 0.013547, 0.013424, 0.013336,
     $     0.013272, 0.013227, 0.0132, 0.013182, 0.013177, 0.013183,
     $     0.0132, 0.013227, 0.013257, 0.013296, 0.01334, 0.013389,
     $     0.013444, 0.013499, 0.013557, 0.013618, 0.013682, 0.013749,
     $     0.013814, 0.01388, 0.013947, 0.014016, 0.014085, 0.01415,
     $     0.014214, 0.014278, 0.014342, 0.014405, 0.014469, 0.014536,
     $     0.014602, 0.014668, 0.014733, 0.014797, 0.01486, 0.014923,
     $     0.014984, 0.015044, 0.015103, 0.015161, 0.015217, 0.015273,
     $     0.015327, 0.015374, 0.015419, 0.015463, 0.015505, 0.015546,
     $     0.015586, 0.015624, 0.01566, 0.015695, 0.015729, 0.015761,
     $     0.015792, 0.015821, 0.015849, 0.015875, 0.015901, 0.015924,
     $     0.015947, 0.015968, 0.015987, 0.016003, 0.016017, 0.016029,
     $     0.01604, 0.01605, 0.016058, 0.016064, 0.016069, 0.016073,
     $     0.016075, 0.016076, 0.016076, 0.016075, 0.016073, 0.016069,
     $     0.016064, 0.016059, 0.016052, 0.016044, 0.016035, 0.016023,
     $     0.016011, 0.015997, 0.015982, 0.015966, 0.015949, 0.015729,
     $     0.01535107, 0.015425, 0.015064, 0.01465, 0.014198, 0.013721,
     $     0.013238, 0.01276, 0.012297, 0.011859, 0.011458, 0.011107,
     $     0.010813, 0.010588/


C     80199, JENDL evaluated data
      nbins(56) = 233
      data ( eg(56,i), i=1, 233) /6663400.0, 6800000.0, 7000000.0,
     $     7200000.0, 7400000.0, 7600000.0, 7800000.0, 8000000.0,
     $     8200000.0, 8400000.0, 8600000.0, 8800000.0, 9000000.0,
     $     9200000.0, 9400000.0, 9600000.0, 9800000.0, 1.0D7, 1.02D7,
     $     1.04D7, 1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7, 1.18D7
     $     , 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7, 1.32D7,
     $     1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7, 1.48D7
     $     , 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7, 1.62D7,
     $     1.64D7, 1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7, 1.78D7
     $     , 1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7, 1.92D7,
     $     1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7, 2.08D7
     $     , 2.1D7, 2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7, 2.22D7,
     $     2.24D7, 2.26D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7, 2.38D7
     $     , 2.4D7, 2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7, 2.52D7,
     $     2.54D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7, 2.68D7
     $     , 2.7D7, 2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7, 2.82D7,
     $     2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7, 2.98D7
     $     , 3.0D7, 3.05D7, 3.1D7, 3.15D7, 3.2D7, 3.25D7, 3.3D7, 3.35D7,
     $     3.4D7, 3.45D7, 3.5D7, 3.55D7, 3.6D7, 3.65D7, 3.7D7, 3.75D7,
     $     3.8D7, 3.85D7, 3.9D7, 3.95D7, 4.0D7, 4.05D7, 4.1D7, 4.15D7,
     $     4.2D7, 4.25D7, 4.3D7, 4.35D7, 4.4D7, 4.45D7, 4.5D7, 4.55D7,
     $     4.6D7, 4.65D7, 4.7D7, 4.75D7, 4.8D7, 4.85D7, 4.9D7, 4.95D7,
     $     5.0D7, 5.05D7, 5.1D7, 5.15D7, 5.2D7, 5.25D7, 5.3D7, 5.35D7,
     $     5.4D7, 5.45D7, 5.5D7, 5.55D7, 5.6D7, 5.65D7, 5.7D7, 5.75D7,
     $     5.8D7, 5.85D7, 5.9D7, 5.95D7, 6.0D7, 6.05D7, 6.1D7, 6.15D7,
     $     6.2D7, 6.25D7, 6.3D7, 6.35D7, 6.4D7, 6.45D7, 6.5D7, 6.55D7,
     $     6.6D7, 6.65D7, 6.7D7, 6.75D7, 6.8D7, 6.85D7, 6.9D7, 6.95D7,
     $     7.0D7, 7.05D7, 7.1D7, 7.15D7, 7.2D7, 7.25D7, 7.3D7, 7.35D7,
     $     7.4D7, 7.45D7, 7.5D7, 7.55D7, 7.6D7, 7.65D7, 7.7D7, 7.75D7,
     $     7.8D7, 7.85D7, 7.9D7, 7.95D7, 8.0D7, 8.5D7, 9.0D7, 9.5D7,
     $     1.0D8, 1.035D8, 1.05D8, 1.1D8, 1.15D8, 1.2D8, 1.25D8, 1.3D8,
     $     1.35D8, 1.4D8, 1.45D8, 1.5D8/

      data ( gdr(56,i), i=1, 233) /0.0, 3.0357D-4, 0.0020935, 0.0045709,
     $     0.0078746, 0.012072, 0.01686, 0.022613, 0.02978, 0.038519,
     $     0.048915, 0.05937, 0.070231, 0.081945, 0.094216, 0.10673,
     $     0.11923, 0.13165, 0.14414, 0.15713, 0.1711, 0.18349, 0.19791,
     $     0.21577, 0.23809, 0.26581, 0.29968, 0.33996, 0.386, 0.43579,
     $     0.48579, 0.5301, 0.56533, 0.58978, 0.59848, 0.59792, 0.58966,
     $     0.57422, 0.55131, 0.52089, 0.48061, 0.43592, 0.39001, 0.34574
     $     , 0.31782, 0.29707, 0.28036, 0.26618, 0.25345, 0.24145,
     $     0.23087, 0.22137, 0.21143, 0.19942, 0.18844, 0.17828, 0.16883
     $     , 0.16008, 0.15206, 0.14481, 0.13837, 0.13277, 0.12806,
     $     0.12419, 0.12106, 0.11874, 0.11705, 0.11567, 0.11465, 0.11371
     $     , 0.11254, 0.11086, 0.10848, 0.10526, 0.10101, 0.09617,
     $     0.090915, 0.085442, 0.079932, 0.074496, 0.06932, 0.064507,
     $     0.060083, 0.05605, 0.052392, 0.049084, 0.046096, 0.043398,
     $     0.040959, 0.038753, 0.036755, 0.034932, 0.033248, 0.031723,
     $     0.030329, 0.029054, 0.027884, 0.026809, 0.02582, 0.024908,
     $     0.024065, 0.023286, 0.022564, 0.021895, 0.021264, 0.020672,
     $     0.02012, 0.019606, 0.019126, 0.018678, 0.018258, 0.017866,
     $     0.017499, 0.017155, 0.016833, 0.016531, 0.016246, 0.01597,
     $     0.015432, 0.014978, 0.014596, 0.014279, 0.01401, 0.013781,
     $     0.013595, 0.013445, 0.013327, 0.013232, 0.01316, 0.013109,
     $     0.013078, 0.013063, 0.013059, 0.013061, 0.013077, 0.013104,
     $     0.013139, 0.01318, 0.013224, 0.013274, 0.013329, 0.013389,
     $     0.01345, 0.013513, 0.013578, 0.013646, 0.013716, 0.013786,
     $     0.013856, 0.013926, 0.013998, 0.01407, 0.014139, 0.014206,
     $     0.014272, 0.014338, 0.014404, 0.01447, 0.014538, 0.014606,
     $     0.014673, 0.01474, 0.014805, 0.01487, 0.014933, 0.014995,
     $     0.015057, 0.015117, 0.015175, 0.015233, 0.015289, 0.015344,
     $     0.015395, 0.01544, 0.015485, 0.015528, 0.015569, 0.015609,
     $     0.015648, 0.015685, 0.015721, 0.015756, 0.015789, 0.01582,
     $     0.01585, 0.015879, 0.015906, 0.015932, 0.015957, 0.01598,
     $     0.016001, 0.016022, 0.016039, 0.016053, 0.016066, 0.016078,
     $     0.016088, 0.016097, 0.016104, 0.01611, 0.016114, 0.016118,
     $     0.01612, 0.016121, 0.016121, 0.016119, 0.016117, 0.016113,
     $     0.016108, 0.016102, 0.016095, 0.016086, 0.016076, 0.016064,
     $     0.016051, 0.016036, 0.016021, 0.016005, 0.015786, 0.015481,
     $     0.01512, 0.014705, 0.01438185, 0.01425, 0.013771, 0.013283,
     $     0.012801, 0.012336, 0.011897, 0.011495, 0.011142, 0.010848,
     $     0.010621/


C     80200, JENDL evaluated data
      nbins(57) = 226
      data ( eg(57,i), i=1, 226) /8029400.0, 8200000.0, 8400000.0,
     $     8600000.0, 8800000.0, 9000000.0, 9200000.0, 9400000.0,
     $     9600000.0, 9800000.0, 1.0D7, 1.02D7, 1.04D7, 1.06D7, 1.08D7,
     $     1.1D7, 1.12D7, 1.14D7, 1.16D7, 1.18D7, 1.2D7, 1.22D7, 1.24D7,
     $     1.26D7, 1.28D7, 1.3D7, 1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7,
     $     1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7
     $     , 1.58D7, 1.6D7, 1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7,
     $     1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7
     $     , 1.88D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7,
     $     2.02D7, 2.04D7, 2.06D7, 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7
     $     , 2.18D7, 2.2D7, 2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7,
     $     2.32D7, 2.34D7, 2.36D7, 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7
     $     , 2.48D7, 2.5D7, 2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7,
     $     2.62D7, 2.64D7, 2.66D7, 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7
     $     , 2.78D7, 2.8D7, 2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7,
     $     2.92D7, 2.94D7, 2.96D7, 2.98D7, 3.0D7, 3.05D7, 3.1D7, 3.15D7,
     $     3.2D7, 3.25D7, 3.3D7, 3.35D7, 3.4D7, 3.45D7, 3.5D7, 3.55D7,
     $     3.6D7, 3.65D7, 3.7D7, 3.75D7, 3.8D7, 3.85D7, 3.9D7, 3.95D7,
     $     4.0D7, 4.05D7, 4.1D7, 4.15D7, 4.2D7, 4.25D7, 4.3D7, 4.35D7,
     $     4.4D7, 4.45D7, 4.5D7, 4.55D7, 4.6D7, 4.65D7, 4.7D7, 4.75D7,
     $     4.8D7, 4.85D7, 4.9D7, 4.95D7, 5.0D7, 5.05D7, 5.1D7, 5.15D7,
     $     5.2D7, 5.25D7, 5.3D7, 5.35D7, 5.4D7, 5.45D7, 5.5D7, 5.55D7,
     $     5.6D7, 5.65D7, 5.7D7, 5.75D7, 5.8D7, 5.85D7, 5.9D7, 5.95D7,
     $     6.0D7, 6.05D7, 6.1D7, 6.15D7, 6.2D7, 6.25D7, 6.3D7, 6.35D7,
     $     6.4D7, 6.45D7, 6.5D7, 6.55D7, 6.6D7, 6.65D7, 6.7D7, 6.75D7,
     $     6.8D7, 6.85D7, 6.9D7, 6.95D7, 7.0D7, 7.05D7, 7.1D7, 7.15D7,
     $     7.2D7, 7.25D7, 7.3D7, 7.35D7, 7.4D7, 7.45D7, 7.5D7, 7.55D7,
     $     7.6D7, 7.65D7, 7.7D7, 7.75D7, 7.8D7, 7.85D7, 7.9D7, 7.95D7,
     $     8.0D7, 8.5D7, 9.0D7, 9.5D7, 1.0D8, 1.05D8, 1.1D8, 1.15D8,
     $     1.1885D8, 1.2D8, 1.25D8, 1.3D8, 1.35D8, 1.4D8, 1.45D8, 1.5D8/

      data ( gdr(57,i), i=1, 226) /0.0, 0.0020645, 0.0086704, 0.017682,
     $     0.029085, 0.042462, 0.055933, 0.070308, 0.085235, 0.10018,
     $     0.11496, 0.1269, 0.13857, 0.15066, 0.16191, 0.17612, 0.19437,
     $     0.21784, 0.24767, 0.28478, 0.32945, 0.37719, 0.42821, 0.47926
     $     , 0.52555, 0.56294, 0.5889, 0.59733, 0.59691, 0.58908,
     $     0.57411, 0.55064, 0.51932, 0.47761, 0.43408, 0.40646, 0.38292
     $     , 0.36292, 0.34483, 0.32704, 0.30934, 0.29158, 0.27336,
     $     0.25935, 0.24946, 0.23757, 0.22421, 0.2116, 0.19968, 0.18855,
     $     0.17823, 0.16871, 0.16004, 0.15228, 0.14542, 0.13947, 0.1345,
     $     0.1305, 0.12741, 0.12509, 0.12339, 0.12205, 0.12054, 0.1188,
     $     0.11649, 0.11346, 0.10957, 0.10477, 0.099431, 0.09374,
     $     0.087917, 0.082158, 0.076615, 0.071391, 0.066505, 0.061989,
     $     0.057876, 0.054115, 0.050715, 0.047645, 0.044872, 0.042365,
     $     0.040096, 0.038039, 0.03617, 0.03447, 0.032898, 0.03144,
     $     0.030102, 0.028873, 0.027741, 0.026697, 0.025732, 0.024839,
     $     0.024011, 0.023243, 0.022529, 0.021864, 0.021245, 0.020662,
     $     0.020112, 0.019598, 0.019116, 0.018664, 0.01824, 0.017842,
     $     0.017468, 0.017116, 0.016786, 0.016475, 0.016183, 0.01562,
     $     0.015141, 0.01474, 0.014405, 0.014128, 0.013892, 0.013693,
     $     0.013532, 0.013404, 0.013304, 0.01323, 0.013174, 0.013137,
     $     0.013117, 0.013112, 0.013116, 0.013128, 0.013148, 0.013181,
     $     0.013222, 0.013267, 0.013316, 0.01337, 0.013428, 0.013491,
     $     0.013554, 0.013619, 0.013686, 0.013756, 0.013827, 0.013897,
     $     0.013968, 0.014039, 0.014111, 0.014183, 0.014252, 0.014318,
     $     0.014384, 0.01445, 0.014515, 0.014585, 0.014654, 0.014722,
     $     0.01479, 0.014856, 0.014922, 0.014986, 0.01505, 0.015112,
     $     0.015173, 0.015233, 0.015291, 0.015348, 0.015404, 0.015458,
     $     0.015507, 0.015553, 0.015596, 0.015639, 0.01568, 0.015718,
     $     0.015755, 0.015791, 0.015825, 0.015857, 0.015888, 0.015918,
     $     0.015946, 0.015973, 0.015999, 0.016023, 0.016045, 0.016067,
     $     0.016086, 0.016105, 0.01612, 0.016133, 0.016144, 0.016154,
     $     0.016163, 0.016168, 0.016173, 0.016177, 0.016179, 0.01618,
     $     0.01618, 0.016178, 0.016176, 0.016172, 0.016167, 0.016161,
     $     0.016154, 0.016146, 0.016137, 0.016126, 0.016114, 0.0161,
     $     0.016085, 0.016069, 0.016052, 0.015831, 0.015526, 0.015163,
     $     0.014747, 0.014291, 0.013812, 0.013323, 0.01295109, 0.01284,
     $     0.012374, 0.011933, 0.01153, 0.011177, 0.010881, 0.010654/


C     80201, JENDL evaluated data
      nbins(58) = 235
      data ( eg(58,i), i=1, 235) /6229400.0, 6400000.0, 6600000.0,
     $     6800000.0, 7000000.0, 7200000.0, 7400000.0, 7600000.0,
     $     7800000.0, 8000000.0, 8200000.0, 8400000.0, 8600000.0,
     $     8800000.0, 9000000.0, 9200000.0, 9400000.0, 9600000.0,
     $     9800000.0, 1.0D7, 1.02D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7,
     $     1.12D7, 1.14D7, 1.16D7, 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7
     $     , 1.28D7, 1.3D7, 1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7,
     $     1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7
     $     , 1.58D7, 1.6D7, 1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7,
     $     1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7
     $     , 1.88D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7,
     $     2.02D7, 2.04D7, 2.06D7, 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7
     $     , 2.18D7, 2.2D7, 2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7,
     $     2.32D7, 2.34D7, 2.36D7, 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7
     $     , 2.48D7, 2.5D7, 2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7,
     $     2.62D7, 2.64D7, 2.66D7, 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7
     $     , 2.78D7, 2.8D7, 2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7,
     $     2.92D7, 2.94D7, 2.96D7, 2.98D7, 3.0D7, 3.05D7, 3.1D7, 3.15D7,
     $     3.2D7, 3.25D7, 3.3D7, 3.35D7, 3.4D7, 3.45D7, 3.5D7, 3.55D7,
     $     3.6D7, 3.65D7, 3.7D7, 3.75D7, 3.8D7, 3.85D7, 3.9D7, 3.95D7,
     $     4.0D7, 4.05D7, 4.1D7, 4.15D7, 4.2D7, 4.25D7, 4.3D7, 4.35D7,
     $     4.4D7, 4.45D7, 4.5D7, 4.55D7, 4.6D7, 4.65D7, 4.7D7, 4.75D7,
     $     4.8D7, 4.85D7, 4.9D7, 4.95D7, 5.0D7, 5.05D7, 5.1D7, 5.15D7,
     $     5.2D7, 5.25D7, 5.3D7, 5.35D7, 5.4D7, 5.45D7, 5.5D7, 5.55D7,
     $     5.6D7, 5.65D7, 5.7D7, 5.75D7, 5.8D7, 5.85D7, 5.9D7, 5.95D7,
     $     6.0D7, 6.05D7, 6.1D7, 6.15D7, 6.2D7, 6.25D7, 6.3D7, 6.35D7,
     $     6.4D7, 6.45D7, 6.5D7, 6.55D7, 6.6D7, 6.65D7, 6.7D7, 6.75D7,
     $     6.8D7, 6.85D7, 6.9D7, 6.95D7, 7.0D7, 7.05D7, 7.1D7, 7.15D7,
     $     7.2D7, 7.25D7, 7.3D7, 7.35D7, 7.4D7, 7.45D7, 7.5D7, 7.55D7,
     $     7.6D7, 7.65D7, 7.7D7, 7.75D7, 7.8D7, 7.85D7, 7.9D7, 7.95D7,
     $     8.0D7, 8.5D7, 9.0D7, 9.5D7, 1.0D8, 1.05D8, 1.1D8, 1.15D8,
     $     1.2D8, 1.25D8, 1.3D8, 1.35D8, 1.3625D8, 1.4D8, 1.45D8, 1.5D8/

      data ( gdr(58,i), i=1, 235) /0.0, 4.4087D-4, 0.0018471, 0.0037713,
     $     0.0062924, 0.0094601, 0.012964, 0.017289, 0.022698, 0.029356,
     $     0.037415, 0.045361, 0.05417, 0.06405, 0.074877, 0.086439,
     $     0.098464, 0.11068, 0.12292, 0.13517, 0.14765, 0.15915,
     $     0.17095, 0.18319, 0.19806, 0.21645, 0.2393, 0.26747, 0.30161,
     $     0.34176, 0.38701, 0.43348, 0.47884, 0.5191, 0.55057, 0.57056,
     $     0.57617, 0.5732, 0.56328, 0.54653, 0.52253, 0.49349, 0.46732,
     $     0.44053, 0.41491, 0.39111, 0.36789, 0.34635, 0.32727, 0.30898
     $     , 0.29049, 0.27018, 0.25422, 0.24372, 0.23161, 0.21864,
     $     0.20633, 0.19467, 0.18377, 0.17362, 0.16424, 0.15562, 0.14784
     $     , 0.14097, 0.13505, 0.13007, 0.12601, 0.12281, 0.12035,
     $     0.11845, 0.11682, 0.11504, 0.11302, 0.11052, 0.10735, 0.10333
     $     , 0.098562, 0.093335, 0.087834, 0.082254, 0.076771, 0.071518,
     $     0.066581, 0.061946, 0.057676, 0.053789, 0.050264, 0.047074,
     $     0.044191, 0.041583, 0.039223, 0.037084, 0.035142, 0.033377,
     $     0.031768, 0.030345, 0.029024, 0.027816, 0.026708, 0.025692,
     $     0.024757, 0.023896, 0.023103, 0.02237, 0.021693, 0.021066,
     $     0.020486, 0.019948, 0.019438, 0.018961, 0.018516, 0.018102,
     $     0.017717, 0.017357, 0.017021, 0.016708, 0.016415, 0.016137,
     $     0.015884, 0.015648, 0.015134, 0.0147, 0.014339, 0.014041,
     $     0.013797, 0.013591, 0.01342, 0.013285, 0.013181, 0.013104,
     $     0.013046, 0.013006, 0.012982, 0.012974, 0.01298, 0.012993,
     $     0.013015, 0.013046, 0.013086, 0.013134, 0.013185, 0.01324,
     $     0.0133, 0.013365, 0.013433, 0.013502, 0.013571, 0.013643,
     $     0.013717, 0.013792, 0.013866, 0.01394, 0.014015, 0.01409,
     $     0.014166, 0.014237, 0.014306, 0.014375, 0.014443, 0.014511,
     $     0.014583, 0.014653, 0.014723, 0.014792, 0.01486, 0.014927,
     $     0.014992, 0.015057, 0.01512, 0.015182, 0.015243, 0.015303,
     $     0.015361, 0.015418, 0.015473, 0.015522, 0.015568, 0.015613,
     $     0.015656, 0.015698, 0.015737, 0.015775, 0.015811, 0.015846,
     $     0.01588, 0.015912, 0.015943, 0.015972, 0.015999, 0.016026,
     $     0.01605, 0.016074, 0.016096, 0.016116, 0.016135, 0.016151,
     $     0.016164, 0.016176, 0.016186, 0.016195, 0.016201, 0.016206,
     $     0.016209, 0.016211, 0.016212, 0.016211, 0.01621, 0.016207,
     $     0.016203, 0.016198, 0.016192, 0.016185, 0.016177, 0.016167,
     $     0.016157, 0.016144, 0.01613, 0.016114, 0.016098, 0.016081,
     $     0.015863, 0.01556, 0.015198, 0.014783, 0.014327, 0.013847,
     $     0.013357, 0.012873, 0.012406, 0.011966, 0.011563, 0.01147425,
     $     0.011208, 0.010913, 0.010686/


C     80202, JENDL evaluated data
      nbins(59) = 227
      data ( eg(59,i), i=1, 227) /7754400.0, 7800000.0, 8000000.0,
     $     8200000.0, 8400000.0, 8600000.0, 8800000.0, 9000000.0,
     $     9200000.0, 9400000.0, 9600000.0, 9800000.0, 1.0D7, 1.02D7,
     $     1.04D7, 1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7, 1.18D7
     $     , 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7, 1.32D7,
     $     1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7, 1.48D7
     $     , 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7, 1.62D7,
     $     1.64D7, 1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7, 1.78D7
     $     , 1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7, 1.92D7,
     $     1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7, 2.08D7
     $     , 2.1D7, 2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7, 2.22D7,
     $     2.24D7, 2.26D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7, 2.38D7
     $     , 2.4D7, 2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7, 2.52D7,
     $     2.54D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7, 2.68D7
     $     , 2.7D7, 2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7, 2.82D7,
     $     2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7, 2.98D7
     $     , 3.0D7, 3.05D7, 3.1D7, 3.15D7, 3.2D7, 3.25D7, 3.3D7, 3.35D7,
     $     3.4D7, 3.45D7, 3.5D7, 3.55D7, 3.6D7, 3.65D7, 3.7D7, 3.75D7,
     $     3.8D7, 3.85D7, 3.9D7, 3.95D7, 4.0D7, 4.05D7, 4.1D7, 4.15D7,
     $     4.2D7, 4.25D7, 4.3D7, 4.35D7, 4.4D7, 4.45D7, 4.5D7, 4.55D7,
     $     4.6D7, 4.65D7, 4.7D7, 4.75D7, 4.8D7, 4.85D7, 4.9D7, 4.95D7,
     $     5.0D7, 5.05D7, 5.1D7, 5.15D7, 5.2D7, 5.25D7, 5.3D7, 5.35D7,
     $     5.4D7, 5.45D7, 5.5D7, 5.55D7, 5.6D7, 5.65D7, 5.7D7, 5.75D7,
     $     5.8D7, 5.85D7, 5.9D7, 5.95D7, 6.0D7, 6.05D7, 6.1D7, 6.15D7,
     $     6.2D7, 6.25D7, 6.3D7, 6.35D7, 6.4D7, 6.45D7, 6.5D7, 6.55D7,
     $     6.6D7, 6.65D7, 6.7D7, 6.75D7, 6.8D7, 6.85D7, 6.9D7, 6.95D7,
     $     7.0D7, 7.05D7, 7.1D7, 7.15D7, 7.2D7, 7.25D7, 7.3D7, 7.35D7,
     $     7.4D7, 7.45D7, 7.5D7, 7.55D7, 7.6D7, 7.65D7, 7.7D7, 7.75D7,
     $     7.8D7, 7.85D7, 7.9D7, 7.95D7, 8.0D7, 8.5D7, 9.0D7, 9.5D7,
     $     1.0D8, 1.05D8, 1.1D8, 1.15D8, 1.2D8, 1.25D8, 1.3D8, 1.35D8,
     $     1.4D8, 1.45D8, 1.5D8/

      data ( gdr(59,i), i=1, 227) /0.0, 5.4008D-5, 0.0033818, 0.0090198,
     $     0.016816, 0.026655, 0.038071, 0.050064, 0.063572, 0.078044,
     $     0.092908, 0.10708, 0.11867, 0.13016, 0.14204, 0.15405,
     $     0.16616, 0.18141, 0.2009, 0.22576, 0.25706, 0.29473, 0.33645,
     $     0.38374, 0.43408, 0.48345, 0.52714, 0.56144, 0.58281, 0.58898
     $     , 0.58664, 0.57659, 0.55846, 0.54033, 0.51969, 0.49385,
     $     0.4669, 0.43987, 0.41278, 0.38851, 0.36649, 0.34574, 0.32519,
     $     0.30229, 0.28021, 0.26371, 0.25307, 0.24002, 0.22654, 0.21373
     $     , 0.20158, 0.19013, 0.17941, 0.16941, 0.16038, 0.1523,
     $     0.14519, 0.13908, 0.13396, 0.12981, 0.12654, 0.12402, 0.12205
     $     , 0.1203, 0.11856, 0.11653, 0.11392, 0.11052, 0.10622,
     $     0.10123, 0.0958, 0.090124, 0.084406, 0.078794, 0.073438,
     $     0.068423, 0.063788, 0.059545, 0.055685, 0.052145, 0.048934,
     $     0.046033, 0.043412, 0.04104, 0.038892, 0.036943, 0.035171,
     $     0.033556, 0.032048, 0.030667, 0.029399, 0.028229, 0.027139,
     $     0.026132, 0.0252, 0.024337, 0.023537, 0.022793, 0.022101,
     $     0.021456, 0.020855, 0.020295, 0.019771, 0.019282, 0.018816,
     $     0.018377, 0.017961, 0.017574, 0.017211, 0.016869, 0.016547,
     $     0.016244, 0.015959, 0.015428, 0.014976, 0.01459, 0.014268,
     $     0.014002, 0.013784, 0.013602, 0.013449, 0.013328, 0.013235,
     $     0.013174, 0.013128, 0.013097, 0.013082, 0.01308, 0.013093,
     $     0.013113, 0.013139, 0.013174, 0.013217, 0.013269, 0.013324,
     $     0.013381, 0.013442, 0.013508, 0.013577, 0.013647, 0.013716,
     $     0.013788, 0.013861, 0.013936, 0.01401, 0.014083, 0.014157,
     $     0.014231, 0.014305, 0.014376, 0.014444, 0.014511, 0.014578,
     $     0.014649, 0.01472, 0.014789, 0.014858, 0.014925, 0.014992,
     $     0.015058, 0.015122, 0.015185, 0.015247, 0.015308, 0.015367,
     $     0.015425, 0.015482, 0.015537, 0.01559, 0.015639, 0.015684,
     $     0.015727, 0.015769, 0.015807, 0.015844, 0.01588, 0.015914,
     $     0.015947, 0.015978, 0.016008, 0.016036, 0.016063, 0.016088,
     $     0.016112, 0.016135, 0.016156, 0.016176, 0.016194, 0.016212,
     $     0.016226, 0.016237, 0.016247, 0.016255, 0.016262, 0.016267,
     $     0.016271, 0.016273, 0.016275, 0.016275, 0.016274, 0.016272,
     $     0.016268, 0.016264, 0.016258, 0.016251, 0.016243, 0.016234,
     $     0.016224, 0.016213, 0.0162, 0.016185, 0.016169, 0.016152,
     $     0.01593, 0.015624, 0.015259, 0.014841, 0.014382, 0.0139,
     $     0.013408, 0.012921, 0.012452, 0.01201, 0.011605, 0.011248,
     $     0.010951, 0.010723/


C     80204, JENDL evaluated data
      nbins(60) = 228
      data ( eg(60,i), i=1, 228) /7495400.0, 7600000.0, 7800000.0,
     $     8000000.0, 8200000.0, 8400000.0, 8600000.0, 8800000.0,
     $     9000000.0, 9200000.0, 9400000.0, 9600000.0, 9800000.0, 1.0D7,
     $     1.02D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7, 1.16D7
     $     , 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.28D7, 1.3D7,
     $     1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7, 1.42D7, 1.44D7, 1.46D7
     $     , 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7, 1.58D7, 1.6D7,
     $     1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7, 1.72D7, 1.74D7, 1.76D7
     $     , 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7,
     $     1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7, 2.06D7
     $     , 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7, 2.18D7, 2.2D7,
     $     2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7, 2.32D7, 2.34D7, 2.36D7
     $     , 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7, 2.48D7, 2.5D7,
     $     2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.64D7, 2.66D7
     $     , 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7,
     $     2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7, 2.96D7
     $     , 2.98D7, 3.0D7, 3.05D7, 3.1D7, 3.15D7, 3.2D7, 3.25D7, 3.3D7,
     $     3.35D7, 3.4D7, 3.45D7, 3.5D7, 3.55D7, 3.6D7, 3.65D7, 3.7D7,
     $     3.75D7, 3.8D7, 3.85D7, 3.9D7, 3.95D7, 4.0D7, 4.05D7, 4.1D7,
     $     4.15D7, 4.2D7, 4.25D7, 4.3D7, 4.35D7, 4.4D7, 4.45D7, 4.5D7,
     $     4.55D7, 4.6D7, 4.65D7, 4.7D7, 4.75D7, 4.8D7, 4.85D7, 4.9D7,
     $     4.95D7, 5.0D7, 5.05D7, 5.1D7, 5.15D7, 5.2D7, 5.25D7, 5.3D7,
     $     5.35D7, 5.4D7, 5.45D7, 5.5D7, 5.55D7, 5.6D7, 5.65D7, 5.7D7,
     $     5.75D7, 5.8D7, 5.85D7, 5.9D7, 5.95D7, 6.0D7, 6.05D7, 6.1D7,
     $     6.15D7, 6.2D7, 6.25D7, 6.3D7, 6.35D7, 6.4D7, 6.45D7, 6.5D7,
     $     6.55D7, 6.6D7, 6.65D7, 6.7D7, 6.75D7, 6.8D7, 6.85D7, 6.9D7,
     $     6.95D7, 7.0D7, 7.05D7, 7.1D7, 7.15D7, 7.2D7, 7.25D7, 7.3D7,
     $     7.35D7, 7.4D7, 7.45D7, 7.5D7, 7.55D7, 7.6D7, 7.65D7, 7.7D7,
     $     7.75D7, 7.8D7, 7.85D7, 7.9D7, 7.95D7, 8.0D7, 8.5D7, 9.0D7,
     $     9.5D7, 1.0D8, 1.05D8, 1.1D8, 1.15D8, 1.2D8, 1.25D8, 1.3D8,
     $     1.35D8, 1.4D8, 1.45D8, 1.5D8/

      data ( gdr(60,i), i=1, 228) /0.0, 2.962D-4, 0.0037393, 0.0087541,
     $     0.015518, 0.024035, 0.033712, 0.04453, 0.057027, 0.070813,
     $     0.085372, 0.098699, 0.11048, 0.12203, 0.13362, 0.14577,
     $     0.15783, 0.17079, 0.18711, 0.20785, 0.2341, 0.26517, 0.30126,
     $     0.34366, 0.39109, 0.44074, 0.48842, 0.52958, 0.56095, 0.57811
     $     , 0.58216, 0.5779, 0.56975, 0.55606, 0.55422, 0.53473,
     $     0.50852, 0.48002, 0.4528, 0.42764, 0.40449, 0.38081, 0.35674,
     $     0.33332, 0.31008, 0.28684, 0.26954, 0.25865, 0.24419, 0.2305,
     $     0.21742, 0.20493, 0.19308, 0.18198, 0.17172, 0.16238, 0.15402
     $     , 0.14667, 0.14042, 0.13521, 0.13097, 0.12755, 0.12481,
     $     0.12258, 0.12058, 0.11853, 0.11613, 0.11316, 0.10953, 0.10491
     $     , 0.099679, 0.093985, 0.088138, 0.082329, 0.076713, 0.071399,
     $     0.066451, 0.061897, 0.057739, 0.053963, 0.050543, 0.04745,
     $     0.044654, 0.042096, 0.03977, 0.03766, 0.035743, 0.033998,
     $     0.032407, 0.030956, 0.029628, 0.02841, 0.02729, 0.02626,
     $     0.025311, 0.024431, 0.023603, 0.022837, 0.022125, 0.021464,
     $     0.020849, 0.020275, 0.019742, 0.019245, 0.018781, 0.018348,
     $     0.017943, 0.017565, 0.017203, 0.016861, 0.01654, 0.016238,
     $     0.015955, 0.015689, 0.015176, 0.014751, 0.014394, 0.014086,
     $     0.013833, 0.013627, 0.013463, 0.013331, 0.013219, 0.013136,
     $     0.01308, 0.013045, 0.013026, 0.013016, 0.013021, 0.013038,
     $     0.013066, 0.013102, 0.013141, 0.013188, 0.013243, 0.013304,
     $     0.013368, 0.013432, 0.0135, 0.013572, 0.013646, 0.013722,
     $     0.013796, 0.013871, 0.013948, 0.014026, 0.014104, 0.01418,
     $     0.014256, 0.014332, 0.014408, 0.014482, 0.01455, 0.014619,
     $     0.014691, 0.014762, 0.014832, 0.014902, 0.01497, 0.015037,
     $     0.015104, 0.015169, 0.015233, 0.015295, 0.015357, 0.015417,
     $     0.015475, 0.015532, 0.015588, 0.015643, 0.015696, 0.015745,
     $     0.015789, 0.015831, 0.01587, 0.015908, 0.015945, 0.01598,
     $     0.016013, 0.016046, 0.016076, 0.016105, 0.016133, 0.01616,
     $     0.016184, 0.016208, 0.01623, 0.01625, 0.01627, 0.016287,
     $     0.016304, 0.016318, 0.016328, 0.016337, 0.016343, 0.016348,
     $     0.016352, 0.016354, 0.016355, 0.016355, 0.016354, 0.016352,
     $     0.016348, 0.016343, 0.016338, 0.016331, 0.016322, 0.016313,
     $     0.016303, 0.016292, 0.016279, 0.016265, 0.016249, 0.016232,
     $     0.016012, 0.015708, 0.015343, 0.014925, 0.014464, 0.013982,
     $     0.013487, 0.012999, 0.012526, 0.012082, 0.011674, 0.011316,
     $     0.011017, 0.010787/


C     82206, JENDL evaluated data
      nbins(61) = 361
      data ( eg(61,i), i=1, 361) /1.0D-5, 100000.0, 776180.0, 1043100.0
     $     , 2691400.0, 3637320.0, 3744240.0, 4637640.0, 4710130.0,
     $     4758730.0, 5211070.0, 5826460.0, 5849390.0, 6397440.0,
     $     7254760.0, 7507750.0, 7586130.0, 7932880.0, 8088190.0,
     $     8150000.0, 8200000.0, 8400000.0, 8600000.0, 8800000.0,
     $     9000000.0, 9063090.0, 9200000.0, 9321320.0, 9400000.0,
     $     9597240.0, 9600000.0, 9800000.0, 1.0D7, 1.02D7, 1.04D7,
     $     1.06D7, 1.06108D7, 1.08D7, 1.1D7, 1.11622D7, 1.12D7,
     $     1.12028D7, 1.14D7, 1.16D7, 1.16582D7, 1.18D7, 1.18745D7,
     $     1.2D7, 1.2162D7, 1.22D7, 1.24D7, 1.25765D7, 1.26D7, 1.26285D7
     $     , 1.27057D7, 1.28D7, 1.29751D7, 1.3D7, 1.32D7, 1.34D7,
     $     1.34484D7, 1.36D7, 1.36713D7, 1.38D7, 1.4D7, 1.40984D7,
     $     1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.48197D7, 1.5D7, 1.51488D7,
     $     1.52D7, 1.54D7, 1.54134D7, 1.55566D7, 1.56D7, 1.5726D7,
     $     1.58D7, 1.6D7, 1.60585D7, 1.62D7, 1.64D7, 1.66D7, 1.68D7,
     $     1.7D7, 1.72D7, 1.74D7, 1.75975D7, 1.76D7, 1.78D7, 1.8D7,
     $     1.82D7, 1.83614D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7, 1.92D7,
     $     1.92416D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.02D7, 2.03476D7,
     $     2.04D7, 2.06D7, 2.08D7, 2.09006D7, 2.1D7, 2.12D7, 2.13572D7,
     $     2.14D7, 2.16D7, 2.18D7, 2.18546D7, 2.2D7, 2.22D7, 2.24D7,
     $     2.2508D7, 2.26D7, 2.28D7, 2.29793D7, 2.3D7, 2.32D7, 2.3214D7,
     $     2.34D7, 2.36D7, 2.37309D7, 2.38D7, 2.4D7, 2.41743D7, 2.42D7,
     $     2.44D7, 2.44281D7, 2.46D7, 2.47049D7, 2.48D7, 2.5D7, 2.52D7,
     $     2.53944D7, 2.54D7, 2.56D7, 2.56936D7, 2.58D7, 2.6D7, 2.62D7,
     $     2.64D7, 2.65177D7, 2.66D7, 2.66782D7, 2.68D7, 2.6952D7, 2.7D7
     $     , 2.72D7, 2.72514D7, 2.73875D7, 2.74D7, 2.76D7, 2.78D7, 2.8D7
     $     , 2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.90616D7, 2.92D7,
     $     2.94D7, 2.96D7, 2.97773D7, 2.98D7, 3.0D7, 3.0138D7, 3.02D7,
     $     3.03585D7, 3.04D7, 3.06D7, 3.08D7, 3.1D7, 3.12D7, 3.14D7,
     $     3.16D7, 3.18D7, 3.2D7, 3.20109D7, 3.21291D7, 3.22D7, 3.24D7,
     $     3.26D7, 3.28D7, 3.3D7, 3.32D7, 3.34D7, 3.34385D7, 3.36D7,
     $     3.38D7, 3.4D7, 3.40616D7, 3.40953D7, 3.41046D7, 3.42D7,
     $     3.44D7, 3.46D7, 3.48D7, 3.5D7, 3.52D7, 3.54D7, 3.56D7, 3.58D7
     $     , 3.6D7, 3.62D7, 3.64D7, 3.66D7, 3.68D7, 3.7D7, 3.72D7,
     $     3.74D7, 3.76D7, 3.77935D7, 3.78D7, 3.8D7, 3.80402D7, 3.82D7,
     $     3.84D7, 3.86D7, 3.86002D7, 3.88D7, 3.88151D7, 3.8864D7, 3.9D7
     $     , 3.91246D7, 3.92D7, 3.93497D7, 3.94D7, 3.95064D7, 3.96D7,
     $     3.98D7, 4.0D7, 4.07674D7, 4.161D7, 4.2D7, 4.34465D7,
     $     4.34619D7, 4.4D7, 4.49649D7, 4.56226D7, 4.57215D7, 4.59749D7,
     $     4.6D7, 4.71025D7, 4.75525D7, 4.8D7, 4.92691D7, 5.0D7,
     $     5.02355D7, 5.10956D7, 5.2D7, 5.25237D7, 5.30189D7, 5.30896D7,
     $     5.4D7, 5.50651D7, 5.56175D7, 5.6D7, 5.6973D7, 5.79714D7,
     $     5.8D7, 5.81856D7, 5.92313D7, 6.0D7, 6.11131D7, 6.113D7, 6.2D7
     $     , 6.22575D7, 6.22712D7, 6.29165D7, 6.4D7, 6.54157D7, 6.6D7,
     $     6.61946D7, 6.65987D7, 6.70447D7, 6.8D7, 7.0D7, 7.01017D7,
     $     7.0136D7, 7.16465D7, 7.2D7, 7.36799D7, 7.39757D7, 7.4D7,
     $     7.48612D7, 7.53D7, 7.6D7, 7.8D7, 7.90938D7, 8.0D7, 8.2D7,
     $     8.20702D7, 8.29242D7, 8.34727D7, 8.36325D7, 8.4D7, 8.6D7,
     $     8.68894D7, 8.8D7, 8.88052D7, 8.93465D7, 9.0D7, 9.06607D7,
     $     9.12425D7, 9.14062D7, 9.2D7, 9.4D7, 9.6D7, 9.64065D7,
     $     9.71047D7, 9.8D7, 9.87545D7, 9.88977D7, 1.0D8, 1.02D8, 1.04D8
     $     , 1.05638D8, 1.06D8, 1.06447D8, 1.0649D8, 1.07153D8, 1.08D8,
     $     1.08647D8, 1.1D8, 1.12D8, 1.14D8, 1.14222D8, 1.1526D8,
     $     1.15468D8, 1.16D8, 1.17246D8, 1.18D8, 1.2D8, 1.22D8,
     $     1.22805D8, 1.24D8, 1.24563D8, 1.24689D8, 1.26D8, 1.2722D8,
     $     1.28D8, 1.3D8, 1.32D8, 1.32617D8, 1.32694D8, 1.32855D8,
     $     1.34D8, 1.35043D8, 1.36D8, 1.38D8, 1.4D8/

      data ( gdr(61,i), i=1, 361) /0.0, 0.0, 9.66654D-4, 0.00129908,
     $     0.00335187, 0.00452992, 0.00466307, 0.00577572, 0.00586599,
     $     0.00592651, 0.00648987, 0.00725627, 0.00728482, 0.00796736,
     $     0.00903506, 0.00935014, 0.00944776, 0.0098796, 0.010073,
     $     0.01015, 0.0184, 0.04094, 0.05219, 0.05635, 0.05864,
     $     0.0600525, 0.0632, 0.0671655, 0.06984, 0.0772126, 0.07732,
     $     0.08525, 0.09411, 0.1044, 0.1162, 0.1297, 0.13049, 0.145,
     $     0.1626, 0.178746, 0.1827, 0.182998, 0.2054, 0.2311, 0.239151,
     $     0.2598, 0.27121, 0.2914, 0.318933, 0.3257, 0.3616, 0.393482,
     $     0.3979, 0.402729, 0.416046, 0.4328, 0.459677, 0.4636, 0.4879,
     $     0.5037, 0.505032, 0.5092, 0.507441, 0.5043, 0.4903, 0.479518,
     $     0.4687, 0.4418, 0.412, 0.3811, 0.377963, 0.3506, 0.328675,
     $     0.3215, 0.2943, 0.292538, 0.274513, 0.2693, 0.25467, 0.2465,
     $     0.226, 0.220403, 0.2075, 0.191, 0.1762, 0.1629, 0.151, 0.1402
     $     , 0.1307, 0.122103, 0.122, 0.1142, 0.1071, 0.1007, 0.0959396,
     $     0.09484, 0.0895, 0.08463, 0.08019, 0.0761, 0.0753013, 0.07235
     $     , 0.06888, 0.0657, 0.06274, 0.06001, 0.0581229, 0.05747,
     $     0.05511, 0.05291, 0.0518663, 0.05086, 0.04895, 0.0475362,
     $     0.04716, 0.04549, 0.04391, 0.0434995, 0.04243, 0.04105,
     $     0.03974, 0.0390641, 0.0385, 0.03734, 0.0363517, 0.03624,
     $     0.0352, 0.0351309, 0.03423, 0.0333, 0.0327205, 0.03242,
     $     0.03158, 0.0308815, 0.03078, 0.03002, 0.0299174, 0.0293,
     $     0.0289302, 0.0286, 0.02794, 0.02732, 0.0267367, 0.02672,
     $     0.02615, 0.0258859, 0.02559, 0.02506, 0.02453, 0.02405,
     $     0.0237956, 0.02362, 0.0234624, 0.02322, 0.0228388, 0.02272,
     $     0.02205, 0.021857, 0.0213554, 0.02131, 0.02065, 0.02012,
     $     0.01974, 0.01947, 0.0193, 0.01922, 0.01921, 0.01926,
     $     0.0192908, 0.01936, 0.01948, 0.01962, 0.0197353, 0.01975,
     $     0.01987, 0.019939, 0.01997, 0.0200176, 0.02003, 0.02007,
     $     0.02009, 0.02009, 0.02008, 0.02005, 0.02, 0.01995, 0.0199,
     $     0.0198973, 0.0198677, 0.01985, 0.01979, 0.01975, 0.01971,
     $     0.01967, 0.01963, 0.0196, 0.0195961, 0.01958, 0.01955,
     $     0.01953, 0.0195207, 0.0195157, 0.0195143, 0.0195, 0.01948,
     $     0.01947, 0.01945, 0.01943, 0.01941, 0.01939, 0.01936, 0.01933
     $     , 0.01931, 0.01929, 0.01926, 0.01923, 0.0192, 0.01917,
     $     0.01914, 0.01911, 0.01909, 0.0190513, 0.01905, 0.01902,
     $     0.019014, 0.01899, 0.01896, 0.01893, 0.01893, 0.01891,
     $     0.0189077, 0.0189004, 0.01888, 0.018855, 0.01884, 0.0188175,
     $     0.01881, 0.018794, 0.01878, 0.01875, 0.01872, 0.0186026,
     $     0.018477, 0.01842, 0.0181865, 0.0181841, 0.0181, 0.0179481,
     $     0.0178471, 0.017832, 0.0177938, 0.01779, 0.0176504, 0.0175947
     $     , 0.01754, 0.0173926, 0.01731, 0.017275, 0.0171491, 0.01702,
     $     0.0169425, 0.0168703, 0.01686, 0.01673, 0.0166007, 0.016535,
     $     0.01649, 0.0163866, 0.0162829, 0.01628, 0.0162572, 0.0161309,
     $     0.01604, 0.0159276, 0.0159259, 0.01584, 0.0158137, 0.0158124,
     $     0.0157472, 0.01564, 0.0155118, 0.01546, 0.0154431, 0.0154084,
     $     0.0153703, 0.01529, 0.01507, 0.0150591, 0.0150554, 0.0148965,
     $     0.01486, 0.0147168, 0.014692, 0.01469, 0.0146204, 0.0145853,
     $     0.01453, 0.01432, 0.0142207, 0.01414, 0.01405, 0.014045,
     $     0.0139847, 0.0139465, 0.0139354, 0.01391, 0.01375, 0.0136962,
     $     0.01363, 0.013569, 0.0135285, 0.01348, 0.0134399, 0.013405,
     $     0.0133952, 0.01336, 0.01322, 0.01309, 0.0130694, 0.0130344,
     $     0.01299, 0.0129519, 0.0129448, 0.01289, 0.0128, 0.01269,
     $     0.0125832, 0.01256, 0.0125397, 0.0125377, 0.0125078, 0.01247,
     $     0.0124439, 0.01239, 0.01232, 0.01222, 0.0122099, 0.0121631,
     $     0.0121537, 0.01213, 0.0120862, 0.01206, 0.01196, 0.01188,
     $     0.0118597, 0.01183, 0.011813, 0.0118092, 0.01177, 0.0117271,
     $     0.0117, 0.01169, 0.01165, 0.0116283, 0.0116255, 0.0116199,
     $     0.01158, 0.0115538, 0.01153, 0.0115, 0.01146/


C     82207, JENDL evaluated data
      nbins(62) = 368
      data ( eg(62,i), i=1, 368) /1.0D-5, 654.0, 1364820.0, 1681130.0,
     $     2803780.0, 3015690.0, 3891350.0, 4143570.0, 4635590.0,
     $     4795630.0, 4821190.0, 5769210.0, 6234410.0, 6491230.0,
     $     6497420.0, 6737790.0, 6800000.0, 7000000.0, 7200000.0,
     $     7400000.0, 7488980.0, 7600000.0, 7800000.0, 8000000.0,
     $     8200000.0, 8400000.0, 8403300.0, 8600000.0, 8800000.0,
     $     9000000.0, 9061770.0, 9200000.0, 9377870.0, 9400000.0,
     $     9600000.0, 9730030.0, 9800000.0, 9978100.0, 1.0D7, 1.00986D7,
     $     1.0137D7, 1.02D7, 1.04D7, 1.06D7, 1.07265D7, 1.08D7, 1.1D7,
     $     1.12D7, 1.14D7, 1.15945D7, 1.16D7, 1.16181D7, 1.1768D7,
     $     1.18D7, 1.19691D7, 1.2D7, 1.22D7, 1.23821D7, 1.24D7, 1.26D7,
     $     1.2691D7, 1.28D7, 1.28986D7, 1.3D7, 1.3057D7, 1.32D7, 1.34D7,
     $     1.36D7, 1.38D7, 1.4D7, 1.4162D7, 1.42D7, 1.43022D7, 1.44D7,
     $     1.46D7, 1.47311D7, 1.47414D7, 1.48D7, 1.4826D7, 1.5D7, 1.52D7
     $     , 1.54D7, 1.56D7, 1.58D7, 1.59787D7, 1.6D7, 1.62D7, 1.64D7,
     $     1.65982D7, 1.66D7, 1.68D7, 1.7D7, 1.70851D7, 1.72D7, 1.74D7,
     $     1.76D7, 1.78D7, 1.7823D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7,
     $     1.88D7, 1.88066D7, 1.9D7, 1.90571D7, 1.92D7, 1.94D7,
     $     1.95067D7, 1.96D7, 1.96827D7, 1.98D7, 2.0D7, 2.02D7, 2.04D7,
     $     2.06D7, 2.08D7, 2.0971D7, 2.1D7, 2.12D7, 2.13469D7, 2.14D7,
     $     2.15575D7, 2.16D7, 2.18D7, 2.2D7, 2.22D7, 2.24D7, 2.26D7,
     $     2.26196D7, 2.28D7, 2.28397D7, 2.3D7, 2.32D7, 2.34D7,
     $     2.34344D7, 2.3567D7, 2.36D7, 2.36975D7, 2.38D7, 2.4D7, 2.42D7
     $     , 2.43603D7, 2.44D7, 2.45018D7, 2.46D7, 2.46277D7, 2.48D7,
     $     2.5D7, 2.52D7, 2.54D7, 2.54278D7, 2.56D7, 2.58D7, 2.6D7,
     $     2.62D7, 2.64D7, 2.64843D7, 2.65907D7, 2.66D7, 2.68D7, 2.7D7,
     $     2.72D7, 2.74D7, 2.76D7, 2.78D7, 2.78698D7, 2.7894D7, 2.8D7,
     $     2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7, 2.92D7, 2.94D7,
     $     2.94753D7, 2.96D7, 2.98D7, 2.99518D7, 3.0D7, 3.02D7, 3.04D7,
     $     3.06D7, 3.08D7, 3.09792D7, 3.1D7, 3.12D7, 3.14D7, 3.16D7,
     $     3.1766D7, 3.18D7, 3.2D7, 3.21057D7, 3.22D7, 3.22414D7,
     $     3.2315D7, 3.24D7, 3.26D7, 3.28D7, 3.3D7, 3.32D7, 3.33036D7,
     $     3.34D7, 3.36D7, 3.38D7, 3.383D7, 3.4D7, 3.42D7, 3.44D7,
     $     3.46D7, 3.48D7, 3.5D7, 3.52D7, 3.54D7, 3.56D7, 3.58D7,
     $     3.58938D7, 3.6D7, 3.62D7, 3.63337D7, 3.64D7, 3.66D7, 3.68D7,
     $     3.68758D7, 3.68766D7, 3.7D7, 3.70721D7, 3.71347D7, 3.72D7,
     $     3.74D7, 3.76D7, 3.78D7, 3.78651D7, 3.8D7, 3.82D7, 3.84D7,
     $     3.86D7, 3.88D7, 3.9D7, 3.92D7, 3.94D7, 3.96D7, 3.98D7,
     $     3.99775D7, 3.99967D7, 4.0D7, 4.07197D7, 4.07658D7, 4.2D7,
     $     4.29578D7, 4.35019D7, 4.37125D7, 4.37628D7, 4.4D7, 4.5442D7,
     $     4.56018D7, 4.6D7, 4.73142D7, 4.8D7, 4.94814D7, 4.98291D7,
     $     5.0D7, 5.00264D7, 5.11355D7, 5.19415D7, 5.2D7, 5.26163D7,
     $     5.27127D7, 5.29748D7, 5.34439D7, 5.4D7, 5.6D7, 5.60844D7,
     $     5.729D7, 5.73248D7, 5.8D7, 6.0D7, 6.00063D7, 6.02921D7,
     $     6.10642D7, 6.11174D7, 6.18029D7, 6.2D7, 6.32717D7, 6.38474D7,
     $     6.4D7, 6.5198D7, 6.6D7, 6.64362D7, 6.8D7, 6.9009D7, 6.94814D7
     $     , 6.96689D7, 6.98495D7, 7.0D7, 7.16029D7, 7.2D7, 7.2445D7,
     $     7.35027D7, 7.4D7, 7.6D7, 7.60712D7, 7.8D7, 7.83843D7,
     $     7.95608D7, 8.0D7, 8.0613D7, 8.2D7, 8.26746D7, 8.33475D7,
     $     8.4D7, 8.4645D7, 8.58316D7, 8.6D7, 8.65219D7, 8.8D7, 9.0D7,
     $     9.03254D7, 9.08248D7, 9.2D7, 9.20817D7, 9.4D7, 9.45469D7,
     $     9.55429D7, 9.6D7, 9.61901D7, 9.8D7, 9.83692D7, 1.0D8, 1.02D8,
     $     1.03144D8, 1.03192D8, 1.03768D8, 1.04D8, 1.06D8, 1.06271D8,
     $     1.08D8, 1.08114D8, 1.1D8, 1.12D8, 1.12619D8, 1.13185D8,
     $     1.14D8, 1.14867D8, 1.16D8, 1.16399D8, 1.18D8, 1.19971D8,
     $     1.2D8, 1.2096D8, 1.21648D8, 1.22D8, 1.24D8, 1.26D8, 1.26356D8
     $     , 1.28D8, 1.29633D8, 1.3D8, 1.30185D8, 1.31301D8, 1.32D8,
     $     1.34D8, 1.34139D8, 1.36D8, 1.38D8, 1.39354D8, 1.4D8/

      data ( gdr(62,i), i=1, 368) /0.0, 4.82229D-7, 0.00100636,
     $     0.00123958, 0.00206737, 0.00222363, 0.0028693, 0.00305527,
     $     0.00341807, 0.00353607, 0.00355492, 0.00425394, 0.00459696,
     $     0.00478633, 0.00479089, 0.00496813, 0.005014, 0.01614,
     $     0.02041, 0.02126, 0.0217331, 0.02233, 0.025, 0.02842, 0.03183
     $     , 0.03533, 0.0353919, 0.03923, 0.04369, 0.04871, 0.0503907,
     $     0.05432, 0.0598866, 0.06061, 0.06771, 0.0728533, 0.07575,
     $     0.0838469, 0.08489, 0.0898991, 0.0919125, 0.0953, 0.1072,
     $     0.1207, 0.13026, 0.1361, 0.1537, 0.1739, 0.1968, 0.221847,
     $     0.2226, 0.225082, 0.246584, 0.2514, 0.277906, 0.283, 0.3172,
     $     0.349393, 0.3527, 0.3882, 0.40312, 0.4216, 0.435527, 0.4502,
     $     0.45618, 0.4714, 0.4836, 0.4853, 0.4768, 0.4598, 0.440569,
     $     0.4362, 0.42166, 0.4083, 0.3783, 0.358178, 0.356646, 0.3481,
     $     0.344123, 0.3188, 0.2912, 0.2657, 0.2425, 0.2215, 0.20461,
     $     0.2027, 0.1858, 0.1707, 0.157317, 0.1572, 0.1452, 0.1344,
     $     0.130167, 0.1247, 0.1159, 0.1081, 0.1011, 0.100328, 0.09461,
     $     0.08879, 0.0835, 0.07867, 0.07426, 0.0741224, 0.07022,
     $     0.069136, 0.06651, 0.06309, 0.0613856, 0.05994, 0.0587107,
     $     0.05702, 0.05432, 0.05182, 0.04949, 0.04731, 0.0453,
     $     0.0436697, 0.0434, 0.04163, 0.0404029, 0.03997, 0.0387511,
     $     0.03843, 0.03696, 0.03558, 0.03428, 0.03305, 0.0319,
     $     0.0317902, 0.0308, 0.0305921, 0.02977, 0.02878, 0.02784,
     $     0.0276845, 0.0270945, 0.02695, 0.0265414, 0.02612, 0.0253,
     $     0.0245, 0.0239045, 0.02376, 0.023447, 0.02315, 0.0230826,
     $     0.02267, 0.02228, 0.02196, 0.02167, 0.0216378, 0.02144,
     $     0.02123, 0.02106, 0.02092, 0.02081, 0.020772, 0.0207242,
     $     0.02072, 0.02065, 0.02059, 0.02054, 0.0205, 0.02047, 0.02043,
     $     0.020416, 0.0204111, 0.02039, 0.02035, 0.02032, 0.02028,
     $     0.02024, 0.02022, 0.02018, 0.02015, 0.0201387, 0.02012,
     $     0.02008, 0.0200572, 0.02005, 0.02002, 0.01999, 0.01996,
     $     0.01993, 0.0199031, 0.0199, 0.01987, 0.01983, 0.01981,
     $     0.0197851, 0.01978, 0.01976, 0.0197494, 0.01974, 0.0197358,
     $     0.0197285, 0.01972, 0.0197, 0.01968, 0.01967, 0.01965,
     $     0.0196448, 0.01964, 0.01963, 0.01961, 0.0196085, 0.0196,
     $     0.01959, 0.01957, 0.01956, 0.01954, 0.01952, 0.0195, 0.01948,
     $     0.01945, 0.01943, 0.0194112, 0.01939, 0.01936, 0.0193399,
     $     0.01933, 0.0193, 0.01927, 0.0192586, 0.0192585, 0.01924,
     $     0.0192292, 0.0192198, 0.01921, 0.01917, 0.01913, 0.0191,
     $     0.0190902, 0.01907, 0.01903, 0.019, 0.01897, 0.01892, 0.01889
     $     , 0.01886, 0.01882, 0.01879, 0.01875, 0.0187145, 0.0187107,
     $     0.01871, 0.0185924, 0.018585, 0.01839, 0.0182537, 0.0181781,
     $     0.0181492, 0.0181423, 0.01811, 0.017921, 0.0179006, 0.01785,
     $     0.0176708, 0.01758, 0.0173935, 0.0173508, 0.01733, 0.0173255,
     $     0.0171403, 0.0170094, 0.017, 0.0169247, 0.016913, 0.0168815,
     $     0.0168255, 0.01676, 0.01653, 0.0165183, 0.0163543, 0.0163496,
     $     0.01626, 0.01603, 0.0160294, 0.0160002, 0.0159225, 0.0159171,
     $     0.0158493, 0.01583, 0.0157146, 0.0156635, 0.01565, 0.0155352,
     $     0.01546, 0.0154134, 0.01525, 0.015148, 0.0151011, 0.0150825,
     $     0.0150648, 0.01505, 0.0148971, 0.01486, 0.0148193, 0.0147241,
     $     0.01468, 0.01449, 0.0144831, 0.0143, 0.0142688, 0.0141746,
     $     0.01414, 0.0140966, 0.014, 0.0139454, 0.0138916, 0.01384,
     $     0.0137911, 0.0137024, 0.01369, 0.0136557, 0.01356, 0.01346,
     $     0.0134319, 0.0133892, 0.01329, 0.0132859, 0.01319, 0.0131624,
     $     0.0131126, 0.01309, 0.0130775, 0.01296, 0.0129432, 0.01287,
     $     0.01276, 0.0127025, 0.0127, 0.0126715, 0.01266, 0.01256,
     $     0.0125476, 0.01247, 0.0124648, 0.01238, 0.01227, 0.012245,
     $     0.0122224, 0.01219, 0.0121594, 0.01212, 0.0121039, 0.01204,
     $     0.0119612, 0.01196, 0.011931, 0.0119105, 0.0119, 0.01185,
     $     0.01178, 0.0117674, 0.01171, 0.0116691, 0.01166, 0.0116563,
     $     0.0116339, 0.01162, 0.01157, 0.0115665, 0.01152, 0.0115,
     $     0.0114728, 0.01146/


C     82208, JENDL evaluated data
      nbins(63) = 369
      data ( eg(63,i), i=1, 369) /1.0D-5, 1.0, 10000.0, 86750.0,
     $     202512.0, 2148130.0, 2541950.0, 3698900.0, 4446840.0,
     $     4511910.0, 5678540.0, 5703070.0, 5777600.0, 6217560.0,
     $     6404780.0, 6633370.0, 7367820.0, 7450000.0, 7500000.0,
     $     7600000.0, 7800000.0, 7866350.0, 7979570.0, 7993320.0,
     $     8000000.0, 8008140.0, 8200000.0, 8400000.0, 8600000.0,
     $     8681320.0, 8800000.0, 9000000.0, 9200000.0, 9400000.0,
     $     9600000.0, 9762160.0, 9800000.0, 1.0D7, 1.02D7, 1.02549D7,
     $     1.0283D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7, 1.12D7, 1.14D7,
     $     1.1434D7, 1.14995D7, 1.16D7, 1.18D7, 1.19147D7, 1.19956D7,
     $     1.2D7, 1.22D7, 1.24D7, 1.26D7, 1.26322D7, 1.26721D7, 1.28D7,
     $     1.28785D7, 1.28954D7, 1.3D7, 1.32D7, 1.34D7, 1.36D7,
     $     1.37767D7, 1.38D7, 1.4D7, 1.41056D7, 1.4177D7, 1.42D7,
     $     1.43912D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7, 1.51084D7, 1.52D7,
     $     1.53816D7, 1.54D7, 1.54336D7, 1.56D7, 1.56091D7, 1.56771D7,
     $     1.58D7, 1.6D7, 1.61942D7, 1.62D7, 1.64D7, 1.66D7, 1.68D7,
     $     1.7D7, 1.72D7, 1.74D7, 1.754D7, 1.76D7, 1.76064D7, 1.78D7,
     $     1.8D7, 1.81851D7, 1.82D7, 1.84D7, 1.86D7, 1.88D7, 1.9D7,
     $     1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.01947D7, 2.02D7,
     $     2.03693D7, 2.04D7, 2.06D7, 2.08D7, 2.1D7, 2.12D7, 2.14D7,
     $     2.16D7, 2.18D7, 2.2D7, 2.21254D7, 2.21938D7, 2.22D7,
     $     2.22437D7, 2.24D7, 2.26D7, 2.26096D7, 2.28D7, 2.3D7, 2.32D7,
     $     2.33315D7, 2.34D7, 2.36D7, 2.38D7, 2.38413D7, 2.4D7,
     $     2.40094D7, 2.42D7, 2.44D7, 2.45547D7, 2.45802D7, 2.46D7,
     $     2.46375D7, 2.48D7, 2.5D7, 2.52D7, 2.52158D7, 2.54D7,
     $     2.54398D7, 2.56D7, 2.58D7, 2.6D7, 2.62D7, 2.63167D7, 2.64D7,
     $     2.6503D7, 2.66D7, 2.68D7, 2.68673D7, 2.7D7, 2.72D7, 2.74D7,
     $     2.76D7, 2.78D7, 2.8D7, 2.82D7, 2.84D7, 2.86D7, 2.87357D7,
     $     2.87708D7, 2.88D7, 2.89253D7, 2.9D7, 2.90099D7, 2.92D7,
     $     2.94D7, 2.96D7, 2.98D7, 3.0D7, 3.02D7, 3.04D7, 3.06D7, 3.08D7
     $     , 3.1D7, 3.12D7, 3.13473D7, 3.14D7, 3.14025D7, 3.16D7, 3.18D7
     $     , 3.2D7, 3.20307D7, 3.22D7, 3.24D7, 3.26D7, 3.27631D7, 3.28D7
     $     , 3.28223D7, 3.2956D7, 3.3D7, 3.30186D7, 3.32D7, 3.34D7,
     $     3.36D7, 3.38D7, 3.4D7, 3.42D7, 3.44D7, 3.46D7, 3.48D7, 3.5D7,
     $     3.52D7, 3.54D7, 3.56D7, 3.58D7, 3.6D7, 3.62D7, 3.64D7, 3.66D7
     $     , 3.68D7, 3.68757D7, 3.7D7, 3.72D7, 3.73196D7, 3.73775D7,
     $     3.74D7, 3.75082D7, 3.76D7, 3.76299D7, 3.78D7, 3.79889D7,
     $     3.8D7, 3.82D7, 3.84D7, 3.86D7, 3.88D7, 3.9D7, 3.91158D7,
     $     3.92D7, 3.94D7, 3.96D7, 3.96031D7, 3.98D7, 4.0D7, 4.07023D7,
     $     4.2D7, 4.22197D7, 4.24543D7, 4.25058D7, 4.26722D7, 4.33226D7,
     $     4.4D7, 4.42436D7, 4.44549D7, 4.46066D7, 4.6D7, 4.8D7,
     $     4.90239D7, 4.93603D7, 5.0D7, 5.05645D7, 5.10166D7, 5.13363D7,
     $     5.1821D7, 5.2D7, 5.27406D7, 5.29696D7, 5.4D7, 5.43748D7,
     $     5.4376D7, 5.6D7, 5.66668D7, 5.76406D7, 5.8D7, 5.95738D7,
     $     5.98109D7, 5.98838D7, 6.0D7, 6.00805D7, 6.0713D7, 6.19363D7,
     $     6.2D7, 6.33196D7, 6.4D7, 6.40041D7, 6.4965D7, 6.6D7,
     $     6.66672D7, 6.8D7, 6.83576D7, 6.85141D7, 6.91707D7, 7.0D7,
     $     7.14692D7, 7.17854D7, 7.2D7, 7.21671D7, 7.30793D7, 7.4D7,
     $     7.6D7, 7.63768D7, 7.8D7, 7.87361D7, 8.0D7, 8.00283D7,
     $     8.09908D7, 8.12276D7, 8.15613D7, 8.2D7, 8.4D7, 8.48079D7,
     $     8.57521D7, 8.6D7, 8.66789D7, 8.8D7, 8.84461D7, 8.91437D7,
     $     9.0D7, 9.2D7, 9.31994D7, 9.4D7, 9.43789D7, 9.51484D7, 9.6D7,
     $     9.65351D7, 9.67548D7, 9.8D7, 1.0D8, 1.02D8, 1.02911D8,
     $     1.03241D8, 1.04D8, 1.04063D8, 1.04344D8, 1.06D8, 1.0618D8,
     $     1.08D8, 1.1D8, 1.10512D8, 1.12D8, 1.1272D8, 1.12775D8, 1.14D8
     $     , 1.14584D8, 1.16D8, 1.18D8, 1.19719D8, 1.2D8, 1.20552D8,
     $     1.21745D8, 1.22D8, 1.24D8, 1.24346D8, 1.26D8, 1.28D8,
     $     1.28328D8, 1.29662D8, 1.29968D8, 1.3D8, 1.31727D8, 1.32D8,
     $     1.34D8, 1.36D8, 1.38D8, 1.38669D8, 1.4D8/

      data ( gdr(63,i), i=1, 369) /0.0, 0.0, 0.0, 3.30582D-5, 7.7172D-5,
     $     8.18595D-4, 9.68672D-4, 0.00140956, 0.00169457, 0.00171937,
     $     0.00216394, 0.00217329, 0.00220169, 0.00236935, 0.00244069,
     $     0.0025278, 0.00280768, 0.002839, 0.004659, 0.0083, 0.01607,
     $     0.0182358, 0.0225716, 0.0231595, 0.02345, 0.0236835, 0.02983,
     $     0.03586, 0.04244, 0.0453858, 0.05, 0.05842, 0.06764, 0.07781,
     $     0.08915, 0.0993763, 0.1019, 0.1162, 0.1323, 0.137095,
     $     0.139607, 0.1505, 0.1711, 0.1945, 0.2208, 0.2505, 0.2837,
     $     0.289711, 0.301618, 0.3207, 0.3612, 0.385682, 0.403798,
     $     0.4048, 0.4505, 0.4965, 0.5405, 0.546703, 0.554453, 0.5799,
     $     0.592169, 0.594831, 0.6115, 0.6323, 0.6412, 0.6373, 0.6234,
     $     0.6216, 0.5963, 0.579017, 0.567686, 0.5641, 0.52915, 0.5276,
     $     0.4894, 0.4513, 0.4145, 0.395314, 0.3799, 0.35052, 0.3477,
     $     0.342543, 0.3183, 0.317025, 0.307725, 0.2917, 0.2676,
     $     0.246505, 0.2459, 0.2265, 0.209, 0.1932, 0.179, 0.1662,
     $     0.1546, 0.147159, 0.1441, 0.143783, 0.1346, 0.126, 0.11876,
     $     0.1182, 0.111, 0.1044, 0.09835, 0.09281, 0.0877, 0.08297,
     $     0.0786, 0.07455, 0.07077, 0.0673398, 0.06725, 0.0644777,
     $     0.06399, 0.06092, 0.05806, 0.05537, 0.05285, 0.05046, 0.04823
     $     , 0.04612, 0.04412, 0.04293, 0.0422968, 0.04224, 0.0418388,
     $     0.04044, 0.03873, 0.0386507, 0.03712, 0.03558, 0.03412,
     $     0.0331981, 0.03273, 0.0314, 0.03011, 0.0298551, 0.0289,
     $     0.028844, 0.02774, 0.02662, 0.0257643, 0.0256263, 0.02552,
     $     0.0253173, 0.02446, 0.0235, 0.02269, 0.0226378, 0.02204,
     $     0.0219372, 0.02153, 0.02114, 0.02086, 0.02067, 0.0206056,
     $     0.02056, 0.0205393, 0.02052, 0.02051, 0.0205235, 0.02055,
     $     0.0206, 0.02065, 0.02068, 0.02069, 0.02069, 0.02066, 0.02062,
     $     0.02057, 0.020536, 0.0205273, 0.02052, 0.0204761, 0.02045,
     $     0.020447, 0.02039, 0.02031, 0.02024, 0.02018, 0.02012,
     $     0.02007, 0.02003, 0.01999, 0.01996, 0.01994, 0.01992,
     $     0.0199053, 0.0199, 0.0198999, 0.01989, 0.01988, 0.01986,
     $     0.0198585, 0.01985, 0.01984, 0.01983, 0.0198218, 0.01982,
     $     0.0198189, 0.0198122, 0.01981, 0.0198081, 0.01979, 0.01976,
     $     0.01974, 0.01972, 0.0197, 0.01968, 0.01965, 0.01962, 0.01959,
     $     0.01957, 0.01955, 0.01952, 0.01949, 0.01947, 0.01943, 0.01941
     $     , 0.01938, 0.01935, 0.01933, 0.0193186, 0.0193, 0.01927,
     $     0.019246, 0.0192345, 0.01923, 0.0192083, 0.01919, 0.0191855,
     $     0.01916, 0.0191317, 0.01913, 0.01909, 0.01906, 0.01901,
     $     0.01898, 0.01894, 0.0189226, 0.01891, 0.01887, 0.01884,
     $     0.0188394, 0.0188, 0.01876, 0.0186524, 0.01846, 0.0184295,
     $     0.0183972, 0.0183901, 0.0183674, 0.0182796, 0.01819,
     $     0.0181537, 0.0181225, 0.0181002, 0.0179, 0.01761, 0.0175115,
     $     0.0174797, 0.01742, 0.0173481, 0.0172912, 0.0172515,
     $     0.0171918, 0.01717, 0.0170608, 0.0170274, 0.01688, 0.0168283,
     $     0.0168281, 0.01661, 0.0165287, 0.0164123, 0.01637, 0.0162198,
     $     0.0161976, 0.0161908, 0.01618, 0.0161714, 0.016104, 0.0159765
     $     , 0.01597, 0.0158237, 0.01575, 0.0157496, 0.0156475, 0.01554,
     $     0.0154791, 0.01536, 0.0153291, 0.0153156, 0.0152597, 0.01519,
     $     0.0150571, 0.015029, 0.01501, 0.014993, 0.0149011, 0.01481,
     $     0.01462, 0.0145856, 0.01444, 0.0143842, 0.01429, 0.0142874,
     $     0.0142, 0.0141787, 0.0141489, 0.01411, 0.01398, 0.013927,
     $     0.0138659, 0.01385, 0.0137985, 0.0137, 0.0136684, 0.0136194,
     $     0.01356, 0.01345, 0.0133776, 0.01333, 0.0133051, 0.0132549,
     $     0.0132, 0.0131675, 0.0131543, 0.01308, 0.01297, 0.01286,
     $     0.0128049, 0.0127851, 0.01274, 0.0127375, 0.0127261, 0.01266,
     $     0.0126527, 0.01258, 0.01247, 0.0124467, 0.01238, 0.012351,
     $     0.0123488, 0.0123, 0.0122765, 0.01222, 0.01213, 0.0120525,
     $     0.01204, 0.0120233, 0.0119876, 0.01198, 0.01192, 0.0119095,
     $     0.01186, 0.0118, 0.0117901, 0.0117501, 0.0117409, 0.01174,
     $     0.0117141, 0.01171, 0.01167, 0.01161, 0.01156, 0.0115398,
     $     0.0115/


C     83209, JENDL evaluated data
      nbins(64) = 230
      data ( eg(64,i), i=1, 230) /3139900.0, 7400000.0, 7600000.0,
     $     7800000.0, 8000000.0, 8200000.0, 8400000.0, 8600000.0,
     $     8800000.0, 9000000.0, 9200000.0, 9400000.0, 9600000.0,
     $     9800000.0, 1.0D7, 1.02D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7,
     $     1.12D7, 1.14D7, 1.16D7, 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7
     $     , 1.28D7, 1.3D7, 1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7,
     $     1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7
     $     , 1.58D7, 1.6D7, 1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7,
     $     1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7
     $     , 1.88D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7,
     $     2.02D7, 2.04D7, 2.06D7, 2.08D7, 2.1D7, 2.12D7, 2.14D7, 2.16D7
     $     , 2.18D7, 2.2D7, 2.22D7, 2.24D7, 2.26D7, 2.28D7, 2.3D7,
     $     2.32D7, 2.34D7, 2.36D7, 2.38D7, 2.4D7, 2.42D7, 2.44D7, 2.46D7
     $     , 2.48D7, 2.5D7, 2.52D7, 2.54D7, 2.56D7, 2.58D7, 2.6D7,
     $     2.62D7, 2.64D7, 2.66D7, 2.68D7, 2.7D7, 2.72D7, 2.74D7, 2.76D7
     $     , 2.78D7, 2.8D7, 2.82D7, 2.84D7, 2.86D7, 2.88D7, 2.9D7,
     $     2.92D7, 2.94D7, 2.96D7, 2.98D7, 3.0D7, 3.05D7, 3.1D7, 3.15D7,
     $     3.2D7, 3.25D7, 3.3D7, 3.35D7, 3.4D7, 3.45D7, 3.5D7, 3.55D7,
     $     3.6D7, 3.65D7, 3.7D7, 3.75D7, 3.8D7, 3.85D7, 3.9D7, 3.95D7,
     $     4.0D7, 4.05D7, 4.1D7, 4.15D7, 4.2D7, 4.25D7, 4.3D7, 4.35D7,
     $     4.4D7, 4.45D7, 4.5D7, 4.55D7, 4.6D7, 4.65D7, 4.7D7, 4.75D7,
     $     4.8D7, 4.85D7, 4.9D7, 4.95D7, 5.0D7, 5.05D7, 5.1D7, 5.15D7,
     $     5.1923D7, 5.2D7, 5.25D7, 5.3D7, 5.35D7, 5.4D7, 5.45D7, 5.5D7,
     $     5.55D7, 5.6D7, 5.65D7, 5.7D7, 5.75D7, 5.8D7, 5.85D7, 5.9D7,
     $     5.95D7, 6.0D7, 6.05D7, 6.1D7, 6.15D7, 6.2D7, 6.25D7, 6.3D7,
     $     6.35D7, 6.4D7, 6.45D7, 6.5D7, 6.55D7, 6.6D7, 6.65D7, 6.7D7,
     $     6.75D7, 6.8D7, 6.85D7, 6.9D7, 6.95D7, 7.0D7, 7.05D7, 7.1D7,
     $     7.15D7, 7.2D7, 7.25D7, 7.3D7, 7.35D7, 7.4D7, 7.45D7, 7.5D7,
     $     7.55D7, 7.6D7, 7.65D7, 7.7D7, 7.75D7, 7.8D7, 7.85D7, 7.9D7,
     $     7.95D7, 8.0D7, 8.5D7, 9.0D7, 9.5D7, 1.0D8, 1.05D8, 1.1D8,
     $     1.15D8, 1.2D8, 1.25D8, 1.3D8, 1.35D8, 1.4D8, 1.45D8, 1.5D8/

      data ( gdr(64,i), i=1, 230) /0.0, 0.0, 0.0016314, 0.0063112,
     $     0.013878, 0.025097, 0.039573, 0.054334, 0.066918, 0.075807,
     $     0.078106, 0.079641, 0.081279, 0.08501, 0.092479, 0.10401,
     $     0.1199, 0.14034, 0.16531, 0.19449, 0.22704, 0.26165, 0.29448,
     $     0.32535, 0.35625, 0.38628, 0.4128, 0.44167, 0.47358, 0.50747,
     $     0.53976, 0.56471, 0.57557, 0.56962, 0.54434, 0.50745, 0.46408
     $     , 0.43711, 0.41124, 0.3865, 0.36352, 0.34242, 0.32303,
     $     0.30605, 0.29106, 0.27364, 0.25755, 0.24249, 0.22869, 0.21598
     $     , 0.20428, 0.19355, 0.18372, 0.17479, 0.16676, 0.15963,
     $     0.15338, 0.14798, 0.14347, 0.13938, 0.13618, 0.13384, 0.13233
     $     , 0.13161, 0.13161, 0.13221, 0.13329, 0.13466, 0.13615,
     $     0.1375, 0.13846, 0.13863, 0.13782, 0.13606, 0.13332, 0.12952,
     $     0.12492, 0.11965, 0.11391, 0.1079, 0.10181, 0.095789,
     $     0.089953, 0.084388, 0.079059, 0.074092, 0.069496, 0.065262,
     $     0.061377, 0.05782, 0.054546, 0.051552, 0.048813, 0.046308,
     $     0.044014, 0.041911, 0.03997, 0.038175, 0.036524, 0.035003,
     $     0.0336, 0.032303, 0.031104, 0.029993, 0.028961, 0.028004,
     $     0.027113, 0.026283, 0.025509, 0.024777, 0.024092, 0.023451,
     $     0.022853, 0.022289, 0.021759, 0.02066, 0.019732, 0.018938,
     $     0.018254, 0.017672, 0.017177, 0.016755, 0.016393, 0.016081,
     $     0.015816, 0.0156, 0.01542, 0.015269, 0.015142, 0.01504,
     $     0.01496, 0.014898, 0.014852, 0.014817, 0.014795, 0.014787,
     $     0.014788, 0.014797, 0.014812, 0.014833, 0.014862, 0.014895,
     $     0.014931, 0.014967, 0.015006, 0.015049, 0.015094, 0.015142,
     $     0.015191, 0.015243, 0.015296, 0.015349, 0.015404, 0.01546,
     $     0.015516, 0.015575, 0.015634, 0.015693, 0.01574291, 0.015752,
     $     0.015811, 0.015869, 0.015926, 0.01598, 0.016029, 0.016078,
     $     0.016125, 0.016171, 0.016217, 0.016261, 0.016304, 0.016346,
     $     0.016387, 0.016427, 0.016465, 0.016503, 0.016539, 0.016573,
     $     0.016607, 0.016639, 0.01667, 0.016699, 0.016727, 0.016753,
     $     0.016774, 0.016795, 0.016814, 0.016832, 0.016848, 0.016863,
     $     0.016877, 0.01689, 0.016901, 0.016911, 0.016919, 0.016926,
     $     0.016932, 0.016937, 0.01694, 0.016942, 0.016943, 0.016942,
     $     0.016941, 0.016937, 0.016931, 0.016924, 0.016916, 0.016906,
     $     0.016896, 0.016884, 0.016871, 0.016858, 0.016843, 0.016827,
     $     0.01681, 0.016792, 0.016551, 0.016221, 0.015841, 0.015403,
     $     0.014927, 0.014425, 0.013916, 0.013413, 0.012923, 0.01246,
     $     0.012041, 0.011673, 0.011366, 0.01113/


C     92235, JENDL evaluated data
      nbins(65) = 61
      data ( eg(65,i), i=1, 61) /5000000.0, 5298000.0, 5500000.0,
     $     6000000.0, 6500000.0, 7000000.0, 7500000.0, 8000000.0,
     $     8500000.0, 9000000.0, 9500000.0, 1.0D7, 1.05D7, 1.1D7, 1.15D7
     $     , 1.2D7, 1.214D7, 1.25D7, 1.3D7, 1.35D7, 1.4D7, 1.45D7, 1.5D7
     $     , 1.55D7, 1.6D7, 1.65D7, 1.7D7, 1.75D7, 1.8D7, 1.85D7, 1.9D7,
     $     1.95D7, 2.0D7, 2.25D7, 2.5D7, 2.75D7, 3.0D7, 3.25D7, 3.5D7,
     $     3.75D7, 4.0D7, 4.5D7, 5.0D7, 5.5D7, 6.0D7, 6.5D7, 7.0D7,
     $     7.5D7, 8.0D7, 8.5D7, 9.0D7, 9.5D7, 1.0D8, 1.05D8, 1.1D8,
     $     1.15D8, 1.2D8, 1.25D8, 1.3D8, 1.35D8, 1.4D8/

      data ( gdr(65,i), i=1, 61) /0.0, 0.0, 0.00534, 0.0133, 0.0207,
     $     0.03098, 0.04428, 0.05996, 0.08435, 0.1195, 0.1911, 0.2885,
     $     0.4097, 0.4786, 0.459, 0.4344, 0.4337, 0.4442, 0.4763, 0.5052
     $     , 0.5075, 0.4757, 0.4207, 0.3587, 0.3008, 0.2514, 0.2111,
     $     0.1787, 0.1527, 0.1319, 0.1149, 0.1011, 0.08962, 0.06064,
     $     0.04434, 0.03558, 0.03043, 0.02719, 0.02499, 0.02341, 0.02223
     $     , 0.02034, 0.0194, 0.01951, 0.01955, 0.01949, 0.01935,
     $     0.01923, 0.01924, 0.0192, 0.01896, 0.01869, 0.01862, 0.01867,
     $     0.01862, 0.01849, 0.01828, 0.01803, 0.01774, 0.01744, 0.01715
     $     /


C     92238, JENDL evaluated data
      nbins(66) = 61
      data ( eg(66,i), i=1, 61) /5000000.0, 5500000.0, 6000000.0,
     $     6153000.0, 6500000.0, 7000000.0, 7500000.0, 8000000.0,
     $     8500000.0, 9000000.0, 9500000.0, 1.0D7, 1.05D7, 1.1D7,
     $     1.128D7, 1.15D7, 1.2D7, 1.25D7, 1.3D7, 1.35D7, 1.4D7, 1.45D7,
     $     1.5D7, 1.55D7, 1.6D7, 1.65D7, 1.7D7, 1.75D7, 1.8D7, 1.85D7,
     $     1.9D7, 1.95D7, 2.0D7, 2.25D7, 2.5D7, 2.75D7, 3.0D7, 3.25D7,
     $     3.5D7, 3.75D7, 4.0D7, 4.5D7, 5.0D7, 5.5D7, 6.0D7, 6.5D7,
     $     7.0D7, 7.5D7, 8.0D7, 8.5D7, 9.0D7, 9.5D7, 1.0D8, 1.05D8,
     $     1.1D8, 1.15D8, 1.2D8, 1.25D8, 1.3D8, 1.35D8, 1.4D8/

      data ( gdr(66,i), i=1, 61) /0.0, 0.00133, 0.01012, 0.00779, 0.0267
     $     , 0.03946, 0.05328, 0.07218, 0.09923, 0.142, 0.2203, 0.3254,
     $     0.438, 0.4825, 0.4736, 0.462, 0.4504, 0.4667, 0.4945, 0.5115,
     $     0.5019, 0.4651, 0.4119, 0.3546, 0.3015, 0.2556, 0.2176,
     $     0.1864, 0.161, 0.1402, 0.1232, 0.109, 0.09718, 0.06635,
     $     0.04865, 0.03893, 0.03312, 0.02941, 0.02686, 0.02501, 0.02362
     $     , 0.02143, 0.02028, 0.02026, 0.0202, 0.02006, 0.01986,
     $     0.01969, 0.01966, 0.01959, 0.01932, 0.01903, 0.01893, 0.01897
     $     , 0.01891, 0.01876, 0.01855, 0.01828, 0.01799, 0.01768,
     $     0.01738/


C     93237, JENDL evaluated data
      nbins(67) = 147
      data ( eg(67,i), i=1, 147) /5300000.0, 5400000.0, 5600000.0,
     $     5800000.0, 6000000.0, 6200000.0, 6400000.0, 6600000.0,
     $     6800000.0, 7000000.0, 7200000.0, 7400000.0, 7600000.0,
     $     7800000.0, 8000000.0, 8200000.0, 8400000.0, 8600000.0,
     $     8800000.0, 9000000.0, 9200000.0, 9400000.0, 9600000.0,
     $     9800000.0, 1.0D7, 1.02D7, 1.04D7, 1.06D7, 1.08D7, 1.1D7,
     $     1.12D7, 1.14D7, 1.16D7, 1.18D7, 1.2D7, 1.22D7, 1.24D7, 1.26D7
     $     , 1.28D7, 1.3D7, 1.32D7, 1.34D7, 1.36D7, 1.38D7, 1.4D7,
     $     1.42D7, 1.44D7, 1.46D7, 1.48D7, 1.5D7, 1.52D7, 1.54D7, 1.56D7
     $     , 1.58D7, 1.6D7, 1.62D7, 1.64D7, 1.66D7, 1.68D7, 1.7D7,
     $     1.72D7, 1.74D7, 1.76D7, 1.78D7, 1.8D7, 1.82D7, 1.84D7, 1.86D7
     $     , 1.88D7, 1.9D7, 1.92D7, 1.94D7, 1.96D7, 1.98D7, 2.0D7, 2.1D7
     $     , 2.2D7, 2.3D7, 2.4D7, 2.5D7, 2.6D7, 2.7D7, 2.8D7, 2.9D7,
     $     3.0D7, 3.1D7, 3.2D7, 3.3D7, 3.4D7, 3.5D7, 3.6D7, 3.7D7, 3.8D7
     $     , 3.9D7, 4.0D7, 4.1D7, 4.2D7, 4.3D7, 4.4D7, 4.5D7, 4.6D7,
     $     4.7D7, 4.8D7, 4.9D7, 5.0D7, 5.1D7, 5.2D7, 5.3D7, 5.4D7, 5.5D7
     $     , 5.6D7, 5.7D7, 5.8D7, 5.9D7, 6.0D7, 6.1D7, 6.2D7, 6.3D7,
     $     6.4D7, 6.5D7, 6.6D7, 6.7D7, 6.8D7, 6.9D7, 7.0D7, 7.1D7, 7.2D7
     $     , 7.3D7, 7.4D7, 7.5D7, 7.6D7, 7.7D7, 7.8D7, 7.9D7, 8.0D7,
     $     8.5D7, 9.0D7, 9.5D7, 1.0D8, 1.05D8, 1.1D8, 1.15D8, 1.2D8,
     $     1.25D8, 1.3D8, 1.35D8, 1.4D8/

      data ( gdr(67,i), i=1, 147) /0.0, 0.0041785, 0.0054983, 0.0070802,
     $     0.0089256, 0.011024, 0.013356, 0.016073, 0.021485, 0.026668,
     $     0.032154, 0.037982, 0.044139, 0.05037, 0.057309, 0.065125,
     $     0.074041, 0.084159, 0.094565, 0.10677, 0.12137, 0.13915,
     $     0.16119, 0.18895, 0.27809, 0.33124, 0.39471, 0.46199, 0.51681
     $     , 0.54196, 0.53307, 0.50569, 0.4789, 0.46282, 0.45934,
     $     0.46639, 0.48528, 0.51454, 0.5396, 0.56164, 0.57931, 0.59215,
     $     0.60136, 0.60946, 0.61747, 0.62435, 0.62555, 0.61546, 0.58992
     $     , 0.54981, 0.50144, 0.45219, 0.40699, 0.36795, 0.33542,
     $     0.30892, 0.2874, 0.27003, 0.25639, 0.24581, 0.23757, 0.23098,
     $     0.22534, 0.22, 0.21437, 0.208, 0.20058, 0.19203, 0.18247,
     $     0.17216, 0.16138, 0.15394, 0.14794, 0.13995, 0.13194,
     $     0.097027, 0.072716, 0.056568, 0.045584, 0.03803, 0.032685,
     $     0.028876, 0.026094, 0.024055, 0.022537, 0.021435, 0.020644,
     $     0.020078, 0.019692, 0.019439, 0.019289, 0.019222, 0.019215,
     $     0.019257, 0.019335, 0.019441, 0.019568, 0.01971, 0.019862,
     $     0.02002, 0.020183, 0.020346, 0.020508, 0.020668, 0.020825,
     $     0.020976, 0.021122, 0.021261, 0.021394, 0.021519, 0.021636,
     $     0.021746, 0.021847, 0.021941, 0.022026, 0.022103, 0.022171,
     $     0.022232, 0.022284, 0.022329, 0.022365, 0.022394, 0.022415,
     $     0.022429, 0.022436, 0.022435, 0.022428, 0.022413, 0.022393,
     $     0.022366, 0.022332, 0.022293, 0.022248, 0.022197, 0.022141,
     $     0.021784, 0.021321, 0.020788, 0.020187, 0.01954, 0.018866,
     $     0.018183, 0.017508, 0.016858, 0.016248, 0.015692, 0.015203/

cABE change @2014/04/25
      getGDRxsec = 0.0d0
      if ( iZ .eq. 1 .and. iN .eq. 0 ) return

      ilo = 0
      ihi = 0
      middle = 0

      if ( Ene .gt. 150.0d0 ) return
cABE end

!     k; counter, l; first isotope number
      k = 0
      Enev = Ene* 1e6
      iA = iZ+ iN

!     Find nucleus
cABE change @2014/05/28, bisection method
      ilo = 1
      ihi = 67

      if (eval(ilo) .gt. iZ*1000+iA .or.
     &    eval(ihi) .lt. iZ*1000+iA) then
         k = 0
         goto 210
      endif

  200 middle = (ilo + ihi) / 2

      if (eval(middle) .eq. iZ*1000+iA)then
         k = middle
         goto 210
      elseif (eval(middle) .lt. iZ*1000+iA) then
         if (ilo .eq. middle) then
            k = 0
            goto 210
         endif
         ilo = middle
         goto 200
      else
         if (ihi .eq. middle) then
            k = 0
            goto 210
         endif
         ihi = middle
         goto 200
      endif

  210 continue
cABE end

!     Evaluation avaialbe
      if ( k /= 0) then
cABE mod @2014/05/28, i -> k
         if ( Enev < eg(k,1)) then
cABE end
           getGDRxsec = 0.
c            goto 100 ! 2014/4/6 Ogawa. Enable isomer generation by low E gamma

cABE add @2014/05/09
         elseif (Enev .ge. eg(k,nbins(k))) then
            getGDRxsec = ( gdr(k,nbins(k)) - gdr(k,nbins(k)-1))
     &                  / ( eg(k,nbins(k)) - eg(k,nbins(k)-1))
     &                  * ( Enev - eg(k,nbins(k))) + gdr(k,nbins(k))
            getGDRxsec = getGDRxsec - sigma_qd(iZ, iN, Ene)
            if (getGDRxsec .lt. 0.0d0 ) getGDRxsec = 0.0d0
cABE end

         else
cABE change @2014/05/28, bisection method
            ilo = 1
            ihi = nbins(k)

  300 middle = (ilo + ihi) / 2

            if (Enev .ge. eg(k,middle) .and.
     &          Enev .lt. eg(k,middle+1)) then
               j = middle
               goto 310
            elseif (eg(k,middle) .lt. Enev) then
               ilo = middle
               goto 300
            else
               ihi = middle
               goto 300
            endif

  310       continue

            getGDRxsec = (gdr(k,j+1) - gdr(k,j)) / (eg(k,j+1) - eg(k,j))
     $                  * (Enev - eg(k,j)) + gdr(k,j)
            getGDRxsec = getGDRxsec - sigma_qd(iZ, iN, Ene)

            if (getGDRxsec .lt. 0.0d0) then
               getGDRxsec = 0.0d0
            endif
cABE end

         end if

!     Evaluation not available, use equation
!     GDR: Lorentzian
      else if ( k == 0) then
100      center = -45.0062+ 81.0222* iA**( -0.0607777)
         width = -1401.18+ 1419.39* iA**( -0.00184425)


         p0 = 2.72795e-06
         p1 = 1.75037
         p2 = -0.501594
         p3 = -0.0105362
         height = p0*(( iZ*1000+ iZ+ iN)**p1)* exp( p2- p3*
     &    shellE(iZ,iN+iZ,2))

         getGDRxsec = height/( 1+(( Ene)**2- center**2)**2/ width**2
     $        /(Ene)**2)/ 1e3

      end if

      return
      end function


! ----------------------------------------------------------------------
!     Modify the decay width of Li-6 in GDR
      subroutine suppressalpha_li6( u)
!     u: incident photon energy (MeV)
! ----------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

      common /exiejn/ nimax
!$OMP THREADPRIVATE(/exiejn/)
      common /std1/ r(70),s(70),sigma,rr(70)
!$OMP THREADPRIVATE(/std1/)
      common /ejectl/ omega(70),ifa(70),ifz(70)

!     implicit none
      double precision u
      double precision sigma

      double precision exen(500)  ! excitation energy (MeV)
      double precision br(30,500) ! branching ration (-)
      integer lines               ! number of energy bins

      integer i, j, k           ! counters
      double precision integral

      double precision linearone

      integral = 0.

!     3006
! (1) 3005    1
! (2) 2005 1001
! (3) 3004    2
! (4) 2004 1002
! (5) 2003 1003
! (6) 1003 2003
! (7) 1002 2004

      lines = 1
      data( exen(k), k=1, 1) / 17.83900 /

!     Branching ratio of 3005
      data( br(1,k), k=1, 1) / 0.00000D+00 /

!     Branching ratio of 2005
      data( br(2,k), k=1, 1) / 0.00000D+00 /

!     Branching ratio of 3004
      data( br(3,k), k=1, 1) / 0.00000D+00 /

!     Branching ratio of 2004
      data( br(4,k), k=1, 1) / 7.66263D-01 /

!     Branching ratio of 2003
      data( br(5,k), k=1, 1) / 0.00000D+00 /

!     Branching ratio of 1003
      data( br(6,k), k=1, 1) / 0.00000D+00 /

!     Branching ratio of 1002
      data( br(7,k), k=1, 1) / 2.33727D-01 /

      do j = 1, lines
!     find corresponding excitation energy
         if ( u >= exen(1)-5. .and. u < exen(1)+5. ) then

            do i = 1, nimax
               if ( ifz(i)*1000+ ifa(i) == 1. ) then
                  if ( r(i) /= 0. ) then
                     integral = integral+ br(1,j)
                     r(i) = sigma* br(1,j)
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 1001. ) then
                  if ( r(i) /= 0. ) then
                     integral = integral+ br(2,j)
                     r(i) = sigma* br(2,j)
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 3004. ) then
                  if ( r(i) /= 0. ) then
                     integral = integral+ br(3,j)
                     r(i) = sigma* br(3,j)
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 1002. ) then
                  if ( r(i) /= 0. ) then
                     integral = integral+ br(4,j)
                     r(i) = sigma* br(4,j)
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 1003. ) then
                  if ( r(i) /= 0. ) then
                     integral = integral+ br(5,j)
                     r(i) = sigma* br(5,j)
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 2003. ) then
                  if ( r(i) /= 0. ) then
                     integral = integral+ br(6,j)
                     r(i) = sigma* br(6,j)
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 2004. ) then
                  if ( r(i) /= 0. ) then
                     integral = integral+ br(7,j)
                     r(i) = sigma* br(7,j)
                  end if

               end if

            end do

            do i = 1, nimax
               r(i) = r(i)/ integral
            end do

         end if

      end do

      end subroutine suppressalpha_li6



! ----------------------------------------------------------------------
!     Modify the decay width of B-10 in GDR
      subroutine suppressalpha_b10( u)
!     u: incident photon energy (MeV)
! ----------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

      common /exiejn/ nimax
!$OMP THREADPRIVATE(/exiejn/)
      common /std1/ r(70),s(70),sigma,rr(70)
!$OMP THREADPRIVATE(/std1/)
      common /ejectl/ omega(70),ifa(70),ifz(70)

!     implicit none
      double precision u
      double precision sigma

      double precision exen(500)  ! excitation energy (MeV)
      double precision br(30,500) ! branching ration (-)
      integer lines               ! number of energy bins

      integer i, j, k           ! counters
      double precision integral

      double precision linearone

      integral = 0.

!      5010
!  (1) 5009    1 o
!  (2) 4009 1001 o
!  (3) 5008    2 oo
!  (4) 4008 1002 o
!  (5) 3008 2002 oo
!  (6) 4007 1003 o
!  (7) 3007 2003 o
!  (8) 4006 1004 not in PHITS
!  (9) 3006 2004 o
! (10) 2006 3004 oo
! (11) 2005 3005 oo
! (12) 2004 3006 o
! (13) 2003 3007 o
! (14) 1003 4007 o
! (15) 1002 4008 not in PHITS

      lines = 245
      data( exen(k), k=1, 245) /
     $     10.78300, 14.07000, 15.26300, 16.19800,
     $     19.03500, 22.48400, 22.94400, 23.12000, 24.09800,
     $     24.30300, 24.73900, 25.35900, 26.00400, 26.19900,
     $     26.31700, 26.79400, 27.29200, 27.46100, 27.93600,
     $     28.24600, 28.61900, 28.97200, 29.05700, 29.60800,
     $     29.85000, 30.28000, 30.69500, 30.83500, 31.23100,
     $     31.51100, 31.81300, 32.00200, 32.15300, 32.57400,
     $     32.82000, 33.09800, 33.15600, 33.70300, 33.79100,
     $     34.33500, 34.44900, 34.74300, 34.94800, 35.14000,
     $     35.33400, 35.69700, 35.95000, 36.06500, 36.35400,
     $     36.83700, 37.12800, 37.33300, 37.52900, 37.58900,
     $     37.92100, 38.14600, 38.19100, 38.65300, 38.81800,
     $     38.88400, 39.13300, 39.45700, 39.60600, 40.00400,
     $     40.59000, 40.77000, 40.83300, 40.94900, 41.11400,
     $     41.25600, 41.53500, 41.67900, 42.05200, 42.15000,
     $     42.25700, 42.48700, 42.72800, 42.75300, 43.01600,
     $     43.15500, 43.38300, 43.52000, 43.86200, 44.01000,
     $     44.16200, 44.32000, 44.42600, 44.75100, 44.93000,
     $     45.02500, 45.12400, 45.34300, 45.65800, 45.91000,
     $     46.26000, 46.28700, 46.57800, 46.79000, 46.88600,
     $     47.07900, 47.25000, 47.34700, 47.48000, 47.69900,
     $     47.89500, 48.10400, 48.17900, 48.36700, 48.68200,
     $     48.86200, 48.98400, 49.17700, 49.21400, 49.36700,
     $     49.62100, 49.75300, 49.78800, 50.08100, 50.23600,
     $     50.42800, 50.56300, 50.86000, 50.92700, 51.11700,
     $     51.28700, 51.44100, 51.61400, 51.81700, 51.98100,
     $     52.19300, 52.21700, 52.35500, 52.53300, 52.66400,
     $     52.77500, 52.92000, 53.09900, 53.19300, 53.35800,
     $     53.53300, 53.75300, 53.95100, 54.17000, 54.37100,
     $     54.50900, 54.62000, 54.65800, 54.87400, 55.06800,
     $     55.16500, 55.40400, 55.45800, 55.65600, 55.73100,
     $     55.99100, 56.13000, 56.27100, 56.30300, 56.58800,
     $     56.66100, 56.72800, 56.87400, 57.15600, 57.24600,
     $     57.40100, 57.43600, 57.64400, 57.69300, 57.77600,
     $     57.98400, 58.25800, 58.44300, 58.50100, 58.51700,
     $     58.59600, 58.74400, 59.01600, 59.08800, 59.13700,
     $     59.36200, 59.50900, 59.65500, 59.74200, 59.91600,
     $     60.01700, 60.05500, 60.20000, 60.23500, 60.56200,
     $     60.64100, 60.67000, 60.94600, 60.99100, 61.07800,
     $     61.37600, 61.49400, 61.62500, 61.71500, 61.77500,
     $     61.90600, 61.99200, 62.10900, 62.41600, 62.51100,
     $     62.59300, 62.76900, 62.87600, 63.03100, 63.12300,
     $     63.28000, 63.47100, 63.56600, 63.67600, 63.76400,
     $     63.90100, 64.01500, 64.20900, 64.27800, 64.40400,
     $     64.54400, 64.65200, 64.74900, 64.85200, 64.96800,
     $     65.12700, 65.13900, 65.32400, 65.52300, 65.77900,
     $     65.86800, 65.92600, 65.99600, 66.09500, 66.32300,
     $     66.56400, 66.60500, 66.78200, 66.91800, 67.10700,
     $     67.30800, 67.37100, 67.42100, 67.45000, 67.57700,
     $     67.69300 /

!     Branching ratio of 5009
      data( br(1,k), k=1, 245) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00 /

!     Branching ratio of 4009
      data( br(2,k), k=1, 245) /
     $     3.26271D-01, 1.67156D-01, 1.37417D-01,
     $     1.22175D-01, 6.86080D-02, 4.28130D-02, 4.09000D-02,
     $     4.01150D-02, 3.78240D-02, 3.81206D-02, 3.81342D-02,
     $     3.71081D-02, 3.51953D-02, 3.45484D-02, 3.40397D-02,
     $     3.25173D-02, 3.07153D-02, 3.01236D-02, 2.85344D-02,
     $     2.75276D-02, 2.62378D-02, 2.49835D-02, 2.46796D-02,
     $     2.28323D-02, 2.20468D-02, 2.08086D-02, 1.95114D-02,
     $     1.91265D-02, 1.80279D-02, 1.72364D-02, 1.64622D-02,
     $     1.59790D-02, 1.56206D-02, 1.46304D-02, 1.40844D-02,
     $     1.34980D-02, 1.33858D-02, 1.24026D-02, 1.22604D-02,
     $     1.14174D-02, 1.12441D-02, 1.08328D-02, 1.05612D-02,
     $     1.03129D-02, 1.00638D-02, 9.60460D-03, 9.28626D-03,
     $     9.14314D-03, 8.77605D-03, 8.15097D-03, 7.79042D-03,
     $     7.55603D-03, 7.33822D-03, 7.27562D-03, 6.92750D-03,
     $     6.70866D-03, 6.66687D-03, 6.24321D-03, 6.10415D-03,
     $     6.04777D-03, 5.84158D-03, 5.59558D-03, 5.48690D-03,
     $     5.20895D-03, 4.82324D-03, 4.71004D-03, 4.67057D-03,
     $     4.60245D-03, 4.50603D-03, 4.42675D-03, 4.27455D-03,
     $     4.20055D-03, 4.01500D-03, 3.96870D-03, 3.91803D-03,
     $     3.81177D-03, 3.70327D-03, 3.69257D-03, 3.57990D-03,
     $     3.52246D-03, 3.42981D-03, 3.37725D-03, 3.24415D-03,
     $     3.18947D-03, 3.13427D-03, 3.07716D-03, 3.03971D-03,
     $     2.92737D-03, 2.86633D-03, 2.83558D-03, 2.80330D-03,
     $     2.73445D-03, 2.64191D-03, 2.57116D-03, 2.47865D-03,
     $     2.47162D-03, 2.39755D-03, 2.34680D-03, 2.32410D-03,
     $     2.27944D-03, 2.24099D-03, 2.21912D-03, 2.18957D-03,
     $     2.14400D-03, 2.10472D-03, 2.06233D-03, 2.04772D-03,
     $     2.01214D-03, 1.95503D-03, 1.92309D-03, 1.90222D-03,
     $     1.87010D-03, 1.86391D-03, 1.83855D-03, 1.79897D-03,
     $     1.77950D-03, 1.77441D-03, 1.73140D-03, 1.70883D-03,
     $     1.68120D-03, 1.66232D-03, 1.62275D-03, 1.61395D-03,
     $     1.58912D-03, 1.56724D-03, 1.54779D-03, 1.52586D-03,
     $     1.50057D-03, 1.48123D-03, 1.45696D-03, 1.45425D-03,
     $     1.43889D-03, 1.41953D-03, 1.40531D-03, 1.39361D-03,
     $     1.37852D-03, 1.36013D-03, 1.35081D-03, 1.33432D-03,
     $     1.31725D-03, 1.29690D-03, 1.27865D-03, 1.25914D-03,
     $     1.24126D-03, 1.22959D-03, 1.22006D-03, 1.21696D-03,
     $     1.19951D-03, 1.18405D-03, 1.17634D-03, 1.15786D-03,
     $     1.15378D-03, 1.13891D-03, 1.13326D-03, 1.11454D-03,
     $     1.10463D-03, 1.09460D-03, 1.09240D-03, 1.07304D-03,
     $     1.06811D-03, 1.06371D-03, 1.05429D-03, 1.03697D-03,
     $     1.03157D-03, 1.02218D-03, 1.02012D-03, 1.00805D-03,
     $     1.00525D-03, 1.00042D-03, 9.88863D-04, 9.74285D-04,
     $     9.64956D-04, 9.62097D-04, 9.61318D-04, 9.57288D-04,
     $     9.50154D-04, 9.37285D-04, 9.33865D-04, 9.31595D-04,
     $     9.21648D-04, 9.15350D-04, 9.09159D-04, 9.05500D-04,
     $     8.98213D-04, 8.94113D-04, 8.92569D-04, 8.86893D-04,
     $     8.85436D-04, 8.73214D-04, 8.70315D-04, 8.69271D-04,
     $     8.59154D-04, 8.57549D-04, 8.54399D-04, 8.43949D-04,
     $     8.39895D-04, 8.35396D-04, 8.32328D-04, 8.30367D-04,
     $     8.26135D-04, 8.23420D-04, 8.19723D-04, 8.10324D-04,
     $     8.07403D-04, 8.04926D-04, 7.99672D-04, 7.96629D-04,
     $     7.92130D-04, 7.89466D-04, 7.85065D-04, 7.80012D-04,
     $     7.77401D-04, 7.74598D-04, 7.72362D-04, 7.68951D-04,
     $     7.66108D-04, 7.61359D-04, 7.59699D-04, 7.56751D-04,
     $     7.53367D-04, 7.50896D-04, 7.48666D-04, 7.46336D-04,
     $     7.43775D-04, 7.40305D-04, 7.40053D-04, 7.36261D-04,
     $     7.32254D-04, 7.27188D-04, 7.25385D-04, 7.24230D-04,
     $     7.22873D-04, 7.21038D-04, 7.16998D-04, 7.12726D-04,
     $     7.11985D-04, 7.09005D-04, 7.06794D-04, 7.03696D-04,
     $     7.00533D-04, 6.99571D-04, 6.98814D-04, 6.98378D-04,
     $     6.96488D-04, 6.94792D-04 /

!     Branching ratio of 5008
      data( br(3,k), k=1, 245) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 4.53731D-08, 8.34383D-08, 2.62732D-07,
     $     1.19779D-06, 3.63154D-05, 6.51501D-05, 7.38644D-05,
     $     1.88941D-04, 3.17658D-04, 4.46324D-04, 8.12439D-04,
     $     8.71764D-04, 1.08607D-03, 1.25818D-03, 1.40306D-03,
     $     1.53833D-03, 1.61344D-03, 1.87423D-03, 2.02338D-03,
     $     2.18783D-03, 2.21446D-03, 2.45987D-03, 2.50073D-03,
     $     2.83901D-03, 2.90787D-03, 3.06114D-03, 3.16278D-03,
     $     3.25184D-03, 3.33754D-03, 3.50760D-03, 3.61989D-03,
     $     3.66407D-03, 3.78112D-03, 3.94613D-03, 4.02923D-03,
     $     4.07154D-03, 4.11598D-03, 4.12719D-03, 4.22411D-03,
     $     4.28575D-03, 4.29777D-03, 4.41277D-03, 4.43876D-03,
     $     4.44959D-03, 4.49492D-03, 4.54567D-03, 4.56703D-03,
     $     4.61716D-03, 4.68917D-03, 4.70898D-03, 4.71741D-03,
     $     4.72612D-03, 4.73794D-03, 4.74568D-03, 4.76356D-03,
     $     4.76797D-03, 4.78328D-03, 4.78671D-03, 4.79149D-03,
     $     4.80110D-03, 4.81317D-03, 4.81371D-03, 4.82348D-03,
     $     4.82737D-03, 4.83554D-03, 4.83624D-03, 4.83552D-03,
     $     4.83249D-03, 4.83061D-03, 4.82635D-03, 4.82271D-03,
     $     4.80880D-03, 4.79980D-03, 4.79403D-03, 4.78884D-03,
     $     4.77955D-03, 4.76599D-03, 4.75609D-03, 4.74224D-03,
     $     4.74132D-03, 4.72981D-03, 4.71886D-03, 4.71321D-03,
     $     4.70217D-03, 4.69121D-03, 4.68508D-03, 4.67586D-03,
     $     4.65993D-03, 4.64555D-03, 4.63204D-03, 4.62718D-03,
     $     4.61595D-03, 4.59559D-03, 4.58381D-03, 4.57512D-03,
     $     4.56131D-03, 4.55869D-03, 4.54810D-03, 4.53037D-03,
     $     4.52076D-03, 4.51817D-03, 4.49592D-03, 4.48380D-03,
     $     4.46925D-03, 4.45863D-03, 4.43450D-03, 4.42897D-03,
     $     4.41325D-03, 4.39831D-03, 4.38427D-03, 4.36802D-03,
     $     4.34860D-03, 4.33304D-03, 4.31289D-03, 4.31059D-03,
     $     4.29743D-03, 4.28034D-03, 4.26811D-03, 4.25767D-03,
     $     4.24388D-03, 4.22688D-03, 4.21783D-03, 4.20169D-03,
     $     4.18446D-03, 4.16330D-03, 4.14454D-03, 4.12430D-03,
     $     4.10530D-03, 4.09245D-03, 4.08209D-03, 4.07860D-03,
     $     4.05874D-03, 4.04067D-03, 4.03162D-03, 4.00914D-03,
     $     4.00402D-03, 3.98545D-03, 3.97830D-03, 3.95410D-03,
     $     3.94117D-03, 3.92779D-03, 3.92480D-03, 3.89809D-03,
     $     3.89119D-03, 3.88500D-03, 3.87182D-03, 3.84728D-03,
     $     3.83966D-03, 3.82642D-03, 3.82345D-03, 3.80609D-03,
     $     3.80206D-03, 3.79506D-03, 3.77821D-03, 3.75656D-03,
     $     3.74243D-03, 3.73807D-03, 3.73688D-03, 3.73077D-03,
     $     3.71984D-03, 3.70030D-03, 3.69519D-03, 3.69174D-03,
     $     3.67656D-03, 3.66685D-03, 3.65724D-03, 3.65152D-03,
     $     3.64005D-03, 3.63351D-03, 3.63105D-03, 3.62183D-03,
     $     3.61950D-03, 3.59931D-03, 3.59449D-03, 3.59273D-03,
     $     3.57558D-03, 3.57280D-03, 3.56738D-03, 3.54923D-03,
     $     3.54192D-03, 3.53374D-03, 3.52811D-03, 3.52446D-03,
     $     3.51658D-03, 3.51148D-03, 3.50456D-03, 3.48656D-03,
     $     3.48092D-03, 3.47611D-03, 3.46579D-03, 3.45971D-03,
     $     3.45077D-03, 3.44549D-03, 3.43676D-03, 3.42661D-03,
     $     3.42142D-03, 3.41576D-03, 3.41127D-03, 3.40440D-03,
     $     3.39866D-03, 3.38898D-03, 3.38553D-03, 3.37930D-03,
     $     3.37216D-03, 3.36683D-03, 3.36201D-03, 3.35695D-03,
     $     3.35136D-03, 3.34379D-03, 3.34324D-03, 3.33487D-03,
     $     3.32607D-03, 3.31493D-03, 3.31104D-03, 3.30855D-03,
     $     3.30559D-03, 3.30153D-03, 3.29250D-03, 3.28283D-03,
     $     3.28117D-03, 3.27434D-03, 3.26923D-03, 3.26210D-03,
     $     3.25470D-03, 3.25243D-03, 3.25063D-03, 3.24959D-03,
     $     3.24507D-03, 3.24098D-03 /

!     Branching ratio of 4008
      data( br(4,k), k=1, 245) /
     $     2.86503D-01, 3.71604D-01, 3.82547D-01,
     $     3.80221D-01, 3.05102D-01, 2.28484D-01, 2.22591D-01,
     $     2.20733D-01, 2.23897D-01, 2.23836D-01, 2.22384D-01,
     $     2.17919D-01, 2.11029D-01, 2.08254D-01, 2.06607D-01,
     $     1.99835D-01, 1.92075D-01, 1.89260D-01, 1.81315D-01,
     $     1.76262D-01, 1.70643D-01, 1.64467D-01, 1.62987D-01,
     $     1.54521D-01, 1.51101D-01, 1.45461D-01, 1.40049D-01,
     $     1.38340D-01, 1.33565D-01, 1.30175D-01, 1.26625D-01,
     $     1.24383D-01, 1.22662D-01, 1.18048D-01, 1.15477D-01,
     $     1.12665D-01, 1.12123D-01, 1.07314D-01, 1.06630D-01,
     $     1.02673D-01, 1.01864D-01, 9.99280D-02, 9.86420D-02,
     $     9.74426D-02, 9.62024D-02, 9.38595D-02, 9.21425D-02,
     $     9.13375D-02, 8.92089D-02, 8.53637D-02, 8.30544D-02,
     $     8.15121D-02, 8.00659D-02, 7.96505D-02, 7.73622D-02,
     $     7.58944D-02, 7.56132D-02, 7.27185D-02, 7.17427D-02,
     $     7.13450D-02, 6.98855D-02, 6.81271D-02, 6.73388D-02,
     $     6.52868D-02, 6.23416D-02, 6.14558D-02, 6.11440D-02,
     $     6.06044D-02, 5.98335D-02, 5.91976D-02, 5.79595D-02,
     $     5.73518D-02, 5.58170D-02, 5.54323D-02, 5.50072D-02,
     $     5.41074D-02, 5.31788D-02, 5.30867D-02, 5.21084D-02,
     $     5.16055D-02, 5.07849D-02, 5.03154D-02, 4.90907D-02,
     $     4.85804D-02, 4.80614D-02, 4.75169D-02, 4.71572D-02,
     $     4.60560D-02, 4.54471D-02, 4.51397D-02, 4.48153D-02,
     $     4.41229D-02, 4.31892D-02, 4.24675D-02, 4.15169D-02,
     $     4.14439D-02, 4.06685D-02, 4.01310D-02, 3.98886D-02,
     $     3.94070D-02, 3.89880D-02, 3.87471D-02, 3.84199D-02,
     $     3.79150D-02, 3.74776D-02, 3.69992D-02, 3.68347D-02,
     $     3.64327D-02, 3.57820D-02, 3.54132D-02, 3.51721D-02,
     $     3.47990D-02, 3.47266D-02, 3.44285D-02, 3.39626D-02,
     $     3.37333D-02, 3.36732D-02, 3.31570D-02, 3.28824D-02,
     $     3.25430D-02, 3.23095D-02, 3.18161D-02, 3.17049D-02,
     $     3.13884D-02, 3.11064D-02, 3.08530D-02, 3.05640D-02,
     $     3.02280D-02, 2.99709D-02, 2.96461D-02, 2.96097D-02,
     $     2.94025D-02, 2.91403D-02, 2.89462D-02, 2.87863D-02,
     $     2.85797D-02, 2.83255D-02, 2.81962D-02, 2.79649D-02,
     $     2.77243D-02, 2.74370D-02, 2.71772D-02, 2.68983D-02,
     $     2.66398D-02, 2.64710D-02, 2.63318D-02, 2.62868D-02,
     $     2.60319D-02, 2.58043D-02, 2.56898D-02, 2.54141D-02,
     $     2.53531D-02, 2.51284D-02, 2.50426D-02, 2.47577D-02,
     $     2.46059D-02, 2.44514D-02, 2.44176D-02, 2.41168D-02,
     $     2.40397D-02, 2.39709D-02, 2.38233D-02, 2.35518D-02,
     $     2.34668D-02, 2.33180D-02, 2.32855D-02, 2.30940D-02,
     $     2.30494D-02, 2.29721D-02, 2.27872D-02, 2.25533D-02,
     $     2.24033D-02, 2.23573D-02, 2.23447D-02, 2.22794D-02,
     $     2.21642D-02, 2.19552D-02, 2.18993D-02, 2.18622D-02,
     $     2.17001D-02, 2.15974D-02, 2.14959D-02, 2.14356D-02,
     $     2.13152D-02, 2.12474D-02, 2.12217D-02, 2.11277D-02,
     $     2.11032D-02, 2.08999D-02, 2.08513D-02, 2.08338D-02,
     $     2.06628D-02, 2.06355D-02, 2.05819D-02, 2.04030D-02,
     $     2.03330D-02, 2.02551D-02, 2.02018D-02, 2.01678D-02,
     $     2.00942D-02, 2.00470D-02, 1.99824D-02, 1.98173D-02,
     $     1.97656D-02, 1.97217D-02, 1.96280D-02, 1.95738D-02,
     $     1.94932D-02, 1.94453D-02, 1.93662D-02, 1.92754D-02,
     $     1.92282D-02, 1.91777D-02, 1.91372D-02, 1.90754D-02,
     $     1.90239D-02, 1.89372D-02, 1.89068D-02, 1.88527D-02,
     $     1.87903D-02, 1.87447D-02, 1.87032D-02, 1.86599D-02,
     $     1.86122D-02, 1.85474D-02, 1.85427D-02, 1.84721D-02,
     $     1.83972D-02, 1.83020D-02, 1.82680D-02, 1.82461D-02,
     $     1.82205D-02, 1.81858D-02, 1.81095D-02, 1.80282D-02,
     $     1.80139D-02, 1.79569D-02, 1.79145D-02, 1.78547D-02,
     $     1.77934D-02, 1.77748D-02, 1.77601D-02, 1.77516D-02,
     $     1.77147D-02, 1.76815D-02 /

!     Branching ratio of 3008
      data( br(5,k), k=1, 245) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 2.55233D-07, 1.98486D-06, 2.14550D-05,
     $     1.12353D-04, 2.93961D-04, 3.58321D-04, 3.95556D-04,
     $     5.71974D-04, 1.25176D-03, 1.52365D-03, 2.16441D-03,
     $     2.59036D-03, 3.12062D-03, 3.72076D-03, 3.88698D-03,
     $     4.94623D-03, 5.37820D-03, 6.12948D-03, 6.78828D-03,
     $     7.01331D-03, 7.60579D-03, 7.97325D-03, 8.33891D-03,
     $     8.52846D-03, 8.67545D-03, 9.02564D-03, 9.22227D-03,
     $     9.42818D-03, 9.47247D-03, 9.83804D-03, 9.89208D-03,
     $     1.01816D-02, 1.02331D-02, 1.03822D-02, 1.04837D-02,
     $     1.05716D-02, 1.06499D-02, 1.07663D-02, 1.08148D-02,
     $     1.08293D-02, 1.08425D-02, 1.07982D-02, 1.07556D-02,
     $     1.07295D-02, 1.06963D-02, 1.06890D-02, 1.06214D-02,
     $     1.05769D-02, 1.05687D-02, 1.04616D-02, 1.04313D-02,
     $     1.04155D-02, 1.03541D-02, 1.02889D-02, 1.02577D-02,
     $     1.01675D-02, 1.00034D-02, 9.94670D-03, 9.92556D-03,
     $     9.89451D-03, 9.84735D-03, 9.81058D-03, 9.73488D-03,
     $     9.70044D-03, 9.60763D-03, 9.58462D-03, 9.55714D-03,
     $     9.49710D-03, 9.43063D-03, 9.42454D-03, 9.35296D-03,
     $     9.31550D-03, 9.25127D-03, 9.21702D-03, 9.11333D-03,
     $     9.07058D-03, 9.02597D-03, 8.97667D-03, 8.94385D-03,
     $     8.83932D-03, 8.77744D-03, 8.74778D-03, 8.71467D-03,
     $     8.64396D-03, 8.54975D-03, 8.47489D-03, 8.37475D-03,
     $     8.36663D-03, 8.27940D-03, 8.22114D-03, 8.19414D-03,
     $     8.14008D-03, 8.09275D-03, 8.06422D-03, 8.02493D-03,
     $     7.96733D-03, 7.91719D-03, 7.85728D-03, 7.83686D-03,
     $     7.78692D-03, 7.70540D-03, 7.65785D-03, 7.62734D-03,
     $     7.58028D-03, 7.57076D-03, 7.53108D-03, 7.47136D-03,
     $     7.44298D-03, 7.43554D-03, 7.36799D-03, 7.33057D-03,
     $     7.28304D-03, 7.25016D-03, 7.18136D-03, 7.16547D-03,
     $     7.11926D-03, 7.07741D-03, 7.03953D-03, 6.99485D-03,
     $     6.94214D-03, 6.90275D-03, 6.85268D-03, 6.84703D-03,
     $     6.81487D-03, 6.77400D-03, 6.74311D-03, 6.71790D-03,
     $     6.68505D-03, 6.64417D-03, 6.62370D-03, 6.58598D-03,
     $     6.54646D-03, 6.50023D-03, 6.45702D-03, 6.41056D-03,
     $     6.36615D-03, 6.33780D-03, 6.31372D-03, 6.30619D-03,
     $     6.26319D-03, 6.22409D-03, 6.20410D-03, 6.15593D-03,
     $     6.14525D-03, 6.10556D-03, 6.09006D-03, 6.03952D-03,
     $     6.01213D-03, 5.98379D-03, 5.97771D-03, 5.92323D-03,
     $     5.90896D-03, 5.89641D-03, 5.86943D-03, 5.82011D-03,
     $     5.80457D-03, 5.77669D-03, 5.77065D-03, 5.73502D-03,
     $     5.72673D-03, 5.71197D-03, 5.67731D-03, 5.63359D-03,
     $     5.60600D-03, 5.59756D-03, 5.59527D-03, 5.58282D-03,
     $     5.56150D-03, 5.52220D-03, 5.51139D-03, 5.50430D-03,
     $     5.47387D-03, 5.45463D-03, 5.43551D-03, 5.42411D-03,
     $     5.40108D-03, 5.38822D-03, 5.38333D-03, 5.36569D-03,
     $     5.36084D-03, 5.32276D-03, 5.31357D-03, 5.31029D-03,
     $     5.27719D-03, 5.27193D-03, 5.26140D-03, 5.22642D-03,
     $     5.21270D-03, 5.19724D-03, 5.18661D-03, 5.17998D-03,
     $     5.16563D-03, 5.15646D-03, 5.14378D-03, 5.11144D-03,
     $     5.10107D-03, 5.09231D-03, 5.07357D-03, 5.06292D-03,
     $     5.04665D-03, 5.03687D-03, 5.02084D-03, 5.00283D-03,
     $     4.99307D-03, 4.98302D-03, 4.97490D-03, 4.96252D-03,
     $     4.95200D-03, 4.93436D-03, 4.92819D-03, 4.91736D-03,
     $     4.90445D-03, 4.89525D-03, 4.88680D-03, 4.87797D-03,
     $     4.86830D-03, 4.85507D-03, 4.85413D-03, 4.83999D-03,
     $     4.82486D-03, 4.80546D-03, 4.79828D-03, 4.79369D-03,
     $     4.78836D-03, 4.78129D-03, 4.76593D-03, 4.74920D-03,
     $     4.74621D-03, 4.73462D-03, 4.72607D-03, 4.71375D-03,
     $     4.70122D-03, 4.69743D-03, 4.69444D-03, 4.69271D-03,
     $     4.68522D-03, 4.67847D-03 /

!     Branching ratio of 4007
      data( br(6,k), k=1, 245) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 1.80728D-03, 3.40243D-02, 3.26454D-02,
     $     3.19925D-02, 2.75733D-02, 2.67841D-02, 2.52963D-02,
     $     2.33695D-02, 2.22462D-02, 2.23326D-02, 2.27900D-02,
     $     2.39350D-02, 2.71622D-02, 2.83137D-02, 3.16998D-02,
     $     3.39195D-02, 3.64318D-02, 3.98605D-02, 4.06824D-02,
     $     4.54246D-02, 4.71735D-02, 4.95057D-02, 5.21531D-02,
     $     5.29067D-02, 5.50870D-02, 5.66663D-02, 5.79079D-02,
     $     5.87140D-02, 5.92648D-02, 6.07584D-02, 6.15241D-02,
     $     6.23388D-02, 6.24777D-02, 6.36503D-02, 6.38116D-02,
     $     6.49165D-02, 6.51656D-02, 6.56815D-02, 6.59814D-02,
     $     6.62167D-02, 6.64255D-02, 6.67839D-02, 6.69419D-02,
     $     6.69795D-02, 6.70618D-02, 6.70566D-02, 6.69609D-02,
     $     6.68296D-02, 6.67128D-02, 6.66701D-02, 6.65678D-02,
     $     6.64714D-02, 6.64511D-02, 6.62033D-02, 6.60580D-02,
     $     6.60024D-02, 6.58302D-02, 6.56027D-02, 6.54956D-02,
     $     6.51844D-02, 6.47087D-02, 6.45625D-02, 6.45159D-02,
     $     6.44073D-02, 6.42521D-02, 6.41181D-02, 6.38756D-02,
     $     6.37387D-02, 6.34124D-02, 6.33304D-02, 6.32428D-02,
     $     6.30529D-02, 6.28801D-02, 6.28588D-02, 6.26508D-02,
     $     6.25354D-02, 6.23544D-02, 6.22307D-02, 6.18922D-02,
     $     6.17323D-02, 6.15760D-02, 6.13960D-02, 6.12751D-02,
     $     6.08723D-02, 6.06460D-02, 6.05235D-02, 6.03982D-02,
     $     6.01410D-02, 5.97938D-02, 5.95239D-02, 5.91668D-02,
     $     5.91395D-02, 5.88462D-02, 5.86302D-02, 5.85292D-02,
     $     5.83278D-02, 5.81431D-02, 5.80401D-02, 5.78917D-02,
     $     5.76530D-02, 5.74467D-02, 5.72337D-02, 5.71590D-02,
     $     5.69811D-02, 5.66741D-02, 5.64997D-02, 5.63823D-02,
     $     5.61964D-02, 5.61620D-02, 5.60172D-02, 5.57863D-02,
     $     5.56693D-02, 5.56383D-02, 5.53696D-02, 5.52240D-02,
     $     5.50494D-02, 5.49241D-02, 5.46479D-02, 5.45831D-02,
     $     5.43985D-02, 5.42274D-02, 5.40679D-02, 5.38824D-02,
     $     5.36664D-02, 5.34966D-02, 5.32787D-02, 5.32540D-02,
     $     5.31112D-02, 5.29296D-02, 5.27980D-02, 5.26877D-02,
     $     5.25456D-02, 5.23688D-02, 5.22764D-02, 5.21111D-02,
     $     5.19356D-02, 5.17241D-02, 5.15346D-02, 5.13311D-02,
     $     5.11410D-02, 5.10156D-02, 5.09131D-02, 5.08793D-02,
     $     5.06868D-02, 5.05133D-02, 5.04248D-02, 5.02085D-02,
     $     5.01596D-02, 4.99793D-02, 4.99094D-02, 4.96780D-02,
     $     4.95540D-02, 4.94263D-02, 4.93982D-02, 4.91459D-02,
     $     4.90802D-02, 4.90220D-02, 4.88973D-02, 4.86679D-02,
     $     4.85965D-02, 4.84719D-02, 4.84443D-02, 4.82820D-02,
     $     4.82445D-02, 4.81792D-02, 4.80229D-02, 4.78250D-02,
     $     4.76976D-02, 4.76585D-02, 4.76478D-02, 4.75924D-02,
     $     4.74944D-02, 4.73183D-02, 4.72713D-02, 4.72403D-02,
     $     4.71055D-02, 4.70199D-02, 4.69350D-02, 4.68848D-02,
     $     4.67832D-02, 4.67256D-02, 4.67039D-02, 4.66237D-02,
     $     4.66028D-02, 4.64278D-02, 4.63856D-02, 4.63703D-02,
     $     4.62186D-02, 4.61943D-02, 4.61461D-02, 4.59845D-02,
     $     4.59205D-02, 4.58477D-02, 4.57976D-02, 4.57656D-02,
     $     4.56961D-02, 4.56513D-02, 4.55900D-02, 4.54320D-02,
     $     4.53816D-02, 4.53386D-02, 4.52469D-02, 4.51936D-02,
     $     4.51143D-02, 4.50672D-02, 4.49894D-02, 4.49003D-02,
     $     4.48540D-02, 4.48046D-02, 4.47650D-02, 4.47044D-02,
     $     4.46533D-02, 4.45673D-02, 4.45369D-02, 4.44824D-02,
     $     4.44184D-02, 4.43714D-02, 4.43286D-02, 4.42836D-02,
     $     4.42340D-02, 4.41666D-02, 4.41617D-02, 4.40885D-02,
     $     4.40110D-02, 4.39133D-02, 4.38783D-02, 4.38559D-02,
     $     4.38296D-02, 4.37942D-02, 4.37157D-02, 4.36308D-02,
     $     4.36159D-02, 4.35559D-02, 4.35112D-02, 4.34477D-02,
     $     4.33821D-02, 4.33620D-02, 4.33461D-02, 4.33369D-02,
     $     4.32969D-02, 4.32606D-02 /

!     Branching ratio of 3007
      data( br(7,k), k=1, 245) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 1.79269D-02, 5.41765D-02, 5.28369D-02,
     $     5.20065D-02, 4.54577D-02, 4.42607D-02, 4.20620D-02,
     $     3.94893D-02, 4.70115D-02, 5.24728D-02, 5.56186D-02,
     $     6.89266D-02, 8.16239D-02, 8.55298D-02, 9.61011D-02,
     $     1.02669D-01, 1.09915D-01, 1.18250D-01, 1.20141D-01,
     $     1.29926D-01, 1.33181D-01, 1.38519D-01, 1.43544D-01,
     $     1.44943D-01, 1.48365D-01, 1.50419D-01, 1.52190D-01,
     $     1.53353D-01, 1.53953D-01, 1.55595D-01, 1.56679D-01,
     $     1.57575D-01, 1.57741D-01, 1.59077D-01, 1.59280D-01,
     $     1.60964D-01, 1.61308D-01, 1.62178D-01, 1.62754D-01,
     $     1.63231D-01, 1.63650D-01, 1.64436D-01, 1.64701D-01,
     $     1.64719D-01, 1.64707D-01, 1.64316D-01, 1.64074D-01,
     $     1.63747D-01, 1.63403D-01, 1.63303D-01, 1.63098D-01,
     $     1.62853D-01, 1.62809D-01, 1.62173D-01, 1.61869D-01,
     $     1.61738D-01, 1.61370D-01, 1.60879D-01, 1.60642D-01,
     $     1.59948D-01, 1.58792D-01, 1.58427D-01, 1.58306D-01,
     $     1.58050D-01, 1.57660D-01, 1.57335D-01, 1.56736D-01,
     $     1.56418D-01, 1.55639D-01, 1.55446D-01, 1.55225D-01,
     $     1.54736D-01, 1.54298D-01, 1.54246D-01, 1.53708D-01,
     $     1.53403D-01, 1.52912D-01, 1.52602D-01, 1.51721D-01,
     $     1.51312D-01, 1.50913D-01, 1.50440D-01, 1.50131D-01,
     $     1.49088D-01, 1.48496D-01, 1.48186D-01, 1.47855D-01,
     $     1.47174D-01, 1.46263D-01, 1.45536D-01, 1.44568D-01,
     $     1.44491D-01, 1.43669D-01, 1.43091D-01, 1.42819D-01,
     $     1.42278D-01, 1.41785D-01, 1.41506D-01, 1.41102D-01,
     $     1.40474D-01, 1.39930D-01, 1.39338D-01, 1.39130D-01,
     $     1.38637D-01, 1.37790D-01, 1.37305D-01, 1.36984D-01,
     $     1.36477D-01, 1.36382D-01, 1.35977D-01, 1.35344D-01,
     $     1.35030D-01, 1.34947D-01, 1.34210D-01, 1.33804D-01,
     $     1.33311D-01, 1.32959D-01, 1.32198D-01, 1.32018D-01,
     $     1.31507D-01, 1.31035D-01, 1.30598D-01, 1.30085D-01,
     $     1.29486D-01, 1.29023D-01, 1.28431D-01, 1.28364D-01,
     $     1.27977D-01, 1.27485D-01, 1.27124D-01, 1.26823D-01,
     $     1.26435D-01, 1.25950D-01, 1.25700D-01, 1.25248D-01,
     $     1.24767D-01, 1.24196D-01, 1.23676D-01, 1.23118D-01,
     $     1.22590D-01, 1.22248D-01, 1.21964D-01, 1.21872D-01,
     $     1.21348D-01, 1.20874D-01, 1.20630D-01, 1.20037D-01,
     $     1.19903D-01, 1.19411D-01, 1.19218D-01, 1.18586D-01,
     $     1.18246D-01, 1.17894D-01, 1.17817D-01, 1.17131D-01,
     $     1.16951D-01, 1.16792D-01, 1.16451D-01, 1.15825D-01,
     $     1.15630D-01, 1.15285D-01, 1.15209D-01, 1.14761D-01,
     $     1.14658D-01, 1.14475D-01, 1.14042D-01, 1.13495D-01,
     $     1.13147D-01, 1.13040D-01, 1.13011D-01, 1.12856D-01,
     $     1.12587D-01, 1.12097D-01, 1.11964D-01, 1.11877D-01,
     $     1.11502D-01, 1.11263D-01, 1.11026D-01, 1.10886D-01,
     $     1.10601D-01, 1.10440D-01, 1.10380D-01, 1.10158D-01,
     $     1.10099D-01, 1.09617D-01, 1.09501D-01, 1.09459D-01,
     $     1.09039D-01, 1.08972D-01, 1.08838D-01, 1.08390D-01,
     $     1.08214D-01, 1.08012D-01, 1.07874D-01, 1.07786D-01,
     $     1.07596D-01, 1.07474D-01, 1.07306D-01, 1.06877D-01,
     $     1.06738D-01, 1.06621D-01, 1.06371D-01, 1.06228D-01,
     $     1.06011D-01, 1.05881D-01, 1.05668D-01, 1.05425D-01,
     $     1.05296D-01, 1.05162D-01, 1.05053D-01, 1.04888D-01,
     $     1.04746D-01, 1.04510D-01, 1.04427D-01, 1.04279D-01,
     $     1.04103D-01, 1.03976D-01, 1.03860D-01, 1.03738D-01,
     $     1.03604D-01, 1.03421D-01, 1.03408D-01, 1.03211D-01,
     $     1.03001D-01, 1.02734D-01, 1.02637D-01, 1.02575D-01,
     $     1.02502D-01, 1.02405D-01, 1.02192D-01, 1.01960D-01,
     $     1.01919D-01, 1.01756D-01, 1.01635D-01, 1.01462D-01,
     $     1.01285D-01, 1.01230D-01, 1.01188D-01, 1.01163D-01,
     $     1.01055D-01, 1.00958D-01 /

!     Branching ratio of 4006
      data( br(8,k), k=1, 245) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     6.36910D-20, 5.83347D-07, 5.11461D-06, 2.02504D-04,
     $     2.37302D-04, 3.75917D-04, 5.32667D-04, 6.72441D-04,
     $     8.25383D-04, 9.58261D-04, 1.58991D-03, 1.92501D-03,
     $     2.21386D-03, 2.26231D-03, 2.62619D-03, 2.67197D-03,
     $     2.93784D-03, 2.98755D-03, 3.07426D-03, 3.11901D-03,
     $     3.15223D-03, 3.17859D-03, 3.25295D-03, 3.38305D-03,
     $     3.45028D-03, 3.63580D-03, 3.99347D-03, 4.21343D-03,
     $     4.35603D-03, 4.53833D-03, 4.59475D-03, 4.97350D-03,
     $     5.18342D-03, 5.22190D-03, 5.68812D-03, 5.85503D-03,
     $     5.93011D-03, 6.22236D-03, 6.52817D-03, 6.64948D-03,
     $     6.93825D-03, 7.45055D-03, 7.63452D-03, 7.70182D-03,
     $     7.80407D-03, 7.95201D-03, 8.06859D-03, 8.28857D-03,
     $     8.37389D-03, 8.57816D-03, 8.62669D-03, 8.68566D-03,
     $     8.81729D-03, 8.97424D-03, 8.98647D-03, 9.13798D-03,
     $     9.21850D-03, 9.36764D-03, 9.43940D-03, 9.63988D-03,
     $     9.71253D-03, 9.79681D-03, 9.88915D-03, 9.95542D-03,
     $     1.01535D-02, 1.02738D-02, 1.03264D-02, 1.03852D-02,
     $     1.05182D-02, 1.07180D-02, 1.08791D-02, 1.10796D-02,
     $     1.10955D-02, 1.12661D-02, 1.13738D-02, 1.14204D-02,
     $     1.15153D-02, 1.15937D-02, 1.16474D-02, 1.17155D-02,
     $     1.18097D-02, 1.18925D-02, 1.20040D-02, 1.20389D-02,
     $     1.21298D-02, 1.22583D-02, 1.23357D-02, 1.23811D-02,
     $     1.24473D-02, 1.24636D-02, 1.25292D-02, 1.26202D-02,
     $     1.26574D-02, 1.26669D-02, 1.27601D-02, 1.28107D-02,
     $     1.28867D-02, 1.29340D-02, 1.30258D-02, 1.30468D-02,
     $     1.31182D-02, 1.31835D-02, 1.32405D-02, 1.33112D-02,
     $     1.34005D-02, 1.34596D-02, 1.35366D-02, 1.35456D-02,
     $     1.35933D-02, 1.36545D-02, 1.37079D-02, 1.37461D-02,
     $     1.37976D-02, 1.38620D-02, 1.38902D-02, 1.39497D-02,
     $     1.40073D-02, 1.40681D-02, 1.41349D-02, 1.42074D-02,
     $     1.42818D-02, 1.43257D-02, 1.43690D-02, 1.43800D-02,
     $     1.44435D-02, 1.45050D-02, 1.45362D-02, 1.46089D-02,
     $     1.46239D-02, 1.46877D-02, 1.47124D-02, 1.47934D-02,
     $     1.48389D-02, 1.48858D-02, 1.48954D-02, 1.49850D-02,
     $     1.50080D-02, 1.50285D-02, 1.50733D-02, 1.51516D-02,
     $     1.51788D-02, 1.52327D-02, 1.52428D-02, 1.53043D-02,
     $     1.53201D-02, 1.53493D-02, 1.54145D-02, 1.54962D-02,
     $     1.55445D-02, 1.55593D-02, 1.55632D-02, 1.55886D-02,
     $     1.56262D-02, 1.57051D-02, 1.57291D-02, 1.57453D-02,
     $     1.58109D-02, 1.58506D-02, 1.58902D-02, 1.59161D-02,
     $     1.59664D-02, 1.59924D-02, 1.60034D-02, 1.60389D-02,
     $     1.60510D-02, 1.61253D-02, 1.61448D-02, 1.61513D-02,
     $     1.62239D-02, 1.62355D-02, 1.62582D-02, 1.63379D-02,
     $     1.63698D-02, 1.64028D-02, 1.64264D-02, 1.64397D-02,
     $     1.64674D-02, 1.64846D-02, 1.65091D-02, 1.65758D-02,
     $     1.65973D-02, 1.66165D-02, 1.66601D-02, 1.66816D-02,
     $     1.67207D-02, 1.67449D-02, 1.67832D-02, 1.68220D-02,
     $     1.68487D-02, 1.68718D-02, 1.68924D-02, 1.69231D-02,
     $     1.69486D-02, 1.69939D-02, 1.70098D-02, 1.70349D-02,
     $     1.70670D-02, 1.70874D-02, 1.71083D-02, 1.71301D-02,
     $     1.71518D-02, 1.71840D-02, 1.71861D-02, 1.72167D-02,
     $     1.72510D-02, 1.73021D-02, 1.73236D-02, 1.73379D-02,
     $     1.73541D-02, 1.73734D-02, 1.74116D-02, 1.74565D-02,
     $     1.74661D-02, 1.74959D-02, 1.75169D-02, 1.75499D-02,
     $     1.75827D-02, 1.75921D-02, 1.75995D-02, 1.76038D-02,
     $     1.76221D-02, 1.76385D-02 /

!     Branching ratio of 3006
      data( br(9,k), k=1, 245) /
     $     3.84365D-01, 1.26957D-01, 9.79471D-02,
     $     8.46339D-02, 4.63329D-02, 2.37761D-02, 2.23374D-02,
     $     2.17632D-02, 1.83934D-02, 1.78313D-02, 1.67907D-02,
     $     1.55009D-02, 1.46250D-02, 1.44039D-02, 1.45471D-02,
     $     1.42709D-02, 1.54498D-02, 1.60157D-02, 1.83669D-02,
     $     2.06375D-02, 2.41220D-02, 2.74333D-02, 2.83990D-02,
     $     3.52157D-02, 3.85387D-02, 4.34921D-02, 4.90706D-02,
     $     5.06969D-02, 5.56125D-02, 5.94736D-02, 6.35274D-02,
     $     6.60284D-02, 6.79876D-02, 7.35437D-02, 7.63844D-02,
     $     7.93264D-02, 7.98791D-02, 8.46374D-02, 8.52857D-02,
     $     8.85410D-02, 8.91307D-02, 9.03867D-02, 9.11517D-02,
     $     9.18533D-02, 9.25725D-02, 9.39171D-02, 9.50959D-02,
     $     9.56818D-02, 9.72564D-02, 1.00237D-01, 1.02014D-01,
     $     1.03268D-01, 1.04530D-01, 1.04886D-01, 1.06651D-01,
     $     1.07752D-01, 1.07944D-01, 1.10082D-01, 1.10794D-01,
     $     1.11103D-01, 1.12143D-01, 1.13407D-01, 1.13975D-01,
     $     1.15436D-01, 1.17495D-01, 1.18098D-01, 1.18307D-01,
     $     1.18679D-01, 1.19276D-01, 1.19759D-01, 1.20667D-01,
     $     1.21102D-01, 1.22253D-01, 1.22539D-01, 1.22880D-01,
     $     1.23616D-01, 1.24237D-01, 1.24303D-01, 1.25020D-01,
     $     1.25424D-01, 1.26059D-01, 1.26388D-01, 1.27227D-01,
     $     1.27585D-01, 1.27903D-01, 1.28309D-01, 1.28537D-01,
     $     1.29298D-01, 1.29706D-01, 1.29908D-01, 1.30159D-01,
     $     1.30665D-01, 1.31283D-01, 1.31824D-01, 1.32573D-01,
     $     1.32640D-01, 1.33327D-01, 1.33714D-01, 1.33906D-01,
     $     1.34261D-01, 1.34582D-01, 1.34749D-01, 1.35025D-01,
     $     1.35415D-01, 1.35750D-01, 1.36160D-01, 1.36316D-01,
     $     1.36650D-01, 1.37255D-01, 1.37587D-01, 1.37805D-01,
     $     1.38164D-01, 1.38225D-01, 1.38511D-01, 1.38929D-01,
     $     1.39129D-01, 1.39181D-01, 1.39663D-01, 1.39955D-01,
     $     1.40266D-01, 1.40496D-01, 1.40927D-01, 1.41031D-01,
     $     1.41277D-01, 1.41484D-01, 1.41669D-01, 1.41896D-01,
     $     1.42134D-01, 1.42315D-01, 1.42527D-01, 1.42550D-01,
     $     1.42695D-01, 1.42872D-01, 1.42992D-01, 1.43110D-01,
     $     1.43260D-01, 1.43444D-01, 1.43534D-01, 1.43686D-01,
     $     1.43864D-01, 1.44064D-01, 1.44258D-01, 1.44467D-01,
     $     1.44675D-01, 1.44792D-01, 1.44889D-01, 1.44921D-01,
     $     1.45105D-01, 1.45262D-01, 1.45352D-01, 1.45557D-01,
     $     1.45604D-01, 1.45744D-01, 1.45813D-01, 1.45987D-01,
     $     1.46092D-01, 1.46216D-01, 1.46237D-01, 1.46409D-01,
     $     1.46466D-01, 1.46507D-01, 1.46597D-01, 1.46769D-01,
     $     1.46818D-01, 1.46914D-01, 1.46939D-01, 1.47082D-01,
     $     1.47107D-01, 1.47165D-01, 1.47286D-01, 1.47437D-01,
     $     1.47526D-01, 1.47553D-01, 1.47560D-01, 1.47606D-01,
     $     1.47684D-01, 1.47826D-01, 1.47873D-01, 1.47897D-01,
     $     1.47996D-01, 1.48066D-01, 1.48141D-01, 1.48175D-01,
     $     1.48264D-01, 1.48316D-01, 1.48331D-01, 1.48388D-01,
     $     1.48406D-01, 1.48534D-01, 1.48559D-01, 1.48569D-01,
     $     1.48669D-01, 1.48682D-01, 1.48720D-01, 1.48807D-01,
     $     1.48828D-01, 1.48870D-01, 1.48893D-01, 1.48906D-01,
     $     1.48939D-01, 1.48960D-01, 1.48995D-01, 1.49045D-01,
     $     1.49071D-01, 1.49084D-01, 1.49100D-01, 1.49112D-01,
     $     1.49130D-01, 1.49145D-01, 1.49172D-01, 1.49203D-01,
     $     1.49218D-01, 1.49233D-01, 1.49242D-01, 1.49257D-01,
     $     1.49283D-01, 1.49309D-01, 1.49315D-01, 1.49325D-01,
     $     1.49347D-01, 1.49355D-01, 1.49358D-01, 1.49360D-01,
     $     1.49369D-01, 1.49376D-01, 1.49377D-01, 1.49386D-01,
     $     1.49405D-01, 1.49418D-01, 1.49430D-01, 1.49433D-01,
     $     1.49436D-01, 1.49441D-01, 1.49452D-01, 1.49462D-01,
     $     1.49461D-01, 1.49463D-01, 1.49464D-01, 1.49471D-01,
     $     1.49468D-01, 1.49467D-01, 1.49466D-01, 1.49465D-01,
     $     1.49460D-01, 1.49456D-01 /

!     Branching ratio of 2006
      data( br(10,k), k=1, 245) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     4.97527D-12, 1.15909D-06, 3.23505D-05, 4.68857D-05,
     $     1.61984D-04, 2.51673D-04, 4.77614D-04, 7.09592D-04,
     $     7.98984D-04, 1.12285D-03, 1.46533D-03, 1.99762D-03,
     $     2.34105D-03, 2.59764D-03, 3.20378D-03, 3.50192D-03,
     $     3.78374D-03, 3.83634D-03, 4.20089D-03, 4.24191D-03,
     $     4.41891D-03, 4.44357D-03, 4.50574D-03, 4.55189D-03,
     $     4.61633D-03, 4.70700D-03, 4.90141D-03, 5.02640D-03,
     $     5.07862D-03, 5.19449D-03, 5.38980D-03, 5.55818D-03,
     $     5.67522D-03, 5.76665D-03, 5.79294D-03, 5.90781D-03,
     $     5.98061D-03, 5.99531D-03, 6.12682D-03, 6.17195D-03,
     $     6.18657D-03, 6.23497D-03, 6.29686D-03, 6.32588D-03,
     $     6.41614D-03, 6.53165D-03, 6.55469D-03, 6.56060D-03,
     $     6.57535D-03, 6.59226D-03, 6.60711D-03, 6.63144D-03,
     $     6.64667D-03, 6.67788D-03, 6.68599D-03, 6.69312D-03,
     $     6.70835D-03, 6.72539D-03, 6.72793D-03, 6.74922D-03,
     $     6.75923D-03, 6.77316D-03, 6.78653D-03, 6.81721D-03,
     $     6.83290D-03, 6.84895D-03, 6.86402D-03, 6.87454D-03,
     $     6.90723D-03, 6.92562D-03, 6.93802D-03, 6.94887D-03,
     $     6.97189D-03, 7.00357D-03, 7.02610D-03, 7.05758D-03,
     $     7.05935D-03, 7.08033D-03, 7.10373D-03, 7.11451D-03,
     $     7.13655D-03, 7.15664D-03, 7.16691D-03, 7.17998D-03,
     $     7.20490D-03, 7.22726D-03, 7.24455D-03, 7.25059D-03,
     $     7.26556D-03, 7.29124D-03, 7.30580D-03, 7.31720D-03,
     $     7.33566D-03, 7.33885D-03, 7.35097D-03, 7.37612D-03,
     $     7.39174D-03, 7.39602D-03, 7.42900D-03, 7.44521D-03,
     $     7.46435D-03, 7.47780D-03, 7.51009D-03, 7.51666D-03,
     $     7.53492D-03, 7.55140D-03, 7.56657D-03, 7.58138D-03,
     $     7.59809D-03, 7.61310D-03, 7.63172D-03, 7.63375D-03,
     $     7.64523D-03, 7.66010D-03, 7.67021D-03, 7.67918D-03,
     $     7.69099D-03, 7.70577D-03, 7.71467D-03, 7.72926D-03,
     $     7.74482D-03, 7.76648D-03, 7.78300D-03, 7.80110D-03,
     $     7.81576D-03, 7.82772D-03, 7.83605D-03, 7.83946D-03,
     $     7.85835D-03, 7.87476D-03, 7.88226D-03, 7.90123D-03,
     $     7.90550D-03, 7.92070D-03, 7.92568D-03, 7.94569D-03,
     $     7.95539D-03, 7.96440D-03, 7.96684D-03, 7.98848D-03,
     $     7.99338D-03, 7.99836D-03, 8.00894D-03, 8.02972D-03,
     $     8.03604D-03, 8.04540D-03, 8.04761D-03, 8.06087D-03,
     $     8.06415D-03, 8.06883D-03, 8.08220D-03, 8.10023D-03,
     $     8.11339D-03, 8.11753D-03, 8.11869D-03, 8.12321D-03,
     $     8.13305D-03, 8.14906D-03, 8.15230D-03, 8.15480D-03,
     $     8.16785D-03, 8.17644D-03, 8.18475D-03, 8.18982D-03,
     $     8.19921D-03, 8.20499D-03, 8.20715D-03, 8.21608D-03,
     $     8.21754D-03, 8.23655D-03, 8.24078D-03, 8.24238D-03,
     $     8.25507D-03, 8.25726D-03, 8.26080D-03, 8.27432D-03,
     $     8.28000D-03, 8.28568D-03, 8.28957D-03, 8.29257D-03,
     $     8.29895D-03, 8.30324D-03, 8.30859D-03, 8.32322D-03,
     $     8.32684D-03, 8.33023D-03, 8.33750D-03, 8.34241D-03,
     $     8.34827D-03, 8.35133D-03, 8.35686D-03, 8.36455D-03,
     $     8.36726D-03, 8.37145D-03, 8.37453D-03, 8.37913D-03,
     $     8.38222D-03, 8.38776D-03, 8.38988D-03, 8.39427D-03,
     $     8.39795D-03, 8.40159D-03, 8.40457D-03, 8.40772D-03,
     $     8.41124D-03, 8.41587D-03, 8.41626D-03, 8.42233D-03,
     $     8.42803D-03, 8.43471D-03, 8.43613D-03, 8.43714D-03,
     $     8.43855D-03, 8.44090D-03, 8.44682D-03, 8.45181D-03,
     $     8.45247D-03, 8.45638D-03, 8.45945D-03, 8.46272D-03,
     $     8.46645D-03, 8.46766D-03, 8.46862D-03, 8.46917D-03,
     $     8.47154D-03, 8.47364D-03 /

!     Branching ratio of 2005
      data( br(11,k), k=1, 245) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00 /

!     Branching ratio of 2004
      data( br(12,k), k=1, 245) /
     $     2.85871D-03, 3.34279D-01, 3.82084D-01,
     $     4.10151D-01, 5.37401D-01, 5.87063D-01, 5.96599D-01,
     $     5.99024D-01, 6.04642D-01, 6.06045D-01, 6.10648D-01,
     $     6.20639D-01, 6.23060D-01, 6.21127D-01, 6.19719D-01,
     $     6.13648D-01, 6.04234D-01, 6.00964D-01, 5.91326D-01,
     $     5.84919D-01, 5.77282D-01, 5.68365D-01, 5.66164D-01,
     $     5.53567D-01, 5.48893D-01, 5.42129D-01, 5.34396D-01,
     $     5.32349D-01, 5.26637D-01, 5.22341D-01, 5.18253D-01,
     $     5.15540D-01, 5.13729D-01, 5.08267D-01, 5.05259D-01,
     $     5.02493D-01, 5.02001D-01, 4.97654D-01, 4.97053D-01,
     $     4.93383D-01, 4.92605D-01, 4.90921D-01, 4.89900D-01,
     $     4.89010D-01, 4.88110D-01, 4.86389D-01, 4.85206D-01,
     $     4.84686D-01, 4.83220D-01, 4.80261D-01, 4.78376D-01,
     $     4.77210D-01, 4.76090D-01, 4.75790D-01, 4.73987D-01,
     $     4.73090D-01, 4.72941D-01, 4.71336D-01, 4.70863D-01,
     $     4.70644D-01, 4.69799D-01, 4.68853D-01, 4.68447D-01,
     $     4.67452D-01, 4.65998D-01, 4.65558D-01, 4.65390D-01,
     $     4.65132D-01, 4.64685D-01, 4.64305D-01, 4.63524D-01,
     $     4.63179D-01, 4.62278D-01, 4.62054D-01, 4.61784D-01,
     $     4.61197D-01, 4.60614D-01, 4.60563D-01, 4.59982D-01,
     $     4.59693D-01, 4.59233D-01, 4.59030D-01, 4.58448D-01,
     $     4.58227D-01, 4.58008D-01, 4.57738D-01, 4.57576D-01,
     $     4.57074D-01, 4.56784D-01, 4.56663D-01, 4.56501D-01,
     $     4.56160D-01, 4.55707D-01, 4.55315D-01, 4.54808D-01,
     $     4.54762D-01, 4.54299D-01, 4.54070D-01, 4.53956D-01,
     $     4.53744D-01, 4.53548D-01, 4.53426D-01, 4.53226D-01,
     $     4.52962D-01, 4.52738D-01, 4.52434D-01, 4.52325D-01,
     $     4.52092D-01, 4.51699D-01, 4.51484D-01, 4.51348D-01,
     $     4.51124D-01, 4.51081D-01, 4.50882D-01, 4.50628D-01,
     $     4.50526D-01, 4.50501D-01, 4.50242D-01, 4.50078D-01,
     $     4.49896D-01, 4.49761D-01, 4.49513D-01, 4.49445D-01,
     $     4.49265D-01, 4.49095D-01, 4.48934D-01, 4.48718D-01,
     $     4.48474D-01, 4.48299D-01, 4.48077D-01, 4.48051D-01,
     $     4.47893D-01, 4.47696D-01, 4.47560D-01, 4.47437D-01,
     $     4.47279D-01, 4.47082D-01, 4.46987D-01, 4.46818D-01,
     $     4.46629D-01, 4.46437D-01, 4.46248D-01, 4.46058D-01,
     $     4.45858D-01, 4.45755D-01, 4.45665D-01, 4.45639D-01,
     $     4.45488D-01, 4.45349D-01, 4.45265D-01, 4.45063D-01,
     $     4.45015D-01, 4.44851D-01, 4.44769D-01, 4.44561D-01,
     $     4.44438D-01, 4.44296D-01, 4.44271D-01, 4.44062D-01,
     $     4.43993D-01, 4.43944D-01, 4.43840D-01, 4.43658D-01,
     $     4.43608D-01, 4.43504D-01, 4.43479D-01, 4.43338D-01,
     $     4.43313D-01, 4.43252D-01, 4.43135D-01, 4.43002D-01,
     $     4.42937D-01, 4.42920D-01, 4.42915D-01, 4.42877D-01,
     $     4.42825D-01, 4.42730D-01, 4.42694D-01, 4.42677D-01,
     $     4.42632D-01, 4.42603D-01, 4.42571D-01, 4.42562D-01,
     $     4.42519D-01, 4.42497D-01, 4.42492D-01, 4.42478D-01,
     $     4.42466D-01, 4.42431D-01, 4.42422D-01, 4.42419D-01,
     $     4.42347D-01, 4.42339D-01, 4.42304D-01, 4.42227D-01,
     $     4.42201D-01, 4.42146D-01, 4.42111D-01, 4.42094D-01,
     $     4.42053D-01, 4.42029D-01, 4.41988D-01, 4.41909D-01,
     $     4.41866D-01, 4.41838D-01, 4.41786D-01, 4.41760D-01,
     $     4.41711D-01, 4.41676D-01, 4.41625D-01, 4.41579D-01,
     $     4.41546D-01, 4.41524D-01, 4.41508D-01, 4.41481D-01,
     $     4.41441D-01, 4.41388D-01, 4.41372D-01, 4.41346D-01,
     $     4.41292D-01, 4.41264D-01, 4.41237D-01, 4.41210D-01,
     $     4.41175D-01, 4.41132D-01, 4.41129D-01, 4.41092D-01,
     $     4.41045D-01, 4.40997D-01, 4.40969D-01, 4.40954D-01,
     $     4.40940D-01, 4.40924D-01, 4.40891D-01, 4.40836D-01,
     $     4.40826D-01, 4.40793D-01, 4.40770D-01, 4.40722D-01,
     $     4.40677D-01, 4.40663D-01, 4.40652D-01, 4.40645D-01,
     $     4.40615D-01, 4.40586D-01 /

!     Branching ratio of 2003
      data( br(13,k), k=1, 245) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 7.81369D-03, 8.73887D-03,
     $     9.29077D-03, 1.10764D-02, 1.14202D-02, 1.21148D-02,
     $     1.26648D-02, 1.27887D-02, 1.27632D-02, 1.26951D-02,
     $     1.25774D-02, 1.25144D-02, 1.25876D-02, 1.30667D-02,
     $     1.34855D-02, 1.40485D-02, 1.49976D-02, 1.52443D-02,
     $     1.65635D-02, 1.69821D-02, 1.76208D-02, 1.84067D-02,
     $     1.86304D-02, 1.91695D-02, 1.95541D-02, 1.99269D-02,
     $     2.01929D-02, 2.04397D-02, 2.12052D-02, 2.16765D-02,
     $     2.21642D-02, 2.22587D-02, 2.33669D-02, 2.35333D-02,
     $     2.44788D-02, 2.47089D-02, 2.52394D-02, 2.55721D-02,
     $     2.58788D-02, 2.62206D-02, 2.68439D-02, 2.73183D-02,
     $     2.75726D-02, 2.83220D-02, 3.00010D-02, 3.10494D-02,
     $     3.17779D-02, 3.24508D-02, 3.26357D-02, 3.35819D-02,
     $     3.41394D-02, 3.42400D-02, 3.52962D-02, 3.56533D-02,
     $     3.58098D-02, 3.63778D-02, 3.70670D-02, 3.73888D-02,
     $     3.82678D-02, 3.95560D-02, 3.99456D-02, 4.00837D-02,
     $     4.03224D-02, 4.06841D-02, 4.09880D-02, 4.16147D-02,
     $     4.19331D-02, 4.27439D-02, 4.29461D-02, 4.31745D-02,
     $     4.36600D-02, 4.41612D-02, 4.42117D-02, 4.47575D-02,
     $     4.50365D-02, 4.54855D-02, 4.57343D-02, 4.64487D-02,
     $     4.67460D-02, 4.70404D-02, 4.73670D-02, 4.75828D-02,
     $     4.82869D-02, 4.86931D-02, 4.88894D-02, 4.91057D-02,
     $     4.95604D-02, 5.01648D-02, 5.06482D-02, 5.13063D-02,
     $     5.13600D-02, 5.19327D-02, 5.23092D-02, 5.24835D-02,
     $     5.28250D-02, 5.31247D-02, 5.33038D-02, 5.35568D-02,
     $     5.39288D-02, 5.42547D-02, 5.46479D-02, 5.47842D-02,
     $     5.51128D-02, 5.56567D-02, 5.59706D-02, 5.61720D-02,
     $     5.64826D-02, 5.65452D-02, 5.68081D-02, 5.71995D-02,
     $     5.73849D-02, 5.74336D-02, 5.78803D-02, 5.81291D-02,
     $     5.84400D-02, 5.86556D-02, 5.91012D-02, 5.92057D-02,
     $     5.95091D-02, 5.97863D-02, 6.00384D-02, 6.03400D-02,
     $     6.06945D-02, 6.09543D-02, 6.12825D-02, 6.13195D-02,
     $     6.15297D-02, 6.17942D-02, 6.19926D-02, 6.21532D-02,
     $     6.23616D-02, 6.26191D-02, 6.27447D-02, 6.29775D-02,
     $     6.32204D-02, 6.34945D-02, 6.37563D-02, 6.40342D-02,
     $     6.43047D-02, 6.44694D-02, 6.46124D-02, 6.46551D-02,
     $     6.48987D-02, 6.51208D-02, 6.52361D-02, 6.55103D-02,
     $     6.55706D-02, 6.57944D-02, 6.58853D-02, 6.61684D-02,
     $     6.63235D-02, 6.64859D-02, 6.65189D-02, 6.68093D-02,
     $     6.68877D-02, 6.69542D-02, 6.70973D-02, 6.73545D-02,
     $     6.74353D-02, 6.75865D-02, 6.76189D-02, 6.78086D-02,
     $     6.78519D-02, 6.79335D-02, 6.81168D-02, 6.83446D-02,
     $     6.84818D-02, 6.85232D-02, 6.85343D-02, 6.86011D-02,
     $     6.87085D-02, 6.89141D-02, 6.89748D-02, 6.90134D-02,
     $     6.91718D-02, 6.92717D-02, 6.93721D-02, 6.94315D-02,
     $     6.95551D-02, 6.96221D-02, 6.96476D-02, 6.97353D-02,
     $     6.97627D-02, 6.99483D-02, 6.99930D-02, 7.00086D-02,
     $     7.01759D-02, 7.02016D-02, 7.02566D-02, 7.04323D-02,
     $     7.04995D-02, 7.05770D-02, 7.06300D-02, 7.06607D-02,
     $     7.07276D-02, 7.07695D-02, 7.08297D-02, 7.09773D-02,
     $     7.10278D-02, 7.10688D-02, 7.11561D-02, 7.12025D-02,
     $     7.12807D-02, 7.13299D-02, 7.14089D-02, 7.14924D-02,
     $     7.15444D-02, 7.15919D-02, 7.16312D-02, 7.16909D-02,
     $     7.17447D-02, 7.18332D-02, 7.18633D-02, 7.19128D-02,
     $     7.19781D-02, 7.20193D-02, 7.20586D-02, 7.20992D-02,
     $     7.21431D-02, 7.22044D-02, 7.22085D-02, 7.22696D-02,
     $     7.23391D-02, 7.24331D-02, 7.24729D-02, 7.24979D-02,
     $     7.25260D-02, 7.25610D-02, 7.26321D-02, 7.27146D-02,
     $     7.27302D-02, 7.27835D-02, 7.28215D-02, 7.28810D-02,
     $     7.29380D-02, 7.29546D-02, 7.29674D-02, 7.29748D-02,
     $     7.30066D-02, 7.30350D-02 /

!     Branching ratio of 1003
      data( br(14,k), k=1, 245) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 1.04226D-03, 1.51260D-03,
     $     2.07422D-03, 5.06418D-03, 5.53191D-03, 6.34877D-03,
     $     7.08677D-03, 7.56452D-03, 7.66995D-03, 7.71637D-03,
     $     8.29095D-03, 9.61002D-03, 1.02068D-02, 1.18810D-02,
     $     1.27733D-02, 1.37018D-02, 1.43206D-02, 1.44286D-02,
     $     1.49131D-02, 1.50823D-02, 1.54459D-02, 1.55518D-02,
     $     1.56234D-02, 1.58440D-02, 1.60159D-02, 1.63305D-02,
     $     1.65323D-02, 1.66826D-02, 1.69432D-02, 1.71909D-02,
     $     1.74858D-02, 1.75527D-02, 1.83300D-02, 1.84486D-02,
     $     1.91573D-02, 1.93454D-02, 1.98758D-02, 2.02578D-02,
     $     2.06490D-02, 2.11064D-02, 2.19617D-02, 2.25860D-02,
     $     2.29169D-02, 2.38097D-02, 2.56560D-02, 2.68725D-02,
     $     2.77212D-02, 2.84494D-02, 2.86592D-02, 2.97172D-02,
     $     3.03765D-02, 3.05047D-02, 3.17375D-02, 3.21862D-02,
     $     3.23660D-02, 3.30340D-02, 3.38771D-02, 3.42801D-02,
     $     3.54006D-02, 3.69835D-02, 3.74658D-02, 3.76335D-02,
     $     3.79431D-02, 3.83722D-02, 3.87394D-02, 3.94922D-02,
     $     3.99028D-02, 4.09336D-02, 4.11947D-02, 4.14714D-02,
     $     4.20444D-02, 4.26556D-02, 4.27192D-02, 4.33716D-02,
     $     4.36942D-02, 4.42068D-02, 4.45218D-02, 4.53789D-02,
     $     4.57404D-02, 4.61035D-02, 4.64712D-02, 4.67269D-02,
     $     4.75339D-02, 4.79893D-02, 4.82239D-02, 4.84591D-02,
     $     4.89560D-02, 4.96227D-02, 5.01215D-02, 5.07987D-02,
     $     5.08488D-02, 5.13945D-02, 5.18041D-02, 5.19876D-02,
     $     5.23567D-02, 5.26814D-02, 5.28707D-02, 5.31215D-02,
     $     5.35255D-02, 5.38815D-02, 5.42570D-02, 5.43849D-02,
     $     5.47076D-02, 5.52347D-02, 5.55378D-02, 5.57385D-02,
     $     5.60441D-02, 5.61044D-02, 5.63442D-02, 5.67306D-02,
     $     5.69262D-02, 5.69779D-02, 5.74142D-02, 5.76377D-02,
     $     5.79182D-02, 5.81094D-02, 5.85299D-02, 5.86227D-02,
     $     5.88942D-02, 5.91399D-02, 5.93617D-02, 5.96099D-02,
     $     5.99014D-02, 6.01260D-02, 6.04112D-02, 6.04432D-02,
     $     6.06220D-02, 6.08476D-02, 6.10128D-02, 6.11455D-02,
     $     6.13161D-02, 6.15233D-02, 6.16300D-02, 6.18202D-02,
     $     6.20144D-02, 6.22496D-02, 6.24561D-02, 6.26767D-02,
     $     6.28761D-02, 6.30104D-02, 6.31191D-02, 6.31548D-02,
     $     6.33543D-02, 6.35318D-02, 6.36180D-02, 6.38265D-02,
     $     6.38724D-02, 6.40430D-02, 6.41051D-02, 6.43203D-02,
     $     6.44316D-02, 6.45414D-02, 6.45663D-02, 6.47863D-02,
     $     6.48398D-02, 6.48889D-02, 6.49930D-02, 6.51822D-02,
     $     6.52415D-02, 6.53424D-02, 6.53640D-02, 6.54903D-02,
     $     6.55206D-02, 6.55708D-02, 6.56925D-02, 6.58462D-02,
     $     6.59456D-02, 6.59763D-02, 6.59846D-02, 6.60265D-02,
     $     6.61010D-02, 6.62346D-02, 6.62684D-02, 6.62920D-02,
     $     6.63969D-02, 6.64625D-02, 6.65263D-02, 6.65655D-02,
     $     6.66402D-02, 6.66819D-02, 6.66983D-02, 6.67581D-02,
     $     6.67728D-02, 6.69007D-02, 6.69312D-02, 6.69421D-02,
     $     6.70449D-02, 6.70616D-02, 6.70921D-02, 6.71974D-02,
     $     6.72395D-02, 6.72827D-02, 6.73128D-02, 6.73324D-02,
     $     6.73739D-02, 6.74006D-02, 6.74361D-02, 6.75304D-02,
     $     6.75576D-02, 6.75819D-02, 6.76350D-02, 6.76656D-02,
     $     6.77108D-02, 6.77369D-02, 6.77798D-02, 6.78295D-02,
     $     6.78551D-02, 6.78831D-02, 6.79058D-02, 6.79397D-02,
     $     6.79658D-02, 6.80116D-02, 6.80280D-02, 6.80571D-02,
     $     6.80882D-02, 6.81117D-02, 6.81333D-02, 6.81558D-02,
     $     6.81795D-02, 6.82124D-02, 6.82147D-02, 6.82508D-02,
     $     6.82877D-02, 6.83366D-02, 6.83528D-02, 6.83636D-02,
     $     6.83767D-02, 6.83942D-02, 6.84322D-02, 6.84705D-02,
     $     6.84774D-02, 6.85044D-02, 6.85244D-02, 6.85508D-02,
     $     6.85781D-02, 6.85862D-02, 6.85926D-02, 6.85962D-02,
     $     6.86118D-02, 6.86256D-02 /

!     Branching ratio of 1002
      data( br(15,k), k=1, 245) /
     $     7.96173D-14, 3.58022D-09, 1.13860D-08,
     $     2.81358D-03, 2.28118D-02, 2.07658D-02, 2.17814D-02,
     $     2.29371D-02, 2.59950D-02, 2.60909D-02, 2.61247D-02,
     $     2.60400D-02, 2.61226D-02, 2.60121D-02, 2.58172D-02,
     $     2.53806D-02, 2.53236D-02, 2.54358D-02, 2.55086D-02,
     $     2.51810D-02, 2.44280D-02, 2.34746D-02, 2.32360D-02,
     $     2.17131D-02, 2.10287D-02, 1.99370D-02, 1.87840D-02,
     $     1.84428D-02, 1.74835D-02, 1.68736D-02, 1.63498D-02,
     $     1.60298D-02, 1.58018D-02, 1.53040D-02, 1.50407D-02,
     $     1.48301D-02, 1.47854D-02, 1.44342D-02, 1.43830D-02,
     $     1.40845D-02, 1.40491D-02, 1.39269D-02, 1.38562D-02,
     $     1.38076D-02, 1.37696D-02, 1.37708D-02, 1.38735D-02,
     $     1.39363D-02, 1.41798D-02, 1.48267D-02, 1.52490D-02,
     $     1.55416D-02, 1.58383D-02, 1.59230D-02, 1.63769D-02,
     $     1.66853D-02, 1.67440D-02, 1.75184D-02, 1.78301D-02,
     $     1.79599D-02, 1.84101D-02, 1.89296D-02, 1.91445D-02,
     $     1.96752D-02, 2.06251D-02, 2.09297D-02, 2.10421D-02,
     $     2.12573D-02, 2.15995D-02, 2.18973D-02, 2.24370D-02,
     $     2.26836D-02, 2.32600D-02, 2.34035D-02, 2.35693D-02,
     $     2.39495D-02, 2.43266D-02, 2.43647D-02, 2.47567D-02,
     $     2.49643D-02, 2.53058D-02, 2.55192D-02, 2.61424D-02,
     $     2.64377D-02, 2.67571D-02, 2.71278D-02, 2.73767D-02,
     $     2.81533D-02, 2.85942D-02, 2.88239D-02, 2.90637D-02,
     $     2.95718D-02, 3.02975D-02, 3.08559D-02, 3.15294D-02,
     $     3.15794D-02, 3.21101D-02, 3.24856D-02, 3.26630D-02,
     $     3.30262D-02, 3.33718D-02, 3.35736D-02, 3.38704D-02,
     $     3.43408D-02, 3.47335D-02, 3.51290D-02, 3.52638D-02,
     $     3.55759D-02, 3.61018D-02, 3.63952D-02, 3.65944D-02,
     $     3.69169D-02, 3.69770D-02, 3.72333D-02, 3.76224D-02,
     $     3.78086D-02, 3.78571D-02, 3.82784D-02, 3.85110D-02,
     $     3.87908D-02, 3.90026D-02, 3.94923D-02, 3.96169D-02,
     $     3.99854D-02, 4.03473D-02, 4.06982D-02, 4.11233D-02,
     $     4.16253D-02, 4.20245D-02, 4.25546D-02, 4.26163D-02,
     $     4.29796D-02, 4.34457D-02, 4.37789D-02, 4.40616D-02,
     $     4.44273D-02, 4.48884D-02, 4.51336D-02, 4.55732D-02,
     $     4.60441D-02, 4.66032D-02, 4.70954D-02, 4.76147D-02,
     $     4.81049D-02, 4.84257D-02, 4.86849D-02, 4.87708D-02,
     $     4.92637D-02, 4.97166D-02, 4.99542D-02, 5.05498D-02,
     $     5.06896D-02, 5.12106D-02, 5.14181D-02, 5.20939D-02,
     $     5.24556D-02, 5.28329D-02, 5.29161D-02, 5.36719D-02,
     $     5.38726D-02, 5.40488D-02, 5.44228D-02, 5.51044D-02,
     $     5.53133D-02, 5.56753D-02, 5.57563D-02, 5.62298D-02,
     $     5.63371D-02, 5.65266D-02, 5.69759D-02, 5.75357D-02,
     $     5.78899D-02, 5.79974D-02, 5.80266D-02, 5.81789D-02,
     $     5.84451D-02, 5.89100D-02, 5.90335D-02, 5.91124D-02,
     $     5.94454D-02, 5.96522D-02, 5.98550D-02, 5.99720D-02,
     $     6.02160D-02, 6.03559D-02, 6.04080D-02, 6.06011D-02,
     $     6.06535D-02, 6.10844D-02, 6.11918D-02, 6.12315D-02,
     $     6.16508D-02, 6.17190D-02, 6.18590D-02, 6.23267D-02,
     $     6.25210D-02, 6.27545D-02, 6.29166D-02, 6.30200D-02,
     $     6.32449D-02, 6.33887D-02, 6.35872D-02, 6.41049D-02,
     $     6.42787D-02, 6.44258D-02, 6.47396D-02, 6.49223D-02,
     $     6.51909D-02, 6.53496D-02, 6.56085D-02, 6.58997D-02,
     $     6.60492D-02, 6.62066D-02, 6.63316D-02, 6.65253D-02,
     $     6.66939D-02, 6.69774D-02, 6.70782D-02, 6.72625D-02,
     $     6.74876D-02, 6.76561D-02, 6.78107D-02, 6.79743D-02,
     $     6.81566D-02, 6.84010D-02, 6.84186D-02, 6.86799D-02,
     $     6.89510D-02, 6.92806D-02, 6.93982D-02, 6.94724D-02,
     $     6.95588D-02, 6.96759D-02, 6.99420D-02, 7.02449D-02,
     $     7.02986D-02, 7.05167D-02, 7.06802D-02, 7.09183D-02,
     $     7.11726D-02, 7.12527D-02, 7.13168D-02, 7.13541D-02,
     $     7.15191D-02, 7.16709D-02 /

      do j = 1, lines-1
         if ( u >= exen(j) .and. u < exen(j+1) ) then

            do i = 1, nimax
               if ( ifz(i)*1000+ ifa(i) == 1. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(1,j), exen(j+1),
     $                    br(1,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 1001. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(2,j), exen(j+1),
     $                    br(2,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 5008. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(3,j), exen(j+1),
     $                    br(3,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 1002. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(4,j), exen(j+1),
     $                    br(4,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 3008. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(5,j), exen(j+1),
     $                    br(5,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 1003. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(6,j), exen(j+1),
     $                    br(6,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 2003. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(7,j), exen(j+1),
     $                    br(7,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 2004. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(9,j), exen(j+1),
     $                    br(9,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 2006. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(10,j), exen(j+1),
     $                    br(10,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 2005. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(11,j), exen(j+1),
     $                    br(11,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 3006. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(12,j), exen(j+1),
     $                    br(12,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 3007. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(13,j), exen(j+1),
     $                    br(13,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 4007. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(14,j), exen(j+1),
     $                    br(14,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               end if

            end do

            do i = 1, nimax
               r(i) = r(i)/ integral
            end do

         end if

      end do

      end subroutine suppressalpha_b10


! ----------------------------------------------------------------------
!     Modify the decay width of C-12 in GDR
      subroutine suppressalpha_c12( u)
!     u: incident photon energy (MeV)
! ----------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

      common /exiejn/ nimax
!$OMP THREADPRIVATE(/exiejn/)
      common /std1/ r(70),s(70),sigma,rr(70)
!$OMP THREADPRIVATE(/std1/)
      common /ejectl/ omega(70),ifa(70),ifz(70)

!     implicit none
      double precision u
      double precision sigma

      double precision exen(500)  ! excitation energy (MeV)
      double precision br(30,500) ! branching ration (-)
      double precision brr        ! branching ration inbetween (-)
      integer lines               ! number of energy bins

      integer i, j, k           ! counters
      double precision integral

      double precision linearone

      integral = 0.

!      6012
!  (1) 6011    1
!  (2) 5011 1001
!  (3) 6010    2
!  (4) 5010 1002
!  (5) 4010 2002
!  (6) 6009    3
!  (7) 5009 1003
!  (8) 4009 2003
!  (9) 3009 3003
! (10) 5008 1004
! (11) 4008 2004
! (12) 3008 3004
! (13) 5007 1005
! (14) 4007 2005
! (15) 3007 3005
! (16) 2007 4005
! (17) 4006 2006
! (18) 3006 3006
! (19) 2006 4006
! (20) 3005 3007
! (21) 2005 4007
! (22) 2004 4008
! (23) 2003 4009
! (24) 1003 5009

      lines = 500
      data( exen(k), k=1, 500) /
     $     19.05300, 20.26200, 22.78700, 23.31000,
     $     24.16000, 26.57600, 26.78300, 27.59100, 28.17700,
     $     28.37500, 29.39500, 29.72600, 30.48300, 30.87200,
     $     31.26100, 31.62900, 31.87100, 32.58100, 33.12700,
     $     33.43000, 33.56100, 34.09800, 34.70100, 34.89700,
     $     35.03900, 35.31300, 35.55700, 35.59600, 35.67700,
     $     35.86400, 36.15900, 36.41500, 36.54700, 36.76600,
     $     37.03600, 37.29600, 37.39400, 37.55700, 37.67500,
     $     38.12800, 38.37300, 38.47200, 38.94300, 38.98000,
     $     39.04700, 39.28300, 39.37800, 39.42000, 39.61000,
     $     39.65900, 39.89400, 39.99100, 40.19100, 40.26900,
     $     40.40800, 40.54900, 40.63200, 40.79600, 40.96400,
     $     41.18000, 41.25500, 41.38500, 41.57400, 41.76300,
     $     41.88100, 42.15300, 42.22300, 42.36900, 42.52700,
     $     42.60500, 42.70000, 42.79900, 42.88100, 42.95000,
     $     43.01400, 43.22700, 43.34100, 43.45200, 43.47600,
     $     43.50200, 43.70500, 43.78900, 43.86600, 43.88500,
     $     43.91000, 43.99100, 44.09700, 44.23800, 44.30900,
     $     44.41800, 44.48300, 44.71300, 44.82400, 44.92400,
     $     44.97300, 45.08700, 45.13400, 45.17200, 45.28100,
     $     45.34900, 45.41400, 45.50900, 45.56000, 45.60500,
     $     45.66900, 45.73600, 45.87800, 45.93700, 45.96100,
     $     46.04100, 46.12100, 46.19200, 46.27200, 46.32800,
     $     46.39000, 46.46100, 46.51300, 46.55300, 46.58000,
     $     46.85300, 46.86800, 46.92500, 46.97700, 47.03600,
     $     47.06400, 47.13300, 47.23800, 47.27600, 47.44800,
     $     47.47400, 47.56400, 47.63200, 47.66100, 47.71300,
     $     47.80500, 47.84600, 47.91400, 48.03400, 48.06100,
     $     48.14200, 48.22200, 48.29000, 48.32900, 48.41400,
     $     48.45500, 48.51300, 48.57000, 48.59300, 48.70900,
     $     48.73700, 48.84900, 48.91300, 48.96600, 48.99800,
     $     49.07100, 49.17500, 49.19500, 49.25700, 49.34100,
     $     49.44600, 49.49300, 49.59200, 49.66500, 49.72100,
     $     49.78400, 49.85200, 49.90300, 49.95700, 50.01900,
     $     50.04500, 50.10300, 50.15800, 50.19500, 50.29000,
     $     50.36200, 50.36500, 50.47200, 50.53100, 50.56100,
     $     50.64900, 50.68900, 50.73400, 50.79300, 50.88300,
     $     50.91700, 51.02200, 51.04400, 51.11600, 51.16300,
     $     51.18700, 51.24600, 51.28500, 51.31100, 51.36200,
     $     51.42800, 51.49200, 51.50400, 51.54400, 51.60800,
     $     51.63100, 51.74200, 51.74900, 51.80900, 51.83400,
     $     51.88200, 51.93100, 51.96600, 52.07000, 52.09700,
     $     52.13200, 52.19900, 52.26700, 52.29100, 52.33300,
     $     52.37800, 52.44700, 52.48100, 52.56100, 52.59100,
     $     52.61500, 52.68600, 52.70900, 52.76300, 52.79900,
     $     52.82600, 52.87400, 52.90300, 52.95800, 53.00300,
     $     53.02900, 53.12100, 53.15500, 53.17900, 53.22300,
     $     53.25600, 53.31000, 53.34600, 53.42000, 53.46700,
     $     53.48700, 53.52200, 53.53600, 53.59200, 53.61800,
     $     53.66500, 53.71000, 53.74800, 53.79400, 53.84100,
     $     53.85200, 53.94900, 53.95400, 54.05200, 54.09500,
     $     54.12700, 54.14600, 54.19500, 54.20200, 54.26600,
     $     54.29500, 54.32800, 54.33800, 54.42600, 54.44200,
     $     54.51200, 54.55500, 54.56500, 54.62000, 54.65900,
     $     54.67300, 54.71000, 54.73800, 54.80500, 54.85000,
     $     54.87800, 54.90500, 54.93200, 54.94500, 55.00300,
     $     55.03800, 55.06400, 55.07900, 55.14100, 55.16100,
     $     55.23700, 55.24900, 55.28200, 55.30200, 55.33900,
     $     55.35000, 55.40700, 55.45500, 55.49600, 55.54600,
     $     55.55100, 55.57300, 55.60500, 55.66600, 55.69700,
     $     55.73000, 55.75200, 55.80700, 55.82500, 55.87400,
     $     55.90700, 55.92800, 55.93800, 55.96700, 55.99100,
     $     56.05500, 56.06900, 56.07100, 56.13000, 56.13800,
     $     56.16200, 56.19200, 56.28200, 56.29600, 56.29900,
     $     56.36600, 56.39000, 56.42900, 56.44300, 56.45000,
     $     56.47600, 56.52600, 56.56000, 56.59500, 56.60900,
     $     56.66100, 56.68400, 56.72400, 56.73500, 56.78500,
     $     56.82000, 56.88200, 56.90300, 56.93600, 56.95300,
     $     56.97600, 56.99600, 57.03000, 57.03400, 57.07700,
     $     57.10800, 57.13600, 57.15600, 57.24900, 57.26600,
     $     57.28700, 57.32100, 57.35200, 57.36100, 57.40000,
     $     57.41300, 57.46700, 57.48400, 57.52900, 57.54500,
     $     57.56300, 57.57800, 57.61000, 57.62200, 57.65300,
     $     57.69300, 57.70200, 57.74600, 57.76700, 57.79000,
     $     57.80800, 57.86300, 57.89300, 57.90400, 57.92600,
     $     57.97900, 58.03000, 58.03800, 58.05300, 58.07400,
     $     58.09900, 58.11900, 58.15100, 58.17200, 58.22400,
     $     58.28500, 58.29100, 58.33100, 58.35600, 58.36400,
     $     58.40600, 58.43400, 58.46100, 58.47700, 58.51100,
     $     58.52100, 58.53100, 58.59000, 58.61300, 58.63400,
     $     58.66100, 58.69000, 58.71700, 58.74000, 58.75500,
     $     58.79200, 58.81600, 58.84400, 58.85600, 58.90800,
     $     58.93300, 58.95200, 58.97500, 58.98900, 59.05000,
     $     59.08400, 59.09400, 59.12300, 59.13200, 59.16900,
     $     59.20000, 59.22900, 59.24100, 59.29200, 59.30100,
     $     59.31300, 59.32900, 59.36700, 59.41300, 59.44100,
     $     59.45800, 59.47000, 59.48700, 59.52800, 59.56000,
     $     59.56900, 59.59600, 59.61200, 59.63100, 59.66100,
     $     59.70900, 59.72000, 59.73100, 59.76000, 59.78400,
     $     59.81600, 59.84200, 59.86900, 59.88400, 59.90200,
     $     59.94200, 59.96600, 59.97100, 60.01200, 60.03200,
     $     60.05300, 60.07200, 60.10800, 60.14700, 60.16600,
     $     60.17100, 60.21600, 60.23700, 60.27000, 60.30400,
     $     60.31900, 60.32900, 60.38400, 60.39100, 60.41500,
     $     60.43700, 60.46600, 60.49100, 60.50400, 60.52200,
     $     60.56100, 60.61000, 60.63000, 60.67000, 60.68000,
     $     60.72100, 60.74200, 60.75700, 60.77700, 60.79500,
     $     60.81200, 60.82900, 60.84000, 60.86400, 60.89300,
     $     60.91400, 60.92900, 60.96800, 60.97400, 61.02300,
     $     61.04700, 61.07500, 61.12000, 61.12900, 61.15600,
     $     61.17700 /

!     Branching ratio of 6011
      data( br(1,k), k=1, 500) /
     $     0.00000D+00, 2.90149D-01, 3.97094D-01,
     $     3.91810D-01, 3.47421D-01, 3.18696D-01, 3.18738D-01,
     $     3.10635D-01, 2.84889D-01, 2.73204D-01, 2.17589D-01,
     $     2.05241D-01, 1.77918D-01, 1.63372D-01, 1.49546D-01,
     $     1.39288D-01, 1.32839D-01, 1.11627D-01, 9.82912D-02,
     $     9.19865D-02, 8.89238D-02, 7.91774D-02, 6.90022D-02,
     $     6.60562D-02, 6.38722D-02, 6.02143D-02, 5.72919D-02,
     $     5.68469D-02, 5.57698D-02, 5.36523D-02, 5.03587D-02,
     $     4.78521D-02, 4.64964D-02, 4.45334D-02, 4.20926D-02,
     $     3.98263D-02, 3.89947D-02, 3.76468D-02, 3.67091D-02,
     $     3.33137D-02, 3.16712D-02, 3.10357D-02, 2.82105D-02,
     $     2.79970D-02, 2.76137D-02, 2.62575D-02, 2.57267D-02,
     $     2.55033D-02, 2.45524D-02, 2.43181D-02, 2.31999D-02,
     $     2.27771D-02, 2.19121D-02, 2.15995D-02, 2.10626D-02,
     $     2.04945D-02, 2.01673D-02, 1.95766D-02, 1.89993D-02,
     $     1.82741D-02, 1.80179D-02, 1.76070D-02, 1.70120D-02,
     $     1.64465D-02, 1.61146D-02, 1.54042D-02, 1.52285D-02,
     $     1.48427D-02, 1.44624D-02, 1.42695D-02, 1.40526D-02,
     $     1.38227D-02, 1.36401D-02, 1.34856D-02, 1.33499D-02,
     $     1.28863D-02, 1.26538D-02, 1.24261D-02, 1.23717D-02,
     $     1.23185D-02, 1.19232D-02, 1.17639D-02, 1.16270D-02,
     $     1.15940D-02, 1.15509D-02, 1.14088D-02, 1.12247D-02,
     $     1.09924D-02, 1.08775D-02, 1.06996D-02, 1.05864D-02,
     $     1.02053D-02, 1.00332D-02, 9.87727D-03, 9.79974D-03,
     $     9.62987D-03, 9.56165D-03, 9.50519D-03, 9.34620D-03,
     $     9.25044D-03, 9.16000D-03, 9.02380D-03, 8.94605D-03,
     $     8.88254D-03, 8.79071D-03, 8.69522D-03, 8.50373D-03,
     $     8.42628D-03, 8.39500D-03, 8.28842D-03, 8.18475D-03,
     $     8.09440D-03, 7.99353D-03, 7.92340D-03, 7.84610D-03,
     $     7.75809D-03, 7.69405D-03, 7.64507D-03, 7.61213D-03,
     $     7.27516D-03, 7.25578D-03, 7.18774D-03, 7.12509D-03,
     $     7.05747D-03, 7.02597D-03, 6.94521D-03, 6.82209D-03,
     $     6.77999D-03, 6.59187D-03, 6.56463D-03, 6.47306D-03,
     $     6.40561D-03, 6.37719D-03, 6.32666D-03, 6.23652D-03,
     $     6.19506D-03, 6.12977D-03, 6.01731D-03, 5.99104D-03,
     $     5.91259D-03, 5.83883D-03, 5.77637D-03, 5.74144D-03,
     $     5.66496D-03, 5.62921D-03, 5.57799D-03, 5.52771D-03,
     $     5.50644D-03, 5.40778D-03, 5.38202D-03, 5.28732D-03,
     $     5.23544D-03, 5.19188D-03, 5.16636D-03, 5.10829D-03,
     $     5.02664D-03, 5.01137D-03, 4.96479D-03, 4.90029D-03,
     $     4.82270D-03, 4.78826D-03, 4.71862D-03, 4.66619D-03,
     $     4.62589D-03, 4.58098D-03, 4.53449D-03, 4.49968D-03,
     $     4.46391D-03, 4.42190D-03, 4.40428D-03, 4.36505D-03,
     $     4.32881D-03, 4.30530D-03, 4.24553D-03, 4.19970D-03,
     $     4.19785D-03, 4.13056D-03, 4.09489D-03, 4.07704D-03,
     $     4.02445D-03, 4.00098D-03, 3.97382D-03, 3.93846D-03,
     $     3.88531D-03, 3.86541D-03, 3.80529D-03, 3.79241D-03,
     $     3.75154D-03, 3.72585D-03, 3.71288D-03, 3.67946D-03,
     $     3.65831D-03, 3.64450D-03, 3.61756D-03, 3.58239D-03,
     $     3.54909D-03, 3.54294D-03, 3.52270D-03, 3.48971D-03,
     $     3.47814D-03, 3.42202D-03, 3.41847D-03, 3.38867D-03,
     $     3.37626D-03, 3.35300D-03, 3.32983D-03, 3.31283D-03,
     $     3.26248D-03, 3.25003D-03, 3.23412D-03, 3.20299D-03,
     $     3.17177D-03, 3.16091D-03, 3.14189D-03, 3.12209D-03,
     $     3.09169D-03, 3.07737D-03, 3.04374D-03, 3.03149D-03,
     $     3.02136D-03, 2.99267D-03, 2.98316D-03, 2.96160D-03,
     $     2.94743D-03, 2.93688D-03, 2.91777D-03, 2.90554D-03,
     $     2.88320D-03, 2.86560D-03, 2.85564D-03, 2.81977D-03,
     $     2.80670D-03, 2.79713D-03, 2.77982D-03, 2.76727D-03,
     $     2.74732D-03, 2.73391D-03, 2.70684D-03, 2.69003D-03,
     $     2.68273D-03, 2.67001D-03, 2.66474D-03, 2.64486D-03,
     $     2.63549D-03, 2.61880D-03, 2.60326D-03, 2.58960D-03,
     $     2.57360D-03, 2.55744D-03, 2.55341D-03, 2.52048D-03,
     $     2.51879D-03, 2.48572D-03, 2.47139D-03, 2.46081D-03,
     $     2.45441D-03, 2.43829D-03, 2.43604D-03, 2.41542D-03,
     $     2.40581D-03, 2.39478D-03, 2.39150D-03, 2.36397D-03,
     $     2.35907D-03, 2.33758D-03, 2.32417D-03, 2.32114D-03,
     $     2.30442D-03, 2.29261D-03, 2.28834D-03, 2.27738D-03,
     $     2.26903D-03, 2.24900D-03, 2.23605D-03, 2.22799D-03,
     $     2.22023D-03, 2.21240D-03, 2.20843D-03, 2.19155D-03,
     $     2.18143D-03, 2.17415D-03, 2.16987D-03, 2.15254D-03,
     $     2.14708D-03, 2.12614D-03, 2.12292D-03, 2.11412D-03,
     $     2.10883D-03, 2.09911D-03, 2.09623D-03, 2.08119D-03,
     $     2.06798D-03, 2.05722D-03, 2.04417D-03, 2.04288D-03,
     $     2.03710D-03, 2.02857D-03, 2.01257D-03, 2.00457D-03,
     $     1.99593D-03, 1.99021D-03, 1.97638D-03, 1.97193D-03,
     $     1.95969D-03, 1.95149D-03, 1.94616D-03, 1.94368D-03,
     $     1.93644D-03, 1.93062D-03, 1.91536D-03, 1.91205D-03,
     $     1.91158D-03, 1.89763D-03, 1.89575D-03, 1.89003D-03,
     $     1.88293D-03, 1.86153D-03, 1.85826D-03, 1.85756D-03,
     $     1.84181D-03, 1.83624D-03, 1.82721D-03, 1.82387D-03,
     $     1.82226D-03, 1.81639D-03, 1.80515D-03, 1.79746D-03,
     $     1.78975D-03, 1.78671D-03, 1.77526D-03, 1.77012D-03,
     $     1.76140D-03, 1.75904D-03, 1.74830D-03, 1.74083D-03,
     $     1.72751D-03, 1.72295D-03, 1.71590D-03, 1.71231D-03,
     $     1.70753D-03, 1.70337D-03, 1.69631D-03, 1.69549D-03,
     $     1.68646D-03, 1.68006D-03, 1.67434D-03, 1.67019D-03,
     $     1.65125D-03, 1.64773D-03, 1.64351D-03, 1.63654D-03,
     $     1.63027D-03, 1.62849D-03, 1.62085D-03, 1.61833D-03,
     $     1.60777D-03, 1.60445D-03, 1.59578D-03, 1.59270D-03,
     $     1.58928D-03, 1.58645D-03, 1.58035D-03, 1.57809D-03,
     $     1.57212D-03, 1.56439D-03, 1.56269D-03, 1.55442D-03,
     $     1.55045D-03, 1.54618D-03, 1.54283D-03, 1.53249D-03,
     $     1.52703D-03, 1.52504D-03, 1.52092D-03, 1.51122D-03,
     $     1.50178D-03, 1.50029D-03, 1.49756D-03, 1.49372D-03,
     $     1.48926D-03, 1.48573D-03, 1.48006D-03, 1.47638D-03,
     $     1.46729D-03, 1.45650D-03, 1.45547D-03, 1.44857D-03,
     $     1.44428D-03, 1.44289D-03, 1.43564D-03, 1.43076D-03,
     $     1.42599D-03, 1.42328D-03, 1.41752D-03, 1.41585D-03,
     $     1.41416D-03, 1.40431D-03, 1.40052D-03, 1.39702D-03,
     $     1.39241D-03, 1.38758D-03, 1.38313D-03, 1.37933D-03,
     $     1.37689D-03, 1.37096D-03, 1.36712D-03, 1.36260D-03,
     $     1.36070D-03, 1.35249D-03, 1.34854D-03, 1.34560D-03,
     $     1.34199D-03, 1.33979D-03, 1.33022D-03, 1.32500D-03,
     $     1.32345D-03, 1.31888D-03, 1.31749D-03, 1.31182D-03,
     $     1.30717D-03, 1.30273D-03, 1.30091D-03, 1.29318D-03,
     $     1.29184D-03, 1.29002D-03, 1.28765D-03, 1.28204D-03,
     $     1.27533D-03, 1.27133D-03, 1.26889D-03, 1.26716D-03,
     $     1.26469D-03, 1.25870D-03, 1.25408D-03, 1.25280D-03,
     $     1.24900D-03, 1.24673D-03, 1.24407D-03, 1.23978D-03,
     $     1.23292D-03, 1.23139D-03, 1.22987D-03, 1.22585D-03,
     $     1.22244D-03, 1.21797D-03, 1.21440D-03, 1.21075D-03,
     $     1.20870D-03, 1.20625D-03, 1.20086D-03, 1.19759D-03,
     $     1.19692D-03, 1.19145D-03, 1.18879D-03, 1.18594D-03,
     $     1.18342D-03, 1.17867D-03, 1.17345D-03, 1.17093D-03,
     $     1.17025D-03, 1.16441D-03, 1.16170D-03, 1.15736D-03,
     $     1.15294D-03, 1.15103D-03, 1.14977D-03, 1.14270D-03,
     $     1.14179D-03, 1.13870D-03, 1.13591D-03, 1.13230D-03,
     $     1.12919D-03, 1.12755D-03, 1.12528D-03, 1.12041D-03,
     $     1.11439D-03, 1.11197D-03, 1.10719D-03, 1.10600D-03,
     $     1.10112D-03, 1.09855D-03, 1.09676D-03, 1.09441D-03,
     $     1.09222D-03, 1.09010D-03, 1.08806D-03, 1.08677D-03,
     $     1.08394D-03, 1.08051D-03, 1.07803D-03, 1.07626D-03,
     $     1.07171D-03, 1.07101D-03, 1.06534D-03, 1.06255D-03,
     $     1.05931D-03, 1.05420D-03, 1.05320D-03, 1.05020D-03,
     $     1.04783D-03 /

!     Branching ratio of 5011
      data( br(2,k), k=1, 500) /
     $     9.99649D-01, 7.09609D-01, 6.02255D-01,
     $     6.02583D-01, 5.35681D-01, 4.01912D-01, 3.98854D-01,
     $     3.75939D-01, 3.38176D-01, 3.22900D-01, 2.55976D-01,
     $     2.42358D-01, 2.12957D-01, 1.97184D-01, 1.82057D-01,
     $     1.70950D-01, 1.63885D-01, 1.39774D-01, 1.24426D-01,
     $     1.17110D-01, 1.13489D-01, 1.02044D-01, 8.98685D-02,
     $     8.63298D-02, 8.36852D-02, 7.92691D-02, 7.57379D-02,
     $     7.52000D-02, 7.38771D-02, 7.12974D-02, 6.72479D-02,
     $     6.41659D-02, 6.24797D-02, 6.00513D-02, 5.69808D-02,
     $     5.41258D-02, 5.30750D-02, 5.13670D-02, 5.01758D-02,
     $     4.58314D-02, 4.37211D-02, 4.29024D-02, 3.92426D-02,
     $     3.89641D-02, 3.84638D-02, 3.66835D-02, 3.59734D-02,
     $     3.56767D-02, 3.44209D-02, 3.41116D-02, 3.26294D-02,
     $     3.20690D-02, 3.09169D-02, 3.05008D-02, 2.97854D-02,
     $     2.90167D-02, 2.85753D-02, 2.77814D-02, 2.70045D-02,
     $     2.60245D-02, 2.56766D-02, 2.51194D-02, 2.43043D-02,
     $     2.35323D-02, 2.30791D-02, 2.21080D-02, 2.18673D-02,
     $     2.13358D-02, 2.08126D-02, 2.05459D-02, 2.02468D-02,
     $     1.99287D-02, 1.96761D-02, 1.94619D-02, 1.92740D-02,
     $     1.86290D-02, 1.83055D-02, 1.79877D-02, 1.79114D-02,
     $     1.78371D-02, 1.72841D-02, 1.70600D-02, 1.68681D-02,
     $     1.68219D-02, 1.67616D-02, 1.65615D-02, 1.63024D-02,
     $     1.59761D-02, 1.58144D-02, 1.55639D-02, 1.54033D-02,
     $     1.48636D-02, 1.46199D-02, 1.43984D-02, 1.42883D-02,
     $     1.40474D-02, 1.39506D-02, 1.38703D-02, 1.36442D-02,
     $     1.35081D-02, 1.33795D-02, 1.31855D-02, 1.30742D-02,
     $     1.29835D-02, 1.28524D-02, 1.27159D-02, 1.24424D-02,
     $     1.23317D-02, 1.22870D-02, 1.21345D-02, 1.19860D-02,
     $     1.18566D-02, 1.17120D-02, 1.16114D-02, 1.15005D-02,
     $     1.13741D-02, 1.12821D-02, 1.12117D-02, 1.11644D-02,
     $     1.06790D-02, 1.06511D-02, 1.05529D-02, 1.04625D-02,
     $     1.03650D-02, 1.03196D-02, 1.02028D-02, 1.00247D-02,
     $     9.96384D-03, 9.69183D-03, 9.65244D-03, 9.52005D-03,
     $     9.42251D-03, 9.38141D-03, 9.30831D-03, 9.17781D-03,
     $     9.11773D-03, 9.02315D-03, 8.86018D-03, 8.82206D-03,
     $     8.70819D-03, 8.60114D-03, 8.51045D-03, 8.45972D-03,
     $     8.34860D-03, 8.29666D-03, 8.22222D-03, 8.14907D-03,
     $     8.11811D-03, 7.97459D-03, 7.93703D-03, 7.79914D-03,
     $     7.72361D-03, 7.66016D-03, 7.62300D-03, 7.53842D-03,
     $     7.41945D-03, 7.39719D-03, 7.32932D-03, 7.23526D-03,
     $     7.12212D-03, 7.07187D-03, 6.97030D-03, 6.89376D-03,
     $     6.83491D-03, 6.76928D-03, 6.70139D-03, 6.65053D-03,
     $     6.59828D-03, 6.53688D-03, 6.51112D-03, 6.45375D-03,
     $     6.40076D-03, 6.36640D-03, 6.27899D-03, 6.21192D-03,
     $     6.20922D-03, 6.11070D-03, 6.05850D-03, 6.03236D-03,
     $     5.95537D-03, 5.92100D-03, 5.88122D-03, 5.82939D-03,
     $     5.75149D-03, 5.72233D-03, 5.63420D-03, 5.61530D-03,
     $     5.55537D-03, 5.51771D-03, 5.49870D-03, 5.44967D-03,
     $     5.41866D-03, 5.39840D-03, 5.35889D-03, 5.30729D-03,
     $     5.25842D-03, 5.24940D-03, 5.21970D-03, 5.17127D-03,
     $     5.15429D-03, 5.07188D-03, 5.06666D-03, 5.02289D-03,
     $     5.00465D-03, 4.97049D-03, 4.93646D-03, 4.91147D-03,
     $     4.83745D-03, 4.81915D-03, 4.79578D-03, 4.75001D-03,
     $     4.70410D-03, 4.68812D-03, 4.66015D-03, 4.63104D-03,
     $     4.58632D-03, 4.56525D-03, 4.51579D-03, 4.49777D-03,
     $     4.48287D-03, 4.44067D-03, 4.42667D-03, 4.39495D-03,
     $     4.37410D-03, 4.35857D-03, 4.33045D-03, 4.31244D-03,
     $     4.27954D-03, 4.25363D-03, 4.23897D-03, 4.18614D-03,
     $     4.16688D-03, 4.15279D-03, 4.12729D-03, 4.10880D-03,
     $     4.07942D-03, 4.05966D-03, 4.01979D-03, 3.99502D-03,
     $     3.98426D-03, 3.96553D-03, 3.95776D-03, 3.92845D-03,
     $     3.91464D-03, 3.89004D-03, 3.86713D-03, 3.84698D-03,
     $     3.82339D-03, 3.79956D-03, 3.79361D-03, 3.74504D-03,
     $     3.74255D-03, 3.69377D-03, 3.67263D-03, 3.65701D-03,
     $     3.64756D-03, 3.62378D-03, 3.62046D-03, 3.59002D-03,
     $     3.57583D-03, 3.55955D-03, 3.55471D-03, 3.51407D-03,
     $     3.50684D-03, 3.47512D-03, 3.45532D-03, 3.45084D-03,
     $     3.42616D-03, 3.40872D-03, 3.40241D-03, 3.38623D-03,
     $     3.37391D-03, 3.34432D-03, 3.32520D-03, 3.31330D-03,
     $     3.30183D-03, 3.29027D-03, 3.28440D-03, 3.25946D-03,
     $     3.24451D-03, 3.23375D-03, 3.22743D-03, 3.20183D-03,
     $     3.19376D-03, 3.16283D-03, 3.15806D-03, 3.14507D-03,
     $     3.13726D-03, 3.12289D-03, 3.11864D-03, 3.09642D-03,
     $     3.07689D-03, 3.06099D-03, 3.04169D-03, 3.03979D-03,
     $     3.03125D-03, 3.01862D-03, 2.99496D-03, 2.98313D-03,
     $     2.97036D-03, 2.96190D-03, 2.94144D-03, 2.93487D-03,
     $     2.91676D-03, 2.90463D-03, 2.89675D-03, 2.89308D-03,
     $     2.88237D-03, 2.87376D-03, 2.85118D-03, 2.84629D-03,
     $     2.84560D-03, 2.82496D-03, 2.82218D-03, 2.81371D-03,
     $     2.80320D-03, 2.77153D-03, 2.76669D-03, 2.76566D-03,
     $     2.74235D-03, 2.73410D-03, 2.72074D-03, 2.71579D-03,
     $     2.71341D-03, 2.70473D-03, 2.68809D-03, 2.67671D-03,
     $     2.66530D-03, 2.66080D-03, 2.64385D-03, 2.63623D-03,
     $     2.62332D-03, 2.61984D-03, 2.60393D-03, 2.59288D-03,
     $     2.57315D-03, 2.56639D-03, 2.55596D-03, 2.55065D-03,
     $     2.54357D-03, 2.53741D-03, 2.52695D-03, 2.52574D-03,
     $     2.51237D-03, 2.50289D-03, 2.49442D-03, 2.48828D-03,
     $     2.46022D-03, 2.45501D-03, 2.44876D-03, 2.43843D-03,
     $     2.42916D-03, 2.42651D-03, 2.41520D-03, 2.41146D-03,
     $     2.39582D-03, 2.39091D-03, 2.37806D-03, 2.37351D-03,
     $     2.36844D-03, 2.36424D-03, 2.35520D-03, 2.35185D-03,
     $     2.34300D-03, 2.33155D-03, 2.32904D-03, 2.31678D-03,
     $     2.31089D-03, 2.30456D-03, 2.29960D-03, 2.28428D-03,
     $     2.27617D-03, 2.27323D-03, 2.26712D-03, 2.25274D-03,
     $     2.23874D-03, 2.23653D-03, 2.23249D-03, 2.22679D-03,
     $     2.22019D-03, 2.21495D-03, 2.20655D-03, 2.20109D-03,
     $     2.18760D-03, 2.17161D-03, 2.17008D-03, 2.15985D-03,
     $     2.15349D-03, 2.15143D-03, 2.14068D-03, 2.13344D-03,
     $     2.12637D-03, 2.12234D-03, 2.11380D-03, 2.11132D-03,
     $     2.10881D-03, 2.09421D-03, 2.08859D-03, 2.08340D-03,
     $     2.07656D-03, 2.06939D-03, 2.06279D-03, 2.05715D-03,
     $     2.05353D-03, 2.04474D-03, 2.03905D-03, 2.03235D-03,
     $     2.02953D-03, 2.01735D-03, 2.01148D-03, 2.00712D-03,
     $     2.00177D-03, 1.99851D-03, 1.98431D-03, 1.97657D-03,
     $     1.97427D-03, 1.96749D-03, 1.96543D-03, 1.95701D-03,
     $     1.95012D-03, 1.94353D-03, 1.94083D-03, 1.92936D-03,
     $     1.92738D-03, 1.92467D-03, 1.92116D-03, 1.91283D-03,
     $     1.90288D-03, 1.89694D-03, 1.89332D-03, 1.89075D-03,
     $     1.88710D-03, 1.87820D-03, 1.87135D-03, 1.86945D-03,
     $     1.86381D-03, 1.86044D-03, 1.85648D-03, 1.85011D-03,
     $     1.83993D-03, 1.83766D-03, 1.83541D-03, 1.82944D-03,
     $     1.82437D-03, 1.81775D-03, 1.81244D-03, 1.80703D-03,
     $     1.80398D-03, 1.80034D-03, 1.79234D-03, 1.78749D-03,
     $     1.78650D-03, 1.77838D-03, 1.77443D-03, 1.77019D-03,
     $     1.76645D-03, 1.75939D-03, 1.75164D-03, 1.74789D-03,
     $     1.74688D-03, 1.73823D-03, 1.73420D-03, 1.72775D-03,
     $     1.72118D-03, 1.71835D-03, 1.71648D-03, 1.70599D-03,
     $     1.70464D-03, 1.70003D-03, 1.69589D-03, 1.69054D-03,
     $     1.68591D-03, 1.68349D-03, 1.68010D-03, 1.67287D-03,
     $     1.66393D-03, 1.66034D-03, 1.65324D-03, 1.65148D-03,
     $     1.64423D-03, 1.64041D-03, 1.63776D-03, 1.63427D-03,
     $     1.63101D-03, 1.62786D-03, 1.62484D-03, 1.62291D-03,
     $     1.61872D-03, 1.61362D-03, 1.60993D-03, 1.60730D-03,
     $     1.60055D-03, 1.59950D-03, 1.59108D-03, 1.58695D-03,
     $     1.58213D-03, 1.57454D-03, 1.57305D-03, 1.56860D-03,
     $     1.56508D-03 /

!     Branching ratio of 6010
      data( br(3,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 1.09251D-14,
     $     1.12251D-10, 7.89811D-07, 2.58175D-06, 7.41501D-05,
     $     1.30995D-04, 1.67227D-04, 2.37554D-04, 3.01981D-04,
     $     3.11609D-04, 3.81761D-04, 4.46793D-04, 5.60098D-04,
     $     6.06611D-04, 6.81858D-04, 7.36753D-04, 1.16235D-03,
     $     1.32129D-03, 1.34491D-03, 1.37490D-03, 1.39245D-03,
     $     1.50885D-03, 1.53796D-03, 1.55621D-03, 1.66216D-03,
     $     1.67038D-03, 1.68498D-03, 1.76601D-03, 2.00209D-03,
     $     2.06346D-03, 2.18122D-03, 2.19983D-03, 2.31535D-03,
     $     2.34406D-03, 2.44543D-03, 2.46913D-03, 2.50691D-03,
     $     2.75549D-03, 2.84704D-03, 2.94752D-03, 3.03309D-03,
     $     3.16227D-03, 3.21789D-03, 3.29324D-03, 3.53630D-03,
     $     3.68137D-03, 3.75383D-03, 3.88490D-03, 3.91727D-03,
     $     4.01494D-03, 4.09705D-03, 4.14878D-03, 4.19697D-03,
     $     4.25663D-03, 4.30109D-03, 4.34181D-03, 4.37305D-03,
     $     4.49602D-03, 4.54780D-03, 4.60152D-03, 4.61925D-03,
     $     4.63261D-03, 4.72866D-03, 4.80589D-03, 4.84930D-03,
     $     4.85856D-03, 4.87022D-03, 4.93599D-03, 5.00494D-03,
     $     5.07145D-03, 5.10163D-03, 5.14664D-03, 5.19733D-03,
     $     5.31341D-03, 5.35249D-03, 5.40196D-03, 5.41916D-03,
     $     5.44592D-03, 5.45487D-03, 5.46884D-03, 5.50498D-03,
     $     5.51854D-03, 5.52982D-03, 5.54864D-03, 5.57316D-03,
     $     5.58476D-03, 5.59968D-03, 5.61346D-03, 5.63240D-03,
     $     5.63844D-03, 5.64067D-03, 5.64979D-03, 5.65658D-03,
     $     5.66089D-03, 5.66439D-03, 5.66601D-03, 5.66709D-03,
     $     5.66763D-03, 5.66748D-03, 5.66722D-03, 5.66693D-03,
     $     5.66884D-03, 5.66975D-03, 5.66953D-03, 5.67036D-03,
     $     5.66947D-03, 5.66884D-03, 5.67807D-03, 5.69724D-03,
     $     5.70087D-03, 5.70673D-03, 5.70664D-03, 5.70445D-03,
     $     5.70166D-03, 5.70021D-03, 5.69723D-03, 5.69197D-03,
     $     5.69021D-03, 5.68504D-03, 5.67384D-03, 5.67183D-03,
     $     5.66561D-03, 5.65744D-03, 5.65030D-03, 5.64569D-03,
     $     5.63566D-03, 5.63011D-03, 5.62245D-03, 5.61750D-03,
     $     5.61594D-03, 5.60025D-03, 5.59994D-03, 5.58503D-03,
     $     5.57462D-03, 5.56624D-03, 5.56077D-03, 5.54827D-03,
     $     5.53031D-03, 5.52670D-03, 5.51530D-03, 5.50420D-03,
     $     5.48631D-03, 5.47817D-03, 5.45987D-03, 5.44666D-03,
     $     5.43867D-03, 5.43008D-03, 5.41761D-03, 5.40783D-03,
     $     5.39692D-03, 5.38468D-03, 5.37954D-03, 5.36809D-03,
     $     5.35693D-03, 5.34916D-03, 5.32915D-03, 5.31589D-03,
     $     5.31530D-03, 5.29517D-03, 5.28253D-03, 5.27583D-03,
     $     5.25573D-03, 5.24626D-03, 5.23575D-03, 5.22248D-03,
     $     5.20342D-03, 5.19548D-03, 5.17137D-03, 5.16636D-03,
     $     5.14921D-03, 5.13768D-03, 5.13177D-03, 5.11781D-03,
     $     5.10835D-03, 5.10198D-03, 5.08995D-03, 5.07505D-03,
     $     5.05985D-03, 5.05696D-03, 5.04723D-03, 5.03186D-03,
     $     5.02624D-03, 5.00079D-03, 4.99935D-03, 4.98534D-03,
     $     4.97984D-03, 4.96867D-03, 4.95673D-03, 4.94829D-03,
     $     4.92432D-03, 4.91782D-03, 4.90936D-03, 4.89436D-03,
     $     4.87929D-03, 4.87377D-03, 4.86407D-03, 4.85358D-03,
     $     4.83844D-03, 4.83059D-03, 4.81210D-03, 4.80509D-03,
     $     4.79956D-03, 4.78292D-03, 4.77806D-03, 4.76563D-03,
     $     4.75718D-03, 4.75082D-03, 4.73955D-03, 4.73339D-03,
     $     4.72128D-03, 4.71101D-03, 4.70496D-03, 4.68434D-03,
     $     4.67674D-03, 4.67143D-03, 4.66174D-03, 4.65424D-03,
     $     4.64162D-03, 4.63357D-03, 4.61614D-03, 4.60485D-03,
     $     4.60003D-03, 4.59157D-03, 4.58851D-03, 4.57513D-03,
     $     4.56885D-03, 4.55788D-03, 4.54707D-03, 4.53838D-03,
     $     4.52737D-03, 4.51605D-03, 4.51358D-03, 4.49051D-03,
     $     4.48931D-03, 4.46599D-03, 4.45571D-03, 4.44813D-03,
     $     4.44373D-03, 4.43202D-03, 4.43033D-03, 4.41530D-03,
     $     4.40863D-03, 4.40107D-03, 4.39870D-03, 4.37757D-03,
     $     4.37370D-03, 4.35673D-03, 4.34662D-03, 4.34422D-03,
     $     4.33095D-03, 4.32151D-03, 4.31811D-03, 4.30915D-03,
     $     4.30236D-03, 4.28643D-03, 4.27562D-03, 4.26887D-03,
     $     4.26236D-03, 4.25584D-03, 4.25281D-03, 4.23900D-03,
     $     4.23084D-03, 4.22469D-03, 4.22126D-03, 4.20665D-03,
     $     4.20191D-03, 4.18381D-03, 4.18095D-03, 4.17311D-03,
     $     4.16835D-03, 4.15954D-03, 4.15693D-03, 4.14332D-03,
     $     4.13221D-03, 4.12254D-03, 4.11064D-03, 4.10945D-03,
     $     4.10418D-03, 4.09683D-03, 4.08236D-03, 4.07495D-03,
     $     4.06711D-03, 4.06190D-03, 4.04879D-03, 4.04450D-03,
     $     4.03277D-03, 4.02487D-03, 4.01993D-03, 4.01758D-03,
     $     4.01068D-03, 4.00499D-03, 3.98982D-03, 3.98651D-03,
     $     3.98604D-03, 3.97205D-03, 3.97016D-03, 3.96446D-03,
     $     3.95746D-03, 3.93677D-03, 3.93353D-03, 3.93284D-03,
     $     3.91731D-03, 3.91173D-03, 3.90264D-03, 3.89946D-03,
     $     3.89785D-03, 3.89186D-03, 3.88031D-03, 3.87256D-03,
     $     3.86454D-03, 3.86133D-03, 3.84947D-03, 3.84416D-03,
     $     3.83493D-03, 3.83241D-03, 3.82086D-03, 3.81285D-03,
     $     3.79858D-03, 3.79373D-03, 3.78607D-03, 3.78213D-03,
     $     3.77680D-03, 3.77219D-03, 3.76438D-03, 3.76346D-03,
     $     3.75356D-03, 3.74637D-03, 3.73987D-03, 3.73519D-03,
     $     3.71354D-03, 3.70957D-03, 3.70468D-03, 3.69692D-03,
     $     3.68973D-03, 3.68765D-03, 3.67864D-03, 3.67564D-03,
     $     3.66312D-03, 3.65918D-03, 3.64878D-03, 3.64508D-03,
     $     3.64094D-03, 3.63749D-03, 3.63012D-03, 3.62736D-03,
     $     3.62028D-03, 3.61112D-03, 3.60907D-03, 3.59898D-03,
     $     3.59414D-03, 3.58887D-03, 3.58477D-03, 3.57225D-03,
     $     3.56541D-03, 3.56290D-03, 3.55782D-03, 3.54572D-03,
     $     3.53406D-03, 3.53221D-03, 3.52878D-03, 3.52394D-03,
     $     3.51824D-03, 3.51368D-03, 3.50637D-03, 3.50159D-03,
     $     3.48984D-03, 3.47593D-03, 3.47457D-03, 3.46557D-03,
     $     3.45995D-03, 3.45813D-03, 3.44870D-03, 3.44240D-03,
     $     3.43632D-03, 3.43273D-03, 3.42510D-03, 3.42286D-03,
     $     3.42063D-03, 3.40747D-03, 3.40233D-03, 3.39762D-03,
     $     3.39154D-03, 3.38509D-03, 3.37908D-03, 3.37394D-03,
     $     3.37060D-03, 3.36242D-03, 3.35710D-03, 3.35094D-03,
     $     3.34831D-03, 3.33688D-03, 3.33138D-03, 3.32723D-03,
     $     3.32224D-03, 3.31919D-03, 3.30599D-03, 3.29865D-03,
     $     3.29648D-03, 3.29022D-03, 3.28830D-03, 3.28034D-03,
     $     3.27373D-03, 3.26749D-03, 3.26494D-03, 3.25405D-03,
     $     3.25215D-03, 3.24958D-03, 3.24618D-03, 3.23812D-03,
     $     3.22840D-03, 3.22252D-03, 3.21895D-03, 3.21642D-03,
     $     3.21286D-03, 3.20423D-03, 3.19753D-03, 3.19566D-03,
     $     3.19004D-03, 3.18669D-03, 3.18273D-03, 3.17645D-03,
     $     3.16642D-03, 3.16413D-03, 3.16186D-03, 3.15584D-03,
     $     3.15084D-03, 3.14421D-03, 3.13886D-03, 3.13332D-03,
     $     3.13025D-03, 3.12656D-03, 3.11839D-03, 3.11347D-03,
     $     3.11245D-03, 3.10413D-03, 3.10006D-03, 3.09576D-03,
     $     3.09191D-03, 3.08465D-03, 3.07672D-03, 3.07286D-03,
     $     3.07183D-03, 3.06277D-03, 3.05854D-03, 3.05186D-03,
     $     3.04501D-03, 3.04202D-03, 3.04002D-03, 3.02898D-03,
     $     3.02758D-03, 3.02279D-03, 3.01839D-03, 3.01264D-03,
     $     3.00768D-03, 3.00510D-03, 3.00151D-03, 2.99375D-03,
     $     2.98408D-03, 2.98016D-03, 2.97237D-03, 2.97043D-03,
     $     2.96248D-03, 2.95834D-03, 2.95543D-03, 2.95157D-03,
     $     2.94805D-03, 2.94470D-03, 2.94140D-03, 2.93928D-03,
     $     2.93465D-03, 2.92904D-03, 2.92500D-03, 2.92212D-03,
     $     2.91465D-03, 2.91350D-03, 2.90415D-03, 2.89957D-03,
     $     2.89423D-03, 2.88572D-03, 2.88403D-03, 2.87898D-03,
     $     2.87503D-03 /

!     Branching ratio of 5010
      data( br(4,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 6.22112D-13, 1.75184D-12,
     $     4.32191D-11, 5.73038D-03, 1.13074D-02, 3.00245D-02,
     $     3.13276D-02, 3.07467D-02, 3.09023D-02, 3.17280D-02,
     $     3.15121D-02, 3.32463D-02, 5.07349D-02, 6.10538D-02,
     $     6.52949D-02, 6.89061D-02, 7.81103D-02, 9.21519D-02,
     $     9.64823D-02, 9.94581D-02, 1.03648D-01, 1.06821D-01,
     $     1.07291D-01, 1.08964D-01, 1.11254D-01, 1.14528D-01,
     $     1.16391D-01, 1.17802D-01, 1.19435D-01, 1.23043D-01,
     $     1.24832D-01, 1.25227D-01, 1.25715D-01, 1.25980D-01,
     $     1.27060D-01, 1.27406D-01, 1.27572D-01, 1.28036D-01,
     $     1.28038D-01, 1.28037D-01, 1.28165D-01, 1.28723D-01,
     $     1.28852D-01, 1.28968D-01, 1.28953D-01, 1.28886D-01,
     $     1.28801D-01, 1.28713D-01, 1.28651D-01, 1.28533D-01,
     $     1.28747D-01, 1.28780D-01, 1.28691D-01, 1.28533D-01,
     $     1.28325D-01, 1.28255D-01, 1.28122D-01, 1.28110D-01,
     $     1.27934D-01, 1.27784D-01, 1.27375D-01, 1.27263D-01,
     $     1.27016D-01, 1.26760D-01, 1.26616D-01, 1.26448D-01,
     $     1.26251D-01, 1.26081D-01, 1.25926D-01, 1.25784D-01,
     $     1.25201D-01, 1.24864D-01, 1.24500D-01, 1.24408D-01,
     $     1.24319D-01, 1.23652D-01, 1.23405D-01, 1.23181D-01,
     $     1.23127D-01, 1.23055D-01, 1.22829D-01, 1.22509D-01,
     $     1.22079D-01, 1.21852D-01, 1.21471D-01, 1.21213D-01,
     $     1.20260D-01, 1.19800D-01, 1.19367D-01, 1.19138D-01,
     $     1.18629D-01, 1.18418D-01, 1.18240D-01, 1.17724D-01,
     $     1.17399D-01, 1.17081D-01, 1.16575D-01, 1.16272D-01,
     $     1.16023D-01, 1.15651D-01, 1.15253D-01, 1.14432D-01,
     $     1.14089D-01, 1.13949D-01, 1.13457D-01, 1.12967D-01,
     $     1.12529D-01, 1.12026D-01, 1.11667D-01, 1.11261D-01,
     $     1.10788D-01, 1.10435D-01, 1.10161D-01, 1.09975D-01,
     $     1.07982D-01, 1.07862D-01, 1.07443D-01, 1.07054D-01,
     $     1.06633D-01, 1.06436D-01, 1.05924D-01, 1.05132D-01,
     $     1.04861D-01, 1.03627D-01, 1.03447D-01, 1.02837D-01,
     $     1.02382D-01, 1.02188D-01, 1.01841D-01, 1.01210D-01,
     $     1.00913D-01, 1.00442D-01, 9.96155D-02, 9.94175D-02,
     $     9.88199D-02, 9.82533D-02, 9.77658D-02, 9.74912D-02,
     $     9.68803D-02, 9.65921D-02, 9.61736D-02, 9.57574D-02,
     $     9.55784D-02, 9.47479D-02, 9.45249D-02, 9.37051D-02,
     $     9.32516D-02, 9.28661D-02, 9.26396D-02, 9.21201D-02,
     $     9.13810D-02, 9.12420D-02, 9.08168D-02, 9.02196D-02,
     $     8.94974D-02, 8.91737D-02, 8.85171D-02, 8.80138D-02,
     $     8.76218D-02, 8.71805D-02, 8.67223D-02, 8.63759D-02,
     $     8.60188D-02, 8.55944D-02, 8.54153D-02, 8.50142D-02,
     $     8.46420D-02, 8.44005D-02, 8.37810D-02, 8.32977D-02,
     $     8.32782D-02, 8.25573D-02, 8.21717D-02, 8.19775D-02,
     $     8.13988D-02, 8.11381D-02, 8.08334D-02, 8.04329D-02,
     $     7.98240D-02, 7.95947D-02, 7.88965D-02, 7.87454D-02,
     $     7.82651D-02, 7.79623D-02, 7.78091D-02, 7.74105D-02,
     $     7.71581D-02, 7.69930D-02, 7.66695D-02, 7.62439D-02,
     $     7.58395D-02, 7.57646D-02, 7.55178D-02, 7.51122D-02,
     $     7.49697D-02, 7.42706D-02, 7.42259D-02, 7.38510D-02,
     $     7.36937D-02, 7.33991D-02, 7.31051D-02, 7.28879D-02,
     $     7.22397D-02, 7.20794D-02, 7.18743D-02, 7.14692D-02,
     $     7.10607D-02, 7.09185D-02, 7.06688D-02, 7.04087D-02,
     $     7.00065D-02, 6.98173D-02, 6.93708D-02, 6.92078D-02,
     $     6.90722D-02, 6.86875D-02, 6.85583D-02, 6.82662D-02,
     $     6.80735D-02, 6.79296D-02, 6.76674D-02, 6.74973D-02,
     $     6.71862D-02, 6.69409D-02, 6.68020D-02, 6.62969D-02,
     $     6.61118D-02, 6.59755D-02, 6.57280D-02, 6.55485D-02,
     $     6.52629D-02, 6.50693D-02, 6.46776D-02, 6.44331D-02,
     $     6.43263D-02, 6.41401D-02, 6.40618D-02, 6.37686D-02,
     $     6.36297D-02, 6.33811D-02, 6.31493D-02, 6.29435D-02,
     $     6.27028D-02, 6.24590D-02, 6.23974D-02, 6.18966D-02,
     $     6.18708D-02, 6.13632D-02, 6.11420D-02, 6.09779D-02,
     $     6.08782D-02, 6.06273D-02, 6.05923D-02, 6.02692D-02,
     $     6.01173D-02, 5.99424D-02, 5.98905D-02, 5.94546D-02,
     $     5.93768D-02, 5.90343D-02, 5.88186D-02, 5.87699D-02,
     $     5.85008D-02, 5.83098D-02, 5.82406D-02, 5.80630D-02,
     $     5.79272D-02, 5.75993D-02, 5.73871D-02, 5.72546D-02,
     $     5.71267D-02, 5.69974D-02, 5.69311D-02, 5.66504D-02,
     $     5.64811D-02, 5.63595D-02, 5.62875D-02, 5.59965D-02,
     $     5.59047D-02, 5.55512D-02, 5.54967D-02, 5.53477D-02,
     $     5.52579D-02, 5.50924D-02, 5.50434D-02, 5.47862D-02,
     $     5.45577D-02, 5.43718D-02, 5.41455D-02, 5.41231D-02,
     $     5.40227D-02, 5.38729D-02, 5.35924D-02, 5.34519D-02,
     $     5.32994D-02, 5.31983D-02, 5.29535D-02, 5.28747D-02,
     $     5.26571D-02, 5.25109D-02, 5.24153D-02, 5.23709D-02,
     $     5.22410D-02, 5.21366D-02, 5.18622D-02, 5.18028D-02,
     $     5.17943D-02, 5.15422D-02, 5.15082D-02, 5.14044D-02,
     $     5.12750D-02, 5.08826D-02, 5.08224D-02, 5.08097D-02,
     $     5.05196D-02, 5.04169D-02, 5.02502D-02, 5.01880D-02,
     $     5.01582D-02, 5.00496D-02, 4.98410D-02, 4.96976D-02,
     $     4.95540D-02, 4.94972D-02, 4.92826D-02, 4.91860D-02,
     $     4.90218D-02, 4.89775D-02, 4.87746D-02, 4.86331D-02,
     $     4.83791D-02, 4.82918D-02, 4.81570D-02, 4.80883D-02,
     $     4.79965D-02, 4.79166D-02, 4.77802D-02, 4.77644D-02,
     $     4.75896D-02, 4.74654D-02, 4.73543D-02, 4.72735D-02,
     $     4.69030D-02, 4.68339D-02, 4.67509D-02, 4.66129D-02,
     $     4.64892D-02, 4.64538D-02, 4.63027D-02, 4.62528D-02,
     $     4.60431D-02, 4.59772D-02, 4.58043D-02, 4.57429D-02,
     $     4.56746D-02, 4.56180D-02, 4.54957D-02, 4.54504D-02,
     $     4.53301D-02, 4.51740D-02, 4.51397D-02, 4.49724D-02,
     $     4.48918D-02, 4.48052D-02, 4.47370D-02, 4.45258D-02,
     $     4.44140D-02, 4.43733D-02, 4.42889D-02, 4.40894D-02,
     $     4.38943D-02, 4.38634D-02, 4.38071D-02, 4.37276D-02,
     $     4.36352D-02, 4.35620D-02, 4.34442D-02, 4.33676D-02,
     $     4.31777D-02, 4.29519D-02, 4.29302D-02, 4.27853D-02,
     $     4.26949D-02, 4.26658D-02, 4.25124D-02, 4.24090D-02,
     $     4.23076D-02, 4.22499D-02, 4.21276D-02, 4.20920D-02,
     $     4.20560D-02, 4.18457D-02, 4.17648D-02, 4.16900D-02,
     $     4.15910D-02, 4.14869D-02, 4.13911D-02, 4.13093D-02,
     $     4.12567D-02, 4.11289D-02, 4.10460D-02, 4.09480D-02,
     $     4.09069D-02, 4.07288D-02, 4.06430D-02, 4.05791D-02,
     $     4.05003D-02, 4.04525D-02, 4.02430D-02, 4.01289D-02,
     $     4.00950D-02, 3.99944D-02, 3.99638D-02, 3.98391D-02,
     $     3.97369D-02, 3.96392D-02, 3.95989D-02, 3.94280D-02,
     $     3.93983D-02, 3.93579D-02, 3.93055D-02, 3.91810D-02,
     $     3.90322D-02, 3.89431D-02, 3.88887D-02, 3.88502D-02,
     $     3.87951D-02, 3.86607D-02, 3.85571D-02, 3.85284D-02,
     $     3.84431D-02, 3.83921D-02, 3.83321D-02, 3.82354D-02,
     $     3.80803D-02, 3.80457D-02, 3.80114D-02, 3.79204D-02,
     $     3.78428D-02, 3.77414D-02, 3.76601D-02, 3.75770D-02,
     $     3.75302D-02, 3.74742D-02, 3.73511D-02, 3.72762D-02,
     $     3.72609D-02, 3.71354D-02, 3.70743D-02, 3.70087D-02,
     $     3.69506D-02, 3.68410D-02, 3.67204D-02, 3.66621D-02,
     $     3.66463D-02, 3.65114D-02, 3.64485D-02, 3.63475D-02,
     $     3.62445D-02, 3.62001D-02, 3.61707D-02, 3.60057D-02,
     $     3.59843D-02, 3.59116D-02, 3.58462D-02, 3.57616D-02,
     $     3.56886D-02, 3.56501D-02, 3.55965D-02, 3.54819D-02,
     $     3.53400D-02, 3.52830D-02, 3.51701D-02, 3.51420D-02,
     $     3.50262D-02, 3.49652D-02, 3.49229D-02, 3.48671D-02,
     $     3.48148D-02, 3.47641D-02, 3.47156D-02, 3.46847D-02,
     $     3.46174D-02, 3.45354D-02, 3.44759D-02, 3.44335D-02,
     $     3.43244D-02, 3.43076D-02, 3.41713D-02, 3.41043D-02,
     $     3.40261D-02, 3.39028D-02, 3.38786D-02, 3.38063D-02,
     $     3.37488D-02 /

!     Branching ratio of 4010
      data( br(5,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 1.10392D-06,
     $     1.32858D-05, 1.76073D-04, 2.09301D-04, 4.53202D-04,
     $     9.73076D-04, 1.36210D-03, 6.67689D-03, 1.22250D-02,
     $     1.37666D-02, 1.42078D-02, 1.58432D-02, 1.75995D-02,
     $     1.86276D-02, 1.92839D-02, 2.03041D-02, 2.09449D-02,
     $     2.10300D-02, 2.11278D-02, 2.13989D-02, 2.16699D-02,
     $     2.19600D-02, 2.20526D-02, 2.23285D-02, 2.28705D-02,
     $     2.37174D-02, 2.40321D-02, 2.44859D-02, 2.47586D-02,
     $     2.55148D-02, 2.60604D-02, 2.62904D-02, 2.71994D-02,
     $     2.72647D-02, 2.73882D-02, 2.78221D-02, 2.79846D-02,
     $     2.80625D-02, 2.84142D-02, 2.84971D-02, 2.87879D-02,
     $     2.89085D-02, 2.91023D-02, 2.91899D-02, 2.93445D-02,
     $     2.94249D-02, 2.94685D-02, 2.96077D-02, 2.97459D-02,
     $     2.98813D-02, 2.99054D-02, 2.99873D-02, 3.00734D-02,
     $     3.01741D-02, 3.02471D-02, 3.04371D-02, 3.04849D-02,
     $     3.05257D-02, 3.06124D-02, 3.06282D-02, 3.06712D-02,
     $     3.06887D-02, 3.07058D-02, 3.07071D-02, 3.07175D-02,
     $     3.06739D-02, 3.06596D-02, 3.06215D-02, 3.05976D-02,
     $     3.05847D-02, 3.04908D-02, 3.04391D-02, 3.04054D-02,
     $     3.03977D-02, 3.03877D-02, 3.03400D-02, 3.02675D-02,
     $     3.01825D-02, 3.01365D-02, 3.00532D-02, 2.99772D-02,
     $     2.97426D-02, 2.96553D-02, 2.95691D-02, 2.95217D-02,
     $     2.94353D-02, 2.94014D-02, 2.93678D-02, 2.92702D-02,
     $     2.92124D-02, 2.91549D-02, 2.90481D-02, 2.89704D-02,
     $     2.89160D-02, 2.88296D-02, 2.87363D-02, 2.85588D-02,
     $     2.84861D-02, 2.84563D-02, 2.83450D-02, 2.82377D-02,
     $     2.81435D-02, 2.80356D-02, 2.79588D-02, 2.78719D-02,
     $     2.77702D-02, 2.76943D-02, 2.76351D-02, 2.75948D-02,
     $     2.71309D-02, 2.71000D-02, 2.70016D-02, 2.69071D-02,
     $     2.68098D-02, 2.67648D-02, 2.66386D-02, 2.64374D-02,
     $     2.63719D-02, 2.60715D-02, 2.60288D-02, 2.58881D-02,
     $     2.57848D-02, 2.57410D-02, 2.56626D-02, 2.55158D-02,
     $     2.54429D-02, 2.53336D-02, 2.51441D-02, 2.50960D-02,
     $     2.49499D-02, 2.48163D-02, 2.47002D-02, 2.46358D-02,
     $     2.44905D-02, 2.44235D-02, 2.43242D-02, 2.42246D-02,
     $     2.41795D-02, 2.39849D-02, 2.39278D-02, 2.37316D-02,
     $     2.36259D-02, 2.35336D-02, 2.34805D-02, 2.33576D-02,
     $     2.31816D-02, 2.31489D-02, 2.30498D-02, 2.29054D-02,
     $     2.27336D-02, 2.26562D-02, 2.25032D-02, 2.23822D-02,
     $     2.22872D-02, 2.21803D-02, 2.20726D-02, 2.19906D-02,
     $     2.19077D-02, 2.18062D-02, 2.17630D-02, 2.16656D-02,
     $     2.15764D-02, 2.15197D-02, 2.13736D-02, 2.12567D-02,
     $     2.12521D-02, 2.10773D-02, 2.09858D-02, 2.09402D-02,
     $     2.08029D-02, 2.07417D-02, 2.06682D-02, 2.05717D-02,
     $     2.04258D-02, 2.03708D-02, 2.02046D-02, 2.01676D-02,
     $     2.00518D-02, 1.99803D-02, 1.99442D-02, 1.98461D-02,
     $     1.97855D-02, 1.97463D-02, 1.96693D-02, 1.95665D-02,
     $     1.94693D-02, 1.94515D-02, 1.93929D-02, 1.92946D-02,
     $     1.92604D-02, 1.90912D-02, 1.90803D-02, 1.89897D-02,
     $     1.89516D-02, 1.88811D-02, 1.88114D-02, 1.87586D-02,
     $     1.85997D-02, 1.85613D-02, 1.85124D-02, 1.84137D-02,
     $     1.83139D-02, 1.82791D-02, 1.82177D-02, 1.81545D-02,
     $     1.80559D-02, 1.80105D-02, 1.79027D-02, 1.78638D-02,
     $     1.78305D-02, 1.77386D-02, 1.77072D-02, 1.76373D-02,
     $     1.75914D-02, 1.75571D-02, 1.74935D-02, 1.74506D-02,
     $     1.73735D-02, 1.73138D-02, 1.72803D-02, 1.71562D-02,
     $     1.71108D-02, 1.70765D-02, 1.70145D-02, 1.69704D-02,
     $     1.69010D-02, 1.68537D-02, 1.67584D-02, 1.66993D-02,
     $     1.66731D-02, 1.66274D-02, 1.66078D-02, 1.65365D-02,
     $     1.65022D-02, 1.64413D-02, 1.63851D-02, 1.63338D-02,
     $     1.62746D-02, 1.62148D-02, 1.61991D-02, 1.60761D-02,
     $     1.60698D-02, 1.59444D-02, 1.58899D-02, 1.58495D-02,
     $     1.58246D-02, 1.57626D-02, 1.57541D-02, 1.56749D-02,
     $     1.56369D-02, 1.55929D-02, 1.55799D-02, 1.54731D-02,
     $     1.54542D-02, 1.53703D-02, 1.53168D-02, 1.53049D-02,
     $     1.52386D-02, 1.51913D-02, 1.51741D-02, 1.51305D-02,
     $     1.50969D-02, 1.50154D-02, 1.49634D-02, 1.49308D-02,
     $     1.48993D-02, 1.48672D-02, 1.48503D-02, 1.47800D-02,
     $     1.47376D-02, 1.47076D-02, 1.46896D-02, 1.46174D-02,
     $     1.45947D-02, 1.45067D-02, 1.44933D-02, 1.44566D-02,
     $     1.44345D-02, 1.43939D-02, 1.43819D-02, 1.43182D-02,
     $     1.42603D-02, 1.42142D-02, 1.41578D-02, 1.41523D-02,
     $     1.41271D-02, 1.40894D-02, 1.40188D-02, 1.39836D-02,
     $     1.39450D-02, 1.39195D-02, 1.38584D-02, 1.38388D-02,
     $     1.37844D-02, 1.37477D-02, 1.37235D-02, 1.37124D-02,
     $     1.36796D-02, 1.36536D-02, 1.35855D-02, 1.35708D-02,
     $     1.35687D-02, 1.35060D-02, 1.34975D-02, 1.34715D-02,
     $     1.34391D-02, 1.33405D-02, 1.33254D-02, 1.33222D-02,
     $     1.32489D-02, 1.32230D-02, 1.31808D-02, 1.31649D-02,
     $     1.31574D-02, 1.31301D-02, 1.30776D-02, 1.30414D-02,
     $     1.30054D-02, 1.29912D-02, 1.29372D-02, 1.29127D-02,
     $     1.28713D-02, 1.28602D-02, 1.28092D-02, 1.27737D-02,
     $     1.27096D-02, 1.26874D-02, 1.26533D-02, 1.26360D-02,
     $     1.26129D-02, 1.25929D-02, 1.25586D-02, 1.25547D-02,
     $     1.25104D-02, 1.24791D-02, 1.24512D-02, 1.24307D-02,
     $     1.23371D-02, 1.23195D-02, 1.22986D-02, 1.22635D-02,
     $     1.22321D-02, 1.22232D-02, 1.21852D-02, 1.21727D-02,
     $     1.21198D-02, 1.21031D-02, 1.20595D-02, 1.20440D-02,
     $     1.20269D-02, 1.20126D-02, 1.19818D-02, 1.19704D-02,
     $     1.19398D-02, 1.19000D-02, 1.18914D-02, 1.18490D-02,
     $     1.18285D-02, 1.18066D-02, 1.17894D-02, 1.17357D-02,
     $     1.17075D-02, 1.16973D-02, 1.16757D-02, 1.16251D-02,
     $     1.15752D-02, 1.15673D-02, 1.15530D-02, 1.15326D-02,
     $     1.15092D-02, 1.14907D-02, 1.14608D-02, 1.14415D-02,
     $     1.13934D-02, 1.13357D-02, 1.13302D-02, 1.12934D-02,
     $     1.12705D-02, 1.12630D-02, 1.12239D-02, 1.11974D-02,
     $     1.11712D-02, 1.11564D-02, 1.11252D-02, 1.11162D-02,
     $     1.11070D-02, 1.10534D-02, 1.10329D-02, 1.10137D-02,
     $     1.09881D-02, 1.09615D-02, 1.09370D-02, 1.09160D-02,
     $     1.09025D-02, 1.08700D-02, 1.08489D-02, 1.08239D-02,
     $     1.08134D-02, 1.07681D-02, 1.07461D-02, 1.07299D-02,
     $     1.07098D-02, 1.06976D-02, 1.06440D-02, 1.06148D-02,
     $     1.06061D-02, 1.05802D-02, 1.05724D-02, 1.05404D-02,
     $     1.05143D-02, 1.04892D-02, 1.04788D-02, 1.04349D-02,
     $     1.04273D-02, 1.04168D-02, 1.04034D-02, 1.03715D-02,
     $     1.03334D-02, 1.03108D-02, 1.02969D-02, 1.02870D-02,
     $     1.02729D-02, 1.02384D-02, 1.02118D-02, 1.02044D-02,
     $     1.01827D-02, 1.01696D-02, 1.01543D-02, 1.01294D-02,
     $     1.00894D-02, 1.00806D-02, 1.00718D-02, 1.00485D-02,
     $     1.00284D-02, 1.00023D-02, 9.98145D-03, 9.96024D-03,
     $     9.94824D-03, 9.93388D-03, 9.90235D-03, 9.88305D-03,
     $     9.87914D-03, 9.84699D-03, 9.83132D-03, 9.81433D-03,
     $     9.79942D-03, 9.77130D-03, 9.74009D-03, 9.72502D-03,
     $     9.72091D-03, 9.68624D-03, 9.67010D-03, 9.64395D-03,
     $     9.61736D-03, 9.60598D-03, 9.59845D-03, 9.55588D-03,
     $     9.55035D-03, 9.53156D-03, 9.51471D-03, 9.49303D-03,
     $     9.47426D-03, 9.46434D-03, 9.45048D-03, 9.42086D-03,
     $     9.38434D-03, 9.36970D-03, 9.34081D-03, 9.33364D-03,
     $     9.30395D-03, 9.28812D-03, 9.27726D-03, 9.26294D-03,
     $     9.24939D-03, 9.23609D-03, 9.22357D-03, 9.21560D-03,
     $     9.19825D-03, 9.17702D-03, 9.16164D-03, 9.15067D-03,
     $     9.12249D-03, 9.11812D-03, 9.08290D-03, 9.06553D-03,
     $     9.04521D-03, 9.01333D-03, 9.00709D-03, 8.98847D-03,
     $     8.97356D-03 /

!     Branching ratio of 6009
      data( br(6,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 5.99728D-26, 3.36300D-25, 2.31515D-15,
     $     7.75881D-15, 1.09494D-14, 1.22178D-14, 1.63236D-14,
     $     1.67297D-14, 2.73814D-14, 7.19379D-14, 7.64302D-14,
     $     1.17829D-13, 1.66537D-13, 1.92564D-13, 2.23691D-13,
     $     2.71040D-13, 5.14876D-13, 5.76192D-13, 6.40878D-13,
     $     7.19949D-13, 8.53655D-13, 1.03198D-12, 1.49796D-12,
     $     1.91542D-12, 2.06857D-12, 2.25720D-12, 2.84727D-12,
     $     2.95625D-12, 3.83324D-12, 3.90802D-12, 4.06622D-12,
     $     4.14957D-12, 4.27861D-12, 4.31140D-12, 5.36155D-12,
     $     1.00948D-11, 1.17079D-11, 1.43964D-11, 1.45686D-11,
     $     1.63242D-11, 2.26435D-11, 3.19528D-11, 3.58408D-11,
     $     4.44461D-11, 5.02041D-11, 5.75778D-11, 5.90462D-11,
     $     6.89349D-11, 7.63241D-11, 8.68674D-11, 9.00489D-11,
     $     1.02673D-10, 1.07733D-10, 1.16528D-10, 1.17977D-10,
     $     1.18173D-10, 1.33734D-10, 1.35608D-10, 1.50060D-10,
     $     1.73258D-10, 2.91943D-10, 3.07644D-10, 3.10329D-10,
     $     4.28647D-10, 4.66772D-10, 5.36708D-10, 6.06224D-10,
     $     6.23602D-10, 6.67647D-10, 7.58138D-10, 8.97215D-10,
     $     9.68525D-10, 9.88123D-10, 1.21029D-09, 1.32903D-09,
     $     1.47089D-09, 1.49459D-09, 1.65353D-09, 1.85656D-09,
     $     2.31585D-09, 2.52027D-09, 2.76181D-09, 2.86091D-09,
     $     2.95413D-09, 3.06392D-09, 3.30759D-09, 3.32376D-09,
     $     4.03847D-09, 4.37463D-09, 4.60611D-09, 4.85988D-09,
     $     6.00907D-09, 6.41491D-09, 6.68859D-09, 8.45195D-09,
     $     9.39558D-09, 9.58407D-09, 1.01499D-08, 1.03009D-08,
     $     1.12777D-08, 1.16447D-08, 1.24243D-08, 1.27098D-08,
     $     1.29406D-08, 1.31086D-08, 1.38624D-08, 1.40587D-08,
     $     1.64411D-08, 1.91938D-08, 1.95730D-08, 2.13195D-08,
     $     2.24626D-08, 2.32629D-08, 2.47374D-08, 3.19483D-08,
     $     3.36439D-08, 3.41149D-08, 3.62334D-08, 4.03404D-08,
     $     4.90902D-08, 5.03612D-08, 5.18941D-08, 5.44087D-08,
     $     5.62469D-08, 5.74202D-08, 5.97687D-08, 6.10254D-08,
     $     6.96495D-08, 7.93889D-08, 7.99867D-08, 8.90936D-08,
     $     9.43058D-08, 9.58755D-08, 1.05085D-07, 1.18407D-07,
     $     1.30467D-07, 1.34420D-07, 1.41836D-07, 1.43705D-07,
     $     1.48329D-07, 1.70123D-07, 1.74776D-07, 1.80204D-07,
     $     1.99297D-07, 2.14465D-07, 2.24704D-07, 2.32991D-07,
     $     2.37438D-07, 2.45510D-07, 2.50989D-07, 2.72729D-07,
     $     2.78089D-07, 2.95919D-07, 3.04513D-07, 3.08961D-07,
     $     3.19617D-07, 3.25060D-07, 3.90357D-07, 4.09792D-07,
     $     4.15718D-07, 4.53180D-07, 4.62189D-07, 4.90320D-07,
     $     5.06009D-07, 5.29749D-07, 5.46558D-07, 6.01009D-07,
     $     6.10505D-07, 6.23386D-07, 6.35552D-07, 6.62546D-07,
     $     6.89231D-07, 7.01816D-07, 7.10598D-07, 7.17169D-07,
     $     7.40687D-07, 7.92788D-07, 8.35756D-07, 8.46089D-07,
     $     8.69678D-07, 8.83922D-07, 8.97900D-07, 9.38573D-07,
     $     1.01248D-06, 1.02334D-06, 1.03308D-06, 1.05955D-06,
     $     1.10788D-06, 1.16088D-06, 1.19827D-06, 1.22613D-06,
     $     1.25321D-06, 1.27921D-06, 1.33204D-06, 1.37716D-06,
     $     1.38367D-06, 1.45214D-06, 1.47916D-06, 1.52391D-06,
     $     1.56016D-06, 1.63519D-06, 1.72388D-06, 1.76273D-06,
     $     1.77360D-06, 1.83846D-06, 1.86625D-06, 1.94780D-06,
     $     2.03502D-06, 2.06117D-06, 2.07627D-06, 2.20471D-06,
     $     2.23577D-06, 2.32765D-06, 2.38078D-06, 2.43454D-06,
     $     2.47767D-06, 2.51857D-06, 2.56823D-06, 2.64905D-06,
     $     2.73302D-06, 2.76247D-06, 2.81543D-06, 2.82798D-06,
     $     2.91711D-06, 2.96623D-06, 2.99234D-06, 3.02299D-06,
     $     3.08978D-06, 3.17099D-06, 3.21954D-06, 3.24530D-06,
     $     3.29710D-06, 3.36140D-06, 3.43679D-06, 3.49696D-06,
     $     3.60850D-06, 3.62406D-06, 3.77029D-06, 3.84753D-06,
     $     3.94630D-06, 4.05793D-06, 4.07630D-06, 4.12846D-06,
     $     4.20846D-06 /

!     Branching ratio of 5009
      data( br(7,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00 /

!     Branching ratio of 4009
      data( br(8,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 9.33579D-07, 1.34966D-04,
     $     2.02967D-02, 4.21122D-02, 4.71611D-02, 5.54527D-02,
     $     5.48121D-02, 5.01659D-02, 4.66681D-02, 4.30683D-02,
     $     4.03341D-02, 3.85588D-02, 3.25892D-02, 2.87807D-02,
     $     2.69810D-02, 2.61043D-02, 2.34703D-02, 2.34578D-02,
     $     2.40036D-02, 2.42883D-02, 2.47728D-02, 2.51340D-02,
     $     2.51892D-02, 2.52329D-02, 2.56351D-02, 2.71682D-02,
     $     2.92297D-02, 3.02170D-02, 3.21313D-02, 3.52513D-02,
     $     3.92659D-02, 4.09786D-02, 4.39456D-02, 4.60878D-02,
     $     5.31388D-02, 5.60420D-02, 5.70071D-02, 6.03348D-02,
     $     6.05541D-02, 6.09417D-02, 6.21272D-02, 6.27099D-02,
     $     6.29415D-02, 6.38311D-02, 6.40304D-02, 6.47626D-02,
     $     6.50525D-02, 6.55454D-02, 6.57596D-02, 6.61500D-02,
     $     6.67493D-02, 6.70096D-02, 6.75906D-02, 6.82501D-02,
     $     6.91295D-02, 6.94172D-02, 6.99840D-02, 7.11142D-02,
     $     7.19999D-02, 7.25622D-02, 7.38803D-02, 7.42185D-02,
     $     7.48326D-02, 7.55403D-02, 7.58376D-02, 7.62199D-02,
     $     7.65631D-02, 7.68410D-02, 7.70492D-02, 7.72519D-02,
     $     7.78215D-02, 7.81463D-02, 7.84195D-02, 7.84541D-02,
     $     7.85094D-02, 7.89214D-02, 7.92026D-02, 7.93934D-02,
     $     7.94371D-02, 7.94935D-02, 7.97548D-02, 8.00073D-02,
     $     8.03017D-02, 8.04413D-02, 8.06380D-02, 8.08178D-02,
     $     8.13053D-02, 8.15786D-02, 8.19040D-02, 8.20140D-02,
     $     8.22790D-02, 8.23838D-02, 8.24914D-02, 8.27711D-02,
     $     8.28986D-02, 8.30033D-02, 8.31004D-02, 8.31854D-02,
     $     8.32373D-02, 8.32728D-02, 8.32889D-02, 8.33238D-02,
     $     8.33319D-02, 8.33339D-02, 8.33194D-02, 8.33080D-02,
     $     8.32983D-02, 8.32846D-02, 8.32722D-02, 8.32555D-02,
     $     8.32318D-02, 8.32109D-02, 8.31925D-02, 8.31789D-02,
     $     8.29149D-02, 8.28897D-02, 8.28263D-02, 8.27606D-02,
     $     8.27017D-02, 8.26756D-02, 8.26460D-02, 8.26105D-02,
     $     8.25949D-02, 8.24359D-02, 8.24134D-02, 8.23457D-02,
     $     8.22975D-02, 8.22766D-02, 8.22378D-02, 8.21502D-02,
     $     8.20961D-02, 8.20257D-02, 8.19033D-02, 8.18644D-02,
     $     8.17425D-02, 8.16383D-02, 8.15409D-02, 8.14885D-02,
     $     8.13608D-02, 8.13042D-02, 8.12140D-02, 8.11427D-02,
     $     8.11036D-02, 8.09392D-02, 8.09018D-02, 8.07241D-02,
     $     8.06239D-02, 8.05276D-02, 8.04733D-02, 8.03428D-02,
     $     8.01478D-02, 8.01120D-02, 8.00043D-02, 7.98657D-02,
     $     7.96797D-02, 7.95934D-02, 7.94307D-02, 7.92920D-02,
     $     7.92006D-02, 7.91040D-02, 7.89928D-02, 7.89015D-02,
     $     7.88111D-02, 7.86906D-02, 7.86381D-02, 7.85164D-02,
     $     7.84051D-02, 7.83363D-02, 7.81534D-02, 7.80124D-02,
     $     7.80068D-02, 7.77904D-02, 7.76713D-02, 7.76109D-02,
     $     7.74231D-02, 7.73395D-02, 7.72346D-02, 7.71042D-02,
     $     7.69243D-02, 7.68489D-02, 7.66304D-02, 7.65784D-02,
     $     7.64140D-02, 7.63133D-02, 7.62625D-02, 7.61138D-02,
     $     7.60246D-02, 7.59675D-02, 7.58592D-02, 7.57152D-02,
     $     7.55722D-02, 7.55459D-02, 7.54599D-02, 7.53098D-02,
     $     7.52583D-02, 7.50149D-02, 7.50009D-02, 7.48701D-02,
     $     7.48182D-02, 7.47194D-02, 7.46192D-02, 7.45392D-02,
     $     7.43055D-02, 7.42490D-02, 7.41774D-02, 7.40360D-02,
     $     7.38916D-02, 7.38387D-02, 7.37428D-02, 7.36449D-02,
     $     7.34981D-02, 7.34295D-02, 7.32619D-02, 7.32022D-02,
     $     7.31488D-02, 7.30062D-02, 7.29616D-02, 7.28555D-02,
     $     7.27844D-02, 7.27311D-02, 7.26283D-02, 7.25608D-02,
     $     7.24398D-02, 7.23451D-02, 7.22915D-02, 7.20940D-02,
     $     7.20224D-02, 7.19658D-02, 7.18655D-02, 7.17939D-02,
     $     7.16809D-02, 7.16071D-02, 7.14497D-02, 7.13513D-02,
     $     7.13060D-02, 7.12266D-02, 7.11957D-02, 7.10755D-02,
     $     7.10149D-02, 7.09131D-02, 7.08168D-02, 7.07298D-02,
     $     7.06256D-02, 7.05193D-02, 7.04920D-02, 7.02774D-02,
     $     7.02662D-02, 7.00455D-02, 6.99490D-02, 6.98787D-02,
     $     6.98358D-02, 6.97259D-02, 6.97109D-02, 6.95753D-02,
     $     6.95096D-02, 6.94331D-02, 6.94097D-02, 6.92196D-02,
     $     6.91858D-02, 6.90331D-02, 6.89373D-02, 6.89156D-02,
     $     6.87930D-02, 6.87045D-02, 6.86716D-02, 6.85900D-02,
     $     6.85263D-02, 6.83739D-02, 6.82765D-02, 6.82146D-02,
     $     6.81540D-02, 6.80913D-02, 6.80588D-02, 6.79219D-02,
     $     6.78409D-02, 6.77828D-02, 6.77494D-02, 6.76085D-02,
     $     6.75640D-02, 6.73861D-02, 6.73590D-02, 6.72852D-02,
     $     6.72406D-02, 6.71584D-02, 6.71340D-02, 6.70025D-02,
     $     6.68850D-02, 6.67901D-02, 6.66721D-02, 6.66606D-02,
     $     6.66069D-02, 6.65309D-02, 6.63813D-02, 6.63058D-02,
     $     6.62230D-02, 6.61684D-02, 6.60372D-02, 6.59953D-02,
     $     6.58764D-02, 6.57961D-02, 6.57440D-02, 6.57201D-02,
     $     6.56480D-02, 6.55914D-02, 6.54435D-02, 6.54116D-02,
     $     6.54070D-02, 6.52695D-02, 6.52510D-02, 6.51930D-02,
     $     6.51229D-02, 6.49127D-02, 6.48800D-02, 6.48731D-02,
     $     6.47116D-02, 6.46538D-02, 6.45588D-02, 6.45237D-02,
     $     6.45070D-02, 6.44464D-02, 6.43280D-02, 6.42474D-02,
     $     6.41665D-02, 6.41346D-02, 6.40135D-02, 6.39567D-02,
     $     6.38615D-02, 6.38362D-02, 6.37184D-02, 6.36378D-02,
     $     6.34898D-02, 6.34379D-02, 6.33576D-02, 6.33167D-02,
     $     6.32626D-02, 6.32158D-02, 6.31363D-02, 6.31270D-02,
     $     6.30229D-02, 6.29481D-02, 6.28811D-02, 6.28310D-02,
     $     6.26036D-02, 6.25600D-02, 6.25087D-02, 6.24252D-02,
     $     6.23480D-02, 6.23261D-02, 6.22328D-02, 6.22019D-02,
     $     6.20700D-02, 6.20281D-02, 6.19188D-02, 6.18799D-02,
     $     6.18368D-02, 6.18011D-02, 6.17230D-02, 6.16942D-02,
     $     6.16178D-02, 6.15176D-02, 6.14957D-02, 6.13882D-02,
     $     6.13353D-02, 6.12793D-02, 6.12357D-02, 6.10994D-02,
     $     6.10275D-02, 6.10013D-02, 6.09446D-02, 6.08148D-02,
     $     6.06855D-02, 6.06644D-02, 6.06268D-02, 6.05728D-02,
     $     6.05114D-02, 6.04629D-02, 6.03838D-02, 6.03325D-02,
     $     6.02064D-02, 6.00507D-02, 6.00359D-02, 5.99379D-02,
     $     5.98764D-02, 5.98561D-02, 5.97512D-02, 5.96797D-02,
     $     5.96089D-02, 5.95691D-02, 5.94837D-02, 5.94589D-02,
     $     5.94340D-02, 5.92877D-02, 5.92310D-02, 5.91776D-02,
     $     5.91062D-02, 5.90330D-02, 5.89652D-02, 5.89064D-02,
     $     5.88688D-02, 5.87784D-02, 5.87192D-02, 5.86498D-02,
     $     5.86207D-02, 5.84936D-02, 5.84316D-02, 5.83861D-02,
     $     5.83306D-02, 5.82962D-02, 5.81472D-02, 5.80654D-02,
     $     5.80407D-02, 5.79680D-02, 5.79462D-02, 5.78558D-02,
     $     5.77822D-02, 5.77101D-02, 5.76811D-02, 5.75556D-02,
     $     5.75341D-02, 5.75040D-02, 5.74656D-02, 5.73737D-02,
     $     5.72639D-02, 5.71988D-02, 5.71586D-02, 5.71300D-02,
     $     5.70896D-02, 5.69900D-02, 5.69134D-02, 5.68923D-02,
     $     5.68295D-02, 5.67914D-02, 5.67468D-02, 5.66739D-02,
     $     5.65573D-02, 5.65314D-02, 5.65058D-02, 5.64371D-02,
     $     5.63782D-02, 5.63014D-02, 5.62401D-02, 5.61776D-02,
     $     5.61427D-02, 5.61004D-02, 5.60072D-02, 5.59499D-02,
     $     5.59382D-02, 5.58434D-02, 5.57967D-02, 5.57457D-02,
     $     5.57013D-02, 5.56180D-02, 5.55241D-02, 5.54786D-02,
     $     5.54659D-02, 5.53611D-02, 5.53121D-02, 5.52328D-02,
     $     5.51524D-02, 5.51180D-02, 5.50952D-02, 5.49656D-02,
     $     5.49491D-02, 5.48925D-02, 5.48412D-02, 5.47752D-02,
     $     5.47176D-02, 5.46875D-02, 5.46448D-02, 5.45532D-02,
     $     5.44400D-02, 5.43947D-02, 5.43054D-02, 5.42832D-02,
     $     5.41918D-02, 5.41419D-02, 5.41080D-02, 5.40634D-02,
     $     5.40211D-02, 5.39793D-02, 5.39400D-02, 5.39150D-02,
     $     5.38604D-02, 5.37929D-02, 5.37447D-02, 5.37104D-02,
     $     5.36211D-02, 5.36071D-02, 5.34951D-02, 5.34395D-02,
     $     5.33744D-02, 5.32714D-02, 5.32513D-02, 5.31912D-02,
     $     5.31432D-02 /

!     Branching ratio of 3009
      data( br(9,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 2.13079D-24, 3.00223D-24,
     $     4.10329D-18, 1.55769D-16, 2.51074D-15, 1.73223D-14,
     $     3.38081D-13, 1.25419D-12, 4.44778D-12, 1.12986D-11,
     $     1.59990D-11, 1.07862D-10, 1.55208D-10, 5.96491D-10,
     $     1.12954D-09, 1.78550D-09, 2.30436D-09, 3.94566D-09,
     $     7.85447D-09, 8.89001D-09, 1.28854D-08, 2.04256D-08,
     $     3.41359D-08, 4.21360D-08, 6.34866D-08, 8.37659D-08,
     $     1.02451D-07, 1.27039D-07, 1.58387D-07, 1.85037D-07,
     $     2.16353D-07, 2.55414D-07, 2.72759D-07, 3.12731D-07,
     $     3.52396D-07, 3.79658D-07, 4.51237D-07, 5.07325D-07,
     $     5.09737D-07, 5.99942D-07, 6.54859D-07, 6.84602D-07,
     $     7.80780D-07, 8.29296D-07, 8.87429D-07, 9.69832D-07,
     $     1.10718D-06, 1.16204D-06, 1.33992D-06, 1.37833D-06,
     $     1.50703D-06, 1.59345D-06, 1.63843D-06, 1.74989D-06,
     $     1.82626D-06, 1.87825D-06, 1.98200D-06, 2.12116D-06,
     $     2.26164D-06, 2.28878D-06, 2.38099D-06, 2.53341D-06,
     $     2.59024D-06, 2.87954D-06, 2.89878D-06, 3.06891D-06,
     $     3.14263D-06, 3.28866D-06, 3.44360D-06, 3.55689D-06,
     $     3.90467D-06, 3.99813D-06, 4.12030D-06, 4.35678D-06,
     $     4.60023D-06, 4.68723D-06, 4.84047D-06, 5.00728D-06,
     $     5.26854D-06, 5.40171D-06, 5.72636D-06, 5.85301D-06,
     $     5.95566D-06, 6.27145D-06, 6.37573D-06, 6.62640D-06,
     $     6.79596D-06, 6.92454D-06, 7.15602D-06, 7.29577D-06,
     $     7.56740D-06, 7.79793D-06, 7.93489D-06, 8.44194D-06,
     $     8.64230D-06, 8.78723D-06, 9.06289D-06, 9.27913D-06,
     $     9.65137D-06, 9.90934D-06, 1.04660D-05, 1.08338D-05,
     $     1.09923D-05, 1.12729D-05, 1.13850D-05, 1.18421D-05,
     $     1.20561D-05, 1.24475D-05, 1.28294D-05, 1.31559D-05,
     $     1.35641D-05, 1.39966D-05, 1.40988D-05, 1.50592D-05,
     $     1.51110D-05, 1.61768D-05, 1.66743D-05, 1.70544D-05,
     $     1.72825D-05, 1.78864D-05, 1.79746D-05, 1.87923D-05,
     $     1.91674D-05, 1.96001D-05, 1.97328D-05, 2.09282D-05,
     $     2.11512D-05, 2.21374D-05, 2.27544D-05, 2.29000D-05,
     $     2.37076D-05, 2.42899D-05, 2.45004D-05, 2.50644D-05,
     $     2.54960D-05, 2.65478D-05, 2.72747D-05, 2.77329D-05,
     $     2.81780D-05, 2.86250D-05, 2.88383D-05, 2.98065D-05,
     $     3.03937D-05, 3.08340D-05, 3.10875D-05, 3.21418D-05,
     $     3.24840D-05, 3.37903D-05, 3.39995D-05, 3.45787D-05,
     $     3.49324D-05, 3.55928D-05, 3.57910D-05, 3.68299D-05,
     $     3.77166D-05, 3.85004D-05, 3.94781D-05, 3.95772D-05,
     $     4.00136D-05, 4.06532D-05, 4.19003D-05, 4.25461D-05,
     $     4.32365D-05, 4.37004D-05, 4.48798D-05, 4.52715D-05,
     $     4.63497D-05, 4.70885D-05, 4.75610D-05, 4.77890D-05,
     $     4.84555D-05, 4.90186D-05, 5.05582D-05, 5.09035D-05,
     $     5.09532D-05, 5.24395D-05, 5.26444D-05, 5.32605D-05,
     $     5.40425D-05, 5.64461D-05, 5.68282D-05, 5.69105D-05,
     $     5.87545D-05, 5.94228D-05, 6.05182D-05, 6.09100D-05,
     $     6.11084D-05, 6.18499D-05, 6.32881D-05, 6.42769D-05,
     $     6.53093D-05, 6.57266D-05, 6.72913D-05, 6.79887D-05,
     $     6.92213D-05, 6.95650D-05, 7.11417D-05, 7.22633D-05,
     $     7.42655D-05, 7.49478D-05, 7.60287D-05, 7.65883D-05,
     $     7.73507D-05, 7.80165D-05, 7.91521D-05, 7.92864D-05,
     $     8.07246D-05, 8.17727D-05, 8.27237D-05, 8.34010D-05,
     $     8.66039D-05, 8.71925D-05, 8.79288D-05, 8.91204D-05,
     $     9.02192D-05, 9.05412D-05, 9.19485D-05, 9.24215D-05,
     $     9.43996D-05, 9.50281D-05, 9.67079D-05, 9.73109D-05,
     $     9.79948D-05, 9.85682D-05, 9.97927D-05, 1.00256D-04,
     $     1.01452D-04, 1.03009D-04, 1.03364D-04, 1.05109D-04,
     $     1.05946D-04, 1.06873D-04, 1.07601D-04, 1.09837D-04,
     $     1.11079D-04, 1.11539D-04, 1.12451D-04, 1.14678D-04,
     $     1.16823D-04, 1.17160D-04, 1.17797D-04, 1.18688D-04,
     $     1.19756D-04, 1.20613D-04, 1.21982D-04, 1.22881D-04,
     $     1.25096D-04, 1.27668D-04, 1.27922D-04, 1.29613D-04,
     $     1.30670D-04, 1.31007D-04, 1.32772D-04, 1.33942D-04,
     $     1.35066D-04, 1.35744D-04, 1.37187D-04, 1.37613D-04,
     $     1.38038D-04, 1.40564D-04, 1.41560D-04, 1.42466D-04,
     $     1.43623D-04, 1.44884D-04, 1.46073D-04, 1.47088D-04,
     $     1.47756D-04, 1.49424D-04, 1.50513D-04, 1.51784D-04,
     $     1.52334D-04, 1.54734D-04, 1.55893D-04, 1.56784D-04,
     $     1.57860D-04, 1.58515D-04, 1.61376D-04, 1.62991D-04,
     $     1.63464D-04, 1.64827D-04, 1.65254D-04, 1.67015D-04,
     $     1.68505D-04, 1.69889D-04, 1.70464D-04, 1.72911D-04,
     $     1.73346D-04, 1.73922D-04, 1.74700D-04, 1.76551D-04,
     $     1.78816D-04, 1.80208D-04, 1.81053D-04, 1.81650D-04,
     $     1.82495D-04, 1.84530D-04, 1.86130D-04, 1.86583D-04,
     $     1.87954D-04, 1.88764D-04, 1.89732D-04, 1.91250D-04,
     $     1.93682D-04, 1.94247D-04, 1.94813D-04, 1.96303D-04,
     $     1.97523D-04, 1.99162D-04, 2.00502D-04, 2.01901D-04,
     $     2.02673D-04, 2.03601D-04, 2.05672D-04, 2.06906D-04,
     $     2.07164D-04, 2.09286D-04, 2.10324D-04, 2.11402D-04,
     $     2.12387D-04, 2.14259D-04, 2.16281D-04, 2.17273D-04,
     $     2.17530D-04, 2.19910D-04, 2.21028D-04, 2.22778D-04,
     $     2.24600D-04, 2.25415D-04, 2.25962D-04, 2.28958D-04,
     $     2.29339D-04, 2.30650D-04, 2.31865D-04, 2.33486D-04,
     $     2.34885D-04, 2.35609D-04, 2.36609D-04, 2.38787D-04,
     $     2.41542D-04, 2.42671D-04, 2.44939D-04, 2.45506D-04,
     $     2.47819D-04, 2.48981D-04, 2.49821D-04, 2.50943D-04,
     $     2.51931D-04, 2.52847D-04, 2.53780D-04, 2.54385D-04,
     $     2.55703D-04, 2.57276D-04, 2.58407D-04, 2.59212D-04,
     $     2.61294D-04, 2.61612D-04, 2.64189D-04, 2.65432D-04,
     $     2.66868D-04, 2.69167D-04, 2.69626D-04, 2.71000D-04,
     $     2.72047D-04 /

!     Branching ratio of 5008
      data( br(10,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 6.43509D-34,
     $     6.41498D-12, 5.52041D-11, 1.54674D-10, 4.09559D-10,
     $     9.43335D-10, 1.60711D-09, 2.31166D-09, 2.94264D-09,
     $     1.78745D-08, 1.92496D-08, 2.50910D-08, 3.11251D-08,
     $     3.88520D-08, 4.27835D-08, 5.32188D-08, 7.84003D-08,
     $     8.85248D-08, 2.73533D-07, 3.06706D-07, 4.17447D-07,
     $     5.29822D-07, 5.94855D-07, 7.49375D-07, 1.56686D-06,
     $     2.40948D-06, 3.39495D-06, 5.58391D-06, 6.76423D-06,
     $     1.07287D-05, 1.38118D-05, 1.71597D-05, 1.88778D-05,
     $     2.36714D-05, 2.56170D-05, 2.93090D-05, 3.48405D-05,
     $     3.79774D-05, 4.68972D-05, 5.19855D-05, 6.42010D-05,
     $     7.06289D-05, 7.75491D-05, 8.14313D-05, 9.12883D-05,
     $     1.07363D-04, 1.10391D-04, 1.19657D-04, 1.37782D-04,
     $     1.56514D-04, 1.65341D-04, 1.81981D-04, 1.96771D-04,
     $     2.11083D-04, 2.27715D-04, 2.42762D-04, 2.54648D-04,
     $     2.66478D-04, 2.82218D-04, 2.89003D-04, 3.04858D-04,
     $     3.19501D-04, 3.28559D-04, 3.52119D-04, 3.73629D-04,
     $     3.74434D-04, 4.06416D-04, 4.22112D-04, 4.29850D-04,
     $     4.53544D-04, 4.64042D-04, 4.77093D-04, 4.95011D-04,
     $     5.23892D-04, 5.34283D-04, 5.66444D-04, 5.73681D-04,
     $     5.96020D-04, 6.09569D-04, 6.16445D-04, 6.36313D-04,
     $     6.48431D-04, 6.56287D-04, 6.72187D-04, 6.94286D-04,
     $     7.14947D-04, 7.18744D-04, 7.31139D-04, 7.52251D-04,
     $     7.59553D-04, 7.97533D-04, 8.00098D-04, 8.20223D-04,
     $     8.29118D-04, 8.44910D-04, 8.59871D-04, 8.71249D-04,
     $     9.06207D-04, 9.14373D-04, 9.24712D-04, 9.46997D-04,
     $     9.69681D-04, 9.77344D-04, 9.90735D-04, 1.00430D-03,
     $     1.02633D-03, 1.03605D-03, 1.05900D-03, 1.06713D-03,
     $     1.07417D-03, 1.09320D-03, 1.10029D-03, 1.11493D-03,
     $     1.12434D-03, 1.13134D-03, 1.14441D-03, 1.15391D-03,
     $     1.17040D-03, 1.18281D-03, 1.18967D-03, 1.21590D-03,
     $     1.22551D-03, 1.23282D-03, 1.24604D-03, 1.25530D-03,
     $     1.26963D-03, 1.27992D-03, 1.29994D-03, 1.31226D-03,
     $     1.31776D-03, 1.32729D-03, 1.33162D-03, 1.34671D-03,
     $     1.35397D-03, 1.36729D-03, 1.37926D-03, 1.39062D-03,
     $     1.40317D-03, 1.41579D-03, 1.41924D-03, 1.44529D-03,
     $     1.44662D-03, 1.47283D-03, 1.48418D-03, 1.49261D-03,
     $     1.49783D-03, 1.51066D-03, 1.51242D-03, 1.52912D-03,
     $     1.53725D-03, 1.54648D-03, 1.54915D-03, 1.57081D-03,
     $     1.57459D-03, 1.59124D-03, 1.60212D-03, 1.60448D-03,
     $     1.61743D-03, 1.62652D-03, 1.62983D-03, 1.63816D-03,
     $     1.64449D-03, 1.66013D-03, 1.66989D-03, 1.67594D-03,
     $     1.68175D-03, 1.68764D-03, 1.69075D-03, 1.70356D-03,
     $     1.71144D-03, 1.71691D-03, 1.72023D-03, 1.73324D-03,
     $     1.73725D-03, 1.75260D-03, 1.75494D-03, 1.76127D-03,
     $     1.76508D-03, 1.77207D-03, 1.77414D-03, 1.78505D-03,
     $     1.79527D-03, 1.80330D-03, 1.81305D-03, 1.81401D-03,
     $     1.81833D-03, 1.82510D-03, 1.83737D-03, 1.84343D-03,
     $     1.85013D-03, 1.85456D-03, 1.86507D-03, 1.86841D-03,
     $     1.87765D-03, 1.88380D-03, 1.88790D-03, 1.88978D-03,
     $     1.89525D-03, 1.89958D-03, 1.91083D-03, 1.91325D-03,
     $     1.91359D-03, 1.92376D-03, 1.92512D-03, 1.92928D-03,
     $     1.93454D-03, 1.95056D-03, 1.95294D-03, 1.95344D-03,
     $     1.96488D-03, 1.96882D-03, 1.97513D-03, 1.97758D-03,
     $     1.97872D-03, 1.98278D-03, 1.99041D-03, 1.99569D-03,
     $     2.00077D-03, 2.00274D-03, 2.01023D-03, 2.01350D-03,
     $     2.01892D-03, 2.02035D-03, 2.02680D-03, 2.03129D-03,
     $     2.03916D-03, 2.04183D-03, 2.04584D-03, 2.04786D-03,
     $     2.05052D-03, 2.05285D-03, 2.05687D-03, 2.05732D-03,
     $     2.06258D-03, 2.06622D-03, 2.06944D-03, 2.07176D-03,
     $     2.08251D-03, 2.08455D-03, 2.08695D-03, 2.09129D-03,
     $     2.09499D-03, 2.09602D-03, 2.10040D-03, 2.10184D-03,
     $     2.10784D-03, 2.10973D-03, 2.11464D-03, 2.11637D-03,
     $     2.11831D-03, 2.11990D-03, 2.12338D-03, 2.12466D-03,
     $     2.12826D-03, 2.13291D-03, 2.13391D-03, 2.13870D-03,
     $     2.14101D-03, 2.14347D-03, 2.14548D-03, 2.15185D-03,
     $     2.15506D-03, 2.15621D-03, 2.15861D-03, 2.16439D-03,
     $     2.17013D-03, 2.17102D-03, 2.17261D-03, 2.17485D-03,
     $     2.17741D-03, 2.17942D-03, 2.18262D-03, 2.18469D-03,
     $     2.18996D-03, 2.19605D-03, 2.19662D-03, 2.20056D-03,
     $     2.20298D-03, 2.20374D-03, 2.20791D-03, 2.21079D-03,
     $     2.21364D-03, 2.21516D-03, 2.21830D-03, 2.21920D-03,
     $     2.22015D-03, 2.22555D-03, 2.22753D-03, 2.22936D-03,
     $     2.23192D-03, 2.23460D-03, 2.23694D-03, 2.23889D-03,
     $     2.24012D-03, 2.24302D-03, 2.24489D-03, 2.24725D-03,
     $     2.24819D-03, 2.25211D-03, 2.25395D-03, 2.25528D-03,
     $     2.25709D-03, 2.25812D-03, 2.26298D-03, 2.26538D-03,
     $     2.26609D-03, 2.26847D-03, 2.26916D-03, 2.27181D-03,
     $     2.27383D-03, 2.27585D-03, 2.27675D-03, 2.28037D-03,
     $     2.28099D-03, 2.28183D-03, 2.28287D-03, 2.28533D-03,
     $     2.28819D-03, 2.28984D-03, 2.29086D-03, 2.29158D-03,
     $     2.29276D-03, 2.29565D-03, 2.29783D-03, 2.29842D-03,
     $     2.30006D-03, 2.30104D-03, 2.30215D-03, 2.30409D-03,
     $     2.30731D-03, 2.30795D-03, 2.30857D-03, 2.31021D-03,
     $     2.31182D-03, 2.31380D-03, 2.31530D-03, 2.31671D-03,
     $     2.31762D-03, 2.31863D-03, 2.32074D-03, 2.32213D-03,
     $     2.32238D-03, 2.32458D-03, 2.32558D-03, 2.32678D-03,
     $     2.32779D-03, 2.32972D-03, 2.33189D-03, 2.33289D-03,
     $     2.33317D-03, 2.33524D-03, 2.33617D-03, 2.33792D-03,
     $     2.33971D-03, 2.34038D-03, 2.34080D-03, 2.34350D-03,
     $     2.34394D-03, 2.34533D-03, 2.34638D-03, 2.34761D-03,
     $     2.34864D-03, 2.34930D-03, 2.35019D-03, 2.35191D-03,
     $     2.35387D-03, 2.35461D-03, 2.35601D-03, 2.35635D-03,
     $     2.35799D-03, 2.35890D-03, 2.35946D-03, 2.36017D-03,
     $     2.36108D-03, 2.36210D-03, 2.36286D-03, 2.36330D-03,
     $     2.36424D-03, 2.36541D-03, 2.36643D-03, 2.36718D-03,
     $     2.36886D-03, 2.36911D-03, 2.37126D-03, 2.37235D-03,
     $     2.37368D-03, 2.37551D-03, 2.37584D-03, 2.37681D-03,
     $     2.37783D-03 /

!     Branching ratio of 4008
      data( br(11,k), k=1, 500) /
     $     8.52019D-06, 6.71249D-06, 1.37061D-05,
     $     1.65824D-05, 2.16001D-05, 8.19258D-02, 8.10113D-02,
     $     7.07902D-02, 7.59394D-02, 8.00483D-02, 8.84541D-02,
     $     8.77001D-02, 8.27734D-02, 8.06466D-02, 8.04965D-02,
     $     8.21723D-02, 8.46290D-02, 9.43511D-02, 9.97125D-02,
     $     1.01377D-01, 1.02462D-01, 1.03869D-01, 1.02890D-01,
     $     1.02021D-01, 1.01392D-01, 1.00192D-01, 9.91872D-02,
     $     9.90315D-02, 9.89889D-02, 9.83624D-02, 9.74104D-02,
     $     9.65802D-02, 9.63669D-02, 9.57855D-02, 9.49819D-02,
     $     9.47909D-02, 9.49129D-02, 9.53071D-02, 9.56814D-02,
     $     9.74799D-02, 9.77157D-02, 9.76740D-02, 9.70327D-02,
     $     9.69779D-02, 9.68773D-02, 9.67059D-02, 9.64865D-02,
     $     9.64127D-02, 9.61838D-02, 9.61396D-02, 9.61324D-02,
     $     9.60716D-02, 9.61194D-02, 9.60689D-02, 9.59628D-02,
     $     9.59479D-02, 9.59801D-02, 9.58990D-02, 9.58419D-02,
     $     9.59552D-02, 9.60745D-02, 9.61847D-02, 9.64361D-02,
     $     9.69650D-02, 9.72913D-02, 9.79517D-02, 9.81011D-02,
     $     9.85880D-02, 9.88763D-02, 9.90503D-02, 9.91569D-02,
     $     9.93223D-02, 9.94093D-02, 9.94959D-02, 9.95373D-02,
     $     9.97778D-02, 9.98199D-02, 9.98806D-02, 9.99395D-02,
     $     9.99603D-02, 1.00029D-01, 1.00056D-01, 1.00045D-01,
     $     1.00040D-01, 1.00034D-01, 1.00045D-01, 1.00085D-01,
     $     1.00119D-01, 1.00143D-01, 1.00201D-01, 1.00294D-01,
     $     1.00521D-01, 1.00544D-01, 1.00601D-01, 1.00636D-01,
     $     1.00641D-01, 1.00634D-01, 1.00648D-01, 1.00677D-01,
     $     1.00673D-01, 1.00667D-01, 1.00722D-01, 1.00806D-01,
     $     1.00829D-01, 1.00878D-01, 1.00927D-01, 1.00945D-01,
     $     1.00943D-01, 1.00941D-01, 1.00976D-01, 1.00994D-01,
     $     1.01001D-01, 1.01009D-01, 1.01016D-01, 1.01024D-01,
     $     1.01035D-01, 1.01042D-01, 1.01048D-01, 1.01052D-01,
     $     1.01214D-01, 1.01236D-01, 1.01264D-01, 1.01311D-01,
     $     1.01331D-01, 1.01337D-01, 1.01420D-01, 1.01584D-01,
     $     1.01619D-01, 1.01769D-01, 1.01783D-01, 1.01813D-01,
     $     1.01828D-01, 1.01834D-01, 1.01843D-01, 1.01888D-01,
     $     1.01934D-01, 1.01968D-01, 1.02009D-01, 1.02033D-01,
     $     1.02113D-01, 1.02157D-01, 1.02197D-01, 1.02213D-01,
     $     1.02255D-01, 1.02265D-01, 1.02290D-01, 1.02326D-01,
     $     1.02354D-01, 1.02386D-01, 1.02444D-01, 1.02503D-01,
     $     1.02511D-01, 1.02529D-01, 1.02532D-01, 1.02545D-01,
     $     1.02565D-01, 1.02564D-01, 1.02558D-01, 1.02594D-01,
     $     1.02608D-01, 1.02614D-01, 1.02600D-01, 1.02615D-01,
     $     1.02642D-01, 1.02683D-01, 1.02687D-01, 1.02689D-01,
     $     1.02680D-01, 1.02686D-01, 1.02692D-01, 1.02705D-01,
     $     1.02707D-01, 1.02700D-01, 1.02687D-01, 1.02699D-01,
     $     1.02699D-01, 1.02724D-01, 1.02715D-01, 1.02705D-01,
     $     1.02680D-01, 1.02664D-01, 1.02659D-01, 1.02657D-01,
     $     1.02660D-01, 1.02654D-01, 1.02631D-01, 1.02633D-01,
     $     1.02617D-01, 1.02595D-01, 1.02582D-01, 1.02576D-01,
     $     1.02561D-01, 1.02549D-01, 1.02528D-01, 1.02518D-01,
     $     1.02491D-01, 1.02485D-01, 1.02461D-01, 1.02437D-01,
     $     1.02424D-01, 1.02392D-01, 1.02392D-01, 1.02371D-01,
     $     1.02366D-01, 1.02345D-01, 1.02314D-01, 1.02301D-01,
     $     1.02281D-01, 1.02267D-01, 1.02246D-01, 1.02227D-01,
     $     1.02209D-01, 1.02199D-01, 1.02184D-01, 1.02162D-01,
     $     1.02142D-01, 1.02120D-01, 1.02072D-01, 1.02050D-01,
     $     1.02037D-01, 1.01984D-01, 1.01976D-01, 1.01939D-01,
     $     1.01911D-01, 1.01890D-01, 1.01858D-01, 1.01856D-01,
     $     1.01840D-01, 1.01815D-01, 1.01797D-01, 1.01754D-01,
     $     1.01737D-01, 1.01731D-01, 1.01719D-01, 1.01703D-01,
     $     1.01664D-01, 1.01644D-01, 1.01589D-01, 1.01548D-01,
     $     1.01532D-01, 1.01505D-01, 1.01503D-01, 1.01458D-01,
     $     1.01438D-01, 1.01404D-01, 1.01363D-01, 1.01342D-01,
     $     1.01306D-01, 1.01265D-01, 1.01262D-01, 1.01182D-01,
     $     1.01178D-01, 1.01102D-01, 1.01067D-01, 1.01043D-01,
     $     1.01032D-01, 1.00993D-01, 1.00987D-01, 1.00935D-01,
     $     1.00917D-01, 1.00903D-01, 1.00896D-01, 1.00816D-01,
     $     1.00800D-01, 1.00732D-01, 1.00698D-01, 1.00688D-01,
     $     1.00638D-01, 1.00602D-01, 1.00590D-01, 1.00552D-01,
     $     1.00526D-01, 1.00466D-01, 1.00418D-01, 1.00388D-01,
     $     1.00360D-01, 1.00332D-01, 1.00325D-01, 1.00271D-01,
     $     1.00239D-01, 1.00210D-01, 1.00198D-01, 1.00132D-01,
     $     1.00108D-01, 1.00021D-01, 1.00006D-01, 9.99637D-02,
     $     9.99374D-02, 9.98877D-02, 9.98727D-02, 9.97967D-02,
     $     9.97483D-02, 9.96956D-02, 9.96289D-02, 9.96220D-02,
     $     9.95937D-02, 9.95596D-02, 9.94838D-02, 9.94430D-02,
     $     9.94023D-02, 9.93753D-02, 9.92994D-02, 9.92733D-02,
     $     9.92051D-02, 9.91590D-02, 9.91337D-02, 9.91204D-02,
     $     9.90821D-02, 9.90478D-02, 9.89531D-02, 9.89321D-02,
     $     9.89291D-02, 9.88415D-02, 9.88295D-02, 9.87950D-02,
     $     9.87530D-02, 9.86379D-02, 9.86190D-02, 9.86149D-02,
     $     9.85258D-02, 9.84933D-02, 9.84399D-02, 9.84241D-02,
     $     9.84147D-02, 9.83779D-02, 9.83061D-02, 9.82606D-02,
     $     9.82097D-02, 9.81887D-02, 9.81146D-02, 9.80833D-02,
     $     9.80246D-02, 9.80078D-02, 9.79321D-02, 9.78797D-02,
     $     9.77896D-02, 9.77601D-02, 9.77114D-02, 9.76855D-02,
     $     9.76492D-02, 9.76182D-02, 9.75665D-02, 9.75600D-02,
     $     9.74953D-02, 9.74459D-02, 9.74000D-02, 9.73687D-02,
     $     9.72176D-02, 9.71910D-02, 9.71560D-02, 9.71043D-02,
     $     9.70536D-02, 9.70383D-02, 9.69706D-02, 9.69477D-02,
     $     9.68543D-02, 9.68248D-02, 9.67459D-02, 9.67178D-02,
     $     9.66856D-02, 9.66586D-02, 9.66019D-02, 9.65803D-02,
     $     9.65280D-02, 9.64617D-02, 9.64461D-02, 9.63689D-02,
     $     9.63323D-02, 9.62910D-02, 9.62593D-02, 9.61654D-02,
     $     9.61108D-02, 9.60905D-02, 9.60524D-02, 9.59590D-02,
     $     9.58720D-02, 9.58582D-02, 9.58313D-02, 9.57940D-02,
     $     9.57479D-02, 9.57106D-02, 9.56517D-02, 9.56126D-02,
     $     9.55170D-02, 9.54066D-02, 9.53953D-02, 9.53215D-02,
     $     9.52753D-02, 9.52607D-02, 9.51855D-02, 9.51361D-02,
     $     9.50911D-02, 9.50618D-02, 9.49987D-02, 9.49798D-02,
     $     9.49615D-02, 9.48512D-02, 9.48070D-02, 9.47673D-02,
     $     9.47192D-02, 9.46669D-02, 9.46170D-02, 9.45744D-02,
     $     9.45460D-02, 9.44747D-02, 9.44284D-02, 9.43765D-02,
     $     9.43535D-02, 9.42537D-02, 9.42066D-02, 9.41698D-02,
     $     9.41278D-02, 9.41018D-02, 9.39915D-02, 9.39277D-02,
     $     9.39094D-02, 9.38592D-02, 9.38431D-02, 9.37755D-02,
     $     9.37173D-02, 9.36647D-02, 9.36432D-02, 9.35511D-02,
     $     9.35346D-02, 9.35133D-02, 9.34836D-02, 9.34127D-02,
     $     9.33255D-02, 9.32716D-02, 9.32390D-02, 9.32161D-02,
     $     9.31846D-02, 9.31104D-02, 9.30517D-02, 9.30348D-02,
     $     9.29834D-02, 9.29534D-02, 9.29172D-02, 9.28620D-02,
     $     9.27748D-02, 9.27539D-02, 9.27329D-02, 9.26778D-02,
     $     9.26348D-02, 9.25757D-02, 9.25270D-02, 9.24753D-02,
     $     9.24479D-02, 9.24145D-02, 9.23391D-02, 9.22948D-02,
     $     9.22854D-02, 9.22091D-02, 9.21717D-02, 9.21336D-02,
     $     9.20983D-02, 9.20322D-02, 9.19615D-02, 9.19265D-02,
     $     9.19176D-02, 9.18326D-02, 9.17926D-02, 9.17329D-02,
     $     9.16704D-02, 9.16420D-02, 9.16228D-02, 9.15206D-02,
     $     9.15084D-02, 9.14654D-02, 9.14242D-02, 9.13689D-02,
     $     9.13214D-02, 9.12975D-02, 9.12641D-02, 9.11904D-02,
     $     9.10972D-02, 9.10587D-02, 9.09811D-02, 9.09617D-02,
     $     9.08834D-02, 9.08443D-02, 9.08155D-02, 9.07769D-02,
     $     9.07444D-02, 9.07144D-02, 9.06827D-02, 9.06619D-02,
     $     9.06165D-02, 9.05619D-02, 9.05231D-02, 9.04957D-02,
     $     9.04226D-02, 9.04113D-02, 9.03194D-02, 9.02750D-02,
     $     9.02238D-02, 9.01387D-02, 9.01215D-02, 9.00695D-02,
     $     9.00308D-02 /

!     Branching ratio of 3008
      data( br(12,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 3.38122D-38,
     $     3.54433D-36, 7.58581D-11, 3.94965D-10, 1.35793D-09,
     $     1.78312D-09, 2.47102D-09, 6.98602D-09, 2.24550D-08,
     $     8.18515D-08, 1.47876D-07, 3.45607D-07, 5.57622D-07,
     $     2.40217D-06, 4.29126D-06, 6.76861D-06, 8.27586D-06,
     $     1.25077D-05, 1.45793D-05, 1.64015D-05, 2.22546D-05,
     $     2.64458D-05, 3.09118D-05, 3.82955D-05, 4.27469D-05,
     $     4.70681D-05, 5.37988D-05, 6.18253D-05, 8.24917D-05,
     $     9.27685D-05, 9.72547D-05, 1.13549D-04, 1.31813D-04,
     $     1.50164D-04, 1.73000D-04, 1.90226D-04, 2.10567D-04,
     $     2.35369D-04, 2.54404D-04, 2.69459D-04, 2.79801D-04,
     $     3.93961D-04, 4.00868D-04, 4.28306D-04, 4.54486D-04,
     $     4.86071D-04, 5.01944D-04, 5.43241D-04, 6.13149D-04,
     $     6.40666D-04, 7.76923D-04, 7.98779D-04, 8.76691D-04,
     $     9.38025D-04, 9.64553D-04, 1.01243D-03, 1.09816D-03,
     $     1.13660D-03, 1.20160D-03, 1.31903D-03, 1.34568D-03,
     $     1.42690D-03, 1.50988D-03, 1.58190D-03, 1.62416D-03,
     $     1.71775D-03, 1.76374D-03, 1.82910D-03, 1.89406D-03,
     $     1.92003D-03, 2.05355D-03, 2.08509D-03, 2.21517D-03,
     $     2.29084D-03, 2.35392D-03, 2.39274D-03, 2.48208D-03,
     $     2.61219D-03, 2.63781D-03, 2.71823D-03, 2.82756D-03,
     $     2.96765D-03, 3.03084D-03, 3.16577D-03, 3.26502D-03,
     $     3.34104D-03, 3.42651D-03, 3.52006D-03, 3.58989D-03,
     $     3.66412D-03, 3.74789D-03, 3.78285D-03, 3.86021D-03,
     $     3.93359D-03, 3.98326D-03, 4.10988D-03, 4.20387D-03,
     $     4.20782D-03, 4.34576D-03, 4.42239D-03, 4.46145D-03,
     $     4.57503D-03, 4.62674D-03, 4.68391D-03, 4.75880D-03,
     $     4.87304D-03, 4.91588D-03, 5.04764D-03, 5.07444D-03,
     $     5.16252D-03, 5.22022D-03, 5.24957D-03, 5.31860D-03,
     $     5.36490D-03, 5.39587D-03, 5.45610D-03, 5.53202D-03,
     $     5.60508D-03, 5.61875D-03, 5.66436D-03, 5.73498D-03,
     $     5.76043D-03, 5.87972D-03, 5.88707D-03, 5.95033D-03,
     $     5.97627D-03, 6.02638D-03, 6.07757D-03, 6.11268D-03,
     $     6.21410D-03, 6.24084D-03, 6.27548D-03, 6.33904D-03,
     $     6.40242D-03, 6.42468D-03, 6.46301D-03, 6.50446D-03,
     $     6.56635D-03, 6.59763D-03, 6.67011D-03, 6.69759D-03,
     $     6.71867D-03, 6.78297D-03, 6.80305D-03, 6.85153D-03,
     $     6.88385D-03, 6.90798D-03, 6.94950D-03, 6.97263D-03,
     $     7.01785D-03, 7.05587D-03, 7.07808D-03, 7.15342D-03,
     $     7.18111D-03, 7.19960D-03, 7.23367D-03, 7.25998D-03,
     $     7.30383D-03, 7.33223D-03, 7.39078D-03, 7.42818D-03,
     $     7.44348D-03, 7.47012D-03, 7.48012D-03, 7.52281D-03,
     $     7.54189D-03, 7.57653D-03, 7.61025D-03, 7.63659D-03,
     $     7.66930D-03, 7.70254D-03, 7.70946D-03, 7.77694D-03,
     $     7.78038D-03, 7.84622D-03, 7.87484D-03, 7.89601D-03,
     $     7.90803D-03, 7.93984D-03, 7.94450D-03, 7.98622D-03,
     $     8.00380D-03, 8.02317D-03, 8.02917D-03, 8.08495D-03,
     $     8.09523D-03, 8.13911D-03, 8.16463D-03, 8.17080D-03,
     $     8.20401D-03, 8.22727D-03, 8.23541D-03, 8.25791D-03,
     $     8.27453D-03, 8.31334D-03, 8.34066D-03, 8.35744D-03,
     $     8.37351D-03, 8.38918D-03, 8.39589D-03, 8.42849D-03,
     $     8.44796D-03, 8.46312D-03, 8.47148D-03, 8.50679D-03,
     $     8.51842D-03, 8.56090D-03, 8.56778D-03, 8.58682D-03,
     $     8.59838D-03, 8.61974D-03, 8.62608D-03, 8.65780D-03,
     $     8.68156D-03, 8.70355D-03, 8.73000D-03, 8.73269D-03,
     $     8.74392D-03, 8.75943D-03, 8.78926D-03, 8.80456D-03,
     $     8.81997D-03, 8.83023D-03, 8.85718D-03, 8.86615D-03,
     $     8.88937D-03, 8.90480D-03, 8.91395D-03, 8.91854D-03,
     $     8.93140D-03, 8.94266D-03, 8.97303D-03, 8.97970D-03,
     $     8.98065D-03, 9.00802D-03, 9.01174D-03, 9.02234D-03,
     $     9.03564D-03, 9.07378D-03, 9.07980D-03, 9.08112D-03,
     $     9.10887D-03, 9.11890D-03, 9.13506D-03, 9.14027D-03,
     $     9.14318D-03, 9.15441D-03, 9.17585D-03, 9.18999D-03,
     $     9.20538D-03, 9.21169D-03, 9.23423D-03, 9.24367D-03,
     $     9.26109D-03, 9.26610D-03, 9.28844D-03, 9.30423D-03,
     $     9.33128D-03, 9.34009D-03, 9.35446D-03, 9.36203D-03,
     $     9.37259D-03, 9.38173D-03, 9.39714D-03, 9.39901D-03,
     $     9.41774D-03, 9.43166D-03, 9.44443D-03, 9.45309D-03,
     $     9.49459D-03, 9.50171D-03, 9.51114D-03, 9.52548D-03,
     $     9.53889D-03, 9.54293D-03, 9.56087D-03, 9.56692D-03,
     $     9.59126D-03, 9.59886D-03, 9.61932D-03, 9.62658D-03,
     $     9.63491D-03, 9.64190D-03, 9.65632D-03, 9.66185D-03,
     $     9.67523D-03, 9.69213D-03, 9.69616D-03, 9.71581D-03,
     $     9.72490D-03, 9.73529D-03, 9.74332D-03, 9.76686D-03,
     $     9.78058D-03, 9.78568D-03, 9.79491D-03, 9.81825D-03,
     $     9.83944D-03, 9.84263D-03, 9.84909D-03, 9.85782D-03,
     $     9.86887D-03, 9.87783D-03, 9.89183D-03, 9.90114D-03,
     $     9.92387D-03, 9.94884D-03, 9.95142D-03, 9.96839D-03,
     $     9.97892D-03, 9.98217D-03, 9.99917D-03, 1.00100D-02,
     $     1.00198D-02, 1.00263D-02, 1.00401D-02, 1.00442D-02,
     $     1.00482D-02, 1.00722D-02, 1.00818D-02, 1.00901D-02,
     $     1.01000D-02, 1.01113D-02, 1.01220D-02, 1.01310D-02,
     $     1.01371D-02, 1.01527D-02, 1.01627D-02, 1.01739D-02,
     $     1.01790D-02, 1.02007D-02, 1.02110D-02, 1.02191D-02,
     $     1.02286D-02, 1.02343D-02, 1.02587D-02, 1.02728D-02,
     $     1.02768D-02, 1.02876D-02, 1.02912D-02, 1.03059D-02,
     $     1.03188D-02, 1.03300D-02, 1.03347D-02, 1.03544D-02,
     $     1.03580D-02, 1.03625D-02, 1.03689D-02, 1.03841D-02,
     $     1.04028D-02, 1.04146D-02, 1.04216D-02, 1.04265D-02,
     $     1.04333D-02, 1.04489D-02, 1.04615D-02, 1.04651D-02,
     $     1.04763D-02, 1.04827D-02, 1.04905D-02, 1.05021D-02,
     $     1.05204D-02, 1.05249D-02, 1.05294D-02, 1.05412D-02,
     $     1.05501D-02, 1.05626D-02, 1.05731D-02, 1.05844D-02,
     $     1.05904D-02, 1.05977D-02, 1.06142D-02, 1.06236D-02,
     $     1.06257D-02, 1.06425D-02, 1.06507D-02, 1.06587D-02,
     $     1.06665D-02, 1.06812D-02, 1.06963D-02, 1.07037D-02,
     $     1.07055D-02, 1.07242D-02, 1.07330D-02, 1.07459D-02,
     $     1.07596D-02, 1.07661D-02, 1.07704D-02, 1.07930D-02,
     $     1.07957D-02, 1.08054D-02, 1.08146D-02, 1.08272D-02,
     $     1.08380D-02, 1.08434D-02, 1.08508D-02, 1.08671D-02,
     $     1.08881D-02, 1.08969D-02, 1.09148D-02, 1.09194D-02,
     $     1.09375D-02, 1.09461D-02, 1.09526D-02, 1.09615D-02,
     $     1.09687D-02, 1.09750D-02, 1.09820D-02, 1.09867D-02,
     $     1.09969D-02, 1.10089D-02, 1.10175D-02, 1.10236D-02,
     $     1.10398D-02, 1.10422D-02, 1.10623D-02, 1.10719D-02,
     $     1.10829D-02, 1.11012D-02, 1.11050D-02, 1.11164D-02,
     $     1.11247D-02 /

!     Branching ratio of 5007
      data( br(13,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00 /

!     Branching ratio of 4007
      data( br(14,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 8.85977D-10,
     $     3.64486D-05, 5.94465D-04, 1.77904D-03, 1.73230D-02,
     $     2.66617D-02, 5.47490D-02, 7.25288D-02, 8.87232D-02,
     $     9.60948D-02, 9.83990D-02, 1.02945D-01, 1.03337D-01,
     $     1.03578D-01, 1.04186D-01, 1.04738D-01, 1.03539D-01,
     $     1.02631D-01, 1.02320D-01, 1.01465D-01, 1.00615D-01,
     $     1.00474D-01, 1.00243D-01, 9.95294D-02, 9.83150D-02,
     $     9.70629D-02, 9.64926D-02, 9.53670D-02, 9.33874D-02,
     $     9.12893D-02, 9.04638D-02, 8.90583D-02, 8.80401D-02,
     $     8.42244D-02, 8.22379D-02, 8.14611D-02, 7.79002D-02,
     $     7.76104D-02, 7.70812D-02, 7.51769D-02, 7.43126D-02,
     $     7.39519D-02, 7.24167D-02, 7.20326D-02, 7.02212D-02,
     $     6.95073D-02, 6.80438D-02, 6.75002D-02, 6.65545D-02,
     $     6.55061D-02, 6.49232D-02, 6.38310D-02, 6.27274D-02,
     $     6.13063D-02, 6.08052D-02, 5.99767D-02, 5.87123D-02,
     $     5.75114D-02, 5.67895D-02, 5.51846D-02, 5.47782D-02,
     $     5.38873D-02, 5.29939D-02, 5.25467D-02, 5.20344D-02,
     $     5.14911D-02, 5.10583D-02, 5.06936D-02, 5.03685D-02,
     $     4.92504D-02, 4.86716D-02, 4.81026D-02, 4.79670D-02,
     $     4.78333D-02, 4.68372D-02, 4.64195D-02, 4.60664D-02,
     $     4.59815D-02, 4.58708D-02, 4.54952D-02, 4.50159D-02,
     $     4.44078D-02, 4.41033D-02, 4.36294D-02, 4.33169D-02,
     $     4.22589D-02, 4.17698D-02, 4.13145D-02, 4.10917D-02,
     $     4.05998D-02, 4.04014D-02, 4.02344D-02, 3.97647D-02,
     $     3.94849D-02, 3.92216D-02, 3.88274D-02, 3.85984D-02,
     $     3.84141D-02, 3.81505D-02, 3.78777D-02, 3.73312D-02,
     $     3.71103D-02, 3.70213D-02, 3.67179D-02, 3.64249D-02,
     $     3.61718D-02, 3.58930D-02, 3.57023D-02, 3.54956D-02,
     $     3.52652D-02, 3.51007D-02, 3.49767D-02, 3.48940D-02,
     $     3.40636D-02, 3.40156D-02, 3.38497D-02, 3.36948D-02,
     $     3.35292D-02, 3.34519D-02, 3.32452D-02, 3.29223D-02,
     $     3.28134D-02, 3.23282D-02, 3.22584D-02, 3.20260D-02,
     $     3.18561D-02, 3.17849D-02, 3.16589D-02, 3.14338D-02,
     $     3.13289D-02, 3.11686D-02, 3.08967D-02, 3.08322D-02,
     $     3.06394D-02, 3.04625D-02, 3.03132D-02, 3.02310D-02,
     $     3.00512D-02, 2.99691D-02, 2.98506D-02, 2.97325D-02,
     $     2.96811D-02, 2.94568D-02, 2.93928D-02, 2.91762D-02,
     $     2.90613D-02, 2.89633D-02, 2.89071D-02, 2.87783D-02,
     $     2.85960D-02, 2.85621D-02, 2.84592D-02, 2.83091D-02,
     $     2.81301D-02, 2.80496D-02, 2.78886D-02, 2.77637D-02,
     $     2.76666D-02, 2.75591D-02, 2.74511D-02, 2.73708D-02,
     $     2.72900D-02, 2.71940D-02, 2.71537D-02, 2.70639D-02,
     $     2.69818D-02, 2.69294D-02, 2.67959D-02, 2.66930D-02,
     $     2.66890D-02, 2.65405D-02, 2.64647D-02, 2.64274D-02,
     $     2.63191D-02, 2.62721D-02, 2.62171D-02, 2.61467D-02,
     $     2.60439D-02, 2.60054D-02, 2.58920D-02, 2.58673D-02,
     $     2.57911D-02, 2.57448D-02, 2.57218D-02, 2.56607D-02,
     $     2.56236D-02, 2.55998D-02, 2.55544D-02, 2.54954D-02,
     $     2.54398D-02, 2.54297D-02, 2.53967D-02, 2.53420D-02,
     $     2.53231D-02, 2.52335D-02, 2.52282D-02, 2.51813D-02,
     $     2.51622D-02, 2.51262D-02, 2.50900D-02, 2.50629D-02,
     $     2.49849D-02, 2.49658D-02, 2.49416D-02, 2.48959D-02,
     $     2.48517D-02, 2.48361D-02, 2.48089D-02, 2.47810D-02,
     $     2.47407D-02, 2.47215D-02, 2.46758D-02, 2.46593D-02,
     $     2.46455D-02, 2.46066D-02, 2.45954D-02, 2.45673D-02,
     $     2.45485D-02, 2.45345D-02, 2.45093D-02, 2.44955D-02,
     $     2.44699D-02, 2.44485D-02, 2.44360D-02, 2.43945D-02,
     $     2.43797D-02, 2.43695D-02, 2.43518D-02, 2.43382D-02,
     $     2.43158D-02, 2.43032D-02, 2.42752D-02, 2.42581D-02,
     $     2.42511D-02, 2.42392D-02, 2.42366D-02, 2.42203D-02,
     $     2.42128D-02, 2.42021D-02, 2.41908D-02, 2.41845D-02,
     $     2.41746D-02, 2.41646D-02, 2.41636D-02, 2.41469D-02,
     $     2.41460D-02, 2.41326D-02, 2.41273D-02, 2.41244D-02,
     $     2.41237D-02, 2.41195D-02, 2.41188D-02, 2.41157D-02,
     $     2.41162D-02, 2.41182D-02, 2.41181D-02, 2.41144D-02,
     $     2.41136D-02, 2.41109D-02, 2.41129D-02, 2.41128D-02,
     $     2.41132D-02, 2.41138D-02, 2.41142D-02, 2.41145D-02,
     $     2.41152D-02, 2.41205D-02, 2.41223D-02, 2.41236D-02,
     $     2.41250D-02, 2.41268D-02, 2.41297D-02, 2.41369D-02,
     $     2.41427D-02, 2.41456D-02, 2.41489D-02, 2.41578D-02,
     $     2.41600D-02, 2.41703D-02, 2.41718D-02, 2.41758D-02,
     $     2.41784D-02, 2.41833D-02, 2.41849D-02, 2.41945D-02,
     $     2.42109D-02, 2.42215D-02, 2.42349D-02, 2.42362D-02,
     $     2.42428D-02, 2.42569D-02, 2.42799D-02, 2.42913D-02,
     $     2.43054D-02, 2.43152D-02, 2.43370D-02, 2.43439D-02,
     $     2.43645D-02, 2.43787D-02, 2.43899D-02, 2.43947D-02,
     $     2.44089D-02, 2.44197D-02, 2.44481D-02, 2.44543D-02,
     $     2.44552D-02, 2.44825D-02, 2.44863D-02, 2.44982D-02,
     $     2.45141D-02, 2.45688D-02, 2.45770D-02, 2.45787D-02,
     $     2.46192D-02, 2.46333D-02, 2.46563D-02, 2.46663D-02,
     $     2.46706D-02, 2.46857D-02, 2.47143D-02, 2.47357D-02,
     $     2.47556D-02, 2.47633D-02, 2.47944D-02, 2.48088D-02,
     $     2.48321D-02, 2.48382D-02, 2.48666D-02, 2.48872D-02,
     $     2.49255D-02, 2.49393D-02, 2.49599D-02, 2.49703D-02,
     $     2.49839D-02, 2.49962D-02, 2.50182D-02, 2.50207D-02,
     $     2.50503D-02, 2.50710D-02, 2.50895D-02, 2.51035D-02,
     $     2.51696D-02, 2.51828D-02, 2.51981D-02, 2.52265D-02,
     $     2.52509D-02, 2.52577D-02, 2.52865D-02, 2.52961D-02,
     $     2.53374D-02, 2.53506D-02, 2.53853D-02, 2.53978D-02,
     $     2.54116D-02, 2.54231D-02, 2.54486D-02, 2.54580D-02,
     $     2.54849D-02, 2.55208D-02, 2.55285D-02, 2.55661D-02,
     $     2.55846D-02, 2.56043D-02, 2.56204D-02, 2.56721D-02,
     $     2.56986D-02, 2.57082D-02, 2.57291D-02, 2.57788D-02,
     $     2.58293D-02, 2.58372D-02, 2.58515D-02, 2.58718D-02,
     $     2.58949D-02, 2.59132D-02, 2.59430D-02, 2.59623D-02,
     $     2.60114D-02, 2.60706D-02, 2.60762D-02, 2.61144D-02,
     $     2.61383D-02, 2.61460D-02, 2.61877D-02, 2.62164D-02,
     $     2.62455D-02, 2.62611D-02, 2.62942D-02, 2.63037D-02,
     $     2.63136D-02, 2.63708D-02, 2.63925D-02, 2.64129D-02,
     $     2.64410D-02, 2.64704D-02, 2.64969D-02, 2.65196D-02,
     $     2.65339D-02, 2.65683D-02, 2.65906D-02, 2.66179D-02,
     $     2.66290D-02, 2.66772D-02, 2.67006D-02, 2.67176D-02,
     $     2.67395D-02, 2.67526D-02, 2.68116D-02, 2.68428D-02,
     $     2.68522D-02, 2.68813D-02, 2.68900D-02, 2.69249D-02,
     $     2.69527D-02, 2.69803D-02, 2.69919D-02, 2.70408D-02,
     $     2.70492D-02, 2.70610D-02, 2.70758D-02, 2.71113D-02,
     $     2.71535D-02, 2.71786D-02, 2.71942D-02, 2.72053D-02,
     $     2.72217D-02, 2.72623D-02, 2.72936D-02, 2.73022D-02,
     $     2.73273D-02, 2.73425D-02, 2.73602D-02, 2.73897D-02,
     $     2.74375D-02, 2.74478D-02, 2.74580D-02, 2.74851D-02,
     $     2.75094D-02, 2.75404D-02, 2.75649D-02, 2.75894D-02,
     $     2.76037D-02, 2.76207D-02, 2.76575D-02, 2.76805D-02,
     $     2.76850D-02, 2.77230D-02, 2.77414D-02, 2.77620D-02,
     $     2.77797D-02, 2.78135D-02, 2.78515D-02, 2.78697D-02,
     $     2.78748D-02, 2.79161D-02, 2.79354D-02, 2.79676D-02,
     $     2.80003D-02, 2.80139D-02, 2.80229D-02, 2.80752D-02,
     $     2.80823D-02, 2.81060D-02, 2.81268D-02, 2.81530D-02,
     $     2.81758D-02, 2.81882D-02, 2.82055D-02, 2.82421D-02,
     $     2.82870D-02, 2.83048D-02, 2.83398D-02, 2.83484D-02,
     $     2.83850D-02, 2.84050D-02, 2.84184D-02, 2.84359D-02,
     $     2.84534D-02, 2.84711D-02, 2.84870D-02, 2.84970D-02,
     $     2.85187D-02, 2.85457D-02, 2.85655D-02, 2.85798D-02,
     $     2.86159D-02, 2.86215D-02, 2.86670D-02, 2.86897D-02,
     $     2.87166D-02, 2.87581D-02, 2.87662D-02, 2.87902D-02,
     $     2.88102D-02 /

!     Branching ratio of 3007
      data( br(15,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 3.62135D-04, 1.19131D-03,
     $     1.56993D-02, 5.88361D-02, 7.76387D-02, 1.91164D-01,
     $     2.17748D-01, 2.74816D-01, 2.92699D-01, 3.05768D-01,
     $     3.15538D-01, 3.19240D-01, 3.15339D-01, 3.10095D-01,
     $     3.07668D-01, 3.04827D-01, 2.98338D-01, 2.85631D-01,
     $     2.81240D-01, 2.77428D-01, 2.71602D-01, 2.67089D-01,
     $     2.66395D-01, 2.64158D-01, 2.60364D-01, 2.53522D-01,
     $     2.48295D-01, 2.44907D-01, 2.40296D-01, 2.33656D-01,
     $     2.27026D-01, 2.24476D-01, 2.20219D-01, 2.17190D-01,
     $     2.05583D-01, 1.99907D-01, 1.97691D-01, 1.87504D-01,
     $     1.86692D-01, 1.85221D-01, 1.79706D-01, 1.77406D-01,
     $     1.76454D-01, 1.72422D-01, 1.71416D-01, 1.66405D-01,
     $     1.64523D-01, 1.60546D-01, 1.59127D-01, 1.56681D-01,
     $     1.53903D-01, 1.52294D-01, 1.49470D-01, 1.46685D-01,
     $     1.43075D-01, 1.41744D-01, 1.39642D-01, 1.36443D-01,
     $     1.33360D-01, 1.31545D-01, 1.27662D-01, 1.26692D-01,
     $     1.24470D-01, 1.22319D-01, 1.21189D-01, 1.19951D-01,
     $     1.18603D-01, 1.17537D-01, 1.16621D-01, 1.15826D-01,
     $     1.13019D-01, 1.11622D-01, 1.10228D-01, 1.09878D-01,
     $     1.09548D-01, 1.07106D-01, 1.06103D-01, 1.05257D-01,
     $     1.05054D-01, 1.04789D-01, 1.03892D-01, 1.02717D-01,
     $     1.01241D-01, 1.00505D-01, 9.93494D-02, 9.85832D-02,
     $     9.60242D-02, 9.48792D-02, 9.38255D-02, 9.32957D-02,
     $     9.21545D-02, 9.16980D-02, 9.13143D-02, 9.02359D-02,
     $     8.95932D-02, 8.89873D-02, 8.80593D-02, 8.75126D-02,
     $     8.70800D-02, 8.64489D-02, 8.57928D-02, 8.45060D-02,
     $     8.39904D-02, 8.37828D-02, 8.30653D-02, 8.23737D-02,
     $     8.17747D-02, 8.11060D-02, 8.06408D-02, 8.01274D-02,
     $     7.95411D-02, 7.91136D-02, 7.87861D-02, 7.85656D-02,
     $     7.62605D-02, 7.61228D-02, 7.56564D-02, 7.52221D-02,
     $     7.47624D-02, 7.45494D-02, 7.39871D-02, 7.31203D-02,
     $     7.28307D-02, 7.15357D-02, 7.13514D-02, 7.07419D-02,
     $     7.02995D-02, 7.01146D-02, 6.97879D-02, 6.92042D-02,
     $     6.89307D-02, 6.85174D-02, 6.78206D-02, 6.76528D-02,
     $     6.71527D-02, 6.66989D-02, 6.63150D-02, 6.61045D-02,
     $     6.56408D-02, 6.54295D-02, 6.51224D-02, 6.48197D-02,
     $     6.46861D-02, 6.41081D-02, 6.39448D-02, 6.33896D-02,
     $     6.30979D-02, 6.28495D-02, 6.27086D-02, 6.23891D-02,
     $     6.19457D-02, 6.18650D-02, 6.16232D-02, 6.12773D-02,
     $     6.08745D-02, 6.06967D-02, 6.03543D-02, 6.00954D-02,
     $     5.99005D-02, 5.96910D-02, 5.94903D-02, 5.93443D-02,
     $     5.92037D-02, 5.90372D-02, 5.89686D-02, 5.88177D-02,
     $     5.86857D-02, 5.86059D-02, 5.84107D-02, 5.82634D-02,
     $     5.82579D-02, 5.80559D-02, 5.79609D-02, 5.79164D-02,
     $     5.77898D-02, 5.77388D-02, 5.76778D-02, 5.76046D-02,
     $     5.75100D-02, 5.74750D-02, 5.73840D-02, 5.73635D-02,
     $     5.73061D-02, 5.72776D-02, 5.72648D-02, 5.72230D-02,
     $     5.72040D-02, 5.71940D-02, 5.71783D-02, 5.71605D-02,
     $     5.71468D-02, 5.71451D-02, 5.71417D-02, 5.71333D-02,
     $     5.71328D-02, 5.71404D-02, 5.71421D-02, 5.71527D-02,
     $     5.71594D-02, 5.71740D-02, 5.71903D-02, 5.71991D-02,
     $     5.72328D-02, 5.72443D-02, 5.72605D-02, 5.72933D-02,
     $     5.73302D-02, 5.73428D-02, 5.73641D-02, 5.73900D-02,
     $     5.74359D-02, 5.74603D-02, 5.75182D-02, 5.75426D-02,
     $     5.75602D-02, 5.76223D-02, 5.76463D-02, 5.77013D-02,
     $     5.77393D-02, 5.77686D-02, 5.78197D-02, 5.78526D-02,
     $     5.79200D-02, 5.79775D-02, 5.80117D-02, 5.81387D-02,
     $     5.81887D-02, 5.82229D-02, 5.82899D-02, 5.83416D-02,
     $     5.84288D-02, 5.84919D-02, 5.86199D-02, 5.87047D-02,
     $     5.87402D-02, 5.88035D-02, 5.88324D-02, 5.89414D-02,
     $     5.89909D-02, 5.90880D-02, 5.91808D-02, 5.92615D-02,
     $     5.93572D-02, 5.94557D-02, 5.94796D-02, 5.96926D-02,
     $     5.97035D-02, 5.99241D-02, 6.00228D-02, 6.00988D-02,
     $     6.01449D-02, 6.02616D-02, 6.02786D-02, 6.04396D-02,
     $     6.05134D-02, 6.05989D-02, 6.06241D-02, 6.08511D-02,
     $     6.08930D-02, 6.10768D-02, 6.11941D-02, 6.12214D-02,
     $     6.13709D-02, 6.14778D-02, 6.15160D-02, 6.16198D-02,
     $     6.16986D-02, 6.18926D-02, 6.20255D-02, 6.21084D-02,
     $     6.21887D-02, 6.22686D-02, 6.23081D-02, 6.24837D-02,
     $     6.25929D-02, 6.26743D-02, 6.27230D-02, 6.29197D-02,
     $     6.29833D-02, 6.32237D-02, 6.32623D-02, 6.33692D-02,
     $     6.34345D-02, 6.35561D-02, 6.35925D-02, 6.37799D-02,
     $     6.39421D-02, 6.40816D-02, 6.42519D-02, 6.42691D-02,
     $     6.43440D-02, 6.44584D-02, 6.46709D-02, 6.47796D-02,
     $     6.48959D-02, 6.49741D-02, 6.51715D-02, 6.52365D-02,
     $     6.54119D-02, 6.55302D-02, 6.56066D-02, 6.56432D-02,
     $     6.57474D-02, 6.58349D-02, 6.60691D-02, 6.61205D-02,
     $     6.61279D-02, 6.63435D-02, 6.63728D-02, 6.64596D-02,
     $     6.65706D-02, 6.69086D-02, 6.69609D-02, 6.69721D-02,
     $     6.72200D-02, 6.73085D-02, 6.74515D-02, 6.75036D-02,
     $     6.75297D-02, 6.76269D-02, 6.78123D-02, 6.79401D-02,
     $     6.80717D-02, 6.81245D-02, 6.83215D-02, 6.84071D-02,
     $     6.85577D-02, 6.85997D-02, 6.87896D-02, 6.89251D-02,
     $     6.91642D-02, 6.92449D-02, 6.93723D-02, 6.94383D-02,
     $     6.95285D-02, 6.96078D-02, 6.97438D-02, 6.97598D-02,
     $     6.99316D-02, 7.00558D-02, 7.01686D-02, 7.02482D-02,
     $     7.06259D-02, 7.06943D-02, 7.07803D-02, 7.09220D-02,
     $     7.10491D-02, 7.10863D-02, 7.12479D-02, 7.13019D-02,
     $     7.15242D-02, 7.15939D-02, 7.17796D-02, 7.18456D-02,
     $     7.19203D-02, 7.19826D-02, 7.21143D-02, 7.21641D-02,
     $     7.22927D-02, 7.24582D-02, 7.24958D-02, 7.26792D-02,
     $     7.27658D-02, 7.28619D-02, 7.29378D-02, 7.31690D-02,
     $     7.32961D-02, 7.33428D-02, 7.34335D-02, 7.36576D-02,
     $     7.38707D-02, 7.39034D-02, 7.39658D-02, 7.40523D-02,
     $     7.41568D-02, 7.42407D-02, 7.43738D-02, 7.44615D-02,
     $     7.46795D-02, 7.49279D-02, 7.49526D-02, 7.51180D-02,
     $     7.52207D-02, 7.52531D-02, 7.54251D-02, 7.55386D-02,
     $     7.56471D-02, 7.57121D-02, 7.58490D-02, 7.58892D-02,
     $     7.59297D-02, 7.61677D-02, 7.62601D-02, 7.63432D-02,
     $     7.64485D-02, 7.65643D-02, 7.66713D-02, 7.67612D-02,
     $     7.68202D-02, 7.69666D-02, 7.70608D-02, 7.71708D-02,
     $     7.72181D-02, 7.74204D-02, 7.75164D-02, 7.75902D-02,
     $     7.76800D-02, 7.77335D-02, 7.79676D-02, 7.80967D-02,
     $     7.81339D-02, 7.82424D-02, 7.82765D-02, 7.84141D-02,
     $     7.85297D-02, 7.86351D-02, 7.86796D-02, 7.88648D-02,
     $     7.88979D-02, 7.89408D-02, 7.89990D-02, 7.91361D-02,
     $     7.93020D-02, 7.94037D-02, 7.94648D-02, 7.95077D-02,
     $     7.95692D-02, 7.97156D-02, 7.98302D-02, 7.98626D-02,
     $     7.99597D-02, 8.00163D-02, 8.00838D-02, 8.01889D-02,
     $     8.03570D-02, 8.03957D-02, 8.04346D-02, 8.05356D-02,
     $     8.06183D-02, 8.07291D-02, 8.08196D-02, 8.09138D-02,
     $     8.09664D-02, 8.10287D-02, 8.11670D-02, 8.12488D-02,
     $     8.12660D-02, 8.14078D-02, 8.14760D-02, 8.15464D-02,
     $     8.16113D-02, 8.17350D-02, 8.18652D-02, 8.19284D-02,
     $     8.19445D-02, 8.20957D-02, 8.21658D-02, 8.22750D-02,
     $     8.23881D-02, 8.24385D-02, 8.24721D-02, 8.26538D-02,
     $     8.26774D-02, 8.27573D-02, 8.28299D-02, 8.29262D-02,
     $     8.30081D-02, 8.30510D-02, 8.31092D-02, 8.32345D-02,
     $     8.33921D-02, 8.34568D-02, 8.35866D-02, 8.36191D-02,
     $     8.37524D-02, 8.38177D-02, 8.38656D-02, 8.39297D-02,
     $     8.39865D-02, 8.40387D-02, 8.40923D-02, 8.41271D-02,
     $     8.42028D-02, 8.42926D-02, 8.43589D-02, 8.44066D-02,
     $     8.45291D-02, 8.45476D-02, 8.47011D-02, 8.47755D-02,
     $     8.48622D-02, 8.50013D-02, 8.50294D-02, 8.51138D-02,
     $     8.51794D-02 /

!     Branching ratio of 2007
      data( br(16,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00 /

!     Branching ratio of 4006
      data( br(17,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     7.95640D-27, 2.03836D-25, 4.07598D-24, 1.93050D-23,
     $     8.11635D-22, 3.63139D-21, 6.00601D-21, 5.96951D-14,
     $     1.06267D-13, 2.64486D-13, 5.01839D-10, 2.95416D-09,
     $     4.26483D-09, 1.59760D-08, 2.14063D-08, 8.59081D-08,
     $     1.26489D-07, 9.17142D-07, 1.04166D-06, 1.24296D-06,
     $     2.09492D-06, 2.62166D-06, 3.19513D-06, 3.75541D-06,
     $     5.26811D-06, 6.48345D-06, 7.63956D-06, 1.90076D-05,
     $     2.50559D-05, 2.75661D-05, 3.14006D-05, 3.26813D-05,
     $     4.32966D-05, 5.03057D-05, 5.63646D-05, 6.12418D-05,
     $     7.06359D-05, 7.69472D-05, 8.34408D-05, 8.77299D-05,
     $     1.10964D-04, 1.19779D-04, 1.31315D-04, 1.37878D-04,
     $     1.41555D-04, 1.65411D-04, 1.98486D-04, 2.12560D-04,
     $     2.15183D-04, 2.18336D-04, 2.46178D-04, 2.72021D-04,
     $     2.93463D-04, 3.03305D-04, 3.18809D-04, 3.47283D-04,
     $     4.11557D-04, 4.33855D-04, 4.68940D-04, 4.81029D-04,
     $     4.98683D-04, 5.04584D-04, 5.16190D-04, 5.46454D-04,
     $     5.57388D-04, 5.67006D-04, 5.86173D-04, 6.12757D-04,
     $     6.25404D-04, 6.42764D-04, 6.59884D-04, 6.86914D-04,
     $     6.97037D-04, 7.01091D-04, 7.18249D-04, 7.33360D-04,
     $     7.45682D-04, 7.59253D-04, 7.68637D-04, 7.79077D-04,
     $     7.91301D-04, 8.00377D-04, 8.07560D-04, 8.12479D-04,
     $     8.79597D-04, 8.85138D-04, 8.99744D-04, 9.15602D-04,
     $     9.30046D-04, 9.36563D-04, 9.72632D-04, 1.03661D-03,
     $     1.05413D-03, 1.11751D-03, 1.12537D-03, 1.14919D-03,
     $     1.16547D-03, 1.17212D-03, 1.18367D-03, 1.20655D-03,
     $     1.21959D-03, 1.23569D-03, 1.26147D-03, 1.26929D-03,
     $     1.29360D-03, 1.31294D-03, 1.33007D-03, 1.33886D-03,
     $     1.35972D-03, 1.36832D-03, 1.38219D-03, 1.40206D-03,
     $     1.41198D-03, 1.44142D-03, 1.45947D-03, 1.49326D-03,
     $     1.50769D-03, 1.52130D-03, 1.52827D-03, 1.54505D-03,
     $     1.56926D-03, 1.57341D-03, 1.58547D-03, 1.61391D-03,
     $     1.63949D-03, 1.65092D-03, 1.67112D-03, 1.69009D-03,
     $     1.71090D-03, 1.73722D-03, 1.75675D-03, 1.77122D-03,
     $     1.78483D-03, 1.80323D-03, 1.81132D-03, 1.82965D-03,
     $     1.84558D-03, 1.85489D-03, 1.87917D-03, 1.90379D-03,
     $     1.90465D-03, 1.94241D-03, 1.95857D-03, 1.96605D-03,
     $     1.98830D-03, 1.99767D-03, 2.01014D-03, 2.02842D-03,
     $     2.05984D-03, 2.06975D-03, 2.10092D-03, 2.10829D-03,
     $     2.12882D-03, 2.14043D-03, 2.14621D-03, 2.16443D-03,
     $     2.17488D-03, 2.18141D-03, 2.19564D-03, 2.21721D-03,
     $     2.23503D-03, 2.23819D-03, 2.24826D-03, 2.26639D-03,
     $     2.27240D-03, 2.30816D-03, 2.31106D-03, 2.33005D-03,
     $     2.33905D-03, 2.35391D-03, 2.36696D-03, 2.37755D-03,
     $     2.41296D-03, 2.42043D-03, 2.42972D-03, 2.45234D-03,
     $     2.47551D-03, 2.48291D-03, 2.49592D-03, 2.50870D-03,
     $     2.53153D-03, 2.54050D-03, 2.56158D-03, 2.56872D-03,
     $     2.57525D-03, 2.59177D-03, 2.59927D-03, 2.61242D-03,
     $     2.62048D-03, 2.62635D-03, 2.63767D-03, 2.64800D-03,
     $     2.66505D-03, 2.67659D-03, 2.68252D-03, 2.70765D-03,
     $     2.71683D-03, 2.72429D-03, 2.73798D-03, 2.74670D-03,
     $     2.75885D-03, 2.76857D-03, 2.78517D-03, 2.79461D-03,
     $     2.79897D-03, 2.80656D-03, 2.81134D-03, 2.82349D-03,
     $     2.82929D-03, 2.84065D-03, 2.84969D-03, 2.86013D-03,
     $     2.87009D-03, 2.87966D-03, 2.88318D-03, 2.90360D-03,
     $     2.90459D-03, 2.92562D-03, 2.93442D-03, 2.94121D-03,
     $     2.94598D-03, 2.95593D-03, 2.95717D-03, 2.97027D-03,
     $     2.97744D-03, 2.98629D-03, 2.98850D-03, 3.00377D-03,
     $     3.00620D-03, 3.01701D-03, 3.02547D-03, 3.02703D-03,
     $     3.03578D-03, 3.04187D-03, 3.04413D-03, 3.04928D-03,
     $     3.05341D-03, 3.06472D-03, 3.07067D-03, 3.07436D-03,
     $     3.07792D-03, 3.08166D-03, 3.08463D-03, 3.09401D-03,
     $     3.10018D-03, 3.10383D-03, 3.10671D-03, 3.11591D-03,
     $     3.11843D-03, 3.12855D-03, 3.12996D-03, 3.13369D-03,
     $     3.13588D-03, 3.13983D-03, 3.14097D-03, 3.14740D-03,
     $     3.15632D-03, 3.16172D-03, 3.16793D-03, 3.16850D-03,
     $     3.17135D-03, 3.17737D-03, 3.18624D-03, 3.19025D-03,
     $     3.19521D-03, 3.19852D-03, 3.20510D-03, 3.20702D-03,
     $     3.21265D-03, 3.21637D-03, 3.21959D-03, 3.22089D-03,
     $     3.22459D-03, 3.22716D-03, 3.23336D-03, 3.23463D-03,
     $     3.23481D-03, 3.24032D-03, 3.24104D-03, 3.24345D-03,
     $     3.24683D-03, 3.25916D-03, 3.26083D-03, 3.26116D-03,
     $     3.26931D-03, 3.27197D-03, 3.27611D-03, 3.27832D-03,
     $     3.27912D-03, 3.28169D-03, 3.28631D-03, 3.29015D-03,
     $     3.29317D-03, 3.29424D-03, 3.29909D-03, 3.30132D-03,
     $     3.30449D-03, 3.30523D-03, 3.30868D-03, 3.31131D-03,
     $     3.31632D-03, 3.31821D-03, 3.32065D-03, 3.32177D-03,
     $     3.32308D-03, 3.32437D-03, 3.32682D-03, 3.32704D-03,
     $     3.33046D-03, 3.33242D-03, 3.33399D-03, 3.33529D-03,
     $     3.34094D-03, 3.34220D-03, 3.34338D-03, 3.34659D-03,
     $     3.34871D-03, 3.34922D-03, 3.35113D-03, 3.35171D-03,
     $     3.35435D-03, 3.35516D-03, 3.35714D-03, 3.35784D-03,
     $     3.35854D-03, 3.35909D-03, 3.36045D-03, 3.36090D-03,
     $     3.36287D-03, 3.36559D-03, 3.36605D-03, 3.36813D-03,
     $     3.36917D-03, 3.37011D-03, 3.37104D-03, 3.37451D-03,
     $     3.37571D-03, 3.37609D-03, 3.37723D-03, 3.37988D-03,
     $     3.38301D-03, 3.38346D-03, 3.38412D-03, 3.38510D-03,
     $     3.38598D-03, 3.38659D-03, 3.38767D-03, 3.38829D-03,
     $     3.39024D-03, 3.39267D-03, 3.39284D-03, 3.39432D-03,
     $     3.39521D-03, 3.39551D-03, 3.39748D-03, 3.39903D-03,
     $     3.40093D-03, 3.40158D-03, 3.40281D-03, 3.40312D-03,
     $     3.40355D-03, 3.40568D-03, 3.40627D-03, 3.40691D-03,
     $     3.40831D-03, 3.40968D-03, 3.41069D-03, 3.41149D-03,
     $     3.41190D-03, 3.41266D-03, 3.41314D-03, 3.41413D-03,
     $     3.41439D-03, 3.41539D-03, 3.41591D-03, 3.41614D-03,
     $     3.41691D-03, 3.41726D-03, 3.41946D-03, 3.42014D-03,
     $     3.42037D-03, 3.42175D-03, 3.42207D-03, 3.42307D-03,
     $     3.42354D-03, 3.42432D-03, 3.42474D-03, 3.42624D-03,
     $     3.42646D-03, 3.42685D-03, 3.42716D-03, 3.42786D-03,
     $     3.42847D-03, 3.42868D-03, 3.42885D-03, 3.42899D-03,
     $     3.42944D-03, 3.43077D-03, 3.43168D-03, 3.43187D-03,
     $     3.43225D-03, 3.43254D-03, 3.43278D-03, 3.43356D-03,
     $     3.43507D-03, 3.43524D-03, 3.43537D-03, 3.43575D-03,
     $     3.43659D-03, 3.43737D-03, 3.43783D-03, 3.43808D-03,
     $     3.43846D-03, 3.43880D-03, 3.43933D-03, 3.43986D-03,
     $     3.43991D-03, 3.44062D-03, 3.44089D-03, 3.44143D-03,
     $     3.44177D-03, 3.44253D-03, 3.44353D-03, 3.44392D-03,
     $     3.44408D-03, 3.44454D-03, 3.44471D-03, 3.44556D-03,
     $     3.44631D-03, 3.44646D-03, 3.44652D-03, 3.44750D-03,
     $     3.44779D-03, 3.44856D-03, 3.44892D-03, 3.44915D-03,
     $     3.44934D-03, 3.44962D-03, 3.44998D-03, 3.45044D-03,
     $     3.45078D-03, 3.45082D-03, 3.45078D-03, 3.45075D-03,
     $     3.45098D-03, 3.45129D-03, 3.45133D-03, 3.45133D-03,
     $     3.45178D-03, 3.45243D-03, 3.45268D-03, 3.45277D-03,
     $     3.45293D-03, 3.45321D-03, 3.45361D-03, 3.45395D-03,
     $     3.45443D-03, 3.45450D-03, 3.45511D-03, 3.45551D-03,
     $     3.45607D-03, 3.45646D-03, 3.45648D-03, 3.45652D-03,
     $     3.45691D-03 /

!     Branching ratio of 3006
      data( br(18,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 1.85412D-22, 1.26292D-16, 1.90626D-14,
     $     1.09917D-07, 9.53948D-07, 6.72410D-06, 9.60380D-06,
     $     1.10046D-05, 1.23873D-05, 1.15031D-04, 2.99873D-03,
     $     4.44209D-03, 5.51204D-03, 7.26262D-03, 8.50920D-03,
     $     8.69122D-03, 9.10504D-03, 1.02100D-02, 1.25443D-02,
     $     1.43859D-02, 1.54722D-02, 1.73148D-02, 2.00135D-02,
     $     2.28332D-02, 2.38997D-02, 2.56703D-02, 2.69341D-02,
     $     3.15498D-02, 3.36583D-02, 3.44596D-02, 3.81095D-02,
     $     3.83964D-02, 3.89131D-02, 4.09671D-02, 4.17756D-02,
     $     4.20950D-02, 4.33628D-02, 4.36638D-02, 4.53412D-02,
     $     4.58937D-02, 4.72085D-02, 4.76529D-02, 4.84429D-02,
     $     4.96010D-02, 5.03435D-02, 5.14501D-02, 5.24890D-02,
     $     5.39704D-02, 5.46465D-02, 5.56199D-02, 5.72953D-02,
     $     5.88932D-02, 5.97782D-02, 6.15001D-02, 6.19386D-02,
     $     6.31805D-02, 6.42899D-02, 6.49830D-02, 6.56465D-02,
     $     6.64376D-02, 6.70458D-02, 6.76051D-02, 6.80423D-02,
     $     6.97655D-02, 7.05301D-02, 7.13563D-02, 7.16101D-02,
     $     7.18156D-02, 7.33316D-02, 7.39901D-02, 7.45100D-02,
     $     7.46331D-02, 7.47930D-02, 7.53759D-02, 7.61714D-02,
     $     7.71129D-02, 7.75802D-02, 7.83428D-02, 7.89072D-02,
     $     8.06899D-02, 8.14214D-02, 8.21068D-02, 8.24748D-02,
     $     8.31951D-02, 8.34710D-02, 8.37097D-02, 8.43639D-02,
     $     8.47340D-02, 8.50758D-02, 8.56379D-02, 8.60017D-02,
     $     8.62573D-02, 8.66535D-02, 8.70685D-02, 8.78077D-02,
     $     8.80918D-02, 8.82048D-02, 8.86184D-02, 8.89969D-02,
     $     8.93113D-02, 8.96538D-02, 8.98871D-02, 9.01405D-02,
     $     9.04251D-02, 9.06298D-02, 9.07853D-02, 9.08894D-02,
     $     9.20865D-02, 9.21717D-02, 9.24148D-02, 9.26532D-02,
     $     9.28848D-02, 9.29901D-02, 9.32935D-02, 9.37797D-02,
     $     9.39333D-02, 9.46653D-02, 9.47663D-02, 9.50866D-02,
     $     9.53131D-02, 9.54069D-02, 9.55718D-02, 9.58776D-02,
     $     9.60325D-02, 9.62439D-02, 9.65882D-02, 9.66821D-02,
     $     9.69643D-02, 9.72028D-02, 9.74116D-02, 9.75223D-02,
     $     9.77745D-02, 9.78829D-02, 9.80480D-02, 9.82031D-02,
     $     9.82803D-02, 9.85717D-02, 9.86602D-02, 9.89522D-02,
     $     9.91017D-02, 9.92406D-02, 9.93170D-02, 9.94971D-02,
     $     9.97619D-02, 9.98101D-02, 9.99549D-02, 1.00168D-01,
     $     1.00429D-01, 1.00549D-01, 1.00777D-01, 1.00965D-01,
     $     1.01101D-01, 1.01249D-01, 1.01401D-01, 1.01520D-01,
     $     1.01637D-01, 1.01787D-01, 1.01851D-01, 1.01997D-01,
     $     1.02127D-01, 1.02205D-01, 1.02402D-01, 1.02547D-01,
     $     1.02552D-01, 1.02753D-01, 1.02851D-01, 1.02897D-01,
     $     1.03033D-01, 1.03089D-01, 1.03159D-01, 1.03242D-01,
     $     1.03344D-01, 1.03387D-01, 1.03497D-01, 1.03524D-01,
     $     1.03606D-01, 1.03651D-01, 1.03674D-01, 1.03748D-01,
     $     1.03789D-01, 1.03815D-01, 1.03861D-01, 1.03922D-01,
     $     1.03982D-01, 1.03993D-01, 1.04027D-01, 1.04088D-01,
     $     1.04107D-01, 1.04188D-01, 1.04192D-01, 1.04231D-01,
     $     1.04243D-01, 1.04266D-01, 1.04289D-01, 1.04309D-01,
     $     1.04360D-01, 1.04369D-01, 1.04381D-01, 1.04399D-01,
     $     1.04415D-01, 1.04422D-01, 1.04435D-01, 1.04447D-01,
     $     1.04456D-01, 1.04459D-01, 1.04468D-01, 1.04468D-01,
     $     1.04470D-01, 1.04466D-01, 1.04459D-01, 1.04446D-01,
     $     1.04436D-01, 1.04428D-01, 1.04415D-01, 1.04405D-01,
     $     1.04380D-01, 1.04359D-01, 1.04345D-01, 1.04292D-01,
     $     1.04270D-01, 1.04256D-01, 1.04227D-01, 1.04205D-01,
     $     1.04167D-01, 1.04135D-01, 1.04078D-01, 1.04039D-01,
     $     1.04023D-01, 1.03995D-01, 1.03979D-01, 1.03926D-01,
     $     1.03903D-01, 1.03852D-01, 1.03804D-01, 1.03760D-01,
     $     1.03711D-01, 1.03661D-01, 1.03647D-01, 1.03533D-01,
     $     1.03527D-01, 1.03406D-01, 1.03352D-01, 1.03309D-01,
     $     1.03282D-01, 1.03216D-01, 1.03207D-01, 1.03112D-01,
     $     1.03067D-01, 1.03014D-01, 1.03000D-01, 1.02868D-01,
     $     1.02844D-01, 1.02740D-01, 1.02669D-01, 1.02653D-01,
     $     1.02567D-01, 1.02507D-01, 1.02485D-01, 1.02427D-01,
     $     1.02383D-01, 1.02270D-01, 1.02196D-01, 1.02151D-01,
     $     1.02107D-01, 1.02063D-01, 1.02040D-01, 1.01943D-01,
     $     1.01880D-01, 1.01836D-01, 1.01808D-01, 1.01700D-01,
     $     1.01667D-01, 1.01541D-01, 1.01521D-01, 1.01467D-01,
     $     1.01434D-01, 1.01373D-01, 1.01355D-01, 1.01262D-01,
     $     1.01174D-01, 1.01105D-01, 1.01022D-01, 1.01013D-01,
     $     1.00977D-01, 1.00916D-01, 1.00811D-01, 1.00759D-01,
     $     1.00702D-01, 1.00663D-01, 1.00570D-01, 1.00539D-01,
     $     1.00458D-01, 1.00403D-01, 1.00365D-01, 1.00347D-01,
     $     1.00298D-01, 1.00257D-01, 1.00150D-01, 1.00126D-01,
     $     1.00123D-01, 1.00024D-01, 1.00011D-01, 9.99713D-02,
     $     9.99182D-02, 9.97482D-02, 9.97229D-02, 9.97176D-02,
     $     9.95982D-02, 9.95566D-02, 9.94907D-02, 9.94639D-02,
     $     9.94515D-02, 9.94067D-02, 9.93229D-02, 9.92621D-02,
     $     9.92027D-02, 9.91792D-02, 9.90884D-02, 9.90495D-02,
     $     9.89827D-02, 9.89643D-02, 9.88814D-02, 9.88202D-02,
     $     9.87118D-02, 9.86749D-02, 9.86187D-02, 9.85900D-02,
     $     9.85512D-02, 9.85164D-02, 9.84555D-02, 9.84486D-02,
     $     9.83710D-02, 9.83174D-02, 9.82697D-02, 9.82359D-02,
     $     9.80762D-02, 9.80469D-02, 9.80111D-02, 9.79454D-02,
     $     9.78910D-02, 9.78755D-02, 9.78091D-02, 9.77871D-02,
     $     9.76970D-02, 9.76690D-02, 9.75949D-02, 9.75687D-02,
     $     9.75393D-02, 9.75147D-02, 9.74628D-02, 9.74433D-02,
     $     9.73891D-02, 9.73195D-02, 9.73042D-02, 9.72312D-02,
     $     9.71970D-02, 9.71596D-02, 9.71287D-02, 9.70321D-02,
     $     9.69822D-02, 9.69641D-02, 9.69285D-02, 9.68388D-02,
     $     9.67514D-02, 9.67384D-02, 9.67141D-02, 9.66806D-02,
     $     9.66407D-02, 9.66089D-02, 9.65585D-02, 9.65256D-02,
     $     9.64406D-02, 9.63461D-02, 9.63369D-02, 9.62730D-02,
     $     9.62338D-02, 9.62215D-02, 9.61539D-02, 9.61082D-02,
     $     9.60631D-02, 9.60382D-02, 9.59867D-02, 9.59718D-02,
     $     9.59561D-02, 9.58659D-02, 9.58324D-02, 9.58023D-02,
     $     9.57615D-02, 9.57165D-02, 9.56764D-02, 9.56436D-02,
     $     9.56225D-02, 9.55708D-02, 9.55379D-02, 9.54968D-02,
     $     9.54800D-02, 9.54096D-02, 9.53765D-02, 9.53515D-02,
     $     9.53182D-02, 9.52995D-02, 9.52120D-02, 9.51673D-02,
     $     9.51544D-02, 9.51129D-02, 9.51002D-02, 9.50518D-02,
     $     9.50127D-02, 9.49765D-02, 9.49601D-02, 9.48954D-02,
     $     9.48839D-02, 9.48691D-02, 9.48497D-02, 9.48049D-02,
     $     9.47519D-02, 9.47201D-02, 9.47011D-02, 9.46877D-02,
     $     9.46667D-02, 9.46168D-02, 9.45782D-02, 9.45676D-02,
     $     9.45367D-02, 9.45190D-02, 9.44982D-02, 9.44641D-02,
     $     9.44083D-02, 9.43963D-02, 9.43843D-02, 9.43537D-02,
     $     9.43259D-02, 9.42901D-02, 9.42617D-02, 9.42333D-02,
     $     9.42160D-02, 9.41965D-02, 9.41542D-02, 9.41282D-02,
     $     9.41230D-02, 9.40785D-02, 9.40580D-02, 9.40358D-02,
     $     9.40155D-02, 9.39762D-02, 9.39354D-02, 9.39160D-02,
     $     9.39112D-02, 9.38672D-02, 9.38473D-02, 9.38135D-02,
     $     9.37783D-02, 9.37636D-02, 9.37539D-02, 9.36993D-02,
     $     9.36912D-02, 9.36647D-02, 9.36431D-02, 9.36155D-02,
     $     9.35925D-02, 9.35792D-02, 9.35619D-02, 9.35270D-02,
     $     9.34845D-02, 9.34673D-02, 9.34335D-02, 9.34251D-02,
     $     9.33884D-02, 9.33711D-02, 9.33587D-02, 9.33424D-02,
     $     9.33259D-02, 9.33100D-02, 9.32954D-02, 9.32863D-02,
     $     9.32668D-02, 9.32441D-02, 9.32254D-02, 9.32116D-02,
     $     9.31795D-02, 9.31748D-02, 9.31347D-02, 9.31154D-02,
     $     9.30927D-02, 9.30592D-02, 9.30526D-02, 9.30331D-02,
     $     9.30161D-02 /

!     Branching ratio of 2006
      data( br(19,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 3.57434D-46,
     $     5.39683D-45, 8.91158D-44, 6.29616D-41, 3.31020D-19,
     $     4.61843D-13, 3.92013D-09, 4.78007D-07, 1.51291D-05,
     $     4.93543D-05, 6.35049D-05, 8.55385D-05, 1.00702D-04,
     $     1.71985D-04, 2.57849D-04, 2.96466D-04, 4.60947D-04,
     $     4.75908D-04, 5.05450D-04, 6.32736D-04, 6.89094D-04,
     $     7.14146D-04, 8.21886D-04, 8.47443D-04, 9.60221D-04,
     $     1.00756D-03, 1.11097D-03, 1.15480D-03, 1.23748D-03,
     $     1.32457D-03, 1.37902D-03, 1.49723D-03, 1.62660D-03,
     $     1.80458D-03, 1.87085D-03, 1.99967D-03, 2.21654D-03,
     $     2.46579D-03, 2.62820D-03, 3.01944D-03, 3.12666D-03,
     $     3.35523D-03, 3.62027D-03, 3.74958D-03, 3.90939D-03,
     $     4.07002D-03, 4.20049D-03, 4.30563D-03, 4.40284D-03,
     $     4.71110D-03, 4.87676D-03, 5.03083D-03, 5.06087D-03,
     $     5.09530D-03, 5.35179D-03, 5.44836D-03, 5.53566D-03,
     $     5.55682D-03, 5.58440D-03, 5.66867D-03, 5.77264D-03,
     $     5.90894D-03, 5.97686D-03, 6.08096D-03, 6.13970D-03,
     $     6.37756D-03, 6.51728D-03, 6.65221D-03, 6.71979D-03,
     $     6.88528D-03, 6.95391D-03, 7.00740D-03, 7.15735D-03,
     $     7.24864D-03, 7.33301D-03, 7.44680D-03, 7.50096D-03,
     $     7.55157D-03, 7.61997D-03, 7.68949D-03, 7.83895D-03,
     $     7.89998D-03, 7.92457D-03, 8.00281D-03, 8.08149D-03,
     $     8.15178D-03, 8.23125D-03, 8.28691D-03, 8.34854D-03,
     $     8.41861D-03, 8.46947D-03, 8.50807D-03, 8.53384D-03,
     $     8.76677D-03, 8.77726D-03, 8.82292D-03, 8.86240D-03,
     $     8.90953D-03, 8.93197D-03, 8.98128D-03, 9.05162D-03,
     $     9.07921D-03, 9.20125D-03, 9.22037D-03, 9.28836D-03,
     $     9.33971D-03, 9.36135D-03, 9.39964D-03, 9.46259D-03,
     $     9.48715D-03, 9.53144D-03, 9.60905D-03, 9.62410D-03,
     $     9.66778D-03, 9.71370D-03, 9.75066D-03, 9.77243D-03,
     $     9.81700D-03, 9.83947D-03, 9.86915D-03, 9.89718D-03,
     $     9.90640D-03, 9.96522D-03, 9.97489D-03, 1.00257D-02,
     $     1.00564D-02, 1.00792D-02, 1.00940D-02, 1.01261D-02,
     $     1.01700D-02, 1.01789D-02, 1.02068D-02, 1.02397D-02,
     $     1.02840D-02, 1.03033D-02, 1.03482D-02, 1.03773D-02,
     $     1.03984D-02, 1.04218D-02, 1.04505D-02, 1.04712D-02,
     $     1.04947D-02, 1.05180D-02, 1.05274D-02, 1.05471D-02,
     $     1.05665D-02, 1.05807D-02, 1.06153D-02, 1.06368D-02,
     $     1.06378D-02, 1.06673D-02, 1.06854D-02, 1.06949D-02,
     $     1.07210D-02, 1.07334D-02, 1.07449D-02, 1.07599D-02,
     $     1.07831D-02, 1.07918D-02, 1.08198D-02, 1.08243D-02,
     $     1.08416D-02, 1.08547D-02, 1.08615D-02, 1.08727D-02,
     $     1.08822D-02, 1.08891D-02, 1.09023D-02, 1.09171D-02,
     $     1.09324D-02, 1.09354D-02, 1.09460D-02, 1.09599D-02,
     $     1.09654D-02, 1.09891D-02, 1.09904D-02, 1.10033D-02,
     $     1.10083D-02, 1.10194D-02, 1.10317D-02, 1.10385D-02,
     $     1.10561D-02, 1.10620D-02, 1.10699D-02, 1.10811D-02,
     $     1.10914D-02, 1.10951D-02, 1.11008D-02, 1.11079D-02,
     $     1.11170D-02, 1.11230D-02, 1.11358D-02, 1.11412D-02,
     $     1.11442D-02, 1.11563D-02, 1.11590D-02, 1.11676D-02,
     $     1.11734D-02, 1.11776D-02, 1.11832D-02, 1.11835D-02,
     $     1.11864D-02, 1.11905D-02, 1.11934D-02, 1.11991D-02,
     $     1.12012D-02, 1.12012D-02, 1.12015D-02, 1.12031D-02,
     $     1.12073D-02, 1.12092D-02, 1.12139D-02, 1.12175D-02,
     $     1.12181D-02, 1.12192D-02, 1.12188D-02, 1.12211D-02,
     $     1.12212D-02, 1.12217D-02, 1.12233D-02, 1.12217D-02,
     $     1.12215D-02, 1.12213D-02, 1.12200D-02, 1.12189D-02,
     $     1.12188D-02, 1.12156D-02, 1.12140D-02, 1.12129D-02,
     $     1.12115D-02, 1.12092D-02, 1.12091D-02, 1.12067D-02,
     $     1.12039D-02, 1.11999D-02, 1.11990D-02, 1.11949D-02,
     $     1.11944D-02, 1.11909D-02, 1.11869D-02, 1.11863D-02,
     $     1.11821D-02, 1.11788D-02, 1.11773D-02, 1.11747D-02,
     $     1.11722D-02, 1.11649D-02, 1.11617D-02, 1.11594D-02,
     $     1.11570D-02, 1.11540D-02, 1.11515D-02, 1.11438D-02,
     $     1.11390D-02, 1.11364D-02, 1.11343D-02, 1.11272D-02,
     $     1.11252D-02, 1.11160D-02, 1.11148D-02, 1.11117D-02,
     $     1.11099D-02, 1.11065D-02, 1.11056D-02, 1.10992D-02,
     $     1.10903D-02, 1.10850D-02, 1.10782D-02, 1.10776D-02,
     $     1.10742D-02, 1.10683D-02, 1.10578D-02, 1.10527D-02,
     $     1.10463D-02, 1.10421D-02, 1.10335D-02, 1.10309D-02,
     $     1.10226D-02, 1.10169D-02, 1.10125D-02, 1.10107D-02,
     $     1.10049D-02, 1.10010D-02, 1.09911D-02, 1.09890D-02,
     $     1.09888D-02, 1.09792D-02, 1.09779D-02, 1.09734D-02,
     $     1.09679D-02, 1.09495D-02, 1.09467D-02, 1.09462D-02,
     $     1.09318D-02, 1.09267D-02, 1.09184D-02, 1.09147D-02,
     $     1.09132D-02, 1.09081D-02, 1.08982D-02, 1.08908D-02,
     $     1.08841D-02, 1.08816D-02, 1.08711D-02, 1.08657D-02,
     $     1.08575D-02, 1.08554D-02, 1.08455D-02, 1.08386D-02,
     $     1.08250D-02, 1.08199D-02, 1.08125D-02, 1.08089D-02,
     $     1.08043D-02, 1.08002D-02, 1.07930D-02, 1.07922D-02,
     $     1.07822D-02, 1.07754D-02, 1.07694D-02, 1.07646D-02,
     $     1.07435D-02, 1.07391D-02, 1.07344D-02, 1.07256D-02,
     $     1.07180D-02, 1.07159D-02, 1.07075D-02, 1.07048D-02,
     $     1.06926D-02, 1.06887D-02, 1.06787D-02, 1.06751D-02,
     $     1.06713D-02, 1.06681D-02, 1.06609D-02, 1.06583D-02,
     $     1.06506D-02, 1.06402D-02, 1.06380D-02, 1.06277D-02,
     $     1.06224D-02, 1.06171D-02, 1.06129D-02, 1.05987D-02,
     $     1.05919D-02, 1.05895D-02, 1.05836D-02, 1.05705D-02,
     $     1.05567D-02, 1.05543D-02, 1.05505D-02, 1.05449D-02,
     $     1.05388D-02, 1.05342D-02, 1.05264D-02, 1.05214D-02,
     $     1.05088D-02, 1.04925D-02, 1.04910D-02, 1.04810D-02,
     $     1.04748D-02, 1.04726D-02, 1.04614D-02, 1.04535D-02,
     $     1.04451D-02, 1.04410D-02, 1.04320D-02, 1.04295D-02,
     $     1.04268D-02, 1.04116D-02, 1.04059D-02, 1.04003D-02,
     $     1.03922D-02, 1.03842D-02, 1.03769D-02, 1.03707D-02,
     $     1.03668D-02, 1.03578D-02, 1.03519D-02, 1.03445D-02,
     $     1.03416D-02, 1.03289D-02, 1.03225D-02, 1.03181D-02,
     $     1.03123D-02, 1.03088D-02, 1.02928D-02, 1.02844D-02,
     $     1.02819D-02, 1.02736D-02, 1.02712D-02, 1.02616D-02,
     $     1.02542D-02, 1.02464D-02, 1.02433D-02, 1.02297D-02,
     $     1.02274D-02, 1.02241D-02, 1.02201D-02, 1.02105D-02,
     $     1.01994D-02, 1.01931D-02, 1.01891D-02, 1.01862D-02,
     $     1.01819D-02, 1.01710D-02, 1.01628D-02, 1.01606D-02,
     $     1.01544D-02, 1.01504D-02, 1.01459D-02, 1.01381D-02,
     $     1.01253D-02, 1.01226D-02, 1.01201D-02, 1.01131D-02,
     $     1.01065D-02, 1.00982D-02, 1.00919D-02, 1.00856D-02,
     $     1.00819D-02, 1.00775D-02, 1.00680D-02, 1.00619D-02,
     $     1.00607D-02, 1.00509D-02, 1.00461D-02, 1.00405D-02,
     $     1.00359D-02, 1.00272D-02, 1.00169D-02, 1.00120D-02,
     $     1.00106D-02, 9.99992D-03, 9.99497D-03, 9.98633D-03,
     $     9.97774D-03, 9.97430D-03, 9.97206D-03, 9.95844D-03,
     $     9.95659D-03, 9.95037D-03, 9.94504D-03, 9.93853D-03,
     $     9.93280D-03, 9.92964D-03, 9.92510D-03, 9.91566D-03,
     $     9.90437D-03, 9.90000D-03, 9.89158D-03, 9.88951D-03,
     $     9.88066D-03, 9.87542D-03, 9.87213D-03, 9.86789D-03,
     $     9.86334D-03, 9.85849D-03, 9.85443D-03, 9.85194D-03,
     $     9.84650D-03, 9.83959D-03, 9.83453D-03, 9.83092D-03,
     $     9.82179D-03, 9.82034D-03, 9.80886D-03, 9.80303D-03,
     $     9.79609D-03, 9.78563D-03, 9.78367D-03, 9.77785D-03,
     $     9.77282D-03 /

!     Branching ratio of 3005
      data( br(20,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00 /

!     Branching ratio of 2005
      data( br(21,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00 /

!     Branching ratio of 2004
      data( br(22,k), k=1, 500) /
     $     1.08272D-05, 9.02133D-06, 2.06704D-05,
     $     5.13439D-03, 1.16600D-01, 1.96830D-01, 1.99781D-01,
     $     2.06357D-01, 1.93551D-01, 1.85809D-01, 1.43917D-01,
     $     1.34028D-01, 1.15352D-01, 1.13545D-01, 1.14838D-01,
     $     1.19449D-01, 1.23989D-01, 1.40340D-01, 1.55140D-01,
     $     1.64574D-01, 1.68514D-01, 1.84146D-01, 2.00386D-01,
     $     2.04898D-01, 2.08649D-01, 2.15639D-01, 2.21659D-01,
     $     2.22617D-01, 2.24666D-01, 2.29249D-01, 2.36192D-01,
     $     2.41526D-01, 2.44251D-01, 2.48071D-01, 2.51356D-01,
     $     2.54553D-01, 2.55722D-01, 2.57601D-01, 2.58923D-01,
     $     2.64400D-01, 2.68057D-01, 2.69614D-01, 2.78140D-01,
     $     2.78908D-01, 2.80322D-01, 2.85439D-01, 2.87192D-01,
     $     2.88057D-01, 2.92391D-01, 2.93570D-01, 2.99294D-01,
     $     3.01661D-01, 3.06249D-01, 3.07977D-01, 3.10922D-01,
     $     3.13264D-01, 3.14674D-01, 3.17519D-01, 3.20350D-01,
     $     3.23678D-01, 3.24723D-01, 3.26462D-01, 3.28373D-01,
     $     3.30340D-01, 3.31566D-01, 3.34432D-01, 3.35151D-01,
     $     3.36557D-01, 3.38043D-01, 3.38762D-01, 3.39666D-01,
     $     3.40622D-01, 3.41439D-01, 3.42139D-01, 3.42809D-01,
     $     3.45193D-01, 3.46560D-01, 3.47927D-01, 3.48224D-01,
     $     3.48544D-01, 3.50988D-01, 3.51814D-01, 3.52642D-01,
     $     3.52848D-01, 3.53119D-01, 3.53888D-01, 3.54979D-01,
     $     3.56502D-01, 3.57293D-01, 3.58554D-01, 3.59252D-01,
     $     3.61973D-01, 3.63304D-01, 3.64424D-01, 3.65033D-01,
     $     3.66470D-01, 3.67078D-01, 3.67549D-01, 3.68946D-01,
     $     3.69879D-01, 3.70798D-01, 3.72205D-01, 3.72936D-01,
     $     3.73608D-01, 3.74605D-01, 3.75674D-01, 3.77948D-01,
     $     3.78906D-01, 3.79298D-01, 3.80631D-01, 3.81970D-01,
     $     3.83171D-01, 3.84541D-01, 3.85512D-01, 3.86600D-01,
     $     3.87858D-01, 3.88787D-01, 3.89504D-01, 3.89991D-01,
     $     3.95028D-01, 3.95317D-01, 3.96363D-01, 3.97316D-01,
     $     3.98368D-01, 3.98862D-01, 4.00036D-01, 4.01800D-01,
     $     4.02426D-01, 4.05336D-01, 4.05767D-01, 4.07235D-01,
     $     4.08333D-01, 4.08801D-01, 4.09640D-01, 4.11155D-01,
     $     4.11854D-01, 4.12983D-01, 4.14969D-01, 4.15437D-01,
     $     4.16838D-01, 4.18175D-01, 4.19319D-01, 4.19964D-01,
     $     4.21389D-01, 4.22063D-01, 4.23034D-01, 4.23969D-01,
     $     4.24367D-01, 4.26265D-01, 4.26727D-01, 4.28564D-01,
     $     4.29588D-01, 4.30453D-01, 4.30961D-01, 4.32115D-01,
     $     4.33744D-01, 4.34049D-01, 4.34981D-01, 4.36241D-01,
     $     4.37786D-01, 4.38477D-01, 4.39875D-01, 4.40933D-01,
     $     4.41726D-01, 4.42594D-01, 4.43521D-01, 4.44222D-01,
     $     4.44942D-01, 4.45792D-01, 4.46148D-01, 4.46948D-01,
     $     4.47691D-01, 4.48174D-01, 4.49414D-01, 4.50364D-01,
     $     4.50403D-01, 4.51811D-01, 4.52579D-01, 4.52967D-01,
     $     4.54128D-01, 4.54651D-01, 4.55260D-01, 4.56046D-01,
     $     4.57213D-01, 4.57664D-01, 4.59022D-01, 4.59315D-01,
     $     4.60255D-01, 4.60846D-01, 4.61145D-01, 4.61921D-01,
     $     4.62408D-01, 4.62726D-01, 4.63337D-01, 4.64125D-01,
     $     4.64887D-01, 4.65028D-01, 4.65492D-01, 4.66255D-01,
     $     4.66522D-01, 4.67799D-01, 4.67876D-01, 4.68559D-01,
     $     4.68840D-01, 4.69373D-01, 4.69913D-01, 4.70314D-01,
     $     4.71493D-01, 4.71788D-01, 4.72167D-01, 4.72903D-01,
     $     4.73647D-01, 4.73910D-01, 4.74376D-01, 4.74861D-01,
     $     4.75598D-01, 4.75950D-01, 4.76789D-01, 4.77096D-01,
     $     4.77354D-01, 4.78083D-01, 4.78318D-01, 4.78867D-01,
     $     4.79233D-01, 4.79508D-01, 4.80011D-01, 4.80328D-01,
     $     4.80909D-01, 4.81373D-01, 4.81638D-01, 4.82590D-01,
     $     4.82938D-01, 4.83194D-01, 4.83655D-01, 4.83992D-01,
     $     4.84532D-01, 4.84890D-01, 4.85631D-01, 4.86095D-01,
     $     4.86300D-01, 4.86656D-01, 4.86796D-01, 4.87349D-01,
     $     4.87615D-01, 4.88081D-01, 4.88521D-01, 4.88906D-01,
     $     4.89365D-01, 4.89834D-01, 4.89948D-01, 4.90903D-01,
     $     4.90952D-01, 4.91921D-01, 4.92345D-01, 4.92656D-01,
     $     4.92843D-01, 4.93323D-01, 4.93390D-01, 4.94000D-01,
     $     4.94287D-01, 4.94613D-01, 4.94712D-01, 4.95548D-01,
     $     4.95698D-01, 4.96362D-01, 4.96774D-01, 4.96867D-01,
     $     4.97387D-01, 4.97757D-01, 4.97892D-01, 4.98236D-01,
     $     4.98498D-01, 4.99126D-01, 4.99532D-01, 4.99786D-01,
     $     5.00031D-01, 5.00280D-01, 5.00404D-01, 5.00937D-01,
     $     5.01254D-01, 5.01482D-01, 5.01612D-01, 5.02153D-01,
     $     5.02324D-01, 5.02984D-01, 5.03085D-01, 5.03359D-01,
     $     5.03524D-01, 5.03827D-01, 5.03916D-01, 5.04388D-01,
     $     5.04797D-01, 5.05130D-01, 5.05537D-01, 5.05577D-01,
     $     5.05758D-01, 5.06017D-01, 5.06516D-01, 5.06766D-01,
     $     5.07038D-01, 5.07217D-01, 5.07649D-01, 5.07787D-01,
     $     5.08174D-01, 5.08433D-01, 5.08601D-01, 5.08679D-01,
     $     5.08911D-01, 5.09094D-01, 5.09576D-01, 5.09681D-01,
     $     5.09695D-01, 5.10140D-01, 5.10200D-01, 5.10386D-01,
     $     5.10613D-01, 5.11293D-01, 5.11398D-01, 5.11420D-01,
     $     5.11932D-01, 5.12114D-01, 5.12410D-01, 5.12520D-01,
     $     5.12572D-01, 5.12762D-01, 5.13130D-01, 5.13381D-01,
     $     5.13632D-01, 5.13731D-01, 5.14105D-01, 5.14277D-01,
     $     5.14565D-01, 5.14643D-01, 5.14999D-01, 5.15244D-01,
     $     5.15686D-01, 5.15839D-01, 5.16074D-01, 5.16193D-01,
     $     5.16351D-01, 5.16487D-01, 5.16718D-01, 5.16745D-01,
     $     5.17044D-01, 5.17255D-01, 5.17444D-01, 5.17582D-01,
     $     5.18205D-01, 5.18323D-01, 5.18461D-01, 5.18690D-01,
     $     5.18896D-01, 5.18955D-01, 5.19203D-01, 5.19284D-01,
     $     5.19630D-01, 5.19739D-01, 5.20021D-01, 5.20122D-01,
     $     5.20232D-01, 5.20324D-01, 5.20524D-01, 5.20597D-01,
     $     5.20794D-01, 5.21048D-01, 5.21104D-01, 5.21373D-01,
     $     5.21505D-01, 5.21644D-01, 5.21753D-01, 5.22093D-01,
     $     5.22269D-01, 5.22333D-01, 5.22471D-01, 5.22786D-01,
     $     5.23099D-01, 5.23150D-01, 5.23240D-01, 5.23369D-01,
     $     5.23515D-01, 5.23630D-01, 5.23817D-01, 5.23937D-01,
     $     5.24235D-01, 5.24597D-01, 5.24631D-01, 5.24859D-01,
     $     5.25001D-01, 5.25047D-01, 5.25288D-01, 5.25453D-01,
     $     5.25615D-01, 5.25706D-01, 5.25898D-01, 5.25954D-01,
     $     5.26011D-01, 5.26340D-01, 5.26466D-01, 5.26585D-01,
     $     5.26745D-01, 5.26909D-01, 5.27059D-01, 5.27189D-01,
     $     5.27272D-01, 5.27471D-01, 5.27601D-01, 5.27755D-01,
     $     5.27819D-01, 5.28099D-01, 5.28235D-01, 5.28335D-01,
     $     5.28457D-01, 5.28533D-01, 5.28865D-01, 5.29046D-01,
     $     5.29101D-01, 5.29264D-01, 5.29312D-01, 5.29513D-01,
     $     5.29675D-01, 5.29835D-01, 5.29900D-01, 5.30177D-01,
     $     5.30225D-01, 5.30291D-01, 5.30375D-01, 5.30577D-01,
     $     5.30816D-01, 5.30957D-01, 5.31044D-01, 5.31106D-01,
     $     5.31194D-01, 5.31411D-01, 5.31578D-01, 5.31623D-01,
     $     5.31758D-01, 5.31840D-01, 5.31936D-01, 5.32092D-01,
     $     5.32344D-01, 5.32399D-01, 5.32454D-01, 5.32600D-01,
     $     5.32728D-01, 5.32892D-01, 5.33023D-01, 5.33156D-01,
     $     5.33230D-01, 5.33321D-01, 5.33519D-01, 5.33641D-01,
     $     5.33666D-01, 5.33867D-01, 5.33965D-01, 5.34074D-01,
     $     5.34167D-01, 5.34343D-01, 5.34542D-01, 5.34638D-01,
     $     5.34665D-01, 5.34883D-01, 5.34986D-01, 5.35152D-01,
     $     5.35320D-01, 5.35391D-01, 5.35438D-01, 5.35708D-01,
     $     5.35742D-01, 5.35861D-01, 5.35967D-01, 5.36103D-01,
     $     5.36220D-01, 5.36283D-01, 5.36371D-01, 5.36559D-01,
     $     5.36790D-01, 5.36882D-01, 5.37063D-01, 5.37107D-01,
     $     5.37292D-01, 5.37394D-01, 5.37462D-01, 5.37552D-01,
     $     5.37638D-01, 5.37724D-01, 5.37803D-01, 5.37853D-01,
     $     5.37963D-01, 5.38098D-01, 5.38195D-01, 5.38264D-01,
     $     5.38442D-01, 5.38470D-01, 5.38692D-01, 5.38802D-01,
     $     5.38931D-01, 5.39132D-01, 5.39171D-01, 5.39288D-01,
     $     5.39381D-01 /

!     Branching ratio of 2003
      data( br(23,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 2.39514D-06,
     $     2.31211D-05, 2.75287D-04, 2.17943D-03, 3.25635D-03,
     $     3.60690D-03, 3.72573D-03, 5.00351D-03, 5.57447D-03,
     $     5.90463D-03, 6.47385D-03, 7.56271D-03, 8.87041D-03,
     $     9.21165D-03, 9.70472D-03, 1.04548D-02, 1.10721D-02,
     $     1.11699D-02, 1.15074D-02, 1.20481D-02, 1.30132D-02,
     $     1.36769D-02, 1.41446D-02, 1.46850D-02, 1.51916D-02,
     $     1.56569D-02, 1.58280D-02, 1.60971D-02, 1.62787D-02,
     $     1.71793D-02, 1.77316D-02, 1.80011D-02, 1.94598D-02,
     $     1.95689D-02, 1.97585D-02, 2.04059D-02, 2.05553D-02,
     $     2.06164D-02, 2.08590D-02, 2.09117D-02, 2.12152D-02,
     $     2.13094D-02, 2.15718D-02, 2.16507D-02, 2.17867D-02,
     $     2.19448D-02, 2.20668D-02, 2.22496D-02, 2.24154D-02,
     $     2.26330D-02, 2.27256D-02, 2.28323D-02, 2.29254D-02,
     $     2.30283D-02, 2.30631D-02, 2.30605D-02, 2.30528D-02,
     $     2.30875D-02, 2.30933D-02, 2.31208D-02, 2.31374D-02,
     $     2.31766D-02, 2.32068D-02, 2.32407D-02, 2.32626D-02,
     $     2.33662D-02, 2.34018D-02, 2.34505D-02, 2.34712D-02,
     $     2.34855D-02, 2.35950D-02, 2.36176D-02, 2.36468D-02,
     $     2.36542D-02, 2.36636D-02, 2.36798D-02, 2.37149D-02,
     $     2.37537D-02, 2.37692D-02, 2.37908D-02, 2.37923D-02,
     $     2.38053D-02, 2.37935D-02, 2.37626D-02, 2.37542D-02,
     $     2.37244D-02, 2.37099D-02, 2.36932D-02, 2.36467D-02,
     $     2.36239D-02, 2.36042D-02, 2.35843D-02, 2.35673D-02,
     $     2.35566D-02, 2.35473D-02, 2.35403D-02, 2.35241D-02,
     $     2.35157D-02, 2.35121D-02, 2.35006D-02, 2.34837D-02,
     $     2.34646D-02, 2.34397D-02, 2.34206D-02, 2.33981D-02,
     $     2.33717D-02, 2.33525D-02, 2.33379D-02, 2.33282D-02,
     $     2.32487D-02, 2.32454D-02, 2.32305D-02, 2.32174D-02,
     $     2.31996D-02, 2.31903D-02, 2.31572D-02, 2.30985D-02,
     $     2.30777D-02, 2.29891D-02, 2.29747D-02, 2.29208D-02,
     $     2.28760D-02, 2.28560D-02, 2.28189D-02, 2.27494D-02,
     $     2.27174D-02, 2.26616D-02, 2.25588D-02, 2.25353D-02,
     $     2.24646D-02, 2.23947D-02, 2.23356D-02, 2.23018D-02,
     $     2.22295D-02, 2.21953D-02, 2.21476D-02, 2.20970D-02,
     $     2.20769D-02, 2.19834D-02, 2.19568D-02, 2.18684D-02,
     $     2.18204D-02, 2.17814D-02, 2.17582D-02, 2.17057D-02,
     $     2.16319D-02, 2.16178D-02, 2.15739D-02, 2.15073D-02,
     $     2.14300D-02, 2.13950D-02, 2.13206D-02, 2.12646D-02,
     $     2.12160D-02, 2.11588D-02, 2.11015D-02, 2.10587D-02,
     $     2.10133D-02, 2.09602D-02, 2.09378D-02, 2.08875D-02,
     $     2.08397D-02, 2.08076D-02, 2.07250D-02, 2.06589D-02,
     $     2.06563D-02, 2.05593D-02, 2.05089D-02, 2.04837D-02,
     $     2.04105D-02, 2.03773D-02, 2.03396D-02, 2.02881D-02,
     $     2.02054D-02, 2.01762D-02, 2.00852D-02, 2.00664D-02,
     $     2.00075D-02, 1.99703D-02, 1.99514D-02, 1.99052D-02,
     $     1.98754D-02, 1.98557D-02, 1.98165D-02, 1.97653D-02,
     $     1.97192D-02, 1.97107D-02, 1.96830D-02, 1.96394D-02,
     $     1.96241D-02, 1.95478D-02, 1.95426D-02, 1.95029D-02,
     $     1.94855D-02, 1.94542D-02, 1.94243D-02, 1.94032D-02,
     $     1.93372D-02, 1.93209D-02, 1.93000D-02, 1.92567D-02,
     $     1.92123D-02, 1.91973D-02, 1.91718D-02, 1.91448D-02,
     $     1.91005D-02, 1.90798D-02, 1.90312D-02, 1.90131D-02,
     $     1.89985D-02, 1.89555D-02, 1.89395D-02, 1.89056D-02,
     $     1.88836D-02, 1.88671D-02, 1.88377D-02, 1.88178D-02,
     $     1.87816D-02, 1.87534D-02, 1.87374D-02, 1.86789D-02,
     $     1.86573D-02, 1.86420D-02, 1.86138D-02, 1.85935D-02,
     $     1.85616D-02, 1.85388D-02, 1.84955D-02, 1.84689D-02,
     $     1.84576D-02, 1.84381D-02, 1.84289D-02, 1.83969D-02,
     $     1.83824D-02, 1.83545D-02, 1.83293D-02, 1.83063D-02,
     $     1.82807D-02, 1.82548D-02, 1.82480D-02, 1.81931D-02,
     $     1.81903D-02, 1.81341D-02, 1.81092D-02, 1.80903D-02,
     $     1.80787D-02, 1.80500D-02, 1.80460D-02, 1.80069D-02,
     $     1.79884D-02, 1.79673D-02, 1.79613D-02, 1.79093D-02,
     $     1.78999D-02, 1.78589D-02, 1.78323D-02, 1.78264D-02,
     $     1.77942D-02, 1.77717D-02, 1.77636D-02, 1.77422D-02,
     $     1.77262D-02, 1.76869D-02, 1.76614D-02, 1.76458D-02,
     $     1.76310D-02, 1.76164D-02, 1.76091D-02, 1.75783D-02,
     $     1.75596D-02, 1.75466D-02, 1.75387D-02, 1.75089D-02,
     $     1.74998D-02, 1.74678D-02, 1.74630D-02, 1.74499D-02,
     $     1.74423D-02, 1.74284D-02, 1.74244D-02, 1.74045D-02,
     $     1.73871D-02, 1.73735D-02, 1.73577D-02, 1.73561D-02,
     $     1.73495D-02, 1.73386D-02, 1.73201D-02, 1.73111D-02,
     $     1.73012D-02, 1.72946D-02, 1.72786D-02, 1.72733D-02,
     $     1.72591D-02, 1.72496D-02, 1.72431D-02, 1.72401D-02,
     $     1.72314D-02, 1.72243D-02, 1.72051D-02, 1.72009D-02,
     $     1.72003D-02, 1.71827D-02, 1.71803D-02, 1.71731D-02,
     $     1.71636D-02, 1.71338D-02, 1.71292D-02, 1.71283D-02,
     $     1.71067D-02, 1.70991D-02, 1.70870D-02, 1.70823D-02,
     $     1.70800D-02, 1.70717D-02, 1.70560D-02, 1.70450D-02,
     $     1.70341D-02, 1.70298D-02, 1.70133D-02, 1.70064D-02,
     $     1.69944D-02, 1.69911D-02, 1.69763D-02, 1.69659D-02,
     $     1.69481D-02, 1.69423D-02, 1.69336D-02, 1.69292D-02,
     $     1.69232D-02, 1.69181D-02, 1.69093D-02, 1.69083D-02,
     $     1.68976D-02, 1.68904D-02, 1.68841D-02, 1.68800D-02,
     $     1.68614D-02, 1.68582D-02, 1.68542D-02, 1.68474D-02,
     $     1.68420D-02, 1.68404D-02, 1.68340D-02, 1.68318D-02,
     $     1.68236D-02, 1.68211D-02, 1.68146D-02, 1.68123D-02,
     $     1.68098D-02, 1.68077D-02, 1.68033D-02, 1.68017D-02,
     $     1.67973D-02, 1.67918D-02, 1.67906D-02, 1.67848D-02,
     $     1.67822D-02, 1.67794D-02, 1.67770D-02, 1.67699D-02,
     $     1.67663D-02, 1.67650D-02, 1.67627D-02, 1.67566D-02,
     $     1.67515D-02, 1.67508D-02, 1.67494D-02, 1.67476D-02,
     $     1.67455D-02, 1.67439D-02, 1.67415D-02, 1.67400D-02,
     $     1.67364D-02, 1.67338D-02, 1.67335D-02, 1.67319D-02,
     $     1.67312D-02, 1.67310D-02, 1.67301D-02, 1.67298D-02,
     $     1.67298D-02, 1.67298D-02, 1.67302D-02, 1.67304D-02,
     $     1.67305D-02, 1.67315D-02, 1.67321D-02, 1.67328D-02,
     $     1.67338D-02, 1.67346D-02, 1.67356D-02, 1.67365D-02,
     $     1.67370D-02, 1.67382D-02, 1.67389D-02, 1.67394D-02,
     $     1.67396D-02, 1.67406D-02, 1.67411D-02, 1.67413D-02,
     $     1.67414D-02, 1.67415D-02, 1.67414D-02, 1.67413D-02,
     $     1.67413D-02, 1.67411D-02, 1.67410D-02, 1.67408D-02,
     $     1.67405D-02, 1.67405D-02, 1.67403D-02, 1.67402D-02,
     $     1.67401D-02, 1.67402D-02, 1.67401D-02, 1.67400D-02,
     $     1.67398D-02, 1.67396D-02, 1.67395D-02, 1.67395D-02,
     $     1.67394D-02, 1.67394D-02, 1.67393D-02, 1.67393D-02,
     $     1.67392D-02, 1.67393D-02, 1.67393D-02, 1.67396D-02,
     $     1.67400D-02, 1.67400D-02, 1.67400D-02, 1.67402D-02,
     $     1.67405D-02, 1.67407D-02, 1.67408D-02, 1.67407D-02,
     $     1.67406D-02, 1.67406D-02, 1.67404D-02, 1.67404D-02,
     $     1.67404D-02, 1.67401D-02, 1.67401D-02, 1.67401D-02,
     $     1.67400D-02, 1.67397D-02, 1.67397D-02, 1.67397D-02,
     $     1.67399D-02, 1.67398D-02, 1.67398D-02, 1.67400D-02,
     $     1.67401D-02, 1.67400D-02, 1.67400D-02, 1.67404D-02,
     $     1.67404D-02, 1.67406D-02, 1.67408D-02, 1.67409D-02,
     $     1.67411D-02, 1.67413D-02, 1.67416D-02, 1.67422D-02,
     $     1.67427D-02, 1.67428D-02, 1.67428D-02, 1.67428D-02,
     $     1.67428D-02, 1.67431D-02, 1.67432D-02, 1.67431D-02,
     $     1.67433D-02, 1.67437D-02, 1.67439D-02, 1.67439D-02,
     $     1.67440D-02, 1.67443D-02, 1.67445D-02, 1.67445D-02,
     $     1.67448D-02, 1.67448D-02, 1.67450D-02, 1.67452D-02,
     $     1.67455D-02, 1.67455D-02, 1.67455D-02, 1.67453D-02,
     $     1.67453D-02 /

!     Branching ratio of 1003
      data( br(24,k), k=1, 500) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     2.30576D-12, 3.01271D-07, 2.17375D-06, 8.54124D-06,
     $     3.00699D-05, 7.85513D-05, 5.74311D-04, 1.32427D-03,
     $     1.71886D-03, 1.86745D-03, 2.56148D-03, 3.51292D-03,
     $     3.90852D-03, 4.22370D-03, 4.92387D-03, 5.62287D-03,
     $     5.73891D-03, 5.96553D-03, 6.54022D-03, 7.45917D-03,
     $     8.25700D-03, 8.62684D-03, 9.25466D-03, 9.98951D-03,
     $     1.07043D-02, 1.09742D-02, 1.14196D-02, 1.17413D-02,
     $     1.30379D-02, 1.39897D-02, 1.44323D-02, 1.67022D-02,
     $     1.68772D-02, 1.71877D-02, 1.81406D-02, 1.84582D-02,
     $     1.85951D-02, 1.91707D-02, 1.93038D-02, 1.98451D-02,
     $     2.00709D-02, 2.05114D-02, 2.06941D-02, 2.10254D-02,
     $     2.13147D-02, 2.14874D-02, 2.18828D-02, 2.22927D-02,
     $     2.27726D-02, 2.29120D-02, 2.31598D-02, 2.34268D-02,
     $     2.36403D-02, 2.37673D-02, 2.40679D-02, 2.41452D-02,
     $     2.42622D-02, 2.44283D-02, 2.44938D-02, 2.46000D-02,
     $     2.46942D-02, 2.47792D-02, 2.48429D-02, 2.49112D-02,
     $     2.50964D-02, 2.52196D-02, 2.53346D-02, 2.53480D-02,
     $     2.53739D-02, 2.55943D-02, 2.56746D-02, 2.57590D-02,
     $     2.57799D-02, 2.58070D-02, 2.58776D-02, 2.59552D-02,
     $     2.60565D-02, 2.61017D-02, 2.61588D-02, 2.61708D-02,
     $     2.62446D-02, 2.62926D-02, 2.63188D-02, 2.63238D-02,
     $     2.63514D-02, 2.63630D-02, 2.63667D-02, 2.63785D-02,
     $     2.63914D-02, 2.64050D-02, 2.64106D-02, 2.63989D-02,
     $     2.64033D-02, 2.64038D-02, 2.64045D-02, 2.64343D-02,
     $     2.64484D-02, 2.64542D-02, 2.64616D-02, 2.64692D-02,
     $     2.64743D-02, 2.64764D-02, 2.64757D-02, 2.64729D-02,
     $     2.64686D-02, 2.64653D-02, 2.64627D-02, 2.64609D-02,
     $     2.64068D-02, 2.63989D-02, 2.63877D-02, 2.63722D-02,
     $     2.63617D-02, 2.63569D-02, 2.63283D-02, 2.62688D-02,
     $     2.62517D-02, 2.61547D-02, 2.61404D-02, 2.60934D-02,
     $     2.60558D-02, 2.60390D-02, 2.60078D-02, 2.59398D-02,
     $     2.59003D-02, 2.58452D-02, 2.57469D-02, 2.57184D-02,
     $     2.56300D-02, 2.55526D-02, 2.54836D-02, 2.54461D-02,
     $     2.53598D-02, 2.53218D-02, 2.52635D-02, 2.52045D-02,
     $     2.51758D-02, 2.50660D-02, 2.50290D-02, 2.49158D-02,
     $     2.48581D-02, 2.48055D-02, 2.47766D-02, 2.47084D-02,
     $     2.46098D-02, 2.45921D-02, 2.45385D-02, 2.44535D-02,
     $     2.43535D-02, 2.43066D-02, 2.42165D-02, 2.41398D-02,
     $     2.40774D-02, 2.40061D-02, 2.39363D-02, 2.38818D-02,
     $     2.38274D-02, 2.37567D-02, 2.37260D-02, 2.36558D-02,
     $     2.35917D-02, 2.35517D-02, 2.34474D-02, 2.33619D-02,
     $     2.33586D-02, 2.32314D-02, 2.31672D-02, 2.31356D-02,
     $     2.30395D-02, 2.29970D-02, 2.29438D-02, 2.28738D-02,
     $     2.27683D-02, 2.27285D-02, 2.26101D-02, 2.25829D-02,
     $     2.25002D-02, 2.24508D-02, 2.24260D-02, 2.23545D-02,
     $     2.23121D-02, 2.22851D-02, 2.22324D-02, 2.21610D-02,
     $     2.20956D-02, 2.20838D-02, 2.20458D-02, 2.19806D-02,
     $     2.19587D-02, 2.18505D-02, 2.18436D-02, 2.17873D-02,
     $     2.17636D-02, 2.17213D-02, 2.16811D-02, 2.16488D-02,
     $     2.15489D-02, 2.15259D-02, 2.14970D-02, 2.14349D-02,
     $     2.13708D-02, 2.13485D-02, 2.13086D-02, 2.12685D-02,
     $     2.12042D-02, 2.11759D-02, 2.11070D-02, 2.10828D-02,
     $     2.10609D-02, 2.10034D-02, 2.09827D-02, 2.09390D-02,
     $     2.09107D-02, 2.08897D-02, 2.08491D-02, 2.08193D-02,
     $     2.07680D-02, 2.07302D-02, 2.07095D-02, 2.06300D-02,
     $     2.06013D-02, 2.05785D-02, 2.05378D-02, 2.05101D-02,
     $     2.04687D-02, 2.04398D-02, 2.03834D-02, 2.03499D-02,
     $     2.03343D-02, 2.03073D-02, 2.02949D-02, 2.02539D-02,
     $     2.02335D-02, 2.01980D-02, 2.01665D-02, 2.01352D-02,
     $     2.01006D-02, 2.00658D-02, 2.00556D-02, 1.99837D-02,
     $     1.99800D-02, 1.99044D-02, 1.98712D-02, 1.98467D-02,
     $     1.98309D-02, 1.97929D-02, 1.97879D-02, 1.97398D-02,
     $     1.97153D-02, 1.96865D-02, 1.96781D-02, 1.96141D-02,
     $     1.96030D-02, 1.95528D-02, 1.95194D-02, 1.95123D-02,
     $     1.94722D-02, 1.94437D-02, 1.94331D-02, 1.94076D-02,
     $     1.93877D-02, 1.93391D-02, 1.93102D-02, 1.92923D-02,
     $     1.92750D-02, 1.92573D-02, 1.92471D-02, 1.92087D-02,
     $     1.91864D-02, 1.91719D-02, 1.91629D-02, 1.91293D-02,
     $     1.91197D-02, 1.90831D-02, 1.90781D-02, 1.90651D-02,
     $     1.90575D-02, 1.90443D-02, 1.90405D-02, 1.90201D-02,
     $     1.89983D-02, 1.89842D-02, 1.89672D-02, 1.89657D-02,
     $     1.89579D-02, 1.89451D-02, 1.89227D-02, 1.89121D-02,
     $     1.88996D-02, 1.88913D-02, 1.88738D-02, 1.88685D-02,
     $     1.88522D-02, 1.88412D-02, 1.88330D-02, 1.88295D-02,
     $     1.88187D-02, 1.88111D-02, 1.87919D-02, 1.87878D-02,
     $     1.87872D-02, 1.87689D-02, 1.87664D-02, 1.87582D-02,
     $     1.87480D-02, 1.87143D-02, 1.87093D-02, 1.87083D-02,
     $     1.86826D-02, 1.86737D-02, 1.86591D-02, 1.86527D-02,
     $     1.86501D-02, 1.86412D-02, 1.86239D-02, 1.86114D-02,
     $     1.86001D-02, 1.85959D-02, 1.85787D-02, 1.85702D-02,
     $     1.85574D-02, 1.85542D-02, 1.85394D-02, 1.85295D-02,
     $     1.85110D-02, 1.85044D-02, 1.84951D-02, 1.84907D-02,
     $     1.84853D-02, 1.84807D-02, 1.84728D-02, 1.84720D-02,
     $     1.84612D-02, 1.84545D-02, 1.84491D-02, 1.84447D-02,
     $     1.84273D-02, 1.84235D-02, 1.84201D-02, 1.84131D-02,
     $     1.84077D-02, 1.84065D-02, 1.84018D-02, 1.84004D-02,
     $     1.83936D-02, 1.83914D-02, 1.83864D-02, 1.83846D-02,
     $     1.83829D-02, 1.83816D-02, 1.83781D-02, 1.83770D-02,
     $     1.83729D-02, 1.83669D-02, 1.83659D-02, 1.83611D-02,
     $     1.83584D-02, 1.83561D-02, 1.83542D-02, 1.83472D-02,
     $     1.83448D-02, 1.83441D-02, 1.83412D-02, 1.83361D-02,
     $     1.83300D-02, 1.83289D-02, 1.83277D-02, 1.83256D-02,
     $     1.83242D-02, 1.83235D-02, 1.83220D-02, 1.83214D-02,
     $     1.83201D-02, 1.83175D-02, 1.83176D-02, 1.83176D-02,
     $     1.83179D-02, 1.83179D-02, 1.83180D-02, 1.83178D-02,
     $     1.83169D-02, 1.83176D-02, 1.83191D-02, 1.83198D-02,
     $     1.83203D-02, 1.83243D-02, 1.83264D-02, 1.83278D-02,
     $     1.83285D-02, 1.83303D-02, 1.83324D-02, 1.83340D-02,
     $     1.83353D-02, 1.83393D-02, 1.83417D-02, 1.83438D-02,
     $     1.83450D-02, 1.83499D-02, 1.83519D-02, 1.83539D-02,
     $     1.83557D-02, 1.83566D-02, 1.83600D-02, 1.83625D-02,
     $     1.83630D-02, 1.83634D-02, 1.83638D-02, 1.83655D-02,
     $     1.83678D-02, 1.83687D-02, 1.83692D-02, 1.83709D-02,
     $     1.83713D-02, 1.83715D-02, 1.83723D-02, 1.83741D-02,
     $     1.83766D-02, 1.83787D-02, 1.83797D-02, 1.83804D-02,
     $     1.83811D-02, 1.83820D-02, 1.83831D-02, 1.83836D-02,
     $     1.83853D-02, 1.83860D-02, 1.83872D-02, 1.83881D-02,
     $     1.83893D-02, 1.83900D-02, 1.83908D-02, 1.83925D-02,
     $     1.83929D-02, 1.83942D-02, 1.83956D-02, 1.83976D-02,
     $     1.83984D-02, 1.83993D-02, 1.84018D-02, 1.84028D-02,
     $     1.84031D-02, 1.84057D-02, 1.84069D-02, 1.84075D-02,
     $     1.84085D-02, 1.84106D-02, 1.84117D-02, 1.84124D-02,
     $     1.84124D-02, 1.84152D-02, 1.84166D-02, 1.84177D-02,
     $     1.84192D-02, 1.84204D-02, 1.84211D-02, 1.84241D-02,
     $     1.84243D-02, 1.84254D-02, 1.84267D-02, 1.84291D-02,
     $     1.84310D-02, 1.84319D-02, 1.84327D-02, 1.84350D-02,
     $     1.84383D-02, 1.84399D-02, 1.84435D-02, 1.84444D-02,
     $     1.84478D-02, 1.84487D-02, 1.84498D-02, 1.84514D-02,
     $     1.84521D-02, 1.84521D-02, 1.84530D-02, 1.84538D-02,
     $     1.84554D-02, 1.84570D-02, 1.84582D-02, 1.84591D-02,
     $     1.84616D-02, 1.84619D-02, 1.84651D-02, 1.84664D-02,
     $     1.84677D-02, 1.84706D-02, 1.84713D-02, 1.84735D-02,
     $     1.84748D-02 /

      do j = 1, lines-1
         if ( u >= exen(j) .and. u < exen(j+1) ) then

            do i = 1, nimax
               if ( ifz(i)*1000+ ifa(i) == 1. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(1,j), exen(j+1),
     $                    br(1,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 1001. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(2,j), exen(j+1),
     $                    br(2,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 6010. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(3,j), exen(j+1),
     $                    br(3,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 1002. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(4,j), exen(j+1),
     $                    br(4,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 4010. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(5,j), exen(j+1),
     $                    br(5,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 6009. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(6,j), exen(j+1),
     $                    br(6,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 1003. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(7,j), exen(j+1),
     $                    br(7,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 2003. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(8,j), exen(j+1),
     $                    br(8,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 3009. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(9,j), exen(j+1),
     $                    br(9,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 5008. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(10,j), exen(j+1),
     $                    br(10,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 4008. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(11,j), exen(j+1),
     $                    br(11,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 3008. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(12,j), exen(j+1),
     $                    br(12,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 5007. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(13,j), exen(j+1),
     $                    br(13,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 4007. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(14,j), exen(j+1),
     $                    br(14,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 3007. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(15,j), exen(j+1),
     $                    br(15,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 2007. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(16,j), exen(j+1),
     $                    br(16,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 4006. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(17,j), exen(j+1),
     $                    br(17,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 3006. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(18,j), exen(j+1),
     $                    br(18,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 2006. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(19,j), exen(j+1),
     $                    br(19,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 3005. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(20,j), exen(j+1),
     $                    br(20,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 2005. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(21,j), exen(j+1),
     $                    br(21,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 2004. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(22,j), exen(j+1),
     $                    br(22,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 2003. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(23,j), exen(j+1),
     $                    br(23,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 1003. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(24,j), exen(j+1),
     $                    br(24,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               end if

            end do

            do i = 1, nimax
               r(i) = r(i)/ integral
            end do

         end if

      end do

      end subroutine suppressalpha_c12


! ----------------------------------------------------------------------
!     Modify the decay width of N-14in GDR
      subroutine suppressalpha_n14( u)
!     u: incident photon energy (MeV)
! ----------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

      common /exiejn/ nimax
!$OMP THREADPRIVATE(/exiejn/)
      common /std1/ r(70),s(70),sigma,rr(70)
!$OMP THREADPRIVATE(/std1/)
      common /ejectl/ omega(70),ifa(70),ifz(70)

!     implicit none
      double precision u
      double precision sigma

      double precision exen(500)  ! excitation energy (MeV)
      double precision br(30,500) ! branching ration (-)
      integer lines               ! number of energy bins

      integer i, j, k           ! counters
      double precision integral

      double precision linearone

      integral = 0.

!      7014
!  (1) 7013     1
!  (2) 6013  1001
!  (3) 7012     2
!  (4) 6012  1002
!  (5) 5012  2002
!  (6) 6011  1003
!  (7) 5011  2003
!  (8) 4011  3003
!  (9) 6010  1004
! (10) 5010  2004
! (11) 4010  3004
! (12) 5009  2005
! (13) 4009  3005
! (14) 3009  4005
! (15) 4008  3006
! (16) 3008  4006
! (17) 4007  3007
! (18) 3007  4007
! (19) 3006  4008
! (20) 2006  5008
! (21) 2005  5009
! (22) 2004  5010
! (23) 2003  5011
! (24) 1003  6011
! (25) 1002  6012

      lines = 150
      data( exen(k), k=1, 150) /
     $     8.30700, 13.20200, 14.96900, 16.03700,
     $     17.48900, 18.40500, 18.86600, 20.06000, 21.12900,
     $     21.38000, 21.40900, 21.88900, 22.37600, 22.69100,
     $     22.97400, 23.43000, 23.76600, 24.52400, 24.88500,
     $     24.94600, 25.44200, 25.80500, 26.12100, 26.23500,
     $     26.65100, 26.73900, 27.25300, 27.28300, 27.59100,
     $     27.83700, 28.17500, 28.36200, 28.42400, 28.72500,
     $     28.79200, 28.83700, 29.08200, 29.16000, 29.22800,
     $     29.45000, 29.72600, 29.77700, 29.94200, 29.99200,
     $     30.21200, 30.31600, 30.45300, 30.52300, 30.71600,
     $     30.75500, 30.85000, 30.94600, 31.22000, 31.35000,
     $     31.40200, 31.66900, 31.86700, 31.92500, 31.96800,
     $     32.15700, 32.19700, 32.27000, 32.52300, 32.55400,
     $     32.61200, 32.73900, 32.80100, 32.90400, 32.99400,
     $     33.06800, 33.22500, 33.33200, 33.40700, 33.44200,
     $     33.48900, 33.60500, 33.69100, 33.72600, 33.83200,
     $     33.87400, 33.90000, 34.03900, 34.12500, 34.14800,
     $     34.25000, 34.30900, 34.32900, 34.43300, 34.48300,
     $     34.54200, 34.59300, 34.66300, 34.69400, 34.75400,
     $     34.89400, 34.98500, 35.02200, 35.12300, 35.15600,
     $     35.27900, 35.32600, 35.45600, 35.46700, 35.54000,
     $     35.59700, 35.65600, 35.71900, 35.75500, 35.82700,
     $     35.87400, 35.94300, 36.03400, 36.10000, 36.11600,
     $     36.16200, 36.26700, 36.27700, 36.35700, 36.39500,
     $     36.46000, 36.61900, 36.66100, 36.69600, 36.76700,
     $     36.78400, 36.81400, 36.87300, 36.91200, 36.97000,
     $     37.05400, 37.07200, 37.15100, 37.22300, 37.26200,
     $     37.27500, 37.35200, 37.46500, 37.51800, 37.54700,
     $     37.56900, 37.59000, 37.67200, 37.73400, 37.73800,
     $     37.83100, 37.86000, 37.89700, 37.98200, 38.02000,
     $     38.05200 /

!     Branching ratio of 7013
      data( br(1,k), k=1, 150) /
     $     0.00000D+00, 2.50427D-01, 1.79207D-01,
     $     1.41033D-01, 1.06638D-01, 8.99762D-02, 8.04656D-02,
     $     6.15176D-02, 4.67546D-02, 4.41720D-02, 4.38831D-02,
     $     3.90679D-02, 3.46730D-02, 3.18648D-02, 2.98670D-02,
     $     2.62153D-02, 2.39573D-02, 1.92060D-02, 1.73213D-02,
     $     1.70485D-02, 1.50824D-02, 1.82960D-02, 1.79806D-02,
     $     1.77780D-02, 1.70288D-02, 1.68340D-02, 1.59916D-02,
     $     1.59504D-02, 1.55521D-02, 1.51593D-02, 1.46238D-02,
     $     1.42937D-02, 1.41692D-02, 1.35491D-02, 1.34230D-02,
     $     1.33403D-02, 1.28785D-02, 1.27448D-02, 1.26298D-02,
     $     1.22600D-02, 1.17767D-02, 1.16905D-02, 1.14122D-02,
     $     1.13276D-02, 1.09321D-02, 1.07355D-02, 1.04891D-02,
     $     1.03661D-02, 1.00259D-02, 9.95208D-03, 9.76478D-03,
     $     9.58841D-03, 9.09945D-03, 8.87885D-03, 8.79244D-03,
     $     8.36031D-03, 8.03442D-03, 7.94383D-03, 7.87808D-03,
     $     7.59826D-03, 7.54071D-03, 7.43234D-03, 7.07065D-03,
     $     7.02807D-03, 6.95010D-03, 6.77914D-03, 6.69908D-03,
     $     6.56324D-03, 6.44257D-03, 6.33899D-03, 6.13405D-03,
     $     5.99780D-03, 5.90579D-03, 5.86363D-03, 5.80762D-03,
     $     5.66900D-03, 5.56922D-03, 5.52943D-03, 5.40558D-03,
     $     5.35761D-03, 5.32858D-03, 5.16969D-03, 5.07484D-03,
     $     5.05032D-03, 4.94415D-03, 4.88163D-03, 4.86116D-03,
     $     4.75428D-03, 4.70499D-03, 4.64573D-03, 4.59455D-03,
     $     4.52782D-03, 4.49891D-03, 4.44203D-03, 4.31279D-03,
     $     4.23238D-03, 4.19896D-03, 4.10785D-03, 4.07961D-03,
     $     3.97664D-03, 3.93815D-03, 3.83348D-03, 3.82454D-03,
     $     3.76831D-03, 3.72421D-03, 3.68068D-03, 3.63375D-03,
     $     3.60799D-03, 3.55749D-03, 3.52515D-03, 3.47701D-03,
     $     3.41492D-03, 3.37093D-03, 3.36027D-03, 3.32950D-03,
     $     3.26115D-03, 3.25486D-03, 3.20542D-03, 3.18172D-03,
     $     3.14129D-03, 3.04328D-03, 3.01788D-03, 2.99685D-03,
     $     2.95512D-03, 2.94524D-03, 2.92815D-03, 2.89430D-03,
     $     2.87263D-03, 2.84061D-03, 2.79332D-03, 2.78342D-03,
     $     2.74110D-03, 2.70337D-03, 2.68319D-03, 2.67651D-03,
     $     2.63726D-03, 2.57870D-03, 2.55077D-03, 2.53609D-03,
     $     2.52515D-03, 2.51459D-03, 2.47431D-03, 2.44267D-03,
     $     2.44071D-03, 2.39610D-03, 2.38273D-03, 2.36587D-03,
     $     2.32734D-03, 2.31004D-03, 2.29593D-03 /

!     Branching ratio of 6013
      data( br(2,k), k=1, 150) /
     $     9.99981D-01, 7.19892D-01, 5.58621D-01,
     $     4.65509D-01, 3.77841D-01, 3.32291D-01, 3.03245D-01,
     $     2.44151D-01, 1.94071D-01, 1.85246D-01, 1.84251D-01,
     $     1.67198D-01, 1.51216D-01, 1.40618D-01, 1.34003D-01,
     $     1.26031D-01, 1.20095D-01, 1.01378D-01, 9.32791D-02,
     $     9.21236D-02, 8.38874D-02, 7.81116D-02, 7.42377D-02,
     $     7.29306D-02, 6.83065D-02, 6.71511D-02, 6.13191D-02,
     $     6.10018D-02, 5.78582D-02, 5.51411D-02, 5.16963D-02,
     $     4.97919D-02, 4.91277D-02, 4.59985D-02, 4.53708D-02,
     $     4.49610D-02, 4.27692D-02, 4.21436D-02, 4.16133D-02,
     $     3.99605D-02, 3.79458D-02, 3.75963D-02, 3.64911D-02,
     $     3.61614D-02, 3.46724D-02, 3.39580D-02, 3.30394D-02,
     $     3.25751D-02, 3.13666D-02, 3.11127D-02, 3.04703D-02,
     $     2.98637D-02, 2.82320D-02, 2.75063D-02, 2.72234D-02,
     $     2.58189D-02, 2.47768D-02, 2.44889D-02, 2.42804D-02,
     $     2.33974D-02, 2.32166D-02, 2.28768D-02, 2.17464D-02,
     $     2.16135D-02, 2.13705D-02, 2.08396D-02, 2.05917D-02,
     $     2.01720D-02, 1.97984D-02, 1.94785D-02, 1.88484D-02,
     $     1.84299D-02, 1.81477D-02, 1.80185D-02, 1.78469D-02,
     $     1.74223D-02, 1.71172D-02, 1.69957D-02, 1.66171D-02,
     $     1.64708D-02, 1.63822D-02, 1.58980D-02, 1.56092D-02,
     $     1.55346D-02, 1.52118D-02, 1.50217D-02, 1.49595D-02,
     $     1.46347D-02, 1.44850D-02, 1.43050D-02, 1.41492D-02,
     $     1.39466D-02, 1.38589D-02, 1.36863D-02, 1.32944D-02,
     $     1.30508D-02, 1.29495D-02, 1.26732D-02, 1.25877D-02,
     $     1.22758D-02, 1.21593D-02, 1.18423D-02, 1.18152D-02,
     $     1.16450D-02, 1.15115D-02, 1.13797D-02, 1.12375D-02,
     $     1.11596D-02, 1.10069D-02, 1.09091D-02, 1.07636D-02,
     $     1.05758D-02, 1.04428D-02, 1.04106D-02, 1.03175D-02,
     $     1.01108D-02, 1.00918D-02, 9.94242D-03, 9.87075D-03,
     $     9.74848D-03, 9.45182D-03, 9.37495D-03, 9.31127D-03,
     $     9.18496D-03, 9.15507D-03, 9.10335D-03, 9.00087D-03,
     $     8.93533D-03, 8.83845D-03, 8.69514D-03, 8.66518D-03,
     $     8.53702D-03, 8.42280D-03, 8.36169D-03, 8.34144D-03,
     $     8.22255D-03, 8.04489D-03, 7.96003D-03, 7.91547D-03,
     $     7.88227D-03, 7.85018D-03, 7.72788D-03, 7.63157D-03,
     $     7.62563D-03, 7.49002D-03, 7.44939D-03, 7.39819D-03,
     $     7.28108D-03, 7.22843D-03, 7.18553D-03 /

!     Branching ratio of 7012
      data( br(3,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 3.33639D-25, 5.61856D-10, 2.97657D-08,
     $     8.41484D-08, 1.92853D-07, 2.59877D-07, 2.88091D-07,
     $     4.41896D-07, 6.02382D-07, 6.46432D-07, 6.77629D-07,
     $     1.05463D-06, 1.16753D-06, 1.41846D-06, 2.47770D-06,
     $     2.60963D-06, 2.85637D-06, 3.49838D-06, 3.80371D-06,
     $     4.61391D-06, 7.48852D-06, 1.00537D-05, 1.28622D-05,
     $     1.62781D-05, 1.78916D-05, 1.85004D-05, 1.92647D-05,
     $     2.15929D-05, 2.34305D-05, 2.42033D-05, 3.72931D-05,
     $     4.11103D-05, 4.33061D-05, 6.21383D-05, 7.47053D-05,
     $     7.77148D-05, 9.05277D-05, 1.03199D-04, 1.06423D-04,
     $     1.25763D-04, 1.33249D-04, 1.48104D-04, 1.91712D-04,
     $     2.17754D-04, 2.26224D-04, 2.43290D-04, 2.79743D-04,
     $     3.00507D-04, 3.10928D-04, 3.57165D-04, 3.67548D-04,
     $     4.03676D-04, 4.17247D-04, 4.54460D-04, 4.58425D-04,
     $     4.79035D-04, 4.95885D-04, 5.11008D-04, 5.44281D-04,
     $     5.55717D-04, 5.76159D-04, 5.88573D-04, 6.07888D-04,
     $     6.33877D-04, 6.58301D-04, 6.65845D-04, 6.83669D-04,
     $     7.17662D-04, 7.20521D-04, 7.42258D-04, 7.53209D-04,
     $     7.72410D-04, 8.25817D-04, 8.39354D-04, 8.51330D-04,
     $     8.76675D-04, 8.83163D-04, 8.92973D-04, 9.11800D-04,
     $     9.22656D-04, 9.38369D-04, 9.69625D-04, 9.75248D-04,
     $     9.96322D-04, 1.01320D-03, 1.02191D-03, 1.02476D-03,
     $     1.04118D-03, 1.07428D-03, 1.09200D-03, 1.10001D-03,
     $     1.10550D-03, 1.11099D-03, 1.13248D-03, 1.16194D-03,
     $     1.16330D-03, 1.19082D-03, 1.19785D-03, 1.20644D-03,
     $     1.22607D-03, 1.23839D-03, 1.24597D-03 /

!     Branching ratio of 6012
      data( br(4,k), k=1, 150) /
     $     0.00000D+00, 2.96415D-02, 2.54692D-01,
     $     2.92902D-01, 4.09774D-01, 4.85346D-01, 5.27528D-01,
     $     5.70127D-01, 6.11137D-01, 6.16738D-01, 6.17258D-01,
     $     6.22713D-01, 6.23806D-01, 6.22899D-01, 6.21663D-01,
     $     6.14101D-01, 6.03003D-01, 5.66296D-01, 5.50416D-01,
     $     5.48122D-01, 5.31460D-01, 5.19256D-01, 5.11375D-01,
     $     5.08630D-01, 4.97676D-01, 4.94423D-01, 4.80295D-01,
     $     4.79640D-01, 4.73258D-01, 4.66112D-01, 4.54767D-01,
     $     4.46988D-01, 4.44026D-01, 4.29288D-01, 4.26288D-01,
     $     4.24324D-01, 4.13653D-01, 4.10606D-01, 4.08020D-01,
     $     4.00052D-01, 3.90324D-01, 3.88661D-01, 3.83345D-01,
     $     3.81725D-01, 3.74054D-01, 3.70169D-01, 3.65668D-01,
     $     3.63515D-01, 3.56734D-01, 3.55197D-01, 3.51340D-01,
     $     3.47817D-01, 3.37482D-01, 3.32735D-01, 3.30856D-01,
     $     3.21275D-01, 3.13865D-01, 3.11796D-01, 3.10292D-01,
     $     3.03885D-01, 3.02562D-01, 3.00043D-01, 2.91529D-01,
     $     2.90515D-01, 2.88654D-01, 2.84527D-01, 2.82587D-01,
     $     2.79239D-01, 2.76357D-01, 2.73808D-01, 2.68547D-01,
     $     2.64989D-01, 2.62540D-01, 2.61405D-01, 2.59886D-01,
     $     2.56069D-01, 2.53296D-01, 2.52185D-01, 2.48719D-01,
     $     2.47358D-01, 2.46534D-01, 2.41953D-01, 2.39198D-01,
     $     2.38485D-01, 2.35381D-01, 2.33525D-01, 2.32918D-01,
     $     2.29712D-01, 2.28227D-01, 2.26419D-01, 2.24892D-01,
     $     2.22861D-01, 2.21974D-01, 2.20210D-01, 2.16165D-01,
     $     2.13631D-01, 2.12565D-01, 2.09661D-01, 2.08759D-01,
     $     2.05449D-01, 2.04204D-01, 2.00794D-01, 2.00500D-01,
     $     1.98661D-01, 1.97204D-01, 1.95770D-01, 1.94240D-01,
     $     1.93390D-01, 1.91714D-01, 1.90635D-01, 1.89009D-01,
     $     1.86896D-01, 1.85402D-01, 1.85040D-01, 1.83983D-01,
     $     1.81613D-01, 1.81395D-01, 1.79672D-01, 1.78837D-01,
     $     1.77402D-01, 1.73875D-01, 1.72952D-01, 1.72184D-01,
     $     1.70660D-01, 1.70298D-01, 1.69672D-01, 1.68420D-01,
     $     1.67618D-01, 1.66424D-01, 1.64637D-01, 1.64262D-01,
     $     1.62650D-01, 1.61202D-01, 1.60423D-01, 1.60164D-01,
     $     1.58637D-01, 1.56323D-01, 1.55202D-01, 1.54613D-01,
     $     1.54175D-01, 1.53749D-01, 1.52124D-01, 1.50833D-01,
     $     1.50753D-01, 1.48927D-01, 1.48379D-01, 1.47688D-01,
     $     1.46095D-01, 1.45376D-01, 1.44789D-01 /

!     Branching ratio of 5012
      data( br(5,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 7.37219D-32, 4.37258D-08, 1.13651D-06,
     $     2.69225D-06, 2.98298D-05, 4.89291D-05, 4.67517D-04,
     $     5.11927D-04, 1.09220D-03, 1.91046D-03, 4.32497D-03,
     $     6.15109D-03, 6.79955D-03, 9.74183D-03, 1.02929D-02,
     $     1.06401D-02, 1.21842D-02, 1.26022D-02, 1.29505D-02,
     $     1.40201D-02, 1.52345D-02, 1.54643D-02, 1.62403D-02,
     $     1.64853D-02, 1.76048D-02, 1.82435D-02, 1.93281D-02,
     $     1.99895D-02, 2.21052D-02, 2.25443D-02, 2.35788D-02,
     $     2.45985D-02, 2.72138D-02, 2.83704D-02, 2.88222D-02,
     $     3.09629D-02, 3.22541D-02, 3.26245D-02, 3.29028D-02,
     $     3.41542D-02, 3.44227D-02, 3.48884D-02, 3.64935D-02,
     $     3.66922D-02, 3.70687D-02, 3.78651D-02, 3.82556D-02,
     $     3.88545D-02, 3.93319D-02, 3.96766D-02, 4.04496D-02,
     $     4.09504D-02, 4.13021D-02, 4.14662D-02, 4.16850D-02,
     $     4.21980D-02, 4.25770D-02, 4.27306D-02, 4.31426D-02,
     $     4.33040D-02, 4.34055D-02, 4.38752D-02, 4.41625D-02,
     $     4.42422D-02, 4.45986D-02, 4.47808D-02, 4.48464D-02,
     $     4.51617D-02, 4.53190D-02, 4.54808D-02, 4.56064D-02,
     $     4.57948D-02, 4.58774D-02, 4.60154D-02, 4.63144D-02,
     $     4.65122D-02, 4.65784D-02, 4.67427D-02, 4.68091D-02,
     $     4.70575D-02, 4.71504D-02, 4.73891D-02, 4.74051D-02,
     $     4.75380D-02, 4.76242D-02, 4.77268D-02, 4.78134D-02,
     $     4.78697D-02, 4.79797D-02, 4.80488D-02, 4.81246D-02,
     $     4.82132D-02, 4.82714D-02, 4.82830D-02, 4.83086D-02,
     $     4.83635D-02, 4.83698D-02, 4.84206D-02, 4.84336D-02,
     $     4.84474D-02, 4.84432D-02, 4.84386D-02, 4.84333D-02,
     $     4.84273D-02, 4.84255D-02, 4.84260D-02, 4.84150D-02,
     $     4.84141D-02, 4.84078D-02, 4.83645D-02, 4.83567D-02,
     $     4.83306D-02, 4.83072D-02, 4.82933D-02, 4.82885D-02,
     $     4.82577D-02, 4.81653D-02, 4.81000D-02, 4.80718D-02,
     $     4.80523D-02, 4.80299D-02, 4.79468D-02, 4.78455D-02,
     $     4.78403D-02, 4.77190D-02, 4.76864D-02, 4.76453D-02,
     $     4.75396D-02, 4.74839D-02, 4.74435D-02 /

!     Branching ratio of 6011
      data( br(6,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 9.95604D-07,
     $     2.38917D-03, 7.85089D-03, 1.81282D-02, 2.07099D-02,
     $     2.12168D-02, 2.61593D-02, 2.81904D-02, 2.91182D-02,
     $     2.92626D-02, 2.90115D-02, 2.87506D-02, 2.70152D-02,
     $     2.69006D-02, 2.56940D-02, 2.46926D-02, 2.47398D-02,
     $     2.54570D-02, 2.57251D-02, 2.70836D-02, 2.73889D-02,
     $     2.75889D-02, 2.84742D-02, 2.87078D-02, 2.88791D-02,
     $     2.92519D-02, 2.93296D-02, 2.93243D-02, 2.92757D-02,
     $     2.92505D-02, 2.90905D-02, 2.89824D-02, 2.87937D-02,
     $     2.86737D-02, 2.83287D-02, 2.82620D-02, 2.81135D-02,
     $     2.79130D-02, 2.74354D-02, 2.71835D-02, 2.70777D-02,
     $     2.65135D-02, 2.62259D-02, 2.61277D-02, 2.60506D-02,
     $     2.57025D-02, 2.56305D-02, 2.55470D-02, 2.53196D-02,
     $     2.52942D-02, 2.52389D-02, 2.51682D-02, 2.51193D-02,
     $     2.51314D-02, 2.50753D-02, 2.51646D-02, 2.53341D-02,
     $     2.54252D-02, 2.54972D-02, 2.55336D-02, 2.55833D-02,
     $     2.57531D-02, 2.58703D-02, 2.59138D-02, 2.60649D-02,
     $     2.61377D-02, 2.61765D-02, 2.64967D-02, 2.66716D-02,
     $     2.67098D-02, 2.68594D-02, 2.69835D-02, 2.70159D-02,
     $     2.72153D-02, 2.72939D-02, 2.74255D-02, 2.75081D-02,
     $     2.76277D-02, 2.76823D-02, 2.78261D-02, 2.81820D-02,
     $     2.83929D-02, 2.85035D-02, 2.87909D-02, 2.88661D-02,
     $     2.91451D-02, 2.92526D-02, 2.95667D-02, 2.95998D-02,
     $     2.97525D-02, 2.98938D-02, 2.99972D-02, 3.01033D-02,
     $     3.01560D-02, 3.02541D-02, 3.03156D-02, 3.04481D-02,
     $     3.06221D-02, 3.07331D-02, 3.07609D-02, 3.08626D-02,
     $     3.10885D-02, 3.11063D-02, 3.12389D-02, 3.13223D-02,
     $     3.14771D-02, 3.19172D-02, 3.20354D-02, 3.21373D-02,
     $     3.23241D-02, 3.23678D-02, 3.24370D-02, 3.25988D-02,
     $     3.26880D-02, 3.28270D-02, 3.30937D-02, 3.31463D-02,
     $     3.33478D-02, 3.35206D-02, 3.36143D-02, 3.36456D-02,
     $     3.38331D-02, 3.42109D-02, 3.44392D-02, 3.45433D-02,
     $     3.46165D-02, 3.46950D-02, 3.49812D-02, 3.52861D-02,
     $     3.53023D-02, 3.56641D-02, 3.57594D-02, 3.58768D-02,
     $     3.61639D-02, 3.63116D-02, 3.64148D-02 /

!     Branching ratio of 5011
      data( br(7,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 2.45050D-07, 4.46781D-05, 6.41049D-05,
     $     2.58140D-03, 1.08366D-02, 1.69139D-02, 2.19030D-02,
     $     2.79156D-02, 3.24356D-02, 3.83286D-02, 3.86509D-02,
     $     3.86272D-02, 3.76475D-02, 3.59285D-02, 3.45815D-02,
     $     3.41510D-02, 3.36782D-02, 3.37684D-02, 3.59067D-02,
     $     3.60678D-02, 3.76536D-02, 3.84603D-02, 3.98438D-02,
     $     4.14755D-02, 4.21803D-02, 4.63840D-02, 4.74249D-02,
     $     4.81252D-02, 5.17870D-02, 5.29682D-02, 5.39898D-02,
     $     5.71997D-02, 6.09428D-02, 6.16862D-02, 6.42474D-02,
     $     6.50758D-02, 6.88755D-02, 7.07751D-02, 7.34617D-02,
     $     7.49098D-02, 7.91824D-02, 8.00438D-02, 8.23708D-02,
     $     8.50552D-02, 9.20614D-02, 9.54335D-02, 9.67875D-02,
     $     1.03763D-01, 1.08648D-01, 1.10053D-01, 1.11084D-01,
     $     1.15418D-01, 1.16283D-01, 1.17748D-01, 1.22393D-01,
     $     1.22923D-01, 1.23905D-01, 1.25911D-01, 1.26868D-01,
     $     1.28393D-01, 1.30487D-01, 1.31882D-01, 1.34359D-01,
     $     1.36454D-01, 1.37761D-01, 1.38328D-01, 1.39062D-01,
     $     1.40750D-01, 1.41952D-01, 1.42433D-01, 1.44327D-01,
     $     1.44953D-01, 1.45329D-01, 1.47160D-01, 1.48245D-01,
     $     1.48537D-01, 1.49830D-01, 1.50545D-01, 1.50790D-01,
     $     1.52013D-01, 1.52611D-01, 1.53290D-01, 1.54151D-01,
     $     1.55082D-01, 1.55470D-01, 1.56173D-01, 1.57809D-01,
     $     1.58919D-01, 1.59353D-01, 1.60893D-01, 1.61342D-01,
     $     1.62937D-01, 1.63520D-01, 1.65054D-01, 1.65177D-01,
     $     1.66015D-01, 1.66624D-01, 1.67261D-01, 1.68179D-01,
     $     1.68589D-01, 1.69361D-01, 1.69842D-01, 1.70489D-01,
     $     1.71323D-01, 1.72039D-01, 1.72236D-01, 1.72706D-01,
     $     1.73661D-01, 1.73750D-01, 1.74455D-01, 1.74769D-01,
     $     1.75296D-01, 1.76576D-01, 1.76912D-01, 1.77194D-01,
     $     1.77880D-01, 1.78066D-01, 1.78363D-01, 1.78892D-01,
     $     1.79239D-01, 1.79743D-01, 1.80618D-01, 1.80783D-01,
     $     1.81497D-01, 1.82141D-01, 1.82493D-01, 1.82611D-01,
     $     1.83315D-01, 1.84476D-01, 1.84962D-01, 1.85232D-01,
     $     1.85436D-01, 1.85624D-01, 1.86388D-01, 1.87124D-01,
     $     1.87165D-01, 1.88042D-01, 1.88303D-01, 1.88632D-01,
     $     1.89350D-01, 1.89774D-01, 1.90068D-01 /

!     Branching ratio of 4011
      data( br(8,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00 /

!     Branching ratio of 6010
      data( br(9,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 3.11066D-44,
     $     9.67215D-41, 1.09260D-10, 1.63732D-10, 4.89701D-10,
     $     3.00200D-09, 3.12767D-09, 3.84916D-09, 2.30912D-08,
     $     8.68153D-08, 7.30975D-07, 9.95726D-07, 1.35476D-06,
     $     2.42311D-06, 2.77064D-06, 3.15567D-06, 4.01641D-06,
     $     4.34276D-06, 4.91617D-06, 8.57556D-06, 9.11608D-06,
     $     1.05741D-05, 1.14367D-05, 1.18381D-05, 1.19652D-05,
     $     1.26990D-05, 1.83962D-05, 2.25271D-05, 2.40147D-05,
     $     2.48815D-05, 2.58933D-05, 3.00652D-05, 4.05099D-05,
     $     4.09083D-05, 4.84742D-05, 5.01380D-05, 5.22589D-05,
     $     5.81352D-05, 6.34693D-05, 6.66906D-05 /

!     Branching ratio of 5010
      data( br(10,k), k=1, 150) /
     $     0.00000D+00, 1.04047D-12, 7.41832D-03,
     $     1.00472D-01, 1.05598D-01, 9.21637D-02, 8.84882D-02,
     $     1.21964D-01, 1.27106D-01, 1.25354D-01, 1.25062D-01,
     $     1.17530D-01, 1.07754D-01, 1.00477D-01, 9.51908D-02,
     $     8.46467D-02, 7.80153D-02, 6.35096D-02, 5.76508D-02,
     $     5.68023D-02, 5.06547D-02, 4.63946D-02, 4.34475D-02,
     $     4.24553D-02, 3.90155D-02, 3.82076D-02, 3.41612D-02,
     $     3.39455D-02, 3.18390D-02, 3.00963D-02, 2.79087D-02,
     $     2.67280D-02, 2.63244D-02, 2.44552D-02, 2.40845D-02,
     $     2.38433D-02, 2.25704D-02, 2.22091D-02, 2.19041D-02,
     $     2.09599D-02, 1.98268D-02, 1.96317D-02, 1.90184D-02,
     $     1.88364D-02, 1.80217D-02, 1.76337D-02, 1.71364D-02,
     $     1.68856D-02, 1.62343D-02, 1.60984D-02, 1.57560D-02,
     $     1.54333D-02, 1.45701D-02, 1.41883D-02, 1.40397D-02,
     $     1.33040D-02, 1.27610D-02, 1.26113D-02, 1.25029D-02,
     $     1.20446D-02, 1.19508D-02, 1.17749D-02, 1.11909D-02,
     $     1.11224D-02, 1.09973D-02, 1.07247D-02, 1.05980D-02,
     $     1.03849D-02, 1.01970D-02, 1.00376D-02, 9.73134D-03,
     $     9.53508D-03, 9.40738D-03, 9.35062D-03, 9.27711D-03,
     $     9.10442D-03, 8.99048D-03, 8.94795D-03, 8.82446D-03,
     $     8.78282D-03, 8.75986D-03, 8.66314D-03, 8.63808D-03,
     $     8.63625D-03, 8.65106D-03, 8.67474D-03, 8.68614D-03,
     $     8.77053D-03, 8.82262D-03, 8.89395D-03, 8.95951D-03,
     $     9.06332D-03, 9.11389D-03, 9.22324D-03, 9.51640D-03,
     $     9.72677D-03, 9.82164D-03, 1.00797D-02, 1.01640D-02,
     $     1.04818D-02, 1.06041D-02, 1.09493D-02, 1.09812D-02,
     $     1.11778D-02, 1.13447D-02, 1.15152D-02, 1.17040D-02,
     $     1.18174D-02, 1.20567D-02, 1.22214D-02, 1.24882D-02,
     $     1.28626D-02, 1.31376D-02, 1.32057D-02, 1.34124D-02,
     $     1.38964D-02, 1.39422D-02, 1.43117D-02, 1.44984D-02,
     $     1.48265D-02, 1.56795D-02, 1.59072D-02, 1.60994D-02,
     $     1.64748D-02, 1.65629D-02, 1.67167D-02, 1.70291D-02,
     $     1.72281D-02, 1.75266D-02, 1.79750D-02, 1.80708D-02,
     $     1.84790D-02, 1.88459D-02, 1.90439D-02, 1.91100D-02,
     $     1.95026D-02, 2.01102D-02, 2.04234D-02, 2.05875D-02,
     $     2.07101D-02, 2.08305D-02, 2.12959D-02, 2.16755D-02,
     $     2.16992D-02, 2.22579D-02, 2.24302D-02, 2.26521D-02,
     $     2.31823D-02, 2.34225D-02, 2.36273D-02 /

!     Branching ratio of 4010
      data( br(11,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 1.29435D-12,
     $     4.76408D-12, 1.47607D-11, 3.18700D-10, 9.07755D-10,
     $     7.96582D-09, 3.15406D-08, 7.12183D-08, 2.49268D-07,
     $     5.23591D-07, 9.13170D-07, 1.17977D-06, 1.64152D-06,
     $     3.36328D-06, 5.18278D-06, 6.05881D-06, 9.30640D-06,
     $     1.09178D-05, 1.20316D-05, 2.00748D-05, 2.76889D-05,
     $     3.01452D-05, 4.38710D-05, 5.42357D-05, 5.81573D-05,
     $     8.19626D-05, 9.56374D-05, 1.13566D-04, 1.30352D-04,
     $     1.55149D-04, 1.66688D-04, 1.90224D-04, 2.49006D-04,
     $     2.89816D-04, 3.06842D-04, 3.53485D-04, 3.69096D-04,
     $     4.27047D-04, 4.48534D-04, 5.06560D-04, 5.11324D-04,
     $     5.42818D-04, 5.66670D-04, 5.90963D-04, 6.15996D-04,
     $     6.30106D-04, 6.57515D-04, 6.74765D-04, 6.99056D-04,
     $     7.29933D-04, 7.51823D-04, 7.57087D-04, 7.72169D-04,
     $     8.08062D-04, 8.11715D-04, 8.42715D-04, 8.58899D-04,
     $     8.89710D-04, 9.91151D-04, 1.02598D-03, 1.05793D-03,
     $     1.13118D-03, 1.15042D-03, 1.18605D-03, 1.26156D-03,
     $     1.31576D-03, 1.40152D-03, 1.53384D-03, 1.56353D-03,
     $     1.69916D-03, 1.82605D-03, 1.89571D-03, 1.91901D-03,
     $     2.05774D-03, 2.25494D-03, 2.34430D-03, 2.39272D-03,
     $     2.42901D-03, 2.46291D-03, 2.59094D-03, 2.68154D-03,
     $     2.68734D-03, 2.81834D-03, 2.85804D-03, 2.90779D-03,
     $     3.01753D-03, 3.06431D-03, 3.10341D-03 /

!     Branching ratio of 5010
      data( br(12,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00 /

!     Branching ratio of 4009
      data( br(13,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 5.27760D-10,
     $     1.76526D-03, 1.99838D-02, 2.61593D-02, 2.69377D-02,
     $     4.15986D-02, 5.56732D-02, 6.36824D-02, 7.02559D-02,
     $     7.97042D-02, 8.72776D-02, 1.11644D-01, 1.22618D-01,
     $     1.24194D-01, 1.34416D-01, 1.37978D-01, 1.41226D-01,
     $     1.42249D-01, 1.45039D-01, 1.44951D-01, 1.44838D-01,
     $     1.44810D-01, 1.44270D-01, 1.42779D-01, 1.40659D-01,
     $     1.39200D-01, 1.38589D-01, 1.35542D-01, 1.34976D-01,
     $     1.34612D-01, 1.32464D-01, 1.31916D-01, 1.31455D-01,
     $     1.30000D-01, 1.27854D-01, 1.27475D-01, 1.26220D-01,
     $     1.25828D-01, 1.23811D-01, 1.22722D-01, 1.21263D-01,
     $     1.20489D-01, 1.18473D-01, 1.18005D-01, 1.16739D-01,
     $     1.15568D-01, 1.12315D-01, 1.10832D-01, 1.10248D-01,
     $     1.07257D-01, 1.04832D-01, 1.04166D-01, 1.03684D-01,
     $     1.01630D-01, 1.01205D-01, 1.00371D-01, 9.75551D-02,
     $     9.72211D-02, 9.66125D-02, 9.52472D-02, 9.46139D-02,
     $     9.35005D-02, 9.24706D-02, 9.15463D-02, 8.97620D-02,
     $     8.85600D-02, 8.77544D-02, 8.73861D-02, 8.68965D-02,
     $     8.56608D-02, 8.47711D-02, 8.44170D-02, 8.32742D-02,
     $     8.28314D-02, 8.25659D-02, 8.10618D-02, 8.01623D-02,
     $     7.99321D-02, 7.89373D-02, 7.83328D-02, 7.81377D-02,
     $     7.71000D-02, 7.66254D-02, 7.60391D-02, 7.55238D-02,
     $     7.48650D-02, 7.45798D-02, 7.40056D-02, 7.26884D-02,
     $     7.18693D-02, 7.15192D-02, 7.05522D-02, 7.02589D-02,
     $     6.91876D-02, 6.87856D-02, 6.76824D-02, 6.75860D-02,
     $     6.69955D-02, 6.65240D-02, 6.60673D-02, 6.55632D-02,
     $     6.52911D-02, 6.47581D-02, 6.44160D-02, 6.38943D-02,
     $     6.32186D-02, 6.27380D-02, 6.26203D-02, 6.22768D-02,
     $     6.15141D-02, 6.14445D-02, 6.08982D-02, 6.06306D-02,
     $     6.01700D-02, 5.90325D-02, 5.87354D-02, 5.84885D-02,
     $     5.80004D-02, 5.78846D-02, 5.76859D-02, 5.72869D-02,
     $     5.70347D-02, 5.66600D-02, 5.60908D-02, 5.59725D-02,
     $     5.54712D-02, 5.50252D-02, 5.47862D-02, 5.47069D-02,
     $     5.42403D-02, 5.35222D-02, 5.31693D-02, 5.29871D-02,
     $     5.28522D-02, 5.27201D-02, 5.22188D-02, 5.18071D-02,
     $     5.17822D-02, 5.12138D-02, 5.10455D-02, 5.08336D-02,
     $     5.03434D-02, 5.01193D-02, 4.99393D-02 /

!     Branching ratio of 3009
      data( br(14,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00 /

!     Branching ratio of 4008
      data( br(15,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 3.37986D-06, 1.71228D-05, 1.89443D-05,
     $     3.99193D-04, 1.16936D-03, 6.08938D-03, 7.51764D-03,
     $     1.70720D-02, 2.20644D-02, 4.07925D-02, 5.21723D-02,
     $     5.36556D-02, 6.32226D-02, 7.13027D-02, 7.62222D-02,
     $     7.79762D-02, 8.47965D-02, 8.82122D-02, 9.93777D-02,
     $     9.98458D-02, 1.04392D-01, 1.11030D-01, 1.17952D-01,
     $     1.21448D-01, 1.22842D-01, 1.29358D-01, 1.30428D-01,
     $     1.31113D-01, 1.35813D-01, 1.36930D-01, 1.37871D-01,
     $     1.40792D-01, 1.45192D-01, 1.45866D-01, 1.47901D-01,
     $     1.48482D-01, 1.51571D-01, 1.53159D-01, 1.54627D-01,
     $     1.55201D-01, 1.56712D-01, 1.57155D-01, 1.58327D-01,
     $     1.58995D-01, 1.61260D-01, 1.62064D-01, 1.62330D-01,
     $     1.63402D-01, 1.64501D-01, 1.64727D-01, 1.64868D-01,
     $     1.65367D-01, 1.65461D-01, 1.65762D-01, 1.66756D-01,
     $     1.66861D-01, 1.67024D-01, 1.67432D-01, 1.67564D-01,
     $     1.67882D-01, 1.67816D-01, 1.67978D-01, 1.68184D-01,
     $     1.68066D-01, 1.67978D-01, 1.67939D-01, 1.67887D-01,
     $     1.67827D-01, 1.67759D-01, 1.67722D-01, 1.67529D-01,
     $     1.67488D-01, 1.67453D-01, 1.67429D-01, 1.67358D-01,
     $     1.67323D-01, 1.67117D-01, 1.67039D-01, 1.66996D-01,
     $     1.66790D-01, 1.66650D-01, 1.66524D-01, 1.66300D-01,
     $     1.66036D-01, 1.65918D-01, 1.65721D-01, 1.65216D-01,
     $     1.64821D-01, 1.64677D-01, 1.64177D-01, 1.64012D-01,
     $     1.63417D-01, 1.63200D-01, 1.62629D-01, 1.62587D-01,
     $     1.62263D-01, 1.62029D-01, 1.61763D-01, 1.61411D-01,
     $     1.61235D-01, 1.60887D-01, 1.60662D-01, 1.60357D-01,
     $     1.59955D-01, 1.59620D-01, 1.59532D-01, 1.59309D-01,
     $     1.58814D-01, 1.58765D-01, 1.58370D-01, 1.58189D-01,
     $     1.57881D-01, 1.57127D-01, 1.56921D-01, 1.56749D-01,
     $     1.56368D-01, 1.56271D-01, 1.56104D-01, 1.55788D-01,
     $     1.55574D-01, 1.55257D-01, 1.54780D-01, 1.54681D-01,
     $     1.54238D-01, 1.53831D-01, 1.53608D-01, 1.53534D-01,
     $     1.53090D-01, 1.52426D-01, 1.52132D-01, 1.51967D-01,
     $     1.51840D-01, 1.51722D-01, 1.51247D-01, 1.50864D-01,
     $     1.50840D-01, 1.50298D-01, 1.50129D-01, 1.49913D-01,
     $     1.49420D-01, 1.49180D-01, 1.48987D-01 /

!     Branching ratio of 3008
      data( br(16,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 6.49015D-91, 5.20716D-90,
     $     1.70547D-46, 1.03712D-21, 3.06385D-21, 1.23870D-19,
     $     1.55082D-18, 1.59833D-16, 1.15451D-15, 1.18881D-14,
     $     1.98603D-13, 7.76109D-13, 1.09363D-12, 2.28437D-12,
     $     2.15186D-11, 2.67954D-11, 1.36343D-10, 2.49108D-10,
     $     6.21260D-10, 3.32737D-09, 4.68492D-09, 6.08656D-09,
     $     9.89749D-09, 1.10474D-08, 1.33244D-08, 1.89930D-08,
     $     2.38060D-08, 3.30317D-08, 5.32093D-08, 5.88238D-08,
     $     9.10263D-08, 1.33541D-07, 1.62833D-07, 1.73766D-07,
     $     2.52836D-07, 4.15500D-07, 5.12074D-07, 5.71047D-07,
     $     6.19448D-07, 6.68577D-07, 8.87760D-07, 1.09357D-06,
     $     1.10836D-06, 1.51321D-06, 1.66897D-06, 1.88938D-06,
     $     2.50787D-06, 2.83850D-06, 3.14595D-06 /

!     Branching ratio of 4007
      data( br(17,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 1.45064D-42,
     $     1.60380D-28, 1.76078D-05, 8.70270D-05, 1.40045D-04,
     $     6.85700D-04, 1.31699D-03, 1.52745D-03, 1.68700D-03,
     $     2.41933D-03, 2.57852D-03, 2.86665D-03, 3.84585D-03,
     $     3.96139D-03, 4.17578D-03, 4.63652D-03, 4.86258D-03,
     $     5.25182D-03, 5.60612D-03, 5.90777D-03, 6.59322D-03,
     $     7.06523D-03, 7.38598D-03, 7.53118D-03, 7.71916D-03,
     $     8.13871D-03, 8.41789D-03, 8.52281D-03, 8.80311D-03,
     $     8.90452D-03, 8.96538D-03, 9.24766D-03, 9.40872D-03,
     $     9.45272D-03, 9.65222D-03, 9.76915D-03, 9.81134D-03,
     $     1.00359D-02, 1.01510D-02, 1.02858D-02, 1.04057D-02,
     $     1.05816D-02, 1.06634D-02, 1.08229D-02, 1.12158D-02,
     $     1.14779D-02, 1.15798D-02, 1.18426D-02, 1.19278D-02,
     $     1.22256D-02, 1.23305D-02, 1.25939D-02, 1.26140D-02,
     $     1.27513D-02, 1.28531D-02, 1.29638D-02, 1.30802D-02,
     $     1.31505D-02, 1.32960D-02, 1.33928D-02, 1.35307D-02,
     $     1.37114D-02, 1.38409D-02, 1.38716D-02, 1.39573D-02,
     $     1.41486D-02, 1.41668D-02, 1.43100D-02, 1.43731D-02,
     $     1.44745D-02, 1.46881D-02, 1.47379D-02, 1.47763D-02,
     $     1.48497D-02, 1.48661D-02, 1.48949D-02, 1.49439D-02,
     $     1.49748D-02, 1.50157D-02, 1.50574D-02, 1.50657D-02,
     $     1.51000D-02, 1.51253D-02, 1.51362D-02, 1.51395D-02,
     $     1.51547D-02, 1.51541D-02, 1.51454D-02, 1.51429D-02,
     $     1.51416D-02, 1.51391D-02, 1.51323D-02, 1.51168D-02,
     $     1.51162D-02, 1.51065D-02, 1.51062D-02, 1.51069D-02,
     $     1.51087D-02, 1.51084D-02, 1.51109D-02 /

!     Branching ratio of 3007
      data( br(18,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 1.58354D-08, 5.34941D-08, 5.04416D-07,
     $     2.46390D-06, 4.48010D-05, 1.11391D-04, 1.51172D-04,
     $     5.03399D-04, 9.13929D-04, 1.05079D-03, 1.15564D-03,
     $     1.64397D-03, 1.75006D-03, 1.94540D-03, 2.65215D-03,
     $     2.74770D-03, 2.93443D-03, 3.39043D-03, 3.64273D-03,
     $     4.10783D-03, 4.56551D-03, 4.97632D-03, 5.96179D-03,
     $     6.70877D-03, 7.26311D-03, 7.52853D-03, 7.89069D-03,
     $     8.79611D-03, 9.47174D-03, 9.74669D-03, 1.05611D-02,
     $     1.08810D-02, 1.10793D-02, 1.21190D-02, 1.27643D-02,
     $     1.29382D-02, 1.37198D-02, 1.41736D-02, 1.43298D-02,
     $     1.51426D-02, 1.55420D-02, 1.60117D-02, 1.64166D-02,
     $     1.69818D-02, 1.72348D-02, 1.77201D-02, 1.88430D-02,
     $     1.95657D-02, 1.98498D-02, 2.05992D-02, 2.08438D-02,
     $     2.17275D-02, 2.20510D-02, 2.29030D-02, 2.29708D-02,
     $     2.34259D-02, 2.37638D-02, 2.41095D-02, 2.44566D-02,
     $     2.46546D-02, 2.50408D-02, 2.52849D-02, 2.56244D-02,
     $     2.60602D-02, 2.63699D-02, 2.64435D-02, 2.66506D-02,
     $     2.71198D-02, 2.71653D-02, 2.75297D-02, 2.76962D-02,
     $     2.79777D-02, 2.86536D-02, 2.88315D-02, 2.89785D-02,
     $     2.92787D-02, 2.93503D-02, 2.94788D-02, 2.97244D-02,
     $     2.98912D-02, 3.01373D-02, 3.04741D-02, 3.05472D-02,
     $     3.08742D-02, 3.11701D-02, 3.13276D-02, 3.13797D-02,
     $     3.16851D-02, 3.20941D-02, 3.22676D-02, 3.23666D-02,
     $     3.24423D-02, 3.25113D-02, 3.27782D-02, 3.29468D-02,
     $     3.29583D-02, 3.32191D-02, 3.33010D-02, 3.34043D-02,
     $     3.36276D-02, 3.37184D-02, 3.37984D-02 /

!     Branching ratio of 3006
      data( br(19,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 5.93570D-19, 5.55023D-08, 2.29011D-06,
     $     2.43118D-05, 3.60235D-04, 1.66045D-03, 1.91342D-03,
     $     8.26447D-03, 1.42003D-02, 1.67931D-02, 1.89410D-02,
     $     2.11803D-02, 2.41959D-02, 3.52166D-02, 3.85826D-02,
     $     3.90867D-02, 4.28081D-02, 4.44762D-02, 4.58034D-02,
     $     4.61107D-02, 4.68703D-02, 4.68245D-02, 4.66082D-02,
     $     4.65817D-02, 4.63295D-02, 4.58789D-02, 4.54002D-02,
     $     4.51295D-02, 4.50010D-02, 4.42919D-02, 4.41503D-02,
     $     4.40547D-02, 4.33906D-02, 4.31896D-02, 4.30061D-02,
     $     4.23495D-02, 4.13241D-02, 4.11268D-02, 4.04560D-02,
     $     4.02447D-02, 3.92021D-02, 3.86638D-02, 3.79544D-02,
     $     3.75884D-02, 3.66498D-02, 3.64495D-02, 3.59372D-02,
     $     3.54731D-02, 3.42692D-02, 3.37462D-02, 3.35416D-02,
     $     3.25020D-02, 3.16708D-02, 3.14412D-02, 3.12754D-02,
     $     3.05652D-02, 3.04174D-02, 3.01295D-02, 2.91466D-02,
     $     2.90292D-02, 2.88144D-02, 2.83341D-02, 2.81094D-02,
     $     2.77177D-02, 2.73578D-02, 2.70394D-02, 2.64202D-02,
     $     2.60068D-02, 2.57296D-02, 2.56027D-02, 2.54340D-02,
     $     2.50102D-02, 2.47050D-02, 2.45837D-02, 2.41975D-02,
     $     2.40482D-02, 2.39585D-02, 2.34561D-02, 2.31554D-02,
     $     2.30781D-02, 2.27430D-02, 2.25410D-02, 2.24753D-02,
     $     2.21277D-02, 2.19679D-02, 2.17723D-02, 2.16013D-02,
     $     2.13808D-02, 2.12853D-02, 2.10945D-02, 2.06570D-02,
     $     2.03845D-02, 2.02692D-02, 1.99525D-02, 1.98557D-02,
     $     1.95026D-02, 1.93704D-02, 1.90094D-02, 1.89781D-02,
     $     1.87852D-02, 1.86323D-02, 1.84831D-02, 1.83197D-02,
     $     1.82312D-02, 1.80579D-02, 1.79468D-02, 1.77792D-02,
     $     1.75629D-02, 1.74096D-02, 1.73721D-02, 1.72634D-02,
     $     1.70230D-02, 1.70011D-02, 1.68292D-02, 1.67457D-02,
     $     1.66024D-02, 1.62505D-02, 1.61588D-02, 1.60828D-02,
     $     1.59324D-02, 1.58967D-02, 1.58354D-02, 1.57129D-02,
     $     1.56351D-02, 1.55198D-02, 1.53457D-02, 1.53095D-02,
     $     1.51555D-02, 1.50184D-02, 1.49450D-02, 1.49206D-02,
     $     1.47770D-02, 1.45575D-02, 1.44504D-02, 1.43948D-02,
     $     1.43536D-02, 1.43134D-02, 1.41608D-02, 1.40370D-02,
     $     1.40295D-02, 1.38581D-02, 1.38073D-02, 1.37434D-02,
     $     1.35964D-02, 1.35296D-02, 1.34760D-02 /

!     Branching ratio of 2006
      data( br(20,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00 /

!     Branching ratio of 2005
      data( br(21,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00 /

!     Branching ratio of 2004
      data( br(22,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     4.32042D-12, 1.44430D-08, 5.60277D-07, 1.46159D-05,
     $     1.70965D-04, 5.85101D-04, 5.11009D-03, 8.26237D-03,
     $     8.79494D-03, 1.43917D-02, 1.98383D-02, 2.58117D-02,
     $     2.82683D-02, 3.83942D-02, 4.06809D-02, 5.38969D-02,
     $     5.46231D-02, 6.19483D-02, 6.86347D-02, 7.79886D-02,
     $     8.32477D-02, 8.51287D-02, 9.42318D-02, 9.60982D-02,
     $     9.73242D-02, 1.03947D-01, 1.05915D-01, 1.07615D-01,
     $     1.13086D-01, 1.20177D-01, 1.21405D-01, 1.25313D-01,
     $     1.26501D-01, 1.32070D-01, 1.34854D-01, 1.38121D-01,
     $     1.39680D-01, 1.44028D-01, 1.45013D-01, 1.47412D-01,
     $     1.49476D-01, 1.55686D-01, 1.58499D-01, 1.59606D-01,
     $     1.65168D-01, 1.69650D-01, 1.70873D-01, 1.71755D-01,
     $     1.75526D-01, 1.76318D-01, 1.77927D-01, 1.83521D-01,
     $     1.84200D-01, 1.85441D-01, 1.88284D-01, 1.89603D-01,
     $     1.91876D-01, 1.93533D-01, 1.95174D-01, 1.98614D-01,
     $     2.00699D-01, 2.02177D-01, 2.02879D-01, 2.03832D-01,
     $     2.06345D-01, 2.08200D-01, 2.08950D-01, 2.11141D-01,
     $     2.12063D-01, 2.12626D-01, 2.15941D-01, 2.17960D-01,
     $     2.18471D-01, 2.20659D-01, 2.21976D-01, 2.22397D-01,
     $     2.24644D-01, 2.25640D-01, 2.26862D-01, 2.27742D-01,
     $     2.28984D-01, 2.29524D-01, 2.30640D-01, 2.33126D-01,
     $     2.34599D-01, 2.35252D-01, 2.36850D-01, 2.37357D-01,
     $     2.39267D-01, 2.40005D-01, 2.42116D-01, 2.42305D-01,
     $     2.43452D-01, 2.44411D-01, 2.45329D-01, 2.46184D-01,
     $     2.46712D-01, 2.47768D-01, 2.48454D-01, 2.49541D-01,
     $     2.50938D-01, 2.51845D-01, 2.52050D-01, 2.52715D-01,
     $     2.54262D-01, 2.54403D-01, 2.55507D-01, 2.56062D-01,
     $     2.57025D-01, 2.59362D-01, 2.59974D-01, 2.60475D-01,
     $     2.61380D-01, 2.61578D-01, 2.61938D-01, 2.62708D-01,
     $     2.63194D-01, 2.63925D-01, 2.64925D-01, 2.65149D-01,
     $     2.66108D-01, 2.66977D-01, 2.67443D-01, 2.67597D-01,
     $     2.68501D-01, 2.69803D-01, 2.70483D-01, 2.70830D-01,
     $     2.71087D-01, 2.71345D-01, 2.72298D-01, 2.72944D-01,
     $     2.72987D-01, 2.74039D-01, 2.74355D-01, 2.74754D-01,
     $     2.75696D-01, 2.76038D-01, 2.76356D-01 /

!     Branching ratio of 2003
      data( br(23,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 1.09365D-13,
     $     1.62272D-12, 5.59022D-12, 9.90885D-12, 1.77607D-11,
     $     6.25488D-11, 1.47084D-10, 2.00480D-10, 1.77762D-09,
     $     3.77095D-09, 5.76514D-09, 4.65706D-07, 2.92405D-06,
     $     4.45597D-06, 1.65685D-05, 2.90782D-05, 3.42568D-05,
     $     6.92151D-05, 9.03371D-05, 1.18294D-04, 1.45119D-04,
     $     1.85524D-04, 2.04839D-04, 2.44016D-04, 3.44212D-04,
     $     4.12821D-04, 4.40796D-04, 5.15313D-04, 5.39161D-04,
     $     6.23551D-04, 6.53846D-04, 7.33017D-04, 7.39465D-04,
     $     7.82674D-04, 8.17340D-04, 8.55554D-04, 8.99619D-04,
     $     9.26692D-04, 9.85621D-04, 1.02713D-03, 1.09194D-03,
     $     1.18346D-03, 1.25293D-03, 1.27012D-03, 1.31993D-03,
     $     1.43502D-03, 1.44601D-03, 1.53398D-03, 1.57542D-03,
     $     1.64498D-03, 1.80604D-03, 1.84638D-03, 1.87873D-03,
     $     1.94217D-03, 1.95685D-03, 1.98244D-03, 2.03074D-03,
     $     2.06165D-03, 2.10636D-03, 2.16792D-03, 2.18092D-03,
     $     2.23724D-03, 2.28701D-03, 2.31328D-03, 2.32195D-03,
     $     2.37233D-03, 2.44169D-03, 2.47290D-03, 2.49050D-03,
     $     2.50404D-03, 2.51684D-03, 2.56840D-03, 2.60739D-03,
     $     2.61004D-03, 2.67484D-03, 2.69650D-03, 2.72518D-03,
     $     2.79397D-03, 2.82596D-03, 2.85403D-03 /

!     Branching ratio of 1003
      data( br(24,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     7.18794D-30, 1.15537D-15, 4.65826D-13, 1.33023D-11,
     $     5.84585D-11, 1.01925D-10, 3.45963D-10, 4.90860D-10,
     $     1.87907D-09, 3.38182D-09, 1.68915D-08, 1.90973D-08,
     $     4.64292D-08, 8.54603D-08, 1.59295D-07, 2.84857D-07,
     $     3.95217D-07, 7.07087D-07, 1.01576D-06, 1.65831D-06,
     $     2.98705D-06, 4.40633D-06, 4.81837D-06, 6.18974D-06,
     $     1.04025D-05, 1.08949D-05, 1.54552D-05, 1.80333D-05,
     $     2.31497D-05, 3.98873D-05, 4.54404D-05, 5.04493D-05,
     $     6.19552D-05, 6.49711D-05, 7.05427D-05, 8.25045D-05,
     $     9.11980D-05, 1.05311D-04, 1.28222D-04, 1.33548D-04,
     $     1.58857D-04, 1.84435D-04, 1.99254D-04, 2.04352D-04,
     $     2.36264D-04, 2.87600D-04, 3.13454D-04, 3.28151D-04,
     $     3.39526D-04, 3.50510D-04, 3.94810D-04, 4.29633D-04,
     $     4.31941D-04, 4.87067D-04, 5.04893D-04, 5.28177D-04,
     $     5.83803D-04, 6.09437D-04, 6.31557D-04 /

!     Branching ratio of 1002
      data( br(25,k), k=1, 150) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 1.35547D-31, 1.78857D-24, 1.78514D-18,
     $     2.19415D-15, 2.47715D-14, 6.83350D-13, 1.85974D-12,
     $     2.12381D-12, 5.87712D-12, 1.39149D-11, 2.90079D-11,
     $     3.77778D-11, 9.66345D-11, 1.15204D-10, 2.66587D-10,
     $     2.77704D-10, 4.13429D-10, 5.62002D-10, 1.00426D-09,
     $     4.09450D-09, 1.01712D-08, 2.79662D-07, 4.76324D-07,
     $     6.75131D-07, 2.86548D-06, 4.08111D-06, 5.34693D-06,
     $     1.10090D-05, 2.09547D-05, 2.32167D-05, 3.17398D-05,
     $     3.47993D-05, 5.12347D-05, 6.12884D-05, 7.73560D-05,
     $     8.69499D-05, 1.21351D-04, 1.29993D-04, 1.53922D-04,
     $     1.81949D-04, 2.80565D-04, 3.32525D-04, 3.53292D-04,
     $     4.55993D-04, 5.23676D-04, 5.43409D-04, 5.58186D-04,
     $     6.23992D-04, 6.38022D-04, 6.63365D-04, 7.55727D-04,
     $     7.67474D-04, 7.89627D-04, 8.37681D-04, 8.61473D-04,
     $     9.00727D-04, 9.35368D-04, 9.63999D-04, 1.03032D-03,
     $     1.07881D-03, 1.11455D-03, 1.13172D-03, 1.15534D-03,
     $     1.21438D-03, 1.25895D-03, 1.27735D-03, 1.33135D-03,
     $     1.35257D-03, 1.36580D-03, 1.43337D-03, 1.47374D-03,
     $     1.48448D-03, 1.53132D-03, 1.55678D-03, 1.56547D-03,
     $     1.60886D-03, 1.62947D-03, 1.65260D-03, 1.67187D-03,
     $     1.69835D-03, 1.70981D-03, 1.73076D-03, 1.77608D-03,
     $     1.80404D-03, 1.81454D-03, 1.84222D-03, 1.85161D-03,
     $     1.88607D-03, 1.89901D-03, 1.93394D-03, 1.93675D-03,
     $     1.95659D-03, 1.97167D-03, 1.98789D-03, 2.00429D-03,
     $     2.01408D-03, 2.03370D-03, 2.04636D-03, 2.06408D-03,
     $     2.08754D-03, 2.10443D-03, 2.10840D-03, 2.11950D-03,
     $     2.14542D-03, 2.14796D-03, 2.16828D-03, 2.17751D-03,
     $     2.19311D-03, 2.22972D-03, 2.23930D-03, 2.24731D-03,
     $     2.26405D-03, 2.26807D-03, 2.27539D-03, 2.28956D-03,
     $     2.29950D-03, 2.31452D-03, 2.33545D-03, 2.34016D-03,
     $     2.36187D-03, 2.38253D-03, 2.39385D-03, 2.39762D-03,
     $     2.41993D-03, 2.45093D-03, 2.46426D-03, 2.47194D-03,
     $     2.47792D-03, 2.48346D-03, 2.50537D-03, 2.51980D-03,
     $     2.52080D-03, 2.54405D-03, 2.55158D-03, 2.56118D-03,
     $     2.58235D-03, 2.59130D-03, 2.59917D-03 /

      do j = 1, lines-1
         if ( u >= exen(j) .and. u < exen(j+1) ) then

            do i = 1, nimax
               if ( ifz(i)*1000+ ifa(i) == 1. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(1,j), exen(j+1),
     $                    br(1,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 1001. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(2,j), exen(j+1),
     $                    br(2,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 7012. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(3,j), exen(j+1),
     $                    br(3,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 1002. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(4,j), exen(j+1),
     $                    br(4,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 5012. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(5,j), exen(j+1),
     $                    br(5,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 1003. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(6,j), exen(j+1),
     $                    br(6,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 2003. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(7,j), exen(j+1),
     $                    br(7,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 4011. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(8,j), exen(j+1),
     $                    br(8,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 6010. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(9,j), exen(j+1),
     $                    br(9,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 2004. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(10,j), exen(j+1),
     $                    br(10,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 4010. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(11,j), exen(j+1),
     $                    br(11,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 4009. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(13,j), exen(j+1),
     $                    br(13,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 3009. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(14,j), exen(j+1),
     $                    br(14,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 3006. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(15,j), exen(j+1),
     $                    br(15,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 3008. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(16,j), exen(j+1),
     $                    br(16,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 3007. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(17,j), exen(j+1),
     $                    br(17,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 4007. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(18,j), exen(j+1),
     $                    br(18,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 4008. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(19,j), exen(j+1),
     $                    br(19,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 2006. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(20,j), exen(j+1),
     $                    br(20,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 5010. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(22,j), exen(j+1),
     $                    br(22,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 5011. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(23,j), exen(j+1),
     $                    br(23,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 6011. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(24,j), exen(j+1),
     $                    br(24,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 6012. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(25,j), exen(j+1),
     $                    br(25,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               end if

            end do

            do i = 1, nimax
               r(i) = r(i)/ integral
            end do

         end if

      end do

      end subroutine suppressalpha_n14


! ----------------------------------------------------------------------
!     Modify the decay width of O-16 in GDR
      subroutine suppressalpha_o16( u)
!     u: incident photon energy (MeV)
! ----------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

      common /exiejn/ nimax
!$OMP THREADPRIVATE(/exiejn/)
      common /std1/ r(70),s(70),sigma,rr(70)
!$OMP THREADPRIVATE(/std1/)
      common /ejectl/ omega(70),ifa(70),ifz(70)

!     implicit none
      double precision u
      double precision sigma

      double precision exen(500)  ! excitation energy (MeV)
      double precision br(30,500) ! branching ration (-)
      integer lines               ! number of energy bins

      integer i, j, k           ! counters
      double precision integral

      double precision linearone

      integral = 0.

!      8016
!  (1) 8015     1
!  (2) 7015  1001
!  (3) 6015  2001
!  (4) 8014     2
!  (5) 7014  1002
!  (6) 6014  2002
!  (7) 7013  1003
!  (8) 6013  2003
!  (9) 6012  2004
! (10) 5012  3004
! (11) 4012  4004
! (12) 6011  2005
! (13) 5011  3005
! (14) 4011  4005
! (15) 6010  2006
! (16) 5010  3006
! (17) 4010  4006
! (18) 5009  3007
! (19) 4009  4007
! (20) 3009  5007
! (21) 4008  4008
! (22) 3008  5008
! (23) 3007  5009
! (24) 3006  5010
! (25) 2006  6010
! (26) 2005  6011
! (27) 2004  6012
! (28) 2003  6013
! (29) 1003  7013
! (30) 1002  7014

      lines = 300
      data( exen(k), k=1, 300) /
     $     14.07900, 19.48900, 22.19300, 22.71600,
     $     24.43600, 24.64000, 25.65000, 26.14000, 27.06300,
     $     27.43700, 28.12200, 28.27500, 28.89700, 29.17200,
     $     29.72500, 29.94700, 30.14600, 30.28500, 30.67600,
     $     31.28900, 31.34500, 31.53800, 31.83500, 31.96800,
     $     32.03800, 32.26200, 32.48200, 32.91300, 33.14000,
     $     33.24300, 33.42100, 33.52500, 33.69200, 33.91700,
     $     34.06400, 34.26600, 34.32500, 34.43400, 34.58300,
     $     34.74300, 34.86000, 35.00600, 35.11700, 35.19300,
     $     35.30000, 35.33100, 35.52100, 35.66500, 35.82000,
     $     35.95000, 35.97200, 36.07400, 36.21400, 36.24100,
     $     36.36600, 36.42000, 36.52700, 36.58600, 36.65000,
     $     36.85700, 36.87800, 36.90400, 36.98500, 37.09200,
     $     37.17400, 37.24800, 37.38500, 37.58600, 37.65000,
     $     37.76100, 37.77800, 37.85400, 37.90500, 37.96200,
     $     38.08500, 38.17700, 38.26700, 38.29800, 38.30800,
     $     38.41000, 38.54600, 38.62600, 38.70600, 38.78100,
     $     38.82700, 38.89500, 38.93000, 39.07800, 39.11500,
     $     39.17400, 39.30200, 39.36000, 39.39000, 39.47300,
     $     39.51400, 39.61300, 39.69300, 39.70300, 39.77500,
     $     39.83500, 39.90500, 39.97500, 40.07900, 40.14900,
     $     40.22400, 40.27100, 40.31100, 40.37800, 40.47600,
     $     40.52000, 40.62000, 40.70600, 40.74700, 40.81600,
     $     40.86700, 40.88700, 40.92000, 40.95800, 41.03800,
     $     41.10500, 41.11300, 41.16200, 41.24300, 41.33600,
     $     41.38400, 41.44400, 41.52400, 41.59800, 41.61800,
     $     41.69400, 41.71100, 41.80100, 41.82700, 41.90700,
     $     41.98200, 42.01300, 42.08400, 42.10900, 42.16000,
     $     42.20400, 42.28100, 42.34200, 42.38900, 42.48400,
     $     42.49000, 42.61800, 42.63500, 42.67400, 42.73800,
     $     42.82700, 42.87600, 42.91800, 42.95000, 42.96800,
     $     43.01800, 43.06400, 43.12400, 43.16200, 43.19400,
     $     43.30900, 43.33600, 43.36600, 43.40600, 43.44700,
     $     43.52200, 43.56400, 43.59000, 43.66800, 43.68600,
     $     43.73100, 43.79100, 43.85300, 43.89100, 43.92500,
     $     44.01300, 44.04700, 44.07100, 44.11500, 44.21000,
     $     44.25800, 44.28700, 44.30900, 44.38800, 44.40300,
     $     44.44500, 44.50400, 44.51200, 44.56600, 44.63200,
     $     44.67300, 44.72500, 44.78300, 44.82500, 44.85000,
     $     44.85800, 44.89700, 44.93700, 44.98700, 45.01100,
     $     45.05200, 45.07500, 45.12900, 45.16200, 45.22800,
     $     45.23500, 45.28000, 45.33700, 45.38400, 45.44900,
     $     45.47600, 45.51400, 45.55900, 45.60000, 45.61500,
     $     45.67600, 45.68300, 45.75100, 45.77000, 45.82600,
     $     45.86500, 45.90200, 45.92700, 45.98000, 45.98900,
     $     46.10800, 46.16100, 46.17300, 46.18600, 46.26900,
     $     46.31300, 46.34300, 46.37700, 46.42000, 46.42600,
     $     46.47700, 46.49600, 46.53600, 46.55500, 46.63400,
     $     46.65800, 46.73400, 46.74000, 46.79200, 46.83400,
     $     46.86100, 46.92500, 46.93700, 46.97700, 47.00800,
     $     47.02600, 47.08900, 47.11900, 47.15600, 47.18200,
     $     47.21100, 47.21900, 47.27900, 47.29800, 47.32600,
     $     47.35800, 47.38000, 47.45800, 47.50900, 47.51500,
     $     47.57100, 47.60800, 47.66400, 47.67500, 47.71400,
     $     47.73700, 47.75600, 47.78400, 47.83300, 47.85100,
     $     47.90800, 47.93100, 47.97300, 48.03600, 48.05800,
     $     48.07900, 48.13700, 48.15400, 48.18200, 48.19300,
     $     48.24000, 48.25100, 48.32600, 48.38400, 48.41900,
     $     48.47400, 48.49200, 48.50800, 48.54900, 48.57600,
     $     48.63500, 48.65000, 48.67800, 48.71500, 48.72800,
     $     48.79600 /

!     Branching ratio of 8015
      data( br(1,k), k=1, 300) /
     $     0.00000D+00, 4.21168D-01, 3.09161D-01, 3.46058D-01,
     $     3.83506D-01, 3.73090D-01, 2.94403D-01, 2.65373D-01,
     $     2.06743D-01, 1.85704D-01, 1.51490D-01, 1.44035D-01,
     $     1.19504D-01, 1.10427D-01, 9.54070D-02, 9.06953D-02,
     $     8.67493D-02, 8.41049D-02, 7.71237D-02, 6.70304D-02,
     $     6.61073D-02, 6.31707D-02, 5.89048D-02, 5.69556D-02,
     $     5.59986D-02, 5.30969D-02, 5.04746D-02, 4.57882D-02,
     $     4.35473D-02, 4.25691D-02, 4.09429D-02, 4.00495D-02,
     $     3.85180D-02, 3.67092D-02, 3.55873D-02, 3.40939D-02,
     $     3.36687D-02, 3.28712D-02, 3.18200D-02, 3.07349D-02,
     $     2.99632D-02, 2.89970D-02, 2.82739D-02, 2.78037D-02,
     $     2.71483D-02, 2.69666D-02, 2.58793D-02, 2.50989D-02,
     $     2.42866D-02, 2.36426D-02, 2.35364D-02, 2.30298D-02,
     $     2.23531D-02, 2.22188D-02, 2.16436D-02, 2.13983D-02,
     $     2.09121D-02, 2.06515D-02, 2.03716D-02, 1.94729D-02,
     $     1.93848D-02, 1.92764D-02, 1.89436D-02, 1.84995D-02,
     $     1.81683D-02, 1.78792D-02, 1.73624D-02, 1.66490D-02,
     $     1.64251D-02, 1.60489D-02, 1.59926D-02, 1.57329D-02,
     $     1.55633D-02, 1.53780D-02, 1.49889D-02, 1.46965D-02,
     $     1.44185D-02, 1.43251D-02, 1.42951D-02, 1.39888D-02,
     $     1.35873D-02, 1.33599D-02, 1.31342D-02, 1.29300D-02,
     $     1.28030D-02, 1.26221D-02, 1.25305D-02, 1.21402D-02,
     $     1.20461D-02, 1.18976D-02, 1.15827D-02, 1.14431D-02,
     $     1.13715D-02, 1.11690D-02, 1.10721D-02, 1.08384D-02,
     $     1.06537D-02, 1.06311D-02, 1.04683D-02, 1.03354D-02,
     $     1.01789D-02, 1.00253D-02, 9.80612D-03, 9.66155D-03,
     $     9.51076D-03, 9.41859D-03, 9.33937D-03, 9.20961D-03,
     $     9.02675D-03, 8.94685D-03, 8.76352D-03, 8.61224D-03,
     $     8.54192D-03, 8.42370D-03, 8.33882D-03, 8.30590D-03,
     $     8.25202D-03, 8.19059D-03, 8.06023D-03, 7.95438D-03,
     $     7.94192D-03, 7.86464D-03, 7.73860D-03, 7.59788D-03,
     $     7.52604D-03, 7.43854D-03, 7.32449D-03, 7.21995D-03,
     $     7.19223D-03, 7.08844D-03, 7.06551D-03, 6.94397D-03,
     $     6.90964D-03, 6.80551D-03, 6.70789D-03, 6.66753D-03,
     $     6.57741D-03, 6.54626D-03, 6.48340D-03, 6.42981D-03,
     $     6.33655D-03, 6.26191D-03, 6.20524D-03, 6.09387D-03,
     $     6.08691D-03, 5.94049D-03, 5.92107D-03, 5.87723D-03,
     $     5.80615D-03, 5.70983D-03, 5.65775D-03, 5.61286D-03,
     $     5.57892D-03, 5.55913D-03, 5.50595D-03, 5.45769D-03,
     $     5.39590D-03, 5.35729D-03, 5.32500D-03, 5.21052D-03,
     $     5.18397D-03, 5.15460D-03, 5.11549D-03, 5.07483D-03,
     $     5.00121D-03, 4.96024D-03, 4.93535D-03, 4.86159D-03,
     $     4.84460D-03, 4.80234D-03, 4.74668D-03, 4.68978D-03,
     $     4.65530D-03, 4.62482D-03, 4.54745D-03, 4.51802D-03,
     $     4.49738D-03, 4.45952D-03, 4.37850D-03, 4.33828D-03,
     $     4.31404D-03, 4.29582D-03, 4.23153D-03, 4.21931D-03,
     $     4.18549D-03, 4.13851D-03, 4.13218D-03, 4.08957D-03,
     $     4.03844D-03, 4.00704D-03, 3.96780D-03, 3.92407D-03,
     $     3.89289D-03, 3.87456D-03, 3.86873D-03, 3.84046D-03,
     $     3.81142D-03, 3.77543D-03, 3.75819D-03, 3.72917D-03,
     $     3.71305D-03, 3.67550D-03, 3.65262D-03, 3.60718D-03,
     $     3.60237D-03, 3.57188D-03, 3.53294D-03, 3.50166D-03,
     $     3.45850D-03, 3.44072D-03, 3.41603D-03, 3.38719D-03,
     $     3.36119D-03, 3.35174D-03, 3.31336D-03, 3.30901D-03,
     $     3.26694D-03, 3.25533D-03, 3.22142D-03, 3.19779D-03,
     $     3.17542D-03, 3.16054D-03, 3.12940D-03, 3.12415D-03,
     $     3.05531D-03, 3.02519D-03, 3.01836D-03, 3.01092D-03,
     $     2.96440D-03, 2.94025D-03, 2.92360D-03, 2.90514D-03,
     $     2.88180D-03, 2.87856D-03, 2.85095D-03, 2.84078D-03,
     $     2.81966D-03, 2.80972D-03, 2.76877D-03, 2.75640D-03,
     $     2.71771D-03, 2.71469D-03, 2.68873D-03, 2.66800D-03,
     $     2.65462D-03, 2.62331D-03, 2.61748D-03, 2.59815D-03,
     $     2.58335D-03, 2.57480D-03, 2.54502D-03, 2.53101D-03,
     $     2.51374D-03, 2.50166D-03, 2.48832D-03, 2.48466D-03,
     $     2.45736D-03, 2.44878D-03, 2.43612D-03, 2.42158D-03,
     $     2.41176D-03, 2.37719D-03, 2.35497D-03, 2.35233D-03,
     $     2.32806D-03, 2.31228D-03, 2.28853D-03, 2.28392D-03,
     $     2.26761D-03, 2.25801D-03, 2.25014D-03, 2.23863D-03,
     $     2.21866D-03, 2.21138D-03, 2.18852D-03, 2.17930D-03,
     $     2.16256D-03, 2.13758D-03, 2.12896D-03, 2.12080D-03,
     $     2.09831D-03, 2.09176D-03, 2.08102D-03, 2.07683D-03,
     $     2.05911D-03, 2.05495D-03, 2.02704D-03, 2.00572D-03,
     $     1.99297D-03, 1.97311D-03, 1.96669D-03, 1.96102D-03,
     $     1.94661D-03, 1.93712D-03, 1.91661D-03, 1.91144D-03,
     $     1.90183D-03, 1.88921D-03, 1.88479D-03, 1.86206D-03 /

!     Branching ratio of 7015
      data( br(2,k), k=1, 300) /
     $     9.99908D-01, 5.78173D-01, 6.89326D-01, 6.52324D-01,
     $     5.61943D-01, 5.40625D-01, 4.16515D-01, 3.75769D-01,
     $     2.96464D-01, 2.68253D-01, 2.21998D-01, 2.11774D-01,
     $     1.78072D-01, 1.65483D-01, 1.44589D-01, 1.38040D-01,
     $     1.32522D-01, 1.28809D-01, 1.18963D-01, 1.04476D-01,
     $     1.03130D-01, 9.88552D-02, 9.26121D-02, 8.97276D-02,
     $     8.83117D-02, 8.40087D-02, 8.01086D-02, 7.31012D-02,
     $     6.97307D-02, 6.82542D-02, 6.57930D-02, 6.44391D-02,
     $     6.20952D-02, 5.93338D-02, 5.76168D-02, 5.53228D-02,
     $     5.46678D-02, 5.34347D-02, 5.18060D-02, 5.01208D-02,
     $     4.89198D-02, 4.74110D-02, 4.62790D-02, 4.55427D-02,
     $     4.45142D-02, 4.42291D-02, 4.25194D-02, 4.12898D-02,
     $     4.00071D-02, 3.89899D-02, 3.88221D-02, 3.80196D-02,
     $     3.69454D-02, 3.67316D-02, 3.58172D-02, 3.54266D-02,
     $     3.46509D-02, 3.42345D-02, 3.37871D-02, 3.23475D-02,
     $     3.22062D-02, 3.20325D-02, 3.14986D-02, 3.07843D-02,
     $     3.02512D-02, 2.97856D-02, 2.89525D-02, 2.78004D-02,
     $     2.74380D-02, 2.68288D-02, 2.67376D-02, 2.63161D-02,
     $     2.60408D-02, 2.57401D-02, 2.51080D-02, 2.46319D-02,
     $     2.41791D-02, 2.40268D-02, 2.39779D-02, 2.34781D-02,
     $     2.28217D-02, 2.24495D-02, 2.20799D-02, 2.17455D-02,
     $     2.15373D-02, 2.12407D-02, 2.10904D-02, 2.04493D-02,
     $     2.02947D-02, 2.00508D-02, 1.95327D-02, 1.93027D-02,
     $     1.91848D-02, 1.88506D-02, 1.86907D-02, 1.83046D-02,
     $     1.79994D-02, 1.79620D-02, 1.76927D-02, 1.74730D-02,
     $     1.72138D-02, 1.69591D-02, 1.65959D-02, 1.63562D-02,
     $     1.61060D-02, 1.59532D-02, 1.58217D-02, 1.56063D-02,
     $     1.53027D-02, 1.51700D-02, 1.48651D-02, 1.46136D-02,
     $     1.44966D-02, 1.43000D-02, 1.41588D-02, 1.41040D-02,
     $     1.40144D-02, 1.39122D-02, 1.36952D-02, 1.35189D-02,
     $     1.34982D-02, 1.33694D-02, 1.31592D-02, 1.29245D-02,
     $     1.28046D-02, 1.26585D-02, 1.24681D-02, 1.22934D-02,
     $     1.22471D-02, 1.20736D-02, 1.20353D-02, 1.18320D-02,
     $     1.17746D-02, 1.16004D-02, 1.14370D-02, 1.13694D-02,
     $     1.12185D-02, 1.11663D-02, 1.10610D-02, 1.09712D-02,
     $     1.08150D-02, 1.06898D-02, 1.05947D-02, 1.04078D-02,
     $     1.03961D-02, 1.01504D-02, 1.01178D-02, 1.00442D-02,
     $     9.92478D-03, 9.76305D-03, 9.67558D-03, 9.60016D-03,
     $     9.54313D-03, 9.50983D-03, 9.42043D-03, 9.33927D-03,
     $     9.23538D-03, 9.17044D-03, 9.11614D-03, 8.92352D-03,
     $     8.87884D-03, 8.82938D-03, 8.76352D-03, 8.69500D-03,
     $     8.57091D-03, 8.50182D-03, 8.45986D-03, 8.33546D-03,
     $     8.30680D-03, 8.23547D-03, 8.14149D-03, 8.04541D-03,
     $     7.98716D-03, 7.93568D-03, 7.80498D-03, 7.75525D-03,
     $     7.72036D-03, 7.65637D-03, 7.51937D-03, 7.45133D-03,
     $     7.41033D-03, 7.37950D-03, 7.27072D-03, 7.25003D-03,
     $     7.19279D-03, 7.11324D-03, 7.10252D-03, 7.03034D-03,
     $     6.94371D-03, 6.89052D-03, 6.82404D-03, 6.74994D-03,
     $     6.69710D-03, 6.66604D-03, 6.65616D-03, 6.60825D-03,
     $     6.55902D-03, 6.49802D-03, 6.46879D-03, 6.41957D-03,
     $     6.39224D-03, 6.32857D-03, 6.28975D-03, 6.21263D-03,
     $     6.20446D-03, 6.15271D-03, 6.08656D-03, 6.03345D-03,
     $     5.96012D-03, 5.92990D-03, 5.88795D-03, 5.83895D-03,
     $     5.79477D-03, 5.77871D-03, 5.71347D-03, 5.70608D-03,
     $     5.63454D-03, 5.61480D-03, 5.55713D-03, 5.51693D-03,
     $     5.47884D-03, 5.45353D-03, 5.40053D-03, 5.39160D-03,
     $     5.27439D-03, 5.22309D-03, 5.21144D-03, 5.19877D-03,
     $     5.11950D-03, 5.07834D-03, 5.04995D-03, 5.01848D-03,
     $     4.97869D-03, 4.97315D-03, 4.92607D-03, 4.90871D-03,
     $     4.87268D-03, 4.85574D-03, 4.78587D-03, 4.76476D-03,
     $     4.69873D-03, 4.69358D-03, 4.64928D-03, 4.61390D-03,
     $     4.59105D-03, 4.53761D-03, 4.52767D-03, 4.49466D-03,
     $     4.46938D-03, 4.45479D-03, 4.40393D-03, 4.38001D-03,
     $     4.35050D-03, 4.32985D-03, 4.30706D-03, 4.30081D-03,
     $     4.25415D-03, 4.23947D-03, 4.21783D-03, 4.19297D-03,
     $     4.17617D-03, 4.11704D-03, 4.07903D-03, 4.07452D-03,
     $     4.03298D-03, 4.00598D-03, 3.96534D-03, 3.95743D-03,
     $     3.92952D-03, 3.91308D-03, 3.89962D-03, 3.87991D-03,
     $     3.84573D-03, 3.83326D-03, 3.79411D-03, 3.77831D-03,
     $     3.74964D-03, 3.70684D-03, 3.69206D-03, 3.67808D-03,
     $     3.63953D-03, 3.62829D-03, 3.60987D-03, 3.60268D-03,
     $     3.57231D-03, 3.56518D-03, 3.51730D-03, 3.48073D-03,
     $     3.45886D-03, 3.42478D-03, 3.41376D-03, 3.40403D-03,
     $     3.37930D-03, 3.36301D-03, 3.32781D-03, 3.31893D-03,
     $     3.30244D-03, 3.28078D-03, 3.27319D-03, 3.23417D-03 /

!     Branching ratio of 6015
      data( br(3,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00 /

!     Branching ratio of 8014
      data( br(4,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 3.65687D-12, 4.01249D-07, 5.70860D-07,
     $     7.16439D-07, 8.21663D-07, 1.07211D-06, 1.79938D-06,
     $     2.37832D-06, 4.31914D-06, 7.15065D-06, 9.49680D-06,
     $     1.02726D-05, 1.22112D-05, 1.36652D-05, 4.66276D-05,
     $     9.53705D-05, 1.29486D-04, 1.97976D-04, 2.31704D-04,
     $     4.72751D-04, 6.26191D-04, 6.98845D-04, 7.83599D-04,
     $     8.05603D-04, 8.52499D-04, 9.07574D-04, 9.57231D-04,
     $     9.89474D-04, 1.03880D-03, 1.07965D-03, 1.10097D-03,
     $     1.13395D-03, 1.14121D-03, 1.18778D-03, 1.21882D-03,
     $     1.30248D-03, 1.34571D-03, 1.35146D-03, 1.41302D-03,
     $     1.49577D-03, 1.51034D-03, 1.55240D-03, 1.56811D-03,
     $     1.59837D-03, 1.61103D-03, 1.62314D-03, 1.66749D-03,
     $     1.67157D-03, 1.67658D-03, 1.69257D-03, 1.74527D-03,
     $     1.77770D-03, 1.80329D-03, 1.84684D-03, 1.90474D-03,
     $     1.92525D-03, 1.95820D-03, 1.96288D-03, 1.98914D-03,
     $     2.00591D-03, 2.02311D-03, 2.05674D-03, 2.08320D-03,
     $     2.10601D-03, 2.11288D-03, 2.11501D-03, 2.13708D-03,
     $     2.16523D-03, 2.17892D-03, 2.19214D-03, 2.20261D-03,
     $     2.21079D-03, 2.21977D-03, 2.22383D-03, 2.24905D-03,
     $     2.25315D-03, 2.25905D-03, 2.27007D-03, 2.27413D-03,
     $     2.27600D-03, 2.28354D-03, 2.28594D-03, 2.29157D-03,
     $     2.29514D-03, 2.29539D-03, 2.29731D-03, 2.29829D-03,
     $     2.30107D-03, 2.30746D-03, 2.31030D-03, 2.31302D-03,
     $     2.31424D-03, 2.31393D-03, 2.31391D-03, 2.31318D-03,
     $     2.31048D-03, 2.30900D-03, 2.30958D-03, 2.30834D-03,
     $     2.30695D-03, 2.30514D-03, 2.30297D-03, 2.30205D-03,
     $     2.30045D-03, 2.29849D-03, 2.29505D-03, 2.29160D-03,
     $     2.29115D-03, 2.28923D-03, 2.28556D-03, 2.28038D-03,
     $     2.27749D-03, 2.27333D-03, 2.26724D-03, 2.26154D-03,
     $     2.25989D-03, 2.25335D-03, 2.25185D-03, 2.24419D-03,
     $     2.24188D-03, 2.23462D-03, 2.22933D-03, 2.22710D-03,
     $     2.22124D-03, 2.21899D-03, 2.21423D-03, 2.20993D-03,
     $     2.20242D-03, 2.19751D-03, 2.19497D-03, 2.18645D-03,
     $     2.18586D-03, 2.17296D-03, 2.17122D-03, 2.16709D-03,
     $     2.16026D-03, 2.15000D-03, 2.14410D-03, 2.13926D-03,
     $     2.13591D-03, 2.13431D-03, 2.12883D-03, 2.12397D-03,
     $     2.11699D-03, 2.11234D-03, 2.10832D-03, 2.09324D-03,
     $     2.08958D-03, 2.08547D-03, 2.08008D-03, 2.07500D-03,
     $     2.06475D-03, 2.05901D-03, 2.05534D-03, 2.04416D-03,
     $     2.04156D-03, 2.03572D-03, 2.02793D-03, 2.01948D-03,
     $     2.01426D-03, 2.00952D-03, 1.99700D-03, 1.99211D-03,
     $     1.98864D-03, 1.98235D-03, 1.96875D-03, 1.96180D-03,
     $     1.95758D-03, 1.95439D-03, 1.94288D-03, 1.94070D-03,
     $     1.93455D-03, 1.92589D-03, 1.92472D-03, 1.91713D-03,
     $     1.90753D-03, 1.90152D-03, 1.89386D-03, 1.88534D-03,
     $     1.87917D-03, 1.87548D-03, 1.87429D-03, 1.86852D-03,
     $     1.86268D-03, 1.85528D-03, 1.85177D-03, 1.84572D-03,
     $     1.84231D-03, 1.83425D-03, 1.82930D-03, 1.81982D-03,
     $     1.81888D-03, 1.81240D-03, 1.80426D-03, 1.79729D-03,
     $     1.78795D-03, 1.78402D-03, 1.77838D-03, 1.77164D-03,
     $     1.76544D-03, 1.76317D-03, 1.75390D-03, 1.75283D-03,
     $     1.74248D-03, 1.73958D-03, 1.73106D-03, 1.72545D-03,
     $     1.71999D-03, 1.71626D-03, 1.70829D-03, 1.70693D-03,
     $     1.68913D-03, 1.68115D-03, 1.67937D-03, 1.67745D-03,
     $     1.66512D-03, 1.65850D-03, 1.65426D-03, 1.64923D-03,
     $     1.64279D-03, 1.64190D-03, 1.63436D-03, 1.63154D-03,
     $     1.62554D-03, 1.62267D-03, 1.61070D-03, 1.60713D-03,
     $     1.59563D-03, 1.59472D-03, 1.58681D-03, 1.58040D-03,
     $     1.57629D-03, 1.56657D-03, 1.56474D-03, 1.55865D-03,
     $     1.55394D-03, 1.55121D-03, 1.54168D-03, 1.53714D-03,
     $     1.53165D-03, 1.52774D-03, 1.52336D-03, 1.52216D-03,
     $     1.51321D-03, 1.51041D-03, 1.50623D-03, 1.50154D-03,
     $     1.49828D-03, 1.48670D-03, 1.47911D-03, 1.47821D-03,
     $     1.47002D-03, 1.46453D-03, 1.45623D-03, 1.45459D-03,
     $     1.44878D-03, 1.44535D-03, 1.44251D-03, 1.43834D-03,
     $     1.43104D-03, 1.42836D-03, 1.41989D-03, 1.41646D-03,
     $     1.41021D-03, 1.40099D-03, 1.39777D-03, 1.39470D-03,
     $     1.38629D-03, 1.38385D-03, 1.37980D-03, 1.37820D-03,
     $     1.37139D-03, 1.36981D-03, 1.35916D-03, 1.35086D-03,
     $     1.34591D-03, 1.33808D-03, 1.33553D-03, 1.33326D-03,
     $     1.32747D-03, 1.32366D-03, 1.31541D-03, 1.31332D-03,
     $     1.30941D-03, 1.30428D-03, 1.30249D-03, 1.29312D-03 /

!     Branching ratio of 7014
      data( br(5,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 1.13857D-12, 3.60883D-11,
     $     2.67946D-02, 4.33020D-02, 1.04705D-01, 1.11446D-01,
     $     1.12847D-01, 1.12582D-01, 1.29268D-01, 1.35463D-01,
     $     1.44431D-01, 1.46238D-01, 1.56771D-01, 1.60410D-01,
     $     1.63250D-01, 1.65344D-01, 1.70822D-01, 1.78649D-01,
     $     1.79833D-01, 1.82351D-01, 1.84364D-01, 1.86151D-01,
     $     1.86626D-01, 1.87556D-01, 1.87983D-01, 1.89372D-01,
     $     1.90435D-01, 1.91169D-01, 1.92577D-01, 1.93256D-01,
     $     1.95315D-01, 1.96748D-01, 1.97238D-01, 1.97586D-01,
     $     1.97612D-01, 1.97830D-01, 1.97836D-01, 1.97547D-01,
     $     1.97211D-01, 1.96970D-01, 1.96866D-01, 1.96680D-01,
     $     1.96519D-01, 1.96427D-01, 1.95906D-01, 1.95476D-01,
     $     1.94915D-01, 1.94445D-01, 1.94362D-01, 1.94074D-01,
     $     1.93672D-01, 1.93643D-01, 1.93125D-01, 1.92867D-01,
     $     1.92326D-01, 1.91957D-01, 1.91525D-01, 1.90118D-01,
     $     1.89960D-01, 1.89764D-01, 1.89141D-01, 1.88298D-01,
     $     1.87719D-01, 1.87207D-01, 1.86269D-01, 1.84800D-01,
     $     1.84339D-01, 1.83499D-01, 1.83366D-01, 1.82829D-01,
     $     1.82442D-01, 1.81984D-01, 1.80946D-01, 1.80185D-01,
     $     1.79405D-01, 1.79124D-01, 1.79033D-01, 1.78112D-01,
     $     1.76865D-01, 1.76096D-01, 1.75323D-01, 1.74566D-01,
     $     1.74115D-01, 1.73415D-01, 1.73047D-01, 1.71410D-01,
     $     1.70996D-01, 1.70327D-01, 1.68860D-01, 1.68189D-01,
     $     1.67839D-01, 1.66889D-01, 1.66410D-01, 1.65248D-01,
     $     1.64302D-01, 1.64182D-01, 1.63325D-01, 1.62607D-01,
     $     1.61773D-01, 1.60900D-01, 1.59612D-01, 1.58725D-01,
     $     1.57785D-01, 1.57205D-01, 1.56718D-01, 1.55908D-01,
     $     1.54738D-01, 1.54217D-01, 1.53029D-01, 1.52015D-01,
     $     1.51536D-01, 1.50738D-01, 1.50148D-01, 1.49917D-01,
     $     1.49536D-01, 1.49098D-01, 1.48176D-01, 1.47407D-01,
     $     1.47315D-01, 1.46753D-01, 1.45823D-01, 1.44756D-01,
     $     1.44205D-01, 1.43518D-01, 1.42603D-01, 1.41758D-01,
     $     1.41531D-01, 1.40669D-01, 1.40476D-01, 1.39451D-01,
     $     1.39156D-01, 1.38251D-01, 1.37381D-01, 1.37024D-01,
     $     1.36212D-01, 1.35928D-01, 1.35348D-01, 1.34848D-01,
     $     1.33969D-01, 1.33255D-01, 1.32697D-01, 1.31600D-01,
     $     1.31531D-01, 1.30070D-01, 1.29876D-01, 1.29432D-01,
     $     1.28705D-01, 1.27697D-01, 1.27142D-01, 1.26663D-01,
     $     1.26293D-01, 1.26078D-01, 1.25496D-01, 1.24954D-01,
     $     1.24256D-01, 1.23814D-01, 1.23443D-01, 1.22103D-01,
     $     1.21787D-01, 1.21436D-01, 1.20963D-01, 1.20468D-01,
     $     1.19567D-01, 1.19062D-01, 1.18751D-01, 1.17821D-01,
     $     1.17605D-01, 1.17061D-01, 1.16334D-01, 1.15589D-01,
     $     1.15134D-01, 1.14730D-01, 1.13691D-01, 1.13292D-01,
     $     1.13010D-01, 1.12492D-01, 1.11374D-01, 1.10811D-01,
     $     1.10471D-01, 1.10214D-01, 1.09297D-01, 1.09122D-01,
     $     1.08635D-01, 1.07954D-01, 1.07862D-01, 1.07233D-01,
     $     1.06476D-01, 1.06007D-01, 1.05418D-01, 1.04758D-01,
     $     1.04284D-01, 1.04004D-01, 1.03914D-01, 1.03479D-01,
     $     1.03032D-01, 1.02473D-01, 1.02205D-01, 1.01750D-01,
     $     1.01496D-01, 1.00900D-01, 1.00536D-01, 9.98017D-02,
     $     9.97232D-02, 9.92262D-02, 9.85887D-02, 9.80733D-02,
     $     9.73564D-02, 9.70601D-02, 9.66466D-02, 9.61607D-02,
     $     9.57203D-02, 9.55597D-02, 9.49057D-02, 9.48311D-02,
     $     9.41077D-02, 9.39069D-02, 9.33174D-02, 9.29031D-02,
     $     9.25111D-02, 9.22491D-02, 9.16976D-02, 9.16043D-02,
     $     9.03692D-02, 8.98236D-02, 8.96990D-02, 8.95639D-02,
     $     8.87105D-02, 8.82643D-02, 8.79552D-02, 8.76112D-02,
     $     8.71762D-02, 8.71155D-02, 8.65971D-02, 8.64052D-02,
     $     8.60052D-02, 8.58164D-02, 8.50344D-02, 8.47954D-02,
     $     8.40463D-02, 8.39876D-02, 8.34810D-02, 8.30743D-02,
     $     8.28114D-02, 8.21931D-02, 8.20776D-02, 8.16932D-02,
     $     8.13974D-02, 8.12261D-02, 8.06273D-02, 8.03441D-02,
     $     7.99922D-02, 7.97463D-02, 7.94738D-02, 7.93989D-02,
     $     7.88371D-02, 7.86594D-02, 7.83975D-02, 7.80953D-02,
     $     7.78908D-02, 7.71671D-02, 7.66984D-02, 7.66429D-02,
     $     7.61275D-02, 7.57911D-02, 7.52831D-02, 7.51839D-02,
     $     7.48329D-02, 7.46258D-02, 7.44558D-02, 7.42060D-02,
     $     7.37712D-02, 7.36122D-02, 7.31108D-02, 7.29082D-02,
     $     7.25395D-02, 7.19860D-02, 7.17941D-02, 7.16121D-02,
     $     7.11082D-02, 7.09604D-02, 7.07183D-02, 7.06238D-02,
     $     7.02230D-02, 7.01291D-02, 6.94933D-02, 6.90060D-02,
     $     6.87128D-02, 6.82556D-02, 6.81073D-02, 6.79760D-02,
     $     6.76417D-02, 6.74214D-02, 6.69427D-02, 6.68217D-02,
     $     6.65966D-02, 6.63000D-02, 6.61958D-02, 6.56587D-02 /

!     Branching ratio of 6014
      data( br(6,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     7.67074D-03, 1.10359D-02, 6.50503D-02, 8.35646D-02,
     $     1.05660D-01, 1.10351D-01, 1.12585D-01, 1.11914D-01,
     $     1.09199D-01, 1.07229D-01, 1.02627D-01, 1.01185D-01,
     $     9.98835D-02, 9.89468D-02, 9.71809D-02, 9.46335D-02,
     $     9.41521D-02, 9.26188D-02, 9.00756D-02, 8.87439D-02,
     $     8.81348D-02, 8.63554D-02, 8.48144D-02, 8.26230D-02,
     $     8.21844D-02, 8.20758D-02, 8.18775D-02, 8.17437D-02,
     $     8.10978D-02, 8.04283D-02, 8.00555D-02, 7.96245D-02,
     $     7.95022D-02, 7.92012D-02, 7.87771D-02, 7.83077D-02,
     $     7.79503D-02, 7.74273D-02, 7.70104D-02, 7.67675D-02,
     $     7.64126D-02, 7.63261D-02, 7.58162D-02, 7.54896D-02,
     $     7.51551D-02, 7.49503D-02, 7.49210D-02, 7.47214D-02,
     $     7.44427D-02, 7.43654D-02, 7.41183D-02, 7.39996D-02,
     $     7.37231D-02, 7.35760D-02, 7.34121D-02, 7.28045D-02,
     $     7.27466D-02, 7.26757D-02, 7.24593D-02, 7.21248D-02,
     $     7.18763D-02, 7.16663D-02, 7.12840D-02, 7.07477D-02,
     $     7.05521D-02, 7.02202D-02, 7.01702D-02, 6.98927D-02,
     $     6.97104D-02, 6.95109D-02, 6.90739D-02, 6.86938D-02,
     $     6.83234D-02, 6.81977D-02, 6.81572D-02, 6.77144D-02,
     $     6.70809D-02, 6.67109D-02, 6.63227D-02, 6.59731D-02,
     $     6.57399D-02, 6.54164D-02, 6.52513D-02, 6.44950D-02,
     $     6.43139D-02, 6.40258D-02, 6.34075D-02, 6.31299D-02,
     $     6.29861D-02, 6.25504D-02, 6.23452D-02, 6.18295D-02,
     $     6.14125D-02, 6.13617D-02, 6.09851D-02, 6.06757D-02,
     $     6.02907D-02, 5.99073D-02, 5.93654D-02, 5.89997D-02,
     $     5.86166D-02, 5.83820D-02, 5.81721D-02, 5.78282D-02,
     $     5.73455D-02, 5.71343D-02, 5.66244D-02, 5.62076D-02,
     $     5.60146D-02, 5.56816D-02, 5.54446D-02, 5.53523D-02,
     $     5.52008D-02, 5.50270D-02, 5.46422D-02, 5.43300D-02,
     $     5.42932D-02, 5.40577D-02, 5.36671D-02, 5.32270D-02,
     $     5.29974D-02, 5.27195D-02, 5.23554D-02, 5.20133D-02,
     $     5.19227D-02, 5.15823D-02, 5.15066D-02, 5.10954D-02,
     $     5.09794D-02, 5.06258D-02, 5.02835D-02, 5.01391D-02,
     $     4.98197D-02, 4.97098D-02, 4.94879D-02, 4.92983D-02,
     $     4.89637D-02, 4.86865D-02, 4.84753D-02, 4.80631D-02,
     $     4.80371D-02, 4.74848D-02, 4.74095D-02, 4.72404D-02,
     $     4.69639D-02, 4.65885D-02, 4.63840D-02, 4.62036D-02,
     $     4.60665D-02, 4.59829D-02, 4.57630D-02, 4.55629D-02,
     $     4.53071D-02, 4.51468D-02, 4.50123D-02, 4.45303D-02,
     $     4.44173D-02, 4.42916D-02, 4.41225D-02, 4.39420D-02,
     $     4.36118D-02, 4.34253D-02, 4.33127D-02, 4.29762D-02,
     $     4.28975D-02, 4.27002D-02, 4.24381D-02, 4.21668D-02,
     $     4.20010D-02, 4.18542D-02, 4.14798D-02, 4.13364D-02,
     $     4.12354D-02, 4.10480D-02, 4.06410D-02, 4.04374D-02,
     $     4.03136D-02, 4.02206D-02, 3.98921D-02, 3.98287D-02,
     $     3.96536D-02, 3.94089D-02, 3.93757D-02, 3.91510D-02,
     $     3.88803D-02, 3.87129D-02, 3.85031D-02, 3.82655D-02,
     $     3.80955D-02, 3.79955D-02, 3.79636D-02, 3.78090D-02,
     $     3.76481D-02, 3.74475D-02, 3.73505D-02, 3.71876D-02,
     $     3.70971D-02, 3.68853D-02, 3.67550D-02, 3.64942D-02,
     $     3.64663D-02, 3.62905D-02, 3.60614D-02, 3.58786D-02,
     $     3.56228D-02, 3.55167D-02, 3.53697D-02, 3.51977D-02,
     $     3.50421D-02, 3.49854D-02, 3.47529D-02, 3.47265D-02,
     $     3.44698D-02, 3.43987D-02, 3.41906D-02, 3.40435D-02,
     $     3.39030D-02, 3.38100D-02, 3.36149D-02, 3.35819D-02,
     $     3.31445D-02, 3.29513D-02, 3.29069D-02, 3.28584D-02,
     $     3.25554D-02, 3.23978D-02, 3.22872D-02, 3.21654D-02,
     $     3.20101D-02, 3.19884D-02, 3.18024D-02, 3.17337D-02,
     $     3.15914D-02, 3.15244D-02, 3.12464D-02, 3.11616D-02,
     $     3.08953D-02, 3.08745D-02, 3.06950D-02, 3.05510D-02,
     $     3.04570D-02, 3.02365D-02, 3.01953D-02, 3.00579D-02,
     $     2.99526D-02, 2.98917D-02, 2.96779D-02, 2.95770D-02,
     $     2.94515D-02, 2.93634D-02, 2.92660D-02, 2.92393D-02,
     $     2.90389D-02, 2.89755D-02, 2.88816D-02, 2.87726D-02,
     $     2.86993D-02, 2.84395D-02, 2.82718D-02, 2.82517D-02,
     $     2.80667D-02, 2.79463D-02, 2.77639D-02, 2.77284D-02,
     $     2.76026D-02, 2.75281D-02, 2.74672D-02, 2.73778D-02,
     $     2.72224D-02, 2.71655D-02, 2.69864D-02, 2.69135D-02,
     $     2.67809D-02, 2.65812D-02, 2.65121D-02, 2.64467D-02,
     $     2.62652D-02, 2.62121D-02, 2.61248D-02, 2.60907D-02,
     $     2.59467D-02, 2.59126D-02, 2.56839D-02, 2.55080D-02,
     $     2.54023D-02, 2.52370D-02, 2.51835D-02, 2.51362D-02,
     $     2.50160D-02, 2.49362D-02, 2.47635D-02, 2.47198D-02,
     $     2.46385D-02, 2.45314D-02, 2.44937D-02, 2.43003D-02 /

!     Branching ratio of 7013
      data( br(7,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 1.67429D-03, 1.67590D-02,
     $     3.54248D-02, 3.60665D-02, 3.24412D-02, 3.11145D-02,
     $     2.61349D-02, 2.41084D-02, 2.06041D-02, 1.94682D-02,
     $     1.85153D-02, 1.78750D-02, 1.61871D-02, 1.37942D-02,
     $     1.35799D-02, 1.28975D-02, 1.19172D-02, 1.14823D-02,
     $     1.12675D-02, 1.06177D-02, 1.00348D-02, 9.13035D-03,
     $     8.76412D-03, 8.64275D-03, 8.47948D-03, 8.38348D-03,
     $     8.42520D-03, 8.32051D-03, 8.22557D-03, 8.09573D-03,
     $     8.05726D-03, 8.04602D-03, 8.00770D-03, 7.94902D-03,
     $     7.90620D-03, 7.93190D-03, 7.98675D-03, 7.99581D-03,
     $     8.03853D-03, 8.03955D-03, 8.06562D-03, 8.07902D-03,
     $     8.07236D-03, 8.07115D-03, 8.07113D-03, 8.10787D-03,
     $     8.18269D-03, 8.21740D-03, 8.28679D-03, 8.32092D-03,
     $     8.40685D-03, 8.44490D-03, 8.48455D-03, 8.64803D-03,
     $     8.66149D-03, 8.67795D-03, 8.72742D-03, 8.80481D-03,
     $     8.87902D-03, 8.94592D-03, 9.07878D-03, 9.28406D-03,
     $     9.36988D-03, 9.51830D-03, 9.54069D-03, 9.67360D-03,
     $     9.75693D-03, 9.84274D-03, 1.00166D-02, 1.01654D-02,
     $     1.02990D-02, 1.03408D-02, 1.03540D-02, 1.05001D-02,
     $     1.07129D-02, 1.08297D-02, 1.09551D-02, 1.10644D-02,
     $     1.11501D-02, 1.12586D-02, 1.13135D-02, 1.15812D-02,
     $     1.16443D-02, 1.17458D-02, 1.19696D-02, 1.20716D-02,
     $     1.21245D-02, 1.22984D-02, 1.23755D-02, 1.25690D-02,
     $     1.27217D-02, 1.27394D-02, 1.28715D-02, 1.29759D-02,
     $     1.31125D-02, 1.32540D-02, 1.34302D-02, 1.35489D-02,
     $     1.36666D-02, 1.37343D-02, 1.37975D-02, 1.38975D-02,
     $     1.40295D-02, 1.40857D-02, 1.42483D-02, 1.43795D-02,
     $     1.44395D-02, 1.45522D-02, 1.46320D-02, 1.46636D-02,
     $     1.47158D-02, 1.47764D-02, 1.49207D-02, 1.50362D-02,
     $     1.50497D-02, 1.51410D-02, 1.52908D-02, 1.54542D-02,
     $     1.55387D-02, 1.56360D-02, 1.57584D-02, 1.58724D-02,
     $     1.59012D-02, 1.60053D-02, 1.60277D-02, 1.61498D-02,
     $     1.61821D-02, 1.62759D-02, 1.63763D-02, 1.64172D-02,
     $     1.64980D-02, 1.65233D-02, 1.65714D-02, 1.66094D-02,
     $     1.66742D-02, 1.67391D-02, 1.67946D-02, 1.68734D-02,
     $     1.68779D-02, 1.69683D-02, 1.69810D-02, 1.70067D-02,
     $     1.70463D-02, 1.70885D-02, 1.71078D-02, 1.71282D-02,
     $     1.71464D-02, 1.71633D-02, 1.71933D-02, 1.72236D-02,
     $     1.72525D-02, 1.72669D-02, 1.72776D-02, 1.73082D-02,
     $     1.73141D-02, 1.73202D-02, 1.73308D-02, 1.73517D-02,
     $     1.73787D-02, 1.73939D-02, 1.74003D-02, 1.74165D-02,
     $     1.74208D-02, 1.74382D-02, 1.74632D-02, 1.74846D-02,
     $     1.74961D-02, 1.75044D-02, 1.75199D-02, 1.75247D-02,
     $     1.75278D-02, 1.75354D-02, 1.75539D-02, 1.75615D-02,
     $     1.75668D-02, 1.75702D-02, 1.75778D-02, 1.75801D-02,
     $     1.75843D-02, 1.75888D-02, 1.75894D-02, 1.76007D-02,
     $     1.76052D-02, 1.76062D-02, 1.76045D-02, 1.76052D-02,
     $     1.76037D-02, 1.76016D-02, 1.76008D-02, 1.75963D-02,
     $     1.75938D-02, 1.75896D-02, 1.75883D-02, 1.75836D-02,
     $     1.75801D-02, 1.75711D-02, 1.75660D-02, 1.75650D-02,
     $     1.75659D-02, 1.75624D-02, 1.75664D-02, 1.75604D-02,
     $     1.75591D-02, 1.75573D-02, 1.75518D-02, 1.75428D-02,
     $     1.75332D-02, 1.75294D-02, 1.75158D-02, 1.75140D-02,
     $     1.74971D-02, 1.74919D-02, 1.74758D-02, 1.74714D-02,
     $     1.74662D-02, 1.74606D-02, 1.74466D-02, 1.74440D-02,
     $     1.74177D-02, 1.74045D-02, 1.74029D-02, 1.74012D-02,
     $     1.73853D-02, 1.73726D-02, 1.73707D-02, 1.73621D-02,
     $     1.73512D-02, 1.73497D-02, 1.73408D-02, 1.73365D-02,
     $     1.73248D-02, 1.73184D-02, 1.72903D-02, 1.72844D-02,
     $     1.72581D-02, 1.72557D-02, 1.72342D-02, 1.72157D-02,
     $     1.72051D-02, 1.71784D-02, 1.71733D-02, 1.71563D-02,
     $     1.71421D-02, 1.71338D-02, 1.71056D-02, 1.70915D-02,
     $     1.70786D-02, 1.70679D-02, 1.70549D-02, 1.70512D-02,
     $     1.70252D-02, 1.70175D-02, 1.70058D-02, 1.69968D-02,
     $     1.69878D-02, 1.69568D-02, 1.69344D-02, 1.69319D-02,
     $     1.69110D-02, 1.68939D-02, 1.68678D-02, 1.68624D-02,
     $     1.68429D-02, 1.68316D-02, 1.68218D-02, 1.68070D-02,
     $     1.67805D-02, 1.67706D-02, 1.67390D-02, 1.67268D-02,
     $     1.67044D-02, 1.66750D-02, 1.66643D-02, 1.66535D-02,
     $     1.66276D-02, 1.66210D-02, 1.66089D-02, 1.66037D-02,
     $     1.65798D-02, 1.65746D-02, 1.65404D-02, 1.65118D-02,
     $     1.64962D-02, 1.64693D-02, 1.64599D-02, 1.64514D-02,
     $     1.64289D-02, 1.64148D-02, 1.63856D-02, 1.63780D-02,
     $     1.63635D-02, 1.63453D-02, 1.63393D-02, 1.63032D-02 /

!     Branching ratio of 6013
      data( br(8,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     1.70555D-02, 2.61023D-02, 5.34139D-02, 5.68246D-02,
     $     7.14238D-02, 8.29703D-02, 9.78742D-02, 9.87911D-02,
     $     9.77181D-02, 9.51186D-02, 8.95324D-02, 8.94624D-02,
     $     9.04333D-02, 9.16302D-02, 9.60993D-02, 9.95064D-02,
     $     9.95031D-02, 9.98882D-02, 1.01938D-01, 1.03221D-01,
     $     1.04048D-01, 1.06801D-01, 1.09311D-01, 1.12601D-01,
     $     1.13490D-01, 1.13754D-01, 1.14043D-01, 1.14251D-01,
     $     1.15687D-01, 1.18342D-01, 1.20813D-01, 1.24768D-01,
     $     1.25979D-01, 1.28175D-01, 1.31117D-01, 1.34204D-01,
     $     1.36460D-01, 1.39293D-01, 1.41422D-01, 1.42820D-01,
     $     1.44747D-01, 1.45286D-01, 1.48483D-01, 1.50728D-01,
     $     1.53468D-01, 1.55454D-01, 1.55768D-01, 1.57447D-01,
     $     1.59648D-01, 1.60048D-01, 1.61784D-01, 1.62532D-01,
     $     1.64026D-01, 1.64835D-01, 1.65695D-01, 1.68469D-01,
     $     1.68743D-01, 1.69081D-01, 1.70131D-01, 1.71935D-01,
     $     1.73156D-01, 1.74211D-01, 1.76129D-01, 1.78978D-01,
     $     1.79949D-01, 1.81611D-01, 1.81863D-01, 1.83055D-01,
     $     1.83833D-01, 1.84681D-01, 1.86457D-01, 1.87787D-01,
     $     1.89023D-01, 1.89431D-01, 1.89562D-01, 1.90901D-01,
     $     1.92641D-01, 1.93627D-01, 1.94636D-01, 1.95564D-01,
     $     1.96166D-01, 1.97034D-01, 1.97481D-01, 1.99736D-01,
     $     2.00254D-01, 2.01080D-01, 2.02893D-01, 2.03717D-01,
     $     2.04142D-01, 2.05366D-01, 2.05945D-01, 2.07334D-01,
     $     2.08405D-01, 2.08535D-01, 2.09468D-01, 2.10218D-01,
     $     2.11100D-01, 2.12076D-01, 2.13328D-01, 2.14204D-01,
     $     2.15067D-01, 2.15557D-01, 2.15971D-01, 2.16629D-01,
     $     2.17517D-01, 2.17895D-01, 2.18855D-01, 2.19615D-01,
     $     2.19942D-01, 2.20476D-01, 2.20846D-01, 2.20986D-01,
     $     2.21213D-01, 2.21468D-01, 2.22015D-01, 2.22444D-01,
     $     2.22493D-01, 2.22799D-01, 2.23286D-01, 2.23803D-01,
     $     2.24061D-01, 2.24362D-01, 2.24738D-01, 2.25077D-01,
     $     2.25164D-01, 2.25480D-01, 2.25548D-01, 2.25910D-01,
     $     2.26010D-01, 2.26304D-01, 2.26661D-01, 2.26792D-01,
     $     2.27064D-01, 2.27153D-01, 2.27330D-01, 2.27477D-01,
     $     2.27728D-01, 2.27993D-01, 2.28241D-01, 2.28599D-01,
     $     2.28619D-01, 2.29018D-01, 2.29069D-01, 2.29177D-01,
     $     2.29338D-01, 2.29524D-01, 2.29610D-01, 2.29682D-01,
     $     2.29752D-01, 2.29804D-01, 2.29901D-01, 2.30017D-01,
     $     2.30131D-01, 2.30188D-01, 2.30231D-01, 2.30354D-01,
     $     2.30377D-01, 2.30400D-01, 2.30439D-01, 2.30502D-01,
     $     2.30554D-01, 2.30571D-01, 2.30575D-01, 2.30574D-01,
     $     2.30571D-01, 2.30590D-01, 2.30622D-01, 2.30619D-01,
     $     2.30602D-01, 2.30580D-01, 2.30500D-01, 2.30462D-01,
     $     2.30434D-01, 2.30377D-01, 2.30240D-01, 2.30167D-01,
     $     2.30121D-01, 2.30086D-01, 2.29955D-01, 2.29929D-01,
     $     2.29854D-01, 2.29741D-01, 2.29725D-01, 2.29647D-01,
     $     2.29513D-01, 2.29418D-01, 2.29288D-01, 2.29130D-01,
     $     2.29007D-01, 2.28932D-01, 2.28908D-01, 2.28786D-01,
     $     2.28657D-01, 2.28490D-01, 2.28408D-01, 2.28267D-01,
     $     2.28186D-01, 2.27994D-01, 2.27874D-01, 2.27670D-01,
     $     2.27651D-01, 2.27502D-01, 2.27320D-01, 2.27150D-01,
     $     2.26924D-01, 2.26822D-01, 2.26674D-01, 2.26492D-01,
     $     2.26322D-01, 2.26259D-01, 2.25996D-01, 2.25965D-01,
     $     2.25664D-01, 2.25579D-01, 2.25325D-01, 2.25163D-01,
     $     2.24997D-01, 2.24882D-01, 2.24635D-01, 2.24592D-01,
     $     2.24040D-01, 2.23784D-01, 2.23729D-01, 2.23666D-01,
     $     2.23267D-01, 2.23043D-01, 2.22900D-01, 2.22725D-01,
     $     2.22494D-01, 2.22461D-01, 2.22190D-01, 2.22086D-01,
     $     2.21862D-01, 2.21754D-01, 2.21291D-01, 2.21157D-01,
     $     2.20701D-01, 2.20664D-01, 2.20340D-01, 2.20073D-01,
     $     2.19898D-01, 2.19477D-01, 2.19398D-01, 2.19129D-01,
     $     2.18920D-01, 2.18798D-01, 2.18363D-01, 2.18155D-01,
     $     2.17908D-01, 2.17727D-01, 2.17522D-01, 2.17465D-01,
     $     2.17041D-01, 2.16908D-01, 2.16705D-01, 2.16481D-01,
     $     2.16323D-01, 2.15756D-01, 2.15379D-01, 2.15333D-01,
     $     2.14920D-01, 2.14641D-01, 2.14209D-01, 2.14124D-01,
     $     2.13820D-01, 2.13639D-01, 2.13489D-01, 2.13268D-01,
     $     2.12880D-01, 2.12737D-01, 2.12281D-01, 2.12095D-01,
     $     2.11753D-01, 2.11242D-01, 2.11063D-01, 2.10892D-01,
     $     2.10424D-01, 2.10290D-01, 2.10063D-01, 2.09973D-01,
     $     2.09588D-01, 2.09496D-01, 2.08885D-01, 2.08401D-01,
     $     2.08114D-01, 2.07653D-01, 2.07502D-01, 2.07367D-01,
     $     2.07023D-01, 2.06793D-01, 2.06298D-01, 2.06172D-01,
     $     2.05935D-01, 2.05623D-01, 2.05514D-01, 2.04940D-01 /

!     Branching ratio of 6012
      data( br(9,k), k=1, 300) /
     $     4.01159D-22, 3.44852D-07, 2.70340D-06, 3.42861D-06,
     $     1.43476D-03, 4.26454D-03, 4.40061D-02, 5.39967D-02,
     $     4.84998D-02, 4.38114D-02, 3.92998D-02, 4.05365D-02,
     $     5.20658D-02, 5.78548D-02, 6.57393D-02, 6.75132D-02,
     $     6.85729D-02, 6.91892D-02, 7.08940D-02, 7.46894D-02,
     $     7.49924D-02, 7.62810D-02, 7.86287D-02, 7.97288D-02,
     $     8.04728D-02, 8.34078D-02, 8.69593D-02, 9.43761D-02,
     $     9.78855D-02, 9.94769D-02, 1.02304D-01, 1.04001D-01,
     $     1.06223D-01, 1.09355D-01, 1.11243D-01, 1.13469D-01,
     $     1.14040D-01, 1.14917D-01, 1.15955D-01, 1.16896D-01,
     $     1.17509D-01, 1.18152D-01, 1.18589D-01, 1.18905D-01,
     $     1.19304D-01, 1.19431D-01, 1.20253D-01, 1.20959D-01,
     $     1.21671D-01, 1.22308D-01, 1.22417D-01, 1.22871D-01,
     $     1.23510D-01, 1.23636D-01, 1.24210D-01, 1.24457D-01,
     $     1.24959D-01, 1.25245D-01, 1.25572D-01, 1.26753D-01,
     $     1.26885D-01, 1.27051D-01, 1.27590D-01, 1.28264D-01,
     $     1.28828D-01, 1.29336D-01, 1.30241D-01, 1.31377D-01,
     $     1.31697D-01, 1.32214D-01, 1.32291D-01, 1.32646D-01,
     $     1.32878D-01, 1.33137D-01, 1.33716D-01, 1.34207D-01,
     $     1.34725D-01, 1.34910D-01, 1.34971D-01, 1.35634D-01,
     $     1.36612D-01, 1.37215D-01, 1.37850D-01, 1.38449D-01,
     $     1.38833D-01, 1.39388D-01, 1.39674D-01, 1.40841D-01,
     $     1.41136D-01, 1.41608D-01, 1.42636D-01, 1.43099D-01,
     $     1.43338D-01, 1.44028D-01, 1.44358D-01, 1.45166D-01,
     $     1.45813D-01, 1.45893D-01, 1.46487D-01, 1.46977D-01,
     $     1.47572D-01, 1.48118D-01, 1.48929D-01, 1.49443D-01,
     $     1.49984D-01, 1.50321D-01, 1.50623D-01, 1.51116D-01,
     $     1.51803D-01, 1.52098D-01, 1.52759D-01, 1.53280D-01,
     $     1.53522D-01, 1.53944D-01, 1.54243D-01, 1.54360D-01,
     $     1.54551D-01, 1.54773D-01, 1.55278D-01, 1.55695D-01,
     $     1.55745D-01, 1.56067D-01, 1.56619D-01, 1.57259D-01,
     $     1.57604D-01, 1.58024D-01, 1.58582D-01, 1.59121D-01,
     $     1.59266D-01, 1.59809D-01, 1.59930D-01, 1.60596D-01,
     $     1.60786D-01, 1.61367D-01, 1.61907D-01, 1.62144D-01,
     $     1.62675D-01, 1.62859D-01, 1.63232D-01, 1.63549D-01,
     $     1.64108D-01, 1.64558D-01, 1.64878D-01, 1.65547D-01,
     $     1.65590D-01, 1.66540D-01, 1.66674D-01, 1.66977D-01,
     $     1.67476D-01, 1.68157D-01, 1.68532D-01, 1.68867D-01,
     $     1.69114D-01, 1.69266D-01, 1.69674D-01, 1.70037D-01,
     $     1.70516D-01, 1.70822D-01, 1.71082D-01, 1.72035D-01,
     $     1.72264D-01, 1.72522D-01, 1.72868D-01, 1.73237D-01,
     $     1.73938D-01, 1.74341D-01, 1.74583D-01, 1.75313D-01,
     $     1.75487D-01, 1.75910D-01, 1.76472D-01, 1.77072D-01,
     $     1.77447D-01, 1.77781D-01, 1.78640D-01, 1.78973D-01,
     $     1.79209D-01, 1.79651D-01, 1.80619D-01, 1.81105D-01,
     $     1.81401D-01, 1.81623D-01, 1.82401D-01, 1.82552D-01,
     $     1.82966D-01, 1.83537D-01, 1.83614D-01, 1.84122D-01,
     $     1.84734D-01, 1.85110D-01, 1.85574D-01, 1.86101D-01,
     $     1.86475D-01, 1.86694D-01, 1.86763D-01, 1.87098D-01,
     $     1.87448D-01, 1.87882D-01, 1.88092D-01, 1.88441D-01,
     $     1.88634D-01, 1.89081D-01, 1.89355D-01, 1.89891D-01,
     $     1.89947D-01, 1.90305D-01, 1.90781D-01, 1.91156D-01,
     $     1.91683D-01, 1.91904D-01, 1.92208D-01, 1.92562D-01,
     $     1.92882D-01, 1.92998D-01, 1.93482D-01, 1.93536D-01,
     $     1.94070D-01, 1.94216D-01, 1.94644D-01, 1.94946D-01,
     $     1.95241D-01, 1.95434D-01, 1.95837D-01, 1.95905D-01,
     $     1.96809D-01, 1.97210D-01, 1.97303D-01, 1.97407D-01,
     $     1.98040D-01, 1.98366D-01, 1.98599D-01, 1.98851D-01,
     $     1.99180D-01, 1.99227D-01, 1.99622D-01, 1.99767D-01,
     $     2.00062D-01, 2.00199D-01, 2.00769D-01, 2.00942D-01,
     $     2.01480D-01, 2.01522D-01, 2.01876D-01, 2.02158D-01,
     $     2.02347D-01, 2.02783D-01, 2.02865D-01, 2.03135D-01,
     $     2.03337D-01, 2.03453D-01, 2.03865D-01, 2.04057D-01,
     $     2.04297D-01, 2.04467D-01, 2.04651D-01, 2.04702D-01,
     $     2.05075D-01, 2.05192D-01, 2.05369D-01, 2.05580D-01,
     $     2.05717D-01, 2.06206D-01, 2.06515D-01, 2.06553D-01,
     $     2.06897D-01, 2.07116D-01, 2.07449D-01, 2.07513D-01,
     $     2.07740D-01, 2.07875D-01, 2.07985D-01, 2.08144D-01,
     $     2.08418D-01, 2.08518D-01, 2.08830D-01, 2.08959D-01,
     $     2.09194D-01, 2.09550D-01, 2.09670D-01, 2.09783D-01,
     $     2.10098D-01, 2.10190D-01, 2.10341D-01, 2.10399D-01,
     $     2.10639D-01, 2.10697D-01, 2.11074D-01, 2.11359D-01,
     $     2.11529D-01, 2.11791D-01, 2.11873D-01, 2.11944D-01,
     $     2.12122D-01, 2.12243D-01, 2.12501D-01, 2.12566D-01,
     $     2.12685D-01, 2.12842D-01, 2.12898D-01, 2.13169D-01 /

!     Branching ratio of 5012
      data( br(10,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 3.85427D-71, 1.28965D-26, 1.87231D-25,
     $     3.32809D-20, 1.30104D-16, 1.09866D-15, 1.74049D-15,
     $     4.36343D-15, 3.56672D-14, 2.09270D-13, 1.24072D-12,
     $     1.45466D-12, 5.62442D-12, 3.96685D-11, 2.33364D-10,
     $     4.56555D-10, 1.13630D-09, 3.20588D-09, 7.87509D-09,
     $     1.00789D-08, 2.37403D-08, 2.84843D-08, 6.95187D-08,
     $     8.84182D-08, 1.73818D-07, 3.22732D-07, 4.09152D-07,
     $     7.16521D-07, 8.69647D-07, 1.30960D-06, 1.81141D-06,
     $     3.19024D-06, 4.90030D-06, 6.70033D-06, 1.23910D-05,
     $     1.28473D-05, 2.68890D-05, 2.96080D-05, 3.66650D-05,
     $     5.07752D-05, 7.69870D-05, 9.57433D-05, 1.14438D-04,
     $     1.30219D-04, 1.39681D-04, 1.68455D-04, 1.98176D-04,
     $     2.42312D-04, 2.73242D-04, 3.00952D-04, 4.12277D-04,
     $     4.41232D-04, 4.74498D-04, 5.20229D-04, 5.68172D-04,
     $     6.59293D-04, 7.12616D-04, 7.46378D-04, 8.49773D-04,
     $     8.73871D-04, 9.33949D-04, 1.01460D-03, 1.09868D-03,
     $     1.15069D-03, 1.19725D-03, 1.31812D-03, 1.36524D-03,
     $     1.39873D-03, 1.46066D-03, 1.59554D-03, 1.66500D-03,
     $     1.70785D-03, 1.74100D-03, 1.86363D-03, 1.88744D-03,
     $     1.95531D-03, 2.05417D-03, 2.06798D-03, 2.16295D-03,
     $     2.28313D-03, 2.35979D-03, 2.45979D-03, 2.57471D-03,
     $     2.66006D-03, 2.71165D-03, 2.72827D-03, 2.80993D-03,
     $     2.89494D-03, 3.00336D-03, 3.05592D-03, 3.14628D-03,
     $     3.19697D-03, 3.31578D-03, 3.38875D-03, 3.53606D-03,
     $     3.55173D-03, 3.65253D-03, 3.77873D-03, 3.88345D-03,
     $     4.02766D-03, 4.08712D-03, 4.17038D-03, 4.26826D-03,
     $     4.35671D-03, 4.38895D-03, 4.51928D-03, 4.53426D-03,
     $     4.67927D-03, 4.71968D-03, 4.83842D-03, 4.92065D-03,
     $     4.99853D-03, 5.05147D-03, 5.16392D-03, 5.18300D-03,
     $     5.43466D-03, 5.54712D-03, 5.57246D-03, 5.59981D-03,
     $     5.77559D-03, 5.86977D-03, 5.93371D-03, 6.00714D-03,
     $     6.10014D-03, 6.11314D-03, 6.22365D-03, 6.26506D-03,
     $     6.35292D-03, 6.39491D-03, 6.57042D-03, 6.62384D-03,
     $     6.79415D-03, 6.80769D-03, 6.92587D-03, 7.02211D-03,
     $     7.08383D-03, 7.23119D-03, 7.25889D-03, 7.35131D-03,
     $     7.42335D-03, 7.46518D-03, 7.61138D-03, 7.68116D-03,
     $     7.76672D-03, 7.82672D-03, 7.89386D-03, 7.91240D-03,
     $     8.05143D-03, 8.09555D-03, 8.16038D-03, 8.23406D-03,
     $     8.28516D-03, 8.46629D-03, 8.58470D-03, 8.59845D-03,
     $     8.72710D-03, 8.81230D-03, 8.94044D-03, 8.96557D-03,
     $     9.05434D-03, 9.10626D-03, 9.14916D-03, 9.21227D-03,
     $     9.32222D-03, 9.36241D-03, 9.48876D-03, 9.53909D-03,
     $     9.63038D-03, 9.76474D-03, 9.81128D-03, 9.85565D-03,
     $     9.97663D-03, 1.00118D-02, 1.00694D-02, 1.00920D-02,
     $     1.01881D-02, 1.02103D-02, 1.03609D-02, 1.04752D-02,
     $     1.05434D-02, 1.06493D-02, 1.06839D-02, 1.07145D-02,
     $     1.07927D-02, 1.08434D-02, 1.09535D-02, 1.09813D-02,
     $     1.10332D-02, 1.11011D-02, 1.11248D-02, 1.12485D-02 /

!     Branching ratio of 4012
      data( br(11,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00 /

!     Branching ratio of 6011
      data( br(12,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     1.26338D-07, 2.52652D-05, 9.50595D-04, 3.23291D-03,
     $     5.68427D-03, 6.80336D-03, 1.27921D-02, 1.45158D-02,
     $     1.58658D-02, 1.68250D-02, 1.95187D-02, 2.69965D-02,
     $     2.81486D-02, 3.14947D-02, 3.62700D-02, 3.87732D-02,
     $     3.98436D-02, 4.28475D-02, 4.53816D-02, 5.03854D-02,
     $     5.29068D-02, 5.40298D-02, 5.58813D-02, 5.68070D-02,
     $     5.86199D-02, 6.01988D-02, 6.09708D-02, 6.18955D-02,
     $     6.21547D-02, 6.28109D-02, 6.36949D-02, 6.46324D-02,
     $     6.53134D-02, 6.62907D-02, 6.70766D-02, 6.75573D-02,
     $     6.82911D-02, 6.84724D-02, 6.95465D-02, 7.02783D-02,
     $     7.09191D-02, 7.14419D-02, 7.15297D-02, 7.19920D-02,
     $     7.26219D-02, 7.27960D-02, 7.33667D-02, 7.36112D-02,
     $     7.41178D-02, 7.43589D-02, 7.46077D-02, 7.54328D-02,
     $     7.55018D-02, 7.55847D-02, 7.58240D-02, 7.60906D-02,
     $     7.63249D-02, 7.65229D-02, 7.68656D-02, 7.72178D-02,
     $     7.73361D-02, 7.74912D-02, 7.75091D-02, 7.76598D-02,
     $     7.77370D-02, 7.78025D-02, 7.78996D-02, 7.80064D-02,
     $     7.80707D-02, 7.80819D-02, 7.80845D-02, 7.81332D-02,
     $     7.81953D-02, 7.81975D-02, 7.82055D-02, 7.81754D-02,
     $     7.81738D-02, 7.81345D-02, 7.81083D-02, 7.79762D-02,
     $     7.79377D-02, 7.78712D-02, 7.77044D-02, 7.76225D-02,
     $     7.75788D-02, 7.74970D-02, 7.74438D-02, 7.73314D-02,
     $     7.72350D-02, 7.72211D-02, 7.71297D-02, 7.70457D-02,
     $     7.69663D-02, 7.68498D-02, 7.66766D-02, 7.65465D-02,
     $     7.64120D-02, 7.63315D-02, 7.62747D-02, 7.61761D-02,
     $     7.60206D-02, 7.59473D-02, 7.57852D-02, 7.56340D-02,
     $     7.55622D-02, 7.54473D-02, 7.53563D-02, 7.53199D-02,
     $     7.52591D-02, 7.51884D-02, 7.50472D-02, 7.49184D-02,
     $     7.49025D-02, 7.48082D-02, 7.46474D-02, 7.44500D-02,
     $     7.43457D-02, 7.42079D-02, 7.40169D-02, 7.38366D-02,
     $     7.37862D-02, 7.35906D-02, 7.35460D-02, 7.33077D-02,
     $     7.32362D-02, 7.30113D-02, 7.27834D-02, 7.26911D-02,
     $     7.24759D-02, 7.23991D-02, 7.22407D-02, 7.21029D-02,
     $     7.18600D-02, 7.16581D-02, 7.14924D-02, 7.11734D-02,
     $     7.11534D-02, 7.07228D-02, 7.06654D-02, 7.05332D-02,
     $     7.03142D-02, 7.00058D-02, 6.98340D-02, 6.96850D-02,
     $     6.95670D-02, 6.94989D-02, 6.93133D-02, 6.91351D-02,
     $     6.89041D-02, 6.87573D-02, 6.86328D-02, 6.81777D-02,
     $     6.80687D-02, 6.79467D-02, 6.77806D-02, 6.76040D-02,
     $     6.72848D-02, 6.71045D-02, 6.69930D-02, 6.66566D-02,
     $     6.65783D-02, 6.63766D-02, 6.61044D-02, 6.58255D-02,
     $     6.56547D-02, 6.55021D-02, 6.51077D-02, 6.49551D-02,
     $     6.48473D-02, 6.46485D-02, 6.42165D-02, 6.39981D-02,
     $     6.38655D-02, 6.37650D-02, 6.34047D-02, 6.33358D-02,
     $     6.31437D-02, 6.28741D-02, 6.28376D-02, 6.25858D-02,
     $     6.22828D-02, 6.20954D-02, 6.18589D-02, 6.15937D-02,
     $     6.14022D-02, 6.12885D-02, 6.12521D-02, 6.10749D-02,
     $     6.08919D-02, 6.06630D-02, 6.05524D-02, 6.03644D-02,
     $     6.02592D-02, 6.00118D-02, 5.98600D-02, 5.95500D-02,
     $     5.95165D-02, 5.93059D-02, 5.90339D-02, 5.88145D-02,
     $     5.85069D-02, 5.83798D-02, 5.82025D-02, 5.79940D-02,
     $     5.78051D-02, 5.77362D-02, 5.74552D-02, 5.74231D-02,
     $     5.71116D-02, 5.70251D-02, 5.67710D-02, 5.65910D-02,
     $     5.64210D-02, 5.63075D-02, 5.60685D-02, 5.60281D-02,
     $     5.54908D-02, 5.52533D-02, 5.51987D-02, 5.51396D-02,
     $     5.47663D-02, 5.45717D-02, 5.44356D-02, 5.42854D-02,
     $     5.40956D-02, 5.40691D-02, 5.38422D-02, 5.37585D-02,
     $     5.35845D-02, 5.35026D-02, 5.31637D-02, 5.30597D-02,
     $     5.27361D-02, 5.27108D-02, 5.24932D-02, 5.23189D-02,
     $     5.22062D-02, 5.19419D-02, 5.18926D-02, 5.17288D-02,
     $     5.16032D-02, 5.15306D-02, 5.12768D-02, 5.11571D-02,
     $     5.10075D-02, 5.09034D-02, 5.07886D-02, 5.07571D-02,
     $     5.05206D-02, 5.04458D-02, 5.03357D-02, 5.02078D-02,
     $     5.01219D-02, 4.98186D-02, 4.96232D-02, 4.96000D-02,
     $     4.93850D-02, 4.92457D-02, 4.90358D-02, 4.89950D-02,
     $     4.88506D-02, 4.87655D-02, 4.86957D-02, 4.85935D-02,
     $     4.84161D-02, 4.83513D-02, 4.81474D-02, 4.80650D-02,
     $     4.79151D-02, 4.76897D-02, 4.76118D-02, 4.75381D-02,
     $     4.73335D-02, 4.72734D-02, 4.71752D-02, 4.71370D-02,
     $     4.69755D-02, 4.69375D-02, 4.66809D-02, 4.64848D-02,
     $     4.63664D-02, 4.61821D-02, 4.61224D-02, 4.60697D-02,
     $     4.59355D-02, 4.58467D-02, 4.56534D-02, 4.56045D-02,
     $     4.55136D-02, 4.53935D-02, 4.53512D-02, 4.51342D-02 /

!     Branching ratio of 5011
      data( br(13,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     1.24186D-05, 6.21164D-05, 1.91219D-02, 3.53050D-02,
     $     1.22232D-01, 1.59621D-01, 2.13636D-01, 2.22716D-01,
     $     2.66874D-01, 2.86459D-01, 3.11716D-01, 3.18504D-01,
     $     3.24015D-01, 3.27091D-01, 3.33047D-01, 3.40086D-01,
     $     3.40417D-01, 3.42310D-01, 3.45158D-01, 3.45078D-01,
     $     3.45149D-01, 3.45103D-01, 3.44583D-01, 3.41659D-01,
     $     3.39650D-01, 3.38409D-01, 3.36100D-01, 3.34844D-01,
     $     3.31259D-01, 3.27216D-01, 3.24479D-01, 3.20680D-01,
     $     3.19616D-01, 3.17501D-01, 3.14942D-01, 3.12523D-01,
     $     3.10851D-01, 3.08506D-01, 3.06606D-01, 3.05434D-01,
     $     3.03671D-01, 3.03216D-01, 3.00359D-01, 2.98232D-01,
     $     2.95792D-01, 2.93845D-01, 2.93522D-01, 2.91731D-01,
     $     2.89232D-01, 2.88676D-01, 2.86613D-01, 2.85724D-01,
     $     2.83900D-01, 2.82952D-01, 2.81931D-01, 2.78373D-01,
     $     2.78017D-01, 2.77575D-01, 2.76185D-01, 2.74078D-01,
     $     2.72439D-01, 2.70980D-01, 2.68290D-01, 2.64498D-01,
     $     2.63223D-01, 2.61071D-01, 2.60749D-01, 2.59130D-01,
     $     2.58084D-01, 2.56952D-01, 2.54557D-01, 2.52632D-01,
     $     2.50789D-01, 2.50171D-01, 2.49973D-01, 2.47860D-01,
     $     2.44951D-01, 2.43281D-01, 2.41558D-01, 2.40003D-01,
     $     2.38981D-01, 2.37550D-01, 2.36821D-01, 2.33537D-01,
     $     2.32748D-01, 2.31491D-01, 2.28770D-01, 2.27542D-01,
     $     2.26907D-01, 2.25018D-01, 2.24126D-01, 2.21926D-01,
     $     2.20174D-01, 2.19961D-01, 2.18395D-01, 2.17114D-01,
     $     2.15547D-01, 2.13994D-01, 2.11803D-01, 2.10339D-01,
     $     2.08811D-01, 2.07880D-01, 2.07053D-01, 2.05699D-01,
     $     2.03802D-01, 2.02971D-01, 2.00984D-01, 1.99359D-01,
     $     1.98606D-01, 1.97311D-01, 1.96387D-01, 1.96028D-01,
     $     1.95439D-01, 1.94764D-01, 1.93279D-01, 1.92074D-01,
     $     1.91932D-01, 1.91028D-01, 1.89535D-01, 1.87854D-01,
     $     1.86978D-01, 1.85918D-01, 1.84528D-01, 1.83222D-01,
     $     1.82876D-01, 1.81573D-01, 1.81283D-01, 1.79712D-01,
     $     1.79268D-01, 1.77909D-01, 1.76596D-01, 1.76042D-01,
     $     1.74810D-01, 1.74384D-01, 1.73522D-01, 1.72782D-01,
     $     1.71477D-01, 1.70397D-01, 1.69571D-01, 1.67949D-01,
     $     1.67847D-01, 1.65669D-01, 1.65373D-01, 1.64708D-01,
     $     1.63622D-01, 1.62148D-01, 1.61347D-01, 1.60643D-01,
     $     1.60108D-01, 1.59784D-01, 1.58929D-01, 1.58149D-01,
     $     1.57150D-01, 1.56522D-01, 1.55995D-01, 1.54105D-01,
     $     1.53662D-01, 1.53169D-01, 1.52506D-01, 1.51803D-01,
     $     1.50519D-01, 1.49796D-01, 1.49359D-01, 1.48057D-01,
     $     1.47754D-01, 1.46995D-01, 1.45989D-01, 1.44951D-01,
     $     1.44318D-01, 1.43759D-01, 1.42334D-01, 1.41789D-01,
     $     1.41406D-01, 1.40697D-01, 1.39163D-01, 1.38396D-01,
     $     1.37931D-01, 1.37582D-01, 1.36346D-01, 1.36109D-01,
     $     1.35452D-01, 1.34536D-01, 1.34412D-01, 1.33574D-01,
     $     1.32568D-01, 1.31947D-01, 1.31172D-01, 1.30299D-01,
     $     1.29676D-01, 1.29309D-01, 1.29193D-01, 1.28627D-01,
     $     1.28041D-01, 1.27311D-01, 1.26959D-01, 1.26368D-01,
     $     1.26040D-01, 1.25273D-01, 1.24802D-01, 1.23862D-01,
     $     1.23761D-01, 1.23129D-01, 1.22307D-01, 1.21652D-01,
     $     1.20737D-01, 1.20359D-01, 1.19835D-01, 1.19223D-01,
     $     1.18670D-01, 1.18469D-01, 1.17647D-01, 1.17554D-01,
     $     1.16650D-01, 1.16401D-01, 1.15672D-01, 1.15160D-01,
     $     1.14672D-01, 1.14350D-01, 1.13677D-01, 1.13564D-01,
     $     1.12075D-01, 1.11426D-01, 1.11277D-01, 1.11115D-01,
     $     1.10112D-01, 1.09595D-01, 1.09234D-01, 1.08838D-01,
     $     1.08336D-01, 1.08266D-01, 1.07670D-01, 1.07451D-01,
     $     1.06998D-01, 1.06786D-01, 1.05911D-01, 1.05645D-01,
     $     1.04815D-01, 1.04751D-01, 1.04196D-01, 1.03752D-01,
     $     1.03463D-01, 1.02787D-01, 1.02662D-01, 1.02243D-01,
     $     1.01923D-01, 1.01739D-01, 1.01093D-01, 1.00790D-01,
     $     1.00414D-01, 1.00150D-01, 9.98609D-02, 9.97816D-02,
     $     9.91897D-02, 9.90035D-02, 9.87278D-02, 9.84084D-02,
     $     9.81948D-02, 9.74431D-02, 9.69632D-02, 9.69057D-02,
     $     9.63795D-02, 9.60404D-02, 9.55291D-02, 9.54301D-02,
     $     9.50807D-02, 9.48744D-02, 9.47064D-02, 9.44609D-02,
     $     9.40363D-02, 9.38818D-02, 9.33978D-02, 9.32016D-02,
     $     9.28464D-02, 9.23156D-02, 9.21336D-02, 9.19624D-02,
     $     9.14897D-02, 9.13523D-02, 9.11276D-02, 9.10404D-02,
     $     9.06751D-02, 9.05890D-02, 9.00173D-02, 8.95836D-02,
     $     8.93258D-02, 8.89262D-02, 8.87985D-02, 8.86862D-02,
     $     8.84029D-02, 8.82151D-02, 8.78126D-02, 8.77116D-02,
     $     8.75245D-02, 8.72795D-02, 8.71937D-02, 8.67591D-02 /

!     Branching ratio of 4011
      data( br(14,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00 /

!     Branching ratio of 6010
      data( br(15,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 6.63050D-46, 9.25281D-26,
     $     4.72673D-24, 2.78040D-20, 3.42725D-20, 5.20943D-19,
     $     7.39283D-19, 7.58664D-19, 2.04241D-17, 3.36530D-17,
     $     1.52864D-16, 2.46552D-16, 3.28068D-16, 4.56031D-16,
     $     1.94330D-15, 4.03445D-15, 4.92088D-14, 2.23588D-13,
     $     4.59499D-13, 5.72287D-13, 5.12463D-11, 7.16439D-11,
     $     7.96459D-11, 3.58532D-10, 4.36625D-10, 4.62886D-10,
     $     5.00420D-10, 5.43530D-10, 7.26214D-10, 8.56761D-10,
     $     8.72823D-10, 2.95343D-09, 6.40638D-09, 9.69226D-09,
     $     1.10228D-08, 1.23024D-08, 1.37644D-08, 1.64883D-08,
     $     1.70533D-08, 1.90243D-08, 1.94641D-08, 2.69919D-08,
     $     2.85431D-08, 3.33870D-08, 6.54322D-08, 8.77393D-08,
     $     1.21755D-07, 1.30229D-07, 1.45646D-07, 1.57643D-07,
     $     2.06541D-07, 3.95761D-07, 9.32898D-07, 1.37209D-06,
     $     1.39397D-06, 1.84942D-06, 1.93075D-06, 2.08367D-06,
     $     2.33732D-06, 2.59021D-06, 2.71466D-06, 2.96993D-06,
     $     3.78937D-06, 4.64672D-06, 5.88046D-06, 7.17439D-06,
     $     8.37551D-06, 9.00340D-06, 9.48165D-06, 1.10340D-05,
     $     1.13758D-05, 1.17554D-05, 1.24977D-05, 1.38869D-05,
     $     1.56719D-05, 1.66714D-05, 1.71955D-05, 1.87969D-05,
     $     1.92340D-05, 2.18267D-05, 2.52783D-05, 2.81818D-05,
     $     2.97913D-05, 3.11423D-05, 3.43576D-05, 3.55696D-05,
     $     3.64214D-05, 3.81399D-05, 4.22388D-05, 4.43738D-05,
     $     4.57720D-05, 4.68322D-05, 5.05606D-05, 5.13602D-05,
     $     5.35209D-05, 5.66163D-05, 5.70446D-05, 6.10078D-05,
     $     6.50549D-05, 6.75200D-05, 7.05653D-05, 7.43520D-05,
     $     7.70343D-05, 7.85703D-05, 7.90567D-05, 8.14558D-05,
     $     8.42201D-05, 8.77216D-05, 8.95382D-05, 9.24772D-05,
     $     9.40911D-05, 9.79360D-05, 1.00417D-04, 1.07083D-04,
     $     1.07982D-04, 1.12413D-04, 1.19154D-04, 1.23631D-04,
     $     1.31190D-04, 1.34165D-04, 1.38046D-04, 1.42409D-04,
     $     1.46294D-04, 1.47710D-04, 1.53768D-04, 1.54446D-04,
     $     1.61260D-04, 1.63146D-04, 1.68711D-04, 1.73848D-04,
     $     1.78484D-04, 1.81340D-04, 1.87144D-04, 1.88106D-04,
     $     2.02508D-04, 2.08841D-04, 2.10531D-04, 2.12345D-04,
     $     2.23457D-04, 2.28739D-04, 2.33591D-04, 2.37983D-04,
     $     2.43557D-04, 2.44348D-04, 2.51829D-04, 2.54486D-04,
     $     2.59672D-04, 2.62020D-04, 2.71711D-04, 2.75230D-04,
     $     2.85173D-04, 2.85917D-04, 2.92291D-04, 2.97342D-04,
     $     3.00867D-04, 3.09089D-04, 3.10639D-04, 3.15862D-04,
     $     3.19807D-04, 3.22100D-04, 3.30441D-04, 3.34345D-04,
     $     3.40133D-04, 3.43930D-04, 3.47959D-04, 3.49051D-04,
     $     3.57814D-04, 3.60735D-04, 3.65008D-04, 3.70859D-04,
     $     3.74315D-04, 3.86944D-04, 3.94909D-04, 3.95904D-04,
     $     4.05663D-04, 4.11444D-04, 4.20274D-04, 4.21950D-04,
     $     4.27897D-04, 4.31503D-04, 4.34375D-04, 4.38551D-04,
     $     4.45771D-04, 4.48418D-04, 4.56781D-04, 4.60360D-04,
     $     4.66906D-04, 4.77934D-04, 4.81718D-04, 4.85185D-04,
     $     4.95879D-04, 4.99306D-04, 5.04664D-04, 5.06651D-04,
     $     5.14770D-04, 5.16776D-04, 5.30903D-04, 5.41525D-04,
     $     5.48472D-04, 5.58868D-04, 5.62138D-04, 5.64987D-04,
     $     5.72156D-04, 5.77156D-04, 5.88542D-04, 5.91420D-04,
     $     5.96705D-04, 6.04021D-04, 6.06724D-04, 6.19450D-04 /

!     Branching ratio of 5010
      data( br(16,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     4.60981D-30, 6.58832D-27, 2.08200D-24, 3.14470D-23,
     $     7.13854D-22, 4.65129D-19, 6.97558D-18, 6.86355D-17,
     $     1.20956D-16, 3.25047D-16, 1.38056D-15, 6.46551D-15,
     $     1.52935D-14, 3.38303D-14, 5.30185D-14, 6.88279D-14,
     $     1.20920D-13, 1.41258D-13, 3.38376D-12, 4.12050D-11,
     $     1.22238D-09, 9.15323D-09, 1.39129D-08, 1.10524D-07,
     $     7.45581D-07, 9.96623D-07, 2.83551D-06, 4.06350D-06,
     $     7.27972D-06, 9.46631D-06, 1.21427D-05, 2.29756D-05,
     $     2.42034D-05, 2.57734D-05, 3.10101D-05, 3.96384D-05,
     $     4.90659D-05, 6.01620D-05, 9.22613D-05, 1.70779D-04,
     $     2.06917D-04, 2.75154D-04, 2.85992D-04, 3.45751D-04,
     $     3.86374D-04, 4.31404D-04, 5.31600D-04, 6.24211D-04,
     $     7.17877D-04, 7.50598D-04, 7.61255D-04, 8.83891D-04,
     $     1.07691D-03, 1.19510D-03, 1.32749D-03, 1.44860D-03,
     $     1.53887D-03, 1.66160D-03, 1.72527D-03, 2.02805D-03,
     $     2.10317D-03, 2.22507D-03, 2.49426D-03, 2.61923D-03,
     $     2.68483D-03, 2.91271D-03, 3.01440D-03, 3.28668D-03,
     $     3.50755D-03, 3.53347D-03, 3.73750D-03, 3.90315D-03,
     $     4.13449D-03, 4.34443D-03, 4.63844D-03, 4.83145D-03,
     $     5.03812D-03, 5.16775D-03, 5.30421D-03, 5.52772D-03,
     $     5.83251D-03, 5.96526D-03, 6.33032D-03, 6.61483D-03,
     $     6.74707D-03, 7.00073D-03, 7.17342D-03, 7.24092D-03,
     $     7.35203D-03, 7.48156D-03, 7.81544D-03, 8.07827D-03,
     $     8.10901D-03, 8.32765D-03, 8.70384D-03, 9.12801D-03,
     $     9.36163D-03, 9.63356D-03, 9.99056D-03, 1.03534D-02,
     $     1.04483D-02, 1.08073D-02, 1.08879D-02, 1.13591D-02,
     $     1.14900D-02, 1.18942D-02, 1.22956D-02, 1.24825D-02,
     $     1.28862D-02, 1.30241D-02, 1.33034D-02, 1.35447D-02,
     $     1.39927D-02, 1.43864D-02, 1.46717D-02, 1.52649D-02,
     $     1.53039D-02, 1.61631D-02, 1.62902D-02, 1.65699D-02,
     $     1.70326D-02, 1.76536D-02, 1.79937D-02, 1.83107D-02,
     $     1.85445D-02, 1.87024D-02, 1.91019D-02, 1.94499D-02,
     $     1.98996D-02, 2.01841D-02, 2.04244D-02, 2.12964D-02,
     $     2.15031D-02, 2.17342D-02, 2.20452D-02, 2.23883D-02,
     $     2.30424D-02, 2.34209D-02, 2.36429D-02, 2.43068D-02,
     $     2.44666D-02, 2.48550D-02, 2.53645D-02, 2.59136D-02,
     $     2.62538D-02, 2.65531D-02, 2.73110D-02, 2.76021D-02,
     $     2.78077D-02, 2.81996D-02, 2.90728D-02, 2.95118D-02,
     $     2.97848D-02, 2.99887D-02, 3.07004D-02, 3.08436D-02,
     $     3.12339D-02, 3.17841D-02, 3.18591D-02, 3.23564D-02,
     $     3.29631D-02, 3.33420D-02, 3.38121D-02, 3.43649D-02,
     $     3.47560D-02, 3.49832D-02, 3.50553D-02, 3.54048D-02,
     $     3.57796D-02, 3.62494D-02, 3.64803D-02, 3.68607D-02,
     $     3.70700D-02, 3.75611D-02, 3.78698D-02, 3.84696D-02,
     $     3.85331D-02, 3.89365D-02, 3.94859D-02, 3.99165D-02,
     $     4.05278D-02, 4.07856D-02, 4.11395D-02, 4.15507D-02,
     $     4.19234D-02, 4.20593D-02, 4.26313D-02, 4.26952D-02,
     $     4.33275D-02, 4.35015D-02, 4.40101D-02, 4.43709D-02,
     $     4.47283D-02, 4.49602D-02, 4.54417D-02, 4.55227D-02,
     $     4.66033D-02, 4.70816D-02, 4.71917D-02, 4.73161D-02,
     $     4.80629D-02, 4.84460D-02, 4.87203D-02, 4.90158D-02,
     $     4.94036D-02, 4.94581D-02, 4.99249D-02, 5.00961D-02,
     $     5.04450D-02, 5.06076D-02, 5.12873D-02, 5.14910D-02,
     $     5.21379D-02, 5.21882D-02, 5.26202D-02, 5.29653D-02,
     $     5.31990D-02, 5.37405D-02, 5.38420D-02, 5.41820D-02,
     $     5.44386D-02, 5.45868D-02, 5.51133D-02, 5.53597D-02,
     $     5.56634D-02, 5.58827D-02, 5.61217D-02, 5.61870D-02,
     $     5.66737D-02, 5.68267D-02, 5.70595D-02, 5.73309D-02,
     $     5.75106D-02, 5.81489D-02, 5.85546D-02, 5.86056D-02,
     $     5.90557D-02, 5.93449D-02, 5.97882D-02, 5.98736D-02,
     $     6.01768D-02, 6.03589D-02, 6.05056D-02, 6.07192D-02,
     $     6.10882D-02, 6.12228D-02, 6.16452D-02, 6.18214D-02,
     $     6.21422D-02, 6.26255D-02, 6.27914D-02, 6.29464D-02,
     $     6.33780D-02, 6.35027D-02, 6.37102D-02, 6.37910D-02,
     $     6.41282D-02, 6.42103D-02, 6.47438D-02, 6.51611D-02,
     $     6.54080D-02, 6.58027D-02, 6.59287D-02, 6.60392D-02,
     $     6.63186D-02, 6.65095D-02, 6.69154D-02, 6.70183D-02,
     $     6.72104D-02, 6.74628D-02, 6.75510D-02, 6.79998D-02 /

!     Branching ratio of 4010
      data( br(17,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 7.42462D-92, 1.53438D-21, 8.82931D-18,
     $     2.57768D-16, 5.49848D-16, 2.46504D-14, 1.72215D-13,
     $     1.31777D-12, 5.20962D-12, 6.41024D-12, 2.97923D-11,
     $     7.33750D-10, 1.13017D-09, 7.35913D-09, 1.42601D-08,
     $     4.15453D-08, 6.89869D-08, 1.21081D-07, 6.13276D-07,
     $     7.04057D-07, 8.30172D-07, 1.31122D-06, 2.15140D-06,
     $     2.93031D-06, 3.70935D-06, 5.33558D-06, 8.36024D-06,
     $     9.53085D-06, 1.18084D-05, 1.21878D-05, 1.39616D-05,
     $     1.52559D-05, 1.68049D-05, 2.06428D-05, 2.39870D-05,
     $     2.76549D-05, 2.90010D-05, 2.94456D-05, 3.41846D-05,
     $     4.09700D-05, 4.51554D-05, 4.94881D-05, 5.37578D-05,
     $     5.65120D-05, 6.08766D-05, 6.32625D-05, 7.43774D-05,
     $     7.74848D-05, 8.27823D-05, 9.57134D-05, 1.02311D-04,
     $     1.05919D-04, 1.16484D-04, 1.22046D-04, 1.36273D-04,
     $     1.48574D-04, 1.50175D-04, 1.62009D-04, 1.72529D-04,
     $     1.85569D-04, 1.99555D-04, 2.21948D-04, 2.37740D-04,
     $     2.55290D-04, 2.66623D-04, 2.76400D-04, 2.93225D-04,
     $     3.19144D-04, 3.31446D-04, 3.61225D-04, 3.89430D-04,
     $     4.03721D-04, 4.28714D-04, 4.47939D-04, 4.55584D-04,
     $     4.68349D-04, 4.83208D-04, 5.14714D-04, 5.41538D-04,
     $     5.44767D-04, 5.64518D-04, 5.97375D-04, 6.35884D-04,
     $     6.56150D-04, 6.82059D-04, 7.17605D-04, 7.51461D-04,
     $     7.60815D-04, 7.97145D-04, 8.05441D-04, 8.50043D-04,
     $     8.63186D-04, 9.04442D-04, 9.44337D-04, 9.61242D-04,
     $     1.00138D-03, 1.01602D-03, 1.04687D-03, 1.07473D-03,
     $     1.12583D-03, 1.16807D-03, 1.20171D-03, 1.27233D-03,
     $     1.27686D-03, 1.37442D-03, 1.38740D-03, 1.41732D-03,
     $     1.46687D-03, 1.53699D-03, 1.57631D-03, 1.61031D-03,
     $     1.63663D-03, 1.65137D-03, 1.69344D-03, 1.73307D-03,
     $     1.78640D-03, 1.82117D-03, 1.85107D-03, 1.96308D-03,
     $     1.99025D-03, 2.02079D-03, 2.06196D-03, 2.10454D-03,
     $     2.18376D-03, 2.22855D-03, 2.25652D-03, 2.34089D-03,
     $     2.36035D-03, 2.40892D-03, 2.47352D-03, 2.53971D-03,
     $     2.58019D-03, 2.61639D-03, 2.71017D-03, 2.74649D-03,
     $     2.77220D-03, 2.81941D-03, 2.92251D-03, 2.97576D-03,
     $     3.00831D-03, 3.03328D-03, 3.12475D-03, 3.14229D-03,
     $     3.19203D-03, 3.26291D-03, 3.27259D-03, 3.33815D-03,
     $     3.41885D-03, 3.46915D-03, 3.53317D-03, 3.60416D-03,
     $     3.65557D-03, 3.68625D-03, 3.69608D-03, 3.74405D-03,
     $     3.79325D-03, 3.85516D-03, 3.88498D-03, 3.93649D-03,
     $     3.96564D-03, 4.03485D-03, 4.07756D-03, 4.16388D-03,
     $     4.17308D-03, 4.23294D-03, 4.30870D-03, 4.37248D-03,
     $     4.46085D-03, 4.49786D-03, 4.55044D-03, 4.61326D-03,
     $     4.67109D-03, 4.69234D-03, 4.77869D-03, 4.78868D-03,
     $     4.88584D-03, 4.91314D-03, 4.99402D-03, 5.05011D-03,
     $     5.10317D-03, 5.13934D-03, 5.21635D-03, 5.22948D-03,
     $     5.40282D-03, 5.48014D-03, 5.49753D-03, 5.51626D-03,
     $     5.63677D-03, 5.70093D-03, 5.74399D-03, 5.79333D-03,
     $     5.85519D-03, 5.86380D-03, 5.93664D-03, 5.96385D-03,
     $     6.02144D-03, 6.04887D-03, 6.16266D-03, 6.19705D-03,
     $     6.30627D-03, 6.31492D-03, 6.39003D-03, 6.45074D-03,
     $     6.48933D-03, 6.58080D-03, 6.59793D-03, 6.65497D-03,
     $     6.69938D-03, 6.72515D-03, 6.81485D-03, 6.85752D-03,
     $     6.90975D-03, 6.94636D-03, 6.98741D-03, 6.99876D-03,
     $     7.08380D-03, 7.11069D-03, 7.15008D-03, 7.19457D-03,
     $     7.22544D-03, 7.33465D-03, 7.40645D-03, 7.41478D-03,
     $     7.49309D-03, 7.54522D-03, 7.62401D-03, 7.63956D-03,
     $     7.69467D-03, 7.72701D-03, 7.75386D-03, 7.79353D-03,
     $     7.86310D-03, 7.88867D-03, 7.96971D-03, 8.00211D-03,
     $     8.06123D-03, 8.14925D-03, 8.17999D-03, 8.20946D-03,
     $     8.29007D-03, 8.31359D-03, 8.35228D-03, 8.36751D-03,
     $     8.43276D-03, 8.44786D-03, 8.55132D-03, 8.63107D-03,
     $     8.67910D-03, 8.75434D-03, 8.77906D-03, 8.80109D-03,
     $     8.85767D-03, 8.89450D-03, 8.97493D-03, 8.99534D-03,
     $     9.03338D-03, 9.08347D-03, 9.10102D-03, 9.19348D-03 /

!     Branching ratio of 5009
      data( br(18,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00 /

!     Branching ratio of 4009
      data( br(19,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 8.61523D-17, 1.52951D-13,
     $     8.33858D-13, 1.50471D-11, 2.20689D-10, 1.70302D-09,
     $     5.14343D-09, 1.55983D-08, 3.28674D-08, 5.40568D-08,
     $     1.35213D-07, 1.91002D-07, 1.52851D-06, 6.04265D-06,
     $     1.99231D-05, 4.53995D-05, 5.16036D-05, 8.79014D-05,
     $     1.65405D-04, 1.85049D-04, 2.96712D-04, 3.57131D-04,
     $     5.02446D-04, 5.96706D-04, 7.12647D-04, 1.17870D-03,
     $     1.23397D-03, 1.30456D-03, 1.53775D-03, 1.86777D-03,
     $     2.13860D-03, 2.39264D-03, 2.88737D-03, 3.63888D-03,
     $     3.88016D-03, 4.30558D-03, 4.37156D-03, 4.66322D-03,
     $     4.85953D-03, 5.08096D-03, 5.56255D-03, 5.92812D-03,
     $     6.29461D-03, 6.42262D-03, 6.46409D-03, 6.89032D-03,
     $     7.47063D-03, 7.81345D-03, 8.15559D-03, 8.47582D-03,
     $     8.66861D-03, 8.95410D-03, 9.10054D-03, 9.70480D-03,
     $     9.85388D-03, 1.00884D-02, 1.05808D-02, 1.07949D-02,
     $     1.09035D-02, 1.11922D-02, 1.13346D-02, 1.16694D-02,
     $     1.19335D-02, 1.19662D-02, 1.21957D-02, 1.23822D-02,
     $     1.25905D-02, 1.27954D-02, 1.31030D-02, 1.33080D-02,
     $     1.35257D-02, 1.36616D-02, 1.37740D-02, 1.39636D-02,
     $     1.42431D-02, 1.43683D-02, 1.46405D-02, 1.48739D-02,
     $     1.49847D-02, 1.51662D-02, 1.53008D-02, 1.53534D-02,
     $     1.54398D-02, 1.55384D-02, 1.57357D-02, 1.59000D-02,
     $     1.59196D-02, 1.60359D-02, 1.62243D-02, 1.64374D-02,
     $     1.65441D-02, 1.66777D-02, 1.68525D-02, 1.70075D-02,
     $     1.70493D-02, 1.72065D-02, 1.72412D-02, 1.74185D-02,
     $     1.74694D-02, 1.76238D-02, 1.77609D-02, 1.78151D-02,
     $     1.79401D-02, 1.79840D-02, 1.80729D-02, 1.81482D-02,
     $     1.82749D-02, 1.83665D-02, 1.84350D-02, 1.85723D-02,
     $     1.85807D-02, 1.87537D-02, 1.87748D-02, 1.88239D-02,
     $     1.89027D-02, 1.90119D-02, 1.90706D-02, 1.91175D-02,
     $     1.91525D-02, 1.91693D-02, 1.92199D-02, 1.92656D-02,
     $     1.93251D-02, 1.93622D-02, 1.93931D-02, 1.95020D-02,
     $     1.95271D-02, 1.95548D-02, 1.95908D-02, 1.96242D-02,
     $     1.96846D-02, 1.97171D-02, 1.97381D-02, 1.98007D-02,
     $     1.98145D-02, 1.98484D-02, 1.98931D-02, 1.99386D-02,
     $     1.99667D-02, 1.99927D-02, 2.00627D-02, 2.00902D-02,
     $     2.01099D-02, 2.01447D-02, 2.02183D-02, 2.02557D-02,
     $     2.02776D-02, 2.02944D-02, 2.03568D-02, 2.03680D-02,
     $     2.03999D-02, 2.04440D-02, 2.04499D-02, 2.04884D-02,
     $     2.05363D-02, 2.05657D-02, 2.06034D-02, 2.06425D-02,
     $     2.06718D-02, 2.06898D-02, 2.06957D-02, 2.07243D-02,
     $     2.07525D-02, 2.07880D-02, 2.08047D-02, 2.08351D-02,
     $     2.08529D-02, 2.08954D-02, 2.09214D-02, 2.09746D-02,
     $     2.09803D-02, 2.10193D-02, 2.10663D-02, 2.11095D-02,
     $     2.11694D-02, 2.11951D-02, 2.12330D-02, 2.12799D-02,
     $     2.13241D-02, 2.13406D-02, 2.14078D-02, 2.14158D-02,
     $     2.14946D-02, 2.15174D-02, 2.15864D-02, 2.16338D-02,
     $     2.16788D-02, 2.17107D-02, 2.17808D-02, 2.17930D-02,
     $     2.19558D-02, 2.20315D-02, 2.20486D-02, 2.20668D-02,
     $     2.21897D-02, 2.22582D-02, 2.23036D-02, 2.23579D-02,
     $     2.24267D-02, 2.24364D-02, 2.25185D-02, 2.25499D-02,
     $     2.26183D-02, 2.26515D-02, 2.27925D-02, 2.28362D-02,
     $     2.29788D-02, 2.29904D-02, 2.30925D-02, 2.31771D-02,
     $     2.32312D-02, 2.33631D-02, 2.33882D-02, 2.34727D-02,
     $     2.35396D-02, 2.35790D-02, 2.37183D-02, 2.37862D-02,
     $     2.38703D-02, 2.39298D-02, 2.39974D-02, 2.40163D-02,
     $     2.41593D-02, 2.42053D-02, 2.42730D-02, 2.43499D-02,
     $     2.44042D-02, 2.45995D-02, 2.47313D-02, 2.47466D-02,
     $     2.48931D-02, 2.49924D-02, 2.51441D-02, 2.51744D-02,
     $     2.52825D-02, 2.53464D-02, 2.53999D-02, 2.54794D-02,
     $     2.56205D-02, 2.56727D-02, 2.58398D-02, 2.59070D-02,
     $     2.60307D-02, 2.62168D-02, 2.62823D-02, 2.63452D-02,
     $     2.65189D-02, 2.65701D-02, 2.66542D-02, 2.66874D-02,
     $     2.68300D-02, 2.68630D-02, 2.70899D-02, 2.72648D-02,
     $     2.73705D-02, 2.75355D-02, 2.75898D-02, 2.76383D-02,
     $     2.77627D-02, 2.78438D-02, 2.80220D-02, 2.80673D-02,
     $     2.81519D-02, 2.82639D-02, 2.83033D-02, 2.85110D-02 /

!     Branching ratio of 3009
      data( br(20,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00 /

!     Branching ratio of 4008
      data( br(21,k), k=1, 300) /
     $     0.00000D+00, 4.11986D-13, 2.46349D-09, 5.40547D-09,
     $     1.02892D-07, 2.24347D-07, 1.85774D-06, 2.33932D-06,
     $     3.02993D-06, 3.35993D-06, 3.89168D-06, 4.05401D-06,
     $     5.89210D-06, 6.93357D-06, 8.41709D-06, 8.78268D-06,
     $     9.14140D-06, 9.38889D-06, 1.00111D-05, 1.12394D-05,
     $     1.13648D-05, 1.18566D-05, 1.27045D-05, 1.30303D-05,
     $     1.31999D-05, 1.37221D-05, 1.42165D-05, 1.52359D-05,
     $     1.58609D-05, 1.61342D-05, 1.66793D-05, 1.71268D-05,
     $     1.82020D-05, 2.32690D-05, 3.28605D-05, 6.88501D-05,
     $     8.70862D-05, 1.35205D-04, 2.31377D-04, 3.81440D-04,
     $     5.19563D-04, 7.20066D-04, 8.92464D-04, 1.01866D-03,
     $     1.20483D-03, 1.26057D-03, 1.61152D-03, 1.87862D-03,
     $     2.15233D-03, 2.37395D-03, 2.41155D-03, 2.58376D-03,
     $     2.82379D-03, 2.87070D-03, 3.09724D-03, 3.19706D-03,
     $     3.39898D-03, 3.51461D-03, 3.64432D-03, 4.07841D-03,
     $     4.12318D-03, 4.17843D-03, 4.35044D-03, 4.57232D-03,
     $     4.74118D-03, 4.89100D-03, 5.16498D-03, 5.58228D-03,
     $     5.72424D-03, 5.98878D-03, 6.03151D-03, 6.22691D-03,
     $     6.36783D-03, 6.53412D-03, 6.92178D-03, 7.23629D-03,
     $     7.56656D-03, 7.68585D-03, 7.72495D-03, 8.13464D-03,
     $     8.71147D-03, 9.06443D-03, 9.41848D-03, 9.75117D-03,
     $     9.95016D-03, 1.02439D-02, 1.03939D-02, 1.10098D-02,
     $     1.11625D-02, 1.14042D-02, 1.19164D-02, 1.21436D-02,
     $     1.22604D-02, 1.25746D-02, 1.27315D-02, 1.31085D-02,
     $     1.34184D-02, 1.34578D-02, 1.37423D-02, 1.39859D-02,
     $     1.42734D-02, 1.45694D-02, 1.50271D-02, 1.53417D-02,
     $     1.56812D-02, 1.58956D-02, 1.60749D-02, 1.63747D-02,
     $     1.68121D-02, 1.70080D-02, 1.74386D-02, 1.78118D-02,
     $     1.79909D-02, 1.82877D-02, 1.85088D-02, 1.85953D-02,
     $     1.87375D-02, 1.89005D-02, 1.92323D-02, 1.95106D-02,
     $     1.95438D-02, 1.97432D-02, 2.00706D-02, 2.04527D-02,
     $     2.06505D-02, 2.09035D-02, 2.12488D-02, 2.15714D-02,
     $     2.16599D-02, 2.20024D-02, 2.20803D-02, 2.24969D-02,
     $     2.26200D-02, 2.30055D-02, 2.33713D-02, 2.35236D-02,
     $     2.38803D-02, 2.40082D-02, 2.42711D-02, 2.45009D-02,
     $     2.49036D-02, 2.52163D-02, 2.54584D-02, 2.59536D-02,
     $     2.59848D-02, 2.66506D-02, 2.67374D-02, 2.69386D-02,
     $     2.72707D-02, 2.77429D-02, 2.80066D-02, 2.82311D-02,
     $     2.84035D-02, 2.84968D-02, 2.87647D-02, 2.90126D-02,
     $     2.93388D-02, 2.95466D-02, 2.97218D-02, 3.03513D-02,
     $     3.04987D-02, 3.06626D-02, 3.08807D-02, 3.10998D-02,
     $     3.15027D-02, 3.17279D-02, 3.18695D-02, 3.22968D-02,
     $     3.23949D-02, 3.26400D-02, 3.29668D-02, 3.33010D-02,
     $     3.35052D-02, 3.36885D-02, 3.41625D-02, 3.43439D-02,
     $     3.44713D-02, 3.47013D-02, 3.51870D-02, 3.54275D-02,
     $     3.55699D-02, 3.56773D-02, 3.60589D-02, 3.61293D-02,
     $     3.63264D-02, 3.66004D-02, 3.66373D-02, 3.68837D-02,
     $     3.71876D-02, 3.73772D-02, 3.76201D-02, 3.78883D-02,
     $     3.80844D-02, 3.82020D-02, 3.82397D-02, 3.84241D-02,
     $     3.86109D-02, 3.88439D-02, 3.89544D-02, 3.91444D-02,
     $     3.92514D-02, 3.95024D-02, 3.96544D-02, 3.99545D-02,
     $     3.99859D-02, 4.01903D-02, 4.04411D-02, 4.06510D-02,
     $     4.09341D-02, 4.10506D-02, 4.12148D-02, 4.14089D-02,
     $     4.15842D-02, 4.16480D-02, 4.19014D-02, 4.19304D-02,
     $     4.22072D-02, 4.22838D-02, 4.25079D-02, 4.26584D-02,
     $     4.27986D-02, 4.28937D-02, 4.30948D-02, 4.31289D-02,
     $     4.35689D-02, 4.37620D-02, 4.38047D-02, 4.38507D-02,
     $     4.41487D-02, 4.43083D-02, 4.44135D-02, 4.45361D-02,
     $     4.46898D-02, 4.47112D-02, 4.48917D-02, 4.49595D-02,
     $     4.51040D-02, 4.51731D-02, 4.54595D-02, 4.55452D-02,
     $     4.58202D-02, 4.58420D-02, 4.60309D-02, 4.61830D-02,
     $     4.62792D-02, 4.65074D-02, 4.65500D-02, 4.66913D-02,
     $     4.68010D-02, 4.68645D-02, 4.70844D-02, 4.71887D-02,
     $     4.73142D-02, 4.74023D-02, 4.75008D-02, 4.75280D-02,
     $     4.77300D-02, 4.77931D-02, 4.78857D-02, 4.79884D-02,
     $     4.80601D-02, 4.83117D-02, 4.84759D-02, 4.84949D-02,
     $     4.86719D-02, 4.87899D-02, 4.89675D-02, 4.90024D-02,
     $     4.91258D-02, 4.91981D-02, 4.92579D-02, 4.93459D-02,
     $     4.94997D-02, 4.95559D-02, 4.97333D-02, 4.98042D-02,
     $     4.99331D-02, 5.01242D-02, 5.01909D-02, 5.02548D-02,
     $     5.04284D-02, 5.04785D-02, 5.05615D-02, 5.05944D-02,
     $     5.07354D-02, 5.07682D-02, 5.09899D-02, 5.11613D-02,
     $     5.12630D-02, 5.14235D-02, 5.14760D-02, 5.15225D-02,
     $     5.16414D-02, 5.17192D-02, 5.18862D-02, 5.19284D-02,
     $     5.20070D-02, 5.21096D-02, 5.21451D-02, 5.23321D-02 /

!     Branching ratio of 3008
      data( br(22,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00 /

!     Branching ratio of 3007
      data( br(23,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 3.63383D-16, 6.66597D-13, 5.06052D-12,
     $     6.00040D-11, 1.19473D-08, 7.99071D-08, 5.14257D-07,
     $     7.83995D-07, 1.68566D-06, 4.24454D-06, 1.03873D-05,
     $     1.81870D-05, 3.27489D-05, 4.82841D-05, 6.13089D-05,
     $     8.35112D-05, 9.08090D-05, 1.45398D-04, 1.99215D-04,
     $     2.69323D-04, 3.40134D-04, 3.53612D-04, 4.21609D-04,
     $     5.34439D-04, 5.59715D-04, 6.97216D-04, 7.68606D-04,
     $     9.34003D-04, 1.04045D-03, 1.16954D-03, 1.66368D-03,
     $     1.71961D-03, 1.78987D-03, 2.01512D-03, 2.32147D-03,
     $     2.56565D-03, 2.79278D-03, 3.24084D-03, 4.00036D-03,
     $     4.27155D-03, 4.77784D-03, 4.85940D-03, 5.23527D-03,
     $     5.50417D-03, 5.81806D-03, 6.54290D-03, 7.12653D-03,
     $     7.72895D-03, 7.94394D-03, 8.01409D-03, 8.74537D-03,
     $     9.76311D-03, 1.03807D-02, 1.10070D-02, 1.16026D-02,
     $     1.19682D-02, 1.25233D-02, 1.28153D-02, 1.40885D-02,
     $     1.44207D-02, 1.49622D-02, 1.61619D-02, 1.67218D-02,
     $     1.70176D-02, 1.78480D-02, 1.82741D-02, 1.93388D-02,
     $     2.02361D-02, 2.03503D-02, 2.11773D-02, 2.18821D-02,
     $     2.27143D-02, 2.35623D-02, 2.48614D-02, 2.57483D-02,
     $     2.67101D-02, 2.73190D-02, 2.78330D-02, 2.86975D-02,
     $     2.99742D-02, 3.05508D-02, 3.18488D-02, 3.29797D-02,
     $     3.35242D-02, 3.44353D-02, 3.51131D-02, 3.53793D-02,
     $     3.58193D-02, 3.63259D-02, 3.73744D-02, 3.82545D-02,
     $     3.83598D-02, 3.89975D-02, 4.00472D-02, 4.12500D-02,
     $     4.18642D-02, 4.26353D-02, 4.36605D-02, 4.45902D-02,
     $     4.48410D-02, 4.57934D-02, 4.60066D-02, 4.71188D-02,
     $     4.74390D-02, 4.84185D-02, 4.93195D-02, 4.96862D-02,
     $     5.05269D-02, 5.08227D-02, 5.14249D-02, 5.19453D-02,
     $     5.28517D-02, 5.35503D-02, 5.40859D-02, 5.51735D-02,
     $     5.52417D-02, 5.66661D-02, 5.68488D-02, 5.72703D-02,
     $     5.79610D-02, 5.89262D-02, 5.94564D-02, 5.99017D-02,
     $     6.02409D-02, 6.04227D-02, 6.09412D-02, 6.14150D-02,
     $     6.20328D-02, 6.24234D-02, 6.27515D-02, 6.39214D-02,
     $     6.41928D-02, 6.44930D-02, 6.48892D-02, 6.52831D-02,
     $     6.60003D-02, 6.63959D-02, 6.66438D-02, 6.73836D-02,
     $     6.75510D-02, 6.79663D-02, 6.85148D-02, 6.90698D-02,
     $     6.94064D-02, 6.97076D-02, 7.04847D-02, 7.07828D-02,
     $     7.09921D-02, 7.13691D-02, 7.21670D-02, 7.25671D-02,
     $     7.28052D-02, 7.29857D-02, 7.36360D-02, 7.37562D-02,
     $     7.40960D-02, 7.45721D-02, 7.46364D-02, 7.50660D-02,
     $     7.55914D-02, 7.59156D-02, 7.63285D-02, 7.67748D-02,
     $     7.70964D-02, 7.72878D-02, 7.73489D-02, 7.76452D-02,
     $     7.79411D-02, 7.83073D-02, 7.84794D-02, 7.87757D-02,
     $     7.89418D-02, 7.93282D-02, 7.95593D-02, 8.00129D-02,
     $     8.00598D-02, 8.03664D-02, 8.07331D-02, 8.10423D-02,
     $     8.14535D-02, 8.16220D-02, 8.18612D-02, 8.21450D-02,
     $     8.24025D-02, 8.24962D-02, 8.28665D-02, 8.29092D-02,
     $     8.33129D-02, 8.34251D-02, 8.37540D-02, 8.39724D-02,
     $     8.41725D-02, 8.43094D-02, 8.45978D-02, 8.46465D-02,
     $     8.52584D-02, 8.55180D-02, 8.55737D-02, 8.56320D-02,
     $     8.60065D-02, 8.62020D-02, 8.63227D-02, 8.64646D-02,
     $     8.66337D-02, 8.66566D-02, 8.68431D-02, 8.69124D-02,
     $     8.70605D-02, 8.71308D-02, 8.74122D-02, 8.74931D-02,
     $     8.77463D-02, 8.77663D-02, 8.79385D-02, 8.80751D-02,
     $     8.81563D-02, 8.83468D-02, 8.83817D-02, 8.84962D-02,
     $     8.85863D-02, 8.86382D-02, 8.88112D-02, 8.88924D-02,
     $     8.89866D-02, 8.90509D-02, 8.91242D-02, 8.91445D-02,
     $     8.92929D-02, 8.93385D-02, 8.94019D-02, 8.94664D-02,
     $     8.95137D-02, 8.96714D-02, 8.97732D-02, 8.97836D-02,
     $     8.98827D-02, 8.99496D-02, 9.00433D-02, 9.00617D-02,
     $     9.01243D-02, 9.01582D-02, 9.01869D-02, 9.02286D-02,
     $     9.02995D-02, 9.03246D-02, 9.04001D-02, 9.04260D-02,
     $     9.04707D-02, 9.05260D-02, 9.05443D-02, 9.05623D-02,
     $     9.06020D-02, 9.06120D-02, 9.06271D-02, 9.06332D-02,
     $     9.06600D-02, 9.06641D-02, 9.06941D-02, 9.07103D-02,
     $     9.07175D-02, 9.07245D-02, 9.07273D-02, 9.07301D-02,
     $     9.07376D-02, 9.07378D-02, 9.07361D-02, 9.07350D-02,
     $     9.07324D-02, 9.07265D-02, 9.07237D-02, 9.07136D-02 /

!     Branching ratio of 3006
      data( br(24,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 4.99006D-37,
     $     4.01863D-36, 1.78447D-35, 1.26411D-30, 3.36568D-28,
     $     1.26976D-26, 2.91349D-26, 3.64506D-26, 5.11164D-25,
     $     5.42574D-23, 3.45670D-22, 1.36612D-21, 4.97870D-21,
     $     9.50866D-21, 2.57128D-20, 5.32927D-20, 3.40949D-19,
     $     7.48828D-18, 8.82975D-17, 1.69103D-15, 1.00426D-14,
     $     1.50989D-14, 2.52343D-13, 5.50434D-13, 2.79155D-12,
     $     7.80843D-12, 8.79501D-12, 1.89625D-11, 3.30039D-11,
     $     5.90712D-11, 9.99322D-11, 2.11139D-10, 3.32247D-10,
     $     5.16842D-10, 6.68111D-10, 8.17176D-10, 1.11623D-09,
     $     1.67646D-09, 1.97968D-09, 2.83730D-09, 3.90183D-09,
     $     4.57270D-09, 6.05868D-09, 7.54885D-09, 8.23531D-09,
     $     9.50917D-09, 1.12141D-08, 1.58675D-08, 2.07795D-08,
     $     2.14301D-08, 2.59708D-08, 3.50934D-08, 4.74367D-08,
     $     5.52238D-08, 6.50483D-08, 7.92657D-08, 9.80921D-08,
     $     1.03014D-07, 1.22403D-07, 1.27009D-07, 1.65788D-07,
     $     1.75623D-07, 2.07324D-07, 2.55778D-07, 2.90146D-07,
     $     3.50797D-07, 3.69292D-07, 4.05890D-07, 4.38217D-07,
     $     5.26634D-07, 6.80241D-07, 7.63695D-07, 1.01689D-06,
     $     1.03886D-06, 1.65573D-06, 1.80603D-06, 2.11365D-06,
     $     2.63529D-06, 3.38821D-06, 3.86199D-06, 4.42483D-06,
     $     4.84673D-06, 5.30386D-06, 6.36691D-06, 7.29376D-06,
     $     8.66515D-06, 9.65208D-06, 1.05772D-05, 1.47345D-05,
     $     1.59167D-05, 1.73394D-05, 1.94228D-05, 2.20001D-05,
     $     2.78194D-05, 3.15257D-05, 3.37748D-05, 4.11857D-05,
     $     4.31715D-05, 4.81221D-05, 5.51532D-05, 6.38061D-05,
     $     6.95949D-05, 7.49943D-05, 9.00858D-05, 9.64529D-05,
     $     1.01147D-04, 1.10516D-04, 1.33402D-04, 1.45930D-04,
     $     1.54077D-04, 1.60295D-04, 1.82961D-04, 1.87752D-04,
     $     2.01074D-04, 2.20567D-04, 2.23284D-04, 2.41490D-04,
     $     2.64723D-04, 2.79743D-04, 2.98859D-04, 3.22172D-04,
     $     3.39121D-04, 3.49168D-04, 3.52387D-04, 3.68225D-04,
     $     3.85645D-04, 4.08277D-04, 4.19664D-04, 4.38810D-04,
     $     4.49575D-04, 4.75567D-04, 4.92495D-04, 5.26110D-04,
     $     5.29735D-04, 5.53431D-04, 5.86897D-04, 6.14169D-04,
     $     6.54318D-04, 6.71832D-04, 6.96399D-04, 7.25720D-04,
     $     7.53035D-04, 7.63177D-04, 8.07003D-04, 8.11985D-04,
     $     8.62567D-04, 8.76826D-04, 9.19366D-04, 9.50072D-04,
     $     9.81363D-04, 1.00204D-03, 1.04599D-03, 1.05352D-03,
     $     1.15730D-03, 1.20569D-03, 1.21698D-03, 1.22986D-03,
     $     1.30928D-03, 1.35181D-03, 1.38258D-03, 1.41672D-03,
     $     1.46272D-03, 1.46929D-03, 1.52611D-03, 1.54742D-03,
     $     1.59171D-03, 1.61275D-03, 1.70342D-03, 1.73120D-03,
     $     1.82261D-03, 1.82991D-03, 1.89357D-03, 1.94585D-03,
     $     1.98185D-03, 2.06730D-03, 2.08365D-03, 2.13908D-03,
     $     2.18169D-03, 2.20659D-03, 2.29675D-03, 2.33995D-03,
     $     2.39357D-03, 2.43301D-03, 2.47660D-03, 2.48862D-03,
     $     2.57937D-03, 2.60837D-03, 2.65300D-03, 2.70519D-03,
     $     2.74042D-03, 2.86832D-03, 2.95220D-03, 2.96287D-03,
     $     3.05808D-03, 3.12086D-03, 3.21909D-03, 3.23834D-03,
     $     3.30745D-03, 3.34947D-03, 3.38368D-03, 3.43411D-03,
     $     3.52277D-03, 3.55557D-03, 3.66022D-03, 3.70447D-03,
     $     3.78615D-03, 3.91125D-03, 3.95493D-03, 3.99617D-03,
     $     4.11218D-03, 4.14605D-03, 4.20303D-03, 4.22546D-03,
     $     4.32021D-03, 4.34343D-03, 4.49659D-03, 4.61926D-03,
     $     4.69275D-03, 4.81220D-03, 4.85086D-03, 4.88501D-03,
     $     4.97227D-03, 5.03237D-03, 5.16193D-03, 5.19521D-03,
     $     5.25784D-03, 5.34093D-03, 5.37018D-03, 5.52167D-03 /

!     Branching ratio of 2006
      data( br(25,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 6.85613D-37, 1.54895D-55, 2.08698D-54,
     $     4.27556D-50, 8.73379D-32, 1.52613D-32, 2.42691D-28,
     $     1.27112D-27, 4.62716D-27, 1.27367D-26, 5.38963D-26,
     $     3.66579D-23, 1.39836D-22, 6.94679D-22, 3.30990D-20,
     $     1.86451D-19, 4.78382D-18, 5.03334D-17, 2.55513D-16,
     $     5.89462D-16, 1.25420D-15, 1.35856D-14, 2.60973D-14,
     $     4.26844D-14, 1.26456D-13, 6.35968D-13, 1.23063D-12,
     $     1.77084D-12, 2.32134D-12, 5.19252D-12, 5.97088D-12,
     $     8.54745D-12, 1.35718D-11, 1.44298D-11, 2.11481D-11,
     $     3.32846D-11, 4.56449D-11, 7.03342D-11, 1.25727D-10,
     $     1.97467D-10, 2.62712D-10, 2.88596D-10, 4.52630D-10,
     $     7.11323D-10, 1.18971D-09, 1.49749D-09, 2.17550D-09,
     $     2.64594D-09, 4.07013D-09, 5.21040D-09, 8.14013D-09,
     $     8.51353D-09, 1.12322D-08, 1.55604D-08, 2.00613D-08,
     $     2.78104D-08, 3.17471D-08, 3.80360D-08, 4.66526D-08,
     $     5.59067D-08, 5.96169D-08, 7.65579D-08, 7.87196D-08,
     $     1.02256D-07, 1.09649D-07, 1.33947D-07, 1.52951D-07,
     $     1.72731D-07, 1.87254D-07, 2.21108D-07, 2.27326D-07,
     $     3.24122D-07, 3.77056D-07, 3.89978D-07, 4.04379D-07,
     $     5.07355D-07, 5.70478D-07, 6.17021D-07, 6.73878D-07,
     $     7.51582D-07, 7.62986D-07, 8.66460D-07, 9.07856D-07,
     $     1.00075D-06, 1.04770D-06, 1.26338D-06, 1.33559D-06,
     $     1.58857D-06, 1.61019D-06, 1.80765D-06, 1.98020D-06,
     $     2.09699D-06, 2.39328D-06, 2.45223D-06, 2.65708D-06,
     $     2.82483D-06, 2.92576D-06, 3.30013D-06, 3.49011D-06,
     $     3.73510D-06, 3.91486D-06, 4.12419D-06, 4.18360D-06,
     $     4.65220D-06, 4.80884D-06, 5.04683D-06, 5.32948D-06,
     $     5.53021D-06, 6.29305D-06, 6.83796D-06, 6.90439D-06,
     $     7.55069D-06, 8.00439D-06, 8.73581D-06, 8.88616D-06,
     $     9.43562D-06, 9.77144D-06, 1.00561D-05, 1.04885D-05,
     $     1.12803D-05, 1.15821D-05, 1.25789D-05, 1.29970D-05,
     $     1.37875D-05, 1.50372D-05, 1.54915D-05, 1.59350D-05,
     $     1.72023D-05, 1.75865D-05, 1.82333D-05, 1.84921D-05,
     $     1.96260D-05, 1.98976D-05, 2.18269D-05, 2.34102D-05,
     $     2.44055D-05, 2.60294D-05, 2.65775D-05, 2.70719D-05,
     $     2.83694D-05, 2.92445D-05, 3.12211D-05, 3.17379D-05,
     $     3.27179D-05, 3.40438D-05, 3.45188D-05, 3.70840D-05 /

!     Branching ratio of 2005
      data( br(26,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00 /

!     Branching ratio of 2004
      data( br(27,k), k=1, 300) /
     $     0.00000D+00, 9.65093D-22, 2.44973D-10, 6.39384D-10,
     $     6.47445D-07, 1.90843D-06, 1.96517D-05, 2.45616D-05,
     $     2.47952D-05, 2.40675D-05, 2.16747D-05, 2.09545D-05,
     $     1.87781D-05, 1.80729D-05, 1.69016D-05, 1.65529D-05,
     $     1.62833D-05, 1.60929D-05, 1.55313D-05, 1.48471D-05,
     $     1.48210D-05, 1.51079D-05, 2.01381D-05, 2.92835D-05,
     $     3.81494D-05, 1.00625D-04, 2.47754D-04, 8.35736D-04,
     $     1.23334D-03, 1.41308D-03, 1.73041D-03, 1.92052D-03,
     $     2.21526D-03, 2.64886D-03, 2.99094D-03, 3.56566D-03,
     $     3.76359D-03, 4.17908D-03, 4.85709D-03, 5.69455D-03,
     $     6.34929D-03, 7.18968D-03, 7.83318D-03, 8.27569D-03,
     $     8.89633D-03, 9.07852D-03, 1.01915D-02, 1.10337D-02,
     $     1.19379D-02, 1.27166D-02, 1.28514D-02, 1.34704D-02,
     $     1.43415D-02, 1.45127D-02, 1.53602D-02, 1.57400D-02,
     $     1.65103D-02, 1.69456D-02, 1.74239D-02, 1.89434D-02,
     $     1.90946D-02, 1.92801D-02, 1.98475D-02, 2.05485D-02,
     $     2.10597D-02, 2.14990D-02, 2.22675D-02, 2.33205D-02,
     $     2.36322D-02, 2.41581D-02, 2.42370D-02, 2.45698D-02,
     $     2.47945D-02, 2.50464D-02, 2.55806D-02, 2.59533D-02,
     $     2.63065D-02, 2.64264D-02, 2.64648D-02, 2.68341D-02,
     $     2.72859D-02, 2.75375D-02, 2.77680D-02, 2.79768D-02,
     $     2.80902D-02, 2.82596D-02, 2.83447D-02, 2.86634D-02,
     $     2.87424D-02, 2.88647D-02, 2.91086D-02, 2.92122D-02,
     $     2.92647D-02, 2.93897D-02, 2.94553D-02, 2.96062D-02,
     $     2.97283D-02, 2.97439D-02, 2.98492D-02, 2.99378D-02,
     $     3.00294D-02, 3.01212D-02, 3.02739D-02, 3.03783D-02,
     $     3.04961D-02, 3.05731D-02, 3.06333D-02, 3.07377D-02,
     $     3.09005D-02, 3.09757D-02, 3.11297D-02, 3.12748D-02,
     $     3.13480D-02, 3.14668D-02, 3.15599D-02, 3.15969D-02,
     $     3.16587D-02, 3.17304D-02, 3.18699D-02, 3.19926D-02,
     $     3.20077D-02, 3.20949D-02, 3.22408D-02, 3.24161D-02,
     $     3.25064D-02, 3.26269D-02, 3.27939D-02, 3.29459D-02,
     $     3.29881D-02, 3.31518D-02, 3.31889D-02, 3.33808D-02,
     $     3.34380D-02, 3.36165D-02, 3.37795D-02, 3.38456D-02,
     $     3.40051D-02, 3.40635D-02, 3.41848D-02, 3.42919D-02,
     $     3.44801D-02, 3.46216D-02, 3.47317D-02, 3.49632D-02,
     $     3.49779D-02, 3.52892D-02, 3.53286D-02, 3.54216D-02,
     $     3.55761D-02, 3.58015D-02, 3.59287D-02, 3.60343D-02,
     $     3.61158D-02, 3.61568D-02, 3.62812D-02, 3.63973D-02,
     $     3.65530D-02, 3.66535D-02, 3.67390D-02, 3.70520D-02,
     $     3.71264D-02, 3.72097D-02, 3.73207D-02, 3.74299D-02,
     $     3.76345D-02, 3.77498D-02, 3.78244D-02, 3.80533D-02,
     $     3.81060D-02, 3.82384D-02, 3.84173D-02, 3.86027D-02,
     $     3.87179D-02, 3.88227D-02, 3.90998D-02, 3.92079D-02,
     $     3.92845D-02, 3.94229D-02, 3.97191D-02, 3.98689D-02,
     $     3.99580D-02, 4.00261D-02, 4.02756D-02, 4.03219D-02,
     $     4.04540D-02, 4.06411D-02, 4.06667D-02, 4.08397D-02,
     $     4.10598D-02, 4.12002D-02, 4.13846D-02, 4.15910D-02,
     $     4.17459D-02, 4.18408D-02, 4.18715D-02, 4.20232D-02,
     $     4.21785D-02, 4.23757D-02, 4.24702D-02, 4.26361D-02,
     $     4.27310D-02, 4.29565D-02, 4.30946D-02, 4.33731D-02,
     $     4.34026D-02, 4.35975D-02, 4.38380D-02, 4.40442D-02,
     $     4.43265D-02, 4.44441D-02, 4.46122D-02, 4.48138D-02,
     $     4.49985D-02, 4.50662D-02, 4.53384D-02, 4.53700D-02,
     $     4.56760D-02, 4.57622D-02, 4.60175D-02, 4.61915D-02,
     $     4.63547D-02, 4.64673D-02, 4.67085D-02, 4.67498D-02,
     $     4.72917D-02, 4.75355D-02, 4.75901D-02, 4.76485D-02,
     $     4.80348D-02, 4.82453D-02, 4.83848D-02, 4.85490D-02,
     $     4.87556D-02, 4.87845D-02, 4.90297D-02, 4.91226D-02,
     $     4.93223D-02, 4.94186D-02, 4.98219D-02, 4.99448D-02,
     $     5.03413D-02, 5.03730D-02, 5.06503D-02, 5.08760D-02,
     $     5.10187D-02, 5.13611D-02, 5.14255D-02, 5.16401D-02,
     $     5.18085D-02, 5.19065D-02, 5.22477D-02, 5.24115D-02,
     $     5.26115D-02, 5.27518D-02, 5.29101D-02, 5.29540D-02,
     $     5.32838D-02, 5.33882D-02, 5.35408D-02, 5.37118D-02,
     $     5.38319D-02, 5.42578D-02, 5.45407D-02, 5.45732D-02,
     $     5.48822D-02, 5.50899D-02, 5.54043D-02, 5.54667D-02,
     $     5.56881D-02, 5.58180D-02, 5.59266D-02, 5.60872D-02,
     $     5.63701D-02, 5.64744D-02, 5.68060D-02, 5.69386D-02,
     $     5.71817D-02, 5.75465D-02, 5.76752D-02, 5.77994D-02,
     $     5.81405D-02, 5.82406D-02, 5.84062D-02, 5.84718D-02,
     $     5.87559D-02, 5.88217D-02, 5.92770D-02, 5.96306D-02,
     $     5.98449D-02, 6.01826D-02, 6.02944D-02, 6.03942D-02,
     $     6.06512D-02, 6.08185D-02, 6.11856D-02, 6.12790D-02,
     $     6.14537D-02, 6.16840D-02, 6.17646D-02, 6.21918D-02 /

!     Branching ratio of 2003
      data( br(28,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     2.81480D-49, 2.93223D-43, 8.68321D-43, 1.17410D-40,
     $     4.38523D-35, 3.04210D-32, 1.20753D-31, 2.51378D-30,
     $     6.14455D-29, 2.55976D-28, 8.20165D-28, 1.35488D-26,
     $     2.77713D-26, 7.43094D-26, 8.51330D-25, 2.72714D-24,
     $     4.51042D-24, 3.87910D-23, 6.13967D-23, 4.54577D-21,
     $     3.28239D-20, 3.98269D-20, 1.64967D-19, 8.30537D-19,
     $     2.80694D-18, 1.51604D-17, 1.18943D-16, 4.47202D-16,
     $     1.42017D-15, 2.86021D-15, 4.89884D-15, 1.19874D-14,
     $     3.83588D-14, 6.27149D-14, 1.74306D-13, 4.01720D-13,
     $     5.83847D-13, 1.09140D-12, 1.82870D-12, 2.23406D-12,
     $     2.99030D-12, 4.42664D-12, 1.00129D-11, 2.02201D-11,
     $     2.18286D-11, 3.57236D-11, 7.58361D-11, 1.75231D-10,
     $     2.56798D-10, 4.20085D-10, 7.79457D-10, 1.28904D-09,
     $     1.46481D-09, 2.37382D-09, 2.64645D-09, 4.51046D-09,
     $     5.20214D-09, 7.85465D-09, 1.14290D-08, 1.32838D-08,
     $     1.83685D-08, 2.05074D-08, 2.54674D-08, 3.06247D-08,
     $     4.19959D-08, 5.34249D-08, 6.37469D-08, 9.14486D-08,
     $     9.35765D-08, 1.50222D-07, 1.59769D-07, 1.84224D-07,
     $     2.33066D-07, 3.21109D-07, 3.81619D-07, 4.41658D-07,
     $     4.93192D-07, 5.24510D-07, 6.21197D-07, 7.23066D-07,
     $     8.78306D-07, 9.90651D-07, 1.09456D-06, 1.55055D-06,
     $     1.67659D-06, 1.82609D-06, 2.04295D-06, 2.28683D-06,
     $     2.79604D-06, 3.11695D-06, 3.33022D-06, 4.04268D-06,
     $     4.22420D-06, 4.70694D-06, 5.41637D-06, 6.23149D-06,
     $     6.78087D-06, 7.30662D-06, 8.82130D-06, 9.47149D-06,
     $     9.95450D-06, 1.08924D-05, 1.31783D-05, 1.44854D-05,
     $     1.53258D-05, 1.59884D-05, 1.85590D-05, 1.90831D-05,
     $     2.06221D-05, 2.29568D-05, 2.32887D-05, 2.56159D-05,
     $     2.86758D-05, 3.07028D-05, 3.34261D-05, 3.66411D-05,
     $     3.90834D-05, 4.05826D-05, 4.10691D-05, 4.34932D-05,
     $     4.60823D-05, 4.94663D-05, 5.11420D-05, 5.40871D-05,
     $     5.57840D-05, 5.99012D-05, 6.25155D-05, 6.79517D-05,
     $     6.85436D-05, 7.24369D-05, 7.75570D-05, 8.19676D-05,
     $     8.83225D-05, 9.10627D-05, 9.50231D-05, 9.98644D-05,
     $     1.04426D-04, 1.06130D-04, 1.13213D-04, 1.14042D-04,
     $     1.22250D-04, 1.24596D-04, 1.31657D-04, 1.36668D-04,
     $     1.41485D-04, 1.44779D-04, 1.51868D-04, 1.53088D-04,
     $     1.69556D-04, 1.77071D-04, 1.78781D-04, 1.80636D-04,
     $     1.92641D-04, 1.99112D-04, 2.03536D-04, 2.08600D-04,
     $     2.15035D-04, 2.15938D-04, 2.23672D-04, 2.26595D-04,
     $     2.32820D-04, 2.35811D-04, 2.48451D-04, 2.52364D-04,
     $     2.65046D-04, 2.66069D-04, 2.75058D-04, 2.82459D-04,
     $     2.87263D-04, 2.98844D-04, 3.01051D-04, 3.08502D-04,
     $     3.14386D-04, 3.17839D-04, 3.30081D-04, 3.36003D-04,
     $     3.43380D-04, 3.48621D-04, 3.54546D-04, 3.56195D-04,
     $     3.68741D-04, 3.72770D-04, 3.78734D-04, 3.85565D-04,
     $     3.90302D-04, 4.07345D-04, 4.18762D-04, 4.20109D-04,
     $     4.32804D-04, 4.41329D-04, 4.54417D-04, 4.57017D-04,
     $     4.66295D-04, 4.71799D-04, 4.76379D-04, 4.83178D-04,
     $     4.95230D-04, 4.99698D-04, 5.13976D-04, 5.19776D-04,
     $     5.30463D-04, 5.46678D-04, 5.52390D-04, 5.57878D-04,
     $     5.73109D-04, 5.77606D-04, 5.85047D-04, 5.87984D-04,
     $     6.00613D-04, 6.03570D-04, 6.23979D-04, 6.40021D-04,
     $     6.49802D-04, 6.65320D-04, 6.70450D-04, 6.75032D-04,
     $     6.86862D-04, 6.94673D-04, 7.11853D-04, 7.16252D-04,
     $     7.24499D-04, 7.35449D-04, 7.39310D-04, 7.59701D-04 /

!     Branching ratio of 1003
      data( br(29,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     3.09637D-40, 1.30133D-21, 1.40228D-18, 5.85810D-16,
     $     1.88318D-15, 1.77033D-14, 2.72868D-14, 1.10738D-13,
     $     3.04871D-13, 7.71719D-13, 4.81828D-12, 1.48419D-11,
     $     4.03644D-11, 5.61402D-11, 6.20789D-11, 1.67086D-10,
     $     6.13740D-10, 1.21642D-09, 2.37637D-09, 4.40661D-09,
     $     6.21568D-09, 1.02464D-08, 1.32196D-08, 3.59678D-08,
     $     4.53837D-08, 6.63315D-08, 1.41683D-07, 1.97621D-07,
     $     2.33842D-07, 3.63200D-07, 4.46753D-07, 7.17205D-07,
     $     1.01329D-06, 1.05611D-06, 1.41181D-06, 1.77530D-06,
     $     2.28064D-06, 2.88525D-06, 4.00217D-06, 4.88994D-06,
     $     5.99235D-06, 6.78091D-06, 7.50756D-06, 8.85954D-06,
     $     1.12113D-05, 1.24334D-05, 1.56431D-05, 1.89491D-05,
     $     2.07155D-05, 2.39757D-05, 2.66449D-05, 2.77533D-05,
     $     2.96591D-05, 3.19764D-05, 3.72667D-05, 4.21793D-05,
     $     4.27961D-05, 4.67138D-05, 5.37348D-05, 6.26799D-05,
     $     6.76577D-05, 7.42666D-05, 8.37892D-05, 9.32510D-05,
     $     9.59289D-05, 1.06650D-04, 1.09184D-04, 1.23315D-04,
     $     1.27625D-04, 1.41627D-04, 1.55871D-04, 1.62064D-04,
     $     1.76937D-04, 1.82423D-04, 1.94046D-04, 2.04582D-04,
     $     2.24164D-04, 2.40633D-04, 2.53908D-04, 2.82463D-04,
     $     2.84341D-04, 3.26016D-04, 3.31794D-04, 3.45346D-04,
     $     3.68524D-04, 4.02673D-04, 4.22374D-04, 4.39728D-04,
     $     4.53318D-04, 4.61033D-04, 4.83036D-04, 5.03795D-04,
     $     5.31778D-04, 5.50101D-04, 5.65903D-04, 6.25617D-04,
     $     6.40277D-04, 6.56876D-04, 6.79488D-04, 7.03164D-04,
     $     7.48069D-04, 7.73972D-04, 7.90348D-04, 8.40822D-04,
     $     8.52732D-04, 8.83001D-04, 9.24367D-04, 9.68201D-04,
     $     9.95721D-04, 1.02082D-03, 1.08777D-03, 1.11439D-03,
     $     1.13344D-03, 1.16885D-03, 1.24782D-03, 1.28919D-03,
     $     1.31460D-03, 1.33410D-03, 1.40594D-03, 1.41986D-03,
     $     1.45954D-03, 1.51684D-03, 1.52475D-03, 1.57875D-03,
     $     1.64656D-03, 1.68964D-03, 1.74549D-03, 1.80877D-03,
     $     1.85530D-03, 1.88329D-03, 1.89229D-03, 1.93647D-03,
     $     1.98232D-03, 2.04059D-03, 2.06885D-03, 2.11775D-03,
     $     2.14548D-03, 2.21145D-03, 2.25230D-03, 2.33524D-03,
     $     2.34412D-03, 2.40185D-03, 2.47560D-03, 2.53767D-03,
     $     2.62456D-03, 2.66116D-03, 2.71329D-03, 2.77585D-03,
     $     2.83365D-03, 2.85496D-03, 2.94214D-03, 2.95223D-03,
     $     3.05067D-03, 3.07844D-03, 3.16101D-03, 3.21867D-03,
     $     3.27348D-03, 3.31077D-03, 3.39026D-03, 3.40384D-03,
     $     3.58416D-03, 3.66498D-03, 3.68322D-03, 3.70293D-03,
     $     3.82963D-03, 3.89723D-03, 3.94294D-03, 3.99520D-03,
     $     4.06106D-03, 4.07025D-03, 4.14833D-03, 4.17760D-03,
     $     4.23964D-03, 4.26926D-03, 4.39290D-03, 4.43061D-03,
     $     4.55120D-03, 4.56081D-03, 4.64456D-03, 4.71272D-03,
     $     4.75646D-03, 4.86083D-03, 4.88051D-03, 4.94642D-03,
     $     4.99800D-03, 5.02809D-03, 5.13369D-03, 5.18431D-03,
     $     5.24680D-03, 5.29088D-03, 5.34045D-03, 5.35419D-03,
     $     5.45780D-03, 5.49077D-03, 5.53929D-03, 5.59444D-03,
     $     5.63264D-03, 5.76853D-03, 5.85836D-03, 5.86887D-03,
     $     5.96747D-03, 6.03322D-03, 6.13306D-03, 6.15278D-03,
     $     6.22277D-03, 6.26401D-03, 6.29824D-03, 6.34883D-03,
     $     6.43775D-03, 6.47050D-03, 6.57442D-03, 6.61623D-03,
     $     6.69277D-03, 6.80756D-03, 6.84771D-03, 6.88618D-03,
     $     6.99209D-03, 7.02312D-03, 7.07427D-03, 7.09441D-03,
     $     7.18073D-03, 7.20080D-03, 7.33845D-03, 7.44525D-03,
     $     7.50984D-03, 7.61154D-03, 7.64502D-03, 7.67487D-03,
     $     7.75162D-03, 7.80192D-03, 7.91196D-03, 7.93998D-03,
     $     7.99236D-03, 8.06154D-03, 8.08583D-03, 8.21369D-03 /

!     Branching ratio of 1002
      data( br(30,k), k=1, 300) /
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 0.00000D+00, 0.00000D+00,
     $     0.00000D+00, 0.00000D+00, 3.33311D-94, 1.06600D-39,
     $     2.34715D-39, 3.33875D-39, 5.22888D-38, 1.17981D-34,
     $     2.82914D-34, 7.05667D-33, 1.39186D-32, 6.26403D-31,
     $     3.78568D-30, 7.04103D-30, 1.05228D-28, 6.58779D-28,
     $     2.39340D-27, 3.81748D-26, 8.81135D-25, 1.75282D-24,
     $     2.82277D-24, 8.01274D-24, 3.78047D-23, 7.40809D-23,
     $     1.12597D-22, 1.59635D-22, 5.03685D-22, 6.12087D-22,
     $     1.29187D-21, 3.90811D-21, 4.42994D-21, 1.06106D-20,
     $     2.82914D-20, 4.89319D-20, 9.57928D-20, 1.88356D-19,
     $     3.05340D-19, 4.00354D-19, 4.35708D-19, 6.52498D-19,
     $     9.74325D-19, 1.58453D-18, 1.97028D-18, 2.82867D-18,
     $     3.46041D-18, 5.40523D-18, 7.00948D-18, 1.14436D-17,
     $     1.20263D-17, 1.64006D-17, 2.37070D-17, 3.15158D-17,
     $     4.57014D-17, 5.31981D-17, 6.57143D-17, 8.39724D-17,
     $     1.06778D-16, 1.17287D-16, 1.74549D-16, 1.83756D-16,
     $     3.32937D-16, 4.02857D-16, 7.82961D-16, 1.25418D-15,
     $     2.03772D-15, 2.82841D-15, 5.55330D-15, 6.20995D-15,
     $     2.52226D-14, 4.58368D-14, 5.28209D-14, 6.24618D-14,
     $     1.42985D-13, 2.16179D-13, 2.94536D-13, 4.03061D-13,
     $     6.41309D-13, 6.83588D-13, 1.19330D-12, 1.41705D-12,
     $     1.91911D-12, 2.19885D-12, 4.10344D-12, 4.77542D-12,
     $     8.50953D-12, 8.84453D-12, 1.20157D-11, 1.50905D-11,
     $     2.01095D-11, 3.22709D-11, 3.51595D-11, 4.64215D-11,
     $     5.43908D-11, 5.93826D-11, 9.03299D-11, 1.04636D-10,
     $     1.22550D-10, 1.45893D-10, 1.68347D-10, 1.74333D-10,
     $     2.21043D-10, 2.37084D-10, 2.83811D-10, 3.68297D-10,
     $     4.10914D-10, 6.03069D-10, 7.20602D-10, 7.55839D-10,
     $     9.56928D-10, 1.07906D-09, 1.43673D-09, 1.49682D-09,
     $     1.74266D-09, 1.94829D-09, 2.08480D-09, 2.27735D-09,
     $     2.61419D-09, 2.74763D-09, 3.20694D-09, 3.58409D-09,
     $     4.35668D-09, 6.04749D-09, 6.59872D-09, 7.05210D-09,
     $     8.87933D-09, 9.35219D-09, 1.03579D-08, 1.07568D-08,
     $     1.22495D-08, 1.30225D-08, 1.62045D-08, 2.00317D-08,
     $     2.21918D-08, 2.67901D-08, 2.81476D-08, 2.92743D-08,
     $     3.20748D-08, 3.54197D-08, 4.17416D-08, 4.34659D-08,
     $     4.69318D-08, 5.23149D-08, 5.42711D-08, 6.31420D-08 /

      do j = 1, lines
         if ( u >= exen(j) .and. u < exen(j+1) ) then

            do i = 1, nimax
               if ( ifz(i)*1000+ ifa(i) == 1. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(1,j), exen(j+1),
     $                    br(1,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 1001. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(2,j), exen(j+1),
     $                    br(2,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 8014. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(4,j), exen(j+1),
     $                    br(4,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 1002. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(5,j), exen(j+1),
     $                    br(5,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 6014. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(6,j), exen(j+1),
     $                    br(6,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 1003. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(7,j), exen(j+1),
     $                    br(7,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 2003. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(8,j), exen(j+1),
     $                    br(8,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 2004. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(9,j), exen(j+1),
     $                    br(9,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 5012. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(10,j), exen(j+1),
     $                    br(10,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 4012. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(11,j), exen(j+1),
     $                    br(11,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 6011. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(12,j), exen(j+1),
     $                    br(12,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 5011. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(13,j), exen(j+1),
     $                    br(13,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 4011. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(14,j), exen(j+1),
     $                    br(14,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 2006. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(15,j), exen(j+1),
     $                    br(15,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 3006. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(16,j), exen(j+1),
     $                    br(16,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 4010. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(17,j), exen(j+1),
     $                    br(17,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 3007. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(18,j), exen(j+1),
     $                    br(18,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 4007. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(19,j), exen(j+1),
     $                    br(19,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 3009. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(20,j), exen(j+1),
     $                    br(20,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 5008. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(22,j), exen(j+1),
     $                    br(22,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 3007. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(23,j), exen(j+1),
     $                    br(23,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 5010. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(24,j), exen(j+1),
     $                    br(24,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 6010. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(25,j), exen(j+1),
     $                    br(25,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 6011. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(26,j), exen(j+1),
     $                    br(26,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 6012. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(27,j), exen(j+1),
     $                    br(27,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 6013. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(28,j), exen(j+1),
     $                    br(28,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 7013. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(29,j), exen(j+1),
     $                    br(29,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               else if ( ifz(i)*1000+ ifa(i) == 7014. ) then
                  if ( r(i) /= 0. ) then
                     brr = linearone( u, exen(j), br(30,j), exen(j+1),
     $                    br(30,j+1))
                     integral = integral+ brr
                     r(i) =  sigma* brr
                  end if

               end if

            end do

            do i = 1, nimax
               r(i) = r(i)/ integral
            end do

         end if

      end do

      end subroutine suppressalpha_o16


!     One-dimensional linear interpolation for (x,y) inbetween (x1,y1) & (x2,y2)
      double precision function linearone(x,x1,y1,x2,y2)

      implicit none
      double precision x1,y1,x2,y2
      double precision x,y

      y = ( y2- y1)/( x2- x1)*( x- x1)+ y1
      linearone = y

      return
      end function linearone


! ----------------------------------------------------------------------
!     Photonuclear pion production
      subroutine pionprod
!     02/18/2014 by S. Noda
! ----------------------------------------------------------------------

      implicit real*8(a-h,o-z)

      include 'param00.inc'

*-----------------------------------------------------------------------

      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr
!$OMP THREADPRIVATE(/const0/)
      common /coodrp/ r(5,nnn), p(6,nnn)
!$OMP THREADPRIVATE(/coodrp/)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)
!$OMP THREADPRIVATE(/coodid/)
      common /pnint/  ipnint
      common /gg005/ xxx, yyy, zzz, uuu, vvv, www, tme, erg,
     &               dls, wgt, vel, dtc,
     &               icl, iii, jjj, kkk, jsu, iap, jgp, ipt,
     &               mtp, iexp, iex, idx,
     &               npa, ncp
!$OMP THREADPRIVATE(/gg005/)

*-----------------------------------------------------------------------
cABE add @2014/06/13
      common /qmdflg/ irqmd
cABE 2016/03/09
      common /muflag/ imuinthit,imubrmhit,imuppdhit,imucaphit,imucapflag
!$OMP THREADPRIVATE(/muflag/)
      common /prmui/ prmui1
cABE 2016/07/13
      common /muflag2/ betaeve,ishadowflag
!$OMP THREADPRIVATE(/muflag2/)

*-----------------------------------------------------------------------
cABE end

      dimension pc(5,2)  ! p in center-of-mass system
      dimension pp(5)    ! momentum and energy of incident photon
      dimension beta(3)
      double precision gamma      ! beta & gamma for Lorentz transform
      integer ideltaid           ! nucleon id to absorb photon
      dimension ecm(3)   ! energy of CM system
      double precision pcm      ! momentum in CM
      double precision theta, phi ! for Lorentz transform
      integer k, l, k2          ! counters
      double precision rand     ! random number
      double precision unirn    ! defined function

      double precision sigma_delta, sigma_nstar ! true pion-production
                                                ! cross section
      double precision sigma_temppion ! temporary pion-production cross section
      double precision spp, stp       ! receivers for above

cABE add @2014/06/13, to calculate (mass^2 + 2*mass*epot)
      double precision epotdl         ! potential energy for nucleon
cABE end

      k = 0
      l = 0
      k2 = 0
      ieo = 6

!-----------------------------------------------------------------------

      e = erg/ 1e3
      pp(1) = uuu* e
      pp(2) = vvv* e
      pp(3) = www* e
      pp(4) = e
      pp(5) = 0.

!     (1) Choose one nucleon at random
      do k2 = 1, 30
         do k = 1, massta
            ideltaid = int( unirn()* massta)+ 1
c     ideltaid = int( unirn()* mstapr)+ 1   ! pick proton
c     ideltaid = int( unirn()*( massta- mstapr))+ mstapr+ 1   ! pick neutron

!     ecm(1): energy of center-of-mass system
cABE change @2014/06/13, p(4,*) is changed energy -> energy(pot included)
            if(irqmd .eq. 0) then
               ecm(1) = sqrt((p(4,ideltaid) + pp(4))**2
     &                       - (p(1,ideltaid) + pp(1))**2
     &                       - (p(2,ideltaid) + pp(2))**2
     &                       - (p(3,ideltaid) + pp(3))**2)
            else
               ecm(1) = sqrt((p(6,ideltaid) + pp(4))**2
     &                       - (p(1,ideltaid) + pp(1))**2
     &                       - (p(2,ideltaid) + pp(2))**2
     &                       - (p(3,ideltaid) + pp(3))**2)
            endif
cABE end

            if ( ecm(1) > 0.938+ 0.135 ) exit
         end do
      end do

      if ( k > massta) then
         write(ieo,*)
     $        "Error: cannot produce pion in photonuclear reaction"
         call parastop( 119)
      end if

cABE add @2014/06/26
      if(irqmd .ne. 0) then
         epotdl = (p(4,ideltaid)**2 - p(5,ideltaid)**2
     &             - p(1,ideltaid)**2 - p(2,ideltaid)**2
     &             - p(3,ideltaid)**2) * 0.5d0 / p(5,ideltaid)
      endif
cABE end

!     (2) Laboratory system (Lab) ==> Center of Mass system (CM)

!     pc(x,1): nuclecon in CM
!     beta(1): parameter beta for x-axis
!     beta(2): parameter beta for y-axis
!     beta(3): parameter beta for z-axis
!     gamma: parameter gamma

cABE change @2014/06/13, p(4,*) is changed energy -> energy(pot included)
      if(irqmd .eq. 0) then
         beta(1) = ( p(1,ideltaid)+ pp(1))/( p(4,ideltaid)+ pp(4))
         beta(2) = ( p(2,ideltaid)+ pp(2))/( p(4,ideltaid)+ pp(4))
         beta(3) = ( p(3,ideltaid)+ pp(3))/( p(4,ideltaid)+ pp(4))
      else
         beta(1) = ( p(1,ideltaid)+ pp(1))/( p(6,ideltaid)+ pp(4))
         beta(2) = ( p(2,ideltaid)+ pp(2))/( p(6,ideltaid)+ pp(4))
         beta(3) = ( p(3,ideltaid)+ pp(3))/( p(6,ideltaid)+ pp(4))
      endif
cABE end

      gamma = 1./ sqrt( 1.- beta(1)* beta(1)- beta(2)* beta(2)- beta(3)*
     $     beta(3))


!     (3) in CM, one nucleon absorbs photon

!     ecm(2): energy of nucleon in CM
!     Nucleon turns into delta
!     Rest mass of delta particles in GeV/c^2
      p(5,ideltaid) = ecm(1)

!     inds: 1 nucleon (proton, neutron)
!           2 delta
!           3 N*
!           4 pion
cABE 2016/07/13
ccABE 2016/03/09
cc      spp = sigma_delta( dble(massta), 0.16, dble(e*1e3),1,1)
cc      stp = sigma_delta( dble(massta), 0.16, dble(e*1e3),2,4)
c      if(imuinthit .ne. 1 ) then
c       spp = sigma_delta( dble(massta), dble(e*1e3),1,1)
c       stp = sigma_delta( dble(massta), dble(e*1e3),2,4)
c      else
c       spp = prmui1*sigma_delta( 1.d0, dble(e*1e3),1,1)
c     &      + (1.d0-prmui1)*sigma_delta( dble(massta), dble(e*1e3),1,1)
c       stp = prmui1*sigma_delta( 1.d0, dble(e*1e3),2,4)
c     &      + (1.d0-prmui1)*sigma_delta( dble(massta), dble(e*1e3),2,4)
c      endif
       spp = sigma_delta( dble(massta), dble(e*1e3),1,1)
       stp = sigma_delta( dble(massta), dble(e*1e3),2,4)

      if ( unirn() < spp /(spp+ stp)) then
         inds(ideltaid) = 2

      else
         inds(ideltaid) = 3

      endif

!     ecm(2) = 1./ 2.*( ecm(1)+( p(5,ideltaid)**2)/ ecm(1))
      ecm(2) = ecm(1)

!     absolute value of nucleon's momentum in CM
      pcm = 0.


!     (4) CM ==> Lab
      pc(1,1) = 0.
      pc(2,1) = 0.
      pc(3,1) = 0.
      pc(4,1) = ecm(2)

!     Lorentz transform
      p(1,ideltaid) = beta(1)* gamma* pc(4,1)
      p(2,ideltaid) = beta(2)* gamma* pc(4,1)
      p(3,ideltaid) = beta(3)* gamma* pc(4,1)
cABE change @2014/06/13, p(4,*) is changed energy -> energy(pot included)
      if(irqmd .eq. 0) then
         p(4,ideltaid) = sqrt(p(1,ideltaid)**2 + p(2,ideltaid)**2
     &                        + p(3,ideltaid)**2 + p(5,ideltaid)**2)
      else
         p(4,ideltaid) = sqrt(p(1,ideltaid)**2 + p(2,ideltaid)**2
     &                        + p(3,ideltaid)**2 + p(5,ideltaid)**2
     &                        + 2.0d0 * p(5,ideltaid) * epotdl)
         p(6,ideltaid) = sqrt(p(1,ideltaid)**2 + p(2,ideltaid)**2
     &                        + p(3,ideltaid)**2 + p(5,ideltaid)**2)
      endif
cABE end

cABE add @2015/04/10
      iavd(ideltaid) = -1

      return
      end subroutine pionprod


! ---
! Port from PICA98
! ---

      function Vcm(Eg,clight,AMp,AMpi)          !Obtain Speed of CMS
      implicit real*8(a-h,o-z)
      Vcm=clight*Eg/(Eg+AMp*clight**2)
      end

      function p1(Eg,pl,v,clight,AMp,AMpi)      !Obtain Photon Moment in CMS
      implicit real*8(a-h,o-z)
      p1=(Eg-pl*v)/(clight*sqrt(1-v**2/clight**2))
      end

      function Ep1(v,clight,AMp,AMpi)           !Obtain Initial Proton Energy in CMS
      implicit real*8(a-h,o-z)
      Ep1=AMp*clight**2/sqrt(1-v**2/clight**2)
      end

      function q1(w,clight,AMp,AMpi)            !Obtain Pion Moment
      implicit real*8(a-h,o-z)
      big=1.001
      small=0.999
!     Non-relativistic Value
      q1=sqrt(2*AMp*AMpi/(AMp+AMpi)* (w-AMpi*clight**2-AMp*clight **2))
      value1=subt(q1,w,clight,AMp,AMpi)
      value2=subt(q1*big,w,clight,AMp,AMpi)
      value3=subt(q1*small,w,clight,AMp,AMpi)
      if(value1.lt.value2.and.value1.lt.value3) return
      if(value2.lt.value3) then
	 times=big
	 value1=value2
      else
	 times=small
	 value1=value3
      endif
      q1=q1*times
      do 10 i=1,5000
	 value2=subt(q1*times,w,clight,AMp,AMpi)
	 if(value2.gt.value1) return
	 value1=value2
	 q1=q1*times
 10   continue

      write(6,*) 'Error in function q1'
      stop
      end

      function subt(q,w,clight,AMp,AMpi)        !Function used in q1(w)
      implicit real*8(a-h,o-z)
      subt=abs(sqrt(AMp**2*clight**4+q**2*clight**2) +sqrt(AMpi**2
     $     *clight**4+q**2*clight**2)-w)
      end

c     **** Subroutine Fshadow is used for calculating Shadowing
c     Parameter **
      subroutine Fshadow(amas)
      implicit real*8(a-h,o-z)
      parameter (Eglow=0.0,Ewid=20.0) !The Lowest Photon Energy, and a
                                      !Step of The Energy
      parameter (nebin=299)
      parameter (nvm=5)         !#Vector Meson
      parameter (maxshad=10000)
      parameter (nal=100)       !#Steps in Alpha Integration
      common /shadow1/amv(nvm),gv(nvm),sigvn(nvm) !Vector Meson
                                                  !Properties
      common /shadow2/ccmps,hbar,pi,ALEM,pii,ck,rc
      dimension sigvmda(nebin)  ! VMD term for nuclear reaction
      dimension sigga(nebin)    ! Sigma for photo-nuclear reaction
      common /shadow3/shadow(nebin) ! Shadowing Effect
      astart=0.5/(nal*1.0)
      astep=1.0/(nal*1.0)
c     **** Set Constant Values ********
      ccmps=2.9979e10           ! Light speed (cm/s)
      hbar=1.0545887e-34/1.602e-13 !Hbar (MeV*s)
      pi=3.141592
      ALEM=1.0/137.0
      pii=1.0/12.0/pi**2*3.0*(1+4+1+4)/9.0
      CK=1.7
      rc=1.3e-13                ! Rc (cm)
      ra=0.93e-13               ! parameter a, for nuclear radius
                                ! sqrt(3/2)*ra*A**(1/3)
C      data amv/769.9,781.9,1019.41,3096.88,3686.0/   !NS 2021.04 comment out for NVIDIA TOOL KIT COMPILL ERR
C      data gv/2.0,23.1,13.2,10.5,30.6/               !NS 2021.04 comment out for NVIDIA TOOL KIT COMPILL ERR
C      data sigvn/22.0,25.0,9.0,2.2,1.3/              !NS 2021.04 comment out for NVIDIA TOOL KIT COMPILL ERR
      do 10 i=1,nvm
         gv(i)=sqrt(4*pi*gv(i))
	 sigvn(i)=sigvn(i)*1.0e-27 !Convert (mb) to (cm**2) if(i.ne.1)
c     sigvn(i)=1/amv(i)**2*(hbar*2*pi*ccmps)**2
 10   continue

c     **** For Photon-Nucleon Reaction ********* **** Calculate First
c     term in Eq(46) ******
      do 20 i=1,nvm
 20      sigvmd=sigvmd+4.0*pi*alem/(gv(i)**2)*sigvn(i)
c     **** Calculate Second term in Eq(46) ******
         u2=amv(3)**2           !Lower Limit of Integration Range
         v2=0.0
         do 40 i=1,maxshad
            du2=u2*0.01
            tmp=0.0
            do 50 ia=1,nal
               alpha=astart+astep*(ia-1)
               sighn=sighnf(sqrt(u2),alpha)
               tmp=tmp+du2*pii/u2*sighn*astep
 50         continue
            v2=v2+tmp
            if(tmp.le.0.0005*v2) goto 60
            u2=u2*1.01
 40      continue
         write(6,*) 'maxshad is too Small'
         stop
 60      siggn=sigvmd+4*pi*alem*ck*v2
c     **** For Photo-Nuclear Reaction *********
         do 65 ie=1,nebin
            eg=eglow+ewid*(ie-1)
            if(eg.le.1000) then
               shadow(ie)=1.0
               goto 65
            endif
c     **** Calculate VMD term *****************
            do 67 i=1,nvm
               u2=amv(i)**2
               ram=2*eg/u2*hbar*ccmps
               tmp1=(1-amas**0.3333*sigvn(i)/(8.0*pi*ra**2)*(amas-1)
     $              /amas*exp(-ra**2*amas**0.6667/(4.0*ram**2)))
               tmp2=amas*sigvn(i)*tmp1
               sigvmda(ie)=sigvmda(ie)+4.0*pi*alem/(gv(i)**2)*tmp2
 67         continue
c     **** Calculate Continuum (second) term ***
            u2=amv(3)**2        !Lower Limit of Integration Range
            v2=0.0
            do 70 i=1,maxshad
               du2=u2*0.01
               ram=2*eg/u2*hbar*ccmps !Lambda
               tmp=0.0
               do 80 ia=1,nal
                  alpha=astart+astep*(ia-1)
                  sighn=sighnf(sqrt(u2),alpha)
                  tmp1=(1-amas**0.3333*sighn/(8.0*pi*ra**2)*(amas-1)
     $                 /amas*exp(-ra**2*amas**0.6667/(4.0*ram**2)))
                  sigha=amas*sighn*tmp1
                  tmp=tmp+du2*pii/u2*sigha*astep
 80            continue
               v2=v2+tmp
               if(tmp.le.0.001*v2) goto 90
               u2=u2*1.01
 70         continue
            write(6,*) 'Maxshad is too Small'
            stop
 90         sigga(ie)=sigvmda(ie)+4*pi*alem*ck*v2
            shadow(ie)=sigga(ie)/(amas*siggn)
 65      continue
         return
         end

c     **** Fucntion SighNf(u**2,alpha) used in Fshadow *********
      function sighnf(u,alpha)
      implicit real*8(a-h,o-z)
      parameter (nvm=5)         !#Vector Meson
      common /shadow1/amv(nvm),gv(nvm),sigvn(nvm) !Vector Meson
                                                  !Properties
      common /shadow2/ccmps,hbar,pi,ALEM,pii,ck,rc
      tmp=1.0/alpha/(1.0-alpha)*4.0/u**2*(hbar*ccmps)**2
      if(tmp.le.rc**2) then
	 sighnf=tmp
      else
	 sighnf=rc**2
      endif
      return
      end

! ---
! END (Port from PICA98)
! ---


!     Photon-incident pion-production cross seciton of delta resonace
!     Port from subroutine resofit in Pic98.for (PICA98)
!     S.Noda 12/13/2013
!
!     amas(input): mass number
!     denst(input): nuclear density in fm^{-3}
!     Eg(input): incoming photon energy (MeV)
!     Ouputs cross section in barn (b)

      double precision function sigma_delta(amas,Eg,nrstart,nrend)
c     subroutine resofit(amas,denst,fittot)

      implicit real*8(a-h,o-z)
      parameter (nreso=4)       !#Resonance
      parameter (nebin=299)     !#Incident Photon Energy
      parameter (Eglow=0.0,Ewid=20.0) !The Lowest Photon Energy, and a
                                      !Step of The Energy
      parameter (Thres=0.1503)  !Threshold of Pi Production in GeV
      parameter (a1=0.091,a2=0.071) !Background parameter
                                    !(mb),(mb*GeV**0.5)
cABE commentout @2014/06/25
cABE      common /shadow3/shadow(nebin) ! Shadowing Effect
cABE end
      dimension Wr(nreso),      !Wr Resonance Energy
     1     Tr(nreso),           !Tr Resonance Width
     2     sigrp(nreso),        !Ir Resonance Intensity from proton
     2     sigrd(nreso),        !Ir Resonance Intensity from deuteron
     2     sigra(nreso),        !Ir Resonance Intensity for average nucleon
     3     l(nreso),jg(nreso),  !l and Jg Angular Moment
     4     qr(nreso),           !Pion Moment at the Resonance Peak
                                !Energy
     5     pr(nreso),           !Photon Moment at the Resonance Peak
                                !Energy
     6     X(nreso),            !Parameter X
     7     Bp(nreso)            !Parameter Bp (for Width Boardening by
c     !P.B.)  dimension shadowC(60),shadowBx(60) !Parameter for
c     Shadowing Effect

      data clight/2.99792e8/     !Light Speed
      data AMp/1.044002737e-14/  !Mass of Nucleon in MeV/c^2 (Proton), 938.3/clight**2
      data AMpi/1.55289705e-15/  !Mass of Pion in MeV/c^2 (Pi+), 139.567/clight**2

c     ***** Set Constant Values *******
c      clight=2.99792e8          !Light Speed
c      AMp=938.3/clight**2       !Mass of Nucleon in MeV/c^2 (Proton)
c      AMpi=139.567/clight**2    !Mass of Pion in MeV/c^2 (Pi+) ***** Set
c     Shadowing Parameter ********* data ShadowC/ 1
c     1.000,1.000,1.000,1.000,1.000, 2 0.990,0.980,0.975,0.970,0.965, 3
c     0.960,0.957,0.955,0.952,0.949, 4 0.946,0.943,0.940,0.937,0.934, 5
c     0.930,0.927,0.926,0.923,0.920, 6 0.916,0.913,0.910,0.906,0.903, 7
c     0.900,0.897,0.893,0.890,0.887, 8 0.885,0.883,0.881,0.879,0.877, 9
c     0.875,0.874,0.872,0.870,0.869, 1 0.868,0.866,0.864,0.863,0.861, 1
c     0.860,0.859,0.858,0.857,0.856, 2 0.855,0.854,0.853,0.852,0.851/
c     data ShadowBx/ 1 0.000, 0.000, 0.000, 0.000, 0.000, 2
c     -0.005,-0.007,-0.010,-0.010,-0.011, 3
c     -0.011,-0.012,-0.014,-0.016,-0.018, 4
c     -0.021,-0.024,-0.026,-0.029,-0.032, 5
c     -0.035,-0.037,-0.039,-0.041,-0.043, 6
c     -0.045,-0.047,-0.049,-0.051,-0.053, 7
c     -0.054,-0.056,-0.058,-0.060,-0.062, 8
c     -0.064,-0.067,-0.070,-0.072,-0.075, 9
c     -0.078,-0.081,-0.083,-0.085,-0.088, 1
c     -0.090,-0.093,-0.095,-0.098,-0.100, 1
c     -0.103,-0.106,-0.109,-0.112,-0.115, 2
c     -0.118,-0.121,-0.124,-0.127,-0.130/ **** Set Resonance Parameter
c      data Wr/1230.0,1505.0,1671.0,1423.0/
      data Tr/122.0,100.0,100.0,66.0/
cABE 2016/03/09, Ir/A for deuteron and average nucleon
      data sigrp/0.415,0.133,0.065,0.047/
      data sigrd/0.3195d0,0.0975d0,0.0029d0,0.0029d0/
      data sigra/0.3195d0,0.045d0,0.002d0,0.002d0/

      data l/1,2,3,1/
cABE 2016/03/09
c      data Jg/1,1,2,1/
      real*8 Jg
      data Jg/1.5d0,1.5d0,2.5d0,0.5d0/
c      data X/150.0,350.0,350.0,350.0/
      data Bp/0.77,0.9,0.9,0.9/

cABE add @2014/07/31
      dimension ishflg(300)
      dimension shadat(300,1000)

      data ishflg/300*0/
      data shadat/300000*0.0/

      save ishflg, shadat
cABE end

cABE 2015/08/10
      common /muflag/ imuinthit,imubrmhit,imuppdhit,imucaphit,imucapflag
!$OMP THREADPRIVATE(/muflag/)
cABE 2016/07/13
      common /muflag2/ betaeve,ishadowflag
!$OMP THREADPRIVATE(/muflag2/)

      sigma_delta = 0.
      if ( Eg < 140. ) return
cABE 2016/03/09
      if( amas .le. 2.d0 ) then
       denst = 0.09d0
      else
       denst = 0.16d0
      endif

cABE 2016/03/09, nuclear medium effects
c      Wr(1)=1230.0+266.67*denst+9.0
      if( amas .le. 1.d0 ) then
       deltM = 0.d0
      else
       deltM = 267.d0 * denst
      endif
      Wr(1)=1230.0+deltM
      Wr(2)=1505.0
      Wr(3)=1671.0
      wr(4)=1423.0
      X(1)=150.0/clight
      X(2)=350.0/clight
      X(3)=350.0/clight
      X(4)=350.0/clight

c      do 5 nr=1,nreso
c 5    X(nr)=X(nr)/clight     !InputData for X is written in MeV,so
c     !convert to MeV/c **** Considering Nuclear Medium Effect ********
c      Wr(1)=Wr(1)+266.67*denst+9.0
c     **** Determine pr and qr *********
         do 10 nr=1,nreso
            tempmin=1000
c            do 20 ie=1,1000
c               Eg=ie*1.0        !Photon Energy in Lab.
               pl=Eg/clight     !Photon Moment in Lab.
               v=Vcm(Eg,clight,AMp,AMpi)        !Speed of C.M.S.
               P=p1(Eg,pl,v,clight,AMp,AMpi)    !Photon Moment in CMS
               Ep=Ep1(v,clight,AMp,AMpi)        !Proton Energy in CMS
               W=p*clight+Ep    !Total Energy in CMS
               temp=abs(Wr(nr)-W)

cABE change @2014/06/26, avoid that Pr(nr) has no parameter if temp .gt. tempmin
               Pr(nr)=P      !Photon Moment at the Resonance in CMS
cABE end

c 20         continue
            qr(nr)=q1(Wr(nr),clight,AMp,AMpi)   !Pion Moment at the Resonance in CMS
 10      continue

c     ***** Calculate Cross Section *************
c         do 30 ie=1,nebin
c            Eg=Eglow+Ewid*(ie-1) !Photon Energy in Lab.
            pl=Eg/clight        !Photon Moment in Lab.
            v=Vcm(Eg,clight,AMp,AMpi)           !Speed of C.M.S.
            P=p1(Eg,pl,v,clight,AMp,AMpi)       !Photon Moment in CMS
            Ep=Ep1(v,clight,AMp,AMpi)           !Proton Energy in CMS
            W=p*clight+Ep       !Total Energy in CMS
            if(w-AMpi*clight**2-AMp*clight**2.lt.0) goto 30 !Threshold
                                                            !Check
            q=q1(w,clight,AMp,AMpi)             !Pion Moment
            do 40 nr=nrstart,nrend
cABE add @2014/06/25, for non-resonant background
             if (nr .le. 4) then
cABE end
cABE 2016/03/09, nuclear medium effects
c               if(nr.eq.1) then
c                  tmpr=Tr(1)*Bp(1)+486.67*denst+17.0
c               else
c                  tmpr=Tr(nr)*Bp(nr)+2033.33*denst-20.0
c               endif
               if( amas .le. 1.d0 ) then
                deltT=0.d0
               else
                if( nr .eq. 1 ) then
                 deltT=487.d0*denst+17.d0
                else
                 deltT=2033.d0*denst-20.d0
                endif
               endif
               tmpr=Tr(nr)*Bp(nr)+deltT

               T=tmpr*(q/qr(nr))**(2*l(nr)+1) *((qr(nr)**2+X(nr)**2)/(q
     $              **2+X(nr)**2))**l(nr)
cABE 2016/03/09
c               Tg=Tr(nr)*(p/pr(nr))**(2*Jg(nr)) *((pr(nr)**2+X(nr)**2)
               Tg=tmpr*(p/pr(nr))**(2*Jg(nr)) *((pr(nr)**2+X(nr)**2)
     $              /(p**2+X(nr)**2))**Jg(nr)
c               tmp1=sigr(nr)*(pr(nr)*qr(nr)/p/q)*Wr(nr)**2*T*Tg !Component
               tmp1=(pr(nr)*qr(nr)/p/q)*Wr(nr)**2*T*Tg !Component
                                                       ! (Bunsi)
               tmp2=(W**2-Wr(nr)**2)**2+Wr(nr)**2*T**2 !Dominantor
                                                       !(Bunbo)

cABE 2016/03/09
c               sigma_delta=sigma_delta+tmp1/tmp2
               if( amas .le. 1.d0 ) then
                sigma_delta=sigma_delta+sigrp(nr)*tmp1/tmp2
               elseif( amas .le. 2.d0 ) then
                sigma_delta=sigma_delta+sigrd(nr)*tmp1/tmp2
               else
                sigma_delta=sigma_delta+sigra(nr)*tmp1/tmp2
               endif

cABE add @2014/06/25, for non-resonant background
cABE change @2014/07/31, to improve cal. time
             elseif (nr .eq. 5) then

               EgGeV = Eg * 1.0d-3

               if (EgGeV .gt. thres) then
c...N.Bianchi,1996 eq(6)
                  bak = (a1 + a2 * EgGeV**(-0.5))
     &                  * (1 - exp(-2.0d0*(EgGeV-thres)))
cABE 2016/07/13
ccABE 2016/03/09, switch to L.B.Bezrukov,proc.ICRC,vol.10,p.245,1979
c                  if( EgGeV .gt. 4.48d0 ) then
c                   bak = 114.3d0 + 1.647d0
c     &                   * dlog(2.d0*0.938272046d0*EgGeV/88.d0)**2.d0
c                   bak = bak * 1.d-3
c                  endif
c...Original function for muon photonuclear interaction
                  if (EgGeV .gt. 3.d0) then
                   dum1 = 0.059679d0
                   dum2 = 0.086457d0
                   dum3 = -0.15393d0
                   bak = (dum1 + dum2 * EgGeV**dum3)
                  endif


                  eshmin = 1.0d+3
                  eshmax = 1.0d+6
                  nsh = 30
                  ddesh = log(eshmax/eshmin) / dble(nsh)

                  ia = idnint(amas)

!$OMP CRITICAL (shadow_crit)
                  if (ishflg(ia) .ne. 1) then
                     do ish = 0, nsh
                        esh = eshmin * exp((dble(ish))*ddesh)
cABE 2016/03/09
c                        call FshadowEg(amas,esh,shadowt)
                        if( ia .eq. 1 ) then
                         shadowt = 1.d0
                        else
                         call FshadowEg(amas,esh,shadowt)
                        endif
                        shadat(ia,ish+1) = shadowt    ! S.Abe 2015/08/18, changed
                     enddo
                     ishflg(ia) = 1
                  endif
!$OMP END CRITICAL (shadow_crit)

cABE 2016/07/13
                  if (Eg .le. eshmin .or.
     &                (imuinthit.eq.1 .and. ishadowflag.eq.-1)) then

                     shadow = 1.0d0

                  elseif (Eg .ge. eshmax) then

                     eshlo = eshmin * exp((dble(nsh-1))*ddesh)
                     eshhi = eshmin * exp((dble(nsh))*ddesh)

                     shadow = (shadat(ia,nsh) - shadat(ia,nsh-1))
     &                       / (eshhi - eshlo) * (Eg - eshhi)
     &                       + shadat(ia,nsh)

                  else

                     ilo = 0
                     ihi = nsh

  100                imid = (ilo + ihi) / 2

                     esh1 = eshmin * exp((dble(imid))*ddesh)
                     esh2 = eshmin * exp((dble(imid+1))*ddesh)

                     if (Eg .ge. esh1 .and. Eg .lt. esh2) then
                        imid = imid + 1      ! S.Abe 2015/08/14, added
                        goto 110
                     elseif (esh1 .lt. Eg) then
                        ilo = imid
                        goto 100
                     else
                        ihi = imid
                        goto 100
                     endif

  110                continue

                     shadow = (shadat(ia,imid+1) - shadat(ia,imid))
     &                       / (esh2 - esh1) * (Eg - esh1)
     &                       + shadat(ia,imid)

                  endif

                  sigma_delta = sigma_delta + bak * shadow

               endif

             endif
cABE end

 40         continue
 30      continue
ccc 30         call Fshadow(amas)  ! Currently Fshadow does not work!!
c  **** Add Non-Resonant Background $ Subtract QDD contribution *******
c         do 50 ie=1,nebin
c$$$            Eg=Eglow+Ewid*(ie-1) !Photon Energy in Lab.
c$$$            EgGeV=Eg*0.001
c$$$            if(EgGeV.gt.thres) then !Add Background
c$$$               bak=(a1+a2*EgGeV**(-0.5))*(1-exp(-2*(EgGeV-thres)))
c$$$c     **** Considering Shadowing Effect only for background term *****
c$$$c     ie1=(ie-1)/5+1 ShadowAx=1.0-ShadowBx(ie1)*log(12.0)
c$$$c     Shadow=ShadowC(ie1)*(ShadowAx+ShadowBx(ie1)*log(amas))
c$$$               sigma_delta=sigma_delta+bak*shadow(ie)
c$$$            endif

c 50      continue

C     Convert from NN cross section (mb) to total cross section (b)
         sigma_delta = sigma_delta/ 1e3* amas

         return
         end


cABE add @2014/06/25
************************************************************************
*                                                                      *
      subroutine FshadowEg(amas,Eg,shadow)
*                                                                      *
*     Subroutine Fshadow is used for calculating Shadowing Parameter   *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*        amas   : mass number                                          *
*        eg     : incoming photon energy (MeV)                         *
*                                                                      *
*     output :                                                         *
*                                                                      *
*        shadow : shadowing parameter                                  *
*                                                                      *
*---- memo ------------------------------------------------------------*
*                                                                      *
*     common /shadow1/ and /shadow2/ is shared with function "sighnf"  *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

      parameter (nvm = 5)                                               ! #Vector Meson
      parameter (maxshad = 10000)
      parameter (nal = 100)                                             ! #Steps in Alpha Integration

      common /shadow1/amv(nvm),gv(nvm),sigvn(nvm)                       ! Vector Meson Properties

      data amv/769.9,781.9,1019.41,3096.88,3686.0/
      data gv/2.0,23.1,13.2,10.5,30.6/
      data sigvn/22.0,25.0,9.0,2.2,1.3/

      common /shadow2/ccmps,hbar,pi,ALEM,pii,ck,rc

      data initialflag /0/
      save initialflag

*-----------------------------------------------------------------------

      shadow = 1.0d0
      if (eg .le. 1.d+3) return

*-----------------------------------------------------------------------
*     Set Constant Values
*-----------------------------------------------------------------------

      astart = 0.5d0 / (nal * 1.0d0)
      astep = 1.0d0 / (nal * 1.0d0)

      ccmps = 2.9979d+10                                                ! Light speed (cm/s)
      hbar = 1.0545887d-34 / 1.602d-13                                  ! Hbar (MeV*s)
      pi = 3.141592d0
      ALEM = 1.0d0 / 137.0d0
      pii = 1.0d0 / 12.0d0 / pi**2 * 3.0d0 * (1d0+4d0+1d0+4d0) / 9.0d0
      CK = 1.7d0
      rc = 1.3d-13                                                      ! Rc (cm)

      ra = 0.93d-13                                                     ! parameter a, for nuclear radius = sqrt(3/2)*ra*A**(1/3)

cABE correct @2015/02/09
      if( initialflag .eq. 0) then
       initialflag = 1
       do 10 i = 1, nvm
        gv(i) = dsqrt(4.0d0 * pi * gv(i))
        sigvn(i) = sigvn(i) * 1.0d-27                                   ! Convert (mb) to (cm**2) if(i.ne.1)
c        sigvn(i)=1/amv(i)**2*(hbar*2*pi*ccmps)**2
   10  continue
      endif

*-----------------------------------------------------------------------
*     For Photon-Nucleon Reaction
*-----------------------------------------------------------------------
*     Calculate First term in Eq(46)
*-----------------------------------------------------------------------

      sigvmd = 0.0d0

      do 20 i = 1, nvm

         sigvmd = sigvmd + 4.0d0 * pi * alem / (gv(i)**2) * sigvn(i)

   20 continue

*-----------------------------------------------------------------------
*     Calculate Second term in Eq(46)
*-----------------------------------------------------------------------

      u2 = amv(3)**2                                                    ! Lower Limit of Integration Range
      v2 = 0.0d0

      do 40 i = 1, maxshad

         du2 = u2 * 0.01d0
         tmp = 0.0d0

         do 50 ia = 1, nal
            alpha = astart + astep * dble(ia-1)
            sighn = sighnf(sqrt(u2),alpha)
            tmp = tmp + du2 * pii / u2 * sighn * astep
   50    continue

         v2 = v2 + tmp

         if (tmp .le. 0.0005d0*v2) goto 60

         u2 = u2 * 1.01d0

   40 continue

*-----------------------------------------------------------------------

      write(6,*) 'maxshad is too Small'
      stop

*-----------------------------------------------------------------------

   60 siggn = sigvmd + 4.0d0 * pi * alem * ck * v2

*-----------------------------------------------------------------------
*     For Photo-Nuclear Reaction
*-----------------------------------------------------------------------
*     Calculate VMD term
*-----------------------------------------------------------------------

      sigvmda = 0.0d0

      do 67 i = 1, nvm

         u2 = amv(i)**2

         ram = 2.0d0 * eg / u2 * hbar * ccmps
         tmp1 = (1 - amas**0.3333 * sigvn(i) / (8.0d0 * pi * ra**2)
     &           * (amas-1.0d0) / amas
     &           * exp(-ra**2 * amas**0.6667 / (4.0d0 * ram**2)))

         tmp2 = amas * sigvn(i) * tmp1

         sigvmda = sigvmda + 4.0d0 * pi * alem / (gv(i)**2) * tmp2

   67 continue

*-----------------------------------------------------------------------
*     Calculate Continuum (second) term
*-----------------------------------------------------------------------

      u2 = amv(3)**2                                                 ! Lower Limit of Integration Range
      v2 = 0.0d0

      do 70 i = 1, maxshad

         du2 = u2 * 0.01d0
         ram = 2.0d0 * eg / u2 * hbar * ccmps                            ! Lambda
         tmp = 0.0d0

         do 80 ia = 1, nal

            alpha = astart + astep * (ia-1)
            sighn = sighnf(sqrt(u2),alpha)

            tmp1 = (1 - amas**0.3333 * sighn / (8.0d0*pi*ra**2)
     $              * (amas-1.0d0) / amas
     $              * exp(-ra**2 * amas**0.6667 / (4.0*ram**2)))
            sigha = amas * sighn * tmp1
            tmp = tmp + du2 * pii / u2 * sigha * astep

   80    continue

         v2 = v2 + tmp

         if (tmp .le. 0.001*v2) goto 90

         u2 = u2 * 1.01d0

   70 continue

*-----------------------------------------------------------------------

      write(6,*) 'Maxshad is too Small'
      stop

*-----------------------------------------------------------------------

   90 continue

      sigga = sigvmda + 4.0d0 * pi * alem * ck * v2
      shadow = sigga / (amas * siggn)

*-----------------------------------------------------------------------
      return
      end
cABE end

************************************************************************
*                                                                      *
      function sigmaNRF(iz, in, ene)
*                                                                      *
*     Nuclear Resonance Fluorescence absorption cross-section          *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*        iz     : atomic number of target nucleus                      *
*        in     : neutron number of target nucleus                     *
*        ene    : incoming photon energy (MeV)                         *
*                                                                      *
*     output :                                                         *
*                                                                      *
*        sigmaNRF :  Cross-section in barns                            *
*                                                                      *
*                                                                      *
************************************************************************
      use levdat
cFURUTA20201007      use levdat1
      use NGSDATAMOD, only : spin, bindeg

      implicit real*8(a-h,o-z)

      include 'err.inc'
      include 'param02.inc'
cFURUTA20201007      common /gamint/ initga
      common /nrfmem/ spis, spgr, levabs, lflgnrf, mpole
!$OMP THREADPRIVATE(/nrfmem/)
      common /pnint/  ipnint

      logical lflg
      logical lflgbrch

      data lflgbrch /.true./

*-----------------------------------------------------------------------

      if(ipnint .le. 1) then  ! save CPU time when NRF is not concerned
       sigmaNRF = 0.d0
       return
      endif

      mpole    = 0
      levabs   = 0
      lflgnrf  = 0
      sigmaNRF = 0.d0
      ifound   = 0
      temp     = 300.d0 ! Target temperature in kelvin, ideally defined in the input.
      spgr     = spin(19, iz * 1000000 + iz + in)
      tarmas   = dble(iz) * rpmass + dble(in) * rnmass - bindeg(iz,in) ! target mass in MeV/c^2
      erecoil  = sqrt(ene**2 + tarmas**2) - tarmas ! Target recoil energy in MeV

*-----------------------------------------------------------------------
*     initialization move to setpar in read00.f
*-----------------------------------------------------------------------

cFURUTA20201007!$OMP CRITICAL (binin2_crit)
cFURUTA20201007         if( initga .eq. 0 ) then
cFURUTA20201007
cFURUTA20201007            initga = initga + 1
cFURUTA20201007            call levdatset(0)
cFURUTA20201007
cFURUTA20201007         end if
cFURUTA20201007!$OMP END CRITICAL (binin2_crit)

*-----------------------------------------------------------------------
*     Search absorption peaks defined in NRF special data
*-----------------------------------------------------------------------

      if((iz .eq. 82 .and. iz + in .eq. 206) .or.
     &   (iz .eq. 82 .and. iz + in .eq. 207) .or.
     &   (iz .eq. 82 .and. iz + in .eq. 208) .or.
     &   (iz .eq. 90 .and. iz + in .eq. 232) .or.
     &   (iz .eq. 92 .and. iz + in .eq. 235) .or.
     &   (iz .eq. 93 .and. iz + in .eq. 237) .or.
     &   (iz .eq. 94 .and. iz + in .eq. 239)  ) then
       call nrflibser(iz, in, ene, ifound, sigmaNRF, erecoil)
       if(ifound .eq. 1) return
      endif

*-----------------------------------------------------------------------
*     Search absorption peaks defined in ENSDF
*-----------------------------------------------------------------------

      nlevel = 0

      call levset(iz+in, iz, lflg)

      if(.not. lflg) goto 900

      labs = 0

      do i = 1, nlevel-1

       if( elevel(i)-ene .gt. 1.d-3) then
           exit
       elseif( elevel(i)-ene .lt. -1.d-3) then
           cycle
       endif

       if(spinpar(i,1) .eq. 0) cycle

       sdifmin = 10000
       do j = 1, nspv(i)
        if(abs(spinpar(i,j) - spgr) .le. sdifmin) then
        sdifmin = abs(spinpar(i,j) - spgr)
        spis    = spinpar(i,j)
        endif
       enddo

       if((spis .eq. 0.d0 .and. sdifmin .eq. 0.d0) .or.                  ! Prohibit 0 -> 0 transition.
     &       minval(nlevdn(i,1:ndch(i))) .ne. 0 ) cycle                  ! Exclude states irrelevant to ground state

       srange = 1.d-5
       SELECT CASE (iz)! expand search for close-to-magic nuclei because their NRF peaks are often wide-spread
        case(6:10,20-2:20+2,28-2:28+2,50-2:50+2,82-2:82+2,126-2:126+2)
        srange = 1.d-4
        case default

        SELECT CASE (in)
         case(6:10,20-2:20+2,28-2:28+2,50-2:50+2,82-2:82+2,126-2:126+2)
         srange = 1.d-4
        end select

       end select

       if( abs(ene - elevel(i) - erecoil) .le. 1.d-5 .and.
     & ( nlevdn(i,ndch(i)) .eq. 0 .or. sdifmin .le. 4.d0 )) then  ! absorption peak width is - 10 eV ~ + 10 eV
        eabs = elevel(i) + erecoil ! absorption peak energy
        labs = i                         ! absorption peak level number
        thal = thalf(i)
        exit
       endif
      enddo

      if(labs .eq. 0) goto 100

      iprgr = parit(0)

      mpole = int(max(1.d0, sdifmin)) * parit(labs) * iprgr
      if(spis .eq. 0.d0 .and. sdifmin .eq. 0.d0) mpole = 2   ! deexcitation from +0 to +0 is usually E2

      if(thal .eq. 0.d0) thal = thalcal(iz, in, labs) ! if ENSDF does not know half-life, estimate t-half based on theory

      gamma = 6.58211928d-22 / thal * Log(2.d0) ! Width in MeV   ! This is OK for eV-order half-life but what about longer half-lives?? check sometime

      if(sum(bratio(labs, 1:ndch(labs))) .eq. 0.d0 ) then ! Workaround for unexpected bug
        if(lflgbrch) then
         write(ErrCha,*)'Warning: Total decay width from NRF level is 0.
     & Z, A, E(MeV) = ', iz, iz+in, ene
         ErrID = 'L:16403/R:sigmaNRF/F:photnucl.f' !W00_003_001
         call ErrWrite(ErrID,ErrCha)
         lflgbrch = .false.
        endif
        gamma0   = gamma
      else
       gamma0  = gamma * max(bratio(labs, ndch(labs)), 1.d-3) ! sometimes branching ratio is missing
     & / sum(bratio(labs, 1:ndch(labs))) ! partial gamma to go to ground state
      endif

cABE 2017/11/24
      if(mpole .eq. 0) then
         sigmaNRF = 0.d0
         goto 100
      endif

      if(temp .eq. 0.d0) then  ! Absolute zero
       sigmaNRF = pi * (2.d0 * spis + 1.d0) / (2.d0 * spgr + 1.d0) /2.d0
     &    * (hbc * 1.d2 / ene) ** 2 * gamma * gamma0 / ( (ene - eabs)**2
     &    + (gamma * 0.5d0)**2 ) /2.d0 ! 1/2 comes from definition of Delta
      elseif(temp .ne. 0.d0) then  ! Finite temprature
       delta = 5.268d-6 * ene * sqrt(temp/ 300.d0 / dble(iz + in) )
       sigmaNRF = sqrt(pi/2.0)
     &    * (2.d0 * abs(spis) + 1.d0) / (4.d0 * abs(spgr) + 2.d0)
     &    * (hbc * 1.d2/ eabs)**2 * gamma0 * gamma / delta
     &    * bwinte(ene, eabs, gamma, delta)
      endif

      levabs = labs ! remember the excitation level

 100  call levUNset

*-----------------------------------------------------------------------
 900  return
      end

************************************************************************
*                                                                      *
      function bwinte(ene, ec, gamma, delta)
*                                                                      *
*     Numerical integration of Breit-Wigner function                   *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*        ec     : Energy of resonance peak (MeV)                       *
*        gamma  : Gamma of Breit-Wigner (MeV)                          *
*                                                                      *
*     output :                                                         *
*                                                                      *
*        bwinte : Breit-Wigner integral  (MeV^-3)                      *
*                                                                      *
*                                                                      *
************************************************************************
      implicit real*8(a-h,o-z)

      bwinte = 0.d0

      do i = -5000, 5000

       e = ec + gamma/1.d2 * i
       dbw = 1.d0/( (e-ec)**2 + (gamma/2.d0)**2 ) *
     & exp( -(e-ene)**2 / ( 2.d0 * delta**2 ) ) * gamma/1.d2

       bwinte = bwinte + dbw
      enddo

      end

************************************************************************
*                                                                      *
      subroutine NRF_ang(axf, ayf, azf, iz, in)
*                                                                      *
*     Nuclear Resonance Fluorescence Gamma emission angle              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*        ax0, ay0, az0  : incoming photon angular cosine               *
*                         usually ax = ay =0, az = 1                   *
*     output :                                                         *
*                                                                      *
*        axf, ayf, azf :  outgoing photon angular cosine               *
*                                                                      *
*                                                                      *
************************************************************************
      use MMBANKMOD ! get spin info  i.e., (spx, spy, spz)

      implicit real*8(a-h,o-z)

      include 'param02.inc'
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /nrfmem/ spis, spgr, levabs, lflgnrf, mpole
!$OMP THREADPRIVATE(/nrfmem/)
      common /nrfdir/ usave(3)
!$OMP THREADPRIVATE(/nrfdir/)

*-----------Noise cut---------------------------------------------
      ax0 = usave(1)  ! renaming for simplicity. usave is incoming direction
      ay0 = usave(2)
      az0 = usave(3)

      if(abs(ax0) .le. 1.d-15) ax0 = 0.d0
      if(abs(ay0) .le. 1.d-15) ay0 = 0.d0
      if(abs(az0) .le. 1.d-15) az0 = 0.d0

      ax = ax0 / sqrt(ax0**2 + ay0**2 + az0**2)
      ay = ay0 / sqrt(ax0**2 + ay0**2 + az0**2)
      az = az0 / sqrt(ax0**2 + ay0**2 + az0**2)

      the  = 0.d0
      phi  = 0.d0
      cos1 = 0.d0
      sin1 = 0.d0
      ww   = 0.d0
      r    = 0.d0
      fac  = 0.d0
      fac1 = 0.d0
      fac2 = 0.d0

*-----------------------------------------------------------------

      if(mod(iz,2) .eq. 0 .and. mod(in,2) .eq. 0) then! even-even
       nsym = 0
      elseif(mod(iz+in,2) .eq. 0) then                ! odd-odd
       nsym = 1
      else                                            ! even-odd
       nsym = 2
      endif

*------------ judge polarization ---------------------------------

!--- change NS 2020.04 del THREADPRIVATE
!      if(sqrt(spx(no)**2 + spy(no)**2 + spz(no)**2) .le. 1.d-8) then ! not polarized
      if(sqrt(spx(no,ipomp+1)**2 + spy(no,ipomp+1)**2 +
     &                         spz(no,ipomp+1)**2) .le. 1.d-8) then ! not polarized
!--- end change NS 2020.04 del THREADPRIVATE
       ipola = 0
      else  ! polarized
       ipola = 1
      endif

*-----Case 1, Even-even nuclei------------------------------------
      IF(nsym .eq. 0) then

       if(ipola .eq. 0) then ! unpolarized photon beam
        if(abs(mpole) .eq. 1) then ! Dipole emission   E1, M1
        r = unirn(dummy)
        cos1 =
     &   (sqrt(16.d0 * r**2 - 16.d0 * r + 5.d0)-4.d0*r+2.d0)**0.333d0
     &  -(sqrt(16.d0 * r**2 - 16.d0 * r + 5.d0)-4.d0*r+2.d0)**(-0.333d0)
        sin1 = sign(sqrt(1.d0 - cos1**2), unirn(dummy) - 0.5d0)
        phi  = 2.d0 * pi * unirn(dummy)
        elseif(abs(mpole) .eq. 2) then ! Quadrepole emission   E2, M2
         do ic = 1, 100
          the  = pi * unirn(dummy)
          cos1 = cos(the)
          sin1 = sin(the)
          ww   = 1.25d0 * sin1 *(1.d0 - 3.d0 * cos1**2 + 4.d0 * cos1**4)
          if(2.5d0 * unirn(dummy) .ge. ww) exit
         enddo
c          sin1  = sign(sqrt(1.d0 - cos1**2), 2.d0 * unirn(dummy) - 1.d0)
           phi  = 2.d0 * pi * unirn(dummy)
        endif
       elseif(ipola .eq. 1) then
        if(abs(mpole) .eq. 1) then ! Dipole emission   mpole = -1 -> E1.  mpole = +1 -> M1.
        do ic = 1, 100
         the  = pi * unirn(dummy)
         cos1 = cos(the)
         sin1 = sin(the)
         phi  = 2.d0 * pi * unirn(dummy)
         ww   = sin1 * (0.75d0 * (1.d0 + cos1**2) + dble(mpole) * 0.75d0
     &         * (1.d0 - cos1**2) * cos(2.d0 * phi))
         if(1.5d0 * unirn(dummy) .ge. ww) exit
        enddo
        elseif(abs(mpole) .eq. 2) then ! Quadrupole emission   mpole = +2 -> E2.  mpole = -2 -> M1.
        do ic = 1, 100
         the  = pi * unirn(dummy)
         cos1 = cos(the)
         sin1 = sin(the)
         phi  = 2.d0 * pi * unirn(dummy)
         ww   = sin1 * (1.25d0 *(1.d0 - 3.d0 * cos1**2 + 4.d0 * cos1**4)
     &               +  1.25d0 *(1.d0 - 5.d0 * cos1**2 + 4.d0 * cos1**4)
     &       * dble(mpole/2) * cos(2.d0 * phi))
         if(2.5d0 * unirn(dummy) .ge. ww) exit
        enddo
        endif
       endif

       if(abs(mpole) .gt. 2) then ! higher polarity emission goes to isotropic
         cos1 = 2 * unirn(dummy) - 1.d0
         sin1 = sign(sqrt(1.d0 - cos1**2), unirn(dummy) - 0.5d0)
          phi  = 2.d0 * pi * unirn(dummy)
       endif

*-----Case 2, Odd-Odd nuclei------------------------------------
      ELSEIF(nsym .eq. 1) then  ! currently, symmetric emission

         cos1 = 2 * unirn(dummy) - 1.d0
         sin1 = sign(sqrt(1.d0 - cos1**2), unirn(dummy) - 0.5d0)
          phi  = 2.d0 * pi * unirn(dummy)


*-----Case 3, Even-Odd nuclei------------------------------------
      ELSEIF(nsym .eq. 2) then

       if(abs(mpole) .eq. 1) then  ! E1, M1 transition
        SelectCase(int(2.d0*spgr))
         Case(7) ! 3.5d0
          SelectCase(int(2.d0*spis))
          Case(9) ! 4.5d0
           fac = 0.3028d0
          Case(7) ! 3.5d0
           fac = -0.4364d0
          Case(5) ! 2.5d0
           fac = 0.1336d0
          EndSelect
         Case(5) ! 2.5d0  ! calculate someday. Isotopic assumption
           fac = 0.d0
         Case(3) ! 1.5d0
          SelectCase(int(2.d0*spis))
          Case(5) ! 2.5d0
           fac = 0.3742d0
          Case(3) ! 1.5d0
           fac = 0.4d0
          Case(1) ! 0.5d0
           fac = 0.d0
          endselect
         Case(1) ! 0.5d0
          SelectCase(int(2.d0*spis))
          Case(3) ! 1.5d0
           fac = 0.5d0
          Case(1) ! 0.5d0
           fac = 0.d0
          endselect
        EndSelect

        do ic = 1, 100
         the  = pi * unirn(dummy)
         cos1 = cos(the)
         sin1 = sin(the)
         phi  = 2.d0 * pi * unirn(dummy)
         ww   = sin1 * (1.d0 + fac**2 * (0.5d0 * (3.d0 * cos1**2 -1)
     &   + dble(mpole/2) * 1.5d0 * (1.d0 - cos1**2) * cos(2.d0 * phi) ))
         if( (1.d0 + fac**2) * unirn(dummy) .ge. ww ) exit
        enddo

       elseif(abs(mpole) .eq. 2) then  ! E2, M2 transition

        SelectCase(int(2.d0*spgr))
         Case(7) ! 3.5d0
          SelectCase(int(2.d0*spis))
          Case(11) ! 5.5d0
           fac1 = -0.4109d0
           fac2 = -0.2237d0
          Case(3) ! 1.5d0
           fac1 = -0.1429d0
           fac2 = 0.d0
          EndSelect
         Case(5) ! 2.5d0 ! calculate someday. Isotropic assumption
           fac1 = 0.d0
           fac2 = 0.d0
         Case(3) ! spgr = 1.5d0  ! spis .eq. 3.5
           fac1 = -0.4676d0
           fac2 = -0.3582d0
         Case(1) ! spgr = 0.5d0  ! spis .eq. 2.5
           fac1 = -0.5345d0
           fac2 = -0.6172d0
        EndSelect

        do ic = 1, 100
         the  = pi * unirn(dummy)
         cos1 = cos(the)
         sin1 = sin(the)
         phi  = 2.d0 * pi * unirn(dummy)
         ww   = sin1 * (1.d0
     &   + fac1**2 * (0.5d0 * (3.d0 * cos1**2 -1)
     &   + dble(mpole/2) * 1.5d0 * (1.d0 - cos1**2) * cos(2.d0 * phi) )
     &   + fac2**2 * (0.125d0 * (35.d0 * cos1**4 - 30.d0 * cos1**2+3.d0)
     &   - dble(mpole/2) / 12.d0 * (-52.5d0 * cos1**4 + 60.d0 * cos1**2
     &   - 7.5d0) * cos(2.d0 * phi) ))
         if( (1.d0 + fac1**2 + fac2**2)* unirn(dummy) .ge. ww ) exit
        enddo

       endif

       if(abs(mpole) .gt. 2) then ! higher polarity emission goes to isotropic
         cos1 = 2 * unirn(dummy) - 1.d0
         sin1 = sign(sqrt(1.d0 - cos1**2), unirn(dummy) - 0.5d0)
          phi  = 2.d0 * pi * unirn(dummy)
       endif

      ENDIF

*--------- Determination of emission angle -----------------------

!--- change NS 2020.04 del THREADPRIVATE
!      dj = (ay * spz(no) - az * spy(no)) + (az * spx(no) - ax * spz(no))
!     &   + (ax * spy(no) - ay * spx(no))    ! judge phi rotation based on outer product of incident vector and spin vector
      dj = (ay * spz(no,ipomp+1) - az * spy(no,ipomp+1))
     &   + (az * spx(no,ipomp+1) - ax * spz(no,ipomp+1))
     &   + (ax * spy(no,ipomp+1) - ay * spx(no,ipomp+1))    ! judge phi rotation based on outer product of incident vector and spin vector
!--- change NS 2020.04 del THREADPRIVATE
       if(dj .lt. 0.d0) phi = phi + pi/2.d0

      if(az .ne. 1.d0) then ! except z beam
       axf = ax * cos1 * cos(phi) - ay * cos1 * sin(phi)
     &  -az * sin1 *
     &  (cos(phi) * ax/sqrt(1-az**2) - sin(phi) * ay/sqrt(1-az**2) )
       ayf = ay * cos1 * cos(phi) + ax * cos1 * sin(phi)
     &  - az * sin1 *
     &  (cos(phi) * ay/sqrt(1-az**2) + sin(phi) * ax/sqrt(1-az**2) )
       azf = az * cos1 + sqrt(ax**2 + ay**2) * sin1 *
     &  sign(1.d0, unirn(dummy) - 0.5d0)
      else ! z beam (y,x components are 0)
       axf = sin1 * cos(phi)
       ayf = sin1 * sin(phi)
       azf = - cos1
      endif

! renormalize vector. photon energy was increased 2024/2/15
      azabs = sqrt(axf**2 + ayf**2 + azf**2)
      if(azabs .eq. 0.d0) then
       axf   = 0.d0
       ayf   = 0.d0
       azf   = 1.d0
      else
       axf   = axf/azabs
       ayf   = ayf/azabs
       azf   = azf/azabs
      endif
*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      function thalcal(iz, in, labs)
*                                                                      *
*     Theoretical calculation of excitation state half life            *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*        iz     : atomic number of target nucleus                      *
*        in     : neutron number of target nucleus                     *
*        labs   : level number                                         *
*                                                                      *
*     output :                                                         *
*                                                                      *
*        thalcal : half life of the excited stat (s)                   *
*                                                                      *
************************************************************************
      use levdat
      use NGSDATAMOD, only : spin

      implicit real*8(a-h,o-z)

      include 'param02.inc'

      thalcal = 1.d-15
      iflgnrf = 1        ! 0 to use this function for non NRF
      tinv = 0.d0        ! S.Abe 2017/11/24

*-----------------------------------------------------------------------
*     Search destination levels
*-----------------------------------------------------------------------
      igrd = 0 ! flag to check if ground is included

      do i = 1, ndch(labs)
       if(nlevdn(labs,i) .eq. 0) then
        sdes = spin(19, iz*1000000 + iz+in)
        igrd = 1
       else
        sdes = spinpar(nlevdn(labs,i),1)
       endif
        l  = int(abs( spinpar(labs,1) - sdes )) ! lambda. Spin shift
       if(l .eq. 0 .and. sdes .eq. 0.d0) then
        l = 2
       elseif(l .eq. 0) then
        l = 1
       endif
        ips = parit( nlevdn(labs,i) ) * parit(labs) ! parity shift. Shift or stay

       if( (-1)**l * ips .ge. 0 )then ! Electric mode
        wsp = 4.4d0 * dble(l + 1) / dble(l) *
     & (3.d0 / dble(l+3))**2 * (2.d0 * sdes + 1.d0) *1.d21
     & * canfact3(elevel(labs), elevel(nlevdn(labs,i)),l, iz + in)
       else      ! Magnetic mode
        wsp = 1.9d0 * dble(l + 1) / dble(l) *
     & (3.d0 / dble(l+3))**2 * (2.d0 * sdes + 1.d0) *1.d21
     & * (1.2d0 * dble(iz + in)**0.33333d0)**(-2)
     & * canfact3(elevel(labs), elevel(nlevdn(labs,i)),l, iz + in)
     & / 50.d0  ! unknown factor to explain U235 and Pu239
       endif

      tinv = tinv + wsp
      enddo

*-----------------------------------------------------------------------
*     Include ground state as a destination
*-----------------------------------------------------------------------

      if(igrd .eq. 0 .and. iflgnrf .eq. 1) then ! include ground as destination
       l = minval(int(abs( spinpar(labs,1:nspv(labs)) -
     &  spin(19, iz*1000000 + iz+in) )))
       if(l .eq. 0 .and. sdes .eq. 0.d0) then
        l = 2
       elseif(l .eq. 0.d0) then
        l = 1
       endif
        ips = parit(0) * parit(labs) ! parity shift. Shift or stay
       if( (-1)**l * ips .ge. 0 )then ! Electric mode
        wsp = 4.4d0 * dble(l + 1) / dble(l) *
     & (3.d0 / dble(l+3))**2 * (2.d0 * sdes + 1.d0) *1.d21
     & * canfact3(elevel(labs), 0.d0, l, iz + in)
       else      ! Magnetic mode
        wsp = 1.9d0 * dble(l + 1) / dble(l) *
     & (3.d0 / dble(l+3))**2 * (2.d0 * sdes + 1.d0) *1.d21
     & * (1.2d0 * dble(iz + in)**0.33333d0)**(-2)
     & * canfact3(elevel(labs), 0.d0, l, iz + in)
     & / 50.d0  ! unknown factor to explain U235 and Pu239
       endif
      tinv = tinv + wsp
      endif

cABE 2017/11/24
c      thalcal = tinv ** (-1)
      if(tinv .ne. 0.d0) thalcal = tinv ** (-1)

*-----------------------------------------------------------------------
      return
      end

************************************************************************
*                                                                      *
      subroutine nrflibser(iz, in, ene, ifound, sigma, erecoil)
*                                                                      *
*     Nuclear Resonance Fluorescence sigma for special isotopes        *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*        iz, in : nuclear charge and neutron number                    *
*        ene    : incoming photon energy (MeV)                         *
*     erecoil   : Recoiling of primary particle (MeV)                  *
*                                                                      *
*     output :                                                         *
*                                                                      *
*        sigma  : calculated cross-section (b)                         *
*        levabs : gamma absorption level. in this subrouitne,          *
*                  value is negative to avoid mixing with gammod levels*
*        ifound : flag to remember if level is found in this subroutine*
*        mpole  : emission multipolarity (1, 2, ...)                   *
*                                                                      *
************************************************************************

      implicit doubleprecision(a-h,o-z)

      include 'param02.inc'

      parameter ( lnrfth232 = 25 )
      parameter ( lnrfu235  = 18 )
      parameter ( lnrfnp237 = 15 )
      parameter ( lnrfp239  = 12 )
      parameter ( lnrfpb206 = 4 )
      parameter ( lnrfpb207 = 2 )
      parameter ( lnrfpb208 = 2 )

      common /NRFTH232/ spnrf3(lnrfth232), csnrf3(lnrfth232),
     & enrf3(lnrfth232), ipnrf3(lnrfth232)
      common /NRFU235/ spnrf1(lnrfu235), csnrf1(lnrfu235),
     & enrf1(lnrfu235), ipnrf1(lnrfu235)
      common /NRFNP237/ spnrf4(lnrfnp237), csnrf4(lnrfnp237),
     & enrf4(lnrfnp237), ipnrf4(lnrfnp237)
      common /NRFP239/ spnrf2(lnrfp239), csnrf2(lnrfp239),
     & enrf2(lnrfp239), ipnrf2(lnrfp239)
      common /NRFPB206/ spnrf5(lnrfpb206), csnrf5(lnrfpb206),
     & enrf5(lnrfpb206), ipnrf5(lnrfpb206)
      common /NRFPB207/ spnrf6(lnrfpb207), csnrf6(lnrfpb207),
     & enrf6(lnrfpb207), ipnrf6(lnrfpb207)
      common /NRFPB208/ spnrf7(lnrfpb208), csnrf7(lnrfpb208),
     & enrf7(lnrfpb208), ipnrf7(lnrfpb208)

      common /nrfmem/ spis, spgr, levabs, lflgnrf, mpole
!$OMP THREADPRIVATE(/nrfmem/)

! Level translation table from enrf* to RIPL.
      integer levtrnspb206(lnrfpb206)
      integer levtrnspb207(lnrfpb207)
      integer levtrnspb208(lnrfpb208)

      data levtrnspb206 /164,190,199,204/
      data levtrnspb207 /122,138/
      data levtrnspb208 /82,103/

      temp     = 300.d0 ! Target temperature in kelvin, ideally defined in the input.
*-----------------------------------------------------------------------

      if(iz .eq. 92 .and. iz + in .eq. 235) then
       do i = 1, lnrfu235
        if(abs(ene - enrf1(i)*1.d-3 - erecoil) .le. 1.d-5) then
         levabs = -i
         ifound = 1
         mpole  = 1
         delta = 5.268d-6 * ene * sqrt(temp/ 300.d0 / dble(iz + in) ) ! in MeV
         sigma = exp(- ((ene - enrf1(i)*1.d-3 - erecoil)**2 / (2.d0 *
     &  delta**2)) ) * csnrf1(i) / (sqrt(2.d0 * pi) * delta * 1.d6)
        return
        endif
       enddo
      elseif(iz .eq. 94 .and. iz + in .eq. 239) then
        do i = 1, lnrfp239
        if(abs(ene - enrf2(i)*1.d-3 - erecoil) .le. 1.d-5) then
         levabs = -i
         ifound = 1
         mpole  = 1
         delta = 5.268d-6 * ene * sqrt(temp/ 300.d0 / dble(iz + in) ) ! in MeV
         sigma = exp(- ((ene - enrf2(i)*1.d-3 - erecoil)**2 / (2.d0
     &   * delta**2)) )  * csnrf2(i) / (sqrt(2.d0 * pi) * delta * 1.d6)
        return
        endif
       enddo
      elseif(iz .eq. 90 .and. iz + in .eq. 232) then
        do i = 1, lnrfth232
        if(abs(ene - enrf3(i)*1.d-3 - erecoil) .le. 1.d-5) then
         levabs = -i
         ifound = 1
         mpole  = 1
         delta = 5.268d-6 * ene * sqrt(temp/ 300.d0 / dble(iz + in) ) ! in MeV
         sigma = exp(- ((ene - enrf3(i)*1.d-3 - erecoil)**2 / (2.d0
     &   * delta**2)) )  * csnrf3(i) / (sqrt(2.d0 * pi) * delta * 1.d6)
        return
        endif
       enddo
      elseif(iz .eq. 93 .and. iz + in .eq. 237) then
        do i = 1, lnrfnp237
        if(abs(ene - enrf4(i)*1.d-3 - erecoil) .le. 1.d-5) then
         levabs = -i
         ifound = 1
         mpole  = 1
         delta = 5.268d-6 * ene * sqrt(temp/ 300.d0 / dble(iz + in) ) ! in MeV
         sigma = exp(- ((ene - enrf4(i)*1.d-3 - erecoil)**2 / (2.d0
     &   * delta**2)) )  * csnrf4(i) / (sqrt(2.d0 * pi) * delta * 1.d6)
        return
        endif
       enddo
      elseif(iz .eq. 82 .and. iz + in .eq. 206) then
        do i = 1, lnrfpb206
        if(abs(ene - enrf5(i)*1.d-3 - erecoil) .le. 1.d-5) then
         levabs = -levtrnspb206(i)
         ifound = 1
         mpole  = 1
         delta = 5.268d-6 * ene * sqrt(temp/ 300.d0 / dble(iz + in) ) ! in MeV
         sigma = exp(- ((ene - enrf5(i)*1.d-3 - erecoil)**2 / (2.d0
     &   * delta**2)) )  * csnrf5(i) / (sqrt(2.d0 * pi) * delta * 1.d6)
        return
        endif
       enddo
      elseif(iz .eq. 82 .and. iz + in .eq. 207) then
        do i = 1, lnrfpb207
        if(abs(ene - enrf6(i)*1.d-3 - erecoil) .le. 1.d-5) then
         levabs = -levtrnspb207(i)
         ifound = 1
         mpole  = 1
         delta = 5.268d-6 * ene * sqrt(temp/ 300.d0 / dble(iz + in) ) ! in MeV
         sigma = exp(- ((ene - enrf6(i)*1.d-3 - erecoil)**2 / (2.d0
     &   * delta**2)) )  * csnrf6(i) / (sqrt(2.d0 * pi) * delta * 1.d6)
        return
        endif
       enddo
      elseif(iz .eq. 82 .and. iz + in .eq. 208) then
        do i = 1, lnrfpb208
        if(abs(ene - enrf7(i)*1.d-3 - erecoil) .le. 1.d-5) then
         levabs = -levtrnspb208(i)
         ifound = 1
         mpole  = 1
         delta = 5.268d-6 * ene * sqrt(temp/ 300.d0 / dble(iz + in) ) ! in MeV
         sigma = exp(- ((ene - enrf7(i)*1.d-3 - erecoil)**2 / (2.d0
     &   * delta**2)) )  * csnrf7(i) / (sqrt(2.d0 * pi) * delta * 1.d6)
        return
        endif
       enddo
      endif



*-----------------------------------------------------------------------
      return
      end



************************************************************************
*                                                                      *
      block data nrflevel
*    NRF level information
*    ipnrf* : parity +1 or -1
*    spnrf* : spin
*    csnrf* : cross section (eV b)
*    ennrf* : energy (MeV)

      implicit doubleprecision(a-h,o-z)

      parameter ( lnrfth232 = 25 )
      parameter ( lnrfu235  = 18 )
      parameter ( lnrfnp237 = 15 )
      parameter ( lnrfp239  = 12 )
      parameter ( lnrfpb206 = 4 )
      parameter ( lnrfpb207 = 2 )
      parameter ( lnrfpb208 = 2 )

      common /NRFTH232/ spnrf3(lnrfth232), csnrf3(lnrfth232),
     & enrf3(lnrfth232), ipnrf3(lnrfth232)
      common /NRFU235/ spnrf1(lnrfu235), csnrf1(lnrfu235),
     & enrf1(lnrfu235), ipnrf1(lnrfu235)
      common /NRFNP237/ spnrf4(lnrfnp237), csnrf4(lnrfnp237),
     & enrf4(lnrfnp237), ipnrf4(lnrfnp237)
      common /NRFP239/ spnrf2(lnrfp239), csnrf2(lnrfp239),
     & enrf2(lnrfp239), ipnrf2(lnrfp239)
      common /NRFPB206/ spnrf5(lnrfpb206), csnrf5(lnrfpb206),
     & enrf5(lnrfpb206), ipnrf5(lnrfpb206)
      common /NRFPB207/ spnrf6(lnrfpb207), csnrf6(lnrfpb207),
     & enrf6(lnrfpb207), ipnrf6(lnrfpb207)
      common /NRFPB208/ spnrf7(lnrfpb208), csnrf7(lnrfpb208),
     & enrf7(lnrfpb208), ipnrf7(lnrfpb208)

c     PHYSICAL REVIEW C 83, 034615 (2011), "Discovery of low-lying E1 and M1 strengths in 232 Th"
      data (ipnrf3(i), spnrf3(i), csnrf3(i), enrf3(i), i= 1, lnrfth232)/
     &   1,  1.0d0, 46.3d0,  2043.7d0,  -1,  1.0d0,  2.1d0,  3060.4d0,
     &   1,  1.0d0, 25.5d0,  2249.5d0,  -1,  1.0d0,  0.9d0,  3395.8d0,
     &   1,  1.0d0, 19.0d0,  2296.3d0,  -1,  1.0d0,  2.2d0,  3607.9d0,
     &   1,  1.0d0,  4.4d0,  2795.2d0,  -1,  1.0d0,  1.9d0,  3626.3d0,
     &   1,  1.0d0,  6.2d0,  2835.0d0,  -1,  1.0d0,  1.1d0,  3639.1d0,
     &   1,  1.0d0,  3.4d0,  2865.5d0,  -1,  1.0d0,  2.0d0,  3731.5d0,
     &   1,  1.0d0,  3.1d0,  2885.8d0,  -1,  1.0d0,  1.6d0,  3742.6d0,
     &   1,  1.0d0,  3.5d0,  2924.6d0,  -1,  1.0d0,  1.3d0,  3752.4d0,
     &   1,  1.0d0,  4.7d0,  2935.6d0,  -1,  1.0d0,  1.5d0,  3820.7d0,
     &   1,  1.0d0,  2.3d0,  2996.3d0,  -1,  1.0d0,  1.6d0,  3920.8d0,
     &   1,  1.0d0,  3.0d0,  3014.4d0,  -1,  1.0d0,  1.1d0,  3935.8d0,
     &   1,  1.0d0,  7.7d0,  3115.1d0,  -1,  1.0d0,  3.0d0,  4002.2d0,
     &   1,  1.0d0,  1.9d0,  3287.1d0/

c     PHYSICAL REVIEW C 83, 041601(R) (2011), "Discrete deexcitations in 235 U below 3 MeV from nuclear resonance fluorescence"
      data (ipnrf1(i), spnrf1(i), csnrf1(i), enrf1(i), i=  1, lnrfu235)/
     &  -1,  3.5d0,  3.0d0,  1656.3d0,  -1,  4.5d0,  2.2d1,  1733.6d0,
     &  -1,  3.5d0,  6.4d0,  1769.3d0,  -1,  4.5d0,  8.9d0,  1815.2d0,
     &  -1,  4.5d0,  5.5d0,  1827.7d0,  -1,  3.5d0,  9.6d0,  1862.4d0,
     &  -1,  3.5d0,  4.6d0,  1973.8d0,  -1,  4.5d0,  6.7d0,  2003.3d0,
     &  -1,  3.5d0,  4.6d0,  2005.9d0,  -1,  3.5d0,  3.0d0,  2010.6d0,
     &  -1,  3.5d0,  3.0d0,  2067.1d0,  -1,  3.5d0,  1.4d0,  2074.2d0,
     &  -1,  3.5d0,  1.1d0,  2086.7d0,  -1,  4.5d0,  3.0d0,  2110.2d0,
     &  -1,  3.5d0,  2.8d0,  2216.1d0,  -1,  3.5d0,  3.6d0,  2416.1d0,
     &  -1,  3.5d0,  2.5d0,  2555.6d0,  -1,  3.5d0,  3.6d0,  2754.7d0/

c     PHYSICAL REVIEW C 82, 054310 (2010), "Nuclear resonance fluorescence of 237 Np"
      data (ipnrf4(i), spnrf4(i), csnrf4(i), enrf4(i), i= 1, lnrfnp237)/
     &   0,  0.0d0,  5.3d0,  1697.8d0,   0,  0.0d0,  7.8d0,  2261.5d0,
     &   0,  0.0d0, 10.6d0,  1728.8d0,   0,  0.0d0,  4.6d0,  2288.2d0,
     &   0,  0.0d0,  4.8d0,  1739.7d0,   0,  0.0d0,  4.9d0,  2375.5d0,
     &   0,  0.0d0,  5.7d0,  1827.1d0,   0,  0.0d0,  2.7d0,  2378.0d0,
     &   0,  0.0d0,  5.6d0,  1861.9d0,   0,  0.0d0,  2.7d0,  2381.5d0,
     &   0,  0.0d0,  6.1d0,  1925.9d0,   0,  0.0d0,  4.4d0,  2402.3d0,
     &   0,  0.0d0,  7.2d0,  2179.5d0,   0,  0.0d0,  3.8d0,  2506.0d0,
     &   0,  0.0d0,  7.3d0,  2252.0d0/

c     PHYSICAL REVIEW C 78, 041601(R) (2008), "Nuclear resonance fluorescence excitations near 2 MeV in 235U and 239Pu"
c     parity, spin is not provided. ggnrf is Is (eV b)
      data (ipnrf2(i), spnrf2(i), csnrf2(i), enrf2(i), i=  1, lnrfp239)/
     &   1,  1.5d0,  8.d0,  2040.25d0, 1,  1.5d0,  5.0d0,  2046.89d0,
     &   1,  1.5d0,  4.d0,  2135.00d0, 1,  1.5d0,  1.3d1,  2143.56d0,
     &   1,  1.5d0,  5.d0,  2150.98d0, 1,  1.5d0,  8.0d0,  2289.02d0,
     &   1,  1.5d0,  1.d1,  2423.48d0, 1,  1.5d0,  9.0d0,  2431.66d0,
     &   1,  1.5d0,  9.d0,  2454.37d0, 1,  1.5d0,  6.0d0,  2460.46d0,
     &   1,  1.5d0,  8.d0,  2464.60d0, 1,  1.5d0,  6.0d0,  2471.07d0/

c     From Prof. Ohgaki 2024/9/1
      data (ipnrf5(i), spnrf5(i), csnrf5(i), enrf5(i), i= 1, lnrfpb206)/
     &  -1,  1.0d0, 2.44d0, 5038.57d0, 0,  1.0d0,0.146d0,  5378.20d0,
     &  -1,  1.0d0,0.171d0, 5471.90d0, 0,  1.0d0,0.203d0,  5525.20d0/

      data (ipnrf6(i), spnrf6(i), csnrf6(i), enrf6(i), i= 1, lnrfpb207)/
     &   1,  1.5d0,0.456d0,  5217.58d0, 1,  1.5d0,0.709d0,  5489.77d0/

      data (ipnrf7(i), spnrf7(i), csnrf7(i), enrf7(i), i= 1, lnrfpb208)/
     &  -1,  1.0d0,1.99d0,  5291.90d0,-1,  1.0d0, 3.09d0,  5511.78d0/

      end
************************************************************************
