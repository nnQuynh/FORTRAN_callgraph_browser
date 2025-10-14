************************************************************************
*                                                                      *
      module ets_art
*                                                                      *
*                                                                      *
*   ETSART                                                             *
*  (Electron Track Structure model for Arbitrary Targets)     *
*                                                                      *
*                                                                      *
************************************************************************
      use ELEDATAMOD
      use MMBANKMOD
      use moddas_material
      use GGMARRAYMOD, only : jemi !  gas/condenced switch
      implicit double precision(a-h, o-z)

      include 'param-physcnst.inc'
      include 'param.inc'
      include 'param00.inc'

*---------- common variables -----------------------------------*
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /eparm/  esmax, esmin, emin(20)

      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /kmat1g/ kmat(kvlmax)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /tscmsg/ ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /celepcc/ enumpcc(kvlmax)

! ielas can be used for ignoring the elastic scattering (developer option)
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common / clustt / nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)
      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE( /clustf/)
      common /dcayp/  adcayp(20), bdcayp(20)
      common / idumc /idum
      common /tstara/ atmrc(10,8)
!$OMP THREADPRIVATE(/tstara/)

! ets parameters
      common / etsminmax / etsmin, etsmax
      common / ele1st / u1st,v1st,w1st
!$OMP THREADPRIVATE(/ele1st/)
      common / ele2nd  / e2nd,u2nd,v2nd,w2nd
!$OMP THREADPRIVATE(/ele2nd/)
      common / etsexe  / dexc_ene
!$OMP THREADPRIVATE(/etsexe/)
      common /etsion/ icoll
!$OMP THREADPRIVATE(/etsion/)
      common / cschng / change, fit
      common /etsart/ bgets(kvlmax), ebgets(kvlmax),
     &                wvets(kvlmax), ewvets(kvlmax)
      common /ets_wvalue/ wvlin
!$OMP THREADPRIVATE(/ets_wvalue/)

*---------------------------------------------------------------*


      double precision, allocatable, private :: xsec_shel(:,:)
!$OMP THREADPRIVATE(xsec_shel)
      double precision, private :: eprob(3)
!$OMP THREADPRIVATE(eprob)
      double precision, public  :: xels,xion,xexc
!$OMP THREADPRIVATE(xels,xion,xexc)
      double precision, public  :: zp ! effective atomic number
!$OMP THREADPRIVATE(zp)

      double precision, private, parameter :: slbnd = 1.d-20   !
      parameter (pi = 3.141592653589793d0 )

      double precision, private :: ele_period(1:108)
!$OMP THREADPRIVATE(ele_period)

      double precision, private :: bexc(100), uexc(100), dnex(100),
     &                             blow(100), denex(100)
!$OMP THREADPRIVATE(bexc,uexc,dnex,blow,denex)

      integer nexc,iex2,isp
      double precision, private :: bg, eplasm,ULave,sumvale,sumdens
      double precision, private :: bimin
!$OMP THREADPRIVATE(nexc,iex2,isp,bg,eplasm,ULave,sumvale,sumdens,bimin)

      integer mexel(100)
      integer msec, iblz000
!$OMP THREADPRIVATE(mexel,msec, iblz000)

! subshel parameters
      data (ele_period(i),i= 1, 108)/
     &    2*1.d0, 8*2.d0,8*3.d0,18*4.d0,18*5.d0,32*6.d0,22*7.d0/
      double precision, private :: shel_priod(1:28)
      data (shel_priod(i),i=1, 28)/
!          K    L      M      N      O      P      Q
     &    1.d0,3*2.d0,5*3.d0,7*4.d0,7*5.d0,4*6.d0,7.d0/


      contains


c=================================================================
      subroutine etsflt_art(sig_macro)
c=================================================================
      implicit real*8(a-h,o-z)


      ec(ibkec+no,ipomp+1) = ec(ibkec+no,ipomp+1) * 1.d+6   ! from (MeV) to (eV)
      e(ibke+no,ipomp+1)   = e(ibke+no,ipomp+1)   * 1.d+6   ! from (MeV) to (eV)
      ene = e(ibke+no,ipomp+1)

c // count the subshell for determination of size of allay
      sig_macro = 0.d0
      nelems = 0
      itzp = 0
       hydro = denh_das(kmat0+mat)
      if(hydro .gt. 0.d0) then
        nelems = 1
        itzp = 1
      endif
      lem   = nint( dnel_das(kmat0+mat) )
      nelems= nelems + lem  ! 20231109 for hydro>0.d0

      if(nelems.ne.0) then
        if( .not. allocated(xsec_shel)) then
            allocate(xsec_shel(1:28,nelems))
            xsec_shel = 0.d0
        else
            deallocate(xsec_shel)
            allocate(xsec_shel(1:28,nelems))
            xsec_shel = 0.d0
        endif
      endif

c // cross section (cm2) //
c--------------------------------
c   elastic scattering
c--------------------------------
      xels = 0.d0
      zp = 0.d0
      call  els_ics_art(ene)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! elastic scattering off for calculate the number of generated electrons
      if(msec.eq.2 .or. ielas.eq.0) then
        xels0= xels
        xels = 0.d0
      endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c--------------------------------
c   ionization
c--------------------------------
      xion = 0.d0
      call  ion_ics_art(ene)

c--------------------------------
c   excitation
c--------------------------------
      xexc = 0.d0
      if(iblz000.ne.iblz1)  call  exc_setup
      iblz000 = iblz1
      if(nexc.ne.0)   call  exc_ics_art(ene)

c // total cross section (sgm) //
! elastice scattering process will be killed when the electron energy is below plasmon peak. 20230422
      if(ene.lt.eplasm .and. (xion+xexc).ne.0.d0) then
         xels = 0.d0
      endif

      sgm = xels + xion + xexc
! avoiding error for termination of electron transport
      if(sgm.eq.0.d0) then
        xels = xels0
        sgm = xels0

      endif

      eprob(1) = xels/ sgm
      eprob(2) = xion/ sgm + eprob(1)
      eprob(3) = xexc/ sgm + eprob(2)

      sig_macro = (xels + xion + xexc) ! = sgm
      ec(ibkec+no,ipomp+1) = ec(ibkec+no,ipomp+1) * 1.d-6    ! from (eV) to (MeV)
      e(ibke+no,ipomp+1)   = e(ibke+ no,ipomp+1)  * 1.d-6    ! from (eV) to (MeV)
      end subroutine


c=================================================================
      subroutine etstrn_art(mark,markp,fpl,nbeta,itmak)
*
*
*       electron transfer
*       and region check
*       modified by K.Niita on 2003/10/12
*
*     input  :
*
*       fpl  : distance
*
*     output :
*
*       mark   : out put code of geom
*       markp  : =1 already check the cell
*       nbeta  : =1,2; reactions or cross 3; stopped
*       itmak  : 0, normal, 1, out of time range
*
c=================================================================
      use MMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'ggsparam.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

            if(e(ibke+no,ipomp+1).gt.etsmin)then
               nbeta = 2
            else
               nbeta = 3
            endif


      if(lev .gt. 0) then
          ll   = lev - 1    ! coordinates to measure the distance to the boundary
        xxx  = udt(1,ll)
        yyy  = udt(2,ll)
        zzz  = udt(3,ll)
        uuu  = udt(4,ll)
        vvv  = udt(5,ll)
        www  = udt(6,ll)
      else
        xxx = x(ibkx+no,ipomp+1)
        yyy = y(ibky+no,ipomp+1)
        zzz = z(ibkz+no,ipomp+1)
        uuu = u(ibku+no,ipomp+1)
        vvv = v(ibkv+no,ipomp+1)
        www = w(ibkw+no,ipomp+1)
      endif

      coincd = coincd * 1.d-4

*-----------------------------------------------------------------------
*        distance to the boundary
*-----------------------------------------------------------------------


               call gomdis(dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                         z(ibkz+no,ipomp+1),
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1),
     &                     mark,markp,nmed(ibknmd+no,ipomp+1),
     &                     iblz(ibkblz+no,ipomp+1))

               if( mark .le. -2 ) then
                   coincd = coincd * 1.d+4
                   return
               endif


*-----------------------------------------------------------------------
*        etectron transport
*-----------------------------------------------------------------------


        if ( fpl .ge. dpr ) then ! (boundary crossing)

               call gomnew(1,dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                     z(ibkz+no,ipomp+1),
     &                     xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                     zc(ibkzc+no,ipomp+1),
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1),
     &                     ec(ibkec+no,ipomp+1),
     &                     nmed(ibknmd+no,ipomp+1),
     &                     iblz(ibkblz+no,ipomp+1),
     &                     mark,markp)

               fpl = dpr

        else ! (normal)

               call gomupr(mark,markp,fpl)

        end if


*-----------------------------------------------------------------------
*        time evolution
*-----------------------------------------------------------------------

               call timtrs(itmak,mark)

*------------------ -----------------------------------------------------
               coincd = coincd * 1.d+4
      return
      end subroutine


c=================================================================
      subroutine etsreac_art
c=================================================================
      use MMBANKMOD
      implicit real*8(a-h,o-z)

      real*8 ran2
      dimension vel(3)
      dimension velp(3),vels(3)


         u1st = 0.d0
         v1st = 0.d0
         w1st = 1.d0
       change = 1.d+5        ! (100 keV)
       ehnum  = 0.d0

      ec(ibkec+no,ipomp+1) = ec(ibkec+no,ipomp+1) * 1.d+6   ! from (MeV) to (eV)
      e(ibke+no,ipomp+1)   = e(ibke+no,ipomp+1)   * 1.d+6   ! from (MeV) to (eV)
      ene = e(ibke+no,ipomp+1)

        if( .not. allocated(xsec_shel)) then
          if(emin(12).gt.etsmin) then
           write(*,'(a80)') 'Warning: etsmin is greater than emin(12).
     & egs is executed after Track Strucutre.'
          endif
          return
        endif
!--- end change NS 2020.04 del THREADPRIVATE
         icoll = 0
         iagr  = 0
         e2nd  = 0.d0
         id_reac = 0

! initialize relaxation number
         lng_rel = 0


c // determiation of e, u, v, w //
         prb   = dble(ran2(idum))
c // elastic scattering (id=1) //
      if(prb.gt.0.d0 .and. prb.le.eprob(1))then
         id_reac = 1      ! S.Abe 2018/02/07
   11 continue
         call els_dcs_art(xthe,xphi)
         call vel_vec(vel,xthe,xphi) !!! subroutines in etsmode.f !!!!

         u1st = vel(1)
         v1st = vel(2)
         w1st = vel(3)
         nclsts     = 1
         dexc_ene = 0.d0

      goto 1000
      endif

c // ionization (id=2) //
      if(prb.gt.eprob(1).and.prb.le.eprob(2))then
        id_reac = 2
        icoll   = 1
        hydro = denh_das(kmat0+mat)
        lem   = nint( dnel_das(kmat0+mat) )
        ssum = 0.d0
        pk = dble(ran2(idum))

        ii = 0
        if(hydro .gt. 0.d0) then ! for 1H
         ii = 1
         i  = 1
         j    = 1
         itz  = 1
         ssum = ssum + xsec_shel(j,i)/xion
         if(pk.le.ssum)   goto 21
        endif

        do i = 1 + ii, lem + ii ! other elements
         do j = 1, 28 ! subshell loop
           ssum = ssum + xsec_shel(j,i)/xion
           itz  = nint( zz_das(kmat(mat)+i-ii) )
           if(pk.le.ssum)   goto 21
         enddo
        enddo

  21    continue
! ion_dcs is different depend on energy
! low energy electron
           call ion_dcs_art(ene,itz,j,bebw2nd)
! 20221128 relativistic effect will be implemented later ....
           call ion_vec_art(itz,j,bebw2nd)
           if(itz.ge.6) call do_atom_relax(itz,j)

      goto 1000
      endif

c // excitation (id=3) //
      if(prb.gt.eprob(2).and.prb.le.eprob(3))then
        id_reac = 3
        icoll   = 1
        call exc_dcs_art(ene,bebw2nd)
        e2nd = bebw2nd
        dexc_ene = bg
        call ion_ang_art(bg,e2nd,pthe,pphi,sthe,sphi)
!!! subroutines in etsmode.f !!!!
        call vel_vec(velp,pthe,pphi)
        call vel_vec(vels,sthe,sphi)
        u1st =  velp(1)
        v1st =  velp(2)
        w1st =  velp(3)
        ec(ibkec+no,ipomp+1) =  e(ibke+no,ipomp+1) - dexc_ene - e2nd
        if(ec(ibkec+no,ipomp+1) .lt. 0.d0)then
! 20230619 epsilon value correction for calculating charged particles
         eleemit = 1.d0
         if(wvlin.gt.0.d0) then
           if(ec(ibkec+no,ipomp+1).gt.wvlin) then
            eleemit = dble(int(ec(ibkec+no,ipomp+1)/wvlin))
            if(eleemit.lt.1.d0) eleemit = 1.d0
           endif
           ehnum = ehnum + eleemit -1.d0 ! 20230802 interaction mark for 'ts_w'
           if(no.eq.1) ehnum = ehnum + 1.d0  ! source electron
         endif
         dexc_ene   =  bg
         e2nd       =  0.d0
! electron generated number

            ec(ibkec+no,ipomp+1) =  e(ibke+no,ipomp+1) - dexc_ene - e2nd
        endif
        u2nd  =  vels(1)
        v2nd  =  vels(2)
        w2nd  =  vels(3)
      endif
 1000 continue
!!!!!!!!! secondary electron killing for range calculation!!!!!!!!!!!!!!!!!!!
        if(msec.eq.1) then
          icoll = 0
          dexc_ene = dexc_ene + e2nd
          e2nd = 0.d0
          lng_rel= 0
        endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c// finalization and data up
! electron stop when energy less than etsmin
      if(ec(ibkec+no,ipomp+1) .le. etsmin*1.d6 )then ! electron absorption
         eleemit = 1.d0
! 20230619 epsilon value correction for calculating charged particles
         if(wvlin.gt.0.d0) then
           if(ec(ibkec+no,ipomp+1).gt.wvlin) then
            eleemit = dble(int(ec(ibkec+no,ipomp+1)/wvlin))
            if(eleemit.lt.1.d0) eleemit = 1.d0
           endif
           ehnum = ehnum + eleemit -1.d0 ! 20230802 interaction mark for 'ts_w'
           if(no.eq.1) ehnum = ehnum + 1.d0  ! source electron
         endif

         dexc_ene     = dexc_ene + ec(ibkec+no,ipomp+1)
! t-interact ets_e-exc
         ec(ibkec+no,ipomp+1) = 0.0d0


! electron stop when energy less than bandgap energy
      else if(eprob(1).eq.1.d0)then ! eprob(1) = 1.d0 means xion and xexc equal zero.
         eleemit = 1.d0
! 20230619 epsilon value correction for calculating charged particles
         if(wvlin.gt.0.d0) then
           if(ec(ibkec+no,ipomp+1).gt.wvlin) then
            eleemit = dble(int(ec(ibkec+no,ipomp+1)/wvlin))
            if(eleemit.lt.1.d0) eleemit = 1.d0
           endif
           ehnum = ehnum + eleemit  ! 20230802 interaction mark for electron counting
         endif
         dexc_ene     = dexc_ene + ec(ibkec+no,ipomp+1)
! t-interact ets_e-exc
         ec(ibkec+no,ipomp+1) = 0.0d0

      endif

      dexc_ene = dexc_ene * 1.d-6 ! eV to MeV for tallying

      id_reac2 = id_reac
      if(id_reac2.eq.3) id_reac2 = 2  ! excitation in ETSART emits an electron. This is defined as the ionization in T-interact.
      if( ityp .eq. 12 .and. id_reac.ne.0) then
         atmrc(6,id_reac2) = atmrc(6,id_reac2) + 1.d0
      elseif( ityp .eq. 13 .and. id_reac.ne.0 ) then
         atmrc(7,id_reac2) = atmrc(7,id_reac2) + 1.d0
      endif

*-----------------------------------------------------------------------
*        data up in bank
*-----------------------------------------------------------------------

      nclsts  =  0

! use a subroutine in etsmode02.f
      call etsdataup02 ! primary electron data up

c // Ionized electrons //
      numpal(12) = 0
      rumpal(12) = 0.d0
      numpal(14) = 0
      rumpal(14) = 0.d0

      wga        = wt(ibkwt+no,ipomp+1)

c // Ionized electrons //
      if(icoll.eq.1)then
! 2nd electron stop when energy less than etsmin
        if(e2nd .le. etsmin*1.d+6 )then ! electron absorption
         eleemit = 1.d0
! 20230619 epsilon value correction for calculating charged particles
         if(wvlin.gt.0.d0) then
           if(e2nd.gt.wvlin) then
            eleemit = dble(int(e2nd/wvlin))
            if(eleemit.lt.1.d0) eleemit = 1.d0
           endif
           ehnum = ehnum + eleemit -1.d0 ! 20230802 interaction mark for 'ts_w'
         endif
           dexc_ene     = dexc_ene + e2nd * 1.d-6
         goto 1002
        endif

         nclsts = nclsts + 1
         iclusts( nclsts ) =  7
         jclusts(  0, nclsts ) = 0
         jclusts(  1, nclsts ) = 0
         jclusts(  2, nclsts ) = 0
         jclusts(  3, nclsts ) = 12
         jclusts(  4, nclsts ) = 0
         jclusts(  5, nclsts ) = -1
         jclusts(  6, nclsts ) = 0
         jclusts(  7, nclsts ) = 11
         jclusts(  8, nclsts ) = 0

         if(e2nd.le.0.0) then
           pel = 0.0
         else
            pel = sqrt( (e2nd*1.d-9)**2 + 2.d0 * e2nd * 1.d-9 *
     &                rmtyp(ityp,ktyp) * 1.d-3 ) ! momentum in GeV
         endif

         qclusts(  0, nclsts ) = 0.d0                      ! impact parameter
         qclusts(  1, nclsts ) = pel * u2nd                ! (GeV/c)
         qclusts(  2, nclsts ) = pel * v2nd                ! (GeV/c)
         qclusts(  3, nclsts ) = pel * w2nd                ! (GeV/c)
         qclusts(  4, nclsts ) = e2nd * 1.d-9 + rmtyp(ityp,ktyp) * 1.d-3! (GeV)
         qclusts(  5, nclsts ) = rmtyp(ityp,ktyp) * 1.d-3  ! (GeV)
         qclusts(  6, nclsts ) = 0.d0                      ! from (eV) to (MeV)
         qclusts(  7, nclsts ) = e2nd * 1.d-6              ! from (eV) to (MeV)
         qclusts(  8, nclsts ) = wga                       ! weight
         qclusts(  9, nclsts ) = 0.d0                      ! (ns)
         qclusts( 10, nclsts ) = 0.d0                      ! (cm)
         qclusts( 11, nclsts ) = 0.d0                      ! (cm)
         qclusts( 12, nclsts ) = 0.d0                      ! (cm)

         numpal(12) = numpal(12) + 1
         rumpal(12) = rumpal(12) + wga
 1002    continue
      endif

! auger electorn dataup copied itsdataup
      if( lng_rel.ne.0) then
      do i = 1, lng_rel ! atomic relaxation
! 20211216 for excitation counting
       if(eng_rel(i).le.etsmin .and. ktp_rel(i).ne.22) then ! 20230105 error etsimin and eng_rel are MeV
         eleemit = 1.d0
! 20230619 epsilon value correction for calculating charged particles
         if(wvlin.gt.0.d0) then
           if(eng_rel(i)*1.d+6 .gt. wvlin) then
            eleemit = dble(int(eng_rel(i)*1.d+6/wvlin))
            if(eleemit.lt.1.d0) eleemit = 1.d0
           endif
           ehnum = ehnum + eleemit ! 20230802 interaction mark for 'ts_w'
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
! deposit energy modify
       dexc_ene = dexc_ene - eng_rel(i)  ! dexc_ene converted to MeV
      enddo
      lng_rel = 0
      endif


c // Pair production //
      if(ityp. eq. 13 .and. ec(ibkec+no,ipomp+1) .eq. etsmin*1.d+6 .and.
     .                      e(ibke+no,ipomp+1)   .gt. etsmin*1.d+6 )then

c // positron decay //
            bdcayp(13) = bdcayp(13) + 1.d0
            adcayp(13) = adcayp(13) + 1.d0

c // prod. photons //

            ranx = dble(ran2(idum))
            rany = dble(ran2(idum))
            ranz = dble(ran2(idum))

            phnu = -0.5D0 + ranx
            phnv = -0.5D0 + rany
            phnw = -0.5D0 + ranz
            phnn = dsqrt(phnu**2 + phnv**2 + phnw**2)

            uphoton = phnu / phnn
            vphoton = phnv / phnn
            wphoton = phnw / phnn

        do i=1,2
         nclsts = nclsts + 1

         iclusts( nclsts ) =  4
         jclusts(  0, nclsts ) = 0
         jclusts(  1, nclsts ) = 0
         jclusts(  2, nclsts ) = 0
         jclusts(  3, nclsts ) = 14
         jclusts(  4, nclsts ) = 0
         jclusts(  5, nclsts ) = 0
         jclusts(  6, nclsts ) = 0
         jclusts(  7, nclsts ) = 22
         jclusts(  8, nclsts ) = 0

         qclusts(  0, nclsts ) = 0.d0                      ! impact parameter
         qclusts(  1, nclsts ) = rmtyp(12,11) * 1.d-3 * uphoton *(-1)**i ! (GeV/c)
         qclusts(  2, nclsts ) = rmtyp(12,11) * 1.d-3 * vphoton *(-1)**i ! (GeV/c)
         qclusts(  3, nclsts ) = rmtyp(12,11) * 1.d-3 * wphoton *(-1)**i ! (GeV/c)
         qclusts(  4, nclsts ) = rmtyp(12,11) * 1.d-3      ! (GeV)
         qclusts(  5, nclsts ) = 0.d0                      ! (GeV)
         qclusts(  6, nclsts ) = 0.d0                      ! from (eV) to (MeV)
         qclusts(  7, nclsts ) = rmtyp(12,11)              ! from (eV) to (MeV)
         qclusts(  8, nclsts ) = wga                       ! weight
         qclusts(  9, nclsts ) = 0.d0                      ! (ns)
         qclusts( 10, nclsts ) = 0.d0                      ! (cm)
         qclusts( 11, nclsts ) = 0.d0                      ! (cm)
         qclusts( 12, nclsts ) = 0.d0                      ! (cm)

         numpal(14) = numpal(14) + 1
         rumpal(14) = rumpal(14) + wga
        enddo


      endif

      if(wvlin .gt. 0.d0 .and. ehnum.ne.0.d0) then
        id_reac2 = id_reac
        if(id_reac2.eq.3) id_reac2 = 2  ! excitation in ETSART emits an electron. This is defined as the ionization in T-interact.
        if( ityp .eq. 12 .and. id_reac.ne.0) then
           atmrc(6,id_reac2) = atmrc(6,id_reac2) + ehnum
        elseif( ityp .eq. 13 .and. id_reac.ne.0 ) then
           atmrc(7,id_reac2) = atmrc(7,id_reac2) + ehnum
        endif
      endif
!--- change NS 2020.04 del THREADPRIVATE
!      ec(ibkec+no) = ec(ibkec+no) * 1.d-6   ! from (eV) to (MeV)
!      e(ibke+no)   = e(ibke+ no)  * 1.d-6   ! from (eV) to (MeV)
      ec(ibkec+no,ipomp+1) = ec(ibkec+no,ipomp+1) * 1.d-6   ! from (eV) to (MeV)
      e(ibke+no,ipomp+1)   = e(ibke+ no,ipomp+1)  * 1.d-6   ! from (eV) to (MeV)
        if( .not. allocated(xsec_shel)) then
        else
            deallocate(xsec_shel)
        endif

      end subroutine


c================================================
      subroutine ion_vec_art(itz,ish,bebw2nd)
c===============================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)
      dimension velp(3),vels(3)


        BI  = pot_elem(ish,itz)
        dexc_ene   =  BI
        e2nd       =  bebw2nd
      call ion_ang_art(BI,e2nd,pthe,pphi,sthe,sphi)
!!! subroutines in etsmode.f !!!!
      call vel_vec(velp,pthe,pphi)
      call vel_vec(vels,sthe,sphi)

        u1st =  velp(1)
        v1st =  velp(2)
        w1st =  velp(3)

      ec(ibkec+no,ipomp+1) =  e(ibke+no,ipomp+1) - dexc_ene - e2nd
      if(ec(ibkec+no,ipomp+1) .lt. 0.d0)then
        dexc_ene   =  BI
        e2nd       =  0.d0
        ec(ibkec+no,ipomp+1) =  e(ibke+no,ipomp+1) - dexc_ene - e2nd
      endif

        u2nd  =  vels(1)
        v2nd  =  vels(2)
        w2nd  =  vels(3)

      end subroutine

c====================================================
      subroutine ion_ang_art(BI,e2nd,pthe,pphi,sthe,sphi)
c====================================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)

      real*8 ran2


      dmc2 = 511.003d+3
      t2   = e(ibke+no,ipomp+1) - e2nd - BI
      wts  = e2nd

      tp   = t2/2.d0/dmc2
      if( t2 .eq. 0.d0 ) then
      p1   = 0.d0
      else
      p1   = wts/t2
      endif
      if(p1.lt.0.d0)then  ! (Takeshi)
      p1   = 0.d0         ! (Takeshi)
      endif               ! (Takeshi)
      p2   = (1.d0 - p1) * tp + 1.d0
      p3   = dsqrt(p1/p2)
      if(t2.ge.100.d0)then
      if(p3.gt.1.d0)then
      pthe = pi * 0.5d0
      else

        if(e(ibke+no,ipomp+1).le.3.d+2)then
        pthe = dasin(p3)        * 1.5d0
        endif

        if(e(ibke+no,ipomp+1).gt.3.d+2 .and.
     .     e(ibke+no,ipomp+1).le.1.d+5)then
        pthe = dasin(p3)        * 3.0d0
        endif

        if(e(ibke+no,ipomp+1).gt.1.d+5 .and.
     .     e(ibke+no,ipomp+1).le.5.d+5)then
        pthe = dasin(p3)        * 2.5d0
        endif

        if(e(ibke+no,ipomp+1).gt.5.d+5)then
        pthe = dasin(p3)        * 1.0d0
        endif

      endif

      else
      pthe = dble(ran2(idum)) * pi / 4.d0
      endif

      wpts = wts/2.d0/dmc2
      s1   = 1.d0 - p1
      if(s1.lt.0.d0)then  ! (Takeshi)
      s1   = 0.d0         ! (Takeshi)
      endif               ! (Takeshi)
      s2   = 1.d0 + wpts
      s3   = dsqrt(s1/s2)
      if(wts.ge.200.d0)then
      sthe = dasin(s3)
      else
       if(wts.le.50.d0)then
       sthe = dble(ran2(idum)) * 1.d0 * pi
       else
       sthe = pi/4.d0 + dble(ran2(idum)) * pi/4.d0
       endif
      endif

      deg  = 180.d0/pi
      pphi = dble(ran2(idum)) * 2.d0 * pi
      if(pphi.gt.pi)then
      sphi = pphi - pi
      else
      sphi = pphi + pi
      endif


      end subroutine

c====================================================
      subroutine els_ics_art(ene)
c====================================================
      implicit real*8(a-h,o-z)

      xels  = 0.d0
      elsum = 0.d0
      sumZi = 0.d0
      sumden= 0.d0
      dm    = 2.94d0

! effective atomic number calculation
       hydro = denh_das(kmat0+mat)
      if(hydro .gt. 0.d0) then
          itz = 1
        elsum = elsum + (dble(itz)**dm)*hydro*dble(itz)
        sumZi = sumZi + hydro*dble(itz)
        sumden= sumden+ hydro
      endif
      lem   = nint( dnel_das(kmat0+mat) )
      do i = 1, lem
        itz  = nint( zz_das(kmat(mat)+i) )
        dens = den_das(kmat(mat)+i)
        elsum = elsum + (dble(itz)**dm)*dens*dble(itz) ! Zi**2.94*(Zi)
        sumZi = sumZi + dens*dble(itz)  ! ‡”Zi
        sumden= sumden+ dens
      enddo
      zp    = (elsum/sumZi)**(1.d0/dm)  ! 2.94ã (Zi**2.94*(Zi/‡”Zi=fi))

      dmc2 = 511.003d+3
      re   = 2.8179d-13
      tp   = ene/dmc2
      bt2  = 1.d0 - 1.d0/(1.d0 + tp)**2

      if(ene.le.5.d+4)then
      etc  = 1.198
      else
      etc  = 1.13 + 3.76 * zp**2 / 137.d0**2 / bt2
      endif

      eta  = etc * 1.7d-5* zp**(2.d0/3.d0) / tp / (tp + 2.d0)
      s1   = pi * re **2 * zp * (zp + 1.d0)
      s2   = (1.d0 - bt2) / bt2**2
      s3   = 1.d0 / eta /(eta + 1.d0)
      cs   = s1 * s2 * s3
      xels = cs*sumden*1.d+24 ! cs (cm2) * sumden(10^24atom/cm3)

      end subroutine


c====================================================
      subroutine els_dcs_art(elsthe,elsphi)
c====================================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)

      dimension dels(0:180),sd(0:180)
      real*8 ran2


      dmc2 = 511.003d+3
      re   = 2.8179d-13


      do iang = 0,180
      the = dble(iang) *pi/180.d0

      tp  = e(ibke+no,ipomp+1)/dmc2
      bt2 = 1.d0 - 1.d0/(1.d0 + tp)**2

      if(e(ibke+no,ipomp+1).le.5.d+4)then
      etc = 1.198
      else
      etc = 1.13 + 3.76 * zp**2 / 137.d0**2 / bt2
      endif
      eta = etc * 1.7d-5* zp**(2.d0/3.d0) / tp / (tp + 2.d0)
      s1  = re **2 * zp * (zp + 1.d0)
      s2  = (1.d0 - bt2) / bt2**2
      s3  = (1.d0 - dcos(the) + 2.d0*eta)**2
      dels(iang) = s1 * s2 / s3
      enddo

      ss = 0.d0
      do i=0,180
      ss = ss + dels(i)
      enddo

      sd(0) = dels(0)/ss
      do i=1,180
      sd(i) = sd(i-1) + dels(i)/ss
      enddo

      pk = dble(ran2(idum))
       if(pk.le.sd(0))   elsthe = 0
      do i=1,180
       if(pk.le.sd(i)) then

        if(e(ibke+no,ipomp+1).le.3.d+2)then
        elsthe = dble(i-ran2(idum)) * pi/180.d0
        endif

        if(e(ibke+no,ipomp+1).gt.3.d+2 .and.
     .     e(ibke+no,ipomp+1).le.1.d+5)then
        elsthe = dble(i-ran2(idum)) * pi/180.d0 * 2.0d0
        endif

        if(e(ibke+no,ipomp+1).gt.1.d+5)then
        elsthe = dble(i-ran2(idum)) * pi/180.d0 * 2.0d0
        endif

        exit
       endif
      enddo
      elsphi = dble(ran2(idum)) * 2.d0 * pi

      end subroutine

c====================================================
      subroutine ion_ics_art(ene)
c====================================================
      implicit real*8(a-h,o-z)

      i = 0
      ii = 0
      xion = 0.d0
      itzp= 0.d0

       hydro = denh_das(kmat0+mat)
      if(hydro .gt. 0.d0) then
        i = 1
        ii = 1
        itzp = 1
        xsec_shel(1,i) = BEBTICS(ene,1,1) * hydro*1.d+24 ! Contribution of 1H
        xion = xion + xsec_shel(1,i)
      endif
      lem   = nint( dnel_das(kmat0+mat) )
      do i = 1 + ii, lem + ii
        itz  = nint( zz_das(kmat(mat)+i-ii) )
        dens = den_das(kmat(mat)+i-ii)
       do j = 1, 28 ! subshell loop
        TICS = BEBTICS(ene,itz,j)*dens*1.d+24  ! BEBTICS (cm2) * dens(10^24atom/cm3)
        xsec_shel(j,i) = xsec_shel(j,i) + TICS
        xion            = xion            + TICS
       enddo
      enddo

      end subroutine

c====================================================
      function BEBTICS(ene,itz,ish)
! Total inelastic cross section(cm2) by binary encounter bethe
! Phys. Rev. A vol. 50(5), pp. 3954-3967 : eq.(57)
c====================================================
      implicit double precision (a-h,o-z)

      parameter (a0=5.29177249d-9)
      parameter(nqmax=10000)  ! T.Sato 2022/08/22
      dimension sdion(nqmax),sd(nqmax)  ! T.Sato 2022/08/22

      BEBTICS = 0.d0

      ts = ene
      bebw2nd=0.d0 ! initialization  ! T.Sato 2022/08/22
      alpha = 1.d0/137.035999679d0
      dmc2  = 511.003d+3
       qt  = 1.d0

! itrn is used in ion_dcs in etsmode.f
      itrn = itz

      sdion(:)=0
      b   = pot_elem(ish,itz)
      uts = dkin_elem(ish,itz)
      ds  = elec_elem(ish,itz)
      if(b.gt.ene) return
      if(b.lt.bg) return  ! avoiding the  binding energy smaller than Band gap
      if(b.eq.0) then
       BEBTICS = 0.d0
       return
      endif

      st  = ts/b; tp = ts/dmc2; bp = b/dmc2; up = uts/dmc2
      bt2 = 1.d0 - 1.d0/(1.d0 + tp)**2
      bb2 = 1.d0 - 1.d0/(1.d0 + bp)**2
      bu2 = 1.d0 - 1.d0/(1.d0 + up)**2

      Q  = 1.d0

      s1 = 4.d0 * pi * a0**2 * alpha**4 * ds
      s2 = (bt2 + bb2 + bu2) * 2.d0 * bp
      s  = s1 / s2
      f1 = 0.5d0*(dlog(bt2/(1.d0-bt2)) - bt2 - dlog(2.d0*bp))
     &    *(1.d0-1.d0/st**2.0)
      f2 = 1.d0-1.d0/st-dlog(st)/(st+1.d0)*(1.d0+2.d0*tp)
     &     /(1.d0+tp/2.d0)**2.d0
      f3 = bp**2.d0/(1.d0+tp/2.d0)**2.d0*(st-1.d0)/2.d0
      BEBTICS = s*(f1+f2+f3)
      end function


c====================================================
      subroutine exc_ics_art(ene)
c====================================================
      implicit real*8(a-h,o-z)
      dimension excpro(nexc)
      real*8 ran2
      dimension Fl(3)
      dimension sexc(2)

      dmc2  = 511.003d+3

      excpro = 0.d0
      sexc = 0.d0
      iex2   = 0
      isp    = 0
      R       = 13.6d0          ! Rydberg Binding Energy
      a0      = 5.29d-11*1.d+2  ! Bohr orbit(cm)
      alpha = 1.d0/137.035999679d0
      xexc = 0.d0
      Q  = 1.d0
      i = 1
      do i =1,nexc + 1
       if(i.ne.nexc + 1) then ! single electron excitation
         BI  = bexc(i)
         Ul  = uexc(i)
         dN  = dnex(i)
         if(BI.gt.ene) cycle
         if(BI.eq.0.d0) cycle
       else             ! plasmon excitation
         BI  = eplasm
         Ul  = ULave
         dN  = sumvale
         if(BI.gt.ene) cycle
       endif

       TICS = 0.d0

! 20230427 Inter band transition using BEB model
        if(i.ne.nexc+1) then
          if(ene.ge.eplasm) goto 1000
          ! when ene>eplasm, then single excitation is prohibited. In this case, TICS=0.d0, when i!=nexc+1.
        endif
! BEB models for excitation
         if(ene.lt.BI) goto 1000
! same formula as BEBTICS function
          b   = BI
          uts = Ul
          ds  = dN
          ts  = ene
          st  = ts/b; tp = ts/dmc2; bp = b/dmc2; up = uts/dmc2
          bt2 = 1.d0 - 1.d0/(1.d0 + tp)**2
          bb2 = 1.d0 - 1.d0/(1.d0 + bp)**2
          bu2 = 1.d0 - 1.d0/(1.d0 + up)**2

          Q  = 1.d0

          s1 = 4.d0 * pi * a0**2 * alpha**4 * ds
          s2 = (bt2 + bb2 + bu2) * 2.d0 * bp
          s  = s1 / s2
          f1 = 0.5d0*(dlog(bt2/(1.d0-bt2)) - bt2 - dlog(2.d0*bp))
     &        *(1.d0-1.d0/st**2.0)
          f2 = 1.d0-1.d0/st-dlog(st)/(st+1.d0)*(1.d0+2.d0*tp)
     &         /(1.d0+tp/2.d0)**2.d0
          f3 = bp**2.d0/(1.d0+tp/2.d0)**2.d0*(st-1.d0)/2.d0
          TICS = s*(f1+f2+f3)
 1000 continue


        if(i.ne.nexc + 1) then
         hydro = denh_das(kmat0+mat)
         if(hydro .gt. 0.d0 .and. i .eq. 1) then
             xexc = xexc + TICS*hydro*1.d+24
             excpro(i) = excpro(i) + TICS*hydro*1.d+24
             sexc(1) = sexc(1) + TICS*hydro*1.d+24
         else
             xexc = xexc + TICS*denex(i)*1.d+24
             excpro(i) = excpro(i) + TICS*denex(i)*1.d+24
             sexc(1) = sexc(1) + TICS*denex(i)*1.d+24
         endif
        else
         xexc = xexc + TICS*sumdens*1.d+24
         sexc(2) = sexc(2) + TICS*sumdens*1.d+24
        endif


      enddo  ! material roop
      if(xexc.eq.0.d0) goto 2000

! determine single electron excitation or plasmon
      pk = dble(ran2(idum))
      if(pk.le.sexc(1)/xexc) then
        isp = 1
      else
        isp = 2
      endif
      if(isp.eq.2) goto 2000

      pk = dble(ran2(idum))
      sumex = 0.d0
      do i = 1,nexc
         sumex = sumex + excpro(i)
         if(pk.le.sumex/xexc) then
          iex2 = i
          goto 2000
         endif
      enddo
 2000 continue
      end subroutine

c====================================================
      subroutine ion_dcs_art(ene,itz,ish,bebw2nd)
c====================================================
      use MMBANKMOD
      implicit real*8(a-h,o-z)

      real*8 ran2

      parameter (a0=5.29177249d-9)
      parameter(nqmax=10000)  ! T.Sato 2022/08/22
      dimension sdion(nqmax),sd(nqmax)  ! T.Sato 2022/08/22

      ts = ene
      bebw2nd=0.0 ! initialization  ! T.Sato 2022/08/22
      alpha = 1.d0/137.035999679d0
      dmc2  = 511.003d+3
       qt  = 1.d0

! itrn is used in ion_dcs in etsmode.f
      itrn = itz

      sdion(:)=0
      b   = pot_elem(ish,itz)
      uts = dkin_elem(ish,itz)
      ds  = elec_elem(ish,itz)
      if(ts.lt.b) return ! T.Kai 2022/08/23

      do 2000 ip=1,nqmax-1  ! T.Sato 2022/08/22
      wi  = 0.d0
      wf  = min(100.0d3,(ts-b)/2.d0) ! T.Sato 2022/08/25, set maximum energy to 100 keV (adjustable)
      hq  = (wf-wi)/dble(nqmax)  ! T.Sato 2022/08/22
      wp  = wi + dble(ip) * hq
      sw  = wp/b
      st  = ts/b; tp = ts/dmc2; bp = b/dmc2; up = uts/dmc2
      bt2 = 1.d0 - 1.d0/(1.d0 + tp)**2
      bb2 = 1.d0 - 1.d0/(1.d0 + bp)**2
      bu2 = 1.d0 - 1.d0/(1.d0 + up)**2

      s1 = 4.d0 * pi * a0**2 * alpha**4 * ds
      s2 = (bt2 + bb2 + bu2) * 2.d0 * bp
      s  = s1 / s2
      w1 = 1.d0/( sw + 1.d0)
      w2 = 1.d0/( st - sw  )
      t1 = 1.d0 + 2.d0 * tp
      t2 = 1.d0 + tp/2.d0
      f1 = (qt - 2.d0)/(st + 1.d0) * (w1 + w2) * t1 / t2**2
      f2 = (2.d0 - qt) * (w1**2 + w2**2 + bp**2/t2**2)
      f3 = qt * ( w1**3 + w2**3 )
      f4 = dlog(bt2/(1.d0-bt2)) - bt2 - dlog(2.d0*bp)
      sdion(ip) = s *( f1 + f2 + f3*f4 )/b
      if(ip.ge.10.and.sdion(ip).lt.sdion(1)*1.0d-5) exit ! not necessary to consider any more

 2000 continue
      nq=ip  ! T.Sato 2022/08/22

      ss = 0.d0
      do i=1,nq
      ss = ss + sdion(i)
      enddo

      sd(1) = sdion(1)/ss
      do i=2,nq
      sd(i) = sd(i-1) + sdion(i)/ss
      enddo

      pk = dble(ran2(idum))
      pk2 = dble(ran2(idum))                                ! T.Kai 2022/08/23
      if(pk.le.sd(1)) then
       bebw2nd = hq*pk2*0.5d0
      else
       do i=1,nq-1
        if(pk.le.sd(i+1)) then
         bebw2nd = dble(i) * hq + hq*(pk2-0.5d0)                ! (old) bebw2nd = dble(i+1) * hq (T.Kai 2022/08/23)
         exit
        endif
       enddo
      endif

      end subroutine
c====================================================
      subroutine exc_dcs_art(ene,bebw2nd)
c====================================================
      use MMBANKMOD
      implicit real*8(a-h,o-z)
      real*8 ran2
      real*8,allocatable :: dsdw(:)
      dimension Fl(3)

      R       = 13.6d0          ! Rydberg Binding Energy
      a0      = 5.29d-11*1.d+2  ! Bohr orbit(cm)
      alpha = 1.d0/137.035999679d0
      dmc2  = 511.003d+3

      if(isp.eq.1) then
       BI  = bexc(iex2) ! band gap for Si
       Ul  = uexc(iex2) ! outer shell
       dN  = dnex(iex2)
      else
       BI  = eplasm
       Ul  = ULave
       dN  = sumvale
      endif
      tp  = ene/BI
      if(tp.lt.1.d0) then
        write(*,*) 'error to determine the excitation shell in etsart.'
      endif
      up  = Ul/BI
      S  = 4.d0*pi*a0**2.d0*dN*(R/BI)**2.d0
      Q  = 1.d0

      wste  = 0.1d0
      nw =  1000000
      allocate(dsdw(nw))

      dsdw = 0.d0
      Fl =0.d0
      adsdw = 0.d0

      do iw=1,nw
       wl = dble(iw)*wste
       wp  = wl/BI
       if((wp+1.d0).gt.tp)   goto 1000
       wmax = (tp-1.d0)/2.d0
       if(wp.gt.wmax)  goto 1000         !loop exit
       if(wp.gt.bimin)  goto 1000         !loop exit
! Phys. Rev. A vol. 50(5), pp. 3954-3967 : eq.(54)
       Fl(2) = (2.d0-Q)/(tp+up+1.d0)
       Fl(1) = - Fl(2)/(tp+1.d0)
       Fl(3) = Q*log(tp)/(tp+up+1.d0)
       do i = 1,3
         fw = 1.d0/(wp+1.d0)**dble(i)
         ftw= 1.d0/(tp-wp)**dble(i)
! Phys. Rev. A vol. 50(5), pp. 3954-3967 : eq.(52) [single differential cross section]
         dsdw(iw) = dsdw(iw) + S*Fl(i)*(fw+ftw)/BI
       enddo
       adsdw = adsdw + dsdw(iw)
       enddo
 1000 continue
      pk = dble(ran2(idum))
      sd = 0.d0
      do iw=1,nw
       sd = sd + dsdw(iw)/adsdw
       if(pk.le.sd) goto 2000
      enddo
 2000 continue
      bebw2nd = dble(iw)*wste

      if(bebw2nd+bg.gt.ene) then
          bebw2nd = ene - bg
      endif

      deallocate(dsdw)
      end subroutine

c====================================================
      subroutine exc_dcs_art_Relative_model(ene,bebw2nd) ! currently not use
c====================================================
      use MMBANKMOD
      implicit real*8(a-h,o-z)
      real*8 ran2

      parameter (a0=5.29177249d-9)
      parameter(nqmax=10000)  ! T.Sato 2022/08/22
      dimension sdion(nqmax),sd(nqmax)  ! T.Sato 2022/08/22

       ts = ene
      bebw2nd=0.0 ! initialization  ! T.Sato 2022/08/22
      alpha = 1.d0/137.035999679d0
      dmc2  = 511.003d+3
       qt  = 1.d0

      sdion(:)=0

      if(isp.eq.1) then
       b  = bexc(iex2) ! band gap for Si
       uts  = uexc(iex2) ! outer shell
       ds  = dnex(iex2)
      else
       b  = eplasm
       uts  = ULave
       ds  = sumvale
      endif


      if(ts.lt.b) return ! T.Kai 2022/08/23

      do 2000 ip=1,nqmax-1  ! T.Sato 2022/08/22
      wi  = 0.d0
      wf  = min(100.0d3,(ts-b)/2.d0) ! T.Sato 2022/08/25, set maximum energy to 100 keV (adjustable)
      hq  = (wf-wi)/dble(nqmax)  ! T.Sato 2022/08/22
      wp  = wi + dble(ip) * hq
      sw  = wp/b
      st  = ts/b; tp = ts/dmc2; bp = b/dmc2; up = uts/dmc2
      bt2 = 1.d0 - 1.d0/(1.d0 + tp)**2
      bb2 = 1.d0 - 1.d0/(1.d0 + bp)**2
      bu2 = 1.d0 - 1.d0/(1.d0 + up)**2

      s1 = 4.d0 * pi * a0**2 * alpha**4 * ds
      s2 = (bt2 + bb2 + bu2) * 2.d0 * bp
      s  = s1 / s2
      w1 = 1.d0/( sw + 1.d0)
      w2 = 1.d0/( st - sw  )
      t1 = 1.d0 + 2.d0 * tp
      t2 = 1.d0 + tp/2.d0
      f1 = (qt - 2.d0)/(st + 1.d0) * (w1 + w2) * t1 / t2**2
      f2 = (2.d0 - qt) * (w1**2 + w2**2 + bp**2/t2**2)
      f3 = qt * ( w1**3 + w2**3 )
      f4 = dlog(bt2/(1.d0-bt2)) - bt2 - dlog(2.d0*bp)
      sdion(ip) = s *( f1 + f2 + f3*f4 )/b
      if(ip.ge.10.and.sdion(ip).lt.sdion(1)*1.0d-5) exit ! not necessary to consider any more

 2000 continue
      nq=ip  ! T.Sato 2022/08/22

      ss = 0.d0
      do i=1,nq
      ss = ss + sdion(i)
      enddo

      sd(1) = sdion(1)/ss
      do i=2,nq
      sd(i) = sd(i-1) + sdion(i)/ss
      enddo

      pk = dble(ran2(idum))
      pk2 = dble(ran2(idum))                                ! T.Kai 2022/08/23
      if(pk.le.sd(1)) then
       bebw2nd = hq*pk2*0.5d0
      else
       do i=1,nq-1
        if(pk.le.sd(i+1)) then
         bebw2nd = dble(i) * hq + hq*(pk2-0.5d0)                ! (old) bebw2nd = dble(i+1) * hq (T.Kai 2022/08/23)
         exit
        endif
       enddo
      endif
      if(bebw2nd+bg.gt.ene) then
          bebw2nd = ene - bg
      endif

      end subroutine
c================================================
      subroutine etsart_db(ierrdb)  ! for debug parameter reading
c================================================
      implicit real*8(a-h,o-z)
      logical   exex

      iblz000 = 0

! initialization
      bexc   = 0.d0
      uexc   = 0.d0
      dnex   = 0.d0
      blow   = 0.d0
      denex  = 0.d0
      eplasm = 0.d0

      msec = 0

! input for secondary electrons and elastic scattering treatment
!      inquire( file = 'secondon.txt', exist = exex )
!      if( exex .eqv. .false. ) then
!       write(*,'(a30)') ' =======  etsart.f ========== '
!       write(*,'(a30)') ' secondon.txt can not be open.'
!       write(*,'(a79)') ' msec = 0 : secnodary electron will generate
!     & and elastic scattering will occur.'
!       write(*,'(a30)') ' ============================ '
!       goto 100
!      endif
!      open(25,file='secondon.txt', status='old')
!        read(25,*) msec
!      close(25)
!      if(msec.eq. 1) then 
!       write(*,*) 'secondary electron will be killed'
!      endif
!      if(msec.eq. 2) then 
!       write(*,*) 'elastic scattering will be ignored'
!      endif
100   continue
      end subroutine

c================================================
      subroutine exc_setup
c================================================
      implicit real*8(a-h,o-z)

      echarg = 4.80320467330d-10  ! electron charge (esu)
      elmass = 9.1093829d-28      ! electron mass (g)
      hbar   = 6.582119569d-16    ! plank constant hbar (eVs)

      hydro = denh_das(kmat0+mat)
      lem   = nint( dnel_das(kmat0+mat) )

      sumdens = 0.d0
      sumvale = 0.d0
      ULave   = 0.d0
      bimin   = 1.d+10
      itz     = 0
      bg  = ebgets( idgr(iblz(ibkblz+no,ipomp+1)) )
      wvlin = ewvets( idgr(iblz(ibkblz+no,ipomp+1)) )
      
      ii = 0
      if(hydro .gt. 0.d0) then
       i = 1
       ii = 1
       mexel(i)  = 1
       bexc(i)   = bg
       uexc(i)   = u_outer(mexel(i))
       ULave     = ULave + u_outer(mexel(i))
       dnex(i)   = valence_ne(mexel(i))
       denex(i)  = hydro
       sumdens   = sumdens + denex(i)
       sumvale   = sumvale + dnex(i)*denex(i)
       if(bimin.gt.bin_outer_low(mexel(i))) then
          bimin  = bin_outer_low(mexel(i))
       endif

      endif
      do i = 1 + ii, lem + ii
       itz = nint( zz_das(kmat(mat)+i-ii))
       mexel(i)  = itz
       bexc(i)   = bg
       uexc(i)   = u_outer(mexel(i))
       ULave     = ULave + u_outer(mexel(i))
       dnex(i)   = valence_ne(mexel(i))
       denex(i)  = den_das(kmat(mat)+i-ii)

       sumdens   = sumdens + denex(i)
       sumvale   = sumvale + dnex(i)*denex(i)

       if(bimin.gt.bin_outer_low(mexel(i))) then
          bimin  = bin_outer_low(mexel(i))
       endif

      enddo
      nexc  = lem + ii

      sumvale = sumvale/sumdens

! plasmon energy and valance electrons
      ULave    = ULave/dble(nexc)
      eledens   = sumdens*sumvale*1.d+24
      eplasm  = hbar*(4.d0*pi*eledens*echarg**2.d0
     &               /elmass)**(0.5d0)
!      write(*,'(a16,f8.2,a3)') 'plasmon energy ',eplasm,' eV'
!      write(*,'(a16,f8.2,a3)') 'average energy ',ULave,' eV'
!      write(*,'(a16,f8.2,a3)') 'BandGap energy ',bg,' eV'

      end subroutine
c================================================
      function u_outer(itz)
c================================================
      implicit real*8(a-h,o-z)
      integer i
      u_outer = 0.d0

      do i = 1,28
       dN  = elec_elem(i,itz)
       if(dN.ne.0.d0) u_outer = dkin_elem(i,itz)
      enddo

      end function

c================================================
      function bin_outer_low(itz)
c================================================
      implicit real*8(a-h,o-z)
      integer i

      bin_outer_low = 0.d0
      peri = ele_period(itz)
      bimin0 = 1.d+10
      do i = 1,28
       dN  = elec_elem(i,itz)
       if(peri.eq.shel_priod(i)) then
        if(bimin0.gt.pot_elem(i,itz).and.pot_elem(i,itz).gt.bg) then ! if binding energy lower than band gap we ignore that process.
          if(dN.ne.0.d0) bimin = pot_elem(i,itz)
        endif
       endif
      enddo
      bin_outer_low = bimin0

      end function

c================================================
      function valence_ne(itz)
c================================================
      implicit real*8(a-h,o-z)
      integer i
       valence_ne = 0
       outele = ele_period(itz) ! outer shell
       do i = 1, 28
         if(outele .eq. shel_priod(i)) then
           valence_ne = valence_ne + elec_elem(i,itz)
         endif
       enddo
      end function
      end module
