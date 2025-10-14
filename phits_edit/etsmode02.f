c=================================================================
c
c  // Monte Carlo Track Stracture Code for low energy electron //
c
c                                          by Takeshi Kai in JAEA
c  for Silicon  by Yuho Hirata
c
c=================================================================
      subroutine etsflt02(sig_macro)
c=================================================================
      use MMBANKMOD
      use usrtalmod
      implicit real*8(a-h,o-z)

      include 'param00.inc'
      include 'param.inc'

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)
      common / etsprob02 / eprob(12),xion(10),xexc(8),prb2,prb3
!$OMP THREADPRIVATE(/etsprob02/)
      common /celepcc/ enumpcc(kvlmax)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)



      dimension vel(3),ss(0:10)

      common / cschng / change, fit
      common / secondon/ m_second
      dimension dnel(1), denh(1), zz(1), a(1), den(1)
      equivalence ( das, dnel, denh, zz, a, den )
      common /kmat1g/ kmat(kvlmax)

      change = 1.d+5        ! (100 keV)
      fit = 2.8851234d-18 / 2.4942487d-18  ! (100 keV) (Liq/Gas)


      ec(ibkec+no,ipomp+1) = ec(ibkec+no,ipomp+1) * 1.d+6   ! from (MeV) to (eV)
      e(ibke+no,ipomp+1)   = e(ibke+no,ipomp+1)   * 1.d+6   ! from (MeV) to (eV)

      ene = e(ibke+no,ipomp+1)
c // cross section (cm2) //

      call  els_ics02(ene,xels)
      call  ion_ics02(ene,xion)
      call  exc_ics02(ene,xexc)
      call  phn_ics02(ene,xphn)
        prb1  = xels

!!!!!! elastic off for e2nd and excitation check !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Si ionization mode => 6
      prb2  = 0.d0
      do ip=1,6
        prb2  = prb2 + xion(ip)
      enddo
! Si excitation
        prb3  = 0.d0
      do ip=1,1
        prb3  = prb3+ xexc(ip)
      enddo
! Si phonon energy loss
        prb4  = xphn

c // total cross section (sgm) //
       sgm      = prb1  +  prb2 + prb3 + prb4
        eprob(1)   =   prb1    /   sgm
        eprob(2)   =   prb2    /   sgm   +   eprob(1)
        eprob(3)   =   prb3    /   sgm   +   eprob(2)
        eprob(4)   =   prb4    /   sgm   +   eprob(3)

c // determination of mean free path (xlm) and flight length (fpl) //
        wanumpcc  = 3.318565377871046d+022
        wenumpcc  = 3.323785188937734d+023
        enumscale = enumpcc(mat)/wenumpcc

! new version 20210915
        sig_macro= sgm * (wanumpcc*enumscale)


      ec(ibkec+no,ipomp+1) = ec(ibkec+no,ipomp+1) * 1.d-6    ! from (eV) to (MeV)
      e(ibke+no,ipomp+1)   = e(ibke+ no,ipomp+1)  * 1.d-6    ! from (eV) to (MeV)

      end


************************************************************************
*                                                                      *
      subroutine etstrn02(mark,markp,fpl,nbeta,itmak)
*                                                                      *
*                                                                      *
*       electron transfer                                              *
*       and region check                                               *
*       modified by K.Niita on 2003/10/12                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       fpl  : distance                                                *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       mark   : out put code of geom                                  *
*       markp  : =1 already check the cell                             *
*       nbeta  : =1,2; reactions or cross 3; stopped                   *
*       itmak  : 0, normal, 1, out of time range                       *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      common / jcomon / nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common / icomon / no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)
      common / tlgeom / iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common / clustt / nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common / clustw / jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common / etsminmax / etsmin, etsmax

*-----------------------------------------------------------------------


            if(e(ibke+no,ipomp+1).gt.etsmin)then
               nbeta = 2
            else
               nbeta = 3
            endif



*-----------------------------------------------------------------------
*        distance to the boundary
*-----------------------------------------------------------------------


               call gomdis(dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                         z(ibkz+no,ipomp+1),
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1),
     &                     mark,markp,nmed(ibknmd+no,ipomp+1),
     &                     iblz(ibkblz+no,ipomp+1))

               if( mark .le. -2 ) return


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

*-----------------------------------------------------------------------

      return
      end


c=================================================================
      subroutine etsreac02
c=================================================================
      use ELEDATAMOD
      use MMBANKMOD
      implicit real*8(a-h,o-z)

      include 'param00.inc'
      include 'param-physcnst.inc'
      include 'param.inc'

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)
      common / clustt / nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE( /clustf/)
      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)

      common /poabs/  iabsms
!$OMP THREADPRIVATE(/poabs/)
      common /dcayp/  adcayp(20), bdcayp(20)

      common / etsprob02 / eprob(12),xion(10),xexc(8),prb2,prb3
!$OMP THREADPRIVATE(/etsprob02/)
      common / ele2nd  / e2nd,u2nd,v2nd,w2nd
!$OMP THREADPRIVATE(/ele2nd/)
      common / etsexe  / dexc_ene
!$OMP THREADPRIVATE(/etsexe/)
      common / ele1st / u1st,v1st,w1st
!$OMP THREADPRIVATE(/ele1st/)

      common / cschng / change, fit

      common / etsminmax / etsmin, etsmax
      common / csion02 / dion(500),si(6,500),xip(6)
      common / cslow / elow(100),se(1,100),sp(1,100)
      common / dcsexc / deexc(500),de(1,500),deprob(1,500),deall
      common / dcsphn / dephn(100),dp(1,100),dpprob(1,500),dpall
      common / heexc / ea1(8)
      common / idumc /idum
      common / secondon/ m_second
      common / Plasmonemit/ nplasm

      real*8 ran2

      dimension vel(3),ss(0:10), eagr(100)

      common /tstara/ atmrc(10,8)
!$OMP THREADPRIVATE(/tstara/)
      common /etsion/ icoll
!$OMP THREADPRIVATE(/etsion/)

      common /etsart/ bgets(kvlmax), ebgets(kvlmax),
     &                wvets(kvlmax), ewvets(kvlmax)
      common /ets_wvalue/ wvlin
!$OMP THREADPRIVATE(/ets_wvalue/)
         u1st = 0.d0
         v1st = 0.d0
         w1st = 1.d0
      ec(ibkec+no,ipomp+1) = ec(ibkec+no,ipomp+1) * 1.d+6   ! from (MeV) to (eV)
      e(ibke+no,ipomp+1)   = e(ibke+no,ipomp+1)   * 1.d+6   ! from (MeV) to (eV)

         icoll = 0
         iagr  = 0
         e2nd  = 0.d0
         id_reac = 0
         ehnum = 0.d0
         wvlin = ewvets( idgr(iblz(ibkblz+no,ipomp+1)) )

! initialize relaxation number
         lng_rel = 0

c // determiation of e, u, v, w //
         prb   = dble(ran2(idum))
c // elastic scattering (id=1) //
      if(prb.gt.0.d0.and.prb.le.eprob(1))then
         id_reac = 1      ! S.Abe 2018/02/07
   11 continue
         call els_dcs02(xthe,xphi)
         call vel_vec02(vel,xthe,xphi)
         u1st = vel(1)
         v1st = vel(2)
         w1st = vel(3)
         nclsts     = 1
         dexc_ene = 0.d0
!20220413 elastic energy loss by momentum transfer

      goto 1000
      endif


c // ionization (id=2) //
      if(prb.gt.eprob(1).and.prb.le.eprob(2))then
         id_reac = 2      ! S.Abe 2018/02/07
         ss(0)  = 0.d0
! plasmon is not ionization
         ss(1)  = ss(0) + xion(1)/prb2
         do ip=2,6
           if(e(ibke+no,ipomp+1).ge.xip(ip))then
           ss(ip) = ss(ip-1) + xion(ip)/prb2
           else
           ss(ip) = ss(ip-1) + 0.d0
           endif
         enddo
         pk = dble(ran2(idum))
         if(pk.gt.0.d0 .and.pk.le.ss(1))   i = 1
         do j=1,5
         if(pk.gt.ss(j).and.pk.le.ss(j+1)) i = j + 1
         enddo
         if(pk.gt.ss(6)) then !ss(6) is sometimes less than 1.0 for rounding error
           do j=0,5
             if(xion(6-j).ne.0.0) then
                i = 6-j
                exit
             endif
           enddo
          endif
! ion_dcs is different depend on energy
! low energy electron
         if( e(ibke+no,ipomp+1) .le. change) then
           call ion_dcs02(i,bebw2nd)
! high energy electron
         else
           if(i.eq.1) then
             ih = 2
           else
             ih = i + 1  ! for ID change for ion_dcs0
           endif
           call ion_dcs0(ih,bebw2nd)
         endif

         call ion_vec02(i,xip,bebw2nd)

         nclsts = 1
         icoll = 1

          if(i.ne.1)then
!! corresponding shell and ID
!    Shell   EADL(ld)   MyOELF(i)
!      K        1          6
!      L1       3          5
!      L2       5          4
!      L3       6          3
!      M1       8          2
! 20220401 Modified subshell ID
              if(i.eq.6) ld = 1  ! K shell
              if(i.eq.5) ld = 3  ! L1 shell
              if(i.eq.4) ld = 5  ! L2 shell
              if(i.eq.3) ld = 6   ! L3 shell
              if(i.eq.2) ld = 8  ! M1 shell (no relaxation)
              call do_atom_relax(14,ld)  ! 14 is atomic number of Si
          endif
       goto 1000
      endif


c// electronic excitation (id=3) //
      if(prb.gt.eprob(2).and.prb.le.eprob(3))then
        id_reac = 3      ! S.Abe 2018/02/07

        call exc_dcs02(e(ibke+no,ipomp+1),1,ea1(1))
! excitation generate a new electron
       if(ea1(1).gt.etsmin * 1.d+6) then
        icoll = 1
        dexc_ene = 1.11  ! band gap energy deposition as a potential
        e2nd = ea1(1) - dexc_ene
        ranx = dble(ran2(idum))
        rany = dble(ran2(idum))
        ranz = dble(ran2(idum))
        vx   = -0.5d0 + ranx
        vy   = -0.5d0 + rany
        vz   = -0.5d0 + ranz
        vv   = dsqrt(vx**2 + vy**2 + vz**2)
        u2nd = vx / vv
        v2nd = vy / vv
        w2nd = vz / vv
c! excitation electron deposits energy at generated position
       else if(ea1(1).le.etsmin * 1.d+6) then
         icoll = 0
         dexc_ene = ea1(1)
         e2nd = 0.0d0

         eleemit = 1.d0
! 20230619 epsilon value correction for calculating charged particles
         if(wvlin.gt.0.d0) then
           if(ea1(1).gt.wvlin) then
            eleemit = dble(int(ea1(1)/wvlin))
            if(eleemit.lt.1.d0) eleemit = 1.d0
           endif
           ehnum = ehnum + eleemit -1.d0
         endif
       endif
        ec(ibkec+no,ipomp+1) = e(ibke+no,ipomp+1) - ea1(1)
! if electron energy was under the plasmon peak, primary electron change the direction
        if(e(ibke+no,ipomp+1).lt.16.7) then
           ranx = dble(ran2(idum))
           rany = dble(ran2(idum))
           ranz = dble(ran2(idum))
           vx   = -0.5d0 + ranx
           vy   = -0.5d0 + rany
           vz   = -0.5d0 + ranz
           vv   = dsqrt(vx**2 + vy**2 + vz**2)
           u1st = vx / vv
           v1st = vy / vv
           w1st = vz / vv
        endif

      goto 1000
      endif

cc // phonon excitation (id=4) //
      if(prb.gt.eprob(3).and.prb.le.eprob(4))then
         id_reac = 6      ! hirata 20230719 changed to same ID as etsmode.f
         icoll = 0
         nclsts       = 1
! phonon interaction will tak at initial positon
         call phn_dcs02(e(ibke+no,ipomp+1),1,pip)
         ec(ibkec+no,ipomp+1) = e(ibke+no,ipomp+1) - pip
         dexc_ene     = pip

       goto 1000
      endif


c// finalization and data up
 1000 continue
! electron stop when energy less than etsmin
      if(ec(ibkec+no,ipomp+1) .le. etsmin * 1.d+6 )then ! electron absorption
         dexc_ene     = dexc_ene + ec(ibkec+no,ipomp+1)

         if(wvlin.gt.0.d0) then
           if(ec(ibkec+no,ipomp+1).gt.wvlin) then
            eleemit = dble(int(ec(ibkec+no,ipomp+1)/wvlin))
            if(eleemit.lt.1.d0) eleemit = 1.d0
           endif
           ehnum = ehnum + eleemit -1.d0
           if(no.eq.1) ehnum = ehnum + 1.d0  ! source electron
         endif

         ec(ibkec+no,ipomp+1) = 0.0d0
      endif
!!!! secondary electron kill for primary electron range calculation
      dexc_ene = dexc_ene * 1.d-6 ! eV to MeV for tallying

c// atmrc increment for [t-interact] 20230719
      if( ityp .eq. 12 ) then
         atmrc(6,id_reac) = atmrc(6,id_reac) + 1.d0
      elseif( ityp .eq. 13 ) then
         atmrc(7,id_reac) = atmrc(7,id_reac) + 1.d0
      endif
! plasmon event
      if(id_reac.eq.2 ) then
        if(i.eq.1) then
          if( ityp .eq. 12 ) then
             atmrc(6,8) = atmrc(6,8) + 1.d0
          elseif( ityp .eq. 13 ) then
             atmrc(7,8) = atmrc(7,8) + 1.d0
          endif
        endif
      endif

*-----------------------------------------------------------------------
*        data up in bank
*-----------------------------------------------------------------------
      nclsts  =  0

      call etsdataup02 ! primary electron data up

c // Ionized electrons //
      numpal(12) = 0
      rumpal(12) = 0.d0
      numpal(14) = 0
      rumpal(14) = 0.d0

      wga        = wt(ibkwt+no,ipomp+1)
      iplasm = 0    ! for plasmon emit mulit electron 20220513
c // Ionized electrons and excited electrons//
      if(icoll.eq.1)then
! plasmon emit multi electrons
        if(id_reac.eq.2 .and. i.eq.1 .and. nplasm .gt. 1) then
               e2nd = e2nd/dble(nplasm)
        endif
 1001  continue
! 2nd electron stop when energy less than etsmin
        if(e2nd .le. etsmin * 1.d+6 )then ! electron absorption
           dexc_ene     = dexc_ene + e2nd * 1.d-6
! 20230619 epsilon value correction for calculating charged particles
         if(wvlin.gt.0.d0) then
           if(e2nd.gt.wvlin) then
            eleemit = dble(int(e2nd/wvlin))
            if(eleemit.lt.1.d0) eleemit = 1.d0
           endif
           ehnum = ehnum + eleemit -1.d0
         endif
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
! plasmon emit multi electrons
        if(id_reac.eq.2 .and. i.eq.1 .and. nplasm .gt. 1) then
           iplasm = iplasm + 1
           if(iplasm.lt.nplasm) goto 1001
        endif

      endif

! auger electorn dataup copied itsdataup
      if( lng_rel.ne.0.0) then
      do i = 1, lng_rel ! atomic relaxation
       
! secondary electron lower than etsmin
        if(ktp_rel(i) .ne. 22 .and. 
     &     eng_rel(i)* 1.d+6  .le. etsmin * 1.d+6 ) then
         if(wvlin.gt.0.d0) then
           if(eng_rel(i)* 1.d+6 .gt.wvlin) then
            eleemit = dble(int(eng_rel(i)* 1.d+6 /wvlin))
            if(eleemit.lt.1.d0) eleemit = 1.d0
           endif
           ehnum = ehnum + eleemit
         endif
         goto 1003
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
 1003    continue
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

c// atmrc increment for [t-interact] 20230719
      if(wvlin.gt.0.d0 .and. ehnum.gt.0.d0) then
        if( ityp .eq. 12 ) then
           atmrc(6,id_reac) = atmrc(6,id_reac) + ehnum
        elseif( ityp .eq. 13 ) then
           atmrc(7,id_reac) = atmrc(7,id_reac) + ehnum
        endif
      endif


      ec(ibkec+no,ipomp+1) = ec(ibkec+no,ipomp+1) * 1.d-6   ! from (eV) to (MeV)
      e(ibke+no,ipomp+1)   = e(ibke+ no,ipomp+1)  * 1.d-6   ! from (eV) to (MeV)

      end

c================================================
      subroutine ion_vec02(i,xip,bebw2nd)
c===============================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)
      dimension velp(3),vels(3),xip(6)
      common / ele1st / u1st,v1st,w1st
!$OMP THREADPRIVATE(/ele1st/)
      common / ele2nd / e2nd,u2nd,v2nd,w2nd
!$OMP THREADPRIVATE(/ele2nd/)
      common / etsexe / dexc_ene
!$OMP THREADPRIVATE(/etsexe/)
      common / etsminmax / etsmin, etsmax

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)

      common / cschng / change, fit
      common / Plasmonemit/ nplasm


      call ion_ang02(i,xip,bebw2nd,pthe,pphi,sthe,sphi)
      call vel_vec02(velp,pthe,pphi)
      call vel_vec02(vels,sthe,sphi)

        u1st =  velp(1)
        v1st =  velp(2)
        w1st =  velp(3)



! Si elf version  bebw2nd = w2nd = elf
! elf distribution is considered as 2nd electron energy distribution......
       if(i.eq.1) then
! plasmnon emit multiple electrons
        dexc_ene   =  1.1d0*dble(nplasm)
        e2nd       =  bebw2nd-dexc_ene
       else
        dexc_ene   =  xip(i)
        e2nd       =  bebw2nd-xip(i)
       endif

      ec(ibkec+no,ipomp+1) =  e(ibke+no,ipomp+1) - dexc_ene - e2nd
      if(ec(ibkec+no,ipomp+1) .lt. 0.d0)then
        dexc_ene   =  xip(i)
        e2nd       =  0.d0
      endif
      ec(ibkec+no,ipomp+1) =  e(ibke+no,ipomp+1) - dexc_ene - e2nd
        u2nd  =  vels(1)
        v2nd  =  vels(2)
        w2nd  =  vels(3)
      end

c================================================
      subroutine exc_vec02    ! currently not used
c===============================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)
      dimension velp(3),vels(3),xip(5)
      common / ele2nd / e2nd,u2nd,v2nd,w2nd
!$OMP THREADPRIVATE(/ele2nd/)
      common / etsexe / dexc_ene
!$OMP THREADPRIVATE(/etsexe/)
      common / etsminmax / etsmin, etsmax

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)

      common / cschng / change, fit

      call exc_ang02(pthe,pphi,sthe,sphi)
      call vel_vec02(velp,pthe,pphi)
      call vel_vec02(vels,sthe,sphi)

      u(ibku+no,ipomp+1)   =  velp(1)
      v(ibkv+no,ipomp+1)   =  velp(2)
      w(ibkw+no,ipomp+1)   =  velp(3)

      u2nd  =  vels(1)
      v2nd  =  vels(2)
      w2nd  =  vels(3)


      end


c====================================================
      subroutine ion_ang02(ip,xip,w2nd,pthe,pphi,sthe,sphi)
c====================================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)
      parameter (pi=3.1415926535897932385d0)
      dimension xip(6)
      common / idumc /idum
      real*8 ran2

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)
      common / angfact/ tangfe(10),tangfi(10),angfe(10),
     &                  angfi(10),iangfe,iangfi


      dmc2 = 511.003d+3
      wts  = w2nd

      if(ip.eq.1) then
        t2   = e(ibke+no,ipomp+1) - wts
! 2021/07/19 plasmon excitation will not change the electron angle
         pthe = 0.0
         pphi = 0.0
! secondary electrons
         if(wts.le.50.d0)then
          sthe = dble(ran2(idum)) * 1.d0 * pi
         else
          sthe = pi/4.d0 + dble(ran2(idum)) * pi/4.d0
         endif
          sphi = dble(ran2(idum)) * 2.d0 * pi
         return
      else
        t2   = e(ibke+no,ipomp+1) - xip(ip) - wts
      endif

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

      pthe = 0.0
      if(t2.ge.100.d0)then
        fang = 1.d0
        if(p3.gt.1.d0)then
         pthe = pi * 0.5d0
        else
!---  scattering angle changed for high energy electron 2021.04
          if(e(ibke+no,ipomp+1).le.tangfi(1)) then
                    fang =  angfi(1)
          else
             do i = 2,iangfi
                 if(e(ibke+no,ipomp+1).gt.tangfi(i-1)
     &        .and. e(ibke+no,ipomp+1).le.tangfi(i)) then
                    fang =  angfi(i)
                endif
             enddo
          endif
          pthe = dasin(p3) * fang
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

      end
c====================================================
      subroutine exc_ang02(pthe,pphi,sthe,sphi)  ! currently not used
c====================================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)
      parameter (pi=3.1415926535897932385d0)
      common / idumc /idum
      real*8 ran2

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)

      common / ele2nd / e2nd,u2nd,v2nd,w2nd
!$OMP THREADPRIVATE(/ele2nd/)
      common / etsexe / dexc_ene
!$OMP THREADPRIVATE(/etsexe/)


      dmc2 = 511.003d+3
      t2   = e(ibke+no,ipomp+1) - dexc_ene - e2nd
      wts  = e2nd

      tp   = t2/2.d0/dmc2
      p1   = wts/t2
      p2   = (1.d0 - p1) * tp + 1.d0
      p3   = dsqrt(p1/p2)
      if(t2.ge.100.d0)then
      pthe = dasin(p3)        * 1.5d0
      else
      pthe = dble(ran2(idum)) * pi / 4.d0
      endif

      wpts = wts/2.d0/dmc2
      s1   = 1.d0 - p1
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


      end

c================================================
      subroutine vel_vec02(vel,the,phi)
c===============================================
      implicit real*8(a-h,o-z)
      dimension vel(3),a(4,4)

      sx = 0.d0
      sy = 0.d0
      sz = 1.d0 ! CM->Lab angular conversion is done in reac.f

      the1 = acos(min(sz,1.d0))
      phi1 = atan2(sy,sx)

      xts = sin(the+the1) * cos(phi1)
      yts = sin(the+the1) * sin(phi1)
      zts = cos(the+the1)

      cs = (1.d0-cos(phi))
      a(1,1) = sx**2 * cs +      cos(phi)
      a(1,2) = sx*sy * cs - sz * sin(phi)
      a(1,3) = sz*sx * cs + sy * sin(phi)
      a(1,4) = 0.d0

      a(2,1) = sx*sy * cs + sz * sin(phi)
      a(2,2) = sy**2 * cs +      cos(phi)
      a(2,3) = sy*sz * cs - sx * sin(phi)
      a(2,4) = 0.d0

      a(3,1) = sz*sx * cs - sy * sin(phi)
      a(3,2) = sy*sz * cs + sx * sin(phi)
      a(3,3) = sz**2 * cs +      cos(phi)
      a(3,4) = 0.d0

      a(4,1) = 0.d0
      a(4,2) = 0.d0
      a(4,3) = 0.d0
      a(4,4) = 1.d0

      vel(1) = a(1,1)*xts + a(1,2)*yts + a(1,3)*zts
      vel(2) = a(2,1)*xts + a(2,2)*yts + a(2,3)*zts
      vel(3) = a(3,1)*xts + a(3,2)*yts + a(3,3)*zts

      sqvel  = sqrt( vel(1)**2 + vel(2)**2 + vel(3)**2 )
      vel(1) = vel(1)/sqvel
      vel(2) = vel(2)/sqvel
      vel(3) = vel(3)/sqvel


      end


c====================================================
      subroutine els_ics02(ene,xels)
c====================================================
      implicit real*8(a-h,o-z)
      parameter (pi=3.1415926535897932385d0)

      dmc2 = 511.003d+3
      z    = 14.0
      re   = 2.8179d-13
      tp   = ene/dmc2
      bt2  = 1.d0 - 1.d0/(1.d0 + tp)**2

      if(ene.le.5.d+4)then
      etc  = 1.198
      else
      etc  = 1.13 + 3.76 * z**2 / 137.d0**2 / bt2
      endif

      eta  = etc * 1.7d-5* z**(2.d0/3.d0) / tp / (tp + 2.d0)
      s1   = pi * re **2 * z * (z + 1.d0)
      s2   = (1.d0 - bt2) / bt2**2
      s3   = 1.d0 / eta /(eta + 1.d0)
      cs   = s1 * s2 * s3
      xels = cs


      end


c====================================================
      subroutine els_dcs02(elsthe,elsphi)
c====================================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)
      parameter (pi=3.1415926535897932385d0)
      dimension dels(0:180),sd(0:180)
      common / idumc /idum
      real*8 ran2

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)
      common / angfact/ tangfe(10),tangfi(10),angfe(10),
     &                  angfi(10),iangfe,iangfi

      dmc2 = 511.003d+3
      zp   = 14.0
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

       fang = 1.0d0
!---  scattering angle changed for high energy electron 2021.04
        if(e(ibke+no,ipomp+1).le.tangfe(1)) then
          fang = angfe(1)
        else
           do ia = 2,iangfe
               if(e(ibke+no,ipomp+1).gt.tangfe(ia-1)
     &      .and. e(ibke+no,ipomp+1).le.tangfe(ia)) then
         fang = angfe(ia)
              endif
           enddo
        endif
        elsthe = dble(i-ran2(idum)) * pi/180.d0*fang

        exit
       endif
      enddo
      elsphi = dble(ran2(idum)) * 2.d0 * pi

      end

c====================================================
      subroutine ion_ics02(ene,xion)
c====================================================
      implicit real*8(a-h,o-z)
      parameter (a0=5.29177249D-9)
      dimension xion(6)
      common / csion02 / dion(500),si(6,500),xip(6)
      common / nline02 / nion,nexc,nvib,nrot,ndea,nphn,nlow

      real(8),allocatable :: ep(:),s(:),dsdwi(:)

      common / cschng / change, fit
      common /DrudeParam/ Ej(6),Aj(6),wj(6)

      dimension y(7),z(7),xip1(7),xion1(7)

      pi    = 3.1415926535897932385d0
      ry    = 27.2116d0 * 0.5d0
      alpha = 1.d0/137.035999679d0
      dmc2   = 511.003d+3

      AUCM = 5.29177D-9 !atomic unit length
      DNCM = 4.9939D+22   !moleculer density /cm3
      em   = 9.1093829D-31   !electron mass kg
      h    = (6.626069D-34)/(2.0*PI)    !planck constant J*s/2pi -> kg*(m/s)^2/2pi
      evJ  = 1.602176634D-19                  ! J to eV  1 eV = 1.602176634x10-19 J
      AU   = ry
! sum rule correction
        corr = 1.0        ! without sum rule correction
!    Drude parameters
      Epl = 16.7
      mode = 6
      allocate(ep(mode),s(mode),dsdwi(mode))

      xip1(1)= 6.52d0  ; y(1)=30.33d0  ; z(1)=1.33d0;
      xip1(2)= 6.55d0  ; y(2)=30.59d0  ; z(2)=0.67d0;
      xip1(3)= 13.63d0 ; y(3)=48.21d0  ; z(3)=2.D0;
      xip1(4)= 107.98d0; y(4)=350.99d0 ; z(4)=4.D0;
      xip1(5)= 108.67d0; y(5)=354.24d0 ; z(5)=2.D0;
      xip1(6)= 151.55d0; y(6)=373.08d0 ; z(6)=2.D0;
      xip1(7)= 1828.5d0; y(7)=2551.0d0 ; z(7)=2.D0;


       S2  = 0.0
        do i=1,mode
          S(i) = 0.D0
          xion(i) = 0.D0
          xip(i) = Ej(i)
        enddo
! escape if the energy less than minimum value in table
       if(ene.lt.dion(1)) goto 2000

      if(ene.le.change)then  ! dion(nion) = change = 100 keV, change is used in water etsmode
      do j=1,6
        xion(j) = 0.d0
       if(j.ne.1 .and.ene.lt.xip(j))then
         xion(j)=0.d0
         goto 1000
       else
         do i=1,nion
         if(ene.eq.dion(i))then
         xion(j) = si(j,i)
         goto 1000
         else
         endif
         enddo
       endif

      do i=1,nion-1
      if(ene.ge.dion(i).and.ene.lt.dion(i+1)) then
         x1 = dion(i)
         x2 = dion(i+1)
         y1 = si(j,i)
         y2 = si(j,i+1)
      else; endif
      enddo

         a1 = (y1-y2)/(x1-x2)
         b1 = y1 - a1 * x1
         xion(j) = a1 * ene + b1
 1000 continue

       if(xion(j).lt.0.0) xion(j)=0.0
      enddo
      else
! KIM formula as water etsmode for high energy electron
        t = ene
        do j=1,7
          xion1(j) = 0.0
        enddo
        do j=1,7
          b = xip1(j); u = y(j); ds = z(j)
          st = t/b   ; su = u/b
          tp = t/dmc2; bp = b/dmc2; up = u/dmc2
          bt2 = 1.D0 - 1.D0/(1.D0 + tp)**2
          bb2 = 1.D0 - 1.D0/(1.D0 + bp)**2
          bu2 = 1.D0 - 1.D0/(1.D0 + up)**2

          s1 = 4.D0 * pi * a0**2 * alpha**4 * ds
          s2 = (bt2 + bb2 + bu2) * 2.D0 * bp
          ss  = s1 / s2
          f1 = 0.5d0 *
     .      (dlog( bt2/(1.d0-bt2) ) - bt2 -dlog(2.D0*bp)) *
     .      (1.d0 - 1.d0/st**2)
          f2 = 1.d0 - 1.d0/st
          f3 = - dlog(st)/(st+1.d0) *
     .      (1.d0 + 2.d0*tp)/(1.d0 + tp/2.d0)**2
          f4 = bp**2/(1.d0 + tp/2.d0)**2 * (st - 1.d0)/2.d0
           xion1(j) = ss * (f1 + f2 + f3 + f4)
        enddo

         xion(6) = fit * xion1(7)                 ! 1s
         xion(5) = fit * xion1(6)                 ! 2s
         xion(4) = fit * (xion1(5) + xion1(4))    ! 2p
         xion(3) = fit * xion1(3)                 ! 3s
         xion(2) = fit * (xion1(2) + xion1(1))    ! 3p

      endif
 2000 continue
      deallocate(ep,s,dsdwi)

      end
c====================================================
      subroutine ion_dcs0(itrn,bebw2nd)
c====================================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)
      parameter (pi=3.1415926535897932385d0)
      parameter (a0=5.29177249d-9)
      dimension sdion(0:10000),sd(0:10000)
      dimension wx(7),wy(7),wz(7)
      common /DrudeParam/ Ej(6),Aj(6),wj(6)
      common / idumc /idum
      real*8 ran2

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)

      alpha = 1.d0/137.035999679d0
      dmc2  = 511.003d+3

      wx(1)= 6.52d0  ; wy(1)=30.33d0  ; wz(1)=1.33d0;
      wx(2)= 6.55d0  ; wy(2)=30.59d0  ; wz(2)=0.67d0;
      wx(3)= 13.63d0 ; wy(3)=48.21d0  ; wz(3)=2.D0;
      wx(4)= 107.98d0; wy(4)=350.99d0 ; wz(4)=4.D0;
      wx(5)= 108.67d0; wy(5)=354.24d0 ; wz(5)=2.D0;
      wx(6)= 151.55d0; wy(6)=373.08d0 ; wz(6)=2.D0;
      wx(7)= 1828.5d0; wy(7)=2551.0d0 ; wz(7)=2.D0;

       ts = e(ibke+no,ipomp+1)
       qt  = 1.d0
       nq = 1000


      do 2000 ip=0,nq

       b  = wx(itrn)  ; uts = wy(itrn)  ; ds = wz(itrn)
      if(ts.lt.b) return
      wi  = 0.d0
      wf  = (ts-b)/2.d0
      hq  = (wf-wi)/nq
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

 2000 continue
 2001 continue
      ss = 0.d0
      do i=0,nq
      ss = ss + sdion(i)
      enddo

      sd(0) = sdion(0)/ss
      do i=1,nq
      sd(i) = sd(i-1) + sdion(i)/ss
      enddo

      pk = dble(ran2(idum))
      if(pk.gt.0.d0 .and.pk.le.sd(0))   bebw2nd = wx(itrn)
      do i=0,nq-1
      if(pk.gt.sd(i).and.pk.le.sd(i+1)) then
       bebw2nd = dble(i) * hq
      endif
      enddo

       if(itrn.eq.1.or.itrn.eq.2) echeck = Ej(2)
       if(itrn.ge.3) echeck = Ej(itrn-1)


      if(bebw2nd.lt.echeck) goto 2001

      end

c====================================================
      subroutine ion_dcs02(itrn,bebw2nd)
c====================================================
      use MMBANKMOD
      implicit real*8(a-h,o-z)
      parameter (a0=5.29177249D-9)
      real(8),allocatable :: ep(:)
      dimension elf(0:5000), sd(0:5000)

      common / cschng / change, fit
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp

      common / etsprob02 / eprob(12),xion(10),xexc(8),prb2,prb3
!$OMP THREADPRIVATE(/etsprob02/)
      common /DrudeParam/ Ej(6),Aj(6),wj(6)

      pi    = 3.1415926535897932385d0
      ry    = 27.2116d0 * 0.5d0
      alpha = 1.d0/137.035999679d0
      dmc2   = 511.003d+3

      AUCM = 5.29177D-9 !atomic unit length
      DNCM = 4.9939D+22   !moleculer density /cm3
      em   = 9.1093829D-31   !electron mass kg
      h    = (6.626069D-34)/(2.0*PI)    !planck constant J*s/2pi -> kg*(m/s)^2/2pi
      evJ  = 1.602176634D-19                  ! J to eV  1 eV = 1.602176634x10-19 J
      AU   = ry
      ene =  e(ibke+no,ipomp+1)
      EV  = ene
! sum rule correction
        corr = 1.0        ! without sum rule correction
!    Drude parameters
      Epl = 16.7d0
      mode = 6
      allocate(ep(mode))
cc //  parameters of dielectric constants //
      if( EV .gt. 5.d+3 ) then
        mesh    =  5000
        emesh   =  1.d0
      else
        mesh    =  5000
        emesh   =  EV / dble(mesh)
      endif

       do j = 0,mesh
        elf(j)  =  0.d0
       enddo

      i = itrn
      do 1000 j = 1, mesh
        WEV = dble(j) * emesh
        if(WEV.ge.EV)  then
            goto 1000            ! the EV is smaller than binding energy then loop( The maxmum energy of WEV is adjusted to maximam Ej)
        endif
        qp = (2.0d0*em*eVJ)**0.5d0/h*(EV**0.5d0+(EV-WEV)**0.5d0)  ! converted J to eV in first term, q is 1/m
        qm = (2.0d0*em*eVJ)**0.5d0/h*(EV**0.5d0-(EV-WEV)**0.5d0)  ! converted J to eV in first term, q is 1/m
        qb = (log10(qp)-log10(qm))/5.D+1
        elf(j) = 0.0
       do qi =log10(qm),log10(qp),qb             ! q integration loop
        q  = 10.0d0**(qi)
        ep(i) = 0.0
          Eq=Ej(i)+(h*q)**2.0d0/(2.0d0*em)/eVJ
          ep(i) = corr*Epl**2.0d0*Aj(i)*wj(i)*WEV/
     &         ((WEV**2.0d0-Eq**2.0d0)**2.0d0+(wj(i)*WEV)**2.0d0)
         if(i.ne.1) then                  ! heviside step function for absorbed edge
            if(WEV.le.Ej(i)) ep(i) = 0.0
         endif
          elf(j) = elf(j) + 1.0d0/(pi*DNCM*AUCM*EV)*ep(i)/q
     &               *(10.0d0**(qi+qb)-10.0d0**qi)  ! using eq.(6) in NIMB 288 (2012) 66-73
       enddo

 1000   CONTINUE
 1001   CONTINUE
        ss  =  0.d0
      do j = 0, mesh
        ss  =  ss + elf(j)
      enddo
        if(ss.eq.0.0) then
          if(EV-Ej(itrn) .gt. 0.0d0) then
             bebw2nd = Ej(itrn)
             deallocate(ep)
             return
          else
          endif
        endif
        sd(0)  =  elf(0) / ss
      do j = 1, mesh
        sd(j)  =  sd(j-1) + elf(j) / ss
      enddo

        pk  =  dble( ran2(idum) )
      do j = 0, mesh-1
        if( pk .gt. sd(j) .and. pk .le. sd(j+1)) then
          w2nd = dble(j) * emesh
          exit
        endif
      enddo
      if(w2nd.lt.Ej(i)) then   ! transfer energy determination error
           goto 1001
      endif
        if( ene - w2nd .lt. 0.d0 ) w2nd = 0.d0
      bebw2nd = w2nd
      deallocate(ep)
      end

c====================================================
      subroutine exc_ics02(ene,xexc)
c====================================================
      implicit real*8 (a-h,o-z)
      parameter (pi=3.1415926535897932385d0)
      parameter (a0=5.29177249d-9) !(cm2)
      dimension xexc(1)
      common / cslow / elow(100),se(1,100),sp(1,100)
      common / dcsexc / deexc(500),de(1,500),deprob(1,500),deall
      common / dcsphn / dephn(100),dp(1,100),dpprob(1,500),dpall

      common / nline02 / nion,nexc,nvib,nrot,ndea,nphn,nlow
      common / heexc / ea1(8)
      dimension dm(7),cs(7)

      common / cschng / change, fit

      ry    = 27.2116d0 * 0.5d0
      dmc2   = 511.003d+3

      xexc(1) = 0.d0
      if(ene.lt.elow(1).or.ene.gt.elow(nlow)) return

      do i = 1, nlow
       if(elow(i-1).lt.ene.and.elow(i).ge.ene) then
          a = (se(1,i)-se(1,i-1))/(elow(i)-elow(i-1))
          b = se(1,i)-a*elow(i)
          xexc(1) = a*ene+b
       endif
      enddo

       if(xexc(1).le.0.0) then
         xexc(1)=0.d0
       endif
       end
c====================================================
      subroutine exc_dcs02(ene,ip,exene)
c====================================================
      implicit real*8 (a-h,o-z)
      parameter (pi=3.1415926535897932385d0)
      parameter (a0=5.29177249d-9) !(cm2)
      common / dcsexc / deexc(500),de(1,500),deprob(1,500),deall
      common / dcsphn / dephn(100),dp(1,100),dpprob(1,500),dpall
      common / nline02 / nion,nexc,nvib,nrot,ndea,nphn,nlow

      real(8),allocatable :: dsdwe(:),dsprobe(:)
      allocate(dsdwe(nexc),dsprobe(nexc))
      do i = 1,nexc
         dsdwe(i) = 0.d0
         dsprobe(i) = 0.d0
      enddo

      EVHZ = 2.41797D+14
      AUHZ = 6.57968D+15
      AU   = 27.2116D0 / 2.D0  !atomic unit energy
      AUCM = 5.29177D-9 !atomic unit length
      DNCM = 4.9939D+22   !moleculer density /cm3
      em   = 9.1093829D-31   !electron mass kg
      h    = (6.626069D-34)/(2.0*PI)    !planck constant J*s/2pi -> kg*(m/s)^2/2pi
      evJ  = 1.602176634D-19                  ! J to eV  1 eV = 1.602176634x10-19 J


        EV   = ene
        dsum = 0.d0
        dmax = 0.d0
        do i = 1, nexc
         WEV = deexc(i)
         if(WEV.gt.EV) exit

         A = WEV/EV
         X1 = (1.D0-A)*DLOG(4.D0/A)
         X2 = -7.D0/4.D0*A
         X3 = A**(1.5D0)
         X4 = -33.D0/32.D0*A**2
         XX = X1 + X2 + X3 + X4
         dsdwe(i) = de(1,i)*XX*1.D0/2.D0/PI/EV/AUCM/DNCM
         if(dsdwe(i).lt.0.d0) dsdwe(i) = 0.d0

         if(dsdwe(i).ne.0.d0) dmax = deexc(i)
         dsum = dsum + dsdwe(i)
         if(i.eq.1) dsprobe(i) = dsdwe(i)
         if(i.gt.1) dsprobe(i) = dsprobe(i-1) + dsdwe(i)
        enddo

! calculate the ELF by Ashley formula 20220531
      pk  =  dble( ran2(idum) )*dsum
      ii = -1
      exene = 0.d0
      do i = 1, nexc
        if(pk .lt.dsprobe(1)) ii = 0
        if( pk .lt. dsprobe(i) .and. pk .ge. dsprobe(i-1)) then
          ii  = i
          exit
        endif
      enddo

      if(pk.ge.dsprobe(nexc)) then
         ii = -1
         exene = dmax
      endif
      if(ii.eq.0) then
        a = deexc(1)/dsprobe(1)
        exene = a*pk
      else if(ii.gt.0) then
        a = (deexc(ii)-deexc(ii-1))/(dsprobe(ii)-dsprobe(ii-1))
        b = deexc(ii)-a*dsprobe(ii)
        exene = a*pk+b
      endif

      deallocate( dsdwe,dsprobe )


      end

c====================================================
      subroutine phn_ics02(ene,xphn)
c====================================================
      implicit real*8 (a-h,o-z)
      common / cslow / elow(100),se(1,100),sp(1,100)
      common / dcsphn / dephn(100),dp(1,100),dpprob(1,500),dpall
      common / nline02 / nion,nexc,nvib,nrot,ndea,nphn,nlow


      ry    = 27.2116d0 * 0.5d0
      dmc2   = 511.003d+3

      xphn = 0.d0
      if(ene.lt.elow(1).or.ene.gt.elow(nlow)) return

      do i = 1, nlow
       if(elow(i-1).lt.ene.and.elow(i).ge.ene) then
          a = (sp(1,i)-sp(1,i-1))/(elow(i)-elow(i-1))
          b = sp(1,i)-a*elow(i)
          xphn = a*ene+b
       endif
      enddo
      if(xphn.lt.0.0) xphn = 0.d0

      end
c====================================================
      subroutine phn_dcs02(ene,ip,phnene)
c====================================================
      implicit real*8 (a-h,o-z)
      parameter (pi=3.1415926535897932385d0)
      parameter (a0=5.29177249d-9) !(cm2)
      common / dcsexc / deexc(500),de(1,500),deprob(1,500),deall
      common / dcsphn / dephn(100),dp(1,100),dpprob(1,500),dpall
      common / nline02 / nion,nexc,nvib,nrot,ndea,nphn,nlow

      real(8),allocatable :: dsdwp(:),dsprobp(:)
      allocate(dsdwp(nphn),dsprobp(nphn))
      do i = 1,nphn
         dsdwp(i) = 0.d0
         dsprobp(i) = 0.d0
      enddo

      EVHZ = 2.41797D+14
      AUHZ = 6.57968D+15
      AU   = 27.2116D0 / 2.D0  !atomic unit energy
      AUCM = 5.29177D-9 !atomic unit length
      DNCM = 4.9939D+22   !moleculer density /cm3
      em   = 9.1093829D-31   !electron mass kg
      h    = (6.626069D-34)/(2.0*PI)    !planck constant J*s/2pi -> kg*(m/s)^2/2pi
      evJ  = 1.602176634D-19                  ! J to eV  1 eV = 1.602176634x10-19 J

! if primary energy is less than the maximum value in table, adjust the probabirity table
        EV   = ene
        dsum = 0.d0
        dmax = 0.d0
        do i = 1, nphn
         WEV = dephn(i)
         if(WEV.gt.EV) exit

         A = WEV/EV
         X1 = (1.D0-A)*DLOG(4.D0/A)
         X2 = -7.D0/4.D0*A
         X3 = A**(1.5D0)
         X4 = -33.D0/32.D0*A**2
         XX = X1 + X2 + X3 + X4
         dsdwp(i) = dp(1,i)*XX*1.D0/2.D0/PI/EV/AUCM/DNCM
         if(dsdwp(i).lt.0.d0) dsdwp(i) = 0.d0

         if(dsdwp(i).ne.0.d0) dmax = dephn(i)
         dsum = dsum + dsdwp(i)
         if(i.eq.1) dsprobp(i) = dsdwp(i)
         if(i.gt.1) dsprobp(i) = dsprobp(i-1) + dsdwp(i)
        enddo

! calculate the ELF by Ashley formula 20220531
      pk  =  dble( ran2(idum) )*dsum
      ii = -1
      phnene = 0.d0
      do i = 1, nphn
        if(pk .lt.dsprobp(1)) ii = 0
        if( pk .lt. dsprobp(i) .and. pk .ge. dsprobp(i-1)) then
          ii  = i
          exit
        endif
      enddo

      if(pk.ge.dsprobp(nphn)) then
         ii = -1
         phnene = dmax
      endif
      if(ii.eq.0) then
        a = dephn(1)/dsprobp(1)
        phnene = a*pk
      else if(ii.gt.0) then
        a = (dephn(ii)-dephn(ii-1))/(dsprobp(ii)-dsprobp(ii-1))
        b = dephn(ii)-a*dsprobp(ii)
        phnene = a*pk+b
      endif

      deallocate( dsdwp,dsprobp )


      end



c================================================
      subroutine silicon_db(ierrdb)
c================================================
      implicit real*8(a-h,o-z)
      common / csion02 / dion(500),si(6,500),xip(6)
      common / cslow / elow(100),se(1,100),sp(1,100)
      common / dcsexc / deexc(500),de(1,500),deprob(1,500),deall
      common / dcsphn / dephn(100),dp(1,100),dpprob(1,500),dpall
      common / nline02 / nion,nexc,nvib,nrot,ndea,nphn,nlow
      common / angfact/ tangfe(10),tangfi(10),angfe(10),
     &                  angfi(10),iangfe,iangfi
      common /DrudeParam/ Ej(6),Aj(6),wj(6)

      common / secondon/ m_second
      common / Plasmonemit/ nplasm

! T.Sato 2016/05/28, for batch.out
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

      logical   exex
! Parameter preparation
        Ej(1) = 16.7d0
        Ej(2) = 13.63d0
        Ej(3) = 107.98d0
        Ej(4) = 108.67d0
        Ej(5) = 151.55d0
        Ej(6) = 1828.5d0

        Aj(1) = 0.827410033d0
        Aj(2) = 0.286524921d0
        Aj(3) = 0.806319615d0
        Aj(4) = 0.755175152d0
        Aj(5) = 1.08870878d0
        Aj(6) = 0.501504064d0

        wj(1) = 2.643317147d0
        wj(2) = 184.4032914d0
        wj(3) = 81.2242336d0
        wj(4) = 82.97694003d0
        wj(5) = 105.4397747d0
        wj(6) = 1149.780696d0

! Hirata angle check 20210830
      iamark = 0
! default values
 100    if(iamark .eq. 0) then
           iangfe = 3
           tangfe(1) = 3.0d+02
           tangfe(2) = 1.0d+05
           angfe(1)  = 1.00
           angfe(2)  = 1.00
           angfe(3)  = 1.00
           iangfi = 3
           tangfi(1) = 1.5d+04
           tangfi(2) = 1.0d+05
           angfi(1)  = 1.2
           angfi(2)  = 2.2
           angfi(3)  = 1.0
        endif

! Hirata angle check 20210830
      m_second = 0
      imark    = 0

! hirata prevent to read the water parameter for Si TS-mode

! Hirata check parameter
      nplasm = 0
      imark    = 0
      imark  = 1
      nplasm = 1

!! Cross section open

      inquire( file = chfn(25)(1:ilfn(25))//'/electron/Si.dat',
     & exist = exex )
      if( exex .eqv. .false. ) then
       ierrdb=1
       return
      endif

      open(25,file=chfn(25)(1:ilfn(25))//'/electron/Si.dat',
     & status='old')

! ionization cross section
      read(25,*)nion
      do i=1,nion
        read(25,*)dion(i),si(1,i),si(2,i),
     &            si(3,i),si(4,i),si(5,i),si(6,i)
      enddo

! excitation OELF
      deall = 0.0
      read(25,*)nexc
      do i=1,nexc
        read(25,'(2e13.4)') deexc(i),de(1,i)
        deall = deall + de(1,i)
      enddo
      deprob(1,1) = de(1,1)/deall
      do i = 2,nexc
       deprob(1,i) = de(1,i)/deall + deprob(1,i-1)
      enddo


! phonon OELF
      dpall =0.0
      read(25,*)nphn
      do i=1,nphn
        read(25,'(2e13.4)') dephn(i),dp(1,i)
        dpall = dpall + dp(1,i)
      enddo
      dpprob(1,1) = dp(1,1)/dpall
      do i = 2,nphn
       dpprob(1,i) = dp(1,i)/dpall + dpprob(1,i-1)
      enddo

! low energy interaction cross section
      read(25,*)nlow
      do i=1,nlow
        read(25,'(3e15.7)') elow(i),se(1,i),sp(1,i)
      enddo

      close(25)


      end

c================================================
      subroutine etsdataup02
c===============================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)

      include 'param00.inc'

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp
!$OMP THREADPRIVATE(/icomon/)
      common / clustt / nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common / clustw / jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common / etsexe / dexc_ene
!$OMP THREADPRIVATE(/etsexe/)
      common / etsminmax / etsmin, etsmax
      common / ele1st / u1st,v1st,w1st
!$OMP THREADPRIVATE(/ele1st/)

         xgv     =   1.d-9                             ! from (eV) to (GeV)
         xmv     =   1.d-6                             ! from (eV) to (MeV)
         erm     =   rmtyp(ityp,ktyp) * 1.d-3          ! electron rest mass from (MeV) to (GeV)
         xme     =   9.10953d-31                       ! electron mass (kg)
         xjl     =   6.24146d+18                       ! from (eV) to (J)
         gvc     =   1.87115653d+18                    ! from (kg m/s) to (GeV/c)


         ej      =   e(ibke+no,ipomp+1) / xjl                  ! energy from (eV) to (J)
         vr      =   dsqrt( 2.D0 * ej / xme )          ! velocity (m/s)
         vx      =   vr  * u1st                        ! x-comp. of velocity (m/s)
         vy      =   vr  * v1st                        ! y-comp. of velocity (m/s)
         vz      =   vr  * w1st                        ! z-comp. of velocity (m/s)
         px      =   xme * vx * gvc                    ! x-comp. of momentum from (kg m/s) to (GeV/c)
         py      =   xme * vy * gvc                    ! y-comp. of momentum from (kg m/s) to (GeV/c)
         pz      =   xme * vz * gvc                    ! z-comp. of momentum from (kg m/s) to (GeV/c)
         etotal  =   dsqrt(px**2 + py**2 + pz**3 + erm**2)

         nclsts = nclsts + 1

         iclusts( nclsts ) =  7
         jclusts(  0, nclsts ) = 0
         jclusts(  1, nclsts ) = 0
         jclusts(  2, nclsts ) = 0
         jclusts(  3, nclsts ) = ityp
         jclusts(  4, nclsts ) = 0
         jclusts(  5, nclsts ) = ichgf(ityp,ktyp)
         jclusts(  6, nclsts ) = 0
         jclusts(  7, nclsts ) = ktyp
         jclusts(  8, nclsts ) = 0


         qclusts(  0, nclsts ) = 0.d0                      ! impact parameter
         qclusts(  1, nclsts ) = px                        ! (GeV/c)
         qclusts(  2, nclsts ) = py                        ! (GeV/c)
         qclusts(  3, nclsts ) = pz                        ! (GeV/c)
         qclusts(  4, nclsts ) = etotal                    ! (GeV)
         qclusts(  5, nclsts ) = rmtyp(ityp,ktyp) * 1.d-3  ! (GeV)
         qclusts(  6, nclsts ) = 0.d0                      ! from (eV) to (MeV)
         qclusts(  7, nclsts ) = ec(ibkec+no,ipomp+1) * xmv! from (eV) to (MeV)
         qclusts(  8, nclsts ) = wt(ibkwt+no,ipomp+1)      ! weight
         qclusts(  9, nclsts ) = 0.d0                      ! (ns)
         qclusts( 10, nclsts ) = 0.d0                      ! (cm)
         qclusts( 11, nclsts ) = 0.d0                      ! (cm)
         qclusts( 12, nclsts ) = 0.d0                      ! (cm)

      return
      end
