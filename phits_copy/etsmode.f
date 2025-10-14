c=================================================================
c
c  // Monte Carlo Track Stracture Code for low energy electron //
c
c                                          by Takeshi Kai in JAEA
c
c=================================================================
      subroutine etsflt(sig_macro)
c=================================================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)

      include 'param00.inc'
      include 'param.inc'

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common / etsprob / eprob(12),xion(5),xexc(6),prb2,prb3
!$OMP THREADPRIVATE(/etsprob/)
      common /celepcc/ enumpcc(kvlmax)

      dimension vel(3),ss(0:10)

      common / cschng / change, fit

      common /regdc/  idrg(kvlmax), idgr(kvmmax)

      change = 1.d+5        ! (100 keV)
      fit = 2.4942487E-18 / 1.9528578E-018


      ec(ibkec+no,ipomp+1) = ec(ibkec+no,ipomp+1) * 1.d+6   ! from (MeV) to (eV)
      e(ibke+no,ipomp+1)   = e(ibke+no,ipomp+1)   * 1.d+6   ! from (MeV) to (eV)
      ene = e(ibke+no,ipomp+1)


c // cross section (cm2) //

      call  els_ics(ene,xels)

      if( ene .gt. 1.d+5 ) then
        call  high_energy( ene, xion, xexc )
        goto 2000
      endif

      call  ion_ics(ene,xion)

      call  exc_ics(ene,xexc)

 2000 continue


      if(e(ibke+no,ipomp+1).gt.4.d0.and.e(ibke+no,ipomp+1).le.13.d0)then
            iprss = 1
            call  dea_ics(iprss,ene,xda1)
            iprss = 2
            call  dea_ics(iprss,ene,xda2)
            iprss = 3
            call  dea_ics(iprss,ene,xda3)

            xda1 = xda1 / 6.d0     ! Y.Matsuya 230818
            xda2 = xda2 / 6.d0     ! Y.Matsuya 230818
            xda3 = xda3 / 6.d0     ! Y.Matsuya 230818

      else
            xda1=0.d0  ;  xda2=0.d0  ;  xda3=0.d0
      endif

      if(e(ibke+no,ipomp+1).le.100.d0)then
            iprss = 1
            call  vib_ics(iprss,ene,xvb1)
            iprss = 2
            call  vib_ics(iprss,ene,xvb2)

            iprss = 1
            call  phn_ics(iprss,ene,xph1)
            iprss = 2
            call  phn_ics(iprss,ene,xph2)

            iprss = 1
            call  rot_ics(iprss,ene,xrt1)
            iprss = 2
            call  rot_ics(iprss,ene,xrt2)
      else
            xvb1=0.d0  ;  xvb2=0.d0
            xph1=0.d0  ;  xph2=0.d0
            xrt1=0.d0  ;  xrt2=0.d0
      endif

        prb1  = xels
        prb2  = 0.d0
      do ip=1,5
        prb2  = prb2 + xion(ip)
      enddo
        prb3  = 0.d0
      do ip=1,6
        prb3  = prb3 + xexc(ip)
      enddo
        prb4  = xda1   ;   prb5  = xda2   ;   prb6  = xda3
        prb7  = xvb1   ;   prb8  = xvb2
        prb9  = xph1   ;   prb10 = xph2
        prb11 = xrt1   ;   prb12 = xrt2


c // total cross section (sgm) //
      sgm      = prb1  +  prb2  +  prb3  + prb4   +  prb5  +  prb6
     .         + prb7  +  prb8  +  prb9  + prb10  +  prb11 +  prb12
        eprob(1)   =   prb1    /   sgm
        eprob(2)   =   prb2    /   sgm   +   eprob(1)
        eprob(3)   =   prb3    /   sgm   +   eprob(2)
        eprob(4)   =   prb4    /   sgm   +   eprob(3)
        eprob(5)   =   prb5    /   sgm   +   eprob(4)
        eprob(6)   =   prb6    /   sgm   +   eprob(5)
        eprob(7)   =   prb7    /   sgm   +   eprob(6)
        eprob(8)   =   prb8    /   sgm   +   eprob(7)
        eprob(9)   =   prb9    /   sgm   +   eprob(8)
        eprob(10)  =   prb10   /   sgm   +   eprob(9)
        eprob(11)  =   prb11   /   sgm   +   eprob(10)
        eprob(12)  =   prb12   /   sgm   +   eprob(11)


c // determination of mean free path (xlm) and flight length (fpl) //
        wanumpcc  = 3.318565377871046d+022
        wenumpcc  = 3.323785188937734d+023
        enumscale = enumpcc(mat)/wenumpcc


        sig_macro   =  sgm * (wanumpcc*enumscale) ! delete /fac (Takeshi 2022/0829)



      ec(ibkec+no,ipomp+1) = ec(ibkec+no,ipomp+1) * 1.d-6    ! from (eV) to (MeV)
      e(ibke+no,ipomp+1)   = e(ibke+ no,ipomp+1)  * 1.d-6    ! from (eV) to (MeV)

      end


************************************************************************
*                                                                      *
      subroutine etstrn(mark,markp,fpl,nbeta,itmak)
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
      include 'param.inc'

*-----------------------------------------------------------------------

      common / icomon / no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common / tlgeom / iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common / etsminmax / etsmin, etsmax
      common /regdc/  idrg(kvlmax), idgr(kvmmax)

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

*-----------------------------------------------------------------------
               coincd = coincd * 1.d+4

      return
      end


c=================================================================
      subroutine etsreac
c=================================================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)

      parameter(nAuger=6) ! T.Sato 2022/09/03 to consider actual energy distribution of Auger electron

      include 'param00.inc'
      include 'param.inc'

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
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

      common / etsprob / eprob(12),xion(5),xexc(6),prb2,prb3
!$OMP THREADPRIVATE(/etsprob/)
      common / ele2nd  / e2nd,u2nd,v2nd,w2nd
!$OMP THREADPRIVATE(/ele2nd/)
      common / etsexe  / dexc_ene
!$OMP THREADPRIVATE(/etsexe/)
      common / ele1st / u1st,v1st,w1st
!$OMP THREADPRIVATE(/ele1st/)

      common / cschng / change, fit

      common / etsminmax / etsmin, etsmax
      common / csion / dion(500),si(5,500),xip(5)
      common / csexc / eexc(500),se(6,500),eip(6)
      common / csdea / edea( 53),sd(3, 53),dip(3)
      common / csvib / evb1( 44),sv(2, 44),vip(2),evb2( 44)
      common / csphn / ephn(308),sp(2,308),pip(2)
      common / csrot / erot(451),sr(2,451),rip(2)
      common / idumc /idum
      real*8 ran2

      dimension vel(3),ss(0:10)

      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)
      common /etsion/ icoll
!$OMP THREADPRIVATE(/etsion/)

      common /etsart/ bgets(kvlmax), ebgets(kvlmax),
     &                wvets(kvlmax), ewvets(kvlmax)
      common /ets_wvalue/ wvlin
!$OMP THREADPRIVATE(/ets_wvalue/)


! T.Sato 2022/09/03 Auger database from EADL
      dimension AugerProb(nAuger),AugerEne(nAuger)
      data AugerProb/1.79E-01,2.95E-01,5.25E-01,5.36E-01,8.27E-01,
     & 9.94E-01/
      data AugerEne/4.79E+02,4.94E+02,4.94E+02,5.09E+02,5.09E+02,
     & 5.09E+02/

         u1st = 0.d0
         v1st = 0.d0
         w1st = 1.d0

      ec(ibkec+no,ipomp+1) = ec(ibkec+no,ipomp+1) * 1.d+6   ! from (MeV) to (eV)
      e(ibke+no,ipomp+1)   = e(ibke+no,ipomp+1)   * 1.d+6   ! from (MeV) to (eV)

         icoll = 0
         iagr  = 0
         e2nd  = 0.d0
         ehnum = 0.d0

         wvlin = ewvets( idgr(iblz(ibkblz+no,ipomp+1)) )

c // determiation of e, u, v, w //  
         prb   = dble(ran2(idum))
c // elastic scattering (id=1) //
      if(prb.le.eprob(1))then
         id_reac = 1      ! S.Abe 2018/02/07
   11 continue
         call els_dcs(xthe,xphi)
         call vel_vec(vel,xthe,xphi)
         u1st = vel(1)
         v1st = vel(2)
         w1st = vel(3)
         nclsts     = 1
         dexc_ene   = 0.d0

c // ionization (id=2) //
      else if(prb.gt.eprob(1).and.prb.le.eprob(2))then
         id_reac = 2      ! S.Abe 2018/02/07
         ss(0)  = 0.d0
         do ip=1,5


           if(e(ibke+no,ipomp+1).ge.xip(ip))then
           ss(ip) = ss(ip-1) + xion(ip)/prb2
           else
           ss(ip) = ss(ip-1) + 0.d0
           endif


         enddo

         pk = dble(ran2(idum))
         if(pk.gt.0.d0 .and.pk.le.ss(1))   i = 1
         do j=1,4
         if(pk.gt.ss(j).and.pk.le.ss(j+1)) i = j + 1
         enddo
         call ion_dcs(i,bebw2nd)
         call ion_vec(i,xip,bebw2nd)

         nclsts = 1
         icoll = 1

c Takeshi (2018/10/18) Auger electron, T.Sato 2022/09/03 to consider energy distribution taken from EADL
          if(i.eq.5)then
           pk = dble(ran2(idum))
           eagr=0.0
           do j=1,nAuger
            if(pk.lt.AugerProb(j)) then
             eagr = AugerEne(j)
             exit
            endif
           enddo
           if(eagr.ne.0.0) then
            dexc_ene = xip(i) - eagr
            ranx = dble(ran2(idum))
            rany = dble(ran2(idum))
            ranz = dble(ran2(idum))
            vx   = -0.5d0 + ranx
            vy   = -0.5d0 + rany
            vz   = -0.5d0 + ranz
            vv   = dsqrt(vx**2 + vy**2 + vz**2)
            uagr = vx / vv
            vagr = vy / vv
            wagr = vz / vv
            iagr = 1
           endif
          endif


c // electronic excitation (id=3) //
      else if(prb.gt.eprob(2).and.prb.le.eprob(3))then
         id_reac = 3      ! S.Abe 2018/02/07
         ss(0)  = 0.d0
         do ip=1,6

           if(e(ibke+no,ipomp+1).ge.eip(ip))then
           ss(ip) = ss(ip-1) + xexc(ip)/prb3
           else
           ss(ip) = ss(ip-1) + 0.d0
           endif

         enddo

         pk = dble(ran2(idum))
         if(pk.gt.0.d0 .and.pk.le.ss(1))   i = 1
         do j=1,5
         if(pk.gt.ss(j).and.pk.le.ss(j+1)) i = j + 1
         enddo

         nclsts       = 1

         ec(ibkec+no,ipomp+1) = e(ibke+no,ipomp+1) - eip(i)
         dexc_ene = e(ibke+no,ipomp+1) - ec(ibkec+no,ipomp+1)


c Takeshi (2018/10/18) collective excidated electron
          if(i.eq.6)then
           dexc_ene = 10.9d0
           id_reac = 2
           e2nd = 21.4d0 - 10.9d0
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
           icoll = 1
          endif


c Takeshi (2022/08/30) Diffuse band excitation
          if(i.eq.5 .and. dble(ran2(idum)).le.0.9d0)then
           dexc_ene = 10.9d0
           id_reac = 2
           e2nd = 14.1d0 - 10.9d0
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
           icoll = 1
          endif


c // dissociative electron attachment; OH(-) (id=4) //
      else if(prb.gt.eprob(3).and.prb.le.eprob(4))then
         id_reac = 4      ! S.Abe 2018/02/07
         nclsts       = 1
         ec(ibkec+no,ipomp+1) = etsmin * 1.d+6 ! (eV)
         dexc_ene     = e(ibke+no,ipomp+1) - ec(ibkec+no,ipomp+1)

c // dissociative electron attachment; O(-) (id=4) //
      else if(prb.gt.eprob(4).and.prb.le.eprob(5))then
         id_reac = 4      ! S.Abe 2018/02/07
         nclsts       = 1
         ec(ibkec+no,ipomp+1) = etsmin * 1.d+6 ! (eV)
         dexc_ene     = e(ibke+no,ipomp+1) - ec(ibkec+no,ipomp+1)

c // dissociative electron attachment; H(-) (id=4) //
      else if(prb.gt.eprob(5).and.prb.le.eprob(6))then
         id_reac = 4      ! S.Abe 2018/02/07
         nclsts       = 1
         ec(ibkec+no,ipomp+1) = etsmin * 1.d+6 ! (eV)
         dexc_ene     = e(ibke+no,ipomp+1) - ec(ibkec+no,ipomp+1)

c // vibration exciation (bending) (id=5) //
      else if(prb.gt.eprob(6).and.prb.le.eprob(7))then
         id_reac = 5      ! S.Abe 2018/02/07
         nclsts       = 1
         ec(ibkec+no,ipomp+1) = e(ibke+no,ipomp+1) - vip(1)
         dexc_ene     = e(ibke+no,ipomp+1) - ec(ibkec+no,ipomp+1)

c // vibration exciation (streching) (id=5) //
      else if(prb.gt.eprob(7).and.prb.le.eprob(8))then
         id_reac = 5      ! S.Abe 2018/02/07
         nclsts       = 1
         ec(ibkec+no,ipomp+1) = e(ibke+no,ipomp+1) - vip(2)
         dexc_ene     = e(ibke+no,ipomp+1) - ec(ibkec+no,ipomp+1)

c // phonon excitation (id=6) //
      else if(prb.gt.eprob(8).and.prb.le.eprob(9))then
         id_reac = 6      ! S.Abe 2018/02/07
         nclsts       = 1
         ec(ibkec+no,ipomp+1) = e(ibke+no,ipomp+1) - pip(1)
         dexc_ene     = e(ibke+no,ipomp+1) - ec(ibkec+no,ipomp+1)

c // phonon excitation (id=6) //
      else if(prb.gt.eprob(9).and.prb.le.eprob(10))then
         id_reac = 6      ! S.Abe 2018/02/07
         nclsts       = 1
         ec(ibkec+no,ipomp+1) = e(ibke+no,ipomp+1) - pip(2)
         dexc_ene     = e(ibke+no,ipomp+1) - ec(ibkec+no,ipomp+1)

c // rotation excitation (id=7) //
      else if(prb.gt.eprob(10).and.prb.le.eprob(11))then
         id_reac = 7      ! S.Abe 2018/02/07
         nclsts       = 1
         ec(ibkec+no,ipomp+1) = e(ibke+no,ipomp+1) - rip(1)
         dexc_ene     = e(ibke+no,ipomp+1) - ec(ibkec+no,ipomp+1)

c // rotation excitation (id=7) //
      else ! if(prb.gt.eprob(11).and.prb.le.eprob(12))then
         id_reac = 7      ! S.Abe 2018/02/07
         nclsts       = 1
         ec(ibkec+no,ipomp+1) = e(ibke+no,ipomp+1) - rip(2)
         dexc_ene     = e(ibke+no,ipomp+1) - ec(ibkec+no,ipomp+1)
      endif


      if(ec(ibkec+no,ipomp+1) .le. etsmin * 1.d+6)then

        if(id_reac .eq. 2) then
          e2nd = 0.d0
        endif

        ec(ibkec+no,ipomp+1) = etsmin * 1.d+6
        dexc_ene     = e(ibke+no,ipomp+1) - ec(ibkec+no,ipomp+1)
! 20230619 epsilon value correction for calculating charged particles
         eleemit = 1.d0
         if(wvlin.gt.0.d0) then
           if(e(ibke+no,ipomp+1).gt.wvlin) then
            eleemit = dble(int(e(ibke+no,ipomp+1)/wvlin))
            if(eleemit.lt.1.d0) eleemit = 1.d0
           endif
           ehnum = ehnum + eleemit -1.d0
           if(no.eq.1) ehnum = ehnum + 1.d0  ! source electron
         endif
      endif

! 20230619 epsilon value correction for calculating charged particles
         eleemit = 1.d0
         if(wvlin.gt.0.d0) then
           if(e2nd.gt.wvlin) then
            eleemit = dble(int(e2nd/wvlin))
            if(eleemit.lt.1.d0) eleemit = 1.d0
           endif
           ehnum = ehnum + eleemit -1.d0
         endif

!--- end change NS 2020.04 del THREADPRIVATE

      dexc_ene = dexc_ene * 1.d-6 ! eV to MeV for tallying

      if( ityp .eq. 12 ) then
         atmrc(6,id_reac) = atmrc(6,id_reac) + 1.d0
      elseif( ityp .eq. 13 ) then
         atmrc(7,id_reac) = atmrc(7,id_reac) + 1.d0
      endif

! plasmon event 20230719 hirata
      if(id_reac.eq.2 .and. i.eq.6) then
        if( ityp .eq. 12 ) then
           atmrc(6,8) = atmrc(6,8) + 1.d0
        elseif( ityp .eq. 13 ) then
           atmrc(7,8) = atmrc(7,8) + 1.d0
        endif
      endif
*-----------------------------------------------------------------------
*        data up in bank
*-----------------------------------------------------------------------

      nclsts  =  0

      call etsdataup

c // Ionized electrons //
      numpal(12) = 0
      rumpal(12) = 0.d0
      numpal(14) = 0
      rumpal(14) = 0.d0

      wga        = wt(ibkwt+no,ipomp+1)

c // Ionized electrons //
      if(icoll.eq.1)then
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

         pel = sqrt( (e2nd*1.d-9)**2 + 2.d0 * e2nd * 1.d-9 *
     &                rmtyp(ityp,ktyp) * 1.d-3 ) ! momentum in GeV

         qclusts(  0, nclsts ) = 0.d0                      ! impact parameter
         qclusts(  1, nclsts ) = pel * u2nd                ! (GeV/c)
         qclusts(  2, nclsts ) = pel * v2nd                ! (GeV/c)
         qclusts(  3, nclsts ) = pel * w2nd                ! (GeV/c)
         qclusts(  4, nclsts ) = e2nd * 1.d-9 + rmtyp(ityp,ktyp) * 1.d-3! (GeV)
         qclusts(  5, nclsts ) = rmtyp(ityp,ktyp) * 1.d-3  ! (GeV)
         qclusts(  6, nclsts ) = 0.d0                      ! from (eV) to (MeV)
         qclusts(  7, nclsts ) = e2nd * 1.d-6              ! from (eV) to (MeV)
         qclusts(  8, nclsts ) = 1.d0                      ! weight change
         qclusts(  9, nclsts ) = 0.d0                      ! (ns)
         qclusts( 10, nclsts ) = 0.d0                      ! (cm)
         qclusts( 11, nclsts ) = 0.d0                      ! (cm)
         qclusts( 12, nclsts ) = 0.d0                      ! (cm)

         numpal(12) = numpal(12) + 1
         rumpal(12) = rumpal(12) + wga
      endif

c Takeshi (2018/10/18) Auger electron
      if(iagr.eq.1)then
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

         pel = sqrt( (eagr*1.d-9)**2 + 2.d0 * eagr * 1.d-9 *
     &                rmtyp(ityp,ktyp) * 1.d-3 ) ! momentum in GeV

         qclusts(  0, nclsts ) = 0.d0                      ! impact parameter
         qclusts(  1, nclsts ) = pel * uagr                ! (GeV/c)
         qclusts(  2, nclsts ) = pel * vagr                ! (GeV/c)
         qclusts(  3, nclsts ) = pel * wagr                ! (GeV/c)
         qclusts(  4, nclsts ) = eagr * 1.d-9 + rmtyp(ityp,ktyp) * 1.d-3! (GeV)
         qclusts(  5, nclsts ) = rmtyp(ityp,ktyp) * 1.d-3  ! (GeV)
         qclusts(  6, nclsts ) = 0.d0                      ! from (eV) to (MeV)
         qclusts(  7, nclsts ) = eagr * 1.d-6              ! from (eV) to (MeV)
         qclusts(  8, nclsts ) = 1.d0                      ! weight change
         qclusts(  9, nclsts ) = 0.d0                      ! (ns)
         qclusts( 10, nclsts ) = 0.d0                      ! (cm)
         qclusts( 11, nclsts ) = 0.d0                      ! (cm)
         qclusts( 12, nclsts ) = 0.d0                      ! (cm)

         numpal(12) = numpal(12) + 1
         rumpal(12) = rumpal(12) + wga

! 20230619 epsilon value correction for calculating charged particles
         eleemit = 1.d0
         if(wvlin.gt.0.d0) then
           if(eagr.gt.wvlin) then
            eleemit = dble(int(eagr/wvlin))
            if(eleemit.lt.1.d0) eleemit = 1.d0
           endif
           ehnum = ehnum + eleemit -1.d0
         endif

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
         qclusts(  8, nclsts ) = 1.d0                      ! weight change
         qclusts(  9, nclsts ) = 0.d0                      ! (ns)
         qclusts( 10, nclsts ) = 0.d0                      ! (cm)
         qclusts( 11, nclsts ) = 0.d0                      ! (cm)
         qclusts( 12, nclsts ) = 0.d0                      ! (cm)

         numpal(14) = numpal(14) + 1
         rumpal(14) = rumpal(14) + wga
        enddo


      endif

      if(wvlin.gt.0.d0 .and. ehnum.gt.0.d0) then
        if( ityp .eq. 12 ) then
           atmrc(6,id_reac) = atmrc(6,id_reac) + ehnum
        elseif( ityp .eq. 13 ) then
           atmrc(7,id_reac) = atmrc(7,id_reac) + ehnum
        endif
      endif

!--- change NS 2020.04 del THREADPRIVATE
!      ec(ibkec+no) = ec(ibkec+no) * 1.d-6   ! from (eV) to (MeV)
!      e(ibke+no)   = e(ibke+ no)  * 1.d-6   ! from (eV) to (MeV)
      ec(ibkec+no,ipomp+1) = ec(ibkec+no,ipomp+1) * 1.d-6   ! from (eV) to (MeV)
      e(ibke+no,ipomp+1)   = e(ibke+ no,ipomp+1)  * 1.d-6   ! from (eV) to (MeV)


      end


c================================================
      subroutine ion_vec(i,xip,bebw2nd)
c===============================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)
      dimension velp(3),vels(3),xip(5)
      common / ele1st / u1st,v1st,w1st
!$OMP THREADPRIVATE(/ele1st/)
      common / ele2nd / e2nd,u2nd,v2nd,w2nd
!$OMP THREADPRIVATE(/ele2nd/)
      common / etsexe / dexc_ene
!$OMP THREADPRIVATE(/etsexe/)
      common / etsminmax / etsmin, etsmax

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)

      common / cschng / change, fit

      call ion_ang(i,xip,bebw2nd,pthe,pphi,sthe,sphi)
      call vel_vec(velp,pthe,pphi)
      call vel_vec(vels,sthe,sphi)


        u1st =  velp(1)
        v1st =  velp(2)
        w1st =  velp(3)



        dexc_ene   =  xip(i)
        e2nd       =  bebw2nd


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
      subroutine exc_vec
c===============================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)
      dimension velp(3),vels(3),xip(5)
      common / ele2nd / e2nd,u2nd,v2nd,w2nd
!$OMP THREADPRIVATE(/ele2nd/)
      common / etsexe / dexc_ene
!$OMP THREADPRIVATE(/etsexe/)
      common / etsminmax / etsmin, etsmax

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)

      common / cschng / change, fit

      call exc_ang(pthe,pphi,sthe,sphi)
      call vel_vec(velp,pthe,pphi)
      call vel_vec(vels,sthe,sphi)

      u(ibku+no,ipomp+1)   =  velp(1)
      v(ibkv+no,ipomp+1)   =  velp(2)
      w(ibkw+no,ipomp+1)   =  velp(3)

      u2nd  =  vels(1)
      v2nd  =  vels(2)
      w2nd  =  vels(3)


      end


c====================================================
      subroutine ion_ang(ip,xip,w2nd,pthe,pphi,sthe,sphi)
c====================================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)
      parameter (pi=3.1415926535897932385d0)
      dimension xip(5)
      common / idumc /idum
      real*8 ran2

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)

      dmc2 = 511.003d+3
      t2   = e(ibke+no,ipomp+1) - xip(ip) - w2nd
      wts  = w2nd

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


      end


c====================================================
      subroutine exc_ang(pthe,pphi,sthe,sphi)
c====================================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)
      parameter (pi=3.1415926535897932385d0)
      common / idumc /idum
      real*8 ran2

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
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
      subroutine vel_vec(vel,the,phi)
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
      subroutine els_ics(ene,xels)
c====================================================
      implicit real*8(a-h,o-z)
      parameter (pi=3.1415926535897932385d0)

      dmc2 = 511.003d+3
      z    = 7.42d0
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
      subroutine els_dcs(elsthe,elsphi)
c====================================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)
      parameter (pi=3.1415926535897932385d0)
      dimension dels(0:180),sd(0:180)
      common / idumc /idum
      real*8 ran2

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)

      dmc2 = 511.003d+3
      zp   = 7.42d0
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

      end


c====================================================
      subroutine high_energy( ene, xion, xexc )
c====================================================
      implicit real*8 (a-h,o-z)


      dimension xion(5), x(5), y(5), z(5)
      dimension xexc(6), xexcg(8), ea(8), dm(7), cs(7)

      dimension wlion(5), wlexc(6)
      dimension scalion(5), scalexc(6)

      pi    = 3.1415926535897932385d0
      ry    = 27.2116d0 * 0.5d0
      alpha = 1.d0/137.035999679d0
      dmc2  = 511.003d+3
      a0    =  5.29177249d-9


      wlion(1) = 1.2980643E-20
      wlion(2) = 1.9812692E-19
      wlion(3) = 3.3105660E-19
      wlion(4) = 1.2827365E-19
      wlion(5) = 2.3101459E-20

      wlexc(1) = 1.4374928E-20
      wlexc(2) = 9.9007440E-21
      wlexc(3) = 1.8622625E-20
      wlexc(4) = 9.0543240E-21
      wlexc(5) = 1.3719891E-19
      wlexc(6) = 1.0701670E-18


        sss = 0.d0
      do i = 1, 5
        sss = sss + wlion(i)
      enddo
      do i = 1, 6
        sss = sss + wlexc(i)
      enddo

      do i = 1, 5
        scalion(i) = wlion(i) / sss
      enddo
      do i = 1, 6
        scalexc(i) = wlexc(i) / sss
      enddo

      x(1) =  12.61D0 ; y(1) =  122.90D0 ; z(1) = 2.D0 ;
      x(2) =  14.73D0 ; y(2) =  118.40D0 ; z(2) = 2.D0 ;
      x(3) =  18.55D0 ; y(3) =   97.39D0 ; z(3) = 2.D0 ;
      x(4) =  32.20D0 ; y(4) =  142.00D0 ; z(4) = 2.D0 ;
      x(5) = 539.70D0 ; y(5) = 1589.50D0 ; z(5) = 2.D0 ;


      ea(1) =  7.4d0 ; dm(1) = 0.099d0  ; cs(1) =   1.25d0 ;
      ea(2) =  9.7d0 ; dm(2) = 0.098d0  ; cs(2) =   1.25d0 ;
      ea(3) = 13.3d0 ; dm(3) = 0.363d0  ; cs(3) =   1.25d0 ;
      ea(4) = 10.0d0 ; dm(4) = 0.041d0  ; cs(4) =   1.25d0 ;
      ea(5) = 11.0d0 ; dm(5) = 0.072d0  ; cs(5) =   1.25d0 ;
      ea(6) = 21.0d0 ; dm(6) = 0.088d0  ; cs(6) = 115.00d0 ;
      ea(7) = 21.0d0 ; dm(7) = 0.0206d0 ; cs(7) =  32.00d0 ;


      t  = ene


* (ioniztion)
      do j = 1, 5

        b   =  x(j)
        u   =  y(j)
        ds  =  z(j)


      if( t .gt. b) then
        st  =  t / b
        su  =  u / b

        tp  =  t / dmc2
        bp  =  b / dmc2
        up  =  u / dmc2

        bt2  =  1.d0 - 1.d0 / ( 1.d0 + tp )**2
        bb2  =  1.d0 - 1.d0 / ( 1.d0 + bp )**2
        bu2  =  1.d0 - 1.d0 / ( 1.d0 + up )**2

        s1  =  4.D0 * pi * a0**2 * alpha**4 * ds
        s2  =  (bt2 + bb2 + bu2) * 2.D0 * bp
        s   =  s1 / s2

        f1 = 0.5d0 *
     .      ( dlog( bt2 / (1.d0 - bt2) ) - bt2 -dlog( 2.D0 * bp ) ) *
     .      ( 1.d0 - 1.d0 / st**2 )

        f2 = 1.d0 - 1.d0 / st

        f3 = - dlog( st ) / ( st + 1.d0 ) *
     .      ( 1.d0 + 2.d0 * tp ) / ( 1.d0 + tp / 2.d0 )**2

        f4 = bp**2 / ( 1.d0 + tp /2.d0 )**2 * ( st - 1.d0 ) / 2.d0

        xion(j) = s * ( f1 + f2 + f3 + f4 )

      else

        xion(j) = 0.d0

      endif
      enddo


        sum  =  0.d0
      do k = 1, 5
        sum  =  sum + xion(k)
      enddo


*(excitation)
      if( t .le. 10.d+4 ) then            ! Kai (2022/10/06)

      do n = 1, 7

      if( t .gt. ea(n) ) then

        alp  =  0.25d0 * ( t / ea(n) - 1.d0 )
        psi  =  1.d0 - dexp( -alp )

        s1   =  4.d0 * pi * a0**2 * ry / t
        s2   =  4.d0 * cs(n) * t / ry
        s3   =  dm(n) * dlog(s2)

        xexcg(n)  =  psi * s1 * s3

      else

        xexcg(n)  =  0.d0

      endif
      enddo

        ea(8)  =  9.d0
        f0c0   =  0.033d0
        om     =  1.d0
        beta   =  2.d0
        v      =  1.d0

        f1  =  4.d0 * pi * a0**2 * ( ry / ea(8) )**2
        f2  =  f0c0 * ( ea(8) / t )**om
        f3  =  ( 1.d0 - ( ea(8) / t )**beta )**v

        ss  =  f1 * f2 * f3

      if( t .gt. ea(8) ) then
        xexcg(8) = ss
      else
        xexcg(8) = 0.d0
      endif

      do k = 1, 8
        sum  =  sum + xexcg(k)
      enddo

      else  ! Kai (2022/10/06)


c Kai (2022/10/06)
      tp  = t / dmc2
      bt2 = 1.d0 - 1.d0/(1.d0 + tp)**2

      a1 = 4.d0 * pi * a0**2 / 137.d0**2 / bt2

      a2 = dlog( bt2 / ( 1.d0 - bt2) ) - bt2

      a3 = 12.3d0 + 1.26d0 * a2

      ss  =  a1 * a3
      sum =  sum + ss
c Kai (2022/10/06)


      endif ! Kai (2022/10/06)


      do k = 1, 5
        xion(k) = sum * scalion(k)
      enddo

      do k = 1, 6
        xexc(k) = sum * scalexc(k)
      enddo


      end


c====================================================
      subroutine ion_ics(ene,xion)
c====================================================
      implicit real*8(a-h,o-z)
      parameter (a0=5.29177249D-9)
      dimension xion(5)
      common / csion / dion(500),si(5,500),xip(5)
      common / nline / nion,nexc,nvib,nrot,ndea,nphn
      common / cschng / change, fit


      do j=1,5

      if(ene.lt.xip(j).or.ene.gt.dion(nion))then
         xion(j)=0.d0
         goto 1000
      else
         do i=1,nion
         if(ene.eq.dion(i))then
         xion(j) = si(j,i)
         goto 1000
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

        xion(j) = fit * xion(j)

      enddo



      end


c====================================================
      subroutine ion_dcs(itrn,bebw2nd)
c====================================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)
      parameter (pi=3.1415926535897932385d0)
      parameter (a0=5.29177249d-9)
      parameter(nqmax=10000)  ! T.Sato 2022/08/22
      dimension sdion(nqmax),sd(nqmax)  ! T.Sato 2022/08/22
      dimension wx(5),wy(5),wz(5)

      common / idumc /idum
      real*8 ran2

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)

      bebw2nd=0.0 ! initialization  ! T.Sato 2022/08/22

      alpha = 1.d0/137.035999679d0
      dmc2  = 511.003d+3


      wx(1)= 10.90d0; wy(1)=  97.39d0; wz(1)=2.d0;
      wx(2)= 13.50d0; wy(2)= 118.40d0; wz(2)=2.d0;
      wx(3)= 17.00d0; wy(3)= 122.90d0; wz(3)=2.d0;
      wx(4)= 26.30d0; wy(4)= 142.00d0; wz(4)=2.d0;
      wx(5)=553.00d0; wy(5)=1589.50d0; wz(5)=2.d0;


       ts = e(ibke+no,ipomp+1)

       qt  = 1.d0

      sdion(:)=0
      b  = wx(itrn)  ; uts = wy(itrn)  ; ds = wz(itrn)
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

      end

c====================================================
      subroutine exc_ics(ene,xexc)
c====================================================
      implicit real*8 (a-h,o-z)
      parameter (pi=3.1415926535897932385d0)
      parameter (a0=5.29177249d-9) !(cm2)
      dimension xexc(6)
      common / csexc / eexc(500),se(6,500),eip(6)
      common / nline / nion,nexc,nvib,nrot,ndea,nphn
      common / cschng / change, fit


      do j=1,6

      if(ene.lt.eip(j).or.ene.gt.eexc(nexc))then
         xexc(j)=0.d0
         goto 1000
      else
         do i=1,nexc
         if(ene.eq.eexc(i))then
         xexc(j) = se(j,i)
         goto 1000
         endif
         enddo
      endif

      do i=1,nexc-1
      if(ene.ge.eexc(i).and.ene.lt.eexc(i+1)) then
         x1 = eexc(i)
         x2 = eexc(i+1)
         y1 = se(j,i)
         y2 = se(j,i+1)
      else; endif
      enddo

         a1 = (y1-y2)/(x1-x2)
         b1 = y1 - a1 * x1
         xexc(j) = a1 * ene + b1
 1000 continue

      xexc(j) = fit * xexc(j)

      enddo

      end


c====================================================
      subroutine dea_ics( iprss, ene, xdea )
c====================================================
      implicit real*8 (a-h,o-z)

      common / csdea / edea( 53), sd(3, 53), dip(3)
      common / nline / nion, nexc, nvib, nrot, ndea, nphn

      if( iprss .eq. 1 ) idea1  =  2 ;  idea2  =  42
      if( iprss .eq. 2 ) idea1  =  3 ;  idea2  =  50
      if( iprss .eq. 3 ) idea1  =  8 ;  idea2  =  29

      if( ene .lt. edea(idea1)  .or.  ene .gt. edea(idea2) ) then
        xdea = 0.d0
        return
      endif

      do i = idea1, idea2
        if( ene .eq. edea(i) ) then
          xdea  =  sd(iprss,i)
          return
        endif
      enddo

      do i = idea1, idea2 - 1
        if( ene .ge. edea(i)  .and.  ene .lt. edea(i+1) ) then
          x1  =  edea(i)
          x2  =  edea(i+1)
          y1  =  sd(iprss,i)
          y2  =  sd(iprss,i+1)
        endif
      enddo

      a1  =  ( y1 - y2 ) / ( x1 - x2 )
      b1  =  y1 - a1 * x1

      xdea  =  a1 * ene + b1


      end


c====================================================
      subroutine vib_ics(iprss,ene,xvib)
c====================================================
      implicit real*8 (a-h,o-z)
      common / csvib / evb1( 44),sv(2, 44),vip(2),evb2( 44)
      common / nline / nion,nexc,nvib,nrot,ndea,nphn
      dimension evb(44)


      do i=1,nvib
      if(iprss.eq.1)then
         evb(i) = evb1(i)
      else
         evb(i) = evb2(i)
      endif
      enddo


      if(ene.lt.evb(1).or.ene.gt.evb(nvib))then
         xvib = 0.d0
         return
      else; endif

      do i=1,nvib
      if(ene.eq.evb(i))then
         xvib = sv(iprss,i)
         return
      else; endif
      enddo

      do i=1,nvib - 1
      if(ene.ge.evb(i).and.ene.lt.evb(i+1)) then
         x1 = evb(i)
         x2 = evb(i+1)
         y1 = sv(iprss,i)
         y2 = sv(iprss,i+1)
      else; endif
      enddo


         a1 = (y1-y2)/(x1-x2)
         b1 = y1 - a1 * x1
         xvib = a1 * ene + b1


      end


c====================================================
      subroutine phn_ics(iprss,ene,xphn)
c====================================================
      implicit real*8 (a-h,o-z)
      common / csphn / ephn(308),sp(2,308),pip(2)
      common / nline / nion,nexc,nvib,nrot,ndea,nphn


      if(ene.lt.pip(iprss).or.ene.gt.ephn(nphn))then
         xphn = 0.d0
         return
      else
         do i=1,nphn
         if(ene.eq.ephn(i))then
         xphn = sp(iprss,i)
         return
         else; endif
         enddo
      endif


      do i=1,nphn - 1
      if(ene.ge.ephn(i).and.ene.lt.ephn(i+1)) then
         x1 = ephn(i)
         x2 = ephn(i+1)
         y1 = sp(iprss,i)
         y2 = sp(iprss,i+1)
      else; endif
      enddo


         a1   = (y1-y2)/(x1-x2)
         b1   = y1 - a1 * x1
         xphn = a1 * ene + b1


      end


c====================================================
      subroutine rot_ics(iprss,ene,xrot)
c====================================================
      implicit real*8 (a-h,o-z)
      common / csrot / erot(451),sr(2,451),rip(2)
      common / nline / nion,nexc,nvib,nrot,ndea,nphn


      if(ene.lt.rip(iprss).or.ene.gt.erot(nrot))then
         xrot = 0.d0
         return
      else
         do i=1,nrot
         if(ene.eq.erot(i))then
         xrot = sr(iprss,i)
         return
         else; endif
         enddo
      endif

      do i=1,nrot - 1
      if(ene.ge.erot(i).and.ene.lt.erot(i+1)) then
         x1 = erot(i)
         x2 = erot(i+1)
         y1 = sr(iprss,i)
         y2 = sr(iprss,i+1)
      else; endif
      enddo

         a1   = (y1-y2)/(x1-x2)
         b1   = y1 - a1 * x1
         xrot = a1 * ene + b1


      end



c================================================
      subroutine water_db(ierrdb)
c================================================
      implicit real*8(a-h,o-z)
      common / csion / dion(500),si(5,500),xip(5)
      common / csexc / eexc(500),se(6,500),eip(6)
      common / csdea / edea( 53),sd(3, 53),dip(3)
      common / csvib / evb1( 44),sv(2, 44),vip(2),evb2( 44)
      common / csphn / ephn(308),sp(2,308),pip(2)
      common / csrot / erot(451),sr(2,451),rip(2)
      common / nline / nion,nexc,nvib,nrot,ndea,nphn

! T.Sato 2016/05/28, for batch.out
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

      logical   exex

      inquire( file = chfn(25)(1:ilfn(25))//'/electron/water.dat',
     & exist = exex )
      if( exex .eqv. .false. ) then
       ierrdb=1
       return
      endif

      open(25,file=chfn(25)(1:ilfn(25))//'/electron/water.dat',
     & status='old')

      read(25,*)
      read(25,*)nion,xip(1),xip(2),xip(3),xip(4),xip(5)
      do i=1,nion
      read(25,*)dion(i),si(1,i),si(2,i),si(3,i),si(4,i),si(5,i)
      enddo

      read(25,*)
      read(25,*)nexc,eip(1),eip(2),eip(3),eip(4),eip(5),eip(6)
      do i=1,nexc
      read(25,*)eexc(i),se(1,i),se(2,i),se(3,i),se(4,i),se(5,i),se(6,i)
      enddo

      read(25,*)
      read(25,*)ndea,dip(1),dip(2),dip(3)
      do i=1,ndea
      read(25,*)edea(i),sd(1,i),sd(2,i),sd(3,i)
      enddo

      read(25,*)
      read(25,*)nvib,vip(1),vip(2)
      do i=1,nvib
      read(25,*)evb1(i),sv(1,i),evb2(i),sv(2,i)
      enddo

      read(25,*)
      read(25,*)nphn,pip(1),pip(2)
      do i=1,nphn
      read(25,*)ephn(i),sp(1,i),sp(2,i)
      enddo

      read(25,*)
      read(25,*)nrot,rip(1),rip(2)
      do i=1,nrot
      read(25,*)erot(i),sr(1,i),sr(2,i)
      enddo


      close(25)


      end




c=================================================================
      function ran2(idum)
c=================================================================
      integer idum,im1,im2,imm1,ia1,ia2,iq1,iq2,ir1,ir2,ntab,ndiv
      real*8 ran2,am,eps,rnmx
      parameter (im1=2147483563,im2=2147483399,am=1./im1,imm1=im1-1,
     *ia1=40014,ia2=40692,iq1=53668,iq2=52774,ir1=12211,ir2=3791,
     *ntab=32,ndiv=1+imm1/ntab,eps=1.2e-7,rnmx=1.-eps)
      integer idum2,j,k,iv(ntab),iy
      save iv,iy,idum2
      data idum2/123456789/, iv/ntab*0/, iy/0/
      real*8 unirn
      integer irandphits

      irandphits = 1
      if(irandphits.eq.1)then
      ran2 = unirn(dummy)
      else

      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=ntab+8,1,-1
          k=idum/iq1
          idum=ia1*(idum-k*iq1)-k*ir1
          if (idum.lt.0) idum=idum+im1
          if (j.le.ntab) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/iq1
      idum=ia1*(idum-k*iq1)-k*ir1
      if (idum.lt.0) idum=idum+im1
      k=idum2/iq2
      idum2=ia2*(idum2-k*iq2)-k*ir2
      if (idum2.lt.0) idum2=idum2+im2
      j=1+iy/ndiv
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+imm1
      ran2=min(am*iy,rnmx)

      endif


      return
      end




c================================================
      subroutine etsdataup
c===============================================
      use MMBANKMOD

      implicit real*8(a-h,o-z)

      include 'param00.inc'

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
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
      common /pnsave/ egs, uus, vvs, wws, wts, tms, nms, nct(3)
!$OMP THREADPRIVATE(/pnsave/)

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
         qclusts(  8, nclsts ) = 1.d0                      ! weight change
         qclusts(  9, nclsts ) = 0.d0                      ! (ns)
         qclusts( 10, nclsts ) = 0.d0                      ! (cm)
         qclusts( 11, nclsts ) = 0.d0                      ! (cm)
         qclusts( 12, nclsts ) = 0.d0                      ! (cm)

         uus  =     u1st
         vvs  =     v1st
         wws  =     w1st
         egs  =    ec(ibkec+no,ipomp+1) * xmv
         wts  =    wt(ibkwt+no,ipomp+1)
         tms  =    tc(ibktc+no,ipomp+1)
         nms  =  name(ibknam+no,ipomp+1)

      return
      end
