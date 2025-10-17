************************************************************************
*                                                                      *
      subroutine egs_magfld(mark,markp,dist,nbeta,itmak,
     &                      a_mag,b_mag,s_mag,p_mag,t_mag,u_mag)
*                                                                      *
*                                                                      *
*       particle transport by d                                        *
*       and region check                                               *
*       modified by Nais on 2016/12/01                                 *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       dist   : distance                                              *
*       a_mag  : magnet gap(cm)                                        *
*       b_mag  : magnet field at pole tip  [kG]                        *
*       s_mag  : speicies of magnet dypole:2, quad:4, sext:6, oct:8    *
*       p_mag  : phase of magnetic field for charge particle           *
*                polarization of neuteron                              *
*       t_mag  : transform id                                          *
*       u_mag  : critical time of time dependent magnetic field        *
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
      use moddas_mesh

      implicit real*8 (a-h,o-z)

      parameter ( rlit = 29.97925d0 )

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /eparm/  esmax, esmin, emin(20)
      common /paraj/  mstz(300), parz(300)
      common /tcntl/  icntl, inucr
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx

*-----------------------------------------------------------------------
      common /ndemax/ dnmax(20)

      include 'include/egs5_h.f'
!     only destep is needed
      common/STACK/
     * egs5e(MXSTACK),      ! Total energy of particle (including rest
     * egs5x(MXSTACK),                                          ! X-po
     * egs5y(MXSTACK),                                          ! Y-po
     * egs5z(MXSTACK),                                          ! Z-po
     * egs5u(MXSTACK),                             ! X-axis direction
     * egs5v(MXSTACK),                             ! Y-axis direction
     * egs5w(MXSTACK),                             ! Z-axis direction
     * egs5uf(MXSTACK),         ! Electric field vectors of polarized
     * egs5vf(MXSTACK),
     * egs5wf(MXSTACK),
     * dnear(MXSTACK),          ! Estimated distance to nearest bounda
     * egs5wt(MXSTACK),                                    ! Particle
     * k1step(MXSTACK),      ! Scat stren to next hinge
     * k1rsd(MXSTACK),       ! Scat stren from hinge to end of step
     * k1init(MXSTACK),      ! Scat of prev hinge end of step
     * time(MXSTACK),
     * deinitial,
     * deresid,
     * denstep,
     * iq(MXSTACK),        ! Particle charge, -1(e-), 0(photons), +1(e
     * ir(MXSTACK),                                      ! Region numb
     * latch(MXSTACK),                               ! Latching variab
     * np,                                         ! Stack pointer ind
     * latchi                        ! Initialization for latch variab

c>ada.2014.12> only destep is needed ..?
!$OMP THREADPRIVATE(/STACK/)
      include 'include/egs5_misc.f'
      real*8 egs5x,egs5y,egs5z,egs5u,egs5v,egs5w,egs5uf,egs5vf,egs5wf,
     $       dnear,egs5wt
      real*8 egs5e,time
      real*8 k1step,k1rsd,k1init
      real*8 deinitial, deresid, denstep
      integer iq,ir,latch,np,latchi
      real*8 k1i,k1r,k1s

      common /egs5cmn3/thard,tinel,tmscat,hardstep,sig,scpow,dedx,
     $                 sig0,ams
!$OMP THREADPRIVATE(/egs5cmn3/)
      common /egs5cmn4/k1i,k1r,k1s
!$OMP THREADPRIVATE(/egs5cmn4/)
      real*8 pstep,dpmfp,gmfp
      common /egs5cmn6/pstep,dpmfp,gmfp
!$OMP THREADPRIVATE(/egs5cmn6/)

      include 'include/egs5_epcont.f'

      include 'include/egs5_useful.f'

      common /egs5cmn10/denstepold,denstepnew,deinit
      real*8 denstepold,denstepnew,deinit
!$OMP THREADPRIVATE(/egs5cmn10/)


*-----------------------------------------------------------------------
      common /celdg/  rhog(kvlmax)

c'Nais added lines for EM 2016/12/01

      include 'include/egs5_bounds.f'

      double precision TEMSTP
      double precision UNPSTE,VNPSTE,WNPSTE
      double precision STPEME,STPBME
      double precision elcx,elcy,elcz
      double precision bmgx,bmgy,bmgz
      double precision ein
      double precision RMSQ

      integer itre,itrm
      integer irl

*-----------------------------------------------------------------------

      common /wwin00/ iwwxyz
!$OMP THREADPRIVATE(/wwin00/)
      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww
      common /wwind0/ dvals, iwmsh, kecho, kgwwp(6), idval
      common /wwxyz/ iwxty(6), iwxnm(6), iwxrg(6),
     &               rwxmi(6), rwxma(6), rwxdl(6),
     &               iwyty(6), iwynm(6), iwyrg(6),
     &               rwymi(6), rwyma(6), rwydl(6),
     &               iwzty(6), iwznm(6), iwzrg(6),
     &               rwzmi(6), rwzma(6), rwzdl(6)

      dimension     idas(1)
      equivalence ( das, idas )
! 2025-03-03 N.SHIGYO treatment as a lost particle in case not coming out even after 1000 times.
! from Y. Nagasawa, MTC.
      data infcheck/0/
      data dprold/1.0e20/
      data zcold/1.0e20/
*-----------------------------------------------------------------------

      STPEME=0.01
      STPBME=0.01
      istcnt = 0

      RMSQ = RM*RM

      justb4=0  ! just before index, introduced by T.Sato 2023/12/29 to make step closer to the boundary

*-----------------------------------------------------------------------
*        transform ( electric field = x,  magnetic field = y )
*-----------------------------------------------------------------------

      itrm = nint( t_mag )
      call trnsuv(0.d0,1.d0,0.d0,bmgx,bmgy,bmgz,itrm)


*-----------------------------------------------------------------------
*        charge state
*-----------------------------------------------------------------------

                  chgp = ctyp

*-----------------------------------------------------------------------
*        constant
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*        dv = fpl : flight length
*-----------------------------------------------------------------------

                  iff = 0

                  dv  = dist

                  esav = e(ibke+no,ipomp+1)
                  tsav = t(ibkt+no,ipomp+1)
                  xsav = x(ibkx+no,ipomp+1)
                  ysav = y(ibky+no,ipomp+1)
                  zsav = z(ibkz+no,ipomp+1)

*-----------------------------------------------------------------------

            if( ityp .lt. 15 ) then

               emint = emin(ityp)

            else

               emint = emin(ityp) * dble( mtyp )

            end if

*-----------------------------------------------------------------------
*        delt1 : maximum flight step
*        delt  : flight step
*-----------------------------------------------------------------------

            if( mstz(111) .eq. 0 ) then
                  delt1 = min( parz(27), dist )
            end if

*-----------------------------------------------------------------------
*        delt1 : maximum flight step
*        delt  : flight step
*            parz(27) : deltm
*            parz(163): deltc for nedisp
*            dist: 1e10
*        dens = rhog(mat): density (g/cm3)
*-----------------------------------------------------------------------

            if( mstz(111) .eq. 1 ) then

                ddeltm = 0.0
                ddeltc = 0.0
               if(mat.ne.0) dens = rhog(mat) ! T.Sato 2023/03/26 avoid error when mat=0

               if( mat .eq. 0 ) dens = 1.0

                ddeltm = parz(27) / dens

                   delt1 = min( ddeltm, dist )




            end if

*-----------------------------------------------------------------------

                  delt  = delt1
	if(ityp.eq.14) delt = 1.0e30 ! T.Sato 2016/03/25, avoid longer range of photon

*-----------------------------------------------------------------------

              !ASTOM 2019/01/15
              x1 = 0
              x2 = 0
              y1 = 0
              y2 = 0
              z1 = 0
              z2 = 0
              r1 = 0
              r2 = 0
  500 continue

              iwwxyz = 0
              jwwxyz = 0
              kwwxyz = 0

              if( jtyp .ne. 0) then
                eie = e(ibke+no,ipomp+1)+RM
              else
                eie = e(ibke+no,ipomp+1)
              endif

              UNPSTE=u(ibku+no,ipomp+1)
              VNPSTE=v(ibkv+no,ipomp+1)
              WNPSTE=w(ibkw+no,ipomp+1)

! T.Sato 2020/06/03, might be used at the end of this subroutine
              uratio=1.0
              vratio=1.0
              wratio=1.0

              irl = idgr(iblz(ibkblz+no,ipomp+1))


              TEMSTP = vacdst

              s_elf = 0.0
              elcx = 0.0
              elcy = 1.0
              elcz = 0.0

              !ASTOM 2019/01/15
              if (s_mag .ge. 0) then
                    s_mgf = b_mag
                    bmgx = 0.0
                    bmgy = 1.0
                    bmgz = 0.0
                    ss_mgf = b_mag
              else
              !use magnetic map data file
                    xxc = x(ibkx+no,ipomp+1)
                    yyc = y(ibky+no,ipomp+1)
                    zzc = z(ibkz+no,ipomp+1)
                    if(s_mag .eq. -1) then
                      call readmap1(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,b_mag,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                    else if(s_mag .eq. -2) then
                       call readmap2(xxc,yyc,zzc,r1,r2,z1,z2,b_mag,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                    else if(s_mag .eq. -3) then
                        call readmap3(xxc,yyc,zzc,
     &                   x1,x2,y1,y2,z1,z2,b_mag,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                    else if(s_mag .eq. -4) then
                       call readmap4(xxc,yyc,zzc,r1,r2,z1,z2,b_mag,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                    end if
!                    s_mgf = b_mag
                    s_mgf = 1.0d0 ! T.Sato 2025/01/24, mgf parameter is duplicately considered in the case of magnetic field map
                    ss_mgf = sqrt(bmgx**2 + bmgy**2 + bmgz**2)
                    bmgx = 0.0
                    bmgy = 1.0
                    bmgz = 0.0
                  end if

              call egs_emstep(TEMSTP,EIE,UNPSTE,VNPSTE,WNPSTE,
     &                        s_elf,elcx,elcy,elcz,
     &                        ss_mgf,bmgx,bmgy,bmgz,
     &                        STPEME,STPBME,RM,RMSQ)



*-----------------------------------------------------------------------
*        mark:     description
*         -1 : outgoing to the void region
*          0 : pass the forward surface
*          2 : reflect surface
*-----------------------------------------------------------------------

            if( mstz(23) .ne. 0 ) then

               if( ( mark .eq. 0 .or. mark .eq. 2 ) .and.
     &               iff .eq. 0 ) then

                  delt = parz(28)
                  iff  = iff + 1

               else if( iff .eq. 1 ) then

                  delt = delt1
                  iff  = iff + 1

               end if

            end if

*-----------------------------------------------------------------------
*        nbeta = 1 : dv > delt,  transport to delt, dv=dv-delt
*        nbeta = 2 : dv < delt,  something happen at dv
*        nbeta = 3 : delt = rng < delt, stopped at rng, dv = dv - delt
*-----------------------------------------------------------------------

               if( dv .gt. delt ) then

                  nbeta = 1
                  dv = dv - delt

               else

                  nbeta = 2
                  delt  = dv

               end if

*-----------------------------------------------------------------------
*        distance to the boundary
*-----------------------------------------------------------------------

! T.Sato 2020/06/02, necessary to reset mark and markp
            markold=mark
            markpold=markp
            if(mstz(159).eq.0.and.mark.ne.2) mark = 1  ! basically reset, but should not reset for complex geometry (2021/01/22) -> revise 148 to 159 on 2023/03/27
            if(mstz(159).eq.0.and.mark.ne.2) markp = 0 ! basically reset, but should not reset for complex geometry (2021/01/22) -> revise 148 to 159 on 2023/03/27

 599         call gomdis(dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1),
     &                  mark,markp,nmed(ibknmd+no,ipomp+1),
     &                  iblz(ibkblz+no,ipomp+1))

            if(mark.le.-3) return ! severe error, need go back to partrs

            if(mark .le. -2) then
             if(markold.eq.-1000) goto 600 ! 2nd time failure
             mark=markold
             markp=markpold
             markold=-1000
             goto 599
            endif

*-----------------------------------------------------------------------
*        WW of xyz mesh
*-----------------------------------------------------------------------

         if( iwwdp .gt. 0 .and. iwmsh .eq. 3 ) then

               dwwt = delt
               if( dwwt .gt. dpr ) dwwt = dpr
               iwctl = 1

            call wwxdis(iwctl,dpx,dwwt,
     &                  x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1),
     &                  iwxnm(1),iwynm(1),iwznm(1),
     &                  das_iwxrg(iwxrg(1)),das_iwyrg(iwyrg(1)),
     &                  das_iwzrg(iwzrg(1)))

            if( dpx .le. dpr ) then

               if( dpx .eq. dpr ) kwwxyz = 1

               dpr = dpx
               jwwxyz = 1

            end if

         end if

*-----------------------------------------------------------------------
*           charge particle
*-----------------------------------------------------------------------

               ecc = e(ibke+no,ipomp+1)
               delt2 = delt
               ec(ibkec+no,ipomp+1) = ecc

*---------------------------------------------------------------------
cc H.Iwase 2014/1/9  (should be called as a subroutine)
*        EGS5 electron transport
*---------------------------------------------------------------------

            if( jtyp .ne. 0 .and. mat .ne. 0 .and.
     &          icntl .ne. 5 ) then


               tstep = min(delt,dpr)
               edep = 0d0
               if( mat .ne. 0 ) edep = tstep * dedx

               thard = thard  - tstep
               tinel = tinel  - tstep
               tmscat = tmscat - tstep

cc H.Iwase 2015/4/2 for keeping the old energy for [t-track]
               denstepold = denstep

               hardstep = thard  * sig
               k1s = tmscat * scpow
               denstep = tinel  * dedx

cc H.Iwase 2015/4/2 for keeping the old energy for [t-track]
               denstepnew = denstep
               deinit    = deinitial

               if(mstz(85).eq.99) then
                  write(93,*)'--------------------------------renewed'
                  write(93,'(a10,f15.9)')'(hardstep)',hardstep
                  write(93,'(a10,f15.9)')'(k1s)'     ,k1s
                  write(93,'(a10,f15.9)')'(denstep)' ,denstep
                  write(93,*)'--------------------------------'
               endif


c' void
          elseif( jtyp .ne. 0 .and. mat .eq. 0 .and.
     &          icntl .ne. 5 ) then

              delt = dpr
              if(temstp.lt.delt) then
                delt = temstp
              endif


           end if

! T.Sato 2023/12/29 particle should go to very close to the boundary before going out from the magenetic field area
       if(justb4.eq.1.and.delt.ge.dpr) then
        justb4=2
        delt=dpr*0.99d0 ! not sure how close should be
       else
        justb4=1
       endif
*-----------------------------------------------------------------------
*        next step or collision or stop
*-----------------------------------------------------------------------

         if( delt .lt. dpr ) then

            if( jtyp .ne. 0 .and. mat .ne. 0 .and.
     &          icntl .ne. 5 ) then

               ein =  e(ibke+no,ipomp+1)
               eout = e(ibke+no,ipomp+1)


               irin = idgr(iblz(ibkblz+no,ipomp+1))
               call egs5ede(irin,ityp,ein,eout,nbeta)

               ec(ibkec+no,ipomp+1) = eout
               ecc = eout

            endif

            if( nbeta .eq. 1 ) then

               goto 800

            else

! 2025-03-03 N.SHIGYO treatment as a lost particle in case not coming out even after 1000 times.
! from Y. Nagasawa, MTC.
              if(((dprold-dpr).lt.1.0e-25.and.(dprold-dpr).gt.-1.0e-25
     &        .and.dpr.lt.1.0e-3
     &        .and.(zc(ibkzc+no,ipomp+1)-zcold.lt.1.0e-14)
     &        .and.(zc(ibkzc+no,ipomp+1)-zcold.gt.-1.0e-14))
     &        .or.(delt.le.0)) then
                infcheck=infcheck+1
                dprold=dpr
                zcold=zc(ibkzc+no,ipomp+1)
                if(infcheck.gt.1000) then
                  write(6,*) '!!! Infinit loop check > 1000 !!!'
                  write(6,*)'mat,ityp,eie,ddeltm,delt1,dist,delt,
     &            dpr,nbet,infchk,zc',
     &            mat, ityp, eie, ddeltm, delt1, dist, delt, dpr, nbeta,
     &            infcheck, zc(ibkzc+no,ipomp+1)
                  mark=-2
                  infcheck=0
                  return
                endif
              else
                infcheck=0
                dprold=dpr
                zcold=zc(ibkzc+no,ipomp+1)
              endif
! to this line

               goto 700

            end if



*-----------------------------------------------------------------------
*        cross surface or stop
*-----------------------------------------------------------------------

         else if( delt .ge. dpr ) then

! 2025-03-03 N.SHIGYO treatment as a lost particle in case not coming out even after 1000 times.
! from Y. Nagasawa, MTC.
            if(((dprold-dpr).lt.1.0e-25.and.(dprold-dpr).gt.-1.0e-25
     &       .and.dpr.lt.1.0e-3
     &       .and.(zc(ibkzc+no,ipomp+1)-zcold.lt.1.0e-14)
     &       .and.(zc(ibkzc+no,ipomp+1)-zcold.gt.-1.0e-14))
     &       .or.(delt.le.0)) then
              infcheck=infcheck+1
              dprold=dpr
              zcold=zc(ibkzc+no,ipomp+1)
              if(infcheck.gt.1000) then
                write(6,*) '!!! Infinit loop check > 1000 !!!'
                write(6,*)'mat,ityp,eie,ddeltm,delt1,dist,delt,
     &           dpr,nbet,infchk,zc',
     &           mat, ityp, eie, ddeltm, delt1, dist, delt, dpr, nbeta,
     &           infcheck, zc(ibkzc+no,ipomp+1)
                mark=-2
                infcheck=0
                return
              endif
            else
              infcheck=0
              dprold=dpr
              zcold=zc(ibkzc+no,ipomp+1)
            endif
! to this line

*-----------------------------------------------------------------------
*           WW of xyz mesh
*-----------------------------------------------------------------------

            if( jwwxyz .eq. 1 .and. kwwxyz .eq. 0 ) then

               nbeta = 2
               delt = dpr
               iwwxyz = 1

               goto 700

            end if

            if( jwwxyz .eq. 1 .and. kwwxyz .eq. 1 ) iwwxyz = 2

*-----------------------------------------------------------------------
*           search new cell
*-----------------------------------------------------------------------

            irold = idgr(iblz(ibkblz+no,ipomp+1))

        uold=u(ibku+no,ipomp+1)
        vold=v(ibkv+no,ipomp+1)
        wold=w(ibkw+no,ipomp+1)

            call gomnew(1,dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                  zc(ibkzc+no,ipomp+1),
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1),
     &                  ec(ibkec+no,ipomp+1),
     &                  nmed(ibknmd+no,ipomp+1),iblz(ibkblz+no,ipomp+1),
     &                  mark,markp)

        if(mark.eq.2) then ! reflection case, remember the ratio
         if(uold.ne.0) uratio=u(ibku+no,ipomp+1)/uold
         if(vold.ne.0) vratio=v(ibkv+no,ipomp+1)/vold
         if(wold.ne.0) wratio=w(ibkw+no,ipomp+1)/wold
        endif

cc H.Iwase 2014/1/9  (photon transport: should be called as subroutine)
*---------------------------------------------------------------------
*           EGS5 photon transport
*---------------------------------------------------------------------

            if( ityp .eq. 14 ) then

               pstep = pstep - dpr
               if( mat .ne. 0 ) dpmfp = max(0.d0,dpmfp-dpr/gmfp)

               medold = medium
               irnew = idgr(iblz(ibkblz+no,ipomp+1))

               if (irnew .ne. irold) then ! Region has changed
                  irl = irnew
                  medium = med(irl)
               end if

               if( mark .eq. -1 ) pstep = 0d0

            endif

*-----------------------------------------------------------------------

               goto 600

         end if

*-----------------------------------------------------------------------
*     go back next step or tally
*-----------------------------------------------------------------------

  800 continue

               call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------
*              call tally if itstep = 1,  ncol = 15
*-----------------------------------------------------------------------


               if( itstep .ne. 0) then

c' tc = t + timd
c' time limit over : itmak = 1
                  call timtrs(itmak,mark)

                  s_elf = 0.0
                  elcx = 0.0
                  elcy = 0.0
                  elcz = 0.0
                  !ASTOM 2019/01/08
                  if (s_mag .ge. 0) then
                     call egs_getbmg(x(ibkx+no,ipomp+1),
     &                           y(ibky+no,ipomp+1),z(ibkz+no,ipomp+1),
     &                      a_mag,b_mag,s_mag,p_mag,t_mag,u_mag,
     %                      s_mgf,bmgx,bmgy,bmgz)
                  else
                    !use magnetic map data file
                    xxc = x(ibkx+no,ipomp+1)
                    yyc = y(ibky+no,ipomp+1)
                    zzc = z(ibkz+no,ipomp+1)
                    if(s_mag .eq. -1) then
                      call readmap1(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,b_mag,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                    else if(s_mag .eq. -2) then
                       call readmap2(xxc,yyc,zzc,r1,r2,z1,z2,b_mag,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                    else if(s_mag .eq. -3) then
                        call readmap3(xxc,yyc,zzc,
     &                   x1,x2,y1,y2,z1,z2,b_mag,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                    else if(s_mag .eq. -4) then
                       call readmap4(xxc,yyc,zzc,r1,r2,z1,z2,b_mag,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                    end if
                  end if

                  call egs_emmov(EIE,chgp,
     &                x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                z(ibkz+no,ipomp+1),
     &                xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                zc(ibkzc+no,ipomp+1),
     &                UNPSTE,VNPSTE,WNPSTE,
     &                u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                w(ibkw+no,ipomp+1),
     &                s_elf,elcx,elcy,elcz,
     &                s_mgf,bmgx,bmgy,bmgz,
     &                RM,RMSQ,mark)

                  if( itmak .eq. 1 ) return

                  esav = ec(ibkec+no,ipomp+1)
                  tsav = tc(ibktc+no,ipomp+1)
                  xsav = xc(ibkxc+no,ipomp+1)
                  ysav = yc(ibkyc+no,ipomp+1)
                  zsav = zc(ibkzc+no,ipomp+1)

                   ncol = 15

                   call analyz(ncol,mark)

                else

                  s_elf = 0.0
                  elcx = 0.0
                  elcy = 0.0
                  elcz = 0.0
                  !ASTOM 2019/01/08
                  if (s_mag .ge. 0) then
                     call egs_getbmg(x(ibkx+no,ipomp+1),
     &                               y(ibky+no,ipomp+1),
     &                               z(ibkz+no,ipomp+1),
     &                      a_mag,b_mag,s_mag,p_mag,t_mag,u_mag,
     &                      s_mgf,bmgx,bmgy,bmgz)
                  else
                    !use magnetic map data file
                    xxc = x(ibkx+no,ipomp+1)
                    yyc = y(ibky+no,ipomp+1)
                    zzc = z(ibkz+no,ipomp+1)
                    if(s_mag .eq. -1) then
                      call readmap1(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,b_mag,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                    else if(s_mag .eq. -2) then
                       call readmap2(xxc,yyc,zzc,r1,r2,z1,z2,b_mag,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                    else if(s_mag .eq. -3) then
                        call readmap3(xxc,yyc,zzc,
     &                   x1,x2,y1,y2,z1,z2,b_mag,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                    else if(s_mag .eq. -4) then
                       call readmap4(xxc,yyc,zzc,r1,r2,z1,z2,b_mag,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                    end if
                  end if

                  call egs_emmov(EIE,chgp,
     &                x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                z(ibkz+no,ipomp+1),
     &                xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                zc(ibkzc+no,ipomp+1),
     &                UNPSTE,VNPSTE,WNPSTE,
     &                u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                w(ibkw+no,ipomp+1),
     &                s_elf,elcx,elcy,elcz,
     &                s_mgf,bmgx,bmgy,bmgz,
     &                RM,RMSQ,mark)

                  end if

*-----------------------------------------------------------------------

                  e(ibke+no,ipomp+1) = ec(ibkec+no,ipomp+1)
                  t(ibkt+no,ipomp+1) = tc(ibktc+no,ipomp+1)
                  x(ibkx+no,ipomp+1) = xc(ibkxc+no,ipomp+1)
                  y(ibky+no,ipomp+1) = yc(ibkyc+no,ipomp+1)
                  z(ibkz+no,ipomp+1) = zc(ibkzc+no,ipomp+1)

                  goto 500

*-----------------------------------------------------------------------
*     update coordinate and momentum
*-----------------------------------------------------------------------

  700 continue

                  call gomupr(mark,markp,delt)


*-----------------------------------------------------------------------

  600 continue

                  esav1 = e(ibke+no,ipomp+1)
                  tsav1 = t(ibkt+no,ipomp+1)
                  xsav1 = x(ibkx+no,ipomp+1)
                  ysav1 = y(ibky+no,ipomp+1)
                  zsav1 = z(ibkz+no,ipomp+1)

c' tc = t + timd
c' time limit over : itmak = 1

                  e(ibke+no,ipomp+1) = esav
                  t(ibkt+no,ipomp+1) = tsav
                  x(ibkx+no,ipomp+1) = xsav
                  y(ibky+no,ipomp+1) = ysav
                  z(ibkz+no,ipomp+1) = zsav

                  call timtrs(itmak,mark)

                  e(ibke+no,ipomp+1) = esav1
                  t(ibkt+no,ipomp+1) = tsav1
                  x(ibkx+no,ipomp+1) = xsav1
                  y(ibky+no,ipomp+1) = ysav1
                  z(ibkz+no,ipomp+1) = zsav1


                  s_elf = 0.0
                  elcx = 0.0
                  elcy = 0.0
                  elcz = 0.0
                  !ASTOM 2019/01/08
                  if (s_mag .ge. 0) then
                     call egs_getbmg(x(ibkx+no,ipomp+1),
     &                               y(ibky+no,ipomp+1),
     &                               z(ibkz+no,ipomp+1),
     &                      a_mag,b_mag,s_mag,p_mag,t_mag,u_mag,
     %                      s_mgf,bmgx,bmgy,bmgz)
                  else
                    !use magnetic map data file
                    xxc = x(ibkx+no,ipomp+1)
                    yyc = y(ibky+no,ipomp+1)
                    zzc = z(ibkz+no,ipomp+1)
                    if(s_mag .eq. -1) then
                      call readmap1(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,b_mag,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                    else if(s_mag .eq. -2) then
                       call readmap2(xxc,yyc,zzc,r1,r2,z1,z2,b_mag,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                    else if(s_mag .eq. -3) then
                        call readmap3(xxc,yyc,zzc,
     &                   x1,x2,y1,y2,z1,z2,b_mag,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                    else if(s_mag .eq. -4) then
                       call readmap4(xxc,yyc,zzc,r1,r2,z1,z2,b_mag,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                    end if
                  end if

                  call egs_emmov(EIE,chgp,
     &                x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                z(ibkz+no,ipomp+1),
     &                xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                zc(ibkzc+no,ipomp+1),
     &                UNPSTE,VNPSTE,WNPSTE,
     &                u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                w(ibkw+no,ipomp+1),
     &                s_elf,elcx,elcy,elcz,
     &                s_mgf,bmgx,bmgy,bmgz,
     &                RM,RMSQ,mark)

                  e(ibke+no,ipomp+1) = esav
                  t(ibkt+no,ipomp+1) = tsav
                  x(ibkx+no,ipomp+1) = xsav
                  y(ibky+no,ipomp+1) = ysav
                  z(ibkz+no,ipomp+1) = zsav

        if(mark.eq.2) then ! reflection case, change the direction again, T.Sato 2020/06/03
         u(ibku+no,ipomp+1)=u(ibku+no,ipomp+1)*uratio
         v(ibkv+no,ipomp+1)=v(ibkv+no,ipomp+1)*vratio
         w(ibkw+no,ipomp+1)=w(ibkw+no,ipomp+1)*wratio
         tmp=sqrt(u(ibku+no,ipomp+1)**2+v(ibkv+no,ipomp+1)**2+
     &                                 w(ibkw+no,ipomp+1)**2)
         u(ibku+no,ipomp+1)=u(ibku+no,ipomp+1)/tmp
         v(ibkv+no,ipomp+1)=v(ibkv+no,ipomp+1)/tmp
         w(ibkw+no,ipomp+1)=w(ibkw+no,ipomp+1)/tmp
        endif

      return
      end

************************************************************************
*                                                                      *
      subroutine egs_elmgfd(mark,markp,dist,nbeta,itmak,
     &                      s_elf,s_mgf,t_elf,t_mgf,e_chs,
     &                      emap_type, mmap_type)
*                                                                      *
*                                                                      *
*       particle transport by d                                        *
*       and region check                                               *
*       modified by Nais on 2016/12/01                                 *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       dist   : distance                                              *
*       s_elf  : electric field (kV/cm)                                *
*       s_mgf  : dipole magnet field  [kG]                             *
*       t_elf  : transform for electirc field                          *
*       t_mgf  : transform for magnetic field                          *
*       e_chs  : charge state for this region                          *
*       emap_type : type of electric data map                          *
*                   -1:xyz-list                                        *
*                   -2:rz-list                                         *
*                   -3:xyz-map                                         *
*                   -4:rz-map                                          *
*       mmap_type : type of magnetic data map                          *
*                   -1:xyz-list                                        *
*                   -2:rz-list                                         *
*                   -3:xyz-map                                         *
*                   -4:rz-map                                          *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       mark   : out put code of geom                                  *
*       markp  : =1 already check the cell                             *
*       nbeta  : =1,2; reactions or cross 3; stopped                   *
*       itmak  : 0, normal, 1, out of time range                       *
*                                                                      *
* AdvanceSoft Hasemi 2019/12/23 modified for electric map                             *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use moddas_mesh

      implicit real*8 (a-h,o-z)

      parameter ( rlit = 29.97925d0 )

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /eparm/  esmax, esmin, emin(20)
      common /paraj/  mstz(300), parz(300)
      common /tcntl/  icntl, inucr
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx

*-----------------------------------------------------------------------
      common /ndemax/ dnmax(20)

      include 'include/egs5_h.f'
!     only destep is needed
      common/STACK/
     * egs5e(MXSTACK),      ! Total energy of particle (including rest
     * egs5x(MXSTACK),                                          ! X-po
     * egs5y(MXSTACK),                                          ! Y-po
     * egs5z(MXSTACK),                                          ! Z-po
     * egs5u(MXSTACK),                             ! X-axis direction
     * egs5v(MXSTACK),                             ! Y-axis direction
     * egs5w(MXSTACK),                             ! Z-axis direction
     * egs5uf(MXSTACK),         ! Electric field vectors of polarized
     * egs5vf(MXSTACK),
     * egs5wf(MXSTACK),
     * dnear(MXSTACK),          ! Estimated distance to nearest bounda
     * egs5wt(MXSTACK),                                    ! Particle
     * k1step(MXSTACK),      ! Scat stren to next hinge
     * k1rsd(MXSTACK),       ! Scat stren from hinge to end of step
     * k1init(MXSTACK),      ! Scat of prev hinge end of step
     * time(MXSTACK),
     * deinitial,
     * deresid,
     * denstep,
     * iq(MXSTACK),        ! Particle charge, -1(e-), 0(photons), +1(e
     * ir(MXSTACK),                                      ! Region numb
     * latch(MXSTACK),                               ! Latching variab
     * np,                                         ! Stack pointer ind
     * latchi                        ! Initialization for latch variab

c>ada.2014.12> only destep is needed ..?
!$OMP THREADPRIVATE(/STACK/)
      include 'include/egs5_misc.f'
      real*8 egs5x,egs5y,egs5z,egs5u,egs5v,egs5w,egs5uf,egs5vf,egs5wf,
     $       dnear,egs5wt
      real*8 egs5e,time
      real*8 k1step,k1rsd,k1init
      real*8 deinitial, deresid, denstep
      integer iq,ir,latch,np,latchi
      real*8 k1i,k1r,k1s

      common /egs5cmn3/thard,tinel,tmscat,hardstep,sig,scpow,dedx,
     $                 sig0,ams
!$OMP THREADPRIVATE(/egs5cmn3/)
      common /egs5cmn4/k1i,k1r,k1s
!$OMP THREADPRIVATE(/egs5cmn4/)
      real*8 pstep,dpmfp,gmfp
      common /egs5cmn6/pstep,dpmfp,gmfp
!$OMP THREADPRIVATE(/egs5cmn6/)

      include 'include/egs5_epcont.f'

      include 'include/egs5_useful.f'

      common /egs5cmn10/denstepold,denstepnew,deinit
      real*8 denstepold,denstepnew,deinit
!$OMP THREADPRIVATE(/egs5cmn10/)


*-----------------------------------------------------------------------
      common /celdg/  rhog(kvlmax)

c'Nais added lines for EM 2016/12/01

      include 'include/egs5_bounds.f'

      double precision TEMSTP
      double precision UNPSTE,VNPSTE,WNPSTE
      double precision EMOVEE
      double precision BFLTMP
      double precision STPEME,STPBME
      double precision elcx,elcy,elcz
      double precision bmgx,bmgy,bmgz
      double precision ein,engup
      double precision RMSQ

      double precision sumemeg

      integer itre,itrm
      integer irl

*-----------------------------------------------------------------------

      common /wwin00/ iwwxyz
!$OMP THREADPRIVATE(/wwin00/)
      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww
      common /wwind0/ dvals, iwmsh, kecho, kgwwp(6), idval
      common /wwxyz/ iwxty(6), iwxnm(6), iwxrg(6),
     &               rwxmi(6), rwxma(6), rwxdl(6),
     &               iwyty(6), iwynm(6), iwyrg(6),
     &               rwymi(6), rwyma(6), rwydl(6),
     &               iwzty(6), iwznm(6), iwzrg(6),
     &               rwzmi(6), rwzma(6), rwzdl(6)

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      STPEME=0.01
      STPBME=0.01

      RMSQ = RM*RM

      sumemeg = 0.0d0



*-----------------------------------------------------------------------
*        transform ( electric field = x,  magnetic field = y )
*-----------------------------------------------------------------------
      itre = nint( t_elf )
      call trnsuv(1.d0,0.d0,0.d0,elcx,elcy,elcz,itre)

      itrm = nint( t_mgf )
      call trnsuv(0.d0,1.d0,0.d0,bmgx,bmgy,bmgz,itrm)


*-----------------------------------------------------------------------
*        charge state
*-----------------------------------------------------------------------

                  chgp = ctyp

*-----------------------------------------------------------------------
*        constant
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*        dv = fpl : flight length
*-----------------------------------------------------------------------

                  iff = 0

                  dv  = dist

                  esav = e(ibke+no,ipomp+1)
                  tsav = t(ibkt+no,ipomp+1)
                  xsav = x(ibkx+no,ipomp+1)
                  ysav = y(ibky+no,ipomp+1)
                  zsav = z(ibkz+no,ipomp+1)

*-----------------------------------------------------------------------

            if( ityp .lt. 15 ) then

               emint = emin(ityp)

            else

               emint = emin(ityp) * dble( mtyp )

            end if

*-----------------------------------------------------------------------
*        delt1 : maximum flight step
*        delt  : flight step
*-----------------------------------------------------------------------

            if( mstz(111) .eq. 0 ) then
                  delt1 = min( parz(27), dist )
            end if

*-----------------------------------------------------------------------
*        delt1 : maximum flight step
*        delt  : flight step
*            parz(27) : deltm
*            parz(163): deltc for nedisp
*            dist: 1e10
*        dens = rhog(mat): density (g/cm3)
*-----------------------------------------------------------------------

            if( mstz(111) .eq. 1 ) then

                ddeltm = 0.0
                ddeltc = 0.0
               if( mat .ne .0 ) dens = rhog(mat) ! S.Abe 2024/11/27 avoid error when mat=0

               if( mat .eq. 0 ) dens = 1.0

                ddeltm = parz(27) / dens

                   delt1 = min( ddeltm, dist )




            end if

*-----------------------------------------------------------------------

                  delt  = delt1
	if(ityp.eq.14) delt = 1.0e30 ! T.Sato 2016/03/25, avoid longer range of photon

*-----------------------------------------------------------------------

              x1 = 0
              x2 = 0
              y1 = 0
              y2 = 0
              z1 = 0
              z2 = 0
              r1 = 0
              r2 = 0
  500 continue

              iwwxyz = 0
              jwwxyz = 0
              kwwxyz = 0



              if( jtyp .ne. 0) then
                eie = e(ibke+no,ipomp+1)+RM
              else
                eie = e(ibke+no,ipomp+1)
              endif

              UNPSTE=u(ibku+no,ipomp+1)
              VNPSTE=v(ibkv+no,ipomp+1)
              WNPSTE=w(ibkw+no,ipomp+1)

! T.Sato 2020/06/03, might be used at the end of this subroutine
              uratio=1.0
              vratio=1.0
              wratio=1.0

              irl = idgr(iblz(ibkblz+no,ipomp+1))


              TEMSTP = vacdst


              elcx0 = elcx
              elcy0 = elcy
              elcz0 = elcz

              if(emap_type .ge. 0) then
                ss_elf = s_elf
              else
                ! use electric field map data
                xxc = x(ibkx+no,ipomp+1)
                yyc = y(ibky+no,ipomp+1)
                zzc = z(ibkz+no,ipomp+1)
                if(emap_type .eq. -1) then
                  call readelcmap1(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,s_elf,
     &                   elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                else if(emap_type .eq. -2) then
                  call readelcmap2(xxc,yyc,zzc,r1,r2,z1,z2,s_elf,
     &                   elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                else if(emap_type .eq. -3) then
                  call readelcmap3(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,s_elf,
     &                   elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                else if(emap_type .eq. -4) then
                  call readelcmap4(xxc,yyc,zzc,r1,r2,z1,z2,s_elf,
     &                   elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                endif
                  ss_elf = sqrt(elcx**2 + elcy**2 + elcz**2)
                  ss_elf = SIGN(ss_elf, s_elf)
                  elcx1 = elcx/ss_elf
                  elcy1 = elcy/ss_elf
                  elcz1 = elcz/ss_elf
                  call trnsuv(elcx1,elcy1,elcz1,elcx,elcy,elcz,itre)
            endif

            bmgx0 = bmgx
            bmgy0 = bmgy
            bmgz0 = bmgz

            if(mmap_type .ge. 0) then
               ss_mgf = s_mgf
            else
               ! use magnetic field map data
                xxc = x(ibkx+no,ipomp+1)
                yyc = y(ibky+no,ipomp+1)
                zzc = z(ibkz+no,ipomp+1)
                if(mmap_type .eq. -1) then
                  call readmap1(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,s_mgf,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                else if(mmap_type .eq. -2) then
                  call readmap2(xxc,yyc,zzc,r1,r2,z1,z2,s_mgf,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                else if(mmap_type .eq. -3) then
                  call readmap3(xxc,yyc,zzc,x1,x2,y1,y2,z1,z2,s_mgf,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                else if(mmap_type .eq. -4) then
                  call readmap4(xxc,yyc,zzc,r1,r2,z1,z2,s_mgf,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                endif
                  ss_mgf = sqrt(bmgx**2 + bmgy**2 + bmgz**2)
                  ss_mgf = SIGN(ss_mgf, s_mgf)
                  bmgx1 = bmgx/ss_mgf
                  bmgy1 = bmgy/ss_mgf
                  bmgz1 = bmgz/ss_mgf
                  call trnsuv(bmgx1,bmgy1,bmgz1,bmgx,bmgy,bmgz,itrm)
            endif

              call egs_emstep(TEMSTP,EIE,UNPSTE,VNPSTE,WNPSTE,
     &                        ss_elf,elcx,elcy,elcz,
     &                        ss_mgf,bmgx,bmgy,bmgz,
     &                        STPEME,STPBME,RM,RMSQ)



*-----------------------------------------------------------------------
*        mark:     description
*         -1 : outgoing to the void region
*          0 : pass the forward surface
*          2 : reflect surface
*-----------------------------------------------------------------------

            if( mstz(23) .ne. 0 ) then

               if( ( mark .eq. 0 .or. mark .eq. 2 ) .and.
     &               iff .eq. 0 ) then

                  delt = parz(28)
                  iff  = iff + 1

               else if( iff .eq. 1 ) then

                  delt = delt1
                  iff  = iff + 1

               end if

            end if

              if(temstp.lt.delt) then
                delt = temstp
              endif

*-----------------------------------------------------------------------
*        nbeta = 1 : dv > delt,  transport to delt, dv=dv-delt
*        nbeta = 2 : dv < delt,  something happen at dv
*        nbeta = 3 : delt = rng < delt, stopped at rng, dv = dv - delt
*-----------------------------------------------------------------------

               if( dv .gt. delt ) then

                  nbeta = 1
                  dv = dv - delt

               else

                  nbeta = 2
                  delt  = dv

               end if



*-----------------------------------------------------------------------
*        distance to the boundary
*-----------------------------------------------------------------------

! T.Sato 2020/06/02, necessary to reset mark and markp
            markold=mark
            markpold=markp
            iteration=0
            if(mark.ne.2) mark = 1  ! basically reset, but should not reset for complext geometry (2021/01/22)
            if(mark.ne.2) markp = 0 ! basically reset, but should not reset for complext geometry (2021/01/22)


 599        call gomdis(dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                      z(ibkz+no,ipomp+1),
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1),
     &                  mark,markp,nmed(ibknmd+no,ipomp+1),
     &                  iblz(ibkblz+no,ipomp+1))

            if(mark.le.-3) then ! T.Sato 2021/01/26, recovery process
             iteration=iteration+1  ! most particle is recovered with iteration=1, but search both direction may be necessary
             if(iteration.ge.1000) then
             return ! severe error, need go back to partrs
             endif
             if(mod(iteration,2).eq.1) then
              posinega=1.0d0
             else
              posinega=-1.0d0
             endif
             x(ibkx+no,ipomp+1)=x(ibkx+no,ipomp+1)+parz(28)*
     &       u(ibku+no,ipomp+1)*iteration*posinega
             y(ibky+no,ipomp+1)=y(ibky+no,ipomp+1)+parz(28)*
     &       v(ibkv+no,ipomp+1)*iteration*posinega
             z(ibkz+no,ipomp+1)=z(ibkz+no,ipomp+1)+parz(28)*
     &       w(ibkw+no,ipomp+1)*iteration*posinega
             mark=0   ! index for boundary crossing
             markp=0  ! check region
             goto 599
            endif

            if(mark .le. -2) then
             if(markold.eq.-1000) goto 600 ! 2nd time failure
             mark=markold
             markp=markpold
             markold=-1000
             goto 599
            endif

*-----------------------------------------------------------------------
*        WW of xyz mesh
*-----------------------------------------------------------------------

         if( iwwdp .gt. 0 .and. iwmsh .eq. 3 ) then

               dwwt = delt
               if( dwwt .gt. dpr ) dwwt = dpr
               iwctl = 1

            call wwxdis(iwctl,dpx,dwwt,
     &                  x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1),
     &                  iwxnm(1),iwynm(1),iwznm(1),
     &                  das_iwxrg(iwxrg(1)),das_iwyrg(iwyrg(1)),
     &                  das_iwzrg(iwzrg(1)))

            if( dpx .le. dpr ) then

               if( dpx .eq. dpr ) kwwxyz = 1

               dpr = dpx
               jwwxyz = 1

            end if

         end if

*-----------------------------------------------------------------------
*           charge particle
*-----------------------------------------------------------------------

               ecc = e(ibke+no,ipomp+1)
               delt2 = delt
               ec(ibkec+no,ipomp+1) = ecc

*---------------------------------------------------------------------
cc H.Iwase 2014/1/9  (should be called as a subroutine)
*        EGS5 electron transport
*---------------------------------------------------------------------

            if( jtyp .ne. 0 .and. mat .ne. 0 .and.
     &          icntl .ne. 5 ) then


               tstep = min(delt,dpr)
               edep = 0d0
               if( mat .ne. 0 ) edep = tstep * dedx

               thard = max(0.d0,thard  - tstep)
               tinel = max(0.d0,tinel  - tstep)
               tmscat = max(0.d0,tmscat - tstep)

cc H.Iwase 2015/4/2 for keeping the old energy for [t-track]
               denstepold = denstep

               hardstep = thard  * sig
               k1s = tmscat * scpow
               denstep = tinel  * dedx

cc H.Iwase 2015/4/2 for keeping the old energy for [t-track]
               denstepnew = denstep
               deinit    = deinitial

               if(mstz(85).eq.99) then
                  write(93,*)'--------------------------------renewed'
                  write(93,'(a10,f15.9)')'(hardstep)',hardstep
                  write(93,'(a10,f15.9)')'(k1s)'     ,k1s
                  write(93,'(a10,f15.9)')'(denstep)' ,denstep
                  write(93,*)'--------------------------------'
               endif


c' void
          elseif( jtyp .ne. 0 .and. mat .eq. 0 .and.
     &          icntl .ne. 5 ) then

              delt = dpr
              if(temstp.lt.delt) then
                delt = temstp
              endif


           end if

*-----------------------------------------------------------------------
*        next step or collision or stop
*-----------------------------------------------------------------------

         if( delt .lt. dpr ) then

            if( jtyp .ne. 0 .and. mat .ne. 0 .and.
     &          icntl .ne. 5 ) then

              call egs_geteng(EMOVEE,chgp,
     &                x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                z(ibkz+no,ipomp+1),
     &                x(ibkx+no,ipomp+1)+delt*UNPSTE,
     &                y(ibky+no,ipomp+1)+delt*VNPSTE,
     &                z(ibkz+no,ipomp+1)+delt*WNPSTE,
     &                ss_elf,elcx,elcy,elcz)


               ein =  e(ibke+no,ipomp+1) + EMOVEE
               eout = e(ibke+no,ipomp+1) + EMOVEE


               irin = idgr(iblz(ibkblz+no,ipomp+1))
               call egs5ede(irin,ityp,ein,eout,nbeta)

               ec(ibkec+no,ipomp+1) = eout
               ecc = eout

c
               sumemeg = sumemeg + EMOVEE
c


            elseif( jtyp .ne. 0 .and. mat .eq. 0 .and.
     &          icntl .ne. 5 ) then

              call egs_geteng(EMOVEE,chgp,
     &                x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                z(ibkz+no,ipomp+1),
     &                x(ibkx+no,ipomp+1)+delt*UNPSTE,
     &                y(ibky+no,ipomp+1)+delt*VNPSTE,
     &                z(ibkz+no,ipomp+1)+delt*WNPSTE,
     &                ss_elf,elcx,elcy,elcz)

                IF ((EIE + EMOVEE.LT.RM)) THEN
                  ec(ibkec+no,ipomp+1) = 0.0d0
                  nbeta = 3
                else
                  ec(ibkec+no,ipomp+1) = EIE + EMOVEE - RM
                endif

            endif

*            write(*,*) 'ec,nbeta=',ec(ibkec+no),nbeta

            if( nbeta .eq. 1 ) then

               goto 800

            else

               goto 700

            end if



*-----------------------------------------------------------------------
*        cross surface or stop
*-----------------------------------------------------------------------

         else if( delt .ge. dpr ) then

*-----------------------------------------------------------------------
*           WW of xyz mesh
*-----------------------------------------------------------------------

            if( jwwxyz .eq. 1 .and. kwwxyz .eq. 0 ) then

                  nbeta = 2
                  delt = dpr
                  iwwxyz = 1

                  goto 700

            end if

            if( jwwxyz .eq. 1 .and. kwwxyz .eq. 1 ) iwwxyz = 2

*-----------------------------------------------------------------------
*           search new cell
*-----------------------------------------------------------------------

        irold = idgr(iblz(ibkblz+no,ipomp+1))
        uold=u(ibku+no,ipomp+1)
        vold=v(ibkv+no,ipomp+1)
        wold=w(ibkw+no,ipomp+1)
            call gomnew(1,dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                  zc(ibkzc+no,ipomp+1),
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1),
     &                  ec(ibkec+no,ipomp+1),
     &                  nmed(ibknmd+no,ipomp+1),iblz(ibkblz+no,ipomp+1),
     &                  mark,markp)
                  engup = EIE

        if(mark.eq.2) then ! reflection case, remember the ratio
         if(uold.ne.0) uratio=u(ibku+no,ipomp+1)/uold
         if(vold.ne.0) vratio=v(ibkv+no,ipomp+1)/vold
         if(wold.ne.0) wratio=w(ibkw+no,ipomp+1)/wold
        endif
                  call egs_upeng(engup,chgp,
     &                x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                z(ibkz+no,ipomp+1),
     &                xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                zc(ibkzc+no,ipomp+1),
     &                ss_elf,elcx,elcy,elcz,
     &                RM,nbeta)

                  if( jtyp .ne. 0) then
                      if(engup.LT.RM) then
                         ec(ibkec+no,ipomp+1) = 0.0
                         nbeta=3
                      else
                         ec(ibkec+no,ipomp+1) = engup - RM
                      endif
                  endif


cc H.Iwase 2014/1/9  (photon transport: should be called as subroutine)
*---------------------------------------------------------------------
*           EGS5 photon transport
*---------------------------------------------------------------------

            if( ityp .eq. 14 ) then

               pstep = pstep - dpr
               if( mat .ne. 0 ) dpmfp = max(0.d0,dpmfp-dpr/gmfp)

               medold = medium
               irnew = idgr(iblz(ibkblz+no,ipomp+1))

               if (irnew .ne. irold) then ! Region has changed
                  irl = irnew
                  medium = med(irl)
               end if

               if( mark .eq. -1 ) pstep = 0d0

            endif

*-----------------------------------------------------------------------


               goto 600

         end if

*-----------------------------------------------------------------------
*     go back next step or tally
*-----------------------------------------------------------------------

  800 continue

               call gomupr(mark,markp,delt)

*-----------------------------------------------------------------------
*              call tally if itstep = 1,  ncol = 15
*-----------------------------------------------------------------------

                  elcx = elcx0
                  elcy = elcy0
                  elcz = elcz0
                  bmgx = bmgx0
                  bmgy = bmgy0
                  bmgz = bmgz0

               if( itstep .ne. 0) then

c' tc = t + timd
c' time limit over : itmak = 1
                  call timtrs(itmak,mark)

                  if(emap_type .ge. 0) then
                    ss_elf = s_elf
                  else
                    ! use electric field map data
                    xxc = x(ibkx+no,ipomp+1)
                    yyc = y(ibky+no,ipomp+1)
                    zzc = z(ibkz+no,ipomp+1)
                    if(emap_type .eq. -1) then
                      call readelcmap1(xxc,yyc,zzc,
     &                   x1,x2,y1,y2,z1,z2,s_elf,
     &                   elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                    else if(emap_type .eq. -2) then
                      call readelcmap2(xxc,yyc,zzc,r1,r2,z1,z2,s_elf,
     &                   elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                    else if(emap_type .eq. -3) then
                      call readelcmap3(xxc,yyc,zzc,
     &                   x1,x2,y1,y2,z1,z2,s_elf,
     &                   elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                    else if(emap_type .eq. -4) then
                      call readelcmap4(xxc,yyc,zzc,r1,r2,z1,z2,s_elf,
     &                   elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                    endif
                     ss_elf = sqrt(elcx**2 + elcy**2 + elcz**2)
                     ss_elf = SIGN(ss_elf, s_elf)
                     elcx1 = elcx/ss_elf
                     elcy1 = elcy/ss_elf
                     elcz1 = elcz/ss_elf
                     call trnsuv(elcx1,elcy1,elcz1,elcx,elcy,elcz,itre)
                  endif

                  if(mmap_type .ge. 0) then
                    ss_mgf = s_mgf
                  else
                  ! use magnetic field map data
                    xxc = x(ibkx+no,ipomp+1)
                    yyc = y(ibky+no,ipomp+1)
                    zzc = z(ibkz+no,ipomp+1)
                    if(mmap_type .eq. -1) then
                     call readmap1(xxc,yyc,zzc,
     &                   x1,x2,y1,y2,z1,z2,s_mgf,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                    else if(mmap_type .eq. -2) then
                     call readmap2(xxc,yyc,zzc,r1,r2,z1,z2,s_mgf,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                    else if(mmap_type .eq. -3) then
                     call readmap3(xxc,yyc,zzc,
     &                   x1,x2,y1,y2,z1,z2,s_mgf,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                    else if(mmap_type .eq. -4) then
                     call readmap4(xxc,yyc,zzc,r1,r2,z1,z2,s_mgf,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                    endif
                     ss_mgf = sqrt(bmgx**2 + bmgy**2 + bmgz**2)
                     ss_mgf = SIGN(ss_mgf, s_mgf)
                     bmgx1 = bmgx/ss_mgf
                     bmgy1 = bmgy/ss_mgf
                     bmgz1 = bmgz/ss_mgf
                     call trnsuv(bmgx1,bmgy1,bmgz1,bmgx,bmgy,bmgz,itrm)
                   endif

                    call egs_emmov(EIE,chgp,
     &                x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                z(ibkz+no,ipomp+1),
     &                xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                zc(ibkzc+no,ipomp+1),
     &                UNPSTE,VNPSTE,WNPSTE,
     &                u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                w(ibkw+no,ipomp+1),
     &                ss_elf,elcx,elcy,elcz,
     &                ss_mgf,bmgx,bmgy,bmgz,
     &                RM,RMSQ,mark)

                  if( itmak .eq. 1 ) return

                  esav = ec(ibkec+no,ipomp+1)
                  tsav = tc(ibktc+no,ipomp+1)
                  xsav = xc(ibkxc+no,ipomp+1)
                  ysav = yc(ibkyc+no,ipomp+1)
                  zsav = zc(ibkzc+no,ipomp+1)

                   ncol = 15

                   call analyz(ncol,mark)

                else

                  if(emap_type .ge. 0) then
                    ss_elf = s_elf
                  else
                    ! use electric field map data
                    xxc = x(ibkx+no,ipomp+1)
                    yyc = y(ibky+no,ipomp+1)
                    zzc = z(ibkz+no,ipomp+1)
                    if(emap_type .eq. -1) then
                      call readelcmap1(xxc,yyc,zzc,
     &                   x1,x2,y1,y2,z1,z2,s_elf,
     &                   elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                    else if(emap_type .eq. -2) then
                      call readelcmap2(xxc,yyc,zzc,r1,r2,z1,z2,s_elf,
     &                   elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                    else if(emap_type .eq. -3) then
                      call readelcmap3(xxc,yyc,zzc,
     &                   x1,x2,y1,y2,z1,z2,s_elf,
     &                   elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                    else if(emap_type .eq. -4) then
                      call readelcmap4(xxc,yyc,zzc,r1,r2,z1,z2,s_elf,
     &                   elcx,elcy,elcz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                    endif
                     ss_elf = sqrt(elcx**2 + elcy**2 + elcz**2)
                     ss_elf = SIGN(ss_elf, s_elf)
                     elcx1 = elcx/ss_elf
                     elcy1 = elcy/ss_elf
                     elcz1 = elcz/ss_elf
                     call trnsuv(elcx1,elcy1,elcz1,elcx,elcy,elcz,itre)
                  endif

                  if(mmap_type .ge. 0) then
                    ss_mgf = s_mgf
                  else
                  ! use magnetic field map data
                    xxc = x(ibkx+no,ipomp+1)
                    yyc = y(ibky+no,ipomp+1)
                    zzc = z(ibkz+no,ipomp+1)
                    if(mmap_type .eq. -1) then
                     call readmap1(xxc,yyc,zzc,
     &                   x1,x2,y1,y2,z1,z2,s_mgf,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                    else if(mmap_type .eq. -2) then
                     call readmap2(xxc,yyc,zzc,r1,r2,z1,z2,s_mgf,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                    else if(mmap_type .eq. -3) then
                     call readmap3(xxc,yyc,zzc,
     &                   x1,x2,y1,y2,z1,z2,s_mgf,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   bxx1y1z1, byx1y1z1, bzx1y1z1,
     &                   bxx1y1z2, byx1y1z2, bzx1y1z2,
     &                   bxx1y2z1, byx1y2z1, bzx1y2z1,
     &                   bxx1y2z2, byx1y2z2, bzx1y2z2,
     &                   bxx2y1z1, byx2y1z1, bzx2y1z1,
     &                   bxx2y1z2, byx2y1z2, bzx2y1z2,
     &                   bxx2y2z1, byx2y2z1, bzx2y2z1,
     &                   bxx2y2z2, byx2y2z2, bzx2y2z2)
                    else if(mmap_type .eq. -4) then
                     call readmap4(xxc,yyc,zzc,r1,r2,z1,z2,s_mgf,
     &                   bmgx,bmgy,bmgz,dxbfbx,dxbfby,dxbfbz,
     &                   dybfbx,dybfby,dybfbz,dzbfbx,dzbfby,dzbfbz,
     &                   brr1z1, bzr1z1, brr1z2, bzr1z2,
     &                   brr2z1, bzr2z1, brr2z2, bzr2z2)
                    endif
                     ss_mgf = sqrt(bmgx**2 + bmgy**2 + bmgz**2)
                     ss_mgf = SIGN(ss_mgf, s_mgf)
                     bmgx1 = bmgx/ss_mgf
                     bmgy1 = bmgy/ss_mgf
                     bmgz1 = bmgz/ss_mgf
                     call trnsuv(bmgx1,bmgy1,bmgz1,bmgx,bmgy,bmgz,itrm)
                   endif

                  call egs_emmov(EIE,chgp,
     &                x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                z(ibkz+no,ipomp+1),
     &                xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                zc(ibkzc+no,ipomp+1),
     &                UNPSTE,VNPSTE,WNPSTE,
     &                u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                w(ibkw+no,ipomp+1),
     &                ss_elf,elcx,elcy,elcz,
     &                ss_mgf,bmgx,bmgy,bmgz,
     &                RM,RMSQ,mark)

                  end if

*-----------------------------------------------------------------------

                  e(ibke+no,ipomp+1) = ec(ibkec+no,ipomp+1)
                  t(ibkt+no,ipomp+1) = tc(ibktc+no,ipomp+1)
                  x(ibkx+no,ipomp+1) = xc(ibkxc+no,ipomp+1)
                  y(ibky+no,ipomp+1) = yc(ibkyc+no,ipomp+1)
                  z(ibkz+no,ipomp+1) = zc(ibkzc+no,ipomp+1)

                  goto 500

*-----------------------------------------------------------------------
*     update coordinate and momentum
*-----------------------------------------------------------------------

  700 continue

                  call gomupr(mark,markp,delt)


*-----------------------------------------------------------------------

  600 continue

                  esav1 = e(ibke+no,ipomp+1)
                  tsav1 = t(ibkt+no,ipomp+1)
                  xsav1 = x(ibkx+no,ipomp+1)
                  ysav1 = y(ibky+no,ipomp+1)
                  zsav1 = z(ibkz+no,ipomp+1)

c' tc = t + timd
c' time limit over : itmak = 1

                  e(ibke+no,ipomp+1) = esav
                  t(ibkt+no,ipomp+1) = tsav
                  x(ibkx+no,ipomp+1) = xsav
                  y(ibky+no,ipomp+1) = ysav
                  z(ibkz+no,ipomp+1) = zsav


                  call timtrs(itmak,mark)


                  e(ibke+no,ipomp+1) = esav1
                  t(ibkt+no,ipomp+1) = tsav1
                  x(ibkx+no,ipomp+1) = xsav1
                  y(ibky+no,ipomp+1) = ysav1
                  z(ibkz+no,ipomp+1) = zsav1

                  call egs_emmov(EIE,chgp,
     &                x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                z(ibkz+no,ipomp+1),
     &                xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                zc(ibkzc+no,ipomp+1),
     &                UNPSTE,VNPSTE,WNPSTE,
     &                u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                w(ibkw+no,ipomp+1),
     &                ss_elf,elcx,elcy,elcz,
     &                ss_mgf,bmgx,bmgy,bmgz,
     &                RM,RMSQ,mark)

                  e(ibke+no,ipomp+1) = esav
                  t(ibkt+no,ipomp+1) = tsav
                  x(ibkx+no,ipomp+1) = xsav
                  y(ibky+no,ipomp+1) = ysav
                  z(ibkz+no,ipomp+1) = zsav


        if(mark.eq.2) then ! reflection case, change the direction again, T.Sato 2020/06/03
         u(ibku+no,ipomp+1)=u(ibku+no,ipomp+1)*uratio
         v(ibkv+no,ipomp+1)=v(ibkv+no,ipomp+1)*vratio
         w(ibkw+no,ipomp+1)=w(ibkw+no,ipomp+1)*wratio
         tmp=sqrt(u(ibku+no,ipomp+1)**2+v(ibkv+no,ipomp+1)**2+
     &                                 w(ibkw+no,ipomp+1)**2)
         u(ibku+no,ipomp+1)=u(ibku+no,ipomp+1)/tmp
         v(ibkv+no,ipomp+1)=v(ibkv+no,ipomp+1)/tmp
         w(ibkw+no,ipomp+1)=w(ibkw+no,ipomp+1)/tmp
        endif

      return
      end

      subroutine uvwnorm(u,v,w,um,vm,wm,ef,bf,bfp)
      implicit none

      double precision u,v,w,um,vm,wm,ef,bf,bfp(3)
      double precision uvwsum,dtnm
      double precision up,vp,wp,dl
      double precision ur,vr,wr,dr
      double precision us,vs,ws,ds
      if(ef.eq.0.0d0.and.bf.ne.0.0d0) then
        dl=(u*bfp(1)+v*bfp(2)+w*bfp(3))/bf
        up=bfp(1)*dl/bf
        vp=bfp(2)*dl/bf
        wp=bfp(3)*dl/bf
        u = u + um
        v = v + vm
        w = w + wm
        dr=(u*bfp(1)+v*bfp(2)+w*bfp(3))/bf
        ur=bfp(1)*dr/bf
        vr=bfp(2)*dr/bf
        wr=bfp(3)*dr/bf
        us=u-ur
        vs=v-vr
        ws=w-wr
        ds=dsqrt(us*us+vs*vs+ws*ws)
        if(ds.eq.0.0d0) then
          u=up
          v=vp
          w=wp
        else
          u=up+dsqrt(1.0d0-dl*dl)*us/ds
          v=vp+dsqrt(1.0d0-dl*dl)*vs/ds
          w=wp+dsqrt(1.0d0-dl*dl)*ws/ds
        endif
      else
        u = u + um
        v = v + vm
        w = w + wm
      endif
      uvwsum = dsqrt(u*u+v*v+w*w)
      u = u/uvwsum
      v = v/uvwsum
      w = w/uvwsum
      return
      end

!-----------------------------------------------------------------------
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12
      subroutine calbfield(unpste,vnpste,wnpste,bfieldp,betadt,cbub)
!
      implicit none

      double precision unpste,vnpste,wnpste,bfieldp(3),betadt,cbub(3)
      double precision c,cb
      c= 2.99792458d0
      cb =c*betadt

      cbub(1)=cb*(vnpste*bfieldp(3)-wnpste*bfieldp(2))
      cbub(2)=cb*(wnpste*bfieldp(1)-unpste*bfieldp(3))
      cbub(3)=cb*(unpste*bfieldp(2)-vnpste*bfieldp(1))
      return
      end
!
      subroutine egs_emstep(TEMSTP,EIE,UNPSTE,VNPSTE,WNPSTE,
     &                      s_elf,elcx,elcy,elcz,
     &                      s_mgf,bmgx,bmgy,bmgz,
     &                      STPEME,STPBME,RM,RMSQ)
!
      implicit none
      double precision TEMSTP,TEMSTPB
      double precision EIE,UNPSTE,VNPSTE,WNPSTE
      double precision s_elf,elcx,elcy,elcz
      double precision s_mgf,bmgx,bmgy,bmgz
      double precision STPEME,STPBME,RM,RMSQ
      double precision EMOVEE,UMOVEE,VMOVEE,WMOVEE
      double precision EFIELDP(3)
      double precision EFLTMP
      double precision BFIELDP(3),CBUB(3),BETADT
      double precision BFLTMP
      double precision PMOVEE,TEMUVW

      TEMSTPB = TEMSTP

      EFIELDP(1)=0.0D0
      EFIELDP(2)=0.0D0
      EFIELDP(3)=0.0D0
      EFLTMP = 0.0d0

c s_elf  : electric field (kV/cm)
c' unit (MV/cm)
      if(s_elf.ne.0.0) then
        EFIELDP(1)=s_elf*elcx*1.0d-3
        EFIELDP(2)=s_elf*elcy*1.0d-3
        EFIELDP(3)=s_elf*elcz*1.0d-3

        EFLTMP = SQRT(EFIELDP(1)*EFIELDP(1)
     &               +EFIELDP(2)*EFIELDP(2)
     &               +EFIELDP(3)*EFIELDP(3))
       end if
!

       BFIELDP(1)=0.0D0
       BFIELDP(2)=0.0D0
       BFIELDP(3)=0.0D0
       BFLTMP  = 0.0d0

c s_mgf  : dipole magnet field  [kG]                             *
c' unit (T)
      if(s_mgf.ne.0.0) then
        BFIELDP(1)=s_mgf*bmgx*1.0d-1
        BFIELDP(2)=s_mgf*bmgy*1.0d-1
        BFIELDP(3)=s_mgf*bmgz*1.0d-1

        BFLTMP = SQRT(BFIELDP(1)*BFIELDP(1)
     &               +BFIELDP(2)*BFIELDP(2)
     &               +BFIELDP(3)*BFIELDP(3))

        BETADT= DSQRT(1.0D0-(RM/EIE)**2)

        call calbfield(UNPSTE,VNPSTE,WNPSTE,BFIELDP,BETADT,
     &                 CBUB)

        PMOVEE = EIE- RMSQ/EIE

        TEMUVW= DSQRT(CBUB(1)*CBUB(1)+CBUB(2)*CBUB(2)
     &               +CBUB(3)*CBUB(3))/PMOVEE
        IF(TEMUVW.GT.0.0D0) THEN
          TEMSTPB=STPBME/TEMUVW
        ENDIF
      endif
!

      IF(EFLTMP.GT.0.0) THEN
        TEMSTP = STPEME*(EIE-RM)/EFLTMP
        TEMSTP = DMAX1(1.0d-4,TEMSTP)
      END IF
      IF(BFLTMP.GT.0.0) THEN
        TEMSTP = DMIN1(TEMSTP,TEMSTPB)
        TEMSTP = DMAX1(1.0d-4,TEMSTP)
      END IF

      return
      end
!
      subroutine egs_emmov(EIE,chgp,x,y,z,xc,yc,zc,
     &                     UNPSTE,VNPSTE,WNPSTE,
     &                     unew,vnew,wnew,
     &                     s_elf,elcx,elcy,elcz,
     &                     s_mgf,bmgx,bmgy,bmgz,
     &                     RM,RMSQ,mark)

      implicit none
      double precision EIE,chgp,x,y,z,xc,yc,zc
      double precision UNPSTE,VNPSTE,WNPSTE
      double precision unew,vnew,wnew
      double precision s_elf,elcx,elcy,elcz
      double precision s_mgf,bmgx,bmgy,bmgz
      double precision RM,RMSQ
      double precision EMOVEE,PMOVEE,SMOVEE,UDCOS
      double precision UMOVEE,VMOVEE,WMOVEE
      double precision EFIELDP(3)
      double precision EFLTMP
      double precision BFIELDP(3),CBUB(3),BETADT
      double precision BFLTMP
      double precision TEMUVW
      double precision TMP
      integer mark ! T.Sato 2023/08/18 avoid infinite loop

      unew = UNPSTE
      vnew = VNPSTE
      wnew = WNPSTE

      if(s_elf.ne.0.0.or.s_mgf.ne.0.0) then

         SMOVEE=sqrt((xc-x)*(xc-x)
     &              +(yc-y)*(yc-y)
     &              +(zc-z)*(zc-z))
         EFIELDP(1)=0.0D0
         EFIELDP(2)=0.0D0
         EFIELDP(3)=0.0D0
         EFLTMP=0.0D0
        if(s_elf.ne.0.0) then
          EFIELDP(1)=s_elf*elcx*1.0d-3
          EFIELDP(2)=s_elf*elcy*1.0d-3
          EFIELDP(3)=s_elf*elcz*1.0d-3
          EFLTMP = SQRT(EFIELDP(1)*EFIELDP(1)
     &                 +EFIELDP(2)*EFIELDP(2)
     &                 +EFIELDP(3)*EFIELDP(3))
        endif
        BFIELDP(1)=0.0D0
        BFIELDP(2)=0.0D0
        BFIELDP(3)=0.0D0
        BFLTMP=0.0D0
        if(s_mgf.ne.0.0) then
          BFIELDP(1)=s_mgf*bmgx*1.0d-1
          BFIELDP(2)=s_mgf*bmgy*1.0d-1
          BFIELDP(3)=s_mgf*bmgz*1.0d-1
          BFLTMP = SQRT(BFIELDP(1)*BFIELDP(1)
     &                 +BFIELDP(2)*BFIELDP(2)
     &                 +BFIELDP(3)*BFIELDP(3))
        endif

        IF(SMOVEE.GT.0.0d0) THEN
          EMOVEE = chgp *( (xc-x)*EFIELDP(1)
     &                   + (yc-y)*EFIELDP(2)
     &                   + (zc-z)*EFIELDP(3))

          IF(EFLTMP.NE.0.0d0.OR.BFLTMP.NE.0.0d0) THEN
            PMOVEE = EIE- RMSQ/EIE
            UDCOS = UNPSTE * EFIELDP(1)
     &            + VNPSTE * EFIELDP(2)
     &            + WNPSTE * EFIELDP(3)
            TMP=1.0D0-(RM/EIE)**2
            if(TMP.lt.0.0) then ! check negative value of dsqrt
             BETADT=0.0
             mark=-1 ! index for lost particle
            else
             BETADT= DSQRT(TMP)
            endif

            call calbfield(UNPSTE,VNPSTE,WNPSTE,BFIELDP,
     &                     BETADT,CBUB)
            UMOVEE = chgp * SMOVEE
     &             * (EFIELDP(1) - UNPSTE*UDCOS+CBUB(1))/PMOVEE
            VMOVEE = chgp * SMOVEE
     &             * (EFIELDP(2) - VNPSTE*UDCOS+CBUB(2))/PMOVEE
            WMOVEE = chgp * SMOVEE
     &             * (EFIELDP(3) - WNPSTE*UDCOS+CBUB(3))/PMOVEE
            call uvwnorm(unew,vnew,wnew,
     &                   umovee,vmovee,wmovee,
     &                   EFLTMP,BFLTMP,BFIELDP)
          endif
        endif
      endif
      return
      end

      subroutine egs_upeng(EIE,chgp,
     &                     x,y,z,xc,yc,zc,
     &                     s_elf,elcx,elcy,elcz,
     &                     RM,nbeta)
      implicit none
      double precision EIE,chgp
      double precision s_elf,elcx,elcy,elcz
      double precision x,y,z,xc,yc,zc
      double precision RM
      double precision EFIELDP(3),EMOVEE
      integer nbeta

      call egs_geteng(EMOVEE,chgp,
     &                x,y,z,xc,yc,zc,
     &                s_elf,elcx,elcy,elcz)

      if(EMOVEE.ne.0.0d0) then
          EIE = EIE + EMOVEE
          IF ((EIE.LT.RM)) THEN
              EIE = RM
              nbeta = 3
          endif
      endif
      return
      end
!
      subroutine egs_geteng(EMOVEE,chgp,
     &                      x,y,z,xc,yc,zc,
     &                      s_elf,elcx,elcy,elcz)

      implicit none
      double precision EMOVEE,chgp
      double precision s_elf,elcx,elcy,elcz
      double precision x,y,z,xc,yc,zc
      double precision EFIELDP(3)
      integer nbeta

      EMOVEE = 0.0d0
      if(s_elf.ne.0.0) then
         EFIELDP(1)=s_elf*elcx*1.0d-3
         EFIELDP(2)=s_elf*elcy*1.0d-3
         EFIELDP(3)=s_elf*elcz*1.0d-3

         EMOVEE = chgp *( (xc-x)*EFIELDP(1)
     &                  + (yc-y)*EFIELDP(2)
     &                  + (zc-z)*EFIELDP(3))

      endif

      return
      end

      subroutine  egs_getbmg(x,y,z,
     &                      a_mag,b_mag,s_mag,p_mag,t_mag,u_mag,
     %                      s_mgf,bmgx,bmgy,bmgz)

      implicit none

*       x,y,z                                                          *
*                : initial position                                    *
*       a_mag  : magnet gap(cm)                                        *
*       b_mag  : magnet field at pole tip  [kG]                        *
*       s_mag  : speicies of magnet dypole:2, quad:4, sext:6, oct:8    *
*       p_mag  : phase of magnet                                       *
*       t_mag  : transform id                                          *

      double precision x,y,z
      double precision a_mag,b_mag,s_mag,p_mag,t_mag,u_mag
      double precision s_mgf,bmgx,bmgy,bmgz
      double precision xx0,yy0,zz0
      double precision s0,bx0,by0,bz0
      double precision xisp,yisp
      complex(kind(0d0))  cb,cbisp

      integer itrs,isp,isp2

      itrs = nint( t_mag )


      call trnsxx(x,y,z,xx0,yy0,zz0,itrs)

      isp = nint( s_mag )
      s_mgf = 0.0

      if(isp.eq.2) then
        s0 = b_mag
        bx0 = 0.0
        by0 = 1.0
        bz0 = 0.0
        s_mgf = b_mag
      elseif(isp.eq.4.or.isp.eq.6.or.isp.eq.8) then
        isp2 = isp/2
        cb = cmplx(xx0,yy0)
        cbisp = cb**(isp2-1)
        yisp = real(cbisp)
        xisp = aimag(cbisp)

        bx0 = b_mag / a_mag * xisp
        by0 = b_mag / a_mag * yisp
        bz0 = 0.0

        s_mgf =  sqrt( bx0*bx0 + by0*by0 + bz0*bz0)
        if(s_mgf.gt.0.0) then
          bx0 = bx0 / s_mgf
          by0 = by0 / s_mgf
          bz0 = bz0 / s_mgf
        else
          bx0 = 0.0
          by0 = 1.0
          bz0 = 0.0
        endif

      endif

      call trnsuv(bx0,by0,bz0,
     &            bmgx,bmgy,bmgz,itrs)

      return
      end

      subroutine clrcmn3()
      real*8
     * thard,tmscat,tinel,hardstep,
     * sig,scpow,dedx,sig0,ams

      common/egs5cmn3/thard,tinel,tmscat,hardstep,sig,scpow,dedx,sig0,
     $                ams
!$OMP THREADPRIVATE(/egs5cmn3/)
        thard = 0.0d0
        tinel = 0.d0
        tmscat = 0.0d0

      return
      end
