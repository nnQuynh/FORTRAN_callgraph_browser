************************************************************************
*                                                                      *
      subroutine ggmsor(x,y,z,u,v,w,nmed,iblz,mark,markp,ici)
*                                                                      *
*       find cell of GG at initial source                              *
*       modified by K.Niita on 2002/05/10                              *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

      include 'err.inc'

*-----------------------------------------------------------------------

      common /inout/  in,io
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)

      common /regdm/  idmg(kvlmax)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /paraj/  mstz(300), parz(300)

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      integer iii0,kkk0
      common /itettal2/ iii0,kkk0
!$OMP THREADPRIVATE(/itettal2/)
*-----------------------------------------------------------------------

               ierr  = 0
               isrr  = 0

         if( mark .eq. 0 .or. mark .eq. 2 .or. markp .eq. 1 ) goto 100

*-----------------------------------------------------------------------

  500    continue

            if( ici .le. 0 ) then

               call initsk

               ih  = 0
               iii = 0 !FURUTA20160607

            else

               ih = idgr( iblz )

            end if

*-----------------------------------------------------------------------

               xxx = x
               yyy = y
               zzz = z

               uuu = u
               vvv = v
               www = w

               jsu  = 0
               levp = 0

               icl  = 0
               lev  = 0

*-----------------------------------------------------------------------
!$OMP CRITICAL (setsuf_crit)

            call setsuf(ih,io,ierr)

!$OMP END CRITICAL (setsuf_crit)
*-----------------------------------------------------------------------

               if( ierr .ne. 0 ) then

                  if( ici .eq. 0 ) then

                     if( isrr .eq. 0 ) then

                        x = x + parz(28) * u
                        y = y + parz(28) * v
                        z = z + parz(28) * w

                        isrr = isrr + 1
                        goto 500

                     else




        ErrCha = ''
        ErrID = 'L:106/R:ggmsor/F:ggs02.f' !E03_003_001
        call ErrWrite(ErrID,ErrCha)

        write(*,'(''Source is not in any cell!'')')
        write(*,'(''batch & history number:'',2i20)') nobch,nocas
        write(*,'(''location (x,y,z)      :'',3es16.8)') x,y,z

                        call parastop( 601 )

                     end if

                  else

                     if( ierr .eq. 1 ) then

cKN 2013/12/25 no defined
                        mark = -2
                        return

                     else if( ierr .eq. 2 ) then

cKN 2013/12/25 double defined
                        mark = -3
                        return

                     end if

                  end if

               end if

*-----------------------------------------------------------------------

            iblz  = idrg(icl)

*-----------------------------------------------------------------------

            if( ici .le. 0 ) then

                  nmed  = idmg(icl)

               if( nmed .gt. 0 ) then

                  nmed = idnm( nmed )

               end if

            else if( ici .gt. 0 .and. icl .ne. ih ) then

                  iblz  = idrg(ih)
                  iblz1 = idrg(ih)
                  iblz2 = idrg(icl)

                  iii0 = iii
                  kkk0 = kkk

                  nmed  = idmg(icl)

                  if( nmed .gt. 0 ) then

                     nmed = idnm( nmed )

                  end if

                  mark = -3
                  return

            end if

*-----------------------------------------------------------------------
*        store initial region ( universe and lattice )
*-----------------------------------------------------------------------

  100    continue

                     iblz1 = iblz
                     iblz2 = iblz

                     iii0 = iii
                     kkk0 = kkk

                     ilev1 = lev
                     ilev2 = lev

*-----------------------------------------------------------------------
*        with level structure
*-----------------------------------------------------------------------

            if( lev .gt. 0 ) then

               if( lat(1,icl) .eq. 0 .or.
     &            abs(lat(1,icl)) .eq. 3) then !FURUTA20150714 TETRA

                     ilev1 = lev

                  do kk = lev - 1, 0, -1

                     k = lev - kk

                     ilat1(1,k) = idrg(nint(udt(7,kk)))
                     ilat1(2,k) = lat(1,nint(udt(7,kk)))
                     ilat1(3,k) = nint(udt( 8,kk))
                     ilat1(4,k) = nint(udt( 9,kk))
                     ilat1(5,k) = nint(udt(10,kk))

                  end do

               else

                     ilev1 = lev + 1

                     ilat1(1,1) = iblz
                     ilat1(2,1) = lat(1,icl)
                     ilat1(3,1) = iii
                     ilat1(4,1) = jjj
                     ilat1(5,1) = kkk

                  do kk = lev - 1, 0, -1

                     k = lev - kk + 1

                     ilat1(1,k) = idrg(nint(udt(7,kk)))
                     ilat1(2,k) = lat(1,nint(udt(7,kk)))
                     ilat1(3,k) = nint(udt( 8,kk))
                     ilat1(4,k) = nint(udt( 9,kk))
                     ilat1(5,k) = nint(udt(10,kk))

                  end do

               end if

                     ilev2 = ilev1

               do k = 1, lev

                  do j = 1, 5

                     ilat2(j,k) = ilat1(j,k)

                  end do

               end do

            end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine ggmpgp(icr,x,y,z,xc,yc,zc,u,v,w,
     &                  nmed,iblz,mark,markp)
*                                                                      *
*       GG geometry                                                    *
*                                                                      *
*       mark:     description                                          *
*         0 : pass the forward surface                                 *
*         1 : collide between present position and the forward surface *
*         2 : reach the reflection surface                             *
*        -1 : outgoing to the void region                              *
*        -2 : error: lost particle                                     *
*        -3 : error: inconsitent after crossing                        *
*        -4 : error: region is differnt after collision                *
*        -5 : error: region is the same after crossing                 *
*                                                                      *
*       modified by K.Niita on 2001/11/29                              *
*       modified by NAIS on 2019/11/19                                 *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      common /regdm/  idmg(kvlmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

      common /gsline/ nowgshow, igsline

      data dps1 / 0.12345678d-07 /

*-----------------------------------------------------------------------

            x0 = x
            y0 = y
            z0 = z

            xc0 = xc
            yc0 = yc
            zc0 = zc

*-----------------------------------------------------------------------

  100 continue

               call ggmprp(icr,x,y,z,xc,yc,zc,u,v,w,
     &                     nmed,iblz,mark,markp,icl0)

      if( nowgshow .ne. 0 ) then

         if( mark .eq. 0 .and. markp .eq. 0 ) then

               icmlnow = 0
               icmlnxt = 0

            if( nowgshow .eq. 2 ) then

               icmlnow = icl0
               icmlnxt = icl

            else if( nowgshow .eq. 1 ) then

               icmlnow = idmg(icl0)
               icmlnxt = idmg(icl)

            end if

            if( icmlnow .eq. icmlnxt ) then

               x = xc
               y = yc
               z = zc

               xc = xc0
               yc = yc0
               zc = zc0

               goto 100

            end if

         end if

  200 continue

               x = x0
               y = y0
               z = z0

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine ggmprp(icr,x,y,z,xc,yc,zc,u,v,w,
     &                  nmed,iblz,mark,markp,icl0)
*                                                                      *
*       GG geometry                                                    *
*                                                                      *
*       mark:     description                                          *
*         0 : pass the forward surface                                 *
*         1 : collide between present position and the forward surface *
*         2 : reach the reflection surface                             *
*        -1 : outgoing to the void region                              *
*        -2 : error: lost particle                                     *
*        -3 : error: inconsitent after crossing                        *
*        -4 : error: region is differnt after collision                *
*        -5 : error: region is the same after crossing                 *
*                                                                      *
*            modified by K.Niita on 2001/11/29                         *
*       Last modified by K.Niita on 2020/04/09                         *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /rcomon/ rcasc

*-----------------------------------------------------------------------

      common/inout/in,io

      common /regdm/  idmg(kvlmax)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

      common /celdb/  idsn(kvlmax), idtn(kvlmax)

*-----------------------------------------------------------------------

      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /regcrs/ icrsflx
      common /tlcost/ costha, uang(3), nsurf
!$OMP THREADPRIVATE(/tlcost/)

*-----------------------------------------------------------------------

      dimension oudt(10,0:mxlv)

*-----------------------------------------------------------------------
cFURUTA20150714 TETRA
      integer iflag
*-----------------------------------------------------------------------

            ierr = 0

            dpr = sqrt( (xc-x)**2 + (yc-y)**2 + (zc-z)**2 )

*-----------------------------------------------------------------------
*        check cell : mark .ne. 0, 2  and markp = 0
*-----------------------------------------------------------------------

            if( mark .ne. 0 .and. mark .ne. 2 .and. markp .eq. 0 ) then

                     ici = 1

                     call ggmsor(x,y,z,u,v,w,nmed,iblz,mark,markp,ici)

                  if( mark .le. -2 ) return

            end if

*-----------------------------------------------------------------------
*           save initial state
*-----------------------------------------------------------------------

                  iii0 = iii
                  jjj0 = jjj
                  kkk0 = kkk
                  icl0 = icl

                  lev0  = lev
                  levp0 = levp

               do k = 0, lev
               do j = 1, 10

                  oudt(j,k) = udt(j,k)

               end do
               end do

*-----------------------------------------------------------------------
*        check minimum distance
*-----------------------------------------------------------------------

                  kdb = 0

            if( lca(icl) .lt. 0 ) call chkcll(icl,3,j)

            call fsurf(icl)

               if( kdb .ne. 0 ) then

                  irror = 1
                  goto 999

               end if

               if(mark.eq.-2)return !FURUTA20150714 TETRA

*-----------------------------------------------------------------------
*        collision occured
*-----------------------------------------------------------------------

      if( dpr .lt. dls ) then

               mark  = 1
               markp = 1

               xxx = xxx + uuu * dpr
               yyy = yyy + vvv * dpr
               zzz = zzz + www * dpr

            do l = 0, lev - 1

               udt(1,l)=udt(1,l)+dpr*udt(4,l)
               udt(2,l)=udt(2,l)+dpr*udt(5,l)
               udt(3,l)=udt(3,l)+dpr*udt(6,l)

            end do

*-----------------------------------------------------------------------
*        crossing boundary
*-----------------------------------------------------------------------

      else

               xc = x + u * dls
               yc = y + v * dls
               zc = z + w * dls

               xxx = xxx + uuu * dls
               yyy = yyy + vvv * dls
               zzz = zzz + www * dls

            do l = 0, lev - 1

               udt(1,l)=udt(1,l)+dls*udt(4,l)
               udt(2,l)=udt(2,l)+dls*udt(5,l)
               udt(3,l)=udt(3,l)+dls*udt(6,l)

            end do

*-----------------------------------------------------------------------

               cs  = icrsflx
               jsu = jap

cFURUTA20150714 TETRA---------------------------------------------------

         if(kkk.le.10000)then
          if(ksu(jsu).ge.0) then
           iflag=0
          else
           iflag=1
          endif
         else
          iflag=0
         endif

*-----------------------------------------------------------------------
*        normal surface
*-----------------------------------------------------------------------

         if( iflag .eq. 0 ) then

               call fcell(io,ierr,cs)

                  mark  = 0
                  icl   = iap
                  if(junf.ne.0)then
                   if(kkk.gt.10000.and.lat(1,icl).eq.0) kkk=10000
                   if(kkk.gt.10000.and.abs(lat(1,icl)).eq.3)iii=jjj
                  endif

               if( ierr .ne. 0 ) then

                  if( ierr .eq. 1 ) irror = 2
                  if( ierr .eq. 2 ) irror = 22

                  goto 999

               end if

               if( kdb .ne. 0 ) then

                  irror = 3
                  goto 999

               end if

*-----------------------------------------------------------------------
*        reflect surface
*-----------------------------------------------------------------------

         else

               call refsuf(cs)

                  mark  = 2

               if( lev .eq. 0 ) then

                  u = uuu
                  v = vvv
                  w = www

               else

                  u = udt(4,0)
                  v = udt(5,0)
                  w = udt(6,0)

               end if

         end if

*-----------------------------------------------------------------------
*        check new cell and store the information
*-----------------------------------------------------------------------

               iblz2 = idrg(icl)
               nmed  = idmg(icl)

            if( nmed .gt. 0 ) then

               nmed = idnm( nmed )

            end if

               iblz = iblz2

               if( nmed .eq. -1 ) mark = -1

            if( mark .eq. 0 .and.
     &          icl  .eq. icl0  .and.
     &          levp .eq. levp0 .and.
     &          lev  .eq. lev0  .and.
     &          iii  .eq. iii0  .and.
     &          jjj  .eq. jjj0  .and.
     &          kkk  .eq. kkk0 ) then

                  if( lev .gt. 0 ) then

                     do k = 0, lev

                        if( nint(udt(7,k)) .ne. nint(oudt(7,k)) )
     &                  goto 5000

                     end do

                  end if

               mark = -5
               return

            end if

 5000       continue

*-----------------------------------------------------------------------

                     ilev2 = lev

            if( lev .gt. 0 ) then

               if( lat(1,icl) .eq. 0 .or.
     &            abs(lat(1,icl)) .eq. 3) then !FURUTA20150714 TETRA

                     ilev2 = lev

                  do kk = lev - 1, 0, -1

                     k = lev - kk

                     ilat2(1,k) = idrg(nint(udt(7,kk)))
                     ilat2(2,k) = lat(1,nint(udt(7,kk)))
                     ilat2(3,k) = nint(udt( 8,kk))
                     ilat2(4,k) = nint(udt( 9,kk))
                     ilat2(5,k) = nint(udt(10,kk))

                  end do

               else

                     ilev2 = lev + 1

                     ilat2(1,1) = iblz2
                     ilat2(2,1) = lat(1,icl)
                     ilat2(3,1) = iii
                     ilat2(4,1) = jjj
                     ilat2(5,1) = kkk

                  do kk = lev - 1, 0, -1

                     k = lev - kk + 1

                     ilat2(1,k) = idrg(nint(udt(7,kk)))
                     ilat2(2,k) = lat(1,nint(udt(7,kk)))
                     ilat2(3,k) = nint(udt( 8,kk))
                     ilat2(4,k) = nint(udt( 9,kk))
                     ilat2(5,k) = nint(udt(10,kk))

                  end do

               end if

            end if

*-----------------------------------------------------------------------

               costha = abs( cs )

               uang(1) = ang(1)
               uang(2) = ang(2)
               uang(3) = ang(3)

               nsurf = iabs( jsu )

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*     error in GG : lost particles : mark = -2
*-----------------------------------------------------------------------

      return

  999    continue

      if( icr .ne. 0 ) then

         write(io,'(/80(''-''))')

         write(io,'(a,f16.0,i10)')
     &            '      ncas  nocas    =', rcasc, nocas

         if( nmed .gt. 0 ) then

            write(io,'( ''    mat ='',i5)') idmn(nmed)

         else

            write(io,'( ''    mat = void'')')

         end if

         write(io,'( ''   iblz1 ='',i7)') iblz1
         write(io,'( ''   iblz2 ='',i7)') iblz2
         write(io,'( ''    mark ='',i5)') mark

         write(io,'(/''    x, y, z  ='',1p3e17.9)') x,y,z
         write(io,'( ''    xc,yc,zc ='',1p3e17.9)') xc,yc,zc
         write(io,'( ''    u, v, w  ='',1p3e17.9)') u,v,w

      if( irror .eq. 1 ) then
         write(io,'(/'' ERROR : cannot find next surface in track'')')
      else if( irror .eq. 2 ) then
         write(io,'(/'' ERROR : crossing surface does not exist in'',
     &               '' this cell ='',i7)') icl
         write(io,'( '' surface ='',i7)') idsn(jsu)
      else if( irror .eq. 22 ) then
         write(io,'(/'' ERROR : zero lattice element hit in fcell.'')')
      else if( irror .eq. 3 ) then
         write(io,'(/'' ERROR : crossing surface does not exist in'',
     &               '' all cell'')')
         write(io,'( '' surface ='',i7)') idsn(jsu)
      end if

      end if

*-----------------------------------------------------------------------

         mark = -2
         return

*-----------------------------------------------------------------------

      end

************************************************************************
*                                                                      *
      subroutine ggmdis(dist,x,y,z,u,v,w,nmed,iblz,mark,markp)
*                                                                      *
*       GG geometry                                                    *
*       detect the distance to the boundary                            *
*                                                                      *
*       mark:     description                                          *
*        -2 : error: lost particle                                     *
*        -6 : error: point is just on the boundary                     *
*                                                                      *
*       modified by K.Niita on 2002/11/21                              *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

      common /paraj/  mstz(300), parz(300)
*-----------------------------------------------------------------------

            icn = 0
            isrr = 0

*-----------------------------------------------------------------------
*        check cell : mark .ne. 0, 2  and markp = 0
*-----------------------------------------------------------------------

            if( mark .ne. 0 .and. mark .ne. 2 .and. markp .eq. 0 ) then

                     ici = 1

                     call ggmsor(x,y,z,u,v,w,nmed,iblz,mark,markp,ici)

                  if( mark .le. -2 ) return

            end if

*-----------------------------------------------------------------------
*        check minimum distance
*-----------------------------------------------------------------------

  100       continue

                  kdb = 0

            if( lca(icl) .lt. 0 ) call chkcll(icl,3,j)

            call fsurf(icl)

               if( kdb .ne. 0 ) then

                  mark = -6
                  return

               end if

               dist = dls

               markp = 1

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine ggmnew(icr,dpr,x,y,z,xc,yc,zc,u,v,w,ec,
     &                  nmed,iblz,mark,markp)
*                                                                      *
*       GG geometry                                                    *
*                                                                      *
*       mark:     description                                          *
*         0 : pass the forward surface                                 *
*         1 : collide between present position and the forward surface *
*         2 : reach the reflection surface                             *
*        -1 : outgoing to the void region                              *
*        -2 : error: lost particle                                     *
*        -3 : error: inconsitent after crossing                        *
*        -4 : error: region is differnt after collision                *
*        -5 : error: region is the same after crossing                 *
*                                                                      *
*       modified by K.Niita on 2002/05/02                              *
*                                                                      *
************************************************************************
      use moddas_region
      use moddas_ggs

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /rcomon/ rcasc

*-----------------------------------------------------------------------

      common/inout/in,io

      common /regdm/  idmg(kvlmax)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

      common /celdb/  idsn(kvlmax), idtn(kvlmax)

*-----------------------------------------------------------------------

      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /regcrs/ icrsflx
      common /tlcost/ costha, uang(3), nsurf
!$OMP THREADPRIVATE(/tlcost/)

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /smireg/ nsreg, isgrc, isgrt, ksmir, ismir

*-----------------------------------------------------------------------

      dimension oudt(10,0:mxlv)

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------
cFURUTA20150714 TETRA
      integer iflag
*-----------------------------------------------------------------------
*        error : collision occured
*-----------------------------------------------------------------------

         if( dpr .lt. dls ) then

            irror = 4
            goto 999

         end if

            ierr = 0

*-----------------------------------------------------------------------
*           save initial state
*-----------------------------------------------------------------------

                  xxx0 = xxx
                  yyy0 = yyy
                  zzz0 = zzz

                  uuu0 = uuu
                  vvv0 = vvv
                  www0 = www

                  iii0 = iii
                  jjj0 = jjj
                  kkk0 = kkk

                  icl0  = icl
                  lev0  = lev
                  levp0 = levp

               do k = 0, lev
               do j = 1, 10

                  oudt(j,k) = udt(j,k)

               end do
               end do

*-----------------------------------------------------------------------
*        crossing boundary
*-----------------------------------------------------------------------

               xc = x + u * dls
               yc = y + v * dls
               zc = z + w * dls

*-----------------------------------------------------------------------

               cs  = icrsflx
               jsu = jap

               xxx = xxx + uuu * dls
               yyy = yyy + vvv * dls
               zzz = zzz + www * dls

            do l = 0, lev - 1

               udt(1,l)=udt(1,l)+dls*udt(4,l)
               udt(2,l)=udt(2,l)+dls*udt(5,l)
               udt(3,l)=udt(3,l)+dls*udt(6,l)

            end do

cFURUTA20150714 TETRA----------------------------------------------------
         if(kkk.le.10000)then
          if(ksu(jsu).ge.0) then
           iflag=0
          else
           iflag=1
          endif
         else
          iflag=0
         endif

*-----------------------------------------------------------------------
*     normal surface
*-----------------------------------------------------------------------

      if( iflag .eq. 0 ) then

               call fcell(io,ierr,cs)

                  mark  = 0
                  icl   = iap
                  if(junf.ne.0)then
                   if(kkk.gt.10000.and.lat(1,icl).eq.0)kkk=10000
                   if(kkk.gt.10000.and.abs(lat(1,icl)).eq.3)iii=jjj !FURUTA20240110
                 endif

                  snt = abs( cs )

*-----------------------------------------------------------------------
*              error check
*-----------------------------------------------------------------------

               if( ierr .ne. 0 ) then

                  if( ierr .eq. 1 ) irror = 2
                  if( ierr .eq. 2 ) irror = 22

                  goto 999

               end if

               if( kdb .ne. 0 ) then

                  irror = 3
                  goto 999

               end if

               if( icl  .eq. icl0  .and.
     &             levp .eq. levp0 .and.
     &             lev  .eq. lev0  .and.
     &             iii  .eq. iii0  .and.
     &             jjj  .eq. jjj0  .and.
     &             kkk  .eq. kkk0 ) then

                     if( lev .gt. 0 ) then

                        do k = 0, lev

                           if( nint(udt(7,k)) .ne. nint(oudt(7,k)) )
     &                     goto 5000

                        end do

                     end if

                  mark = -5
                  return

               end if

 5000          continue

*-----------------------------------------------------------------------
*        Super Mirror for neutron,
*        e > 10 eV and  sin > 0.001 sharp cut off
*-----------------------------------------------------------------------

         if( nsreg .gt. 0 .and. ityp .eq. 2 .and.
     &     ( ec .le. 10.e-6 .or. snt .gt. 1.e-3 ) ) then

*-----------------------------------------------------------------------

                     iblz2 = idrg(icl)
                     ilev2 = lev

            if( lev .gt. 0 ) then

               if( lat(1,icl) .eq. 0 .or.
     &            abs(lat(1,icl)) .eq. 3) then !FURUTA20150714 TETRA

                     ilev2 = lev

                  do kk = lev - 1, 0, -1

                     k = lev - kk

                     ilat2(1,k) = idrg(nint(udt(7,kk)))
                     ilat2(2,k) = lat(1,nint(udt(7,kk)))
                     ilat2(3,k) = nint(udt( 8,kk))
                     ilat2(4,k) = nint(udt( 9,kk))
                     ilat2(5,k) = nint(udt(10,kk))

                  end do

               else

                     ilev2 = lev + 1

                     ilat2(1,1) = iblz2
                     ilat2(2,1) = lat(1,icl)
                     ilat2(3,1) = iii
                     ilat2(4,1) = jjj
                     ilat2(5,1) = kkk

                  do kk = lev - 1, 0, -1

                     k = lev - kk + 1

                     ilat2(1,k) = idrg(nint(udt(7,kk)))
                     ilat2(2,k) = lat(1,nint(udt(7,kk)))
                     ilat2(3,k) = nint(udt( 8,kk))
                     ilat2(4,k) = nint(udt( 9,kk))
                     ilat2(5,k) = nint(udt(10,kk))

                  end do

               end if

            end if

*-----------------------------------------------------------------------
*           check region and Q value
*-----------------------------------------------------------------------

                  qin = 4.d0 * 3.141592653589793d0 * snt
     &                * sqrt( ec / 8.180425d-8 )

                     idsm = isgrt
                     jdsm = -1

                     kdsm = ksmir
                     ldsm = 0

            do ir = 1, nsreg

                     ldsm  = ldsm + 1
                     smm   = das_ksmir(kdsm+ldsm)
                     ldsm  = ldsm + 1
                     sr0   = das_ksmir(kdsm+ldsm)
                     ldsm  = ldsm + 1
                     sqc   = das_ksmir(kdsm+ldsm)
                     ldsm  = ldsm + 1
                     sam   = das_ksmir(kdsm+ldsm)
                     ldsm  = ldsm + 1
                     swm   = das_ksmir(kdsm+ldsm)

                     ic0 = 0
                     ic1 = 0
                     rpr = sr0

                  if( qin .gt. sqc ) then

                     rpr = 0.5d0 * sr0
     &                   * ( 1.d0 - tanh( ( qin - smm * sqc ) / swm ) )
     &                   * ( 1.d0 - sam * ( qin - sqc ) )

                  end if

               if( unirn(dummy) .lt. rpr ) ic0 = 1

                     jdsm = jdsm + 1
                     ntrn = idas_isgrt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_isgrt(idsm+jdsm)

                     jj = 0

               if( ic0 .ne. 0 ) then

                  do k = 1, ntrn

                     call tregck(iblz1,ilev1,ilat1,
     &                           mtrn,idas_isgrt(idsm+jdsm+1),jj,ic1)

                  end do

               end if

                     jdsm = jdsm + mtrn

                     jdsm = jdsm + 1
                     ntrn = idas_isgrt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_isgrt(idsm+jdsm)

               if( ic1 .ne. 0 ) then

                     jj  = 0

                  do k = 1, ntrn

                     call tregck(iblz2,ilev2,ilat2,
     &                           mtrn,idas_isgrt(idsm+jdsm+1),jj,ic2)

                     if( ic2 .ne. 0 ) goto 400

                  end do

               end if

                     jdsm = jdsm + mtrn

            end do

               goto 6000

*-----------------------------------------------------------------------
*           super mirror, restore initial state
*-----------------------------------------------------------------------

  400       continue


                  xxx = xxx0
                  yyy = yyy0
                  zzz = zzz0

                  uuu = uuu0
                  vvv = vvv0
                  www = www0

                  iii = iii0
                  jjj = jjj0
                  kkk = kkk0

                  icl  = icl0
                  lev  = lev0
                  levp = levp0

               do k = 0, lev
               do j = 1, 10

                  udt(j,k) = oudt(j,k)

               end do
               end do

               cs  = icrsflx
               jsu = jap

               xxx = xxx + uuu * dls
               yyy = yyy + vvv * dls
               zzz = zzz + www * dls

            do l = 0, lev - 1

               udt(1,l)=udt(1,l)+dls*udt(4,l)
               udt(2,l)=udt(2,l)+dls*udt(5,l)
               udt(3,l)=udt(3,l)+dls*udt(6,l)

            end do

*-----------------------------------------------------------------------
*           super mirror, reflection
*-----------------------------------------------------------------------

               call refsuf(cs)

                  mark  = 2

               if( lev .eq. 0 ) then

                  u = uuu
                  v = vvv
                  w = www

               else

                  u = udt(4,0)
                  v = udt(5,0)
                  w = udt(6,0)

               end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*     reflect surface
*-----------------------------------------------------------------------

      else

               call refsuf(cs)

                  mark  = 2

               if( lev .eq. 0 ) then

                  u = uuu
                  v = vvv
                  w = www

               else

                  u = udt(4,0)
                  v = udt(5,0)
                  w = udt(6,0)

               end if

      end if

*-----------------------------------------------------------------------
*        store the lattice information
*-----------------------------------------------------------------------

                     ilev2 = lev

            if( lev .gt. 0 ) then

               if( lat(1,icl) .eq. 0 .or.
     &            abs(lat(1,icl)) .eq. 3) then !FURUTA20150714 TETRA

                     ilev2 = lev

                  do kk = lev - 1, 0, -1

                     k = lev - kk

                     ilat2(1,k) = idrg(nint(udt(7,kk)))
                     ilat2(2,k) = lat(1,nint(udt(7,kk)))
                     ilat2(3,k) = nint(udt( 8,kk))
                     ilat2(4,k) = nint(udt( 9,kk))
                     ilat2(5,k) = nint(udt(10,kk))

                  end do

               else

                     ilev2 = lev + 1

                     ilat2(1,1) = iblz2
                     ilat2(2,1) = lat(1,icl)
                     ilat2(3,1) = iii
                     ilat2(4,1) = jjj
                     ilat2(5,1) = kkk

                  do kk = lev - 1, 0, -1

                     k = lev - kk + 1

                     ilat2(1,k) = idrg(nint(udt(7,kk)))
                     ilat2(2,k) = lat(1,nint(udt(7,kk)))
                     ilat2(3,k) = nint(udt( 8,kk))
                     ilat2(4,k) = nint(udt( 9,kk))
                     ilat2(5,k) = nint(udt(10,kk))

                  end do

               end if

            end if

*-----------------------------------------------------------------------
*        new cell information
*-----------------------------------------------------------------------

 6000    continue

               iblz2 = idrg(icl)
               nmed  = idmg(icl)

            if( nmed .gt. 0 ) then

               nmed = idnm( nmed )

            end if

               iblz = iblz2

               if( nmed .eq. -1 ) mark = -1

*-----------------------------------------------------------------------

               costha = abs( cs )

               uang(1) = ang(1)
               uang(2) = ang(2)
               uang(3) = ang(3)

               nsurf = iabs( jsu )

*-----------------------------------------------------------------------
*     error in GG : lost particles : mark = -2
*-----------------------------------------------------------------------

      return

  999    continue

      if( icr .ne. 0 ) then

         write(io,'(/80(''-''))')

         write(io,'(a,f16.0,i10)')
     &            '      ncas  nocas    =', rcasc, nocas

         if( nmed .gt. 0 ) then

            write(io,'( ''    mat ='',i5)') idmn(nmed)

         else

            write(io,'( ''    mat = void'')')

         end if

         write(io,'( ''   iblz1 ='',i7)') iblz1
         write(io,'( ''   iblz2 ='',i7)') iblz2
         write(io,'( ''    mark ='',i5)') mark

         write(io,'(/''    x, y, z  ='',1p3e17.9)') x,y,z
         write(io,'( ''    xc,yc,zc ='',1p3e17.9)') xc,yc,zc
         write(io,'( ''    u, v, w  ='',1p3e17.9)') u,v,w

      if( irror .eq. 0 ) then
         write(io,'(/'' ERROR : cannot find cell after collision'')')
      else if( irror .eq. 4 ) then
         write(io,'(/'' ERROR : delt is smaller than dls'')')
      else if( irror .eq. 1 ) then
         write(io,'(/'' ERROR : cannot find next surface in track'')')
      else if( irror .eq. 2 ) then
         write(io,'(/'' ERROR : crossing surface does not exist in'',
     &               '' this cell ='',i7)') icl
         write(io,'( '' surface ='',i7)') idsn(jsu)
      else if( irror .eq. 22 ) then
         write(io,'(/'' ERROR : zero lattice element hit in fcell.'')')
      else if( irror .eq. 3 ) then
         write(io,'(/'' ERROR : crossing surface does not exist in'',
     &               '' all cell'')')
         write(io,'( '' surface ='',i7)') idsn(jsu)
      end if

      end if

*-----------------------------------------------------------------------

         mark = -2
         return

*-----------------------------------------------------------------------

      end

