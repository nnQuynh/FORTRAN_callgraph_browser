************************************************************************
*                                                                      *
      subroutine cggsor(x,y,z,u,v,w,nmed,iblz,mark,markp,ici)
*                                                                      *
*       find region for CG                                             *
*       modified by K.Niita on 2001/11/29                              *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      real*4 gms(1)
      equivalence (das,gms)

*-----------------------------------------------------------------------
c
      common/inout/in,io
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)

      common/ark/ngy,nll,nbb,nzy(3),xd(3)
      common/arar/nby,nlev,nar,nq,iaw,iay,nf,nx1(3)
!$OMP THREADPRIVATE(/arar/)
      common/gomloc/kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     &              kkr2,knsr,kvoll,nadd,ldata,ltma,lfpd,numr,irtru,
     &              numb,nir,kbiz,kbcz
!$OMP THREADPRIVATE(/gomloc/)
      common/parem/xb(3),wb(3),wp(3),xp(3),rin,rout,pinf,dist,jr,idbg,
     &             irprim,nasc,lsurf,nbo,lri,lro,kloop,loop,itype,noa
!$OMP THREADPRIVATE(/parem/)

      common/cggstr/ lgms

      integer blzold
      common/orgi/dist0,markg,nmedg,nblz,blzold,irpold
!$OMP THREADPRIVATE(/orgi/)
      common/mgomv/mus,muz,ll,ipret,iflow,iect,nlo,igx
!$OMP THREADPRIVATE(/mgomv/)
      common/tape/intt,iot,iout,iou2,idm(4),ioecg
!$OMP THREADPRIVATE(/tape/)
      common/ss/s1
!$OMP THREADPRIVATE(/ss/)
      common/repeat/ip(20)

      common /paraj/  mstz(300), parz(300)

      common/flag/ierror
!$OMP THREADPRIVATE(/flag/)

      common/tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)

      dimension xa(3),dc(3)

*-----------------------------------------------------------------------

               ierr  = 0
               isrr  = 0

         if( mark .eq. 0 .or. mark .eq. 2 .or. markp .eq. 1 ) goto 100

*-----------------------------------------------------------------------

            if( ici .gt. 0 ) then

               ih = iblz

            end if

*-----------------------------------------------------------------------

  500    continue

                  dist = 0.0

                  xa(1) = x
                  xa(2) = y
                  xa(3) = z

                  dc(1) = u
                  dc(2) = v
                  dc(3) = w

                  ll  = 0
                  nl  = 0
                  nlu = 0

                  irprim = 0
                  i001 = ip(12)
                  i002 = i001 + 6 * nlev + 4 - 1
                  if(nlev.gt.0) call resetz(gms(lgms),i001,i002)

                  ip4 = ip(4)
                  ip5 = ip(5)
                  ip7 = ip(7)
                  ip8 = ip(8)
                  ip0 = ip(10)

                  ierror = 0

            call cali(nl,nlu,xa,dc,ip,
     &                gms(lgms-1+ip4),gms(lgms-1+ip5),gms(lgms-1+ip7),
     &                gms(lgms-1+ip8),gms(lgms-1+ip0),gms(lgms-1+kfpd),
     &                gms(lgms-1+kma),gms(lgms-1+klcr),
     &                gms(lgms),gms(lgms))

                  irpold = irprim
                  nasc   = -1

*-----------------------------------------------------------------------
*           error after cali
*-----------------------------------------------------------------------

               if( ierror .eq. 1 .or. irprim .le. 0 ) then

                  if( ici .eq. 0 ) then

                     if( isrr .eq. 0 ) then

                        x = x + parz(28) * u
                        y = y + parz(28) * v
                        z = z + parz(28) * w

                        isrr = isrr + 1
                        goto 500

                     else

                        write(io,61) irprim,ierror
   61                   format(/' *** error message from s.cggsor ***'
     &                  /' illegal geometry condition was found.'
     &                  /' irprim =',i5,'  ierror =',i5)
                        call parastop( 827 )

                     end if

                  else

                     mark = -2
                     return

                  end if

               end if

*-----------------------------------------------------------------------

            iblz = nregno(gms(lgms),irprim)

*-----------------------------------------------------------------------

            if( ici .le. 0 ) then

                  nmed  = medno(gms(lgms),irprim)

                  if( nmed .eq. 0 )      nmed = -1
                  if( nmed .eq. kvlmax ) nmed =  0

            else if( ici .gt. 0 .and. iblz .ne. ih ) then

                  iblz2 = iblz
                  iblz1 = ih
                  iblz  = ih

                  mark = -3
                  return

            end if

*-----------------------------------------------------------------------
*        store initial region
*-----------------------------------------------------------------------

  100    continue

                  iblz1 = iblz
                  iblz2 = iblz

                  ilev1 = 0
                  ilev2 = 0

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine cggprp(icr,x,y,z,xc,yc,zc,u,v,w,
     &                  nmed,iblz,mark,markp)
*                                                                      *
*       CG geometry                                                    *
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
*       modified by K.Niita on 2001/11/30                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      real*4 gms(1)
      equivalence (das,gms)

*-----------------------------------------------------------------------

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /rcomon/ rcasc

*-----------------------------------------------------------------------

      common/flag/ierror
!$OMP THREADPRIVATE(/flag/)

      common/ark/ngy,nll,nbb,nzy(3),xd(3)
      common/arar/nby,nlev,nar,nq,iaw,iay,nf,nx1(3)
!$OMP THREADPRIVATE(/arar/)
      common/gomloc/kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     &              kkr2,knsr,kvoll,nadd,ldata,ltma,lfpd,numr,irtru,
     &              numb,nir,kbiz,kbcz
!$OMP THREADPRIVATE(/gomloc/)

      common/parem/xb(3),wb(3),wp(3),xp(3),rin,rout,pinf,dist,jr,idbg,
     &             irprim,nasc,lsurf,nbo,lri,lro,kloop,loop,itype,noa
!$OMP THREADPRIVATE(/parem/)

      common/cggstr/ lgms

      integer blzold
      common/orgi/dist0,markg,nmedg,nblz,blzold,irpold
!$OMP THREADPRIVATE(/orgi/)
      common/mgomv/mus,muz,ll,ipret,iflow,iect,nlo,igx
!$OMP THREADPRIVATE(/mgomv/)
      common/tape/intt,iot,iout,iou2,idm(4),ioecg
!$OMP THREADPRIVATE(/tape/)
      common/ss/s1
!$OMP THREADPRIVATE(/ss/)
      common/repeat/ip(20)

*-----------------------------------------------------------------------

      common/inout/in,io

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /regcrs/ icrsflx
      common /tlcost/ costha, uang(3), nsurf
!$OMP THREADPRIVATE(/tlcost/)

*-----------------------------------------------------------------------

      dimension xa(3), dc(3)
      dimension un(3)

*-----------------------------------------------------------------------
*           initial values
*-----------------------------------------------------------------------

            ierr = 0
            ierror = 0

            etath = dsqrt( (xc-x)**2 + (yc-y)**2 + (zc-z)**2 )

*-----------------------------------------------------------------------
*        check region : mark .ne. 0 and markp = 0
*-----------------------------------------------------------------------

            if( mark .ne. 0 .and. mark .ne. 2 .and. markp .eq. 0 ) then

                  ici = 1

               call cggsor(x,y,z,u,v,w,nmed,iblz,mark,markp,ici)

                  if( mark .le. -2 ) return

            end if

*-----------------------------------------------------------------------
*           for pilot
*-----------------------------------------------------------------------

               dist0 = dist + etath

               call reset(gms(lgms),gms(lgms),nl,nlu)

               ip4 = ip(4)
               ip5 = ip(5)
               ip7 = ip(7)
               ip8 = ip(8)
               ip0 = ip(10)

               xa(1) = xb(1)
               xa(2) = xb(2)
               xa(3) = xb(3)

*-----------------------------------------------------------------------

            call pilot(nl,nlu,xa,ip,
     &                 gms(lgms-1+ip4),gms(lgms-1+ip5),gms(lgms-1+ip7),
     &                 gms(lgms-1+ip8),gms(lgms-1+ip0),gms(lgms-1+kfpd),
     &                 gms(lgms-1+kma),gms(lgms-1+klcr),
     &                 gms(lgms),gms(lgms),distd)

*-----------------------------------------------------------------------

            if( irprim .le. 0 .or. ierror .ne. 0 ) then

               write(io,'(/''** ERROR: after pilot'')')
               goto 999

            end if

*-----------------------------------------------------------------------

               xc = x + u * s1
               yc = y + v * s1
               zc = z + w * s1

               nmed  = medno(gms(lgms),irprim)
               iblz2 = nregno(gms(lgms),irprim)
               mark  = markg

                  if( nmed .eq. 0 )      nmed = -1
                  if( nmed .eq. kvlmax ) nmed =  0
                  if( nmed .eq. -1 )     mark = -1

               iblz = iblz2

*-----------------------------------------------------------------------
*        for tally
*-----------------------------------------------------------------------

            if( markg .eq. 0 .and. icrsflx .ne. 0 ) then

               jrtmp = jr
               jr    = irpold

               call norml(gms(lgms-1+kma),gms(lgms-1+kfpd),
     &                    gms(lgms-1+klcr),gms(lgms-1+knbd),
     &                    gms(lgms-1+kkr1),gms(lgms-1+kkr2),un)

               jr = jrtmp

               xcd = xc - x
               ycd = yc - y
               zcd = zc - z

               absv = sqrt( xcd**2 + ycd**2 + zcd**2 )
               absu = sqrt( un(1)**2 + un(2)**2 + un(3)**2 )

               costha = ( un(1) * xcd + un(2) * ycd + un(3) * zcd )
     &                / ( absu * absv )

               costha = abs(costha)

               uang(1) = un(1) / absu
               uang(2) = un(2) / absu
               uang(3) = un(3) / absu

               nsurf = iabs( lsurf )

            end if

               irpold = irprim

*-----------------------------------------------------------------------
*     region error : mark = -4, -5
*-----------------------------------------------------------------------

            if( mark .eq. 1 .and. iblz1 .ne. iblz2 ) then

               mark = -4
               return

            else if( ( mark .eq. 0 .or. mark .eq. -1 )
     &               .and. ( iblz1 .eq. iblz2 ) ) then

               mark = -5
               return

            end if

*-----------------------------------------------------------------------
*     error in CG : lost particles : mark = -2
*-----------------------------------------------------------------------

      return

  999    continue

      if( icr .ne. 0 ) then

         write(io,'(/a,f16.0,i10)')
     &            '      ncas  nocas    =', rcasc, nocas

         if( nmed .gt. 0 ) then

            write(io,'( ''    mat ='',i5)') idmn(nmed)

         else

            write(io,'( ''    mat = void'')')

         end if

         write(io,'( ''   iblz1 ='',i5)') iblz1
         write(io,'( ''   iblz2 ='',i5)') iblz2
         write(io,'( ''    mark ='',i5)') mark

         write(io,'(/''    x, y, z  ='',1p3e17.9)') x,y,z
         write(io,'( ''    xc,yc,zc ='',1p3e17.9)') xc,yc,zc

         write(io,'(/a,1pe15.8,a,e15.8)')
     &                '    dist0 =',dist0,'  etath =',etath

         write(io,'(/73(''-''))')

         write(io,'(/a,i5,a,i5)') '  error condition: irprim=',irprim,
     &                          '  ierror=',ierror

      end if

*-----------------------------------------------------------------------

         mark = -2
         return

*-----------------------------------------------------------------------

      end

************************************************************************
*                                                                      *
      subroutine cggdis(dist,x,y,z,u,v,w,
     &                  nmed,iblz,mark,markp)
*                                                                      *
*       get distance up to boundary                                    *
*       modified by K.Niita on 2002/05/10                              *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /paraj/  mstz(300), parz(300)

*-----------------------------------------------------------------------

      save iblz0, nmed0, mark0
      save xc0, yc0, zc0
!$OMP THREADPRIVATE(iblz0, nmed0, mark0,xc0, yc0, zc0)
*-----------------------------------------------------------------------

            if( ( mark .eq. 0 .or. mark .eq. 2 ) .and.
     &            mstz(23) .eq. 0 ) then

               xc = x + u * parz(28)
               yc = y + v * parz(28)
               zc = z + w * parz(28)

               call cggprp(1,x,y,z,xc,yc,zc,u,v,w,
     &                     nmed,iblz,mark,markp)

               if( mark .ne. 1 ) goto 999

               x = x + u * parz(28)
               y = y + v * parz(28)
               z = z + w * parz(28)

               markp = 0

            end if

               iblz2 = iblz
               nmed2 = nmed

               xc = x + u * 1.0d+10
               yc = y + v * 1.0d+10
               zc = z + w * 1.0d+10

               call cggprp(1,x,y,z,xc,yc,zc,u,v,w,
     &                     nmed,iblz,mark,markp)

               if( mark .ne. 0 .and.
     &             mark .ne. 2 .and.
     &             mark .ne. -1 ) goto 999

               dist = sqrt( ( x - xc )**2
     &                    + ( y - yc )**2
     &                    + ( z - zc )**2 )

               xc0 = xc
               yc0 = yc
               zc0 = zc

               mark0 = mark
               iblz0 = iblz
               nmed0 = nmed

               iblz = iblz2
               nmed = nmed2

               mark  = 1
               markp = 0

*-----------------------------------------------------------------------

         return

*-----------------------------------------------------------------------

  999    continue

            if( mark .gt .-2 ) mark = -2

         return


************************************************************************
*                                                                      *
      entry cggnew(icr,dpr,x1,y1,z1,xc1,yc1,zc1,u1,v1,w1,e1,
     &             nmed1,iblz1,mark1,markp1)
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

            xc1 = xc0
            yc1 = yc0
            zc1 = zc0

            iblz1 = iblz0
            nmed1 = nmed0

            mark1  = mark0
            markp1 = 0

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine setcg(io,jo,ierr)
*                                                                      *
*       setup CG geometry                                              *
*       modified by K.Niita on 19/06/2000                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_character

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      real*4 gms(1)
      equivalence (das,gms)

*-----------------------------------------------------------------------

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      common /ccggg/  icgg
      common /cggmm/  ngstar, ngfini, ngfin0

      common /tcntl/  icntl, inucr
      common /inpec/  ititl, ipara, ibody, iregn, llarr, itby, itar
      common /inggs/  iog, igcel, ioa, igsuf, iob, igtrs

      common /regda/  ichl(kvlmax), chsm(kvlmax), ichmx, iod
      character       chsm*10
      common /regdu/  iuni(kvlmax)
      common /regdm/  idmg(kvlmax)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

      common/flag/ierror
!$OMP THREADPRIVATE(/flag/)

      common/cggstr/ lgms

      common/gomloc/kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     &              kkr2,knsr,kvoll,nadd,ldata,ltma,lfpd,numr,irtru,
     &              numb,nir,kbiz,kbcz
!$OMP THREADPRIVATE(/gomloc/)

      common /taliin/ rsouin, nzztin, nrgnin

      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200

      common/tape/inttcg,iotcg,ioutcg,iou2cg,idmggcg(4),ioecg
!$OMP THREADPRIVATE(/tape/)

*-----------------------------------------------------------------------

      character chtit*60

      dimension ipva(4)
      dimension bval(50)
      character chbd*3

      character chlw*200

      character chdf*80
      character cblan*200

      dimension iddm(10)

      dimension ibck(kvmmax)
      data ibck / kvmmax*0 /

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

               if( icgg .ne. 0 ) return

               ierr  = 0

*-----------------------------------------------------------------------

            ngstar = mmmax
            lgms   = mmmax * 2 + 1

*-----------------------------------------------------------------------

            mcmx = ( mdas  - 1 ) * 8 + 1
               call moddas_allocate_cha(MAX_NUM_CHRG, chrg)
               mci = 0

*-----------------------------------------------------------------------

         do i = 1, 200

            cblan(i:i) = ' '

         end do

*-----------------------------------------------------------------------
*     file = 20 is temporary file
*-----------------------------------------------------------------------

            iog = 20
            open(iog,form='formatted',status='scratch')

*-----------------------------------------------------------------------
*     open temporary file( 15 ) or MARSPF.IN for CG
*-----------------------------------------------------------------------

            iot = 15

         if( icntl .ne. 4 ) then

            open(iot,form='formatted',status='scratch')

         else if( icntl .eq. 4 ) then

            open(iot,file = chfn(4), form='formatted',status='unknown')

         end if

*-----------------------------------------------------------------------
*     write geometry for CG package on iot = 15 from itby = 18
*                                                    itar = 19
*-----------------------------------------------------------------------
*        write body infomation
*-----------------------------------------------------------------------

               rewind itby
               read(itby) chtit
               write(iot,'(a60)') chtit

               read(itby) ipva

               if( ipva(2) .gt. 1 ) then

                  ipv2 = 1

               else

                  ipv2 = 0

               end if

               write(iot,'(4i10)') ipva(1), ipv2, ipva(3), ipva(4)

            do k = 1, ibody

                  read(itby) chbd, ibnum, ibva, ( bval(i), i = 1, ibva )

                  ibck(k) = ibnum

               if( ipva(3) .eq. 0 ) then

                  if( ibva .eq. 4 ) then

                     write(iot,'(a3,3x,4(1p1e15.7))')
     &                    chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 8 ) then

                     write(iot,'(a3,3x,
     &                          4(1p1e15.7)/6x,4(1p1e15.7))')
     &                    chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 12 ) then

                     write(iot,'(a3,3x,
     &                          4(1p1e15.7)/6x,4(1p1e15.7),
     &                                     /6x,4(1p1e15.7))')
     &                    chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 30 ) then

                     write(iot,'(a3,3x,
     &                          4(1p1e15.7),5(/6x,4(1p1e15.7)))')
     &                    chbd, ( bval(i), i = 1, 24 )

                     write(iot,'(6x,6i10)')
     &                    ( nint( bval(i) ), i = 25, 30 )

                  else

                     write(iot,'(a3,3x,
     &                          4(1p1e15.7),5(/6x,4(1p1e15.7)))')
     &                    chbd, ( bval(i), i = 1, ibva )

                  end if

               else if( ipva(3) .gt. 0 ) then

                  if( ibva .eq. 4 ) then

                     write(iot,'(a3,i5,3x,4(1p1e15.7))')
     &                    chbd, ibnum, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 8 ) then

                     write(iot,'(a3,i5,3x,
     &                          4(1p1e15.7)/11x,4(1p1e15.7))')
     &                    chbd, ibnum, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 12 ) then

                     write(iot,'(a3,i5,3x,
     &                          4(1p1e15.7)/11x,4(1p1e15.7),
     &                                     /11x,4(1p1e15.7))')
     &                    chbd, ibnum, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 30 ) then

                     write(iot,'(a3,i5,3x,
     &                          4(1p1e15.7),5(/11x,4(1p1e15.7)))')
     &                    chbd, ibnum, ( bval(i), i = 1, 24 )

                     write(iot,'(11x,6i10)')
     &                    ( nint( bval(i) ), i = 25, 30 )

                  else

                     write(iot,'(a3,i5,3x,
     &                          4(1p1e15.7),5(/11x,4(1p1e15.7)))')
     &                    chbd, ibnum, ( bval(i), i = 1, ibva )

                  end if

               end if

            end do

                     write(iot,'(''end'')')

*-----------------------------------------------------------------------
*        write region definition
*-----------------------------------------------------------------------

                  irerr = 0
                  rewind iod

                  ioe = 22
                  open(ioe,status='scratch',form='unformatted')

            do i = 1, iregn

                  rewind ioe
                  read(iod) (chrg(mci+k:mci+k),k=1,ichl(i))

                  ild1 = 6
                  ild0 = 80 - ild1

                  ilrm = ichl(i)
                  isqd = 0
                  isrm = 0

  430             isqd = isqd + 1

                  if( ilrm .le. ild0 ) then

                     if( isqd .eq. 1 ) then

                           ildf = ilrm

                        do k = 1, ilrm

                           chdf(k:k) = chrg(mci+k+isrm:mci+k+isrm)

                        end do

                     else

                           write(ioe) ilrm
                           write(ioe)
     &                     (chrg(mci+k+isrm:mci+k+isrm),k=1,ilrm)

                     end if

                  else

                     do k = ild0, 1, -1

                        if( chrg(mci+k+isrm:mci+k+isrm) .eq. ' ' )
     &                  goto 420

                     end do

  420                k2 = k

                     if( isqd .eq. 1 ) then

                           ildf = k2 - 1

                        do k = 1, k2 - 1

                           chdf(k:k) = chrg(mci+k+isrm:mci+k+isrm)

                        end do

                     else

                           write(ioe) k2 - 1
                           write(ioe)
     &                     (chrg(mci+k+isrm:mci+k+isrm),k=1,k2-1)

                     end if

                     isrm = isrm + k2
                     ilrm = ilrm - k2

                     goto 430

                  end if

                     write(iot,'(a3,3x,200a1)') chsm(i),
     &                    (chdf(j:j),j=1,ildf)

               if( isqd .gt. 1 ) then

                     rewind ioe

                  do k = 2, isqd

                     read(ioe) ildf
                     read(ioe) (chdf(j:j),j=1,ildf)

                     write(iot,'(200a1)') (cblan(j:j),j=1,ild1),
     &                    (chdf(j:j),j=1,ildf)

                  end do

               end if

*-----------------------------------------------------------------------
*           check body number
*-----------------------------------------------------------------------

            ic = 0

  450       ic = ic + 1

            if( ic .gt. ichl(i) ) goto 460

               if( chrg(mci+ic:mci+ic) .eq. ' ' ) goto 450

               if( ( chrg(mci+ic:mci+ic) .eq. 'o' .or.
     &               chrg(mci+ic:mci+ic) .eq. 'O' ) .and.
     &             ( chrg(mci+ic+1:mci+ic+1) .eq. 'r' .or.
     &               chrg(mci+ic+1:mci+ic+1) .eq. 'R' ) ) then

                  ic = ic + 1
                  goto 450

               end if

                     call snum(chrg(mci+1:mci+1),ic,ichl(i),ic2,
     &                         cvvv,ierrt)

                     ic = ic2 - 1

                     ibnm = abs( nint( cvvv ) )

                     do j = 1, ibody

                        if( ibnm .eq. ibck(j) ) goto 450

                     end do

                  irerr = irerr + 1

                     goto 450

  460       continue

*-----------------------------------------------------------------------

            end do

                  close(ioe)

                  write(iot,'(''end'')')

                  if( irerr .ne. 0 ) goto 999

*-----------------------------------------------------------------------
*        write region number
*-----------------------------------------------------------------------

                  nrg = ( iregn - 1 ) / 10 + 1

               do i = 1, nrg

                  nj = ( i - 1 ) * 10 + 1
                  nk = min( i * 10, iregn )

                  write(iot,'(10i6)') ( idrg(j), j = nj, nk )

               end do

*-----------------------------------------------------------------------
*        write universe number
*-----------------------------------------------------------------------

               do i = 1, nrg

                  nj = ( i - 1 ) * 10 + 1
                  nk = min( i * 10, iregn )

                  write(iot,'(10i6)') ( iuni(j), j = nj, nk )

               end do

*-----------------------------------------------------------------------
*        write material number
*-----------------------------------------------------------------------

               do i = 1, nrg

                  nj = ( i - 1 ) * 10 + 1
                  nk = min( i * 10, iregn )

                  do j = nj, nk

                        k = j - nj + 1

                     if( idmg(j) .gt. 0 ) then

                        iddm(k) = idnm(idmg(j))

                     else if( idmg(j) .eq. 0 ) then

                        iddm(k) = kvlmax

                     else if( idmg(j) .eq. -1 ) then

                        iddm(k) = 0

                     end if

                  end do

                     kf = nk - nj + 1

                  write(iot,'(10i6)') ( iddm(j), j = 1, kf )

               end do

*-----------------------------------------------------------------------
*        write array information for CG package
*-----------------------------------------------------------------------

            if( llarr .gt. 0 ) then

               rewind itar

               do j = 1, llarr

                  read(itar,'(i6,200a1)')
     &                 i3, ( chlw(i:i), i = 1, i3 )
                  write(iot,'(200a1)') ( chlw(i:i), i = 1, i3 )

               end do

            else

                  write(iot,'(''     0'')')

            end if

*-----------------------------------------------------------------------
*        end for MARSPF out
*-----------------------------------------------------------------------

            if( icntl .eq. 4 ) then

               close(iot)
               return

            end if

*-----------------------------------------------------------------------
*     jomin reads all geometry data - combinatorial + array
*        read input data from iot = 15
*-----------------------------------------------------------------------

            rewind iot

            call scanof

*-----------------------------------------------------------------------
*        file = 16, 17 are temporary file for CG
*-----------------------------------------------------------------------

            n16 = 16
            n17 = 17

            open(n16,form='unformatted',status='scratch')
            open(n17,form='unformatted',status='scratch')

*-----------------------------------------------------------------------

            call jomin(gms(lgms),gms(lgms),1,ngeom,iot,iog,iog,n16,n17,
     &                 0,ndum,nmost,1,1,mdas)

               mmmax  = ngstar + ( ngeom + mod(ngeom,2) ) / 2
               ngfini = mmmax
               ngfin0 = mmmax

*-----------------------------------------------------------------------

         if( ngeom .gt. mdas ) then

               write(io,'(/''** Memory ERROR : at the end of CG,''/
     &                    ''** memory exceeds mdas'',/
     &                    ''** total memory: mdas ='',i9,/
     &                    ''** end   of GG memory ='',i9,/
     &                    ''** start of GG memory ='',i9,/
     &               ''<<<  Please extend mdas in param.inc >>>'')')
     &                      mdas, ngfini, ngstar

               ErrCha = ''
               MsgID = 'L:1099/R:setcg/F:marscg.f'
               call ErrWrite(MsgID, ErrCha)
               write(jo,'(/''** Memory ERROR : at the end of CG,''/
     &                    ''** memory exceeds mdas'',/
     &                    ''** total memory: mdas ='',i9,/
     &                    ''** end   of GG memory ='',i9,/
     &                    ''** start of GG memory ='',i9,/
     &               ''<<<  Please extend mdas in param.inc >>>'')')
     &                      mdas, ngfini, ngstar

               goto 999

         end if

         if( ierror .ne. 0 ) then

            write(io,1001)
            ErrCha = ''
            MsgID = 'L:1117/R:setcg/F:marscg.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1001)
1001        format(/' Error in setcg: job terminateted',
     &              ' check above lists')

            goto 999

         end if

         write(iog,'(/79(''-''))')

*-----------------------------------------------------------------------
*        close temporary files
*-----------------------------------------------------------------------

            close(iot)
            close(n16)
            close(n17)

*-----------------------------------------------------------------------

            call scanon
            call rcrdln(72,lenold)

*-----------------------------------------------------------------------
*     data for tally ( nrgn : number of geometry regions )
*-----------------------------------------------------------------------

         nrgnin = irtru

*-----------------------------------------------------------------------
*     check nrgnin and iregn
*-----------------------------------------------------------------------

      if( iregn .ne. nrgnin ) then

            write(io,'(''***ERROR: CG array, '',
     &                 ''Sorry array is not crrect.''/
     &                 '' iregn ='',i4,''  irtru ='',i4)')
     &            iregn, nrgnin

            ErrCha = ''
            MsgID = 'L:1160/R:setcg/F:marscg.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,'(''***ERROR: CG array, '',
     &                 ''Sorry array is not crrect.''/
     &                 '' iregn ='',i4,''  irtru ='',i4)')
     &            iregn, nrgnin

            goto 999

      end if

      goto 998

*-----------------------------------------------------------------------

  999 continue

            ierr = 1

*-----------------------------------------------------------------------

  998 continue

            if( me .ne. 0 ) then

                  close( iog )

            else

               if( ipva(2) .eq. 0 .and. ierr .eq. 0 ) then

                  close( iog )
                  open(iog,form='formatted',status='scratch')

               end if

                  endfile( iog )

            end if

*-----------------------------------------------------------------------
*     set file id to 6
*-----------------------------------------------------------------------

         inttcg  = 6
         iotcg   = 6
         ioutcg  = 6
         iout2cg = 6
         ioecg   = 6

*-----------------------------------------------------------------------

      call moddas_deallocate_cha(chrg)
      return
      end

