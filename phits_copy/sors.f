************************************************************************
*                                                                      *
      subroutine sors(iresample)
*                                                                      *
*        generate a source particle.                                   *
*        last modified by K.Niita on 2008/01/28                        *
*                                                                      *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use TETRAMOD, only : tetrasors !FURUTA20210324
      use mod_ompparallel
      use COSMICMOD ! T.Sato 2020/12/15
C for noncritical_ovly12 strict_rsouin
!$    use mod_rsouin
      use moddas_mesh
      use moddas_region
      use moddas_source

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'err.inc'
      include 'param-physcnst.inc' ! T.Sato 2023/07/23

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /mpi00/  npe, me
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /inout/  in,io
      common /tcntl/  icntl, inucr
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /iradkk/ randkk,irskip
      common /taliin/ rsouin, nzztin, nrgnin
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /eparm/  esmax, esmin, emin(20)
      common /ndemax/ dnmax(20)
      common /neulo/  reutn
      common /paraj/  mstz(300), parz(300)
      common /rcomon/ rcasc

      common /enginit/ engini
!$OMP THREADPRIVATE(/enginit/)

      common /srsph/  pmphs
!$OMP THREADPRIVATE(/srsph/)
      common /mttmc/  smttc(kvlmax), mttcn, mttc1(kvlmax), mttc2(kvlmax)
      common /wparm/  swtm(20), wc1(20), wc2(20) !FURUTA
      common /wparm0/ wc01(20), wc02(20)         !FURUTA
!$OMP THREADPRIVATE(/wparm0/)

      common /isorbia/ isbias(isrc) ! T.Sato, source bias method for xyz-mesh source

*-----------------------------------------------------------------------

      common /isocor/ iscorr, itcorr, imlwt(isrc)
      common /isomul/ smlwt(isrc), totfact, imsrc
      common /isorfs/ ispfs(isrc), rspfn, rspfz, ispfn
      common /isorst/ jstyp(isrc), istyp(isrc), inkf0(isrc), lstyp(isrc)
!$OMP THREADPRIVATE(/isorst/)
      common /isorsp/ sx0(isrc), sy0(isrc), sz0(isrc), sx1(isrc),
     &                sy1(isrc), sz1(isrc), sr0(isrc), se0(isrc),
     &                sdir(isrc), srx(isrc), sry(isrc), swem(isrc),
     &                sphi(isrc), sdom(isrc), swt0(isrc)
      common /isorsg/ seg0(isrc),seg1(isrc),seg2(isrc),seg3(isrc),
     &                set0(isrc),set1(isrc),set2(isrc),set3(isrc),
     &                jetyp(isrc),jptyp(isrc)
      common /isorsf/ isorf(isrc),lsfile(isrc), sfile(isrc)
      character sfile*100
      common /isorsm/ stm0(isrc),stmw(isrc),stmc(isrc),stmd(isrc),
     &                jttyp(isrc),jttpn(isrc)
      common /isorsn/ sr1(isrc), sr2(isrc), isrn(isrc)
      common /isorsr/ nsmx(isrc), nsrn(isrc), nsrc(isrc)
      common /isorsc/ isort(isrc,4), rsort(isrc,13)
      common /isorsd/ isdmp(isrc,0:30), jsdmp(isrc,0:30)
      common /isorsa/ sfactor(isrc)
      common /isorpn/ ssx(isrc), ssy(isrc), ssz(isrc)
!$OMP THREADPRIVATE(/isorpn/)
      common /isodct/ sdl0(isrc), sdl1(isrc), sdl2(isrc), sdpf(isrc),
     &                sdxw(isrc), sdyw(isrc), sdrd(isrc), sdebg(isrc),
     &                sdsxp(isrc), sdsxn(isrc), sdsyp(isrc),
     &                sdsyn(isrc), dnorm(isrc), prs1(isrc), prs2(isrc),
     &                prs3(isrc), prs4(isrc), psxp(isrc), psxn(isrc),
     &                psyp(isrc), psyn(isrc), sdls(isrc), sdrs(isrc),
     &                sdxs(isrc), sdys(isrc)

      common /isorsl/ nglp(isrc), ngli(isrc), ngla(isrc), ngfl(isrc),
     &                nglc(isrc), nglw(isrc)

      common /isousr/ jsusr(isrc,0:30)

      common /isoucn/ csphi,csth,snphi,snth

      common /isbeam/ sx2(isrc),sy2(isrc),
     &       sxmrad1(isrc),sxmrad2(isrc),symrad1(isrc),symrad2(isrc)   ! T.Sato, beam emittance source


*-----------------------------------------------------------------------

      common /isoras/ sag1(isrc),sag2(isrc),jatyp(isrc),jqtyp(isrc)
      common /isospg/ spg1(isrc),spg2(isrc) ! T.Sato 2021/02/26

*-----------------------------------------------------------------------

      common /isormd/ istdd(isrc),istcc(isrc),
     &                isxtp(isrc),isinx(isrc),istxx(isrc),
     &                isytp(isrc),isiny(isrc),istyy(isrc),
     &                isztp(isrc),isinz(isrc),istzz(isrc),
     &                sxmin(isrc),sxmax(isrc),sxdel(isrc),
     &                symin(isrc),symax(isrc),sydel(isrc),
     &                szmin(isrc),szmax(isrc),szdel(isrc)

*-----------------------------------------------------------------------

      common /isosuf/ issuf(isrc), iscut(isrc), isvct(isrc,8),
     &                issfd(isrc), isvfd(isrc,8),
     &                ivsfd(isrc), ivvfd(isrc,8),
     &                dvsfd(isrc,4), dvvfd(isrc,8,4)

      common /issufd/ isxdf(isrc), isydf(isrc), iszdf(isrc),
     &                isdef(isrc), icdef(isrc),
     &                ixdef(isrc,5), iydef(isrc,5), izdef(isrc,5),
     &                vxpos(isrc,5), vypos(isrc,5), vzpos(isrc,5),
     &                sxpos(isrc), sypos(isrc), szpos(isrc), ssrad(isrc)

*-----------------------------------------------------------------------

! T.Sato 2023/07/17, no cell region is regarded as lost particle
      common /cgerr/  nlost, ilost, igerr, icger, ncger, nrecover

      data ksou0 /0/
      save ksou0
!$OMP THREADPRIVATE(ksou0)  ! T.Sato 2020/12/26
      dimension ksou1(isrc)
      data ksou1 /isrc*0/
      save ksou1
!$OMP THREADPRIVATE(ksou1)
      dimension ksou1node(isrc)
      data ksou1node /isrc*0/ !FURUTA
      save ksou1node     !FURUTA

      data iserr /0/
      save iserr
!$OMP THREADPRIVATE(iserr)
      dimension     idas(1)
      equivalence ( das, idas )

      dimension dmpd(30)

      dimension pmlwt(isrc), pmlpr(isrc)
      save      pmlwt, pmlpr
!$OMP THREADPRIVATE(pmlwt, pmlpr)
      dimension idorr(isrc*100)
      save      idorr
!$OMP THREADPRIVATE(idorr)
      data kduct /0/
      save kduct
!$OMP THREADPRIVATE(kduct)
      data rsouin /0.0d0/ !FURUTA
*-----------------------------------------------------------------------
      common /bnkmem/ maxbnk, maxbn2, rtrckflp
      real(8) dmpmulti
      integer idmpmode,ibchjmp,idmpjmp
      common /stat2/ dmpmulti,idmpmode,ibchjmp,idmpjmp(2)
      integer nmulti
      real(8) rrr
*-----------------------------------------------------------------------
      common /itetsor/ ksoutnode,itetreg(isrc)
      real(8) x0(3)

C for noncritical_ovly12
      real(8), save :: rsouin_local
!$OMP THREADPRIVATE(rsouin_local)
      common /cntsrc/icntsrc(3,isrc) ! initial counter

! T.Sato 2020/12/07 for cosmic-ray source
      common /cosmicint/ipcosmic(isrc),icenv(isrc)
      common /cosmicani/ie511(isrc),annihratio(isrc) ! T.Sato 2023/07/23

      common /ibchsor/ ibcsn(isrc), ibsor(isrc,isrc)
      common /talout/ itall
      dimension bchson(isrc)
      integer, save :: ibfirst = 0
!$OMP THREADPRIVATE(ibfirst)
      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

      integer ipompl
      ipompl = ipomp+1

*-----------------------------------------------------------------------
*     initialization of normarization
*-----------------------------------------------------------------------
      iresample=0 ! T.Sato 2023/07/17 (=1) resampling is required
      iengbin=0 ! T.Sato 2020/12/06, (=1) determine energy first, then angle
      izeroe=0  ! T.Sato 2020/12/14, (=1) energy will be zero at the end

      if( ksou0 .eq. 0 ) then

            pmlpr(1) = 1.d0

         if( imsrc .gt. 1 .and. iscorr .eq. 0 ) then

            if( totfact .gt. 0 ) then

                  sek = 0.d0
               do i = 1, imsrc
                  sek = sek + smlwt(i)
               end do

                  pmlwt(1) = smlwt(1) / sek

               do i = 2, imsrc
                  pmlwt(i) = pmlwt(i-1) + smlwt(i) / sek
               end do

               do i = 1, imsrc
                  if( smlwt(i) .gt. 0.d0 ) then
                     pmlpr(i) = 1.d0
                  else
                     pmlpr(i) = 0.d0
                  end if
               end do

            else

                  sek = 0.d0
                  nek = 0
               do i = 1, imsrc
                  sek = sek + smlwt(i)
                  if( smlwt(i) .gt. 0.d0 ) nek = nek + 1
               end do

               do i = 1, imsrc
                  pmlpr(i) = smlwt(i) / sek * dble( nek )
               end do

                  sek = dble( nek )

                  if( smlwt(1) .gt. 0.d0 ) then
                     pmlwt(1) = 1.d0 / sek
                  else
                     pmlwt(1) = 0.d0
                  end if

               do i = 2, imsrc
                  if( smlwt(i) .gt. 0.d0 ) then
                     pmlwt(i) = pmlwt(i-1) + 1.d0 / sek
                  else
                     pmlwt(i) = pmlwt(i-1)
                  end if
               end do

            end if

*-----------------------------------------------------------------------

         else if( imsrc .gt. 1 .and. iscorr .ne. 0 ) then

               do i = 1, imsrc
                     pmlpr(i) = 1.d0
               end do

                  icj = 0
                  icm = 1

               do i = 1, itcorr

                     idorr(i) = icm
                     icj = icj + 1

                  if( icj .eq. imlwt(icm) ) then

                     icm = icm + 1
                     icj = 0

                  end if

               end do

         end if

         ksou0 = 1   ! S.Abe 2020/11/30

      end if

*-----------------------------------------------------------------------
*     normalization for itall = 4
*-----------------------------------------------------------------------
      if ( itall.eq.4 .and. nobch.gt.ibfirst ) then
        ibfirst = ibfirst + 1
        ibchson = 0
        do i = 1, imsrc
          bchson(i) = 1.d0
          if ( ibcsn(i).gt.0 ) then
            itt = 0
            do k = 1, ibcsn(i)
              if ( nobch.eq.ibsor(k,i) ) itt = itt + 1
            end do
            if ( itt.eq.0 ) then
              bchson(i) = 0.d0
            else
              ibchson = 1
            end if
          else
            ibchson = 1
          end if
        end do
        if ( ibchson.eq.0 ) then
!$OMP CRITICAL (nomultisource)
          write(ErrCha,'(''no multi-source is assigned to batch'',I4)')
     &    nobch
          ErrID = 'L:326/R:sors/F:sors.f'
          call ErrWrite(ErrID,ErrCha)
          call parastop( 999 )
!$OMP END CRITICAL (nomultisource)
        end if

        if( totfact .gt. 0 ) then
          sek = 0.d0
          do i = 1, imsrc
            sek = sek + bchson(i) * smlwt(i)
          end do
          pmlwt(1) = bchson(1) * smlwt(1) / sek
          do i = 2, imsrc
            pmlwt(i) = pmlwt(i-1) + bchson(i) * smlwt(i) / sek
          end do
          do i = 1, imsrc
            if( bchson(i)*smlwt(i) .gt. 0.d0 ) then
              pmlpr(i) = 1.d0
            else
              pmlpr(i) = 0.d0
            end if
          end do
        else
          sek = 0.d0
          nek = 0
          do i = 1, imsrc
            sek = sek + bchson(i) * smlwt(i)
            if( bchson(i)*smlwt(i) .gt. 0.d0 ) nek = nek + 1
          end do
          do i = 1, imsrc
            pmlpr(i) = bchson(i) * smlwt(i) / sek * dble( nek )
          end do
          sek = dble( nek )
          if( bchson(1)*smlwt(1) .gt. 0.d0 ) then
            pmlwt(1) = 1.d0 / sek
          else
            pmlwt(1) = 0.d0
          end if
          do i = 2, imsrc
            if( bchson(i)*smlwt(i) .gt. 0.d0 ) then
              pmlwt(i) = pmlwt(i-1) + 1.d0 / sek
            else
              pmlwt(i) = pmlwt(i-1)
            end if
          end do
        end if
      end if

*-----------------------------------------------------------------------
*     start of source
*-----------------------------------------------------------------------

        no    = 1
        nomax = 1
        ifirst = 0     !FURUTA20150515
C for NONCRITICAL_OVLY12
        rsouin_local = 0.0d0

*-----------------------------------------------------------------------

        nj = 0

 7000 continue

        nj = nj + 1
        nmulti = 1

        itetpos(nj,ipompl)=0   !FURUTA20241120
        itetposa(nj,ipompl)=0  !FURUTA20241120

        if(nj.gt.maxbnk)then
         write(6,'(/''**** Fatal ERROR: lack of bank'',
     &                       '' space ****''/
     &                       ''     please increase maxbnk'')')
         call parastop( 837 )
        endif
c------------------
*-----------------------------------------------------------------------
*     multi-source:  j
*           kduct : outer region of the duct source,
*                   repeat the same source
*-----------------------------------------------------------------------

 4000 continue

               j = 1

      if( imsrc .gt. 1 ) then

         if( kduct .eq. 0 ) then

            if( iscorr .eq. 0 ) then

                  prb = unirn(dummy)

               do k = 1, imsrc

                  if( pmlwt(k) .gt. prb ) goto 5000

               end do

                  k = imsrc

 5000             j = k

*-----------------------------------------------------------------------
*       multi correlation source
*-----------------------------------------------------------------------

            else if( iscorr .ne. 0 ) then

               if( nj .gt. itcorr ) then

                  nj = itcorr
                  goto 6000

               end if

                  j = idorr(nj)

            end if

*-----------------------------------------------------------------------

         else

               j = kduct

         end if

      end if

*-----------------------------------------------------------------------
*     dump source routine
* FURUTA20210726
*-----------------------------------------------------------------------
      if(jstyp(j).eq.17)then
       iduct=1
       jduct=1
       wtnorm=1.0d0
       wtduct=1.0d0
!$OMP CRITICAL (dump_sors)
       if(ibchjmp.ne.0)then ! No further dump source
        nmulti=0
       else
        call dumpsors(j,nj,nmulti,pmlpr,wtnorm,wtspfs,wtduct)
       endif
!$OMP END CRITICAL (dump_sors)
       goto 5899 ! T.Sato 2021/12/07 for considering angle data when u,v,w is not read from dump data
      endif

*-----------------------------------------------------------------------
*     number of try for finding region
*-----------------------------------------------------------------------

            itsrg = 0

 2000 continue

            iduct = 1
            jduct = 1

            wtnorm = 1.0d0
            wtduct = 1.0d0

*-----------------------------------------------------------------------

            t(ibkt+nj,ipompl) = 0.0d0

             ncnt(ibknct+1,nj,ipompl) = icntsrc(1,j) ! T.Sato 2019/05/19
             ncnt(ibknct+2,nj,ipompl) = icntsrc(2,j) ! T.Sato 2019/05/19
             ncnt(ibknct+3,nj,ipompl) = icntsrc(3,j) ! T.Sato 2019/05/19

             do itt =1, 3
               if ( nj.eq.1 ) then
                 nctmxsr(itt) = icntsrc(itt,j)
               else if ( icntsrc(itt,j).gt.nctmxsr(itt) ) then
                 nctmxsr(itt) = icntsrc(itt,j)
               end if
             end do

            name0 = 1

*-----------------------------------------------------------------------
*     surface sorce
*-----------------------------------------------------------------------

      if( jstyp(j) .eq. 26 ) then

               istry = 0
  111          istry = istry + 1
               if( istry .gt. nsmx(j) ) then
                  write(io,'(/'' Error : cut selection is failed,'',
     &                        '' in source.''/
     &                        '' max. trial ntmax = '',i5)') nsmx(j)
                  write(*,'(/'' Error : cut selection is failed,'',
     &                        '' in source.''/
     &                        '' max. trial ntmax = '',i5)') nsmx(j)
                  call parastop( 889 )
               end if

*-----------------------------------------------------------------------

         if( issfd(j) .eq. 2 ) then

               sxval = dvsfd(j,1) + parz(28)

            if( isdef(j) .eq. -1 .or. icdef(j) .eq. -1 ) then

                  thet = unirn(dummy) * 2.0 * pi
                    cs = cos(thet)
                    sn = sin(thet)
                    r  = unirn(dummy)
                    r  = ssrad(j) * dmax1( r, unirn(dummy) )
                  x(ibkx+nj,ipompl) = sxval
                  y(ibky+nj,ipompl) = r * cs + sypos(j)
                  z(ibkz+nj,ipompl) = r * sn + szpos(j)

               do k = 1, isydf(j)
                  if( iydef(j,k) .gt. 0 .and.
     &                y(ibky+nj,ipompl) .lt. vypos(j,k) ) goto 111
                  if( iydef(j,k) .lt. 0 .and.
     &                y(ibky+nj,ipompl) .gt. vypos(j,k) ) goto 111
               end do
               do k = 1, iszdf(j)
                  if( izdef(j,k) .gt. 0 .and.
     &                z(ibkz+nj,ipompl) .lt. vzpos(j,k) ) goto 111
                  if( izdef(j,k) .lt. 0 .and.
     &                z(ibkz+nj,ipompl) .gt. vzpos(j,k) ) goto 111
               end do

            else
                  ymin = vypos(j,1)
                  ymax = vypos(j,2)
                  zmin = vzpos(j,1)
                  zmax = vzpos(j,2)
                  x(ibkx+nj,ipompl) = sxval
                  y(ibky+nj,ipompl) =
     &                          unirn(dummy) * ( ymax - ymin ) + ymin
                  z(ibkz+nj,ipompl) =
     &                          unirn(dummy) * ( zmax - zmin ) + zmin

               if( icdef(j) .eq. 1 .or. isdef(j) .eq. 1 ) then
                  dist = sqrt( ( y(ibky+nj,ipompl) - sypos(j) )**2
     &                       + ( z(ibkz+nj,ipompl) - szpos(j) )**2 )
                  if( dist. le. ssrad(j) ) goto 111
               end if

            end if

               xcom = 1.0d0
               ycom = 0.0d0
               zcom = 0.0d0

*-----------------------------------------------------------------------
         else if( issfd(j) .eq. 3 ) then

               syval = dvsfd(j,1) + parz(28)

            if( isdef(j) .eq. -1 .or. icdef(j) .eq. -1 ) then

                  thet = unirn(dummy) * 2.0 * pi
                    cs = cos(thet)
                    sn = sin(thet)
                    r  = unirn(dummy)
                    r  = ssrad(j) * dmax1( r, unirn(dummy) )
                  x(ibkx+nj,ipompl) = r * cs + sxpos(j)
                  y(ibky+nj,ipompl) = syval
                  z(ibkz+nj,ipompl) = r * sn + szpos(j)

               do k = 1, isxdf(j)
                  if( ixdef(j,k) .gt. 0 .and.
     &                x(ibkx+nj,ipompl) .lt. vxpos(j,k) ) goto 111
                  if( ixdef(j,k) .lt. 0 .and.
     &                x(ibkx+nj,ipompl) .gt. vxpos(j,k) ) goto 111
               end do
               do k = 1, iszdf(j)
                  if( izdef(j,k) .gt. 0 .and.
     &                z(ibkz+nj,ipompl) .lt. vzpos(j,k) ) goto 111
                  if( izdef(j,k) .lt. 0 .and.
     &                z(ibkz+nj,ipompl) .gt. vzpos(j,k) ) goto 111
               end do

            else
                  xmin = vxpos(j,1)
                  xmax = vxpos(j,2)
                  zmin = vzpos(j,1)
                  zmax = vzpos(j,2)
                  x(ibkx+nj,ipompl) =
     &                         unirn(dummy) * ( xmax - xmin ) + xmin
                  y(ibky+nj,ipompl) = syval
                  z(ibkz+nj,ipompl) =
     &                         unirn(dummy) * ( zmax - zmin ) + zmin

               if( icdef(j) .eq. 1 .or. isdef(j) .eq. 1 ) then
                  dist = sqrt( ( x(ibkx+nj,ipompl) - sxpos(j) )**2
     &                       + ( z(ibkz+nj,ipompl) - szpos(j) )**2 )
                  if( dist. le. ssrad(j) ) goto 111
               end if

            end if

               xcom = 0.0d0
               ycom = 1.0d0
               zcom = 0.0d0

*-----------------------------------------------------------------------
         else if( issfd(j) .eq. 4 ) then

               szval = dvsfd(j,1) + parz(28)

            if( isdef(j) .eq. -1 .or. icdef(j) .eq. -1 ) then

                  thet = unirn(dummy) * 2.0 * pi
                    cs = cos(thet)
                    sn = sin(thet)
                    r  = unirn(dummy)
                    r  = ssrad(j) * dmax1( r, unirn(dummy) )
                  x(ibkx+nj,ipompl) = r * cs + sxpos(j)
                  y(ibky+nj,ipompl) = r * sn + sypos(j)
                  z(ibkz+nj,ipompl) = szval

               do k = 1, isxdf(j)
                  if( ixdef(j,k) .gt. 0 .and.
     &                x(ibkx+nj,ipompl) .lt. vxpos(j,k) ) goto 111
                  if( ixdef(j,k) .lt. 0 .and.
     &                x(ibkx+nj,ipompl) .gt. vxpos(j,k) ) goto 111
               end do
               do k = 1, isydf(j)
                  if( iydef(j,k) .gt. 0 .and.
     &                y(ibky+nj,ipompl) .lt. vypos(j,k) ) goto 111
                  if( iydef(j,k) .lt. 0 .and.
     &                y(ibky+nj,ipompl) .gt. vypos(j,k) ) goto 111
               end do

            else
                  xmin = vxpos(j,1)
                  xmax = vxpos(j,2)
                  ymin = vypos(j,1)
                  ymax = vypos(j,2)
                  x(ibkx+nj,ipompl) =
     &                        unirn(dummy) * ( xmax - xmin ) + xmin
                  y(ibky+nj,ipompl) =
     &                        unirn(dummy) * ( ymax - ymin ) + ymin
                  z(ibkz+nj,ipompl) = szval

               if( icdef(j) .eq. 1 .or. isdef(j) .eq. 1 ) then
                  dist = sqrt( ( x(ibkx+nj,ipompl) - sxpos(j) )**2
     &                       + ( y(ibky+nj,ipompl) - sypos(j) )**2 )
                  if( dist. le. ssrad(j) ) goto 111
               end if

            end if

               xcom = 0.0d0
               ycom = 0.0d0
               zcom = 1.0d0

*-----------------------------------------------------------------------
         else if( issfd(j) .ge. 5 .and. issfd(j) .le. 9 ) then

               ra = ssrad(j) + parz(28)

               cx = 1.0 - 2.0 * unirn(dummy)
               sx = sqrt( 1.0 - cx**2 )
               ph = unirn(dummy) * 2.0 * pi

               x(ibkx+nj,ipompl)  = sxpos(j) + ra * sx * cos( ph )
               y(ibky+nj,ipompl)  = sypos(j) + ra * sx * sin( ph )
               z(ibkz+nj,ipompl)  = szpos(j) + ra * cx

               xcom = sx * cos( ph )
               ycom = sx * sin( ph )
               zcom = cx

               do k = 1, isxdf(j)
                  if( ixdef(j,k) .gt. 0 .and.
     &                x(ibkx+nj,ipompl) .lt. vxpos(j,k) ) goto 111
                  if( ixdef(j,k) .lt. 0 .and.
     &                x(ibkx+nj,ipompl) .gt. vxpos(j,k) ) goto 111
               end do
               do k = 1, isydf(j)
                  if( iydef(j,k) .gt. 0 .and.
     &                y(ibky+nj,ipompl) .lt. vypos(j,k) ) goto 111
                  if( iydef(j,k) .lt. 0 .and.
     &                y(ibky+nj,ipompl) .gt. vypos(j,k) ) goto 111
               end do
               do k = 1, iszdf(j)
                  if( izdef(j,k) .gt. 0 .and.
     &                z(ibkz+nj,ipompl) .lt. vzpos(j,k) ) goto 111
                  if( izdef(j,k) .lt. 0 .and.
     &                z(ibkz+nj,ipompl) .gt. vzpos(j,k) ) goto 111
               end do

*-----------------------------------------------------------------------
         else if( issfd(j) .eq. 10 .or. issfd(j) .eq. 13 ) then

               ra = ssrad(j) + parz(28)

               thet = unirn(dummy) * 2.0 * pi
                 cs = cos(thet)
                 sn = sin(thet)

               x(ibkx+nj,ipompl) = unirn(dummy)
     &                    * ( vxpos(j,2) - vxpos(j,1) ) + vxpos(j,1)
               y(ibky+nj,ipompl) = ra * cs + sypos(j)
               z(ibkz+nj,ipompl) = ra * sn + szpos(j)

               xcom = 0.0
               ycom = cs
               zcom = sn

               do k = 1, isydf(j)
                  if( iydef(j,k) .gt. 0 .and.
     &                y(ibky+nj,ipompl) .lt. vypos(j,k) ) goto 111
                  if( iydef(j,k) .lt. 0 .and.
     &                y(ibky+nj,ipompl) .gt. vypos(j,k) ) goto 111
               end do
               do k = 1, iszdf(j)
                  if( izdef(j,k) .gt. 0 .and.
     &                z(ibkz+nj,ipompl) .lt. vzpos(j,k) ) goto 111
                  if( izdef(j,k) .lt. 0 .and.
     &                z(ibkz+nj,ipompl) .gt. vzpos(j,k) ) goto 111
               end do

*-----------------------------------------------------------------------
         else if( issfd(j) .eq. 11 .or. issfd(j) .eq. 14 ) then

               ra = ssrad(j) + parz(28)

               thet = unirn(dummy) * 2.0 * pi
                 cs = cos(thet)
                 sn = sin(thet)

               x(ibkx+nj,ipompl) = ra * cs + sxpos(j)
               y(ibky+nj,ipompl) = unirn(dummy)
     &                    * ( vypos(j,2) - vypos(j,1) ) + vypos(j,1)
               z(ibkz+nj,ipompl) = ra * sn + szpos(j)

               xcom = cs
               ycom = 0.0
               zcom = sn

               do k = 1, isxdf(j)
                  if( ixdef(j,k) .gt. 0 .and.
     &                x(ibkx+nj,ipompl) .lt. vxpos(j,k) ) goto 111
                  if( ixdef(j,k) .lt. 0 .and.
     &                x(ibkx+nj,ipompl) .gt. vxpos(j,k) ) goto 111
               end do
               do k = 1, iszdf(j)
                  if( izdef(j,k) .gt. 0 .and.
     &                z(ibkz+nj,ipompl) .lt. vzpos(j,k) ) goto 111
                  if( izdef(j,k) .lt. 0 .and.
     &                z(ibkz+nj,ipompl) .gt. vzpos(j,k) ) goto 111
               end do

*-----------------------------------------------------------------------
         else if( issfd(j) .eq. 12 .or. issfd(j) .eq. 15 ) then

               ra = ssrad(j) + parz(28)

               thet = unirn(dummy) * 2.0 * pi
                 cs = cos(thet)
                 sn = sin(thet)

               x(ibkx+nj,ipompl) = ra * cs + sxpos(j)
               y(ibky+nj,ipompl) = ra * sn + sypos(j)
               z(ibkz+nj,ipompl) = unirn(dummy)
     &                    * ( vzpos(j,2) - vzpos(j,1) ) + vzpos(j,1)

               xcom = cs
               ycom = sn
               zcom = 0.0

               do k = 1, isxdf(j)
                  if( ixdef(j,k) .gt. 0 .and.
     &                x(ibkx+nj,ipompl) .lt. vxpos(j,k) ) goto 111
                  if( ixdef(j,k) .lt. 0 .and.
     &                x(ibkx+nj,ipompl) .gt. vxpos(j,k) ) goto 111
               end do
               do k = 1, isydf(j)
                  if( iydef(j,k) .gt. 0 .and.
     &                y(ibky+nj,ipompl) .lt. vypos(j,k) ) goto 111
                  if( iydef(j,k) .lt. 0 .and.
     &                y(ibky+nj,ipompl) .gt. vypos(j,k) ) goto 111
               end do

         end if

*-----------------------------------------------------------------------

                  scosth = zcom
                  rt2 = xcom**2 + ycom**2
               if( rt2 .eq. 0.0d0 ) then
                  ssinth = 0.0d0
                  scosphi= 1.0d0
                  ssinphi= 0.0d0
               else
                  rt = sqrt(rt2)
                  ssinth  = rt
                  scosphi = xcom / rt
                  ssinphi = ycom / rt
               end if


*-----------------------------------------------------------------------
*        x(ibkx+1), y(ibky+1), z(ibkz+1), wtnorm
*-----------------------------------------------------------------------

      else if( jstyp(j) .eq. 1 .or. jstyp(j) .eq. 4 ) then

            wtnorm = swt0(j)

         if( sr0(j) .le. 0.0d0 ) then

            x(ibkx+nj,ipompl) = sx0(j)
            y(ibky+nj,ipompl) = sy0(j)
            z(ibkz+nj,ipompl) =
     &               unirn(dummy) * ( sz1(j) - sz0(j) ) + sz0(j)

         else

            thet = unirn(dummy) * 2.0 * pi
            cs   = cos(thet)
            sn   = sin(thet)

            if( sr0(j) .ge. sr1(j) ) then
  106          r = sr1(j) + ( sr0(j) - sr1(j) ) * unirn(dummy)
               if( unirn(dummy) .gt.  r / sr0(j) ) goto 106
            else
               r = unirn(dummy)
               r = sr0(j) * dmax1( r, unirn(dummy) )
            end if

            x(ibkx+nj,ipompl) = r * sn + sx0(j)
            y(ibky+nj,ipompl) = r * cs + sy0(j)
            z(ibkz+nj,ipompl) =
     &                   unirn(dummy) * ( sz1(j) - sz0(j) ) + sz0(j)

         end if

*-----------------------------------------------------------------------

      else if( jstyp(j) .eq. 2 .or. jstyp(j) .eq. 5 ) then

            wtnorm = swt0(j)

            x(ibkx+nj,ipompl) =
     &                     unirn(dummy) * ( sx1(j) - sx0(j) ) + sx0(j)
            y(ibky+nj,ipompl) =
     &                     unirn(dummy) * ( sy1(j) - sy0(j) ) + sy0(j)
            z(ibkz+nj,ipompl) =
     &                     unirn(dummy) * ( sz1(j) - sz0(j) ) + sz0(j)

*-----------------------------------------------------------------------

      else if( jstyp(j) .eq. 3 .or. jstyp(j) .eq. 6 ) then

            wtnorm = swt0(j)

            xsg = sx1(j) / 2.35482d+0
            ysg = sy1(j) / 2.35482d+0
            zsg = sz1(j) / 2.35482d+0

            x(ibkx+nj,ipompl) = xsg * gaurn(dummy) + sx0(j)
            y(ibky+nj,ipompl) = ysg * gaurn(dummy) + sy0(j)
            z(ibkz+nj,ipompl) = zsg * gaurn(dummy) + sz0(j)

*-----------------------------------------------------------------------

      else if( jstyp(j) .eq. 13 .or. jstyp(j) .eq. 14 ) then

            wtnorm = swt0(j)

            rsg = sr1(j) / 2.35482d+0

            rrx = gaurn(dummy)
            rry = gaurn(dummy)

            rrm = sqrt( rrx**2 + rry**2 )
            rrm = rsg * rrm

            thet = unirn(dummy) * 2.0 * pi
            cs   = cos(thet)
            sn   = sin(thet)

            x(ibkx+nj,ipompl) = rrm * cs + sx0(j)
            y(ibky+nj,ipompl) = rrm * sn + sy0(j)
            z(ibkz+nj,ipompl) =
     &               unirn(dummy) * ( sz1(j) - sz0(j) ) + sz0(j)

*-----------------------------------------------------------------------

      else if( jstyp(j) .eq. 7 .or. jstyp(j) .eq. 8 ) then

            wtnorm = swt0(j)

            ax = dble(isrn(j)+1)/2.0/dble(isrn(j))/sx1(j)**(isrn(j)+1)
            bx = dble(isrn(j)+1)/2.0/dble(isrn(j))/sx1(j)

            ay = dble(isrn(j)+1)/2.0/dble(isrn(j))/sy1(j)**(isrn(j)+1)
            by = dble(isrn(j)+1)/2.0/dble(isrn(j))/sy1(j)

            az = dble(isrn(j)+1)/2.0/dble(isrn(j))/sz1(j)**(isrn(j)+1)
            bz = dble(isrn(j)+1)/2.0/dble(isrn(j))/sz1(j)

  300       xr = ( 1.0 - 2.0 * unirn(dummy) ) * sx1(j)
            if( unirn(dummy) .lt. ax*xr**isrn(j)/bx ) goto 300

  301       yr = ( 1.0 - 2.0 * unirn(dummy) ) * sy1(j)
            if( unirn(dummy) .lt. ay*yr**isrn(j)/by ) goto 301

  302       zr = ( 1.0 - 2.0 * unirn(dummy) ) * sy1(j)
            if( unirn(dummy) .lt. az*zr**isrn(j)/bz ) goto 302

            x(ibkx+nj,ipompl) = xr + sx0(j)
            y(ibky+nj,ipompl) = yr + sy0(j)
            z(ibkz+nj,ipompl) = zr + sz0(j)

*-----------------------------------------------------------------------

      else if( jstyp(j) .eq. 15 .or. jstyp(j) .eq. 16 ) then

            wtnorm = swt0(j)

            ax = dble(isrn(j)+1)/2.0/dble(isrn(j))/sr1(j)**(isrn(j)+1)
            bx = dble(isrn(j)+1)/2.0/dble(isrn(j))/sr1(j)

            ay = dble(isrn(j)+1)/2.0/dble(isrn(j))/sr1(j)**(isrn(j)+1)
            by = dble(isrn(j)+1)/2.0/dble(isrn(j))/sr1(j)

  303       xr = ( 1.0 - 2.0 * unirn(dummy) ) * sr1(j)
            if( unirn(dummy) .lt. ax*xr**isrn(j)/bx ) goto 303

  304       yr = ( 1.0 - 2.0 * unirn(dummy) ) * sr1(j)
            if( unirn(dummy) .lt. ay*yr**isrn(j)/by ) goto 304

            rrm = sqrt( xr**2 + yr**2 )

            thet = unirn(dummy) * 2.0 * pi
            cs   = cos(thet)
            sn   = sin(thet)

            x(ibkx+nj,ipompl) = rrm * cs + sx0(j)
            y(ibky+nj,ipompl) = rrm * sn + sy0(j)
            z(ibkz+nj,ipompl) =
     &                     unirn(dummy) * ( sz1(j) - sz0(j) ) + sz0(j)

*-----------------------------------------------------------------------
*cKN 2011/04/05 corrected by Sato.

      else if( jstyp(j) .eq. 9 .or. jstyp(j) .eq. 10 ) then

            wtnorm = swt0(j)

            rh = sr2(j)**2

         if( rh .gt. 0.0 ) then
  100       ra = sr1(j) + ( sr2(j) - sr1(j) ) * unirn(dummy)
            if( unirn(dummy) .gt. ra**2 / rh ) goto 100
         else
            ra = sr1(j)
         end if

            cx = 1.0 - 2.0 * unirn(dummy)
            sx = sqrt( 1.0 - cx**2 )
            ph = unirn(dummy) * 2.0 * pi

            x(ibkx+nj,ipompl)  = sx0(j) + ra * sx * cos( ph )
            y(ibky+nj,ipompl)  = sy0(j) + ra * sx * sin( ph )
            z(ibkz+nj,ipompl)  = sz0(j) + ra * cx

c *********************************************************************
c Modified by T.Sato for isotopic irradiation without using weight 2012/4/4
c *********************************************************************
         if( sdir(j) .le. -2.5 ) then
 1001       xd = (unirn(dummy)-0.5) * 2.0 * sr2(j)
            yd = (unirn(dummy)-0.5) * 2.0 * sr2(j)

            if (xd**2 + yd**2 .gt. sr2(j)**2) goto 1001

            zd = sr1(j)

c 2012/04/03 by T.Sato, revised 2020/12/14
c           cx is cosine of zenith angle, phi is azimuthal angle
            if(icenv(j).ge.1) then ! terrestrial cosmic-ray mode
             jin = j
             call engbin(eo,ew0,jin,ie) ! energy must be determined before determining angle
             if(ipcosmic(j).eq.33.and.ie.eq.ie511(j)) then
              if(unirn(dummy).gt.1.0d0/annihratio(j))
     &        eo=rstms(12)*1.0d3
             endif
             iengbin=1
             cx = getCosmicAng(jin,ie)
            else  ! isotropic distribution
             cx = 1.0 - 2.0 * unirn(dummy)
            endif
            sx = sqrt( 1.0 - cx**2 )
            phi  =2.0*pi*unirn(dummy)

            x(ibkx+nj,ipompl) = xd*cx*cos(phi) - yd*sin(phi)
     &                 + zd*sx*cos(phi) + sx0(j)
            y(ibky+nj,ipompl) = xd*cx*sin(phi) + yd*cos(phi)
     &                 + zd*sx*sin(phi) + sy0(j)
            z(ibkz+nj,ipompl) =-xd*sx + zd*cx + sz0(j)

            u(ibku+nj,ipompl) = -sx*cos(phi)
            v(ibkv+nj,ipompl) = -sx*sin(phi)
            w(ibkw+nj,ipompl) = -cx

            if(w(ibkw+nj,ipompl).lt.sag1(j).or.
     &      w(ibkw+nj,ipompl).gt.sag2(j)) izeroe=1  ! zero-energy index

            phideg=phi*180.0d0/acos(-1.0d0) ! convert rad to degree
            if(spg1(j).lt.spg2(j)) then ! normal case
             if(phideg.lt.spg1(j).or.phideg.gt.spg2(j)) izeroe=1 ! zero-energy index, T.Sato 2021/02/26
            else  ! inverse case
             if(phideg.lt.spg1(j).and.phideg.gt.spg2(j)) izeroe=1 ! zero-energy index
            endif

            if(izeroe.eq.1.and.isbias(j).eq.1) then ! do not allow to produce dead particle, T.Sato 2021/03/03
             izeroe=0
             goto 1001
            endif

! 2017/11/21 by T.Sato to move source location on the sphere
            sintmp=sqrt(xd**2+yd**2)/sr1(j)
            costmp=sqrt(1.0-sintmp**2)
            dtmp = sr1(j)*(1.0-costmp)

            x(ibkx+nj,ipompl) =
     &                   x(ibkx+nj,ipompl) + dtmp*u(ibku+nj,ipompl)
            y(ibky+nj,ipompl) =
     &                   y(ibky+nj,ipompl) + dtmp*v(ibkv+nj,ipompl)
            z(ibkz+nj,ipompl) =
     &                   z(ibkz+nj,ipompl) + dtmp*w(ibkw+nj,ipompl)

         end if
c **************  End of Modification by T.Sato  ***********************

*-----------------------------------------------------------------------

      else if( jstyp(j) .eq. 11 .or. jstyp(j) .eq. -11) then ! T.Sato 2020/11/29, s-type=-11 is used for e-type

            wtnorm = swt0(j)

      if(swem(j).ne.0) then ! original mode

  101    continue

            xp =  1. - 2. * unirn(dummy)
            xq =  1. - 2. * unirn(dummy)

            if( xp**2 + xq**2 .gt. 1.0 ) goto 101

  102    continue

            yp =  1. - 2. * unirn(dummy)
            yq =  1. - 2. * unirn(dummy)

            if( yp**2 + yq**2 .gt. 1.0 ) goto 102

            xpo = sqrt( swem(j) * sx1(j) ) * xp
            xqo = sqrt( swem(j) / sx1(j) ) * xq
            ypo = sqrt( swem(j) * sy1(j) ) * yp
            yqo = sqrt( swem(j) / sy1(j) ) * yq

      else ! gaussian distribution, T.Sato 2020/09/28
            xpo = sx1(j) * gaurn(dummy)
            xqo = sxmrad1(j) * gaurn(dummy)
            ypo = sy1(j) * gaurn(dummy)
            yqo = symrad1(j) * gaurn(dummy)
      endif

            xcm = xpo * cos(srx(j)) - xqo * sin(srx(j)) + sx2(j)
            xpm = xpo * sin(srx(j)) + xqo * cos(srx(j)) + sxmrad2(j)
            ycm = ypo * cos(sry(j)) - yqo * sin(sry(j)) + sy2(j)
            ypm = ypo * sin(sry(j)) + yqo * cos(sry(j)) + symrad2(j)

            x(ibkx+nj,ipompl) = sx0(j) + xcm
            y(ibky+nj,ipompl) = sy0(j) + ycm
            z(ibkz+nj,ipompl) =
     &                 unirn(dummy) * ( sz1(j) - sz0(j) ) + sz0(j)

            vz = sdir(j)
            vx = vz * tan( xpm / 1000.0 )
            vy = vz * tan( ypm / 1000.0 )

            av = sqrt( 1.0 + vx**2 + vy**2 )

            u(ibku+nj,ipompl) = vx / av
            v(ibkv+nj,ipompl) = vy / av
            w(ibkw+nj,ipompl) = vz / av

*-----------------------------------------------------------------------

      else if( jstyp(j) .eq. 12 ) then

         if( ksou1(j) .eq. 0 .and. ksou1node(j) .eq. 0 ) then !FURUTA

            if( irskip .ne. 0 ) then

                  jrskip = iabs( irskip )

               do i = 1, jrskip

  887             continue

                  read(isorf(j),*,iostat=ios)
     &                 xp, xq, yp, yq, eo, wt000, pz000
                  if( ios .eq. -1 ) goto 888

                  goto 889
  888             rewind isorf(j)
                  goto 887
  889             continue

               end do

            end if

         end if

  987             continue

                  read(isorf(j),*,iostat=ios)
     &                 xp, xq, yp, yq, eo, wt000, pz000
                  if( ios .eq. -1 ) goto 988

                  goto 989
  988             rewind isorf(j)
                  goto 987
  989             continue

            x(ibkx+nj,ipompl) = xp + sx0(j)
            y(ibky+nj,ipompl) = yp + sy0(j)
            z(ibkz+nj,ipompl) = sz0(j)

            vz = sdir(j)
            vx = vz * tan( xq / 1000 )
            vy = vz * tan( yq / 1000 )

            av = sqrt( 1.0 + vx**2 + vy**2 )

            u(ibku+nj,ipompl) = vx / av
            v(ibkv+nj,ipompl) = vy / av
            w(ibkw+nj,ipompl) = vz / av

            e(ibke+nj,ipompl) = sqrt((938.28)**2 + (eo*1000.0)**2)
     &                                          - 938.28
            wtnorm = wt000


*-----------------------------------------------------------------------

      else if( jstyp(j) .eq. 18 .or. jstyp(j) .eq. 19 ) then


               znor = sqrt( sx1(j)**2 + sy1(j)**2 + sz1(j)**2 )

               if( znor .eq. 0.0d0 ) then
                  write(ErrCha,*) ' (sx1, sy1, sz1) should not be zero'
                  ErrID = 'L:1188/R:sors/F:sors.f' !E03_104_001
                  call ErrWrite(ErrID,ErrCha)
                  stop
               end if

                  xx0 = sx1(j) / znor
                  yy0 = sy1(j) / znor
                  zz0 = sz1(j) / znor

                  csth = zz0
                  rt2 = xx0**2 + yy0**2

               if( rt2 .eq. 0.0d0 ) then

                  snth = 0.0d0
                  csphi= 1.0d0
                  snphi= 0.0d0

               else

                  rt = sqrt(rt2)
                  snth  = rt
                  csphi = xx0 / rt
                  snphi = yy0 / rt

               end if


  114          rrr = ( sr1(j) - sr0(j) ) * unirn(dummy)
               if( sr1(j) * unirn(dummy) .gt. sr0(j) + rrr ) goto 114
               rr0 = rrr + sr0(j)

               rr1 = rr0 * sin( sr2(j) * pi / 180.d0 )
               zz0 = rr0 * cos( sr2(j) * pi / 180.d0 )

               thet = unirn(dummy) * 2.0 * pi
               xx0 = rr1 * cos(thet)
               yy0 = rr1 * sin(thet)

               pab = sqrt( xx0**2 + yy0**2 + zz0**2 )

               if( pab .gt. 0.0d0 ) then

                  adc = xx0 / pab
                  bdc = yy0 / pab
                  gdc = zz0 / pab

               else

                  adc = 0.0
                  bdc = 0.0
                  gdc = 1.0

               end if

                  t1  = csth  * adc + snth  * gdc
                  xx1 = csphi * t1  - snphi * bdc
                  yy1 = snphi * t1  + csphi * bdc
                  zz1 = csth  * gdc - snth  * adc

            x(ibkx+nj,ipompl) = xx1 * pab + sx0(j)
            y(ibky+nj,ipompl) = yy1 * pab + sy0(j)
            z(ibkz+nj,ipompl) = zz1 * pab + sz0(j)

            wtnorm = swt0(j)

*-----------------------------------------------------------------------

      else if( jstyp(j) .eq. 20 .or. jstyp(j) .eq. 21 ) then

            v1x = sx1(j)
            v1y = sy1(j)
            v1z = sz1(j)

            v2x = sr0(j)
            v2y = sr1(j)
            v2z = sr2(j)

            v3x = srx(j)
            v3y = sry(j)
            v3z = swem(j)

            v3h = sqrt( v3x**2 + v3y**2 + v3z**2 )

            vtx = v1x + v2x
            vty = v1y + v2y
            vtz = v1z + v2z
            vtl = sqrt( vtx**2 + vty**2 + vtz**2 )

            rp1 = unirn(dummy)
            rp2 = unirn(dummy)

         if( rp2 .gt. 1.d0 - rp1 ) then

            rp1 = 1.d0 - rp1
            rp2 = 1.d0 - rp2

         end if

            ptx = ( v1x * rp1 + v2x * rp2 )
            pty = ( v1y * rp1 + v2y * rp2 )
            ptz = ( v1z * rp1 + v2z * rp2 )


         if( sdl0(j) .gt. 0.d0 ) then

  592       rp3 = exprnf(dummy) / sdl0(j)
            if( rp3 .gt. v3h ) goto 592

         else

cKN 2017/02/22 typo!!!
            rp3 = v3h * unirn(dummy)

         end if

            x(ibkx+nj,ipompl) = ptx + rp3 * v3x / v3h + sx0(j)
            y(ibky+nj,ipompl) = pty + rp3 * v3y / v3h + sy0(j)
            z(ibkz+nj,ipompl) = ptz + rp3 * v3z / v3h + sz0(j)

            wtnorm = swt0(j)


*-----------------------------------------------------------------------

      else if( jstyp(j) .eq. 22 .or. jstyp(j) .eq. 23 ) then

            wtnorm = swt0(j)

            istdg=istdd(j)
            istdc=istcc(j)

         if(isbias(j).eq.0) then ! normal bias mode

            random = unirn(dummy)

            inxyz=iabs(isinz(j)*isiny(j)*isinx(j))
            call dichotomy(klm,inxyz,random,das(istdc+1))
            k=(klm-1)/iabs(isinx(j)*isiny(j))+1
            l=(klm-1-(k-1)*iabs(isinx(j)*isiny(j)))/iabs(isinx(j))+1
            m=klm-(k-1)*iabs(isinx(j)*isiny(j))-(l-1)*iabs(isinx(j))

! 123        continue

         elseif(isbias(j).eq.1) then  ! Weight bias mode, T.Sato 2020/08/09

            do i=1,100000
               k=int(unirn(dummy)*iabs(isinz(j)))+1
               l=int(unirn(dummy)*iabs(isiny(j)))+1
               m=int(unirn(dummy)*iabs(isinx(j)))+1
               klm=(k-1)*iabs(isinx(j)*isiny(j))+(l-1)*iabs(isinx(j))+m
               if(das(istdg+klm).ne.0.0) then  ! data exist
                  wtnorm = wtnorm*das(istdg+klm)
                  exit
               endif
            enddo

         else  ! equal distribution mode, T.Sato 2020/08/09

!$OMP CRITICAL (stype22_crit)
            inxyz=iabs(isinz(j)*isiny(j)*isinx(j))
            random = unirn(dummy)*das(istdc+inxyz)
            call dichotomy(klm,inxyz,random,das(istdc+1))
            k=(klm-1)/iabs(isinx(j)*isiny(j))+1
            l=(klm-1-(k-1)*iabs(isinx(j)*isiny(j)))/iabs(isinx(j))+1
            m=klm-(k-1)*iabs(isinx(j)*isiny(j))-(l-1)*iabs(isinx(j))


            wtnorm = wtnorm*das(istdg+klm)

            do iklm=klm,inxyz
               das(istdc+iklm)=das(istdc+iklm)-1.0d0
            enddo

            if(das(istdc+inxyz).lt.0.01) then ! no data anymore, reset das(istdc)
                  if(das(istdg+1).eq.0.0) then
                     das(istdc+1) = 0.0
                  else
                     das(istdc+1) = 1.0d0
                  endif

                  do k1 = 2, inxyz

                     if(das(istdg+k1).eq.0.0) then
                        das(istdc+k1)=das(istdc+k1-1)
                     else
                        das(istdc+k1)=das(istdc+k1-1)+1.0d0
                     endif

                  end do

            endif

!$OMP END CRITICAL (stype22_crit)

         endif

            rx1 = das_istxx(istxx(j)+m-1)
            rx2 = das_istxx(istxx(j)+m)
            ry1 = das_istyy(istyy(j)+l-1)
            ry2 = das_istyy(istyy(j)+l)
            rz1 = das_istzz(istzz(j)+k-1)
            rz2 = das_istzz(istzz(j)+k)

            x(ibkx+nj,ipompl) = unirn(dummy) * ( rx2 - rx1 ) + rx1
            y(ibky+nj,ipompl) = unirn(dummy) * ( ry2 - ry1 ) + ry1
            z(ibkz+nj,ipompl) = unirn(dummy) * ( rz2 - rz1 ) + rz1

*-----------------------------------------------------------------------

      else if( jstyp(j) .eq. 17 ) then
cFURUTA20210726----
       ErrCha='Something wrong. jstyp(j)=17 was detected'
       ErrID = 'L:1401/R:sors/F:sors.f'
       call ErrWrite(ErrID,ErrCha)

*-----------------------------------------------------------------------
*     TETRA source : jstyp(j) = 24 mono = 25 spectra
*-----------------------------------------------------------------------
      else if( jstyp(j) .eq. 24 .or. jstyp(j) .eq. 25 ) then

cFURUTA20210324 move to sors_init

         do
          icl=idgr(itetreg(j))
          call tetrasors(icl,x0,ierr)
          if(ierr.gt.0)then

           ErrCha = ''
           ErrID = 'L:1417/R:sors/F:sors.f' !E80_002_001
           call ErrWrite(ErrID,ErrCha)

           write(*,'(''*** TETRA SOURCE ERROR: '',
     &          ''source REG is not found in tetrahedrons'')')
           write(*,'(''itetreg ='',i7)')itetreg(j)
           call parastop( 602 )
          endif
          x(ibkx+nj,ipompl) = x0(1)
          y(ibky+nj,ipompl) = x0(2)
          z(ibkz+nj,ipompl) = x0(3)
          ici=0
          mark=1
          markp=0
          call gomsor(x(ibkx+nj,ipompl),y(ibky+nj,ipompl),
     &                z(ibkz+nj,ipompl),
     &                  u(ibku+nj,ipompl),v(ibkv+nj,ipompl),
     &                  w(ibkw+nj,ipompl),
     &                  nmed(ibknmd+nj,ipompl),iblz(ibkblz+nj,ipompl),
     &                  mark,markp,ici)
          if(iblz(ibkblz+nj,ipompl).eq.itetreg(j))exit
          itsrg=itsrg+1
          if(itsrg.gt.nsmx(j))then

          ErrCha = ''
          ErrID = 'L:1442/R:sors/F:sors.f' !E80_002_002
          call ErrWrite(ErrID,ErrCha)

           write(io,'(/''*** TETRA SOURCE ERROR: '',
     &          ''region selection is failed in source.''/
     &                     '' max. trial ntmax = '',i5)') nsmx(j)
           write(io,'(''itetreg ='',i7)')itetreg(j)

           write(*,'(''*** TETRA SOURCE ERROR: '',
     &          ''region selection is failed in source.''/
     &                     '' max. trial ntmax = '',i5)') nsmx(j)
           write(*,'(''itetreg ='',i7)')itetreg(j)

               call parastop( 889 )

          endif
         enddo

*-----------------------------------------------------------------------
*     user defined source : jstyp(j) = 100
*-----------------------------------------------------------------------

      else if( jstyp(j) .eq. 100 ) then

         if( ksou1(j) .eq. 0 .and. ksou1node(j) .eq. 0 ) then !FURUTA

            if( irskip .ne. 0 ) then

                  jrskip = iabs( irskip )

               do i = 1, jrskip

                  call usrsors(ux,uy,uz,uu,uv,uw,ue,uwt,utim,inam,iukf,
     &                         int1,int2,int3,usx,usy,usz)

               end do

            end if

         end if
            call initusrsors(ux,uy,uz,uu,uv,uw,ue,uwt,utim,inam,iukf,
     &        int1,int2,int3,usx,usy,usz) !FURUTA20200206
            call usrsors(ux,uy,uz,uu,uv,uw,ue,uwt,utim,inam,iukf,
     &                   int1,int2,int3,usx,usy,usz)
*-----------------------------------------------------------------------

            if( jsusr(j,1) .eq. 0 ) then

               inkf0(j) = iukf
               istyp(j) = kftp( iukf )

            end if

*-----------------------------------------------------------------------

            if( jsusr(j,2) .gt. 0 ) then
               ux = unirn(dummy) * ( sx1(j) - sx0(j) ) + sx0(j)
            end if

            if( jsusr(j,3) .gt. 0 ) then
               uy = unirn(dummy) * ( sy1(j) - sy0(j) ) + sy0(j)
            end if

            if( jsusr(j,4) .gt. 0 ) then
               uz = unirn(dummy) * ( sz1(j) - sz0(j) ) + sz0(j)
            end if

               x(ibkx+nj,ipompl) = ux
               y(ibky+nj,ipompl) = uy
               z(ibkz+nj,ipompl) = uz

*-----------------------------------------------------------------------

            if( sdir(j) .gt. 1000.0 ) then

               u(ibku+nj,ipompl) = uu
               v(ibkv+nj,ipompl) = uv
               w(ibkw+nj,ipompl) = uw

            end if

*-----------------------------------------------------------------------

            if( jsusr(j, 8) .eq. 0 ) se0(j) = ue
            if( jsusr(j, 9) .eq. 0 ) wtnorm = uwt
            if( jsusr(j,10) .eq. 0 ) t(ibkt+nj,ipompl) = utim

            if(ncnt(ibknct+1,nj,ipompl).eq.0)
     &              ncnt(ibknct+1,nj,ipompl) = int1 ! T.Sato 2019/06/25, use initial value if defined
            if(ncnt(ibknct+2,nj,ipompl).eq.0)
     &              ncnt(ibknct+2,nj,ipompl) = int2 ! T.Sato 2019/06/25, use initial value if defined
            if(ncnt(ibknct+3,nj,ipompl).eq.0)
     &              ncnt(ibknct+3,nj,ipompl) = int3 ! T.Sato 2019/06/25, use initial value if defined

            if( jsusr(j,14) .eq. 0 ) ssx(j) = usx
            if( jsusr(j,15) .eq. 0 ) ssy(j) = usy
            if( jsusr(j,16) .eq. 0 ) ssz(j) = usz
            if( jsusr(j,17) .eq. 0 ) name0  = inam

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*     u(ibku+1), v(ibkv+1), w(ibkw+1)
*-----------------------------------------------------------------------
 5899 continue ! T.Sato 2021/12/07 for determining source direction when dir is specified
            aw0 = 1.d0

      if( abs(jstyp(j)) .eq. 11 .or. jstyp(j) .eq. 12 ) then ! T.Sato 2021/01/01

*-----------------------------------------------------------------------
*     Spontaneus Fission Source
*-----------------------------------------------------------------------

      else if( ispfs(j) .ne. 0 ) then

            cx = 1.0 - 2.0 * unirn(dummy)
            sx = sqrt( 1.0 - cx**2 )
            ph = unirn(dummy) * 2.0 * pi

            u(ibku+nj,ipompl)  = sx * cos( ph )
            v(ibkv+nj,ipompl)  = sx * sin( ph )
            w(ibkw+nj,ipompl)  = cx

*-----------------------------------------------------------------------

      else if( jstyp(j) .eq. 9 .or. jstyp(j) .eq. 10) then

*-----------------------------------------------------------------------
       if(sdir(j) .gt. -2.5) then  ! T.Sato 2024/09/07 add condition for excluding iso source because source condition has alrady been determined

         if( sdir(j) .gt. 100.0 .and. sdir(j). lt. 250.0) then ! all

            cx = 1.0 - 2.0 * unirn(dummy)
            sx = sqrt( 1.0 - cx**2 )
            ph = unirn(dummy) * 2.0 * pi

            u(ibku+nj,ipompl)  = sx * cos( ph )
            v(ibkv+nj,ipompl)  = sx * sin( ph )
            w(ibkw+nj,ipompl)  = cx


         elseif( sdir(j) .eq. 1.0d0 ) then  ! dir = 1.0

            u(ibku+nj,ipompl)  = ( x(ibkx+nj,ipompl) - sx0(j) ) / ra
            v(ibkv+nj,ipompl)  = ( y(ibky+nj,ipompl) - sy0(j) ) / ra
            w(ibkw+nj,ipompl)  = ( z(ibkz+nj,ipompl) - sz0(j) ) / ra

         else if( sdir(j) .eq. -1.0d0 ) then  ! dir = -1.0

            u(ibku+nj,ipompl)  = -( x(ibkx+nj,ipompl) - sx0(j) ) / ra
            v(ibkv+nj,ipompl)  = -( y(ibky+nj,ipompl) - sy0(j) ) / ra
            w(ibkw+nj,ipompl)  = -( z(ibkz+nj,ipompl) - sz0(j) ) / ra

         else if( sdir(j) .ge. 1.5 .and. sdir(j) .lt. 2.5 ) then ! +all

            ux    = ( x(ibkx+nj,ipompl) - sx0(j) ) / ra
            uy    = ( y(ibky+nj,ipompl) - sy0(j) ) / ra
            costh = ( z(ibkz+nj,ipompl) - sz0(j) ) / ra
            rt2   = ux**2 + uy**2

            if( rt2 .eq. 0.0d0 ) then

               sinth = 0.0d0
               cosphi= 1.0d0
               sinphi= 0.0d0

            else

               rt = sqrt(rt2)
               sinth  = rt
               cosphi = ux / rt
               sinphi = uy / rt

            end if

  103       cx = unirn(dummy)
            if( unirn(dummy) .gt. cx ) goto 103

            sx = sqrt( 1.0 - cx**2 )
            ph = unirn(dummy) * 2.0 * pi

            adc = sx * cos( ph )
            bdc = sx * sin( ph )
            gdc = cx

            t1   = costh  * adc + sinth  * gdc
            adcc = cosphi * t1  - sinphi * bdc
            bdcc = sinphi * t1  + cosphi * bdc
            gdcc = costh  * gdc - sinth  * adc

            u(ibku+nj,ipompl) = adcc
            v(ibkv+nj,ipompl) = bdcc
            w(ibkw+nj,ipompl) = gdcc

         else if( sdir(j) .le. -1.5 .and. sdir(j) .gt. -2.5 ) then  ! -all

            ux    = -( x(ibkx+nj,ipompl) - sx0(j) ) / ra
            uy    = -( y(ibky+nj,ipompl) - sy0(j) ) / ra
            costh = -( z(ibkz+nj,ipompl) - sz0(j) ) / ra

            rt2   = ux**2 + uy**2

            if( rt2 .eq. 0.0d0 ) then

               sinth = 0.0d0
               cosphi= 1.0d0
               sinphi= 0.0d0

            else

               rt = sqrt(rt2)
               sinth  = rt
               cosphi = ux / rt
               sinphi = uy / rt

            end if

  105       cx = unirn(dummy)
            if( unirn(dummy) .gt. cx**2 ) goto 105

            wtnorm = wtnorm * 2./3. / cx

            sx = sqrt( 1.0 - cx**2 )
            ph = unirn(dummy) * 2.0 * pi

            adc = sx * cos( ph )
            bdc = sx * sin( ph )
            gdc = cx

            t1   = costh  * adc + sinth  * gdc
            adcc = cosphi * t1  - sinphi * bdc
            bdcc = sinphi * t1  + cosphi * bdc
            gdcc = costh  * gdc - sinth  * adc

            u(ibku+nj,ipompl) = adcc
            v(ibkv+nj,ipompl) = bdcc
            w(ibkw+nj,ipompl) = gdcc

         else  ! new option (speific angle) introduced in version 3.35

          if( sdir(j) .lt. 350.0 .and. sdir(j) .gt. 250 ) then ! adata is specified, T.Sato 2024/07/20

            jin = j

            call angbin(wdir,aw0,jin)

            wtnorm = wtnorm * aw0 ! change weight

          else
            wdir = sdir(j)
          endif

          if(wdir.lt.-1.0d0) wdir=-1.0d0
          if(wdir.gt.1.0d0) wdir=1.0d0

          if( wdir.gt.0.0 ) then ! go outside

            ux    = ( x(ibkx+nj,ipompl) - sx0(j) ) / ra
            uy    = ( y(ibky+nj,ipompl) - sy0(j) ) / ra
            costh = ( z(ibkz+nj,ipompl) - sz0(j) ) / ra
            rt2   = ux**2 + uy**2

            if( rt2 .eq. 0.0d0 ) then

               sinth = 0.0d0
               cosphi= 1.0d0
               sinphi= 0.0d0

            else

               rt = sqrt(rt2)
               sinth  = rt
               cosphi = ux / rt
               sinphi = uy / rt

            end if

            sx = sqrt( 1.0 - cx**2 )
            ph = unirn(dummy) * 2.0 * pi

            adc = sx * cos( ph )
            bdc = sx * sin( ph )
            gdc = cx

            t1   = costh  * adc + sinth  * gdc
            adcc = cosphi * t1  - sinphi * bdc
            bdcc = sinphi * t1  + cosphi * bdc
            gdcc = costh  * gdc - sinth  * adc

            u(ibku+nj,ipompl) = adcc
            v(ibkv+nj,ipompl) = bdcc
            w(ibkw+nj,ipompl) = gdcc

          else ! go inside

            ux    = -( x(ibkx+nj,ipompl) - sx0(j) ) / ra
            uy    = -( y(ibky+nj,ipompl) - sy0(j) ) / ra
            costh = -( z(ibkz+nj,ipompl) - sz0(j) ) / ra
            rt2   = ux**2 + uy**2

            if( rt2 .eq. 0.0d0 ) then

               sinth = 0.0d0
               cosphi= 1.0d0
               sinphi= 0.0d0

            else

               rt = sqrt(rt2)
               sinth  = rt
               cosphi = ux / rt
               sinphi = uy / rt

            end if

            sx = sqrt( 1.0 - cx**2 )
            ph = unirn(dummy) * 2.0 * pi

            adc = sx * cos( ph )
            bdc = sx * sin( ph )
            gdc = cx

            t1   = costh  * adc + sinth  * gdc
            adcc = cosphi * t1  - sinphi * bdc
            bdcc = sinphi * t1  + cosphi * bdc
            gdcc = costh  * gdc - sinth  * adc

            u(ibku+nj,ipompl) = adcc
            v(ibkv+nj,ipompl) = bdcc
            w(ibkw+nj,ipompl) = gdcc

          end if

         end if

       end if

*-----------------------------------------------------------------------

      else if( sdir(j) .lt. 1000.0 ) then
*-----------------------------------------------------------------------

            if( sdir(j) .lt. 100.0 ) then

               cdir = min(  1.0d0, sdir(j) )
               cdir = max( -1.0d0, cdir )

               wdir = cdir

            else if( sdir(j) .lt. 350.0 .and. sdir(j) .gt. 250 ) then

               jin = j

               call angbin(wdir,aw0,jin)

            else

               wdir = 2.0d0 * unirn(dummy) - 1.d+0

            end if

            if( sphi(j) .le. -1000.0 .or. (sdir(j) .gt. 100.0 .and. ! sphi=-1000 when phi=all (S.H. 2018.10.2)
     &       sdir(j) .lt. 250.0)) then ! 2016/10/27 Ogawa. Accept phi definition in case dir = data

               phi = 2.d+0 * unirn(dummy) * pi

            else

               phi = sphi(j) / 180.0 * pi

            end if

               costh  = wdir
               sinth  = sqrt( 1.d+0 - wdir**2 )
               cosphi = cos( phi )
               sinphi = sin( phi )

*-----------------------------------------------------------------------
*        sdom: normal
*-----------------------------------------------------------------------

            if( sdom(j) .gt. 0.0 ) then

               gdc = 1.0d0
     &             + ( cos( sdom(j) / 180.0 * pi ) - 1.0 )
     &             * unirn(dummy)

               ph = 2.d+0 * unirn(dummy) * pi
               dr = sqrt( 1.d+0 - gdc**2 )

               adc = dr * cos( ph )
               bdc = dr * sin( ph )

            else if( sdom(j) .gt. -1.5 .and. sdom(j) .lt. -0.5 ) then

  104          gdc = unirn(dummy)
               if( unirn(dummy) .gt. gdc ) goto 104

               ph = 2.d+0 * unirn(dummy) * pi
               dr = sqrt( 1.d+0 - gdc**2 )

               adc = dr * cos( ph )
               bdc = dr * sin( ph )

*-----------------------------------------------------------------------
*        sdom: duct source for jstyp = 1, 2, 4, 5
*-----------------------------------------------------------------------

            else if( sdom(j) .gt. -10.5 .and. sdom(j) .lt. -9.5 .and.
     &             ( jstyp(j) .eq. 1 .or. jstyp(j) .eq. 4 .or.
     &               jstyp(j) .eq. 2 .or. jstyp(j) .eq. 5 ) ) then

*-----------------------------------------------------------------------

                     fac1 = 1.d0

                     costh  = 1.0d0
                     sinth  = 0.0d0
                     cosphi = 1.0d0
                     sinphi = 0.0d0

                     xp0 = x(ibkx+nj,ipompl)
                     yp0 = y(ibky+nj,ipompl)
                     zp0 = z(ibkz+nj,ipompl)

               if( jstyp(j) .eq. 1 .or. jstyp(j) .eq. 4 ) then

                     thet = unirn(dummy) * 2.0 * pi
                     cs   = cos(thet)
                     sn   = sin(thet)

                     ard  = sdrd(j)

               else if( jstyp(j) .eq. 2 .or. jstyp(j) .eq. 5 ) then

                     xd0 = ( sx0(j) + sx1(j) ) / 2.0
                     yd0 = ( sy0(j) + sy1(j) ) / 2.0

                     ard = ( sdxw(j) * ( abs(  sdyw(j) / 2.0 - yp0 )
     &                              + abs( -sdyw(j) / 2.0 - yp0 ) )
     &                     + sdyw(j) * ( abs(  sdxw(j) / 2.0 - xp0 )
     &                              + abs( -sdxw(j) / 2.0 - xp0 ) ) )
     &                   / ( 2.0 * ( sdxw(j) + sdyw(j) ) )

               end if

               if( unirn(dummy) .lt. sdpf(j) ) then

                     fac0 = 1.d0 / sdpf(j) * sdl1(j)**2 / sdl2(j)**2

                  if( jstyp(j) .eq. 1 .or. jstyp(j) .eq. 4 ) then

                     r = unirn(dummy)
                     r = sdrd(j) * dmax1( r, unirn(dummy) )

                     xp1 = r * cs + sx0(j)
                     yp1 = r * sn + sy0(j)

                  else if( jstyp(j) .eq. 2 .or. jstyp(j) .eq. 5 ) then

                     xp1 = xd0 + sdxw(j) * unirn(dummy) - sdxw(j) / 2.d0
                     yp1 = yd0 + sdyw(j) * unirn(dummy) - sdyw(j) / 2.d0

                  end if

                     zp1 = sdl2(j) + z(ibkz+1,ipompl)

               else

                     fac0 = 1.d0 / ( 1.d0 - sdpf(j) )
     &                    * ( sdl2(j)**2 - sdl1(j)**2 ) / sdl2(j)**2

                  if( nglp(j) .le. 0 ) then

                     dlng = sdl1(j)
     &                    + ( sdl2(j) - sdl1(j) ) * unirn(dummy)
                     facl = 1.d0

                  else

                     random = unirn(dummy)
                     igs =  1

                     do ig = 1, nglp(j)

                        if( random. le. rgl(nglc(j)+ig) ) then

                           igs = ig
                           goto 3000

                        end if

                     end do

 3000                continue

                     dlng = slmin(ngli(j)+igs)
     &                    + ( slmax(ngla(j)+igs)
     &                      - slmin(ngli(j)+igs) )
     &                    * unirn(dummy)

                     facl = rgw(nglw(j)+igs)

                  end if

                     zp1 = dlng + z(ibkz+1,ipompl)

                  if( jstyp(j) .eq. 1 .or. jstyp(j) .eq. 4 ) then

                        xp1 = sdrd(j) * cs + sx0(j)
                        yp1 = sdrd(j) * sn + sy0(j)

                        dvp =  abs( sdrd(j) - xp0 * cs - yp0 * sn )

                  else if( jstyp(j) .eq. 2 .or. jstyp(j) .eq. 5 ) then

                        pram = unirn(dummy)

                     if( pram .le. psxp(j) ) then

                        xp1  = xd0 + sdxw(j) / 2.0
                        yp1  = yd0 + sdyw(j) * unirn(dummy)
     &                             - sdyw(j) / 2.0
                        fac1 = prs1(j)

                        dvp = abs( sdxw(j) / 2.0 - xp0 )

                     else if( pram .le. psyp(j) ) then

                        xp1  = xd0 + sdxw(j) * unirn(dummy)
     &                             - sdxw(j) / 2.0
                        yp1  = yd0 + sdyw(j) / 2.d0
                        fac1 = prs2(j)

                        dvp = abs( sdyw(j) / 2.0 - yp0 )

                     else if( pram .le. psxn(j) ) then

                        xp1  = xd0 - sdxw(j) / 2.0
                        yp1  = yd0 + sdyw(j) * unirn(dummy)
     &                             - sdyw(j) / 2.0
                        fac1 = prs3(j)

                        dvp = abs( sdxw(j) / 2.0 + xp0 )

                     else

                        xp1  = xd0 + sdxw(j) * unirn(dummy)
     &                             - sdxw(j) / 2.0
                        yp1  = yd0 - sdyw(j) / 2.d0
                        fac1 = prs4(j)

                        dvp = abs( sdyw(j) / 2.0 + yp0 )

                     end if

                  end if

                     fac0 = fac0 * facl / dlng**3
     &                    * 2.d0 * sdl1(j)**2 * sdl2(j)**2
     &                    / ( sdl1(j) + sdl2(j) )

                     fac0 = fac0 * dvp / ard

               end if

                  adc = xp1 - xp0
                  bdc = yp1 - yp0
                  gdc = zp1 - zp0

                  rnor = sqrt( adc**2 + bdc**2 + gdc**2 )

                  adc = adc / rnor
                  bdc = bdc / rnor
                  gdc = gdc / rnor

                  xd1 = adc / gdc * sdl0(j) + x(ibkx+1,ipompl) - xd0
                  yd1 = bdc / gdc * sdl0(j) + y(ibky+1,ipompl) - yd0

                  xd2 = adc / gdc * sdls(j) + x(ibkx+1,ipompl) - xd0
                  yd2 = bdc / gdc * sdls(j) + y(ibky+1,ipompl) - yd0

               if( jstyp(j) .eq. 1 .or. jstyp(j) .eq. 4 ) then

                  if( sqrt( xd1**2 + yd1**2 ) .gt. sdrd(j) ) goto 2000
                  if( sqrt( xp0**2 + yp0**2 ) .gt. sdrd(j) ) iduct = 0

               else if( jstyp(j) .eq. 2 .or. jstyp(j) .eq. 5 ) then

                     if( abs( xd1 ) .gt. sdxw(j) / 2.0 ) goto 2000
                     if( abs( yd1 ) .gt. sdyw(j) / 2.0 ) goto 2000

                     if( abs( xp0 ) .gt. sdxw(j) / 2.0 .or.
     &                   abs( yp0 ) .gt. sdyw(j) / 2.0 ) iduct = 0

               end if

               if( sdls(j) .gt. 0.0d0 ) then

                  if( sdrs(j) .gt. 0.0d0 ) then

                     if( sqrt( xd2**2 + yd2**2 ) .gt. sdrs(j) )
     &                  jduct = 0

                  else

                     if( abs( xd2 ) .gt. sdxs(j) / 2.0 .or.
     &                   abs( yd2 ) .gt. sdys(j) / 2.0 ) jduct = 0

                  end if

               end if

*-----------------------------------------------------------------------

               if( sdebg(j) .eq. 0.0d0 ) then

                  wtduct = fac0 * fac1

               end if

*-----------------------------------------------------------------------
*        sdom(j) = 0.0
*-----------------------------------------------------------------------

            else

               adc = 0.0
               bdc = 0.0
               gdc = 1.0

            end if

*-----------------------------------------------------------------------

               t1   = costh  * adc + sinth  * gdc
               adcc = cosphi * t1  - sinphi * bdc
               bdcc = sinphi * t1  + cosphi * bdc
               gdcc = costh  * gdc - sinth  * adc

               u(ibku+nj,ipompl) = adcc
               v(ibkv+nj,ipompl) = bdcc
               w(ibkw+nj,ipompl) = gdcc

               uabs = sqrt( u(ibku+nj,ipompl)**2
     &                    + v(ibkv+nj,ipompl)**2
     &                    + w(ibkw+nj,ipompl)**2 )

               u(ibku+nj,ipompl) = u(ibku+nj,ipompl) / uabs
               v(ibkv+nj,ipompl) = v(ibkv+nj,ipompl) / uabs
               w(ibkw+nj,ipompl) = w(ibkw+nj,ipompl) / uabs

*-----------------------------------------------------------------------

         if( jstyp(j) .eq. 26 ) then

               adc = u(ibku+nj,ipompl)
               bdc = v(ibkv+nj,ipompl)
               gdc = w(ibkw+nj,ipompl)

               t1   = scosth  * adc + ssinth  * gdc
               adcc = scosphi * t1  - ssinphi * bdc
               bdcc = ssinphi * t1  + scosphi * bdc
               gdcc = scosth  * gdc - ssinth  * adc

               u(ibku+nj,ipompl) = adcc
               v(ibkv+nj,ipompl) = bdcc
               w(ibkw+nj,ipompl) = gdcc

         end if

*-----------------------------------------------------------------------

      end if

      if(jstyp(j).eq.17) goto 5900 ! T.Sato 2021/12/07

c  T.Sato for Special correlation source (iscorr=2 and 3), 2013/11/15
            if(iscorr.eq.2) then ! all particles are emitted from the same location
             if(nj.ge.2) then ! all particles are generated at the same location
              x(ibkx+nj,ipompl) = x(ibkx+1,ipompl)
              y(ibky+nj,ipompl) = y(ibky+1,ipompl)
              z(ibkz+nj,ipompl) = z(ibkz+1,ipompl)
             endif
            endif
            if(iscorr.eq.3) then ! 2 particles are going to opposite direction
             if(nj.ge.2) then
              x(ibkx+nj,ipompl) = x(ibkx+1,ipompl)
              y(ibky+nj,ipompl) = y(ibky+1,ipompl)
              z(ibkz+nj,ipompl) = z(ibkz+1,ipompl)
              u(ibku+nj,ipompl) = -u(ibku+1,ipompl)
              v(ibkv+nj,ipompl) = -v(ibkv+1,ipompl)
              w(ibkw+nj,ipompl) = -w(ibkw+1,ipompl)
             endif
            endif
c *** End of special correlation source ********

*-----------------------------------------------------------------------
*     weight factor
*-----------------------------------------------------------------------

               wtnorm = wtnorm * aw0 * pmlpr(j)

*-----------------------------------------------------------------------
*     e(ibke+1)
*-----------------------------------------------------------------------

         if(iengbin.ne.1) ew0 = 1.0d0 ! T.Sato 2020/12/07 (iengbin=1 indicates eo & ew0 have has been already determined)

         if( ispfs(j) .ne. 0 ) then

         else
     &   if( jstyp(j) .eq.  4 .or.
     &       jstyp(j) .eq.  5 .or.
     &       jstyp(j) .eq.  6 .or.
     &       jstyp(j) .eq.  8 .or.
     &       jstyp(j) .eq. 10 .or.
     &       jstyp(j) .eq.-11 .or.  ! T.Sato 2020/11/29
     &       jstyp(j) .eq. 14 .or.
     &       jstyp(j) .eq. 16 .or.
     &       jstyp(j) .eq. 19 .or.
     &       jstyp(j) .eq. 21 .or.
     &       jstyp(j) .eq. 23 .or.
     &      ( jstyp(j) .eq. 26 .and. jetyp(j) .gt. 0 ) .or.
     &       jstyp(j) .eq. 25 .or.
     &     ( jstyp(j) .eq. 100 .and. jetyp(j) .gt. 0 ) ) then

               jin = j
               if(iengbin.ne.1) then ! T.Sato 2020/12/07 (iengbin=1 indicates eo & ew0 have has been already determined)
                if(icenv(j).le.0) then ! conventional mode
                 call engbin(eo,ew0,jin,idummy)
                 if(sdir(j).lt.-3.5) then ! plane cosmic-ray mode, cos(theta) biased
 1003             wtmp = 1.0d0 - 2.0d0*unirn(dummy)
                  if(unirn(dummy).gt.abs(wtmp)) goto 1003
                  if(wtmp.lt.sag1(j).or.wtmp.gt.sag2(j)) izeroe=1  ! zero-energy index
                  if(izeroe.eq.1.and.isbias(j).eq.1) then ! do not allow to produce dead particle, T.Sato 2021/03/03
                   izeroe=0
                   goto 1003
                  endif
                 endif
                else  ! terrestrial cosmic-ray mode, check angular limitation
                 call engbin(eo,ew0,jin,ie)
 1002            wtmp = -getCosmicAng(jin,ie)  ! w = -cx = -getCosmicAng(jin,ie)
                 if(sdir(j).lt.-3.5) then ! plane cosmic-ray mode, cos(theta) biased
                  if(unirn(dummy).gt.abs(wtmp)) izeroe=1 ! not necessary to be re-sampling, T.Sato 2024/01/18
                 endif
                 if(wtmp.lt.sag1(j).or.wtmp.gt.sag2(j)) izeroe=1  ! zero-energy index
                 if(izeroe.eq.1.and.isbias(j).eq.1) then ! do not allow to produce dead particle, T.Sato 2021/03/03
                  izeroe=0
                  goto 1002
                 endif
                 if(ipcosmic(j).eq.33.and.ie.eq.ie511(j)) then
                  if(unirn(dummy).gt.1.0d0/annihratio(j))
     &            eo=rstms(12)*1.0d3
                 endif
                endif
                if(sdir(j).lt.-3.5) then ! plane cosmic-ray mode, use wtmp for determing direction
                 cx = wtmp
                 sx = sqrt( 1.0 - cx**2 )
                 ph = unirn(dummy) * 2.0 * pi
                 u(ibku+nj,ipompl)  = sx * cos( ph )
                 v(ibkv+nj,ipompl)  = sx * sin( ph )
                 w(ibkw+nj,ipompl)  = cx
                endif
            endif


               e(ibke+nj,ipompl) = eo
               wtnorm = wtnorm * ew0

         else if( jstyp(j) .ne. 12 ) then !FURUTA20210726

               e(ibke+nj,ipompl) = se0(j)

         end if


*-----------------------------------------------------------------------
*     nty(ibknty+1), nkf(ibknkf+1), name(ibknam+1),
*-----------------------------------------------------------------------

      if( ispfs(j) .eq. 0 ) then

         nty(ibknty+nj,ipompl)  = istyp(j)
         nkf(ibknkf+nj,ipompl)  = inkf0(j)
         name(ibknam+nj,ipompl) = name0

*-----------------------------------------------------------------------
*     Spontaneus Fission Source
*-----------------------------------------------------------------------

      else

         nty(ibknty+nj,ipompl)  = 2
         nkf(ibknkf+nj,ipompl)  = 2112
         name(ibknam+nj,ipompl) = name0

      end if

*-----------------------------------------------------------------------
*     for nucleus, input energy is MeV / nuleon
*-----------------------------------------------------------------------

      if( nty(ibknty+nj,ipompl) .ge. 15 ) then

         e(ibke+nj,ipompl) = e(ibke+nj,ipompl) *
     &                       ibryf(nty(ibknty+nj,ipompl),
     &                       nkf(ibknkf+nj,ipompl))

      end if

      if(izeroe.eq.1) then
       e(ibke+nj,ipompl)=0.0  ! T.Sato 2020/12/14, energy is forced to be 0
       nty(ibknty+nj,ipompl)=14 ! forced to be photon because some particle decays
       nkf(ibknkf+nj,ipompl)=22 ! forced to be photon because some particle decays
      endif

*-----------------------------------------------------------------------
*     t(ibkt+1)
*-----------------------------------------------------------------------


               tw0 = 1.0d0

         if( jttyp(j) .eq.   1 .or.
     &       jttyp(j) .eq.   2 .or.
     &       jttyp(j) .eq.   3 .or.
     &       jttyp(j) .eq.   4 .or.
     &       jttyp(j) .eq.   5 .or.
     &       jttyp(j) .eq.   6 .or.
     &       jttyp(j) .eq. 100 .or.
     &     ( jstyp(j) .eq. 100 .and. jttyp(j) .gt. 0 ) ) then

               eo = e(ibke+nj,ipompl)

               jin = j

               call timbin(to,tw0,eo,jin)

               t(ibkt+nj,ipompl) = to
               wtnorm = wtnorm * tw0

         end if

*-----------------------------------------------------------------------
*     initial flag
*-----------------------------------------------------------------------

      if( ksou1(j) .eq. 0 ) ksou1(j) = 1

*-----------------------------------------------------------------------
*     transformation
*-----------------------------------------------------------------------

         if( isort(j,4) .ne. 0 .and. (iscorr .le. 1 .or.    ! 2018/09/04 Ogawa. Apply if source positions are uncorrelated
     &   (iscorr .ge. 2 .and. j .eq. 1) )) then             ! do not apply second and subsequent if source positions is that of the first one

            call trnsxv(x(ibkx+nj,ipompl),y(ibky+nj,ipompl),
     &                  z(ibkz+nj,ipompl),
     &                  ux,uy,uz,isort(j,4))

            call trnsuv(u(ibku+nj,ipompl),v(ibkv+nj,ipompl),
     &                  w(ibkw+nj,ipompl),
     &                  uu,uv,uw,isort(j,4))

               x(ibkx+nj,ipompl) = ux
               y(ibky+nj,ipompl) = uy
               z(ibkz+nj,ipompl) = uz

               u(ibku+nj,ipompl) = uu
               v(ibkv+nj,ipompl) = uv
               w(ibkw+nj,ipompl) = uw

         end if

*-----------------------------------------------------------------------
*     find initial region
*-----------------------------------------------------------------------

      if( icntl .eq. 0 .or. icntl .eq. 5 .or.
     &    icntl .eq. 6 .or.
     &    icntl .eq. 14 .or.
     &    icntl .eq. 15 ) then

               isrr  = 0

  700    continue

               ici   = 0
               mark  = 1
               markp = 0

            call gomsor(x(ibkx+nj,ipompl),y(ibky+nj,ipompl),
     &                  z(ibkz+nj,ipompl),
     &                  u(ibku+nj,ipompl),v(ibkv+nj,ipompl),
     &                  w(ibkw+nj,ipompl),
     &                  nmed(ibknmd+nj,ipompl),iblz(ibkblz+nj,ipompl),
     &                  mark,markp,ici)

*-----------------------------------------------------------------------
*           check the point on boundary
*-----------------------------------------------------------------------

               markp = 1

            call gomdis(dpr,
     &                  x(ibkx+nj,ipompl),y(ibky+nj,ipompl),
     &                  z(ibkz+nj,ipompl),
     &                  u(ibku+nj,ipompl),v(ibkv+nj,ipompl),
     &                  w(ibkw+nj,ipompl),
     &                  mark,markp,
     &                  nmed(ibknmd+nj,ipompl),iblz(ibkblz+nj,ipompl))

               markp = 0

            if( mark .eq. -6 ) then

               if( isrr .eq. 0 ) then

                  x(ibkx+nj,ipompl) = x(ibkx+nj,ipompl) +
     &                                     u(ibku+nj,ipompl) * parz(28)
                  y(ibky+nj,ipompl) = y(ibky+nj,ipompl) +
     &                                     v(ibkv+nj,ipompl) * parz(28)
                  z(ibkz+nj,ipompl) = z(ibkz+nj,ipompl) +
     &                                     w(ibkw+nj,ipompl) * parz(28)

                  isrr = isrr + 1
                  goto 700

               else

                  ErrCha = ''
                  ErrID = 'L:2376/R:sors/F:sors.f' !E03_106_001
                  call ErrWrite(ErrID,ErrCha)

        if( npe .gt. 0 ) then ! T.Sato 2019/01/07, for MPI
         write(*,'( '' my ip       = '',i3)') me
        endif

        write(*,'(/
     &  ''** lost particle in sors, source is not in any cell'')')
        write(*,*)  x(ibkx+nj,ipompl), y(ibky+nj,ipompl),
     &              z(ibkz+nj,ipompl)

            ilost = ilost + 1

            if( ilost .ge. nlost ) then

             write(ErrCha,'(''Calculation is terminated!!!'')')
             ErrID = 'L:2393/R:sors/F:sors.f'
             call ErrWrite(ErrID,ErrCha)
             write(ErrCha,'(''number of'',
     &       '' lost particles exceeds nlost (='',i8,'')'')') nlost
             ErrID = 'L:2397/R:sors/F:sors.f'
             call ErrWrite(ErrID,ErrCha)
             call parastop( 602 )

            else ! resampling is required

             iresample=1
             return

            endif


               end if

            end if

*-----------------------------------------------------------------------
*        region selection ( maxmum trial is ntmax )
*-----------------------------------------------------------------------

            if( nsrn(j) .gt. 0 .and. jstyp(j) .ne. 12
     &           .and. jstyp(j) .ne. 24 !FURUTA20180115
     &           .and. jstyp(j) .ne. 25 !FURUTA20180115
     &           .and. icntl .ne. 14    !T.Sato 2025/03/05, reg should be ignored in volume calculation
     &           ) then

                  idsm = iaddress_nsrn(j)
                  jdsm = 0

                  jj  = 0

                  jdsm = jdsm + 1
                  ntrn = idas_nsrn(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_nsrn(idsm+jdsm)

               do ii = 1, ntrn

                  call tregck(iblz1,ilev1,ilat1,
     &                        mtrn,idas_nsrn(idsm+jdsm+1),jj,icc)

                  if( icc .ne. 0 ) goto 500

               end do

               itsrg = itsrg + 1

            if( itsrg .gt. nsmx(j) ) then
              ErrCha = ''
              ErrID = 'L:2446/R:sors/F:sors.f' !E03_107_001
              call ErrWrite(ErrID,ErrCha)

               write(io,'(/'' Error : region selection failed,'',
     &                     '' in source.''/
     &                     '' max. trial ntmax = '',i5)') nsmx(j)

               write(*,'(/'' Error : region selection failed,'',
     &                     '' in source.''/
     &                     '' max. trial ntmax = '',i5)') nsmx(j)


               call parastop( 889 )

            end if

            goto 2000

  500       continue

         end if

      end if

*-----------------------------------------------------------------------
*     phase of this event
*-----------------------------------------------------------------------

            pmphs = 2.0d0 * pi * unirn(dummy)

*-----------------------------------------------------------------------
*     sorce number j
*-----------------------------------------------------------------------

            nsos(ibksos+nj,ipompl) = j

*-----------------------------------------------------------------------
*     charge state
*-----------------------------------------------------------------------

         if( lstyp(j) .gt. -1000.d0 ) then

            nzst(ibkzst+nj,ipompl) = lstyp(j)

         else

            nzst(ibkzst+nj,ipompl) = ichgf(istyp(j),inkf0(j))

         end if

*-----------------------------------------------------------------------
*     initial spin,  ssr = 0 ; no defined
*-----------------------------------------------------------------------

             ssr = ssx(j)**2 + ssy(j)**2 + ssz(j)**2

         if( ssr .gt. 1.d-8 ) then

            spx(ibkspx+nj,ipompl) = ssx(j) / ssr
            spy(ibkspy+nj,ipompl) = ssy(j) / ssr
            spz(ibkspz+nj,ipompl) = ssz(j) / ssr

         else

            spx(ibkspx+nj,ipompl) = 0.d0
            spy(ibkspy+nj,ipompl) = 0.d0
            spz(ibkspz+nj,ipompl) = 0.d0

         end if

*-----------------------------------------------------------------------
*        wt(ibkwt+nj)
*-----------------------------------------------------------------------

            wt(ibkwt+nj,ipompl) = wtnorm * wtduct

*-----------------------------------------------------------------------
*     Spontaneus Fission Source
*-----------------------------------------------------------------------

            wtspfs = 1.0d0

         if( ispfs(j) .ne. 0 ) then

            ichag = abs(inkf0(j)) / 1000000
            imass = abs(inkf0(j)) - ichag * 1000000
            izaid = ichag * 1000 + imass

            nCf252option = 0

            nneut = nSmpSpNuDistData(izaid, nCf252option)

*-----------------------------------------------------------------------
*           if neutron is zero, add source weight for ispfs=1
*-----------------------------------------------------------------------

            if( nneut .eq. 0 ) then

               if( ispfs(j) .eq. 1 ) then

C for NONCRITICAL_OVLY12
                  rsouin_local = rsouin_local
     &                         + wtnorm * wtspfs * wtduct
     &                           / sfactor(j) / abs( totfact )

               end if
!$OMP CRITICAL (rspfz_crit)
                  rspfz  = rspfz + wtnorm * wtspfs * wtduct
     &                   / sfactor(j) / abs( totfact )
!$OMP END CRITICAL (rspfz_crit)

                  goto 4000

            else if( nneut .lt. 0 ) then

               write(io,'(/'' Error : proj is not listed in the'',
     &                     '' SPFS table. ZAID = '',i5)') izaid
               call parastop( 789 )

            else
!$OMP ATOMIC
               rspfn = rspfn + nneut

            end if

*-----------------------------------------------------------------------
*              add number of neutron on source weight for ispfs=2
*-----------------------------------------------------------------------

               if( ispfs(j) .eq. 2 ) then

                  wtspfs = nneut

               end if

*-----------------------------------------------------------------------
*           energy of spontaneus fission neutron
*           and store SP fission neutrons in bank for nneut > 1
*-----------------------------------------------------------------------

                  ePart = 0.0d0

                  ewatt = SmpWatt(ePart, izaid)
                  e(ibke+nj,ipompl) = ewatt

*-----------------------------------------------------------------------

            if( nneut .gt. 1 ) then

                  nj0 = nj + nneut - 1

               do kk = 2, nneut

                   k = kk + nj - 1

                  ewatt = SmpWatt(ePart, izaid)
                  e(ibke+k,ipompl) = ewatt

                  t(ibkt+k,ipompl)      = t(ibkt+nj,ipompl)
                  name(ibknam+k,ipompl) = name(ibknam+nj,ipompl)
                  nty(ibknty+k,ipompl)  = nty(ibknty+nj,ipompl)
                  nkf(ibknkf+k,ipompl)  = nkf(ibknkf+nj,ipompl)

                  x(ibkx+k,ipompl) = x(ibkx+nj,ipompl)
                  y(ibky+k,ipompl) = y(ibky+nj,ipompl)
                  z(ibkz+k,ipompl) = z(ibkz+nj,ipompl)

                  cx = 1.0 - 2.0 * unirn(dummy)
                  sx = sqrt( 1.0 - cx**2 )
                  ph = unirn(dummy) * 2.0 * pi

                  u(ibku+k,ipompl)  = sx * cos( ph )
                  v(ibkv+k,ipompl)  = sx * sin( ph )
                  w(ibkw+k,ipompl)  = cx

                  wt(ibkwt+k,ipompl)      = wt(ibkwt+nj,ipompl)
                  iblz(ibkblz+k,ipompl)   = iblz(ibkblz+nj,ipompl)
                  nmed(ibknmd+k,ipompl)   = nmed(ibknmd+nj,ipompl)

                  spx(ibkspx+k,ipompl)    = spx(ibkspx+nj,ipompl)
                  spy(ibkspy+k,ipompl)    = spy(ibkspy+nj,ipompl)
                  spz(ibkspz+k,ipompl)    = spz(ibkspz+nj,ipompl)

                  nzst(ibkzst+k,ipompl)    = nzst(ibkzst+nj,ipompl)
                  nsos(ibksos+k,ipompl)    = nsos(ibksos+nj,ipompl)

                  ncnt(ibknct+1,k,ipompl) = ncnt(ibknct+1,nj,ipompl)
                  ncnt(ibknct+2,k,ipompl) = ncnt(ibknct+2,nj,ipompl)
                  ncnt(ibknct+3,k,ipompl) = ncnt(ibknct+3,nj,ipompl)
                  itetpos(ibtetpos+k,ipompl)=itetpos(ibtetpos+nj,ipompl)

               end do

                  nj = nj0

            end if

         end if

 5800    continue

         if( ksou1node(j) .eq. 0) ksou1node(j) = 1

 5900    continue

*-----------------------------------------------------------------------
*     rsouin : normalization of source with sfactor and totfact
*              dnorm for special duct
*-----------------------------------------------------------------------

               kduct = 0
        if(jstyp(j) .eq. 17 .and. idmpmode .eq. 1)then

C for NONCRITICAL_OVLY12
         rsouin_local = rsouin_local + 1.0d0
        else
c------------------
         if( sdom(j) .gt. -10.5 .and. sdom(j) .lt. -9.5 .and.
     &       sdebg(j) .eq. 0.0d0 ) then

            if( iduct .eq. 1 ) then
C for NONCRITICAL_OVLY12
               rsouin_local  = rsouin_local
     &                       + wtnorm * wtspfs / dnorm(j)
     &                         / sfactor(j) / abs( totfact )
            else

               kduct = j

            end if

         else
               rsouin_local  = rsouin_local
     &                       + wtnorm * wtspfs * wtduct
     &                         / sfactor(j) / abs( totfact )

         end if !FURUTA20150515

        end if

*-----------------------------------------------------------------------
*        special for duct source ( slit rejection )
*-----------------------------------------------------------------------

         if( jduct .eq. 0 ) goto 4000

*-----------------------------------------------------------------------

         if( iscorr .ne. 0 ) goto 7000

*-----------------------------------------------------------------------
         if(nmulti.eq.0)then
          nomax=0
          no=1
          return
         endif

*-----------------------------------------------------------------------
*     start initial particle
*-----------------------------------------------------------------------

 6000 continue

            if( jstyp(j) .ne. 100 ) nomax = nj

*-----------------------------------------------------------------------

            ec(ibkec+1,ipompl) = e(ibke+1,ipompl)
            tc(ibktc+1,ipompl) = t(ibkt+1,ipompl)
            xc(ibkxc+1,ipompl) = x(ibkx+1,ipompl)
            yc(ibkyc+1,ipompl) = y(ibky+1,ipompl)
            zc(ibkzc+1,ipompl) = z(ibkz+1,ipompl)

            ityp = nty(ibknty+1,ipompl)
            ktyp = nkf(ibknkf+1,ipompl)
            jtyp = ichgf(ityp,ktyp,ipompl)
            mtyp = ibryf(ityp,ktyp,ipompl)
            rtyp = rmtyp(ityp,ktyp,ipompl)

            oldwt = wt(ibkwt+1,ipompl)
            mat   = nmed(ibknmd+1,ipompl)

            if( ityp .eq. 2 .and.
     &          e(ibke+1,ipompl) .le. max(emin(2),dnmax(2)) ) then
!$OMP ATOMIC
                 reutn = reutn + wtspfs

            end if

*-----------------------------------------------------------------------
*        initial weight
*-----------------------------------------------------------------------

         do k = 1, 19

            wc01(k) = wc1(k) * wt(ibkwt+1,ipompl)
            wc02(k) = wc2(k) * wt(ibkwt+1,ipompl)

         end do

*-----------------------------------------------------------------------
*        keep initial cell importance: forced collisions flag
*-----------------------------------------------------------------------

      if( icntl .eq. 0 .or. icntl .eq. 5 .or.
     &    icntl .eq. 6 .or.
     &    icntl .eq. 14 .or.
     &    icntl .eq. 15 ) then

         do k = 1, nomax

            wtin(ibkwin+k,ipompl) =
     &                  aimp(nty(ibknty+k,ipompl),iblz1,ilev1,ilat1,ii1)
            wtnz(ibkwnz+k,ipompl) =
     &                  aimp(nty(ibknty+k,ipompl),iblz1,ilev1,ilat1,ii1)

            nfcs(ibknfc+k,ipompl) = 0

            if( wtin(ibkwin+k,ipompl) .le. 0.0 ) then

               write(io,'(/'' Error : importance is zero at'',
     &                     '' initial position in source.''/
     &                     '' region = '',i5)') iblz1

               call parastop( 888 )

            end if

         end do

      end if

*-----------------------------------------------------------------------
*           mat time change
*-----------------------------------------------------------------------

            if( mttcn .gt. 0 .and. mat .gt. 0 ) then

               do k = 1, mttcn

                  if( mttc1(k) .eq. idmn(mat) .and.
     &                smttc(k) .lt. t(ibkt+1,ipompl) ) then

                     if( mttc2(k) .gt. 0 ) then

                        nmed(ibknmd+1,ipompl) = idnm(mttc2(k))

                     else

                        nmed(ibknmd+1,ipompl) = mttc2(k)

                     end if

                        mat = nmed(ibknmd+1,ipompl)

                  end if

               end do

            end if

*-----------------------------------------------------------------------
*           write source information on output
*-----------------------------------------------------------------------

            if( rcasc .le. dble(mstz(55)) ) then

               if( rcasc .lt. 2.0d0 ) then

                  write(io,'(/''***** source information upto'',i4,
     &            '' histories *****'')') mstz(55)

                  write(io,'(''  no ityp    reg  mat'',
     &            ''    energy  '',
     &            ''    x       '',
     &            ''    y       '',
     &            ''    z       '',
     &            ''    u       '',
     &            ''    v       '',
     &            ''    w       '',
     &            ''    wt      '',
     &            ''  time (ns)'')')  ! MATSUDA 2023.05.10

               end if

                  do k = 1, nomax

                  write(io,'(i4,i5,i7,i5,9(1p1e12.4))')
     &               nint(rcasc),
     &               nty(ibknty+k,ipompl),
     &               iblz(ibkblz+k,ipompl),
     &               idmn(nmed(ibknmd+k,ipompl)),
     &               e(ibke+k,ipompl),
     &               x(ibkx+k,ipompl),
     &               y(ibky+k,ipompl),
     &               z(ibkz+k,ipompl),
     &               u(ibku+k,ipompl),
     &               v(ibkv+k,ipompl),
     &               w(ibkw+k,ipompl),
     &               wt(ibkwt+k,ipompl),
     &               t(ibkwt+k,ipompl)  ! MATSUDA 2023.05.10

                  end do

            end if

*-----------------------------------------------------------------------
*     save initial energy for adjoint mode
*-----------------------------------------------------------------------
               engini = e(ibke+1,ipompl)

*-----------------------------------------------------------------------
C for NONCRITICAL_OVLY12
!$OMP ATOMIC
        rsouin = rsouin + rsouin_local
C for STRICT_RSOUIN
!$      rsouin_arr(nocas) = rsouin_local
      return
      end


************************************************************************
*                                                                      *
      subroutine engbin(eo,wo,j,igs)
*                                                                      *
*        purpose : determine the incident energy from e-group          *
*        last modified by T.Sato on 2020/12/15                         *
*                                                                      *
************************************************************************
      use moddas_source

      implicit real*8 (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      common /isorsg/ seg0(isrc),seg1(isrc),seg2(isrc),seg3(isrc),
     &                set0(isrc),set1(isrc),set2(isrc),set3(isrc),
     &                jetyp(isrc),jptyp(isrc)

      common /isorse/ ngrp(isrc), ngei(isrc), ngea(isrc), ngfe(isrc),
     &                ngft(isrc), ngll(isrc), ngpi(isrc), ngpw(isrc)

*-----------------------------------------------------------------------
*     given by the data
*-----------------------------------------------------------------------

      if( jetyp(j) .eq.  1 .or. jetyp(j) .eq. 11 .or.
     &    jetyp(j) .eq.  4 .or. jetyp(j) .eq. 14 .or.
     &    jetyp(j) .eq. 21 .or. jetyp(j) .eq. 31 .or.
     &    jetyp(j) .eq. 24 .or. jetyp(j) .eq. 34 .or.
     &    jetyp(j) .eq.  8 .or. jetyp(j) .eq. 18 .or.
     &    jetyp(j) .eq.  9 .or. jetyp(j) .eq. 19 .or.
     &    jetyp(j) .eq. 20 .or.                      ! S.H. (2016.12.30)
     &    jetyp(j) .eq. 28 .or. jetyp(j) .eq. 29 .or.
     &    jetyp(j) .eq. 22 .or. jetyp(j) .eq. 32 .or.
     &    jetyp(j) .eq. 23 .or. jetyp(j) .eq. 33 .or.
     &    jetyp(j) .eq. 25 .or. jetyp(j) .eq. 26 .or. ! T.Sato 2020/12/06 cosmic-ray source
c
     &    jetyp(j) .eq.  5 .or. jetyp(j) .eq. 15 .or.
     &    jetyp(j) .eq.  6 .or. jetyp(j) .eq. 16 .or.
     &    jetyp(j) .eq.  3 .or. jetyp(j) .eq.  7 ) then

*-----------------------------------------------------------------------

            random = unirn(dummy)
            igs =  1

         do ig = 1, ngrp(j)

            if( random. le. rfe(ngft(j)+ig) ) then

               igs = ig
               goto 1000

            end if

         end do

 1000    continue

            wo = pwt(ngpw(j)+igs)

*-----------------------------------------------------------------------

         if( egmax(ngea(j)+igs) .eq. egmin(ngei(j)+igs) ) then

               eo =   egmin(ngei(j)+igs)

         else if( ngll(j) .gt. 0 ) then

               eo =   egmin(ngei(j)+igs)
     &            + ( egmax(ngea(j)+igs)
     &            -   egmin(ngei(j)+igs) ) * unirn(dummy)

         else

               eo = egmin(ngei(j)+igs)
     &            * exp( log( egmax(ngea(j)+igs) / egmin(ngei(j)+igs) )
     &            * unirn(dummy) )

         end if

*-----------------------------------------------------------------------
*     gaussian distribution
*-----------------------------------------------------------------------

      else if( jetyp(j) .eq. 2 .or. jetyp(j) .eq. 12 ) then

*-----------------------------------------------------------------------

            xeg = seg1(j) / 2.35482d+0

  100       eo = xeg * gaurn(dummy) + seg0(j)

            if( eo .lt. seg2(j) .or. eo .gt. seg3(j) ) goto 100

            wo = 1.0d0

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*     wave length to energy
*-----------------------------------------------------------------------

      if( jetyp(j) .eq. 11 .or. jetyp(j) .eq. 12 .or.
     &    jetyp(j) .eq. 15 .or. jetyp(j) .eq. 16 .or.
     &    jetyp(j) .eq. 18 .or. jetyp(j) .eq. 19 .or.
     &    jetyp(j) .eq. 32 .or. jetyp(j) .eq. 33 .or.
c
     &    jetyp(j) .eq. 31 .or. jetyp(j) .eq. 34 ) then

         eo = 8.1804250d-8 / eo**2

      end if

*-----------------------------------------------------------------------

      return
      end

*-----------------------------------------------------------------------

************************************************************************
*                                                                      *
      subroutine angbin(ao,wo,j)
*                                                                      *
*        purpose : determine the angle from a-group                    *
*        last modified by K.Niita on 2005/11/08                        *
*                                                                      *
************************************************************************
      use moddas_source

      implicit real*8 (a-h,o-z)

      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /isoraa/ narp(isrc), naei(isrc), naea(isrc), nafe(isrc),
     &                naft(isrc), nall(isrc), napi(isrc), napw(isrc)


      common /isoras/ sag1(isrc),sag2(isrc),jatyp(isrc),jqtyp(isrc)

*-----------------------------------------------------------------------

      common /isowtt/ iggos(isrc)

            iggos(j) = 1

*-----------------------------------------------------------------------
*     given by the data
*-----------------------------------------------------------------------

      if( jatyp(j) .eq.  1 .or. jatyp(j) .eq. 11 .or.
     &    jatyp(j) .eq.  4 .or. jatyp(j) .eq. 14 .or.
     &    jatyp(j) .eq.  5 .or. jatyp(j) .eq. 15 .or.
     &    jatyp(j) .eq.  6 .or. jatyp(j) .eq. 16 ) then

*-----------------------------------------------------------------------

            random = unirn(dummy)
            igs =  1

         do ig = 1, narp(j)

            if( random. le. rfa(naft(j)+ig) ) then

               igs = ig
               goto 1000

            end if

         end do

 1000    continue

            iggos(j) = igs

            wo = pat(napw(j)+igs)

         if( jatyp(j) .lt. 10 ) then

            ao =   agmin(naei(j)+igs)
     &         + ( agmax(naea(j)+igs)
     &         -   agmin(naei(j)+igs) ) * unirn(dummy)

         else

            gamin = cos( agmin(naei(j)+igs) / 180.0d0 * pi )
            gamax = cos( agmax(naea(j)+igs) / 180.0d0 * pi )

            ao = gamin + ( gamax - gamin ) * unirn(dummy)

         end if

*-----------------------------------------------------------------------

      end if

      if( ao .lt. -1.d0 .or. ao .gt. 1.d0 )  call parastop( 889 )

*-----------------------------------------------------------------------

      return
      end

*-----------------------------------------------------------------------

************************************************************************
*                                                                      *
      subroutine timbin(to,wo,eo,j)
*                                                                      *
*        purpose : determine the time from t-group                     *
*        last modified by K.Niita on 2015/02/09                        *
*                                                                      *
************************************************************************
      use moddas_source

      implicit real*8 (a-h,o-z)

      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /isorsm/ stm0(isrc),stmw(isrc),stmc(isrc),stmd(isrc),
     &                jttyp(isrc),jttpn(isrc)

      common /isortt/ ntrp(isrc), ntei(isrc), ntea(isrc), ntfe(isrc),
     &                ntft(isrc), ntll(isrc), ntpi(isrc), ntpw(isrc)

      common /isorft/ isll(isrc), lsfz(isrc), srfz(isrc)
      character srfz*200


      common /isorts/ stg1(isrc),stg2(isrc),jotyp(isrc)

*-----------------------------------------------------------------------
*     given by the data
*-----------------------------------------------------------------------

      if( jttyp(j) .eq. 1 .or. jttyp(j) .eq. 2 ) then

               npch = int( unirn(dummy) * jttpn(j) )
               if( npch .eq. jttpn(j) ) npch = npch - 1

               pkcn = stm0(j) + npch * stmd(j)

            if( jttyp(j) .eq. 1 ) then

               to = pkcn + ( 2.0 * unirn(dummy) - 1.0 ) * stmw(j) / 2.0

            else if( jttyp(j) .eq. 2 ) then

  133          widtt = stmw(j) / 2.35482d+0 * gaurn(dummy)
               if( abs(widtt) .gt. stmc(j) ) goto 133

               to = pkcn + widtt

            end if

               wo = 1.d0

*-----------------------------------------------------------------------

      else if( jttyp(j) .eq.  3 .or. jttyp(j) .eq. 4 .or.
     &         jttyp(j) .eq.  5 .or. jttyp(j) .eq. 6  ) then

*-----------------------------------------------------------------------

            random = unirn(dummy)
            igs =  1

         do ig = 1, ntrp(j)

            if( random. le. rft(ntft(j)+ig) ) then

               igs = ig
               goto 200

            end if

         end do

  200    continue

            wo = ptt(ntpw(j)+igs)

            to =   tgmin(ntei(j)+igs)
     &         + ( tgmax(ntea(j)+igs)
     &         -   tgmin(ntei(j)+igs) ) * unirn(dummy)

*-----------------------------------------------------------------------

      else if( jttyp(j) .eq.  100 ) then

            t1 = stg1(j)
            t2 = stg2(j)
            e0 = eo
            t0 = 0.d0
            ll = isll(j)

         call tdis01(t0,e0,t1,t2,ll)

            to = t0

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine tdis01(t0,e0,t1,t2,ll)
*                                                                      *
*        purpose : determine the time from the special parametrization *
*        last modified by K.Niita on 2015/02/10                        *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      data ifirst /0/
      save ifirst

*-----------------------------------------------------------------------

         ee  = e0 * 1.d+6
         tm1 = t1 * 1.d-3
         tm2 = t2 * 1.d-3

         dta = tm2 - tm1

*-----------------------------------------------------------------------

         rr  = 0.404   - 0.290   * exp( -2.78d-4 * ee )
         tz0 = 2.27d-2 + 2.03    * ee**(-0.460)
         g1  = 2.95d-2 + 0.905   * ee**( 0.343)
         g2  = 6.78d-2 + 9.77d-2 * ee**( 0.447)
         s1  = 6.80d-3 + 0.658   * ee**(-0.468)
         s2  = 3.15d-2 + 1.71    * ee**(-0.476)

         tz1 = tz0 + g1 * s2**2
         tz2 = tz0 + g2 * s2**2

*-----------------------------------------------------------------------

         if( ifirst .eq. -1 ) then

             write(*,'(3(1p1g13.5))') rr
             write(*,'(3(1p1g13.5))') s1,s2
             write(*,'(3(1p1g13.5))') g1,g2
             write(*,'(3(1p1g13.5))') tz0*1.e3, tz1*1.e3, tz2*1.e3
             ifirst = 1

         end if

*-----------------------------------------------------------------------

  100    continue

            ran1 = unirn(dummy)
            ran2 = unirn(dummy)

            tt = tm1 + dta * ran1
            dt = tt - tz0

         if( tt .lt. tz0 ) then

               ff = exp( -0.5d0 * ( dt / s1 )**2 )

         else

            if( tt .lt. tz1 ) then

               f1 = exp( -0.5d0 * ( dt / s2 )**2 )

            else

               f1 = exp( 0.5d0 * ( g1 * s2 )**2 - g1 * dt )

            end if

            if( tt .lt. tz2 ) then

               f2 = exp( -0.5d0 * ( dt / s2 )**2 )

            else

               f2 = exp( 0.5d0 * ( g2 * s2 )**2 - g2 * dt )

            end if

               ff = ( 1.d0 - rr ) * f1 + rr * f2

         end if

         if( ran2 .gt. ff ) goto 100

*-----------------------------------------------------------------------

            t0 = tt * 1.d+3

*-----------------------------------------------------------------------

      return
      end

cFURUTA20200206 Not to use unstandard FORTRAN function isNaN
cFURUTA20200206 Instead initusrsors was introduced to avoid uninit vals
c$$$************************************************************************
c$$$*                                                                      *
c$$$      subroutine checkusrsors(ux,uy,uz,uu,uv,uw,ue,uwt,utim,inam,iukf,
c$$$     &                   int1,int2,int3,usx,usy,usz,iflag)
c$$$*                                                                      *
c$$$*        purpose : check variables defined in usrsors are not NaN      *
c$$$*        last modified by T.Furuta on 2019/04/23                        *
c$$$*                                                                      *
c$$$************************************************************************
c$$$
c$$$      implicit none
c$$$      real(8) ux,uy,uz,uu,uv,uw,ue,uwt,utim,usx,usy,usz
c$$$      integer inam,iukf,int1,int2,int3
c$$$      integer iflag
c$$$      logical flags
c$$$
c$$$      iflag=0
c$$$      flags=isnan(ux).or.isnan(uy).or.isnan(uz)
c$$$     &     .or.isnan(uu).or.isnan(uv).or.isnan(uw)
c$$$     &     .or.isnan(ue).or.isnan(uwt).or.isnan(utim)
c$$$     &     .or.isnan(usx).or.isnan(usy).or.isnan(usz)
c$$$      if(flags)iflag=1
c$$$
c$$$      return
c$$$      end

************************************************************************
*                                                                      *
      subroutine initusrsors(ux,uy,uz,uu,uv,uw,ue,uwt,utim,inam,iukf,
     &                   int1,int2,int3,usx,usy,usz)
*                                                                      *
*        purpose : initialize variables before defined in usrsors      *
*        replaced with checkusrsors by T.Furuta on 2020/02/06          *
*                                                                      *
************************************************************************

      implicit none
      real(8) ux,uy,uz,uu,uv,uw,ue,uwt,utim,usx,usy,usz
      integer inam,iukf,int1,int2,int3

      ux = 0.0
      uy = 0.0
      uz = 0.0

      uu = 0.0
      uv = 0.0
      uw = 1.0

      ue = 3000.0

      uwt = 1.0
      utim = 0.0
      inam = 1

      iukf = 2212

      int1 = 0
      int2 = 0
      int3 = 0

      usx = 0.d0
      usy = 0.d0
      usz = 0.d0

      return
      end

************************************************************************
*                                                                      *
      subroutine init_sors
*                                                                      *
*        initialize source variable before OpenMP loop                 *
*        last modified by T.Furuta on 2021/03/24                       *
*                                                                      *
************************************************************************
      use TETRAMOD,only:tetrasorsinit
      implicit real*8 (a-h,o-z)
*-----------------------------------------------------------------------
      include 'param.inc'
*-----------------------------------------------------------------------
      common /isorst/ jstyp(isrc), istyp(isrc), inkf0(isrc), lstyp(isrc)
!$OMP THREADPRIVATE(/isorst/)
      common /itetsor/ ksoutnode,itetreg(isrc)
*-----------------------------------------------------------------------
cFURUTA20210324 tetrasorsinit moved from sors
      do j=1,isrc
       if( jstyp(j) .eq. 24 .or. jstyp(j) .eq. 25 )then
        if(ksoutnode.eq.0)then
         call tetrasorsinit(kvlmax)
         ksoutnode=1
        endif
       endif
      enddo
      return
      end

************************************************************************
*                                                                      *
      subroutine dumpsors(j,nj,nmulti,pmlpr,wtnorm,wtspfs,wtduct)
*                                                                      *
*        subroutine for dump source                                    *
*        last modified by T.Furuta on 2021/07/26                       *
*                                                                      *
************************************************************************
      use moddas_region
      use MMBANKMOD !FURUTA
      use mod_ompparallel
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'err.inc'

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /mpi00/  npe, me
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /inout/  in,io
      common /tcntl/  icntl, inucr
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /iradkk/ randkk,irskip
      common /taliin/ rsouin, nzztin, nrgnin
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /paraj/  mstz(300), parz(300)

      common /srsph/  pmphs
!$OMP THREADPRIVATE(/srsph/)

*-----------------------------------------------------------------------

      common /isocor/ iscorr, itcorr, imlwt(isrc)
      common /isorst/ jstyp(isrc), istyp(isrc), inkf0(isrc), lstyp(isrc)
!$OMP THREADPRIVATE(/isorst/)
      common /isorsp/ sx0(isrc), sy0(isrc), sz0(isrc), sx1(isrc),
     &                sy1(isrc), sz1(isrc), sr0(isrc), se0(isrc),
     &                sdir(isrc), srx(isrc), sry(isrc), swem(isrc),
     &                sphi(isrc), sdom(isrc), swt0(isrc)
      common /isorsg/ seg0(isrc),seg1(isrc),seg2(isrc),seg3(isrc),
     &                set0(isrc),set1(isrc),set2(isrc),set3(isrc),
     &                jetyp(isrc),jptyp(isrc)
      common /isorsf/ isorf(isrc),lsfile(isrc), sfile(isrc)
      character sfile*100
      common /isorsm/ stm0(isrc),stmw(isrc),stmc(isrc),stmd(isrc),
     &                jttyp(isrc),jttpn(isrc)
      common /isorsr/ nsmx(isrc), nsrn(isrc), nsrc(isrc)
      common /isorsc/ isort(isrc,4), rsort(isrc,13)
      common /isorsd/ isdmp(isrc,0:30), jsdmp(isrc,0:30)
      common /isorpn/ ssx(isrc), ssy(isrc), ssz(isrc)
!$OMP THREADPRIVATE(/isorpn/)

*-----------------------------------------------------------------------

      dimension ksou1(isrc)
      data ksou1 /isrc*0/
      save ksou1
!$OMP THREADPRIVATE(ksou1)
      dimension ksou1node(isrc)
      data ksou1node /isrc*0/ !FURUTA
      save ksou1node     !FURUTA

      data iserr /0/
      save iserr
!$OMP THREADPRIVATE(iserr)
      dimension     idas(1)
      equivalence ( das, idas )

      dimension dmpd(30)
      dimension pmlpr(isrc)
*-----------------------------------------------------------------------
      common /bnkmem/ maxbnk, maxbn2, rtrckflp
      real(8) dmpmulti
      integer idmpmode,ibchjmp,idmpjmp
      common /stat2/ dmpmulti,idmpmode,ibchjmp,idmpjmp(2)
      integer nmulti
      real(8) dmpd_old(30),wdump,dmul,rrr
      save dmpd_old
      integer jpsf
      common /stat3/ jpsf
      integer idata
      real(4) fdata(7)
      real(8) thetadata,phidata
*-----------------------------------------------------------------------

      common /cntsrc/icntsrc(3,isrc) ! initial counter

      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

      integer ipompl
      ipompl = ipomp+1

*-----------------------------------------------------------------------
*     start of source
*-----------------------------------------------------------------------

      ifirst = 0                !FURUTA20150515
      wtnorm=1.0d0
      wtspfs=1.0d0
      wtduct=1.0d0

*-----------------------------------------------------------------------
*     weight setting for multiplication by dmpmulti
*     roussian roulette for number of decimal places in dmpmulti
* FURUTA20150515
*-----------------------------------------------------------------------
      if(dmpmulti.gt.0.0)then
       nmulti=int(dmpmulti)
       dmul=dmpmulti-dble(nmulti)
       rrr = unirn(dummy)
       if(rrr.lt.dmul)then
        nmulti=nmulti+1
       endif
       if(dmpmulti.ge.1.0)then
        wdump=1.0d0/dble(nmulti)
       else
        wdump=1.0d0/dmpmulti
       endif
      endif

*-----------------------------------------------------------------------
*     number of try for finding region
*-----------------------------------------------------------------------

      itsrg = 0

 2000 continue

      wtnorm = 1.0d0

      t(ibkt+nj,ipompl) = 0.0d0
      ncnt(ibknct+1,nj,ipompl) = icntsrc(1,j) ! T.Sato 2019/05/19
      ncnt(ibknct+2,nj,ipompl) = icntsrc(2,j) ! T.Sato 2019/05/19
      ncnt(ibknct+3,nj,ipompl) = icntsrc(3,j) ! T.Sato 2019/05/19
      do itt =1, 3
       if ( nj.eq.1 ) then
        nctmxsr(itt) = icntsrc(itt,j)
       else if ( icntsrc(itt,j).gt.nctmxsr(itt) ) then
        nctmxsr(itt) = icntsrc(itt,j)
       end if
      end do

      name0 = 1

      if(ibchjmp.ne.0)then      ! No further dump source
       nmulti=0
       goto 5900
      endif

*-----------------------------------------------------------------------
*     idmpmode = 0:  event by event statistical processing is OFF
*-----------------------------------------------------------------------
      if(idmpmode.eq.0)then     !FURUTA20150515

       if( ksou1(j) .eq. 0 .and. ksou1node(j) .eq. 0 ) then !FURUTA

        if( irskip .ne. 0 ) then

         jrskip = iabs( irskip )

         do i = 1, jrskip

 787      continue

          if( isdmp(j,0) .gt. 0 ) then

           if(jpsf.eq.0)then    !FURUTA20201110
            read(isorf(j),iostat=ios,err=790)
     &           ( dmpd(isdmp(j,k)), k = 1, abs( isdmp(j,0) ) )
            if( ios .eq. -1 ) goto 788
           else                 !FURUTA20201110
            read(isorf(j),iostat=ios,err=790)
     &           idata,(fdata(k),k=1,4+jpsf)
            if( ios .eq. -1) goto 788
           endif

          else

           read(isorf(j),'(30(1p1d24.15))',iostat=ios,err=790)
     &          ( dmpd(isdmp(j,k)), k = 1, abs( isdmp(j,0) ) )
           if( ios .eq. -1 ) goto 788

          end if

          goto 789
 788      rewind isorf(j)
          goto 787
 790      continue

          iserr = iserr + 1
          write( *,'(''Error in dump file no ='',i5)') iserr
          write(io,'(''Error in dump file no ='',i5)') iserr

 789      continue

         end do

        end if

       end if

 687   continue

       if( isdmp(j,0) .gt. 0 ) then

        if(jpsf.eq.0)then       !FURUTA20201110
         read(isorf(j),iostat=ios,err=690)
     &        ( dmpd(isdmp(j,k)), k = 1, abs( isdmp(j,0) ) )
         if( ios .eq. -1 ) goto 688
        else                    !FURUTA20201110
         read(isorf(j),iostat=ios,err=690)
     &        idata,(fdata(k),k=1,4+jpsf)
         if( ios .eq. -1) goto 688
         dmpd(1)=dble(idata)
         dmpd(2)=dble(fdata(1))
         dmpd(3)=dble(fdata(2))
         thetadata=dble(fdata(3))
         phidata=dble(fdata(4))
         dmpd(5)=sin(fdata(3))*cos(fdata(4))
         dmpd(6)=sin(fdata(3))*sin(fdata(4))
         dmpd(7)=cos(fdata(3))
         dmpd(8)=dble(fdata(5))
         if(jpsf.eq.2)dmpd(9)=dble(fdata(6))
        endif

       else

        read(isorf(j),*,iostat=ios,err=690) ! T.Sato 2025/04/08 change to free format
     &       ( dmpd(isdmp(j,k)), k = 1, abs( isdmp(j,0) ) )
        if( ios .eq. -1 ) goto 688

       end if

       goto 689
 688   if(dmpmulti.eq.0.0)then  !FURUTA20150515
        rewind isorf(j)
        goto 687
       else

        ErrCha = ''
        ErrID = 'L:3640/R:dumpsors/F:sors.f' !E03_105_001
        call ErrWrite(ErrID,ErrCha)

        write(*,'(/
     &     ''** ERROR in sours, dump source file ended''/
     &     '' at NOBCH ='',i5,'' NOCAS='',i5)')
     &       nobch,nocas
        write(io,'(/
     &     ''** ERROR in sours, dump source file ended''/
     &     '' at NOBCH ='',i5,'' NOCAS='',i5)')
     &       nobch,nocas
        call parastop( 711 )
       endif
c------------------
 690   continue

       iserr = iserr + 1
       write(ErrCha,'(''Error in dump file no ='',i5)') iserr
       ErrID = 'L:3658/R:dumpsors/F:sors.f' !E03_105_002
       call ErrWrite(ErrID,ErrCha)
       call ErrWriteIO(ErrID,ErrCha,io)

 689   continue

*-----------------------------------------------------------------------
*     idmpmode = 1:  event by event statistical processing is ON
* FURUTA20150515
*-----------------------------------------------------------------------

      else

       if(irskip .ne. 0)then
        ErrCha = ''
        ErrID = 'L:3673/R:dumpsors/F:sors.f' !E03_105_003
        call ErrWrite(ErrID,ErrCha)

        write(*,'(/
     &     ''** ERROR in sours, irskip should be 0''/
     &     '' in event dump source (idmpmode=1)'')')
        write(io,'(/
     &     ''** ERROR in sours, irskip should be 0''/
     &     '' in event dump source (idmpmode=1)'')')
        call parastop( 710 )
       endif

       if(ifirst.gt.0.or.ksou1node(j).eq.0)then

 587    continue
        if( isdmp(j,0) .gt. 0 ) then
         read(isorf(j),iostat=ios,err=590)
     &        ( dmpd(isdmp(j,k)), k = 1, abs( isdmp(j,0) ) )
         if( ios .eq. -1 ) goto 588
        else
         read(isorf(j),'(30(1p1d24.15))',iostat=ios,err=590)
     &        ( dmpd(isdmp(j,k)), k = 1, abs( isdmp(j,0) ) )
         if( ios .eq. -1 ) goto 588
        end if
        goto 589
 588    nj=nj-1
        ibchjmp=-1
        goto 5900
 590    continue
        iserr = iserr + 1
        write(ErrCha,'(''Error in dump file no ='',i5)') iserr
        ErrID = 'L:3704/R:dumpsors/F:sors.f' !E03_105_004
        call ErrWrite(ErrID,ErrCha)
        call ErrWriteIO(ErrID,ErrCha,io)
 589    continue
        if(ksou1node(j).ne.0)then
         if(dmpd(18).ne.dmpd_old(18)
     &        .or.dmpd(19).ne.dmpd_old(19))then
          if(dmpd(19).ne.dmpd_old(19))ibchjmp=1 ! Jumper for next batch
          nj=nj-1
          dmpd_old=dmpd
          goto 5900
         endif
        endif
        dmpd_old=dmpd

       else

        dmpd=dmpd_old

       endif

       ifirst=1

      endif

*-----------------------------------------------------------------------

      wtnorm = swt0(j)
      if(nmulti.eq.0) goto 5800 !FURUTA20150515
      e(ibke+nj,ipompl) = se0(j)

*-----------------------------------------------------------------------

      if( jsdmp(j,1)  .gt. 0 ) then
       inkf0(j) = nint( dmpd(1) )
       istyp(j) = kftp(inkf0(j))
      end if

      if( jsdmp(j,2)  .gt. 0 ) then
       x(ibkx+nj,ipompl) = dmpd(2)
      else
       x(ibkx+nj,ipompl) =
     &      unirn(dummy) * ( sx1(j) - sx0(j) ) + sx0(j)
      end if

      if( jsdmp(j,3)  .gt. 0 ) then
       y(ibky+nj,ipompl) = dmpd(3)
      else
       y(ibky+nj,ipompl) =
     &      unirn(dummy) * ( sy1(j) - sy0(j) ) + sy0(j)
      end if

      if( jsdmp(j,4)  .gt. 0 ) then
       z(ibkz+nj,ipompl) = dmpd(4)
      else
       z(ibkz+nj,ipompl) =
     &      unirn(dummy) * ( sz1(j) - sz0(j) ) + sz0(j)
      end if

      if( jsdmp(j,5)  .gt. 0 ) vx = dmpd(5)
      if( jsdmp(j,6)  .gt. 0 ) vy = dmpd(6)
      if( jsdmp(j,7)  .gt. 0 ) vz = dmpd(7)
      if( jsdmp(j,8)  .gt. 0 ) e(ibke+nj,ipompl) = dmpd(8)
      if( jsdmp(j,9)  .gt. 0 ) wtnorm = dmpd(9)
      if( jsdmp(j,10) .gt. 0 ) t(ibkt+nj,ipompl) = dmpd(10)
      if( jsdmp(j,11) .gt. 0 )
     &     ncnt(ibknct+1,nj,ipompl) = nint(dmpd(11))
      if( jsdmp(j,12) .gt. 0 )
     &     ncnt(ibknct+2,nj,ipompl) = nint(dmpd(12))
      if( jsdmp(j,13) .gt. 0 )
     &     ncnt(ibknct+3,nj,ipompl) = nint(dmpd(13))
      if( jsdmp(j,14) .gt. 0 ) ssx(j) = dmpd(14)
      if( jsdmp(j,15) .gt. 0 ) ssy(j) = dmpd(15)
      if( jsdmp(j,16) .gt. 0 ) ssz(j) = dmpd(16)

      if( jsdmp(j,17) .gt. 0 ) name0 = nint(dmpd(17))

      if( sdir(j) .gt. 1000.0 ) then

       av = sqrt( vx**2 + vy**2 + vz**2 )

       u(ibku+nj,ipompl) = vx / av
       v(ibkv+nj,ipompl) = vy / av
       w(ibkw+nj,ipompl) = vz / av

      else ! temporary determine direction to avoid error, actual direction will be determined later, T.Sato 2021/12/07

       u(ibku+nj,ipompl) = 0.0
       v(ibkv+nj,ipompl) = 0.0
       w(ibkw+nj,ipompl) = 1.0

      end if


*-----------------------------------------------------------------------
*     u(ibku+1), v(ibkv+1), w(ibkw+1)
*-----------------------------------------------------------------------

      aw0 = 1.d0

*-----------------------------------------------------------------------
*     weight factor
*-----------------------------------------------------------------------

      wtnorm = wtnorm * aw0 * pmlpr(j)

*-----------------------------------------------------------------------
*     e(ibke+1)
*-----------------------------------------------------------------------

      ew0 = 1.0d0

      if( jetyp(j) .gt. 0 )then

       jin = j

       call engbin(eo,ew0,jin,idummy)

       e(ibke+nj,ipompl) = eo
       wtnorm = wtnorm * ew0

      end if


*-----------------------------------------------------------------------
*     nty(ibknty+1), nkf(ibknkf+1), name(ibknam+1),
*-----------------------------------------------------------------------

      nty(ibknty+nj,ipompl)  = istyp(j)
      nkf(ibknkf+nj,ipompl)  = inkf0(j)
      name(ibknam+nj,ipompl) = name0

*-----------------------------------------------------------------------
*     for nucleus, input energy is MeV / nuleon
*-----------------------------------------------------------------------

      if( nty(ibknty+nj,ipompl) .ge. 15 ) then

       e(ibke+nj,ipompl) = e(ibke+nj,ipompl) *
     &      ibryf(nty(ibknty+nj,ipompl),
     &      nkf(ibknkf+nj,ipompl))

      end if

*-----------------------------------------------------------------------
*     t(ibkt+1)
*-----------------------------------------------------------------------


      tw0 = 1.0d0

      if( jttyp(j) .eq.   1 .or.
     &     jttyp(j) .eq.   2 .or.
     &     jttyp(j) .eq.   3 .or.
     &     jttyp(j) .eq.   4 .or.
     &     jttyp(j) .eq.   5 .or.
     &     jttyp(j) .eq.   6 .or.
     &     jttyp(j) .eq. 100 .or.
     &     ( jstyp(j) .eq. 17  .and. jttyp(j) .gt. 0 ) .or.
     &     ( jstyp(j) .eq. 100 .and. jttyp(j) .gt. 0 ) ) then

       eo = e(ibke+nj,ipompl)
       jin = j

       call timbin(to,tw0,eo,jin)

       t(ibkt+nj,ipompl) = to
       wtnorm = wtnorm * tw0

      end if

*-----------------------------------------------------------------------
*     initial flag
*-----------------------------------------------------------------------

      if( ksou1(j) .eq. 0 ) ksou1(j) = 1

*-----------------------------------------------------------------------
*     transformation
*-----------------------------------------------------------------------

      if( isort(j,4) .ne. 0 .and. (iscorr .le. 1 .or. ! 2018/09/04 Ogawa. Apply if source positions are uncorrelated
     &     (iscorr .ge. 2 .and. j .eq. 1) )) then ! do not apply second and subsequent if source positions is that of the first one

       call trnsxv(x(ibkx+nj,ipompl),y(ibky+nj,ipompl),
     &      z(ibkz+nj,ipompl),
     &      ux,uy,uz,isort(j,4))

       call trnsuv(u(ibku+nj,ipompl),v(ibkv+nj,ipompl),
     &      w(ibkw+nj,ipompl),
     &      uu,uv,uw,isort(j,4))

       x(ibkx+nj,ipompl) = ux
       y(ibky+nj,ipompl) = uy
       z(ibkz+nj,ipompl) = uz

       u(ibku+nj,ipompl) = uu
       v(ibkv+nj,ipompl) = uv
       w(ibkw+nj,ipompl) = uw

      end if

*-----------------------------------------------------------------------
*     find initial region
*-----------------------------------------------------------------------

      if( icntl .eq. 0 .or. icntl .eq. 5 .or.
     &     icntl .eq. 6 .or.
     &     icntl .eq. 14 .or.
     &     icntl .eq. 15 ) then

       isrr  = 0

 700   continue

       ici   = 0
       mark  = 1
       markp = 0

       call gomsor(x(ibkx+nj,ipompl),y(ibky+nj,ipompl),
     &      z(ibkz+nj,ipompl),
     &      u(ibku+nj,ipompl),v(ibkv+nj,ipompl),
     &      w(ibkw+nj,ipompl),
     &      nmed(ibknmd+nj,ipompl),iblz(ibkblz+nj,ipompl),
     &      mark,markp,ici)

*-----------------------------------------------------------------------
*           check the point on boundary
*-----------------------------------------------------------------------

       markp = 1

       call gomdis(dpr,
     &      x(ibkx+nj,ipompl),y(ibky+nj,ipompl),
     &      z(ibkz+nj,ipompl),
     &      u(ibku+nj,ipompl),v(ibkv+nj,ipompl),
     &      w(ibkw+nj,ipompl),
     &      mark,markp,
     &      nmed(ibknmd+nj,ipompl),iblz(ibkblz+nj,ipompl))

       markp = 0

       if( mark .eq. -6 ) then

        if( isrr .eq. 0 ) then

         x(ibkx+nj,ipompl) = x(ibkx+nj,ipompl) +
     &        u(ibku+nj,ipompl) * parz(28)
         y(ibky+nj,ipompl) = y(ibky+nj,ipompl) +
     &        v(ibkv+nj,ipompl) * parz(28)
         z(ibkz+nj,ipompl) = z(ibkz+nj,ipompl) +
     &        w(ibkw+nj,ipompl) * parz(28)

         isrr = isrr + 1
         goto 700

        else

         ErrCha = ''
         ErrID = 'L:3963/R:dumpsors/F:sors.f' !E03_106_001
         call ErrWrite(ErrID,ErrCha)

         write(*,'(/
     &    ''** ERROR in sors, source is not in any cell'')')
         write(*,*)  x(ibkx+nj,ipompl), y(ibky+nj,ipompl),
     &        z(ibkz+nj,ipompl)

         write(io,'(/
     &     ''** ERROR in sors, source is not in any cell'')')
         write(io,*) x(ibkx+nj,ipompl), y(ibky+nj,ipompl),
     &        z(ibkz+nj,ipompl)
         call parastop( 602 )

        end if

       end if

*-----------------------------------------------------------------------
*        region selection ( maxmum trial is ntmax )
*-----------------------------------------------------------------------

       if( nsrn(j) .gt. 0 .and. jstyp(j) .ne. 12
     &      .and. jstyp(j) .ne. 24 !FURUTA20180115
     &      .and. jstyp(j) .ne. 25 !FURUTA20180115
     &      ) then

        idsm = iaddress_nsrn(j)
        jdsm = 0

        jj  = 0

        jdsm = jdsm + 1
        ntrn = idas_nsrn(idsm+jdsm)
        jdsm = jdsm + 1
        mtrn = idas_nsrn(idsm+jdsm)

        do ii = 1, ntrn

         call tregck(iblz1,ilev1,ilat1,
     &        mtrn,idas_nsrn(idsm+jdsm+1),jj,icc)

         if( icc .ne. 0 ) goto 500

        end do

        itsrg = itsrg + 1

        if( itsrg .gt. nsmx(j) ) then
         ErrCha = ''
         ErrID = 'L:4013/R:dumpsors/F:sors.f' !E03_107_001
         call ErrWrite(ErrID,ErrCha)

         write(io,'(/'' Error : region selection failed,'',
     &    '' in source.''/
     &    '' max. trial ntmax = '',i5)') nsmx(j)

         write(*,'(/'' Error : region selection failed,'',
     &    '' in source.''/
     &    '' max. trial ntmax = '',i5)') nsmx(j)


         call parastop( 889 )

        end if

        goto 2000

 500    continue

       end if

      end if

*-----------------------------------------------------------------------
*     phase of this event
*-----------------------------------------------------------------------

      pmphs = 2.0d0 * pi * unirn(dummy)

*-----------------------------------------------------------------------
*     sorce number j
*-----------------------------------------------------------------------

      nsos(ibksos+nj,ipompl) = j

*-----------------------------------------------------------------------
*     charge state
*-----------------------------------------------------------------------

      if( lstyp(j) .gt. -1000.d0 ) then

       nzst(ibkzst+nj,ipompl) = lstyp(j)

      else

       nzst(ibkzst+nj,ipompl) = ichgf(istyp(j),inkf0(j))

      end if

*-----------------------------------------------------------------------
*     initial spin,  ssr = 0 ; no defined
*-----------------------------------------------------------------------

      ssr = ssx(j)**2 + ssy(j)**2 + ssz(j)**2

      if( ssr .gt. 1.d-8 ) then

       spx(ibkspx+nj,ipompl) = ssx(j) / ssr
       spy(ibkspy+nj,ipompl) = ssy(j) / ssr
       spz(ibkspz+nj,ipompl) = ssz(j) / ssr

      else

       spx(ibkspx+nj,ipompl) = 0.d0
       spy(ibkspy+nj,ipompl) = 0.d0
       spz(ibkspz+nj,ipompl) = 0.d0

      end if

*-----------------------------------------------------------------------
*        wt(ibkwt+nj)
*-----------------------------------------------------------------------

      wt(ibkwt+nj,ipompl) = wtnorm * wtduct

*-----------------------------------------------------------------------
*     nmulti : multiplicity for dump source
* FURUTA20150515
*-----------------------------------------------------------------------
      if(dmpmulti.gt.0.0)then
       wt(ibkwt+nj,ipompl)=wt(ibkwt+nj,ipompl)*wdump
       if(nmulti.gt.1)then
        if(nj+nmulti-1.gt.maxbnk)then
         write(6,'(/''**** Fatal ERROR: lack of bank'',
     &     '' space ****''/
     &     ''     please increase maxbnk'')')
         call parastop( 837 )
        endif
        do k=nj+1,nj+nmulti-1
         nkf(ibknkf+k,ipompl)=nkf(ibknkf+nj,ipompl)
         x(ibkx+k,ipompl)=x(ibkx+nj,ipompl)
         y(ibky+k,ipompl)=y(ibky+nj,ipompl)
         z(ibkz+k,ipompl)=z(ibkz+nj,ipompl)
         u(ibku+k,ipompl)=u(ibku+nj,ipompl)
         v(ibkv+k,ipompl)=v(ibkv+nj,ipompl)
         w(ibkw+k,ipompl)=w(ibkw+nj,ipompl)
         e(ibke+k,ipompl)=e(ibke+nj,ipompl)
         wt(ibkwt+k,ipompl)=wt(ibkwt+nj,ipompl)
         t(ibkt+k,ipompl)=t(ibkt+nj,ipompl)
         ncnt(ibknct+1,k,ipompl)=ncnt(ibknct+1,nj,ipompl)
         ncnt(ibknct+2,k,ipompl)=ncnt(ibknct+2,nj,ipompl)
         ncnt(ibknct+3,k,ipompl)=ncnt(ibknct+3,nj,ipompl)
         spx(ibkspx+k,ipompl)=spx(ibkspx+nj,ipompl)
         spy(ibkspy+k,ipompl)=spy(ibkspy+nj,ipompl)
         spz(ibkspz+k,ipompl)=spz(ibkspz+nj,ipompl)
         name(ibknam+k,ipompl)=name(ibknam+nj,ipompl)
         nty(ibknty+k,ipompl)=nty(ibknty+nj,ipompl)
         iblz(ibkblz+k,ipompl)=iblz(ibkblz+nj,ipompl)
         nmed(ibknmd+k,ipompl)=nmed(ibknmd+nj,ipompl)
         nzst(ibkzst+k,ipompl)=nzst(ibkzst+nj,ipompl)
         nsos(ibksos+k,ipompl)=nsos(ibksos+nj,ipompl)
         itetpos(ibtetpos+k,ipompl)=itetpos(ibtetpos+nj,ipompl)

        enddo
        nj=nj+nmulti-1
       endif
      endif

 5800 continue

      if( ksou1node(j) .eq. 0) ksou1node(j) = 1

*-----------------------------------------------------------------------
*     next dump souce for idmpmode=1
* FURUTA20150515
*-----------------------------------------------------------------------
      if( idmpmode .eq. 1 )then
       nj=nj+1
       goto 2000
      endif

 5900 continue

      return
      end subroutine dumpsors
