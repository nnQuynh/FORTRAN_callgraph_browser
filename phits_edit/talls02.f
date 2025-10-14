************************************************************************
*                                                                      *
      subroutine tsufreg(ncol,m,np,nr,ne,na,nt,nm,mr,kr,eb,ab,tb,tr,
     &                   trEVENT)
*                                                                      *
*       region crossing surface tally of current spectrum and flux.    *
*       last modified by K.Niita on 2015/10/30                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*             10 : geometry boundary crossing                          *
*             12 : termination by escape or leakage                    *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use partmod, only: itmxpt, itpan, itpat, jtpat ! frtati 2021/10/05
*-----------------------------------------------------------------------

      implicit double precision(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /tlcost/ costha, uang(3), nsurf
!$OMP THREADPRIVATE(/tlcost/)
      common /paraj/  mstz(300), parz(300)

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)

      common /tall13/ itrcn(itlmax), itrcr(itlmax), itrca(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /tall73/ itenclo(itlmax), itangform(itlmax)
      dimension ud(3)

      common /tall82/ itcnth(9,itlmax)

      common /tall50/ itlmt(itlmax), ite2l(itlmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

*-----------------------------------------------------------------------

      dimension   kr(mr+1)  ! T.Sato 2022/10/28
      dimension   eb(ne+1)
      dimension   ab(na+1)
      dimension   tb(nt+1)
      dimension   tr(np,ne,na,nt,nr,nm,2)
      dimension   trEVENT(np,ne,na,nt,nr,nm)       !OBINATA(2012.7.11): as Ct
      real(8),allocatable,save:: tr0(:,:,:,:,:,:)  !OBINATA(2012.7.11): as C


      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt
      dimension facm(6)
      dimension dmpd(30)

      common /stat / istdev, irestart, ireschk
      common /cparm/ maxbch,maxcas

      common /egs5cmn10/denstepold,denstepnew,deinit
      real*8            denstepold,denstepnew,deinit
!$OMP THREADPRIVATE(/egs5cmn10/)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

*-----------------------------------------------------------------------
      integer idmpomp !FURUTA20150427
      common /idmpomp0/idmpomp !FURUTA20150427
*-----------------------------------------------------------------------

      data grab,hgrab / 0.00349d+0, 0.001745d+0 /

*-----------------------------------------------------------------------
* check of history counter
*-----------------------------------------------------------------------
      ihistcount = 0 ! history counter index, T.Sato 2022/12/20 for speed up
         if ( ncol.eq.0 .or. ncol.eq.4 ) then
           do i = 1, 3
             if( itcnth(i,m) .eq. 1 ) then
               if( ncntmx(i) .lt. itcnth(i*2+2,m) .or.
     &             ncntmx(i) .gt. itcnth(i*2+3,m) ) then
                 ihistcount = 1
               end if
             end if
           end do
         end if

*-----------------------------------------------------------------------
*        source particle or end of batch ( in case of istdev = 2 )
*
* OBINATA(2012.7.11): change tr(,,,,,3) to trEVENT(,,,,)
*-----------------------------------------------------------------------

         if (( ncol .eq. 0 .or. ncol .eq. 4 )
     &                              .and. istdev .eq. 2) then
           if ((nocas.gt.1.or.ncol.eq.0) .and. ihistcount.ne.1 ) then

             tr(:,:,:,:,:,:,1) = tr(:,:,:,:,:,:,1)
     &                         + trEVENT(:,:,:,:,:,:)
             tr(:,:,:,:,:,:,2) = tr(:,:,:,:,:,:,2)
     &                         + trEVENT(:,:,:,:,:,:) ** 2

! sumover
            call tsufreg_sumover(m,1,
     &                 np,    ne,    na,   nt,  nr,   nm,    trEVENT)



           end if

           trEVENT(:,:,:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------
*        end of batch ( in case of istdev = 1 )
*
* OBINATA(2012.7.11): modificate for thread parallel
*-----------------------------------------------------------------------

         if ( ncol .eq. 0 .and. istdev .eq. 1) then
!$OMP MASTER
             allocate( tr0(np,ne,na,nt,nr,nm) )
             tr0(:,:,:,:,:,:) = 0.d0

!$OMP END MASTER
!$OMP BARRIER
!$OMP CRITICAL (tsufreg_crit_ist1)
             tr0(:,:,:,:,:,:) = tr0(:,:,:,:,:,:) + trEVENT(:,:,:,:,:,:)
!$OMP END CRITICAL (tsufreg_crit_ist1)
!$OMP BARRIER
!$OMP MASTER
             tr(:,:,:,:,:,:,1) = tr(:,:,:,:,:,:,1)
     &                         + tr0(:,:,:,:,:,:) / maxcas
             tr(:,:,:,:,:,:,2) = tr(:,:,:,:,:,:,2)
     &                         + ( tr0(:,:,:,:,:,:) / maxcas ) ** 2

! sumover
            call tsufreg_sumover(m,maxcas,
     &                 np,    ne,    na,   nt,  nr,   nm,    tr0)

             deallocate( tr0 )
!$OMP END MASTER

           trEVENT(:,:,:,:,:,:) = 0

         end if


*-----------------------------------------------------------------------
*        check of ncol and backward ( it is impossible )
*-----------------------------------------------------------------------

            if( ncol .ne. 10 .and. ncol .ne. 12 ) return

            if( itout(m) .eq. 4 .or. itout(m) .eq. 7 ) return

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncnt(ibknct+i,no,ipomp+1) .lt. itcnt(i*2+2,m) .or.
     &                ncnt(ibknct+i,no,ipomp+1) .gt. itcnt(i*2+3,m) )
     &                      return

               end if

            end do

*-----------------------------------------------------------------------
*        material for LET
*-----------------------------------------------------------------------

               if( itlmt(m) .gt. 0 ) then

                  lmat = idnm( itlmt(m) )

               else if( itlmt(m) .eq. 0 ) then

                  lmat = mat

               else

                  lmat = -idnm( -itlmt(m) )

               end if

*-----------------------------------------------------------------------
*        check of particles
*-----------------------------------------------------------------------

            call pcheck(m,np,ityp,ktyp,jtyp,ipn,ips)

               if( ipn .eq. 0 ) return

*-----------------------------------------------------------------------
cc H.Iwase 2015/4/2 for electrons
cc                  change the step energy to coutinuous energy
*-----------------------------------------------------------------------

      if( ityp.eq.12 .or. ityp.eq.13 )then

         e1 =  e(ibke  +no,ipomp+1) - deinit  + denstepold
         e2 = ec(ibkec +no,ipomp+1) - deinit  + denstepnew

         if( e1 .lt. eb(1)    ) return
         if( e2 .ge. eb(ne+1) ) return

         tparti = abs(t(ibkt+no,ipomp+1))
         tpartf = abs(tc(ibktc+no,ipomp+1))
         if( tparti .ge. tb(nt+1) ) return
         if( tpartf .lt. tb(1) ) return

          e(ibke +no,ipomp+1) = e1
         ec(ibkec+no,ipomp+1) = e2

      endif

*-----------------------------------------------------------------------

      if( iMeVperu.eq.1 .and. ityp.ge.15 .and. ityp.le.19 ) then
           ebm = ktyp - ktyp / 1000000 * 1000000
      else
           ebm = 1.d0
      end if

*-----------------------------------------------------------------------
*        check of energy
*-----------------------------------------------------------------------

      if (ite2l(m) .eq. 0) then   ! ccse 2022/08/31

         if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &       itout(m) .eq. 10) then  ! T.Sato 2021/05/05

               if( ec(ibkec+no,ipomp+1) .lt. eb(1)*ebm ) goto 999
               if( ec(ibkec+no,ipomp+1) .ge. eb(ne+1)*ebm ) goto 999

            do i = 1, ne

               if( ec(ibkec+no,ipomp+1) .ge. eb(i)*ebm .and.
     &             ec(ibkec+no,ipomp+1) .lt. eb(i+1)*ebm ) goto 30

            end do

   30          ie = i

         else

               ie = 1

         end if

      else if (ite2l(m) .eq. 1) then   ! convert to energy to LET

         if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &       itout(m) .eq. 10) then    ! T.Sato 2021/05/05

            call dedxas(ec(ibkec+no,ipomp+1),dedx,lmat,
     &                  ityp,ktyp,jtyp,rtyp)
            dedx = dedx /10.0d0

            if( dedx .lt. eb(1) ) goto 999
            if( dedx .ge. eb(ne+1) ) goto 999

            do i = 1, ne
               if( dedx .ge. eb(i) .and.
     &             dedx .lt. eb(i+1) ) goto 35
            end do

   35       ie = i

         else

            ie = 1

         end if
      end if

*-----------------------------------------------------------------------
*         FM factor
*-----------------------------------------------------------------------

               icli = idgr(iblz1)

               call fmfac(m,icli,ec(ibkec+no,ipomp+1),oldwt,facm)

*-----------------------------------------------------------------------
*        check of angle
*-----------------------------------------------------------------------

         if( itout(m) .ge. 8 ) then ! T.Sato 2021/05/05

            if( itangform(m).ge.1 .and. itangform(m).le.3 ) then

               dis = ( xc(ibkxc+no,ipomp+1) - x(ibkx+no,ipomp+1) )**2
     &             + ( yc(ibkyc+no,ipomp+1) - y(ibky+no,ipomp+1) )**2
     &             + ( zc(ibkzc+no,ipomp+1) - z(ibkz+no,ipomp+1) )**2
               dis = dsqrt(dis)

               ud(1) = ( xc(ibkxc+no,ipomp+1)
     &                  - x(ibkx+no, ipomp+1) ) / dis
               ud(2) = ( yc(ibkyc+no,ipomp+1)
     &                  - y(ibky+no, ipomp+1) ) / dis
               ud(3) = ( zc(ibkzc+no,ipomp+1)
     &                  - z(ibkz+no, ipomp+1) ) / dis

               cst = ud(itangform(m))

            else

               cst = costha

            endif

               if( itaty(m) .gt. 0 ) then

                  if( cst .lt. ab(1) ) goto 999
                  if( cst .gt. ab(na+1) ) goto 999

               else

                  if( cst .lt. cos( ab(na+1) / 180.d0 * pi ) ) goto 999
                  if( cst .gt. cos( ab(1) / 180.d0 * pi ) ) goto 999

               end if

            do i = 1, na

               if( itaty(m) .gt. 0 ) then

                  if( cst .ge. ab(i) .and.
     &                cst .le. ab(i+1) ) goto 40

               else

                  if( cst .ge. cos( ab(i+1) / 180.d0 * pi ) .and.
     &                cst .le. cos( ab(i) / 180.d0 * pi ) ) goto 40

               end if

            end do

   40          ia = i

         else

               ia = 1

         end if

*-----------------------------------------------------------------------
*           check of time
*-----------------------------------------------------------------------

               tpart = abs(tc(ibktc+no,ipomp+1))

               if( tpart .lt. tb(1) ) goto 999
               if( tpart .ge. tb(nt+1) ) goto 999

*-----------------------------------------------------------------------
*           time
*-----------------------------------------------------------------------

            do i = 1, nt

               if( tpart .ge. tb(i) .and.
     &             tpart .lt. tb(i+1) ) goto 50

            end do

   50          it = i

*-----------------------------------------------------------------------
*        surface crossing flux tally
*-----------------------------------------------------------------------

         if( itout(m) .eq. 1
     &   .or.itout(m).eq.10.or.itout(m).eq.11) then ! T.Sato 2021/05/05
               costhe = abs( costha )

               if( costhe .le. grab ) costhe = hgrab

               wt1 = oldwt / costhe
               wt2 = wt1**2

*-----------------------------------------------------------------------
*        surface crossing current tally
*-----------------------------------------------------------------------

         else

               wt1 = oldwt
               wt2 = oldwt**2

         end if

*-----------------------------------------------------------------------
*        check of region
*-----------------------------------------------------------------------

                     j = 0

      do 1000 ir = 1, nr

                     j = j + 1
                     ntrn = kr(j)
                     j = j + 1
                     mtrn = kr(j)

                     jj  = 0

                  do k = 1, ntrn

                     call tregck(iblz1,ilev1,ilat1,
     &                           mtrn,kr(j+1),jj,ic1)

                  end do

                     j = j + mtrn

                     j = j + 1
                     ntrn = kr(j)
                     j = j + 1
                     mtrn = kr(j)

               if( ic1 .ne. 0 ) then

                     jj  = 0

                  do k = 1, ntrn

                     call tregck(iblz2,ilev2,ilat2,
     &                           mtrn,kr(j+1),jj,ic2)

                     if( ic2 .ne. 0 ) goto 20

                  end do

               end if

                     goto 500

   20       continue

*-----------------------------------------------------------------------

            do im = 1, nm
            do ip = 1, ipn



!OBINATA(2012.7.11): Ct = Ct + wi
               trEVENT(ips(ip),ie,ia,it,ir,im) =
     &         trEVENT(ips(ip),ie,ia,it,ir,im) + wt1 * facm(im)


            end do
            end do

*-----------------------------------------------------------------------
*        dump data on file
*-----------------------------------------------------------------------

         if( itmdp(m,0) .ne. 0 ) then
               dmpd(1)  = dble( ktyp )
               dmpd(2)  = xc(ibkxc+no,ipomp+1) +
     &                          u(ibku+no,ipomp+1) * parz(28)
               dmpd(3)  = yc(ibkyc+no,ipomp+1) +
     &                          v(ibkv+no,ipomp+1) * parz(28)
               dmpd(4)  = zc(ibkzc+no,ipomp+1) +
     &                          w(ibkw+no,ipomp+1) * parz(28)
               dmpd(5)  = u(ibku+no,ipomp+1)
               dmpd(6)  = v(ibkv+no,ipomp+1)
               dmpd(7)  = w(ibkw+no,ipomp+1)
               dmpd(8)  = ec(ibkec+no,ipomp+1)
               dmpd(9)  = wt(ibkwt+no,ipomp+1)
               dmpd(10) = abs(tc(ibktc+no,ipomp+1))
               dmpd(11) = ncnt(ibknct+1,no,ipomp+1)
               dmpd(12) = ncnt(ibknct+2,no,ipomp+1)
               dmpd(13) = ncnt(ibknct+3,no,ipomp+1)
               dmpd(14) = spx(ibkspx+no,ipomp+1)
               dmpd(15) = spy(ibkspy+no,ipomp+1)
               dmpd(16) = spz(ibkspz+no,ipomp+1)
               dmpd(17) = name(ibknam+no,ipomp+1)
               dmpd(18) = nocas
               dmpd(19) = nobch
               dmpd(20) = no
               dmpd(21) = ctyp

               if( ityp .ge. 15 )
     &         dmpd(8) = dmpd(8) / ibryf(ityp,ktyp)

            if(idmpomp.ne.0)then       !FURUTA20150427
             call writeompfile(m,dmpd) !FURUTA20150427
            else                       !FURUTA20150427

             if( itmdp(m,0) .gt. 0 ) then

               write(itmdf(m))
     &         ( dmpd(itmdp(m,k)), k = 1, abs( itmdp(m,0) ) )

              else

               write(itmdf(m),'(30(1p1d24.15))')
     &         ( dmpd(itmdp(m,k)), k = 1, abs( itmdp(m,0) ) )

              end if

            endif !FURUTA20150427

         end if

*-----------------------------------------------------------------------

  500 continue

               j = j + mtrn

 1000 continue

*-----------------------------------------------------------------------
cc H.Iwase 2015/4/2 for electrons
cc                  change the continuous energy back to the step energy
cc T.Sato 2017/2/14, use "goto 999" to come here
 999  if( ityp.eq.12 .or. ityp.eq.13 )then

          e(ibke +no,ipomp+1) =  e(ibke +no,ipomp+1) +
     &                                     deinit  - denstepold
         ec(ibkec+no,ipomp+1) = ec(ibkec+no,ipomp+1) +
     &                                     deinit  - denstepnew


      endif


*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine psufreg(m,np,nr,ne,na,nt,nm,mr,kr,ar,eb,ab,tb,
     &                   tr,idasa)
*                                                                      *
*       output of region crossing surface tally of                     *
*       current spectrum and flux.                                     *
*       last modified by K.Niita on 2015/10/30                         *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall13/ itrcn(itlmax), itrcr(itlmax), itrca(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80
      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall21/ rtfac(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /istcut/ ist_cut, ist_bat

      common /fact01/ facmax(itlmax) ! kitamura23/03/31

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      common /talout/ itall

      character fname*100, fnume*3

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /mpi00/ npe, me

      common /tall50/ itlmt(itlmax), ite2l(itlmax)
*-----------------------------------------------------------------------

      dimension   kr(mr+1)  ! T.Sato 2022/10/28
      dimension   eb(ne+1)
      dimension   ew(ne)
      dimension   ab(na+1)
      dimension   tb(nt+1)
      dimension   aw(na)
      dimension   tw(nt)
      dimension   ar(nr)
      dimension   tr(np,ne,na,nt,nr,nm,2)

*-----------------------------------------------------------------------

      dimension tott(np,2) ! frtati 2021/10/05 6 -> np

      character hsunit(16)*32

      data hsunit( 1) / '[1/cm^2/source]                 '/
      data hsunit( 2) / '[1/cm^2/MeV/source]             '/
      data hsunit( 3) / '[1/cm^2/Lethargy/source]        '/
      data hsunit( 4) / '[1/cm^2/sr/source]              '/
      data hsunit( 5) / '[1/cm^2/MeV/sr/source]          '/
      data hsunit( 6) / '[1/cm^2/Lethargy/sr/source]     '/
      data hsunit(11) / '[1/cm^2/nsec/source]            '/
      data hsunit(12) / '[1/cm^2/MeV/nsec/source]        '/
      data hsunit(13) / '[1/cm^2/Lethargy/nsec/source]   '/
      data hsunit(14) / '[1/cm^2/sr/nsec/source]         '/
      data hsunit(15) / '[1/cm^2/MeV/sr/nsec/source]     '/
      data hsunit(16) / '[1/cm^2/Lethargy/sr/nsec/source]'/

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      character cha*1
      data cha /"'"/

      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31

      character rpa*1
      data rpa /'}'/

      character dum1*10000
      character dum2*10000
      character dum3*10000

      character cname*7
      character dname*7

      character aname*3

      include 'samepage_include/samepage000.inc'

      character yen*1
      yen  = char(92)

      include 'samepage_include/samepage001.inc'

*-----------------------------------------------------------------------

      if ( iMeVperu.eq. 1 ) then
         hsunit( 2) = '[1/cm^2/(MeV/n)/source]         '
         hsunit( 5) = '[1/cm^2/(MeV/n)/sr/source]      '
         hsunit(12) = '[1/cm^2/(MeV/n)/nsec/source]    '
         hsunit(15) = '[1/cm^2/(MeV/n)/sr/nsec/source] '
      end if

      if( ite2l(m) .eq. 1 ) then   ! convert to energy to LET
         hsunit( 2) = '[1/cm^2/(keV/um)/source]        '
         hsunit( 5) = '[1/cm^2/(keV/um)/sr/source]     '
         hsunit(12) = '[1/cm^2/(keV/um)/nsec/source]   '
         hsunit(15) = '[1/cm^2/(keV/um)/sr/nsec/source]'
      end if

*-----------------------------------------------------------------------
      np_mxang = np
      if ( itmxang(m).ne.0 .and. iabs(itmxang(m)).lt.np ) then
        np_mxang = iabs(itmxang(m))
      end if

*-----------------------------------------------------------------------
*        current or flux
*-----------------------------------------------------------------------

            if( itout(m) .eq. 1 ) then

               cname = '  Flux '
               dname = '  flux '

            else if( itout(m) .eq. 2 ) then

               cname = 'Current'
               dname = 'current'

            else if( itout(m) .eq. 3 ) then

               cname = ' F-Curr'
               dname = ' f-curr'

            else if( itout(m) .eq. 4 ) then

               cname = ' B-Curr'
               dname = ' b-curr'

            else if( itout(m) .eq. 5 ) then

               cname = ' O-Curr'
               dname = ' o-curr'

            else if( itout(m) .eq. 6 ) then

               cname = 'OF-Curr'
               dname = 'of-curr'

            else if( itout(m) .eq. 7 ) then

               cname = 'OB-Curr'
               dname = 'ob-curr'

            else if( itout(m) .eq. 8 ) then

               cname = ' A-Curr'
               dname = ' a-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            else if( itout(m) .eq. 9 ) then

               cname = 'OA-Curr'
               dname = 'oa-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            else if( itout(m) .eq. 10 ) then ! T.Sato 2021/05/05

               cname = ' A-Flux'
               dname = ' a-flux'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            else if( itout(m) .eq. 11 ) then ! T.Sato 2021/05/05

               cname = 'OA-Flux'
               dname = 'oa-flux'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               call wtpname(m,np,chp,chq)

*-----------------------------------------------------------------------
*        c1 : nomalization for source
*-----------------------------------------------------------------------

            if ( rsouin .gt. 0d0 ) then ! S.H. (2016.7.26)

               c1 = 1.0d+0 / rsouin

            else

               c1 = 0d0

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 1 : /cm^2/source
*                    = 2 : /cm^2/MeV/source
*                    = 3 : /cm^2/Lethargy/source
*                    = 4 : /cm^2/source/SR
*                    = 5 : /cm^2/MeV/source/SR
*                    = 6 : /cm^2/Lethargy/source/SR
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  4 .or.
     &          itunt(m) .eq. 11 .or. itunt(m) .eq. 14 ) then

               do i = 1, ne

                  ew(i) = 1.d+0

               end do
               ew_sum = 1.0d0

            else if( itunt(m) .eq.  2 .or. itunt(m) .eq.  5 .or.
     &               itunt(m) .eq. 12 .or. itunt(m) .eq. 15 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &               itunt(m) .eq. 13 .or. itunt(m) .eq. 16 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

            end if

*-----------------------------------------------------------------------
            aw_sum = 0.0d0

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  2 .or.
     &          itunt(m) .eq.  3 .or. itunt(m) .eq. 11 .or.
     &          itunt(m) .eq. 12 .or. itunt(m) .eq. 13 ) then

               do i = 1, na

                  aw(i) = 1.d+0

               end do
               aw_sum = 1.0d0

            else
               

               do i = 1, na

                  if( itaty(m) .gt. 0 ) then

                     aw(i) = ( ab(i+1) - ab(i) ) * 2.0 * pi

                  else

                     aw(i) = ( cos( ab(i) / 180.d0 * pi )
     &                       - cos( ab(i+1) / 180.d0 * pi ) ) * 2.0 * pi

                  end if
                  aw_sum = aw_sum + aw(i)

               end do

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 11, 12, 13, 14 : /.../nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .gt. 10 ) then

               do i = 1, nt

                  tw(i) = tb(i+1) - tb(i)

               end do
               tw_sum = tb(nt+1) - tb(1)

            else

               do i = 1, nt

                  tw(i) = 1.d+0

               end do
               tw_sum = 1.0d0

            end if

*-----------------------------------------------------------------------
*        relative error and  unit conversion
*-----------------------------------------------------------------------
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

            ar_sum = sum(ar(1:nr))

            facmax(m) = 1.d0
            if( rtfac(m) .lt. 0.d0 ) then
               facmax(m) = 0.d0

               do 101 im = 1, nm
               do 101 ir = 1, nr
               do 101 ia = 1, na
               do 101 it = 1, nt
               do 101 ie = 1, ne
               do 101 ip = 1, np

                  if( tr(ip,ie,ia,it,ir,im,1) .gt. 0.d0 ) then

                     fmaxfc = tr(ip,ie,ia,it,ir,im,1)
     &                               / ar(ir) / ew(ie) / aw(ia) / tw(it)

                     if( facmax(m) .lt. fmaxfc) facmax(m) = fmaxfc

                  end if

  101          continue

            end if

            do 100 im = 1, nm
            do 100 ir = 1, nr
            do 100 ia = 1, na
            do 100 it = 1, nt
            do 100 ie = 1, ne
            do 100 ip = 1, np

               if( tr(ip,ie,ia,it,ir,im,1) .gt. 0.d0 ) then

                  call calc_stdev(m,Xa,sigx,
     &                            tr(ip,ie,ia,it,ir,im,1),
     &                            tr(ip,ie,ia,it,ir,im,2),
     &             abs(rtfac(m)/facmax(m))/ar(ir)/ew(ie)/aw(ia)/tw(it))

                  tr(ip,ie,ia,it,ir,im,1) = Xa
                  tr(ip,ie,ia,it,ir,im,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

               else

                  isdz = 1
                  tr(ip,ie,ia,it,ir,im,2) = 0.0

               end if

! sumover
               fact_in = abs(rtfac(m)/facmax(m))
               call psufreg_sumover_stdev(0,m,ip,ie,ia,it,ir,im,
     &                fact_in,ew(ie),aw(ia),tw(it),ar(ir),
     &                ew_sum,aw_sum,tw_sum,ar_sum)


  100       continue

            if( nobch .gt. ist_bat ) then
            if( rtstd(m) .gt. 0.d0 .and.
     &          stdm .gt. 0.d0 .and. isdz .eq. 0 .and.
     &          stdm .lt. rtstd(m) ) itstd(m) = 0
            end if

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do iax = 1, itfln(m)

         if( itmdp(m,0) .eq. 0 ) then

            if( ( itall .eq. 2 .and. nobch .lt. maxbch )
     &           .or. itall .eq. 4 ) then

               call mk_2dnumfn(ctfln(m,iax),fname,itfll(m,iax),
     &              nobch,maxbch,npe)

            else

               fname = ctfln(m,iax)

            end if

         else

            if( ( itall .eq. 2 .and. nobch .lt. maxbch )
     &           .or. itall .eq. 4 ) then

               call mk_2dnumfn(ctfln(m,iax),fname,itfll(m,iax),
     &              nobch,maxbch,npe)

            else

               fname = ctfln(m,iax)(1:itfll(m,iax))

            end if

         end if

            iot = 31
            open(iot, file = fname, status = 'unknown' )



*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

               call tcrsech(iot,m,iax,1)

*-----------------------------------------------------------------------
*        energy axis  ( LET axis )
*-----------------------------------------------------------------------

      include 'samepage_include/samepage002_pmerat.inc'
      include 'samepage_include/samepagedum.inc'
      include 'samepage_include/samepagechp_pmerat.inc'
      include 'samepage_include/samepageseti.inc'

               inum = 0

         if( itaxs(m,iax) .eq. 1 .or.
     &       itaxs(m,iax) .eq. 14 ) then


            nmstepi = 1
            do imi = 1, nm
            do iri = 1, nr, nrstepi

             idum=(imi-1)*nr+iri


            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ir = iri
               ia = iai
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)



      l = 0
      dum3(l+1:l+9) = '#   no. ='
      l = l + 9
      write(dum3(l+1:l+5),'(i5)') inum
      l = l + 5
      if(iloopmode .ne. 0) then
        dum3(l+1:l+11) = '   '//chp(ipi)(3:10)
        l = l + 11
      endif
      if(iloopmode .ne. 2) then
        dum3(l+1:l+9) = '   reg = '
        l = l + 9
        dum3(l+1:l+lng1_n(idum)) = dum1_n(idum)(1:lng1_n(idum))
        l = l + lng1_n(idum)
        dum3(l+1:l+3) = ' - '
        l = l + 3
        dum3(l+1:l+lng2_n(idum)) = dum2_n(idum)(1:lng2_n(idum))
        l = l + lng2_n(idum)
      endif

      write(iot,'(a)') dum3(1:l)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,imi)

               if( itout(m) .gt. 7 ) then
                if(nastepi .eq. 1) then
                  write(iot,'(
     &            "#   ia  =",i3,/
     &            "#   ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+1)
                else
                  write(iot,'(
     &            "#   ia  =",i3," - ",i3/
     &            "#   ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, ia+nastepi-1,
     &                     aname, ab(ia), ab(ia+nastepi)
                endif
               end if

               if( ittty(m) .ne. 0 ) then
                if(ntstepi .eq. 1) then
                  write(iot,'(
     &            "#  it ="i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+1)
                else
                  write(iot,'(
     &            "#  it ="i3," - ",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, it+ntstepi-1,
     &                     tb(it), tb(it+ntstepi)
                endif
               end if

               if( itaxl(m) .eq. 0 ) then

                if ( itaxs(m,iax) .eq. 1 ) then   ! ccse 2022/09/30

                if ( iMeVperu.eq. 1 ) then
                  write(iot,'(/"x: Energy [MeV/n]")')
                else
                  write(iot,'(/"x: Energy [MeV]")')
                end if

                else if( itaxs(m,iax) .eq. 14 ) then
                 write(iot,'(/"x: LET [keV/um]")')
                end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                           cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 13 .or. itunt(m) .eq. 16 .or.
     &             itety(m) .eq.  3 .or. itety(m) .eq.  5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

       include 'samepage_include/peatrm_e.inc'

             if(iloopmode .eq. 0) then

                  l = 0
                  write(dum3(l+1:l+1),'(a1)') cha
                  l = l + 1
                  dum3(l+1:l+5) = 'no. ='
                  l = l + 5
                  write(dum3(l+1:l+5),'(i5)') inum
                  l = l + 5
                  dum3(l+1:l+9) = '   reg = '
                  l = l + 9
                  dum3(l+1:l+lng1_n(idum))=dum1_n(idum)(1:lng1_n(idum))
                  l = l + lng1_n(idum)
                  dum3(l+1:l+3) = ' - '
                  l = l + 3
                  dum3(l+1:l+lng2_n(idum))=dum2_n(idum)(1:lng2_n(idum))
                  l = l + lng2_n(idum)

               if( itout(m) .gt. 7 ) then
                  dum3(l+1:l+11) = ',   ang  = '
                  l = l + 11
                  write(dum3(l+1:l+3),'(i3)') ia
                  l = l + 3
               end if

               if( ittty(m) .ne. 0 ) then
                  dum3(l+1:l+9) = ',   t  = '
                  l = l + 9
                  write(dum3(l+1:l+3),'(i3)') it
                  l = l + 3
               end if

                  write(dum3(l+1:l+1),'(a1)') cha
                  l = l + 1

                  write(iot,'(600a1)') (dum3(i:i),i=1,l)

          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                changelsub(3) = " "

                changelsub(4) =
     &           '   reg = '//dum1_n(idum)(1:lng1_n(idum))//' - '//
     &                     dum2_n(idum)(1:lng2_n(idum))

                write(changelsub(5),'(",  ang =",i3)') ia
                write(changelsub(6),'(",  ie =",i3)') ie
                write(changelsub(7),'(",  it =",i3)') it
                write(changelsub(8),'(a1)') cha

               if( ittty(m) .eq. 0 ) then
                  changelsub(7) = " "
               end if

               if( itout(m) .lt. 8 ) then
                  changelsub(5) = " "
               end if
               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then
               else
                  changelsub(6) = " "
               end if

                changelsub(6) = " "

               if(iloopmode .eq. 2) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(7) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))
                write(iot,'(/a)') trim(angeltitle)
             end if

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

                  areasum =0.0d0
                  do i=iri,iri+nrstepi-1
                    areasum = areasum + ar(i)
                  end do
                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  region surface crossing"/
     &                     "  area  &=&",1pe13.4," [cm^2]")')
     &                     yen, areasum

               if( itout(m) .gt. 7 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        reg axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 2 ) then

             nmstepi = 1
            do imi = 1, nm, nmstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ie = iei
               ia = iai
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 ) then
                  if(nestepi .eq. 1) then
                    write(iot,'("#   no. =",i3,3x,
     &              "ie  =",i3/
     &              "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                          inum, ie,
     &                          eb(ie), eb(ie+nestepi)
                  else
                    write(iot,'("#   no. =",i3,3x,
     &              "ie  =",i3," - ",i3/
     &              "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                          inum, ie,ie+nestepi-1,
     &                          eb(ie), eb(ie+nestepi)
                  end if
               else if( itout(m) .eq. 8 .or.
     &                  itout(m) .eq.10 ) then


                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,"ia  =",i3,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iei, ia,
     &                        eb(ie), eb(ie+nestepi),
     &                        aname, ab(ia), ab(ia+nastepi)

               else if( itout(m) .le. 7 ) then

                  write(iot,'("#   no. =",i3)')
     &                        inum

               else if( itout(m) .eq. 9 .or.
     &                  itout(m) .eq. 11) then  ! T.Sato 2021/05/05

                  if(nastepi .eq. 1) then
                    write(iot,'("#   no. =",i3,3x,
     &              "ia  =",i3,/
     &              "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ia, aname,
     &                         ab(ia), ab(ia+nastepi)
                  else
                    write(iot,'("#   no. =",i3,3x,
     &              "ia  =",i3," - ",i3/
     &              "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ia,iai+nastepi-1,aname,
     &                         ab(ia), ab(ia+nastepi)
                 end if
               end if

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     tb(it), tb(it+ntstepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Serial Num. of Region")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname,hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

       include 'samepage_include/peatrm_r.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,",  e =",
     &                     1pe12.4,a1)')
     &                     cha, inum, eb(ie+1), cha
                  else
                     write(iot,'(/a1,"no. =",i3,",  e =",i3,
     &                     ",  t =",i3,a1)')
     &                     cha, inum, ie, it, cha
                  end if

               else if( itout(m) .le. 7 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,a1)')
     &                     cha, inum, cha
                  else
                     write(iot,'(/a1,"no. =",i3,",  t =",i3,a1)')
     &                     cha, inum, it, cha
                  end if

               else if( itout(m) .eq. 8 .or.
     &                  itout(m) .eq.10) then  ! T.Sato 2021/05/05

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,",  e =",
     &                     1pe12.4,",   ",a3," =",1pe12.4,a1)')
     &                     cha, inum, eb(ie+1), aname, ab(ia+1), cha
                  else
                     write(iot,'(/a1,"no. =",i3,",  e =",i3,
     &                     ",   ",a3," =",i3,",  t =",i3,a1)')
     &                     cha, inum, ie, aname, ia, it, cha
                  end if

               else if( itout(m) .eq. 9 .or.
     &                  itout(m) .eq.11) then ! T.Sato 2021/05/05

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,",   ",a3," =",
     &                     1pe12.4,a1)')
     &                     cha, inum, aname, ab(ia+1), cha
                  else
                     write(iot,'(/a1,"no. =",i3,",   ",a3," =",
     &                     i3,",  t =",i3,a1)')
     &                     cha, inum, aname, ia, it, cha
                  end if

               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,",  e =",
     &                     1pe12.4,",  mset ="i3,a1)')
     &                     cha, inum, eb(ie+1), itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,",  e =",i3,
     &                     ",  t =",i3,",  mset ="i3,a1)')
     &                     cha, inum, ie, it, itmnt(m,im), cha
                  end if

               else if( itout(m) .le. 7 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,",  mset ="i3,a1)')
     &                     cha, inum, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,",  t =",i3,
     &                               ",  mset ="i3,a1)')
     &                     cha, inum, it, itmnt(m,im), cha
                  end if

               else if( itout(m) .eq. 8 .or.
     &                  itout(m) .eq.10) then  ! T.Sato 2021/05/05

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,",  e =",
     &                     1pe12.4,",   ",a3," =",1pe12.4,
     &                     ",  mset ="i3,a1)')
     &                     cha, inum, eb(ie+1), aname, ab(ia+1),
     &                     itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,",  e =",i3,
     &                     ",   ",a3," =",i3,",  t =",i3,
     &                     ",  mset ="i3,a1)')
     &                     cha, inum, ie, aname, ia, it,
     &                     itmnt(m,im), cha
                  end if

               else if( itout(m) .eq. 9 .or.
     &                  itout(m) .eq.11) then ! T.Sato 2021/05/05


                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,",   ",a3," =",
     &                     1pe12.4,",  mset ="i3,a1)')
     &                     cha, inum, aname, ab(ia+1), itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,",   ",a3," =",
     &                     i3,",  t =",i3,",  mset ="i3,a1)')
     &                     cha, inum, aname, ia, it, itmnt(m,im), cha
                  end if

               end if

            end if

          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                changelsub(4) =
     &           '   reg = '//dum1_n(idum)(1:lng1_n(idum))//
     &                      ' - '//dum2_n(idum)(1:lng2_n(idum))

                write(changelsub(5),'(",  ang =",i3)') ia
                write(changelsub(6),'(",  ie =",i3)') ie
                write(changelsub(7),'(",  it =",i3)') it
                write(changelsub(8),'(a1)') cha

               if( ittty(m) .eq. 0 ) then
                  changelsub(7) = " "
               end if

               if( itout(m) .lt. 8 ) then
                  changelsub(5) = " "
               end if
               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then
               else
                  changelsub(6) = " "
               end if

                changelsub(4) = " "

               if(iloopmode .eq. 2) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(7) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  region surface crossing")') yen

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif


               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then  ! T.Sato 2021/05/05

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+nestepi)
                end if
               end if

               if( itout(m) .ge. 8 ) then

                     write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        angle axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 .or.
     &            itaxs(m,iax) .eq. 10 ) then

                  j = 0
                  igm = ( mmmax - 1 ) * 2 + 1

            nmstepi = 1
            do imi = 1, nm
            do iri = 1, nr, nrstepi

               idum = (imi-1)*nr + iri

            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ir = iri
               ie = iei
               it = iti
               ip = ipi


               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)



      l = 0
      dum3(l+1:l+9) = '#   no. ='
      l = l + 9
      write(dum3(l+1:l+5),'(i5)') inum
      l = l + 5
      if(iloopmode .ne. 0) then
        dum3(l+1:l+11) = '   '//chp(ipi)(3:10)
        l = l + 11
      endif
      if(iloopmode .ne. 2) then
        dum3(l+1:l+9) = '   reg = '
        l = l + 9
        dum3(l+1:l+lng1_n(idum)) = dum1_n(idum)(1:lng1_n(idum))
        l = l + lng1_n(idum)
        dum3(l+1:l+3) = ' - '
        l = l + 3
        dum3(l+1:l+lng2_n(idum)) = dum2_n(idum)(1:lng2_n(idum))
        l = l + lng2_n(idum)
      endif

      write(iot,'(a)') dum3(1:l)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

               if( itout(m) .eq. 8 .or.
     &                  itout(m) .eq.10) then  ! T.Sato 2021/05/05

                  write(iot,'("#   ie  =",i3,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ie, eb(ie), eb(ie+nestepi)

               end if

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     tb(it), tb(it+ntstepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  if( itaxs(m,iax) .eq. 8 ) then
                     write(iot,'(/"x: cos(",a1,"theta)")') yen
                  else
                     write(iot,'(/"x: ",a1,"theta  [deg]")') yen
                  end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

       include 'samepage_include/peatrm_a.inc'

             if(iloopmode .eq. 0) then

                  l = 0
                  write(dum3(l+1:l+1),'(a1)') cha
                  l = l + 1
                  dum3(l+1:l+5) = 'no. ='
                  l = l + 5
                  write(dum3(l+1:l+5),'(i5)') inum
                  l = l + 5
                  dum3(l+1:l+9) = '   reg = '
                  l = l + 9
                  dum3(l+1:l+lng1_n(idum)) =
     &                dum1_n(idum)(1:lng1_n(idum))
                  l = l + lng1_n(idum)
                  dum3(l+1:l+3) = ' - '
                  l = l + 3
                  dum3(l+1:l+lng2_n(idum)) =
     &                dum2_n(idum)(1:lng2_n(idum))
                  l = l + lng2_n(idum)

               if( itout(m) .eq. 8 .or.
     &                  itout(m) .eq.10) then  ! T.Sato 2021/05/05
                  dum3(l+1:l+10) = ',    e  = '
                  l = l + 10
                  write(dum3(l+1:l+3),'(i3)') ie
                  l = l + 3
               end if
               if( ittty(m) .ne. 0 ) then
                  dum3(l+1:l+9) = ',   t  = '
                  l = l + 9
                  write(dum3(l+1:l+3),'(i3)') it
                  l = l + 3
               end if

                  write(dum3(l+1:l+1),'(a1)') cha
                  l = l + 1

                  write(iot,'(600a1)') (dum3(i:i),i=1,l)

          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                changelsub(3) = " "

                changelsub(4) =
     &           '   reg = '//dum1_n(idum)(1:lng1_n(idum))//' - '//
     &                        dum2_n(idum)(1:lng2_n(idum))

                write(changelsub(5),'(",  ang =",i3)') ia
                write(changelsub(6),'(",  ie =",i3)') ie
                write(changelsub(7),'(",  it =",i3)') it
                write(changelsub(8),'(a1)') cha

               if( ittty(m) .eq. 0 ) then
                  changelsub(7) = " "
               end if

               if( itout(m) .lt. 8 ) then
                  changelsub(5) = " "
               end if
               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then
               else
                  changelsub(6) = " "
               end if

                changelsub(5) = " "

               if(iloopmode .eq. 2) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(7) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))
                write(iot,'(/a)') trim(angeltitle)
             end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

                  areasum = 0.0d0
                  do i = iri,iri+nrstepi-1
                    areasum = areasum + ar(i)
                  end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  region surface crossing"/
     &                     "  area  &=&",1pe13.4," [cm^2]")')
     &                     yen, areasum

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then  ! T.Sato 2021/05/05

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(iei), eb(iei+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+nestepi)
                end if
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        time axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9 ) then

                  j = 0
                  igm = ( mmmax - 1 ) * 2 + 1

            nmstepi = 1
            do imi = 1, nm
            do iri = 1, nr, nrstepi

              idum = (imi - 1)*nr + iri

            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
               im = imi
               ir = iri
               ie = iei
               ia = iai
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)



      l = 0
      dum3(l+1:l+9) = '#   no. ='
      l = l + 9
      write(dum3(l+1:l+5),'(i5)') inum
      l = l + 5
      if(iloopmode .ne. 0) then
        dum3(l+1:l+11) = '   '//chp(ipi)(3:10)
        l = l + 11
      endif
      if(iloopmode .ne. 2) then
        dum3(l+1:l+9) = '   reg = '
        l = l + 9
        dum3(l+1:l+lng1_n(idum)) = dum1_n(idum)(1:lng1_n(idum))
        l = l + lng1_n(idum)
        dum3(l+1:l+3) = ' - '
        l = l + 3
        dum3(l+1:l+lng2_n(idum)) = dum2_n(idum)(1:lng2_n(idum))
        l = l + lng2_n(idum)
      endif

      write(iot,'(a)') dum3(1:l)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then  ! T.Sato 2021/05/05

                  if(nestepi .eq. 1) then
                    write(iot,'("#   ie  =",i3,/
     &              "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                       ie, eb(ie), eb(ie+1)
                  else
                    write(iot,'("#   ie  =",i3," - ",i3/
     &              "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                       ie, ie+nestepi, eb(ie), eb(ie+nestepi)
                  endif

               end if

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     aname, ab(ia), ab(ia+nastepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Time [nsec]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( ittty(m) .eq. 3 .or. ittty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

       include 'samepage_include/peatrm_t.inc'

             if(iloopmode .eq. 0) then

                  l = 0
                  write(dum3(l+1:l+1),'(a1)') cha
                  l = l + 1
                  dum3(l+1:l+5) = 'no. ='
                  l = l + 5
                  write(dum3(l+1:l+5),'(i5)') inum
                  l = l + 5
                  dum3(l+1:l+9) = '   reg = '
                  l = l + 9
                  dum3(l+1:l+lng1_n(idum))=dum1_n(idum)(1:lng1_n(idum))
                  l = l + lng1_n(idum)
                  dum3(l+1:l+3) = ' - '
                  l = l + 3
                  dum3(l+1:l+lng2_n(idum))=dum2_n(idum)(1:lng2_n(idum))
                  l = l + lng2_n(idum)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then  ! T.Sato 2021/05/05

                  dum3(l+1:l+10) = ',    e  = '
                  l = l + 10
                  write(dum3(l+1:l+3),'(i3)') ie
                  l = l + 3
               end if
               if( itout(m) .ge. 8 ) then
                  dum3(l+1:l+10) = ',  ang  = '
                  l = l + 10
                  write(dum3(l+1:l+3),'(i3)') ia
                  l = l + 3
               end if

                  write(dum3(l+1:l+1),'(a1)') cha
                  l = l + 1

                  write(iot,'(600a1)') (dum3(i:i),i=1,l)

          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                changelsub(3) = " "

                changelsub(4) =
     &           '   reg = '//dum1_n(idum)(1:lng1_n(idum))//' - '//
     &                        dum2_n(idum)(1:lng2_n(idum))

                write(changelsub(5),'(",  ang =",i3)') ia
                write(changelsub(6),'(",  ie =",i3)') ie
                write(changelsub(7),'(",  it =",i3)') it
                write(changelsub(8),'(a1)') cha

               if( ittty(m) .eq. 0 ) then
                  changelsub(7) = " "
               end if

               if( itout(m) .lt. 8 ) then
                  changelsub(5) = " "
               end if
               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then
               else
                  changelsub(6) = " "
               end if

                changelsub(7) = " "

               if(iloopmode .eq. 2) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(7) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))
                write(iot,'(/a)') trim(angeltitle)
             end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

                  areasum = 0.0d0
                  do i=iri,iri+nrstepi-1
                    areasum = areasum + ar(i)
                  end do
                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  region surface crossing"/
     &                     "  area  &=&",1pe13.4," [cm^2]")')
     &                     yen, areasum

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if
               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then  ! T.Sato 2021/05/05


                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+nestepi)
                end if
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do

         end if

*-----------------------------------------------------------------------

            call prestart(m,iot) !OBINATA(2012.6.13)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

      end do

*-----------------------------------------------------------------------

      include 'samepage_include/samepage999.inc'

      return
      end


************************************************************************
*                                                                      *
      subroutine tsufrz(ncol,m,np,nr,nz,ne,na,nt,nm,
     &                  rm,zm,eb,ab,tb,tr,tz,
     &                  trEVENT,tzEVENT,
     &                  itrmax,itrmin,itzmax,itzmin)
*                                                                      *
*       surface crossing tally of current spectrum and flux            *
*       for r-z scoring mesh                                           *
*       last modified by K.Niita on 2015/10/30                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*             10 : geometry boundary crossing                          *
*             11 : termination by energy cut-off                       *
*             12 : termination by escape or leakage                    *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
*             15 : sequential transport only for tally                 *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use partmod, only: itmxpt, itpan, itpat, jtpat ! frtati 2021/10/05
*-----------------------------------------------------------------------

      implicit double precision( a-h, o-z )

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter( rlit = 29.97925d0 )
      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)

      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)

*-----------------------------------------------------------------------

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall08/ rtrx0(itlmax), rtry0(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)
      common /tall73/ itenclo(itlmax), itangform(itlmax)

      common /tall82/ itcnth(9,itlmax)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)

      common /tall50/ itlmt(itlmax), ite2l(itlmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

*-----------------------------------------------------------------------

      dimension   rm(nr+1)
      dimension   zm(nz+1)
      dimension   eb(ne+1)
      dimension   ab(na+1)
      dimension   tb(nt+1)

      dimension   tr(np,ne,na,nt,(nr+1)*nz,nm,2)

      dimension   tz(np,ne,na,nt,nr*(nz+1),nm,2)

      dimension   trEVENT(np,ne,na,nt,(nr+1)*nz,nm)    !OBINATA(2012.7.11): as Ct

      dimension   tzEVENT(np,ne,na,nt,nr*(nz+1),nm)    !OBINATA(2012.7.11): as Ct

      real(8),allocatable,save:: tr0(:,:,:,:,:,:) !OBINATA(2012.7.11): as C
      real(8),allocatable,save:: tz0(:,:,:,:,:,:) !OBINATA(2012.7.11): as C

      dimension   itrmax(7),itrmin(7)
      dimension   itzmax(7),itzmin(7)

*-----------------------------------------------------------------------

      dimension ud(3)

      data dmax  /1.0d+19/
      data dmax0 /1.0d+18/

*-----------------------------------------------------------------------

      data grab,hgrab / 0.00349d+0, 0.001745d+0 /

      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt

*-----------------------------------------------------------------------

      dimension facm(6)

*-----------------------------------------------------------------------

      common /stat / istdev, irestart, ireschk
      common /cparm/ maxbch,maxcas

*-----------------------------------------------------------------------

      common /paraj/ mstz(300), parz(300)

      common /egs5cmn10/denstepold,denstepnew,deinit
      real*8            denstepold,denstepnew,deinit
!$OMP THREADPRIVATE(/egs5cmn10/)

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

*------------------------------------------------------------------------

      common /spred/ nspred, nwsprd, nedisp, itstep, ndedx


*-----------------------------------------------------------------------

         irf(ir,iz) = ir + ( iz - 1 ) * ( nr + 1 )
         izf(ir,iz) = iz + ( ir - 1 ) * ( nz + 1 )

*-----------------------------------------------------------------------

      if (istdev .eq. 2) then

        call readitrminmax7(itrmin,itrmax,(/ np,ne,na,nt,nr+1,nz,nm/),
     &                      mnpr,mner,mnar,mntr,mnrr,mnzr,mnmr,
     &                      mxpr,mxer,mxar,mxtr,mxrr,mxzr,mxmr)
        call readitrminmax7(itzmin,itzmax,(/ np,ne,na,nt,nr,nz+1,nm/),
     &                      mnpz,mnez,mnaz,mntz,mnrz,mnzz,mnmz,
     &                      mxpz,mxez,mxaz,mxtz,mxrz,mxzz,mxmz)


      endif

*-----------------------------------------------------------------------
* check of history counter
*-----------------------------------------------------------------------
      ihistcount = 0 ! history counter index, T.Sato 2022/12/20 for speed up
         if ( ncol.eq.0 .or. ncol.eq.4 ) then
           do i = 1, 3
             if( itcnth(i,m) .eq. 1 ) then
               if( ncntmx(i) .lt. itcnth(i*2+2,m) .or.
     &             ncntmx(i) .gt. itcnth(i*2+3,m) ) then
                 ihistcount = 1
               end if
             end if
           end do
         end if

*-----------------------------------------------------------------------
*        source particle or end of batch ( in case of istdev = 2 )
*
* OBINATA(2012.7.11): change tr(,,,,,3) to trEVENT(,,,,)
*-----------------------------------------------------------------------

         if (( ncol .eq. 0 .or. ncol .eq. 4 )
     &                              .and. istdev .eq. 2) then
           if ((nocas.gt.1.or.ncol.eq.0) .and. ihistcount.ne.1 ) then

             do im = mnmr,mxmr
             do iz = mnzr,mxzr
             do ir = mnrr,mxrr
             do ie = mner,mxer
               irz = irf(ir,iz)
               tr(:,ie,:,:,irz,im,1) = tr(:,ie,:,:,irz,im,1)
     &                              + trEVENT(:,ie,:,:,irz,im)
               tr(:,ie,:,:,irz,im,2) = tr(:,ie,:,:,irz,im,2)
     &                              + trEVENT(:,ie,:,:,irz,im) ** 2
             enddo
             enddo
             enddo
             enddo

             do im = mnmz,mxmz
             do iz = mnzz,mxzz
             do ir = mnrz,mxrz
             do ie = mnez,mxez
               irz = izf(ir,iz)
               tz(:,ie,:,:,irz,im,1) = tz(:,ie,:,:,irz,im,1)
     &                              + tzEVENT(:,ie,:,:,irz,im)
               tz(:,ie,:,:,irz,im,2) = tz(:,ie,:,:,irz,im,2)
     &                              + tzEVENT(:,ie,:,:,irz,im) ** 2
             enddo
             enddo
             enddo
             enddo


             call tsufrz_sumover(m,1,
     &                   np, ne, na, nt, nr, nz,  nm, trEVENT, tzEVENT)


           end if

           do im = mnmr,mxmr
           do iz = mnzr,mxzr
           do ir = mnrr,mxrr
           do ie = mner,mxer
             irz = irf(ir,iz)
             trEVENT(:,ie,:,:,irz,im) = 0
           enddo
           enddo
           enddo
           enddo

           do im = mnmz,mxmz
           do iz = mnzz,mxzz
           do ir = mnrz,mxrz
           do ie = mnez,mxez
             irz = izf(ir,iz)
             tzEVENT(:,ie,:,:,irz,im) = 0
           enddo
           enddo
           enddo
           enddo

           call resetitrminmax
     &          (itrmin,itrmax,7,(/np,ne,na,nt,nr+1,nz,nm/))
           call resetitrminmax
     &          (itzmin,itzmax,7,(/np,ne,na,nt,nr,nz+1,nm/))

         end if

*-----------------------------------------------------------------------
*        end of batch ( in case of istdev = 1 )
*
* OBINATA(2012.5.29): modificate for thread parallel
*-----------------------------------------------------------------------

         if ( ncol .eq. 0 .and. istdev .eq. 1) then
!$OMP MASTER
             allocate( tr0(np,ne,na,nt,(nr+1)*nz,nm) )
             allocate( tz0(np,ne,na,nt,nr*(nz+1),nm) )


             tr0(:,:,:,:,:,:) = 0.d0
             tz0(:,:,:,:,:,:) = 0.d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP CRITICAL (tsufrz_crit_ist1_r)
             tr0(:,:,:,:,:,:) = tr0(:,:,:,:,:,:) + trEVENT(:,:,:,:,:,:)
!$OMP END CRITICAL (tsufrz_crit_ist1_r)
!$OMP CRITICAL (tsufrz_crit_ist1_z)
             tz0(:,:,:,:,:,:) = tz0(:,:,:,:,:,:) + tzEVENT(:,:,:,:,:,:)
!$OMP END CRITICAL (tsufrz_crit_ist1_z)
!$OMP BARRIER
!$OMP MASTER
             tr(:,:,:,:,:,:,1) = tr(:,:,:,:,:,:,1)
     &                         + tr0(:,:,:,:,:,:) / maxcas
             tr(:,:,:,:,:,:,2) = tr(:,:,:,:,:,:,2)
     &                         + ( tr0(:,:,:,:,:,:) / maxcas ) ** 2
             tz(:,:,:,:,:,:,1) = tz(:,:,:,:,:,:,1)
     &                         + tz0(:,:,:,:,:,:) / maxcas
             tz(:,:,:,:,:,:,2) = tz(:,:,:,:,:,:,2)
     &                         + ( tz0(:,:,:,:,:,:) / maxcas ) ** 2

             call tsufrz_sumover(m,maxcas,
     &                   np, ne, na, nt, nr, nz,  nm, tr0, tz0)

             deallocate( tr0 )
             deallocate( tz0 )
!$OMP END MASTER

           trEVENT(:,:,:,:,:,:) = 0
           tzEVENT(:,:,:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------

         small = parz(28) * 10.d0

*-----------------------------------------------------------------------
*        check of ncol
*-----------------------------------------------------------------------

         if( ncol .lt. 10 ) return

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncnt(ibknct+i,no,ipomp+1) .lt. itcnt(i*2+2,m) .or.
     &                ncnt(ibknct+i,no,ipomp+1) .gt. itcnt(i*2+3,m) )
     &                return

               end if

            end do

*-----------------------------------------------------------------------
*        material for LET
*-----------------------------------------------------------------------

               if( itlmt(m) .gt. 0 ) then

                  lmat = idnm( itlmt(m) )

               else if( itlmt(m) .eq. 0 ) then

                  lmat = mat

               else

                  lmat = -idnm( -itlmt(m) )

               end if

*-----------------------------------------------------------------------
*        check of particles
*-----------------------------------------------------------------------

            call pcheck(m,np,ityp,ktyp,jtyp,ipn,ips)

               if( ipn .eq. 0 ) return

*-----------------------------------------------------------------------
cc H.Iwase 2015/4/2 for electrons
cc                  change the step energy to coutinuous energy
*-----------------------------------------------------------------------

      if( ityp.eq.12 .or. ityp.eq.13 )then

         e1 =  e(ibke  +no,ipomp+1) - deinit  + denstepold
         e2 = ec(ibkec +no,ipomp+1) - deinit  + denstepnew

        if (ite2l(m) .eq. 0) then   ! ccse 2022/08/31

         if( e1 .lt. eb(1)    ) return
         if( e2 .ge. eb(ne+1) ) return

        else if (ite2l(m) .eq. 1) then   ! energy to LTE
            call dedxas( e1,dedx,lmat,
     &                  ityp,ktyp,jtyp,rtyp)
            dedxl = dedx /10.0d0
            call dedxas( e2,dedx,lmat,
     &                  ityp,ktyp,jtyp,rtyp)
            dedxh = dedx /10.0d0
            if ( dedxl .gt. dedxh ) then   ! hight low swap
               dedx  = dedxh
               dedxh = dedxl
               dedxl = dedx
            end if

            if( dedxl .lt. eb(1)    ) return
            if( dedxh .ge. eb(ne+1) ) return
        end if

         tparti = abs(t(ibkt+no,ipomp+1))
         tpartf = abs(tc(ibktc+no,ipomp+1))
         if( tparti .ge. tb(nt+1) ) return
         if( tpartf .lt. tb(1) ) return

          e(ibke +no,ipomp+1) = e1
         ec(ibkec+no,ipomp+1) = e2

      endif

*-----------------------------------------------------------------------

      if( iMeVperu.eq.1 .and. ityp.ge.15 .and. ityp.le.19 ) then
           ebm = ktyp - ktyp / 1000000 * 1000000
      else
           ebm = 1.d0
      end if

*-----------------------------------------------------------------------
*        transform positions
*-----------------------------------------------------------------------

            call trnsxx(x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  xxa,yya,zza,itmtr(m,4))

            call trnsxx(xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                  zc(ibkzc+no,ipomp+1),
     &                  xxc,yyc,zzc,itmtr(m,4))

*-----------------------------------------------------------------------
*        check z mesh and r mesh region
*-----------------------------------------------------------------------

               x0 = rtrx0(m)
               y0 = rtry0(m)

            if( zza+small .lt. zm(1) .and.
     &          zzc+small .lt. zm(1) ) goto 999

            if( zza-small .ge. zm(nz+1) .and.
     &          zzc-small .ge. zm(nz+1) ) goto 999

               dis0 = sqrt( ( xxa - x0 )**2
     &                    + ( yya - y0 )**2 )

               dis1 = sqrt( ( xxc - x0 )**2
     &                    + ( yyc - y0 )**2 )

            if( dis0 .lt. rm(1) .and.
     &          dis1 .lt. rm(1) ) goto 999

            if( dis0 .ge. rm(nr+1) .and.
     &          dis1 .ge. rm(nr+1) ) then

               aa =   yyc - yya
               bb = - xxc + xxa
               cc = - aa * xxa - bb * yya

               if( aa**2 + bb**2 .ne. 0.0d0 ) then

                  dd = abs( aa * x0 + bb * y0 + cc )
     &               / sqrt( aa**2 + bb**2 )

                  if( dd .ge. rm(nr+1) ) goto 999

               end if

            end if

*-----------------------------------------------------------------------
*        check of energy
*-----------------------------------------------------------------------

         if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &       itout(m) .eq.10) then  ! T.Sato 2021/05/05

         if (ite2l(m) .eq. 0) then   ! ccse 2022/08/31
            if( e(ibke+no,ipomp+1)   .lt. eb(1)*ebm ) goto 999
            if( ec(ibkec+no,ipomp+1) .ge. eb(ne + 1)*ebm ) goto 999

         else if (ite2l(m) .eq. 1) then   ! convert to energy to LET
            call dedxas( e(ibke+no,ipomp+1),dedx,lmat,
     &                  ityp,ktyp,jtyp,rtyp)
            dedxl = dedx /10.0d0
            call dedxas(ec(ibkec+no,ipomp+1),dedx,lmat,
     &                  ityp,ktyp,jtyp,rtyp)
            dedxh = dedx /10.0d0
            if ( dedxl .gt. dedxh ) then   ! hight low swap
               dedx  = dedxh
               dedxh = dedxl
               dedxl = dedx
            end if
            if( dedxl .lt. eb(1) ) goto 999
            if( dedxh .ge. eb(ne + 1) ) goto 999
         end if

         end if

*-----------------------------------------------------------------------
*        check of time
*-----------------------------------------------------------------------

               tpart = abs(t(ibkt+no,ipomp+1))

               if( tpart  .ge. tb(nt+1) ) goto 999
               if( abs(tc(ibktc+no,ipomp+1)) .lt. tb(1) ) goto 999

*-----------------------------------------------------------------------

            icli = idgr(iblz1)

*-----------------------------------------------------------------------
*        distance and unit vector ud(i)
*-----------------------------------------------------------------------

            dis = ( xxc - xxa )**2
     &          + ( yyc - yya )**2
     &          + ( zzc - zza )**2

            if( dis .le. small**2 ) goto 999

            dis = sqrt(dis)

            ud(1) = ( xxc - xxa ) / dis
            ud(2) = ( yyc - yya ) / dis
            ud(3) = ( zzc - zza ) / dis

*-----------------------------------------------------------------------
*     initial position, energy and initial range rng
*-----------------------------------------------------------------------

            tot = 0.0d0

            xpp = xxa
            ypp = yya
            zpp = zza

            se  = e(ibke+no,ipomp+1)
            ee  = ec(ibkec+no,ipomp+1)

            if( jtyp .ne. 0 .and. mat .gt. 0 )
     &          call rainge(se,rng,mat,ityp,ktyp,jtyp,rtyp)

*-----------------------------------------------------------------------
*     initial z-position
*-----------------------------------------------------------------------

               izm = 0
               izc = 0

            if( ud(3) .gt. 0.0d0 ) izk =  1
            if( ud(3) .eq. 0.0d0 ) izk =  0
            if( ud(3) .lt. 0.0d0 ) izk = -1

            if( izk .ge. 0 ) then

               do i = 1, nz + 1
                  if( zm(i) .ge. zpp - small ) then
                         izm = i - 1
                         izc = i - 1
                         goto 38
                  end if
               end do

            else

               izm = nz + 2
               izc = nz + 1

               do i = 1, nz + 1
                  if( zm(i) .gt. zpp + small ) then
                         izm = i
                         izc = i - 1
                         goto 38
                  end if
               end do

            end if

   38       continue

*-----------------------------------------------------------------------
*     initial r-position
*-----------------------------------------------------------------------

               uvec = ( xxa - x0 ) * ud(1) + ( yya - y0 ) * ud(2)

            if( uvec .ge. 0.0d0 ) then

               do i = 1, nr + 1
                  if( rm(i) .gt. dis0 + small ) then
                        irc = i - 1
                        goto 47
                  end if
               end do

            else

               do i = 1, nr + 1
                  if( rm(i) .ge. dis0 - small ) then
                        irc = i - 1
                        goto 47
                  end if
               end do

            end if

                  irc = nr + 1

   47          continue

*-----------------------------------------------------------------------
*     loop for finding mesh and booking upto total distance
*-----------------------------------------------------------------------

   50 continue

*-----------------------------------------------------------------------
*      calculation is finished
*-----------------------------------------------------------------------

         if( tot .ge. dis ) goto 999

*-----------------------------------------------------------------------
*        check of energy
*-----------------------------------------------------------------------

         if (ite2l(m) .eq. 0) then   ! ccse 2022/08/31

         if( ( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &         itout(m) .eq.10 ) .and.  ! T.Sato 2021/05/05
     &       se .lt. eb(1)*ebm ) goto 999
         if( ( ( itout(m) .ge. 5 .and. itout(m) .le. 7 ) .or.
     &           itout(m) .eq. 9 .or. itout(m) .eq. 11) .and.
     &       se .lt. 0.0d0 ) goto 999

         else if (ite2l(m) .eq. 1) then   ! convert to energy to LET
            call dedxas( se,dedx,lmat,ityp,ktyp,jtyp,rtyp)
            dedx = dedx /10.0d0
            if( ( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &            itout(m) .eq.10 ) .and.
     &            dedx .lt. eb(1) ) goto 999
            if( ( ( itout(m) .ge. 5 .and. itout(m) .le. 7 ) .or.
     &              itout(m) .eq. 9 .or. itout(m) .eq. 11) .and.
     &          dedx .lt. 0.0d0 ) goto 999
         end if

*-----------------------------------------------------------------------
*        dr : distance to the nearest r mesh
*-----------------------------------------------------------------------

               dr  = dmax
               idr = 0

               a  = ud(1)**2 + ud(2)**2

         if( a .gt. 0.0d0 ) then

               b = 2.0d0 * ( ud(1) * ( xpp - x0 )
     &                     + ud(2) * ( ypp - y0 ) )

               ir1 = irc
               ir2 = irc + 1

              if( b .ge. 0.0 ) ir1 = irc + 1

            do 51 i = max( 1, ir1 ), min( nr + 1, ir2 )

               if( rm(i) .eq. 0.0d0 ) goto 51

               c = ( xpp - x0 )**2
     &           + ( ypp - y0 )**2 - rm(i)**2

               dchk = b**2 - 4.0d0 * a * c

               if( dchk .lt. 0.0d0 ) goto 51

               t1 = ( -b + sqrt(dchk) ) / ( a * 2.0 )
               t2 = ( -b - sqrt(dchk) ) / ( a * 2.0 )

               if( t1 .lt. small ) t1 = dmax
               if( t2 .lt. small ) t2 = dmax

               ttr = min(t1,t2)

               if( ttr .lt. dr ) then

                  dr  = ttr
                  idr = i

               end if

   51       continue

         end if

*-----------------------------------------------------------------------
*        dz : distance to the nearest z mesh
*-----------------------------------------------------------------------

               dz = dmax

         if( izk .ne. 0 ) then
            if( izm + izk .ge. 1 .and.
     &          izm + izk .le. nz + 1 ) then

               dz = abs( ( zm(izm+izk) - zpp ) / ud(3) )

            end if
         end if

*-----------------------------------------------------------------------
*        which boundary is the nearlist
*-----------------------------------------------------------------------

            if( dr .gt. dz ) then

               dd = dz
               jz = 1

            else

               dd = dr
               jz = 0

            end if

            if( dd .gt. dmax0 ) goto 999

*-----------------------------------------------------------------------
*        recalculate rng for ncol = 11 and nedisp
*-----------------------------------------------------------------------

            if( mat .gt. 0 .and. jtyp .ne. 0 .and. nedisp .ne. 0 .and.
     &           (  ityp .ne. 12 .and. ityp .ne. 13 ) ) then


               call rainge(se,rng,mat,ityp,ktyp,jtyp,rtyp)

               if( ncol .eq. 11 )  dis = tot + rng

            end if

*-----------------------------------------------------------------------
*        propagate position upto the boundary or the final point
*-----------------------------------------------------------------------

            if( tot + dd .ge. dis - small ) then

               if( jz .eq. 1 ) goto 999
               if( tot + dd .ge. dis + small ) goto 999

               dd  = dis - tot
               tot = dis
               xpp = xxc
               ypp = yyc
               zpp = zzc

            else

               tot = tot + dd
               xpp = xpp + dd * ud(1)
               ypp = ypp + dd * ud(2)
               zpp = zpp + dd * ud(3)

            end if

*-----------------------------------------------------------------------
*        final energy ee and range rng
*-----------------------------------------------------------------------

            if( jtyp .ne. 0 .and. mat .gt. 0 ) then

               delt = dd

               call ecol(ee,delt,se,rng,
     &                   mat,ityp,ktyp,jtyp,rtyp)

               rng = max( 0.0d0, rng - delt )
               ee  = max( 0.0d0, ee )
               ee  = min( se, ee )

               if (ite2l(m) .eq. 1) then   ! ccse 2022/08/31
                  call dedxas( ee,dedx,lmat,ityp,ktyp,jtyp,rtyp)
                  dedx = dedx /10.0d0
               end if

            end if

*-----------------------------------------------------------------------
*        time evolution
*-----------------------------------------------------------------------

            if( t(ibkt+no,ipomp+1) .gt. 0.d0 ) then

               ekin = ( se + ee ) / 2.0
               dist = dd
               timd = 0.0

               if( ekin .gt. 0.1 .or. rtyp .eq. 0.0d0 ) then

                  timd = dist * ( ekin + rtyp )
     &                 / sqrt( ekin * ( ekin + 2.0 * rtyp ) )
     &                 / rlit

               else if( ekin .gt. 0.0 ) then

                  timd = dist * sqrt( rtyp / 2.0 / ekin ) / rlit

               end if

                  tpart = tpart + timd

            end if

*-----------------------------------------------------------------------
*     booking
*     for z-crossing jz = 1, r-crossing jz = 0
*-----------------------------------------------------------------------

            if( jz .eq. 1 ) then

               iz = izm + izk
               ir = irc

            else if( jz .eq. 0 ) then

               iz = izc
               ir = idr

            end if

*-----------------------------------------------------------------------
*        booking the crossing
*-----------------------------------------------------------------------

         if( ( ( itenclo(m) .eq. 1 ) .and.
     &         ( ir .ge. 1 .and. ir .le. nr + 1 .and.
     &           iz .ge. 1 .and. iz .le. nz + 1 ) )
     &       .or.
     &       ( ( itenclo(m) .ne. 1 .and. jz .eq. 1 ) .and.
     &         ( ir .ge. 1 .and. ir .lt. nr + 1 .and.
     &           iz .ge. 1 .and. iz .le. nz + 1 ) )
     &       .or.
     &       ( ( itenclo(m) .ne. 1 .and. jz .eq. 0 ) .and.
     &         ( ir .ge. 1 .and. ir .le. nr + 1 .and.
     &           iz .ge. 1 .and. iz .lt. nz + 1 ) ) ) then

*-----------------------------------------------------------------------
*        final energy cell
*-----------------------------------------------------------------------

            if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &          itout(m) .eq.10) then  ! T.Sato 2021/05/05

                     ie = 0

               if (ite2l(m) .eq. 0) then   ! ccse 2022/08/31
               do i = ne, 1, -1

                  if( ee .ge. eb(i)*ebm .and.
     &                ee .lt. eb(i+1)*ebm ) then

                     ie = i

                     goto 40

                  end if

               end do
               else if (ite2l(m) .eq. 1) then   ! convert to energy to LET
               do i = ne, 1, -1
                  if( dedx .ge. eb(i) .and.
     &                dedx .lt. eb(i+1) ) then
                     ie = i
                     goto 40
                  end if
               end do
               end if

   40          continue

*-----------------------------------------------------------------------
                     call fmfac(m,icli,ee,oldwt,facm)

*-----------------------------------------------------------------------
*        omni current case
*-----------------------------------------------------------------------

            else

               ie = 1

               facm(:) = 1.0d0

            end if

*-----------------------------------------------------------------------
*        tally
*-----------------------------------------------------------------------

            if( ie .gt. 0 ) then

                     ia = 1

*-----------------------------------------------------------------------
*              time
*-----------------------------------------------------------------------

                  do i = 1, nt

                     if( tpart .ge. tb(i) .and.
     &                   tpart .lt. tb(i+1) ) goto 60

                  end do

                     goto 61

   60                it = i

*-----------------------------------------------------------------------
*     identify ir, iz for enclosure mode
*        nrpt =  0 : normal mode
*             = -1 : save only backward
*             =  1 : save only forward
*             =  2 : save both direction, backward id is saved as ixrpt
*-----------------------------------------------------------------------

               nrpt = 0

               if( itenclo(m) .eq. 1 ) then

                  nrpt = 1
                  tlwrpt = 0.d0

                  if( jz . eq. 0 ) then
                     smu = 1.0
                     dxy = ( xpp - x0 ) * ud(1)
     &                   + ( ypp - y0 ) * ud(2)
                     rxy = ( xpp - x0 )**2
     &                   + ( ypp - y0 )**2
                     if( rxy .gt. 0.0 ) smu = dxy / sqrt(rxy)
                  endif

                  if( itout(m) .eq. 1 .or. itout(m) .eq. 2 .or.
     &                itout(m) .eq. 5 .or. itout(m) .ge. 8) then ! T.Sato 2021/05/05

                     nrpt = 2
                     iarpt = ia
                     irrpt = ir
                     izrpt = iz

                     if( jz .eq. 1 ) then
                        iz = izc + izk
                        izrpt = izc
                     elseif( jz .eq. 0 ) then
                        if( smu .ge. 0.d0 ) then
                           irrpt = ir - 1
                        else
                           ir = ir - 1
                        endif
                     endif

                     if( irrpt .lt. 1 .or. irrpt .ge. nr + 1 .or.
     &                   izrpt .lt. 1 .or. izrpt .ge. nz + 1 ) nrpt = 1

                  elseif( itout(m) .eq. 3 .or. itout(m) .eq. 6 ) then

                     if( jz .eq. 1 ) then
                        iz = izc + izk
                     elseif( jz .eq. 0 ) then
                        if( smu .ge. 0.d0 ) then
                           ir = ir
                        else
                           ir = ir - 1
                        endif
                     endif

                  elseif( itout(m) .eq. 4 .or. itout(m) .eq. 7 ) then

                     if( jz .eq. 1 ) then
                        iz  = izc
                     elseif( jz .eq. 0 ) then
                        if( smu .ge. 0.d0 ) then
                           ir = ir - 1
                        else
                           ir = ir
                        endif
                     endif

                  endif

                  if( ir .lt. 1 .or. ir .ge. nr + 1 .or.
     &                iz .lt. 1 .or. iz .ge. nz + 1 ) then
                        if( nrpt .eq. 2 ) then
                           nrpt = -1
                           ir = irrpt
                           iz = izrpt
                        else
                           goto 55
                        endif
                  endif

               endif

*-----------------------------------------------------------------------
*              surface crossing omni or spectrum current
*-----------------------------------------------------------------------

                  if( itout(m) .eq. 2 .or. itout(m) .eq. 5 ) then

                     tlw = oldwt
                     if( nrpt .eq. 2 ) tlwrpt = tlw   ! S.Abe 2018/10/18

*-----------------------------------------------------------------------
*              surface crossing flux spectrum
*-----------------------------------------------------------------------

                  else if( itout(m) .eq. 1) then

                     if( jz .eq. 1 ) then

                        smu = ud(3)

                     else

                        smu = 1.0
                        dxy = ( xpp - x0 ) * ud(1)
     &                      + ( ypp - y0 ) * ud(2)
                        rxy = ( xpp - x0 )**2
     &                      + ( ypp - y0 )**2

                        if( rxy .gt. 0.0 ) smu = dxy / sqrt(rxy)

                     end if

                     amu = abs( smu )

                     if( amu .le. grab ) amu = hgrab

                     tlw = oldwt / amu
                     if( nrpt .eq. 2 ) tlwrpt = tlw   ! S.Abe 2018/10/18

*-----------------------------------------------------------------------
*              surface crossing with angle
*-----------------------------------------------------------------------

                  else if( itout(m) .ge. 8) then  ! T.Sato 2021/05/05

                     if( itangform(m) .ge. 1 .and.
     &                   itangform(m) .le. 3 ) then

                        smu = ud(itangform(m))

                     else

                        if( jz .eq. 1 ) then

                           smu = ud(3)

                        else

                           smu = 1.0
                           dxy = ( xpp - x0 ) * ud(1)
     &                         + ( ypp - y0 ) * ud(2)
                           rxy = ( xpp - x0 )**2
     &                         + ( ypp - y0 )**2

                           if( rxy .gt. 0.0 ) smu = dxy / sqrt(rxy)

                        end if

                     endif

*-----------------------------------------------------------------------

                     tlw = 0.d0
                     tlwrpt = 0.d0

*-----------------------------------------------------------------------

                     if( nrpt .eq. 2 ) then

                        if( itangform(m).ge.1 .and.
     &                      itangform(m).le.3 ) then
                           amu = smu
                        else
                           amu = -1.d0 * abs( smu )
                        endif

                        if( itaty(m) .gt. 0 ) then
                           if( amu .lt. ab(1) .or.
     &                         amu .gt. ab(na+1) ) then
                                 nrpt = 1
                                 goto 44
                           endif
                        else
                           if( amu .lt. cos( ab(na+1)/180.d0*pi ).or.
     &                         amu .gt. cos( ab(1)/180.d0*pi ) ) then
                              nrpt = 1
                              goto 44
                           endif
                        endif

                        do i = 1, na

                           if( itaty(m) .gt. 0 ) then
                              if( amu .ge. ab(i) .and.
     &                            amu .le. ab(i+1) )
     &                         goto 43
                           else
                              if( amu .ge. cos( ab(i+1)/180.d0*pi )
     &                           .and.
     &                            amu .le. cos( ab(i)/180.d0*pi ) )
     &                         goto 43
                           end if

                        end do

   43                   iarpt = i

                        tlwrpt = oldwt

                     endif

   44                continue

*-----------------------------------------------------------------------

                     if( itangform(m).ge.1 .and.
     &                   itangform(m).le.3 ) then
                        amu = smu
                     else
                        if( nrpt .eq. 0 ) then
                           amu = smu
                        elseif( nrpt .eq. -1 ) then
                           amu = -1.d0 * abs( smu )
                        else
                           amu = abs( smu )
                        endif
                     endif

                     if( itaty(m) .gt. 0 ) then

                        if( smu .lt. ab(1) ) goto 46
                        if( smu .gt. ab(na+1) ) goto 46

                     else

                        if( smu .lt. cos( ab(na+1) / 180.d0 * pi ) )
     &                     goto 46
                        if( smu .gt. cos( ab(1) / 180.d0 * pi ) )
     &                     goto 46

                     end if

                     do i = 1, na

                        if( itaty(m) .gt. 0 ) then

                           if( smu .ge. ab(i) .and.
     &                         smu .le. ab(i+1) ) goto 45

                        else

                           if( smu .ge. cos( ab(i+1) / 180.d0 * pi )
     &                        .and.
     &                         smu .le. cos( ab(i) / 180.d0 * pi ) )
     &                         goto 45

                           end if

                     end do

   45                ia = i

                     tlw = oldwt

*-----------------------------------------------------------------------

   46                continue

                     if( tlw .le. 0.d0 ) then

                        if( tlwrpt .gt. 0.d0 ) then

                           nrpt = -1
                           ia = iarpt
                           ir = irrpt
                           iz = izrpt
                           tlw = tlwrpt

                        else

                           goto 55

                        endif

                     endif

                     if(itout(m) .ge. 10) then ! T.Sato 2021/05/05
                        amu = abs( smu )
                        if( amu .le. grab ) amu = hgrab
                        tlw = tlw / amu
                        tlwrpt = tlwrpt / amu  ! not sure the purpose of tlwrpt
                     endif

*-----------------------------------------------------------------------
*              surface crossing forward or back ward current
*-----------------------------------------------------------------------

                  else if( itout(m) .eq. 3 .or. itout(m) .eq. 4 .or.
     &                     itout(m) .eq. 6 .or. itout(m) .eq. 7 ) then

                     if( nrpt .eq. 0 ) then

                        if( jz .eq. 1 ) then

                           smu = ud(3)

                        else

                           smu = ( xpp - x0 ) * ud(1)
     &                         + ( ypp - y0 ) * ud(2)

                        end if

                           tlw = 0.0

                        if( ( itout(m).eq.3 .or. itout(m).eq.6 ) .and.
     &                        smu .gt. 0 ) tlw = oldwt

                        if( ( itout(m).eq.4 .or. itout(m).eq.7 ) .and.
     &                        smu .lt. 0 ) tlw = oldwt

                     else

                        tlw = oldwt

                     endif

                  end if

*-----------------------------------------------------------------------
*                    booking
*-----------------------------------------------------------------------

               if( nrpt .eq. 0 ) nrpt = 1
               mrpt = iabs(nrpt)

               do irpt = 1, mrpt

                  if( irpt .eq. 2 ) then
                     ia = iarpt
                     ir = irrpt
                     iz = izrpt
                     tlw = tlwrpt
                  endif

*-----------------------------------------------------------------------

                  if( itenclo(m) .ne. 1 .and. jz .eq. 1 ) then
                     do im = 1, nm
                     do ip = 1, ipn



!OBINATA(2012.7.11): Ct = Ct + xi.wi
                        tzEVENT(ips(ip),ie,ia,it,izf(ir,iz),im) =
     &                  tzEVENT(ips(ip),ie,ia,it,izf(ir,iz),im)
     &                  + tlw * facm(im)

                     end do
                     end do

                     if (istdev .eq. 2) then

                       call setitrmin(itzmin,2,7,(/ie,ia,it,ir,iz,nm/))
                       call setitrmax(itzmax,2,7,(/ie,ia,it,ir,iz,nm/))


                     endif

                  else

                     do im = 1, nm
                     do ip = 1, ipn



!OBINATA(2012.7.11): Ct = Ct + xi.wi
                        trEVENT(ips(ip),ie,ia,it,irf(ir,iz),im) =
     &                  trEVENT(ips(ip),ie,ia,it,irf(ir,iz),im)
     &                  + tlw * facm(im)

                     end do
                     end do

                     if (istdev .eq. 2) then

                       call setitrmin(itrmin,2,7,(/ie,ia,it,ir,iz,nm/))
                       call setitrmax(itrmax,2,7,(/ie,ia,it,ir,iz,nm/))

                     endif
                  end if

               enddo   ! S.Abe 2018/10/18

   61          continue

            end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*     next position
*     for z-crossing jz = 1, r-crossing jz = 0
*-----------------------------------------------------------------------

   55 continue

            if( jz .eq. 1 ) then

                  izm = izm + izk
                  izc = izc + izk

            else if( jz .eq. 0 ) then

               if( idr .gt. irc ) then

                  irc = idr

               else

                  irc = idr - 1

               end if

            end if

            se  = ee

*-----------------------------------------------------------------------

      goto 50

*-----------------------------------------------------------------------
cc H.Iwase 2015/4/2 for electrons
cc                  change the continuous energy back to the step energy
cc T.Sato 2017/2/14, use "goto 999" to come here

 999    if( ityp.eq.12 .or. ityp.eq.13 )then

          e(ibke +no,ipomp+1) =  e(ibke +no,ipomp+1) +
     &                                      deinit  - denstepold
         ec(ibkec+no,ipomp+1) = ec(ibkec+no,ipomp+1) +
     &                                      deinit  - denstepnew


      endif

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine psufrz(m,np,nr,nz,ne,na,nt,nm,rm,zm,eb,ab,tb,
     &                  tr,tz,idasa)
*                                                                      *
*       output r-z scoring mesh crossing current & flux                *
*       last modified by K.Niita on 2003/04/15                         *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)
      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80
      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall21/ rtfac(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /istcut/ ist_cut, ist_bat

      common /fact02/ facmxr(itlmax), facmxz(itlmax) ! kitamura23/03/31

*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      common /talout/ itall

      character fname*100, fnume*3

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /mpi00/ npe, me

      common /tall50/ itlmt(itlmax), ite2l(itlmax)
*-----------------------------------------------------------------------

      dimension   rm(nr+1)
      dimension   zm(nz+1)
      dimension   eb(ne+1)
      dimension   ab(na+1)
      dimension   tb(nt+1)
      dimension   ew(ne)
      dimension   aw(na)
      dimension   tw(nt)

      dimension   tr(np,ne,na,nt,(nr+1)*nz,nm,2)
      dimension   tz(np,ne,na,nt,nr*(nz+1),nm,2)


*-----------------------------------------------------------------------

      dimension tott(np,2) ! frtati 2021/10/05 6 -> np

      character hsunit(16)*32

      data hsunit( 1) / '[1/cm^2/source]                 '/
      data hsunit( 2) / '[1/cm^2/MeV/source]             '/
      data hsunit( 3) / '[1/cm^2/Lethargy/source]        '/
      data hsunit( 4) / '[1/cm^2/sr/source]              '/
      data hsunit( 5) / '[1/cm^2/MeV/sr/source]          '/
      data hsunit( 6) / '[1/cm^2/Lethargy/sr/source]     '/
      data hsunit(11) / '[1/cm^2/nsec/source]            '/
      data hsunit(12) / '[1/cm^2/MeV/nsec/source]        '/
      data hsunit(13) / '[1/cm^2/Lethargy/nsec/source]   '/
      data hsunit(14) / '[1/cm^2/sr/nsec/source]         '/
      data hsunit(15) / '[1/cm^2/MeV/sr/nsec/source]     '/
      data hsunit(16) / '[1/cm^2/Lethargy/sr/nsec/source]'/

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      character cha*1
      data cha /"'"/

      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31

*-----------------------------------------------------------------------

      character rpa*1
      data rpa /'}'/

      character cname*7
      character dname*7
      character aname*3

      character yen*1

      real(8),allocatable :: ar_r(:),ar_z(:)

      include 'samepage_include/samepage000.inc'


*-----------------------------------------------------------------------

         irf(ir,iz) = ir + ( iz - 1 ) * ( nr + 1 )
         izf(ir,iz) = max(iz,1) + ( ir - 1 ) * ( nz + 1 )

*-----------------------------------------------------------------------
*        set mesh area
*-----------------------------------------------------------------------

               az(ir) = pi * ( rm(ir+1)**2 - rm(ir)**2 )

               ar(ir,iz) = 2.0 * pi * rm(ir)
     &                     * ( zm(iz+1) - zm(iz) )

*-----------------------------------------------------------------------

      include 'samepage_include/samepage001.inc'

! sumover
      az_sum = 0.0d0
      do ir=1,nr
        az_sum = az_sum + az(ir)
      enddo
      allocate (ar_r(nz),ar_z(nr+1))
      ar_r(:) = 0.0d0
      ar_z(:) = 0.0d0
      do iz = 1,nz
        do ir = 1,nr+1
          ar_r(iz) = ar_r(iz) + ar(ir,iz)
          ar_z(ir) = ar_z(ir) + ar(ir,iz)
        enddo
      enddo

      yen  = char(92)
      igsh = 0

*-----------------------------------------------------------------------

      if ( iMeVperu.eq. 1 ) then
         hsunit( 2) = '[1/cm^2/(MeV/n)/source]         '
         hsunit( 5) = '[1/cm^2/(MeV/n)/sr/source]      '
         hsunit(12) = '[1/cm^2/(MeV/n)/nsec/source]    '
         hsunit(15) = '[1/cm^2/(MeV/n)/sr/nsec/source] '
      end if

      if ( ite2l(m) .eq. 1 ) then   ! convert to energy to LET
         hsunit( 2) = '[1/cm^2/(keV/um)/source]        '
         hsunit( 5) = '[1/cm^2/(keV/um)/sr/source]     '
         hsunit(12) = '[1/cm^2/(keV/um)/nsec/source]   '
         hsunit(15) = '[1/cm^2/(keV/um)/sr/nsec/source]'
      end if
*-----------------------------------------------------------------------
      np_mxang = np
      if ( itmxang(m).ne.0 .and. iabs(itmxang(m)).lt.np ) then
        np_mxang = iabs(itmxang(m))
      end if

*-----------------------------------------------------------------------
*        current or flux
*-----------------------------------------------------------------------

            if( itout(m) .eq. 1 ) then

               cname = '  Flux '
               dname = '  flux '

            else if( itout(m) .eq. 2 ) then

               cname = 'Current'
               dname = 'current'

            else if( itout(m) .eq. 3 ) then

               cname = ' F-Curr'
               dname = ' f-curr'

            else if( itout(m) .eq. 4 ) then

               cname = ' B-Curr'
               dname = ' b-curr'

            else if( itout(m) .eq. 5 ) then

               cname = ' O-Curr'
               dname = ' o-curr'

            else if( itout(m) .eq. 6 ) then

               cname = 'OF-Curr'
               dname = 'of-curr'

            else if( itout(m) .eq. 7 ) then

               cname = 'OB-Curr'
               dname = 'ob-curr'

            else if( itout(m) .eq. 8 ) then

               cname = ' A-Curr'
               dname = ' a-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            else if( itout(m) .eq. 9 ) then

               cname = 'OA-Curr'
               dname = 'oa-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            else if( itout(m) .eq. 10 ) then ! T.Sato 2021/05/05

               cname = ' A-Flux'
               dname = ' a-flux'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            else if( itout(m) .eq. 11 ) then ! T.Sato 2021/05/05

               cname = 'OA-Flux'
               dname = 'oa-flux'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               call wtpname(m,np,chp,chq)

*-----------------------------------------------------------------------
*        c1 : nomalization for source
*-----------------------------------------------------------------------

            if ( rsouin .gt. 0d0 ) then ! S.H. (2016.7.26)

               c1 = 1.0d+0 / rsouin

            else

               c1 = 0d0

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 1 : /cm^2/source
*                    = 2 : /cm^2/MeV/source
*                    = 3 : /cm^2/Lethargy/source
*                    = 4 : /cm^2/source/SR
*                    = 5 : /cm^2/MeV/source/SR
*                    = 6 : /cm^2/Lethargy/source/SR
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  4 .or.
     &          itunt(m) .eq. 11 .or. itunt(m) .eq. 14 ) then

               do i = 1, ne

                  ew(i) = 1.d+0

               end do
               ew_sum = 1.0d0

            else if( itunt(m) .eq.  2 .or. itunt(m) .eq.  5 .or.
     &               itunt(m) .eq. 12 .or. itunt(m) .eq. 15 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &               itunt(m) .eq. 13 .or. itunt(m) .eq. 16 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

            end if

*-----------------------------------------------------------------------

            aw_sum = 0.0d0

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  2 .or.
     &          itunt(m) .eq.  3 .or. itunt(m) .eq. 11 .or.
     &          itunt(m) .eq. 12 .or. itunt(m) .eq. 13 ) then

               do i = 1, na

                  aw(i) = 1.d+0

               end do
               aw_sum = 1.0d0

            else

               do i = 1, na

                  if( itaty(m) .gt. 0 ) then

                     aw(i) = ( ab(i+1) - ab(i) ) * 2.0 * pi
                  else

                     aw(i) = ( cos( ab(i) / 180.d0 * pi )
     &                       - cos( ab(i+1) / 180.d0 * pi ) ) * 2.0 * pi

                  end if
                  aw_sum = aw_sum +  aw(i)

               end do

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 11, 12, 13, 14 : /.../nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .gt. 10 ) then

               do i = 1, nt

                  tw(i) = tb(i+1) - tb(i)

               end do
               tw_sum = tb(nt+1) - tb(1)

            else

               do i = 1, nt

                  tw(i) = 1.d+0

               end do
               tw_sum = 1.0d0

            end if

*-----------------------------------------------------------------------
*        relative error and  unit conversion
*-----------------------------------------------------------------------
*        z-crossing
*-----------------------------------------------------------------------

                  cmax = 0.0
                  cmin = 1.e+33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1


            facmxz(m) = 1.d0
            if( rtfac(m) .lt. 0.d0 ) then
               facmxz(m) = 0.d0

               do 101 im = 1, nm
               do 101 iz = 1, nz + 1
               do 101 ir = 1, nr
               do 101 it = 1, nt
               do 101 ia = 1, na
               do 101 ie = 1, ne
               do 101 ip = 1, np

                  if( tz(ip,ie,ia,it,izf(ir,iz),im,1) .gt. 0.0d0 ) then

                     fmaxfc = tz(ip,ie,ia,it,izf(ir,iz),im,1)
     &                               / az(ir) / ew(ie) / aw(ia) / tw(it)

                    if( facmxz(m) .lt. fmaxfc) facmxz(m) = fmaxfc

                  end if

  101          continue

            end if

            do 100 im = 1, nm
            do 100 iz = 1, nz + 1
            do 100 ir = 1, nr
            do 100 it = 1, nt
            do 100 ia = 1, na
            do 100 ie = 1, ne
            do 100 ip = 1, np

               if( tz(ip,ie,ia,it,izf(ir,iz),im,1) .gt. 0.0d0 ) then

                  call calc_stdev(m,Xa,sigx,
     &                         tz(ip,ie,ia,it,izf(ir,iz),im,1),
     &                         tz(ip,ie,ia,it,izf(ir,iz),im,2),
     &          abs(rtfac(m)/facmxz(m))/az(ir)/ew(ie)/aw(ia)/tw(it))

                  tz(ip,ie,ia,it,izf(ir,iz),im,1) = Xa
                  tz(ip,ie,ia,it,izf(ir,iz),im,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tz(ip,ie,ia,it,izf(ir,iz),im,1) .gt. cmax )
     &                          cmax = tz(ip,ie,ia,it,izf(ir,iz),im,1)

                  if( tz(ip,ie,ia,it,izf(ir,iz),im,1) .lt. cmin )
     &                          cmin = tz(ip,ie,ia,it,izf(ir,iz),im,1)

               else

                  isdz = 1
                  tz(ip,ie,ia,it,izf(ir,iz),im,2) = 0.0

               end if

! sumover
               fact_in = abs(rtfac(m)/facmxz(m))
               call psufrz_sumover_tz_stdev(0,m,ip,ie,ia,it,ir,iz,im,
     &                fact_in,ew(ie),aw(ia),tw(it),az(ir),
     &                ew_sum,aw_sum,tw_sum,az_sum)


  100       continue

*-----------------------------------------------------------------------
*        r-crossing
*-----------------------------------------------------------------------

            facmxr(m) = 1.d0
            if( rtfac(m) .lt. 0.d0 ) then
               facmxr(m) = 0.d0

               do 201 im = 1, nm
               do 201 iz = 1, nz
               do 201 ir = 1, nr + 1
               do 201 it = 1, nt
               do 201 ia = 1, na
               do 201 ie = 1, ne
               do 201 ip = 1, np

                  if( tr(ip,ie,ia,it,irf(ir,iz),im,1) .gt. 0.0d0 ) then

                     fmaxfc = tr(ip,ie,ia,it,irf(ir,iz),im,1)
     &                            / ar(ir,iz) / ew(ie) / aw(ia) / tw(it)

                    if( facmxr(m) .lt. fmaxfc) facmxr(m) = fmaxfc

                  end if

  201          continue

            end if

            do 200 im = 1, nm
            do 200 iz = 1, nz
            do 200 ir = 1, nr + 1
            do 200 it = 1, nt
            do 200 ia = 1, na
            do 200 ie = 1, ne
            do 200 ip = 1, np

               if( tr(ip,ie,ia,it,irf(ir,iz),im,1) .gt. 0.d0 ) then

                call calc_stdev(m,Xa,sigx,
     &                       tr(ip,ie,ia,it,irf(ir,iz),im,1),
     &                       tr(ip,ie,ia,it,irf(ir,iz),im,2),
     &        abs(rtfac(m)/facmxr(m))/ar(ir,iz)/ew(ie)/aw(ia)/tw(it))

                  tr(ip,ie,ia,it,irf(ir,iz),im,1) = Xa
                  tr(ip,ie,ia,it,irf(ir,iz),im,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

               else

                  isdz = 1
                  tr(ip,ie,ia,it,irf(ir,iz),im,2) = 0.0

               end if
! sumover
               fact_in = abs(rtfac(m)/facmxr(m))
               call psufrz_sumover_stdev(0,m,ip,ie,ia,it,ir,iz,im,
     &                fact_in,ew(ie),aw(ia),tw(it),ar(ir,iz),
     &                ew_sum,aw_sum,tw_sum,ar_r(iz),ar_z(ir))


  200       continue

            if( nobch .gt. ist_bat ) then
            if( rtstd(m) .gt. 0.d0 .and.
     &          stdm .gt. 0.d0 .and. isdz .eq. 0 .and.
     &          stdm .lt. rtstd(m) ) itstd(m) = 0
            end if

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do iax = 1, itfln(m)

!OBINATA(2012.7.12): output *.err
         noe = 1
         if ( any( itaxs(m,iax) .eq. (/ 7 /) )
     &       .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &     noe = 2

         do ioe = 1, noe

         if( itall .eq. 2 .and. nobch .lt. maxbch ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.7.12): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)(1:itfll(m,iax))//'.'//fnume
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         else if ( itall.eq.4 ) then
            write(fnume,'(i3.3)') nobch
            if ( npe.gt.1 ) write (fnume,'(i3.3)') nobch/(npe-1)
            if ( ioe.eq.1 ) then
              call mk_2dnumfn(ctfln(m,iax),fname,itfll(m,iax),
     &             nobch,maxbch,npe)
            else
              call mk_2dnerfn(ctfln(m,iax),fname,itfll(m,iax),fnume)
            end if

         else

!OBINATA(2012.7.12): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         end if

            iot = 31
            open(iot, file = fname, status = 'unknown' )

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

               call tcrsech(iot,m,iax,1)

*-----------------------------------------------------------------------
*        z-crossing
*        energy axis  ( LET axis )
*-----------------------------------------------------------------------

      include 'samepage_include/samepage002_pmeatr1z1.inc'
      include 'samepage_include/samepagechp_pmeatrz.inc'
      include 'samepage_include/samepageseti.inc'

            nrstepi_0 = nrstepi
            nzstepi_0 = nzstepi

               inum = 0

         if( itaxs(m,iax) .eq. 1 .or.
     &       itaxs(m,iax) .eq. 14 ) then

            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nzstepi, nrstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            do imi = 1, nm
            do izi = 1, nz + 1, nzstepi
            do iri = 1, nr, nrstepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
              im = imi
              iz = izi
              ir = iri
              ia = iai
              it = iti
              ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:"
     &           ,1p1e14.7,1p1e14.7)')
     &           facmxr(m),facmxz(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif


                if(nzstepi .eq. 1. and. nrstepi .eq. 1) then
                  write(iot,'("#   no. =",i3,3x,
     &            "iz surf =",i3,3x,
     &            "ir  =",i3,/
     &            "#   zmesh = ",1p1e13.4,/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iz, ir,
     &                        zm(iz),
     &                        rm(ir), rm(ir+nrstepi)
                else
                  write(iot,'("#   no. =",i3,3x,
     &            "iz surf =",i3," - ",i3,3x,
     &            "ir  =",i3," - ",i3/
     &            "#   zmesh = ",1p1e13.4,/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  inum, iz, iz+nzstepi-1,ir,ir+nrstepi-1,
     &                        zm(iz), ! zm(iz+nzstepi), T.Sato 2024/11/16, not necessary to write higher z
     &                        rm(ir), rm(ir+nrstepi)
                endif

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .ge. 8 ) then
                 if(nastepi .eq. 1) then
                   write(iot,'(
     &             "#  ia ="i3/
     &             "# ",a3," = (",
     &                          1p1e13.4,"  -",1p1e13.4,"  )")')
     &                         ia, aname, ab(ia), ab(ia+nastepi)
                  else
                   write(iot,'(
     &             "#  ia ="i3," - ",i3/
     &             "# ",a3," = (",
     &                          1p1e13.4,"  -",1p1e13.4,"  )")')
     &                         ia,ia+nastepi-1, aname,
     &                         ab(ia), ab(ia+nastepi)
                 endif
               end if
               if( ittty(m) .ne. 0 ) then
                 if(ntstepi .eq. 1) then
                   write(iot,'(
     &             "#  it ="i3/
     &             "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
                 else
                   write(iot,'(
     &             "#  it ="i3," - ",i3/
     &             "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, it+ntstepi-1,tb(it), tb(it+ntstepi)
                 endif
               end if

               if( itaxl(m) .eq. 0 ) then

                if ( itaxs(m,iax) .eq. 1 ) then   ! ccse 2022/09/30

                if ( iMeVperu.eq. 1 ) then
                  write(iot,'(/"x: Energy [MeV/n]")')
                else
                  write(iot,'(/"x: Energy [MeV]")')
                end if

                else if ( itaxs(m,iax) .eq. 14 ) then
                 write(iot,'(/"x: LET [keV/um]")')
                end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 13 .or. itunt(m) .eq. 16 .or.
     &             itety(m) .eq.  3 .or. itety(m) .eq.  5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatrzm_e.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 7 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,a1)')
     &                        cha, inum, iz, ir, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,",  it =",i3,a1)')
     &                        cha, inum, iz, ir, it, cha
                  end if

               else

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ia =",i3,a1)')
     &                        cha, inum, iz, ir, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ia =",i3,",  it =",i3,a1)')
     &                        cha, inum, iz, ir, ia, it, cha
                  end if

               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 7 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, iz, ir, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,",  it =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, iz, ir, it, itmnt(m,im), cha
                  end if

               else

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ia =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, iz, ir, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ia =",i3,",  it =",i3,
     &                        ",  mset =",i3,a1)')
     &                     cha, inum, iz, ir, ia, it, itmnt(m,im), cha
                  end if

               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir =",i3)') ir
                write(changelsub(5),'(",  iz surf =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(7) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------
                  areasum = 0.0d0
                  do i=iri,iri+nrstepi-1
                    areasum = areasum + az(i)
                  end do

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  z surface crossing"/
     &                        "  z     &=&",1pe13.4," [cm]"/
     &                        "  area  &=&",1pe13.4," [cm^2]"/
     &                        "  rmin  &=&",1pe13.4," [cm]"/
     &                        "  rmax  &=&",1pe13.4," [cm]")')
     &                        yen,
     &                        zm(iz), areasum, rm(ir), rm(ir+nrstepi)
               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif


               if( itout(m) .gt. 7 ) then
                  write(iot,'(
     &                        "  amin  &=&",1pe13.4,/
     &                        "  amax  &=&",1pe13.4)')
     &                        ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        z-crossing
*        r axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 6 ) then

            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nzstepi, nrstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            do imi = 1, nm
            do izi = 1, nz + 1, nzstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
              im = imi
              iz = izi
              ie = iei
              ia = iai
              it = iti
              ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:"
     &           ,1p1e14.7,1p1e14.7)')
     &           facmxr(m),facmxz(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif


                if(nzstepi .eq. 1) then
                  write(iot,'("#   no. =",i3,3x,
     &            "iz surf =",i3)')
     &                        inum, iz
                else
                  write(iot,'("#   no. =",i3,3x,
     &            "iz surf =",i3," - ",i3)')
     &                        inum, iz,iz+nzstepi-1
                end if
               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05
                if(nestepi .eq. 1) then
                  write(iot,'(
     &            "#  ie ="i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(iei), eb(iei+nestepi)
                else
                  write(iot,'(
     &            "#  ie ="i3," - ",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, iei+nestepi-1,
     &                        eb(ie), eb(ie+nestepi)
                end if
               end if
               if( itout(m) .ge. 8 ) then
                if(nastepi .eq. 1) then
                  write(iot,'(
     &            "#  ia ="i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
                else
                  write(iot,'(
     &            "#  ia ="i3," - ",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, ia+nastepi-1,aname,
     &                        ab(ia), ab(ia+nastepi)
                end if
               end if
               if( ittty(m) .ne. 0 ) then
                if(ntstepi .eq. 1) then
                  write(iot,'(
     &            "#  it ="i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
                else
                  write(iot,'(
     &            "#  it ="i3," -",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, it+ntstepi-1, tb(it), tb(it+ntstepi)
                end if
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: r [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itrty(m) .eq. 3 .or. itrty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatrzm_r.inc'
      voll = samewtt(1)

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",   ie =",i3,a1)')
     &                     cha, inum, iz, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",   ie =",i3,",   it =",i3,a1)')
     &                     cha, inum, iz, ie, it, cha
                  end if

               else if( itout(m) .le. 7 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,a1)')
     &                     cha, inum, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,",   it =",i3,a1)')
     &                     cha, inum, iz, it, cha
                  end if

               else if( itout(m) .eq. 8 .or.
     &                  itout(m) .eq.10) then ! T.Sato 2021/05/05


                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iz, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",   it =",i3,a1)')
     &                     cha, inum, iz, ie, ia, it, cha
                  end if

               else if( itout(m) .eq. 9 .or.
     &                  itout(m) .eq.11) then ! T.Sato 2021/05/05

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iz, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ia =",i3,",   it =",i3,a1)')
     &                     cha, inum, iz, ia, it, cha
                  end if

               end if

*-----------------------------------------------------------------------
            else

               if( itout(m) .le. 4 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",   ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",   ie =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, it, itmnt(m,im), cha
                  end if

               else if( itout(m) .le. 7 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, it, itmnt(m,im), cha
                  end if

               else if( itout(m) .eq. 8 .or.
     &                  itout(m) .eq.10) then ! T.Sato 2021/05/05


                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, ia, it, itmnt(m,im), cha
                  end if

               else if( itout(m) .eq. 9 .or.
     &                  itout(m) .eq.11) then ! T.Sato 2021/05/05

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ia =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ia, it, itmnt(m,im), cha
                  end if

               end if

            end if

          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir =",i3)') ir
                write(changelsub(5),'(",  iz surf =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(4) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  z surface crossing"/
     &                     "  z     &=&",1pe13.4," [cm]")')
     &                     yen, zm(iz)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif


               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05


                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+nestepi)
                end if
               end if

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if


               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'(
     &                     " tot area  &=&",1pe13.4," [cm^2]")')
     &                     samewtt(1)

                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        z-crossing
*        z axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 5 ) then


            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nzstepi, nrstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            do imi = 1, nm
            do iei = 1, ne, nestepi
            do iri = 1, nr, nrstepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
            im = imi
            ie = iei
            ir = iri
            ia = iai
            it = iti
            ip = ipi


               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:"
     &           ,1p1e14.7,1p1e14.7)')
     &           facmxr(m),facmxz(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif


                  write(iot,'("#   no. =",i3,3x,"ir  =",i3/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ir,
     &                        rm(ir), rm(ir+nrstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05
                if(nestepi .eq. 1) then
                  write(iot,'(
     &            "#  ie ="i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(iei), eb(iei+nestepi)
                else
                  write(iot,'(
     &            "#  ie ="i3," - ",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                   ie,ie+nestepi-1, eb(ie), eb(ie+nestepi)
                endif
               end if
               if( itout(m) .ge. 8 ) then
                if(nastepi .eq. 1) then
                  write(iot,'(
     &            "#  ia ="i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
                else
                  write(iot,'(
     &            "#  ia ="i3," - ",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &             ia, ia+nastepi-1,name, ab(ia), ab(ia+nastepi)
                end if
               end if
               if( ittty(m) .ne. 0 ) then
                if(ntstepi .eq. 1) then
                  write(iot,'(
     &            "#  it ="i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
                else
                  write(iot,'(
     &            "#  it ="i3," - ",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, it+ntstepi-1 ,tb(it), tb(it+ntstepi)
                end if
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itzty(m) .eq. 3 .or. itzty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatrzm_z.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ir =",i3,a1)')
     &                     cha, inum, ie, ir, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ir =",i3,",   it =",i3,a1)')
     &                     cha, inum, ie, ir, it, cha
                  end if

               else if( itout(m) .le. 7 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,a1)')
     &                     cha, inum, ir, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, it, cha
                  end if

               else if( itout(m) .eq. 8 .or.
     &                  itout(m) .eq.10) then ! T.Sato 2021/05/05


                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie ="i3,
     &                     ",  ia ="i3,
     &                     ",  ir =",i3,a1)')
     &                     cha, inum, ie, ia, ir, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie ="i3,
     &                     ",  ia ="i3,
     &                     ",  ir =",i3,",   it =",i3,a1)')
     &                     cha, inum, ie, ia, ir, it, cha
                  end if

               else if( itout(m) .eq. 9 .or.
     &                  itout(m) .eq.11) then ! T.Sato 2021/05/05

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia ="i3,
     &                     ",   ir =",i3,a1)')
     &                     cha, inum, ia, ir, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia ="i3,
     &                     ",   ir =",i3,",   it =",i3,a1)')
     &                     cha, inum, ia, ir, it, cha
                  end if

               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ir =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ir, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ir =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ir, it, itmnt(m,im), cha
                  end if

               else if( itout(m) .le. 7 ) then

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, it, itmnt(m,im), cha
                  end if

               else if( itout(m) .eq. 8 .or.
     &                  itout(m) .eq.10) then ! T.Sato 2021/05/05

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie ="i3,
     &                     ",  ia ="i3,
     &                     ",  ir =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, ir, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie ="i3,
     &                     ",  ia ="i3,
     &                     ",  ir =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, ir, it, itmnt(m,im), cha
                  end if

               else if( itout(m) .eq. 9 .or.
     &                  itout(m) .eq.11) then ! T.Sato 2021/05/05

                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia ="i3,
     &                     ",   ir =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, ir, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia ="i3,
     &                     ",   ir =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, ir, it, itmnt(m,im), cha
                  end if

               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir =",i3)') ir
                write(changelsub(5),'(",  iz surf =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(5) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------
                  areasum = 0.0d0
                  do i=iri,iri+nrstepi-1
                    areasum = areasum + az(i)
                  end do
                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  z surface crossing"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]")')
     &                     yen, areasum,
     &                     rm(ir), rm(ir+nrstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif


               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if


               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05


                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+nestepi)
                end if
               end if
               if( itout(m) .ge. 8 ) then

                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        z-crossing
*        angle axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 .or.
     &            itaxs(m,iax) .eq. 10 ) then

            if(nrstepi .gt. nr) then
               nrstepi = nr
            end if
            do imi = 1, nm
            do izi = 1, nz + 1, nzstepi
            do iri = 1, nr, nrstepi
            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
            im = imi
            iz = izi
            ir = iri
            ie = iei
            it = iti
            ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:"
     &           ,1p1e14.7,1p1e14.7)')
     &           facmxr(m),facmxz(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif


                  write(iot,'("#   no. =",i3,3x,
     &            "iz surf =",i3,3x,
     &            "ir  =",i3/
     &            "#   zmesh = ",1p1e13.4,/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iz, ir,
     &                        zm(iz),
     &                        rm(ir), rm(ir+nrstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)


               if( itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05

                if(nestepi .eq. 1) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
                else
                  write(iot,'(
     &            "#  ie =",i3," - ",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                 ie, ie+nestepi-1,eb(ie), eb(ie+nestepi)
                end if
               end if
               if( ittty(m) .ne. 0 ) then
                if(ntstepi .eq. 1) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        it, tb(it), tb(it+ntstepi)
                else
                  write(iot,'(
     &            "#  it =",i3," - ",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                 it, it+ntstepi-1, tb(it), tb(it+ntstepi)
                end if
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  if( itaxs(m,iax) .eq. 8 ) then
                     write(iot,'(/"x: cos(",a1,"theta)")') yen
                  else
                     write(iot,'(/"x: ",a1,"theta  [deg]")') yen
                  end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatrzm_a.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .eq. 8 .or. itout(m) .eq.10) then ! T.Sato 2021/05/05

                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ie =",i3,a1)')
     &                        cha, inum, iz, ir, ie, cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ie =",i3,",   it =",i3,a1)')
     &                        cha, inum, iz, ir, ie, it, cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq. 11) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,a1)')
     &                        cha, inum, iz, ir, cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,",   it =",i3,a1)')
     &                        cha, inum, iz, ir, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ie =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, iz, ir, ie, itmnt(m,im), cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ie =",i3,",   it =",i3,
     &                        ",  mset =",i3,a1)')
     &                     cha, inum, iz, ir, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9  .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, iz, ir, itmnt(m,im), cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,",   it =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, iz, ir, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir =",i3)') ir
                write(changelsub(5),'(",  iz surf =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(6) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------
                  areasum = 0.0d0
                  do i=iri,iri+nrstepi-1
                    areasum = areasum + az(i)
                  end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  z surface crossing"/
     &                     "  z     &=&",1pe13.4," [cm]"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     zm(iz), areasum, rm(ir), rm(ir+nrstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if


               if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4,/
     &                     "  emax  &=&",1pe13.4)')
     &                     eb(ie), eb(ie+nestepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                     write(iot,'("e:")')
            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        z-crossing
*        time axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9 ) then

            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nzstepi, nrstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            do imi = 1, nm
            do izi = 1, nz + 1, nzstepi
            do iri = 1, nr, nrstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
            im = imi
            iz = izi
            ir = iri
            ie = iei
            ia = iai
            ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:"
     &           ,1p1e14.7,1p1e14.7)')
     &           facmxr(m),facmxz(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif


                  write(iot,'("#   no. =",i3,3x,
     &            "iz surf =",i3,3x,
     &            "ir  =",i3/
     &            "#   zmesh = ",1p1e13.4,/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iz, ir,
     &                        zm(iz),
     &                        rm(ir), rm(ir+nrstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05
                 if(nestepi .eq. 1) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        iei, eb(iei), eb(iei+nestepi)
                 else
                  write(iot,'(
     &            "#  ie =",i3," - ",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                 ie, ie+nestepi-1, eb(ie), eb(ie+nestepi)
                 endif
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Time [nsec]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( ittty(m) .eq. 3 .or. ittty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatrzm_t.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ie =",i3,a1)')
     &                        cha, inum, iz, ir, ie, cha
               else if( itout(m) .le. 7 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,a1)')
     &                        cha, inum, iz, ir, cha
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ie =",i3,
     &                        ",  ia =",i3,a1)')
     &                        cha, inum, iz, ir, ie, ia, cha
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ia =",i3,a1)')
     &                        cha, inum, iz, ir, ia, cha
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ie =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, iz, ir, ie, itmnt(m,im), cha
               else if( itout(m) .le. 7 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, iz, ir, itmnt(m,im), cha
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ie =",i3,
     &                        ",  ia =",i3,
     &                        ",  mset =",i3,a1)')
     &                     cha, inum, iz, ir, ie, ia, itmnt(m,im), cha
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  iz surf =",i3,
     &                        ",  ir =",i3,
     &                        ",  ia =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, iz, ir, ia, itmnt(m,im), cha
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir =",i3)') ir
                write(changelsub(5),'(",  iz surf =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(8) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  areasum = 0.0d0
                  do i=iri,iri+nrstepi-1
                     areasum = areasum + az(i)
                  end do
                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  z surface crossing"/
     &                     "  z     &=&",1pe13.4," [cm]"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     zm(iz), areasum, rm(ir), rm(ir+nrstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05

                  write(iot,'(
     &                     "  emin  &=&",1pe13.4,/
     &                     "  emax  &=&",1pe13.4)')
     &                     eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

         end if

*-----------------------------------------------------------------------
*        r-crossing
*        energy axis  ( LET axis )
*-----------------------------------------------------------------------
      if(nz.ge.1) then ! T.Sato 2023/04/13 ! in the case of nz = 0, no need to output r-crossing data
         if( itaxs(m,iax) .eq. 1 .or.
     &       itaxs(m,iax) .eq. 14 ) then

            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nzstepi, nrstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if

          nsame_1 = 1
          nr_1 = 1
          if(rm(1) .le. 0.0) then
            nr_1 =2
            if(iloopmode .eq. 6) then
              nrstepi = nrstepi - 1
              nsame_1 = 2
            end if
          end if
            do imi = 1, nm
            do iri = nr_1, nr + 1, nrstepi
            do izi = 1, nz, nzstepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
            im = imi
            ir = iri
            iz = izi
            ia = iai
            it = iti
            ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:"
     &           ,1p1e14.7,1p1e14.7)')
     &           facmxr(m),facmxz(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif


                  write(iot,'("#   no. =",i3,3x,
     &            "ir surf =",i3,3x,
     &            "iz  =",i3/
     &            "#   rmesh = ",1p1e13.4,/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ir, iz,
     &                        rm(ir),
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#    ia =",i3/
     &            "#   ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#    it =",i3/
     &            "#     t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                if ( itaxs(m,iax) .eq. 1 ) then   ! ccse 2022/09/30

                  write(iot,'(/"x: Energy [MeV]")')
                else if ( itaxs(m,iax) .eq. 14 ) then
                  write(iot,'(/"x: LET [keV/um]")')
                end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 13 .or. itunt(m) .eq. 16 .or.
     &             itety(m) .eq.  3 .or. itety(m) .eq.  5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatrzm_e_r_cross.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ir, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, iz, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,a1)')
     &                     cha, inum, ir, iz, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, iz, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, it, itmnt(m,im), cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir surf =",i3)') ir
                write(changelsub(5),'(",  iz =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(7) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------
                 areasum = 0.0d0
                 do i=izi,izi+nzstepi-1
                   areasum = areasum + ar(iri,i)
                 end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  r surface crossing"/
     &                     "  r     &=&",1pe13.4," [cm]"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                   rm(ir), areasum, zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if


               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        r-crossing
*        r axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 6 ) then

            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nzstepi, nrstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            do imi = 1, nm
            do iei = 1, ne, nestepi
            do izi = 1, nz, nzstepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
            im = imi
            ie = iei
            iz = izi
            ia = iai
            it = iti
            ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:"
     &           ,1p1e14.7,1p1e14.7)')
     &           facmxr(m),facmxz(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "iz  =",i3/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iz,
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05

                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: r [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itzty(m) .eq. 3 .or. itzty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatrzm_r_r.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ie, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ie, iz, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, iz, it, cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ie, ia, iz, cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ie, ia, iz, it, cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ia, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ia, iz, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, iz, itmnt(m,im), cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, iz, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir surf =",i3)') ir
                write(changelsub(5),'(",  iz =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(4) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------
                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  r surface crossing"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05


                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(iei), eb(iei+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+nestepi)
                end if
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'(
     &                     " tot area  &=&",1pe13.4," [cm^2]")')
     &                     samewtt(1)

                     write(iot,'("e:")')
            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        r-crossing
*        z axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 5 ) then

            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nzstepi, nrstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if

          nsame_1 = 1
          nr_1 = 1
          if(rm(1) .le. 0.0) then
            nr_1 =2
            if(iloopmode .eq. 6) then
              nrstepi = nrstepi - 1
              nsame_1 = 2
            end if
          end if

            do imi = 1, nm
            do iri = nr_1, nr + 1, nrstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
            im = imi
            ir = iri
            ie = iei
            ia = iai
            it = iti
            ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:"
     &           ,1p1e14.7,1p1e14.7)')
     &           facmxr(m),facmxz(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif


                  write(iot,'("#   no. =",i3,3x,
     &            "ir surf =",i3/
     &            "#   rmesh = ",1p1e13.4)')
     &                        inum, ir, rm(ir)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)


               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05

                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                          ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itrty(m) .eq. 3 .or. itrty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatrzm_z_r_cross.inc'
*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ie =",i3,a1)')
     &                     cha, inum, ir, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ie =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",    ir surf =",i3,a1)')
     &                     cha, inum, ir, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",    ir surf =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, it, cha
                  end if
               else if( itout(m) .eq. 8  .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ie =",i3,
     &                     ",   ia =",i3,a1)')
     &                     cha, inum, ir, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ie =",i3,
     &                     ",   ia =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ia =",i3,a1)')
     &                     cha, inum, ir, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ia =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ie =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",    ir surf =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",    ir surf =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ie =",i3,
     &                     ",   ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ie =",i3,
     &                     ",   ia =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, ie, ia, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir surf =",i3,
     &                     ",   ia =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir surf =",i3)') ir
                write(changelsub(5),'(",  iz =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(5) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  r surface crossing"/
     &                     "  r     &=&",1pe13.4," [cm]")')
     &                     yen, rm(ir)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if


               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05


                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+nestepi)
                end if
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'(
     &                     " tot area  &=&",1pe13.4," [cm^2]")')
     &                     samewtt(1)

                     write(iot,'("e:")')
            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        r-crossing
*        angle axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 .or.
     &            itaxs(m,iax) .eq. 10 ) then

            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nzstepi, nrstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
          nsame_1 = 1
          nr_1 = 1
          if(rm(1) .le. 0.0) then
            nr_1 =2
            if(iloopmode .eq. 6) then
              nrstepi = nrstepi - 1
              nsame_1 = 2
            end if
          end if

            do imi = 1, nm
            do iri = nr_1, nr + 1, nrstepi
            do izi = 1, nz, nzstepi
            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
            im = imi
            ir = iri
            iz = izi
            ie = iei
            it = iti
            ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:"
     &           ,1p1e14.7,1p1e14.7)')
     &           facmxr(m),facmxz(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ir surf =",i3,3x,
     &            "iz  =",i3/
     &            "#   rmesh = ",1p1e13.4,/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ir, iz,
     &                        rm(ir),
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)


               if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  if( itaxs(m,iax) .eq. 8 ) then
                     write(iot,'(/"x: cos(",a1,"theta)")') yen
                  else
                     write(iot,'(/"x: ",a1,"theta  [deg]")') yen
                  end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatrzm_a_r_cross.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,a1)')
     &                     cha, inum, ir, iz, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, iz, ie, it, cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ir, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, iz, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir surf =",i3)') ir
                write(changelsub(5),'(",  iz =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(6) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                 areasum = 0.0
                 do i=izi,izi+nzstepi-1
                   areasum = areasum + ar(iri,i)
                 end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  r surface crossing"/
     &                     "  r     &=&",1pe13.4," [cm]"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     rm(ir), areasum, zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if


               if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4,/
     &                     "  emax  &=&",1pe13.4)')
     &                     eb(ie), eb(ie+nestepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if
                     write(iot,'("e:")')
            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        r-crossing
*        time axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9 ) then

            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nzstepi, nrstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
          nsame_1 = 1
          nr_1 = 1
          if(rm(1) .le. 0.0) then
            nr_1 =2
            if(iloopmode .eq. 6) then
              nrstepi = nrstepi - 1
              nsame_1 = 2
            end if
          end if

            do imi = 1, nm
            do iri = nr_1, nr + 1, nrstepi
            do izi = 1, nz, nzstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
            im = imi
            ir = iri
            iz = izi
            ie = iei
            ia = iai
            ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:"
     &           ,1p1e14.7,1p1e14.7)')
     &           facmxr(m),facmxz(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif


                  write(iot,'("#   no. =",i3,3x,
     &            "ir surf =",i3,3x,
     &            "iz  =",i3/
     &            "#   rmesh = ",1p1e13.4,/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ir, iz,
     &                        rm(ir),
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Time [nsec]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( ittty(m) .eq. 3 .or. ittty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatrzm_t_r_cross.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,a1)')
     &                     cha, inum, ir, iz, ie, cha
               else if( itout(m) .le. 7 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ir, iz, cha
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,
     &                     ",  ia = ",i3,a1)')
     &                     cha, inum, ir, iz, ie, ia, cha
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,a1)')
     &                     cha, inum, ir, iz, ia, cha
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ie, itmnt(m,im), cha
               else if( itout(m) .le. 7 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, itmnt(m,im), cha
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,
     &                     ",  ia = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ie, ia, itmnt(m,im), cha
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir surf =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ia, itmnt(m,im), cha
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ir surf =",i3)') ir
                write(changelsub(5),'(",  iz =",i3)') iz
                write(changelsub(6),'(",  ia =",i3)') ia
                write(changelsub(7),'(",  ie =",i3)') ie
                write(changelsub(8),'(",  it =",i3)') it
                write(changelsub(9),'(a1)') cha

                changelsub(8) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(6) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(6) = " "
                  changelsub(7) = " "
                else if( itout(m) .le. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(7) = " "
                end if
               if(iloopmode .eq. 6) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------
                  areasum = 0.0d0
                  do i=izi,izi+nzstepi-1
                    areasum = areasum + ar(iri,i)
                  end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  r surface crossing"/
     &                     "  r     &=&",1pe13.4," [cm]"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     rm(ir), areasum, zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4,/
     &                     "  emax  &=&",1pe13.4)')
     &                     eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------

         end if
      end if ! T.Sato 2023/04/13, end of nz.ge.1 condition

*-----------------------------------------------------------------------

            call prestart(m,iot) !OBINATA(2012.7.12)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

      end do
      end do

*-----------------------------------------------------------------------
      include 'samepage_include/samepage999.inc'

      deallocate (ar_r,ar_z)

      return
      end


************************************************************************
*                                                                      *
      subroutine tsufxyz(ncol,m,np,nx,ny,nz,ne,na,nt,nm,
     &                   xm,ym,zm,eb,ab,tb,tr,trEVENT,
     &                   itrmax,itrmin)
*                                                                      *
*       surface crossing tally of current spectrum and flux            *
*       for xyz scoring mesh                                           *
*       last modified by K.Niita on 2015/10/30                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*             10 : geometry boundary crossing                          *
*             11 : termination by energy cut-off                       *
*             12 : termination by escape or leakage                    *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
*             15 : sequential transport only for tally                 *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use partmod, only: itmxpt, itpan, itpat, jtpat ! frtati 2021/10/05
*-----------------------------------------------------------------------

      implicit double precision( a-h, o-z )

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter( rlit = 29.97925d0 )
      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)

      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)

*-----------------------------------------------------------------------

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)
      common /tall73/ itenclo(itlmax), itangform(itlmax)

      common /tall82/ itcnth(9,itlmax)

      common /tall50/ itlmt(itlmax), ite2l(itlmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

*-----------------------------------------------------------------------

      dimension xm(nx+1)
      dimension ym(ny+1)
      dimension zm(nz+1)
      dimension eb(ne+1)
      dimension ab(na+1)
      dimension tb(nt+1)

      dimension tr(np,ne,na,nt,nx*ny*(nz+1),nm,2)
      dimension trEVENT(np,ne,na,nt,nx*ny*(nz+1),nm) !OBINATA(2012.5.29): as Ct
      real(8),allocatable,save:: tr0(:,:,:,:,:,:)   !OBINATA(2012.5.29): as C
      dimension itrmax(8),itrmin(8)


*-----------------------------------------------------------------------

      dimension ud(3)

      data dmax  /1.0d+19/
      data dmax0 /1.0d+18/

*-----------------------------------------------------------------------

      data grab,hgrab / 0.00349d+0, 0.001745d+0 /

      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt

*-----------------------------------------------------------------------

      dimension facm(6)

*-----------------------------------------------------------------------

      common /stat / istdev, irestart, ireschk
      common /cparm/ maxbch,maxcas

*-----------------------------------------------------------------------

      common /paraj/ mstz(300), parz(300)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)

      common /egs5cmn10/denstepold,denstepnew,deinit
      real*8            denstepold,denstepnew,deinit
!$OMP THREADPRIVATE(/egs5cmn10/)

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

*------------------------------------------------------------------------

      common /spred/ nspred, nwsprd, nedisp, itstep, ndedx

*-----------------------------------------------------------------------

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny

*-----------------------------------------------------------------------

      if (istdev .eq. 2) then

        call readitrminmax8(itrmin,itrmax,(/np,ne,na,nt,nx,ny,nz+1,nm/),
     &                      mnp,mne,mna,mnt,mnx,mny,mnz,mnm,
     &                      mxp,mxe,mxa,mxt,mxx,mxy,mxz,mxm)

      endif

*-----------------------------------------------------------------------
* check of history counter
*-----------------------------------------------------------------------
      ihistcount = 0 ! history counter index, T.Sato 2022/12/20 for speed up
         if ( ncol.eq.0 .or. ncol.eq.4 ) then
           do i = 1, 3
             if( itcnth(i,m) .eq. 1 ) then
               if( ncntmx(i) .lt. itcnth(i*2+2,m) .or.
     &             ncntmx(i) .gt. itcnth(i*2+3,m) ) then
                 ihistcount = 1
               end if
             end if
           end do
         end if

*-----------------------------------------------------------------------
*        source particle or end of batch ( in case of istdev = 2 )
*
* OBINATA(2012.7.11): change tr(,,,,,3) to trEVENT(,,,,)
*-----------------------------------------------------------------------

         if (( ncol .eq. 0 .or. ncol .eq. 4 )
     &                              .and. istdev .eq. 2) then
           if ((nocas.gt.1.or.ncol.eq.0) .and. ihistcount.ne.1 ) then

             do im = mnm,mxm
             do iz = mnz,mxz
             do iy = mny,mxy
             do ix = mnx,mxx
               ixyz = icf(ix,iy,iz)
               do ie = mne,mxe
                 tr(:,ie,:,:,ixyz,im,1) = tr(:,ie,:,:,ixyz,im,1)
     &                                  + trEVENT(:,ie,:,:,ixyz,im)
                 tr(:,ie,:,:,ixyz,im,2) = tr(:,ie,:,:,ixyz,im,2)
     &                                  + trEVENT(:,ie,:,:,ixyz,im) ** 2

               enddo
             enddo
             enddo
             enddo
             enddo

             call tsufxyz_sumover(m,1,
     &                   np, ne, na, nt, nx, ny, nz, nm,  trEVENT)

           end if

           do im = mnm,mxm
           do iz = mnz,mxz
           do iy = mny,mxy
           do ix = mnx,mxx
             ixyz = icf(ix,iy,iz)
             do ie = mne,mxe
               trEVENT(:,ie,:,:,ixyz,im) = 0
             enddo
           enddo
           enddo
           enddo
           enddo

           call resetitrminmax(itrmin,itrmax,8,
     &                         (/np,ne,na,nt,nx,ny,nz+1,nm/))

         end if

*-----------------------------------------------------------------------
*        end of batch ( in case of istdev = 1 )
*
* OBINATA(2012.7.11): modificate for thread parallel
*-----------------------------------------------------------------------

         if ( ncol .eq. 0 .and. istdev .eq. 1) then
!$OMP MASTER
             allocate( tr0(np,ne,na,nt,nx*ny*(nz+1),nm) )
             tr0(:,:,:,:,:,:) = 0.d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP CRITICAL (tsufxyz_crit_ist1)
             tr0(:,:,:,:,:,:) = tr0(:,:,:,:,:,:) + trEVENT(:,:,:,:,:,:)
!$OMP END CRITICAL (tsufxyz_crit_ist1)
!$OMP BARRIER
!$OMP MASTER
             tr(:,:,:,:,:,:,1) = tr(:,:,:,:,:,:,1)
     &                         + tr0(:,:,:,:,:,:) / maxcas
             tr(:,:,:,:,:,:,2) = tr(:,:,:,:,:,:,2)
     &                         + ( tr0(:,:,:,:,:,:) / maxcas ) ** 2

              call tsufxyz_sumover(m,maxcas,
     &                   np, ne, na, nt, nx, ny, nz, nm,  tr0)

             deallocate( tr0 )
!$OMP END MASTER

           trEVENT(:,:,:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------

         small = parz(28) * 10.d0

*-----------------------------------------------------------------------
*        check of ncol
*-----------------------------------------------------------------------

         if( ncol .lt. 10 ) return

*-----------------------------------------------------------------------

            if( ncol .eq. 10 .or. ncol .eq. 12 ) then

               bsmall = small

            else

               bsmall = 0.0

            end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncnt(ibknct+i,no,ipomp+1) .lt. itcnt(i*2+2,m) .or.
     &                ncnt(ibknct+i,no,ipomp+1) .gt. itcnt(i*2+3,m) )
     &                 return

               end if

            end do

*-----------------------------------------------------------------------
*        material for LET
*-----------------------------------------------------------------------

               if( itlmt(m) .gt. 0 ) then

                  lmat = idnm( itlmt(m) )

               else if( itlmt(m) .eq. 0 ) then

                  lmat = mat

               else

                  lmat = -idnm( -itlmt(m) )

               end if

*-----------------------------------------------------------------------
*        check of particles
*-----------------------------------------------------------------------

            call pcheck(m,np,ityp,ktyp,jtyp,ipn,ips)

               if( ipn .eq. 0 ) return

*-----------------------------------------------------------------------
cc H.Iwase 2015/4/2 for electrons
cc                  change the step energy to coutinuous energy
*-----------------------------------------------------------------------

      if( ityp.eq.12 .or. ityp.eq.13 )then

         e1 =  e(ibke  +no,ipomp+1) - deinit  + denstepold
         e2 = ec(ibkec +no,ipomp+1) - deinit  + denstepnew

        if (ite2l(m) .eq. 0) then   ! ccse 2022/08/31

         if( e1 .lt. eb(1)    ) return
         if( e2 .ge. eb(ne+1) ) return

        else if (ite2l(m) .eq. 1) then   ! energy to LTE
            call dedxas(e1,dedx,lmat,
     &                  ityp,ktyp,jtyp,rtyp)
            dedxl = dedx /10.0d0
            call dedxas(e2,dedx,lmat,
     &                  ityp,ktyp,jtyp,rtyp)
            dedxh = dedx /10.0d0
            if ( dedxl .gt. dedxh ) then   ! hight low swap
               dedx  = dedxh
               dedxh = dedxl
               dedxl = dedx
            end if

            if( dedxl .lt. eb(1)    ) return
            if( dedxh .ge. eb(ne+1) ) return
        end if

         tparti = abs(t(ibkt+no,ipomp+1))
         tpartf = abs(tc(ibktc+no,ipomp+1))
         if( tparti .ge. tb(nt+1) ) return
         if( tpartf .lt. tb(1) ) return

          e(ibke +no,ipomp+1) = e1
         ec(ibkec+no,ipomp+1) = e2

      endif

*-----------------------------------------------------------------------

      if( iMeVperu.eq.1 .and. ityp.ge.15 .and. ityp.le.19 ) then
           ebm = ktyp - ktyp / 1000000 * 1000000
      else
           ebm = 1.d0
      end if

*-----------------------------------------------------------------------
*        transform positions
*-----------------------------------------------------------------------

            call trnsxx(x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  xxa,yya,zza,itmtr(m,4))

            call trnsxx(xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                  zc(ibkzc+no,ipomp+1),
     &                  xxc,yyc,zzc,itmtr(m,4))

*-----------------------------------------------------------------------
*        check position : out of rainge
*-----------------------------------------------------------------------

            if( xxa+small .lt. xm(1) .and.
     &          xxc+small .lt. xm(1) ) goto 999

            if( yya+small .lt. ym(1) .and.
     &          yyc+small .lt. ym(1) ) goto 999

            if( zza+small .lt. zm(1) .and.
     &          zzc+small .lt. zm(1) ) goto 999

            if( xxa-small .ge. xm(nx+1) .and.
     &          xxc-small .ge. xm(nx+1) ) goto 999

            if( yya-small .ge. ym(ny+1) .and.
     &          yyc-small .ge. ym(ny+1) ) goto 999

            if( zza-small .ge. zm(nz+1) .and.
     &          zzc-small .ge. zm(nz+1) ) goto 999

*-----------------------------------------------------------------------
*           check of energy
*-----------------------------------------------------------------------

            if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &          itout(m) .eq.10) then ! T.Sato 2021/05/05

            if (ite2l(m) .eq. 0) then   ! ccse 2022/08/31
               if( e(ibke+no,ipomp+1)   .lt. eb(1)*ebm ) goto 999
               if( ec(ibkec+no,ipomp+1) .ge. eb(ne + 1)*ebm ) goto 999

            else if (ite2l(m) .eq. 1) then   ! convert to energy to LET
               call dedxas( e(ibke+no,ipomp+1),dedx,lmat,
     &                  ityp,ktyp,jtyp,rtyp)
               dedxl = dedx /10.0d0
               call dedxas(ec(ibkec+no,ipomp+1),dedx,lmat,
     &                     ityp,ktyp,jtyp,rtyp)
               dedxh = dedx /10.0d0
               if ( dedxl .gt. dedxh ) then   ! hight low swap
                  dedx  = dedxh
                  dedxh = dedxl
                  dedxl = dedx
               end if
               if( dedxl .lt. eb(1) ) goto 999
               if( dedxh .ge. eb(ne + 1) ) goto 999
            end if

            end if

*-----------------------------------------------------------------------

            icli = idgr(iblz1)

*-----------------------------------------------------------------------
*           check of time
*-----------------------------------------------------------------------

               tpart = abs(t(ibkt+no,ipomp+1))

               if( tpart  .ge. tb(nt+1) ) goto 999
               if( abs(tc(ibktc+no,ipomp+1)) .lt. tb(1) ) goto 999

*-----------------------------------------------------------------------
*        distance and unit vector ud(i)
*-----------------------------------------------------------------------

            dis = ( xxc - xxa )**2
     &          + ( yyc - yya )**2
     &          + ( zzc - zza )**2

            if( dis .le. small**2 ) goto 999

            dis = sqrt(dis)

            ud(1) = ( xxc - xxa ) / dis
            ud(2) = ( yyc - yya ) / dis
            ud(3) = ( zzc - zza ) / dis

*-----------------------------------------------------------------------
*     initial position, energy and initial range
*-----------------------------------------------------------------------

            tot = 0.0d0

            se  = e(ibke+no,ipomp+1)
            ee  = ec(ibkec+no,ipomp+1)

            xpp = xxa
            ypp = yya
            zpp = zza

            if( jtyp .ne. 0 .and. mat .gt. 0 )
     &          call rainge(se,rng,mat,ityp,ktyp,jtyp,rtyp)

*-----------------------------------------------------------------------
*     initial xyz position
*-----------------------------------------------------------------------

               ixm = 0
               ixc = 0

            if( ud(1) .gt. 0.0d0 ) ixk =  1
            if( ud(1) .eq. 0.0d0 ) ixk =  0
            if( ud(1) .lt. 0.0d0 ) ixk = -1

            if( ixk .ge. 0 ) then

               do i = 1, nx + 1
                  if( xm(i) .ge. xpp - small ) then
                         ixm = i - 1
                         ixc = i - 1
                         goto 32
                  end if
               end do

            else

               ixm = nx + 2
               ixc = nx + 1

               do i = 1, nx + 1
                  if( xm(i) .gt. xpp + small ) then
                         ixm = i
                         ixc = i - 1
                         goto 32
                  end if
               end do

            end if

   32       continue

*-----------------------------------------------------------------------

               iym = 0
               iyc = 0

            if( ud(2) .gt. 0.0d0 ) iyk =  1
            if( ud(2) .eq. 0.0d0 ) iyk =  0
            if( ud(2) .lt. 0.0d0 ) iyk = -1

            if( iyk .ge. 0 ) then

               do i = 1, ny + 1
                  if( ym(i) .ge. ypp - small ) then
                         iym = i - 1
                         iyc = i - 1
                         goto 35
                  end if
               end do

            else

               iym = ny + 2
               iyc = ny + 1

               do i = 1, ny + 1
                  if( ym(i) .gt. ypp + small ) then
                         iym = i
                         iyc = i - 1
                         goto 35
                  end if
               end do

            end if

   35       continue

*-----------------------------------------------------------------------

               izm = 0
               izc = 0

            if( ud(3) .gt. 0.0d0 ) izk =  1
            if( ud(3) .eq. 0.0d0 ) izk =  0
            if( ud(3) .lt. 0.0d0 ) izk = -1

            if( izk .ge. 0 ) then

               do i = 1, nz + 1
                  if( zm(i) .ge. zpp - small ) then
                         izm = i - 1
                         izc = i - 1
                         goto 38
                  end if
               end do

            else

               izm = nz + 2
               izc = nz + 1

               do i = 1, nz + 1
                  if( zm(i) .gt. zpp + small ) then
                         izm = i
                         izc = i - 1
                         goto 38
                  end if
               end do

            end if

   38       continue

*-----------------------------------------------------------------------
*     loop for finding mesh and booking upto total distance
*-----------------------------------------------------------------------

   50 continue

*-----------------------------------------------------------------------
*        check of energy
*-----------------------------------------------------------------------

         if (ite2l(m) .eq. 0) then   ! ccse 2022/08/31

         if( ( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &         itout(m) .eq. 10) .and.
     &       se .lt. eb(1)*ebm ) goto 999
         if( ( ( itout(m) .ge. 5 .and. itout(m) .le. 7 ) .or.
     &           itout(m) .eq. 9 .or. itout(m) .eq. 11 ) .and.
     &       se .lt. 0.0d0 ) goto 999

         else if (ite2l(m) .eq. 1) then   ! convert to energy to LET
            call dedxas( se,dedx,lmat,ityp,ktyp,jtyp,rtyp)
            dedx = dedx /10.0d0
            if( ( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &            itout(m) .eq.10 ) .and.
     &            dedx .lt. eb(1) ) goto 999
            if( ( ( itout(m) .ge. 5 .and. itout(m) .le. 7 ) .or.
     &              itout(m) .eq. 9 .or. itout(m) .eq. 11) .and.
     &          dedx .lt. 0.0d0 ) goto 999
         end if

*-----------------------------------------------------------------------
*        dx,dy,dz : distance to the nearest xyz mesh
*-----------------------------------------------------------------------

            dx = dmax
            dy = dmax
            dz = dmax

         if( ixk .ne. 0 ) then
            if( ixm + ixk .ge. 1 .and.
     &         ixm + ixk .le. nx + 1 ) then

               dx = abs( ( xm(ixm+ixk) - xpp ) / ud(1) )

            end if
         end if

         if( iyk .ne. 0 ) then
            if( iym + iyk .ge. 1 .and.
     &          iym + iyk .le. ny + 1 ) then

               dy = abs( ( ym(iym+iyk) - ypp ) / ud(2) )

            end if
         end if

         if( izk .ne. 0 ) then
            if( izm + izk .ge. 1 .and.
     &          izm + izk .le. nz + 1 ) then

               dz = abs( ( zm(izm+izk) - zpp ) / ud(3) )

            end if
         end if

*-----------------------------------------------------------------------
*        which boundary is the nearlist
*-----------------------------------------------------------------------

               dd = dx
               jj = 1

            if( dy .lt. dd ) then

               dd = dy
               jj = 2

            end if

            if( dz .lt. dd ) then

               dd = dz
               jj = 3

            end if

            if( dd .gt. dmax0 ) goto 999

*-----------------------------------------------------------------------
*        recalculate rng for ncol = 11 and nedisp
*-----------------------------------------------------------------------

            if( mat .gt. 0 .and. jtyp .ne. 0 .and. nedisp .ne. 0 .and.
     &           (  ityp .ne. 12 .and. ityp .ne. 13 ) ) then


               call rainge(se,rng,mat,ityp,ktyp,jtyp,rtyp)

               if( ncol .eq. 11 )  dis = tot + rng

            end if

*-----------------------------------------------------------------------
*         calculation is finished
*-----------------------------------------------------------------------

            if( tot + dd .ge. dis - small ) goto 999

*-----------------------------------------------------------------------
*        propagate position upto the boundary or the final point
*-----------------------------------------------------------------------

               tot = tot + dd
               xpp = xpp + dd * ud(1)
               ypp = ypp + dd * ud(2)
               zpp = zpp + dd * ud(3)

*-----------------------------------------------------------------------
*        final energy ee and range rng
*-----------------------------------------------------------------------

            if( jtyp .ne. 0 .and. mat .gt. 0 ) then

               delt = dd

               call ecol(ee,delt,se,rng,
     &                   mat,ityp,ktyp,jtyp,rtyp)

               rng = max( 0.0d0, rng - delt )
               ee  = max( 0.0d0, ee )
               ee  = min( se, ee )

               if (ite2l(m) .eq. 1) then   ! ccse 2022/08/31
                  call dedxas( ee,dedx,lmat,ityp,ktyp,jtyp,rtyp)
                  dedx = dedx /10.0d0
               end if

            end if

*-----------------------------------------------------------------------
*        time evolution
*-----------------------------------------------------------------------

            if( t(ibkt+no,ipomp+1) .gt. 0.d0 ) then

               ekin = ( se + ee ) / 2.0
               dist = dd
               timd = 0.0

               if( ekin .gt. 0.1 .or. rtyp .eq. 0.0d0 ) then

                  timd = dist * ( ekin + rtyp )
     &                 / sqrt( ekin * ( ekin + 2.0 * rtyp ) )
     &                 / rlit

               else if( ekin .gt. 0.0 ) then

                  timd = dist * sqrt( rtyp / 2.0 / ekin ) / rlit

               end if

                  tpart = tpart + timd

            end if

*-----------------------------------------------------------------------
*           booking position and next position
*-----------------------------------------------------------------------

            if( jj .eq. 1 ) then

               ix  = ixm + ixk
               iy  = iyc
               iz  = izc

            else if( jj .eq. 2 ) then

               ix  = ixc
               iy  = iym + iyk
               iz  = izc

            else if( jj .eq. 3 ) then

               ix  = ixc
               iy  = iyc
               iz  = izm + izk

            end if

*-----------------------------------------------------------------------
*        booking
*-----------------------------------------------------------------------

         if( ( itenclo(m) .eq. 1 .and.
     &         ix .ge. 1 .and. ix .le. nx + 1 .and.
     &         iy .ge. 1 .and. iy .le. ny + 1 .and.
     &         iz .ge. 1 .and. iz .le. nz + 1 ) .or.
     &       ( itenclo(m) .ne. 1 .and. jj .eq. 3 .and.
     &         ix .ge. 1 .and. ix .lt. nx + 1 .and.
     &         iy .ge. 1 .and. iy .lt. ny + 1 .and.
     &         iz .ge. 1 .and. iz .le. nz + 1 ) ) then

*-----------------------------------------------------------------------
*           final energy cell
*-----------------------------------------------------------------------

            if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &          itout(m) .eq.10) then ! T.Sato 2021/05/05

                     ie = 0

               if (ite2l(m) .eq. 0) then   ! ccse 2022/08/31
               do i = ne, 1, -1

                  if( ee .ge. eb(i)*ebm .and.
     &                ee .lt. eb(i+1)*ebm ) then

                     ie = i

                     goto 40

                  end if

               end do
               else if (ite2l(m) .eq. 1) then   ! convert to energy to LET
               do i = ne, 1, -1
                  if( dedx .ge. eb(i) .and.
     &                dedx .lt. eb(i+1) ) then
                     ie = i
                     goto 40
                  end if
               end do
               end if

   40          continue

*-----------------------------------------------------------------------
                     call fmfac(m,icli,ee,oldwt,facm)

*-----------------------------------------------------------------------
*           omni current case
*-----------------------------------------------------------------------

            else

               facm(:) = 1.0d0

               ie = 1

            end if

*-----------------------------------------------------------------------
*           tally
*-----------------------------------------------------------------------

            if( ie .gt. 0 ) then

                     ia = 1

*-----------------------------------------------------------------------
*              time
*-----------------------------------------------------------------------

                  do i = 1, nt

                     if( tpart .ge. tb(i) .and.
     &                   tpart .lt. tb(i+1) ) goto 60

                  end do

                     goto 61

   60                it = i

*-----------------------------------------------------------------------
*     identify ix, iy, iz for enclosure mode
*        nrpt =  0 : normal mode
*             = -1 : save only backward
*             =  1 : save only forward
*             =  2 : save both direction, backward id is saved as ixrpt
*-----------------------------------------------------------------------

               nrpt = 0

               if( itenclo(m) .eq. 1 ) then

                  nrpt = 1
                  tlwrpt = 0.d0

                  if( itout(m) .eq. 1 .or. itout(m) .eq. 2 .or.
     &                itout(m) .eq. 5 .or. itout(m) .ge. 8) then ! T.Sato 2021/05/05

                     nrpt = 2            ! backward wo rpt ni settei
                     iarpt = ia
                     ixrpt = ix
                     iyrpt = iy
                     izrpt = iz

                     if( jj .eq. 1 ) then
                        ix = ixc + ixk
                        ixrpt = ixc
                     elseif( jj .eq. 2 ) then
                        iy = iyc + iyk
                        iyrpt = iyc
                     elseif( jj .eq. 3 ) then
                        iz = izc + izk
                        izrpt = izc
                     endif

                     if( ixrpt .lt. 1 .or. ixrpt .ge. nx + 1 .or.
     &                   iyrpt .lt. 1 .or. iyrpt .ge. ny + 1 .or.
     &                   izrpt .lt. 1 .or. izrpt .ge. nz + 1 ) nrpt = 1

                  elseif( itout(m) .eq. 3 .or. itout(m) .eq. 6 ) then

                     if( jj .eq. 1 ) then
                        ix = ixc + ixk
                     elseif( jj .eq. 2 ) then
                        iy = iyc + iyk
                     elseif( jj .eq. 3 ) then
                        iz = izc + izk
                     endif

                  elseif( itout(m) .eq. 4 .or. itout(m) .eq. 7 ) then

                     if( jj .eq. 1 ) then
                        ix  = ixc
                     elseif( jj .eq. 2 ) then
                        iy  = iyc
                     elseif( jj .eq. 3 ) then
                        iz  = izc
                     endif

                  endif

                  if( ix .lt. 1 .or. ix .ge. nx + 1 .or.
     &                iy .lt. 1 .or. iy .ge. ny + 1 .or.
     &                iz .lt. 1 .or. iz .ge. nz + 1 ) then
                        if( nrpt .eq. 2 ) then
                           nrpt = -1
                           ix = ixrpt
                           iy = iyrpt
                           iz = izrpt
                        else
                           goto 55
                        endif
                  endif

               endif

*-----------------------------------------------------------------------
*              surface crossing omni or spectrum current
*-----------------------------------------------------------------------

                  if( itout(m) .eq. 2 .or. itout(m) .eq. 5 ) then

                     tlw = oldwt
                     if( nrpt .eq. 2 ) tlwrpt = tlw   ! S.Abe 2018/10/18

*-----------------------------------------------------------------------
*              surface crossing flux spectrum
*-----------------------------------------------------------------------

                  else if( itout(m) .eq. 1) then

                     smu = abs( ud(jj) )   ! S.Abe 2018/10/18, ud(3) -> ud(jj)

                     if( smu .le. grab ) smu = hgrab

                     tlw = oldwt / smu
                     if( nrpt .eq. 2 ) tlwrpt = tlw   ! S.Abe 2018/10/18

*-----------------------------------------------------------------------
*              surface crossing with angle
*-----------------------------------------------------------------------

                  else if( itout(m) .ge. 8 ) then  ! T.Sato 2021/05/05

                     tlw = 0.d0
                     tlwrpt = 0.d0

*-----------------------------------------------------------------------

                     if( nrpt .eq. 2 ) then

                        if( itangform(m) .ge. 1 .and.
     &                      itangform(m) .le. 3 ) then
                           smu = ud(itangform(m))
                        else
                           smu = -1.d0 * abs( ud(jj) )
                        endif

                        if( itaty(m) .gt. 0 ) then
                           if( smu .lt. ab(1) .or.
     &                         smu .gt. ab(na+1) ) then
                              nrpt = 1
                              goto 44
                           endif
                        else
                           if( smu .lt. cos( ab(na+1)/180.d0*pi ) .or.
     &                         smu .gt. cos( ab(1)/180.d0*pi ) ) then
                              nrpt = 1
                              goto 44
                           endif
                        endif

                        do i = 1, na
                           if( itaty(m) .gt. 0 ) then
                              if( smu .ge. ab(i) .and.
     &                            smu .le. ab(i+1) ) goto 43
                           else
                              if( smu .ge. cos( ab(i+1) / 180.d0 * pi )
     &                           .and.
     &                            smu .le. cos( ab(i) / 180.d0 * pi ) )
     &                            goto 43
                           end if
                        end do

   43                   iarpt = i

                        tlwrpt = oldwt

                     endif

   44                continue

*-----------------------------------------------------------------------

                        if( itangform(m) .ge. 1 .and.
     &                      itangform(m) .le. 3 ) then
                           smu = ud(itangform(m))
                        else
                           if( nrpt .eq. 0 ) then
                              smu = ud(jj)
                           elseif( nrpt .eq. -1 ) then
                              smu = -1.d0 * abs( ud(jj) )
                           else
                              smu = abs( ud(jj) )
                           endif
                        endif

                        if( itaty(m) .gt. 0 ) then

                           if( smu .lt. ab(1) ) goto 46
                           if( smu .gt. ab(na+1) ) goto 46

                        else

                           if( smu .lt.
     &                         cos( ab(na+1) / 180.d0 * pi ) )
     &                         goto 46
                           if( smu .gt. cos( ab(1) / 180.d0 * pi ) )
     &                         goto 46

                        end if

                     do i = 1, na

                        if( itaty(m) .gt. 0 ) then

                           if( smu .ge. ab(i) .and.
     &                         smu .le. ab(i+1) ) goto 45

                        else

                           if( smu .ge. cos( ab(i+1) / 180.d0 * pi )
     &                        .and.
     &                         smu .le. cos( ab(i) / 180.d0 * pi ) )
     &                         goto 45

                        end if

                     end do

   45                ia = i

                     tlw = oldwt

*-----------------------------------------------------------------------

   46                continue

                     if( tlw .le. 0.d0 ) then

                        if( tlwrpt .gt. 0.d0 ) then

                           nrpt = -1
                           ia = iarpt
                           ix = ixrpt
                           iy = iyrpt
                           iz = izrpt
                           tlw = tlwrpt

                        else

                           goto 55

                        endif

                     endif

                     if(itout(m) .ge. 10) then ! T.Sato 2021/05/05
                        amu = abs( smu )
                        if( amu .le. grab ) amu = hgrab
                        tlw = tlw / amu
                        tlwrpt = tlwrpt / amu  ! not sure the purpose of tlwrpt
                     endif

*-----------------------------------------------------------------------
*              surface crossing forward or back ward current
*-----------------------------------------------------------------------

                  else if( itout(m) .eq. 3 .or. itout(m) .eq. 4 .or.
     &                     itout(m) .eq. 6 .or. itout(m) .eq. 7 ) then

                     if( nrpt .eq. 0 ) then

                        smu = ud(jj)   ! S.Abe 2018/10/18, ud(3) -> ud(jj)

                        tlw = 0.0

                        if( ( itout(m).eq.3 .or. itout(m).eq.6 ) .and.
     &                        smu .gt. 0 ) tlw = oldwt

                        if( ( itout(m).eq.4 .or. itout(m).eq.7 ) .and.
     &                        smu .lt. 0 ) tlw = oldwt

                     else

                        tlw = oldwt

                     endif

                  end if

*-----------------------------------------------------------------------
*              booking
*-----------------------------------------------------------------------

               if( nrpt .eq. 0 ) nrpt = 1   ! normal case
               mrpt = iabs(nrpt)

               do irpt = 1, mrpt

                  if( irpt .eq. 2 ) then
                     ia = iarpt
                     ix = ixrpt
                     iy = iyrpt
                     iz = izrpt
                     tlw = tlwrpt
                  endif

*-----------------------------------------------------------------------


                  do im = 1, nm
                  do ip = 1, ipn



!OBINATA(2012.7.11): Ct = Ct + xi.wi
         trEVENT(ips(ip),ie,ia,it,icf(ix,iy,iz),im) =
     &   trEVENT(ips(ip),ie,ia,it,icf(ix,iy,iz),im) + tlw * facm(im)


                  end do
                  end do

                 if (istdev .eq. 2) then

cKN 2015/10/30 ??
                   call setitrmin(itrmin,2,8,(/ie,ia,it,ix,iy,iz, 1/))
                   call setitrmax(itrmax,2,8,(/ie,ia,it,ix,iy,iz,nm/))


                 endif

               enddo   ! S.Abe 2018/10/18

   61          continue

            end if

         end if

*-----------------------------------------------------------------------
*     next position
*-----------------------------------------------------------------------

   55 continue

            if( jj .eq. 1 ) then

               ixm = ixm + ixk
               ixc = ixc + ixk

            else if( jj .eq. 2 ) then

               iym = iym + iyk
               iyc = iyc + iyk

            else if( jj .eq. 3 ) then

               izm = izm + izk
               izc = izc + izk

            end if

               se  = ee

*-----------------------------------------------------------------------

      goto 50

*-----------------------------------------------------------------------
cc H.Iwase 2015/4/2 for electrons
cc                  change the continuous energy back to the step energy
cc T.Sato 2017/2/14, use "goto 999" to come here

 999   if( ityp.eq.12 .or. ityp.eq.13 )then

          e(ibke +no,ipomp+1) =  e(ibke +no,ipomp+1) +
     &                                  deinit  - denstepold
         ec(ibkec+no,ipomp+1) = ec(ibkec+no,ipomp+1) +
     &                                  deinit  - denstepnew

      endif


*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine psufxyz(m,np,nx,ny,nz,ne,na,nt,nm,xm,ym,zm,eb,ab,tb,
     &                   tr,igsh,idasa)
*                                                                      *
*       output xyz scoring mesh z surfacecrossing current & flux       *
*       last modified by K.Niita on 2015/10/30                         *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
      common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &                rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)
      common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &                rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)
      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80
      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall21/ rtfac(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /istcut/ ist_cut, ist_bat

      common /fact01/ facmax(itlmax) ! kitamura23/03/31

*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      common /talout/ itall

      character fname*100, fnume*3

*-----------------------------------------------------------------------
      character gfnam*100
      common /tall60/ itger(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /mpi00/ npe, me

      common /tall50/ itlmt(itlmax), ite2l(itlmax)
*-----------------------------------------------------------------------

      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)
      dimension   eb(ne+1)
      dimension   ab(na+1)
      dimension   tb(nt+1)
      dimension   ew(ne)
      dimension   aw(na)
      dimension   tw(nt)

      integer,allocatable :: ixyz(:)
      dimension   tr(np,ne,na,nt,nx*ny*(nz+1),nm,2)

*-----------------------------------------------------------------------

      dimension tott(np,2) ! frtati 2021/10/05 6 -> np

      character hsunit(16)*32

      data hsunit( 1) / '[1/cm^2/source]                 '/
      data hsunit( 2) / '[1/cm^2/MeV/source]             '/
      data hsunit( 3) / '[1/cm^2/Lethargy/source]        '/
      data hsunit( 4) / '[1/cm^2/sr/source]              '/
      data hsunit( 5) / '[1/cm^2/MeV/sr/source]          '/
      data hsunit( 6) / '[1/cm^2/Lethargy/sr/source]     '/
      data hsunit(11) / '[1/cm^2/nsec/source]            '/
      data hsunit(12) / '[1/cm^2/MeV/nsec/source]        '/
      data hsunit(13) / '[1/cm^2/Lethargy/nsec/source]   '/
      data hsunit(14) / '[1/cm^2/sr/nsec/source]         '/
      data hsunit(15) / '[1/cm^2/MeV/sr/nsec/source]     '/
      data hsunit(16) / '[1/cm^2/Lethargy/sr/nsec/source]'/

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      character cha*1
      data cha /"'"/

      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31

*-----------------------------------------------------------------------

      character dc2*4

      character rpa*1
      data rpa /'}'/

      character cname*7
      character dname*7
      character aname*3

      character yen*1

      real(8),allocatable :: ax_x(:),ax_y(:)

      include 'samepage_include/samepage000.inc'


*-----------------------------------------------------------------------

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny

*-----------------------------------------------------------------------
*        set mesh area
*-----------------------------------------------------------------------

               ax(ix,iy) = ( xm(ix+1) - xm(ix) )
     &                   * ( ym(iy+1) - ym(iy) )

      include 'samepage_include/samepage001.inc'

*-----------------------------------------------------------------------

      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

*-----------------------------------------------------------------------
! sumover

      allocate (ax_x(ny),ax_y(nx))
      ax_x(:) = 0.0d0
      ax_y(:) = 0.0d0
      do iy = 1,ny
        do ix = 1,nx
          ax_x(iy) = ax_x(iy) + ax(ix,iy)
          ax_y(ix) = ax_y(ix) + ax(ix,iy)
        enddo
      enddo

      if ( iMeVperu.eq. 1 ) then
         hsunit( 2) = '[1/cm^2/(MeV/n)/source]         '
         hsunit( 5) = '[1/cm^2/(MeV/n)/sr/source]      '
         hsunit(12) = '[1/cm^2/(MeV/n)/nsec/source]    '
         hsunit(15) = '[1/cm^2/(MeV/n)/sr/nsec/source] '
      end if

      if ( ite2l(m) .eq. 1 ) then   ! convert to energy to LET
         hsunit( 2) = '[1/cm^2/(keV/um)/source]        '
         hsunit( 5) = '[1/cm^2/(keV/um)/sr/source]     '
         hsunit(12) = '[1/cm^2/(keV/um)/nsec/source]   '
         hsunit(15) = '[1/cm^2/(keV/um)/sr/nsec/source]'
      end if
*-----------------------------------------------------------------------
      np_mxang = np
      if ( itmxang(m).ne.0 .and. iabs(itmxang(m)).lt.np ) then
        np_mxang = iabs(itmxang(m))
      end if

*-----------------------------------------------------------------------
*     for igshow
*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

            npg = np
            neg = ne
            nag = na
            ntg = nt

         else

            npg = 1
            neg = 1
            nag = 1
            ntg = 1

         end if

*-----------------------------------------------------------------------
*        current or flux
*-----------------------------------------------------------------------

            if( itout(m) .eq. 1 ) then

               cname = '  Flux '
               dname = '  flux '

            else if( itout(m) .eq. 2 ) then

               cname = 'Current'
               dname = 'current'

            else if( itout(m) .eq. 3 ) then

               cname = ' F-Curr'
               dname = ' f-curr'

            else if( itout(m) .eq. 4 ) then

               cname = ' B-Curr'
               dname = ' b-curr'

            else if( itout(m) .eq. 5 ) then

               cname = ' O-Curr'
               dname = ' o-curr'

            else if( itout(m) .eq. 6 ) then

               cname = 'OF-Curr'
               dname = 'of-curr'

            else if( itout(m) .eq. 7 ) then

               cname = 'OB-Curr'
               dname = 'ob-curr'

            else if( itout(m) .eq. 8 ) then

               cname = ' A-Curr'
               dname = ' a-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            else if( itout(m) .eq. 9 ) then

               cname = 'OA-Curr'
               dname = 'oa-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               call wtpname(m,np,chp,chq)

*-----------------------------------------------------------------------
*           itunt(m) = 1 : /cm^2/source
*                    = 2 : /cm^2/MeV/source
*                    = 3 : /cm^2/Lethargy/source
*                    = 4 : /cm^2/source/SR
*                    = 5 : /cm^2/MeV/source/SR
*                    = 6 : /cm^2/Lethargy/source/SR
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  4 .or.
     &          itunt(m) .eq. 11 .or. itunt(m) .eq. 14 ) then

               do i = 1, ne

                  ew(i) = 1.d+0

               end do
               ew_sum = 1.0d0

            else if( itunt(m) .eq.  2 .or. itunt(m) .eq.  5 .or.
     &               itunt(m) .eq. 12 .or. itunt(m) .eq. 15 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &               itunt(m) .eq. 13 .or. itunt(m) .eq. 16 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

            end if

*-----------------------------------------------------------------------

            aw_sum = 0.0d0
            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  2 .or.
     &          itunt(m) .eq.  3 .or. itunt(m) .eq. 11 .or.
     &          itunt(m) .eq. 12 .or. itunt(m) .eq. 13 ) then

               do i = 1, na

                  aw(i) = 1.d+0

               end do
               aw_sum = 1.0d0

            else

               do i = 1, na

                  if( itaty(m) .gt. 0 ) then

                     aw(i) = ( ab(i+1) - ab(i) ) * 2.0 * pi

                  else

                     aw(i) = ( cos( ab(i) / 180.d0 * pi )
     &                       - cos( ab(i+1) / 180.d0 * pi ) ) * 2.0 * pi

                  end if
                  aw_sum = aw_sum + aw(i)

               end do

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 11, 12, 13, 14 : /.../nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .gt. 10 ) then

               do i = 1, nt

                  tw(i) = tb(i+1) - tb(i)

               end do
               tw_sum = tb(nt+1) - tb(1)

            else

               do i = 1, nt

                  tw(i) = 1.d+0

               end do
               tw_sum = 1.0d0

            end if

*-----------------------------------------------------------------------
*        relative error and  unit conversion
*        c1 : nomalization for source
*-----------------------------------------------------------------------
*        z-crossing
*-----------------------------------------------------------------------

                  cmax = 0.0
                  cmin = 1.e+33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

         if( igsh .eq. 0 ) then

            if ( rsouin .gt. 0d0 ) then ! S.H. (2016.7.26)

               c1 = 1.0d+0 / rsouin

            else

               c1 = 0d0

            end if


            facmax(m) = 1.d0
            if( rtfac(m) .lt. 0.d0 ) then
               facmax(m) = 0.d0

               do 101 im = 1, nm
               do 101 iz = 1, nz + 1
               do 101 ix = 1, nx
               do 101 iy = 1, ny
               do 101 ie = 1, ne
               do 101 ia = 1, na
               do 101 ip = 1, np
               do 101 it = 1, nt

                  if(tr(ip,ie,ia,it,icf(ix,iy,iz),im,1) .gt. 0.d0 ) then

                     fmaxfc = tr(ip,ie,ia,it,icf(ix,iy,iz),im,1)
     &                            / ax(ix,iy) / ew(ie) / aw(ia) / tw(it)

                    if( facmax(m) .lt. fmaxfc) facmax(m) = fmaxfc

                  end if

  101          continue

            end if

            do 100 im = 1, nm
            do 100 iz = 1, nz + 1
            do 100 ix = 1, nx
            do 100 iy = 1, ny
            do 100 ie = 1, ne
            do 100 ia = 1, na
            do 100 ip = 1, np
            do 100 it = 1, nt

               if( tr(ip,ie,ia,it,icf(ix,iy,iz),im,1) .gt. 0.d0 ) then

                call calc_stdev(m,Xa,sigx,
     &                          tr(ip,ie,ia,it,icf(ix,iy,iz),im,1),
     &                          tr(ip,ie,ia,it,icf(ix,iy,iz),im,2),
     &           abs(rtfac(m)/facmax(m))/ax(ix,iy)/ew(ie)/aw(ia)/tw(it))

                  tr(ip,ie,ia,it,icf(ix,iy,iz),im,1) = Xa
                  tr(ip,ie,ia,it,icf(ix,iy,iz),im,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,ia,it,icf(ix,iy,iz),im,1) .gt. cmax )
     &                     cmax = tr(ip,ie,ia,it,icf(ix,iy,iz),im,1)

                  if( tr(ip,ie,ia,it,icf(ix,iy,iz),im,1) .lt. cmin )
     &                     cmin = tr(ip,ie,ia,it,icf(ix,iy,iz),im,1)

               else

                  isdz = 1
                  tr(ip,ie,ia,it,icf(ix,iy,iz),im,2) = 0.0

               end if

! sumover
               fact_in = abs(rtfac(m)/facmax(m))
               call psufxyz_sumover_stdev(0,m,ip,ie,ia,it,ix,iy,iz,im,
     &            fact_in,ew(ie),aw(ia),tw(it),ax(ix,iy),
     &            ew_sum,aw_sum,tw_sum,ax_x(iy),ax_y(ix))


  100       continue

            if( nobch .gt. ist_bat ) then
            if( rtstd(m) .gt. 0.d0 .and.
     &          stdm .gt. 0.d0 .and. isdz .eq. 0 .and.
     &          stdm .lt. rtstd(m) ) itstd(m) = 0
            end if

         end if

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do 900 iax = 1, itfln(m)

!OBINATA(2012.7.12): output *.err
         noe = 1
         if ( any( itaxs(m,iax) .eq. (/ 7 /) )
     &       .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &     noe = 2

         do 900 ioe = 1, noe

         if( igsh .ne. 0 .and.
     &     ( itaxs(m,iax) .lt. 7 .or.
     &       ittwo(m) .eq. 4 .or. ittwo(m) .eq. 5 ) ) goto 900

         if( itall .eq. 2 .and. nobch .lt. maxbch .and.
     &       igsh .eq. 0 ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.7.2): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)(1:itfll(m,iax))//'.'//fnume
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         else if ( itall.eq.4 .and. igsh.eq.0 ) then
            write(fnume,'(i3.3)') nobch
            if ( npe.gt.1 ) write (fnume,'(i3.3)') nobch/(npe-1)
            if ( ioe.eq.1 ) then
              call mk_2dnumfn(ctfln(m,iax),fname,itfll(m,iax),
     &             nobch,maxbch,npe)
            else
              call mk_2dnerfn(ctfln(m,iax),fname,itfll(m,iax),fnume)
            end if

         else

!OBINATA(2012.7.2): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         end if

            iot = 31
            open(iot, file = fname, status = 'unknown' )

*-----------------------------------------------------------------------
            gfnam = ctfln(m,iax)
            igfmn = itfll(m,iax)
            igser = itger(m)

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

               call tcrsech(iot,m,iax,1)

*-----------------------------------------------------------------------
*        z-crossing
*        energy axis  ( LET axis )
*-----------------------------------------------------------------------

      include 'samepage_include/samepage002_pmeatxyz1.inc'
      include 'samepage_include/samepagechp_pmeatxyz.inc'
      include 'samepage_include/samepageseti.inc'

            nxstepi_0 = nxstepi
            nystepi_0 = nystepi
            nzstepi_0 = nzstepi

         if( itaxs(m,iax) .eq. 1 .or.
     &       itaxs(m,iax) .eq. 14 ) then

               inum = 0

            if(nxstepi .gt. nx) then
               nxstepi = nx
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nystepi .gt. ny) then
               nystepi = ny
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if

            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz + 1, nzstepi
            do ixi = 1, nx, nxstepi
            do iyi = 1, ny, nystepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               ix = ixi
               iy = iyi
               ia = iai
               it = iti
               ip = ipi


               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "iz surf =",i3,3x,
     &            "ix  =",i3,3x,
     &            "iy  =",i3/
     &            "#   zmesh = ",1p1e13.4,/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iz, ix, iyi,
     &                        zm(iz),
     &                        xm(ix), xm(ix+nxstepi),
     &                        ym(iy), ym(iy+nystepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)


               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                if ( itaxs(m,iax) .eq. 1 ) then   ! ccse 2022/09/30

                if ( iMeVperu.eq. 1 ) then
                  write(iot,'(/"x: Energy [MeV/n]")')
                else
                  write(iot,'(/"x: Energy [MeV]")')
                end if

                else if ( itaxs(m,iax) .eq. 14 ) then
                 write(iot,'(/"x: LET [keV/um]")')
                end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 13 .or. itunt(m) .eq. 16 .or.
     &             itety(m) .eq.  3 .or. itety(m) .eq.  5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatxyzm_e.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, iz, ix, iy, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, ix, iy, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iz, ix, iy, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, ix, iy, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, iy, it, itmnt(m,im), cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, iy, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, iy, ia, it,
     &                     itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz surf =",i3)') iz
                write(changelsub(7),'(",  ia =",i3)') ia
                write(changelsub(8),'(",  ie =",i3)') ie
                write(changelsub(9),'(",  it =",i3)') it
                write(changelsub(10),'(a1)') cha

                changelsub(8) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(9) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(7) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(7) = " "
                  changelsub(8) = " "
                else if( itout(m) .eq. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(8) = " "
                end if
               if(iloopmode .eq. 3) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 4) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(8) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(9) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))//trim(changelsub(10))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------
                  areasum = 0.0d0
                  do ii=1,iyi,iyi+nystepi-1
                    do i=1,ixi,ixi+nxstepi-1
                        areasum = areasum + ax(i,ii)
                    end do
                  end do
                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  z surface crossing"/
     &                     "  z     &=&",1pe13.4," [cm]"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     zm(iz), areasum,
     &            xm(ix), xm(ix+nxstepi), ym(iy), ym(iy+nystepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do
            end do

            nxstepi = nxstepi_0
            nystepi = nystepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        z-crossing
*        x axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 3 ) then

               inum = 0

            if(nxstepi .gt. nx) then
               nxstepi = nx
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nystepi .gt. ny) then
               nystepi = ny
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz + 1, nzstepi
            do iyi = 1, ny, nystepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               iy = iyi
               ie = iei
               ia = iai
               it = iti
               ip = ipi


               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "iz surf =",i3,3x
     &            "iy  =",i3/
     &            "#   zmesh = ",1p1e13.4/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iz, iy, zm(iz),
     &                        ym(iy), ym(iy+nystepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: x [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itxty(m) .eq. 3 .or. itxty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatxyzm_x.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, iz, iy, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, iy, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  iy =",i3, a1)')
     &                     cha, inum, iz, iy, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, iy, it, cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iz, iy, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, iy, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  iy =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iz, iy, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  iy =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, iy, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, iy, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, iy, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  iy =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, iy, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, iy, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, iy, ie, ia, it,
     &                     itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  iy =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, iy, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  iy =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, iy, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz surf =",i3)') iz
                write(changelsub(7),'(",  ia =",i3)') ia
                write(changelsub(8),'(",  ie =",i3)') ie
                write(changelsub(9),'(",  it =",i3)') it
                write(changelsub(10),'(a1)') cha

                changelsub(4) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(9) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(7) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(7) = " "
                  changelsub(8) = " "
                else if( itout(m) .eq. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(8) = " "
                end if
               if(iloopmode .eq. 3) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 4) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(8) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(9) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))//trim(changelsub(10))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------
                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  z surface crossing"/
     &                     "  z     &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]")')
     &                     yen, zm(izi), ym(iyi), ym(iyi+nystepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if


               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+nestepi)
                end if
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'(
     &                     " tot area  &=&",1pe13.4," [cm^2]")')
     &                     samewtt(1)

                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do
            end do

            nxstepi = nxstepi_0
            nystepi = nystepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        z-crossing
*        y axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 4 ) then

               inum = 0

            if(nxstepi .gt. nx) then
               nxstepi = nx
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nystepi .gt. ny) then
               nystepi = ny
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz + 1, nzstepi
            do ixi = 1, nx, nxstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               ix = ixi
               ie = iei
               ia = iai
               it = iti
               ip = ipi


               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "iz surf =",i3,3x,
     &            "ix  =",i3/
     &            "#   zmesh = ",1p1e13.4/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iz, ix, zm(iz),
     &                        xm(ix), xm(ix+nxstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)


               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: y [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( ityty(m) .eq. 3 .or. ityty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatxyzm_y.inc'
      voll = samewtt(1)

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, iz, ix, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, ix, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3, a1)')
     &                     cha, inum, iz, ix, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, ix, it, cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iz, ix, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, ix, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iz, ix, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, ix, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, ie, ia, it,
     &                     itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz surf =",i3)') iz
                write(changelsub(7),'(",  ia =",i3)') ia
                write(changelsub(8),'(",  ie =",i3)') ie
                write(changelsub(9),'(",  it =",i3)') it
                write(changelsub(10),'(a1)') cha

                changelsub(5) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(9) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(7) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(7) = " "
                  changelsub(8) = " "
                else if( itout(m) .eq. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(8) = " "
                end if
               if(iloopmode .eq. 3) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 4) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(8) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(9) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))//trim(changelsub(10))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------


                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  z surface crossing"/
     &                     "  z     &=&",1pe13.4," [cm]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]")')
     &                     yen, zm(iz), xm(ix), xm(ix+nxstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+nestepi)
                end if
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'(
     &                     " tot area  &=&",1pe13.4," [cm^2]")')
     &                     voll

                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do
            end do

            nxstepi = nxstepi_0
            nystepi = nystepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        z-crossing
*        z axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 5 ) then

               inum = 0

            if(nxstepi .gt. nx) then
               nxstepi = nx
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nystepi .gt. ny) then
               nystepi = ny
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do ixi = 1, nx, nxstepi
            do iyi = 1, ny, nystepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ie = iei
               ia = iai
               ix = ixi
               iy = iyi
               it = iti
               ip = ipi


               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ix  =",i3,3x
     &            "iy  =",i3,/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ix, iy,
     &          xm(ix), xm(ix+nxstepi), ym(iy), ym(iy+nystepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)


               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itzty(m) .eq. 3 .or. itzty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatxyzm_z.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, ie, ix, iy, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, ix, iy, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, ix, iy, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iy, it, cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, ie, ia, ix, iy, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, ia, ix, iy, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, ia, ix, iy, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, ia, ix, iy, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ix, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ix, iy, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, ix, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, ix, iy, it,
     &                     itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, ix, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, ix, iy, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz  surf =",i3)') iz
                write(changelsub(7),'(",  ia =",i3)') ia
                write(changelsub(8),'(",  ie =",i3)') ie
                write(changelsub(9),'(",  it =",i3)') it
                write(changelsub(10),'(a1)') cha

                changelsub(6) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(9) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(7) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(7) = " "
                  changelsub(8) = " "
                else if( itout(m) .eq. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(8) = " "
                end if
               if(iloopmode .eq. 3) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 4) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(8) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(9) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))//trim(changelsub(10))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  areasum = 0.0d0
                  do ii=iyi,iyi+nystepi-1
                    do i=ixi,ixi+nxstepi-1
                      areasum = areasum + ax(i,ii)
                    end do
                  end do
                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  z surface crossing"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]")')
     &                     yen, areasum,
     &               xm(ix), xm(ix+nxstepi), ym(iy), ym(iy+nystepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if


               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+nestepi)
                end if
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if
                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do
            end do

            nxstepi = nxstepi_0
            nystepi = nystepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        z-crossing
*        angle axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 .or.
     &            itaxs(m,iax) .eq. 10 ) then

               inum = 0

            if(nxstepi .gt. nx) then
               nxstepi = nx
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nystepi .gt. ny) then
               nystepi = ny
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz + 1, nzstepi
            do ixi = 1, nx, nxstepi
            do iyi = 1, ny, nystepi
            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               ix = ixi
               iy = iyi
               ie = iei
               it = iti
               ip = ipi


               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "iz surf =",i3,3x,
     &            "ix  =",i3,3x,
     &            "iy  =",i3/
     &            "#   zmesh = ",1p1e13.4,/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iz, ix, iy,
     &                        zm(iz),
     &             xm(ix), xm(ix+nxstepi), ym(iy), ym(iy+nystepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)


               if( itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(iei+nestepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  if( itaxs(m,iax) .eq. 8 ) then
                     write(iot,'(/"x: cos(",a1,"theta)")') yen
                  else
                     write(iot,'(/"x: ",a1,"theta  [deg]")') yen
                  end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatxyzm_a.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, iz, ix, iy, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, ix, iy, ie, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, iz, ix, iy, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, ix, iy, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, iy, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, iy, ie, it,
     &                     itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, iy, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz surf =",i3)') iz
                write(changelsub(7),'(",  ia =",i3)') ia
                write(changelsub(8),'(",  ie =",i3)') ie
                write(changelsub(9),'(",  it =",i3)') it
                write(changelsub(10),'(a1)') cha

                changelsub(7) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(9) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(7) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(7) = " "
                  changelsub(8) = " "
                else if( itout(m) .eq. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(8) = " "
                end if
               if(iloopmode .eq. 3) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 4) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(8) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(9) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))//trim(changelsub(10))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  areasum = 0.0d0
                  do ii=iyi,iyi+nystepi-1
                    do i=ixi,ixi+nxstepi-1
                      areasum = areasum + ax(i,ii)
                    end do
                  end do
                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  z surface crossing"/
     &                     "  z     &=&",1pe13.4," [cm]"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     zm(iz), areasum,
     &           xm(ix), xm(ix+nxstepi), ym(iy), ym(iy+nystepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .eq. 8 ) then
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4,/
     &                     "  emax  &=&",1pe13.4)')
     &                     eb(ie), eb(ie+nestepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if
                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do
            end do

            nxstepi = nxstepi_0
            nystepi = nystepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        z-crossing
*        t axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9 ) then

               inum = 0

            if(nxstepi .gt. nx) then
               nxstepi = nx
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nystepi .gt. ny) then
               nystepi = ny
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz + 1, nzstepi
            do ixi = 1, nx, nxstepi
            do iyi = 1, ny, nystepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               ix = ixi
               iy = iyi
               ie = iei
               ia = iai
               ip = ipi



               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "iz surf =",i3,3x,
     &            "ix  =",i3,3x,
     &            "iy  =",i3/
     &            "#   zmesh = ",1p1e13.4,/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iz, ix, iy,
     &                        zm(iz),
     &            xm(ix), xm(ix+nxstepi), ym(iy), ym(iy+nystepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        iei, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        iai, aname, ab(ia), ab(ia+nastepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Time [nsec]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( ittty(m) .eq. 3 .or. ittty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatxyzm_t.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, iz, ix, iy, ie, cha
               else if( itout(m) .le. 7 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, iz, ix, iy, cha
               else if( itout(m) .eq. 8 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,",  ia =",i3,a1)')
     &                     cha, inum, iz, ix, iy, ie, ia, cha
               else if( itout(m) .eq. 9 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  ia =",i3,a1)')
     &                     cha, inum, iz, ix, iy, ia, cha
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, iy, ie, itmnt(m,im), cha
               else if( itout(m) .le. 7 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, iy, itmnt(m,im), cha
               else if( itout(m) .eq. 8 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ix, iy, ie, ia,
     &                     itmnt(m,im), cha
               else if( itout(m) .eq. 9 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  ia =",i3,
     &                     ",  mset =",i3, a1)')
     &                     cha, inum, iz, ix, iy, ia, itmnt(m,im), cha
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz =",i3)') iz
                write(changelsub(7),'(",  ia =",i3)') ia
                write(changelsub(8),'(",  ie =",i3)') ie
                write(changelsub(9),'(",  it =",i3)') it
                write(changelsub(10),'(a1)') cha

                changelsub(9) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(9) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(7) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(7) = " "
                  changelsub(8) = " "
                else if( itout(m) .eq. 8 ) then
                else if( itout(m) .eq. 9 ) then
                  changelsub(8) = " "
                end if
               if(iloopmode .eq. 3) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 4) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(8) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(9) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))//trim(changelsub(10))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------
                  areasum = 0.0d0
                  do ii = iyi, iyi+nystepi-1
                    do i = ixi, ixi+nxstepi-1
                      areasum = areasum + ax(i,ii)
                    end do
                  end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  z surface crossing"/
     &                     "  z     &=&",1pe13.4," [cm]"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     zm(izi), areasum,
     &           xm(ix), xm(ix+nxstepi), ym(iy), ym(iy+nystepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if


               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4,/
     &                     "  emax  &=&",1pe13.4)')
     &                     eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                     write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do
            end do

            nxstepi = nxstepi_0
            nystepi = nystepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        z-crossing
*        xy axis ( matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 7 ) then

               if( ittwo(m) .eq. 1 ) then

                  dc2 = 'h2: '

               else if( ittwo(m) .eq. 2 ) then

                  dc2 = 'hd: '

               else if( ittwo(m) .eq. 3 ) then

                  dc2 = 'hc: '

               else if( ittwo(m) .eq. 6 ) then

                  dc2 = 'hd2:'

               else if( ittwo(m) .eq. 7 ) then

                  dc2 = 'hc2:'

               end if

*-----------------------------------------------------------------------

               inum = 0

            do im = 1, nm
            do iz = 1, nz + 1
            do ip = 1, npg
            do ie = 1, neg
            do ia = 1, nag
            do it = 1, ntg

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if ( ip.gt.np_mxang ) then
                 write(iot,'( " SKIPPAGE:")')
                 inum = inum - 1
               end if

*-----------------------------------------------------------------------

                  write(iot,'("#   no. =",i3,3x,
     &            "iz surf =",i3/
     &            "#   zmesh = ",1p1e13.4)')
     &                        inum, iz, zm(iz)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+1)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+1)
               end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, iz, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ie =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,a1)')
     &                     cha, inum, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, it, cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iz, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iz, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ie =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, ia, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz surf =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: x [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: y [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = ( ym(ny+1) - ym(1) )
     &                  / ( xm(nx+1) - xm(1) )
                  xfac  = 0.9
                  afac  = 0.8
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( ( ( ittwo(m) .ge. 2 .and. ittwo(m) .le. 3 ) .or.
     &               ( ittwo(m) .ge. 6 .and. ittwo(m) .le. 7 ) ) .and.
     &             cmin .gt. 0.0 .and. cmax .gt. cmin ) then

                if (ioe .eq. 1 ) then
                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax
                else
                  write(iot,'( "set: c3[1.0e-4] c4[1.0]")')
                end if
                  write(iot,'( "p: cmin[c3] cmax[c4]")')
                  write(iot,'( "p: dmin(1e-31)")')

                  if( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &                 write(iot,'( "p: zlog")')

               end if

*-----------------------------------------------------------------------

               xmin = xm(1)
               xmax = xm(nx+1)
               ymin = ym(1)
               ymax = ym(ny+1)

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') xmin, xmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') ymin, ymax

               if( itanl(m) .gt. 0 .and. ioe .eq. 1 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               else if( itanl(m) .gt. 0 .and. ioe .ne. 1 ) then

                   call terrang(iot,m)

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

            end if

*-----------------------------------------------------------------------

                  write(iot,'("#  ny = ",i3,"   nx = ",i3)')
     &                         ny, nx

            if( ittwo(m) .ne. 4 ) then

                  write(iot,'( "# ( ( data(x,y), x = 1, nx ),",
     &                         " y = ny, 1, -1 )")')

            end if

*-----------------------------------------------------------------------

            if( ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) .and.
     &            igsh .eq. 0 ) then

                  write(iot,'(/a4," y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7," ; x = ",
     &            1p1g14.7," to ",1p1g14.7," by ",
     &            1p1g14.7," ;")') dc2,
     &            ym(ny) + rtydl(m)/2.0, ym(1) + rtydl(m)/2.0, rtydl(m),
     &            xm(1) + rtxdl(m)/2.0, xm(nx) + rtxdl(m)/2.0, rtxdl(m)

               write(iot,'(1p10e11.3)')
     &         ( ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), ix = 1, nx ),
     &             iy = ny, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# x          y        ",
     &                      "  flux       r.err")')

               do iy = 1, ny
               do ix = 1, nx

                  write(iot,'(1p3e11.3,0pf8.4)')
     &               xm(ix)  + rtxdl(m)/2.0,
     &               ym(iy)  + rtydl(m)/2.0,
     &               tr(ip,ie,ia,it,icf(ix,iy,iz),im,1),
     &               tr(ip,ie,ia,it,icf(ix,iy,iz),im,2)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

                  write(iot,'( "#   x = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7,/
     &                         "#   y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            xm(1) + rtxdl(m)/2.0, xm(nx) + rtxdl(m)/2.0, rtxdl(m),
     &            ym(1) + rtydl(m)/2.0, ym(ny) + rtydl(m)/2.0, rtydl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'y/x',( xm(ix) + rtxdl(m)/2.0, ix = 1, nx )

               do iy = ny, 1, -1

                  write(iot,'(1p1000e11.3)')
     &            ym(iy) + rtydl(m)/2.0,
     &            ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), ix = 1, nx )

               end do

            end if

*-----------------------------------------------------------------------
*        gshow
*-----------------------------------------------------------------------

         if( itgsh(m) .ne. 0 .and.
     &     ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) ) then

               write(iot,'(/"#",78("-"))')
               write(iot,'( "# gshow")')
               write(iot,'( "#",78("-"))')
               write(iot,'(/"p: legs[c5*0.875]")')

               zval = zm(iz)
               none = 1
               iaxs = 1
               iuni = itgsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nx+1,ny+1,none,xm,ym,zval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

         end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.05] form[c1/0.05] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.05] form[c1/0.05] ",
     &"nosp afac[c5*0.625] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if( itazl(m) .eq. 0 ) then

         write(iot,'("y: ",a7,1x,a32)') cname, hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  z surface crossing"/
     &                     "  z     &=&",1pe13.4," [cm]"/
     &                     "  part  &=&  ",a8)')
     &                     yen, zm(iz), chq(ip)
               if( itout(m) .le. 4 .or. itout(m) .eq. 8 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+1)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+1)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+1)
                end if
               end if

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,it)
               end if

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
               end if
                     write(iot,'("e:")')

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') xmin, xmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') ymin, ymax

      end if

*-----------------------------------------------------------------------

            end do
            end do
            end do
            end do
            end do
            end do

         end if

*-----------------------------------------------------------------------

            call prestart(m,iot) !OBINATA(2012.7.12)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

  900 continue

*-----------------------------------------------------------------------
      deallocate( ixyz )
      include 'samepage_include/samepage999.inc'

      deallocate (ax_x,ax_y)

      return
      end


************************************************************************
*                                                                      *
      subroutine psufxyz_rpp(m,np,nx,ny,nz,ne,na,nt,nm,xm,ym,zm,eb,ab,
     &                   tb,tr,igsh,idasa)
*                                                                      *
*       output xyz scoring mesh current & flux for enclosure mode      *
*       last modified by S.Abe on 2018/10/18                           *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
      common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &                rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)
      common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &                rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)
      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80
      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall21/ rtfac(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /tall73/ itenclo(itlmax), itangform(itlmax)

      common /istcut/ ist_cut, ist_bat

      common /fact01/ facmax(itlmax) ! kitamura23/03/31

*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      common /talout/ itall

      character fname*100, fnume*3

*-----------------------------------------------------------------------
      character gfnam*100
      common /tall60/ itger(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /tall50/ itlmt(itlmax), ite2l(itlmax)
*-----------------------------------------------------------------------

      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)
      dimension   eb(ne+1)
      dimension   ab(na+1)
      dimension   tb(nt+1)
      dimension   ew(ne)
      dimension   aw(na)
      dimension   tw(nt)
      integer,allocatable :: ixyz(:)
cKN 2015/10/30 ! 2023/2/17 Ogawa. (nz+1)th is not used in enclosure mode but for the sake of coherence with tsufxyz, +1 remains
      dimension   tr(np,ne,na,nt,nx*ny*(nz+1),nm,2)

*-----------------------------------------------------------------------

      dimension tott(np,2) ! frtati 2021/10/05 6 -> np

      character hsunit(16)*32

      data hsunit( 1) / '[1/cm^2/source]                 '/
      data hsunit( 2) / '[1/cm^2/MeV/source]             '/
      data hsunit( 3) / '[1/cm^2/Lethargy/source]        '/
      data hsunit( 4) / '[1/cm^2/sr/source]              '/
      data hsunit( 5) / '[1/cm^2/MeV/sr/source]          '/
      data hsunit( 6) / '[1/cm^2/Lethargy/sr/source]     '/
      data hsunit(11) / '[1/cm^2/nsec/source]            '/
      data hsunit(12) / '[1/cm^2/MeV/nsec/source]        '/
      data hsunit(13) / '[1/cm^2/Lethargy/nsec/source]   '/
      data hsunit(14) / '[1/cm^2/sr/nsec/source]         '/
      data hsunit(15) / '[1/cm^2/MeV/sr/nsec/source]     '/
      data hsunit(16) / '[1/cm^2/Lethargy/sr/nsec/source]'/

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      character cha*1
      data cha /"'"/

cfrtati 2021/10/05 chl, chm moved to partmod

cfrtati 2021/10/05 6 -> itmxpt
      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31


*-----------------------------------------------------------------------

      character dc2*4

      character rpa*1
      data rpa /'}'/

      character cname*7
      character dname*7
      character aname*3

      character yen*1

      real(8),allocatable :: ax_x(:,:),ax_y(:,:),ax_z(:,:)

      include 'samepage_include/samepage000.inc'

*-----------------------------------------------------------------------

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny

*-----------------------------------------------------------------------
*        set mesh area for rpp
*-----------------------------------------------------------------------

      ax(ix,iy,iz) = 2.d0 * ( xm(ix+1)-xm(ix) ) * ( ym(iy+1)-ym(iy) )
     &             + 2.d0 * ( ym(iy+1)-ym(iy) ) * ( zm(iz+1)-zm(iz) )
     &             + 2.d0 * ( zm(iz+1)-zm(iz) ) * ( xm(ix+1)-xm(ix) )

*-----------------------------------------------------------------------
      include 'samepage_include/samepage001.inc'

      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

*-----------------------------------------------------------------------

! sumover

      allocate (ax_x(ny,nz),ax_y(nx,nz),ax_z(nx,ny))
      ax_x(:,:) = 0.0d0
      ax_y(:,:) = 0.0d0
      ax_z(:,:) = 0.0d0
      do iz = 1,nz
        do iy = 1,ny
          do ix = 1,nx
            ax_x(iy,iz) = ax_x(iy,iz) + ax(ix,iy,iz)
            ax_y(ix,iz) = ax_y(ix,iz) + ax(ix,iy,iz)
            ax_z(ix,iy) = ax_z(ix,iy) + ax(ix,iy,iz)
          enddo
        enddo
      enddo

      if ( iMeVperu.eq. 1 ) then
         hsunit( 2) = '[1/cm^2/(MeV/n)/source]         '
         hsunit( 5) = '[1/cm^2/(MeV/n)/sr/source]      '
         hsunit(12) = '[1/cm^2/(MeV/n)/nsec/source]    '
         hsunit(15) = '[1/cm^2/(MeV/n)/sr/nsec/source] '
      end if

      if ( ite2l(m) .eq. 1 ) then   ! convert to energy to LET
         hsunit( 2) = '[1/cm^2/(keV/um)/source]        '
         hsunit( 5) = '[1/cm^2/(keV/um)/sr/source]     '
         hsunit(12) = '[1/cm^2/(keV/um)/nsec/source]   '
         hsunit(15) = '[1/cm^2/(keV/um)/sr/nsec/source]'
      end if
*-----------------------------------------------------------------------
      np_mxang = np
      if ( itmxang(m).ne.0 .and. iabs(itmxang(m)).lt.np ) then
        np_mxang = iabs(itmxang(m))
      end if

*-----------------------------------------------------------------------
*     for igshow
*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

            npg = np
            neg = ne
            nag = na
            ntg = nt

         else

            npg = 1
            neg = 1
            nag = 1
            ntg = 1

         end if

*-----------------------------------------------------------------------
*        current or flux
*-----------------------------------------------------------------------

            if( itout(m) .eq. 1 ) then

               cname = '  Flux '
               dname = '  flux '

            else if( itout(m) .eq. 2 ) then

               cname = 'Current'
               dname = 'current'

            else if( itout(m) .eq. 3 ) then

               cname = ' F-Curr'
               dname = ' f-curr'

            else if( itout(m) .eq. 4 ) then

               cname = ' B-Curr'
               dname = ' b-curr'

            else if( itout(m) .eq. 5 ) then

               cname = ' O-Curr'
               dname = ' o-curr'

            else if( itout(m) .eq. 6 ) then

               cname = 'OF-Curr'
               dname = 'of-curr'

            else if( itout(m) .eq. 7 ) then

               cname = 'OB-Curr'
               dname = 'ob-curr'

            else if( itout(m) .eq. 8 ) then

               cname = ' A-Curr'
               dname = ' a-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            else if( itout(m) .eq. 9 ) then

               cname = 'OA-Curr'
               dname = 'oa-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            else if( itout(m) .eq. 10 ) then ! T.Sato 2021/05/05

               cname = ' A-Flux'
               dname = ' a-flux'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            else if( itout(m) .eq. 11 ) then ! T.Sato 2021/05/05

               cname = 'OA-Flux'
               dname = 'oa-flux'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               call wtpname(m,np,chp,chq)

*-----------------------------------------------------------------------
*           itunt(m) = 1 : /cm^2/source
*                    = 2 : /cm^2/MeV/source
*                    = 3 : /cm^2/Lethargy/source
*                    = 4 : /cm^2/source/SR
*                    = 5 : /cm^2/MeV/source/SR
*                    = 6 : /cm^2/Lethargy/source/SR
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  4 .or.
     &          itunt(m) .eq. 11 .or. itunt(m) .eq. 14 ) then

               do i = 1, ne

                  ew(i) = 1.d+0

               end do
               ew_sum = 1.0d0

            else if( itunt(m) .eq.  2 .or. itunt(m) .eq.  5 .or.
     &               itunt(m) .eq. 12 .or. itunt(m) .eq. 15 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &               itunt(m) .eq. 13 .or. itunt(m) .eq. 16 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

            end if

*-----------------------------------------------------------------------

            aw_sum = 0.0d0
            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  2 .or.
     &          itunt(m) .eq.  3 .or. itunt(m) .eq. 11 .or.
     &          itunt(m) .eq. 12 .or. itunt(m) .eq. 13 ) then

               do i = 1, na

                  aw(i) = 1.d+0

               end do
               aw_sum = 1.0d0

            else

               do i = 1, na

                  if( itaty(m) .gt. 0 ) then

                     aw(i) = ( ab(i+1) - ab(i) ) * 2.0 * pi

                  else

                     aw(i) = ( cos( ab(i) / 180.d0 * pi )
     &                       - cos( ab(i+1) / 180.d0 * pi ) ) * 2.0 * pi

                  end if
                  aw_sum = aw_sum + aw(i)

               end do

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 11, 12, 13, 14 : /.../nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .gt. 10 ) then

               do i = 1, nt

                  tw(i) = tb(i+1) - tb(i)

               end do
               tw_sum = tb(nt+1) - tb(1)

            else

               do i = 1, nt

                  tw(i) = 1.d+0

               end do
               tw_sum = 1.0d0

            end if

*-----------------------------------------------------------------------
*        relative error and  unit conversion
*        c1 : nomalization for source
*-----------------------------------------------------------------------
*        rpp-crossing
*-----------------------------------------------------------------------

                  cmax = 0.0
                  cmin = 1.e+33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

         if( igsh .eq. 0 ) then

            if ( rsouin .gt. 0d0 ) then ! S.H. (2016.7.26)

               c1 = 1.0d+0 / rsouin

            else

               c1 = 0d0

            end if

            facmax(m) = 1.d0
            if( rtfac(m) .lt. 0.d0 ) then
               facmax(m) = 0.d0

               do 101 im = 1, nm
               do 101 iz = 1, nz
               do 101 ix = 1, nx
               do 101 iy = 1, ny
               do 101 ie = 1, ne
               do 101 ia = 1, na
               do 101 ip = 1, np
               do 101 it = 1, nt

                  if(tr(ip,ie,ia,it,icf(ix,iy,iz),im,1) .gt. 0.d0 ) then

                     fmaxfc = tr(ip,ie,ia,it,icf(ix,iy,iz),im,1)
     &                         / ax(ix,iy,iz) / ew(ie) / aw(ia) / tw(it)

                    if( facmax(m) .lt. fmaxfc) facmax(m) = fmaxfc

                  end if

  101          continue

            end if

            do 100 im = 1, nm
            do 100 iz = 1, nz
            do 100 ix = 1, nx
            do 100 iy = 1, ny
            do 100 ie = 1, ne
            do 100 ia = 1, na
            do 100 ip = 1, np
            do 100 it = 1, nt

               if( tr(ip,ie,ia,it,icf(ix,iy,iz),im,1) .gt. 0.d0 ) then

                call calc_stdev(m,Xa,sigx,
     &                          tr(ip,ie,ia,it,icf(ix,iy,iz),im,1),
     &                          tr(ip,ie,ia,it,icf(ix,iy,iz),im,2),
     &        abs(rtfac(m)/facmax(m))/ax(ix,iy,iz)/ew(ie)/aw(ia)/tw(it))

                  tr(ip,ie,ia,it,icf(ix,iy,iz),im,1) = Xa
                  tr(ip,ie,ia,it,icf(ix,iy,iz),im,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,ia,it,icf(ix,iy,iz),im,1) .gt. cmax )
     &                     cmax = tr(ip,ie,ia,it,icf(ix,iy,iz),im,1)

                  if( tr(ip,ie,ia,it,icf(ix,iy,iz),im,1) .lt. cmin )
     &                     cmin = tr(ip,ie,ia,it,icf(ix,iy,iz),im,1)

               else

                  isdz = 1
                  tr(ip,ie,ia,it,icf(ix,iy,iz),im,2) = 0.0

               end if

! sumover
               fact_in = abs(rtfac(m)/facmax(m))
               call psufxyz_sumover_rpp_stdev(0,m,
     &            ip,ie,ia,it,ix,iy,iz,im,
     &            fact_in,ew(ie),aw(ia),tw(it),ax(ix,iy,iz),
     &            ew_sum,aw_sum,tw_sum,
     &            ax_x(iy,iz),ax_y(ix,iz),ax_z(ix,iy))


  100       continue

            if( nobch .gt. ist_bat ) then
            if( rtstd(m) .gt. 0.d0 .and.
     &          stdm .gt. 0.d0 .and. isdz .eq. 0 .and.
     &          stdm .lt. rtstd(m) ) itstd(m) = 0
            end if

         end if

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do 900 iax = 1, itfln(m)

!OBINATA(2012.7.12): output *.err
         noe = 1
         if ( any( itaxs(m,iax) .eq. (/ 7 /) )
     &       .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &     noe = 2

         do 900 ioe = 1, noe

         if( igsh .ne. 0 .and.
     &     ( itaxs(m,iax) .lt. 7 .or.
     &       ittwo(m) .eq. 4 .or. ittwo(m) .eq. 5 ) ) goto 900

         if( itall .eq. 2 .and. nobch .lt. maxbch .and.
     &       igsh .eq. 0 ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.7.2): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)(1:itfll(m,iax))//'.'//fnume
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         else

!OBINATA(2012.7.2): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         end if

            iot = 31
            open(iot, file = fname, status = 'unknown' )

*-----------------------------------------------------------------------
            gfnam = ctfln(m,iax)
            igfmn = itfll(m,iax)
            igser = itger(m)

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

               call tcrsech(iot,m,iax,1)

*-----------------------------------------------------------------------
*        rpp-crossing
*        energy axis  ( LET axis )
*-----------------------------------------------------------------------

      include 'samepage_include/samepage002_pmeatxyz.inc'
      include 'samepage_include/samepagechp_pmeatxyz.inc'
      include 'samepage_include/samepageseti.inc'

            nxstepi_0 = nxstepi
            nystepi_0 = nystepi
            nzstepi_0 = nzstepi

         if( itaxs(m,iax) .eq. 1 .or.
     &       itaxs(m,iax) .eq. 14 ) then

               inum = 0

            if(nxstepi .gt. nx) then
               nxstepi = nx
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nystepi .gt. ny) then
               nystepi = ny
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz, nzstepi
            do ixi = 1, nx, nxstepi
            do iyi = 1, ny, nystepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               ix = ixi
               iy = iyi
               ia = iai
               it = iti
               ip = ipi


               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif


                  write(iot,'("#   no. =",i3,3x,
     &            "ix  =",i3,3x,
     &            "iy  =",i3/
     &            "iz  =",i3,3x,
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ix, iy, iz,
     &                        xm(ix), xm(ix+nxstepi),
     &                        ym(iy), ym(iy+nystepi),
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                if ( itaxs(m,iax) .eq. 1 ) then   ! ccse 2022/09/30

                if ( iMeVperu.eq. 1 ) then
                  write(iot,'(/"x: Energy [MeV/n]")')
                else
                  write(iot,'(/"x: Energy [MeV]")')
                end if

                else if ( itaxs(m,iax) .eq. 14 ) then
                 write(iot,'(/"x: LET [keV/um]")')
                end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 13 .or. itunt(m) .eq. 16 .or.
     &             itety(m) .eq.  3 .or. itety(m) .eq.  5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatxyzm_e_rpp.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ix, iy, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iy, iz, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, it, itmnt(m,im), cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ia, it,
     &                     itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz =",i3)') iz
                write(changelsub(7),'(",  ia =",i3)') ia
                write(changelsub(8),'(",  ie =",i3)') ie
                write(changelsub(9),'(",  it =",i3)') it
                write(changelsub(10),'(a1)') cha

                changelsub(8) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(9) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(7) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(7) = " "
                  changelsub(8) = " "
                else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then
                else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then
                  changelsub(8) = " "
                end if
               if(iloopmode .eq. 3) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 4) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(8) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(9) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))//trim(changelsub(10))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------


                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, ax(ix,iy,iz), xm(ix), xm(ix+nxstepi),
     &                ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do
            end do

            nxstepi = nxstepi_0
            nystepi = nystepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        rpp-crossing
*        x axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 3 ) then

               inum = 0

            if(nxstepi .gt. nx) then
               nxstepi = nx
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nystepi .gt. ny) then
               nystepi = ny
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz, nzstepi
            do iyi = 1, ny, nystepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               iy = iyi
               ie = iei
               ia = iai
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "iy  =",i3,3x
     &            "iz  =",i3/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iy, iz,
     &                 ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: x [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itxty(m) .eq. 3 .or. itxty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatxyzm_x_rpp.inc'
      voll = samewtt(1)

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, iy, iz, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  it =",i3,a1)')
     &                     cha, inum, iy, iz, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3, a1)')
     &                     cha, inum, iy, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,",  it =",i3,a1)')
     &                     cha, inum, iy, iz, it, cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iy, iz, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iy, iz, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iy, iz, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iy, iz, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, iz, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, iz, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8  .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, iz, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, iz, ie, ia, it,
     &                     itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9  .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, iz, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, iz, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz =",i3)') iz
                write(changelsub(7),'(",  ia =",i3)') ia
                write(changelsub(8),'(",  ie =",i3)') ie
                write(changelsub(9),'(",  it =",i3)') it
                write(changelsub(10),'(a1)') cha

                changelsub(4) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(9) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(7) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(7) = " "
                  changelsub(8) = " "
                else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then
                else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then
                  changelsub(8) = " "
                end if
               if(iloopmode .eq. 3) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 4) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(8) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(9) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))//trim(changelsub(10))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------
                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &             yen, ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05cKN 2018/01/22

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+nestepi)
                end if
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'(
     &                     " tot area  &=&",1pe13.4," [cm^2]")')
     &                     voll

                     write(iot,'("e:")')
            end do

            end do
            end do
            end do
            end do
            end do
            end do

            nxstepi = nxstepi_0
            nystepi = nystepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        rpp-crossing
*        y axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 4 ) then

               inum = 0

            if(nxstepi .gt. nx) then
               nxstepi = nx
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nystepi .gt. ny) then
               nystepi = ny
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz, nzstepi
            do ixi = 1, nx, nxstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               ix = ixi
               ie = iei
               ia = iai
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif


                  write(iot,'("#   no. =",i3,3x,
     &            "ix  =",i3,3x,
     &            "iz  =",i3/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ix, iz,
     &                        xm(ix), xm(ix+nxstepi),
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: y [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( ityty(m) .eq. 3 .or. ityty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatxyzm_y_rpp.inc'
      voll = samewtt(1)

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, ix, iz, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iz, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3, a1)')
     &                     cha, inum, ix, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iz, it, cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, ix, iz, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iz, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, ix, iz, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iz, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iz, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iz, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iz, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iz, ie, ia, it,
     &                     itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iz, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iz, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz =",i3)') iz
                write(changelsub(7),'(",  ia =",i3)') ia
                write(changelsub(8),'(",  ie =",i3)') ie
                write(changelsub(9),'(",  it =",i3)') it
                write(changelsub(10),'(a1)') cha

                changelsub(5) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(9) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(7) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(7) = " "
                  changelsub(8) = " "
                else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then
                else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then
                  changelsub(8) = " "
                end if
               if(iloopmode .eq. 3) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 4) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(8) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(9) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))//trim(changelsub(10))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------
                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &            yen, xm(ix), xm(ix+nxstepi), zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+nestepi)
                end if
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'(
     &                     " tot area  &=&",1pe13.4," [cm^2]")')
     &                     voll

                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do
            end do

            nxstepi = nxstepi_0
            nystepi = nystepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        rpp-crossing
*        z axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 5 ) then

               inum = 0

            if(nxstepi .gt. nx) then
               nxstepi = nx
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nystepi .gt. ny) then
               nystepi = ny
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do ixi = 1, nx, nxstepi
            do iyi = 1, ny, nystepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ie = iei
               ia = iai
               ix = ixi
               iy = iyi
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ix  =",i3,3x
     &            "iy  =",i3,/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ix, iy,
     &                  xm(ix), xm(ix+nxstepi), ym(iy), ym(iy+nystepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itzty(m) .eq. 3 .or. itzty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatxyzm_z_rpp.inc'
      voll = samewtt(1)

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, ie, ix, iy, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, ix, iy, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, ix, iy, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iy, it, cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, ie, ia, ix, iy, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, ia, ix, iy, it, cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, ia, ix, iy, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, ia, ix, iy, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ix, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ix, iy, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, ix, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, ix, iy, it,
     &                     itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, ix, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ia =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, ix, iy, it, itmnt(m,im), cha
                  end if
               end if
            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz =",i3)') iz
                write(changelsub(7),'(",  ia =",i3)') ia
                write(changelsub(8),'(",  ie =",i3)') ie
                write(changelsub(9),'(",  it =",i3)') it
                write(changelsub(10),'(a1)') cha

                changelsub(6) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(9) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(7) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(7) = " "
                  changelsub(8) = " "
                else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then
                else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then
                  changelsub(8) = " "
                end if
               if(iloopmode .eq. 3) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 4) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(8) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(9) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))//trim(changelsub(10))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                 xm(ix), xm(ix+nxstepi), ym(iy), ym(iy+nystepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+nestepi)
                end if
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'(
     &                     " tot area  &=&",1pe13.4," [cm^2]")')
     &                     voll

                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do
            end do

            nxstepi = nxstepi_0
            nystepi = nystepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        rpp-crossing
*        angle axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 .or.
     &            itaxs(m,iax) .eq. 10 ) then

               inum = 0

            if(nxstepi .gt. nx) then
               nxstepi = nx
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nystepi .gt. ny) then
               nystepi = ny
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz, nzstepi
            do iyi = 1, ny, nystepi
            do ixi = 1, nx, nxstepi
            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               iy = iyi
               ix = ixi
               ie = iei
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif


                  write(iot,'("#   no. =",i3,3x,
     &            "ix  =",i3,3x,
     &            "iy  =",i3,3x,
     &            "iz  =",i3/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ix, iy, iz,
     &                        xm(ix), xm(ix+nxstepi),
     &                        ym(iy), ym(iy+nystepi),
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  if( itaxs(m,iax) .eq. 8 ) then
                     write(iot,'(/"x: cos(",a1,"theta)")') yen
                  else
                     write(iot,'(/"x: ",a1,"theta  [deg]")') yen
                  end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

      include 'samepage_include/peatxyzm_a_rpp.inc'

*-----------------------------------------------------------------------
          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ie, it, cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ix, iy, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, iy, iz, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ie, it,
     &                     itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, it, itmnt(m,im), cha
                  end if
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz =",i3)') iz
                write(changelsub(7),'(",  ia =",i3)') ia
                write(changelsub(8),'(",  ie =",i3)') ie
                write(changelsub(9),'(",  it =",i3)') it
                write(changelsub(10),'(a1)') cha

                changelsub(7) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(9) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(7) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(7) = " "
                  changelsub(8) = " "
                else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then
                else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then
                  changelsub(8) = " "
                end if
               if(iloopmode .eq. 3) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 4) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(8) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(9) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))//trim(changelsub(10))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                 yen, ax(ix,iy,iz), xm(ix), xm(ix+nxstepi),
     &                 ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4,/
     &                     "  emax  &=&",1pe13.4)')
     &                     eb(ie), eb(ie+nestepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if
                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do
            end do

            nxstepi = nxstepi_0
            nystepi = nystepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        rpp-crossing
*        t axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9 ) then

               inum = 0

            if(nxstepi .gt. nx) then
               nxstepi = nx
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nystepi .gt. ny) then
               nystepi = ny
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nxstepi, nystepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do izi = 1, nz, nzstepi
            do iyi = 1, ny, nystepi
            do ixi = 1, nx, nxstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
               im = imi
               iz = izi
               iy = iyi
               ix = ixi
               ie = iei
               ia = iai
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ix =",i3,3x,
     &            "iy  =",i3,3x,
     &            "iz  =",i3/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ix, iy, iz,
     &                        xm(ix), xm(ix+nxstepi),
     &                        ym(iy), ym(iy+nystepi),
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Time [nsec]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( ittty(m) .eq. 3 .or. ittty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatxyzm_t_rpp.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ie, cha
               else if( itout(m) .le. 7 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ix, iy, iz, cha
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  ia =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ie, ia, cha
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,",  ia =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ia, cha
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ie, itmnt(m,im), cha
               else if( itout(m) .le. 7 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, itmnt(m,im), cha
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ie, ia,
     &                     itmnt(m,im), cha
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,",  ia =",i3,
     &                     ",  mset =",i3, a1)')
     &                     cha, inum, ix, iy, iz, ia, itmnt(m,im), cha
               end if

            end if
          else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                 changelsub(3) = " " ! T.Sato 2024/11/16, avoid NULL
                else
                  write(changelsub(3),'(",  mset ="i4)') itmnt(m,im)
                end if
                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz =",i3)') iz
                write(changelsub(7),'(",  ia =",i3)') ia
                write(changelsub(8),'(",  ie =",i3)') ie
                write(changelsub(9),'(",  it =",i3)') it
                write(changelsub(10),'(a1)') cha

                changelsub(9) = " "

                if( ittty(m) .eq. 0 ) then
                  changelsub(9) = " "
                end if
                if( itout(m) .le. 4 ) then
                  changelsub(7) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(7) = " "
                  changelsub(8) = " "
                else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then
                else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then
                  changelsub(8) = " "
                end if
               if(iloopmode .eq. 3) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 4) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 5) then
                  changelsub(6) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9 ) then
                  changelsub(7) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(8) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(9) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))//trim(changelsub(10))
                write(iot,'(/a)') trim(angeltitle)
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                   yen, ax(ix,iy,iz), xm(ix), xm(ix+nxstepi),
     &                  ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4,/
     &                     "  emax  &=&",1pe13.4)')
     &                     eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                     write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do
            end do

            nxstepi = nxstepi_0
            nystepi = nystepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        rpp-crossing
*        xy axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 7 ) then

               if( ittwo(m) .eq. 1 ) then

                  dc2 = 'h2: '

               else if( ittwo(m) .eq. 2 ) then

                  dc2 = 'hd: '

               else if( ittwo(m) .eq. 3 ) then

                  dc2 = 'hc: '

               else if( ittwo(m) .eq. 6 ) then

                  dc2 = 'hd2:'

               else if( ittwo(m) .eq. 7 ) then

                  dc2 = 'hc2:'

               end if

*-----------------------------------------------------------------------

               inum = 0

            do im = 1, nm
            do iz = 1, nz
            do ip = 1, npg
            do ie = 1, neg
            do ia = 1, nag
            do it = 1, ntg

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if ( ip.gt.np_mxang ) then
                 write(iot,'( " SKIPPAGE:")')
                 inum = inum - 1
               end if

*-----------------------------------------------------------------------

                  write(iot,'("#   no. =",i3,3x,
     &            "iz =",i3/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iz, zm(iz), zm(iz+1)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+1)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+1)
               end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, iz, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, it, cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iz, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iz, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iz, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ie, ia, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: x [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: y [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = ( ym(ny+1) - ym(1) )
     &                  / ( xm(nx+1) - xm(1) )
                  xfac  = 0.9
                  afac  = 0.8
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( ( ( ittwo(m) .ge. 2 .and. ittwo(m) .le. 3 ) .or.
     &               ( ittwo(m) .ge. 6 .and. ittwo(m) .le. 7 ) ) .and.
     &             cmin .gt. 0.0 .and. cmax .gt. cmin ) then

                if (ioe .eq. 1 ) then
                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax
                else
                  write(iot,'( "set: c3[1.0e-4] c4[1.0]")')
                end if
                  write(iot,'( "p: cmin[c3] cmax[c4]")')
                  write(iot,'( "p: dmin(1e-31)")')

                  if( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &                 write(iot,'( "p: zlog")')

               end if

*-----------------------------------------------------------------------

               xmin = xm(1)
               xmax = xm(nx+1)
               ymin = ym(1)
               ymax = ym(ny+1)

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') xmin, xmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') ymin, ymax

               if( itanl(m) .gt. 0 .and. ioe .eq. 1 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               else if( itanl(m) .gt. 0 .and. ioe .ne. 1 ) then

                   call terrang(iot,m)

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

            end if

*-----------------------------------------------------------------------

                  write(iot,'("#  ny = ",i3,"   nx = ",i3)')
     &                         ny, nx

            if( ittwo(m) .ne. 4 ) then

                  write(iot,'( "# ( ( data(x,y), x = 1, nx ),",
     &                         " y = ny, 1, -1 )")')

            end if

*-----------------------------------------------------------------------

            if( ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) .and.
     &            igsh .eq. 0 ) then

                  write(iot,'(/a4," y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7," ; x = ",
     &            1p1g14.7," to ",1p1g14.7," by ",
     &            1p1g14.7," ;")') dc2,
     &            ym(ny) + rtydl(m)/2.0, ym(1) + rtydl(m)/2.0, rtydl(m),
     &            xm(1) + rtxdl(m)/2.0, xm(nx) + rtxdl(m)/2.0, rtxdl(m)

               write(iot,'(1p10e11.3)')
     &         ( ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), ix = 1, nx ),
     &             iy = ny, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# x          y        ",
     &                      "  flux       r.err")')

               do iy = 1, ny
               do ix = 1, nx

                  write(iot,'(1p3e11.3,0pf8.4)')
     &               xm(ix)  + rtxdl(m)/2.0,
     &               ym(iy)  + rtydl(m)/2.0,
     &               tr(ip,ie,ia,it,icf(ix,iy,iz),im,1),
     &               tr(ip,ie,ia,it,icf(ix,iy,iz),im,2)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

                  write(iot,'( "#   x = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7,/
     &                         "#   y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            xm(1) + rtxdl(m)/2.0, xm(nx) + rtxdl(m)/2.0, rtxdl(m),
     &            ym(1) + rtydl(m)/2.0, ym(ny) + rtydl(m)/2.0, rtydl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'y/x',( xm(ix) + rtxdl(m)/2.0, ix = 1, nx )

               do iy = ny, 1, -1

                  write(iot,'(1p1000e11.3)')
     &            ym(iy) + rtydl(m)/2.0,
     &            ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), ix = 1, nx )

               end do

            end if

*-----------------------------------------------------------------------
*        gshow
*-----------------------------------------------------------------------

         if( itgsh(m) .ne. 0 .and.
     &     ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) ) then

               write(iot,'(/"#",78("-"))')
               write(iot,'( "# gshow")')
               write(iot,'( "#",78("-"))')
               write(iot,'(/"p: legs[c5*0.875]")')

               zval = ( zm(iz) + zm(iz+1) ) / 2.0d0
               none = 1
               iaxs = 1
               iuni = itgsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nx+1,ny+1,none,xm,ym,zval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

         end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.05] form[c1/0.05] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.05] form[c1/0.05] ",
     &"nosp afac[c5*0.625] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if( itazl(m) .eq. 0 ) then

         write(iot,'("y: ",a7,1x,a32)') cname, hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]"/
     &                     "  part  &=&  ",a8)')
     &                     yen, zm(iz), zm(iz+1), chq(ip)
               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+1)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+1)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+1)
                end if
               end if

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,it)
               end if

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
               end if
                     write(iot,'("e:")')

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') xmin, xmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') ymin, ymax

      end if

*-----------------------------------------------------------------------

            end do
            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        rpp-crossing
*        yz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 11 ) then

               if( ittwo(m) .eq. 1 ) then

                  dc2 = 'h2: '

               else if( ittwo(m) .eq. 2 ) then

                  dc2 = 'hd: '

               else if( ittwo(m) .eq. 3 ) then

                  dc2 = 'hc: '

               else if( ittwo(m) .eq. 6 ) then

                  dc2 = 'hd2:'

               else if( ittwo(m) .eq. 7 ) then

                  dc2 = 'hc2:'

               end if

*-----------------------------------------------------------------------

               inum = 0

            do im = 1, nm
            do ix = 1, nx
            do ip = 1, npg
            do ie = 1, neg
            do ia = 1, nag
            do it = 1, ntg

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if ( ip.gt.np_mxang ) then
                 write(iot,'( " SKIPPAGE:")')
                 inum = inum - 1
               end if

*-----------------------------------------------------------------------

                  write(iot,'("#   no. =",i3,3x,
     &            "ix =",i3/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ix, xm(ix), xm(ix+1)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+1)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+1)
               end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, ix, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,a1)')
     &                     cha, inum, ix, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, it, cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, ix, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, ix, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, ix, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, ie, ia, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ix =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ix, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: y [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = ( ym(ny+1) - ym(1) )
     &                  / ( zm(nz+1) - zm(1) )
                  xfac  = 0.9
                  afac  = 0.8
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( ( ( ittwo(m) .ge. 2 .and. ittwo(m) .le. 3 ) .or.
     &               ( ittwo(m) .ge. 6 .and. ittwo(m) .le. 7 ) ) .and.
     &             cmin .gt. 0.0 .and. cmax .gt. cmin ) then

                if (ioe .eq. 1 ) then
                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax
                else
                  write(iot,'( "set: c3[1.0e-4] c4[1.0]")')
                end if
                  write(iot,'( "p: cmin[c3] cmax[c4]")')
                  write(iot,'( "p: dmin(1e-31)")')

                  if( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &                 write(iot,'( "p: zlog")')

               end if

*-----------------------------------------------------------------------

               zmin = zm(1)
               zmax = zm(nz+1)
               ymin = ym(1)
               ymax = ym(ny+1)

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') ymin, ymax

               if( itanl(m) .gt. 0 .and. ioe .eq. 1 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               else if( itanl(m) .gt. 0 .and. ioe .ne. 1 ) then

                   call terrang(iot,m)

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

            end if

*-----------------------------------------------------------------------

                  write(iot,'("#  ny = ",i3,"   nz = ",i3)')
     &                         ny, nz

            if( ittwo(m) .ne. 4 ) then

                  write(iot,'( "# ( ( data(z,y), z = 1, nz ),",
     &                         " y = ny, 1, -1 )")')

            end if

*-----------------------------------------------------------------------

            if( ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) .and.
     &            igsh .eq. 0 ) then

                  write(iot,'(/a4," y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7," ; x = ",
     &            1p1g14.7," to ",1p1g14.7," by ",
     &            1p1g14.7," ;")') dc2,
     &            ym(ny) + rtydl(m)/2.0, ym(1) + rtydl(m)/2.0, rtydl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

               write(iot,'(1p10e11.3)')
     &         ( ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), iz = 1, nz ),
     &             iy = ny, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# y          z        ",
     &                      "  flux       r.err")')

               do iz = 1, nz
               do iy = 1, ny

                  write(iot,'(1p3e11.3,0pf8.4)')
     &               ym(iy)  + rtydl(m)/2.0,
     &               zm(iz)  + rtzdl(m)/2.0,
     &               tr(ip,ie,ia,it,icf(ix,iy,iz),im,1),
     &               tr(ip,ie,ia,it,icf(ix,iy,iz),im,2)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

                  write(iot,'( "#   y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7,/
     &                         "#   z = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            ym(1) + rtydl(m)/2.0, ym(ny) + rtydl(m)/2.0, rtydl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'y/z',( zm(iz) + rtzdl(m)/2.0, iz = 1, nz )

               do iy = ny, 1, -1

                  write(iot,'(1p1000e11.3)')
     &            ym(iy) + rtydl(m)/2.0,
     &            ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), iz = 1, nz )

               end do

            end if

*-----------------------------------------------------------------------
*        gshow
*-----------------------------------------------------------------------

         if( itgsh(m) .ne. 0 .and.
     &     ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) ) then

               write(iot,'(/"#",78("-"))')
               write(iot,'( "# gshow")')
               write(iot,'( "#",78("-"))')
               write(iot,'(/"p: legs[c5*0.875]")')

               xval = ( xm(ix) + xm(ix+1) ) / 2.0d0
               none = 1
               iaxs = 2
               iuni = itgsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,ny+1,none,zm,ym,xval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

         end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.05] form[c1/0.05] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.05] form[c1/0.05] ",
     &"nosp afac[c5*0.625] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if( itazl(m) .eq. 0 ) then

         write(iot,'("y: ",a7,1x,a32)') cname, hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  part  &=&  ",a8)')
     &                     yen, xm(ix), xm(ix+1), chq(ip)
               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+1)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+1)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+1)
                end if
               end if

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,it)
               end if

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
               end if
                     write(iot,'("e:")')

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') ymin, ymax

      end if

*-----------------------------------------------------------------------

            end do
            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        rpp-crossing
*        xz axis ( matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 12 ) then

               if( ittwo(m) .eq. 1 ) then

                  dc2 = 'h2: '

               else if( ittwo(m) .eq. 2 ) then

                  dc2 = 'hd: '

               else if( ittwo(m) .eq. 3 ) then

                  dc2 = 'hc: '

               else if( ittwo(m) .eq. 6 ) then

                  dc2 = 'hd2:'

               else if( ittwo(m) .eq. 7 ) then

                  dc2 = 'hc2:'

               end if

*-----------------------------------------------------------------------

               inum = 0

            do im = 1, nm
            do iy = 1, ny
            do ip = 1, npg
            do ie = 1, neg
            do ia = 1, nag
            do it = 1, ntg

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if ( ip.gt.np_mxang ) then
                 write(iot,'( " SKIPPAGE:")')
                 inum = inum - 1
               end if

*-----------------------------------------------------------------------

                  write(iot,'("#   no. =",i3,3x,
     &            "iy =",i3/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iy, ym(iy), ym(iy+1)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+1)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+1)
               end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, iy, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,",  it =",i3,a1)')
     &                     cha, inum, iy, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, iy, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, iy, it, cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iy, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iy, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, iy, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, iy, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, ie, ia, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iy =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iy, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: x [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = ( xm(nx+1) - xm(1) )
     &                  / ( zm(nz+1) - zm(1) )
                  xfac  = 0.9
                  afac  = 0.8
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( ( ( ittwo(m) .ge. 2 .and. ittwo(m) .le. 3 ) .or.
     &               ( ittwo(m) .ge. 6 .and. ittwo(m) .le. 7 ) ) .and.
     &             cmin .gt. 0.0 .and. cmax .gt. cmin ) then

                if (ioe .eq. 1 ) then
                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax
                else
                  write(iot,'( "set: c3[1.0e-4] c4[1.0]")')
                end if
                  write(iot,'( "p: cmin[c3] cmax[c4]")')
                  write(iot,'( "p: dmin(1e-31)")')

                  if( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &                 write(iot,'( "p: zlog")')

               end if

*-----------------------------------------------------------------------

               zmin = zm(1)
               zmax = zm(nz+1)
               xmin = xm(1)
               xmax = xm(nx+1)

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') xmin, xmax

               if( itanl(m) .gt. 0 .and. ioe .eq. 1 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               else if( itanl(m) .gt. 0 .and. ioe .ne. 1 ) then

                   call terrang(iot,m)

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

            end if

*-----------------------------------------------------------------------

                  write(iot,'("#  nx = ",i3,"   nz = ",i3)')
     &                         nx, nz

            if( ittwo(m) .ne. 4 ) then

                  write(iot,'( "# ( ( data(z,x), z = 1, nz ),",
     &                         " x = nx, 1, -1 )")')

            end if

*-----------------------------------------------------------------------

            if( ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) .and.
     &            igsh .eq. 0 ) then

                  write(iot,'(/a4," y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7," ; x = ",
     &            1p1g14.7," to ",1p1g14.7," by ",
     &            1p1g14.7," ;")') dc2,
     &            xm(nx) + rtxdl(m)/2.0, xm(1) + rtxdl(m)/2.0, rtxdl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

               write(iot,'(1p10e11.3)')
     &         ( ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), iz = 1, nz ),
     &             ix = nx, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# x          z        ",
     &                      "  flux       r.err")')

               do iz = 1, nz
               do ix = 1, nx

                  write(iot,'(1p3e11.3,0pf8.4)')
     &               xm(ix)  + rtxdl(m)/2.0,
     &               zm(iz)  + rtzdl(m)/2.0,
     &               tr(ip,ie,ia,it,icf(ix,iy,iz),im,1),
     &               tr(ip,ie,ia,it,icf(ix,iy,iz),im,2)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

                  write(iot,'( "#   y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7,/
     &                         "#   z = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            xm(1) + rtxdl(m)/2.0, xm(nx) + rtxdl(m)/2.0, rtxdl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'x/z',( zm(iz) + rtzdl(m)/2.0, iz = 1, nz )

               do ix = nx, 1, -1

                  write(iot,'(1p1000e11.3)')
     &            xm(ix) + rtxdl(m)/2.0,
     &            ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), iz = 1, nz )

               end do

            end if

*-----------------------------------------------------------------------
*        gshow
*-----------------------------------------------------------------------

         if( itgsh(m) .ne. 0 .and.
     &     ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) ) then

               write(iot,'(/"#",78("-"))')
               write(iot,'( "# gshow")')
               write(iot,'( "#",78("-"))')
               write(iot,'(/"p: legs[c5*0.875]")')

               yval = ( ym(iy) + ym(iy+1) ) / 2.0d0
               none = 1
               iaxs = 2
               iuni = itgsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,nx+1,none,zm,xm,yval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

         end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.05] form[c1/0.05] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.05] form[c1/0.05] ",
     &"nosp afac[c5*0.625] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if( itazl(m) .eq. 0 ) then

         write(iot,'("y: ",a7,1x,a32)') cname, hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  part  &=&  ",a8)')
     &                     yen, ym(iy), ym(iy+1), chq(ip)
               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+1)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+1)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+1)
                end if
               end if

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,it)
               end if

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
               end if
                     write(iot,'("e:")')

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') xmin, xmax

      end if

*-----------------------------------------------------------------------

            end do
            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------

         end if          ! itaxis

*-----------------------------------------------------------------------

            call prestart(m,iot) !OBINATA(2012.7.12)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

  900 continue

*-----------------------------------------------------------------------
      deallocate( ixyz )
      include 'samepage_include/samepage999.inc'

      deallocate (ax_x,ax_y,ax_z)

      return
      end


************************************************************************
*                                                                      *
      subroutine psufrz_rcc(m,np,nr,nz,ne,na,nt,nm,rm,zm,eb,ab,tb,
     &                      tr,tz,idasa)
*                                                                      *
*       output r-z scoring mesh current & flux for enclosure mode      *
*       last modified by S.Abe on 2018/10/18                           *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)
      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80
      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall21/ rtfac(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /tall73/ itenclo(itlmax), itangform(itlmax)

      common /istcut/ ist_cut, ist_bat

      common /fact01/ facmax(itlmax) ! kitamura23/03/31

*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      common /talout/ itall

      character fname*100, fnume*3

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /tall50/ itlmt(itlmax), ite2l(itlmax)
*-----------------------------------------------------------------------

      dimension   rm(nr+1)
      dimension   zm(nz+1)
      dimension   eb(ne+1)
      dimension   ab(na+1)
      dimension   tb(nt+1)
      dimension   ew(ne)
      dimension   aw(na)
      dimension   tw(nt)

cKN 2015/10/30 ! 2023/2/17 Ogawa. (nz+1)th and (nr+1)th are not used in enclosure mode but for the sake of coherence with tsufrz, +1 remains
      dimension   tr(np,ne,na,nt,(nr+1)*nz,nm,2)
      dimension   tz(np,ne,na,nt,nr*(nz+1),nm,2)


*-----------------------------------------------------------------------

      dimension tott(np,2) ! frtati 2021/10/05 6 -> np

      character hsunit(16)*32

      data hsunit( 1) / '[1/cm^2/source]                 '/
      data hsunit( 2) / '[1/cm^2/MeV/source]             '/
      data hsunit( 3) / '[1/cm^2/Lethargy/source]        '/
      data hsunit( 4) / '[1/cm^2/sr/source]              '/
      data hsunit( 5) / '[1/cm^2/MeV/sr/source]          '/
      data hsunit( 6) / '[1/cm^2/Lethargy/sr/source]     '/
      data hsunit(11) / '[1/cm^2/nsec/source]            '/
      data hsunit(12) / '[1/cm^2/MeV/nsec/source]        '/
      data hsunit(13) / '[1/cm^2/Lethargy/nsec/source]   '/
      data hsunit(14) / '[1/cm^2/sr/nsec/source]         '/
      data hsunit(15) / '[1/cm^2/MeV/sr/nsec/source]     '/
      data hsunit(16) / '[1/cm^2/Lethargy/sr/nsec/source]'/

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      character cha*1
      data cha /"'"/

      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31

*-----------------------------------------------------------------------

      character rpa*1
      data rpa /'}'/

      character cname*7
      character dname*7
      character aname*3

      character yen*1

      character dc2*4

      real(8),allocatable :: ar_r(:),ar_z(:)

      include 'samepage_include/samepage000.inc'

*-----------------------------------------------------------------------

         irf(ir,iz) = ir + ( iz - 1 ) * ( nr + 1 )

*-----------------------------------------------------------------------
*        set mesh area ... two circules and two side walls
*-----------------------------------------------------------------------

         ar(ir,iz) = 2.d0 * pi * ( rm(ir+1)**2 - rm(ir)**2 )
     &             + 2.d0 * pi * ( rm(ir+1) + rm(ir) )
     &               * ( zm(iz+1) - zm(iz) )

*-----------------------------------------------------------------------

      include 'samepage_include/samepage001.inc'

      allocate (ar_r(nz),ar_z(nr))
      ar_r(:) = 0.0d0
      ar_z(:) = 0.0d0
      do iz=1,nz
        do ir=1,nr
          ar_r(iz) = ar_r(iz) + ar(ir,iz)
          ar_z(ir) = ar_z(ir) + ar(ir,iz)
        enddo
      enddo

      yen  = char(92)
      igsh = 0

*-----------------------------------------------------------------------

      if ( iMeVperu.eq. 1 ) then
         hsunit( 2) = '[1/cm^2/(MeV/n)/source]         '
         hsunit( 5) = '[1/cm^2/(MeV/n)/sr/source]      '
         hsunit(12) = '[1/cm^2/(MeV/n)/nsec/source]    '
         hsunit(15) = '[1/cm^2/(MeV/n)/sr/nsec/source] '
      end if

      if ( ite2l(m) .eq. 1 ) then   ! convert to energy to LET
         hsunit( 2) = '[1/cm^2/(keV/um)/source]        '
         hsunit( 5) = '[1/cm^2/(keV/um)/sr/source]     '
         hsunit(12) = '[1/cm^2/(keV/um)/nsec/source]   '
         hsunit(15) = '[1/cm^2/(keV/um)/sr/nsec/source]'
      end if
*-----------------------------------------------------------------------
      np_mxang = np
      if ( itmxang(m).ne.0 .and. iabs(itmxang(m)).lt.np ) then
        np_mxang = iabs(itmxang(m))
      end if

*-----------------------------------------------------------------------
*     for igshow
*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

            npg = np
            neg = ne
            nag = na
            ntg = nt

         else

            npg = 1
            neg = 1
            nag = 1
            ntg = 1

         end if

*-----------------------------------------------------------------------
*        current or flux
*-----------------------------------------------------------------------

            if( itout(m) .eq. 1 ) then

               cname = '  Flux '
               dname = '  flux '

            else if( itout(m) .eq. 2 ) then

               cname = 'Current'
               dname = 'current'

            else if( itout(m) .eq. 3 ) then

               cname = ' F-Curr'
               dname = ' f-curr'

            else if( itout(m) .eq. 4 ) then

               cname = ' B-Curr'
               dname = ' b-curr'

            else if( itout(m) .eq. 5 ) then

               cname = ' O-Curr'
               dname = ' o-curr'

            else if( itout(m) .eq. 6 ) then

               cname = 'OF-Curr'
               dname = 'of-curr'

            else if( itout(m) .eq. 7 ) then

               cname = 'OB-Curr'
               dname = 'ob-curr'

            else if( itout(m) .eq. 8 ) then

               cname = ' A-Curr'
               dname = ' a-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            else if( itout(m) .eq. 9 ) then

               cname = 'OA-Curr'
               dname = 'oa-curr'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            else if( itout(m) .eq. 10 ) then ! T.Sato 2021/05/05

               cname = ' A-Flux'
               dname = ' a-flux'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            else if( itout(m) .eq. 11 ) then ! T.Sato 2021/05/05

               cname = 'OA-Flux'
               dname = 'oa-flux'

               if( itaty(m) .gt. 0 ) then
                  aname = 'cos'
               else
                  aname = 'the'
               end if

            end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               call wtpname(m,np,chp,chq)

*-----------------------------------------------------------------------
*        c1 : nomalization for source
*-----------------------------------------------------------------------

            if ( rsouin .gt. 0d0 ) then ! S.H. (2016.7.26)

               c1 = 1.0d+0 / rsouin

            else

               c1 = 0d0

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 1 : /cm^2/source
*                    = 2 : /cm^2/MeV/source
*                    = 3 : /cm^2/Lethargy/source
*                    = 4 : /cm^2/source/SR
*                    = 5 : /cm^2/MeV/source/SR
*                    = 6 : /cm^2/Lethargy/source/SR
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  4 .or.
     &          itunt(m) .eq. 11 .or. itunt(m) .eq. 14 ) then

               do i = 1, ne

                  ew(i) = 1.d+0

               end do
               ew_sum = 1.0d0

            else if( itunt(m) .eq.  2 .or. itunt(m) .eq.  5 .or.
     &               itunt(m) .eq. 12 .or. itunt(m) .eq. 15 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &               itunt(m) .eq. 13 .or. itunt(m) .eq. 16 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

            end if

*-----------------------------------------------------------------------

            aw_sum = 0.0d0
            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  2 .or.
     &          itunt(m) .eq.  3 .or. itunt(m) .eq. 11 .or.
     &          itunt(m) .eq. 12 .or. itunt(m) .eq. 13 ) then

               do i = 1, na

                  aw(i) = 1.d+0

               end do
               aw_sum = 1.0d0

            else

               do i = 1, na

                  if( itaty(m) .gt. 0 ) then

                     aw(i) = ( ab(i+1) - ab(i) ) * 2.0 * pi

                  else

                     aw(i) = ( cos( ab(i) / 180.d0 * pi )
     &                       - cos( ab(i+1) / 180.d0 * pi ) ) * 2.0 * pi

                  end if
                  aw_sum = aw_sum + aw(i)

               end do

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 11, 12, 13, 14 : /.../nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .gt. 10 ) then

               do i = 1, nt

                  tw(i) = tb(i+1) - tb(i)

               end do
               tw_sum = tb(nt+1) - tb(1)

            else

               do i = 1, nt

                  tw(i) = 1.d+0

               end do
               tw_sum = 1.0d0

            end if

*-----------------------------------------------------------------------
*        relative error and  unit conversion
*-----------------------------------------------------------------------
*        rcc-crossing
*-----------------------------------------------------------------------

                  cmax = 0.0
                  cmin = 1.e+33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

            facmax(m) = 1.d0
            if( rtfac(m) .lt. 0.d0 ) then
               facmax(m) = 0.d0

               do 201 im = 1, nm
               do 201 iz = 1, nz
               do 201 ir = 1, nr
               do 201 it = 1, nt
               do 201 ia = 1, na
               do 201 ie = 1, ne
               do 201 ip = 1, np

                  if( tr(ip,ie,ia,it,irf(ir,iz),im,1) .gt. 0.d0 ) then

                     fmaxfc = tr(ip,ie,ia,it,irf(ir,iz),im,1)
     &                       / ar(ir,iz) / ew(ie) / aw(ia) / tw(it)

                    if( facmax(m) .lt. fmaxfc) facmax(m) = fmaxfc

                  end if

  201          continue

            end if

            do 200 im = 1, nm
            do 200 iz = 1, nz
            do 200 ir = 1, nr
            do 200 it = 1, nt
            do 200 ia = 1, na
            do 200 ie = 1, ne
            do 200 ip = 1, np

               if( tr(ip,ie,ia,it,irf(ir,iz),im,1) .gt. 0.d0 ) then

                call calc_stdev(m,Xa,sigx,
     &                       tr(ip,ie,ia,it,irf(ir,iz),im,1),
     &                       tr(ip,ie,ia,it,irf(ir,iz),im,2),
     &        abs(rtfac(m)/facmax(m))/ar(ir,iz)/ew(ie)/aw(ia)/tw(it))

                  tr(ip,ie,ia,it,irf(ir,iz),im,1) = Xa
                  tr(ip,ie,ia,it,irf(ir,iz),im,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,ia,it,irf(ir,iz),im,1) .gt. cmax )
     &                     cmax = tr(ip,ie,ia,it,irf(ir,iz),im,1)

                  if( tr(ip,ie,ia,it,irf(ir,iz),im,1) .lt. cmin )
     &                     cmin = tr(ip,ie,ia,it,irf(ir,iz),im,1)

               else

                  isdz = 1
                  tr(ip,ie,ia,it,irf(ir,iz),im,2) = 0.0

               end if

! sumover
               fact_in = abs(rtfac(m)/facmax(m))
               call psufrz_sumover_stdev(0,m,ip,ie,ia,it,ir,iz,im,
     &                fact_in,ew(ie),aw(ia),tw(it),ar(ir,iz),
     &                ew_sum,aw_sum,tw_sum,ar_r(iz),ar_z(ir))


  200       continue

            if( nobch .gt. ist_bat ) then
            if( rtstd(m) .gt. 0.d0 .and.
     &          stdm .gt. 0.d0 .and. isdz .eq. 0 .and.
     &          stdm .lt. rtstd(m) ) itstd(m) = 0
            end if

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do iax = 1, itfln(m)

!OBINATA(2012.7.12): output *.err
         noe = 1
         if ( any( itaxs(m,iax) .eq. (/ 7 /) )
     &       .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &     noe = 2

         do ioe = 1, noe

         if( itall .eq. 2 .and. nobch .lt. maxbch ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.7.12): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)(1:itfll(m,iax))//'.'//fnume
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         else

!OBINATA(2012.7.12): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         end if

            iot = 31
            open(iot, file = fname, status = 'unknown' )

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

               call tcrsech(iot,m,iax,1)

*-----------------------------------------------------------------------
*        rcc-crossing
*        energy axis  ( LET axis )
*-----------------------------------------------------------------------
      include 'samepage_include/samepage002_pmeatrz.inc'
      include 'samepage_include/samepagechp_pmeatrz.inc'
      include 'samepage_include/samepageseti.inc'

            nrstepi_0 = nrstepi
            nzstepi_0 = nzstepi

               inum = 0

         if( itaxs(m,iax) .eq. 1 .or.
     &       itaxs(m,iax) .eq. 14 ) then

            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iri = 1, nr, nrstepi
            do izi = 1, nz, nzstepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ir = iri
               iz = izi
               ia = iai
               it = iti
               ip = ipi


               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif


                  write(iot,'("#   no. =",i3,3x,
     &            "ir  =",i3,3x,
     &            "iz  =",i3/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ir, iz,
     &                        rm(ir), rm(ir+nrstepi),
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#    ia =",i3/
     &            "#   ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#    it =",i3/
     &            "#     t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                if ( itaxs(m,iax) .eq. 1 ) then   ! ccse 2022/09/30

                  write(iot,'(/"x: Energy [MeV]")')

                else if ( itaxs(m,iax) .eq. 14 ) then
                 write(iot,'(/"x: LET [keV/um]")')
                end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 13 .or. itunt(m) .eq. 16 .or.
     &             itety(m) .eq.  3 .or. itety(m) .eq.  5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatrzm_e_rcc.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ir, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, iz, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,a1)')
     &                     cha, inum, ir, iz, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, iz, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, it, itmnt(m,im), cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if
             else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                  changelsub(3) = " "
                else
                  write(changelsub(3),'(",  mset ="i3)') itmnt(m,im)
                end if
               if( na .le. 1 ) then
                  changelsub(4) = " "
                else
                  write(changelsub(4),'(",  ia =",i3)') ia
                end if
                write(changelsub(5),'(",  (ir,iz) = (",i3,",",i3,")")')
     &                                                 ir, iz
                if(iloopmode .eq. 5) then
                   write(changelsub(5),'(",  ir =",i3)') ir
                end if
                if(iloopmode .eq. 6) then
                   write(changelsub(5),'(",  iz =",i3)') iz
                end if
                write(changelsub(6),'(",  ie =",i3)') ie
                if( ittty(m) .eq. 0 ) then
                  changelsub(7) = " "
                else
                  write(changelsub(7),'(",  it =",i3)') it
                end if
                write(changelsub(8),'(a1)') cha
c
                if( itout(m) .le. 4 ) then
                  changelsub(4) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(4) = " "
                  changelsub(6) = " "
                else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then
                else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then
                  changelsub(6) = " "
                end if

                changelsub(6) = " "
c
                if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                  changelsub(4) = " "
                end if
                if( ittty(m) .eq. 0 ) then
                  changelsub(7) = " "
                end if
                if(iloopmode .eq. 1) then
                  changelsub(6) = " "
                end if
                if(iloopmode .eq. 10) then
                  changelsub(7) = " "
                end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))
                write(iot,'(/a)') trim(angeltitle)
             end if
*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------
                 areasum = 0.0d0
                 do ii=izi,izi+nzstepi-1
                   do i=iri,iri+nrstepi-1
                     areasum = areasum + ar(i,ii)
                   end do
                 end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rcc surface crossing"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &       areasum, rm(ir), rm(ir+nrstepi), zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        rcc-crossing
*        r axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 6 ) then

            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iei = 1, ne, nestepi
            do izi = 1, nz, nzstepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ie = iei
               iz = izi
               ia = iai
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

*-----------------------------------------------------------------------

               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif


                  write(iot,'("#   no. =",i3,3x,
     &            "iz  =",i3/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, iz,
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05

                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: r [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itzty(m) .eq. 3 .or. itzty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatrzm_r_rcc.inc'
      voll = samewtt(1)

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ie, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ie, iz, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, iz, it, cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ie, ia, iz, cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ie, ia, iz, it, cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ia, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ia, iz, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, iz, itmnt(m,im), cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, iz, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, iz, it, itmnt(m,im), cha
                  end if
               end if

            end if
             else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                  changelsub(3) = " "
                else
                  write(changelsub(3),'(",  mset ="i3)') itmnt(m,im)
                end if
               if( na .le. 1 ) then
                  changelsub(4) = " "
                else
                  write(changelsub(4),'(",  ia =",i3)') ia
                end if
                write(changelsub(5),'(",  (ir,iz) = (",i3,",",i3,")")')
     &                                                 ir, iz
                if(iloopmode .eq. 5) then
                   write(changelsub(5),'(",  ir =",i3)') ir
                end if
                if(iloopmode .eq. 6) then
                   write(changelsub(5),'(",  iz =",i3)') iz
                end if
                write(changelsub(6),'(",  ie =",i3)') ie
                if( ittty(m) .eq. 0 ) then
                  changelsub(7) = " "
                else
                  write(changelsub(7),'(",  it =",i3)') it
                end if
                write(changelsub(8),'(a1)') cha
c
                if( itout(m) .le. 4 ) then
                  changelsub(4) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(4) = " "
                  changelsub(6) = " "
                else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then
                else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then
                  changelsub(6) = " "
                end if

                if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                  changelsub(4) = " "
                end if
                if( ittty(m) .eq. 0 ) then
                  changelsub(7) = " "
                end if
                if(iloopmode .eq. 1) then
                  changelsub(6) = " "
                end if
                if(iloopmode .eq. 10) then
                  changelsub(7) = " "
                end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))
                write(iot,'(/a)') trim(angeltitle)

              end if
*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------


                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rcc surface crossing"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05


                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+nestepi)
                end if
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'(
     &                     " tot area  &=&",1pe13.4," [cm^2]")')
     &                     voll

                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        rcc-crossing
*        z axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 5 ) then

            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iri = 1, nr, nrstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ir = iri
               ie = iei
               ia = iai
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ir  =",i3/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ir,
     &                        rm(ir), rm(ir+nrstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05

                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                          ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itrty(m) .eq. 3 .or. itrty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatrzm_z_rcc.inc'
      voll = samewtt(1)

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",   ie =",i3,a1)')
     &                     cha, inum, ir, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",   ie =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",    ir =",i3,a1)')
     &                     cha, inum, ir, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",    ir =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, it, cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",   ie =",i3,
     &                     ",   ia =",i3,a1)')
     &                     cha, inum, ir, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",   ie =",i3,
     &                     ",   ia =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",   ia =",i3,a1)')
     &                     cha, inum, ir, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",   ia =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",   ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",   ie =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",    ir =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",    ir =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",   ie =",i3,
     &                     ",   ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",   ie =",i3,
     &                     ",   ia =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, ie, ia, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",   ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",   ir =",i3,
     &                     ",   ia =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if
             else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                  changelsub(3) = " "
                else
                  write(changelsub(3),'(",  mset ="i3)') itmnt(m,im)
                end if
               if( na .le. 1 ) then
                  changelsub(4) = " "
                else
                  write(changelsub(4),'(",  ia =",i3)') ia
                end if
                write(changelsub(5),'(",  (ir,iz) = (",i3,",",i3,")")')
     &                                                 ir, iz
                if(iloopmode .eq. 5) then
                   write(changelsub(5),'(",  ir =",i3)') ir
                end if
                if(iloopmode .eq. 6) then
                   write(changelsub(5),'(",  iz =",i3)') iz
                end if
                write(changelsub(6),'(",  ie =",i3)') ie
                if( ittty(m) .eq. 0 ) then
                  changelsub(7) = " "
                else
                  write(changelsub(7),'(",  it =",i3)') it
                end if
                write(changelsub(8),'(a1)') cha
c
                if( itout(m) .le. 4 ) then
                  changelsub(4) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(4) = " "
                  changelsub(6) = " "
                else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then
                else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then
                  changelsub(6) = " "
                end if

                if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                  changelsub(4) = " "
                end if
                if( ittty(m) .eq. 0 ) then
                  changelsub(7) = " "
                end if
                if(iloopmode .eq. 1) then
                  changelsub(6) = " "
                end if
                if(iloopmode .eq. 10) then
                  changelsub(7) = " "
                end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))
                write(iot,'(/a)') trim(angeltitle)
              end if
*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rcc surface crossing"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]")')
     &                     yen, rm(ir), rm(ir+nrstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05


                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+1)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+nestepi)
                end if
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'(
     &                     " tot area  &=&",1pe13.4," [cm^2]")')
     &                     voll

                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        rcc-crossing
*        angle axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 .or.
     &            itaxs(m,iax) .eq. 10 ) then

            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iri = 1, nr, nrstepi
            do izi = 1, nz, nzstepi
            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ir = iri
               iz = izi
               ie = iei
               it = iti
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)
*-----------------------------------------------------------------------
               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ir  =",i3,3x,
     &            "iz  =",i3/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ir, iz,
     &                        rm(ir), rm(ir+nrstepi),
     &                        zm(iz), zm(iz+nzstepi)


               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)



               if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  if( itaxs(m,iax) .eq. 8 ) then
                     write(iot,'(/"x: cos(",a1,"theta)")') yen
                  else
                     write(iot,'(/"x: ",a1,"theta  [deg]")') yen
                  end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatrzm_a_rcc.inc'

*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,a1)')
     &                     cha, inum, ir, iz, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, iz, ie, it, cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ir, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,",   it =",i3,a1)')
     &                     cha, inum, ir, iz, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,",   it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, it, itmnt(m,im), cha
                  end if
               end if

            end if
             else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                  changelsub(3) = " "
                else
                  write(changelsub(3),'(",  mset ="i3)') itmnt(m,im)
                end if
               if( na .le. 1 ) then
                  changelsub(4) = " "
                else
                  write(changelsub(4),'(",  ia =",i3)') ia
                end if
                write(changelsub(5),'(",  (ir,iz) = (",i3,",",i3,")")')
     &                                                 ir, iz
                if(iloopmode .eq. 5) then
                   write(changelsub(5),'(",  ir =",i3)') ir
                end if
                if(iloopmode .eq. 6) then
                   write(changelsub(5),'(",  iz =",i3)') iz
                end if
                write(changelsub(6),'(",  ie =",i3)') ie
                if( ittty(m) .eq. 0 ) then
                  changelsub(7) = " "
                else
                  write(changelsub(7),'(",  it =",i3)') it
                end if
                write(changelsub(8),'(a1)') cha
c
                if( itout(m) .le. 4 ) then
                  changelsub(4) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(4) = " "
                  changelsub(6) = " "
                else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then
                else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then
                  changelsub(6) = " "
                end if

                changelsub(4) = " "
c
                if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                  changelsub(4) = " "
                end if
                if( ittty(m) .eq. 0 ) then
                  changelsub(7) = " "
                end if
                if(iloopmode .eq. 1) then
                  changelsub(6) = " "
                end if
                if(iloopmode .eq. 10) then
                  changelsub(7) = " "
                end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))
                write(iot,'(/a)') trim(angeltitle)

              end if
*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

                 areasum = 0.0d0
                 do ii=izi,izi+nzstepi-1
                   do i=iri,iri+nrstepi-1
                     areasum = areasum + ar(i,ii)
                   end do
                 end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rcc surface crossing"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &      areasum, rm(ir), rm(ir+nrstepi), zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4,/
     &                     "  emax  &=&",1pe13.4)')
     &                     eb(ie), eb(ie+nestepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if
                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        r-crossing
*        time axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9 ) then

            if(nzstepi .gt. nz) then
               nzstepi = nz
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            if(nrstepi .gt. nr) then
               nrstepi = nr
               nsame = max(nrstepi, nzstepi, nastepi,
     &                  ntstepi, npstepi, nestepi)
               np_mxang =min(np_mxang, nsame)
            end if
            nmstepi = 1
            do imi = 1, nm, nmstepi
            do iri = 1, nr, nrstepi
            do izi = 1, nz, nzstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
               im = imi
               ir = iri
               iz = izi
               ie = iei
               ia = iai
               ip = ipi

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if(iloopmode .ne. 0) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

*-----------------------------------------------------------------------

                  write(iot,'("#   no. =",i3,3x,
     &            "ir  =",i3,3x,
     &            "iz  =",i3/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ir, iz,
     &                        rm(ir), rm(ir+nrstepi),
     &                        zm(iz), zm(iz+nzstepi)

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05

                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Time [nsec]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: ",a7,1x,a32)')
     &                            cname, hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( ittty(m) .eq. 3 .or. ittty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

      include 'samepage_include/peatrzm_t_rcc.inc'
*-----------------------------------------------------------------------

          if(iloopmode .eq. 0) then

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,a1)')
     &                     cha, inum, ir, iz, ie, cha
               else if( itout(m) .le. 7 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ir, iz, cha
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,
     &                     ",  ia = ",i3,a1)')
     &                     cha, inum, ir, iz, ie, ia, cha
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,a1)')
     &                     cha, inum, ir, iz, ia, cha
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ie, itmnt(m,im), cha
               else if( itout(m) .le. 7 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, itmnt(m,im), cha
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ie = ",i3,
     &                     ",  ia = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ie, ia, itmnt(m,im), cha
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ir =",i3,
     &                     ",  iz =",i3,
     &                     ",  ia = ",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ir, iz, ia, itmnt(m,im), cha
               end if

            end if

             else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                  changelsub(3) = " "
                else
                  write(changelsub(3),'(",  mset ="i3)') itmnt(m,im)
                end if
               if( na .le. 1 ) then
                  changelsub(4) = " "
                else
                  write(changelsub(4),'(",  ia =",i3)') ia
                end if
                write(changelsub(5),'(",  (ir,iz) = (",i3,",",i3,")")')
     &                                                 ir, iz
                if(iloopmode .eq. 5) then
                   write(changelsub(5),'(",  ir =",i3)') ir
                end if
                if(iloopmode .eq. 6) then
                   write(changelsub(5),'(",  iz =",i3)') iz
                end if
                write(changelsub(6),'(",  ie =",i3)') ie
                if( ittty(m) .eq. 0 ) then
                  changelsub(7) = " "
                else
                  write(changelsub(7),'(",  it =",i3)') it
                end if
                write(changelsub(8),'(a1)') cha
c
                if( itout(m) .le. 4 ) then
                  changelsub(4) = " "
                else if( itout(m) .le. 7 ) then
                  changelsub(4) = " "
                  changelsub(6) = " "
                else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then
                else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then
                  changelsub(6) = " "
                end if

                changelsub(7) = " "
c
                if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                  changelsub(4) = " "
                end if
                if( ittty(m) .eq. 0 ) then
                  changelsub(7) = " "
                end if
                if(iloopmode .eq. 1) then
                  changelsub(6) = " "
                end if
                if(iloopmode .eq. 10) then
                  changelsub(7) = " "
                end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))
                write(iot,'(/a)') trim(angeltitle)
               end if
*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------
                 areasum = 0.0d0
                 do ii=izi,izi+nzstepi-1
                   do i=iri,iri+nrstepi-1
                     areasum = areasum + ar(i,ii)
                   end do
                 end do

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rcc surface crossing"/
     &                     "  area  &=&",1pe13.4," [cm^2]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &     areasum, rm(ir), rm(ir+nrstepi), zm(iz), zm(iz+nzstepi)

               if(iloopmode .ne. 0) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05

                  write(iot,'(
     &                     "  emin  &=&",1pe13.4,/
     &                     "  emax  &=&",1pe13.4)')
     &                     eb(ie), eb(ie+nestepi)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

            nrstepi = nrstepi_0
            nzstepi = nzstepi_0

*-----------------------------------------------------------------------
*        rcc-crossing
*        rz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 13 ) then

               if( ittwo(m) .eq. 1 ) then

                  dc2 = 'h2: '

               else if( ittwo(m) .eq. 2 ) then

                  dc2 = 'hd: '

               else if( ittwo(m) .eq. 3 ) then

                  dc2 = 'hc: '

               else if( ittwo(m) .eq. 6 ) then

                  dc2 = 'hd2:'

               else if( ittwo(m) .eq. 7 ) then

                  dc2 = 'hc2:'

               end if

*-----------------------------------------------------------------------

               inum = 0

            do im = 1, nm
            do ip = 1, npg
            do ie = 1, neg
            do ia = 1, nag
            do it = 1, ntg

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if ( ip.gt.np_mxang ) then
                 write(iot,'( " SKIPPAGE:")')
                 inum = inum - 1
               end if

*-----------------------------------------------------------------------

               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05

                  write(iot,'(
     &            "#  ie =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ie, eb(ie), eb(ie+1)
               end if
               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &            "#  ia =",i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+1)
               end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

            if( itmlp(m) .eq. 0 ) then

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, it, cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,a1)')
     &                     cha, inum, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ", it =",i3,a1)')
     &                     cha, inum, it, cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, ie, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, ia, it, cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,a1)')
     &                     cha, inum, ia, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,",  it =",i3,a1)')
     &                     cha, inum, ia, it, cha
                  end if
               end if

*-----------------------------------------------------------------------

            else

               if( itout(m) .le. 4 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .le. 7 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 8 .or. itout(m) .eq.10 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ia, it, itmnt(m,im), cha
                  end if
               else if( itout(m) .eq. 9 .or. itout(m) .eq.11 ) then  ! T.Sato 2021/05/05
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ia =",i3,",  it =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ia, it, itmnt(m,im), cha
                  end if
               end if

            end if

*-----------------------------------------------------------------------

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: r [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = ( rm(nr+1) - rm(1) )
     &                  / ( zm(nz+1) - zm(1) )
                  xfac  = 0.9
                  afac  = 0.8
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( ( ( ittwo(m) .ge. 2 .and. ittwo(m) .le. 3 ) .or.
     &               ( ittwo(m) .ge. 6 .and. ittwo(m) .le. 7 ) ) .and.
     &             cmin .gt. 0.0 .and. cmax .gt. cmin ) then

                if (ioe .eq. 1 ) then
                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax
                else
                  write(iot,'( "set: c3[1.0e-4] c4[1.0]")')
                end if
                  write(iot,'( "p: cmin[c3] cmax[c4]")')
                  write(iot,'( "p: dmin(1e-31)")')

                  if( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &                 write(iot,'( "p: zlog")')

               end if

*-----------------------------------------------------------------------

               zmin = zm(1)
               zmax = zm(nz+1)
               rmin = rm(1)
               rmax = rm(nr+1)

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') rmin, rmax

               if( itanl(m) .gt. 0 .and. ioe .eq. 1 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               else if( itanl(m) .gt. 0 .and. ioe .ne. 1 ) then

                   call terrang(iot,m)

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

            end if

*-----------------------------------------------------------------------

                  write(iot,'("#  nr = ",i3,"   nz = ",i3)')
     &                         nr, nz

            if( ittwo(m) .ne. 4 ) then

                  write(iot,'( "# ( ( data(z,r), z = 1, nz ),",
     &                         " r = nr, 1, -1 )")')

            end if

*-----------------------------------------------------------------------

            if( ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) .and.
     &            igsh .eq. 0 ) then

                  write(iot,'(/a4," y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7," ; x = ",
     &            1p1g14.7," to ",1p1g14.7," by ",
     &            1p1g14.7," ;")') dc2,
     &            rm(nr) + rtrdl(m)/2.0, rm(1) + rtrdl(m)/2.0, rtrdl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

               write(iot,'(1p10e11.3)')
     &         ( ( tr(ip,ie,ia,it,irf(ir,iz),im,ioe),
     &              iz = 1, nz ), ir = nr, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# r          z        ",
     &                      "  flux       r.err")')

               do iz = 1, nz
               do ir = 1, nr

                  write(iot,'(1p3e11.3,0pf8.4)')
     &               rm(ir)  + rtrdl(m)/2.0,
     &               zm(iz)  + rtzdl(m)/2.0,
     &               tr(ip,ie,ia,it,irf(ir,iz),im,1),
     &               tr(ip,ie,ia,it,irf(ir,iz),im,2)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

                  write(iot,'( "#   r = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7,/
     &                         "#   z = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            rm(1) + rtrdl(m)/2.0, rm(nr) + rtrdl(m)/2.0, rtrdl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'r/z',( zm(iz) + rtzdl(m)/2.0, iz = 1, nz )

               do ir = nr, 1, -1

                  write(iot,'(1p1000e11.3)')
     &            rm(ir) + rtrdl(m)/2.0,
     &            ( tr(ip,ie,ia,it,irf(ir,iz),im,ioe), iz = 1, nz )

               end do

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.05] form[c1/0.05] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.05] form[c1/0.05] ",
     &"nosp afac[c5*0.625] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if( itazl(m) .eq. 0 ) then

         write(iot,'("y: ",a7,1x,a32)') cname, hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

                  write(iot,'(/"wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  rpp surface crossing"/
     &                     "  part  &=&  ",a8)')
     &                     yen, chq(ip)
               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq.10) then ! T.Sato 2021/05/05


                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+1)
                  else
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     eb(ie), eb(ie+1)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     eb(ie), eb(ie+1)
                end if
               end if

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,it)
               end if

               if( itout(m) .ge. 8 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
               end if
                     write(iot,'("e:")')

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') rmin, rmax

      end if

*-----------------------------------------------------------------------

            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------

            call prestart(m,iot) !OBINATA(2012.7.12)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

      end do
      end do

*-----------------------------------------------------------------------

      include 'samepage_include/samepage999.inc'

      deallocate (ar_r,ar_z)

      return
      end

************************************************************************
*                                                                      *
* sumover subroutine group                                             *
*                                                                      *
************************************************************************
*
************************************************************************
*                                                                      *
      subroutine tsufreg_sumover(m,maxcas,
     &                 np,    ne,    na,   nt,  nr,   nm,    tr0x)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0

      implicit double precision (a-h,o-z)
      dimension   tr0x(np,ne,na,nt,nr,nm)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)


C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)

      do iax=1,itaxn(m)

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if
      if(n_tr_sum > 0) then

        call tsufreg_sumover_sub(maxcas,itaxs(m,iax),
     &     np,    ne,    na,   nt,    nr,   nm,     tr0x,
     &     itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &     ittnm_sum(m,iax),itrcn_sum(m,iax),
     &     itmst_sum(m,iax),tr_sum)

        endif

      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine tsufreg_sumover_sub(maxcas,itaxs_in,
     &                np,    ne,    na,    nt,    nr,    nm,    tr0,
     &                np_sum,ne_sum,na_sum,nt_sum,nr_sum,nm_sum,tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      dimension   tr0(np,ne,na,nt,nr,nm)
      dimension   tr_sum(np_sum,ne_sum,na_sum,nt_sum,nr_sum,nm_sum,2)

      if(itaxs_in == 1 .or. itaxs_in == 14) then  ! energ
        do im = 1, nm
        do ir = 1, nr
        do it = 1, nt
        do ia = 1, na
        do ip = 1, np
          tr0_sum = 0.0d0
          do ie = 1, ne
            tr0_sum = tr0_sum + tr0(ip,ie,ia,it,ir,im) / maxcas
          end do
          tr_sum(ip,1,ia,it,ir,im,1) = tr_sum(ip,1,ia,it,ir,im,1)
     &                            + tr0_sum
          tr_sum(ip,1,ia,it,ir,im,2) = tr_sum(ip,1,ia,it,ir,im,2)
     &                            + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
!
      elseif(itaxs_in == 2) then    ! reg
        do im = 1, nm
        do it = 1, nt
        do ia = 1, na
        do ie = 1, ne
        do ip = 1, np
          tr0_sum = 0.0d0
          do ir = 1, nr
            tr0_sum = tr0_sum + tr0(ip,ie,ia,it,ir,im) / maxcas
          end do
          tr_sum(ip,ie,ia,it,1,im,1) = tr_sum(ip,ie,ia,it,1,im,1)
     &                            + tr0_sum
          tr_sum(ip,ie,ia,it,1,im,2) = tr_sum(ip,ie,ia,it,1,im,2)
     &                            + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 8 .or. itaxs_in == 10) then   ! ang
        do im = 1, nm
        do ir = 1, nr
        do it = 1, nt
        do ie = 1, ne
        do ip = 1, np
          tr0_sum = 0.0d0
          do ia = 1, na
            tr0_sum = tr0_sum + tr0(ip,ie,ia,it,ir,im) / maxcas
          end do
          tr_sum(ip,ie,1,it,ir,im,1) = tr_sum(ip,ie,1,it,ir,im,1)
     &                            + tr0_sum
          tr_sum(ip,ie,1,it,ir,im,2) = tr_sum(ip,ie,1,it,ir,im,2)
     &                            + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 9) then   ! time
        do im = 1, nm
        do ir = 1, nr
        do ia = 1, na
        do ie = 1, ne
        do ip = 1, np
          tr0_sum = 0.0d0
          do it = 1, nt
            tr0_sum = tr0_sum + tr0(ip,ie,ia,it,ir,im) / maxcas
          end do
          tr_sum(ip,ie,ia,1,ir,im,1) = tr_sum(ip,ie,ia,1,ir,im,1)
     &                            + tr0_sum
          tr_sum(ip,ie,ia,1,ir,im,2) = tr_sum(ip,ie,ia,1,ir,im,2)
     &                            + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do

       end if
       return
       end


************************************************************************
*                                                                      *
      subroutine tsufrz_sumover(m,maxcas,
     &                   np, ne, na, nt, nr, nz,  nm, tr0x, tz0x)
*                                                                      *
*     tr_sum(np,ne,na,nt,(nr+1)*nz,nm,2)                                       *
*     tz_sum(np,ne,na,nt,nr*(nz+1),nm,2)                                       *
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0

      implicit double precision (a-h,o-z)
      dimension   tr0x(np,ne,na,nt,(nr+1)*nz,nm)
      dimension   tz0x(np,ne,na,nt,nr*(nz+1),nm)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)


C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:),tz_sum(:)

      do iax=1,itaxn(m)

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$          call GET_TZ_HEAD_POINTER0_SUM(tz_sum,m,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
            call GET_TZ_HEAD_POINTER_SUM(tz_sum,m,iax)
C for nonshared_tally option
!$       end if
      if(n_tr_sum > 0) then

        call tsufrz_sumover_sub(maxcas,itaxs(m,iax),
     &     np, ne, na, nt, nr, nz, nm,  tr0x, tz0x,
     &     itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &     ittnm_sum(m,iax),itrnm_sum(m,iax),itznm_sum(m,iax),
     &     itmst_sum(m,iax),tr_sum,tz_sum)

        endif

      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine tsufrz_sumover_sub(maxcas,itaxs_in,
     &           np, ne, na, nt, nr, nz, nm,  tr0, tz0,
     &           np_sum, ne_sum, na_sum, nt_sum, nr_sum, nz_sum,
     &           nm_sum, tr_sum,tz_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      dimension   tr0(np,ne,na,nt,(nr+1)*nz,nm)
      dimension   tz0(np,ne,na,nt,nr*(nz+1),nm)
      dimension   tr_sum(np_sum,ne_sum,na_sum,nt_sum,(nr_sum+1)*nz_sum,
     &                   nm_sum,2)
      dimension   tz_sum(np_sum,ne_sum,na_sum,nt_sum,nr_sum*(nz_sum+1),
     &                   nm_sum,2)

      icf(ir,iz)     = ir + ( iz - 1 ) * (nr + 1)
      icf_sum(ir,iz) = ir + ( iz - 1 ) * (nr_sum +1)
      izf(ir,iz)     = iz + ( ir - 1 ) * (nz + 1)
      izf_sum(ir,iz) = iz + ( ir - 1 ) * (nz_sum +1)

      if(itaxs_in == 1 .or. itaxs_in == 14) then  ! energ
        do im = 1, nm
        do iz = 1, nz
        do ir = 1, nr+1
        do it = 1, nt
        do ia = 1, na
        do ip = 1, np
          tr0_sum = 0.0d0
          do ie = 1, ne
            tr0_sum = tr0_sum + tr0(ip,ie,ia,it,icf(ir,iz),im) / maxcas
          end do
          tr_sum(ip,1,ia,it,icf_sum(ir,iz),im,1) =
     &           tr_sum(ip,1,ia,it,icf_sum(ir,iz),im,1) + tr0_sum
          tr_sum(ip,1,ia,it,icf_sum(ir,iz),im,2) =
     &          tr_sum(ip,1,ia,it,icf_sum(ir,iz),im,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do

        do im = 1, nm
        do iz = 1, nz+1
        do ir = 1, nr
        do it = 1, nt
        do ia = 1, na
        do ip = 1, np
          tz0_sum = 0.0d0
          do ie = 1, ne
            tz0_sum = tz0_sum + tz0(ip,ie,ia,it,izf(ir,iz),im) / maxcas
          end do
          tz_sum(ip,1,ia,it,izf_sum(ir,iz),im,1) =
     &           tz_sum(ip,1,ia,it,izf_sum(ir,iz),im,1) + tz0_sum
          tz_sum(ip,1,ia,it,izf_sum(ir,iz),im,2) =
     &           tz_sum(ip,1,ia,it,izf_sum(ir,iz),im,2) + tz0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do
!
      elseif(itaxs_in == 6) then    ! r
        do im = 1, nm
        do iz = 1, nz
        do it = 1, nt
        do ia = 1, na
        do ie = 1, ne
        do ip = 1, np
          tr0_sum = 0.0d0
          do ir = 1, nr+1
            tr0_sum = tr0_sum + tr0(ip,ie,ia,it,icf(ir,iz),im) / maxcas
          end do
          tr_sum(ip,ie,ia,it,icf_sum(1,iz),im,1) =
     &           tr_sum(ip,ie,ia,it,icf_sum(1,iz),im,1) + tr0_sum
          tr_sum(ip,ie,ia,it,icf_sum(1,iz),im,2) =
     &           tr_sum(ip,ie,ia,it,icf_sum(1,iz),im,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do

        do im = 1, nm
        do iz = 1, nz+1
        do it = 1, nt
        do ia = 1, na
        do ie = 1, ne
        do ip = 1, np
          tz0_sum = 0.0d0
          do ir = 1, nr
            tz0_sum = tz0_sum + tz0(ip,ie,ia,it,izf(ir,iz),im) / maxcas
          end do
          tz_sum(ip,ie,ia,it,izf_sum(1,iz),im,1) =
     &           tz_sum(ip,ie,ia,it,izf_sum(1,iz),im,1) + tz0_sum
          tz_sum(ip,ie,ia,it,izf_sum(1,iz),im,2) =
     &           tz_sum(ip,ie,ia,it,izf_sum(1,iz),im,2) + tz0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do


      elseif(itaxs_in == 5) then    ! z
        do im = 1, nm
        do ir = 1, nr+1
        do it = 1, nt
        do ia = 1, na
        do ie = 1, ne
        do ip = 1, np
          tr0_sum = 0.0d0
          do iz = 1, nz
            tr0_sum = tr0_sum + tr0(ip,ie,ia,it,icf(ir,iz),im) / maxcas
          end do
          tr_sum(ip,ie,ia,it,icf_sum(ir,1),im,1) =
     &           tr_sum(ip,ie,ia,it,icf_sum(ir,1),im,1) + tr0_sum
          tr_sum(ip,ie,ia,it,icf_sum(ir,1),im,2) =
     &           tr_sum(ip,ie,ia,it,icf_sum(ir,1),im,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do

        do im = 1, nm
        do ir = 1, nr
        do it = 1, nt
        do ia = 1, na
        do ie = 1, ne
        do ip = 1, np
          tz0_sum = 0.0d0
          do iz = 1, nz+1
            tz0_sum = tz0_sum + tz0(ip,ie,ia,it,izf(ir,iz),im) / maxcas
          end do
          tz_sum(ip,ie,ia,it,izf_sum(ir,1),im,1) =
     &           tz_sum(ip,ie,ia,it,izf_sum(ir,1),im,1) + tz0_sum
          tz_sum(ip,ie,ia,it,icf_sum(ir,1),im,2) =
     &           tz_sum(ip,ie,ia,it,izf_sum(ir,1),im,2) + tz0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 9) then   ! time
        do im = 1, nm
        do iz = 1, nz
        do ir = 1, nr+1
        do ia = 1, na
        do ie = 1, ne
        do ip = 1, np
          tr0_sum = 0.0d0
          do it = 1, nt
            tr0_sum = tr0_sum + tr0(ip,ie,ia,it,icf(ir,iz),im) / maxcas
          end do
          tr_sum(ip,ie,ia,1,icf_sum(ir,iz),im,1) =
     &           tr_sum(ip,ie,ia,1,icf_sum(ir,iz),im,1) + tr0_sum
          tr_sum(ip,ie,ia,1,icf_sum(ir,iz),im,2) =
     &           tr_sum(ip,ie,ia,1,icf_sum(ir,iz),im,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do

        do im = 1, nm
        do iz = 1, nz+1
        do ir = 1, nr
        do ia = 1, na
        do ie = 1, ne
        do ip = 1, np
          tz0_sum = 0.0d0
          do it = 1, nt
            tz0_sum = tz0_sum + tz0(ip,ie,ia,it,izf(ir,iz),im) / maxcas
          end do
          tz_sum(ip,ie,ia,1,izf_sum(ir,iz),im,1) =
     &           tz_sum(ip,ie,ia,1,izf_sum(ir,iz),im,1) + tz0_sum
          tz_sum(ip,ie,ia,1,izf_sum(ir,iz),im,2) =
     &           tz_sum(ip,ie,ia,1,izf_sum(ir,iz),im,2) + tz0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 8 .or. itaxs_in == 10) then   ! ang
        do im = 1, nm
        do iz = 1, nz
        do ir = 1, nr+1
        do it = 1, nt
        do ie = 1, ne
        do ip = 1, np
          tr0_sum = 0.0d0
          do ia = 1, na
            tr0_sum = tr0_sum + tr0(ip,ie,ia,it,icf(ir,iz),im) / maxcas
          end do
          tr_sum(ip,ie,1,it,icf_sum(ir,iz),im,1) =
     &           tr_sum(ip,ie,1,it,icf_sum(ir,iz),im,1) + tr0_sum
          tr_sum(ip,ie,1,it,icf_sum(ir,iz),im,2) =
     &           tr_sum(ip,ie,1,it,icf_sum(ir,iz),im,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do

        do im = 1, nm
        do iz = 1, nz+1
        do ir = 1, nr
        do it = 1, nt
        do ie = 1, ne
        do ip = 1, np
          tz0_sum = 0.0d0
          do ia = 1, na
            tz0_sum = tz0_sum + tz0(ip,ie,ia,it,izf(ir,iz),im) / maxcas
          end do
          tz_sum(ip,ie,1,it,izf_sum(ir,iz),im,1) =
     &           tz_sum(ip,ie,1,it,izf_sum(ir,iz),im,1) + tz0_sum
          tz_sum(ip,ie,1,it,izf_sum(ir,iz),im,2) =
     &          tz_sum(ip,ie,1,it,izf_sum(ir,iz),im,2) + tz0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do

       end if
       return
       end

************************************************************************
*                                                                      *
      subroutine tsufxyz_sumover(m,maxcas,
     &                   np, ne, na, nt, nx, ny, nz, nm,  tr0x)
*                                                                      *
*     tr_sum(np,ne,na,nt,nx*ny*(nz+1),nm,2)                                   *
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0

      implicit double precision (a-h,o-z)
      dimension   tr0x(np,ne,na,nt,nx*ny*(nz+1),nm)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)


C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)

      do iax=1,itaxn(m)

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if
      if(n_tr_sum > 0) then

        call tsufxyz_sumover_sub(maxcas,itaxs(m,iax),
     &     np, ne, na,  nt, nx, ny, nz, nm,    tr0x,
     &     itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &     ittnm_sum(m,iax),
     &     itxnm_sum(m,iax),itynm_sum(m,iax),
     &     itznm_sum(m,iax),itmst_sum(m,iax),
     &     tr_sum)

        endif

      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine tsufxyz_sumover_sub(maxcas,itaxs_in,
     &           np, ne, na,  nt, nx, ny, nz, nm,    tr0,
     &           np_sum, ne_sum, na_sum, nt_sum,
     &           nx_sum, ny_sum, nz_sum, nm_sum, 
     &           tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      dimension   tr0(np,ne,na,nt,nx*ny*(nz+1),nm)
      dimension   tr_sum(np_sum,ne_sum,na_sum,nt_sum,
     &                   nx_sum*ny_sum*(nz_sum+1),nm_sum,2)

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny
      icf_sum(ix,iy,iz) = ix + ( iy - 1 ) * nx_sum 
     &                  + ( iz - 1 ) * nx_sum * ny_sum

      if(itaxs_in == 1 .or. itaxs_in == 14 ) then  ! energ
        do im = 1, nm
        do iz = 1, nz+1
        do iy = 1, ny
        do ix = 1, nx
        do it = 1, nt
        do ia = 1, na
        do ip = 1, np
          tr0_sum = 0.0d0
          do ie = 1, ne
            tr0_sum = tr0_sum
     &               + tr0(ip,ie,ia,it,icf(ix,iy,iz),im) / maxcas
          end do
          tr_sum(ip,1,ia,it,icf_sum(ix,iy,iz),im,1) =
     &         tr_sum(ip,1,ia,it,icf_sum(ix,iy,iz),im,1) + tr0_sum
          tr_sum(ip,1,ia,it,icf_sum(ix,iy,iz),im,2) =
     &        tr_sum(ip,1,ia,it,icf_sum(ix,iy,iz),im,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do
        end do
!
      elseif(itaxs_in == 3) then    ! x
        do im = 1, nm
        do iz = 1, nz+1
        do iy = 1, ny
        do it = 1, nt
        do ia = 1, na
        do ie = 1, ne
        do ip = 1, np
          tr0_sum = 0.0d0
          do ix = 1, nx
            tr0_sum = tr0_sum
     &              + tr0(ip,ie,ia,it,icf(ix,iy,iz),im) / maxcas
          end do
          tr_sum(ip,ie,ia,it,icf_sum(1,iy,iz),im,1) =
     &         tr_sum(ip,ie,ia,it,icf_sum(1,iy,iz),im,1) + tr0_sum
          tr_sum(ip,ie,ia,it,icf_sum(1,iy,iz),im,2) =
     &         tr_sum(ip,ie,ia,it,icf_sum(1,iy,iz),im,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 4) then    ! y
        do im = 1, nm
        do iz = 1, nz+1
        do ix = 1, nx
        do it = 1, nt
        do ia = 1, na
        do ie = 1, ne
        do ip = 1, np
          tr0_sum = 0.0d0
          do iy = 1, ny
            tr0_sum = tr0_sum
     &              + tr0(ip,ie,ia,it,icf(ix,iy,iz),im) / maxcas
          end do
          tr_sum(ip,ie,ia,it,icf_sum(ix,1,iz),im,1) =
     &         tr_sum(ip,ie,ia,it,icf_sum(ix,1,iz),im,1) + tr0_sum
          tr_sum(ip,ie,ia,it,icf_sum(ix,1,iz),im,2) =
     &         tr_sum(ip,ie,ia,it,icf_sum(ix,1,iz),im,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 5) then    ! z
        do im = 1, nm
        do iy = 1, ny
        do ix = 1, nx
        do it = 1, nt
        do ia = 1, na
        do ie = 1, ne
        do ip = 1, np
          tr0_sum = 0.0d0
          do iz = 1, nz+1
            tr0_sum = tr0_sum
     &              + tr0(ip,ie,ia,it,icf(ix,iy,iz),im) / maxcas
          end do
          tr_sum(ip,ie,ia,it,icf_sum(ix,iy,1),im,1) =
     &         tr_sum(ip,ie,ia,it,icf_sum(ix,iy,1),im,1) + tr0_sum
          tr_sum(ip,ie,ia,it,icf_sum(ix,iy,1),im,2) =
     &         tr_sum(ip,ie,ia,it,icf_sum(ix,iy,1),im,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 9) then   ! time
        do im = 1, nm
        do iz = 1, nz+1
        do iy = 1, ny
        do ix = 1, nx
        do ia = 1, na
        do ie = 1, ne
        do ip = 1, np
          tr0_sum = 0.0d0
          do it = 1, nt
            tr0_sum = tr0_sum
     &              + tr0(ip,ie,ia,it,icf(ix,iy,iz),im) / maxcas
          end do
          tr_sum(ip,ie,ia,1,icf_sum(ix,iy,iz),im,1) =
     &           tr_sum(ip,ie,ia,1,icf_sum(ix,iy,iz),im,1) + tr0_sum
          tr_sum(ip,ie,ia,1,icf_sum(ix,iy,iz),im,2) =
     &          tr_sum(ip,ie,ia,1,icf_sum(ix,iy,iz),im,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 8 .or. itaxs_in == 10) then   ! ang
        do im = 1, nm
        do iz = 1, nz+1
        do iy = 1, ny
        do ix = 1, nx
        do it = 1, nt
        do ie = 1, ne
        do ip = 1, np
          tr0_sum = 0.0d0
          do ia = 1, na
            tr0_sum = tr0_sum
     &              + tr0(ip,ie,ia,it,icf(ix,iy,iz),im) / maxcas
          end do
          tr_sum(ip,ie,1,it,icf_sum(ix,iy,iz),im,1) =
     &           tr_sum(ip,ie,1,it,icf_sum(ix,iy,iz),im,1) + tr0_sum
          tr_sum(ip,ie,1,it,icf_sum(ix,iy,iz),im,2) =
     &          tr_sum(ip,ie,1,it,icf_sum(ix,iy,iz),im,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do
        end do

       end if
       return
       end


************************************************************************
*                                                                      *
      subroutine psufreg_sumover_stdev(mode,m,ip,ie,ia,it,ir,im,
     &                fact_in,ew,aw,tw,ar,ew_sum,aw_sum,tw_sum,ar_sum)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)
     
      do iax=1,itaxn(m)

        if(mode == 0) then

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if

        else

          n_tr_sum = italsize_sum(m,iax)
          tr_sum => trRES_sum(irestalm_sum(m,iax):)

        endif

        if(n_tr_sum > 0) then
             sum_fact = fact_in/ar/ew/aw/tw
          if((itaxs(m,iax) == 1 .or. itaxs(m,iax) == 14) .and.
     &        ie == 1) then
             sum_fact = sum_fact * ew / ew_sum
          else if(itaxs(m,iax) == 2 .and. ir == 1) then
             sum_fact = sum_fact * ar / ar_sum
          else if(itaxs(m,iax) == 9 .and. it == 1) then
             sum_fact = sum_fact * tw / tw_sum
          else if((itaxs(m,iax) == 8 .or. itaxs(m,iax) == 10) .and.
     &        ia == 1) then
             sum_fact = sum_fact * aw / aw_sum
          else
            sum_fact = 0.0d0
          endif
          if(sum_fact /= 0.0d0) then
            call psufreg_sumover_stdev_sub(mode,m,sum_fact,
     &         ip,ie,ia,it,ir,im,
     &         itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &         ittnm_sum(m,iax),
     &         itrcn_sum(m,iax),itmst_sum(m,iax),tr_sum)
          endif

        endif
      enddo

      return
      end


************************************************************************
*                                                                      *
      subroutine psufreg_sumover_stdev_ntf(mode,m,ntf,
     &                ip,ie,ia,it,ir,im,
     &                fact_in,ew,aw,tw,ar,ew_sum,aw_sum,tw_sum,ar_sum)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)
     
      do iax=1,itaxn(m)

        if(mode == 0) then

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM_NTF(tr_sum,m,ntf,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM_NTF(tr_sum,m,ntf,iax)
C for nonshared_tally option
!$       end if

        else

          ntfbase = mtalsize_sum(m) * (ntf-1)
          n_tr_sum = italsize_sum(m,iax)
          tr_sum => trRES_sum(irestalm_sum(m,iax)+ntfbase:)

        endif

        if(n_tr_sum > 0) then
             sum_fact = fact_in/ar/ew/aw/tw
          if((itaxs(m,iax) == 1 .or. itaxs(m,iax) == 14) .and.
     &        ie == 1) then
             sum_fact = sum_fact * ew / ew_sum
          else if(itaxs(m,iax) == 2 .and. ir == 1) then
             sum_fact = sum_fact * ar / ar_sum
          else if(itaxs(m,iax) == 9 .and. it == 1) then
             sum_fact = sum_fact * tw / tw_sum
          else if((itaxs(m,iax) == 8 .or. itaxs(m,iax) == 10) .and.
     &        ia == 1) then
             sum_fact = sum_fact * aw / aw_sum
          else
            sum_fact = 0.0d0
          endif
          if(sum_fact /= 0.0d0) then
            call psufreg_sumover_stdev_sub(mode,m,sum_fact,
     &         ip,ie,ia,it,ir,im,
     &         itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &         ittnm_sum(m,iax),
     &         itrcn_sum(m,iax),itmst_sum(m,iax),tr_sum)
          endif

        endif
      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine psufreg_sumover_stdev_sub(mode,m,sum_fact,
     &               ip,ie,ia,it,ir,im,
     &               np_sum,ne_sum,na_sum,nt_sum,nr_sum,nm_sum,tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      dimension   tr_sum(np_sum,ne_sum,na_sum,nt_sum,nr_sum,nm_sum,2)
      
      if(tr_sum(ip,ie,ia,it,ir,im,1) > 0.0) then
        if(mode == 0) then
          call calc_stdev(m,Xa,sigx,
     &                  tr_sum(ip,ie,ia,it,ir,im,1),
     &                  tr_sum(ip,ie,ia,it,ir,im,2),
     &                  sum_fact)
        else
          sum_fact_r = 1.0d0 / sum_fact
          call invert_stdev(m,Xa,sigx,
     &                  tr_sum(ip,ie,ia,it,ir,im,1),
     &                  tr_sum(ip,ie,ia,it,ir,im,2),
     &                  sum_fact_r)
        endif

        tr_sum(ip,ie,ia,it,ir,im,1) = Xa
        tr_sum(ip,ie,ia,it,ir,im,2) = sigx

      else

        tr_sum(ip,ie,ia,it,ir,im,2) = 0.0

      end if

      return
      end

************************************************************************
*                                                                      *
      subroutine  psufreg_sumover_getput(mode,m,iax,
     &    npstepi,nestepi,nastepi,ntstepi,nrstepi,nmstepi,
     &    ipi,iei,iai,iti,iri,imi,maxtott,tott_sum)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      real(8) :: tott_sum(maxtott,2)

      real(8),allocatable :: tott_in(:,:)
      character :: chin*200

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

      real(8),pointer :: tr_sum(:)

      if(mode == 0) then
 
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if

         call psufreg_sumover_getput_sub(mode,
     &    npstepi,nestepi,nastepi,ntstepi,nrstepi,nmstepi,
     &    ipi,iei,iai,iti,iri,imi,maxtott,tott_sum,
     &    itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &    ittnm_sum(m,iax),itrcn_sum(m,iax),itmst_sum(m,iax),tr_sum)


      else

        tr_sum =>  trRES_sum(irestalm_sum(m,iax):)

        nsame = npstepi * nestepi * nastepi * ntstepi * nrstepi
     &        * nmstepi
        allocate(tott_in(nsame,2))

        read(mode,'(a)') chin
        read(mode,'(26x,1000(1pe13.4,0pf8.4))')
     &  (tott_in(i,1),tott_in(i,2),i=1,nsame)

        call psufreg_sumover_getput_sub(mode,
     &    npstepi,nestepi,nastepi,ntstepi,nrstepi,nmstepi,
     &    ipi,iei,iai,iti,iri,imi,nsame,tott_in,
     &    itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &    ittnm_sum(m,iax),itrcn_sum(m,iax),itmst_sum(m,iax),tr_sum)

        deallocate(tott_in)

      endif

      return
      end

************************************************************************
*                                                                      *
      subroutine psufreg_sumover_getput_sub(mode,
     &    npstepi,nestepi,nastepi,ntstepi,nrstepi,nmstepi,
     &    ipi,iei,iai,iti,iri,imi,maxtott,tott_sum,
     &    np_sum,ne_sum,na_sum,nt_sum,nr_sum,nm_sum,tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      real(8) :: tr_sum(np_sum,ne_sum,na_sum,nt_sum,nr_sum,nm_sum,2)

      real(8) :: tott_sum(maxtott,2)

        i = 0
        do imloop=1,nmstepi
          do irloop=1,nrstepi
            do itloop=1,ntstepi
             do ialoop=1,nastepi
              do ieloop=1,nestepi
                do iploop=1,npstepi
                  i = i + 1
                  do k=1,2
                    if(mode == 0) then
                      tott_sum(i,k) =
     &                tr_sum(ipi+iploop-1,iei+ieloop-1,iai+ialoop-1,
     &                iti+itloop-1,iri+irloop-1,
     &                imi+imloop-1,k)
                    else
                      tr_sum(ipi+iploop-1,iei+ieloop-1,iai+ialoop-1,
     &                       iti+itloop-1,iri+irloop-1,
     &                       imi+imloop-1,k)
     &                = tott_sum(i,k)
                    endif
                  enddo
                end do
              end do
             end do
            end do
          end do
        end do

      return
      end


************************************************************************
*                                                                      *
      subroutine psufrz_sumover_stdev(mode,m,ip,ie,ia,it,ir,iz,im,
     &                fact_in,ew,aw,tw,vl_in,ew_sum,aw_sum,tw_sum,
     &                vl_r_in,vl_z_in)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)
     
      vl   = vl_in
      vl_r = vl_r_in
      vl_z = vl_z_in
      if(vl == 0.0d0) then
        vl = 1.0d0
      endif
      if(vl_r == 0.0d0) then
        vl_r = 1.0d0
      endif
      if(vl_z == 0.0d0) then
        vl_z = 1.0d0
      endif

      do iax=1,itaxn(m)

        if(mode == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_2_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            n_tr_sum = italsize_2_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if
       else

          n_tr_sum = italsize_2_sum(m,iax)
          tr_sum => trRES_sum(irestalm_sum(m,iax):)

       endif

        if(n_tr_sum > 0) then
           sum_fact = fact_in/vl/aw/ew/tw
          if((itaxs(m,iax) == 1 .or. itaxs(m,iax) == 14 ) .and.
     &                               ie == 1) then
             sum_fact = sum_fact * ew / ew_sum
          else if(itaxs(m,iax) == 5 .and. iz == 1) then
             sum_fact = sum_fact * vl / vl_z
          else if(itaxs(m,iax) == 6 .and. ir == 1) then
             sum_fact = sum_fact * vl / vl_r
          else if(itaxs(m,iax) == 9 .and. it == 1) then
             sum_fact = sum_fact * tw / tw_sum
          else if((itaxs(m,iax) == 8 .or. itaxs(m,iax) == 10) .and.
     &             ia == 1) then
             sum_fact = sum_fact * aw / aw_sum
          else
            sum_fact = 0.0d0
          endif
          if(sum_fact /= 0.0d0) then
            call psufrz_sumover_stdev_sub(mode,m,sum_fact,
     &         ip,ie,ia,it,ir,iz,im,
     &         itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &         ittnm_sum(m,iax),itrnm_sum(m,iax),itznm_sum(m,iax),
     &         itmst_sum(m,iax),tr_sum)
          endif

        endif
      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine psufrz_sumover_stdev_ntf(mode,m,ntf,
     &                ip,ie,ia,it,ir,iz,im,
     &                fact_in,ew,aw,tw,vl_in,ew_sum,aw_sum,tw_sum,
     &                vl_r_in,vl_z_in)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)
 
      vl   = vl_in
      vl_r = vl_r_in
      vl_z = vl_z_in
      if(vl == 0.0d0) then
        vl = 1.0d0
      endif
      if(vl_r == 0.0d0) then
        vl_r = 1.0d0
      endif
      if(vl_z == 0.0d0) then
        vl_z = 1.0d0
      endif

      do iax=1,itaxn(m)

        if(mode == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_2_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM_NTF(tr_sum,m,ntf,iax)
!$       else
            n_tr_sum = italsize_2_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM_NTF(tr_sum,m,ntf,iax)
C for nonshared_tally option
!$       end if
       else

          ntfbase =  mtalsize_sum(m) * (ntf-1)
          n_tr_sum = italsize_2_sum(m,iax)
          tr_sum => trRES_sum(irestalm_sum(m,iax)+ntfbase:)

       endif

        if(n_tr_sum > 0) then
           sum_fact = fact_in/vl/aw/ew/tw
          if((itaxs(m,iax) == 1 .or. itaxs(m,iax) == 14 ) .and.
     &                               ie == 1) then
             sum_fact = sum_fact * ew / ew_sum
          else if(itaxs(m,iax) == 5 .and. iz == 1) then
             sum_fact = sum_fact * vl / vl_z
          else if(itaxs(m,iax) == 6 .and. ir == 1) then
             sum_fact = sum_fact * vl / vl_r
          else if(itaxs(m,iax) == 9 .and. it == 1) then
             sum_fact = sum_fact * tw / tw_sum
          else if((itaxs(m,iax) == 8 .or. itaxs(m,iax) == 10) .and.
     &             ia == 1) then
             sum_fact = sum_fact * aw / aw_sum
          else
            sum_fact = 0.0d0
          endif
          if(sum_fact /= 0.0d0) then
            call psufrz_sumover_stdev_sub(mode,m,sum_fact,
     &         ip,ie,ia,it,ir,iz,im,
     &         itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &         ittnm_sum(m,iax),itrnm_sum(m,iax),itznm_sum(m,iax),
     &         itmst_sum(m,iax),tr_sum)
          endif

        endif
      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine psufrz_sumover_stdev_sub(mode,m,sum_fact,
     &                   ip,ie,ia,it,ir,iz,im,
     &                   np_sum,ne_sum,na_sum,nt_sum,nr_sum,nz_sum,
     &                   nm_sum,tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      dimension   tr_sum(np_sum,ne_sum,na_sum,nt_sum,(nr_sum+1)*nz_sum,
     &                   nm_sum,2)
         icf_sum(ir,iz) = ir + ( iz - 1 ) * (nr_sum + 1)

      
      if(tr_sum(ip,ie,ia,it,icf_sum(ir,iz),im,1) > 0.0) then
        if(mode == 0) then
          call calc_stdev(m,Xa,sigx,
     &                  tr_sum(ip,ie,ia,it,icf_sum(ir,iz),im,1),
     &                  tr_sum(ip,ie,ia,it,icf_sum(ir,iz),im,2),
     &                  sum_fact)
        else
          sum_fact_r = 1.0d0 / sum_fact
          call invert_stdev(m,Xa,sigx,
     &                  tr_sum(ip,ie,ia,it,icf_sum(ir,iz),im,1),
     &                  tr_sum(ip,ie,ia,it,icf_sum(ir,iz),im,2),
     &                  sum_fact_r)
        endif

        tr_sum(ip,ie,ia,it,icf_sum(ir,iz),im,1) = Xa
        tr_sum(ip,ie,ia,it,icf_sum(ir,iz),im,2) = sigx

      else
        tr_sum(ip,ie,ia,it,icf_sum(ir,iz),im,2) = 0.0

      end if

      return
      end

************************************************************************
*                                                                      *
      subroutine psufrz_sumover_tz_stdev(mode,m,ip,ie,ia,it,ir,iz,im,
     &                fact_in,ew,aw,tw,ar,ew_sum,aw_sum,tw_sum,
     &                ar_sum)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tz_sum(:)
     
      do iax=1,itaxn(m)

        if(mode == 0) then

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TZ_HEAD_POINTER0_SUM(tz_sum,m,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TZ_HEAD_POINTER_SUM(tz_sum,m,iax)
C for nonshared_tally option
!$       end if

       else

          n_tr_sum = italsize_sum(m,iax) - italsize_2_sum(m,iax)
          tz_sum =>
     &      trRES_sum(irestalm_sum(m,iax)+italsize_2_sum(m,iax):)

       endif

        if(n_tr_sum > 0) then
             sum_fact = fact_in/ar/aw/ew/tw
          if((itaxs(m,iax) == 1 .or. itaxs(m,iax) == 14 ) .and.
     &                               ie == 1) then
             sum_fact = sum_fact * ew / ew_sum
          else if(itaxs(m,iax) == 5 .and. iz == 1) then

          else if(itaxs(m,iax) == 6 .and. ir == 1) then
             sum_fact = sum_fact * ar / ar_sum
          else if(itaxs(m,iax) == 9 .and. it == 1) then
             sum_fact = sum_fact * tw / tw_sum
          else if((itaxs(m,iax) == 8 .or. itaxs(m,iax) == 10) .and.
     &             ia == 1) then
             sum_fact = sum_fact * aw / aw_sum
          else
            sum_fact = 0.0d0
          endif
          if(sum_fact /= 0.0d0) then
            call psufrz_sumover_sub_tz_stdev(mode,m,sum_fact,
     &         ip,ie,ia,it,ir,iz,im,
     &         itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &         ittnm_sum(m,iax),itrnm_sum(m,iax),itznm_sum(m,iax),
     &         itmst_sum(m,iax),tz_sum)
          endif

        endif
      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine psufrz_sumover_tz_stdev_ntf(mode,m,ntf,
     &                ip,ie,ia,it,ir,iz,im,
     &                fact_in,ew,aw,tw,ar,ew_sum,aw_sum,tw_sum,
     &                ar_sum)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tz_sum(:)
     
      do iax=1,itaxn(m)

        if(mode == 0) then

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TZ_HEAD_POINTER0_SUM_NTF(tz_sum,m,ntf,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TZ_HEAD_POINTER_SUM_NTF(tz_sum,m,ntf,iax)
C for nonshared_tally option
!$       end if

       else

          ntfbase = mtalsize_sum(m) * (ntf-1)
          n_tr_sum = italsize_sum(m,iax) - italsize_2_sum(m,iax)
          tz_sum =>
     &      trRES_sum(irestalm_sum(m,iax)+italsize_2_sum(m,iax)
     &               +ntfbase:)

       endif

        if(n_tr_sum > 0) then
             sum_fact = fact_in/ar/aw/ew/tw
          if((itaxs(m,iax) == 1 .or. itaxs(m,iax) == 14 ) .and.
     &                               ie == 1) then
             sum_fact = sum_fact * ew / ew_sum
          else if(itaxs(m,iax) == 5 .and. iz == 1) then

          else if(itaxs(m,iax) == 6 .and. ir == 1) then
             sum_fact = sum_fact * ar / ar_sum
          else if(itaxs(m,iax) == 9 .and. it == 1) then
             sum_fact = sum_fact * tw / tw_sum
          else if((itaxs(m,iax) == 8 .or. itaxs(m,iax) == 10) .and.
     &             ia == 1) then
             sum_fact = sum_fact * aw / aw_sum
          else
            sum_fact = 0.0d0
          endif
          if(sum_fact /= 0.0d0) then
            call psufrz_sumover_sub_tz_stdev(mode,m,sum_fact,
     &         ip,ie,ia,it,ir,iz,im,
     &         itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &         ittnm_sum(m,iax),itrnm_sum(m,iax),itznm_sum(m,iax),
     &         itmst_sum(m,iax),tz_sum)
          endif

        endif
      enddo

      return
      end


************************************************************************
*                                                                      *
      subroutine psufrz_sumover_sub_tz_stdev(mode,m,sum_fact,
     &                   ip,ie,ia,it,ir,iz,im,
     &                   np_sum,ne_sum,na_sum,nt_sum,nr_sum,nz_sum,
     &                   nm_sum,tz_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      dimension   tz_sum(np_sum,ne_sum,na_sum,nt_sum,nr_sum*(nz_sum+1),
     &                   nm_sum,2)
         izf_sum(ir,iz) = iz + ( ir - 1 ) * (nz_sum + 1)

      
      if(tz_sum(ip,ie,ia,it,izf_sum(ir,iz),im,1) > 0.0) then
        if(mode == 0) then
          call calc_stdev(m,Xa,sigx,
     &                  tz_sum(ip,ie,ia,it,izf_sum(ir,iz),im,1),
     &                  tz_sum(ip,ie,ia,it,izf_sum(ir,iz),im,2),
     &                  sum_fact)
        else
          sum_fact_r = 1.0d0 / sum_fact
          call invert_stdev(m,Xa,sigx,
     &                  tz_sum(ip,ie,ia,it,izf_sum(ir,iz),im,1),
     &                  tz_sum(ip,ie,ia,it,izf_sum(ir,iz),im,2),
     &                  sum_fact_r)
        endif

        tz_sum(ip,ie,ia,it,izf_sum(ir,iz),im,1) = Xa
        tz_sum(ip,ie,ia,it,izf_sum(ir,iz),im,2) = sigx

      else

        tz_sum(ip,ie,ia,it,izf_sum(ir,iz),im,2) = 0.0

      end if

      return
      end

************************************************************************
*                                                                      *
      subroutine  psufrz_sumover_getput(mode,m,iax,
     &    npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &    ipi,iei,iai,iti,iri,izi,imi,maxtott,tott_sum)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      real(8) :: tott_sum(maxtott,2)

      real(8),allocatable :: tott_in(:,:)
      character :: chin*200


C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

      real(8),pointer :: tr_sum(:)

      if(mode == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if
         call psufrz_sumover_getput_sub(mode,
     &    npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &    ipi,iei,iai,iti,iri,izi,imi,maxtott,tott_sum,
     &    itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &    ittnm_sum(m,iax),itrnm_sum(m,iax),itznm_sum(m,iax),
     &    itmst_sum(m,iax),tr_sum)

      else

        tr_sum => trRES_sum(irestalm_sum(m,iax):)

        nsame = npstepi * nestepi * nastepi * ntstepi
     &        * (nrstepi + 1) * nzstepi * nmstepi
        allocate(tott_in(nsame,2))

        read(mode,'(a)') chin
        read(mode,'(26x,1000(1pe13.4,0pf8.4))')
     &  (tott_in(i,1),tott_in(i,2),i=1,nsame)

        call psufrz_sumover_getput_sub(mode,
     &    npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &    ipi,iei,iai,iti,iri,izi,imi,nsame,tott_in,
     &    itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &    ittnm_sum(m,iax),itrnm_sum(m,iax),itznm_sum(m,iax),
     &    itmst_sum(m,iax),tr_sum)

       deallocate(tott_in)

      endif

      return
      end

************************************************************************
*                                                                      *
      subroutine psufrz_sumover_getput_sub(mode,
     &    npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &    ipi,iei,iai,iti,iri,izi,imi,maxtott,tott_sum,
     &    np_sum,ne_sum,na_sum,nt_sum,nr_sum,nz_sum,nm_sum,tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      real(8) :: tr_sum(np_sum,ne_sum,na_sum,nt_sum,(nr_sum+1)*nz_sum,
     &                  nm_sum,2)

      real(8) :: tott_sum(maxtott,2)

      icf_sum(ir,iz) = ir + ( iz - 1 ) * (nr_sum +1)

        i = 0
        do imloop=1,nmstepi
          do izloop=1,nzstepi
           do irloop=1,nrstepi
            do itloop=1,ntstepi
             do ialoop=1,nastepi
              do ieloop=1,nestepi
                do iploop=1,npstepi
                  i = i + 1
                  do k=1,2
                    if(mode == 0) then
                      tott_sum(i,k) =
     &                tr_sum(ipi+iploop-1,iei+ieloop-1,
     &                iai+ialoop-1,iti+itloop-1,
     &                icf_sum(iri+irloop-1,izi+izloop-1),
     &                imi+imloop-1,k)
                    else
                      tr_sum(ipi+iploop-1,iei+ieloop-1,
     &                iai+ialoop-1,iti+itloop-1,
     &                icf_sum(iri+irloop-1,izi+izloop-1),
     &                imi+imloop-1,k)
     &                = tott_sum(i,k)
                    endif
                  enddo
                end do
              end do
            end do
          end do
        end do
       end do
      end do

      return
      end

************************************************************************
      subroutine  psufrz_sumover_getput_z(mode,m,iax,
     &    npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &    ipi,iei,iai,iti,iri,izi,imi,maxtott,tott_sum)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      real(8) :: tott_sum(maxtott,2)

      real(8),allocatable :: tott_in(:,:)
      character :: chin*200

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

      real(8),pointer :: tz_sum(:)

      if(mode == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TZ_HEAD_POINTER0_SUM(tz_sum,m,iax)
!$       else
            call GET_TZ_HEAD_POINTER_SUM(tz_sum,m,iax)
C for nonshared_tally option
!$       end if

         call psufrz_sumover_getput_sub_z(mode,
     &     npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &     ipi,iei,iai,iti,iri,izi,imi,maxtott,tott_sum,
     &     itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &     ittnm_sum(m,iax),itrnm_sum(m,iax),itznm_sum(m,iax),
     &     itmst_sum(m,iax),tz_sum)
      else

        tz_sum =>
     &     trRES_sum(irestalm_sum(m,iax)+italsize_2_sum(m,iax):)

        nsame = npstepi * nestepi * nastepi * ntstepi
     &        * nrstepi * (nzstepi + 1) * nmstepi
        allocate(tott_in(nsame,2))

        read(mode,'(a)') chin
        read(mode,'(26x,1000(1pe13.4,0pf8.4))')
     &  (tott_in(i,1),tott_in(i,2),i=1,nsame)

         call psufrz_sumover_getput_sub_z(mode,
     &     npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &     ipi,iei,iai,iti,iri,izi,imi,nsame,tott_in,
     &     itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &     ittnm_sum(m,iax),itrnm_sum(m,iax),itznm_sum(m,iax),
     &     itmst_sum(m,iax),tz_sum)

        deallocate(tott_in)

      endif

      return
      end

************************************************************************
*                                                                      *
      subroutine psufrz_sumover_getput_sub_z(mode,
     &    npstepi,nrstepi,nzstepi,nestepi,nastepi,nmstepi,ntstepi,
     &    ipi,iei,iai,iti,iri,izi,imi,maxtott,tott_sum,
     &    np_sum,ne_sum,na_sum,nt_sum,nr_sum,nz_sum,nm_sum,tz_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      real(8) :: tz_sum(np_sum,ne_sum,na_sum,nt_sum,nr_sum*(nz_sum+1),
     &                  nm_sum,2)

      real(8) :: tott_sum(maxtott,2)

      izf_sum(ir,iz) = iz + ( ir - 1 ) * (nz_sum +1)

        i = 0
        do imloop=1,nmstepi
          do izloop=1,nzstepi
           do irloop=1,nrstepi
            do itloop=1,ntstepi
             do ialoop=1,nastepi
              do ieloop=1,nestepi
                do iploop=1,npstepi
                  i = i + 1
                  do k=1,2
                    if(mode == 0) then
                      tott_sum(i,k) =
     &                tz_sum(ipi+iploop-1,iei+ieloop-1,
     &                iai+ialoop-1,iti+itloop-1,
     &                izf_sum(iri+irloop-1,izi+izloop-1),
     &                imi+imloop-1,k)
                    else
                      tz_sum(ipi+iploop-1,iei+ieloop-1,
     &                iai+ialoop-1,iti+itloop-1,
     &                izf_sum(iri+irloop-1,izi+izloop-1),
     &                imi+imloop-1,k)
     &                = tott_sum(i,k)
                    endif
                  enddo
                end do
              end do
            end do
          end do
        end do
       end do
      end do

      return
      end

************************************************************************
*                                                                      *
      subroutine psufxyz_sumover_stdev(mode,m,ip,ie,ia,it,ix,iy,iz,im,
     &            fact_in,ew,aw,tw,vl,
     &            ew_sum,aw_sum,tw_sum,vl_x,vl_y)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)
     
      do iax=1,itaxn(m)

        if(mode == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if
        else

          n_tr_sum = italsize_sum(m,iax)
          tr_sum => trRES_sum(irestalm_sum(m,iax):)

        endif

        if(n_tr_sum > 0) then
           sum_fact = fact_in/vl/ew/tw

          if((itaxs(m,iax) == 1 .or. itaxs(m,iax) == 14).and.
     &                                                ie == 1) then
             sum_fact = sum_fact * ew / ew_sum
          else if(itaxs(m,iax) == 3 .and. ix == 1) then
             sum_fact = sum_fact * vl / vl_x
          else if(itaxs(m,iax) == 4 .and. iy == 1) then
             sum_fact = sum_fact * vl / vl_y
          else if(itaxs(m,iax) == 5 .and. iz == 1) then

          else if(itaxs(m,iax) == 9 .and. it == 1) then
             sum_fact = sum_fact * tw / tw_sum
          else if((itaxs(m,iax) == 8 .or. itaxs(m,iax) == 10).and.
     &                                               ia == 1) then
             sum_fact = sum_fact * aw / aw_sum
          else
            sum_fact = 0.0d0
          endif
          if(sum_fact /= 0.0d0) then
            call psufxyz_sumover_stdev_sub(mode,m,sum_fact,
     &         ip,ie,ia,it,ix,iy,iz,im,
     &         itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &         ittnm_sum(m,iax),
     &         itxnm_sum(m,iax),itynm_sum(m,iax),itznm_sum(m,iax),
     &         itmst_sum(m,iax),tr_sum)
          endif

        endif
      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine psufxyz_sumover_stdev_ntf(mode,m,ntf,
     &            ip,ie,ia,it,ix,iy,iz,im,
     &            fact_in,ew,aw,tw,vl,
     &            ew_sum,aw_sum,tw_sum,vl_x,vl_y)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)
     
      do iax=1,itaxn(m)

        if(mode == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM_NTF(tr_sum,m,ntf,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM_NTF(tr_sum,m,ntf,iax)
C for nonshared_tally option
!$       end if
        else

          ntfbase = mtalsize_sum(m) * (ntf-1)
          n_tr_sum = italsize_sum(m,iax)
          tr_sum => trRES_sum(irestalm_sum(m,iax)+ntfbase:)

        endif

        if(n_tr_sum > 0) then
           sum_fact = fact_in/vl/ew/tw

          if((itaxs(m,iax) == 1 .or. itaxs(m,iax) == 14).and.
     &                                                ie == 1) then
             sum_fact = sum_fact * ew / ew_sum
          else if(itaxs(m,iax) == 3 .and. ix == 1) then
             sum_fact = sum_fact * vl / vl_x
          else if(itaxs(m,iax) == 4 .and. iy == 1) then
             sum_fact = sum_fact * vl / vl_y
          else if(itaxs(m,iax) == 5 .and. iz == 1) then

          else if(itaxs(m,iax) == 9 .and. it == 1) then
             sum_fact = sum_fact * tw / tw_sum
          else if((itaxs(m,iax) == 8 .or. itaxs(m,iax) == 10).and.
     &                                               ia == 1) then
             sum_fact = sum_fact * aw / aw_sum
          else
            sum_fact = 0.0d0
          endif
          if(sum_fact /= 0.0d0) then
            call psufxyz_sumover_stdev_sub(mode,m,sum_fact,
     &         ip,ie,ia,it,ix,iy,iz,im,
     &         itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &         ittnm_sum(m,iax),
     &         itxnm_sum(m,iax),itynm_sum(m,iax),itznm_sum(m,iax),
     &         itmst_sum(m,iax),tr_sum)
          endif

        endif
      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine psufxyz_sumover_rpp_stdev(mode,m,
     &            ip,ie,ia,it,ix,iy,iz,im,
     &            fact_in,ew,aw,tw,vl,
     &            ew_sum,aw_sum,tw_sum,vl_x,vl_y,vl_z)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)
     
      do iax=1,itaxn(m)

        if(mode == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if

        else

          n_tr_sum = italsize_sum(m,iax)
          tr_sum => trRES_sum(irestalm_sum(m,iax):)

        endif

        if(n_tr_sum > 0) then
          sum_fact = fact_in/vl/ew/aw/tw
          if((itaxs(m,iax) == 1 .or. itaxs(m,iax) == 14).and.
     &                                                ie == 1) then
             sum_fact = sum_fact * ew / ew_sum
          else if(itaxs(m,iax) == 3 .and. ix == 1) then
             sum_fact = sum_fact * vl / vl_x
          else if(itaxs(m,iax) == 4 .and. iy == 1) then
             sum_fact = sum_fact * vl / vl_y
          else if(itaxs(m,iax) == 5 .and. iz == 1) then
             sum_fact = sum_fact * vl / vl_z
          else if(itaxs(m,iax) == 9 .and. it == 1) then
             sum_fact = sum_fact * tw / tw_sum
          else if((itaxs(m,iax) == 8 .or. itaxs(m,iax) == 10).and.
     &                                               ia == 1) then
             sum_fact = sum_fact * aw /aw_sum
          else
            sum_fact = 0.0d0
          endif
          if(sum_fact /= 0.0d0) then
            call psufxyz_sumover_stdev_sub(mode,m,sum_fact,
     &         ip,ie,ia,it,ix,iy,iz,im,
     &         itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &         ittnm_sum(m,iax),
     &         itxnm_sum(m,iax),itynm_sum(m,iax),itznm_sum(m,iax),
     &         itmst_sum(m,iax),tr_sum)
          endif

        endif
      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine psufxyz_sumover_rpp_stdev_ntf(mode,m,ntf,
     &            ip,ie,ia,it,ix,iy,iz,im,
     &            fact_in,ew,aw,tw,vl,
     &            ew_sum,aw_sum,tw_sum,vl_x,vl_y,vl_z)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)
     
      do iax=1,itaxn(m)

        if(mode == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM_NTF(tr_sum,m,ntf,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM_NTF(tr_sum,m,ntf,iax)
C for nonshared_tally option
!$       end if

        else

          n_tr_sum = italsize_sum(m,iax)
          ntfbase = mtalsize_sum(m) * (ntf-1)
          tr_sum => trRES_sum(irestalm_sum(m,iax)+ntfbase:)

        endif

        if(n_tr_sum > 0) then
          sum_fact = fact_in/vl/ew/aw/tw
          if((itaxs(m,iax) == 1 .or. itaxs(m,iax) == 14).and.
     &                                                ie == 1) then
             sum_fact = sum_fact * ew / ew_sum
          else if(itaxs(m,iax) == 3 .and. ix == 1) then
             sum_fact = sum_fact * vl / vl_x
          else if(itaxs(m,iax) == 4 .and. iy == 1) then
             sum_fact = sum_fact * vl / vl_y
          else if(itaxs(m,iax) == 5 .and. iz == 1) then
             sum_fact = sum_fact * vl / vl_z
          else if(itaxs(m,iax) == 9 .and. it == 1) then
             sum_fact = sum_fact * tw / tw_sum
          else if((itaxs(m,iax) == 8 .or. itaxs(m,iax) == 10).and.
     &                                               ia == 1) then
             sum_fact = sum_fact * aw /aw_sum
          else
            sum_fact = 0.0d0
          endif
          if(sum_fact /= 0.0d0) then
            call psufxyz_sumover_stdev_sub(mode,m,sum_fact,
     &         ip,ie,ia,it,ix,iy,iz,im,
     &         itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &         ittnm_sum(m,iax),
     &         itxnm_sum(m,iax),itynm_sum(m,iax),itznm_sum(m,iax),
     &         itmst_sum(m,iax),tr_sum)
          endif

        endif
      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine psufxyz_sumover_stdev_sub(mode,m,sum_fact,
     &                   ip,ie,ia,it,ix,iy,iz,im,
     &                   np_sum,ne_sum,na_sum,nt_sum,nx_sum,ny_sum,
     &                   nz_sum,nm_sum,tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      dimension   tr_sum(np_sum,ne_sum,na_sum,nt_sum,
     &                   nx_sum*ny_sum*(nz_sum+1),nm_sum,2)
      
      icf_sum(ix,iy,iz) = ix + ( iy - 1 ) * nx_sum
     &                  + ( iz - 1 ) * nx_sum * ny_sum

      if(tr_sum(ip,ie,ia,it,icf_sum(ix,iy,iz),im,1) > 0.0) then
        if(mode == 0) then
          call calc_stdev(m,Xa,sigx,
     &                  tr_sum(ip,ie,ia,it,icf_sum(ix,iy,iz),im,1),
     &                  tr_sum(ip,ie,ia,it,icf_sum(ix,iy,iz),im,2),
     &                  sum_fact)

        else
          sum_fact_r = 1.0d0 / sum_fact
          call invert_stdev(m,Xa,sigx,
     &                  tr_sum(ip,ie,ia,it,icf_sum(ix,iy,iz),im,1),
     &                  tr_sum(ip,ie,ia,it,icf_sum(ix,iy,iz),im,2),
     &                  sum_fact_r)
        endif
        tr_sum(ip,ie,ia,it,icf_sum(ix,iy,iz),im,1) = Xa
        tr_sum(ip,ie,ia,it,icf_sum(ix,iy,iz),im,2) = sigx

      else

        tr_sum(ip,ie,ia,it,icf_sum(ix,iy,iz),im,2) = 0.0

      end if

      return
      end


************************************************************************
* test                                                                     *
      subroutine  psufxyz_sumover_get(m,iax,
     &    npstepi,nestepi,nastepi,ntstepi,nxstepi,nystepi,nzstepi,
     &    nmstepi,
     &    ipi,iei,iai,iti,ixi,iyi,izi,imi,maxtott,tott)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0

      implicit double precision (a-h,o-z)

      parameter (ndim_p = 8)

      real(8) :: tott(maxtott,2)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

      dimension :: nstep(ndim_p),istart(ndim_p),ibase(ndim_p),
     &             n_sum_l(ndim_p),n_sum(ndim_p)

      real(8),pointer :: tr_sum(:)

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if

      ntf = 1
      ndim       = ndim_p
      istart(:)  = 1
      n_sum_l(:) = 1

      call set_integer_8(nstep,
     &    npstepi,nestepi,nastepi,ntstepi,nxstepi,nystepi,nzstepi,
     &    nmstepi)

      call set_integer_8(ibase,
     &      ipi,iei,iai,iti,ixi,iyi,izi,imi)

      call set_integer_8(n_sum,
     &    itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &    ittnm_sum(m,iax),
     &    itxnm_sum(m,iax),itynm_sum(m,iax),itznm_sum(m,iax)+1,
     &    itmst_sum(m,iax))

      call tr_sum_to_tott(ntf, ndim, nstep, istart, ibase,
     &                    n_sum_l,n_sum, maxtott, tott, tr_sum)
     
      return
      end



************************************************************************
*                                                                      *
      subroutine  psufxyz_sumover_getput(mode,m,iax,
     &    npstepi,nestepi,nastepi,ntstepi,nxstepi,nystepi,nzstepi,
     &    nmstepi,
     &    ipi,iei,iai,iti,ixi,iyi,izi,imi,maxtott,tott_sum)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      real(8) :: tott_sum(maxtott,2)

      real(8),allocatable :: tott_in(:,:)
      character :: chin*200

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

      real(8),pointer :: tr_sum(:)

      if(mode == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if

        call psufxyz_sumover_getput_sub(mode,
     &    npstepi,nestepi,nastepi,ntstepi,nxstepi,nystepi,nzstepi,
     &    nmstepi,
     &    ipi,iei,iai,iti,ixi,iyi,izi,imi,maxtott,tott_sum,
     &    itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &    ittnm_sum(m,iax),
     &    itxnm_sum(m,iax),itynm_sum(m,iax),itznm_sum(m,iax),
     &    itmst_sum(m,iax),tr_sum)


      else

          tr_sum => trRES_sum(irestalm_sum(m,iax):)

        nsame = npstepi * nestepi * nastepi * ntstepi
     &        * nxstepi * nystepi * (nzstepi + 1) * nmstepi
        allocate(tott_in(nsame,2))

        read(mode,'(a)') chin
        read(mode,'(26x,1000(1pe13.4,0pf8.4))')
     &  (tott_in(i,1),tott_in(i,2),i=1,nsame)

          call psufxyz_sumover_getput_sub(mode,
     &    npstepi,nestepi,nastepi,ntstepi,nxstepi,nystepi,nzstepi,
     &    nmstepi,
     &    ipi,iei,iai,iti,ixi,iyi,izi,imi,nsame,tott_in,
     &    itpan_sum(m,iax),itenm_sum(m,iax),itanm_sum(m,iax),
     &    ittnm_sum(m,iax),
     &    itxnm_sum(m,iax),itynm_sum(m,iax),itznm_sum(m,iax),
     &    itmst_sum(m,iax),tr_sum)

        deallocate(tott_in)

      endif


      return
      end

************************************************************************
*                                                                      *
      subroutine psufxyz_sumover_getput_sub(mode,
     &    npstepi,nestepi,nastepi,ntstepi,nxstepi,nystepi,nzstepi,
     &    nmstepi,
     &    ipi,iei,iai,iti,ixi,iyi,izi,imi,maxtott,tott_sum,
     &    np_sum,ne_sum,na_sum,nt_sum,nx_sum,ny_sum,nz_sum,nm_sum,
     &    tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      real(8) :: tr_sum(np_sum,ne_sum,na_sum,nt_sum,
     &                  nx_sum*ny_sum*(nz_sum+1),
     &                  nm_sum,2)

      real(8) :: tott_sum(maxtott,2)

      icf_sum(ix,iy,iz) = ix + ( iy - 1 ) * nx_sum
     &                  + ( iz - 1 ) * nx_sum * ny_sum

        i = 0
        do imloop=1,nmstepi
         do izloop=1,nzstepi
          do iyloop=1,nystepi
           do ixloop=1,nxstepi
            do itloop=1,ntstepi
              do ialoop=1,nastepi
               do ieloop=1,nestepi
                do iploop=1,npstepi
                  i = i + 1
                  do k=1,2
                   if(mode == 0) then
                      tott_sum(i,k) =
     &              tr_sum(ipi+iploop-1,iei+ieloop-1,iai+ialoop-1,
     &              iti+itloop-1,
     &              icf_sum(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),
     &              imi+imloop-1,k)
                   else
                    tr_sum(ipi+iploop-1,iei+ieloop-1,iai+ialoop-1,
     &              iti+itloop-1,
     &              icf_sum(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),
     &              imi+imloop-1,k)
     &              = tott_sum(i,k)
                   endif
                  enddo
                end do
              end do
             end do
            end do
          end do
         end do
        end do
      end do

      return
      end
