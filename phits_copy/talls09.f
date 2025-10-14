************************************************************************
*                                                                      *
      subroutine tadjntreg(ncol,m,nl,lt,np,nr,mr,ne,nm,nt,kr,eb,tb,tr,
     &                    trEVENT)
*                                                                      *
*       adjoint tally in region mesh                                   *
*       last modified by K. Niita on 2016/01/22                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*              0 : end of batch                                        *
*              4 : source particle                                     *
*                                                                      *
*              9 : termination by time cut-off                         *
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
!$    use omp_lib
*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

      parameter( rlit = 29.97925d0 )

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)

*-----------------------------------------------------------------------

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)

*-----------------------------------------------------------------------
!sumover
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)

      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall66/ itadm(itlmax), rtade(itlmax), rtadw(itlmax)

      common /tall82/ itcnth(9,itlmax)

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)

      common /enginit/ engini
!$OMP THREADPRIVATE(/enginit/)

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension   kr(mr)
      dimension   eb(ne+1)
      dimension   tb(nt+1)
      dimension   tr(np,ne,nt,nr,nm,2)
      dimension   trEVENT(np,ne,nt,nr,nm)       !OBINATA(2012.5.29): as Ct
      real(8),allocatable,save:: tr0(:,:,:,:,:) !OBINATA(2012.5.29): as C

      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt
      dimension facm(6)

      common /stat / istdev, irestart, ireschk
      common /cparm/ maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      common /egs5cmn10/denstepold,denstepnew,deinit
      real*8            denstepold,denstepnew,deinit
!$OMP THREADPRIVATE(/egs5cmn10/)


      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

*-----------------------------------------------------------------------
cKN 2016/01/22  initial multi-source number
*               energy bin
*-----------------------------------------------------------------------

         jsos = nsos(ibksos+1,ipomp+1)

         eaini = engini
         eamin = rtade(m)
         eamax = rtadw(m)

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
* OBINATA(2012.5.29): change tr(,,,,,3) to trEVENT(,,,,)
*-----------------------------------------------------------------------

         if (( ncol .eq. 0 .or. ncol .eq. 4 )
     &                              .and. istdev .eq. 2) then
           if ((nocas.gt.1.or.ncol.eq.0) .and. ihistcount.ne.1 ) then

             tr(:,:,:,:,:,1) = tr(:,:,:,:,:,1) + trEVENT(:,:,:,:,:)
             tr(:,:,:,:,:,2) = tr(:,:,:,:,:,2) + trEVENT(:,:,:,:,:) ** 2
! sumover
             call ttracreg_sumover(m,1,
     &                   np,    ne,    nt,    nr,    nm,    trEVENT)

           end if

           trEVENT(:,:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------
*        end of batch ( in case of istdev = 1 )
*
* OBINATA(2012.5.29): modificate for thread parallel
*-----------------------------------------------------------------------

         if ( ncol .eq. 0 .and. istdev .eq. 1) then
!$OMP MASTER
             allocate( tr0(np,ne,nt,nr,nm) )
             tr0(:,:,:,:,:) = 0.d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP CRITICAL (tadjntreg_crit_ist1)
             tr0(:,:,:,:,:) = tr0(:,:,:,:,:) + trEVENT(:,:,:,:,:)
!$OMP END CRITICAL (tadjntreg_crit_ist1)
!$OMP BARRIER
!$OMP MASTER
             tr(:,:,:,:,:,1) = tr(:,:,:,:,:,1) + tr0(:,:,:,:,:) / maxcas
             tr(:,:,:,:,:,2) = tr(:,:,:,:,:,2)
     &                     + ( tr0(:,:,:,:,:) / maxcas ) ** 2

c sumover
c tr0(np,ne,nt,nr,nm)
             call ttracreg_sumover(m,maxcas,
     &                   np,    ne,    nt,    nr,    nm,    tr0)

             deallocate( tr0 )
!$OMP END MASTER

           trEVENT(:,:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------
*        check of ncol
*-----------------------------------------------------------------------

         if( ncol .lt. 9 ) return

*-----------------------------------------------------------------------
*        check of energy
*-----------------------------------------------------------------------

            if( e(ibke+no,ipomp+1)   .lt. eamin ) goto 999
            if( ec(ibkec+no,ipomp+1) .ge. eamax ) goto 999

            if( eaini .lt. eb(1) ) goto 999
            if( eaini .ge. eb(ne+1) ) goto 999

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncnt(ibknct+i,no,ipomp+1) .lt. itcnt(i*2+2,m) .or.
     &                ncnt(ibknct+i,no,ipomp+1) .gt. itcnt(i*2+3,m) )
     &                   return

               end if

            end do

*-----------------------------------------------------------------------
*        check of mat
*-----------------------------------------------------------------------

            if( nl .gt. 0 ) then

                  do i = 1, nl

                     if( itmcn(m) .gt. 0 .and.
     &                   idmn(mat) .eq. lt(i) ) goto 502
                     if( itmcn(m) .lt. 0 .and.
     &                   idmn(mat) .eq. lt(i) ) return

                  end do

                     if( itmcn(m) .gt. 0 ) return

            end if

  502          continue

*-----------------------------------------------------------------------
*        check of particles
*-----------------------------------------------------------------------

               if( jsos .gt. np ) return

*-----------------------------------------------------------------------
*           check of time
*-----------------------------------------------------------------------

               tparti = abs(t(ibkt+no,ipomp+1))
               tpartf = abs(tc(ibktc+no,ipomp+1))

               if( tparti .ge. tb(nt+1) ) goto 999
               if( tpartf .lt. tb(1) ) goto 999

*-----------------------------------------------------------------------
*           check region
*-----------------------------------------------------------------------

                  jj = 0

               do ii = 1, nr

                  call tregck(iblz1,ilev1,ilat1,mr,kr,jj,icc)

                  if( icc .ne. 0 ) goto 501

              end do

                  goto 999

  501         continue

*-----------------------------------------------------------------------
*           track length and initial time
*-----------------------------------------------------------------------

               tlngth = sqrt( ( xc(ibkxc+no,ipomp+1) -
     &                                       x(ibkx+no,ipomp+1) )**2
     &                      + ( yc(ibkyc+no,ipomp+1) -
     &                                       y(ibky+no,ipomp+1) )**2
     &                      + ( zc(ibkzc+no,ipomp+1) -
     &                                       z(ibkz+no,ipomp+1) )**2 )

*-----------------------------------------------------------------------
*           initial and final energy
*-----------------------------------------------------------------------

                  se  = eaini
                  ee  = eaini

               do i = 1, ne

                  if( se .ge. eb(i) .and.
     &                se .lt. eb(i+1) ) goto 30

               end do

   30             ie1 = min( i, ne )

               do i = ne, 1, -1

                  if( ee .ge. eb(i) .and.
     &                ee .lt. eb(i+1) ) goto 40

               end do

   40             ie2 = max( i, 1 )

*-----------------------------------------------------------------------
*           rrl : track length and rrr : charge particles correction
*-----------------------------------------------------------------------

               rrl = tlngth
               rrr = tlngth
               rrt = 0.0d0
               tpc = tparti

               rrl0 = rrl
               rrr0 = rrr
               rrt0 = rrt
               tpi0 = tpc

*-----------------------------------------------------------------------
*     check of region
*-----------------------------------------------------------------------

                  jj = 0

*-----------------------------------------------------------------------

      do 100 ii = 1, nr

*-----------------------------------------------------------------------

               call tregck(iblz1,ilev1,ilat1,mr,kr,jj,icc)

               if( icc .eq. 0 ) goto 100

               ir = ii

               icli = idgr(iblz1)

*-----------------------------------------------------------------------

               tpi = tpi0
               rrl = rrl0
               rrr = rrr0
               rrt = rrt0

*-----------------------------------------------------------------------
*        tally
*-----------------------------------------------------------------------

            do 270 ie = ie1, ie2, -1

                     emie = max( ee, eb(ie) )
                     emiu = min( se, eb(ie+1) )
                     erg  = ( emiu + emie ) / 2.0

                     call fmfac(m,icli,erg,oldwt,facm)

*-----------------------------------------------------------------------
*              initial and final time
*-----------------------------------------------------------------------

               dst = rrl
               tmd = 0.0d0
               CC  = 0.0d0

               if( erg .gt. 0.1 .or. rtyp .eq. 0.0d0 ) then

                  tmd = dst * ( erg + rtyp )
     &                / sqrt( erg * ( erg + 2.0 * rtyp ) )
     &                / rlit
                  CC  = 1.d0 / ( erg + rtyp )
     &                * sqrt( erg * ( erg + 2.0 * rtyp ) )
     &                * rlit

               else if( erg .gt. 0.0 ) then

                  tmd = dst * sqrt( rtyp / 2.0 / erg ) / rlit
                  CC = 1.d0 / sqrt( rtyp / 2.0 / erg ) * rlit

               end if

               tpf = tpi + tmd

               do 271 it = 1, nt

                 if( tpf .lt. tb(it) .or. tpi .ge. tb(it+1) ) goto 271

                 tin = max( tpi, tb(it) )
                 tfn = min( tpf, tb(it+1) )
                 tmd = max( 0.0d0, tfn - tin )

                 dst = tmd * CC

*-----------------------------------------------------------------------

                  tlv = oldwt * dst
                  ip = jsos

                  tlv = tlv * (eb(ne+1)-eb(1)) / (eamax-eamin) ! T.Sato 2020/09/06

                  do im = 1, nm

                     trEVENT(ip,ie,it,ir,im) =
     &               trEVENT(ip,ie,it,ir,im) + tlv * facm(im)

                  end do

*-----------------------------------------------------------------------

  271          continue

                  rrr = rrr - rrl
                  tpi = tpf

  270       continue

*-----------------------------------------------------------------------

  100 continue

*-----------------------------------------------------------------------

  999 continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine padjntreg(m,np,nr,mr,ne,nm,nt,kr,eb,tb,tr,
     &                    nvl,ivl,rvl,
     &                    nx,ny,nz,xm,ym,zm,igsh,idasa)
*                                                                      *
*       output the adjoint tally in region mesh                        *
*       last modified by K. Niita on 2016/01/22                        *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      common /stat / istdev, irestart, ireschk
      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
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
      common /tall29/ itrsh(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)

      common /istcut/ ist_cut, ist_bat

      common /fact01/ facmax(itlmax) ! kitamura23/03/31

*-----------------------------------------------------------------------

      common /volreg/ dvol(kvlmax)

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

*-----------------------------------------------------------------------

      dimension   kr(mr)
      dimension   eb(ne+1)
      dimension   ew(ne)
      dimension   tb(nt+1)
      dimension   tw(nt)
      dimension   vl(nr)
      dimension   lr(nr)
      dimension   tr(np,ne,nt,nr,nm,2)
      dimension   ivl(nvl)
      dimension   rvl(nvl)

      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)
      dimension   val(nr)

      integer,allocatable :: ixyz(:)

*-----------------------------------------------------------------------

      dimension tott(np,2) ! frtati 2021/10/05 6 -> np

      character hsunit(14)*29

      data hsunit( 1) / '[1/cm^2/source]              '/
      data hsunit( 2) / '[1/cm^2/MeV/source]          '/
      data hsunit( 3) / '[1/cm^2/Lethargy/source]     '/
      data hsunit( 4) / '[cm/source]                  '/
      data hsunit(11) / '[1/cm^2/nsec/source]         '/
      data hsunit(12) / '[1/cm^2/nsec/MeV/source]     '/
      data hsunit(13) / '[1/cm^2/nsec/Lethargy/source]'/
      data hsunit(14) / '[cm/nsec/source]             '/

      character cha*1
      data cha /"'"/

cfrtati 2021/10/05 chl, chm moved to partmod

      character chp(6)*10
      character chq(6)*8


      character rpa*1
      data rpa /'}'/
      character yen*1

      include 'samepage_include/samepage000.inc'

      include 'samepage_include/samepage001.inc'



      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

*-----------------------------------------------------------------------
*     for igshow
*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

            npg = np
            neg = ne
            ntg = nt

         else

            npg = 1
            neg = 1
            ntg = 1

         end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               chp(1) = 'y(msos-1  '
               chp(2) = 'y(msos-2  '
               chp(3) = 'y(msos-3  '
               chp(4) = 'y(msos-4  '
               chp(5) = 'y(msos-5  '
               chp(6) = 'y(msos-6  '
               chq(1) = '1       '
               chq(2) = '2       '
               chq(3) = '3       '
               chq(4) = '4       '
               chq(5) = '5       '
               chq(6) = '6       '

*-----------------------------------------------------------------------
*           itunt(m) = 1,4 11,14: /cm^2/source
*                                 /cm^2/sec/source
*                    = 2   12   : /cm^2/MeV/source
*                                 /cm^2/MeV/nsec/source
*                    = 3   13   : /cm^2/Lethargy/source
*                                 /cm^2/Lethargy/nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  4 .or.
     &          itunt(m) .eq. 11 .or. itunt(m) .eq. 14 ) then

               do i = 1, ne

                  ew(i) = 1.d+0

               end do
               ew_sum = 1.0

            else if( itunt(m) .eq. 2 .or. itunt(m) .eq. 12 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) -eb(1)

            else if( itunt(m) .eq. 3 .or. itunt(m) .eq. 13 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log(eb(ne+1) -eb(1))

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
               tw_sum = 1.0

            end if

*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

               call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)
               vl_sum = sum(vl(1:nr))

*-----------------------------------------------------------------------
*        ( unit = 4 or 14  ; vol = 1.0 )
*-----------------------------------------------------------------------

            if( itunt(m) .eq. 4 .or. itunt(m) .eq. 14 ) then

               do ir = 1, nr

                  vl(ir) = 1.0d0

               end do
               vl_sum = 1.0

            end if

*-----------------------------------------------------------------------
*        relative error and  unit conversion
*        c1 : nomalization for source
*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

                  cmax = 0.0
                  cmin = 1.e+33
                  dnon = 1.e-33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

                  c1 = 1.0d+0 / rsouin

            facmax(m) = 1.d0
            if( rtfac(m) .lt. 0.d0 ) then
               facmax(m) = 0.d0

               do 101 im = 1, nm
               do 101 ir = 1, nr
               do 101 it = 1, nt
               do 101 ie = 1, ne
               do 101 ip = 1, np

                  if( tr(ip,ie,it,ir,im,1) .gt. 0.d0 ) then

                     fmaxfc = tr(ip,ie,it,ir,im,1)
     &                                   / vl(ir) / ew(ie) / tw(it)

                     if( facmax(m) .lt. fmaxfc) facmax(m) = fmaxfc

                  end if

  101          continue

            end if

            do 100 im = 1, nm
            do 100 ir = 1, nr
            do 100 it = 1, nt
            do 100 ie = 1, ne
            do 100 ip = 1, np

               if( tr(ip,ie,it,ir,im,1) .gt. 0.d0 ) then

                  call calc_stdev(m,Xa,sigx,
     &                            tr(ip,ie,it,ir,im,1),
     &                            tr(ip,ie,it,ir,im,2),
     &            abs(rtfac(m)/facmax(m))/vl(ir)/ew(ie)/tw(it))

                  tr(ip,ie,it,ir,im,1) = Xa
                  tr(ip,ie,it,ir,im,2) = sigx

                  if( sigx .gt. stdm ) stdm = sigx
                  if( tr(ip,ie,it,ir,im,1) .gt. cmax )
     &                        cmax = tr(ip,ie,it,ir,im,1)

                  if( tr(ip,ie,it,ir,im,1) .lt. cmin )
     &                        cmin = tr(ip,ie,it,ir,im,1)

               else

                  isdz = 1
                  tr(ip,ie,it,ir,im,2) = 0.0

               end if

! sumover
               fact_in = abs(rtfac(m)/facmax(m))
               call ptracreg_sumover_stdev(0,m,ip,ie,it,ir,im,
     &                           fact_in,ew(ie),tw(it),vl(ir),
     &                                  ew_sum,tw_sum,vl_sum)

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

         if( igsh .ne. 0 .and.
     &     ( itaxs(m,iax) .lt. 7 .or. itrsh(m) .eq. 0 ) ) goto 900

         if( (itall .eq. 2 .and. nobch .lt. maxbch .and. igsh .eq. 0)
     &        .or. ( itall .eq. 4 .and. igsh .eq. 0 ) ) then

            call mk_2dnumfn(ctfln(m,iax),fname,itfll(m,iax),
     &           nobch,maxbch,npe)

         else

            fname = ctfln(m,iax)

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

               call tadjech(iot,m,iax,1)


      include 'samepage_include/samepage002_petrm.inc'
      include 'samepage_include/samepagechp_petrm.inc'
      include 'samepage_include/samepageseti.inc'

*-----------------------------------------------------------------------
*        energy axis
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 1 ) then

               inum = 0

            do imi = 1, nm
            do iri = 1, nr, nrstepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ir = iri
               it = iti
               ip = ipi


               inum = inum + 1

               ireg = lr(ir)

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

               write(iot,'("#   no. =",i3,3x,"reg =",i7)')
     &                     inum, ireg
               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3,/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Energy [MeV]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Flux ",a29)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itunt(m) .eq. 3 .or. itunt(m) .eq. 13 .or.
     &             itety(m) .eq. 3 .or. itety(m) .eq. 5 ) then

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

       include 'samepage_include/flux_petrm_erg.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( itmlp(m) .eq. 0 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  reg  =",i7,a1)')
     &                        cha, inum, ireg, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  reg  =",i7,",  it =",i3,a1)')
     &                        cha, inum, ireg, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  reg =",i7,
     &                        ",  mset =",i4,a1)')
     &                        cha, inum, ireg, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  reg =",i7,
     &                        ",  mset =",i3,",  it =",i3,a1)')
     &                        cha, inum, ireg, itmnt(m,im), it, cha
                  end if
               end if
             else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                write(changelsub(4),'(",  reg  =",i7)') ireg
                write(changelsub(5),'(",  it =",i3)') it
                write(changelsub(6),'(",  ie =",i3)') ie
                write(changelsub(7),'(a1)') cha
               changelsub(6) = " "
               if( itmlp(m) .eq. 0 ) then
                  changelsub(3) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 2) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(6) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))
                write(iot,'(/a)') trim(angeltitle)
             end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]")')
     &                     yen, vl(ir)
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
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

*-----------------------------------------------------------------------
*        reg axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 2 ) then

               inum = 0

            do imi = 1, nm
            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3/
     &            "#  ie =",i3,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie,
     &                        eb(ie), eb(ie+nestepi)
               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3,/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Serial Num. of Region")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Flux ",a29)') hsunit(itunt(m))

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

       include 'samepage_include/flux_petrm_reg.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( itmlp(m) .eq. 0 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,",  ie =",
     &                     i3,a1)')
     &                     cha, inum, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,",  ie =",
     &                     i3,",  it =",i3,a1)')
     &                     cha, inum, ie, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,",  ie =",
     &                     i3,",  mset =",i3,a1)')
     &                     cha, inum, ie, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,",  ie =",
     &                     i3,",  mset =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, itmnt(m,im), it, cha
                  end if
               end if

             else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                write(changelsub(4),'(",  reg  =",i7)') ireg
                write(changelsub(5),'(",  it =",i3)') it
                write(changelsub(6),'(",  ie =",i3)') ie
                write(changelsub(7),'(a1)') cha
               changelsub(4) = " "
               if( itmlp(m) .eq. 0 ) then
                  changelsub(3) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 2) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(6) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))
                write(iot,'(/a)') trim(angeltitle)
             end if


               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     yen, eb(ie), eb(ie+nestepi)
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
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

*-----------------------------------------------------------------------
*        time axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 11 ) then

               inum = 0

            do imi = 1, nm
            do iri = 1, nr, nrstepi
            do iei = 1, ne, nestepi
            do ipi = 1, np, npstepi
               im = imi
               ir = iri
               ie = iei
               ip = ipi

               inum = inum + 1

               ireg = lr(ir)

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

               write(iot,'("#   no. =",i3,3x,"reg =",i7)')
     &                     inum, ireg
               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)
                  write(iot,'(
     &            "#  ie =",i3,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ie, eb(ie), eb(ie+nestepi)

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Time [nsec]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Flux ",a29)') hsunit(itunt(m))

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

       include 'samepage_include/flux_petrm_time.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( itmlp(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  reg  =",i7,",  ie =",i3,a1)')
     &                        cha, inum, ireg, ie, cha
               else
                     write(iot,'(/a1,"no. =",i3,
     &                        ",  reg =",i7,
     &                        ",  mset =",i3,",  ie =",i3,a1)')
     &                        cha, inum, ireg, itmnt(m,im), ie, cha
               end if

             else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                write(changelsub(3),'(",  mset =",i4)') itmnt(m,im)
                write(changelsub(4),'(",  reg  =",i7)') ireg
                write(changelsub(5),'(",  it =",i3)') it
                write(changelsub(6),'(",  ie =",i3)') ie
                write(changelsub(7),'(a1)') cha
               changelsub(5) = " "
               if( itmlp(m) .eq. 0 ) then
                  changelsub(3) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 2) then
                  changelsub(4) = " "
               end if
               if(iloopmode .eq. 10) then
                  changelsub(5) = " "
               end if
               if(iloopmode .eq. 1) then
                  changelsub(6) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))
                write(iot,'(/a)') trim(angeltitle)
             end if



               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]")')
     &                     yen, vl(ir)
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                     write(iot,'("e:")')

            end do

            end do
            end do
            end do

*-----------------------------------------------------------------------
*        xy axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 7 ) then

*-----------------------------------------------------------------------

               inum = 0

            do im = 1, nm
            do iz = 1, nz
            do ip = 1, npg
            do ie = 1, neg
            do it = 1, ntg

               zval = ( zm(iz) + zm(iz+1) ) / 2.0d0

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

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "iz  =",i3,3x,
     &            "msos = ",a8,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = ",1p1e13.4)')
     &                     inum, ie, iz, chq(ip),
     &                     eb(ie), eb(ie+1), zval
               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3,/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+1)
               end if

*-----------------------------------------------------------------------

               if( itmlp(m) .eq. 0 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ie, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, iz, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, iz, itmnt(m,im), it, cha
                  end if
               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

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

         if( igsh .eq. 0 ) then

               do i = 1, nr

                  val(i) = tr(ip,ie,it,i,im,1)

               end do

                  izlog = 1
                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

            if( izlog .eq. 0 ) then

               do i = 1, nr

                  if( val(i) .gt. dnon .and. cmin .ne. cmax ) then

                     val(i) = ( val(i) - cmin ) / ( cmax - cmin )
                     val(i) = max(0.0d0,val(i))
                     val(i) = min(1.0d0,val(i))

                  else

                     val(i) = -1.0

                  end if

               end do

            else

               do i = 1, nr

                  if( val(i) .gt. dnon .and. cmin .ne. cmax ) then

                     val(i) = log10(val(i)/cmin) / log10(cmax/cmin)
                     val(i) = max(0.0d0,val(i))
                     val(i) = min(1.0d0,val(i))

                  else

                     val(i) = -1.0d0

                  end if

               end do

            end if

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

               if( cmin .gt. 0.0 .and. cmax .gt. cmin ) then

                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax

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

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

*-----------------------------------------------------------------------
*        rshow
*-----------------------------------------------------------------------

               write(iot,'(/"#",78("-"))')
               write(iot,'( "# rshow")')
               write(iot,'( "#",78("-"))')
               write(iot,'(/"p: legs[c5*0.875]")')

               none = 1
               iaxs = 1
               iuni = itrsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

            if( igsh .eq. 0 ) then

               call gshow(1,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nx+1,ny+1,none,xm,ym,zval,ixyz(1),
     &                    nr,mr,kr,val,itmtr(m,4),itglt(m))

            else

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nx+1,ny+1,none,xm,ym,zval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. igsh .eq. 0 .and.
     &    cmin .gt. 0.0 .and. cmax .gt. cmin ) then

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

         write(iot,'("y: Flux ",a29)') hsunit(itunt(m))

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
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  z     &=&",1pe13.4," [cm]"/
     &                     "  msos &=& ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     zval, chq(ip)
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
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

      if( inolg .eq. 1 .and. igsh .eq. 0 .and.
     &    cmin .gt. 0.0 .and. cmax .gt. cmin ) then

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

*-----------------------------------------------------------------------
*        yz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 ) then

*-----------------------------------------------------------------------

               inum = 0

            do im = 1, nm
            do ix = 1, nx
            do ip = 1, npg
            do ie = 1, neg
            do it = 1, ntg

               xval = ( xm(ix) + xm(ix+1) ) / 2.0d0

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

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "ix  =",i3,3x,
     &            "msos = ",a8,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   x = ",1p1e13.4)')
     &                     inum, ie, ix, chq(ip),
     &                     eb(ie), eb(ie+1), xval
               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3,/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+1)
               end if

*-----------------------------------------------------------------------

               if( itmlp(m) .eq. 0 ) then
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,a1)')
     &                     cha, inum, ie, ix, cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, ix, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ix, itmnt(m,im), cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  mset =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, ix, itmnt(m,im), it, cha
                  end if
               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

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

         if( igsh .eq. 0 ) then

               do i = 1, nr

                  val(i) = tr(ip,ie,it,i,im,1)

               end do

                  izlog = 1
                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

            if( izlog .eq. 0 ) then

               do i = 1, nr

                  if( val(i) .gt. dnon .and. cmin .ne. cmax ) then

                     val(i) = ( val(i) - cmin ) / ( cmax - cmin )
                     val(i) = max(0.0d0,val(i))
                     val(i) = min(1.0d0,val(i))

                  else

                     val(i) = -1.0

                  end if

               end do

            else

               do i = 1, nr

                  if( val(i) .gt. dnon .and. cmin .ne. cmax ) then

                     val(i) = log10(val(i)/cmin) / log10(cmax/cmin)
                     val(i) = max(0.0d0,val(i))
                     val(i) = min(1.0d0,val(i))

                  else

                     val(i) = -1.0d0

                  end if

               end do

            end if

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

               if( cmin .gt. 0.0 .and. cmax .gt. cmin ) then

                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax

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

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

*-----------------------------------------------------------------------
*        rshow
*-----------------------------------------------------------------------

               write(iot,'(/"#",78("-"))')
               write(iot,'( "# rshow")')
               write(iot,'( "#",78("-"))')
               write(iot,'(/"p: legs[c5*0.875]")')

               none = 1
               iaxs = 2
               iuni = itrsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

            if( igsh .eq. 0 ) then

               call gshow(1,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,ny+1,none,zm,ym,xval,ixyz(1),
     &                    nr,mr,kr,val,itmtr(m,4),itglt(m))

            else

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,ny+1,none,zm,ym,xval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. igsh .eq. 0 .and.
     &    cmin .gt. 0.0 .and. cmax .gt. cmin ) then

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

         write(iot,'("y: Flux ",a29)') hsunit(itunt(m))

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
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  x     &=&",1pe13.4," [cm]"/
     &                     "  msos &=& ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     xval, chq(ip)
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
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

      if( inolg .eq. 1 .and. igsh .eq. 0 .and.
     &    cmin .gt. 0.0 .and. cmax .gt. cmin ) then

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

*-----------------------------------------------------------------------
*        xz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9 ) then

*-----------------------------------------------------------------------

               inum = 0

            do im = 1, nm
            do iy = 1, ny
            do ip = 1, npg
            do ie = 1, neg
            do it = 1, ntg

               yval = ( ym(iy) + ym(iy+1) ) / 2.0d0

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

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "iy  =",i3,3x,
     &            "msos = ",a8,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = ",1p1e13.4)')
     &                        inum, ie, iy, chq(ip),
     &                        eb(ie), eb(ie+1), yval
               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3,/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+1)
               end if

*-----------------------------------------------------------------------

               if( itmlp(m) .eq. 0 ) then
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, ie, iy, cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, iy, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, iy, itmnt(m,im), cha
                  else
                        write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, iy, itmnt(m,im), it, cha
                  end if
               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

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

         if( igsh .eq. 0 ) then

               do i = 1, nr

                  val(i) = tr(ip,ie,it,i,im,1)

               end do

                  izlog = 1
                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

            if( izlog .eq. 0 ) then

               do i = 1, nr

                  if( val(i) .gt. dnon .and. cmin .ne. cmax ) then

                     val(i) = ( val(i) - cmin ) / ( cmax - cmin )
                     val(i) = max(0.0d0,val(i))
                     val(i) = min(1.0d0,val(i))

                  else

                     val(i) = -1.0

                  end if

               end do

            else

               do i = 1, nr

                  if( val(i) .gt. dnon .and. cmin .ne. cmax ) then

                     val(i) = log10(val(i)/cmin) / log10(cmax/cmin)
                     val(i) = max(0.0d0,val(i))
                     val(i) = min(1.0d0,val(i))

                  else

                     val(i) = -1.0d0

                  end if

               end do

            end if

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

               if( cmin .gt. 0.0 .and. cmax .gt. cmin ) then

                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax

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

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

*-----------------------------------------------------------------------
*        rshow
*-----------------------------------------------------------------------

               write(iot,'(/"#",78("-"))')
               write(iot,'( "# rshow")')
               write(iot,'( "#",78("-"))')
               write(iot,'(/"p: legs[c5*0.875]")')

               none = 1
               iaxs = 3
               iuni = itrsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

            if( igsh .eq. 0 ) then

               call gshow(1,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,nx+1,none,zm,xm,yval,ixyz(1),
     &                    nr,mr,kr,val,itmtr(m,4),itglt(m))

            else

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,nx+1,none,zm,xm,yval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. igsh .eq. 0 .and.
     &    cmin .gt. 0.0 .and. cmax .gt. cmin ) then

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

         write(iot,'("y: Flux ",a29)') hsunit(itunt(m))

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
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  y     &=&",1pe13.4," [cm]"/
     &                     "  msos &=& ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     yval, chq(ip)
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
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

      if( inolg .eq. 1 .and. igsh .eq. 0 .and.
     &    cmin .gt. 0.0 .and. cmax .gt. cmin ) then

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

*-----------------------------------------------------------------------

         end if

            call prestart(m,iot) !OBINATA(2012.6.13)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

*-----------------------------------------------------------------------

  900 continue

*-----------------------------------------------------------------------
      deallocate( ixyz )
      include 'samepage_include/samepage999.inc'


      return
      end


************************************************************************
*                                                                      *

      subroutine tadjntrz(ncol,m,nl,lt,np,nr,nz,ne,nm,nt,na,
     &                   rm,zm,ab,eb,tb,tr,
     &                   trEVENT,itrmax,itrmin)
*                                                                      *
*       adjoint tally in r-z scoring mesh                              *
*       last modified by K.Niita on 2016/01/22                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*              9 : termination by time cut-off                         *
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

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)

*-----------------------------------------------------------------------

      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall08/ rtrx0(itlmax), rtry0(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /tall66/ itadm(itlmax), rtade(itlmax), rtadw(itlmax)

      common /tall82/ itcnth(9,itlmax)

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)

      common /enginit/ engini
!$OMP THREADPRIVATE(/enginit/)

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension rm(nr+1)
      dimension zm(nz+1)

      dimension ab(na+1)

      dimension eb(ne+1)
      dimension tb(nt+1)

      dimension tr(np,ne,nt,nr*nz,na,nm,2)

      dimension   trEVENT(np,ne,nt,nr*nz,na,nm)   !cKN 2014/12/23
      real(8),allocatable,save:: tr0(:,:,:,:,:,:) !OBINATA(2012.6.13): as C

      dimension itrmax(7),itrmin(7)

      data pi /3.1415926535897932384d0/

*-----------------------------------------------------------------------

      dimension ud(3)

      data dmax  /1.0d+19/
      data dmax0 /1.0d+18/

      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt
      dimension facm(6)

*-----------------------------------------------------------------------

      common /paraj/ mstz(300), parz(300)

*-----------------------------------------------------------------------

      common /stat / istdev, irestart, ireschk
      common /cparm/ maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

*-----------------------------------------------------------------------

         icf(ir,iz) = ir + ( iz - 1 ) * nr

*-----------------------------------------------------------------------

         small = parz(28) * 10.d0

*-----------------------------------------------------------------------
cKN 2016/01/22  initial multi-source number
*               energy bin
*-----------------------------------------------------------------------

         jsos = nsos(ibksos+1,ipomp+1)

         eaini = engini
         eamin = rtade(m)
         eamax = rtadw(m)

*-----------------------------------------------------------------------

      if (istdev .eq. 2) then

        call readitrminmax7(itrmin,itrmax,(/ np,ne,nt,nr,nz,na,nm/),
     &                      mnp,mne,mnt,mnr,mnz,mna,mnm,
     &                      mxp,mxe,mxt,mxr,mxz,mxa,mxm)

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
*OBINATA(2012.6.13)
*-----------------------------------------------------------------------

         if (( ncol .eq. 0 .or. ncol .eq. 4 )
     &                              .and. istdev .eq. 2) then
           if ((nocas.gt.1.or.ncol.eq.0) .and. ihistcount.ne.1 ) then

             do im = mnm,mxm
             do ia = mna,mxa
             do iz = mnz,mxz
             do ir = mnr,mxr
               irz = icf(ir,iz)
               tr(:,mne:mxe,:,irz,ia,im,1)
     &           = tr(:,mne:mxe,:,irz,ia,im,1)
     &           + trEVENT(:,mne:mxe,:,irz,ia,im)
               tr(:,mne:mxe,:,irz,ia,im,2)
     &           = tr(:,mne:mxe,:,irz,ia,im,2)
     &           + trEVENT(:,mne:mxe,:,irz,ia,im) ** 2
             enddo
             enddo
             enddo
             enddo

! sumover
             call ttracrz_sumover(m,1,
     &                   np, ne, nt, nr, nz, na, nm, trEVENT)

           end if


           do im = mnm,mxm
           do ia = mna,mxa
           do iz = mnz,mxz
           do ir = mnr,mxr
             irz = icf(ir,iz)
             trEVENT(:,mne:mxe,:,irz,ia,im) = 0
           enddo
           enddo
           enddo
           enddo

           call resetitrminmax(itrmin,itrmax,7,(/np,ne,nt,nr,nz,na,nm/))

         end if

*-----------------------------------------------------------------------
*        end of batch ( in case of istdev = 1 )
*OBINATA(2012.6.13)
*-----------------------------------------------------------------------

         if ( ncol .eq. 0 .and. istdev .eq. 1) then
!$OMP MASTER
             allocate( tr0(np,ne,nt,nr*nz,na,nm) )
             tr0(:,:,:,:,:,:) = 0.d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP CRITICAL (tadjntrz_crit)
             tr0(:,:,:,:,:,:) = tr0(:,:,:,:,:,:) + trEVENT(:,:,:,:,:,:)
!$OMP END CRITICAL (tadjntrz_crit)
!$OMP BARRIER
!$OMP MASTER
             tr(:,:,:,:,:,:,1) = tr(:,:,:,:,:,:,1)
     &                         + tr0(:,:,:,:,:,:) / maxcas
             tr(:,:,:,:,:,:,2) = tr(:,:,:,:,:,:,2)
     &                     + ( tr0(:,:,:,:,:,:) / maxcas ) ** 2

! sumover
             call ttracrz_sumover(m,maxcas,
     &                   np, ne, nt, nr, nz, na, nm, tr0)

             deallocate( tr0 )
!$OMP END MASTER

           trEVENT(:,:,:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------
*        check of ncol
*-----------------------------------------------------------------------

         if( ncol .lt. 9 ) return

*-----------------------------------------------------------------------
*        check of energy
*-----------------------------------------------------------------------

            if( e(ibke+no,ipomp+1)   .lt. eamin ) goto 999
            if( ec(ibkec+no,ipomp+1) .ge. eamax ) goto 999

            if( eaini .lt. eb(1) ) goto 999
            if( eaini .ge. eb(ne+1) ) goto 999

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncnt(ibknct+i,no,ipomp+1) .lt. itcnt(i*2+2,m) .or.
     &                ncnt(ibknct+i,no,ipomp+1) .gt. itcnt(i*2+3,m) )
     &                                         return

               end if

            end do

*-----------------------------------------------------------------------
*        check of mat
*-----------------------------------------------------------------------

            if( nl .gt. 0 ) then

                  do i = 1, nl

                     if( itmcn(m) .gt. 0 .and.
     &                   idmn(mat) .eq. lt(i) ) goto 502
                     if( itmcn(m) .lt. 0 .and.
     &                   idmn(mat) .eq. lt(i) ) return

                  end do

                     if( itmcn(m) .gt. 0 ) return

            end if

  502          continue

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

            if( zza-small .lt. zm(1) .and.
     &          zzc-small .lt. zm(1) ) goto 999

            if( zza+small .ge. zm(nz+1) .and.
     &          zzc+small .ge. zm(nz+1) ) goto 999

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

         if( itaty(m) .ne. 0 ) then

               xfa = xxa - x0
               yfa = yya - y0

               xfc = xxc - x0
               yfc = yyc - y0

               aveg = xfa * yfc - yfa * xfc

         end if

*-----------------------------------------------------------------------
*        check of particles
*-----------------------------------------------------------------------

               if( jsos .gt. np ) return

*-----------------------------------------------------------------------
*           check of time
*-----------------------------------------------------------------------

               tparti = abs(t(ibkt+no,ipomp+1))
               tpartf = abs(tc(ibktc+no,ipomp+1))

               if( tparti .ge. tb(nt+1) ) goto 999
               if( tpartf .lt. tb(1) ) goto 999

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
*     initial position, time, energy and initial range rng
*-----------------------------------------------------------------------

            tot = 0.0d0
            tpt = tparti

            se  = eaini
            xpp = xxa
            ypp = yya
            zpp = zza

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
                  if( zm(i) .gt. zpp + small ) then
                         izm = i - 1
                         izc = i - 1
                         goto 38
                  end if
               end do

            else

               izm = nz + 2
               izc = nz + 1

               do i = 1, nz + 1
                  if( zm(i) .ge. zpp - small ) then
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
*     initial a-position
*-----------------------------------------------------------------------

               iac = 1

*-----------------------------------------------------------------------
*     loop for finding mesh and booking upto total distance
*-----------------------------------------------------------------------

      do 50

*-----------------------------------------------------------------------
*      calculation is finished
*-----------------------------------------------------------------------

         if( tot .ge. dis ) return

*-----------------------------------------------------------------------
*        check of angle
*-----------------------------------------------------------------------

            if( itaty(m) .ne. 0 ) then

               dis0 = sqrt( ( xpp - x0 )**2
     &                    + ( ypp - y0 )**2 )

               if( dis0 .gt. 0.0d0 ) then
                     ang0 = acos( ( xpp - x0 ) / dis0 )
               else
                     ang0 = 0.0d0
               end if

               if( ypp - y0 .lt. 0.0d0 ) ang0 = 2.d0 * pi - ang0

               if( ang0 .eq. 0.d0 .and. aveg .lt. 0.d0 )
     &               ang0 = 2.d0 * pi

               if( ang0 .eq. 2.0 * pi .and. aveg .gt. 0.d0 )
     &               ang0 = 0.d0

               if( itaty(m) .lt. 0.0d0 ) then
                     ang0 = ang0 / pi * 180.d0
               end if

               do i = 1, na
                  if( aveg .ge. 0.d0 ) then
                     if( ang0 .ge. ab(i) - small .and.
     &                   ang0 .lt. ab(i+1) - small ) then
                           iac = i
                           goto 49
                     end if
                  else if( aveg .lt. 0.d0 ) then
                     if( ang0 .gt. ab(i) + small .and.
     &                   ang0 .le. ab(i+1) + small ) then
                           iac = i
                           goto 49
                     end if
                  end if
               end do

                  iac = -1

   49          continue

            end if

*-----------------------------------------------------------------------
*        dr : distance to the nearest r mesh
*-----------------------------------------------------------------------

               dr = dmax

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
*        da : distance to the nearest a mesh
*-----------------------------------------------------------------------

            da = dmax * 2.d0

         if( itaty(m) .ne. 0 .and.
     &       iac .gt. 0 .and. aveg .ne. 0 ) then

            if( aveg .gt. 0.d0 ) then
               iaa = iac + 1
            else if( aveg .lt. 0.d0 ) then
               iaa = iac
            end if

            thet = ab(iaa)
            if( itaty(m) .lt. 0.0d0 ) thet = thet / 180.d0 * pi

            xsn = -sin( thet )
            ysn =  cos( thet )

            xpc = xpp - x0
            ypc = ypp - y0

            agn1 = ( xfc * xsn + yfc * ysn )
            agn2 = ( xpc * xsn + ypc * ysn )

            agnm = agn1 * agn2

            if( agnm .lt. -small .or. abs(agn1) .lt. small ) then

               da = sqrt( ( xpp - xxc )**2
     &                  + ( ypp - yyc )**2
     &                  + ( zpp - zzc )**2 )
     &              * abs( agn2 ) / ( abs( agn1 ) + abs( agn2 ) )

            end if

         end if

*-----------------------------------------------------------------------
*        which boundary is the nearlist
*-----------------------------------------------------------------------

            if( da .lt. dr .and. da .lt. dz ) then

               dd = da
               jz = -1

            else if( dr .gt. dz ) then

               dd = dz
               jz = 1

            else

               dd = dr
               jz = 0

            end if

            if( dd .gt. dmax0 ) return

*-----------------------------------------------------------------------
*        propagate position upto the boundary or the final point
*-----------------------------------------------------------------------

            if( tot + dd .ge. dis - small ) then

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

               ee = se

*-----------------------------------------------------------------------
*        time evolution
*-----------------------------------------------------------------------

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

                  tpc = tpt
                  tpt = tpt + timd

*-----------------------------------------------------------------------
*        booking
*-----------------------------------------------------------------------

            iz = izc
            ir = irc
            ia = iac

         if( ir .ge. 1 .and. ir .lt. nr + 1 .and.
     &       iz .ge. 1 .and. iz .lt. nz + 1 .and.
     &       ia .ge. 1 .and. ia .lt. na + 1 ) then

*-----------------------------------------------------------------------
*           initial and final energy cell
*-----------------------------------------------------------------------

               do i = 1, ne

                  if( se .ge. eb(i) .and.
     &                se .lt. eb(i+1) ) goto 30

               end do

   30             ie1 = min( i, ne )

               do i = ne, 1, -1

                  if( ee .ge. eb(i) .and.
     &                ee .lt. eb(i+1) ) goto 40

               end do

   40             ie2 = max( i, 1 )

*-----------------------------------------------------------------------
*           rrl : track length and rrr : charge particles correction
*-----------------------------------------------------------------------

                  rrl = dd
                  rrr = dd
                  rrt = 0.0d0

                  tpi = tpc

*-----------------------------------------------------------------------
*        tally
*-----------------------------------------------------------------------

            do 270 ie = ie1, ie2, -1

                     emie = max( ee, eb(ie) )
                     emiu = min( se, eb(ie+1) )
                     erg  = ( emiu + emie ) / 2.0

                     call fmfac(m,icli,erg,oldwt,facm)

*-----------------------------------------------------------------------
*              initial and final time
*-----------------------------------------------------------------------

               dst = rrl
               tmd = 0.0d0
               CC  = 0.0d0

               if( erg .gt. 0.1 .or. rtyp .eq. 0.0d0 ) then

                  tmd = dst * ( erg + rtyp )
     &                / sqrt( erg * ( erg + 2.0 * rtyp ) )
     &                / rlit
                  CC  = 1.d0 / ( erg + rtyp )
     &                * sqrt( erg * ( erg + 2.0 * rtyp ) )
     &                * rlit

               else if( erg .gt. 0.0 ) then

                  tmd = dst * sqrt( rtyp / 2.0 / erg ) / rlit
                  CC = 1.d0 / sqrt( rtyp / 2.0 / erg ) * rlit

               end if

               tpf = tpi + tmd

               do 271 it = 1, nt

                 if( tpf .lt. tb(it) .or. tpi .ge. tb(it+1) ) goto 271

                 tin = max( tpi, tb(it) )
                 tfn = min( tpf, tb(it+1) )
                 tmd = max( 0.0d0, tfn - tin )

                 dst = tmd * CC

*-----------------------------------------------------------------------

                  tlv = oldwt * dst
                  ip = jsos

                  tlv = tlv * (eb(ne+1)-eb(1)) / (eamax-eamin) ! T.Sato 2020/09/06

                  do im = 1, nm

                trEVENT(ip,ie,it,icf(ir,iz),ia,im) =
     &          trEVENT(ip,ie,it,icf(ir,iz),ia,im) + tlv * facm(im)

                  end do

*-----------------------------------------------------------------------

  271          continue

                  rrr = rrr - rrl
                  tpi = tpf

  270       continue

         if (istdev .eq. 2) then

           call setitrmin(itrmin,2,7,(/ie2, 1,ir,iz, 1, 1/))
           call setitrmax(itrmax,2,7,(/ie1,it,ir,iz,ia,nm/))


         endif

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*     next position
*     for z-crossing jz = 1, r-crossing jz = 0, a-crossing jz = -1
*-----------------------------------------------------------------------

            if( jz .eq. -1 ) then

                  iac = iac

            else if( jz .eq. 1 ) then

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

   50 continue

*-----------------------------------------------------------------------

  999 continue

*-----------------------------------------------------------------------

      end


************************************************************************

      subroutine padjntrz(m,np,nr,nz,na,ne,nm,nt,rm,zm,ab,eb,
     &                   tb,tr,idasa)

*       output r-z scoring mesh adjoint tally                          *
*       last modified by K.Niita on 2016/01/22                         *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

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

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

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

*-----------------------------------------------------------------------

      dimension   rm(nr+1)
      dimension   zm(nz+1)
      dimension   eb(ne+1)
      dimension   ab(na+1)
      dimension   ew(ne)
      dimension   tb(nt+1)
      dimension   tw(nt)

      dimension   aw(na)
      dimension   tr(np,ne,nt,nr*nz,na,nm,2)

*-----------------------------------------------------------------------

      dimension tott(np,2) ! frtati 2021/10/05 6 -> np

      character hsunit(14)*29

      data hsunit( 1) / '[1/cm^2/source]              '/
      data hsunit( 2) / '[1/cm^2/MeV/source]          '/
      data hsunit( 3) / '[1/cm^2/Lethargy/source]     '/
      data hsunit( 4) / '[cm/source]                  '/
      data hsunit(11) / '[1/cm^2/nsec/source]         '/
      data hsunit(12) / '[1/cm^2/nsec/MeV/source]     '/
      data hsunit(13) / '[1/cm^2/nsec/Lethargy/source]'/
      data hsunit(14) / '[cm/nsec/source]             '/

      character cha*1
      data cha /"'"/

cfrtati 2021/10/05 chl, chm moved to partmod

      character dc2*4
      character chp(6)*10
      character chq(6)*8


*-----------------------------------------------------------------------

      data pi/3.14159265d+0/

      character aname*3

      character rpa*1
      data rpa /'}'/
      character yen*1

      real(8), allocatable :: vl_r(:),vl_z(:)

      include 'samepage_include/samepage000.inc'

*-----------------------------------------------------------------------
*        set mesh volume ( unit = 4, 14 ; vol = 1.0 )
*-----------------------------------------------------------------------

        vl(ir,iz) = dble( itunt(m)/4  +  itunt(m)/14
     &                  -(itunt(m)/4) * (itunt(m)/10) )
     &            - dble( itunt(m)/4  +  itunt(m)/14
     &                  -(itunt(m)/4) * (itunt(m)/10) - 1 )
     &            * pi * ( rm(ir+1)**2 - rm(ir)**2 )
     &                 * ( zm(iz+1) - zm(iz) )

*-----------------------------------------------------------------------

      icf(ir,iz) = ir + ( iz - 1 ) * nr

*-----------------------------------------------------------------------

      include 'samepage_include/samepage001.inc'

      allocate (vl_r(nz),vl_z(nr))
      vl_r(:) = 0.0d0
      vl_z(:) = 0.0d0
      do iz =1 ,nz
         do ir=1,nr
           vl_r(iz) = vl_r(iz) + vl(ir,iz)
           vl_z(ir) = vl_z(ir) + vl(ir,iz)
         enddo
       enddo


      yen  = char(92)
      igsh = 0

*-----------------------------------------------------------------------
*     if angel is specified
*-----------------------------------------------------------------------
               if( itaty(m) .gt. 0 ) then
                  aname = 'rad'
               else
                  aname = 'deg'
               end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               chp(1) = 'y(msos-1  '
               chp(2) = 'y(msos-2  '
               chp(3) = 'y(msos-3  '
               chp(4) = 'y(msos-4  '
               chp(5) = 'y(msos-5  '
               chp(6) = 'y(msos-6  '
               chq(1) = '1       '
               chq(2) = '2       '
               chq(3) = '3       '
               chq(4) = '4       '
               chq(5) = '5       '
               chq(6) = '6       '

*-----------------------------------------------------------------------
*        c1 : nomalization for source
*-----------------------------------------------------------------------

               c1 = 1.0d+0 / rsouin

*-----------------------------------------------------------------------
*           itunt(m) = 1,4 11,14: /cm^2/source
*                                 /cm^2/sec/source
*                    = 2   12   : /cm^2/MeV/source
*                                 /cm^2/MeV/nsec/source
*                    = 3   13   : /cm^2/Lethargy/source
*                                 /cm^2/Lethargy/nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  4 .or.
     &          itunt(m) .eq. 11 .or. itunt(m) .eq. 14 ) then

               do i = 1, ne

                  ew(i) = 1.d+0

               end do
               ew_sum = 1.0d0

            else if( itunt(m) .eq. 2 .or. itunt(m) .eq. 12 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else if( itunt(m) .eq. 3 .or. itunt(m) .eq. 13 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

            end if

*-----------------------------------------------------------------------

            if( na .le. 1 ) then

                  aw(1) = 1.d+0
                  aw_sum = 1.0d0

            else

               aw_sum = 0.0d0
               do i = 1, na

                  if( itaty(m) .gt. 0 ) then

                     aw(i) = ( ab(i+1) - ab(i) ) / 2.0 / pi

                  else

                     aw(i) = ( ab(i+1) - ab(i) ) / 360.0

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
               tw_sum = tb(nt+1) -tb(1)

            else

               do i = 1, nt

                  tw(i) = 1.d+0

               end do
                tw_sum = 1.0d0

            end if

*-----------------------------------------------------------------------
*        relative error and  unit conversion
*-----------------------------------------------------------------------

                  cmax = 0.0
                  cmin = 1.e+33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1


            facmax(m) = 1.d0
            if( rtfac(m) .lt. 0.d0 ) then
               facmax(m) = 0.d0

               do 101 im = 1, nm
               do 101 ia = 1, na
               do 101 iz = 1, nz
               do 101 ir = 1, nr
               do 101 it = 1, nt
               do 101 ie = 1, ne
               do 101 ip = 1, np

                  if( tr(ip,ie,it,icf(ir,iz),ia,im,1) .gt. 0.d0 ) then

                     fmaxfc = tr(ip,ie,it,icf(ir,iz),ia,im,1)
     &                            / vl(ir,iz) / aw(ia) / ew(ie) / tw(it)

                     if( facmax(m) .lt. fmaxfc) facmax(m) = fmaxfc

                  end if

  101          continue

            end if

            do 100 im = 1, nm
            do 100 ia = 1, na
            do 100 iz = 1, nz
            do 100 ir = 1, nr
            do 100 it = 1, nt
            do 100 ie = 1, ne
            do 100 ip = 1, np

               if( tr(ip,ie,it,icf(ir,iz),ia,im,1) .gt. 0.d0 ) then

!OBINATA(2012.7.11): new method.
                  call calc_stdev(m,Xa,sigx,
     &                            tr(ip,ie,it,icf(ir,iz),ia,im,1),
     &                            tr(ip,ie,it,icf(ir,iz),ia,im,2),
     &           abs(rtfac(m)/facmax(m))/vl(ir,iz)/aw(ia)/ew(ie)/tw(it))

                  tr(ip,ie,it,icf(ir,iz),ia,im,1) = Xa
                  tr(ip,ie,it,icf(ir,iz),ia,im,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,it,icf(ir,iz),ia,im,1) .gt. cmax )
     &                          cmax = tr(ip,ie,it,icf(ir,iz),ia,im,1)

                  if( tr(ip,ie,it,icf(ir,iz),ia,im,1) .lt. cmin )
     &                          cmin = tr(ip,ie,it,icf(ir,iz),ia,im,1)

               else

                  isdz = 1
                  tr(ip,ie,it,icf(ir,iz),ia,im,2) = 0.0

               end if

! sumover
               fact_in = abs(rtfac(m)/facmax(m))
               call ptracrz_sumover_stdev(0,m,
     &                   ip,ie,it,ir,iz,ia,im,
     &                   fact_in,ew(ie),tw(it),vl(ir,iz),aw(ia),
     &                   ew_sum,tw_sum,vl_r(iz),vl_z(ir),aw_sum)

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

!OBINATA(2012.7.2): output *.err
        noe = 1
        if ( itaxs(m,iax) .eq. 10 .and.
     &       ittwo(m) .ne. 4 .and. icntl .ne. 8 ) noe = 2

        do ioe = 1, noe

         if( itall .eq. 2 .and. nobch .lt. maxbch ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.7.2): output *.err
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

               call tadjech(iot,m,iax,1)


      include 'samepage_include/samepage002_petrzam.inc'
      include 'samepage_include/samepagechp_petrzam.inc'
      include 'samepage_include/samepageseti.inc'

*-----------------------------------------------------------------------
*        angle axis
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 12 .or. itaxs(m,iax) .eq. 13 ) then

               inum = 0

            do imi = 1, nm
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

               if(iloopmode .ge. 1) then
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

                  write(iot,'(
     &            "#  ie =",i3,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ie, eb(ie), eb(ie+nestepi)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3,/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  if( itaty(m) .gt. 0.0 ) then

                     write(iot,'(/"x: Angle [radian]")')

                  else

                     write(iot,'(/"x: Angle [degree]")')

                  end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Flux ",a29)') hsunit(itunt(m))

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

       include 'samepage_include/flux_petrzam_a.inc'

*-----------------------------------------------------------------------
             if(iloopmode .eq. 0) then
               if( itmlp(m) .eq. 0 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  (ir,iz) = (",
     &                     i3,",",i3,"),  ie =",i3,a1)')
     &                     cha, inum, ir, iz, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  (ir,iz) = (",
     &                     i3,",",i3,"),  ie =",i3,
     &                     ",  it =",i3,a1)')
     &                     cha, inum, ir, iz, ie, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  (ir,iz) = (",i3,",",i3,")",
     &                     ",  mset =",i3,",  ie =",i3,a1)')
     &                     cha, inum, ir, iz, itmnt(m,im), ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  (ir,iz) = (",i3,",",i3,")",
     &                     ",  mset =",i3,",  ie =",i3,
     &                     ",  it =",i3,a1)')
     &                     cha, inum, ir, iz, itmnt(m,im), ie, it, cha
                  end if
               end if

             else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                  changelsub(3) = " "
                else
                  write(changelsub(3),'(",  mset =",i3)') itmnt(m,im)
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


               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                     volsum =0.0d0
                     do iii = iai, iai+nastepi-1
                     do ii  = izi, izi+nzstepi-1
                     do i   = iri, iri+nrstepi-1
                       volsum = volsum + vl(i,ii) * aw(iii)
                     end do
                     end do
                     end do

                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, volsum,
     &                  rm(ir), rm(ir+nrstepi), zm(iz), zm(iz+nzstepi)
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
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

*-----------------------------------------------------------------------
*        energy axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 1 ) then

               inum = 0

            do imi = 1, nm
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

               if(iloopmode .ge. 1) then
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
               if( na .gt. 1 ) then
                  write(iot,'(
     &            "#  ia ="i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3,/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstep)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Energy [MeV]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Flux ",a29)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if
               if( itunt(m) .eq. 3 .or. itunt(m) .eq. 13 .or.
     &             itety(m) .eq. 3 .or. itety(m) .eq. 5 ) then

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

       include 'samepage_include/flux_petrzam_erg.inc'

*-----------------------------------------------------------------------

              if(iloopmode .eq. 0) then
               if( itmlp(m) .eq. 0 ) then
                  if( ittty(m) .eq. 0 ) then
                     if( na .le. 1 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  (ir,iz) = (",
     &                        i3,",",i3,")",a1)')
     &                        cha, inum, ir, iz, cha
                     else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  (ir,iz) = (",
     &                        i3,",",i3,")",",  ia =",i3,a1)')
     &                        cha, inum, ir, iz, ia, cha
                     end if
                  else
                     if( na .le. 1 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  (ir,iz) = (",
     &                        i3,",",i3,"),  it =",i3,a1)')
     &                        cha, inum, ir, iz, it, cha
                     else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  (ir,iz) = (",
     &                        i3,",",i3,"),  it =",i3,",  ia =",
     &                        i3,a1)')
     &                        cha, inum, ir, iz, it, ia, cha
                     end if
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     if( na .le. 1 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  (ir,iz) = (",i3,",",i3,")",
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, ir, iz, itmnt(m,im), cha
                     else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  (ir,iz) = (",i3,",",i3,")",
     &                        ",  ia =",i3,",  mset =",i3,a1)')
     &                        cha, inum, ir, iz, ia, itmnt(m,im), cha
                     end if
                  else
                     if( na .le. 1 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  (ir,iz) = (",i3,",",i3,")",
     &                        ",  mset =",i3,",  it =",i3,a1)')
     &                        cha, inum, ir, iz, itmnt(m,im), it, cha
                     else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  (ir,iz) = (",i3,",",i3,")",
     &                        ",  ia =", i3,
     &                        ",  mset =",i3,",  it =",i3,a1)')
     &                        cha, inum, ir, iz, ia,
     &                        itmnt(m,im), it, cha
                     end if
                  end if
               end if
             else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                  changelsub(3) = " "
                else
                  write(changelsub(3),'(",  mset =",i3)') itmnt(m,im)
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


               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                     volsum =0.0d0
                     do iii = iai, iai+nastepi-1
                     do ii  = izi, izi+nzstepi-1
                     do i   = iri, iri+nrstepi-1
                       volsum = volsum + vl(i,ii) * aw(iii)
                     end do
                     end do
                     end do

                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, volsum,
     &                  rm(ir), rm(ir+nrstepi), zm(iz), zm(iz+nzstepi)
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if
               if( na .gt. 1 ) then
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

*-----------------------------------------------------------------------
*        r axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 6 ) then

               inum = 0

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
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "iz  =",i3,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie, iz,
     &                        eb(ie), eb(ie+nestepi),
     &                        zm(iz), zm(iz+nzstepi)
               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)
               if( na .gt. 1 ) then
                  write(iot,'(
     &            "#  ia ="i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3,/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: r [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Flux ",a29)') hsunit(itunt(m))

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

       include 'samepage_include/flux_petrzam_r.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( itmlp(m) .eq. 0 ) then
                  if( ittty(m) .eq. 0 ) then
                     if( na .le. 1 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,
     &                        ",  iz =",i3,a1)')
     &                        cha, inum, ie, iz, cha
                     else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,
     &                        ",  iz =",i3,",  ia =",i3,a1)')
     &                        cha, inum, ie, iz, ia, cha
                     end if
                  else
                     if( na .le. 1 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,
     &                        ",  iz =",i3,",  it =",i3,a1)')
     &                        cha, inum, ie, iz, it, cha
                     else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,
     &                        ",  iz =",i3,",  ia =",i3,
     &                        ",  it =",i3,a1)')
     &                        cha, inum, ie, iz, ia, it, cha
                     end if
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     if( na .le. 1 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,
     &                        ",  iz =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, ie, iz, itmnt(m,im), cha
                     else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,
     &                        ",  iz =",i3,",  ia =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, ie, iz, ia, itmnt(m,im), cha
                     end if
                  else
                     if( na .le. 1 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,
     &                        ",  iz =",i3,
     &                        ",  mset =",i3,",  it =",i3,a1)')
     &                        cha, inum, ie, iz, itmnt(m,im), it, cha
                     else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,
     &                        ",  iz =",i3,",  ia =",i3,
     &                        ",  mset =",i3,",  it =",i3,a1)')
     &                        cha, inum, ie, iz, ia,
     &                        itmnt(m,im), it, cha
                     end if
                  end if
               end if

             else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                  changelsub(3) = " "
                else
                  write(changelsub(3),'(",  mset =",i3)') itmnt(m,im)
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
                write(changelsub(5),'(",  iz =",i3)') iz
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


               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                  eb(ie), eb(ie+nestepi), zm(iz), zm(iz+nzstepi)
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)

               end if
               if( na .gt. 1 ) then
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

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 5 ) then

               inum = 0

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
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "ir  =",i3,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie, ir,
     &                        eb(ie), eb(ie+nestepi),
     &                        rm(ir), rm(ir+nrstepi)
               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)
               if( na .gt. 1 ) then
                  write(iot,'(
     &            "#  ia ="i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+nastepi)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3,/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Flux ",a29)') hsunit(itunt(m))

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

       include 'samepage_include/flux_petrzam_z.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( itmlp(m) .eq. 0 ) then
                  if( ittty(m) .eq. 0 ) then
                     if( na .le. 1 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,
     &                        ",  ir =",i3,a1)')
     &                        cha, inum, ie, ir, cha
                     else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,
     &                        ",  ir =",i3,",  ia =",i3,a1)')
     &                        cha, inum, ie, ir, ia, cha
                     end if
                  else
                     if( na .le. 1 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,
     &                        ",  ir =",i3,",  it =",i3,a1)')
     &                        cha, inum, ie, ir, it, cha
                     else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,
     &                        ",  ir =",i3,",  ia =",i3,
     &                        ",  it =",i3,a1)')
     &                        cha, inum, ie, ir, ia, it, cha
                     end if
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     if( na .le. 1 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,
     &                        ",  ir =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, ie, ir, itmnt(m,im), cha
                     else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,
     &                        ",  ir =",i3,",  ia =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, ie, ir, ia, itmnt(m,im), cha
                     end if
                  else
                     if( na .le. 1 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,
     &                        ",  ir =",i3,
     &                        ",  mset =",i3,",  it =",i3,a1)')
     &                        cha, inum, ie, ir, itmnt(m,im), it, cha
                     else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,
     &                        ",  ir =",i3,",  ia =",i3,
     &                        ",  mset =",i3,",  it =",i3,a1)')
     &                        cha, inum, ie, ir, ia,
     &                        itmnt(m,im), it, cha
                     end if
                  end if
               end if
             else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                  changelsub(3) = " "
                else
                  write(changelsub(3),'(",  mset =",i3)') itmnt(m,im)
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
                write(changelsub(5),'(",  ir =",i3)') ir
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


               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                  eb(ie), eb(ie+nestepi), rm(ir), rm(ir+nrstepi)
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if
               if( na .gt. 1 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi1)
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

*-----------------------------------------------------------------------
*        time axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 11 ) then

               inum = 0

            do imi = 1, nm
            do iri = 1, nr, nrstepi
            do izi = 1, nz, nzstepi
            do iai = 1, na, nastepi
            do iei = 1, ne, nestepi
            do ipi = 1, np, npstepi
               im = imi
               ir = iri
               iz = izi
               ia = iai
               ie = iei
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

               if(iloopmode .ge. 1) then
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
               if( na .gt. 1 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "#   a = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, ab(ia), ab(ia+nastepi)
               end if
                  write(iot,'(
     &            "#  ie =",i3,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ie, eb(ie), eb(ie+nestepi)

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Time [nsec]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Flux ",a29)') hsunit(itunt(m))

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

       include 'samepage_include/flux_petrzam_time.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( itmlp(m) .eq. 0 ) then
                  if( na .le. 1 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  (ir,iz) = (",
     &                     i3,",",i3,"),  ie =",i3,a1)')
     &                     cha, inum, ir, iz, ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  (ir,iz) = (",
     &                     i3,",",i3,"),  ia =",i3,
     &                     ",  ie =",i3,a1)')
     &                     cha, inum, ir, iz, ia, ie, cha
                  end if
               else
                  if( na .le. 1 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  (ir,iz) = (",i3,",",i3,")",
     &                     ",  mset =",i3,",  ie =",i3,a1)')
     &                     cha, inum, ir, iz, itmnt(m,im), ie, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  (ir,iz) = (",i3,",",i3,")",
     &                     ",  ia =",i3,
     &                     ",  mset =",i3,",  ie =",i3,a1)')
     &                     cha, inum, ir, iz, ia, itmnt(m,im), ie, cha
                  end if
               end if

             else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                  changelsub(3) = " "
                else
                  write(changelsub(3),'(",  mset =",i3)') itmnt(m,im)
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


               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                     volsum =0.0d0
                     do iii = iai, iai+nastepi-1
                     do ii  = izi, izi+nzstepi-1
                     do i   = iri, iri+nrstepi-1
                       volsum = volsum + vl(i,ii) * aw(iii)
                     end do
                     end do
                     end do

                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, volsum,
     &                  rm(ir), rm(ir+nrstepi), zm(iz), zm(iz+nzstepi)
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if
               if( na .gt. 1 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        rz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 10 ) then

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
            do ip = 1, np
            do ie = 1, ne
            do ia = 1, na
            do it = 1, nt

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

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "msos = ",a8,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie, chq(ip),
     &                        eb(ie), eb(ie+1)
               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)
               if( na .gt. 1 ) then
                  write(iot,'(
     &            "#  ia ="i3/
     &            "# ",a3," = (",
     &                         1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        ia, aname, ab(ia), ab(ia+1)
               end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3,/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+1)
               end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               if( itmlp(m) .eq. 0 ) then
                  if( ittty(m) .eq. 0 ) then
                     if( na .le. 1 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,a1)')
     &                        cha, inum, ie, cha
                     else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ia =",i3,",  ie =",i3,a1)')
     &                        cha, inum, ia, ie, cha
                     end if
                  else
                     if( na .le. 1 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,",  it =",i3,a1)')
     &                        cha, inum, ie, it, cha
                     else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ia =",i3,
     &                        ",  ie =",i3,",  it =",i3,a1)')
     &                        cha, inum, ia, ie, it, cha
                     end if
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     if( na .le. 1 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, ie, itmnt(m,im), cha
                     else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ia =",i3,
     &                        ",  ie =",i3,
     &                        ",  mset =",i3,a1)')
     &                        cha, inum, ia, ie, itmnt(m,im), cha
                     end if
                  else
                     if( na .le. 1 ) then
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ie =",i3,
     &                        ",  mset =",i3,",  it =",i3,a1)')
     &                        cha, inum, ie, itmnt(m,im), it, cha
                     else
                        write(iot,'(/a1,"no. =",i3,
     &                        ",  ia =",i3,
     &                        ",  ie =",i3,
     &                        ",  mset =",i3,",  it =",i3,a1)')
     &                        cha, inum, ia, ie, itmnt(m,im), it, cha
                     end if
                  end if
               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

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

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

                  write(iot,'(/a4," y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7," ; x = ",
     &            1p1g14.7," to ",1p1g14.7," by ",
     &            1p1g14.7," ;")') dc2,
     &            rm(nr) + rtrdl(m)/2.0, rm(1) + rtrdl(m)/2.0, rtrdl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

!OBINATA(2012.7.2): output *.err
               write(iot,'(1p10e11.3)')
     &        ( ( tr(ip,ie,it,icf(ir,iz),ia,im,ioe),
     &            iz = 1, nz ), ir = nr, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# r          z        ",
     &                      "  flux       r.err")')

               do iz = 1, nz
               do ir = 1, nr

                  write(iot,'(1p3e11.3,0pf8.4)')
     &               rm(ir)  + rtrdl(m)/2.0,
     &               zm(iz)  + rtzdl(m)/2.0,
     &            tr(ip,ie,it,icf(ir,iz),ia,im,1),
     &            tr(ip,ie,it,icf(ir,iz),ia,im,2)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

                  write(iot,'("#   r = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7/
     &                        "#   z = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            rm(1) + rtrdl(m)/2.0, rm(nr) + rtrdl(m)/2.0, rtrdl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'r/z',( zm(iz) + rtzdl(m)/2.0, iz = 1, nz )

               do ir = nr, 1, -1

!OBINATA(2012.7.2): output *.err
                  write(iot,'(1p1000e11.3)')
     &            rm(ir) + rtrdl(m)/2.0,
     &            ( tr(ip,ie,it,icf(ir,iz),ia,im,ioe), iz = 1, nz )

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

      if ( ioe .eq. 2 ) then

          write(iot,'("y: Relative Error")')

      else if( itazl(m) .eq. 0 ) then

          write(iot,'("y: Flux ",a29)') hsunit(itunt(m))

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
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  msos &=& ",a8)')
     &                     yen, eb(ie), eb(ie+1), chq(ip)
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if
               if( na .gt. 1 ) then
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

            call prestart(m,iot) !OBINATA(2012.6.13)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

      end do
      end do

*-----------------------------------------------------------------------
      include 'samepage_include/samepage999.inc'

! sumover
      deallocate (vl_r,vl_z)

      return
      end


************************************************************************
*                                                                      *
      subroutine tadjntxyz(ncol,m,nl,lt,np,nx,ny,nz,ne,
     &                    nm,nt,xm,ym,zm,eb,tb,tr,trEVENT,
     &                    itrmax,itrmin)
*                                                                      *
*       adjoint tally in xyz scoring mesh                              *
*       last modified by K.Niita on 2016/01/22                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*              9 : termination by time cut-off                         *
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
      include 'err.inc'

      parameter( rlit = 29.97925d0 )

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)

*-----------------------------------------------------------------------
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)

      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall66/ itadm(itlmax), rtade(itlmax), rtadw(itlmax)

      common /tall82/ itcnth(9,itlmax)

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)

      common /enginit/ engini
!$OMP THREADPRIVATE(/enginit/)

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension xm(nx+1)
      dimension ym(ny+1)
      dimension zm(nz+1)
      dimension eb(ne+1)
      dimension tb(nt+1)
      dimension tr(np,ne,nt,nx*ny*nz,nm,2)
      dimension trEVENT(np,ne,nt,nx*ny*nz,nm)   !OBINATA(2012.6.13): as Ct
      real(8),allocatable,save:: tr0(:,:,:,:,:) !OBINATA(2012.6.13): as C
      dimension itrmax(7),itrmin(7)

*-----------------------------------------------------------------------

      dimension ud(3)

      data dmax  /1.0d+19/
      data dmax0 /1.0d+18/

      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt
      dimension facm(6)

*-----------------------------------------------------------------------

      common /paraj/ mstz(300), parz(300)

*-----------------------------------------------------------------------

      common /stat / istdev, irestart, ireschk
      common /cparm/ maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      common /egs5cmn10/denstepold,denstepnew,deinit
      real*8            denstepold,denstepnew,deinit
!$OMP THREADPRIVATE(/egs5cmn10/)

*-----------------------------------------------------------------------

      common /spred/ nspred, nwsprd, nedisp, itstep, ndedx

      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

*-----------------------------------------------------------------------

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny

*-----------------------------------------------------------------------

         small = parz(28) * 10.d0
         rng1 = 0.0d0 !FURUTA

*-----------------------------------------------------------------------
cKN 2016/01/22  initial multi-source number
*               energy bin
*-----------------------------------------------------------------------

         jsos = nsos(ibksos+1,ipomp+1)

         eaini = engini
         eamin = rtade(m)
         eamax = rtadw(m)

*-----------------------------------------------------------------------

      if (istdev .eq. 2) then

        call readitrminmax7(itrmin,itrmax,(/ np,ne,nt,nx,ny,nz,nm/),
     &                      mnp,mne,mnt,mnx,mny,mnz,mnm,
     &                      mxp,mxe,mxt,mxx,mxy,mxz,mxm)

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
*OBINATA(2012.6.13)
*-----------------------------------------------------------------------

         if (( ncol .eq. 0 .or. ncol .eq. 4 )
     &                              .and. istdev .eq. 2) then
           if ((nocas.gt.1.or.ncol.eq.0) .and. ihistcount.ne.1 ) then

             do im = mnm,mxm
             do iz = mnz,mxz
             do iy = mny,mxy
             do ix = mnx,mxx
               ixyz = icf(ix,iy,iz)
               tr(:,mne:mxe,:,ixyz,im,1)
     &           = tr(:,mne:mxe,:,ixyz,im,1)
     &           + trEVENT(:,mne:mxe,:,ixyz,im)
               tr(:,mne:mxe,:,ixyz,im,2)
     &           = tr(:,mne:mxe,:,ixyz,im,2)
     &           + trEVENT(:,mne:mxe,:,ixyz,im) ** 2
             enddo
             enddo
             enddo
             enddo

! sumover
              call ttracxyz_sumover(m,1,
     &                   np, ne, nt, nx, ny, nz, nm, trEVENT)

           end if

           do im = mnm,mxm
           do iz = mnz,mxz
           do iy = mny,mxy
           do ix = mnx,mxx
             trEVENT(:,mne:mxe,:,icf(ix,iy,iz),im) = 0
           enddo
           enddo
           enddo
           enddo

           call resetitrminmax(itrmin,itrmax,7,(/np,ne,nt,nx,ny,nz,nm/))

         end if

*-----------------------------------------------------------------------
*        end of batch ( in case of istdev = 1 )
*OBINATA(2012.6.13)
*-----------------------------------------------------------------------

         if ( ncol .eq. 0 .and. istdev .eq. 1) then
!$OMP MASTER
             allocate( tr0(np,ne,nt,nx*ny*nz,nm) )
             tr0(:,:,:,:,:) = 0.d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP CRITICAL (tadjntxyz_crit)
             tr0(:,:,:,:,:) = tr0(:,:,:,:,:) + trEVENT(:,:,:,:,:)
!$OMP END CRITICAL (tadjntxyz_crit)
!$OMP BARRIER
!$OMP MASTER
             tr(:,:,:,:,:,1) = tr(:,:,:,:,:,1)
     &                       + tr0(:,:,:,:,:) / maxcas
             tr(:,:,:,:,:,2) = tr(:,:,:,:,:,2)
     &                     + ( tr0(:,:,:,:,:) / maxcas ) ** 2

! sumover
              call ttracxyz_sumover(m,maxcas,
     &                   np, ne, nt, nx, ny, nz, nm, tr0)

             deallocate( tr0 )
!$OMP END MASTER

           trEVENT(:,:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------
*        check of ncol
*-----------------------------------------------------------------------

         if( ncol .lt. 9 ) return

*-----------------------------------------------------------------------
*        check of energy
*-----------------------------------------------------------------------

            if( e(ibke+no,ipomp+1)   .lt. eamin ) goto 999
            if( ec(ibkec+no,ipomp+1) .ge. eamax ) goto 999

            if( eaini .lt. eb(1) ) goto 999
            if( eaini .ge. eb(ne+1) ) goto 999

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
*        check of mat
*-----------------------------------------------------------------------

            if( nl .gt. 0 ) then

                  do i = 1, nl

                     if( itmcn(m) .gt. 0 .and.
     &                   idmn(mat) .eq. lt(i) ) goto 502
                     if( itmcn(m) .lt. 0 .and.
     &                   idmn(mat) .eq. lt(i) ) return

                  end do

                     if( itmcn(m) .gt. 0 ) return

            end if

  502          continue

*-----------------------------------------------------------------------
*        check of particles
*-----------------------------------------------------------------------

               if( jsos .gt. np ) return

*-----------------------------------------------------------------------
*           check of time
*-----------------------------------------------------------------------

               tparti = abs(t(ibkt+no,ipomp+1))
               tpartf = abs(tc(ibktc+no,ipomp+1))

               if( tparti .ge. tb(nt+1) ) goto 999
               if( tpartf .lt. tb(1) ) goto 999

*-----------------------------------------------------------------------

            icli = idgr(iblz1)

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

            if( xxa-small .lt. xm(1) .and.
     &          xxc-small .lt. xm(1) ) goto 999

            if( yya-small .lt. ym(1) .and.
     &          yyc-small .lt. ym(1) ) goto 999

            if( zza-small .lt. zm(1) .and.
     &          zzc-small .lt. zm(1) ) goto 999

            if( xxa+small .ge. xm(nx+1) .and.
     &          xxc+small .ge. xm(nx+1) ) goto 999

            if( yya+small .ge. ym(ny+1) .and.
     &          yyc+small .ge. ym(ny+1) ) goto 999

            if( zza+small .ge. zm(nz+1) .and.
     &          zzc+small .ge. zm(nz+1) ) goto 999

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
            tpt = tparti

            se  = eaini
            xpp = xxa
            ypp = yya
            zpp = zza

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
                  if( xm(i) .gt. xpp + small ) then
                         ixm = i - 1
                         ixc = i - 1
                         goto 32
                  end if
               end do

            else

               ixm = nx + 2
               ixc = nx + 1

               do i = 1, nx + 1
                  if( xm(i) .ge. xpp - small ) then
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
                  if( ym(i) .gt. ypp + small ) then
                         iym = i - 1
                         iyc = i - 1
                         goto 35
                  end if
               end do

            else

               iym = ny + 2
               iyc = ny + 1

               do i = 1, ny + 1
                  if( ym(i) .ge. ypp - small ) then
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
                  if( zm(i) .gt. zpp + small ) then
                         izm = i - 1
                         izc = i - 1
                         goto 38
                  end if
               end do

            else

               izm = nz + 2
               izc = nz + 1

               do i = 1, nz + 1
                  if( zm(i) .ge. zpp - small ) then
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

      iloop=0 ! T.Sato 2021/11/06 avoid infinite loop
      do 50 while ( .not. tot .ge. dis )
       iloop=iloop+1
       if(iloop.gt.1000000) then
        write(ErrCha,'("warning: iloop is too large")')
        ErrID = 'L:6062/R:tadjntxyz/F:talls09.f'
        call ErrWrite(ErrID,ErrCha)
        exit
       endif

*-----------------------------------------------------------------------
*     calculation is finished
*-----------------------------------------------------------------------


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

               tt = dd

            if( dd .gt. dmax0 ) return

*-----------------------------------------------------------------------
*        propagate position upto the boundary or the final point
*-----------------------------------------------------------------------

            if( tot + dd .ge. dis - small ) then

               dd  = dis - tot
               tot = dis
               xpp = xxc
               ypp = yyc
               zpp = zzc

            else

               tot = tot + dd
               xpp = xpp + tt * ud(1)
               ypp = ypp + tt * ud(2)
               zpp = zpp + tt * ud(3)

            end if

*-----------------------------------------------------------------------

               ee = se

*-----------------------------------------------------------------------
*        time evolution
*-----------------------------------------------------------------------

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

                  tpc = tpt
                  tpt = tpt + timd

*-----------------------------------------------------------------------
*        booking
*-----------------------------------------------------------------------

         if( ixc .ge. 1 .and. ixc .lt. nx + 1 .and.
     &       iyc .ge. 1 .and. iyc .lt. ny + 1 .and.
     &       izc .ge. 1 .and. izc .lt. nz + 1 ) then

*-----------------------------------------------------------------------
*           initial and final energy cell
*-----------------------------------------------------------------------

            do i = 1, ne

               if( se .ge. eb(i) .and.
     &             se .lt. eb(i+1) ) goto 30

            end do

   30          ie1 = min( i, ne )

            do i = ne, 1, -1

               if( ee .ge. eb(i) .and.
     &             ee .lt. eb(i+1) ) goto 40

            end do

   40          ie2 = max( i, 1 )

*-----------------------------------------------------------------------
*           rrl : track length and rrr : charge particles correction
*-----------------------------------------------------------------------

                  rrl = dd
                  rrr = dd
                  rrt = 0.0d0

                  tpi = tpc

*-----------------------------------------------------------------------
*        tally
*-----------------------------------------------------------------------

            do 270 ie = ie1, ie2, -1

               emie = max( ee, eb(ie) )
               emiu = min( se, eb(ie+1) )
               erg  = ( emiu + emie ) / 2.0

               call fmfac(m,icli,erg,oldwt,facm)

*-----------------------------------------------------------------------
*              initial and final time
*-----------------------------------------------------------------------

               dst = rrl
               tmd = 0.0d0
               CC  = 0.0d0

               if( erg .gt. 0.1 .or. rtyp .eq. 0.0d0 ) then

                  tmd = dst * ( erg + rtyp )
     &                / sqrt( erg * ( erg + 2.0 * rtyp ) )
     &                / rlit
                  CC  = 1.d0 / ( erg + rtyp )
     &                * sqrt( erg * ( erg + 2.0 * rtyp ) )
     &                * rlit

               else if( erg .gt. 0.0 ) then

                  tmd = dst * sqrt( rtyp / 2.0 / erg ) / rlit
                  CC = 1.d0 / sqrt( rtyp / 2.0 / erg ) * rlit

               end if

               tpf = tpi + tmd

               do 271 it = 1, nt

                  if( tpf .lt. tb(it) .or. tpi .ge. tb(it+1) ) goto 271

                  tin = max( tpi, tb(it) )
                  tfn = min( tpf, tb(it+1) )
                  tmd = max( 0.0d0, tfn - tin )

                  dst = tmd * CC

*-----------------------------------------------------------------------

                  tlv = oldwt * dst
                  ip = jsos

                  tlv = tlv * (eb(ne+1)-eb(1)) / (eamax-eamin) ! T.Sato 2020/09/06

               do im = 1, nm

                  trEVENT(ip,ie,it,icf(ixc,iyc,izc),im) =
     &            trEVENT(ip,ie,it,icf(ixc,iyc,izc),im)
     &                      + tlv * facm(im)

               end do

*-----------------------------------------------------------------------

  271          continue

               rrr = rrr - rrl
               tpi = tpf

  270       continue

         if (istdev .eq. 2) then

cKN 2015/12/01 ??
           call setitrmin(itrmin,2,7,(/ie2, 1,ixc,iyc,izc, 1/))
           call setitrmax(itrmax,2,7,(/ie1,it,ixc,iyc,izc,nm/))

         endif

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*     next position
*-----------------------------------------------------------------------

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

  50  continue

*-----------------------------------------------------------------------

  999 continue

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine padjntxyz(m,np,nl,lt,
     &                    nx,ny,nz,ne,nm,nt,xm,ym,zm,eb,tb,tr,
     &                    igsh,idasa)
*                                                                      *
*       output xyz scoring mesh adjoint tally                          *
*       last modified by K.Niita on 2016/01/22                         *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
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

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall21/ rtfac(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)

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

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)
      dimension   eb(ne+1)
      dimension   ew(ne)
      dimension   tb(nt+1)
      dimension   tw(nt)
      dimension   tr(np,ne,nt,nx*ny*nz,nm,2)

      integer,allocatable :: ixyz(:)

*-----------------------------------------------------------------------

      dimension tott(np,2) ! frtati 2021/10/05 6 -> np

      character hsunit(14)*29

      data hsunit( 1) / '[1/cm^2/source]              '/
      data hsunit( 2) / '[1/cm^2/MeV/source]          '/
      data hsunit( 3) / '[1/cm^2/Lethargy/source]     '/
      data hsunit( 4) / '[cm/source]                  '/
      data hsunit(11) / '[1/cm^2/nsec/source]         '/
      data hsunit(12) / '[1/cm^2/nsec/MeV/source]     '/
      data hsunit(13) / '[1/cm^2/nsec/Lethargy/source]'/
      data hsunit(14) / '[cm/nsec/source]             '/

      character cha*1
      data cha /"'"/

cfrtati 2021/10/05 chl, chm moved to partmod

      character dc2*4

      character chp(6)*10
      character chq(6)*8


*-----------------------------------------------------------------------

      character rpa*1
      data rpa /'}'/
      character yen*1

      real(8),allocatable :: vl_x(:,:),vl_y(:,:),vl_z(:,:)

      include 'samepage_include/samepage000.inc'

*-----------------------------------------------------------------------

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny

*-----------------------------------------------------------------------
*        set mesh volume ( unit = 4, 14 ; vol = 1.0 )
*-----------------------------------------------------------------------

         vl(ix,iy,iz) = dble( itunt(m)/4  +  itunt(m)/14
     &                      -(itunt(m)/4) * (itunt(m)/10) )
     &                - dble( itunt(m)/4  +  itunt(m)/14
     &                      -(itunt(m)/4) * (itunt(m)/10) - 1 )
     &                * vls(nl,lt,itmcn(m),itvm(m),itmtr(m,4),
     &                      xm(ix),xm(ix+1),
     &                      ym(iy),ym(iy+1),
     &                      zm(iz),zm(iz+1))

*-----------------------------------------------------------------------
      include 'samepage_include/samepage001.inc'

! sumover
      allocate (vl_x(ny,nz),vl_y(nx,nz),vl_z(nx,ny))
      vl_x(:,:) = 0.0d0
      vl_y(:,:) = 0.0d0
      vl_z(:,:) = 0.0d0
      do iz = 1,nz
         do iy = 1,ny
           do ix = 1,nx
             vl_x(iy,iz) = vl_x(iy,iz) + vl(ix,iy,iz)
             vl_y(ix,iz) = vl_y(ix,iz) + vl(ix,iy,iz)
             vl_z(ix,iy) = vl_z(ix,iy) + vl(ix,iy,iz)
           enddo
          enddo
      enddo

      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

*-----------------------------------------------------------------------
*     for igshow
*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

            npg = np
            neg = ne
            ntg = nt

         else

            npg = 1
            neg = 1
            ntg = 1

         end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               chp(1) = 'y(msos-1  '
               chp(2) = 'y(msos-2  '
               chp(3) = 'y(msos-3  '
               chp(4) = 'y(msos-4  '
               chp(5) = 'y(msos-5  '
               chp(6) = 'y(msos-6  '
               chq(1) = '1       '
               chq(2) = '2       '
               chq(3) = '3       '
               chq(4) = '4       '
               chq(5) = '5       '
               chq(6) = '6       '

*-----------------------------------------------------------------------
*           itunt(m) = 1,4 11,14: /cm^2/source
*                                 /cm^2/sec/source
*                    = 2   12   : /cm^2/MeV/source
*                                 /cm^2/MeV/nsec/source
*                    = 3   13   : /cm^2/Lethargy/source
*                                 /cm^2/Lethargy/nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  4 .or.
     &          itunt(m) .eq. 11 .or. itunt(m) .eq. 14 ) then

               do i = 1, ne

                  ew(i) = 1.d+0

               end do
               ew_sum = 1.0d0

            else if( itunt(m) .eq. 2 .or. itunt(m) .eq. 12 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else if( itunt(m) .eq. 3 .or. itunt(m) .eq. 13 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

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

                  cmax = 0.0
                  cmin = 1.e+33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

         if( igsh .eq. 0 ) then

               c1 = 1.0d+0 / rsouin

            facmax(m) = 1.d0
            if( rtfac(m) .lt. 0.d0 ) then
               facmax(m) = 0.d0

               do 101 im = 1, nm
               do 101 iz = 1, nz
               do 101 iy = 1, ny
               do 101 ix = 1, nx
               do 101 it = 1, nt
               do 101 ie = 1, ne
               do 101 ip = 1, np

                  if( tr(ip,ie,it,icf(ix,iy,iz),im,1) .gt. 0.d0 ) then

                     fmaxfc = tr(ip,ie,it,icf(ix,iy,iz),im,1)
     &                                 / vl(ix,iy,iz) / ew(ie) / tw(it)

                     if( facmax(m) .lt. fmaxfc) facmax(m) = fmaxfc

                  end if

  101          continue

            end if

            do 100 im = 1, nm
            do 100 iz = 1, nz
            do 100 iy = 1, ny
            do 100 ix = 1, nx
            do 100 it = 1, nt
            do 100 ie = 1, ne
            do 100 ip = 1, np

               if( tr(ip,ie,it,icf(ix,iy,iz),im,1) .gt. 0.d0 ) then
!OBINATA(2012.7.11): new method.
                  call calc_stdev(m,Xa,sigx,
     &                            tr(ip,ie,it,icf(ix,iy,iz),im,1),
     &                            tr(ip,ie,it,icf(ix,iy,iz),im,2),
     &             abs(rtfac(m)/facmax(m))/vl(ix,iy,iz)/ew(ie)/tw(it))

                  tr(ip,ie,it,icf(ix,iy,iz),im,1) = Xa
                  tr(ip,ie,it,icf(ix,iy,iz),im,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,it,icf(ix,iy,iz),im,1) .gt. cmax )
     &                     cmax = tr(ip,ie,it,icf(ix,iy,iz),im,1)

                  if( tr(ip,ie,it,icf(ix,iy,iz),im,1) .lt. cmin )
     &                     cmin = tr(ip,ie,it,icf(ix,iy,iz),im,1)

               else

                  isdz = 1
                  tr(ip,ie,it,icf(ix,iy,iz),im,2) = 0.d+0

               end if
! sumover
               fact_in = abs(rtfac(m)/facmax(m))
               call ptracxyz_sumover_stdev(0,m,ip,ie,it,ix,iy,iz,im,
     &              fact_in,ew(ie),tw(it),vl(ix,iy,iz),
     &              ew_sum,tw_sum,vl_x(iy,iz),vl_y(ix,iz),vl_z(ix,iy))


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

!OBINATA(2012.7.2): output *.err
        noe = 1
        if ( any( itaxs(m,iax) .eq. (/ 7, 8, 9 /) )
     &      .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &    noe = 2

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

               call tadjech(iot,m,iax,1)

      include 'samepage_include/samepage002_petxyzm.inc'
      include 'samepage_include/samepagechp_petxyzm.inc'
      include 'samepage_include/samepageseti.inc'

*-----------------------------------------------------------------------
*        energy axis
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 1 ) then

               inum = 0

            do imi = 1, nm
            do ixi = 1, nx, nxstepi
            do iyi = 1, ny, nystepi
            do izi = 1, nz, nzstepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ix = ixi
               iy = iyi
               iz = izi
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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ix  =",i3,3x,
     &            "iy  =",i3,3x,
     &            "iz  =",i3,/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ix, iy, iz,
     &                        xm(ix), xm(ix+nxstepi),
     &                        ym(iy), ym(iy+nystepi),
     &                        zm(iz), zm(iz+nzstepi)
               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3,/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Energy [MeV]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Flux ",a29)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itunt(m) .eq. 3 .or. itunt(m) .eq. 13 .or.
     &             itety(m) .eq. 3 .or. itety(m) .eq. 5 ) then

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

       include 'samepage_include/flux_petxyzm_erg.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( itmlp(m) .eq. 0 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &               ", (ix,iy,iz) = (",i3,",",
     &               i3,",",i3,")",a1)')
     &               cha, inum, ix, iy, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &               ", (ix,iy,iz) = (",i3,",",
     &               i3,",",i3,"),  it =",i3,a1)')
     &               cha, inum, ix, iy, iz, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &               ", (ix,iy,iz) = (",i3,
     &               ",",i3,",",i3,")",
     &               ",  mset =",i3,a1)')
     &               cha, inum, ix, iy, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &               ", (ix,iy,iz) = (",i3,
     &               ",",i3,",",i3,")",
     &               ",  mset =",i3,",  it =",i3,a1)')
     &               cha, inum, ix, iy, iz, itmnt(m,im), it, cha
                  end if
               end if

             else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                  changelsub(3) = " "
                else
                  write(changelsub(3),'(",  mset =",i3)') itmnt(m,im)
                end if

               write(changelsub(4),'(", (ix,iy,iz) = (",i3,",",
     &               i3,",",i3,")")') ix, iy, iz
                changelsub(5) = " "
                changelsub(6) = " "
                if(iloopmode .eq. 3) then
                  changelsub(4) = " "
                  write(changelsub(5),'(",  iy =",i3)') iy
                  write(changelsub(6),'(",  iz =",i3)') iz
                end if
                if(iloopmode .eq. 4) then
                  write(changelsub(4),'(",  ix =",i3)') ix
                  changelsub(5) = " "
                  write(changelsub(6),'(",  iz =",i3)') iz
                end if
                if(iloopmode .eq. 5) then
                  write(changelsub(4),'(",  ix =",i3)') ix
                  write(changelsub(5),'(",  iy =",i3)') iy
                  changelsub(6) = " "
                end if
                write(changelsub(7),'(",  ie =",i3)') ie
                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                else
                  write(changelsub(8),'(",  it =",i3)') it
                end if
                write(changelsub(9),'(a1)') cha
c
                  changelsub(7) = " "
c
                if(iloopmode .eq. 1) then
                  changelsub(7) = " "
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

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

                     volsum =0.0d0
                     do iii = izi, izi+nzstepi-1
                     do ii  = iyi, iyi+nystepi-1
                     do i   = ixi, ixi+nxstepi-1
                       volsum = volsum + vl(i,ii,iii)
                     end do
                     end do
                     end do

                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, volsum, xm(ix), xm(ix+nxstepi),
     &                  ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi)
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
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

*-----------------------------------------------------------------------
*        x axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 3 ) then

               inum = 0

            do imi = 1, nm
            do iei = 1, ne, nestepi
            do iyi = 1, ny, nystepi
            do izi = 1, nz, nzstepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ie = iei
               iy = iyi
               iz = izi
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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "iy  =",i3,3x,
     &            "iz  =",i3,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie, iy, iz,
     &                        eb(ie), eb(ie+nestepi),
     &                        ym(iy), ym(iy+nystepi),
     &                        zm(iz), zm(iz+nzstepi)
               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3,/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: x [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Flux ",a29)') hsunit(itunt(m))

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

       include 'samepage_include/flux_petxyzm_x.inc'

*-----------------------------------------------------------------------
              if(iloopmode .eq. 0) then
               if( itmlp(m) .eq. 0 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ie, iy, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, iy, iz, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                   cha, inum, ie, iy, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iy =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,",  it =",i3,a1)')
     &                   cha, inum, ie, iy, iz, itmnt(m,im), it, cha
                  end if
               end if

             else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                  changelsub(3) = " "
                else
                  write(changelsub(3),'(",  mset =",i3)') itmnt(m,im)
                end if

                changelsub(4) = " "
                write(changelsub(5),'(",  iy =",i3)') iy
                write(changelsub(6),'(",  iz =",i3)') iz
                if(iloopmode .eq. 3) then
                end if
                if(iloopmode .eq. 4) then
                  changelsub(5) = " "
                end if
                if(iloopmode .eq. 5) then
                  changelsub(6) = " "
                end if
                write(changelsub(7),'(",  ie =",i3)') ie
                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                else
                  write(changelsub(8),'(",  it =",i3)') it
                end if
                write(changelsub(9),'(a1)') cha
c
                if(iloopmode .eq. 1) then
                  changelsub(7) = " "
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


               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, eb(ie), eb(ie+nestepi),
     &                  ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi)
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
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

*-----------------------------------------------------------------------
*        y axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 4 ) then

               inum = 0

            do imi = 1, nm
            do iei = 1, ne, nestepi
            do ixi = 1, nx, nxstepi
            do izi = 1, nz, nzstepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               im = imi
               ie = iei
               ix = ixi
               iz = izi
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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "ix  =",i3,3x,
     &            "iz  =",i3,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie, ix, iz,
     &                        eb(ie), eb(ie+nestepi),
     &                        xm(ix), xm(ix+nxstepi),
     &                        zm(iz), zm(iz+nzstepi)
               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3,/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: y [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Flux ",a29)') hsunit(itunt(m))

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

       include 'samepage_include/flux_petxyzm_y.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( itmlp(m) .eq. 0 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ie, ix, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, ix, iz, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                   cha, inum, ie, ix, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,",  it =",i3,a1)')
     &                   cha, inum, ie, ix, iz, itmnt(m,im), it, cha
                  end if
               end if

             else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                  changelsub(3) = " "
                else
                  write(changelsub(3),'(",  mset =",i3)') itmnt(m,im)
                end if

                write(changelsub(4),'(",  ix =",i3)') ix
                changelsub(5) = " "
                write(changelsub(6),'(",  iz =",i3)') iz
                if(iloopmode .eq. 3) then
                  changelsub(4) = " "
                end if
                if(iloopmode .eq. 4) then
                end if
                if(iloopmode .eq. 5) then
                  changelsub(6) = " "
                end if
                write(changelsub(7),'(",  ie =",i3)') ie
                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                else
                  write(changelsub(8),'(",  it =",i3)') it
                end if
                write(changelsub(9),'(a1)') cha
c
                if(iloopmode .eq. 1) then
                  changelsub(7) = " "
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


               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, eb(ie), eb(ie+1),
     &                  xm(ix), xm(ix+nxstepi), zm(iz), zm(iz+nzstepi)
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
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

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 5 ) then

               inum = 0

            do imi = 1, nm
            do iei = 1, ne, nestepi
            do ixi = 1, nx, nxstepi
            do iyi = 1, ny, nystepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               ie = iei
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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "ix  =",i3,3x,
     &            "iy  =",i3,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie, ix, iy,
     &                        eb(ie), eb(ie+nestepi),
     &                        xm(ix), xm(ix+nxstepi),
     &                        ym(iy), ym(iy+nystepi)
               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3,/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+ntstepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z[cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Flux ",a29)') hsunit(itunt(m))

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

       include 'samepage_include/flux_petxyzm_z.inc'

*-----------------------------------------------------------------------

              if(iloopmode .eq. 0) then
               if( itmlp(m) .eq. 0 ) then
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
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                  cha, inum, ie, ix, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,",  it =",i3, a1)')
     &                  cha, inum, ie, ix, iy, itmnt(m,im), it, cha
                  end if
               end if

             else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                  changelsub(3) = " "
                else
                  write(changelsub(3),'(",  mset =",i3)') itmnt(m,im)
                end if

                write(changelsub(4),'(",  ix =",i3)') ix
                write(changelsub(5),'(",  iy =",i3)') iy
                changelsub(6) = " "
                if(iloopmode .eq. 3) then
                  changelsub(4) = " "
                end if
                if(iloopmode .eq. 4) then
                  changelsub(5) = " "
                end if
                if(iloopmode .eq. 5) then
                end if
                write(changelsub(7),'(",  ie =",i3)') ie
                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                else
                  write(changelsub(8),'(",  it =",i3)') it
                end if
                write(changelsub(9),'(a1)') cha
c
                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                end if
                if(iloopmode .eq. 1) then
                  changelsub(7) = " "
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

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]")')
     &                     yen, eb(ie), eb(ie+nestepi),
     &                  xm(ix), xm(ix+nxstepi), ym(iy), ym(iy+nystepi)
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
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

*-----------------------------------------------------------------------
*        time axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 11 ) then

               inum = 0

            do imi = 1, nm
            do ixi = 1, nx, nxstepi
            do iyi = 1, ny, nystepi
            do izi = 1, nz, nzstepi
            do iei = 1, ne, nestepi
            do ipi = 1, np, npstepi
               im = imi
               ix = ixi
               iy = iyi
               iz = izi
               ie = iei
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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ix  =",i3,3x,
     &            "iy  =",i3,3x,
     &            "iz  =",i3,/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ix, iy, iz,
     &                        xm(ix), xm(ix+nxstepi),
     &                        ym(iy), ym(iy+nystepi),
     &                        zm(iz), zm(iz+nzstepi)
               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)
                  write(iot,'(
     &            "#  ie =",i3,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ie, eb(ie), eb(ie+nestepi)

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Time [nsec]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Flux ",a29)') hsunit(itunt(m))

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

       include 'samepage_include/flux_petxyzm_time.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( itmlp(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &               ", (ix,iy,iz) = (",i3,
     &               ",",i3,",",i3,"),  ie =",i3,a1)')
     &               cha, inum, ix, iy, iz, ie, cha
               else
                     write(iot,'(/a1,"no. =",i3,
     &               ", (ix,iy,iz) = (",i3,
     &               ",",i3,",",i3,")",
     &               ",  mset =",i3,",  ie =",i3,a1)')
     &               cha, inum, ix, iy, iz, itmnt(m,im), ie, cha
               end if

             else
                write(changelsub(1),'(a1,"no. =",i3)') cha,inum
                write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
                if( itmlp(m) .eq. 0 ) then
                  changelsub(3) = " "
                else
                  write(changelsub(3),'(",  mset =",i3)') itmnt(m,im)
                end if

                write(changelsub(4),'(", (ix,iy,iz) = (",i3,",",
     &               i3,",",i3,")")') ix, iy, iz
                changelsub(5) = " "
                changelsub(6) = " "
                if(iloopmode .eq. 3) then
                  changelsub(4) = " "
                  write(changelsub(5),'(",  iy =",i3)') iy
                  write(changelsub(6),'(",  iz =",i3)') iz
                end if
                if(iloopmode .eq. 4) then
                  write(changelsub(4),'(",  ix =",i3)') ix
                  changelsub(5) = " "
                  write(changelsub(6),'(",  iz =",i3)') iz
                end if
                if(iloopmode .eq. 5) then
                  write(changelsub(4),'(",  ix =",i3)') ix
                  write(changelsub(5),'(",  iy =",i3)') iy
                  changelsub(6) = " "
                end if
                write(changelsub(7),'(",  ie =",i3)') ie
                if( ittty(m) .eq. 0 ) then
                  changelsub(8) = " "
                else
                  write(changelsub(8),'(",  it =",i3)') it
                end if
                write(changelsub(9),'(a1)') cha
c
                  changelsub(8) = " "
c
                if(iloopmode .eq. 1) then
                  changelsub(7) = " "
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

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                     volsum =0.0d0
                     do iii = izi, izi+nzstepi-1
                     do ii  = iyi, iyi+nystepi-1
                     do i   = ixi, ixi+nxstepi-1
                       volsum = volsum + vl(i,ii,iii)
                     end do
                     end do
                     end do

                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, volsum, xm(ix), xm(ix+nxstepi),
     &                  ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi)
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
               end if
                  write(iot,'(
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     eb(ie), eb(ie+nestepi)
                     write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
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
            do it = 1, ntg

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

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "iz  =",i3,3x,
     &            "msos = ",a8,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie, iz, chq(ip),
     &                        eb(ie), eb(ie+1),
     &                        zm(iz), zm(iz+1)
               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3,/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+1)
               end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               if( itmlp(m) .eq. 0 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,a1)')
     &                     cha, inum, ie, iz, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, iz, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, iz, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iz =",i3,
     &                     ",  mset =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, iz, itmnt(m,im), it, cha
                  end if
               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

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

!OBINATA(2012.7.2): output *.err
               write(iot,'(1p10e11.3)')
     &         ( ( tr(ip,ie,it,icf(ix,iy,iz),im,ioe),
     &             ix = 1, nx ), iy = ny, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# x          y        ",
     &                      "  flux       r.err")')

               do iy = 1, ny
               do ix = 1, nx

                  write(iot,'(1p3e13.4,0pf8.4)')
     &              xm(ix)  + rtxdl(m)/2.0,
     &              ym(iy)  + rtydl(m)/2.0,
     &              (tr(ip,ie,it,icf(ix,iy,iz),im,k),k=1,2)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

                  write(iot,'("#   x = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7/
     &                        "#   y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            xm(1) + rtxdl(m)/2.0, xm(nx) + rtxdl(m)/2.0, rtxdl(m),
     &            ym(1) + rtydl(m)/2.0, ym(ny) + rtydl(m)/2.0, rtydl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'y/x',( xm(ix) + rtxdl(m)/2.0, ix = 1, nx )

               do iy = ny, 1, -1

!OBINATA(2012.7.2): output *.err
                  write(iot,'(1p1000e11.3)')
     &            ym(iy) + rtydl(m)/2.0,
     &            ( tr(ip,ie,it,icf(ix,iy,iz),im,ioe), ix = 1, nx )

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

      if ( ioe .eq. 2 ) then

          write(iot,'("y: Relative Error")')

      else if( itazl(m) .eq. 0 ) then

          write(iot,'("y: Flux ",a29)') hsunit(itunt(m))

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
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]"/
     &                     "  msos &=& ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     zm(iz), zm(iz+1), chq(ip)
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,it)
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

*-----------------------------------------------------------------------
*        yz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 ) then

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
            do it = 1, ntg

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

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "ix  =",i3,3x,
     &            "msos = ",a8,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie, ix, chq(ip),
     &                        eb(ie), eb(ie+1),
     &                        xm(ix), xm(ix+1)
               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3,/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+1)
               end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               if( itmlp(m) .eq. 0 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,a1)')
     &                     cha, inum, ie, ix, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, ix, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, ix, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  ix =",i3,
     &                     ",  mset =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, ix, itmnt(m,im), it, cha
                  end if
               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

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


               write(iot,'(/"#   no. =",i3,3x,
     &         "ie  =",i3,3x,
     &         "ix  =",i3,3x/
     &         "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, ie, ix,
     &                     eb(ie), eb(ie+1),
     &                     xm(ix), xm(ix+1)

*-----------------------------------------------------------------------

            if( ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) .and.
     &            igsh .eq. 0 ) then

                  write(iot,'(/a4," y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7," ; x = ",
     &            1p1g14.7," to ",1p1g14.7," by ",
     &            1p1g14.7," ;")') dc2,
     &            ym(ny) + rtydl(m)/2.0, ym(1) + rtydl(m)/2.0, rtydl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

!OBINATA(2012.7.2): output *.err
               write(iot,'(1p10e11.3)')
     &         ( ( tr(ip,ie,it,icf(ix,iy,iz),im,ioe),
     &                iz = 1, nz ), iy = ny, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# y          z        ",
     &                      "  flux       r.err")')

               do iz = 1, nz
               do iy = 1, ny

                  write(iot,'(1p3e13.4,0pf8.4)')
     &              ym(iy)  + rtydl(m)/2.0,
     &              zm(iy)  + rtzdl(m)/2.0,
     &              (tr(ip,ie,it,icf(ix,iy,iz),im,k),k=1,2)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

                  write(iot,'("#   y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7/
     &                        "#   z = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            ym(1) + rtydl(m)/2.0, ym(ny) + rtydl(m)/2.0, rtydl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'y/z',( zm(iz) + rtzdl(m)/2.0, iz = 1, nz )

               do iy = ny, 1, -1

!OBINATA(2012.7.2): output *.err
                  write(iot,'(1p1000e11.3)')
     &            ym(iy) + rtydl(m)/2.0,
     &            ( tr(ip,ie,it,icf(ix,iy,iz),im,ioe), iz = 1, nz )

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

      if ( ioe .eq. 2 ) then

          write(iot,'("y: Relative Error")')

      else if( itazl(m) .eq. 0 ) then

          write(iot,'("y: Flux ",a29)') hsunit(itunt(m))

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
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  msos &=& ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     xm(ix), xm(ix+1), chq(ip)
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
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

*-----------------------------------------------------------------------
*        xz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9 ) then

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
            do it = 1, ntg

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

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "iy  =",i3,3x,
     &            "msos = ",a8,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie, iy, chq(ip),
     &                        eb(ie), eb(ie+1),
     &                        ym(iy), ym(iy+1)
               if( itmlp(m) .gt. 0 )
     &            write(iot,'("#   mset =",i3)') itmnt(m,im)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3,/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     it, tb(it), tb(it+1)
               end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               if( itmlp(m) .eq. 0 ) then
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iy =",i3,a1)')
     &                     cha, inum, ie, iy, cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iy =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, iy, it, cha
                  end if
               else
                  if( ittty(m) .eq. 0 ) then
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,a1)')
     &                     cha, inum, ie, iy, itmnt(m,im), cha
                  else
                     write(iot,'(/a1,"no. =",i3,
     &                     ",  ie =",i3,
     &                     ",  iy =",i3,
     &                     ",  mset =",i3,",  it =",i3,a1)')
     &                     cha, inum, ie, iy, itmnt(m,im), it, cha
                  end if
               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

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

!OBINATA(2012.7.2): output *.err
               write(iot,'(1p10e11.3)')
     &         ( ( tr(ip,ie,it,icf(ix,iy,iz),im,ioe),
     &                iz = 1, nz ), ix = nx, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# x          z        ",
     &                      "  flux       r.err")')

               do iz = 1, nz
               do ix = 1, nx

                  write(iot,'(1p3e11.3,0pf8.4)')
     &               xm(ix)  + rtxdl(m)/2.0,
     &               zm(iz)  + rtzdl(m)/2.0,
     &               tr(ip,ie,it,icf(ix,iy,iz),im,1),
     &               tr(ip,ie,it,icf(ix,iy,iz),im,2)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

                  write(iot,'("#   x = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7/
     &                        "#   z = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            xm(1) + rtxdl(m)/2.0, xm(nx) + rtxdl(m)/2.0, rtxdl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'x/z',( zm(iz) + rtzdl(m)/2.0, iz = 1, nz )

               do ix = nx, 1, -1

!OBINATA(2012.7.2): output *.err
                  write(iot,'(1p1000e11.3)')
     &            xm(ix) + rtxdl(m)/2.0,
     &            ( tr(ip,ie,it,icf(ix,iy,iz),im,ioe), iz = 1, nz )

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
               iaxs = 3
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

      if ( ioe .eq. 2 ) then

          write(iot,'("y: Relative Error")')

      else if( itazl(m) .eq. 0 ) then

          write(iot,'("y: Flux ",a29)') hsunit(itunt(m))

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
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  msos &=& ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     ym(iy), ym(iy+1), chq(ip)
               if( itmlp(m) .gt. 0 ) then
                     write(iot,'(
     &                     "  mset  &=& ",i3)')
     &                     itmnt(m,im)
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

*-----------------------------------------------------------------------

         end if

            call prestart(m,iot) !OBINATA(2012.6.13)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

  900    continue

*-----------------------------------------------------------------------
         deallocate( ixyz )
      include 'samepage_include/samepage999.inc'

! sumover
      deallocate (vl_x,vl_y,vl_z)

      return
      end
