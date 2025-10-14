************************************************************************
*                                                                      *
      subroutine setmd(iom,jom,ierr)
*                                                                      *
*       setup material and nuclear data                                *
*       Last modified by K.Niita on 2009/10/05                         *
*                                                                      *
************************************************************************
      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      use dedx_file
      use moddas_material

      implicit real*8 (a-h,o-z)

      parameter ( imfnmax=1000) ! T.Sato 2022/10/30

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------

      common /nnmode/ nmode
      common /ndemax/ dnmax(20)
      common /kmat1h/ kmatg(kvlmax)
      common /eparm/  esmax, esmin, emin(20)
      common /tparm/  tmax(20)
      common /matmd/  nmlow, nmintl, nmfinl
      common /kmat1a/ mxmat, mxmat0, mxnel
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /regdm/  idmg(kvlmax)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)
      common /kmat1j/ intum
      common /paraj/  mstz(300), parz(300)
      common /tmpreg/ dtmp(kvlmax)
      common /pwtreg/ dpwt(kvlmax)
      common /kmat1o/ iom1, iom2, iom3
      common /kmat1p/ kmout
      common /brsmsg/ mnbrs, icbrs, mbbrs(kvlmax), cbrem(49)
      common /paran/  icfn(100), ilfn(100), chfn(100)
      character       chfn*200

      common /wtcntl/ iwt, icimp(20), ifcls(20), iwwin(20), ircls(20)

      common /inggs/  iog, igcel, ioa, igsuf, iob, igtrs

      common /mpi00/  npe, me

      common /egsemi/ iegsemi, iegsout

      dimension ix(3)
      character ht*10


*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall33/ itln(itlmax,2), itli(itlmax,2), itlr(itlmax,2),
     &                rtdm(itlmax,2)
      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)
      common /fmcard/ ifm, ifmi(imfnmax,3)

      common /dpnmaxcom/ dpnmax ! photon nuclear library maximum energy

      common /  tscreg   / ntscell(kvlmax)

*-----------------------------------------------------------------------

               ierr = 0
               call cputime(2)

               mxa = igcel

*-----------------------------------------------------------------------
*     mode of neutron, photon and electron, proton
*-----------------------------------------------------------------------

            do i = 1, mipt

               in_erg(i) = 8
               kpt(i)    = 0

            end do

               ides = mstz(37)

c             dpnmax=150.0 !!!!! temporary !!!!!!!
             if ( dpnmax.gt.emin(14) ) then
               ispn = 1
             end if

            if( dnmax( 2) .gt. emin( 2) ) kpt(1) = 1
            if( dnmax(14) .gt. emin(14) ) kpt(2) = 1
! T.Sato 2021/08/31 not necessary to read electron library for negs=1
            if( dnmax(12) .gt. emin(12) .and. mstz(85).le.0) kpt(3) = 1
            if( dnmax( 1) .gt. emin( 1) ) kpt(9) = 1
            if( dnmax(15) .gt. emin(15) ) kpt(31) = 1
            if( dnmax(18) .gt. emin(18) ) kpt(34) = 1

*-----------------------------------------------------------------------

! Emission particle type. 1:n, 9:p, 31:d, 32:t, 33:3he, 34:a
            if( ( kpt(9) .eq. 1 .and. dnmax(1) .gt. 150.0 ) .or.
     &          ( kpt(1) .eq. 1 .and. dnmax(2) .gt. 150.0 ) ) then

               kpt(20) = 1
               kpt(21) = 1

            end if

            if( ( kpt(9) .eq. 1 .and. dnmax(1) .gt. 20.0 ) .or.
     &          ( kpt(1) .eq. 1 .and. dnmax(2) .gt. 20.0 ) ) then

               kpt(31) = 1
               kpt(32) = 1
               kpt(33) = 1
               kpt(34) = 1

            end if

*-----------------------------------------------------------------------

            if( dnmax(2) .gt. 21.0 .and. kpt(1) .eq. 1 .and.
     &          kpt(9) .eq. 0 ) kpt(9) = 2
            if( dnmax(1) .gt. 21.0 .and. kpt(9) .eq. 1 .and.
     &          kpt(1) .eq. 0 ) kpt(1) = 2

*-----------------------------------------------------------------------
*        DPA library in tally
*-----------------------------------------------------------------------

               do k = 1, itnm
               do j = 1, 2

                  if( itln(k,j) .gt. 0 ) then

                     if( rtdm(k,j) .lt. 0.0 ) then

                        rtdm(k,j) = dnmax(j)

                     end if

                     if( rtdm(k,j) .gt. emin(j) ) then

                        if( j .eq. 1 ) then

                           if( kpt(9) .eq. 0 ) kpt(9) = 2

                        else

                           if( kpt(1) .eq. 0 ) kpt(1) = 2

                        end if

                     end if

                  end if

               end do
               end do

*-----------------------------------------------------------------------
*        FM card
*-----------------------------------------------------------------------

         if( ifm .gt. 0 ) then

*-----------------------------------------------------------------------

            do j = 1, ifm

                  k = ifmi(j,1)

                  if( kpt(k) .eq. 0 ) then

                     kpt(k) = 2

                  end if

            end do

*-----------------------------------------------------------------------

            do m = 1, itnm

               if( itmlp(m) .gt. 0 ) then

                  do k = 1, itmlp(m)

                     do 55 l = 1, itmpn(m,k)

                        ityp = itmpt(m,k,l,1)

                        if( ityp .eq. 1 ) then
                           ipim = 9
                        else if( ityp .eq. 2 ) then
                           ipim = 1
                        else if( ityp .eq. 14 ) then
                           ipim = 2
                        else
                           goto 55
                        end if

                        if( rtmme(m,k) .lt. 0.0 ) then

                           rtmme(m,k) = dnmax(ityp)

                        end if

   55                continue

                  end do

               end if

            end do

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*        mode summary
*-----------------------------------------------------------------------

               if( kpt(2) .eq. 0 .and. kpt(3) .eq. 0 ) ides = 1

               nmode = 0

            if( kpt(1) .ne. 0 ) nmode = nmode + 1
            if( kpt(2) .ne. 0 ) nmode = nmode + 1
            if( kpt(3) .ne. 0 ) nmode = nmode + 1

            if( kpt(2) .ne. 0 .and. kpt(3) .eq. 0 .and.
     &          ides .eq. 0 ) nmode = nmode + 1

            if( kpt(9) .ne. 0 ) nmode = nmode + 1

            if( kpt(31) .ne. 0 ) nmode = nmode + 1
            if( kpt(32) .ne. 0 ) nmode = nmode + 1
            if( kpt(33) .ne. 0 ) nmode = nmode + 1
            if( kpt(34) .ne. 0 ) nmode = nmode + 1

            if( ispn .ne. 0 .and. kpt(2) .eq. 0 ) ispn = 0
            if( ispn .ne. 0 ) nmode = nmode + 1

*-----------------------------------------------------------------------
*        nmode = 0, return
*-----------------------------------------------------------------------

            if( nmode .eq. 0 ) then

               nmlow  = 0
               goto 900

            end if

*-----------------------------------------------------------------------
*        cut-off energy
*-----------------------------------------------------------------------

            ecf(1)  = emin( 2)
            ecf(2)  = emin(14)
            ecf(3)  = emin(12)
            ecf(4)  = emin(6)

         do i = 5, 8
            ecf(i)  = emin(11)
         end do

            ecf(9)  = emin(1)

         do i = 10, 19
            ecf(i)  = emin(11)
         end do

            ecf(20)  = emin(3)
            ecf(21)  = emin(4)
            ecf(22)  = emin(8)
            ecf(23)  = emin(9)
            ecf(24)  = emin(9)

         do i = 25, 30
            ecf(i)  = emin(11)
         end do

            ecf(31)  = emin(15)
            ecf(32)  = emin(16)
            ecf(33)  = emin(17)
            ecf(34)  = emin(18)

*-----------------------------------------------------------------------

cKN??? 20.0-0.1 for neutron case, it is stopped !!

            if( ecf(1) .gt. 0.01 ) ecf(1) = 0.01d0

            if( ecf(2) .lt. 0.001 )  ecf(2) = 0.001d0
            if( ecf(3) .lt. 0.001 )  ecf(3) = 0.001d0

            elc(1) = ecf(1)
            elc(2) = ecf(2)
            elc(3) = ecf(3)

            emin(12) = ecf(3)
            emin(13) = ecf(3)
            emin(14) = ecf(2)

            parz(12) = ecf(3)
            parz(13) = ecf(3)
            parz(14) = ecf(2)

            if( kpt(2) .eq. 0 ) ecf(2) = 0.001d0
            if( kpt(3) .eq. 0 ) ecf(3) = 0.001d0



*-----------------------------------------------------------------------
*        nmode is not zero; we need weight cut-off
*-----------------------------------------------------------------------

               iwt = iwt + 1

*-----------------------------------------------------------------------
*        open temporary file for information
*-----------------------------------------------------------------------

               iom1 = 23
               open(iom1,form='formatted',status='scratch')

*-----------------------------------------------------------------------
*     some unused default values
*-----------------------------------------------------------------------
*     mxt  = 1 : no time dependent of temperature
*     mcal = 0 : no mgopt card
*      igm = 0 : for mgopt card
*      img = 0 : for mgopt card
*     naw  = 0 : no awtab card
*     nxsc = 0 : no xs card
*   npikmt = 0 : no photo-production bias card
*     ntal = 0 : no tally card
*    npert = 0 : no perturbation card
*    ipert = 0 : no perturbation card
*   itotnu = 0 : no total fission card
*    itfxs = 0 : no fission turnoff card
*      nsr = 0 : for source
*      nsa = 0 : for source
*    ksdef = 0 : for source
*      mct = 0 : for prdmp card
*-----------------------------------------------------------------------

                mxt = 1
               mcal = 0
                igm = 0
                img = 0
                naw = 0
               nxsc = 0
             npikmt = 0
               ntal = 0
              npert = 0
              ipert = 0
             itotnu = 0

              itfxs = 0
              if( mstz(113) .eq. 0 ) itfxs = 1

                mct = 0

               nfer = 0
              lfatl = 0

                nsr = 0
                nsa = 0
              ksdef = 0

               itty = 5
               jtty = 6

               nkxs = 0
               kc8  = 0

               kcy  = 0
               ikz  = 0
               mjss = 0
               nilw = 0

              krflg = 0
              its30 = 0
               ntop = 0
               nstp = 0
               nwng = 0
              xunru = 0.
             istern = 0
                mkc = 0
                lmb = 0
                lxs = 0
                ngp = 0
               ntyn = 0
               mpan = 0
                iet = 0
               ipsc = 0
               ixre = 0
              ixcos = 0
               nter = 0
                irt = 0

*-----------------------------------------------------------------------
*      setup some parameters
*-----------------------------------------------------------------------

            jovr(1) = 1
            jovr(2) = 0
            jovr(3) = 1
            jovr(4) = 1
            jovr(5) = 0

            do i = 1, mink

               ink(i) = 0

            end do

*-----------------------------------------------------------------------
*        zero set
*-----------------------------------------------------------------------

            do i = 1, mipt
               dxw(i,1) = 0.
               dxw(i,2) = 0.
               dxw(i,3) = 0.
               nww(i) = 0
               ndx(i) = 0
               jgm(i) = 0
               nch(i) = 0
               mgegbt(i) = 0
            do j = 1, 5
            do k = 1, mxdx
               dxx(i,j,k) = 0.
            end do
            end do
            end do

            do i = 1, 5
               rnb(i) = 0.
            end do

            do i = 0, 50
               thgf(i) = 0.
            end do

            do i = 1, mtop
               rkt(i)  = 0.
            end do

            do i = 1, mbng
               rka(i) = 0.
            end do

            do i = 1, maxi
              calph(i) = 0.
            end do

*-----------------------------------------------------------------------
*     cut-off time
*-----------------------------------------------------------------------

               tco(1) = tmax( 2) / 10.0d0
               tco(2) = tmax(14) / 10.0d0
               tco(3) = tmax(12) / 10.0d0

*-----------------------------------------------------------------------
*     set debug information
*-----------------------------------------------------------------------

            do i = 1, 30

               dbcn(i) = 0.0d0

            end do

               dbcn( 5) = 600.
               dbcn( 9) = 1.e-4
               dbcn(10) = 100.
               dbcn(17) = parz(158)
               dbcn(18) = parz(159)

*-----------------------------------------------------------------------
*        set phys. parameters
*-----------------------------------------------------------------------

               emx(1)  = dnmax(2)
               emcf(1) = parz(151)
                         if( emcf(1) .lt. 0.0 ) goto 998
               iunr    = mstz(36)
                         if( iunr .ne. 0 ) iunr = 100000
               dnb     = parz(152)

               emx(2)  = dnmax(14)
               emcf(2) = max( parz(157), 0.001d0 )

               nocoh   = mstz(38)
                         if( nocoh .ne. 0 ) nocoh = 1

               emx(3)  = dnmax(12)
                         if( kpt(3) .eq. 0 ) emx(3) = emx(2)
cT.Sato 20140917, add icntl check on 2019/05/16
                         if( kpt(3) .ne. 0 .and. iegsemi.eq.0 .and.
     &              emx(3) .gt. 10000. and. mstz(1).eq.0) goto 997
               emcf(3) = 100.0
               iphot   = mstz(39)
                         if( kpt(2) .eq. 0 ) iphot = 1
               ibad    = mstz(40)
               istrg   = mstz(41)
                         if( kpt(3) .eq. 0 ) istrg = 0
               bnum    = parz(153)
                         if( abs(bnum-1.0d0) .le. 0.00001 ) bnum = 1.0d0
               xnum    = max(0.0d0,parz(154))
               rnok    = max(0.0d0,parz(155))
                         if( abs(parz(156)-1.0d0) .gt. 0.0001)
     &         enum    = max(0.0d0,parz(156))
                         if( bnum .ne. 0. ) enum = 1. / abs(bnum)
               numb    = mstz(42)

               if( bnum .ne. 0. .and. numb .ne. 0 ) goto 996

            do i = 4, 34
               emx(i) = 0.d0
            end do

! frtati 2024/06/06 bugfix MeV/u -> MeV
!            emx(1)=max(dnmax(1),dnmax(2),dnmax(15),dnmax(16),dnmax(17),
!     &      dnmax(18),dnmax(19),dpnmax) ! highest energy library used
            emx(1)=max(dnmax(1),dnmax(2),dnmax(15)*2.d0,dnmax(16)*3.d0,
     &      dnmax(17)*3.d0,dnmax(18)*4.d0,dnmax(19),dpnmax) ! highest energy library used

            emx(9)=emx(1)

               ex = emx(1)

            if( ( kpt(3) .ne. 0 .or. kpt(2) .ne. 0 ) .and.
     &            ides .eq. 0 ) ex = max(ex,emx(3))

            do mp = 4, mipt
               if( kpt(mp) .ne. 0 ) ex = max(ex,emx(mp))
            end do

            do mp = 4, mipt
               if( kpt(mp) .ne. 0 .and. emx(mp) .eq. 0.d0 )
     &         emx(mp) = ex
            end do

*-----------------------------------------------------------------------
*        brems bias
*-----------------------------------------------------------------------

            do i = 1, 49

               if( cbrem(i) .ne. 0.0d0 ) then

                  bbrem(i) = 1.0 / cbrem(i)

               else

                  bbrem(i) = 0.0d0

               end if

            end do

*-----------------------------------------------------------------------
*        setup input data
*-----------------------------------------------------------------------

               nmlow  = 0
               nmintl = mmmax

               lkxs = mmmax * ndp2

               nmat1 = mxmat0
               nmat  = mxmat0
               mxe1  = 0
               indt  = 0
               mix   = 0
               mnnm  = 0

         do m = 1, mxmat0

                     nel  = nint( das_kmatg(kmatg(m)+1) )
                     inlb = nint( das_kmatg(kmatg(m)+4) )
                     iplb = nint( das_kmatg(kmatg(m)+5) )
                     ielb = nint( das_kmatg(kmatg(m)+6) )
                     imts = nint( das_kmatg(kmatg(m)+8) )
                     know = nint( das_kmatg(kmatg(m)+9) )

                     iulb = nint( das_kmatg(kmatg(m)+10) )
                     ihlb = nint( das_kmatg(kmatg(m)+11) )

                     mix  = mix + nel
                     mnnm = max( mnnm, nel )

*-----------------------------------------------------------------------
*        dedxfile read
*-----------------------------------------------------------------------
               do ii=1,20
                idedx(ii) = nint( das_kmatg(kmatg(m)+11+ii) )
               enddo

               call char_int(2,dedx_filename(m),idedx)


               if( index(dedx_filename(m),' ') /= 0) then
                  call dedx_file_read(m)
               endif
*-----------------------------------------------------------------------
               do 100 l = 1, nel

                     izia = nint( das_kmatg(kmatg(m)+(l-1)*3+32) )
                     libi = nint( das_kmatg(kmatg(m)+(l-1)*3+34) )

                  do k = 1, 2*mxe1, 2

                     if( das(lkxs+k  ) .eq. izia .and.
     &                   das(lkxs+k+1) .eq. libi ) goto 100

                  end do

                     das(lkxs+2*mxe1+1) = izia
                     das(lkxs+2*mxe1+2) = libi

                     mxe1 = mxe1 + 1

  100          continue

            if( inlb .ne. 0 ) then

               do 101 l = 1, nel

                     izia = nint( das_kmatg(kmatg(m)+(l-1)*3+32) )
                     libi = nint( das_kmatg(kmatg(m)+(l-1)*3+34) )

                  do k = 1, 2*mxe1, 2

                     if( das(lkxs+k  ) .eq. izia .and.
     &                   das(lkxs+k+1) .eq. inlb ) goto 101

                  end do

                     das(lkxs+2*mxe1+1) = izia
                     das(lkxs+2*mxe1+2) = inlb
                     mxe1 = mxe1 + 1

  101          continue

            end if

            if( iplb .ne. 0 ) then

               do 102 l = 1, nel

                     izia = nint( das_kmatg(kmatg(m)+(l-1)*3+32) )
                     libi = nint( das_kmatg(kmatg(m)+(l-1)*3+34) )

                  do k = 1, 2*mxe1, 2

                     if( das(lkxs+k  ) .eq. izia .and.
     &                   das(lkxs+k+1) .eq. iplb ) goto 102

                  end do

                     das(lkxs+2*mxe1+1) = izia
                     das(lkxs+2*mxe1+2) = iplb
                     mxe1 = mxe1 + 1

  102          continue

            end if

            if( ielb .ne. 0 ) then

               do 103 l = 1, nel

                     izia = nint( das_kmatg(kmatg(m)+(l-1)*3+32) )
                     libi = nint( das_kmatg(kmatg(m)+(l-1)*3+34) )

                  do k = 1, 2*mxe1, 2

                     if( das(lkxs+k  ) .eq. izia .and.
     &                   das(lkxs+k+1) .eq. ielb ) goto 103

                  end do

                     das(lkxs+2*mxe1+1) = izia
                     das(lkxs+2*mxe1+2) = ielb
                     mxe1 = mxe1 + 1

  103          continue

            end if

            if( imts .ne. 0 ) then

                     indt = indt + imts

               do 104 l = 1, imts

                     ix(1) = nint( das_kmatg(kmatg(m)+know+(l-1)*3+1) )
                     ix(2) = nint( das_kmatg(kmatg(m)+know+(l-1)*3+2) )
                     ix(3) = nint( das_kmatg(kmatg(m)+know+(l-1)*3+3) )

                     call zaid(2,ht,ix)

                  kk = 0
                  ll = 0

                  do i = 1, 5
                     kk = kk * 40
     &                  + index('abcdefghijklmnopqrstuvwxyz0123456789,',
     &                    ht(i:i) ) + 1

                     ll = ll * 40
     &                  + index('abcdefghijklmnopqrstuvwxyz0123456789.',
     &                    ht(i+5:i+5) ) + 1

                  end do

                  do k = 1, 2*mxe1, 2

                     if( das(lkxs+k  ) .eq. kk .and.
     &                   das(lkxs+k+1) .eq. ll ) goto 104

                  end do

                     das(lkxs+2*mxe1+1) = kk
                     das(lkxs+2*mxe1+2) = ll
                     mxe1 = mxe1 + 1

  104          continue

            end if

         end do

*-----------------------------------------------------------------------

            mxe1 = mxe1 * max( 1, nmode )

*-----------------------------------------------------------------------
*        for electron and positron
*-----------------------------------------------------------------------

            nee = 0
            if( kpt(3) .ne. 0 .or. kpt(2) .ne. 0 .and. ides .eq. 0 )
     &      nee = in_erg(3) * log(emx(3)/ecf(3)) / log(two) + two

            d = in_erg(3)
            efac = half**(one/d)

*-----------------------------------------------------------------------
*        for proton
*-----------------------------------------------------------------------

            neep = 0

         if( kpt(9) .ne. 0 ) then

            neep = in_erg(9) * log(emx(9)/ecf(9)) / log(two) + two

         end if

*-----------------------------------------------------------------------

            mnnm = mnnm * ( kpt(1) + kpt(2) )
            ke   = max( kpt(3) , kpt(2) * ( 1 - ides ) )

*-----------------------------------------------------------------------
*     set das
*-----------------------------------------------------------------------

            mtasks = 1

            lrho = nmintl
            ltmp = lrho + mxa * mtasks
            lfme = ltmp + mxa * mxt
            ltth = lfme + mix * mtasks
            lawc = ltth + mxt
            lawt = lawc + mix
            lawn = lawt + naw
            ltbt = lawn + mxe1
            lden = ltbt + mxe1
            lesa = lden + mxa
            lpik = lesa + mxe1
            lgwt = lpik + npikmt * mtasks
            ldpt = lgwt + mxa * kpt(1) * kpt(2)
            lrpt = ldpt + 3 * npert * mnnm
            leba = lrpt + ipert
            leee = leba + nee * nmat1 * mtop
            lpru = leee + nee
            lpkn = lpru + nee * nmat1
            lfim = lpkn + nee * nmat1
            lrng = lfim + ( mipt + 1 ) * mxa
            ldrs = lrng + nee * nmat1
            lqav = ldrs + nee * nmat1
            lasp = lqav + nee * nmat1
            lear = lasp + nee * nmat1 * ( 1 - istrg )
            lqcn = lear + nee * nmat1 * ( 1 - istrg )
            lflc = lqcn + nee * nmat1 * ( 1 - istrg )
            lx85 = lflc + nee * nmat1 * ( 1 - istrg )
            legg = lx85 + nee * nmat1 * 10
            ledg = legg + nee * nmat1 * maxi
            leek = ledg + ke  * nmat1
            lwwk = leek + ke  * nmat1
            lpxr = lwwk + ke  * nmat1
            lxnm = lpxr + nee * nmat1
            lpbr = lxnm + ke  * nmat1
            lpbt = lpbr + nee * nmat1
            lebt = lpbt + nee * nmat1
            lebd = lebt + nee * nmat1 * mtop
            lfst = lebd + nee * nmat1 * mtop * nint( bbrem(1) )
            lftt = lfst + nee * nmat1 * nint( bbrem(1) )
            lech = lftt + nee * nmat1 * nint( bbrem(1) )
            lrtc = lech + mpng * mwng * ( nee / mstp + 1) * nmat1 * ke
            ltgp = lrtc !FURUTA + 15 * mxe1 * mtasks
            lemaxt = ltgp !FURUTA + mxe1 * mtasks
            lpnt = lemaxt + mxe1
            lqav0 = lpnt + nmat1
            lrng0 = lqav0 + nee * nmat1
            intm = lrng0 + nee * nmat1
            lintm = intm

*-----------------------------------------------------------------------

            intm = intm * ndp2

            lmat = intm
            llxd = lmat + mxa
            lnmt = llxd + nmat1 * mipt

            ljem = lnmt + nmat1 + 1
            lnsb = ljem + nmat1 * ke
            ljco = lnsb + nmat1 * ke
            liza = ljco + nmat1 * ke
            lizn = liza + mix
            llmn = lizn + mix
            lkmm = llmn + mix
            ljmd = lkmm + mix
            lnpq = ljmd + nmat1 + 1
            lmbi = lnpq + nmat1
            lkdr = lmbi + nmat1 * ke
            llme = lkdr + mxe1
            lkmt = llme + mipt * mix
            lixl = lkmt + 3 * indt
            lnty = lixl + 3 * mxe1
            lixc = lnty + mxe1
            lkxd = lixc + 61 * mxe1
            lkaw = lkxd + mxe1
            lipa = lkaw + naw
            ljmt = lipa + mxa + 1
            ljxs = ljmt + indt
            lnxs = ljxs + 32 * mxe1
            llfc = lnxs + 16 * mxe1
            llmt = llfc + mxa + 2
            lngm = llmt + mix
            lnht = lngm + mxe1
            lipt = lnht + mxe1
            lipb = lipt + 48 * ntal
            lnpt = lipb + (2+2*npkey)*npert
            ljpt = lnpt + npert + 1
            lktp = ljpt + 18 * ntal
            lkxs = lktp + mipt * ntal
            lktc = lkxs + mxe1
            ljss = lktc !FURUTA + 2 * mxe1 * mtasks
            lixs = ljss + mjss + nilw
            mnax = lixs + mixs * maxsec * mxe1

*-----------------------------------------------------------------------

            nrtc_bank=mxe1*mtasks !FURUTA
            nktc_bank=mxe1*mtasks !FURUTA
            ntgp_bank=mxe1*mtasks !FURUTA
            kpik = lpik

*-----------------------------------------------------------------------

            mnax = ( mnax + mod(mnax,2) ) / ndp2

            mmmax  = mnax + 1
            nmfinl = mnax + 1
            nmlow  = nmfinl - nmintl

            if( mnax .gt. mdas ) then

               write(iom,'(/
     &               ''<<< Memory ERROR : at the end of setup MD >>>''/
     &               ''*  memory exceeds mdas'',/
     &               ''*     start of MD memory ='',i9,/
     &               ''*       end of MD memory ='',i9,/
     &               ''*  total Material memory ='',i9,/
     &               ''*     total memory: mdas ='',i9,/
     &               ''<<<  Please extend mdas in param.inc >>>''/)')
     &                      nmintl, nmfinl, nmlow, mdas

               write(jom,'(/
     &               ''<<< Memory ERROR : at the end of setup MD >>>''/
     &               ''*  memory exceeds mdas'',/
     &               ''*     start of MD memory ='',i9,/
     &               ''*       end of MD memory ='',i9,/
     &               ''*  total Material memory ='',i9,/
     &               ''*     total memory: mdas ='',i9,/
     &               ''<<<  Please extend mdas in param.inc >>>''/)')
     &                      nmintl, nmfinl, nmlow, mdas

               goto 999

            end if
*-----------------------------------------------------------------------
       call ALLOCATE_GGMARRAY(mxa,mtasks,mxt,mix,naw,mxe1,npikmt,
     & kpt(1),kpt(2),nee,nmat1,istrg,ke,nint(bbrem(1)),mipt,indt,
     & mpng,mwng,mtop,maxi,ipert,ntal,npert,npkey,mnnm)
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*        initialize of dimensions
*-----------------------------------------------------------------------

            do i = nmintl + 1, mnax

               das(i) = 0.0

            end do


*-----------------------------------------------------------------------
*     set initial value
*-----------------------------------------------------------------------

            das(ljmd)   = 1
            das(ljmd+1) = 1
            do i = 1, mxmat0 + 1
               jmd(i) = das(ljmd + i - 1)
            end do

            do i = 1, mxe1
               kxs(i) = nint( das(mmmax*ndp2 + i) )
               nty(i) = nint( das(lnty + i) )
               kxd(i) = nint( das(lkxd + i) )
               do j = 1, 32
                  jxs(j,i) = nint( das(ljxs + j + 32*(i-1)))
               end do
               do j = 1, 16
                  nxs(j,i) = nint( das(lnxs + j + 16*(i-1)))
               end do
            end do

            do i = 1, mxa
               lfcl(i) = nint( das(llfc + i) )
            end do

            do i = 1, mix
               lmt(i) = nint( das(llmt + i) )
            end do


*-----------------------------------------------------------------------
*     set importance = 1.0 ( 0.0 for outer void ) for temporary
*-----------------------------------------------------------------------

         fim(mipt+1,1:mxa) = 0.d0
         do i = 1, mipt
         do j = 1, mxa

            if( idmg(j) .eq. -1 ) then

                  fim(i,j) = 0.0d0

            else

               if( kpt(i) .eq. 1 ) then

                  fim(i,j) = 1.0d0

               else

                  fim(i,j) = 0.0d0

               end if

            end if

         end do
         end do

*-----------------------------------------------------------------------
*     set cross section directory file name
*-----------------------------------------------------------------------

            xsdir  = chfn(7)(1:ilfn(7))
            hdpath = ' '
            hdpth  = ' '

            m = ichar(' ') + 256 * (ichar(' ') + 256 * ichar(' '))

         do i = 1, mxe1
         do j = 4, 60

            ixc(j,i) = m

         end do
         end do

*-----------------------------------------------------------------------
*     set up data
*-----------------------------------------------------------------------
         tth(1) = 0.0d0 
         do i = 1, mxa

            mat(i) = idmg(i)
            if( mat(i) .eq. -1 ) mat(i) = 0

            rho(i) = denr(i)
            tmp(i) = dtmp(i)

            if( kpt(1) * kpt(2) .ne. 0 ) gwt(i) = dpwt(i)

         end do

            m1 = ichar(' ') + 256 * (ichar(' ') + 256 * ichar(' '))
            m2 = ichar(' ') + 256 * (ichar(' ') + 256 * ichar('p'))
            m3 = ichar(' ') + 256 * (ichar(' ') + 256 * ichar('e'))
            m4 = ichar(' ') + 256 * (ichar(' ') + 256 * ichar('u'))
            m9 = ichar(' ') + 256 * (ichar(' ') + 256 * ichar('h'))

         do i = 1, nmat1

            lxd(1,i) = m1
            lxd(2,i) = m2
            lxd(3,i) = m3
            lxd(4,i) = m4
            lxd(9,i) = m9

            pnt(1) = huge

         end do


*-----------------------------------------------------------------------
*        material ID
*-----------------------------------------------------------------------

            nmat = nmat1

         do m = 1, mxmat0

            nmt(m) = idmn(m)
            if( nmt(m) .eq. -1 ) nmt(m) = 0

         end do

*-----------------------------------------------------------------------
*        check nlib, plib and elib
*-----------------------------------------------------------------------

         do m = 1, mxmat0

               inlb = nint( das_kmatg(kmatg(m)+4) )
               iplb = nint( das_kmatg(kmatg(m)+5) )
               ielb = nint( das_kmatg(kmatg(m)+6) )
               iulb = nint( das_kmatg(kmatg(m)+10) )
               ihlb = nint( das_kmatg(kmatg(m)+11) )

               ichn = inlb / 256**2
               ichp = iplb / 256**2
               iche = ielb / 256**2
               ichu = iulb / 256**2
               ichh = ihlb / 256**2

            if( inlb .ne. 0 .and. char(ichn) .ne. ' ' .and.
     &          char(ichn) .ne. 'c' .and. char(ichn) .ne. 'd' .and.
     &          char(ichn) .ne. 'y' .and. char(ichn) .ne. 'm' ) then

                  write(iom,'(/''Error : default nlib is wrong'')')
                  write(iom,'( ''        in material '',i4)') idmn(m)
                  write(jom,'(/''Error : default nlib is wrong'')')
                  write(jom,'( ''        in material '',i4)') idmn(m)
                  ierr = ierr + 1

            end if

            if( iplb .ne. 0 .and. char(ichp) .ne. ' ' .and.
     &          char(ichp) .ne. 'p' .and. char(ichp) .ne. 'g' ) then

                  write(iom,'(/''Error : default plib is wrong'')')
                  write(iom,'( ''        in material '',i4)') idmn(m)
                  write(jom,'(/''Error : default plib is wrong'')')
                  write(jom,'( ''        in material '',i4)') idmn(m)
                  ierr = ierr + 1

            end if

            if( ielb .ne. 0 .and. char(iche) .ne. ' ' .and.
     &          char(iche) .ne. 'e' ) then

                  write(iom,'(/''Error : default elib is wrong'')')
                  write(iom,'( ''        in material '',i4)') idmn(m)
                  write(jom,'(/''Error : default elib is wrong'')')
                  write(jom,'( ''        in material '',i4)') idmn(m)
                  ierr = ierr + 1

            end if

            if( iulb .ne. 0 .and. char(ichu) .ne. ' ' .and.
     &          char(ichu) .ne. 'u' ) then

                  write(iom,'(/''Error : default pnlib is wrong'')')
                  write(iom,'( ''        in material '',i4)') idmn(m)
                  write(jom,'(/''Error : default pnlib is wrong'')')
                  write(jom,'( ''        in material '',i4)') idmn(m)
                  ierr = ierr + 1

            end if

            if( ihlb .ne. 0 .and. char(ichh) .ne. ' ' .and.
     &          char(ichh) .ne. 'h' ) then

                  write(iom,'(/''Error : default hlib is wrong'')')
                  write(iom,'( ''        in material '',i4)') idmn(m)
                  write(jom,'(/''Error : default hlib is wrong'')')
                  write(jom,'( ''        in material '',i4)') idmn(m)
                  ierr = ierr + 1

            end if

            if( inlb .ne. 0 .and. char(ichn) .eq. 'y' .and.
     &          das_intum(intum+m) .gt. 0.0d0 ) then

                  write(iom,'(/''Error : dosimetry table is used '',
     &                         ''in material of cell'')')
                  write(iom,'( ''        in material '',i4)') idmn(m)
                  write(jom,'(/''Error : dosimetry table is used '',
     &                         ''in material of cell'')')
                  write(jom,'( ''        in material '',i4)') idmn(m)
                  ierr = ierr + 1

            end if

         end do

            if( ierr .ne. 0 ) goto 999

*-----------------------------------------------------------------------
*        input jemi, nsb, jcond, lxd
*-----------------------------------------------------------------------

         do m = 1, mxmat0

               igas = nint( das_kmatg(kmatg(m)+2) )
               istp = nint( das_kmatg(kmatg(m)+3) )
               inlb = nint( das_kmatg(kmatg(m)+4) )
               iplb = nint( das_kmatg(kmatg(m)+5) )
               ielb = nint( das_kmatg(kmatg(m)+6) )
               icnd = nint( das_kmatg(kmatg(m)+7) )
               iulb = nint( das_kmatg(kmatg(m)+10) )
               ihlb = nint( das_kmatg(kmatg(m)+11) )

            if( iplb .ne. 0 ) lxd(2,m) = iplb
            if( ielb .ne. 0 ) lxd(3,m) = ielb

cfrtati 2021/12/17 extensions for high-energy libraries
            if ( dnmax(2).gt.20.d0 ) then
              lxd(1,m) = mstz(150)
            end if
            lxd(4,m) = mstz(151)
            lxd(9,m) = mstz(149)
            lxd(31,m) = mstz(152)
            lxd(32,m) = mstz(153)
            lxd(33,m) = mstz(154)
            lxd(34,m) = mstz(155)

cfrtati 2021/12/17 inlb, iulb, ihlb are prior to high-energy default
            if( inlb .ne. 0 ) lxd(1,m) = inlb
            if( iulb .ne. 0 ) lxd(4,m) = iulb
            if( ihlb .ne. 0 ) lxd(9,m) = ihlb

            if( nee .gt. 0 ) then

               jemi(m) = igas
                nsb(m) = istp
              jcond(m) = icnd

            elseif( minval(ntscell(1:kvlmax)) .ne. 0 ) then ! 2022/7/21 ITSART is used.

               jemi(m) = igas

            end if

         end do

*-----------------------------------------------------------------------
*        nucleus information
*-----------------------------------------------------------------------

               mix = 0

         do m = 1, mxmat0

               nel = nint( das_kmatg(kmatg(m)+1) )

               nen = 0

            do l = 1, nel

               izia = nint( das_kmatg(kmatg(m)+(l-1)*3+32) )
               dens =       das_kmatg(kmatg(m)+(l-1)*3+33)
               libi = nint( das_kmatg(kmatg(m)+(l-1)*3+34) )

               if( dens .ne. 0. ) then

                  mix = mix + 1
                  nen = nen + 1

                  iza(mix) = izia
                  if( ispn .ne. 0) izn(mix) = izia
                  kmm(mix) = libi
                  fme(mix) = dens

               end if

            end do

               if( nen .eq. 0 ) goto 995

                  npq(m) = nen
                  jmd(1+m+1) = jmd(1+m) + npq(m)

! Not necessary to check duplicate isotope because they are automatically merged, T.Sato 2024/01/24
!            do i = jmd(1+m), jmd(1+m+1) - 1
!            do j = i+1, jmd(1+m+1) - 1
!
!               if( iza(i) .eq. iza(j) ) then
!
!                  write(iom,'(/''$** warning : duplicate isotope'',i8,
!     &                         '' within mateial '',i4)')
!     &            iza(i), idmn(m)
!
!               end if
!
!            end do
!            end do

         end do

*-----------------------------------------------------------------------
*        set s(a,b)
*-----------------------------------------------------------------------

               jndt = 0

         do m = 1, mxmat0

               imts = nint( das_kmatg(kmatg(m)+8) )
               know = nint( das_kmatg(kmatg(m)+9) )

            do l = 1, imts

               jndt = jndt + 1

               jmt(jndt) = idmn(m)

               kmt(1,jndt) = nint( das_kmatg(kmatg(m)+know+(l-1)*3+1) )
               kmt(2,jndt) = nint( das_kmatg(kmatg(m)+know+(l-1)*3+2) )
               kmt(3,jndt) = nint( das_kmatg(kmatg(m)+know+(l-1)*3+3) )

            end do

         end do

*-----------------------------------------------------------------------
*        brems bias
*-----------------------------------------------------------------------

            do i = 1, nmat1 * ke

               mbi(i) = nint( das(lmbi+i) )

            end do

            do i = 1, mnbrs

               mbi(i) = -mbbrs(i)

            end do

*-----------------------------------------------------------------------
*     set up the mateials
*-----------------------------------------------------------------------

               mxe = 0

            call setmat(iom,ierr)

               if( ierr .ne. 0 ) goto 900

*-----------------------------------------------------------------------
*     set up cross section tables
*-----------------------------------------------------------------------

            lxss = 0
            lexs = lxss

            call setxst(iom,jom,ierr) ! S.H. added jom (2017.1.4)

               if( ierr .ne. 0 ) goto 900

            nmfinl = mmmax
            nmlow  = nmfinl - nmintl

            if( mmmax .gt. mdas ) then

               write(iom,'(/
     &               ''<<< Memory ERROR : at the end of SETXST >>>''/
     &               ''*  memory exceeds mdas'',/
     &               ''*     start of MD memory ='',i9,/
     &               ''*       end of MD memory ='',i9,/
     &               ''*  total Material memory ='',i9,/
     &               ''*     total memory: mdas ='',i9,/
     &               ''<<<  Please extend mdas in param.inc >>>''/)')
     &                      nmintl, nmfinl, nmlow, mdas

               write(jom,'(/
     &               ''<<< Memory ERROR : at the end of SETXST >>>''/
     &               ''*  memory exceeds mdas'',/
     &               ''*     start of MD memory ='',i9,/
     &               ''*       end of MD memory ='',i9,/
     &               ''*  total Material memory ='',i9,/
     &               ''*     total memory: mdas ='',i9,/
     &               ''<<<  Please extend mdas in param.inc >>>''/)')
     &                      nmintl, nmfinl, nmlow, mdas

               goto 999

            end if

*-----------------------------------------------------------------------
*     mark any dr tables that have identical energy group structures.
*-----------------------------------------------------------------------

      do 50 ie=1,mxe
   50    nxs(16,ie)=ie
      do 80 ie=1,mxe
         if(nty(ie).ne.2.or.nxs(16,ie).gt.1000)go to 80
         nxs(16,ie)=1000+ie
      do 70 i=ie+1,mxe
         if(nty(i).ne.2.or.nxs(16,i).gt.1000.or.
     &    nxs(3,i).ne.nxs(3,ie))go to 70
      do 60 j=1,nxs(3,ie)
   60    if(xss(jxs(1,i)+j-1).ne.xss(jxs(1,ie)+j-1))go to 70
         nxs(16,i)=nxs(16,ie)
   70 continue
   80 continue

*-----------------------------------------------------------------------
*     initialize cross-section indexes.
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*     set maximum energy of library for data max
*-----------------------------------------------------------------------
      call datamaxsummary()

*-----------------------------------------------------------------------
*     set non-elastic cross-section maximum for high-energy event mode
*-----------------------------------------------------------------------
      call SET_LIBCSMAX()

*-----------------------------------------------------------------------

      goto 999

  995 continue
         ierr = 1
         write(iom,'(/''Error : no non-zero fractions. in '',i5)')
     &         idmn(m)
         write(jom,'(/''Error : no non-zero fractions. in '',i5)')
     &         idmn(m)
         goto 999

  996 continue
         ierr = 1
         write(iom,'(/''Error : only one kind of '',
     &                ''brems biasing allowed'')')
         write(jom,'(/''Error : only one kind of '',
     &                ''brems biasing allowed'')')
         goto 999

  997 continue
         ierr = 1
         write(iom,'(/''Error : electron dmax is > 10000.0'')')
         write(jom,'(/''Error : electron dmax is > 10000.0'')')
         goto 999

  998 continue
         ierr = 1
         write(iom,'(/''Error : neutron emcf is less than zero'')')
         write(jom,'(/''Error : neutron emcf is less than zero'')')
         goto 999

*-----------------------------------------------------------------------

  999 continue

            if( me .ne. 0 ) then

                  close( iom1 )

            else

               if( kmout .eq. 0 .and. ierr .eq. 0 ) then

                  close( iom1 )
                  open(iom1,form='formatted',status='scratch')

               end if

                  endfile( iom1 )

            end if

*-----------------------------------------------------------------------

  900 continue

               call cputime(2)

      return
      end

