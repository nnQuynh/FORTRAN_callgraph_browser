
************************************************************************
*                                                                      *
      subroutine check_tlet(m,iax,jsn,jsi,dsin,idsi,ill,ilf,ierr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************
      use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05
      use moddas
      use moddas_tally

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'param01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

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

      common /tall08/ rtrx0(itlmax), rtry0(itlmax)

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

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall21/ rtfac(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)

      common /tall28/ itgsh(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall31/ iterl(itlmax)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall37/ itcnt(9,itlmax)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)

      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall49/ itglt(itlmax)
      common /tall50/ itlmt(itlmax), ite2l(itlmax)   ! CCSE 2022.08.31

!OBINATA: for resfile
      common /tall59/ irfll(itlmax), crfln(itlmax), itrff(itlmax)
      character crfln*100
      character irfile*100

      dimension idas(mdas*2)
      equivalence ( das, idas )


*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character title*80
      character angelp*200
      character cxtxt*200
      character cytxt*200
      character cztxt*200

*-----------------------------------------------------------------------

      dimension lschn(32), ischn(32)
      character schan(32)*8

      data icsu / 32 /

      data ( schan(i), i = 1, 32 ) /
     &    'mesh    ','part    ','l-type  ','unit    ','axis    ',
     &    'file    ','title   ','angel   ','2d-type ','factor  ',
     &    'x-txt   ','y-txt   ','z-txt   ','gshow   ','rshow   ',
     &    'iechrl  ','material','epsout  ','ctmin(1)','ctmax(1)',
     &    'ctmin(2)','ctmax(2)','ctmin(3)','ctmax(3)','resol   ',
     &    'width   ','trcl    ','*trcl   ','gslat   ','letmat  ',
     &    'volmat  ','resfile '/

      data ( lschn(i), i = 1, 32 ) /
     &     4,         4,         6,         4,         4,
     &     4,         5,         5,         7,         6,
     &     5,         5,         5,         5,         5,
     &     6,         8,         6,         8,         8,
     &     8,         8,         8,         8,         5,
     &     5,         4,         5,         5,         6,
     &     6,         7/

*-----------------------------------------------------------------------

cfrtati 2021/10/05 6 -> mxpart
      dimension iptyp(mxpart), ipnkf(mxpart)
      dimension imtyp(mxpart,mxpart), imnkf(mxpart,mxpart)
      dimension jstyp(mxpart), jnkf0(mxpart)
      dimension icount(9)

      common /subtra/ isubt, ipsub(mxpart)   ! kitamura22/03/31

      dimension iaxis(6)
      character ifile(6)*100
      dimension lfile(6)

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      character dkam*9

      character tname*7
      data      tname /'[t-let]'/

      logical deqn4

*-----------------------------------------------------------------------

      character rglnrf*200, rglnech*200

      character*10 cepn(100)
      dimension lepn(100)
      character*20 ctl(2)

*-----------------------------------------------------------------------

      dimension jptyp(6), jpnkf(6)
      dimension kptyp(6,6), kpnkf(6,6)
      dimension vtrs(13)

*-----------------------------------------------------------------------

            ierr  = 0
            jpn   = 0

            iunt  = 1
            imate = 0
            jmate = 1
            inpat = 0
            inaxi = 0
            infil = 0
            idtyp = 3

            langel = 0
            lxtxt  = 0
            lytxt  = 0
            lztxt  = 0
            lgshow = 0
            lrshow = 0
            iechrl = 72
            matvol = 9
            ieps   = 0
            igkst  = 0
            idtt   = 0
            ktrs   = 0
            igslt  = 1
            letmat = 0
            lrfile = 0 !OBINATA
            irfflg = 0 !OBINATA

            icount(1) = 0
            icount(2) = 0
            icount(3) = 0
            icount(4) = -9999
            icount(5) =  9999
            icount(6) = -9999
            icount(7) =  9999
            icount(8) = -9999
            icount(9) =  9999

            rfact = 1.0

            ireso = 1
            width = 0.5

         do i = 1, icsu

            ischn(i) = 0

         end do

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( jpn  .eq. 3 ) goto 800

               if( iskip .ne. 0 ) goto 140

  150 continue

               if( jpn  .eq. 3 ) goto 800

*-----------------------------------------------------------------------
*        head of the section
*-----------------------------------------------------------------------

            if( chlw(i1:i1) .eq. '[' .and.
     &          chcm(i1:i4) .eq. tname ) then

               goto 140

            end if

*-----------------------------------------------------------------------
*        end of the section
*-----------------------------------------------------------------------

            if( chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 800

            end if

*-----------------------------------------------------------------------
*        identify the parameters
*-----------------------------------------------------------------------

            icl = i1

  200    continue

            chlc = chlw
            call chcomp(chlc,icl,i3,i5)

         do i = 1, icsu

            il = icl + lschn(i) - 1

            if( chlc(icl:il) .eq. schan(i)(1:lschn(i)) ) goto 100

         end do

               goto 800

*-----------------------------------------------------------------------
*        read value of parameters
*-----------------------------------------------------------------------

  100    continue

               ipm = i

               ischn( ipm ) = ischn( ipm ) + 1

            if( ipm .ne.  5 .and. ipm .ne. 6 .and.
     &          ipm .ne.  2 .and.
     &          ischn( ipm ) .gt. 1 ) goto 989

               ic = inumc(chlw,il+1,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 997

               icl = inumc(chlw,ic,i3,';') - 1

*-----------------------------------------------------------------------
*        mesh = region, r-z, or xyz
*-----------------------------------------------------------------------

         if( ipm .eq. 1 ) then

            if( chlw(ic:ic+2) .eq. 'reg' ) then

               imesh = 1

            else if( chlw(ic:ic+2) .eq. 'r-z' ) then

               imesh = 2

            else if( chlw(ic:ic+2) .eq. 'xyz' ) then

               imesh = 3

            else

               goto 998

            end if

            if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

            if( imesh .eq. 1 ) then

               call tregion0(1,jsn,jsi,dsin,idsi,ill,ilf,
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                      ntrn,mtrn,ndsm,nvol,ivl,irvl,0,
     &                      rglnrf)

                  if( jpn  .eq. 3 ) goto 800
                  if( ntrn .lt. -1 ) goto 996

            else if( imesh .eq. 2 ) then

               call trzmesh(jsn,jsi,dsin,idsi,ill,ilf,
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                      rzx0,rzy0,
     &                      irtp,inr,rmin,rmax,rdel,istrg,
     &                      iztp,inz,zmin,zmax,zdel,istzg)

                  if( jpn  .eq. 3 ) goto 800
                  if( irtp .lt. 0 ) goto 996
                  if( iztp .lt. 0 ) goto 996

            else if( imesh .eq. 3 ) then

               call txymesh(jsn,jsi,dsin,idsi,ill,ilf,
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                      ixtp,inx,xmin,xmax,xdel,istxg,
     &                      iytp,iny,ymin,ymax,ydel,istyg,
     &                      iztp,inz,zmin,zmax,zdel,istzg)

                  if( jpn  .eq. 3 ) goto 800
                  if( ixtp .lt. 0 ) goto 996
                  if( iytp .lt. 0 ) goto 996
                  if( iztp .lt. 0 ) goto 996

            end if

                  goto 150

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

         else if( ipm .eq. 2 ) then

  400       continue

               ic = jnumc(chlw,ic,icl)

               if( ic .gt. icl ) then

                  icl = jnumc(chlw,icl+2,i3)

                  if( icl .le. i3 ) goto 200

                  goto 140

               end if

*-----------------------------------------------------------------------

            call rdpname(ic,icl,chlw,istyp,inkf0,jstyp,jnkf0,ierr)

               if( ierr .eq. 994 ) goto 994
               if( ierr .eq. 998 ) goto 997

*-----------------------------------------------------------------------

               inpat = inpat + 1


                  iptyp(inpat) = istyp
                  ipnkf(inpat) = inkf0
                  ipsub(inpat) = isubt  ! kitamura22/03/31

               if( istyp .lt. 0 ) then

                  do i = 1, -istyp

                     imtyp(inpat,i) = jstyp(i)
                     imnkf(inpat,i) = jnkf0(i)

                  end do

               end if

               goto 400

*-----------------------------------------------------------------------
*        material
*-----------------------------------------------------------------------

         else if( ipm .eq. 17 ) then

            if( chlw(ic:ic+2) .eq. 'all' ) then

                  imate = 0

            else if( deqn4( chlw(ic:ic) ) ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 979

                  imate = nint( cvvv )

               if( imate .lt. 0 ) then

                  imate = -imate
                  jmate = -1

               end if

               if( imate .eq. 0 ) goto 979

                  nsmte = 1
                  call moddas_allocate_int(imate, ismte_temporary)

                  if( mmmax .gt. mdas ) goto 950

  151                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                        if( jpn  .eq. 3 ) goto 979

                        if( iskip .ne. 0 ) goto 151

                  ic = i1

               do k = 1, imate

                  if( ic .gt. i3 ) then

  152                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( jpn  .eq. 3 ) goto 979

                     if( iskip .ne. 0 ) goto 152

                     ic = i1

                  end if

                     ic = jnumc(chlw,ic,i3)

                     call snum(chlw,ic,i3,ic2,cvvv,ierr)

                     if( ierr .ne. 0 ) goto 979

                     matei = nint( cvvv )

                     ic = ic2

                     ismte_temporary(nsmte-1+k) = matei

               end do
               call moddas_deallocate_int(ismte_temporary)

            else

               goto 979

            end if

*-----------------------------------------------------------------------
*        let mesh
*-----------------------------------------------------------------------

         else if( ipm .eq. 3 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

               ietp = nint( cvvv )

               call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     'l',ietp,ine,emin,emax,edel,isteg)


               goto 150

*-----------------------------------------------------------------------
*        transform
*-----------------------------------------------------------------------

         else if( ipm .eq. 27 .or. ipm .eq. 28 ) then

                  if( ipm .eq. 28 ) ktrs = 1

               call ttrans(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     ic,ktrs,igkst,idtt,vtrs)


               goto 150

*-----------------------------------------------------------------------
*        unit
*-----------------------------------------------------------------------

         else if( ipm .eq. 4 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               iunt = nint( cvvv )

               if( iunt .lt. 1 .or.
     &             iunt .gt. 14 ) goto 988 ! T.Sato 2017/11/21

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        axis
*-----------------------------------------------------------------------

         else if( ipm .eq. 5 ) then

  600       inaxi = inaxi + 1

            if( inaxi .gt. 6 ) goto 991

            if( chlw(ic:ic+2) .eq. 'let' ) then

               iaxis(inaxi) = 1
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               iaxis(inaxi) = 2
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+1) .eq. 'xy' .or.
     &               chlw(ic:ic+1) .eq. 'yx' ) then

               iaxis(inaxi) = 7
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic+1) .eq. 'yz' .or.
     &               chlw(ic:ic+1) .eq. 'zy' ) then

               iaxis(inaxi) = 8
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic+1) .eq. 'zx' .or.
     &               chlw(ic:ic+1) .eq. 'xz' ) then

               iaxis(inaxi) = 9
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic+1) .eq. 'rz' .or.
     &               chlw(ic:ic+1) .eq. 'zr' ) then

               iaxis(inaxi) = 10
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic) .eq. 'x' ) then

               iaxis(inaxi) = 3
               ic =jnumc(chlw,ic+2,icl)

            else if( chlw(ic:ic) .eq. 'y' ) then

               iaxis(inaxi) = 4
               ic =jnumc(chlw,ic+2,icl)

            else if( chlw(ic:ic) .eq. 'z' ) then

               iaxis(inaxi) = 5
               ic =jnumc(chlw,ic+2,icl)

            else if( chlw(ic:ic) .eq. 'r' ) then

               iaxis(inaxi) = 6
               ic =jnumc(chlw,ic+2,icl)

            else

               goto 992

            end if

               if( ic .le. icl ) goto 600

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        file name
*-----------------------------------------------------------------------

         else if( ipm .eq. 6 ) then

  700          infil = infil + 1

               if( infil .gt. 6 ) goto 990

               icf = min( inumc(chlw,ic,icl,' ') - 1, icl )

               lfile(infil) = icf - ic + 1
               ifile(infil)(1:icf-ic+1) = chin(ic:icf)

               ic = jnumc(chlw,icf+2,icl)

               if( ic .le. icl ) goto 700

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        restart file name
*-----------------------------------------------------------------------
         else if( ipm .eq. 32 ) then

               icf = min( inumc(chlw,ic,icl,' ') - 1, icl )

               lrfile = icf-ic+1

               irfile(1:lrfile) = chin(ic:icf)

               irfflg = 1

*-----------------------------------------------------------------------
*        title
*-----------------------------------------------------------------------

         else if( ipm .eq. 7 ) then

               ict = min( ic + 79, i2 )

               title = chin(ic:ict)

               titll = ict - ic + 1

*-----------------------------------------------------------------------
*        angel parameters
*-----------------------------------------------------------------------

         else if( ipm .eq. 8 ) then

               ict = min( ic + 199, i2 )

               angelp = chin(ic:ict)

               langel = ict - ic + 1

*-----------------------------------------------------------------------
*        x-txt
*-----------------------------------------------------------------------

         else if( ipm .eq. 11 ) then

               ict = min( ic + 199, i2 )

               cxtxt = chin(ic:ict)

               lxtxt = ict - ic + 1

*-----------------------------------------------------------------------
*        y-txt
*-----------------------------------------------------------------------

         else if( ipm .eq. 12 ) then

               ict = min( ic + 199, i2 )

               cytxt = chin(ic:ict)

               lytxt = ict - ic + 1

*-----------------------------------------------------------------------
*        z-txt
*-----------------------------------------------------------------------

         else if( ipm .eq. 13 ) then

               ict = min( ic + 199, i2 )

               cztxt = chin(ic:ict)

               lztxt = ict - ic + 1

*-----------------------------------------------------------------------
*        gshow
*-----------------------------------------------------------------------

         else if( ipm .eq. 14 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               lgshow = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        gslat
*-----------------------------------------------------------------------

         else if( ipm .eq. 29 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               igslt = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        rshow
*-----------------------------------------------------------------------

         else if( ipm .eq. 15 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               lrshow = nint( cvvv )

               if( lrshow .le. 0 ) then

                  lrshow = 0

                  icl = jnumc(chlw,icl+2,i3)
                  if( icl .le. i3 ) goto 200

               end if

               if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

               call txymesh(jsn,jsi,dsin,idsi,ill,ilf,
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                      ixtp,inx,xmin,xmax,xdel,istxg,
     &                      iytp,iny,ymin,ymax,ydel,istyg,
     &                      iztp,inz,zmin,zmax,zdel,istzg)

                  if( jpn  .eq. 3 ) goto 800
                  if( ixtp .lt. 0 ) goto 996
                  if( iytp .lt. 0 ) goto 996
                  if( iztp .lt. 0 ) goto 996

                  goto 150

*-----------------------------------------------------------------------
*        reg echo length
*-----------------------------------------------------------------------

         else if( ipm .eq. 16 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               iechrl = nint( cvvv )

               if( iechrl .lt. 40 ) goto 997

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        2d-type
*-----------------------------------------------------------------------

         else if( ipm .eq. 9 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               idtyp = nint( cvvv )

               if( idtyp .lt. 1 .or. idtyp .gt. 7 ) goto 983

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        factor
*-----------------------------------------------------------------------

         else if( ipm .eq. 10 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               rfact = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        epsout
*-----------------------------------------------------------------------

         else if( ipm .eq. 18 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ieps = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        resolution
*-----------------------------------------------------------------------

         else if( ipm .eq. 25 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ireso = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        width
*-----------------------------------------------------------------------

         else if( ipm .eq. 26 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               width = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        counter
*-----------------------------------------------------------------------

         else if( ipm .ge. 19 .and. ipm .le. 24 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               icount((ipm-19+2)/2) = 1
               icount( ipm-19+4)    = nint( cvvv )

               if( abs( nint(cvvv) ) .gt. 9999 ) goto 978

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        letmat
*-----------------------------------------------------------------------

         else if( ipm .eq. 30 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               letmat = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        volmat
*-----------------------------------------------------------------------

         else if( ipm .eq. 31 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               matvol = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------

         end if

            goto 140

*-----------------------------------------------------------------------
*     check
*-----------------------------------------------------------------------
  800 continue

      iec = 0

      call check_mesh(m,imesh,iec,cepn,lepn,ierr)

      if ( itmsh(m) .eq. 1 ) then

        call check_reg(m,rglnrf,iec,cepn,lepn,ierr)

      else if ( itmsh(m) .eq. 2 ) then

        call check_x0y0(m,rzx0,rzy0,iec,cepn,lepn,ierr)

        call check_type('r',iec,cepn,lepn,
     &                  itrty(m),rtrma(m),rtrmi(m),itrnm(m),
     &                  irtp, rmax, rmin, inr)

        call check_type('z',iec,cepn,lepn,
     &                  itzty(m),rtzma(m),rtzmi(m),itznm(m),
     &                  iztp, zmax, zmin, inz)

      else if ( itmsh(m) .eq. 3 ) then

        call check_type('x',iec,cepn,lepn,
     &                  itxty(m),rtxma(m),rtxmi(m),itxnm(m),
     &                  ixtp, xmax, xmin, inx)

        call check_type('y',iec,cepn,lepn,
     &                  ityty(m),rtyma(m),rtymi(m),itynm(m),
     &                  iytp, ymax, ymin, iny)

        call check_type('z',iec,cepn,lepn,
     &                  itzty(m),rtzma(m),rtzmi(m),itznm(m),
     &                  iztp, zmax, zmin, inz)

      end if

      call check_type('l',iec,cepn,lepn,
     &                itety(m),rtema(m),rtemi(m),itenm(m),
     &                ietp, emax, emin, ine)

      call check_unit(m,iunt,iec,cepn,lepn,ierr)

      call check_axis(m,iaxis,inaxi,iec,cepn,lepn,ierr)

      do i = 1, itaxn(m)

        if ( any( itaxs(m,i) .eq. (/ 7, 8, 9, 10 /))) then
          call check_2dtype(m,idtyp,iec,cepn,lepn,ierr)
        end if

      end do

      call check_factor(m,rfact,iec,cepn,lepn,ierr)

      call check_part(m,iptyp,ipnkf,inpat,iec,cepn,lepn,ierr)

      call check_letmat(m,letmat,iec,cepn,lepn,ierr)

*-----------------------------------------------------------------------

      call judge_tall_check(iec,cepn,lepn,dsin(jsn),idsi(jsn),ierr)

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  950    write(dkam,'(i9)') mdas
         m_err = 'Total tally storage number exceeds mdas ='//dkam
         ErrCha = ''
         ErrID = 'L:1020/R:check_tlet/F:restlet.f'
         goto 999

  976    m_err = 'LET mesh is Lethergy but emin is zero'
         ErrCha = ''
         ErrID = 'L:1025/R:check_tlet/F:restlet.f'
         goto 999

  978    m_err = 'counter should be from -9999 to 9999'
         ErrCha = ''
         ErrID = 'L:1030/R:check_tlet/F:restlet.f'
         goto 999

  979    m_err = 'Description of material parameter is wrong in tally '
     &            //tname
         ErrCha = ''
         ErrID = 'L:1036/R:check_tlet/F:restlet.f'
         goto 999

  983    m_err = '2D-type should be 1-7 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1041/R:check_tlet/F:restlet.f'
         goto 999

  986    write(dkam,'(i5)') itlmax
         m_err = 'Total tally number exceeds itlmax = '//dkam
         ErrCha = ''
         ErrID = 'L:1047/R:check_tlet/F:restlet.f'
         goto 999

  987    m_err = 'Unknown parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1052/R:check_tlet/F:restlet.f'
         goto 999

  988    m_err = 'Unit should be 1 - 12 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1057/R:check_tlet/F:restlet.f'
         goto 999

  989    m_err = 'Double definition of the parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1062/R:check_tlet/F:restlet.f'
         goto 999

  990    m_err = 'Too many file name in tally '//tname
         ErrCha = ''
         ErrID = 'L:1067/R:check_tlet/F:restlet.f'
         goto 999

  991    m_err = 'Too many axis in tally '//tname
         ErrCha = ''
         ErrID = 'L:1072/R:check_tlet/F:restlet.f'
         goto 999

  992    m_err = 'Unknown axis name in tally '//tname
         ErrCha = ''
         ErrID = 'L:1077/R:check_tlet/F:restlet.f'
         goto 999

  993    m_err = 'In this line, [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:1082/R:check_tlet/F:restlet.f'
         goto 999

  994    m_err = 'Name of particle is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1087/R:check_tlet/F:restlet.f'
         goto 999

  995    m_err = 'Number of particles is larger than 6 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1092/R:check_tlet/F:restlet.f'
         goto 999

  996    m_err = 'Description of mesh is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1097/R:check_tlet/F:restlet.f'
         goto 999

  997    m_err = 'Description of parameter is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1102/R:check_tlet/F:restlet.f'
         goto 999

  998    m_err = 'Unknown mesh parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1107/R:check_tlet/F:restlet.f'
         goto 999

*-----------------------------------------------------------------------

  999 continue

         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1

*-----------------------------------------------------------------------

      return

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_tlet(m,iax,ierr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        use RESTALMOD
        use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
      use moddas_tally

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'
        include 'param01.inc'

*-----------------------------------------------------------------------

        common /mpi00/ npe, me

*-----------------------------------------------------------------------

        common /talout/ itall
        common /talmm/  nmmax, lmmax, itlmx

*-----------------------------------------------------------------------

        common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------

        common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
        common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                  rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
        common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &                  rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)
        common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &                  rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)
        common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                  rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)
        common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                  rtemi(itlmax), rtema(itlmax), rtedl(itlmax)
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
        common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
        common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

        dimension     idas(mdas*2)
        equivalence ( das, idas )

*-----------------------------------------------------------------------

        character dsin(0:9,2)*200
        dimension idsi(0:9,2)
        dimension ill(0:9,2), ilf(0:9,2)
        dimension jsn(2), jsi(2)
        dimension ierrs(2)

*-----------------------------------------------------------------------

        ierr = 0
        noe  = 1

        if ( any( itaxs(m,iax) .eq. (/ 7, 8, 9, 10 /) )
     &      .and. ittwo(m) .ne. 4 ) noe = 2

*-----------------------------------------------------------------------
*   open restart file
*-----------------------------------------------------------------------
        call open_resfile(m,noe,jsn,jsi,dsin,idsi,ill,ilf,newtall,ierr)

        if ( newtall .ne. 0 ) goto 900  !! it's new tally
        if ( ierr    .ne. 0 ) goto 900

*-----------------------------------------------------------------------
*   check tally
*-----------------------------------------------------------------------
        do 800 ioe = 1, noe

        if(ireschk.eq.0) then  ! T.Sato 2013/10/19
        call check_tlet(m,iax,jsn(ioe),jsi(ioe),
     &                 dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                 ierr)
        endif

        if (ierr .ne. 0 ) goto 900

*-----------------------------------------------------------------------

        rewind(jsi(ioe))
        ill(jsn(ioe),ioe) = 1

*-----------------------------------------------------------------------

        idas0 = nmmax
        idas1 = lmmax
        if (itmsh(m) .eq. 1 ) then
          idas2 = idas1 + itenm(m)
          idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
          idas4 = idas2 + itrgn(m)
     &          + ( itrgn(m) + mod(itrgn(m),2) ) / 2
          idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
          idasa = idas4
        else if (itmsh(m) .eq. 3 ) then
          idas2 = ( idas1 + itenm(m) - 1 ) * 2 + 1
          idasa = idas1 + itenm(m)
        else if (itmsh(m) .eq. 3 ) then
          idasa = idas1 + itenm(m)
        end if

*-----------------------------------------------------------------------
*   read tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call read_letreg(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itpan(m),itrgn(m),itrgm(m),itenm(m),
     &      idas_itreg(itreg(m)),
     &      das_iterg(iterg(m)),
     &      trRES(irestalm(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        else if ( itmsh(m) .eq. 2 ) then

          call read_letrz(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itpan(m),itrnm(m),itznm(m),itenm(m),
     &      das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),das_iterg(iterg(m)),
     &      trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 3 ) then

          call read_letxyz(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itpan(m),
     &      itxnm(m),itynm(m),itznm(m),itenm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)),
     &      das_iterg(iterg(m)),
     &      trRES(irestalm(m)))

        end if

  800   continue
*-----------------------------------------------------------------------
*   restore tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call restore_letreg(m,
     &      itpan(m),itrgn(m),itrgm(m),itenm(m),
     &      idas_itreg(itreg(m)),
     &      das_iterg(iterg(m)),
     &      trRES(irestalm(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        else if ( itmsh(m) .eq. 2 ) then

          call restore_letrz(m,
     &      itpan(m),itrnm(m),itznm(m),itenm(m),
     &      das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),das_iterg(iterg(m)),
     &      trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 3 ) then

          call restore_letxyz(m,
     &      itpan(m),
     &      itmtn(m),ismte(itmtt(m)),
     &      itxnm(m),itynm(m),itznm(m),itenm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)),
     &      das_iterg(iterg(m)),
     &      trRES(irestalm(m)))

        end if

  900   continue
*-----------------------------------------------------------------------
*   close restart file
*-----------------------------------------------------------------------

        do ioe = 1, noe

          close(jsi(ioe))

        end do

*-----------------------------------------------------------------------
  999   continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_letreg(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                       np,nr,mr,ne,kr,eb,tr,
     &                       nvl,ivl,rvl,
     &                       nx,ny,nz,xm,ym,zm)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'
        include 'param01.inc'

*-----------------------------------------------------------------------

        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)

        common /tall21/ rtfac(itlmax)
        common /fact01/ facmax(itlmax)
        common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------

        dimension   kr(mr)
        dimension   eb(ne+1)
        dimension   ew(ne)
        dimension   vl(nr)
        dimension   lr(nr)
        dimension   tr(np,ne,nr,2)

        dimension   ivl(nvl)
        dimension   rvl(nvl)

        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

        character chin*200, chlw*200, chcm*200

        real(8),allocatable:: totnorm_nn(:) ! T.Sato 2025/03/13

! sumover dummy
        dimension tott_sum(1,2)
        data nsame /1/

*-----------------------------------------------------------------------

        character dc2*20

*-----------------------------------------------------------------------
*        LET axis
*-----------------------------------------------------------------------
         if( itaxs(m,iax) .eq. 1 ) then

           dc2  = '#let-lwrlet-upr'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31
           nr_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_per.inc'

           do 190 ir = 1, nr, nrstepi
           do 190 ip = 1, np, npstepi

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             nestepi = 1
             do ie = 1, ne, nestepi


               read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &              dmm0, dmm1,
     &      ((((tr(ip+iploop-1,ie+ieloop-1,
     &               ir+irloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      irloop=1,nrstepi)

             end do

! sumover
            ie = 1
            call pletreg_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nrstepi,
     &      ip,ie,ir,nsame,tott_sum)

! T.Sato 2025/03/12 for re-normalization
          allocate (totnorm_nn(nsame))
          read(jsi,'(26x,1000(1pe13.4,8x))')
     &    (totnorm_nn(i),i=1,nsame)
          do ie = 1, ne
           do i=1,nsame
            if(nrstepi.ne.1) then ! samepage = reg
             tr(ip,ie,i,1)=tr(ip,ie,i,1)*totnorm_nn(i)
            elseif(npstepi.ne.1) then ! samepage = part
             tr(i,ie,ir,1)=tr(i,ie,ir,1)*totnorm_nn(i)
            else  ! only one normalization factor (nsame = 1)
             tr(ip,ie,ir,1)=tr(ip,ie,ir,1)*totnorm_nn(i)
            endif
           enddo
          enddo
          deallocate (totnorm_nn)

  190     continue

*-----------------------------------------------------------------------
*        reg axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 2 ) then

           dc2  = '#numregvolume'
           ldc2 = 13

           facmx = 1.d0  ! kitamura23/03/31
           nr_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_per.inc'

           do 290 ie = 1, ne, nestepi
           do 290 ip = 1, np, npstepi

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             nrstepi = 1
             do ir = 1, nr, nrstepi


               read(jsi,'(i5,1x,i7,1pe13.4,1000(1pe13.4,0pf8.4))')
     &              idmm0, idmm1, dmm2,
     &      ((((tr(ip+iploop-1,ie+ieloop-1,
     &               ir+irloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      irloop=1,nrstepi)

             end do

! sumover
            ir = 1
            call pletreg_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nrstepi,
     &      ip,ie,ir,nsame,tott_sum)


  290     continue

*-----------------------------------------------------------------------
         end if

  900 continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_letrz(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                      np,nr,nz,ne,rm,zm,eb,tr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'
        include 'param01.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)

        common /tall21/ rtfac(itlmax)
        common /fact01/ facmax(itlmax)
        common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------

        dimension   rm(nr+1)
        dimension   zm(nz+1)
        dimension   eb(ne+1)
        dimension   ew(ne)
        dimension   tr(np,ne,nr,nz,2)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

        character chin*200, chlw*200, chcm*200

        real(8),allocatable:: totnorm_nn(:) ! T.Sato 2025/03/13

! sumover dummy
        dimension tott_sum(1,2)
        data nsame /1/

*-----------------------------------------------------------------------

        character dc2*20

*-----------------------------------------------------------------------
*        LET axis
*-----------------------------------------------------------------------
         if( itaxs(m,iax) .eq. 1 ) then

           dc2  = '#let-lwrlet-upr'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31
           nz_0 = 0
           nr_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_perz.inc'

           do 190 ir = 1, nr, nrstepi
           do 190 iz = 1, nz, nzstepi
           do 190 ip = 1, np, npstepi
             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             nestepi = 1
             do ie = 1, ne, nestepi


               read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &              dmm0, dmm1,
     &      (((((tr(ip+iploop-1,ie+ieloop-1,
     &               ir+irloop-1,iz+izloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)

             end do

! sumover
            ie = 1
            call pletrz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nrstepi,nzstepi,
     &      ip,ie,ir,iz,nsame,tott_sum)

! T.Sato 2025/03/12 for re-normalization
          allocate (totnorm_nn(nsame))
          read(jsi,'(26x,1000(1pe13.4,8x))')
     &    (totnorm_nn(i),i=1,nsame)
          do ie = 1, ne
           do i=1,nsame
            if(nrstepi.ne.1) then ! samepage = r
             tr(ip,ie,i,iz,1)=tr(ip,ie,i,iz,1)*totnorm_nn(i)
            elseif(nzstepi.ne.1) then ! samepage = z
             tr(ip,ie,ir,i,1)=tr(ip,ie,ir,i,1)*totnorm_nn(i)
            elseif(npstepi.ne.1) then ! samepage = part
             tr(i,ie,ir,iz,1)=tr(i,ie,ir,iz,1)*totnorm_nn(i)
            else  ! only one normalization factor (nsame = 1)
             tr(ip,ie,ir,iz,1)=tr(ip,ie,ir,iz,1)*totnorm_nn(i)
            endif
           enddo
          enddo
          deallocate (totnorm_nn)

  190     continue

*-----------------------------------------------------------------------
*        r axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 6 ) then

           dc2  = '#r-lowerr-upper'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31
           nz_0 = 0
           nr_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_perz.inc'

           do 290 ie = 1, ne, nestepi
           do 290 iz = 1, nz, nzstepi
           do 290 ip = 1, np, npstepi
             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             nrstepi = 1
             do ir = 1, nr, nrstepi


               read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &              dmm0, dmm1,
     &      (((((tr(ip+iploop-1,ie+ieloop-1,
     &               ir+irloop-1,iz+izloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)

             end do

! sumover
            ir = 1
            call pletrz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nrstepi,nzstepi,
     &      ip,ie,ir,iz,nsame,tott_sum)

  290     continue

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 5 ) then

           dc2  = '#z-lowerz-upper'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31
           nz_0 = 0
           nr_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_perz.inc'

           do 390 ie = 1, ne, nestepi
           do 390 ir = 1, nr, nrstepi
           do 390 ip = 1, np, npstepi

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             nzstepi = 1
             do iz = 1, nz, nzstepi


               read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &              dmm0, dmm1,
     &      (((((tr(ip+iploop-1,ie+ieloop-1,
     &               ir+irloop-1,iz+izloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)

             end do

! sumover
            iz = 1
            call pletrz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nrstepi,nzstepi,
     &      ip,ie,ir,iz,nsame,tott_sum)

  390     continue

*-----------------------------------------------------------------------
*        rz axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 10 ) then

           call sel_dc2(m,dc2,ldc2)
           if( ittwo(m) .eq. 4 ) then
             dc2  = '#rzletr.err'
             ldc2 = 11
           else if( ittwo(m) .eq. 5 ) then
             dc2  = 'r/z'
             ldc2 = 3
           end if

           facmx = 1.d0  ! kitamura23/03/31

           do 490 ip = 1, np
           do 490 ie = 1, ne

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

*-----------------------------------------------------------------------
             if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               read(jsi,'(1p10e11.3)')
     &         ( ( tr(ip,ie,ir,iz,ioe), iz = 1, nz ), ir = nr, 1, -1 )

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 4 ) then

               do iz = 1, nz
               do ir = 1, nr

                 read(jsi,'(1p3e11.3,0pf8.4)')
     &                dmm0, dmm1,
     &                tr(ip,ie,ir,iz,1), tr(ip,ie,ir,iz,2)

               end do
               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 5 ) then

               do ir = nr, 1, -1

                  read(jsi,'(1p1000e11.3)')
     &            dmm1,
     &            ( tr(ip,ie,ir,iz,ioe), iz = 1, nz )

               end do

*-----------------------------------------------------------------------
             end if

  490     continue

*-----------------------------------------------------------------------
         end if

  900 continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_letxyz(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                       np,nx,ny,nz,ne,xm,ym,zm,eb,tr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'
        include 'param01.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)

        common /tall21/ rtfac(itlmax)
        common /fact01/ facmax(itlmax)
        common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------

        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)
        dimension   eb(ne+1)
        dimension   ew(ne)
        dimension   tr(np,ne,nx*ny*nz,2)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

        character chin*200, chlw*200, chcm*200

        real(8),allocatable:: totnorm_nn(:) ! T.Sato 2025/03/13

! sumover dummy
        dimension tott_sum(1,2)
        data nsame /1/
*-----------------------------------------------------------------------

        character dc2*20

*-----------------------------------------------------------------------

        icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny

        npg = np
        neg = ne

*-----------------------------------------------------------------------
*        LET axis
*-----------------------------------------------------------------------
         if( itaxs(m,iax) .eq. 1 ) then

           dc2  = '#let-lwrlet-upr'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_pexyz.inc'

           do 190 ix = 1, nx, nxstepi
           do 190 iy = 1, ny, nystepi
           do 190 iz = 1, nz, nzstepi
           do 190 ip = 1, np, npstepi

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             nestepi = 1
             do ie = 1, ne, nestepi


               read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &              dmm0, dmm1,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,
     &       icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &        k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi)

             end do

! sumover
            ie = 1
             call pletxyz_sumover_getput(jsi,m,iax,
     &       npstepi,nestepi,nxstepi,nystepi,nzstepi,
     &       ip,ie,ix,iy,iz,maxtott,tott_sum)

! T.Sato 2025/03/12 for re-normalization
          allocate (totnorm_nn(nsame))
          read(jsi,'(26x,1000(1pe13.4,8x))')
     &    (totnorm_nn(i),i=1,nsame)
          do ie = 1, ne
           do i=1,nsame
            if(nxstepi.ne.1) then ! samepage = x
      tr(ip,ie,icf(i,iy,iz),1)=tr(ip,ie,icf(i,iy,iz),1)*totnorm_nn(i)
            elseif(nxstepi.ne.1) then ! samepage = y
      tr(ip,ie,icf(ix,i,iz),1)=tr(ip,ie,icf(ix,i,iz),1)*totnorm_nn(i)
            elseif(nzstepi.ne.1) then ! samepage = z
      tr(ip,ie,icf(ix,iy,i),1)=tr(ip,ie,icf(ix,iy,i),1)*totnorm_nn(i)
            elseif(npstepi.ne.1) then ! samepage = part
      tr(i,ie,icf(ix,iy,iz),1)=tr(i,ie,icf(ix,iy,iz),1)*totnorm_nn(i)
            else  ! only one normalization factor (nsame = 1)
      tr(ip,ie,icf(ix,iy,iz),1)=tr(ip,ie,icf(ix,iy,iz),1)*totnorm_nn(i)
            endif
           enddo
          enddo
          deallocate (totnorm_nn)

  190     continue

*-----------------------------------------------------------------------
*        x axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 3 ) then

           dc2  = '#x-lowerx-upper'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_pexyz.inc'

           do 290 ie = 1, ne, nestepi
           do 290 iy = 1, ny, nystepi
           do 290 iz = 1, nz, nzstepi
           do 290 ip = 1, np, npstepi

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             nxstepi = 1
             do ix = 1, nx, nxstepi


               read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &              dmm0, dmm1,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,
     &       icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &        k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi)

             end do

! sumover
            ix = 1
             call pletxyz_sumover_getput(jsi,m,iax,
     &       npstepi,nestepi,nxstepi,nystepi,nzstepi,
     &       ip,ie,ix,iy,iz,maxtott,tott_sum)

  290     continue

*-----------------------------------------------------------------------
*        y axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 4 ) then

           dc2  = '#y-lowery-upper'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_pexyz.inc'

           do 390 ie = 1, ne, nestepi
           do 390 ix = 1, nx, nxstepi
           do 390 iz = 1, nz, nzstepi

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             nystepi = 1
             do iy = 1, ny, nystepi


               read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &              dmm0, dmm1,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,
     &       icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &        k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi)

             end do

! sumover
            iy = 1
             call pletxyz_sumover_getput(jsi,m,iax,
     &       npstepi,nestepi,nxstepi,nystepi,nzstepi,
     &       ip,ie,ix,iy,iz,maxtott,tott_sum)

  390     continue

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 5 ) then

           dc2  = '#z-lowerz-upper'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_pexyz.inc'

           do 490 ie = 1, ne, nestepi
           do 490 ix = 1, nx, nxstepi
           do 490 iy = 1, ny, nystepi
           do 490 ip = 1, np, npstepi

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             nzstepi = 1
             do iz = 1, nz, nzstepi


               read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &              dmm0, dmm1,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,
     &       icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &        k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi)

             end do

! sumover
            iz = 1
             call pletxyz_sumover_getput(jsi,m,iax,
     &       npstepi,nestepi,nxstepi,nystepi,nzstepi,
     &       ip,ie,ix,iy,iz,maxtott,tott_sum)

  490     continue


*-----------------------------------------------------------------------
*        xy axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 7 ) then

           call sel_dc2(m,dc2,ldc2)
           if( ittwo(m) .eq. 4 ) then
             dc2  = '#xyletr.err'
             ldc2 = 11
           else if( ittwo(m) .eq. 5 ) then
             dc2  = 'y/x'
             ldc2 = 3
           end if

           facmx = 1.d0  ! kitamura23/03/31

           do 590 iz = 1, nz
           do 590 ip = 1, npg
           do 590 ie = 1, neg

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

*-----------------------------------------------------------------------
             if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               read(jsi,'(1p10e11.3)')
     &         ( ( tr(ip,ie,icf(ix,iy,iz),ioe),
     &             ix = 1, nx ), iy = ny, 1, -1 )

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 4 ) then

               do iy = 1, ny
               do ix = 1, nx

                 read(jsi,'(1p3e13.4,0pf8.4)')
     &                dmm0, dmm1,
     &                (tr(ip,ie,icf(ix,iy,iz),k),k=1,2)

               end do
               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 5 ) then

               do iy = ny, 1, -1

                  read(jsi,'(1p1000e11.3)')
     &            dmm1,
     &            ( tr(ip,ie,icf(ix,iy,iz),ioe), ix = 1, nx )

               end do

*-----------------------------------------------------------------------
             end if

  590     continue

*-----------------------------------------------------------------------
*        yz axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 8 ) then

           call sel_dc2(m,dc2,ldc2)
           if( ittwo(m) .eq. 4 ) then
             dc2  = '#yzletr.err'
             ldc2 = 11
           else if( ittwo(m) .eq. 5 ) then
             dc2  = 'y/z'
             ldc2 = 3
           end if

           facmx = 1.d0  ! kitamura23/03/31

           do 690 ix = 1, nx
           do 690 ip = 1, npg
           do 690 ie = 1, neg

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

*-----------------------------------------------------------------------
             if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               read(jsi,'(1p10e11.3)')
     &         ( ( tr(ip,ie,icf(ix,iy,iz),ioe),
     &             iz = 1, nz ), iy = ny, 1, -1 )

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 4 ) then

               do iz = 1, nz
               do iy = 1, ny

                 read(jsi,'(1p3e13.4,0pf8.4)')
     &                dmm0, dmm1,
     &                (tr(ip,ie,icf(ix,iy,iz),k),k=1,2)

               end do
               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 5 ) then

               do iy = ny, 1, -1

                  read(jsi,'(1p1000e11.3)')
     &            dmm1,
     &            ( tr(ip,ie,icf(ix,iy,iz),ioe), iz = 1, nz )

               end do

*-----------------------------------------------------------------------
             end if

  690     continue

*-----------------------------------------------------------------------
*        xz axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 9 ) then

           call sel_dc2(m,dc2,ldc2)
           if( ittwo(m) .eq. 4 ) then
             dc2  = '#xzletr.err'
             ldc2 = 11
           else if( ittwo(m) .eq. 5 ) then
             dc2  = 'x/z'
             ldc2 = 3
           end if

           facmx = 1.d0  ! kitamura23/03/31

           do 790 iy = 1, ny
           do 790 ip = 1, npg
           do 790 ie = 1, neg

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

*-----------------------------------------------------------------------
             if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               read(jsi,'(1p10e11.3)')
     &         ( ( tr(ip,ie,icf(ix,iy,iz),ioe),
     &             iz = 1, nz ), ix = nx, 1, -1 )

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 4 ) then

               do iz = 1, nz
               do ix = 1, nx

                 read(jsi,'(1p3e11.3,0pf8.4)')
     &                dmm0, dmm1,
     &                (tr(ip,ie,icf(ix,iy,iz),k),k=1,2)

               end do
               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 5 ) then

               do ix = nx, 1, -1

                  read(jsi,'(1p1000e11.3)')
     &            dmm1,
     &            ( tr(ip,ie,icf(ix,iy,iz),ioe), iz = 1, nz )

               end do

*-----------------------------------------------------------------------
             end if

  790     continue

*-----------------------------------------------------------------------
         end if

  900 continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_letreg(m,
     &                          np,nr,mr,ne,kr,eb,tr,
     &                          nvl,ivl,rvl,
     &                          nx,ny,nz,xm,ym,zm)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'
        include 'param01.inc'

*-----------------------------------------------------------------------

        common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall21/ rtfac(itlmax)

        common /fact01/ facmax(itlmax)  ! kitamura23/03/31

*-----------------------------------------------------------------------

        dimension   kr(mr)
        dimension   eb(ne+1)
        dimension   ew(ne)
        dimension   vl(nr)
        dimension   lr(nr)
        dimension   tr(np,ne,nr,2)
        dimension   ivl(nvl)
        dimension   rvl(nvl)

        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)

*-----------------------------------------------------------------------

            if( itunt(m) .eq. 1 .or. itunt(m) .eq. 2 .or.
     &          itunt(m) .eq. 7 .or. itunt(m) .eq. 8 .or.
     &          itunt(m) .eq.13 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else
     &      if( itunt(m) .eq. 3 .or. itunt(m) .eq. 4 .or.
     &          itunt(m) .eq. 9 .or. itunt(m) .eq. 10.or.
     &          itunt(m) .eq.14 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

            else

               do i = 1, ne

                  ew(i) = 1.d0

               end do
               ew_sum = 1.0d0

            end if

*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

               call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)
               
               vl_sum= sum(vl(:))

*-----------------------------------------------------------------------
*        ( unit = 1-6 ; vol = 1.0 )
*-----------------------------------------------------------------------

            if( itunt(m) .le. 6 .or. itunt(m) .ge. 13) then ! T.Sato 2025/03/13

               do ir = 1, nr

                  vl(ir) = 1.0d0

               end do
               vl_sum = 1.0d0

            end if

*-----------------------------------------------------------------------

        if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

        do 100 ir = 1, nr
        do 100 ie = 1, ne
        do 100 ip = 1, np

           call invert_stdev(m,A,B,
     &                     tr(ip,ie,ir,1),
     &                     tr(ip,ie,ir,2),
     &      vl(ir)*ew(ie)/abs(rtfac(m)/facmax(m)))

           tr(ip,ie,ir,1) = A
           tr(ip,ie,ir,2) = B

! sumover
           fact_in = abs(rtfac(m)/facmax(m))
           call pletreg_sumover_stdev(1,m,ip,ie,ir,
     &          fact_in,ew(ie),vl(ir),ew_sum,vl_sum)


  100   continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_letrz(m,
     &                         np,nr,nz,ne,rm,zm,eb,tr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'
        include 'param01.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall21/ rtfac(itlmax)

        common /fact01/ facmax(itlmax)  ! kitamura23/03/31

*-----------------------------------------------------------------------

        dimension   rm(nr+1)
        dimension   zm(nz+1)
        dimension   eb(ne+1)
        dimension   ew(ne)
        dimension   tr(np,ne,nr,nz,2)
        real(8),allocatable :: vl_r(:),vl_z(:)

*-----------------------------------------------------------------------

        vl(ir,iz) = dble( 1 - itunt(m) / 7 )
     &            + dble( itunt(m) / 7 )
     &            * pi * ( rm(ir+1)**2 - rm(ir)**2 )
     &                 * ( zm(iz+1) - zm(iz) )

        pi = 2.d0 * asin(1.d0)

*-----------------------------------------------------------------------
! sumover
        allocate (vl_r(nz),vl_z(nr))
        vl_r(:) = 0.0d0
        vl_z(:) = 0.0d0
        do iz = 1, nz
          do ir = 1, nr
            vl_r(iz) = vl_r(iz) + vl(ir,iz)
            vl_z(ir) = vl_z(ir) + vl(ir,iz)
          enddo
        enddo
         

            if( itunt(m) .eq. 1 .or. itunt(m) .eq. 2 .or.
     &          itunt(m) .eq. 7 .or. itunt(m) .eq. 8 .or.
     &          itunt(m) .eq.13 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else
     &      if( itunt(m) .eq. 3 .or. itunt(m) .eq. 4 .or.
     &          itunt(m) .eq. 9 .or. itunt(m) .eq. 10.or.
     &          itunt(m) .eq.14 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

            else

               do i = 1, ne

                  ew(i) = 1.d0

               end do
               ew_sum = 1.0d0

            end if

*-----------------------------------------------------------------------

            if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

            do 100 iz = 1, nz
            do 100 ir = 1, nr
            do 100 ie = 1, ne
            do 100 ip = 1, np

              if(itunt(m).ge.13) then  ! T.Sato 2025/03/13

               call invert_stdev(m,A,B,
     &                         tr(ip,ie,ir,iz,1),
     &                         tr(ip,ie,ir,iz,2),
     &          ew(ie)/abs(rtfac(m)/facmax(m)))

              else

               call invert_stdev(m,A,B,
     &                         tr(ip,ie,ir,iz,1),
     &                         tr(ip,ie,ir,iz,2),
     &          vl(ir,iz)*ew(ie)/abs(rtfac(m)/facmax(m)))

              endif

               tr(ip,ie,ir,iz,1) = A
               tr(ip,ie,ir,iz,2) = B

! sumover
               fact_in = abs(rtfac(m)/facmax(m))
               call pletrz_sumover_stdev(1,m,ip,ie,ir,iz,
     &              fact_in,ew(ie),vl(ir,iz),ew_sum,vl_r(iz),vl_z(ir))

  100       continue
        deallocate (vl_r,vl_z)

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_letxyz(m,
     &                       np,nl,lt,nx,ny,nz,ne,xm,ym,zm,eb,tr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'
        include 'param01.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall21/ rtfac(itlmax)
        common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
        common /tall34/ itvm(itlmax)
        common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

        common /fact01/ facmax(itlmax)  ! kitamura23/03/31

*-----------------------------------------------------------------------

        dimension   lt(nl)
        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)
        dimension   eb(ne+1)
        dimension   ew(ne)
        dimension   tr(np,ne,nx*ny*nz,2)
        real(8),allocatable :: vl_x(:,:),vl_y(:,:),vl_z(:,:)

*-----------------------------------------------------------------------

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny

*-----------------------------------------------------------------------
*        set mesh volume ( unit = 1-6 ; vol = 1.0 )
*-----------------------------------------------------------------------

         vl(ix,iy,iz) = dble( 1 - itunt(m) / 7 )
     &                + dble( itunt(m) / 7 )
     &                * vls(nl,lt,itmcn(m),itvm(m),itmtr(m,4),
     &                      xm(ix),xm(ix+1),
     &                      ym(iy),ym(iy+1),
     &                      zm(iz),zm(iz+1))

*-----------------------------------------------------------------------
*           itunt(m) = 1 : cm/(keV/um)/source
*                    = 2 : MeV/(keV/um)/source
*                    = 3 : cm/ln(keV/um)/source
*                    = 4 : MeV/ln(keV/um)/source
*                    = 5 : cm/source
*                    = 6 : MeV/source
*                    = 7 : 1/cm^2/(keV/um)/source
*                    = 8 : MeV/cm^3/(keV/um)/source
*                    = 9 : 1/cm^2/ln(keV/um)/source
*                    =10 : MeV/cm^3/ln(keV/um)/source
*                    =11 : 1/cm^2/source
*                    =12 : MeV/cm^3/source
*-----------------------------------------------------------------------
! sumover
            allocate (vl_x(ny,nz),vl_y(nx,nz),vl_z(nx,ny))
            vl_x(:,:) = 0.0d0
            vl_y(:,:) = 0.0d0
            vl_z(:,:) = 0.0d0
            do iz= 1, nz
              do iy= 1, ny
                do ix= 1, nx
                  vl_x(iy,iz) = vl_x(iy,iz) + vl(ix,iy,iz)
                  vl_y(ix,iz) = vl_y(ix,iz) + vl(ix,iy,iz)
                  vl_z(ix,iy) = vl_z(ix,iy) + vl(ix,iy,iz)
                enddo
              enddo
            enddo

            if( itunt(m) .eq. 1 .or. itunt(m) .eq. 2 .or.
     &          itunt(m) .eq. 7 .or. itunt(m) .eq. 8 .or.
     &          itunt(m) .eq.13 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else
     &      if( itunt(m) .eq. 3 .or. itunt(m) .eq. 4 .or.
     &          itunt(m) .eq. 9 .or. itunt(m) .eq. 10.or.
     &          itunt(m) .eq.14 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

            else

               do i = 1, ne

                  ew(i) = 1.d0

               end do
               ew_sum = 1.0d0

            end if

*-----------------------------------------------------------------------

            if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

            do 100 iz = 1, nz
            do 100 iy = 1, ny
            do 100 ix = 1, nx
            do 100 ie = 1, ne
            do 100 ip = 1, np

              if(itunt(m).ge.13) then  ! T.Sato 2025/03/13

               call invert_stdev(m,A,B,
     &                         tr(ip,ie,icf(ix,iy,iz),1),
     &                         tr(ip,ie,icf(ix,iy,iz),2),
     &          ew(ie)/abs(rtfac(m)/facmax(m)))

              else

               call invert_stdev(m,A,B,
     &                         tr(ip,ie,icf(ix,iy,iz),1),
     &                         tr(ip,ie,icf(ix,iy,iz),2),
     &          vl(ix,iy,iz)*ew(ie)/abs(rtfac(m)/facmax(m)))

              endif

               tr(ip,ie,icf(ix,iy,iz),1) = A
               tr(ip,ie,icf(ix,iy,iz),2) = B

! sumover
               fact_in = abs(rtfac(m)/facmax(m))
               call pletxyz_sumover_stdev(1,m,ip,ie,ix,iy,iz,
     &                fact_in,ew(ie),vl(ix,iy,iz),
     &                ew_sum,vl_x(iy,iz),vl_y(ix,iz),vl_z(ix,iy))


  100       continue
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
      deallocate (vl_x,vl_y,vl_z)

      end subroutine


