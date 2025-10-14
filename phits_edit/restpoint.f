************************************************************************
*                                                                      *

      subroutine check_tpoint(m,iax,jsn,jsi,dsin,idsi,ill,ilf,ierr)
*                                                                      *
*       modified by K.Niita on 2015/11/30                              *
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************
      use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'param01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
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

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall21/ rtfac(itlmax)

      common /tall35/ iteps(itlmax)
      common /tall37/ itcnt(9,itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

!OBINATA: for resfile
      common /tall59/ irfll(itlmax), crfln(itlmax), itrff(itlmax)
      character crfln*100
      character irfile*100

      common /tall62/ itpon(itlmax), rtpon(itlmax,20,4)

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

      dimension lschn(35), ischn(35)
      character schan(35)*8

      data icsu / 22 /

      data ( schan(i), i = 1, 22 ) /
     &    'point   ','part    ','e-type  ','unit    ','axis    ',
     &    'file    ','title   ','factor  ','x-txt   ','y-txt   ',
     &    'epsout  ','ctmin(1)','ctmax(1)','ctmin(2)','ctmax(2)',
     &    'ctmin(3)','ctmax(3)','multipli','t-type  ','resfile ',
     &    'angel   ','ring    '/

      data ( lschn(i), i = 1, 22 ) /
     &     5,         4,         6,         4,         4,
     &     4,         5,         6,         5,         5,
     &     6,         8,         8,         8,         8,
     &     8,         8,         8,         6,         7,
     &     5,         4/

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

      character tname*9
      data      tname /'[t-point]'/

*-----------------------------------------------------------------------

      dimension jptyp(6), jpnkf(6)
      dimension kptyp(6,6), kpnkf(6,6)
      dimension ktln(6), ktli(6), ktls(6), kmst(6), knpat(6)
      dimension dkmax(6)
      dimension imst(6), kimst(6,6)
      dimension imtinf(4) ! frtati 2023/12/07

*-----------------------------------------------------------------------

      dimension stpon(20,4)
      dimension irsq(10)

*-----------------------------------------------------------------------

      character*10 cepn(100)
      dimension lepn(100)
      character*20 ctl(2)

*-----------------------------------------------------------------------

            jpn   = 0
            ierr  = 0

            iunt  = 1
            nrsq  = 0
            npont = -1
            itpnt = -1
            mtcn  = 0

            irsq(1) = 1
            irsq(2) = 2
            irsq(3) = 3
            irsq(4) = 4

         do i = 5, 10

            irsq(i) = 0

         end do

            iunt  = 1
            inpat = 0
            inaxi = 0
            infil = 0
            jmul  = 0

            langel = 0
            lxtxt  = 0
            lytxt  = 0
            ieps   = 0
            ittp   = 0
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

         do i = 1, icsu

            ischn(i) = 0

         end do

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( jpn  .ne. 0 ) goto 800
               if( jpn  .eq. 3 ) goto 800

               if( iskip .ne. 0 ) goto 140

  150 continue

               if( ierr .ne. 0 ) goto 800
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
     &          ipm .ne.  2 .and. ipm .ne. 18 .and.
     &          ipm .ne. 20 .and.
     &          ischn( ipm ) .gt. 1 ) goto 989

               ic = inumc(chlw,il+1,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 997

               icl = inumc(chlw,ic,i3,';') - 1

*-----------------------------------------------------------------------
*        point = number of point estimator, max is 20 in one tally
*-----------------------------------------------------------------------

         if( ipm .eq. 1 ) then

               if( itpnt .gt. 0 ) goto 965

               itpnt = 1

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               npont = nint( cvvv )

               if( npont .gt. 20 ) goto 967

*-----------------------------------------------------------------------

  141 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( jpn  .ne. 0 ) goto 800
               if( jpn  .eq. 3 ) goto 800

               if( iskip .ne. 0 ) goto 141

*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  102       if( ic .gt. i3 ) goto 101

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               irsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic) .eq. 'x' ) then

               irsq( mrsq + 1 ) = 1
               ic = ic +1

            else if( chlw(ic:ic) .eq. 'y' ) then

               irsq( mrsq + 1 ) = 2
               ic = ic + 1

            else if( chlw(ic:ic) .eq. 'z' ) then

               irsq( mrsq + 1 ) = 3
               ic = ic + 1

            else if( chlw(ic:ic+1) .eq. 'r0' ) then

               irsq( mrsq + 1 ) = 4
               ic = ic + 2

            else

               goto 101

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 102

  101       continue

            if( mrsq .gt. 0 ) then

                  ipox = 0
                  ipoy = 0
                  ipoz = 0
                  ipor = 0

               do k = 1, mrsq
                  if( irsq(k) .eq. 1 ) ipox = ipox + 1
                  if( irsq(k) .eq. 2 ) ipoy = ipoy + 1
                  if( irsq(k) .eq. 3 ) ipoz = ipoz + 1
                  if( irsq(k) .eq. 4 ) ipor = ipor + 1
               end do

                  if( ipox .ne. 1 .or. ipoy .ne. 1 .or.
     &                ipoz .ne. 1 .or. ipor .ne. 1 ) goto 998

                  nrsq = mrsq

                  goto 141

            else

                  nrsq = 4

            end if

         end if

*-----------------------------------------------------------------------

               mtcn = mtcn + 1

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( irsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( irsq(k) .gt. 0 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               stpon(mtcn,irsq(k)) = cvvv

            end if

         end do

            if( mtcn .lt. npont ) goto 141
            if( mtcn .eq. npont ) goto 140

*-----------------------------------------------------------------------
*        ring = number of ring estimator, max is 20 in one tally
*-----------------------------------------------------------------------

         else if( ipm .eq. 22 ) then

               if( itpnt .gt. 0 ) goto 965

               itpnt = 2

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               npont = nint( cvvv )

               if( npont .gt. 20 ) goto 967

*-----------------------------------------------------------------------

  142 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( jpn  .ne. 0 ) goto 800
               if( jpn  .eq. 3 ) goto 800

               if( iskip .ne. 0 ) goto 142

*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  112       if( ic .gt. i3 ) goto 111

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               irsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+3) .eq. 'axis' ) then

               irsq( mrsq + 1 ) = 1
               ic = ic + 4

            else if( chlw(ic:ic+1) .eq. 'ar' ) then

               irsq( mrsq + 1 ) = 2
               ic = ic + 2

            else if( chlw(ic:ic+1) .eq. 'rr' ) then

               irsq( mrsq + 1 ) = 3
               ic = ic + 2

            else if( chlw(ic:ic+1) .eq. 'r0' ) then

               irsq( mrsq + 1 ) = 4
               ic = ic + 2

            else

               goto 111

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 112

  111       continue

            if( mrsq .gt. 0 ) then

                  ipox = 0
                  ipoa = 0
                  ipod = 0
                  ipor = 0

               do k = 1, mrsq
                  if( irsq(k) .eq. 1 ) ipox = ipox + 1
                  if( irsq(k) .eq. 2 ) ipoa = ipoa + 1
                  if( irsq(k) .eq. 3 ) ipod = ipod + 1
                  if( irsq(k) .eq. 4 ) ipor = ipor + 1
               end do

                  if( ipox .ne. 1 .or. ipoa .ne. 1 .or.
     &                ipod .ne. 1 .or. ipor .ne. 1 ) goto 998

                  nrsq = mrsq

                  goto 142

            else

                  nrsq = 4

            end if

         end if

*-----------------------------------------------------------------------

               mtcn = mtcn + 1

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( irsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( irsq(k) .eq. 1 ) then

               if( chlw(ic:ic) .eq. 'x' ) then

                  stpon(mtcn,irsq(k)) = 1.d0

               else if( chlw(ic:ic) .eq. 'y' ) then

                  stpon(mtcn,irsq(k)) = 2.d0

               else if( chlw(ic:ic) .eq. 'z' ) then

                  stpon(mtcn,irsq(k)) = 3.d0

               else

                  goto 964

               end if

                  ic2 = ic + 1

            else if( irsq(k) .gt. 1 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               stpon(mtcn,irsq(k)) = cvvv

            end if

         end do

            if( mtcn .lt. npont ) goto 142
            if( mtcn .eq. npont ) goto 140

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

               if( istyp .ne. 2 .and. istyp .ne. 14 ) goto 995
               if( inpat .gt. 2 ) goto 995

                  iptyp(inpat) = istyp
                  ipnkf(inpat) = inkf0
                  ipsub(inpat) = isubt  ! kitamura22/03/31

               goto 400

*-----------------------------------------------------------------------
*        energy mesh
*-----------------------------------------------------------------------

         else if( ipm .eq. 3 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

               ietp = nint( cvvv )

               call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     'e',ietp,ine,emin,emax,edel,isteg)

               if( ierr .ne. 0 ) return

               goto 150

*-----------------------------------------------------------------------
*        time mesh
*-----------------------------------------------------------------------

         else if( ipm .eq. 19 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

               ittp = nint( cvvv )

               call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     't',ittp,int,tmin,tmax,tdel,isttg)

               if( ierr .ne. 0 ) return

               goto 150

*-----------------------------------------------------------------------
*        unit
*-----------------------------------------------------------------------

         else if( ipm .eq. 4 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               iunt = nint( cvvv )

               if( iunt .lt. 1 .or.
     &           ( iunt .gt. 3 .and. iunt .lt. 11 ) .or.
     &             iunt .gt. 13 ) goto 988

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        axis
*-----------------------------------------------------------------------

         else if( ipm .eq. 5 ) then

  600       inaxi = inaxi + 1

            if( inaxi .gt. 2 ) goto 991

            if( chlw(ic:ic+2) .eq. 'eng' ) then

               iaxis(inaxi) = 1
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic) .eq. 't' ) then

               iaxis(inaxi) = 2
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

               if( infil .gt. 2 ) goto 990

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

         else if( ipm .eq. 20 ) then

               icf = min( inumc(chlw,ic,icl,' ') - 1, icl )

               lrfile = icf-ic+1

               irfile(1:lrfile) = chin(ic:icf)

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

         else if( ipm .eq. 21 ) then

               ict = min( ic + 199, i2 )

               angelp = chin(ic:ict)

               langel = ict - ic + 1

*-----------------------------------------------------------------------
*        x-txt
*-----------------------------------------------------------------------

         else if( ipm .eq. 9 ) then

               ict = min( ic + 199, i2 )

               cxtxt = chin(ic:ict)

               lxtxt = ict - ic + 1

*-----------------------------------------------------------------------
*        y-txt
*-----------------------------------------------------------------------

         else if( ipm .eq. 10 ) then

               ict = min( ic + 199, i2 )

               cytxt = chin(ic:ict)

               lytxt = ict - ic + 1

*-----------------------------------------------------------------------
*        factor
*-----------------------------------------------------------------------

         else if( ipm .eq. 8 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               rfact = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        epsout
*-----------------------------------------------------------------------

         else if( ipm .eq. 11 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ieps = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        counter
*-----------------------------------------------------------------------

         else if( ipm .ge. 12 .and. ipm .le. 17 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               icount((ipm-20+2)/2) = 1
               icount( ipm-20+4)    = nint( cvvv )

               if( abs( nint(cvvv) ) .gt. 9999 ) goto 978

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        multiplier
*-----------------------------------------------------------------------

         else if( ipm .eq. 18 ) then

               jmul = jmul + 1

               if( jmul .gt. 6 ) goto 970

            if( chlw(ic:ic+2) .eq. 'all' ) then

               jtln = 1
               jtal = 1

            else

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               jtln = nint( cvvv )
               jtal = 0

               if( jtln .lt. 0 ) goto 971

            end if

cfrtati 2023/12/07 imtinf added
               call tmultipl(jsn,jsi,dsin,idsi,ill,ilf,
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                      jtln,jtal,jnpat,jptyp,jpnkf,dlmax,
     &                      jtli,jtls,mmst,imst,imtinf)

                  if( ierr .ne. 0 ) return

                  ktln(jmul)  = jtln
                  ktli(jmul)  = jtli
                  ktls(jmul)  = jtls
                  kmst(jmul)  = mmst
                  dkmax(jmul) = dlmax

                  knpat(jmul) = jnpat

               do j = 1, mmst

                  kimst(jmul,j) = imst(j)

               end do

               do j = 1, jnpat

                  kptyp(jmul,j) = jptyp(j)
                  kpnkf(jmul,j) = jpnkf(j)

               end do

                  if( jpn  .eq. 3 ) goto 800

               goto 150

*-----------------------------------------------------------------------

         end if

            goto 140

  800 continue

*-----------------------------------------------------------------------
*     check
*-----------------------------------------------------------------------

      iec = 0

      call check_point(m,npont,itpnt,stpon,iec,cepn,lepn,ierr)

      call check_part(m,iptyp,ipnkf,inpat,iec,cepn,lepn,ierr)

      call check_type('e',iec,cepn,lepn,
     &                itety(m),rtema(m),rtemi(m),itenm(m),
     &                ietp, emax, emin, ine)

      call check_type('t',iec,cepn,lepn,
     &                ittty(m),rttma(m),rttmi(m),ittnm(m),
     &                ittp, tmax, tmin, int)

      call check_unit(m,iunt,iec,cepn,lepn,ierr)

      call check_factor(m,rfact,iec,cepn,lepn,ierr)

      call check_axis(m,iaxis,inaxi,iec,cepn,lepn,ierr)

*-----------------------------------------------------------------------

      call judge_tall_check(iec,cepn,lepn,dsin(jsn),idsi(jsn),ierr)

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  964    m_err = 'Axis should be x, y, z'//tname
         ErrCha = ''
         ErrID = 'L:992/R:check_tpoint/F:restpoint.f'
         goto 999

  965    m_err = 'only one point or ring in a tally.'//tname
         ErrCha = ''
         ErrID = 'L:997/R:check_tpoint/F:restpoint.f'
         goto 999

  967    m_err = 'max point number is 20 in one tall section.'//tname
         ErrCha = ''
         ErrID = 'L:1002/R:check_tpoint/F:restpoint.f'
         goto 999

  968    m_err = 'mset number is inconsistent.'//tname
         ErrCha = ''
         ErrID = 'L:1007/R:check_tpoint/F:restpoint.f'
         goto 999

  969    m_err = 'number of mset should be the same.'//tname
         ErrCha = ''
         ErrID = 'L:1012/R:check_tpoint/F:restpoint.f'
         goto 999

  970    m_err = 'def of multiplier should be less than 7 '//tname
         ErrCha = ''
         ErrID = 'L:1017/R:check_tpoint/F:restpoint.f'
         goto 999

  971    m_err = 'number of multiplier is negative'//tname
         ErrCha = ''
         ErrID = 'L:1022/R:check_tpoint/F:restpoint.f'
         goto 999

  978    m_err = 'counter should be from -9999 to 9999'
         ErrCha = ''
         ErrID = 'L:1027/R:check_tpoint/F:restpoint.f'
         goto 999

  985    m_err = 'Unit is Lethargy but energy mesh points are negative'
         ErrCha = ''
         ErrID = 'L:1032/R:check_tpoint/F:restpoint.f'
         goto 999

  986    write(dkam,'(i5)') itlmax
         m_err = 'Total tally number exceeds itlmax = '//dkam
         ErrCha = ''
         ErrID = 'L:1038/R:check_tpoint/F:restpoint.f'
         goto 999

  987    m_err = 'Unknown parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1043/R:check_tpoint/F:restpoint.f'
         goto 999

  988    m_err = 'Unit should be 1,2,3,11,12,13 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1048/R:check_tpoint/F:restpoint.f'
         goto 999

  989    m_err = 'Double definition of the parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1053/R:check_tpoint/F:restpoint.f'
         goto 999

  990    m_err = 'Too many file name in tally '//tname
         ErrCha = ''
         ErrID = 'L:1058/R:check_tpoint/F:restpoint.f'
         goto 999

  991    m_err = 'Too many axis in tally '//tname
         ErrCha = ''
         ErrID = 'L:1063/R:check_tpoint/F:restpoint.f'
         goto 999

  992    m_err = 'Unknown axis name in tally '//tname
         ErrCha = ''
         ErrID = 'L:1068/R:check_tpoint/F:restpoint.f'
         goto 999

  993    m_err = 'In this line, [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:1073/R:check_tpoint/F:restpoint.f'
         goto 999

  994    m_err = 'Name of particle is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1078/R:check_tpoint/F:restpoint.f'
         goto 999

  995    m_err = 'Particles is only neutron, photon in tally '//tname
         ErrCha = ''
         ErrID = 'L:1083/R:check_tpoint/F:restpoint.f'
         goto 999

  997    m_err = 'Description of parameter is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1088/R:check_tpoint/F:restpoint.f'
         goto 999

  998    m_err = 'Unknown point parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1093/R:check_tpoint/F:restpoint.f'
         goto 999

*-----------------------------------------------------------------------

  999 continue

         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1

         write(ErrCha,*) 'Error: ' // m_err
         ErrID = 'L:1105/R:check_tpoint/F:restpoint.f' !E53_001_001
         call ErrWrite(ErrID,ErrCha)

*-----------------------------------------------------------------------

      return

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_tpoint(m,iax,ierr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        use RESTALMOD
        use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05
      use moddas_mesh

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        common /mpi00/ npe, me

*-----------------------------------------------------------------------

        common /talmm/  nmmax, lmmax, itlmx

*-----------------------------------------------------------------------

        common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------

        common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                  rtemi(itlmax), rtema(itlmax), rtedl(itlmax)
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
        common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                  rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
        common /tall21/ rtfac(itlmax)
        common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
        common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)
        common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

        common /tall59/ irfll(itlmax), crfln(itlmax), itrff(itlmax)
        character crfln*100

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
        ioe  = 1

*-----------------------------------------------------------------------
*   open restart file
*-----------------------------------------------------------------------

        call open_resfile(m,noe,jsn,jsi,dsin,idsi,ill,ilf,newtall,ierr)

        if ( newtall .ne. 0 ) goto 900  !! it's new tally
        if ( ierr    .ne. 0 ) goto 900

*-----------------------------------------------------------------------
*   check tally
*-----------------------------------------------------------------------

        if(ireschk.eq.0) then  ! T.Sato 2013/10/19
        call check_tpoint(m,iax,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    ierr)
        endif

        if (ierr .ne. 0 ) goto 900

*-----------------------------------------------------------------------

        rewind(jsi(ioe))
        ill(jsn(ioe),ioe) = 1

*-----------------------------------------------------------------------

        idas1 = lmmax
        idas2 = idas1 + itenm(m)
        idas3 = idas2 + ittnm(m)

*-----------------------------------------------------------------------
*   read tally
*-----------------------------------------------------------------------

          call read_point(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itmsh(m),itenm(m),itmst(m),ittnm(m),
     &                    das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))


*-----------------------------------------------------------------------
*   restore tally
*-----------------------------------------------------------------------

          call restore_point(m,
     &                    itpan(m),itmsh(m),itenm(m),itmst(m),ittnm(m),
     &                    das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

  900   continue
*-----------------------------------------------------------------------
*   close restart file
*-----------------------------------------------------------------------

          close(jsi(ioe))

*-----------------------------------------------------------------------
  999   continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_point(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                      np,nr,ne,nm,nt,eb,tb,tr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)

        common /cparm/  maxbch,maxcas
        common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
        common /restart/ resc2(itlmax), resc3(itlmax)

        common /tall21/ rtfac(itlmax)
        common /fact01/ facmax(itlmax)
        common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------

      dimension   eb(ne+1)
      dimension   ew(ne)
      dimension   tb(nt+1)
      dimension   tw(nt)
      dimension   tr(np,ne,nt,nr,nm,2)
! sumover
      dimension   tott_sum(1,2)
      data nsame /1/
*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        character chin*200, chlw*200, chcm*200
        character dc2*15

*-----------------------------------------------------------------------
*        energy axis
*-----------------------------------------------------------------------

        if( itaxs(m,iax) .eq. 1 ) then
          !'#e-lowere-upper'

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nr_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petrm.inc'

          nmstepi = 1
          do 190 ir = 1, nr, nrstepi
          do 190 im = 1, nm, nmstepi
          do 190 it = 1, nt, ntstepi
          do 190 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#e-lowere-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nestepi = 1
            do ie = 1, ne, nestepi



              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &                  dmm0, dmm1,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &      ir+irloop-1,im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),imloop=1,nmstepi)

            end do

! sumover
            ie = 1
            call ptpointp_sumover_getput(jsi,m,iax,
     &           npstepi,nestepi,ntstepi,nrstepi,nmstepi,
     &           ip,ie,it,ir,im,nsame,tott_sum)

  190     continue

*-----------------------------------------------------------------------
*        time axis
*-----------------------------------------------------------------------

        else if( itaxs(m,iax) .eq. 11 ) then
          !'#t-lowert-upper'

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nr_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petrm.inc'

          nmstepi = 1
          do 490 ir = 1, nr, nrstepi
          do 490 im = 1, nm, nmstepi
          do 490 ie = 1, ne, nestpei
          do 490 ip = 1, np, npstepi
            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                       '#t-lowert-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            ntstepi = 1
            do it = 1, nt, ntstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &                  dmm0, dmm1,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &      ir+irloop-1,im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),imloop=1,nmstepi)

            end do

! sumover
            it = 1
            call ptpointp_sumover_getput(jsi,m,iax,
     &           npstepi,nestepi,ntstepi,nrstepi,nmstepi,
     &           ip,ie,it,ir,im,nsame,tott_sum)

  490     continue

*-----------------------------------------------------------------------
        end if

  800 continue

*-----------------------------------------------------------------------

  900 continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_point(m,
     &                         np,nr,ne,nm,nt,eb,tb,tr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall21/ rtfac(itlmax)

        common /cparm/  maxbch,maxcas
        common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
        common /restart/ resc2(itlmax), resc3(itlmax)

        common /fact01/ facmax(itlmax)  ! kitamura23/03/31

*-----------------------------------------------------------------------

      dimension   eb(ne+1)
      dimension   ew(ne)
      dimension   tb(nt+1)
      dimension   tw(nt)
      dimension   tr(np,ne,nt,nr,nm,2)

*-----------------------------------------------------------------------
*           itunt(m) = 1,  11   : /cm^2/source
*                                 /cm^2/sec/source
*                    = 2   12   : /cm^2/MeV/source
*                                 /cm^2/MeV/nsec/source
*                    = 3   13   : /cm^2/Lethargy/source
*                                 /cm^2/Lethargy/nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  1 .or. itunt(m) .eq. 11 ) then

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
*           itunt(m) = 11, 12, 13 : /.../nsec/source
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

            if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

            do 100 im = 1, nm
            do 100 ir = 1, nr
            do 100 it = 1, nt
            do 100 ie = 1, ne
            do 100 ip = 1, np

                  call invert_stdev(m,A,B,
     &                             tr(ip,ie,it,ir,im,1),
     &                             tr(ip,ie,it,ir,im,2),
     &              ew(ie)*tw(it)/abs(rtfac(m)/facmax(m)))

                  tr(ip,ie,it,ir,im,1) = A
                  tr(ip,ie,it,ir,im,2) = B

! sumover
                 fact_in = abs(rtfac(m)/facmax(m))
                 call ptpointp_sumover_stdev(1,m,ip,ie,it,ir,im,
     &                fact_in,ew(ie),tw(it),ew_sum,tw_sum)


  100       continue

*-----------------------------------------------------------------------

      end subroutine
