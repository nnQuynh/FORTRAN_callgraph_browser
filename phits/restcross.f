
************************************************************************
*                                                                      *
      subroutine check_tcross(m,iax,jsn,jsi,dsin,idsi,ill,ilf,ierr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************
      use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'param01.inc' ! frtati 2021/10/05
      include 'err.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)

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

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /tall08/ rtrx0(itlmax), rtry0(itlmax)

      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall13/ itrcn(itlmax), itrcr(itlmax), itrca(itlmax)
      common /tall26/ itrcs(itlmax), itrcc(itlmax), itrss(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall24/ itrcm(itlmax), itrcg(itlmax)
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)
      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall19/ itsmn(itlmax), itstm(itlmax)
      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall21/ rtfac(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall31/ iterl(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall37/ itcnt(9,itlmax)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)

*-----------------------------------------------------------------------
      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

*-----------------------------------------------------------------------

      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)
      common /tall49/ itglt(itlmax)

!OBINATA: for resfile
      common /tall59/ irfll(itlmax), crfln(itlmax), itrff(itlmax)
      character crfln*100
      character irfile*100

      common /tall60/ itger(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /tall50/ itlmt(itlmax), ite2l(itlmax)   ! CCSE 2022.08.31

      common/nzerocom/nzerocheck ! T.Sato 2018/3/6, allow nz = 0 for [t-cross]

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
CCSE added for LET parameter, 37 -> 39 (2022.08.31) >>>>>

      dimension lschn(39), ischn(39)
      character schan(39)*8

      data icsu / 39 /

      data ( schan(i), i = 1, 39 ) /
     &    'mesh    ','part    ','e-type  ','unit    ','axis    ',
     &    'file    ','title   ','angel   ','output  ','2d-type ',
     &    'factor  ','a-type  ','x-txt   ','y-txt   ','z-txt   ',
     &    'gshow   ','epsout  ','ctmin(1)','ctmax(1)','ctmin(2)',
     &    'ctmax(2)','ctmin(3)','ctmax(3)','resol   ','width   ',
     &    't-type  ','trcl    ','*trcl   ','dump    ','gslat   ',
     &    'resfile ','ginfo   ','multipli','stdcut  ','sangel  ',
     &    'enclos  ','iangform','letmat  ','eng2let '/

      data ( lschn(i), i = 1, 39 ) /
     &     4,         4,         6,         4,         4,
     &     4,         5,         5,         6,         7,
     &     6,         6,         5,         5,         5,
     &     5,         6,         8,         8,         8,
     &     8,         8,         8,         5,         5,
     &     6,         4,         5,         4,         5,
     &     7,         5,         8,         6,         6,
     &     6,         8,         6,         7/

*-----------------------------------------------------------------------
cfrtati 2021/10/05 6 -> mxpart
      dimension iptyp(mxpart), ipnkf(mxpart)
      dimension imtyp(mxpart,mxpart), imnkf(mxpart,mxpart)
      dimension jstyp(mxpart), jnkf0(mxpart)

      common /subtra/ isubt, ipsub(mxpart)  ! kitamura22/03/31

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
      data      tname /'[t-cross]'/

      dimension icount(9)
      dimension vtrs(13)

      dimension jdump(30)

*-----------------------------------------------------------------------

      dimension jptyp(6), jpnkf(6)
      dimension kptyp(6,6), kpnkf(6,6)
      dimension ktln(6), ktli(6), ktls(6), kmst(6), knpat(6)
      dimension dkmax(6)
      dimension imst(6), kimst(6,6)
      dimension imtinf(4) ! frtati 2023/12/07

*-----------------------------------------------------------------------


      character*10 cepn(100)
      dimension lepn(100)
      character*20 ctl(2)

*-----------------------------------------------------------------------

            ierr  = 0
            jpn   = 0

            iunt  = 1
            imate = 0
            jmate = 1
            iout  = 1
            inpat = 0
            inaxi = 0
            infil = 0
            idtyp = 3   ! S.Abe 2019/03/29, 2 -> 3
            jmul  = 0

            langel = 0
            lxtxt  = 0
            lytxt  = 0
            lztxt  = 0
            lgshow = 0
            ieps   = 0
            ietp   = 0
            ittp   = 0
            iatp   = 0
            jatp   = 0
            igkst  = 0
            idtt   = 0
            ktrs   = 0
            idump  = 0
            igslt  = 1
            letmat = 0  ! CCSE 2022/08/31
            ieng2let= 0 ! CCSE 2022/08/31 Flag to convert energy to LET
            lrfile = 0 !OBINATA
            irfflg = 0 !OBINATA
            ntrcn  = 0 !OBINATA

            amin = 0.0
            amax = 0.0
            adel = 0.0 !FURUTA

            icount(1) = 0
            icount(2) = 0
            icount(3) = 0
            icount(4) = -9999
            icount(5) =  9999
            icount(6) = -9999
            icount(7) =  9999
            icount(8) = -9999
            icount(9) =  9999

            ireso = 1
            width = 0.5

            rfact = 1.0

         do i = 1, icsu

            ischn(i) = 0

         end do

         do i = 1, 30
            jdump(i) = 0
         end do

            stdcut = -1.0

            ienclo = 0   ! S.Abe 2018/10/26

            iangform = 0   ! S.Abe 2019/11/26

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

            if( ipm .ne. 5 .and. ipm .ne. 6 .and.
     &          ipm .ne. 2 .and. ipm .ne. 33 .and.
     &          ischn( ipm ) .gt. 1 ) goto 989

               ic = inumc(chlw,il+1,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 997

               icl = inumc(chlw,ic,i3,';') - 1

*-----------------------------------------------------------------------
*        mesh = region,  r-z, xyz
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

               !! serach line 'reg='
  160          call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &              jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                  if( jpn  .eq. 3 ) goto 800
                  if( iskip .ne. 0 ) goto 160

                  icr = inumc(chlw,1,i3,'=') + 1
                  call onum(chlw,icr,i3,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 800

                  ntrcn  = nint( cvvv )

                  do j = 1, ntrcn + 1
                    call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
                    if( ierr .ne. 0 ) goto 800
                    if( jpn  .eq. 3 ) goto 800
                  end do
                  call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

            else if( imesh .eq. 2 ) then
               nzerocheck=1  ! allow nz = 0, T.Sato 2018/03/06
               call trzmesh(jsn,jsi,dsin,idsi,ill,ilf,
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                      rzx0,rzy0,
     &                      irtp,inr,rmin,rmax,rdel,istrg,
     &                      iztp,inz,zmin,zmax,zdel,istzg)
               nzerocheck=0  ! allow nz = 0, T.Sato 2018/03/06
                  if( jpn  .eq. 3 ) goto 800
                  if( irtp .lt. 0 ) goto 996
                  if( iztp .lt. 0 ) goto 996

            else if( imesh .eq. 3 ) then
               nzerocheck=1  ! allow nz = 0, T.Sato 2018/03/06
               call txymesh(jsn,jsi,dsin,idsi,ill,ilf,
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                      ixtp,inx,xmin,xmax,xdel,istxg,
     &                      iytp,iny,ymin,ymax,ydel,istyg,
     &                      iztp,inz,zmin,zmax,zdel,istzg)
               nzerocheck=0  ! allow nz = 0, T.Sato 2018/03/06

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


               goto 150

*-----------------------------------------------------------------------
*        angle mesh
*-----------------------------------------------------------------------

         else if( ipm .eq. 12 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

               iatp = nint( cvvv )
               jatp = abs( iatp )

               call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     'a',jatp,ina,amin,amax,adel,istag)


               goto 150

*-----------------------------------------------------------------------
*        time mesh
*-----------------------------------------------------------------------

         else if( ipm .eq. 26 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

               ittp = nint( cvvv )

               call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     't',ittp,intt,tmin,tmax,tdel,isttg)


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
     &           ( iunt .gt. 6 .and. iunt .lt. 11 ) .or.
     &             iunt .gt. 16 ) goto 988

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        axis
*-----------------------------------------------------------------------

         else if( ipm .eq. 5 ) then

  600       inaxi = inaxi + 1

            if( inaxi .gt. 6 ) goto 991

            if( chlw(ic:ic+2) .eq. 'eng' ) then

               iaxis(inaxi) = 1
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+2) .eq. 'let' ) then

               iaxis(inaxi) = 14
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               iaxis(inaxi) = 2
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+2) .eq. 'the' ) then

               iaxis(inaxi) = 10
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+2) .eq. 'cos' ) then

               iaxis(inaxi) = 8
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+1) .eq. 'xy' .or.
     &               chlw(ic:ic+1) .eq. 'yx' ) then

               iaxis(inaxi) = 7
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic+1) .eq. 'yz' .or.
     &               chlw(ic:ic+1) .eq. 'zy' ) then

               iaxis(inaxi) = 11
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic+1) .eq. 'zx' .or.
     &               chlw(ic:ic+1) .eq. 'xz' ) then

               iaxis(inaxi) = 12
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic+1) .eq. 'rz' .or.
     &               chlw(ic:ic+1) .eq. 'zr' ) then

               iaxis(inaxi) = 13
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

            else if( chlw(ic:ic) .eq. 't' ) then

               iaxis(inaxi) = 9
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
         else if( ipm .eq. 31 ) then

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

!-----------------------------------------------------------------------
!        sangel parameters
!-----------------------------------------------------------------------

         else if( ipm .eq. 35 ) then

            itsanf = itsanf + 1

            call read_sangel(jsn,jsi,dsin,idsi,ill,ilf,
     &                       jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr
     &                      ,ic ,icl ,tname)
            if( ierr .ne. 0 ) return

*-----------------------------------------------------------------------
*        x-txt
*-----------------------------------------------------------------------

         else if( ipm .eq. 13 ) then

               ict = min( ic + 199, i2 )

               cxtxt = chin(ic:ict)

               lxtxt = ict - ic + 1

*-----------------------------------------------------------------------
*        y-txt
*-----------------------------------------------------------------------

         else if( ipm .eq. 14 ) then

               ict = min( ic + 199, i2 )

               cytxt = chin(ic:ict)

               lytxt = ict - ic + 1

*-----------------------------------------------------------------------
*        z-txt
*-----------------------------------------------------------------------

         else if( ipm .eq. 15 ) then

               ict = min( ic + 199, i2 )

               cztxt = chin(ic:ict)

               lztxt = ict - ic + 1

*-----------------------------------------------------------------------
*        gshow
*-----------------------------------------------------------------------

         else if( ipm .eq. 16 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               lgshow = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        gslat
*-----------------------------------------------------------------------

         else if( ipm .eq. 30 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               igslt = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        epsout
*-----------------------------------------------------------------------

         else if( ipm .eq. 17 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ieps = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        resolution
*-----------------------------------------------------------------------

         else if( ipm .eq. 24 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ireso = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        width
*-----------------------------------------------------------------------

         else if( ipm .eq. 25 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               width = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        counter
*-----------------------------------------------------------------------

         else if( ipm .ge. 18 .and. ipm .le. 23 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               icount((ipm-18+2)/2) = 1
               icount( ipm-18+4)    = nint( cvvv )

               if( abs( nint(cvvv) ) .gt. 9999 ) goto 978

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        output
*-----------------------------------------------------------------------

         else if( ipm .eq. 9 ) then

            if( chlw(ic:ic+3) .eq. 'flux' ) then

               iout = 1

            else if( chlw(ic:ic+6) .eq. 'current' ) then

               iout = 2

            else if( chlw(ic:ic+5) .eq. 'f-curr' ) then

               iout = 3

            else if( chlw(ic:ic+5) .eq. 'b-curr' ) then

               iout = 4

            else if( chlw(ic:ic+5) .eq. 'o-curr' ) then

               iout = 5

            else if( chlw(ic:ic+6) .eq. 'of-curr' ) then

               iout = 6

            else if( chlw(ic:ic+6) .eq. 'ob-curr' ) then

               iout = 7

            else if( chlw(ic:ic+5) .eq. 'a-curr' ) then

               iout = 8

            else if( chlw(ic:ic+6) .eq. 'oa-curr' ) then

               iout = 9

            else if( chlw(ic:ic+5) .eq. 'a-flux' ) then

               iout = 10

            else if( chlw(ic:ic+6) .eq. 'oa-flux' ) then

               iout = 11

            else

               goto 984

            end if

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        2d-type
*-----------------------------------------------------------------------

         else if( ipm .eq. 10 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               idtyp = nint( cvvv )

               if( idtyp .lt. 1 .or. idtyp .gt. 7 ) goto 983

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        factor
*-----------------------------------------------------------------------

         else if( ipm .eq. 11 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               rfact = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        dump
*-----------------------------------------------------------------------

         else if( ipm .eq. 29 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               idump = nint( cvvv )

            if( idump .ne. 0 ) then

  151                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                        if( jpn  .eq. 3 ) goto 979

                        if( iskip .ne. 0 ) goto 151

                  ic = i1

               do k = 1, abs( idump )

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

                     jdump(k) = nint( cvvv )

                     if( jdump(k) .gt. 20 .or.
     &                   jdump(k) .le.  0 ) goto 977

                     ic = ic2

               end do

            end if

*-----------------------------------------------------------------------
*        infog
*-----------------------------------------------------------------------

         else if( ipm .eq. 32 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               infog = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        stdcut
*-----------------------------------------------------------------------

         else if( ipm .eq. 34 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               stdcut = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        multiplier
*-----------------------------------------------------------------------

         else if( ipm .eq. 33 ) then

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
*        Parameter for detection area setting
*-----------------------------------------------------------------------

         else if( ipm .eq. 36 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ienclo = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        Parameter for the angle formed from
*-----------------------------------------------------------------------

         else if( ipm .eq. 37 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               iangform = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

CCSE added for LET parameter, 37 -> 39 (2022.08.31) >>>>>
*-----------------------------------------------------------------------
*        letmat
*-----------------------------------------------------------------------

         else if( ipm .eq. 38 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               letmat = nint( cvvv )


               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        eng2let
*-----------------------------------------------------------------------

         else if( ipm .eq. 39 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ieng2let = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------

         end if

            goto 140

  800 continue

*-----------------------------------------------------------------------
*     check
*-----------------------------------------------------------------------

      iec = 0

      call check_mesh(m,imesh,iec,cepn,lepn,ierr)

      if ( itmsh(m) .eq. 1 ) then

        if ( itrcn(m) .ne. ntrcn ) then
          iec = iec + 1
          cepn(iec) = 'reg'
          lepn(iec) = 3
        end if

        call check_type('e',iec,cepn,lepn,
     &                  itety(m),rtema(m),rtemi(m),itenm(m),
     &                  ietp, emax, emin, ine)

      else if ( itmsh(m) .eq. 2 ) then

        call check_x0y0(m,rzx0,rzy0,iec,cepn,lepn,ierr)

        call check_type('r',iec,cepn,lepn,
     &                  itrty(m),rtrma(m),rtrmi(m),itrnm(m),
     &                  irtp, rmax, rmin, inr)

        call check_type('z',iec,cepn,lepn,
     &                  itzty(m),rtzma(m),rtzmi(m),itznm(m),
     &                  iztp, zmax, zmin, inz)
          call check_type('e',iec,cepn,lepn,
     &                    itety(m),rtema(m),rtemi(m),itenm(m),
     &                    ietp, emax, emin, ine)

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

        call check_type('e',iec,cepn,lepn,
     &                  itety(m),rtema(m),rtemi(m),itenm(m),
     &                  ietp, emax, emin, ine)

      end if

      call check_unit(m,iunt,iec,cepn,lepn,ierr)

      call check_axis(m,iaxis,inaxi,iec,cepn,lepn,ierr)

      if ( itaxs(m,iax) .eq. 9 ) then !! time axis
        call check_type('t',iec,cepn,lepn,
     &                  ittty(m),rttma(m),rttmi(m),ittnm(m),
     &                  ittp, tmax, tmin, intt)
      end if

      call check_2dtype(m,idtyp,iec,cepn,lepn,ierr)

      call check_factor(m,rfact,iec,cepn,lepn,ierr)

      call check_part(m,iptyp,ipnkf,inpat,iec,cepn,lepn,ierr)

      call check_output(m,iout,iec,cepn,lepn,ierr)

      call check_enclos(m,ienclo,iec,cepn,lepn,ierr)

      call check_angform(m,iangform,iec,cepn,lepn,ierr)

*-----------------------------------------------------------------------

      call judge_tall_check(iec,cepn,lepn,dsin(jsn),idsi(jsn),ierr)

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  970    m_err = 'def of multiplier should be less than 7 '//tname
         ErrCha = ''
         ErrID = 'L:1284/R:check_tcross/F:restcross.f'
         goto 999

  971    m_err = 'number of multiplier is negative'//tname
         ErrCha = ''
         ErrID = 'L:1289/R:check_tcross/F:restcross.f'
         goto 999

  968    m_err = 'mset number is inconsistent.'//tname
         ErrCha = ''
         ErrID = 'L:1294/R:check_tcross/F:restcross.f'
         goto 999

  969    m_err = 'number of mset should be the same.'//tname
         ErrCha = ''
         ErrID = 'L:1299/R:check_tcross/F:restcross.f'
         goto 999


  976    m_err = 'dump is available only mesh = reg.'
         ErrCha = ''
         ErrID = 'L:1305/R:check_tcross/F:restcross.f'
         goto 999

  977    m_err = 'dump id should be 1 - 20.'
         ErrCha = ''
         ErrID = 'L:1310/R:check_tcross/F:restcross.f'
         goto 999

  978    m_err = 'counter should be from -9999 to 9999'
         ErrCha = ''
         ErrID = 'L:1315/R:check_tcross/F:restcross.f'
         goto 999

  979    m_err = 'Description of dump parameter is wrong in tally '
     &            //tname
         ErrCha = ''
         ErrID = 'L:1321/R:check_tcross/F:restcross.f'
         goto 999

  983    m_err = '2D-type should be 1-7 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1326/R:check_tcross/F:restcross.f'
         goto 999

  984    m_err = 'Unknown output parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1331/R:check_tcross/F:restcross.f'
         goto 999

  985    m_err = 'Unit is Lethargy but energy mesh points are negative'
         ErrCha = ''
         ErrID = 'L:1336/R:check_tcross/F:restcross.f'
         goto 999

  986    write(dkam,'(i5)') itlmax
         m_err = 'Total tally number exceeds itlmax = '//dkam
         ErrCha = ''
         ErrID = 'L:1342/R:check_tcross/F:restcross.f'
         goto 999

  987    m_err = 'Unknown parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1347/R:check_tcross/F:restcross.f'
         goto 999

  988    m_err = 'Unit should be 1, 2, 3 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1352/R:check_tcross/F:restcross.f'
         goto 999

  989    m_err = 'Double definition of the parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1357/R:check_tcross/F:restcross.f'
         goto 999

  990    m_err = 'Too many file name in tally '//tname
         ErrCha = ''
         ErrID = 'L:1362/R:check_tcross/F:restcross.f'
         goto 999

  991    m_err = 'Too many axis in tally '//tname
         ErrCha = ''
         ErrID = 'L:1367/R:check_tcross/F:restcross.f'
         goto 999

  992    m_err = 'Unknown axis name in tally '//tname
         ErrCha = ''
         ErrID = 'L:1372/R:check_tcross/F:restcross.f'
         goto 999

  993    m_err = 'In this line, [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:1377/R:check_tcross/F:restcross.f'
         goto 999

  994    m_err = 'Name of particle is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1382/R:check_tcross/F:restcross.f'
         goto 999

  995    m_err = 'Number of particles is larger than 6 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1387/R:check_tcross/F:restcross.f'
         goto 999

  996    m_err = 'Description of mesh is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1392/R:check_tcross/F:restcross.f'
         goto 999

  997    m_err = 'Description of parameter is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1397/R:check_tcross/F:restcross.f'
         goto 999

  998    m_err = 'Unknown mesh parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1402/R:check_tcross/F:restcross.f'
         goto 999

*-----------------------------------------------------------------------

  999 continue

         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1

         write(ErrCha,*) 'Error: ' // m_err
         ErrID = 'L:1414/R:check_tcross/F:restcross.f' !E52_001_001
         call ErrWrite(ErrID,ErrCha)

*-----------------------------------------------------------------------

      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine read_tcross(m,iax,ierr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        use RESTALMOD
        use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region

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
        common /tall13/ itrcn(itlmax), itrcr(itlmax), itrca(itlmax)
        common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                  rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
        common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                  rtami(itlmax), rtama(itlmax), rtadl(itlmax)
        common /tall26/ itrcs(itlmax), itrcc(itlmax), itrss(itlmax)

*-----------------------------------------------------------------------
      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /tall73/ itenclo(itlmax), itangform(itlmax)

*-----------------------------------------------------------------------

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

        if ( any( itaxs(m,iax) .eq. (/ 7 /) )
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

        if(ireschk.eq.0) then ! T.Sato 2013/10/19
        call check_tcross(m,iax,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    ierr)
        endif

        if (ierr .ne. 0 ) goto 900

*-----------------------------------------------------------------------

        rewind(jsi(ioe))
        ill(jsn(ioe),ioe) = 1

*-----------------------------------------------------------------------

        idas0 = nmmax
        idas1 = lmmax
        idas2 = idas1 + itenm(m)
        idas3 = idas2 + itanm(m)
        idas4 = ( idas3 + ittnm(m) - 1 ) * 2 + 1
        idasa = idas3 + ittnm(m)

        idasz = itenm(m) * itpan(m) * 2
     &        * ( itrnm(m) + 1 ) * itznm(m)
     &        * itanm(m)
     &        * ittnm(m)
     &        * itmst(m)

! sumover
        idasz_sum = itenm(m) * itpan(m) * 2
     &        * ( itrnm(m) + 1 ) * itznm(m)
     &        * itanm(m)
     &        * ittnm(m)
     &        * itmst(m)


*-----------------------------------------------------------------------
*   read tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call read_crsreg(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itrcn(m),itenm(m),itanm(m),ittnm(m),
     &                    itmst(m),
     &                    itrss(m),
     &                    idas_itrcc(itrcc(m)),das_itrca(itrca(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))


        else if ( itmsh(m) .eq. 2 ) then

         if( itenclo(m) .eq. 1 ) then

          call read_crsrz_rcc(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itrnm(m),itznm(m),itenm(m),itanm(m),
     &                    ittnm(m),itmst(m),
     &                    das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))
!     &                    trRES(irestalm(m)+idasz))

         else

          call read_crsrz(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itrnm(m),itznm(m),itenm(m),itanm(m),
     &                    ittnm(m),itmst(m),
     &                    das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)),
     &                    trRES(irestalm(m)+idasz))

         endif

        else if ( itmsh(m) .eq. 3 ) then

         if( itenclo(m) .eq. 1 ) then

          call read_crsxyz_rpp(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itxnm(m),itynm(m),itznm(m),
     &                    itenm(m),itanm(m),ittnm(m),itmst(m),
     &                    das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                    das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

         else

          call read_crsxyz(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itxnm(m),itynm(m),itznm(m),
     &                    itenm(m),itanm(m),ittnm(m),itmst(m),
     &                    das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                    das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

         endif

        end if

  800   continue
*-----------------------------------------------------------------------
*   restore tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call restore_crsreg(m,
     &                    itpan(m),itrcn(m),itenm(m),itanm(m),ittnm(m),
     &                    itmst(m),
     &                    itrss(m),
     &                    idas_itrcc(itrcc(m)),das_itrca(itrca(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 2 ) then

         if( itenclo(m) .eq. 1 ) then

          call restore_crsrz_rcc(m,
     &                    itpan(m),itrnm(m),itznm(m),itenm(m),itanm(m),
     &                    ittnm(m),itmst(m),
     &                    das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)),
     &                    trRES(irestalm(m)+idasz))

         else

          call restore_crsrz(m,
     &                    itpan(m),itrnm(m),itznm(m),itenm(m),itanm(m),
     &                    ittnm(m),itmst(m),
     &                    das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)),
     &                    trRES(irestalm(m)+idasz))

         endif

        else if ( itmsh(m) .eq. 3 ) then

         if( itenclo(m) .eq. 1 ) then

          call restore_crsxyz_rpp(m,
     &                    itpan(m),itxnm(m),itynm(m),itznm(m),
     &                    itenm(m),itanm(m),ittnm(m),itmst(m),
     &                    das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                    das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

         else

          call restore_crsxyz(m,
     &                    itpan(m),itxnm(m),itynm(m),itznm(m),
     &                    itenm(m),itanm(m),ittnm(m),itmst(m),
     &                    das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                    das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

         endif

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

      subroutine read_crsreg(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                       np,nr,ne,na,nt,nm,mr,kr,ar,eb,ab,tb,tr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)

        common /tall21/ rtfac(itlmax)
        common /fact01/ facmax(itlmax)
        common /stat / istdev, irestart, ireschk

        dimension   kr(mr)
        dimension   eb(ne+1)
        dimension   ew(ne)
        dimension   ab(na+1)
        dimension   tb(nt+1)
        dimension   aw(na)
        dimension   tw(nt)
        dimension   ar(nr)

        dimension   tr(np,ne,na,nt,nr,nm,2)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        character chin*200, chlw*200, chcm*200

! sumover dummy
        dimension :: tott_sum(1,2)
        data nsame /1/

*-----------------------------------------------------------------------
*        energy axis  ( LET axis )
*-----------------------------------------------------------------------

        if( itaxs(m,iax) .eq. 1 .or.
     &      itaxs(m,iax) .eq. 14 ) then

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nr_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrm.inc'

          do 190 im = 1, nm
          do 190 ir = 1, nr, nrstepi
          do 190 ia = 1, na, nastepi
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
     &             dmm0, dmm1,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             ir+irloop-1,im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi)

            end do

! sumover
            ie = 1
            call psufreg_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nmstepi,
     &      ip,ie,ia,it,ir,im,nsame,tott_sum)

  190     continue

*-----------------------------------------------------------------------
*        reg axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 2 ) then

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nr_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrm.inc'

          do 290 im = 1, nm
          do 290 ie = 1, ne, nestepi
          do 290 ia = 1, na, nastepi
          do 290 it = 1, nt, ntstepi
          do 290 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#numarea',8,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nrstepi = 1
            do ir = 1, nr, nrstepi



              read(jsi,'(1x,i5,2x,1pe13.4,
     &                   1000(1pe13.4,0pf8.4))')
     &             idmm0, dmm1,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             ir+irloop-1,im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi)

            end do

! sumover
            ir = 1
            call psufreg_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nmstepi,
     &      ip,ie,ia,it,ir,im,nsame,tott_sum)


  290     continue

*-----------------------------------------------------------------------
*        time axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 9 ) then

          facmx = 1.d0  ! kitamura23/03/31

           nm_0 = 0
           nr_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrm.inc'

          do 390 im = 1, nm
          do 390 ir = 1, nr, nrstepi
          do 390 ie = 1, ne, nestepi
          do 390 ia = 1, na, nastepi
          do 390 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#t-lowert-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            ntstepi = 1
            do it = 1, nt, ntstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             ir+irloop-1,im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi)

            end do

! sumover
            it = 1
            call psufreg_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nmstepi,
     &      ip,ie,ia,it,ir,im,nsame,tott_sum)

  390     continue

*-----------------------------------------------------------------------
*        angle axis ! T.Sato 2024/12/17
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 8 .or. itaxs(m,iax) .eq. 10 ) then

          facmx = 1.d0  ! kitamura23/03/31

           nm_0 = 0
           nr_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrm.inc'

          do 490 im = 1, nm
          do 490 ir = 1, nr, nrstepi
          do 490 ie = 1, ne, nestepi
          do 490 it = 1, nt, ntstepi
          do 490 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#a-lowera-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nastepi = 1
            do ia = 1, na, nastepi

              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             ir+irloop-1,im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi)

            end do

! sumover
            ia = 1
            call psufreg_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nmstepi,
     &      ip,ie,ia,it,ir,im,nsame,tott_sum)

  490     continue

*-----------------------------------------------------------------------
         end if

  900 continue

      end subroutine



************************************************************************
*                                                                      *
      subroutine read_crsrz(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                      np,nr,nz,ne,na,nt,nm,rm,zm,eb,ab,tb,
     &                      tr,tz)

*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
        common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                  rtami(itlmax), rtama(itlmax), rtadl(itlmax)

        common /tall21/ rtfac(itlmax)
        common /fact02/ facmxr(itlmax), facmxz(itlmax)  ! kitamura23/03/31
        common /stat / istdev, irestart, ireschk

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

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        character chin*200, chlw*200, chcm*200

! sumover dummy
         dimension tott_sum(1,2)
         data nsame /1/

*-----------------------------------------------------------------------

         irf(ir,iz) = ir + ( iz - 1 ) * ( nr + 1 )
         izf(ir,iz) = iz + ( ir - 1 ) * ( nz + 1 )

*-----------------------------------------------------------------------
*        energy axis  ( LET axis )
*-----------------------------------------------------------------------

        if( itaxs(m,iax) .eq. 1 .or.
     &      itaxs(m,iax) .eq. 14 ) then

           nm_0 = 0
           nz_0 = 1
           nr_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrzm.inc'

*-----------------------------------------------------------------------
*          z-crossing
*-----------------------------------------------------------------------

          facmr = 1.d0
          facmz = 1.d0

          do 180 im = 1, nm
          do 180 iz = 1, nz + 1, nzstepi
          do 180 ir = 1, nr, nrstepi
          do 180 ia = 1, na, nastepi
          do 180 it = 1, nt, ntstepi
          do 180 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa2(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                             facmr,facmz,ierr)

            facmxr(m) = facmr
            facmxz(m) = facmz

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#e-lowere-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nestepi = 1
            do ie = 1, ne, nestepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      (((((((tz(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             izf(ir+irloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)

            end do

! sumover
            ie = 1
            call psufrz_sumover_getput_z(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &      ip,ie,ia,it,ir,iz,im,nsame,tott_sum)

  180     continue
*-----------------------------------------------------------------------
*          r-crossing
*-----------------------------------------------------------------------

          facmr = 1.d0
          facmz = 1.d0

           nm_0 = 0
           nz_0 = 0
           nr_0 = 1
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrzm.inc'

          nr_1 = 1
          if(rm(1) .le. 0.0) then
            nr_1 =2
            if(iloopmode .eq. 6) then
              nrstepi = nrstepi - 1
            end if
          end if
          do im = 1, nm
          do ir = nr_1, nr + 1, nrstepi
          do 190 iz = 1, nz, nzstepi
          do 190 ia = 1, na, nastepi
          do 190 it = 1, nt, ntstepi
          do 190 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa2(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                             facmr,facmz,ierr)

            facmxr(m) = facmr
            facmxz(m) = facmz

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#e-lowere-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nestepi = 1
            do ie = 1, ne, nestepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      (((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             irf(ir+irloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)

            end do

! sumover
            ie = 1
            call psufrz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &      ip,ie,ia,it,ir,iz,im,nsame,tott_sum)


  190     continue
          end do
          end do

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 5 ) then

           nm_0 = 0
           nz_0 = 1
           nr_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrzm.inc'

*-----------------------------------------------------------------------
*          z-crossing
*-----------------------------------------------------------------------

          facmr = 1.d0
          facmz = 1.d0

          do 280 im = 1, nm
          do 280 ie = 1, ne, nestepi
          do 280 ir = 1, nr, nrstepi
          do 280 ia = 1, na, nastepi
          do 280 it = 1, nt, ntstepi
          do 280 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa2(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                             facmr,facmz,ierr)

            facmxr(m) = facmr
            facmxz(m) = facmz

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#zsurfaceposition',17,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nzstepi = 1
            do iz = 1, nz + 1, nzstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      (((((((tz(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             izf(ir+irloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)

            end do

! sumover
            iz = 1
            call psufrz_sumover_getput_z(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &      ip,ie,ia,it,ir,iz,im,nsame,tott_sum)


  280     continue
*-----------------------------------------------------------------------
*          r-crossing
*-----------------------------------------------------------------------

          facmr = 1.d0
          facmz = 1.d0

           nm_0 = 0
           nz_0 = 0
           nr_0 = 1
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrzm.inc'

          nr_1 = 1
          if(rm(1) .le. 0.0) then
            nr_1 =2
            if(iloopmode .eq. 6) then
              nrstepi = nrstepi - 1
            end if
          end if
          do im = 1, nm
          do ir = nr_1, nr + 1, nrstepi
          do 290 ie = 1, ne, nestepi
          do 290 ia = 1, na, nastepi
          do 290 it = 1, nt, ntstepi
          do 290 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa2(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                             facmr,facmz,ierr)

            facmxr(m) = facmr
            facmxz(m) = facmz

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#z-lowerz-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nzstepi = 1
            do iz = 1, nz, nzstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      (((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             irf(ir+irloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)


            end do

! sumover
            iz = 1
            call psufrz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &      ip,ie,ia,it,ir,iz,im,nsame,tott_sum)

  290     continue
          end do
          end do

*-----------------------------------------------------------------------
*        r axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 6 ) then

*-----------------------------------------------------------------------
*        z-crossing
*-----------------------------------------------------------------------

          facmr = 1.d0
          facmz = 1.d0

           nm_0 = 0
           nz_0 = 1
           nr_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrzm.inc'

          do 380 im = 1, nm
          do 380 iz = 1, nz + 1, nzstepi
          do 380 ie = 1, ne, nestepi
          do 380 ia = 1, na, nastepi
          do 380 it = 1, nt, ntstepi
          do 380 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa2(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                             facmr,facmz,ierr)

            facmxr(m) = facmr
            facmxz(m) = facmz

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#r-lowerr-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nrstepi = 1
            do ir = 1, nr, nrstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      (((((((tz(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             izf(ir+irloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)

            end do

! sumover
            ir = 1
            call psufrz_sumover_getput_z(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &      ip,ie,ia,it,ir,iz,im,nsame,tott_sum)

  380     continue
*-----------------------------------------------------------------------
*          r-crossing
*-----------------------------------------------------------------------

          facmr = 1.d0
          facmz = 1.d0

           nm_0 = 0
           nz_0 = 0
           nr_0 = 1
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrzm.inc'

          do 390 im = 1, nm
          do 390 ie = 1, ne, nestepi
          do 390 iz = 1, nz, nzstepi
          do 390 ia = 1, na, nastepi
          do 390 it = 1, nt, ntstepi
          do 390 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa2(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                             facmr,facmz,ierr)

            facmxr(m) = facmr
            facmxz(m) = facmz

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#rsurfaceposition',17,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nrstepi = 1
            do ir = 1, nr + 1, nrstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      (((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             irf(ir+irloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)

            end do

! sumover
            ir = 1
            call psufrz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &      ip,ie,ia,it,ir,iz,im,nsame,tott_sum)

  390     continue

*-----------------------------------------------------------------------
*        angle axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 8 .or.
     &            itaxs(m,iax) .eq. 10 ) then

*-----------------------------------------------------------------------
         if(( itaxs(m,iax) .eq. 8  .and. itaty(m) .gt. 0 ) .or.
     &      ( itaxs(m,iax) .eq. 10 .and. itaty(m) .lt. 0 )) then

               iai = 1
               iaf = na
               iad = 1

         else if(( itaxs(m,iax) .eq. 8  .and. itaty(m) .lt. 0 ) .or.
     &           ( itaxs(m,iax) .eq. 10 .and. itaty(m) .gt. 0 )) then

               iai = na
               iaf = 1
               iad = -1

         end if

*-----------------------------------------------------------------------
*        z-crossing
*-----------------------------------------------------------------------

          facmr = 1.d0
          facmz = 1.d0

           nm_0 = 0
           nz_0 = 1
           nr_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrzm.inc'

          do 480 im = 1, nm
          do 480 iz = 1, nz + 1, nzstepi
          do 480 ir = 1, nr, nrstepi
          do 480 ie = 1, ne, nestepi
          do 480 it = 1, nt, ntstepi
          do 480 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa2(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                             facmr,facmz,ierr)

            facmxr(m) = facmr
            facmxz(m) = facmz

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#a-lowera-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            do ia = iai, iaf, iad

              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      (((((((tz(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             izf(ir+irloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)

            end do

! sumover
            ia= 1
            call psufrz_sumover_getput_z(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &      ip,ie,ia,it,ir,iz,im,nsame,tott_sum)

  480     continue
*-----------------------------------------------------------------------
*          r-crossing
*-----------------------------------------------------------------------

          facmr = 1.d0
          facmz = 1.d0

           nm_0 = 0
           nz_0 = 0
           nr_0 = 1
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrzm.inc'

          nr_1 = 1
          if(rm(1) .le. 0.0) then
            nr_1 =2
            if(iloopmode .eq. 6) then
              nrstepi = nrstepi - 1
            end if
          end if
          do im = 1, nm
          do ir = nr_1, nr + 1, nrstepi
          do 490 iz = 1, nz, nzstepi
          do 490 ie = 1, ne, nestepi
          do 490 it = 1, nt, ntstepi
          do 490 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa2(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                             facmr,facmz,ierr)

            facmxr(m) = facmr
            facmxz(m) = facmz

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#a-lowera-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            do ia = iai, iaf, iad

              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      (((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             irf(ir+irloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)

            end do

! sumover
            ia= 1
            call psufrz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &      ip,ie,ia,it,ir,iz,im,nsame,tott_sum)

  490     continue
          end do
          end do

*-----------------------------------------------------------------------
*        time axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 9 ) then

           nm_0 = 0
           nz_0 = 1
           nr_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrzm.inc'

*-----------------------------------------------------------------------
*          z-crossing
*-----------------------------------------------------------------------

          facmr = 1.d0
          facmz = 1.d0

          do 580 im = 1, nm
          do 580 iz = 1, nz + 1, nzstepi
          do 580 ir = 1, nr, nrstepi
          do 580 ie = 1, ne, nestepi
          do 580 ia = 1, na, nastepi
          do 580 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa2(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                             facmr,facmz,ierr)

            facmxr(m) = facmr
            facmxz(m) = facmz

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#t-lowert-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

             ntstepi = 1
            do it = 1, nt, ntstepi



              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      (((((((tz(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             izf(ir+irloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)

            end do

! sumover
            it= 1
            call psufrz_sumover_getput_z(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &      ip,ie,ia,it,ir,iz,im,nsame,tott_sum)

  580     continue
*-----------------------------------------------------------------------
*          r-crossing
*-----------------------------------------------------------------------

          facmr = 1.d0
          facmz = 1.d0

           nm_0 = 0
           nz_0 = 0
           nr_0 = 1
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrzm.inc'

          nr_1 = 1
          if(rm(1) .le. 0.0) then
            nr_1 =2
            if(iloopmode .eq. 6) then
              nrstepi = nrstepi - 1
            end if
          end if
          do im = 1, nm
          do ir = nr_1, nr + 1, nrstepi
          do 590 iz = 1, nz, nzstepi
          do 590 ie = 1, ne, nestepi
          do 590 ia = 1, na, nastepi
          do 590 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa2(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                             facmr,facmz,ierr)

            facmxr(m) = facmr
            facmxz(m) = facmz

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#t-lowert-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            ntstepi = 1
            do it = 1, nt, ntstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      (((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             irf(ir+irloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)


            end do

! sumover
            it= 1
            call psufrz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &      ip,ie,ia,it,ir,iz,im,nsame,tott_sum)

  590     continue
          end do
          end do

*-----------------------------------------------------------------------
         end if

  900 continue

      end subroutine



************************************************************************
*                                                                      *
      subroutine read_crsrz_rcc(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                      np,nr,nz,ne,na,nt,nm,rm,zm,eb,ab,tb,
     &                      tr)
! nais
!    &                      tz)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
        common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                  rtami(itlmax), rtama(itlmax), rtadl(itlmax)

        common /tall21/ rtfac(itlmax)
        common /fact01/ facmax(itlmax)
        common /stat / istdev, irestart, ireschk

        dimension   rm(nr+1)
        dimension   zm(nz+1)
        dimension   eb(ne+1)
        dimension   ab(na+1)
        dimension   tb(nt+1)
        dimension   ew(ne)
        dimension   aw(na)
        dimension   tw(nt)

        dimension   tr(np,ne,na,nt,(nr+1)*nz,nm,2)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        character chin*200, chlw*200, chcm*200
        character dc2*15

! sumover dummy
        dimension tott_sum(1,2)
        data nsame /1/

*-----------------------------------------------------------------------

         irf(ir,iz) = ir + ( iz - 1 ) * ( nr + 1 )

*-----------------------------------------------------------------------
*        rcc-crossing
*        energy axis  ( LET axis )
*-----------------------------------------------------------------------

        if( itaxs(m,iax) .eq. 1 .or.
     &      itaxs(m,iax) .eq. 14 ) then

*-----------------------------------------------------------------------

          facmx = 1.d0  ! kitamura23/03/31

           nm_0 = 0
           nz_0 = 0
           nr_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrzm.inc'

          do 190 im = 1, nm
          do 190 ir = 1, nr, nrstepi
          do 190 iz = 1, nz, nzstepi
          do 190 ia = 1, na, nastepi
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
     &             dmm0, dmm1,
     &      (((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             irf(ir+irloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)

            end do

! sumover
            ie= 1
            call psufrz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &      ip,ie,ia,it,ir,iz,im,nsame,tott_sum)

  190     continue

*-----------------------------------------------------------------------
*        rcc-crossing
*        z axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 5 ) then

*-----------------------------------------------------------------------

          facmx = 1.d0  ! kitamura23/03/31

           nm_0 = 0
           nz_0 = 0
           nr_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrzm.inc'

          do 290 im = 1, nm
          do 290 ir = 1, nr, nrstepi
          do 290 ie = 1, ne, nestepi
          do 290 ia = 1, na, nastepi
          do 290 it = 1, nt, ntstepi
          do 290 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#z-lowerz-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nzstepi = 1
            do iz = 1, nz, nzstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      (((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             irf(ir+irloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)

            end do

! sumover
            iz= 1
            call psufrz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &      ip,ie,ia,it,ir,iz,im,nsame,tott_sum)


  290     continue

*-----------------------------------------------------------------------
*        rcc-crossing
*        r axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 6 ) then

           nm_0 = 0
           nz_0 = 0
           nr_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrzm.inc'

*-----------------------------------------------------------------------
*          r-crossing
*-----------------------------------------------------------------------

          facmx = 1.d0  ! kitamura23/03/31

          do 390 im = 1, nm
          do 390 ie = 1, ne, nestepi
          do 390 iz = 1, nz, nzstepi
          do 390 ia = 1, na, nastepi
          do 390 it = 1, nt, ntstepi
          do 390 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#r-lowerr-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nrstepi = 1
            do ir = 1, nr, nrstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      (((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             irf(ir+irloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)

            end do

! sumover
            ir= 1
            call psufrz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &      ip,ie,ia,it,ir,iz,im,nsame,tott_sum)


  390     continue

*-----------------------------------------------------------------------
*        rcc-crossing
*        angle axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 .or.
     &            itaxs(m,iax) .eq. 10 ) then

*-----------------------------------------------------------------------

         if(( itaxs(m,iax) .eq. 8  .and. itaty(m) .gt. 0 ) .or.
     &      ( itaxs(m,iax) .eq. 10 .and. itaty(m) .lt. 0 )) then

               iai = 1
               iaf = na
               iad = 1

         else if(( itaxs(m,iax) .eq. 8  .and. itaty(m) .lt. 0 ) .or.
     &           ( itaxs(m,iax) .eq. 10 .and. itaty(m) .gt. 0 )) then

               iai = na
               iaf = 1
               iad = -1

         end if

*-----------------------------------------------------------------------

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nz_0 = 0
           nr_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrzm.inc'

          do 490 im = 1, nm
          do 490 ir = 1, nr, nrstepi
          do 490 iz = 1, nz, nzstepi
          do 490 ie = 1, ne, nestepi
          do 490 it = 1, nt, ntstepi
          do 490 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#a-lowera-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            do ia = iai, iaf, iad


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      (((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             irf(ir+irloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)

            end do

! sumover
            ia= 1
            call psufrz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &      ip,ie,ia,it,ir,iz,im,nsame,tott_sum)

  490     continue

*-----------------------------------------------------------------------
*        time axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9 ) then

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nz_0 = 0
           nr_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatrzm.inc'

*-----------------------------------------------------------------------

          do 590 im = 1, nm
          do 590 ir = 1, nr, nrstepi
          do 590 iz = 1, nz, nzstepi
          do 590 ie = 1, ne, nestepi
          do 590 ia = 1, na, nastepi
          do 590 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#t-lowert-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            ntstepi = 1
            do it = 1, nt, ntstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      (((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             irf(ir+irloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)

            end do

! sumover
            it= 1
            call psufrz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &      ip,ie,ia,it,ir,iz,im,nsame,tott_sum)

  590     continue

*-----------------------------------------------------------------------
*        rz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 13 ) then

          call sel_dc2(m,dc2,ldc2)
          if( ittwo(m) .eq. 4 ) then
            dc2 = '#rzfluxr.err'
            ldc2 = 12
          else if( ittwo(m) .eq. 5 ) then
            dc2 = 'r/z'
            ldc2 = 3
          end if

*-----------------------------------------------------------------------

          facmx = 1.d0  ! kitamura23/03/31

          do 690 im = 1, nm
          do 690 ip = 1, np
          do 690 ie = 1, ne
          do 690 ia = 1, na
          do 690 it = 1, nt

            call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                       dc2,ldc2,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               read(jsi,'(1p10e11.3)')
     &          ((tr(ip,ie,ia,it,irf(ir,iz),im,ioe), iz = 1,nz ),
     &                                            ir = nr,1,-1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

              do iz = 1, nz
              do ir = 1, nr

               read(jsi,'(1p3e11.3,0pf8.4)')
     &           dmm0, dmm1,
     &           tr(ip,ie,ia,it,irf(ir,iz),im,1),
     &           tr(ip,ie,ia,it,irf(ir,iz),im,2)

              end do
              end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

              do ir = nr, 1, -1

                read(jsi,'(1p1000e11.3)')
     &            dmm0,
     &            ( tr(ip,ie,ia,it,irf(ir,iz),im,ioe), iz = 1, nz )

              end do

*-----------------------------------------------------------------------
            end if

  690     continue

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
  900 continue

      end subroutine



************************************************************************
*                                                                      *
      subroutine read_crsxyz(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                       np,nx,ny,nz,ne,na,nt,nm,xm,ym,zm,eb,ab,tb,
     &                       tr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
        common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                  rtami(itlmax), rtama(itlmax), rtadl(itlmax)

        common /tall21/ rtfac(itlmax)
        common /fact01/ facmax(itlmax)
        common /stat / istdev, irestart, ireschk

        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)
        dimension   eb(ne+1)
        dimension   ab(na+1)
        dimension   tb(nt+1)
        dimension   ew(ne)
        dimension   aw(na)
        dimension   tw(nt)

        dimension   tr(np,ne,na,nt,nx*ny*(nz+1),nm,2)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        character chin*200, chlw*200, chcm*200
        character dc2*15
!
! sumover dummy
        dimension tott_sum(1,2)
        data nsame /1/

*-----------------------------------------------------------------------

        data igsh / 0 /

*-----------------------------------------------------------------------

        icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny

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
*        energy axis  ( LET axis )
*-----------------------------------------------------------------------
        if( itaxs(m,iax) .eq. 1 .or.
     &      itaxs(m,iax) .eq. 14 ) then

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nz_0 = 1
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatxyzm.inc'

          do 190 im = 1, nm
          do 190 iz = 1, nz + 1, nzstepi
          do 190 ix = 1, nx, nxstepi
          do 190 iy = 1, ny, nystepi
          do 190 ia = 1, na, nastepi
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
     &           dmm0, dmm1,

     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi)

            end do

! sumover
            ie= 1
            call psufxyz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nxstepi,nystepi,nzstepi,
     &      nmstepi,
     &      ip,ie,ia,it,ix,iy,iz,im,nsame,tott_sum)

  190     continue

*-----------------------------------------------------------------------
*        x axis
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 3 ) then

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nz_0 = 1
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatxyzm.inc'


          do 290 im = 1, nm
          do 290 iz = 1, nz + 1, nzstepi
          do 290 iy = 1, ny, nystepi
          do 290 ie = 1, ne, nestepi
          do 290 ia = 1, na, nastepi
          do 290 it = 1, nt, ntstepi
          do 290 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#x-lowerx-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nxstepi = 1
            do ix = 1, nx, nxstepi

              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &           dmm0, dmm1,

     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi)

            end do

! sumover
            ix= 1
            call psufxyz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nxstepi,nystepi,nzstepi,
     &      nmstepi,
     &      ip,ie,ia,it,ix,iy,iz,im,nsame,tott_sum)

  290     continue

*-----------------------------------------------------------------------
*        y axis
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 4 ) then

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nz_0 = 1
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatxyzm.inc'

          do 390 im = 1, nm, nmstepi
          do 390 iz = 1, nz + 1, nzstepi
          do 390 ix = 1, nx, nxstepi
          do 390 ie = 1, ne, nestepi
          do 390 ia = 1, na, nastepi
          do 390 it = 1, nt, ntstepi
          do 390 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#y-lowery-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nystepi = 1
            do iy = 1, ny, nystepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &           dmm0, dmm1,
     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi)

            end do

! sumover
            iy= 1
            call psufxyz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nxstepi,nystepi,nzstepi,
     &      nmstepi,
     &      ip,ie,ia,it,ix,iy,iz,im,nsame,tott_sum)

  390     continue

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 5 ) then

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nz_0 = 1
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatxyzm.inc'

          do 490 im = 1, nm
          do 490 ie = 1, ne, nestepi
          do 490 ia = 1, na, nastepi
          do 490 ix = 1, nx, nxstepi
          do 490 iy = 1, ny, nystepi
          do 490 it = 1, nt, ntstepi
          do 490 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#zsurfaceposition',17,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nzstepi = 1
            do iz = 1, nz + 1, nzstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &           dmm0, dmm1,
     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi)

            end do

! sumover
            iz= 1
            call psufxyz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nxstepi,nystepi,nzstepi,
     &      nmstepi,
     &      ip,ie,ia,it,ix,iy,iz,im,nsame,tott_sum)

  490     continue

*-----------------------------------------------------------------------
*        t axis
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 9 ) then

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nz_0 = 1
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatxyzm.inc'

          do 590 im = 1, nm
          do 590 iz = 1, nz + 1, nzstepi
          do 590 ix = 1, nx, nxstepi
          do 590 iy = 1, ny, nystepi
          do 590 ie = 1, ne, nestepi
          do 590 ia = 1, na, nastepi
          do 590 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#t-lowert-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            ntstepi = 1
            do it = 1, nt, ntstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &           dmm0, dmm1,
     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi)

            end do

! sumover
            call psufxyz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nxstepi,nystepi,nzstepi,
     &      nmstepi,
     &      ip,ie,ia,it,ix,iy,iz,im,nsame,tott_sum)

  590     continue

*-----------------------------------------------------------------------
*        xy axis
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 7 ) then

          call sel_dc2(m,dc2,ldc2)
          if( ittwo(m) .eq. 4 ) then
            dc2 = '#xyfluxr.err'
            ldc2 = 12
          else if( ittwo(m) .eq. 5 ) then
            dc2 = 'y/x'
            ldc2 = 3
          end if

          facmx = 1.d0  ! kitamura23/03/31

          do 690 im = 1, nm
          do 690 iz = 1, nz + 1
          do 690 ip = 1, npg
          do 690 ie = 1, neg
          do 690 ia = 1, nag
          do 690 it = 1, ntg

            call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

*-----------------------------------------------------------------------
            if ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               read(jsi,'(1p10e11.3)')
     &         ( ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), ix = 1, nx ),
     &             iy = ny, 1, -1 )

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 4 ) then

              do iy = 1, ny
              do ix = 1, nx

                read(jsi,'(1p3e11.3,0pf8.4)')
     &             dmm0, dmm1,
     &             tr(ip,ie,ia,it,icf(ix,iy,iz),im,1),
     &             tr(ip,ie,ia,it,icf(ix,iy,iz),im,2)

              end do
              end do

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 5 ) then

              do iy = ny, 1, -1

                read(jsi,'(1p1000e11.3)')
     &          dmm0,
     &          ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), ix = 1, nx )

              end do

*-----------------------------------------------------------------------
            end if

  690     continue

*-----------------------------------------------------------------------
        end if

  900   continue

      end subroutine



************************************************************************
*                                                                      *
      subroutine read_crsxyz_rpp(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                       np,nx,ny,nz,ne,na,nt,nm,xm,ym,zm,eb,ab,tb,
     &                       tr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
        common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                  rtami(itlmax), rtama(itlmax), rtadl(itlmax)

        common /tall21/ rtfac(itlmax)
        common /fact01/ facmax(itlmax)
        common /stat / istdev, irestart, ireschk

        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)
        dimension   eb(ne+1)
        dimension   ab(na+1)
        dimension   tb(nt+1)
        dimension   ew(ne)
        dimension   aw(na)
        dimension   tw(nt)

        dimension   tr(np,ne,na,nt,nx*ny*(nz+1),nm,2)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        character chin*200, chlw*200, chcm*200
        character dc2*15

! sumover dummy
        dimension tott_sum(1,2)
        data nsame /1/
*-----------------------------------------------------------------------

        data igsh / 0 /

*-----------------------------------------------------------------------

        icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny

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
*        energy axis  ( LET axis )
*-----------------------------------------------------------------------

        if( itaxs(m,iax) .eq. 1 .or.
     &      itaxs(m,iax) .eq. 14 ) then

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatxyzm.inc'

          do 190 im = 1, nm
          do 190 iz = 1, nz, nzstepi
          do 190 ix = 1, nx, nxstepi
          do 190 iy = 1, ny, nystepi
          do 190 ia = 1, na, nastepi
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
     &           dmm0, dmm1,
     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi)

            end do

! sumover
            ie= 1
            call psufxyz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nxstepi,nystepi,nzstepi,
     &      nmstepi,
     &      ip,ie,ia,it,ix,iy,iz,im,nsame,tott_sum)

  190     continue

*-----------------------------------------------------------------------
*        x axis
*-----------------------------------------------------------------------

        else if( itaxs(m,iax) .eq. 3 ) then

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatxyzm.inc'

          do 290 im = 1, nm
          do 290 iz = 1, nz, nzstepi
          do 290 iy = 1, ny, nystepi
          do 290 ie = 1, ne, nestepi
          do 290 ia = 1, na, nastepi
          do 290 it = 1, nt, ntstepi
          do 290 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#x-lowerx-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nxstepi = 1
            do ix = 1, nx, nxstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &           dmm0, dmm1,
     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi)

            end do

! sumover
            ix= 1
            call psufxyz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nxstepi,nystepi,nzstepi,
     &      nmstepi,
     &      ip,ie,ia,it,ix,iy,iz,im,nsame,tott_sum)

  290     continue

*-----------------------------------------------------------------------
*        y axis
*-----------------------------------------------------------------------

        else if( itaxs(m,iax) .eq. 4 ) then

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatxyzm.inc'

          do 390 im = 1, nm
          do 390 iz = 1, nz, nzstepi
          do 390 ix = 1, nx, nxstepi
          do 390 ie = 1, ne, nestepi
          do 390 ia = 1, na, nastepi
          do 390 it = 1, nt, ntstepi
          do 390 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#y-lowery-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nystepi = 1
            do iy = 1, ny, nystepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &           dmm0, dmm1,
     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi)

            end do

! sumover
            iy= 1
            call psufxyz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nxstepi,nystepi,nzstepi,
     &      nmstepi,
     &      ip,ie,ia,it,ix,iy,iz,im,nsame,tott_sum)

  390     continue

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------

        else if( itaxs(m,iax) .eq. 5 ) then

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatxyzm.inc'

          do 490 im = 1, nm
          do 490 ie = 1, ne, nestepi
          do 490 ia = 1, na, nastepi
          do 490 ix = 1, nx, nxstepi
          do 490 iy = 1, ny, nystepi
          do 490 it = 1, nt, ntstepi
          do 490 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#z-lowerz-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nzstepi = 1
            do iz = 1, nz, nzstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &           dmm0, dmm1,
     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi)

            end do

! sumover
            iz= 1
            call psufxyz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nxstepi,nystepi,nzstepi,
     &      nmstepi,
     &      ip,ie,ia,it,ix,iy,iz,im,nsame,tott_sum)

  490     continue

*-----------------------------------------------------------------------
*        t axis
*-----------------------------------------------------------------------

        else if( itaxs(m,iax) .eq. 9 ) then

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           na_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_peatxyzm.inc'

          do 590 im = 1, nm
          do 590 iz = 1, nz, nzstepi
          do 590 ix = 1, nx, nxstepi
          do 590 iy = 1, ny, nystepi
          do 590 ie = 1, ne, nestepi
          do 590 ia = 1, na, nastepi
          do 590 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#t-lowert-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            ntstepi = 1
            do it = 1, nt, ntstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &           dmm0, dmm1,
     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,
     &             ia+ialoop-1,it+itloop-1,
     &             icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &             im,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi)

            end do

! sumover
            it= 1
            call psufxyz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,nastepi,ntstepi,nxstepi,nystepi,nzstepi,
     &      nmstepi,
     &      ip,ie,ia,it,ix,iy,iz,im,nsame,tott_sum)

  590     continue

*-----------------------------------------------------------------------
*        xy axis
*-----------------------------------------------------------------------

        else if( itaxs(m,iax) .eq. 7 ) then

          call sel_dc2(m,dc2,ldc2)
          if( ittwo(m) .eq. 4 ) then
            dc2 = '#xyfluxr.err'
            ldc2 = 12
          else if( ittwo(m) .eq. 5 ) then
            dc2 = 'y/x'
            ldc2 = 3
          end if

          facmx = 1.d0  ! kitamura23/03/31

          do 690 im = 1, nm
          do 690 iz = 1, nz + 1
          do 690 ip = 1, npg
          do 690 ie = 1, neg
          do 690 ia = 1, nag
          do 690 it = 1, ntg

            call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

*-----------------------------------------------------------------------

            if ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               read(jsi,'(1p10e11.3)')
     &         ( ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), ix = 1, nx ),
     &             iy = ny, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

              do iy = 1, ny
              do ix = 1, nx


              end do
              end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

              do iy = ny, 1, -1

                read(jsi,'(1p1000e11.3)')
     &          dmm0,
     &          ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), ix = 1, nx )

              end do

*-----------------------------------------------------------------------

            end if

  690     continue

*-----------------------------------------------------------------------
*        yz axis
*-----------------------------------------------------------------------

        else if( itaxs(m,iax) .eq. 11 ) then

          call sel_dc2(m,dc2,ldc2)
          if( ittwo(m) .eq. 4 ) then
            dc2 = '#yzfluxr.err'
            ldc2 = 12
          else if( ittwo(m) .eq. 5 ) then
            dc2 = 'y/z'
            ldc2 = 3
          end if

          facmx = 1.d0  ! kitamura23/03/31

          do 790 im = 1, nm
          do 790 ix = 1, nx
          do 790 ip = 1, npg
          do 790 ie = 1, neg
          do 790 ia = 1, nag
          do 790 it = 1, ntg

            call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

*-----------------------------------------------------------------------

            if ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               read(jsi,'(1p10e11.3)')
     &         ( ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), iz = 1, nz ),
     &             iy = ny, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

              do iz = 1, nz
              do iy = 1, ny

                read(jsi,'(1p3e11.3,0pf8.4)')
     &             dmm0, dmm1,
     &             tr(ip,ie,ia,it,icf(ix,iy,iz),im,1),
     &             tr(ip,ie,ia,it,icf(ix,iy,iz),im,2)

              end do
              end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

              do iy = ny, 1, -1

                read(jsi,'(1p1000e11.3)')
     &          dmm0,
     &          ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), iz = 1, nz )

              end do

*-----------------------------------------------------------------------

            end if

  790     continue

*-----------------------------------------------------------------------
*        zx axis
*-----------------------------------------------------------------------

        else if( itaxs(m,iax) .eq. 12 ) then

          call sel_dc2(m,dc2,ldc2)
          if( ittwo(m) .eq. 4 ) then
            dc2 = '#xzfluxr.err'
            ldc2 = 12
          else if( ittwo(m) .eq. 5 ) then
            dc2 = 'x/z'
            ldc2 = 3
          end if

          facmx = 1.d0  ! kitamura23/03/31

          do 890 im = 1, nm
          do 890 iy = 1, ny
          do 890 ip = 1, npg
          do 890 ie = 1, neg
          do 890 ia = 1, nag
          do 890 it = 1, ntg

            call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

*-----------------------------------------------------------------------

            if ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               read(jsi,'(1p10e11.3)')
     &         ( ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), iz = 1, nz ),
     &             ix = nx, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

              do iz = 1, nz
              do ix = 1, nx

                read(jsi,'(1p3e11.3,0pf8.4)')
     &             dmm0, dmm1,
     &             tr(ip,ie,ia,it,icf(ix,iy,iz),im,1),
     &             tr(ip,ie,ia,it,icf(ix,iy,iz),im,2)

              end do
              end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

              do ix = nx, 1, -1

                read(jsi,'(1p1000e11.3)')
     &          dmm0,
     &          ( tr(ip,ie,ia,it,icf(ix,iy,iz),im,ioe), iz = 1, nz )

              end do

*-----------------------------------------------------------------------

            end if

  890     continue

*-----------------------------------------------------------------------

        end if

*-----------------------------------------------------------------------
  900   continue

      end subroutine



************************************************************************
*                                                                      *
      subroutine restore_crsreg(m,
     &                       np,nr,ne,na,nt,nm,mr,kr,ar,eb,ab,tb,tr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall21/ rtfac(itlmax)
        common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                  rtami(itlmax), rtama(itlmax), rtadl(itlmax)

        common /fact01/ facmax(itlmax)  ! kitamura23/03/31

        dimension   kr(mr)
        dimension   eb(ne+1)
        dimension   ew(ne)
        dimension   ab(na+1)
        dimension   tb(nt+1)
        dimension   aw(na)
        dimension   tw(nt)
        dimension   ar(nr)

        dimension   tr(np,ne,na,nt,nr,nm,2)

*-----------------------------------------------------------------------

        pi = 2.d0 * asin(1.d0)

*-----------------------------------------------------------------------
*           itunt(m) = 1 : /cm^2/source
*                    = 2 : /cm^2/MeV/source
*                    = 3 : /cm^2/Lethargy/source
*                    = 4 : /cm^2/source/SR
*                    = 5 : /cm^2/MeV/source/SR
*                    = 6 : /cm^2/Lethargy/source/SR
*-----------------------------------------------------------------------

         if( itunt(m) .eq.  1 .or. itunt(m) .eq.  4 .or.
     &       itunt(m) .eq. 11 .or. itunt(m) .eq. 14 ) then

            ew(1:ne) = 1.d+0
            ew_sum = 1.0d0

         else if( itunt(m) .eq.  2 .or. itunt(m) .eq.  5 .or.
     &            itunt(m) .eq. 12 .or. itunt(m) .eq. 15 ) then

            ew(1:ne) = eb(2:ne+1) - eb(1:ne)
            ew_sum = eb(ne+1) - eb(1)

         else if( itunt(m) .eq.  3 .or. itunt(m) .eq.  6 .or.
     &            itunt(m) .eq. 13 .or. itunt(m) .eq. 16 ) then

            ew(1:ne) = log( eb(2:ne+1) / eb(1:ne) )
            ew_sum = log( eb(ne+1) / eb(1) )

         end if

*-----------------------------------------------------------------------

         if( itunt(m) .eq.  1 .or. itunt(m) .eq.  2 .or.
     &       itunt(m) .eq.  3 .or. itunt(m) .eq. 11 .or.
     &       itunt(m) .eq. 12 .or. itunt(m) .eq. 13 ) then

            do i = 1, na

             aw(i) = 1.d+0

            end do
            aw_sum = 1.0d0

         else

            aw_sum = 0.0d0
            do i = 1, na

               if( itaty(m) .gt. 0 ) then

                  aw(i) = ( ab(i+1) - ab(i) ) * 2.0 * pi

               else

                  aw(i) = ( cos( ab(i) / 180.d0 * pi )
     &                    - cos( ab(i+1) / 180.d0 * pi ) ) * 2.0 * pi

               end if
               aw_sum = aw_sum + aw(i)

            end do

         end if


*-----------------------------------------------------------------------
*           itunt(m) = 11, 12, 13, 14 : /.../nsec/source
*-----------------------------------------------------------------------

        if( itunt(m) .gt. 10 ) then

          tw(1:nt) = tb(2:nt+1) - tb(1:nt)
          tw_sum = tb(nt+1) - tb(1)

        else

          tw(1:nt) = 1.d+0
          tw_sum = 1.0d0

        end if

        ar_sum = sum(ar(:))

*-----------------------------------------------------------------------

        if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

        do 100 im = 1, nm
        do 100 ir = 1, nr
        do 100 it = 1, nt
        do 100 ia = 1, na
        do 100 ie = 1, ne
        do 100 ip = 1, np

          call invert_stdev(m,A,B,
     &                      tr(ip,ie,ia,it,ir,im,1),
     &                      tr(ip,ie,ia,it,ir,im,2),
     &       ar(ir)*ew(ie)*aw(ia)*tw(it)/abs(rtfac(m)/facmax(m)))

          tr(ip,ie,ia,it,ir,im,1) = A
          tr(ip,ie,ia,it,ir,im,2) = B

! sumover
          fact_in = abs(rtfac(m)/facmax(m))
          call psufreg_sumover_stdev(1,m,ip,ie,ia,it,ir,im,
     &        fact_in,ew(ie),aw(ia),ar(ir),tw(it),
     &        ew_sum,aw_sum,ar_sum,tw_sum)


  100 continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_crsrz(m,
     &                         np,nr,nz,ne,na,nt,nm,rm,zm,eb,ab,tb,
     &                         tr,tz)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall21/ rtfac(itlmax)
        common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                  rtami(itlmax), rtama(itlmax), rtadl(itlmax)

        common /fact02/ facmxr(itlmax), facmxz(itlmax) ! kitamura23/03/31

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

         real(8), allocatable :: ar_r(:), ar_z(:)

*-----------------------------------------------------------------------

        parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

         irf(ir,iz) = ir + ( iz - 1 ) * ( nr + 1 )
         izf(ir,iz) = iz + ( ir - 1 ) * ( nz + 1 )

*-----------------------------------------------------------------------
*        set mesh area
*-----------------------------------------------------------------------

               az(ir) = pi * ( rm(ir+1)**2 - rm(ir)**2 )

               ar(ir,iz) = 2.0 * pi * rm(ir)
     &                     * ( zm(iz+1) - zm(iz) )

*-----------------------------------------------------------------------
*           itunt(m) = 1 : /cm^2/source
*                    = 2 : /cm^2/MeV/source
*                    = 3 : /cm^2/Lethargy/source
*                    = 4 : /cm^2/source/SR
*                    = 5 : /cm^2/MeV/source/SR
*                    = 6 : /cm^2/Lethargy/source/SR
*-----------------------------------------------------------------------
! sumover
            az_sum = 0.0d0
            do ir=1,nr
              az_sum = az_sum + az(ir)
            enddo
            allocate (ar_r(nz), ar_z(nr+1))
            ar_r(:) = 0.0d0
            ar_z(:) = 0.0d0
            do iz=1,nz
              do ir = 1,nr+1
                ar_r(iz) = ar_r(iz) + ar(ir,iz)
                ar_z(ir) = ar_z(ir) + ar(ir,iz)
              enddo
            enddo

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

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  2 .or.
     &          itunt(m) .eq.  3 .or. itunt(m) .eq. 11 .or.
     &          itunt(m) .eq. 12 .or. itunt(m) .eq. 13 ) then

               do i = 1, na

                  aw(i) = 1.d+0

               end do
               aw_sum = 1.0d0

            else

               aw_sum = 0.0d0
               do i = 1, na

                  if( itaty(m) .gt. 0 ) then

                     aw(i) = ( ab(i+1) - ab(i) ) * 2.0 * pi

                  else

                     aw(i) = ( cos( ab(i) / 180.d0 * pi )
     &                       - cos( ab(i+1) / 180.d0 * pi ) ) * 2.0 * pi

                  end if
                  aw_sum = a_sum + aw(i)

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
*        z-crossing
*-----------------------------------------------------------------------

            if( facmxz(m) .eq. 0.d0) facmxz(m) = 1.d0  ! kitamura23/03/31

            do 100 im = 1, nm
            do 100 iz = 1, nz + 1
            do 100 ir = 1, nr
            do 100 it = 1, nt
            do 100 ia = 1, na
            do 100 ie = 1, ne
            do 100 ip = 1, np

              call invert_stdev(m,A,B,
     &                          tz(ip,ie,ia,it,izf(ir,iz),im,1),
     &                          tz(ip,ie,ia,it,izf(ir,iz),im,2),
     &           az(ir)*ew(ie)*aw(ia)*tw(it)/abs(rtfac(m)/facmxz(m)))

              tz(ip,ie,ia,it,izf(ir,iz),im,1) = A
              tz(ip,ie,ia,it,izf(ir,iz),im,2) = B

! sumover
              fact_in = abs(rtfac(m)/facmxz(m))
              call psufrz_sumover_tz_stdev(1,m,ip,ie,ia,it,ir,iz,im,
     &                fact_in,ew(ie),aw(ia),tw(it),az(ir),
     &                ew_sum,aw_sum,tw_sum,az_sum)

  100       continue

*-----------------------------------------------------------------------
*        r-crossing
*-----------------------------------------------------------------------

            if( facmxr(m) .eq. 0.d0) facmxr(m) = 1.d0  ! kitamura23/03/31

            do 200 im = 1, nm
            do 200 iz = 1, nz
            do 200 ir = 1, nr + 1
            do 200 it = 1, nt
            do 200 ia = 1, na
            do 200 ie = 1, ne
            do 200 ip = 1, np

              call invert_stdev(m,A,B,
     &                          tr(ip,ie,ia,it,irf(ir,iz),im,1),
     &                          tr(ip,ie,ia,it,irf(ir,iz),im,2),
     &           ar(ir,iz)*ew(ie)*aw(ia)*tw(it)/abs(rtfac(m)/facmxr(m)))

              tr(ip,ie,ia,it,irf(ir,iz),im,1) = A
              tr(ip,ie,ia,it,irf(ir,iz),im,2) = B

! sumover
              fact_in = abs(rtfac(m)/facmxz(m))
              call psufrz_sumover_stdev(1,m,ip,ie,ia,it,ir,iz,im,
     &             fact_in,ew(ie),aw(ia),tw(it),ar(ir,iz),
     &             ew_sum,aw_sum,tw_sum,ar_r(iz),ar_z(ir))


  200       continue

      deallocate (ar_r, ar_z)

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_crsrz_rcc(m,
     &                         np,nr,nz,ne,na,nt,nm,rm,zm,eb,ab,tb,
     &                         tr,tz)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall21/ rtfac(itlmax)
        common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                  rtami(itlmax), rtama(itlmax), rtadl(itlmax)

        common /fact01/ facmax(itlmax)  ! kitamura23/03/31

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

        real(8), allocatable :: ar_r(:),ar_z(:)

*-----------------------------------------------------------------------

        parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

         irf(ir,iz) = ir + ( iz - 1 ) * ( nr + 1 )

*-----------------------------------------------------------------------
*        set mesh area ... two circules and two side walls
*-----------------------------------------------------------------------

         ar(ir,iz) = 2.d0 * pi * ( rm(ir+1)**2 - rm(ir)**2 )
     &             + 2.d0 * pi * ( rm(ir+1) + rm(ir) )
     &               * ( zm(iz+1) - zm(iz) )

*-----------------------------------------------------------------------
*           itunt(m) = 1 : /cm^2/source
*                    = 2 : /cm^2/MeV/source
*                    = 3 : /cm^2/Lethargy/source
*                    = 4 : /cm^2/source/SR
*                    = 5 : /cm^2/MeV/source/SR
*                    = 6 : /cm^2/Lethargy/source/SR
*-----------------------------------------------------------------------
! sumover
         allocate (ar_r(nz),ar_z(nr))
         ar_r(:) = 0.0d0
         ar_z(:) = 0.0d0
         do iz= 1, nz
           do ir = 1, nr
             ar_r(iz) = ar_r(iz) + ar(ir,iz)
             ar_z(ir) = ar_z(ir) + ar(ir,iz)
           enddo
         enddo

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

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  2 .or.
     &          itunt(m) .eq.  3 .or. itunt(m) .eq. 11 .or.
     &          itunt(m) .eq. 12 .or. itunt(m) .eq. 13 ) then

               do i = 1, na

                  aw(i) = 1.d+0

               end do
               aw_sum = 1.0d0

            else

               aw_sum = 0.0d0
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
               tw_sum = 1.d+0

            end if

*-----------------------------------------------------------------------
*        rcc-crossing
*-----------------------------------------------------------------------

            if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

            do 200 im = 1, nm
            do 200 iz = 1, nz
            do 200 ir = 1, nr
            do 200 it = 1, nt
            do 200 ia = 1, na
            do 200 ie = 1, ne
            do 200 ip = 1, np

              call invert_stdev(m,A,B,
     &                          tr(ip,ie,ia,it,irf(ir,iz),im,1),
     &                          tr(ip,ie,ia,it,irf(ir,iz),im,2),
     &           ar(ir,iz)*ew(ie)*aw(ia)*tw(it)/abs(rtfac(m)/facmax(m)))

              tr(ip,ie,ia,it,irf(ir,iz),im,1) = A
              tr(ip,ie,ia,it,irf(ir,iz),im,2) = B

! sumover
              fact_in = abs(rtfac(m)/facmax(m))
              call psufrz_sumover_stdev(1,m,ip,ie,ia,it,ir,iz,im,
     &              fact_in,ew(ie),aw(ia),tw(it),ar(ir,iz),
     &              ew_sum,aw_sum,tw_sum,ar_r(iz),ar_z(ir))


  200       continue

      deallocate (ar_r,ar_z)


      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_crsxyz(m,
     &                       np,nx,ny,nz,ne,na,nt,nm,xm,ym,zm,eb,ab,tb,
     &                       tr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall21/ rtfac(itlmax)
        common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                  rtami(itlmax), rtama(itlmax), rtadl(itlmax)

        common /fact01/ facmax(itlmax)  ! kitamura23/03/31

        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)
        dimension   eb(ne+1)
        dimension   ab(na+1)
        dimension   tb(nt+1)
        dimension   ew(ne)
        dimension   aw(na)
        dimension   tw(nt)
        dimension   tr(np,ne,na,nt,nx*ny*(nz+1),nm,2)

        real(8),allocatable :: ax_x(:),ax_y(:)

*-----------------------------------------------------------------------

        icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny

        ax(ix,iy) = ( xm(ix+1) - xm(ix) )
     &            * ( ym(iy+1) - ym(iy) )

        pi = acos(-1.d0)

*-----------------------------------------------------------------------
*           itunt(m) = 1 : /cm^2/source
*                    = 2 : /cm^2/MeV/source
*                    = 3 : /cm^2/Lethargy/source
*                    = 4 : /cm^2/source/SR
*                    = 5 : /cm^2/MeV/source/SR
*                    = 6 : /cm^2/Lethargy/source/SR
*-----------------------------------------------------------------------
! sumover
            allocate (ax_x(ny),ax_y(nx))
            ax_x(:) = 0.0d0
            ax_y(:) = 0.0d0
            do iy=1,ny
              do ix=1,nx
                ax_x(iy) = ax_x(iy)  + ax(ix,iy)
                ax_y(ix) = ax_y(ix)  + ax(ix,iy)
              enddo
            enddo


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

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  2 .or.
     &          itunt(m) .eq.  3 .or. itunt(m) .eq. 11 .or.
     &          itunt(m) .eq. 12 .or. itunt(m) .eq. 13 ) then

               do i = 1, na

                  aw(i) = 1.d+0

               end do
               aw_sum = 1.0d0

            else

               aw_sum = 0.0d0
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
*        z-crossing
*-----------------------------------------------------------------------

        if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

          do 100 im = 1, nm
          do 100 iz = 1, nz + 1
          do 100 ix = 1, nx
          do 100 iy = 1, ny
          do 100 ie = 1, ne
          do 100 ia = 1, na
          do 100 ip = 1, np
          do 100 it = 1, nt

            call invert_stdev(m,A,B,
     &                        tr(ip,ie,ia,it,icf(ix,iy,iz),im,1),
     &                        tr(ip,ie,ia,it,icf(ix,iy,iz),im,2),
     &         ax(ix,iy)*ew(ie)*aw(ia)*tw(it)/abs(rtfac(m)/facmax(m)))


            tr(ip,ie,ia,it,icf(ix,iy,iz),im,1) = A
            tr(ip,ie,ia,it,icf(ix,iy,iz),im,2) = B

! sumover
            fact_in = abs(rtfac(m)/facmax(m))
            call psufxyz_sumover_stdev(1,m,
     &            ip,ie,ia,it,ix,iy,iz,im,
     &            fact_in,ew(ie),aw(ia),tw(it),ax(ix,iy),
     &            ew_sum,aw_sum,tw_sum,ax_x(iy),ax_y(ix))

  100     continue

      deallocate (ax_x,ax_y)

      end subroutine

************************************************************************
*                                                                      *
      subroutine restore_crsxyz_rpp(m,
     &                       np,nx,ny,nz,ne,na,nt,nm,xm,ym,zm,eb,ab,tb,
     &                       tr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall21/ rtfac(itlmax)
        common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                  rtami(itlmax), rtama(itlmax), rtadl(itlmax)

        common /fact01/ facmax(itlmax)  ! kitamura23/03/31

        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)
        dimension   eb(ne+1)
        dimension   ab(na+1)
        dimension   tb(nt+1)
        dimension   ew(ne)
        dimension   aw(na)
        dimension   tw(nt)
        dimension   tr(np,ne,na,nt,nx*ny*(nz+1),nm,2)

        real(8), allocatable :: ax_x(:,:), ax_y(:,:), ax_z(:,:)
*-----------------------------------------------------------------------

        icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny

*-----------------------------------------------------------------------
*        set mesh area for rpp
*-----------------------------------------------------------------------

      ax(ix,iy,iz) = 2.d0 * ( xm(ix+1)-xm(ix) ) * ( ym(iy+1)-ym(iy) )
     &             + 2.d0 * ( ym(iy+1)-ym(iy) ) * ( zm(iz+1)-zm(iz) )
     &             + 2.d0 * ( zm(iz+1)-zm(iz) ) * ( xm(ix+1)-xm(ix) )

        pi = acos(-1.d0)

*-----------------------------------------------------------------------
*           itunt(m) = 1 : /cm^2/source
*                    = 2 : /cm^2/MeV/source
*                    = 3 : /cm^2/Lethargy/source
*                    = 4 : /cm^2/source/SR
*                    = 5 : /cm^2/MeV/source/SR
*                    = 6 : /cm^2/Lethargy/source/SR
*-----------------------------------------------------------------------
       allocate (ax_x(ny,nz),ax_y(nx,nz),ax_z(nx,ny))
       ax_x(:,:) = 0.0d0
       ax_y(:,:) = 0.0d0
       ax_y(:,:) = 0.0d0
       do iz = 1,nz
         do iy = 1,ny
           do ix = 1,nx
             ax_x(iy,iz) = ax_x(iy,iz) + ax(ix,iy,iz)
             ax_y(ix,iz) = ax_y(ix,iz) + ax(ix,iy,iz)
             ax_z(ix,iy) = ax_z(ix,iy) + ax(ix,iy,iz)
           enddo
          enddo
       enddo

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

            if( itunt(m) .eq.  1 .or. itunt(m) .eq.  2 .or.
     &          itunt(m) .eq.  3 .or. itunt(m) .eq. 11 .or.
     &          itunt(m) .eq. 12 .or. itunt(m) .eq. 13 ) then

               do i = 1, na

                  aw(i) = 1.d+0

               end do
               aw_sum = 1.0d0

            else

               aw_sum = 0.0d0
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
*        rpp-crossing
*-----------------------------------------------------------------------

          if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

          do 100 im = 1, nm
          do 100 iz = 1, nz
          do 100 ix = 1, nx
          do 100 iy = 1, ny
          do 100 ie = 1, ne
          do 100 ia = 1, na
          do 100 ip = 1, np
          do 100 it = 1, nt

            call invert_stdev(m,A,B,
     &                        tr(ip,ie,ia,it,icf(ix,iy,iz),im,1),
     &                        tr(ip,ie,ia,it,icf(ix,iy,iz),im,2),
     &       ax(ix,iy,iz)*ew(ie)*aw(ia)*tw(it)/abs(rtfac(m)/facmax(m)))

            tr(ip,ie,ia,it,icf(ix,iy,iz),im,1) = A
            tr(ip,ie,ia,it,icf(ix,iy,iz),im,2) = B

! sumover
            fact_in = abs(rtfac(m)/facmax(m))
            call psufxyz_sumover_rpp_stdev(1,m,
     &            ip,ie,ia,it,ix,iy,iz,im,
     &            fact_in,ew(ie),aw(ia),tw(it),ax(ix,iy,iz),
     &            ew_sum,aw_sum,tw_sum,
     &            ax_x(iy,iz),ax_y(ix,iz),ax_z(ix,iy))

  100     continue

       deallocate (ax_x,ax_y,ax_z)


      end subroutine
