
************************************************************************
*                                                                      *
      subroutine check_tdeposit2(m,iax,jsn,jsi,dsin,idsi,ill,ilf,ierr)
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
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

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

      common /tall31/ iterl(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall37/ itcnt(9,itlmax)

      common /tall50/ itlmt(itlmax), ite2l(itlmax)   ! CCSE 2022.08.31
      common /tall53/ itdfn(itlmax,2)

      common /tall54/ itdfn2(itlmax,2)
      common /tall55/ itlmt2(itlmax)

      common /tall56/ itety2(itlmax), itenm2(itlmax), iterg2(itlmax),
     &                rtemi2(itlmax), rtema2(itlmax), rtedl2(itlmax)
!OBINATA: for resfile
      common /tall59/ irfll(itlmax), crfln(itlmax), itrff(itlmax)
      character crfln*100
      character irfile*100

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

      dimension lschn(28), ischn(28)
      character schan(28)*8

      data icsu / 28 /

      data ( schan(i), i = 1, 28 ) /
     &    'mesh    ','part    ','e1-type ','e2-type ','unit    ',
     &    'axis    ','file    ','title   ','angel   ','2d-type ',
     &    'factor  ','x-txt   ','y-txt   ','z-txt   ','iechrl  ',
     &    'epsout  ','ctmin(1)','ctmax(1)','ctmin(2)','ctmax(2)',
     &    'ctmin(3)','ctmax(3)','letmat1 ','letmat2 ','dedxfnc1',
     &    'dedxfnc2','t-type  ','resfile '/

      data ( lschn(i), i = 1, 28 ) /
     &     4,         4,         7,         7,         4,
     &     4,         4,         5,         5,         7,
     &     6,         5,         5,         5,         6,
     &     6,         8,         8,         8,         8,
     &     8,         8,         7,         7,         8,
     &     8,         6,         7/

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

      character tname*12
      data      tname /'[t-deposit2]'/

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
            inpat = 0
            inaxi = 0
            infil = 0
            idtyp = 3
            ietp1  = 0
            ietp2  = 0

            langel = 0
            lxtxt  = 0
            lytxt  = 0
            lztxt  = 0
            iechrl = 72
            ieps   = 0
            ittp   = 0

            letmat1= 0
            letmat2= 0
            idxfn1 = 0
            idxfn2 = 0

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

            if( ipm .ne.  6 .and. ipm .ne. 7 .and.
     &          ipm .ne.  2 .and.
     &          ischn( ipm ) .gt. 1 ) goto 989

               ic = inumc(chlw,il+1,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 997

               icl = inumc(chlw,ic,i3,';') - 1

*-----------------------------------------------------------------------
*        mesh = region
*-----------------------------------------------------------------------

         if( ipm .eq. 1 ) then

            if( chlw(ic:ic+2) .eq. 'reg' ) then

               imesh = 1

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
*        eng1,2 mesh
*-----------------------------------------------------------------------

         else if( ipm .eq. 3 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

               ietp1 = nint( cvvv )

               call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     'e',ietp1,ine1,emin1,emax1,edel1,isteg1)


               goto 150

         else if( ipm .eq. 4 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

               ietp2 = nint( cvvv )

               call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     'e',ietp2,ine2,emin2,emax2,edel2,isteg2)


               goto 150

*-----------------------------------------------------------------------
*        time mesh
*-----------------------------------------------------------------------

         else if( ipm .eq. 27 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

               ittp = nint( cvvv )

               call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     't',ittp,int,tmin,tmax,tdel,isttg)


               goto 150

*-----------------------------------------------------------------------
*        unit
*-----------------------------------------------------------------------

         else if( ipm .eq. 5 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               iunt = nint( cvvv )

               if( iunt .lt. 1 .or.
     &             iunt .gt. 2 ) goto 988

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        axis
*-----------------------------------------------------------------------

         else if( ipm .eq. 6 ) then

  600       inaxi = inaxi + 1

            if( inaxi .gt. 6 ) goto 991

            if( chlw(ic:ic+3) .eq. 'eng1' ) then

               iaxis(inaxi) = 1
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+3) .eq. 'eng2' ) then

               iaxis(inaxi) = 2
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+2) .eq. 'e12' ) then

               iaxis(inaxi) = 3
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic+2) .eq. 'e21' ) then

               iaxis(inaxi) = 4
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic+3) .eq. 't-e1' ) then

               iaxis(inaxi) = 5
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+3) .eq. 'e1-t' ) then

               iaxis(inaxi) = 6
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+3) .eq. 't-e2' ) then

               iaxis(inaxi) = 7
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+3) .eq. 'e2-t' ) then

               iaxis(inaxi) = 8
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic) .eq. 't' ) then

               iaxis(inaxi) = 9
               ic =jnumc(chlw,ic+1,icl)

            else

               goto 992

            end if

               if( ic .le. icl ) goto 600

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        file name
*-----------------------------------------------------------------------

         else if( ipm .eq. 7 ) then

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
         else if( ipm .eq. 28 ) then

               icf = min( inumc(chlw,ic,icl,' ') - 1, icl )

               lrfile = icf-ic+1

               irfile(1:lrfile) = chin(ic:icf)

               irfflg = 1

*-----------------------------------------------------------------------
*        title
*-----------------------------------------------------------------------

         else if( ipm .eq. 8 ) then

               ict = min( ic + 79, i2 )

               title = chin(ic:ict)

               titll = ict - ic + 1

*-----------------------------------------------------------------------
*        angel parameters
*-----------------------------------------------------------------------

         else if( ipm .eq. 9 ) then

               ict = min( ic + 199, i2 )

               angelp = chin(ic:ict)

               langel = ict - ic + 1

*-----------------------------------------------------------------------
*        x-txt
*-----------------------------------------------------------------------

         else if( ipm .eq. 12 ) then

               ict = min( ic + 199, i2 )

               cxtxt = chin(ic:ict)

               lxtxt = ict - ic + 1

*-----------------------------------------------------------------------
*        y-txt
*-----------------------------------------------------------------------

         else if( ipm .eq. 13 ) then

               ict = min( ic + 199, i2 )

               cytxt = chin(ic:ict)

               lytxt = ict - ic + 1

*-----------------------------------------------------------------------
*        z-txt
*-----------------------------------------------------------------------

         else if( ipm .eq. 14 ) then

               ict = min( ic + 199, i2 )

               cztxt = chin(ic:ict)

               lztxt = ict - ic + 1

*-----------------------------------------------------------------------
*        dedxfnc1,2
*-----------------------------------------------------------------------

         else if( ipm .eq. 25 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               idxfn1 = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

         else if( ipm .eq. 26 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               idxfn2 = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        reg echo length
*-----------------------------------------------------------------------

         else if( ipm .eq. 15 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               iechrl = nint( cvvv )

               if( iechrl .lt. 40 ) goto 997

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
*        epsout
*-----------------------------------------------------------------------

         else if( ipm .eq. 16 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ieps = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        counter
*-----------------------------------------------------------------------

         else if( ipm .ge. 17 .and. ipm .le. 22 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               icount((ipm-17+2)/2) = 1
               icount( ipm-17+4)    = nint( cvvv )

               if( abs( nint(cvvv) ) .gt. 9999 ) goto 978

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        letmat1,2
*-----------------------------------------------------------------------

         else if( ipm .eq. 23 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               letmat1 = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

         else if( ipm .eq. 24 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               letmat2 = nint( cvvv )

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

      call check_reg(m,rglnrf,iec,cepn,lepn,ierr)

      call check_unit(m,iunt,iec,cepn,lepn,ierr)

      call check_axis(m,iaxis,inaxi,iec,cepn,lepn,ierr)

      do i = 1, itaxn(m)

        if ( any( itaxs(m,i) .eq. (/ 3, 4, 5, 6, 7, 8 /))) then
          call check_2dtype(m,idtyp,iec,cepn,lepn,ierr)
        end if

      end do

      call check_e1type(m,iec,cepn,lepn,ietp1,emax1,emin1,ine1)
      call check_e2type(m,iec,cepn,lepn,ietp2,emax2,emin2,ine2)

      call check_type('t',iec,cepn,lepn,
     &                ittty(m),rttma(m),rttmi(m),ittnm(m),
     &                ittp, tmax, tmin, int)

      call check_factor(m,rfact,iec,cepn,lepn,ierr)

      call check_part(m,iptyp,ipnkf,inpat,iec,cepn,lepn,ierr)

      call check_letmat1(m,letmat1,iec,cepn,lepn,ierr)
      call check_letmat2(m,letmat2,iec,cepn,lepn,ierr)

      call check_dedxfnc1(m,idxfn1,iec,cepn,lepn,ierr)
      call check_dedxfnc2(m,idxfn2,iec,cepn,lepn,ierr)

*-----------------------------------------------------------------------

      call judge_tall_check(iec,cepn,lepn,dsin(jsn),idsi(jsn),ierr)

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  950    write(dkam,'(i9)') mdas
         m_err = 'Total tally storage number exceeds mdas ='//dkam
         ErrCha = ''
         ErrID = 'L:800/R:check_tdeposit2/F:restdeposit2.f'
         goto 999

  975    m_err = 'number of reg should be two'
         ErrCha = ''
         ErrID = 'L:805/R:check_tdeposit2/F:restdeposit2.f'
         goto 999

  976    m_err = 'dedxfnc1(2) should be 0, 1, 2'
         ErrCha = ''
         ErrID = 'L:810/R:check_tdeposit2/F:restdeposit2.f'
         goto 999

  978    m_err = 'counter should be from -9999 to 9999'
         ErrCha = ''
         ErrID = 'L:815/R:check_tdeposit2/F:restdeposit2.f'
         goto 999

  983    m_err = '2D-type should be 1-7 in tally '//tname
         ErrCha = ''
         ErrID = 'L:820/R:check_tdeposit2/F:restdeposit2.f'
         goto 999

  986    write(dkam,'(i5)') itlmax
         m_err = 'Total tally number exceeds itlmax = '//dkam
         ErrCha = ''
         ErrID = 'L:826/R:check_tdeposit2/F:restdeposit2.f'
         goto 999

  987    m_err = 'Unknown parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:831/R:check_tdeposit2/F:restdeposit2.f'
         goto 999

  988    m_err = 'Unit should be 1 - 2 in tally '//tname
         ErrCha = ''
         ErrID = 'L:836/R:check_tdeposit2/F:restdeposit2.f'
         goto 999

  989    m_err = 'Double definition of the parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:841/R:check_tdeposit2/F:restdeposit2.f'
         goto 999

  990    m_err = 'Too many file name in tally '//tname
         ErrCha = ''
         ErrID = 'L:846/R:check_tdeposit2/F:restdeposit2.f'
         goto 999

  991    m_err = 'Too many axis in tally '//tname
         ErrCha = ''
         ErrID = 'L:851/R:check_tdeposit2/F:restdeposit2.f'
         goto 999

  992    m_err = 'Unknown axis name in tally '//tname
         ErrCha = ''
         ErrID = 'L:856/R:check_tdeposit2/F:restdeposit2.f'
         goto 999

  993    m_err = 'In this line, [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:861/R:check_tdeposit2/F:restdeposit2.f'
         goto 999

  994    m_err = 'Name of particle is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:866/R:check_tdeposit2/F:restdeposit2.f'
         goto 999

  995    m_err = 'Number of particles is larger than 6 in tally '//tname
         ErrCha = ''
         ErrID = 'L:871/R:check_tdeposit2/F:restdeposit2.f'
         goto 999

  996    m_err = 'Description of mesh is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:876/R:check_tdeposit2/F:restdeposit2.f'
         goto 999

  997    m_err = 'Description of parameter is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:881/R:check_tdeposit2/F:restdeposit2.f'
         goto 999

  998    m_err = 'Unknown mesh parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:886/R:check_tdeposit2/F:restdeposit2.f'
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
      subroutine read_tdeposit2(m,iax,ierr)
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

        common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                  rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
        common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)
        common /tall56/ itety2(itlmax), itenm2(itlmax), iterg2(itlmax),
     &                  rtemi2(itlmax), rtema2(itlmax), rtedl2(itlmax)

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

        if ( any( itaxs(m,iax) .eq. (/ 3, 4, 5, 6, 7, 8 /) )
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
        call check_tdeposit2(m,iax,jsn(ioe),jsi(ioe),
     &                   dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                   ierr)
        endif

        if (ierr .ne. 0 ) goto 900

*-----------------------------------------------------------------------

        rewind(jsi(ioe))
        ill(jsn(ioe),ioe) = 1

*-----------------------------------------------------------------------

          idas0 = nmmax
          idas1 = lmmax
          idas2 = idas1 + ittnm(m)
          idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
          idasa = idas2 + itrgn(m)
     &          + ( itrgn(m) + mod(itrgn(m),2) ) / 2

*-----------------------------------------------------------------------
*   read tally
*-----------------------------------------------------------------------

          call read_deposit2reg(m,iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itpan(m),itrgn(m),itrgm(m),itenm(m),itenm2(m),
     &      ittnm(m),idas_itreg(itreg(m)),
     &      das_iterg(iterg(m)),das_iterg2(iterg2(m)),
     &      das_ittrg(ittrg(m)),
     &      trRES(irestalm(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)))

  800   continue
*-----------------------------------------------------------------------
*   restore tally
*-----------------------------------------------------------------------

          call restore_deposit2reg(m,
     &      itpan(m),itrgn(m),itrgm(m),itenm(m),itenm2(m),
     &      ittnm(m),idas_itreg(itreg(m)),
     &      das_iterg(iterg(m)),das_iterg2(iterg2(m)),
     &      das_ittrg(ittrg(m)),
     &      trRES(irestalm(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)))

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
      subroutine read_deposit2reg(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                            np,nr,mr,ne1,ne2,nt,kr,eb1,eb2,tb,tr,
     &                            nvl,ivl,rvl)
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

        common /tall21/ rtfac(itlmax)
        common /fact01/ facmax(itlmax)
        common /stat / istdev, irestart, ireschk

        dimension   kr(mr)
        dimension   eb1(ne1+1)
        dimension   eb2(ne2+1)
        dimension   tb(nt+1)
        dimension   tw(nt)
        dimension   tr(np,0:ne1,0:ne2,nt,2)
        dimension   vl(nr)
        dimension   lr(nr)
        dimension   ivl(nvl)
        dimension   rvl(nvl)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        character dc2*20

*-----------------------------------------------------------------------
*        eng1 axis for deposit energy distribution
*-----------------------------------------------------------------------
         if( itaxs(m,iax) .eq. 1 ) then

           dc2  = '#e-lowere-upper'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31
           nt_0 = 0
           ne2_0 = 0
           ne1_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_pe1e2t.inc'

           do 190 ie2 = 1, ne2, ne2stepi
           do 190 it  = 1, nt, ntstepi
           do 190 ip  = 1, np, npstepi

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             ne1stepi = 1
             do ie1 = 1, ne1, ne1stepi


               read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &              dmm0, dmm1,
     &      (((((tr(ip+iploop-1,ie1+ie1loop-1,ie2+ie2loop-1,
     &              it+itloop-1,k),k=1,2),
     &      iploop=1,npstepi),ie1loop=1,ne1stepi),ie2loop=1,ne2stepi),
     &      itloop=1,ntstepi)

             end do

! sumover
            ie1= 1
            call  pdpst2reg_sumover_getput(jsi,m,iax,
     &            npstepi,ne1stepi,ne2stepi,ntstepi,
     &            ip,ie1,ie2,it,nsame,tott_sum)


  190     continue

*-----------------------------------------------------------------------
*        eng2 axis for deposit energy distribution
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 2 ) then

           dc2  = '#e-lowere-upper'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31
           nt_0 = 0
           ne2_0 = 0
           ne1_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_pe1e2t.inc'

           do 290 ie1 = 1, ne1, ne1stepi
           do 290 it  = 1, nt, ntstepi
           do 290 ip  = 1, np, npstepi

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             ne2stepi = 1
             do ie2 = 1, ne2, ne2stepi


               read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &              dmm0, dmm1,
     &      (((((tr(ip+iploop-1,ie1+ie1loop-1,ie2+ie2loop-1,
     &              it+itloop-1,k),k=1,2),
     &      iploop=1,npstepi),ie1loop=1,ne1stepi),ie2loop=1,ne2stepi),
     &      itloop=1,ntstepi)

             end do

! sumover
            ie2= 1
            call  pdpst2reg_sumover_getput(jsi,m,iax,
     &            npstepi,ne1stepi,ne2stepi,ntstepi,
     &            ip,ie1,ie2,it,nsame,tott_sum)


  290     continue

*-----------------------------------------------------------------------
*        e12 axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 3 ) then

           call sel_dc2(m,dc2,ldc2)
           if( ittwo(m) .eq. 4 ) then
             dc2 = '#e1e2numberr.err'
             ldc2 = 16
           else if( ittwo(m) .eq. 5 ) then
             dc2 = 'y/x'
             ldc2 = 3
           end if
*-----------------------------------------------------------------------

           facmx = 1.d0  ! kitamura23/03/31

           do 390 ip = 1, np
           do 390 it = 1, nt

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
     &         ( ( tr(ip,ie1,ie2,it,ioe),
     &             ie1 = 1, ne1 ), ie2 = ne2, 1, -1 )

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 4 ) then

               do ie2 = 1, ne2
               do ie1 = 1, ne1

                 read(jsi,'(1p3e13.4,0pf8.4)')
     &                dmm0, dmm1,
     &                (tr(ip,ie1,ie2,it,k),k=1,2)

               end do
               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 5 ) then

               do ie2 = ne2, 1, -1

                  read(jsi,'(1p1000e11.3)')
     &            dmm0,
     &            ( tr(ip,ie1,ie2,it,ioe), ie1 = 1, ne1 )

               end do

             end if

  390     continue

*-----------------------------------------------------------------------
*        e21 axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 4 ) then

           call sel_dc2(m,dc2,ldc2)
           if( ittwo(m) .eq. 4 ) then
             dc2 = '#e2e1numberr.err'
             ldc2 = 16
           else if( ittwo(m) .eq. 5 ) then
             dc2 = 'y/x'
             ldc2 = 3
           end if
*-----------------------------------------------------------------------

           facmx = 1.d0  ! kitamura23/03/31

           do 490 ip = 1, np
           do 490 it = 1, nt

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
     &         ( ( tr(ip,ie1,ie2,it,ioe),
     &             ie2 = 1, ne2 ), ie1 = ne1, 1, -1 )

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 4 ) then

               do ie1 = 1, ne1
               do ie2 = 1, ne2

                 read(jsi,'(1p3e13.4,0pf8.4)')
     &                dmm0, dmm1,
     &                (tr(ip,ie1,ie2,it,k),k=1,2)

               end do
               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 5 ) then

               do ie1 = ne1, 1, -1

                  read(jsi,'(1p1000e11.3)')
     &            dmm0,
     &            ( tr(ip,ie1,ie2,it,ioe), ie2 = 1, ne2 )

               end do

             end if

  490     continue

*-----------------------------------------------------------------------
*        time axis for deposit energy distribution
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 9 ) then

           dc2  = '#t-lowert-upper'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31
           nt_0 = 0
           ne2_0 = 0
           ne1_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_pe1e2t.inc'

           do 590 ie1 = 1, ne1, ne1stepi
           do 590 ie2 = 1, ne2, ne2stepi
           do 590 ip = 1, np, npstepi

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             ntstepi = 1
             do it = 1, nt, ntstepi


               read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &              dmm0, dmm1,
     &      (((((tr(ip+iploop-1,ie1+ie1loop-1,ie2+ie2loop-1,
     &              it+itloop-1,k),k=1,2),
     &      iploop=1,npstepi),ie1loop=1,ne1stepi),ie2loop=1,ne2stepi),
     &      itloop=1,ntstepi)

             end do

! sumover
            it = 1
            call  pdpst2reg_sumover_getput(jsi,m,iax,
     &            npstepi,ne1stepi,ne2stepi,ntstepi,
     &            ip,ie1,ie2,it,nsame,tott_sum)


  590     continue

*-----------------------------------------------------------------------
*        t-e1 axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 5 ) then

           call sel_dc2(m,dc2,ldc2)
           if( ittwo(m) .eq. 4 ) then
             dc2 = '#te1numberr.err'
             ldc2 = 15
           else if( ittwo(m) .eq. 5 ) then
             dc2 = 'y/x'
             ldc2 = 3
           end if
*-----------------------------------------------------------------------

           facmx = 1.d0  ! kitamura23/03/31

           do 690 ip  = 1, np
           do 690 ie2 = 1, ne2

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
     &         ( ( tr(ip,ie1,ie2,it,ioe),
     &             it = 1, nt ), ie1 = ne1, 1, -1 )

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 4 ) then

               do ie1 = 1, ne1
               do it  = 1, nt

                 read(jsi,'(1p3e13.4,0pf8.4)')
     &                dmm0, dmm1,
     &                (tr(ip,ie1,ie2,it,k),k=1,2)

               end do
               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 5 ) then

               do ie1 = ne1, 1, -1

                  read(jsi,'(1p1000e11.3)')
     &            dmm0,
     &            ( tr(ip,ie1,ie2,it,ioe), it = 1, nt )

               end do

             end if

  690     continue

*-----------------------------------------------------------------------
*        e1-t axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 6 ) then

           call sel_dc2(m,dc2,ldc2)
           if( ittwo(m) .eq. 4 ) then
             dc2 = '#e1tnumberr.err'
             ldc2 = 15
           else if( ittwo(m) .eq. 5 ) then
             dc2 = 'y/x'
             ldc2 = 3
           end if
*-----------------------------------------------------------------------

           facmx = 1.d0  ! kitamura23/03/31

           do 790 ip  = 1, np
           do 790 ie2 = 1, ne2

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
     &         ( ( tr(ip,ie1,ie2,it,ioe),
     &             ie1 = 1, ne1 ), it = nt, 1, -1 )

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 4 ) then

               do it  = 1, nt
               do ie1 = 1, ne1

                 read(jsi,'(1p3e13.4,0pf8.4)')
     &                dmm0, dmm1,
     &                (tr(ip,ie1,ie2,it,k),k=1,2)

               end do
               end do
*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 5 ) then

               do it = nt, 1, -1

                  read(jsi,'(1p1000e11.3)')
     &            dmm0,
     &            ( tr(ip,ie1,ie2,it,ioe), ie1 = 1, ne1 )

               end do

             end if

  790     continue

*-----------------------------------------------------------------------
*        t-e2 axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 7 ) then

           call sel_dc2(m,dc2,ldc2)
           if( ittwo(m) .eq. 4 ) then
             dc2 = '#te2numberr.err'
             ldc2 = 15
           else if( ittwo(m) .eq. 5 ) then
             dc2 = 'y/x'
             ldc2 = 3
           end if
*-----------------------------------------------------------------------

           facmx = 1.d0  ! kitamura23/03/31

           do 890 ip  = 1, np
           do 890 ie1 = 1, ne1

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
     &         ( ( tr(ip,ie1,ie2,it,ioe),
     &             it = 1, nt ), ie2 = ne2, 1, -1 )

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 4 ) then

               do ie2 = 1, ne2
               do it  = 1, nt

                 read(jsi,'(1p3e13.4,0pf8.4)')
     &                dmm0, dmm1,
     &                (tr(ip,ie1,ie2,it,k),k=1,2)

               end do
               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 5 ) then

               do ie2 = ne2, 1, -1

                  read(jsi,'(1p1000e11.3)')
     &            dmm0,
     &            ( tr(ip,ie1,ie2,it,ioe), it = 1, nt )

               end do

             end if

  890     continue

*-----------------------------------------------------------------------
*        e2-t axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 8 ) then

           call sel_dc2(m,dc2,ldc2)
           if( ittwo(m) .eq. 4 ) then
             dc2 = '#e2tnumberr.err'
             ldc2 = 15
           else if( ittwo(m) .eq. 5 ) then
             dc2 = 'y/x'
             ldc2 = 3
           end if
*-----------------------------------------------------------------------

           facmx = 1.d0  ! kitamura23/03/31

           do 990 ip  = 1, np
           do 990 ie1 = 1, ne1

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
     &         ( ( tr(ip,ie1,ie2,it,ioe),
     &             ie2 = 1, ne2 ), it = nt, 1, -1 )

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 4 ) then

               do it  = 1, nt
               do ie2 = 1, ne2

                 read(jsi,'(1p3e13.4,0pf8.4)')
     &                dmm0, dmm1,
     &                (tr(ip,ie1,ie2,it,k),k=1,2)

               end do
               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 5 ) then

               do it = nt, 1, -1

                  read(jsi,'(1p1000e11.3)')
     &            dmm0,
     &            ( tr(ip,ie1,ie2,it,ioe), ie2 = 1, ne2 )

               end do

             end if

  990     continue

*-----------------------------------------------------------------------
         end if

  900 continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_deposit2reg(m,
     &                            np,nr,mr,ne1,ne2,nt,kr,eb1,eb2,tb,tr,
     &                            nvl,ivl,rvl)
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
     &                  itout(itlmax), ittwo(itlmax)

        common /tall21/ rtfac(itlmax)

        common /fact01/ facmax(itlmax)  ! kitamura23/03/31

        dimension   kr(mr)
        dimension   eb1(ne1+1)
        dimension   eb2(ne2+1)
        dimension   tb(nt+1)
        dimension   tw(nt)
        dimension   tr(np,0:ne1,0:ne2,nt,2)
        dimension   vl(nr)
        dimension   lr(nr)
        dimension   ivl(nvl)
        dimension   rvl(nvl)

*-----------------------------------------------------------------------
*           itunt(m) = 1 : 1/source
*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

               call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)
               vl_sum = sum(vl(:))

*-----------------------------------------------------------------------
*        ( unit = 1, 2 ; vol = 1.0 )
*-----------------------------------------------------------------------

            if( itunt(m) .eq. 1 ) then

               do ir = 1, nr

                  vl(ir) = 1.0d0

               end do
               vl_sum = 1.0d0

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 2 : /nsec
*-----------------------------------------------------------------------

            if( itunt(m) .eq. 2 ) then

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

            do 100 ie1 = 1, ne1
            do 100 ie2 = 1, ne2
            do 100 it  = 1, nt
            do 100 ip  = 1, np

                  call invert_deposit_stdev(m,A,B,
     &                            tr(ip,ie1,ie2,it,1),
     &                            tr(ip,ie1,ie2,it,2),
     &                vl(1)/abs(rtfac(m)/facmax(m)),ip)

                  tr(ip,ie1,ie2,it,1) = A
                  tr(ip,ie1,ie2,it,2) = B

! sumover
                  fact_in = abs(rtfac(m)/facmax(m))
                  call pdpst2reg_sumover_stdev(1,m,ip,ie1,ie2,it,
     &                1,1,
     &                fact_in,1.0d0,1.0d0,vl(1),1.0d0,vl(1))


  100       continue

      end subroutine


