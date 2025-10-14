
************************************************************************
*                                                                      *
      subroutine check_theat(m,iax,jsn,jsi,dsin,idsi,ill,ilf,ierr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************
      use moddas
      use moddas_tally

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'param01.inc' ! frtati 2021/10/05
      include 'err.inc'

*-----------------------------------------------------------------------

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

      dimension lschn(35), ischn(35)
      character schan(35)*8

      data icsu / 35 /

      data ( schan(i), i = 1, 35 ) /
     &    'mesh    ','part    ','output  ','unit    ','axis    ',
     &    'file    ','title   ','angel   ','2d-type ','factor  ',
     &    'material','x-txt   ','y-txt   ','z-txt   ','gshow   ',
     &    'rshow   ','iechrl  ','electron','volmat  ','epsout  ',
     &    'e-type  ','deposit ','ctmin(1)','ctmax(1)','ctmin(2)',
     &    'ctmax(2)','ctmin(3)','ctmax(3)','resol   ','width   ',
     &    'trcl    ','*trcl   ','gslat   ','err2d'   ,'resfile '/

      data ( lschn(i), i = 1, 35 ) /
     &     4,         4,         6,         4,         4,
     &     4,         5,         5,         7,         6,
     &     8,         5,         5,         5,         5,
     &     5,         6,         8,         6,         6,
     &     6,         7,         8,         8,         8,
     &     8,         8,         8,         5,         5,
     &     4,         5,         5,         5,         7/

*-----------------------------------------------------------------------

cfrtati 2021/10/05 6 -> mxpart
      dimension iptyp(mxpart), ipnkf(mxpart)
      dimension imtyp(mxpart,mxpart), imnkf(mxpart,mxpart)
      dimension jstyp(mxpart), jnkf0(mxpart)

      common /subtra/ isubt, ipsub(mxpart)   ! kitamura22/03/31

      dimension iaxis(6)
      character ifile(6)*100
      dimension lfile(6)


*-----------------------------------------------------------------------

      logical deqn4

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      character dkam*9

      character tname*8
      data      tname /'[t-heat]'/

      dimension icount(9)
      dimension vtrs(13)

      character rglnrf*200, rglnech*200

*-----------------------------------------------------------------------

      character*10 cepn(100)
      dimension lepn(100)
      character*20 ctl(2)

*-----------------------------------------------------------------------

            ierr  = 0
            jpn   = 0

            iout  = 2
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
            idepo  = 0
            igkst  = 0
            idtt   = 0
            ktrs   = 0
            igslt  = 1
            iter2  = 0
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

            ielct  = 0
            imesh  = 0

            ireso = 1
            width = 0.5

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

            if( ipm .ne. 5 .and. ipm .ne. 6 .and.
     &          ipm .ne. 2 .and.
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

         else if( ipm .eq. 11 ) then

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
*        energy mesh
*-----------------------------------------------------------------------

         else if( ipm .eq. 21 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

               ietp = nint( cvvv )

               call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     'e',ietp,ine,emin,emax,edel,isteg)


               goto 150

*-----------------------------------------------------------------------
*        transform
*-----------------------------------------------------------------------

         else if( ipm .eq. 31 .or. ipm .eq. 32 ) then

                  if( ipm .eq. 32 ) ktrs = 1

               call ttrans(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     ic,ktrs,igkst,idtt,vtrs)


               goto 150

*-----------------------------------------------------------------------
*        output
*-----------------------------------------------------------------------

         else if( ipm .eq. 3 ) then

            if( chlw(ic:ic+2) .eq. 'all' ) then

               iout = 1

            else if( chlw(ic:ic+5) .eq. 'simple' ) then

               iout = 2

            else if( chlw(ic:ic+3) .eq. 'heat' ) then

               iout = 3

            else if( chlw(ic:ic+10) .eq. 'deposit_all' .or.
     &               chlw(ic:ic+10) .eq. 'deposit-all' ) then

               iout = 4

            else if( chlw(ic:ic+13) .eq. 'deposit_simple' .or.
     &               chlw(ic:ic+13) .eq. 'deposit-simple' ) then

               iout = 5

            else if( chlw(ic:ic+11) .eq. 'deposit_heat' .or.
     &               chlw(ic:ic+11) .eq. 'deposit-heat' ) then

               iout = 6

            else

               goto 984

            end if

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        unit
*-----------------------------------------------------------------------

         else if( ipm .eq. 4 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               iunt = nint( cvvv )
               if( iunt .lt. 0 .or. iunt .gt. 3 ) goto 988

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        axis
*-----------------------------------------------------------------------

         else if( ipm .eq. 5 ) then

  600       inaxi = inaxi + 1

            if( inaxi .gt. 6 ) goto 991

            if( chlw(ic:ic+2) .eq. 'reg' ) then

               iaxis(inaxi) = 2
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+2) .eq. 'eng' ) then

               iaxis(inaxi) = 1
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
         else if( ipm .eq. 35 ) then

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
*        gshow
*-----------------------------------------------------------------------

         else if( ipm .eq. 15 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               lgshow = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        gslat
*-----------------------------------------------------------------------

         else if( ipm .eq. 33 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               igslt = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        rshow
*-----------------------------------------------------------------------

         else if( ipm .eq. 16 ) then

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

         else if( ipm .eq. 17 ) then

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
*        volmat
*-----------------------------------------------------------------------

         else if( ipm .eq. 19 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               matvol = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        epsout
*-----------------------------------------------------------------------

         else if( ipm .eq. 20 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ieps = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        resolution
*-----------------------------------------------------------------------

         else if( ipm .eq. 29 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ireso = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        width
*-----------------------------------------------------------------------

         else if( ipm .eq. 30 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               width = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        electron heat
*-----------------------------------------------------------------------

         else if( ipm .eq. 18 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ielct = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        deposit
*-----------------------------------------------------------------------

         else if( ipm .eq. 22 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               idepo = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        counter
*-----------------------------------------------------------------------

         else if( ipm .ge. 23 .and. ipm .le. 28 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               icount((ipm-23+2)/2) = 1
               icount( ipm-23+4)    = nint( cvvv )

               if( abs( nint(cvvv) ) .gt. 9999 ) goto 978

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        gshow
*-----------------------------------------------------------------------

         else if( ipm .eq. 34 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               iter2 = nint( cvvv )

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

      call check_unit(m,iunt,iec,cepn,lepn,ierr)

      call check_axis(m,iaxis,inaxi,iec,cepn,lepn,ierr)

      call check_output(m,iout,iec,cepn,lepn,ierr)

      if ( any( itout(m) .eq. (/ 4, 5, 6 /) ) ) then !! output is deposit
        call check_type('e',iec,cepn,lepn,
     &                  itety(m),rtema(m),rtemi(m),itenm(m),
     &                  ietp, emax, emin, ine)
      end if

      call check_2dtype(m,idtyp,iec,cepn,lepn,ierr)

      call check_factor(m,rfact,iec,cepn,lepn,ierr)

      call check_part(m,iptyp,ipnkf,inpat,iec,cepn,lepn,ierr)

*-----------------------------------------------------------------------

      call judge_tall_check(iec,cepn,lepn,dsin(jsn),idsi(jsn),ierr)

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  950    write(dkam,'(i9)') mdas
         m_err = 'Total tally storage number exceeds mdas ='//dkam
         ErrCha = ''
         ErrID = 'L:1052/R:check_theat/F:restheat.f'
         goto 999

  978    m_err = 'counter should be from -9999 to 9999'
         ErrCha = ''
         ErrID = 'L:1057/R:check_theat/F:restheat.f'
         goto 999

  979    m_err = 'Description of material parameter is wrong in tally '
     &            //tname
         ErrCha = ''
         ErrID = 'L:1063/R:check_theat/F:restheat.f'
         goto 999

  982    m_err = 'all particle is not available in tally '//tname
         ErrCha = ''
         ErrID = 'L:1068/R:check_theat/F:restheat.f'
         goto 999

  983    m_err = '2D-type should be 1-7 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1073/R:check_theat/F:restheat.f'
         goto 999

  984    m_err = 'Unknown output parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1078/R:check_theat/F:restheat.f'
         goto 999

  986    write(dkam,'(i5)') itlmax
         m_err = 'Total tally number exceeds itlmax = '//dkam
         ErrCha = ''
         ErrID = 'L:1084/R:check_theat/F:restheat.f'
         goto 999

  987    m_err = 'Unknown parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1089/R:check_theat/F:restheat.f'
         goto 999
  988    m_err = 'Unit should be 0, 1, 2, 3 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1093/R:check_theat/F:restheat.f'
         goto 999

  989    m_err = 'Double definition of the parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1098/R:check_theat/F:restheat.f'
         goto 999

  990    m_err = 'Too many file name in tally '//tname
         ErrCha = ''
         ErrID = 'L:1103/R:check_theat/F:restheat.f'
         goto 999

  991    m_err = 'Too many axis in tally '//tname
         ErrCha = ''
         ErrID = 'L:1108/R:check_theat/F:restheat.f'
         goto 999

  992    m_err = 'Unknown axis name in tally '//tname
         ErrCha = ''
         ErrID = 'L:1113/R:check_theat/F:restheat.f'
         goto 999

  993    m_err = 'In this line, [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:1118/R:check_theat/F:restheat.f'
         goto 999

  994    m_err = 'Name of particle is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1123/R:check_theat/F:restheat.f'
         goto 999

  995    m_err = 'Number of particles is larger than 6 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1128/R:check_theat/F:restheat.f'
         goto 999

  996    m_err = 'Description of mesh is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1133/R:check_theat/F:restheat.f'
         goto 999

  997    m_err = 'Description of parameter is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1138/R:check_theat/F:restheat.f'
         goto 999

  998    m_err = 'Unknown mesh parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1143/R:check_theat/F:restheat.f'
         goto 999

*-----------------------------------------------------------------------

  999 continue

         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1

         call ErrWrite(ErrID, ErrCha)
         write(*,*) 'Error: ' // m_err

*-----------------------------------------------------------------------

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_theat(m,iax,ierr)
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
        common /tall17/ itndy(itlmax)
        common /tall18/ ithet(itlmax)
        common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
        common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

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
        call check_theat(m,iax,jsn(ioe),jsi(ioe),
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
        idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
        idas3 = idas1 + itrgn(m) + ( itrgn(m) + mod(itrgn(m),2) ) / 2
        idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1

*-----------------------------------------------------------------------
*   read tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call read_hetreg(m,iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndy(m),itout(m),
     &      itpan(m),itrgn(m),itrgm(m),
     &      idas_itreg(itreg(m)),itenm(m),das_iterg(iterg(m)),
     &      trRES(irestalm(m)),das(ithet(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        else if ( itmsh(m) .eq. 2 ) then

          call read_hetrz(m,iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndy(m),itout(m),
     &      itpan(m),itrnm(m),itznm(m),
     &      das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &      itenm(m),das_iterg(iterg(m)),
     &      trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 3 ) then

          call read_hetxyz(m,iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndy(m),itout(m),
     &      itpan(m),itmtn(m),ismte(itmtt(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)),
     &      itenm(m),das_iterg(iterg(m)),
     &      trRES(irestalm(m)),das(ithet(m)))

        end if

  800   continue
*-----------------------------------------------------------------------
*   restore tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call restore_hetreg(m,
     &      itndy(m),itout(m),
     &      itpan(m),itrgn(m),itrgm(m),
     &      idas_itreg(itreg(m)),itenm(m),das_iterg(iterg(m)),
     &      trRES(irestalm(m)),das(ithet(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        else if ( itmsh(m) .eq. 2 ) then

          call restore_hetrz(m,
     &      itndy(m),itout(m),
     &      itpan(m),itrnm(m),itznm(m),
     &      das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &      itenm(m),das_iterg(iterg(m)),
     &      trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 3 ) then

          call restore_hetxyz(m,
     &      itndy(m),itout(m),
     &      itpan(m),itmtn(m),ismte(itmtt(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)),
     &      itenm(m),das_iterg(iterg(m)),
     &      trRES(irestalm(m)),das(ithet(m)))

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
      subroutine read_hetreg(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                       nd,ns,np,nr,mr,kr,ne,eb,tr,
     &                       rabs,nvl,ivl,rvl,
     &                       nx,ny,nz,xm,ym,zm)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

*-----------------------------------------------------------------------

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)

        common /tall21/ rtfac(itlmax)
        common /fact01/ facmax(itlmax)
        common /stat / istdev, irestart, ireschk

        dimension   kr(mr)
        dimension   vl(nr)
        dimension   lr(nr)
        dimension   eb(ne+1)
        dimension   tr(nd,nr,0:ne,2)

        dimension   rabs(5)
        dimension   ivl(nvl)
        dimension   rvl(nvl)
        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)

*-----------------------------------------------------------------------

        dimension kl(7,7)
        dimension ks(7)

        data ks / 1, 5, 7, 6, 5, 6, 3/
        data ( kl(1,i),i=1,1) /1/
        data ( kl(2,i),i=1,5) /1,2,3,4,5/
        data ( kl(3,i),i=1,7) /1,6,7,8,9,10,11/
        data ( kl(4,i),i=1,6) /6,12,13,14,15,16/
        data ( kl(5,i),i=1,5) /7,17,18,19,20/
        data ( kl(6,i),i=1,6) /21,22,23,24,25,26/
        data ( kl(7,i),i=1,3) /27,28,29/

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        character chin*200, chlw*200, chcm*200
        character dc2*15

!sumover
        real(8),allocatable :: tott_sum(:,:)


        maxtott = max(7,np+1)
        allocate (tott_sum(maxtott,2))
*-----------------------------------------------------------------------
*        initialization
*-----------------------------------------------------------------------

        if( ns .eq. 3 .or. ns .eq. 6 ) then

           ids = 1
           idf = 1
           idw = 1

        else if( ns .eq. 2 .or. ns .eq. 5 ) then

           ids = 2
           idf = 3
           idw = 7

        else if( ns .eq. 1 .or. ns .eq. 4 ) then

           ids = 2
           idf = 7
           idw = 7

        end if

*-----------------------------------------------------------------------
*        eng axis for deposit energy distributioon
*-----------------------------------------------------------------------
        if( itaxs(m,iax) .eq. 1 ) then

          dc2  = '#e-lowere-upper'
          ldc2 = 15

          facmx = 1.d0  ! kitamura23/03/31

          do 190 ir = 1, nr
          do 190 idd = ids, idf

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            do ie = 1, ne

              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &             ((tr(kl(idd,k),ir,ie,j),j=1,2),k=1,ks(idd))

            end do

! sumover
            read(jsi,'(a)') chin
            read(jsi,'(26x,7(1pe13.4,0pf8.4))')
     &             ((tott_sum(k,j),j=1,2),k=1,ks(idd))

            itype_het = 1
            nsame = ks(idd)
            ie = 1
            call restheatreg_sumover_put(m,iax,
     &           kl,ks,np,itype_het,idd,
     &           ik,ir,ie,maxtott,nsame,tott_sum)

*-----------------------------------------------------------------------
            if( ( idd .eq. 5 .or. idd .eq. 6 ) .and. np .gt. 0 ) then

              if ( idd .eq. 5 ) then
                k0 = 29
                kk = 7
              else
                k0 = 29 + np
                kk = 11
              end if

              call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                          dc2,ldc2,ierr)

              if ( ierr  .ne. 0 ) goto 900
              if ( jpn   .eq. 3 ) goto 900

              do ie = 1, ne

                read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &               dmm0, dmm1,
     &               (tr(kk,ir,ie,j),j=1,2),
     &               ((tr(k0+k,ir,ie,j),j=1,2),k=1,np)

              end do

! sumover
            read(jsi,'(a)') chin
            read(jsi,'(26x,1000(1pe13.4,0pf8.4))')
     &             (tott_sum(1,j),j=1,2),
     &             ((tott_sum(k+1,j),j=1,2),k=1,np)

            if ( idd .eq. 5 ) then
              itype_het = 2
            else
              itype_het = 3
            endif
            nsame = np+1
            ie = 1
            call restheatreg_sumover_put(m,iax,
     &           kl,ks,np,itype_het,idd,
     &           ik,ir,ie,maxtott,nsame,tott_sum)


            end if

  190     continue

*-----------------------------------------------------------------------
*        reg axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 2 ) then

          dc2  = '#numregvolume'
          ldc2 = 13

          facmx = 1.d0  ! kitamura23/03/31

          do 290 idd = ids, idf

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            do ir = 1, nr

              read(jsi,'(i5,1x,i7,1pe13.4,7(1pe13.4,0pf8.4))')
     &             idmm0, idmm1, dmm2,
     &             ((tr(kl(idd,k),ir,0,j),j=1,2),k=1,ks(idd))

            end do

! sumover
            read(jsi,'(a)') chin
            read(jsi,'(26x,7(1pe13.4,0pf8.4))')
     &             ((tott_sum(k,j),j=1,2),k=1,ks(idd))

            itype_het = 1
            nsame = ks(idd)
            ie = 0
            ir = 1
            call restheatreg_sumover_put(m,iax,
     &           kl,ks,np,itype_het,idd,
     &           ik,ir,ie,maxtott,nsame,tott_sum)


*-----------------------------------------------------------------------
            if( ( idd .eq. 5 .or. idd .eq. 6 ) .and. np .gt. 0 ) then

              if ( idd .eq. 5 ) then
                k0 = 29
                kk = 7
              else
                k0 = 29 + np
                kk = 11
              end if

              call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                          dc2,ldc2,ierr)

              if ( ierr  .ne. 0 ) goto 900
              if ( jpn   .eq. 3 ) goto 900

              do ir = 1, nr

                read(jsi,
     &               '(i5,1x,i7,1pe13.4,1000(1pe13.4,0pf8.4))')
     &               idmm0, idmm1, dmm2,
     &               (tr(kk,ir,0,j),j=1,2),
     &                ((tr(k0+k,ir,0,j),j=1,2),k=1,np)

              end do

! sumover
              read(jsi,'(a)') chin
              read(jsi,'(26x,1000(1pe13.4,0pf8.4))')
     &               (tott_sum(1,j),j=1,2),
     &               ((tott_sum(k+1,j),j=1,2),k=1,np)

              if ( idd .eq. 5 ) then
                itype_het = 2
              else
                itype_het = 3
              endif
              nsame = np+1
              ir = 1
              ie = 0
              call restheatreg_sumover_put(m,iax,
     &             kl,ks,np,itype_het,idd,
     &             ik,ir,ie,maxtott,nsame,tott_sum)

            end if

  290     continue

*-----------------------------------------------------------------------
         end if

  900 continue
      
      deallocate (tott_sum)

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_hetrz(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                      nd,ns,np,nr,nz,rm,zm,ne,eb,tr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

*-----------------------------------------------------------------------

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

        dimension   rm(nr+1)
        dimension   zm(nz+1)
        dimension   tr(nd,nr,nz,0:ne,2)

        dimension   rabs(5)
        dimension   eb(ne+1)

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        dimension kl(7,7)
        dimension ks(7)

        data ks / 1, 5, 7, 6, 5, 6, 3/
        data ( kl(1,i),i=1,1) /1/
        data ( kl(2,i),i=1,5) /1,2,3,4,5/
        data ( kl(3,i),i=1,7) /1,6,7,8,9,10,11/
        data ( kl(4,i),i=1,6) /6,12,13,14,15,16/
        data ( kl(5,i),i=1,5) /7,17,18,19,20/
        data ( kl(6,i),i=1,6) /21,22,23,24,25,26/
        data ( kl(7,i),i=1,3) /27,28,29/

        dimension ict(7)
        data ict /1,6,7,8,9,10,11/

*-----------------------------------------------------------------------

        character chin*200, chlw*200, chcm*200
        character dc2*15

!sumover
        real(8),allocatable :: tott_sum(:,:)

        maxtott = max(7,np+1)
        allocate (tott_sum(maxtott,2))

*-----------------------------------------------------------------------
*        initialization
*-----------------------------------------------------------------------

        if( ns .eq. 3 .or. ns .eq. 6 ) then

           ids = 1
           idf = 1
           idw = 1

        else if( ns .eq. 2 .or. ns .eq. 5 ) then

           ids = 2
           idf = 3
           idw = 7

        else if( ns .eq. 1 .or. ns .eq. 4 ) then

           ids = 2
           idf = 7
           idw = 7

        end if

*-----------------------------------------------------------------------
*        eng axis for deposit energy distributioon
*-----------------------------------------------------------------------
        if( itaxs(m,iax) .eq. 1 ) then

          dc2  = '#e-lowere-upper'
          ldc2 = 15

          facmx = 1.d0  ! kitamura23/03/31

          do 190 ir = 1, nr
          do 190 iz = 1, nz
          do 190 idd = ids, idf

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            do ie = 1, ne

              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &              dmm0, dmm1,
     &              ((tr(kl(idd,k),ir,iz,ie,j),j=1,2),k=1,ks(idd))

            end do

! sumover
            read(jsi,'(a)') chin
            read(jsi,'(26x,7(1pe13.4,0pf8.4))')
     &             ((tott_sum(k,j),j=1,2),k=1,ks(idd))

            itype_het = 1
            nsame = ks(idd)
            ie = 1
            call restheatrz_sumover_put(m,iax,
     &           kl,ks,np,itype_het,idd,
     &           ik,ir,iz,ie,maxtott,nsame,tott_sum)


*-----------------------------------------------------------------------
            if( ( idd .eq. 5 .or. idd .eq. 6 ) .and. np .gt. 0 ) then

              if ( idd .eq. 5 ) then
                k0 = 29
                kk = 7
              else
                k0 = 29 + np
                kk = 11
              end if

              call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                          dc2,ldc2,ierr)

              if ( ierr  .ne. 0 ) goto 900
              if ( jpn   .eq. 3 ) goto 900

              do ie = 1, ne

                read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &               dmm0, dmm1,
     &               (tr(kk,ir,iz,ie,j),j=1,2),
     &               ((tr(k0+k,ir,iz,ie,j),j=1,2),k=1,np)

              end do

! sumover
            read(jsi,'(a)') chin
            read(jsi,'(26x,7(1pe13.4,0pf8.4))')
     &             (tott_sum(1,j),j=1,2),
     &             ((tott_sum(k+1,j),j=1,2),k=1,np)

            if ( idd .eq. 5 ) then
              itype_het = 2
            else
              itype_het = 3
            endif
            nsame = np+1
            ie = 1
            call restheatrz_sumover_put(m,iax,
     &           kl,ks,np,itype_het,idd,
     &           ik,ir,iz,ie,maxtott,nsame,tott_sum)

            end if

  190     continue

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 5 ) then

          dc2  = '#z-lowerz-upper'
          ldc2 = 15

          facmx = 1.d0  ! kitamura23/03/31

          do 290 ir = 1, nr
          do 290 idd = ids, idf

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            do iz = 1, nz

              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &              dmm0, dmm1,
     &              ((tr(kl(idd,k),ir,iz,0,j),j=1,2),k=1,ks(idd))

            end do

! sumover
            read(jsi,'(a)') chin
            read(jsi,'(26x,7(1pe13.4,0pf8.4))')
     &             ((tott_sum(k,j),j=1,2),k=1,ks(idd))

            itype_het = 1
            nsame = ks(idd)
            iz = 1
            ie = 0
            call restheatrz_sumover_put(m,iax,
     &           kl,ks,np,itype_het,idd,
     &           ik,ir,iz,ie,maxtott,nsame,tott_sum)


*-----------------------------------------------------------------------
            if( ( idd .eq. 5 .or. idd .eq. 6 ) .and. np .gt. 0 ) then

              if ( idd .eq. 5 ) then
                k0 = 29
                kk = 7
              else
                k0 = 29 + np
                kk = 11
              end if

              call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                          dc2,ldc2,ierr)

              if ( ierr  .ne. 0 ) goto 900
              if ( jpn   .eq. 3 ) goto 900

              do iz = 1, nz

                read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &               dmm0, dmm1,
     &               (tr(kk,ir,iz,0,j),j=1,2),
     &               ((tr(k0+k,ir,iz,0,j),j=1,2),k=1,np)

              end do

! sumover
            read(jsi,'(a)') chin
            read(jsi,'(26x,7(1pe13.4,0pf8.4))')
     &             (tott_sum(1,j),j=1,2),
     &             ((tott_sum(k+1,j),j=1,2),k=1,np)

            if ( idd .eq. 5 ) then
              itype_het = 2
            else
              itype_het = 3
            endif
            nsame = np+1
            iz = 1
            ie = 0
            call restheatrz_sumover_put(m,iax,
     &           kl,ks,np,itype_het,idd,
     &           ik,ir,iz,ie,maxtott,nsame,tott_sum)

            end if

  290     continue

*-----------------------------------------------------------------------
*        r axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 6 ) then

          dc2  = '#r-lowerr-upper'
          ldc2 = 15

          facmx = 1.d0  ! kitamura23/03/31

          do 390 iz = 1, nz
          do 390 idd = ids, idf

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            do ir = 1, nr

              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &              dmm0, dmm1,
     &              ((tr(kl(idd,k),ir,iz,0,j),j=1,2),k=1,ks(idd))

            end do

! sumover
            read(jsi,'(a)') chin
            read(jsi,'(26x,7(1pe13.4,0pf8.4))')
     &             ((tott_sum(k,j),j=1,2),k=1,ks(idd))

            itype_het = 1
            nsame = ks(idd)
            ir = 1
            ie = 0
            call restheatrz_sumover_put(m,iax,
     &           kl,ks,np,itype_het,idd,
     &           ik,ir,iz,ie,maxtott,nsame,tott_sum)


*-----------------------------------------------------------------------
            if( ( idd .eq. 5 .or. idd .eq. 6 ) .and. np .gt. 0 ) then

              if ( idd .eq. 5 ) then
                k0 = 29
                kk = 7
              else
                k0 = 29 + np
                kk = 11
              end if

              call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                          dc2,ldc2,ierr)

              if ( ierr  .ne. 0 ) goto 900
              if ( jpn   .eq. 3 ) goto 900

              do ir = 1, nr

                read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &               dmm0, dmm1,
     &               (tr(kk,ir,iz,0,j),j=1,2),
     &               ((tr(k0+k,ir,iz,0,j),j=1,2),k=1,np)

              end do

! sumover
            read(jsi,'(a)') chin
            read(jsi,'(26x,7(1pe13.4,0pf8.4))')
     &             (tott_sum(1,j),j=1,2),
     &             ((tott_sum(k+1,j),j=1,2),k=1,np)

            if ( idd .eq. 5 ) then
              itype_het = 2
            else
              itype_het = 3
            endif
            nsame = np+1
            ir = 1
            ie = 0
            call restheatrz_sumover_put(m,iax,
     &           kl,ks,np,itype_het,idd,
     &           ik,ir,iz,ie,maxtott,nsame,tott_sum)


            end if

  390     continue

*-----------------------------------------------------------------------
*        rz axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 10 ) then

          call sel_dc2(m,dc2,ldc2)
          if( ittwo(m) .eq. 4 ) then
            dc2 = '#rzheatr.err'
            ldc2 = 12
          else if( ittwo(m) .eq. 5 ) then
            dc2 = 'r/z'
            ldc2 = 3
          end if

*-----------------------------------------------------------------------

          facmx = 1.d0  ! kitamura23/03/31

          do 490 ip = 1, idw

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
            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               read(jsi,'(1p10e11.3)')
     &         ( ( tr(ict(ip),ir,iz,0,ioe), iz = 1, nz ), ir = nr,1,-1 )

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 4 ) then

              do iz = 1, nz
              do ir = 1, nr

                read(jsi,'(1p3e11.3,0pf8.4)')
     &             dmm0, dmm1,
     &             tr(ict(ip),ir,iz,0,1), tr(ict(ip),ir,iz,0,2)

              end do
              end do

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 5 ) then

              do ir = nr, 1, -1

                read(jsi,'(1p1000e11.3)')
     &            dmm0,
     &            ( tr(ict(ip),ir,iz,0,ioe), iz = 1, nz )

              end do

            end if

  490     continue

*-----------------------------------------------------------------------
         end if

  900 continue

      deallocate (tott_sum)

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_hetxyz(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                       nd,ns,np,nl,lt,nx,ny,nz,xm,ym,zm,ne,eb,
     &                       tr,rabs)
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

*-----------------------------------------------------------------------

        dimension   lt(nl)
        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)
        dimension   tr(nd,nx,ny,nz,0:ne,2)

        dimension   rabs(5)
        dimension   eb(ne+1)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        dimension kl(7,7)
        dimension ks(7)

        data ks / 1, 5, 7, 6, 5, 6, 3/
        data ( kl(1,i),i=1,1) /1/
        data ( kl(2,i),i=1,5) /1,2,3,4,5/
        data ( kl(3,i),i=1,7) /1,6,7,8,9,10,11/
        data ( kl(4,i),i=1,6) /6,12,13,14,15,16/
        data ( kl(5,i),i=1,5) /7,17,18,19,20/
        data ( kl(6,i),i=1,6) /21,22,23,24,25,26/
        data ( kl(7,i),i=1,3) /27,28,29/

        dimension ict(7)
        data ict /1,6,7,8,9,10,11/

*-----------------------------------------------------------------------

        character chin*200, chlw*200, chcm*200
        character dc2*15

!sumover
        real(8),allocatable :: tott_sum(:,:)

        maxtott = max(7,np+1)
        allocate (tott_sum(maxtott,2))

*-----------------------------------------------------------------------
*        initialization
*-----------------------------------------------------------------------

        if( ns .eq. 3 .or. ns .eq. 6 ) then

           ids = 1
           idf = 1
           idw = 1

        else if( ns .eq. 2 .or. ns .eq. 5 ) then

           ids = 2
           idf = 3
           idw = 7

        else if( ns .eq. 1 .or. ns .eq. 4 ) then

           ids = 2
           idf = 7
           idw = 7

        end if

*-----------------------------------------------------------------------
*        eng axis for deposit energy distributioon
*-----------------------------------------------------------------------
        if( itaxs(m,iax) .eq. 1 ) then

          dc2  = '#e-lowere-upper'
          ldc2 = 15

          facmx = 1.d0  ! kitamura23/03/31

          do 190 ix = 1, nx
          do 190 iy = 1, ny
          do 190 iz = 1, nz
          do 190 idd = ids, idf

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            do ie = 1, ne

              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &             ((tr(kl(idd,k),ix,iy,iz,ie,j),j=1,2),k=1,ks(idd))

            end do

! sumover
            read(jsi,'(a)') chin
            read(jsi,'(26x,7(1pe13.4,0pf8.4))')
     &             ((tott_sum(k,j),j=1,2),k=1,ks(idd))

            itype_het = 1
            nsame = ks(idd)
            ie = 1
            call restheatxyz_sumover_put(m,iax,
     &           kl,ks,np,itype_het,idd,
     &           ik,ix,iy,iz,ie,maxtott,nsame,tott_sum)

*-----------------------------------------------------------------------
            if( ( idd .eq. 5 .or. idd .eq. 6 ) .and. np .gt. 0 ) then

              if ( idd .eq. 5 ) then
                k0 = 29
                kk = 7
              else
                k0 = 29 + np
                kk = 11
              end if

              call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                          dc2,ldc2,ierr)

              if ( ierr  .ne. 0 ) goto 900
              if ( jpn   .eq. 3 ) goto 900

              do ie = 1, ne

                read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &               dmm0, dmm1,
     &               (tr(kk,ix,iy,iz,ie,j),j=1,2),
     &               ((tr(k0+np+k,ix,iy,iz,ie,j),j=1,2),k=1,np)

              end do

! sumover
            read(jsi,'(a)') chin
            read(jsi,'(26x,7(1pe13.4,0pf8.4))')
     &             (tott_sum(1,j),j=1,2),
     &             ((tott_sum(k+1,j),j=1,2),k=1,np)

            if ( idd .eq. 5 ) then
              itype_het = 2
            else
              itype_het = 3
            endif
            nsame = np+1
            ie = 1
            call restheatxyz_sumover_put(m,iax,
     &           kl,ks,np,itype_het,idd,
     &           ik,ix,iy,iz,ie,maxtott,nsame,tott_sum)


            end if

  190     continue

*-----------------------------------------------------------------------
*        x axis
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 3 ) then

          dc2  = '#x-lowerx-upper'
          ldc2 = 15

          facmx = 1.d0  ! kitamura23/03/31

          do 290 iy = 1, ny
          do 290 iz = 1, nz
          do 290 idd = ids, idf

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            do ix = 1, nx

              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &             ((tr(kl(idd,k),ix,iy,iz,0,j),j=1,2),k=1,ks(idd))

            end do

! sumover
            read(jsi,'(a)') chin
            read(jsi,'(26x,7(1pe13.4,0pf8.4))')
     &             ((tott_sum(k,j),j=1,2),k=1,ks(idd))

            itype_het = 1
            nsame = ks(idd)
            ix = 1
            ie = 0
            call restheatxyz_sumover_put(m,iax,
     &           kl,ks,np,itype_het,idd,
     &           ik,ix,iy,iz,ie,maxtott,nsame,tott_sum)


*-----------------------------------------------------------------------
            if( ( idd .eq. 5 .or. idd .eq. 6 ) .and. np .gt. 0 ) then

              if ( idd .eq. 5 ) then
                k0 = 29
                kk = 7
              else
                k0 = 29 + np
                kk = 11
              end if

              call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                          dc2,ldc2,ierr)

              if ( ierr  .ne. 0 ) goto 900
              if ( jpn   .eq. 3 ) goto 900

              do ix = 1, nx

                read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &               dmm0, dmm1,
     &               (tr(kk,ix,iy,iz,0,j),j=1,2),
     &               ((tr(k0+np+k,ix,iy,iz,0,j),j=1,2),k=1,np)

              end do

! sumover
            read(jsi,'(a)') chin
            read(jsi,'(26x,7(1pe13.4,0pf8.4))')
     &             (tott_sum(1,j),j=1,2),
     &             ((tott_sum(k+1,j),j=1,2),k=1,np)

            if ( idd .eq. 5 ) then
              itype_het = 2
            else
              itype_het = 3
            endif
            nsame = np+1
            ix = 1
            ie = 0
            call restheatxyz_sumover_put(m,iax,
     &           kl,ks,np,itype_het,idd,
     &           ik,ix,iy,iz,ie,maxtott,nsame,tott_sum)


            end if

  290     continue

*-----------------------------------------------------------------------
*        y axis
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 4 ) then

          dc2  = '#y-lowery-upper'
          ldc2 = 15

          facmx = 1.d0  ! kitamura23/03/31

          do 390 ix = 1, nx
          do 390 iz = 1, nz
          do 390 idd = ids, idf

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            do iy = 1, ny

              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &             ((tr(kl(idd,k),ix,iy,iz,0,j),j=1,2),k=1,ks(idd))

            end do

! sumover
            read(jsi,'(a)') chin
            read(jsi,'(26x,7(1pe13.4,0pf8.4))')
     &             ((tott_sum(k,j),j=1,2),k=1,ks(idd))

            itype_het = 1
            nsame = ks(idd)
            iy = 1
            ie = 0
            call restheatxyz_sumover_put(m,iax,
     &           kl,ks,np,itype_het,idd,
     &           ik,ix,iy,iz,ie,maxtott,nsame,tott_sum)


*-----------------------------------------------------------------------
            if( ( idd .eq. 5 .or. idd .eq. 6 ) .and. np .gt. 0 ) then

              if ( idd .eq. 5 ) then
                k0 = 29
                kk = 7
              else
                k0 = 29 + np
                kk = 11
              end if

              call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                          dc2,ldc2,ierr)

              if ( ierr  .ne. 0 ) goto 900
              if ( jpn   .eq. 3 ) goto 900

              do iy = 1, ny

                read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &               dmm0, dmm1,
     &               (tr(kk,ix,iy,iz,0,j),j=1,2),
     &               ((tr(k0+np+k,ix,iy,iz,0,j),j=1,2),k=1,np)

              end do

! sumover
            read(jsi,'(a)') chin
            read(jsi,'(26x,7(1pe13.4,0pf8.4))')
     &             (tott_sum(1,j),j=1,2),
     &             ((tott_sum(k+1,j),j=1,2),k=1,np)

            if ( idd .eq. 5 ) then
              itype_het = 2
            else
              itype_het = 3
            endif
            nsame = np+1
            iy = 1
            ie = 0
            call restheatxyz_sumover_put(m,iax,
     &           kl,ks,np,itype_het,idd,
     &           ik,ix,iy,iz,ie,maxtott,nsame,tott_sum)


            end if

  390     continue

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 5 ) then

          dc2  = '#z-lowerz-upper'
          ldc2 = 15

          facmx = 1.d0  ! kitamura23/03/31

          do 490 ix = 1, nx
          do 490 iy = 1, ny
          do 490 idd = ids, idf

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            do iz = 1, nz

              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &             ((tr(kl(idd,k),ix,iy,iz,0,j),j=1,2),k=1,ks(idd))

            end do

! sumover
            read(jsi,'(a)') chin
            read(jsi,'(26x,7(1pe13.4,0pf8.4))')
     &             ((tott_sum(k,j),j=1,2),k=1,ks(idd))

            itype_het = 1
            nsame = ks(idd)
            iz = 1
            ie = 0
            call restheatxyz_sumover_put(m,iax,
     &           kl,ks,np,itype_het,idd,
     &           ik,ix,iy,iz,ie,maxtott,nsame,tott_sum)


*-----------------------------------------------------------------------
            if( ( idd .eq. 5 .or. idd .eq. 6 ) .and. np .gt. 0 ) then

              if ( idd .eq. 5 ) then
                k0 = 29
                kk = 7
              else
                k0 = 29 + np
                kk = 11
              end if

              call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                          dc2,ldc2,ierr)

              if ( ierr  .ne. 0 ) goto 900
              if ( jpn   .eq. 3 ) goto 900

              do iz = 1, nz

                read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &               dmm0, dmm1,
     &               (tr(kk,ix,iy,iz,0,j),j=1,2),
     &               ((tr(k0+np+k,ix,iy,iz,0,j),j=1,2),k=1,np)

              end do

! sumover
            read(jsi,'(a)') chin
            read(jsi,'(26x,7(1pe13.4,0pf8.4))')
     &             (tott_sum(1,j),j=1,2),
     &             ((tott_sum(k+1,j),j=1,2),k=1,np)

            if ( idd .eq. 5 ) then
              itype_het = 2
            else
              itype_het = 3
            endif
            nsame = np+1
            iz = 1
            ie = 0
            call restheatxyz_sumover_put(m,iax,
     &           kl,ks,np,itype_het,idd,
     &           ik,ix,iy,iz,ie,maxtott,nsame,tott_sum)


            end if

  490     continue

*-----------------------------------------------------------------------
*        xy axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 7 ) then

          call sel_dc2(m,dc2,ldc2)
          if( ittwo(m) .eq. 4 ) then
            dc2 = '#xyheatr.err'
            ldc2 = 12
          else if( ittwo(m) .eq. 5 ) then
            dc2 = 'y/x'
            ldc2 = 3
          end if

*-----------------------------------------------------------------------

          facmx = 1.d0  ! kitamura23/03/31

          do 590 iz = 1, nz
          do 590 ip = 1, idw

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
            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

              read(jsi,'(1p10e11.3)')
     &        ( ( tr(ict(ip),ix,iy,iz,0,ioe),
     &                        ix = 1, nx ), iy = ny, 1, -1 )

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 4 ) then

              do iy = 1, ny
              do ix = 1, nx

                read(jsi,'(1p3e11.3,0pf8.4)')
     &             dmm0, dmm1,
     &             tr(ict(ip),ix,iy,iz,0,1), tr(ict(ip),ix,iy,iz,0,2)

              end do
              end do

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 5 ) then

              do iy = ny, 1, -1

                read(jsi,'(1p1000e11.3)')
     &          dmm0,
     &          ( tr(ict(ip),ix,iy,iz,0,ioe), ix = 1, nx )

              end do

            end if

  590     continue

*-----------------------------------------------------------------------
*        yz axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 8 ) then

          call sel_dc2(m,dc2,ldc2)
          if( ittwo(m) .eq. 4 ) then
            dc2 = '#yzheatr.err'
            ldc2 = 12
          else if( ittwo(m) .eq. 5 ) then
            dc2 = 'y/z'
            ldc2 = 3
          end if

*-----------------------------------------------------------------------

          facmx = 1.d0  ! kitamura23/03/31

          do 690 ix = 1, nx
          do 690 ip = 1, idw

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
            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

              read(jsi,'(1p10e11.3)')
     &        ( ( tr(ict(ip),ix,iy,iz,0,ioe),
     &                        iz = 1, nz ), iy = ny, 1, -1 )

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 4 ) then

              do iz = 1, nz
              do iy = 1, ny

                read(jsi,'(1p3e11.3,0pf8.4)')
     &             dmm0, dmm1,
     &             tr(ict(ip),ix,iy,iz,0,1), tr(ict(ip),ix,iy,iz,0,2)

              end do
              end do

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 5 ) then

              do iy = ny, 1, -1

                read(jsi,'(1p1000e11.3)')
     &          dmm0,
     &          ( tr(ict(ip),ix,iy,iz,0,ioe), iz = 1, nz )

              end do

            end if

  690     continue

*-----------------------------------------------------------------------
*        xz axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 9 ) then

          call sel_dc2(m,dc2,ldc2)
          if( ittwo(m) .eq. 4 ) then
            dc2 = '#xzheatr.err'
            ldc2 = 12
          else if( ittwo(m) .eq. 5 ) then
            dc2 = 'x/z'
            ldc2 = 3
          end if

*-----------------------------------------------------------------------

          facmx = 1.d0  ! kitamura23/03/31

          do 790 iy = 1, ny
          do 790 ip = 1, idw

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
            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

              read(jsi,'(1p10e11.3)')
     &        ( ( tr(ict(ip),ix,iy,iz,0,ioe),
     &                       iz = 1, nz ), ix = nx, 1, -1 )

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 4 ) then

              do iz = 1, nz
              do ix = 1, nx

                read(jsi,'(1p3e11.3,0pf8.4)')
     &             dmm0, dmm1,
     &             tr(ict(ip),ix,iy,iz,0,1), tr(ict(ip),ix,iy,iz,0,2)

              end do
              end do

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 5 ) then

              do ix = nx, 1, -1

                read(jsi,'(1p1000e11.3)')
     &          dmm0,
     &          ( tr(ict(ip),ix,iy,iz,0,ioe), iz = 1, nz )

              end do

            end if

  790     continue

*-----------------------------------------------------------------------
         end if

  900 continue

      deallocate (tott_sum)

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_hetreg(m,
     &                          nd,ns,np,nr,mr,kr,ne,eb,tr,
     &                          rabs,nvl,ivl,rvl,
     &                          nx,ny,nz,xm,ym,zm)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        dimension   kr(mr)
        dimension   vl(nr)
        dimension   lr(nr)
        dimension   eb(ne+1)
        dimension   tr(nd,nr,0:ne,2)
        dimension   rabs(5)
        dimension   ivl(nvl)
        dimension   rvl(nvl)
        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)
        real(8) cfac
        real(8),parameter:: c2gy=1.602d-10
        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)

        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)

        common /tall21/ rtfac(itlmax)
        common /fact01/ facmax(itlmax)

*-----------------------------------------------------------------------
*           itunt(m) = 0 ; cfac = c2gy
*-----------------------------------------------------------------------
        if( itunt(m) .eq. 0)then
          cfac = c2gy
        else
          cfac = 1.0d0
        endif

*-----------------------------------------------------------------------

        call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)

*-----------------------------------------------------------------------

        if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/3

        vl_sum = 0.0d0
        do ir=1,nr
          vl_sum = vl_sum + vl(ir)
        enddo

        do 100 ir = 1, nr
        do 100 ik = 1, nd
        do 100 ie = 0, ne

          call invert_heat_stdev(m,A,B,
     &                           tr(ik,ir,ie,1),
     &                           tr(ik,ir,ie,2),
     &    vl(ir)/cfac/abs(rtfac(m)/facmax(m)))

          tr(ik,ir,ie,1) = A
          tr(ik,ir,ie,2) = B

! sumover
          do iax = 1, itaxn(m)

            if(itaxs(m,iax) == 1 .and. (ie == 0 .or. ie == 1)) then  ! energ
              fact_sum = vl(ir)/cfac/abs(rtfac(m)/facmax(m))
              call restore_heatreg_sumover(m,iax,
     &             ik,ir,ie,fact_sum,nd,nr,1)           ! nr = 1
            else if(itaxs(m,iax) == 2 .and. ir == 1) then  ! reg
              fact_sum = vl_sum/cfac/abs(rtfac(m)/facmax(m))
              call restore_heatreg_sumover(m,iax,
     &             ik,ir,ie,fact_sum,nd,1,ne)           ! nr = 1
            endif
          enddo

  100 continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_hetrz(m,
     &                         nd,ns,np,nr,nz,rm,zm,ne,eb,tr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        dimension   rm(nr+1)
        dimension   zm(nz+1)
        dimension   tr(nd,nr,nz,0:ne,2)
        dimension   rabs(5)
        dimension   eb(ne+1)
        real(8) cfac
        real(8),parameter:: c2gy=1.602d-10
        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)

        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)

        common /tall21/ rtfac(itlmax)
        common /fact01/ facmax(itlmax)

        real(8),allocatable :: vl_z(:),vl_r(:)

*-----------------------------------------------------------------------

        vl(ir,iz) = pi * ( rm(ir+1)**2 - rm(ir)**2 )
     &                 * ( zm(iz+1) - zm(iz) )

        pi = 2.0 * asin(1.0)

        allocate (vl_z(nr),vl_r(nz))
        vl_z(:) =0.0d0
        vl_r(:) = 0.0d0

        do ir=1,nr
          do iz=1,nz
             vl_z(ir) = vl_z(ir) + vl(ir,iz)
             vl_r(iz) = vl_r(iz) + vl(ir,iz)
          enddo
        enddo 

*-----------------------------------------------------------------------
*           itunt(m) = 0 ; cfac = c2gy
*-----------------------------------------------------------------------
        if( itunt(m) .eq. 0)then
          cfac = c2gy
        else
          cfac = 1.0d0
        endif

*-----------------------------------------------------------------------

        if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

        do 100 ir = 1, nr
        do 100 iz = 1, nz
        do 100 ik = 1, nd
        do 100 ie = 0, ne
          call invert_heat_stdev(m,A,B,
     &                           tr(ik,ir,iz,ie,1),
     &                           tr(ik,ir,iz,ie,2),
     &        vl(ir,iz)/cfac/abs(rtfac(m)/facmax(m)))

          tr(ik,ir,iz,ie,1) = A
          tr(ik,ir,iz,ie,2) = B

! sumover
          do iax = 1, itaxn(m)

            if(itaxs(m,iax) == 1 .and. (ie == 0 .or. ie == 1)) then  ! energ
              fact_sum = vl(ir,iz)/cfac/abs(rtfac(m)/facmax(m))
              call restore_heatrz_sumover(m,iax,
     &             ik,ir,iz,ie,fact_sum,nd,nr,nz,1)           ! ne = 1
            else if(itaxs(m,iax) == 5 .and. iz == 1) then  ! reg
              fact_sum = vl_z(ir)/cfac/abs(rtfac(m)/facmax(m))
              call restore_heatrz_sumover(m,iax,
     &             ik,ir,iz,ie,fact_sum,nd,nr,1,ne)           ! nz = 1
            else if(itaxs(m,iax) == 6 .and. ir == 1) then  ! reg
              fact_sum = vl_r(iz)/cfac/abs(rtfac(m)/facmax(m))
              call restore_heatrz_sumover(m,iax,
     &             ik,ir,iz,ie,fact_sum,nd,1,nz,ne)           ! nr = 1
            endif
          enddo

  100 continue
      deallocate (vl_z,vl_r)

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_hetxyz(m,
     &                          nd,ns,np,nl,lt,nx,ny,nz,xm,ym,zm,ne,eb,
     &                          tr,rabs)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
        common /tall34/ itvm(itlmax)
        common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

        common /tall21/ rtfac(itlmax)
        common /fact01/ facmax(itlmax)
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)

*-----------------------------------------------------------------------

        dimension   lt(nl)
        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)
        dimension   tr(nd,nx,ny,nz,0:ne,2)
        dimension   rabs(5)
        dimension   eb(ne+1)
        real(8) cfac
        real(8),parameter:: c2gy=1.602d-10
        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)

        real(8),allocatable :: vl_x(:,:),vl_y(:,:),vl_z(:,:)

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

         vl(ix,iy,iz) = vls(nl,lt,itmcn(m),itvm(m),itmtr(m,4),
     &                      xm(ix),xm(ix+1),
     &                      ym(iy),ym(iy+1),
     &                      zm(iz),zm(iz+1))

        allocate (vl_x(ny,nz),vl_y(nx,nz),vl_z(nx,ny))
        vl_x(:,:) = 0.0d0
        vl_y(:,:) = 0.0d0
        vl_z(:,:) = 0.0d0
        do iz=1,nz
          do iy=1,ny
            do ix=1,nx
              vl_x(iy,iz) = vl_x(iy,iz) + vl(ix,iy,iz)
              vl_y(ix,iz) = vl_y(ix,iz) + vl(ix,iy,iz)
              vl_z(ix,iy) = vl_z(ix,iy) + vl(ix,iy,iz)
            enddo
          enddo
        enddo 


*-----------------------------------------------------------------------
*           itunt(m) = 0 ; cfac = c2gy
*-----------------------------------------------------------------------
        if( itunt(m) .eq. 0)then
          cfac = c2gy
        else
          cfac = 1.0d0
        endif

*-----------------------------------------------------------------------

        if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

        do 100 ix = 1, nx
        do 100 iy = 1, ny
        do 100 iz = 1, nz
        do 100 ik = 1, nd
        do 100 ie = 0, ne
          call invert_heat_stdev(m,A,B,
     &                           tr(ik,ix,iy,iz,ie,1),
     &                           tr(ik,ix,iy,iz,ie,2),
     &      vl(ix,iy,iz)/cfac/abs(rtfac(m)/facmax(m)))

          tr(ik,ix,iy,iz,ie,1) = A
          tr(ik,ix,iy,iz,ie,2) = B

! sumover
          do iax = 1, itaxn(m)

            if(itaxs(m,iax) == 1 .and. (ie == 0 .or. ie == 1)) then  ! energ
              fact_sum = vl(ix,iy,iz)/cfac/abs(rtfac(m)/facmax(m))
              call restore_heatxyz_sumover(m,iax,
     &             ik,ix,iy,iz,ie,fact_sum,nd,nx,ny,nz,1)           ! ne = 1
            else if(itaxs(m,iax) == 3 .and. ix == 1) then  ! reg
              fact_sum = vl_x(iy,iz)/cfac/abs(rtfac(m)/facmax(m))
              call restore_heatxyz_sumover(m,iax,
     &             ik,ix,iy,iz,ie,fact_sum,nd,1,ny,nx,ne)           ! nx = 1
            else if(itaxs(m,iax) == 4 .and. iy == 1) then  ! reg
              fact_sum = vl_y(ix,iz)/cfac/abs(rtfac(m)/facmax(m))
              call restore_heatxyz_sumover(m,iax,
     &             ik,ix,iy,iz,ie,fact_sum,nd,nx,1,nz,ne)           ! ny = 1
            else if(itaxs(m,iax) == 5 .and. iz == 1) then  ! reg
              fact_sum = vl_z(ix,iy)/cfac/abs(rtfac(m)/facmax(m))
              call restore_heatxyz_sumover(m,iax,
     &             ik,ix,iy,iz,ie,fact_sum,nd,nx,ny,1,ne)           ! nz = 1
            endif
          enddo

  100 continue

        deallocate (vl_x,vl_y,vl_z)

      end subroutine


***********************************************************************
*                                                                     *
* sumover subroutine group                                            *
*                                                                     *
***********************************************************************

***********************************************************************
*                                                                     *
* heatreg                                                             *
*                                                                     *
***********************************************************************

***********************************************************************
*                                                                     *
      subroutine  restheatreg_sumover_put(m,iax,
     &    kl,ks,np,itype_het,idd,
     &    ik,ir,ie,maxtott,nsame,tott_sum)
*                                                                     *
***********************************************************************

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

      dimension kl(7,7)
      dimension ks(7)

      real(8),pointer :: tr_sum(:)


      real(8) :: tott_sum(maxtott,2)

      tr_sum => trRES_sum(irestalm_sum(m,iax):)

      call restheatreg_sumover_put_sub(
     &    kl,ks,np,itype_het,idd,
     &    ik,ir,ie,maxtott,nsame,tott_sum,
     &    itndy_sum(m,iax),itrgn_sum(m,iax),itenm_sum(m,iax),
     &    tr_sum)

      return
      end

************************************************************************
*                                                                      *
      subroutine restheatreg_sumover_put_sub(
     &    kl,ks,np,itype_het,idd,
     &    ik,ir,ie,maxtott,nsame,tott_sum,
     &    nd_sum,nr_sum,ne_sum,tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      dimension kl(7,7)
      dimension ks(7)

      real(8) :: tr_sum(nd_sum,nr_sum,0:ne_sum,2)

      real(8) :: tott_sum(maxtott,2)

      if(itype_het == 1) then
        do k = 1, ks(idd)
           tr_sum(kl(idd,k),ir,ie,1) = tott_sum(k,1)
           tr_sum(kl(idd,k),ir,ie,2) = tott_sum(k,2)
        end do
      else if(itype_het == 2) then
        tr_sum(7,ir,ie,1) = tott_sum(1,1)
        tr_sum(7,ir,ie,2) = tott_sum(1,2)
        do k=1,np
          tr_sum(29+k,ir,ie,1) = tott_sum(k+1,1)
          tr_sum(29+k,ir,ie,2) = tott_sum(k+1,2)
        enddo
      else if(itype_het == 3) then
        tr_sum(11,ir,ie,1) = tott_sum(1,1)
        tr_sum(11,ir,ie,2) = tott_sum(1,2)
        do k=1,np
           tr_sum(29+np+k,ir,ie,1) = tott_sum(k+1,1)
           tr_sum(29+np+k,ir,ie,2) = tott_sum(k+1,2)
        enddo
      endif

      return

      end

************************************************************************
*                                                                      *
      subroutine restore_heatreg_sumover(m,iax,
     &           ik,ir,ie,fact_sum,nd,nr,ne)
*                                                                      *
************************************************************************
      
      use RESTALMOD

      implicit double precision (a-h,o-z)


C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)

      tr_sum =>  trRES_sum(irestalm_sum(m,iax):)
      
      call restore_heatreg_sumover_sub(m,
     &     ik,ir,ie,fact_sum,
     &     nd,nr,ne,tr_sum)

      return
      end


************************************************************************
*                                                                      *
      subroutine restore_heatreg_sumover_sub(m,
     &            ik,ir,ie,fact_sum,
     &            nd_sum,nr_sum,ne_sum,tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      dimension   tr_sum(nd_sum,nr_sum,0:ne_sum,2)

          call invert_heat_stdev(m,A,B,
     &                           tr_sum(ik,ir,ie,1),
     &                           tr_sum(ik,ir,ie,2),
     &                           fact_sum)

          tr_sum(ik,ir,ie,1) = A
          tr_sum(ik,ir,ie,2) = B

      return
      end

***********************************************************************
*                                                                     *
* heatrz                                                              *
*                                                                     *
***********************************************************************

***********************************************************************
*                                                                     *
      subroutine  restheatrz_sumover_put(m,iax,
     &    kl,ks,np,itype_het,idd,
     &    ik,ir,iz,ie,maxtott,nsame,tott_sum)
*                                                                     *
***********************************************************************

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

      dimension kl(7,7)
      dimension ks(7)

      real(8),pointer :: tr_sum(:)


      real(8) :: tott_sum(maxtott,2)

      tr_sum => trRES_sum(irestalm_sum(m,iax):)

      call restheatrz_sumover_put_sub(
     &    kl,ks,np,itype_het,idd,
     &    ik,ir,iz,ie,maxtott,nsame,tott_sum,
     &    itndy_sum(m,iax),itrnm_sum(m,iax),itznm_sum(m,iax),
     &    itenm_sum(m,iax),
     &    tr_sum)

      return
      end

************************************************************************
*                                                                      *
      subroutine restheatrz_sumover_put_sub(
     &    kl,ks,np,itype_het,idd,
     &    ik,ir,iz,ie,maxtott,nsame,tott_sum,
     &    nd_sum,nr_sum,nz_sum,ne_sum,tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      dimension kl(7,7)
      dimension ks(7)

      real(8) :: tr_sum(nd_sum,nr_sum,nz_sum,0:ne_sum,2)

      real(8) :: tott_sum(maxtott,2)

      if(itype_het == 1) then
        do k = 1, ks(idd)
           tr_sum(kl(idd,k),ir,iz,ie,1) = tott_sum(k,1)
           tr_sum(kl(idd,k),ir,iz,ie,2) = tott_sum(k,2)
        end do
      else if(itype_het == 2) then
        tr_sum(7,ir,iz,ie,1) = tott_sum(1,1)
        tr_sum(7,ir,iz,ie,2) = tott_sum(1,2)
        do k=1,np
          tr_sum(29+k,ir,iz,ie,1) = tott_sum(k+1,1)
          tr_sum(29+k,ir,iz,ie,2) = tott_sum(k+1,2)
        enddo
      else if(itype_het == 3) then
        tr_sum(11,ir,iz,ie,1) = tott_sum(1,1)
        tr_sum(11,ir,iz,ie,2) = tott_sum(1,2)
        do k=1,np
           tr_sum(29+np+k,ir,iz,ie,1) = tott_sum(k+1,1)
           tr_sum(29+np+k,ir,iz,ie,2) = tott_sum(k+1,2)
        enddo
      endif

      return

      end

************************************************************************
*                                                                      *
      subroutine restore_heatrz_sumover(m,iax,
     &           ik,ir,iz,ie,fact_sum,nd,nr,nz,ne)
*                                                                      *
************************************************************************
      
      use RESTALMOD

      implicit double precision (a-h,o-z)


C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)

      tr_sum =>  trRES_sum(irestalm_sum(m,iax):)
      
      call restore_heatrz_sumover_sub(m,
     &     ik,ir,iz,ie,fact_sum,
     &     nd,nr,nz,ne,tr_sum)

      return
      end


************************************************************************
*                                                                      *
      subroutine restore_heatrz_sumover_sub(m,
     &            ik,ir,iz,ie,fact_sum,
     &            nd_sum,nr_sum,nz_sum,ne_sum,tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      dimension   tr_sum(nd_sum,nr_sum,nz_sum,0:ne_sum,2)

          call invert_heat_stdev(m,A,B,
     &                           tr_sum(ik,ir,iz,ie,1),
     &                           tr_sum(ik,ir,iz,ie,2),
     &                           fact_sum)

          tr_sum(ik,ir,iz,ie,1) = A
          tr_sum(ik,ir,iz,ie,2) = B

      return
      end

***********************************************************************
*                                                                     *
* heatxyz                                                              *
*                                                                     *
***********************************************************************

***********************************************************************
*                                                                     *
      subroutine  restheatxyz_sumover_put(m,iax,
     &    kl,ks,np,itype_het,idd,
     &    ik,ix,iy,iz,ie,maxtott,nsame,tott_sum)
*                                                                     *
***********************************************************************

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

      dimension kl(7,7)
      dimension ks(7)

      real(8),pointer :: tr_sum(:)


      real(8) :: tott_sum(maxtott,2)

      tr_sum => trRES_sum(irestalm_sum(m,iax):)

      call restheatxyz_sumover_put_sub(
     &    kl,ks,np,itype_het,idd,
     &    ik,ix,iy,iz,ie,maxtott,nsame,tott_sum,
     &    itndy_sum(m,iax),itxnm_sum(m,iax),itynm_sum(m,iax),
     &    itznm_sum(m,iax),itenm_sum(m,iax),
     &    tr_sum)

      return
      end

************************************************************************
*                                                                      *
      subroutine restheatxyz_sumover_put_sub(
     &    kl,ks,np,itype_het,idd,
     &    ik,ix,iy,iz,ie,maxtott,nsame,tott_sum,
     &    nd_sum,nx_sum,ny_sum,nz_sum,ne_sum,tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      dimension kl(7,7)
      dimension ks(7)

      real(8) :: tr_sum(nd_sum,nx_sum,ny_sum,nz_sum,0:ne_sum,2)

      real(8) :: tott_sum(maxtott,2)

      if(itype_het == 1) then
        do k = 1, ks(idd)
           tr_sum(kl(idd,k),ix,iy,iz,ie,1) = tott_sum(k,1)
           tr_sum(kl(idd,k),ix,iy,iz,ie,2) = tott_sum(k,2)
        end do
      else if(itype_het == 2) then
        tr_sum(7,ix,iy,iz,ie,1) = tott_sum(1,1)
        tr_sum(7,ix,iy,iz,ie,2) = tott_sum(1,2)
        do k=1,np
          tr_sum(29+k,ix,iy,iz,ie,1) = tott_sum(k+1,1)
          tr_sum(29+k,ix,iy,iz,ie,2) = tott_sum(k+1,2)
        enddo
      else if(itype_het == 3) then
        tr_sum(11,ix,iy,iz,ie,1) = tott_sum(1,1)
        tr_sum(11,ix,iy,iz,ie,2) = tott_sum(1,2)
        do k=1,np
           tr_sum(29+np+k,ix,iy,iz,ie,1) = tott_sum(k+1,1)
           tr_sum(29+np+k,ix,iy,iz,ie,2) = tott_sum(k+1,2)
        enddo
      endif

      return

      end


************************************************************************
*                                                                      *
      subroutine restore_heatxyz_sumover(m,iax,
     &           ik,ix,iy,iz,ie,fact_sum,nd,nx,ny,nz,ne)
*                                                                      *
************************************************************************
      
      use RESTALMOD

      implicit double precision (a-h,o-z)


C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)

      tr_sum =>  trRES_sum(irestalm_sum(m,iax):)
      
      call restore_heatxyz_sumover_sub(m,
     &     ik,ix,iy,iz,ie,fact_sum,
     &     nd,nx,ny,nz,ne,tr_sum)

      return
      end


************************************************************************
*                                                                      *
      subroutine restore_heatxyz_sumover_sub(m,
     &            ik,ix,iy,iz,ie,fact_sum,
     &            nd_sum,nx_sum,ny_sum,nz_sum,ne_sum,tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      dimension   tr_sum(nd_sum,nx_sum,ny_sum,nz_sum,0:ne_sum,2)

          call invert_heat_stdev(m,A,B,
     &                           tr_sum(ik,ix,iy,iz,ie,1),
     &                           tr_sum(ik,ix,iy,iz,ie,2),
     &                           fact_sum)

          tr_sum(ik,ix,iy,iz,ie,1) = A
          tr_sum(ik,ix,iy,iz,ie,2) = B

      return
      end
