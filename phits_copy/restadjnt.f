************************************************************************
*                                                                      *
      subroutine check_tadjnt(m,iax,jsn,jsi,dsin,idsi,ill,ilf,ierr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************
      use moddas
      use moddas_tally

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'param01.inc'
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

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

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

      dimension lschn(38), ischn(38)
      character schan(38)*8

      data icsu / 38 /

      data ( schan(i), i = 1, 38 ) /
     &    'mesh    ','part    ','e-type  ','unit    ','axis    ',
     &    'file    ','title   ','angel   ','2d-type ','factor  ',
     &    'x-txt   ','y-txt   ','z-txt   ','gshow   ','rshow   ',
     &    'iechrl  ','material','volmat  ','epsout  ','ctmin(1)',
     &    'ctmax(1)','ctmin(2)','ctmax(2)','ctmin(3)','ctmax(3)',
     &    'resol   ','width   ','multipli','t-type  ','trcl    ',
     &    '*trcl   ','gslat   ','resfile ','ginfo   ','a-type  ',
     &    'e-smin  ','e-smax  ','m-source'/

      data ( lschn(i), i = 1, 38 ) /
     &     4,         4,         6,         4,         4,
     &     4,         5,         5,         7,         6,
     &     5,         5,         5,         5,         5,
     &     6,         8,         6,         6,         8,
     &     8,         8,         8,         8,         8,
     &     5,         5,         8,         6,         4,
     &     5,         5,         7,         5,         6,
     &     6,         6,         8/

*-----------------------------------------------------------------------

cfrtati 2021/10/05 6 -> mxpart
      dimension iptyp(mxpart), ipnkf(mxpart)
      dimension imtyp(mxpart,mxpart), imnkf(mxpart,mxpart)
      dimension jstyp(mxpart), jnkf0(mxpart)
      dimension icount(9)

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

      character tname*11
      data      tname /'[t-adjoint]'/

      logical deqn4

      character rglnrf*200, rglnech*200

*-----------------------------------------------------------------------

      dimension jptyp(6), jpnkf(6)
      dimension kptyp(6,6), kpnkf(6,6)
      dimension ktln(6), ktli(6), ktls(6), kmst(6), knpat(6)
      dimension dkmax(6)
      dimension imst(6), kimst(6,6)
      dimension vtrs(13)
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
            inpat = 0
            inaxi = 0
            infil = 0
            idtyp = 3
            jmul  = 0

            langel = 0
            lxtxt  = 0
            lytxt  = 0
            lztxt  = 0
            lgshow = 0
            lrshow = 0
            iechrl = 72
            matvol = 9
            ieps   = 0
            ittp   = 0
            igkst  = 0
            idtt   = 0
            ktrs   = 0
            igslt  = 1
            lrfile = 0

            iatp   = 0
            jatp   = 0

            amin = 0.0
            amax = 0.0
            adel = 0.0

            infog = 0

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
     &          ipm .ne.  2 .and. ipm .ne. 28 .and.
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
*        time mesh
*-----------------------------------------------------------------------

         else if( ipm .eq. 29 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

               ittp = nint( cvvv )

               call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     't',ittp,int,tmin,tmax,tdel,isttg)


               goto 150

*-----------------------------------------------------------------------
*        angle mesh
*-----------------------------------------------------------------------

         else if( ipm .eq. 35 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

               iatp = nint( cvvv )
               jatp = abs( iatp )

               call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     'a',jatp,ina,amin,amax,adel,istag)

               if( ierr .ne. 0 ) return

               goto 150

*-----------------------------------------------------------------------
*        transform
*-----------------------------------------------------------------------

         else if( ipm .eq. 30 .or. ipm .eq. 31 ) then

                  if( ipm .eq. 31 ) ktrs = 1

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
     &           ( iunt .gt. 5 .and. iunt .lt. 11 ) .or.
     &             iunt .gt. 14 ) goto 988

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

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               iaxis(inaxi) = 2
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+2) .eq. 'rad' ) then

               iaxis(inaxi) = 12
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+2) .eq. 'deg' ) then

               iaxis(inaxi) = 13
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

            else if( chlw(ic:ic) .eq. 't' ) then

               iaxis(inaxi) = 11
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

         else if( ipm .eq. 33 ) then

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

         else if( ipm .eq. 32 ) then

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
*        volmat
*-----------------------------------------------------------------------

         else if( ipm .eq. 18 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               matvol = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        epsout
*-----------------------------------------------------------------------

         else if( ipm .eq. 19 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ieps = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        resolution
*-----------------------------------------------------------------------

         else if( ipm .eq. 26 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ireso = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        width
*-----------------------------------------------------------------------

         else if( ipm .eq. 27 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               width = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        counter
*-----------------------------------------------------------------------

         else if( ipm .ge. 20 .and. ipm .le. 25 ) then

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

         else if( ipm .eq. 28 ) then

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
*        infog
*-----------------------------------------------------------------------

         else if( ipm .eq. 34 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               infog = nint( cvvv )

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

      call check_type('e',iec,cepn,lepn,
     &                itety(m),rtema(m),rtemi(m),itenm(m),
     &                ietp, emax, emin, ine)

      call check_type('t',iec,cepn,lepn,
     &                ittty(m),rttma(m),rttmi(m),ittnm(m),
     &                ittp, tmax, tmin, int)

      call check_type('a',iec,cepn,lepn,
     &                itaty(m),rtama(m),rtami(m),itanm(m),
     &                iatp, amax, amin, ina)

      call check_2dtype(m,idtyp,iec,cepn,lepn,ierr)

      call check_factor(m,rfact,iec,cepn,lepn,ierr)


*-----------------------------------------------------------------------

      call judge_tall_check(iec,cepn,lepn,dsin(jsn),idsi(jsn),ierr)

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  950    write(dkam,'(i9)') mdas
         m_err = 'Total tally storage number exceeds mdas ='//dkam
         ErrCha = ''
         ErrID = 'L:1118/R:check_tadjnt/F:restadjnt.f'
         goto 999

  968    m_err = 'mset number is inconsistent.'//tname
         ErrCha = ''
         ErrID = 'L:1123/R:check_tadjnt/F:restadjnt.f'
         goto 999

  969    m_err = 'number of mset should be the same.'//tname
         ErrCha = ''
         ErrID = 'L:1128/R:check_tadjnt/F:restadjnt.f'
         goto 999

  970    m_err = 'def of multiplier should be less than 7 '//tname
         ErrCha = ''
         ErrID = 'L:1133/R:check_tadjnt/F:restadjnt.f'
         goto 999

  971    m_err = 'number of multiplier is negative'//tname
         ErrCha = ''
         ErrID = 'L:1138/R:check_tadjnt/F:restadjnt.f'
         goto 999

  978    m_err = 'counter should be from -9999 to 9999'
         ErrCha = ''
         ErrID = 'L:1143/R:check_tadjnt/F:restadjnt.f'
         goto 999

  979    m_err = 'Description of material parameter is wrong in tally '
     &            //tname
         ErrCha = ''
         ErrID = 'L:1149/R:check_tadjnt/F:restadjnt.f'
         goto 999

  983    m_err = '2D-type should be 1-7 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1154/R:check_tadjnt/F:restadjnt.f'
         goto 999

  985    m_err = 'Unit is Lethargy but energy mesh points are negative'
         ErrCha = ''
         ErrID = 'L:1159/R:check_tadjnt/F:restadjnt.f'
         goto 999

  986    write(dkam,'(i5)') itlmax
         m_err = 'Total tally number exceeds itlmax = '//dkam
         ErrCha = ''
         ErrID = 'L:1165/R:check_tadjnt/F:restadjnt.f'
         goto 999

  987    m_err = 'Unknown parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1170/R:check_tadjnt/F:restadjnt.f'
         goto 999

  988    m_err = 'Unit should be 1,2,3,4,11,12,13,14 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1175/R:check_tadjnt/F:restadjnt.f'
         goto 999

  989    m_err = 'Double definition of the parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1180/R:check_tadjnt/F:restadjnt.f'
         goto 999

  990    m_err = 'Too many file name in tally '//tname
         ErrCha = ''
         ErrID = 'L:1185/R:check_tadjnt/F:restadjnt.f'
         goto 999

  991    m_err = 'Too many axis in tally '//tname
         ErrCha = ''
         ErrID = 'L:1190/R:check_tadjnt/F:restadjnt.f'
         goto 999

  992    m_err = 'Unknown axis name in tally '//tname
         ErrCha = ''
         ErrID = 'L:1195/R:check_tadjnt/F:restadjnt.f'
         goto 999

  993    m_err = 'In this line, [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:1200/R:check_tadjnt/F:restadjnt.f'
         goto 999

  994    m_err = 'Name of particle is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1205/R:check_tadjnt/F:restadjnt.f'
         goto 999

  995    m_err = 'Number of particles is larger than 6 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1210/R:check_tadjnt/F:restadjnt.f'
         goto 999

  996    m_err = 'Description of mesh is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1215/R:check_tadjnt/F:restadjnt.f'
         goto 999

  997    m_err = 'Description of parameter is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1220/R:check_tadjnt/F:restadjnt.f'
         goto 999

  998    m_err = 'Unknown mesh parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1225/R:check_tadjnt/F:restadjnt.f'
         goto 999

*-----------------------------------------------------------------------

  999 continue

         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1

         call ErrWrite(ErrID, ErrCha)
         write(*,*) 'Error: ' // m_err

*-----------------------------------------------------------------------

      return

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_tadjnt(m,iax,ierr)
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
        common /tall08/ rtrx0(itlmax), rtry0(itlmax)
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
        common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                  rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
        common /tall21/ rtfac(itlmax)
        common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
        common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)
        common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

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
        call check_tadjnt(m,iax,jsn(ioe),jsi(ioe),
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
        idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1

*-----------------------------------------------------------------------
*   read tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call read_tadjreg(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itrgn(m),itrgm(m),itenm(m),itmst(m),
     &                    ittnm(m),
     &                    idas_itreg(itreg(m)),
     &                    das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)),
     &                    itrnv(m),idas(itriv(m)),das(itrrv(m)))


        else if ( itmsh(m) .eq. 2 ) then

          call read_tadjrz(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itrnm(m),itznm(m),itenm(m),itmst(m),
     &                    ittnm(m),itanm(m),
     &                    das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &                    das_itarg(itarg(m)),das_iterg(iterg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))


        else if ( itmsh(m) .eq. 3 ) then

          call read_tadjxyz(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itmtn(m),ismte(itmtt(m)),
     &                    itxnm(m),itynm(m),itznm(m),itenm(m),itmst(m),
     &                    ittnm(m),
     &                    das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                    das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

        end if

  800   continue
*-----------------------------------------------------------------------
*   restore tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call restore_tadjreg(m,
     &                    itpan(m),itrgn(m),itrgm(m),itenm(m),itmst(m),
     &                    ittnm(m),
     &                    idas_itreg(itreg(m)),
     &                    das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)),
     &                    itrnv(m),idas(itriv(m)),das(itrrv(m)) )

        else if ( itmsh(m) .eq. 2 ) then

          call restore_tadjrz(m,
     &                    itpan(m),itrnm(m),itznm(m),itenm(m),itmst(m),
     &                    ittnm(m),itanm(m),
     &                    das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &                    das_itarg(itarg(m)),das_iterg(iterg(m)),
     &                    das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 3 ) then

          call restore_tadjxyz(m,
     &                    itpan(m),itmtn(m),ismte(itmtt(m)),
     &                    itxnm(m),itynm(m),itznm(m),itenm(m),itmst(m),
     &                    ittnm(m),
     &                    das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                    das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

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
      subroutine read_tadjreg(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                        np,nr,mr,ne,nm,nt,kr,eb,tb,tr,
     &                        nvl,ivl,rvl)
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
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
        common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
        character ctfln*100

        common /tall21/ rtfac(itlmax)

        common /cparm/  maxbch,maxcas
        common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
        common /restart/ resc2(itlmax), resc3(itlmax)

        common /fact01/ facmax(itlmax)
        common /stat / istdev, irestart, ireschk

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
*        energy axis
*-----------------------------------------------------------------------
        if( itaxs(m,iax) .eq. 1 ) then

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nr_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petrm.inc'

          do 190 im = 1, nm, nmstepi
          do 190 ir = 1, nr, nrstepi
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


              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &             dmm0,dmm1,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &      ir+irloop-1,im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),imloop=1,nmstepi)

            end do

! sumover
            ie = 1
            call ptracreg_sumover_getput(jsi,m,iax,                ! ie = 1
     &        npstepi,nestepi,ntstepi,nrstepi,nmstepi,
     &        ip,ie,it,ir,im,nsame,tott_sum)

  190     continue

*-----------------------------------------------------------------------
*        reg axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 2 ) then

           facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nr_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petrm.inc'

           do 290 im = 1, nm, nmstepi
           do 290 ie = 1, ne, nestepi
           do 290 it = 1, nt, ntstepi
           do 290 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#numregvolume',13,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             nrstepi = 1
             do ir = 1, nr, nrstepi


               read(jsi,'(i5,1x,i7,1pe13.4,7(1pe13.4,0pf8.4))')
     &              idmm0, idmm1, dmm2,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &      ir+irloop-1,im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),imloop=1,nmstepi)

             end do

! sumover
            ir = 1
            call ptracreg_sumover_getput(jsi,m,iax,                ! ir = 1
     &        npstepi,nestepi,ntstepi,nrstepi,nmstepi,
     &        ip,ie,it,ir,im,nsame,tott_sum)

  290     continue

*-----------------------------------------------------------------------
*        time axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 11 ) then

           facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nr_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petrm.inc'

           do 390 im = 1, nm, nmstepi
           do 390 ir = 1, nr, nrstepi
           do 390 ie = 1, ne, nestepi
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


              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &      ir+irloop-1,im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),imloop=1,nmstepi)

             end do

! sumover
            it = 1
            call ptracreg_sumover_getput(jsi,m,iax,                ! it = 1
     &        npstepi,nestepi,ntstepi,nrstepi,nmstepi,
     &        ip,ie,it,ir,im,nsame,tott_sum)

  390     continue

*-----------------------------------------------------------------------
         end if

  900 continue

      end subroutine



************************************************************************
*                                                                      *
      subroutine read_tadjrz(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                       np,nr,nz,ne,nm,nt,na,rm,zm,ab,eb,tb,tr)
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
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
        common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
        character ctfln*100

        common /tall21/ rtfac(itlmax)

        common /cparm/  maxbch,maxcas
        common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
        common /restart/ resc2(itlmax), resc3(itlmax)

        common /fact01/ facmax(itlmax)
        common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------

        dimension rm(nr+1)
        dimension zm(nz+1)
        dimension eb(ne+1)
        dimension ab(na+1)
        dimension ew(ne)
        dimension tb(nt+1)
        dimension tw(nt)

        dimension aw(na)
        dimension tr(np,ne,nt,nr*nz,na,nm,2)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        character chin*200, chlw*200, chcm*200
        character dc2*15

! sumover dummy
        dimension tott_sum(1,2)
        data nsame/1/

*-----------------------------------------------------------------------

      icf(ir,iz) = ir + ( iz - 1 ) * nr

*-----------------------------------------------------------------------
*        angle axis
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 12 .or. itaxs(m,iax) .eq. 13 ) then
          !'#a-lowera-upper'

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           na_0 = 0
           nz_0 = 0
           nr_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petrzam.inc'

          do 690 im = 1, nm, nmstepi
          do 690 ir = 1, nr, nrstepi
          do 690 iz = 1, nz, nzstepi
          do 690 ie = 1, ne, nestepi
          do 690 it = 1, nt, ntstepi
          do 690 ip = 1, np, npstepi

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


              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &       icf(ir+irloop-1,iz+izloop-1),
     &      ia+ialoop-1,im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi),ialoop=1,nastepi),
     &      imloop=1,nmstepi)

            end do

! sumover
            ia = 1
            call ptracrz_sumover_getput(jsi,m,iax,
     &        npstepi,nestepi,ntstepi,nrstepi,nzstepi,nastepi,nmstepi,
     &        ip,ie,it,ir,iz,ia,im,nsame,tott_sum)

  690     continue

*-----------------------------------------------------------------------
*        energy axis
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 1 ) then
          !'#e-lowere-upper'

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           na_0 = 0
           nz_0 = 0
           nr_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petrzam.inc'

          nmstepi = 1
          do 190 im = 1, nm, nmstepi
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


              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &       icf(ir+irloop-1,iz+izloop-1),
     &      ia+ialoop-1,im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi),ialoop=1,nastepi),
     &      imloop=1,nmstepi)

            end do

! sumover
            ie = 1
            call ptracrz_sumover_getput(jsi,m,iax,
     &        npstepi,nestepi,ntstepi,nrstepi,nzstepi,nastepi,nmstepi,
     &        ip,ie,it,ir,iz,ia,im,nsame,tott_sum)

  190     continue

*-----------------------------------------------------------------------
*        r axis
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 6 ) then
          !'#r-lowerr-upper'

          facmx = 1.d0  ! kitamura23/03/31

           nm_0 = 0
           na_0 = 0
           nz_0 = 0
           nr_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petrzam.inc'

          do 290 im = 1, nm, nmstepi
          do 290 ie = 1, ne, nestepi
          do 290 iz = 1, nz, nzstepi
          do 290 ia = 1, na, nastepi
          do 290 it = 1, nt, ntstepi
          do 290 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                       '#r-lowerr-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nrstepi = 1
            do ir = 1, nr, nrstepi


              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &       icf(ir+irloop-1,iz+izloop-1),
     &      ia+ialoop-1,im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi),ialoop=1,nastepi),
     &      imloop=1,nmstepi)

            end do

! sumover
            ir = 1
            call ptracrz_sumover_getput(jsi,m,iax,
     &        npstepi,nestepi,ntstepi,nrstepi,nzstepi,nastepi,nmstepi,
     &        ip,ie,it,ir,iz,ia,im,nsame,tott_sum)

  290     continue

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 5 ) then
          !'#z-lowerz-upper'

          facmx = 1.d0  ! kitamura23/03/31

           nm_0 = 0
           na_0 = 0
           nz_0 = 0
           nr_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petrzam.inc'

          do 390 im = 1, nm, nmstepi
          do 390 ie = 1, ne, nestepi
          do 390 ir = 1, nr, nrstepi
          do 390 ia = 1, na, nastepi
          do 390 it = 1, nt, ntstepi
          do 390 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                       '#z-lowerz-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nzstepi = 1
            do iz = 1, nz, nzstepi


              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &       icf(ir+irloop-1,iz+izloop-1),
     &      ia+ialoop-1,im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi),ialoop=1,nastepi),
     &      imloop=1,nmstepi)

            end do

! sumover
            iz = 1
            call ptracrz_sumover_getput(jsi,m,iax,
     &        npstepi,nestepi,ntstepi,nrstepi,nzstepi,nastepi,nmstepi,
     &        ip,ie,it,ir,iz,ia,im,nsame,tott_sum)

  390     continue

*-----------------------------------------------------------------------
*        time axis
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 11 ) then
          !'#t-lowert-upper'

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           na_0 = 0
           nz_0 = 0
           nr_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petrzam.inc'

          do 490 im = 1, nm, nmstepi
          do 490 ir = 1, nr, nrstepi
          do 490 iz = 1, nz, nzstepi
          do 490 ia = 1, na, nastepi
          do 490 ie = 1, ne, nestepi
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


              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &       icf(ir+irloop-1,iz+izloop-1),
     &      ia+ialoop-1,im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi),ialoop=1,nastepi),
     &      imloop=1,nmstepi)

            end do

! sumover
            it = 1
            call ptracrz_sumover_getput(jsi,m,iax,
     &        npstepi,nestepi,ntstepi,nrstepi,nzstepi,nastepi,nmstepi,
     &        ip,ie,it,ir,iz,ia,im,nsame,tott_sum)

  490     continue

*-----------------------------------------------------------------------
*        rz axis (matrix)
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 10 ) then

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

          do 590 im = 1, nm
          do 590 ip = 1, np
          do 590 ie = 1, ne
          do 590 ia = 1, na
          do 590 it = 1, nt

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
     &          ((tr(ip,ie,it,icf(ir,iz),ia,im,ioe), iz = 1,nz ),
     &                                            ir = nr,1,-1 )

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 4 ) then

              do iz = 1, nz
              do ir = 1, nr

               read(jsi,'(1p3e11.3,0pf8.4)')
     &           dmm0, dmm1,
     &           tr(ip,ie,it,icf(ir,iz),ia,im,1),
     &           tr(ip,ie,it,icf(ir,iz),ia,im,2)

              end do
              end do

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 5 ) then

              do ir = nr, 1, -1

                read(jsi,'(1p1000e11.3)')
     &            dmm0,
     &            ( tr(ip,ie,it,icf(ir,iz),ia,im,ioe), iz = 1, nz )

              end do

*-----------------------------------------------------------------------
            end if

  590     continue


*-----------------------------------------------------------------------
        end if

  800 continue

*-----------------------------------------------------------------------

  900 continue

      end subroutine



************************************************************************
*                                                                      *
      subroutine read_tadjxyz(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                       np,nl,lt,nx,ny,nz,ne,nm,nt,xm,ym,zm,
     &                       eb,tb,tr)
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
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
        common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
        character ctfln*100

        common /tall21/ rtfac(itlmax)
        common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
        common /tall34/ itvm(itlmax)
        common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

        common /cparm/  maxbch,maxcas
        common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
        common /restart/ resc2(itlmax), resc3(itlmax)

        common /fact01/ facmax(itlmax)
        common /stat / istdev, irestart, ireschk

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

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        character chin*200, chlw*200, chcm*200
        character dc2*15

! sumover dummy
        dimension :: tott_sum(1,2)
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
            ntg = nt

         else

            npg = 1
            neg = 1
            ntg = 1

         end if

*-----------------------------------------------------------------------
*        energy axis
*-----------------------------------------------------------------------
        if( itaxs(m,iax) .eq. 1 ) then
          !'#e-lowere-upper'

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petxyzm.inc'

          do 190 im = 1, nm, nmstepi
          do 190 ix = 1, nx, nxstepi
          do 190 iy = 1, ny, nystepi
          do 190 iz = 1, nz, nzstepi
          do 190 it = 1, nt, ntstepi
          do 190 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                       '#e-lowere-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 999
            if ( jpn   .eq. 3 ) goto 999

            nestepi = 1
            do ie = 1, ne, nestepi


              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &           dmm0, dmm1,
     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &       icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &      im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi),
     &      imloop=1,nmstepi)

            end do

! sumover
            ie = 1
            call ptracxyz_sumover_getput(jsi,m,iax,
     &        npstepi,nestepi,ntstepi,nxstepi,nystepi,nzstepi,nmstepi,
     &        ip,ie,it,ix,iy,iz,im,nsame,tott_sum)

  190     continue

*-----------------------------------------------------------------------
*        x axis
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 3 ) then
          !'#x-lowerx-upper'

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petxyzm.inc'

          do 290 im = 1, nm, nmstepi
          do 290 ie = 1, ne, nestepi
          do 290 iy = 1, ny, nystepi
          do 290 iz = 1, nz, nzstepi
          do 290 it = 1, nt, ntstepi
          do 290 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                       '#x-lowerx-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 999
            if ( jpn   .eq. 3 ) goto 999

            nxstepi = 1
            do ix = 1, nx, nxstepi


              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &           dmm0, dmm1,
     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &       icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &      im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi),
     &      imloop=1,nmstepi)

            end do

! sumover
            ix = 1
            call ptracxyz_sumover_getput(jsi,m,iax,
     &        npstepi,nestepi,ntstepi,nxstepi,nystepi,nzstepi,nmstepi,
     &        ip,ie,it,ix,iy,iz,im,nsame,tott_sum)

  290     continue

*-----------------------------------------------------------------------
*        y axis
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 4 ) then
          !'#y-lowery-upper'

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petxyzm.inc'

          do 390 im = 1, nm, nmstepi
          do 390 ie = 1, ne, nestepi
          do 390 ix = 1, nx, nxstepi
          do 390 iz = 1, nz, nzstepi
          do 390 it = 1, nt, ntstepi
          do 390 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                       '#y-lowery-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 999
            if ( jpn   .eq. 3 ) goto 999

            nystepi = 1
            do iy = 1, ny, nystepi


              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &           dmm0, dmm1,
     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &       icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &      im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi),
     &      imloop=1,nmstepi)

            end do

! sumover
            iy = 1
            call ptracxyz_sumover_getput(jsi,m,iax,
     &        npstepi,nestepi,ntstepi,nxstepi,nystepi,nzstepi,nmstepi,
     &        ip,ie,it,ix,iy,iz,im,nsame,tott_sum)

  390     continue

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 5 ) then
          !'#z-lowerz-upper'

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petxyzm.inc'

          do 490 im = 1, nm, nmstepi
          do 490 ie = 1, ne, nestepi
          do 490 ix = 1, nx, nxstepi
          do 490 iy = 1, ny, nystepi
          do 490 it = 1, nt, ntstepi
          do 490 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                       '#z-lowerz-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 999
            if ( jpn   .eq. 3 ) goto 999

            nzstepi = 1
            do iz = 1, nz, nzstepi


              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &           dmm0, dmm1,
     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &       icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &      im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi),
     &      imloop=1,nmstepi)

            end do

! sumover
            iz = 1
            call ptracxyz_sumover_getput(jsi,m,iax,
     &        npstepi,nestepi,ntstepi,nxstepi,nystepi,nzstepi,nmstepi,
     &        ip,ie,it,ix,iy,iz,im,nsame,tott_sum)

  490     continue

*-----------------------------------------------------------------------
*        time axis
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 11 ) then
          !'#t-lowert-upper'

          facmx = 1.d0  ! kitamura23/03/31
           nm_0 = 0
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petxyzm.inc'

          do 590 im = 1, nm, nmstepi
          do 590 ie = 1, ne, nestepi
          do 590 ix = 1, nx, nxstepi
          do 590 iy = 1, ny, nystepi
          do 590 iz = 1, nz, nzstepi
          do 590 ip = 1, np, npstepi

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                       '#t-lowert-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 999
            if ( jpn   .eq. 3 ) goto 999

            ntstepi = 1
            do it = 1, nt, ntstepi


              read(jsi,'(1p2e13.4,7(1pe13.4,0pf8.4))')
     &           dmm0, dmm1,
     &      ((((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &       icf(ix+ixloop-1,iy+iyloop-1,iz+izloop-1),
     &      im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi),
     &      imloop=1,nmstepi)

            end do

! sumover
            it = 1
            call ptracxyz_sumover_getput(jsi,m,iax,
     &        npstepi,nestepi,ntstepi,nxstepi,nystepi,nzstepi,nmstepi,
     &        ip,ie,it,ix,iy,iz,im,nsame,tott_sum)

  590     continue

*-----------------------------------------------------------------------
*        xy axis (matrix)
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

*-----------------------------------------------------------------------

          facmx = 1.d0  ! kitamura23/03/31

          do 690 im = 1, nm
          do 690 iz = 1, nz
          do 690 ip = 1, npg
          do 690 ie = 1, neg
          do 690 it = 1, ntg

            call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                       dc2,ldc2,ierr)

            if ( ierr  .ne. 0 ) goto 999
            if ( jpn   .eq. 3 ) goto 999

*-----------------------------------------------------------------------
            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

              read(jsi,'(1p10e11.3)')
     &        ( ( tr(ip,ie,it,icf(ix,iy,iz),im,ioe),
     &            ix = 1, nx ), iy = ny, 1, -1 )

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 4 ) then

              do iy = 1, ny
              do ix = 1, nx

                read(jsi,'(1p3e13.4,0pf8.4)')
     &            dmm0, dmm1,
     &            (tr(ip,ie,it,icf(ix,iy,iz),im,k),k=1,2)

              end do
              end do

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 5 ) then

              do iy = ny, 1, -1

                read(jsi,'(1p1000e11.3)')
     &            dmm0,
     &            ( tr(ip,ie,it,icf(ix,iy,iz),im,ioe), ix = 1, nx )

              end do

*-----------------------------------------------------------------------
            end if

  690     continue

*-----------------------------------------------------------------------
*        yz axis (matrix)
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 8 ) then

          call sel_dc2(m,dc2,ldc2)
          if( ittwo(m) .eq. 4 ) then
            dc2 = '#yzfluxr.err'
            ldc2 = 12
          else if( ittwo(m) .eq. 5 ) then
            dc2 = 'y/z'
            ldc2 = 3
          end if

*-----------------------------------------------------------------------

          facmx = 1.d0  ! kitamura23/03/31

          do 790 im = 1, nm
          do 790 ix = 1, nx
          do 790 ip = 1, npg
          do 790 ie = 1, neg
          do 790 it = 1, ntg

            call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                       dc2,ldc2,ierr)

            if ( ierr  .ne. 0 ) goto 999
            if ( jpn   .eq. 3 ) goto 999

*-----------------------------------------------------------------------
            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

              read(jsi,'(1p10e11.3)')
     &        ( ( tr(ip,ie,it,icf(ix,iy,iz),im,ioe),
     &               iz = 1, nz ), iy = ny, 1, -1 )

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 4 ) then

              do iz = 1, nz
              do iy = 1, ny

                read(jsi,'(1p3e13.4,0pf8.4)')
     &            dmm0, dmm1,
     &            (tr(ip,ie,it,icf(ix,iy,iz),im,k),k=1,2)

              end do
              end do

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 5 ) then

              do iy = ny, 1, -1

                read(jsi,'(1p1000e11.3)')
     &          dmm0,
     &          ( tr(ip,ie,it,icf(ix,iy,iz),im,ioe), iz = 1, nz )

              end do

*-----------------------------------------------------------------------
            end if

  790     continue

*-----------------------------------------------------------------------
*        xz axis (matrix)
*-----------------------------------------------------------------------
        else if( itaxs(m,iax) .eq. 9 ) then

          call sel_dc2(m,dc2,ldc2)
          if( ittwo(m) .eq. 4 ) then
            dc2 = '#xzfluxr.err'
            ldc2 = 12
          else if( ittwo(m) .eq. 5 ) then
            dc2 = 'x/z'
            ldc2 = 3
          end if

*-----------------------------------------------------------------------

          facmx = 1.d0  ! kitamura23/03/31

          do 890 im = 1, nm
          do 890 iy = 1, ny
          do 890 ip = 1, npg
          do 890 ie = 1, neg
          do 890 it = 1, ntg

            call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                       dc2,ldc2,ierr)

            if ( ierr  .ne. 0 ) goto 999
            if ( jpn   .eq. 3 ) goto 999

*-----------------------------------------------------------------------
            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

              read(jsi,'(1p10e11.3)')
     &        ( ( tr(ip,ie,it,icf(ix,iy,iz),im,ioe),
     &               iz = 1, nz ), ix = nx, 1, -1 )

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 4 ) then

              do iz = 1, nz
              do ix = 1, nx

                read(jsi,'(1p3e11.3,0pf8.4)')
     &             dmm0, dmm1,
     &             tr(ip,ie,it,icf(ix,iy,iz),im,1),
     &             tr(ip,ie,it,icf(ix,iy,iz),im,2)

              end do
              end do

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 5 ) then

              do ix = nx, 1, -1

                read(jsi,'(1p1000e11.3)')
     &          dmm0,
     &          ( tr(ip,ie,it,icf(ix,iy,iz),im,ioe), iz = 1, nz )

              end do

*-----------------------------------------------------------------------
            end if

  890     continue

*-----------------------------------------------------------------------
        end if

  999 continue

      end subroutine



************************************************************************
*                                                                      *
      subroutine restore_tadjreg(m,
     &                           np,nr,mr,ne,nm,nt,kr,eb,tb,tr,
     &                           nvl,ivl,rvl)
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
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)

        common /tall21/ rtfac(itlmax)

        common /cparm/  maxbch,maxcas
        common /restart/ resc2(itlmax), resc3(itlmax)

        common /fact01/ facmax(itlmax)  ! kitamura23/03/31

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

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*           itunt(m) = 1,4 11,14: /cm^2/source
*                                 /cm^2/sec/source
*                    = 2   12   : /cm^2/MeV/source
*                                 /cm^2/MeV/nsec/source
*                    = 3   13   : /cm^2/Lethargy/source
*                                 /cm^2/Lethargy/nsec/source
*-----------------------------------------------------------------------

         if( itunt(m) .eq.  1 .or. itunt(m) .eq.  4 .or.
     &       itunt(m) .eq. 11 .or. itunt(m) .eq. 14 ) then

            ew(1:ne) = 1.d+0
            ew_sum = 1.0d0

         else if( itunt(m) .eq. 2 .or. itunt(m) .eq. 12 ) then

            ew(1:ne) = eb(2:ne+1) - eb(1:ne)
            ew_sum  = eb(ne+1) -eb(1)

         else if( itunt(m) .eq. 3 .or. itunt(m) .eq. 13 ) then

            ew(1:ne) = log( eb(2:ne+1) / eb(1:ne) )
            ew_sum = log(eb(ne+1) /eb(1))

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

*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

        call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)
        vl_sum = sum(vl(1:nr))

*-----------------------------------------------------------------------
*        ( unit = 4 or 14  ; vol = 1.0 )
*-----------------------------------------------------------------------

        if( itunt(m) .eq. 4 .or. itunt(m) .eq. 14 ) then

          vl(1:nr) = 1.0d0
          vl_sum = 1.0d0

        end if

*-----------------------------------------------------------------------

        if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

        do 100 im = 1, nm
        do 100 ir = 1, nr
        do 100 it = 1, nt
        do 100 ie = 1, ne
        do 100 ip = 1, np

          call invert_stdev(m,A,B,
     &                      tr(ip,ie,it,ir,im,1),
     &                      tr(ip,ie,it,ir,im,2),
     &       vl(ir)*ew(ie)*tw(it)/abs(rtfac(m)/facmax(m)))

          tr(ip,ie,it,ir,im,1) = A
          tr(ip,ie,it,ir,im,2) = B

! sumover
          fact_in = abs(rtfac(m)/facmax(m))
          call ptracreg_sumover_stdev(1,m,ip,ie,it,ir,im,
     &                           fact_in,ew(ie),tw(it),vl(ir),
     &                                  ew_sum,tw_sum,vl_sum)

  100 continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_tadjrz(m,
     &                       np,nr,nz,ne,nm,nt,na,rm,zm,ab,eb,tb,tr)
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
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
        common /tall21/ rtfac(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

        common /cparm/  maxbch,maxcas
        common /restart/ resc2(itlmax), resc3(itlmax)

        common /fact01/ facmax(itlmax)  ! kitamura23/03/31

*-----------------------------------------------------------------------

        dimension rm(nr+1)
        dimension zm(nz+1)
        dimension eb(ne+1)
        dimension ab(na+1)
        dimension ew(ne)
        dimension tb(nt+1)
        dimension tw(nt)

        dimension aw(na)
        dimension tr(np,ne,nt,nr*nz,na,nm,2)

! sumover
      real(8),allocatable :: vl_r(:),vl_z(:)

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

        pi = 2.d0 * asin(1.d0)

*-----------------------------------------------------------------------
*           itunt(m) = 1,4 11,14: /cm^2/source
*                                 /cm^2/sec/source
*                    = 2   12   : /cm^2/MeV/source
*                                 /cm^2/MeV/nsec/source
*                    = 3   13   : /cm^2/Lethargy/source
*                                 /cm^2/Lethargy/nsec/source
*-----------------------------------------------------------------------
! sumovver
      allocate (vl_r(nz),vl_z(nr))
      vl_r(:) = 0.0d0
      vl_z(:) = 0.0d0
      do iz=1,nz
         do ir=1,nr
           vl_r(iz) = vl_r(iz) + vl(ir,iz)
           vl_z(ir) = vl_z(ir) + vl(ir,iz)
         enddo
      enddo

        if( itunt(m) .eq.  1 .or. itunt(m) .eq.  4 .or.
     &      itunt(m) .eq. 11 .or. itunt(m) .eq. 14 ) then

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

               do i = 1, na

                  if( itaty(m) .gt. 0 ) then

                     aw(i) = ( ab(i+1) - ab(i) ) / 2.0 / pi
                     if(i == 1) then
                       aw_sum = ( ab(na+1) - ab(1) ) / 2.0 / pi
                     endif

                  else

                     aw(i) = ( ab(i+1) - ab(i) ) / 360.0
                     if(i == 1) then
                       aw_sum = ( ab(na+1) - ab(1) ) / 360.0
                     endif

                  end if

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

        if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

        do 100 im = 1, nm
        do 100 ia = 1, na
        do 100 iz = 1, nz
        do 100 ir = 1, nr
        do 100 it = 1, nt
        do 100 ie = 1, ne
        do 100 ip = 1, np

          call invert_stdev(m,A,B,
     &                      tr(ip,ie,it,icf(ir,iz),ia,im,1),
     &                      tr(ip,ie,it,icf(ir,iz),ia,im,2),
     &        vl(ir,iz)*aw(ia)*ew(ie)*tw(it)/abs(rtfac(m)/facmax(m)))

          tr(ip,ie,it,icf(ir,iz),ia,im,1) = A
          tr(ip,ie,it,icf(ir,iz),ia,im,2) = B

! sumover
          fact_in = abs(rtfac(m)/facmax(m))
          call ptracrz_sumover_stdev(1,m,ip,ie,it,ir,iz,ia,im,
     &                   fact_in,ew(ie),tw(it),vl(ir,iz),aw(ia),
     &                   ew_sum,tw_sum,vl_r(iz),vl_z(ir),aw_sum)

  100 continue

! sumover
      deallocate (vl_r,vl_z)

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_tadjxyz(m,
     &                       np,nl,lt,nx,ny,nz,ne,nm,nt,xm,ym,zm,
     &                       eb,tb,tr)
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
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
        common /tall21/ rtfac(itlmax)
        common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
        common /tall34/ itvm(itlmax)
        common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

        common /cparm/  maxbch,maxcas
        common /restart/ resc2(itlmax), resc3(itlmax)

        common /fact01/ facmax(itlmax)  ! kitamura23/03/31

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

! sumover
      real(8), allocatable :: vl_x(:,:),vl_y(:,:),vl_z(:,:)

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
* restore tally
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*           itunt(m) = 1,4 11,14: /cm^2/source
*                                 /cm^2/sec/source
*                    = 2   12   : /cm^2/MeV/source
*                                 /cm^2/MeV/nsec/source
*                    = 3   13   : /cm^2/Lethargy/source
*                                 /cm^2/Lethargy/nsec/source
*-----------------------------------------------------------------------

! sumovver
      allocate (vl_x(ny,nz),vl_y(nx,nz),vl_z(nx,ny))
      vl_x(:,:) = 0.0d0
      vl_y(:,:) = 0.0d0
      vl_z(:,:) = 0.0d0
      do iz=1,nz
        do iy=1,ny
          do ix =1,nx
            vl_x(iy,iz) = vl_x(iy,iz) + vl(ix,iy,iz)
            vl_y(ix,iz) = vl_y(ix,iz) + vl(ix,iy,iz)
            vl_z(ix,iy) = vl_z(ix,iy) + vl(ix,iy,iz)
          enddo
        enddo
      enddo

         if( itunt(m) .eq.  1 .or. itunt(m) .eq.  4 .or.
     &       itunt(m) .eq. 11 .or. itunt(m) .eq. 14 ) then

            ew(1:ne) = 1.d+0
            ew_sum = 1.0d0

         else if( itunt(m) .eq. 2 .or. itunt(m) .eq. 12 ) then

            ew(1:ne) = eb(2:ne+1) - eb(1:ne)
            ew_sum = eb(ne+1) - eb(1)

         else if( itunt(m) .eq. 3 .or. itunt(m) .eq. 13 ) then

            ew(1:ne) = log( eb(2:ne+1) / eb(1:ne) )
            ew_sum = log( eb(ne+1) / eb(1) )

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

*-----------------------------------------------------------------------

        if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

        do 100 im = 1, nm
        do 100 iz = 1, nz
        do 100 iy = 1, ny
        do 100 ix = 1, nx
        do 100 it = 1, nt
        do 100 ie = 1, ne
        do 100 ip = 1, np

          call invert_stdev(m,A,B,
     &                      tr(ip,ie,it,icf(ix,iy,iz),im,1),
     &                      tr(ip,ie,it,icf(ix,iy,iz),im,2),
     &        vl(ix,iy,iz)*ew(ie)*tw(it)/abs(rtfac(m)/facmax(m)))

          tr(ip,ie,it,icf(ix,iy,iz),im,1) = A
          tr(ip,ie,it,icf(ix,iy,iz),im,2) = B

! sumover
          fact_in = abs(rtfac(m)/facmax(m))
          call ptracxyz_sumover_stdev(1,m,ip,ie,it,ix,iy,iz,im,
     &             fact_in,ew(ie),tw(it),vl(ix,iy,iz),
     &             ew_sum,tw_sum,vl_x(iy,iz),vl_y(ix,iz),vl_z(ix,iy))

  100 continue

! sumover
      deallocate (vl_x,vl_y,vl_z)

      end subroutine

