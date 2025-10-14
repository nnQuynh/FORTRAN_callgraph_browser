************************************************************************
*                                                                      *
      subroutine check_twwg(m,iax,jsn,jsi,dsin,idsi,ill,ilf,ierr)
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

      dimension lschn(35), ischn(35)
      character schan(35)*8

      data icsu / 35 /

      data ( schan(i), i = 1, 35 ) /
     &    'mesh    ','part    ','e-type  ','unit    ','axis    ',
     &    'file    ','title   ','angel   ','2d-type ','factor  ',
     &    'x-txt   ','y-txt   ','z-txt   ','gshow   ','rshow   ',
     &    'iechrl  ','material','volmat  ','epsout  ','ctmin(1)',
     &    'ctmax(1)','ctmin(2)','ctmax(2)','ctmin(3)','ctmax(3)',
     &    'resol   ','width   ','multipli','t-type  ','trcl    ',
     &    '*trcl   ','gslat   ','resfile ','ginfo   ','a-type  '/

      data ( lschn(i), i = 1, 35 ) /
     &     4,         4,         6,         4,         4,
     &     4,         5,         5,         7,         6,
     &     5,         5,         5,         5,         5,
     &     6,         8,         6,         6,         8,
     &     8,         8,         8,         8,         8,
     &     5,         5,         8,         6,         4,
     &     5,         5,         7,         5,         6/

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
      data      tname /'[t-wwg]'/

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

      integer,allocatable :: mtetreg(:)

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

            else if( chlw(ic:ic+2) .eq. 'xyz' ) then

               imesh = 3

            else if( chlw(ic:ic+2) .eq. 'tet' ) then
cFURUTA20240110
               imesh = 4

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

            else if( imesh .eq. 3 ) then

               call txymesh(jsn,jsi,dsin,idsi,ill,ilf,
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                      ixtp,inx,xmin,xmax,xdel,istxg,
     &                      iytp,iny,ymin,ymax,ydel,istyg,
     &                      iztp,inz,zmin,zmax,zdel,istzg)

                  if( ierr .ne. 0 ) return
                  if( jpn  .eq. 3 ) goto 800
                  if( ixtp .lt. 0 ) goto 996
                  if( iytp .lt. 0 ) goto 996
                  if( iztp .lt. 0 ) goto 996

            else if( imesh .eq. 4 ) then

               call moddas_allocate_int(MAX_NUM_MTRG,mtetreg)

               call ttetmesh0(jsn,jsi,dsin,idsi,ill,ilf,
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                      MAX_NUM_MTRG,mtetreg)

                  if( jpn  .eq. 3 ) goto 800

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

               goto 939

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

               if( iunt .ne. 1 ) goto 938

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        axis
*-----------------------------------------------------------------------

         else if( ipm .eq. 5 ) then

  600       inaxi = inaxi + 1

            if( inaxi .gt. 6 ) goto 991

            if( chlw(ic:ic+2) .eq. 'wwg' ) then

               iaxis(inaxi) = 7
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+2) .eq. 'eng' ) then

               iaxis(inaxi) = 1
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               iaxis(inaxi) = 2
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+2) .eq. 'tet' ) then

               iaxis(inaxi) = 11
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+1) .eq. 'xy' .or.
     &               chlw(ic:ic+1) .eq. 'yx' ) then

               iaxis(inaxi) = 3
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic+1) .eq. 'yz' .or.
     &               chlw(ic:ic+1) .eq. 'zy' ) then

               iaxis(inaxi) = 4
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic+1) .eq. 'zx' .or.
     &               chlw(ic:ic+1) .eq. 'xz' ) then

               iaxis(inaxi) = 5
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic) .eq. 't' ) then

               iaxis(inaxi) = 6
               ic =jnumc(chlw,ic+2,icl)

            else if( chlw(ic:ic) .eq. 'x' ) then

               iaxis(inaxi) = 8
               ic =jnumc(chlw,ic+2,icl)

            else if( chlw(ic:ic) .eq. 'y' ) then

               iaxis(inaxi) = 9
               ic =jnumc(chlw,ic+2,icl)

            else if( chlw(ic:ic) .eq. 'z' ) then

               iaxis(inaxi) = 10
               ic =jnumc(chlw,ic+2,icl)

            else

               goto 937

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

      else if ( itmsh(m) .eq. 4 ) then

        call check_tet(m,MAX_NUM_ITREG,mtetreg,
     &      iec,cepn,lepn,ierr)
        call moddas_deallocate_int(mtetreg)

      end if

      call check_unit(m,iunt,iec,cepn,lepn,ierr)

      call check_axis(m,iaxis,inaxi,iec,cepn,lepn,ierr)

      call check_type('e',iec,cepn,lepn,
     &                itety(m),rtema(m),rtemi(m),itenm(m),
     &                ietp, emax, emin, ine)

      call check_type('t',iec,cepn,lepn,
     &                ittty(m),rttma(m),rttmi(m),ittnm(m),
     &                ittp, tmax, tmin, int)

      call check_factor(m,rfact,iec,cepn,lepn,ierr)

      call check_part(m,iptyp,ipnkf,inpat,iec,cepn,lepn,ierr)

*-----------------------------------------------------------------------

      call judge_tall_check(iec,cepn,lepn,dsin(jsn),idsi(jsn),ierr)

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  933    m_err = 'multiplier set should be one in tally '//tname
         ErrCha = ''
         ErrID = 'L:1082/R:check_twwg/F:restwwg.f'
         goto 999

  934    m_err = 'emin should be zero in tally '//tname
         ErrCha = ''
         ErrID = 'L:1087/R:check_twwg/F:restwwg.f'
         goto 999

  935    m_err = 'ne or nt shoud be one in tally '//tname
         ErrCha = ''
         ErrID = 'L:1092/R:check_twwg/F:restwwg.f'
         goto 999

  936    m_err = 'gshow cannot be used in tally '//tname
         ErrCha = ''
         ErrID = 'L:1097/R:check_twwg/F:restwwg.f'
         goto 999

  937    m_err = 'axis should be reg,eng,t,xy,yz,zx in tally '//tname
         ErrCha = ''
         ErrID = 'L:1102/R:check_twwg/F:restwwg.f'
         goto 999

  938    m_err = 'unit should be 1 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1107/R:check_twwg/F:restwwg.f'
         goto 999

  939    m_err = 'angle mesh cannot be used in tally '//tname
         ErrCha = ''
         ErrID = 'L:1112/R:check_twwg/F:restwwg.f'
         goto 999

  940    m_err = 'mesh=reg should be chosen in tally '//tname
         ErrCha = ''
         ErrID = 'L:1117/R:check_twwg/F:restwwg.f'
         goto 999

*-----------------------------------------------------------------------

  950    write(dkam,'(i9)') mdas
         m_err = 'Total tally storage number exceeds mdas ='//dkam
         ErrCha = ''
         ErrID = 'L:1125/R:check_twwg/F:restwwg.f'
         goto 999

  968    m_err = 'mset number is inconsistent.'//tname
         ErrCha = ''
         ErrID = 'L:1130/R:check_twwg/F:restwwg.f'
         goto 999

  969    m_err = 'number of mset should be the same.'//tname
         ErrCha = ''
         ErrID = 'L:1135/R:check_twwg/F:restwwg.f'
         goto 999

  970    m_err = 'def of multiplier should be less than 7 '//tname
         ErrCha = ''
         ErrID = 'L:1140/R:check_twwg/F:restwwg.f'
         goto 999

  971    m_err = 'number of multiplier is negative'//tname
         ErrCha = ''
         ErrID = 'L:1145/R:check_twwg/F:restwwg.f'
         goto 999

  978    m_err = 'counter should be from -9999 to 9999'
         ErrCha = ''
         ErrID = 'L:1150/R:check_twwg/F:restwwg.f'
         goto 999

  979    m_err = 'Description of material parameter is wrong in tally '
     &            //tname
         ErrCha = ''
         ErrID = 'L:1156/R:check_twwg/F:restwwg.f'
         goto 999

  983    m_err = '2D-type should be 1-7 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1161/R:check_twwg/F:restwwg.f'
         goto 999

  985    m_err = 'Unit is Lethargy but energy mesh points are negative'
         ErrCha = ''
         ErrID = 'L:1166/R:check_twwg/F:restwwg.f'
         goto 999

  986    write(dkam,'(i5)') itlmax
         m_err = 'Total tally number exceeds itlmax = '//dkam
         ErrCha = ''
         ErrID = 'L:1172/R:check_twwg/F:restwwg.f'
         goto 999

  987    m_err = 'Unknown parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1177/R:check_twwg/F:restwwg.f'
         goto 999

  988    m_err = 'Unit should be 1,2,3,4,11,12,13,14 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1182/R:check_twwg/F:restwwg.f'
         goto 999

  989    m_err = 'Double definition of the parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1187/R:check_twwg/F:restwwg.f'
         goto 999

  990    m_err = 'Too many file name in tally '//tname
         ErrCha = ''
         ErrID = 'L:1192/R:check_twwg/F:restwwg.f'
         goto 999

  991    m_err = 'Too many axis in tally '//tname
         ErrCha = ''
         ErrID = 'L:1197/R:check_twwg/F:restwwg.f'
         goto 999

  992    m_err = 'Unknown axis name in tally '//tname
         ErrCha = ''
         ErrID = 'L:1202/R:check_twwg/F:restwwg.f'
         goto 999

  993    m_err = 'In this line, [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:1207/R:check_twwg/F:restwwg.f'
         goto 999

  994    m_err = 'Name of particle is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1212/R:check_twwg/F:restwwg.f'
         goto 999

  995    m_err = 'Number of particles is larger than 6 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1217/R:check_twwg/F:restwwg.f'
         goto 999

  996    m_err = 'Description of mesh is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1222/R:check_twwg/F:restwwg.f'
         goto 999

  997    m_err = 'Description of parameter is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1227/R:check_twwg/F:restwwg.f'
         goto 999

  998    m_err = 'Unknown mesh parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1232/R:check_twwg/F:restwwg.f'
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
      subroutine read_twwg(m,iax,ierr)
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
        include 'err.inc'

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

! T.Sato 2024/03/18 weighted history counter ID
      common /tall92/ ichnum(itlmax),chbias(itlmax),pedest(itlmax)
     &               ,ictnum(itlmax),ctbias(itlmax),ictidx(itlmax)

        dimension     idas(mdas*2)
        equivalence ( das, idas )

*-----------------------------------------------------------------------

        character dsin(0:9,2)*200
        dimension idsi(0:9,2)
        dimension ill(0:9,2), ilf(0:9,2)
        dimension jsn(2), jsi(2)
        dimension ierrs(2)

*-----------------------------------------------------------------------

      character m_err*200

*-----------------------------------------------------------------------

        ierr = 0
        noe  = 1
        if(ichnum(m).ne.0.or.ictnum(m).ne.0) then ! T.Sato 2024/03/24
         ErrCha = 'Warning: biased factor in [t-wwg] determined by '//
     &   '(history) counter is initialized even in the restart '//
     &   'calculation'
         ErrID = 'L:1345/R:read_twwg/F:restwwg.f'
         call ErrWrite(ErrID, ErrCha)
        endif

        if(itaxs(m,iax).eq.7)then !FURUTA20240111
         goto 998
        endif

cKN 2018/01/08
        if ( any( itaxs(m,iax) .eq. (/ 3, 4, 5 /) )
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
        call check_twwg(m,iax,jsn(ioe),jsi(ioe),
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

          call read_wwgreg(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itrgn(m),itrgm(m),itenm(m),itmst(m),
     &                    ittnm(m),
     &                    idas_itreg(itreg(m)),
     &                    das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)),
     &                    itrnv(m),idas(itriv(m)),das(itrrv(m)))

        else if ( itmsh(m) .eq. 3 ) then

          call read_wwgxyz(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itmtn(m),ismte(itmtt(m)),
     &                    itxnm(m),itynm(m),itznm(m),itenm(m),itmst(m),
     &                    ittnm(m),
     &                    das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                    das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

        else if( itmsh(m) .eq. 4 ) then

          call read_wwgtet(m,iax,ioe,jsn(ioe),jsi(ioe),
     &                    dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                    itpan(m),itrgn(m),itenm(m),itmst(m),ittnm(m),
     &                    das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

        end if

  800   continue
*-----------------------------------------------------------------------
*   restore tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call restore_wwgreg(m,
     &                    itpan(m),itrgn(m),itrgm(m),itenm(m),itmst(m),
     &                    ittnm(m),
     &                    idas_itreg(itreg(m)),
     &                    das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)),
     &                    itrnv(m),idas(itriv(m)),das(itrrv(m)) )

        else if ( itmsh(m) .eq. 3 ) then

          call restore_wwgxyz(m,
     &                    itpan(m),itmtn(m),ismte(itmtt(m)),
     &                    itxnm(m),itynm(m),itznm(m),itenm(m),itmst(m),
     &                    ittnm(m),
     &                    das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                    das_itzrg(itzrg(m)),
     &                    das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &                    trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 4 ) then

          call restore_wwgtet(m,
     &                    itpan(m),itrgn(m),itrgm(m),
     &                    itenm(m),itmst(m),ittnm(m),
     &                    idas_itreg(itreg(m)),
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

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  998    m_err = '1st axis of T-WWG shold not be wwg for restart'
         ErrCha = ''
         ErrID = 'L:1483/R:read_twwg/F:restwwg.f'
         goto 999

*-----------------------------------------------------------------------
  999   continue
         ierr  = 1

         call ErrWrite(ErrID, ErrCha)
         write(*,*) 'Error: ' // m_err

*-----------------------------------------------------------------------

      return

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_wwgreg(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                        np,nr,mr,ne,nm,nt,kr,eb,tb,tr,
     &                        nvl,ivl,rvl)
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

           nm_0 = 0
           nr_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petrm.inc'

          nmstepi = 1
          do 190 im = 1, nm, nmstepi
          do 190 ir = 1, nr, nrstepi
          do 190 it = 1, nt, ntstepi
          do 190 ip = 1, np, npstepi

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#e-lowere-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nestepi = 1
            do ie = 1, ne, nestepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &      ir+irloop-1,im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),imloop=1,nmstepi)

            end do

! sumover
            ie = 1
            call pwwgreg_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,ntstepi,nrstepi,nmstepi,
     &      ip,ie,it,ir,im,nsame,tott_sum)


  190     continue

*-----------------------------------------------------------------------
*        reg axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 2 ) then

           nm_0 = 0
           nr_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petrm.inc'

           nmstepi = 1
           do 290 im = 1, nm, nmstepi
           do 290 ie = 1, ne, nestepi
           do 290 it = 1, nt, ntstepi
           do 290 ip = 1, np, npstepi

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#numregvolume',13,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             nrstepi = 1
             do ir = 1, nr, nrstepi


               read(jsi,'(i5,1x,i7,1pe13.4,1000(1pe13.4,0pf8.4))')
     &              idmm0, idmm1, dmm2,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &      ir+irloop-1,im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),imloop=1,nmstepi)

             end do

! sumover
            ir = 1
            call pwwgreg_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,ntstepi,nrstepi,nmstepi,
     &      ip,ie,it,ir,im,nsame,tott_sum)

  290     continue

*-----------------------------------------------------------------------
*        time axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 6 ) then

           nm_0 = 0
           nr_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petrm.inc'

           nmstepi = 1
           do 390 im = 1, nm, nmstepi
           do 390 ir = 1, nr, nrstepi
           do 390 ie = 1, ne, nestepi
           do 390 ip = 1, np, npstepi

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#t-lowert-upper',15,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             ntstepi = 1
             do it = 1, nt, ntstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &      ir+irloop-1,im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),imloop=1,nmstepi)

             end do

! sumover
            it = 1
            call pwwgreg_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,ntstepi,nrstepi,nmstepi,
     &      ip,ie,it,ir,im,nsame,tott_sum)

  390     continue

*-----------------------------------------------------------------------
         end if

  900 continue

      end subroutine

************************************************************************
*                                                                      *
      subroutine read_wwgtet(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                        np,nr,ne,nm,nt,eb,tb,tr)
*                                                                      *
*   m: the tally number, index of ital.                                *
*                                                                      *
*     Last Modified by T.Furuta on 2024/01/10                          *
*                                                                      *
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

*-----------------------------------------------------------------------

        dimension   eb(ne+1)
        dimension   ew(ne)
        dimension   tb(nt+1)
        dimension   tr(np,ne,nt,nr,nm,2)

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

cnais 2023/01/31
           nm_0 = 0
           nr_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petrm.inc'

          nmstepi = 1
          do 190 im = 1, nm, nmstepi
          do 190 ir = 1, nr, nrstepi
          do 190 it = 1, nt, ntstepi
          do 190 ip = 1, np, npstepi

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#e-lowere-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 900
            if ( jpn   .eq. 3 ) goto 900

            nestepi = 1
            do ie = 1, ne, nestepi

cfrtati 2021/10/05
c              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
c     &             dmm0,dmm1,
c     &             ((tr(ip,ie,it,ir,im,k),k=1,2),ip=1,np)

              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0,dmm1,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &      ir+irloop-1,im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),imloop=1,nmstepi)

            end do

! sumover
            ie = 1
            call pwwgreg_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,ntstepi,nrstepi,nmstepi,
     &      ip,ie,it,ir,im,nsame,tott_sum)

  190     continue

*-----------------------------------------------------------------------
*        tet axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 11 ) then

cnais 2023/01/31
           nm_0 = 0
           nr_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petrm.inc'

           nmstepi = 1
           do 290 im = 1, nm, nmstepi
           do 290 ie = 1, ne, nestepi
           do 290 it = 1, nt, ntstepi
           do 290 ip = 1, np, npstepi

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#numtetravolume',15,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             nrstepi = 1
             do ir = 1, nr, nrstepi

cfrtati 2021/10/05
c               read(jsi,'(i8,1x,i8,1pe13.4,1000(1pe13.4,0pf8.4))')
c     &              idmm0, idmm1, dmm2,
c     &              ((tr(ip,ie,it,ir,im,k),k=1,2),ip=1,np)

               read(jsi,'(i8,1x,i8,1pe13.4,1000(1pe13.4,0pf8.4))')
     &              idmm0, idmm1, dmm2,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &      ir+irloop-1,im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),imloop=1,nmstepi)

             end do

! sumover
            ir = 1
            call pwwgreg_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,ntstepi,nrstepi,nmstepi,
     &      ip,ie,it,ir,im,nsame,tott_sum)

  290     continue

*-----------------------------------------------------------------------
*        time axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 6 ) then

           nm_0 = 0
           nr_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petrm.inc'

           nmstepi = 1
           do 390 im = 1, nm, nmstepi
           do 390 ir = 1, nr, nrstepi
           do 390 ie = 1, ne, nestepi
           do 390 ip = 1, np, npstepi

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#t-lowert-upper',15,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             ntstepi = 1
             do it = 1, nt, ntstepi

cfrtati 2021/10/05
c              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
c     &             dmm0, dmm1,
c     &             ((tr(ip,ie,it,ir,im,k),k=1,2),ip=1,np)

              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             dmm0, dmm1,
     &      ((((((tr(ip+iploop-1,ie+ieloop-1,it+itloop-1,
     &      ir+irloop-1,im+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),imloop=1,nmstepi)

             end do

! sumover
            it = 1
            call pwwgreg_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,ntstepi,nrstepi,nmstepi,
     &      ip,ie,it,ir,im,nsame,tott_sum)

  390     continue

*-----------------------------------------------------------------------
         end if

  900 continue

      end subroutine



************************************************************************
*                                                                      *
      subroutine read_wwgxyz(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
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

           nm_0 = 0
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petxyzm.inc'

          nmstepi = 1
          do 190 im = 1, nm, nmstepi
          do 190 ix = 1, nx, nxstepi
          do 190 iy = 1, ny, nystepi
          do 190 iz = 1, nz, nzstepi
          do 190 it = 1, nt, ntstepi
          do 190 ip = 1, np, npstepi

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                       '#e-lowere-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 999
            if ( jpn   .eq. 3 ) goto 999

            nestepi = 1
            do ie = 1, ne, nestepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
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
            call pwwgxyz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,ntstepi,nxstepi,nystepi,nzstepi,nmstepi,
     &      ip,ie,it,ix,iy,iz,im,nsame,tott_sum)

  190     continue

*-----------------------------------------------------------------------
*        x axis
*-----------------------------------------------------------------------
cKN 2023/12/26
c       else if( itaxs(m,iax) .eq. 3 ) then
        else if( itaxs(m,iax) .eq. 8 ) then
          !'#x-lowerx-upper'

           nm_0 = 0
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petxyzm.inc'

          nmstepi = 1
          do 290 im = 1, nm, nmstepi
          do 290 ie = 1, ne, nestepi
          do 290 iy = 1, ny, nystepi
          do 290 iz = 1, nz, nzstepi
          do 290 it = 1, nt, ntstepi
          do 290 ip = 1, np, npstepi

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                       '#x-lowerx-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 999
            if ( jpn   .eq. 3 ) goto 999

            nxstepi = 1
            do ix = 1, nx, nxstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
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
            call pwwgxyz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,ntstepi,nxstepi,nystepi,nzstepi,nmstepi,
     &      ip,ie,it,ix,iy,iz,im,nsame,tott_sum)

  290     continue

*-----------------------------------------------------------------------
*        y axis
*-----------------------------------------------------------------------
cKN 2023/12/26
c       else if( itaxs(m,iax) .eq. 4 ) then
        else if( itaxs(m,iax) .eq. 9 ) then
          !'#y-lowery-upper'

           nm_0 = 0
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petxyzm.inc'

          nmstepi = 1
          do 390 im = 1, nm, nmstepi
          do 390 ie = 1, ne, nestepi
          do 390 ix = 1, nx, nxstepi
          do 390 iz = 1, nz, nzstepi
          do 390 it = 1, nt, ntstepi
          do 390 ip = 1, np, npstepi

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                       '#y-lowery-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 999
            if ( jpn   .eq. 3 ) goto 999

            nystepi = 1
            do iy = 1, ny, nystepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
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
            call pwwgxyz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,ntstepi,nxstepi,nystepi,nzstepi,nmstepi,
     &      ip,ie,it,ix,iy,iz,im,nsame,tott_sum)

  390     continue

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------
cKN 2023/12/26
c       else if( itaxs(m,iax) .eq. 5 ) then
        else if( itaxs(m,iax) .eq. 10 ) then
          !'#z-lowerz-upper'

           nm_0 = 0
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petxyzm.inc'

          nmstepi = 1
          do 490 im = 1, nm, nmstepi
          do 490 ie = 1, ne, nestepi
          do 490 ix = 1, nx, nxstepi
          do 490 iy = 1, ny, nystepi
          do 490 it = 1, nt, ntstepi
          do 490 ip = 1, np, npstepi

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                       '#z-lowerz-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 999
            if ( jpn   .eq. 3 ) goto 999

            nzstepi = 1
            do iz = 1, nz, nzstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
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
            call pwwgxyz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,ntstepi,nxstepi,nystepi,nzstepi,nmstepi,
     &      ip,ie,it,ix,iy,iz,im,nsame,tott_sum)

  490     continue

*-----------------------------------------------------------------------
*        time axis
*-----------------------------------------------------------------------
cKN 2023/12/26
c       else if( itaxs(m,iax) .eq. 11 ) then
        else if( itaxs(m,iax) .eq. 6 ) then
          !'#t-lowert-upper'

           nm_0 = 0
           nz_0 = 0
           ny_0 = 0
           nx_0 = 0
           nt_0 = 0
           ne_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_petxyzm.inc'

          nmstepi = 1
          do 590 im = 1, nm, nmstepi
          do 590 ie = 1, ne, nestepi
          do 590 ix = 1, nx, nxstepi
          do 590 iy = 1, ny, nystepi
          do 590 iz = 1, nz, nzstepi
          do 590 ip = 1, np, npstepi

            call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                       '#t-lowert-upper',15,ierr)

            if ( ierr  .ne. 0 ) goto 999
            if ( jpn   .eq. 3 ) goto 999

            ntstepi = 1
            do it = 1, nt, ntstepi


              read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
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
            call pwwgxyz_sumover_getput(jsi,m,iax,
     &      npstepi,nestepi,ntstepi,nxstepi,nystepi,nzstepi,nmstepi,
     &      ip,ie,it,ix,iy,iz,im,nsame,tott_sum)

  590     continue

*-----------------------------------------------------------------------
*        xy axis (matrix)
*-----------------------------------------------------------------------
cKN 2023/12/26
c       else if( itaxs(m,iax) .eq. 7 ) then
        else if( itaxs(m,iax) .eq. 3 ) then

          call sel_dc2(m,dc2,ldc2)
          if( ittwo(m) .eq. 4 ) then
            dc2 = '#xyfluxr.err'
            ldc2 = 12
          else if( ittwo(m) .eq. 5 ) then
            dc2 = 'y/x'
            ldc2 = 3
          end if

*-----------------------------------------------------------------------

          do 690 im = 1, nm
          do 690 iz = 1, nz
          do 690 ip = 1, npg
          do 690 ie = 1, neg
          do 690 it = 1, ntg

            call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)
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
cKN 2023/12/26
c       else if( itaxs(m,iax) .eq. 8 ) then
        else if( itaxs(m,iax) .eq. 4 ) then

          call sel_dc2(m,dc2,ldc2)
          if( ittwo(m) .eq. 4 ) then
            dc2 = '#yzfluxr.err'
            ldc2 = 12
          else if( ittwo(m) .eq. 5 ) then
            dc2 = 'y/z'
            ldc2 = 3
          end if

*-----------------------------------------------------------------------

          do 790 im = 1, nm
          do 790 ix = 1, nx
          do 790 ip = 1, npg
          do 790 ie = 1, neg
          do 790 it = 1, ntg

            call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)
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
cKN 2023/12/26
c       else if( itaxs(m,iax) .eq. 9 ) then
        else if( itaxs(m,iax) .eq. 5 ) then

          call sel_dc2(m,dc2,ldc2)
          if( ittwo(m) .eq. 4 ) then
            dc2 = '#xzfluxr.err'
            ldc2 = 12
          else if( ittwo(m) .eq. 5 ) then
            dc2 = 'x/z'
            ldc2 = 3
          end if

*-----------------------------------------------------------------------

          do 890 im = 1, nm
          do 890 iy = 1, ny
          do 890 ip = 1, npg
          do 890 ie = 1, neg
          do 890 it = 1, ntg

            call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)
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
      subroutine restore_wwgreg(m,
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

        common /tall21/ rtfac(itlmax)

        common /cparm/  maxbch,maxcas
        common /restart/ resc2(itlmax), resc3(itlmax)

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
*           itunt(m) = 1   /cm^2/source
*-----------------------------------------------------------------------

            ew(1:ne) = 1.d+0
            ew_sum = 1.d0

*-----------------------------------------------------------------------
*           itunt(m) = 1
*-----------------------------------------------------------------------

          tw(1:nt) = 1.d+0
          tw_sum = 1.0d0

*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

        call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)

        vl_sum = sum(vl(:))

*-----------------------------------------------------------------------

        do 100 im = 1, nm
        do 100 ir = 1, nr
        do 100 it = 1, nt
        do 100 ie = 1, ne
        do 100 ip = 1, np

          call invert_stdev(m,A,B,
     &                      tr(ip,ie,it,ir,im,1),
     &                      tr(ip,ie,it,ir,im,2),
     &                      vl(ir)*ew(ie)*tw(it)/rtfac(m))

          tr(ip,ie,it,ir,im,1) = A
          tr(ip,ie,it,ir,im,2) = B

! sumover
         fact_in = rtfac(m)
         call pwwgreg_sumover_stdev(1,m,ip,ie,it,ir,im,
     &        fact_in,ew(ie),tw(it),vl(ir),ew_sum,tw_sum,vl_sum)

  100 continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_wwgtet(m,
     &                           np,nr,mr,ne,nm,nt,kr,eb,tb,tr)
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

        common /cparm/  maxbch,maxcas
        common /restart/ resc2(itlmax), resc3(itlmax)

*-----------------------------------------------------------------------

        dimension   kr(mr)
        dimension   eb(ne+1)
        dimension   ew(ne)
        dimension   vl(nr)
        dimension   lr(nr)
        dimension   tb(nt+1)
        dimension   tw(nt)
        dimension   tr(np,ne,nt,nr,nm,2)

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*           itunt(m) = 1   /cm^2/source
*-----------------------------------------------------------------------

            ew(1:ne) = 1.d+0
            ew_sum = 1.0d0

*-----------------------------------------------------------------------
*           itunt(m) = 1
*-----------------------------------------------------------------------

          tw(1:nt) = 1.d+0
          tw_sum = 1.0d0

*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

        call ttetvl(mr,kr,nr,vl,lr)
        vl_sum = sum(vl(:))

*-----------------------------------------------------------------------

        do 100 im = 1, nm
        do 100 ir = 1, nr
        do 100 it = 1, nt
        do 100 ie = 1, ne
        do 100 ip = 1, np

          call invert_stdev(m,A,B,
     &                      tr(ip,ie,it,ir,im,1),
     &                      tr(ip,ie,it,ir,im,2),
     &                      vl(ir)*ew(ie)*tw(it)/rtfac(m))

          tr(ip,ie,it,ir,im,1) = A
          tr(ip,ie,it,ir,im,2) = B

! sumover
         fact_in = rtfac(m)
         call pwwgreg_sumover_stdev(1,m,ip,ie,it,ir,im,
     &        fact_in,ew(ie),tw(it),vl(ir),ew_sum,tw_sum,vl_sum)

  100 continue

      end subroutine


************************************************************************
cKN 2018/01/08
*                                                                      *
      subroutine restore_wwgxyz(m,
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
        common /tall21/ rtfac(itlmax)
        common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
        common /tall34/ itvm(itlmax)
        common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

        common /cparm/  maxbch,maxcas
        common /restart/ resc2(itlmax), resc3(itlmax)

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
        real(8),allocatable :: vl_x(:,:),vl_y(:,:),vl_z(:,:)

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
! sumover
       allocate (vl_x(ny,nz),vl_y(nx,nz),vl_z(nx,ny))
       vl_x(:,:) = 0.0d0
       vl_y(:,:) = 0.0d0
       vl_z(:,:) = 0.0d0
       do iz = 1, nz
         do iy = 1, ny
           do ix = 1, nx
             vl_x(iy,iz) = vl_x(iy,iz) + vl(ix,iy,iz)
             vl_y(ix,iz) = vl_y(ix,iz) + vl(ix,iy,iz)
             vl_z(ix,iy) = vl_z(ix,iy) + vl(ix,iy,iz)
           enddo
         enddo
       enddo

          ew(1:ne) = 1.d+0
          ew_sum = 1.0d0

          tw(1:nt) = 1.d+0
          tw_sum = 1.0d0

*-----------------------------------------------------------------------

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
     &                      vl(ix,iy,iz)*ew(ie)*tw(it)/rtfac(m) )

          tr(ip,ie,it,icf(ix,iy,iz),im,1) = A
          tr(ip,ie,it,icf(ix,iy,iz),im,2) = B

! sumover
              fact_in = rtfac(m)
              call pwwgxyz_sumover_stdev(1,m,ip,ie,it,ix,iy,iz,im,
     &             fact_in,ew(ie),tw(it),vl(ix,iy,iz),
     &             ew_sum,tw_sum,vl_x(iy,iz),vl_y(ix,iz),vl_z(ix,iy))

  100 continue

      deallocate (vl_x,vl_y,vl_z)

      end subroutine


