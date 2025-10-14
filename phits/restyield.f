
************************************************************************
*                                                                      *
      subroutine check_tyield(m,iax,jsn,jsi,dsin,idsi,ill,ilf,ierr)
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
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)

!OBINATA: for resfile
      common /tall59/ irfll(itlmax), crfln(itlmax), itrff(itlmax)
      character crfln*100
      character irfile*100

      dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
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

      dimension lschn(36), ischn(36)
      character schan(36)*8

      data icsu / 36 /

      data ( schan(i), i = 1, 36 ) /
     &    'mesh    ','special ','mother  ','nucleus ','axis    ',
     &    'file    ','title   ','angel   ','unit    ','info    ',
     &    '2d-type ','factor  ','material','x-txt   ','y-txt   ',
     &    'z-txt   ','gshow   ','rshow   ','ndata   ','iechrl  ',
     &    'volmat  ','epsout  ','ctmin(1)','ctmax(1)','ctmin(2)',
     &    'ctmax(2)','ctmin(3)','ctmax(3)','resol   ','width   ',
     &    'trcl    ','*trcl   ','part    ','gslat   ','output  ',
     &    'resfile '/

      data ( lschn(i), i = 1, 36 ) /
     &     4,         7,         6,         7,         4,
     &     4,         5,         5,         4,         4,
     &     7,         6,         8,         5,         5,
     &     5,         5,         5,         5,         6,
     &     6,         6,         8,         8,         8,
     &     8,         8,         8,         5,         5,
     &     4,         5,         4,         5,         6,
     &     7/

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

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
      data      tname /'[t-yield]'/

      logical deqn1
      logical deqn4
      logical dcom2

      dimension ikzz(maxpt,maxnt), iknn(maxpt,maxnt)

      dimension icount(9)
      dimension vtrs(13)
cfrtati 2021/10/05 6 -> mxpart
      dimension iptyp(mxpart), ipnkf(mxpart)
      dimension imtyp(mxpart,mxpart), imnkf(mxpart,mxpart)
      dimension jstyp(mxpart), jnkf0(mxpart)

      common /subtra/ isubt, ipsub(mxpart)   ! kitamura22/03/31

*-----------------------------------------------------------------------

      character element(104)*3,cnuc*3

      data element/
     & 'h  ','he ','li ','be ','b  ','c  ','n  ','o  ',
     & 'f  ','ne ','na ','mg ','al ','si ','p  ','s  ',
     & 'cl ','ar ','k  ','ca ','sc ','ti ','v  ','cr ',
     & 'mn ','fe ','co ','ni ','cu ','zn ','ga ','ge ',
     & 'as ','se ','br ','kr ','rb ','sr ','y  ','zr ',
     & 'nb ','mo ','tc ','ru ','rh ','pd ','ag ','cd ',
     & 'in ','sn ','sb ','te ','i  ','xe ','cs ','ba ',
     & 'la ','ce ','pr ','nd ','pm ','sm ','eu ','gd ',
     & 'tb ','dy ','ho ','er ','tm ','yb ','lu ','hf ',
     & 'ta ','w  ','re ','os ','ir ','pt ','au ','hg ',
     & 'tl ','pb ','bi ','po ','at ','rn ','fr ','ra ',
     & 'ac ','th ','pa ','u  ','np ','pu ','am ','cm ',
     & 'bk ','cf ','es ','fm ','md ','no ','lr ','ku '/

*-----------------------------------------------------------------------

      character rglnrf*200, rglnech*200

      character*10 cepn(100)
      dimension lepn(100)
      character*20 ctl(2)

*-----------------------------------------------------------------------

      integer,allocatable :: mtetreg(:)

*-----------------------------------------------------------------------


            ierr  = 0
            jpn   = 0

            inmat = 0
            inaxi = 0
            infil = 0
            ispec = 0
            imoth = 0
            jmoth = 1
            inucl = 0
            imate = 0
            jmate = 1
            imasi = 0
            iunt  = 1
            info  = 0
            idtyp = 3
            inpat = 0
            iout  = 0

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

            ndata  = 0

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
     &          ischn( ipm ) .gt. 1 ) goto 989

               ic = inumc(chlw,il+1,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 997

               icl = inumc(chlw,ic,i3,';') - 1

*-----------------------------------------------------------------------
*        mesh = region, r-z, xyz, or tet
*-----------------------------------------------------------------------

         if( ipm .eq. 1 ) then

            if( chlw(ic:ic+2) .eq. 'reg' ) then

               imesh = 1

            else if( chlw(ic:ic+2) .eq. 'r-z' ) then

               imesh = 2

            else if( chlw(ic:ic+2) .eq. 'xyz' ) then

               imesh = 3

            else if( chlw(ic:ic+2) .eq. 'tet' ) then
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

         else if( ipm .eq. 33 ) then

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
*        yield special
*-----------------------------------------------------------------------

         else if( ipm .eq. 2 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ispec = nint( cvvv )

               if( ispec .lt. 0 ) goto 988

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        mother
*-----------------------------------------------------------------------

         else if( ipm .eq. 3 ) then

            if( chlw(ic:ic+2) .eq. 'all' ) then

               imoth = 0

            else if( deqn4( chlw(ic:ic) ) ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               imoth = nint( cvvv )

               if( imoth .lt. 0 ) then

                  imoth = -imoth
                  jmoth = -1

               end if

               if( imoth .eq. 0 ) goto 982

                  nsmat = 1
                  call moddas_allocate_int(imoth, ismat_temporary)

                  if( mmmax .gt. mdas ) goto 950

  141                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                        if( jpn  .eq. 3 ) goto 982

                        if( iskip .ne. 0 ) goto 141

                  ic = i1

               do k = 1, imoth

                  if( ic .gt. i3 ) then

  142                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( jpn  .eq. 3 ) goto 982

                     if( iskip .ne. 0 ) goto 142

                     ic = i1

                  end if

                        ic = jnumc(chlw,ic,i3)

                        isa = 0
                        isn = 0

                        ica = 0
                        icb = 0
                        icm = 0
                        icn = 0

                     do i = ic, i3

                        if( chlw(i:i) .ge. 'a' .and.
     &                      chlw(i:i) .le. 'z' ) then

                           isa = isa + 1

                           if( isa .eq. 1 ) ica = i

                           icb = i

                        else if( deqn1( chlw(i:i) ) ) then

                           isn = isn + 1

                           if( isn .eq. 1 ) icm = i

                           icn = i

                        else if( dcom2( chlw(i:i) ) .or.
     &                           i .eq. i3 ) then

                           icd = i
                           goto 502

                        else

                           goto 982

                        end if

                     end do

                        icd = i3

  502                continue

                     if( ica .eq. 0 .or. icb .eq. 0 ) goto 982
                     if( icb-ica .lt. 0 .or. icb-ica .gt. 1 ) goto 982

                     cnuc = chlw(ica:icb)//'  '

                     do j = 1, 104

                        if( cnuc(1:3) .eq. element(j)(1:3) ) then

                           icha = j

                           goto 452

                        end if

                     end do

                           goto 982

  452                continue

                     if( icha .gt. 104 )  goto 982

                  if( icm .eq. 0 .or. icn .eq. 0 ) then

                     ismat_temporary(nsmat-1+k) = icha * 1000

                  else

                     if( isn .gt. 3 ) goto 982

                     read(chlw(icm:icn),'(i5)') masi

                     if( masi .lt. icha ) goto 982
                     if( masi-icha .gt. maxnt ) goto 982

                     ismat_temporary(nsmat-1+k) = icha * 1000 + masi

                  end if

                     ic = icd + 1

               end do
               call moddas_deallocate_int(ismat_temporary)

            else

               goto 982

            end if

*-----------------------------------------------------------------------
*        nucleus
*-----------------------------------------------------------------------

         else if( ipm .eq. 4 ) then

            if( chlw(ic:ic+2) .eq. 'all' ) then

               inucl = 0
               imasi = 0

            else if( deqn1( chlw(ic:ic) ) ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               inucl = nint( cvvv )

               if( inucl .le. 0 ) goto 981

                  nnucl = 1
                  call moddas_allocate_int(inucl, isnuc_temporary)

                  if( mmmax .gt. mdas ) goto 950

  143                 call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                        if( jpn  .eq. 3 ) goto 981

                        if( iskip .ne. 0 ) goto 143

                  ic = i1

               do k = 1, inucl

                  if( ic .gt. i3 ) then

  144                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( jpn  .eq. 3 ) goto 981

                     if( iskip .ne. 0 ) goto 144

                     ic = i1

                  end if

                        ic = jnumc(chlw,ic,i3)

                        isa = 0
                        isn = 0

                        ica = 0
                        icb = 0
                        icm = 0
                        icn = 0

                     do i = ic, i3

                        if( chlw(i:i) .ge. 'a' .and.
     &                      chlw(i:i) .le. 'z' ) then

                           isa = isa + 1

                           if( isa .eq. 1 ) ica = i

                           icb = i

                        else if( deqn1( chlw(i:i) ) ) then

                           isn = isn + 1

                           if( isn .eq. 1 ) icm = i

                           icn = i

                        else if( dcom2( chlw(i:i) ) .or.
     &                           i .eq. i3 ) then

                           icd = i
                           goto 500

                        else

                           goto 981

                        end if

                     end do

                        icd = i3

  500                continue

                     if( ica .eq. 0 .or. icb .eq. 0 ) goto 981
                     if( icb-ica .lt. 0 .or. icb-ica .gt. 1 ) goto 981

                     cnuc = chlw(ica:icb)//'  '

                     do j = 1, 104

                        if( cnuc(1:3) .eq. element(j)(1:3) ) then

                           icha = j

                           goto 450

                        end if

                     end do

                           goto 981

  450                continue

                     if( icha .gt. 104 )  goto 981

                  if( icm .eq. 0 .or. icn .eq. 0 ) then

                     isnuc_temporary(nnucl-1+k) = icha * 1000

                  else

                     if( isn .gt. 3 ) goto 981

                     read(chlw(icm:icn),'(i5)') masi

                     if( masi .lt. icha ) goto 981
                     if( masi-icha .gt. maxnt ) goto 981
                     if( masi .eq. icha ) goto 975

                     isnuc_temporary(nnucl-1+k) = icha * 1000 + masi

                     imasi = imasi + 1

                  end if

                     ic = icd + 1

               end do
               call moddas_deallocate_int(isnuc_temporary)

            else

               goto 981

            end if

*-----------------------------------------------------------------------
*        material
*-----------------------------------------------------------------------

         else if( ipm .eq. 13 ) then

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
*        axis
*-----------------------------------------------------------------------

         else if( ipm .eq. 5 ) then

  600       inaxi = inaxi + 1

            if( inaxi .gt. 6 ) goto 991

            if( chlw(ic:ic+3) .eq. 'mass' ) then

               iaxis(inaxi) = 1
               ic =jnumc(chlw,ic+5,icl)

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               iaxis(inaxi) = 2
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+2) .eq. 'tet' ) then

               iaxis(inaxi) = 14
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+5) .eq. 'dchain' ) then

               iaxis(inaxi) = 13
               ic =jnumc(chlw,ic+7,icl)

            else if( chlw(ic:ic+1) .eq. 'xy' .or.
     &               chlw(ic:ic+1) .eq. 'yx' ) then

               iaxis(inaxi) = 9
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic+1) .eq. 'yz' .or.
     &               chlw(ic:ic+1) .eq. 'zy' ) then

               iaxis(inaxi) = 10
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic+1) .eq. 'zx' .or.
     &               chlw(ic:ic+1) .eq. 'xz' ) then

               iaxis(inaxi) = 11
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic+1) .eq. 'zr' .or.
     &               chlw(ic:ic+1) .eq. 'rz' ) then

               iaxis(inaxi) = 12
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

            else if( chlw(ic:ic+5) .eq. 'charge' ) then

               iaxis(inaxi) = 7
               ic =jnumc(chlw,ic+7,icl)

            else if( chlw(ic:ic+4) .eq. 'chart' ) then

               iaxis(inaxi) = 8
               ic =jnumc(chlw,ic+6,icl)

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
         else if( ipm .eq. 36 ) then

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

         else if( ipm .eq. 14 ) then

               ict = min( ic + 199, i2 )

               cxtxt = chin(ic:ict)

               lxtxt = ict - ic + 1

*-----------------------------------------------------------------------
*        y-txt
*-----------------------------------------------------------------------

         else if( ipm .eq. 15 ) then

               ict = min( ic + 199, i2 )

               cytxt = chin(ic:ict)

               lytxt = ict - ic + 1

*-----------------------------------------------------------------------
*        z-txt
*-----------------------------------------------------------------------

         else if( ipm .eq. 16 ) then

               ict = min( ic + 199, i2 )

               cztxt = chin(ic:ict)

               lztxt = ict - ic + 1

*-----------------------------------------------------------------------
*        gshow
*-----------------------------------------------------------------------

         else if( ipm .eq. 17 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               lgshow = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        gslat
*-----------------------------------------------------------------------

         else if( ipm .eq. 34 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               igslt = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        output
*-----------------------------------------------------------------------

         else if( ipm .eq. 35 ) then

            if( chlw(ic:ic+6) .eq. 'product' ) then

               iout = 0

            else if( chlw(ic:ic+5) .eq. 'cutoff' ) then

               iout = 1

            else

               goto 997

            end if

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        rshow
*-----------------------------------------------------------------------

         else if( ipm .eq. 18 ) then

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
*        transform
*-----------------------------------------------------------------------

         else if( ipm .eq. 31 .or. ipm .eq. 32 ) then

                  if( ipm .eq. 32 ) ktrs = 1

               call ttrans(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     ic,ktrs,igkst,idtt,vtrs)


               goto 150

*-----------------------------------------------------------------------
*        reg echo length
*-----------------------------------------------------------------------

         else if( ipm .eq. 20 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               iechrl = nint( cvvv )

               if( iechrl .lt. 40 ) goto 997

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        unit
*-----------------------------------------------------------------------

         else if( ipm .eq. 9 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               iunt = nint( cvvv )

               if( iunt .lt. 1 .or. iunt .gt. 2 ) goto 980

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        info
*-----------------------------------------------------------------------

         else if( ipm .eq. 10 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               info = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        2d-type
*-----------------------------------------------------------------------

         else if( ipm .eq. 11 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               idtyp = nint( cvvv )

               if( idtyp .lt. 1 .or. idtyp .gt. 7 ) goto 983

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        factor
*-----------------------------------------------------------------------

         else if( ipm .eq. 12 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               rfact = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        volmat
*-----------------------------------------------------------------------

         else if( ipm .eq. 21 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               matvol = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        epsout
*-----------------------------------------------------------------------

         else if( ipm .eq. 22 ) then

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
*        counter
*-----------------------------------------------------------------------

         else if( ipm .ge. 23 .and. ipm .le. 28 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               icount((ipm-23+2)/2) = 1
               icount( ipm-23+4)    = nint( cvvv )

               if( abs( nint(cvvv) ) .gt. 9999 ) goto 976

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        ndata
*-----------------------------------------------------------------------

         else if( ipm .eq. 19 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ndata = nint( cvvv )

               if( ndata .lt. 0 .or. ndata .gt. 3 ) goto 978

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

      else if ( itmsh(m) .eq. 4 ) then

        call check_tet(m,MAX_NUM_ITREG,mtetreg,
     &      iec,cepn,lepn,ierr)
        call moddas_deallocate_int(mtetreg)

      end if

      call check_unit(m,iunt,iec,cepn,lepn,ierr)

      call check_axis(m,iaxis,inaxi,iec,cepn,lepn,ierr)

      call check_output(m,iout,iec,cepn,lepn,ierr)

      do i = 1, itaxn(m)

        if ( any( itaxs(m,i) .eq. (/ 8, 9, 10, 11, 12 /))) then
          call check_2dtype(m,idtyp,iec,cepn,lepn,ierr)
        end if

      end do

      call check_part(m,iptyp,ipnkf,inpat,iec,cepn,lepn,ierr)

*-----------------------------------------------------------------------

      call judge_tall_check(iec,cepn,lepn,dsin(jsn),idsi(jsn),ierr)

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  976    m_err = 'counter should be from -9999 to 9999'
         ErrCha = ''
         ErrID = 'L:1385/R:check_tyield/F:restyield.f'
         goto 999

  950    write(dkam,'(i9)') mdas
         m_err = 'Total tally storage number exceeds mdas ='//dkam
         ErrCha = ''
         ErrID = 'L:1391/R:check_tyield/F:restyield.f'
         goto 999

  980    m_err = 'unit should be 1 or 2 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1396/R:check_tyield/F:restyield.f'
         goto 999

  978    m_err = 'ndata should be 0 or 1 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1401/R:check_tyield/F:restyield.f'
         goto 999

  979    m_err = 'Description of material parameter is wrong in tally '
     &            //tname
         ErrCha = ''
         ErrID = 'L:1407/R:check_tyield/F:restyield.f'
         goto 999

  975    m_err = '1H or Z=A cannot be specified for nucleus in tally '
     &            //tname
         ErrCha = ''
         ErrID = 'L:1413/R:check_tyield/F:restyield.f'
         goto 999

  981    m_err = 'Description of nucleus parameter is wrong in tally '
     &            //tname
         ErrCha = ''
         ErrID = 'L:1419/R:check_tyield/F:restyield.f'
         goto 999

  982    m_err = 'Description of mother parameter is wrong in tally '
     &            //tname
         ErrCha = ''
         ErrID = 'L:1425/R:check_tyield/F:restyield.f'
         goto 999

  983    m_err = '2D-type should be 1-7 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1430/R:check_tyield/F:restyield.f'
         goto 999

  986    write(dkam,'(i5)') itlmax
         m_err = 'Total tally number exceeds itlmax = '//dkam
         ErrCha = ''
         ErrID = 'L:1436/R:check_tyield/F:restyield.f'
         goto 999

  987    m_err = 'Unknown parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1441/R:check_tyield/F:restyield.f'
         goto 999

  988    m_err = 'special should be positive in tally '//tname
         ErrCha = ''
         ErrID = 'L:1446/R:check_tyield/F:restyield.f'
         goto 999

  989    m_err = 'Double definition of the parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1451/R:check_tyield/F:restyield.f'
         goto 999

  990    m_err = 'Too many file name in tally '//tname
         ErrCha = ''
         ErrID = 'L:1456/R:check_tyield/F:restyield.f'
         goto 999

  991    m_err = 'Too many axis in tally '//tname
         ErrCha = ''
         ErrID = 'L:1461/R:check_tyield/F:restyield.f'
         goto 999

  992    m_err = 'Unknown axis name in tally '//tname
         ErrCha = ''
         ErrID = 'L:1466/R:check_tyield/F:restyield.f'
         goto 999

  993    m_err = 'In this line, [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:1471/R:check_tyield/F:restyield.f'
         goto 999

  994    m_err = 'Name of particle is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1476/R:check_tyield/F:restyield.f'
         goto 999

  995    m_err = 'Number of particles is larger than 6 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1481/R:check_tyield/F:restyield.f'
         goto 999

  996    m_err = 'Description of mesh is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1486/R:check_tyield/F:restyield.f'
         goto 999

  997    m_err = 'Description of parameter is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1491/R:check_tyield/F:restyield.f'
         goto 999

  998    m_err = 'Unknown mesh parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1496/R:check_tyield/F:restyield.f'
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
      subroutine read_tyield(m,iax,ierr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************

        use RESTALMOD
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
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
        common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                  itnun(itlmax), itnuc(itlmax),
     &                  itndz(itlmax), itndn(itlmax),
     &                  itnkz(itlmax), itnkn(itlmax)
        common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)
        common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
        common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)
        common /tall83/ itnzn(itlmax), itndm(itlmax) ! frtati 2022/02/18

        dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
        equivalence ( das, idas )

*-----------------------------------------------------------------------

        character dsin(0:9,2)*200
        dimension idsi(0:9,2)
        dimension ill(0:9,2), ilf(0:9,2)
        dimension jsn(2), jsi(2)
        dimension ierrs(2)

*-----------------------------------------------------------------------
        common /redufmt/ iredufmt(itlmax) !FURUTA20200615
*-----------------------------------------------------------------------

        ierr = 0
        noe  = 1

        if ( any( itaxs(m,iax) .eq. (/ 8, 9, 10, 11, 12, 13 /) )
     &      .and. ittwo(m) .ne. 4 ) noe = 2

        if(itaxs(m,iax).eq.13.and.iredufmt(m).eq.1) noe=1 !FURUTA20200615

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
        call check_tyield(m,iax,jsn(ioe),jsi(ioe),
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
        idas2 = idas1 + ( maxnt + maxpt ) * 2
        if( itmsh(m) .eq. 1 ) then
          idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
          idas4 = idas2 + itrgn(m)
     &          + ( itrgn(m) + mod(itrgn(m),2) ) / 2
          idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
        else if ( itmsh(m) .eq. 3 ) then
          idas3 = ( idas2 + itnfn(m) * itnfn(m) - 1 ) * 2 + 1
        end if

*-----------------------------------------------------------------------
*   read tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call read_yieldreg(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndz(m),itndn(m),itndm(m),
     &      itrgn(m),itrgm(m),itnun(m),
     &      idas_itreg(itreg(m)),
     &      isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &      trRES(irestalm(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        else if ( itmsh(m) .eq. 2 ) then

          call read_yieldrz(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndz(m),itndn(m),itndm(m),itnfr(m),itnfz(m),
     &      itrnm(m),itznm(m),itnun(m),
     &      das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &      isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &      trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 3 ) then

          call read_yieldxyz(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndz(m),itndn(m),itndm(m),itnfn(m),
     &      itmtn(m),ismte(itmtt(m)),
     &      itxnm(m),itynm(m),itznm(m),itnun(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)),
     &      isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &      trRES(irestalm(m)))

        else if( itmsh(m) .eq. 4 ) then

          call read_yieldtet(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndz(m),itndn(m),itndm(m),itrgn(m),itnun(m),
     &      isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &      trRES(irestalm(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        end if

  800   continue
*-----------------------------------------------------------------------
*   restore tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call restore_yieldreg(m,
     &      itndz(m),itndn(m),itndm(m),
     &      itrgn(m),itrgm(m),itnun(m),
     &      idas_itreg(itreg(m)),
     &      isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &      trRES(irestalm(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        else if ( itmsh(m) .eq. 2 ) then

          call restore_yieldrz(m,
     &      itndz(m),itndn(m),itndm(m),itnfr(m),itnfz(m),
     &      itrnm(m),itznm(m),itnun(m),
     &      das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &      isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &      trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 3 ) then

          call restore_yieldxyz(m,
     &      itndz(m),itndn(m),itndm(m),itnfn(m),
     &      itmtn(m),ismte(itmtt(m)),
     &      itxnm(m),itynm(m),itznm(m),itnun(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)),
     &      isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &      trRES(irestalm(m)))

        else if( itmsh(m) .eq. 4 ) then

          call restore_yieldtet(m,
     &      itndz(m),itndn(m),itndm(m),
     &      itrgn(m),itrgm(m),itnun(m),
     &      isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &      trRES(irestalm(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      idas_itreg(itreg(m)),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

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
      subroutine read_yieldreg(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                         mz,mn,mm, ! frtati 2022/02/18 added mm
     &                         nr,mr,nn,kr,nt,ikzz,iknn,
     &                         tr,nvl,ivl,rvl,
     &                         nx,ny,nz,xm,ym,zm)
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
        common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
        common /tall83/ itnzn(itlmax), itndm(itlmax) ! frtati 2022/02/18

        common /tall21/ rtfac(itlmax)
        common /fact01/ facmax(itlmax)
        common /stat / istdev, irestart, ireschk
        character ctfln*100

        dimension   kr(mr)
        dimension   nt(nn)
        dimension   ikzz(maxpt,maxnt), iknn(maxpt,maxnt)
        dimension   tm(maxnt+maxpt,2)
        dimension   vl(nr)
        dimension   lr(nr)
        dimension   tr(nr,mz,mn,0:mm,2) ! frtati 2022/02/18 2 -> mm
        dimension   ivl(nvl)
        dimension   rvl(nvl)

        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)
        dimension   val(nr)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)
        dimension a0_tr(maxnt) ! frtati 2022/02/18

        character chin*200, chlw*200, chcm*200, ctem*2

*-----------------------------------------------------------------------

        character dc2*200

*-----------------------------------------------------------------------

        character elmnt(104)*3

        data elmnt /
     &    'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ',
     &    'Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ',
     &    'Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

*-----------------------------------------------------------------------
        common /redufmt/ iredufmt(itlmax) !FURUTA20200615
*-----------------------------------------------------------------------
*        dchain axis
*-----------------------------------------------------------------------
         if( itaxs(m,iax) .eq. 13 ) then

          if(iredufmt(m).eq.0)then !FURUTA20200615

              !! seek to next isotope production
              jpn = 0
  110         call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$#',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
              if ( ierr  .ne. 0  ) goto 900
              if ( jpn   .eq. 3  ) goto 900

              if ( iskip .ne. 0  ) goto 110
              if ( i2    .lt. 27 ) goto 110

              if (chin(9:27) .eq. ' isotope production') then
               il = 0
              elseif (chin(9:27) .eq. ' 1st metastable iso') then
               il = 1
              elseif (chin(9:27) .eq. ' 2nd metastable iso') then
               il = 2
              else
               goto 110
              endif

*-----------------------------------------------------------------------

              read(chin(1:5),'(i5)') iz
              IF(il .eq. 0) then
               read(chin(39:45),'(i3,1x,i3)') n3, n4
              ELSE
               read(chin(54:60),'(i3,1x,i3)') n3, n4
              ENDIF


              call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! reg

              do ir = 1, nr

               if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                 read(jsi,'(i7,1p12e11.3)')
     &                lr(ir), ( tr(ir,iz,i,il,ioe), i = n3, n4 )
               else
                 read(jsi,'(i7,1p12e11.3)')
     &           lr(ir),  ( a0_tr(i), i = n3, n4 )
                 do i = n3, n4
                   if( a0_tr(i).ne.0.d0 ) then
                     tr(ir,igetiznm(iz,i,il,m),1,0,ioe) = a0_tr(i)
                   end if
                 end do
               end if

              end do

              call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$#',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)



           goto 110

*-----------------------------------------------------------------------
          else

           call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#numnucleusidyieldr.err',13,ierr)

           if ( ierr  .ne. 0 ) goto 900
           if ( jpn   .eq. 3 ) goto 900

           do
            read(jsi,*)ir,nucleusID,yield,e_yield
            if(ir.eq.0)exit
            iz=nucleusID/10000
            ia=(nucleusID-iz*10000)/10
            in=ia-iz
            il=mod(nucleusID,10)
            if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
             tr(ir,iz,in,il,1)=yield
             tr(ir,iz,in,il,2)=e_yield
            else
             tr(ir,igetiznm(iz,in,il,m),1,0,1)=yield
             tr(ir,igetiznm(iz,in,il,m),1,0,2)=e_yield
            end if
           enddo

          endif
*-----------------------------------------------------------------------

  190 continue

*-----------------------------------------------------------------------
*        reg axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 2 ) then

           dc2  = '#numregvolume'
           ldc2 = 13

           if( nn .eq. 0 ) then
              nc = 1
           else
              nc = nn
           end if

          facmx = 1.d0  ! kitamura23/03/31

           do 290 ic = 1, nc

             if( nn .gt. 0 ) then
                iz = nt(ic) / 1000
                ia = nt(ic) - iz * 1000
                if( ia .gt. 0 ) then
                   in = ia - iz
                end if
             end if

*-----------------------------------------------------------------------

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             do ir = 1, nr

               read(jsi,'(i5,1x,i7,1pe13.4,1pe13.4,0pf8.4)')
     &              idmm0, idmm1, dmm2, sek, ser

               if( nn .eq. 0 ) then

                  do kn = 1, mn
                  do kz = 1, mz
                  do il = 0, mm ! frtati 2022/02/18 2 -> mm

                     lz  = ikzz( kz, kn )
                     ln  = iknn( kz, kn )

                     IF(il .eq. 0) then
                      tr(ir,lz,ln,il,1) = sek
                      tr(ir,lz,ln,il,2) = ser
                     ELSE
                      tr(ir,lz,ln,il,1) = 0.d0
                      tr(ir,lz,ln,il,2) = 0.d0
                     ENDIF

                  end do
                  end do
                  end do

               else if( ia .eq. 0 ) then

                  do kn = 1, mn
                  do il = 0, mm ! frtati 2022/02/18 2 -> mm

                     lz  = ikzz( iz, kn )
                     ln  = iknn( iz, kn )

                     IF(il .eq. 0) then
                      tr(ir,lz,ln,il,1) = sek
                      tr(ir,lz,ln,il,2) = ser
                     ELSE
                      tr(ir,lz,ln,il,1) = 0.d0
                      tr(ir,lz,ln,il,2) = 0.d0
                     ENDIF

                  end do
                  end do

               else

                  do il = 0, mm ! frtati 2022/02/18 2 -> mm
                     lz  = ikzz( iz, in )
                     ln  = iknn( iz, in )

                     IF(il .eq. 0) then
                      tr(ir,lz,ln,il,1) = sek
                      tr(ir,lz,ln,il,2) = ser
                     ELSE
                      tr(ir,lz,ln,il,1) = 0.d0
                      tr(ir,lz,ln,il,2) = 0.d0
                     ENDIF
                  end do

               end if

             end do

  290     continue

*-----------------------------------------------------------------------
*        mass axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 1 ) then

           dc2  = '#massnumberr.err'
           ldc2 = 16

           if( nn .eq. 0 ) then
              nc = 1
           else
              nc = nn
           end if

           facmx = 1.d0  ! kitamura23/03/31

           do 390 ir = 1, nr
           do 390 ic = 1, nc

             if( nn .gt. 0 ) iz = nt(ic) / 1000

             do i = 1, maxnt + maxpt

                tm(i,1) = 0.0d+0
                tm(i,2) = 0.0d+0

             end do

*-----------------------------------------------------------------------

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
             read(chin(11:i2),'(i3,1x,i3)') im, jm
             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             if( im .eq. 1 ) then

                read(jsi,'(3x,f4.1,2x,1pe13.4,0pf8.4)') dmm0, dmm1, dmm2

             end if

             do i = im, jm

               read(jsi,'(3x,i4,2x,1pe13.4,0pf8.4)')
     &                       idmm0, tm(i,1), tm(i,2)

             end do

*-----------------------------------------------------------------------

            if( nn .eq. 0 .and. itfln(m) .eq. 1 ) then

               do ln = 1, mn
               do il = 0, mm ! frtati 2022/02/18 2 -> mm

                  lz = 1
                  im = ln

                  IF(il .eq. 0) then
                   tr(ir,lz,ln,il,1) = tm(im,1)
                   tr(ir,lz,ln,il,2) = tm(im,2)
                  ELSE
                   tr(ir,lz,ln,il,1) = 0.d0
                   tr(ir,lz,ln,il,2) = 0.d0
                  ENDIF

               end do
               end do

            else if( nn .eq. 0 .and. itfln(m) .gt. 1 ) then

               do kn = 1, mn
               do kz = 1, mz
               do il = 0, mm ! frtati 2022/02/18 2 -> mm

                  lz = ikzz( kz, kn )
                  ln = iknn( kz, kn )

                  im = kz + kn

                  IF(il .eq. 0) then
                   tr(ir,lz,ln,il,1) = tm(im,1)
                   tr(ir,lz,ln,il,2) = tm(im,2)
                  ELSE
                   tr(ir,lz,ln,il,1) = 0.d0
                   tr(ir,lz,ln,il,2) = 0.d0
                  ENDIF

               end do
               end do
               end do

            else

               do kn = 1, mn
               do il = 0, mm ! frtati 2022/02/18 2 -> mm

                  lz = ikzz( iz, kn )
                  ln = iknn( iz, kn )

                  im = iz + kn

                  IF(il .eq. 0) then
                   tr(ir,lz,ln,il,1) = tm(im,1)
                   tr(ir,lz,ln,il,2) = tm(im,2)
                  ELSE
                   tr(ir,lz,ln,il,1) = 0.d0
                   tr(ir,lz,ln,il,2) = 0.d0
                  ENDIF

               end do
               end do

            end if

  390     continue

*-----------------------------------------------------------------------
*        charge axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 7 ) then

           dc2  = '#chargenumberr.err'
           ldc2 = 18

           facmx = 1.d0  ! kitamura23/03/31

           do 490 ir = 1, nr

             do i = 1, maxpt

                tm(i,1) = 0.0d+0
                tm(i,2) = 0.0d+0

             end do

*-----------------------------------------------------------------------

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
             read(chin(11:i2),'(i3,1x,i3)') im, jm
             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             if( im .eq. 1 ) then

                read(jsi,'(3x,f4.1,2x,1pe13.4,0pf8.4)') dmm0, dmm1, dmm2

             end if

             do i = im, jm

               read(jsi,'(3x,i4,2x,1pe13.4,0pf8.4)')
     &                       idmm0, tm(i,1), tm(i,2)

             end do

*-----------------------------------------------------------------------

            do kz = 1, mz
            do kn = 1, mn
            do il = 0, mm ! frtati 2022/02/18 2 -> mm

                  iz  = ikzz( kz, kn )
                  in  = iknn( kz, kn )

                  IF(il .eq. 0) then
                   tr(ir,iz,in,il,1) = tm(iz,1)
                   tr(ir,iz,in,il,2) = tm(iz,2)
                  ELSE
                   tr(ir,iz,in,il,1) = 0.d0
                   tr(ir,iz,in,il,2) = 0.d0
                  ENDIF

            end do
            end do
            end do

  490     continue

*-----------------------------------------------------------------------
*        chart axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 8 ) then

           if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then
             dc2 = 'hc:y='
             ldc2 = 5
           else if( ittwo(m) .eq. 4 ) then
             dc2 = '#znmassnumberr.err'
             ldc2 = 18
           else if( ittwo(m) .eq. 5 ) then
             dc2 = 'z/n'
             ldc2 = 3
           end if

           if( itnzn(m).ne.0 ) tr(:,:,1,0,ioe) = 0.d0 ! frtati 2022/03/18

*-----------------------------------------------------------------------

           facmx = 1.d0  ! kitamura23/03/31

           do 590

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)
             if( jsn.eq.-1 ) goto 900 ! T.Sato 2023/04/13, no more data
            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
             read(chin(16:17),'(A2)') ctem ! after reading all data, it goes to '# icmax...' line
             IF(ctem .ne. 'ir') goto 900
             read(chin(21:24),'(i3)') ir
             read(chin(34:34),'(i1)') il
             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
             read(chin(17:i2),'(i3,1x,i3)') icmax, inmax
             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)
             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

*-----------------------------------------------------------------------
             if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               do i = icmax+2, 1, -1

                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                 read(jsi,'(1p10e11.3)')
     &             ( tr(ir,i,l,il,ioe), l = 1, inmax+2 )
                else
                 read(jsi,'(1p10e11.3)')
     &             ( a0_tr(l), l = 1, inmax+2 )
                 do l = 1, inmax+2
                   if( a0_tr(l).ne.0.d0 ) then
                      tr(ir,igetiznm(i,l,il,m),1,0,ioe) = a0_tr(l)
                   end if
                 end do
                end if

               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 4 ) then

               do i0 = 1, icmax
               do j0 = 1, inmax

                 call readl(jsn,jsi,dsin,idsi,ill,ilf,'',
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
                 if ( ierr .ne. 0  ) goto 590
                 if ( jpn  .eq. 3  ) goto 590
                 if ( i1   .eq. i3 ) goto 590

                 read(chin(1:i2),'(3i5,1pe13.4,0pf8.4)')
     &                i, j, ij, A, B

                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                 tr(ir,i,j,il,1) = A
                 tr(ir,i,j,il,2) = B
                else if( A.ne.0.d0 ) then
                 tr(ir,igetiznm(i,j,il,m),1,0,1) = A
                 tr(ir,igetiznm(i,j,il,m),1,0,2) = B
                end if

               end do
               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 5 ) then

              if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
               do i = icmax+2, 1, -1
                  read(jsi,'(1p1000e11.3)')
     &            dmm0,
     &            ( tr(ir,i,l,il,ioe), l = 1, inmax+2 )

               end do
              else
               do i = icmax+2, 1, -1
                  read(jsi,'(1p1000e11.3)') dmm0,
     &            ( a0_tr(l), l = 1, inmax+2 )
                 do l = 1, inmax+2
                   if( a0_tr(l).ne.0.d0 ) then
                     tr(ir,igetiznm(i,l,il,m),1,0,ioe) = a0_tr(l)
                   end if
                 end do
               end do
              end if

*-----------------------------------------------------------------------
             end if

  590     continue

*-----------------------------------------------------------------------
         end if

  900 continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_yieldtet(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                         mz,mn,mm, ! frtati 2022/02/18 added mm
     &                         nr,nn,nt,ikzz,iknn,
     &                         tr,
     &                         nx,ny,nz,xm,ym,zm)
*                                                                      *
*   m: the tally number, index of ital.                                *
*                                                                      *
*     Created by T.Furuta on 2019/01/21
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'
        include 'param01.inc'

*-----------------------------------------------------------------------

        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
        common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
        common /tall83/ itnzn(itlmax), itndm(itlmax) ! frtati 2022/02/18

        common /tall21/ rtfac(itlmax)
        common /fact01/ facmax(itlmax)
        common /stat / istdev, irestart, ireschk
        character ctfln*100

        dimension   nt(nn)
        dimension   ikzz(maxpt,maxnt), iknn(maxpt,maxnt)
        dimension   tm(maxnt+maxpt,2)
        dimension   vl(nr)
        dimension   lr(nr)
        dimension   tr(nr,mz,mn,0:mm,2) ! frtati 2022/02/18 2 -> mm

        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)
        dimension   val(nr)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)
        dimension a0_tr(maxnt) ! frtati 2022/02/18

        character chin*200, chlw*200, chcm*200, ctem*2

*-----------------------------------------------------------------------

        character dc2*200

*-----------------------------------------------------------------------

        character elmnt(104)*3

        data elmnt /
     &    'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ',
     &    'Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ',
     &    'Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

*-----------------------------------------------------------------------
        common /redufmt/ iredufmt(itlmax) !FURUTA20200615
*-----------------------------------------------------------------------
*        dchain axis
*-----------------------------------------------------------------------
         if( itaxs(m,iax) .eq. 13 ) then

          if(iredufmt(m).eq.0)then !FURUTA20200615

              !! seek to next isotope production
              jpn = 0
  110         call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$#',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
              if ( ierr  .ne. 0  ) goto 900
              if ( jpn   .eq. 3  ) goto 900

              if ( iskip .ne. 0  ) goto 110
              if ( i2    .lt. 27 ) goto 110

              if (chin(9:27) .eq. ' isotope production') then
               il = 0
              elseif (chin(9:27) .eq. ' 1st metastable iso') then
               il = 1
              elseif (chin(9:27) .eq. ' 2nd metastable iso') then
               il = 2
              else
               goto 110
              endif

*-----------------------------------------------------------------------

              read(chin(1:5),'(i5)') iz
              IF(il .eq. 0) then
               read(chin(39:45),'(i3,1x,i3)') n3, n4
              ELSE
               read(chin(54:60),'(i3,1x,i3)') n3, n4
              ENDIF


              call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! reg

              do ir = 1, nr

               if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                 read(jsi,'(i7,1p12e11.3)')
     &                lr(ir), ( tr(ir,iz,i,il,ioe), i = n3, n4 )
               else
                 read(jsi,'(i7,1p12e11.3)') lr(ir),
     &           ( tr(ir,igetiznm(iz,i,il,m),1,0,ioe), i = n3, n4 )
               end if

              end do

              call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$#',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)



           goto 110

*-----------------------------------------------------------------------
          else

           call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#numnucleusidyieldr.err',13,ierr)

           if ( ierr  .ne. 0 ) goto 900
           if ( jpn   .eq. 3 ) goto 900

           do
            read(jsi,*)ir,nucleusID,yield,e_yield
            if(ir.eq.0)exit
            iz=nucleusID/10000
            ia=(nucleusID-iz*10000)/10
            in=ia-iz
            il=mod(nucleusID,10)
            if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
             tr(ir,iz,in,il,1)=yield
             tr(ir,iz,in,il,2)=e_yield
            else
             tr(ir,igetiznm(iz,in,il,m),1,0,1)=yield
             tr(ir,igetiznm(iz,in,il,m),1,0,2)=e_yield
            end if
           enddo

          endif
*-----------------------------------------------------------------------

  190 continue

*-----------------------------------------------------------------------
*        tet axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 14 ) then

           dc2  = '#numtetravolume'
           ldc2 = 13

           if( nn .eq. 0 ) then
              nc = 1
           else
              nc = nn
           end if

          facmx = 1.d0  ! kitamura23/03/31

           do 290 ic = 1, nc

             if( nn .gt. 0 ) then
                iz = nt(ic) / 1000
                ia = nt(ic) - iz * 1000
                if( ia .gt. 0 ) then
                   in = ia - iz
                end if
             end if

*-----------------------------------------------------------------------

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             do ir = 1, nr

               read(jsi,'(i8,1x,i8,1pe13.4,1pe13.4,0pf8.4)')
     &              idmm0, idmm1, dmm2, sek, ser

               if( nn .eq. 0 ) then

                  do kn = 1, mn
                  do kz = 1, mz
                  do il = 0, mm ! frtati 2022/02/18 2 -> mm

                     lz  = ikzz( kz, kn )
                     ln  = iknn( kz, kn )

                     IF(il .eq. 0) then
                      tr(ir,lz,ln,il,1) = sek
                      tr(ir,lz,ln,il,2) = ser
                     ELSE
                      tr(ir,lz,ln,il,1) = 0.d0
                      tr(ir,lz,ln,il,2) = 0.d0
                     ENDIF

                  end do
                  end do
                  end do

               else if( ia .eq. 0 ) then

                  do kn = 1, mn
                  do il = 0, mm ! frtati 2022/02/18 2 -> mm

                     lz  = ikzz( iz, kn )
                     ln  = iknn( iz, kn )

                     IF(il .eq. 0) then
                      tr(ir,lz,ln,il,1) = sek
                      tr(ir,lz,ln,il,2) = ser
                     ELSE
                      tr(ir,lz,ln,il,1) = 0.d0
                      tr(ir,lz,ln,il,2) = 0.d0
                     ENDIF

                  end do
                  end do

               else

                  do il = 0, mm ! frtati 2022/02/18 2 -> mm
                     lz  = ikzz( iz, in )
                     ln  = iknn( iz, in )

                     IF(il .eq. 0) then
                      tr(ir,lz,ln,il,1) = sek
                      tr(ir,lz,ln,il,2) = ser
                     ELSE
                      tr(ir,lz,ln,il,1) = 0.d0
                      tr(ir,lz,ln,il,2) = 0.d0
                     ENDIF
                  end do

               end if

             end do

  290     continue

*-----------------------------------------------------------------------
*        mass axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 1 ) then

           dc2  = '#massnumberr.err'
           ldc2 = 16

           if( nn .eq. 0 ) then
              nc = 1
           else
              nc = nn
           end if

           facmx = 1.d0  ! kitamura23/03/31

           do 390 ir = 1, nr
           do 390 ic = 1, nc

             if( nn .gt. 0 ) iz = nt(ic) / 1000

             do i = 1, maxnt + maxpt

                tm(i,1) = 0.0d+0
                tm(i,2) = 0.0d+0

             end do

*-----------------------------------------------------------------------

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
             read(chin(11:i2),'(i3,1x,i3)') im, jm
             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             if( im .eq. 1 ) then

                read(jsi,'(3x,f4.1,2x,1pe13.4,0pf8.4)') dmm0, dmm1, dmm2

             end if

             do i = im, jm

               read(jsi,'(3x,i4,2x,1pe13.4,0pf8.4)')
     &                       idmm0, tm(i,1), tm(i,2)

             end do

*-----------------------------------------------------------------------

            if( nn .eq. 0 .and. itfln(m) .eq. 1 ) then

               do ln = 1, mn
               do il = 0, mm ! frtati 2022/02/18 2 -> mm

                  lz = 1
                  im = ln

                  IF(il .eq. 0) then
                   tr(ir,lz,ln,il,1) = tm(im,1)
                   tr(ir,lz,ln,il,2) = tm(im,2)
                  ELSE
                   tr(ir,lz,ln,il,1) = 0.d0
                   tr(ir,lz,ln,il,2) = 0.d0
                  ENDIF

               end do
               end do

            else if( nn .eq. 0 .and. itfln(m) .gt. 1 ) then

               do kn = 1, mn
               do kz = 1, mz
               do il = 0, mm ! frtati 2022/02/18 2 -> mm

                  lz = ikzz( kz, kn )
                  ln = iknn( kz, kn )

                  im = kz + kn

                  IF(il .eq. 0) then
                   tr(ir,lz,ln,il,1) = tm(im,1)
                   tr(ir,lz,ln,il,2) = tm(im,2)
                  ELSE
                   tr(ir,lz,ln,il,1) = 0.d0
                   tr(ir,lz,ln,il,2) = 0.d0
                  ENDIF

               end do
               end do
               end do

            else

               do kn = 1, mn
               do il = 0, mm ! frtati 2022/02/18 2 -> mm

                  lz = ikzz( iz, kn )
                  ln = iknn( iz, kn )

                  im = iz + kn

                  IF(il .eq. 0) then
                   tr(ir,lz,ln,il,1) = tm(im,1)
                   tr(ir,lz,ln,il,2) = tm(im,2)
                  ELSE
                   tr(ir,lz,ln,il,1) = 0.d0
                   tr(ir,lz,ln,il,2) = 0.d0
                  ENDIF

               end do
               end do

            end if

  390     continue

*-----------------------------------------------------------------------
*        charge axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 7 ) then

           dc2  = '#chargenumberr.err'
           ldc2 = 18

           facmx = 1.d0  ! kitamura23/03/31

           do 490 ir = 1, nr

             do i = 1, maxpt

                tm(i,1) = 0.0d+0
                tm(i,2) = 0.0d+0

             end do

*-----------------------------------------------------------------------

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
             read(chin(11:i2),'(i3,1x,i3)') im, jm
             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             if( im .eq. 1 ) then

                read(jsi,'(3x,f4.1,2x,1pe13.4,0pf8.4)') dmm0, dmm1, dmm2

             end if

             do i = im, jm

               read(jsi,'(3x,i4,2x,1pe13.4,0pf8.4)')
     &                       idmm0, tm(i,1), tm(i,2)

             end do

*-----------------------------------------------------------------------

            do kz = 1, mz
            do kn = 1, mn
            do il = 0, mm ! frtati 2022/02/18 2 -> mm

                  iz  = ikzz( kz, kn )
                  in  = iknn( kz, kn )

                  IF(il .eq. 0) then
                   tr(ir,iz,in,il,1) = tm(iz,1)
                   tr(ir,iz,in,il,2) = tm(iz,2)
                  ELSE
                   tr(ir,iz,in,il,1) = 0.d0
                   tr(ir,iz,in,il,2) = 0.d0
                  ENDIF

            end do
            end do
            end do

  490     continue

*-----------------------------------------------------------------------
*        chart axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 8 ) then

           if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then
             dc2 = 'hc:y='
             ldc2 = 5
           else if( ittwo(m) .eq. 4 ) then
             dc2 = '#znmassnumberr.err'
             ldc2 = 18
           else if( ittwo(m) .eq. 5 ) then
             dc2 = 'z/n'
             ldc2 = 3
           end if

           if( itnzn(m).ne.0 ) tr(:,:,1,0,ioe) = 0.d0 ! frtati 2022/03/18

*-----------------------------------------------------------------------

           facmx = 1.d0  ! kitamura23/03/31

           do 590

             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
             read(chin(16:17),'(A2)') ctem ! after reading all data, it goes to '# icmax...' line
             IF(ctem .ne. 'ir') goto 900
             read(chin(21:24),'(i3)') ir
             read(chin(34:34),'(i1)') il
             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
             read(chin(17:i2),'(i3,1x,i3)') icmax, inmax
             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)
             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

*-----------------------------------------------------------------------
             if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               do i = icmax+2, 1, -1

                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                 read(jsi,'(1p10e11.3)')
     &             ( tr(ir,i,l,il,ioe), l = 1, inmax+2 )
                else
                 read(jsi,'(1p10e11.3)')
     &             ( a0_tr(l), l = 1, inmax+2 )
                 do l = 1, inmax+2
                   if( a0_tr(l).ne.0.d0 ) then
                      tr(ir,igetiznm(i,l,il,m),1,0,ioe) = a0_tr(l)
                   end if
                 end do
                end if

               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 4 ) then

               do i0 = 1, icmax
               do j0 = 1, inmax

                 call readl(jsn,jsi,dsin,idsi,ill,ilf,'',
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
                 if ( ierr .ne. 0  ) goto 590
                 if ( jpn  .eq. 3  ) goto 590
                 if ( i1   .eq. i3 ) goto 590

                 read(chin(1:i2),'(3i5,1pe13.4,0pf8.4)')
     &                i, j, ij, A, B

                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                 tr(ir,i,j,il,1) = A
                 tr(ir,i,j,il,2) = B
                else if( A.ne.0.d0 ) then
                 tr(ir,igetiznm(i,j,il,m),1,0,1) = A
                 tr(ir,igetiznm(i,j,il,m),1,0,2) = B
                end if

               end do
               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 5 ) then

              if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
               do i = icmax+2, 1, -1
                  read(jsi,'(1p1000e11.3)')
     &            dmm0,
     &            ( tr(ir,i,l,il,ioe), l = 1, inmax+2 )

               end do
              else
               do i = icmax+2, 1, -1
                  read(jsi,'(1p1000e11.3)') dmm0,
     &            ( a0_tr(l), l = 1, inmax+2 )
                 do l = 1, inmax+2
                   if( a0_tr(l).ne.0.d0 ) then
                     tr(ir,igetiznm(i,l,il,m),1,0,ioe) = a0_tr(l)
                   end if
                 end do
               end do
              end if

*-----------------------------------------------------------------------
             end if

  590     continue

*-----------------------------------------------------------------------
         end if

  900 continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_yieldrz(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                        mz,mn,mm,nfr,nfz, ! frtati 2022/02/18 added mm
     &                        nr,nz,nn,rm,zm,nt,ikzz,iknn,
     &                        tr)

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
        common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
        common /tall83/ itnzn(itlmax), itndm(itlmax) ! frtati 2022/02/18

        common /tall21/ rtfac(itlmax)
        common /fact01/ facmax(itlmax)
        common /stat / istdev, irestart, ireschk
        character ctfln*100

        dimension   rm(nr+1)
        dimension   zm(nz+1)
        dimension   nt(nn)
        dimension   ikzz(maxpt,maxnt), iknn(maxpt,maxnt)
        dimension   tm(maxnt+maxpt,2)
        dimension   fm(nfz,nfr)
        dimension   tr(nr,nz,mz,mn,0:mm,2) ! frtati 2022/02/18 2 -> mm

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)
        dimension a0_tr(maxnt) ! frtati 2022/02/18

        character chin*200, chlw*200, chcm*200, ctem*2

*-----------------------------------------------------------------------

        character dc2*20
        character dc3*200

*-----------------------------------------------------------------------

        character elmnt(104)*3

        data elmnt /
     &    'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ',
     &    'Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ',
     &    'Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

        character chau*8
        character cha*1
        data cha /"'"/

        character rpa*1
        data rpa /'}'/
        character yen*1

*-----------------------------------------------------------------------

        iseek = 1
        inum  = 0

        if( nn .eq. 0 ) then
           nc = 1
        else
           nc = nn
        end if

*-----------------------------------------------------------------------
*        r axis
*-----------------------------------------------------------------------
         if( itaxs(m,iax) .eq. 6 ) then

           dc2  = '#r-lowerr-upper'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31

           do 190 jz = 1, nz
           do 190 ic = 1, nc

             if (iseek.eq.1) inum = inum + 1

             if( nn .gt. 0 ) then
                iz = nt(ic) / 1000
                ia = nt(ic) - iz * 1000
                if( ia .gt. 0 ) then
                   in = ia - iz
                   call chname(idum,ia,iz,chau)
                end if
             end if

             if( nn .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''iz  ='',i3,3x,''Z = all'')')
     &                      inum, jz
             else if( ia .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''iz  ='',i3,3x,''Z = '',i3,'' : '',a3)')
     &                      inum, jz, iz, elmnt(iz)
             else
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''iz  ='',i3,3x,a8)')
     &                      inum, jz, chau
             end if

*-----------------------------------------------------------------------

             if ( iseek .eq. 1 ) then
               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900
             end if

             if ( chin(1:i2).ne.dc3(1:i2) ) then
               iseek = 0
               goto 190
             else
               iseek = 1
             end if

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             do jr = 1, nr

               read(jsi,'(1p2e13.4,1pe13.4,0pf8.4)')
     &              dmm0, dmm1, sek, ser

               if( nn .eq. 0 ) then

                  do kn = 1, mn
                  do kz = 1, mz
                  do il = 0, mm ! frtati 2022/02/18 2 -> mm

                     lz  = ikzz( kz, kn )
                     ln  = iknn( kz, kn )

                  IF(il .eq. 0) then
                   tr(jr,jz,lz,ln,il,1) = sek
                   tr(jr,jz,lz,ln,il,2) = ser
                  ELSE
                   tr(jr,jz,lz,ln,il,1) = 0.d0
                   tr(jr,jz,lz,ln,il,2) = 0.d0
                  ENDIF

                  end do
                  end do
                  end do

               else if( ia .eq. 0 ) then

                  do kn = 1, mn
                  do il = 0, mm ! frtati 2022/02/18 2 -> mm

                     lz  = ikzz( iz, kn )
                     ln  = iknn( iz, kn )

                  IF(il .eq. 0) then
                   tr(jr,jz,lz,ln,il,1) = sek
                   tr(jr,jz,lz,ln,il,2) = ser
                  ELSE
                   tr(jr,jz,lz,ln,il,1) = 0.d0
                   tr(jr,jz,lz,ln,il,2) = 0.d0
                  ENDIF

                  end do
                  end do

               else

                  do il = 0, mm ! frtati 2022/02/18 2 -> mm
                     lz  = ikzz( iz, in )
                     ln  = iknn( iz, in )

                  IF(il .eq. 0) then
                   tr(jr,jz,lz,ln,il,1) = sek
                   tr(jr,jz,lz,ln,il,2) = ser
                  ELSE
                   tr(jr,jz,lz,ln,il,1) = 0.d0
                   tr(jr,jz,lz,ln,il,2) = 0.d0
                  ENDIF

                  end do

               end if

             end do

  190    continue

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 5 ) then

           dc2  = '#z-lowerz-upper'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31

           do 290 jr = 1, nr
           do 290 ic = 1, nc

             if (iseek.eq.1) inum = inum + 1

             if( nn .gt. 0 ) then
                iz = nt(ic) / 1000
                ia = nt(ic) - iz * 1000
                if( ia .gt. 0 ) then
                   in = ia - iz
                   call chname(idum,ia,iz,chau)
                end if
             end if

             if( nn .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''ir  ='',i3,3x,''Z = all'')')
     &                      inum, jr
             else if( ia .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''ir  ='',i3,3x,''Z = '',i3,'' : '',a3)')
     &                      inum, jr, iz, elmnt(iz)
             else
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''ir  ='',i3,3x,a8)')
     &                      inum, jr, chau
             end if

*-----------------------------------------------------------------------

             if ( iseek .eq. 1 ) then
               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900
             end if

             if ( chin(1:i2).ne.dc3(1:i2) ) then
               iseek = 0
               goto 290
             else
               iseek = 1
             end if

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             do jz = 1, nz

               read(jsi,'(1p2e13.4,1pe13.4,0pf8.4)')
     &              dmm0, dmm1, sek, ser

             if( nn .eq. 0 ) then

                do kn = 1, mn
                do kz = 1, mz
                do il = 0, mm ! frtati 2022/02/18 2 -> mm

                   lz  = ikzz( kz, kn )
                   ln  = iknn( kz, kn )

                  IF(il .eq. 0) then
                   tr(jr,jz,lz,ln,il,1) = sek
                   tr(jr,jz,lz,ln,il,2) = ser
                  ELSE
                   tr(jr,jz,lz,ln,il,1) = 0.d0
                   tr(jr,jz,lz,ln,il,2) = 0.d0
                  ENDIF

                end do
                end do
                end do

             else if( ia .eq. 0 ) then

                do kn = 1, mn
                do il = 0, mm ! frtati 2022/02/18 2 -> mm

                   lz  = ikzz( iz, kn )
                   ln  = iknn( iz, kn )

                  IF(il .eq. 0) then
                   tr(jr,jz,lz,ln,il,1) = sek
                   tr(jr,jz,lz,ln,il,2) = ser
                  ELSE
                   tr(jr,jz,lz,ln,il,1) = 0.d0
                   tr(jr,jz,lz,ln,il,2) = 0.d0
                  ENDIF

                end do
                end do

             else

                do il = 0, mm ! frtati 2022/02/18 2 -> mm
                   lz  = ikzz( iz, in )
                   ln  = iknn( iz, in )

                  IF(il .eq. 0) then
                   tr(jr,jz,lz,ln,il,1) = sek
                   tr(jr,jz,lz,ln,il,2) = ser
                  ELSE
                   tr(jr,jz,lz,ln,il,1) = 0.d0
                   tr(jr,jz,lz,ln,il,2) = 0.d0
                  ENDIF

                end do

             end if

             end do

  290    continue

*-----------------------------------------------------------------------
*        mass axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 1 ) then

           dc2  = '#massnumberr.err'
           ldc2 = 16

           facmx = 1.d0  ! kitamura23/03/31

           do 390 jr = 1, nr
           do 390 jz = 1, nz
           do 390 ic = 1, nc

             if (iseek.eq.1) inum = inum + 1

             if( nn .gt. 0 ) iz = nt(ic) / 1000

             do i = 1, maxnt + maxpt
                tm(i,1) = 0.0d+0
                tm(i,2) = 0.0d+0
             end do

             if( nn .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''ir ='',i3,3x,''iz ='',i3,3x,
     &          ''Z = all'')')
     &                      inum, jr, jz
             else
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''ir ='',i3,3x,''iz ='',i3,3x,''Z = '',i3,'' : '',a3)')
     &                      inum, jr, jz, iz, elmnt(iz)
             end if

*-----------------------------------------------------------------------

             if ( iseek .eq. 1 ) then
               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900
             end if

             if ( chin(1:i2).ne.dc3(1:i2) ) then
               iseek = 0
               goto 390
             else
               iseek = 1
             end if

             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# r ='
             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# z ='
             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
             read(chin(11:i2),'(i3,1x,i3)') im, jm
             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             if( im .eq. 1 ) then

                read(jsi,'(3x,f4.1,2x,1pe13.4,0pf8.4)') dmm0, dmm1, dmm2

             end if

             do i = im, jm

               read(jsi,'(3x,i4,2x,1pe13.4,0pf8.4)')
     &                       idmm0, tm(i,1), tm(i,2)

             end do

*-----------------------------------------------------------------------

             if( nn .eq. 0 .and. itfln(m) .eq. 1 ) then

                do ln = 1, mn
                do il = 0, mm ! frtati 2022/02/18 2 -> mm

                   lz = 1
                   im = ln

                  IF(il .eq. 0) then
                   tr(jr,jz,lz,ln,il,1) = tm(im,1)
                   tr(jr,jz,lz,ln,il,2) = tm(im,2)
                  ELSE
                   tr(jr,jz,lz,ln,il,1) = 0.d0
                   tr(jr,jz,lz,ln,il,2) = 0.d0
                  ENDIF

                end do
                end do

             else if( nn .eq. 0 .and. itfln(m) .gt. 1 ) then

                do kn = 1, mn
                do kz = 1, mz
                do il = 0, mm ! frtati 2022/02/18 2 -> mm

                   lz = ikzz( kz, kn )
                   ln = iknn( kz, kn )

                   im = iz + kn

                  IF(il .eq. 0) then
                   tr(jr,jz,lz,ln,il,1) = tm(im,1)
                   tr(jr,jz,lz,ln,il,2) = tm(im,2)
                  ELSE
                   tr(jr,jz,lz,ln,il,1) = 0.d0
                   tr(jr,jz,lz,ln,il,2) = 0.d0
                  ENDIF

                end do
                end do
                end do

             else

                do kn = 1, mn
                do il = 0, mm ! frtati 2022/02/18 2 -> mm

                   lz = ikzz( iz, kn )
                   ln = iknn( iz, kn )

                   im = iz + kn

                  IF(il .eq. 0) then
                   tr(jr,jz,lz,ln,il,1) = tm(im,1)
                   tr(jr,jz,lz,ln,il,2) = tm(im,2)
                  ELSE
                   tr(jr,jz,lz,ln,il,1) = 0.d0
                   tr(jr,jz,lz,ln,il,2) = 0.d0
                  ENDIF

                end do
                end do

             end if

  390 continue

*-----------------------------------------------------------------------
*        charge axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 7 ) then

           dc2  = '#chargenumberr.err'
           ldc2 = 18

           facmx = 1.d0  ! kitamura23/03/31

           do 490 jr = 1, nr
           do 490 jz = 1, nz

             if (iseek.eq.1) inum = inum + 1

             do i = 1, maxpt
                tm(i,1) = 0.0d+0
                tm(i,2) = 0.0d+0
             end do

             write(dc3,'(''#   no. ='',i3,3x,
     &       ''ir ='',i3,3x,''iz ='',i3,3x)')
     &                   inum, jr, jz

*-----------------------------------------------------------------------

             if ( iseek .eq. 1 ) then
               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900
             end if

             if ( chin(1:i2).ne.dc3(1:i2) ) then
               iseek = 0
               goto 490
             else
               iseek = 1
             end if

             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# r ='
             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# z ='
             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
             read(chin(11:i2),'(i3,1x,i3)') im, jm
             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             if( im .eq. 1 ) then

                read(jsi,'(3x,f4.1,2x,1pe13.4,0pf8.4)') dmm0, dmm1, dmm2

             end if

             do i = im, jm

               read(jsi,'(3x,i4,2x,1pe13.4,0pf8.4)')
     &                       idmm0, tm(i,1), tm(i,2)

             end do

*-----------------------------------------------------------------------

            do kz = 1, mz
            do kn = 1, mn
            do il = 0, mm ! frtati 2022/02/18 2 -> mm

                  iz  = ikzz( kz, kn )
                  in  = iknn( kz, kn )

                  IF(il .eq. 0) then
                   tr(jr,jz,iz,in,il,1) = tm(iz,1)
                   tr(jr,jz,iz,in,il,2) = tm(iz,2)
                  ELSE
                   tr(jr,jz,iz,in,il,1) = 0.d0
                   tr(jr,jz,iz,in,il,2) = 0.d0
                  ENDIF

            end do
            end do
            end do

  490     continue

*-----------------------------------------------------------------------
*        chart axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 8 ) then

           if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then
             dc2 = 'hc:y='
             ldc2 = 5
           else if( ittwo(m) .eq. 4 ) then
             dc2 = '#znmassnumberr.err'
             ldc2 = 18
           else if( ittwo(m) .eq. 5 ) then
             dc2 = 'z/n'
             ldc2 = 3
           end if

           if( itnzn(m).ne.0 ) tr(:,:,:,1,0,ioe) = 0.d0 ! frtati 2022/03/18

*-----------------------------------------------------------------------


           facmx = 1.d0  ! kitamura23/03/3

           do 590



*-----------------------------------------------------------------------

               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)
               if( jsn.eq.-1 ) goto 900 ! T.Sato 2023/04/13
            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='

               read(chin(16:17),'(A2)') ctem ! after reading all data, it goes to '# icmax...' line
               IF(ctem .ne. 'ir') goto 900
               read(chin(20:22),'(i3)') jr
               read(chin(31:33),'(i3)') jz
               read(chin(44:44),'(i1)') il



             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
             read(chin(17:i2),'(i3,1x,i3)') icmax, inmax
             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

*-----------------------------------------------------------------------
             if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               do i = icmax+2, 1, -1

                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                 read(jsi,'(1p10e11.3)')
     &             ( tr(jr,jz,i,l,il,ioe), l = 1, inmax+2 )
                else
                 read(jsi,'(1p10e11.3)')
     &             ( a0_tr(l), l = 1, inmax+2 )
                 do l = 1, inmax+2
                   if( a0_tr(l).ne.0.d0 ) then
                      tr(jr,jz,igetiznm(i,l,il,m),1,0,ioe) = a0_tr(l)
                   end if
                 end do
                end if

               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 4 ) then

               do i0 = 1, icmax
               do j0 = 1, inmax

                 call readl(jsn,jsi,dsin,idsi,ill,ilf,'',
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
                 if ( ierr .ne. 0  ) goto 590
                 if ( jpn  .eq. 3  ) goto 590
                 if ( i1   .eq. i3 ) goto 590

                 read(chin(1:i2),'(3i5,1pe13.4,0pf8.4)')
     &                i, j, ij, A, B

                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                 tr(jr,jz,i,j,il,1) = A
                 tr(jr,jz,i,j,il,2) = B
                else if( A.ne.0.d0 ) then
                 tr(jr,jz,igetiznm(i,j,il,m),1,0,1) = A
                 tr(jr,jz,igetiznm(i,j,il,m),1,0,2) = B
                end if

               end do
               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 5 ) then

              if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
               do i = icmax+2, 1, -1
                  read(jsi,'(1p1000e11.3)')
     &            dmm0,
     &            ( tr(jr,jz,i,l,il,ioe), l = 1, inmax+2 )

               end do
              else
               do i = icmax+2, 1, -1
                  read(jsi,'(1p1000e11.3)') dmm0,
     &            ( a0_tr(l), l = 1, inmax+2 )
                 do l = 1, inmax+2
                   if( a0_tr(l).ne.0.d0 ) then
                     tr(jr,jz,igetiznm(i,l,il,m),1,0,ioe) = a0_tr(l)
                   end if
                 end do
               end do
              end if

*-----------------------------------------------------------------------
             end if

  590     continue

*-----------------------------------------------------------------------
*        rz axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 12 ) then

           call sel_dc2(m,dc2,ldc2)
           if( ittwo(m) .eq. 4 ) then
             dc2 = '#rznumber'
             ldc2 = 9
           else if( ittwo(m) .eq. 5 ) then
             dc2 = 'r/z'
             ldc2 = 3
           end if

*-----------------------------------------------------------------------

           facmx = 1.d0  ! kitamura23/03/31

           do 690 ic = 1, nc

             if (iseek.eq.1) inum = inum + 1

             if( nn .gt. 0 ) then
               iz = nt(ic) / 1000
               ia = nt(ic) - iz * 1000
               if( ia .gt. 0 ) in = ia - iz
               call chname(idum,ia,iz,chau)
             end if

             if( nn .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,''Z = all'')')
     &                      inum
             else if( ia .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''Z = '',i3,'' : '',a3)')
     &                      inum, iz, elmnt(iz)
             else
                write(dc3,'(''#   no. ='',i3,3x,a8)')
     &                      inum, chau
             end if

*-----------------------------------------------------------------------

             if ( iseek .eq. 1 ) then
               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900
             end if

             if ( chin(1:i2).ne.dc3(1:i2) ) then
               iseek = 0
               goto 690
             else
               iseek = 1
             end if

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                      dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

*-----------------------------------------------------------------------
             if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               read(jsi,'(1p10e11.3)')
     &           ( ( fm(jz,jr), jz = 1, nz ), jr = nr, 1, -1 )

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 4 ) then

               do jz = 1, nz
               do jr = 1, nr

                 read(jsi,'(1p10e11.3)')
     &             dmm0, dmm1, fm(jz,jr)

               end do
               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 5 ) then

               do jr = nr, 1, -1

                  read(jsi,'(1p1000e11.3)')
     &              dmm0, ( fm(jz,jr), jz = 1, nz )

               end do

*-----------------------------------------------------------------------
             end if

             do 691 jz = 1, nz
             do 691 jr = 1, nr

                if( nn .eq. 0 ) then

                   do kn = 1, mn
                   do kz = 1, mz
                   do il = 0, mm ! frtati 2022/02/18 2 -> mm

                      lz  = ikzz( kz, kn )
                      ln  = iknn( kz, kn )

                     IF(il .eq. 0) then
                      tr(jr,jz,lz,ln,il,ioe) = fm(jz,jr)
                     ELSE
                      tr(jr,jz,lz,ln,il,ioe) = 0.d0
                     ENDIF

                   end do
                   end do
                   end do

                else if( ia .eq. 0 ) then

                   do kn = 1, mn
                   do il = 0, mm ! frtati 2022/02/18 2 -> mm

                      lz  = ikzz( iz, kn )
                      ln  = iknn( iz, kn )

                     IF(il .eq. 0) then
                      tr(jr,jz,lz,ln,il,ioe) = fm(jz,jr)
                     ELSE
                      tr(jr,jz,lz,ln,il,ioe) = 0.d0
                     ENDIF

                   end do
                   end do

                else
                  do il = 0, mm ! frtati 2022/02/18 2 -> mm

                      lz  = ikzz( iz, in )
                      ln  = iknn( iz, in )

                     IF(il .eq. 0) then
                      tr(jr,jz,lz,ln,il,ioe) = fm(jz,jr)
                     ELSE
                      tr(jr,jz,lz,ln,il,ioe) = 0.d0
                     ENDIF
                  end do

                end if

  691     continue
  690     continue

*-----------------------------------------------------------------------
         end if

  900 continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_yieldxyz(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                       mz,mn,mm,nf,nl,lt,! frtati 2022/02/18 added mm
     &                       nx,ny,nz,nn,xm,ym,zm,nt,ikzz,iknn,
     &                       tr)
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
        common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
        common /tall83/ itnzn(itlmax), itndm(itlmax) ! frtati 2022/02/18

        common /tall21/ rtfac(itlmax)
        common /fact01/ facmax(itlmax)
        common /stat / istdev, irestart, ireschk
        character ctfln*100

        dimension   lt(nl)
        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)
        dimension   nt(nn)
        dimension   ikzz(maxpt,maxnt), iknn(maxpt,maxnt)
        dimension   tm(maxpt+maxnt,2)
        dimension   fm(nf,nf)
        dimension   tr(nx,ny,nz,mz,mn,0:mm,2) ! frtati 2022/02/18 2 -> mm

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)
        dimension a0_tr(maxnt) ! frtati 2022/02/18

        character chin*200, chlw*200, chcm*200, ctem*2

*-----------------------------------------------------------------------

        character dc2*20
        character dc3*200

*-----------------------------------------------------------------------

        character elmnt(104)*3

        data elmnt /
     &    'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ',
     &    'Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ',
     &    'Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

        character chau*8
        character cha*1
        data cha /"'"/

        character rpa*1
        data rpa /'}'/
        character yen*1

*-----------------------------------------------------------------------
        common /redufmt/ iredufmt(itlmax) !FURUTA20200615

        iseek = 1
        inum  = 0

        if( nn .eq. 0 ) then
           nc = 1
        else
           nc = nn
        end if

*-----------------------------------------------------------------------
*     dchain axis
*-----------------------------------------------------------------------
         if( itaxs(m,iax) .eq. 13 ) then

          if(iredufmt(m).eq.0)then !FURUTA20200615

              !! seek to next isotope production
              jpn = 0
  110         call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$#',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
              if ( ierr  .ne. 0  ) goto 900
              if ( jpn   .eq. 3  ) goto 900

              if ( iskip .ne. 0  ) goto 110
              if ( i2    .lt. 27 ) goto 110

              if (chin(9:27) .eq. ' isotope production') then
               il = 0
              elseif (chin(9:27) .eq. ' 1st metastable iso') then
               il = 1
              elseif (chin(9:27) .eq. ' 2nd metastable iso') then
               il = 2
              else
               goto 110
              endif

*-----------------------------------------------------------------------

              read(chin(1:5),'(i5)') iz
              IF(il .eq. 0) then
               read(chin(39:45),'(i3,1x,i3)') n3, n4
              ELSE
               read(chin(54:60),'(i3,1x,i3)') n3, n4
              ENDIF

              call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

              do jx=1,nx
              do jy=1,ny
              do jz=1,nz

                 read(jsi,'(3i7,1p12e11.3)')
     &              jx0,jy0,jz0,(tr(jx,jy,jz,iz,i,il,ioe),i=n3,n4)

              end do
              enddo
              enddo

              call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$#',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

           goto 110

*-----------------------------------------------------------------------
          else

           call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        '#numnucleusidyieldr.err',13,ierr)

           if ( ierr  .ne. 0 ) goto 900
           if ( jpn   .eq. 3 ) goto 900

           do
            read(jsi,*)ir,nucleusID,yield,e_yield
            if(ir.eq.0)exit
            jz=(ir-1)/(nx*ny)+1
            jy=(ir-(jz-1)*nx*ny-1)/nx+1
            jx=mod(ir-1,nx)+1
            iz=nucleusID/10000
            ia=(nucleusID-iz*10000)/10
            in=ia-iz
            il=mod(nucleusID,10)
            if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
             tr(jx,jy,jz,iz,in,il,1)=yield
             tr(jx,jy,jz,iz,in,il,2)=e_yield
            else
             tr(jx,jy,jz,igetiznm(iz,in,il,m),1,0,1)=yield
             tr(jx,jy,jz,igetiznm(iz,in,il,m),1,0,2)=e_yield
            end if
           enddo

          endif
*-----------------------------------------------------------------------
*        x axis
*-----------------------------------------------------------------------
         elseif( itaxs(m,iax) .eq. 3 ) then

           dc2  = '#x-lowerx-upper'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31

           do 190 jy = 1, ny
           do 190 jz = 1, nz
           do 190 ic = 1, nc

             if (iseek.eq.1) inum = inum + 1

             if( nn .gt. 0 ) then
                iz = nt(ic) / 1000
                ia = nt(ic) - iz * 1000
                if( ia .gt. 0 ) then
                   in = ia - iz
                   call chname(idum,ia,iz,chau)
                end if
             end if

             if( nn .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''iy  ='',i3,3x,
     &          ''iz  ='',i3,3x,''Z = all'')')
     &                      inum, jy, jz
             else if( ia .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''iy  ='',i3,3x,
     &          ''iz  ='',i3,3x,''Z = '',i3,'' : '',a3)')
     &                      inum, jy, jz, iz, elmnt(iz)
             else
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''iy  ='',i3,3x,
     &          ''iz  ='',i3,3x,a8)')
     &                      inum, jy, jz, chau
             end if

*-----------------------------------------------------------------------

             if ( iseek .eq. 1 ) then
               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900
             end if

             if ( chin(1:i2).ne.dc3(1:i2) ) then
               iseek = 0
               goto 190
             else
               iseek = 1
             end if

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             do jx = 1, nx

               read(jsi,'(1p2e13.4,1pe13.4,0pf8.4)')
     &              dmm0, dmm1, sek, ser

               if( nn .eq. 0 ) then

                  do kn = 1, mn
                  do kz = 1, mz
                  do il = 0, mm ! frtati 2022/02/18 2 -> mm

                     lz  = ikzz( kz, kn )
                     ln  = iknn( kz, kn )

                    IF(il .eq. 0) then
                     tr(jx,jy,jz,lz,ln,il,1) = sek
                     tr(jx,jy,jz,lz,ln,il,2) = ser
                    ELSE
                     tr(jx,jy,jz,lz,ln,il,1) = 0.d0
                     tr(jx,jy,jz,lz,ln,il,2) = 0.d0
                    ENDIF

                  end do
                  end do
                  end do

               else if( ia .eq. 0 ) then

                  do kn = 1, mn
                  do il = 0, mm ! frtati 2022/02/18 2 -> mm

                     lz  = ikzz( iz, kn )
                     ln  = iknn( iz, kn )

                    IF(il .eq. 0) then
                     tr(jx,jy,jz,lz,ln,il,1) = sek
                     tr(jx,jy,jz,lz,ln,il,2) = ser
                    ELSE
                     tr(jx,jy,jz,lz,ln,il,1) = 0.d0
                     tr(jx,jy,jz,lz,ln,il,2) = 0.d0
                    ENDIF

                  end do
                  end do

               else

                  do il = 0, mm ! frtati 2022/02/18 2 -> mm
                     lz  = ikzz( iz, in )
                     ln  = iknn( iz, in )

                    IF(il .eq. 0) then
                     tr(jx,jy,jz,lz,ln,il,1) = sek
                     tr(jx,jy,jz,lz,ln,il,2) = ser
                    ELSE
                     tr(jx,jy,jz,lz,ln,il,1) = 0.d0
                     tr(jx,jy,jz,lz,ln,il,2) = 0.d0
                    ENDIF
                  end do

               end if

             end do

  190      continue

*-----------------------------------------------------------------------
*        y axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 4 ) then

           dc2  = '#y-lowery-upper'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31

           do 290 jx = 1, nx
           do 290 jz = 1, nz
           do 290 ic = 1, nc

             if (iseek.eq.1) inum = inum + 1

             if( nn .gt. 0 ) then
                iz = nt(ic) / 1000
                ia = nt(ic) - iz * 1000
                if( ia .gt. 0 ) then
                   in = ia - iz
                   call chname(idum,ia,iz,chau)
                end if
             end if

             if( nn .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''ix  ='',i3,3x,
     &          ''iz  ='',i3,3x,''Z = all'')')
     &                      inum, jx, jz
             else if( ia .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''ix  ='',i3,3x,
     &          ''iz  ='',i3,3x,''Z = '',i3,'' : '',a3)')
     &                      inum, jx, jz, iz, elmnt(iz)
             else
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''ix  ='',i3,3x,
     &          ''iz  ='',i3,3x,a8)')
     &                      inum, jx, jz, chau
             end if

*-----------------------------------------------------------------------

             if ( iseek .eq. 1 ) then
               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900
             end if

             if ( chin(1:i2).ne.dc3(1:i2) ) then
               iseek = 0
               goto 290
             else
               iseek = 1
             end if

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             do jy = 1, ny

               read(jsi,'(1p2e13.4,1pe13.4,0pf8.4)')
     &              dmm0, dmm1, sek, ser

               if( nn .eq. 0 ) then

                  do kn = 1, mn
                  do kz = 1, mz
                  do il = 0, mm ! frtati 2022/02/18 2 -> mm

                     lz  = ikzz( kz, kn )
                     ln  = iknn( kz, kn )

                    IF(il .eq. 0) then
                     tr(jx,jy,jz,lz,ln,il,1) = sek
                     tr(jx,jy,jz,lz,ln,il,2) = ser
                    ELSE
                     tr(jx,jy,jz,lz,ln,il,1) = 0.d0
                     tr(jx,jy,jz,lz,ln,il,2) = 0.d0
                    ENDIF

                  end do
                  end do
                  end do

               else if( ia .eq. 0 ) then

                  do kn = 1, mn
                  do il = 0, mm ! frtati 2022/02/18 2 -> mm

                     lz  = ikzz( iz, kn )
                     ln  = iknn( iz, kn )

                    IF(il .eq. 0) then
                     tr(jx,jy,jz,lz,ln,il,1) = sek
                     tr(jx,jy,jz,lz,ln,il,2) = ser
                    ELSE
                     tr(jx,jy,jz,lz,ln,il,1) = 0.d0
                     tr(jx,jy,jz,lz,ln,il,2) = 0.d0
                    ENDIF

                  end do
                  end do

               else

                  do il = 0, mm ! frtati 2022/02/18 2 -> mm
                     lz  = ikzz( iz, in )
                     ln  = iknn( iz, in )

                    IF(il .eq. 0) then
                     tr(jx,jy,jz,lz,ln,il,1) = sek
                     tr(jx,jy,jz,lz,ln,il,2) = ser
                    ELSE
                     tr(jx,jy,jz,lz,ln,il,1) = 0.d0
                     tr(jx,jy,jz,lz,ln,il,2) = 0.d0
                    ENDIF
                  end do

               end if

             end do

  290      continue

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 5 ) then

           dc2  = '#z-lowerz-upper'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31

           do 390 jx = 1, nx
           do 390 jy = 1, ny
           do 390 ic = 1, nc

             if (iseek.eq.1) inum = inum + 1

             if( nn .gt. 0 ) then
                iz = nt(ic) / 1000
                ia = nt(ic) - iz * 1000
                if( ia .gt. 0 ) then
                   in = ia - iz
                   call chname(idum,ia,iz,chau)
                end if
             end if

             if( nn .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''ix  ='',i3,3x,
     &          ''iy  ='',i3,3x,''Z = all'')')
     &                      inum, jx, jy
             else if( ia .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''ix  ='',i3,3x,
     &          ''iy  ='',i3,3x,''Z = '',i3,'' : '',a3)')
     &                      inum, jx, jy, iz, elmnt(iz)
             else
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''ix  ='',i3,3x,
     &          ''iy  ='',i3,3x,a8)')
     &                      inum, jx, jy, chau
             end if

*-----------------------------------------------------------------------

             if ( iseek .eq. 1 ) then
               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900
             end if

             if ( chin(1:i2).ne.dc3(1:i2) ) then
               iseek = 0
               goto 390
             else
               iseek = 1
             end if

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             do jz = 1, nz

               read(jsi,'(1p2e13.4,1pe13.4,0pf8.4)')
     &              dmm0, dmm1, sek, ser

               if( nn .eq. 0 ) then

                  do kn = 1, mn
                  do kz = 1, mz
                  do il = 0, mm ! frtati 2022/02/18 2 -> mm

                     lz  = ikzz( kz, kn )
                     ln  = iknn( kz, kn )

                    IF(il .eq. 0) then
                     tr(jx,jy,jz,lz,ln,il,1) = sek
                     tr(jx,jy,jz,lz,ln,il,2) = ser
                    ELSE
                     tr(jx,jy,jz,lz,ln,il,1) = 0.d0
                     tr(jx,jy,jz,lz,ln,il,2) = 0.d0
                    ENDIF

                  end do
                  end do
                  end do

               else if( ia .eq. 0 ) then

                  do kn = 1, mn
                  do il = 0, mm ! frtati 2022/02/18 2 -> mm

                     lz  = ikzz( iz, kn )
                     ln  = iknn( iz, kn )

                    IF(il .eq. 0) then
                     tr(jx,jy,jz,lz,ln,il,1) = sek
                     tr(jx,jy,jz,lz,ln,il,2) = ser
                    ELSE
                     tr(jx,jy,jz,lz,ln,il,1) = 0.d0
                     tr(jx,jy,jz,lz,ln,il,2) = 0.d0
                    ENDIF

                  end do
                  end do

               else

                  do il = 0, mm ! frtati 2022/02/18 2 -> mm
                     lz  = ikzz( iz, in )
                     ln  = iknn( iz, in )

                    IF(il .eq. 0) then
                     tr(jx,jy,jz,lz,ln,il,1) = sek
                     tr(jx,jy,jz,lz,ln,il,2) = ser
                    ELSE
                     tr(jx,jy,jz,lz,ln,il,1) = 0.d0
                     tr(jx,jy,jz,lz,ln,il,2) = 0.d0
                    ENDIF
                  end do

               end if

             end do

  390      continue

*-----------------------------------------------------------------------
*        mass axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 1 ) then

           dc2  = '#massnumberr.err'
           ldc2 = 16

           facmx = 1.d0  ! kitamura23/03/31

           do 490 jx = 1, nx
           do 490 jy = 1, ny
           do 490 jz = 1, nz
           do 490 ic = 1, nc

             if (iseek.eq.1) inum = inum + 1

             if( nn .gt. 0 ) iz = nt(ic) / 1000

             do i = 1, maxnt + maxpt
                tm(i,1) = 0.0d+0
                tm(i,2) = 0.0d+0
             end do

             if( nn .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''ix ='',i3,3x,''iy ='',i3,3x,''iz ='',i3,3x,
     &          ''Z = all'')')
     &                      inum, jx, jy, jz
             else
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''ix ='',i3,3x,''iy ='',i3,3x,''iz ='',i3,3x,
     &          ''Z = '',i3,'' : '',a3)')
     &                      inum, jx, jy, jz, iz, elmnt(iz)
             end if

*-----------------------------------------------------------------------

             if ( iseek .eq. 1 ) then
               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900
             end if

             if ( chin(1:i2).ne.dc3(1:i2) ) then
               iseek = 0
               goto 490
             else
               iseek = 1
             end if

             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# x ='
             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# y ='
             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# z ='

             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
             read(chin(11:i2),'(i3,1x,i3)') im, jm

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             if( im .eq. 1 ) then

                read(jsi,'(3x,f4.1,2x,1pe13.4,0pf8.4)') dmm0, dmm1, dmm2

             end if

             do i = im, jm

               read(jsi,'(3x,i4,2x,1pe13.4,0pf8.4)')
     &                       idmm0, tm(i,1), tm(i,2)

             end do

*-----------------------------------------------------------------------

            if( nn .eq. 0 .and. itfln(m) .eq. 1 ) then

               do ln = 1, mn
               do il = 0, mm ! frtati 2022/02/18 2 -> mm

                  lz = 1
                  im = ln

                  IF(il .eq. 0) then
                   tr(jx,jy,jz,lz,ln,il,1) = tm(im,1)
                   tr(jx,jy,jz,lz,ln,il,2) = tm(im,2)
                  ELSE
                   tr(jx,jy,jz,lz,ln,il,1) = 0.d0
                   tr(jx,jy,jz,lz,ln,il,2) = 0.d0
                  ENDIF

               end do
               end do

            else if( nn .eq. 0 .and. itfln(m) .gt. 1 ) then

               do kn = 1, mn
               do kz = 1, mz
               do il = 0, mm ! frtati 2022/02/18 2 -> mm

                  lz = ikzz( kz, kn )
                  ln = iknn( kz, kn )

                  im = kz + kn

                  IF(il .eq. 0) then
                   tr(jx,jy,jz,lz,ln,il,1) = tm(im,1)
                   tr(jx,jy,jz,lz,ln,il,2) = tm(im,2)
                  ELSE
                   tr(jx,jy,jz,lz,ln,il,1) = 0.d0
                   tr(jx,jy,jz,lz,ln,il,2) = 0.d0
                  ENDIF

               end do
               end do
               end do

            else

               do kn = 1, mn
               do il = 0, mm ! frtati 2022/02/18 2 -> mm

                  lz = ikzz( iz, kn )
                  ln = iknn( iz, kn )

                  im = iz + kn

                  IF(il .eq. 0) then
                   tr(jx,jy,jz,lz,ln,il,1) = tm(im,1)
                   tr(jx,jy,jz,lz,ln,il,2) = tm(im,2)
                  ELSE
                   tr(jx,jy,jz,lz,ln,il,1) = 0.d0
                   tr(jx,jy,jz,lz,ln,il,2) = 0.d0
                  ENDIF

               end do
               end do

            end if

  490     continue

*-----------------------------------------------------------------------
*        charge axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 7 ) then

           dc2  = '#chargenumberr.err'
           ldc2 = 18

           iseek = 1

           facmx = 1.d0  ! kitamura23/03/31

           do 590 jx = 1, nx
           do 590 jy = 1, ny
           do 590 jz = 1, nz

             if (iseek.eq.1) inum = inum + 1

             do i = 1, maxpt
                tm(i,1) = 0.0d+0
                tm(i,2) = 0.0d+0
             end do

             write(dc3,'(''#   no. ='',i3,3x,
     &       ''ix ='',i3,3x,''iy ='',i3,3x,''iz ='',i3,3x)')
     &                   inum, jx, jy, jz

*-----------------------------------------------------------------------

             if ( iseek .eq. 1 ) then
               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900
             end if

             if ( chin(1:i2).ne.dc3(1:i2) ) then
               iseek = 0
               goto 590
             else
               iseek = 1
             end if

             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# x ='
             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# y ='
             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# z ='
             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
             read(chin(11:i2),'(i3,1x,i3)') im, jm
             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             if( im .eq. 1 ) then

                read(jsi,'(3x,f4.1,2x,1pe13.4,0pf8.4)') dmm0, dmm1, dmm2

             end if

             do i = im, jm

               read(jsi,'(3x,i4,2x,1pe13.4,0pf8.4)')
     &                       idmm0, tm(i,1), tm(i,2)

             end do

*-----------------------------------------------------------------------

            do kz = 1, mz
            do kn = 1, mn
            do il = 0, mm ! frtati 2022/02/18 2 -> mm

                  iz  = ikzz( kz, kn )
                  in  = iknn( kz, kn )

                  IF(il .eq. 0) then
                   tr(jx,jy,jz,iz,in,il,1) = tm(iz,1)
                   tr(jx,jy,jz,iz,in,il,2) = tm(iz,2)
                  ELSE
                   tr(jx,jy,jz,iz,in,il,1) = 0.d0
                   tr(jx,jy,jz,iz,in,il,2) = 0.d0
                  ENDIF

            end do
            end do
            end do

  590     continue

*-----------------------------------------------------------------------
*        chart axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 8 ) then

           if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then
             dc2 = 'hc:y='
             ldc2 = 5
           else if( ittwo(m) .eq. 4 ) then
             dc2 = '#znmassnumberr.err'
             ldc2 = 18
           else if( ittwo(m) .eq. 5 ) then
             dc2 = 'z/n'
             ldc2 = 3
           end if

           if( itnzn(m).ne.0 ) tr(:,:,:,:,1,0,ioe) = 0.d0 ! frtati 2022/03/18

*-----------------------------------------------------------------------


           facmx = 1.d0  ! kitamura23/03/31

           do 690



*-----------------------------------------------------------------------

               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)
             if( jsn.eq.-1 ) goto 900 ! T.Sato 2023/04/13, no more data
            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='


               read(chin(16:17),'(A2)') ctem ! after reading all data, it goes to '# icmax...' line
               IF(ctem .ne. 'ix') goto 900
               read(chin(20:22),'(i3)') jx
               read(chin(30:32),'(i3)') jy
               read(chin(41:43),'(i3)') jz
               read(chin(54:54),'(i1)') il

             call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
             read(chin(17:i2),'(i3,1x,i3)') icmax, inmax
             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

*-----------------------------------------------------------------------
             if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               do i = icmax+2, 1, -1

                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                 read(jsi,'(1p10e11.3)')
     &             ( tr(jx,jy,jz,i,l,il,ioe), l = 1, inmax+2 )
                else
                 read(jsi,'(1p10e11.3)')
     &             ( a0_tr(l), l = 1, inmax+2 )
                 do l = 1, inmax+2
                   if( a0_tr(l).ne.0.d0 ) then
                     tr(jx,jy,jz,igetiznm(i,l,il,m),1,0,ioe) = a0_tr(l)
                   end if
                 end do
                end if

               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 4 ) then

               do i0 = 1, icmax
               do j0 = 1, inmax

                 call readl(jsn,jsi,dsin,idsi,ill,ilf,'',
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
                 if ( ierr .ne. 0  ) goto 690
                 if ( jpn  .eq. 3  ) goto 690
                 if ( i1   .eq. i3 ) goto 690

                 read(chin(1:i2),'(3i5,1pe13.4,0pf8.4)')
     &                i, j, ij, A, B

                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                 tr(jx,jy,jz,i,j,il,1) = A
                 tr(jx,jy,jz,i,j,il,2) = B
                else if( A.ne.0.d0 ) then
                 tr(jx,jy,jz,igetiznm(i,j,il,m),1,0,1) = A
                 tr(jx,jy,jz,igetiznm(i,j,il,m),1,0,2) = B
                end if

               end do
               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 5 ) then

              if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
               do i = icmax+2, 1, -1
                  read(jsi,'(1p1000e11.3)')
     &            dmm0,
     &            ( tr(jx,jy,jz,i,l,il,ioe), l = 1, inmax+2 )

               end do
              else
               do i = icmax+2, 1, -1
                  read(jsi,'(1p1000e11.3)') dmm0,
     &            ( a0_tr(l), l = 1, inmax+2 )
                 do l = 1, inmax+2
                   if( a0_tr(l).ne.0.d0 ) then
                     tr(jx,jy,jz,igetiznm(i,l,il,m),1,0,ioe) = a0_tr(l)
                   end if
                 end do
               end do
              end if

*-----------------------------------------------------------------------
             end if

  690     continue

*-----------------------------------------------------------------------
*        xy axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 9 ) then

           call sel_dc2(m,dc2,ldc2)
           if( ittwo(m) .eq. 4 ) then
             dc2 = '#xynumber'
             ldc2 = 9
           else if( ittwo(m) .eq. 5 ) then
             dc2 = 'y/x'
             ldc2 = 3
           end if

*-----------------------------------------------------------------------

           facmx = 1.d0  ! kitamura23/03/31

           do 790 jz = 1, nz
           do 790 ic = 1, nc

             if (iseek.eq.1) inum = inum + 1

             if( nn .gt. 0 ) then
               iz = nt(ic) / 1000
               ia = nt(ic) - iz * 1000
               if( ia .gt. 0 ) then
                 in = ia - iz
                 call chname(idum,ia,iz,chau)
               end if
             end if

             if( nn .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''iz  ='',i3,3x,''Z = all'')')
     &                      inum, jz
             else if( ia .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''iz  ='',i3,3x,''Z = '',i3,'' : '',a3)')
     &                      inum, jz, iz, elmnt(iz)
             else
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''iz  ='',i3,3x,a8)')
     &                      inum, jz, chau
             end if

*-----------------------------------------------------------------------

             if ( iseek .eq. 1 ) then
               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900
             end if

             if ( chin(1:i2).ne.dc3(1:i2) ) then
               iseek = 0
               goto 790
             else
               iseek = 1
             end if

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

*-----------------------------------------------------------------------
             if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               read(jsi,'(1p10e11.3)')
     &           ( ( fm(jx,jy), jx = 1, nx ), jy = ny, 1, -1 )

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 4 ) then

               do jy = 1, ny
               do jx = 1, nx

                 read(jsi,'(1p10e11.3)')
     &             dmm0, dmm1, fm(jx,jy)

               end do
               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 5 ) then

               do jy = ny, 1, -1

                  read(jsi,'(1p1000e11.3)')
     &              dmm0, ( fm(jx,jy), jx = 1, nx )

               end do

*-----------------------------------------------------------------------
             end if

             do 791 jx = 1, nx
             do 791 jy = 1, ny

                if( nn .eq. 0 ) then

                   do kn = 1, mn
                   do kz = 1, mz
                   do il = 0, mm ! frtati 2022/02/18 2 -> mm

                      lz  = ikzz( kz, kn )
                      ln  = iknn( kz, kn )

                      IF(il .eq. 0) then
                       tr(jx,jy,jz,lz,ln,il,ioe) = fm(jx,jy)
                      ELSE
                       tr(jx,jy,jz,lz,ln,il,ioe) = 0.d0
                      ENDIF

                   end do
                   end do
                   end do

                else if( ia .eq. 0 ) then

                   do kn = 1, mn
                   do il = 0, mm ! frtati 2022/02/18 2 -> mm

                      lz  = ikzz( iz, kn )
                      ln  = iknn( iz, kn )

                      IF(il .eq. 0) then
                       tr(jx,jy,jz,lz,ln,il,ioe) = fm(jx,jy)
                      ELSE
                       tr(jx,jy,jz,lz,ln,il,ioe) = 0.d0
                      ENDIF

                   end do
                   end do

                else

                  do il = 0, mm ! frtati 2022/02/18 2 -> mm
                      lz  = ikzz( iz, in )
                      ln  = iknn( iz, in )

                      IF(il .eq. 0) then
                       tr(jx,jy,jz,lz,ln,il,ioe) = fm(jx,jy)
                      ELSE
                       tr(jx,jy,jz,lz,ln,il,ioe) = 0.d0
                      ENDIF
                  end do

                end if

  791     continue
  790     continue

*-----------------------------------------------------------------------
*        yz axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 10 ) then

           call sel_dc2(m,dc2,ldc2)
           if( ittwo(m) .eq. 4 ) then
             dc2 = '#yznumber'
             ldc2 = 9
           else if( ittwo(m) .eq. 5 ) then
             dc2 = 'y/z'
             ldc2 = 3
           end if

*-----------------------------------------------------------------------

           facmx = 1.d0  ! kitamura23/03/31

           do 890 jx = 1, nx
           do 890 ic = 1, nc

             if (iseek.eq.1) inum = inum + 1

             if( nn .gt. 0 ) then
               iz = nt(ic) / 1000
               ia = nt(ic) - iz * 1000
               if( ia .gt. 0 ) then
                 in = ia - iz
                 call chname(idum,ia,iz,chau)
               end if
             end if

             if( nn .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''ix  ='',i3,3x,''Z = all'')')
     &                      inum, jx
             else if( ia .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''ix  ='',i3,3x,''Z = '',i3,'' : '',a3)')
     &                      inum, jx, iz, elmnt(iz)
             else
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''ix  ='',i3,3x,a8)')
     &                      inum, jx, chau
             end if

*-----------------------------------------------------------------------

             if ( iseek .eq. 1 ) then
               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900
             end if

             if ( chin(1:i2).ne.dc3(1:i2) ) then
               iseek = 0
               goto 890
             else
               iseek = 1
             end if

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

*-----------------------------------------------------------------------
             if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               read(jsi,'(1p10e11.3)')
     &           ( ( fm(jz,jy), jz = 1, nz ), jy = ny, 1, -1 )

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 4 ) then

               do jz = 1, nz
               do jy = 1, ny

                 read(jsi,'(1p10e11.3)')
     &             dmm0, dmm1, fm(jz,jy)

               end do
               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 5 ) then

               do jy = ny, 1, -1

                  read(jsi,'(1p1000e11.3)')
     &              dmm0, ( fm(jz,jy), jz = 1, nz )

               end do

*-----------------------------------------------------------------------
             end if

             do 891 jz = 1, nz
             do 891 jy = 1, ny

                if( nn .eq. 0 ) then

                   do kn = 1, mn
                   do kz = 1, mz
                   do il = 0, mm ! frtati 2022/02/18 2 -> mm

                      lz  = ikzz( kz, kn )
                      ln  = iknn( kz, kn )

                      IF(il .eq. 0) then
                       tr(jx,jy,jz,lz,ln,il,ioe) = fm(jz,jy)
                      ELSE
                       tr(jx,jy,jz,lz,ln,il,ioe) = 0.d0
                      ENDIF

                   end do
                   end do
                   end do

                else if( ia .eq. 0 ) then

                   do kn = 1, mn
                   do il = 0, mm ! frtati 2022/02/18 2 -> mm

                      lz  = ikzz( iz, kn )
                      ln  = iknn( iz, kn )

                      IF(il .eq. 0) then
                       tr(jx,jy,jz,lz,ln,il,ioe) = fm(jz,jy)
                      ELSE
                       tr(jx,jy,jz,lz,ln,il,ioe) = 0.d0
                      ENDIF

                   end do
                   end do

                else

                  do il = 0, mm ! frtati 2022/02/18 2 -> mm
                      lz  = ikzz( iz, in )
                      ln  = iknn( iz, in )

                      IF(il .eq. 0) then
                       tr(jx,jy,jz,lz,ln,il,ioe) = fm(jz,jy)
                      ELSE
                       tr(jx,jy,jz,lz,ln,il,ioe) = 0.d0
                      ENDIF

                  end do

                end if

  891     continue
  890     continue

*-----------------------------------------------------------------------
*        xz axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 11 ) then

           call sel_dc2(m,dc2,ldc2)
           if( ittwo(m) .eq. 4 ) then
             dc2 = '#xznumber'
             ldc2 = 9
           else if( ittwo(m) .eq. 5 ) then
             dc2 = 'x/z'
             ldc2 = 3
           end if

*-----------------------------------------------------------------------

           facmx = 1.d0  ! kitamura23/03/31

           do 990 jy = 1, ny
           do 990 ic = 1, nc

             if (iseek.eq.1) inum = inum + 1

             if( nn .gt. 0 ) then
               iz = nt(ic) / 1000
               ia = nt(ic) - iz * 1000
               if( ia .gt. 0 ) then
                 in = ia - iz
                 call chname(idum,ia,iz,chau)
               end if
             end if

             if( nn .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''iy  ='',i3,3x,''Z = all'')')
     &                      inum, jy
             else if( ia .eq. 0 ) then
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''iy  ='',i3,3x,''Z = '',i3,'' : '',a3)')
     &                      inum, jy, iz, elmnt(iz)
             else
                write(dc3,'(''#   no. ='',i3,3x,
     &          ''iy  ='',i3,3x,a8)')
     &                      inum, jy, chau
             end if

*-----------------------------------------------------------------------

             if ( iseek .eq. 1 ) then
               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr) ! '# no ='
               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900
             end if

             if ( chin(1:i2).ne.dc3(1:i2) ) then
               iseek = 0
               goto 990
             else
               iseek = 1
             end if

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                        dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

*-----------------------------------------------------------------------
             if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               read(jsi,'(1p10e11.3)')
     &           ( ( fm(jz,jx), jz = 1, nz ), jx = nx, 1, -1 )

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 4 ) then

               do jz = 1, nz
               do jx = 1, nx

                 read(jsi,'(1p10e11.3)')
     &             dmm0, dmm1, fm(jz,jx)

               end do
               end do

*-----------------------------------------------------------------------
             else if( ittwo(m) .eq. 5 ) then

               do jx = nx, 1, -1

                  read(jsi,'(1p1000e11.3)')
     &              dmm0, ( fm(jz,jx), jz = 1, nz )

               end do

*-----------------------------------------------------------------------
             end if

             do 991 jz = 1, nz
             do 991 jx = 1, nx

                if( nn .eq. 0 ) then

                   do kn = 1, mn
                   do kz = 1, mz
                   do il = 0, mm ! frtati 2022/02/18 2 -> mm

                      lz  = ikzz( kz, kn )
                      ln  = iknn( kz, kn )

                      IF(il .eq. 0) then
                       tr(jx,jy,jz,lz,ln,il,ioe) = fm(jz,jx)
                      ELSE
                       tr(jx,jy,jz,lz,ln,il,ioe) = 0.d0
                      ENDIF

                   end do
                   end do
                   end do

                else if( ia .eq. 0 ) then

                   do kn = 1, mn
                   do il = 0, mm ! frtati 2022/02/18 2 -> mm

                      lz  = ikzz( iz, kn )
                      ln  = iknn( iz, kn )

                      IF(il .eq. 0) then
                       tr(jx,jy,jz,lz,ln,il,ioe) = fm(jz,jx)
                      ELSE
                       tr(jx,jy,jz,lz,ln,il,ioe) = 0.d0
                      ENDIF

                   end do
                   end do

                else

                  do il = 0, mm ! frtati 2022/02/18 2 -> mm
                      lz  = ikzz( iz, in )
                      ln  = iknn( iz, in )

                      IF(il .eq. 0) then
                       tr(jx,jy,jz,lz,ln,il,ioe) = fm(jz,jx)
                      ELSE
                       tr(jx,jy,jz,lz,ln,il,ioe) = 0.d0
                      ENDIF

                  end do

                end if

  991     continue
  990     continue

*-----------------------------------------------------------------------
         end if

  900 continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_yieldreg(m,
     &                         mz,mn,mm, ! frtati 2022/02/18 added mm
     &                         nr,mr,nn,kr,nt,ikzz,iknn,
     &                         tr,nvl,ivl,rvl,
     &                         nx,ny,nz,xm,ym,zm)
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

        dimension   kr(mr)
        dimension   nt(nn)
        dimension   ikzz(maxpt,maxnt), iknn(maxpt,maxnt)
        dimension   tm(maxnt+maxpt,2)
        dimension   vl(nr)
        dimension   lr(nr)
        dimension   tr(nr,mz,mn,0:mm,2) ! frtati 2022/05/02 2 -> mm
        dimension   ivl(nvl)
        dimension   rvl(nvl)

        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)

*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

        call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)

        if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

        do 100 ir = 1, nr

           if( itunt(m) .eq. 1 ) then

              cc = 1.d0 / abs(rtfac(m)/facmax(m))

           else if( itunt(m) .eq. 2 ) then

              cc = vl(ir) / abs(rtfac(m)/facmax(m))

           end if

           do 100 iz = 1, mz
           do 100 in = 1, mn
           do 100 il = 0, mm ! frtati 2022/02/18 2 -> mm

                 call invert_stdev(m,A,B,
     &                           tr(ir,iz,in,il,1),
     &                           tr(ir,iz,in,il,2),
     &                           cc)

                 tr(ir,iz,in,il,1) = A
                 tr(ir,iz,in,il,2) = B

  100   continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_yieldtet(m,
     &                         mz,mn,mm,nr,mr,nn,nt,ikzz,iknn, ! frtati 2022/02/18 added mm
     &                         tr,
     &                         nx,ny,nz,kr,xm,ym,zm)
*                                                                      *
*   m: the tally number, index of ital.                                *
*                                                                      *
*     Modifeid by T.Furuta on 2025/01/17                               *
*                                                                      *
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

        dimension   kr(mr)
        dimension   nt(nn)
        dimension   ikzz(maxpt,maxnt), iknn(maxpt,maxnt)
        dimension   tm(maxnt+maxpt,2)
        dimension   vl(nr)
        dimension   lr(nr)
        dimension   tr(nr,mz,mn,0:mm,2) ! frtati 2022/02/18 2 -> mm

        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)

*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

        call ttetvl(mr,kr,nr,vl,lr)

        if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

        do 100 ir = 1, nr

           if( itunt(m) .eq. 1 ) then

              cc = 1.d0 / abs(rtfac(m)/facmax(m))

           else if( itunt(m) .eq. 2 ) then

              cc = vl(ir) / abs(rtfac(m)/facmax(m))

           end if

           do 100 iz = 1, mz
           do 100 in = 1, mn
           do 100 il = 0, mm ! frtati 2022/02/18 2 -> mm

                 call invert_stdev(m,A,B,
     &                           tr(ir,iz,in,il,1),
     &                           tr(ir,iz,in,il,2),
     &                           cc)

                 tr(ir,iz,in,il,1) = A
                 tr(ir,iz,in,il,2) = B

  100   continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_yieldrz(m,
     &                        mz,mn,mm,nfr,nfz, ! frtati 2022/02/18 added mm
     &                        nr,nz,nn,rm,zm,nt,ikzz,iknn,
     &                        tr)
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
        dimension   nt(nn)
        dimension   ikzz(maxpt,maxnt), iknn(maxpt,maxnt)
        dimension   tm(maxnt+maxpt,2)
        dimension   fm(nfz,nfr)
        dimension   tr(nr,nz,mz,mn,0:mm,2) ! frtati 2022/02/18 2 -> mm

*-----------------------------------------------------------------------
*        set mesh volume
*-----------------------------------------------------------------------

        vl(jr,jz) = pi * ( rm(jr+1)**2 - rm(jr)**2 )
     &                 * ( zm(jz+1) - zm(jz) )

        pi = 2.d0 * asin(1.d0)

*-----------------------------------------------------------------------

        if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

        do 100 jr = 1, nr
        do 100 jz = 1, nz

          if( itunt(m) .eq. 1 ) then

            cc = 1.d0 / abs(rtfac(m)/facmax(m))

          else if( itunt(m) .eq. 2 ) then

            cc = vl(jr,jz) / abs(rtfac(m)/facmax(m))

          end if

          do 100 iz = 1, mz
          do 100 in = 1, mn
          do 100 il = 0, mm ! frtati 2022/02/18 2 -> mm

             call invert_stdev(m,A,B,
     &                       tr(jr,jz,iz,in,il,1),
     &                       tr(jr,jz,iz,in,il,2),
     &                       cc)

             tr(jr,jz,iz,in,il,1) = A
             tr(jr,jz,iz,in,il,2) = B

  100   continue


      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_yieldxyz(m,
     &                       mz,mn,mm,nf,nl,lt, ! frtati 2022/02/18 added mm
     &                       nx,ny,nz,nn,xm,ym,zm,nt,ikzz,iknn,
     &                       tr)
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

        dimension   lt(nl)
        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)
        dimension   nt(nn)
        dimension   ikzz(maxpt,maxnt), iknn(maxpt,maxnt)
        dimension   tm(maxpt+maxnt,2)
        dimension   fm(nf,nf)
        dimension   tr(nx,ny,nz,mz,mn,0:mm,2) ! frtati 2022/02/18 2 -> mm

*-----------------------------------------------------------------------
*        set mesh volume
*-----------------------------------------------------------------------

         vl(ix,iy,iz) = vls(nl,lt,itmcn(m),itvm(m),itmtr(m,4),
     &                      xm(ix),xm(ix+1),
     &                      ym(iy),ym(iy+1),
     &                      zm(iz),zm(iz+1))

*-----------------------------------------------------------------------

            if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

            do 100 jx = 1, nx
            do 100 jy = 1, ny
            do 100 jz = 1, nz

               if( itunt(m) .eq. 1 ) then

                  cc = 1.d0 / abs(rtfac(m)/facmax(m))

               else if( itunt(m) .eq. 2 ) then

                  cc =  vl(jx,jy,jz) / abs(rtfac(m)/facmax(m))

               end if

               do 100 iz = 1, mz
               do 100 in = 1, mn
               do 100 il = 0, mm ! frtati 2022/02/18 2 -> mm

                  call invert_stdev(m,A,B,
     &                            tr(jx,jy,jz,iz,in,il,1),
     &                            tr(jx,jy,jz,iz,in,il,2),
     &                            cc)

                  tr(jx,jy,jz,iz,in,il,1) = A
                  tr(jx,jy,jz,iz,in,il,2) = B

  100       continue

      end subroutine

