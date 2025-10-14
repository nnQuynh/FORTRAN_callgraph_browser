
************************************************************************
*                                                                      *
      subroutine check_tdpa(m,iax,jsn,jsi,dsin,idsi,ill,ilf,ierr)
*                                                                      *
*   m: the tally number, index of ital.                                *
************************************************************************
      use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05
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


      common /tall08/ rtrx0(itlmax), rtry0(itlmax)

      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)

      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall17/ itndy(itlmax)

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall21/ rtfac(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)

      common /tall28/ itgsh(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall31/ iterl(itlmax)

      common /tall33/ itln(itlmax,2), itli(itlmax,2), itlr(itlmax,2),
     &                rtdm(itlmax,2)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall37/ itcnt(9,itlmax)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall49/ itglt(itlmax)

!OBINATA: for resfile
      common /tall59/ irfll(itlmax), crfln(itlmax), itrff(itlmax)
      character crfln*100
      character irfile*100

      integer :: itdpa
      integer :: iteth
      common /tall80/ itdpa(itlmax)
      common /tall81/ iteth(itlmax)


      dimension idas(mdas*2)
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

*-----------------------------------------------------------------------

      character title*80
      character angelp*200
      character cxtxt*200
      character cytxt*200
      character cztxt*200

*-----------------------------------------------------------------------
      dimension lschn(42), ischn(42)
      character schan(42)*8

      data icsu / 42 /

      data ( schan(i), i = 1, 42 ) /
     &    'mesh    ','part    ','mother  ','axis    ','file    ',
     &    'title   ','angel   ','unit    ','2d-type ','output  ',
     &    'factor  ','material','x-txt   ','y-txt   ','z-txt   ',
     &    'gshow   ','rshow   ','iechrl  ','library ','volmat  ',
     &    'epsout  ','ctmin(1)','ctmax(1)','ctmin(2)','ctmax(2)',
     &    'ctmin(3)','ctmax(3)','resol   ','width   ','trcl    ',
     &    '*trcl   ','gslat   ','resfile ','ginfo   ','bmpout  ',
     &    'vtkout  ','vtkfmt  ','stdcut  ','sangel  ','foamout ',
     &    'idpa    ','ith     '/

      data ( lschn(i), i = 1, 42 ) /
     &     4,         4,         6,         4,         4,
     &     5,         5,         4,         7,         6,
     &     6,         8,         5,         5,         5,
     &     5,         5,         6,         7,         6,
     &     6,         8,         8,         8,         8,
     &     8,         8,         5,         5,         4,
     &     5,         5,         7,         5,         6,
     &     6,         6,         6,         6,         7,
     &     4,         3/

*-----------------------------------------------------------------------


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

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      character dkam*9

      character tname*7
      data      tname /'[t-dpa]'/

      logical deqn1
      logical deqn4
      logical dcom2

      dimension icount(9)
      dimension vtrs(13)

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
            inpat = 0
            inaxi = 0
            infil = 0
            imoth = 0
            jmoth = 1
            imate = 0
            jmate = 1
            iunt  = 1
            iout  = 1
            idtyp = 3

            jtlnn = 0
            jtlin = 0
            jtlrn = 0
            jtlnp = 0
            jtlip = 0
            jtlrp = 0
            dlmxn = -1.0
            dlmxp = -1.0

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

            if( ipm .ne. 4 .and. ipm .ne. 5 .and. ipm .ne. 19 .and.
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

         else if( ipm .eq. 12 ) then

            if( chlw(ic:ic+2) .eq. 'all' ) then

                  imate = 0

            else if( deqn4( chlw(ic:ic) ) ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 981

                  imate = nint( cvvv )

               if( imate .lt. 0 ) then

                  imate = -imate
                  jmate = -1

               end if

               if( imate .eq. 0 ) goto 981

                  nsmte = 1
                  call moddas_allocate_int(imate, ismte_temporary)

                  if( mmmax .gt. mdas ) goto 950

  151                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                        if( jpn  .eq. 3 ) goto 981

                        if( iskip .ne. 0 ) goto 151

                  ic = i1

               do k = 1, imate

                  if( ic .gt. i3 ) then

  152                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( jpn  .eq. 3 ) goto 981

                     if( iskip .ne. 0 ) goto 152

                     ic = i1

                  end if

                     ic = jnumc(chlw,ic,i3)

                     call snum(chlw,ic,i3,ic2,cvvv,ierr)

                     if( ierr .ne. 0 ) goto 981

                     matei = nint( cvvv )

                     ic = ic2

                     ismte_temporary(nsmte-1+k) = matei

               end do
               call moddas_deallocate_int(ismte_temporary)

            else

               goto 981

            end if

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

  145                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                        if( jpn  .eq. 3 ) goto 982

                        if( iskip .ne. 0 ) goto 145

                  ic = i1

               do k = 1, imoth

                  if( ic .gt. i3 ) then

  146                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( jpn  .eq. 3 ) goto 982

                     if( iskip .ne. 0 ) goto 146

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
*        axis
*-----------------------------------------------------------------------

         else if( ipm .eq. 4 ) then

  600       inaxi = inaxi + 1

            if( inaxi .gt. 6 ) goto 991

            if( chlw(ic:ic+2) .eq. 'reg' ) then

               iaxis(inaxi) = 1
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+2) .eq. 'tet' ) then

               iaxis(inaxi) = 10
               ic =jnumc(chlw,ic+4,icl)

            else if( chlw(ic:ic+1) .eq. 'xy' .or.
     &               chlw(ic:ic+1) .eq. 'yx' ) then

               iaxis(inaxi) = 6
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic+1) .eq. 'yz' .or.
     &               chlw(ic:ic+1) .eq. 'zy' ) then

               iaxis(inaxi) = 7
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic+1) .eq. 'zx' .or.
     &               chlw(ic:ic+1) .eq. 'xz' ) then

               iaxis(inaxi) = 8
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic+1) .eq. 'zr' .or.
     &               chlw(ic:ic+1) .eq. 'rz' ) then

               iaxis(inaxi) = 9
               ic =jnumc(chlw,ic+3,icl)

            else if( chlw(ic:ic) .eq. 'x' ) then

               iaxis(inaxi) = 2
               ic =jnumc(chlw,ic+2,icl)

            else if( chlw(ic:ic) .eq. 'y' ) then

               iaxis(inaxi) = 3
               ic =jnumc(chlw,ic+2,icl)

            else if( chlw(ic:ic) .eq. 'z' ) then

               iaxis(inaxi) = 4
               ic =jnumc(chlw,ic+2,icl)

            else if( chlw(ic:ic) .eq. 'r' ) then

               iaxis(inaxi) = 5
               ic =jnumc(chlw,ic+2,icl)

            else

               goto 992

            end if

               if( ic .le. icl ) goto 600

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        output
*-----------------------------------------------------------------------

         else if( ipm .eq. 10 ) then

            if( chlw(ic:ic+2) .eq. 'dpa' ) then

               iout = 1

            else if( chlw(ic:ic+5) .eq. 'simple' ) then

               iout = 2

            else if( chlw(ic:ic+2) .eq. 'all' ) then

               iout = 3

            else

               goto 984

            end if

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

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
*        file name
*-----------------------------------------------------------------------

         else if( ipm .eq. 5 ) then

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

               irfflg = 1

*-----------------------------------------------------------------------
*        title
*-----------------------------------------------------------------------

         else if( ipm .eq. 6 ) then

               ict = min( ic + 79, i2 )

               title = chin(ic:ict)

               titll = ict - ic + 1

*-----------------------------------------------------------------------
*        angel parameters
*-----------------------------------------------------------------------

         else if( ipm .eq. 7 ) then

               ict = min( ic + 199, i2 )

               angelp = chin(ic:ict)

               langel = ict - ic + 1

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

         else if( ipm .eq. 32 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               igslt = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        rshow
*-----------------------------------------------------------------------

         else if( ipm .eq. 17 ) then

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

         else if( ipm .eq. 18 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               iechrl = nint( cvvv )

               if( iechrl .lt. 40 ) goto 997

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        unit
*-----------------------------------------------------------------------

         else if( ipm .eq. 8 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               iunt = nint( cvvv )

               if( iunt .lt. 1 .or. iunt .gt. 2 ) goto 980

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

         else if( ipm .eq. 11 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               rfact = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        volmat
*-----------------------------------------------------------------------

         else if( ipm .eq. 20 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               matvol = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        epsout
*-----------------------------------------------------------------------

         else if( ipm .eq. 21 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ieps = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        resolution
*-----------------------------------------------------------------------

         else if( ipm .eq. 28 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ireso = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        width
*-----------------------------------------------------------------------

         else if( ipm .eq. 29 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               width = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        counter
*-----------------------------------------------------------------------

         else if( ipm .ge. 22 .and. ipm .le. 27 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               icount((ipm-22+2)/2) = 1
               icount( ipm-22+4)    = nint( cvvv )

               if( abs( nint(cvvv) ) .gt. 9999 ) goto 978

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        library
*-----------------------------------------------------------------------

         else if( ipm .eq. 19 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               jtln = nint( cvvv )

               if( jtln .lt. 0 ) goto 971

               call tdpalib(jsn,jsi,dsin,idsi,ill,ilf,
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                      icnp,ienp,dlmax,jtln,jtli,jtlr)


               if( icnp .eq. 1 ) then

                     jtlnp = jtln
                     jtlip = jtli
                     jtlrp = jtlr

                  if( ienp .gt. 0 ) then

                     dlmxp = dlmax

                  end if

               else

                     jtlnn = jtln
                     jtlin = jtli
                     jtlrn = jtlr

                  if( ienp .gt. 0 ) then

                     dlmxn = dlmax

                  end if

               end if

                  if( jpn  .eq. 3 ) goto 800

               goto 150

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

      do i = 1, itaxn(m)

        if ( any( itaxs(m,i) .eq. (/ 6, 7, 8, 9 /))) then
          call check_2dtype(m,idtyp,iec,cepn,lepn,ierr)
        end if

      end do

      call check_factor(m,rfact,iec,cepn,lepn,ierr)

      call check_part(m,iptyp,ipnkf,inpat,iec,cepn,lepn,ierr)

      call check_output(m,iout,iec,cepn,lepn,ierr)

*-----------------------------------------------------------------------

      call judge_tall_check(iec,cepn,lepn,dsin(jsn),idsi(jsn),ierr)

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  950    write(dkam,'(i9)') mdas
         m_err = 'Total tally storage number exceeds mdas ='//dkam
         ErrCha = ''
         ErrID = 'L:1281/R:check_tdpa/F:restdpa.f'
         goto 999

  978    m_err = 'counter should be from -9999 to 9999'
         ErrCha = ''
         ErrID = 'L:1286/R:check_tdpa/F:restdpa.f'
         goto 999

  971    m_err = 'number of library is negative'//tname
         ErrCha = ''
         ErrID = 'L:1291/R:check_tdpa/F:restdpa.f'
         goto 999

  980    m_err = 'unit should be 1 or 2 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1296/R:check_tdpa/F:restdpa.f'
         goto 999

  981    m_err = 'Description of material parameter is wrong in tally '
     &            //tname
         ErrCha = ''
         ErrID = 'L:1302/R:check_tdpa/F:restdpa.f'
         goto 999

  982    m_err = 'Description of mother parameter is wrong in tally '
     &            //tname
         ErrCha = ''
         ErrID = 'L:1308/R:check_tdpa/F:restdpa.f'
         goto 999

  983    m_err = '2D-type should be 1-7 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1313/R:check_tdpa/F:restdpa.f'
         goto 999

  984    m_err = 'Unknown output parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1318/R:check_tdpa/F:restdpa.f'
         goto 999

  986    write(dkam,'(i5)') itlmax
         m_err = 'Total tally number exceeds itlmax = '//dkam
         ErrCha = ''
         ErrID = 'L:1324/R:check_tdpa/F:restdpa.f'
         goto 999

  987    m_err = 'Unknown parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1329/R:check_tdpa/F:restdpa.f'
         goto 999

  989    m_err = 'Double definition of the parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1334/R:check_tdpa/F:restdpa.f'
         goto 999

  990    m_err = 'Too many file name in tally '//tname
         ErrCha = ''
         ErrID = 'L:1339/R:check_tdpa/F:restdpa.f'
         goto 999

  991    m_err = 'Too many axis in tally '//tname
         ErrCha = ''
         ErrID = 'L:1344/R:check_tdpa/F:restdpa.f'
         goto 999

  992    m_err = 'Unknown axis name in tally '//tname
         ErrCha = ''
         ErrID = 'L:1349/R:check_tdpa/F:restdpa.f'
         goto 999

  993    m_err = 'In this line, [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:1354/R:check_tdpa/F:restdpa.f'
         goto 999

  994    m_err = 'Name of particle is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1359/R:check_tdpa/F:restdpa.f'
         goto 999

  995    m_err = 'Number of particles is larger than 6 in tally '//tname
         ErrCha = ''
         ErrID = 'L:1364/R:check_tdpa/F:restdpa.f'
         goto 999

  996    m_err = 'Description of mesh is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1369/R:check_tdpa/F:restdpa.f'
         goto 999

  997    m_err = 'Description of parameter is wrong in tally '//tname
         ErrCha = ''
         ErrID = 'L:1374/R:check_tdpa/F:restdpa.f'
         goto 999

  998    m_err = 'Unknown mesh parameter in tally '//tname
         ErrCha = ''
         ErrID = 'L:1379/R:check_tdpa/F:restdpa.f'
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
      subroutine read_tdpa(m,iax,ierr)
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

        if ( any( itaxs(m,iax) .eq. (/ 6, 7, 8, 9 /) )
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
        call check_tdpa(m,iax,jsn(ioe),jsi(ioe),
     &                 dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &                 ierr)
        endif

        if (ierr .ne. 0 ) goto 900

*-----------------------------------------------------------------------

        rewind(jsi(ioe))
        ill(jsn(ioe),ioe) = 1

*-----------------------------------------------------------------------

        if (itmsh(m) .eq. 1 ) then
          idas0 = nmmax
          idas1 = lmmax
          idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
          idas3 = idas1 + itrgn(m)
     &          + ( itrgn(m) + mod(itrgn(m),2) ) / 2
          idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1
          idasa = idas3
        else if (itmsh(m) .eq. 2 ) then
          idas0 = nmmax
          idas1 = ( lmmax - 1 ) * 2 + 1
          idasa = lmmax
        else if (itmsh(m) .eq. 3 ) then
          idas0 = nmmax
          idasa = lmmax
        end if

*-----------------------------------------------------------------------
*   read tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call read_dpareg(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndy(m),itout(m),
     &      itpan(m),itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &      trRES(irestalm(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        else if ( itmsh(m) .eq. 2 ) then

          call read_dparz(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndy(m),itout(m),
     &      itpan(m),itrnm(m),itznm(m),
     &      das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &      trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 3 ) then

          call read_dpaxyz(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndy(m),itout(m),
     &      itpan(m),itmtn(m),ismte(itmtt(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)),
     &      trRES(irestalm(m)))

        else if( itmsh(m) .eq. 4 ) then

          call read_dpatet(m,
     &      iax,ioe,jsn(ioe),jsi(ioe),
     &      dsin(0,ioe),idsi(0,ioe),ill(0,ioe),ilf(0,ioe),
     &      itndy(m),itout(m),itpan(m),itrgn(m),
     &      trRES(irestalm(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        end if

  800   continue
*-----------------------------------------------------------------------
*   restore tally
*-----------------------------------------------------------------------

        if( itmsh(m) .eq. 1 ) then

          call restore_dpareg(m,
     &      itndy(m),itout(m),
     &      itpan(m),itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &      trRES(irestalm(m)),
     &      itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)))

        else if ( itmsh(m) .eq. 2 ) then

          call restore_dparz(m,
     &      itndy(m),itout(m),
     &      itpan(m),itrnm(m),itznm(m),
     &      das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &      trRES(irestalm(m)))

        else if ( itmsh(m) .eq. 3 ) then

          call restore_dpaxyz(m,
     &      itndy(m),itout(m),
     &      itpan(m),itmtn(m),ismte(itmtt(m)),
     &      itxnm(m),itynm(m),itznm(m),
     &      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),das_itzrg(itzrg(m)),
     &      trRES(irestalm(m)))

        else if( itmsh(m) .eq. 4 ) then

          call restore_dpatet(m,
     &      itndy(m),itout(m),itpan(m),itrgn(m),itrgm(m),
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
      subroutine read_dpareg(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                       nd,ns,np,nr,mr,kr,tr,nvl,ivl,rvl,
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
        dimension   vl(nr)
        dimension   lr(nr)
        dimension   tr(nd,np,nr,2)
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

*-----------------------------------------------------------------------

        character dc2*20

*-----------------------------------------------------------------------
*        reg axis
*-----------------------------------------------------------------------
         if( itaxs(m,iax) .eq. 1 ) then

           dc2  = '#numregvolume'
           ldc2 = 13

           facmx = 1.d0  ! kitamura23/03/31

           if ( ns .eq .1 ) then

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

                read(jsi,'(1x,i5,2x,i5,1pe13.4,1000(1pe13.4,0pf8.4))')
     &               idmm0, idmm1, dmm2,
     &               ((tr(1,ip,ir,k),k=1,2),ip=1,np)

             end do

! sumover
!            ir = 1
!            read(jsi,'(a)') chin
!            read(jsi,'(26x,1000(1pe13.4,0pf8.4))')
!     &               ((tr_sum(1,ip,ir,k),k=1,2),ip=1,np)

           else

             do 190 ip = 1, np

               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                           dc2,ldc2,ierr)

               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900

               do ir = 1, nr

                 read(jsi,'(1x,i5,2x,i5,1pe13.4,6(1pe13.4,0pf8.4))')
     &                idmm0, idmm1, dmm2,
     &                ((tr(ik,ip,ir,k),k=1,2),ik=1,5)

               end do

! sumover
!            ir = 1
!            read(jsi,'(a)') chin
!            read(jsi,'(26x,6(1pe13.4,0pf8.4))')
!     &                ((tr_sum(ik,ip,ir,k),k=1,2),ik=1,5)


               if ( ns .eq .3 ) then

                 call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                             dc2,ldc2,ierr)

                 if ( ierr  .ne. 0 ) goto 900
                 if ( jpn   .eq. 3 ) goto 900

                 do ir = 1, nr

                    read(jsi,
     &                  '(1x,i5,2x,i5,1pe13.4,6(1pe13.4,0pf8.4))')
     &                   idmm0, idmm1, dmm2,
     &                   ((tr(ik,ip,ir,k),k=1,2),ik=6,10)

                 end do

               end if

  190        continue

           end if

*-----------------------------------------------------------------------
         end if

  900 continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_dpatet(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                       nd,ns,np,nr,tr,
     &                       nx,ny,nz,xm,ym,zm)
*                                                                      *
*   m: the tally number, index of ital.                                *
*                                                                      *
*     Created by T.Furuta on 2019/01/21                                *
*                                                                      *
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

        dimension   vl(nr)
        dimension   lr(nr)
        dimension   tr(nd,np,nr,2)

        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

        character chin*200, chlw*200, chcm*200

*-----------------------------------------------------------------------

        character dc2*20

*-----------------------------------------------------------------------
*        tet axis
*-----------------------------------------------------------------------
         if( itaxs(m,iax) .eq. 10 ) then

           dc2  = '#numtetravolume'
           ldc2 = 13

           facmx = 1.d0  ! kitamura23/03/31

           if ( ns .eq .1 ) then

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

                read(jsi,'(1x,i8,2x,i8,1pe13.4,1000(1pe13.4,0pf8.4))')
     &               idmm0, idmm1, dmm2,
     &               ((tr(1,ip,ir,k),k=1,2),ip=1,np)

             end do

! sumover
!            ir = 1
!            read(jsi,'(a)') chin
!            read(jsi,'(32x,1000(1pe13.4,0pf8.4))')
!     &               ((tr_sum(1,ip,ir,k),k=1,2),ip=1,np)


           else

             do 190 ip = 1, np

               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                           dc2,ldc2,ierr)

               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900

               do ir = 1, nr

                 read(jsi,'(1x,i8,2x,i8,1pe13.4,6(1pe13.4,0pf8.4))')
     &                idmm0, idmm1, dmm2,
     &                ((tr(ik,ip,ir,k),k=1,2),ik=1,5)

               end do

! sumover
!               ir = 1
!               read(jsi,'(a)') chin
!               read(jsi,'(32x,6(1pe13.4,0pf8.4))')
!     &                ((tr_sum(ik,ip,ir,k),k=1,2),ik=1,5)


               if ( ns .eq .3 ) then

                 call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                             dc2,ldc2,ierr)

                 if ( ierr  .ne. 0 ) goto 900
                 if ( jpn   .eq. 3 ) goto 900

                 do ir = 1, nr

                    read(jsi,
     &                  '(1x,i8,2x,i8,1pe13.4,6(1pe13.4,0pf8.4))')
     &                   idmm0, idmm1, dmm2,
     &                   ((tr(ik,ip,ir,k),k=1,2),ik=6,10)

                 end do

               end if

  190        continue

           end if

*-----------------------------------------------------------------------
         end if

  900 continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_dparz(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                      nd,ns,np,nr,nz,rm,zm,tr)
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
        dimension   tr(nd,np,nr,nz,2)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

        character chin*200, chlw*200, chcm*200

*-----------------------------------------------------------------------

        character dc2*20

*-----------------------------------------------------------------------
*        r axis
*-----------------------------------------------------------------------
         if( itaxs(m,iax) .eq. 5 ) then

           dc2  = '#r-lowerr-upper'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31


           if ( ns .eq .1 ) then

           nr_0 = 0
           nz_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_prz.inc'

            do iz = 1, nz, nzstepi
            do ip = 1, np, npstepi
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
     &               dmm0, dmm1,
     &      ((((tr(1,ip+iploop-1,ir+irloop-1,
     &               iz+izloop-1,k),k=1,2),
     &               iploop=1,npstepi),irloop=1,nrstepi),
     &               izloop=1,nzstepi)

             end do

! sumover
!               ir = 1
!               read(jsi,'(a)') chin
!               read(jsi,'(26x,1000(1pe13.4,0pf8.4))')
!     &      ((((tr_sum(1,ip+iploop-1,ir+irloop-1,
!     &               iz+izloop-1,k),k=1,2),
!     &               iploop=1,npstepi),irloop=1,nrstepi),
!     &               izloop=1,nzstepi)


            end do
            end do

           else

            do iz = 1, nz

             do ip = 1, np

               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                           dc2,ldc2,ierr)

               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900

               do ir = 1, nr

                 read(jsi,'(1p2e13.4,6(1pe13.4,0pf8.4))')
     &                dmm0, dmm1,
     &                ((tr(ik,ip,ir,iz,k),k=1,2),ik=1,5)

               end do

! sumover
!               ir = 1
!               read(jsi,'(a)') chin
!               read(jsi,'(26x,6(1pe13.4,0pf8.4))')
!     &                ((tr_sum(ik,ip,ir,iz,k),k=1,2),ik=1,5)


               if ( ns .eq .3 ) then

                 call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                             dc2,ldc2,ierr)

                 if ( ierr  .ne. 0 ) goto 900
                 if ( jpn   .eq. 3 ) goto 900

                 do ir = 1, nr

                    read(jsi,'(1p2e13.4,6(1pe13.4,0pf8.4))')
     &                   dmm0, dmm1,
     &                   ((tr(ik,ip,ir,iz,k),k=1,2),ik=6,10)

                 end do

! sumover
!               ir = 1
!               read(jsi,'(a)') chin
!               read(jsi,'(26x,6(1pe13.4,0pf8.4))')
!     &                   ((tr_sum(ik,ip,ir,iz,k),k=1,2),ik=6,10)


               end if

             end do
            end do

           end if

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 4 ) then

           dc2  = '#z-lowerz-upper'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31


           if ( ns .eq .1 ) then

           nr_0 = 0
           nz_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_prz.inc'

            do ir = 1, nr, nrstepi
            do ip = 1, np, npstepi
             call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx
            write(*,*)' DPA R-Z facmax Z-IN ',facmax(m),'',m

             call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                         dc2,ldc2,ierr)

             if ( ierr  .ne. 0 ) goto 900
             if ( jpn   .eq. 3 ) goto 900

             nzstepi = 1
             do iz = 1, nz, nzstepi


                read(jsi,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &               dmm0, dmm1,
     &      ((((tr(1,ip+iploop-1,ir+irloop-1,
     &               iz+izloop-1,k),k=1,2),
     &               iploop=1,npstepi),irloop=1,nrstepi),
     &               izloop=1,nzstepi)

             end do

! sumover
!               iz = 1
!               read(jsi,'(a)') chin
!               read(jsi,'(26x,1000(1pe13.4,0pf8.4))')
!     &      ((((tr_sum(1,ip+iploop-1,ir+irloop-1,
!     &               iz+izloop-1,k),k=1,2),
!     &               iploop=1,npstepi),irloop=1,nrstepi),
!     &               izloop=1,nzstepi)

            end do
            end do

           else

            do ir = 1, nr

             do ip = 1, np

               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                           dc2,ldc2,ierr)

               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900

               do iz = 1, nz

                 read(jsi,'(1p2e13.4,6(1pe13.4,0pf8.4))')
     &                dmm0, dmm1,
     &                ((tr(ik,ip,ir,iz,k),k=1,2),ik=1,5)

               end do

! sumover
!               iz = 1
!               read(jsi,'(a)') chin
!               read(jsi,'(26x,6(1pe13.4,0pf8.4))')
!     &                ((tr_sum(ik,ip,ir,iz,k),k=1,2),ik=1,5)


               if ( ns .eq .3 ) then

                 call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                             dc2,ldc2,ierr)

                 if ( ierr  .ne. 0 ) goto 900
                 if ( jpn   .eq. 3 ) goto 900

                 do iz = 1, nz

                    read(jsi,'(1p2e13.4,6(1pe13.4,0pf8.4))')
     &                   dmm0, dmm1,
     &                   ((tr(ik,ip,ir,iz,k),k=1,2),ik=6,10)

                 end do

! sumover
!               iz = 1
!               read(jsi,'(a)') chin
!               read(jsi,'(26x,6(1pe13.4,0pf8.4))')
!     &                   ((tr_sum(ik,ip,ir,iz,k),k=1,2),ik=6,10)


               end if

             end do
            end do

           end if


*-----------------------------------------------------------------------
*        rz axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 9 ) then

           call sel_dc2(m,dc2,ldc2)
           if( ittwo(m) .eq. 4 ) then
             dc2  = '#rzdpar.err'
             ldc2 = 11
           else if( ittwo(m) .eq. 5 ) then
             dc2  = 'r/z'
             ldc2 = 3
           end if

           facmx = 1.d0  ! kitamura23/03/31

           do 390 ip = 1, np
           do 390 ik = 1, nd

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
     &         ( ( tr(ik,ip,ir,iz,ioe), iz = 1, nz ), ir = nr, 1, -1 )

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 4 ) then

               do iz = 1, nz
               do ir = 1, nr

                  read(jsi,'(1p3e11.3,0pf8.4)') dmm0, dmm1,
     &               tr(ik,ip,ir,iz,1), tr(ik,ip,ir,iz,2)

               end do
               end do

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 5 ) then

               do ir = nr, 1, -1

                  read(jsi,'(1p1000e11.3)') dmm0,
     &            ( tr(ik,ip,ir,iz,ioe), iz = 1, nz )

               end do

*-----------------------------------------------------------------------
            end if

  390    continue

*-----------------------------------------------------------------------
         end if

  900 continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_dpaxyz(m,iax,ioe,jsn,jsi,dsin,idsi,ill,ilf,
     &                       nd,ns,np,nl,lt,
     &                       nx,ny,nz,xm,ym,zm,tr)
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

        dimension   lt(nl)
        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)
        dimension   tr(nd,np,nx,ny,nz,2)

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

        character chin*200, chlw*200, chcm*200

*-----------------------------------------------------------------------

        character dc2*20

*-----------------------------------------------------------------------

            npg = np
            ndg = nd

*-----------------------------------------------------------------------
*        x axis
*-----------------------------------------------------------------------
         if( itaxs(m,iax) .eq. 2 ) then

           dc2  = '#x-lowerx-upper'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31


           if ( ns .eq .1 ) then

           nx_0 = 0
           ny_0 = 0
           nz_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_pxyz.inc'

           do iy = 1, ny, nystepi
           do iz = 1, nz, nzstepi
           do ip = 1, np, npstepi
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
     &               dmm0, dmm1,
     &      (((((tr(1,ip+iploop-1,ix+ixloop-1,iy+iyloop-1,
     &                iz+izloop-1,k),k=1,2),
     &                iploop=1,npstepi),ixloop=1,nxstepi),
     &                iyloop=1,nystepi),izloop=1,nzstepi)

             end do

! sumover
!               ix = 1
!               read(jsi,'(a)') chin
!               read(jsi,'(26x,1000(1pe13.4,0pf8.4))')
!     &      (((((tr_sum(1,ip+iploop-1,ix+ixloop-1,iy+iyloop-1,
!     &                iz+izloop-1,k),k=1,2),
!     &                iploop=1,npstepi),ixloop=1,nxstepi),
!     &                iyloop=1,nystepi),izloop=1,nzstepi)

           end do
           end do
           end do

           else

           do iy = 1, ny
           do iz = 1, nz

             do ip = 1, np

               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                           dc2,ldc2,ierr)

               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900

               do ix = 1, nx

                 read(jsi,'(1p2e13.4,6(1pe13.4,0pf8.4))')
     &                dmm0, dmm1,
     &                ((tr(ik,ip,ix,iy,iz,k),k=1,2),ik=1,5)

               end do

! sumover
!               ix = 1
!               read(jsi,'(a)') chin
!               read(jsi,'(26x,6(1pe13.4,0pf8.4))')
!     &                ((tr_sum(ik,ip,ix,iy,iz,k),k=1,2),ik=1,5)


               if ( ns .eq .3 ) then

                 call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                             dc2,ldc2,ierr)

                 if ( ierr  .ne. 0 ) goto 900
                 if ( jpn   .eq. 3 ) goto 900

                 do ix = 1, nx

                    read(jsi,'(1p2e13.4,6(1pe13.4,0pf8.4))')
     &                   dmm0, dmm1,
     &                   ((tr(ik,ip,ix,iy,iz,k),k=1,2),ik=6,10)

                 end do

! sumover
!               ix = 1
!               read(jsi,'(a)') chin
!               read(jsi,'(26x,6(1pe13.4,0pf8.4))')
!     &              ((tr_sum(ik,ip,ix,iy,iz,k),k=1,2),ik=6,10)


               end if

             end do
            end do
            end do

           end if


*-----------------------------------------------------------------------
*        y axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 3 ) then

           dc2  = '#y-lowery-upper'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31

           if ( ns .eq .1 ) then

           nx_0 = 0
           ny_0 = 0
           nz_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_pxyz.inc'

           do ix = 1, nx, nxstepi
           do iz = 1, nz, nzstepi
           do ip = 1, np, npstepi
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
     &               dmm0, dmm1,
     &      (((((tr(1,ip+iploop-1,ix+ixloop-1,iy+iyloop-1,
     &                iz+izloop-1,k),k=1,2),
     &                iploop=1,npstepi),ixloop=1,nxstepi),
     &                iyloop=1,nystepi),izloop=1,nzstepi)

             end do

! sumover
!               iy = 1
!               read(jsi,'(a)') chin
!               read(jsi,'(26x,1000(1pe13.4,0pf8.4))')
!     &      (((((tr_sum(1,ip+iploop-1,ix+ixloop-1,iy+iyloop-1,
!     &                iz+izloop-1,k),k=1,2),
!     &                iploop=1,npstepi),ixloop=1,nxstepi),
!     &                iyloop=1,nystepi),izloop=1,nzstepi)

            end do
            end do
            end do

           else
            do ix = 1, nx
            do iz = 1, nz

             do ip = 1, np

               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                           dc2,ldc2,ierr)

               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900

               do iy = 1, ny

                 read(jsi,'(1p2e13.4,6(1pe13.4,0pf8.4))')
     &                dmm0, dmm1,
     &                ((tr(ik,ip,ix,iy,iz,k),k=1,2),ik=1,5)

               end do

! sumover
!               iy = 1
!               read(jsi,'(a)') chin
!               read(jsi,'(26x,6(1pe13.4,0pf8.4))')
!     &                ((tr_sum(ik,ip,ix,iy,iz,k),k=1,2),ik=1,5)


               if ( ns .eq .3 ) then

                 call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                             dc2,ldc2,ierr)

                 if ( ierr  .ne. 0 ) goto 900
                 if ( jpn   .eq. 3 ) goto 900

                 do iy = 1, ny

                    read(jsi,'(1p2e13.4,6(1pe13.4,0pf8.4))')
     &                   dmm0, dmm1,
     &                   ((tr(ik,ip,ix,iy,iz,k),k=1,2),ik=6,10)

                 end do

! sumover
!               iy = 1
!               read(jsi,'(a)') chin
!               read(jsi,'(26x,6(1pe13.4,0pf8.4))')
!     &              ((tr_sum(ik,ip,ix,iy,iz,k),k=1,2),ik=6,10)


               end if

             end do
            end do
            end do

           end if


*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 4 ) then

           dc2  = '#z-lowerz-upper'
           ldc2 = 15

           facmx = 1.d0  ! kitamura23/03/31

           if ( ns .eq .1 ) then

           nx_0 = 0
           ny_0 = 0
           nz_0 = 0
           np_0 = 0
      include 'samepage_include/samepagestepi_pxyz.inc'

           do ix = 1, nx, nxstepi
           do iy = 1, ny, nystepi
           do ip = 1, np, npstepi
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
     &               dmm0, dmm1,
     &      (((((tr(1,ip+iploop-1,ix+ixloop-1,iy+iyloop-1,
     &                iz+izloop-1,k),k=1,2),
     &                iploop=1,npstepi),ixloop=1,nxstepi),
     &                iyloop=1,nystepi),izloop=1,nzstepi)

             end do

! sumover
!               iz = 1
!               read(jsi,'(a)') chin
!               read(jsi,'(26x,1000(1pe13.4,0pf8.4))')
!     &      (((((tr_sum(1,ip+iploop-1,ix+ixloop-1,iy+iyloop-1,
!     &                iz+izloop-1,k),k=1,2),
!     &                iploop=1,npstepi),ixloop=1,nxstepi),
!     &                iyloop=1,nystepi),izloop=1,nzstepi)


            end do
            end do
            end do

           else

            do ix = 1, nx
            do iy = 1, ny

             do ip = 1, np

               call seek_resfile_np(jsn,jsi,dsin,idsi,ill,ilf,jpn,ierr)

            if ( irestart .eq. 1 .and. rtfac(m) .lt. 0.d0 )
     &      call seek_resfile_fa(jsn,jsi,dsin,idsi,ill,ilf,jpn,facmx,
     &                          ierr)

            facmax(m) = facmx

               call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                           dc2,ldc2,ierr)

               if ( ierr  .ne. 0 ) goto 900
               if ( jpn   .eq. 3 ) goto 900

               do iz = 1, nz

                 read(jsi,'(1p2e13.4,6(1pe13.4,0pf8.4))')
     &                dmm0, dmm1,
     &                ((tr(ik,ip,ix,iy,iz,k),k=1,2),ik=1,5)

               end do

! sumover
!               iz = 1
!               read(jsi,'(a)') chin
!               read(jsi,'(26x,6(1pe13.4,0pf8.4))')
!     &                ((tr_sum(ik,ip,ix,iy,iz,k),k=1,2),ik=1,5)


               if ( ns .eq .3 ) then

                 call seek_resfile(jsn,jsi,dsin,idsi,ill,ilf,jpn,
     &                             dc2,ldc2,ierr)

                 if ( ierr  .ne. 0 ) goto 900
                 if ( jpn   .eq. 3 ) goto 900

                 do iz = 1, nz

                    read(jsi,'(1p2e13.4,6(1pe13.4,0pf8.4))')
     &                   dmm0, dmm1,
     &                   ((tr(ik,ip,ix,iy,iz,k),k=1,2),ik=6,10)

                 end do

! sumover
!               iz = 1
!               read(jsi,'(a)') chin
!               read(jsi,'(26x,6(1pe13.4,0pf8.4))')
!     &                   ((tr_sum(ik,ip,ix,iy,iz,k),k=1,2),ik=6,10)

               end if

             end do
            end do
            end do

           end if


*-----------------------------------------------------------------------
*        xy axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 6 ) then

           call sel_dc2(m,dc2,ldc2)
           if( ittwo(m) .eq. 4 ) then
             dc2  = '#xydpar.err'
             ldc2 = 11
           else if( ittwo(m) .eq. 5 ) then
             dc2  = 'y/x'
             ldc2 = 3
           end if

           facmx = 1.d0  ! kitamura23/03/31

           do 490 iz = 1, nz
           do 490 ip = 1, npg
           do 490 ik = 1, ndg

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
     &         (( tr(ik,ip,ix,iy,iz,ioe), ix = 1, nx ), iy = ny, 1, -1 )

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 4 ) then

               do iy = 1, ny
               do ix = 1, nx

                  read(jsi,'(1p3e11.3,0pf8.4)') dmm0, dmm1,
     &               tr(ik,ip,ix,iy,iz,1), tr(ik,ip,ix,iy,iz,2)

               end do
               end do

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 5 ) then

               do iy = ny, 1, -1

                  read(jsi,'(1p1000e11.3)') dmm0,
     &            ( tr(ik,ip,ix,iy,iz,ioe), ix = 1, nx )

               end do

*-----------------------------------------------------------------------
            end if

  490    continue

*-----------------------------------------------------------------------
*        yz axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 7 ) then

           call sel_dc2(m,dc2,ldc2)
           if( ittwo(m) .eq. 4 ) then
             dc2  = '#zydpar.err'
             ldc2 = 11
           else if( ittwo(m) .eq. 5 ) then
             dc2  = 'y/z'
             ldc2 = 3
           end if

           do 590 ix = 1, nx
           do 590 ip = 1, npg
           do 590 ik = 1, ndg

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
     &         (( tr(ik,ip,ix,iy,iz,ioe), iz = 1, nz ), iy = ny, 1, -1 )

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 4 ) then

               do iy = 1, ny
               do iz = 1, nz

                  read(jsi,'(1p3e11.3,0pf8.4)') dmm0, dmm1,
     &               tr(ik,ip,ix,iy,iz,1), tr(ik,ip,ix,iy,iz,2)

               end do
               end do

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 5 ) then

               do iy = ny, 1, -1

                  read(jsi,'(1p1000e11.3)') dmm0,
     &            ( tr(ik,ip,ix,iy,iz,ioe), iz = 1, nz )

               end do

*-----------------------------------------------------------------------
            end if

  590    continue

*-----------------------------------------------------------------------
*        yz axis (matrix)
*-----------------------------------------------------------------------
         else if( itaxs(m,iax) .eq. 8 ) then

           call sel_dc2(m,dc2,ldc2)
           if( ittwo(m) .eq. 4 ) then
             dc2  = '#zxdpar.err'
             ldc2 = 11
           else if( ittwo(m) .eq. 5 ) then
             dc2  = 'x/z'
             ldc2 = 3
           end if

           facmx = 1.d0  ! kitamura23/03/31

           do 690 iy = 1, ny
           do 690 ip = 1, npg
           do 690 ik = 1, ndg

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
     &         (( tr(ik,ip,ix,iy,iz,ioe), iz = 1, nz ), ix = nx, 1, -1 )

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 4 ) then

               do ix = 1, nx
               do iz = 1, nz

                  read(jsi,'(1p3e11.3,0pf8.4)') dmm0, dmm1,
     &               tr(ik,ip,ix,iy,iz,1), tr(ik,ip,ix,iy,iz,2)

               end do
               end do

*-----------------------------------------------------------------------
            else if( ittwo(m) .eq. 5 ) then

               do ix = nx, 1, -1

                  read(jsi,'(1p1000e11.3)') dmm0,
     &            ( tr(ik,ip,ix,iy,iz,ioe), iz = 1, nz )

               end do

*-----------------------------------------------------------------------
            end if

  690    continue

*-----------------------------------------------------------------------
         end if


  900 continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_dpareg(m,
     &                          nd,ns,np,nr,mr,kr,tr,nvl,ivl,rvl,
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
        dimension   vl(nr)
        dimension   lr(nr)
        dimension   tr(nd,np,nr,2)
        dimension   ivl(nvl)
        dimension   rvl(nvl)

        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)

*-----------------------------------------------------------------------

            npg = np
            ndg = nd

*-----------------------------------------------------------------------
*           itunt(m) = 1 : DPA* 1.E-24/source
*                    = 2 : DPA/source
*-----------------------------------------------------------------------

            if( itunt(m) .eq. 1 ) then

               dnor = 1.0

            else if( itunt(m) .eq. 2 ) then

               dnor = 1.d+24

            end if

*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

               call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)

*-----------------------------------------------------------------------
*        relative error and  unit conversion
*        c1 : nomalization for source
*-----------------------------------------------------------------------

            if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

            do 100 ir = 1, nr
            do 100 ip = 1, np
            do 100 ik = 1, nd

                  call invert_stdev(m,A,B,
     &                            tr(ik,ip,ir,1),
     &                            tr(ik,ip,ir,2),

     &             dnor*vl(ir)/abs(rtfac(m)/facmax(m)))

                  tr(ik,ip,ir,1) = A
                  tr(ik,ip,ir,2) = B

  100   continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_dpatet(m,
     &                          nd,ns,np,nr,mr,tr,
     &                          nx,ny,nz,kr,xm,ym,zm)
*                                                                      *
*   m: the tally number, index of ital.                                *
*                                                                      *
*     Modified by T.Furuta on 2025/01/17                               *
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

*-----------------------------------------------------------------------

        dimension   kr(mr)
        dimension   vl(nr)
        dimension   lr(nr)
        dimension   tr(nd,np,nr,2)

        dimension   xm(nx+1)
        dimension   ym(ny+1)
        dimension   zm(nz+1)

*-----------------------------------------------------------------------

            npg = np
            ndg = nd

*-----------------------------------------------------------------------
*           itunt(m) = 1 : DPA* 1.E-24/source
*                    = 2 : DPA/source
*-----------------------------------------------------------------------

            if( itunt(m) .eq. 1 ) then

               dnor = 1.0

            else if( itunt(m) .eq. 2 ) then

               dnor = 1.d+24

            end if

*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

               call ttetvl(mr,kr,nr,vl,lr)

*-----------------------------------------------------------------------
*        relative error and  unit conversion
*        c1 : nomalization for source
*-----------------------------------------------------------------------

            if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

            do 100 ir = 1, nr
            do 100 ip = 1, np
            do 100 ik = 1, nd

                  call invert_stdev(m,A,B,
     &                            tr(ik,ip,ir,1),
     &                            tr(ik,ip,ir,2),
     &             dnor*vl(ir)/abs(rtfac(m)/facmax(m)))

                  tr(ik,ip,ir,1) = A
                  tr(ik,ip,ir,2) = B

  100   continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_dparz(m,
     &                         nd,ns,np,nr,nz,rm,zm,tr)
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
        dimension   tr(nd,np,nr,nz,2)

*-----------------------------------------------------------------------
*        set mesh volume
*-----------------------------------------------------------------------

        vl(ir,iz) = pi * ( rm(ir+1)**2 - rm(ir)**2 )
     &                 * ( zm(iz+1) - zm(iz) )

        pi = 2.d0 * asin(1.d0)

*-----------------------------------------------------------------------
*           itunt(m) = 1 : DPA* 1.E-24/source
*                    = 2 : DPA/source
*-----------------------------------------------------------------------

            if( itunt(m) .eq. 1 ) then

               dnor = 1.0

            else if( itunt(m) .eq. 2 ) then

               dnor = 1.d+24

            end if

*-----------------------------------------------------------------------
*        relative error and  unit conversion
*-----------------------------------------------------------------------

            if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

            do 100 iz = 1, nz
            do 100 ir = 1, nr
            do 100 ip = 1, np
            do 100 ik = 1, nd

                  call invert_stdev(m,A,B,
     &                            tr(ik,ip,ir,iz,1),
     &                            tr(ik,ip,ir,iz,2),
     &             vl(ir,iz)*dnor/abs(rtfac(m)/facmax(m)))

                  tr(ik,ip,ir,iz,1) = A
                  tr(ik,ip,ir,iz,2) = B


  100       continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine restore_dpaxyz(m,
     &                          nd,ns,np,nl,lt,
     &                          nx,ny,nz,xm,ym,zm,tr)
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
        dimension   tr(nd,np,nx,ny,nz,2)

*-----------------------------------------------------------------------
*        set mesh volume
*-----------------------------------------------------------------------

         vl(ix,iy,iz) = vls(nl,lt,itmcn(m),itvm(m),itmtr(m,4),
     &                      xm(ix),xm(ix+1),
     &                      ym(iy),ym(iy+1),
     &                      zm(iz),zm(iz+1))

*-----------------------------------------------------------------------

        npg = np
        ndg = nd

*-----------------------------------------------------------------------
*           itunt(m) = 1 : DPA* 1.E-24/source
*                    = 2 : DPA/source
*-----------------------------------------------------------------------

            if( itunt(m) .eq. 1 ) then

               dnor = 1.0

            else if( itunt(m) .eq. 2 ) then

               dnor = 1.d+24

            end if

*-----------------------------------------------------------------------
*        relative error and  unit conversion
*        c1 : nomalization for source
*-----------------------------------------------------------------------

            if( facmax(m) .eq. 0.d0) facmax(m) = 1.d0  ! kitamura23/03/31

            do 100 iz = 1, nz
            do 100 iy = 1, ny
            do 100 ix = 1, nx
            do 100 ip = 1, np
            do 100 ik = 1, nd

               call invert_stdev(m,A,B,
     &                         tr(ik,ip,ix,iy,iz,1),
     &                         tr(ik,ip,ix,iy,iz,2),
     &          vl(ix,iy,iz)*dnor/abs(rtfac(m)/facmax(m)))

               tr(ik,ip,ix,iy,iz,1) = A
               tr(ik,ip,ix,iy,iz,2) = B

  100       continue

*-----------------------------------------------------------------------

      end subroutine


