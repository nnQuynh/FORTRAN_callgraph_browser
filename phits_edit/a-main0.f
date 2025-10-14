************************************************************************
*                                                                      *
      subroutine a_main0(icang,iindat)
*                                                                      *
*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

************************************************************************
*                                                                      *
*         ANGEL         Graphic Program                                *
                        parameter ( version = 4.51 )
*                       2005/11/11 last reviced                        *
*         Copyright (C) Koji NIITA  1993-2005. All rights reserved.    *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /stwhb/  stwhb(0:13,3,0:255)
      common /stwrl/  stwrl4(3,4)
      common /stwhbk/ stwhbk(3)

*-----------------------------------------------------------------------

      common /box/  bxws, bxwl, bxss, bxsl, bxds, bxdl, bxls, bxll
      common /suf/  ftht, ftss, ftks, ftsv,
     &              yplu, ypld, ypku, ypkd, yplh, ypls
      common /con/  cm, dd
      common /pap/  a4w, a4h, wmg, hmg, wct, hct

      common /frm/  xal, yal

*-----------------------------------------------------------------------

      common /wtval1/ strl0, strh0, strb0
      common /wtval2/ clc(0:inig,3), fts(0:inig),
     &                stw(0:inig), sth(0:inig), stb(0:inig),
     &                xps(0:inig), yps(0:inig),
     &                xpp(0:inig), ypp(0:inig),
     &                ixs(0:inig), iys(0:inig), jfn(0:inig)

*-----------------------------------------------------------------------

      real(8),allocatable:: daxy(:,:),dax(:),day(:)

*-----------------------------------------------------------------------

      common /wtval3/ strla(0:inig), strha(0:inig), strba(0:inig)
      common /wtval7/ strhtp, strbtp, strhbt, strbbt, strhall
      common /wtval8/ xpsi(0:inig), ypsi(0:inig), strl(0:inig),
     &                xpsl(0:inig), xpsr(0:inig), rind(0:inig)

      common /fdfn/ ifd(0:13),ikan1,ikan2,ikan3,ikan4,
     &              idsm(9), idln(6), idxe, idye, ibd1, ibd2, ibd3,
     &              ibd4, ibd5, ibd6, ibd7, ibd8, ibd9, ibd10, ibd11

      common /lbox/  bcb(3), bcl(3), bcs(3), bds, bdd, bdc

      common /bnbox/ xog0, yog0, rog0, sxg0, syg0,
     &               bx01, bx02, by01, by02, icb

*-----------------------------------------------------------------------

      common /h2/ xdd,ydd,
     &            xymax,xymin,clus,cmax,cmin,cxymax,cxymin,xyabs,
     &            dmax,dmin,ixnm,iynm

      common /h3/ smax(3), smin(3)

      common /h5/ nhc, mhc(100)
      common /h6/ dhgh(100), dwih(100)

*-----------------------------------------------------------------------

      common /error/ m_err, l_err, k_err
      character m_err*200

*-----------------------------------------------------------------------

      character     today1*11, timestr1*5, today*100, timestr*100
      character     infnm(ichrl)*1, infn(ichrl)*1, avers*5
      common /dfil1/today, timestr, today1, timestr1, infnm, infn, avers
      common /dfil2/jtoday, jtime, jnm, mpage, iname, javers

*-----------------------------------------------------------------------

      character iindat(200)*1

      character iofn(ichrl)*1
      character indat(ichrl)*1

*-----------------------------------------------------------------------

      common /rval1/ cval(mxcval), aval(mxcval)

*-----------------------------------------------------------------------

      character dum(ichrl)*1
      character lum(ichrl)*1
      character ccmt(ichrl)*1

*-----------------------------------------------------------------------

      character   jikn(3)*18, funn(3)*18,
     &            nenn(3)*18, tuki(3)*18, niti(3)*18

      character   month*2

      dimension   ijikn(3,2), ifunn(3,2),
     &            inenn(3,2), ituki(3,2), initi(3,2)

      data (ijikn(1,i),i=1,2) / 142, 158 /
      data (ijikn(2,i),i=1,2) /  59, 126 /
      data (ijikn(3,i),i=1,2) / 187, 254 /

      data (ifunn(1,i),i=1,2) / 149, 170 /
      data (ifunn(2,i),i=1,2) /  74,  44 /
      data (ifunn(3,i),i=1,2) / 202, 172 /

      data (inenn(1,i),i=1,2) / 148,  78 /
      data (inenn(2,i),i=1,2) /  71,  47 /
      data (inenn(3,i),i=1,2) / 199, 175 /

      data (ituki(1,i),i=1,2) / 140, 142 /
      data (ituki(2,i),i=1,2) /  55, 110 /
      data (ituki(3,i),i=1,2) / 183, 238 /

      data (initi(1,i),i=1,2) / 147, 250 /
      data (initi(2,i),i=1,2) /  49,  49 /
      data (initi(3,i),i=1,2) / 198, 252 /

      dimension itopi(8), iboti(8)

      data itopi / 92,  51,  55,  55,  92,  48,  48,  49 /
      data iboti / 92,  51,  55,  55,  92,  48,  48,  48 /

*-----------------------------------------------------------------------

      character mcm(6,ichrl)*1
      dimension imcm(6)

*-----------------------------------------------------------------------

      dimension xlta(numtic), xsta(numtic), ylta(numtic), ysta(numtic)
      dimension nxtch(numtic), nytch(numtic)
      character xcta(numtic,30)*1,ycta(numtic,30)*1

      character dsjht*200
      character dsjhe*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      character xtt(ichrl)*1,ytt(ichrl)*1,htt(ichrl)*1

      character xttn(ichrl)*1
      character yttn(ichrl)*1
      character cmch(-5:7)*5, tnch(-9:3)*7
      dimension icmc(-5:7), itnc(-9:3)
      character scch(4)*4
      dimension iscc(4)

*-----------------------------------------------------------------------

      character dpin(ipsm)*200
      dimension idpi(ipsm)
      dimension psbbx(ipsm,4)
      dimension xyps(ipsm,6)
      dimension ixps(ipsm), iyps(ipsm)

*-----------------------------------------------------------------------

      dimension wxys(mc,9),iwx(mc),iwy(mc),iwn(mc),iwf(mc),iwb(mc)
      dimension ibox(mc),cbox(mc,3,3)

      dimension axys(mc,11),iax(mc),ian(mc),iaf(mc)

      dimension bxys(mc,10)
      dimension cxys(mc,13)
      dimension pxys(mc,14)
      dimension oxys(mc,16)
      dimension rxys(mc,16)
      dimension sxys(mc,15)

      dimension iycm(mc,8),rycm(mc,7)
      dimension rcut(mc,5)
      dimension iwd2(mc)

*-----------------------------------------------------------------------

      dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

      dimension  clal(3), clax(3), cltx(3), clnm(3), cltl(3), cllg(3),
     &           clpa(3), clcn(3), clfr(3), clms(3), clgb(3), clgl(3),
     &           clgs(3), clpb(3), clpl(3), clps(3), clnb(3), clnl(3),
     &           clns(3), clhd(3), clbg(3), clin(3)

      dimension  rcol(3)
      dimension  rcoll(3), rcolc(3), rcolb(3)

*-----------------------------------------------------------------------

      logical   exex

*-----------------------------------------------------------------------

      character devnul*13
      data devnul/ ' >/dev/null &'/

*-----------------------------------------------------------------------

      character,allocatable:: chag(:,:)

      character yen*1
      character tub*1

*-----------------------------------------------------------------------

      common /cmap/ ndis, icrev, cmap
      character cmap*99

*-----------------------------------------------------------------------

      allocate(chag(inig,0:ichrl))

*-----------------------------------------------------------------------

            yen = char(92)
            tub = char(9)

            scch(1)  = 'min'
            iscc(1)  = 3
            scch(2)  = 'hour'
            iscc(2)  = 4
            scch(3)  = 'day'
            iscc(3)  = 3
            scch(4)  = 'year'
            iscc(4)  = 4

            cmch(7)  = 'nm'
            icmc(7)  = 2
            cmch(4)  = yen//'mu m'
            icmc(4)  = 5
            cmch(1)  = 'mm'
            icmc(1)  = 2
            cmch(-2) = 'm'
            icmc(-2) = 1
            cmch(-5) = 'km'
            icmc(-5) = 2

            tnch(3)  = 'psec'
            itnc(3)  = 4
            tnch(-3) = yen//'mu sec'
            itnc(-3) = 7
            tnch(-6) = 'msec'
            itnc(-6) = 4
            tnch(-9) = 'sec'
            itnc(-9) = 3

*-----------------------------------------------------------------------
*     definition of special character from integer
*-----------------------------------------------------------------------

            dsjht = ' '
            dsjhe = ' '

         do k = 1, 3
         do i = 1, 8

            jikn(k)(i:i) = char( itopi(i) )
            funn(k)(i:i) = char( itopi(i) )
            nenn(k)(i:i) = char( itopi(i) )
            tuki(k)(i:i) = char( itopi(i) )
            niti(k)(i:i) = char( itopi(i) )

         end do
         end do

         do k = 1, 3
         do i = 11, 18

            j = i - 10

            jikn(k)(i:i) = char( iboti(j) )
            funn(k)(i:i) = char( iboti(j) )
            nenn(k)(i:i) = char( iboti(j) )
            tuki(k)(i:i) = char( iboti(j) )
            niti(k)(i:i) = char( iboti(j) )

         end do
         end do

         do k = 1, 3
         do j = 1, 2

            i = 8 + j

            jikn(k)(i:i) = char( ijikn(k,j) )
            funn(k)(i:i) = char( ifunn(k,j) )
            nenn(k)(i:i) = char( inenn(k,j) )
            tuki(k)(i:i) = char( ituki(k,j) )
            niti(k)(i:i) = char( initi(k,j) )

         end do
         end do

*-----------------------------------------------------------------------
*     data and time
*-----------------------------------------------------------------------

               call date_a_time(iyer,imon,iday,ihor,imin,isec)

                     write(today1(1:2),'(i2.2)') iday

                     if( imon .eq.  1 ) today1(4:6) = 'Jan'
                     if( imon .eq.  2 ) today1(4:6) = 'Feb'
                     if( imon .eq.  3 ) today1(4:6) = 'Mar'
                     if( imon .eq.  4 ) today1(4:6) = 'Apr'
                     if( imon .eq.  5 ) today1(4:6) = 'May'
                     if( imon .eq.  6 ) today1(4:6) = 'Jun'
                     if( imon .eq.  7 ) today1(4:6) = 'Jul'
                     if( imon .eq.  8 ) today1(4:6) = 'Aug'
                     if( imon .eq.  9 ) today1(4:6) = 'Sep'
                     if( imon .eq. 10 ) today1(4:6) = 'Oct'
                     if( imon .eq. 11 ) today1(4:6) = 'Nov'
                     if( imon .eq. 12 ) today1(4:6) = 'Dec'

                     write(today1(8:11),'(i4.4)') iyer

                     today1(3:3) = '-'
                     today1(7:7) = '-'

                     write(timestr1(1:2),'(i2.2)') ihor
                     write(timestr1(4:5),'(i2.2)') imin

                     timestr1(3:3) = ':'

*-----------------------------------------------------------------------
*     initialization of random number from date
*-----------------------------------------------------------------------

         if( icang .eq. 0 ) then

               iseed1 = isec * 500 + ihor * 100 + imin
               iseed2 = imin * 500 + ihor * 100 + isec

               call rmarin(iseed1,iseed2)

         end if

*-----------------------------------------------------------------------
*     input file number constant
*-----------------------------------------------------------------------

            jsn = 0

*-----------------------------------------------------------------------

      m_err = 'Something Wrong'
      ErrCha = ''
      ErrID = 'L:369/R:a_main0/F:a-main0.f'
      l_err = 0
      k_err = jsn+1

      write(avers,'(f5.2)') version
      javers = 5

*-----------------------------------------------------------------------
*     CONSTANTS
*-----------------------------------------------------------------------

*        CM : [ 28.346457 ] CM  to Point
*        DD : [ 0.24 ]      DPI to Point

*-----------------------------------------------------------------------

         cm = 28.346457
         dd = 0.24

*-----------------------------------------------------------------------
*     WRITE A4 PAPER AND MARGIN
*-----------------------------------------------------------------------

*        A4W : [ 21.00 (20.950:old) cm ] A4 wide size
*        A4H : [ 29.70 (29.628:old) cm ] A4 hight size

*        WMG : [  1.000  cm ] left and right margin
*        HMG : [  1.414  cm ] bottom and top margin

*        WCT : [  0.050  cm ] left correction of position
*        HCT : [ -0.050  cm ] bottom correction of position

*-----------------------------------------------------------------------
         a4w = 21.00  * cm
         a4h = 29.70  * cm

         wmg = 1.0    * cm
         hmg = 1.414  * cm

         wct = 0.05   * cm
         hct =-0.05   * cm

*-----------------------------------------------------------------------
*     SUFFIXES
*-----------------------------------------------------------------------

*        FTHT : [ 0.7   ] Actual Font hight relative to Font Size
*        FTSS : [ 0.65  ] Suffix Font Size relative to Parent Font Size
*        FTKS : [ 0.333 ] FURIGANA Suffix Font Size rel to Parent Font
*        FTSV : [ 0.50  ] Vector Font Size relative to Parent Font Size
*        YPLU : [ 0.434 ] Upper Suffix Y Position relative to P.F.S
*        YPLD : [-0.189 ] Lower Suffix Y Position relative to P.F.S
*        YPKU : [ 0.850 ] FURIGANA Upper Suffix Y Position rel to P.F.S
*        YPKD : [-0.450 ] FURIGANA Lower Suffix Y Position rel to P.F.S
*        YPLH : [ 0.200 ] Hat Y Position relative to Parent Font Size
*        YPLS : [ 0.150 ] Sum and Frac Y Position relative to
*                         Parent Font Size

*-----------------------------------------------------------------------

         ftht =  0.7
         ftss =  0.65
         ftks =  0.333
         ftsv =  0.5
         yplu =  0.434
         ypld = -0.189
         ypku =  0.900
         ypkd = -0.450
         yplh =  0.200
         ypls =  0.150

*-----------------------------------------------------------------------
*     BOX CONSTANTS
*-----------------------------------------------------------------------

*        BXWS =  0.25      : box small spacing
*        BXWL =  0.50      : box large spacing
*        BXSS =  0.15      : small obal box corner and shadow small size
*        BXSL =  0.25      : large obal box corner and shadow large size
*        BXDS =  0.10      : small double line distance
*        BXDL =  0.20      : large double line distance
*        BXLS =  0.028     : thin line width
*        BXLL =  0.056     : thick line width

*-----------------------------------------------------------------------

         bxws =  0.25
         bxwl =  0.50
         bxss =  0.15
         bxsl =  0.25
         bxds =  0.10
         bxdl =  0.20
         bxls =  0.028
         bxll =  0.056

*-----------------------------------------------------------------------
*     BOUNDING BOX INITIAL VALUE AND MARGIN
*-----------------------------------------------------------------------

         bxi1 =  10000.
         bxi2 = -10000.
         byi1 =  10000.
         byi2 = -10000.

************************************************************************
*                                                                      *

      if( icang .eq. 0 ) then

      write(6,'(/
     &'' ANGEL         Graphic Program Version '',a5/
     &'' Copyright (C) Koji NIITA  1993-2005. All rights reserved.''
     & /)') avers

      end if

*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
*     read input information from 5 and interpret it
*-----------------------------------------------------------------------

         if( icang. eq. 0 ) then

            read(5,'(10000a1)', iostat = ios ) (indat(i),i=1,icolm)
            if( ios .eq. -1 ) goto 48

         else

            do i = 1, ichrl

               indat(i) = ' '

            end do

            do i = 1, 200

               indat(i) = iindat(i)

            end do

         end if

*-----------------------------------------------------------------------

         call inpdat(indat,infn,iname,iexp,
     &               iangb,iver2,inpage,ifpage,
     &               iofn,ioname,ierr)

            if( ierr .ne. 0 ) then

               m_err = 'Input information is wrong'
               ErrCha = ''
               ErrID = 'L:523/R:a_main0/F:a-main0.f'
               l_err = 1
               k_err = jsn+1
               goto 999

            end if

*-----------------------------------------------------------------------
* (1) input file name
*-----------------------------------------------------------------------

            if( iname .eq. 0 ) then

               m_err = 'Input File Name Error.'
               ErrCha = ''
               ErrID = 'L:538/R:a_main0/F:a-main0.f'
               l_err = 1
               k_err = jsn+1
               goto 999

            end if

               do 57 i = 1, iname

                  dsin(jsn+1)(i:i) = infn(i)

   57          continue

                  idsi(jsn+1) = iname

*-----------------------------------------------------------------------
*     input file dsin exist ?
*-----------------------------------------------------------------------

         if( icang. eq. 0 ) then

            inquire( file = dsin(jsn+1), exist = exex )

            if( exex .eqv. .false. ) then
               m_err = 'File Name Error. (File does not exist)'
               ErrCha = ''
               ErrID = 'L:564/R:a_main0/F:a-main0.f'
               l_err = 1
               k_err = jsn+1
               goto 999
            end if

*-----------------------------------------------------------------------
*        open first input file for alone
*        or 31 is present file from PHITS
*-----------------------------------------------------------------------

            call openf(jsi,jsn,dsin)

         else

               jsn = jsn + 1
               jsi = jsn + 30

         end if

*-----------------------------------------------------------------------
* (2) read control number for angelbat for unix
*     and open 'angelbat' file for shell script

*     [-a]

*        iangb = 0 ; default
*        iangb = 1 ; make angelbat

*-----------------------------------------------------------------------

            if( iangb .eq. 1 ) then

               open(jab,file='angelbat',status='unknown')

            end if

*-----------------------------------------------------------------------
* (3) for angel version 1.70
*
*     [-v12]
*
*        iver2 = 0 ; default
*        iver2 = 1 ; ver1 to ver2

*-----------------------------------------------------------------------
* (4,5) initial and final page number
*
*     [-f 3 -t 5] ; -f first page, -t terminate page
*
*-----------------------------------------------------------------------

            if( inpage .lt. 0 ) inpage = 1
            if( ifpage .lt. 0 ) ifpage = 100000

*-----------------------------------------------------------------------
* (6) output file name
*
*     [-o filename]
*
*-----------------------------------------------------------------------
* (7) expand input file
*
*     [-exp]
*
*        iexp = 0 ; default
*        iexp = 1 ; expand input file on ****.new
*-----------------------------------------------------------------------

            if( ioname .eq. 0 ) then

               iofile = 0

            else

               iofile = 1

            end if

*-----------------------------------------------------------------------
*     END OF READ INPUT.ANNG
*-----------------------------------------------------------------------

            goto 49

   48             m_err = 'Something wrong in input file input.ang '//
     &                    'or the number of input data mismatches.'
                  ErrCha = ''
                  ErrID = 'L:652/R:a_main0/F:a-main0.f'
                  l_err = 0
                  k_err = jsn
                  goto 999

   49       continue

*-----------------------------------------------------------------------
*     DETERMINE THE OUTPUT FILE NAME : DEFAULT IS *****.eps
*     OPEN FILE FOR 'eps' FILE
*-----------------------------------------------------------------------

               do 52 i = iname , 1, -1
                  if( infn(i) .eq. '.' ) goto 53
   52          continue

               jfname = iname
               goto 56

   53          if( i .eq. 1 ) then
                  jfname = iname
                  goto 56
               end if

               jfname = i - 1

   56          continue

*-----------------------------------------------------------------------

         if( iofile .eq. 0 ) then

            do 58 i = 1, jfname

               dsjht(i:i) = infn(i)

   58       continue

               i = jfname + 1

               dsjht(i:i+3) = '.eps'

               ioname = jfname + 4

         else if( iofile .eq. 1 ) then

            do 59 i = 1, ioname

               dsjht(i:i) = iofn(i)

   59       continue

         end if


            open( jht, file = dsjht, status = 'unknown' )


*-----------------------------------------------------------------------
*     expand input file on jht
*-----------------------------------------------------------------------

      if( iexp .ne. 0 ) then

            do i = 1, jfname

               dsjhe(i:i) = infn(i)

            end do

               i = jfname + 1

               dsjhe(i:i+3) = '.new'

               iename = jfname + 4

               open( jhs, file = dsjhe, status = 'unknown' )

            call inpexp(jsi,jsn,dsin,idsi,ierr)

               close(jhs)

               if( ierr .ne. 0 ) goto 999

               call closef(jsi,jsn)

               do i = 1, iename

                  dsin(jsn+1)(i:i) = dsjhe(i:i)

               end do

                  idsi(jsn+1) = iename

               call openf(jsi,jsn,dsin)

      end if

*-----------------------------------------------------------------------
*     LANGAGE CONTROL NUMBER
*-----------------------------------------------------------------------

*        IFON = 0 ; ENGLISH ONLY
*        IFON > 0 ; JAPANESE AND ENGLISH

*        IFON = 1 ; Shift JIS
*        IFON = 2 ; Normal JIS
*        IFON = 3 ; EUC JIS

*-----------------------------------------------------------------------
*        INITIAL VALUE OF JAPANESE KANJI CODE
*-----------------------------------------------------------------------

            ifon = -1

*-----------------------------------------------------------------------
*     File Name Convert ; \ -> \\, ( -> \(, ) -> \)
*-----------------------------------------------------------------------

               jnm = 0

               do 620 j = 1, iname

                     jnm = jnm + 1

                  if( infn(j) .eq. yen .or.
     &                infn(j) .eq. '(' .or.
     &                infn(j) .eq. ')' ) then

                     infnm(jnm)   = yen
                     infnm(jnm+1) = infn(j)
                     jnm = jnm + 1

                  else

                     infnm(jnm) = infn(j)

                  end if

  620          continue


                     call jpncode(infnm,jnm,ifon)
                     call jpnprep(infnm,jnm,ifon)


*-----------------------------------------------------------------------
*     OPEN WORK FILES
*-----------------------------------------------------------------------

                  open(jol,form='unformatted',status='scratch')
                  open(jil,form='unformatted',status='scratch')
                  open(jhc,form='unformatted',status='scratch')
                  open(jhl,form='unformatted',status='scratch')
                  open(jhy,form='unformatted',status='scratch')

                  open(jhb,form='unformatted',status='scratch')
                  open(jhd,form='unformatted',status='scratch')
                  open(jhp,form='unformatted',status='scratch')
                  open(jhq,form='unformatted',status='scratch')
                  open(jhr,form='unformatted',status='scratch')


************************************************************************

*-----------------------------------------------------------------------
*     INITIAL VALUE FOR INCLUDE FILE
*-----------------------------------------------------------------------

         do 30 i = 0, 9

            ill(i) = 0
            ilf(i) = 10000000

   30    continue


*-----------------------------------------------------------------------
*     FOR ANGEL VERSION 1.70

*        IVER2 = 0 ; DEFAULT
*        IVER2 = 1 ; VER1 TO VER2

*-----------------------------------------------------------------------

         if( iver2 .eq. 1 ) then

                call ver1to2(jsi,infn,jfname,jsn)

         else

             read(jsi,'(10000a1)',iostat=ios) ( dum(ic), ic = 1, icolm )
             if( ios .eq. -1 ) goto 141

               call chlow(dum,lum)

            if( lum(1) .eq. 'v' .and. lum(2) .eq. 'e' .and.
     &          lum(3) .eq. 'r' .and.
     &        ( lum(4) .eq. '1' .or. lum(4) .eq. '0' ) ) then


                call ver1to2(jsi,infn,jfname,jsn)

            else

                rewind(jsi)

            end if

         end if

*-----------------------------------------------------------------------
*     SOME CONSTANTS BEFORE DO LOOP OF MULTI PAGE
*-----------------------------------------------------------------------

            ierr  = 0

*-----------------------------------------------------------------------

            ibsw = 0
            idbg = 0

*-----------------------------------------------------------------------

            ipdo  = 0
            jpage = 1
            lpage = 1

            kpg = 0
            jpn = 0

            smin(1) = 3.0
            smin(2) = 1.0
            smin(3) = 1.0
            smax(1) = 1.0
            smax(2) = 1.0
            smax(3) = 1.0

************************************************************************
*                                                                      *
*     DO LOOP FOR MULTI PAGES                                          *
*                                                                      *
************************************************************************

 5000    continue

*-----------------------------------------------------------------------

               ipdo  = ipdo + 1

            if( kpg .eq. 1 ) then

               jpage = jpage + 1
               lpage = lpage + 1
               kpg   = 0

            end if

*-----------------------------------------------------------------------
*        END OF THE PROCEDURE
*-----------------------------------------------------------------------


            if( ipdo .gt. jpage ) goto 1000


*-----------------------------------------------------------------------
*        SKIP PAGE ACCORDING TO INITIAL PAGE CONTROL
*-----------------------------------------------------------------------

            if( jpage .lt. inpage ) then

               jpn = 3

            else if( jpage .gt. ifpage ) then

               jpn = 4

            end if


************************************************************************
*                                                                      *
*     DEFAULT VALUE                                                    *
*                                                                      *
*         SCAL    = 1.0                                                *
*                                                                      *
*         XORG    = 0.0                                                *
*         YORG    = 0.0                                                *
*                                                                      *
*         XFAC    = 1.0                                                *
*                                                                      *
*         IBSW    = 0                                                  *
*                                                                      *
*         IROT    = 1                                                  *
*         NOMS    = 0      ; NO MESSAGES FILENAME AND DATE             *
*         NOFR    = 0      ; NO FRAME                                  *
*         FRAM    = 0      ; WRITE FRAME                               *
*         FORM    = 0.75                                               *
*         AFAC    = 1.0                                                *
*         I()LOG  = 0      ; LINEAR IS DEFAULT                         *
*         ISPAC   = 0      ; SPACE                                     *
*         I()TIC  = 0      ; 0-> NORMAL, 1-> LARGE, -1-> SMALL TIC     *
*         ITIC    = 1      ; 1-> INSIDE, -1-> OUTSIDE                  *
*         CLUS    = 1.0    ; STORENGTH OF CLUSTER PLOT                 *
*         NOLG    = 1      ; WRITE LEGEND IF EXIST                     *
*                            0-> NO LEGEND IF DEFAULT POSITION         *
*         NOTL    = 1      ; 0-> NO TITLE IF EXIST                     *
*         NOCM    = 1      ; 0-> NO COMMENTS IF EXIST                  *
*         NOXT    = 1      ; 0-> NO X-AXIS TEXT                        *
*         NOXN    = 1      ; 0-> NO X-AXIS NUMBER                      *
*         NOYT    = 1      ; 0-> NO Y-AXIS TEXT                        *
*         NOYN    = 1      ; 0-> NO Y-AXIS NUMBER                      *
*                                                                      *
************************************************************************
*-----------------------------------------------------------------------
*     INITIAL VALUSE OF LANGAGE
*-----------------------------------------------------------------------

         if( lpage .gt. 1 ) then

            ifon = -1

         end if

            idat = -1

*-----------------------------------------------------------------------
*     BASIC FONT
*-----------------------------------------------------------------------

            ibfon = -1

*-----------------------------------------------------------------------

      ibgc  =  0

      isec  =  0

*-----------------------------------------------------------------------

      noms  =  0
      nofr  =  0

      irot  =  1

      ipsd  =  0

*-----------------------------------------------------------------------
*     MASSAGE FLAG AND DEFAULT MASSAGE

*           1: UPPER LEFT  2: UPPER CENTER  3: UPPER RIGHT
*           4: DOWN  LEFT  5: DOWN  CENTER  6: DOWN  RIGHT

*-----------------------------------------------------------------------

      imcm(1) = -2
      imcm(2) = -1
      imcm(3) = -2
      imcm(4) = -1
      imcm(5) = -1
      imcm(6) = -1

*-----------------------------------------------------------------------
*  INITIALIZATION OF THE CONTROL NUMBER
*-----------------------------------------------------------------------

      idhl  = 0
      idhc  = 0

      xfac  =  1.0
      form  =  0.75
      afac  =  1.0

      ispac =  0
      ispax =  0
      ispay =  0

      icut  =  0
      ncut  =  8
      icct  =  0
      iccr  =  1

      clus  =  1.0
      ipd2  =  0
      ipds  =  0
      ipdc  =  0

      ixtic  = 0
      iytic  = 0
      jxtic  = 2
      jytic  = 2
      itic   = 1

      cmax  = -r1max
      cmin  =  r1max

      dmax  = -r1max
      dmin  =  r1max

      do 130 i = 1, mc
         rcut(i,1) =  0.0
         rcut(i,2) = -2.0
         rcut(i,3) =  1.0
         rcut(i,4) =  1.0
         rcut(i,5) =  0.0
         iwd2(i) = 4
  130 continue

      clmo = -r1max

      clal(1) = -r1max
      clax(1) = -r1max
      cltx(1) = -r1max
      clnm(1) = -r1max
      cltl(1) = -r1max
      cllg(1) = -r1max
      clpa(1) = -r1max
      clcn(1) = -r1max
      clfr(1) = -r1max
      clms(1) = -r1max
      clgb(1) = -r1max
      clgl(1) = -r1max
      clgs(1) = -r1max
      clpb(1) = -r1max
      clpl(1) = -r1max
      clps(1) = -r1max
      clnb(1) = -r1max
      clnl(1) = -r1max
      clns(1) = -r1max
      clhd(1) = -r1max
      clbg(1) = -r1max
      clin(1) = -r1max

      clal(2) = 1.0
      clax(2) = 1.0
      cltx(2) = 1.0
      clnm(2) = 1.0
      cltl(2) = 1.0
      cllg(2) = 1.0
      clpa(2) = 1.0
      clcn(2) = 1.0
      clfr(2) = 1.0
      clms(2) = 1.0
      clgb(2) = 1.0
      clgl(2) = 1.0
      clgs(2) = 1.0
      clpb(2) = 1.0
      clpl(2) = 1.0
      clps(2) = 1.0
      clnb(2) = 1.0
      clnl(2) = 1.0
      clns(2) = 1.0
      clhd(2) = 1.0
      clbg(2) = 1.0
      clin(2) = 1.0

      clal(3) = 1.0
      clax(3) = 1.0
      cltx(3) = 1.0
      clnm(3) = 1.0
      cltl(3) = 1.0
      cllg(3) = 1.0
      clpa(3) = 1.0
      clcn(3) = 1.0
      clfr(3) = 1.0
      clms(3) = 1.0
      clgb(3) = 1.0
      clgl(3) = 1.0
      clgs(3) = 1.0
      clpb(3) = 1.0
      clpl(3) = 1.0
      clps(3) = 1.0
      clnb(3) = 1.0
      clnl(3) = 1.0
      clns(3) = 1.0
      clhd(3) = 1.0
      clbg(3) = 1.0
      clin(3) = 1.0

      ilbox = 0
      ipbox = 0
      isbox = 0

      ih2fs = 0

      icmyy  = 0
      itnyy  = 0
      isccg  = 0

      icmcg  = 0
      itncg  = 0

      xmul  = 0.0d0
      ymul  = 0.0d0
      zmul  = 0.0d0

      cmap=""
      ndis=0
      icrev=0

*-----------------------------------------------------------------------
*     INITIALIZATION OF THE FIRST DEFINITION OF OPERATORS
*-----------------------------------------------------------------------

      do 132 i = 1, 9

         idsm(i) = 0

  132 continue

      do 133 i = 1, 6

         idln(i) = 0

  133 continue

         idxe = 0
         idye = 0

         ibd1  = 0
         ibd2  = 0
         ibd3  = 0
         ibd4  = 0
         ibd5  = 0
         ibd6  = 0
         ibd7  = 0
         ibd8  = 0
         ibd9  = 0
         ibd10 = 0
         ibd11 = 0

*-----------------------------------------------------------------------
*     constants
*-----------------------------------------------------------------------

      do 81 i = 1, mxcval

         cval(i) = -r1max
         aval(i) = -r1max

   81 continue


************************************************************************
*        MULTI-GRAPH IN ONE PAGE
************************************************************************

         izdo = 0
         izcc = 1

 9010    izdo = izdo + 1

         if( izdo .gt. izcc ) goto 9000

         if( izdo .gt. 1 ) then

            xorg = xorgn
            yorg = yorgn
            angz = angzn
            noxt = noxtn
            noxn = noxnn
            noyt = noytn
            noyn = noynn

         else

            xorg = 0.0
            yorg = 0.0
            angz = -r1max
            noxt = 1
            noxn = 1
            noyt = 1
            noyn = 1

            angzp  = -r1max
            iangzp = 0

         end if

************************************************************************

*-----------------------------------------------------------------------
*     INITIALIZATION OF THE FIRST DEFINITION OF OPERATORS FOR FONTS
*-----------------------------------------------------------------------

      do 131 i = 0, 13

         ifd(i) = 0

  131 continue

         ikan1 = 0
         ikan2 = 0
         ikan3 = 0
         ikan4 = 0

*-----------------------------------------------------------------------
*     BASIC TIC FONT
*-----------------------------------------------------------------------

         itfon = -1

************************************************************************

      scal  =  1.0

      atxs  =  1.0

      iylog =  0
      ixlog =  0
      izlog =  0

      ixnum = 1
      iynum = 1
      ixtxt = 1
      iytxt = 1

      xtxp  = 0.5
      ytxp  = 0.5
      tlxp  = 0.5

      nolg  = 1
      nopa  = 0
      nocn  = 0
      notl  = 1
      nocm  = 1
      notc  = 1
      notn  = 1
      notf  = 1

      ixstd = 0
      iystd = 0
      ixltd = 0
      iyltd = 0

      noaf  = 1

      mdecx = -10
      mdecy = -10

      xstv  = 0.0
      xltv  = 0.0
      ystv  = 0.0
      yltv  = 0.0

      noned = 0
      noner = 0
      ncom  = 0

      ibmap = 0

      nhc   = 0

      do i = 1, 100

         mhc(i)  = 0
         dhgh(i) = 0.0
         dwih(i) = 0.0

      end do

      nhl   = 0
      iw    = 0
      iv    = 0
      ia    = 0
      ib    = 0
      ibc   = 0
      ipc   = 0
      isc   = 0
      irc   = 0
      ioc   = 0
      ihtl  = 0
      ihtt  = 0


      ix    = 0
      iy    = 0
      ixtt  = 0
      iytt  = 0

      ilbox = 0
      ipbox = 0
      isbox = 0




      xmul  = 0.0d0
      ymul  = 0.0d0
      zmul  = 0.0d0

      ymax  = -r0max
      ymin  =  r0max
      xmax  = -r0max
      xmin  =  r0max

      ymin2 =  r0max
      xmin2 =  r0max

      xpmin =  r1max
      xpmax = -r1max
      ypmin =  r1max
      ypmax = -r1max

      regx  = -r1max
      regy  = -r1max
      regs  =  1.0
      regd  =  1.0

      parx  = -r1max
      pary  = -r1max
      pars  =  1.0

      conx  = -r1max
      cony  = -r1max
      cons  =  1.0

      rlptl = -r1max
      sybw  = -r1max

*----------------------------------------------------------------------*
cKN 2024/01/24
      erwd  = 1.0
      ixexp = 1
      iyexp = 1
*----------------------------------------------------------------------*

      inps  = 0

      do 80 i = 1, ipsm

         xyps(i,1) = -r1max
         xyps(i,2) = -r1max
         xyps(i,3) = -r1max
         xyps(i,4) = -r1max
         xyps(i,5) = -r1max
         xyps(i,6) = -r1max

         ixps(i) = 1
         iyps(i) = 1

   80 continue

      ipdc  =  0


************************************************************************
*                                                                      *
*     START OF READ DATA AND MULTI RECORD CONTROL                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
*     READ ONE LINE from JSI
*-----------------------------------------------------------------------

  140       ill(jsn) = ill(jsn) + 1

*-----------------------------------------------------------------------

               read(jsi,'(10000a1)', iostat = ios ) (dum(ic),ic=1,icolm)
cFURUTA20131227 gfortran gives ios .ne. -1 for the second read after EOF
               if( ios .ne. 0 ) goto 251

               call chlow(dum,lum)

*-----------------------------------------------------------------------

  240       continue

               if( ill(jsn) .gt. ilf(jsn) ) goto 251

*-----------------------------------------------------------------------

               goto 250

  251          continue

               if( jpn .eq. 0 ) then

                     if( jsn .le. 1 ) then

                        goto 141

                     else if( jsn .gt. 1 ) then

                        call closef(jsi,jsn)

                        goto 140

                     end if

               else if( jpn .ne. 0 ) then

                     if( jsn .le. 1 ) then

                        goto 1000

                     else

                        call closef(jsi,jsn)

                        goto 140

                     end if

               end if

  250          continue

*-----------------------------------------------------------------------
*     SKIP THE BLANK
*-----------------------------------------------------------------------

               k = 1

            do 145 l = 1, icolm

               if(lum(l) .ne. ' ' .and. lum(l) .ne. tub ) goto 146

  145       continue

               goto 140

  146          k = l


*-----------------------------------------------------------------------
*     INCLFL: INCLUDE FILE ONLY FOR JPN = 0
*-----------------------------------------------------------------------

            if( jpn .ne. 4 .and.
     &          lum(k)//lum(k+1)//lum(k+2)//lum(k+3)//lum(k+4) .eq.
     &          'infl:' ) then

                  k = k + 4

               call inclf(jsn,jsi,dsin,dum,k,icolm,ill,ilf,idsi,ierr)

                  if( ierr .ne. 0 ) goto 999

               goto 140

            end if

*-----------------------------------------------------------------------
*     Q: STOP READING ALL THE FILE
*-----------------------------------------------------------------------

            if( lum(k)//lum(k+1) .eq. 'q:' ) then

               if( jpn .eq. 0 ) then

                  jpn = 2

                  goto 141

               else

                  goto 1000

               end if

            end if

*-----------------------------------------------------------------------
*     QP: STOP READING THE PART OF THIS PAGE UNTIL NEWPAGE:
*-----------------------------------------------------------------------

            if( jpn .eq. 0 .and.
     &          lum(k)//lum(k+1)//lum(k+2) .eq. 'qp:' ) then

                  kpg = 1
                  jpn = 1

                  goto 141

            end if


*-----------------------------------------------------------------------
*     SKIPPAGE: SKIP THE PART UNTIL NEWPAGE:
*-----------------------------------------------------------------------

            if( jpn .eq. 0 .and.
     &          lum(k  )//lum(k+1)//lum(k+2)//lum(k+3)//lum(k+4)//
     &          lum(k+5)//lum(k+6)//lum(k+7)//lum(k+8) .eq.
     &          'skippage:' ) then

                  kpg = 1
                  jpn = 1

                  goto 141

            end if

*-----------------------------------------------------------------------
*     NEWPAGE: START THE NEWPAGE
*-----------------------------------------------------------------------

            if( lum(k  )//lum(k+1)//lum(k+2)//lum(k+3)//lum(k+4)//
     &          lum(k+5)//lum(k+6)//lum(k+7) .eq.
     &          'newpage:' ) then

               if( jpn .eq. 0 ) then

                  kpg = 1
                  jpn = 0

                  goto 141

               else if( jpn .eq. 1 ) then

                  kpg = 0
                  jpn = 0

                  goto 140

               else if( jpn .eq. 3 .or. jpn .eq. 4 ) then

                  kpg = 0
                  jpn = 0

                  jpage = jpage + 1

                  goto 5000

               end if



            end if


************************************************************************
*     WITH SKIP CASE ( Q: QP: SKIPPAGE: )
************************************************************************

      if( jpn .ne. 0 ) goto 140

************************************************************************

*-----------------------------------------------------------------------
*     Z: MULTI-GRAPH
*-----------------------------------------------------------------------

            if( lum(k)//lum(k+1) .eq. 'z:' ) then

               call zmult(lum,k+2,xorgn,yorgn,angzn,
     &                    noxtn,noxnn,noytn,noynn,ierr)


               if(ierr.ne.0) then
                  l_err = ill(jsn)
                  k_err = jsn
                  goto 999
               end if

               izcc = izcc + 1

               if( izcc .gt. mc ) then
                  m_err =
     &            'Number of Records is Larger Than mc Records (MAX).'
                  ErrCha = ''
                  ErrID = 'L:1617/R:a_main0/F:a-main0.f'
                  l_err = ill(jsn)
                  k_err = jsn
                  goto 999
               end if

               goto 141

            end if

*-----------------------------------------------------------------------
*           RESET PARAMETERS OF SECTION PAPER VERSION : SECP, SECL
*-----------------------------------------------------------------------

               if( isec .ne. 0 .and. izcc .gt. 1 ) then

                  ispac = 0
                  ispax = 0
                  ispay = 0
                  xfac  = 1.0
                  form  = 0.75
                  itic  = 1
                  afac  = 1.0
                  ixltd = 0
                  iyltd = 0
                  ixstd = 0
                  iystd = 0

                  isec = 0

               end if


************************************************************************
*                                                                      *
* ### READ X:,Y:,H:,H2:,HC:,W:,P:,AW:,A: LINE AND READ DATA            *
*          Q: END OF DATA                                              *
*                                                                      *
*     IX, IY, IHTL : NUMBER OF LINE FOR X:,Y:' '                       *
*                                                                      *
*     XTT(ICHRL), YTT(ICHRL), TITLE OF AXIS                            *
*                                                                      *
*     HTT(ICHRL) ; TITLE OF FIGURE                                     *
*                                                                      *
*        IXTT,IYTT,IHTT : LENGTH OF TITLES                             *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
*     TITLE OF X
*-----------------------------------------------------------------------

         if( lum(k)//lum(k+1) .eq. 'x:' ) then

                  ix = ix + 1

                  call atitle(ill,jsi,jsn,dum,lum,k+2,xtt,ixtt,
     &                        ifon,iend)

                  if( ixtt .eq. 0 ) then

                     xtt(1) = ' '
                     ixtt   =  1

                  end if

                  if( iend .ne. 0 ) then

                     if( jpn .eq. 0 ) then

                        goto 141

                     else if( jpn .ne. 0 ) then

                        goto 1000

                     end if

                  end if

               goto 140

         end if

*-----------------------------------------------------------------------
*     TITLE OF Y
*-----------------------------------------------------------------------

         if( lum(k)//lum(k+1) .eq. 'y:' ) then

                  iy = iy + 1

                  call atitle(ill,jsi,jsn,dum,lum,k+2,ytt,iytt,
     &                        ifon,iend)

                  if( iytt .eq. 0 ) then

                     ytt(1) = ' '
                     iytt   =  1

                  end if

                  if( iend .ne. 0 ) then

                     if( jpn .eq. 0 ) then

                        goto 141

                     else if( jpn .ne. 0 ) then

                        goto 1000

                     end if

                  end if

               goto 140

         end if

*-----------------------------------------------------------------------
*     TITLE OF GRAPH
*-----------------------------------------------------------------------

         if( lum(k) .eq. "'" ) THEN

                  call atitle(ill,jsi,jsn,dum,lum,k+1,htt,ihtt,
     &                        ifon,iend)

                     if( ihtt .ne. 0 ) then

                        if( htt(ihtt) .ne. "'" ) THEN

                           ihtt = 0

                        else

                           ihtt = ihtt - 1

                        end if

                     end if

                  if( iend .ne. 0 ) then

                     if( jpn .eq. 0 ) then

                        goto 141

                     else if( jpn .ne. 0 ) then

                        goto 1000

                     end if

                  end if

               goto 140

         end if

*-----------------------------------------------------------------------
*     MASSAGE OF THE CORNER
*-----------------------------------------------------------------------

         if( izdo .eq. 1 .and.
     &       lum(k)  .eq.'m' .and. lum(k+1).eq.'s' .and.
     &      (lum(k+2).eq.'u' .or.  lum(k+2).eq.'d').and.
     &      (lum(k+3).eq.'l' .or.  lum(k+3).eq.'c'.or.
     &       lum(k+3).eq.'r').and.
     &       lum(k+4).eq.':') then

                  call atitle(ill,jsi,jsn,dum,lum,k+5,ccmt,icmm,
     &                        ifon,iend)


            if( lum(k+2)//lum(k+3) .eq.'ul' ) then

                  imcp = 1
                  imcm(imcp) = icmm

               do 601 i = 1, icmm

                  mcm(imcp,i) = ccmt(i)

  601          continue

            else if( lum(k+2)//lum(k+3) .eq.'uc' ) then

                  imcp = 2
                  imcm(imcp) = icmm

               do 602 i = 1, icmm

                  mcm(imcp,i) = ccmt(i)

  602          continue

            else if( lum(k+2)//lum(k+3) .eq.'ur' ) then

                  imcp = 3
                  imcm(imcp) = icmm

               do 603 i = 1, icmm

                  mcm(imcp,i) = ccmt(i)

  603          continue

            else if( lum(k+2)//lum(k+3) .eq.'dl' ) then

                  imcp = 4
                  imcm(imcp) = icmm

               do 604 i = 1, icmm

                  mcm(imcp,i) = ccmt(i)

  604          continue

            else if( lum(k+2)//lum(k+3) .eq.'dc' ) then

                  imcp = 5
                  imcm(imcp) = icmm

               do 605 i = 1, icmm

                  mcm(imcp,i) = ccmt(i)

  605          continue

            else if( lum(k+2)//lum(k+3) .eq.'dr' ) then

                  imcp = 6
                  imcm(imcp) = icmm

               do 606 i = 1, icmm

                  mcm(imcp,i) = ccmt(i)

  606          continue

            end if


                  if( iend .ne. 0 ) then

                     if( jpn .eq. 0 ) then

                        goto 141

                     else if( jpn .ne. 0 ) then

                        goto 1000

                     end if

                  end if

               goto 140

         end if


*-----------------------------------------------------------------------
*     LINE PLOT H:
*-----------------------------------------------------------------------

         if( lum(k)//lum(k+1) .eq. 'h:' ) then

*-----------------------------------------------------------------------
*           MULTI LINE DESCRIPTION
*-----------------------------------------------------------------------

               in = k + 2

               call seqlin(dum,in,jsi,jsn,ill,ilf,ierr)
               call chlow(dum,lum)

               if( ierr .ne. 0 ) goto 999

*-----------------------------------------------------------------------

            call hone(dum,lum,k+2,jsi,jsn,ill,ilf,dsin,idsi,
     &                ymax,ymin,ymin2,
     &                xmax,xmin,xmin2,noned,noner,ncom,ierr,jol,jil,
     &                iycm,rycm,ifon)

               if( ierr .ne. 0 ) goto 999

               goto 240

         end if

*-----------------------------------------------------------------------
*     3-D CLUSTER PLOT AND CONTOUR PLOT
*-----------------------------------------------------------------------

         if( ( lum(k) .eq. 'h' .and. lum(k+2) .eq. ':' .and.
     &       ( lum(k+1) .eq. 'c' .or.
     &         lum(k+1) .eq. '2' .or.
     &         lum(k+1) .eq. 'd') ) .or.
     &       ( lum(k) .eq. 'h' .and. lum(k+3) .eq. ':' .and.
     &     ( ( lum(k+1) .eq. 'c' .and. lum(k+2) .eq. '2' ) .or.
     &       ( lum(k+1) .eq. 'd' .and. lum(k+2) .eq. '2' ) ) ) ) then

                  nh2  = 0

                  if( lum(k+1) .eq. 'd' .and. lum(k+2) .eq. ':' )
     &                                    mhc2 = 1
                  if( lum(k+1) .eq. '2' ) mhc2 = 2
                  if( lum(k+1) .eq. 'c' .and. lum(k+2) .eq. ':' )
     &                                    mhc2 = 3
                  if( lum(k+1) .eq. 'd' .and. lum(k+2) .eq. '2' )
     &                                    mhc2 = 4
                  if( lum(k+1) .eq. 'c' .and. lum(k+2) .eq. '2' )
     &                                    mhc2 = 5

               call prehtwo(lum,k+3,jsn,ill,
     &                      xfin,yfin,xint,yint,ixy,ic,
     &                      ierr)
                  if( ierr .ne. 0 ) goto 999

               allocate(daxy(ixnm,iynm),dax(ixnm),day(iynm))

               call htwo(dum,lum,k+3,jsi,jsn,ill,ilf,dsin,idsi,
     &                   ymax,ymin,xmax,xmin,
     &                   nhl,nh2,ierr,mhc2,ipds,ipdc,clhd,
     &                   izlog,daxy,dax,day,
     &                   xfin,yfin,xint,yint,ixy,ic,
     &                   ih2fs)

                  if( ierr .ne. 0 ) goto 999

*-----------------------------------------------------------------------
*              WRITE CONTOUR ON JOL
*-----------------------------------------------------------------------

               if( nh2 .gt. 0 ) then

                  call contjol(jol,icut,ncut,rcut,icct,iccr,
     &                         iwd2,ipd2,noned,izlog,ibmap,
     &                         daxy,ixnm,iynm,dax,day,ierr)
                  if( ierr .ne. 0 ) goto 999

               end if

*-----------------------------------------------------------------------
*           END OF CONTOUR PLOT
*-----------------------------------------------------------------------

               deallocate(daxy,dax,day)

               goto 240

         end if

*-----------------------------------------------------------------------
*     bitmap with clipping or without
*-----------------------------------------------------------------------

         if( lum(k) .eq. 'h' .and. lum(k+1) .eq. 'b' .and.
     &       lum(k+2) .eq. ':' ) then

               call bmap(dum,lum,k+3,jsi,jsn,ill,ilf,dsin,idsi,
     &                   ibmap,ymax,ymin,xmax,xmin,ierr)

                  if( ierr .ne. 0 ) goto 999

               goto 240

         end if

*-----------------------------------------------------------------------
*     INPS: INCLUDE PS FILE
*-----------------------------------------------------------------------

         if( lum(k)//lum(k+1)//lum(k+2)//lum(k+3)//lum(k+4) .eq.
     &       'inps:' ) then

                  k = k + 4

               call inclp(inps,dpin,idpi,ill,jsn,dum,lum,k,ierr,
     &                    psbbx,xyps,ixps,iyps)

                  if( ierr .ne. 0 ) goto 999

               goto 140

         end if

*-----------------------------------------------------------------------
*     SET: SET CONSTANTS C1, C2, ....
*-----------------------------------------------------------------------

         if( lum(k)//lum(k+1)//lum(k+2)//lum(k+3) .eq.
     &       'set:' ) then

                  k = k + 3

               call setcv(lum,k,icolm,dsin,idsi,ill,jsn,ierr)
                  if( ierr .ne. 0 ) goto 999

               goto 140

         end if

*-----------------------------------------------------------------------
*     PARAMETERS
*-----------------------------------------------------------------------

         if(lum(k)//lum(k+1) .eq. 'p:' ) then

            call para(dum,lum,k+2,form,xpmin,xpmax,ypmin,ypmax,
     &                afac,iylog,ixlog,atxs,
     &                ispac,ispax,ispay,icut,rcut,ncut,
     &                ixtic,iytic,nocm,notc,notn,notf,
     &                itic,iwd2,ipd2,ipds,ipdc,
     &                scal,noms,nolg,notl,ibsw,
     &                jxtic,jytic,iccr,
     &                regx,regy,regs,regd,irot,xorg,yorg,xfac,nofr,
     &                noxt,noxn,noyt,noyn,ierr,izdo,icct,ibfon,itfon,
     &                clal,clax,cltx,clnm,cltl,cllg,clfr,clms,clhd,
     &                clbg,clin,clmo,izlog,isec,noaf,angz,
     &                ixstd,iystd,ixltd,iyltd,xsd,xld,ysd,yld,
     &                mdecx,mdecy,xstv,xltv,ystv,yltv,rlptl,sybw,idbg,
     &                ixnum,iynum,ixtxt,iytxt,xtxp,ytxp,tlxp,ifon,idat,
     &                clgb,clgl,clgs,ilbox,
     &                nopa,parx,pary,pars,clpa,ipbox,clpb,clpl,clps,
     &                nocn,conx,cony,cons,clcn,isbox,clnb,clnl,clns,
     &                xmul,ymul,zmul,icmcg,itncg,icmyy,itnyy,isccg,
     &                erwd,ixexp,iyexp,
cKN 2024/01/24
     &                ih2fs)

            if(ierr.ne.0) then
               l_err = ill(jsn)
               k_err = jsn
               goto 999
            end if


               goto 140

         end if

*-----------------------------------------------------------------------
*     TABLE COMMENT
*-----------------------------------------------------------------------

         if( lum(k)//lum(k+1)//lum(k+2)//lum(k+3)//lum(k+4) .eq.
     &       'wtab:' ) then

               iv = iv + 1

                  if( iv .eq. 1 )
     &            open(jwt,form='unformatted',status='scratch')

            call wctb(dum,lum,k+5,iv,ierr,ifon,ill,ilf,jsi,jsn,iend)

                  if( ierr .ne. 0 ) goto 999

                  if( iend .ne. 0 ) then

                     if( jpn .eq. 0 ) then

                        goto 141

                     else if( jpn .ne. 0 ) then

                        goto 1000

                     end if

                  end if

               goto 140

         end if

*-----------------------------------------------------------------------
*     TEXT COMMENT
*-----------------------------------------------------------------------

         if(lum(k)//lum(k+1)//lum(k+2) .eq. 'wt:' ) then

               iw = iw + 1

               if( iw .eq. 1 )
     &         open(jhm,form='unformatted',status='scratch')

                     if( iw .gt. mc ) then
                        l_err = ill(jsn)
                        k_err = jsn
      m_err = 'Number of Comment is larger than mc. Please try output =
     & 2 or gshow = 1'
                        ErrCha = ''
                        ErrID = 'L:2114/R:a_main0/F:a-main0.f'
                        goto 999
                     end if

            call wcmt(dum,lum,k+3,iw,wxys,iwx,iwy,iwn,ierr,ifon,
     &                ill,ilf,jsi,jsn,iend,iwf,iwb,ibox,cbox)

                     if( ierr .ne. 0 ) then
                        l_err = ill(jsn)
                        k_err = jsn
                        goto 999
                     end if

                  if( iend .ne. 0 ) then

                     if( jpn .eq. 0 ) then

                        goto 141

                     else if( jpn .ne. 0 ) then

                        goto 1000

                     end if

                  end if

               goto 140

         end if

*-----------------------------------------------------------------------
*     ONE LINE COMMENT
*-----------------------------------------------------------------------

         if(lum(k)//lum(k+1) .eq. 'w:' ) then

               iw = iw + 1

               if( iw .eq. 1 )
     &         open(jhm,form='unformatted',status='scratch')

                     if( iw .gt. mc ) then
                        l_err = ill(jsn)
                        k_err = jsn
      m_err = 'Number of Comment is larger than mc. Please try output =
     & 2 or gshow = 1'
                        ErrCha = ''
                        ErrID = 'L:2162/R:a_main0/F:a-main0.f'
                        goto 999
                     end if

            call wcmm(dum,lum,k+2,iw,wxys,iwx,iwy,iwn,ierr,ifon,
     &                ill,ilf,jsi,jsn,iend,iwf,iwb)

                     if( ierr .ne. 0 ) then
                        l_err = ill(jsn)
                        k_err = jsn
                        goto 999
                     end if

                  if( iend .ne. 0 ) then

                     if( jpn .eq. 0 ) then

                        goto 141

                     else if( jpn .ne. 0 ) then

                        goto 1000

                     end if

                  end if

               goto 140

         end if

*-----------------------------------------------------------------------
*     ARROW AND COMMENT
*-----------------------------------------------------------------------

         if(lum(k)//lum(k+1)//lum(k+2) .eq. 'aw:' ) then

               ia = ia + 1

               if( ia .eq. 1 )
     &         open(jha,form='unformatted',status='scratch')

               call arrcm(dum,lum,k+3,ia,axys,iax,ian,ierr,ifon,
     &                    ill,ilf,jsi,jsn,iend,iaf)

                     if( ierr .ne. 0 ) then
                        l_err = ill(jsn)
                        k_err = jsn
                        goto 999
                     end if

                  if( iend .ne. 0 ) then

                     if( jpn .eq. 0 ) then

                        goto 141

                     else if( jpn .ne. 0 ) then

                        goto 1000

                     end if

                  end if

               goto 140

         end if

*-----------------------------------------------------------------------
*     ARROW ONLY
*-----------------------------------------------------------------------

         if(lum(k)//lum(k+1) .eq.'a:' ) then

               ib = ib + 1

               call arrow(dum,lum,k+2,ill,ilf,jsn,ib,bxys,ierr)

                  if( ierr. ne. 0 ) then
                     l_err = ill(jsn)
                     k_err = jsn
                     goto 999
                  end if

               goto 140

         end if

*-----------------------------------------------------------------------
*     ARROW BOUND ONLY
*-----------------------------------------------------------------------

         if(lum(k)//lum(k+1)//lum(k+2) .eq. 'ab:' ) then

               ibc = ibc + 1

               call arrowc(dum,lum,k+3,ill,ilf,jsn,ibc,cxys,ierr)

                  if( ierr .ne. 0 ) then
                     l_err = ill(jsn)
                     k_err = jsn
                     goto 999
                  end if

               goto 140

         end if

*-----------------------------------------------------------------------
*     POLYGON
*-----------------------------------------------------------------------

         if(lum(k)//lum(k+1)//lum(k+2)//lum(k+3)//lum(k+4) .eq.
     &      'polg:' ) then

               ipc = ipc + 1

               call polg(dum,lum,k+5,ill,ilf,jsn,ipc,pxys,ierr)

                  if( ierr .ne. 0 ) then
                     l_err = ill(jsn)
                     k_err = jsn
                     goto 999
                  end if

               goto 140

         end if

*-----------------------------------------------------------------------
*     STAR
*-----------------------------------------------------------------------

         if(lum(k)//lum(k+1)//lum(k+2)//lum(k+3)//lum(k+4) .eq.
     &      'star:' ) then

               isc = isc + 1

               call star(dum,lum,k+5,ill,ilf,jsn,isc,sxys,ierr)

                  if( ierr .ne. 0 ) then
                     l_err = ill(jsn)
                     k_err = jsn
                     goto 999
                  end if

               goto 140

         end if

*-----------------------------------------------------------------------
*     RIBBON
*-----------------------------------------------------------------------

         if(lum(k)//lum(k+1)//lum(k+2)//lum(k+3)//lum(k+4) .eq.
     &      'ribn:' ) then

               irc = irc + 1

               call ribn(dum,lum,k+5,ill,ilf,jsn,irc,rxys,ierr)

                  if( ierr .ne. 0 ) then
                     l_err = ill(jsn)
                     k_err = jsn
                     goto 999
                  end if

               goto 140

         end if

*-----------------------------------------------------------------------
*     BOX
*-----------------------------------------------------------------------

         if(lum(k)//lum(k+1)//lum(k+2)//lum(k+3) .eq.
     &      'box:' ) then

               ioc = ioc + 1

               call rbox(dum,lum,k+4,ill,ilf,jsn,ioc,oxys,ierr)

                  if( ierr .ne. 0 ) then
                     l_err = ill(jsn)
                     k_err = jsn
                     goto 999
                  end if

               goto 140

         end if

************************************************************************


         goto 140


************************************************************************
*     END OF READ SOURCE FILE
************************************************************************


  141    continue


*-----------------------------------------------------------------------

            xmuc = 1.0d0
            ymuc = 1.0d0
            xmut = 1.0d0
            ymut = 1.0d0

*-----------------------------------------------------------------------
*     sec -> min(scmn), hour(schr), day(scdy), year(scyr)
*-----------------------------------------------------------------------

       if( isccg .gt. 0 ) then

                do ii = 1, ixtt
                   xttn(ii) = xtt(ii)
                end do
          do i = 1, ixtt-4
              if( xttn(i)//xttn(i+1)//xttn(i+2)//xttn(i+3)//xttn(i+4)
     &            .eq. '[sec]' ) then
                do j = 1, iscc(isccg)
                   xtt(i+j) = scch(isccg)(j:j)
                end do
                do k = i + 4, ixtt
                   xtt(i+iscc(isccg)+k-i-3) = xttn(k)
                end do
                   ixtt = iscc(isccg) + ixtt - 3
                   if( isccg .eq. 1 ) tfact = 60.d0
                   if( isccg .eq. 2 ) tfact = 3600.d0
                   if( isccg .eq. 3 ) tfact = 3600.d0*24.0
                   if( isccg .eq. 4 ) tfact = 3600.d0*24.0*365.0
                   xmuc = 1.d0/tfact
             end if
          end do

       end if

*-----------------------------------------------------------------------

       if( isccg .eq. 0 ) then

               if( xpmin .gt. r0max ) then
                  xxmax = xmax
               else
                  xxmax = xpmax
               end if

                do ii = 1, ixtt
                   xttn(ii) = xtt(ii)
                end do
          do i = 1, ixtt-4
              if( xttn(i)//xttn(i+1)//xttn(i+2)//xttn(i+3)//xttn(i+4)
     &            .eq. '[sec]' ) then
             if( xxmax .ge. 600.0 ) then
                   if( xxmax .lt. 3600.0*10.0 ) then
                      icccg = 1
                   else if( xxmax .lt. 3600.0*24.0*10.0 ) then
                      icccg = 2
                   else if( xxmax .lt. 3600.0*24.0*365.0*10.0 ) then
                      icccg = 3
                   else
                      icccg = 4
                   end if

                do j = 1, iscc(icccg)
                   xtt(i+j) = scch(icccg)(j:j)
                end do
                do k = i + 4, ixtt
                   xtt(i+iscc(icccg)+k-i-3) = xttn(k)
                end do
                   ixtt = iscc(icccg) + ixtt - 3
                   if( icccg .eq. 1 ) tfact = 60.d0
                   if( icccg .eq. 2 ) tfact = 3600.d0
                   if( icccg .eq. 3 ) tfact = 3600.d0*24.0
                   if( icccg .eq. 4 ) tfact = 3600.d0*24.0*365.0
                   xmuc = 1.d0/tfact
             end if
             end if
          end do

       end if

*-----------------------------------------------------------------------
*     cm -> nm(cmnm), um(cmum), mm(cmmm), mt(cmmt), Km(cmkm)
*-----------------------------------------------------------------------

       if( icmcg .ne. 0 ) then

                do ii = 1, ixtt
                   xttn(ii) = xtt(ii)
                end do
          do i = 1, ixtt-3
              if( xttn(i)//xttn(i+1)//xttn(i+2)//xttn(i+3)
     &            .eq. '[cm]' ) then
                do j = 1, icmc(icmcg)
                   xtt(i+j) = cmch(icmcg)(j:j)
                end do
                do k = i + 3, ixtt
                   xtt(i+icmc(icmcg)+k-i-2) = xttn(k)
                end do
                   ixtt = icmc(icmcg) + ixtt - 2
                   xmuc = 10.d0**(icmcg)
             end if
          end do

*-----------------------------------------------------------------------

                do ii = 1, iytt
                   yttn(ii) = ytt(ii)
                end do
          do i = 1, iytt-3
              if( yttn(i)//yttn(i+1)//yttn(i+2)//yttn(i+3)
     &            .eq. '[cm]' ) then
                do j = 1, icmc(icmcg)
                   ytt(i+j) = cmch(icmcg)(j:j)
                end do
                do k = i + 3, iytt
                   ytt(i+icmc(icmcg)+k-i-2) = yttn(k)
                end do
                   iytt = icmc(icmcg) + iytt - 2
                   ymuc = 10.d0**(icmcg)
             end if
          end do

*-----------------------------------------------------------------------

       if( icmyy .ne. 0 ) then

                do ii = 1, ixtt
                   xttn(ii) = xtt(ii)
                end do
          do l = 1, ixtt-2
              if( xttn(i)//xttn(i+1)//xttn(i+2) .eq. '/cm' ) then
                   ipww = 1
                if( xttn(l+4) .eq. '2' ) ipww = 2
                if( xttn(l+4) .eq. '3' ) ipww = 3
                do j = 1, icmc(icmcg)
                   xtt(l+j) = cmch(icmcg)(j:j)
                end do
                do k = l + 3, ixtt
                   xtt(l+icmc(icmcg)+k-l-2) = xttn(k)
                end do
                   ixtt = icmc(icmcg) + ixtt - 2
                   xmuc = 1.d0/(10.d0**(icmcg))**(ipww)
             end if
          end do

*-----------------------------------------------------------------------

                do ii = 1, iytt
                   yttn(ii) = ytt(ii)
                end do
          do l = 1, iytt-2
              if( yttn(i)//yttn(i+1)//yttn(i+2) .eq. '/cm' ) then
                   ipww = 1
                if( yttn(l+4) .eq. '2' ) ipww = 2
                if( yttn(l+4) .eq. '3' ) ipww = 3
                do j = 1, icmc(icmcg)
                   ytt(l+j) = cmch(icmcg)(j:j)
                end do
                do k = l + 3, iytt
                   ytt(l+icmc(icmcg)+k-l-2) = yttn(k)
                end do
                   iytt = icmc(icmcg) + iytt - 2
                   ymuc = 1.d0/(10.d0**(icmcg))**(ipww)
             end if
          end do

*-----------------------------------------------------------------------

             do ii = 1, ixtt
                xttn(ii) = xtt(ii)
             end do
          do l = 1, ixtt-1
              if( xttn(i)//xttn(i+1) .eq. 'cm' ) then
                do j = 1, icmc(icmcg)
                   xtt(l+j-1) = cmch(icmcg)(j:j)
                end do
                do k = l + 2, ixtt
                   xtt(l+icmc(icmcg)-1+k-l-1) = xttn(k)
                end do
                   ixtt = icmc(icmcg) + ixtt - 2
                   xmuc = 10.d0**(icmcg)
             end if
          end do

*-----------------------------------------------------------------------

             do ii = 1, iytt
                yttn(ii) = ytt(ii)
             end do
          do l = 1, iytt-1
              if( yttn(i)//yttn(i+1)
     &            .eq. 'cm' ) then
                do j = 1, icmc(icmcg)
                   ytt(l+j-1) = cmch(icmcg)(j:j)
                end do
                do k = l + 2, iytt
                   ytt(l+icmc(icmcg)-1+k-l-1) = yttn(k)
                end do
                   iytt = icmc(icmcg) + iytt - 2
                   ymuc = 10.d0**(icmcg)
             end if
          end do

       end if

       end if

*-----------------------------------------------------------------------
*      nsec -> psec(nsps), usec(nsus), msec(nsms), sec(nssc)
*-----------------------------------------------------------------------

       if( itncg .ne. 0 ) then

                do ii = 1, ixtt
                   xttn(ii) = xtt(ii)
                end do
          do i = 1, ixtt-5
              if( xttn(i)//xttn(i+1)//xttn(i+2)//xttn(i+3)//xttn(i+4)
     &//xttn(i+5)
     &            .eq. '[nsec]' ) then
                do j = 1, itnc(itncg)
                   xtt(i+j) = tnch(itncg)(j:j)
                end do
                do k = i + 5, ixtt
                   xtt(i+itnc(itncg)+k-i-4) = xttn(k)
                end do
                   ixtt = itnc(itncg) + ixtt - 4
                   xmut = 10.d0**(itncg)
             end if
          end do

*-----------------------------------------------------------------------

                do ii = 1, iytt
                   yttn(ii) = ytt(ii)
                end do
          do i = 1, iytt-5
              if( yttn(i)//yttn(i+1)//yttn(i+2)//yttn(i+3)//yttn(i+4)
     &//yttn(i+5)
     &            .eq. '[nsec]' ) then
                do j = 1, itnc(itncg)
                   ytt(i+j) = tnch(itncg)(j:j)
                end do
                do k = i + 5, iytt
                   ytt(i+itnc(itncg)+k-i-4) = yttn(k)
                end do
                   iytt = itnc(itncg) + iytt - 4
                   ymut = 10.d0**(itncg)
             end if
          end do

*-----------------------------------------------------------------------

       if( itnyy .ne. 0 ) then

                do ii = 1, ixtt
                   xttn(ii) = xtt(ii)
                end do
          do l = 1, ixtt-4
              if( xttn(i)//xttn(i+1)//xttn(i+2)//xttn(i+3)//xttn(i+4)
     &            .eq. '/nsec' ) then
                do j = 1, itnc(itncg)
                   xtt(l+j) = tnch(itncg)(j:j)
                end do
                do k = l + 5, ixtt
                   xtt(l+itnc(itncg)+k-l-4) = xttn(k)
                end do
                   ixtt = itnc(itncg) + ixtt - 4
                   xmut = 1.d0/10.d0**(itncg)
             end if
          end do

*-----------------------------------------------------------------------

                do ii = 1, iytt
                   yttn(ii) = ytt(ii)
                end do
          do l = 1, iytt-4
              if( yttn(i)//yttn(i+1)//yttn(i+2)//yttn(i+3)//yttn(i+4)
     &            .eq. '/nsec' ) then
                do j = 1, itnc(itncg)
                   ytt(l+j) = tnch(itncg)(j:j)
                end do
                do k = l + 5, iytt
                   ytt(l+itnc(itncg)+k-l-4) = yttn(k)
                end do
                   iytt = itnc(itncg) + iytt - 4
                   ymut = 1.d0/10.d0**(itncg)
             end if
          end do

       end if

       end if

*-----------------------------------------------------------------------

          if( xmuc * xmut .ne. 1.0d0 ) then
             xmul = xmuc * xmut
          end if

          if( ymuc * ymut .ne. 1.0d0 ) then
             ymul = ymuc * ymut
          end if

*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*     macro parameters
*-----------------------------------------------------------------------

         if( notc .eq. 0 ) then

              itic=0
              noxn=0
              noyn=0

         end if

         if( notn .eq. 0 ) then

              itic=0
              noxn=0
              noyn=0
              noxt=0
              noyt=0
              notl=0

         end if

         if( notf .eq. 0 ) then

              itic=0
              noxn=0
              noyn=0
              noxt=0
              noyt=0
              noaf=0
              notl=0

         end if

*-----------------------------------------------------------------------
*     COLOR OF AXIS, TEXT, NUMBER, TITLE, LEGEND, COMMENT AND LINE
*           DETERMINE BY CLAL( )
*           EXCEPT FOR COLOR OF REGION, COLOR CLUSTER
*-----------------------------------------------------------------------

               if( clal(1) .gt. -r0max ) then

                 clax(1) = clal(1)
                 cltx(1) = clal(1)
                 clnm(1) = clal(1)
                 cltl(1) = clal(1)
                 cllg(1) = clal(1)
                 clpa(1) = clal(1)
                 clfr(1) = clal(1)
                 clms(1) = clal(1)

                 clax(2) = clal(2)
                 cltx(2) = clal(2)
                 clnm(2) = clal(2)
                 cltl(2) = clal(2)
                 cllg(2) = clal(2)
                 clpa(2) = clal(2)
                 clfr(2) = clal(2)
                 clms(2) = clal(2)

                 clax(3) = clal(3)
                 cltx(3) = clal(3)
                 clnm(3) = clal(3)
                 cltl(3) = clal(3)
                 cllg(3) = clal(3)
                 clpa(3) = clal(3)
                 clfr(3) = clal(3)
                 clms(3) = clal(3)

               end if

*-----------------------------------------------------------------------
*     CHANGE ALL COLOR TO MONOCHROME
*       COLOR OF AXIS, TEXT, NUMBER, TITLE, LEGEND, COMMENT AND LINE
*       AND COLOR OF REGION
*-----------------------------------------------------------------------

            if( clmo .gt. -r0max ) then

                 clax(1) = -r1max
                 cltx(1) = -r1max
                 clnm(1) = -r1max
                 cltl(1) = -r1max
                 cllg(1) = -r1max
                 clpa(1) = -r1max
                 clfr(1) = -r1max
                 clms(1) = -r1max

            end if

*-----------------------------------------------------------------------
*### X-AXIS AND TIC RECORD ###
*-----------------------------------------------------------------------

               if( ispac .eq.  1 ) ispax = 1
               if( ispac .eq.  1 ) ispay = 1

               if( xpmin .lt.  r0max ) xmin = xpmin
               if( xpmax .gt. -r0max ) xmax = xpmax

               if( xmin .gt. r9max .and. xmax .lt. -r9max ) then

                  if( ixlog .eq. 0 ) then

                     xmin = 0.0
                     xmax = 1.0

                  else

                     xmin = 1.0
                     xmax = 100.0

                  end if

               end if

               if( xmin .gt. r9max .and. xmax .ge. -r9max ) then

                  if( ixlog .eq. 0 ) then

                     xmin = xmax / 2.0

                  else

                     xmin = xmax / 100.0

                  end if

               end if

               if( xmin .le. r9max .and. xmax .lt. -r9max ) then

                  if( ixlog .eq. 0 ) then

                     xmax = xmin * 2.0

                  else

                     xmax = xmin * 100.0

                  end if

               end if

               if( xmin .eq. xmax ) then

                  if( xmin .eq. 0.0 .and. ixlog .eq. 0 ) then

                     xmin = -1.0
                     xmax =  1.0

                  else if( xmin .eq. 0.0 .and. ixlog .eq. 1 ) then

                     xmin = 0.1
                     xmax = 10.0

                  else if( ixlog .eq. 0 ) then

                     xmin = xmin * 0.5
                     xmax = xmax * 1.5

                  else

                     xmin = xmin / 10.0
                     xmax = xmax * 10.0

                  end if

               end if

               if( ixlog .eq. 1 ) then

                  if( xmin .le. 0.0 .and.
     &                xmin2 .gt. 0.0 .and.
     &                xmin2 .le. r9max ) then

                     xmin = xmin2 * 0.5

                  end if

               end if

               if( xmax .lt. xmin ) then

                  m_err = 'XMAX is smaller than XMIN, Please CHeck '//
     &                    'XMAX( ) and XMIN( )'
                  ErrCha = ''
                  ErrID = 'L:2865/R:a_main0/F:a-main0.f'
                  l_err = ill(jsn)
                  k_err = jsn
                  goto 999

               end if


               xmaxd = 30.0

               if( ixlog .eq. 1 )then
                  if( log10( xmax / xmin ) .gt. xmaxd ) then

                  xmin = 10.d0**( log10( xmax) - xmaxd )


                  end if
               endif

               call xyrmax(xmin,xmax,ixlog,ispax,ierr,
     &                     xpmin,xpmax)
cKN 2024/01/24

               if( ierr .ne. 0 ) then
                  m_err = 'X-Axis is LOG but Xmin is Negative.'
                  ErrCha = ''
                  ErrID = 'L:2891/R:a_main0/F:a-main0.f'
                  goto 999
               end if


               call xytic(ixlog,xmin,xmax,ixtic,xld,xsd,
     &                    nxlta,nxsta,xlta,xsta,xcta,nxtch,
     &                    idecx,ixstd,ixltd,mdecx,xstv,xltv,ierr,
     &                    xmul)


               if( ierr .ne. 0 ) then
                  m_err = 'Error in Subroutine XYTIC, Sorry.'
                  ErrCha = ''
                  ErrID = 'L:2905/R:a_main0/F:a-main0.f'
                  goto 999
               end if


*-----------------------------------------------------------------------
*### Y-AXIS ###
*-----------------------------------------------------------------------

               if( ypmin .lt.  r0max ) ymin = ypmin
               if( ypmax .gt. -r0max ) ymax = ypmax

               if( ymin .gt. r9max .and. ymax .lt. -r9max ) then

                  if( iylog .eq. 0 ) then

                     ymin = 0.0
                     ymax = 1.0

                  else

                     ymin = 1.0
                     ymax = 100.0

                  end if

               end if

               if( ymin .gt. r9max .and. ymax .ge. -r9max ) then

                  if( iylog .eq. 0 ) then

                     ymin = ymax / 2.0

                  else

                     ymin = ymax / 100.0

                  end if

               end if

               if( ymin .le. r9max .and. ymax .lt. -r9max ) then

                  if( iylog .eq. 0 ) then

                     ymax = ymin * 2.0

                  else

                     ymax = ymin * 100.0

                  end if

               end if

               if( ymin .eq. ymax ) then

                  if( ymin .eq. 0.0 .and. iylog .eq. 0 ) then

                     ymin = -1.0
                     ymax =  1.0

                  else if( ymin .eq. 0.0 .and. iylog .eq. 1 ) then

                     ymin = 0.1
                     ymax = 10.0

                  else if( iylog .eq. 0 ) then

                     ymin = ymin * 0.5
                     ymax = ymax * 1.5

                  else

                     ymin = ymin / 10.0
                     ymax = ymax * 10.0

                  end if

               end if

               if( iylog .eq. 1 ) then

                  if( ymin .le. 0.0 .and.
     &                ymin2 .gt. 0.0 .and.
     &                ymin2 .le. r9max ) then

                     ymin = ymin2 * 0.5

                  end if

               end if

            if( ymax .lt. ymin ) then

               if( iylog .eq. 1 .and. ymax .eq. 0.0 ) then

                     ymin = 0.1
                     ymax = 10.0

               else

                  m_err = 'YMAX is smaller than YMIN, Please CHeck '//
     &                    'YMAX( ) and YMIN( )'
                  ErrCha = ''
                  ErrID = 'L:3011/R:a_main0/F:a-main0.f'
                  l_err = ill(jsn)
                  k_err = jsn
                  goto 999

               end if

            end if


               ymaxd = 80.0

               if( iylog .eq. 1 )then
                 if( log10( ymax / ymin ) .gt. ymaxd ) then

                  ymin = 10.d0**( log10( ymax) - ymaxd )


                endif
               end if

               call xyrmax(ymin,ymax,iylog,ispay,ierr,
     &                     ypmin,ypmax)
cKN 2024/01/24

               if(ierr.ne.0) then
                  m_err = 'Y-Axis is LOG but Ymin is Negative.'
                  ErrCha = ''
                  ErrID = 'L:3039/R:a_main0/F:a-main0.f'
                  goto 999
               end if

               call xytic(iylog,ymin,ymax,iytic,yld,ysd,
     &                    nylta,nysta,ylta,ysta,ycta,nytch,
     &                    idecy,iystd,iyltd,mdecy,ystv,yltv,ierr,
     &                    ymul)


               if(ierr.ne.0) then
                  m_err = 'Error in Subroutine XYTIC, Sorry.'
                  ErrCha = ''
                  ErrID = 'L:3052/R:a_main0/F:a-main0.f'
                  goto 999
               end if



************************************************************************

*     START OF WRITING PS FILE

************************************************************************

*-----------------------------------------------------------------------
*     JAPANESE KANJI CODE
*-----------------------------------------------------------------------

            if( ifon .lt. 0 ) ifon = 0
            if( idat .lt. 0 ) idat = 0

*-----------------------------------------------------------------------
*     BASIC FONT AND TIC FONT
*-----------------------------------------------------------------------

            if( ibfon .lt. 0 ) ibfon = 0
            if( itfon .lt. 0 ) itfon = ibfon

*-----------------------------------------------------------------------
*        KANJI CODE FOR JAPANESE DATE
*-----------------------------------------------------------------------

            if( idat .eq. 0 ) then

                     today(1:11) = today1(1:11)

                     jtoday = 11

                     timestr(1:5) = timestr1(1:5)

                     jtime = 5

            else if( idat .eq. 1 .or.
     &               idat .eq. 2 .or.
     &               idat .eq. 3 ) then

                     write(month,'(i2)') imon
                     write(today( 1: 4),'(i4.4)') iyer
                     write(today(43:44),'(i2)') iday

                     today( 5:22) = nenn(idat)(1:18)
                     today(23:24) = month(1:2)
                     today(25:42) = tuki(idat)(1:18)
                     today(45:62) = niti(idat)(1:18)

                     jtoday = 62

                     write(timestr( 1: 2),'(i2)') ihor
                     write(timestr(21:22),'(i2)') imin

                     timestr( 3:20) = jikn(idat)(1:18)
                     timestr(23:40) = funn(idat)(1:18)

                     jtime = 40

            end if

*-----------------------------------------------------------------------
*     SPECIAL DEFAULT VALUE FOR COURIER FONT

*           FONT SIZE FOR AXIS TEXT   : FSTL
*           FONT SIZE FOR TICS NUMBER : FSTC
*           FONT SIZE FOR COMMENTS    : FSCM
*           WIDTH OF FRAME            : IDWF
*           WIDTH OF AXIS FRAME       : DWAF
*           WIDTH OF AXIS TCS         : DWAT

*-----------------------------------------------------------------------

               if( ibfon .ge. 9 .or. itfon .ge. 9 ) then

                  fstl = 17.0
                  fstc = 17.0
                  fscm = 17.0

                  idwf = 1
                  dwaf = 2.0
                  dwat = 2.0

               else

                  fstl = 24.0
                  fstc = 24.0
                  fscm = 17.0

                  idwf = 3
                  dwaf = 5.0
                  dwat = 4.0

               end if

*-----------------------------------------------------------------------
*        TMUP = Position Up of Minus in Tic Number : 0.08 * FNS
*-----------------------------------------------------------------------

               if( itfon .le. 3 ) then

                  tmup =  0.08

               else if( itfon .le. 8 ) then

                  tmup =  0.05

               else if( itfon .le. 12 ) then

                  tmup =  0.03

               else

                  tmup =  0.00

               end if

*-----------------------------------------------------------------------
*     WRITE HEADER AND BASIC OPERATORS FOR EPS
*-----------------------------------------------------------------------

         if( lpage .eq. 1 .and. izdo .eq. 1 ) then

            call writehd(jht,idbg,ifon)

         end if

*-----------------------------------------------------------------------
*     WRITE PAGE NUMBER
*-----------------------------------------------------------------------

         if( izdo .eq. 1 ) then

               if( idbg .eq. 1 )
     &         write(jht,'(/''%'',71(''*''),/
     &                      ''%     New Page Starts Here''
     &                    ,/''%'',71(''*''))')

               mpage = jpage

               write(jht,'(/''%%Page: '',i3,'' '',i3)') mpage, mpage
               write(jht,'( ''%%PageBoundingBox: (atend)'')')

*-----------------------------------------------------------------------
*cKN 2023/08/23

            call writehd1(jht,idbg,ifon)

*-----------------------------------------------------------------------

            if( ibsw .eq. 1 ) then

               if( idbg .eq. 1 )  write(jht,'()')

               write(jht,'(''gs'')')

            end if


            if( ifon .ne. 0 ) then

               if( idbg .eq. 1 )  write(jht,'()')

               write(jht,'(
     &         ''/ftnm {12 dict begin /FontName} N /fnm1 '',
     &         ''{N /FontType 0 N /WMode 0 N''/
     &         ''/FMapType 3 N /FontMatrix matrix N '',
     &         ''/Encoding [0 1] N /FDepVector [ } N''/
     &         ''/fend {FontName currentdict end definefont pop} N'')')

            end if


            if( mpage .eq. inpage ) then

               irot0 = irot

            end if

         end if


*-----------------------------------------------------------------------
*     FRAME ORIGIN FOR PORTRAIT OR LANDSCAPE
*-----------------------------------------------------------------------

         if( izdo .eq. 1 ) then

*-----------------------------------------------------------------------
*     INITIALIZE OF THE BOUNDING BOX
*-----------------------------------------------------------------------

               bx02 = 0.0
               bx01 = a4w
               by02 = 0.0
               by01 = a4h

               icb  = 0

*-----------------------------------------------------------------------
*     LANDSCAPE
*-----------------------------------------------------------------------

            if( irot .eq. 1 ) then

               wid  = a4h - 2.0 * hmg
               hgt  = a4w - 2.0 * wmg

               xog0 = wmg + wct + hgt
               yog0 = hmg + hct
               rog0 = 90.0
               sxg0 = 1.0
               syg0 = 1.0

               if( idbg .eq. 1 )
     &         write(jht,'(/''%'',71(''-''),/
     &                      ''% Frame Origin for Landscape''
     &                    ,/''%'',71(''-'')/)')

*-----------------------------------------------------------------------
*     PORTRAIT
*-----------------------------------------------------------------------

            else if( irot .eq. 0 ) then

               wid  = a4w - 2.0 * wmg
               hgt  = a4h - 2.0 * hmg

               xog0 = wmg + wct
               yog0 = hmg + hct
               rog0 = 0.0
               sxg0 = 1.0
               syg0 = 1.0

               if( idbg .eq. 1 )
     &         write(jht,'(/''%'',71(''-''),/
     &                      ''% Frame Origin for Portrait''
     &                    ,/''%'',71(''-'')/)')

*-----------------------------------------------------------------------
*     SMALL LANDSCAPE

*           1 / SQRT(2) SCALE DOWN OF LANDSCAPE
*                       scals = 1 / SQRT(2) = 0.70711
*-----------------------------------------------------------------------

            else if( irot .eq. -1 ) then

               wid  = a4h - 2.0 * hmg
               hgt  = a4w - 2.0 * wmg

               xog0 = wmg + wct
               yog0 = hmg + hct
               rog0 = 0.0
               sxg0 = 0.70711
               syg0 = 0.70711

               if( idbg .eq. 1 )
     &         write(jht,'(/''%'',71(''-''),/
     &                      ''% Frame Origin for Small Landscape''
     &                    ,/''%'',71(''-'')/)')

            end if

*-----------------------------------------------------------------------

            write(jht,'(2g14.5,'' TR '',g14.5,'' rotate '',/
     &                  2g14.5,'' scale'')')
     &                  xog0, yog0, rog0, sxg0, syg0

         end if

*-----------------------------------------------------------------------
*     BACK GROUND AND FRAME

*           IBGC = 1 ( Back Ground On)
*                  0 ( Back Ground Off)
*-----------------------------------------------------------------------

      if(    izdo .eq. 1 .and.
     &   ( ( clmo    .le. -r0max .and. clbg(1) .gt. 0.5 ) .or.
     &     ( clbg(1) .gt. -r0max .and. clbg(1) .le. 0.5 ) .or.
     &       nofr .eq. 0 ) ) then

*----------------------------------------------------------------------*

            if( idbg .eq. 1 )
     &      write(jht,'(/''%'',71(''-''),/
     &                   ''% Back Ground and Frame''
     &                 ,/''%'',71(''-'')/)')

*----------------------------------------------------------------------*
*     BACK GROUND COLOR OR GRAYSCALE
*----------------------------------------------------------------------*

         if( ( clmo    .le. -r0max .and. clbg(1) .gt. 0.0 ) .or.
     &       ( clbg(1) .gt. -r0max .and. clbg(1) .le. 0.0 ) ) then

            ibgc = 1

            write(jht,'(3f7.3,'' sc'')') clbg

            write(jht,'(
     &      ''0 0 M '',g14.5,'' 0 L '',2g14.5,'' L''/
     &      ''0 '',g14.5,'' L cp fl''
     &      )') wid, wid, hgt, hgt

            call bbox(1,0,0.d0,0.d0,0.d0)
            call bbox(1,0,wid,hgt,0.d0)
            icb = 1

         end if

*----------------------------------------------------------------------*
*     WRITE FRAME LINE

*           Default Line Width of Frame is IDWF = 3 dd

*----------------------------------------------------------------------*

         if( nofr .eq. 0 ) then

            if( clfr(1) .lt. -r0max ) clfr(1) = -2.0

            write(jht,'(3f7.3,'' sc sd0 '',i1.1,'' dd lw'')')
     &      clfr, idwf

            write(jht,'(
     &      ''0 0 M '',g14.5,'' 0 L '',2g14.5,'' L''/
     &      ''0 '',g14.5,'' L cp st''
     &      )') wid, wid, hgt, hgt

            icb = 0
            call bbox(1,0,0.d0,0.d0,0.d0)
            call bbox(1,0,wid,hgt,0.d0)
            icb = 1

         end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*        WRITE FRAME MESSAGE
*-----------------------------------------------------------------------

               if( noms .eq. 0 .and .izdo .eq. 1 ) then

                  call wfrmg(jht,idbg,ibfon,clms,imcm,wid,hgt,mcm,
     &                       clal,clmo,ifon,idat)

               end if


*-----------------------------------------------------------------------

************************************************************************
*     START MALTI-GRAPH IN ONE PAGE
************************************************************************

               if( idbg .eq. 1 )
     &         write(jht,'(/''%'',71(''=''),/
     &                      ''% Start No.'',I2,
     &                      '' Graph in Page '',I2
     &                    ,/''%'',71(''=''))') izdo, mpage


*-----------------------------------------------------------------------
*     Font definition
*-----------------------------------------------------------------------

            if( ifon .ne. 0 .and. izdo .gt. 1 ) then

               if( idbg .eq. 1 )  write(jht,'()')

               write(jht,'(
     &         ''/ftnm {12 dict begin /FontName} N /fnm1 '',
     &         ''{N /FontType 0 N /WMode 0 N''/
     &         ''/FMapType 3 N /FontMatrix matrix N '',
     &         ''/Encoding [0 1] N /FDepVector [ } N''/
     &         ''/fend {FontName currentdict end definefont pop} N'')')

            end if


*-----------------------------------------------------------------------
*     Rotout
*-----------------------------------------------------------------------

            if( iangzp .eq. 1 ) then

               if( idbg .eq. 1 )  write(jht,'()')

               write(jht,'(''rotout0'')')

               call bbox(1,1,0.d0,0.d0,-angzp)

            end if

*-----------------------------------------------------------------------
*     Scale
*-----------------------------------------------------------------------

            if( idbg .eq. 1 )
     &      write(jht,'(/''%'',71(''-''),/
     &                   ''% Scale''
     &                 ,/''%'',71(''-'')/)')


               write(jht,'(2g14.5,'' scale'')') scal, scal

               sxg0 = sxg0 * scal
               syg0 = syg0 * scal

*-----------------------------------------------------------------------
*     Axis Origin

*           XAL ; Default Length of X Axis is 14.0 cm

*           Origin from Frame corner
*              ( 4.0, 3.5 ) cm for Portrait
*              ( 4.5, 3.0 ) cm for Landscape

*-----------------------------------------------------------------------

               if( idbg .eq. 1 )
     &         write(jht,'(/''%'',71(''-''),/
     &                      ''% Axis Origin and Axis Size''
     &                    ,/''%'',71(''-'')/)')


               xal = 14.0  * xfac
               yal = xal * form

*-----------------------------------------------------------------------

            if( izdo .eq. 1 ) then

                  xalp = xal
                  yalp = yal

               if( irot .eq. 1 .or. irot .eq. -1 ) then

                  xog1 = 4.5
                  yog1 = 3.0

               else if( irot .eq. 0) then

                  xog1 = 4.0
                  yog1 = 3.5

               end if

                  xog2 = xalp * xorg
                  yog2 = yalp * yorg

            else if( izdo .gt. 1 ) then

                  xog1 = 0.0
                  yog1 = 0.0

                  xog2 = xalp * xorg / scal
                  yog2 = yalp * yorg / scal

            end if

            if( iangzp .eq. 1 ) then

                  dcs = cos( angzp / 180.0 * pi )
                  dsn = sin( angzp / 180.0 * pi )

                  xog3 = xog2 * dcs - yog2 * dsn
                  yog3 = xog2 * dsn + yog2 * dcs

                  xog2 = xog3
                  yog2 = yog3

            end if

                  xog1 = xog1 + xog2
                  yog1 = yog1 + yog2

*-----------------------------------------------------------------------

               write(jht,'(''/xor '',g14.5,'' cm N '',
     &                     ''/yor '',g14.5,'' cm N'')') xog1, yog1

*-----------------------------------------------------------------------

                  xalp = xal
                  yalp = yal

            if( angz .gt. -r0max ) then

                  iangz = 1

               if( iangzp .eq. 1 ) then

                  angz = angz + angzp

               end if

            else

               if( iangzp .eq. 1 ) then

                  angz = angzp
                  iangz = 1

               else

                  iangz = 0

               end if

            end if

                  angzp  = angz
                  iangzp = iangz

*-----------------------------------------------------------------------

               if( idbg .eq. 1 )  write(jht,'()')

            if( iangz .eq. 0 ) then

               write(jht,'(''xor yor TR'')')

                  call bbox(2,1,xog1,yog1,0.d0)

            else if( iangz .eq. 1 ) then

               write(jht,'(''xor yor '',g14.5,'' rotin0'')') angz

                  call bbox(2,1,xog1,yog1,angz)

            end if


               if( idbg .eq. 1 )  write(jht,'()')

               write(jht,'(''/xal '',g14.5,'' cm N '',
     &                     ''/yal '',g14.5,'' cm N'')') xal, yal

*-----------------------------------------------------------------------
*     Color Inside the Axis
*-----------------------------------------------------------------------

            if( ( clmo    .le. -r0max .and. clin(1) .gt. 0.0 ) .or.
     &          ( clin(1) .gt. -r0max .and. clin(1) .le. 0.0 ) .or.
     &            ibgc .eq. 1 ) then

               if( idbg .eq. 1 )
     &         write(jht,'(/''%'',71(''-''),/
     &                      ''% Color Inside the Axis''
     &                    ,/''%'',71(''-'')/)')

               if( ( clmo    .le. -r0max .and. clin(1) .gt. 0.0 ) .or.
     &             ( clin(1) .gt. -r0max .and. clin(1) .le. 0.0 ) ) then

                  write(jht,'(3f7.3,'' sc ax fl'')')  clin

               else

                  clin(1) = -1.0

                  write(jht,'(3f7.3,'' sc ax fl'')')  clin

               end if


                  call bbox(3,0,0.d0,0.d0,0.d0)
                  call bbox(3,0,1.d0,0.d0,0.d0)
                  call bbox(3,0,0.d0,1.d0,0.d0)
                  call bbox(3,0,1.d0,1.d0,0.d0)

            end if


*-----------------------------------------------------------------------
*     WRITE COLOR ON JHT READ FROM JHC
*-----------------------------------------------------------------------

      if( nhc .gt. 0 ) then

            if( ixlog .ne. 0 .or. iylog .ne. 0 ) then

               write(jht,'(/,''showpage'')')

               m_err = 'Please use LINEAR scal for Color Cluster Plot'
               ErrCha = ''
               ErrID = 'L:3648/R:a_main0/F:a-main0.f'
               goto 999

            end if

            rewind(jhc)

*-----------------------------------------------------------------------

      do jj = 1, nhc

            if( idbg .eq. 1 )
     &      write(jht,'(/''%'',71(''-''),/
     &                   ''% Color Cluster Plot''
     &                 ,/''%'',71(''-'')/)')

         if( idhc .eq. 0 ) then

            write(jht,'(
     &               ''/hs {/yps S N /xps S N'',/
     &               ''      np xps yps m xps wih add yps l '',
     &               ''xps wih add yps hgh add l '',/
     &               ''      xps yps hgh add l cp fl } N'')')

            if( idbg .eq. 1 )  write(jht,'()')

            idhc = idhc + 1

         end if

            read(jhc) wid0, hgt0

               wid  = zonep(xmin,xmax,wid0+xmin,ixlog)
               hgt  = zonep(ymin,ymax,hgt0+ymin,iylog)
               wid2 = zonep(xmin,xmax,wid0/2.0+xmin,ixlog)
               hgt2 = zonep(ymin,ymax,hgt0/2.0+ymin,iylog)

            write(jht,'(''/hgh '',g13.5,'' N /wih '',g13.5,'' N'')')
     &                     hgt, wid

            if( idbg .eq. 1 )  write(jht,'()')

            write(jht,'(''gs ax cl np '')')

            if( idbg .eq. 1 )  write(jht,'()')

*-----------------------------------------------------------------------

            do 401 ihc = 1, mhc(nhc)

               read(jhc) xpo0,ypo0,rcol

               if( xpo0 - wid0 / 2.0 .lt. xmax .and.
     &             xpo0 + wid0 / 2.0 .gt. xmin .and.
     &             ypo0 - hgt0 / 2.0 .lt. ymax .and.
     &             ypo0 + hgt0 / 2.0 .gt. ymin ) then

                     xpo = zonep(xmin,xmax,xpo0,ixlog) - wid2
                     ypo = zonep(ymin,ymax,ypo0,iylog) - hgt2

                     if( clmo .gt. -r0max .and. rcol(1) .gt. 0.0 )
     &               call ctomo(rcol,clmo)

                     write(jht,'(1p2e14.6,0p3f7.3,'' sc hs'')')
     &                       xpo,ypo,rcol

                     xp1 = xpo
                     xp2 = xpo + wid
                     yp1 = ypo
                     yp2 = ypo + hgt

                     xpm = max(0.0d0,xp1)
                     xp1 = min(1.0d0,xpm)
                     xpm = max(0.0d0,xp2)
                     xp2 = min(1.0d0,xpm)
                     ypm = max(0.0d0,yp1)
                     yp1 = min(1.0d0,ypm)
                     ypm = max(0.0d0,yp2)
                     yp2 = min(1.0d0,ypm)

                     call bbox(3,0,xp1,yp1,0.d0)
                     call bbox(3,0,xp1,yp2,0.d0)
                     call bbox(3,0,xp2,yp1,0.d0)
                     call bbox(3,0,xp2,yp2,0.d0)

               end if

  401       continue

            if( idbg .eq. 1 )  write(jht,'()')

            write(jht,'(''gr'')')

      end do
      end if


*-----------------------------------------------------------------------
*     WRITE CLUSTER ON JHT READ FROM JHL
*-----------------------------------------------------------------------

      if( nhl .gt. 0 ) then

            if( ixlog .ne. 0 .or. iylog .ne. 0 ) then

               write(jht,'(/,''showpage'')')

               m_err = 'Please use LINEAR scal for Cluster Plot'
               ErrCha = ''
               ErrID = 'L:3757/R:a_main0/F:a-main0.f'
               goto 999

            end if

*-----------------------------------------------------------------------

            if( idbg .eq. 1 )
     &      write(jht,'(/''%'',71(''-''),/
     &                   ''% Cluster Plot''
     &                 ,/''%'',71(''-'')/)')

         if( idhl .eq. 0 ) then

            write(jht,'(
     &               ''/hc {/hgh S N /wih S N'',/
     &               ''     /yps S hgh 2 div sub N '',
     &               ''/xps S wih 2 div sub N'',/
     &               ''      np xps yps m xps wih add yps l '',
     &               ''xps wih add yps hgh add l '',/
     &               ''      xps yps hgh add l cp fl } N'')')

            if( idbg .eq. 1 )  write(jht,'()')

            idhl = idhl + 1

         end if

            write(jht,'(''gs ax cl np'')')

            if( idbg .eq. 1 )  write(jht,'()')

*-----------------------------------------------------------------------

            rewind(jhl)

            do 400 ihc = 1, nhl

                  read(jhl) xpo,ypo,wid,hgt,rcol

                  wid = wid * clus
                  hgt = hgt * clus

               if( xpo - wid / 2.0 .lt. xmax .and.
     &             xpo + wid / 2.0 .gt. xmin .and.
     &             ypo - hgt / 2.0 .lt. ymax .and.
     &             ypo + hgt / 2.0 .gt. ymin ) then


                     xpo = zonep(xmin,xmax,xpo,ixlog)
                     ypo = zonep(ymin,ymax,ypo,iylog)
                     wid = zonep(xmin,xmax,wid+xmin,ixlog)
                     hgt = zonep(ymin,ymax,hgt+ymin,iylog)

                     if( clmo .gt. -r0max .and. rcol(1) .gt. 0.0 )
     &               call ctomo(rcol,clmo)

                     write(jht,'(1p4e12.4,0p3f7.3,'' sc hc'')')
     &                       xpo,ypo,wid,hgt,rcol

                     xp1 = xpo - wid / 2.0
                     xp2 = xpo + wid / 2.0
                     yp1 = ypo - hgt / 2.0
                     yp2 = ypo + hgt / 2.0

                     xpm = max(0.0d0,xp1)
                     xp1 = min(1.0d0,xpm)
                     xpm = max(0.0d0,xp2)
                     xp2 = min(1.0d0,xpm)
                     ypm = max(0.0d0,yp1)
                     yp1 = min(1.0d0,ypm)
                     ypm = max(0.0d0,yp2)
                     yp2 = min(1.0d0,ypm)

                     call bbox(3,0,xp1,yp1,0.d0)
                     call bbox(3,0,xp1,yp2,0.d0)
                     call bbox(3,0,xp2,yp1,0.d0)
                     call bbox(3,0,xp2,yp2,0.d0)

               end if

  400       continue

            if( idbg .eq. 1 )  write(jht,'()')

            write(jht,'(''gr'')')

      end if


*-----------------------------------------------------------------------
*     WRITE REGION ON JHT READ FROM JIL  ; NORMAL REGION
*-----------------------------------------------------------------------

         if( noner .gt. 0 ) then

                  if( idbg .eq. 1 )
     &            write(jht,'(/''%'',71(''-''),/
     &                         ''% Colored Region before Axis''
     &                       ,/''%'',71(''-'')/)')


                  rewind(jil)

                  call whone(1,jil,jht,noner,
     &                       ixlog,iylog,xmax,xmin,ymax,ymin,
     &                       ierr,clal,clmo,idbg,xfac,rlptl,sybw,erwd)
cKN 2024/01/24

                  if( ierr .ne .0 ) goto 999

         end if


*-----------------------------------------------------------------------
*     write bitmap on jht read from jhb, jhd and jhp
*-----------------------------------------------------------------------

      if( ibmap .gt. 0 ) then

            if( ixlog .ne. 0 .or. iylog .ne. 0 ) then

               write(jht,'(/,''showpage'')')

               m_err = 'Please use LINEAR scal for Bit Map hb:'
               ErrCha = ''
               ErrID = 'L:3883/R:a_main0/F:a-main0.f'
               goto 999

            end if

            rewind(jhd)
            rewind(jhb)
            rewind(jhp)
            rewind(jhr)

*-----------------------------------------------------------------------

      do kk = 1, ibmap

            if( idbg .eq. 1 )
     &      write(jht,'(/''%'',71(''-''),/
     &                   ''% Bit Map Plot''
     &                 ,/''%'',71(''-'')/)')

         if( idhc .eq. 0 ) then

            write(jht,'(
     &               ''/hs {/yps S N /xps S N'',/
     &               ''      np xps yps m xps wih add yps l '',
     &               ''xps wih add yps hgh add l '',/
     &               ''      xps yps hgh add l cp fl } N'')')

            if( idbg .eq. 1 )  write(jht,'()')

            idhc = idhc + 1

         end if

            read(jhd) iclip, nfrmm, nclmm, nbmap, npath,
     &                ihsb, iline, ilizt, rcolc, rcolb

                     if( clmo .gt. -r0max .and. rcolc(1) .gt. 0.0 )
     &               call ctomo(rcolc,clmo)

                     if( clmo .gt. -r0max .and. rcolb(1) .gt. 0.0 )
     &               call ctomo(rcolb,clmo)

*-----------------------------------------------------------------------
*        no frame, no cliping
*-----------------------------------------------------------------------

         if( iclip .eq. 0 .and. nfrmm .eq. 0 ) then

               if( idbg .eq. 1 )  write(jht,'()')

               write(jht,'(''gs ax cl np'')')

               if( idbg .eq. 1 )  write(jht,'()')

         end if

*-----------------------------------------------------------------------
*        frame
*-----------------------------------------------------------------------

         if( nfrmm .gt. 0 ) then

               if( idbg .eq. 1 )  write(jht,'()')

               write(jht,'(''gs ax cl np gs'')')

               if( idbg .eq. 1 )  write(jht,'()')

            do jj = 1, nfrmm

               read(jhb) xpo0, ypo0

                     xpo = zonep(xmin,xmax,xpo0,ixlog)
                     ypo = zonep(ymin,ymax,ypo0,iylog)

                  if( jj .eq. 1 ) then

                     write(jht,'(1p2e14.6,'' m'')') xpo,ypo

                  else

                     write(jht,'(1p2e14.6,'' l'')') xpo,ypo

                  end if

                     xp1 = xpo
                     yp1 = ypo

                     xpm = max(0.0d0,xp1)
                     xp1 = min(1.0d0,xpm)
                     ypm = max(0.0d0,yp1)
                     yp1 = min(1.0d0,ypm)

                     call bbox(3,0,xp1,yp1,0.d0)

            end do

               if( idbg .eq. 1 )  write(jht,'()')

                     write(jht,'(''cp cl np'')')

         end if

*-----------------------------------------------------------------------
*        clipping
*-----------------------------------------------------------------------

         if( iclip .ne. 0 ) then

            if( nfrmm .eq. 0 ) then

               if( idbg .eq. 1 )  write(jht,'()')

               write(jht,'(''gs ax cl np gs'')')

            end if

               if( idbg .eq. 1 )  write(jht,'()')

            do jj = 1, nclmm

               read(jhb) xpo0, ypo0

                     xpo = zonep(xmin,xmax,xpo0,ixlog)
                     ypo = zonep(ymin,ymax,ypo0,iylog)

                  if( jj .eq. 1 ) then

                     write(jht,'(1p2e14.6,'' m'')') xpo,ypo

                  else

                     write(jht,'(1p2e14.6,'' l'')') xpo,ypo

                  end if

                  if( nfrmm .eq. 0 ) then

                     xp1 = xpo
                     yp1 = ypo

                     xpm = max(0.0d0,xp1)
                     xp1 = min(1.0d0,xpm)
                     ypm = max(0.0d0,yp1)
                     yp1 = min(1.0d0,ypm)

                     call bbox(3,0,xp1,yp1,0.d0)

                  end if

            end do

               if( idbg .eq. 1 )  write(jht,'()')

                  if( iline .eq. 0 ) then

                     if( rcolb(1) .lt. -r0max ) then

                        write(jht,'(''cp cl np'')')

                     else

                        write(jht,'(''cp cl gs '',3f7.3,
     &                              '' sc fl np'')') rcolb

                     end if

                  else

                     if( rcolb(1) .lt. -r0max ) then

                        write(jht,'(''cp cl gs np'')')

                     else

                        write(jht,'(''cp cl gs '',3f7.3,
     &                              '' sc fl gs np'')') rcolb

                     end if

                  end if

               if( idbg .eq. 1 )  write(jht,'()')

         end if

*-----------------------------------------------------------------------

         if( npath .gt. 0 ) then

            do jj = 1, npath

               read(jhp) nptmm, rcol

                     if( clmo .gt. -r0max .and. rcol(1) .gt. 0.0 )
     &               call ctomo(rcol,clmo)

                  write(jht,'(3f7.3,'' sc'')') rcol

               do ii = 1, nptmm

                  read(jhp) xpo0, ypo0

                  xpo = zonep(xmin,xmax,xpo0,ixlog)
                  ypo = zonep(ymin,ymax,ypo0,iylog)

                  if( ii .eq. 1 ) then

                     write(jht,'(1p2e14.6,'' m'')') xpo, ypo

                  else if( ii .eq. nptmm ) then

                     write(jht,'(1p2e14.6,'' l cp fl'')') xpo, ypo

                  else

                     write(jht,'(1p2e14.6,'' l'')') xpo, ypo

                  end if

               end do

                  if( idbg .eq. 1 )  write(jht,'()')

               if( iclip .eq. 0 .and. nfrmm .eq. 0 ) then

                  xp1 = xpo
                  yp1 = ypo

                  xpm = max(0.0d0,xp1)
                  xp1 = min(1.0d0,xpm)
                  ypm = max(0.0d0,yp1)
                  yp1 = min(1.0d0,ypm)

                  call bbox(3,0,xp1,yp1,0.d0)

               end if

            end do

         end if

*-----------------------------------------------------------------------

         if( nbmap .gt. 0 ) then

            do ii = 1, nbmap

               read(jhr) nbtmm, wid0, hgt0

                  wid  = zonep(xmin,xmax,wid0+xmin,ixlog)
                  hgt  = zonep(ymin,ymax,hgt0+ymin,iylog)
                  wid2 = zonep(xmin,xmax,wid0/2.0+xmin,ixlog)
                  hgt2 = zonep(ymin,ymax,hgt0/2.0+ymin,iylog)

                  write(jht,'(''/hgh '',g13.5,'' N /wih '',
     &                        g13.5,'' N'')') hgt, wid

            do jj = 1, nbtmm

                  read(jhr) xpo0,ypo0,rcol

               if( xpo0 - wid0 / 2.0 .lt. xmax .and.
     &             xpo0 + wid0 / 2.0 .gt. xmin .and.
     &             ypo0 - hgt0 / 2.0 .lt. ymax .and.
     &             ypo0 + hgt0 / 2.0 .gt. ymin ) then

                     xpo = zonep(xmin,xmax,xpo0,ixlog) - wid2
                     ypo = zonep(ymin,ymax,ypo0,iylog) - hgt2

                     if( ihsb .eq. 0 ) then
                           rcol(1) = rcol(1) + 1.0d0
                     else

                        if( rcol(1) .gt. 0.0 ) then
                           rcol(1) = 2.5 * ( 1.0 - rcol(1) ) + 1.0
                        else
                           rcol(1) = -2.0 - rcol(1)
                        end if

                           rcol(1) = chue( rcol(1) )

                     end if

                     if( clmo .gt. -r0max .and. rcol(1) .gt. 0.0 )
     &               call ctomo(rcol,clmo)

                  if(  ihsb .eq. 1 ) then

                     write(jht,'(1p2e14.6,0p3f7.3,'' sc hs'')')
     &                       xpo,ypo,rcol

                  else

                     write(jht,'(1p2e14.6,0p3f7.3,'' sr hs'')')
     &                       xpo,ypo,rcol

                  end if

                  if( iclip .eq. 0 .and. nfrmm .eq. 0 ) then

                     xp1 = xpo
                     xp2 = xpo + wid
                     yp1 = ypo
                     yp2 = ypo + hgt

                     xpm = max(0.0d0,xp1)
                     xp1 = min(1.0d0,xpm)
                     xpm = max(0.0d0,xp2)
                     xp2 = min(1.0d0,xpm)
                     ypm = max(0.0d0,yp1)
                     yp1 = min(1.0d0,ypm)
                     ypm = max(0.0d0,yp2)
                     yp2 = min(1.0d0,ypm)

                     call bbox(3,0,xp1,yp1,0.d0)
                     call bbox(3,0,xp1,yp2,0.d0)
                     call bbox(3,0,xp2,yp1,0.d0)
                     call bbox(3,0,xp2,yp2,0.d0)

                  end if

               end if

            end do

            end do

         end if

*-----------------------------------------------------------------------

            if( idbg .eq. 1 )  write(jht,'()')

            if( iclip .ne. 0 ) then

               if( iline .eq. 0 ) then

                  if( rcolb(1) .lt. -r0max ) then

                     write(jht,'(''gr gr'')')

                  else

                     write(jht,'(''gr gr gr'')')

                  end if

               else

                  if( rcolb(1) .lt. -r0max ) then

                     write(jht,'(''gr '',3f7.3,
     &                           '' sc sd0'',i3,'' dd lw st gr gr'')')
     &                           rcolc, ilizt

                  else

                     write(jht,'(''gr gr '',3f7.3,
     &                           '' sc sd0'',i3,'' dd lw st gr gr'')')
     &                           rcolc, ilizt

                  end if

               end if

            else

                  write(jht,'(''gr'')')

            end if

      end do
      end if


*-----------------------------------------------------------------------
*     WRITE AXIS FRAME, TICS AND AXIS TEXT
*-----------------------------------------------------------------------

                  call wftat(jht,idbg,noaf,dwaf,afac,clax,dwat,
     &                       itic,jxtic,jytic,nxlta,nxsta,nylta,nysta,
     &                       idecx,idecy,ixexp,iyexp,
cKN 2024/01/24
     &                       txl,txr,tyd,tyu,noxn,noyn,
     &                       ifon,itfon,clnm,ixlog,iylog,
     &                       ixnum,iynum,nxtch,nytch,xcta,ycta,
     &                       xlta,xsta,ylta,ysta,tmup,
     &                       xtt,ixtt,ytt,iytt,
     &                       fstl,fstc,ibfon,cltx,atxs,ix,iy,
     &                       ixtxt,xtxp,iytxt,ytxp,noxt,noyt)

*-----------------------------------------------------------------------
*     WRITE HONE ON JHT READ FROM JOL  ; NORMAL LINES
*-----------------------------------------------------------------------

            if( noned .gt. 0 ) then

                  if( idbg .eq. 1 )
     &            write(jht,'(/''%'',71(''-''),/
     &                         ''% Lines, Markers and Error bars''
     &                       ,/''%'',71(''-'')/)')


                  rewind(jol)

                  call whone(0,jol,jht,noned,
     &                       ixlog,iylog,xmax,xmin,ymax,ymin,
     &                       ierr,clal,clmo,idbg,xfac,rlptl,sybw,erwd)
cKN 2024/01/24

                  if( ierr .ne .0 ) goto 999

            end if


*-----------------------------------------------------------------------
*     WRITE BOX
*-----------------------------------------------------------------------

            if( ioc .gt. 0 ) then

                  if( idbg .eq. 1 )
     &            write(jht,'(/''%'',71(''-''),/
     &                         ''% Write Boxes BOX:''
     &                       ,/''%'',71(''-'')/)')


                 do 475 i = 1, ioc

                     ioib = i

                     call wsbox(jht,idbg,clal,clmo,oxys,ioib,
     &                          xmin,xmax,ymin,ymax,ixlog,iylog)

  475            continue

            end if

*-----------------------------------------------------------------------
*     WRITE POLYGON
*-----------------------------------------------------------------------

            if( ipc .gt. 0 ) then

                  if( idbg .eq. 1 )
     &            write(jht,'(/''%'',71(''-''),/
     &                         ''% Write Polygons POLG:''
     &                       ,/''%'',71(''-'')/)')


                 do 472 i = 1, ipc

                     iplg = i

                     call wpolg(jht,idbg,clal,clmo,pxys,iplg,
     &                          xmin,xmax,ymin,ymax,ixlog,iylog)

  472            continue

            end if

*-----------------------------------------------------------------------
*     WRITE STAR
*-----------------------------------------------------------------------

            if( isc .gt. 0 ) then

                  if( idbg .eq. 1 )
     &            write(jht,'(/''%'',71(''-''),/
     &                         ''% Write Stars STAR:''
     &                       ,/''%'',71(''-'')/)')


                 do 473 i = 1, isc

                     isar = i

                     call wstar(jht,idbg,clal,clmo,sxys,isar,
     &                          xmin,xmax,ymin,ymax,ixlog,iylog)

  473            continue

            end if

*-----------------------------------------------------------------------
*     WRITE RIBBON
*-----------------------------------------------------------------------

            if( irc .gt. 0 ) then

                  if( idbg .eq. 1 )
     &            write(jht,'(/''%'',71(''-''),/
     &                         ''% Write Ribbon RIBN:''
     &                       ,/''%'',71(''-'')/)')


                 do 474 i = 1, irc

                     irib = i

                     call wribn(jht,idbg,clal,clmo,rxys,irib,
     &                          xmin,xmax,ymin,ymax,ixlog,iylog)

  474            continue

            end if

*-----------------------------------------------------------------------
*     WRITE ARROW BOUND   AB:
*-----------------------------------------------------------------------

            if( ibc .gt. 0 ) then

                  if( idbg .eq. 1 )
     &            write(jht,'(/''%'',71(''-''),/
     &                         ''% Arrow Bound AB:''
     &                       ,/''%'',71(''-'')/)')


                 do 471 i = 1, ibc

                        xi = zonep(xmin,xmax,cxys(i,1),ixlog)
                        xf = zonep(xmin,xmax,cxys(i,3),ixlog)
                        yi = zonep(ymin,ymax,cxys(i,2),iylog)
                        yf = zonep(ymin,ymax,cxys(i,4),iylog)

                        rcoll(1) = cxys(i,5)
                        rcoll(2) = cxys(i,6)
                        rcoll(3) = cxys(i,7)

                        rcolb(1) = cxys(i,11)
                        rcolb(2) = cxys(i,12)
                        rcolb(3) = cxys(i,13)

                        dd  = cxys(i,8)
                        rnl = cxys(i,9)
                        ara = cxys(i,10)

                     call warrw02(jht,idbg,clal,clmo,rcoll,rcolb,
     &                            xi,yi,xf,yf,xl,yl,dd,rnl,ara)

  471            continue

            end if

*-----------------------------------------------------------------------
*     WRITE ARROW ONLY   A:
*-----------------------------------------------------------------------

            if( ib .gt. 0 ) then

                  if( idbg .eq. 1 )
     &            write(jht,'(/''%'',71(''-''),/
     &                         ''% Arrow Only A:''
     &                       ,/''%'',71(''-'')/)')


                 do 470 i = 1, ib

                        xi = zonep(xmin,xmax,bxys(i,1),ixlog)
                        xf = zonep(xmin,xmax,bxys(i,3),ixlog)
                        yi = zonep(ymin,ymax,bxys(i,2),iylog)
                        yf = zonep(ymin,ymax,bxys(i,4),iylog)

                        rcol(1) = bxys(i,5)
                        rcol(2) = bxys(i,6)
                        rcol(3) = bxys(i,7)

                        dd   = bxys(i, 8)
                        rnl  = bxys(i, 9)
                        ara  = bxys(i,10)


                     call warrw01(0,jht,idbg,clal,clmo,
     &                            xi,yi,xf,yf,xl,yl,rcol,dd,rnl,ara)

  470            continue

            end if

*-----------------------------------------------------------------------
*     WRITE ARROW OF ARROW COMMENT  AW:
*-----------------------------------------------------------------------

                  acxl0 = 0.2 / xal
                  acxd0 = 0.1 / xal

            if( ia .gt. 0 ) then

                  if( idbg .eq. 1 )
     &            write(jht,'(/''%'',71(''-''),/
     &                         ''% Arrow of Arrow Comment AW:''
     &                       ,/''%'',71(''-'')/)')

               do 461 i = 1, ia

                        xi = zonep(xmin,xmax,axys(i,1),ixlog)
                        xf = zonep(xmin,xmax,axys(i,3),ixlog)
                        yi = zonep(ymin,ymax,axys(i,2),iylog)
                        yf = zonep(ymin,ymax,axys(i,4),iylog)

                        rcol(1) = axys(i,7)
                        rcol(2) = axys(i,8)
                        rcol(3) = axys(i,9)

                        dd   = axys(i,10)
                        rnl  = 1.0

                      if( iax(i) .eq. 0 ) then

                        if( axys(i,1) .lt. axys(i,3) ) then

                           iax(i) = 3

                         else

                           iax(i) = 1

                         end if

                      end if

                           acxl = acxl0 * axys(i,6)
                           acxd = acxd0 * axys(i,6)

                      if( iax(i) .eq. 1 ) then

                           xl = xi + acxl
                           xc = xl + acxd

                      else

                           xl = xi - acxl
                           xc = xl - acxd

                      end if

                           axys(i,5) = xc

                           ara = axys(i,11)

                     call warrw01(1,jht,idbg,clal,clmo,
     &                            xi,yi,xf,yf,xl,yi,rcol,dd,rnl,ara)

                           if( idbg .eq. 1 )  write(jht,'()')

  461          continue

            end if

*-----------------------------------------------------------------------
*     START OF COMMENTS
*-----------------------------------------------------------------------

*           Default Tic Font Size is 24 pt

*           YMARG = Margin of Text from number or axis : 1.0 * TFS

*-----------------------------------------------------------------------
*     WRITE TITLE OF GRAPH
*-----------------------------------------------------------------------

      if( ihtt .gt. 0 .and. notl .eq. 1 ) then

            if( idbg .eq. 1 )
     &      write(jht,'(/''%'',71(''-''),/
     &                   ''% Write Title of Graph''
     &                 ,/''%'',71(''-'')/)')

               iftgr = ibfon

               tfs = fstl * afac

               if( cltl(1) .lt. -r0max ) cltl(1) = -2.0

               tymag = 1.0

               xpo = tlxp * xal * cm
               ypo = tyu + tfs * tymag

               ixc  = 2
               iyc  = 1

*-----------------------------------------------------------------------

               call wtext(jht,chag,idbg,clal,clmo,0,xpo,ypo,
     &                    htt,ihtt,ixc,iyc,0,0.d0,iftgr,ifon,
     &                    0,indt,tfs,cltl)

                     if( idbg .eq. 1 )  write(jht,'()')

      end if

*-----------------------------------------------------------------------
*     WRITE ARROW COMMENT
*-----------------------------------------------------------------------

      if( ia .gt. 0 ) then

            if( idbg .eq. 1 )
     &      write(jht,'(/''%'',71(''-''),/
     &                   ''% Write Arrow Comment AW:''
     &                 ,/''%'',71(''-'')/)')

               rewind(jha)

        do 465 i = 1, ia

               if( ian(i) .eq. 0 ) goto 465

                  tfs = fscm * axys(i,6)

               rcol(1) = axys(i,7)
               rcol(2) = axys(i,8)
               rcol(3) = axys(i,9)

               if( clal(1) .gt. -r0max ) then
                  rcol(1) = clal(1)
                  rcol(2) = clal(2)
                  rcol(3) = clal(3)
               end if

               if( clmo .gt. -r0max .and. rcol(1) .gt. 0.0 )
     &            rcol(1) = -2.0
               if( rcol(1) .lt. -r0max ) rcol(1) = -2.0

               xpo = axys(i,5) * xal * cm
               ypo = zonep(ymin,ymax,axys(i,2),iylog) * yal * cm

               ixc = iax(i)
               iyc = 2

               read(jha) (ccmt(k),k=1,ian(i))

*-----------------------------------------------------------------------

               if( iaf(i) .lt. 0 ) then

                  ifarc = ibfon

               else

                  ifarc = iaf(i)

               end if

*-----------------------------------------------------------------------

                  call wtext(jht,chag,idbg,clal,clmo,0,xpo,ypo,
     &                       ccmt,ian(i),ixc,iyc,0,0.d0,ifarc,ifon,
     &                       0,indt,tfs,rcol)

                     if( idbg .eq. 1 )  write(jht,'()')

  465   continue

      end if

*-----------------------------------------------------------------------
*     WRITE VALUES OF THE PARAMETERS
*-----------------------------------------------------------------------

            if( nopa .eq. 1 ) then

                  call wparn(jht,idbg,'A',0.d0,
     &                       xmin,xmax,ixlog,ymin,ymax,iylog,ifon,
     &                       ibfon,fscm,clal,clmo,
     &                       parx,pary,pars,ipbox,clpa,clpb,clpl,clps)

            end if

*-----------------------------------------------------------------------
*     WRITE VALUES OF THE PARAMETERS
*-----------------------------------------------------------------------

            if( nocn .eq. 1 ) then

                  call wparn(jht,idbg,'C',0.5d0,
     &                       xmin,xmax,ixlog,ymin,ymax,iylog,ifon,
     &                       ibfon,fscm,clal,clmo,
     &                       conx,cony,cons,isbox,clcn,clnb,clnl,clns)

            end if

*-----------------------------------------------------------------------
*     WRITE LEGEND COMMENTS, BOX, LINES AND SYMBOLS
*-----------------------------------------------------------------------

            if( ncom .gt. 0 .and. nolg .eq. 1 ) then

                  call wlgnd(jht,idbg,txr,regx,regy,regs,regd,
     &                       ncom,ilbox,
     &                       xmin,xmax,ixlog,ymin,ymax,iylog,
     &                       ibfon,fscm,cllg,clal,clmo,
     &                       ifon,clgb,clgl,clgs,xfac,iycm,rycm,
     &                       rlptl,sybw,erwd)
cKN 2024/01/24

            end if

*-----------------------------------------------------------------------
*        default x coordinate for comment
*-----------------------------------------------------------------------

            if( ixlog .eq. 0 ) then

               rgxm = ( txr / cm / xal + 0.4 / xal )
     &              * ( xmax - xmin ) + xmin

            else if( ixlog .eq. 1 ) then

               rgxm = 10**( ( txr / cm / xal + 0.4 / xal )
     &                    * ( log10(xmax) - log10(xmin) )
     &                      + log10(xmin) )

            end if

            if( iylog .eq. 0 ) then

               rgym = ( tyu / cm / yal + 0.4 / yal )
     &              * ( ymax - ymin ) + ymin

            else if( iylog .eq. 1 ) then

               rgym = 10**( ( tyu / cm / yal + 0.4 / yal )
     &                    * ( log10(ymax) - log10(ymin) )
     &                      + log10(ymin) )

            end if

*-----------------------------------------------------------------------
*     WRITE COMMENTS
*-----------------------------------------------------------------------

            if( iw .gt. 0 .and. nocm .eq. 1 ) then

               call commt0(jht,idbg,ifon,
     &                     iw,wxys,iwx,iwy,iwn,iwf,iwb,ibox,cbox,
     &                     fscm,clal,clmo,xmax,xmin,ymax,ymin,
     &                     rgxm,rgym,
     &                     ixlog,iylog,ibfon)

            end if


*-----------------------------------------------------------------------
*     WRITE TABLE COMMENT
*-----------------------------------------------------------------------

            if( iv .gt. 0 ) then

               call commt1(jht,idbg,ifon,iv,
     &                     fscm,clal,clmo,xmax,xmin,ymax,ymin,
     &                     rgxm,rgym,
     &                     ixlog,iylog,ibfon)

            end if


*-----------------------------------------------------------------------
*     WRITE INCLUDE PS FILE
*-----------------------------------------------------------------------

            if( inps .gt. 0 ) then

               call winclp(jht,idbg,inps,
     &                     dpin,idpi,
     &                     psbbx,xyps,ixps,iyps,
     &                     xmax,xmin,ymax,ymin,ixlog,iylog)


            end if

************************************************************************
*     END OF MULTI GRAPH
************************************************************************

      rewind(jol)
      rewind(jil)
      rewind(jhc)
      rewind(jhl)
      rewind(jhy)

      rewind(jhq)
      rewind(jhd)
      rewind(jhb)
      rewind(jhp)
      rewind(jhr)

      if( iv .gt. 0 ) close(jwt)
      if( iw .gt. 0 ) close(jhm)
      if( ia .gt. 0 ) close(jha)

      goto 9010

 9000 continue

************************************************************************

*-----------------------------------------------------------------------
*     SUMMARY AND SHOW BOUNDING BOX

*        BBOX +- 3 AND +- 0.025 MARGIN

*        BBMG = 0.025
*        BBMG = 0.0

*-----------------------------------------------------------------------

               bbmg = 0.025
               bbmg = 0.0

               bx01 = bx01 - 3
               bx02 = bx02 + 3
               by01 = by01 - 3
               by02 = by02 + 3

               blng = max( bx02-bx01, by02-by01 ) * bbmg

               bx01 = bx01 - blng
               bx02 = bx02 + blng
               by01 = by01 - blng
               by02 = by02 + blng


               bxi1 = min( bxi1, bx01 )
               bxi2 = max( bxi2, bx02 )
               byi1 = min( byi1, by01 )
               byi2 = max( byi2, by02 )


         if( ibsw .eq. 1 ) then

               if( idbg .eq. 1 )
     &         write(jht,'(/''%'',71(''-''),/
     &                      ''% Show Bounding Box''
     &                    ,/''%'',71(''-'')/)')

               write(jht,'(''gr 1 1 1 sc sd0 1 dd lw'',2I5,'' M'',2I5,
     &                     '' L'',2I5,'' L'',2I5,'' L cp st'')')
     &         nint(bx01),
     &         nint(by01),
     &         nint(bx01),
     &         nint(by02),
     &         nint(bx02),
     &         nint(by02),
     &         nint(bx02),
     &         nint(by01)

         end if

*-----------------------------------------------------------------------

               if( idbg .eq. 1 ) then

                  write(jht,'(/''%'',71(''-''),/
     &                         ''      showpage''
     &                       ,/''%'',71(''-'')/)')

               else

                  write(jht,'(''showpage'')')

               end if

*-----------------------------------------------------------------------
*     SUMMARY OF ENDING PAGE
*-----------------------------------------------------------------------

               write(jht,'(''%%PageTrailer'')')

*-----------------------------------------------------------------------
*     PAGE BOUNDING BOX
*-----------------------------------------------------------------------

               write(jht,'(''%%PageBoundingBox:'',4i9)')
     &         nint(bx01),
     &         nint(by01),
     &         nint(bx02),
     &         nint(by02)

               if( idbg .eq. 1 )  write(jht,'()')

************************************************************************
*     END OF MULTI PAGE
************************************************************************


         goto 5000


************************************************************************

*-----------------------------------------------------------------------

 1000    continue

*-----------------------------------------------------------------------

           if( jpn .eq. 1 .or. jpn .eq. 3 .or. jpn .eq. 4 )
     &     lpage = lpage - 1

*-----------------------------------------------------------------------
*     ENDING OF PS FILE
*-----------------------------------------------------------------------

            if( idbg .eq. 1 )
     &      write(jht,'(/''%'',71(''*''),/
     &                   ''%     Ending of PS file''
     &                 ,/''%'',71(''*'')/)')

            write(jht,'(''%%Trailer'')')
            write(jht,'(''%%Pages: '',i3)') lpage

            write(jht,'(''%%BoundingBox:'',4i9)')
     &      nint(bxi1),
     &      nint(byi1),
     &      nint(bxi2),
     &      nint(byi2)

            write(jht,'(''%%EOF'')')

*-----------------------------------------------------------------------
*     WRITE 'angelbat' FOR UNIX SHELL SCRIPT
*-----------------------------------------------------------------------

         if( iangb .ge. 1 ) then

            if( irot0 .eq. 0 .or. irot0 .eq. -1 ) then

               write(jab,'('' angel-gp '',200a1)')
     &                       (dsjht(i:i), i=1,ioname),
     &                       (devnul(i:i), i=1,13)

            else if( irot0 .eq. 1 ) then

               write(jab,'('' angel-gs '',200a1)')
     &                       (dsjht(i:i), i=1,ioname),
     &                       (devnul(i:i), i=1,13)

            end if

               write(jab,*) ' rm angelbat'

         end if

               goto 1100

*-----------------------------------------------------------------------
*     END OF ERROR
*-----------------------------------------------------------------------

  999 open(jer,file='error.ang',status='unknown')

      write(jer,'('' ***** Error Message in ANGEL *****''/)')
      write(  6,'('' ***** Error Message in ANGEL *****''/)')
      write(jer,*) dsin(k_err)(1:idsi(k_err)),l_err,':'
      write(  6,*) dsin(k_err)(1:idsi(k_err)),l_err,':'

            icf = 200

         do 910 i = 200, 1, -1

            if( m_err(i:i) .ne. ' ' ) goto 911

  910    continue

  911       icf = i

      call ErrWrite(ErrID, ErrCha)
      write(jer,'('' ERROR = '',200a1)') ( m_err(i:i), i=1,icf )
      write(  6,'('' ERROR = '',200a1)') ( m_err(i:i), i=1,icf )

*-----------------------------------------------------------------------

 1100 continue

*-----------------------------------------------------------------------

      close(jer)

      close(jol)
      close(jil)
      close(jhc)
      close(jhl)
      close(jhy)

      close(jsi)
      close(jht)

      close(jhd)
      close(jhb)
      close(jhp)
      close(jhq)
      close(jhr)

      if( iv .gt. 0 ) close(jwt)
      if( iw .gt. 0 ) close(jhm)
      if( ia .gt. 0 ) close(jha)

      if( iangb .ge. 1 )  close(jab)

      do i = jsn, 1, -1

          call closef(jsi,jsn)

      end do

*-----------------------------------------------------------------------

      deallocate(chag)
      return
      end


