************************************************************************
*                                                                      *
      subroutine read00(jsi,ivers,ierr)
*                                                                      *
*       read input data                                                *
*       initialization of PHITS                                        *
*       modified by K.Niita on 2010/07/27                              *
*                                                                      *
*       output :                                                       *
*           jsi      unit of reading file                              *
*         ivers      version of input file                             *
*          ierr      error flag                                        *
*                                                                      *
************************************************************************
      use dedx_file
      use TDCHAINMOD, only: itdc,tdchain,talldcinit !FURUTA20200522
      use udm_Utility
      use moddas
      use liboutmod, only: libout ! frtati 2022/12/28
      use CHARVARMOD, only: Set_Charpara,ErrLine_Adjust, filnm,irwt

      use cvaloutmod ! S.H. 2023.10.27
      use t4dtrack_mod, only: read_t4dtrack

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

      include 'err.inc'

*-----------------------------------------------------------------------

      common /error/  m_err, l_err, k_err
      character       m_err*200

*-----------------------------------------------------------------------

      common /parai/  ipsq(400)
      common /bnkmem/ maxbnk, maxbn2, rtrckflp

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /dmpfil/ idmpf

      common /impmsg/ iimpn, isimp, iswct, isstr, maxip,
     &                mnimp(7,0:20),
     &                kfimp(7), inimc(7), inimt(7)

      common /fclmsg/ ifcln, isfcl, maxrg,
     &                mnfcl(6,0:20), mrfcl(kvlmax),
     &                kfcls(6), inflc(6), inflt(6)

      common /rclmsg/ ircln, isrcl, maxrr,
     &                mnrcl(6,0:20), mrrcl(kvlmax),
     &                krcls(6), lrcls(6), inrlc(6), inrlt(6),
     &                irpem(6,2), irman(6), irmct(6), jsmat(6)

      common /wwindp/ wupn, wsurvn, mxspln, mwhere, mvoww
      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww

      common /magreg/ nmreg, ingrc, ingrt, kmags, magtin, imgusr
      common /tmtreg/ ntmrg, intmc, intmt, ktime
      common /smireg/ nsreg, isgrc, isgrt, ksmir, ismir

      common /elmgrg/ nelctf, nereg, inerc, inert, kelcs

      common /splreg/ isptn, npreg(6), mnspt(6,0:20),
     &                ipgrc(6), ipgrt(6), ksplt(6), isplt(6),
     &                ispct(6,9), ispem(6,2)
      common /splrge/ espem(6,2)

      common /multipl/ imltp,iimlt(multmax),inmlt(multmax),
C MATSUDA 2024.11.25 (lagrange: ilmlt)
     &  ilmlt(multmax),
     &  idmlt(multmax),ismlt(multmax),impan(multmax),impat(multmax,6,3),
     &                 jmpat(multmax,6,6,2), imdfl(multmax)

      common /delesg/ rdels(kvlmax), mndel, ndels(kvlmax)
      common /delreg/ delm(kvlmax), kdelt

      common /tscreg/ ntscell(kvlmax)  ! ntcscell=0:no track structure, 1:water

      common /impreg/ dimp(kvlmax)
      common /volmsg/ rvols(kvlmax), mnvol, nvols(kvlmax)
      common /volreg/ dvol(kvlmax)
      common /brsmsg/ mnbrs, icbrs, mbbrs(kvlmax), cbrem(49)
      common /tmpmsg/ rtmps(kvlmax), mntmp, ntmps(kvlmax)
      common /tmpreg/ dtmp(kvlmax)
      common /pwtmsg/ rpwts(kvlmax), mnpwt, npwts(kvlmax)
      common /pwtreg/ dpwt(kvlmax)
      common /regdd/  ivolm, iimpo
      common /regdn/  irden
      common /regim/  iimps
      common /paraj/  mstz(300), parz(300) ! S.H. added (2017.8.11)
      common /wtcntl/ iwt, icimp(20), ifcls(20), iwwin(20), ircls(20)
      common /cggmm/  ngstar, ngfini, ngfin0
      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)
      common /mttmc/  smttc(kvlmax), mttcn, mttc1(kvlmax), mttc2(kvlmax)
      common /mtnmc/  smtnc(kvlmax), dmtnc(kvlmax,2),
     &                mtncn, mtnc(kvlmax,2), nmtnc(kvlmax,2)
      character dmtnc*80
      common /mtreg/  smtrg(kvlmax), dmtrg(kvlmax),
     &                mtrgn, mtrg(kvlmax,2), nmtrg(kvlmax)
      character dmtrg*80

      common /elreg/  elarg(kvlmax,4), mlrgn, melrg(kvlmax,2), ielusr

      common /fgdata/ ifrgd, ifgm, ifgdf(kvlmax,5), frgfl(kvlmax)
      character frgfl*200

      common /regcrs/ icrsflx
      common /fsmul/  ridfs, iidfs

      common /tscmsg/ ktsc(kvlmax), mntsc, ntsc(kvlmax)
      data mntsc /0/

      data iidfs /0/
      data ridfs /-1.0d0/
      data icrsflx /0/
      data idmpf /0/
      data ifrgd /0/

      data ncntc /0,0,0/
      data ngstar /0/
      data ngfini /0/
      data ngfin0 /0/

      data ivolm /0/
      data iimpo /0/
      data irden /0/
      data mnpwt /0/
      data mtncn /0/
      data mtrgn /0/
      data mlrgn /0/
      data mttcn /0/
      data mtreg /0/
      data imltp /0/

      data cbrem / 49 * 0.0d0 /
      data nmtnc / kvlmax*0, kvlmax*0 /

*-----------------------------------------------------------------------

      common /htitl/  iclgt(100), ctitl(100)
      character       ctitl*200
      common /inpec/  ititl, ipara, ibody, iregn, llarr, itby, itar
      common /inggs/  iog, igcel, ioa, igsuf, iob, igtrs
      common /tcntl/  icntl, inucr
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /regdu/  iuni(kvlmax)
      common /regde/  ichp(kvlmax), ilat(kvlmax), idct(kvlmax)
      common /regdf/  ioc, ifilt, ifil(kvlmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /kmat1a/ mxmat, mxmat0, mxnel
      common /kmat1p/ kmout
      common /ivoll/ ivoll

*-----------------------------------------------------------------------

      common /gamint/ initga

      data initga /0/

*-----------------------------------------------------------------------

      data idmn(0) /0/

      data kmout /0/

      data mxmat / 0 /
      data iuni / kvlmax*0 /
      data ilat / kvlmax*0 /
      data idct / kvlmax*0 /
      data ifil / kvlmax*0 /
      data ifilt / 0 /

      data mrfcl / kvlmax*0 /
      data mrrcl / kvlmax*0 /

      data delm / kvlmax*1.d+10 /

      data ntscell / kvlmax*0 /

      data dimp / kvlmax*1.d0 /
      data dvol / kvlmax*1.d0 /
      data dtmp / kvlmax*2.585d-8 /   ! T.Sato 2022/03/20 to adjust temparature in JENDL-4.0
      data dpwt / kvlmax*-1.0d0 /

      data idgr / kvmmax*0 /
      data idnm / kvmmax*0 /

      data iwt   / 0 /
      data icimp / 20*0 /
      data ifcls / 20*0 /
      data ircls / 20*0 /
      data iwwin / 20*0 /

      data maxip / 0 /
      data maxrg / 0 /
      data maxrr / 0 /
      data maxww / 0 /

      data iswct / 0 /
      data isimp / 0 /
      data iimps / 0 /
      data iimpn / 0 /
      data ifcln / 0 /
      data ircln / 0 /
      data isptn / 0 /
      data iwwdp / 0 /
      data iswwp / 0 /

      data junf / 0 /

*-----------------------------------------------------------------------

      common /wwind0/ dvals, iwmsh, kecho, kgwwp(6), idval
      data iwmsh / 0 /

*-----------------------------------------------------------------------
      common /wwbiasn/ eenwb(6,100), iwbdp, mnwbp(6,0:20), kfwbp(6),
     &                 inwbc(6), ienwb(6), inwbt(6), iswbp, maxwb
      data iwbdp / 0 /
      data iswbp / 0 /
      data maxwb / 0 /

*-----------------------------------------------------------------------

      common /igsherr/ igsher, icl01, icl02
      data igsher / 0 /

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      logical   exex

      data ititl /0/
      data ipara /0/
      data ibody /0/
      data iregn /0/
      data llarr /0/

      data igcel /0/
      data igsuf /0/
      data igtrs /0/

      data itnm  /0/
      data idcnm /0/
      dimension ild(0:9)

      data nereg /0/

      data nmreg /0/
      data ntmrg /0/
      data nsreg /0/
      data mnvol /0/
      data mndel /0/
      data mntmp /0/
      data mnbrs /-1/
      data mmmax /1/

      data ndels / kvlmax*0 /
      data rdels / kvlmax*1.0d0 /
      data nvols / kvlmax*0 /
      data rvols / kvlmax*1.0d0 /
      data ntmps / kvlmax*0 /
      data npwts / kvlmax*0 /
      data mbbrs / kvlmax*0 /
      data rtmps / kvlmax*2.585d-8 /  ! T.Sato 2022/03/20 to adjust temparature in JENDL-4.0
      data rpwts / kvlmax*-1.0d0 /
      data iodw  /201/

*-----------------------------------------------------------------------
      common /canat/  ianat
      data ianat /0/

*-----------------------------------------------------------------------
      common /ggcell/ icells, iobo
      character chcfg*100
      data iobo /0/

      common /mpi00/ npe, me
      character chprojall*200
      character chme*5

*-----------------------------------------------------------------------
      itdc = 0 !FURUTA20200522
*-----------------------------------------------------------------------

            ierr = 0
            ivoll = 0

            ireadpara = 0 ! S.H. added (2017.8.14)
            ifileflag = 0 !FURUTA20230208

*-----------------------------------------------------------------------
*     input file number constant
*-----------------------------------------------------------------------

            jsn = 0

               dsin(0) = 'Error Line'
               idsi(0) = 12
               m_err = ' Error !!'
               ErrCha = ''  ! Initialize Error message
               ErrID = 'L:316/R:read00/F:read00.f'
               l_err = 1
               k_err = 0

               SOURCE_Eflag = 0

*-----------------------------------------------------------------------
      call moddas_initialize()

*-----------------------------------------------------------------------
*     read first non comment line of unit 5
*-----------------------------------------------------------------------

  514 continue

            read(jsi,'(a200)', iostat = ios ) chin
            if( ios .eq. -1 ) goto 998

                  call chlngt(chin,200,i1,i2)
                  chlw = chin
                  call chcaps(chlw,i1,i2,i3,'#!$')
                  chcm = chlw
                  call chcomp(chcm,i1,i3,i4)

*-----------------------------------------------------------------------

            if( i1 .eq. 0 .and. i2 .eq. 0 ) then
               iskip = 1
            else if( index('#!$',chin(i1:i1)) .ne. 0 ) then
               iskip = 2
            else
               iskip = 0
            end if

            if( iskip .ne. 0 ) goto 514

*-----------------------------------------------------------------------
*     check unit 5 file
*     if first line is started as 'file =' -> this is input file name
*                                  other   -> read unit jsi = 5
*-----------------------------------------------------------------------

            if( chcm(i1:i1+4) .ne. 'file=' ) then

               rewind jsi

               goto 101

            end if

*-----------------------------------------------------------------------
*     get input file name from file =
*-----------------------------------------------------------------------

               ifileflag = 1 !FURUTA20230208

               ic1 = inumc(chlw,i1,i3,'=') + 1
               ic1 = jnumc(chlw,ic1,i3)
               ic2 = min( i3, inumc(chlw,ic1+1,i3,' ') )

               iname = ic2 - ic1 + 1

               do i = 1, iname

                  dsin(jsn+1)(i:i) = chin(ic1+i-1:ic1+i-1)

               end do

               do i = iname + 1, 200

                  dsin(jsn+1)(i:i) = ' '

               end do

                  idsi(jsn+1) = iname

*-----------------------------------------------------------------------
*     input file dsin exist ?
*     open first input file
*-----------------------------------------------------------------------

            inquire( file = dsin(jsn+1), exist = exex )

            if( exex .eqv. .false. ) then

               m_err = 'Input File Name Error. File not exist'//
     &                 ' ->> '//dsin(jsn+1)(1:iname)
               ErrCha = ''
               ErrID = 'L:404/R:read00/F:read00.f'
               l_err = 1
               k_err = jsn+1

               dsin(k_err) = 'Error Line'
               idsi(k_err) = 12

               goto 999

            end if

               call openf(jsi,jsn,dsin)

*-----------------------------------------------------------------------
*     Relpace set Character variable
*-----------------------------------------------------------------------
 101           call set_Charpara(jsi,ifileflag,dsin(jsn),ierr)
               if(ierr.ne.0)return

*-----------------------------------------------------------------------
*     check input file : whether [  ] is exist or not
*-----------------------------------------------------------------------

  100    continue

            ivers = 0

            read(jsi, '(a10)', iostat = ios ) chin
            if( ios .eq. -1 ) goto 200

            call chlngt(chin,10,i1,i2)

            if( i1.ne.0 ) then
              if( chin(i1:i1) .ne. '[' ) goto 100
            end if

            ivers = 1

  200    continue

            rewind jsi

*-----------------------------------------------------------------------
*        old version
*-----------------------------------------------------------------------

            if( ivers .eq. 0 ) then

               iot = 6
               ict = 0

               call read01(jsi,iot,ict,ierr)

               return

            end if

*-----------------------------------------------------------------------
*     read input file
*     initial value for include file
*-----------------------------------------------------------------------

         do i = 0, 9

            ill(i) = 0
            ilf(i) = 10000000
            ild(i) = 0

         end do

            jpn  = 0


************************************************************************
*                                                                      *
*     start of read sections                                           *
*                                                                      *
************************************************************************

C S.H. added for reading [parameters] at first (2017.8.11)
*-----------------------------------------------------------------------
 141        continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 999
               if( jpn .eq. 3 ) then
                   jsn = 1
                   goto 995 ! [parameters] section is mondatory.
               endif
               if( iskip .ne. 0 ) goto 141

*-----------------------------------------------------------------------
*        [parameters]
*-----------------------------------------------------------------------

            if( chcm(i1:i1+11) .eq. '[parameters]' ) then

               if( chcm(i1+12:i1+14) .eq. 'off' ) then

                  jpn = 2
                  goto 141

               end if

               call param(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  goto 142

            end if
            goto 141

 142        continue

C for the case that [parameters] is the last section
            if ( jsi.eq.5 .and. jpn.eq.3 ) then
             if ( ifileflag .eq. 0 ) then ! input is given by redirection
              rewind jsi
             else
              jsi = 31
              jsn = jsn + 1
              if ( irwt .eq. 0 ) then ! $RWT = 0, open the file of jsn=1
               open(unit=jsi,file=trim(dsin(jsn)),status='unknown')
              else ! $RWT > 0, open the rewrite file
               open(unit=jsi,file=trim(filnm),status='unknown')
              end if
             end if
            else
             do while ( jsi .gt. 31 )
              call closef(jsi,jsn)
             end do
             rewind jsi
            end if

            do i = 0, 9
               ill(i) = 0
               ilf(i) = 10000000
               ild(i) = 0
            end do
            jpn  = 0
            if ( mstz(1) .eq. 16 ) then
               call anainp(jsi,jsn,ivers,ierr)
               return
            end if

            if ( mstz(161) .ne. 0 ) then
               call allocate_cvalout
               imemorycval = 1
            end if

*-----------------------------------------------------------------------
C S.H. added for reading [parameters] at first (2017.8.11)

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue


            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 999
               if( jpn .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

  240 continue

*-----------------------------------------------------------------------
*        [title]
*-----------------------------------------------------------------------

            if( chcm(i1:i1+6) .eq. '[title]' ) then

               if( chcm(i1+7:i1+9) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call title(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [parameters]
*-----------------------------------------------------------------------
c skip reading [parameters] here (S.H. revised 2017.8.11)

            else if( chcm(i1:i1+11) .eq. '[parameters]' ) then

               if( chcm(i1+12:i1+14) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

c error when [parameters] are read two times (S.H. added 2017.8.14)
                  ireadpara = ireadpara + 1
                  if ( ireadpara .ge. 2 ) goto 996

 143              continue
                  call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
                  if( ierr .ne. 0 ) goto 999
                  if( jpn .eq. 3 ) return

! --- add NS 2020.04 ---
! for DEBUG_BOUNDARY
                  if (i1 .ne. 0) then
                  if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then
                     jpn = 1
                     goto 240
                  else
                     goto 143
                  end if
! for DEBUG_BOUNDARY
                  else
                    goto 143
                  endif

*-----------------------------------------------------------------------
*        [material]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+9) .eq. '[material]' ) then

               if( chcm(i1+10:i1+12) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call mater(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)


                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [body]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 5) .eq. '[body]' ) then

               if( chcm(i1+6:i1+8) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call bodyin(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [region]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 7) .eq. '[region]' ) then

               if( chcm(i1+8:i1+10) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call region(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [cell]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 5) .eq. '[cell]' ) then

               if( chcm(i1+6:i1+8) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call gcell(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [surface]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 8) .eq. '[surface]' ) then

               if( chcm(i1+9:i1+11) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call gsurf(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [transform]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+10) .eq. '[transform]' ) then

               if( chcm(i1+11:i1+13) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call gtrans(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [array]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 6) .eq. '[array]' ) then

               if( chcm(i1+7:i1+9) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call arrayg(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [source]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 7) .eq. '[source]' ) then

               if( chcm(i1+8:i1+10) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

c [source] set check flag
               SOURCE_Eflag = 1


c delete the extracted mulit-souce sub-section file
               if( npe .le. 1 ) then
                  chprojall = 'risrc-extract.tmp'
               else
                  iorder = aint(log10(real(npe))) + 1
                  if ( iorder .lt. 3) iorder = 3
                  write(chme,'(i5.5)') me
                  chprojall = 'risrc-extract.tmp'//chme(6-iorder:5)
               endif

               inquire(file=chprojall,exist=exex)
               if ( exex ) then
                  open (80,file=chprojall) ! 80 <- jsi+1 (y.sakaki,2023/01)
                  close(80,status='delete')
               end if

c proj=all, make multi-source sub-section data(file)
               npcunt = -1   ! input: flag for proj=all, output: total incident particles
               call sours(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,npcunt,ierr)

               if( ierr .ne. 0 ) goto 999

               if( jpn .eq. 3 ) return

               if ( npcunt .gt. 0 ) then     ! multi-source sub-section
                  dsin(jsn+1) = chprojall
                  idsi(jsn+1) = 17
                  ill (jsn+1) = 0
                  ilf (jsn+1) = 100000000
                  call openf(jsi,jsn,dsin)    ! expanded multi-source sub-section

                  npcunt = 0   ! input: flag for ignore proj=all, output: total incident particles
                  call sours(jsn,jsi,dsin,idsi,ill,ilf,
     &                       jpn,chin,chlw,chcm,i1,i2,i3,i4,npcunt,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

c switch the file to be read to the original input file.
                  call closef(jsi,jsn)
               end if

c delete the extracted mulit-souce sub-section file
               inquire(file=chprojall,exist=exex)
               if ( exex ) then
                  open (80,file=chprojall) ! 80 <- jsi+1 (y.sakaki,2023/01)
                  close(80,status='delete')
               end if

               call caltot()

                  goto 240

*-----------------------------------------------------------------------
*        [magnetic field]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+14) .eq. '[magneticfield]' ) then

               if( chcm(i1+15:i1+17) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call mgnet(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [electro magnetic field]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+21) .eq. '[electromagneticfield]' ) then

               if( chcm(i1+22:i1+24) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call elmgf(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [charge state]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+12) .eq. '[chargestate]' ) then

               if( chcm(i1+13:i1+15) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if


                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [Delta Ray]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+9) .eq. '[deltaray]' ) then

               if( chcm(i1+10:i1+12) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call delte(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [Track Structure]  ! T.Sato 2017/5/09
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+15) .eq. '[trackstructure]' ) then

               if( chcm(i1+16:i1+18) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call trackst(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [volume]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+7) .eq. '[volume]' ) then

               if( chcm(i1+8:i1+10) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if
               ivoll = 1

               call volum(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [importance]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+11) .eq. '[importance]' ) then

               if( chcm(i1+12:i1+14) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call impot(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

! T.Sato 2024/06/21, maxbnk is automatically changed to 1000000 if not defined
               itmp=0
               do kk = 1, ipara
                if( ipsq(kk) .eq. 34 ) itmp=1
               end do
               if(itmp.eq.0) then
                maxbnk=1000000
                mstz( 34 ) = 1000000 ! in case overwrite later
               endif
               if(mstz(85).eq.1 .and. mstz(23).eq.1) then   ! When [importace] is used in egs, over-division occurs if igchk=1.
                write(*,'(a75)') "*** Warning: negs=1 and [importace]
     & are used, igchk is changed to igchk=0."
                 mstz(23) = 0
               endif
                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [forced collisions]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+17) .eq. '[forcedcollisions]' ) then

               if( chcm(i1+18:i1+20) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call foccl(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [repeated collisions]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+19) .eq. '[repeatedcollisions]' ) then

               if( chcm(i1+20:i1+22) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call repcl(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [splitting]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+10) .eq. '[splitting]' ) then

               if( chcm(i1+11:i1+13) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call split(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [super mirror]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+12) .eq. '[supermirror]' ) then

               if( chcm(i1+13:i1+15) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call supmir(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [temperature]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+12) .eq. '[temperature]' ) then

               if( chcm(i1+13:i1+15) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call tempe(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [bbrem]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+10) .eq. '[bremsbias]' ) then

               if( chcm(i1+11:i1+13) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call bremb(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [photon weight]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+13) .eq. '[photonweight]' ) then

               if( chcm(i1+14:i1+16) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call photw(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [counter]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 8) .eq. '[counter]' ) then

               if( chcm(i1+9:i1+11) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call counts(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [timer]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 6) .eq. '[timer]' ) then

               if( chcm(i1+7:i1+9) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call timers(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [matnamecolor]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+13) .eq. '[matnamecolor]' ) then

               if( chcm(i1+14:i1+16) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call matnc(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [mattimechange]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+14) .eq. '[mattimechange]' ) then

               if( chcm(i1+15:i1+17) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call mattc(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [elastic option]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+14) .eq. '[elasticoption]' ) then

               if( chcm(i1+15:i1+17) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call elasop(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [regname]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+8) .eq. '[regname]' ) then

               if( chcm(i1+9:i1+11) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call regnm(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [fragdata]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+9) .eq. '[fragdata]' ) then

               if( chcm(i1+10:i1+12) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call frgdat(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [WW Bias]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+7) .eq. '[wwbias]' ) then

               if( chcm(i1+8:i1+10) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call wwbias(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [weightwindow]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+13) .eq. '[weightwindow]' ) then

               if( chcm(i1+14:i1+16) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call wwind(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

! T.Sato 2024/06/21, maxbnk is automatically changed to 1000000 if not defined
               itmp=0
               do kk = 1, ipara
                if( ipsq(kk) .eq. 34 ) itmp=1
               end do
               if(itmp.eq.0) then
                maxbnk=1000000
                mstz( 34 ) = 1000000 ! in case overwrite later
               endif
                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [multiplier]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+11) .eq. '[multiplier]' ) then

               if( chcm(i1+12:i1+14) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call multip(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [datamax]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+8) .eq. '[datamax]' ) then

               if( chcm(i1+9:i1+11) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call ndatmax(jsn,jsi,dsin,idsi,ill,ilf,
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [libout] ! frtati 2022/12/28
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+7) .eq. '[libout]' ) then

               if( chcm(i1+8:i1+10) .eq. 'off' .or.
     &             mstz(1).ne.1 .or. mstz(15).ne.100 ) then

                  jpn = 2
                  goto 140

               end if

               call libout(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [anatally]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+9) .eq. '[anatally]' ) then

               if( chcm(i1+10:i1+12) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call read_anatal(jsn,jsi,dsin,idsi,ill,ilf,
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [User Defined Interaction]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+23) .eq.
     &                                 '[userdefinedinteraction]' ) then

               if( chcm(i1+24:i1+26) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call read_udinteract(jsn,jsi,dsin,idsi,ill,ilf,
     &                              jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [User Defined Particle]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+20) .eq. '[userdefinedparticle]' ) then

               if( chcm(i1+21:i1+23) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call read_udpart(jsn,jsi,dsin,idsi,ill,ilf,
     &                          jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*        [t-4Dtrack]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+10) .eq. '[t-4dtrack]' ) then

               if( chcm(i1+11:i1+13) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call read_t4dtrack(jsn,jsi,dsin,idsi,ill,ilf,
     &                          jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .eq. 900 ) then ! t-4dtrack skipped
                   jpn = 2
                   goto 140
                  endif

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

            end if ! S.H. added (2017.8.11)

*-----------------------------------------------------------------------
* For tally
*-----------------------------------------------------------------------
           if ( mstz(1).ne.14 .and. mstz(1).ne.15 ) then ! for icntl.ne.14,15

*-----------------------------------------------------------------------
*        [t-deposit2]
*-----------------------------------------------------------------------

            if( chcm(i1:i1+ 11) .eq. '[t-deposit2]' ) then ! S.H. revised (2017.8.11)

               if( chcm(i1+12:i1+14) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call tdepsit2(jsn,jsi,dsin,idsi,ill,ilf,
     &                       jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-deposit]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 10) .eq. '[t-deposit]' ) then

               if( chcm(i1+11:i1+13) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call tdeposit(jsn,jsi,dsin,idsi,ill,ilf,
     &                       jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-point]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 8) .eq. '[t-point]' ) then

               if( chcm(i1+9:i1+11) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call tpoint(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-volume]
*-----------------------------------------------------------------------
c skip reading [t-volume] here (S.H. revised 2017.8.11)

            else if( chcm(i1:i1+ 9) .eq. '[t-volume]' ) then

                  jpn = 2
                  goto 140

*-----------------------------------------------------------------------
*        [t-wwbg]
*-----------------------------------------------------------------------
c skip reading [t-wwbg] here (S.H. revised 2017.8.11)

            else if( chcm(i1:i1+ 7) .eq. '[t-wwbg]' ) then

                  jpn = 2
                  goto 140

*-----------------------------------------------------------------------
*        [t-wwg]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 6) .eq. '[t-wwg]' ) then

               if( chcm(i1+7:i1+9) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call twwg(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-sed]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 6) .eq. '[t-sed]' ) then

               if( chcm(i1+7:i1+9) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call tsed(jsn,jsi,dsin,idsi,ill,ilf,
     &                   jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-let]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 6) .eq. '[t-let]' ) then

               if( chcm(i1+7:i1+9) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call tlet(jsn,jsi,dsin,idsi,ill,ilf,
     &                   jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-track]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 8) .eq. '[t-track]' ) then

               if( chcm(i1+9:i1+11) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call ttrack(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-adjoint]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+10) .eq. '[t-adjoint]' ) then

               if( chcm(i1+11:i1+13) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call tadjnt(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-cross]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 8) .eq. '[t-cross]' ) then

               if( chcm(i1+9:i1+11) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call tcross(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-yield]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 8) .eq. '[t-yield]' ) then

               if( chcm(i1+9:i1+11) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call tyield(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-dchain]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 9) .eq. '[t-dchain]' ) then

               if( chcm(i1+10:i1+12) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if


               if(itdc.eq.0) then
                call talldcinit(itlmax) !FURUTA20200522
                if(mstz(35).ne.1) then
                 write(ErrCha,'("Warning: jmout is set to 1 ",
     &           "because [t-dchain] is defined")')
                 ErrID = 'L:1823/R:read00/F:read00.f'
                 call ErrWrite(ErrID,ErrCha)
                 mstz(35)=1 ! jmout should be 1 when [t-dchain] is defined
                endif
               endif

               call tdchain(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                    jpnd   = 0
                  if( jpn .eq. 3 ) then
                     jpn   = 0
                    jpnd   = 3
                  end if

              iodw = 201
              write(iodw,'("[end]")')
             rewind(iodw)

  145 continue

            call readl(jsn,iodw,dsin,idsi,ild,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 999
               if( jpn .eq. 3 ) return
            if( iskip .ne. 0 ) goto 145

            if( chcm(i1:i1+ 11) .eq. '[t-track]off' ) then

               call ttrack(jsn,iodw,dsin,idsi,ild,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
                    if( ierr .ne. 0 ) goto 999
                    if( jpn .eq. 3 ) return
               backspace(jsi,err=146)  ! T.Sato 2024/11/19, avoid "without [end] or other section behind [t-dchain]" issue
               ill(jsn) = ill(jsn) - 1
 146           continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               itnm = itnm + 1

            else
                  goto 997
            end if

            if( jpnd .eq. 3 ) return

            close(iodw)

                  goto 240


*-----------------------------------------------------------------------
*        [t-heat]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 7) .eq. '[t-heat]' ) then

               if( chcm(i1+8:i1+10) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call theat(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-star]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 7) .eq. '[t-star]' ) then

               if( chcm(i1+8:i1+10) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               itnam=1  ! T.Sato 2018/2/15, [t-star] or [t-interact]

               call tstar(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,itnam,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-interact] ! Equivalent to [t-star], T.Sato 2018/2/15
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+11) .eq. '[t-interact]' ) then

               if( chcm(i1+12:i1+14) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               itnam=2  ! T.Sato 2018/2/15, [t-star] or [t-interact]

               call tstar(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,itnam,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-time]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 7) .eq. '[t-time]' ) then

               if( chcm(i1+8:i1+10) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call ttime(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-dpa]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 6) .eq. '[t-dpa]' ) then

               if( chcm(i1+7:i1+9) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call tdpa(jsn,jsi,dsin,idsi,ill,ilf,
     &                   jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-product]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 10) .eq. '[t-product]' ) then

               if( chcm(i1+11:i1+13) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call tprodct(jsn,jsi,dsin,idsi,ill,ilf,
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-gshow]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 8) .eq. '[t-gshow]' ) then

               if( chcm(i1+9:i1+11) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call tgshow(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-rshow]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 8) .eq. '[t-rshow]' ) then

               if( chcm(i1+9:i1+11) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call trshow(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-3dshow]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 9) .eq. '[t-3dshow]' ) then

               if( chcm(i1+10:i1+12) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call tdshow(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

*-----------------------------------------------------------------------
*        [t-userdefined]
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 14) .eq. '[t-userdefined]' ) then

               if( chcm(i1+15:i1+17) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call tusrdf(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

            end if

           end if ! for icntl.ne.14,15 (2017.8.11)

*-----------------------------------------------------------------------
*        [t-volume]
*-----------------------------------------------------------------------
           if ( mstz(1).eq.14 ) then ! for icntl.eq.14

            if( chcm(i1:i1+ 9) .eq. '[t-volume]' ) then

               if( chcm(i1+10:i1+12) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call tvolume(jsn,jsi,dsin,idsi,ill,ilf,
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

c skip reading the other tally sections here (S.H. revised 2017.8.11)
            else if( chcm(i1:i1+ 2) .eq. '[t-' ) then

                  jpn = 2
                  goto 140

            end if

           end if

*-----------------------------------------------------------------------
*        [t-wwbg]
*-----------------------------------------------------------------------
           if ( mstz(1).eq.15 ) then ! for icntl.eq.15

            if( chcm(i1:i1+ 7) .eq. '[t-wwbg]' ) then

               if( chcm(i1+8:i1+10) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call twwbg(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  if( ierr .ne. 0 ) goto 999

                  if( jpn .eq. 3 ) return

                  goto 240

c skip reading the other tally sections here (S.H. revised 2017.8.11)
            else if( chcm(i1:i1+ 2) .eq. '[t-' ) then

                  jpn = 2
                  goto 140

            end if

           end if

*-----------------------------------------------------------------------
*        [end]
*-----------------------------------------------------------------------

            if( chcm(i1:i1+ 4) .eq. '[end]' ) then ! S.H. revised (2017.8.11)

               if( chcm(i1+5:i1+7) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

                  return

*-----------------------------------------------------------------------
*        [old] : for old tally input
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 4) .eq. '[old]' ) then

               iot = 6
               ict = 1

               call read01(jsi,iot,ict,ierr)

               ivers = 0

               return

*-----------------------------------------------------------------------
*        [oldgo] : for old input and go
*-----------------------------------------------------------------------

            else if( chcm(i1:i1+ 6) .eq. '[oldgo]' ) then

               iot = 6
               ict = 0

               call read01(jsi,iot,ict,ierr)

               ivers = 1

               return

*-----------------------------------------------------------------------
*        else
*-----------------------------------------------------------------------

            else

                  if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) goto 997

                  goto 140

            end if

*-----------------------------------------------------------------------
*     end of error
*-----------------------------------------------------------------------
  951 continue

            m_err = 'Double definition of [T-DCHAIN] is not allowed.'
            ErrCha = ''
            ErrID = 'L:2237/R:read00/F:read00.f'
            l_err = ill(jsn)
            k_err = jsn

            goto 999
*-----------------------------------------------------------------------
  995 continue

            m_err = '[parameters] section is missing.'
            ErrCha = ''
            ErrID = 'L:2247/R:read00/F:read00.f'
            l_err = ill(jsn)
            k_err = jsn

            goto 999
*-----------------------------------------------------------------------
c error when [parameters] are read two times (S.H. added 2017.8.14)

  996 continue

            m_err = '[parameters] section appears twice.'
            ErrCha = ''
            ErrID = 'L:2259/R:read00/F:read00.f'
            l_err = ill(jsn)
            k_err = jsn

            goto 999

*-----------------------------------------------------------------------

  997 continue

            m_err = 'Unknown section name'
            ErrCha = ''
            ErrID = 'L:2271/R:read00/F:read00.f'
            l_err = ill(jsn)
            k_err = jsn

            goto 999

*-----------------------------------------------------------------------

  998 continue

            dsin(0) = 'Error Line'
            idsi(0) = 12
            m_err = 'There is nothing in the normal input'
            ErrCha = ''
            ErrID = 'L:2285/R:read00/F:read00.f'
            l_err = 1
            k_err = 0

            goto 999

*-----------------------------------------------------------------------

  999 continue

         write(*,'(/" ***** Error Message from Input File *****"/)')

C ----- set char parameter ----------------
         call ErrLine_Adjust(dsin(k_err),l_err)
         write(*,*) trim(dsin(k_err)),l_err,':'
C -----------------------------------------


            icf = 200

         do 910 i = 200, 1, -1

            if( m_err(i:i) .ne. ' ' ) goto 911

  910    continue

  911       icf = i

         call ErrWrite(ErrID, ErrCha)
         write(*,'(" error = ",200A1/)') ( m_err(i:i), i=1,icf )

            ierr = 1

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine setbnk(io,jo,ierr)
*                                                                      *
*       set bank memory array                                          *
*       modified by K.Niita on 2002/02/05                              *
*                                                                      *
************************************************************************
      use MMBANKMOD,only:maxbnk_bank !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /bnkmem/ maxbnk, maxbn2, rtrckflp

      common /cggmm/  ngstar, ngfini, ngfin0
      common /talmm/  nmmax, lmmax, itlmx
      common /bnkmm/  mbmax, mbfin, mbtfin
      common /matmm/  nmhigh, nminth, nmfinh

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall19/ itsmn(itlmax), itstm(itlmax)

*-----------------------------------------------------------------------

            ierr = 0

*-----------------------------------------------------------------------

            mnax = nmmax * 2
            maxbnk_bank=maxbnk !FURUTA

            ibkblz = 0 !FURUTA mnax
            ibknam = 0 !FURUTA ibkblz + maxbnk
            ibknmd = 0 !FURUTA ibknam + maxbnk
            ibknty = 0 !FURUTA ibknmd + maxbnk
            ibknkf = 0 !FURUTA ibknty + maxbnk
            ibknfc = 0 !FURUTA ibknkf + maxbnk
            ibknct = 0 !FURUTA ibknfc + maxbnk
            ibksos = 0 !FURUTA ibknct + maxbn2
            ibkzst = 0 !FURUTA ibksos + maxbnk * 3

            mnax = ibkzst + maxbnk
            mnax = ( mnax + mod(mnax,2) ) / 2

            ibkwt  = 0 !FURUTA mnax
            ibku   = 0 !FURUTA ibkwt  + maxbnk
            ibkv   = 0 !FURUTA ibku   + maxbnk
            ibkw   = 0 !FURUTA ibkv   + maxbnk
            ibke   = 0 !FURUTA ibkw   + maxbnk
            ibkec  = 0 !FURUTA ibke   + maxbnk
            ibkx   = 0 !FURUTA ibkec  + maxbnk
            ibkxc  = 0 !FURUTA ibkx   + maxbnk
            ibky   = 0 !FURUTA ibkxc  + maxbnk
            ibkyc  = 0 !FURUTA ibky   + maxbnk
            ibkz   = 0 !FURUTA ibkyc  + maxbnk
            ibkzc  = 0 !FURUTA ibkz   + maxbnk
            ibkt   = 0 !FURUTA ibkzc  + maxbnk
            ibktc  = 0 !FURUTA ibkt   + maxbnk
            ibkwin = 0 !FURUTA ibktc  + maxbnk
            ibkwnz = 0 !FURUTA ibkwin + maxbnk
            ibkxfc = 0 !FURUTA ibkwnz + maxbnk
            ibkspx = 0 !FURUTA ibkxfc + maxbnk
            ibkspy = 0 !FURUTA ibkspx + maxbnk
            ibkspz = 0 !FURUTA ibkspy + maxbnk

            mnax = ibkspz + maxbnk
            mnax = mnax * 2

            iakblz = 0 !FURUTA mnax
            iaknam = 0 !FURUTA iakblz + maxbn2
            iaknmd = 0 !FURUTA iaknam + maxbn2
            iaknty = 0 !FURUTA iaknmd + maxbn2
            iaknkf = 0 !FURUTA iaknty + maxbn2
            iaknfc = 0 !FURUTA iaknkf + maxbn2
            iaknct = 0 !FURUTA iaknfc + maxbn2
            iaksos = 0 !FURUTA iaknct + maxbn2
            iakzst = 0 !FURUTA iaksos + maxbn2 * 3

            mnax = iakzst + maxbn2
            mnax = ( mnax + mod(mnax,2) ) / 2

            iakwt  = 0 !FURUTA mnax
            iaku   = 0 !FURUTA iakwt  + maxbn2
            iakv   = 0 !FURUTA iaku   + maxbn2
            iakw   = 0 !FURUTA iakv   + maxbn2
            iake   = 0 !FURUTA iakw   + maxbn2
            iakx   = 0 !FURUTA iake   + maxbn2
            iaky   = 0 !FURUTA iakx   + maxbn2
            iakz   = 0 !FURUTA iaky   + maxbn2
            iakt   = 0 !FURUTA iakz   + maxbn2
            iakwin = 0 !FURUTA iakt   + maxbn2
            iakwnz = 0 !FURUTA iakwin + maxbn2
            iakxfc = 0 !FURUTA iakwnz + maxbn2
            iakspx = 0 !FURUTA iakxfc + maxbn2
            iakspy = 0 !FURUTA iakspx + maxbn2
            iakspz = 0 !FURUTA iakspy + maxbn2
            ibtetpos = 0
            ibtetposa= 0

            mnax  = iakspz !+ maxbn2
            mbfin = mnax
            mbmax = mnax - nmmax

         if( mnax .gt. mmmax ) then

            mmmax = mnax
            itlmx = 0

         else

            mmmax = mmmax
            itlmx = itlmx - mbmax

         end if

            mbtfin = mmmax

            mmmax = mmmax + 1


*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine setpar(io,jo,ivers,ierr)
*                                                                      *
*       set parameters and initial values                              *
*       last modified by K.Niita and S. Hashimoto on 2011/08/02        *
*                                                                      *
************************************************************************
!$    use omp_lib
      use dedx_file
      use MEMBANKMOD,only:mxmat_bank,maxnl_bank,mxnel_bank,msumnel_bank,
     &     mdbatima,dbcutoff,mdbpseud !FURUTA20160126
      use TETRAMOD,only:itetra,ntetsurf,ntetelem !FURUTA20160908
      use neutrino_mod,only:ntrnore
      use NGSDATAMOD, only : gstbl, weitn, fstbl
      use ELEDATAMOD, only : ion_potential_table, allocate_frac_ichem,
     & ichem, set_atomic_relaxation_table,
     & electron_kinetic_energy_table

      use NDATA2MOD
      use levdenmod, only : sstbl
      use fission_mod, only : ifiss
      use moddas
      use moddas_character
      use moddas_material
      use moddas_tally
      use liboutmod, only: ilibpart, elibreadmax ! frtati 2022/12/28
      use ets_art, only : etsart_db ! hirata 2023/10/31

      use t4dtrack_mod, only: nt4domp

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param01.inc'

      include 'err.inc'

*-----------------------------------------------------------------------

      common /mpi00/ npe, me
      common /mpi03/ filhd, nfihd
      character filhd*100

*-----------------------------------------------------------------------

      common /irndm/ irndmode,idmprijk
      common /ccggg/  icgg
      common /ggmes/  iggcm
      common /bnkmem/ maxbnk, maxbn2, rtrckflp
      common /inggs/  iog, igcel, ioa, igsuf, iob, igtrs
      common /celdb/  idsn(kvlmax), idtn(kvlmax)
      common /celda/  deng(kvlmax)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)

      common /celdg/  rhog(kvlmax)

      common /inpec/  ititl, ipara, ibody, iregn, llarr, itby, itar
      common /htitl/  iclgt(100), ctitl(100)
      character       ctitl*200

      common /regdm/  idmg(kvlmax)
      common /regda/  ichl(kvlmax), chsm(kvlmax), ichmx, iod
      character       chsm*10
      common /regde/  ichp(kvlmax), ilat(kvlmax), idct(kvlmax)
      common /regdb/  nrsq, irsq(10)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /regdd/  ivolm, iimpo
      common /verjam/ versn, lastr, iyeav, imonv, idayv
      common /startf/ iday0,imon0,iyer0,ihor0,imin0,isec0
      common /parai/  ipsq(400)
      common /paraj/  mstz(300), parz(300)
      common /parak/  icnu, icdf(400), icdl(400), chnm(400)
      character       chnm*8
      common /paran/  icfn(100), ilfn(100), chfn(100)
      character       chfn*200

      logical   exex

      common /delesg/ rdels(kvlmax), mndel, ndels(kvlmax)
      common /delreg/ delm(kvlmax), kdelt

! T.Sato for Track structure simulation
      common /tscmsg/ ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /tscreg/ ntscell(kvlmax)  ! ntcscell=0:no track structure, 1:water

      common /volmsg/ rvols(kvlmax), mnvol, nvols(kvlmax)
      common /volreg/ dvol(kvlmax)

      common /tmpmsg/ rtmps(kvlmax), mntmp, ntmps(kvlmax)
      common /tmpreg/ dtmp(kvlmax)
      common /pwtmsg/ rpwts(kvlmax), mnpwt, npwts(kvlmax)
      common /pwtreg/ dpwt(kvlmax)

      common /impreg/ dimp(kvlmax)
      common /wtcntl/ iwt, icimp(20), ifcls(20), iwwin(20), ircls(20)

      common /brsmsg/ mnbrs, icbrs, mbbrs(kvlmax), cbrem(49)
      common /regdn/  irden

*-----------------------------------------------------------------------

      common /mttmc/  smttc(kvlmax), mttcn, mttc1(kvlmax), mttc2(kvlmax)
      common /mtnmc/  smtnc(kvlmax), dmtnc(kvlmax,2),
     &                mtncn, mtnc(kvlmax,2), nmtnc(kvlmax,2)
      character dmtnc*80

      common /mtnmcl/ dmhsb(-1:kvlmax,4),
     &                dmtnm(-1:kvlmax), dmtcl(-1:kvlmax),
     &                nmtnm(-1:kvlmax), nmtcl(-1:kvlmax)
      character dmtnm*80, dmtcl*30

      common /elreg/  elarg(kvlmax,4), mlrgn, melrg(kvlmax,2), ielusr

      common /magreg/ nmreg, ingrc, ingrt, kmags, magtin, imgusr
      common /elmgrg/ nelctf, nereg, inerc, inert, kelcs

      parameter ( icnm = 36 )
      parameter ( iccx = 15 )
      parameter ( imfnmax=1000) ! T.Sato 2022/10/30

      common /coldef/ coldf(3,icnm), lcoln(icnm), colnm(icnm), icolr(30)
      character colnm*15

*-----------------------------------------------------------------------

      common /rcmini/ rcmin(20)
      common /ndemax/ dnmax(20)
      common /dpnmaxcom/ dpnmax ! frtati 2021/12/17

*-----------------------------------------------------------------------

      common /wwindp/ wupn, wsurvn, mxspln, mwhere, mvoww

*-----------------------------------------------------------------------

c S.H. xorshift (2020.2.6)
      parameter (period_LCG=2d0**48)
      common /randnbit/ bitrseed
      common /randn/ nrandgen
      common /iradkk/ randkk,irskip
      common /eparm/  esmax, esmin, emin(20)
      common /tparm/  tmax(20)
      common /wparm/  swtm(20), wc1(20), wc2(20) !FURUTA

      common /kmat1a/ mxmat, mxmat0, mxnel
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

      common /kmat1g/ kmat(kvlmax)
      common /kmat1h/ kmatg(kvlmax)
      common /kmat1i/ kmatc(kvlmax)
      common /kmat1j/ intum

      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax)
      common /kmat1ka/ kmathd(kvlmax), kmathe(kvlmax)

      common /kmat1o/ iom1, iom2, iom3
      common /kmat1p/ kmout
      common /kmatad/ matadd
      common /matmm/  nmhigh, nminth, nmfinh
      common /xgeosm/ ksig(kvlmax)
      common /xinels/ ksige, ksign
      common /sparcmn/ khb, kh1, kh2, kh3, khe(kvlmax), khp,
     &                 kdf, kdg, kro
      common /argcns/ kcar, kczs, kcze

      common /natnuc/ natnn(maxpt), natnm(maxpt,10), patnn(maxpt,10)
      common /geosig/ geosig(250)

      common /cparm/  maxbch,maxcas
      common /bparm/  andt,jevap,npidk
      common /ngcut/  incut, igcut, ipcut
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx
      common /iaoru/  iflagAorU !FURUTA20200219 E/A <-> E/u
      common /spsgn/  nspsgn
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /emode/  emodem, ge1, ge2, iemode
      common /cugnon/ icugn
      common /tcntl/  icntl, inucr
      common /engch/  ejamnu, ejampi, eisobar, eqmdnu, eqmdmn, ejamqmd
      common /taliin/ rsouin, nzztin, nrgnin
      common /talout/ itall
      common /talsav/ iptall
      common /cgerr/  nlost, ilost, igerr, icger, ncger, nrecover
      common /cgstr/  novp, nrovp(3,1000)
      common /clionprd/  lionprd
      common /gravit/ grav(3), igrav
      common /crshi/  bplus, icrhi, ijudg, imadj, iqmax
      common /fsmul/  ridfs, iidfs
      common /sabad/  isaba
      common /sansd/  isans
      common /voxel/  ivoxel
      common /pnint/  ipnint
      common /gmuppd/ igmuppd ! S.Abe 2019/11/07
      common /cincl/  inclg, inclv
      common /cincelf/incelf
      common /ceinc/  einclmin, einclmax, eielfmin, eielfmax
      common /cidwba/ idwba
      common /smmflg/ ismm, ifbm
      common /ckurotama/  dsck
      common /ccxsm/  icxsni, icxspi
      common /qmdflg/ irqmd
      common /gem/    ngem ! 2017/4/28 Ogawa. GEM version
      common /muint/ imuint, imubrm, imuppd, imucap
      common /emurng/ emumin, emumax
      common /prmui/ prmui1
      integer italsh
      common /talsh/ italsh
      data italsh / 0 /
      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall30/ itnda(itlmax)

      common /tall40/ rtorg(itlmax,3), rteye(itlmax,3), rtlit(itlmax,3),
     &                rtwin(itlmax,3), rtbrt(itlmax,2), ithvn(itlmax),
     &                itwin(itlmax,2), itbox(itlmax), itmir(itlmax),
     &                rtbox(itlmax,5,10), rtout(itlmax), rthet(itlmax),
     &                itlin(itlmax), itshd(itlmax), itgxs(itlmax)

      common /tall41/ rdmax

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /ivmth/  icnt14
      data icnt14 / 0 /

      common /tall69/ itnms(itlmax), trmsh(itlmax,10),
     &                tzmsh(itlmax,10), tfmsh(itlmax,10),
     &                trmpo(itlmax,10), tzmpo(itlmax,10)
      common /isorsc/ isort(isrc,4), rsort(isrc,13)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall68/ itstp(itlmax), itmth(itlmax), rtvr0(itlmax),
     &                rtvx0(itlmax), rtvy0(itlmax), rtvz0(itlmax),
     &                rtvx1(itlmax), rtvy1(itlmax), rtvz1(itlmax)
      data eps / 0.0000000001234d0 /

      common /tall79/ itallech

      common /gsline/ nowgshow, igsline

      common /isomul/ smlwt(isrc), totfact, imsrc
      common /isorst/ jstyp(isrc), istyp(isrc), inkf0(isrc), lstyp(isrc)
!$OMP THREADPRIVATE(/isorst/)
      common /isorori/ jstypori(isrc) ! T.Sato, original jstyp written in input file
      common /isorsp/ sx0(isrc), sy0(isrc), sz0(isrc), sx1(isrc),
     &                sy1(isrc), sz1(isrc), sr0(isrc), se0(isrc),
     &                sdir(isrc), srx(isrc), sry(isrc), swem(isrc),
     &                sphi(isrc), sdom(isrc), swt0(isrc)
      common /isorsn/ sr1(isrc), sr2(isrc), isrn(isrc)
      common /isorsm/ stm0(isrc),stmw(isrc),stmc(isrc),stmd(isrc),
     &                jttyp(isrc),jttpn(isrc)
      common /isorsa/ sfactor(isrc)

      common /fmcard/ ifm, ifmi(imfnmax,3)


      common /electh/ ieleh
      common /dircha/ idirch
      common /rngmem/ nrange, krnge, krngn !FURUTA
      common /estrag/ kesta, kestz, kestd
      common /atima01/ katima
      common /dumpall/ idumpall
      common /dcyd/ rtdcy(50), kfdcy(50), ndcy

      common /infprint/ infout

      common /adjoint/ iadjnt
      common /adjenrg/ adjemax

      common /timecut/ timeout

*-----------------------------------------------------------------------
      common /istcut/ ist_cut, ist_bat

*-----------------------------------------------------------------------
      common /nwwbias/ iwwbias

*-----------------------------------------------------------------------
      common /expnatur/ inatur

*-----------------------------------------------------------------------
      common /ndtmax/ indmp, indmn, indmu, indmd, indma,
     & nucdxp(500), nucdxn(500), nucdxu(500), nucdxd(500), nucdxa(500),
     & matdxp(500), matdxn(500), matdxu(500), matdxd(500), matdxa(500)
      common /ddtmax/ dmxdxp(500), dmxdxn(500),
     &                dmxdxu(500), dmxdxd(500), dmxdxa(500)

      common /ggcell/ icells, iobo

*-----------------------------------------------------------------------
      common /egsemi/ iegsemi, iegsout, ipegs

*-----------------------------------------------------------------------

! T.Sato, 2015/03/10, introduce free parameter for Moliere 1st
      common /aspara/aspara1,aspara2

cMIURA 2016.09.30 tally output unit, the energy per nucleon. [MeV/n]
      integer           iMeVperu
      common /cMeVperu/ iMeVperu

! T.Sato 2017/05/11, Track structure mode
      common /etsminmax/ etsmin,etsmax
      common / tsxcl / itsxcl
      common / ruthang/ rathe1, rathe2

! T.Sato 2018/03/06, Kerma control
      common /kermon/ikerman,ikermap

      common /fixcharge/ifixchg     !T.Sato 2019/02/17

      character chme*5

      logical deqn1

      dimension ix(3)
      dimension nbbrs(kvlmax)

      dimension tcol(3)
      character c1*1, c2*1
      character lum(200)*1

      dimension bval(50)
      dimension ipva(4)
      character chbd*3
      character chtit*60

*-----------------------------------------------------------------------

      character filnm*100

      character mtm*6

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

      common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------
cFURUTA20150714 TETRA
      integer nlat3,itetcl(10),itetvol,itetauto,itgchk
      common /tetf1/ nlat3,itetcl,itetvol,itetauto,itgchk
*-----------------------------------------------------------------------
      common /pnmul/ pnimul
      common /gmmul/ gmumul ! S.Abe 2019/11/07

*-----------------------------------------------------------------------
      common /celepcc/ enumpcc(kvlmax)

*-----------------------------------------------------------------------
      common /reslet/ irlet

*-----------------------------------------------------------------------
      common /multipl/ imltp,iimlt(multmax),inmlt(multmax),
C MATSUDA 2024.11.25 (lagrange: ilmlt)
     &  ilmlt(multmax),
     &  idmlt(multmax),ismlt(multmax),impan(multmax),impat(multmax,6,3),
     &                 jmpat(multmax,6,6,2), imdfl(multmax)

      character idcar*3  ! T.Sato 2018/01/31
      character mltfl*200

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      common /scincom/iscinful  ! use SCINFUL or not
      common /scinswt/iswitch
!$OMP THREADPRIVATE(/scinswt/)

      common / tsminmax  / tsmax
      common /etsart/ bgets(kvlmax), ebgets(kvlmax),
     &                wvets(kvlmax), ewvets(kvlmax)
      common /qmdscm/ qmdscm_h0, qmdscm_d, qmdscm_rcls, iqmdscm

      common /geomemcom/ igeomem  ! T.Sato 2024/12/25

*-----------------------------------------------------------------------
cKN 2024/03/26

      common /bnkmew/ iwwbnk

*-----------------------------------------------------------------------
cKN 2018/01/30

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      character jdcar*5
      character element(104)*2

      data element/
     & 'H_','He','Li','Be','B_','C_','N_','O_',
     & 'F_','Ne','Na','Mg','Al','Si','P_','S_',
     & 'Cl','Ar','K_','Ca','Sc','Ti','V_','Cr',
     & 'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge',
     & 'As','Se','Br','Kr','Rb','Sr','Y_','Zr',
     & 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',
     & 'In','Sn','Sb','Te','I_','Xe','Cs','Ba',
     & 'La','Ce','Pr','Nd','Pm','Sm','Eu','Gd',
     & 'Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf',
     & 'Ta','W_','Re','Os','Ir','Pt','Au','Hg',
     & 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra',
     & 'Ac','Th','Pa','U_','Np','Pu','Am','Cm',
     & 'Bk','Cf','Es','Fm','Md','No','Lr','Ku'/

      dimension remsh(2)

      common /trwwbg/ itrwwbg(4),rtrwwbg(13)

cKN 2018/02/12
*-----------------------------------------------------------------------

            ierr  = 0

      ifixchg=mstz(132) ! T.Sato 2019/02/18
      iscinful=mstz(141) ! T.Sato 2020/02/16
      tsmax=parz(196)   ! T.Sato 2020/03/18
      igeomem=mstz(167)  ! T.Sato 2024/12/25
      KBCflg = 0
      if(tsmax.gt.1.0d-3) call KURBUCcheck(KBCflg)

      iswitch=0
      if(iscinful.ne.0) call CX_init ! Initializaiton of SCINFUL cross section, T.Sato 2020/02/16

*-----------------------------------------------------------------------
*     open output data file ( unit = 28 )
*-----------------------------------------------------------------------

             if( me .eq. 0 ) then

                  io = 28

               if( mstz(82) .ne. -1 ) then   ! infout = mstz( 82 )

                  open(io, file = chfn(6), status = 'unknown' )

               else

                  open(io, form='formatted',status='scratch')

               end if

             end if

*-----------------------------------------------------------------------
*        initialize nuclear structure table
*-----------------------------------------------------------------------

               call gstbl ! ground state table
               call sstbl ! statistical nuclear strucutre table
               call fstbl(ierr) ! fission barrier table
               call igamma1set ! discrete level data table FURUTA20201007
! igmma1set required almost all cases (igamma!=0 and jevap!=0)

               call ion_potential_table(ierr) ! ionization potential table used by track structure and EBITEM
               call set_atomic_relaxation_table(ierr)
               if(mntsc.ne.0) call electron_kinetic_energy_table ! if track structure is on, call kinetic energy table
               if(.not.allocated(ichem)) call allocate_frac_ichem ! material composition array

*-----------------------------------------------------------------------
*     minimum and maxmum energy for energy loss
*-----------------------------------------------------------------------

                  esmax = parz(170)
                  esmin = parz(162)
                  ! 2023/9/15 Ogawa. ITSART uses ATIMA DEDX at 3 MEV. If esmin too high, ITSART crashes.
                  ! 2024/4/8  Ogawa. ITSART infinite loop if emin is below valid cross section
                  if(mntsc.ne.0) esmin = max(1.d-5, min(2.9d0, esmin))

*-----------------------------------------------------------------------
*     set parameters from input and check the values
*-----------------------------------------------------------------------

! T.Sato 2024/02/18, check nfcseg parameter
      if(mstz(162).le.0) then
       ErrCha = 'Warning: nfcseg should be greater than 0'//
     & ' and it is changed to 1'
       MsgID = 'L:2953/R:setpar/F:read00.f'
       call ErrWrite(MsgID, ErrCha)
       mstz(162)=1
      elseif(mstz(162).gt.1000) then
       ErrCha = 'Warning: nfcseg should be less than 1000'//
     & ' and it is changed to 1000'
       MsgID = 'L:2959/R:setpar/F:read00.f'
       call ErrWrite(MsgID, ErrCha)
       mstz(162)=1000
      endif

c S.H. xorshift (2020.2.6)
      nrandgen = mstz(139)
      itimrand = mstz(140)

      if( nrandgen .eq. 0 ) then ! Linear Congruential Generator
         if ( parz(195) .eq. 0d0 ) then
            bitrseed = 6647299061401.d0 ! default value of bitrseed
         else
            bitrseed = parz(195)
            write(*,'(a)')
     &           '*** Warning: bitrseed is not used when nrandgen = 0'
         end if

         if ( itimrand .eq. 0 ) then
            randkk = dabs( parz(21) )
            if ( randkk .ge. period_LCG ) then
               ErrCha = ''
               MsgID = 'L:2981/R:setpar/F:read00.f'
               call ErrWrite(MsgID, ErrCha)
               write(jo,'(" **** Error : absolute value of rseed "
     &,1p1e25.16e3," must be less than "
     &,1p1e16.8," when nrandgen = 0")') randkk, period_LCG
               write(io,'(" **** Error : absolute value of rseed "
     &,1p1e25.16e3," must be less than "
     &,1p1e16.8," when nrandgen = 0")') randkk, period_LCG
               goto 999
            end if
         else                   ! time dependent initial random number
            iseed = isec0*10000 + ihor0*100 + imin0
            if(iseed / 2 * 2 .eq. iseed ) iseed = iseed + 101
            randkk = dble( iseed )
            if ( parz(21) .ne. 0d0 )
     &           write(*,'(a)')
     &      '*** Warning: specified rseed is not used when itimrand = 1'
         end if

      else ! xorshift
         if ( itimrand.eq.0 ) then
            if ( parz(195) .eq. 0d0 ) then
               if ( parz(21) .eq. 0d0 ) then
                  bitrseed = 6647299061401.d0 ! default value of bitrseed
               else
                  bitrseed = parz(21) ! rseed is used
               end if
               randkk = bitrseed
            else
               bitrseed = parz(195)
               if ( parz(21) .eq. 0d0 ) then
                  randkk = bitrseed
               else
                  randkk = parz(21)
                  write(*,'(a,a)')
     &                 '*** Warning: specified rseed is not used ',
     &                 'when nrandgen=1 and bitrseed is specified'
               end if
            end if

         else ! time dependent initial random number
            iseed = isec0*10000 + ihor0*100 + imin0
            if(iseed / 2 * 2 .eq. iseed ) iseed = iseed + 101
            bitrseed = dble( iseed )
            if ( parz(195) .eq. 0d0 ) then
               if ( parz(21) .eq. 0d0 ) then
                  randkk = bitrseed
               else
                  randkk = parz(21)
                  write(*,'(a)')
     &      '*** Warning: specified rseed is not used when itimrand = 1'
               end if
            else
               write(*,'(a)')
     &   '*** Warning: specified bitrseed is not used when itimrand = 1'
               if ( parz(21) .eq. 0d0 ) then
                  randkk = bitrseed
               else
                  randkk = parz(21)
                  write(*,'(a)')
     &      '*** Warning: specified rseed is not used when itimrand = 1'
               end if
            end if
         end if
      end if

                  irskip = mstz( 2 )

*-----------------------------------------------------------------------
*           cutoff energy
*-----------------------------------------------------------------------

               do i = 1, 20

                  emin(i) = parz(i)

               end do

*-----------------------------------------------------------------------
*           Track structure cutoffenergy  T.Sato 2017/05/11
*-----------------------------------------------------------------------
               etsmin=parz(191)
               etsmax=parz(192)
               if(etsmax.lt.1.0d-3) then
                write(ErrCha,'("etsmax should be greater than 1 keV,",
     &          " and automatically set to etsmax = 1.0d-3")')
                ErrID = 'L:3067/R:setpar/F:read00.f' !W04_008_001
                call ErrWrite(ErrID,ErrCha)
                etsmax = 1.0d-3
               endif

*-----------------------------------------------------------------------
*              check charged particle minimum energy
*-----------------------------------------------------------------------
cFURUTA20160422 Warning message when emin changed

                  if( emin(1) .lt. esmin )then

                   ErrCha = ''
                   ErrID = 'L:3080/R:setpar/F:read00.f' !W04_008_002
                   call ErrWrite(ErrID,ErrCha)

                   emin(1) = esmin
                   write(*,'(a,i2,a,1p1d14.7,a)')
     &                  '*** Warning: emin(',1,') is adjusted to',
     &                  emin(1),
     &                  ' by esmin'
                  endif
               do i = 3, 10

                  if( ichgf(i,0) .ne. 0 ) then

                     if( emin(i) .lt. esmin ) then

                      ErrCha = ''
                      ErrID = 'L:3096/R:setpar/F:read00.f' !W04_008_003
                      call ErrWrite(ErrID,ErrCha)

                      emin(i) = esmin
                      write(*,'(a,i2,a,1p1d14.7,a)')
     &                  '*** Warning: emin(',i,') is adjusted to',
     &                  emin(i),
     &                  ' by esmin'
                     endif

                  end if

               end do

               if( emin(11) .lt. esmin * 2.0 )then

                ErrCha = ''
                ErrID = 'L:3113/R:setpar/F:read00.f' !W04_008_004
                call ErrWrite(ErrID,ErrCha)

                emin(11) = esmin * 2.0
                write(*,'(a,i2,a,1p1d14.7,a)')
     &               '*** Warning: emin(',11,') is adjusted to',
     &               emin(11),
     &               ' by esmin'
               endif


! Check cut-off energy of electron, positron, and photon, T.Sato 2018/07/06
       if(emin(12).lt.1.0d-3) then
        ErrCha = ''
        MsgID = 'L:3127/R:setpar/F:read00.f'
        call ErrWrite(MsgID, ErrCha)
        write(jo,'("emin(12)=",es12.3,
     &  " is too small, and changed to 1.0E-3")') emin(12)
        emin(12)=1.0d-3
       endif
       if(emin(13).lt.1.0d-3) then
        ErrCha = ''
        MsgID = 'L:3135/R:setpar/F:read00.f'
        call ErrWrite(MsgID, ErrCha)
        write(jo,'("emin(13)=",es12.3,
     &  " is too small, and changed to 1.0E-3")') emin(13)
        emin(13)=1.0d-3
       endif
       if(emin(14).lt.1.0d-3) then
        ErrCha = ''
        MsgID = 'L:3143/R:setpar/F:read00.f'
        call ErrWrite(MsgID, ErrCha)
        write(jo,'("emin(14)=",es12.3,
     &  " is too small, and changed to 1.0E-3")') emin(14)
        emin(14)=1.0d-3
       endif


cFURUTA20160422---------------------------------------------------------
               do i=15,19
                if(emin(i).lt.esmin)then

                 ErrCha = ''
                 ErrID = 'L:3156/R:setpar/F:read00.f' !W04_008_005
                 call ErrWrite(ErrID,ErrCha)

                 emin(i)=esmin
                 write(*,'(a,i2,a,1p1d14.7,a)')
     &                '*** Warning: emin(',i,') is adjusted to',
     &                emin(i),
     &                ' by esmin'
                endif
               enddo
c-----------------------------------------------------------------------

               do i = 1, 20

                  parz(i) = emin(i)

               end do

*-----------------------------------------------------------------------
*           reaction cutoff energy
*-----------------------------------------------------------------------

               do i = 1, 20

                  rcmin(i) = parz(i)

                  if(i.ne.2.and.i.ne.12.and.i.ne.13.and.i.ne.14
     &            .and.i.ne.15 )  ! T.Sato 2024/05/25, deuteron cmin should be below 1 MeV/n
     &            rcmin(i)=max(1.0d0,rcmin(i)) ! should be above 1 MeV, T.Sato 2019/08/01

                  k = 110 + i

                  do l = 1, ipara

                     if( k .eq. -ipsq(l) ) rcmin(i) = parz(k)

                  end do

               end do

                  rcmin(13) = rcmin(12)

*-----------------------------------------------------------------------
*           nuclear data maxmum energy
*-----------------------------------------------------------------------

               do i = 1, 20

                  dnmax(i) = parz(i)

! T.Sato 2021/10/03 For icntl = 1, dnmax is basically 0 except for neutrons and photons
                  if(mstz(1).eq.1.and.i.ne.2.and.i.ne.14) dnmax(i)=0.0

cKN 2016/08/09 I do not know why this is here
c                 if( i .eq. 2 ) dnmax(i) = parz(132)

                  k = 130 + i

                  do l = 1, ipara

                     if( k .eq. -ipsq(l) ) dnmax(i) = parz(k)

                  end do

               end do

                  dnmax(13) = dnmax(12)

                  dpnmax = parz(203) ! frtati 2021/12/17

*-----------------------------------------------------------------------
*           cutoff time and weight
*-----------------------------------------------------------------------

               do i = 1, 20

                  k = 28 + i
                  tmax(i) = parz(k)

                  k = 48 + i
                  swtm(i) = parz(k)

                  k = 68 + i
                  wc1(i) = parz(k)
                  if( wc1(i) .lt. 0.0 ) wc1(i) = -wc1(i) * swtm(i)

                  k = 88 + i
                  wc2(i) = parz(k)

               if( wc2(i) .le. -99.d0 ) then
                  wc2(i) = 0.5d0 * wc1(i)
               end if

                  do l = 1, ipara

                     if( k .eq. -ipsq(l) ) goto 500

                  end do

                  wc2(i) = wc1(i) * wc2(i)

                  goto 510

  500             continue

                  if( wc2(i) .lt. 0.0 ) wc2(i) = -wc2(i) * swtm(i)

  510             continue

                  if( wc2(i) .gt. wc1(i) ) then

                     write(io,'(" **** Error : in wc1 and wc2",
     &                          " definition:",
     &                          "  it should be wc1 > wc2.")')
                     write(io,'(" ityp =",i3,"  wc1 =",e13.5,
     &                          "  wc2 =",e13.5)')
     &                         i, wc1(i), wc2(i)

                     ErrCha = ''
                     MsgID = 'L:3275/R:setpar/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'(" **** Error : in wc1 and wc2",
     &                          " definition:",
     &                          "  it should be wc1 > wc2.")')
                     write(jo,'(" ityp =",i3,"  wc1 =",e13.5,
     &                          "  wc2 =",e13.5)')
     &                         i, wc1(i), wc2(i)

                     goto 999

                  end if

               end do

*-----------------------------------------------------------------------
*        weight window parameters
*-----------------------------------------------------------------------

                  wupn   = parz( 160 )
                  wsurvn = parz( 161 )
                  mxspln = mstz( 47 )
                  mwhere = mstz( 48 )
                  mvoww  = mstz( 83 )

                  iewwd = 0

               if( wupn .lt. 2.0d0 ) then

                  write(io,'(" **** Error : wupn should be larger ",
     &                       "than 2, in Weight Window parameter")')

                  ErrCha = ''
                  MsgID = 'L:3308/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(" **** Error : wupn should be larger ",
     &                       "than 2, in Weight Window parameter")')

                     iewwd = iewwd + 1

               end if

               if( wsurvn .gt. -100.d0 ) then

                  if( wsurvn .le. 1.d0 .or. wsurvn .gt. wupn ) then

                     write(io,'(" **** Error : wsurvn should be ",
     &                       "1 < wsurvn < wupn")')

                     ErrCha = ''
                     MsgID = 'L:3325/R:setpar/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'(" **** Error : wsurvn should be ",
     &                       "1 < wsurvn < wupn")')

                     iewwd = iewwd + 1

                  end if

               end if

               if( wsurvn .le. -100.d0 ) wsurvn = 0.6 * wupn

               if( mxspln .lt. 2 ) then

                  write(io,'(" **** Error : mxspln should be larger ",
     &                       "than 1, in Weight Window parameter")')

                  ErrCha = ''
                  MsgID = 'L:3344/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(" **** Error : mxspln should be larger ",
     &                       "than 1, in Weight Window parameter")')

                     iewwd = iewwd + 1

               end if

               if( mwhere .ne. 1 .and. mwhere .ne. 0 .and.
     &             mwhere .ne. -1 ) then

                  write(io,'(" **** Error : mwherre should be ",
     &                       "-1, 0, 1, in Weight Window parameter")')

                  ErrCha = ''
                  MsgID = 'L:3360/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(" **** Error : mwherre should be ",
     &                       "-1, 0, 1, in Weight Window parameter")')

                     iewwd = iewwd + 1

               end if

                  if( iewwd .gt. 0 ) goto 999

*-----------------------------------------------------------------------
cKN 2024/03/26
                  iwwbnk = mstz( 163 )

*-----------------------------------------------------------------------

                  maxbnk = mstz( 34 )
                  maxbn2 = maxbnk / 2
                  rtrckflp = parz(202)
                  maxcas = mstz( 3 )
                  maxbch = mstz( 4 )

*-----------------------------------------------------------------------

                  istdev = mstz( 73 ) !OBINATA(2012.6.13)

                  italsh = mstz( 144 )
                  if ( mstz(1).eq.13 .or. mstz(1).eq.17 ) italsh = 1 !S.H.(2020.11.5)

                  ireschk = mstz( 80 ) !T.Sato(2013.10.19)

                  if(istdev .lt. 0) then

                    istdev = -istdev

                  else if ( istdev .eq. 0 ) then
!$                  if( omp_get_max_threads() .gt. 1 .and.
!$   &                  italsh .eq. 1 ) then
!$                    istdev = 1
!$                  else
                      istdev = 2
!$                  end if
                  end if

*-----------------------------------------------------------------------

                  incut  = mstz( 17 )
                  igcut  = mstz( 18 )
                  ipcut  = mstz( 30 )

                  iggcm  = mstz( 32 )

                  npidk  = mstz( 5 )
                  lvlopt = mstz( 11 )
                  igamma = mstz( 12 )

                  andt   = parz( 22 )

                  ejamnu  = parz( 23 )
                  ejampi  = parz( 24 )
                  eisobar = parz( 25 )
                  eqmdnu  = parz( 110 )
                  eqmdmn  = parz( 164 )
                  ejamqmd = parz( 169 )

                  jevap  = mstz( 6 )
                  ntrnore = mstz( 131 )
                  nspred = iabs( mstz( 7 ) )
                  nspsgn = isign( 1, mstz( 7 ) )
! T.Sato, 2015/03/10, introduce free parameter for Moliere 1st
                  aspara1 = parz (183 )
                  aspara2 = parz (184 )

                  nwsprd = mstz( 14 )
                  nedisp = mstz( 51 )
                  itstep = mstz( 52 )
                  ndedx  = mstz( 58 )
cFURUTA20200219 E/A <-> E/u
                  if(ndedx.lt.0)then
                   iflagAorU=1
                   ndedx=abs(ndedx)
                  else
                   iflagAorU=0
                  endif
                  irlet  = mstz( 127 )

cFURUTA 20160126 ATIMA database by Wada
                  mdbatima = mstz( 115 )
                  dbcutoff = parz( 187 )

                  ielas = mstz( 8 )
                  icugn = mstz( 13 )

                  icntl = mstz( 1 )
*-----------------------------------------------------------------------
               if( icntl .eq. 5 .or. icntl .eq. 14 ) then
                  parz(27)  = 100000.d0
                  parz(163) = 100000.d0
               end if

*-----------------------------------------------------------------------

                  inucr = mstz( 15 )

                  itall  = mstz( 26 )
                  iptall = mstz( 33 )

                  if ( itall.eq.4 .and. istdev.eq.1 ) istdev = 2 !frtati 2021/03/06

                  icasc = mstz( 9 ) + 1

                  inclg = mstz( 71 )
                  inclv = mstz( 72 )
                  einclmin = parz ( 176 )
                  einclmax = parz ( 177 )

                  incelf = mstz( 76 )
                  eielfmin = parz ( 178 )
                  eielfmax = parz ( 179 )

                  idwba = mstz( 77 )

                  infout = mstz( 82 )

! T.Sato 2021/10/03 default of iMeVperU = 1 for icntl = 1
                 if(icntl.eq.1) then
                  imevperu = 1
                 else
                  imevperu = 0
                 endif
                 do kk = 1, ipara
                  if( ipsq(kk) .eq. 123   ) imevperu = mstz( 123 ) ! iMeVperu is changed only when explicitly defined
                 end do

               if( mstz(6) .eq. 0 ) then

                  jevap = 0

               else

                  jevap = 1

               end if


               if( mstz(6) .eq. 1 .and. mstz(10) .eq. 1 ) then

                  iqstep = 3

               else if( mstz(6) .eq. 1 ) then

                  iqstep = 5

               else if( mstz(6) .eq. 2 ) then

                  iqstep = 4

               else if( mstz(6) .eq. 3 ) then

                  iqstep = 6

               end if

               ngem = max(1, mstz(126)) ! 2017/4/28 Ogawa. GEM version. 0: recommended version(1), 1: Ver1, 2: Ver2

                  nlost = mstz( 21 )
                  ilost = 0

                  igerr = mstz( 22 )
                  nrecover = mstz(103)
                  icger = 0
                  ncger = 0

                  novp  = 0

                  do i = 1, 1000

                     nrovp(1,i) = 0
                     nrovp(2,i) = 0
                     nrovp(3,i) = 0

                  end do

                  matadd = mstz( 27 )
                  kmout  = mstz( 43 )
                  ieleh  = mstz( 46 )
                  idirch = mstz( 49 )

                  lionprd = mstz( 53 )

                  icrhi  = mstz( 54 )
                  bplus  = parz( 168 )
                  ijudg  = mstz( 56 )
                  imadj  = mstz( 57 )
                  iqmax  = mstz( 64 )

                  if( icrhi .gt. 3 .or. icrhi .lt. 0 ) icrhi = 1
                  if( ijudg .gt. 1 .or. ijudg .lt. 0 ) ijudg = 1
                  if( imadj .gt. 1 .or. imadj .lt. 0 ) imadj = 1

                  imgusr = mstz( 62 )
                  ielusr = mstz( 63 )

                  iidfs  = mstz( 66 )
                  isaba  = mstz( 67 )

                  ivoxel = mstz( 68 )
               if( ivoxel .lt. 0 .or. ivoxel .ge. 4 ) then !FURUTA20200515

                  write(io,'(/" Error : ivoxel should be",
     &            " 0,  1, 2 or 3.")')

                  ErrCha = ''
                  MsgID = 'L:3574/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/" Error : ivoxel should be",
     &            " 0,  1, 2 or 3.")')

                  goto 999

               end if

                  icells = mstz( 142 )
               if( icells .lt. 0 .or. icells .ge. 4 ) then

                  write(io,'(/" Error : icells should be",
     &            " 0,  1, 2 or 3.")')

                  ErrCha = ''
                  MsgID = 'L:3590/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/" Error : icells should be",
     &            " 0,  1, 2 or 3.")')

                  goto 999

               end if

                  nelctf = mstz( 70 )

                  ifbm = mstz( 81 ) ! T.Ogawa 2014/4/16
                  ismm = mstz( 74 )
                  if ( ismm .eq. 1 ) iqmax = 75

                  dsck = parz( 175 )
                  icxsni = mstz( 75 )
                  if( icxsni .gt. 2 .or. icxsni .lt. 0 ) icxsni = 0

                  icxspi = mstz( 117 )

                  irndmode = mstz( 78 )
                  idmprijk = mstz( 166 )
!$                if (idmprijk.eq.1) then
!$                 write(ErrCha,'("*** Warning ! idmprijk=1 is not",
!$   &                  " compatible with OpenMP enabled executable.",
!$   &                  " Ignored.")')
!$                 ErrID=''
!$                 call ErrWrite(ErrID,ErrCha)
!$                endif

                  irqmd = mstz( 79 )

                  ifiss = mstz( 143 )
               if( ifiss .lt. 0 .or. ifiss .ge. 3 ) then
                  write(io,'(/" Error : ifission should be",
     &            " 0,  1, or 2.")')
                  ErrCha = ''
                  MsgID = 'L:3628/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/" Error : ifission should be",
     &            " 0,  1, or 2.")')
                  goto 999
               elseif( ifiss .eq. 0 ) then ! default treatment
                  ifiss = 1
               end if

                  iqmdscm = mstz( 146 )
                  qmdscm_h0   = parz( 198 )
                  qmdscm_d    = parz( 199 )
                  qmdscm_rcls = parz( 200 )

*-----------------------------------------------------------------------

                  timeout = parz(190)

*-----------------------------------------------------------------------

                  iadjnt = mstz(114)
                  adjemax = parz(189)

      if(iadjnt.eq.1) then
       write(*,'("Adjoint mode for photon:",
     & " please see phits/utility/adjoint/adjoint for instruction")')
      elseif(iadjnt.eq.2) then
       write(*,'("Adjoint mode for charged particles:",
     & " please see phits/utility/adjoint/charged for instruction")')
       rcmin(:)=1.0e9 ! all nuclear reactions are prohibited
      endif

*-----------------------------------------------------------------------

                  inatur  = mstz(118)

*-----------------------------------------------------------------------

                  iwwbias  = mstz(125)

*-----------------------------------------------------------------------

                  ist_cut  = mstz(134)
                  ist_bat  = mstz(136)

*-----------------------------------------------------------------------

                  itallech = mstz(135)

*-----------------------------------------------------------------------
                  igsline = mstz(138)

*-----------------------------------------------------------------------
                  isans = mstz(156)  ! T.Sato 2022/03/02 change from 149 to 156

*-----------------------------------------------------------------------
                  mdbpseud = mstz(157)

*-----------------------------------------------------------------------
c Ogawa 2023/06/02  KURBUC exclusive mode
                  itsxcl = mstz(160)
                  rathe1  = parz(205)
                  rathe2  = parz(206)

*-----------------------------------------------------------------------
*           photonuclear reaction GDR part
*-----------------------------------------------------------------------

                  ipnint = mstz( 69 )


*-----------------------------------------------------------------------
*           multiplying factor for photonuclear interaction CS
*-----------------------------------------------------------------------

                  pnimul = parz( 186 )

*-----------------------------------------------------------------------
*           photon-induced muon pair production
*-----------------------------------------------------------------------

                  igmuppd = mstz( 137 )

*-----------------------------------------------------------------------
*           multiplying factor for photon-induced muon pair production CS
*-----------------------------------------------------------------------

                  gmumul = parz( 194 )

*-----------------------------------------------------------------------
*           muon interaction
*-----------------------------------------------------------------------
                  imuint = mstz( 84 )
                  emumin = parz( 180 )
                  emumax = parz( 181 )
                  prmui1 = parz( 188 )   ! S.Abe 2016/03/09
                  imucap = mstz( 104 )
                  imubrm = mstz( 109 )
                  imuppd = mstz( 110 )

*-----------------------------------------------------------------------
*           TETRA parameter
*-----------------------------------------------------------------------
                  itgchk = mstz(158)
                  itetvol = mstz(105)
                  if(abs(itetra).ne.1)then
                   ntetsurf = mstz(121)
                   ntetelem = mstz(122)
                  else
                   if(mstz(121).ne.ntetsurf)then
                    write(*,'("*** TETRA WARNING: ",
     &                   "ntetsurf is adjusted from ",i8,
     &                   " to ",i8,
     &                   ", which is written in Tetra.bin")')
     &                   mstz(121),ntetsurf
                    mstz(121)=ntetsurf
                   endif
                   if(mstz(122).ne.ntetelem)then
                    write(*,'("*** TETRA WARNING: ",
     &                   "ntetelem is adjusted from ",i8,
     &                   " to ",i8,
     &                   ", which is written in Tetra.bin")')
     &                   mstz(122),ntetelem
                    mstz(122)=ntetelem
                   endif
                  endif

*-----------------------------------------------------------------------

                  if( mstz( 65 ) .eq. 0 ) then

                     kfdcy(19) = 0

                  end if

*-----------------------------------------------------------------------
*              dump all
*-----------------------------------------------------------------------

                  idumpall = mstz( 60 )

                  if(idumpall.ne.0.and.nt4domp.gt.0)then
                   write(ErrCha,'('' **** Error ''
     &                  '' T-4Dtrack is not compatible''
     &                  '' with idumpall option for OpenMP'')')
                   ErrID = 'L:3773/R:setpar/F:read00.f' !W02_001_001
                   call ErrWrite(ErrID,ErrCha)
                   goto 999
                  endif

*-----------------------------------------------------------------------
*              e-mode; igamma = 2, nevap = 3, ! 2013/04/17 Ogawa
*                      emcnf = 10000, emcpf = 10000   :2009/06/02 KN
*-----------------------------------------------------------------------

                  iemode = mstz( 59 )

                  if( iemode .ne. 0 .and. igamma .eq. 0) then

                   write(ErrCha,'('' **** Warning ''
     &                 '' igamma = 0 is not allowed when ''
     &                 '' iemode != 0. igamma has been changed to 2'')')
                   ErrID = 'L:3790/R:setpar/F:read00.f' !W02_001_001
                   call ErrWrite(ErrID,ErrCha)

                     igamma = 2        ! 2013/04/17 Ogawa
                     mstz( 12 ) = 2    ! 2013/04/17 Ogawa

                     jevap  = 3
                     mstz( 6 ) = 3

                     if( emin(2) .lt. 1.d-11 ) then ! 2025/01/14 Ogawa Default neutron cutoff was 10^-10 eV even though JENDL minimum is 10^-11 eV.
                        emin(2) = 1.d-11
                        parz(2) = emin(2)
                     end if

                     parz( 151 ) = 10000.0
                     parz( 157 ) = 10000.0

                  end if

                  ge1 = parz( 172 )
                  ge2 = parz( 173 )

                  emodem = parz( 193 )

*-----------------------------------------------------------------------
*  emin & dmax adjustment (iegs & nucdata) T.Sato 2017/06/09, S.H. 2020.3.11
*-----------------------------------------------------------------------
                    dnmax( 2 )  = parz( 132 )
!       if(mstz(128).eq.1) then ! nucdata=1 for JENDL
       if(mstz(128).eq.1.or.(icntl.ge.7.and.icntl.le.11).or.
     &   icntl.eq.13.or.icntl.eq.17) then ! nucdata=1, icntl = 7-11, 13, 17
       else ! nucdata=0, S.H. added 2020.3.11
                 km2 = 0
                 kd2 = 0
                 kcm2 = 0
                 do kk = 1, ipara
                     if( ipsq(kk) .eq. -2   ) km2 = 1 ! emin(2)
                     if( ipsq(kk) .eq. -132 ) kd2 = 1 ! dmax(2)
                     if( ipsq(kk) .eq. -112 ) kcm2 = 1 ! cmin(2)
                 end do
                 if( km2 .eq. 0 ) then
                    parz( 2 ) = 1.0d-3
                    emin( 2 ) = 1.0d-3
                    if ( kcm2 .eq. 0 ) rcmin(2 ) = 1.0d-3
                 end if
                 if( kd2 .eq. 0 ) then
                    parz( 132 ) = parz( 2 )
                    dnmax( 2 )  = parz( 2 )
                 end if
       endif

                    dnmax( 14 ) = parz( 144 )

       if(mstz(85).eq.-1) then ! iegs=-1 for original photon
!       else if(mstz(85).eq.0) then ! negs=0, S.H. added 2020.3.11
       else if(mstz(85).eq.0.or.(icntl.ge.7.and.icntl.le.11).or.
     & icntl.eq.13.or.icntl.eq.17) then  ! for icntl = 7-11, 13, 17, negs is set to 0 to reduce the computational time
                 km14 = 0
                 kd14 = 0
                 kcm14 = 0
                 do kk = 1, ipara
                     if( ipsq(kk) .eq. -14  ) km14 = 1
                     if( ipsq(kk) .eq. -144 ) kd14 = 1
                     if( ipsq(kk) .eq. -124 ) kcm14 = 1
                 end do
                 if( km14 .eq. 0 ) then
                    parz( 14 ) = 1d9
                    emin( 14 ) = 1d9
                    if ( kcm14 .eq. 0 ) rcmin(14 ) = 1d9
                 end if
                 if( kd14 .eq. 0 ) then
                    parz( 144 ) = 1d9
                    dnmax( 14 ) = 1d9
                 end if
       endif

       if(icntl.eq.5.or.icntl.eq.6) then ! for icntl=5,6, negs should be -1, T.Sato 2018/10/02
        mstz(85)=-1
        parz( 12 ) = 0.001d0
        emin( 12 ) = 0.001d0
        parz( 13 ) = 0.001d0
        emin( 13 ) = 0.001d0
        parz( 14 ) = 0.001d0
        emin( 14 ) = 0.001d0
        parz( 142 ) = 10000.d0
        dnmax( 12 ) = 10000.d0
        parz( 143 ) = 10000.d0
        dnmax( 13 ) = 10000.d0
        parz( 144 ) = 10000.d0
        dnmax( 14 ) = 10000.d0
       end if

       if(mstz(85).le.0.and.emin(12).le.1000.0d0) then
          mstz(37)=0 ! T.Sato, 2017/08/09 photon produce electron
          parz(201) = 3.0d0 ! hirata 2023/05/26 extending 'xsmemory' for reading electron XS.
       endif

       if( icntl.eq.1 .and. inucr.eq.100 ) then ! frtati 2022/12/28
         if( dnmax(ilibpart).eq.0.d0 .or.
     &       ilibpart.eq.2.and.dnmax(2).eq.20.d0 ) then
           dnmax(ilibpart) = elibreadmax
         end if
       end if

*-----------------------------------------------------------------------
*  warning for dmax(2) setting in event generator Y. Iwamoto 2019/04/02
*-----------------------------------------------------------------------
       if( (iemode.eq.1 .or. iemode.eq.2) .and. parz(132).gt.20.0) then
	            ErrCha = ''
                    ErrID = 'L:3899/R:setpar/F:read00.f' !W02_001_001
                    call ErrWrite(ErrID,ErrCha)
                   write(*,'("e-mode dose not work well ",
     &              "when dmax(2) is above 20 MeV.")')
       end if

*-----------------------------------------------------------------------
*           EGS
*-----------------------------------------------------------------------

                  iegsemi = max(0,mstz( 85 )) ! T.Sato 2017/06/09 mstz(85)=-1 for original photon
                  iegsout = mstz( 102 )
! T.Sato 2015/07/25 ! PEGS run option
                  ipegs= mstz( 106 )

! T.Sato 2016/02/29, do not use EGS5 for icntl <> 0
                  if(icntl .ne. 0) then
                   iegsemi = 0
                   iegsout = 0
                   ipegs = 0
                  endif

! T.Sato 2016/03/02, for ipegs >= 1, iegsout should not be 0
                  if(ipegs .ge. 1 .and. iegsout.eq.0) iegsout = 1

! T.Sato 2014/08/29, Check input parameter for EGS
                  if(mstz(86).le.-1.or.mstz(86).ge.2) then
                   write(io,'(" Error in EGS parameter!! ",
     &             "iedgfl (=",i5,") should be 0 or 1")') mstz(86)
                   ErrCha = ''
                   MsgID = 'L:3929/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'(" Error in EGS parameter!! ",
     &             "iedgfl (=",i5,") should be 0 or 1")') mstz(86)
                   goto 999
                  endif

                  if(mstz(87).le.-1.or.mstz(87).ge.2) then
                   write(io,'(" Error in EGS parameter!! ",
     &             "iauger (=",i5,") should be 0 or 1")') mstz(87)
                   ErrCha = ''
                   MsgID = 'L:3940/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'(" Error in EGS parameter!! ",
     &             "iauger (=",i5,") should be 0 or 1")') mstz(87)
                   goto 999
                  endif

                  if(mstz(88).le.-1.or.mstz(88).ge.2) then
                   write(io,'(" Error in EGS parameter!! ",
     &             "iraylr (=",i5,") should be 0 or 1")') mstz(88)
                   ErrCha = ''
                   MsgID = 'L:3951/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'(" Error in EGS parameter!! ",
     &             "iraylr (=",i5,") should be 0 or 1")') mstz(88)
                   goto 999
                  endif

                  if(mstz(89).le.-1.or.mstz(89).ge.1) then
                   write(io,'(" Error in EGS parameter!! ",
     &             "lpolar (=",i5,") should be 0 at this moment")')
     &             mstz(89)
                   ErrCha = ''
                   MsgID = 'L:3963/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'(" Error in EGS parameter!! ",
     &             "lpolar (=",i5,") should be 0 at this moment")')
     &             mstz(89)
                   goto 999
                  endif

                  if(mstz(90).le.-1.or.mstz(90).ge.2) then
                   write(io,'(" Error in EGS parameter!! ",
     &             "incohr (=",i5,") should be 0 or 1")') mstz(90)
                   ErrCha = ''
                   MsgID = 'L:3975/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'(" Error in EGS parameter!! ",
     &             "incohr (=",i5,") should be 0 or 1")') mstz(90)
                   goto 999
                  endif

                  if(mstz(91).le.-1.or.mstz(91).ge.2) then
                   write(io,'(" Error in EGS parameter!! ",
     &             "iprofr (=",i5,") should be 0 or 1")') mstz(91)
                   ErrCha = ''
                   MsgID = 'L:3986/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'(" Error in EGS parameter!! ",
     &             "iprofr (=",i5,") should be 0 or 1")') mstz(91)
                   goto 999
                  endif

                  if(mstz(92).le.-1.or.mstz(92).ge.2) then
                   write(io,'(" Error in EGS parameter!! ",
     &             "impacr (=",i5,") should be 0 or 1")') mstz(92)
                   ErrCha = ''
                   MsgID = 'L:3997/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'(" Error in EGS parameter!! ",
     &             "impacr (=",i5,") should be 0 or 1")') mstz(92)
                   goto 999
                  endif

                  if(mstz(95).le.-1.or.mstz(95).ge.2) then
                   write(io,'(" Error in EGS parameter!! ",
     &             "ieispl (=",i5,") should be 0 or 1")') mstz(95)
                   ErrCha = ''
                   MsgID = 'L:4008/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'(" Error in EGS parameter!! ",
     &             "ieispl (=",i5,") should be 0 or 1")') mstz(95)
                   goto 999
                  endif

                  if(mstz(97).le.-1.or.mstz(97).ge.2) then
                   write(io,'(" Error in EGS parameter!! ",
     &             "ibrdst (=",i5,") should be 0 or 1")') mstz(97)
                   ErrCha = ''
                   MsgID = 'L:4019/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'(" Error in EGS parameter!! ",
     &             "ibrdst (=",i5,") should be 0 or 1")') mstz(97)
                   goto 999
                  endif

                  if(mstz(98).le.-1.or.mstz(98).ge.2) then
                   write(io,'(" Error in EGS parameter!! ",
     &             "iprdst (=",i5,") should be 0 or 1")') mstz(98)
                   ErrCha = ''
                   MsgID = 'L:4030/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'(" Error in EGS parameter!! ",
     &             "iprdst (=",i5,") should be 0 or 1")') mstz(98)
                   goto 999
                  endif

                  if(mstz(99).le.-1.or.mstz(99).ge.2) then
                   write(io,'(" Error in EGS parameter!! ",
     &             "iphter (=",i5,") should be 0 or 1")') mstz(99)
                   ErrCha = ''
                   MsgID = 'L:4041/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'(" Error in EGS parameter!! ",
     &             "iphter (=",i5,") should be 0 or 1")') mstz(99)
                   goto 999
                  endif

                  if(mstz(100).le.-1.or.mstz(100).ge.2) then
                   write(io,'(" Error in EGS parameter!! ",
     &             "ibound (=",i5,") should be 0 or 1")') mstz(100)
                   ErrCha = ''
                   MsgID = 'L:4052/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'(" Error in EGS parameter!! ",
     &             "ibound (=",i5,") should be 0 or 1")') mstz(100)
                   goto 999
                  endif

                  if(mstz(101).le.-1.or.mstz(101).ge.3) then
                   write(io,'(" Error in EGS parameter!! ",
     &             "iaprim (=",i5,") should be 0 - 2")') mstz(101)
                   ErrCha = ''
                   MsgID = 'L:4063/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'(" Error in EGS parameter!! ",
     &             "iaprim (=",i5,") should be 0 - 2")') mstz(101)
                   goto 999
                  endif
! T.Sato 2014/08/29, End of check parameter

               if( iegsemi .eq. 0 ) then

                     kfdcy(20) = 0

               else

                     mstz( 27 ) = 1
                     matadd = mstz( 27 )

                     km12 = 0
                     km13 = 0
                     km14 = 0
                     kd12 = 0
                     kd13 = 0
                     kd14 = 0
                     kcm12 = 0
                     kcm13 = 0
                     kcm14 = 0

                 do kk = 1, ipara
                     if( ipsq(kk) .eq. -12  ) km12 = 1
                     if( ipsq(kk) .eq. -13  ) km13 = 1
                     if( ipsq(kk) .eq. -14  ) km14 = 1
                     if( ipsq(kk) .eq. -142 ) kd12 = 1
                     if( ipsq(kk) .eq. -143 ) kd13 = 1
                     if( ipsq(kk) .eq. -144 ) kd14 = 1
                     if( ipsq(kk) .eq. -122  ) kcm12 = 1
                     if( ipsq(kk) .eq. -123  ) kcm13 = 1
                     if( ipsq(kk) .eq. -124  ) kcm14 = 1
                 end do

                 if( km12 .eq. 0 ) then
                    parz( 12 ) = 0.1d0
                    emin( 12 ) = 0.1d0
                    if ( kcm12 .eq. 0 ) rcmin(12 ) = 0.1d0
                 end if
                 if( km13 .eq. 0 ) then
                    parz( 13 ) = 0.1d0
                    emin( 13 ) = 0.1d0
                    if ( kcm13 .eq. 0 ) rcmin(13 ) = 0.1d0
                 end if
                 if( km14 .eq. 0 ) then
                    parz( 14 ) = 0.001d0
                    emin( 14 ) = 0.001d0
                    if ( kcm14 .eq. 0 ) rcmin(14 ) = 0.001d0
                 end if
                 if( kd12 .eq. 0 ) then
                    if( mstz(85) .eq. 2 ) then   ! y.sakaki 2022/12 negs=2. egs high energy mode (10TeV)
                       parz( 142 ) = 10000000.d0
                       dnmax( 12 ) = 10000000.d0
                    else
                       parz( 142 ) = 1000.d0
                       dnmax( 12 ) = 1000.d0
                    end if
                 end if
                 if( kd13 .eq. 0 ) then
                    if( mstz(85) .eq. 2 ) then   ! y.sakaki 2022/12 negs=2. egs high energy mode (10TeV)
                       parz( 143 ) = 10000000.d0
                       dnmax( 13 ) = 10000000.d0
                    else
                       parz( 143 ) = 1000.d0
                       dnmax( 13 ) = 1000.d0
                    end if
                 end if
                 if( kd14 .eq. 0 ) then
                    if( mstz(85) .eq. 2 ) then   ! y.sakaki 2022/12 negs=2. egs high energy mode (10TeV)
                       parz( 144 ) = 10000000.d0
                       dnmax( 14 ) = 10000000.d0
                    else
                       parz( 144 ) = 1000.d0
                       dnmax( 14 ) = 1000.d0
                    end if
                 end if

               end if

! Check emin(12)-(14) consistency ! T.Sato 2017/06/25
       if(emin(12).gt.emin(13)) then
        ErrCha = ''
        MsgID = 'L:4150/R:setpar/F:read00.f'
        call ErrWrite(MsgID, ErrCha)
        write(jo,*) 'emin(12) should be equal to emin(13)'
        write(jo,'("emin(12) is adjusted to emin(13) =",es12.3)')
     &  emin(13)
        emin(12)=emin(13)
        if ( kcm12 .eq. 0 ) rcmin(12) = emin(12) ! T.Sato 2019/01/16
       endif
       if(emin(13).gt.emin(12)) then
        ErrCha = ''
        MsgID = 'L:4160/R:setpar/F:read00.f'
        call ErrWrite(MsgID, ErrCha)
        write(jo,*) 'emin(13) should be equal to emin(12)'
        write(jo,'("emin(13) is adjusted to emin(12) =",es12.3)')
     &  emin(12)
        emin(13)=emin(12)
        if ( kcm13 .eq. 0 ) rcmin(13) = emin(13) ! T.Sato 2019/01/16
       endif
       if(emin(14).gt.emin(12)) then
        ErrCha = ''
        MsgID = 'L:4170/R:setpar/F:read00.f'
        call ErrWrite(MsgID, ErrCha)
        write(jo,*) 'emin(14) should be less than emin(12)'
        write(jo,'("emin(14) is adjusted to emin(12) =",es12.3)')
     &  emin(12)
        emin(14)=emin(12)
        if ( kcm14 .eq. 0 ) rcmin(14) = emin(14) ! T.Sato 2019/01/16
       endif
! T.Sato 2023/03/30 warning for too high emin(12) for EGS mode
        if(mstz(85).ge.1.and.emin(12).gt.10.0.and.icntl.eq.0) then ! add condition of icntl
         ErrCha = ''
         MsgID = 'L:4181/R:setpar/F:read00.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,*)'emin(12-13) < 10 is recommended when negs >= 1'
        endif

! Ogawa. 2024/6/3
! These warnings are intended to avoid distortion of transport phenomena.
! Neutral particles, which can travel further even at low E, should be transported until threshold.
! High cut-off saves little CPU time but has various disadvantages such as cutoff photons are not scored by T-deposit, capture gammas are disregarded, etc.
      if(icntl.eq.0.and.mstz(85) .ne. 0 .and.
     & abs(emin(14) - 1.d-3) .gt. 1.d-5 ) then
       write(ErrCha,'("*** Warning! emin(14) is changed. Transport ",
     & "physics is distorted. Default emin(14) is encouraged")')
       ErrID = 'L:4194/R:setpar/F:read00.f'
       call ErrWrite(ErrID,ErrCha)
      endif
      if (icntl.eq.0.and.mstz(128).eq.1.and.
     & abs(emin(2) - 1.d-11).gt.1.d-13) then
       write(ErrCha,'("*** Warning! emin(2) is changed. Transport ",
     & "physics is distorted. Default emin(2) is encouraged")')
       ErrID = 'L:4201/R:setpar/F:read00.f'
       call ErrWrite(ErrID,ErrCha)
      endif

! Check ikerman, ikermap, T.Sato 2018/03/06
      if(mstz(129).eq.1) then
       ikerman=1 ! do not use kerma for neutron
      elseif(mstz(129).eq.2) then
       ikerman=2 ! use kerma for neutron
      else
       if(iemode.eq.0) then ! non-event generator
        ikerman=2  ! use kerma for neutron for non-EG mode
       else
        ikerman=1  ! do not use kerma for EG mode
       endif
      endif

      if(mstz(130).eq.1) then
       ikermap=1 ! do not use kerma for photon
      elseif(mstz(130).eq.2) then
       ikermap=2 ! use kerma for photon
      else
       if(mstz(85).ge.1) then ! EGS mode y.sakaki 2022/12
        ikermap=1  ! do not use kerma for photon in EGS mode
       elseif(emin(12).lt.10.0) then ! transport electron in non-EGS mode
        ikermap=1  ! do not use kerma for photon in original-PHITS mode
        if(icntl.eq.0) then
                       ErrCha = ''
                       MsgID = 'L:4229/R:setpar/F:read00.f'
                       call ErrWrite(MsgID, ErrCha)
                   write(jo,*) 'Warning: negs = 1 or 2 is recommended ',
     &  'for electron transport'
        end if
       else
        ikermap=2  ! use kerma for photon
       endif
      endif


*-----------------------------------------------------------------------
*           gravity
*-----------------------------------------------------------------------

                  grav(1) = parz( 165 )
                  grav(2) = parz( 166 )
                  grav(3) = parz( 167 )

                  igrav = 0

                  if( grav(1) .ne. 0.d0 .or.
     &                grav(2) .ne. 0.d0 .or.
     &                grav(3) .ne. 0.d0 ) then

                     agrav = sqrt( grav(1)**2 + grav(2)**2
     &                           + grav(3)**2 )

                     grav(1) = grav(1) / agrav
                     grav(2) = grav(2) / agrav
                     grav(3) = grav(3) / agrav

                     igrav = 1

                  end if

*-----------------------------------------------------------------------
*           nspred = 3 is not complete
*-----------------------------------------------------------------------

            if( nspred .eq. 3 ) then

                  write(io,'(/" Error : nspred = 3 is not",
     &            " complete now, Sorry !!")')

                  ErrCha = ''
                  MsgID = 'L:4275/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/" Error : nspred = 3 is not",
     &            " complete now, Sorry !!")')

                  goto 999

            end if

*-----------------------------------------------------------------------
*           delt0 = deltc / 10.0 always
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*           check cross section directory file ( unit = 7 )
*-----------------------------------------------------------------------

            if( icfn(7) .ne. 0 ) then

               inquire( file = chfn(7), exist = exex )

               if( exex .eqv. .false. ) then

                  write(io,'(/" Error : input data file for",
     &            " cross section directory does not exist."/
     &            " file name = ",200a1)')
     &            ( chfn(7)(i:i),i=1, ilfn(7) )

                  ErrCha = ''
                  MsgID = 'L:4304/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/" Error : input data file for",
     &            " cross section directory does not exist."/
     &            " file name = ",200a1)')
     &            ( chfn(7)(i:i),i=1, ilfn(7) )

                  goto 999

               end if

            end if

*-----------------------------------------------------------------------
*           ivoxel = 1
*-----------------------------------------------------------------------

            if( ivoxel .eq. 1 ) then

               inquire( file = chfn(18), exist = exex )

               if( exex .eqv. .false. ) then

                  write(io,'(/" Error : voxel data file",
     &            " does not exist."/
     &            " file name = ",200a1)')
     &            ( chfn(18)(i:i),i=1, ilfn(18) )

                  ErrCha = ''
                  MsgID = 'L:4333/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/" Error : voxel data file",
     &            " does not exist."/
     &            " file name = ",200a1)')
     &            ( chfn(18)(i:i),i=1, ilfn(18) )

                  goto 999

               end if

            end if

*-----------------------------------------------------------------------
*           icells = 1
*-----------------------------------------------------------------------

            if( icells .eq. 1 ) then

               inquire( file = chfn(19), exist = exex )

               if( exex .eqv. .false. ) then

                  write(io,'(/" Error : gcell data file",
     &            " does not exist."/
     &            " file name = ",200a1)')
     &            ( chfn(19)(i:i),i=1, ilfn(19) )

                  ErrCha = ''
                  MsgID = 'L:4362/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/" Error : gcell data file",
     &            " does not exist."/
     &            " file name = ",200a1)')
     &            ( chfn(19)(i:i),i=1, ilfn(19) )

                  goto 999

               end if

            end if

*-----------------------------------------------------------------------
*        nuclear reaction case
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*           cgview files
*-----------------------------------------------------------------------

            if( icntl .eq. 2 ) then

               if( icfn(2) .eq. 0 ) icfn(2) = 2
               if( icfn(3) .eq. 0 ) icfn(3) = 2

            end if

*-----------------------------------------------------------------------
*           marspf files
*-----------------------------------------------------------------------

            if( icntl .eq. 4 ) then

               if( icfn(4) .eq. 0 ) icfn(4) = 2

            end if

*-----------------------------------------------------------------------

            if( icntl .eq. 2 .or. icntl .eq. 4 ) then

               if( iregn .eq. 0 ) then

                  write(io,'(" **** Error : in CG definition")')
                  write(io,'(" You should choose CG input.")')

                  ErrCha = ''
                  MsgID = 'L:4411/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(" **** Error : in CG definition")')
                  write(jo,'(" You should choose CG input.")')

                  goto 999

               end if

            end if

*-----------------------------------------------------------------------
*     choice of CG (region) or GG (cell)
*-----------------------------------------------------------------------

            icgg = 0


            if( ( iregn .ne. 0 .and. ibody .eq. 0 ) .or.
     &          ( iregn .eq. 0 .and. ibody .ne. 0 ) ) then

               write(io,'(" **** Error : in CG definition")')
               write(io,'(" [region] or [body] is missing.")')

               ErrCha = ''
               MsgID = 'L:4436/R:setpar/F:read00.f'
               call ErrWrite(MsgID, ErrCha)
               write(jo,'(" **** Error : in CG definition")')
               write(jo,'(" [region] or [body] is missing.")')

               goto 999

            end if

            if( ( igcel .ne. 0 .and. igsuf .eq. 0 ) .or.
     &          ( igcel .eq. 0 .and. igsuf .ne. 0 ) ) then

               write(io,'(" **** Error : in GG definition")')
               write(io,'(" [cell] or [surface] is missing.")')

               ErrCha = ''
               MsgID = 'L:4452/R:setpar/F:read00.f'
               call ErrWrite(MsgID, ErrCha)
               write(jo,'(" **** Error : in GG definition")')
               write(jo,'(" [cell] or [surface] is missing.")')

               goto 999

            end if

            if( iregn .eq. 0 .and. igcel .eq. 0 ) then

               write(io,'(" **** Error : in CG, GG definition")')
               write(io,'(" both [region] and [cell] are missing.")')

               ErrCha = ''
               MsgID = 'L:4467/R:setpar/F:read00.f'
               call ErrWrite(MsgID, ErrCha)
               write(jo,'(" **** Error : in CG, GG definition")')
               write(jo,'(" both [region] and [cell] are missing.")')

               goto 999

            end if


            if( iregn .ne. 0 .and. abs(idumpall) .eq. 1 ) then

               write(io,'(" **** Error : dumpall=1 but CG")')
               write(io,'(" dumpall is only available in GG.")')

               ErrCha = ''
               MsgID = 'L:4483/R:setpar/F:read00.f'
               call ErrWrite(MsgID, ErrCha)
               write(jo,'(" **** Error : dumpall=1 but CG")')
               write(jo,'(" dumpall is only available in GG.")')

               goto 999

            end if

*-----------------------------------------------------------------------
*        GG case : iregn <-- igcel
*-----------------------------------------------------------------------

         if( igcel .ne. 0 ) then

               icgg = 1
               iregn = igcel

         else

               igcel = iregn

         end if

*-----------------------------------------------------------------------
*           rsous : number of total history
*-----------------------------------------------------------------------

            rsous = dble( maxcas ) * dble( maxbch )

            if( dble( irskip ) .ge. rsous ) then






          ErrCha = ''
          ErrID = 'L:4521/R:setpar/F:read00.f' !E02_002_001
          call ErrWriteIO(ErrID,ErrCha,io)
          call ErrWriteIO(ErrID,ErrCha,jo)

          write(io,*) '**** Error : irskip >= maxcas * maxbch'
          write(io,*) 'maxcas is intger overflow'
     &    ,'limit 2,147,483,647'
          write(io,*) ' maxcas =', maxcas
          write(io,*) ' maxbch =', maxbch
          write(io,*) ' total  =', rsous
          write(io,*) ' irskip =', irskip

          ErrCha = ''
          MsgID = 'L:4534/R:setpar/F:read00.f'
          call ErrWrite(MsgID, ErrCha)
          write(jo,*) '**** Error : irskip >= maxcas * maxbch'
          write(jo,*) 'maxcas is intger overflow'
     &    ,'limit 2,147,483,647'
          write(jo,*) ' maxcas =', maxcas
          write(jo,*) ' maxbch =', maxbch
          write(jo,*) ' total  =', rsous
          write(jo,*) ' irskip =', irskip

               goto 999

            end if

*-----------------------------------------------------------------------
*           check of input data file for photon emission
*-----------------------------------------------------------------------









*-----------------------------------------------------------------------
*     check of input dump file for dumpall calculation
*-----------------------------------------------------------------------

      if( icntl .eq. 12 ) then

*-----------------------------------------------------------------------

         if( npe .le. 1 ) then

               inquire( file = chfn(15), exist = exex )

               if( exex .eqv. .false. ) then

                  write(io,'(/" Error : input dump file for",
     &            " dumpall calculation does not exist."/
     &            " file name = ",200a1)')
     &            ( chfn(15)(i:i),i=1, ilfn(15) )

                  ErrCha = ''
                  MsgID = 'L:4580/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/" Error : input dump file for",
     &            " dumpall calculation does not exist."/
     &            " file name = ",200a1)')
     &            ( chfn(15)(i:i),i=1, ilfn(15) )

                  goto 999

               end if

               if( idumpall .eq. 1 ) then

                  open(9, file = chfn(15),
     &                    form='unformatted',status = 'unknown' )

                  read(9,iostat=ios,err=999) ncol

               else if( idumpall .eq. -1 ) then

                  open(9, file = chfn(15),
     &                    form='formatted',status = 'unknown' )

                  read(9,*,iostat=ios,err=999) ncol

               end if

                  if( ios .eq. -1 ) then

                     write(io,'(/" Error : input dump file for",
     &               " dumpall is empty."/
     &               " file name = ",200a1)')
     &               ( chfn(15)(i:i),i=1, ilfn(15) )

                     ErrCha = ''
                     MsgID = 'L:4615/R:setpar/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'(/" Error : input dump file for",
     &               " dumpall is empty."/
     &               " file name = ",200a1)')
     &               ( chfn(15)(i:i),i=1, ilfn(15) )

                     goto 999

                  end if

                     rewind 9

*-----------------------------------------------------------------------
*        mpi case
*-----------------------------------------------------------------------

         else if( me .gt. 0 ) then

            if( mstz(61) .eq. 0 ) then

                  filnm = filhd(1:nfihd) // chfn(15)(1:ilfn(15))

                  inquire( file = filnm, exist = exex )

                  if( exex .eqv. .false. ) then

                     write(io,'(/" Error : input dump file for",
     &               " dumpall calculation does not exist."/
     &               " file name = ",200a1)')
     &               ( filnm(i:i),i=1,nfihd+ilfn(15) )

                     ErrCha = ''
                     MsgID = 'L:4648/R:setpar/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'(/" Error : input dump file for",
     &               " dumpall calculation does not exist."/
     &               " file name = ",200a1)')
     &               ( filnm(i:i),i=1,nfihd+ilfn(15) )

                     goto 999

                  end if

               if( idumpall .eq. 1 ) then

                  open(9, file = filnm(1:nfihd+ilfn(15)),
     &                    form='unformatted',status = 'unknown' )

                  read(9,iostat=ios,err=999) ncol

               else if( idumpall .eq. -1 ) then

                  open(9, file = filnm(1:nfihd+ilfn(15)),
     &                    form='formatted',status = 'unknown' )

                  read(9,*,iostat=ios,err=999) ncol

               end if

                  if( ios .eq. -1 ) then

                     write(io,'(/" Error : input dump file for",
     &               " dumpall calculation is empty."/
     &               " file name = ",200a1)')
     &               ( filnm(i:i),i=1,nfihd+ilfn(15) )

                     ErrCha = ''
                     MsgID = 'L:4683/R:setpar/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'(/" Error : input dump file for",
     &               " dumpall calculation is empty."/
     &               " file name = ",200a1)')
     &               ( filnm(i:i),i=1,nfihd+ilfn(15) )

                     goto 999

                  end if

                     rewind 9

            else if( mstz(61) .eq. 1 ) then

                  iorder = aint(log10(real(npe))) + 1
                  if ( iorder .lt. 3) iorder = 3

                  write(chme,'(i5.5)') me

                  filnm = filhd(1:nfihd) // chfn(15)(1:ilfn(15))
     &                  //'.' // chme(6-iorder:5)

                  inquire( file = filnm, exist = exex )

                  if( exex .eqv. .false. ) then

                     write(io,'(/" Error : input dump file for",
     &               " dumpall calculation does not exist."/
     &               " file name = ",200a1)')
     &               ( filnm(i:i),i=1,nfihd+ilfn(15)+1+iorder )

                     ErrCha = ''
                     MsgID = 'L:4716/R:setpar/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'(/" Error : input dump file for",
     &               " dumpall calculation does not exist."/
     &               " file name = ",200a1)')
     &               ( filnm(i:i),i=1,nfihd+ilfn(15)+1+iorder )

                     goto 999

                  end if

               if( idumpall .eq. 1 ) then

                  open(9, file = filnm(1:nfihd+ilfn(15)+1+iorder),
     &                    form='unformatted',status = 'unknown' )

                  read(9,iostat=ios,err=999) ncol

               else if( idumpall .eq. -1 ) then

                  open(9, file = filnm(1:nfihd+ilfn(15)+1+iorder),
     &                    form='formatted',status = 'unknown' )

                  read(9,*,iostat=ios,err=999) ncol

               end if

                  if( ios .eq. -1 ) then

                     write(io,'(/" Error : input dump file for",
     &               " dumpall calculation is empty."/
     &               " file name = ",200a1)')
     &               ( filnm(i:i),i=1,nfihd+ilfn(15)+1+iorder )

                     ErrCha = ''
                     MsgID = 'L:4751/R:setpar/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'(/" Error : input dump file for",
     &               " dumpall calculation is empty."/
     &               " file name = ",200a1)')
     &               ( filnm(i:i),i=1,nfihd+ilfn(15)+1+iorder )

                     goto 999

                  end if

                     rewind 9

            else if( mstz(61) .eq. 2 ) then

                  filnm = chfn(15)(1:ilfn(15))

                  inquire( file = filnm, exist = exex )

                  if( exex .eqv. .false. ) then

                     write(io,'(/" Error : input dump file for",
     &               " dumpall calculation does not exist."/
     &               " file name = ",200a1)')
     &               ( filnm(i:i),i=1,ilfn(15) )

                     ErrCha = ''
                     MsgID = 'L:4778/R:setpar/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'(/" Error : input dump file for",
     &               " dumpall calculation does not exist."/
     &               " file name = ",200a1)')
     &               ( filnm(i:i),i=1,ilfn(15) )

                     goto 999

                  end if

               if( idumpall .eq. 1 ) then

                  open(9, file = filnm(1:ilfn(15)),
     &                    form='unformatted',status = 'unknown' )

                  read(9,iostat=ios,err=999) ncol

               else if( idumpall .eq. -1 ) then

                  open(9, file = filnm(1:ilfn(15)),
     &                    form='formatted',status = 'unknown' )

                  read(9,*,iostat=ios,err=999) ncol

               end if

                  if( ios .eq. -1 ) then

                     write(io,'(/" Error : input dump file for",
     &               " dumpall calculation is empty."/
     &               " file name = ",200a1)')
     &               ( filnm(i:i),i=1,ilfn(15) )

                     ErrCha = ''
                     MsgID = 'L:4813/R:setpar/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'(/" Error : input dump file for",
     &               " dumpall calculation is empty."/
     &               " file name = ",200a1)')
     &               ( filnm(i:i),i=1,ilfn(15) )

                     goto 999

                  end if

                     rewind 9

            else if( mstz(61) .eq. 3 ) then

                  iorder = aint(log10(real(npe))) + 1
                  if ( iorder .lt. 3) iorder = 3

                  write(chme,'(i5.5)') me

                  filnm = chfn(15)(1:ilfn(15))
     &                  //'.' // chme(6-iorder:5)

                  inquire( file = filnm, exist = exex )

                  if( exex .eqv. .false. ) then

                     write(io,'(/" Error : input dump file for",
     &               " dumpall calculation does not exist."/
     &               " file name = ",200a1)')
     &               ( filnm(i:i),i=1,ilfn(15)+1+iorder )

                     ErrCha = ''
                     MsgID = 'L:4846/R:setpar/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'(/" Error : input dump file for",
     &               " dumpall calculation does not exist."/
     &               " file name = ",200a1)')
     &               ( filnm(i:i),i=1,ilfn(15)+1+iorder )

                     goto 999

                  end if

               if( idumpall .eq. 1 ) then

                  open(9, file = filnm(1:ilfn(15)+1+iorder),
     &                    form='unformatted',status = 'unknown' )

                  read(9,iostat=ios,err=999) ncol

               else if( idumpall .eq. -1 ) then

                  open(9, file = filnm(1:ilfn(15)+1+iorder),
     &                    form='formatted',status = 'unknown' )

                  read(9,*,iostat=ios,err=999) ncol

               end if

                  if( ios .eq. -1 ) then

                     write(io,'(/" Error : input dump file for",
     &               " dumpall calculation is empty."/
     &               " file name = ",200a1)')
     &               ( filnm(i:i),i=1,ilfn(15)+1+iorder )

                     ErrCha = ''
                     MsgID = 'L:4881/R:setpar/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'(/" Error : input dump file for",
     &               " dumpall calculation is empty."/
     &               " file name = ",200a1)')
     &               ( filnm(i:i),i=1,ilfn(15)+1+iorder )

                     goto 999

                  end if

                     rewind 9

            end if

         end if

      end if

*-----------------------------------------------------------------------
*     open ncut and gcut files
*-----------------------------------------------------------------------

      if( npe .le. 1 ) then

            if( incut .ne. 0 ) then

C S.H. revised for Dump-Restart on 2014/5/7
               if ( irestart .eq. 1 ) then

                  open(12, file = chfn(12),
     &                 form='unformatted',status = 'unknown',
     &                 access = 'append' )

               else

                  open(12, file = chfn(12),
     &                 form='unformatted',status = 'unknown' )

               end if

            end if

            if( igcut .ne. 0 ) then

               if ( irestart .eq. 1 ) then

                  open(13, file = chfn(13),
     &                 form='unformatted',status = 'unknown',
     &                 access = 'append' )

               else

                  open(13, file = chfn(13),
     &                 form='unformatted',status = 'unknown' )

               end if

            end if

            if( ipcut .ne. 0 ) then

               if ( irestart .eq. 1 ) then

                  open(10, file = chfn(10),
     &                 form='unformatted',status = 'unknown',
     &                 access = 'append' )

               else

                  open(10, file = chfn(10),
     &                 form='unformatted',status = 'unknown' )

               end if

            end if

*-----------------------------------------------------------------------

         if( icntl .ne. 12 ) then

            if( idumpall .eq. 1 ) then

C S.H. revised for Dump-Restart on 2014/5/7
               if ( irestart .eq. 1 ) then

                  open(9, file = chfn(15),
     &                 form='unformatted',status = 'unknown',
     &                 access = 'append' )

               else

                  open(9, file = chfn(15),
     &                 form='unformatted',status = 'unknown' )

               end if

            end if

            if( idumpall .eq. -1 ) then

               if ( irestart .eq. 1 ) then

                  open(9, file = chfn(15),
     &                 form='formatted',status = 'unknown',
     &                 access = 'append' )

               else

                  open(9, file = chfn(15),
     &                 form='formatted',status = 'unknown' )

               end if

            end if

         end if

*-----------------------------------------------------------------------
*     mpi case
*-----------------------------------------------------------------------

      else if( me .gt. 0 ) then

         if( icntl .ne. 12 ) then

            if( mstz(61) .eq. 0 ) then

                  filnm = filhd(1:nfihd) // chfn(15)(1:ilfn(15))

               if( idumpall .eq. 1 ) then

C S.H. revised for Dump-Restart on 2014/5/7
                  if ( irestart .eq. 1 ) then

                     open(9, file = filnm(1:nfihd+ilfn(15)),
     &                    form='unformatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(9, file = filnm(1:nfihd+ilfn(15)),
     &                    form='unformatted',status = 'unknown' )

                  end if

               else if( idumpall .eq. -1 ) then

                  if ( irestart .eq. 1 ) then

                     open(9, file = filnm(1:nfihd+ilfn(15)),
     &                    form='formatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(9, file = filnm(1:nfihd+ilfn(15)),
     &                    form='formatted',status = 'unknown' )

                  end if

               end if

            else if( mstz(61) .eq. 1 ) then

                  iorder = aint(log10(real(npe))) + 1
                  if ( iorder .lt. 3) iorder = 3

                  write(chme,'(i5.5)') me

                  filnm = filhd(1:nfihd) // chfn(15)(1:ilfn(15))
     &                  //'.' // chme(6-iorder:5)

               if( idumpall .eq. 1 ) then

                  if ( irestart .eq. 1 ) then

                     open(9, file = filnm(1:nfihd+ilfn(15)+1+iorder),
     &                    form='unformatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(9, file = filnm(1:nfihd+ilfn(15)+1+iorder),
     &                    form='unformatted',status = 'unknown' )

                  end if

               else if( idumpall .eq. -1 ) then

                  if ( irestart .eq. 1 ) then

                     open(9, file = filnm(1:nfihd+ilfn(15)+1+iorder),
     &                    form='formatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(9, file = filnm(1:nfihd+ilfn(15)+1+iorder),
     &                    form='formatted',status = 'unknown' )

                  end if

               end if

            else if( mstz(61) .eq. 2 ) then

                  filnm = chfn(15)(1:ilfn(15))

               if( idumpall .eq. 1 ) then

                  if ( irestart .eq. 1 ) then

                     open(9, file = filnm(1:ilfn(15)),
     &                    form='unformatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(9, file = filnm(1:ilfn(15)),
     &                    form='unformatted',status = 'unknown' )

                  end if

               else if( idumpall .eq. -1 ) then

                  if ( irestart .eq. 1 ) then

                     open(9, file = filnm(1:ilfn(15)),
     &                    form='formatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(9, file = filnm(1:ilfn(15)),
     &                    form='formatted',status = 'unknown' )

                  end if

               end if

            else if( mstz(61) .eq. 3 ) then

                  iorder = aint(log10(real(npe))) + 1
                  if ( iorder .lt. 3) iorder = 3

                  write(chme,'(i5.5)') me

                  filnm = chfn(15)(1:ilfn(15))
     &                  //'.' // chme(6-iorder:5)

               if( idumpall .eq. 1 ) then

                  if ( irestart .eq. 1 ) then

                     open(9, file = filnm(1:ilfn(15)+1+iorder),
     &                    form='unformatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(9, file = filnm(1:ilfn(15)+1+iorder),
     &                    form='unformatted',status = 'unknown' )

                  end if

               else if( idumpall .eq. -1 ) then

                  if ( irestart .eq. 1 ) then

                     open(9, file = filnm(1:ilfn(15)+1+iorder),
     &                    form='formatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(9, file = filnm(1:ilfn(15)+1+iorder),
     &                    form='formatted',status = 'unknown' )

                  end if

               end if

            end if


         end if

*-----------------------------------------------------------------------

            if( incut .ne. 0 ) then

               if( mstz(28) .eq. 0 ) then

                  filnm = filhd(1:nfihd) // chfn(12)(1:ilfn(12))

                  if ( irestart .eq. 1 ) then

                     open(12, file = filnm(1:nfihd+ilfn(12)),
     &                    form='unformatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(12, file = filnm(1:nfihd+ilfn(12)),
     &                    form='unformatted',status = 'unknown' )

                  end if

               else if( mstz(28) .eq. 1 ) then

                  iorder = aint(log10(real(npe))) + 1
                  if ( iorder .lt. 3) iorder = 3

                  write(chme,'(i5.5)') me

                  filnm = filhd(1:nfihd) // chfn(12)(1:ilfn(12))
     &                  //'.' // chme(6-iorder:5)

                  if ( irestart .eq. 1 ) then

                     open(12, file = filnm(1:nfihd+ilfn(12)+1+iorder),
     &                    form='unformatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(12, file = filnm(1:nfihd+ilfn(12)+1+iorder),
     &                    form='unformatted',status = 'unknown' )

                  end if

               else if( mstz(28) .eq. 2 ) then

                  filnm = chfn(12)(1:ilfn(12))

                  if ( irestart .eq. 1 ) then

                     open(12, file = filnm(1:ilfn(12)),
     &                    form='unformatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(12, file = filnm(1:ilfn(12)),
     &                    form='unformatted',status = 'unknown' )

                  end if

               else if( mstz(28) .eq. 3 ) then

                  iorder = aint(log10(real(npe))) + 1
                  if ( iorder .lt. 3) iorder = 3

                  write(chme,'(i5.5)') me

                  filnm = chfn(12)(1:ilfn(12))
     &                  //'.' // chme(6-iorder:5)

                  if ( irestart .eq. 1 ) then

                     open(12, file = filnm(1:ilfn(12)+1+iorder),
     &                    form='unformatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(12, file = filnm(1:ilfn(12)+1+iorder),
     &                    form='unformatted',status = 'unknown' )

                  end if

               end if

            end if

*-----------------------------------------------------------------------

            if( igcut .ne. 0 ) then

               if( mstz(29) .eq. 0 ) then

                  filnm = filhd(1:nfihd) // chfn(13)(1:ilfn(13))

                  if ( irestart .eq. 1 ) then

                     open(13, file = filnm(1:nfihd+ilfn(13)),
     &                    form='unformatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(13, file = filnm(1:nfihd+ilfn(13)),
     &                    form='unformatted',status = 'unknown' )

                  end if

               else if( mstz(29) .eq. 1 ) then

                  iorder = aint(log10(real(npe))) + 1
                  if ( iorder .lt. 3) iorder = 3

                  write(chme,'(i5.5)') me

                  filnm = filhd(1:nfihd) // chfn(13)(1:ilfn(13))
     &                  //'.' // chme(6-iorder:5)

                  if ( irestart .eq. 1 ) then

                     open(13, file = filnm(1:nfihd+ilfn(13)+1+iorder),
     &                    form='unformatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(13, file = filnm(1:nfihd+ilfn(13)+1+iorder),
     &                    form='unformatted',status = 'unknown' )

                  end if

               else if( mstz(29) .eq. 2 ) then

                  filnm = chfn(13)(1:ilfn(13))

                  if ( irestart .eq. 1 ) then

                     open(13, file = filnm(1:ilfn(13)),
     &                    form='unformatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(13, file = filnm(1:ilfn(13)),
     &                    form='unformatted',status = 'unknown' )

                  end if

               else if( mstz(29) .eq. 3 ) then

                  iorder = aint(log10(real(npe))) + 1
                  if ( iorder .lt. 3) iorder = 3

                  write(chme,'(i5.5)') me

                  filnm = chfn(13)(1:ilfn(13))
     &                  //'.' // chme(6-iorder:5)

                  if ( irestart .eq. 1 ) then

                     open(13, file = filnm(1:ilfn(13)+1+iorder),
     &                    form='unformatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(13, file = filnm(1:ilfn(13)+1+iorder),
     &                    form='unformatted',status = 'unknown' )

                  end if

               end if

            end if

*-----------------------------------------------------------------------

            if( ipcut .ne. 0 ) then

               if( mstz(31) .eq. 0 ) then

                  filnm = filhd(1:nfihd) // chfn(10)(1:ilfn(10))

                  if ( irestart .eq. 1 ) then

                     open(10, file = filnm(1:nfihd+ilfn(10)),
     &                    form='unformatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(10, file = filnm(1:nfihd+ilfn(10)),
     &                    form='unformatted',status = 'unknown' )

                  end if

               else if( mstz(31) .eq. 1 ) then

                  iorder = aint(log10(real(npe))) + 1
                  if ( iorder .lt. 3) iorder = 3

                  write(chme,'(i5.5)') me

                  filnm = filhd(1:nfihd) // chfn(10)(1:ilfn(10))
     &                  //'.' // chme(6-iorder:5)

                  if ( irestart .eq. 1 ) then

                     open(10, file = filnm(1:nfihd+ilfn(10)+1+iorder),
     &                    form='unformatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(10, file = filnm(1:nfihd+ilfn(10)+1+iorder),
     &                    form='unformatted',status = 'unknown' )

                  end if

               else if( mstz(31) .eq. 2 ) then

                  filnm = chfn(10)(1:ilfn(10))

                  if ( irestart .eq. 1 ) then

                     open(10, file = filnm(1:ilfn(10)),
     &                    form='unformatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(10, file = filnm(1:ilfn(10)),
     &                    form='unformatted',status = 'unknown' )

                  end if

               else if( mstz(31) .eq. 3 ) then

                  iorder = aint(log10(real(npe))) + 1
                  if ( iorder .lt. 3) iorder = 3

                  write(chme,'(i5.5)') me

                  filnm = chfn(10)(1:ilfn(10))
     &                  //'.' // chme(6-iorder:5)

                  if ( irestart .eq. 1 ) then

                     open(10, file = filnm(1:ilfn(10)+1+iorder),
     &                    form='unformatted',status = 'unknown',
     &                    access = 'append' )

                  else

                     open(10, file = filnm(1:ilfn(10)+1+iorder),
     &                    form='unformatted',status = 'unknown' )

                  end if

               end if

            end if

      end if

*-----------------------------------------------------------------------
*     Check material
*-----------------------------------------------------------------------

               irerr = 0

            if( mxmat .eq. 0 ) then

                  write(io,'(/"*** Error no [material] section")')
                  ErrCha = ''
                  MsgID = 'L:5445/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/"*** Error no [material] section")')

                  goto 999

            end if

            do i = 1, iregn

               if( idmg(i) .gt. 0 ) then

                  if( idnm(idmg(i)) .eq. 0 ) then

                     write(io,'("*** ERROR : undefined ",
     &                   "material ID number =",i5,
     &                   " in region/cell = ",i6)') idmg(i), idrg(i)

                     ErrCha = ''
                     MsgID = 'L:5464/R:setpar/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'("*** ERROR : undefined ",
     &                   "material ID number =",i5,
     &                   " in region/cell = ",i6)') idmg(i), idrg(i)

                     irerr = irerr + 1

                  end if

               end if

            end do

               if( irerr .gt. 0 ) then

                  write(io,'(/"*** Error in [region/cell] section",
     &                        " listed above, stop!!"/
     &                        "    Number of errors is ",i3/)') irerr

                  ErrCha = ''
                  MsgID = 'L:5485/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/"*** Error in [region/cell] section",
     &                        " listed above, stop!!"/
     &                        "    Number of errors is ",i3/)') irerr

                  goto 999

               end if

*-----------------------------------------------------------------------
*        read material data
*-----------------------------------------------------------------------

                     nmhigh  = 0
                     nminth  = mmmax
                     mxmat0  = mxmat

                     rewind iom1

               do m = 1, mxmat

                     read(iom1) nel
                     read(iom1) denh, libh
                     read(iom1) igas, istp, inlb, iplb, ielb, icnd
                     read(iom1) iulb, ihlb
                     read(iom1) imts
                     read(iom1) (idedx(ii),ii=1,20)

                     know = nel * 4 + 33

                  num_kmat_add = 33 + nel*4 + imts*3
                  call moddas_reallocate_dbl(
     &                    kvlmax, m, num_kmat_add, kmat, das_kmat)
!
                  das_kmat(kmat(m)+1) = dble( nel )
                  das_kmat(kmat(m)+2) = denh
                  das_kmat(kmat(m)+3) = dble( libh )
                  das_kmat(kmat(m)+4) = dble( igas )
                  das_kmat(kmat(m)+5) = dble( istp )
                  das_kmat(kmat(m)+6) = dble( inlb )
                  das_kmat(kmat(m)+7) = dble( iplb )
                  das_kmat(kmat(m)+8) = dble( ielb )
                  das_kmat(kmat(m)+9) = dble( icnd )
                  das_kmat(kmat(m)+10) = dble( imts )
                  das_kmat(kmat(m)+11) = dble( know )
                  das_kmat(kmat(m)+12) = dble( iulb )
                  das_kmat(kmat(m)+13) = dble( ihlb )
!
                  do ii=1,20
                     das_kmat(kmat(m)+13+ii) = dble( idedx(ii) )
                  enddo

                  do l = 1, nel
                     read(iom1) icha, masi, denst, libi
                     das_kmat(kmat(m)+(l-1)*4+34) = dble( icha )
                     das_kmat(kmat(m)+(l-1)*4+35) = dble( masi )
                     das_kmat(kmat(m)+(l-1)*4+36) = denst
                     das_kmat(kmat(m)+(l-1)*4+37) = dble( libi )
                  end do

                  do l = 1, imts
                     read(iom1) ix(1), ix(2), ix(3)
                     das_kmat(kmat(m)+know+(l-1)*3+1) = dble( ix(1) )
                     das_kmat(kmat(m)+know+(l-1)*3+2) = dble( ix(2) )
                     das_kmat(kmat(m)+know+(l-1)*3+3) = dble( ix(3) )
                  end do

               end do

*-----------------------------------------------------------------------
*        check materials with different density
*        add new material number
*-----------------------------------------------------------------------

         if(mstz(145).eq.1) then ! T.Sato 2020/09/11
           if(iregn.gt.3000) then
             write(*,'("Number of cell ",i10, " is so large that",
     &       " it is better to set ichkmat = 0")') iregn
           endif
         endif

            do i = 1, iregn - 1

               if( idmg(i) .gt. 0 ) then

                     n = idnm( idmg(i) )

                  do j = i + 1, iregn

                     if( idmg(j) .gt. 0 ) then

                        if( idmg(i) .eq. idmg(j) .and.
     &                      abs( deng(i) - deng(j) ) .gt. 1.d-8 )
     &                  then

                              mxmat = mxmat + 1
                              m = mxmat

                              do k = 1, 99999
                                    if( idnm(k) .eq. 0 ) goto 335
                              end do
                              write(io,'(/
     &                        "*** Error; material number is lack,",
     &                        " for the additional material ",
     &                        " with differnt density")')
                              ErrCha = ''
                              MsgID = 'L:5592/R:setpar/F:read00.f'
                              call ErrWrite(MsgID, ErrCha)
                              write(jo,'(/
     &                        "*** Error; material number is lack,",
     &                        " for the additional material ",
     &                        " with differnt density")')
                                 goto 999
  335                         continue

                              idmg(j) = k
         if(mstz(145).eq.1) then ! T.Sato 2020/09/11
          do jsearch = i+2, iregn
           if( idmg(i) .eq. idmg(jsearch) .and.
     &         abs( deng(j) - deng(jsearch) ) .le. 1.d-8 ) then
              idmg(jsearch) = idmg(j)
           end if
          end do
         endif
                              idnm(k) = m

                           if( matadd .eq. 0 ) then

                              idmn(m) = idmg(i)

                              write(io,'(/
     &                        "*** Warning; same material but",
     &                        " density is different"/
     &                        " Please use matadd = 1",2i6)')
     &                        idmg(i), idmg(j)

                           else

                              idmn(m) = k

                           end if


                           nel1 = nint( das_kmat(kmat(n)+ 1) )
                           imt1 = nint( das_kmat(kmat(n)+10) )
                           know = nint( das_kmat(kmat(n)+11) )

                           num_kmat_add = 33 + nel1*4 + imt1*3
                           call moddas_reallocate_dbl(
     &                        kvlmax, m, num_kmat_add, kmat, das_kmat)

                           das_kmat(kmat(m)+ 1) = das_kmat(kmat(n)+ 1)
                           das_kmat(kmat(m)+ 2) = das_kmat(kmat(n)+ 2)
                           das_kmat(kmat(m)+ 3) = das_kmat(kmat(n)+ 3)
                           das_kmat(kmat(m)+ 4) = das_kmat(kmat(n)+ 4)
                           das_kmat(kmat(m)+ 5) = das_kmat(kmat(n)+ 5)
                           das_kmat(kmat(m)+ 6) = das_kmat(kmat(n)+ 6)
                           das_kmat(kmat(m)+ 7) = das_kmat(kmat(n)+ 7)
                           das_kmat(kmat(m)+ 8) = das_kmat(kmat(n)+ 8)
                           das_kmat(kmat(m)+ 9) = das_kmat(kmat(n)+ 9)
                           das_kmat(kmat(m)+10) = das_kmat(kmat(n)+10)
                           das_kmat(kmat(m)+11) = das_kmat(kmat(n)+11)
                           das_kmat(kmat(m)+12) = das_kmat(kmat(n)+12)
                           das_kmat(kmat(m)+13) = das_kmat(kmat(n)+13)

                           do ii=1,20
                              das_kmat(kmat(m)+13+ii)
     &                           = das_kmat(kmat(n)+13+ii)
                           enddo

                           if( matadd .eq. 1 ) then
                              write(iom1) nint( das_kmat(kmat(m)+ 1) )
                              write(iom1)       das_kmat(kmat(m)+ 2),
     &                                    nint( das_kmat(kmat(m)+ 3) )
                              write(iom1) nint( das_kmat(kmat(m)+ 4) ),
     &                                    nint( das_kmat(kmat(m)+ 5) ),
     &                                    nint( das_kmat(kmat(m)+ 6) ),
     &                                    nint( das_kmat(kmat(m)+ 7) ),
     &                                    nint( das_kmat(kmat(m)+ 8) ),
     &                                    nint( das_kmat(kmat(m)+ 9) )
                              write(iom1) nint( das_kmat(kmat(m)+12) ),
     &                                    nint( das_kmat(kmat(m)+13) )
                              write(iom1) nint( das_kmat(kmat(m)+10) )
                              write(iom1) nint( das_kmat(kmat(m)+14) ),
     &                                    nint( das_kmat(kmat(m)+15) ),
     &                                    nint( das_kmat(kmat(m)+16) ),
     &                                    nint( das_kmat(kmat(m)+17) ),
     &                                    nint( das_kmat(kmat(m)+18) ),
     &                                    nint( das_kmat(kmat(m)+19) ),
     &                                    nint( das_kmat(kmat(m)+20) ),
     &                                    nint( das_kmat(kmat(m)+21) ),
     &                                    nint( das_kmat(kmat(m)+22) ),
     &                                    nint( das_kmat(kmat(m)+23) ),
     &                                    nint( das_kmat(kmat(m)+24) ),
     &                                    nint( das_kmat(kmat(m)+25) ),
     &                                    nint( das_kmat(kmat(m)+26) ),
     &                                    nint( das_kmat(kmat(m)+27) ),
     &                                    nint( das_kmat(kmat(m)+28) ),
     &                                    nint( das_kmat(kmat(m)+29) ),
     &                                    nint( das_kmat(kmat(m)+30) ),
     &                                    nint( das_kmat(kmat(m)+31) ),
     &                                    nint( das_kmat(kmat(m)+32) ),
     &                                    nint( das_kmat(kmat(m)+33) )
                           end if

                           do l = 1, nel1
                              das_kmat(kmat(m)+(l-1)*4+34)
     &                           = das_kmat(kmat(n)+(l-1)*4+34)
                              das_kmat(kmat(m)+(l-1)*4+35)
     &                           = das_kmat(kmat(n)+(l-1)*4+35)
                              das_kmat(kmat(m)+(l-1)*4+36)
     &                           = das_kmat(kmat(n)+(l-1)*4+36)
                              das_kmat(kmat(m)+(l-1)*4+37)
     &                           = das_kmat(kmat(n)+(l-1)*4+37)
                              if( matadd .eq. 1 ) then
                                 write(iom1)
     &                              nint( das_kmat(kmat(m)
     &                                             +(l-1)*4+34) ),
     &                              nint( das_kmat(kmat(m)
     &                                             +(l-1)*4+35) ),
     &                                    das_kmat(kmat(m)
     &                                             +(l-1)*4+36),
     &                              nint( das_kmat(kmat(m)
     &                                             +(l-1)*4+37) )
                              end if
                           end do

                           do l = 1, imt1
                              das_kmat(kmat(m)+know+(l-1)*3+1)
     &                           = das_kmat(kmat(n)+know+(l-1)*3+1)
                              das_kmat(kmat(m)+know+(l-1)*3+2)
     &                           = das_kmat(kmat(n)+know+(l-1)*3+2)
                              das_kmat(kmat(m)+know+(l-1)*3+3)
     &                           = das_kmat(kmat(n)+know+(l-1)*3+3)
                              if( matadd .eq. 1 ) then
                                 write(iom1)
     &                              nint( das_kmat(kmat(m)
     &                                             +know+(l-1)*3+1) ),
     &                              nint( das_kmat(kmat(m)
     &                                             +know+(l-1)*3+2) ),
     &                              nint( das_kmat(kmat(m)
     &                                             +know+(l-1)*3+3) )
                              end if
                           end do

                        end if
                     end if
                  end do
               end if
            end do

            if( matadd .eq. 1 ) mxmat0 = mxmat

*-----------------------------------------------------------------------
                     call dedx_file_allocate(mxmat) ! move to here, T.Sato 2020/10/22
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        check special nuclei ID like 1011
*-----------------------------------------------------------------------

         do m = 1, mxmat

               nel = nint( das_kmat(kmat(m)+1) )

            do l = 1, nel

               icha = nint( das_kmat(kmat(m)+(l-1)*4+34) )
               masi = nint( das_kmat(kmat(m)+(l-1)*4+35) )
               if( icha .eq. 1 .and. masi .gt. 3 ) then

                  write(io,'(/"**** Warning : Strange H isotope is",
     &            " appeared."/
     &            "             : We change this nucleus to"
     &            " 1001,1H Hydrogen,"/
     &            "             : for High Energy Calculation."/
     &            "             : Strange H =",i7/)') icha*1000+masi

                  write(6,'(/"**** Warning : Strange H isotope is",
     &            " appeared."/
     &            "             : We change this nucleus to"
     &            " 1001,1H Hydrogen,"/
     &            "             : for High Energy Calculation."/
     &            "             : Strange H =",i7/)') icha*1000+masi

                  das_kmat(kmat(m)+(l-1)*4+35) = 1
               end if

            end do

         end do

*-----------------------------------------------------------------------
*        rewrite information and expand the natural nuclei
*-----------------------------------------------------------------------

         call readnatural ! T.Sato 2023/12/26 to read natural abundance data from natural_abundance.dat

               iom2 = 24
               iom3 = 25

               open(iom2,form='unformatted',status='scratch')
               open(iom3,form='unformatted',status='scratch')

            do m = 1, mxmat

                        dnel  = das_kmat(kmat(m)+1)
                        denh  = das_kmat(kmat(m)+2)
                        nel1 = nint( dnel )
                        nel2 = 0

               do l = 1, nel1

                        dicha  = das_kmat(kmat(m)+(l-1)*4+34)
                        dmasi  = das_kmat(kmat(m)+(l-1)*4+35)
                        denst  = das_kmat(kmat(m)+(l-1)*4+36)

               if( denst .ne. 0.0d0 ) then

                        nz = nint( dicha )
                        na = nint( dmasi )

                  if(nz.eq.0) then ! T.Sato 2024/12/12, probably O (oxygen) is written as 0 (zero)
                   write(ErrCha,'("More than two numbers are defined",
     &             " as [material] fraction. If you define O (oxygen)",
     &             ", you may mistakenly entered the number 0 (zero)",
     &             " instead of the letter 0")')
                   MsgID = 'L:5813/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   goto 999
                  endif

                  if( na .eq. 0 .and. natnn(nz) .ne. 0 ) then

                     if( denst .gt. 0.0 ) then

                        do k = 1, natnn(nz)

                           if( nz .eq. 1 .and. natnm(nz,k) .eq. 1 ) then

                              denh = denh + denst * patnn(nz,k) / 100.0

                           else

                              nel2 = nel2 + 1

                              dmasim = dble( natnm(nz,k) )
                              denstm = denst * patnn(nz,k) / 100.0

                              write(iom2) dicha, dmasim, denstm

                           end if

                        end do

                     else

                           sek = 0.0

                        do k = 1, natnn(nz)

                            nn = natnm(nz,k) - nz
                           sek = sek + patnn(nz,k) * weitn(nz,nn)

                        end do

                        do k = 1, natnn(nz)
                           if( nz .eq. 1 .and. natnm(nz,k) .eq. 1 ) then
                              denh = denh + denst / sek
     &                                    * patnn(nz,k) * weitn(1,0)
                           else

                              nel2 = nel2 + 1
                                nn = natnm(nz,k) - nz

                              dmasim = dble( natnm(nz,k) )
                              denstm = denst / sek
     &                               * patnn(nz,k) * weitn(nz,nn)

                              write(iom2) dicha, dmasim, denstm

                           end if

                        end do

                     end if

                  else if( na .eq. 0 .and. natnn(nz) .eq. 0 ) then

                   write(io,'("*** ERROR : Nucleus Z =",i3
     &              " is not a natural element.",/12x,
     &              "Please also specify the mass number.")') nz

                   ErrCha = ''
                   MsgID = 'L:5880/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'("*** ERROR : Nucleus Z =",i3
     &              " is not a natural element.",/12x,
     &              "Please also specify the mass number.")') nz

                   goto 999

C S.H. added the error below to stop calc. with . (2017.7.3)
                  else if( ( ndedx .eq. 1 .or. ndedx .eq. 3 )
     &                       .and. nz .gt. 97 ) then

                   write(io,'("*** ERROR : Nucleus Z =",i3,
     &              " is included.",/12x,"When Z > 92, ",
     &              "ATIMA (ndedx = 1 or 3) cannot be used.",
     &              /12x,"Please set ndedx = 2 in [parameters].")') nz

                   ErrCha = ''
                   MsgID = 'L:5898/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'("*** ERROR : Nucleus Z =",i3,
     &              " is included.",/12x,"When Z > 92, ",
     &              "ATIMA (ndedx = 1 or 3) cannot be used.",
     &              /12x,"Please set ndedx = 2 in [parameters].")') nz

                   goto 999

                  else

                        nel2 = nel2 + 1
                        write(iom2) dicha, dmasi, denst

                  end if

               end if
               end do

                        nel = nel2

                        write(iom3) dble(nel), denh

                  if( nel .gt. 0 ) then

                        rewind iom2

                     do i = 1, nel

                         read(iom2) dicha, dmasi, denst
                        write(iom3) dicha, dmasi, denst

                     end do

                        rewind iom2

                  end if

            end do

                     mmmax = nminth

                     rewind iom3

               call moddas_deallocate_dbl(das_kmat)

               call moddas_allocate_dbl(mxmat, dnel_das)
               call moddas_allocate_dbl(mxmat, denh_das)

                     mxnel = 0

               do m = 1, mxmat

                     read(iom3) dnel, denh

                     dnel_das(kmat0+m) = dnel
                     denh_das(kmat0+m) = denh

                     nel0 = nint( dnel )

                     if( nel0 .gt. mxnel ) mxnel = nel0

                     call moddas_reallocate_dbl(
     &                       kvlmax, m, nel0, kmat, zz_das)
                     call moddas_reallocate_dbl(
     &                       kvlmax, m, nel0, kmat, a_das)
                     call moddas_reallocate_dbl(
     &                       kvlmax, m, nel0, kmat, den_das)

                  do l = 1, nel0

                     read(iom3) dicha, dmasi, denst

                     zz_das(kmat(m)+l) = dicha
                     if(dmasi.gt.400) dmasi = dmasi - 400  ! S.H. for meta-stable (2024.1.24)
                     a_das(kmat(m)+l) = dmasi
                     den_das(kmat(m)+l) = denst
                  end do

               end do

                     call moddas_allocate_dbl(mxmat, das_intum)

               do m = 1, mxmat

                     das_intum(intum+m) = -1.0d0

               end do


            close( iom2 )
            close( iom3 )

*-----------------------------------------------------------------------
*        check unused material
*-----------------------------------------------------------------------

            do i = 1, iregn

               if( idmg(i) .gt. 0 ) then

                     m = idnm(idmg(i))

                     das_intum(intum+m) = 1.0d0

               end if

            end do

*-----------------------------------------------------------------------
*     renormalization of density to particle density for unused material
*-----------------------------------------------------------------------

            do m = 1, mxmat

               if( das_intum(intum+m) .lt. 0 ) then
                        nel1 = nint( dnel_das(kmat0+m) )
                        sek = denh_das(kmat0+m)
                        sej = denh_das(kmat0+m)

                     do l = 1, nel1

                        sek = sek + den_das(kmat(m)+l)
                        sej = sej + den_das(kmat(m)+l)*zz_das(kmat(m)+l)

                     end do

                  if( sek .lt. 0 ) then

                        rnrm = -1.0d0

                        denh_das(kmat0+m) = denh_das(m)/weitn(1,0)*rnrm

                     do l = 1, nel1

                              nz  = nint( zz_das(kmat(m)+l) )
                              na  = nint( a_das(kmat(m)+l) )

                        if( na .gt. 0 ) then

                              nn  = na - nz
                              den_das(kmat(m)+l) = den_das(kmat(m)+l)
     &                                             /weitn(nz,nn)*rnrm

                        else if( na .eq. 0 .and. natnn(nz) .ne. 0 ) then

                              sek = 0.0

                           do k = 1, natnn(nz)

                              nn  = natnm(nz,k) - nz
                              sek = sek + den_das(kmat(m)+l)
     &                                    /weitn(nz,nn)*rnrm
     &                                    *patnn(nz,k)/100.0
                           end do

                              den_das(kmat(m)+l) = sek

                        end if

                     end do

                        nel1 = nint( dnel_das(kmat0+m) )
                        sek = denh_das(kmat0+m)
                        sej = denh_das(kmat0+m)

                     do l = 1, nel1

                        sek = sek + den_das(kmat(m)+l)
                        sej = sej + den_das(kmat(m)+l)*zz_das(kmat(m)+l)

                     end do

                  end if

                        denm(m) = sek
                        denc(m) = sej

               end if

            end do


*-----------------------------------------------------------------------
*        set deng(i) for CG case with irden = 0
*-----------------------------------------------------------------------

         if( icgg .eq. 0 .and. irden .eq. 0 ) then

            do i = 1, iregn

                        deng(i) = 0.0d0

               if( idmg(i) .gt. 0 ) then

                        m    = idnm(idmg(i))
                        nel1 = nint( dnel_das(kmat0+m) )
                        sek  = denh_das(kmat0+m)

                     do l = 1, nel1

                        sek = sek + den_das(kmat(m)+l)

                     end do

                        deng(i) = sek

               end if

            end do

         end if

*-----------------------------------------------------------------------
*        renormalization of density to particle density for used mat
*-----------------------------------------------------------------------

               irerr = 0

            do i = 1, iregn

                        denr(i) = 0.0d0

               if( idmg(i) .gt. 0 ) then

                        m    = idnm(idmg(i))
                        nel1 = nint( dnel_das(kmat0+m) )
                        sek  = denh_das(kmat0+m)
                        sej  = denh_das(kmat0+m)

                     do l = 1, nel1

                        sek = sek + den_das(kmat(m)+l)
                        sej = sej + den_das(kmat(m)+l)*zz_das(kmat(m)+l)

                     end do

                  if( sek .eq. 0.0 ) then

                        write(io,'("*** ERROR : material is ",
     &                      "empty MAT = ",i5)') m

                        ErrCha = ''
                        MsgID = 'L:6141/R:setpar/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'("*** ERROR : material is ",
     &                      "empty MAT = ",i5)') m

                        irerr = irerr + 1

                  else if( sek .gt. 0.0 .and. deng(i) .gt. 0.0 ) then

                        rnrm = deng(i) / sek

                        denh_das(kmat0+m) = denh_das(kmat0+m)*rnrm

                     do l = 1, nel1

                        den_das(kmat(m)+l) = den_das(kmat(m)+l)*rnrm

                     end do

                  else if( sek .gt. 0.0 .and. deng(i) .lt. 0.0 ) then

                        sek = denh_das(kmat0+m)*weitn(1,0)

                     do l = 1, nel1

                        nz  = nint( zz_das(kmat(m)+l) )
                        nn  = nint( a_das(kmat(m)+l)
     &                              - zz_das(kmat(m)+l) )
                        sek = sek + den_das(kmat(m)+l)*weitn(nz,nn)

                     end do

                        rnrm = - deng(i) / sek

                        denh_das(kmat0+m) = denh_das(kmat0+m)*rnrm

                     do l = 1, nel1

                        den_das(kmat(m)+l) = den_das(kmat(m)+l)*rnrm

                     end do

                  else if( sek .lt. 0.0 .and. deng(i) .gt. 0.0 ) then

                        sek = denh_das(kmat0+m)/weitn(1,0)

                     do l = 1, nel1

                        nz  = nint( zz_das(kmat(m)+l) )
                        nn  = nint( a_das(kmat(m)+l)
     &                              - zz_das(kmat(m)+l) )
                        sek = sek + den_das(kmat(m)+l)/weitn(nz,nn)

                     end do

                        rnrm = deng(i) / sek

                        denh_das(kmat0+m) = denh_das(kmat0+m)
     &                                      /weitn(1,0)*rnrm

                     do l = 1, nel1

                        nz  = nint( zz_das(kmat(m)+l) )
                        nn  = nint( a_das(kmat(m)+l)
     &                              - zz_das(kmat(m)+l) )
                        den_das(kmat(m)+l) = den_das(kmat(m)+l)
     &                                       /weitn(nz,nn)*rnrm

                     end do

                  else if( sek .lt. 0.0 .and. deng(i) .lt. 0.0 ) then

                        rnrm = - deng(i) / sek

                        denh_das(kmat0+m) = denh_das(kmat0+m)
     &                                      /weitn(1,0)*rnrm

                     do l = 1, nel1

                        nz  = nint( zz_das(kmat(m)+l) )
                        nn  = nint( a_das(kmat(m)+l)
     &                              - zz_das(kmat(m)+l) )
                        den_das(kmat(m)+l) = den_das(kmat(m)+l)
     &                                       /weitn(nz,nn)*rnrm

                     end do

                  end if

                        sek  = denh_das(kmat0+m)
                        sej  = denh_das(kmat0+m)

                     do l = 1, nel1

                        sek = sek + den_das(kmat(m)+l)
                        sej = sej + den_das(kmat(m)+l)
     &                              *zz_das(kmat(m)+l)

                     end do

                        denr(i) = sek
                        denm(m) = sek
                        denc(m) = sej

               end if

            end do

               if( irerr .gt. 0 ) then

                  write(io,'(/"*** Error in [material] section",
     &                        " listed above, stop!!"/
     &                        "    Number of errors is ",i3/)') irerr

                  ErrCha = ''
                  MsgID = 'L:6256/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/"*** Error in [material] section",
     &                        " listed above, stop!!"/
     &                        "    Number of errors is ",i3/)') irerr

                  goto 999

               end if

*-----------------------------------------------------------------------
*           store the information on kmat() for high energy reaction
*
*                    das(kmat(m)+1) = dnel
*                    das(kmat(m)+2) = denh
*
*                    das(kmat(m)+(l-1)*3+3) = dicha
*                    das(kmat(m)+(l-1)*3+4) = dmasi
*                    das(kmat(m)+(l-1)*3+5) = denst
*
*                    das(kmatd(m)+1) = dmaxp1
*                    das(kmatd(m)+2) = dmaxp2
*                    das(kmatd(m)+3) = dmaxn1
*                    das(kmatd(m)+4) = dmaxn2

*                    das(kmatd(m)+(l-1)*2+11) = dmax(1)
*                    das(kmatd(m)+(l-1)*2+12) = dmax(2)
*
*-----------------------------------------------------------------------
cKN 2016/08/01 store dmax for proton and neutron
*-----------------------------------------------------------------------

            do m = 1, mxmat

                     nel = nint( dnel_das(kmat0+m) )
                     denh = denh_das(kmat0+m)

                     num_kmat_add = 10 + nel*5
                     call moddas_reallocate_dbl(
     &                  kvlmax, m, num_kmat_add, kmatd, das_kmatd)
                     num_kmath_add = 6*mxmat
                     call moddas_reallocate_dbl(
     &                  kvlmax, m, num_kmath_add, kmathd, das_kmathd)

                        dmaxp1 = 0.0d0
                        dmaxp2 = 1.d+30
                        dmaxn1 = 0.0d0
                        dmaxn2 = 1.d+30
                        dmaxu1 = 0.0d0
                        dmaxu2 = 1.d+30
                        dmaxd1 = 0.0d0
                        dmaxd2 = 1.d+30
                        dmaxa1 = 0.0d0
                        dmaxa2 = 1.d+30

               do l = 0, nel

                     if( l .eq. 0 ) then

                        iz = 1
                        ia = 1

                     else

                        iz = nint( zz_das(kmat(m)+l) )
                        ia = nint( a_das(kmat(m)+l) )

                     end if

                        das_kmatd(kmatd(m)+(l-1)*5+11) = dnmax(1)
                        das_kmatd(kmatd(m)+(l-1)*5+12) = dnmax(2)
                        das_kmatd(kmatd(m)+(l-1)*5+13) = dpnmax
                        das_kmatd(kmatd(m)+(l-1)*5+14) = dnmax(15)
                        das_kmatd(kmatd(m)+(l-1)*5+15) = dnmax(18)

*-----------------------------------------------------------------------

                        ihmn = 0

                  do i = 1, indmp

                           imm = matdxp(i)

                     if( imm .eq. 0 .or. imm .eq. idmn(m) ) then

                           izi = nucdxp(i) / 1000
                           iai = nucdxp(i) - izi * 1000

                        if( izi .eq.  0 .or.
     &                    ( iai .eq.  0 .and. izi .eq. iz ) .or.
     &                    ( iai .eq. ia .and. izi .eq. iz ) ) then

                           das_kmatd(kmatd(m)+(l-1)*5+11) = dmxdxp(i)

                           if( l .gt. 0 ) then
                           if( dmxdxp(i) .gt. dmaxp1 )
     &                         dmaxp1 = dmxdxp(i)
                           if( dmxdxp(i) .lt. dmaxp2 )
     &                         dmaxp2 = dmxdxp(i)
                           end if

                           ihmn = ihmn + 1

                        end if

                     end if

                  end do

                        if( ihmn .eq. 0 .and.
     &                    ( l .ne. 0 .or. denh .gt. 0.0d0 ) ) then

                           if( dnmax(1) .gt. dmaxp1 )
     &                         dmaxp1 = dnmax(1)
                           if( dnmax(1) .lt. dmaxp2 )
     &                         dmaxp2 = dnmax(1)

                        end if

*-----------------------------------------------------------------------

                        ihmn = 0

                  do i = 1, indmn

                           imm = matdxn(i)

                     if( imm .eq. 0 .or. imm .eq. idmn(m) ) then

                           izi = nucdxn(i) / 1000
                           iai = nucdxn(i) - izi * 1000

                        if( izi .eq.  0 .or.
     &                    ( iai .eq.  0 .and. izi .eq. iz ) .or.
     &                    ( iai .eq. ia .and. izi .eq. iz ) ) then

                           das_kmatd(kmatd(m)+(l-1)*5+12) = dmxdxn(i)

                           if( l .gt. 0 ) then
                           if( dmxdxn(i) .gt. dmaxn1 )
     &                         dmaxn1 = dmxdxn(i)
                           if( dmxdxn(i) .lt. dmaxn2 )
     &                         dmaxn2 = dmxdxn(i)
                           end if

                           ihmn = ihmn + 1

                        end if

                     end if

                  end do

                        if( ihmn .eq. 0 .and.
     &                    ( l .ne. 0 .or. denh .gt. 0.0d0 ) ) then

                           if( dnmax(2) .gt. dmaxn1 )
     &                         dmaxn1 = dnmax(2)
                           if( dnmax(2) .lt. dmaxn2 )
     &                         dmaxn2 = dnmax(2)

                        end if

*-----------------------------------------------------------------------
cfrtati 2021/12/17 data max for photo-nuclear library
                        ihmn = 0

                  do i = 1, indmu

                           imm = matdxu(i)

                     if( imm .eq. 0 .or. imm .eq. idmn(m) ) then

                           izi = nucdxu(i) / 1000
                           iai = nucdxu(i) - izi * 1000

                        if( izi .eq.  0 .or.
     &                    ( iai .eq.  0 .and. izi .eq. iz ) .or.
     &                    ( iai .eq. ia .and. izi .eq. iz ) ) then

                            das_kmatd(kmatd(m)+(l-1)*5+13) = dmxdxu(i)

                           if( l .gt. 0 ) then
                           if( dmxdxu(i) .gt. dmaxu1 )
     &                         dmaxu1 = dmxdxu(i)
                           if( dmxdxu(i) .lt. dmaxu2 )
     &                         dmaxu2 = dmxdxu(i)
                           end if

                           ihmn = ihmn + 1

                        end if

                     end if

                  end do

                        if( ihmn .eq. 0 .and.
     &                    ( l .ne. 0 .or. denh .gt. 0.0d0 ) ) then

                           if( dpnmax .gt. dmaxu1 )
     &                         dmaxu1 = dpnmax
                           if( dpnmax .lt. dmaxu2 )
     &                         dmaxu2 = dpnmax

                        end if

*-----------------------------------------------------------------------
cfrtati 2021/12/17 data max for deuteron library
                        ihmn = 0

                  do i = 1, indmd

                           imm = matdxd(i)

                     if( imm .eq. 0 .or. imm .eq. idmn(m) ) then

                           izi = nucdxd(i) / 1000
                           iai = nucdxd(i) - izi * 1000

                        if( izi .eq.  0 .or.
     &                    ( iai .eq.  0 .and. izi .eq. iz ) .or.
     &                    ( iai .eq. ia .and. izi .eq. iz ) ) then

                            das_kmatd(kmatd(m)+(l-1)*5+14) = dmxdxd(i)

                           if( l .gt. 0 ) then
                           if( dmxdxd(i) .gt. dmaxd1 )
     &                         dmaxd1 = dmxdxd(i)
                           if( dmxdxd(i) .lt. dmaxd2 )
     &                         dmaxd2 = dmxdxd(i)
                           end if

                           ihmn = ihmn + 1

                        end if

                     end if

                  end do

                        if( ihmn .eq. 0 .and.
     &                    ( l .ne. 0 .or. denh .gt. 0.0d0 ) ) then

                           if( dnmax(15) .gt. dmaxd1 )
     &                         dmaxd1 = dnmax(15)
                           if( dnmax(15) .lt. dmaxd2 )
     &                         dmaxd2 = dnmax(15)

                        end if

*-----------------------------------------------------------------------
cfrtati 2021/12/17 data max for alpha library
                        ihmn = 0

                  do i = 1, indma

                           imm = matdxa(i)

                     if( imm .eq. 0 .or. imm .eq. idmn(m) ) then

                           izi = nucdxa(i) / 1000
                           iai = nucdxa(i) - izi * 1000

                        if( izi .eq.  0 .or.
     &                    ( iai .eq.  0 .and. izi .eq. iz ) .or.
     &                    ( iai .eq. ia .and. izi .eq. iz ) ) then

                           das_kmatd(kmatd(m)+(l-1)*5+15) = dmxdxa(i)

                           if( l .gt. 0 ) then
                           if( dmxdxa(i) .gt. dmaxa1 )
     &                         dmaxa1 = dmxdxa(i)
                           if( dmxdxd(i) .lt. dmaxd2 )
     &                         dmaxa2 = dmxdxa(i)
                           end if

                           ihmn = ihmn + 1

                        end if

                     end if

                  end do

                        if( ihmn .eq. 0 .and.
     &                    ( l .ne. 0 .or. denh .gt. 0.0d0 ) ) then

                           if( dnmax(18) .gt. dmaxa1 )
     &                         dmaxa1 = dnmax(18)
                           if( dnmax(18) .lt. dmaxa2 )
     &                         dmaxa2 = dnmax(18)

                        end if

*-----------------------------------------------------------------------

               end do

                        das_kmatd(kmatd(m)+1) = dmaxp1
                        das_kmatd(kmatd(m)+2) = dmaxp2
                        das_kmatd(kmatd(m)+3) = dmaxn1
                        das_kmatd(kmatd(m)+4) = dmaxn2

                        das_kmathd(kmathd(m)+1) = dmaxu1
                        das_kmathd(kmathd(m)+2) = dmaxu2
                        das_kmathd(kmathd(m)+3) = dmaxd1
                        das_kmathd(kmathd(m)+4) = dmaxd2
                        das_kmathd(kmathd(m)+5) = dmaxa1
                        das_kmathd(kmathd(m)+6) = dmaxa2

            end do

*-----------------------------------------------------------------------
*           store the information on kmatg()
*           material = i, nuclei = j
*
*              nel  = nint( das(kmatg(i)+1) )
*              igas = nint( das(kmatg(i)+2) )
*              istp = nint( das(kmatg(i)+3) )
*              inlb = nint( das(kmatg(i)+4) )
*              iplb = nint( das(kmatg(i)+5) )
*              ielb = nint( das(kmatg(i)+6) )
*              icnd = nint( das(kmatg(i)+7) )
*              imts = nint( das(kmatg(i)+8) )
*              know = nint( das(kmatg(i)+9) )
*              iulb = nint( das(kmatg(i)+10) )
*              ihlb = nint( das(kmatg(i)+11) )
*              do ii=1,20
*               idedx(ii) = nint( das(kmatg(i)+11+ii) )
*              enddo
*
*              izia = nint( das(kmatg(i)+(j-1)*3+32) )
*              den  =       das(kmatg(i)+(j-1)*3+33)
*              libi = nint( das(kmatg(i)+(j-1)*3+34) )
*
*              ix(1) = nint( das(kmatg(i)+know+(j-1)*3+1) )
*              ix(2) = nint( das(kmatg(i)+know+(j-1)*3+2) )
*              ix(3) = nint( das(kmatg(i)+know+(j-1)*3+3) )
*
*-----------------------------------------------------------------------
cKN 2016/07/31 expand natural nucleus

                  iom2 = 24
                  iom3 = 25

                  open(iom2,form='unformatted',status='scratch')
                  open(iom3,form='unformatted',status='scratch')

                  rewind iom1

               do m = 1, mxmat0

                     read(iom1) nel
                     read(iom1) denh, libh
                     read(iom1) igas, istp, inlb, iplb, ielb, icnd
                     read(iom1) iulb, ihlb
                     read(iom1) imts
                     read(iom1) (idedx(ii),ii=1,20)

                     nel2 = 0

               do ll = 1, nel

                     read(iom1) nz, na, denst, libi

               if( denst .ne. 0.0d0 ) then

                  if( na .eq. 0 .and. natnn(nz) .ne. 0 .and.
     &                inatur .ne. 0 ) then ! T.Sato 2023/08/30 carbon exception is removed

                     if( denst .gt. 0.0 ) then

                        do k = 1, natnn(nz)

                           if( nz .eq. 1 .and. natnm(nz,k) .eq. 1 ) then

                              denh = denh + denst * patnn(nz,k) / 100.0

                           else

                              nel2 = nel2 + 1

                              masim  = natnm(nz,k)
                              denstm = denst * patnn(nz,k) / 100.0

                              if( nz .eq. 1 ) then
                                 write(iom2) nz, masim, denstm, libh
                              else
                                 write(iom2) nz, masim, denstm, libi
                              end if

                           end if

                        end do

                     else if( denst .le. 0.0 ) then

                           sek = 0.0

                        do k = 1, natnn(nz)

                            nn = natnm(nz,k) - nz
                           sek = sek + patnn(nz,k) * weitn(nz,nn)

                        end do

                        do k = 1, natnn(nz)

                           if( nz .eq. 1 .and. natnm(nz,k) .eq. 1 ) then

                              denh = denh + denst / sek
     &                                    * patnn(nz,k) * weitn(1,0)
                           else

                              nel2 = nel2 + 1
                              nn = natnm(nz,k) - nz

                              masim  = natnm(nz,k)
                              denstm = denst / sek
     &                               * patnn(nz,k) * weitn(nz,nn)

                              if( nz .eq. 1 ) then
                                 write(iom2) nz, masim, denstm, libh
                              else
                                 write(iom2) nz, masim, denstm, libi
                              end if

                           end if

                        end do

                     end if

                  else if( na .eq. 0 .and. natnn(nz) .eq. 0 .and.
     &                     inatur .ne. 0 ) then ! T.Sato 2023/08/30 carbon exception is removed

                   write(io,'("*** ERROR : Nucleus Z =",i3
     &              " is not a natural element.",/12x,
     &              "Please also specify the mass number.")') nz

                   ErrCha = ''
                   MsgID = 'L:6698/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'("*** ERROR : Nucleus Z =",i3
     &              " is not a natural element.",/12x,
     &              "Please also specify the mass number.")') nz

                   goto 999

C S.H. added the error below to stop calc. with . (2017.7.3)
                  else if( ( ndedx .eq. 1 .or. ndedx .eq. 3 )
     &                       .and. nz .gt. 97 ) then

                   write(io,'("*** ERROR : Nucleus Z =",i3,
     &              " is included.",/12x,"When Z > 92, ",
     &              "ATIMA (ndedx = 1 or 3) cannot be used.",
     &              /12x,"Please set ndedx = 2 in [parameters].")') nz

                   ErrCha = ''
                   MsgID = 'L:6716/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'("*** ERROR : Nucleus Z =",i3,
     &              " is included.",/12x,"When Z > 92, ",
     &              "ATIMA (ndedx = 1 or 3) cannot be used.",
     &              /12x,"Please set ndedx = 2 in [parameters].")') nz

                   goto 999

                  else

                        nel2 = nel2 + 1
                        write(iom2) nz, na, denst, libi

                  end if

            end if
            end do

               if( imts .gt. 0 ) then
                  do idum = 1, imts
                     read(iom1) idum1, idum2, idum3
                  enddo
               endif

                        nel = nel2

                        write(iom3) nel, denh

                  if( nel .gt. 0 ) then

                        rewind iom2

                     do i = 1, nel

                         read(iom2) nz, na, denst, libi
                        write(iom3) nz, na, denst, libi

                     end do

                        rewind iom2

                  end if

            end do

                  rewind iom3

*-----------------------------------------------------------------------

                     rewind iom1

               do m = 1, mxmat0

                     read(iom1) nel0
                     read(iom1) denh, libh
                     read(iom1) igas, istp, inlb, iplb, ielb, icnd
                     read(iom1) iulb, ihlb
                     read(iom1) imts
                     read(iom1) (idedx(ii),ii=1,20)

*-----------------------------------------------------------------------
cKN 2016/07/31 expand natural nucleus

                     read(iom3) nel, denh

*-----------------------------------------------------------------------


                     nelp = 0
                     nel3 = nel

                  if( denh .ne. 0.0d0 ) then

                     nel  = nel + 1
                     nelp = 1

                  end if

                     if( nel .gt. mxnel ) mxnel = nel

                     know = nel * 3 + 31

                  num_kmat_add = 31 + nel*3 + imts*3
                  call moddas_reallocate_dbl(
     &                    kvlmax, m, num_kmat_add, kmatg, das_kmatg)

                     das_kmatg(kmatg(m)+1) = dble( nel )
                     das_kmatg(kmatg(m)+2) = dble( igas )
                     das_kmatg(kmatg(m)+3) = dble( istp )
                     das_kmatg(kmatg(m)+4) = dble( inlb )
                     das_kmatg(kmatg(m)+5) = dble( iplb )
                     das_kmatg(kmatg(m)+6) = dble( ielb )
                     das_kmatg(kmatg(m)+7) = dble( icnd )
                     das_kmatg(kmatg(m)+8) = dble( imts )
                     das_kmatg(kmatg(m)+9) = dble( know )
                     das_kmatg(kmatg(m)+10) = dble( iulb )
                     das_kmatg(kmatg(m)+11) = dble( ihlb )

                  do ii=1, 20
                     das_kmatg(kmatg(m)+11+ii) = dble( idedx(ii) )
                  enddo

                  if( nelp .eq. 1 ) then
                     l = 1
                     das_kmatg(kmatg(m)+(l-1)*3+32) = 1001.0d0
                     das_kmatg(kmatg(m)+(l-1)*3+33) = denh
                     das_kmatg(kmatg(m)+(l-1)*3+34) = dble( libh )
                  end if

                  do ll = 1, nel0
                     read(iom1) icha, masi, denst, libi
                  end do

                  do ll = 1, nel3

                     l = ll + nelp

*-----------------------------------------------------------------------
cKN 2016/07/31 expand natural nucleus

                     read(iom3) icha, masi, denst, libi

*-----------------------------------------------------------------------

                     das_kmatg(kmatg(m)+(l-1)*3+32) = dble( icha )
     &                                                *1000.0
     &                                                + dble( masi )
                     das_kmatg(kmatg(m)+(l-1)*3+33) = denst
                     das_kmatg(kmatg(m)+(l-1)*3+34) = dble( libi )

                  end do

                  do l = 1, imts

                     read(iom1) ix(1), ix(2), ix(3)

                     das_kmatg(kmatg(m)+know+(l-1)*3+1) = dble( ix(1) )
                     das_kmatg(kmatg(m)+know+(l-1)*3+2) = dble( ix(2) )
                     das_kmatg(kmatg(m)+know+(l-1)*3+3) = dble( ix(3) )

                  end do

               end do

*-----------------------------------------------------------------------
cKN 2016/07/31 expand natural nucleus

            close( iom2 )
            close( iom3 )

*-----------------------------------------------------------------------
cKN 2016/08/01 store dmax for proton and neutron
*-----------------------------------------------------------------------

*                    das(kmate(m)+1) = dmaxp1
*                    das(kmate(m)+2) = dmaxp2
*                    das(kmate(m)+3) = dmaxn1
*                    das(kmate(m)+4) = dmaxn2
*
*                    das(kmate(m)+(l-1)*2+11) = dmax(1)
*                    das(kmate(m)+(l-1)*2+12) = dmax(2)
*
*-----------------------------------------------------------------------

            do m = 1, mxmat

                     nel = nint( das_kmatg(kmatg(m)+1) )

                     num_kmat_add = 10 + nel*5
                     call moddas_reallocate_dbl(
     &                  kvlmax, m, num_kmat_add, kmate, das_kmate)
                     num_kmath_add = 6*mxmat
                     call moddas_reallocate_dbl(
     &                  kvlmax, m, num_kmath_add, kmathe, das_kmathe)

                        dmaxp1 = 0.0d0
                        dmaxp2 = 1.d+30
                        dmaxn1 = 0.0d0
                        dmaxn2 = 1.d+30
                        dmaxu1 = 0.0d0
                        dmaxu2 = 1.d+30
                        dmaxd1 = 0.0d0
                        dmaxd2 = 1.d+30
                        dmaxa1 = 0.0d0
                        dmaxa2 = 1.d+30

               do l = 1, nel

                        iza = nint( das_kmatg(kmatg(m)+(l-1)*3+32) )
                        iz  = iza / 1000
                        ia  = iza - iz * 1000

                        das_kmate(kmate(m)+(l-1)*5+11) = dnmax(1)
                        das_kmate(kmate(m)+(l-1)*5+12) = dnmax(2)
                        das_kmate(kmate(m)+(l-1)*5+13) = dpnmax
                        das_kmate(kmate(m)+(l-1)*5+14) = dnmax(15)
                        das_kmate(kmate(m)+(l-1)*5+15) = dnmax(18)

                        ihmn = 0

                  do i = 1, indmp

                           imm = matdxp(i)

                     if( imm .eq. 0 .or. imm .eq. idmn(m) ) then

                           izi = nucdxp(i) / 1000
                           iai = nucdxp(i) - izi * 1000

                        if( izi .eq.  0 .or.
     &                    ( iai .eq.  0 .and. izi .eq. iz ) .or.
     &                    ( iai .eq. ia .and. izi .eq. iz ) ) then

                           das_kmate(kmate(m)+(l-1)*5+11) = dmxdxp(i)

                           if( dmxdxp(i) .gt. dmaxp1 )
     &                         dmaxp1 = dmxdxp(i)
                           if( dmxdxp(i) .lt. dmaxp2 )
     &                         dmaxp2 = dmxdxp(i)

                           ihmn = ihmn + 1

                        end if

                     end if

                  end do

                        if( ihmn .eq. 0 ) then

                           if( dnmax(1) .gt. dmaxp1 )
     &                         dmaxp1 = dnmax(1)
                           if( dnmax(1) .lt. dmaxp2 )
     &                         dmaxp2 = dnmax(1)

                        end if

                        ihmn = 0

                  do i = 1, indmn

                           imm = matdxn(i)

                     if( imm .eq. 0 .or. imm .eq. idmn(m) ) then

                           izi = nucdxn(i) / 1000
                           iai = nucdxn(i) - izi * 1000

                        if( izi .eq.  0 .or.
     &                    ( iai .eq.  0 .and. izi .eq. iz ) .or.
     &                    ( iai .eq. ia .and. izi .eq. iz ) ) then

                           das_kmate(kmate(m)+(l-1)*5+12) = dmxdxn(i)
                           if( dmxdxn(i) .gt. dmaxn1 )
     &                         dmaxn1 = dmxdxn(i)
                           if( dmxdxn(i) .lt. dmaxn2 )
     &                         dmaxn2 = dmxdxn(i)

                           ihmn = ihmn + 1

                        end if

                     end if

                  end do

                        if( ihmn .eq. 0 ) then

                           if( dnmax(2) .gt. dmaxn1 )
     &                         dmaxn1 = dnmax(2)
                           if( dnmax(2) .lt. dmaxn2 )
     &                         dmaxn2 = dnmax(2)

                        end if

cfrtati 2021/12/17 data max for photo-nuclear library
                        ihmn = 0

                  do i = 1, indmu

                           imm = matdxu(i)

                     if( imm .eq. 0 .or. imm .eq. idmn(m) ) then

                           izi = nucdxu(i) / 1000
                           iai = nucdxu(i) - izi * 1000

                        if( izi .eq.  0 .or.
     &                    ( iai .eq.  0 .and. izi .eq. iz ) .or.
     &                    ( iai .eq. ia .and. izi .eq. iz ) ) then

                           das_kmate(kmate(m)+(l-1)*5+13) = dmxdxu(i)
                           if( dmxdxu(i) .gt. dmaxu1 )
     &                         dmaxu1 = dmxdxu(i)
                           if( dmxdxu(i) .lt. dmaxu2 )
     &                         dmaxu2 = dmxdxu(i)

                           ihmn = ihmn + 1

                        end if

                     end if

                  end do

                        if( ihmn .eq. 0 ) then

                           if( dpnmax .gt. dmaxu1 )
     &                         dmaxu1 = dpnmax
                           if( dpnmax .lt. dmaxu2 )
     &                         dmaxu2 = dpnmax

                        end if

cfrtati 2021/12/17 data max for deuteron library
                        ihmn = 0

                  do i = 1, indmd

                           imm = matdxd(i)

                     if( imm .eq. 0 .or. imm .eq. idmn(m) ) then

                           izi = nucdxd(i) / 1000
                           iai = nucdxd(i) - izi * 1000

                        if( izi .eq.  0 .or.
     &                    ( iai .eq.  0 .and. izi .eq. iz ) .or.
     &                    ( iai .eq. ia .and. izi .eq. iz ) ) then

                           das_kmate(kmate(m)+(l-1)*5+14) = dmxdxd(i)

                           if( dmxdxd(i) .gt. dmaxd1 )
     &                         dmaxd1 = dmxdxd(i)
                           if( dmxdxd(i) .lt. dmaxd2 )
     &                         dmaxd2 = dmxdxd(i)

                           ihmn = ihmn + 1

                        end if

                     end if

                  end do
                        if( ihmn .eq. 0 ) then

                           if( dnmax(15) .gt. dmaxd1 )
     &                         dmaxd1 = dnmax(15)
                           if( dnmax(15) .lt. dmaxd2 )
     &                         dmaxd2 = dnmax(15)

                        end if

cfrtati 2021/12/17 data max for alpha library
                        ihmn = 0

                  do i = 1, indma

                           imm = matdxa(i)

                     if( imm .eq. 0 .or. imm .eq. idmn(m) ) then

                           izi = nucdxa(i) / 1000
                           iai = nucdxa(i) - izi * 1000

                        if( izi .eq.  0 .or.
     &                    ( iai .eq.  0 .and. izi .eq. iz ) .or.
     &                    ( iai .eq. ia .and. izi .eq. iz ) ) then

                           das_kmate(kmate(m)+(l-1)*5+15) = dmxdxa(i)
                           if( dmxdxa(i) .gt. dmaxa1 )
     &                         dmaxa1 = dmxdxa(i)
                           if( dmxdxa(i) .lt. dmaxa2 )
     &                         dmaxa2 = dmxdxa(i)

                           ihmn = ihmn + 1

                        end if

                     end if

                  end do

                        if( ihmn .eq. 0 ) then

                           if( dnmax(18) .gt. dmaxa1 )
     &                         dmaxa1 = dnmax(18)
                           if( dnmax(18) .lt. dmaxa2 )
     &                         dmaxa2 = dnmax(18)

                        end if

               end do

                        das_kmate(kmate(m)+1) = dmaxp1
                        das_kmate(kmate(m)+2) = dmaxp2
                        das_kmate(kmate(m)+3) = dmaxn1
                        das_kmate(kmate(m)+4) = dmaxn2

                        das_kmathe(kmathe(m)+1) = dmaxu1
                        das_kmathe(kmathe(m)+2) = dmaxu2
                        das_kmathe(kmathe(m)+3) = dmaxd1
                        das_kmathe(kmathe(m)+4) = dmaxd2
                        das_kmathe(kmathe(m)+5) = dmaxa1
                        das_kmathe(kmathe(m)+6) = dmaxa2

            end do

*-----------------------------------------------------------------------
*           store the information on kmatc() for ECHO
*           material = i, nuclei = j
*
*              nel  = nint( das(kmatc(i)+ 1) )
*              denh =       das(kmatc(i)+ 2)
*              libh = nint( das(kmatc(i)+ 3) )
*
*              igas = nint( das(kmatc(i)+ 4) )
*              istp = nint( das(kmatc(i)+ 5) )
*              inlb = nint( das(kmatc(i)+ 6) )
*              iplb = nint( das(kmatc(i)+ 7) )
*              ielb = nint( das(kmatc(i)+ 8) )
*              icnd = nint( das(kmatc(i)+ 9) )
*              imts = nint( das(kmatc(i)+10) )
*              know = nint( das(kmatc(i)+11) )
*              iulb = nint( das(kmatc(i)+12) )
*              ihlb = nint( das(kmatc(i)+13) )
*              do ii=1,20
*               idedx(ii) = nint( das(kmatg(i)+13+ii) )
*              enddo
*
*              iz   = nint( das(kmatc(i)+(j-1)*4+34) )
*              ia   = nint( das(kmatc(i)+(j-1)*4+35) )
*              den  =       das(kmatc(i)+(j-1)*4+36)
*              libi = nint( das(kmatc(i)+(j-1)*4+37) )
*
*              ix(1) = nint( das(kmatc(i)+know+(j-1)*3+1) )
*              ix(2) = nint( das(kmatc(i)+know+(j-1)*3+2) )
*              ix(3) = nint( das(kmatc(i)+know+(j-1)*3+3) )
*
*-----------------------------------------------------------------------
cKN 2016/07/31 expand natural nucleus

                  iom2 = 24
                  iom3 = 25

                  open(iom2,form='unformatted',status='scratch')
                  open(iom3,form='unformatted',status='scratch')

                  rewind iom1

               do m = 1, mxmat0

                     read(iom1) nel
                     read(iom1) denh, libh
                     read(iom1) igas, istp, inlb, iplb, ielb, icnd
                     read(iom1) iulb, ihlb
                     read(iom1) imts
                     read(iom1) (idedx(ii),ii=1,20)

                     nel2 = 0

                  do ll = 1, nel

                     read(iom1) nz, na, denst, libi

                  if( na .eq. 0 .and. natnn(nz) .ne. 0 .and.
     &                inatur .ge. 2 ) then  ! T.Sato 2023/08/31

                     if( denst .gt. 0.0 ) then

                        do k = 1, natnn(nz)

                           if( nz .eq. 1 .and. natnm(nz,k) .eq. 1 ) then

                              denh = denh + denst * patnn(nz,k) / 100.0

                           else

                              nel2 = nel2 + 1

                              masim  = natnm(nz,k)
                              denstm = denst * patnn(nz,k) / 100.0

                              write(iom2) nz, masim, denstm, libi

                           end if

                        end do

                     else

                           sek = 0.0

                        do k = 1, natnn(nz)

                            nn = natnm(nz,k) - nz
                           sek = sek + patnn(nz,k) * weitn(nz,nn)

                        end do

                        do k = 1, natnn(nz)
                           if( nz .eq. 1 .and. natnm(nz,k) .eq. 1 ) then

                              denh = denh + denst / sek
     &                                    * patnn(nz,k) * weitn(1,0)
                           else

                              nel2 = nel2 + 1
                                nn = natnm(nz,k) - nz

                              masim  = natnm(nz,k)
                              denstm = denst / sek
     &                               * patnn(nz,k) * weitn(nz,nn)

                              write(iom2) nz, masim, denstm, libi

                           end if

                        end do

                     end if

                  else if( na .eq. 0 .and. natnn(nz) .eq. 0 .and.
     &                     inatur .ge. 2 ) then   ! T.Sato 2023/08/31

                   write(io,'("*** ERROR : Nucleus Z =",i3
     &              " is not a natural element.",/12x,
     &              "Please also specify the mass number.")') nz

                   ErrCha = ''
                   MsgID = 'L:7247/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'("*** ERROR : Nucleus Z =",i3
     &              " is not a natural element.",/12x,
     &              "Please also specify the mass number.")') nz

                   goto 999

C S.H. added the error below to stop calc. with . (2017.7.3)
                  else if( ( ndedx .eq. 1 .or. ndedx .eq. 3 )
     &                       .and. nz .gt. 97 ) then

                   write(io,'("*** ERROR : Nucleus Z =",i3,
     &              " is included.",/12x,"When Z > 92, ",
     &              "ATIMA (ndedx = 1 or 3) cannot be used.",
     &              /12x,"Please set ndedx = 2 in [parameters].")') nz

                   ErrCha = ''
                   MsgID = 'L:7265/R:setpar/F:read00.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'("*** ERROR : Nucleus Z =",i3,
     &              " is included.",/12x,"When Z > 92, ",
     &              "ATIMA (ndedx = 1 or 3) cannot be used.",
     &              /12x,"Please set ndedx = 2 in [parameters].")') nz

                   goto 999

                  else

                        nel2 = nel2 + 1
                        write(iom2) nz, na, denst, libi

                  end if

               end do

               if( imts .gt. 0 ) then
                  do idum = 1, imts
                     read(iom1) idum1, idum2, idum3
                  enddo
               endif

                        nel = nel2

                        write(iom3) nel, denh

                  if( nel .gt. 0 ) then

                        rewind iom2

                     do i = 1, nel

                         read(iom2) nz, na, denst, libi
                        write(iom3) nz, na, denst, libi

                     end do

                        rewind iom2

                  end if

            end do

                  rewind iom3

*-----------------------------------------------------------------------

                     rewind iom1

               do m = 1, mxmat0

                     read(iom1) nel0
                     read(iom1) denh, libh
                     read(iom1) igas, istp, inlb, iplb, ielb, icnd
                     read(iom1) iulb, ihlb
                     read(iom1) imts
                     read(iom1) (idedx(ii),ii=1,20)

*-----------------------------------------------------------------------
cKN 2016/07/31 expand natural nucleus

                     read(iom3) nel, denh

*-----------------------------------------------------------------------


                     know = nel * 4 + 33

                     num_kmat_add = 33 + nel*4 + imts*3
                     call moddas_reallocate_dbl(
     &                       kvlmax, m, num_kmat_add, kmatc, das_kmatc)

                     das_kmatc(kmatc(m)+ 1) = dble( nel )
                     das_kmatc(kmatc(m)+ 2) = denh
                     das_kmatc(kmatc(m)+ 3) = dble( libh )
                     das_kmatc(kmatc(m)+ 4) = dble( igas )
                     das_kmatc(kmatc(m)+ 5) = dble( istp )
                     das_kmatc(kmatc(m)+ 6) = dble( inlb )
                     das_kmatc(kmatc(m)+ 7) = dble( iplb )
                     das_kmatc(kmatc(m)+ 8) = dble( ielb )
                     das_kmatc(kmatc(m)+ 9) = dble( icnd )
                     das_kmatc(kmatc(m)+10) = dble( imts )
                     das_kmatc(kmatc(m)+11) = dble( know )
                     das_kmatc(kmatc(m)+12) = dble( iulb )
                     das_kmatc(kmatc(m)+13) = dble( ihlb )

                     do ii=1,20
                        das_kmatc(kmatc(m)+13+ii) = dble( idedx(ii) )
                     enddo

                  do ll = 1, nel0
                     read(iom1) icha, masi, denst, libi
                  end do

                  do l = 1, nel

*-----------------------------------------------------------------------
cKN 2016/07/31 expand natural nucleus

                     read(iom3) icha, masi, denst, libi

*-----------------------------------------------------------------------

                     das_kmatc(kmatc(m)+(l-1)*4+34) = dble( icha )
                     das_kmatc(kmatc(m)+(l-1)*4+35) = dble( masi )
                     das_kmatc(kmatc(m)+(l-1)*4+36) = denst
                     das_kmatc(kmatc(m)+(l-1)*4+37) = dble( libi )

                  end do

                  do l = 1, imts

                     read(iom1) ix(1), ix(2), ix(3)

                     das_kmatc(kmatc(m)+know+(l-1)*3+1) = dble( ix(1) )
                     das_kmatc(kmatc(m)+know+(l-1)*3+2) = dble( ix(2) )
                     das_kmatc(kmatc(m)+know+(l-1)*3+3) = dble( ix(3) )

                  end do

               end do

               close( iom1 )

*-----------------------------------------------------------------------
cKN 2016/07/31 expand natural nucleus

            close( iom2 )
            close( iom3 )

*-----------------------------------------------------------------------
*     renormalization of density to particle density
*     for mstz(35) jmout = 1, 2 for echo
*-----------------------------------------------------------------------

      if( mstz(35) .ge. 1 ) then

*-----------------------------------------------------------------------
*           for unused material
*-----------------------------------------------------------------------

            do m = 1, mxmat0

               if( das_intum(intum+m) .lt. 0 ) then
                     nel1 = nint( das_kmatc(kmatc(m)+1) )
                        sek = das_kmatc(kmatc(m)+2)

                     do l = 1, nel1

                        sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)

                     end do

*-----------------------------------------------------------------------

                  if( sek .lt. 0 .and. mstz(35) .eq. 1 ) then

                        sek = das_kmatc(kmatc(m)+2) / weitn(1,0)

                     do l = 1, nel1

                        nz  = nint( das_kmatc(kmatc(m)+(l-1)*4+34) )
                        na  = nint( das_kmatc(kmatc(m)+(l-1)*4+35) )

                        if( na .gt. 0 ) then

                              nn  = na - nz
                              sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)
     &                                    / weitn(nz,nn)

                        else if( na .eq. 0 .and. natnn(nz) .ne. 0 ) then

                           do k = 1, natnn(nz)

                              nn  = natnm(nz,k) - nz
                              sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)
     &                                    / weitn(nz,nn)
     &                                    * patnn(nz,k) / 100.0
                           end do

                        end if

                     end do

                        rnrm = 1.0d0 / sek

                        das_kmatc(kmatc(m)+2) = das_kmatc(kmatc(m)+2)
     &                                          / weitn(1,0) * rnrm

                     do l = 1, nel1

                        nz  = nint( das_kmatc(kmatc(m)+(l-1)*4+34) )
                        na  = nint( das_kmatc(kmatc(m)+(l-1)*4+35) )

                        if( na .gt. 0 ) then

                              nn  = na - nz
                              das_kmatc(kmatc(m)+(l-1)*4+36)
     &                           = das_kmatc(kmatc(m)+(l-1)*4+36)
     &                             / weitn(nz,nn) * rnrm

                        else if( na .eq. 0 .and. natnn(nz) .ne. 0 ) then

                              sek = 0.0

                           do k = 1, natnn(nz)

                              nn  = natnm(nz,k) - nz
                              sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)
     &                                    / weitn(nz,nn) * rnrm
     &                                    * patnn(nz,k) / 100.0
                           end do

                              das_kmatc(kmatc(m)+(l-1)*4+36) = sek

                        end if
                     end do

*-----------------------------------------------------------------------

                  else if( sek .gt. 0 .and. mstz(35) .eq. 2 ) then

                        sek = das_kmatc(kmatc(m)+2) * weitn(1,0)

                     do l = 1, nel1

                        nz  = nint( das_kmatc(kmatc(m)+(l-1)*4+34) )
                        na  = nint( das_kmatc(kmatc(m)+(l-1)*4+35) )

                        if( na .gt. 0 ) then

                              nn  = na - nz
                              sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)
     &                                    * weitn(nz,nn)

                        else if( na .eq. 0 .and. natnn(nz) .ne. 0 ) then

                           do k = 1, natnn(nz)

                              nn  = natnm(nz,k) - nz
                              sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)
     &                                    * weitn(nz,nn)
     &                                    * patnn(nz,k) / 100.0
                           end do

                        end if

                     end do

                        rnrm = -1.0d0 / sek

                        das_kmatc(kmatc(m)+2) = das_kmatc(kmatc(m)+2)
     &                                          * weitn(1,0) * rnrm

                     do l = 1, nel1

                        nz  = nint( das_kmatc(kmatc(m)+(l-1)*4+34) )
                        na  = nint( das_kmatc(kmatc(m)+(l-1)*4+35) )

                        if( na .gt. 0 ) then

                              nn  = na - nz
                              das_kmatc(kmatc(m)+(l-1)*4+36)
     &                           = das_kmatc(kmatc(m)+(l-1)*4+36)
     &                             * weitn(nz,nn) * rnrm

                        else if( na .eq. 0 .and. natnn(nz) .ne. 0 ) then

                              sek = 0.0

                           do k = 1, natnn(nz)

                              nn  = natnm(nz,k) - nz
                              sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)
     &                                    * weitn(nz,nn) * rnrm
     &                                    * patnn(nz,k) / 100.0
                           end do

                              das_kmatc(kmatc(m)+(l-1)*4+36) = sek

                        end if
                     end do

                  end if

*-----------------------------------------------------------------------

               end if
            end do

*-----------------------------------------------------------------------
*        renormalization of density to particle density
*-----------------------------------------------------------------------

            do i = 1, iregn

               if( idmg(i) .gt. 0 ) then

                     m = idnm(idmg(i))

                     nel1 = nint( das_kmatc(kmatc(m)+1) )
                     sek = das_kmatc(kmatc(m)+2)

                     do l = 1, nel1

                        sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)

                     end do

*-----------------------------------------------------------------------

                  if( mstz(35) .eq. 1 ) then

                  if( sek .gt. 0.0 .and. deng(i) .gt. 0.0 ) then

                        rnrm = deng(i) / sek

                        das_kmatc(kmatc(m)+2) = das(kmatc(m)+2) * rnrm

                     do l = 1, nel1

                        das_kmatc(kmatc(m)+(l-1)*4+36)
     &                     = das_kmatc(kmatc(m)+(l-1)*4+36) * rnrm

                     end do

                  else if( sek .gt. 0.0 .and. deng(i) .lt. 0.0 ) then

                        sek = das_kmatc(kmatc(m)+2) * weitn(1,0)

                     do l = 1, nel1

                        nz  = nint( das_kmatc(kmatc(m)+(l-1)*4+34) )
                        na  = nint( das_kmatc(kmatc(m)+(l-1)*4+35) )

                        if( na .gt. 0 ) then

                              nn  = na - nz
                              sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)
     &                                    * weitn(nz,nn)

                        else if( na .eq. 0 .and. natnn(nz) .ne. 0 ) then

                           do k = 1, natnn(nz)

                              nn  = natnm(nz,k) - nz
                              sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)
     &                                    * weitn(nz,nn)
     &                                    * patnn(nz,k) / 100.0

                           end do

                        end if

                     end do

                        rnrm = - deng(i) / sek

                        das_kmatc(kmatc(m)+2) = das_kmatc(kmatc(m)+2)
     &                                          * rnrm

                     do l = 1, nel1

                        das_kmatc(kmatc(m)+(l-1)*4+36)
     *                     = das_kmatc(kmatc(m)+(l-1)*4+36) * rnrm

                     end do

                  else if( sek .lt. 0.0 .and. deng(i) .gt. 0.0 ) then

                        sek = das_kmatc(kmatc(m)+2) / weitn(1,0)

                     do l = 1, nel1

                        nz  = nint( das_kmatc(kmatc(m)+(l-1)*4+34) )
                        na  = nint( das_kmatc(kmatc(m)+(l-1)*4+35) )

                        if( na .gt. 0 ) then

                              nn  = na - nz
                              sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)
     &                                    / weitn(nz,nn)

                        else if( na .eq. 0 .and. natnn(nz) .ne. 0 ) then

                           do k = 1, natnn(nz)

                              nn  = natnm(nz,k) - nz
                              sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)
     &                                    / weitn(nz,nn)
     &                                    * patnn(nz,k) / 100.0

                           end do

                        end if

                     end do

                        rnrm = deng(i) / sek

                        das(kmatc(m)+2) = das_kmatc(kmatc(m)+2)
     &                                    / weitn(1,0) * rnrm

                     do l = 1, nel1

                        nz  = nint( das_kmatc(kmatc(m)+(l-1)*4+34) )
                        na  = nint( das_kmatc(kmatc(m)+(l-1)*4+35) )

                        if( na .gt. 0 ) then

                              nn  = na - nz
                              das_kmatc(kmatc(m)+(l-1)*4+36)
     &                           = das_kmatc(kmatc(m)+(l-1)*4+36)
     &                             / weitn(nz,nn) * rnrm

                        else if( na .eq. 0 .and. natnn(nz) .ne. 0 ) then

                              sek = 0.0

                           do k = 1, natnn(nz)

                              nn  = natnm(nz,k) - nz
                              sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)
     &                                    / weitn(nz,nn) * rnrm
     &                                    * patnn(nz,k) / 100.0

                           end do

                              das_kmatc(kmatc(m)+(l-1)*4+36) = sek

                        end if

                     end do

                  else if( sek .lt. 0.0 .and. deng(i) .lt. 0.0 ) then

                        rnrm = - deng(i) / sek

                        das_kmatc(kmatc(m)+2) = das_kmatc(kmatc(m)+2)
     &                                          / weitn(1,0) * rnrm

                     do l = 1, nel1

                        nz  = nint( das_kmatc(kmatc(m)+(l-1)*4+34) )
                        na  = nint( das_kmatc(kmatc(m)+(l-1)*4+35) )

                        if( na .gt. 0 ) then

                              nn  = na - nz
                              das_kmatc(kmatc(m)+(l-1)*4+36)
     &                           = das_kmatc(kmatc(m)+(l-1)*4+36)
     &                             / weitn(nz,nn) * rnrm

                        else if( na .eq. 0 .and. natnn(nz) .ne. 0 ) then

                              sek = 0.0

                           do k = 1, natnn(nz)

                              nn  = natnm(nz,k) - nz
                              sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)
     &                                    / weitn(nz,nn) * rnrm
     &                                    * patnn(nz,k) / 100.0

                           end do

                              das_kmatc(kmatc(m)+(l-1)*4+36) = sek

                        end if
                     end do

                  end if
                  end if

*-----------------------------------------------------------------------

                  if( mstz(35) .eq. 2 ) then

                  if( sek .lt. 0.0 .and. deng(i) .lt. 0.0 ) then

                        rnrm = deng(i) / sek

                        das_kmatc(kmatc(m)+2) = das_kmatc(kmatc(m)+2)
     &                                          * rnrm

                     do l = 1, nel1

                        das_kmatc(kmatc(m)+(l-1)*4+36)
     &                     = das_kmatc(kmatc(m)+(l-1)*4+36) * rnrm

                     end do

                  else if( sek .lt. 0.0 .and. deng(i) .gt. 0.0 ) then

                        sek = das_kmatc(kmatc(m)+2) / weitn(1,0)

                     do l = 1, nel1

                        nz  = nint( das_kmatc(kmatc(m)+(l-1)*4+34) )
                        na  = nint( das_kmatc(kmatc(m)+(l-1)*4+35) )

                        if( na .gt. 0 ) then

                              nn  = na - nz
                              sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)
     &                                    / weitn(nz,nn)

                        else if( na .eq. 0 .and. natnn(nz) .ne. 0 ) then

                           do k = 1, natnn(nz)

                              nn  = natnm(nz,k) - nz
                              sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)
     &                                    / weitn(nz,nn)
     &                                    * patnn(nz,k) / 100.0

                           end do

                        end if

                     end do

                        rnrm = - deng(i) / sek

                        das_kmatc(kmatc(m)+2) = das_kmatc(kmatc(m)+2)
     &                                          * rnrm

                     do l = 1, nel1

                        das_kmatc(kmatc(m)+(l-1)*4+36)
     &                     = das_kmatc(kmatc(m)+(l-1)*4+36) * rnrm

                     end do

                  else if( sek .gt. 0.0 .and. deng(i) .lt. 0.0 ) then

                        sek = das_kmatc(kmatc(m)+2) * weitn(1,0)

                     do l = 1, nel1

                        nz  = nint( das_kmatc(kmatc(m)+(l-1)*4+34) )
                        na  = nint( das_kmatc(kmatc(m)+(l-1)*4+35) )

                        if( na .gt. 0 ) then

                              nn  = na - nz
                              sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)
     &                                    * weitn(nz,nn)

                        else if( na .eq. 0 .and. natnn(nz) .ne. 0 ) then

                           do k = 1, natnn(nz)

                              nn  = natnm(nz,k) - nz
                              sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)
     &                                    * weitn(nz,nn)
     &                                    * patnn(nz,k) / 100.0

                           end do

                        end if

                     end do

                        rnrm = deng(i) / sek

                        das_kmatc(kmatc(m)+2) = das_kmatc(kmatc(m)+2)
     &                                          * weitn(1,0) * rnrm

                     do l = 1, nel1

                        nz  = nint( das_kmatc(kmatc(m)+(l-1)*4+34) )
                        na  = nint( das_kmatc(kmatc(m)+(l-1)*4+35) )

                        if( na .gt. 0 ) then

                              nn  = na - nz
                              das_kmatc(kmatc(m)+(l-1)*4+36)
     &                           = das_kmatc(kmatc(m)+(l-1)*4+36)
     &                             * weitn(nz,nn) * rnrm

                        else if( na .eq. 0 .and. natnn(nz) .ne. 0 ) then

                              sek = 0.0

                           do k = 1, natnn(nz)

                              nn  = natnm(nz,k) - nz
                              sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)
     &                                    * weitn(nz,nn) * rnrm
     &                                    * patnn(nz,k) / 100.0

                           end do

                              das_kmatc(kmatc(m)+(l-1)*4+36) = sek

                        end if

                     end do

                  else if( sek .gt. 0.0 .and. deng(i) .gt. 0.0 ) then

                        rnrm = - deng(i) / sek

                        das_kmatc(kmatc(m)+2) = das_kmatc(kmatc(m)+2)
     &                                          * weitn(1,0) * rnrm

                     do l = 1, nel1

                        nz  = nint( das_kmatc(kmatc(m)+(l-1)*4+34) )
                        na  = nint( das_kmatc(kmatc(m)+(l-1)*4+35) )

                        if( na .gt. 0 ) then

                              nn  = na - nz
                              das_kmatc(kmatc(m)+(l-1)*4+36)
     &                           = das_kmatc(kmatc(m)+(l-1)*4+36)
     &                             * weitn(nz,nn) * rnrm

                        else if( na .eq. 0 .and. natnn(nz) .ne. 0 ) then

                              sek = 0.0

                           do k = 1, natnn(nz)

                              nn  = natnm(nz,k) - nz
                              sek = sek + das_kmatc(kmatc(m)+(l-1)*4+36)
     &                                    * weitn(nz,nn) * rnrm
     &                                    * patnn(nz,k) / 100.0

                           end do

                              das_kmatc(kmatc(m)+(l-1)*4+36) = sek

                        end if
                     end do

                  end if
                  end if

*-----------------------------------------------------------------------

               end if
            end do

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*        store the weight density rhog
*-----------------------------------------------------------------------
c Takeshi Kai 2017/06/02 (enumpcc(m))

         pmass = 1.6726219d-24

         do m = 1, mxmat

               sekm = 0.0d0
               arho = 0.0d0
               erho = 0.0d0

               apcc = 0.0
               epcc = 0.0

               lem   = nint( dnel_das(kmat0+m) )
               hydro = denh_das(kmat0+m)

            if( hydro .gt. 0.0d0 ) then

               zin = 1.0
               ain = 1.0
               itz = nint( zin )
               itn = nint( ain - zin )

               rhom = hydro * weitn(itz,itn)
               arho = rhom / (dble(itz+itn) * pmass)  ! Takeshi Kai
               erho = arho *  dble(itz)               ! Takeshi Kai

               sekm = sekm + rhom
               apcc = apcc + arho                     ! Takeshi Kai
               epcc = epcc + erho                     ! Takeshi Kai

            end if

            do i = 1, lem

               zin = zz_das(kmat(m)+i)
               ain = a_das(kmat(m)+i)
               itz = nint( zin )
               itn = nint( ain - zin )

               rhom = den_das(kmat(m)+i)*weitn(itz,itn)
               arho = rhom / (dble(itz+itn) * pmass)  ! Takeshi Kai
               erho = arho *  dble(itz)               ! Takeshi Kai

               sekm = sekm + rhom
               apcc = apcc + arho                     ! Takeshi Kai
               epcc = epcc + erho                     ! Takeshi Kai

            end do

               rhog(m) = sekm
               enumpcc(m) = epcc                      ! Takeshi kai

         end do

*-----------------------------------------------------------------------
*        nzztin : maximum atomic number to be stored nuclear yield
*
*        sigg (ksig(m)+(l-1)*2+1) : geometrical cross section
*        siggm(ksig(m)+(l-1)*2+2) : maximum cross section for proton
*-----------------------------------------------------------------------

               nzztin = 0

         do m = 1, mxmat

               nel1 = nint( dnel_das(kmat0+m) )
               call moddas_reallocate_dbl(kvlmax, m, nel1, ksig, sigg)
               call moddas_reallocate_dbl(kvlmax, m, nel1, ksig, siggm)

            do l = 0, nel1

               if( l .eq. 0 )  then

                  nz   = 1
                  na   = 1
                  den = denh_das(kmat0+m)

               else

                  nz  = nint( zz_das(kmat(m)+l) )
                  na  = nint( a_das(kmat(m)+l) )
                  den = den_das(kmat(m)+l)

               end if

               if( nz .gt. nzztin ) nzztin = nz

                  sigme = 0.0
                  ityp  = 1

            if( den .gt. 0.0d0 ) then

               do i = 1, 100

                  ein = 10.0 + dble(i)

                  if( l .eq. 0 ) then

                        kf1 = 2212
                        kf2 = 2212

                        call sigjam(kf1,kf2,ein,sigt,sigel,signe)

                  else

                        call sigrc(ityp,ein,na,nz,sigt,signe,sigel)

                  end if

                  if( ielas .le. 1 ) then
                     sigte = signe
                  else
                     sigte = sigt
                  end if

                  if( sigte .gt. sigme ) sigme = sigte

               end do

            end if

               if( l .gt. 0 ) then

                  sigg(ksig(m)+l) = den * geosig(na)
                  siggm(ksig(m)+l) = den * sigme

               end if

            end do

         end do

*-----------------------------------------------------------------------
*        memory for sigge(l) = das(ksige+l),
*                   signn(l) = das(ksign+l)
*-----------------------------------------------------------------------

               ksige = 0!FURUTA mmmax
               ksign = 0!FURUTA ksige + mxnel

               mxnel_bank=mxnel !FURUTA

*-----------------------------------------------------------------------
*        memory for range of charge particles
*        common /rngmem/ nrange, krnge, krngn, initr
*-----------------------------------------------------------------------

               nrange = 2048

               krnge = 0!FURUTA mmmax
               krngn = 0!FURUTA krnge + nrange
               mxmat_bank=mxmat !FURUTA

*-----------------------------------------------------------------------
*        memory for energy straggling
*        common /estrag/ kesta, kestz, kestd
*-----------------------------------------------------------------------

            if( nedisp .ne. 0 ) then

               kesta = 0!FURUTA mmmax
               kestz = 0!FURUTA kesta + mxmat
               kestd = 0!FURUTA kestz + mxmat

            end if

*-----------------------------------------------------------------------
*        memory for ATIMA
*        common /atima01/ katima
*-----------------------------------------------------------------------

            if( nedisp .eq. 10 .or. nspred .eq. 10 .or.
     &           ndedx .eq. 1  .or.  ndedx .eq. 3 ) then   ! S.Abe 2016/08/08

               katima = 0 !FURUTA mmmax

            end if

*-----------------------------------------------------------------------
*        memory for Delta Ray
*-----------------------------------------------------------------------

            if( mndel .gt. 0 ) then

               call moddas_allocate_dbl(mxmat, edns)

            end if

*-----------------------------------------------------------------------
*        memory for range of heavy ion for SPAR
*        common /sparcmn/ khb, kh1, kh2, kh3, khe(kvlmax), khp
*                         kdf, kdg, kro
*-----------------------------------------------------------------------

         if( ndedx .le. 3 ) then   ! S.Abe 2016/08/08

               khb = 0!FURUTA mmmax
               kh1 = 0!FURUTA khb + mxmat
               kh2 = 0!FURUTA kh1 + mxmat
               kh3 = 0!FURUTA kh2 + mxmat


               maxnl = 0
               mmm0 = 0 !FURUTA
            do m = 1, mxmat

                  khe(m) = mmm0 !FURUTA mmmax

                  nel = nint( dnel_das(kmat0+m) )
                  mmm0=mmm0+nel !FURUTA

               if( nel .gt. maxnl ) maxnl = nel

            end do

            msumnel_bank=mmm0 !FURUTA

                  maxnl = maxnl + 1

               khp = 0!FURUTA mmmax
               kdf = 0!FURUTA khp + 341 * mxmat + 1
               kdg = 0!FURUTA kdf + maxnl
               kro = 0!FURUTA kdg + maxnl

               maxnl_bank=maxnl !FURUTA

         end if

*-----------------------------------------------------------------------
*        memory for Coulomb spreading
*        common /argcns/ kcar, kczs, kcze
*        common /argcns/ arg(),zsmcs(),zemcs()
*-----------------------------------------------------------------------

            if( nspred .ne. 0 ) then

               kcar = 0!FURUTA mmmax
               kczs = 0!FURUTA kcar + mxmat
               kcze = 0!FURUTA kczs + mxmat


            end if

*-----------------------------------------------------------------------
*        memory summary of material of high energy part
*-----------------------------------------------------------------------

               mmmax  = mmmax + 1
               nmfinh = mmmax
               nmhigh = mmmax - nminth

               if( mmmax .gt. mdas ) then

                  write(io,'(/
     &            "<<< Memory ERROR : at material section",
     &                              " for high energy part"/
     &            "*  memory exceeds mdas",/
     &            "*  start of memory =",i9,/
     &            "*    end of memory =",i9,/
     &            "*     total memory =",i9,/
     &            "*       total mdas =",i9,/
     &            "<<<  Please extend mdas in param.inc >>>")')
     &                   nminth, mmmax, nmhigh, mdas

                  ErrCha = ''
                  MsgID = 'L:8185/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/
     &            "<<< Memory ERROR : at material section",
     &                              " for high energy part"/
     &            "*  memory exceeds mdas",/
     &            "*  start of memory =",i9,/
     &            "*    end of memory =",i9,/
     &            "*     total memory =",i9,/
     &            "*       total mdas =",i9,/
     &            "<<<  Please extend mdas in param.inc >>>")')
     &                   nminth, mmmax, nmhigh, mdas

                  goto 999

               end if

*-----------------------------------------------------------------------
*        modify region definition
*-----------------------------------------------------------------------
*           skip double space
*-----------------------------------------------------------------------

         if( icgg .eq. 0 ) then

*-----------------------------------------------------------------------
               mcmx = ( mdas / 2 - 1 ) * 8 + 1
               call moddas_allocate_cha(MAX_NUM_CHRG, chrg)
               call moddas_allocate_cha(MAX_NUM_CHRG, chtm)
               mci = 0
               mct = 0

               rewind iod

               ioe = 22

               open(ioe,status='scratch',form='unformatted')

*-----------------------------------------------------------------------

               do i = 1, iregn

                     read(iod) (chrg(mci+k:mci+k),k=1,ichl(i))
                     chrg(mci+ichl(i)+1:mci+ichmx) = ' '

                     l = 0

                  do k = 1, ichl(i)

                     if( chrg(mci+k:mci+k) .ne. ' ' .or.
     &                   chrg(mci+k+1:mci+k+1) .ne. ' ' ) then

                        l = l + 1
                        chtm(mct+l:mct+l) = chrg(mci+k:mci+k)

                     end if

                  end do

cFURUTA20131226 optimization bug? in gfortran 4.8 can be fixed
cFURUTA20140213 but it does not work in BX900
                  do k = 1, l

                        chrg(mci+k:mci+k) = chtm(mct+k:mct+k)

                  end do

                        ichl(i) = l

*-----------------------------------------------------------------------
*           add + on number
*-----------------------------------------------------------------------

                     l = 0

                     if( deqn1( chrg(mci+1:mci+1) ) ) then

                        l = l + 1
                        chtm(mct+l:mct+l) = '+'

                     end if

                        l = l + 1
                        chtm(mct+l:mct+l) = chrg(mci+1:mci+1)

                  do k = 2, ichl(i)

                     if( chrg(mci+k-1:mci+k-1) .eq. ' ' .and.
     &                   deqn1( chrg(mci+k:mci+k) ) ) then

                        l = l + 1
                        chtm(mct+l:mct+l) = '+'

                     end if

                        l = l + 1
                        chtm(mct+l:mct+l) = chrg(mci+k:mci+k)

                  end do

cFURUTA20131226 optimization bug? in gfortran 4.8 can be fixed
cFURUTA20140213 but it does not work in BX900
                  do k = 1, l

                        chrg(mci+k:mci+k) = chtm(mct+k:mct+k)

                  end do

                        ichl(i) = l

                     write(ioe) (chrg(mci+k:mci+k),k=1,ichl(i))

               end do

*-----------------------------------------------------------------------

               rewind iod
               rewind ioe

               do i = 1, iregn

                     read(ioe)  (chrg(mci+k:mci+k),k=1,ichl(i))
                     write(iod) (chrg(mci+k:mci+k),k=1,ichl(i))

               end do

               close(ioe)

         end if

*-----------------------------------------------------------------------
*        summary for Delta Ray
*-----------------------------------------------------------------------

         if( mndel .gt. 0 ) then

                  iverr = 0

            do i = 1, mndel

               if( idgr( ndels(i) ) .eq. 0 ) then

                  write(io,'(/" **** Error in [Delta Ray], reg =",
     &            i7," is not defined in [region/cell] ",
     &                                 "section")')
     &            ndels(i)

                  ErrCha = ''
                  MsgID = 'L:8333/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/" **** Error in [Delta Ray], reg =",
     &            i7," is not defined in [region/cell] ",
     &                                 "section")')
     &            ndels(i)

                     iverr = iverr + 1

               else

                      delm( idgr( ndels(i) ) ) = rdels(i)

               end if

            end do

                  if( iverr .gt. 0 ) then

                     write(io,'(" **** stop!!, ",
     &               "number of above error is ",i3)') iverr

                     ErrCha = ''
                     MsgID = 'L:8356/R:setpar/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'(" **** stop!!, ",
     &               "number of above error is ",i3)') iverr

                     goto 999

                  end if

*-----------------------------------------------------------------------

            do m = 1, mxmat

                  sekm = 0.0d0

                  lem   = nint( dnel_das(kmat0+m) )
                  hydro = denh_das(kmat0+m)

               if( hydro .gt. 0.0d0 ) sekm = sekm + hydro

               do i = 1, lem

                  zin = zz_das(kmat(m)+i)
                  sekm = sekm + den_das(kmat(m)+i)*zin

               end do

                  edns(kdelt+m) = sekm

            end do


*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*        summary for Track Structure  T.Sato 2017/05/09
*-----------------------------------------------------------------------

           do i = 1, mntsc
             if(ktsc(i).eq.1)then
               call pts_input        ! read database
               call cts_input        ! read database
               exit
             endif
           enddo

         if( mntsc .gt. 0 ) then

                  iverr = 0
                  ierrdb = 0

            matID=0
            matID2=0
            matID1=0
            do i = 1, mntsc  ! check mID parameter
             if(ktsc(i).eq.0) then
              continue
             elseif(ktsc(i).eq.1) then
              if(matID.ne.1)  call water_db(ierrdb) ! read database
              matID=1
chirata etsmode02 20210915
             elseif(ktsc(i).eq.-1) then
               mat0 = idnm(idmg(idgr(ntsc(i))))
               if(ichem(mat0,1).eq.22014 .and. ichem(mat0,2).eq.0) then
                 if(matID2.ne.2) call silicon_db(ierrdb) ! read database
                 matID2=2
               elseif(ichem(mat0,1).eq.1 .and. ichem(mat0,2).eq.0) then
                    if(matID.ne.1) call water_db(ierrdb) ! read database
                    matID=1
               else
                 if(matID1.ne.1) call etsart_db(ierrdb) ! read database
                 matID1=1
               endif
             elseif(ktsc(i).eq.-100) then
              if(matID.ne.1) call water_db(ierrdb) ! read database
              matID=1
             else
              write(io,'(/" **** Error in [Track Structure], mID =",
     &        i7," is not available")') ktsc(i)
              ErrCha = ''
              MsgID = 'L:8438/R:setpar/F:read00.f'
              call ErrWrite(MsgID, ErrCha)
              write(jo,'(/" **** Error in [Track Structure], mID =",
     &        i7," is not available")') ktsc(i)
              iverr = iverr + 1
             endif
            enddo
chirata etsmode02 20211213
            matID = matID + matID2+matID1

            if(ierrdb.ne.0) then ! error in reading data
             if(ierrdb.eq.1) then
              write(io,'(/" **** Error in opening track-structure",
     &        " database, file does not exist")')
              write(io,*) chfn(25)(1:ilfn(25)),'electron/water.dat'
              ErrCha = ''
              MsgID = 'L:8454/R:setpar/F:read00.f'
              call ErrWrite(MsgID, ErrCha)
              write(jo,'(/" **** Error in opening track-structure",
     &        " database, file does not exist")')
              write(jo,*) chfn(25)(1:ilfn(25)),'electron/water.dat'
             endif
             iverr = iverr + 1
            endif

            if(matID.ge.1) then ! track structure mode, check EGS parameter
             if(iegsemi.ne.1.and.icntl.eq.0) then ! EGS should be used
              write(io,'(/"Warning :: Please be sure that track ",
     &       " structure mode is on but secondary electrsons are not ",
     &       " transported because EGS is not used "
     &       "(i.e. consider negs = 1 or 2)")')
              ErrCha = ''
              MsgID = 'L:8470/R:setpar/F:read00.f'
              call ErrWrite(MsgID, ErrCha)
              write(jo,'(/"Warning :: Please be sure that track ",
     &       " structure mode is on but secondary electrsons are not ",
     &       " transported because EGS is not used "
     &       "(i.e. consider negs = 1 or 2)")')
             endif
             if(emin(12).gt.1.0001d-3.and.icntl.eq.0) then
              write(io,'(/"Warning :: Please be sure that track ",
     &       " structure mode is on but emin(12) is higher than 1 keV",
     &       " (i.e. consider emin(12) & emin(13) = 1.0d-3)")')
              ErrCha = ''
              MsgID = 'L:8482/R:setpar/F:read00.f'
              call ErrWrite(MsgID, ErrCha)
              write(jo,'(/"Warning :: Please be sure that track ",
     &       " structure mode is on but emin(12) is higher than 1 keV",
     &       " (i.e. consider emin(12) & emin(13) = 1.0d-3)")')
             endif
            endif

            do i = 1, mntsc

             if( idgr( ntsc(i) ) .eq. 0 ) then

               write(io,'(/" **** Error in [Track Structure], reg =",
     &         i7," is not defined in [region/cell] ",
     &                                 "section")')
     &         ntsc(i)

               ErrCha = ''
               MsgID = 'L:8500/R:setpar/F:read00.f'
               call ErrWrite(MsgID, ErrCha)
               write(jo,'(/" **** Error in [Track Structure], reg =",
     &         i7," is not defined in [region/cell] ",
     &                                 "section")')
     &         ntsc(i)

                     iverr = iverr + 1

             else

c hirata 20211213 for Si track structute (ntscell > 1)
               if(abs(ktsc(i)).le.1) then
                    ntscell( idgr( ntsc(i) ) )
     &                     = ktsc(i)**(KBCflg+1)*(-1)**KBCflg ! if KURBUC.f is absent and ntscell = 1. Convert to ->  -1
               else
                    ntscell( idgr( ntsc(i) ) ) = ktsc(i)
               endif
               ebgets(idgr( ntsc(i) ) ) = 0.d0
               if(bgets(i).ne.0.d0) then
                    ebgets( idgr( ntsc(i) ) ) = bgets(i)
               endif

               ewvets(idgr( ntsc(i) ) ) = 0.d0
               if(wvets(i).ne.0.d0) then
                    ewvets( idgr( ntsc(i) ) ) = wvets(i)
               endif

             end if

            end do

                  if( iverr .gt. 0 ) then

                     write(io,'(" **** stop!!, ",
     &               "number of above error is ",i3)') iverr

                     ErrCha = ''
                     MsgID = 'L:8538/R:setpar/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'(" **** stop!!, ",
     &               "number of above error is ",i3)') iverr

                     goto 999

                  end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*        summary for volume
*-----------------------------------------------------------------------

         if( mnvol .gt. 0 ) then

                  iverr = 0

            do i = 1, mnvol

               if( idgr( nvols(i) ) .eq. 0 ) then

                  write(io,'(/" **** Error in [volume], reg =",
     &            i7," is not defined in [region/cell] ",
     &                                 "section")')
     &            nvols(i)

                  ErrCha = ''
                  MsgID = 'L:8569/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/" **** Error in [volume], reg =",
     &            i7," is not defined in [region/cell] ",
     &                                 "section")')
     &            nvols(i)

                     iverr = iverr + 1

               else

                      dvol( idgr( nvols(i) ) ) = rvols(i)

               end if

            end do

                  if( iverr .gt. 0 ) then

                     write(io,'(" **** stop!!, ",
     &               "number of above error is ",i3)') iverr

                     ErrCha = ''
                     MsgID = 'L:8592/R:setpar/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'(" **** stop!!, ",
     &               "number of above error is ",i3)') iverr

                     goto 999

                  end if

         end if

         if( ivolm .gt. 0 .and. mstz(19) .eq. 0 ) then

               mnvol = iregn

            do i = 1, iregn

                     nvols(i) = idrg(i)
                     rvols(i) = dvol(i)

            end do

         end if

*-----------------------------------------------------------------------
*        summary for temperature
*-----------------------------------------------------------------------

         if( mntmp .gt. 0 ) then

                  iverr = 0

            do i = 1, mntmp

               if( idgr( ntmps(i) ) .eq. 0 ) then

                  write(io,'(/" **** Error in [temperature], reg =",
     &            i7," is not defined in [region/cell] ",
     &                                 "section")')
     &            ntmps(i)

                  ErrCha = ''
                  MsgID = 'L:8634/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/" **** Error in [temperature], reg =",
     &            i7," is not defined in [region/cell] ",
     &                                 "section")')
     &            ntmps(i)

                     iverr = iverr + 1

               else

                      dtmp( idgr( ntmps(i) ) ) = rtmps(i)

               end if

            end do

                  if( iverr .gt. 0 ) then

                     write(io,'(" **** stop!!, ",
     &               "number of above error is ",i3)') iverr

                     ErrCha = ''
                     MsgID = 'L:8657/R:setpar/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'(" **** stop!!, ",
     &               "number of above error is ",i3)') iverr

                     goto 999

                  end if

         end if

*-----------------------------------------------------------------------
*        summary for brems bias
*-----------------------------------------------------------------------

         if( mnbrs .gt. -1 ) then

                  iverr = 0

            if( mnbrs .eq. 0 ) then

                     mnbrs = 0

                  do i = 1, mxmat

                     if( das_intum(intum+m) .gt. 0 ) then

                        mnbrs = mnbrs + 1

                        mbbrs(mnbrs) = idmn(i)

                     end if

                  end do

            else

               if( icbrs .gt. 0 ) then

                     jcbrs = 0

                  do i = 1, mnbrs

                     if( idnm( mbbrs(i) ) .eq. 0 ) then

                        write(io,'(/" **** Error in [brems bias],",
     &                        " mat =",i6,
     &                        " is not defined in [material] ",
     &                        "section")')
     &                        mbbrs(i)

                        ErrCha = ''
                        MsgID = 'L:8709/R:setpar/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(/" **** Error in [brems bias],",
     &                        " mat =",i6,
     &                        " is not defined in [material] ",
     &                        "section")')
     &                        mbbrs(i)

                           iverr = iverr + 1

                     end if

                     if( das_intum(intum+idnm(mbbrs(i))) .gt. 0 ) then

                        jcbrs = jcbrs + 1
                        nbbrs(jcbrs) = mbbrs(i)

                     end if

                  end do

                     if( iverr .gt. 0 ) then

                        write(io,'(" **** stop!!, ",
     &                  "number of above error is ",i3)') iverr

                        ErrCha = ''
                        MsgID = 'L:8736/R:setpar/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **** stop!!, ",
     &                  "number of above error is ",i3)') iverr

                        goto 999

                     end if

               else

                     jcbrs = 0

                  do 522 i = 1, mxmat

                     do j = 1, mnbrs

                        if( mbbrs(j) .eq. idmn(i) ) goto 522

                     end do

                     if( das_intum(intum+i) .gt. 0 ) then

                        jcbrs = jcbrs + 1
                        nbbrs(jcbrs) = idmn(i)

                     end if

  522             continue

               end if

                     mnbrs = jcbrs

                  do i = 1, mnbrs

                     mbbrs(i) = nbbrs(i)

                  end do

            end if

         end if

*-----------------------------------------------------------------------
*        summary for photon weight
*-----------------------------------------------------------------------

         if( mnpwt .gt. 0 ) then

                  iverr = 0

            do i = 1, mnpwt

               if( idgr( npwts(i) ) .eq. 0 ) then

                  write(io,'(/" **** Error in [photon weight], reg =",
     &            i7," is not defined in [region/cell] ",
     &                                 "section")')
     &            npwts(i)

                  ErrCha = ''
                  MsgID = 'L:8798/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/" **** Error in [photon weight], reg =",
     &            i7," is not defined in [region/cell] ",
     &                                 "section")')
     &            npwts(i)

                     iverr = iverr + 1

               else

                      dpwt( idgr( npwts(i) ) ) = rpwts(i)

               end if

            end do

                  if( iverr .gt. 0 ) then

                     write(io,'(" **** stop!!, ",
     &               "number of above error is ",i3)') iverr

                     ErrCha = ''
                     MsgID = 'L:8821/R:setpar/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'(" **** stop!!, ",
     &               "number of above error is ",i3)') iverr

                     goto 999

                  end if

         end if

*-----------------------------------------------------------------------
*        summary of mat name and color
*-----------------------------------------------------------------------
                  nmtnm(-1) = 10
                  dmtnm(-1)(1:10) = 'outer void'
                  nmtcl(-1) = 5
                  dmtcl(-1)(1:5) = 'white'
                  dmhsb(-1,4) = 1.0d0

                  nmtnm(0) = 4
                  dmtnm(0)(1:4) = 'void'
                  nmtcl(0) = 9
                  dmtcl(0)(1:9) = 'lightgray'
                  dmhsb(0,4) = 1.0d0

                  iclr = 0

            do i = 1, mxmat

                  write(mtm,'(i6)') idmn(i)

               do k = 1, 6

                  if( mtm(k:k) .ne. ' ' ) goto 121

               end do

                  k = 1

  121             inm = k

                  nmtnm(i) = 6 - k + 1
                  dmtnm(i)(1:nmtnm(i)) = mtm(inm:6)

                  iclr = iclr + 1
                  if( iclr .gt. 30 ) iclr = 1

                  icll = icolr(iclr)

                  nmtcl(i) = lcoln(icll)
                  dmtcl(i)(1:lcoln(icll)) = colnm(icll)(1:lcoln(icll))

                  dmhsb(i,4) = 1.0d0

            end do

*-----------------------------------------------------------------------

         if( mttcn .gt. 0 ) then

               imnc = 0

            do i = 1, mttcn

                  k = mttc1(i)

                  if( k .gt. 0 ) then

                     j = idnm(k)

                     if( j .eq. 0 ) then

                      write(io,'(/" **** Error in [mat time change],",
     &                        " mat =",i6,
     &                        " is not defined in [material]")')
     &                        k

                      ErrCha = ''
                      MsgID = 'L:8900/R:setpar/F:read00.f'
                      call ErrWrite(MsgID, ErrCha)
                      write(jo,'(/" **** Error in [mat time change],",
     &                        " mat =",i6,
     &                        " is not defined in [material]")')
     &                        k

                           imnc = imnc + 1

                     end if

                  else

                      write(io,'(/" **** Error in [mat time change],",
     &                        " void cannot be defined"
     &                        " as the initial mat")')

                      ErrCha = ''
                      MsgID = 'L:8918/R:setpar/F:read00.f'
                      call ErrWrite(MsgID, ErrCha)
                      write(jo,'(/" **** Error in [mat time change],",
     &                        " void cannot be defined"
     &                        " as the initial mat")')

                           imnc = imnc + 1

                  end if

            end do

                     if( imnc .gt. 0 ) then

                        write(io,'(" **** stop!!, ",
     &                  "number of above error is ",i3)') imnc

                        ErrCha = ''
                        MsgID = 'L:8936/R:setpar/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **** stop!!, ",
     &                  "number of above error is ",i3)') imnc

                        goto 999

                     end if

         end if

*-----------------------------------------------------------------------

         if( mtncn .gt. 0 ) then

               imnc = 0

            do i = 1, mtncn

               do k = mtnc(i,1), mtnc(i,2)

                  if( k .gt. 0 ) then

                     j = idnm(k)

                     if( j .eq. 0 ) then

                       write(io,'(/" **** Error in [mat name color],",
     &                        " mat =",i6,
     &                        " is not defined in [material]")')
     &                        k

                       ErrCha = ''
                       MsgID = 'L:8969/R:setpar/F:read00.f'
                       call ErrWrite(MsgID, ErrCha)
                       write(jo,'(/" **** Error in [mat name color],",
     &                        " mat =",i6,
     &                        " is not defined in [material]")')
     &                        k

                           imnc = imnc + 1

                     else

                        if( nmtnc(i,1) .gt. 0 ) then
                           nmtnm(j) = nmtnc(i,1)
                           dmtnm(j)(1:nmtnm(j)) = dmtnc(i,1)(1:nmtnm(j))
                        end if
                        if( nmtnc(i,2) .gt. 0 ) then
                           nmtcl(j) = nmtnc(i,2)
                           dmtcl(j)(1:nmtcl(j)) = dmtnc(i,2)(1:nmtcl(j))
                        end if

                     end if

                           dmhsb(j,4) = smtnc(i)

                  else

                        if( nmtnc(i,1) .gt. 0 ) then
                           nmtnm(0) = nmtnc(i,1)
                           dmtnm(0)(1:nmtnm(0)) = dmtnc(i,1)(1:nmtnm(0))
                        end if
                        if( nmtnc(i,2) .gt. 0 ) then
                           nmtcl(0) = nmtnc(i,2)
                           dmtcl(0)(1:nmtcl(0)) = dmtnc(i,2)(1:nmtcl(0))
                        end if

                           dmhsb(0,4) = smtnc(i)

                  end if

               end do

            end do

                     if( imnc .gt. 0 ) then

                        write(io,'(" **** stop!!, ",
     &                  "number of above error is ",i3)') imnc

                        ErrCha = ''
                        MsgID = 'L:9018/R:setpar/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **** stop!!, ",
     &                  "number of above error is ",i3)') imnc

                        goto 999

                     end if

         end if

*-----------------------------------------------------------------------
*           get HSB number
*-----------------------------------------------------------------------

            do i = -1, mxmat

                  icc = 0
                  ic  = 0
                  c1  = '['
                  c2  = ']'
                  iclm = nmtcl(i) + 2
                  lum(1)    = c1
                  lum(iclm) = c2

                  do k = 1, nmtcl(i)

                     lum(k+1) = dmtcl(i)(k:k)

                  end do

               call dcols(icc,lum,
     &                    ic,tcol,ierr,iclm,c1,c2)

                  if( ierr .ne. 0 ) then

                        write(io,'(" **** Error: ",
     &                  "color definition is wrong in ",
     &                  "[mat name color] of mat = ",i5)') idmn(i)

                        ErrCha = ''
                        MsgID = 'L:9059/R:setpar/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **** Error: ",
     &                  "color definition is wrong in ",
     &                  "[mat name color] of mat = ",i5)') idmn(i)

                        goto 999

                  end if

                  if( tcol(1) .gt. 0.0 )
     &               tcol(1) = 1.0 - ( tcol(1) - 1.0 ) / 2.5

                  if( tcol(1) .lt. 0.0 )
     &               tcol(1) = -2.0 - tcol(1)

                  dmhsb(i,1) = tcol(1)
                  dmhsb(i,2) = tcol(2)
                  dmhsb(i,3) = tcol(3)

            end do

*-----------------------------------------------------------------------
*        summary of elastic option
*-----------------------------------------------------------------------

         if( mlrgn .gt. 0 ) then

               imnc = 0

            do i = 1, mlrgn

               do k = melrg(i,1), melrg(i,2)

                  if( k .gt. 0 ) then

                     j = idgr(k)

                     if( j .eq. 0 ) then

                      write(io,'(/" **** Error in [elastic option],",
     &                        " reg =",i7,
     &                        " is not defined in [cell]")')
     &                        k

                      ErrCha = ''
                      MsgID = 'L:9105/R:setpar/F:read00.f'
                      call ErrWrite(MsgID, ErrCha)
                      write(jo,'(/" **** Error in [elastic option],",
     &                        " reg =",i7,
     &                        " is not defined in [cell]")')
     &                        k

                           imnc = imnc + 1

                     end if

                  else

                      write(io,'(/" **** Error in [elastic option],",
     &                        " reg =",i7,
     &                        " is not defined in [cell]")')
     &                        k

                      ErrCha = ''
                      MsgID = 'L:9124/R:setpar/F:read00.f'
                      call ErrWrite(MsgID, ErrCha)
                      write(jo,'(/" **** Error in [elastic option],",
     &                        " reg =",i7,
     &                        " is not defined in [cell]")')
     &                        k

                           imnc = imnc + 1

                  end if

               end do

            end do

                     if( imnc .gt. 0 ) then

                        write(io,'(" **** stop!!, ",
     &                  "number of above error is ",i3)') imnc

                        ErrCha = ''
                        MsgID = 'L:9145/R:setpar/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **** stop!!, ",
     &                  "number of above error is ",i3)') imnc

                        goto 999

                     end if

         end if

*-----------------------------------------------------------------------
*     set tally ndata data
*-----------------------------------------------------------------------

         ndatp = 0
         ndatn = 0
         ndatd = 0
         ndata = 0
         ndatg = 0

         mmdatp = 0
         mmdatn = 0
         mmdatd = 0
         mmdata = 0
         mmdatg = 0

         iemmap = 0
         iemman = 0
         iemmad = 0
         iemmaa = 0
         iemmag = 0

         itallo = 0

         indatp(:) = 0
         indatn(:) = 0
         indatd(:) = 0
         indata(:) = 0
         indatg(:) = 0

*-----------------------------------------------------------------------

               jyield = 0
               jndata = 0

*-----------------------------------------------------------------------

      do 400 mm = 1, itnm

         if( ( ital(mm) .ne. 3 .and. ital(mm) .ne. 16 ) .or.
     &       ( itnda(mm) .ne. 2 .and. itnda(mm) .ne. 3 ) ) goto 400

               nl = itmtn(mm)
               ml = itmtt(mm)

*-----------------------------------------------------------------------

         do 401 m = 1, mxmat

*-----------------------------------------------------------------------
*           material choice in t-yield
*-----------------------------------------------------------------------

            if( nl .gt. 0 ) then

                  do i = 1, nl

                     if( itmcn(mm) .gt. 0 .and.
     &                   m .eq. ismte(ml+i-1) ) goto 402
                     if( itmcn(mm) .lt. 0 .and.
     &                   m .eq. ismte(ml+i-1) ) goto 401

                  end do

                     if( itmcn(mm) .gt. 0 ) goto 401

            end if

  402          continue

*-----------------------------------------------------------------------
*        check of each nucleus
*-----------------------------------------------------------------------

               nel = nint( das_kmatg(kmatg(m)+1) )

         do l = 1, nel

              iza = nint( das_kmatg(kmatg(m)+(l-1)*3+32) )
               iz  = iza / 1000
               ia  = iza - iz * 1000

               demaxp = das_kmate(kmate(m)+(l-1)*5+11)
               demaxn = das_kmate(kmate(m)+(l-1)*5+12)
               demaxd = das_kmate(kmate(m)+(l-1)*5+14)
               demaxa = das_kmate(kmate(m)+(l-1)*5+15)
               demaxg = das_kmate(kmate(m)+(l-1)*5+13)

               ipyield = 0
               inyield = 0
               idyield = 0
               iayield = 0
               igyield = 0

               jdcar(1:2) = element(iz)
               write(jdcar(3:5),'(i3.3)') ia

*-----------------------------------------------------------------------
*           read data for proton
*-----------------------------------------------------------------------

                  jdcar(1:2) = element(iz)
                  write(jdcar(3:5),'(i3.3)') ia
                  mltfl = chfn(27)(1:ilfn(27))//jdcar//'-y-p.dat'

                  do kk = 1, ndatp
                     if( iza .eq. indatp(kk) ) goto 403
                  end do

*-----------------------------------------------------------------------

                  inquire( file = mltfl, exist = exex )

               if( exex .eqv. .false. ) then

                  if ( mstz(165) .ge. 2 )
     &                 write(*,'("ndata does not exist ",a60)') mltfl

               else

                  if( itnda(mm) .eq. 2 .or. itnda(mm) .eq. 3 ) then

                     ipyield = 1
                     if( itallo .eq. 0 ) itallo = itallo + 1

                  end if

               end if

*-----------------------------------------------------------------------

         if( ipyield .eq. 1 ) then

                  ndatp = ndatp + 1
                  indatp(ndatp) = iza

                  mdatp = 0
                  ldatp = 0
                  idatp = 2
                  mskip = 1

*-----------------------------------------------------------------------

            iotmp = 15
            open(iotmp, file = mltfl, status = 'old' )

               jsn = 1
               jsi = iotmp
               dsin(1) = mltfl
               idsi(1) = ilfn(27)+11
               ill(1) = 0
               ilf(1) = 10000000
               jpn = 0

 141        continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 996
               if( jpn .eq. 3 ) goto 53
               if( iskip .ne. 0 ) goto 141

            if( chcm(i1:i1+4) .eq. '#elab' ) then

               goto 141

            else if( chcm(i1:i1+3) .eq. '#zap' ) then

               if( mskip .eq. 0 ) then
                  if( ldatp .gt. iemmap ) iemmap = ldatp
               end if

               ic = inumc(chlw,i1+5,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               icl = inumc(chlw,ic,i3,'l') - 1

               call onum(chlw,ic,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 996

                  irzap = nint( cvvv )

                  mskip = 0
                  mdatp = mdatp + 1
                  ldatp = 0

            else
                  ic = i1
               do i = 1, 2
                  call snum(chlw,ic,i3,ic2,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 996
                  remsh(i) = cvvv
                  ic = ic2
               end do

               if( mskip .eq. 0 ) ldatp = ldatp + 1

            end if

                  goto 141

   53       close(iotmp)

               if( mskip .eq. 0 ) then
                  if( ldatp .gt. iemmap ) iemmap = ldatp
               end if

*-----------------------------------------------------------------------

         end if

  403    continue

*-----------------------------------------------------------------------
*           read data for neutron
*-----------------------------------------------------------------------

                  jdcar(1:2) = element(iz)
                  write(jdcar(3:5),'(i3.3)') ia
                  mltfl = chfn(27)(1:ilfn(27))//jdcar//'-y-n.dat'

                  do kk = 1, ndatn
                     if( iza .eq. indatn(kk) ) goto 404
                  end do

*-----------------------------------------------------------------------

                  inquire( file = mltfl, exist = exex )

               if( exex .eqv. .false. ) then

                  if ( mstz(165) .ge. 2 )
     &                 write(*,'("ndata does not exist ",a60)') mltfl

               else

                  if( itnda(mm) .eq. 2 .or. itnda(mm) .eq. 3 ) then

                     inyield = 1
                     if( itallo .eq. 0 ) itallo = itallo + 1

                  end if

               end if

*-----------------------------------------------------------------------

         if( inyield .eq. 1 ) then

                  ndatn = ndatn + 1
                  indatn(ndatn) = iza

                  mdatn = 0
                  ldatn = 0
                  idatn = 2
                  mskip = 1

*-----------------------------------------------------------------------

            iotmp = 15
            open(iotmp, file = mltfl, status = 'old' )

               jsn = 1
               jsi = iotmp
               dsin(1) = mltfl
               idsi(1) = ilfn(27)+11
               ill(1) = 0
               ilf(1) = 10000000
               jpn = 0

 142        continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 996
               if( jpn .eq. 3 ) goto 54
               if( iskip .ne. 0 ) goto 142

            if( chcm(i1:i1+4) .eq. '#elab' ) then

               goto 142

            else if( chcm(i1:i1+3) .eq. '#zap' ) then

               if( mskip .eq. 0 ) then
                  if( ldatn .gt. iemman ) iemman = ldatn
               end if

               ic = inumc(chlw,i1+5,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               icl = inumc(chlw,ic,i3,'l') - 1

               call onum(chlw,ic,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 996

                  irzap = nint( cvvv )

                  mskip = 0
                  mdatn = mdatn + 1
                  ldatn = 0

            else
                  ic = i1
               do i = 1, 2
                  call snum(chlw,ic,i3,ic2,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 996
                  remsh(i) = cvvv
                  ic = ic2
               end do

               if( mskip .eq. 0 ) ldatn = ldatn + 1

            end if

                  goto 142

   54       close(iotmp)

               if( mskip .eq. 0 ) then
                  if( ldatn .gt. iemman ) iemman = ldatn
               end if

*-----------------------------------------------------------------------

         end if

  404    continue

*-----------------------------------------------------------------------
*           read data for deuteron
*-----------------------------------------------------------------------

                  jdcar(1:2) = element(iz)
                  write(jdcar(3:5),'(i3.3)') ia
                  mltfl = chfn(27)(1:ilfn(27))//jdcar//'-y-d.dat'

                  do kk = 1, ndatd
                     if( iza .eq. indatd(kk) ) goto 405
                  end do

*-----------------------------------------------------------------------

                  inquire( file = mltfl, exist = exex )

               if( exex .eqv. .false. ) then

                  if ( mstz(165) .ge. 2 )
     &                 write(*,'("ndata does not exist ",a60)') mltfl

               else

                  if( itnda(mm) .eq. 2 .or. itnda(mm) .eq. 3 ) then

                     idyield = 1
                     if( itallo .eq. 0 ) itallo = itallo + 1

                  end if

               end if

*-----------------------------------------------------------------------

         if( idyield .eq. 1 ) then

                  ndatd = ndatd + 1
                  indatd(ndatd) = iza

                  mdatd = 0
                  ldatd = 0
                  idatd = 2
                  mskip = 1

*-----------------------------------------------------------------------

            iotmp = 15
            open(iotmp, file = mltfl, status = 'old' )

               jsn = 1
               jsi = iotmp
               dsin(1) = mltfl
               idsi(1) = ilfn(27)+11
               ill(1) = 0
               ilf(1) = 10000000
               jpn = 0

 143        continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 996
               if( jpn .eq. 3 ) goto 55
               if( iskip .ne. 0 ) goto 143

            if( chcm(i1:i1+4) .eq. '#elab' ) then

               goto 143

            else if( chcm(i1:i1+3) .eq. '#zap' ) then

               if( mskip .eq. 0 ) then
                  if( ldatd .gt. iemmad ) iemmad = ldatd
               end if

               ic = inumc(chlw,i1+5,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               icl = inumc(chlw,ic,i3,'l') - 1

               call onum(chlw,ic,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 996

                  irzap = nint( cvvv )

                  mskip = 0
                  mdatd = mdatd + 1
                  ldatd = 0

            else
                  ic = i1
               do i = 1, 2
                  call snum(chlw,ic,i3,ic2,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 996
                  remsh(i) = cvvv
                  ic = ic2
               end do

               if( mskip .eq. 0 ) ldatd = ldatd + 1

            end if

                  goto 143

   55       close(iotmp)

               if( mskip .eq. 0 ) then
                  if( ldatd .gt. iemmad ) iemmad = ldatd
               end if

*-----------------------------------------------------------------------

         end if

  405    continue

*-----------------------------------------------------------------------
*           read data for alpha
*-----------------------------------------------------------------------

                  jdcar(1:2) = element(iz)
                  write(jdcar(3:5),'(i3.3)') ia
                  mltfl = chfn(27)(1:ilfn(27))//jdcar//'-y-a.dat'

                  do kk = 1, ndata
                     if( iza .eq. indata(kk) ) goto 406
                  end do

*-----------------------------------------------------------------------

                  inquire( file = mltfl, exist = exex )

               if( exex .eqv. .false. ) then

                  if ( mstz(165) .ge. 2 )
     &                 write(*,'("ndata does not exist ",a60)') mltfl

               else

                  if( itnda(mm) .eq. 2 .or. itnda(mm) .eq. 3 ) then

                     iayield = 1
                     if( itallo .eq. 0 ) itallo = itallo + 1

                  end if

               end if

*-----------------------------------------------------------------------

         if( iayield .eq. 1 ) then

                  ndata = ndata + 1
                  indata(ndata) = iza

                  mdata = 0
                  ldata = 0
                  idata = 2
                  mskip = 1

*-----------------------------------------------------------------------

            iotmp = 15
            open(iotmp, file = mltfl, status = 'old' )

               jsn = 1
               jsi = iotmp
               dsin(1) = mltfl
               idsi(1) = ilfn(27)+11
               ill(1) = 0
               ilf(1) = 10000000
               jpn = 0

 144        continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 996
               if( jpn .eq. 3 ) goto 56
               if( iskip .ne. 0 ) goto 144

            if( chcm(i1:i1+4) .eq. '#elab' ) then

               goto 144

            else if( chcm(i1:i1+3) .eq. '#zap' ) then

               if( mskip .eq. 0 ) then
                  if( ldata .gt. iemmaa ) iemmaa = ldata
               end if

               ic = inumc(chlw,i1+5,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               icl = inumc(chlw,ic,i3,'l') - 1

               call onum(chlw,ic,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 996

                  irzap = nint( cvvv )

                  mskip = 0
                  mdata = mdata + 1
                  ldata = 0

            else
                  ic = i1
               do i = 1, 2
                  call snum(chlw,ic,i3,ic2,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 996
                  remsh(i) = cvvv
                  ic = ic2
               end do

               if( mskip .eq. 0 ) ldata = ldata + 1

            end if

            goto 144

 56         close(iotmp)

               if( mskip .eq. 0 ) then
                  if( ldata .gt. iemmaa ) iemmaa = ldata
               end if

*-----------------------------------------------------------------------

         end if

 406     continue

*-----------------------------------------------------------------------
*           read data for gamma (photon)
*-----------------------------------------------------------------------

                  jdcar(1:2) = element(iz)
                  write(jdcar(3:5),'(i3.3)') ia
                  mltfl = chfn(27)(1:ilfn(27))//jdcar//'-y-g.dat'

                  do kk = 1, ndatg
                     if( iza .eq. indatg(kk) ) goto 407
                  end do

*-----------------------------------------------------------------------

                  inquire( file = mltfl, exist = exex )

               if( exex .eqv. .false. ) then

                  if ( mstz(165) .ge. 2 )
     &                 write(*,'("ndata does not exist ",a60)') mltfl

               else

                  if( itnda(mm) .eq. 2 .or. itnda(mm) .eq. 3 ) then

                     igyield = 1
                     if( itallo .eq. 0 ) itallo = itallo + 1

                  end if

               end if

*-----------------------------------------------------------------------

         if( igyield .eq. 1 ) then

                  ndatg = ndatg + 1
                  indatg(ndatg) = iza

                  mdatg = 0
                  ldatg = 0
                  idatg = 2
                  mskip = 1

*-----------------------------------------------------------------------

            iotmp = 15
            open(iotmp, file = mltfl, status = 'old' )

               jsn = 1
               jsi = iotmp
               dsin(1) = mltfl
               idsi(1) = ilfn(27)+11
               ill(1) = 0
               ilf(1) = 10000000
               jpn = 0

 145        continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 996
               if( jpn .eq. 3 ) goto 57
               if( iskip .ne. 0 ) goto 145

            if( chcm(i1:i1+4) .eq. '#elab' ) then

               goto 145

            else if( chcm(i1:i1+3) .eq. '#zap' ) then

               if( mskip .eq. 0 ) then
                  if( ldatg .gt. iemmag ) iemmag = ldatg
               end if

               ic = inumc(chlw,i1+5,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               icl = inumc(chlw,ic,i3,'l') - 1

               call onum(chlw,ic,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 996

                  irzap = nint( cvvv )

                  mskip = 0
                  mdatg = mdatg + 1
                  ldatg = 0

            else
                  ic = i1
               do i = 1, 2
                  call snum(chlw,ic,i3,ic2,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 996
                  remsh(i) = cvvv
                  ic = ic2
               end do

               if( mskip .eq. 0 ) ldatg = ldatg + 1

            end if

            goto 145

 57         close(iotmp)

               if( mskip .eq. 0 ) then
                  if( ldatg .gt. iemmag ) iemmag = ldatg
               end if

*-----------------------------------------------------------------------

         end if

 407     continue

*-----------------------------------------------------------------------

               if( ipyield.eq.1 .and. mdatp .gt. mmdatp ) mmdatp = mdatp
               if( inyield.eq.1 .and. mdatn .gt. mmdatn ) mmdatn = mdatn
               if( idyield.eq.1 .and. mdatd .gt. mmdatd ) mmdatd = mdatd
               if( iayield.eq.1 .and. mdata .gt. mmdata ) mmdata = mdata
               if( igyield.eq.1 .and. mdatg .gt. mmdatg ) mmdatg = mdatg

*-----------------------------------------------------------------------

c S.H. moved the following block from outside of the 401 do-loop (2021.12.16)
            if( itnda(mm) .eq. 2 .or.  itnda(mm) .eq. 3 ) then
                  jndata = 1
               if( ipyield .ne. 0 ) jyield = 1
               if( inyield .ne. 0 ) jyield = 1
               if( idyield .ne. 0 ) jyield = 1
               if( iayield .ne. 0 ) jyield = 1
               if( igyield .ne. 0 ) jyield = 1
            end if

         end do

  401    continue

  400 continue


            if( jndata .eq. 1 .and. jyield .eq. 0 ) then

                     write(io,'(/" Warning : no cross section",
     &    " data exists in the folder below inspite of ndata = 2,3")')
                     write(io,*) chfn(27)(1:ilfn(27))

                     ErrCha = ''
                     MsgID = 'L:9867/R:setpar/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'(/" Warning : no cross section",
     &    " data exists in the folder below inspite of ndata = 2,3")')
                     write(jo,*) chfn(27)(1:ilfn(27))

            end if

*-----------------------------------------------------------------------
*     re-read tally ndata data
*-----------------------------------------------------------------------

      if( itallo .ne. 0 ) then

         call ALLOCATE_NDATA2


         ndatp0 = 0
         ndatn0 = 0
         ndatd0 = 0
         ndata0 = 0
         ndatg0 = 0

*-----------------------------------------------------------------------

      do 410 mm = 1, itnm

         if( ital(mm) .ne. 3 .or.
     &     ( itnda(mm) .ne. 2 .and. itnda(mm) .ne. 3 ) ) goto 410

               nl = itmtn(mm)
               ml = itmtt(mm)

*-----------------------------------------------------------------------

         do 411 m = 1, mxmat

*-----------------------------------------------------------------------
*           material choice in t-yield
*-----------------------------------------------------------------------

            if( nl .gt. 0 ) then

                  do i = 1, nl

                     if( itmcn(mm) .gt. 0 .and.
     &                   m .eq. ismte(ml+i-1) ) goto 412
                     if( itmcn(mm) .lt. 0 .and.
     &                   m .eq. ismte(ml+i-1) ) goto 411

                  end do

                     if( itmcn(mm) .gt. 0 ) goto 411

            end if

  412          continue

*-----------------------------------------------------------------------
*        check of each nucleus
*-----------------------------------------------------------------------

               nel = nint( das_kmatg(kmatg(m)+1) )

         do l = 1, nel

               iza = nint( das_kmatg(kmatg(m)+(l-1)*3+32) )
               iz  = iza / 1000
               ia  = iza - iz * 1000

               demaxp = das_kmate(kmate(m)+(l-1)*5+11)
               demaxn = das_kmate(kmate(m)+(l-1)*5+12)
               demaxd = das_kmate(kmate(m)+(l-1)*5+14)
               demaxa = das_kmate(kmate(m)+(l-1)*5+15)
               demaxg = das_kmate(kmate(m)+(l-1)*5+13)

               ipyield = 0
               inyield = 0
               idyield = 0
               iayield = 0
               igyield = 0

               jdcar(1:2) = element(iz)
               write(jdcar(3:5),'(i3.3)') ia

*-----------------------------------------------------------------------
*           read data for proton
*-----------------------------------------------------------------------

                  jdcar(1:2) = element(iz)
                  write(jdcar(3:5),'(i3.3)') ia
                  mltfl = chfn(27)(1:ilfn(27))//jdcar//'-y-p.dat'

                  do kk = 1, ndatp0
                     if( iza .eq. indatp(kk) ) goto 413
                  end do

*-----------------------------------------------------------------------

               if( itnda(mm) .eq. 2 .or. itnda(mm) .eq. 3 ) then
                  inquire( file = mltfl, exist = exex )
                  if( exex .eqv. .true. ) ipyield = 1
               end if

*-----------------------------------------------------------------------

         if( ipyield .eq. 1 ) then

                  ndatp0 = ndatp0 + 1

                  mdatp = 0
                  ldatp = 0
                  idatp = 2
                  mskip = 1

*-----------------------------------------------------------------------

            iotmp = 15
            open(iotmp, file = mltfl, status = 'old' )

            if ( mstz(165) .ge. 1 )
     &           write(*,'("Read ndata from ",a60)')  mltfl

               jsn = 1
               jsi = iotmp
               dsin(1) = mltfl
               idsi(1) = ilfn(27)+11
               ill(1) = 0
               ilf(1) = 10000000
               jpn = 0

 151        continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 996
               if( jpn .eq. 3 ) goto 63
               if( iskip .ne. 0 ) goto 151

            if( chcm(i1:i1+4) .eq. '#elab' ) then

               goto 151

            else if( chcm(i1:i1+3) .eq. '#zap' ) then

               if( mskip .eq. 0 ) then
                  kndatp(ndatp0,mdatp) = irzap
                  kldatp(ndatp0,mdatp) = irlip
                  kedatp(ndatp0,mdatp) = ldatp
                  kidatp(ndatp0,mdatp) = idatp
               end if

               ic = inumc(chlw,i1+5,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               icl = inumc(chlw,ic,i3,'l') - 1

               call onum(chlw,ic,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 996

                  irzap = nint( cvvv )

                  mskip = 0
                  mdatp = mdatp + 1
                  ldatp = 0

                  ic = inumc(chlw,icl,i3,'=') + 1
                  ic = jnumc(chlw,ic,i3)
                  icl = inumc(chlw,ic,i3,'i') - 1

                  call onum(chlw,ic,icl,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 998

                  irlip = nint( cvvv )

                  ic = inumc(chlw,icl,i3,'=') + 1
                  ic = jnumc(chlw,ic,i3)

                  call onum(chlw,ic,i3,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 998

                  idatp = nint( cvvv )

            else
                  ic = i1
               do i = 1, 2
                  call snum(chlw,ic,i3,ic2,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 996
                  remsh(i) = cvvv
                  ic = ic2
               end do

               if( mskip .eq. 0 ) then

                  ldatp = ldatp + 1
                   cxproton(ndatp0,mdatp,1,ldatp) = remsh(1)
                   cxproton(ndatp0,mdatp,2,ldatp) = remsh(2)

               end if

            end if

                  goto 151

   63       close(iotmp)

               if( mskip .eq. 0 ) then
                  kndatp(ndatp0,mdatp) = irzap
                  kldatp(ndatp0,mdatp) = irlip
                  kedatp(ndatp0,mdatp) = ldatp
                  kidatp(ndatp0,mdatp) = idatp
               end if

*-----------------------------------------------------------------------

            jdatp(ndatp0) = mdatp

         end if

  413    continue

*-----------------------------------------------------------------------
*          read data for neutron
*-----------------------------------------------------------------------

                  jdcar(1:2) = element(iz)
                  write(jdcar(3:5),'(i3.3)') ia
                  mltfl = chfn(27)(1:ilfn(27))//jdcar//'-y-n.dat'

                  do kk = 1, ndatn0
                     if( iza .eq. indatn(kk) ) goto 414
                  end do

*-----------------------------------------------------------------------

               if( itnda(mm) .eq. 2 .or. itnda(mm) .eq. 3 ) then
                  inquire( file = mltfl, exist = exex )
                  if( exex .eqv. .true. ) inyield = 1
               end if

*-----------------------------------------------------------------------

         if( inyield .eq. 1 ) then

                  ndatn0 = ndatn0 + 1

                  mdatn = 0
                  ldatn = 0
                  idatn = 2
                  mskip = 1

*-----------------------------------------------------------------------

            iotmp = 15
            open(iotmp, file = mltfl, status = 'old' )

            if ( mstz(165) .ge. 1 )
     &           write(*,'("Read ndata from ",a60)')  mltfl

               jsn = 1
               jsi = iotmp
               dsin(1) = mltfl
               idsi(1) = ilfn(27)+11
               ill(1) = 0
               ilf(1) = 10000000
               jpn = 0

 152        continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 996
               if( jpn .eq. 3 ) goto 64
               if( iskip .ne. 0 ) goto 152

            if( chcm(i1:i1+4) .eq. '#elab' ) then

               goto 152

            else if( chcm(i1:i1+3) .eq. '#zap' ) then

               if( mskip .eq. 0 ) then
                  kndatn(ndatn0,mdatn) = irzap
                  kldatn(ndatn0,mdatn) = irlip
                  kedatn(ndatn0,mdatn) = ldatn
                  kidatn(ndatn0,mdatn) = idatn
               end if

               ic = inumc(chlw,i1+5,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               icl = inumc(chlw,ic,i3,'l') - 1

               call onum(chlw,ic,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 996

                  irzap = nint( cvvv )

                  mskip = 0
                  mdatn = mdatn + 1
                  ldatn = 0

                  ic = inumc(chlw,icl,i3,'=') + 1
                  ic = jnumc(chlw,ic,i3)
                  icl = inumc(chlw,ic,i3,'i') - 1

                  call onum(chlw,ic,icl,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 998

                  irlip = nint( cvvv )

                  ic = inumc(chlw,icl,i3,'=') + 1
                  ic = jnumc(chlw,ic,i3)

                  call onum(chlw,ic,i3,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 998

                  idatn = nint( cvvv )

            else
                  ic = i1
               do i = 1, 2
                  call snum(chlw,ic,i3,ic2,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 996
                  remsh(i) = cvvv
                  ic = ic2
               end do

               if( mskip .eq. 0 ) then

                  ldatn = ldatn + 1
                  cxneutron(ndatn0,mdatn,1,ldatn) = remsh(1)
                  cxneutron(ndatn0,mdatn,2,ldatn) = remsh(2)

               end if

            end if

                  goto 152

   64       close(iotmp)

               if( mskip .eq. 0 ) then
                  kndatn(ndatn0,mdatn) = irzap
                  kldatn(ndatn0,mdatn) = irlip
                  kedatn(ndatn0,mdatn) = ldatn
                  kidatn(ndatn0,mdatn) = idatn
               end if

*-----------------------------------------------------------------------

            jdatn(ndatn0) = mdatn

         end if

  414    continue

*-----------------------------------------------------------------------
*          read data for deuteron
*-----------------------------------------------------------------------

                  jdcar(1:2) = element(iz)
                  write(jdcar(3:5),'(i3.3)') ia
                  mltfl = chfn(27)(1:ilfn(27))//jdcar//'-y-d.dat'

                  do kk = 1, ndatd0
                     if( iza .eq. indatd(kk) ) goto 415
                  end do

*-----------------------------------------------------------------------

               if( itnda(mm) .eq. 2 .or. itnda(mm) .eq. 3 ) then
                  inquire( file = mltfl, exist = exex )
                  if( exex .eqv. .true. ) idyield = 1
               end if

*-----------------------------------------------------------------------

         if( idyield .eq. 1 ) then

                  ndatd0 = ndatd0 + 1

                  mdatd = 0
                  ldatd = 0
                  idatd = 2
                  mskip = 1

*-----------------------------------------------------------------------

            iotmp = 15
            open(iotmp, file = mltfl, status = 'old' )

            if ( mstz(165) .ge. 1 )
     &           write(*,'("Read ndata from ",a60)')  mltfl

               jsn = 1
               jsi = iotmp
               dsin(1) = mltfl
               idsi(1) = ilfn(27)+11
               ill(1) = 0
               ilf(1) = 10000000
               jpn = 0

 153        continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 996
               if( jpn .eq. 3 ) goto 65
               if( iskip .ne. 0 ) goto 153

            if( chcm(i1:i1+4) .eq. '#elab' ) then

               goto 153

            else if( chcm(i1:i1+3) .eq. '#zap' ) then

               if( mskip .eq. 0 ) then
                  kndatd(ndatd0,mdatd) = irzap
                  kldatd(ndatd0,mdatd) = irlip
                  kedatd(ndatd0,mdatd) = ldatd
                  kidatd(ndatd0,mdatd) = idatd
               end if

               ic = inumc(chlw,i1+5,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               icl = inumc(chlw,ic,i3,'l') - 1

               call onum(chlw,ic,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 996

                  irzap = nint( cvvv )

                  mskip = 0
                  mdatd = mdatd + 1
                  ldatd = 0

                  ic = inumc(chlw,icl,i3,'=') + 1
                  ic = jnumc(chlw,ic,i3)
                  icl = inumc(chlw,ic,i3,'i') - 1

                  call onum(chlw,ic,icl,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 998

                  irlip = nint( cvvv )

                  ic = inumc(chlw,icl,i3,'=') + 1
                  ic = jnumc(chlw,ic,i3)

                  call onum(chlw,ic,i3,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 998

                  idatd = nint( cvvv )

            else
                  ic = i1
               do i = 1, 2
                  call snum(chlw,ic,i3,ic2,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 996
                  remsh(i) = cvvv
                  ic = ic2
               end do

               if( mskip .eq. 0 ) then

                  ldatd = ldatd + 1
                  cxdeuteron(ndatd0,mdatd,1,ldatd) = remsh(1)
                  cxdeuteron(ndatd0,mdatd,2,ldatd) = remsh(2)

               end if

            end if

                  goto 153

   65       close(iotmp)

               if( mskip .eq. 0 ) then
                  kndatd(ndatd0,mdatd) = irzap
                  kldatd(ndatd0,mdatd) = irlip
                  kedatd(ndatd0,mdatd) = ldatd
                  kidatd(ndatd0,mdatd) = idatd
               end if

*-----------------------------------------------------------------------

            jdatd(ndatd0) = mdatd

         end if

  415    continue

*-----------------------------------------------------------------------
*          read data for alpha
*-----------------------------------------------------------------------

                  jdcar(1:2) = element(iz)
                  write(jdcar(3:5),'(i3.3)') ia
                  mltfl = chfn(27)(1:ilfn(27))//jdcar//'-y-a.dat'

                  do kk = 1, ndata0
                     if( iza .eq. indata(kk) ) goto 416
                  end do

*-----------------------------------------------------------------------

               if( itnda(mm) .eq. 2 .or. itnda(mm) .eq. 3 ) then
                  inquire( file = mltfl, exist = exex )
                  if( exex .eqv. .true. ) iayield = 1
               end if

*-----------------------------------------------------------------------

         if( iayield .eq. 1 ) then

                  ndata0 = ndata0 + 1

                  mdata = 0
                  ldata = 0
                  idata = 2
                  mskip = 1

*-----------------------------------------------------------------------

            iotmp = 15
            open(iotmp, file = mltfl, status = 'old' )

            if ( mstz(165) .ge. 1 )
     &           write(*,'("Read ndata from ",a60)')  mltfl

               jsn = 1
               jsi = iotmp
               dsin(1) = mltfl
               idsi(1) = ilfn(27)+11
               ill(1) = 0
               ilf(1) = 10000000
               jpn = 0

 154        continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 996
               if( jpn .eq. 3 ) goto 66
               if( iskip .ne. 0 ) goto 154

            if( chcm(i1:i1+4) .eq. '#elab' ) then

               goto 154

            else if( chcm(i1:i1+3) .eq. '#zap' ) then

               if( mskip .eq. 0 ) then
                  kndata(ndata0,mdata) = irzap
                  kldata(ndata0,mdata) = irlip
                  kedata(ndata0,mdata) = ldata
                  kidata(ndata0,mdata) = idata
               end if

               ic = inumc(chlw,i1+5,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               icl = inumc(chlw,ic,i3,'l') - 1

               call onum(chlw,ic,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 996

                  irzap = nint( cvvv )

                  mskip = 0
                  mdata = mdata + 1
                  ldata = 0

                  ic = inumc(chlw,icl,i3,'=') + 1
                  ic = jnumc(chlw,ic,i3)
                  icl = inumc(chlw,ic,i3,'i') - 1

                  call onum(chlw,ic,icl,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 998

                  irlip = nint( cvvv )

                  ic = inumc(chlw,icl,i3,'=') + 1
                  ic = jnumc(chlw,ic,i3)

                  call onum(chlw,ic,i3,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 998

                  idata = nint( cvvv )

            else
                  ic = i1
               do i = 1, 2
                  call snum(chlw,ic,i3,ic2,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 996
                  remsh(i) = cvvv
                  ic = ic2
               end do

               if( mskip .eq. 0 ) then

                  ldata = ldata + 1
                  cxalpha(ndata0,mdata,1,ldata) = remsh(1)
                  cxalpha(ndata0,mdata,2,ldata) = remsh(2)

               end if

            end if

                  goto 154

 66         close(iotmp)

               if( mskip .eq. 0 ) then
                  kndata(ndata0,mdata) = irzap
                  kldata(ndata0,mdata) = irlip
                  kedata(ndata0,mdata) = ldata
                  kidata(ndata0,mdata) = idata
               end if

*-----------------------------------------------------------------------

            jdata(ndata0) = mdata

         end if

 416     continue

*-----------------------------------------------------------------------
*          read data for gamma (photon)
*-----------------------------------------------------------------------

                  jdcar(1:2) = element(iz)
                  write(jdcar(3:5),'(i3.3)') ia
                  mltfl = chfn(27)(1:ilfn(27))//jdcar//'-y-g.dat'

                  do kk = 1, ndatg0
                     if( iza .eq. indatg(kk) ) goto 417
                  end do

*-----------------------------------------------------------------------

               if( itnda(mm) .eq. 2 .or. itnda(mm) .eq. 3 ) then
                  inquire( file = mltfl, exist = exex )
                  if( exex .eqv. .true. ) igyield = 1
               end if

*-----------------------------------------------------------------------

         if( igyield .eq. 1 ) then

                  ndatg0 = ndatg0 + 1

                  mdatg = 0
                  ldatg = 0
                  idatg = 2
                  mskip = 1

*-----------------------------------------------------------------------

            iotmp = 15
            open(iotmp, file = mltfl, status = 'old' )

            if ( mstz(165) .ge. 1 )
     &           write(*,'("Read ndata from ",a60)')  mltfl

               jsn = 1
               jsi = iotmp
               dsin(1) = mltfl
               idsi(1) = ilfn(27)+11
               ill(1) = 0
               ilf(1) = 10000000
               jpn = 0

 155           continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 996
               if( jpn .eq. 3 ) goto 67
               if( iskip .ne. 0 ) goto 155

            if( chcm(i1:i1+4) .eq. '#elab' ) then

               goto 155

            else if( chcm(i1:i1+3) .eq. '#zap' ) then

               if( mskip .eq. 0 ) then
                  kndatg(ndatg0,mdatg) = irzap
                  kldatg(ndatg0,mdatg) = irlip
                  kedatg(ndatg0,mdatg) = ldatg
                  kidatg(ndatg0,mdatg) = idatg
               end if

               ic = inumc(chlw,i1+5,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)
               icl = inumc(chlw,ic,i3,'l') - 1

               call onum(chlw,ic,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 996

                  irzap = nint( cvvv )

                  mskip = 0
                  mdatg = mdatg + 1
                  ldatg = 0

                  ic = inumc(chlw,icl,i3,'=') + 1
                  ic = jnumc(chlw,ic,i3)
                  icl = inumc(chlw,ic,i3,'i') - 1

                  call onum(chlw,ic,icl,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 998

                  irlip = nint( cvvv )

                  ic = inumc(chlw,icl,i3,'=') + 1
                  ic = jnumc(chlw,ic,i3)

                  call onum(chlw,ic,i3,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 998

                  idatg = nint( cvvv )

            else
                  ic = i1
               do i = 1, 2
                  call snum(chlw,ic,i3,ic2,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 996
                  remsh(i) = cvvv
                  ic = ic2
               end do

               if( mskip .eq. 0 ) then

                  ldatg = ldatg + 1
                  cxgamma(ndatg0,mdatg,1,ldatg) = remsh(1)
                  cxgamma(ndatg0,mdatg,2,ldatg) = remsh(2)

               end if

            end if

                  goto 155

 67         close(iotmp)

               if( mskip .eq. 0 ) then
                  kndatg(ndatg0,mdatg) = irzap
                  kldatg(ndatg0,mdatg) = irlip
                  kedatg(ndatg0,mdatg) = ldatg
                  kidatg(ndatg0,mdatg) = idatg
               end if

*-----------------------------------------------------------------------

            jdatg(ndatg0) = mdatg

         end if

 417     continue

*-----------------------------------------------------------------------

         end do

  411    continue
  410 continue

*-----------------------------------------------------------------------

          icontp = 0
      if( icontp .eq. 1 ) then
             write(*,*) 'proton', ndatp, mmdatp
         if( ndatp .gt. 0 ) then
          do i = 1, ndatp
             write(*,*) i, jdatp(i), indatp(i)
             do j = 1, jdatp(i)
                   write(*,*) i,j,kedatp(i,j), kndatp(i,j), kldatp(i,j),
     &                        kidatp(i,j)
             end do
          end do
         end if

             write(*,*) 'neutron', ndatn, mmdatn
         if( ndatn .gt. 0 ) then
          do i = 1, ndatn
             write(*,*) i, jdatn(i), indatn(i)
             do j = 1, jdatn(i)
                   write(*,*) i,j,kedatn(i,j), kndatn(i,j), kldatn(i,j),
     &                        kidatn(i,j)
             end do
          end do
         end if

             write(*,*) 'deuteron', ndatd, mmdatd
         if( ndatd .gt. 0 ) then
          do i = 1, ndatd
             write(*,*) i, jdatd(i), indatd(i)
             do j = 1, jdatd(i)
                   write(*,*) i,j,kedatd(i,j), kndatd(i,j), kldatd(i,j),
     &                        kidatd(i,j)
             end do
          end do
         end if

             write(*,*) 'alpha', ndata, mmdata
         if( ndata .gt. 0 ) then
          do i = 1, ndata
             write(*,*) i, jdata(i), indata(i)
             do j = 1, jdata(i)
                   write(*,*) i,j,kedata(i,j), kndata(i,j), kldata(i,j),
     &                        kidata(i,j)
             end do
          end do
         end if

             write(*,*) 'gamma', ndatg, mmdatg
         if( ndatg .gt. 0 ) then
          do i = 1, ndatg
             write(*,*) i, jdatg(i), indatg(i)
             do j = 1, jdatg(i)
                   write(*,*) i,j,kedatg(i,j), kndatg(i,j), kldatg(i,j),
     &                        kidatg(i,j)
             end do
          end do
         end if

             stop 888

      end if

*-----------------------------------------------------------------------

      end if


*-----------------------------------------------------------------------
*     set tally multipliers
*-----------------------------------------------------------------------

         ifm = 0

*-----------------------------------------------------------------------

      do 300 m = 1, itnm

         if( itmlp(m) .le. 0 ) goto 300

*-----------------------------------------------------------------------
*        check of multiplier -200 to -299
*-----------------------------------------------------------------------

      do it = 1, itmst(m)
            nn = it * 2 - 1
      do k = 1, itmlp(m)
      do j = 1, itmln(m,k)
                  jjin = mltp(1+nn,itmli(m,k)/13+j)
                  jjfn = mltp(1+nn+1,itmli(m,k)/13+j)
         do jj = jjin, jjfn, 2
         if( nint(slib(jj)) .eq. 2 ) then
         if( nint(slib(jj-2)) .eq. 5 .and.
     &       nint(slib(jj-1)) .ge. -299 .and.
     &       nint(slib(jj-1)) .le. -200 ) then

            idtal = nint(slib(jj-1))
            write(idcar,'(i3)') abs(idtal) ! T.Sato 2018/01/31
            if(idtal.lt.0) then ! Basically negative
             mltfl = chfn(26)(1:ilfn(26))//'/m'//idcar//'.inp'
            else
             mltfl = chfn(26)(1:ilfn(26))//'/p'//idcar//'.inp'
            endif

         do kk = 1, imltp

            idmlp = idmlt(kk)

            if( idmlp .eq. idtal ) goto 51

         end do

*-----------------------------------------------------------------------

            inquire( file = mltfl, exist = exex )

               if( exex .eqv. .false. ) then

                  write(io,'(/" Error : input data file for",
     &            " multiplier does not exist."/
     &            " file name = ",200a1)')
     &            ( mltfl(i:i),i=1, ilfn(26)+9 )

                  ErrCha = ''
                  MsgID = 'L:10763/R:setpar/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/" Error : input data file for",
     &            " multiplier does not exist."/
     &            " file name = ",200a1)')
     &            ( mltfl(i:i),i=1, ilfn(26)+9 )

                  goto 999

               end if

*-----------------------------------------------------------------------

         iotmp = 15
         open(iotmp, file = mltfl, status = 'old' )

               jsn = 1
               jsi = iotmp
               dsin(1) = mltfl
               idsi(1) = ilfn(26)+9
               ill(1) = 0
               ilf(1) = 10000000
               jpn = 0

 140        continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 998
               if( jpn .eq. 3 ) goto 52
               if( iskip .ne. 0 ) goto 140

  240       continue

            if( chcm(i1:i1+11) .eq. '[multiplier]' ) then

               if( chcm(i1+12:i1+14) .eq. 'off' ) then

                  jpn = 2
                  goto 140

               end if

               call multip(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)

                  imdfl(imltp) = 1

                  if( ierr .ne. 0 ) goto 998
                  if( jpn .eq. 3 ) goto 52

                  goto 240

            end if

                  goto 140

   52    close(iotmp)

*-----------------------------------------------------------------------

   51    continue

         end if
         end if
         end do
      end do
      end do
      end do


*-----------------------------------------------------------------------
*        fmcard library : multiplier
*-----------------------------------------------------------------------

      do it = 1, itmst(m)

            nn = it * 2 - 1

         do k = 1, itmlp(m)

            do 50 l = 1, itmpn(m,k)

               ityp = itmpt(m,k,l,1)
               ktyp = itmpt(m,k,l,2)

               if( ityp .eq. 1 ) then
                  ipim = 9
               else if( ityp .eq. 2 ) then
                  ipim = 1
               else if( ityp .eq. 14 ) then
                  ipim = 2
               else
                  goto 50
               end if

               do j = 1, itmln(m,k)

                     jjin = mltp(1+nn,itmli(m,k)/13+j)
                     jjfn = mltp(1+nn+1,itmli(m,k)/13+j)

                  do jj = jjin, jjfn, 2

                        if( nint(slib(jj)) .eq. 5 ) then

                           imat = nint(slib(jj+1))

                        else if( nint(slib(jj)) .eq. 6 ) then

                           lmat = imat
                           kmtt = nint(slib(jj+1))

                           ifm = ifm + 1
                           ifmi(ifm,1) = ipim
                           ifmi(ifm,2) = lmat
                           ifmi(ifm,3) = kmtt

                        else if( nint(slib(jj)) .eq. 7 ) then

                           jmat = nint(slib(jj+1))

                        else if( nint(slib(jj)) .eq. 8 ) then

                           lmat = jmat
                           kmtt = 1

                           ifm = ifm + 1
                           ifmi(ifm,1) = ipim
                           ifmi(ifm,2) = lmat
                           ifmi(ifm,3) = kmtt

                        end if

                  end do

               end do

   50       continue

         end do

      end do

*-----------------------------------------------------------------------

  300 continue

*-----------------------------------------------------------------------
*     icntl = 14 for t-volume
*-----------------------------------------------------------------------

      if( icntl .eq. 14 ) then

            m = 0

         do mm = 1, itnm

            if( ital(mm) .eq. 21 ) m = mm

         end do

            if( m .eq. 0 ) then

               write(io,'(" **** icntl=14 but no [t-volume] ")')
               ErrCha = ''
               MsgID = 'L:10929/R:setpar/F:read00.f'
               call ErrWrite(MsgID, ErrCha)
               write(jo,'(" **** icntl=14 but no [t-volume] ")')

               goto 999

            end if

*-----------------------------------------------------------------------

            do j = 1, 6

               isort(j,1) = 0
               isort(j,2) = 0
               isort(j,3) = 0
               isort(j,4) = 0

            end do

*-----------------------------------------------------------------------
*       define source for icntl = 14 and [t-volume]
*-----------------------------------------------------------------------


      if( itmth(m) .eq. 0 ) then

        if( itstp(m) .eq. 1 ) then

              imsrc = 1
              totfact = rtvr0(m)**2 * 3.1415926535897932d0
              j = 1
              smlwt(j) = 1.d0
              swt0(j) = 1.d0
              sfactor(j) = 1.d0

              se0(j) = 100.0
              istyp(j) = 2
              inkf0(j) = 2112

              jstyp(j) = 9
              jstypori(j) = 9
              sx0(j) = rtvx0(m)
              sy0(j) = rtvy0(m)
              sz0(j) = rtvz0(m)
              sr1(j) = rtvr0(m) + eps
              sr2(j) = rtvr0(m) + eps
              sdir(j) = -2.0
              sphi(j) = -1000.0
              sdom(j) = -1000.0
              jttyp(j) = -10000

        else if( itstp(m) .eq. 2 ) then

              imsrc = 6
              totfact = ( rtvx1(m) - rtvx0(m) )
     &                * ( rtvy1(m) - rtvy0(m) )
     &                + ( rtvy1(m) - rtvy0(m) )
     &                * ( rtvz1(m) - rtvz0(m) )
     &                + ( rtvz1(m) - rtvz0(m) )
     &                * ( rtvx1(m) - rtvx0(m) )
              totfact = totfact * 2.d0 / 6.d0

            do j = 1, 6

              swt0(j) = 1.d0
              se0(j) = 100.0
              istyp(j) = 2
              inkf0(j) = 2112
              jstyp(j) = 2
              jstypori(j) = 2
              sdom(j) = -1000.0
              jttyp(j) = -10000
              sfactor(j) = 1.d0

            end do

              smlwt(1) = ( rtvy1(m) - rtvy0(m) )
     &                 * ( rtvz1(m) - rtvz0(m) )
              smlwt(2) = smlwt(1)
              smlwt(3) = ( rtvz1(m) - rtvz0(m) )
     &                 * ( rtvx1(m) - rtvx0(m) )
              smlwt(4) = smlwt(3)
              smlwt(5) = ( rtvx1(m) - rtvx0(m) )
     &                 * ( rtvy1(m) - rtvy0(m) )
              smlwt(6) = smlwt(5)

              sx0(1) = rtvx0(m) - eps
              sx1(1) = rtvx0(m) + eps
              sy0(1) = rtvy0(m) - eps
              sy1(1) = rtvy1(m) + eps
              sz0(1) = rtvz0(m) - eps
              sz1(1) = rtvz1(m) + eps
              sdir(1) = 0.d0
              sphi(1) = 0.d0

              sx0(2) = rtvx1(m) - eps
              sx1(2) = rtvx1(m) + eps
              sy0(2) = rtvy0(m) - eps
              sy1(2) = rtvy1(m) + eps
              sz0(2) = rtvz0(m) - eps
              sz1(2) = rtvz1(m) + eps
              sdir(2) = 0.d0
              sphi(2) = 180.d0

              sx0(3) = rtvx0(m) - eps
              sx1(3) = rtvx1(m) + eps
              sy0(3) = rtvy0(m) - eps
              sy1(3) = rtvy0(m) + eps
              sz0(3) = rtvz0(m) - eps
              sz1(3) = rtvz1(m) + eps
              sdir(3) = 0.d0
              sphi(3) = 90.d0

              sx0(4) = rtvx0(m) - eps
              sx1(4) = rtvx1(m) + eps
              sy0(4) = rtvy1(m) - eps
              sy1(4) = rtvy1(m) + eps
              sz0(4) = rtvz0(m) - eps
              sz1(4) = rtvz1(m) + eps
              sdir(4) = 0.d0
              sphi(4) = -90.d0

              sx0(5) = rtvx0(m) - eps
              sx1(5) = rtvx1(m) + eps
              sy0(5) = rtvy0(m) - eps
              sy1(5) = rtvy1(m) + eps
              sz0(5) = rtvz0(m) - eps
              sz1(5) = rtvz0(m) + eps
              sdir(5) = 1.d0
              sphi(5) = 0.d0

              sx0(6) = rtvx0(m) - eps
              sx1(6) = rtvx1(m) + eps
              sy0(6) = rtvy0(m) - eps
              sy1(6) = rtvy1(m) + eps
              sz0(6) = rtvz1(m) - eps
              sz1(6) = rtvz1(m) + eps
              sdir(6) = -1.d0
              sphi(6) = 0.d0

        end if

      else if( itmth(m) .ne. 0 ) then

             icnt14 = 1

         if( itstp(m) .eq. 1 ) then

              imsrc = 1
              totfact = 4.d0/3.d0* rtvr0(m)**3 * 3.1415926535897932d0
              j = 1
              smlwt(j) = 1.d0
              swt0(j) = 1.d0
              sfactor(j) = 1.d0

              se0(j) = 100.0
              istyp(j) = 2
              inkf0(j) = 2112

              jstyp(j) = 9
              jstypori(j) = 9
              sx0(j) = rtvx0(m)
              sy0(j) = rtvy0(m)
              sz0(j) = rtvz0(m)
              sr1(j) = 0.d0
              sr2(j) = rtvr0(m)
              sdir(j) = 1.d0
              sphi(j) = -1000.0
              sdom(j) = -1000.0
              jttyp(j) = -10000

         else if( itstp(m) .eq. 2 ) then

              imsrc = 1
              totfact = ( rtvx1(m) - rtvx0(m) )
     &                * ( rtvy1(m) - rtvy0(m) )
     &                * ( rtvz1(m) - rtvz0(m) )
              totfact = totfact

              j = 1
              swt0(j) = 1.d0
              se0(j) = 100.0
              istyp(j) = 2
              inkf0(j) = 2112
              jstyp(j) = 2
              jstypori(j) = 2
              jttyp(j) = -10000
              sfactor(j) = 1.d0
              smlwt(j) = 1.d0
              sdir(j) = 1.d0
              sphi(j) = -1000.0
              sdom(j) = -1000.0

              sx0(j) = rtvx0(m)
              sx1(j) = rtvx1(m)
              sy0(j) = rtvy0(m)
              sy1(j) = rtvy1(m)
              sz0(j) = rtvz0(m)
              sz1(j) = rtvz1(m)

         end if

      end if

      end if


*-----------------------------------------------------------------------
*     icntl = 15 for t-wwbg
*-----------------------------------------------------------------------

      if( icntl .eq. 15 ) then

            m = 0

         do mm = 1, itnm

            if( ital(mm) .eq. 22 ) m = mm

         end do

            if( m .eq. 0 ) then

               write(io,'(" **** icntl=15 but no [t-wwbg] ")')
               ErrCha = ''
               MsgID = 'L:11154/R:setpar/F:read00.f'
               call ErrWrite(MsgID, ErrCha)
               write(jo,'(" **** icntl=15 but no [t-wwbg] ")')

               goto 999

            end if

*-----------------------------------------------------------------------
*       define source for icntl = 15 and [t-wwbg]
*-----------------------------------------------------------------------

               icnt14 = 1

               imsrc = 1
               totfact = 1.d0

               j = 1
               swt0(j) = 1.d0
               se0(j) = 100.0
               istyp(j) = 2
               inkf0(j) = 2112
               jstyp(j) = 1
               jstypori(j) = 1
               jttyp(j) = -10000
               sfactor(j) = 1.d0
               smlwt(j) = 1.d0
               sdir(j) = 200.d0
               sphi(j) = -1000.0
               sdom(j) = -1000.0

               sx0(j) = 0.d0
               sy0(j) = 0.d0

               sz0(j) = -tzmpo(m,itnms(m))
               sz1(j) =  tzmpo(m,itnms(m))

               sr0(j) =  trmpo(m,itnms(m))
               sr1(j) = 0.d0      ! S.Abe 2020/08/05, change "trmpo(m,itnms(m))" to "0.d0"

               isort(j,1) = itrwwbg(1)
               isort(j,2) = itrwwbg(2)
               isort(j,3) = itrwwbg(3)
               isort(j,4) = itrwwbg(4)

            do k = 1, 13

               rsort(j,k) = rtrwwbg(k)

            end do

      end if


*-----------------------------------------------------------------------

      call moddas_deallocate_cha(chrg)
      call moddas_deallocate_cha(chtm)

      return

*-----------------------------------------------------------------------

  996 continue

         write(*,'(/" ***** Error from Yield Files in XS *****")')
         write(*,*) dsin(k_err)(1:idsi(k_err)),l_err,':'

         goto 997

*-----------------------------------------------------------------------

  998 continue

         write(*,'(/" ***** Error from Multiplier Files *****")')
         write(*,*) dsin(k_err)(1:idsi(k_err)),l_err,':'

         goto 997

*-----------------------------------------------------------------------

  997 continue

            icf = 200

         do 910 i = 200, 1, -1

            if( m_err(i:i) .ne. ' ' ) goto 911

  910    continue

  911       icf = i

         call ErrWrite(ErrID, ErrCha)
         write(*,'(" error = ",200A1/)') ( m_err(i:i), i=1,icf )

*-----------------------------------------------------------------------

  999 continue

         ierr = 1

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine setegs(io,jo,ierr)
*                                                                      *
*       write material data for EGS                                    *
*       modified by K.Niita on 2014/08/21                              *
*                                                                      *
************************************************************************
      use moddas_material

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      include 'include/egs5_h.f'

*-----------------------------------------------------------------------

      common /kmat1a/ mxmat, mxmat0, mxnel
      common /kmat1g/ kmat(kvlmax)
      common /kmat1j/ intum
      common /celdg/  rhog(kvlmax)

      common /ndemax/ dnmax(20)

*-----------------------------------------------------------------------

       integer med_p2e(kvlmax)
       common /egs5cmn7/med_p2e

      real*8 emaxend
      common /egs5cmn1/emaxend

! 2015/5/28 T.Sato, Datapath common
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn, ilfn

*-----------------------------------------------------------------------

         ierr = 0

*-----------------------------------------------------------------------
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

               lemmx = 0

            do m = 1, mxmat
               lem = nint( dnel_das(kmat0+m) )
               if( lem .gt. lemmx ) lemmx = lem
            end do

*-----------------------------------------------------------------------
* Obtain the maximum calculation energy, emaxend,

         emaxend = max( dnmax(12), dnmax(14) )

*-----------------------------------------------------------------------

         open(566, file=chfn(23)(1:ilfn(23))//'.tmp', status='unknown')

               nmedium = 0

         do m = 1, mxmat

               nmedium = nmedium + 1  ! write only use to EGS5

               med_p2e(m) = nmedium   ! PHITS-mat to EGS-med list

               dnel = dnel_das(kmat0+m)
               denh = denh_das(kmat0+m)
               nel1  = nint( dnel )

               if(denh .ne. 0) then          ! if H is includ
                                             ! then add it
                  denst = denh
                  na = 1
                  nz = 1
                     write(566,679) nmedium, -rhog(m), 1, nz, na,
     $                              denst
               endif

               do l = 1, nel1

                  dicha = zz_das(kmat(m)+l)
                  dmasi = a_das(kmat(m)+l)
                  denst = den_das(kmat(m)+l)

                  nz = nint( dicha )
                  na = nint( dmasi )

                  if(denh .ne. 0) then
                     write(566,679) nmedium, -rhog(m), l+1, nz, na,
     $                              denst
                  else
                     write(566,679) nmedium, -rhog(m), l,   nz, na,
     $                              denst
                  endif

               enddo

         enddo

         close(566)

 679     format(i4,es15.7,3i5,es15.7)
c> ada.check
c  "lemmx"   "MXEPERMED"
c  "nmedium" is may be same material "mxmat" = MXMED
c  "l=nel1" is may be same number of  "mxnel" = MXNEL or MXEPERMED
c< ada.check
*-----------------------------------------------------------------------
*EGS5 end
*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine setpag(io,jo,ierr)
*                                                                      *
*       set parameters after geometry set up                           *
*       modified by K.Niita on 2011/02/03                              *
*                                                                      *
************************************************************************

      use fragdatamod
      use mod_counter, only:iTYinevt,iTYiregn,iTYmxmat !FURUTA20200119
      use moddas
      use moddas_character
      use moddas_fragdata
      use moddas_region
      use moddas_variance_reduction
      use tetramod, only: nelemtot
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /paraj/  mstz(300), parz(300)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /inpec/  ititl, ipara, ibody, iregn, llarr, itby, itar
      common /kmat1a/ mxmat, mxmat0, mxnel

      common /impreg/ dimp(kvlmax)
      common /regdd/  ivolm, iimpo
      common /regim/  iimps

      common /impmsg/ iimpn, isimp, iswct, isstr, maxip,
     &                mnimp(7,0:20),
     &                kfimp(7), inimc(7), inimt(7)

      common /fclmsg/ ifcln, isfcl, maxrg,
     &                mnfcl(6,0:20), mrfcl(kvlmax),
     &                kfcls(6), inflc(6), inflt(6)

      common /rclmsg/ ircln, isrcl, maxrr,
     &                mnrcl(6,0:20), mrrcl(kvlmax),
     &                krcls(6), lrcls(6), inrlc(6), inrlt(6),
     &                irpem(6,2), irman(6), irmct(6), jsmat(6)

      common /wwindp/ wupn, wsurvn, mxspln, mwhere, mvoww
      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww

*-----------------------------------------------------------------------

      common /wwind0/ dvals, iwmsh, kecho, kgwwp(6), idval
      common /wwxyz/ iwxty(6), iwxnm(6), iwxrg(6),
     &               rwxmi(6), rwxma(6), rwxdl(6),
     &               iwyty(6), iwynm(6), iwyrg(6),
     &               rwymi(6), rwyma(6), rwydl(6),
     &               iwzty(6), iwznm(6), iwzrg(6),
     &               rwzmi(6), rwzma(6), rwzdl(6)
      common /wwtrs/ iwwtr(6,4), rwwtr(6,13)
      common /wbxyz/ iwbxty(6), iwbxnm(6), iwbxrg(6),
     &               rwbxmi(6), rwbxma(6), rwbxdl(6),
     &               iwbyty(6), iwbynm(6), iwbyrg(6),
     &               rwbymi(6), rwbyma(6), rwbydl(6),
     &               iwbzty(6), iwbznm(6), iwbzrg(6),
     &               rwbzmi(6), rwbzma(6), rwbzdl(6)
      common /wbtrs/ iwbtr(6,4), rwbtr(6,13)

*-----------------------------------------------------------------------

      common /magreg/ nmreg, ingrc, ingrt, kmags, magtin, imgusr
      common /tmtreg/ ntmrg, intmc, intmt, ktime
      common /smireg/ nsreg, isgrc, isgrt, ksmir, ismir

      common /elmgrg/ nelctf, nereg, inerc, inert, kelcs

      common /splreg/ isptn, npreg(6), mnspt(6,0:20),
     &                ipgrc(6), ipgrt(6), ksplt(6), isplt(6),
     &                ispct(6,9), ispem(6,2)
      common /splrge/ espem(6,2)

      common /inggs/  iog, igcel, ioa, igsuf, iob, igtrs
      common /celdb/  idsn(kvlmax), idtn(kvlmax)

      common /isomul/ smlwt(isrc), totfact, imsrc
      common /isorsr/ nsmx(isrc), nsrn(isrc), nsrc(isrc)
      common /isorsc/ isort(isrc,4), rsort(isrc,13)

      common /wtcntl/ iwt, icimp(20), ifcls(20), iwwin(20), ircls(20)

      common /fgdata/ ifrgd, ifgm, ifgdf(kvlmax,5), frgfl(kvlmax)
      character frgfl*200

      common /ptname/ pname(20), ipln(20)
      character       pname*8

      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)

      common /regcrs/ icrsflx

*-----------------------------------------------------------------------
      common /wwbiasn/ eenwb(6,100), iwbdp, mnwbp(6,0:20), kfwbp(6),
     &                 inwbc(6), ienwb(6), inwbt(6), iswbp, maxwb

      common /nwwbias/ iwwbias

*-----------------------------------------------------------------------

      common /isosuf/ issuf(isrc), iscut(isrc), isvct(isrc,8),
     &                issfd(isrc), isvfd(isrc,8),
     &                ivsfd(isrc), ivvfd(isrc,8),
     &                dvsfd(isrc,4), dvvfd(isrc,8,4)

      common /issufd/ isxdf(isrc), isydf(isrc), iszdf(isrc),
     &                isdef(isrc), icdef(isrc),
     &                ixdef(isrc,5), iydef(isrc,5), izdef(isrc,5),
     &                vxpos(isrc,5), vypos(isrc,5), vzpos(isrc,5),
     &                sxpos(isrc), sypos(isrc), szpos(isrc), ssrad(isrc)

      dimension bval(50)
      dimension iscte(8)

      common /isorst/ jstyp(isrc), istyp(isrc), inkf0(isrc), lstyp(isrc)
!$OMP THREADPRIVATE(/isorst/)

*-----------------------------------------------------------------------

      dimension     idas(1)
      equivalence ( das, idas )


      character dmm*1
      character chbd*3

*-----------------------------------------------------------------------

      dimension ipar(0:20)
      dimension jpar(0:20)
      dimension klev(0:20)

      logical   exex

*-----------------------------------------------------------------------
cKN 2018/10/29 jstyp->kstyp, jnkf0->knkf0, istyp->mstyp, inkf0->lnkf0

      integer kstyp(6), knkf0(6)
      character chin*200, chlw*200, chcm*200
      integer ic,nch,mstyp,lnkf0,ierrtmp,ichg,imas

C S.Hashimoto revised for the new [Frag Data]. (2014.12.26)
      integer   jsn,jpn,iskip
      character dsin(0:9)*200
      integer   idsi(0:9), ill(0:9), ilf(0:9)
C S.Hashimoto added a new flag for [Frag Data]. (2015.8.25)
      common /initfgdata/ init_frag ! S.H. corrected. (2016.9.25)
      integer   init_frag

*-----------------------------------------------------------------------
      common /mtnreg/ dmhsg(-1:kvlmax),
     &                nmtng(-1:kvlmax), dmtng(-1:kvlmax)
      character dmtng*80

      common /mtreg/  smtrg(kvlmax), dmtrg(kvlmax),
     &                mtrgn, mtrg(kvlmax,2), nmtrg(kvlmax)
      character dmtrg*80

      character mtm*6

*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr
      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

      common /tall40/ rtorg(itlmax,3), rteye(itlmax,3), rtlit(itlmax,3),
     &                rtwin(itlmax,3), rtbrt(itlmax,2), ithvn(itlmax),
     &                itwin(itlmax,2), itbox(itlmax), itmir(itlmax),
     &                rtbox(itlmax,5,10), rtout(itlmax), rthet(itlmax),
     &                itlin(itlmax), itshd(itlmax), itgxs(itlmax)

      common /tall41/ rdmax

      common /ccggg/  icgg

      common /regdm/  idmg(kvlmax)
      common /regda/  ichl(kvlmax), chsm(kvlmax), ichmx, iod
      character       chsm*10
      common /regde/  ichp(kvlmax), ilat(kvlmax), idct(kvlmax)
      common /celda/  deng(kvlmax)

      common /rdouti/ nrdout, nrdi(100)
      data nrdout / 0 /

      integer italsh
      common /talsh/ italsh

      integer,allocatable :: jdata(:),idas_temporary(:)

*-----------------------------------------------------------------------

         ierr = 0

         jpn = 0

*-----------------------------------------------------------------------
*     icntl = 11, 3dshow,  add big sphere for outer void
*-----------------------------------------------------------------------

      if( icntl .eq. 11 .or. icntl .eq. 14 .or. icntl .eq. 15 ) then

               rdmax = 1.0

            do i = 1, itnm

               if( ital(i) .eq. 11 .or.
     &             ital(i) .eq. 21 .or. ital(i) .eq. 22 ) then

                  if( rtout(i) .gt. rdmax ) rdmax = rtout(i)

               end if

            end do

         if( rdmax .gt. 1.0 ) then

*-----------------------------------------------------------------------
*        for GG
*-----------------------------------------------------------------------

            if( icgg .eq. 1 ) then

               rewind  ioa

               do i = 1, igsuf

                  read(ioa) idrf, idtr, idsf, igkst,
     &                    ( bval(j), j = 1, igkst )

               end do

                     idrf    = 0
                     idtr    = 0
                     idsf    = 5
                     igkst   = 1
                     bval(1) = rdmax

                  write(ioa) idrf, idtr, idsf, igkst,
     &                     ( bval(j), j = 1, igkst )

                  igsuf = igsuf + 1
                  idsn(igsuf) = kvmmax

                  ichmx = ichmx + 10
                  mcmx = ( mdas / 2 - 1 ) * 8 + 1
       if(mcmx.lt.0.0) then ! T.Sato 2020/09/17
         write(*,'("mdas is too large, you have to activate ",
     &   "integer*8 option in param.inc")')
         goto 999
       endif

               call moddas_allocate_cha(MAX_NUM_CHRG, chrg)
                  mci = 0

                  rewind iod

                  ioe = 22

                  open(ioe,status='scratch',form='unformatted')

               do i = 1, igcel

                     read(iod) (chrg(mci+k:mci+k),k=1,ichp(i))

                  if( idmg(i) .eq. -1 ) then

                     nrdout = nrdout + 1
                     nrdi(nrdout) = idrg(i)

                     idmg(i) = 0

                     do k = ichp(i), ichl(i) + 1, -1

                        chrg(mci+k+14:mci+k+14) = chrg(mci+k:mci+k)

                     end do

                        chrg(mci+ichl(i)+3:mci+ichl(i)+4) = ') '

                     do k = ichl(i), 1, -1

                        chrg(mci+k+2:mci+k+2) = chrg(mci+k:mci+k)

                     end do

                        chrg(mci+1:mci+2) = '( '

                        write(chrg(mci+ichl(i)+5:
     &                             mci+ichl(i)+14),'(i10)') -kvmmax

                        ichl(i) = ichl(i) + 14
                        ichp(i) = ichp(i) + 14


                  end if

                     write(ioe) (chrg(mci+k:mci+k),k=1,ichp(i))

               end do

                     igcel = igcel + 1
                     iregn = igcel
                     if(igcel.gt.kvlmax) goto 901
                     ichl(igcel)  = 10
                     ichp(igcel)  = 10
                     idmg(igcel)  = -1
                     idrg(igcel)  = kvmmax
                     idgr(kvmmax) = igcel

                     write(chrg(mci+1:mci+10),'(i10)') kvmmax
                     write(ioe) (chrg(mci+k:mci+k),k=1,ichp(igcel))

                  rewind iod
                  rewind ioe

               do i = 1, igcel

                     read(ioe) (chrg(mci+k:mci+k),k=1,ichp(i))
                    write(iod) (chrg(mci+k:mci+k),k=1,ichp(i))

               end do

                  close(ioe)

*-----------------------------------------------------------------------
*        for CG
*-----------------------------------------------------------------------

            else

               rewind itby

               read(itby) chtit
               read(itby) ipva

               do k = 1, ibody

                  read(itby) chbd, ibnum, ibva, ( bval(i), i = 1, ibva )

               end do

                  ibody = ibody + 1
                  chbd  = 'sph'
                  ibnum = 9999
                  ibva  = 4
                  bval(1) = 0.0d0
                  bval(2) = 0.0d0
                  bval(3) = 0.0d0
                  bval(4) = rdmax

                  write(itby) chbd, ibnum, ibva,
     &                      ( bval(i), i = 1, ibva )

                  ichmx = ichmx + 10
                  mcmx = ( mdas / 2 - 1 ) * 8 + 1
       if(mcmx.lt.0.0) then ! T.Sato 2020/09/17
         write(*,'("mdas is too large, you have to activate ",
     &   "integer*8 option in param.inc")')
         goto 999
       endif

                  rewind iod

                  ioe = 22
                  open(ioe,status='scratch',form='unformatted')

               do i = 1, iregn

                     read(iod) (chrg(mci+k:mci+k),k=1,ichl(i))

                  if( idmg(i) .eq. -1 ) then

                     idmg(i) = 0

                        write(chrg(mci+ichl(i)+1:
     &                             mci+ichl(i)+10),'(i10)') 9999

                        ichl(i) = ichl(i) + 10

                  end if

                     write(ioe) (chrg(mci+k:mci+k),k=1,ichl(i))

               end do

                     iregn = iregn + 1
                     igcel = iregn

                     ichl(iregn)  = 10
                     idmg(iregn)  = -1
                     idrg(iregn)  = 9999
                     idgr(9999) = iregn
                     chsm(iregn)  = 'tmp'
                     deng(iregn)  = 0.0

                     write(chrg(mci+1:mci+10),'(i10)') -9999
                     write(ioe) (chrg(mci+k:mci+k),k=1,ichl(iregn))

                  rewind iod
                  rewind ioe

               do i = 1, igcel

                     read(ioe) (chrg(mci+k:mci+k),k=1,ichl(i))
                    write(iod) (chrg(mci+k:mci+k),k=1,ichl(i))

               end do

                  close(ioe)

            end if

*-----------------------------------------------------------------------

         end if

      end if

*-----------------------------------------------------------------------
*        summary of reg name
*-----------------------------------------------------------------------
                  nmtng(-1) = 10
                  dmtng(-1)(1:10) = 'outer void'
                  dmhsg(-1) = 1.0d0

                  nmtng(0) = 4
                  dmtng(0)(1:4) = 'void'
                  dmhsg(0) = 1.0d0

                  iclr = 0

            do i = 1, igcel

                  write(mtm,'(i6)') idrg(i)

               do k = 1, 6

                  if( mtm(k:k) .ne. ' ' ) goto 131

               end do

                  k = 1

  131             inm = k

                  nmtng(i) = 6 - k + 1
                  dmtng(i)(1:nmtng(i)) = mtm(inm:6)

                  dmhsg(i) = 1.0d0

            end do

*-----------------------------------------------------------------------

         if( mtrgn .gt. 0 ) then

            do i = 1, mtrgn

               do k = mtrg(i,1), mtrg(i,2)

                  if( k .gt. 0 ) then

                     j = idgr(k)

                     if( j .eq. 0 ) then

                        goto 444

                     else

                        if( nmtrg(i) .gt. 0 ) then
                           nmtng(j) = nmtrg(i)
                           dmtng(j)(1:nmtng(j)) = dmtrg(i)(1:nmtng(j))
                        end if

                     end if

                           dmhsg(j) = smtrg(i)

                  else if( k .eq. 0 ) then

                        if( nmtrg(i) .gt. 0 ) then
                           nmtng(0) = nmtrg(i)
                           dmtng(0)(1:nmtng(0)) = dmtrg(i)(1:nmtng(0))
                        end if

                           dmhsg(0) = smtrg(i)

                  end if

  444             continue

               end do

            end do

         end if

*-----------------------------------------------------------------------
*        check transform id in source
*-----------------------------------------------------------------------

            do j = 1, imsrc

               if( isort(j,1) .gt. 0 ) then

                     l = 0

                  do k = 1, igtrs

                     if( isort(j,3) .eq. idtn(k) ) l = k

                  end do

                  if( l .gt. 0 ) then

                     isort(j,4) = l

                  else

                     write(io,'("*Error: in [source],"
     &               " trcl =",i7," is not defined ",
     &               "in [transform]")') isort(j,3)

                     ErrCha = ''
                     MsgID = 'L:11926/R:setpag/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'("*Error: in [source],"
     &               " trcl =",i7," is not defined ",
     &               "in [transform]")') isort(j,3)

                     ierr = ierr + 1
                     goto 999

                  end if

               end if

            end do


*-----------------------------------------------------------------------
*     setting fragment XS data by reading frag data files
*-----------------------------------------------------------------------

      if( ifrgd .gt. 0 ) then

                  mmmin = mmmax

                  call moddas_allocate_int(ifrgd+1, ifgs01_nei)
                  call moddas_allocate_int(ifrgd+1, ifgs02_kne)
                  call moddas_allocate_int(ifrgd+1, ifgs03_kxs)
                  call moddas_allocate_int(ifrgd+1, ifgs04_neo)
                  call moddas_allocate_int(ifrgd+1, ifgs05_kef)
                  call moddas_allocate_int(ifrgd+1, ifgs06_nag)
                  call moddas_allocate_int(ifrgd+1, ifgs07_kaf)
                  call moddas_allocate_int(ifrgd+1, ifgs08_nfrg)
                  call moddas_allocate_int(ifrgd+1, ifgs09_kim)
                  call moddas_allocate_int(ifrgd+1, ifgs10_ks0)

*-----------------------------------------------------------------------
C do-loop for frag data files
         do i = 1, ifrgd

            if( ifgdf(i,1) .gt. 0 ) then

               inquire( file = frgfl(i), exist = exex )

               if( exex .eqv. .false. ) goto 300

               iot = 31
               open(iot, file = frgfl(i), status = 'old' )

               jsn  = 0
               ill(jsn) = 0
               ilf(jsn) = 10000000


C projectile
               iskip = 100
               do iii=1,100000000 ! optimization bug, do while sentence cannot be used here
                if(iskip.eq.0) exit
                iskip = 0
                call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &               jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
               end do
               if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 302
               call rdpname(i1,i3,chlw,mstyp,lnkf0,kstyp,knkf0,ierrtmp)
               if ( ierrtmp .ne. 0 ) goto 302

               if ( mstyp .eq. 11 ) then ! for old format: Z*1000+A

                if ( lnkf0 .eq. 1 ) then ! neutron
                   iproj = 2112
                else if ( lnkf0 .eq. 1001 ) then ! proton
                   iproj = 2212
                else if ( (lnkf0 .gt. 1001 .and. lnkf0 .lt. 2112)
     &                  .or. (lnkf0 .ge. 3001 .and. lnkf0 .lt. 3112)
     &                  .or. (lnkf0 .ge. 4001) ) then ! nucleus
                   ichg = abs( lnkf0 / 1000 )
                   imas = abs( lnkf0 - lnkf0 / 1000 * 1000 )
                   if ( ichg .ge. 1000 .or. imas .ge. 1000) goto 303
                   iproj = ichg * 1000000 + imas
                else ! for kf-code when ityp=11 (other particles)
                   iproj = lnkf0
                end if

               else ! for kf-code
                iproj = lnkf0
               end if

               if( iproj .ne. ifgdf(i,2) ) goto 304


C target
               iskip = 100
               do while ( iskip .ne. 0 )
                iskip = 0
                call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &               jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
               end do
               if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 305
               call rdpname(i1,i3,chlw,mstyp,lnkf0,kstyp,knkf0,ierrtmp)
               if ( ierrtmp .ne. 0 ) goto 305

               if ( mstyp .eq. 11 ) then ! for old format: Z*1000+A

                if ( lnkf0 .eq. 1001 ) then ! proton
                   itarg = 2212
                else if ( (lnkf0 .gt. 1001 .and. lnkf0 .lt. 2112)
     &                  .or. (lnkf0 .ge. 3001 .and. lnkf0 .lt. 3112)
     &                  .or. (lnkf0 .ge. 4001) ) then ! nucleus
                   ichg = abs( lnkf0 / 1000 )
                   imas = abs( lnkf0 - lnkf0 / 1000 * 1000 )
                   if ( ichg .ge. 1000 .or. imas .ge. 1000) goto 306
                   itarg = ichg * 1000000 + imas
                else ! for kf-code when ityp=11 (other particles)
                   iproj = lnkf0
                end if

               else ! for kf-code
                itarg = lnkf0
               end if

               if( itarg .ne. ifgdf(i,3) ) goto 307


C nei: number of incident energy bin
               iskip = 100
               do while ( iskip .ne. 0 )
                iskip = 0
                call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &               jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
               end do
               if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 308
               call onum(chlw,i1,i3,prn,ierrtmp)
               if ( ierrtmp .ne. 0 ) goto 308
               nei = idnint( prn )
               kne = ifgs02_kne(i)
               ifgs01_nei(i) = nei
               call moddas_reallocate_dbl(
     &                 ifrgd+1, i, nei+1, ifgs02_kne, frgne)


C ein(nei+1): data points of incident energy (MeV/n)
               j = 0
               do while ( j .lt. nei+1 )

                call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
                if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 309
                ic = i1

                do while ( ic .le. i3 )
                 call snum(chlw,ic,i3,ic2,prn,ierrtmp)
                 if ( ierrtmp .ne. 0 ) goto 309
                 ic = ic2
                 j = j+1
                 frgne(kne+j) = prn
                end do

               end do

               if ( j .ne. nei+1 ) goto 310

               kxs = ifgs03_kxs(i)
               call moddas_reallocate_dbl(
     &                 ifrgd+1, i, nei+1, ifgs03_kxs, frgxs)


C sigtot(nei+1): total reaction cross sections (mb)
               j = 0
               do while ( j .lt. nei+1 )

                call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
                if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 311
                ic = i1

                do while ( ic .le. i3 )
                 call snum(chlw,ic,i3,ic2,prn,ierrtmp)
                 if ( ierrtmp .ne. 0 ) goto 311
                 ic = ic2
                 j = j+1
                 frgxs(kxs+j) = prn
                end do

               end do

               if ( j .ne. nei+1 ) goto 312

C When sigtot is not given (= 0mb), total cross section model is used.
               do inei = 1, nei+1
                if ( frgxs(kxs+inei) .le. 0d0 ) then
                 ein = frgne(kne+inei)
                 call txsmodel(iproj,itarg,ein,signe)
                 frgxs(kxs+inei) = signe * 1d3 ! b -> mb
                end if
               end do


C neo: number of XS data for energy of fragment particles
C nef: number of energy mesh point of fragment particles (= neo+1 or neo)
               iskip = 100
               do while ( iskip .ne. 0 )
                iskip = 0
                call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &               jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
               end do
               if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 313

               if (chlw(i1:i1+4).eq.'model') then
                  ifgdf(i,5) = 1
                  exit
               end if

               call onum(chlw,i1,i3,prn,ierrtmp)
               if ( ierrtmp .ne. 0 ) goto 313
               neo = idnint( prn )
               ifgs04_neo(i) = neo

               if ( ifgdf(i,1).eq.4 .and. neo.le.0 ) goto 314

               if ( neo .gt. 0 ) then
                  if ( ifgdf(i,1) .eq. 5 ) then
                     nef = neo
                  else
                     nef = neo + 1
                  end if
               else
                  nef = neo
               end if


C frgef(nef): energy mesh points of fragment particles (MeV)
               if ( nef .gt. 0 ) then

                kef = ifgs05_kef(i)
                call moddas_reallocate_dbl(
     &                  ifrgd+1, i, nef+1, ifgs05_kef, frgef)

                j = 0
                do while ( j .lt. nef )

                 call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
                 if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 315
                 ic = i1

                 do while ( ic .le. i3 )
                  call snum(chlw,ic,i3,ic2,prn,ierrtmp)
                  if ( ierrtmp .ne. 0 ) goto 315
                  ic = ic2
                  j = j+1
                  frgef(kef+j) = prn
                 end do

                end do

                if ( j .ne. nef ) then
                  if ( ifgdf(i,1) .eq. 5 ) then
                     goto 317
                  else
                     goto 316
                  end if
                end if

                if ( ifgdf(i,1).eq.4 .and. frgef(kef+1).gt.0d0 ) then
                 do j = nef, 1, -1
                  frgef(kef+j+1) = frgef(kef+j)
                 end do
                 frgef(kef+1) = 0d0
                 neo = nef
                 ifgs04_neo(i) = neo
                 iextene = 1
                else
                 iextene = 0
                end if

C nef < 0: for discrete cross sections
               else if( nef .lt. 0 ) then

                kef = ifgs05_kef(i)
                call moddas_reallocate_dbl(
     &                  ifrgd+1, i, iabs(nef), ifgs05_kef, frgef)

                j = 0
                do while ( j .lt. iabs(nef) )

                 call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
                 if ( ierrtmp .ne. 0 ) goto 315
                 ic = i1

                 do while ( ic .le. i3 )
                  call snum(chlw,ic,i3,ic2,prn,ierrtmp)
                  if ( ierrtmp .ne. 0 ) goto 315
                  ic = ic2
                  j = j+1
                  frgef(kef+j) = prn
                 end do

                end do

                if ( j .ne. iabs(nef) ) goto 317

               end if


C nag: number of XS data for angle of fragment particles
C naf: number of angle mesh point of fragment particles (= |nag|+1 or |nag|)
               iskip = 100
               do while ( iskip .ne. 0 )
                iskip = 0
                call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &               jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
               end do
               if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 318
               call onum(chlw,i1,i3,prn,ierrtmp)
               if ( ierrtmp .ne. 0 ) goto 318
               nag = idnint( prn )
               ifgs06_nag(i) = nag

               if ( ifgdf(i,1).eq.4 .and. nag.eq.0 ) goto 314

               if ( ifgdf(i,1) .eq. 5 ) then
                  naf = nag
               else
                  if ( nag .gt. 0 ) then
                     naf = nag + 1
                  else if ( nag .lt. 0 ) then
                     naf = nag - 1
                  else
                     naf = nag
                  end if
               end if

C frgaf(naf): angular mesh points (rad(when nag>0), degree(when nag<0))
               if( iabs(naf) .gt. 0 ) then

                kaf = ifgs07_kaf(i)
                call moddas_reallocate_dbl(
     &                  ifrgd+1, i, iabs(naf)+2, ifgs07_kaf, frgaf)

                j = 0
                do while ( j .lt. iabs(naf) )

                 call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
                 if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 319
                 ic = i1

                 do while ( ic .le. i3 )
                  call snum(chlw,ic,i3,ic2,prn,ierrtmp)
                  if ( ierrtmp .ne. 0 ) goto 319
                  ic = ic2
                  j = j+1
                  if ( naf .gt. 0 ) then
                   frgaf(kaf+j) = prn
                  else if ( naf .lt. 0 ) then
                   frgaf(kaf+j) = prn * pi / 180d0
                  end if
                 end do

                end do

                if ( j .ne. iabs(naf) ) then
                  if ( ifgdf(i,1) .eq. 5 ) then
                     goto 341
                  else
                     goto 320
                  end if
                end if

                if ( ifgdf(i,1).eq.4 .and. frgaf(kaf+1).gt.0d0 ) then
                 do j = iabs(naf), 1, -1
                  frgaf(kaf+j+1) = frgaf(kaf+j)
                 end do
                 frgaf(kaf+1) = 0d0
                 if ( naf .gt. 0 ) then
                  naf = naf + 1
                  nag = naf
                 else if ( naf .lt. 0 ) then
                  naf = naf - 1
                  nag = naf
                 end if
                 ifgs06_nag(i) = nag
                 iexta0 = 1
                else
                 iexta0 = 0
                end if

                if ( ifgdf(i,1).eq.4
     &               .and. frgaf(kaf+iabs(naf)).lt.pi ) then
                 frgaf(kaf+iabs(naf)+1) = pi
                 nag = naf
                 ifgs06_nag(i) = nag
                 iextapi = 1
                else
                 iextapi = 0
                end if

               end if


C nfrg: number of fragment particles
               iskip = 100
               do while ( iskip .ne. 0 )
                iskip = 0
                call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &               jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
               end do
               if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 321
               call onum(chlw,i1,i3,prn,ierrtmp)
               if ( ierrtmp .ne. 0 ) goto 321
               nfrg = idnint( prn )
               ifgs08_nfrg(i) = nfrg
               kim = ifgs09_kim(i)
               call moddas_reallocate_int(
     &                 ifrgd+1, i, nfrg, ifgs09_kim, ifrgm)


C ifrgm(nfrg): kind of fragment particles
               j = 0
               do while ( j .lt. nfrg )

                call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
                if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 322
                ic = i1

                do while ( ic .le. i3 )
                 call rdpname(ic,i3,chlw,mstyp,lnkf0,kstyp,knkf0,
     &                  ierrtmp)
                 if ( ierrtmp .ne. 0 ) goto 322
                 j = j+1
                 ifrgm(kim+j) = lnkf0
                end do

               end do

               if ( j .ne. nfrg ) goto 323

               do j=1, nfrg

C for old format: Z*1000+A
                if ( ifrgm(kim+j) .eq. 1 ) then ! neutron
                 ifrgm(kim+j) = 2112

                else if ( ifrgm(kim+j) .eq. 1001 ) then ! proton
                 ifrgm(kim+j) = 2212

                else if ( (ifrgm(kim+j) .gt. 1001
     &                .and. ifrgm(kim+j) .lt. 2112)
     &                .or. (ifrgm(kim+j) .ge. 3001
     &                .and. ifrgm(kim+j) .lt. 3112)
     &                .or. (ifrgm(kim+j) .ge. 4001
     &                .and. ifrgm(kim+j) .lt. 1000001) ) then ! nucleus
                 ichg = abs( ifrgm(kim+j) / 1000 )
                 imas = abs( ifrgm(kim+j) - ifrgm(kim+j)/1000*1000 )
                 if ( ichg .ge. 1000 .or. imas .ge. 1000) goto 324
                 ifrgm(kim+j) = ichg*1000000 + imas

                end if

               end do

*-----------------------------------------------------------------------

               ks1 = 0
               ks2 = 0
               ks3 = 0

               ks0 = ifgs10_ks0(i)
               call moddas_reallocate_int(
     &                 ifrgd+1, i, nei+2, ifgs10_ks0, ifrge_ksf)
               call moddas_reallocate_int(
     &                 ifrgd+1, i, nei+2, ifgs10_ks0, ifrge_ks1)
               call moddas_reallocate_int(
     &                 ifrgd+1, i, nei+2, ifgs10_ks0, ifrge_ks2)
               call moddas_reallocate_int(
     &                 ifrgd+1, i, nei+2, ifgs10_ks0, ifrge_ks3)
               num_frgsf = nfrg
               num_frgee = 0
               num_frgdd = 0
               num_frgdx = 0

C case 1
               if ( nef .gt. 0 .and. naf .ne. 0 ) then

                num_frgdd = nfrg*nef*iabs(naf)

C case 2
               else if ( nef .eq. 0 .and. naf .ne. 0 ) then

                num_frgee = nfrg*2
                num_frgdx = nfrg*iabs(naf)

C case 3
               else if ( nef .eq. 0 .and. naf .eq. 0 ) then

                num_frgee = nfrg*2

C case 4
               else if ( nef .gt. 0 .and. naf .eq. 0 ) then

                num_frgdx = nfrg*nef

C case 5
               else if ( nef .lt. 0 .and. naf .ne. 0 ) then

                num_frgdd = nfrg*iabs(nef)*iabs(naf)

C case 6
               else if ( nef .lt. 0 .and. naf .eq. 0 ) then

                num_frgdx = nfrg*iabs(nef)

               end if

               if ( ifgdf(i,1).eq.4 .and. iexta0 .eq. 1 ) then
                if ( naf .gt. 0 ) then
                 naf = naf - 1
                else if ( naf .lt. 0 ) then
                 naf = naf + 1
                end if
               end if

*-----------------------------------------------------------------------
C do-loop for incident energy mesh points
            if( i > 1 ) then
               ifrge_ksf(ks0+1) = ifrge_ksf(ks0)
               ifrge_ks1(ks0+1) = ifrge_ks1(ks0)
               ifrge_ks2(ks0+1) = ifrge_ks2(ks0)
               ifrge_ks3(ks0+1) = ifrge_ks3(ks0)
            end if
            do m = 1, nei + 1

             mmlm = ifrge_ksf(ks0+m)
             call moddas_reallocate_dbl(
     &               ifgs10_ks0(i+1), ks0+m, nfrg, ifrge_ksf, frgsf)

C frgsf((nei+1)*nfrg): production cross sections of fragment particles (mb)
             j = 0
             do while ( j .lt. nfrg )
              iskip = 100
              do while ( iskip .ne. 0 )
               iskip = 0
               call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &              jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
              end do
              if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 325
              ic = i1

              do while ( ic .le. i3 )
               call snum(chlw,ic,i3,ic2,prn,ierrtmp)
               if ( ierrtmp .ne. 0 ) goto 325
               ic = ic2
               j = j+1
               frgsf(mmlm+j) = prn
              end do

             end do

             if ( j .ne. nfrg ) goto 326

             mmlm = 0

             ks1 = ifrge_ks1(ks0+m)
             if( num_frgee > 0 ) then
                call moddas_reallocate_dbl(
     &                  ifgs10_ks0(i+1), ks0+m,
     &                  num_frgee, ifrge_ks1, frgee)
             end if

             ks2 = ifrge_ks2(ks0+m)
             if( num_frgdd > 0 ) then
                call moddas_reallocate_dbl(
     &                  ifgs10_ks0(i+1), ks0+m,
     &                  num_frgdd, ifrge_ks2, frgdd)
             end if

             ks3 = ifrge_ks3(ks0+m)
             if( num_frgdx > 0 ) then
                call moddas_reallocate_dbl(
     &                  ifgs10_ks0(i+1), ks0+m,
     &                  num_frgdx, ifrge_ks3, frgdx)
             end if

C case 1
C frgdd((nei+1)*nfrg*neo*nag): double differential cross sections (mb/MeV/sr)
             if ( neo .gt. 0 .and. nag .ne. 0 ) then

              nfrgdata = nfrg*neo*iabs(nag)
              if ( ifgdf(i,1) .eq. 4 ) then

               if ( iextene .eq. 1 ) then
                neotmp = neo - 1
               else
                neotmp = neo
               end if
               if ( iexta0 .eq. 1 ) then
                nagtmp = iabs(nag) - 1
               else
                nagtmp = iabs(nag)
               end if
               if ( iextapi .eq. 1 ) then
                nagtmp = nagtmp - 1
               else
                nagtmp = nagtmp
               end if

               nfrgdata = nfrg*neotmp*nagtmp

              end if

              j = 0
              do while ( j .lt. nfrgdata )
               iskip = 100
               do while ( iskip .ne. 0 )
                iskip = 0
                call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &               jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
               end do
               if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 327
               ic = i1

               do while ( ic .le. i3 )
                call snum(chlw,ic,i3,ic2,prn,ierrtmp)
                if ( ierrtmp .ne. 0 ) goto 327
                ic = ic2
                if ( j.eq.0 .and. prn.lt.0d0 ) goto 328
                if ( prn .ge. 0d0 ) then
                 j = j+1
                 frgdd(mmlm+ks2+j) = prn
                else
                 n0repeat = idnint(dabs(prn))
                 do i0repeat = 1, n0repeat
                  j = j+1
                  frgdd(mmlm+ks2+j) = frgdd(mmlm+ks2+j-1)
                 end do
                end if
               end do

              end do

              if ( j .ne. nfrgdata ) goto 329


              if ( ifgdf(i,1) .eq. 4 ) then

               if ( iextene .eq. 1 ) then

                nfrgdata = nfrgdata+nfrg*nagtmp
                infrgdata = 0
                do infrg = nfrg, 1, -1
                 do ineo = neotmp, 1, -1
                  do inag = nagtmp, 1, -1
                   infrgdata = infrgdata + 1
                   frgdd(mmlm+ks2+nfrgdata-infrgdata+1)
     &               = frgdd(mmlm+ks2+nfrgdata-infrg*nagtmp-infrgdata+1)
                  end do
                 end do
                 do inag = nagtmp, 1, -1
                  infrgdata = infrgdata + 1
                  frgdd(mmlm+ks2+nfrgdata-infrgdata+1)
     &                 = frgdd(mmlm+ks2+nfrgdata+nagtmp-infrgdata+1)
                 end do
                end do

               end if

               if ( iexta0 .eq. 1 ) then

                nfrgdata = nfrgdata+nfrg*neo
                infrgdata = 0
                do infrg = nfrg, 1, -1
                 do ineo = neo, 1, -1
                  do inag = nagtmp, 1, -1
                   infrgdata = infrgdata + 1
                   frgdd(mmlm+ks2+nfrgdata-infrgdata+1)
     &                  = frgdd(mmlm+ks2+nfrgdata
     &                  +(1-infrg)*neo-ineo-infrgdata+1)
                  end do
                  infrgdata = infrgdata + 1
                  frgdd(mmlm+ks2+nfrgdata-infrgdata+1)
     &                 = frgdd(mmlm+ks2+nfrgdata-infrgdata+1+1)
                 end do
                end do

                nagtmp = nagtmp +1

               end if

               if ( iextapi .eq. 1 ) then

                nfrgdata = nfrgdata+nfrg*neo
                infrgdata = 0
                do infrg = nfrg, 1, -1
                 do ineo = neo, 1, -1
                  infrgdata = infrgdata + 1
                  frgdd(mmlm+ks2+nfrgdata-infrgdata+1)
     &                 = frgdd(mmlm+ks2+nfrgdata
     &                 +(1-infrg)*neo-ineo-infrgdata+1)
                  do inag = nagtmp, 1, -1
                   infrgdata = infrgdata + 1
                   frgdd(mmlm+ks2+nfrgdata-infrgdata+1)
     &                  = frgdd(mmlm+ks2+nfrgdata
     &                  +(1-infrg)*neo-ineo-infrgdata+1+1)
                  end do
                 end do
                end do

               end if

              end if


C case 2
C frgee((nei+1)*nfrg*2): mean value and dispersion for Gaussian (MeV)
             else if ( neo .eq. 0 .and. nag .ne. 0 ) then

              j = 0
              do while ( j .lt. nfrg*2 )
               iskip = 100
               do while ( iskip .ne. 0 )
                iskip = 0
                call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &               jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
               end do
               if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 330
               ic = i1

               do while ( ic .le. i3 )
                call snum(chlw,ic,i3,ic2,prn,ierrtmp)
                if( ierrtmp .ne. 0 ) goto 330
                ic = ic2
                j = j+1
                frgee(mmlm+ks1+j) = prn
               end do

              end do

              if ( j .ne. nfrg*2 ) goto 331

C and frgdx((nei+1)*nfrg*nag): angular differential cross sections (mb/sr)
              j = 0
              do while ( j .lt. nfrg*iabs(nag) )
               iskip = 100
               do while ( iskip .ne. 0 )
                iskip = 0
                call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &               jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
               end do
               if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 332
               ic = i1

               do while ( ic .le. i3 )
                call snum(chlw,ic,i3,ic2,prn,ierrtmp)
                if( ierrtmp .ne. 0 ) goto 332
                ic = ic2
                if ( j.eq.0 .and. prn.lt.0d0 ) goto 328
                if ( prn .ge. 0d0 ) then
                 j = j+1
                 frgdx(mmlm+ks3+j) = prn
                else
                 n0repeat = idnint(dabs(prn))
                 do i0repeat = 1, n0repeat
                  j = j+1
                  frgdx(mmlm+ks3+j) = frgdx(mmlm+ks3+j-1)
                 end do
                end if
               end do

              end do

              if ( j .ne. nfrg*iabs(nag) ) goto 333


C case 3
C frgee((nei+1)*nfrg*2): mean value and dispersion for Gaussian (MeV)
C angular distribution is not given (isotropic is assumed)
             else if ( neo .eq. 0 .and. nag .eq. 0 ) then

              j = 0
              do while ( j .lt. nfrg*2 )
               iskip = 100
               do while ( iskip .ne. 0 )
                iskip = 0
                call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &               jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
               end do
               if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 330
               ic = i1

               do while ( ic .le. i3 )
                call snum(chlw,ic,i3,ic2,prn,ierrtmp)
                if( ierrtmp .ne. 0 ) goto 330
                ic = ic2
                j = j+1
                frgee(mmlm+ks1+j) = prn
               end do

              end do

              if ( j .ne. nfrg*2 ) goto 331


C case 4
C frgdx((nei+1)*nfrg*neo): energy differential cross sections (mb/MeV)
C angular distribution is not given (isotropic is assumed)
             else if ( neo .gt. 0 .and. nag .eq. 0 ) then

              j = 0
              do while ( j .lt. nfrg*neo )
               iskip = 100
               do while ( iskip .ne. 0 )
                iskip = 0
                call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &               jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
               end do
               if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 334
               ic = i1

               do while ( ic .le. i3 )
                call snum(chlw,ic,i3,ic2,prn,ierrtmp)
                if( ierrtmp .ne. 0 ) goto 334
                ic = ic2
                if ( j.eq.0 .and. prn.lt.0d0 ) goto 328
                if ( prn .ge. 0d0 ) then
                 j = j+1
                 frgdx(mmlm+ks3+j) = prn
                else
                 n0repeat = idnint(dabs(prn))
                 do i0repeat = 1, n0repeat
                  j = j+1
                  frgdx(mmlm+ks3+j) = frgdx(mmlm+ks3+j-1)
                 end do
                end if
               end do

              end do

              if ( j .ne. nfrg*neo ) goto 335


C case 5
C frgdd((nei+1)*nfrg*neo*nag): angular differential cross sections
C of discrete type (mb/sr)
             else if ( neo .lt. 0 .and. nag .ne. 0 ) then

              j = 0
              do while ( j .lt. nfrg*iabs(neo)*iabs(nag) )
               iskip = 100
               do while ( iskip .ne. 0 )
                iskip = 0
                call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &               jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
               end do
               if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 332
               ic = i1

               do while ( ic .le. i3 )
                call snum(chlw,ic,i3,ic2,prn,ierrtmp)
                if( ierrtmp .ne. 0 ) goto 332
                ic = ic2
                if ( j.eq.0 .and. prn.lt.0d0 ) goto 328
                if ( prn .ge. 0d0 ) then
                 j = j+1
                 frgdd(mmlm+ks2+j) = prn
                else
                 n0repeat = idnint(dabs(prn))
                 do i0repeat = 1, n0repeat
                  j = j+1
                  frgdd(mmlm+ks2+j) = frgdd(mmlm+ks2+j-1)
                 end do
                end if
               end do

              end do

              if ( j .ne. nfrg*iabs(neo)*iabs(nag) ) goto 336


C case 6
C frgdx((nei+1)*nfrg*neo): discrete cross sections (mb)
C angular distribution is not given (isotropic is assumed)
             else if ( neo .lt. 0 .and. nag .eq. 0 ) then

              j = 0
              do while ( j .lt. nfrg*iabs(neo) )
               iskip = 100
               do while ( iskip .ne. 0 )
                iskip = 0
                call readl(jsn,iot,dsin,idsi,ill,ilf,'#!$',
     &               jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierrtmp)
               end do
               if( ierrtmp.ne.0 .or. jpn.eq.3 ) goto 337
               ic = i1

               do while ( ic .le. i3 )
                call snum(chlw,ic,i3,ic2,prn,ierrtmp)
                if( ierrtmp .ne. 0 ) goto 337
                ic = ic2
                if ( j.eq.0 .and. prn.lt.0d0 ) goto 328
                if ( prn .ge. 0d0 ) then
                 j = j+1
                 frgdx(mmlm+ks3+j) = prn
                else
                 n0repeat = idnint(dabs(prn))
                 do i0repeat = 1, n0repeat
                  j = j+1
                  frgdx(mmlm+ks3+j) = frgdx(mmlm+ks3+j-1)
                 end do
                end if
               end do

              end do

              if ( j .ne. nfrg*iabs(neo) ) goto 338

             else

              goto 339

             end if


            end do

*-----------------------------------------------------------------------

            goto 399

*-----------------------------------------------------------------------
*        write normalized data for check
*-----------------------------------------------------------------------

            icheck = 0

            if( icheck .eq. 1 ) then

               iot = 6

               write(iot,'(i5)') ifgdf(i,2)
               write(iot,'(i5)') ifgdf(i,3)

               write(iot,'(i5)') nei
               write(iot,'(10(1pe11.4))')
     &            ( frgne(kne+j), j = 1, nei + 1 )
               write(iot,'(10(1pe11.4))')
     &            ( frgxs(kxs+j), j = 1, nei + 1 )

               write(iot,'(i5)') neo
               if( neo .gt. 0 ) then
                  write(iot,'(10(1pe11.4))')
     &                 ( frgef(kef+j), j = 1, neo + 1 )
               else if( neo .lt. 0 ) then
                  write(iot,'(10(1pe11.4))')
     &                 ( frgef(kef+j), j = 1, iabs(neo) )
               end if

               write(iot,'(i5)') nag
               if( nag .ne. 0 ) then
                  write(iot,'(10(1pe11.4))')
     &                 ( frgaf(kaf+j), j = 1, iabs(nag) + 1 )
               end if

               write(iot,'(i5)') nfrg
               write(iot,'(10(i6))')
     &              ( ifrgm(kim+j), j = 1, nfrg )

               do m = 1, nei + 1
                  write(iot,'()')
                  ksf = ifrge_ksf(ks0+m)
                  kk1 = ifrge_ks1(ks0+m)
                  kk2 = ifrge_ks2(ks0+m)
                  kk3 = ifrge_ks3(ks0+m)

                  write(iot,'(10(1pe11.4))')
     &                 (frgsf(ksf+j), j =1, nfrg )

                  if( neo .eq. 0 ) then
                     do k = 1, nfrg
                        write(iot,'(2(1pe11.4))')
     &                       frgee(kk1+(k-1)*2+1),
     &                       frgee(kk1+(k-1)*2+2)
                     end do
                  end if

                  if ( neo .ne. 0 .and. nag .ne. 0 ) then

                     do k = 1, nfrg
                     do l = 1, iabs(neo)
                        write(iot,'(10(1pe11.4))')
     &         ( frgdd(kk2+(k-1)*iabs(neo)*iabs(nag)+(l-1)*iabs(nag)+j),
     &            j=1,iabs(nag) )
                     end do
                     end do

                  else if ( neo .eq. 0 .and. nag .ne. 0 ) then

                     do k = 1, nfrg
                        write(iot,'(10(1pe11.4))')
     &                  ( frgdx(kk3+(k-1)*iabs(nag)+j), j=1,iabs(nag) )
                     end do

                  else if ( neo .ne. 0 .and. nag .eq. 0 ) then

                     do k = 1, nfrg
                        write(iot,'(10(1pe11.4))')
     &                  ( frgdx(kk3+(k-1)*iabs(neo)+j), j=1,iabs(neo) )
                     end do

                  end if

               end do

            end if

*-----------------------------------------------------------------------

            goto 399

*-----------------------------------------------------------------------

 300        continue
 1300       format('** Error : when opening ',i2,'-th frag data file,'
     &           /a, ' cannot be found.')
            write(io,1300) i, frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:12952/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1300) i, frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            exit

 302        continue
 1302       format('** Error : description of projectile is wrong'
     &           /'in ',a)
            write(io,1302) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:12963/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1302) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 303        continue
 1303       format('** Error : projectile is very exotic'
     &           /'in ',a)
            write(io,1303) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:12974/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1303) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 304        continue
 1304       format('** Error : projectile in ',a
     &           /'is not consistent with proj in [frag data].')
            write(io,1304) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:12985/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1304) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 305        continue
 1305       format('** Error : description of target is wrong'
     &           /'in ',a)
            write(io,1305) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:12996/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1305) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 306        continue
 1306       format('** Error : target is very exotic'
     &           /'in ',a)
            write(io,1306) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13007/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1306) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 307        continue
 1307       format('** Error : target in ',a
     &           /'is not consistent with targ in [frag data].')
            write(io,1307) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13018/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1307) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 308        continue
 1308       format('** Error : description of nei is wrong'
     &           /'in ',a)
            write(io,1308) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13029/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1308) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 309        continue
 1309       format('** Error : description of ein data is wrong'
     &           /'in ',a)
            write(io,1309) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13040/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1309) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 310        continue
 1310       format('** Error : the number of ein data in ',a
     &           /'is not consistent with nei+1.')
            write(io,1310) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13051/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1310) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 311        continue
 1311       format('** Error : description of totxs data is wrong'
     &           /'in ',a)
            write(io,1311) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13062/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1311) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 312        continue
 1312       format('** Error : the number of totxs data in ',a
     &           /'is not consistent with nei+1.')
            write(io,1312) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13073/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1312) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 313        continue
 1313       format('** Error : description of neo is wrong'
     &           /'in ',a)
            write(io,1313) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13084/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1313) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 314        continue
 1314       format('** Error : in ',a,
     &           /'when neo=<0 or nag=0, opt=4 cannot be used.')
            write(io,1314) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13095/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1314) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 315        continue
 1315       format('** Error : description of eout data is wrong'
     &           /'in ',a)
            write(io,1315) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13106/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1315) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 316        continue
 1316       format('** Error : the number of eout data in ',a
     &           /'is not consistent with neo+1.')
            write(io,1316) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13117/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1316) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 317        continue
 1317       format('** Error : the number of eout data in ',a
     &           /'is not consistent with |neo|.')
            write(io,1317) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13128/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1317) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 318        continue
 1318       format('** Error : description of nag is wrong'
     &           /'in ',a)
            write(io,1318) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13139/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1318) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 319        continue
 1319       format('** Error : description of angle data is wrong'
     &           /'in ',a)
            write(io,1319) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13150/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1319) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 320        continue
 1320       format('** Error : the number of angle data in ',a
     &           /'is not consistent with |nag|+1.')
            write(io,1320) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13161/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1320) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 321        continue
 1321       format('** Error : description of nfrg is wrong'
     &           /'in ',a)
            write(io,1321) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13172/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1321) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 322        continue
 1322       format('** Error : description of frag data is wrong'
     &           /'in ',a)
            write(io,1322) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13183/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1322) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 323        continue
 1323       format('** Error : the number of frag data in ',a
     &           /'is not consistent with nfrg.')
            write(io,1323) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13194/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1323) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 324        continue
 1324       format('** Error : nucleus of frag data is very exotic'
     &           /'in ',a)
            write(io,1324) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13205/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1324) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 325        continue
 1325       format('** Error : description of proxs data is wrong',
     &           /'for ',i5,'-th incident energy in ',a)
            write(io,1325) m,frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13216/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1325) m,frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 326        continue
 1326       format('** Error : the number of proxs data '
     &           ,'is not consistent with nfrg'
     &           /'for',i5,'-th incident energy in ',a)
            write(io,1326) m,frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13228/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1326) m,frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 327        continue
 1327       format('** Error : description of DDX data is wrong',
     &           /'for ',i5,'-th incident energy in ',a)
            write(io,1327) m,frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13239/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1327) m,frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 328        continue
 1328       format('** Error : the first value of XS data is negative',
     &           /'for ',i5,'-th incident energy in ',a)
            write(io,1328) m,frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13250/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1328) m,frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 329        continue
 1329       format('** Error : the number of DDX data '
     &           ,'is not consistent with nfrg*neo*|nag|'
     &           /'for ',i5,'-th incident energy in ',a)
            write(io,1329) m,frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13262/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1329) m,frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 330        continue
 1330       format('** Error : description of mean and width values '
     &           ,'for Gaussian is wrong',
     &           /'for ',i5,'-th incident energy in ',a)
            write(io,1330) m,frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13274/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1330) m,frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 331        continue
 1331       format('** Error : the number of mean and width values '
     &           ,'for Gaussian is not consistent with nfrg*2'
     &           /'for ',i5,'-th incident energy in ',a)
            write(io,1331) m,frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13286/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1331) m,frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 332        continue
 1332       format('** Error : description of angular differential XS '
     &           ,'is wrong'
     &           /'for ',i5,'-th incident energy in ',a)
            write(io,1332) m,frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13298/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1332) m,frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 333        continue
 1333       format('** Error : the number of angular differential XS '
     &           ,'is not consistent with nfrg*|nag|'
     &           /'for ',i5,'-th incident energy in ',a)
            write(io,1333) m,frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13310/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1333) m,frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 334        continue
 1334       format('** Error : description of energy differential XS '
     &           ,'is wrong'
     &           /'for ',i5,'-th incident energy in ',a)
            write(io,1334) m,frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13322/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1334) m,frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 335        continue
 1335       format('** Error : the number of energy differential XS '
     &           ,'is not consistent with nfrg*neo'
     &           /'for ',i5,'-th incident energy in ',a)
            write(io,1335) m,frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13334/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1335) m,frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 336        continue
 1336       format('** Error : the number of angular differential XS '
     &           ,'is not consistent with nfrg*|neo|*|nag|'
     &           /'for ',i5,'-th incident energy in ',a)
            write(io,1336) m,frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13346/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1336) m,frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 337        continue
 1337       format('** Error : description of discrete XS is wrong'
     &           /'for ',i5,'-th incident energy in ',a)
            write(io,1337) m,frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13357/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1337) m,frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 338        continue
 1338       format('** Error : the number of discrete XS '
     &           ,'is not consistent with nfrg*|neo|'
     &           /'for ',i5,'-th incident energy in ',a)
            write(io,1338) m,frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13369/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1338) m,frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 339        continue
 1339       format('** Error : combination of neo=',i5,' and nag=',i5
     &           ,'is not allowed'
     &           /'for ',i5,'-th incident energy in ',a)
            write(io,1339) neo,nag,m,frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13381/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1339) neo,nag,m,frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399

 341        continue
 1341       format('** Error : the number of angle data in ',a
     &           /'is not consistent with |nag|.')
            write(io,1341) frgfl(i)(1:ifgdf(i,4))
            ErrCha = ''
            MsgID = 'L:13392/R:setpag/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,1341) frgfl(i)(1:ifgdf(i,4))
            ierr = ierr + 1
            goto 399


*-----------------------------------------------------------------------

 399        close(iot)

            end if

         end do

*-----------------------------------------------------------------------

* record fragdata cross sections when opt=5
         call FragData_record(io,jo,ierr)

         call FragData_cumulative(io,jo,ierr)

*-----------------------------------------------------------------------


         init_frag = 0 ! flag for [frag data] S.Hashimoto (2015.8.25)

      end if

*-----------------------------------------------------------------------
*     splitting
*-----------------------------------------------------------------------

      if( isptn .gt. 0 ) then

            iwt = iwt + 1

         do ll = 1, isptn

*-----------------------------------------------------------------------
*              exchange reg1 and reg2
*-----------------------------------------------------------------------

               if( npreg(ll) .lt. 0 ) then

                        npreg(ll) = - npreg(ll)

                        kdsm = ipgrc(ll)
                        ldsm = -1
                        idsm = -1

                  do i = 1, npreg(ll)

                        ldsm = ldsm + 1
                        ntr1 = idas_ipgrc(kdsm+ldsm)
                        ldsm = ldsm + 1
                        mtr1 = idas_ipgrc(kdsm+ldsm)

                        mdsm = 1
                        call moddas_allocate_int(
     &                          mtr1, idas_ipgrc_temporary1)

                     do k = 1, mtr1

                        ldsm = ldsm + 1
                        idas_ipgrc_temporary1( mdsm + k - 1 )
     &                     = idas_ipgrc(kdsm+ldsm)

                     end do

                        ldsm = ldsm + 1
                        ntr2 = idas_ipgrc(kdsm+ldsm)
                        ldsm = ldsm + 1
                        mtr2 = idas_ipgrc(kdsm+ldsm)

                        ndsm = 1
                        call moddas_allocate_int(
     &                          mtr2, idas_ipgrc_temporary2)

                     do k = 1, mtr2

                        ldsm = ldsm + 1
                        idas_ipgrc_temporary2( ndsm + k - 1 )
     &                     = idas_ipgrc(kdsm+ldsm)

                     end do

                        idsm = idsm + 1
                        idas_ipgrc(kdsm+idsm) = ntr2
                        idsm = idsm + 1
                        idas_ipgrc(kdsm+idsm) = mtr2

                     do k = 1, mtr2

                        idsm = idsm + 1
                        idas_ipgrc(kdsm+idsm)
     &                     = idas_ipgrc_temporary2( ndsm + k - 1 )

                     end do

                        idsm = idsm + 1
                        idas_ipgrc(kdsm+idsm) = ntr1
                        idsm = idsm + 1
                        idas_ipgrc(kdsm+idsm) = mtr1

                     do k = 1, mtr1

                        idsm = idsm + 1
                        idas_ipgrc(kdsm+idsm)
     &                     = idas_ipgrc_temporary1( mdsm + k - 1 )

                     end do
                     call moddas_deallocate_int(idas_ipgrc_temporary1)
                     call moddas_deallocate_int(idas_ipgrc_temporary2)

                  end do

               end if

*-----------------------------------------------------------------------

                  mmmin = mmmax

                  kdsm = ipgrc(ll)
                  ldsm = -1

                  call moddas_reallocate_int(
     &                    6, ll, MAX_NUM_IPGRT, ipgrt, idas_ipgrt)
                  idsm  = ipgrt(ll)
                  jdsm  = -1

         do 221 m = 1, npreg(ll) * 2

                  ldsm = ldsm + 1
                  ntrn = idas_ipgrc(kdsm+ldsm)
                  ldsm = ldsm + 1
                  mtrn = idas_ipgrc(kdsm+ldsm)
                  kssm = kdsm + ldsm + 1
                  ldsm = ldsm + mtrn

                     igm = idsm + jdsm + 3

                  call tregion2(io,jo,ierrm,ntrn,mtrn,idas_ipgrc(kssm)
     &                          ,igm,MAX_NUM_IPGRT,idas_ipgrt)

                     if( ierrm .ne. 0 ) then

                        write(io,'(" **** above error in ",i3,
     &                  "-th splitting region")') ( m - 1 ) / 2 + 1

                        ErrCha = ''
                        MsgID = 'L:13543/R:setpag/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **** above error in ",i3,
     &                  "-th splitting region")') ( m - 1 ) / 2 + 1

                        ierr = ierr + 1
                        goto 221

                     end if

                  jdsm = jdsm + 1
                  idas_ipgrt(idsm+jdsm) = ntrn

                  jdsm = jdsm + 1
                  idas_ipgrt(idsm+jdsm) = mtrn

                  jdsm = jdsm + mtrn

  221    continue

               if( ierr .ne. 0 ) goto 999

               if( jdsm > MAX_NUM_INWWT ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.setpag@read00.f'
     &                    //' ?dimension over idas_ipgrt?'
     &                    //' jdsm > MAX_NUM_IPGRT'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_IPGRT@moddas.f=',MAX_NUM_IPGRT,')'
                  ErrID = 'L:13572/R:setpag/F:read00.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(
     &                 6, ll, jdsm+1, ipgrt, idas_ipgrt)

*-----------------------------------------------------------------------

         end do

      end if

*-----------------------------------------------------------------------
*       repeated collisions cannot be used with forced collisions
*       at the same time
*-----------------------------------------------------------------------

      if( ircln .ne. 0 .and. ifcln .ne. 0 ) then

          write(io,'(" **** Error : repeated collisions cannot ",
     &    "be used with forced collisions at the same time")')

          ErrCha = ''
          MsgID = 'L:13596/R:setpag/F:read00.f'
          call ErrWrite(MsgID, ErrCha)
          write(jo,'(" **** Error : repeated collisions cannot ",
     &    "be used with forced collisions at the same time")')

          goto 999

      end if

*-----------------------------------------------------------------------
*     repeated collisions
*-----------------------------------------------------------------------

      if( ircln .ne. 0 ) then

                  mmmin = mmmax

*-----------------------------------------------------------------------

         do k = 1, ircln

*-----------------------------------------------------------------------
*           check particle species
*-----------------------------------------------------------------------

            do m = 1, mnrcl(k,20)

                  j = mnrcl(k,m)

                  if( j .eq. 12 .or. j .eq. 13 ) then

                     write( *,'("** Warning : at ",i2,
     &                          "-th [repeated collision] section."/
     &                          "   the particle ",a8,
     &                          " is not available.")')
     &                          k, pname(j)

                     mnrcl(k,m) = 0

                  end if

                  ircls(j) = ircls(j) + 1

                  if( ircls(j) .gt. 1 ) then

                     write(io,'("** Error : at ",i2,
     &                          "-th [repeated collision] section."/
     &                          "   the particle ",a8,
     &                          " is defined in duplicate.")')
     &                          k, pname(j)

                     ErrCha = ''
                     MsgID = 'L:13648/R:setpag/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'("** Error : at ",i2,
     &                          "-th [repeated collision] section."/
     &                          "   the particle ",a8,
     &                          " is defined in duplicate.")')
     &                          k, pname(j)

                     ierr = ierr + 1

                  end if

                  iwt = iwt + 1

            end do

*-----------------------------------------------------------------------
*           region check
*-----------------------------------------------------------------------

                  kdsm = inrlc(k)
                  ldsm = 0

                  call moddas_reallocate_int(
     &                    6, k, MAX_NUM_INRLT, inrlt, idas_inrlt)
                  idsm = inrlt(k)
                  jdsm = 0

            do 101 m = 1, mnrcl(k,0)

                  ldsm = ldsm + 1
                  ntrn = idas_inrlc(kdsm+ldsm)
                  ldsm = ldsm + 1
                  mtrn = idas_inrlc(kdsm+ldsm)
                  kssm = kdsm + ldsm + 1
                  ldsm = ldsm + mtrn

               do l = 1, mtrn

                  ntrkg = idas_inrlc(kssm+l-1)

                  if( ntrkg .eq. 6000000 ) then

                        write(io,'(" **Error at ",i3,
     &                  "-th repeated cell in ",i2,
     &                  "-th [repeated collision] section."/
     &                  "  (all) cannot be used")') m, k

                        ErrCha = ''
                        MsgID = 'L:13697/R:setpag/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **Error at ",i3,
     &                  "-th repeated cell in ",i2,
     &                  "-th [repeated collision] section."/
     &                  "  (all) cannot be used")') m, k

                        ierr = ierr + 1

                  end if

               end do

               if( idas_inrlc(kssm) .eq. 1000000 ) then

                     ntrn = idas_inrlc( kssm + 1 )
                     mtrn = mtrn - 3
                     kssm = kssm + 2

               end if

                     igm = idsm + jdsm + 3

                  call tregion2(io,jo,ierrm,ntrn,mtrn,idas_inrlc(kssm),
     &                          igm,MAX_NUM_INRLT,idas_inrlt)

                     if( ierrm .ne. 0 ) then

                        write(io,'(" **** above error at ",i3,
     &                  "-th repeated cell in ",i2,
     &                  "-th [repeated collision] section.")') m, k

                        ErrCha = ''
                        MsgID = 'L:13730/R:setpag/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **** above error at ",i3,
     &                  "-th repeated cell in ",i2,
     &                  "-th [repeated collision] section.")') m, k

                        ierr = ierr + 1
                        goto 101

                     end if

                  jdsm = jdsm + 1
                  idas_inrlt(idsm+jdsm) = ntrn

                  jdsm = jdsm + 1
                  idas_inrlt(idsm+jdsm) = mtrn

                  jdsm = jdsm + mtrn

  101       continue

               if( ierr .ne. 0 ) goto 999

               if( jdsm > MAX_NUM_INRLT ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.setpag@read00.f ?dimension over idas_inrlt?'
     &                    //' jdsm > MAX_NUM_INRLT'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_INRLT@moddas.f=',MAX_NUM_INRLT,')'
                  ErrID = 'L:13759/R:setpag/F:read00.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(6, k, jdsm+1, inrlt, idas_inrlt)

*-----------------------------------------------------------------------
*           check the simple region
*-----------------------------------------------------------------------

                  jdsm = 0

            do 701 m = 1, mnrcl(k,0)

                  jdsm = jdsm + 1
                  ntrn = idas_inrlt(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_inrlt(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                     j  = 0

            do 501 ir = 1, ntrn

                     j = j + 1
                     nreg = idas_inrlt(kssm+j-1)

                  if( nreg .gt. 0 .and. nreg .lt. 1000000 ) then

                        mrrcl(idgr(nreg)) = 1

                        goto 501

                  else if( nreg .eq. 6000000 ) then

                        ierr = ierr + 1

                        goto 501

                  end if

*-----------------------------------------------------------------------

                        kpar = 0

                     do i = 0, 20

                        ipar(i) = 0
                        jpar(i) = 0
                        klev(i) = 0

                     end do

  601          continue

               if( nreg .gt. 0 .and. nreg .lt. 1000000 ) then

                        ipar(kpar) = ipar(kpar) + 1

                        mrrcl(idgr(nreg)) = 1

               else if( nreg .lt. 0 ) then

                        kpar = kpar + 1
                        jpar(kpar) = -nreg

               else if( nreg .gt. 3000000 .and. nreg .lt. 4000000 ) then

                        knum = nreg - 3000000

                        kpar = kpar + 1
                        jpar(kpar) = knum
                        klev(kpar) = 1

               else if( nreg .gt. 4000000 .and. nreg .lt. 5000000 ) then

                        kpar = kpar + 1
                        jpar(kpar) = -1
                        ipar(kpar) = -1
                        klev(kpar) = 0

                        nreg = nreg - 4000000

                        j = j + 1
                        jlat = idas_inrlt(kssm+j-1)

                     do ll = 1, jlat

                        do mm = 1, 6

                           j = j + 1

                        end do

                     end do

               end if

*-----------------------------------------------------------------------

  661          continue

               if( ipar(kpar) .eq. jpar(kpar) ) then

                     ipar(kpar) = 0
                     jpar(kpar) = 0

                  if( klev(kpar) .gt. 0 ) then

                     klev(kpar) = 0

                  end if

                     kpar = kpar - 1
                     ipar(kpar) = ipar(kpar) + 1

                  if( kpar .eq. 0 ) goto 501

                     goto 661

               end if

                     j = j + 1
                     nreg = idas_inrlt(kssm+j-1)

                     goto 601

  501       continue
  701       continue

         end do

*-----------------------------------------------------------------------

               if( ierr .ne. 0 ) goto 999

*-----------------------------------------------------------------------
*           storage space for analysis
*-----------------------------------------------------------------------

               do j = 1, ircln

                     idsm = inrlt(j)
                     jdsm = 0
                     ii = 0

                  do i = 1, mnrcl(j,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inrlt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inrlt(idsm+jdsm)
                     jdsm = jdsm + mtrn

                     ii = ii + ntrn

                  end do

                  if( ii .gt. maxrr ) maxrr = ii

               end do

               isrcl = 0

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*     forced collision
*-----------------------------------------------------------------------

      if( ifcln .ne. 0 ) then

                  mmmin = mmmax

*-----------------------------------------------------------------------

         do k = 1, ifcln

*-----------------------------------------------------------------------
*           check particle species
*-----------------------------------------------------------------------

            do m = 1, mnfcl(k,20)

                  j = mnfcl(k,m)

                  if( j .eq. 12 .or. j .eq. 13 ) then

                     write( *,'("** Warning : at ",i2,
     &                          "-th [forced collision] section."/
     &                          "   the particle ",a8,
     &                          " is not available.")')
     &                          k, pname(j)

                     mnfcl(k,m) = 0

                  end if

                  ifcls(j) = ifcls(j) + 1

                  if( ifcls(j) .gt. 1 ) then

                     write(io,'("** Error : at ",i2,
     &                          "-th [forced collision] section."/
     &                          "   the particle ",a8,
     &                          " is defined in duplicate.")')
     &                          k, pname(j)

                     ErrCha = ''
                     MsgID = 'L:13971/R:setpag/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'("** Error : at ",i2,
     &                          "-th [forced collision] section."/
     &                          "   the particle ",a8,
     &                          " is defined in duplicate.")')
     &                          k, pname(j)

                     ierr = ierr + 1

                  end if

                  iwt = iwt + 1

            end do

*-----------------------------------------------------------------------
*           region check
*-----------------------------------------------------------------------

                  kdsm = inflc(k)
                  ldsm = 0

                  call moddas_reallocate_int(
     &                    6, k, MAX_NUM_INFLT, inflt, idas_inflt)
                  idsm = inflt(k)
                  jdsm = 0

            do 100 m = 1, mnfcl(k,0)

                  ldsm = ldsm + 1
                  ntrn = idas_inflc(kdsm+ldsm)
                  ldsm = ldsm + 1
                  mtrn = idas_inflc(kdsm+ldsm)
                  kssm = kdsm + ldsm + 1
                  ldsm = ldsm + mtrn

               do l = 1, mtrn

                  ntrkg = idas_inflc(kssm+l-1)
                  if( ntrkg .eq. 6000000 ) then

                        write(io,'(" **Error at ",i3,
     &                  "-th forced cell in ",i2,
     &                  "-th [forced collision] section."/
     &                  "  (all) cannot be used")') m, k

                        ErrCha = ''
                        MsgID = 'L:14019/R:setpag/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **Error at ",i3,
     &                  "-th forced cell in ",i2,
     &                  "-th [forced collision] section."/
     &                  "  (all) cannot be used")') m, k

                        ierr = ierr + 1

                  end if

               end do

               if( idas_inflc(kssm) .eq. 1000000 ) then
                     ntrn = idas_inflc( kssm + 1 )
                     mtrn = mtrn - 3
                     kssm = kssm + 2

               end if

                     igm = idsm + jdsm + 3

                  call tregion2(io,jo,ierrm,ntrn,mtrn,idas_inflc(kssm)
     &                          ,igm, MAX_NUM_INFLT,idas_inflt)

                     if( ierrm .ne. 0 ) then

                        write(io,'(" **** above error at ",i3,
     &                  "-th forced cell in ",i2,
     &                  "-th [forced collision] section.")') m, k

                        ErrCha = ''
                        MsgID = 'L:14051/R:setpag/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **** above error at ",i3,
     &                  "-th forced cell in ",i2,
     &                  "-th [forced collision] section.")') m, k

                        ierr = ierr + 1
                        goto 100

                     end if

                  jdsm = jdsm + 1
                  idas_inflt(idsm+jdsm) = ntrn

                  jdsm = jdsm + 1
                  idas_inflt(idsm+jdsm) = mtrn

                  jdsm = jdsm + mtrn

  100       continue

               if( ierr .ne. 0 ) goto 999

               if( jdsm > MAX_NUM_INFLT ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.setpag@read00.f ?dimension over idas_inflt?'
     &                    //' jdsm > MAX_NUM_INFLT'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_INFLT@moddas.f=',MAX_NUM_INFLT,')'
                  ErrID = 'L:14080/R:setpag/F:read00.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(6, k, jdsm+1, inflt, idas_inflt)

*-----------------------------------------------------------------------
*           check the simple region
*-----------------------------------------------------------------------

                  jdsm = 0

            do 700 m = 1, mnfcl(k,0)

                  jdsm = jdsm + 1
                  ntrn = idas_inflt(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_inflt(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                     j  = 0

            do 500 ir = 1, ntrn

                     j = j + 1
                     nreg = idas_inflt(kssm+j-1)

                  if( nreg .gt. 0 .and. nreg .lt. 1000000 ) then

                        mrfcl(idgr(nreg)) = 1

                        goto 500

                  else if( nreg .eq. 6000000 ) then

                        ierr = ierr + 1

                        goto 500

                  end if

*-----------------------------------------------------------------------

                        kpar = 0

                     do i = 0, 20

                        ipar(i) = 0
                        jpar(i) = 0
                        klev(i) = 0

                     end do

  600          continue

               if( nreg .gt. 0 .and. nreg .lt. 1000000 ) then

                        ipar(kpar) = ipar(kpar) + 1

                        mrfcl(idgr(nreg)) = 1

               else if( nreg .lt. 0 ) then

                        kpar = kpar + 1
                        jpar(kpar) = -nreg

               else if( nreg .gt. 3000000 .and. nreg .lt. 4000000 ) then

                        knum = nreg - 3000000

                        kpar = kpar + 1
                        jpar(kpar) = knum
                        klev(kpar) = 1

               else if( nreg .gt. 4000000 .and. nreg .lt. 5000000 ) then

                        kpar = kpar + 1
                        jpar(kpar) = -1
                        ipar(kpar) = -1
                        klev(kpar) = 0

                        nreg = nreg - 4000000

                        j = j + 1
                        jlat = idas_inflt(kssm+j-1)

                     do ll = 1, jlat

                        do mm = 1, 6

                           j = j + 1

                        end do

                     end do

               end if

*-----------------------------------------------------------------------

  660          continue

               if( ipar(kpar) .eq. jpar(kpar) ) then

                     ipar(kpar) = 0
                     jpar(kpar) = 0

                  if( klev(kpar) .gt. 0 ) then

                     klev(kpar) = 0

                  end if

                     kpar = kpar - 1
                     ipar(kpar) = ipar(kpar) + 1

                  if( kpar .eq. 0 ) goto 500

                     goto 660

               end if

                     j = j + 1
                     nreg = idas_inflt(kssm+j-1)

                     goto 600

  500       continue
  700       continue

         end do

*-----------------------------------------------------------------------

               if( ierr .ne. 0 ) goto 999

*-----------------------------------------------------------------------
*           storage space for analysis
*-----------------------------------------------------------------------

               do j = 1, ifcln

                     idsm = inflt(j)
                     jdsm = 0
                     ii = 0

                  do i = 1, mnfcl(j,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inflt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inflt(idsm+jdsm)
                     jdsm = jdsm + mtrn

                     ii = ii + ntrn

                  end do

                  if( ii .gt. maxrg ) maxrg = ii

               end do

               isfcl = 0
               call moddas_allocate_dbl3(20, 2, maxrg, wtxfcl)
               wtxfcl_pointer => wtxfcl(:,1,1)

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*     weight window
*-----------------------------------------------------------------------

      if( iwwdp .ne. 0 ) then

            iwt = iwt + 1

         do k = 1, iwwdp

*-----------------------------------------------------------------------
*           check particle species
*-----------------------------------------------------------------------

            do m = 1, mnwwp(k,20)

                  j = mnwwp(k,m)

                  iwwin(j) = iwwin(j) + 1

                  if( iwwin(j) .gt. 1 ) then

                     write(io,'("** Error : at ",i2,
     &                          "-th [weight window] section."/
     &                          "   the particle ",a8,
     &                          " is defined in duplicate.")')
     &                          k, pname(j)

                     ErrCha = ''
                     MsgID = 'L:14280/R:setpag/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'("** Error : at ",i2,
     &                          "-th [weight window] section."/
     &                          "   the particle ",a8,
     &                          " is defined in duplicate.")')
     &                          k, pname(j)

                     ierr = ierr + 1

                  end if

            end do

*-----------------------------------------------------------------------
*           region check
*-----------------------------------------------------------------------

         if( iwmsh .eq. 1 ) then

                  kdsm = inwwc(k)
                  ldsm = 0

                  call moddas_reallocate_int(
     &                    6, k, MAX_NUM_INWWT, inwwt, idas_inwwt)
                  idsm = inwwt(k)
                  jdsm = 0

            do 120 m = 1, mnwwp(k,0)

                  ldsm = ldsm + 1
                  ntrn = idas_inwwc(kdsm+ldsm)
                  ldsm = ldsm + 1
                  mtrn = idas_inwwc(kdsm+ldsm)
                  kssm = kdsm + ldsm + 1
                  ldsm = ldsm + mtrn

               do l = 1, mtrn

                  ntrkg = idas_inwwc(kssm+l-1)

                  if( ntrkg .eq. 6000000 ) then

                        write(io,'(" **Error at ",i3,
     &                  "-th weight window cell in ",i2,
     &                  "-th [weight window] section."/
     &                  "  (all) cannot be used")') m, k

                        ErrCha = ''
                        MsgID = 'L:14329/R:setpag/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **Error at ",i3,
     &                  "-th weight window cell in ",i2,
     &                  "-th [weight window] section."/
     &                  "  (all) cannot be used")') m, k

                        ierr = ierr + 1

                  end if

               end do

               if( idas_inwwc(kssm) .eq. 1000000 ) then
                     ntrn = idas_inwwc( kssm + 1 )
                     mtrn = mtrn - 3
                     kssm = kssm + 2

               end if

                     igm = idsm + jdsm + 3

                  call tregion2(io,jo,ierrm,ntrn,mtrn,idas_inwwc(kssm)
     &                          ,igm,MAX_NUM_INWWT,idas_inwwt)

                     if( ierrm .ne. 0 ) then

                        write(io,'(" **** above error at ",i3,
     &                  "-th cell in ",i2,
     &                  "-th [weight window] section.")') m, k

                        ErrCha = ''
                        MsgID = 'L:14361/R:setpag/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **** above error at ",i3,
     &                  "-th cell in ",i2,
     &                  "-th [weight window] section.")') m, k

                        ierr = ierr + 1
                        goto 120

                     end if

                  jdsm = jdsm + 1
                  idas_inwwt(idsm+jdsm) = ntrn

                  jdsm = jdsm + 1
                  idas_inwwt(idsm+jdsm) = mtrn

                  jdsm = jdsm + mtrn

  120       continue

               if( ierr .ne. 0 ) goto 999

                  if( jdsm > MAX_NUM_INWWT ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.setpag@read00.f'
     &                       //' ?dimension over idas_inwwt?'
     &                       //' jdsm > MAX_NUM_INWWT'
     &                    ,' (jdsm=',jdsm,')'
     &                    ,' (MAX_NUM_INWWT@moddas.f=',MAX_NUM_INWWT,')'
                     ErrID = 'L:14391/R:setpag/F:read00.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                  call moddas_reduce_int(
     &                    6, k, jdsm+1, inwwt, idas_inwwt)

         elseif( iwmsh .eq. 4 ) then

                  kdsm = inwwc(k)
                  ldsm = kdsm+1
                  ndata = mnwwp(k,0)
                  mtrn = 3+idas_inwwc(ldsm+2)
                  ldsm = ldsm+mtrn

                  call moddas_allocate_int(MAX_NUM_INWWC,idas_temporary)
                  call moddas_allocate_int(nelemtot,jdata)

                  call ttetmesh4(io,jo,ntrn,
     &                 mtrn,idas_inwwc(inwwc(k)+1),
     &                 MAX_NUM_INWWC,idas_temporary,
     &                 ndata,idas_inwwc(ldsm+1),
     &                 k,iwwdp,jdata,ierr)

                  if( mtrn > MAX_NUM_INWWC ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.setpag@read00.f'
     &                       //' ?dimension over idas_inwwc?'
     &                       //' mtrn > MAX_NUM_INWWC'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_INWWC@moddas.f=',MAX_NUM_INWWC,')'
                     ErrID = 'L:14422/R:setpag/F:read00.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                  ntrn0=idas_temporary(4+idas_temporary(3))

                  call moddas_reallocate_int(
     &                    6, k, mtrn+ntrn0, inwwt, idas_inwwt)

                  idsm = inwwt(k)
                  jdsm = idsm+1

                  idas_inwwt(jdsm:jdsm+mtrn-1)
     &                 =idas_temporary(1:mtrn)
                  idas_inwwt(jdsm+mtrn:jdsm+mtrn+ntrn0-1)
     &                 =jdata(1:ntrn0)

                  call moddas_deallocate_int(jdata)
                  call moddas_deallocate_int(idas_temporary)

         end if

*-----------------------------------------------------------------------

         end do

               if( ierr .ne. 0 ) goto 999

*-----------------------------------------------------------------------
*           storage space for analysis
*-----------------------------------------------------------------------

         if( iwmsh .eq. 1 ) then

               do j = 1, iwwdp

                     idsm = inwwt(j)
                     jdsm = 0
                     ii = 0

                  do i = 1, mnwwp(j,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inwwt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inwwt(idsm+jdsm)
                     jdsm = jdsm + mtrn

                     ii = ii + ntrn

                  end do

                  if( ii .gt. maxww ) maxww = ii

               end do

         else if( iwmsh .eq. 3 ) then

               do j = 1, iwwdp

                     inx = iwxnm(j)
                     iny = iwynm(j)
                     inz = iwznm(j)

                     ii = inx * iny * inz + 1

                  if( ii .gt. maxww ) maxww = ii

               end do

         else if( iwmsh .eq. 4 ) then

               do j = 1, iwwdp

                  if( mnwwp(j,0) .gt. maxww ) maxww = mnwwp(j,0)+1

               end do

         end if

               iswwp = 0
               call moddas_allocate_dbl3(20, 14, maxww, wtxwwp)
               wtxwwp_pointer => wtxwwp(:,1,1)

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*     WW Bias
*-----------------------------------------------------------------------

      if( iwwbias .ne. 0 ) then

*-----------------------------------------------------------------------

      if( iwbdp .ne. 0 .and. iwwdp .ne. 0 ) then

                    jerr = 0

                  if( iwbdp .ne. iwwdp ) jerr = jerr + 1

             do k = 1, iwbdp

                  if( mnwbp(k, 0) .ne. mnwwp(k, 0) ) jerr = jerr + 1
                  if( mnwbp(k,20) .ne. mnwwp(k,20) ) jerr = jerr + 1

               do l = 1, mnwbp(k,20)

                  if( mnwbp(k,l) .ne. mnwwp(k,l) ) jerr = jerr + 1

               end do

               do l = 1, ienwb(k)

                  if( eenwb(k,l) .ne. eenww(k,l) ) jerr = jerr + 1

               end do

               if( iwmsh .eq. 3)then

                if( iwbxty(k) .ne. iwxty(k) ) jerr = jerr + 1
                if( iwbyty(k) .ne. iwyty(k) ) jerr = jerr + 1
                if( iwbzty(k) .ne. iwzty(k) ) jerr = jerr + 1
                if( iwbxnm(k) .ne. iwxnm(k) ) jerr = jerr + 1
                if( iwbynm(k) .ne. iwynm(k) ) jerr = jerr + 1
                if( iwbznm(k) .ne. iwznm(k) ) jerr = jerr + 1
                if( rwbxdl(k) .ne. rwxdl(k) ) jerr = jerr + 1
                if( rwbydl(k) .ne. rwydl(k) ) jerr = jerr + 1
                if( rwbzdl(k) .ne. rwzdl(k) ) jerr = jerr + 1
                if( rwbxmi(k) .ne. rwxmi(k) ) jerr = jerr + 1
                if( rwbymi(k) .ne. rwymi(k) ) jerr = jerr + 1
                if( rwbzmi(k) .ne. rwzmi(k) ) jerr = jerr + 1
                if( rwbxma(k) .ne. rwxma(k) ) jerr = jerr + 1
                if( rwbyma(k) .ne. rwyma(k) ) jerr = jerr + 1
                if( rwbzma(k) .ne. rwzma(k) ) jerr = jerr + 1

               elseif( iwmsh .eq. 4)then

                if( idas_inwbc(inwbc(k)+1) .ne. idas_inwwc(inwwc(k)+1) )
     &               jerr = jerr + 1

                  kdsm = inwbc(k)
                  ldsm = kdsm+1
                  ndata = mnwbp(k,0)
                  mtrn = 3+idas_inwbc(ldsm+2)
                  ldsm = ldsm+mtrn

                  call moddas_allocate_int(MAX_NUM_INWBC,idas_temporary)
                  call moddas_allocate_int(nelemtot,jdata)

                  call ttetmesh4(io,jo,ntrn,
     &                 mtrn,idas_inwbc(inwbc(k)+1),
     &                 MAX_NUM_INWBC,idas_temporary,
     &                 ndata,idas_inwbc(ldsm+1),
     &                 k,iwbdp,jdata,ierr)

                  if( mtrn > MAX_NUM_INWBC ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.setpag@read00.f'
     &                       //' ?dimension over idas_inwbc?'
     &                       //' mtrn > MAX_NUM_INWBC'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_INWBC@moddas.f=',MAX_NUM_INWBC,')'
                     ErrID = 'L:14586/R:setpag/F:read00.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                   ntrn0=idas_temporary(4+idas_temporary(3))

                  call moddas_reallocate_int(
     &                    6, k, mtrn+ntrn0, inwbt, idas_inwbt)

                  idsm = inwbt(k)
                  jdsm = idsm+1

                  idas_inwbt(jdsm:jdsm+mtrn-1)
     &                 =idas_temporary(1:mtrn)
                  idas_inwbt(jdsm+mtrn:jdsm+mtrn+ntrn0-1)
     &                 =jdata(1:ntrn0)

                  call moddas_deallocate_int(jdata)
                  call moddas_deallocate_int(idas_temporary)

               endif

             end do

               if( jerr .ne. 0 ) then

                     write(io,'("** Error : Mismatch in",
     &                          " [weight window] and [WW bias].")')

                     ErrCha = ''
                     MsgID = 'L:14616/R:setpag/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'("** Error : Mismatch in",
     &                          " [weight window] and [WW bias].")')

                     ierr = ierr + 1

               end if

               if( ierr .ne. 0 ) goto 999

*-----------------------------------------------------------------------

         do j = 1, iwwdp

                  if( iwmsh .eq. 1)then
                   nloop=kvlmax
                  elseif( iwmsh .eq. 3)then
                   nloop=iwxnm(j)*iwynm(j)*iwznm(j)
                  elseif( iwmsh .eq. 4)then
                   nloop=mnwwp(j,0)
                  endif

                  ice = abs( ienww(j) )
                  if( ice .le. 0 ) ice = 1

            do ll = 1, ice

                     kdsm = kfwwp(j)
                     kdsn = kfwbp(j)

               do m = 1, mnwwp(j,0)

                if( iwmsh .eq. 4)then
                 mm=idas_inwwt(inwwt(j)+1+m)
                 nn=idas_inwbt(inwbt(j)+1+m)
                else
                 mm=m
                 nn=m
                endif

                     wwt = das_kfwwp(kdsm+(ll-1)*nloop+mm-1)
                     wbt = das_kfwbp(kdsn+(ll-1)*nloop+nn-1)

                  if( wbt .gt. 0.0d0 ) then

                     das_kfwwp(kdsm+(ll-1)*nloop+m-1) = wwt / wbt

                  else

                     write(io,'("** Error : WW Bias shoule be",
     &                          " greater than zero")')

                     ErrCha = ''
                     MsgID = 'L:14670/R:setpag/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'("** Error : WW Bias shoule be",
     &                          " greater than zero")')

                     ierr = ierr + 1

                  end if

               end do

            end do

         end do

               if( ierr .ne. 0 ) goto 999

*-----------------------------------------------------------------------

      end if
      end if

*-----------------------------------------------------------------------
*     importance
*-----------------------------------------------------------------------

      if( iimpn .ne. 0 ) then

                  mmmin = mmmax

*-----------------------------------------------------------------------

         do k = 1, iimpn

*-----------------------------------------------------------------------
*           check particle species
*-----------------------------------------------------------------------

            do m = 1, mnimp(k,20)

                  j = mnimp(k,m)

                  icimp(j) = icimp(j) + 1

                  if( icimp(j) .gt. 1 ) then

                     write(io,'("** Error : at ",i2,
     &                          "-th [importance] section."/
     &                          "   the particle ",a8,
     &                          " is defined in duplicate.")')
     &                          k, pname(j)

                     ErrCha = ''
                     MsgID = 'L:14723/R:setpag/F:read00.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'("** Error : at ",i2,
     &                          "-th [importance] section."/
     &                          "   the particle ",a8,
     &                          " is defined in duplicate.")')
     &                          k, pname(j)

                     ierr = ierr + 1

                  end if

            end do

*-----------------------------------------------------------------------
*           region check
*-----------------------------------------------------------------------

                  kdsm = inimc(k)
                  ldsm = 0

                  call moddas_reallocate_int(
     &                    7, k, MAX_NUM_INIMT, inimt, idas_inimt)
                  idsm = inimt(k)
                  jdsm = 0

            do 110 m = 1, mnimp(k,0)

                  ldsm = ldsm + 1
                  ntrn = idas_inimc(kdsm+ldsm)
                  ldsm = ldsm + 1
                  mtrn = idas_inimc(kdsm+ldsm)
                  kssm = kdsm + ldsm + 1
                  ldsm = ldsm + mtrn

               do l = 1, mtrn

                  ntrkg = idas_inimc(kssm+l-1)

                  if( ntrkg .eq. 6000000 ) then

                        write(io,'(" **Error at ",i3,
     &                  "-th importance cell in ",i2,
     &                  "-th [importance] section."/
     &                  "  (all) cannot be used")') m, k

                        ErrCha = ''
                        MsgID = 'L:14770/R:setpag/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **Error at ",i3,
     &                  "-th importance cell in ",i2,
     &                  "-th [importance] section."/
     &                  "  (all) cannot be used")') m, k

                        ierr = ierr + 1

                  end if

               end do

               if( idas_inimc(kssm) .eq. 1000000 ) then
                     ntrn = idas_inimc( kssm + 1 )
                     mtrn = mtrn - 3
                     kssm = kssm + 2

               end if

                     igm = idsm + jdsm + 3

                  call tregion2(io,jo,ierrm,ntrn,mtrn,idas_inimc(kssm)
     &                          ,igm,MAX_NUM_INIMT,idas_inimt)

                     if( ierrm .ne. 0 ) then

                        write(io,'(" **** above error at ",i3,
     &                  "-th cell in ",i2,
     &                  "-th [importance] section.")') m, k

                        ErrCha = ''
                        MsgID = 'L:14802/R:setpag/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **** above error at ",i3,
     &                  "-th cell in ",i2,
     &                  "-th [importance] section.")') m, k

                        ierr = ierr + 1
                        goto 110

                     end if

                  jdsm = jdsm + 1
                  idas_inimt(idsm+jdsm) = ntrn

                  jdsm = jdsm + 1
                  idas_inimt(idsm+jdsm) = mtrn

                  jdsm = jdsm + mtrn

  110       continue

               if( ierr .ne. 0 ) goto 999

                  jjsm = jdsm
                  jssm = idsm + jjsm + 3

*-----------------------------------------------------------------------
*           add undefined region : imp = 1.0
*-----------------------------------------------------------------------

                     iadn = 0
                     ilev = 0

            do 520 i = 1, iregn

                     ireg = idrg(i)
                     jdsm = 0

               do 510 m = 1, mnimp(k,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inimt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inimt(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     jj = 0

                  do ir = 1, ntrn

                     call tregck(ireg,ilev,ilat,mtrn
     &                           ,idas_inimt(kssm),jj,icc)

                     if( icc .ne. 0 ) goto 520

                  end do

  510          continue

                     iadn = iadn + 1
                     idas_inimt(jssm+iadn-1) = ireg

  520       continue

                  jdsm = jjsm + 1
                  idas_inimt(idsm+jdsm) = iadn
                  jdsm = jdsm + 1
                  idas_inimt(idsm+jdsm) = iadn
                  jdsm = jdsm + iadn

                  if( jdsm > MAX_NUM_INIMT ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.setpag@read00.f'
     &                       //' ?dimension over idas_inimt?'
     &                       //' jdsm > MAX_NUM_INIMT'
     &                    ,' (jdsm=',jdsm,')'
     &                    ,' (MAX_NUM_INIMT@moddas.f=',MAX_NUM_INIMT,')'
                     ErrID = 'L:14880/R:setpag/F:read00.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                  call moddas_reduce_int(
     &                    7, k, jdsm+1, inimt, idas_inimt)

                  idsm = kfimp(k) + mnimp(k,0)
                  das_kfimp(idsm) = 1.0d0

*-----------------------------------------------------------------------

         end do

               if( ierr .ne. 0 ) goto 999

               iimps = iimpn
               iwt   = iwt + 1

*-----------------------------------------------------------------------
*     importance from [region] section
*-----------------------------------------------------------------------

      else if( iimpo .eq. 0 ) then

               iwt = iwt

*-----------------------------------------------------------------------

      else if( iimpo .ne. 0 ) then

            do i = 1, iregn

               if( abs( dimp(i) - 1.d0 ) .gt. 1.0d-8 ) goto 130

            end do

               goto 160

  130       continue

            do i = 1, iregn

               if( dimp(i) .ne. dimp(1) ) goto 170

            end do

               goto 160

  170       continue

               do j = 1, 19

                  icimp(j)  = 1

               end do

               iwt = iwt + 1

  160    continue

*-----------------------------------------------------------------------
*        set importance from [region] section
*-----------------------------------------------------------------------

         if( iwt .gt. 0 ) then

                  iimpn = 1

                  k = iimpn

                  mnimp(k,20) = 19

               do m = 1, mnimp(k,20)

                  mnimp(k,m) = m

               end do

                  mnimp(k,0) = iregn

                  call moddas_reallocate_int(
     &                    7, k, MAX_NUM_INIMC, inimc, idas_inimc)
                  idsm = inimc(k)
                  jdsm = 0

               do m = 1, mnimp(k,0)

                  mtrn = 1

                  jdsm = jdsm + 1
                  idas_inimc(idsm+jdsm) = mtrn
                  jdsm = jdsm + 1
                  idas_inimc(idsm+jdsm) = mtrn
                  jdsm = jdsm + 1
                  idas_inimc(idsm+jdsm) = idrg(m)

               end do

               call moddas_reduce_int(7, k, jdsm+1, inimc, idas_inimc)

               call moddas_reallocate_dbl(
     &                 7, k, mnimp(k,0)+1, kfimp, das_kfimp)
                  kdsm = kfimp(k)


               do m = 1, mnimp(k,0)

                  das_kfimp(kdsm-1+m) = dimp(m)

               end do

                  call moddas_reallocate_int(
     &                    7, k, MAX_NUM_INIMT, inimt, idas_inimt)
                  idsm = inimt(k)
                  jdsm = 0

               do m = 1, mnimp(k,0)

                  mtrn = 1

                  jdsm = jdsm + 1
                  idas_inimt(idsm+jdsm) = mtrn
                  jdsm = jdsm + 1
                  idas_inimt(idsm+jdsm) = mtrn
                  jdsm = jdsm + 1
                  idas_inimt(idsm+jdsm) = idrg(m)

               end do

                  iadn = 0

                  jdsm = jdsm + 1
                  idas_inimt(idsm+jdsm) = iadn
                  jdsm = jdsm + 1
                  idas_inimt(idsm+jdsm) = iadn

                  if( jdsm > MAX_NUM_INIMT ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.setpag@read00.f'
     &                       //' ?dimension over idas_inimt?'
     &                       //' jdsm > MAX_NUM_INIMT'
     &                    ,' (jdsm=',jdsm,')'
     &                    ,' (MAX_NUM_INIMT@moddas.f=',MAX_NUM_INIMT,')'
                     ErrID = 'L:15024/R:setpag/F:read00.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                  call moddas_reduce_int(
     &                    7, k, jdsm+1, inimt, idas_inimt)

                  idsm = kfimp(k) + mnimp(k,0)
                  das_kfimp(idsm) = 1.0d0

         end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*     set one for undefined particles
*-----------------------------------------------------------------------

      if( iwt .gt. 0 ) then

                  iadp = 0
                  k = iimpn + 1

            do j = 1, 19

               if( icimp(j) .eq. 0 ) then

                  iadp = iadp + 1
                  mnimp(k,iadp) = j

               end if

            end do

                  mnimp(k,20) = iadp
                  mnimp(k,0)  = 0

               if( iadp .gt. 0 ) then

                  mtrn = iregn

               else

                  mtrn = 0

               end if

                  call moddas_reallocate_int(
     &                    7, k, MAX_NUM_INIMT, inimt, idas_inimt)

                  idsm = inimt(k)
                  jdsm = 0

                  jdsm = jdsm + 1
                  idas_inimt(idsm+jdsm) = mtrn

                  jdsm = jdsm + 1
                  idas_inimt(idsm+jdsm) = mtrn

               do l = 1, mtrn

                  jdsm = jdsm + 1
                  idas_inimt(idsm+jdsm) = idrg(l)

               end do

                  if( jdsm > MAX_NUM_INIMT ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.setpag@read00.f'
     &                       //' ?dimension over idas_inimt?'
     &                       //' jdsm > MAX_NUM_INIMT'
     &                    ,' (jdsm=',jdsm,')'
     &                    ,' (MAX_NUM_INIMT@moddas.f=',MAX_NUM_INIMT,')'
                     ErrID = 'L:15099/R:setpag/F:read00.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                  call moddas_reduce_int(
     &                    7, k, jdsm+1, inimt, idas_inimt)

               call moddas_reallocate_dbl(7, k, 1, kfimp, das_kfimp)
                  idsm = kfimp(k)

                  das_kfimp(idsm-1+1) = 1.0d0

*-----------------------------------------------------------------------
*           storage space for analysis
*-----------------------------------------------------------------------

               do j = 1, iimpn + 1

                     idsm = inimt(j)
                     jdsm = 0
                     ii = 0

                  do i = 1, mnimp(j,0) + 1

                     jdsm = jdsm + 1
                     ntrn = idas_inimt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inimt(idsm+jdsm)
                     jdsm = jdsm + mtrn

                     ii = ii + ntrn

                  end do

                  if( ii .gt. maxip ) maxip = ii

               end do

               isimp = 0
               call moddas_allocate_dbl3(20, 7, maxip, wtximp)
               wtximp_pointer => wtximp(:,1,1)

               iswct = 0
               call moddas_allocate_dbl3(20, 8, maxip, wtxcut)
               wtxcut_pointer => wtxcut(:,1,1)

      end if

*-----------------------------------------------------------------------

      if( iwt .gt. 0 ) then

               isnum = max(maxww, max(maxip, maxrg))*11
               call moddas_allocate_dbl(isnum, das_isstr)

      end if

*-----------------------------------------------------------------------
*     counter
*-----------------------------------------------------------------------

      do k = 1, 3
      if( ncntc(k) .eq. 1 ) then

                  mmmin = mmmax

                  kdsm = incrc(k)
                  ldsm = 0

                  call moddas_reallocate_int(
     &                    3, k, MAX_NUM_INCRT, incrt, idas_incrt)
                  idsm     = incrt(k)
                  jdsm     = 0

         do 210 m = 1, ncreg(k)

                  ldsm = ldsm + 1
                  ntrn = idas_incrc(kdsm+ldsm)
                  ldsm = ldsm + 1
                  mtrn = idas_incrc(kdsm+ldsm)
                  kssm = kdsm + ldsm + 1
                  ldsm = ldsm + mtrn

                     igm = idsm + jdsm + 3

                  call tregion2(io,jo,ierrm,ntrn,mtrn,idas_incrc(kssm)
     &                          ,igm,MAX_NUM_INCRT,idas_incrt)

                     if( ierrm .ne. 0 ) then

                        write(io,'(" **** above error in ",i1,
     &                  " counter and ",i3,"-th region")') k, m

                        ErrCha = ''
                        MsgID = 'L:15193/R:setpag/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **** above error in ",i1,
     &                  " counter and ",i3,"-th region")') k, m

                        ierr = ierr + 1
                        goto 210

                     end if

                  jdsm = jdsm + 1
                  idas_incrt(idsm+jdsm) = ntrn

                  jdsm = jdsm + 1
                  idas_incrt(idsm+jdsm) = mtrn

                  jdsm = jdsm + mtrn

  210    continue

               if( ierr .ne. 0 ) goto 999

                  if( jdsm > MAX_NUM_INCRT ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.setpag@read00.f'
     &                       //' ?dimension over idas_incrt?'
     &                       //' jdsm > MAX_NUM_INCRT'
     &                    ,' (jdsm=',jdsm,')'
     &                    ,' (MAX_NUM_INCRT@moddas.f=',MAX_NUM_INCRT,')'
                     ErrID = 'L:15222/R:setpag/F:read00.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                  call moddas_reduce_int(
     &                    3, k, jdsm+1, incrt, idas_incrt)

      end if
      end do

*-----------------------------------------------------------------------
*     super mirror
*-----------------------------------------------------------------------

      if( nsreg .ne. 0 ) then

*-----------------------------------------------------------------------
*              exchange reg1 and reg2
*-----------------------------------------------------------------------

               if( nsreg .lt. 0 ) then

                        nsreg = - nsreg

                        kdsm = isgrc
                        ldsm = -1
                        idsm = -1

                  do i = 1, nsreg

                        ldsm = ldsm + 1
                        ntr1 = idas_isgrc(kdsm+ldsm)
                        ldsm = ldsm + 1
                        mtr1 = idas_isgrc(kdsm+ldsm)

                        mdsm = 1
                        call moddas_allocate_int(
     &                          mtr1, idas_isgrc_temporary1)

                     do k = 1, mtr1

                        ldsm = ldsm + 1
                        idas_isgrc_temporary1( mdsm + k - 1 )
     &                     = idas_isgrc(kdsm+ldsm)

                     end do

                        ldsm = ldsm + 1
                        ntr2 = idas_isgrc(kdsm+ldsm)
                        ldsm = ldsm + 1
                        mtr2 = idas_isgrc(kdsm+ldsm)

                        ndsm = 1
                        call moddas_allocate_int(
     &                          mtr2, idas_isgrc_temporary2)

                     do k = 1, mtr2

                        ldsm = ldsm + 1
                        idas_isgrc_temporary2( ndsm + k - 1 )
     &                     = idas_isgrc(kdsm+ldsm)

                     end do

                        idsm = idsm + 1
                        idas_isgrc(kdsm+idsm) = ntr2
                        idsm = idsm + 1
                        idas_isgrc(kdsm+idsm) = mtr2

                     do k = 1, mtr2

                        idsm = idsm + 1
                        idas_isgrc(kdsm+idsm)
     &                     = idas_isgrc_temporary2( ndsm + k - 1 )

                     end do

                        idsm = idsm + 1
                        idas_isgrc(kdsm+idsm) = ntr1
                        idsm = idsm + 1
                        idas_isgrc(kdsm+idsm) = mtr1

                     do k = 1, mtr1

                        idsm = idsm + 1
                        idas_isgrc(kdsm+idsm)
     &                     = idas_isgrc_temporary1( mdsm + k - 1 )

                     end do
                     call moddas_deallocate_int(idas_isgrc_temporary1)
                     call moddas_deallocate_int(idas_isgrc_temporary2)

                  end do

               end if

*-----------------------------------------------------------------------

                  mmmin = mmmax

                  kdsm = isgrc
                  ldsm = -1

                  isgrt = 1
                  iaddress_region(:) = 1
                  call moddas_reallocate_int(
     &                    2, 1, MAX_NUM_ISGRT
     &                    , iaddress_region, idas_isgrt)
                  idsm  = isgrt
                  jdsm  = -1

         do 220 m = 1, nsreg * 2

                  ldsm = ldsm + 1
                  ntrn = idas_isgrc(kdsm+ldsm)
                  ldsm = ldsm + 1
                  mtrn = idas_isgrc(kdsm+ldsm)
                  kssm = kdsm + ldsm + 1
                  ldsm = ldsm + mtrn

                     igm = idsm + jdsm + 3

                  call tregion2(io,jo,ierrm,ntrn,mtrn,idas_isgrc(kssm)
     &                          ,igm,MAX_NUM_ISGRT,idas_isgrt)

                     if( ierrm .ne. 0 ) then

                        write(io,'(" **** above error in ",i3,
     &                  "-th super mirror region")') ( m - 1 ) / 2 + 1

                        ErrCha = ''
                        MsgID = 'L:15353/R:setpag/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **** above error in ",i3,
     &                  "-th super mirror region")') ( m - 1 ) / 2 + 1

                        ierr = ierr + 1
                        goto 220

                     end if

                  jdsm = jdsm + 1
                  idas_isgrt(idsm+jdsm) = ntrn

                  jdsm = jdsm + 1
                  idas_isgrt(idsm+jdsm) = mtrn

                  jdsm = jdsm + mtrn

  220    continue

               if( ierr .ne. 0 ) goto 999

               if( jdsm > MAX_NUM_ISGRT ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.setpag@read00.f'
     &                    //' ?dimension over idas_isgrt?'
     &                    //' jdsm > MAX_NUM_ISGRT'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_ISGRT@moddas.f=',MAX_NUM_ISGRT,')'
                  ErrID = 'L:15382/R:setpag/F:read00.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(
     &                 2, 1, jdsm+1, iaddress_region, idas_isgrt)

               icrsflx = 1

      end if

*-----------------------------------------------------------------------
*     timer
*-----------------------------------------------------------------------

      if( ntmrg .ne. 0 ) then

                  mmmin = mmmax

                  kdsm = intmc
                  ldsm = 0

                  intmt = 0
                  iaddress_region(:) = 0
                  call moddas_reallocate_int(
     &                    2, 1, MAX_NUM_INTMT
     &                    , iaddress_region, idas_intmt)
                  idsm  = intmt
                  jdsm  = 0

         do 230 m = 1, ntmrg

                  ldsm = ldsm + 1
                  ntrn = idas_intmc(kdsm+ldsm)
                  ldsm = ldsm + 1
                  mtrn = idas_intmc(kdsm+ldsm)
                  kssm = kdsm + ldsm + 1
                  ldsm = ldsm + mtrn

                     igm = idsm + jdsm + 3

                  call tregion2(io,jo,ierrm,ntrn,mtrn,idas_intmc(kssm)
     &                          ,igm,MAX_NUM_INTMT,idas_intmt)

                     if( ierrm .ne. 0 ) then

                        write(io,'(" **** above error in ",i3,
     &                  "-th timer")') m

                        ErrCha = ''
                        MsgID = 'L:15432/R:setpag/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **** above error in ",i3,
     &                  "-th timer")') m

                        ierr = ierr + 1
                        goto 230

                     end if

                  jdsm = jdsm + 1
                  idas_intmt(idsm+jdsm) = ntrn

                  jdsm = jdsm + 1
                  idas_intmt(idsm+jdsm) = mtrn

                  jdsm = jdsm + mtrn

  230    continue

               if( ierr .ne. 0 ) goto 999

               if( jdsm > MAX_NUM_INTMT ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.setpag@read00.f'
     &                    //' ?dimension over idas_intmt?'
     &                    //' jdsm > MAX_NUM_INTMT'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_INTMT@moddas.f=',MAX_NUM_INTMT,')'
                  ErrID = 'L:15461/R:setpag/F:read00.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(
     &                 2, 1, jdsm+1, iaddress_region, idas_intmt)

      end if

*-----------------------------------------------------------------------
*     magnetic field
*-----------------------------------------------------------------------

      if( nmreg .ne. 0 .and. mstz(14) .ne. 0 ) then

                  mmmin = mmmax

                  kdsm = ingrc
                  ldsm = 0

                  ingrt = 0
                  iaddress_region(:) = 0
                  call moddas_reallocate_int(
     &                    2, 1, MAX_NUM_INGRT
     &                    , iaddress_region, idas_ingrt)
                  idsm  = ingrt
                  jdsm  = 0

         do 200 m = 1, nmreg

                  ldsm = ldsm + 1
                  ntrn = idas_ingrc(kdsm+ldsm)
                  ldsm = ldsm + 1
                  mtrn = idas_ingrc(kdsm+ldsm)
                  kssm = kdsm + ldsm + 1
                  ldsm = ldsm + mtrn

                     igm = idsm + jdsm + 3

                  call tregion2(io,jo,ierrm,ntrn,mtrn,idas_ingrc(kssm)
     &                          ,igm,MAX_NUM_INGRT,idas_ingrt)

                     if( ierrm .ne. 0 ) then

                        write(io,'(" **** above error in ",i3,
     &                  "-th magnet field")') m

                        ErrCha = ''
                        MsgID = 'L:15509/R:setpag/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **** above error in ",i3,
     &                  "-th magnet field")') m

                        ierr = ierr + 1
                        goto 200

                     end if

                  jdsm = jdsm + 1
                  idas_ingrt(idsm+jdsm) = ntrn

                  jdsm = jdsm + 1
                  idas_ingrt(idsm+jdsm) = mtrn

                  jdsm = jdsm + mtrn

  200    continue

               if( ierr .ne. 0 ) goto 999

                  if( jdsm > MAX_NUM_INGRT ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.setpag@read00.f'
     &                       //' ?dimension over idas_ingrt?'
     &                       //' jdsm > MAX_NUM_INGRT'
     &                    ,' (jdsm=',jdsm,')'
     &                    ,' (MAX_NUM_INGRT@moddas.f=',MAX_NUM_INGRT,')'
                     ErrID = 'L:15538/R:setpag/F:read00.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                  call moddas_reduce_int(
     &                    2, 1, jdsm+1, iaddress_region, idas_ingrt)

*-----------------------------------------------------------------------

                  kdsm = kmags
                  ldsm = 0
                  ltrc = 0

               do m = 1, nmreg

                  ldsm    = ldsm + 6
                  it_mag  = nint( das_kmags(kdsm+ldsm) )

                  if( it_mag .gt. 0 ) then

                        l = 0

                     do k = 1, igtrs

                        if( it_mag .eq. idtn(k) ) l = k

                     end do

                     if( l .gt. 0 ) then

                        das_kmags(kdsm+ldsm-1) = l

                     else

                        ierr = 1
                        write(io,'("*** Error: in [Magnetic Field] ",
     &                  "trcl =",i7," is not defined ",
     &                  "in [transform]")') it_mag

                        ErrCha = ''
                        MsgID = 'L:15578/R:setpag/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'("*** Error: in [Magnetic Field] ",
     &                  "trcl =",i7," is not defined ",
     &                  "in [transform]")') it_mag

                        goto 999

                     end if

                  end if

                  ldsm    = ldsm + 1

               end do

      end if

*-----------------------------------------------------------------------
*     electro magnetic field
*-----------------------------------------------------------------------

      if( nereg .ne. 0 .and. mstz(70) .ne. 0 ) then

                  mmmin = mmmax

                  kdsm = inerc
                  ldsm = 0

                  inert = 0
                  iaddress_region(:) = 0
                  call moddas_reallocate_int(
     &                    2, 1, MAX_NUM_INERT
     &                    , iaddress_region, idas_inert)
                  idsm  = inert
                  jdsm  = 0

         do 201 m = 1, nereg

                  ldsm = ldsm + 1
                  ntrn = idas_inerc(kdsm+ldsm)
                  ldsm = ldsm + 1
                  mtrn = idas_inerc(kdsm+ldsm)
                  kssm = kdsm + ldsm + 1
                  ldsm = ldsm + mtrn

                     igm = idsm + jdsm + 3

                  call tregion2(io,jo,ierrm,ntrn,mtrn,idas_inerc(kssm),
     &                          igm,MAX_NUM_INERT,idas_inert)

                     if( ierrm .ne. 0 ) then

                        write(io,'(" **** above error in ",i3,
     &                  "-th electro magnetic field")') m

                        ErrCha = ''
                        MsgID = 'L:15635/R:setpag/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **** above error in ",i3,
     &                  "-th electro magnetic field")') m

                        ierr = ierr + 1
                        goto 201

                     end if

                  jdsm = jdsm + 1
                  idas_inert(idsm+jdsm) = ntrn

                  jdsm = jdsm + 1
                  idas_inert(idsm+jdsm) = mtrn

                  jdsm = jdsm + mtrn

  201    continue

               if( ierr .ne. 0 ) goto 999

               if( jdsm > MAX_NUM_INERT ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.setpag@read00.f'
     &                    //' ?dimension over idas_inert?'
     &                    //' jdsm > MAX_NUM_INERT'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_INERT@moddas.f=',MAX_NUM_INERT,')'
                  ErrID = 'L:15664/R:setpag/F:read00.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(
     &                 2, 1, jdsm+1, iaddress_region, idas_inert)

*-----------------------------------------------------------------------

                  kdsm = kelcs
                  ldsm = 0
                  ltrc = 0

               do m = 1, nereg

                     ldsm    = ldsm + 4
                     it_elf  = nint( das_kelcs(kdsm+ldsm) )

                  if( it_elf .gt. 0 ) then

                        l = 0

                     do k = 1, igtrs

                        if( it_elf .eq. idtn(k) ) l = k

                     end do

                     if( l .gt. 0 ) then

                        das_kelcs(kdsm+ldsm-1) = l

                     else

                        ierr = 1
                        write(io,
     &                  '("*** Error: in [Electro Magnetic Field] ",
     &                  "trcle =",i5," is not defined ",
     &                  "in [transform]")') it_elf

                        ErrCha = ''
                        MsgID = 'L:15705/R:setpag/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,
     &                  '("*** Error: in [Electro Magnetic Field] ",
     &                  "trcle =",i5," is not defined ",
     &                  "in [transform]")') it_elf

                        goto 999

                     end if

                  end if

                     ldsm    = ldsm + 2
                     it_mgf  = nint( das_kelcs(kdsm+ldsm) )

                  if( it_mgf .gt. 0 ) then

                        l = 0

                     do k = 1, igtrs

                        if( it_mgf .eq. idtn(k) ) l = k

                     end do

                     if( l .gt. 0 ) then

                        das_kelcs(kdsm+ldsm-1) = l

                     else

                        ierr = 1
                        write(io,
     &                  '("*** Error: in [Electro Magnetic Field] ",
     &                  "trclm =",i7," is not defined ",
     &                  "in [transform]")') it_mgf

                        ErrCha = ''
                        MsgID = 'L:15744/R:setpag/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,
     &                  '("*** Error: in [Electro Magnetic Field] ",
     &                  "trclm =",i7," is not defined ",
     &                  "in [transform]")') it_mgf

                        goto 999

                     end if

                  end if

                  ldsm    = ldsm + 3

               end do

      end if

*-----------------------------------------------------------------------
*     region in [source]
*-----------------------------------------------------------------------

          do j = 1, imsrc

            if( nsrn(j) .gt. 0 ) then

                  mmmin = mmmax

                  kdsm = nsrc(j)
                  ldsm = 0

                  idsm  = iaddress_nsrn(j)
                  call moddas_reallocate_int(
     &                    isrc, j, MAX_NUM_NSRN, iaddress_nsrn
     &                    , idas_nsrn)
                  jdsm  = 0

                  ldsm = ldsm + 1
                  ntrn = idas_nsrc(kdsm+ldsm)
                  ldsm = ldsm + 1
                  mtrn = idas_nsrc(kdsm+ldsm)
                  kssm = kdsm + ldsm + 1
                  ldsm = ldsm + mtrn

                     igm = idsm + jdsm + 3

                  call tregion2(io,jo,ierrm,ntrn,mtrn,idas_nsrc(kssm)
     &                         ,igm,MAX_NUM_NSRN,idas_nsrn)

                     if( ierrm .ne. 0 ) then

                        write(io,'(" **** above error in ",
     &                  "reg = in [source]")')

                        ErrCha = ''
                        MsgID = 'L:15800/R:setpag/F:read00.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(" **** above error in ",
     &                  "reg = in [source]")')

                        ierr = ierr + 1
                        goto 999

                     end if

                  jdsm = jdsm + 1
                  idas_nsrn(idsm+jdsm) = ntrn

                  jdsm = jdsm + 1
                  idas_nsrn(idsm+jdsm) = mtrn

                  jdsm = jdsm + mtrn

               if( jdsm > MAX_NUM_NSRN ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.setpag@read00.f ?dimension over idas_nsrn?'
     &                    //' jdsm > MAX_NUM_NSRN'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_NSRN@moddas.f=',MAX_NUM_NSRN,')'
                  ErrID = 'L:15824/R:setpag/F:read00.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(
     &                 isrc, j, jdsm+1, iaddress_nsrn, idas_nsrn)

            end if

         end do

*-----------------------------------------------------------------------
*     surface source
*-----------------------------------------------------------------------

         do j = 1, imsrc
         if( jstyp(j) .eq. 26 ) then

                  icsf = 0
               do i = 1, 8
                  iscte(i) = 0
               end do

            rewind ioa

            do i = 1, igsuf

               read(ioa) idrf, idtr, idsf, igkst,
     &                  ( bval(k), k = 1, igkst )

               if( issuf(j) .eq. idsn(i) ) then

                  if( idsf .lt. 2 .or. idsf .gt. 15 ) then
                     write(io,'(/"*** Error, surface ",i7,
     &               " in surface source should be plane,",
     &               " sphere or cylinder")') issuf(j)
                     ierr = ierr + 1
                     goto 999
                  end if

                     issfd(j) = idsf
                     ivsfd(j) = igkst
                  do k = 1, igkst
                     dvsfd(j,k) = bval(k)
                  end do

                  icsf = 1

               else

                  do l = 1, iscut(j)
                  if( abs(isvct(j,l)) .eq. idsn(i) ) then

                     isvfd(j,l) = idsf
                     ivvfd(j,l) = igkst
                  do k = 1, igkst
                     dvvfd(j,l,k) = bval(k)
                  end do
                     iscte(l) = 1
                  end if
                  end do

               end if

            end do

            if( icsf .eq. 0 ) then
                  write(jo,'(/"*** Error, surface ",i7,
     &            " in surface source is not defined.")') issuf(j)
                  ierr = ierr + 1
                  goto 999
            end if

            do k = 1, iscut(j)
               if( iscte(k) .eq. 0 ) then
                  write(jo,'(/"*** Error, cut surface ",i7,
     &            " in surface source is not defined.")')
     &            abs(isvct(j,k))
                  ierr = ierr + 1
                  goto 999
               end if
            end do

         end if
         end do

*-----------------------------------------------------------------------
*        check and evaluate the surface and cut
*-----------------------------------------------------------------------

         do j = 1, imsrc
         if( jstyp(j) .eq. 26 ) then

*-----------------------------------------------------------------------

            if( issfd(j) .eq. 2 ) then
                  iss = 0
                  icc = 0
                  sxval = dvsfd(j,1)
                  isxdf(j) = 0
                  isydf(j) = 0
                  iszdf(j) = 0
                  isdef(j) = 0
                  icdef(j) = 0

               do k = 1, iscut(j)

*-----------------------------------------------------------------------
                  if( isvfd(j,k) .eq. 3 ) then
                     isydf(j) = isydf(j) + 1
                     iydef(j,isydf(j)) = isvct(j,k)
                     vypos(j,isydf(j)) = dvvfd(j,k,1)

*-----------------------------------------------------------------------
                  else if( isvfd(j,k) .eq. 4 ) then
                     iszdf(j) = iszdf(j) + 1
                     izdef(j,iszdf(j)) = isvct(j,k)
                     vzpos(j,iszdf(j)) = dvvfd(j,k,1)

*-----------------------------------------------------------------------
                  else if( isvfd(j,k) .ge. 5 .and.
     &                     isvfd(j,k) .le. 9 ) then
                     iss = iss + 1

                     if( iss .gt. 1 ) then
                        write(jo,'(/"*** Error, cut surface for ",i7,
     &                  " ,sphere or cylinder should be used once.")')
     &                  issuf(j)
                        ierr = ierr + 1
                        goto 999
                     end if
                     if( isvfd(j,k) .eq. 5 ) then
                              ssval = dvvfd(j,k,1)
                              sxmin = -ssval
                              sxmax =  ssval
                              sxpos(j) = 0.0
                              sypos(j) = 0.0
                              szpos(j) = 0.0

                     else if( isvfd(j,k) .eq. 6 ) then
                              ssval = dvvfd(j,k,4)
                              sxmin = dvvfd(j,k,1)-ssval
                              sxmax = dvvfd(j,k,1)+ssval
                              sxpos(j) = dvvfd(j,k,1)
                              sypos(j) = dvvfd(j,k,2)
                              szpos(j) = dvvfd(j,k,3)

                     else if( isvfd(j,k) .ge. 7 .and.
     &                              isvfd(j,k) .le. 9 ) then
                              ssval = dvvfd(j,k,2)
                           if( isvfd(j,k) .eq. 7 ) then
                              sxmin = dvvfd(j,k,1)-ssval
                              sxmax = dvvfd(j,k,1)+ssval
                              sxpos(j) = dvvfd(j,k,1)
                              sypos(j) = 0.0
                              szpos(j) = 0.0
                           else if( isvfd(j,k) .eq. 8 ) then
                              sxmin = -ssval
                              sxmax =  ssval
                              sxpos(j) = 0.0
                              sypos(j) = dvvfd(j,k,1)
                              szpos(j) = 0.0
                           else if( isvfd(j,k) .eq. 9 ) then
                              sxmin = -ssval
                              sxmax =  ssval
                              sxpos(j) = 0.0
                              sypos(j) = 0.0
                              szpos(j) = dvvfd(j,k,1)
                           end if
                     end if

*-----------------------------------------------------------------------

                     if( sxval .ge. sxmax .or.
     &                      sxval .le. sxmin ) then
                        if( isvct(j,k) .lt. 0 ) then
                           write(jo,'(/"*** Error, cut surface and",
     &                     " source surface is not overlap.")')
                           ierr = ierr + 1
                           goto 999
                        else
                           isdef(j) = 0
                        end if
                     else
                        if( isvct(j,k) .lt. 0 ) then
                           isdef(j) = -1
                        else
                           isdef(j) =  1
                        end if
                           ssrad(j) = sqrt( ssval**2 -
     &                                    ( sxval - sxpos(j) )**2 )
                     end if

*-----------------------------------------------------------------------

                  else if( isvfd(j,k) .ge. 10 .and.
     &                     isvfd(j,k) .le. 15 ) then
                     iss = iss + 1

                     if( iss .gt. 1 ) then
                        write(jo,'(/"*** Error, cut surface for ",i7,
     &                  " sphere or cylinder should be used once.")')
     &                  issuf(j)
                        ierr = ierr + 1
                        goto 999
                     end if

                     if( isvfd(j,k) .eq. 10 ) then
                              ssval = dvvfd(j,k,3)
                              sxpos(j) = 0.0
                              sypos(j) = dvvfd(j,k,1)
                              szpos(j) = dvvfd(j,k,2)
                              ssrad(j) = ssval
                     else if( isvfd(j,k) .eq. 11 ) then
                              ssval = dvvfd(j,k,3)
                              sxmin = dvvfd(j,k,1)-ssval
                              sxmax = dvvfd(j,k,1)+ssval
                              sxpos(j) = dvvfd(j,k,1)
                              sypos(j) = 0.0
                              szpos(j) = dvvfd(j,k,2)
                              ssrad(j) = ssval
                     else if( isvfd(j,k) .eq. 12 ) then
                              ssval = dvvfd(j,k,3)
                              sxmin = dvvfd(j,k,1)-ssval
                              sxmax = dvvfd(j,k,1)+ssval
                              sxpos(j) = dvvfd(j,k,1)
                              sypos(j) = dvvfd(j,k,2)
                              szpos(j) = 0.0
                              ssrad(j) = ssval
                     else if( isvfd(j,k) .eq. 13 ) then
                              ssval = dvvfd(j,k,1)
                              sxpos(j) = 0.0
                              sypos(j) = 0.0
                              szpos(j) = 0.0
                              ssrad(j) = ssval
                     else if( isvfd(j,k) .eq. 14 ) then
                              ssval = dvvfd(j,k,1)
                              sxmin = -ssval
                              sxmax =  ssval
                              sxpos(j) = 0.0
                              sypos(j) = 0.0
                              szpos(j) = 0.0
                              ssrad(j) = ssval
                     else if( isvfd(j,k) .eq. 15 ) then
                              ssval = dvvfd(j,k,1)
                              sxmin = -ssval
                              sxmax =  ssval
                              sxpos(j) = 0.0
                              sypos(j) = 0.0
                              szpos(j) = 0.0
                              ssrad(j) = ssval
                     end if

*-----------------------------------------------------------------------

                     if( isvfd(j,k) .eq. 10 .or.
     &                   isvfd(j,k) .eq. 13 ) then
                        if( isvct(j,k) .lt. 0 ) then
                           icdef(j) = -1
                        else
                           icdef(j) =  1
                        end if
                     else if( sxval .ge. sxmax .or.
     &                        sxval .le. sxmin ) then
                        if( isvct(j,k) .lt. 0 ) then
                           write(jo,'(/"*** Error, cut surface and",
     &                     " source surface is not overlap.")')
                           ierr = ierr + 1
                           goto 999
                        else
                           icdef(j) =  0
                        end if
                     else
                           icdef(j) = 0
                        if( isvfd(j,k) .eq. 11 .or.
     &                      isvfd(j,k) .eq. 14 ) then
                           zval = sqrt( ssrad(j)**2
     &                              - ( sxval - sxpos(j) )**2 )
                           iszdf(j) = iszdf(j) + 1
                           izdef(j,iszdf(j)) = 1
                           vzpos(j,iszdf(j)) = -zval - szpos(j)
                           iszdf(j) = iszdf(j) + 1
                           izdef(j,iszdf(j)) = -1
                           vzpos(j,iszdf(j)) =  zval - szpos(j)
                        else if( isvfd(j,k) .eq. 12 .or.
     &                           isvfd(j,k) .eq. 15 ) then
                           yval = sqrt( ssrad(j)**2
     &                              - ( sxval - sxpos(j) )**2 )
                           isydf(j) = isydf(j) + 1
                           iydef(j,isydf(j)) = 1
                           vypos(j,isydf(j)) = -yval - sypos(j)
                           isydf(j) = isydf(j) + 1
                           iydef(j,isydf(j)) = -1
                           vypos(j,isydf(j)) =  yval - sypos(j)
                        end if
                     end if

*-----------------------------------------------------------------------

                  else
                     write(jo,'(/"*** Error, cut surface for ",i7,
     &               " surface source shoud be py, pz, s or cx.")')
     &               issuf(j)
                     ierr = ierr + 1
                     goto 999
                  end if

               end do

*-----------------------------------------------------------------------

               do k = 1, iscut(j)

                           ymmin =  1.d+18
                           ymmax = -1.d+18
                     do  l = 1, isydf(j)
                        if( iydef(j,l) .gt. 0 .and.
     &                      vypos(j,l) .gt. ymmax ) ymmax = vypos(j,l)
                        if( iydef(j,l) .lt. 0 .and.
     &                      vypos(j,l) .lt. ymmin ) ymmin = vypos(j,l)
                     end do
                           isydf(j) = 0
                        if( ymmax .gt. -1.d+18 ) then
                           isydf(j) = isydf(j) + 1
                           iydef(j,isydf(j)) =  1
                           vypos(j,isydf(j)) = ymmax
                        end if
                        if( ymmin .lt. 1.d+18 ) then
                           isydf(j) = isydf(j) + 1
                           iydef(j,isydf(j)) = -1
                           vypos(j,isydf(j)) = ymmin
                        end if

                     if( isydf(j) .eq. 2 .and. iydef(j,1) .eq. -1 ) then
                           vtemp = vypos(j,1)
                           iydef(j,1) =  1
                           vypos(j,1) = vypos(j,2)
                           iydef(j,2) = -1
                           vypos(j,2) = vtemp
                     end if

                           zmmin =  1.d+18
                           zmmax = -1.d+18
                     do  l = 1, iszdf(j)
                        if( izdef(j,l) .gt. 0 .and.
     &                      vzpos(j,l) .gt. zmmax ) zmmax = vzpos(j,l)
                        if( izdef(j,l) .lt. 0 .and.
     &                      vzpos(j,l) .lt. zmmin ) zmmin = vzpos(j,l)
                     end do
                           iszdf(j) = 0
                        if( zmmax .gt. -1.d+18 ) then
                           iszdf(j) = iszdf(j) + 1
                           izdef(j,iszdf(j)) =  1
                           vzpos(j,iszdf(j)) = zmmax
                        end if
                        if( zmmin .lt. 1.d+18 ) then
                           iszdf(j) = iszdf(j) + 1
                           izdef(j,iszdf(j)) = -1
                           vzpos(j,iszdf(j)) = zmmin
                        end if

                     if( iszdf(j) .eq. 2 .and. izdef(j,1) .eq. -1 ) then
                           vtemp = vzpos(j,1)
                           izdef(j,1) =  1
                           vzpos(j,1) = vzpos(j,2)
                           izdef(j,2) = -1
                           vzpos(j,2) = vtemp
                     end if

               end do

*-----------------------------------------------------------------------

               if( isdef(j) .ne. -1 .and. icdef(j) .ne. -1 ) then
                  if( isydf(j) .ne. 2 .or. iszdf(j) .ne. 2 ) then
                     write(jo,'(/"*** Error, cut surface for ",i7,
     &               " surface source is poor.")')
     &               issuf(j)
                     ierr = ierr + 1
                     goto 999
                  end if
               end if

*-----------------------------------------------------------------------
*-----------------------------------------------------------------------

            else if( issfd(j) .eq. 3 ) then

                  iss = 0
                  icc = 0
                  syval = dvsfd(j,1)
                  isxdf(j) = 0
                  isydf(j) = 0
                  iszdf(j) = 0
                  isdef(j) = 0
                  icdef(j) = 0

               do k = 1, iscut(j)

*-----------------------------------------------------------------------
                  if( isvfd(j,k) .eq. 2 ) then
                     isxdf(j) = isxdf(j) + 1
                     ixdef(j,isxdf(j)) = isvct(j,k)
                     vxpos(j,isxdf(j)) = dvvfd(j,k,1)

*-----------------------------------------------------------------------
                  else if( isvfd(j,k) .eq. 4 ) then
                     iszdf(j) = iszdf(j) + 1
                     izdef(j,iszdf(j)) = isvct(j,k)
                     vzpos(j,iszdf(j)) = dvvfd(j,k,1)

*-----------------------------------------------------------------------
                  else if( isvfd(j,k) .ge. 5 .and.
     &                     isvfd(j,k) .le. 9 ) then
                     iss = iss + 1

                     if( iss .gt. 1 ) then
                        write(jo,'(/"*** Error, cut surface for ",i7,
     &                  " ,sphere or cylinder should be used once.")')
     &                  issuf(j)
                        ierr = ierr + 1
                        goto 999
                     end if
                     if( isvfd(j,k) .eq. 5 ) then
                              ssval = dvvfd(j,k,1)
                              symin = -ssval
                              symax =  ssval
                              sxpos(j) = 0.0
                              sypos(j) = 0.0
                              szpos(j) = 0.0

                     else if( isvfd(j,k) .eq. 6 ) then
                              ssval = dvvfd(j,k,4)
                              symin = dvvfd(j,k,2)-ssval
                              symax = dvvfd(j,k,2)+ssval
                              sxpos(j) = dvvfd(j,k,1)
                              sypos(j) = dvvfd(j,k,2)
                              szpos(j) = dvvfd(j,k,3)

                     else if( isvfd(j,k) .ge. 7 .and.
     &                              isvfd(j,k) .le. 9 ) then
                              ssval = dvvfd(j,k,2)
                           if( isvfd(j,k) .eq. 7 ) then
                              symin = -ssval
                              symax =  ssval
                              sxpos(j) = dvvfd(j,k,1)
                              sypos(j) = 0.0
                              szpos(j) = 0.0
                           else if( isvfd(j,k) .eq. 8 ) then
                              symin = dvvfd(j,k,1)-ssval
                              symax = dvvfd(j,k,1)+ssval
                              sxpos(j) = 0.0
                              sypos(j) = dvvfd(j,k,1)
                              szpos(j) = 0.0
                           else if( isvfd(j,k) .eq. 9 ) then
                              symin = -ssval
                              symax =  ssval
                              sxpos(j) = 0.0
                              sypos(j) = 0.0
                              szpos(j) = dvvfd(j,k,1)
                           end if
                     end if

*-----------------------------------------------------------------------

                     if( syval .ge. symax .or.
     &                      syval .le. symin ) then
                        if( isvct(j,k) .lt. 0 ) then
                           write(jo,'(/"*** Error, cut surface and",
     &                     " source surface is not overlap.")')
                           ierr = ierr + 1
                           goto 999
                        else
                           isdef(j) = 0
                        end if
                     else
                        if( isvct(j,k) .lt. 0 ) then
                           isdef(j) = -1
                        else
                           isdef(j) =  1
                        end if
                           ssrad(j) = sqrt( ssval**2 -
     &                                    ( syval - sypos(j) )**2 )
                     end if

*-----------------------------------------------------------------------

                  else if( isvfd(j,k) .ge. 10 .and.
     &                     isvfd(j,k) .le. 15 ) then
                     iss = iss + 1

                     if( iss .gt. 1 ) then
                        write(jo,'(/"*** Error, cut surface for ",i7,
     &                  " sphere or cylinder should be used once.")')
     &                  issuf(j)
                        ierr = ierr + 1
                        goto 999
                     end if

                     if( isvfd(j,k) .eq. 10 ) then
                              ssval = dvvfd(j,k,3)
                              symin = dvvfd(j,k,1)-ssval
                              symax = dvvfd(j,k,1)+ssval
                              sxpos(j) = 0.0
                              sypos(j) = dvvfd(j,k,1)
                              szpos(j) = dvvfd(j,k,2)
                              ssrad(j) = sval
                     else if( isvfd(j,k) .eq. 11 ) then
                              ssval = dvvfd(j,k,3)
                              sxpos(j) = dvvfd(j,k,1)
                              sypos(j) = 0.0
                              szpos(j) = dvvfd(j,k,2)
                              ssrad(j) = ssval
                     else if( isvfd(j,k) .eq. 12 ) then
                              ssval = dvvfd(j,k,3)
                              symin = dvvfd(j,k,1)-ssval
                              symax = dvvfd(j,k,1)+ssval
                              sxpos(j) = dvvfd(j,k,1)
                              sypos(j) = dvvfd(j,k,2)
                              szpos(j) = 0.0
                              ssrad(j) = ssval

                     else if( isvfd(j,k) .eq. 13 ) then
                              ssval = dvvfd(j,k,1)
                              symin = -ssval
                              symax =  ssval
                              sxpos(j) = 0.0
                              sypos(j) = 0.0
                              szpos(j) = 0.0
                              ssrad(j) = ssval
                     else if( isvfd(j,k) .eq. 14 ) then
                              ssval = dvvfd(j,k,1)
                              sxpos(j) = 0.0
                              sypos(j) = 0.0
                              szpos(j) = 0.0
                              ssrad(j) = ssval
                     else if( isvfd(j,k) .eq. 15 ) then
                              ssval = dvvfd(j,k,1)
                              symin = -ssval
                              symax =  ssval
                              sxpos(j) = 0.0
                              sypos(j) = 0.0
                              szpos(j) = 0.0
                              ssrad(j) = ssval
                     end if

*-----------------------------------------------------------------------

                     if( isvfd(j,k) .eq. 11 .or.
     &                   isvfd(j,k) .eq. 14 ) then
                        if( isvct(j,k) .lt. 0 ) then
                           icdef(j) = -1
                        else
                           icdef(j) =  1
                        end if
                     else if( syval .ge. symax .or.
     &                        syval .le. symin ) then
                        if( isvct(j,k) .lt. 0 ) then
                           write(jo,'(/"*** Error, cut surface and",
     &                     " source surface is not overlap.")')
                           ierr = ierr + 1
                           goto 999
                        else
                           icdef(j) =  0
                        end if
                     else
                           icdef(j) = 0
                        if( isvfd(j,k) .eq. 10 .or.
     &                      isvfd(j,k) .eq. 13 ) then
                           zval = sqrt( ssrad(j)**2
     &                              - ( syval - sypos(j) )**2 )
                           iszdf(j) = iszdf(j) + 1
                           izdef(j,iszdf(j)) = 1
                           vzpos(j,iszdf(j)) = -zval - szpos(j)
                           iszdf(j) = iszdf(j) + 1
                           izdef(j,iszdf(j)) = -1
                           vzpos(j,iszdf(j)) =  zval - szpos(j)
                        else if( isvfd(j,k) .eq. 12 .or.
     &                           isvfd(j,k) .eq. 15 ) then
                           xval = sqrt( ssrad(j)**2
     &                              - ( syval - sypos(j) )**2 )
                           isxdf(j) = isxdf(j) + 1
                           ixdef(j,isxdf(j)) = 1
                           vxpos(j,isxdf(j)) = -xval - sxpos(j)
                           isxdf(j) = isxdf(j) + 1
                           ixdef(j,isxdf(j)) = -1
                           vxpos(j,isxdf(j)) =  xval - sxpos(j)
                        end if
                     end if

*-----------------------------------------------------------------------

                  else
                     write(jo,'(/"*** Error, cut surface for ",i7,
     &               " surface source shoud be px, pz, s or c.")')
     &               issuf(j)
                     ierr = ierr + 1
                     goto 999
                  end if

               end do

*-----------------------------------------------------------------------

               do k = 1, iscut(j)

                           xmmin =  1.d+18
                           xmmax = -1.d+18
                     do  l = 1, isxdf(j)
                        if( ixdef(j,l) .gt. 0 .and.
     &                      vxpos(j,l) .gt. xmmax ) xmmax = vxpos(j,l)
                        if( ixdef(j,l) .lt. 0 .and.
     &                      vxpos(j,l) .lt. xmmin ) xmmin = vxpos(j,l)
                     end do
                           isxdf(j) = 0
                        if( xmmax .gt. -1.d+18 ) then
                           isxdf(j) = isxdf(j) + 1
                           ixdef(j,isxdf(j)) =  1
                           vxpos(j,isxdf(j)) = xmmax
                        end if
                        if( xmmin .lt. 1.d+18 ) then
                           isxdf(j) = isxdf(j) + 1
                           ixdef(j,isxdf(j)) = -1
                           vxpos(j,isxdf(j)) = xmmin
                        end if

                     if( isxdf(j) .eq. 2 .and. ixdef(j,1) .eq. -1 ) then
                           vtemp = vxpos(j,1)
                           ixdef(j,1) =  1
                           vxpos(j,1) = vxpos(j,2)
                           ixdef(j,2) = -1
                           vxpos(j,2) = vtemp
                     end if

                           zmmin =  1.d+18
                           zmmax = -1.d+18
                     do  l = 1, iszdf(j)
                        if( izdef(j,l) .gt. 0 .and.
     &                      vzpos(j,l) .gt. zmmax ) zmmax = vzpos(j,l)
                        if( izdef(j,l) .lt. 0 .and.
     &                      vzpos(j,l) .lt. zmmin ) zmmin = vzpos(j,l)
                     end do
                           iszdf(j) = 0
                        if( zmmax .gt. -1.d+18 ) then
                           iszdf(j) = iszdf(j) + 1
                           izdef(j,iszdf(j)) =  1
                           vzpos(j,iszdf(j)) = zmmax
                        end if
                        if( zmmin .lt. 1.d+18 ) then
                           iszdf(j) = iszdf(j) + 1
                           izdef(j,iszdf(j)) = -1
                           vzpos(j,iszdf(j)) = zmmin
                        end if

                     if( iszdf(j) .eq. 2 .and. izdef(j,1) .eq. -1 ) then
                           vtemp = vzpos(j,1)
                           izdef(j,1) =  1
                           vzpos(j,1) = vzpos(j,2)
                           izdef(j,2) = -1
                           vzpos(j,2) = vtemp
                     end if

               end do

*-----------------------------------------------------------------------

               if( isdef(j) .ne. -1 .and. icdef(j) .ne. -1 ) then
                  if( isxdf(j) .ne. 2 .or. iszdf(j) .ne. 2 ) then
                     write(jo,'(/"*** Error, cut surface for ",i7,
     &               " surface source is poor.")')
     &               issuf(j)
                     ierr = ierr + 1
                     goto 999
                  end if
               end if

*-----------------------------------------------------------------------
*-----------------------------------------------------------------------

            else if( issfd(j) .eq. 4 ) then

                  iss = 0
                  icc = 0
                  szval = dvsfd(j,1)
                  isxdf(j) = 0
                  isydf(j) = 0
                  iszdf(j) = 0
                  isdef(j) = 0
                  icdef(j) = 0

               do k = 1, iscut(j)

*-----------------------------------------------------------------------
                  if( isvfd(j,k) .eq. 2 ) then
                     isxdf(j) = isxdf(j) + 1
                     ixdef(j,isxdf(j)) = isvct(j,k)
                     vxpos(j,isxdf(j)) = dvvfd(j,k,1)

*-----------------------------------------------------------------------
                  else if( isvfd(j,k) .eq. 3 ) then
                     isydf(j) = isydf(j) + 1
                     iydef(j,isydf(j)) = isvct(j,k)
                     vypos(j,isydf(j)) = dvvfd(j,k,1)

*-----------------------------------------------------------------------
                  else if( isvfd(j,k) .ge. 5 .and.
     &                     isvfd(j,k) .le. 9 ) then
                     iss = iss + 1

                     if( iss .gt. 1 ) then
                        write(jo,'(/"*** Error, cut surface for ",i7,
     &                  " ,sphere or cylinder should be used once.")')
     &                  issuf(j)
                        ierr = ierr + 1
                        goto 999
                     end if
                     if( isvfd(j,k) .eq. 5 ) then
                              ssval = dvvfd(j,k,1)
                              szmin = -ssval
                              szmax =  ssval
                              sxpos(j) = 0.0
                              sypos(j) = 0.0
                              szpos(j) = 0.0

                     else if( isvfd(j,k) .eq. 6 ) then
                              ssval = dvvfd(j,k,4)
                              szmin = dvvfd(j,k,3)-ssval
                              szmax = dvvfd(j,k,3)+ssval
                              sxpos(j) = dvvfd(j,k,1)
                              sypos(j) = dvvfd(j,k,2)
                              szpos(j) = dvvfd(j,k,3)

                     else if( isvfd(j,k) .ge. 7 .and.
     &                              isvfd(j,k) .le. 9 ) then
                              ssval = dvvfd(j,k,2)
                           if( isvfd(j,k) .eq. 7 ) then
                              szmin = -ssval
                              szmax =  ssval
                              sxpos(j) = dvvfd(j,k,1)
                              sypos(j) = 0.0
                              szpos(j) = 0.0
                           else if( isvfd(j,k) .eq. 8 ) then
                              szmin = -ssval
                              szmax =  ssval
                              sxpos(j) = 0.0
                              sypos(j) = dvvfd(j,k,1)
                              szpos(j) = 0.0
                           else if( isvfd(j,k) .eq. 9 ) then
                              szmin = dvvfd(j,k,1)-ssval
                              szmax = dvvfd(j,k,1)+ssval
                              sxpos(j) = 0.0
                              sypos(j) = 0.0
                              szpos(j) = dvvfd(j,k,1)
                           end if
                     end if

*-----------------------------------------------------------------------

                     if( szval .ge. szmax .or.
     &                      szval .le. szmin ) then
                        if( isvct(j,k) .lt. 0 ) then
                           write(jo,'(/"*** Error, cut surface and",
     &                     " source surface is not overlap.")')
                           ierr = ierr + 1
                           goto 999
                        else
                           isdef(j) = 0
                        end if
                     else
                        if( isvct(j,k) .lt. 0 ) then
                           isdef(j) = -1
                        else
                           isdef(j) =  1
                        end if
                           ssrad(j) = sqrt( ssval**2 -
     &                                    ( szval - szpos(j) )**2 )
                     end if

*-----------------------------------------------------------------------

                  else if( isvfd(j,k) .ge. 10 .and.
     &                     isvfd(j,k) .le. 15 ) then
                     iss = iss + 1

                     if( iss .gt. 1 ) then
                        write(jo,'(/"*** Error, cut surface for ",i7,
     &                  " ,sphere or cylinder should be used once.")')
     &                  issuf(j)
                        ierr = ierr + 1
                        goto 999
                     end if

                     if( isvfd(j,k) .eq. 10 ) then
                              ssval = dvvfd(j,k,3)
                              szmin = dvvfd(j,k,1)-ssval
                              szmax = dvvfd(j,k,1)+ssval
                              sxpos(j) = 0.0
                              sypos(j) = dvvfd(j,k,1)
                              szpos(j) = dvvfd(j,k,2)
                              ssrad(j) = sval
                     else if( isvfd(j,k) .eq. 11 ) then
                              ssval = dvvfd(j,k,3)
                              szmin = dvvfd(j,k,1)-ssval
                              szmax = dvvfd(j,k,1)+ssval
                              sxpos(j) = dvvfd(j,k,1)
                              sypos(j) = 0.0
                              szpos(j) = dvvfd(j,k,2)
                              ssrad(j) = ssval
                     else if( isvfd(j,k) .eq. 12 ) then
                              ssval = dvvfd(j,k,3)
                              sxpos(j) = dvvfd(j,k,1)
                              sypos(j) = dvvfd(j,k,2)
                              szpos(j) = 0.0
                              ssrad(j) = ssval

                     else if( isvfd(j,k) .eq. 13 ) then
                              ssval = dvvfd(j,k,1)
                              szmin = -ssval
                              szmax =  ssval
                              sxpos(j) = 0.0
                              sypos(j) = 0.0
                              szpos(j) = 0.0
                              ssrad(j) = ssval
                     else if( isvfd(j,k) .eq. 14 ) then
                              ssval = dvvfd(j,k,1)
                              szmin = -ssval
                              szmax =  ssval
                              sxpos(j) = 0.0
                              sypos(j) = 0.0
                              szpos(j) = 0.0
                              ssrad(j) = ssval
                     else if( isvfd(j,k) .eq. 15 ) then
                              ssval = dvvfd(j,k,1)
                              sxpos(j) = 0.0
                              sypos(j) = 0.0
                              szpos(j) = 0.0
                              ssrad(j) = ssval
                     end if

*-----------------------------------------------------------------------

                     if( isvfd(j,k) .eq. 12 .or.
     &                   isvfd(j,k) .eq. 15 ) then
                        if( isvct(j,k) .lt. 0 ) then
                           icdef(j) = -1
                        else
                           icdef(j) =  1
                        end if
                     else if( szval .ge. szmax .or.
     &                        szval .le. szmin ) then
                        if( isvct(j,k) .lt. 0 ) then
                           write(jo,'(/"*** Error, cut surface and",
     &                     " source surface is not overlap.")')
                           ierr = ierr + 1
                           goto 999
                        else
                           icdef(j) =  0
                        end if
                     else
                           icdef(j) = 0
                        if( isvfd(j,k) .eq. 10 .or.
     &                      isvfd(j,k) .eq. 13 ) then
                           yval = sqrt( ssrad(j)**2
     &                              - ( szval - szpos(j) )**2 )
                           isydf(j) = isydf(j) + 1
                           iydef(j,isydf(j)) = 1
                           vypos(j,isydf(j)) = -yval - sypos(j)
                           isydf(j) = isydf(j) + 1
                           iydef(j,isydf(j)) = -1
                           vypos(j,isydf(j)) =  yval - sypos(j)
                        else if( isvfd(j,k) .eq. 11 .or.
     &                           isvfd(j,k) .eq. 14 ) then
                           xval = sqrt( ssrad(j)**2
     &                              - ( szval - szpos(j) )**2 )
                           isxdf(j) = isxdf(j) + 1
                           ixdef(j,isxdf(j)) = 1
                           vxpos(j,isxdf(j)) = -xval - sxpos(j)
                           isxdf(j) = isxdf(j) + 1
                           ixdef(j,isxdf(j)) = -1
                           vxpos(j,isxdf(j)) =  xval - sxpos(j)
                        end if
                     end if

*-----------------------------------------------------------------------

                  else
                     write(jo,'(/"*** Error, cut surface for ",i7,
     &               " surface source shoud be px, py, s or c.")')
     &               issuf(j)
                     ierr = ierr + 1
                     goto 999
                  end if

               end do

*-----------------------------------------------------------------------

               do k = 1, iscut(j)

                           xmmin =  1.d+18
                           xmmax = -1.d+18
                     do  l = 1, isxdf(j)
                        if( ixdef(j,l) .gt. 0 .and.
     &                      vxpos(j,l) .gt. xmmax ) xmmax = vxpos(j,l)
                        if( ixdef(j,l) .lt. 0 .and.
     &                      vxpos(j,l) .lt. xmmin ) xmmin = vxpos(j,l)
                     end do
                           isxdf(j) = 0
                        if( xmmax .gt. -1.d+18 ) then
                           isxdf(j) = isxdf(j) + 1
                           ixdef(j,isxdf(j)) =  1
                           vxpos(j,isxdf(j)) = xmmax
                        end if
                        if( xmmin .lt. 1.d+18 ) then
                           isxdf(j) = isxdf(j) + 1
                           ixdef(j,isxdf(j)) = -1
                           vxpos(j,isxdf(j)) = xmmin
                        end if

                     if( isxdf(j) .eq. 2 .and. ixdef(j,1) .eq. -1 ) then
                           vtemp = vxpos(j,1)
                           ixdef(j,1) =  1
                           vxpos(j,1) = vxpos(j,2)
                           ixdef(j,2) = -1
                           vxpos(j,2) = vtemp
                     end if
                           ymmin =  1.d+18
                           ymmax = -1.d+18
                     do  l = 1, isydf(j)
                        if( iydef(j,l) .gt. 0 .and.
     &                      vypos(j,l) .gt. ymmax ) ymmax = vypos(j,l)
                        if( iydef(j,l) .lt. 0 .and.
     &                      vypos(j,l) .lt. ymmin ) ymmin = vypos(j,l)
                     end do
                           isydf(j) = 0
                        if( ymmax .gt. -1.d+18 ) then
                           isydf(j) = isydf(j) + 1
                           iydef(j,isydf(j)) =  1
                           vypos(j,isydf(j)) = ymmax
                        end if
                        if( ymmin .lt. 1.d+18 ) then
                           isydf(j) = isydf(j) + 1
                           iydef(j,isydf(j)) = -1
                           vypos(j,isydf(j)) = ymmin
                        end if

                     if( isydf(j) .eq. 2 .and. iydef(j,1) .eq. -1 ) then
                           vtemp = vypos(j,1)
                           iydef(j,1) =  1
                           vypos(j,1) = vypos(j,2)
                           iydef(j,2) = -1
                           vypos(j,2) = vtemp
                     end if

               end do

*-----------------------------------------------------------------------

               if( isdef(j) .ne. -1 .and. icdef(j) .ne. -1 ) then
                  if( isxdf(j) .ne. 2 .or. isydf(j) .ne. 2 ) then
                     write(jo,'(/"*** Error, cut surface for ",i7,
     &               " surface source is poor.")')
     &               issuf(j)
                     ierr = ierr + 1
                     goto 999
                  end if
               end if

*-----------------------------------------------------------------------
*-----------------------------------------------------------------------

            else if( issfd(j) .ge. 5 .and. issfd(j) .le. 9 ) then

                  iss = 0
                  icc = 0
                  isxdf(j) = 0
                  isydf(j) = 0
                  iszdf(j) = 0
                  isdef(j) = 0
                  icdef(j) = 0

               if( issfd(j) .eq. 5 ) then
                  sxpos(j) = 0.0
                  sypos(j) = 0.0
                  szpos(j) = 0.0
                  ssrad(j) = dvsfd(j,1)

               else if( issfd(j) .eq. 6 ) then
                  sxpos(j) = dvsfd(j,1)
                  sypos(j) = dvsfd(j,2)
                  szpos(j) = dvsfd(j,3)
                  ssrad(j) = dvsfd(j,4)

               else if( issfd(j) .eq. 7 ) then
                  sxpos(j) = dvsfd(j,1)
                  sypos(j) = 0.0
                  szpos(j) = 0.0
                  ssrad(j) = dvsfd(j,2)

               else if( issfd(j) .eq. 8 ) then
                  sxpos(j) = 0.0
                  sypos(j) = dvsfd(j,1)
                  szpos(j) = 0.0
                  ssrad(j) = dvsfd(j,2)

               else if( issfd(j) .eq. 9 ) then
                  sxpos(j) = 0.0
                  sypos(j) = 0.0
                  szpos(j) = dvsfd(j,1)
                  ssrad(j) = dvsfd(j,2)

               end if

*-----------------------------------------------------------------------

               do k = 1, iscut(j)

                  if( isvfd(j,k) .eq. 2 ) then
                     isxdf(j) = isxdf(j) + 1
                     ixdef(j,isxdf(j)) = isvct(j,k)
                     vxpos(j,isxdf(j)) = dvvfd(j,k,1)

                  else if( isvfd(j,k) .eq. 3 ) then
                     isydf(j) = isydf(j) + 1
                     iydef(j,isydf(j)) = isvct(j,k)
                     vypos(j,isydf(j)) = dvvfd(j,k,1)

                  else if( isvfd(j,k) .eq. 4 ) then
                     iszdf(j) = iszdf(j) + 1
                     izdef(j,iszdf(j)) = isvct(j,k)
                     vzpos(j,iszdf(j)) = dvvfd(j,k,1)

                  else
                     write(jo,'(/"*** Error, cut surface for ",i7,
     &               " surface source shoud be px, py, pz.")')
     &               issuf(j)
                     ierr = ierr + 1
                     goto 999
                  end if

               end do

*-----------------------------------------------------------------------

                     xmmin =  1.d+18
                     xmmax = -1.d+18
               do  l = 1, isxdf(j)
                  if( ixdef(j,l) .gt. 0 .and.
     &                vxpos(j,l) .gt. xmmax ) xmmax = vxpos(j,l)
                  if( ixdef(j,l) .lt. 0 .and.
     &                vxpos(j,l) .lt. xmmin ) xmmin = vxpos(j,l)
               end do

                     isxdf(j) = 0
                  if( xmmax .gt. -1.d+18 ) then
                     isxdf(j) = isxdf(j) + 1
                     ixdef(j,isxdf(j)) =  1
                     vxpos(j,isxdf(j)) = xmmax
                  end if
                  if( xmmin .lt. 1.d+18 ) then
                     isxdf(j) = isxdf(j) + 1
                     ixdef(j,isxdf(j)) = -1
                     vxpos(j,isxdf(j)) = xmmin
                  end if

               if( isxdf(j) .eq. 2 .and. ixdef(j,1) .eq. -1 ) then
                     vtemp = vxpos(j,1)
                     ixdef(j,1) =  1
                     vxpos(j,1) = vxpos(j,2)
                     ixdef(j,2) = -1
                     vxpos(j,2) = vtemp
               end if

                     ymmin =  1.d+18
                     ymmax = -1.d+18
               do  l = 1, isydf(j)
                  if( iydef(j,l) .gt. 0 .and.
     &                vypos(j,l) .gt. ymmax ) ymmax = vypos(j,l)
                  if( iydef(j,l) .lt. 0 .and.
     &                vypos(j,l) .lt. ymmin ) ymmin = vypos(j,l)
               end do
                     isydf(j) = 0
                  if( ymmax .gt. -1.d+18 ) then
                     isydf(j) = isydf(j) + 1
                     iydef(j,isydf(j)) =  1
                     vypos(j,isydf(j)) = ymmax
                  end if
                  if( ymmin .lt. 1.d+18 ) then
                     isydf(j) = isydf(j) + 1
                     iydef(j,isydf(j)) = -1
                     vypos(j,isydf(j)) = ymmin
                  end if

               if( isydf(j) .eq. 2 .and. iydef(j,1) .eq. -1 ) then
                     vtemp = vypos(j,1)
                     iydef(j,1) =  1
                     vypos(j,1) = vypos(j,2)
                     iydef(j,2) = -1
                     vypos(j,2) = vtemp
               end if

                     zmmin =  1.d+18
                     zmmax = -1.d+18
               do  l = 1, iszdf(j)
                  if( izdef(j,l) .gt. 0 .and.
     &                vzpos(j,l) .gt. zmmax ) zmmax = vzpos(j,l)
                  if( izdef(j,l) .lt. 0 .and.
     &                vzpos(j,l) .lt. zmmin ) zmmin = vzpos(j,l)
               end do
                     iszdf(j) = 0
                  if( zmmax .gt. -1.d+18 ) then
                     iszdf(j) = iszdf(j) + 1
                     izdef(j,iszdf(j)) =  1
                     vzpos(j,iszdf(j)) = zmmax
                  end if
                  if( zmmin .lt. 1.d+18 ) then
                     iszdf(j) = iszdf(j) + 1
                     izdef(j,iszdf(j)) = -1
                     vzpos(j,iszdf(j)) = zmmin
                  end if

               if( iszdf(j) .eq. 2 .and. izdef(j,1) .eq. -1 ) then
                     vtemp = vzpos(j,1)
                     izdef(j,1) =  1
                     vzpos(j,1) = vzpos(j,2)
                     izdef(j,2) = -1
                     vzpos(j,2) = vtemp
               end if

*-----------------------------------------------------------------------
*-----------------------------------------------------------------------

            else if( issfd(j) .ge. 10 .and. issfd(j) .le. 15 ) then

                  iss = 0
                  icc = 0
                  isxdf(j) = 0
                  isydf(j) = 0
                  iszdf(j) = 0
                  isdef(j) = 0
                  icdef(j) = 0

               if( issfd(j) .eq. 10 ) then
                  sxpos(j) = 0.0
                  sypos(j) = dvsfd(j,1)
                  szpos(j) = dvsfd(j,2)
                  ssrad(j) = dvsfd(j,3)

               else if( issfd(j) .eq. 11 ) then
                  sxpos(j) = dvsfd(j,1)
                  sypos(j) = 0.0
                  szpos(j) = dvsfd(j,2)
                  ssrad(j) = dvsfd(j,3)

               else if( issfd(j) .eq. 12 ) then
                  sxpos(j) = dvsfd(j,1)
                  sypos(j) = dvsfd(j,2)
                  szpos(j) = 0.0
                  ssrad(j) = dvsfd(j,3)

               else if( issfd(j) .ge. 13 .and.
     &                  issfd(j) .le. 15) then
                  sxpos(j) = 0.0
                  sypos(j) = 0.0
                  szpos(j) = 0.0
                  ssrad(j) = dvsfd(j,1)

               end if

*-----------------------------------------------------------------------

               do k = 1, iscut(j)

                  if( isvfd(j,k) .eq. 2 ) then
                     isxdf(j) = isxdf(j) + 1
                     ixdef(j,isxdf(j)) = isvct(j,k)
                     vxpos(j,isxdf(j)) = dvvfd(j,k,1)

                  else if( isvfd(j,k) .eq. 3 ) then
                     isydf(j) = isydf(j) + 1
                     iydef(j,isydf(j)) = isvct(j,k)
                     vypos(j,isydf(j)) = dvvfd(j,k,1)

                  else if( isvfd(j,k) .eq. 4 ) then
                     iszdf(j) = iszdf(j) + 1
                     izdef(j,iszdf(j)) = isvct(j,k)
                     vzpos(j,iszdf(j)) = dvvfd(j,k,1)

                  else
                     write(jo,'(/"*** Error, cut surface for ",i7,
     &               " surface source shoud be px, py, pz.")')
     &               issuf(j)
                     ierr = ierr + 1
                     goto 999
                  end if

               end do

*-----------------------------------------------------------------------

                     xmmin =  1.d+18
                     xmmax = -1.d+18
               do  l = 1, isxdf(j)
                  if( ixdef(j,l) .gt. 0 .and.
     &                vxpos(j,l) .gt. xmmax ) xmmax = vxpos(j,l)
                  if( ixdef(j,l) .lt. 0 .and.
     &                vxpos(j,l) .lt. xmmin ) xmmin = vxpos(j,l)
               end do
                     isxdf(j) = 0
                  if( xmmax .gt. -1.d+18 ) then
                     isxdf(j) = isxdf(j) + 1
                     ixdef(j,isxdf(j)) =  1
                     vxpos(j,isxdf(j)) = xmmax
                  end if
                  if( xmmin .lt. 1.d+18 ) then
                     isxdf(j) = isxdf(j) + 1
                     ixdef(j,isxdf(j)) = -1
                     vxpos(j,isxdf(j)) = xmmin
                  end if

               if( isxdf(j) .eq. 2 .and. ixdef(j,1) .eq. -1 ) then
                     vtemp = vxpos(j,1)
                     ixdef(j,1) =  1
                     vxpos(j,1) = vxpos(j,2)
                     ixdef(j,2) = -1
                     vxpos(j,2) = vtemp
               end if

               if( ( issfd(j) .eq. 10 .or. issfd(j) .eq. 13 )
     &               .and. isxdf(j) .ne. 2 ) then
                  write(jo,'(/"*** Error, cut surface for ",i7,
     &            " surface source shoud have 2 px cut.")')
     &            issuf(j)
                  ierr = ierr + 1
                  goto 999
               end if

                     ymmin =  1.d+18
                     ymmax = -1.d+18
               do  l = 1, isydf(j)
                  if( iydef(j,l) .gt. 0 .and.
     &                vypos(j,l) .gt. ymmax ) ymmax = vypos(j,l)
                  if( iydef(j,l) .lt. 0 .and.
     &                vypos(j,l) .lt. ymmin ) ymmin = vypos(j,l)
               end do
                     isydf(j) = 0
                  if( ymmax .gt. -1.d+18 ) then
                     isydf(j) = isydf(j) + 1
                     iydef(j,isydf(j)) =  1
                     vypos(j,isydf(j)) = ymmax
                  end if
                  if( ymmin .lt. 1.d+18 ) then
                     isydf(j) = isydf(j) + 1
                     iydef(j,isydf(j)) = -1
                     vypos(j,isydf(j)) = ymmin
                  end if

               if( isydf(j) .eq. 2 .and. iydef(j,1) .eq. -1 ) then
                     vtemp = vypos(j,1)
                     iydef(j,1) =  1
                     vypos(j,1) = vypos(j,2)
                     iydef(j,2) = -1
                     vypos(j,2) = vtemp
               end if

               if( ( issfd(j) .eq. 11 .or.issfd(j) .eq. 14 )
     &               .and. isydf(j) .ne. 2 ) then
                  write(jo,'(/"*** Error, cut surface for ",i7,
     &            " surface source shoud have 2 py cut.")')
     &            issuf(j)
                  ierr = ierr + 1
                  goto 999
               end if

                     zmmin =  1.d+18
                     zmmax = -1.d+18
               do  l = 1, iszdf(j)
                  if( izdef(j,l) .gt. 0 .and.
     &                vzpos(j,l) .gt. zmmax ) zmmax = vzpos(j,l)
                  if( izdef(j,l) .lt. 0 .and.
     &                vzpos(j,l) .lt. zmmin ) zmmin = vzpos(j,l)
               end do
                     iszdf(j) = 0
                  if( zmmax .gt. -1.d+18 ) then
                     iszdf(j) = iszdf(j) + 1
                     izdef(j,iszdf(j)) =  1
                     vzpos(j,iszdf(j)) = zmmax
                  end if
                  if( zmmin .lt. 1.d+18 ) then
                     iszdf(j) = iszdf(j) + 1
                     izdef(j,iszdf(j)) = -1
                     vzpos(j,iszdf(j)) = zmmin
                  end if

               if( iszdf(j) .eq. 2 .and. izdef(j,1) .eq. -1 ) then
                     vtemp = vzpos(j,1)
                     izdef(j,1) =  1
                     vzpos(j,1) = vzpos(j,2)
                     izdef(j,2) = -1
                     vzpos(j,2) = vtemp
               end if

               if( ( issfd(j) .eq. 12 .or.issfd(j) .eq. 15 )
     &               .and. iszdf(j) .ne. 2 ) then
                  write(jo,'(/"*** Error, cut surface for ",i7,
     &            " surface source shoud have 2 pz cut.")')
     &            issuf(j)
                  ierr = ierr + 1
                  goto 999
               end if

*-----------------------------------------------------------------------

            end if

         end if
         end do


*-----------------------------------------------------------------------
*     summary of collision
*-----------------------------------------------------------------------

                  iTYinevt=inevt   !FURUTA20210119
                  iTYiregn=iregn+1 !FURUTA20210119
                  iTYmxmat=mxmat   !FURUTA20210119

*-----------------------------------------------------------------------

      call moddas_deallocate_cha(chrg)
         return


*-----------------------------------------------------------------------

 901     continue
 1901    format('** Error : # of cells exceeds kvlmax = ',i9,
     &        /a, ' Increase kvlmax.')
         write(io,1901) kvlmax
         ErrCha = ''
         MsgID = 'L:17164/R:setpag/F:read00.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,1901) kvlmax
         ierr = ierr + 1
         goto 999

*-----------------------------------------------------------------------

  910    write(io,'(/
     &         "<<< Memory ERROR : at the timer >>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                mmmin, mmmax, mmmax-mmmin, mdas

         ierr = ierr + 1
         goto 999

*-----------------------------------------------------------------------

  920    write(io,'(/
     &         "<<< Memory ERROR : at frag data >>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                mmmin, mmmax, mmmax-mmmin, mdas

         ErrCha = ''
         MsgID = 'L:17198/R:setpag/F:read00.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,'(/
     &         "<<< Memory ERROR : at frag data >>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                mmmin, mmmax, mmmax-mmmin, mdas

         ierr = ierr + 1
         goto 999

*-----------------------------------------------------------------------

  930    write(io,'(/
     &         "<<< Memory ERROR : at the counter field >>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                mmmin, mmmax, mmmax-mmmin, mdas

         ErrCha = ''
         MsgID = 'L:17226/R:setpag/F:read00.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,'(/
     &         "<<< Memory ERROR : at the counter field >>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                mmmin, mmmax, mmmax-mmmin, mdas

         ierr = ierr + 1
         goto 999

*-----------------------------------------------------------------------

  940    write(io,'(/
     &         "<<< Memory ERROR : at summary part of collision "/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                mmmin, mmmax, mmmax-mmmin, mdas

         ErrCha = ''
         MsgID = 'L:17254/R:setpag/F:read00.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,'(/
     &         "<<< Memory ERROR : at summary part of collision "/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                mmmin, mmmax, mmmax-mmmin, mdas

         ierr = ierr + 1
         goto 999

*-----------------------------------------------------------------------

  952    write(io,'(/
     &         "<<< Memory ERROR : at the splitting >>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                mmmin, mmmax, mmmax-mmmin, mdas

         ErrCha = ''
         MsgID = 'L:17282/R:setpag/F:read00.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,'(/
     &         "<<< Memory ERROR : at the splitting >>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                mmmin, mmmax, mmmax-mmmin, mdas

         ierr = ierr + 1
         goto 999

*-----------------------------------------------------------------------

  951    write(io,'(/
     &         "<<< Memory ERROR : at the super mirror >>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                mmmin, mmmax, mmmax-mmmin, mdas

         ErrCha = ''
         MsgID = 'L:17310/R:setpag/F:read00.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,'(/
     &         "<<< Memory ERROR : at the super mirror >>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                mmmin, mmmax, mmmax-mmmin, mdas

         ierr = ierr + 1
         goto 999

*-----------------------------------------------------------------------

  950    write(io,'(/
     &         "<<< Memory ERROR : at the magnetic field >>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                mmmin, mmmax, mmmax-mmmin, mdas

         ErrCha = ''
         MsgID = 'L:17338/R:setpag/F:read00.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,'(/
     &         "<<< Memory ERROR : at the magnetic field >>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                mmmin, mmmax, mmmax-mmmin, mdas

         ierr = ierr + 1
         goto 999

*-----------------------------------------------------------------------

  953    write(io,'(/
     &         "<<< Memory ERROR : at the electro magnetic field >>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                mmmin, mmmax, mmmax-mmmin, mdas

         ErrCha = ''
         MsgID = 'L:17366/R:setpag/F:read00.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,'(/
     &         "<<< Memory ERROR : at the electro magnetic field >>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                mmmin, mmmax, mmmax-mmmin, mdas

         ierr = ierr + 1
         goto 999

*-----------------------------------------------------------------------

  960    write(io,'(/
     &         "<<< Memory ERROR : at ",
     &         "reg = in [source] >>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                mmmin, mmmax, mmmax-mmmin, mdas

         ErrCha = ''
         MsgID = 'L:17395/R:setpag/F:read00.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,'(/
     &         "<<< Memory ERROR : at ",
     &         "reg = in [source] >>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                mmmin, mmmax, mmmax-mmmin, mdas

         ierr = ierr + 1
         goto 999

*-----------------------------------------------------------------------

  971    write(io,'(/
     &         "<<< Memory ERROR : at ",i3,
     &         "-th repeated cell in ",i2,
     &         "-th [repeated collision] section.>>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                m, k, mmmin, mmmax, mmmax-mmmin, mdas

         ErrCha = ''
         MsgID = 'L:17426/R:setpag/F:read00.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,'(/
     &         "<<< Memory ERROR : at ",i3,
     &         "-th repeated cell in ",i2,
     &         "-th [repeated collision] section.>>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                m, k, mmmin, mmmax, mmmax-mmmin, mdas

         ierr = ierr + 1
         goto 999

*-----------------------------------------------------------------------

  970    write(io,'(/
     &         "<<< Memory ERROR : at ",i3,
     &         "-th forced cell in ",i2,
     &         "-th [forced collision] section.>>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                m, k, mmmin, mmmax, mmmax-mmmin, mdas

         ErrCha = ''
         MsgID = 'L:17458/R:setpag/F:read00.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,'(/
     &         "<<< Memory ERROR : at ",i3,
     &         "-th forced cell in ",i2,
     &         "-th [forced collision] section.>>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                m, k, mmmin, mmmax, mmmax-mmmin, mdas

         ierr = ierr + 1
         goto 999

*-----------------------------------------------------------------------

  980    write(io,'(/
     &         "<<< Memory ERROR : at ",i3,
     &         "-th importance cell in ",i2,
     &         "-th [importance] section.>>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                m, k, mmmin, mmmax, mmmax-mmmin, mdas

         ErrCha = ''
         MsgID = 'L:17490/R:setpag/F:read00.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,'(/
     &         "<<< Memory ERROR : at ",i3,
     &         "-th importance cell in ",i2,
     &         "-th [importance] section.>>>"/
     &         "*  memory exceeds mdas",/
     &         "*  start of memory =",i9,/
     &         "*    end of memory =",i9,/
     &         "*     total memory =",i9,/
     &         "*       total mdas =",i9,/
     &         "<<<  Please extend mdas in param.inc >>>")')
     &                m, k, mmmin, mmmax, mmmax-mmmin, mdas

         ierr = ierr + 1
         goto 999

*-----------------------------------------------------------------------

  999 continue
      return
      end

************************************************************************
*                                                                      *
      subroutine echoi(io,jo,ivers,ierr,ierrg)
*                                                                      *
*       echo input parameters                                          *
*       modified by K.Niita on 2010/12/23                              *
*                                                                      *
************************************************************************
      use TDCHAINMOD, only : tdchech
      use COSMICMOD
      use dedx_file
      use ELEDATAMOD, only : ichem, frac, form_chemi
      use moddas
      use moddas_character
      use moddas_material
      use moddas_mesh
      use moddas_multiplier
      use moddas_region
      use moddas_repeated_collisions
      use moddas_source
      use cvaloutmod ! S.H. 2023.10.27
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param01.inc'
      include 'risrcparam.inc'
      include 'err.inc'
*     risrcparam.inc: include parameters
*     integer nuclinmax, nuclpumax, maxchain, chainmax
*     parameter ( nuclinmax =    100 ) nuclides including daughter
*     parameter ( nuclpumax =  1,000 ) Not used
*     parameter ( maxchain  =     23 ) from Number of chain (.NDX)
*     parameter ( chainmax  =  4,050 ) from Number of linear chain (.NDX)
*-----------------------------------------------------------------------

      parameter ( ibmt = 40 )

*-----------------------------------------------------------------------

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      common /ccggg/  icgg
      common /cggmm/  ngstar, ngfini, ngfin0
      common /matmm/  nmhigh, nminth, nmfinh
      common /matmd/  nmlow, nmintl, nmfinl

      common /inggs/  iog, igcel, ioa, igsuf, iob, igtrs
      common /celda/  deng(kvlmax)
      common /celdb/  idsn(kvlmax), idtn(kvlmax)

      common /nnmode/ nmode

      common /inpec/  ititl, ipara, ibody, iregn, llarr, itby, itar
      common /htitl/  iclgt(100), ctitl(100)
      character       ctitl*200

      common /regdu/  iuni(kvlmax)
      common /regde/  ichp(kvlmax), ilat(kvlmax), idct(kvlmax)
      common /regdm/  idmg(kvlmax)
      common /regda/  ichl(kvlmax), chsm(kvlmax), ichmx, iod
      character       chsm*10
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /regdb/  nrsq, irsq(10)
      common /regdl/  ilike(kvlmax)
      common /verjam/ versn, lastr, iyeav, imonv, idayv
      common /startf/ iday0,imon0,iyer0,ihor0,imin0,isec0
      common /parai/  ipsq(400)
      common /paraj/  mstz(300), parz(300)
      common /parak/  icnu, icdf(400), icdl(400), chnm(400)
      character       chnm*8
      common /paral/  lpcn(400), lpcr(400), pcmn(400), pcmr(400)
      character       pcmn*70, pcmr*70  !OBINATA(2012.6.13): change size
      common /paran/  icfn(100), ilfn(100), chfn(100)
      character       chfn*200
      common /userp/  idam(100), rdam(100)
      common /ngcut/  incut, igcut, ipcut

      common /kmatad/ matadd
      common /kmat1a/ mxmat, mxmat0, mxnel
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

      common /kmat1h/ kmatg(kvlmax)

      common /kmat1i/ kmatc(kvlmax)
      common /kmat1j/ intum
      common /kmat1o/ iom1, iom2, iom3
      common /kmat1p/ kmout

      character cnd*1, chd*1

*-----------------------------------------------------------------------

      common /isocor/ iscorr, itcorr, imlwt(isrc)
      common /isomul/ smlwt(isrc), totfact, imsrc
      common /isorfs/ ispfs(isrc), rspfn, rspfz, ispfn
      common /isorst/ jstyp(isrc), istyp(isrc), inkf0(isrc), lstyp(isrc)
!$OMP THREADPRIVATE(/isorst/)
      common /isorori/ jstypori(isrc) ! T.Sato, original jstyp written in input file
      common /isorbia/ isbias(isrc)   ! T.Sato, source bias method for xyz-mesh source
      common /isorsg/ seg0(isrc),seg1(isrc),seg2(isrc),seg3(isrc),
     &                set0(isrc),set1(isrc),set2(isrc),set3(isrc),
     &                jetyp(isrc),jptyp(isrc)
      common /isorsp/ sx0(isrc), sy0(isrc), sz0(isrc), sx1(isrc),
     &                sy1(isrc), sz1(isrc), sr0(isrc), se0(isrc),
     &                sdir(isrc), srx(isrc), sry(isrc), swem(isrc),
     &                sphi(isrc), sdom(isrc), swt0(isrc)
      common /isbeam/ sx2(isrc),sy2(isrc),
     &       sxmrad1(isrc),sxmrad2(isrc),symrad1(isrc),symrad2(isrc)   ! T.Sato, beam emittance source
      common /isorsm/ stm0(isrc),stmw(isrc),stmc(isrc),stmd(isrc),
     &                jttyp(isrc),jttpn(isrc)
      common /isorsn/ sr1(isrc), sr2(isrc), isrn(isrc)
      common /isorsf/ isorf(isrc),lsfile(isrc), sfile(isrc)
      character sfile*100
      common /isorsr/ nsmx(isrc), nsrn(isrc), nsrc(isrc)
      common /isorsc/ isort(isrc,4), rsort(isrc,13)
      common /isorsd/ isdmp(isrc,0:30), jsdmp(isrc,0:30)
      common /isorsa/ sfactor(isrc)
      common /isorpn/ ssx(isrc), ssy(isrc), ssz(isrc)
!$OMP THREADPRIVATE(/isorpn/)
      common /isodct/ sdl0(isrc), sdl1(isrc), sdl2(isrc), sdpf(isrc),
     &                sdxw(isrc), sdyw(isrc), sdrd(isrc), sdebg(isrc),
     &                sdsxp(isrc), sdsxn(isrc), sdsyp(isrc),
     &                sdsyn(isrc), dnorm(isrc), prs1(isrc), prs2(isrc),
     &                prs3(isrc), prs4(isrc), psxp(isrc), psxn(isrc),
     &                psyp(isrc), psyn(isrc), sdls(isrc), sdrs(isrc),
     &                sdxs(isrc), sdys(isrc)

      common /isorse/ ngrp(isrc), ngei(isrc), ngea(isrc), ngfe(isrc),
     &                ngft(isrc), ngll(isrc), ngpi(isrc), ngpw(isrc)

      common /isorfx/ isnm(isrc), lsfx(isrc), srfx(isrc)
      character srfx*200
      common /isousr/ jsusr(isrc,0:30)


      common /isorsl/ nglp(isrc), ngli(isrc), ngla(isrc), ngfl(isrc),
     &                nglc(isrc), nglw(isrc)

*-----------------------------------------------------------------------

      common /isosuf/ issuf(isrc), iscut(isrc), isvct(isrc,8),
     &                issfd(isrc), isvfd(isrc,8),
     &                ivsfd(isrc), ivvfd(isrc,8),
     &                dvsfd(isrc,4), dvvfd(isrc,8,4)

*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma

*-----------------------------------------------------------------------

      common /elmgrg/ nelctf, nereg, inerc, inert, kelcs

      common /magreg/ nmreg, ingrc, ingrt, kmags, magtin, imgusr
      common /tmtreg/ ntmrg, intmc, intmt, ktime
      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)
      common /smireg/ nsreg, isgrc, isgrt, ksmir, ismir

      common /splreg/ isptn, npreg(6), mnspt(6,0:20),
     &                ipgrc(6), ipgrt(6), ksplt(6), isplt(6),
     &                ispct(6,9), ispem(6,2)
      common /splrge/ espem(6,2)

      common /multipl/ imltp,iimlt(multmax),inmlt(multmax),
C MATSUDA 2024.11.25 (lagrange: ilmlt)
     &  ilmlt(multmax),
     &  idmlt(multmax),ismlt(multmax),impan(multmax),impat(multmax,6,3),
     &                 jmpat(multmax,6,6,2), imdfl(multmax)

C MATSUDA 2024.11.25 (multplf: file name)
      common /multplf/ imltf(multmax),lmltfile(multmax),mltfile(multmax)
      character mltfile*100

*-----------------------------------------------------------------------

      common /talmm/  nmmax, lmmax, itlmx
      common /bnkmm/  mbmax, mbfin, mbtfin

      common /tall19/ itsmn(itlmax), itstm(itlmax)

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

      common /tall08/ rtrx0(itlmax), rtry0(itlmax)

      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)
      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)

      common /tall13/ itrcn(itlmax), itrcr(itlmax), itrca(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall37/ itcnt(9,itlmax)
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)

      common /tall82/ itcnth(9,itlmax)

! T.Sato 2024/03/18 weighted history counter ID
      common /tall92/ ichnum(itlmax),chbias(itlmax),pedest(itlmax)
     &               ,ictnum(itlmax),ctbias(itlmax),ictidx(itlmax)
*-----------------------------------------------------------------------

      common /impmsg/ iimpn, isimp, iswct, isstr, maxip,
     &                mnimp(7,0:20),
     &                kfimp(7), inimc(7), inimt(7)

      common /fclmsg/ ifcln, isfcl, maxrg,
     &                mnfcl(6,0:20), mrfcl(kvlmax),
     &                kfcls(6), inflc(6), inflt(6)

*-----------------------------------------------------------------------
      common /rclmsg/ ircln, isrcl, maxrr,
     &                mnrcl(6,0:20), mrrcl(kvlmax),
     &                krcls(6), lrcls(6), inrlc(6), inrlt(6),
     &                irpem(6,2), irman(6), irmct(6), jsmat(6)
      common /rclemm/ erpem(6,2)

      common /isoras/ sag1(isrc),sag2(isrc),jatyp(isrc),jqtyp(isrc)
      common /cosmicint/ipcosmic(isrc),icenv(isrc)
      common /cosmicreal/solarmod(isrc),rigid(isrc),
     &        depatom(isrc),ground(isrc),environ(isrc)
      common /isospg/ spg1(isrc),spg2(isrc) ! T.Sato 2021/02/26

      character chauu*8
      character chaur(7)*8
      character elmnts(104)*3
      character choutall*80

      data elmnts /
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

      common /wwindp/ wupn, wsurvn, mxspln, mwhere, mvoww
      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww

      common /wwind0/ dvals, iwmsh, kecho, kgwwp(6), idval
      common /wwxyz/ iwxty(6), iwxnm(6), iwxrg(6),
     &               rwxmi(6), rwxma(6), rwxdl(6),
     &               iwyty(6), iwynm(6), iwyrg(6),
     &               rwymi(6), rwyma(6), rwydl(6),
     &               iwzty(6), iwznm(6), iwzrg(6),
     &               rwzmi(6), rwzma(6), rwzdl(6)
      common /wwtrs/ iwwtr(6,4), rwwtr(6,13)
      common /wbxyz/ iwbxty(6), iwbxnm(6), iwbxrg(6),
     &               rwbxmi(6), rwbxma(6), rwbxdl(6),
     &               iwbyty(6), iwbynm(6), iwbyrg(6),
     &               rwbymi(6), rwbyma(6), rwbydl(6),
     &               iwbzty(6), iwbznm(6), iwbzrg(6),
     &               rwbzmi(6), rwbzma(6), rwbzdl(6)
      common /wbtrs/ iwbtr(6,4), rwbtr(6,13)

      common /delesg/ rdels(kvlmax), mndel, ndels(kvlmax)
      common /delreg/ delm(kvlmax), kdelt

! T.Sato for Track structure simulation
      common /tscmsg/ ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /etsart/ bgets(kvlmax), ebgets(kvlmax),
     &                wvets(kvlmax), ewvets(kvlmax)

      common /volmsg/ rvols(kvlmax), mnvol, nvols(kvlmax)
      common /regdd/  ivolm, iimpo
      common /regim/  iimps
      common /regdn/  irden
      common /tmpmsg/ rtmps(kvlmax), mntmp, ntmps(kvlmax)
      common /brsmsg/ mnbrs, icbrs, mbbrs(kvlmax), cbrem(49)
      common /pwtmsg/ rpwts(kvlmax), mnpwt, npwts(kvlmax)
      common /pwtreg/ dpwt(kvlmax)

      common /impreg/ dimp(kvlmax)
      common /volreg/ dvol(kvlmax)

      common /ptname/ pname(20), ipln(20)
      character       pname*8

      common /mtnmcl/ dmhsb(-1:kvlmax,4),
     &                dmtnm(-1:kvlmax), dmtcl(-1:kvlmax),
     &                nmtnm(-1:kvlmax), nmtcl(-1:kvlmax)
      character dmtnm*80, dmtcl*30
      common /mttmc/  smttc(kvlmax), mttcn, mttc1(kvlmax), mttc2(kvlmax)
      common /mtnmc/  smtnc(kvlmax), dmtnc(kvlmax,2),
     &                mtncn, mtnc(kvlmax,2), nmtnc(kvlmax,2)
      character dmtnc*80

      common /mtreg/  smtrg(kvlmax), dmtrg(kvlmax),
     &                mtrgn, mtrg(kvlmax,2), nmtrg(kvlmax)
      character dmtrg*80

      common /mtnreg/ dmhsg(-1:kvlmax),
     &                nmtng(-1:kvlmax), dmtng(-1:kvlmax)
      character dmtng*80

      common /elreg/  elarg(kvlmax,4), mlrgn, melrg(kvlmax,2), ielusr

      common /fgdata/ ifrgd, ifgm, ifgdf(kvlmax,5), frgfl(kvlmax)
      character frgfl*200

      character nctemp*200
      character nrtem1*7
      character nrtem2*7

      common /dumpall/ idumpall
      common /voxel/  ivoxel
      common /pnint/  ipnint

      common /ggcell/ icells, iobo

      common /infprint/ infout

      common /risource/ chalct(nuclinmax,isrc), act000(nuclinmax,isrc),
     &                  activy(nuclinmax,isrc), thf(nuclinmax,isrc),
     &                  decayt(isrc)
      real  chalct
      double precision  act000, activy, thf, decayt

C MATSUDA 2017.05.29 (iaugers)
C MATSUDA 2018.08.15 (icharacterx)
C MATSUDA 2019.05.08 (iannih)
      common /risrc00/ niorg(isrc), nicur(isrc), norm(isrc),
     &                 iaugers(isrc), icharacterx(isrc), iannih(isrc),
     &                 aclow(isrc)
      integer  niorg, nicur, norm, iaugers, icharacterx, iannih
      double precision  aclow
      common /risrc01/ normfact(isrc), asfsum(isrc)
      double precision normfact, asfsum
      common /risrc02/ smlwt2(isrc), totfact2
      double precision karival, smlwt2

cfrtati 2021/12/17 added [data max] for photo-nuclear, deuteron, alpha
      common /nntmax/ dmxdxx(6,500),indmm,ipdmm(6,6),ipdnn(6),ipdpt(6),
     &                nucdxx(6,500),matdxx(6,500)
      common /ndtmax/ indmp, indmn, indmu, indmd, indma,
     & nucdxp(500), nucdxn(500), nucdxu(500), nucdxd(500), nucdxa(500),
     & matdxp(500), matdxn(500), matdxu(500), matdxd(500), matdxa(500)
      common /ddtmax/ dmxdxp(500), dmxdxn(500),
     &                dmxdxu(500), dmxdxd(500), dmxdxa(500)

      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax)

*-----------------------------------------------------------------------
      common /wwbiasn/ eenwb(6,100), iwbdp, mnwbp(6,0:20), kfwbp(6),
     &                 inwbc(6), ienwb(6), inwbt(6), iswbp, maxwb
      common /nwwbias/ iwwbias

*-----------------------------------------------------------------------

      common /isormd/ istdd(isrc),istcc(isrc),
     &                isxtp(isrc),isinx(isrc),istxx(isrc),
     &                isytp(isrc),isiny(isrc),istyy(isrc),
     &                isztp(isrc),isinz(isrc),istzz(isrc),
     &                sxmin(isrc),sxmax(isrc),sxdel(isrc),
     &                symin(isrc),symax(isrc),sydel(isrc),
     &                szmin(isrc),szmax(isrc),szdel(isrc)

*-----------------------------------------------------------------------
      common /cntech/ icnech(3,0:maxcntr+1) ! S.H. set maxcntr+'non' (2022.3.24)

*-----------------------------------------------------------------------
      common /ibchsor/ ibcsn(isrc), ibsor(isrc,isrc)
      common /stat / istdev, irestart, ireschk
      character fmtt*20

*-----------------------------------------------------------------------

      character dum1*10000
      character dum2*10000
      character dum3*10000

      character asfil*100

      logical   exex

      character chtit*60

      dimension ipva(4)
      dimension bval(50)
      character chbd*3
      character chau*8
      character chat*8
      character chav*8

      character chin*200
      character jobtl*13
      data jobtl /'[ Job Title ]'/
      character cblan*200
      character cmins*200

      character chdf*80

      logical deqn1
      logical deqn5

      dimension ix(3)
      character ht*10
      character hs*80

*-----------------------------------------------------------------------

      character sname*9
      character rname*51
      character qname*16

      character uname(20)*8

      dimension jmat(20)
      dimension rmat(20)

*-----------------------------------------------------------------------

      dimension ibck(kvmmax)
      data ibck / kvmmax*0 /

      dimension     idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )



*-----------------------------------------------------------------------
      common /canat/  ianat

*-----------------------------------------------------------------------
      character chaus*8
      character chaum*6

      character       elmnt(104)*2

      data elmnt /
     &    ' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
     &    'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca',
     &    'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     &    'Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr',
     &    'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     &    'Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd',
     &    'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     &    'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     &    'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     &    'Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     &    'Md','No','Lr','Ku' /

*-----------------------------------------------------------------------

      character chsf(50)*3

      data ( chsf(i), i = 1, ibmt ) /
     &         'p  ','px ','py ','pz ','so ',
     &         's  ','sx ','sy ','sz ','c/x',
     &         'c/y','c/z','cx ','cy ','cz ',
     &         'k/x','k/y','k/z','kx ','ky ',
     &         'kz ','sq ','gq ','tx ','ty ',
     &         'tz ','x  ','y  ','z  ','box',
     &         'rpp','sph','rcc','rec','ell',
     &         'trc','wed','arb','rhp','hex'/

      character chss*6

      character chsn*7
      character chsc*1

      character chtr*13

      dimension vtrs(13)

*-----------------------------------------------------------------------

      character dmpc(20)*4
      data dmpc / '  kf','   x','   y','   z','   u','   v','   w',
     &            '   e','  wt','  tm','  c1','  c2','  c3',
     &            '  sx','  sy','  sz','    ','    ','    ','    '/

*-----------------------------------------------------------------------
      common /stat2/ dmpmulti,idmpmode,ibchjmp,idmpjmp(2)
      common /stat3/ jpsf

*-----------------------------------------------------------------------
      common /itetsor/ ksoutnode,itetreg(isrc)

*-----------------------------------------------------------------------
      character chemID*200
      character ckam*3

*-----------------------------------------------------------------------
      character yen*1
      yen = char(92)

*-----------------------------------------------------------------------

            lib0 = ichar(' ') + 256 * ( ichar(' ')
     &                        + 256 * ichar(' ') )

*-----------------------------------------------------------------------

            ierr = ierrg

         do i = 1, 200

            cblan(i:i) = ' '
            cmins(i:i) = '_'

         end do

*-----------------------------------------------------------------------
*    open temporary file io = 29
*-----------------------------------------------------------------------

               ioi = 29

               open(ioi,form='formatted',status='scratch')

               iot = io

*-----------------------------------------------------------------------
*     header
*-----------------------------------------------------------------------

      if( ivers .ne. 0 .and. ierrg .eq. 0 ) then

      write(iot,'(10x," ",57("_")/10x,"|",57x,"|")')

      write(iot,'(10x,
     &"|       _/_/_/_/                                         ",
     &" |"/10x,
     &"|      _/      _/                  _/_/_/_/_/_/          ",
     &" |"/10x,
     &"|     _/      _/  _/      _/   _/      _/      _/_/_/_/_/",
     &" |"/10x,
     &"|    _/_/_/_/    _/      _/   _/      _/      _/         ",
     &" |"/10x,
     &"|   _/          _/_/_/_/_/   _/      _/       _/_/_/_/   ",
     &" |"/10x,
     &"|  _/          _/      _/   _/      _/              _/   ",
     &" |"/10x,
     &"| _/          _/      _/   _/      _/      _/_/_/_/_/    ",
     &" |")
     &')

      write(iot,'(10x,
     &"|                                                        ",
     &" |"/10x,
     &"|                                                        ",
     &" |"/10x,
     &"|       Particle and Heavy Ion Transport code System     ",
     &" |"/10x,
     &"|                      Version =",
     &                                  f7.3,"                  ",
     &" |")') versn

      if( infout .eq. 7 .or. infout .eq. 8 ) then

      write(iot,'(10x,
     &"|                       developed by                      |"/10x,
     &"|                                                         |"/10x,
     &"|  Tatsuhiko SATO, Yosuke IWAMOTO, Shintaro HASHIMOTO,    |"/10x,
     &"|    Tatsuhiko OGAWA, Takuya FURUTA, Shinichiro ABE,      |"/10x,
     &"|    Takeshi KAI, Norihiro MATSUDA, Yusuke MATSUYA,       |"/10x,
     &"|      Yuho HIRATA, Takuya SEKIKAWA, Pi-En TSAI,          |"/10x,
     &"|                Hunter RATLIFF (JAEA),                   |"/10x,
     &"|                                                         |"/10x,
     &"|  Hiroshi IWASE, Yasuhito SAKAKI, Kenta SUGIHARA (KEK),  |"/10x,
     &"|                                                         |"/10x,
     &"|           Nobuhiro SHIGYO (Kyushu University),          |"/10x,
     &"|                                                         |"/10x,
     &"|      Lembit SIHVER (Technische Universitat Wien), and   |"/10x,
     &"|                                                         |"/10x,
     &"|                     Koji NIITA (RIST)                   |")')
      end if

      write(iot,'(10x,
     &"|                                                        ",
     &" |"/10x,
     &"|                 Last Revised  ",
     &                 i4,"-",i2.2,"-",i2.2,"               ",
     &" |"
     &                  )') iyeav,imonv, idayv

      write(iot,'(10x,"|",57("_"),"|")')

      if( npe .gt. 0 ) then

      write(iot,'(10x,
     &"|                                                        ",
     &" |"/10x,
     &"|             This is a Parallel Version by MPI          ",
     &" |"/10x,
     &"|         Number of Executable PE is ( Total - 1 )       ",
     &" |"/10x,
     &"|                                                        ",
     &" |"/10x,
     &"|                    Total PE =",
     &                                i4,"                      ",
     &" |"
     &                  )') npe


      write(iot,'(10x,"|",57("_"),"|")')

      end if

      end if

*-----------------------------------------------------------------------
*     old version to new version
*-----------------------------------------------------------------------

      if( ivers .eq. 0 ) then

      write(iot,'(10x," ",57("_")/10x,"|",57x,"|")')

      write(iot,'(10x,
     &"|                   Congratulations !!                   ",
     &" |"/10x,
     &"|                                                        ",
     &" |"/10x,
     &"|            Your input file is successfully             ",
     &" |"/10x,
     &"|            transformed to new input file below.        ",
     &" |"
     &                  )')

      write(iot,'(10x,
     &"|                                                        ",
     &" |"/10x,
     &"|",13x,"Created Date = ",
     &                      i4,"-",i2.2,"-",i2.2,19x,"|"/10x,
     &"|",13x,"        Time = ",
     &                      i2.2,"h ",i2.2,"m ",i2.2,19x,"|")')
     &                      iyer0,imon0,iday0,
     &                      ihor0,imin0,isec0

      write(iot,'(10x,
     &"|                                                        ",
     &" |"/10x,
     &"|      Please replace below < fort.***  >                ",
     &" |"/10x,
     &"|      by the file name with the path name,              ",
     &" |"/10x,
     &"|      which is linked to < fort.*** > in your shell     ",
     &" |"
     &                  )')

      write(iot,'(10x,"|",57("_"),"|")')

      end if

*-----------------------------------------------------------------------
*     job title
*-----------------------------------------------------------------------

      if( ititl .gt. 0 .and. ivers .ne. 0 .and. icntl .ne. 3 .and.
     &    ierrg .eq. 0 ) then

            lmaxc = 0

         do i = 1, ititl

            if( iclgt(i) .gt. lmaxc ) lmaxc = iclgt(i)

         end do

            lmaxa = max( lmaxc, 14 )

            if( lmaxa + 6 .gt. 74 ) then
               istr = 1
            else
               istr = max( 1, 34 - lmaxa / 2 )
            end if

            lmins = max( 1, ( lmaxc - 13 ) / 2 + 2 )
            lminl = max( 1, lmaxa - lmins - 13 + 2 )

            write(iot,'(/200a1)')
     &         (cblan(j:j),j=1,istr),      ' ',
     &         (cmins(k:k),k=1,lmins),
     &         ' ',(jobtl(j:j),j=1,13),' ',
     &         (cmins(k:k),k=1,lminl),     ' '

            write(iot,'(200a1)')
     &         (cblan(j:j),j=1,istr),      '|',
     &         (cblan(k:k),k=1,lmaxa+4),   '|'

         do i = 1, ititl

            ienr = lmaxa - iclgt(i) + 2

            write(iot,'(200a1)')
     &         (cblan(j:j),j=1,istr),      '|',' ',' ',
     &         (ctitl(i)(k:k),k=1,iclgt(i)),
     &         (cblan(j:j),j=1,ienr),      '|'

         end do

            write(iot,'(200a1)')
     &         (cblan(j:j),j=1,istr),      '|',
     &         (cmins(k:k),k=1,lmaxa+4),   '|'

      end if

*-----------------------------------------------------------------------
*     starting date
*-----------------------------------------------------------------------

      if( ivers .ne. 0 .and. icntl .ne. 3 .and. ierrg .eq. 0 ) then

         write(iot,'(/24x,"Starting Date = ",
     &                      i4,"-",i2.2,"-",i2.2)')
     &                      iyer0,imon0,iday0
         write(iot,'( 24x,"Starting Time = ",
     &                      i2.2,"h ",i2.2,"m ",i2.2)')
     &                      ihor0,imin0,isec0

      end if

*-----------------------------------------------------------------------
*     input echo or print geometry error
*-----------------------------------------------------------------------


      if( infout .eq. 1 .or. infout .eq. 4 .or. infout .eq. 5 .or.
     &    infout .eq. 7 .or. infout .eq. 8 .or. ierrg. ne. 0 ) then
         iot = io
      else
         iot = ioi
      end if

*-----------------------------------------------------------------------

      if( ivers .ne. 0 .and. icntl .ne. 3 .and. ierrg .eq. 0 ) then

         write(iot,'(/">>> Input Echo >>>",61("="))')

      else if( ivers .ne. 0 .and. icntl .ne. 3 .and. ierrg .ne. 0 ) then

         write(io,'(/">>> Print Geometry Errors >>>")')
         ErrCha = ''
         MsgID = 'L:18288/R:echoi/F:read00.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,'(/">>> Print Geometry Errors >>>")')

      end if

*-----------------------------------------------------------------------
*     title
*-----------------------------------------------------------------------

      if( ititl .gt. 0 .and. ierrg .eq. 0 ) then

         write(iot,'(/"[ Title ]")')

         do i = 1, ititl

            write(iot,'(200a1)') (ctitl(i)(k:k),k=1,iclgt(i))

         end do

      end if

*-----------------------------------------------------------------------
*     parameters
*-----------------------------------------------------------------------

      if( ierrg .eq. 0 ) then

*-----------------------------------------------------------------------

         write(iot,'(/"[ Parameters ]")')

         if( ipara .gt. 0 ) then

            do k = 1, ipara

                  j = ipsq(k)

               do i = 1, icnu

                  if( j .eq. icdf(i) ) goto 100

               end do

*-----------------------------------------------------------------------

               if( abs( j ) .lt. 1000 ) then

                  if( j .gt. 0 ) then
                     write(iot,'(1x,"mstz(",i3,") =",i11)')
     &               j, mstz(j)

                  else
                     write(iot,'(1x,"parz(",i3,") =",1p1g16.9)')
     &               -j, parz(-j)

                  end if

               else

                  if( j .gt. 0 ) then

                     jj = j - 1000
                     write(iot,'(1x,"idam(",i3,") =",i11)')
     &               jj, idam(jj)

                  else

                     jj = -j - 1000
                     write(iot,'(1x,"rdam(",i3,") =",1p1g16.9)')
     &               jj, rdam(jj)

                  end if

               end if

               goto 400

*-----------------------------------------------------------------------

  100          continue
!OBINATA(2012.6.13): formats are modified for change size of pcmn, pcmr.
               if( j .gt. 0 ) then

                if(i.ne.326 .and. (i.lt.351.or.i.gt.357)) then ! frtati 2021/12/17 for lib

                     write(iot,'(1x,a8," =",i12,4x," # ",70a1)')
     &               chnm(i), mstz(j), (pcmn(j)(l:l),l=1,lpcn(j))

                else if(i.ge.351.and.i.le.357) then

                     ckam = char(mod(mstz(j),256))//
     &               char(mod(mstz(j)/256,256))//char(mstz(j)/65536)
                     write(iot,'(1x,a8," =",a12,4x," # ",70a1)')
     &               chnm(i), ckam, (pcmn(j)(l:l),l=1,lpcn(j))

                else
                     write(iot,'("$ ",a8," =",i12,4x," # ",70a1)')
     &               chnm(i), mstz(j), (pcmn(j)(l:l),l=1,lpcn(j))

                endif

               else

c S.H. xorshift (2020.2.6)
                  if( -j .eq. 195 ) then

                  write(iot,'(1x,a8," =",1x,b64.64,/28x,"# ",70a1)')
     &               chnm(i), parz(-j), (pcmr(-j)(l:l),l=1,lpcr(-j))

                  else if( -j .eq. 21 ) then

                  write(iot,'(1x,a8," =",1p1e25.16e3,4x,"# ",70a1)')
     &               chnm(i), parz(-j), (pcmr(-j)(l:l),l=1,lpcr(-j))

                  else if( parz(-j) .eq. 0.0 ) then

                     write(iot,'(1x,a8," =",f12.1,4x," # ",70a1)')
     &               chnm(i), parz(-j), (pcmr(-j)(l:l),l=1,lpcr(-j))

                  else

                     write(iot,'(1x,a8," =",1p1g16.9," # ",70a1)')
     &               chnm(i), parz(-j), (pcmr(-j)(l:l),l=1,lpcr(-j))

                  end if

               end if

  400          continue

            end do

         end if

*-----------------------------------------------------------------------
*        file name
*-----------------------------------------------------------------------

               if( icntl .ne. 2 )  icfn( 2) = 0
               if( icntl .ne. 2 )  icfn( 3) = 0
               if( icntl .ne. 4 )  icfn( 4) = 0
               if( ipcut .ne. 0 )  icfn(10) = 1
               if( incut .ne. 0 )  icfn(12) = 1
               if( igcut .ne. 0 )  icfn(13) = 1
               icfn(14) = 0
               if( ivers .eq. 0 )  icfn( 6) = 1

               if( idumpall .ne. 0 ) icfn(15) = 1

               if( ivoxel .ne. 0 .and. ivoxel .ne. 3 ) icfn(18) = 1
               if( icells .ne. 0 .and. icells .ne. 3 ) icfn(19) = 1


         do k = 1, 100

            if( icfn(k) .ne. 0 ) then

                  ilbl = max( 1, 16 - ilfn(k) )

               if( k .eq. 1 ) then ! T.Sato 2016/06/08

                  chin = '# (D=PhitsPath) PHITS install folder name'
                  ilcm = 41

               elseif( k .eq. 2 ) then

                  chin = '# (D=cgview.in) cgview input file name'
                  ilcm = 38

               else if( k .eq. 3 ) then

                  chin = '# (D=cgview.set) cgview setting file name'
                  ilcm = 41

               else if( k .eq. 4 ) then

                  chin = '# (D=marspf.in) marspf input file name'
                  ilcm = 38

               else if( k .eq. 6 ) then

                  chin = '# (D=phits.out) general output file name'
                  ilcm = 40

               else if( k .eq. 7 ) then

      chin='# (D=c:/phits/data/xsdir.jnd) nuclear data input file name'
                  ilcm = 58

               else if( k .eq. 11 ) then

                  chin = '# (D=compas.out) COMPAS output file name'
                  ilcm = 40

               else if( k .eq. 10 ) then

                  chin = '# (D=fort.10) cutoff proton file name'
                  ilcm = 37

               else if( k .eq. 12 ) then

                  chin = '# (D=fort.12) cutoff neutron file name'
                  ilcm = 38

               else if( k .eq. 13 ) then

                  chin = '# (D=fort.13) cutoff photon file name'
                  ilcm = 37



               else if( k .eq. 15 ) then

                  chin = '# (D=dumpall.dat) dumpall file name'
                  ilcm = 35

               else if( k .eq. 18 ) then

                  chin = '# (D=voxel.bin) binary voxel data'
                  ilcm = 33

               else if( k .eq. 19 ) then

                  chin = '# (D=gcell.bin) binary gcell data'
                  ilcm = 33

               else if( k .eq. 20 ) then

           chin='# (D=c:/phits/XS/egs) EGS library data folder name'
            ilcm = 50

               else if( k .eq. 21 ) then

           chin='# (D=c:/phits/dchain-sp/data) dchain data folder name'
            ilcm = 53

               else if( k .eq. 24 ) then

                  chin='# (D=c:/phits/data) DECDC2 data folder name'
                  ilcm = 43

               else if( k .eq. 25 ) then
           chin='# (D=c:/phits/XS/tra) TS data folder name'
                  ilcm = 41

               else if( k .eq. 26 ) then

           chin='# (D=c:/phits/data/multiplier) multiplier folder name'
                  ilcm = 53

               else if( k .eq. 27 ) then

           chin='# (D=c:/phits/XS/yield) yield data files folder name'
                  ilcm = 52

               else if( k .eq. 28 ) then

           chin='# (D=c:/phits/data/aama.dat) aama parameter file name'
                  ilcm = 53

               else if( k .eq. 29 ) then
           chin='# (D=c:/phits/data/dedx) dedx folder name'
                  ilcm = 40

               else if( k .eq. 30 ) then
                  chin = '# (D=tetra.bin) binary tetra data'
                  ilcm = 33

               end if

               if( k .lt. 10 ) then

                  write(iot,'(1x,"file(",i1,")  = ",200a1)')
     &               k, (chfn(k)(l:l),l=1,ilfn(k)),
     &                  (cblan(l:l), l=1,ilbl),
     &                  (chin(l:l), l=1,ilcm)

               else

                  write(iot,'(1x,"file(",i2,") = ",200a1)')
     &               k, (chfn(k)(l:l),l=1,ilfn(k)),
     &                  (cblan(l:l), l=1,ilbl),
     &                  (chin(l:l), l=1,ilcm)

               end if

            end if

         end do

*-----------------------------------------------------------------------

      if( mstz(45) .gt. 0 ) then

               write(iot,'(/"#------- default values of parameters",
     &         " ----------------------------------------#")')

            do i = 1, icnu

               do k = 1, ipara

                  j = ipsq(k)

                  if( j .eq. icdf(i) ) goto 410

               end do

                  j = icdf(i)

               if( j .gt. 0 ) then

                     write(iot,'(1x,a8," =",i12,4x," # ",70a1)')
     &               chnm(i), mstz(j), (pcmn(j)(l:l),l=1,lpcn(j))

               else

c S.H. xorshift (2020.2.6)
                  if( -j .eq. 195 ) then

                  write(iot,'(1x,a8," =",1x,b64.64,/28x,"# ",70a1)')
     &               chnm(i), parz(-j), (pcmr(-j)(l:l),l=1,lpcr(-j))

                  else if( -j .eq. 21 ) then

                  write(iot,'(1x,a8," =",1p1e25.16e3,4x,"# ",70a1)')
     &               chnm(i), parz(-j), (pcmr(-j)(l:l),l=1,lpcr(-j))

                  else if( parz(-j) .eq. 0.0 ) then

                     write(iot,'(1x,a8," =",f12.1,4x," # ",70a1)')
     &               chnm(i), parz(-j), (pcmr(-j)(l:l),l=1,lpcr(-j))

                  else

                     write(iot,'(1x,a8," =",1p1g16.9," # ",70a1)')
     &               chnm(i), parz(-j), (pcmr(-j)(l:l),l=1,lpcr(-j))

                  end if

               end if

  410          continue

            end do

*-----------------------------------------------------------------------

               write(iot,'(/"#------- default filename -----------",
     &         "-----------------------------------------#")')

         do k = 1, 100

            if( icfn(k) .eq. 0 ) then

                  ilbl = max( 1, 16 - ilfn(k) )
                  ilcm = 0

               if( k .eq. 2 ) then

                  chin = '# (D=cgview.in) cgview input file name'
                  ilcm = 38

               else if( k .eq. 3 ) then

                  chin = '# (D=cgview.set) cgview setting file name'
                  ilcm = 41

               else if( k .eq. 4 ) then

                  chin = '# (D=marspf.in) marspf input file name'
                  ilcm = 38

               else if( k .eq. 6 ) then

                  chin = '# (D=phits.out) general output file name'
                  ilcm = 40

               else if( k .eq. 7 ) then

                  chin = '# (D=xdirs) nuclear data input file name'
                  ilcm = 40

               else if( k .eq. 11 ) then

                  chin = '# (D=nuclcal.out) n-reac. output file name'
                  ilcm = 42

               else if( k .eq. 10 ) then

                  chin = '# (D=fort.10) cutoff proton file name'
                  ilcm = 37

               else if( k .eq. 12 ) then

                  chin = '# (D=fort.12) cutoff neutron file name'
                  ilcm = 38

               else if( k .eq. 13 ) then

                  chin = '# (D=fort.13) cutoff photon file name'
                  ilcm = 36



               else if( k .eq. 15 ) then

                  chin = '# (D=dumpall.dat) dump all data file name'
                  ilcm = 40

               else if( k .eq. 18 ) then

                  chin = '# (D=voxel.bin) voxel binary data file name'
                  ilcm = 44

               else if( k .eq. 21 ) then

           chin='# (D=c:/phits/dchain-sp/data) dchain data folder name'
            ilcm = 53

               else if( k .eq. 24 ) then

                  chin='# (D=c:/phits/data) DECDC2 data folder name'
                  ilcm = 43

               else if( k .eq. 26 ) then

           chin='# (D=c:/phits/data/multiplier) multiplier folder name'
                  ilcm = 53

               else if( k .eq. 27 ) then

           chin='# (D=c:/phits/XS/yield) yield data files folder name'
                  ilcm = 52

               else if( k .eq. 30 ) then

                  chin = '# (D=tetra.bin) tetra binary data file name'
                  ilcm = 44

               end if

               if( ilcm .gt. 0 ) then

                  if( k .lt. 10 ) then

                     write(iot,'(1x,"file(",i1,")  = ",200a1)')
     &                  k, (chfn(k)(l:l),l=1,ilfn(k)),
     &                     (cblan(l:l), l=1,ilbl),
     &                     (chin(l:l), l=1,ilcm)

                  else

                     write(iot,'(1x,"file(",i2,") = ",200a1)')
     &                  k, (chfn(k)(l:l),l=1,ilfn(k)),
     &                     (cblan(l:l), l=1,ilbl),
     &                     (chin(l:l), l=1,ilcm)

                  end if

               end if

            end if

         end do

      end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*     source
*-----------------------------------------------------------------------

      if( ierrg .eq. 0 ) then

         write(iot,'(/"[ Source ]")')

          do j = 1, imsrc

            if( jetyp(j) .ne. 28 .and. jetyp(j) .ne. 29 .and.
     &          jetyp(j) .ne. 25 .and. jetyp(j) .ne. 26) then  ! T.Sato 2020/12/06
              if( j .eq. imsrc )
     &        write(iot,'("  totfact = ",1p1g12.5,3x,
     &        " # (D=1.0) global factor")') totfact
            else
              write(iot,'("  totfact = ",1p1g12.5,3x,
     &        " # (D=1.0) global factor")') totfact2
              write(iot,'("$ totfact = ",1p1g12.5,3x,
     &        " # revised global factor by e-type = 25,26,28 or 29")')
     &        totfact
                exit
            end if

          end do


            if( iscorr .ne. 0 ) then

               write(iot,'("   iscorr = ",i3,12x,
     &         " # (D=0) multi correlation source")') iscorr

               write(iot,'("            ",3x,12x,
     &         " # total number of corr. source = ",i3)') itcorr

            end if

      do j = 1, imsrc

         if( imsrc .gt. 1 .and. iscorr .eq. 0 ) then

            if( jetyp(j) .ne. 28 .and. jetyp(j) .ne. 29 .and.
     &          jetyp(j) .ne. 25 .and. jetyp(j) .ne. 26) then  ! T.Sato 2020/12/06
               write(iot,'(/" <Source> = ",1p1g12.5,3x,
     &         " # weight of this sub-source")') smlwt(j)
           else
               write(iot,'(/" <Source> = ",1p1g12.5,3x,
     &         " # weight of this sub-source")') smlwt2(j)
               write(iot,'("$<Source> = ",1p1g12.5,3x,
     &         " # revised <source> by e-type=25,26,28 or 29")')
     &         smlwt(j)
           end if

         else if( imsrc .gt. 1 .and. iscorr .ne. 0 ) then

               write(iot,'(/" <Source> = ",i3,12x,
     &         " # number of this sub-source")') imlwt(j)

         end if

*-----------------------------------------------------------------------

            if( istyp(j) .eq. 19 ) then

               iz = ichgf(istyp(j),inkf0(j))
               ia = ibryf(istyp(j),inkf0(j))

               call chname(idum,ia,iz,chau)

               rname = chau(1:8)
     &                 //'       # kind of incident nucleus'

            else if( istyp(j) .ne. 11 .and. istyp(j) .gt. 0 ) then

               rname = pname(istyp(j))(1:8)
     &                 //'       # kind of incident particle'

            else if( istyp(j) .eq. 11 ) then

               call jamname(inkf0(j),0,0,qname)
               call kfcname(inkf0(j),9,sname)

               rname = sname(1:9)
     &                 //'      # '//qname(1:8)
     &                 //' : kind of incident particle'

            end if

*-----------------------------------------------------------------------

         if( jstyp(j) .eq. 1 ) then

            write(iot,'("   s-type = ",i3,12x,
     &      " # cylindrical source")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if

          if( ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:# of SF, 2:neutron" )')
     &         ispfs(j)

          end if

            write(iot,'("       e0 = ",1p1g12.5,3x,
     &      " # energy of beam [MeV/n]")') se0(j)


            write(iot,'("       r0 = ",1p1g12.5,3x,
     &      " # radius [cm]")') sr0(j)

            if( sr0(j) .gt. sr1(j) )
     &      write(iot,'("       r1 = ",1p1g12.5,3x,
     &      " # minimum radius [cm]")') sr1(j)

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of x-axis [cm]")') sx0(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of y-axis [cm]")') sy0(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # minimum position of z-axis [cm]")') sz0(j)

            write(iot,'("       z1 = ",1p1g12.5,3x,
     &      " # maximum position of z-axis [cm]")') sz1(j)

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

C S.H. added all option for phi (2018.10.2)
            else if( dabs(sphi(j)+1000.0) .lt. 1d-9 ) then

               write(iot,'("      phi =   all",10x,
     &         " # azimuthal angle of beam [random]")')

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if


            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 4 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # cylindrical source")')
     &            jstypori(j)  ! T.Sato 2017/07/04

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            write(iot,'("       r0 = ",1p1g12.5,3x,
     &      " # radius [cm]")') sr0(j)

            if( sr0(j) .gt. sr1(j) )
     &      write(iot,'("       r1 = ",1p1g12.5,3x,
     &      " # minimum radius [cm]")') sr1(j)

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of x-axis [cm]")') sx0(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of y-axis [cm]")') sy0(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # minimum position of z-axis [cm]")') sz0(j)

            write(iot,'("       z1 = ",1p1g12.5,3x,
     &      " # maximum position of z-axis [cm]")') sz1(j)

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 2 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # rectangular-solid source")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            write(iot,'("       e0 = ",1p1g12.5,3x,
     &      " # energy of beam [MeV/n]")') se0(j)

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # minimum position of x-axis [cm]")') sx0(j)

            write(iot,'("       x1 = ",1p1g12.5,3x,
     &      " # maximum position of x-axis [cm]")') sx1(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # minimum position of y-axis [cm]")') sy0(j)

            write(iot,'("       y1 = ",1p1g12.5,3x,
     &      " # maximum position of y-axis [cm]")') sy1(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # minimum position of z-axis [cm]")') sz0(j)

            write(iot,'("       z1 = ",1p1g12.5,3x,
     &      " # maximum position of z-axis [cm]")') sz1(j)

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 5 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # rectangular-solid source")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # minimum position of x-axis [cm]")') sx0(j)

            write(iot,'("       x1 = ",1p1g12.5,3x,
     &      " # maximum position of x-axis [cm]")') sx1(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # minimum position of y-axis [cm]")') sy0(j)

            write(iot,'("       y1 = ",1p1g12.5,3x,
     &      " # maximum position of y-axis [cm]")') sy1(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # minimum position of z-axis [cm]")') sz0(j)

            write(iot,'("       z1 = ",1p1g12.5,3x,
     &      " # maximum position of z-axis [cm]")') sz1(j)

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 3 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # Gaussian distribution source")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            write(iot,'("       e0 = ",1p1g12.5,3x,
     &      " # energy of beam [MeV/n]")') se0(j)

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of x-axis [cm]")') sx0(j)

            write(iot,'("       x1 = ",1p1g12.5,3x,
     &      " # Full Width at Half Maximum of x-axis [cm]")') sx1(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of y-axis [cm]")') sy0(j)

            write(iot,'("       y1 = ",1p1g12.5,3x,
     &      " # Full Width at Half Maximum of y-axis [cm]")') sy1(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # center position of z-axis [cm]")') sz0(j)

            write(iot,'("       z1 = ",1p1g12.5,3x,
     &      " # Full Width at Half Maximum of z-axis [cm]")') sz1(j)

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 6 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # Gaussian source")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of x-axis [cm]")') sx0(j)

            write(iot,'("       x1 = ",1p1g12.5,3x,
     &      " # Full Width at Half Maximum of x-axis [cm]")') sx1(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of y-axis [cm]")') sy0(j)

            write(iot,'("       y1 = ",1p1g12.5,3x,
     &      " # Full Width at Half Maximum of y-axis [cm]")') sy1(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # center position of z-axis [cm]")') sz0(j)

            write(iot,'("       z1 = ",1p1g12.5,3x,
     &      " # Full Width at Half Maximum of z-axis [cm]")') sz1(j)

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 7 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # parabola source")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            write(iot,'("       rn = ",i3,12x,
     &      " # (D=2) pawer of parabola")') isrn(j)

            write(iot,'("       e0 = ",1p1g12.5,3x,
     &      " # energy of beam [MeV/n]")') se0(j)

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of x-axis [cm]")') sx0(j)

            write(iot,'("       x1 = ",1p1g12.5,3x,
     &      " # Radius of x-axis [cm]")') sx1(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of y-axis [cm]")') sy0(j)

            write(iot,'("       y1 = ",1p1g12.5,3x,
     &      " # Radius of y-axis [cm]")') sy1(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of z-axis [cm]")') sz0(j)

            write(iot,'("       z1 = ",1p1g12.5,3x,
     &      " # Radius of z-axis [cm]")') sz1(j)

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 8 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # parabola source")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            write(iot,'("       rn = ",i3,12x,
     &      " # (D=2) pawer of parabola")') isrn(j)

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of x-axis [cm]")') sx0(j)

            write(iot,'("       x1 = ",1p1g12.5,3x,
     &      " # Radius of x-axis [cm]")') sx1(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of y-axis [cm]")') sy0(j)

            write(iot,'("       y1 = ",1p1g12.5,3x,
     &      " # Radius of y-axis [cm]")') sy1(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of z-axis [cm]")') sz0(j)

            write(iot,'("       z1 = ",1p1g12.5,3x,
     &      " # Radius of z-axis [cm]")') sz1(j)

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 9 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # spherical distribution source")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            write(iot,'("       e0 = ",1p1g12.5,3x,
     &      " # energy of beam [MeV/n]")') se0(j)

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of x-axis [cm]")') sx0(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of y-axis [cm]")') sy0(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of z-axis [cm]")') sz0(j)

            write(iot,'("       r1 = ",1p1g12.5,3x,
     &      " # inner radius [cm]")') sr1(j)

            write(iot,'("       r2 = ",1p1g12.5,3x,
     &      " # outer radius [cm]")') sr2(j)

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # direction of beam [isotropic]")')

            else if( sdir(j) .gt. 0.0 .and. sdir(j) .lt. 1.5 ) then

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # direction of beam is outer normal to sphere")')
     &            sdir(j)

            else if( sdir(j) .lt. 0.0 .and. sdir(j) .gt. -1.5 ) then

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # direction of beam is inner normal to sphere")')
     &            sdir(j)

            else if( sdir(j) .ge. 1.5 .and. sdir(j) .lt. 2.5 ) then

               write(iot,'("      dir =  +all",10x,
     &         " # outer direction with cos distribution")')

            else if( sdir(j) .le. -1.5 .and. sdir(j) .gt. -2.5 ) then

               write(iot,'("      dir =  -all",10x,
     &         " # inner direction with cos dis. and cos^2 bias")')

            else if( sdir(j) .le. -2.5 .and. sdir(j) .gt. -3.5 ) then

               write(iot,'("      dir =  iso ",10x,
     &         " # inner direction with uniform dis. by analog")')

               write(iot,'("      ag1 = ",1p1g12.5,3x,
     &         " # minimum cutoff cosine")') sag1(j)

               write(iot,'("      ag2 = ",1p1g12.5,3x,
     &         " # maximum cutoff cosine")') sag2(j)

               write(iot,'("      pg1 = ",1p1g12.5,3x,
     &         " # minimum azimuth angle in degree")') spg1(j)

               write(iot,'("      pg2 = ",1p1g12.5,3x,
     &         " # maximum azimuth angle in degree")') spg2(j)

               if(isbias(j).ge.1) then ! T.Sato 2020/08/09
                 write(iot,'("   isbias = ",i3, 12x,
     &           " # (D=0) 0:photon with no energy, 1:resampling" )')
     &           isbias(j)
               endif

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 10 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # spherical source")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of x-axis [cm]")') sx0(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of y-axis [cm]")') sy0(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of z-axis [cm]")') sz0(j)

            write(iot,'("       r1 = ",1p1g12.5,3x,
     &      " # inner radius [cm]")') sr1(j)

            write(iot,'("       r2 = ",1p1g12.5,3x,
     &      " # outer radius [cm]")') sr2(j)

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # direction of beam [isotropic]")')

            else if( sdir(j) .gt. 0.0 .and. sdir(j) .lt. 1.5 ) then

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # direction of beam is outer normal to sphere")')
     &            sdir(j)

            else if( sdir(j) .lt. 0.0 .and. sdir(j) .gt. -1.5 ) then

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # direction of beam is inner normal to sphere")')
     &            sdir(j)

            else if( sdir(j) .ge. 1.5 .and. sdir(j) .lt. 2.5 ) then

               write(iot,'("      dir =  +all",10x,
     &         " # outer direction with cos distribution")')

            else if( sdir(j) .le. -1.5 .and. sdir(j) .gt. -2.5 ) then

               write(iot,'("      dir =  -all",10x,
     &         " # inner direction with cos dis. and cos^2 bias")')

            else if( sdir(j) .le. -2.5 .and. sdir(j) .gt. -3.5 ) then

               write(iot,'("      dir =  iso ",10x,
     &         " # inner direction with uniform dis. by analog")')

               write(iot,'("      ag1 = ",1p1g12.5,3x,
     &         " # minimum zenith angle in cosine")') sag1(j)

               write(iot,'("      ag2 = ",1p1g12.5,3x,
     &         " # maximum zenith angle in cosine")') sag2(j)

               write(iot,'("      pg1 = ",1p1g12.5,3x,
     &         " # minimum azimuth angle in degree")') spg1(j)

               write(iot,'("      pg2 = ",1p1g12.5,3x,
     &         " # maximum azimuth angle in degree")') spg2(j)

               if(isbias(j).ge.1) then ! T.Sato 2020/08/09
                 write(iot,'("   isbias = ",i3, 12x,
     &           " # (D=0) 0:photon with no energy, 1:resampling" )')
     &           isbias(j)
               endif

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 11 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # ellipsoid source")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            write(iot,'("       e0 = ",1p1g12.5,3x,
     &      " # energy of beam [MeV/n]")') se0(j)

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center of x-axis [cm]")') sx0(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center of y-axis [cm]")') sy0(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # minimum position of z-axis [cm]")') sz0(j)

            write(iot,'("       z1 = ",1p1g12.5,3x,
     &      " # maximum position of z-axis [cm]")') sz1(j)

            write(iot,'("       rx = ",1p1g12.5,3x,
     &      " # phase space angle of x-axis [rad]")') srx(j)

            write(iot,'("       ry = ",1p1g12.5,3x,
     &      " # phase space angle of y-axis [rad]")') sry(j)

         if(swem(j).ne.0.0) then ! T.Sato 2020/09/21

            write(iot,'("      wem = ",1p1g12.5,3x,
     &      " # emittance [cm*mrad]")') swem(j)

            write(iot,'("       x1 = ",1p1g12.5,3x,
     &      " # max radius / max angle of x-axis [cm/mrad]")') sx1(j)

            write(iot,'("       y1 = ",1p1g12.5,3x,
     &      " # max radius / max angle of y-axis [cm/mrad]")') sy1(j)

         else  ! Gaussian distribution

            write(iot,'("       x1 = ",1p1g12.5,3x,
     &      " # Sigma of x-corrdinate Gaussian [cm]")') sx1(j)

            write(iot,'("   xmrad1 = ",1p1g12.5,3x,
     &      " # Sigma of x-angle Gaussian [mrad]")') sxmrad1(j)

            write(iot,'("       y1 = ",1p1g12.5,3x,
     &      " # Sigma of y-corrdinate Gaussian [cm]")') sy1(j)

            write(iot,'("   ymrad1 = ",1p1g12.5,3x,
     &      " # Sigma of y-angle Gaussian [mrad]")') symrad1(j)

         end if

            write(iot,'("       x2 = ",1p1g12.5,3x,
     &      " # Center of phase space x-corrdinate [cm]")') sx2(j)

            write(iot,'("   xmrad2 = ",1p1g12.5,3x,
     &      " # Center of phase space x-angle [mrad]")') sxmrad2(j)

            write(iot,'("       y2 = ",1p1g12.5,3x,
     &      " # Center of phase space y-corrdinate [cm]")') sx2(j)

            write(iot,'("   ymrad2 = ",1p1g12.5,3x,
     &      " # Center of phase space y-angle [mrad]")') symrad2(j)

            write(iot,'("      dir = ",1p1g12.5,3x,
     &      " # z-direction of beam [cosine]")') sdir(j)

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 12 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # external source with Decay-Turtle output")')
     &            jstypori(j) ! T.Sato 2017/07/14

            asfil = '  # file name of decay-turtle out put'

            msfile = max( 13, lsfile(j) )

            write(iot,'("     file =  ",100a1)')
     &                  ( sfile(j)(i:i), i = 1, msfile ),
     &                  ( asfil(i:i), i = 1, 38 )

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center of x-axis [cm]")') sx0(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center of y-axis [cm]")') sy0(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center of z-axis [cm]")') sz0(j)

            write(iot,'("      dir = ",1p1g12.5,3x,
     &      " # z-direction of beam [cosine]")') sdir(j)

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 13 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # Gaussian distribution in XY plane source")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            write(iot,'("       e0 = ",1p1g12.5,3x,
     &      " # energy of beam [MeV/n]")') se0(j)

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of x-axis [cm]")') sx0(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of y-axis [cm]")') sy0(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # (D=0.0) minimum position of z-axis [cm]")') sz0(j)

            write(iot,'("       z1 = ",1p1g12.5,3x,
     &      " # (D=0.0) maximum position of z-axis [cm]")') sz1(j)

            write(iot,'("       r1 = ",1p1g12.5,3x,
     &      " # Full Width at Half Maximum of Gaussian [cm]")') sr1(j)

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 14 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # Gaussian distribution in XY plane source")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of x-axis [cm]")') sx0(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of y-axis [cm]")') sy0(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # (D=0.0) minimum position of z-axis [cm]")') sz0(j)

            write(iot,'("       z1 = ",1p1g12.5,3x,
     &      " # (D=0.0) maximum position of z-axis [cm]")') sz1(j)

            write(iot,'("       r1 = ",1p1g12.5,3x,
     &      " # Full Width at Half Maximum of Gaussian [cm]")') sr1(j)

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 15 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # parabola in XY plane source")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            write(iot,'("       rn = ",i3,12x,
     &      " # (D=2) pawer of parabola")') isrn(j)

            write(iot,'("       e0 = ",1p1g12.5,3x,
     &      " # energy of beam [MeV/n]")') se0(j)

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of x-axis [cm]")') sx0(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of y-axis [cm]")') sy0(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # (D=0.0) minmum position of z-axis [cm]")') sz0(j)

            write(iot,'("       z1 = ",1p1g12.5,3x,
     &      " # (D=0.0) maximum position of z-axis [cm]")') sz1(j)

            write(iot,'("       r1 = ",1p1g12.5,3x,
     &      " # Radius of parabola [cm]")') sr1(j)

            if( sdom(j) .gt. -1000.0 ) then

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 16 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # parabola source in XY plane")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            write(iot,'("       rn = ",i3,12x,
     &      " # (D=2) pawer of parabola")') isrn(j)

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of x-axis [cm]")') sx0(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # (D=0.0) center position of y-axis [cm]")') sy0(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # (D=0.0) minimum position of z-axis [cm]")') sz0(j)

            write(iot,'("       z1 = ",1p1g12.5,3x,
     &      " # (D=0.0) maximum position of z-axis [cm]")') sz1(j)

            write(iot,'("       r1 = ",1p1g12.5,3x,
     &      " # Radius of parabola [cm]")') sr1(j)

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 18 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # cone source")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            write(iot,'("       e0 = ",1p1g12.5,3x,
     &      " # energy of beam [MeV/n]")') se0(j)

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # (D=0.0) top position of x-axis [cm]")') sx0(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # (D=0.0) top position of y-axis [cm]")') sy0(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # (D=0.0) top position of z-axis [cm]")') sz0(j)

            write(iot,'("       x1 = ",1p1g12.5,3x,
     &      " # x-vector from top to bottom")') sx1(j)

            write(iot,'("       y1 = ",1p1g12.5,3x,
     &      " # y-vector from top to bottom")') sy1(j)

            write(iot,'("       z1 = ",1p1g12.5,3x,
     &      " # z-vector from top to bottom")') sz1(j)

            write(iot,'("       r0 = ",1p1g12.5,3x,
     &      " # (D=0.0) Distance from top to upper suf [cm]")') sr0(j)

            write(iot,'("       r1 = ",1p1g12.5,3x,
     &      " # Distance from top to lower suface [cm]")') sr1(j)

            write(iot,'("       r2 = ",1p1g12.5,3x,
     &      " # Angle of cone [degree]")') sr2(j)

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 19 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # cone source")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # (D=0.0) top position of x-axis [cm]")') sx0(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # (D=0.0) top position of y-axis [cm]")') sy0(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # (D=0.0) top position of z-axis [cm]")') sz0(j)

            write(iot,'("       x1 = ",1p1g12.5,3x,
     &      " # x-vector from top to bottom")') sx1(j)

            write(iot,'("       y1 = ",1p1g12.5,3x,
     &      " # y-vector from top to bottom")') sy1(j)

            write(iot,'("       z1 = ",1p1g12.5,3x,
     &      " # z-vector from top to bottom")') sz1(j)

            write(iot,'("       r0 = ",1p1g12.5,3x,
     &      " # (D=0.0) Distance from top to upper suf [cm]")') sr0(j)

            write(iot,'("       r1 = ",1p1g12.5,3x,
     &      " # Distance from top to lower suface [cm]")') sr1(j)

            write(iot,'("       r2 = ",1p1g12.5,3x,
     &      " # Angle of cone [degree]")') sr2(j)

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 20 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # wedge source")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1: # of SF, 2:neutron" )')
     &         ispfs(j)

          end if

            write(iot,'("       e0 = ",1p1g12.5,3x,
     &      " # energy of beam [MeV/n]")') se0(j)

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # x of top position")') sx0(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # y of top position")') sy0(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # z of top position")') sz0(j)

            write(iot,'("       x1 = ",1p1g12.5,3x,
     &      " # x of the first triangle corner")') sx1(j)

            write(iot,'("       y1 = ",1p1g12.5,3x,
     &      " # y of the first triangle corner")') sy1(j)

            write(iot,'("       z1 = ",1p1g12.5,3x,
     &      " # z of the first triangle corner")') sz1(j)

            write(iot,'("       x2 = ",1p1g12.5,3x,
     &      " # x of the second triangle corner")') sr0(j)

            write(iot,'("       y2 = ",1p1g12.5,3x,
     &      " # y of the second triangle corner")') sr1(j)

            write(iot,'("       z2 = ",1p1g12.5,3x,
     &      " # z of the second triangle corner")') sr2(j)

            write(iot,'("       x3 = ",1p1g12.5,3x,
     &      " # x of the bottom position")') srx(j)

            write(iot,'("       y3 = ",1p1g12.5,3x,
     &      " # y of the bottom position")') sry(j)

            write(iot,'("       z3 = ",1p1g12.5,3x,
     &      " # z of the bottom position")') swem(j)

            write(iot,'("      exa = ",1p1g12.5,3x,
     &      " # exp(-ax): a>0.0, 0.0:uniform")') sdl0(j)

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 21 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # wedge source")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1: # of SF, 2:neutron" )')
     &         ispfs(j)

          end if

            write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # x of top position")') sx0(j)

            write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # y of top position")') sy0(j)

            write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # z of top position")') sz0(j)

            write(iot,'("       x1 = ",1p1g12.5,3x,
     &      " # x of the first triangle corner")') sx1(j)

            write(iot,'("       y1 = ",1p1g12.5,3x,
     &      " # y of the first triangle corner")') sy1(j)

            write(iot,'("       z1 = ",1p1g12.5,3x,
     &      " # z of the first triangle corner")') sz1(j)

            write(iot,'("       x2 = ",1p1g12.5,3x,
     &      " # x of the second triangle corner")') sr0(j)

            write(iot,'("       y2 = ",1p1g12.5,3x,
     &      " # y of the second triangle corner")') sr1(j)

            write(iot,'("       z2 = ",1p1g12.5,3x,
     &      " # z of the second triangle corner")') sr2(j)

            write(iot,'("       x3 = ",1p1g12.5,3x,
     &      " # x of the bottom position")') srx(j)

            write(iot,'("       y3 = ",1p1g12.5,3x,
     &      " # y of the bottom position")') sry(j)

            write(iot,'("       z3 = ",1p1g12.5,3x,
     &      " # z of the bottom position")') swem(j)

            write(iot,'("      exa = ",1p1g12.5,3x,
     &      " # exp(-ax): a>0.0, 0.0:uniform")') sdl0(j)

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 22 .or. jstyp(j) .eq. 23 ) then

          if( jstyp(j) .eq. 22  ) then
            write(iot,'("   s-type = ",i3,8x,
     &      "     # mesh-weighted source")')
     &            jstypori(j) ! T.Sato 2017/07/14

            write(iot,'("       e0 = ",1p1g12.5,3x,
     &      " # energy of beam [MeV/n]")') se0(j)

          else if( jstyp(j) .eq. 23  ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # mesh-weighted source with energy spectrum")')
     &            jstyp(j)

          end if

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

            if(isbias(j).ge.1) then ! T.Sato 2020/08/09
               write(iot,'("   isbias = ",i3, 12x,
     &      " # Source bias method, 0:number, 1:weight, 2:equal" )')
     &      isbias(j)
            endif

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh only")')

                  call echmty(0,iot,'x',isxtp(j),isinx(j),
     &                        sxdel(j),sxmin(j),sxmax(j),istxx(j)
     &                        ,abs(isinx(j)+1),das_istxx(1))  ! T.Sato 2025/02/02, bug fix
                  call echmty(0,iot,'y',isytp(j),isiny(j),
     &                        sydel(j),symin(j),symax(j),istyy(j)
     &                        ,abs(isiny(j)+1),das_istyy(1))  ! T.Sato 2025/02/02, bug fix
                  call echmty(0,iot,'z',isztp(j),isinz(j),
     &                        szdel(j),szmin(j),szmax(j),istzz(j)
     &                        ,abs(isinz(j)+1),das_istzz(1))  ! T.Sato 2025/02/02, bug fix

                  write(iot,'("# relative strength : ",
     &            "((((x,y,z),ix=1,nx),iy=1,ny),iz=1,nz) ")')

                  inxyz = iabs(isinx(j)*isiny(j)*isinz(j))
                  write(iot,'(1p7e11.3)') ! S.H. changed 10e->7e (2017.2.9)
     &            ( das(istdd(j)+inxyz+k),k=1,inxyz )

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 24 .or. jstyp(j) .eq. 25 ) then

          write(iot,'("   s-type = ",i3,8x,
     &      "     # tetra-mesh source")')
     &         jstypori(j)

          write(iot,'("   tetreg =  ",i5, 9x,
     &      " # tetra-mesh source region")') itetreg(j)

          if( jstyp(j) .eq. 24  ) then

            write(iot,'("       e0 = ",1p1g12.5,3x,
     &      " # energy of beam [MeV/n]")') se0(j)

          endif

            write(iot,'("     proj =  ",a51)') rname

          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------
cKN 2018/10/29 for surface source

         else if( jstyp(j) .eq. 26 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # surface source")')
     &            jstypori(j)  ! T.Sato 2017/07/04

            write(iot,'("     proj =  ",a51)') rname

          if( se0(j) .gt. 0.0 ) then

            write(iot,'("       e0 = ",1p1g12.5,3x,
     &      " # energy of beam [MeV/n]")') se0(j)

          endif
          if( lstyp(j) .gt. -1000.d0 ) then

            write(iot,'("     izst =  ",i3, 11x,
     &      " # charge state of source particle")') lstyp(j)

          end if
          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # direction to surface by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # direction to surface [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # direction to surface [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            if( swt0(j) .ne. 1.0d0 ) then

               write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

            end if

*-----------------------------------------------------------------------

            write(iot,'("      suf = ",i3,8x,
     &      "     # surface number")') issuf(j)

          if( iscut(j) .gt. 0 ) then

                jk = 0
             do kk = 1, iscut(j)

                if(      iabs(isvct(j,kk)) .lt. 10 ) then
                   write(dum1(jk+1:jk+3),'(i3)') isvct(j,kk)
                   jk = jk + 3
                else if( iabs(isvct(j,kk)) .lt. 100 ) then
                   write(dum1(jk+1:jk+4),'(i4)') isvct(j,kk)
                   jk = jk + 4
                else if( iabs(isvct(j,kk)) .lt. 1000 ) then
                   write(dum1(jk+1:jk+5),'(i5)') isvct(j,kk)
                   jk = jk + 5
                else if( iabs(isvct(j,kk)) .lt. 10000 ) then
                   write(dum1(jk+1:jk+6),'(i6)') isvct(j,kk)
                   jk = jk + 6
                else if( iabs(isvct(j,kk)) .lt. 100000 ) then
                   write(dum1(jk+1:jk+7),'(i7)') isvct(j,kk)
                   jk = jk + 7
                else
                   write(dum1(jk+1:jk+8),'(i8)') isvct(j,kk)
                   jk = jk + 8
                end if

             end do

            write(iot,'("      cut = ",64(a1),8x,
     &      "     # cut surfaces definition")') (dum1(kk:kk),kk=1,jk)

          end if

cKN 2018/10/29 for surface source
*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 17 ) then

            write(iot,'("   s-type = ",i3,8x,
     &      "     # external source with PHITS dump file")')
     &            jstypori(j) ! T.Sato 2017/07/14

            asfil = '  # file name of dump file'

            msfile = max( 13, lsfile(j) )

            write(iot,'("     file =  ",100a1)')
     &                  ( sfile(j)(i:i), i = 1, msfile ),
     &                  ( asfil(i:i), i = 1, 38 )

            if(jpsf.eq.0)then !FURUTA20201110

             write(iot,'("     dump =",i5,11x,
     &        " # number of dumped data",
     &        " <0: ascii, >0: binary")') isdmp(j,0)

             write(iot,'(13x,30(i4))')
     &            ( isdmp(j,k), k = 1, abs( isdmp(j,0) ) )

             write(iot,'("# dump data  ",30(a4))')
     &            ( dmpc(isdmp(j,k)), k = 1, abs( isdmp(j,0) ) )
             write(iot,'("   idmpmode =",i3,13x,
     &         " # (D=1) 0: off, 1: consider event number of",
     &         " dump source")') idmpmode

             if(idmpmode.eq.1)then
              write(iot,'("#  totfact is ignored with idmpmode=1")')
             endif

             write(iot,'("  dmpmulti =",1p1g12.5,3x,
     &         " # multiplication factor for dump source")') dmpmulti
c------------------
            else              !FURUTA20201110
             write(iot,'("     jpsf =",i5,11x,
     &        " # (D=0) 0: off, >0: special option for",
     &        "  Phase Space File")') jpsf
            endif

            if( jsdmp(j,1)  .eq. 0 ) then

                  write(iot,'("     proj =  ",a51)') rname
               if( lstyp(j) .gt. -1000.d0 ) then

                  write(iot,'("     izst =  ",i3, 11x,
     &            " # charge state of source particle")') lstyp(j)

               end if

          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if
            end if

            if( jsdmp(j,2)  .eq. 0 )
     &      write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # x-coordinate [cm]")') sx0(j)

            if( jsdmp(j,3)  .eq. 0 )
     &      write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # y-coordinate [cm]")') sy0(j)

            if( jsdmp(j,4)  .eq. 0 )
     &      write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # z-coordinate [cm]")') sz0(j)

            if( jsdmp(j,5) .eq. 0 .or. jsdmp(j,6) .eq. 0 .or.
     &          jsdmp(j,7) .eq. 0 ) then

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            end if

            if( jsdmp(j,8)  .eq. 0 .and. se0(j) .gt. 0.0 )
     &      write(iot,'("       e0 = ",1p1g12.5,3x,
     &      " # energy of beam [MeV/n]")') se0(j)

            if( jsdmp(j,9)  .eq. 0 )
     &         write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

*-----------------------------------------------------------------------

         else if( jstyp(j) .eq. 100 ) then

            write(iot,'("   s-type = ",i3,12x,
     &      " # user defined source.")')
     &            jstypori(j) ! T.Sato 2017/07/14

            if( jsusr(j,1)  .eq. 1 ) then
            write(iot,'("     proj =  ",a51)') rname
               if( lstyp(j) .gt. -1000.d0 ) then

                  write(iot,'("     izst =  ",i3, 11x,
     &            " # charge state of source particle")') lstyp(j)

               end if

          if(ispfs(j) .ne. 0 ) then

            write(iot,'("    ispfs =  ",i2, 12x,
     &      " # Spontaneus Fission Source, 1:neutron, 2: # of SF" )')
     &         ispfs(j)

          end if
            end if

            if( jsusr(j,2)  .eq. 1 )
     &      write(iot,'("       x0 = ",1p1g12.5,3x,
     &      " # x-coordinate [cm]")') sx0(j)

            if( jsusr(j,3)  .eq. 1 )
     &      write(iot,'("       y0 = ",1p1g12.5,3x,
     &      " # y-coordinate [cm]")') sy0(j)

            if( jsusr(j,4)  .eq. 1 )
     &      write(iot,'("       z0 = ",1p1g12.5,3x,
     &      " # z-coordinate [cm]")') sz0(j)

            if( jsusr(j,5) .eq. 1 .or. jsusr(j,6) .eq. 1 .or.
     &          jsusr(j,7) .eq. 1 ) then

            if( sdir(j) .gt. 250.0 ) then

               write(iot,'("      dir =  data",10x,
     &         " # z-direction of beam by below data")')

               call echsrsa(j,iot,ierr)

            else if( sdir(j) .gt. 100.0 ) then

               write(iot,'("      dir =   all",10x,
     &         " # z-direction of beam [isotropic]")')

            else

               write(iot,'("      dir = ",1p1g12.5,3x,
     &         " # z-direction of beam [cosine]")') sdir(j)

            end if

            if( sphi(j) .gt. -1000.0 ) then

               write(iot,'("      phi = ",1p1g12.5,3x,
     &         " # azimuthal angle of beam [degree]")') sphi(j)

            end if

            if( sdom(j) .gt. -1000.0 ) then

               write(iot,'("      dom = ",1p1g12.5,3x,
     &         " # spread of beam direction [degree]")') sdom(j)

            end if

            end if

            if( jsusr(j,8)  .eq. 1 .and. se0(j) .gt. 0.0 )
     &      write(iot,'("       e0 = ",1p1g12.5,3x,
     &      " # energy of beam [MeV/n]")') se0(j)

            if( jsusr(j,9)  .eq. 1 )
     &         write(iot,'("      wgt = ",1p1g12.5,3x,
     &         " # weight of the particle")') swt0(j)

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*     source energy group
*-----------------------------------------------------------------------

      if( jstyp(j) .eq.   4 .or.
     &    jstyp(j) .eq.   5 .or.
     &    jstyp(j) .eq.   6 .or.
     &    jstyp(j) .eq.   8 .or.
     &    jstyp(j) .eq.  10 .or.
     &    jstyp(j) .eq.  14 .or.
     &    jstyp(j) .eq.  16 .or.
     &    jstyp(j) .eq.  19 .or.
     &    jstyp(j) .eq.  21 .or.
     &    jstyp(j) .eq.  23 .or.
     &    jstyp(j) .eq.  25 .or.
     &  ( jstyp(j) .eq.  26 .and. jetyp(j) .gt. 0 ) .or.

     &  ( jstyp(j) .eq.  17 .and. jetyp(j) .gt. 0 ) .or.
     &  ( jstyp(j) .eq. 100 .and. jetyp(j) .gt. 0 ) ) then

         if( jetyp(j) .eq. 1 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # energy distribution given by data")') jetyp(j)

            write(iot,'("       ne = ",i4,11x,
     &      " # number of energy and weight"/27x,
     &      " # ne>0: lin, ne<0: log interpolation"/
     &      " #  data = ( e(i), w(i), i = 1, ne ), emax")')
     &      ngrp(j) * ngll(j)

            write(iot,'(1p6e13.5)') (egmin(ngei(j)+i),
     &                               fegrp(ngfe(j)+i),i=1,ngrp(j)),
     &                               egmax(ngea(j)+ngrp(j))

         else if( jetyp(j) .eq. 11 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # wave length [A] distribution given by data")')
     &         jetyp(j)

            write(iot,'("       ne = ",i4,11x,
     &      " # number of wave length [A] and weight"/
     &      " # ne>0: lin, ne<0: log interpolation"/
     &      " #  data = ( e(i), w(i), i = 1, ne ), emax")')
     &      ngrp(j) * ngll(j)

            write(iot,'(1p6e13.5)') (egmin(ngei(j)+i),
     &                               fegrp(ngfe(j)+i),i=1,ngrp(j)),
     &                               egmax(ngea(j)+ngrp(j))

         else if( jetyp(j) .eq. 4 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # energy distribution given by data")') jetyp(j)

            write(iot,'("       ne = ",i4,11x,
     &      " # number of energy and weight"/27x,
     &      " # ne>0: lin, ne<0: log interpolation"/
     &      " #  data = ( e(i), w(i), i = 1, ne ), emax")')
     &      ngrp(j) * ngll(j)

            write(iot,'(1p6e13.5)') (egmin(ngei(j)+i),
     &                               fegrp(ngfe(j)+i),i=1,ngrp(j)),
     &                               egmax(ngea(j)+ngrp(j))

            write(iot,'("   p-type = ",i3,12x,
     &      " # generate weight type, 0:equal, 1: given by data")')
     &         jptyp(j)

          if( jptyp(j) .eq. 1 ) then
            write(iot,'(1p6e13.5)') (prw(ngpi(j)+i),i=1,ngrp(j))
          end if

         else if( jetyp(j) .eq. 14 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # wave length [A] distribution given by data")')
     &           jetyp(j)

            write(iot,'("       ne = ",i4,11x,
     &      " # number of wave length [A] and weight"/
     &      " # ne>0: lin, ne<0: log interpolation"/
     &      " #  data = ( e(i), w(i), i = 1, ne ), emax")')
     &      ngrp(j) * ngll(j)

            write(iot,'(1p6e13.5)') (egmin(ngei(j)+i),
     &                               fegrp(ngfe(j)+i),i=1,ngrp(j)),
     &                               egmax(ngea(j)+ngrp(j))

            write(iot,'("   p-type = ",i3,12x,
     &      " # generate weight type, 0:equal, 1: given by data")')
     &         jptyp(j)

          if( jptyp(j) .eq. 1 ) then
            write(iot,'(1p6e13.5)') (prw(ngpi(j)+i),i=1,ngrp(j))
          end if

*-----------------------------------------------------------------------

         else if( jetyp(j) .eq. 8 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # pointwise energies given by data")') jetyp(j)

            write(iot,'("       ne = ",i4,11x,
     &      " # number of energy and weight"/27x,
     &      " #  data = ( e(i), w(i), i = 1, ne )")')
     &      ngrp(j)

            write(iot,'(1p6e13.5)') (egmin(ngei(j)+i),
     &                               fegrp(ngfe(j)+i),i=1,ngrp(j))

         else if( jetyp(j) .eq. 18 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # pointwise wave length [A] given by data")')
     &         jetyp(j)

            write(iot,'("       ne = ",i4,11x,
     &      " # number of wave length [A] and weight"/
     &      " #  data = ( e(i), w(i), i = 1, ne )")')
     &      ngrp(j)

            write(iot,'(1p6e13.5)') (egmin(ngei(j)+i),
     &                               fegrp(ngfe(j)+i),i=1,ngrp(j))

         else if( jetyp(j) .eq. 9 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # pointwise energies given by data")') jetyp(j)

            write(iot,'("       ne = ",i4,11x,
     &      " # number of energy and weight"/27x,
     &      " #  data = ( e(i), w(i), i = 1, ne )")')
     &      ngrp(j)

            write(iot,'(1p6e13.5)') (egmin(ngei(j)+i),
     &                               fegrp(ngfe(j)+i),i=1,ngrp(j))

            write(iot,'("   p-type = ",i3,12x,
     &      " # generate weight type, 0:equal, 1: given by data")')
     &         jptyp(j)

          if( jptyp(j) .eq. 1 ) then
            write(iot,'(1p6e13.5)') (prw(ngpi(j)+i),i=1,ngrp(j))
          end if

         else if( jetyp(j) .eq. 19 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # pointwize wave length [A] given by data")')
     &           jetyp(j)

            write(iot,'("       ne = ",i4,11x,
     &      " # number of wave length [A] and weight"/
     &      " #  data = ( e(i), w(i), i = 1, ne )")')
     &      ngrp(j)

            write(iot,'(1p6e13.5)') (egmin(ngei(j)+i),
     &                               fegrp(ngfe(j)+i),i=1,ngrp(j))

            write(iot,'("   p-type = ",i3,12x,
     &      " # generate weight type, 0:equal, 1: given by data")')
     &         jptyp(j)

          if( jptyp(j) .eq. 1 ) then
            write(iot,'(1p6e13.5)') (prw(ngpi(j)+i),i=1,ngrp(j))
          end if

*-----------------------------------------------------------------------

         else if( jetyp(j) .eq. 22 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # emin and emax energies given by data")') jetyp(j)

            write(iot,'("       ne = ",i4,11x,
     &      " # number of energy and weight"/27x,
     &      " #  data = ( emin(i), emax(i), w(i), i = 1, ne )")')
     &      ngrp(j)

            write(iot,'(1p6e13.5)') (egmin(ngei(j)+i),egmax(ngea(j)+i),
     &                               fegrp(ngfe(j)+i),i=1,ngrp(j))

         else if( jetyp(j) .eq. 32 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # min and max wave length [A] given by data")')
     &         jetyp(j)

            write(iot,'("       ne = ",i4,11x,
     &      " # number of wave length [A] and weight"/
     &      " #  data = ( emin(i), emax(i), w(i), i = 1, ne )")')
     &      ngrp(j)

            write(iot,'(1p6e13.5)') (egmin(ngei(j)+i),egmax(ngea(j)+i),
     &                               fegrp(ngfe(j)+i),i=1,ngrp(j))

         else if( jetyp(j) .eq. 23 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # emin and emax energies given by data")') jetyp(j)

            write(iot,'("       ne = ",i4,11x,
     &      " # number of energy and weight"/27x,
     &      " #  data = ( emin(i), emax(i), w(i), i = 1, ne )")')
     &      ngrp(j)

            write(iot,'(1p6e13.5)') (egmin(ngei(j)+i),egmax(ngea(j)+i),
     &                               fegrp(ngfe(j)+i),i=1,ngrp(j))

            write(iot,'("   p-type = ",i3,12x,
     &      " # generate weight type, 0:equal, 1: given by data")')
     &         jptyp(j)

          if( jptyp(j) .eq. 1 ) then
            write(iot,'(1p6e13.5)') (prw(ngpi(j)+i),i=1,ngrp(j))
          end if

         else if( jetyp(j) .eq. 33 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # min and max wave length [A] given by data")')
     &           jetyp(j)

            write(iot,'("       ne = ",i4,11x,
     &      " # number of wave length [A] and weight"/
     &      " #  data = ( emin(i), emax(i), w(i), i = 1, ne )")')
     &      ngrp(j)

            write(iot,'(1p6e13.5)') (egmin(ngei(j)+i),egmax(ngea(j)+i),
     &                               fegrp(ngfe(j)+i),i=1,ngrp(j))

            write(iot,'("   p-type = ",i3,12x,
     &      " # generate weight type, 0:equal, 1: given by data")')
     &         jptyp(j)

          if( jptyp(j) .eq. 1 ) then
            write(iot,'(1p6e13.5)') (prw(ngpi(j)+i),i=1,ngrp(j))
          end if

*-----------------------------------------------------------------------
         else if( jetyp(j) .eq. 21 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # energy distribution given by data")') jetyp(j)

            write(iot,'("       ne = ",i4,11x,
     &      " # number of energy and dN/dE"/27x,
     &      " # ne>0: lin, ne<0: log interpolation"/
     &      " #  data = ( e(i), dN/dE(i), i = 1, ne ), emax")')
     &      ngrp(j) * ngll(j)

            write(iot,'(1p6e13.5)') (egmin(ngei(j)+i),
     &                               fegrp(ngfe(j)+i),i=1,ngrp(j)),
     &                               egmax(ngea(j)+ngrp(j))

         else if( jetyp(j) .eq. 31 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # wave length [A] distribution given by data")')
     &         jetyp(j)

            write(iot,'("       ne = ",i4,11x,
     &      " # number of wave length [A] and dN/dA"/
     &      " # ne>0: lin, ne<0: log interpolation"/
     &      " #  data = ( e(i), dN/dA(i), i = 1, ne ), emax")')
     &      ngrp(j) * ngll(j)

            write(iot,'(1p6e13.5)') (egmin(ngei(j)+i),
     &                               fegrp(ngfe(j)+i),i=1,ngrp(j)),
     &                               egmax(ngea(j)+ngrp(j))

         else if( jetyp(j) .eq. 24 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # energy distribution given by data")') jetyp(j)

            write(iot,'("       ne = ",i4,11x,
     &      " # number of energy and dN/dE"/27x,
     &      " # ne>0: lin, ne<0: log interpolation"/
     &      " #  data = ( e(i), w(i), i = 1, ne ), emax")')
     &      ngrp(j) * ngll(j)

            write(iot,'(1p6e13.5)') (egmin(ngei(j)+i),
     &                               fegrp(ngfe(j)+i),i=1,ngrp(j)),
     &                               egmax(ngea(j)+ngrp(j))

            write(iot,'("   p-type = ",i3,12x,
     &      " # generate weight type, 0:equal, 1: given by data")')
     &         jptyp(j)

          if( jptyp(j) .eq. 1 ) then
            write(iot,'(1p6e13.5)') (prw(ngpi(j)+i),i=1,ngrp(j))
          end if

         else if( jetyp(j) .eq. 34 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # wave length [A] distribution given by data")')
     &           jetyp(j)

            write(iot,'("       ne = ",i4,11x,
     &      " # number of wave length [A] and dN/dA"/
     &      " # ne>0: lin, ne<0: log interpolation"/
     &      " #  data = ( e(i), w(i), i = 1, ne ), emax")')
     &      ngrp(j) * ngll(j)

            write(iot,'(1p6e13.5)') (egmin(ngei(j)+i),
     &                               fegrp(ngfe(j)+i),i=1,ngrp(j)),
     &                               egmax(ngea(j)+ngrp(j))

            write(iot,'("   p-type = ",i3,12x,
     &      " # generate weight type, 0:equal, 1: given by data")')
     &         jptyp(j)

          if( jptyp(j) .eq. 1 ) then
            write(iot,'(1p6e13.5)') (prw(ngpi(j)+i),i=1,ngrp(j))
          end if

*-----------------------------------------------------------------------

         else if( jetyp(j) .eq. 2 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # Gaussian energy distribution")') jetyp(j)

            write(iot,'("      eg0 = ",1p1g12.5,3x,
     &      " # center of Gaussian [MeV]")') seg0(j)

            write(iot,'("      eg1 = ",1p1g12.5,3x,
     &      " # Full Width at Half Maximum of Gaussian [MeV]")')
     &         seg1(j)

            write(iot,'("      eg2 = ",1p1g12.5,3x,
     &      " # minimum cutoff energy of Gaussian [MeV]")') seg2(j)

            write(iot,'("      eg3 = ",1p1g12.5,3x,
     &      " # maximum cutoff energy of Gaussian [MeV]")') seg3(j)

            if( seg3(j) .le. seg2(j) ) then

               write(io,'(/"Error: eg3 should be > eg2")')
               ErrCha = ''
               MsgID = 'L:21445/R:echoi/F:read00.f'
               call ErrWrite(MsgID, ErrCha)
               write(jo,'(/"Error: eg3 should be > eg2")')
               stop 844

            end if

         else if( jetyp(j) .eq. 12 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # Gaussian wave lenght [A] distribution")') jetyp(j)

            write(iot,'("      eg0 = ",1p1g12.5,3x,
     &      " # center of Gaussian [A]")') seg0(j)

            write(iot,'("      eg1 = ",1p1g12.5,3x,
     &      " # Full Width at Half Maximum of Gaussian [A]")') seg1(j)

            write(iot,'("      eg2 = ",1p1g12.5,3x,
     &      " # minimum cutoff wave length of Gaussian [A]")') seg2(j)

            write(iot,'("      eg3 = ",1p1g12.5,3x,
     &      " # maximum cutoff wave length of Gaussian [A]")') seg3(j)

            if( seg3(j) .le. seg2(j) ) then

               write(io,'(/"Error: eg3 should be > eg2")')
               ErrCha = ''
               MsgID = 'L:21473/R:echoi/F:read00.f'
               call ErrWrite(MsgID, ErrCha)
               write(jo,'(/"Error: eg3 should be > eg2")')
               stop 844

            end if

         else if( jetyp(j) .eq. 3 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # Maxwellian generated by probability")') jetyp(j)

            write(iot,'("      et0 = ",1p1g12.5,3x,
     &      " # temperature of Maxwellian [MeV]")') set0(j)

            write(iot,'("      et1 = ",1p1g12.5,3x,
     &      " # minimum cutoff energy of Maxwellian [MeV]")') set1(j)

            write(iot,'("      et2 = ",1p1g12.5,3x,
     &      " # maximum cutoff energy of Maxwellian [MeV]")') set2(j)

            write(iot,'("      et3 = ",1p1g12.5,3x,
     &      " # power index of Maxwellian")') set3(j)

            write(iot,'("       nm = ",i4,11x,
     &      " # (D=-200) number of energy mesh"/27x,
     &      " # nm>0: lin, nm<0: log mesh and interpolation")')
     &      ngrp(j) * ngll(j)

         else if( jetyp(j) .eq. 7 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # Maxwellian generated by p-type")') jetyp(j)

            write(iot,'("      et0 = ",1p1g12.5,3x,
     &      " # temperature of Maxwellian [MeV]")') set0(j)

            write(iot,'("      et1 = ",1p1g12.5,3x,
     &      " # minimum cutoff energy of Maxwellian [MeV]")') set1(j)

            write(iot,'("      et2 = ",1p1g12.5,3x,
     &      " # maximum cutoff energy of Maxwellian [MeV]")') set2(j)

            write(iot,'("      et3 = ",1p1g12.5,3x,
     &      " # power index of Maxwellian")') set3(j)

            write(iot,'("       nm = ",i4,11x,
     &      " # (D=-200) number of energy mesh"/27x,
     &      " # nm>0: lin, nm<0: log mesh and interpolation")')
     &      ngrp(j) * ngll(j)

            write(iot,'("   p-type = ",i3,12x,
     &      " # generate weight type, 0:equal, 1: given by data")')
     &         jptyp(j)

          if( jptyp(j) .eq. 1 ) then
            write(iot,'(1p6e13.5)') (prw(ngpi(j)+i),i=1,ngrp(j))
          end if

         else if( jetyp(j) .eq. 5 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # Functional energy distribution")') jetyp(j)

            write(iot,'("     f(x) = ",1024a1)')
     &                 ( srfx(j)(k:k), k = 1, lsfx(j) )

            write(iot,'("      eg1 = ",1p1g12.5,3x,
     &      " # minimum cutoff energy [MeV]")') seg1(j)

            write(iot,'("      eg2 = ",1p1g12.5,3x,
     &      " # maximum cutoff energy [MeV]")') seg2(j)

            write(iot,'("       nm = ",i4,11x,
     &      " # number of energy mesh"/27x,
     &      " # nm>0: lin, nm<0: log mesh and interpolation")')
     &      ngrp(j) * ngll(j)

         else if( jetyp(j) .eq. 6 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # Functional energy distribution")') jetyp(j)

            write(iot,'("     f(x) = ",1024a1)')
     &                 ( srfx(j)(k:k), k = 1, lsfx(j) )

            write(iot,'("      eg1 = ",1p1g12.5,3x,
     &      " # minimum cutoff energy [MeV]")') seg1(j)

            write(iot,'("      eg2 = ",1p1g12.5,3x,
     &      " # maximum cutoff energy [MeV]")') seg2(j)

            write(iot,'("       nm = ",i4,11x,
     &      " # number of energy mesh"/27x,
     &      " # nm>0: lin, nm<0: log mesh and interpolation")')
     &      ngrp(j) * ngll(j)

            write(iot,'("   p-type = ",i3,12x,
     &      " # generate weight type, 0:equal, 1: given by data")')
     &         jptyp(j)

          if( jptyp(j) .eq. 1 ) then
            write(iot,'(1p6e13.5)') (prw(ngpi(j)+i),i=1,ngrp(j))
          end if

         else if( jetyp(j) .eq. 15 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # Functional wave length [A]  distribution")') jetyp(j)

            write(iot,'("     f(x) = ",1024a1)')
     &                 ( srfx(j)(k:k), k = 1, lsfx(j) )

            write(iot,'("      eg1 = ",1p1g12.5,3x,
     &      " # minimum cutoff wave length [A]")') seg1(j)

            write(iot,'("      eg2 = ",1p1g12.5,3x,
     &      " # maximum cutoff wave length [A]")') seg2(j)

            write(iot,'("       nm = ",i4,11x,
     &      " # number of wave length mesh"/27x,
     &      " # nm>0: lin, nm<0: log mesh and interpolation")')
     &      ngrp(j) * ngll(j)

         else if( jetyp(j) .eq. 16 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # Functional wave length [A]  distribution")') jetyp(j)

            write(iot,'("     f(x) = ",1024a1)')
     &                 ( srfx(j)(k:k), k = 1, lsfx(j) )

            write(iot,'("      eg1 = ",1p1g12.5,3x,
     &      " # minimum cutoff wave length [A]")') seg1(j)

            write(iot,'("      eg2 = ",1p1g12.5,3x,
     &      " # maximum cutoff wave length [A]")') seg2(j)

            write(iot,'("       nm = ",i4,11x,
     &      " # number of wave length mesh"/27x,
     &      " # nm>0: lin, nm<0: log mesh and interpolation")')
     &      ngrp(j) * ngll(j)

            write(iot,'("   p-type = ",i3,12x,
     &      " # generate weight type, 0:equal, 1: given by data")')
     &         jptyp(j)

          if( jptyp(j) .eq. 1 ) then
            write(iot,'(1p6e13.5)') (prw(ngpi(j)+i),i=1,ngrp(j))
          end if

*-----------------------------------------------------------------------
         else if( jetyp(j) .eq. 28 .or. jetyp(j) .eq. 29 ) then
             nk0 = niorg(j)

             write(iot,'("     norm = ",i3,12x,
     &       " # (D=0) normalization factor =0:(/sec) =1:(/source)")')
     &       norm(j)

             write(iot,'("    dtime = ",1pe13.5,2x,
     &       " # (D=-10.0) cooling time: sec(+) or half_life(-)")')
     &       decayt(j)

           if( nk0 .ne. 0 ) then  ! MATSUDA 2017.05.29
             write(iot,'("$   dtime =   0.0000",7x,
     &       " # if you want to use like e-type = 22",
     &       " dtime must be 0.0")')
           end if

             write(iot,'("   actlow = ",1pe13.5,2x,
     &       " # lower limit of activity (Bq)")')
     &       aclow(j)

           if( inkf0(j) .eq. 11 ) then
             write(iot,'("  iaugers = ",i3,12x,
     &       " # (D=0) =0 with, 1 w/o, or 2 only Auger electron")')
     &       iaugers(j)
           end if

           if( inkf0(j) .eq. 22 ) then
             write(iot,'(" icharctx = ",i3,12x,
     &       " # (D=0) =0 with, 1 w/o, or 2 only character X-rays")')
     &       icharacterx(j)
           end if

           if( inkf0(j) .eq. 22 ) then
             write(iot,'("   iannih = ",i3,12x,
     &       " # (D=0) =0 with, or 1 w/o annihilation photons")')
     &       iannih(j)
           end if

             write(iot,'("   e-type = ",i3,12x,
     &       " # RI source")')
     &       jetyp(j)

*-----------------------------------------------------------------------
*     Read nuclides list
*-----------------------------------------------------------------------
           if( nk0 .ne. 0 ) then
             write(iot,'("       ni = ",i3,12x,
     &       " # number of registered nuclide"/27x,
     &       " #  data = ( nuclide(i), activity(i), i = 1, ni )")')
     &       nk0

             nk1 = nicur(j)
             do n = 1, nk1
               call ddc2echo(iot,n,j,ierr)
             end do

             if( asfsum(j) .gt. 0.0d+00 ) then
               write(iot,'("*** warning ***  Spontaneous fision: "
     &                    , 1pe13.5, " (fission/sec)")') asfsum(j)
             end if

! MATSUDA 2017.05.29
               write(iot,'("$      ne = ",i4,11x,
     &         " # number of energy and weight"/27x,
     &         " #  data = ( emin(i), emax(i), w(i), i = 1, ne )")')
     &         ngrp(j)

             if( jetyp(j) .eq. 28 ) then
               write(iot,'("$",1p6e13.5)') (egmin(ngei(j)+i),
     &               egmax(ngea(j)+i),fegrp(ngfe(j)+i),i=1,ngrp(j))
             else  ! jetyp(j) .eq. 29
               write(iot,'("$",1p6e13.5)') (egmin(ngei(j)+i),
     &               egmax(ngea(j)+i),fegrp(ngfe(j)+i),i=1,ngrp(j))
               write(iot,'("$  p-type = ",i3,12x,
     &         " # generate weight type, 0:equal, 1: given by data"
     &         )')  jptyp(j)

               if( jptyp(j) .eq. 1 ) then
                 write(iot,'("$",1p6e13.5)')
     &                                    (prw(ngpi(j)+i),i=1,ngrp(j))
               end if
             end if

           else

               write(iot,'("       ne = ",i4,11x,
     &         " # number of energy and weight"/27x,
     &         " #  data = ( emin(i), emax(i), w(i), i = 1, ne )")')
     &         ngrp(j)

             if( jetyp(j) .eq. 28 ) then
               write(iot,'(1p6e13.5)') (egmin(ngei(j)+i),
     &               egmax(ngea(j)+i),fegrp(ngfe(j)+i),i=1,ngrp(j))

             else  ! jetyp(j) .eq. 29
               write(iot,'(1p6e13.5)') (egmin(ngei(j)+i),
     &               egmax(ngea(j)+i),fegrp(ngfe(j)+i),i=1,ngrp(j))
               write(iot,'("   p-type = ",i3,12x,
     &         " # generate weight type, 0:equal, 1: given by data"
     &         )')  jptyp(j)

               if( jptyp(j) .eq. 1 ) then
                 write(iot,'(1p6e13.5)') (prw(ngpi(j)+i),i=1,ngrp(j))
               end if

             end if

           end if

*-----------------------------------------------------------------------
         else if( jetyp(j) .eq. 20 ) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # energy distribution given by tally output")') jetyp(j)

            asfil = '  # file name of tally output'
            msfile = max( 13, lsfile(j) )
            write(iot,'("     file =  ",100a1)')
     &                  ( sfile(j)(i:i), i = 1, msfile ),
     &                  ( asfil(i:i), i = 1, 29 )

            write(iot,'("#      ne = ",i4,11x,
     &      " # number of energy bin"/
     &      "#  e-lower      e-upper      fluence")')
     &      ngrp(j)

            do idata = 1, ngrp(j)
             write(iot,'("#",1p3e13.5)')
     &              egmin(ngei(j)+idata),egmax(ngea(j)+idata)
     &              ,fegrp(ngfe(j)+idata)
            end do

*-----------------------------------------------------------------------
         else if( jetyp(j) .eq. 25 .or. jetyp(j) .eq. 26) then

            write(iot,'("   e-type = ",i3,12x,
     &      " # Cosmic-ray source")') jetyp(j)

            write(iot,'("      eg1 = ",1p1g12.5,3x,
     &      " # minimum cutoff energy [MeV]")') seg1(j)

            write(iot,'("      eg2 = ",1p1g12.5,3x,
     &      " # maximum cutoff energy [MeV]")') seg2(j)

            write(iot,'("       nm = ",i4,11x,
     &      " # number of energy mesh"/27x,
     &      " # nm>0: lin, nm<0: log mesh and interpolation")')
     &      ngrp(j) * ngll(j)

            if(jetyp(j) .eq. 25) then
             write(iot,'("#   e-type = 21",12x,
     &       " # energy distribution given by data")')
            else
             write(iot,'("#   e-type = 24",12x,
     &       " # energy distribution given by data")')
            endif

            write(iot,'("#      ne = ",i4,11x,
     &      " # number of energy bin"/
     &      "#  e-lower      fluence(/cm2/s/(MeV/n))")')
     &      ngrp(j) * ngll(j)

            do ig = 1, ngrp(j)
             write(iot,'("#",1p2e13.5)')
     &       egmin(ngei(j)+ig),
     &       fegrp(ngfe(j)+ig)/(egmax(ngea(j)+ig)-egmin(ngei(j)+ig))
            end do

       write(iot,'("    icenv = ",i3,12x,
     & " # 0>: Terrestiral GCR, D=0: Free-space GCR, ",
     & "<0:Free-space SEP")') icenv(j)

       write(iot,'(" solarmod = ",1p1g12.5,3x,
     & " # Solar activity W-index")') solarmod(j)

       write(iot,'("    rigid = ",1p1g12.5,3x,
     & " # Cut-off rigidity in GV")') rigid(j)

       write(iot,'("  depatom = ",1p1g12.5,3x,
     & " # Atmospheric depth in g/cm2")') depatom(j)

       write(iot,'("  enviorn = ",1p1g12.5,3x,
     & " # Surrounding environment (icenv=2-4) or",
     & " SPE index (icenv=-1)")') environ(j)

*-----------------------------------------------------------------------
         else

            write(io,'(/"*** Error : e-type is wrong in [source]")')
            ErrCha = ''
            MsgID = 'L:21816/R:echoi/F:read00.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,'(/"*** Error : e-type is wrong in [source]")')
            ierr = ierr + 1

         end if

      end if

*-----------------------------------------------------------------------
*     source time parameters
*-----------------------------------------------------------------------

               call echsrst(j,iot,ierr)

*-----------------------------------------------------------------------
*     source region selection
*-----------------------------------------------------------------------

         if( nsrn(j) .gt. 0 .and. jstyp(j) .ne. 12 ) then

               write(iot,'(" # region selection of source")')

                     idsm = iaddress_nsrn(j)
                     jdsm = 0

                     jdsm = jdsm + 1
                     ntrn = idas_nsrn(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_nsrn(idsm+jdsm)

                     kdsm = nsrc(j)
                     ldsm = 0

                     ldsm = ldsm + 1
                     ntrc = idas_nsrc(kdsm+ldsm)
                     ldsm = ldsm + 1
                     mtrc = idas_nsrc(kdsm+ldsm)

                     idas1 = mmmax
                     idas2 = ( idas1 + ntrn - 1 ) * 2 + 1
                     idas3 = idas2 + ntrn

                  call echrg(1,72,0,0,iot,mtrc,idas_nsrc(kdsm+ldsm+1),
     &                       ntrn,mtrn,idas_nsrn(idsm+jdsm+1),
     &                       das(idas1),idas(idas2),
     &                       1,ivl,rvl,
     &                       idas3)

               write(iot,'("    ntmax = ",i7,8x,
     &         " # (D=1000) maximum trial number of region selection"
     &         )') nsmx(j)

         end if

*-----------------------------------------------------------------------
*        transform of source
*-----------------------------------------------------------------------

         if( isort(j,1) .gt. 0 ) then

               if( isort(j,1) .eq. 1 ) then
                     write(iot,'("     trcl = ",i6,9x,
     &               " # transform ID")') isort(j,3)

               else

                  if( isort(j,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rsort(j,k), k = 1, 12 ),
     &               nint( rsort(j,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rsort(j,k), k = 1, 12 ),
     &               nint( rsort(j,13) )

                  end if

               end if

         end if

*-----------------------------------------------------------------------
*        spin
*-----------------------------------------------------------------------

             ssr = ssx(j)**2 + ssy(j)**2 + ssz(j)**2

         if( ssr .gt. 0.d0 ) then

               write(iot,'("       sx = ",1p1g12.5,3x,
     &         " # (D=0.0) spin x component")') ssx(j)
               write(iot,'("       sy = ",1p1g12.5,3x,
     &         " # (D=0.0) spin y component")') ssy(j)
               write(iot,'("       sz = ",1p1g12.5,3x,
     &         " # (D=0.0) spin z component")') ssz(j)

         end if

*-----------------------------------------------------------------------
*        global factor
*-----------------------------------------------------------------------

         if( sfactor(j) .ne. 1.0d0 ) then

               write(iot,'("   factor = ",1p1g12.5,3x,
     &         " # (D=1.0) factor for each source")') sfactor(j)

         end if

*-----------------------------------------------------------------------
*        special for duct
*-----------------------------------------------------------------------

         if( ( jstyp(j) .eq. 1 .or.
     &         jstyp(j) .eq. 4 .or.
     &         jstyp(j) .eq. 2 .or.
     &         jstyp(j) .eq. 5 ) .and.
     &       ( sdom(j) .gt. -10.5 .and. sdom(j) .lt. -9.5 ) ) then

               write(iot,'(/"      dl0 = ",1p1g12.5,3x,
     &         " # distance to entrance")') sdl0(j)
               write(iot,'("      dl1 = ",1p1g12.5,3x,
     &         " # start point of duct")') sdl1(j)
               write(iot,'("      dl2 = ",1p1g12.5,3x,
     &         " # end point of duct")') sdl2(j)
               write(iot,'("      dpf = ",1p1g12.5,3x,
     &         " # (D=dl1**2/dl2**2) prob of end 0<dpf<1")') sdpf(j)

            if( jstyp(j) .eq. 1 .or. jstyp(j) .eq. 4 ) then

               write(iot,'("      drd = ",1p1g12.5,3x,
     &         " # radius of cylinder duct")') sdrd(j)

            else if( jstyp(j) .eq. 2 .or. jstyp(j) .eq. 5 ) then

               write(iot,'("      dxw = ",1p1g12.5,3x,
     &         " # x-width of square duct")') sdxw(j)
               write(iot,'("      dyw = ",1p1g12.5,3x,
     &         " # y-width of square duct")') sdyw(j)

            end if

            if( sdls(j) .gt. 0.0d0 ) then

                  write(iot,'("      dls = ",1p1g12.5,3x,
     &            " # slit positive")') sdls(j)

               if( sdrs(j) .gt. 0.0d0 ) then

                  write(iot,'("      drs = ",1p1g12.5,3x,
     &            " # slit radius")') sdrs(j)

               else

                  write(iot,'("      dxs = ",1p1g12.5,3x,
     &            " # slit x-length")') sdxs(j)
                  write(iot,'("      dys = ",1p1g12.5,3x,
     &            " # slit y-length")') sdys(j)

               end if

            end if

            if( sdsxp(j) .ne. 1.0d0 .or.
     &          sdsxn(j) .ne. 1.0d0 .or.
     &          sdsyp(j) .ne. 1.0d0 .or.
     &          sdsyn(j) .ne. 1.0d0 ) then

               write(iot,'("     dsxp = ",1p1g12.5,3x,
     &         " # x-positive side weight")') sdsxp(j)
               write(iot,'("     dsxn = ",1p1g12.5,3x,
     &         " # x-negative side weight")') sdsxn(j)
               write(iot,'("     dsyp = ",1p1g12.5,3x,
     &         " # y-positive side weight")') sdsyp(j)
               write(iot,'("     dsyn = ",1p1g12.5,3x,
     &         " # y-negative side weight")') sdsyn(j)

            end if

            if( nglp(j) .gt. 0 ) then

               write(iot,'(/"       nl = ",i4,11x,
     &         " # number of length and weight"/
     &         " #  data = ( l(i), lw(i), i = 1, nl ), lmax")')
     &         nglp(j)

               write(iot,'(1p6e13.5)') (slmin(ngli(j)+i),
     &                                  flgrp(ngfl(j)+i),i=1,nglp(j)),
     &                                  slmax(ngla(j)+nglp(j))

            end if

            if( sdebg(j) .ne. 0.0d0 ) then

               write(iot,'(/"     debg = ",1p1g12.5,3x,
     &         " # debug mode on, 0-> off")') sdebg(j)

            end if

         end if

*-----------------------------------------------------------------------
*        ibatch for itall = 4
*-----------------------------------------------------------------------
         if( mstz(26).eq.4 ) then
           if( ibcsn(j).eq.0 ) then
             write(iot,'("   ibatch =   all"11x,
     &       "# (D=all) activated batch number for itall = 4")')

           else if( ibcsn(j) .gt. 0 ) then
             write(iot,'("   ibatch =")',advance='no')
             itt1 = ibcsn(j)/10
             itt2 = mod(ibcsn(j),10)
             if( itt1.gt.0 ) then
               do k = 0, itt1-1
                 write(iot,'(10(i6))') ( ibsor(i+10*k,j),i = 1, 10 )
               end do
             end if
             if ( itt2.ne.0 ) then
               write(fmtt,'("("I0"i6)")') itt2
               write(iot,fmtt) ( ibsor(i+10*itt1,j),i = 1,itt2 )
             end if
           end if
         end if

*-----------------------------------------------------------------------

      end do
      end if

*-----------------------------------------------------------------------
*     material
*-----------------------------------------------------------------------
*              nel  = nint( das(kmatc(i)+ 1) )
*              denh =       das(kmatc(i)+ 2)
*              libh = nint( das(kmatc(i)+ 3) )
*
*              igas = nint( das(kmatc(i)+ 4) )
*              istp = nint( das(kmatc(i)+ 5) )
*              inlb = nint( das(kmatc(i)+ 6) )
*              iplb = nint( das(kmatc(i)+ 7) )
*              ielb = nint( das(kmatc(i)+ 8) )
*              icnd = nint( das(kmatc(i)+ 9) )
*              imts = nint( das(kmatc(i)+10) )
*              know = nint( das(kmatc(i)+11) )
*              iulb = nint( das(kmatc(i)+12) )
*              ihlb = nint( das(kmatc(i)+13) )
*              do ii=1,20
*               idedx(ii) = nint( das(kmatg(i)+13+ii) )
*              enddo
*
*              iz   = nint( das(kmatc(i)+(j-1)*4+34) )
*              ia   = nint( das(kmatc(i)+(j-1)*4+35) )
*              den  =       das(kmatc(i)+(j-1)*4+36)
*              libi = nint( das(kmatc(i)+(j-1)*4+37) )
*
*              ix(1) = nint( das(kmatc(i)+know+(j-1)*3+1) )
*              ix(2) = nint( das(kmatc(i)+know+(j-1)*3+2) )
*              ix(3) = nint( das(kmatc(i)+know+(j-1)*3+3) )
*-----------------------------------------------------------------------

      if( mxmat .gt. 0 .and. ierrg .eq. 0 ) then

         write(iot,'(/"[ Material ]")')


         do i = 1, mxmat0

               nel  = nint( das_kmatc(kmatc(i)+1) )
               denh =       das_kmatc(kmatc(i)+2)
               libh = nint( das_kmatc(kmatc(i)+3) )

                     hs = ' '
                     j = idmn(i)

                  if( mstz(16) .le. 1 ) then

                     if( j .lt. 10 ) then

                        write(hs,'("MAT[ ",i1," ]")') j
                        k = 9

                     else if( j .lt. 100 ) then

                        write(hs,'("MAT[ ",i2," ]")') j
                        k = 10

                     else if( j .lt. 1000 ) then

                        write(hs,'("MAT[ ",i3," ]")') j
                        k = 11

                     else if( j .lt. 10000 ) then

                        write(hs,'("MAT[ ",i4," ]")') j
                        k = 12

                     else if( j .lt. 100000 ) then

                        write(hs,'("MAT[ ",i5," ]")') j
                        k = 13

                     else if( j .lt. 1000000 ) then

                        write(hs,'("MAT[ ",i6," ]")') j
                        k = 14

                     end if

                  else

                     k = 6

                     if( j .lt. 10 ) then

                        write(hs,'("m",i1)') j

                     else if( j .lt. 100 ) then

                        write(hs,'("m",i2)') j

                     else if( j .lt. 1000 ) then

                        write(hs,'("m",i3)') j

                     else if( j .lt. 10000 ) then

                        write(hs,'("m",i4)') j

                     else if( j .lt. 100000 ) then

                        write(hs,'("m",i5)') j
                        k = 7

                     else if( j .lt. 1000000 ) then

                        write(hs,'("m",i6)') j
                        k = 8

                     end if

                  end if


*-----------------------------------------------------------------------
*           keyword parameters
*-----------------------------------------------------------------------

               igas = nint( das_kmatc(kmatc(i)+4) )
               istp = nint( das_kmatc(kmatc(i)+5) )
               inlb = nint( das_kmatc(kmatc(i)+6) )
               iplb = nint( das_kmatc(kmatc(i)+7) )
               ielb = nint( das_kmatc(kmatc(i)+8) )
               icnd = nint( das_kmatc(kmatc(i)+9) )
               iulb = nint( das_kmatc(kmatc(i)+12) )
               ihlb = nint( das_kmatc(kmatc(i)+13) )

*-----------------------------------------------------------------------
*           *) idedx <= dedx file read ggm01.f:setmd
*-----------------------------------------------------------------------

            if( igas .ne. 0 ) then

               ht = ' '
               write(ht,'("gas=",i4)') igas
               call chcomp(ht,1,8,l)
               hs(k+2:k+l+1)=ht
               k = k + l + 2

            end if

            if( istp .ne. 0 ) then

               ht = ' '
               write(ht,'("estep=",i4)') istp
               call chcomp(ht,1,10,l)
               hs(k+2:k+l+1)=ht
               k = k + l + 2

            end if

            if( inlb .ne. 0 ) then

               ich3 =   inlb / 256**2
               ich2 = ( inlb - ich3 * 256**2 ) / 256
               ich1 =   inlb - ich3 * 256**2 - ich2 * 256

               ht = ' '
               write(ht,'("nlib=",3a1)')
     &         char(ich1),char(ich2),char(ich3)
               call chcomp(ht,1,8,l)
               hs(k+2:k+l+1)=ht
               k = k + l + 2

            end if

            if( iplb .ne. 0 ) then

               ich3 =   iplb / 256**2
               ich2 = ( iplb - ich3 * 256**2 ) / 256
               ich1 =   iplb - ich3 * 256**2 - ich2 * 256

               ht = ' '
               write(ht,'("plib=",3a1)')
     &         char(ich1),char(ich2),char(ich3)
               call chcomp(ht,1,8,l)
               hs(k+2:k+l+1)=ht
               k = k + l + 2

            end if

            if( ielb .ne. 0 ) then

               ich3 =   ielb / 256**2
               ich2 = ( ielb - ich3 * 256**2 ) / 256
               ich1 =   ielb - ich3 * 256**2 - ich2 * 256

               ht = ' '
               write(ht,'("elib=",3a1)')
     &         char(ich1),char(ich2),char(ich3)
               call chcomp(ht,1,8,l)
               hs(k+2:k+l+1)=ht
               k = k + l + 2

            end if

            if( iulb .ne. 0 ) then

               ich3 =   iulb / 256**2
               ich2 = ( iulb - ich3 * 256**2 ) / 256
               ich1 =   iulb - ich3 * 256**2 - ich2 * 256

               ht = ' '
               write(ht,'("pnlib=",3a1)')
     &         char(ich1),char(ich2),char(ich3)
               call chcomp(ht,1,9,l)
               hs(k+2:k+l+1)=ht
               k = k + l + 2

            end if

            if( ihlb .ne. 0 ) then

               ich3 =   ihlb / 256**2
               ich2 = ( ihlb - ich3 * 256**2 ) / 256
               ich1 =   ihlb - ich3 * 256**2 - ich2 * 256

               ht = ' '
               write(ht,'("hlib=",3a1)')
     &         char(ich1),char(ich2),char(ich3)
               call chcomp(ht,1,8,l)
               hs(k+2:k+l+1)=ht
               k = k + l + 2

            end if

            if( icnd .ne. 0 ) then

               ht = ' '
               write(ht,'("cond=",i4)') icnd
               call chcomp(ht,1,9,l)
               hs(k+2:k+l+1)=ht
               k = k + l + 2

            end if

               write(iot,'(80a1)') (hs(m:m),m=1,k)

*-----------------------------------------------------------------------

            if( mstz(16) .eq. 0 ) then

               if( denh .ne. 0.0d0 ) then

                     chau = '   1H'

                  if( libh .ne. lib0 ) then

                     ich3 =   libh / 256**2
                     ich2 = ( libh - ich3 * 256**2 ) / 256
                     ich1 =   libh - ich3 * 256**2 - ich2 * 256

                     write(iot,'(8x,a5,".",3a1,3x,1p1e15.7)')
     &               chau, char(ich1),char(ich2),char(ich3), denh

                  else

                     write(iot,'(8x,a5,7x,1p1e15.7)')
     &                           chau, denh

                  end if

               end if

               do j = 1, nel

                  iz   = nint( das_kmatc(kmatc(i)+(j-1)*4+34) )
                  ia   = nint( das_kmatc(kmatc(i)+(j-1)*4+35) )
                  den  =       das_kmatc(kmatc(i)+(j-1)*4+36)
                  libi = nint( das_kmatc(kmatc(i)+(j-1)*4+37) )
                  chat = elmnt(iz)

                  ich3 =   libi / 256**2
                  ich2 = ( libi - ich3 * 256**2 ) / 256
                  ich1 =   libi - ich3 * 256**2 - ich2 * 256

                  if( ia .gt. 0 ) then

                     if( chat(1:1) .eq. ' ' ) then

                        write(chau,'(i4,a1)') ia, chat(2:2)

                     else

                        write(chau,'(i3,a2)') ia, chat(1:2)

                     end if

                  else

                     if( chat(1:1) .eq. ' ' ) then

                        write(chau,'(4x,a1)') chat(2:2)

                     else

                        write(chau,'(3x,a2)') chat(1:2)

                     end if

                  end if

                  if( libi .ne. lib0 ) then

                     write(iot,'(8x,a5,".",3a1,3x,1p1e15.7)')
     &               chau, char(ich1),char(ich2),char(ich3), den

                  else

                     write(iot,'(8x,a5,7x,1p1e15.7)')
     &                           chau, den

                  end if

               end do

*-----------------------------------------------------------------------

            else if( mstz(16) .eq. 1 ) then

               if( denh .ne. 0.0d0 ) then

                     chau = ' H-1'

                     if( libh .ne. lib0 ) then

                        ich3 =   libh / 256**2
                        ich2 = ( libh - ich3 * 256**2 ) / 256
                        ich1 =   libh - ich3 * 256**2 - ich2 * 256

                        write(iot,'(9x,a4,".",3a1,3x,1p1e15.7)')
     &                  chau(1:4), char(ich1),char(ich2),char(ich3),
     &                  denh

                     else

                        write(iot,'(9x,a4,7x,1p1e15.7)')
     &                    chau(1:4), denh

                     end if

               end if

               do j = 1, nel

                  iz   = nint( das_kmatc(kmatc(i)+(j-1)*4+34) )
                  ia   = nint( das_kmatc(kmatc(i)+(j-1)*4+35) )
                  den  =       das_kmatc(kmatc(i)+(j-1)*4+36)
                  libi = nint( das_kmatc(kmatc(i)+(j-1)*4+37) )
                  chau = elmnt(iz)

                  ich3 =   libi / 256**2
                  ich2 = ( libi - ich3 * 256**2 ) / 256
                  ich1 =   libi - ich3 * 256**2 - ich2 * 256

                  if( ia .eq. 0 ) then

                     if( libi .ne. lib0 ) then

                        write(iot,'(11x,a2,
     &                             ".",3a1,3x,1p1e15.7)')
     &               chau(1:2), char(ich1),char(ich2),char(ich3), den

                     else

                        write(iot,'(11x,a2,7x,1p1e15.7)')
     &                    chau(1:2), den

                     end if

                  else if( ia .lt. 10 ) then

                     if( libi .ne. lib0 ) then

                        write(iot,'(9x,a2,"-",i1,
     &                             ".",3a1,3x,1p1e15.7)')
     &                  chau(1:2), ia, char(ich1),char(ich2),char(ich3),
     &                  den

                     else

                        write(iot,'(9x,a2,"-",i1,7x,1p1e15.7)')
     &                    chau(1:2), ia, den

                     end if

                  else if( ia .lt. 100 ) then

                     if( libi .ne. lib0 ) then

                        write(iot,'(8x,a2,"-",i2,
     &                             ".",3a1,3x,1p1e15.7)')
     &                  chau(1:2), ia, char(ich1),char(ich2),char(ich3),
     &                  den

                     else

                        write(iot,'(8x,a2,"-",i2,7x,1p1e15.7)')
     &                    chau(1:2), ia, den

                     end if

                  else

                     if( libi .ne. lib0 ) then

                        write(iot,'(7x,a2,"-",i3,
     &                             ".",3a1,3x,1p1e15.7)')
     &                  chau(1:2), ia, char(ich1),char(ich2),char(ich3),
     &                  den

                     else

                        write(iot,'(7x,a2,"-",i3,7x,1p1e15.7)')
     &                    chau(1:2), ia, den

                     end if

                  end if

               end do

*-----------------------------------------------------------------------

            else if( mstz(16) .eq. 2 ) then

               if( denh .ne. 0.0d0 ) then

                  if( libh .ne. lib0 ) then

                     ich3 =   libh / 256**2
                     ich2 = ( libh - ich3 * 256**2 ) / 256
                     ich1 =   libh - ich3 * 256**2 - ich2 * 256

                     write(iot,'(7x,i6,".",3a1,3x,1p1e15.7)')
     &                  1001, char(ich1),char(ich2),char(ich3), denh

                  else

                     write(iot,'(7x,i6,7x,1p1e15.7)')
     &                          1001, denh

                  end if

               end if

               do j = 1, nel

                  iz   = nint( das_kmatc(kmatc(i)+(j-1)*4+34) )
                  ia   = nint( das_kmatc(kmatc(i)+(j-1)*4+35) )
                  den  =       das_kmatc(kmatc(i)+(j-1)*4+36)
                  libi = nint( das_kmatc(kmatc(i)+(j-1)*4+37) )
                  iza  = iz * 1000 + ia

                  ich3 =   libi / 256**2
                  ich2 = ( libi - ich3 * 256**2 ) / 256
                  ich1 =   libi - ich3 * 256**2 - ich2 * 256

                  if( libi .ne. lib0 ) then

                     write(iot,'(7x,i6,".",3a1,3x,1p1e15.7)')
     &               iza, char(ich1),char(ich2),char(ich3), den

                  else

                     write(iot,'(7x,i6,7x,1p1e15.7)')
     &                          iza, den

                  end if

               end do

            end if

*-----------------------------------------------------------------------
*           s(a,b)
*-----------------------------------------------------------------------

               imts = nint( das_kmatc(kmatc(i)+10) )
               know = nint( das_kmatc(kmatc(i)+11) )

            if( imts .gt. 0 ) then

                        hs = ' '
                        j = idmn(i)

                     k = 8

                     if( j .lt. 10 ) then

                        write(hs,'("mt",i1)') j

                     else if( j .lt. 100 ) then

                        write(hs,'("mt",i2)') j

                     else if( j .lt. 1000 ) then

                        write(hs,'("mt",i3)') j

                     else if( j .lt. 10000 ) then

                        write(hs,'("mt",i4)') j

                     else if( j .lt. 100000 ) then

                        write(hs,'("mt",i5)') j
                        k = 9

                     else if( j .lt. 1000000 ) then

                        write(hs,'("mt",i6)') j
                        k = 10

                     end if

                     l = 0

               do j = 1, imts

                     ix(1) = nint( das_kmatc(kmatc(i)+know+(j-1)*3+1) )
                     ix(2) = nint( das_kmatc(kmatc(i)+know+(j-1)*3+2) )
                     ix(3) = nint( das_kmatc(kmatc(i)+know+(j-1)*3+3) )

                     call zaid(2,ht,ix)

                     write(hs(k:k+10),'(a10)') ht
                     k = k + 11
                     l = l + 1

                  if( j .eq. imts .or. mod(j,5) .eq. 0 ) then

                     write(iot,'(80a1)') (hs(m:m),m=1,k)

                     hs = ' '
                     k = 8
                     l = 0

                  end if

               end do

            end if

*-----------------------------------------------------------------------
            if( index(dedx_filename(i),' ') /= 0) then
             write(iot,'("dedxfile = "a)')
     &        dedx_filename(i)(1:len_trim(dedx_filename(i)))
            endif
*-----------------------------------------------------------------------
            kk = 1
            do while( frac(i,kk) .gt. 0)
             if(kk .eq. 1) write(iot,*) "chem = "
             call form_chemi(chemID, len, ichem(i,kk),ierr)
             call chcptl(chemID,1,len)
             write(iot,*) chemID(1:len), frac(i,kk)
             kk = kk + 1
            enddo
*-----------------------------------------------------------------------
         end do

      end if

*-----------------------------------------------------------------------
*     body
*-----------------------------------------------------------------------

      if( ibody .gt. 0 ) then

               rewind itby
               read(itby) chtit

               call chlngt(chtit,60,i1,i2)

         write(iot,'(/"[ Body ]  ",60a1)') (chtit(k:k),k=1,i2)

               read(itby) ipva

               write(iot,'( " idbg =",i2," ;  ibod =",i2,
     &                                     " ;  naz = ",i2 )')
     &         ipva(2), ipva(3), ipva(4)


            if( ipva(3) .ne. 0 ) then

               write(iot,'( "  num   sym     def")')

            else

               write(iot,'( "sym     def")')

            end if

            do k = 1, ibody

                  read(itby) chbd, ibnum, ibva, ( bval(i), i = 1, ibva )

                  ibck(k) = ibnum

               if( ipva(3) .eq. 0 ) then

                  if( ibva .eq. 4 ) then

                     write(iot,'(a3,3x,4(1p1e15.7))')
     &                    chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 6 ) then

                     write(iot,'(a3,3x,
     &                          2(1p1e15.7)/6x,2(1p1e15.7)
     &                                     /6x,2(1p1e15.7))')
     &                    chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 7 ) then

                     write(iot,'(a3,3x,
     &                          3(1p1e15.7)/6x,4(1p1e15.7))')
     &                    chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 8 ) then

                     write(iot,'(a3,3x,
     &                          3(1p1e15.7)/6x,3(1p1e15.7)
     &                                     /6x,2(1p1e15.7))')
     &                    chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 9 .and. chbd .ne. 'tor' ) then

                     write(iot,'(a3,3x,
     &                          2(1p1e15.7)/6x,2(1p1e15.7)
     &                                     /6x,2(1p1e15.7)
     &                                     /6x,3(1p1e15.7))')
     &                    chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 9 .and. chbd .eq. 'tor' ) then

                     write(iot,'(a3,3x,
     &                          3(1p1e15.7)/6x,3(1p1e15.7)
     &                                     /6x,3(1p1e15.7))')
     &                    chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 10 ) then

                     write(iot,'(a3,3x,
     &                          3(1p1e15.7)/6x,3(1p1e15.7)
     &                                     /6x,4(1p1e15.7))')
     &                    chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 12 ) then

                     write(iot,'(a3,3x,
     &                          3(1p1e15.7)/6x,3(1p1e15.7),
     &                                     /6x,3(1p1e15.7)
     &                                     /6x,3(1p1e15.7))')
     &                    chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 30 ) then

                     write(iot,'(a3,3x,
     &                          3(1p1e15.7),7(/6x,3(1p1e15.7)))')
     &                    chbd, ( bval(i), i = 1, 24 )

                     write(iot,'(7x,6i7)')
     &                    ( nint( bval(i) ), i = 25, 30 )

                  else

                     write(iot,'(a3,3x,
     &                          4(1p1e15.7),5(/6x,4(1p1e15.7)))')
     &                    chbd, ( bval(i), i = 1, ibva )

                  end if

               else if( ipva(3) .gt. 0 ) then

                  if( ibva .eq. 4 ) then

                     write(iot,'(i5,3x,a3,3x,4(1p1e15.7))')
     &                    ibnum, chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 6 ) then

                     write(iot,'(i5,3x,a3,3x,
     &                          2(1p1e15.7)/14x,2(1p1e15.7)
     &                                     /14x,2(1p1e15.7))')
     &                    ibnum, chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 7 ) then

                     write(iot,'(i5,3x,a3,3x,
     &                          3(1p1e15.7)/14x,4(1p1e15.7))')
     &                    ibnum, chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 8 ) then

                     write(iot,'(i5,3x,a3,3x,
     &                          3(1p1e15.7)/14x,3(1p1e15.7)
     &                                     /14x,2(1p1e15.7))')
     &                    ibnum, chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 9 .and. chbd .ne. 'tor' ) then

                     write(iot,'(i5,3x,a3,3x,
     &                          2(1p1e15.7)/14x,2(1p1e15.7)
     &                                     /14x,2(1p1e15.7)
     &                                     /14x,3(1p1e15.7))')
     &                    ibnum, chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 9 .and. chbd .eq. 'tor' ) then

                     write(iot,'(i5,3x,a3,3x,
     &                          3(1p1e15.7)/14x,3(1p1e15.7)
     &                                     /14x,3(1p1e15.7))')
     &                    ibnum, chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 10 ) then

                     write(iot,'(i5,3x,a3,3x,
     &                          3(1p1e15.7)/14x,3(1p1e15.7)
     &                                     /14x,4(1p1e15.7))')
     &                    ibnum, chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 12 ) then

                     write(iot,'(i5,3x,a3,3x,
     &                          3(1p1e15.7)/14x,3(1p1e15.7)
     &                                     /14x,3(1p1e15.7)
     &                                     /14x,3(1p1e15.7))')
     &                    ibnum, chbd, ( bval(i), i = 1, ibva )

                  else if( ibva .eq. 30 ) then

                     write(iot,'(i5,3x,a3,3x,
     &                          3(1p1e15.7),7(/14x,3(1p1e15.7)))')
     &                    ibnum, chbd, ( bval(i), i = 1, 24 )

                     write(iot,'(15x,6i7)')
     &                    ( nint( bval(i) ), i = 25, 30 )

                  else

                     write(iot,'(i5,3x,a3,3x,
     &                          4(1p1e15.7),5(/14x,4(1p1e15.7)))')
     &                    ibnum, chbd, ( bval(i), i = 1, ibva )

                  end if

               end if

            end do

         if( icntl .ne. 2 ) close(itby)

      end if

*-----------------------------------------------------------------------
*     region
*-----------------------------------------------------------------------

               jimpo = 0
               iovid = 0

      if( iregn .gt. 0 .and. icgg .eq. 0 ) then

         write(iot,'(/"[ Region ]")')

               if( mstz(20) .eq. 1 .and.
     &             iimpo .gt. 0 .and.
     &             iimps .eq. 0 )  jimpo = 1

               iuniv = 0

            do i = 1, nrsq

               if( irsq(i) .eq. 5 ) iuniv = 1

            end do

            if( irden .eq. 0 ) then

               if( mstz(19) .eq. 0 ) then

                  if( jimpo .eq. 0 .and. iuniv .eq. 0 ) ityr = 1
                  if( jimpo .ne. 0 .and. iuniv .eq. 0 ) ityr = 2
                  if( jimpo .eq. 0 .and. iuniv .eq. 1 ) ityr = 3
                  if( jimpo .ne. 0 .and. iuniv .eq. 1 ) ityr = 4

               else if( mstz(19) .eq. 1 ) then

                  if( jimpo .eq. 0 .and. iuniv .eq. 0 ) ityr = 5
                  if( jimpo .ne. 0 .and. iuniv .eq. 0 ) ityr = 6
                  if( jimpo .eq. 0 .and. iuniv .eq. 1 ) ityr = 7
                  if( jimpo .ne. 0 .and. iuniv .eq. 1 ) ityr = 8

               end if

            else

               if( mstz(19) .eq. 0 ) then

                  if( jimpo .eq. 0 .and. iuniv .eq. 0 ) ityr = 9
                  if( jimpo .ne. 0 .and. iuniv .eq. 0 ) ityr = 10
                  if( jimpo .eq. 0 .and. iuniv .eq. 1 ) ityr = 11
                  if( jimpo .ne. 0 .and. iuniv .eq. 1 ) ityr = 12

               else if( mstz(19) .eq. 1 ) then

                  if( jimpo .eq. 0 .and. iuniv .eq. 0 ) ityr = 13
                  if( jimpo .ne. 0 .and. iuniv .eq. 0 ) ityr = 14
                  if( jimpo .eq. 0 .and. iuniv .eq. 1 ) ityr = 15
                  if( jimpo .ne. 0 .and. iuniv .eq. 1 ) ityr = 16

               end if

            end if

*-----------------------------------------------------------------------

            if( ityr .eq. 1 ) then

               ild1 = 20

               write(iot,'(
     &         "  num   mat    sym  def")')

            else if( ityr .eq. 2 ) then

               ild1 = 31

               write(iot,'(
     &         "  num   mat   imp         sym  def")')

            else if( ityr .eq. 3 ) then

               ild1 = 25

               write(iot,'(
     &         "  num   mat  uni    sym  def")')

            else if( ityr .eq. 4 ) then

               ild1 = 36

               write(iot,'(
     &         "  num   mat  uni   imp         sym  def")')

            else if( ityr .eq. 5 ) then

               ild1 = 31

               write(iot,'(
     &         "  num   mat   vol         sym  def")')

            else if( ityr .eq. 6 ) then

               ild1 = 44

               write(iot,'(
     &         "  num   mat   imp          vol         ",
     &         "sym  def")')

            else if( ityr .eq. 7 ) then

               ild1 = 36

               write(iot,'(
     &         "  num   mat  uni   vol         sym  def")')

            else if( ityr .eq. 8 ) then

               ild1 = 55

               write(iot,'(
     &         "  num   mat  uni   imp          vol         ",
     &         "sym  def")')

            else if( ityr .eq. 9 ) then

               ild1 = 31

               write(iot,'(
     &         "  num   mat   den         sym  def")')

            else if( ityr .eq. 10 ) then

               ild1 = 44

               write(iot,'(
     &         "  num   mat   den          imp         ",
     &         "sym  def")')

            else if( ityr .eq. 11 ) then

               ild1 = 36

               write(iot,'(
     &         "  num   mat  uni   den         sym  def")')

            else if( ityr .eq. 12 ) then

               ild1 = 55

               write(iot,'(
     &         "  num   mat  uni   den          imp         ",
     &         "sym  def")')

            else if( ityr .eq. 13 ) then

               ild1 = 44

               write(iot,'(
     &         "  num   mat   den          vol         ",
     &         "sym  def")')

            else if( ityr .eq. 14 ) then

               ild1 = 57

               write(iot,'(
     &         "  num   mat   den          imp          ",
     &         "vol         sym  def")')

            else if( ityr .eq. 15 ) then

               ild1 = 55

               write(iot,'(
     &         "  num   mat  uni   den          vol         ",
     &         "sym  def")')

            else if( ityr .eq. 16 ) then

               ild1 = 68

               write(iot,'(
     &         "  num   mat  uni   den          imp          ",
     &         "vol         sym  def")')

            end if

               ild0 = 80 - ild1

*-----------------------------------------------------------------------

            mcmx = ( mdas / 2 - 1 ) * 8 + 1
       if(mcmx.lt.0.0) then ! T.Sato 2020/09/17
         write(*,'("mdas is too large, you have to activate ",
     &   "integer*8 option in param.inc")')
         return
       endif
               call moddas_allocate_cha(MAX_NUM_CHRG, chrg)
            mci = 0

            rewind iod

            ioe = 22
            open(ioe,status='scratch',form='unformatted')

*-----------------------------------------------------------------------

         irerr = 0

         do i = 1, iregn

               rewind ioe
               read(iod) (chrg(mci+k:mci+k),k=1,ichl(i))
               chrg(mci+ichl(i)+1:mci+ichmx) = ' '

               if( idmg(i) .eq. -1 ) iovid = iovid + 1

*-----------------------------------------------------------------------

                  ilrm = ichl(i)
                  isqd = 0
                  isrm = 0

  430             isqd = isqd + 1

                  if( ilrm .le. ild0 ) then

                     if( isqd .eq. 1 ) then

                           ildf = ilrm

cFURUTA20131226 optimization bug? in gfortran 4.8 can be fixed
cFURUTA20140213 but it does not work in BX900
c$$$                        chdf(1:ilrm)=chrg(mci+1+isrm:mci+ilrm+isrm)
                        do k = 1, ilrm

                           chdf(k:k) = chrg(mci+k+isrm:mci+k+isrm)

                        end do

                     else

                           write(ioe) ilrm
                           write(ioe)
     &                     (chrg(mci+k+isrm:mci+k+isrm),k=1,ilrm)

                     end if

                  else

                     do k = ild0, 1, -1

                        if( chrg(mci+k+isrm:mci+k+isrm) .eq. ' ' )
     &                  goto 420

                     end do

  420                k2 = k

                     if( isqd .eq. 1 ) then

                           ildf = k2 - 1

cFURUTA20131226 optimization bug? in gfortran 4.8 can be fixed
cFURUTA20140213 but it does not work in BX900
c$$$                        chdf(1:k2-1)=chrg(mci+1+isrm:mci+k2-1+isrm)
                        do k = 1, k2 - 1

                           chdf(k:k) = chrg(mci+k+isrm:mci+k+isrm)

                        end do

                     else

                           write(ioe) k2 - 1
                           write(ioe)
     &                     (chrg(mci+k+isrm:mci+k+isrm),k=1,k2-1)

                     end if

                     isrm = isrm + k2
                     ilrm = ilrm - k2

                     goto 430

                  end if

*-----------------------------------------------------------------------

            if( ityr .eq. 1 ) then

               write(iot,'(i5,1x,i5,4x,a3,2x,200a1)')
     &                         idrg(i), idmg(i), chsm(i),
     &                         (chdf(j:j),j=1,ildf)

            else if( ityr .eq. 2 ) then

               if( dimp(i) .ne. 0.0d0 ) then

                  write(iot,'(i5,1x,i5,1x,1p1g13.6,1x,a3,2x,200a1)')
     &                       idrg(i), idmg(i), dimp(i), chsm(i),
     &                       (chdf(j:j),j=1,ildf)

               else

                  write(iot,'(i5,1x,i5,1x,f9.5,5x,a3,2x,200a1)')
     &                       idrg(i), idmg(i), dimp(i), chsm(i),
     &                       (chdf(j:j),j=1,ildf)

               end if

            else if( ityr .eq. 3 ) then

               write(iot,'(i5,1x,i5,i5,4x,a3,2x,200a1)')
     &                         idrg(i), idmg(i), iuni(i), chsm(i),
     &                         (chdf(j:j),j=1,ildf)

            else if( ityr .eq. 4 ) then

               if( dimp(i) .ne. 0.0d0 ) then

                  write(iot,'(i5,1x,i5,i5,1x,1p1g13.6,1x,a3,2x,200a1)')
     &                  idrg(i), idmg(i), iuni(i), dimp(i), chsm(i),
     &                  (chdf(j:j),j=1,ildf)

               else

                  write(iot,'(i5,1x,i5,i5,1x,f9.5,5x,a3,2x,200a1)')
     &                  idrg(i), idmg(i), iuni(i), dimp(i), chsm(i),
     &                  (chdf(j:j),j=1,ildf)

               end if

            else if( ityr .eq. 5 ) then

                  write(iot,'(i5,1x,i5,1x,1p1g13.6,1x,a3,2x,200a1)')
     &                       idrg(i), idmg(i), dvol(i), chsm(i),
     &                       (chdf(j:j),j=1,ildf)

            else if( ityr .eq. 6 ) then

                  write(iot,'(i5,1x,i5,1x,2(1p1g13.6),1x,a3,2x,200a1)')
     &                       idrg(i), idmg(i), dimp(i),
     &                       dvol(i), chsm(i),
     &                       (chdf(j:j),j=1,ildf)

            else if( ityr .eq. 7 ) then

                  write(iot,'(i5,1x,i5,i5,1x,1p1g13.6,1x,a3,2x,200a1)')
     &                    idrg(i), idmg(i), iuni(i), dvol(i), chsm(i),
     &                    (chdf(j:j),j=1,ildf)

            else if( ityr .eq. 8 ) then

                  write(iot,'(i5,1x,i5,i5,1x,2(1p1g13.6),1x,a3,2x,
     &                        200a1)')
     &                    idrg(i), idmg(i), iuni(i),
     &                    dimp(i), dvol(i), chsm(i),
     &                    (chdf(j:j),j=1,ildf)

            else if( ityr .eq. 9 ) then

               if( deng(i) .ne. 0.0d0 ) then

                  write(iot,'(i5,1x,i5,1x,1p1g13.6,1x,a3,2x,200a1)')
     &                       idrg(i), idmg(i), deng(i), chsm(i),
     &                       (chdf(j:j),j=1,ildf)

               else

                  write(iot,'(i5,1x,i5,1x,f9.5,5x,a3,2x,200a1)')
     &                       idrg(i), idmg(i), deng(i), chsm(i),
     &                       (chdf(j:j),j=1,ildf)

               end if

            else if( ityr .eq. 10 ) then

                  write(iot,'(i5,1x,i5,1x,2(1p1g13.6),1x,a3,2x,200a1)')
     &                       idrg(i), idmg(i), deng(i),
     &                       dimp(i), chsm(i),
     &                       (chdf(j:j),j=1,ildf)

            else if( ityr .eq. 11 ) then

               if( deng(i) .ne. 0.0d0 ) then

                  write(iot,'(i5,1x,i5,i5,1x,1p1g13.6,1x,a3,2x,200a1)')
     &                  idrg(i), idmg(i), iuni(i), deng(i), chsm(i),
     &                  (chdf(j:j),j=1,ildf)

               else

                  write(iot,'(i5,1x,i5,i5,1x,f9.5,5x,a3,2x,200a1)')
     &                  idrg(i), idmg(i), iuni(i), deng(i), chsm(i),
     &                  (chdf(j:j),j=1,ildf)

               end if

            else if( ityr .eq. 12 ) then

                  write(iot,'(i5,1x,i5,i5,1x,2(1p1g13.6),1x,a3,2x,
     &                        200a1)')
     &                    idrg(i), idmg(i), iuni(i),
     &                    deng(i), dimp(i), chsm(i),
     &                    (chdf(j:j),j=1,ildf)

            else if( ityr .eq. 13 ) then

                  write(iot,'(i5,1x,i5,1x,2(1p1g13.6),1x,a3,2x,200a1)')
     &                       idrg(i), idmg(i), deng(i),
     &                       dvol(i), chsm(i),
     &                       (chdf(j:j),j=1,ildf)

            else if( ityr .eq. 14 ) then

                  write(iot,'(i5,1x,i5,1x,3(1p1g13.6),1x,a3,2x,200a1)')
     &                       idrg(i), idmg(i), deng(i),
     &                       dimp(i), dvol(i), chsm(i),
     &                       (chdf(j:j),j=1,ildf)

            else if( ityr .eq. 15 ) then

                  write(iot,'(i5,1x,i5,i5,1x,2(1p1g13.6),1x,a3,2x,
     &                        200a1)')
     &                    idrg(i), idmg(i), iuni(i),
     &                    deng(i), dvol(i), chsm(i),
     &                    (chdf(j:j),j=1,ildf)

            else if( ityr .eq. 16 ) then

                  write(iot,'(i5,1x,i5,i5,1x,3(1p1g13.6),1x,a3,2x,
     &                        200a1)')
     &                    idrg(i), idmg(i), iuni(i),
     &                    deng(i), dimp(i), dvol(i), chsm(i),
     &                    (chdf(j:j),j=1,ildf)

            end if

            if( isqd .gt. 1 ) then

                  rewind ioe

               do k = 2, isqd

                  read(ioe) ildf
                  read(ioe) (chdf(j:j),j=1,ildf)

                  write(iot,'(200a1)') (cblan(j:j),j=1,ild1),
     &                    (chdf(j:j),j=1,ildf)

               end do

            end if

*-----------------------------------------------------------------------
*           check body number
*-----------------------------------------------------------------------

            ic = 0

  450       ic = ic + 1

            if( ic .gt. ichl(i) ) goto 460

               if( chrg(mci+ic:mci+ic) .eq. ' ' ) goto 450

               if( ( chrg(mci+ic:mci+ic) .eq. 'o' .or.
     &               chrg(mci+ic:mci+ic) .eq. 'O' ) .and.
     &             ( chrg(mci+ic+1:mci+ic+1) .eq. 'r' .or.
     &               chrg(mci+ic+1:mci+ic+1) .eq. 'R' ) ) then

                  ic = ic + 1
                  goto 450

               end if

                     call snum(chrg(mci+1:mci+1),ic,ichl(i),ic2,
     &                         cvvv,ierrt)

                     ic = ic2 - 1

                     ibnm = abs( nint( cvvv ) )

                     do j = 1, ibody

                        if( ibnm .eq. ibck(j) ) goto 450

                     end do

                  write(io,'(/"* Error: above body number = [",
     &               i5," ] is not defined in [body] section"/
     &               )') ibnm

                  ErrCha = ''
                  MsgID = 'L:23302/R:echoi/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/"* Error: above body number = [",
     &               i5," ] is not defined in [body] section"/
     &               )') ibnm

                  irerr = irerr + 1

                     goto 450

  460       continue

*-----------------------------------------------------------------------

         end do

               close(ioe)

               if( iovid .eq. 0 ) then

                  write(io,'(/"*** Error in [region] section",
     &                         "    there is no outer void.")')
                  ErrCha = ''
                  MsgID = 'L:23325/R:echoi/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/"*** Error in [region] section",
     &                         "    there is no outer void.")')

                  irerr = irerr + 1

               end if

               if( irerr .gt. 0 ) then

                  write(io,'(/"*** Error in [region] section",
     &                         " listed above,"/
     &                         "    Number of errors is ",i3)') irerr
                  ErrCha = ''
                  MsgID = 'L:23340/R:echoi/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/"*** Error in [region] section",
     &                         " listed above,"/
     &                         "    Number of errors is ",i3)') irerr

                  ierr = ierr + irerr

               end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*     array
*-----------------------------------------------------------------------

      if( llarr .gt. 0 ) then

         write(iot,'(/"[ Array ]")')

               rewind itar

               do j = 1, llarr

                  read(itar,'(i6,200a1)')
     &                 i2, ( chin(i:i), i = 1, i2 )
                  write(iot,'(200a1)') ( chin(i:i), i = 1, i2 )

               end do

               close( itar )

      end if

*-----------------------------------------------------------------------
*     cell
*-----------------------------------------------------------------------

      if( icgg .eq. 1 ) then
         if( ierrg .ne. 0 ) then
            ioq = iot
            iot = io
         end if

            write(iot,'(/"[ Cell ]")')

*-----------------------------------------------------------------------
*        icells = 0, 1, 2  :  no echo for cell
*-----------------------------------------------------------------------

         if( icells .le. 2 ) then

            write(iot,'("$ no echo for icells = 0, 1, 2")')

         else if( icells .gt. 2 ) then

*-----------------------------------------------------------------------

               ild1 = 30
               ild0 = 72 - ild1

               icerr = 0
               iovid = 0

*-----------------------------------------------------------------------

            mcmx = ( mdas / 2 - 1 ) * 8 + 1
       if(mcmx.lt.0.0) then ! T.Sato 2020/09/17
         write(*,'("mdas is too large, you have to activate ",
     &   "integer*8 option in param.inc")')
         return
       endif
               call moddas_allocate_cha(MAX_NUM_CHRG, chrg)
            mci = 0

            rewind iod

            ioe = 22
            open(ioe,status='scratch',form='unformatted')

*-----------------------------------------------------------------------

         do i = 1, igcel

               read(iod) (chrg(mci+k:mci+k),k=1,ichp(i))
               chrg(mci+ichp(i)+1:mci+ichmx) = ' '

*-----------------------------------------------------------------------

               if( ilike(i) .ne. 0 ) then

                  rewind ioe
                  write(ioe) (chrg(mci+k:mci+k),k=ichl(i)+2,ichp(i))

                     ichp(i) = ichp(i) - ichl(i) - 1
                     ichl(i) = 0

                  rewind ioe
                  read(ioe) (chrg(mci+k:mci+k),k=1,ichp(i))
                  chrg(mci+ichp(i)+1:mci+ichmx) = ' '

               end if

*-----------------------------------------------------------------------

               rewind ioe

                  ilrm = ichl(i)
                  ilrm = ichp(i)

                  isqd = 0
                  isrm = 0

  431             isqd = isqd + 1

                     do k = 1, ild0-3

                        if( chrg(mci+k+isrm:mci+k+isrm+2) .eq. ' #(' )
     &                                                       goto 421
                     end do

                     do k = ild0, 2, -1

                        if( chrg(mci+k+isrm-1:mci+k+isrm) .eq. ') '
     &                  .and. k+isrm .lt. ichl(i) ) goto 421

                     end do

                  if( ilrm .le. ild0 ) then

                     if( isqd .eq. 1 ) then

                           ildf = ilrm

cFURUTA20131226 optimization bug? in gfortran 4.8 can be fixed
cFURUTA20140213 but it does not work in BX900
c$$$                        chdf(1:ilrm)=chrg(mci+1+isrm:mci+ilrm+isrm)
                        do k = 1, ilrm

                           chdf(k:k) = chrg(mci+k+isrm:mci+k+isrm)

                        end do

                     else

                           write(ioe) ilrm
                           write(ioe)
     &                     (chrg(mci+k+isrm:mci+k+isrm),k=1,ilrm)

                     end if

                        goto 422

                  end if

                     do k = 1, ild0-2

                        if( chrg(mci+k+isrm:mci+k+isrm+2) .eq. ' : ' )
     &                                                       goto 421
                     end do

                     do k = ild0, 1, -1

                        if( chrg(mci+k+isrm:mci+k+isrm) .eq. ' ' )
     &                  goto 421

                     end do

  421                k2 = k

                     if( isqd .eq. 1 ) then

                           ildf = k2 - 1

cFURUTA20131226 optimization bug? in gfortran 4.8 can be fixed
cFURUTA20140213 but it does not work in BX900
c$$$                        chdf(1:k2-1)=chrg(mci+1+isrm:mci+k2-1+isrm)
                        do k = 1, k2 - 1

                           chdf(k:k) = chrg(mci+k+isrm:mci+k+isrm)

                        end do

                     else

                           write(ioe) k2 - 1
                           write(ioe)
     &                     (chrg(mci+k+isrm:mci+k+isrm),k=1,k2-1)

                     end if

                     isrm = isrm + k2
                     ilrm = ilrm - k2

                     goto 431

  422             continue

*-----------------------------------------------------------------------

            if( idmg(i) .gt. 0 ) then

               if( matadd .ne. 0 ) then

                  matnum = idmg(i)

               else

                  matnum = idmn(idnm(idmg(i)))

               end if

            else

                  matnum = idmg(i)

            end if

                     chsn(1:7) = '       '

                  if( idrg(i) .lt. 10 ) then

                     write(chsn(2:2),'(i1)') idrg(i)

                  else if( idrg(i) .lt. 100 ) then

                     write(chsn(2:3),'(i2)') idrg(i)

                  else if( idrg(i) .lt. 1000 ) then

                     write(chsn(2:4),'(i3)') idrg(i)

                  else if( idrg(i) .lt. 10000 ) then

                     write(chsn(2:5),'(i4)') idrg(i)

                  else if( idrg(i) .lt. 100000 ) then

                     write(chsn(2:6),'(i5)') idrg(i)

                  else if( idrg(i) .lt. 1000000 ) then

                     write(chsn(2:7),'(i6)') idrg(i)

                  else

                     write(chsn(1:7),'(i7)') idrg(i)

                  end if

            if( ilike(i) .ne. 0 ) then

               write(iot,'(a7,3x,"like ",i6," but",5x,200a1)')
     &            chsn, ilike(i),
     &            (chdf(j:j),j=1,ildf)

            else if( idmg(i) .gt. 0 ) then

               write(iot,'(a7,1x,i6,1p1e15.7,2x,200a1)') ! T.Sato 2020/09/17
     &            chsn, matnum, deng(i),
     &            (chdf(j:j),j=1,ildf)


            else

               if( idmg(i) .eq. -1 ) iovid = iovid + 1

               write(iot,'(a7,1x,i6,15x,2x,200a1)')      ! T.Sato 2020/09/17
     &                         chsn, matnum,
     &                         (chdf(j:j),j=1,ildf)

            end if

            if( isqd .gt. 1 ) then

                  rewind ioe

               do k = 2, isqd

                  read(ioe) ildf
                  read(ioe) (chdf(j:j),j=1,ildf)

                  write(iot,'(200a1)') (cblan(j:j),j=1,ild1),
     &                    (chdf(j:j),j=1,ildf)

               end do

            end if

*-----------------------------------------------------------------------
*           check surface number and parentheses
*-----------------------------------------------------------------------

               ic = 0

               img = 0
               ipr = 0
               isf = 0
               ibk = 0
               k   = 0

  451       ic = ic + 1

               if( ic .gt. ichl(i) .and. ipr .ne. 0 ) then

                  write(io,'(/"* Error: in above definition,",
     &               " parentheses are used incorrectly"/ )')
                  ErrCha = ''
                  MsgID = 'L:23650/R:echoi/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/"* Geometry error occurs!,",
     &       " see [cell] section in the input echo (phits.out)"/)')


                  icerr = icerr + 1

               end if

            if( ic .gt. ichl(i) ) goto 461

               ks = k

               k = index('():#',chrg(mci+ic:mci+ic))

! Nais_2024 >>>
               k_all = 0
               if(k == 4) then
                  if(mci + ic + 3 < ichl(i)) then
                    k_all = index(chrg(mci+ic:),'#all ')
                  else
                    k_all = index(chrg(mci+ic:),'#all')
                  end if
                  if(k_all == 1) then
! #all: get cells
                     choutall =' '
                     ichout = 0
                     do i_all = 1,igcel
                       if(i /= i_all) then
                         if(idmg(i_all) == -1) then
                          else
                            if(iuni(i_all) == 0) then
                               ichout=ichout+1
                               write(choutall((ichout-1)*8+1:ichout*8),
     &                                       '("#",i0," ")') idrg(i_all)
                               if(ichout == 8) then
                                 write(iot,'(a,a)') '$ #all: ',choutall
                                 ichout = 0
                                 choutall = ''
                               endif
                             end if
                          end if
                        end if
                      end do
                      if(ichout >= 1) then
                        write(iot,'(a,a)') '$ #all: ',choutall
                        ichout = 0
                        choutall = ''
                      endif
                    end if
                  end if
! Nais_2024 <<<
!
               if( deqn5( chrg(mci+ic:mci+ic) ) ) k = 5

               if( k .eq. 1 ) then

                  ipr = ipr + 1

               else if( k .eq. 2 ) then

                  ipr = ipr - 1

               end if

               if( k .eq. 5 .and. img .eq. 0 ) then

                     img = 1
                     ini = ic
                     if( ks .eq. 4 ) ibk = 1

               end if


! Nais_2024 >>>
! #all
              if( k == 4 .and. k_all == 1) then
                 ibk = 0

                 go to 452

!                   ibnm = 0

                    do j = 1, igcel
                      if( ibnm .eq. idrg(j) ) then
                         ibk = 0
                        goto 452
                      end if
                    end do

                    write(io,'(/"* Error: above cell number = [",
     &                    i6," ] is not defined in [cell] section"/
     &                    )') ibnm
                    ErrCha = ''
                    MsgID = 'L:23745/R:echoi/F:read00.f'
                    call ErrWrite(MsgID, ErrCha)
                    write(jo,'(/"* Error: above cell number = [",
     &                    i6," ] is not defined in [cell] section"/
     &                     )') ibnm

                    icerr = icerr + 1

  452            continue

                 if(mci + ic + 3 < ichl(i)) then
                    ic = ic + 4
                 else
                   ic = ic + 3
                 end if

                 goto 451

              end if
! Nais_2024 <<<

               if( img .eq. 1 .and.
     &           ( ( k .eq. 5 .and. ic .eq. ichl(i) ) .or.
     &             ( k .ne. 5 ) ) ) then

                     img = 0

                  if( k .eq. 5 ) then

                     ifi = ic

                  else

                     ifi = ic - 1

                  end if

                     isf = isf + 1


                     call onum(chrg(1:ifi),ini,ifi,cvvv,ierrt)

                     ibnm = abs( nint( cvvv ) )

                     if( ibk .eq. 0 ) goto 454
                     if( ibk .eq. 1 ) goto 453

               end if

                  goto 451

  453          continue

! Nais_2024 >>>
               if(k == 4 .and. k_all == 0) then
! Nais_2024 <<<

                     do j = 1, igcel

! Nais_2024 >>>
!                        if( ibnm .eq. idrg(j) ) goto 454

                        if( ibnm .eq. idrg(j) ) then
!                          write(*,*) 'k,j,ibnm=',k,j,ibnm
                          go to 454
                        end if
! Nais_2024 <<<

                     end do

                  write(io,'(/"* Error: above cell number = [",
     &               i6," ] is not defined in [cell] section"/
     &               )') ibnm
                  ErrCha = ''
                  MsgID = 'L:23819/R:echoi/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/"* Error: above cell number = [",
     &               i6," ] is not defined in [cell] section"/
     &               )') ibnm

                  icerr = icerr + 1

                  goto 454

! Nais_2024 >>>
               end if
! Nais_2024 <<<

  454          continue

                  ibk = 0
                  goto 451

  461       continue

               if( isf .eq. 0 .and. ilike(i) .eq. 0 ) then

                  write(io,'(/"* Error: in above definition,",
     &               " there is no surface"/)')
                  ErrCha = ''
                  MsgID = 'L:23845/R:echoi/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/"* Geometry error occurs!,",
     &       " see [cell] section in the input echo (phits.out)"/)')

                  icerr = icerr + 1

               end if

*-----------------------------------------------------------------------

         end do

               close(ioe)

               if( iovid .eq. 0 ) then

                  write(io,'(/"*** Error in [cell] section",
     &                         "    there is no outer void.")')
                  ErrCha = ''
                  MsgID = 'L:23865/R:echoi/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/"*** Error in [cell] section",
     &                         "    there is no outer void.")')

                  icerr = icerr + 1

               end if

               if( icerr .gt. 0 ) then

                  write(io,'(/"*** Error in [cell] section",
     &                         " listed above, stop!!"/
     &                         "    Number of errors is ",i3)') icerr
                  ErrCha = ''
                  MsgID = 'L:23880/R:echoi/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/"*** Error in [cell] section",
     &                         " listed above, stop!!"/
     &                         "    Number of errors is ",i3)') icerr

                  ierr = ierr + icerr

               end if

*-----------------------------------------------------------------------

         if( ierrg .ne. 0 ) then
            iot = ioq
         end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*     surface
*-----------------------------------------------------------------------

      if( icgg .eq. 1 ) then

         if( ierrg .ne. 0 ) then
            ioq = iot
            iot = io
         end if

            write(iot,'(/"[ Surface ]")')

*-----------------------------------------------------------------------

            rewind  ioa

            iserr = 0

         do i = 1, igsuf

            read(ioa) idrf, idtr, idsf, igkst, ( bval(j), j = 1, igkst )

                  if( idrf .eq. -1 ) then
                     chsc = '*'
                  else if( idrf .eq. -2 ) then
                     chsc = '+'
                  else
                     chsc = ' '
                  end if

                     chsn(1:1) = chsc

                  if( idsn(i) .lt. 10 ) then

                     write(chsn(2:2),'(i1)') idsn(i)
                     chsn(3:7) = '     '

                  else if( idsn(i) .lt. 100 ) then

                     write(chsn(2:3),'(i2)') idsn(i)
                     chsn(4:7) = '    '

                  else if( idsn(i) .lt. 1000 ) then

                     write(chsn(2:4),'(i3)') idsn(i)
                     chsn(5:7) = '   '

                  else if( idsn(i) .lt. 10000 ) then

                     write(chsn(2:5),'(i4)') idsn(i)
                     chsn(6:7) = '  '

                  else if( idsn(i) .lt. 100000 ) then

                     write(chsn(2:6),'(i5)') idsn(i)
                     chsn(7:7) = ' '

                  else if( idsn(i) .lt. 1000000 ) then

                     write(chsn(2:7),'(i6)') idsn(i)

                  else

                     write(chsn(1:7),'(i7)') idsn(i)

                  end if

               if( idtr .ne. 0 ) then

                  write(chss,'(i6)') idtr

               else

                  chss = '      '

               end if

                  if( igkst .le. 3 ) then

                     write(iot,'(a7,1x,a6,1x,a3,1x,3(1p1e15.7))')
     &                        chsn, chss, chsf(idsf),
     &                        ( bval(j), j = 1, igkst )

                  else if( igkst .le. 6 ) then

                     write(iot,'(a7,1x,a6,1x,a3,1x,3(1p1e15.7)
     &                           /19x,3(1p1e15.7))')
     &                        chsn, chss, chsf(idsf),
     &                        ( bval(j), j = 1, igkst )

                  else if( igkst .le. 9 ) then

                     write(iot,'(a7,1x,a6,1x,a3,1x,3(1p1e15.7)
     &                           /19x,3(1p1e15.7)
     &                           /19x,3(1p1e15.7))')
     &                        chsn, chss, chsf(idsf),
     &                        ( bval(j), j = 1, igkst )

                  else if( igkst .le. 12 ) then

                     write(iot,'(a7,1x,a6,1x,a3,1x,3(1p1e15.7)
     &                           /19x,3(1p1e15.7)
     &                           /19x,3(1p1e15.7)
     &                           /19x,3(1p1e15.7))')
     &                        chsn, chss, chsf(idsf),
     &                        ( bval(j), j = 1, igkst )

                  else

                     write(iot,'(a7,1x,a6,1x,a3,1x,3(1p1e15.7),
     &                           7(/19x,3(1p1e15.7)))')
     &                        chsn, chss, chsf(idsf),
     &                        ( bval(j), j = 1, igkst )

                  end if

*-----------------------------------------------------------------------

               if( idtr .ne. 0 ) then

                  do j = 1, igtrs

                     if( idtr .eq. idtn(j) ) goto 457

                  end do

                  write(io,'(/"* Error: above transform number = [",
     &               i7," ] is not defined in [transform] section"/
     &               )') idtr
                  ErrCha = ''
                  MsgID = 'L:24035/R:echoi/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/"* Error: above transform number = [",
     &               i7," ] is not defined in [transform] section"/
     &               )') idtr

                  iserr = iserr + 1

  457             continue

               end if

         end do

            close( ioa )

*-----------------------------------------------------------------------

               if( iserr .gt. 0 ) then

                  write(io,'(/"*** Error in [surface] section",
     &                         " listed above, stop!!"/
     &                         "    Number of errors is ",i3)') iserr
                  ErrCha = ''
                  MsgID = 'L:24059/R:echoi/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/"*** Error in [surface] section",
     &                         " listed above, stop!!"/
     &                         "    Number of errors is ",i3)') iserr

                  ierr = ierr + iserr

               end if

*-----------------------------------------------------------------------

         if( ierrg .ne. 0 ) then
            iot = ioq
         end if

      end if

*-----------------------------------------------------------------------
*     transform
*-----------------------------------------------------------------------

      if( igtrs .gt. 0 ) then

            write(iot,'(/"[ Transform ]")')

               innm = 0

               rewind  iob

            do i = 1, igtrs

               read(iob) itrs, ( vtrs(j), j = 1, 13 )

                  write(chsn,'(i7)') idtn(i)

                  call chlngt(chsn,7,ic1,ic2)

                  if( idtn(i) .gt. 999999 ) innm = 1

               if( itrs .eq. 0 ) then

                  if( idtn(i) .le. 999999 ) then

                     chtr = 'tr'//chsn(ic1:ic2)
                     chtr(4+ic2-ic1:13) = ' '

                     write(iot,'(a13,3(1p1e15.7),
     &                           3(/13x,3(1p1e15.7)),i5)')
     &               chtr, ( vtrs(j), j = 1, 12 ), nint( vtrs(13) )

                  end if

               else

                  if( idtn(i) .le. 999999 ) then

                     chtr = '*tr'//chsn(ic1:ic2)
                     chtr(5+ic2-ic1:13) = ' '

                     write(iot,'(a13,3(1p1e15.7),
     &                           3(/13x,3(1p1e15.7)),i5)')
     &               chtr, ( vtrs(j), j = 1, 12 ), nint( vtrs(13) )

                  end if

               end if

            end do

*-----------------------------------------------------------------------

            if( innm .eq. 1 ) then

               write(iot,'(/"$ The following transforms are",
     &              " defined in cell, lat, fill, source or tally")')

               rewind  iob

            do i = 1, igtrs

               read(iob) itrs, ( vtrs(j), j = 1, 13 )

                  write(chsn,'(i7)') idtn(i)

                  call chlngt(chsn,7,ic1,ic2)

               if( itrs .eq. 0 ) then

                  if( idtn(i) .gt. 999999 ) then

                     chtr = '$ tr'//chsn(ic1:ic2)
                     chtr(6+ic2-ic1:13) = ' '

                     write(iot,'(a13,3(1p1e15.7),
     &                           3(/"$",12x,3(1p1e15.7)),i5)')
     &               chtr, ( vtrs(j), j = 1, 12 ), nint( vtrs(13) )

                  end if

               else

                  if( idtn(i) .gt. 999999 ) then

                     chtr = '$ *tr'//chsn(ic1:ic2)
                     chtr(7+ic2-ic1:13) = ' '

                     write(iot,'(a13,3(1p1e15.7),
     &                           3(/"$",12x,3(1p1e15.7)),i5)')
     &               chtr, ( vtrs(j), j = 1, 12 ), nint( vtrs(13) )

                  end if

               end if

            end do

            end if

            close( iob )

      end if

*-----------------------------------------------------------------------
*     Delta Ray
*-----------------------------------------------------------------------

      if( mndel .gt. 0 .and. ierrg .eq. 0 ) then

            write(iot,'(/"[ Delta Ray ]")')

            write(iot,'("    reg        del")')

         do i = 1, mndel

            write(iot,'(i7,3x,1p1g15.7)') ndels(i), rdels(i)

         end do

      end if

*-----------------------------------------------------------------------
*     Delta Ray
*-----------------------------------------------------------------------

      if( mntsc .gt. 0 .and. ierrg .eq. 0 ) then

            write(iot,'(/"[ Track Structure ]")')
            mebg = 0
            mwvl = 0
         do i = 1, mntsc
           if(bgets(i).gt.0.d0) mebg = 1
           if(wvets(i).gt.0.d0) mwvl = 1
         enddo
         if(mebg.eq.0 .and. mwvl.eq.0) then
            write(iot,'("    reg    mID")')
         else
            write(iot,'("    reg    mID       eBG       wvl")')
         endif

         do i = 1, mntsc
           if(mebg.eq.0 .and. mwvl.eq.0) then
            write(iot,'(2i7)') ntsc(i), ktsc(i)
           else
            write(iot,'(2i7,2f10.4)')  ntsc(i), ktsc(i),
     &                                bgets(i), wvets(i)
           endif
         end do

      end if

*-----------------------------------------------------------------------
*     volume
*-----------------------------------------------------------------------

      if( mnvol .gt. 0 .and. ( mstz(19) .eq. 0 .or. icgg .eq. 1 ) .and.
     &    ierrg .eq. 0 ) then

            write(iot,'(/"[ Volume ]")')

            write(iot,'("    reg        vol")')

         do i = 1, mnvol

            write(iot,'(i7,3x,1p1g15.7)') nvols(i), rvols(i)

         end do

      end if

*-----------------------------------------------------------------------
*     temperature
*-----------------------------------------------------------------------

      if( mntmp .gt. 0 .and. ierrg .eq. 0 ) then

            write(iot,'(/"[ Temperature ]")')

            write(iot,'("    reg        tmp")')

         do i = 1, mntmp

            write(iot,'(i7,3x,1p1g15.7)') ntmps(i), rtmps(i)

         end do

      end if

*-----------------------------------------------------------------------
*     brems bias
*-----------------------------------------------------------------------

      if( mnbrs .gt. 0 .and. ierrg .eq. 0 ) then

            write(iot,'(/"[ Brems  Bias ]")')

            write(iot,'(" material = ",i4)')
     &         mnbrs

               k = 0

            do j = 1, mnbrs

               k = k + 1

               jmat(k) = mbbrs(j)

               if( k / 10 * 10 .eq. k .or. j .eq. mnbrs ) then

                  write(iot,'(11x,10(i5))') ( jmat(i),i = 1, k )

                  k = 0

               end if

            end do

                  write(iot,'("    num        bias")')

                     nnm = 0

               do i = 1, 49

                     nnm = nnm + 1

                  if( nnm .eq. 2 .and.
     &                cbrem(i) .ne. cbrem(i-1) ) then

                        write(iot,'(5x,i2,4x,1p1g15.7)') i-1, cbrem(i-1)
                        nnm = 1

                  else if( nnm .gt. 2 .and.
     &                     cbrem(i) .ne. cbrem(i-1) ) then

                     if( nnm .eq. 3 ) then

                        write(iot,'(5x,i2,4x,1p1g15.7)') i-2, cbrem(i-2)
                        write(iot,'(5x,i2,4x,1p1g15.7)') i-1, cbrem(i-1)

                     else

                        j = i - nnm + 1

                        write(iot,'(2x,"{",i2,"-",i2,"}",
     &                              2x,1p1g15.7)') j, i-1, cbrem(i-1)

                     end if

                        nnm = 1

                  end if

                  if( i .eq. 49 ) then

                     if( nnm .le. 3 ) then

                        do j = 1, nnm

                           write(iot,'(5x,i2,4x,1p1g15.7)')
     &                           i-j+1, cbrem(i-j+1)

                        end do

                     else

                        j = i - nnm + 1

                        write(iot,'(2x,"{",i2,"-",i2,"}",
     &                              2x,1p1g15.7)') j, i, cbrem(i)

                     end if

                  end if

               end do

      end if

*-----------------------------------------------------------------------
*     photon weight
*-----------------------------------------------------------------------

      if( mnpwt .gt. 0 .and. ierrg .eq. 0 ) then

            write(iot,'(/"[ Photon  Weight ]")')

            write(iot,'("    reg        pwt")')

         do i = 1, mnpwt

            write(iot,'(i7,3x,1p1g15.7)') npwts(i), rpwts(i)

         end do

      end if

*-----------------------------------------------------------------------
*     importance
*-----------------------------------------------------------------------

      if( jimpo .ne. 1 .and. ierrg .eq. 0 ) then

         if( iimpo .gt. 0 .and. iimps .eq. 0 ) then

                  write(iot,'(/"[ Importance ]")')

                  write(iot,'("   part = ",a8)') pname(20)
                  write(iot,'("    reg        imp")')

               do i = 1, iregn

                  write(iot,'(i7,3x,1p1e15.7)') idrg(i), dimp(i)

               end do

*-----------------------------------------------------------------------

         else if( iimpn .gt. 0 ) then

            do j = 1, iimpn

                  write(iot,'(/"[ Importance ]")')

                  nump = mnimp(j,20)

               if( nump .lt. 19 ) then

                  do k = 1, nump

                     uname(k) = pname(mnimp(j,k))(1:8)

                  end do

               else

                     nump = 1
                     uname(1) = pname(20)(1:8)

               end if

                     write(iot,'("   part = ",a8,19(1x,a8))')
     &               ( uname(i), i = 1, nump )

*-----------------------------------------------------------------------

                     igm = ( mmmax - 1 ) * 2 + 1
                     lngmax = 0

                     idsm = inimc(j)
                     jdsm = 0

               do m = 1, mnimp(j,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inimc(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inimc(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     call echrg2(mtrn,idas_inimc(kssm)
     &                           ,dum1,lng1,icmb,igm)

                     if( lng1 .gt. lngmax ) lngmax = lng1

               end do

                     lblk = lngmax
                     dum2(1:6) = '   reg'

                   do i = 1, lblk
                     dum2(i+6:i+6) = ' '
                   end do

                     dum2(lblk+6+1:lblk+6+7) =
     &                          '    imp'

                     write(iot,'(600a1)') (dum2(i:i),i=1,lblk+6+7)

*-----------------------------------------------------------------------

                     idsm = inimc(j)
                     jdsm = 0

                     kdsm = kfimp(j)

               do m = 1, mnimp(j,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inimc(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inimc(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     rimp = das_kfimp(kdsm-1+m)

                     call echrg2(mtrn,idas_inimc(kssm)
     &                           ,dum1,lng1,icmb,igm)

                     dum2(1:3) = '   '

                  do i = 1, lng1
                     dum2(i+3:i+3) = dum1(i:i)
                  end do

                     lngb = lngmax - lng1

                  do i = 1, lngb
                     dum2(i+lng1+3:i+lng1+3) = ' '
                  end do

                     write(dum2(lngb+lng1+4:lngb+lng1+22),
     &                         '(3x,1p1e15.7)') rimp

                     write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+22)

               end do

            end do

         end if

      end if

*-----------------------------------------------------------------------
*     Data Max
*-----------------------------------------------------------------------

      if(  indmm .gt. 0 .and. ierrg .eq. 0 ) then

         do i = 1, indmm

                  write(iot,'(/"[ Data Max ]")')

               write(iot,'("   part = ",6A10)')
     &         (pname(ipdmm(indmm,ip)),ip=1,ipdpt(indmm))

                  write(iot,'("    mat     nucleus     dmax")')

               do j = 1, ipdnn(i)

                        imm = matdxx(i,j)

                     if( imm .eq. 0 ) then

                        chaum = '   all'

                     else

                        write(chaum,'(i6)') imm

                     end if

                        iaz = nucdxx(i,j)
                        iz = iaz / 1000
                        ia = iaz - iz * 1000

                     if( iaz .eq. 0 ) then

                        chaus = 'all'

                     else if( ia .eq. 0 ) then

                        chaus = elmnt(iz)

                     else

                        call chname(idum,ia,iz,chau)

                        chaus = chau

                     end if

                  write(iot,'(1x,a6,5x,a8,1x,1p1g15.7)')
     &                       chaum, chaus, dmxdxx(i,j)

               end do

*-----------------------------------------------------------------------

         end do

      end if

*-----------------------------------------------------------------------
*     Anatally
*-----------------------------------------------------------------------

      if(  ianat .gt. 0 ) then
         call anatal_echo(iot)
      end if

*-----------------------------------------------------------------------
*     weight window
*-----------------------------------------------------------------------

      if( iwwdp .gt. 0 .and. ierrg .eq. 0 ) then

         do j = 1, iwwdp

                  write(iot,'(/"[ Weight Window ]")')

*-----------------------------------------------------------------------

            if( iwmsh .eq. 1 ) then

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is region-wise")')

            else if( iwmsh .eq. 3 ) then

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

                  call echmty(0,iot,'x',iwxty(j),iwxnm(j),
     &                        rwxdl(j),rwxmi(j),rwxma(j),iwxrg(j)
     &                        ,abs(iwxnm(j)+1),das_iwxrg(iwxrg(j)))
                  call echmty(0,iot,'y',iwyty(j),iwynm(j),
     &                        rwydl(j),rwymi(j),rwyma(j),iwyrg(j)
     &                        ,abs(iwynm(j)+1),das_iwyrg(iwyrg(j)))
                  call echmty(0,iot,'z',iwzty(j),iwznm(j),
     &                        rwzdl(j),rwzmi(j),rwzma(j),iwzrg(j)
     &                        ,abs(iwznm(j)+1),das_iwzrg(iwzrg(j)))

*-----------------------------------------------------------------------

            if( iwwtr(j,1) .gt. 0 ) then

               if( iwwtr(j,1) .eq. 1 ) then

                     write(iot,'("     trcl = ",i6,9x,
     &               " # transformation ID")') iwwtr(j,3)

               else

                  if( iwwtr(j,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rwwtr(j,k), k = 1, 12 ),
     &               nint( rwwtr(j,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rwwtr(j,k), k = 1, 12 ),
     &               nint( rwwtr(j,13) )

                  end if

               end if

            end if

*-----------------------------------------------------------------------

               if( kecho .ne. 0 ) then

                  write(iot,'( "     echo =    1",11x,
     &            " # (D=0); only input, =1; all meshes ")')

               else

                  write(iot,'( "     echo =    0",11x,
     &            " # (D=0); only input, =1; all meshes ")')

               end if

               if( idval .ne. 0 ) then

                  write(iot,'( "     dval =",1p1e14.5,2x,
     &           " # (D=1.0); default WW value for undefined mesh ")')
     &            dvals

               end if

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 4 ) then

                  write(iot,'("     mesh =  tet",11x,
     &            " # mesh type is tetra geometry")')

                  ntrn=idas_inwwc(inwwc(j)+3)

                  call echtet(iot,4+ntrn,idas_inwwc(inwwc(j)+1))

            end if

*-----------------------------------------------------------------------

                  nump = mnwwp(j,20)

               if( nump .lt. 19 ) then

                  do k = 1, nump

                     uname(k) = pname(mnwwp(j,k))(1:8)

                  end do

               else

                     nump = 1
                     uname(1) = pname(20)(1:8)

               end if

                     write(iot,'(/"   part = ",a8,19(1x,a8))')
     &               ( uname(i), i = 1, nump )

*-----------------------------------------------------------------------

                     ien = abs( ienww(j) )

               if( ien .gt. 0 ) then

                  if( ienww(j) .gt. 0 ) then

                        write(iot,'("    eng =",i3)') ien

                  else

                        write(iot,'("    tim =",i3)') ien

                  end if

                        k   = 0

                  do i = 1, ien

                        k = k + 1

                        rmat(k) = eenww(j,i)

                        if( k / 5 * 5 .eq. k .or. k .eq. ien ) then

                           write(iot,'(6x,5(1p1e14.5))')
     &                          ( rmat(l),l = 1, k )

                           k = 0

                        end if

                  end do

               else

                     ien = 1

               end if

*-----------------------------------------------------------------------

               if( iwmsh .eq. 1 ) then

                     igm = ( mmmax - 1 ) * 2 + 1
                     lngmax = 0

                     idsm = inwwc(j)
                     jdsm = 0

                  do m = 1, mnwwp(j,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inwwc(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inwwc(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     call echrg2(mtrn,idas_inwwc(kssm)
     &                           ,dum3,lng1,icmb,igm)

                     if( lng1 .gt. lngmax ) lngmax = lng1

                  end do

                     lblk = lngmax
                     dum1(1:6) = '   reg'

                  do i = 1, lblk
                     dum1(i+6:i+6) = ' '
                  end do

                     lbin = lblk + 6

*-----------------------------------------------------------------------

               else if( iwmsh .eq. 3 ) then

                     dum1(1:6) = '   xyz'
                     lni = 6

                  if( iwxnm(j) .lt. 10 ) then
                     laix = 1
                  else if( iwxnm(j) .lt. 100 ) then
                     laix = 2
                  else if( iwxnm(j) .lt. 1000 ) then
                     laix = 3
                  else if( iwxnm(j) .lt. 10000 ) then
                     laix = 4
                  else if( iwxnm(j) .lt. 100000 ) then
                     laix = 5
                  end if

                  if( iwynm(j) .lt. 10 ) then
                     laiy = 1
                  else if( iwynm(j) .lt. 100 ) then
                     laiy = 2
                  else if( iwynm(j) .lt. 1000 ) then
                     laiy = 3
                  else if( iwynm(j) .lt. 10000 ) then
                     laiy = 4
                  else if( iwynm(j) .lt. 100000 ) then
                     laiy = 5
                  end if

                  if( iwznm(j) .lt. 10 ) then
                     laiz = 1
                  else if( iwznm(j) .lt. 100 ) then
                     laiz = 2
                  else if( iwznm(j) .lt. 1000 ) then
                     laiz = 3
                  else if( iwznm(j) .lt. 10000 ) then
                     laiz = 4
                  else if( iwznm(j) .lt. 100000 ) then
                     laiz = 5
                  end if

                     lbin = laix + laiy + laiz + 3

                  do i = 1, lbin
                     dum1(lni+i:lni+i) = ' '
                  end do

                     lbin = lbin + lni

*-----------------------------------------------------------------------

               else if( iwmsh .eq. 4 ) then

                     lngmax = 8
                     lblk = lngmax
                     dum1(1:6) = '   tet'

                   do i = 1, lblk
                     dum1(i+6:i+6) = ' '
                   end do

                     lbin = lngmax + 6

               end if

*-----------------------------------------------------------------------

               ienm = ( ien - 1 ) / 5 + 1

            do ll = 1, ienm

               li = ( ll - 1 ) * 5 + 1
               lf = min( li + 4, ien )

               ld = lf - li + 1

               lbb = lbin

               do jj = 1, ld

                     nn = li + jj - 1

                  if( nn .lt. 10 ) then

                     dum1(lbb+1:lbb+8) = '      ww'
                     write(dum1(lbb+9:lbb+9),'(i1)') nn
                     dum1(lbb+10:lbb+14) = '     '

                  else

                     dum1(lbb+1:lbb+7) = '     ww'
                     write(dum1(lbb+8:lbb+9),'(i2)') nn
                     dum1(lbb+10:lbb+14) = '     '

                  end if

                     lbb = lbb + 14

               end do

                     write(iot,'(/600a1)') (dum1(i:i),i=1,lbb)

*-----------------------------------------------------------------------

            if( iwmsh .eq. 1 ) then

                     idsm = inwwc(j)
                     jdsm = 0

                     kdsm = kfwwp(j)

               do m = 1, mnwwp(j,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inwwc(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inwwc(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     call echrg2(mtrn,idas_inwwc(kssm)
     &                           ,dum3,lng1,icmb,igm)

                     dum2(1:3) = '   '

                  do i = 1, lng1
                     dum2(i+3:i+3) = dum3(i:i)
                  end do

                     lngb = lngmax - lng1

                  do i = 1, lngb + 2
                     dum2(i+lng1+3:i+lng1+3) = ' '
                  end do

                     lni = lngb + lng1 + 3 + 2 + 1

                  do jj = 1, ld

                     nn = li + jj - 1

                     wwt = das_kfwwp(kdsm+(nn-1)*kvlmax+m-1)

                     write(dum2(lni:lni+14),
     &                         '(1p1e14.5)') wwt

                     lni = lni + 14

                  end do

                     write(iot,'(600a1)') (dum2(i:i),i=1,lni)

               end do

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 3 ) then

            if( kecho .eq. 0 ) then

                     idsm = inwwc(j)
                     ixyz = iwxnm(j) * iwynm(j) * iwznm(j)

               do m = 1, mnwwp(j,0)

                     ivx = idas_inwwc(idsm+(m-1)*3+1-1)
                     ivy = idas_inwwc(idsm+(m-1)*3+2-1)
                     ivz = idas_inwwc(idsm+(m-1)*3+3-1)

                     dum2(1:3) = '  ('
                     lni = 3 + 1

                  if( laix .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') ivx
                  else if( laix .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') ivx
                  else if( laix .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') ivx
                  else if( laix .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') ivx
                  else if( laix .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') ivx
                  end if

                     lni = lni + laix
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( laiy .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') ivy
                  else if( laiy .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') ivy
                  else if( laiy .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') ivy
                  else if( laiy .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') ivy
                  else if( laiy .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') ivy
                  end if

                     lni = lni + laiy
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( laiz .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') ivz
                  else if( laiz .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') ivz
                  else if( laiz .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') ivz
                  else if( laiz .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') ivz
                  else if( laiz .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') ivz
                  end if

                     lni = lni + laiz
                     dum2(lni:lni) = ')'
                     lni = lni + 1

                  do jj = 1, ld

                     nn = li + jj - 1

                     wwt = das_kfwwp(kfwwp(j)+(nn-1)*ixyz+m-1)

                     write(dum2(lni:lni+14),
     &                         '(1p1e14.5)') wwt

                     lni = lni + 14

                  end do

                     write(iot,'(600a1)') (dum2(i:i),i=1,lni)

               end do

*-----------------------------------------------------------------------

            else if( kecho .ne. 0 ) then

                     inx = iwxnm(j)
                     iny = iwynm(j)
                     inz = iwznm(j)
                     ixyz = inx * iny * inz

               do jx = 1, inx
               do jy = 1, iny
               do jz = 1, inz

                     dum2(1:3) = '  ('
                     lni = 3 + 1

                  if( laix .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jx
                  else if( laix .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jx
                  else if( laix .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jx
                  else if( laix .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jx
                  else if( laix .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jx
                  end if

                     lni = lni + laix
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( laiy .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jy
                  else if( laiy .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jy
                  else if( laiy .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jy
                  else if( laiy .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jy
                  else if( laiy .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jy
                  end if

                     lni = lni + laiy
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( laiz .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jz
                  else if( laiz .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jz
                  else if( laiz .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jz
                  else if( laiz .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jz
                  else if( laiz .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jz
                  end if

                     lni = lni + laiz
                     dum2(lni:lni) = ')'
                     lni = lni + 1

                  icf = jz + inz * ( jy - 1 )+ iny * inz * ( jx - 1 )

                  do jj = 1, ld

                     nn = li + jj - 1

                     wwt = das(kgwwp(j)+(nn-1)*ixyz+icf-1)

                     write(dum2(lni:lni+14),
     &                         '(1p1e14.5)') wwt

                     lni = lni + 14

                  end do

                     write(iot,'(600a1)') (dum2(i:i),i=1,lni)

               end do
               end do
               end do

*-----------------------------------------------------------------------

            end if

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 4 ) then

                     jdsm = inwwc(j)
                     kdsm = kfwwp(j)
                     ntrn = idas_inwwc(jdsm+3)

                   do m = 1, mnwwp(j,0)

                     lng1=8
                     write(dum1(1:lng1),'(i8)')
     &                    idas_inwwc(jdsm+4+ntrn+m)

                     dum2(1:3) = '   '

                     do i = 1, lng1
                      dum2(i+3:i+3) = dum1(i:i)
                     end do

                     lngb = lngmax - lng1

                     do i = 1, lngb + 2
                      dum2(i+lng1+3:i+lng1+3) = ' '
                     end do

                     lni = lngb + lng1 + 3 + 2 + 1

                     do jj = 1, ld

                      nn = li + jj - 1

                      wwt = das_kfwwp(kdsm+(nn-1)*mnwwp(j,0)+m-1)

                      write(dum2(lni:lni+14),
     &                     '(1p1e14.5)') wwt

                      lni = lni + 14

                     enddo

                     write(iot,'(600a1)') (dum2(i:i),i=1,lni)

                   enddo

            end if

*-----------------------------------------------------------------------

            end do

         end do

      end if

*-----------------------------------------------------------------------
*     WW Bais
*-----------------------------------------------------------------------

      if( iwbdp .gt. 0 .and. ierrg .eq. 0 ) then

         do j = 1, iwbdp

               if( iwwbias .eq. 0 .or. iwwdp .gt. 0 ) then

                  write(iot,'(/"[ WW Bias ] off")')

               else

                  write(iot,'(/"[ WW Bias ]")')

               end if

            if( iwmsh .eq. 1 ) then

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is region-wise")')

            else if( iwmsh .eq. 3 ) then

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

                  call echmty(0,iot,'x',iwbxty(j),iwbxnm(j),
     &                        rwbxdl(j),rwbxmi(j),rwbxma(j),iwbxrg(j)
     &                        ,abs(iwbxnm(j)+1),das_iwbxrg(iwbxrg(j)))
                  call echmty(0,iot,'y',iwbyty(j),iwbynm(j),
     &                        rwbydl(j),rwbymi(j),rwbyma(j),iwbyrg(j)
     &                        ,abs(iwbynm(j)+1),das_iwbyrg(iwbyrg(j)))
                  call echmty(0,iot,'z',iwbzty(j),iwbznm(j),
     &                        rwbzdl(j),rwbzmi(j),rwbzma(j),iwbzrg(j)
     &                        ,abs(iwbznm(j)+1),das_iwbzrg(iwbzrg(j)))

*-----------------------------------------------------------------------

            if( iwbtr(j,1) .gt. 0 ) then

               if( iwbtr(j,1) .eq. 1 ) then

                     write(iot,'("     trcl = ",i6,9x,
     &               " # transformation ID")') iwbtr(j,3)

               else

                  if( iwbtr(j,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rwbtr(j,k), k = 1, 12 ),
     &               nint( rwbtr(j,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rwbtr(j,k), k = 1, 12 ),
     &               nint( rwbtr(j,13) )

                  end if

               end if

            end if

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 4 ) then

                  write(iot,'("     mesh =  tet",11x,
     &            " # mesh type is tetra geometry")')

                  ntrn=idas_inwbc(inwbc(j)+3)

                  call echtet(iot,4+ntrn,idas_inwbc(inwbc(j)+1))

            end if

*-----------------------------------------------------------------------

                  nump = mnwbp(j,20)

               if( nump .lt. 19 ) then

                  do k = 1, nump

                     uname(k) = pname(mnwbp(j,k))(1:8)

                  end do

               else

                     nump = 1
                     uname(1) = pname(20)(1:8)

               end if

                     write(iot,'("   part = ",a8,19(1x,a8))')
     &               ( uname(i), i = 1, nump )

*-----------------------------------------------------------------------

                     ien = abs( ienwb(j) )

               if( ien .gt. 0 ) then

                  if( ienwb(j) .gt. 0 ) then

                        write(iot,'("    eng =",i3)') ien

                  else

                        write(iot,'("    tim =",i3)') ien

                  end if

                        k   = 0

                  do i = 1, ien

                        k = k + 1

                        rmat(k) = eenwb(j,i)

                        if( k / 5 * 5 .eq. k .or. k .eq. ien ) then

                           write(iot,'(6x,5(1p1e14.5))')
     &                          ( rmat(l),l = 1, k )

                           k = 0

                        end if

                  end do

               else

                     ien = 1

               end if

*-----------------------------------------------------------------------

               if( iwmsh .eq. 1 ) then

                     igm = ( mmmax - 1 ) * 2 + 1
                     lngmax = 0

                     idsm = inwbc(j)
                     jdsm = 0

                  do m = 1, mnwbp(j,0)

                     jdsm = jdsm + 1
!<-20220126murofushi update
!--                     ntrn = idas(idsm+jdsm)
                     ntrn = idas_inwbc(idsm+jdsm)
                     jdsm = jdsm + 1
!<-20220126murofushi update
!--                     mtrn = idas(idsm+jdsm)
                     mtrn = idas_inwbc(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

!<-20220126murofushi update
!--                     call echrg2(mtrn,idas(kssm),dum3,lng1,icmb,igm)
                     call echrg2(mtrn,idas_inwbc(kssm)
     &                           ,dum3,lng1,icmb,igm)

                     if( lng1 .gt. lngmax ) lngmax = lng1

                  end do

                     lblk = lngmax
                     dum1(1:6) = '   reg'

                  do i = 1, lblk
                     dum1(i+6:i+6) = ' '
                  end do

                     lbin = lblk + 6

*-----------------------------------------------------------------------

               else if( iwmsh .eq. 3 ) then

                     dum1(1:6) = '   xyz'
                     lni = 6

                  if( iwbxnm(j) .lt. 10 ) then
                     laix = 1
                  else if( iwbxnm(j) .lt. 100 ) then
                     laix = 2
                  else if( iwbxnm(j) .lt. 1000 ) then
                     laix = 3
                  else if( iwbxnm(j) .lt. 10000 ) then
                     laix = 4
                  else if( iwbxnm(j) .lt. 100000 ) then
                     laix = 5
                  end if

                  if( iwbynm(j) .lt. 10 ) then
                     laiy = 1
                  else if( iwbynm(j) .lt. 100 ) then
                     laiy = 2
                  else if( iwbynm(j) .lt. 1000 ) then
                     laiy = 3
                  else if( iwbynm(j) .lt. 10000 ) then
                     laiy = 4
                  else if( iwbynm(j) .lt. 100000 ) then
                     laiy = 5
                  end if

                  if( iwbznm(j) .lt. 10 ) then
                     laiz = 1
                  else if( iwbznm(j) .lt. 100 ) then
                     laiz = 2
                  else if( iwbznm(j) .lt. 1000 ) then
                     laiz = 3
                  else if( iwbznm(j) .lt. 10000 ) then
                     laiz = 4
                  else if( iwbznm(j) .lt. 100000 ) then
                     laiz = 5
                  end if

                     lbin = laix + laiy + laiz + 3

                  do i = 1, lbin
                     dum1(lni+i:lni+i) = ' '
                  end do

                     lbin = lbin + lni

*-----------------------------------------------------------------------

               else if( iwmsh .eq. 4 ) then

                     lngmax = 8
                     lblk = lngmax
                     dum1(1:6) = '   tet'

                   do i = 1, lblk
                     dum1(i+6:i+6) = ' '
                   end do

                     lbin = lngmax + 6

               end if

*-----------------------------------------------------------------------

               ienm = ( ien - 1 ) / 5 + 1

            do ll = 1, ienm

               li = ( ll - 1 ) * 5 + 1
               lf = min( li + 4, ien )

               ld = lf - li + 1

               lbb = lbin

               do jj = 1, ld

                     nn = li + jj - 1

                  if( nn .lt. 10 ) then

                     dum1(lbb+1:lbb+8) = '     wwb'
                     write(dum1(lbb+9:lbb+9),'(i1)') nn
                     dum1(lbb+10:lbb+13) = '    '

                  else

                     dum1(lbb+1:lbb+7) = '    wwb'
                     write(dum1(lbb+8:lbb+9),'(i2)') nn
                     dum1(lbb+10:lbb+13) = '    '

                  end if

                     lbb = lbb + 13

               end do

                     write(iot,'(/600a1)') (dum1(i:i),i=1,lbb)

*-----------------------------------------------------------------------

            if( iwmsh .eq. 1 ) then

                     idsm = inwbc(j)
                     jdsm = 0

                     kdsm = kfwbp(j)

               do m = 1, mnwbp(j,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inwbc(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inwbc(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     call echrg2(mtrn,idas_inwbc(kssm)
     &                           ,dum3,lng1,icmb,igm)

                     dum2(1:3) = '   '

                  do i = 1, lng1
                     dum2(i+3:i+3) = dum3(i:i)
                  end do

                     lngb = lngmax - lng1

                  do i = 1, lngb + 2
                     dum2(i+lng1+3:i+lng1+3) = ' '
                  end do

                     lni = lngb + lng1 + 3 + 2 + 1

                  do jj = 1, ld

                     nn = li + jj - 1

                     wwt = das_kfwbp(kdsm+(nn-1)*kvlmax+m-1)

                     write(dum2(lni:lni+14),
     &                         '(1p1e14.5)') wwt

                     lni = lni + 14

                  end do

                     write(iot,'(600a1)') (dum2(i:i),i=1,lni)

               end do

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 3 ) then

                     inx = iwbxnm(j)
                     iny = iwbynm(j)
                     inz = iwbznm(j)
                     ixyz = inx * iny * inz

               do jx = 1, inx
               do jy = 1, iny
               do jz = 1, inz

                     dum2(1:3) = '  ('
                     lni = 3 + 1

                  if( laix .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jx
                  else if( laix .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jx
                  else if( laix .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jx
                  else if( laix .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jx
                  else if( laix .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jx
                  end if

                     lni = lni + laix
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( laiy .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jy
                  else if( laiy .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jy
                  else if( laiy .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jy
                  else if( laiy .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jy
                  else if( laiy .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jy
                  end if

                     lni = lni + laiy
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( laiz .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jz
                  else if( laiz .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jz
                  else if( laiz .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jz
                  else if( laiz .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jz
                  else if( laiz .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jz
                  end if

                     lni = lni + laiz
                     dum2(lni:lni) = ')'
                     lni = lni + 1

                  icf = jz + inz * ( jy - 1 )+ iny * inz * ( jx - 1 )

                  do jj = 1, ld

                     nn = li + jj - 1

                     wwt = das_kfwbp(kfwbp(j)+(nn-1)*ixyz+icf-1)

                     write(dum2(lni:lni+14),
     &                         '(1p1e14.5)') wwt

                     lni = lni + 14

                  end do

                     write(iot,'(600a1)') (dum2(i:i),i=1,lni)

               end do
               end do
               end do

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 4 ) then

                     jdsm = inwbc(j)
                     kdsm = kfwbp(j)
                     ntrn = idas_inwbc(jdsm+3)

                   do m = 1, mnwbp(j,0)

                     lng1=8
                     write(dum1(1:lng1),'(i8)')
     &                    idas_inwbc(jdsm+4+ntrn+m)

                     dum2(1:3) = '   '

                     do i = 1, lng1
                      dum2(i+3:i+3) = dum1(i:i)
                     end do

                     lngb = lngmax - lng1

                     do i = 1, lngb + 2
                      dum2(i+lng1+3:i+lng1+3) = ' '
                     end do

                     lni = lngb + lng1 + 3 + 2 + 1

                     do jj = 1, ld

                      nn = li + jj - 1

                      wwt = das_kfwbp(kdsm+(nn-1)*mnwbp(j,0)+m-1)

                      write(dum2(lni:lni+14),
     &                     '(1p1e14.5)') wwt

                      lni = lni + 14

                     enddo

                     write(iot,'(600a1)') (dum2(i:i),i=1,lni)

                   enddo

            end if

*-----------------------------------------------------------------------

            end do

         end do

      end if

*-----------------------------------------------------------------------
*     iwwbias = 1, iwbdp = 0 and iwwdp > 0
*-----------------------------------------------------------------------

      if( iwwbias .ne. 0 .and. iwbdp .eq. 0 .and. iwwdp .ne.0 .and.
     &    ierrg .eq. 0 ) then

         do j = 1, iwwdp

                  write(iot,'(/"[ WW Bias ] off")')

                  nump = mnwwp(j,20)

               if( nump .lt. 19 ) then

                  do k = 1, nump

                     uname(k) = pname(mnwwp(j,k))(1:8)

                  end do

               else

                     nump = 1
                     uname(1) = pname(20)(1:8)

               end if

                     write(iot,'("   part = ",a8,19(1x,a8))')
     &               ( uname(i), i = 1, nump )

*-----------------------------------------------------------------------

                     ien = abs( ienww(j) )

               if( ien .gt. 0 ) then

                  if( ienww(j) .gt. 0 ) then

                        write(iot,'("    eng =",i3)') ien

                  else

                        write(iot,'("    tim =",i3)') ien

                  end if

                        k   = 0

                  do i = 1, ien

                        k = k + 1

                        rmat(k) = eenww(j,i)

                        if( k / 5 * 5 .eq. k .or. k .eq. ien ) then

                           write(iot,'(6x,5(1p1e14.5))')
     &                          ( rmat(l),l = 1, k )

                           k = 0

                        end if

                  end do

               else

                     ien = 1

               end if

*-----------------------------------------------------------------------

               if( iwmsh .eq. 1 ) then

                     igm = ( mmmax - 1 ) * 2 + 1
                     lngmax = 0

                     idsm = inwwc(j)
                     jdsm = 0

                  do m = 1, mnwwp(j,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inwwc(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inwwc(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     call echrg2(mtrn,idas_inwwc(kssm)
     &                           ,dum3,lng1,icmb,igm)

                     if( lng1 .gt. lngmax ) lngmax = lng1

                  end do

                     lblk = lngmax
                     dum1(1:6) = '   reg'

                   do i = 1, lblk
                     dum1(i+6:i+6) = ' '
                   end do

                     lbin = lblk + 6

*-----------------------------------------------------------------------

               else if( iwmsh .eq. 3 ) then

                     dum1(1:6) = '   xyz'
                     lni = 6

                  if( iwxnm(j) .lt. 10 ) then
                     laix = 1
                  else if( iwxnm(j) .lt. 100 ) then
                     laix = 2
                  else if( iwxnm(j) .lt. 1000 ) then
                     laix = 3
                  else if( iwxnm(j) .lt. 10000 ) then
                     laix = 4
                  else if( iwxnm(j) .lt. 100000 ) then
                     laix = 5
                  end if

                  if( iwynm(j) .lt. 10 ) then
                     laiy = 1
                  else if( iwynm(j) .lt. 100 ) then
                     laiy = 2
                  else if( iwynm(j) .lt. 1000 ) then
                     laiy = 3
                  else if( iwynm(j) .lt. 10000 ) then
                     laiy = 4
                  else if( iwynm(j) .lt. 100000 ) then
                     laiy = 5
                  end if

                  if( iwznm(j) .lt. 10 ) then
                     laiz = 1
                  else if( iwznm(j) .lt. 100 ) then
                     laiz = 2
                  else if( iwznm(j) .lt. 1000 ) then
                     laiz = 3
                  else if( iwznm(j) .lt. 10000 ) then
                     laiz = 4
                  else if( iwznm(j) .lt. 100000 ) then
                     laiz = 5
                  end if

                     lbin = laix + laiy + laiz + 3

                  do i = 1, lbin
                     dum1(lni+i:lni+i) = ' '
                  end do

                     lbin = lbin + lni

*-----------------------------------------------------------------------

               else if( iwmsh .eq. 4 ) then

                     lngmax = 8
                     lblk = lngmax
                     dum1(1:6) = '   tet'

                   do i = 1, lblk
                     dum1(i+6:i+6) = ' '
                   end do

                     lbin = lngmax + 6

               end if

*-----------------------------------------------------------------------

               ienm = ( ien - 1 ) / 5 + 1

            do ll = 1, ienm

               li = ( ll - 1 ) * 5 + 1
               lf = min( li + 4, ien )

               ld = lf - li + 1

               lbb = lbin

               do jj = 1, ld

                     nn = li + jj - 1

                  if( nn .lt. 10 ) then

                     dum1(lbb+1:lbb+8) = '     wwb'
                     write(dum1(lbb+9:lbb+9),'(i1)') nn
                     dum1(lbb+10:lbb+13) = '    '

                  else

                     dum1(lbb+1:lbb+7) = '    wwb'
                     write(dum1(lbb+8:lbb+9),'(i2)') nn
                     dum1(lbb+10:lbb+13) = '    '

                  end if

                     lbb = lbb + 13

               end do

                     write(iot,'(/600a1)') (dum1(i:i),i=1,lbb)

*-----------------------------------------------------------------------

            if( iwmsh .eq. 1 ) then

                     idsm = inwwc(j)
                     jdsm = 0

                     kdsm = kfwwp(j)

               do m = 1, mnwwp(j,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inwwc(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inwwc(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     call echrg2(mtrn,idas_inwwc(kssm)
     &                           ,dum3,lng1,icmb,igm)

                     dum2(1:3) = '   '

                  do i = 1, lng1
                     dum2(i+3:i+3) = dum3(i:i)
                  end do

                     lngb = lngmax - lng1

                  do i = 1, lngb + 2
                     dum2(i+lng1+3:i+lng1+3) = ' '
                  end do

                     lni = lngb + lng1 + 3 + 2 + 1

                  do jj = 1, ld

                     nn = li + jj - 1

                     wwt = 1.0d0

                     write(dum2(lni:lni+14),
     &                         '(1p1e14.5)') wwt

                     lni = lni + 14

                  end do

                     write(iot,'(600a1)') (dum2(i:i),i=1,lni)

               end do

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 3 ) then

                     inx = iwxnm(j)
                     iny = iwynm(j)
                     inz = iwznm(j)
                     ixyz = inx * iny * inz

               do jx = 1, inx
               do jy = 1, iny
               do jz = 1, inz

                     dum2(1:3) = '  ('
                     lni = 3 + 1

                  if( laix .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jx
                  else if( laix .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jx
                  else if( laix .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jx
                  else if( laix .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jx
                  else if( laix .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jx
                  end if

                     lni = lni + laix
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( laiy .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jy
                  else if( laiy .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jy
                  else if( laiy .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jy
                  else if( laiy .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jy
                  else if( laiy .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jy
                  end if

                     lni = lni + laiy
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( laiz .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jz
                  else if( laiz .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jz
                  else if( laiz .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jz
                  else if( laiz .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jz
                  else if( laiz .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jz
                  end if

                     lni = lni + laiz
                     dum2(lni:lni) = ')'
                     lni = lni + 1

                  icf = jz + inz * ( jy - 1 )+ iny * inz * ( jx - 1 )

                  do jj = 1, ld

                     nn = li + jj - 1

                     wwt = das(kgwwp(j)+(nn-1)*ixyz+icf-1)

                     write(dum2(lni:lni+14),
     &                         '(1p1e14.5)') wwt

                     lni = lni + 14

                  end do

                     write(iot,'(600a1)') (dum2(i:i),i=1,lni)

               end do
               end do
               end do

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 4 ) then

                     jdsm = inwwc(j)
                     kdsm = kfwwp(j)

                   do m = 1, mnwwp(j,0)

                     lng1=8
                     write(dum1(1:lng1),'(i8)')
     &                    idas_inwwc(jdsm+1+m)

                     dum2(1:3) = '   '

                     do i = 1, lng1
                      dum2(i+3:i+3) = dum1(i:i)
                     end do

                     lngb = lngmax - lng1

                     do i = 1, lngb + 2
                      dum2(i+lng1+3:i+lng1+3) = ' '
                     end do

                     lni = lngb + lng1 + 3 + 2 + 1

                     do jj = 1, ld

                      nn = li + jj - 1

                      wwt = das_kfwwp(kdsm+(nn-1)*mnwwp(j,0)+m-1)

                      write(dum2(lni:lni+14),
     &                     '(1p1e14.5)') wwt

                      lni = lni + 14

                     enddo

                     write(iot,'(600a1)') (dum2(i:i),i=1,lni)

                   enddo

            end if

*-----------------------------------------------------------------------

            end do

         end do

      end if

*-----------------------------------------------------------------------
*     repeated collisions
*-----------------------------------------------------------------------

         if( ircln .gt. 0 .and. ierrg .eq. 0 ) then

            do j = 1, ircln

                  write(iot,'(/"[ Repeated Collisions ]")')

                  nump = mnrcl(j,20)

               if( nump .lt. 19 ) then

                  do k = 1, nump

                     uname(k) = pname(mnrcl(j,k))(1:8)

                  end do

               else

                     nump = 1
                     uname(1) = pname(20)(1:8)

               end if

                     write(iot,'("   part = ",a8,19(1x,a8))')
     &               ( uname(i), i = 1, nump )

*-----------------------------------------------------------------------

               if( irman(j) .eq. 0 ) then

                  write(iot,'("   mother =  all"11x,
     &            " # (D=all) number of specific mother")')

               else if( irman(j) .gt. 0 ) then

                  write(iot,'("   mother = ",i4,11x,
     &            " # (D=all) number of specific mother")')
     &            irman(j) * irmct(j)

                     k = 0

                  do jj = 1, irman(j)

                     k = k + 1

                     iaz = ismat_jsmat( jsmat(j) + jj - 1 )

                     iz = iaz / 1000
                     ia = iaz - iz * 1000

                     if( ia .eq. 0 ) then

                        chaur(k) = elmnts(iz)

                     else

                        call chname(idum,ia,iz,chauu)

                        chaur(k) = chauu

                     end if

                     if( k / 7 * 7 .eq. k .or. jj .eq. irman(j) ) then

                        write(iot,'(12x,7(a8))') ( chaur(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               if( irpem(j,1) .gt. 0 ) then
                     write(iot,'("   emin = ",1p1g15.7)') erpem(j,1)
               end if

               if( irpem(j,2) .gt. 0 ) then
                     write(iot,'("   emax = ",1p1g15.7)') erpem(j,2)
               end if

*-----------------------------------------------------------------------

                     igm = ( mmmax - 1 ) * 2 + 1
                     lngmax = 0

                     idsm = inrlc(j)
                     jdsm = 0

               do m = 1, mnrcl(j,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inrlc(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inrlc(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     call echrg2(mtrn,idas_inrlc(kssm)
     &                           ,dum1,lng1,icmb,igm)

                     if( lng1 .gt. lngmax ) lngmax = lng1

               end do

                     lblk = lngmax
                     dum2(1:6) = '   reg'

                  do i = 1, lblk
                     dum2(i+6:i+6) = ' '
                  end do

                     lblk = lblk + 6

                     dum2(lblk+1:lblk+11) = '     n-coll'

                  iflag_evap = 0
                  if( allocated(idas_lrcls) ) then
                     if( j < 6 ) then
                        if( lrcls(j+1)-lrcls(j) > 0 ) then
                           iflag_evap = 1
                        end if
                     else
                        if( idas_lrcls(lrcls(j)) > 0 ) then
                           iflag_evap = 1
                        end if
                     end if
                  end if
                  if( iflag_evap == 1 ) then

                     lblk = lblk + 11
                     dum2(lblk+1:lblk+11) = '     n-evap'

                  end if

                     write(iot,'(600a1)') (dum2(i:i),i=1,lblk+11)

*-----------------------------------------------------------------------

                     idsm = inrlc(j)
                     jdsm = 0

                     kdsm = krcls(j)
                     ldsm = lrcls(j)

               do m = 1, mnrcl(j,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inrlc(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inrlc(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     ircls = idas_krcls(kdsm-1+m)
                     call echrg2(mtrn,idas_inrlc(kssm)
     &                           ,dum1,lng1,icmb,igm)

                     dum2(1:3) = '   '

                  do i = 1, lng1
                     dum2(i+3:i+3) = dum1(i:i)
                  end do

                     lngb = lngmax - lng1

                  do i = 1, lngb
                     dum2(i+lng1+3:i+lng1+3) = ' '
                  end do

                     lngb = lngb + lng1 + 3

                     write(dum2(lngb+1:lngb+11),'(6x,i5)') ircls

                  if( iflag_evap == 1 ) then
                     iecls = idas_lrcls(ldsm-1+m)

                     lngb = lngb + 11
                     write(dum2(lngb+1:lngb+11),'(6x,i5)') iecls

                  end if

                     write(iot,'(600a1)') (dum2(i:i),i=1,lngb+11)

               end do

            end do

         end if

*-----------------------------------------------------------------------
*     forced collisions
*-----------------------------------------------------------------------

         if( ifcln .gt. 0 .and. ierrg .eq. 0 ) then

            do j = 1, ifcln

                  write(iot,'(/"[ Forced Collisions ]")')

                  nump = mnfcl(j,20)

               if( nump .lt. 19 ) then

                  do k = 1, nump

                     uname(k) = pname(mnfcl(j,k))(1:8)

                  end do

               else

                     nump = 1
                     uname(1) = pname(20)(1:8)

               end if

                     write(iot,'("   part = ",a8,19(1x,a8))')
     &               ( uname(i), i = 1, nump )

*-----------------------------------------------------------------------

                     igm = ( mmmax - 1 ) * 2 + 1
                     lngmax = 0

                     idsm = inflc(j)
                     jdsm = 0

               do m = 1, mnfcl(j,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inflc(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inflc(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     call echrg2(mtrn,idas_inflc(kssm)
     &                           ,dum1,lng1,icmb,igm)

                     if( lng1 .gt. lngmax ) lngmax = lng1

               end do

                     lblk = lngmax
                     dum2(1:6) = '   reg'

                   do i = 1, lblk
                     dum2(i+6:i+6) = ' '
                   end do

                     dum2(lblk+6+1:lblk+6+7) =
     &                          '    fcl'

                     write(iot,'(600a1)') (dum2(i:i),i=1,lblk+6+7)

*-----------------------------------------------------------------------

                     idsm = inflc(j)
                     jdsm = 0

                     kdsm = kfcls(j)

               do m = 1, mnfcl(j,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inflc(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inflc(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     rfcls = das_kfcls(kdsm-1+m)

                     call echrg2(mtrn,idas_inflc(kssm)
     &                           ,dum1,lng1,icmb,igm)

                     dum2(1:3) = '   '

                  do i = 1, lng1
                     dum2(i+3:i+3) = dum1(i:i)
                  end do

                     lngb = lngmax - lng1

                  do i = 1, lngb
                     dum2(i+lng1+3:i+lng1+3) = ' '
                  end do

                     write(dum2(lngb+lng1+4:lngb+lng1+22),
     &                         '(3x,1p1g15.7)') rfcls

                     write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+22)

               end do

            end do

         end if

*-----------------------------------------------------------------------
*     splitting
*-----------------------------------------------------------------------

      if( isptn .gt. 0 .and. ierrg .eq. 0 ) then

      do ll = 1, isptn

         write(iot,'(/"[ Splitting ]")')

*-----------------------------------------------------------------------

                  nump = mnspt(ll,20)

               if( nump .lt. 19 ) then

                  do k = 1, nump

                     uname(k) = pname(mnspt(ll,k))(1:8)

                  end do

               else

                     nump = 1
                     uname(1) = pname(20)(1:8)

               end if

                     write(iot,'("     part = ",a8,19(1x,a8))')
     &               ( uname(i), i = 1, nump )

*-----------------------------------------------------------------------

               if( ispem(ll,1) .ne. 0 ) then

                  write(iot,'("     emin =",1p1g15.7)')
     &                         espem(ll,1)

               end if

               if( ispem(ll,2) .ne. 0 ) then

                  write(iot,'("     emax =",1p1g15.7)')
     &                         espem(ll,2)

               end if

*-----------------------------------------------------------------------

            do m = 1, 3

               if( ispct(ll,m) .ne. 0 ) then

                  write(iot,'(" ctmin(",i1,") = ",i5)')
     &                         m, ispct(ll,m*2+2)
                  write(iot,'(" ctmax(",i1,") = ",i5)')
     &                         m, ispct(ll,m*2+3)

               end if

            end do

*-----------------------------------------------------------------------

                  lgmax1 = 0
                  lgmax2 = 0

                  igm = ( mmmax - 1 ) * 2 + 2
                  idsm = ipgrc(ll)
                  jdsm = -1

            do m = 1, npreg(ll)

                  jdsm = jdsm + 1
                  ntrn = idas_ipgrc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_ipgrc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  call echrg2(mtrn,idas_ipgrc(kssm),dum1,lng1,icmb,igm)

                  if( lng1 .gt. lgmax1 ) lgmax1 = lng1

                  jdsm = jdsm + 1
                  ntrn = idas_ipgrc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_ipgrc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  call echrg2(mtrn,idas_ipgrc(kssm),dum1,lng1,icmb,igm)

                  if( lng1 .gt. lgmax2 ) lgmax2 = lng1

            end do

                  dum2(1:5) = '  non'
                  l = 5
                  dum2(l+1:l+9) = '     r-in'
                  l = l + 9
                  lblk = lgmax1
                do i = 1, lblk
                  dum2(i+l:i+l) = ' '
                end do
                  l = l + lblk
                  dum2(l+1:l+6) =
     &                       ' r-out'
                  l = l + 6
                  lblk = lgmax2
                do i = 1, lblk
                  dum2(i+l:i+l) = ' '
                end do
                  l = l + lblk
                  dum2(l+1:l+9) =
     &                       '   factor'
                  l = l + 9

                  write(iot,'(600a1)') (dum2(i:i),i=1,l)

*-----------------------------------------------------------------------

                  idsm = ipgrc(ll)
                  jdsm = -1

                  kdsm = ksplt(ll)
                  ldsm = 0

            do j = 1, npreg(ll)

                  jdsm = jdsm + 1
                  ntrn = idas_ipgrc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_ipgrc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  ldsm  = ldsm + 1
                  sfc   = das_ksplt(kdsm+ldsm)

                  call echrg2(mtrn,idas_ipgrc(kssm),dum1,lng1,icmb,igm)

                     l = 2
                     dum2(1:l) = '  '
                     write(dum2(l+1:l+3),'(i3)') j
                     l = l + 3
                     dum2(l+1:l+5) = '     '
                     l = l + 5
                  do i = 1, lng1
                     dum2(i+l:i+l) = dum1(i:i)
                  end do
                     l = l + lng1
                     lngb = lgmax1 - lng1
                  do i = 1, lngb
                     dum2(i+l:i+l) = ' '
                  end do
                     l = l + lngb

                  jdsm = jdsm + 1
                  ntrn = idas_ipgrc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_ipgrc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  call echrg2(mtrn,idas_ipgrc(kssm),dum1,lng1,icmb,igm)

                     dum2(l+1:l+5) = '     '
                     l = l + 5
                  do i = 1, lng1
                     dum2(i+l:i+l) = dum1(i:i)
                  end do
                     l = l + lng1
                     lngb = lgmax2 - lng1
                  do i = 1, lngb
                     dum2(i+l:i+l) = ' '
                  end do
                     l = l + lngb

                  write(dum2(l+1:l+20),'(5x,1p1g15.7)') sfc
                     l = l + 20

                     write(iot,'(600a1)') (dum2(i:i),i=1,l)

            end do

      end do

      end if

*-----------------------------------------------------------------------
*     super mirror
*-----------------------------------------------------------------------

      if( nsreg .gt. 0 .and. ierrg .eq. 0 ) then

         write(iot,'(/"[ Super Mirror ]")')

                  lgmax1 = 0
                  lgmax2 = 0

                  igm = ( mmmax - 1 ) * 2 + 2
                  idsm = isgrc
                  jdsm = -1

            do m = 1, nsreg

                  jdsm = jdsm + 1
                  ntrn = idas_isgrc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_isgrc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  call echrg2(mtrn,idas_isgrc(kssm),dum1,lng1,icmb,igm)

                  if( lng1 .gt. lgmax1 ) lgmax1 = lng1

                  jdsm = jdsm + 1
                  ntrn = idas_isgrc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_isgrc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  call echrg2(mtrn,idas_isgrc(kssm),dum1,lng1,icmb,igm)

                  if( lng1 .gt. lgmax2 ) lgmax2 = lng1

            end do

                  dum2(1:5) = '  non'
                  l = 5
                  dum2(l+1:l+9) = '     r-in'
                  l = l + 9
                  lblk = lgmax1
                do i = 1, lblk
                  dum2(i+l:i+l) = ' '
                end do
                  l = l + lblk
                  dum2(l+1:l+6) =
     &                       ' r-out'
                  l = l + 6
                  lblk = lgmax2
                do i = 1, lblk
                  dum2(i+l:i+l) = ' '
                end do
                  l = l + lblk
                  dum2(l+1:l+10) =
     &                       ' Mm       '
                  l = l + 10
                  dum2(l+1:l+10) =
     &                       ' R0       '
                  l = l + 10
                  dum2(l+1:l+10) =
     &                       ' Qc       '
                  l = l + 10
                  dum2(l+1:l+10) =
     &                       ' Am       '
                  l = l + 10
                  dum2(l+1:l+10) =
     &                       ' Wm       '
                  l = l + 10

                  write(iot,'(600a1)') (dum2(i:i),i=1,l)

*-----------------------------------------------------------------------

                  idsm = isgrc
                  jdsm = -1

                  kdsm = ksmir
                  ldsm = 0

            do j = 1, nsreg

                  jdsm = jdsm + 1
                  ntrn = idas_isgrc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_isgrc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  ldsm  = ldsm + 1
                  smm   = das_ksmir(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  srm   = das_ksmir(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  sqc   = das_ksmir(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  sam   = das_ksmir(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  swm   = das_ksmir(kdsm+ldsm)

                  call echrg2(mtrn,idas_isgrc(kssm),dum1,lng1,icmb,igm)

                     l = 2
                     dum2(1:l) = '  '
                     write(dum2(l+1:l+3),'(i3)') j
                     l = l + 3
                     dum2(l+1:l+5) = '     '
                     l = l + 5
                  do i = 1, lng1
                     dum2(i+l:i+l) = dum1(i:i)
                  end do
                     l = l + lng1
                     lngb = lgmax1 - lng1
                  do i = 1, lngb
                     dum2(i+l:i+l) = ' '
                  end do
                     l = l + lngb

                  jdsm = jdsm + 1
                  ntrn = idas_isgrc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_isgrc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  call echrg2(mtrn,idas_isgrc(kssm),dum1,lng1,icmb,igm)

                     dum2(l+1:l+5) = '     '
                     l = l + 5
                  do i = 1, lng1
                     dum2(i+l:i+l) = dum1(i:i)
                  end do
                     l = l + lng1
                     lngb = lgmax2 - lng1
                  do i = 1, lngb
                     dum2(i+l:i+l) = ' '
                  end do
                     l = l + lngb

                  write(dum2(l+1:l+53),'(3x,5f10.5)')
     &               smm, srm, sqc, sam, swm
                     l = l + 53

                     write(iot,'(600a1)') (dum2(i:i),i=1,l)

            end do

      end if

*-----------------------------------------------------------------------
*     magnetic field
*-----------------------------------------------------------------------

      if( nmreg .gt. 0 .and. mstz(14) .ne. 0 .and. ierrg .eq. 0 ) then

         write(iot,'(/"[ Magnetic  Field ]")')

                  lngmax = 0

                  igm = ( mmmax - 1 ) * 2 + 1
                  idsm = ingrc
                  jdsm = 0

            do m = 1, nmreg

                  jdsm = jdsm + 1
                  ntrn = idas_ingrc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_ingrc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  call echrg2(mtrn,idas_ingrc(kssm),dum1,lng1,icmb,igm)

                  if( lng1 .gt. lngmax ) lngmax = lng1

            end do

                  lblk = lngmax
                  dum2(1:6) = '   reg'

                do i = 1, lblk
                  dum2(i+6:i+6) = ' '
                end do

*-----------------------------------------------------------------------

                  kdsm = kmags
                  ldsm = 0
                  lphs = 0
                  ltrc = 0
                  ltim = 0

               do m = 1, nmreg

                  ldsm   = ldsm + 4
                  p_mag  = das_kmags(kdsm+ldsm)
                  ldsm   = ldsm + 2
                  t_mag  = das_kmags(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  u_mag  = das_kmags(kdsm+ldsm)

                  if( p_mag .gt. -10000.0d0 ) lphs = lphs + 1
                  if( nint( t_mag ) .ne. 0 )  ltrc = ltrc + 1
                  if( u_mag .gt. -1.0d+9 )    ltim = ltim + 1

               end do

*-----------------------------------------------------------------------

            if( ltim .eq. 0 ) then

               if( ltrc .eq. 0 .and. lphs .eq. 0 ) then

                  dum2(lblk+6+1:lblk+6+26) =
     &            ' typ    gap            mgf'

                  write(iot,'(600a1)') (dum2(i:i),i=1,lblk+6+26)

               else if( ltrc .eq. 0 .and. lphs .gt. 0 ) then

                  dum2(lblk+6+1:lblk+6+41) =
     &            ' typ    gap            mgf          polar'

                  write(iot,'(600a1)') (dum2(i:i),i=1,lblk+6+41)

               else if( ltrc .gt. 0 .and. lphs .eq. 0 ) then

                  dum2(lblk+6+1:lblk+6+41) =
     &            ' typ    gap            mgf           trcl'

                  write(iot,'(600a1)') (dum2(i:i),i=1,lblk+6+41)

               else if( ltrc .gt. 0 .and. lphs .gt. 0 ) then

                  dum2(lblk+6+1:lblk+6+51) =
     &       ' typ    gap            mgf           trcl     polar'

                  write(iot,'(600a1)') (dum2(i:i),i=1,lblk+6+51)

               end if

*-----------------------------------------------------------------------

            else if( ltim .gt. 0 ) then

               if( ltrc .eq. 0 .and. lphs .eq. 0 ) then

                  dum2(lblk+6+1:lblk+6+41) =
     &            ' typ    gap            mgf'//
     &            '           time'

                  write(iot,'(600a1)') (dum2(i:i),i=1,lblk+6+41)

               else if( ltrc .eq. 0 .and. lphs .gt. 0 ) then

                  dum2(lblk+6+1:lblk+6+57) =
     &            ' typ    gap            mgf            polar'//
     &            '          time'

                  write(iot,'(600a1)') (dum2(i:i),i=1,lblk+6+57)

               else if( ltrc .gt. 0 .and. lphs .eq. 0 ) then

                  dum2(lblk+6+1:lblk+6+50) =
     &            ' typ    gap            mgf           trcl'//
     &            '     time'

                  write(iot,'(600a1)') (dum2(i:i),i=1,lblk+6+50)

               else if( ltrc .gt. 0 .and. lphs .gt. 0 ) then

                  dum2(lblk+6+1:lblk+6+65) =
     &     ' typ    gap            mgf           trcl     polar'//
     &            '          time'

                  write(iot,'(600a1)') (dum2(i:i),i=1,lblk+6+65)

               end if

            end if

*-----------------------------------------------------------------------

                  idsm = ingrc
                  jdsm = 0

                  kdsm = kmags
                  ldsm = 0

            do m = 1, nmreg

                  jdsm = jdsm + 1
                  ntrn = idas_ingrc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_ingrc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  ldsm  = ldsm + 1
                  a_mag = das_kmags(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  b_mag = das_kmags(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  s_mag = das_kmags(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  p_mag = das_kmags(kdsm+ldsm)
                  ldsm  = ldsm + 2
                  t_mag = das_kmags(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  u_mag = das_kmags(kdsm+ldsm)

                  call echrg2(mtrn,idas_ingrc(kssm),dum1,lng1,icmb,igm)

                  dum2(1:3) = '   '

               do i = 1, lng1
                  dum2(i+3:i+3) = dum1(i:i)
               end do

                  lngb = lngmax - lng1

               do i = 1, lngb
                  dum2(i+lng1+3:i+lng1+3) = ' '
               end do

                     chbd = 'non'

*-----------------------------------------------------------------------

            if( ltim .eq. 0 ) then

               if( ltrc .eq. 0 .and. lphs .eq. 0 ) then

                  write(dum2(lngb+lng1+4:lngb+lng1+41),
     &                      '(3x,i3,2x,1p2g15.7)')
     &               nint(s_mag), a_mag, b_mag

                  write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+41)

               else if( ltrc .eq. 0 .and. lphs .gt. 0 ) then

                  if( p_mag .gt. -10000.0d0 ) then

                     write(dum2(lngb+lng1+4:lngb+lng1+56),
     &                         '(3x,i3,2x,1p2g15.7,1p1g15.7)')
     &                  nint(s_mag), a_mag, b_mag, p_mag

                     write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+56)

                  else

                     write(dum2(lngb+lng1+4:lngb+lng1+47),
     &                         '(3x,i3,2x,1p2g15.7,3x,a3)')
     &                  nint(s_mag), a_mag, b_mag, chbd

                     write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+47)

                  end if

               else if( ltrc .gt. 0 .and. lphs .eq. 0 ) then

                  write(dum2(lngb+lng1+4:lngb+lng1+47),
     &                      '(3x,i3,2x,1p2g15.7,2x,i4)')
     &               nint(s_mag), a_mag, b_mag, nint(t_mag)

                  write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+47)

               else if( ltrc .gt. 0 .and. lphs .gt. 0 ) then

                  if( p_mag .gt. -10000.0d0 ) then

                     write(dum2(lngb+lng1+4:lngb+lng1+64),
     &                         '(3x,i3,2x,1p2g15.7,2x,i4,2x,1p1g15.7)')
     &                  nint(s_mag), a_mag, b_mag, nint(t_mag), p_mag

                     write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+64)

                  else

                     write(dum2(lngb+lng1+4:lngb+lng1+55),
     &                         '(3x,i3,2x,1p2g15.7,2x,i4,5x,a3)')
     &                  nint(s_mag), a_mag, b_mag, nint(t_mag), chbd

                     write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+55)

                  end if

               end if

*-----------------------------------------------------------------------

            else if( ltim .gt. 0 .and. u_mag .gt. -1.0d+9 ) then

               if( ltrc .eq. 0 .and. lphs .eq. 0 ) then

                  write(dum2(lngb+lng1+4:lngb+lng1+56),
     &                      '(3x,i3,2x,1p3g15.7)')
     &               nint(s_mag), a_mag, b_mag, u_mag

                  write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+56)

               else if( ltrc .eq. 0 .and. lphs .gt. 0 ) then

                  if( p_mag .gt. -10000.0d0 ) then

                     write(dum2(lngb+lng1+4:lngb+lng1+71),
     &                         '(3x,i3,2x,1p2g15.7,1p2g15.7)')
     &                  nint(s_mag), a_mag, b_mag, p_mag, u_mag

                     write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+71)

                  else

                     write(dum2(lngb+lng1+4:lngb+lng1+71),
     &                         '(3x,i3,2x,1p2g15.7,3x,a3,9x,1p1g15.7)')
     &                  nint(s_mag), a_mag, b_mag, chbd, u_mag

                     write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+71)

                  end if

               else if( ltrc .gt. 0 .and. lphs .eq. 0 ) then

                  write(dum2(lngb+lng1+4:lngb+lng1+64),
     &                      '(3x,i3,2x,1p2g15.7,2x,i4,2x,1p1g15.7)')
     &               nint(s_mag), a_mag, b_mag, nint(t_mag), u_mag

                  write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+64)

               else if( ltrc .gt. 0 .and. lphs .gt. 0 ) then

                  if( p_mag .gt. -10000.0d0 ) then

                     write(dum2(lngb+lng1+4:lngb+lng1+79),
     &                         '(3x,i3,2x,1p2g15.7,2x,i4,2x,1p2g15.7)')
     &                  nint(s_mag), a_mag, b_mag, nint(t_mag),
     &                  p_mag, u_mag

                     write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+79)

                  else

                     write(dum2(lngb+lng1+4:lngb+lng1+79),
     &                         '(3x,i3,2x,1p2g15.7,2x,i4,5x,a3,
     &                           9x,1p1g15.7)')
     &                  nint(s_mag), a_mag, b_mag, nint(t_mag),
     &                  chbd, u_mag

                     write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+79)

                  end if

               end if

*-----------------------------------------------------------------------

            else if( ltim .gt. 0 .and. u_mag .le. -1.0d+9 ) then

               if( ltrc .eq. 0 .and. lphs .eq. 0 ) then

                  write(dum2(lngb+lng1+4:lngb+lng1+47),
     &                      '(3x,i3,2x,1p2g15.7,3x,3a)')
     &               nint(s_mag), a_mag, b_mag, chbd

                  write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+47)

               else if( ltrc .eq. 0 .and. lphs .gt. 0 ) then

                  if( p_mag .gt. -10000.0d0 ) then

                     write(dum2(lngb+lng1+4:lngb+lng1+62),
     &                         '(3x,i3,2x,1p3g15.7,3x,3a)')
     &                  nint(s_mag), a_mag, b_mag, p_mag, chbd

                     write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+62)

                  else

                     write(dum2(lngb+lng1+4:lngb+lng1+62),
     &                         '(3x,i3,2x,1p2g15.7,3x,a3,9x,3x,3a)')
     &                  nint(s_mag), a_mag, b_mag, chbd, chbd

                     write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+62)

                  end if

               else if( ltrc .gt. 0 .and. lphs .eq. 0 ) then

                  write(dum2(lngb+lng1+4:lngb+lng1+55),
     &                      '(3x,i3,2x,1p2g15.7,2x,i4,2x,3x,3a)')
     &               nint(s_mag), a_mag, b_mag, nint(t_mag), chbd

                  write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+55)

               else if( ltrc .gt. 0 .and. lphs .gt. 0 ) then

                  if( p_mag .gt. -10000.0d0 ) then

                     write(dum2(lngb+lng1+4:lngb+lng1+70),
     &                         '(3x,i3,2x,1p2g15.7,2x,i4,2x,1p1g15.7,
     &                           3x,3a)')
     &                  nint(s_mag), a_mag, b_mag, nint(t_mag),
     &                  p_mag, chbd

                     write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+70)

                  else

                     write(dum2(lngb+lng1+4:lngb+lng1+70),
     &                         '(3x,i3,2x,1p2g15.7,2x,i4,5x,a3,
     &                           9x,3x,3a)')
     &                  nint(s_mag), a_mag, b_mag, nint(t_mag),
     &                  chbd, chbd

                     write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+70)

                  end if

               end if

            end if

*-----------------------------------------------------------------------

            end do

      end if

*-----------------------------------------------------------------------
*     electro magnetic field
*-----------------------------------------------------------------------

      if( nereg .gt. 0 .and. mstz(70) .ne. 0 .and. ierrg .eq. 0 ) then

         write(iot,'(/"[ Electro Magnetic Field ]")')

                  lngmax = 0

                  igm = ( mmmax - 1 ) * 2 + 1
                  idsm = inerc
                  jdsm = 0

            do m = 1, nereg

                  jdsm = jdsm + 1
                  ntrn = idas_inerc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_inerc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  call echrg2(mtrn,idas_inerc(kssm),dum1,lng1,icmb,igm)

                  if( lng1 .gt. lngmax ) lngmax = lng1

            end do

                  lblk = lngmax
                  dum2(1:6) = '   reg'

                do i = 1, lblk
                  dum2(i+6:i+6) = ' '
                end do

*-----------------------------------------------------------------------

                  dum2(lblk+6+1:lblk+6+45) =
     &            '  elf            mgf            trcle'//
     &            '   trclm'

                  write(iot,'(600a1)') (dum2(i:i),i=1,lblk+6+45)

*-----------------------------------------------------------------------

                  idsm = inerc
                  jdsm = 0

                  kdsm = kelcs
                  ldsm = 0

            do m = 1, nereg

                  jdsm = jdsm + 1
                  ntrn = idas_inerc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_inerc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  ldsm  = ldsm + 1
                  s_elf = das_kelcs(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  s_mgf = das_kelcs(kdsm+ldsm)
                  ldsm  = ldsm + 2
                  t_elf = das_kelcs(kdsm+ldsm)
                  ldsm  = ldsm + 2
                  t_mgf = das_kelcs(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  emap_type = das_kelcs(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  mmap_type = das_kelcs(kdsm+ldsm)
                  ldsm  = ldsm + 1
                  a_elmg = das_kelcs(kdsm+ldsm) ! gap

                  call echrg2(mtrn,idas_inerc(kssm),dum1,lng1,icmb,igm)

                  dum2(1:3) = '   '

               do i = 1, lng1
                  dum2(i+3:i+3) = dum1(i:i)
               end do

                  lngb = lngmax - lng1

               do i = 1, lngb
                  dum2(i+lng1+3:i+lng1+3) = ' '
               end do

*-----------------------------------------------------------------------

                  write(dum2(lngb+lng1+4:lngb+lng1+51),
     &                      '(2x,1p2g15.7,2i8)')
     &               s_elf, s_mgf, nint(t_elf), nint(t_mgf)

                  write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+51)

*-----------------------------------------------------------------------

            end do

      end if

*-----------------------------------------------------------------------
*     timer
*-----------------------------------------------------------------------

      if( ierrg .eq. 0 .and. ntmrg .ne. 0 ) then

            write(iot,'(/"[ Timer ]")')

                  lngmax = 0

                  igm = ( mmmax - 1 ) * 2 + 1
                  idsm = intmc
                  jdsm = 0

            do m = 1, ntmrg

                  jdsm = jdsm + 1
                  ntrn = idas_intmc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_intmc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  call echrg2(mtrn,idas_intmc(kssm),dum1,lng1,icmb,igm)

                  if( lng1 .gt. lngmax ) lngmax = lng1

            end do

                  lblk = lngmax
                  dum2(1:6) = '   reg'

                do i = 1, lblk
                  dum2(i+6:i+6) = ' '
                end do

                  dum2(lblk+6+1:lblk+6+29) =
     &                       '   in     out    coll     ref'

                  write(iot,'(600a1)') (dum2(i:i),i=1,lblk+6+29)

*-----------------------------------------------------------------------

                  idsm = intmc
                  jdsm = 0

                  kdsm = ktime
                  ldsm = 0

            do m = 1, ntmrg

                  jdsm = jdsm + 1
                  ntrn = idas_intmc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_intmc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  ldsm   = ldsm + 1
                  inin   = idas_ktime(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  inout  = idas_ktime(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  incol  = idas_ktime(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  inref  = idas_ktime(kdsm+ldsm)

                  call echrg2(mtrn,idas_intmc(kssm),dum1,lng1,icmb,igm)

                  dum2(1:5) = '     '

               do i = 1, lng1
                  dum2(i+5:i+5) = dum1(i:i)
               end do

                  lngb = lngmax - lng1

               do i = 1, lngb
                  dum2(i+lng1+5:i+lng1+5) = ' '
               end do

                  write(dum2(lngb+lng1+6:lngb+lng1+6+32),
     &                      '(4i8)')
     &               inin, inout, incol, inref

                  write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+6+32)

            end do

      end if

*-----------------------------------------------------------------------
*     counter
*-----------------------------------------------------------------------

      if( ierrg .eq. 0 .and.
     &  ( ncntc(1) .eq. 1 .or. ncntc(2) .eq. 1 .or.
     &    ncntc(3) .eq. 1 ) ) then

               icont = 0

            do k = 1, itnm

               do j = 1, 3

                  if( itcnt(j,k) .ne. 0 ) icont = icont + 1

               end do

                  idump = itmdp(k,0)

               if( idump .ne. 0 ) then

                  do j = 1, abs( idump )

                     if( itmdp(k,j) .ge. 11 .and.
     &                   itmdp(k,j) .le. 13 ) icont = icont + 1

                  end do

               end if

               if( ital(k) .eq. 20 ) then
                  icont = icont + 1
               endif

            end do

            if( idumpall .ne. 0 ) icont = icont + 1

         iconth = 0
         do k = 1, itnm
           do j = 1, 3
             if( itcnth(j,k) .ne. 0 ) iconth = iconth + 1
           end do
         end do
         if ( iconth.gt.0 .and. istdev.ne.2 ) then
           write(jo,*)
     &     'warning: history counter must be used with istdev = 2'
           iconth = 0
           do j = 1, 3
             itcnth(j,:) = 0
           end do
         end if

! T.Sato 2024/03/19 counter bias for [t-wwg]
         iconwei=0
         do k = 1, itnm
          if(ichnum(k).ne.0.or.ictnum(k).ne.0) iconwei=iconwei+1
         enddo

         if( icont.gt.0.or.iconth.gt.0.or.iconwei.gt.0 ) then ! frtati 2021/03/05

            write(iot,'(/"[ Counter ]")')

         else

            write(iot,'(/"[ Counter ] off")')

         end if

*-----------------------------------------------------------------------

         do k = 1, 3
         if( ncntc(k) .eq. 1 ) then

            write(iot,'("  counter =",i2)') k

*-----------------------------------------------------------------------

               call echcprt(1,iot,k,icpan,icpat)

*-----------------------------------------------------------------------

                  lngmax = 0

                  igm = ( mmmax - 1 ) * 2 + 1
                  idsm = incrc(k)
                  jdsm = 0

            do m = 1, ncreg(k)

                  jdsm = jdsm + 1
                  ntrn = idas_incrc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_incrc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  call echrg2(mtrn,idas_incrc(kssm),dum1,lng1,icmb,igm)

                  if( lng1 .gt. lngmax ) lngmax = lng1

            end do

                  lblk = lngmax
                  dum2(1:6) = '   reg'

                do i = 1, lblk
                  dum2(i+6:i+6) = ' '
                end do

                icn = 0
                do i = 1, icnech(k,0)

                  if( icnech(k,i) .eq. 1 ) then
                    cycle
                  elseif( icnech(k,i) .eq. 0 ) then ! S.H. added non case (2022.3.24)
                    cycle
                  elseif( icnech(k,i) .eq. 2 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '      in'
                  elseif( icnech(k,i) .eq. 3 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '     out'
                  elseif( icnech(k,i) .eq. 4 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '    coll'
                  elseif( icnech(k,i) .eq. 5 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '     ref'
                  elseif( icnech(k,i) .eq. 6 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '    fiss'
                  elseif( icnech(k,i) .eq. 7 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '    elst'
                  elseif( icnech(k,i) .eq. 8 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '    iels'
                  elseif( icnech(k,i) .eq. 9 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '    nucl'
                  elseif( icnech(k,i) .eq. 10 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '    dcay'
                  elseif( icnech(k,i) .eq. 11 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '    atom'
                  elseif( icnech(k,i) .eq. 12 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '    delt'
                  elseif( icnech(k,i) .eq. 13 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '    fluo'
                  elseif( icnech(k,i) .eq. 14 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '    auge'
                  elseif( icnech(k,i) .eq. 15 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '    brem'
                  elseif( icnech(k,i) .eq. 16 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '    phel'
                  elseif( icnech(k,i) .eq. 17 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '    cmpt'
                  elseif( icnech(k,i) .eq. 18 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '    pprd'
                  elseif( icnech(k,i) .eq. 19 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '    anih'
                  elseif( icnech(k,i) .eq. 20 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '    msct'
                  elseif( icnech(k,i) .eq. 21 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '    raly'
                  elseif( icnech(k,i) .eq. 22 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '    knoe'
                  elseif( icnech(k,i) .eq. 23 ) then
                    dum2(lblk+6+icn:lblk+6+icn+8) = '   ndata'
                  endif

                  if( icnech(k,i) .ne. 1 ) icn = icn + 8   ! except 'reg'

                enddo

                  write(iot,'(600a1)') (dum2(i:i),i=1,lblk+5+icn)
*-----------------------------------------------------------------------

                  idsm = incrc(k)
                  jdsm = 0

                  kdsm = kcont(k)
                  ldsm = 0

            do m = 1, ncreg(k)

                  jdsm = jdsm + 1
                  ntrn = idas_incrc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_incrc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  ldsm   = ldsm + 1
                  inin   = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  inout  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  incol  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  inref  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  infis  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  inels  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  iniel  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  inncr  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  indcy  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  inato  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  indlr  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  influ  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  inaug  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  inbrm  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  inphe  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  incmp  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  inppd  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  inanh  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  inmst  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  inray  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  inkoe  = idas_kcont(kdsm+ldsm)
                  ldsm   = ldsm + 1
                  indat  = idas_kcont(kdsm+ldsm)

                  call echrg2(mtrn,idas_incrc(kssm),dum1,lng1,icmb,igm)

                  dum2(1:5) = '     '

               do i = 1, lng1
                  dum2(i+5:i+5) = dum1(i:i)
               end do

                  lngb = lngmax - lng1

               do i = 1, lngb
                  dum2(i+lng1+5:i+lng1+5) = ' '
               end do

                icn = 0
                do i = 1, icnech(k,0)

                  if( icnech(k,i) .eq. 1 ) then
                    cycle
                  elseif( icnech(k,i) .eq. 0 ) then ! S.H. added non case (2022.3.24)
                    cycle
                  elseif( icnech(k,i) .eq. 2 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') inin
                  elseif( icnech(k,i) .eq. 3 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') inout
                  elseif( icnech(k,i) .eq. 4 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') incol
                  elseif( icnech(k,i) .eq. 5 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') inref
                  elseif( icnech(k,i) .eq. 6 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') infis
                  elseif( icnech(k,i) .eq. 7 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') inels
                  elseif( icnech(k,i) .eq. 8 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') iniel
                  elseif( icnech(k,i) .eq. 9 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') inncr
                  elseif( icnech(k,i) .eq. 10 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') indcy
                  elseif( icnech(k,i) .eq. 11 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') inato
                  elseif( icnech(k,i) .eq. 12 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') indlr
                  elseif( icnech(k,i) .eq. 13 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') influ
                  elseif( icnech(k,i) .eq. 14 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') inaug
                  elseif( icnech(k,i) .eq. 15 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') inbrm
                  elseif( icnech(k,i) .eq. 16 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') inphe
                  elseif( icnech(k,i) .eq. 17 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') incmp
                  elseif( icnech(k,i) .eq. 18 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') inppd
                  elseif( icnech(k,i) .eq. 19 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') inanh
                  elseif( icnech(k,i) .eq. 20 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') inmst
                  elseif( icnech(k,i) .eq. 21 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') inray
                  elseif( icnech(k,i) .eq. 22 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') inkoe
                  elseif( icnech(k,i) .eq. 23 ) then
                    write(dum2(lngb+lng1+6+icn:lngb+lng1+6+icn+8),
     &              '(i8)') indat
                  endif

                  if( icnech(k,i) .ne. 1 ) icn = icn + 8   ! except 'reg'

                enddo

                  write(iot,'(600a1)') (dum2(i:i),i=1,lngb+lng1+6+icn)

            end do

         end if
         end do

*-----------------------------------------------------------------------
*        if there is no tally with counter, we neglect counter process
*-----------------------------------------------------------------------

         if( icont.eq.0.and.iconth.eq.0.and.iconwei.eq.0 ) then ! frtati 2021/03/05, T.Sato 2024/03/19

            do j = 1, 3

               ncntc(j) = 0

            end do

               write(iot,'("# Warning: there is no tally with ",
     &                     "counter.")')

         end if

      end if

*-----------------------------------------------------------------------
*     Mat Name Color
*-----------------------------------------------------------------------

      if( ierrg .eq. 0 .and. mtncn .gt. 0 ) then

            write(iot,'(/"[ Mat Name Color ]")')

                  inmax = 0
                  icmax = 0

            do i = 1, mxmat

                  inp = 0
                  ins = 0
               do k = 1, nmtnm(i)
                  if( dmtnm(i)(k:k) .eq. ' ' ) ins = 2
                  if( dmtnm(i)(k:k) .eq. '{' ) inp = inp + 1
                  if( dmtnm(i)(k:k) .eq. '}' ) inp = inp + 1
               end do

                  ics = 0
               do k = 1, nmtcl(i)
                  if( dmtcl(i)(k:k) .eq. ' ' ) ics = 2
               end do

                  nmtmg = nmtnm(i) + ins + inp
                  nmtcg = nmtcl(i) + ics

               if( nmtmg .gt. inmax ) inmax = nmtmg
               if( nmtcg .gt. icmax ) icmax = nmtcg

            end do

                  nctemp(1:6) = '   mat'
                  ic = 6

               do k = 1, inmax
                  nctemp(ic+k:ic+k) = ' '
               end do

                  ic = ic + inmax
                  nctemp(ic+1:ic+4) = 'name'
                  ic = ic + 4

                  nctemp(ic+1:ic+9) = '     size'
                  ic = ic + 9

               do k = 1, icmax
                  nctemp(ic+k:ic+k) = ' '
               end do

                  ic = ic + icmax
                  nctemp(ic+1:ic+5) = 'color'
                  ic = ic + 5

                  nctemp(ic+1:ic+20) = '   # HSB color value'
                  ic = ic + 20

               write(iot,'(200a1)') (nctemp(k:k),k=1,ic)

*-----------------------------------------------------------------------

            do i = 0, mxmat

                  inp = 0
                  ins = 0
               do k = 1, nmtnm(i)
                  if( dmtnm(i)(k:k) .eq. ' ' ) ins = 2
                  if( dmtnm(i)(k:k) .eq. '{' ) inp = inp + 1
                  if( dmtnm(i)(k:k) .eq. '}' ) inp = inp + 1
               end do

                  ics = 0
               do k = 1, nmtcl(i)
                  if( dmtcl(i)(k:k) .eq. ' ' ) ics = 2
               end do

                  nmtmg = nmtnm(i) + ins + inp
                  nmtcg = nmtcl(i) + ics

               write(nctemp(1:6),'(i6)') idmn(i)

                  ic = 6
                  icc = 4 + inmax - nmtmg

               do k = 1, icc
                  nctemp(ic+k:ic+k) = ' '
               end do

                  ic = ic + icc

               if( ins .gt. 0 ) then

                  nctemp(ic+1:ic+1) = '{'
                  ic = ic + 1

               end if

                        ii = 0
                  do jj = 1, nmtnm(i)
                     if( dmtnm(i)(jj:jj) .eq. '{' .or.
     &                   dmtnm(i)(jj:jj) .eq. '}' ) then
                        ii = ii + 1
                        nctemp(ic+ii:ic+ii) = yen
                        ii = ii + 1
                        nctemp(ic+ii:ic+ii) = dmtnm(i)(jj:jj)
                     else
                        ii = ii + 1
                        nctemp(ic+ii:ic+ii) = dmtnm(i)(jj:jj)
                     end if
                  end do

                  ic = ic + ii

               if( ins .gt. 0 ) then

                  nctemp(ic+1:ic+1) = '}'
                  ic = ic + 1

               end if

               write(nctemp(ic+1:ic+9),'(f9.2)') dmhsb(i,4)

                  ic = ic + 9

*-----------------------------------------------------------------------

                  icc = 5 + icmax - nmtcg

               do k = 1, icc
                  nctemp(ic+k:ic+k) = ' '
               end do

                  ic = ic + icc

               if( ics .gt. 0 ) then

                  nctemp(ic+1:ic+1) = '{'
                  ic = ic + 1

               end if

                  nctemp(ic+1:ic+nmtcl(i)) = dmtcl(i)(1:nmtcl(i))

                  ic = ic + nmtcl(i)

               if( ics .gt. 0 ) then

                  nctemp(ic+1:ic+1) = '}'
                  ic = ic + 1

               end if

*-----------------------------------------------------------------------

                  nctemp(ic+1:ic+6) = '   # {'
                  ic = ic + 6

                  write(nctemp(ic+1:ic+21),'(3f7.3)')
     &            dmhsb(i,1), dmhsb(i,2), dmhsb(i,3)

                  ic = ic + 21

                  nctemp(ic+1:ic+2) = ' }'
                  ic = ic + 2

*-----------------------------------------------------------------------

               write(iot,'(200a1)') (nctemp(k:k),k=1,ic)

            end do

      end if

*-----------------------------------------------------------------------
*     Mat Time Change
*-----------------------------------------------------------------------

      if( ierrg .eq. 0 .and. mttcn .gt. 0 ) then

            write(iot,'(/"[ Mat Time Change ]")')

                  nctemp(1:6) = '   mat'
                  ic = 6
                  nctemp(ic+1:ic+15) = '       time    '
                  ic = ic + 15
                  nctemp(ic+1:ic+6) = 'change'
                  ic = ic + 6

               write(iot,'(200a1)') (nctemp(k:k),k=1,ic)

            do i = 1, mttcn

               write(iot,'(i6,1p1g15.7,i6)')
     &                mttc1(i), smttc(i), mttc2(i)

            end do

      end if

*-----------------------------------------------------------------------
*     Elastic option
*-----------------------------------------------------------------------

      if( ierrg .eq. 0 .and. mlrgn .gt. 0 ) then

            write(iot,'(/"[ Elastic Option ]")')

                  irmax = 0

            do i = 1, mlrgn

                  write(nrtem1,'(i7)') melrg(i,1)
                  do k = 1, 7
                     if( nrtem1(k:k) .ne. ' ' ) goto 510
                  end do
  510             isl1 = 8 - k

               if( melrg(i,1) .eq. melrg(i,2) ) then

                  isll = isl1

               else

                  write(nrtem2,'(i7)') melrg(i,2)
                  do k = 1, 7
                     if( nrtem2(k:k) .ne. ' ' ) goto 520
                  end do
  520             isl2 = 8 - k

                  isll = isl1 + isl2 + 3

               end if

                  if( isll .gt. irmax ) irmax = isll

            end do

                  ic = 0
               do k = 1, irmax
                  nctemp(ic+k:ic+k) = ' '
               end do

                  ic = ic + irmax
                  nctemp(ic+1:ic+3) = 'reg'
                  ic = ic + 3

                  nctemp(ic+1:ic+15) = '      c1       '
                  ic = ic + 15
                  nctemp(ic+1:ic+15) = '      c2       '
                  ic = ic + 15
                  nctemp(ic+1:ic+15) = '      c3       '
                  ic = ic + 15
                  nctemp(ic+1:ic+15) = '      c4       '
                  ic = ic + 15

               write(iot,'(200a1)') (nctemp(k:k),k=1,ic)

*-----------------------------------------------------------------------

            do i = 1, mlrgn

                  write(nrtem1,'(i7)') melrg(i,1)
                  do k = 1, 7
                     if( nrtem1(k:k) .ne. ' ' ) goto 610
                  end do
  610             isl1 = 8 - k

               if( melrg(i,1) .eq. melrg(i,2) ) then

                  ic = irmax - isl1 + 3

                  do k = 1, ic
                     nctemp(k:k) = ' '
                  end do

                  do k = 1, isl1
                     nctemp(ic+k:ic+k) = nrtem1(7-isl1+k:7-isl1+k)
                  end do

                  ic = ic + isl1

               else

                  write(nrtem2,'(i7)') melrg(i,2)
                  do k = 1, 7
                     if( nrtem2(k:k) .ne. ' ' ) goto 620
                  end do
  620             isl2 = 8 - k

                  isll = isl1 + isl2 + 3

                  ic = irmax - isll + 3

                  do k = 1, ic
                     nctemp(k:k) = ' '
                  end do

                     nctemp(ic+1:ic+1) = '{'
                     ic = ic + 1

                  do k = 1, isl1
                     nctemp(ic+k:ic+k) = nrtem1(7-isl1+k:7-isl1+k)
                  end do
                     ic = ic + isl1

                     nctemp(ic+1:ic+1) = '-'
                     ic = ic + 1

                  do k = 1, isl2
                     nctemp(ic+k:ic+k) = nrtem2(7-isl2+k:7-isl2+k)
                  end do
                     ic = ic + isl2

                     nctemp(ic+1:ic+1) = '}'
                     ic = ic + 1

               end if

               write(nctemp(ic+1:ic+15),'(1p1g15.4)') elarg(i,1)
                  ic = ic + 15
               write(nctemp(ic+1:ic+15),'(1p1g15.4)') elarg(i,2)
                  ic = ic + 15
               write(nctemp(ic+1:ic+15),'(1p1g15.4)') elarg(i,3)
                  ic = ic + 15
               write(nctemp(ic+1:ic+15),'(1p1g15.4)') elarg(i,4)
                  ic = ic + 15

*-----------------------------------------------------------------------

               write(iot,'(200a1)') (nctemp(k:k),k=1,ic)

            end do

      end if

*-----------------------------------------------------------------------
*     Reg Name
*-----------------------------------------------------------------------

      if( ierrg .eq. 0 .and. mtrgn .gt. 0 ) then

            write(iot,'(/"[ Reg Name ]")')

                  inmax = 0
                  irmax = 0

            do i = 1, mtrgn

                  write(nrtem1,'(i7)') mtrg(i,1)
                  do k = 1, 7
                     if( nrtem1(k:k) .ne. ' ' ) goto 110
                  end do
  110             isl1 = 8 - k

               if( mtrg(i,1) .eq. mtrg(i,2) ) then

                  isll = isl1

               else

                  write(nrtem2,'(i7)') mtrg(i,2)
                  do k = 1, 7
                     if( nrtem2(k:k) .ne. ' ' ) goto 120
                  end do
  120             isl2 = 8 - k

                  isll = isl1 + isl2 + 3

               end if

                  if( isll .gt. irmax ) irmax = isll

                     inp = 0
                     ins = 0
                  do k = 1, nmtrg(i)
                     if( dmtrg(i)(k:k) .eq. ' ' ) ins = 2
                     if( dmtrg(i)(k:k) .eq. '{' ) inp = inp + 1
                     if( dmtrg(i)(k:k) .eq. '}' ) inp = inp + 1
                  end do

                     nmtmg = nmtrg(i) + ins + inp

                  if( nmtmg .gt. inmax ) inmax = nmtmg

            end do

                  ic = 0
               do k = 1, irmax
                  nctemp(ic+k:ic+k) = ' '
               end do

                  ic = ic + irmax
                  nctemp(ic+1:ic+3) = 'reg'
                  ic = ic + 3

               do k = 1, inmax
                  nctemp(ic+k:ic+k) = ' '
               end do

                  ic = ic + inmax
                  nctemp(ic+1:ic+4) = 'name'
                  ic = ic + 4

                  nctemp(ic+1:ic+9) = '     size'
                  ic = ic + 9

               write(iot,'(200a1)') (nctemp(k:k),k=1,ic)

*-----------------------------------------------------------------------

            do i = 1, mtrgn

                  write(nrtem1,'(i7)') mtrg(i,1)
                  do k = 1, 7
                     if( nrtem1(k:k) .ne. ' ' ) goto 210
                  end do
  210             isl1 = 8 - k

               if( mtrg(i,1) .eq. mtrg(i,2) ) then

                  ic = irmax - isl1 + 3

                  do k = 1, ic
                     nctemp(k:k) = ' '
                  end do

                  do k = 1, isl1
                     nctemp(ic+k:ic+k) = nrtem1(7-isl1+k:7-isl1+k)
                  end do

                  ic = ic + isl1

               else

                  write(nrtem2,'(i7)') mtrg(i,2)
                  do k = 1, 7
                     if( nrtem2(k:k) .ne. ' ' ) goto 220
                  end do
  220             isl2 = 8 - k

                  isll = isl1 + isl2 + 3

                  ic = irmax - isll + 3

                  do k = 1, ic
                     nctemp(k:k) = ' '
                  end do

                     nctemp(ic+1:ic+1) = '{'
                     ic = ic + 1

                  do k = 1, isl1
                     nctemp(ic+k:ic+k) = nrtem1(7-isl1+k:7-isl1+k)
                  end do
                     ic = ic + isl1

                     nctemp(ic+1:ic+1) = '-'
                     ic = ic + 1

                  do k = 1, isl2
                     nctemp(ic+k:ic+k) = nrtem2(7-isl2+k:7-isl2+k)
                  end do
                     ic = ic + isl2

                     nctemp(ic+1:ic+1) = '}'
                     ic = ic + 1

               end if

                  inp = 0
                  ins = 0
               do k = 1, nmtrg(i)
                  if( dmtrg(i)(k:k) .eq. ' ' ) ins = 2
                  if( dmtrg(i)(k:k) .eq. '{' ) inp = inp + 1
                  if( dmtrg(i)(k:k) .eq. '}' ) inp = inp + 1
               end do

                  nmtmg = nmtrg(i) + ins + inp

                  icc = inmax - nmtmg + 4

                  do k = 1, icc
                     nctemp(ic+k:ic+k) = ' '
                  end do

                     ic = ic + icc

               if( ins .gt. 0 ) then

                  nctemp(ic+1:ic+1) = '{'
                  ic = ic + 1

               end if

                        ii = 0
                  do jj = 1, nmtrg(i)
                     if( dmtrg(i)(jj:jj) .eq. '{' .or.
     &                   dmtrg(i)(jj:jj) .eq. '}' ) then
                        ii = ii + 1
                        nctemp(ic+ii:ic+ii) = yen
                        ii = ii + 1
                        nctemp(ic+ii:ic+ii) = dmtrg(i)(jj:jj)
                     else
                        ii = ii + 1
                        nctemp(ic+ii:ic+ii) = dmtrg(i)(jj:jj)
                     end if
                  end do

                  ic = ic + ii

               if( ins .gt. 0 ) then

                  nctemp(ic+1:ic+1) = '}'
                  ic = ic + 1

               end if

               write(nctemp(ic+1:ic+9),'(f9.2)') smtrg(i)

                  ic = ic + 9

*-----------------------------------------------------------------------

               write(iot,'(200a1)') (nctemp(k:k),k=1,ic)

            end do

      end if

*-----------------------------------------------------------------------
*     frag data
*-----------------------------------------------------------------------

      if( ierrg .eq. 0 .and. ifrgd .gt. 0 ) then

            write(iot,'(/"[ Frag Data ]")')

            write(iot,'("  opt    proj        targ         file")')

            do i = 1, ifrgd

C S.Hashimoto revised for the new [Frag Data]. (2014.12.26)
C projectile
               ityp = kftp(ifgdf(i,2))
               if ( ityp .eq. 19 ) then

                  iz   = ifgdf(i,2) / 1000000
                  ia   = ifgdf(i,2) - ifgdf(i,2) / 1000000 * 1000000
                  chat = elmnt(iz)

                  if( ia .gt. 0 ) then

                     if( chat(1:1) .eq. ' ' ) then

                        write(chau,'(i4,a1)') ia, chat(2:2)

                     else

                        write(chau,'(i3,a2)') ia, chat(1:2)

                     end if

                  else

                     if( chat(1:1) .eq. ' ' ) then

                        write(chau,'(4x,a1)') chat(2:2)

                     else

                        write(chau,'(3x,a2)') chat(1:2)

                     end if

                  end if

               else if ( ityp .le. 18 .and. ityp .ne. 11 ) then

                  write(chau,'(a8)') pname(ityp)

               else if ( ityp .eq. 11 .or. ityp .eq. 20 ) then

                  write(chau,'(i8)') ifgdf(i,2)

               end if


C target
               ityp = kftp(ifgdf(i,3))
               if ( ityp .eq. 19 ) then

                  iz   = ifgdf(i,3) / 1000000
                  ia   = ifgdf(i,3) - ifgdf(i,3) / 1000000 * 1000000
                  chat = elmnt(iz)

                  if( ia .gt. 0 ) then

                     if( chat(1:1) .eq. ' ' ) then

                        write(chav,'(i4,a1)') ia, chat(2:2)

                     else

                        write(chav,'(i3,a2)') ia, chat(1:2)

                     end if

                  else

                     if( chat(1:1) .eq. ' ' ) then

                        write(chav,'(4x,a1)') chat(2:2)

                     else

                        write(chav,'(3x,a2)') chat(1:2)

                     end if

                  end if

               else if ( ityp .le. 18 .and. ityp .ne. 11 ) then

                  write(chav,'(a8)') pname(ityp)

               else if ( ityp .eq. 11 .or. ityp .eq. 20 ) then

                  write(chav,'(i8)') ifgdf(i,2)

               end if

               write(iot,'(i5,4x,a8,4x,a8,5x,a)')
     &         ifgdf(i,1), chau, chav, frgfl(i)(1:ifgdf(i,4))

            end do

      end if

*-----------------------------------------------------------------------
*     Multiplier
*-----------------------------------------------------------------------

      if( ierrg .eq. 0 .and. imltp .gt. 0 ) then

         do kk = 1, imltp

            if( imdfl(kk) .eq. 0 .or. infout .eq. 8 ) then

c           write(iot,'(/"[ Multiplier ]")')
c
c           call echprt(1,iot,kk,impan,impat,jmpat,multmax,6,6)  ! T.Sato 2024/12/05, change from 100 to multmax
c
*-----------------------------------------------------------------------
c
c           write(iot,'("  number = ",i4)') idmlt(kk)
c
C MATSUDA 2024.11.25 (options: xlin, ylog, xlog, and ylin)
c           if( iimlt(kk) .eq. 1 ) then
c              write(iot,'("   interpolation = lin")')
c c         else if( iimlt(kk) .lt. -1 ) then
c           else if( iimlt(kk) .eq. -1 ) then
c              write(iot,'("   interpolation = log")')
c           else if( iimlt(kk) .eq. -2 ) then
c              write(iot,'("   interpolation = glow")')
c           else if( iimlt(kk) .eq.  2 ) then
c              write(iot,'("   interpolation = ghigh")')
c           else if( iimlt(kk) .eq.  3 ) then
c              write(iot,'("   interpolation = xlin")')
c           else if( iimlt(kk) .eq. -3 ) then
c              write(iot,'("   interpolation = xlog")')
c           else if( iimlt(kk) .eq.  4 ) then
c              write(iot,'("   interpolation = ylog")')
c           else if( iimlt(kk) .eq. -4 ) then
c              write(iot,'("   interpolation = ylin")')
c           end if
c
C MATSUDA 2024.11.25 (lagrange)
C           write(iot,'("  lagrange = ",i4)') ilmlt(kk)
c
c           write(iot,'("  ne = ",i4)') inmlt(kk)
c
c              igr = inmlt(kk)
c              igm = ismlt(kk)
c
c           do i = 1, igr
c
c              write(iot,'(1p2e13.5)')
c    &            gmsh_ismlt(igm-1+2*i-1), gmsh_ismlt(igm-1+2*i)
c
c           end do

               write(iot,'()')

            call multech(iot,kk,1,1)

            end if

         end do

      end if

*-----------------------------------------------------------------------
*     tally
*-----------------------------------------------------------------------

      if( itnm .gt. 0 .and. ierrg .eq. 0 ) then

         do m = 1, itnm

               write(iot,'()')

*-----------------------------------------------------------------------
*           t-track tally
*-----------------------------------------------------------------------

            if( ital(m) .eq. 1 ) then
              if( ital(m+1) .eq. 16 ) then

              else
                 call trckech(iot,m,0,0)
              end if

*-----------------------------------------------------------------------
*           t-adjoint tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 19 ) then

               call tadjech(iot,m,0,0)

*-----------------------------------------------------------------------
*           t-cross tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 2 ) then

               call tcrsech(iot,m,0,0)

*-----------------------------------------------------------------------
*           t-yield tally
*-----------------------------------------------------------------------


            else if( ital(m) .eq. 3 ) then
              if( ital(m+2) .eq. 16 ) then

              else
                 call tyilech(iot,m,0,0)
              end if



*-----------------------------------------------------------------------
*           t-dchain tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 16 ) then

               call tdchech(iot,m,0,0)

*-----------------------------------------------------------------------
*           t-heat tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 4 ) then

               call thetech(iot,m,0,0)

*-----------------------------------------------------------------------
*           t-star tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 5 ) then

               call tstaech(iot,m,0,0)

*-----------------------------------------------------------------------
*           t-time tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 6 ) then

               call ttimech(iot,m,0,0)

*-----------------------------------------------------------------------
*           t-dpa tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 7 ) then

               call tdpaech(iot,m,0,0)

*-----------------------------------------------------------------------
*           t-product tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 8 ) then

               call tproech(iot,m,0,0)

*-----------------------------------------------------------------------
*           t-gshow tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 9 ) then

               call tgshech(iot,m,0,0)

*-----------------------------------------------------------------------
*           t-rshow tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 10 ) then

               call trshech(iot,m,0,0)

*-----------------------------------------------------------------------
*           t-3dshow tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 11 ) then

               call tdshech(iot,m,0,0)

*-----------------------------------------------------------------------
*           LET tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 12 ) then

               call letech(iot,m,0,0)

*-----------------------------------------------------------------------
*           deposit tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 13 ) then

               call depstech(iot,m,0,0)

*-----------------------------------------------------------------------
*           deposit2 tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 14 ) then

               call dps2tech(iot,m,0,0)

*-----------------------------------------------------------------------
*           SED tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 15 ) then

               call sedech(iot,m,0,0)

*-----------------------------------------------------------------------
*           point tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 17 ) then

               call tponech(iot,m,0,0)

*-----------------------------------------------------------------------
*           wwg tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 18 ) then

               call twwgech(iot,m,0,0)

*-----------------------------------------------------------------------
*           user defined tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 20 ) then

               call usrdfech(iot,m,0,0)

*-----------------------------------------------------------------------
*           volume tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 21 ) then

               call tvolech(iot,m,0,0)

*-----------------------------------------------------------------------
*           wwbg  tally
*-----------------------------------------------------------------------

            else if( ital(m) .eq. 22 ) then

               call twbgech(iot,m,0,0)

            end if

*-----------------------------------------------------------------------

         end do

      end if

*-----------------------------------------------------------------------
*     end of input echo
*-----------------------------------------------------------------------

         if( ierrg .eq. 0 ) then

            if( ivers .ne. 0 .and. icntl .ne. 3 ) then

               write(iot,'(/"[END] of Input Echo <<<",56("="))')

            else

               write(iot,'(/"[END]")')

            end if

         end if

*-----------------------------------------------------------------------

*-----------------------------------------------------------------------

      if( infout .eq. 2 .or. infout .eq. 4 .or. infout .eq. 6 .or.
     &    infout .eq. 7 .or. infout .eq. 8 ) then
         iot = io
      else
         iot = ioi
      end if

*-----------------------------------------------------------------------
*     output c-values
*-----------------------------------------------------------------------

      if ( mstz(161) .ne. 0 ) then
         write(iot,'(/
     &              "$ c-values output")')
         do icvalue = 1, icvalout
            write(iot,'(3x,a1,i3,a11,1p1e16.6,a8,i9,a4,a)')
     &           'c',idint(cvalout(1,icvalue)),' was set to'
     &           ,cvalout(2,icvalue),' at line'
     &           ,idint(cvalout(3,icvalue)),' of '
     &           ,cfncvalout(icvalue)(1:idint(cvalout(4,icvalue)))
         end do
         call deallocate_cvalout
      end if

*-----------------------------------------------------------------------
*     write geometry memory report
*-----------------------------------------------------------------------

            mgtot = ngfini - ngstar

         if( ngfini - ngstar .gt. 0 .and. ierrg .eq. 0 ) then

                  write(iot,'(/
     &                   "<<< Report of real Geometry memory >>>"/
     &                   "*           GG/CG memory =",i9,/
     &                   "*     total memory: mdas =",i9,/
     &                   "*                percent =",f9.2," %")')
     &                     ngfini - ngstar,
     &                     mdas,
     &                     dble(ngfini-ngstar)/dble(mdas)*100.

            if( icgg .eq. 1 .and. icntl .eq. 3 ) then

                  write(iot,'(/
     &                   "<<< Report of temporary GG memory >>>"/
     &                   "*    temporary GG memory =",i9,/
     &                   "*     total memory: mdas =",i9,/
     &                   "*              ( percent =",f9.2," % )")')
     &                     ngfin0 - ngfini,
     &                     mdas,
     &                     dble(ngfin0-ngfini)/dble(mdas)*100.

            end if

         end if

*-----------------------------------------------------------------------
*     write material memory of high energy part
*-----------------------------------------------------------------------

         if( nmhigh .gt. 0 ) then

            write(iot,'(/"<<< Report of material memory",
     &                   " for high energy>>>"/
     &                   "*   high energy material =",i9,/
     &                   "*     total memory: mdas =",i9,/
     &                   "*                percent =",f9.2," %")')
     &                    nmhigh,
     &                    mdas,
     &                    dble(nmhigh)/dble(mdas)*100.

         end if

*-----------------------------------------------------------------------
*     write material memory of low energy part
*-----------------------------------------------------------------------

         if( nmlow .gt. 0 ) then

            write(iot,'(/"<<< Report of material memory",
     &                   " for low energy>>>"/
     &                   "*    low energy material =",i9,/
     &                   "*     total memory: mdas =",i9,/
     &                   "*                percent =",f9.2," %")')
     &                    nmlow,
     &                    mdas,
     &                    dble(nmlow)/dble(mdas)*100.

         end if

*-----------------------------------------------------------------------
*     write tally memory report
*-----------------------------------------------------------------------

               mttot = 0

         if( itnm .gt. 0 .and. ierrg .eq. 0 ) then

            do i = 1, itnm

               mttot = mttot + itsmn(i)

            end do

            write(iot,'(/"<<< Report of real tally memory >>>"/
     &                   "*      real tally memory =",i9,/
     &                   "*     total memory: mdas =",i9,/
     &                   "*                percent =",f9.2," %")')
     &                    mttot,
     &                    mdas,
     &                    dble(mttot)/dble(mdas)*100.

            if( itlmx .gt. 0 ) then

            write(iot,'(/"<<< Report of temporary tally memory >>>"/
     &                   "* temporary tally memory =",i9,/
     &                   "*     total memory: mdas =",i9,/
     &                   "*                percent =",f9.2," %")')
     &                    itlmx,
     &                    mdas,
     &                    dble(itlmx)/dble(mdas)*100.
            end if

         end if

*-----------------------------------------------------------------------
*     write bank memory
*-----------------------------------------------------------------------

         if( mbmax .gt. 0 .and. ierrg .eq. 0 ) then

            write(iot,'(/"<<< Report of bank memory >>>"/
     &                   "*            bank memory =",i9,/
     &                   "*     total memory: mdas =",i9,/
     &                   "*                percent =",f9.2," %")')
     &                    mbmax,
     &                    mdas,
     &                    dble(mbmax)/dble(mdas)*100.

         end if

*-----------------------------------------------------------------------
*     write other memory report
*-----------------------------------------------------------------------

            motot = mmmax - mgtot - mttot - itlmx - mbmax
     &            - nmhigh - nmlow

         if( motot .gt. 0 .and. ierrg .eq. 0 ) then

            write(iot,'(/"<<< Report of other memory >>>"/
     &                   "*           other memory =",i9,/
     &                   "*     total memory: mdas =",i9,/
     &                   "*                percent =",f9.2," %")')
     &                    motot,
     &                    mdas,
     &                    dble(motot)/dble(mdas)*100.

         end if

*-----------------------------------------------------------------------
*     write total memory report
*-----------------------------------------------------------------------

         if( ierrg .eq. 0 ) then

            write(iot,'(/"<<< Report of total memory >>>"/
     &                   "*      used total memory =",i9,/
     &                   "*     total memory: mdas =",i9,/
     &                   "*                percent =",f9.2," %"/
     &                   "*                ---------------------")')
     &                    mmmax,
     &                    mdas,
     &                    dble(mmmax)/dble(mdas)*100.

         end if

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*     write warnings from CG or GG
*-----------------------------------------------------------------------

               rewind iog

               lin = 0

  300       continue

            read(iog,'(a200)', iostat = ios ) chin
            if( ios .eq. -1 ) goto 301

               lin = lin + 1

            if( lin .eq. 1 ) then

               if( icgg .eq. 1 ) then

                  write(io,'(//"<<< Message from GG set up >>>"/)')
                  ErrCha = ''
                  MsgID = 'L:28743/R:echoi/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(//"<<< Message from GG set up >>>"/)')

               else

                  write(io,'(//"<<< Message from CG set up >>>"/)')
                  ErrCha = ''
                  MsgID = 'L:28751/R:echoi/F:read00.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(//"<<< Message from CG set up >>>"/)')

               end if

            end if

               call chlngt(chin,200,i1,i2)

               write(io,'(200a1)') (chin(i:i),i=1,i2)
               ErrCha = ''
               MsgID = 'L:28763/R:echoi/F:read00.f'
               call ErrWrite(MsgID, ErrCha)
               write(jo,'(200a1)') (chin(i:i),i=1,i2)

               goto 300

  301       close( iog )

*-----------------------------------------------------------------------
*     write warnings from MAT
*-----------------------------------------------------------------------

         if( nmode .ne. 0 .and. kmout .eq. 1 ) then

               rewind iom1

               lin = 0

  310       continue

            read(iom1,'(a200)', iostat = ios ) chin
            if( ios .eq. -1 ) goto 311

               lin = lin + 1

            if( lin .eq. 1 ) then

                  write(io,'(//"<<< Message from Material and",
     &                          " Data Libraries >>>")')


            end if

               call chlngt(chin,200,i1,i2)

               write(io,'(200a1)') (chin(i:i),i=1,i2)

               goto 310

  311       close( iom1 )

         end if

         close( ioi )

      call moddas_deallocate_cha(chrg)
      return
      end


************************************************************************
*                                                                      *
      subroutine echsrsa(j,iot,ierr)
*                                                                      *
*       echo dir in [source] with data                                 *
*       modified by K.Niita on 2005/11/08                              *
*                                                                      *
************************************************************************
      use moddas_source

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /isorfs/ ispfs(isrc), rspfn, rspfz, ispfn
      common /isorst/ jstyp(isrc), istyp(isrc), inkf0(isrc), lstyp(isrc)
!$OMP THREADPRIVATE(/isorst/)
      common /isoraa/ narp(isrc), naei(isrc), naea(isrc), nafe(isrc),
     &                naft(isrc), nall(isrc), napi(isrc), napw(isrc)

      common /isorfa/ isnn(isrc), lsfy(isrc), srfy(isrc)
      character srfy*200


      common /isoras/ sag1(isrc),sag2(isrc),jatyp(isrc),jqtyp(isrc)

*-----------------------------------------------------------------------
*     source energy group
*-----------------------------------------------------------------------

      if( jatyp(j) .eq.   1 .or.
     &    jatyp(j) .eq.   4 .or.
     &    jatyp(j) .eq.   5 .or.
     &    jatyp(j) .eq.   6 .or.
     &    jatyp(j) .eq.  11 .or.
     &    jatyp(j) .eq.  14 .or.
     &    jatyp(j) .eq.  15 .or.
     &    jatyp(j) .eq.  16 .or.
     &  ( jstyp(j) .eq.  17 .and. jatyp(j) .gt. 0 ) .or.
     &  ( jstyp(j) .eq. 100 .and. jatyp(j) .gt. 0 ) ) then

         if( jatyp(j) .eq. 1 ) then

            write(iot,'("   a-type = ",i3,12x,
     &      " # angular distribution given by cosine")') jatyp(j)

            write(iot,'("       na = ",i4,11x,
     &      " # number of cosine and weight"/27x,
     &      " # na>0: lin, na<0: log interpolation"/
     &      " #  data = ( a(i), w(i), i = 1, na ), amax")')
     &      narp(j) * nall(j)

            write(iot,'(1p6e13.5)') (agmin(naei(j)+i),
     &                               fagrp(nafe(j)+i),i=1,narp(j)),
     &                               agmax(naea(j)+narp(j))

         else if( jatyp(j) .eq. 11 ) then

            write(iot,'("   a-type = ",i3,12x,
     &      " # angular distribution given by theta(degree)")')
     &         jatyp(j)

            write(iot,'("       na = ",i4,11x,
     &      " # number of theta and weight"/
     &      " #  data = ( a(i), w(i), i = 1, na ), amax")')
     &      narp(j)

            write(iot,'(1p6e13.5)') (agmin(naei(j)+i),
     &                               fagrp(nafe(j)+i),i=1,narp(j)),
     &                               agmax(naea(j)+narp(j))

         else if( jatyp(j) .eq. 4 ) then

            write(iot,'("   a-type = ",i3,12x,
     &      " # angular distribution given by cosine")') jatyp(j)

            write(iot,'("       na = ",i4,11x,
     &      " # number of cosine and weight"/27x,
     &      " # na>0: lin, na<0: log interpolation"/
     &      " #  data = ( a(i), w(i), i = 1, na ), amax")')
     &      narp(j) * nall(j)

            write(iot,'(1p6e13.5)') (agmin(naei(j)+i),
     &                               fagrp(nafe(j)+i),i=1,narp(j)),
     &                               agmax(naea(j)+narp(j))

            write(iot,'("   q-type = ",i3,12x,
     &      " # generate weight type, 0:equal, 1: given by data")')
     &         jqtyp(j)

          if( jqtyp(j) .eq. 1 ) then
            write(iot,'(1p6e13.5)') (paw(napi(j)+i),i=1,narp(j))
          end if

         else if( jatyp(j) .eq. 14 ) then

            write(iot,'("   a-type = ",i3,12x,
     &      " # angular distribution given by theta(degree)")')
     &           jatyp(j)

            write(iot,'("       na = ",i4,11x,
     &      " # number of theta and weight"/
     &      " #  data = ( a(i), w(i), i = 1, na ), amax")')
     &      narp(j)

            write(iot,'(1p6e13.5)') (agmin(naei(j)+i),
     &                               fagrp(nafe(j)+i),i=1,narp(j)),
     &                               agmax(naea(j)+narp(j))

            write(iot,'("   q-type = ",i3,12x,
     &      " # generate weight type, 0:equal, 1: given by data")')
     &         jqtyp(j)

          if( jqtyp(j) .eq. 1 ) then
            write(iot,'(1p6e13.5)') (paw(napi(j)+i),i=1,narp(j))
          end if

         else if( jatyp(j) .eq. 5 ) then

            write(iot,'("   a-type = ",i3,12x,
     &      " # Functional angular distribution by cosine")') jatyp(j)

            write(iot,'("     g(x) = ",1024a1)')
     &                 ( srfy(j)(k:k), k = 1, lsfy(j) )

            write(iot,'("      ag1 = ",1p1g12.5,3x,
     &      " # minimum cutoff cosine")') sag1(j)

            write(iot,'("      ag2 = ",1p1g12.5,3x,
     &      " # maximum cutoff cosine")') sag2(j)

            write(iot,'("       nn = ",i4,11x,
     &      " # number of cosine mesh"/27x,
     &      " # nn>0: lin, nn<0: log mesh and interpolation")')
     &      narp(j) * nall(j)

         else if( jatyp(j) .eq. 6 ) then

            write(iot,'("   a-type = ",i3,12x,
     &      " # Functional angular distribution by theta")') jatyp(j)

            write(iot,'("     g(x) = ",1024a1)')
     &                 ( srfy(j)(k:k), k = 1, lsfy(j) )

            write(iot,'("      ag1 = ",1p1g12.5,3x,
     &      " # minimum cutoff theta(degree)")') sag1(j)

            write(iot,'("      ag2 = ",1p1g12.5,3x,
     &      " # maximum cutoff theta(degree)")') sag2(j)

            write(iot,'("       nn = ",i4,11x,
     &      " # number of theta mesh"/27x,
     &      " # nn>0: lin, nn<0: log mesh and interpolation")')
     &      narp(j) * nall(j)

            write(iot,'("   q-type = ",i3,12x,
     &      " # generate weight type, 0:equal, 1: given by data")')
     &         jqtyp(j)

          if( jqtyp(j) .eq. 1 ) then
            write(iot,'(1p6e13.5)') (paw(napi(j)+i),i=1,narp(j))
          end if

         else if( jatyp(j) .eq. 15 ) then

            write(iot,'("   a-type = ",i3,12x,
     &      " # Functional angular distribution by cosine")') jatyp(j)

            write(iot,'("     g(x) = ",1024a1)')
     &                 ( srfy(j)(k:k), k = 1, lsfy(j) )

            write(iot,'("      ag1 = ",1p1g12.5,3x,
     &      " # minimum cutoff cosine")') sag1(j)

            write(iot,'("      ag2 = ",1p1g12.5,3x,
     &      " # maximum cutoff cosine")') sag2(j)

            write(iot,'("       nn = ",i4,11x,
     &      " # number of cosine mesh"/27x,
     &      " # nn>0: lin, nn<0: log mesh and interpolation")')
     &      narp(j) * nall(j)

         else if( jatyp(j) .eq. 16 ) then

            write(iot,'("   a-type = ",i3,12x,
     &      " # Functional angular distribution by theta")') jatyp(j)

            write(iot,'("     g(x) = ",1024a1)')
     &                 ( srfy(j)(k:k), k = 1, lsfy(j) )

            write(iot,'("      ag1 = ",1p1g12.5,3x,
     &      " # minimum cutoff theta(degree)")') sag1(j)

            write(iot,'("      ag2 = ",1p1g12.5,3x,
     &      " # maximum cutoff theta(degree)")') sag2(j)

            write(iot,'("       nn = ",i4,11x,
     &      " # number of theta mesh"/27x,
     &      " # nn>0: lin, nn<0: log mesh and interpolation")')
     &      narp(j) * nall(j)

            write(iot,'("   q-type = ",i3,12x,
     &      " # generate weight type, 0:equal, 1: given by data")')
     &         jqtyp(j)

          if( jqtyp(j) .eq. 1 ) then
            write(iot,'(1p6e13.5)') (paw(napi(j)+i),i=1,narp(j))
          end if

         else

          write(ErrCha,'("*** Error : a-type is wrong in [source]")')
          ErrID = 'L:29028/R:echsrsa/F:read00.f' !E03_004_001
          call ErrWrite(ErrID,ErrCha)
          ierr = ierr + 1

         end if

      end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine echsrst(j,iot,ierr)
*                                                                      *
*       echo time in [source] with data                                *
*       modified by K.Niita on 2015/02/09                              *
*                                                                      *
************************************************************************
      use moddas_source

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /isorst/ jstyp(isrc), istyp(isrc), inkf0(isrc), lstyp(isrc)
!$OMP THREADPRIVATE(/isorst/)

      common /isorsm/ stm0(isrc),stmw(isrc),stmc(isrc),stmd(isrc),
     &                jttyp(isrc),jttpn(isrc)

      common /isortt/ ntrp(isrc), ntei(isrc), ntea(isrc), ntfe(isrc),
     &                ntft(isrc), ntll(isrc), ntpi(isrc), ntpw(isrc)

      common /isorft/ isll(isrc), lsfz(isrc), srfz(isrc)
      character srfz*200


      common /isorts/ stg1(isrc),stg2(isrc),jotyp(isrc)

*-----------------------------------------------------------------------

      if( jttyp(j) .eq.   1 .or.
     &    jttyp(j) .eq.   2 .or.
     &    jttyp(j) .eq.   3 .or.
     &    jttyp(j) .eq.   4 .or.
     &    jttyp(j) .eq.   5 .or.
     &    jttyp(j) .eq.   6 .or.
     &    jttyp(j) .eq. 100 .or.
     &  ( jstyp(j) .eq.  17 .and. jttyp(j) .gt. 0 ) .or.
     &  ( jstyp(j) .eq. 100 .and. jttyp(j) .gt. 0 ) ) then

*-----------------------------------------------------------------------

         if( jttyp(j) .eq. 1 .or. jttyp(j) .eq. 2 ) then

               write(iot,'("   t-type = ",i3,12x,
     &         " # (D=0) 0:t=0, 1:rectangular, 2:Gaussian")') jttyp(j)

               write(iot,'("       t0 = ",1p1g12.5,3x,
     &         " # (D=0.0) center time of first pulse (ns)")') stm0(j)

               write(iot,'("       tw = ",1p1g12.5,3x,
     &         " # FWHM of one pulse [ns]")') stmw(j)

               write(iot,'("       tn = ",i4,11x,
     &         " # number of time pulse")') jttpn(j)

               write(iot,'("       td = ",1p1g12.5,3x,
     &         " # time between pulses (ns)")') stmd(j)

            if( jttyp(j) .eq. 2 ) then

               write(iot,'("       tc = ",1p1g12.5,3x,
     &         " # (D=10*tw) cut off time Gaussian pulse (ns)")')
     &               stmc(j)

            end if

         else if( jttyp(j) .eq. 3 ) then

            write(iot,'("   t-type = ",i3,12x,
     &      " # time distribution given by data")') jttyp(j)

            write(iot,'("      ntt = ",i4,11x,
     &      " # number of time and weight"/27x,
     &      " # ntt>0: lin, ntt<0: log interpolation"/
     &      " #  data = ( t(i), w(i), i = 1, ntt ), tmax")')
     &      ntrp(j) * ntll(j)

            write(iot,'(1p6e13.5)') (tgmin(ntei(j)+i),
     &                               ftgrp(ntfe(j)+i),i=1,ntrp(j)),
     &                               tgmax(ntea(j)+ntrp(j))

         else if( jttyp(j) .eq. 4 ) then

            write(iot,'("   t-type = ",i3,12x,
     &      " # time distribution given by data")') jttyp(j)

            write(iot,'("      ntt = ",i4,11x,
     &      " # number of time and weight"/27x,
     &      " # ntt>0: lin, ntt<0: log interpolation"/
     &      " #  data = ( t(i), w(i), i = 1, ntt ), tmax")')
     &      ntrp(j) * ntll(j)

            write(iot,'(1p6e13.5)') (tgmin(ntei(j)+i),
     &                               ftgrp(ntfe(j)+i),i=1,ntrp(j)),
     &                               tgmax(ntea(j)+ntrp(j))

            write(iot,'("   o-type = ",i3,12x,
     &      " # generate weight type, 0:equal, 1: given by data")')
     &         jotyp(j)

          if( jotyp(j) .eq. 1 ) then
            write(iot,'(1p6e13.5)') (ptw(ntpi(j)+i),i=1,ntrp(j))
          end if

         else if( jttyp(j) .eq. 5 ) then

            write(iot,'("   t-type = ",i3,12x,
     &      " # Functional time distribution")') jttyp(j)

            write(iot,'("     h(x) = ",1024a1)')
     &                 ( srfz(j)(k:k), k = 1, lsfz(j) )

            write(iot,'("      tg1 = ",1p1g12.5,3x,
     &      " # minimum cutoff time")') stg1(j)

            write(iot,'("      tg2 = ",1p1g12.5,3x,
     &      " # maximum cutoff time")') stg2(j)

            write(iot,'("       ll = ",i4,11x,
     &      " # number of time mesh"/27x,
     &      " # ll>0: lin, ll<0: log mesh and interpolation")')
     &      ntrp(j) * ntll(j)

         else if( jttyp(j) .eq. 6 ) then

            write(iot,'("   t-type = ",i3,12x,
     &      " # Functional time distribution")') jttyp(j)

            write(iot,'("     h(x) = ",1024a1)')
     &                 ( srfz(j)(k:k), k = 1, lsfz(j) )

            write(iot,'("      tg1 = ",1p1g12.5,3x,
     &      " # minimum cutoff time")') stg1(j)

            write(iot,'("      tg2 = ",1p1g12.5,3x,
     &      " # maximum cutoff time")') stg2(j)

            write(iot,'("       ll = ",i4,11x,
     &      " # number of time mesh"/27x,
     &      " # ll>0: lin, ll<0: log mesh and interpolation")')
     &      ntrp(j) * ntll(j)

            write(iot,'("   o-type = ",i3,12x,
     &      " # generate weight type, 0:equal, 1: given by data")')
     &         jotyp(j)

          if( jotyp(j) .eq. 1 ) then
            write(iot,'(1p6e13.5)') (ptw(ntpi(j)+i),i=1,ntrp(j))
          end if

         else if( jttyp(j) .eq. 100 ) then

            write(iot,'("   t-type = ",i3,12x,
     &      " # Special time distribution")') jttyp(j)

            write(iot,'("      tg1 = ",1p1g12.5,3x,
     &      " # minimum cutoff time")') stg1(j)

            write(iot,'("      tg2 = ",1p1g12.5,3x,
     &      " # maximum cutoff time")') stg2(j)

         else

          write(ErrCha,'("*** Error : t-type is wrong in [source]")')
          ErrID = 'L:29210/R:echsrst/F:read00.f' !E03_005_001
          call ErrWrite(ErrID,ErrCha)
          ierr = ierr + 1

         end if

      end if

*-----------------------------------------------------------------------

      return
      end


