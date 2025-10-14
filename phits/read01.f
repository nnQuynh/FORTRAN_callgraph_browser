************************************************************************
*                                                                      *
      subroutine read01(in,io,ict,ierr)
*                                                                      *
*       read old input data                                            *
*       modified by K.Niita on 2005/11/02                              *
*                                                                      *
*                                                                      *
************************************************************************
      use partmod ! frtati 2021/10/05
      use moddas
      use moddas_character
      use moddas_mesh
      use moddas_region
      use moddas_source
      use moddas_tally

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param01.inc'

      parameter ( ibmt = 23 )

*-----------------------------------------------------------------------

      common /kmat1a/ mxmat, mxmat0, mxnel
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

      common /kmat1o/ iom1, iom2, iom3

*-----------------------------------------------------------------------

      common /inpec/  ititl, ipara, ibody, iregn, llarr, itby, itar
      common /htitl/  iclgt(100), ctitl(100)
      character       ctitl*200

      common /parai/  ipsq(400)
      common /paraj/  mstz(300), parz(300)

      common /regdu/  iuni(kvlmax)
      common /regdm/  idmg(kvlmax)
      common /regda/  ichl(kvlmax), chsm(kvlmax), ichmx, iod
      character       chsm*10
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /regdb/  nrsq, irsq(10)

      common /volmsg/ rvols(kvlmax), mnvol, nvols(kvlmax)
      common /regdd/  ivolm, iimpo
      common /impreg/ dimp(kvlmax)

*-----------------------------------------------------------------------

      common /isomul/ smlwt(isrc), totfact, imsrc
      common /isorst/ jstyp(isrc), istyp(isrc), inkf0(isrc), lstyp(isrc)
!$OMP THREADPRIVATE(/isorst/)
      common /isorsg/ seg0(isrc),seg1(isrc),seg2(isrc),seg3(isrc),
     &                set0(isrc),set1(isrc),set2(isrc),set3(isrc),
     &                jetyp(isrc),jptyp(isrc)
      common /isorsp/ sx0(isrc), sy0(isrc), sz0(isrc), sx1(isrc),
     &                sy1(isrc), sz1(isrc), sr0(isrc), se0(isrc),
     &                sdir(isrc), srx(isrc), sry(isrc), swem(isrc),
     &                sphi(isrc), sdom(isrc), swt0(isrc)
      common /isorsf/ isorf(isrc),lsfile(isrc), sfile(isrc)
      character sfile*100

      common /isorse/ ngrp(isrc), ngei(isrc), ngea(isrc), ngfe(isrc),
     &                ngft(isrc), ngll(isrc), ngpi(isrc), ngpw(isrc)

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

      common /tall17/ itndy(itlmax)
      common /tall18/ ithet(itlmax)

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall21/ rtfac(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall26/ itrcs(itlmax), itrcc(itlmax), itrss(itlmax)
      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall31/ iterl(itlmax)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      common /regcrs/ icrsflx

*-----------------------------------------------------------------------

      character hl*80, hl2*80
      dimension italxx(8), imesxx(6)

      data italxx /8*0/

      dimension ntcr(kvlmax,2)
      dimension stca(kvlmax)
      dimension ntrg(kvlmax)
      dimension ntyg(kvlmax)
      dimension nthg(kvlmax)
      dimension nttg(kvlmax)
      dimension steg(kvlmax)
      dimension strg(kvlmax)
      dimension stzr(kvlmax)
      dimension stxg(kvlmax)
      dimension styg(kvlmax)
      dimension stzg(kvlmax)
      dimension strs(kvlmax)
      dimension stzs(kvlmax)

      dimension lydn(0:100), lydz(100), lyda(100)
      dimension ikzz(maxpt,maxnt), iknn(maxpt,maxnt)

      dimension iptyp(6), ipnkf(6)

      data iptyp /   1,   2,   3,   5,   6,   7/
      data ipnkf /2212,2112, 211,-211,  13, -13/

*-----------------------------------------------------------------------

      character dam*1
      character ia*80

      character chin*200, chlw*200

      character chtit*60
      dimension ipva(4)
      dimension irby(10)

      dimension bval(50)

      dimension ibva(50), iblg(50)
      character chbd(50)*3

      dimension ibck(kvmmax)
      data ibck / kvmmax*0 /

      data ( chbd(i), i = 1, ibmt ) /
     &         'arb','sph','rcc','rec','trc',
     &         'ell','box','wed','rpp','gel',
     &         'tor','qua','bpp','wpp','p  ',
     &         'px ','py ','pz ','ps ','c  ',
     &         'cx ','cy ','cz '/

      data ( iblg(i), i = 1, ibmt ) /
     &           3,    3,    3,    3,    3,
     &           3,    3,    3,    3,    3,
     &           3,    3,    3,    3,    2,
     &           2,    2,    2,    2,    2,
     &           2,    2,    2/

      data ( ibva(i), i = 1, ibmt ) /
     &          30,    4,    7,   12,    8,
     &           7,   12,   12,    6,   12,
     &           9,   10,    9,    9,    4,
     &           1,    1,    1,    7,    7,
     &           4,    4,    4/

*-----------------------------------------------------------------------

      data irreg /0/
      data iyreg /0/
      data ihreg /0/
      data ikreg /0/

      data iengm  /0/
      data iengn  /0/
      data ksunit /-1/
      data kmunit /-1/
      data kpsurf /3/
      data kpmesh /3/
      data nrsurf /0/
      data nzsurf /0/
      data kyield /0/
      data lyield /0/
      data lydn /101*0/
      data nbsurf /0/
      data nrmshr /0/
      data nzmshr /0/
      data nxmshx /0/
      data nymshx /0/
      data nzmshx /0/

*-----------------------------------------------------------------------

      include 'err.inc'

*-----------------------------------------------------------------------

      imsrc   = 1
      totfact = 1.d0
      smlwts  = 1.d0

*-----------------------------------------------------------------------

            mci  = ( mmmax - 1 ) * 8 + 1
            mcmx = ( mdas  - 1 ) * 8 + 1 - mci

            iod = 26

            open(iod,status='scratch',form='unformatted')

            ichmx = 0

*-----------------------------------------------------------------------
*        old tally
*-----------------------------------------------------------------------

         if( ict .eq. 1 ) goto 9000

*-----------------------------------------------------------------------
*        icntl
*-----------------------------------------------------------------------

               ipara = ipara + 1
               ipsq(ipara) = 1

*-----------------------------------------------------------------------
*     card 1 : main title up to 80a
*-----------------------------------------------------------------------

         read(in,3) ia
    3    format(a80)

               call chnumb(ia,80,i1,i2)

               ititl = 1
               iclgt(ititl) = i2
               ctitl(ititl)(1:i2) = ia(1:i2)

*-----------------------------------------------------------------------
*     card 2 : sub title up to 80a
*-----------------------------------------------------------------------

         read(in,3)  ia

               call chnumb(ia,80,i1,i2)

               ititl = 2
               iclgt(ititl) = i2
               ctitl(ititl)(1:i2) = ia(1:i2)

*-----------------------------------------------------------------------
*     card 3 : random number
*-----------------------------------------------------------------------

         read(in,*) randkk,irskip

               parz( 21 ) = randkk

               ipara = ipara + 1
               ipsq(ipara) = -21

               mstz( 2 ) = irskip

               ipara = ipara + 1
               ipsq(ipara) = 2

*-----------------------------------------------------------------------
*     card 4 : emax,emin(1),emin(2),mxmat,maxcas,maxbch,n1col
*-----------------------------------------------------------------------

         read(in,*) emax,emin1,emin2,mxmat,maxcas,maxbch,n1col

               mstz( 3 ) = maxcas

               ipara = ipara + 1
               ipsq(ipara) = 3

               mstz( 4 ) = maxbch

               ipara = ipara + 1
               ipsq(ipara) = 4

               parz( 1 ) = emin1

               ipara = ipara + 1
               ipsq(ipara) = -1

               parz( 2 ) = emin2

               ipara = ipara + 1
               ipsq(ipara) = -2

               ipara = ipara + 1
               ipsq(ipara) = -14

*-----------------------------------------------------------------------
*     card 5 : nquit,neutp,nbertp,npowr2,npidk,nhstp
*-----------------------------------------------------------------------

         read(in,*) nquit,neutp,nbertp,npowr2,npidk,nhstp

            if( npidk .ne. mstz( 5 ) ) then

               mstz( 5 ) = npidk

               ipara = ipara + 1
               ipsq(ipara) = 5

            end if

*-----------------------------------------------------------------------
*     card 6 : andt,ctofe,nbogus,nspred,nwsprd,nseudo
*-----------------------------------------------------------------------

         read(in,*) andt,ctofe,nbogus,nspred,nwsprd,nseudo

            if( andt .ne. parz( 22 ) ) then

               parz( 22 ) = andt

               ipara = ipara + 1
               ipsq(ipara) = -22

            end if

            if( nbogus .eq. 0 ) then

               mstz( 6 ) = nbogus

               ipara = ipara + 1
               ipsq(ipara) = 6

            end if

               mstz( 7 ) = nspred

               ipara = ipara + 1
               ipsq(ipara) = 7

            if( nwsprd .ne. mstz( 14 ) ) then

               mstz( 14 ) = nwsprd

               ipara = ipara + 1
               ipsq(ipara) = 14

            end if

*-----------------------------------------------------------------------
*     card 7 : ielas,icasc,iqstep,lvlopt,igamma
*-----------------------------------------------------------------------

         read(in,*) ielas,icasc,iqstep,lvlopt,igamma

            if( icasc .lt. -100 ) then

               icasc = - icasc - 100
               icugn = 2

            else if( icasc .lt. 0 ) then

               icasc = - icasc
               icugn = 1

            else

               icugn = 0

            end if

               mstz( 13 ) = icugn

               ipara = ipara + 1
               ipsq(ipara) = 13

            if( ielas .gt. 10 .and. icasc .gt. 10 ) then

               mstz( 15 ) = ielas - 10

               ipara = ipara + 1
               ipsq(ipara) = 15

            else

               mstz( 8 ) = ielas

               ipara = ipara + 1
               ipsq(ipara) = 8

            end if


            if( icasc .gt. 10 ) then

               mstz( 1 ) = 1

               icasc = icasc - 10

            else

               mstz( 1 ) = 0

            end if

            if( icasc .eq. 1 ) then


            else if( icasc .eq. 2 ) then

               mstz( 9 ) = 1

               ipara = ipara + 1
               ipsq(ipara) = 9

               parz( 25 ) = 1000.0

               ipara = ipara + 1
               ipsq(ipara) = -25

            else if( icasc .eq. 3 ) then

               parz( 23 ) = 1000.0
               parz( 24 ) = 1000.0

               ipara = ipara + 1
               ipsq(ipara) = -23
               ipara = ipara + 1
               ipsq(ipara) = -24

               icasc = 1

            else if( icasc .eq. 4 ) then

               parz( 23 ) = 500.0
               parz( 24 ) = 500.0

               ipara = ipara + 1
               ipsq(ipara) = -23
               ipara = ipara + 1
               ipsq(ipara) = -24

               icasc = 1

            else if( icasc .eq. 5 ) then

               parz( 23 ) = 150.0
               parz( 24 ) = 150.0

               ipara = ipara + 1
               ipsq(ipara) = -23
               ipara = ipara + 1
               ipsq(ipara) = -24

               icasc = 1

            else if( icasc .eq. 6 ) then

               parz( 23 ) = 0.0
               parz( 24 ) = 0.0

               ipara = ipara + 1
               ipsq(ipara) = -23
               ipara = ipara + 1
               ipsq(ipara) = -24

               icasc = 1

            else

               write(io,'(''*** ERROR : value of icasc is wrong'',
     &                    '' in old input. icasc = '',i3)') icasc
               goto 888

            end if

         if( nbogus .ne. 0 ) then

            if( iqstep .eq. 2 ) then

               mstz( 6 ) = 1

               ipara = ipara + 1
               ipsq(ipara) = 6

            else if( iqstep .eq. 3 ) then

               mstz( 6 ) = 1

               ipara = ipara + 1
               ipsq(ipara) = 6

               mstz( 10 ) = 1

               ipara = ipara + 1
               ipsq(ipara) = 10

            else if( iqstep .eq. 4 ) then

               mstz( 6 ) = 2

               ipara = ipara + 1
               ipsq(ipara) = 6

            else if( iqstep .eq. 5 ) then

               mstz( 6 ) = 1

               ipara = ipara + 1
               ipsq(ipara) = 6

            else if( iqstep .eq. 6 ) then

               mstz( 6 ) = 3

               ipara = ipara + 1
               ipsq(ipara) = 6

            end if

         end if

            if( lvlopt .ne. mstz(11) ) then

               mstz( 11 ) = lvlopt

               ipara = ipara + 1
               ipsq(ipara) = 11

            end if

               mstz( 12 ) = igamma

               ipara = ipara + 1
               ipsq(ipara) = 12

*-----------------------------------------------------------------------

               ipara = ipara + 1
               ipsq(ipara) = 16

               ipara = ipara + 1
               ipsq(ipara) = 17

               ipara = ipara + 1
               ipsq(ipara) = 18

*-----------------------------------------------------------------------
*     card 8a : ( denh(m), nel1, m = 1, mxmat )
*          8b : ( zz(l,m),a(l,m),den(l,m), l = 1, nel1 )
*-----------------------------------------------------------------------

               iom1 = 23

               open(iom1,form='unformatted',status='scratch')

               libh = 0
               igas = 0
               istp = 0
               inlb = 0
               iplb = 0
               ielb = 0
               icnd = 0
               imts = -1
               libi  = 0

         do m = 1, mxmat

               read(in,*) denh, nel1

               write(iom1) nel1
               write(iom1) denh, libh
               write(iom1) igas, istp, inlb, iplb, ielb, icnd
               write(iom1) imts

               idnm(m) = m
               idmn(m) = m

            do l = 1, nel1

               read(in,*) zzin, ain, denin

               write(iom1) nint(zzin), nint(ain), denin, libi

            end do

         end do

*-----------------------------------------------------------------------
*        open temporary binary file : unit 18
*-----------------------------------------------------------------------

            itby = 18

            open(itby,form='unformatted',status='scratch')

*-----------------------------------------------------------------------
*     card 9 : title for body
*-----------------------------------------------------------------------

            read(in,'(a60)') chtit
            write(itby) chtit

*-----------------------------------------------------------------------
*     card 10 : parameters for body
*-----------------------------------------------------------------------

            read(in,*) ipva(1), ipva(2), ipva(3), ipva(4)

               if( ipva(2) .gt. 0 ) ipva(2) = 2

            write(itby) ipva

            if( ipva(3) .eq. 0 ) then

               nbdy    = 2
               irby(1) = 1
               irby(2) = 100

            else

               nbdy    = 3
               irby(1) = 1
               irby(2) = 2
               irby(3) = 100

            end if

*-----------------------------------------------------------------------
*     card 11 : body information
*        read one line from in upto [ end ]
*-----------------------------------------------------------------------

  140 continue

            read(in,'(a200)') chin

                  call chlngt(chin,200,i1,i2)

                  chlw = chin
                  call chcaps(chlw,i1,i2,i3,'#!$')

                  if( chlw(i1:i1+2) .eq. 'end' ) goto 150

*-----------------------------------------------------------------------

               ibody = ibody + 1

               ic2 = i1

         do k = 1, nbdy

               ic = jnumc(chlw,ic2,i3)

            if( irby(k) .eq. 1 ) then

               do i = 1, ibmt

                  il = ic + iblg(i) - 1

                  if( chlw(ic:il) .eq. chbd(i)(1:iblg(i)) ) goto 301

               end do

                  goto 999

  301          ipm = i
               ic2 = il + 1

            else if( irby(k) .eq. 2 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               ibnum = nint( cvvv )


            else if( irby(k) .eq. 100 ) then

               do i = 1, ibva(ipm)

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  bval(i) = cvvv

                  ic = ic2

                  if( i .lt. ibva(ipm) .and. ic .ge. i3 ) then

                        read(in,'(a200)') chin

                        call chlngt(chin,200,i1,i2)

                        chlw = chin
                        call chcaps(chlw,i1,i2,i3,'#!$')

                     ic = 1

                  end if

               end do

            end if

         end do

               if( ipva(3) .eq. 0 ) ibnum = ibody

               if( ibnum .le. 0 .or. ibnum .gt. 9999 ) goto 999

               if( ibck(ibnum) .ne. 0 ) goto 999
               ibck(ibnum) = ibody

               write(itby) chbd(ipm), ibnum,
     &                    ibva(ipm), ( bval(i), i = 1, ibva(ipm) )

               goto 140

*-----------------------------------------------------------------------

  999 continue

            write(io,'(''*** Error : body description is wrong'',
     &                 '' in old input'')')
            write(io,'(a200)') chin

            goto 888

*-----------------------------------------------------------------------

  150 continue

*-----------------------------------------------------------------------
*     card 12 : region definition
*-----------------------------------------------------------------------

               nrsq    = 2
               irsq(1) = 4
               irsq(2) = 100

*-----------------------------------------------------------------------

               call moddas_allocate_cha(MAX_NUM_CHRG, chrg)
               mci = 0
  160 continue

            read(in,'(a200)') chin

                  call chlngt(chin,200,i1,i2)

                  chlw = chin
                  call chcaps(chlw,i1,i2,i3,'#!$')


                  if( chlw(i1:i1+2) .eq. 'end' ) then

                     ichmx = max(ichmx,ichl(iregn))

                     write(iod) (chrg(mci+l:mci+l),l=1,ichl(iregn))

                     goto 170

                  end if

*-----------------------------------------------------------------------

               iregn = iregn + 1

               if( iregn .gt. kvlmax ) goto 998

               ic2  = i1

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( irsq(k) .eq. 4 ) then

               if( chlw(ic:ic) .lt. 'a' .or.
     &             chlw(ic:ic) .gt. 'z' .or.
     &             chlw(ic:ic+2) .eq. 'or ' ) then

                  iregn = iregn - 1

                  if( ichl(iregn) + i3 - ic + 2 .gt. mcmx ) goto 990

                  ichs = ichl(iregn) + 1
                  chrg(mci+ichs:mci+ichs) = ' '
                  ichc = i3 - ic + 1
                  ichl(iregn) = ichl(iregn) + ichc + 1

                  do i = 1, ichc

                     j = i + ichs
                     chrg(mci+j:mci+j) = chin(ic+i-1:ic+i-1)

                  end do

                     goto 160

               else

                  ic2 = inumc(chlw,ic,i3,' ')
                  chsm(iregn) = chin(ic:ic2-1)

                  if( iregn .gt. 1 ) then

                     ichmx = max(ichmx,ichl(iregn-1))

                     write(iod) (chrg(mci+l:mci+l),l=1,ichl(iregn-1))

                  end if

               end if

            else if( irsq(k) .eq. 100 ) then

               ichl(iregn) = i3 - ic + 1

               do i = 1, ichl(iregn)

                  chrg(mci+i:mci+i) = chin(ic+i-1:ic+i-1)

               end do

            end if

         end do

            goto 160

*-----------------------------------------------------------------------

  998 continue

         write(io,'('' ** ERROR : Number of region exceeds'',
     &              '' kvlmax in old input'')')
         goto 888

*-----------------------------------------------------------------------

  170 continue
         if( ichmx > MAX_NUM_CHRG ) then
            write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &           'sub.read01@read01.f'
     &              //' ?dimension over chrg?'
     &              //' ichmx > MAX_NUM_CHRG'
     &           ,' (ichmx=',ichmx,')'
     &           ,' (MAX_NUM_CHRG@moddas.f=',MAX_NUM_CHRG,')'
            ErrID = 'L:889/R:read01/F:read01.f'
            call ErrWrite(ErrID,ErrCha)
         endif
         call moddas_deallocate_cha(chrg)

*-----------------------------------------------------------------------
*     card 13 : region number definition
*-----------------------------------------------------------------------

         read(in,*) ( idrg(i),i = 1, iregn )

            do i = 1, iregn

               idgr(idrg(i)) = i

            end do

*-----------------------------------------------------------------------
*     card 14 : universe definition
*-----------------------------------------------------------------------

         read(in,*) ( iuni(i),i = 1, iregn )

               iuniv = 0

            do i = 1, iregn

               if( iuni(i) .ne. 0 ) iuniv = 1

            end do

            if( iuniv .gt. 0 ) then

               nrsq    = 5
               irsq(1) = 1
               irsq(2) = 2
               irsq(3) = 4
               irsq(4) = 5
               irsq(5) = 100

            end if

*-----------------------------------------------------------------------
*     card 15 : material definition
*-----------------------------------------------------------------------

         read(in,*) ( idmg(i),i = 1, iregn )

            do i = 1, iregn

               if( idmg(i) .eq. 0 )    idmg(i) = -1
               if( idmg(i) .eq. 1000 ) idmg(i) = 0

            end do

*-----------------------------------------------------------------------
*     card 16 - 19 : array
*        write array infomation on the temporary file 29
*-----------------------------------------------------------------------

            read(in,'(a200)') chin

                  call chlngt(chin,200,i1,i2)

                  call snum(chin,i1,i2,i3,cvvv,ierr)

                  if( ierr .ne. 0 ) then

                     write(io,'(''*** ERROR : in card 16'',
     &                          '' in old input'')')
                     goto 888

                  end if

               ic16 = nint( cvvv )

                  if( ic16 .ne. 0 .and. iuniv .eq. 0 ) then

                     write(io,'(''*** ERROR : universe is not'',
     &                          '' defined in old input'')')
                     goto 888

                  end if

         if( ic16 .ne. 0 ) then

               itar = 19
               open(itar,form='formatted',status='scratch')

  185          continue

               llarr = llarr + 1

               write(itar,'(i6,200a1)') i2, ( chin(i:i), i = 1, i2 )

               read(in,'(a200)') chin

                  call chlngt(chin,200,i1,i2)

                  chlw = chin
                  call chcaps(chlw,i1,i2,i3,'#!$')

               if(      chlw(i1:i1+9) .eq. 'importance' ) then

                  iimpo = 1
                  goto 180

               else if( chlw(i1:i1+9) .eq. 'no-importa' ) then

                  iimpo = 0
                  goto 180

               else if( chlw(i1:i1+5) .eq. 'source' ) then

                  iimpo = 0
                  backspace in
                  goto 180

               else

                  goto 185

               end if

         else

               read(in,'(a200)') chin

                  call chlngt(chin,200,i1,i2)

                  chlw = chin
                  call chcaps(chlw,i1,i2,i3,'#!$')

               if(      chlw(i1:i1+9) .eq. 'importance' ) then

                  iimpo = 1

               else if( chlw(i1:i1+9) .eq. 'no-importa' ) then

                  iimpo = 0

               else

                  iimpo = 0
                  backspace in

               end if

                  goto 180

         end if

*-----------------------------------------------------------------------

  180 continue

*-----------------------------------------------------------------------
*     card 20 : importance
*-----------------------------------------------------------------------

         if( iimpo .gt. 0 ) then

            read(in,*) ( dimp(i),i = 1, iregn )

            if( iuniv .eq. 0 ) then

               nrsq    = 5
               irsq(1) = 1
               irsq(2) = 2
               irsq(3) = 3
               irsq(4) = 4
               irsq(5) = 100

            else

               nrsq    = 6
               irsq(1) = 1
               irsq(2) = 2
               irsq(3) = 3
               irsq(4) = 4
               irsq(4) = 5
               irsq(5) = 100

            end if

         end if

*-----------------------------------------------------------------------
*     card 21 : source type
*-----------------------------------------------------------------------

         read(in,'(a200)') chin

            call chlngt(chin,200,i1,i2)

            chlw = chin
            call chcaps(chlw,i1,i2,i3,'#!$')

            if( chlw(i1:i1+5) .ne. 'source' .and.
     &          chlw(i1:i1+3) .ne. 'beam' ) then

               write(io,'('' Error : source identification label'',
     &                    '' was not found in old input'')')

               goto 888

            end if

               i4 = inumc(chlw,i1,i3,' ')
               i5 = jnumc(chlw,i4,i3)

               call snum(chlw,i5,i3,ic2,vvv,ierr)

               if( ierr .ne. 0 ) goto 997

               jstyp(1) = nint( vvv )

               if( jstyp(1) .lt. 1 .or.
     &           ( jstyp(1) .gt. 5 .and. jstyp(1) .lt. 11 ) .or.
     &             jstyp(1) .gt. 12 ) goto 997

*-----------------------------------------------------------------------
*     card 22 : source information
*-----------------------------------------------------------------------

         if( jstyp(1) .eq. 1 ) then

            read(in,*) sr0(1), sz0(1), sz1(1), se0(1), stip0,
     &                 stl, sdir(1)

         else if( jstyp(1) .eq. 2 ) then

            read(in,*) sx0(1), sx1(1), sy0(1), sy1(1), sz0(1), sz1(1),
     &                 se0(1), stip0, stl, sdir(1)

         else if( jstyp(1) .eq. 3 ) then

            read(in,*) sx0(1), sx1(1), sy0(1), sy1(1), sz0(1), sz1(1),
     &                 se0(1), stip0, stl, sdir(1)

         else if( jstyp(1) .eq. 4 ) then

            read(in,*) sr0(1), sz0(1), sz1(1), se0(1), stip0,
     &                 stl, sdir(1)

         else if( jstyp(1) .eq. 5 ) then

            read(in,*) sx0(1), sx1(1), sy0(1), sy1(1), sz0(1), sz1(1),
     &                 se0(1), stip0, stl, sdir(1)

         else if( jstyp(1) .eq. 11 ) then

            read(in,*) sx0(1), sx1(1), sy0(1), sy1(1), sz0(1), sz1(1),
     &                 srx(1), sry(1), swem(1), se0(1), stip0,
     &                 stl, sdir(1)

         else if( jstyp(1) .eq. 12 ) then

            read(in,*) sx0(1), sy0(1), sz0(1), stip0, stl, sdir(1)

            lsfile(1) = 12
            sfile(1)  = 'd-turtle.dat'

         end if

            istyp(1) = nint( stip0 ) + 1
            inkf0(1) = kfft( istyp(1) )
            lstyp(1) = -10000

            sphi(1) = -10000.0
            sdom(1) = -10000.0
            swt0(1) = 1.0d0

*-----------------------------------------------------------------------
*     card 23 : source energy group
*-----------------------------------------------------------------------

         if( jstyp(1) .eq. 4 .or. jstyp(1) .eq. 5 ) then

            read(in,*) ngrp(1)

               call moddas_reallocate_dbl(isrc, 1, ngrp(1), ngei, egmin)
               call moddas_reallocate_dbl(isrc, 1, ngrp(1), ngea, egmax)
               call moddas_reallocate_dbl(isrc, 1, ngrp(1), ngfe, fegrp)
               call moddas_reallocate_dbl(isrc, 1, ngrp(1), ngft, rfe)

            read(in,*) ( egmin(ngei(1)+i),
     &                   fegrp(ngfe(1)+i), i = 1, ngrp(1) ),
     &                   egmax(ngea(1)+ngrp(1))

            do i = 1, ngrp(1) - 1

               egmax(ngea(1)+i) = egmin(ngei(1)+i+1)

            end do

               jetyp(1) = 1

         end if

*-----------------------------------------------------------------------
*        open temporary file for tally
*-----------------------------------------------------------------------

 9000       continue

            isc = 15
            open(isc,status='scratch',form='formatted')

*-----------------------------------------------------------------------
*     card 24 : read tally options
*-----------------------------------------------------------------------

  881    read(in,'(a)', iostat = ios ) hl
         if( ios .eq. -1 ) goto 520

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 881
            call chnumb(hl,10,icl1,icl2)
            if( icl1 .eq. 0 ) goto 881

*-----------------------------------------------------------------------

         if( hl(icl1:icl2) .eq. 'tally*' ) then

            read(hl(11:80),*) italxx

*-----------------------------------------------------------------------

         else if( hl(icl1:icl2) .eq. 'all tally' .or.
     &            hl(icl1:icl2) .eq. 'all-tally' ) then

            do i = 1, 8
               italxx(i) = 1
            end  do

*-----------------------------------------------------------------------

         else if( hl(icl1:icl2) .eq. 'tally' ) then


  110          call chnumb(hl(11:80),70,icl3,icl4)

               icl3 = icl3 + 10
               icl4 = icl4 + 10

            if( hl(icl3:icl4) .eq. 'track' ) then

               italxx(1) = 1

            else if( hl(icl3:icl4) .eq. 'surface' ) then

               italxx(2) = 1
               italxx(3) = 1
               italxx(4) = 1

            else if( hl(icl3:icl4) .eq. 'surface-spectrum' .or.
     &               hl(icl3:icl4) .eq. 'surface-spect' ) then

               italxx(2) =1

            else if( hl(icl3:icl4) .eq. 'surface-flux' ) then

               italxx(3) = 1

            else if( hl(icl3:icl4) .eq. 'surface-current' ) then

               italxx(4) = 1

            else if( hl(icl3:icl4) .eq. 'yield' ) then

               italxx(5) = 1

            else if( hl(icl3:icl4) .eq. 'yield-special' ) then

               italxx(5) = 2

            else if( hl(icl3:icl4) .eq. 'yield-select' ) then

               italxx(5) = 3

            else if( hl(icl3:icl4) .eq. 'yield-select-special' .or.
     &               hl(icl3:icl4) .eq. 'yield-special-select' ) then

               italxx(5) = 4

            else if( hl(icl3:icl4) .eq.'heat' .or.
     &               hl(icl3:icl4) .eq.'heating' ) then

               italxx(6) = 1

            else if( hl(icl3:icl4) .eq. 'bndsurf-current' ) then

               italxx(7) = 1

            else if( hl(icl3:icl4) .eq. 'bndsurf-flux' ) then

               italxx(8) = 1

            else if( hl(icl3:icl4) .eq. 'bndsurf' ) then

               italxx(7) = 1
               italxx(8) = 1

            end if

*-----------------------------------------------------------------------

  882       read(in,'(a)', iostat = ios ) hl
            if( ios .eq. -1 ) goto 804

               if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 882
               call chnumb(hl,10,icl3,icl4)

*-----------------------------------------------------------------------

            if( icl3 .eq. 0 .or.
     &          hl(icl3:icl4) .ne. 'tally' ) then

               goto 201

            else

               goto 110

            end if

*-----------------------------------------------------------------------

         else

            goto 510

         end if

*-----------------------------------------------------------------------
*     card 25 - 36 : read tally options in region mesh
*-----------------------------------------------------------------------

  200    read(in,'(a)',end=500) hl

  201    continue

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 200
            call chnumb(hl,10,icl1,icl2)
            if( icl1 .eq. 0 ) goto 200

*-----------------------------------------------------------------------
*     card 25 : energy bin
*-----------------------------------------------------------------------

         if( hl(icl1:icl2) .eq. 'energy' ) then

            read(hl(11:80),*) neng

            icc = 0

  210       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 221

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 210
            call chnumb(hl,10,icl3,icl4)
            if( hl(icl3:icl4) .ne. 'energy' ) goto 220
            write(isc,'(a)') hl(11:80)
            goto 210

  221       icc = 1

  220       rewind isc

            read(isc,*,iostat=ios) ( steg(i), i = 1, neng + 1 )
            if( ios .eq. -1 ) goto 225
  225       rewind isc

            do i = 1, neng
               if( steg(i+1) .le. steg(i) ) then

                  write(io,'(''energy bin must be the '',
     &                       ''ascending order.''/
     &                       '' bin='',i4,''   energy='',1p2e13.4)')
     &            i, steg(i+1), steg(i)

                  goto 888

               end if
            end do

            iengm = 1

            if( icc .eq. 1 ) goto 200

*-----------------------------------------------------------------------
*     card 26 : region volume
*-----------------------------------------------------------------------

         else if( hl(icl1:icl2) .eq. 'volume' ) then

            write(isc,'(a)') hl(11:80)

            icc = 0

  240       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 251

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 240
            call chnumb(hl,10,icl3,icl4)
            if(hl(icl3:icl4).ne.'volume') goto 250
            write(isc,'(a)') hl(11:80)
            goto 240

  251       icc = 1

  250       rewind isc

            do i = 1, iregn
               rvols(i) = -1.0
            end do

            read(isc,*,iostat=ios) (rvols(i),i=1,iregn)
            if( ios .eq. -1 ) goto 255
  255       rewind isc

            do i = 1, iregn

               if( rvols(i) .gt. 0.0 ) then
                  mnvol = mnvol + 1
                  nvols(i) = i
               end if

            end do

            if( icc .eq. 1 ) goto 200

*-----------------------------------------------------------------------
*     card 27 : r crossing surface
*-----------------------------------------------------------------------

         else if(hl(icl1:icl2) .eq. 'r-surf' ) then

            read(hl(11:80),*) nrsurf

            icc = 0

  290       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 299

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 290
            call chnumb(hl,10,icl3,icl4)
            if(hl(icl3:icl4).ne.'r-surf') goto 300
            write(isc,'(a)') hl(11:80)
            goto 290

  299       icc = 1

  300       rewind isc

            read(isc,*) (strg(i),i=1,nrsurf)
            rewind isc

            do i = 1, nrsurf - 1
               if( strg(i+1) .le. strg(i) ) then

                  write(io,'(''r-surf bin must be the '',
     &                       ''ascending order.''/
     &                       '' bin='',i4,''   r-surf='',1p2e13.4)')
     &            i, strg(i+1), strg(i)

                  goto 888

               end if
            end do

            if( icc .eq. 1 ) goto 200

*-----------------------------------------------------------------------
*     card 28 : z crossing surface
*-----------------------------------------------------------------------

         else if( hl(icl1:icl2) .eq. 'z-surf' ) then

            read(hl(11:80),*) nzsurf

            icc = 0

  340       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 349

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 340
            call chnumb(hl,10,icl3,icl4)
            if(hl(icl3:icl4).ne.'z-surf') goto 350
            write(isc,'(a)') hl(11:80)
            goto 340

  349       icc = 1

  350       rewind isc

            read(isc,*) (stzg(i),i=1,nzsurf)
            rewind isc

            do i = 1, nzsurf - 1
               if( stzg(i+1) .le. stzg(i) ) then

                  write(io,'(''z-surf bin must be the '',
     &                       ''ascending order.''/
     &                       '' bin='',i4,''   z-surf='',1p2e13.4)')
     &            i, stzg(i+1), stzg(i)

                  goto 888

               end if
            end do

            if( icc .eq. 1 ) goto 200

*-----------------------------------------------------------------------
*     card 29 : region mesh
*-----------------------------------------------------------------------

         else if( hl(icl1:icl2) .eq. 'region' ) then

            read(hl(11:80),*) nrreg

            nd = nrreg
            if( nrreg .lt. 0 ) nd = 2

            icc = 0

  370       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 378

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 370
            call chnumb(hl,10,icl3,icl4)
            if(hl(icl3:icl4).ne.'region') goto 380
            write(isc,'(a)') hl(11:80)
            goto 370

  378       icc = 1

  380       rewind isc

            read(isc,*) (ntrg(i),i=1,nd)
            rewind isc

            if(nrreg.lt.0) then
               i01=ntrg(1)
               i02=ntrg(2)
               nrreg=i02-i01+1
               do i=2,nrreg
                  ntrg(i)=i01+i-1
               end do
            end if

            irreg = 1

            if( icc .eq. 1 ) goto 200

*-----------------------------------------------------------------------
*     card 30 : yield region mesh
*-----------------------------------------------------------------------

         else if( hl(icl1:icl2) .eq. 'y-region' ) then

            read(hl(11:80),*) nyreg

            nd = nyreg
            if( nyreg .lt. 0 ) nd = 2

            icc = 0

  371       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 379

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 371
            call chnumb(hl,10,icl3,icl4)
            if(hl(icl3:icl4).ne.'y-region') goto 381
            write(isc,'(a)') hl(11:80)
            goto 371

  379       icc = 1

  381       rewind isc

            read(isc,*) (ntyg(i),i=1,nd)
            rewind isc

            if(nyreg.lt.0) then
               i01=ntyg(1)
               i02=ntyg(2)
               nyreg=i02-i01+1
               do i=2,nyreg
                  ntyg(i)=i01+i-1
               end do
            end if

            iyreg = 1

            if( icc .eq. 1 ) goto 200

*-----------------------------------------------------------------------
*     card 31 : heat region mesh
*-----------------------------------------------------------------------

         else if( hl(icl1:icl2) .eq. 'h-region' ) then

            read(hl(11:80),*) nhreg

            nd = nhreg
            if( nhreg .lt. 0 ) nd = 2

            icc = 0

  372       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 376

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 372
            call chnumb(hl,10,icl3,icl4)
            if(hl(icl3:icl4).ne.'h-region') goto 382
            write(isc,'(a)') hl(11:80)
            goto 372

  376       icc = 1

  382       rewind isc

            read(isc,*) (nthg(i),i=1,nd)
            rewind isc

            if(nhreg.lt.0) then
               i01=nthg(1)
               i02=nthg(2)
               nhreg=i02-i01+1
               do i=2,nhreg
                  nthg(i)=i01+i-1
               end do
            end if

            ihreg = 1

            if( icc .eq. 1 ) goto 200

*-----------------------------------------------------------------------
*     card 32 : track rength region mesh
*-----------------------------------------------------------------------

         else if( hl(icl1:icl2) .eq. 't-region' ) then

            read(hl(11:80),*) ntreg

            nd = ntreg
            if( ntreg .lt. 0 ) nd = 2

            icc = 0

  373       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 384

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 373
            call chnumb(hl,10,icl3,icl4)
            if(hl(icl3:icl4).ne.'t-region') goto 383
            write(isc,'(a)') hl(11:80)
            goto 373

  384       icc = 1

  383       rewind isc

            read(isc,*) (nttg(i),i=1,nd)
            rewind isc

            if(ntreg.lt.0) then
               i01=nttg(1)
               i02=nttg(2)
               ntreg=i02-i01+1
               do i=2,ntreg
                  nttg(i)=i01+i-1
               end do
            end if

            ikreg = 1

            if( icc .eq. 1 ) goto 200

*-----------------------------------------------------------------------
*     card 33 : sp-unit
*-----------------------------------------------------------------------

         else if( hl(icl1:icl2) .eq. 'sp-unit' ) then

            read(hl(11:80),*) ksunit

            if( ksunit .lt. 0 .or. ksunit .gt. 2 ) then

                  write(io,'(''ksunit should be 0, 1, 2 in old input''
     &                     )')

                  goto 888

            end if

            goto 200

*-----------------------------------------------------------------------
*     card 34 : pt-surf
*-----------------------------------------------------------------------

         else if(hl(icl1:icl2).eq.'pt-surf') then

            read(hl(11:80),*) kpsurf

            if( kpsurf .ne. 2 .and.
     &          kpsurf .ne. 3 .and.
     &          kpsurf .ne. 5 ) kpsurf = 3

            goto 200

*-----------------------------------------------------------------------
*     card 35 : crossiong region
*-----------------------------------------------------------------------

         else if( hl(icl1:icl2) .eq. 'b-surf' ) then

            read(hl(11:80),*) nbsurf

            icc = 0

  390       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 389

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 390
            call chnumb(hl,10,icl3,icl4)
            if(hl(icl3:icl4).ne.'b-surf') goto 391
            write(isc,'(a)') hl(11:80)
            goto 390

  389       icc = 1

  391       rewind isc

            read(isc,*) (ntcr(i,1),ntcr(i,2),i=1,nbsurf)
            rewind isc

            if( icc .eq. 1 ) goto 200

*-----------------------------------------------------------------------
*     card 36 : area for crossing region
*-----------------------------------------------------------------------

         else if( hl(icl1:icl2) .eq. 'area' ) then

            write(isc,'(a)') hl(11:80)

            icc = 0

  392       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 394

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 392
            call chnumb(hl,10,icl3,icl4)
            if(hl(icl3:icl4).ne.'area') goto 393
            write(isc,'(a)') hl(11:80)
            goto 392

  394       icc = 1

  393       rewind isc

            do i = 1, nbsurf
               stca(i) = -1.0
            end do

            read(isc,*,iostat=ios) (stca(i),i=1,nbsurf)
            if( ios .eq. -1 ) goto 455
  455       rewind isc

            do i = 1, nbsurf

               if( stca(i) .lt. 0.0 ) then
                  stca(i) = 1.0
               end if

            end do

            if( icc .eq. 1 ) goto 200

*-----------------------------------------------------------------------
*     card ** : yield special
*-----------------------------------------------------------------------

         else if(hl(icl1:icl2).eq.'n-yield') then

            read(hl(11:80),*) kyield

            goto 200

*-----------------------------------------------------------------------
*     card ** : yield select mother
*-----------------------------------------------------------------------

         else if(hl(icl1:icl2).eq.'s-yield') then

            read(hl(11:80),*) lyield

               if( lyield .le. 0 ) then

                  write(io,'(''  **** s-yield is wrong,'',
     &                       '' should be greater than 0'')')
                  goto 888

               end if

               iclyd = 0

               icc = 0

  633          read(in,'(a)',iostat=ios) hl
               if( ios .eq. -1 ) goto 636

               if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 633
               call chnumb(hl,10,icl3,icl4)
               if( hl(icl3:icl4) .ne. 's-yield' ) goto 634

               iclyd = iclyd + 1

               write(isc,'(a)') hl(11:80)
               goto 633

  636          icc = 1

  634          rewind isc

               if( iclyd .ne. lyield + 1 ) then

                  write(io,'(
     &                '' **** description of s-yield is wrong'')')
                  goto 888

               end if

               read(isc,*) ( lydn(i), i = 1, lyield )

                  lyzi = 0

               do i = 1, lyield

                  lyzi = lyzi + lydn(i-1)

                  read(isc,*,iostat=ios)
     &            ( lydz(lyzi+j), lyda(lyzi+j), j = 1, lydn(i) )
                  if( ios .eq. -1 ) goto 635

               end do

  635          rewind isc

               if( icc .eq. 1 ) goto 200

*-----------------------------------------------------------------------

         else

            goto 500

         end if

            goto 201

*-----------------------------------------------------------------------
*     summary and check of region mesh tally
*-----------------------------------------------------------------------

  500 continue

*-----------------------------------------------------------------------
*        t-track of region mesh
*-----------------------------------------------------------------------

            if( italxx(1) .ne. 0 .and.
     &          ikreg .eq. 0 .and. irreg .eq. 0 ) then

                  write(io,'(
     &                ''t-region is not defined'')')
                  goto 888

            end if

            if( italxx(1) .ne. 0 .and.
     &          iengm .eq. 0 ) then

                  write(io,'(
     &                ''energy mesh is not defined'')')
                  goto 888

            end if

            if( italxx(1) .ne. 0 .and.
     &          ikreg .eq. 0 .and. irreg .ne. 0 ) then

                  ntreg = nrreg

               do i = 1, nrreg

                  nttg(i) = ntrg(i)

               end do

            end if

*-----------------------------------------------------------------------

         if( italxx(1) .ne. 0 ) then

                  itnm = itnm + 1

                  ital( itnm ) = 1

                  itsmn(itnm) = 0
                  itstm(itnm) = 0

                  ittll(itnm) = 24
                  ittle(itnm) = '[t-track] in region mesh'

                  itanl(itnm) = 0
                  itaxl(itnm) = 0
                  itayl(itnm) = 0
                  itazl(itnm) = 0
                  itgsh(itnm) = 0
                  itrsh(itnm) = 0
                  itmtn(itnm) = 0
                  itmcn(itnm) = 0
                  iterl(itnm) = 72
                  rtfac(itnm) = 1.0
                  itvm(itnm)  = 9
                  iteps(itnm) = 0

                  itmsh(itnm) = 1

               if( ntreg .eq. 0 ) then

                  ntreg = 1
                  nttg(1) = 50000

               end if

                  itrgn(itnm) = ntreg
                  itrnv(itnm) = 0

                  itrgm(itnm) = ntreg
                  call moddas_reallocate_int(
     &                    itlmax, itnm, ntreg, itreg, idas_itreg)

                  itsmn(itnm) = itsmn(itnm)
     &                        + ( ntreg + mod(ntreg,2) ) / 2


               do i = 1, ntreg

                  idas_itreg( itreg(itnm) + i - 1 ) = nttg(i)

               end do

               if( kpsurf .eq. 2 ) then

                  itpan(itnm) = 2

               else if( kpsurf .eq. 3 ) then

                  itpan(itnm) = 4

               else if( kpsurf .eq. 5 ) then

                  itpan(itnm) = 6

               else

                  itpan(itnm) = 2

               end if

               do i = 1, 6

                  itpat(itnm,i,1) = iptyp(i)
                  itpat(itnm,i,2) = ipnkf(i)

               end do

                  itety(itnm) = 1
                  itenm(itnm) = neng
                  rtedl(itnm) = 0.0
                  rtemi(itnm) = steg(neng+1)
                  rtema(itnm) = steg(1)

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, neng+1, iterg, das_iterg)

                  itsmn(itnm) = itsmn(itnm) + neng + 1


               do i = 1, neng + 1

                  das_iterg( iterg(itnm) + i - 1 ) = steg(i)

               end do

                  itunt(itnm) = ksunit + 1

                  ittwo(itnm) = 3

                  itaxn(itnm) = 1
                  itaxs(itnm,1) = 1

                  itfln(itnm) = 1
                  itfll(itnm,1) = 7
                  ctfln(itnm,1) = 'fort.32'

         end if

*-----------------------------------------------------------------------
*        t-cross  ( currennt spectrum )
*-----------------------------------------------------------------------

            if( italxx(2) .ne. 0 .and.
     &          iengm .eq. 0 ) then

                  write(io,'(
     &                ''energy mesh is not defined'')')
                  goto 888

            end if

            if( italxx(2) .ne. 0 .and.
     &        ( nrsurf .le. 1 .or. nzsurf .le. 1 ) ) then

                  write(io,'(
     &                ''r-z mesh definitions are wrong or missing'')')
                  goto 888

            end if

*-----------------------------------------------------------------------

         if( italxx(2) .ne. 0 ) then

                  itnm = itnm + 1

                  ital( itnm ) = 2

                  itsmn(itnm) = 0
                  itstm(itnm) = 0

                  ittll(itnm) = 21
                  ittle(itnm) = '[t-cross] in r-z mesh'

                  itanl(itnm) = 0
                  itaxl(itnm) = 0
                  itayl(itnm) = 0
                  itazl(itnm) = 0
                  itgsh(itnm) = 0
                  itrsh(itnm) = 0
                  itmtn(itnm) = 0
                  itmcn(itnm) = 0
                  iterl(itnm) = 72
                  rtfac(itnm) = 1.0
                  itvm(itnm)  = 9
                  iteps(itnm) = 0

                  itmsh(itnm) = 2

                  rtrx0(itnm) = 0.0
                  rtry0(itnm) = 0.0

                  inr = nrsurf - 1

                  itrty(itnm) = 1
                  itrnm(itnm) = inr
                  rtrdl(itnm) = 0.0
                  rtrmi(itnm) = strg(1)
                  rtrma(itnm) = strg(inr+1)

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, inr+1, itrrg, das_itrrg)

                  itsmn(itnm) = itsmn(itnm) + inr + 1


               do i = 1, inr + 1

                  das_itrrg( itrrg(itnm) + i - 1 ) = strg(i)

               end do

                  inz = nzsurf - 1

                  itzty(itnm) = 1
                  itznm(itnm) = inz
                  rtzdl(itnm) = 0.0
                  rtzmi(itnm) = stzg(1)
                  rtzma(itnm) = stzg(inz+1)

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, inz+1, itzrg, das_itzrg)

                  itsmn(itnm) = itsmn(itnm) + inz + 1


               do i = 1, inz + 1

                  das_itzrg( itzrg(itnm) + i - 1 ) = stzg(i)

               end do

               if( kpsurf .eq. 2 ) then

                  itpan(itnm) = 2

               else if( kpsurf .eq. 3 ) then

                  itpan(itnm) = 4

               else if( kpsurf .eq. 5 ) then

                  itpan(itnm) = 6

               else

                  itpan(itnm) = 2

               end if

               do i = 1, 6

                  itpat(itnm,i,1) = iptyp(i)
                  itpat(itnm,i,2) = ipnkf(i)

               end do

                  itety(itnm) = 1
                  itenm(itnm) = neng
                  rtedl(itnm) = 0.0
                  rtemi(itnm) = steg(1)
                  rtema(itnm) = steg(neng+1)

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, neng+1, iterg, das_iterg)

                  itsmn(itnm) = itsmn(itnm) + neng + 1


               do i = 1, neng + 1

                  das_iterg( iterg(itnm) + i - 1 ) = steg(i)

               end do

                  itunt(itnm) = ksunit + 1

                  ittwo(itnm) = 3

                  itaxn(itnm) = 1
                  itaxs(itnm,1) = 1

                  itfln(itnm) = 1
                  itfll(itnm,1) = 7
                  ctfln(itnm,1) = 'fort.36'

                  itout(itnm) = 2

         end if

*-----------------------------------------------------------------------
*        t-cross  ( flux spectrum )
*-----------------------------------------------------------------------

            if( italxx(3) .ne. 0 .and.
     &          iengm .eq. 0 ) then

                  write(io,'(
     &                ''energy mesh is not defined'')')
                  goto 888

            end if

            if( italxx(3) .ne. 0 .and.
     &        ( nrsurf .le. 1 .or. nzsurf .le. 1 ) ) then

                  write(io,'(
     &                ''r-z mesh definitions are wrong or missing'')')
                  goto 888

            end if

*-----------------------------------------------------------------------

         if( italxx(3) .ne. 0 ) then

                  itnm = itnm + 1

                  ital( itnm ) = 2

                  itsmn(itnm) = 0
                  itstm(itnm) = 0

                  ittll(itnm) = 21
                  ittle(itnm) = '[t-cross] in r-z mesh'

                  itanl(itnm) = 0
                  itaxl(itnm) = 0
                  itayl(itnm) = 0
                  itazl(itnm) = 0
                  itgsh(itnm) = 0
                  itrsh(itnm) = 0
                  itmtn(itnm) = 0
                  itmcn(itnm) = 0
                  iterl(itnm) = 72
                  rtfac(itnm) = 1.0
                  itvm(itnm)  = 9
                  iteps(itnm) = 0

                  itmsh(itnm) = 2

                  rtrx0(itnm) = 0.0
                  rtry0(itnm) = 0.0

                  inr = nrsurf - 1

                  itrty(itnm) = 1
                  itrnm(itnm) = inr
                  rtrdl(itnm) = 0.0
                  rtrmi(itnm) = strg(1)
                  rtrma(itnm) = strg(inr+1)

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, inr+1, itrrg, das_itrrg)

                  itsmn(itnm) = itsmn(itnm) + inr + 1


               do i = 1, inr + 1

                  das_itrrg( itrrg(itnm) + i - 1 ) = strg(i)

               end do

                  inz = nzsurf - 1

                  itzty(itnm) = 1
                  itznm(itnm) = inz
                  rtzdl(itnm) = 0.0
                  rtzmi(itnm) = stzg(1)
                  rtzma(itnm) = stzg(inz+1)

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, inz+1, itzrg, das_itzrg)

                  itsmn(itnm) = itsmn(itnm) + inz + 1


               do i = 1, inz + 1

                  das_itzrg( itzrg(itnm) + i - 1 ) = stzg(i)

               end do

               if( kpsurf .eq. 2 ) then

                  itpan(itnm) = 2

               else if( kpsurf .eq. 3 ) then

                  itpan(itnm) = 4

               else if( kpsurf .eq. 5 ) then

                  itpan(itnm) = 6

               else

                  itpan(itnm) = 2

               end if

               do i = 1, 6

                  itpat(itnm,i,1) = iptyp(i)
                  itpat(itnm,i,2) = ipnkf(i)

               end do

                  itety(itnm) = 1
                  itenm(itnm) = neng
                  rtedl(itnm) = 0.0
                  rtemi(itnm) = steg(1)
                  rtema(itnm) = steg(neng+1)

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, neng+1, iterg, das_iterg)

                  itsmn(itnm) = itsmn(itnm) + neng + 1


               do i = 1, neng + 1

                  das_iterg( iterg(itnm) + i - 1 ) = steg(i)

               end do

                  itunt(itnm) = ksunit + 1

                  ittwo(itnm) = 3

                  itaxn(itnm) = 1
                  itaxs(itnm,1) = 1

                  itfln(itnm) = 1
                  itfll(itnm,1) = 7
                  ctfln(itnm,1) = 'fort.37'

                  itout(itnm) = 1

         end if

*-----------------------------------------------------------------------
*        t-cross  ( omni current )
*-----------------------------------------------------------------------

            if( italxx(4) .ne. 0 .and.
     &        ( nrsurf .le. 1 .or. nzsurf .le. 1 ) ) then

                  write(io,'(
     &                ''r-z mesh definitions are wrong or missing'')')
                  goto 888

            end if

*-----------------------------------------------------------------------

         if( italxx(4) .ne. 0 ) then

                  itnm = itnm + 1

                  ital( itnm ) = 2

                  itsmn(itnm) = 0
                  itstm(itnm) = 0

                  ittll(itnm) = 21
                  ittle(itnm) = '[t-cross] in r-z mesh'

                  itanl(itnm) = 0
                  itaxl(itnm) = 0
                  itayl(itnm) = 0
                  itazl(itnm) = 0
                  itgsh(itnm) = 0
                  itrsh(itnm) = 0
                  itmtn(itnm) = 0
                  itmcn(itnm) = 0
                  iterl(itnm) = 72
                  rtfac(itnm) = 1.0
                  itvm(itnm)  = 9
                  iteps(itnm) = 0

                  itmsh(itnm) = 2

                  rtrx0(itnm) = 0.0
                  rtry0(itnm) = 0.0

                  inr = nrsurf - 1

                  itrty(itnm) = 1
                  itrnm(itnm) = inr
                  rtrdl(itnm) = 0.0
                  rtrmi(itnm) = strg(1)
                  rtrma(itnm) = strg(inr+1)

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, inr+1, itrrg, das_itrrg)

                  itsmn(itnm) = itsmn(itnm) + inr + 1


               do i = 1, inr + 1

                  das_itrrg( itrrg(itnm) + i - 1 ) = strg(i)

               end do

                  inz = nzsurf - 1

                  itzty(itnm) = 1
                  itznm(itnm) = inz
                  rtzdl(itnm) = 0.0
                  rtzmi(itnm) = stzg(1)
                  rtzma(itnm) = stzg(inz+1)

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, inz+1, itzrg, das_itzrg)

                  itsmn(itnm) = itsmn(itnm) + inz + 1


               do i = 1, inz + 1

                  das_itzrg( itzrg(itnm) + i - 1 ) = stzg(i)

               end do

               if( kpsurf .eq. 2 ) then

                  itpan(itnm) = 2

               else if( kpsurf .eq. 3 ) then

                  itpan(itnm) = 4

               else if( kpsurf .eq. 5 ) then

                  itpan(itnm) = 6

               else

                  itpan(itnm) = 2

               end if

               do i = 1, 6

                  itpat(itnm,i,1) = iptyp(i)
                  itpat(itnm,i,2) = ipnkf(i)

               end do

                  itety(itnm) = 1
                  itenm(itnm) = 1
                  rtedl(itnm) = 0.0
                  rtemi(itnm) = 0.0
                  rtema(itnm) = 1.0d+10

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, neng+1, iterg, das_iterg)

                  itsmn(itnm) = itsmn(itnm) + neng + 1


                  das_iterg( iterg(itnm)     ) = rtemi(itnm)
                  das_iterg( iterg(itnm) + 1 ) = rtema(itnm)

                  itunt(itnm) = 1
                  ittwo(itnm) = 3

                  itaxn(itnm) = 1
                  itaxs(itnm,1) = 5

                  itfln(itnm) = 1
                  itfll(itnm,1) = 7
                  ctfln(itnm,1) = 'fort.38'

                  itout(itnm) = 5

         end if

*-----------------------------------------------------------------------
*        t-yield
*-----------------------------------------------------------------------

            if( italxx(5) .ne. 0 .and. lyield .gt. 9 ) then

                  write(io,'(
     &                ''mother group should be less than 10'')')
                  goto 888

            end if

            if( italxx(5) .ne. 0 .and.
     &          iyreg .eq. 0 .and. irreg .eq. 0 ) then

                  write(io,'(
     &                ''y-region is not defined'')')
                  goto 888

            end if

            if( italxx(5) .ne. 0 .and.
     &          iyreg .eq. 0 .and. irreg .ne. 0 ) then

                  nyreg = nrreg

               do i = 1, nrreg

                  ntyg(i) = ntrg(i)

               end do

            end if

*-----------------------------------------------------------------------

         if( italxx(5) .ne. 0 ) then

               if( lyield .eq. 0 ) then

                  iydo  = 1

               else

                  iydo  = lyield
                  lyzi  = 0

               end if

*-----------------------------------------------------------------------

            do l = 1, iydo

                  itnm = itnm + 1

                  ital( itnm ) = 3

                  itsmn(itnm) = 0
                  itstm(itnm) = 0

                  ittll(itnm) = 24
                  ittle(itnm) = '[t-yield] in region mesh'

                  itanl(itnm) = 0
                  itaxl(itnm) = 0
                  itayl(itnm) = 0
                  itazl(itnm) = 0
                  itgsh(itnm) = 0
                  itrsh(itnm) = 0
                  itmtn(itnm) = 0
                  itmcn(itnm) = 0
                  iterl(itnm) = 72
                  rtfac(itnm) = 1.0
                  itvm(itnm)  = 9
                  iteps(itnm) = 0

                  itspc(itnm) = kyield

                  imoth = lydn(l)

                  itman(itnm) = imoth

*-----------------------------------------------------------------------

               if( imoth .gt. 0 ) then

                     call moddas_reallocate_int(
     &                       itlmax, itnm, imoth, itmat, ismat)

                     itsmn(itnm) = itsmn(itnm)
     &                           + ( imoth + mod(imoth,2) ) / 2


                     lyzi = lyzi + lydn(l-1)

                  do i = 1, imoth

                     izma = lydz(lyzi+i) * 1000 + lyda(lyzi+i)

                     ismat( itmat(itnm) + i - 1 ) = izma

                  end do

               end if

*-----------------------------------------------------------------------

                  itnun(itnm) = 0

                  itmsh(itnm) = 1

               if( nyreg .eq. 0 ) then

                  nyreg = 1
                  ntyg(1) = 50000

               end if

                  itrgn(itnm) = nyreg
                  itrnv(itnm) = 0

                  itrgm(itnm) = nyreg

                  call moddas_reallocate_int(
     &                    itlmax, itnm, nyreg, itreg, idas_itreg)
                  itsmn(itnm) = itsmn(itnm)
     &                        + ( nyreg + mod(nyreg,2) ) / 2


               do i = 1, nyreg

                  idas_itreg( itreg(itnm) + i - 1 ) = ntyg(i)

               end do

                  itaxn(itnm)   = 1
                  itaxs(itnm,1) = 13

                  itfln(itnm) = 1
                  itfll(itnm,1) = 7

                  write(dam,'(i1)') l-1
                  ctfln(itnm,1) = 'fort.4'//dam

                  itunt(itnm) = 1
                  itout(itnm) = 0
                  ittwo(itnm) = 3

                  itnfn(itnm) = 1
                  itnfr(itnm) = 1
                  itnfz(itnm) = 1

                  itndz(itnm) = maxpt
                  itndn(itnm) = maxnt

               do jz = 1, maxpt
               do jn = 1, maxnt

                  ikzz(jz,jn) = jz
                  iknn(jz,jn) = jn

               end do
               end do

                  maxik = maxpt * maxnt

                  itnkz(itnm) = ( mmmax - 1 ) * 2 + 1
                  mmmax = mmmax + ( maxik + mod(maxik,2) ) / 2

                  itsmn(itnm) = itsmn(itnm)
     &                        + ( maxik + mod(maxik,2) ) / 2

                  if( mmmax .gt. mdas ) goto 950

                  k = 0

               do j = 1, maxnt
               do i = 1, maxpt

                  k = k + 1

                  idas( itnkz(itnm) + k - 1 ) = ikzz(i,j)

               end do
               end do

                  itnkn(itnm) = ( mmmax - 1 ) * 2 + 1
                  mmmax = mmmax + ( maxik + mod(maxik,2) ) / 2

                  itsmn(itnm) = itsmn(itnm)
     &                        + ( maxik + mod(maxik,2) ) / 2

                  if( mmmax .gt. mdas ) goto 950

                  k = 0

               do j = 1, maxnt
               do i = 1, maxpt

                  k = k + 1

                  idas( itnkn(itnm) + k - 1 ) = iknn(i,j)

               end do
               end do

*-----------------------------------------------------------------------

            end do

         end if

*-----------------------------------------------------------------------
*        t-heat
*-----------------------------------------------------------------------

            if( italxx(6) .ne. 0 .and.
     &          ihreg .eq. 0 .and. irreg .eq. 0 ) then

                  write(io,'(
     &                ''h-region is not defined'')')
                  goto 888

            end if

            if( italxx(6) .ne. 0 .and.
     &          ihreg .eq. 0 .and. irreg .ne. 0 ) then

                  nhreg = nrreg

               do i = 1, nrreg

                  nthg(i) = ntrg(i)

               end do

            end if

*-----------------------------------------------------------------------

         if( italxx(6) .ne. 0 ) then

                  itnm = itnm + 1

                  ital( itnm ) = 4

                  itsmn(itnm) = 0
                  itstm(itnm) = 0

                  ittll(itnm) = 23
                  ittle(itnm) = '[t-heat] in region mesh'

                  itanl(itnm) = 0
                  itaxl(itnm) = 0
                  itayl(itnm) = 0
                  itazl(itnm) = 0
                  itgsh(itnm) = 0
                  itrsh(itnm) = 0
                  itmtn(itnm) = 0
                  itmcn(itnm) = 0
                  iterl(itnm) = 72
                  rtfac(itnm) = 1.0
                  itvm(itnm)  = 9
                  iteps(itnm) = 0

                  itmsh(itnm) = 1

               if( nhreg .eq. 0 ) then

                  nhreg = 1
                  nthg(1) = 50000

               end if

                  itrgn(itnm) = nhreg
                  itrnv(itnm) = 0

                  itrgm(itnm) = nhreg
                  call moddas_reallocate_int(
     &                    itlmax, itnm, nhreg, itreg, idas_itreg)

                  itsmn(itnm) = itsmn(itnm)
     &                        + ( nhreg + mod(nhreg,2) ) / 2


               do i = 1, nhreg

                  idas_itreg( itreg(itnm) + i - 1 ) = nthg(i)

               end do

                  itpan(itnm) = 0
                  itunt(itnm) = 1
                  ittwo(itnm) = 3

                  itaxn(itnm)   = 1
                  itaxs(itnm,1) = 2

                  itfln(itnm) = 1
                  itfll(itnm,1) = 7
                  ctfln(itnm,1) = 'fort.56'

                  itout(itnm) = 1
                  itndy(itnm) = 28

         end if

*-----------------------------------------------------------------------
*        t-cross  ( currennt spectrum by boundary crossing)
*-----------------------------------------------------------------------

            if( italxx(7) .ne. 0 .and.
     &          iengm .eq. 0 ) then

                  write(io,'(
     &                ''energy mesh is not defined'')')
                  goto 888

            end if

            if( italxx(7) .ne. 0 .and. nbsurf .eq. 1 ) then

                  write(io,'(
     &                ''b-surf is not defined'')')
                  goto 888

            end if

*-----------------------------------------------------------------------

         if( italxx(7) .ne. 0 ) then

                  itnm = itnm + 1

                  ital( itnm ) = 2

                  itsmn(itnm) = 0
                  itstm(itnm) = 0

                  ittll(itnm) = 24
                  ittle(itnm) = '[t-cross] in region mesh'

                  itanl(itnm) = 0
                  itaxl(itnm) = 0
                  itayl(itnm) = 0
                  itazl(itnm) = 0
                  itgsh(itnm) = 0
                  itrsh(itnm) = 0
                  itmtn(itnm) = 0
                  itmcn(itnm) = 0
                  iterl(itnm) = 72
                  rtfac(itnm) = 1.0
                  itvm(itnm)  = 9
                  iteps(itnm) = 0

                  itmsh(itnm) = 1

                  itrcn(itnm) = nbsurf


                  itsmn(itnm) = itsmn(itnm) + nbsurf * 4

                  call moddas_reallocate_int(
     &                    itlmax, itnm, 6*nbsurf, itrcr, idas_itrcr)
                  itrcs(itnm) = 6 * nbsurf
                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, nbsurf, itrca, das_itrca)

               do i = 1, nbsurf

                  idas_itrcr( itrcr(itnm) + 6 * ( i - 1 )     ) = 1
                  idas_itrcr( itrcr(itnm) + 6 * ( i - 1 ) + 1 ) = 1
                  idas_itrcr( itrcr(itnm) + 6 * ( i - 1 ) + 2 )
     &                      = ntcr(i,1)
                  idas_itrcr( itrcr(itnm) + 6 * ( i - 1 ) + 3 ) = 1
                  idas_itrcr( itrcr(itnm) + 6 * ( i - 1 ) + 4 ) = 1
                  idas_itrcr( itrcr(itnm) + 6 * ( i - 1 ) + 5 )
     &                      = ntcr(i,2)

               end do

               do i = 1, nbsurf

                  das_itrca( itrca(itnm) + i - 1 ) = stca(i)

               end do

               if( kpsurf .eq. 2 ) then

                  itpan(itnm) = 2

               else if( kpsurf .eq. 3 ) then

                  itpan(itnm) = 4

               else if( kpsurf .eq. 5 ) then

                  itpan(itnm) = 6

               else

                  itpan(itnm) = 2

               end if

               do i = 1, 6

                  itpat(itnm,i,1) = iptyp(i)
                  itpat(itnm,i,2) = ipnkf(i)

               end do

                  itety(itnm) = 1
                  itenm(itnm) = neng
                  rtedl(itnm) = 0.0
                  rtemi(itnm) = steg(1)
                  rtema(itnm) = steg(neng+1)

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, neng+1, iterg, das_iterg)

                  itsmn(itnm) = itsmn(itnm) + neng + 1


               do i = 1, neng + 1

                  das_iterg( iterg(itnm) + i - 1 ) = steg(i)

               end do

                  itunt(itnm) = ksunit + 1

                  ittwo(itnm) = 3

                  itaxn(itnm) = 1
                  itaxs(itnm,1) = 1

                  itfln(itnm) = 1
                  itfll(itnm,1) = 7
                  ctfln(itnm,1) = 'fort.71'

                  itout(itnm) = 2

         end if

*-----------------------------------------------------------------------
*        t-cross  ( flux spectrum by boundary crossing)
*-----------------------------------------------------------------------

            if( italxx(8) .ne. 0 .and.
     &          iengm .eq. 0 ) then

                  write(io,'(
     &                ''energy mesh is not defined'')')
                  goto 888

            end if

            if( italxx(8) .ne. 0 .and. nbsurf .eq. 1 ) then

                  write(io,'(
     &                ''b-surf is not defined'')')
                  goto 888

            end if

*-----------------------------------------------------------------------

         if( italxx(8) .ne. 0 ) then

                  itnm = itnm + 1

                  ital( itnm ) = 2

                  itsmn(itnm) = 0
                  itstm(itnm) = 0

                  ittll(itnm) = 24
                  ittle(itnm) = '[t-cross] in region mesh'

                  itanl(itnm) = 0
                  itaxl(itnm) = 0
                  itayl(itnm) = 0
                  itazl(itnm) = 0
                  itgsh(itnm) = 0
                  itrsh(itnm) = 0
                  itmtn(itnm) = 0
                  itmcn(itnm) = 0
                  iterl(itnm) = 72
                  rtfac(itnm) = 1.0
                  itvm(itnm)  = 9
                  iteps(itnm) = 0

                  itmsh(itnm) = 1

                  itrcn(itnm) = nbsurf


                  itsmn(itnm) = itsmn(itnm) + nbsurf * 4

                  call moddas_reallocate_int(
     &                    itlmax, itnm, 6*nbsurf, itrcr, idas_itrcr)
                  itrcs(itnm) = 6 * nbsurf
                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, nbsurf, itrca, das_itrca)

               do i = 1, nbsurf

                  idas_itrcr( itrcr(itnm) + 6 * ( i - 1 )     ) = 1
                  idas_itrcr( itrcr(itnm) + 6 * ( i - 1 ) + 1 ) = 1
                  idas_itrcr( itrcr(itnm) + 6 * ( i - 1 ) + 2 )
     &                      = ntcr(i,1)
                  idas_itrcr( itrcr(itnm) + 6 * ( i - 1 ) + 3 ) = 1
                  idas_itrcr( itrcr(itnm) + 6 * ( i - 1 ) + 4 ) = 1
                  idas_itrcr( itrcr(itnm) + 6 * ( i - 1 ) + 5 )
     &                      = ntcr(i,2)

               end do

               do i = 1, nbsurf

                  das_itrca( itrca(itnm) + i - 1 ) = stca(i)

               end do

               if( kpsurf .eq. 2 ) then

                  itpan(itnm) = 2

               else if( kpsurf .eq. 3 ) then

                  itpan(itnm) = 4

               else if( kpsurf .eq. 5 ) then

                  itpan(itnm) = 6

               else

                  itpan(itnm) = 2

               end if

               do i = 1, 6

                  itpat(itnm,i,1) = iptyp(i)
                  itpat(itnm,i,2) = ipnkf(i)

               end do

                  itety(itnm) = 1
                  itenm(itnm) = neng
                  rtedl(itnm) = 0.0
                  rtemi(itnm) = steg(1)
                  rtema(itnm) = steg(neng+1)

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, neng+1, iterg, das_iterg)

                  itsmn(itnm) = itsmn(itnm) + neng + 1


               do i = 1, neng + 1

                  das_iterg( iterg(itnm) + i - 1 ) = steg(i)

               end do

                  itunt(itnm) = ksunit + 1

                  ittwo(itnm) = 3

                  itaxn(itnm) = 1
                  itaxs(itnm,1) = 1

                  itfln(itnm) = 1
                  itfll(itnm,1) = 7
                  ctfln(itnm,1) = 'fort.72'

                  itout(itnm) = 1

         end if

*-----------------------------------------------------------------------
*     read information of scoring mesh tally
*-----------------------------------------------------------------------

  510 continue

*-----------------------------------------------------------------------
*     card 37
*-----------------------------------------------------------------------

         if(     hl(icl1:icl2).eq.'mesh*'   ) then

            read(hl(11:80),*,err=803,iostat=ios) imesxx
            if( ios .eq. -1 ) goto 803

         else if(hl(icl1:icl2).eq.'all-mesh') then

            imesxx(1)=1
            imesxx(2)=1
            imesxx(3)=1
            imesxx(4)=1
            imesxx(5)=1
            imesxx(6)=1

         else if(hl(icl1:icl2).eq.'mesh'    ) then

   10       call chnumb(hl(11:80),70,icl3,icl4)
            icl3=icl3+10
            icl4=icl4+10
            if(     hl(icl3:icl4).eq.'rz-track'     ) then
                 imesxx(1)=1
            else if(hl(icl3:icl4).eq.'rz-heat'      ) then
                 imesxx(2)=1
            else if(hl(icl3:icl4).eq.'rz-current'   ) then
                 imesxx(3)=1
            else if(hl(icl3:icl4).eq.'rz-flux'      ) then
                 imesxx(4)=1
            else if(hl(icl3:icl4).eq.'rz-all'       ) then
                 imesxx(1)=1
                 imesxx(2)=1
                 imesxx(3)=1
                 imesxx(4)=1
            else if(hl(icl3:icl4).eq.'xyz-track'    ) then
                 imesxx(5)=1
            else if(hl(icl3:icl4).eq.'xyz-heat'     ) then
                 imesxx(6)=1
            else if(hl(icl3:icl4).eq.'xyz-all'      ) then
                 imesxx(5)=1
                 imesxx(6)=1
            end if

  570       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 803
            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 570
            call chnumb(hl,10,icl3,icl4)
            if(icl3.eq.0) goto 570
            if(hl(icl3:icl4) .eq. 'mesh') goto 10

            goto 590

         end if

*-----------------------------------------------------------------------
*     read rest of scoring mesh data( card 38-47 )
*-----------------------------------------------------------------------

  580    read(in,'(a)',iostat=ios) hl
         if( ios .eq. -1 ) goto 600

  590       rewind isc

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 580
            call chnumb(hl,10,icl1,icl2)
            if(icl1.eq.0) goto 580

*-----------------------------------------------------------------------
*     card 38
*-----------------------------------------------------------------------

         if(hl(icl1:icl2).eq.'r-rzmesh') then

            read(hl(11:80),*,err=803,iostat=ios) xbaser,ybaser
            if( ios .eq. -1 ) goto 803

  602       read(in,'(a)',end=691) hl

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 602
            call chnumb(hl,10,icl3,icl4)
            if(icl3.eq.0) goto 602

            if(hl(icl3:icl4).eq.'r-rzmesh') then
               read(hl(11:80),*,err=803,iostat=ios) nrmshr
            if( ios .eq. -1 ) goto 803
            else
               goto 803
            endif

  690       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 691
            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 690
            call chnumb(hl,10,icl3,icl4)
            if(icl3.eq.0) goto 690

            if(hl(icl3:icl4).eq.'r-rzmesh') then
               hl2=hl(11:80)
               call divcha(isc,hl2,70)
               goto 690
            endif

            ird = 0
            goto 692
  691       ird = 1
  692       rewind isc

            if(imesxx(1)+imesxx(2).gt.0) then

               mrmshr = nrmshr

               read(isc,*,err=803,iostat=ios) rdat1
               if( ios .eq. -1 ) goto 803
               strg(1) = rdat1

                     ic=1
               do i=1, nrmshr
                     read(isc,*,err=803,iostat=ios) idiv
                     if( ios .eq. -1 ) goto 803
                     read(isc,*,err=803,iostat=ios) rdat2
                     if( ios .eq. -1 ) goto 803
                     drdat=(rdat2-rdat1)/idiv
                  do j=1,idiv
                     ic=ic+1
                     strg(ic)=strg(ic-1)+drdat
                  end do
                     rdat1=rdat2
               end do

               nrmshr = ic

            end if

            if( ird .eq. 1 ) goto 580
            goto 590

*-----------------------------------------------------------------------
*     card 39
*-----------------------------------------------------------------------

         else if(hl(icl1:icl2).eq.'z-rzmesh') then

            read(hl(11:80),*,err=803,iostat=ios) nzmshr
            if( ios .eq. -1 ) goto 803

  400       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 402

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 400
            call chnumb(hl,10,icl3,icl4)
            if(icl3.eq.0) goto 400

            if(hl(icl3:icl4).eq.'z-rzmesh') then
               hl2=hl(11:80)
               call divcha(isc,hl2,70)
               goto 400
            endif

            ird = 0
            goto 401
  402       ird = 1
  401       rewind isc

            if(imesxx(1)+imesxx(2).gt.0) then

               mzmshr = nzmshr

               read(isc,*,err=803,iostat=ios) rdat1
               if( ios .eq. -1 ) goto 803
               stzr(1) = rdat1

                     ic=1
               do i=1,nzmshr
                     read(isc,*,err=803,iostat=ios) idiv
                     if( ios .eq. -1 ) goto 803
                     read(isc,*,err=803,iostat=ios) rdat2
                     if( ios .eq. -1 ) goto 803
                     drdat=(rdat2-rdat1)/idiv
                  do j=1,idiv
                     ic=ic+1
                     stzr(ic)=stzr(ic-1)+drdat
                  end do
                     rdat1=rdat2
               end do

               nzmshr = ic

            end if

            if( ird .eq. 1 ) goto 580
            goto 590

*-----------------------------------------------------------------------
*     card 40
*-----------------------------------------------------------------------

         else if(hl(icl1:icl2).eq.'x-xyzmesh') then

            read(hl(11:80),*,err=803,iostat=ios) nxmshx
            if( ios .eq. -1 ) goto 803

            ird = 0
  410       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 412

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 410
            call chnumb(hl,10,icl3,icl4)
            if(icl3.eq.0) goto 410

            if(hl(icl3:icl4).eq.'x-xyzmesh') then
               hl2=hl(11:80)
               call divcha(isc,hl2,70)
               goto 410
            endif

            goto 411
  412       ird = 1
  411       rewind isc

            if(imesxx(5)+imesxx(6).gt.0) then

               mxmshx = nxmshx

               read(isc,*,err=803,iostat=ios) rdat1
               if( ios .eq. -1 ) goto 803
               stxg(1) = rdat1

                     ic=1
               do i=1,nxmshx
                     read(isc,*,err=803,iostat=ios) idiv
                     if( ios .eq. -1 ) goto 803
                     read(isc,*,err=803,iostat=ios) rdat2
                     if( ios .eq. -1 ) goto 803
                     drdat=(rdat2-rdat1)/idiv
                  do j=1,idiv
                     ic=ic+1
                     stxg(ic)=stxg(ic-1)+drdat
                  end do
                     rdat1=rdat2
               end do

               nxmshx = ic

            end if

            if( ird .eq. 1 ) goto 580
            goto 590

*-----------------------------------------------------------------------
*     card 41
*-----------------------------------------------------------------------

         else if(hl(icl1:icl2).eq.'y-xyzmesh') then

            read(hl(11:80),*,err=803,iostat=ios) nymshx
            if( ios .eq. -1 ) goto 803

            ird = 0
  420       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 422

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 420
            call chnumb(hl,10,icl3,icl4)
            if(icl3.eq.0) goto 420

            if(hl(icl3:icl4).eq.'y-xyzmesh') then
               hl2=hl(11:80)
               call divcha(isc,hl2,70)
               goto 420
            endif

            goto 421
  422       ird = 1
  421       rewind isc

            if(imesxx(5)+imesxx(6).gt.0) then

               mymshx = nymshx

               read(isc,*,err=803,iostat=ios) rdat1
               if( ios .eq. -1 ) goto 803
               styg(1) = rdat1

                     ic=1
               do i=1,nymshx
                     read(isc,*,err=803,iostat=ios) idiv
                     if( ios .eq. -1 ) goto 803
                     read(isc,*,err=803,iostat=ios) rdat2
                     if( ios .eq. -1 ) goto 803
                     drdat=(rdat2-rdat1)/idiv
                  do j=1,idiv
                     ic=ic+1
                     styg(ic)=styg(ic-1)+drdat
                  end do
                     rdat1=rdat2
               end do

               nymshx = ic

            end if

            if( ird .eq. 1 ) goto 580
            goto 590

*-----------------------------------------------------------------------
*     card 42
*-----------------------------------------------------------------------

         else if(hl(icl1:icl2).eq.'z-xyzmesh') then

            read(hl(11:80),*,err=803,iostat=ios) nzmshx
            if( ios .eq. -1 ) goto 803

            ird = 0
  430       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 432

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 430
            call chnumb(hl,10,icl3,icl4)
            if(icl3.eq.0) goto 430

            if(hl(icl3:icl4).eq.'z-xyzmesh') then
               hl2=hl(11:80)
               call divcha(isc,hl2,70)
               goto 430
            endif

            goto 431
  432       ird = 1
  431       rewind isc

            if(imesxx(5)+imesxx(6).gt.0) then

               mzmshx = nzmshx

               read(isc,*,err=803,iostat=ios) rdat1
               if( ios .eq. -1 ) goto 803
               stzg(1) = rdat1

                     ic=1
               do i=1,nzmshx
                     read(isc,*,err=803,iostat=ios) idiv
                     if( ios .eq. -1 ) goto 803
                     read(isc,*,err=803,iostat=ios) rdat2
                     if( ios .eq. -1 ) goto 803
                     drdat=(rdat2-rdat1)/idiv
                  do j=1,idiv
                     ic=ic+1
                     stzg(ic)=stzg(ic-1)+drdat
                  end do
                     rdat1=rdat2
               end do

               nzmshx = ic

            end if

            if( ird .eq. 1 ) goto 580
            goto 590

*-----------------------------------------------------------------------
*     card 43
*-----------------------------------------------------------------------

         else if(hl(icl1:icl2).eq.'rsurf-mesh') then

            read(hl(11:80),*,err=803,iostat=ios) xbases,ybases
            if( ios .eq. -1 ) goto 803

            ird = 0
  603       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 792

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 603
            call chnumb(hl,10,icl3,icl4)
            if(icl3.eq.0) goto 603

            if(hl(icl3:icl4).eq.'rsurf-mesh') then
               read(hl(11:80),*,err=803,iostat=ios) nrsrfm
               if( ios .eq. -1 ) goto 803
            else
               goto 803
            endif

  790       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 792
            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 790
            call chnumb(hl,10,icl3,icl4)
            if(icl3.eq.0) goto 790

            if(hl(icl3:icl4).eq.'rsurf-mesh') then
               hl2=hl(11:80)
               call divcha(isc,hl2,70)
               goto 790
            endif

            goto 791
  792       ird = 1
  791       rewind isc

            if(imesxx(3)+imesxx(4).gt.0) then

               mrsrfm = nrsrfm

               read(isc,*,err=803,iostat=ios) rdat1
               if( ios .eq. -1 ) goto 803
               strs(1) = rdat1

                     ic=1
               do i=1, nrsrfm
                     read(isc,*,err=803,iostat=ios) idiv
                     if( ios .eq. -1 ) goto 803
                     read(isc,*,err=803,iostat=ios) rdat2
                     if( ios .eq. -1 ) goto 803
                     drdat=(rdat2-rdat1)/idiv
                  do j=1,idiv
                     ic=ic+1
                     strs(ic)=strs(ic-1)+drdat
                  end do
                     rdat1=rdat2
               end do

               nrsrfm = ic

            end if

            if( ird .eq. 1 ) goto 580
            goto 590

*-----------------------------------------------------------------------
*     card 44
*-----------------------------------------------------------------------

         else if(hl(icl1:icl2).eq.'zsurf-mesh') then

            read(hl(11:80),*,err=803,iostat=ios) nzsrfm
            if( ios .eq. -1 ) goto 803

            ird = 0
  440       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 442

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 440
            call chnumb(hl,10,icl3,icl4)
            if(icl3.eq.0) goto 440

            if(hl(icl3:icl4).eq.'zsurf-mesh') then
               hl2=hl(11:80)
               call divcha(isc,hl2,70)
               goto 440
            endif

            goto 441
  442       ird = 1
  441       rewind isc

            if(imesxx(3)+imesxx(4).gt.0) then

               mzsrfm = nzsrfm

               read(isc,*,err=803,iostat=ios) rdat1
               if( ios .eq. -1 ) goto 803
               stzs(1) = rdat1

                     ic=1
               do i=1,nzsrfm
                     read(isc,*,err=803,iostat=iios) idiv
                     if( ios .eq. -1 ) goto 803
                     read(isc,*,err=803,iostat=ios) rdat2
                     if( ios .eq. -1 ) goto 803
                     drdat=(rdat2-rdat1)/idiv
                  do j=1,idiv
                     ic=ic+1
                     stzs(ic)=stzs(ic-1)+drdat
                  end do
                     rdat1=rdat2
               end do

               nzsrfm = ic

            end if

            if( ird .eq. 1 ) goto 580
            goto 590

*-----------------------------------------------------------------------
*     card 45 : energy bin
*-----------------------------------------------------------------------

         elseif( hl(icl1:icl2) .eq. 'en-mesh' ) then

            read(hl(11:80),*) neng

            ird = 0
  710       read(in,'(a)',iostat=ios) hl
            if( ios .eq. -1 ) goto 721

            if( hl(1:1) .eq. '#' .or. hl(1:1) .eq. '%' ) goto 710
            call chnumb(hl,10,icl3,icl4)

            if( hl(icl3:icl4) .ne. 'en-mesh' ) goto 720
            write(isc,'(a)') hl(11:80)
            goto 710

            goto 720
  721       ird = 1
  720       rewind isc

            read(isc,*,iostat=ios) ( steg(i), i = 1, neng + 1 )
            if( ios .eq. -1 ) goto 725
  725       rewind isc

            do i = 1, neng
               if( steg(i+1) .le. steg(i) ) then

                  write(io,'(''energy bin must be the '',
     &                       ''ascending order.''/
     &                       '' bin='',i4,''   energy='',1p2e13.4)')
     &            i, steg(i+1), steg(i)

                  goto 888

               end if
            end do

            iengn = 1

            if( ird .eq. 1 ) goto 580
            goto 590

*-----------------------------------------------------------------------
*     card 46
*-----------------------------------------------------------------------

         else if(hl(icl1:icl2).eq.'unit-mesh') then

            read(hl(11:80),*,err=803,iostat=ios) kmunit
            if( ios .eq. -1 ) goto 803

            if( kmunit .lt. 0 .or. ksunit .gt. 2 ) then

                  write(io,'(''kmunit should be 0, 1, 2 in old input''
     &                     )')

                  goto 888

            end if

            goto 580

*-----------------------------------------------------------------------
*     card 47
*-----------------------------------------------------------------------

         else if(hl(icl1:icl2).eq.'pt-mesh') then

            read(hl(11:80),*,err=803,iostat=ios) kpmesh
            if( ios .eq. -1 ) goto 803

            goto 580

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------

  600 continue

*-----------------------------------------------------------------------
*     t-track of r-z mesh
*-----------------------------------------------------------------------

            if( imesxx(1) .ne. 0 .and.
     &        ( nrmshr .le. 1 .or. nzmshr .le. 1 ) ) then

                  write(io,'(
     &                ''r-rzmesh or z-rzmesh is not defined'')')
                  goto 888

            end if

            if( imesxx(1) .ne. 0 .and.
     &          iengn .eq. 0 ) then

                  write(io,'(
     &                ''energy mesh is not defined'')')
                  goto 888

            end if

*-----------------------------------------------------------------------

         if( imesxx(1) .ne. 0 ) then

                  inr = nrmshr - 1
                  inz = nzmshr - 1

                  itnm = itnm + 1

                  ital( itnm ) = 1

                  itsmn(itnm) = 0
                  itstm(itnm) = 0

                  ittll(itnm) = 21
                  ittle(itnm) = '[t-track] in r-z mesh'

                  itanl(itnm) = 0
                  itaxl(itnm) = 0
                  itayl(itnm) = 0
                  itazl(itnm) = 0
                  itgsh(itnm) = 0
                  itrsh(itnm) = 0
                  itmtn(itnm) = 0
                  itmcn(itnm) = 0
                  iterl(itnm) = 72
                  rtfac(itnm) = 1.0
                  itvm(itnm)  = 9
                  iteps(itnm) = 0

                  itmsh(itnm) = 2

                  rtrx0(itnm) = xbaser
                  rtry0(itnm) = ybaser

               if( mrmshr .eq. 1 ) then

                  itrty(itnm) = 2
                  itrnm(itnm) = inr

                  rtrmi(itnm) = strg(1)
                  rtrma(itnm) = strg(inr+1)
                  rtrdl(itnm) = ( strg(inr+1) - strg(1) ) / inr

               else

                  itrty(itnm) = 1
                  itrnm(itnm) = inr

                  rtrmi(itnm) = strg(1)
                  rtrma(itnm) = strg(inr+1)
                  rtrdl(itnm) = 0.0

               end if

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, inr+1, itrrg, das_itrrg)

                  itsmn(itnm) = itsmn(itnm) + inr + 1


                  do i = 1, inr + 1

                     das_itrrg( itrrg(itnm) + i - 1 ) = strg(i)

                  end do

               if( mzmshr .eq. 1 ) then

                  itzty(itnm) = 2
                  itznm(itnm) = inz

                  rtzmi(itnm) = stzr(1)
                  rtzma(itnm) = stzr(inz+1)
                  rtzdl(itnm) = ( stzr(inz+1) - stzr(1) ) / inz

               else

                  itzty(itnm) = 1
                  itznm(itnm) = inz

                  rtzmi(itnm) = stzr(1)
                  rtzma(itnm) = stzr(inz+1)
                  rtzdl(itnm) = 0.0

               end if

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, inz+1, itzrg, das_itzrg)

                  itsmn(itnm) = itsmn(itnm) + inz + 1


                  do i = 1, inz + 1

                     das_itzrg( itzrg(itnm) + i - 1 ) = stzr(i)

                  end do

               if( kpmesh .eq. 2 ) then

                  itpan(itnm) = 2

               else if( kpmesh .eq. 3 ) then

                  itpan(itnm) = 4

               else if( kpmesh .eq. 5 ) then

                  itpan(itnm) = 6

               else

                  itpan(itnm) = 2

               end if

               do i = 1, 6

                  itpat(itnm,i,1) = iptyp(i)
                  itpat(itnm,i,2) = ipnkf(i)

               end do

                  itety(itnm) = 1
                  itenm(itnm) = neng
                  rtedl(itnm) = 0.0
                  rtemi(itnm) = steg(1)
                  rtema(itnm) = steg(neng+1)

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, neng+1, iterg, das_iterg)

                  itsmn(itnm) = itsmn(itnm) + neng + 1


               do i = 1, neng + 1

                  das_iterg( iterg(itnm) + i - 1 ) = steg(i)

               end do

                  itunt(itnm) = kmunit + 1

                  ittwo(itnm) = 3

                  itaxn(itnm) = 1
                  itaxs(itnm,1) = 1

                  itfln(itnm) = 1
                  itfll(itnm,1) = 7
                  ctfln(itnm,1) = 'fort.73'

         end if

*-----------------------------------------------------------------------

         if( imesxx(2) .ne. 0 ) then

                  inr = nrmshr - 1
                  inz = nzmshr - 1

                  itnm = itnm + 1

                  ital( itnm ) = 4

                  itsmn(itnm) = 0
                  itstm(itnm) = 0

                  ittll(itnm) = 20
                  ittle(itnm) = '[t-heat] in r-z mesh'

                  itanl(itnm) = 0
                  itaxl(itnm) = 0
                  itayl(itnm) = 0
                  itazl(itnm) = 0
                  itgsh(itnm) = 0
                  itrsh(itnm) = 0
                  itmtn(itnm) = 0
                  itmcn(itnm) = 0
                  iterl(itnm) = 72
                  rtfac(itnm) = 1.0
                  itvm(itnm)  = 9
                  iteps(itnm) = 0

                  itmsh(itnm) = 2

                  rtrx0(itnm) = xbaser
                  rtry0(itnm) = ybaser

               if( mrmshr .eq. 1 ) then

                  itrty(itnm) = 2
                  itrnm(itnm) = inr

                  rtrmi(itnm) = strg(1)
                  rtrma(itnm) = strg(inr+1)
                  rtrdl(itnm) = ( strg(inr+1) - strg(1) ) / inr

               else

                  itrty(itnm) = 1
                  itrnm(itnm) = inr

                  rtrmi(itnm) = strg(1)
                  rtrma(itnm) = strg(inr+1)
                  rtrdl(itnm) = 0.0

               end if

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, inr+1, itrrg, das_itrrg)

                  itsmn(itnm) = itsmn(itnm) + inr + 1


                  do i = 1, inr + 1

                     das_itrrg( itrrg(itnm) + i - 1 ) = strg(i)

                  end do

               if( mzmshr .eq. 1 ) then

                  itzty(itnm) = 2
                  itznm(itnm) = inz

                  rtzmi(itnm) = stzr(1)
                  rtzma(itnm) = stzr(inz+1)
                  rtzdl(itnm) = ( stzr(inz+1) - stzr(1) ) / inz

               else

                  itzty(itnm) = 1
                  itznm(itnm) = inz

                  rtzmi(itnm) = stzr(1)
                  rtzma(itnm) = stzr(inz+1)
                  rtzdl(itnm) = 0.0

               end if

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, inz+1, itzrg, das_itzrg)

                  itsmn(itnm) = itsmn(itnm) + inz + 1


                  do i = 1, inz + 1

                     das_itzrg( itzrg(itnm) + i - 1 ) = stzr(i)

                  end do

                  itpan(itnm) = 0

                  itunt(itnm) = 1

                  ittwo(itnm) = 3

                  itaxn(itnm) = 1
                  itaxs(itnm,1) = 5

                  itfln(itnm) = 1
                  itfll(itnm,1) = 7
                  ctfln(itnm,1) = 'fort.76'

                  itout(itnm) = 1
                  itndy(itnm) = 28

         end if

*-----------------------------------------------------------------------

         if( imesxx(3) .ne. 0 ) then

                  inr = nrsrfm - 1
                  inz = nzsrfm - 1

                  itnm = itnm + 1

                  ital( itnm ) = 2

                  itsmn(itnm) = 0
                  itstm(itnm) = 0

                  ittll(itnm) = 21
                  ittle(itnm) = '[t-cross] in r-z mesh'

                  itanl(itnm) = 0
                  itaxl(itnm) = 0
                  itayl(itnm) = 0
                  itazl(itnm) = 0
                  itgsh(itnm) = 0
                  itrsh(itnm) = 0
                  itmtn(itnm) = 0
                  itmcn(itnm) = 0
                  iterl(itnm) = 72
                  rtfac(itnm) = 1.0
                  itvm(itnm)  = 9
                  iteps(itnm) = 0

                  itmsh(itnm) = 2

                  rtrx0(itnm) = xbases
                  rtry0(itnm) = ybases

               if( mrsrfm .eq. 1 ) then

                  itrty(itnm) = 2
                  itrnm(itnm) = inr

                  rtrmi(itnm) = strs(1)
                  rtrma(itnm) = strs(inr+1)
                  rtrdl(itnm) = ( strs(inr+1) - strs(1) ) / inr

               else

                  itrty(itnm) = 1
                  itrnm(itnm) = inr

                  rtrmi(itnm) = strs(1)
                  rtrma(itnm) = strs(inr+1)
                  rtrdl(itnm) = 0.0

               end if

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, inr+1, itrrg, das_itrrg)

                  itsmn(itnm) = itsmn(itnm) + inr + 1


                  do i = 1, inr + 1

                     das_itrrg( itrrg(itnm) + i - 1 ) = strs(i)

                  end do

               if( mzsrfm .eq. 1 ) then

                  itzty(itnm) = 2
                  itznm(itnm) = inz

                  rtzmi(itnm) = stzs(1)
                  rtzma(itnm) = stzs(inz+1)
                  rtzdl(itnm) = ( stzs(inz+1) - stzs(1) ) / inz

               else

                  itzty(itnm) = 1
                  itznm(itnm) = inz

                  rtzmi(itnm) = stzs(1)
                  rtzma(itnm) = stzs(inz+1)
                  rtzdl(itnm) = 0.0

               end if

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, inz+1, itzrg, das_itzrg)

                  itsmn(itnm) = itsmn(itnm) + inz + 1


                  do i = 1, inz + 1

                     das_itzrg( itzrg(itnm) + i - 1 ) = stzs(i)

                  end do

               if( kpmesh .eq. 2 ) then

                  itpan(itnm) = 2

               else if( kpmesh .eq. 3 ) then

                  itpan(itnm) = 4

               else if( kpmesh .eq. 5 ) then

                  itpan(itnm) = 6

               else

                  itpan(itnm) = 2

               end if

               do i = 1, 6

                  itpat(itnm,i,1) = iptyp(i)
                  itpat(itnm,i,2) = ipnkf(i)

               end do

                  itety(itnm) = 1
                  itenm(itnm) = neng
                  rtedl(itnm) = 0.0
                  rtemi(itnm) = steg(1)
                  rtema(itnm) = steg(neng+1)

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, neng+1, iterg, das_iterg)

                  itsmn(itnm) = itsmn(itnm) + neng + 1


               do i = 1, neng + 1

                  das_iterg( iterg(itnm) + i - 1 ) = steg(i)

               end do

                  itunt(itnm) = kmunit + 1

                  ittwo(itnm) = 3

                  itaxn(itnm) = 1
                  itaxs(itnm,1) = 1

                  itfln(itnm) = 1
                  itfll(itnm,1) = 7
                  ctfln(itnm,1) = 'fort.79'

                  itout(itnm) = 2

         end if

*-----------------------------------------------------------------------

         if( imesxx(4) .ne. 0 ) then

                  inr = nrsrfm - 1
                  inz = nzsrfm - 1

                  itnm = itnm + 1

                  ital( itnm ) = 2

                  itsmn(itnm) = 0
                  itstm(itnm) = 0

                  ittll(itnm) = 21
                  ittle(itnm) = '[t-cross] in r-z mesh'

                  itanl(itnm) = 0
                  itaxl(itnm) = 0
                  itayl(itnm) = 0
                  itazl(itnm) = 0
                  itgsh(itnm) = 0
                  itrsh(itnm) = 0
                  itmtn(itnm) = 0
                  itmcn(itnm) = 0
                  iterl(itnm) = 72
                  rtfac(itnm) = 1.0
                  itvm(itnm)  = 9
                  iteps(itnm) = 0

                  itmsh(itnm) = 2

                  rtrx0(itnm) = xbases
                  rtry0(itnm) = ybases

               if( mrsrfm .eq. 1 ) then

                  itrty(itnm) = 2
                  itrnm(itnm) = inr

                  rtrmi(itnm) = strs(1)
                  rtrma(itnm) = strs(inr+1)
                  rtrdl(itnm) = ( strs(inr+1) - strs(1) ) / inr

               else

                  itrty(itnm) = 1
                  itrnm(itnm) = inr

                  rtrmi(itnm) = strs(1)
                  rtrma(itnm) = strs(inr+1)
                  rtrdl(itnm) = 0.0

               end if

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, inr+1, itrrg, das_itrrg)

                  itsmn(itnm) = itsmn(itnm) + inr + 1


                  do i = 1, inr + 1

                     das_itrrg( itrrg(itnm) + i - 1 ) = strs(i)

                  end do

               if( mzsrfm .eq. 1 ) then

                  itzty(itnm) = 2
                  itznm(itnm) = inz

                  rtzmi(itnm) = stzs(1)
                  rtzma(itnm) = stzs(inz+1)
                  rtzdl(itnm) = ( stzs(inz+1) - stzs(1) ) / inz

               else

                  itzty(itnm) = 1
                  itznm(itnm) = inz

                  rtzmi(itnm) = stzs(1)
                  rtzma(itnm) = stzs(inz+1)
                  rtzdl(itnm) = 0.0

               end if

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, inz+1, itzrg, das_itzrg)

                  itsmn(itnm) = itsmn(itnm) + inz + 1


                  do i = 1, inz + 1

                     das_itzrg( itzrg(itnm) + i - 1 ) = stzs(i)

                  end do

               if( kpmesh .eq. 2 ) then

                  itpan(itnm) = 2

               else if( kpmesh .eq. 3 ) then

                  itpan(itnm) = 4

               else if( kpmesh .eq. 5 ) then

                  itpan(itnm) = 6

               else

                  itpan(itnm) = 2

               end if

               do i = 1, 6

                  itpat(itnm,i,1) = iptyp(i)
                  itpat(itnm,i,2) = ipnkf(i)

               end do

                  itety(itnm) = 1
                  itenm(itnm) = neng
                  rtedl(itnm) = 0.0
                  rtemi(itnm) = steg(1)
                  rtema(itnm) = steg(neng+1)

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, neng+1, iterg, das_iterg)

                  itsmn(itnm) = itsmn(itnm) + neng + 1


               do i = 1, neng + 1

                  das_iterg( iterg(itnm) + i - 1 ) = steg(i)

               end do

                  itunt(itnm) = kmunit + 1

                  ittwo(itnm) = 3

                  itaxn(itnm) = 1
                  itaxs(itnm,1) = 1

                  itfln(itnm) = 1
                  itfll(itnm,1) = 7
                  ctfln(itnm,1) = 'fort.80'

                  itout(itnm) = 1

         end if

*-----------------------------------------------------------------------

         if( imesxx(5) .ne. 0 ) then

                  inx = nxmshx - 1
                  iny = nymshx - 1
                  inz = nzmshx - 1

                  itnm = itnm + 1

                  ital( itnm ) = 1

                  itsmn(itnm) = 0
                  itstm(itnm) = 0

                  ittll(itnm) = 21
                  ittle(itnm) = '[t-track] in xyz mesh'

                  itanl(itnm) = 0
                  itaxl(itnm) = 0
                  itayl(itnm) = 0
                  itazl(itnm) = 0
                  itgsh(itnm) = 0
                  itrsh(itnm) = 0
                  itmtn(itnm) = 0
                  itmcn(itnm) = 0
                  iterl(itnm) = 72
                  rtfac(itnm) = 1.0
                  itvm(itnm)  = 9
                  iteps(itnm) = 0

                  itmsh(itnm) = 3

               if( mxmshx .eq. 1 ) then

                  itxty(itnm) = 2
                  itxnm(itnm) = inx

                  rtxmi(itnm) = stxg(1)
                  rtxma(itnm) = stxg(inx+1)
                  rtxdl(itnm) = ( styg(inx+1) - stxg(1) ) / inx

               else

                  itxty(itnm) = 1
                  itxnm(itnm) = inx

                  rtxmi(itnm) = stxg(1)
                  rtxma(itnm) = stxg(inx+1)
                  rtxdl(itnm) = 0.0

               end if

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, inx+1, itxrg, das_itxrg)

                  itsmn(itnm) = itsmn(itnm) + inx + 1


                  do i = 1, inx + 1

                     das_itxrg( itxrg(itnm) + i - 1 ) = stxg(i)

                  end do

               if( mymshx .eq. 1 ) then

                  ityty(itnm) = 2
                  itynm(itnm) = iny

                  rtymi(itnm) = styg(1)
                  rtyma(itnm) = styg(iny+1)
                  rtydl(itnm) = ( styg(iny+1) - styg(1) ) / iny

               else

                  ityty(itnm) = 1
                  itynm(itnm) = iny

                  rtymi(itnm) = styg(1)
                  rtyma(itnm) = styg(iny+1)
                  rtydl(itnm) = 0.0

               end if

                  ityrg(itnm) = mmmax
                  mmmax = mmmax + iny + 1
                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, iny+1, ityrg, das_ityrg)

                  itsmn(itnm) = itsmn(itnm) + iny + 1


                  do i = 1, iny + 1

                     das_ityrg( ityrg(itnm) + i - 1 ) = styg(i)

                  end do

               if( mzmshx .eq. 1 ) then

                  itzty(itnm) = 2
                  itznm(itnm) = inz

                  rtzmi(itnm) = stzg(1)
                  rtzma(itnm) = stzg(inz+1)
                  rtzdl(itnm) = ( stzg(inz+1) - stzg(1) ) / inz

               else

                  itzty(itnm) = 1
                  itznm(itnm) = inz

                  rtzmi(itnm) = stzg(1)
                  rtzma(itnm) = stzg(inz+1)
                  rtzdl(itnm) = 0.0

               end if

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, inz+1, itzrg, das_itzrg)

                  itsmn(itnm) = itsmn(itnm) + inz + 1


                  do i = 1, inz + 1

                     das_itzrg( itzrg(itnm) + i - 1 ) = stzg(i)

                  end do

               if( kpmesh .eq. 2 ) then

                  itpan(itnm) = 2

               else if( kpmesh .eq. 3 ) then

                  itpan(itnm) = 4

               else if( kpmesh .eq. 5 ) then

                  itpan(itnm) = 6

               else

                  itpan(itnm) = 2

               end if

               do i = 1, 6

                  itpat(itnm,i,1) = iptyp(i)
                  itpat(itnm,i,2) = ipnkf(i)

               end do

                  itety(itnm) = 1
                  itenm(itnm) = neng
                  rtedl(itnm) = 0.0
                  rtemi(itnm) = steg(1)
                  rtema(itnm) = steg(neng+1)

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, neng+1, iterg, das_iterg)

                  itsmn(itnm) = itsmn(itnm) + neng + 1


               do i = 1, neng + 1

                  das_iterg( iterg(itnm) + i - 1 ) = steg(i)

               end do

                  itunt(itnm) = kmunit + 1

                  ittwo(itnm) = 3

                  itaxn(itnm) = 1
                  itaxs(itnm,1) = 1

                  itfln(itnm) = 1
                  itfll(itnm,1) = 7
                  ctfln(itnm,1) = 'fort.81'

         end if

*-----------------------------------------------------------------------

         if( imesxx(6) .ne. 0 ) then

                  inx = nxmshx - 1
                  iny = nymshx - 1
                  inz = nzmshx - 1

                  itnm = itnm + 1

                  ital( itnm ) = 4

                  itsmn(itnm) = 0
                  itstm(itnm) = 0

                  ittll(itnm) = 20
                  ittle(itnm) = '[t-heat] in xyz mesh'

                  itanl(itnm) = 0
                  itaxl(itnm) = 0
                  itayl(itnm) = 0
                  itazl(itnm) = 0
                  itgsh(itnm) = 0
                  itrsh(itnm) = 0
                  itmtn(itnm) = 0
                  itmcn(itnm) = 0
                  iterl(itnm) = 72
                  rtfac(itnm) = 1.0
                  itvm(itnm)  = 9
                  iteps(itnm) = 0

                  itmsh(itnm) = 3

               if( mxmshx .eq. 1 ) then

                  itxty(itnm) = 2
                  itxnm(itnm) = inx

                  rtxmi(itnm) = stxg(1)
                  rtxma(itnm) = stxg(inx+1)
                  rtxdl(itnm) = ( styg(inx+1) - stxg(1) ) / inx

               else

                  itxty(itnm) = 1
                  itxnm(itnm) = inx

                  rtxmi(itnm) = stxg(1)
                  rtxma(itnm) = stxg(inx+1)
                  rtxdl(itnm) = 0.0

               end if

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, inx+1, itxrg, das_itxrg)

                  itsmn(itnm) = itsmn(itnm) + inx + 1


                  do i = 1, inx + 1

                     das_itxrg( itxrg(itnm) + i - 1 ) = stxg(i)

                  end do

               if( mymshx .eq. 1 ) then

                  ityty(itnm) = 2
                  itynm(itnm) = iny

                  rtymi(itnm) = styg(1)
                  rtyma(itnm) = styg(iny+1)
                  rtydl(itnm) = ( styg(iny+1) - styg(1) ) / iny

               else

                  ityty(itnm) = 1
                  itynm(itnm) = iny

                  rtymi(itnm) = styg(1)
                  rtyma(itnm) = styg(iny+1)
                  rtydl(itnm) = 0.0

               end if

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, iny+1, ityrg, das_ityrg)

                  itsmn(itnm) = itsmn(itnm) + iny + 1


                  do i = 1, iny + 1

                     das_ityrg( ityrg(itnm) + i - 1 ) = styg(i)

                  end do

               if( mzmshx .eq. 1 ) then

                  itzty(itnm) = 2
                  itznm(itnm) = inz

                  rtzmi(itnm) = stzg(1)
                  rtzma(itnm) = stzg(inz+1)
                  rtzdl(itnm) = ( stzg(inz+1) - stzg(1) ) / inz

               else

                  itzty(itnm) = 1
                  itznm(itnm) = inz

                  rtzmi(itnm) = stzg(1)
                  rtzma(itnm) = stzg(inz+1)
                  rtzdl(itnm) = 0.0

               end if

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, inz+1, itzrg, das_itzrg)

                  itsmn(itnm) = itsmn(itnm) + inz + 1


                  do i = 1, inz + 1

                     das_itzrg( itzrg(itnm) + i - 1 ) = stzg(i)

                  end do

                  itpan(itnm) = 0

                  itety(itnm) = 1
                  itenm(itnm) = neng
                  rtedl(itnm) = 0.0
                  rtemi(itnm) = steg(1)
                  rtema(itnm) = steg(neng+1)

                  call moddas_reallocate_dbl(
     &                    itlmax, itnm, neng+1, iterg, das_iterg)

                  itsmn(itnm) = itsmn(itnm) + neng + 1


               do i = 1, neng + 1

                  das_iterg( iterg(itnm) + i - 1 ) = steg(i)

               end do

                  itunt(itnm) = 1

                  ittwo(itnm) = 3

                  itaxn(itnm) = 1
                  itaxs(itnm,1) = 5

                  itfln(itnm) = 1
                  itfll(itnm,1) = 7
                  ctfln(itnm,1) = 'fort.86'

                  itout(itnm) = 1
                  itndy(itnm) = 28

         end if

*-----------------------------------------------------------------------

         close(isc)

  520 return

*-----------------------------------------------------------------------

  803 continue

         write(io,'(
     &        ''tally description is wrong in scoring mesh'')')
         goto 888

*-----------------------------------------------------------------------

  804 continue

         write(io,'(
     &        ''tally description is wrong in region mesh'')')
         goto 888

*-----------------------------------------------------------------------

  950 continue

         write(io,'(
     &        ''Total tally storage number exceeds mdas'',i9)')
     &         mdas
         goto 888

*-----------------------------------------------------------------------

  990 continue

         write(io,'(''Sequential line of definitin is too long'')')
         goto 888

*-----------------------------------------------------------------------

  997 continue

         write(io,'(''Error : illegal source option was detected'',
     &              '' in old input'')')
         goto 888

*-----------------------------------------------------------------------

  888 continue

         ierr = 888

*-----------------------------------------------------------------------

      return
      end


