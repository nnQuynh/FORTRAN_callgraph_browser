************************************************************************
*                                                                      *
      module TDCHAINMOD
*                                                                      *
*     Modules for tdchain to increase number of mesh                   *
*
*     Following subroutines are called externally
*       tdchain
*       pdchreg2
*       tdchech
*       talldcinit
*       talldcfin
*
*     Created by T. Furuta on 2020/05/22
*     Last modified 2020/05/22 by T. Furuta
*                                                                      *
************************************************************************
      integer :: itdc,itdcreg,iregitm,itdcregitm
      integer :: maxtdcreg,maxregitm,maxtdcregitm,maxreg
      integer,parameter :: maxitm = 100, nientry = 53, nrentry = 10
      integer,parameter :: ncentry = 2 ! T.Sato 2025/02/12
      integer,parameter :: maxnmtcf = 10
      integer,parameter :: zero = 0

      integer,allocatable :: itnm2tdc(:),itdc2reg(:)

      character(8),allocatable :: tglst(:),tglsttmp(:)
      real(8),allocatable :: tgnlt(:),tgnlttmp(:)

      real(8),allocatable :: rtgvlm(:),rtgvlmtmp(:)
      integer,allocatable :: itgcel(:),itgceltmp(:)
      integer,allocatable :: itgnum(:),itgnumtmp(:)
      integer,allocatable :: ireg2itm(:),ireg2itmtmp(:)
      real(8),allocatable :: fflux(:),ffluxtmp(:)

      real(8),allocatable :: tgnzas(:),tgnzastmp(:)
      real(8),allocatable :: ctgnnds(:),ctgnndstmp(:)
      character(8),allocatable :: htgnzas(:),htgnzastmp(:)
      integer,allocatable :: jreg2itm(:),jreg2itmtmp(:)

      integer,allocatable :: itgcell(:),itgcelltmp(:)
      integer,allocatable :: itglist(:),itglisttmp(:)
      real(8),allocatable :: rtgvol(:),rtgvoltmp(:)
      character(200),allocatable :: dcrgm(:),dcrgmtmp(:)

      integer,parameter :: n1next=1000,n2next=1000,n3next=1000

*------------------------------------------------------------------------
      contains

************************************************************************
      subroutine talldcinit(n)
************************************************************************
      integer,intent(in) :: n
      itdcreg=0
      itdcregitm=0
      allocate( itnm2tdc(n),itdc2reg(n) )
      call allocate_talldc1
      call allocate_talldc2
      call allocate_talldc3
      call allocate_talldc4

      end subroutine talldcinit

************************************************************************
      subroutine talldcfin
************************************************************************
      deallocate( itnm2tdc,itdc2reg )
      call deallocate_talldc1
      call deallocate_talldc2
      call deallocate_talldc3
      call deallocate_talldc4

      end subroutine talldcfin

************************************************************************
      subroutine allocate_talldc1
************************************************************************
      maxtdcregitm=n3next
      allocate( tglst(maxtdcregitm),tgnlt(maxtdcregitm) )
      end subroutine allocate_talldc1

************************************************************************
      subroutine allocate_talldc2
************************************************************************
      maxtdcreg=n2next
      allocate( itgcel(maxtdcreg),rtgvlm(maxtdcreg),itgnum(maxtdcreg) )
      allocate( ireg2itm(maxtdcreg) )
      allocate( fflux(maxtdcreg) )
      end subroutine allocate_talldc2

************************************************************************
      subroutine allocate_talldc3
************************************************************************
      maxregitm=n2next
      allocate( tgnzas(maxregitm),ctgnnds(maxregitm) )
      allocate( htgnzas(maxregitm) )
      tgnzas(:) = 0.0d0
      end subroutine allocate_talldc3

************************************************************************
      subroutine allocate_talldc4
************************************************************************
      maxreg=n1next
      allocate( itgcell(maxreg),itglist(maxreg),rtgvol(maxreg) )
      allocate( dcrgm(maxreg) )
      allocate( jreg2itm(maxreg) )
      rtgvol(1:maxreg)=1.0
      end subroutine allocate_talldc4

************************************************************************
      subroutine deallocate_talldc1
************************************************************************
      deallocate( tglst,tgnlt )
      end subroutine deallocate_talldc1

************************************************************************
      subroutine deallocate_talldc2
************************************************************************
      deallocate( itgcel,rtgvlm,itgnum )
      deallocate( ireg2itm )
      deallocate( fflux )
      end subroutine deallocate_talldc2

************************************************************************
      subroutine deallocate_talldc3
************************************************************************
      deallocate( tgnzas,ctgnnds,htgnzas )
      end subroutine deallocate_talldc3

************************************************************************
      subroutine deallocate_talldc4
************************************************************************
      deallocate( itgcell,itglist,rtgvol,dcrgm,jreg2itm )
      end subroutine deallocate_talldc4

************************************************************************
      subroutine reallocate_talldc1(n)
************************************************************************
      implicit none
      integer,intent(in) :: n
      integer maxtdcregitm0
      allocate( tglsttmp(maxtdcregitm) )
      allocate( tgnlttmp(maxtdcregitm) )
      tglsttmp(1:maxtdcregitm)=tglst(1:maxtdcregitm)
      tgnlttmp(1:maxtdcregitm)=tgnlt(1:maxtdcregitm)
      deallocate( tglst,tgnlt )
      maxtdcregitm0=maxtdcregitm
      do while(maxtdcregitm.lt.n)
       maxtdcregitm=maxtdcregitm+n3next
      enddo
      allocate( tglst(maxtdcregitm),tgnlt(maxtdcregitm) )
      tglsttmp(1:maxtdcregitm0)=tglst(1:maxtdcregitm0)
      tgnlttmp(1:maxtdcregitm0)=tgnlt(1:maxtdcregitm0)
      deallocate( tglsttmp,tgnlttmp )
      end subroutine reallocate_talldc1

************************************************************************
      subroutine reallocate_talldc2(n)
************************************************************************
      implicit none
      integer,intent(in) :: n
      integer maxtdcreg0
      allocate( itgceltmp(maxtdcreg) )
      allocate( rtgvlmtmp(maxtdcreg) )
      allocate( itgnumtmp(maxtdcreg) )
      allocate( ireg2itmtmp(maxtdcreg) )
      allocate( ffluxtmp(maxtdcreg) )
      itgceltmp(1:maxtdcreg)=itgcel(1:maxtdcreg)
      rtgvlmtmp(1:maxtdcreg)=rtgvlm(1:maxtdcreg)
      itgnumtmp(1:maxtdcreg)=itgnum(1:maxtdcreg)
      ireg2itmtmp(1:maxtdcreg)=ireg2itm(1:maxtdcreg)
      ffluxtmp(1:maxtdcreg)=fflux(1:maxtdcreg)
      deallocate( itgcel,rtgvlm,itgnum,ireg2itm,fflux )
      maxtdcreg0=maxtdcreg
      do while(maxtdcreg.lt.n)
       maxtdcreg=maxtdcreg+n2next
      enddo
      allocate( itgcel(maxtdcreg),rtgvlm(maxtdcreg),itgnum(maxtdcreg) )
      allocate( ireg2itm(maxtdcreg) )
      allocate( fflux(maxtdcreg) )
      itgcel(1:maxtdcreg0)=itgceltmp(1:maxtdcreg0)
      rtgvlm(1:maxtdcreg0)=rtgvlmtmp(1:maxtdcreg0)
      itgnum(1:maxtdcreg0)=itgnumtmp(1:maxtdcreg0)
      ireg2itm(1:maxtdcreg0)=ireg2itmtmp(1:maxtdcreg0)
      fflux(1:maxtdcreg0)=ffluxtmp(1:maxtdcreg0)
      deallocate( itgceltmp,rtgvlmtmp,itgnumtmp )
      deallocate( ireg2itmtmp )
      deallocate( ffluxtmp )
      if(n.gt.maxtdcregitm)call reallocate_talldc1(n) !FURUTA20200630
      end subroutine reallocate_talldc2

************************************************************************
      subroutine reallocate_talldc3(n)
************************************************************************
      implicit none
      integer,intent(in) :: n
      integer maxregitm0
      allocate( tgnzastmp(maxregitm),ctgnndstmp(maxregitm) )
      allocate( htgnzastmp(maxregitm) )
      tgnzastmp(1:maxregitm)=tgnzas(1:maxregitm)
      ctgnndstmp(1:maxregitm)=ctgnnds(1:maxregitm)
      htgnzastmp(1:maxregitm)=htgnzas(1:maxregitm)
      deallocate( tgnzas,ctgnnds,htgnzas )
      maxregitm0=maxregitm
      do while(maxregitm.lt.n)
       maxregitm=maxregitm+n2next
      enddo
      allocate( tgnzas(maxregitm),ctgnnds(maxregitm) )
      allocate( htgnzas(maxregitm) )
      tgnzas(1:maxregitm0)=tgnzastmp(1:maxregitm0)
      ctgnnds(1:maxregitm0)=ctgnndstmp(1:maxregitm0)
      htgnzas(1:maxregitm0)=htgnzastmp(1:maxregitm0)
      deallocate( tgnzastmp,ctgnndstmp,htgnzastmp )
      if(n.gt.maxtdcregitm)call reallocate_talldc1(n) !FURUTA20200630
      end subroutine reallocate_talldc3

************************************************************************
      subroutine reallocate_talldc4(n)
************************************************************************
      implicit none
      integer,intent(in) :: n
      integer maxreg0
      allocate( itgcelltmp(maxreg),itglisttmp(maxreg) )
      allocate( rtgvoltmp(maxreg) )
      allocate( dcrgmtmp(maxreg) )
      allocate( jreg2itmtmp(maxreg))
      itgcelltmp(1:maxreg)=itgcell(1:maxreg)
      itglisttmp(1:maxreg)=itglist(1:maxreg)
      rtgvoltmp(1:maxreg)=rtgvol(1:maxreg)
      dcrgmtmp(1:maxreg)=dcrgm(1:maxreg)
      jreg2itmtmp(1:maxreg)=jreg2itm(1:maxreg)
      deallocate( itgcell,itglist,rtgvol,dcrgm,jreg2itm )
      maxreg0=maxreg
      do while(maxreg.lt.n)
       maxreg=maxreg+n1next
      enddo
      allocate( itgcell(maxreg),itglist(maxreg),rtgvol(maxreg) )
      allocate( dcrgm(maxreg) )
      allocate( jreg2itm(maxreg))
      itgcell(1:maxreg0)=itgcelltmp(1:maxreg0)
      itglist(1:maxreg0)=itglisttmp(1:maxreg0)
      rtgvol(1:maxreg0)=rtgvoltmp(1:maxreg0)
      rtgvol(maxreg0+1:maxreg)=1.0
      dcrgm(1:maxreg0)=dcrgmtmp(1:maxreg0)
      jreg2itm(1:maxreg0)=jreg2itmtmp(1:maxreg0)
      deallocate( itgcelltmp,itglisttmp,rtgvoltmp,dcrgmtmp,jreg2itmtmp )
      if(n.gt.maxtdcregitm)call reallocate_talldc1(n) !FURUTA20200630
      if(n.gt.maxtdcreg)call reallocate_talldc2(n)    !FURUTA20200630
      if(n.gt.maxregitm)call reallocate_talldc3(n)    !FURUTA20200630
      end subroutine reallocate_talldc4

************************************************************************
*                                                                      *
      subroutine tdchain(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [t-dchain] section of input files                         *
*       based on [t-yield] section and additionally-make [t-dchain]    *
*            modified by N.Matsuda on 2015/01/16                       *
*       last modified by T.Furuta on 2020/05/05                        *
*                                                                      *
************************************************************************
      use sangelmod, only: itsanf
      use anatallymod, only: nanatalRead ! S.H. 2022.1.21
      use partmod ! frtati 2021/10/05
      use moddas
      use moddas_mesh
      use moddas_region
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
*     for mesh
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

      common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
      common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &                rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)
      common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &                rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)

*    &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)

      common /tall08/ rtrx0(itlmax), rtry0(itlmax)

      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)

      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax),
     &                ihxll(itlmax), ihxle(itlmax),
     &                ihdll(itlmax), ihdle(itlmax),
     &                ihfll(itlmax), ihfle(itlmax),
     &                ihsll(itlmax), ihsle(itlmax)
      character ittle*80, ihxle*80, ihdle*80, ihfle*80, ihsle*80
      integer hxsll, hdcll, hfsll, hssll

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall19/ itsmn(itlmax), itstm(itlmax)

*    &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)

      common /tall21/ rtfac(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)

      common /tall28/ itgsh(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall30/ itnda(itlmax)
      common /tall31/ iterl(itlmax)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall36/ itdpo(itlmax)

      common /tall37/ itcnt(9,itlmax)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)

*    &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
*    &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall49/ itglt(itlmax)

      common /tall52/ itprd(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /paraj/  mstz(300), parz(300) ! S.H. 2022.1.21

!OBINATA: for resfile
      common /tall59/ irfll(itlmax), crfln(itlmax), itrff(itlmax)
      character crfln*100
      character irfile*100
      common /talldc2/ itarget
      common /talldc/ itdci(nientry,itlmax), rtdci(nrentry,itlmax),
     &               tbina(itlmax,maxitm), tbinb(itlmax,maxitm),
     &               beampw(itlmax,maxitm),
     &               tmina(itlmax,maxitm), tminb(itlmax,maxitm),
     &               aclst(itlmax,maxitm), bqlst(itlmax,maxitm)
      common /talldc3/ ctdci(ncentry,itlmax) ! T.Sato, character variable in [t-dchain]
      character ctdci*200

CCSE added for dchain parameter, mtscore (2018.07.31) >>>>>
      common /tall77/ itscore(itlmax)
CCSE added for dchain parameter, mtscore (2018.07.31) <<<<<

      common /tall83/ itnzn(itlmax), itndm(itlmax) ! frtati 2022/02/18

      common /rdfisparm/ irdfisyd !H.Ratliff 2021.05.19

      common /paran/  icfn(100), ilfn(100), chfn(100)
      character       chfn*200
      integer icfn,ilfn
      character line*256

      dimension ancza(maxitm)
      character hancza(maxitm)*8
      character hl*8
      character chaia3*3
      character chaia2*2
      character chaia1*1
      character chaim*1
      character aclst*8


      common /redufmt/ iredufmt(itlmax) !FURUTA20200601

      integer itetp, iottp, iancza, itgnum2



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
      character hnxslib*80
      character hdcylib*80
      character hnfylib*80
      character hsfylib*80

      real*8 radlev,ergthfa,ergfafu,chrlvth

*-----------------------------------------------------------------------
      parameter ( nkwarga = 109 )
CCSE added for dchain parameter (2018.07.31) >>>>>
      dimension lschn(nkwarga), ischn(nkwarga)
      character schan(nkwarga)*8

      data icsu / nkwarga /

      data ( schan(i), i = 1, nkwarga ) /
     &    'mesh    ','special ','mother  ','nucleus ','axis    ',
     &    'file    ','title   ','angel   ','unit    ','info    ',
     &    '2d-type ','factor  ','material','x-txt   ','y-txt   ',
     &    'z-txt   ','gshow   ','rshow   ','ndata   ','iechrl  ',
     &    'volmat  ','epsout  ','ctmin(1)','ctmax(1)','ctmin(2)',
     &    'ctmax(2)','ctmin(3)','ctmax(3)','resol   ','width   ',
     &    'trcl    ','*trcl   ','part    ','gslat   ','output  ',
     &    'resfile ','imode   ','jmode   ','itstep  ','itout   ',
     &    'idivs   ','iregon  ','inmtcf  ','ichain  ','itdecs  ',
     &    'itdecn  ','isomtr  ','ifisyd  ','ifisye  ','iyild   ',
     &    'iggrp   ','ibetap  ','acmin   ','istabl  ','igsdef  ',
     &    'iprtb1  ','iprtb2  ','rprtb2  ','iprtb3  ','igsorg  ',
     &    'amp     ','ebeam   ','prodnp  ','timeevo ','outtime ',
     &    'ac-list ','target  ','dversion','stdcut  ','sangel  ',
     &    'mtscore ','iertdcho','inxslib ','idcylib ','infylib ',
     &    'hnxslib ','hdcylib ','hnfylib ','iwrtchn ','chrlvth ',
     &    'iwrchdt ','iwrchss ','ixsrall ','idosecf ',
     &    'ipltmode','ipltaxis','foamout ','foamvals','iredufmt',
     &    'irdonce ','imtcard ','imtcnum ','imtcmeta','thmatnd ',
     &    'idosunit','iangbpwr','iphtout ','isfylib ','radlev  ',
     &    'ergthfa ','ergfafu ','infyerg ','hsfylib ','ilchain ',
     &    'mxnuclei','itriton ','iaonucl ','aonucl  ','aoreg   '/

      data ( lschn(i), i = 1, nkwarga ) /
     &     4,         7,         6,         7,         4,
     &     4,         5,         5,         4,         4,
     &     7,         6,         8,         5,         5,
     &     5,         5,         5,         5,         6,
     &     6,         6,         8,         8,         8,
     &     8,         8,         8,         5,         5,
     &     4,         5,         4,         5,         6,
     &     7,         5,         5,         6,         5,
     &     5,         6,         6,         6,         6,
     &     6,         6,         6,         6,         5,
     &     5,         6,         5,         6,         6,
     &     6,         6,         6,         6,         6,
     &     3,         5,         6,         7,         7,
     &     7,         6,         8,         6,         6,
     &     7,         8,         7,         7,         7,
     &     7,         7,         7,         7,         7,
     &     7,         7,         7,         7,
     &     8,         8,         7,         8,         8,
     &     7,         7,         7,         8,         7,
     &     8,         8,         7,         7,         6,
     &     7,         7,         7,         7,         7,
     &     8,         7,         7,         6,         5/
CCSE added for dchain parameter (2018.07.31) <<<<<

*-----------------------------------------------------------------------



*-----------------------------------------------------------------------

      dimension iaxis(6)
      data ( iaxis(i), i = 1, 1 ) /
     &    13/
      dimension iaxdc(6)
      data ( iaxdc(i), i = 1, 6 ) /
     &    15,        16,        17,        18,        19,
     &    20/
      character ifile(6)*100
      dimension lfile(6)
      character filnm*100

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      character dkam*9

      character tnamey*9
      character tnamed*10
      character tnamet*9
      data      tnamey /'[t-yield]'/
      data      tnamed /'[t-dchain]'/
      data      tnamet /'[t-track]'/

      character dtrack*11
      data      dtrack /'DCTrack.dat'/

      character yield*5
      data      yield /'.dyld'/
      character dcout*5
      data      dcout /'.dout'/
      character track*5
      data      track /'.dtrk'/

      character aonucl*200, aoreg*200 ! T.Sato 2025/02/12

      dimension temsh(4,maxitm), itemsh(2,maxitm)

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

      common /subtra/ isubt, ipsub(mxpart)  ! kitamura22/03/31

*-----------------------------------------------------------------------

      character element(104)*3,cnuc*3,elmnt(104)*3

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

      data elmnt /
     &    ' H ','He ','Li ','Be ',' B ',' C ',' N ',' O ',' F ','Ne ',
     &    'Na ','Mg ','Al ','Si ',' P ',' S ','Cl ','Ar ',' K ','Ca ',
     &    'Sc ','Ti ',' V ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ',' Y ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ',' I ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ',' W ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ',' U ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

*-----------------------------------------------------------------------

      itdc=itdc+1            !FURUTA20200520
      itdc2reg(itdc)=itdcreg !FURUTA20200522
      iregitm=0              !FURUTA20200522
*-----------------------------------------------------------------------
*     default values of [T-Yield] parameters
*-----------------------------------------------------------------------

            ierr  = 0

            inmat = 0
            inaxi = 1
           inaxdc = 6
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
            iout  = 1  ! MATSUDA 2017.03.21 ! 0: product, 1: cutoff

            langel = 0
            itsanf = 0
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
            igslt  = -1
            lrfile = 0 !OBINATA
            irfflg = 0 !OBINATA
            imxnuc = -1 !frtati 2022/03/25

            icount(1) = 0
            icount(2) = 0
            icount(3) = 0
            icount(4) = -9999
            icount(5) =  9999
            icount(6) = -9999
            icount(7) =  9999
            icount(8) = -9999
            icount(9) =  9999
            ndata  = 2  ! S.H. 2024/12/19, default value is change to 2.

            ireso = 1
            width = 0.5

            rfact = 1.0

         do i = 1, icsu

            ischn(i) = 0

         end do

         do i = 1, maxreg !FURUTA20200522
            rtgvol(i) = 1.0
         end do

ccse 2023/03/29, zero set itgvoll
         itgvoll = 0

            stdcut = -1.0

*-----------------------------------------------------------------------
*     default values of [T-Dchain] parameters
*-----------------------------------------------------------------------

*     Card[1] title
            titll = 0
*     Card[2] primary control data
            imode = 2
            jmode = 2
*     Card[3] first calculation condition data
           idcstp = 1
           idcout = 1
            idivs = 50
           iregon = 1
           inmtcf = 1
           ichain = 100
           itdecs = 1
           itdecn = 1
           isomtr = 2
           ifisyd = 0
           ifisye = 0
*     Card[4] output condition data
            iyild = 2
            iggrp = 3
           ibetap = 1
            acmin = 0.0  ! T.Sato 2025/01/31 change 1e-20 to 0.0 (down to 1e-10 of activity)
           istabl = 0
           igsdef = 0   ! T.Sato 2025/02/10 change 1 to 0
           iprtb1 = 1
           iprtb2 = 1
           rprtb2 = 10
           iprtb3 = 0
           igsorg = 0   ! T.Sato 2025/02/10 change 1 to 0
*     Card[5] spallation data
c iwamoto 2012/12/10 amp=source/sec
              amp = 1.0
             ampp = 1.0
            ebeam = 3.0
           prodnp = 1.0
*     Card[6] irra. / cooling info.
            itetp = 0
*     Card[7] output time
            iottp = 0
*     Card[8] pick up nuclide
           iancza = 0
*     Card[10] neutron flux
           ifilnm = 0

            klst=0
         iversion = 1

CCSE added for r-z/xyz parameter (2018.07.31) >>>>>
          mtscore = 0
CCSE added for r-z/xyz parameter (2018.07.31) <<<<<

*     Card[9] toggle of error propagation in DCHAIN
         iertdcho = 1
*     Card[3b] data library selection
         inxslib  = 100
         idcylib  = 5
         infylib  = 17 ! H.Ratliff 2021.04.07 - re-enable fission
         write(hnxslib,"('custom-lib',70x)")
         write(hdcylib,"('custom-lib',70x)")
         write(hnfylib,"('custom-lib',70x)")
         write(hsfylib,"('custom-lib',70x)")
         hxsll = 10
         hdcll = 10
         hfsll = 10
         hssll = 10
*     Card[4] output options for chain data file
         iwrtchn  = 1
         chrlvth  = -1
         iwrchdt  = 0
         iwrchss  = 0
*     Card[3] calculation option for loading all cross sections into memory
          ixsrall = 1
*     Card[4] output options for dose rate coefficient
          idosecf = 1
*     Card[9d] output plot options for mesh geometries
        ipltmode = 0
        ipltaxis = 1
*     Card[9d] output plot options for tetrahedral geometries
        foamout  = 0
        foamvals = 0
        iredufmt0 = 1
        irdonce = 1
        imtcard = 1
        imtcnum = 0
        imtcmeta = 2
        thmatnd = 1.d-6
        idosunit = 333
        iangbpwr = 1
        iphtout = 0  ! T.Sato 2025/01/31 change from 1 (gamma spectrum) to 0 (no output)
        itriton = -1  ! T.Sato 2025/02/08 tritium-production cross section (0: consider, 1:ignore when HE library avaiable, 2: always ignore)
        iaonucl = 0   ! T.Sato 2025/02/12 for angelout_nuclides
        laonucl = 0 ! T.Sato 2025/02/12 for character length of angelout_nuclides data
        laoreg = 0  ! T.Sato 2025/02/12 for character length of angelout_region data
       isfylib = 1
        ergthfa = 0.1
        ergfafu = 10.0
        radlev  = 1.0e-3 ! T.Sato 2025/02/02 change from 0 to 1.0 Bq
        infyerg = 4
        ilchain = 100

! H.Ratliff 2021.05.19 - determines if ifisyd manually specified in T-Dchain or left at default
        irdfisyd = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 800

               if( iskip .ne. 0 ) goto 140

  150 continue

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 800
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

!-----------------------------------------------------------------------
!        sumtally sub section
!-----------------------------------------------------------------------
      if ( chlw(i1:i3) .eq. 'sumtally start' ) then
       isumtalRead = 1
       call read_sumtal_dchain(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
       if( jpn == 3 ) goto 800   ! 'q:'
       if( jpn == 1 ) goto 800   ! end of the section
       go to 140
      end if

*-----------------------------------------------------------------------
*        anatally
*-----------------------------------------------------------------------

      if ( chlw(i1:i3) .eq. 'anatally start' ) then

         if (  mstz(1) .eq. 17 ) then
            nanatalRead(1+itnm) = 1
            nanatalRead(1+itnm+1) = 1
            nanatalRead(1+itnm+2) = 1
            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &           jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
            if ( chlw(i1:i3) .ne. 'anatally end' ) then
               ierr = 1
            end if
         else

            istbat = 1

            call stbatm(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                  ist_nm,isd_nm,kst_nm)

            if( ierr .ne. 0 ) return
            if( jpn == 3 ) goto 800   ! 'q:'
            if( jpn == 1 ) goto 800   ! end of the section

         end if ! S.H. (2020.2.25)
            goto 140

      end if

!-----------------------------------------------------------------------

  200    continue

            chlc = chlw
            call chcomp(chlc,icl,i3,i5)

         do i = 1, icsu

            il = icl + lschn(i) - 1

            if( chlc(icl:il) .eq. schan(i)(1:lschn(i)) ) goto 100

         end do

               goto 987

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
*        mesh = region only! (r-z, xyz or tet) 2020.05.05
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

*    &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
*    &                      ntrn,mtrn,ndsm,nvol,ivl,irvl,0)
               ndsm = 1
               call moddas_allocate_int(
     &                 MAX_NUM_ITREG, idas_itreg_temporary)

               call tregdc(1,jsn,jsi,dsin,idsi,ill,ilf,
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                      ntrn,mtrn,ndsm,nvol,ivl,irvl,0,
     &                      idcrg,MAX_NUM_ITREG,idas_itreg_temporary)

               if( mtrn > MAX_NUM_ITREG ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.tdchain@tallsm4.f'
     &                    //' ?dimension over idas_itreg_temporary?'
     &                    //' mtrn > MAX_NUM_ITREG'
     &                 ,' (mtrn=',mtrn,')'
     &                 ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                  ErrID = 'L:930/R:tdchain/F:tallsm4.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

                  if( ierr .ne. 0 ) return
                  if( jpn  .eq. 3 ) goto 800
                  if( ntrn .lt. -1 ) goto 996

            else if( imesh .eq. 2 ) then

CCSE change for mesh=r-z parameter (2017.11.30) >>>>>
               call trzmeshdc(jsn,jsi,dsin,idsi,ill,ilf,
     &                       jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                       rzx0,rzy0,
     &                       irtp,inr,rmin,rmax,rdel,istrg,
     &                       iztp,inz,zmin,zmax,zdel,istzg,
     &                       idcrg)

CCSE add for mesh=r-z parameter (2018.07.31) >>>>>
               idcrg = idcrg - 1
CCSE add for mesh=r-z parameter (2018.07.31) <<<<<
               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 800
               if( irtp .lt. 0 ) goto 996
               if( iztp .lt. 0 ) goto 996
CCSE change for mesh=r-z parameter (2017.11.30) <<<<<

            else if( imesh .eq. 3 ) then

CCSE change for mesh=xyz parameter (2017.11.30) >>>>>
               call txymeshdc(jsn,jsi,dsin,idsi,ill,ilf,
     &                       jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                       ixtp,inx,xmin,xmax,xdel,istxg,
     &                       iytp,iny,ymin,ymax,ydel,istyg,
     &                       iztp,inz,zmin,zmax,zdel,istzg,
     &                       idcrg)

CCSE add for mesh=r-z parameter (2018.07.31) >>>>>
               idcrg = idcrg - 1
CCSE add for mesh=r-z parameter (2018.07.31) <<<<<
               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 800
               if( ixtp .lt. 0 ) goto 996
               if( iytp .lt. 0 ) goto 996
               if( iztp .lt. 0 ) goto 996
CCSE change for mesh=xyz parameter (2017.11.30) <<<<<

            else if( imesh .eq. 4 ) then

               call moddas_allocate_int(
     &                 MAX_NUM_ITREG, idas_itreg_temporary)

               call ttetmeshdc(jsn,jsi,dsin,idsi,ill,ilf,
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                      mtrn,ntrn,
     &                      MAX_NUM_ITREG,idas_itreg_temporary,
     &                      idcrg)

               if( mtrn > MAX_NUM_ITREG ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.ttrack@tallsm2.f'
     &                    //' ?dimension over idas_itreg_temporary?'
     &                    //' mtrn > MAX_NUM_ITREG'
     &                 ,' (mtrn=',mtrn,')'
     &                 ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                  ErrID = 'L:995/R:tdchain/F:tallsm4.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 800

             end if

                  goto 150

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

         else if( ipm .eq. 33 ) then
*           particle name for DCHAIN mode is allowed only 'all'.

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
               if( isubt .eq. 1 )  goto 994   ! kitamura22/03/31

*-----------------------------------------------------------------------

               inpat = inpat + 1

               if( inpat .gt. 6 ) goto 995

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
*           yield special for [T-Yield] is allowed.

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
*           all option of mother for DCHAIN mode is allowed (def=all).

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

                  nsmat = itmat(itnm+1)
                  call moddas_reallocate_int(
     &                    itlmax, itnm+1, imoth, itmat, ismat)
                  if( mmmax .gt. mdas ) goto 950

  141                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                        if( ierr .ne. 0 ) return
                        if( jpn  .eq. 3 ) goto 982

                        if( iskip .ne. 0 ) goto 141

                  ic = i1

               do k = 1, imoth

                  if( ic .gt. i3 ) then

  142                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
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

                     ismat(nsmat-1+k) = icha * 1000

                  else

                     if( isn .gt. 3 ) goto 982

                     read(chlw(icm:icn),'(i5)') masi

                     if( masi .lt. icha ) goto 982
                     if( masi-icha .gt. maxnt ) goto 982

                     ismat(nsmat-1+k) = icha * 1000 + masi

                  end if

                     ic = icd + 1

               end do

            else

               goto 982

            end if

*-----------------------------------------------------------------------
*        nucleus
*-----------------------------------------------------------------------

         else if( ipm .eq. 4 ) then
*           nucleus for DCHAIN mode is not allowed.

            if( chlw(ic:ic+2) .eq. 'all' ) then

               inucl = 0
               imasi = 0

            else if( deqn1( chlw(ic:ic) ) ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               inucl = nint( cvvv )

               if( inucl .le. 0 ) goto 981

                  nnucl = itnuc(itnm+1)
                  call moddas_reallocate_int(
     &                    itlmax, itnm+1, inucl, itnuc, isnuc)
                  if( mmmax .gt. mdas ) goto 950

  143                 call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                        if( ierr .ne. 0 ) return
                        if( jpn  .eq. 3 ) goto 981

                        if( iskip .ne. 0 ) goto 143

                  ic = i1

               do k = 1, inucl

                  if( ic .gt. i3 ) then

  144                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
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

                     isnuc(nnucl-1+k) = icha * 1000

                  else

                     if( isn .gt. 3 ) goto 981

                     read(chlw(icm:icn),'(i5)') masi

                     if( masi .lt. icha ) goto 981
                     if( masi-icha .gt. maxnt ) goto 981
                     if( masi .eq. icha ) goto 975

                     isnuc(nnucl-1+k) = icha * 1000 + masi

                     imasi = imasi + 1

                  end if

                     ic = icd + 1

               end do

            else

               goto 981

            end if

*-----------------------------------------------------------------------
*        material
*-----------------------------------------------------------------------

         else if( ipm .eq. 13 ) then
*           material for DCHAIN mode is allowed.
*           but, this option is for r-z, and xyz mesh.

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

                  nsmte = itmtt(itnm+1)
                  call moddas_reallocate_int(
     &                    itlmax, itnm+1, imate, itmtt, ismte)
                  if( mmmax .gt. mdas ) goto 950

  151                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                        if( ierr .ne. 0 ) return
                        if( jpn  .eq. 3 ) goto 979

                        if( iskip .ne. 0 ) goto 151

                  ic = i1

               do k = 1, imate

                  if( ic .gt. i3 ) then

  152                call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) goto 979

                     if( iskip .ne. 0 ) goto 152

                     ic = i1

                  end if

                     ic = jnumc(chlw,ic,i3)

                     call snum(chlw,ic,i3,ic2,cvvv,ierr)

                     if( ierr .ne. 0 ) goto 979

                     matei = nint( cvvv )

                     ic = ic2

                     ismte(nsmte-1+k) = matei

               end do

            else

               goto 979

            end if

*-----------------------------------------------------------------------
*        axis
*-----------------------------------------------------------------------

         else if( ipm .eq. 5 ) then
*           axis is modified for DCHAIN mode.
*           but, now is automatically defined 'dchain' for [T-Yield].
               inaxi = 1
               iaxis(inaxi) = 13

  600       inaxdc = inaxdc + 1

            if( inaxdc .gt. 6 ) goto 991

            if( chlw(ic:ic+2) .eq. 'all' ) then
               iaxdc(inaxdc) = 14
               ic =jnumc(chlw,ic+4,icl)
            else if( chlw(ic:ic+2) .eq. 'cla' ) then
               iaxdc(inaxdc) = 15
               ic =jnumc(chlw,ic+4,icl)
            else if( chlw(ic:ic+2) .eq. 'gsd' ) then
               iaxdc(inaxdc) = 16
               ic =jnumc(chlw,ic+4,icl)
            else if( chlw(ic:ic+2) .eq. 'gso' ) then
               iaxdc(inaxdc) = 17
               ic =jnumc(chlw,ic+4,icl)
            else if( chlw(ic:ic+2) .eq. 'tba' ) then
               iaxdc(inaxdc) = 18
               ic =jnumc(chlw,ic+4,icl)
            else if( chlw(ic:ic+2) .eq. 'tbh' ) then
               iaxdc(inaxdc) = 19
               ic =jnumc(chlw,ic+4,icl)
            else if( chlw(ic:ic+2) .eq. 'yld' ) then
               iaxdc(inaxdc) = 20
               ic =jnumc(chlw,ic+4,icl)
            else
               goto 992
            end if

               if( ic .le. icl ) goto 600

*        check 'all' option
            do i = 1, inaxdc
            if ( iaxdc(i) .eq. 14 ) then
               inaxdc = 6
               do ia = 1, inaxdc
                  iaxdc(ia) = 14 + ia
               end do
            end if
            end do

*        check 'gsd' and 'gso' option
            do i = 1, inaxdc
            if ( iaxdc(i) .eq. 16 .or. iaxdc(i) .eq. 17 ) then
               imode = 2
               if ( iaxdc(i) .eq. 16 ) then
                  igsdef = 1
               else
                  igsorg = 1
               end if
            end if
            end do

*        check 'cla' option
            do i = 1, inaxdc
            if ( iaxdc(i) .eq. 15) then
               iprtb1 = 1
            end if
            end do

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        file name
*-----------------------------------------------------------------------

         else if( ipm .eq. 6 ) then
*           file name is modified for DCHAIN mode.

  700          infil = infil + 1

               if( infil .gt. 2 ) goto 990

               icf = min( inumc(chlw,ic,icl,' ') - 1, icl )

               lfile(infil) = icf - ic + 1
               ifile(infil)(1:icf-ic+1) = chin(ic:icf)

               jjj = lfile(infil) - 4
               if( jjj .gt. 0 ) then
                  if( ifile(infil)(jjj:jjj+4)
     &                   .eq. '.dyld' .or.
     &                ifile(infil)(jjj:jjj+4)
     &                   .eq. '.dtrk' .or.
     &                ifile(infil)(jjj:jjj+4)
     &                   .eq. '.dout' ) then
                         goto 963
                  end if
               end if

               ic = jnumc(chlw,icf+2,icl)

               if( ic .le. icl ) goto 700

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        restart file name
*-----------------------------------------------------------------------
         else if( ipm .eq. 36 ) then

               call found_param_resfile(itnm,
     &                                  irfflg,ic,icl,chlw,chin,
     &                                  lrfile, irfile)

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

         else if( ipm .eq. 70 ) then

            itsanf = itsanf + 1

            call read_sangel(jsn,jsi,dsin,idsi,ill,ilf,
     &                       jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr
     &                      ,ic ,icl ,'[t-dchain]')
            if( ierr .ne. 0 ) return

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
*           output for DCHAIN mode is allowed only 'product'.

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

                  if( ierr .ne. 0 ) return
                  if( jpn  .eq. 3 ) goto 800
                  if( ixtp .lt. 0 ) goto 996
                  if( iytp .lt. 0 ) goto 996
                  if( iztp .lt. 0 ) goto 996

                  goto 150

*-----------------------------------------------------------------------
*        transform
*-----------------------------------------------------------------------

         else if( ipm .eq. 31 .or. ipm .eq. 32 ) then
*           transform for DCHAIN mode is allowed.
*           but, this option is for r-z, and xyz mesh.

                  if( ipm .eq. 32 ) ktrs = 1

               call ttrans(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     ic,ktrs,igkst,idtt,vtrs)

               if( ierr .ne. 0 ) return

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
*           unit is variable for axis of DCHAIN mode.

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
*           info for DCHAIN mode is allowed.

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
*           factor for DCHAIN mode is allowed.

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               rfact = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        volmat
*-----------------------------------------------------------------------

         else if( ipm .eq. 21 ) then
*           volmat for DCHAIN mode is not allowed.
*           because, this option is for xyz mesh.

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
*        stdcut
*-----------------------------------------------------------------------

         else if( ipm .eq. 69 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               stdcut = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        ndata
*-----------------------------------------------------------------------

         else if( ipm .eq. 19 ) then
*           ndata for DCHAIN mode is allowed.

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ndata = nint( cvvv )

               if( ndata .lt. 0 .or. ndata .gt. 3 ) goto 978

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        The other DCHAIN option (control) ipm = 37 to 67
*-----------------------------------------------------------------------
*        initialization
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*        imode     in Card[2]-1
*-----------------------------------------------------------------------

         else if( ipm .eq. 37 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 922

               imode = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        jmode     in Card[2]-2
*-----------------------------------------------------------------------

         else if( ipm .eq. 38 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 923

               jmode = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        itstep    in Card[3]-1
*-----------------------------------------------------------------------

         else if( ipm .eq. 39 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 924

               idcstp = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        itout     in Card[3]-2
*-----------------------------------------------------------------------

         else if( ipm .eq. 40 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 925

               idcout = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        idivs     in Card[3]-3
*-----------------------------------------------------------------------

         else if( ipm .eq. 41 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 926

               idivs = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        iregon    in Card[3]-4
*-----------------------------------------------------------------------

         else if( ipm .eq. 42 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 927

               iregon = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        inmtcf    in Card[3]-5
*-----------------------------------------------------------------------

         else if( ipm .eq. 43 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 928

               inmtcf = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        ichain    in Card[3]-6
*-----------------------------------------------------------------------

         else if( ipm .eq. 44 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 929

               ichain = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        itdecs    in Card[3]-7
*-----------------------------------------------------------------------

         else if( ipm .eq. 45 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 930

               itdecs = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        itdecn    in Card[3]-8
*-----------------------------------------------------------------------

         else if( ipm .eq. 46 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 931

               itdecn = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        isomtr    in Card[3]-9
*-----------------------------------------------------------------------

         else if( ipm .eq. 47 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 932

               isomtr = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        ifisyd    in Card[3]-10
*-----------------------------------------------------------------------

         else if( ipm .eq. 48 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 933

               ifisyd = nint( cvvv )
               irdfisyd = 1

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        ifisye    in Card[3]-11
*-----------------------------------------------------------------------

         else if( ipm .eq. 49 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 934

               ifisye = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        iyild     in Card[4]-1
*-----------------------------------------------------------------------

         else if( ipm .eq. 50 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 935

               iyild = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        iggrp     in Card[4]-2
*-----------------------------------------------------------------------

         else if( ipm .eq. 51 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 936

               iggrp = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        ibetap    in Card[4]-3
*-----------------------------------------------------------------------

         else if( ipm .eq. 52 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 937

               ibetap = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        acmin     in Card[4]-4
*-----------------------------------------------------------------------

         else if( ipm .eq. 53 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 938

               acmin = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        istabl    in Card[4]-5
*-----------------------------------------------------------------------

         else if( ipm .eq. 54 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 939

               istabl = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        igsdef    in Card[4]-6
*-----------------------------------------------------------------------

         else if( ipm .eq. 55 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 940

               igsdef = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        iprtb1    in Card[4]-7
*-----------------------------------------------------------------------

         else if( ipm .eq. 56 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 941

               iprtb1 = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        iprtb2    in Card[4]-8
*-----------------------------------------------------------------------

         else if( ipm .eq. 57 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 942

               iprtb2 = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        rprtb2    in Card[4]-9
*-----------------------------------------------------------------------

         else if( ipm .eq. 58 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 943

               rprtb2 = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        iprtb3    in Card[4]-10
*-----------------------------------------------------------------------

         else if( ipm .eq. 59 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 944

               iprtb3 = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        igsorg    in Card[4]-11
*-----------------------------------------------------------------------

         else if( ipm .eq. 60 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 945

               igsorg = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        amp       in Card[5]-1
*-----------------------------------------------------------------------

         else if( ipm .eq. 61 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 946

               amp = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        ebeam     in Card[5]-2
*-----------------------------------------------------------------------

         else if( ipm .eq. 62 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 947

               ebeam = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        prodnp    in Card[5]-3
*-----------------------------------------------------------------------

         else if( ipm .eq. 63 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 948

               prodnp = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        tbin and beampw    in Card[6]  irra. / cooling info.
*-----------------------------------------------------------------------

         else if( ipm .eq. 64 ) then


               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 952

               itetp = nint( cvvv )

               call gettev(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     1,itetp,temsh,itemsh)

               if( ierr .ne. 0 ) goto 952
               if ( itetp .ne. idcstp ) idcstp = itetp

               if( ierr .ne. 0 ) return
               goto 140

*-----------------------------------------------------------------------
*        tmin    in Card[7]  output time
*-----------------------------------------------------------------------

         else if( ipm .eq. 65 ) then


               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 953

               iottp = nint( cvvv )

               call gettev(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     2,iottp,temsh,itemsh)

               if( ierr .ne. 0 ) goto 953
               if ( iottp .ne. idcout ) idcout = iottp

               if( ierr .ne. 0 ) return
               goto 140

*-----------------------------------------------------------------------
*        ancza or hancza    in Card[8]  pick up nuclide
*-----------------------------------------------------------------------

         else if( ipm .eq. 66 ) then


         m_err = 'ac-list is disabled function for NEW DCHAIN. '//tnamed
         ErrCha = ''
         ErrID = 'L:2472/R:tdchain/F:tallsm4.f'
           goto 999


*---
*---
! T.Sato 2015/8/25, delete tgnzas & ctgnnds from parameter list
               goto 140

*-----------------------------------------------------------------------
*        hhnmtcf    in Card[9]  [T-Yield] file name
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*        hregcmm    in Card[10] and [10a]  special calc. info.
*-----------------------------------------------------------------------

         else if( ipm .eq. 67 ) then

                   itarget = 0
               call onum(chlw,ic,icl,cvvv,ierr)
               if( ierr .ne. 0 ) return
               itarget = nint( cvvv )
               if (itarget .eq. 0) goto 911

            klst = 0

  163       call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &           jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 800

               if( chlw(i1:i1+2) .eq. 'non' .or.
     &             chlw(i1:i1+2) .eq. 'reg' .or.
     &             chlw(i1:i1+2) .eq. 'vol' ) goto 163
               if( iskip .ne. 0 ) goto 163

               ic = i1

*           for number of target
164          call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 955
               klst = klst + 1

               if(klst.gt.maxreg)call reallocate_talldc4(klst) !FURUTA20200522

               itdcreg=itdcreg+1             !FURUTA20200522
               if(itdcreg.gt.maxtdcreg)call reallocate_talldc2(itdcreg) !FURUTA20200522
               ireg2itm(itdcreg)=itdcregitm  !FURUTA20200522
               jreg2itm(klst)=iregitm        !FURUTA20200522


               ic = jnumc(chlw,ic2,i3)
               if( ic .gt. i3 ) goto 955

*           for cell No.
               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 955

                  itgcell(klst) = nint(cvvv)

               ic = jnumc(chlw,ic2,i3)
               if( ic .gt. i3 ) goto 955

*           for volume
               call onum(chlw,ic,i3,cvvv,ierr)
               if( ierr .ne. 0 ) goto 955

                  rtgvol(klst) = cvvv

                  if( rtgvol(klst) .gt. 0.0d0 ) itgvoll = 1


  165       call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &           jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 955

               if( iskip .ne. 0 ) goto 165

                  ic = i1

               if( chlw(ic:ic+6) .eq. 'tg-list' ) then
                  ic = inumc(chlw,ic,i3,'=') + 1
                  ic = jnumc(chlw,ic,i3)

                  if( ic .gt. i3 ) goto 955

               else
                  goto 955
               end if

*           for number of material
                     call onum(chlw,ic,i3,cvvv,ierr)
                        if( ierr .ne. 0 ) goto 955

                        itglist(klst) = nint( cvvv )
                        itgnum2 = itglist(klst)
*---
*---
*           for material and ratio
! T.Sato 2015/8/25, delete tgnzas & ctgnnds from parameter list
               call getndc(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     2,itgnum2,klst,hancza,ancza)

                 if( ierr .ne. 0 ) return
                 if( jpn  .eq. 3 ) goto 955

  911 continue
* conversion

  172    call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &        jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

            if( ierr .ne. 0 ) return
            if( jpn  .eq. 3 ) goto 800

            if( iskip .ne. 0 ) goto 172

               ic = i1
            if( chlw(ic:ic) .ne. '[' .and. chlw(ic:ic) .lt. 'a' ) then
               goto 164
            else
               goto 150
            end if

*-----------------------------------------------------------------------
*        version    DCHAIN-SP2001 or DCHAIN-SP2015
*-----------------------------------------------------------------------

         else if( ipm .eq. 68 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 957

               iversion = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

CCSE added for dchain parameter (2018.07.31) >>>>>
*-----------------------------------------------------------------------
*         mat number for DCHAIN
*-----------------------------------------------------------------------

         else if( ipm .eq. 71 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 957

               mtscore = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

CCSE added for dchain parameter (2018.07.31) <<<<<

c H.Ratliff START 2019.07.30, add uncertainty propagation, library selection,
c                             and decay chain info file output options

*-----------------------------------------------------------------------
*        iertdcho    in Card[9]-5
*-----------------------------------------------------------------------

         else if( ipm .eq. 72 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 821

               iertdcho = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200


*-----------------------------------------------------------------------
*        inxslib    in Card[3b]-1
*-----------------------------------------------------------------------

         else if( ipm .eq. 73 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 822

               inxslib = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200


*-----------------------------------------------------------------------
*        idcylib    in Card[3b]-2
*-----------------------------------------------------------------------

         else if( ipm .eq. 74 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 823

               idcylib = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200


*-----------------------------------------------------------------------
*        infylib    in Card[3b]-3
*-----------------------------------------------------------------------

         else if( ipm .eq. 75 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 824

               infylib = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        hnxslib    in Card[3b]-4
*-----------------------------------------------------------------------

         else if( ipm .eq. 76 ) then

               ict = min( ic + 79, i2 )

               hnxslib = chin(ic:ict)

               hxsll = ict - ic + 1


*-----------------------------------------------------------------------
*        hdcylib    in Card[3b]-5
*-----------------------------------------------------------------------

         else if( ipm .eq. 77 ) then

               ict = min( ic + 79, i2 )

               hdcylib = chin(ic:ict)

               hdcll = ict - ic + 1


*-----------------------------------------------------------------------
*        hnfylib    in Card[3b]-6
*-----------------------------------------------------------------------

         else if( ipm .eq. 78 ) then

               ict = min( ic + 79, i2 )

               hnfylib = chin(ic:ict)

               hfsll = ict - ic + 1



*-----------------------------------------------------------------------
*        iwrtchn    in Card[4]-12
*-----------------------------------------------------------------------

         else if( ipm .eq. 79 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 828

               iwrtchn = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200


*-----------------------------------------------------------------------
*        chrlvth    in Card[4]-13
*-----------------------------------------------------------------------

         else if( ipm .eq. 80 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 829

               chrlvth = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200


*-----------------------------------------------------------------------
*        iwrchdt    in Card[4]-14
*-----------------------------------------------------------------------

         else if( ipm .eq. 81 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 2830

               iwrchdt = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200


*-----------------------------------------------------------------------
*        iwrchss    in Card[4]-15
*-----------------------------------------------------------------------

         else if( ipm .eq. 82 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 831

               iwrchss = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200


c H.Ratliff START 2020.04.17

*-----------------------------------------------------------------------
*        ixsrall    in Card[3]-7
*-----------------------------------------------------------------------

         else if( ipm .eq. 83 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 832

               ixsrall = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        idosecf    in Card[4]-16
*-----------------------------------------------------------------------

         else if( ipm .eq. 84 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 833

               idosecf = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        ipltmode    in Card[9d]-1
*-----------------------------------------------------------------------

         else if( ipm .eq. 85 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 834

               ipltmode = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        ipltaxis    in Card[9d]-2
*-----------------------------------------------------------------------

         else if( ipm .eq. 86 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 835

               ipltaxis = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200


c H.Ratliff START 2020.05.29
*-----------------------------------------------------------------------
*        foamout    in Card[9d]-1
*-----------------------------------------------------------------------

         else if( ipm .eq. 87 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 836

               foamout = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        foamvals    in Card[9d]-2
*-----------------------------------------------------------------------

         else if( ipm .eq. 88 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 837

               foamvals = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
c H.Ratliff START 2020.06.05
*-----------------------------------------------------------------------
*        iredufmt    in Card[8]-7
*-----------------------------------------------------------------------

         else if( ipm .eq. 89 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 838

               iredufmt0 = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
c H.Ratliff START 2020.06.17
*-----------------------------------------------------------------------
*        irdonce    in Card[8]-7
*-----------------------------------------------------------------------

         else if( ipm .eq. 90 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 839

               irdonce = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
c H.Ratliff START 2020.11.16
*-----------------------------------------------------------------------
*        imtcard    in Card[4]-17
*-----------------------------------------------------------------------

         else if( ipm .eq. 91 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 840

               imtcard = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
*-----------------------------------------------------------------------
*        imtcnum    in Card[4]-18
*-----------------------------------------------------------------------

         else if( ipm .eq. 92 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 841

               imtcnum = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
*-----------------------------------------------------------------------
*        imtcmeta    in Card[4]-19
*-----------------------------------------------------------------------

         else if( ipm .eq. 93 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 842

               imtcmeta = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
*-----------------------------------------------------------------------
*        thmatnd    in Card[4]-20
*-----------------------------------------------------------------------

         else if( ipm .eq. 94 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 843

               thmatnd = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
c H.Ratliff START 2020.12.18
*-----------------------------------------------------------------------
*        idosunit    in Card[4]-16b
*-----------------------------------------------------------------------

         else if( ipm .eq. 95 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 844

               idosunit = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
*-----------------------------------------------------------------------
*        iangbpwr    in Card[10]-2
*-----------------------------------------------------------------------

         else if( ipm .eq. 96 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 845

               iangbpwr = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
c H.Ratliff START 2021.02.17
*-----------------------------------------------------------------------
*        iphtout    in Card[10]-1
*-----------------------------------------------------------------------

         else if( ipm .eq. 97 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 846

               iphtout = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
c H.Ratliff START 2021.04.07
*-----------------------------------------------------------------------
*        isfylib    in Card[3b]-7
*-----------------------------------------------------------------------

         else if( ipm .eq. 98 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 847

               isfylib = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
*-----------------------------------------------------------------------
*        radlev    in Card[10]-3
*-----------------------------------------------------------------------

         else if( ipm .eq. 99 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 848

               radlev = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
*-----------------------------------------------------------------------
*        ergthfa    in Card[3b]-8
*-----------------------------------------------------------------------

         else if( ipm .eq. 100 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 849

               ergthfa = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
*-----------------------------------------------------------------------
*        ergfafu    in Card[3b]-9
*-----------------------------------------------------------------------

         else if( ipm .eq. 101 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 1850

               ergfafu = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
*-----------------------------------------------------------------------
*        infyerg    in Card[3b]-10
*-----------------------------------------------------------------------

         else if( ipm .eq. 102 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 851

               infyerg = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
*-----------------------------------------------------------------------
*        hsfylib    in Card[3b]-11
*-----------------------------------------------------------------------

         else if( ipm .eq. 103 ) then

               ict = min( ic + 79, i2 )

               hsfylib = chin(ic:ict)

               hssll = ict - ic + 1

c H.Ratliff START 2021.04.26
*-----------------------------------------------------------------------
*        ilchain    in Card[3]-12
*-----------------------------------------------------------------------

         else if( ipm .eq. 104 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 852

               ilchain = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        mxnuclei
*-----------------------------------------------------------------------
         else if( ipm .eq. 105 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               imxnuc = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        itriton
*-----------------------------------------------------------------------
         else if( ipm .eq. 106 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 853

               itriton = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        iaonucl
*-----------------------------------------------------------------------
         else if( ipm .eq. 107 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 854

               iaonucl = nint( cvvv )

               if(iaonucl.ge.14) goto 854

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        aonucl
*-----------------------------------------------------------------------

         else if( ipm .eq. 108 ) then

               ict = min( ic + 199, i2 )

               aonucl = chin(ic:ict)

               laonucl = ict - ic + 1

*-----------------------------------------------------------------------
*        aoreg
*-----------------------------------------------------------------------

         else if( ipm .eq. 109 ) then

               ict = min( ic + 199, i2 )

               aoreg = chin(ic:ict)

               laoreg = ict - ic + 1

*-----------------------------------------------------------------------

         end if

            goto 140

*-----------------------------------------------------------------------
*     check
*-----------------------------------------------------------------------

  800 continue

      if( ischn(7) .eq. 0 ) then
*        title
         if( imesh .eq. 1 ) then
*        mesh for DCHAIN mode is allowed only region mesh.

            title = '[t-yield] in region mesh'
            titll = 24

         else if( imesh .eq. 2 ) then

            title = '[t-yield] in r-z mesh'
            titll = 21

         else if( imesh .eq. 3 ) then

            title = '[t-yield] in xyz mesh'
            titll = 21

         else if( imesh .eq. 4 ) then
            title = '[t-yield] in tet mesh'
            titll = 21

         end if

      end if

      if( ischn(1) .eq. 0 ) then

         m_err = 'Mesh is not defined in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:3272/R:tdchain/F:tallsm4.f'
         goto 999

      else if( inaxi .eq. 0 ) then

         m_err = 'Axis is not defined in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:3279/R:tdchain/F:tallsm4.f'
         goto 999

      else if( infil .eq. 0 ) then

         m_err = 'File is not defined in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:3286/R:tdchain/F:tallsm4.f'
         goto 999

*        axis is modified, so This option is deleted.
      else if( inaxi .ne. infil ) then

        m_err = 'Number of Axis and File is different in tally '//tnamed
        ErrCha = ''
        ErrID = 'L:3294/R:tdchain/F:tallsm4.f'
         goto 999

      else if( itsanf > 0 ) then

         m_err ='sangel parameter is invalid in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:3301/R:tdchain/F:tallsm4.f'
         goto 999

      end if

      do i = 1, inaxi

         if( lrshow .eq. 0 ) then

            if( ( iaxis(i) .eq. 9 .or. iaxis(i) .eq. 10 .or.
     &            iaxis(i) .eq. 11 ) .and. imesh .ne. 3 ) then

               m_err = 'Axis is xy, yz, or zx, '//
     &                 'but mesh is not xyz in tally '//tnamed
               ErrCha = ''
               ErrID = 'L:3316/R:tdchain/F:tallsm4.f'
               goto 999

            end if

         end if

         if( iaxis(i) .eq. 2 .and. imesh .ne. 1 ) then

            m_err = 'Axis is reg but mesh is not reg in tally '//tnamed
            ErrCha = ''
            ErrID = 'L:3327/R:tdchain/F:tallsm4.f'
            goto 999

CCSE chg for mesh=r-z (2018.07.31) >>>>>
         else if( iaxis(i) .eq. 13 .and. imesh .eq. 2 ) then

            m_err = 'Axis is dchain but mesh is r-z in tally '
     &              //tnamed
            ErrCha = ''
            ErrID = 'L:3336/R:tdchain/F:tallsm4.f'
            goto 999
CCSE chg for mesh=r-z (2018.07.31) <<<<<

         else if( iaxis(i) .eq. 13 .and. iunt .ne. 1 ) then

            m_err = 'Axis is dchain but unit is not 1 in tally '
     &              //tnamed
            ErrCha = ''
            ErrID = 'L:3345/R:tdchain/F:tallsm4.f'
            goto 999

         else if( iaxis(i) .eq. 3 .and. imesh .ne. 3 ) then

            m_err = 'Axis is x but mesh is not xyz in tally '//tnamed
            ErrCha = ''
            ErrID = 'L:3352/R:tdchain/F:tallsm4.f'
            goto 999

         else if( iaxis(i) .eq. 4 .and. imesh .ne. 3 ) then

            m_err = 'Axis is y but mesh is not xyz in tally '//tnamed
            ErrCha = ''
            ErrID = 'L:3359/R:tdchain/F:tallsm4.f'
            goto 999

         else if( iaxis(i) .eq. 5 .and. imesh .eq. 1 ) then

            m_err = 'Axis is z but mesh is reg in tally '//tnamed
            ErrCha = ''
            ErrID = 'L:3366/R:tdchain/F:tallsm4.f'
            goto 999

         else if( iaxis(i) .eq. 6 .and. imesh .ne. 2 ) then

            m_err = 'Axis is r but mesh is not r-z in tally '//tnamed
            ErrCha = ''
            ErrID = 'L:3373/R:tdchain/F:tallsm4.f'
            goto 999

         else if( ( iaxis(i) .eq. 7 .or. iaxis(i) .eq. 8 ) .and.
     &              inucl .gt. 0 ) then

            m_err = 'Axis is charge or chart'//
     &              ' but nucleus is restricted in tally '//tnamed
            ErrCha = ''
            ErrID = 'L:3382/R:tdchain/F:tallsm4.f'
            goto 999

         else if( ( iaxis(i) .eq. 1 .or. iaxis(i) .eq. 13 ) .and.
     &              imasi .gt. 0 ) then

            m_err = 'Axis is dchain or mass,'//
     &              ' but isotope is restricted in tally '//tnamed
            ErrCha = ''
            ErrID = 'L:3391/R:tdchain/F:tallsm4.f'
            goto 999

         else if( iaxis(i) .eq. 12 .and. imesh .ne. 2 ) then

            m_err = 'Axis is rz,'//
     &              ' but mesh is not r-z in tally '//tnamed
            ErrCha = ''
            ErrID = 'L:3399/R:tdchain/F:tallsm4.f'
            goto 999

*    &           ( irtp .ne. 2 .and. irtp .ne. 4 ) ) then
*
*    &              ' but r mesh type is not 2 or 4 in tally '//tnamed
         else if( iaxis(i) .eq. 12 .and.
     &           ( iztp .ne. 2 .and. iztp .ne. 4 ) ) then

            m_err = 'Axis is rz,'//
     &              ' but z mesh type is not 2 or 4 in tally '//tnamed
            ErrCha = ''
            ErrID = 'L:3411/R:tdchain/F:tallsm4.f'
            goto 999

         else if( ( iaxis(i) .eq. 9 .or. iaxis(i) .eq. 11 ) .and.
     &           ( ixtp .ne. 2 .and. ixtp .ne. 4 ) ) then

            m_err = 'Axis is xy or zx,'//
     &              ' but x mesh type is not 2 or 4 in tally '//tnamed
            ErrCha = ''
            ErrID = 'L:3420/R:tdchain/F:tallsm4.f'
            goto 999

         else if( ( iaxis(i) .eq. 10 .or. iaxis(i) .eq. 11 ) .and.
     &           ( iztp .ne. 2 .and. iztp .ne. 4 ) ) then

            m_err = 'Axis is yz or zx,'//
     &              ' but z mesh type is not 2 or 4 in tally '//tnamed
            ErrCha = ''
            ErrID = 'L:3429/R:tdchain/F:tallsm4.f'
            goto 999

         else if( ( iaxis(i) .eq. 10 .or. iaxis(i) .eq. 9 ) .and.
     &           ( iytp .ne. 2 .and. iytp .ne. 4 ) ) then

            m_err = 'Axis is xy or yz,'//
     &              ' but y mesh type is not 2 or 4 in tally '//tnamed
            ErrCha = ''
            ErrID = 'L:3438/R:tdchain/F:tallsm4.f'
            goto 999

         else if( iaxis(i) .eq. 1 .and. inucl .gt. 1 ) then

            do k = 1, inucl
               do j = k + 1, inucl

                  if( isnuc(nnucl-1+k)/1000 .eq.
     &                isnuc(nnucl-1+j)/1000 ) then

                     m_err = 'Axis is mass,'//
     &               ' and there is the same nucleus in tally '//tnamed
                     ErrCha = ''
                     ErrID = 'L:3452/R:tdchain/F:tallsm4.f'
                     goto 999

                  end if

               end do
            end do

         end if

      end do

*     Card[2]
      if ( imode .lt. 0 .or. imode .gt. 2 ) then
         m_err = 'imode is required 0, 1 or 2'
         ErrCha = ''
         ErrID = 'L:3468/R:tdchain/F:tallsm4.f'
         goto 999
      end if
      if ( jmode .lt. -1 .or. jmode .gt. 2) then
         m_err = 'jmode is required -1, 0, 1 or 2'
         ErrCha = ''
         ErrID = 'L:3474/R:tdchain/F:tallsm4.f'
         goto 999
      end if

*     Card[3]
      if ( idcstp .gt. maxitm .or. idcout .gt. maxitm ) then
         m_err = 'the number of irradiation/cooling time bins was'//
     &           ' exceeded the limit.'
         ErrCha = ''
         ErrID = 'L:3483/R:tdchain/F:tallsm4.f'
         goto 999
      end if
      if ( iregon .gt. maxreg ) then !FURUTA20200522
         m_err = 'the number of processing many regions exceeded'//
     &           ' the limit.'
         ErrCha = ''
         ErrID = 'L:3490/R:tdchain/F:tallsm4.f'
         goto 999
      end if
      if ( inmtcf .gt. maxnmtcf ) then
         m_err = 'the number of [T-YIELD] files exceeded the'//
     &           ' limit.'
         ErrCha = ''
         ErrID = 'L:3497/R:tdchain/F:tallsm4.f'
         goto 999
      end if
      if ( jmode .ne. 0 .and. jmode .ne. 2 .and.
     &     itdecs .eq. 1 ) itdecs = 0
      if ( jmode .ne. 1 .and. jmode .ne. 2 .and.
     &     itdecn .eq. 1 ) itdecn = 0

*     Card[4]
      if ( iregon .le. 1 ) iprtb3 = 0
      if ( imode .le. 1 ) then
         iggrp = 0
         ibetap = 0
         igsdef = 0
         igsorg = 0
      end if
      if ( iggrp .eq. 0 ) then
         ibetap = 0
         igsdef = 0
         igsorg = 0
      end if

*     Card[5]
      if (( jmode .eq. 0 .or. jmode .eq. 2 ) .and.
     &      amp. eq. zero ) then
         m_err = 'beam current of spallation was less than 0.'
         ErrCha = ''
         ErrID = 'L:3524/R:tdchain/F:tallsm4.f'
         goto 999
c Ratliff 2021.01.27 comment out incorrect jmode statement
      else if ( jmode .eq. -1 ) then
         amp = 0
         ebeam = 0
      end if
      if ( prodnp .le. zero) prodnp = 1
      if ( jmode .le. 0) prodnp = 0

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

                  itnm = itnm + 1

*        for [T-Yield]: 3, [T-Track]: 1 and [T-Dchain]: 16.
                  if( itnm+2 .gt. itlmax ) goto 986

                  ital( itnm ) = 3
                  ital( itnm + 2 ) = 16

                  itsmn(itnm) = 0
                  itstm(itnm) = 0

                  itxnm(itnm) = 0
                  itynm(itnm) = 0
                  itznm(itnm) = 0

* ---

* ---
                  itsmn(itnm+2) = 0
                  itstm(itnm+2) = 0

                  itxnm(itnm+2) = 0
                  itynm(itnm+2) = 0
                  itznm(itnm+2) = 0

                  itnm2tdc(itnm)=itdc !FURUTA20200522

*-----------------------------------------------------------------------
*     create DCTrack and DCYield tallies
*-----------------------------------------------------------------------

            iodw = 201
            open( iodw, status = 'scratch', action = 'readwrite' )
            if ( iversion .eq. 0 ) then
              write(iodw, 503)
            else if ( iversion .eq. 2 ) then  ! user defined energy bin, T.Sato 2024/11/19
             iodc = 22  ! temporary used
             open(iodc,file=chfn(1)(1:ilfn(1))//
     &       '/data/dchain_EnGroup.dat', status='old')
             do
              read(iodc, '(A)', end=504) line
              write(iodw, '(A)') line
             end do
 504         close(iodc)
            else
              write(iodw, 501)
            end if

  501 format('[ T-TRACK ] off'
     &/'     part =  neutron'
     &/'   e-type =    1'
     &/'       ne = 1968'
     &/'   1.00001E-11 3.00000E-09 5.00000E-09 6.90000E-09 1.00000E-08'
     &/'   1.50000E-08 2.00000E-08 2.50000E-08 3.00000E-08 3.50000E-08'
     &/'   4.20000E-08 5.00000E-08 5.80000E-08 6.70000E-08 7.70000E-08'
     &/'   8.00000E-08 9.50000E-08 1.00000E-07 1.15000E-07 1.34000E-07'
     &/'   1.40000E-07 1.46370E-07 1.53030E-07 1.60000E-07 1.69710E-07'
     &/'   1.80000E-07 1.89000E-07 1.98810E-07 2.09140E-07 2.20000E-07'
     &/'   2.33580E-07 2.48000E-07 2.63510E-07 2.80000E-07 3.00000E-07'
     &/'   3.14500E-07 3.20000E-07 3.34660E-07 3.50000E-07 3.69930E-07'
     &/'   3.91000E-07 4.00000E-07 4.13990E-07 4.33000E-07 4.49680E-07'
     &/'   4.67010E-07 4.85000E-07 5.00000E-07 5.19620E-07 5.31580E-07'
     &/'   5.40000E-07 5.66960E-07 5.95280E-07 6.25000E-07 6.53150E-07'
     &/'   6.82560E-07 7.05000E-07 7.41550E-07 7.80000E-07 7.90000E-07'
     &/'   8.19450E-07 8.50000E-07 8.60000E-07 8.76425E-07 9.10000E-07'
     &/'   9.30000E-07 9.50000E-07 9.72000E-07 9.86000E-07 9.96000E-07'
     &/'   1.02000E-06 1.03500E-06 1.04500E-06 1.07100E-06 1.08000E-06'
     &/'   1.09700E-06 1.11000E-06 1.12300E-06 1.15000E-06 1.17000E-06'
     &/'   1.20206E-06 1.23500E-06 1.26708E-06 1.30000E-06 1.33750E-06'
     &/'   1.37000E-06 1.40456E-06 1.44000E-06 1.47500E-06 1.50000E-06'
     &/'   1.54434E-06 1.59000E-06 1.62951E-06 1.67000E-06 1.71197E-06'
     &/'   1.75500E-06 1.79700E-06 1.84000E-06 1.85539E-06 1.88446E-06'
     &/'   1.93000E-06 1.97449E-06 2.02000E-06 2.05961E-06 2.10000E-06'
     &/'   2.13000E-06 2.18531E-06 2.24205E-06 2.30027E-06 2.36000E-06'
     &/'   2.38237E-06 2.42171E-06 2.48503E-06 2.55000E-06 2.60000E-06'
     &/'   2.65932E-06 2.72000E-06 2.76792E-06 2.83799E-06 2.90983E-06'
     &/'   2.98349E-06 3.05902E-06 3.13733E-06 3.21763E-06 3.30000E-06'
     &/'   3.38075E-06 3.46633E-06 3.55408E-06 3.64405E-06 3.73630E-06'
     &/'   3.83088E-06 3.92786E-06 4.00000E-06 4.12925E-06 4.23378E-06'
     &/'   4.34096E-06 4.45085E-06 4.56353E-06 4.67905E-06 4.79750E-06'
     &/'   4.91895E-06 5.04348E-06 5.08568E-06 5.12824E-06 5.17115E-06'
     &/'   5.21443E-06 5.25806E-06 5.30206E-06 5.34643E-06 5.39117E-06'
     &/'   5.43628E-06 5.48177E-06 5.52765E-06 5.57390E-06 5.62055E-06'
     &/'   5.66758E-06 5.71501E-06 5.76283E-06 5.81106E-06 5.85968E-06'
     &/'   5.90872E-06 5.95816E-06 6.00802E-06 6.05830E-06 6.10900E-06'
     &/'   6.16012E-06 6.21166E-06 6.26365E-06 6.31606E-06 6.36891E-06'
     &/'   6.42221E-06 6.47595E-06 6.53014E-06 6.58479E-06 6.63989E-06'
     &/'   6.69545E-06 6.75148E-06 6.80798E-06 6.86495E-06 6.92240E-06'
     &/'   6.98033E-06 7.03874E-06 7.09764E-06 7.15703E-06 7.21692E-06'
     &/'   7.27732E-06 7.33821E-06 7.39962E-06 7.46154E-06 7.52398E-06'
     &/'   7.58695E-06 7.65043E-06 7.71445E-06 7.77901E-06 7.84410E-06'
     &/'   7.90975E-06 7.97594E-06 8.04268E-06 8.10998E-06 8.17785E-06'
     &/'   8.24628E-06 8.31529E-06 8.38487E-06 8.45504E-06 8.52579E-06'
     &/'   8.59713E-06 8.66908E-06 8.74162E-06 8.81477E-06 8.88854E-06'
     &/'   8.96292E-06 9.03792E-06 9.11355E-06 9.18981E-06 9.26671E-06'
     &/'   9.34426E-06 9.42245E-06 9.50130E-06 9.58081E-06 9.66099E-06'
     &/'   9.74183E-06 9.82335E-06 9.90555E-06 9.98845E-06 1.00720E-05'
     &/'   1.01563E-05 1.02413E-05 1.03270E-05 1.04134E-05 1.05006E-05'
     &/'   1.05884E-05 1.06770E-05 1.07664E-05 1.08565E-05 1.09473E-05'
     &/'   1.10389E-05 1.11313E-05 1.12245E-05 1.13184E-05 1.14131E-05'
     &/'   1.15086E-05 1.16049E-05 1.17020E-05 1.18000E-05 1.18987E-05'
     &/'   1.19983E-05 1.20987E-05 1.21999E-05 1.23020E-05 1.24049E-05'
     &/'   1.25088E-05 1.26134E-05 1.27190E-05 1.28254E-05 1.29327E-05'
     &/'   1.30410E-05 1.31501E-05 1.32601E-05 1.33711E-05 1.34830E-05'
     &/'   1.35958E-05 1.37096E-05 1.38243E-05 1.39400E-05 1.40566E-05'
     &/'   1.41743E-05 1.42929E-05 1.44125E-05 1.45331E-05 1.46547E-05'
     &/'   1.47774E-05 1.49010E-05 1.50257E-05 1.51514E-05 1.52782E-05'
     &/'   1.54061E-05 1.55350E-05 1.56650E-05 1.57961E-05 1.59283E-05'
     &/'   1.60616E-05 1.61960E-05 1.63315E-05 1.64682E-05 1.66060E-05'
     &/'   1.67449E-05 1.68851E-05 1.70264E-05 1.71688E-05 1.73125E-05'
     &/'   1.74574E-05 1.76035E-05 1.77508E-05 1.78993E-05 1.80491E-05'
     &/'   1.82001E-05 1.83524E-05 1.85060E-05 1.86609E-05 1.88170E-05'
     &/'   1.89745E-05 1.91333E-05 1.92934E-05 1.94548E-05 1.96176E-05'
     &/'   1.97818E-05 1.99473E-05 2.01143E-05 2.02826E-05 2.04523E-05'
     &/'   2.06234E-05 2.07960E-05 2.09701E-05 2.11455E-05 2.13225E-05'
     &/'   2.15009E-05 2.16808E-05 2.18623E-05 2.20452E-05 2.22297E-05'
     &/'   2.24157E-05 2.26033E-05 2.27924E-05 2.29832E-05 2.31755E-05'
     &/'   2.33694E-05 2.35650E-05 2.37622E-05 2.39610E-05 2.41615E-05'
     &/'   2.43637E-05 2.45676E-05 2.47732E-05 2.49805E-05 2.51895E-05'
     &/'   2.54003E-05 2.56129E-05 2.58272E-05 2.60434E-05 2.62613E-05'
     &/'   2.64810E-05 2.67026E-05 2.69261E-05 2.71514E-05 2.73786E-05'
     &/'   2.76077E-05 2.78388E-05 2.80717E-05 2.83066E-05 2.85435E-05'
     &/'   2.87824E-05 2.90232E-05 2.92661E-05 2.95110E-05 2.97579E-05'
     &/'   3.00069E-05 3.02581E-05 3.05113E-05 3.07666E-05 3.10240E-05'
     &/'   3.12837E-05 3.15454E-05 3.18094E-05 3.20756E-05 3.23440E-05'
     &/'   3.26147E-05 3.28876E-05 3.31628E-05 3.34403E-05 3.37201E-05'
     &/'   3.40023E-05 3.42869E-05 3.45738E-05 3.48631E-05 3.51548E-05'
     &/'   3.54490E-05 3.57457E-05 3.60448E-05 3.63464E-05 3.66506E-05'
     &/'   3.69573E-05 3.72665E-05 3.75784E-05 3.78929E-05 3.82099E-05'
     &/'   3.85297E-05 3.88521E-05 3.91772E-05 3.95051E-05 3.98357E-05'
     &/'   4.01690E-05 4.05051E-05 4.08441E-05 4.11859E-05 4.15305E-05'
     &/'   4.18781E-05 4.22285E-05 4.25819E-05 4.29382E-05 4.32975E-05'
     &/'   4.36598E-05 4.40252E-05 4.43936E-05 4.47651E-05 4.51397E-05'
     &/'   4.55174E-05 4.58983E-05 4.62824E-05 4.66697E-05 4.70603E-05'
     &/'   4.74541E-05 4.78512E-05 4.82516E-05 4.86554E-05 4.90625E-05'
     &/'   4.94731E-05 4.98871E-05 5.03046E-05 5.07255E-05 5.11500E-05'
     &/'   5.15780E-05 5.20096E-05 5.24449E-05 5.28837E-05 5.33263E-05'
     &/'   5.37725E-05 5.42225E-05 5.46762E-05 5.51338E-05 5.55951E-05'
     &/'   5.60604E-05 5.65295E-05 5.70025E-05 5.74795E-05 5.79605E-05'
     &/'   5.84456E-05 5.89346E-05 5.94278E-05 5.99251E-05 6.04266E-05'
     &/'   6.09322E-05 6.14421E-05 6.19563E-05 6.24747E-05 6.29975E-05'
     &/'   6.35247E-05 6.40563E-05 6.45923E-05 6.51329E-05 6.56779E-05'
     &/'   6.62275E-05 6.67817E-05 6.73405E-05 6.79040E-05 6.84723E-05'
     &/'   6.90453E-05 6.96231E-05 7.02057E-05 7.07932E-05 7.13856E-05'
     &/'   7.19829E-05 7.25853E-05 7.31927E-05 7.38052E-05 7.44228E-05'
     &/'   7.50456E-05 7.56736E-05 7.63068E-05 7.69454E-05 7.75893E-05'
     &/'   7.82385E-05 7.88932E-05 7.95534E-05 8.02191E-05 8.08904E-05'
     &/'   8.15673E-05 8.22499E-05 8.29382E-05 8.36322E-05 8.43321E-05'
     &/'   8.50378E-05 8.57494E-05 8.64669E-05 8.71905E-05 8.79202E-05'
     &/'   8.86559E-05 8.93978E-05 9.01459E-05 9.09002E-05 9.16609E-05'
     &/'   9.24279E-05 9.32014E-05 9.39813E-05 9.47677E-05 9.55608E-05'
     &/'   9.63604E-05 9.71668E-05 9.79799E-05 9.87998E-05 9.96266E-05'
     &/'   1.00460E-04 1.01301E-04 1.02149E-04 1.03003E-04 1.03865E-04'
     &/'   1.04734E-04 1.05611E-04 1.06495E-04 1.07386E-04 1.08284E-04'
     &/'   1.09191E-04 1.10104E-04 1.11026E-04 1.11955E-04 1.12892E-04'
     &/'   1.13836E-04 1.14789E-04 1.15750E-04 1.16718E-04 1.17695E-04'
     &/'   1.18680E-04 1.19673E-04 1.20674E-04 1.21684E-04 1.22702E-04'
     &/'   1.23729E-04 1.24765E-04 1.25809E-04 1.26862E-04 1.27923E-04'
     &/'   1.28994E-04 1.30073E-04 1.31161E-04 1.32259E-04 1.33366E-04'
     &/'   1.34482E-04 1.35607E-04 1.36742E-04 1.37886E-04 1.39040E-04'
     &/'   1.40204E-04 1.41377E-04 1.42560E-04 1.43753E-04 1.44956E-04'
     &/'   1.46169E-04 1.47392E-04 1.48625E-04 1.49869E-04 1.51123E-04'
     &/'   1.52388E-04 1.53663E-04 1.54949E-04 1.56246E-04 1.57553E-04'
     &/'   1.58871E-04 1.60201E-04 1.61542E-04 1.62893E-04 1.64256E-04'
     &/'   1.65631E-04 1.67017E-04 1.68415E-04 1.69824E-04 1.71245E-04'
     &/'   1.72678E-04 1.74123E-04 1.75580E-04 1.77049E-04 1.78531E-04'
     &/'   1.80025E-04 1.81531E-04 1.83050E-04 1.84582E-04 1.86127E-04'
     &/'   1.87685E-04 1.89255E-04 1.90839E-04 1.92436E-04 1.94046E-04'
     &/'   1.95670E-04 1.97307E-04 1.98958E-04 2.00623E-04 2.02302E-04'
     &/'   2.03995E-04 2.05702E-04 2.07423E-04 2.09159E-04 2.10910E-04'
     &/'   2.12674E-04 2.14454E-04 2.16249E-04 2.18058E-04 2.19883E-04'
     &/'   2.21723E-04 2.23578E-04 2.25449E-04 2.27336E-04 2.29238E-04'
     &/'   2.31157E-04 2.33091E-04 2.35042E-04 2.37008E-04 2.38992E-04'
     &/'   2.40992E-04 2.43008E-04 2.45042E-04 2.47092E-04 2.49160E-04'
     &/'   2.51245E-04 2.53348E-04 2.55468E-04 2.57605E-04 2.59761E-04'
     &/'   2.61935E-04 2.64127E-04 2.66337E-04 2.68566E-04 2.70813E-04'
     &/'   2.73079E-04 2.75365E-04 2.77669E-04 2.79992E-04 2.82335E-04'
     &/'   2.84698E-04 2.87080E-04 2.89483E-04 2.91905E-04 2.94348E-04'
     &/'   2.96811E-04 2.99295E-04 3.01799E-04 3.04325E-04 3.06872E-04'
     &/'   3.09439E-04 3.12029E-04 3.14640E-04 3.17273E-04 3.19928E-04'
     &/'   3.22605E-04 3.25305E-04 3.28027E-04 3.30772E-04 3.33540E-04'
     &/'   3.36331E-04 3.39145E-04 3.41983E-04 3.44845E-04 3.47731E-04'
     &/'   3.50641E-04 3.53575E-04 3.56534E-04 3.59517E-04 3.62526E-04'
     &/'   3.65560E-04 3.68618E-04 3.71703E-04 3.74814E-04 3.77950E-04'
     &/'   3.81113E-04 3.84302E-04 3.87518E-04 3.90761E-04 3.94031E-04'
     &/'   3.97328E-04 4.00653E-04 4.04006E-04 4.07387E-04 4.10796E-04'
     &/'   4.14233E-04 4.17700E-04 4.21195E-04 4.24719E-04 4.28274E-04'
     &/'   4.31858E-04 4.35471E-04 4.39115E-04 4.42790E-04 4.46495E-04'
     &/'   4.50232E-04 4.53999E-04 4.57798E-04 4.61629E-04 4.65492E-04'
     &/'   4.69388E-04 4.73316E-04 4.77276E-04 4.81270E-04 4.85298E-04'
     &/'   4.89359E-04 4.93454E-04 4.97583E-04 5.01747E-04 5.05945E-04'
     &/'   5.10179E-04 5.14449E-04 5.18754E-04 5.23095E-04 5.27472E-04'
     &/'   5.31886E-04 5.36337E-04 5.40825E-04 5.45351E-04 5.49914E-04'
     &/'   5.54516E-04 5.59156E-04 5.63835E-04 5.68554E-04 5.73311E-04'
     &/'   5.78109E-04 5.82947E-04 5.87825E-04 5.92744E-04 5.97704E-04'
     &/'   6.02706E-04 6.07749E-04 6.12835E-04 6.17963E-04 6.23135E-04'
     &/'   6.28349E-04 6.33607E-04 6.38909E-04 6.44256E-04 6.49647E-04'
     &/'   6.55083E-04 6.60565E-04 6.66093E-04 6.71667E-04 6.77287E-04'
     &/'   6.82955E-04 6.88670E-04 6.94433E-04 7.00244E-04 7.06104E-04'
     &/'   7.12013E-04 7.17971E-04 7.23979E-04 7.30037E-04 7.36146E-04'
     &/'   7.42307E-04 7.48518E-04 7.54782E-04 7.61098E-04 7.67467E-04'
     &/'   7.73889E-04 7.80365E-04 7.86896E-04 7.93481E-04 8.00120E-04'
     &/'   8.06816E-04 8.13568E-04 8.20376E-04 8.27241E-04 8.34163E-04'
     &/'   8.41144E-04 8.48182E-04 8.55280E-04 8.62437E-04 8.69654E-04'
     &/'   8.76932E-04 8.84270E-04 8.91670E-04 8.99131E-04 9.06655E-04'
     &/'   9.14242E-04 9.21893E-04 9.29607E-04 9.37386E-04 9.45231E-04'
     &/'   9.53140E-04 9.61116E-04 9.69159E-04 9.77269E-04 9.85447E-04'
     &/'   9.93694E-04 1.00201E-03 1.01039E-03 1.01885E-03 1.02738E-03'
     &/'   1.03597E-03 1.04464E-03 1.05338E-03 1.06220E-03 1.07109E-03'
     &/'   1.08005E-03 1.08909E-03 1.09820E-03 1.10739E-03 1.11666E-03'
     &/'   1.12600E-03 1.13542E-03 1.14493E-03 1.15451E-03 1.16417E-03'
     &/'   1.17391E-03 1.18373E-03 1.19364E-03 1.20363E-03 1.21370E-03'
     &/'   1.22386E-03 1.23410E-03 1.24442E-03 1.25484E-03 1.26534E-03'
     &/'   1.27593E-03 1.28661E-03 1.29737E-03 1.30823E-03 1.31918E-03'
     &/'   1.33021E-03 1.34135E-03 1.35257E-03 1.36389E-03 1.37530E-03'
     &/'   1.38681E-03 1.39842E-03 1.41012E-03 1.42192E-03 1.43382E-03'
     &/'   1.44582E-03 1.45791E-03 1.47012E-03 1.48242E-03 1.49482E-03'
     &/'   1.50733E-03 1.51994E-03 1.53266E-03 1.54549E-03 1.55842E-03'
     &/'   1.57146E-03 1.58461E-03 1.59787E-03 1.61124E-03 1.62473E-03'
     &/'   1.63832E-03 1.65203E-03 1.66586E-03 1.67980E-03 1.69386E-03'
     &/'   1.70803E-03 1.72232E-03 1.73673E-03 1.75127E-03 1.76592E-03'
     &/'   1.78070E-03 1.79560E-03 1.81063E-03 1.82578E-03 1.84106E-03'
     &/'   1.85646E-03 1.87200E-03 1.88766E-03 1.90346E-03 1.91939E-03'
     &/'   1.93545E-03 1.95165E-03 1.96798E-03 1.98445E-03 2.00105E-03'
     &/'   2.01780E-03 2.03468E-03 2.05171E-03 2.06888E-03 2.08619E-03'
     &/'   2.10365E-03 2.12125E-03 2.13900E-03 2.15690E-03 2.17495E-03'
     &/'   2.19315E-03 2.21151E-03 2.23001E-03 2.24867E-03 2.26749E-03'
     &/'   2.28647E-03 2.30560E-03 2.32489E-03 2.34435E-03 2.36397E-03'
     &/'   2.38375E-03 2.40370E-03 2.42381E-03 2.44409E-03 2.46455E-03'
     &/'   2.48517E-03 2.50597E-03 2.52693E-03 2.54808E-03 2.56940E-03'
     &/'   2.59090E-03 2.61259E-03 2.63445E-03 2.65649E-03 2.67872E-03'
     &/'   2.70114E-03 2.72374E-03 2.74654E-03 2.76952E-03 2.79270E-03'
     &/'   2.81607E-03 2.83963E-03 2.86339E-03 2.88735E-03 2.91151E-03'
     &/'   2.93588E-03 2.96045E-03 2.98522E-03 3.01020E-03 3.03539E-03'
     &/'   3.06079E-03 3.08641E-03 3.11223E-03 3.13828E-03 3.16454E-03'
     &/'   3.19102E-03 3.21772E-03 3.24465E-03 3.27180E-03 3.29918E-03'
     &/'   3.32679E-03 3.35463E-03 3.38270E-03 3.41101E-03 3.43955E-03'
     &/'   3.46833E-03 3.49736E-03 3.52662E-03 3.55613E-03 3.58589E-03'
     &/'   3.61590E-03 3.64616E-03 3.67667E-03 3.70743E-03 3.73846E-03'
     &/'   3.76974E-03 3.80129E-03 3.83310E-03 3.86518E-03 3.89752E-03'
     &/'   3.93014E-03 3.96302E-03 3.99619E-03 4.02963E-03 4.06335E-03'
     &/'   4.09735E-03 4.13164E-03 4.16621E-03 4.20107E-03 4.23623E-03'
     &/'   4.27168E-03 4.30743E-03 4.34347E-03 4.37982E-03 4.41647E-03'
     &/'   4.45343E-03 4.49069E-03 4.52827E-03 4.56616E-03 4.60438E-03'
     &/'   4.64291E-03 4.68176E-03 4.72094E-03 4.76044E-03 4.80028E-03'
     &/'   4.84045E-03 4.88095E-03 4.92180E-03 4.96298E-03 5.00451E-03'
     &/'   5.04639E-03 5.08862E-03 5.13120E-03 5.17414E-03 5.21744E-03'
     &/'   5.26110E-03 5.30513E-03 5.34952E-03 5.39429E-03 5.43943E-03'
     &/'   5.48494E-03 5.53084E-03 5.57713E-03 5.62380E-03 5.67086E-03'
     &/'   5.71831E-03 5.76616E-03 5.81442E-03 5.86307E-03 5.91214E-03'
     &/'   5.96161E-03 6.01150E-03 6.06180E-03 6.11253E-03 6.16368E-03'
     &/'   6.21526E-03 6.26727E-03 6.31971E-03 6.37260E-03 6.42592E-03'
     &/'   6.47970E-03 6.53392E-03 6.58860E-03 6.64373E-03 6.69933E-03'
     &/'   6.75539E-03 6.81192E-03 6.86892E-03 6.92640E-03 6.98436E-03'
     &/'   7.04281E-03 7.10174E-03 7.16117E-03 7.22110E-03 7.28152E-03'
     &/'   7.34246E-03 7.40390E-03 7.46586E-03 7.52833E-03 7.59133E-03'
     &/'   7.65486E-03 7.71891E-03 7.78351E-03 7.84864E-03 7.91432E-03'
     &/'   7.98055E-03 8.04733E-03 8.11467E-03 8.18258E-03 8.25105E-03'
     &/'   8.32009E-03 8.38972E-03 8.45993E-03 8.53072E-03 8.60211E-03'
     &/'   8.67409E-03 8.74668E-03 8.81987E-03 8.89367E-03 8.96810E-03'
     &/'   9.04315E-03 9.11882E-03 9.19513E-03 9.27207E-03 9.34966E-03'
     &/'   9.42790E-03 9.50680E-03 9.58635E-03 9.66657E-03 9.74746E-03'
     &/'   9.82903E-03 9.91128E-03 9.99422E-03 1.00779E-02 1.01622E-02'
     &/'   1.02472E-02 1.03330E-02 1.04194E-02 1.05066E-02 1.05946E-02'
     &/'   1.06832E-02 1.07726E-02 1.08628E-02 1.09537E-02 1.10453E-02'
     &/'   1.11378E-02 1.12310E-02 1.13249E-02 1.14197E-02 1.15153E-02'
     &/'   1.16116E-02 1.17088E-02 1.18068E-02 1.19056E-02 1.20052E-02'
     &/'   1.21057E-02 1.22070E-02 1.23091E-02 1.24121E-02 1.25160E-02'
     &/'   1.26207E-02 1.27263E-02 1.28328E-02 1.29402E-02 1.30485E-02'
     &/'   1.31577E-02 1.32678E-02 1.33788E-02 1.34908E-02 1.36037E-02'
     &/'   1.37175E-02 1.38323E-02 1.39481E-02 1.40648E-02 1.41825E-02'
     &/'   1.43012E-02 1.44208E-02 1.45415E-02 1.46632E-02 1.47859E-02'
     &/'   1.49096E-02 1.50344E-02 1.51602E-02 1.52871E-02 1.54150E-02'
     &/'   1.55440E-02 1.56741E-02 1.58052E-02 1.59375E-02 1.60709E-02'
     &/'   1.62053E-02 1.63409E-02 1.64777E-02 1.66156E-02 1.67546E-02'
     &/'   1.68948E-02 1.70362E-02 1.71788E-02 1.73225E-02 1.74675E-02'
     &/'   1.76136E-02 1.77610E-02 1.79097E-02 1.80595E-02 1.82107E-02'
     &/'   1.83630E-02 1.85167E-02 1.86717E-02 1.88279E-02 1.89855E-02'
     &/'   1.91443E-02 1.93045E-02 1.94661E-02 1.96290E-02 1.97932E-02'
     &/'   1.99589E-02 2.01259E-02 2.02943E-02 2.04641E-02 2.06354E-02'
     &/'   2.08081E-02 2.09822E-02 2.11578E-02 2.13348E-02 2.15133E-02'
     &/'   2.16934E-02 2.18749E-02 2.20580E-02 2.22425E-02 2.24287E-02'
     &/'   2.26164E-02 2.28056E-02 2.29965E-02 2.31889E-02 2.33830E-02'
     &/'   2.35786E-02 2.37759E-02 2.39749E-02 2.41755E-02 2.43778E-02'
     &/'   2.45818E-02 2.47875E-02 2.49950E-02 2.52041E-02 2.54150E-02'
     &/'   2.56277E-02 2.58422E-02 2.60584E-02 2.62765E-02 2.64963E-02'
     &/'   2.67181E-02 2.69417E-02 2.70000E-02 2.71671E-02 2.73945E-02'
     &/'   2.76237E-02 2.78548E-02 2.80879E-02 2.83230E-02 2.85000E-02'
     &/'   2.85600E-02 2.87990E-02 2.90400E-02 2.92830E-02 2.95280E-02'
     &/'   2.97751E-02 3.00243E-02 3.02755E-02 3.05289E-02 3.07844E-02'
     &/'   3.10420E-02 3.13017E-02 3.15637E-02 3.18278E-02 3.20942E-02'
     &/'   3.23627E-02 3.26335E-02 3.29066E-02 3.31820E-02 3.34597E-02'
     &/'   3.37397E-02 3.40220E-02 3.43067E-02 3.45938E-02 3.48833E-02'
     &/'   3.51752E-02 3.54695E-02 3.57663E-02 3.60656E-02 3.63674E-02'
     &/'   3.66718E-02 3.69786E-02 3.72881E-02 3.76001E-02 3.79148E-02'
     &/'   3.82320E-02 3.85520E-02 3.88746E-02 3.91999E-02 3.95279E-02'
     &/'   3.98587E-02 4.01922E-02 4.05286E-02 4.08677E-02 4.12097E-02'
     &/'   4.15546E-02 4.19023E-02 4.22529E-02 4.26065E-02 4.29631E-02'
     &/'   4.33226E-02 4.36851E-02 4.40507E-02 4.44193E-02 4.47910E-02'
     &/'   4.51658E-02 4.55438E-02 4.59249E-02 4.63092E-02 4.66967E-02'
     &/'   4.70875E-02 4.74815E-02 4.78788E-02 4.82795E-02 4.86835E-02'
     &/'   4.90909E-02 4.95017E-02 4.99159E-02 5.03336E-02 5.07548E-02'
     &/'   5.11796E-02 5.16078E-02 5.20397E-02 5.24752E-02 5.29143E-02'
     &/'   5.33571E-02 5.38036E-02 5.42538E-02 5.47078E-02 5.51656E-02'
     &/'   5.56273E-02 5.60928E-02 5.65622E-02 5.70355E-02 5.75128E-02'
     &/'   5.79941E-02 5.84793E-02 5.89687E-02 5.94622E-02 5.99598E-02'
     &/'   6.04615E-02 6.09675E-02 6.14777E-02 6.19921E-02 6.25109E-02'
     &/'   6.30340E-02 6.35614E-02 6.40933E-02 6.46297E-02 6.51705E-02'
     &/'   6.57159E-02 6.62658E-02 6.68203E-02 6.73795E-02 6.79433E-02'
     &/'   6.85119E-02 6.90852E-02 6.96633E-02 7.02463E-02 7.08341E-02'
     &/'   7.14268E-02 7.20246E-02 7.26273E-02 7.32350E-02 7.38479E-02'
     &/'   7.44658E-02 7.50890E-02 7.57173E-02 7.63509E-02 7.69899E-02'
     &/'   7.76341E-02 7.82838E-02 7.89389E-02 7.95000E-02 7.95994E-02'
     &/'   8.02655E-02 8.09372E-02 8.16145E-02 8.22975E-02 8.25000E-02'
     &/'   8.29861E-02 8.36806E-02 8.43808E-02 8.50869E-02 8.57990E-02'
     &/'   8.65169E-02 8.72409E-02 8.79710E-02 8.87071E-02 8.94494E-02'
     &/'   9.01980E-02 9.09528E-02 9.17139E-02 9.24814E-02 9.32552E-02'
     &/'   9.40356E-02 9.48225E-02 9.56160E-02 9.64161E-02 9.72230E-02'
     &/'   9.80366E-02 9.88569E-02 9.96842E-02 1.00518E-01 1.01360E-01'
     &/'   1.02208E-01 1.03063E-01 1.03925E-01 1.04795E-01 1.05672E-01'
     &/'   1.06556E-01 1.07448E-01 1.08347E-01 1.09254E-01 1.10168E-01'
     &/'   1.11090E-01 1.12020E-01 1.12957E-01 1.13902E-01 1.14855E-01'
     &/'   1.15816E-01 1.16786E-01 1.17763E-01 1.18748E-01 1.19742E-01'
     &/'   1.20744E-01 1.21754E-01 1.22773E-01 1.23801E-01 1.24837E-01'
     &/'   1.25881E-01 1.26935E-01 1.27997E-01 1.29068E-01 1.30148E-01'
     &/'   1.31237E-01 1.32336E-01 1.33443E-01 1.34560E-01 1.35686E-01'
     &/'   1.36821E-01 1.37966E-01 1.39120E-01 1.40285E-01 1.41459E-01'
     &/'   1.42642E-01 1.43836E-01 1.45040E-01 1.46253E-01 1.47477E-01'
     &/'   1.48711E-01 1.49956E-01 1.51211E-01 1.52476E-01 1.53752E-01'
     &/'   1.55038E-01 1.56336E-01 1.57644E-01 1.58963E-01 1.60294E-01'
     &/'   1.61635E-01 1.62988E-01 1.64351E-01 1.65727E-01 1.67114E-01'
     &/'   1.68512E-01 1.69922E-01 1.71344E-01 1.72778E-01 1.74224E-01'
     &/'   1.75682E-01 1.77152E-01 1.78634E-01 1.80129E-01 1.81636E-01'
     &/'   1.83156E-01 1.84689E-01 1.86235E-01 1.87793E-01 1.89364E-01'
     &/'   1.90949E-01 1.92547E-01 1.94158E-01 1.95783E-01 1.97421E-01'
     &/'   1.99073E-01 2.00739E-01 2.02419E-01 2.04113E-01 2.05821E-01'
     &/'   2.07543E-01 2.09280E-01 2.11031E-01 2.12797E-01 2.14578E-01'
     &/'   2.16374E-01 2.18184E-01 2.20010E-01 2.21851E-01 2.23708E-01'
     &/'   2.25580E-01 2.27467E-01 2.29371E-01 2.31290E-01 2.33226E-01'
     &/'   2.35178E-01 2.37146E-01 2.39130E-01 2.41131E-01 2.43149E-01'
     &/'   2.45184E-01 2.47235E-01 2.49304E-01 2.51390E-01 2.53494E-01'
     &/'   2.55615E-01 2.57754E-01 2.59911E-01 2.62086E-01 2.64279E-01'
     &/'   2.66491E-01 2.68721E-01 2.70970E-01 2.73237E-01 2.75524E-01'
     &/'   2.77829E-01 2.80154E-01 2.82499E-01 2.84863E-01 2.87246E-01'
     &/'   2.89650E-01 2.92074E-01 2.94518E-01 2.96983E-01 2.97200E-01'
     &/'   2.98500E-01 2.99468E-01 3.01974E-01 3.04501E-01 3.07049E-01'
     &/'   3.09618E-01 3.12209E-01 3.14822E-01 3.17456E-01 3.20113E-01'
     &/'   3.22792E-01 3.25493E-01 3.28217E-01 3.30963E-01 3.33733E-01'
     &/'   3.36525E-01 3.39342E-01 3.42181E-01 3.45045E-01 3.47932E-01'
     &/'   3.50844E-01 3.53780E-01 3.56740E-01 3.59725E-01 3.62735E-01'
     &/'   3.65771E-01 3.68832E-01 3.71918E-01 3.75030E-01 3.78169E-01'
     &/'   3.81333E-01 3.84524E-01 3.87742E-01 3.90987E-01 3.94259E-01'
     &/'   3.97558E-01 4.00885E-01 4.04239E-01 4.07622E-01 4.11033E-01'
     &/'   4.14473E-01 4.17941E-01 4.21438E-01 4.24965E-01 4.28521E-01'
     &/'   4.32107E-01 4.35723E-01 4.39369E-01 4.43046E-01 4.46754E-01'
     &/'   4.50492E-01 4.54262E-01 4.58063E-01 4.61896E-01 4.65762E-01'
     &/'   4.69659E-01 4.73589E-01 4.77552E-01 4.81548E-01 4.85578E-01'
     &/'   4.89642E-01 4.93739E-01 4.97871E-01 5.02037E-01 5.06238E-01'
     &/'   5.10474E-01 5.14746E-01 5.19054E-01 5.23397E-01 5.27777E-01'
     &/'   5.32193E-01 5.36647E-01 5.41138E-01 5.45666E-01 5.50232E-01'
     &/'   5.54837E-01 5.59480E-01 5.64161E-01 5.68882E-01 5.73643E-01'
     &/'   5.78443E-01 5.83284E-01 5.88165E-01 5.93087E-01 5.98050E-01'
     &/'   6.03054E-01 6.08101E-01 6.13189E-01 6.18321E-01 6.23495E-01'
     &/'   6.28712E-01 6.33973E-01 6.39279E-01 6.44628E-01 6.50022E-01'
     &/'   6.55462E-01 6.60947E-01 6.66478E-01 6.72055E-01 6.77679E-01'
     &/'   6.83350E-01 6.89068E-01 6.94834E-01 7.00649E-01 7.06512E-01'
     &/'   7.12424E-01 7.18386E-01 7.24398E-01 7.30459E-01 7.36572E-01'
     &/'   7.42736E-01 7.48951E-01 7.55218E-01 7.61538E-01 7.67911E-01'
     &/'   7.74337E-01 7.80817E-01 7.87351E-01 7.93939E-01 8.00583E-01'
     &/'   8.07282E-01 8.14038E-01 8.20850E-01 8.27719E-01 8.34646E-01'
     &/'   8.41630E-01 8.48673E-01 8.55775E-01 8.62936E-01 8.70157E-01'
     &/'   8.77439E-01 8.84781E-01 8.92185E-01 8.99651E-01 9.07180E-01'
     &/'   9.14771E-01 9.22426E-01 9.30145E-01 9.37928E-01 9.45777E-01'
     &/'   9.53692E-01 9.61672E-01 9.69720E-01 9.77834E-01 9.86017E-01'
     &/'   9.94268E-01 1.00259E+00 1.01098E+00 1.01944E+00 1.02797E+00'
     &/'   1.03657E+00 1.04524E+00 1.05399E+00 1.06281E+00 1.07171E+00'
     &/'   1.08067E+00 1.08972E+00 1.09884E+00 1.10803E+00 1.11730E+00'
     &/'   1.12665E+00 1.13608E+00 1.14559E+00 1.15518E+00 1.16484E+00'
     &/'   1.17459E+00 1.18442E+00 1.19433E+00 1.20432E+00 1.21440E+00'
     &/'   1.22456E+00 1.23481E+00 1.24514E+00 1.25556E+00 1.26607E+00'
     &/'   1.27667E+00 1.28735E+00 1.29812E+00 1.30898E+00 1.31994E+00'
     &/'   1.33098E+00 1.34212E+00 1.35335E+00 1.36468E+00 1.37610E+00'
     &/'   1.38761E+00 1.39922E+00 1.41093E+00 1.42274E+00 1.43465E+00'
     &/'   1.44665E+00 1.45876E+00 1.47096E+00 1.48327E+00 1.49569E+00'
     &/'   1.50820E+00 1.52082E+00 1.53355E+00 1.54638E+00 1.55932E+00'
     &/'   1.57237E+00 1.58553E+00 1.59880E+00 1.61218E+00 1.62567E+00'
     &/'   1.63927E+00 1.65299E+00 1.66682E+00 1.68077E+00 1.69483E+00'
     &/'   1.70902E+00 1.72332E+00 1.73774E+00 1.75228E+00 1.76694E+00'
     &/'   1.78173E+00 1.79664E+00 1.81168E+00 1.82684E+00 1.84212E+00'
     &/'   1.85754E+00 1.87308E+00 1.88876E+00 1.90456E+00 1.92050E+00'
     &/'   1.93657E+00 1.95278E+00 1.96912E+00 1.98560E+00 2.00221E+00'
     &/'   2.01896E+00 2.03586E+00 2.05290E+00 2.07008E+00 2.08740E+00'
     &/'   2.10487E+00 2.12248E+00 2.14024E+00 2.15815E+00 2.17621E+00'
     &/'   2.19442E+00 2.21278E+00 2.23130E+00 2.24997E+00 2.26880E+00'
     &/'   2.28779E+00 2.30693E+00 2.32624E+00 2.34570E+00 2.36533E+00'
     &/'   2.38513E+00 2.40508E+00 2.42521E+00 2.44550E+00 2.46597E+00'
     &/'   2.48660E+00 2.50741E+00 2.52840E+00 2.54955E+00 2.57089E+00'
     &/'   2.59240E+00 2.61410E+00 2.63597E+00 2.65803E+00 2.68027E+00'
     &/'   2.70270E+00 2.72532E+00 2.74812E+00 2.77112E+00 2.79431E+00'
     &/'   2.81769E+00 2.84127E+00 2.86505E+00 2.88902E+00 2.91320E+00'
     &/'   2.93758E+00 2.96216E+00 2.98695E+00 3.01194E+00 3.03715E+00'
     &/'   3.06256E+00 3.08819E+00 3.11403E+00 3.14009E+00 3.16637E+00'
     &/'   3.19286E+00 3.21958E+00 3.24652E+00 3.27369E+00 3.30109E+00'
     &/'   3.32871E+00 3.35657E+00 3.38465E+00 3.41298E+00 3.44154E+00'
     &/'   3.47034E+00 3.49938E+00 3.52866E+00 3.55819E+00 3.58796E+00'
     &/'   3.61799E+00 3.64826E+00 3.67879E+00 3.70958E+00 3.74062E+00'
     &/'   3.77192E+00 3.80349E+00 3.83532E+00 3.86741E+00 3.89977E+00'
     &/'   3.93241E+00 3.96531E+00 3.99850E+00 4.03196E+00 4.06570E+00'
     &/'   4.09972E+00 4.13403E+00 4.16862E+00 4.20350E+00 4.23868E+00'
     &/'   4.27415E+00 4.30992E+00 4.34598E+00 4.38235E+00 4.41902E+00'
     &/'   4.45600E+00 4.49329E+00 4.53089E+00 4.56880E+00 4.60704E+00'
     &/'   4.64559E+00 4.68446E+00 4.72367E+00 4.76319E+00 4.80305E+00'
     &/'   4.84325E+00 4.88378E+00 4.92464E+00 4.96585E+00 5.00741E+00'
     &/'   5.04931E+00 5.09156E+00 5.13417E+00 5.17714E+00 5.22046E+00'
     &/'   5.26414E+00 5.30820E+00 5.35261E+00 5.39741E+00 5.44257E+00'
     &/'   5.48812E+00 5.53404E+00 5.58035E+00 5.62705E+00 5.67414E+00'
     &/'   5.72162E+00 5.76950E+00 5.81778E+00 5.86646E+00 5.91555E+00'
     &/'   5.96506E+00 6.01497E+00 6.06531E+00 6.11606E+00 6.16724E+00'
     &/'   6.21885E+00 6.27089E+00 6.32337E+00 6.37628E+00 6.42964E+00'
     &/'   6.48344E+00 6.53770E+00 6.59241E+00 6.64757E+00 6.70320E+00'
     &/'   6.75929E+00 6.81586E+00 6.87289E+00 6.93041E+00 6.98840E+00'
     &/'   7.04688E+00 7.10585E+00 7.16531E+00 7.22527E+00 7.28574E+00'
     &/'   7.34670E+00 7.40818E+00 7.47018E+00 7.53269E+00 7.59572E+00'
     &/'   7.65928E+00 7.72338E+00 7.78801E+00 7.85318E+00 7.91890E+00'
     &/'   7.98516E+00 8.05198E+00 8.11936E+00 8.18731E+00 8.25582E+00'
     &/'   8.32491E+00 8.39457E+00 8.46482E+00 8.53565E+00 8.60708E+00'
     &/'   8.67910E+00 8.75173E+00 8.82497E+00 8.89882E+00 8.97328E+00'
     &/'   9.04837E+00 9.12409E+00 9.20044E+00 9.27744E+00 9.35507E+00'
     &/'   9.43335E+00 9.51229E+00 9.59190E+00 9.67216E+00 9.75310E+00'
     &/'   9.83472E+00 9.91701E+00 1.00000E+01 1.00837E+01 1.01681E+01'
     &/'   1.02532E+01 1.03390E+01 1.04255E+01 1.05127E+01 1.06007E+01'
     &/'   1.06894E+01 1.07788E+01 1.08690E+01 1.09600E+01 1.10517E+01'
     &/'   1.11442E+01 1.12374E+01 1.13315E+01 1.14263E+01 1.15219E+01'
     &/'   1.16183E+01 1.17156E+01 1.18136E+01 1.19125E+01 1.20122E+01'
     &/'   1.21127E+01 1.22140E+01 1.23162E+01 1.24193E+01 1.25232E+01'
     &/'   1.26280E+01 1.27337E+01 1.28402E+01 1.29477E+01 1.30560E+01'
     &/'   1.31653E+01 1.32755E+01 1.33866E+01 1.34986E+01 1.36116E+01'
     &/'   1.37254E+01 1.38403E+01 1.39561E+01 1.40729E+01 1.41907E+01'
     &/'   1.43094E+01 1.44292E+01 1.45499E+01 1.46717E+01 1.47944E+01'
     &/'   1.49182E+01 1.50431E+01 1.51690E+01 1.52959E+01 1.54239E+01'
     &/'   1.55530E+01 1.56831E+01 1.58144E+01 1.59467E+01 1.60801E+01'
     &/'   1.62147E+01 1.63504E+01 1.64872E+01 1.66252E+01 1.67643E+01'
     &/'   1.69046E+01 1.70460E+01 1.71887E+01 1.73325E+01 1.74776E+01'
     &/'   1.76238E+01 1.77713E+01 1.79200E+01 1.80700E+01 1.82212E+01'
     &/'   1.83737E+01 1.85274E+01 1.86825E+01 1.88388E+01 1.89964E+01'
     &/'   1.91554E+01 1.93157E+01 1.94773E+01 2.00001E+01'
     &/'     unit =    1'
     &/'     axis =  dchain') !FURUTA20200601

  503 format('[ T-TRACK ] off'
     &/'     part =  neutron'
     &/'   e-type =    1'
     &/'       ne =  175'
     &/'   1.00001E-11 1.00001E-07 4.13994E-07 5.31578E-07 6.82560E-07'
     &/'   8.76425E-07 1.12535E-06 1.44498E-06 1.85539E-06 2.38237E-06'
     &/'   3.05902E-06 3.92786E-06 5.04348E-06 6.47595E-06 8.31529E-06'
     &/'   1.06770E-05 1.37096E-05 1.76035E-05 2.26033E-05 2.90232E-05'
     &/'   3.72665E-05 4.78512E-05 6.14421E-05 7.88932E-05 1.01301E-04'
     &/'   1.30073E-04 1.67017E-04 2.14454E-04 2.75364E-04 3.53575E-04'
     &/'   4.53999E-04 5.82947E-04 7.48518E-04 9.61116E-04 1.23410E-03'
     &/'   1.58461E-03 2.03468E-03 2.24867E-03 2.48517E-03 2.61259E-03'
     &/'   2.74654E-03 3.03539E-03 3.35463E-03 3.70744E-03 4.30742E-03'
     &/'   5.53084E-03 7.10174E-03 9.11882E-03 1.05946E-02 1.17088E-02'
     &/'   1.50344E-02 1.93045E-02 2.18749E-02 2.35786E-02 2.41755E-02'
     &/'   2.47875E-02 2.60584E-02 2.70001E-02 2.85011E-02 3.18278E-02'
     &/'   3.43067E-02 4.08677E-02 4.63092E-02 5.24752E-02 5.65622E-02'
     &/'   6.73794E-02 7.19981E-02 7.94987E-02 8.25034E-02 8.65169E-02'
     &/'   9.80366E-02 1.11090E-01 1.16786E-01 1.22773E-01 1.29068E-01'
     &/'   1.35686E-01 1.42642E-01 1.49956E-01 1.57644E-01 1.65727E-01'
     &/'   1.74224E-01 1.83156E-01 1.92547E-01 2.02419E-01 2.12797E-01'
     &/'   2.23708E-01 2.35177E-01 2.47235E-01 2.73237E-01 2.87247E-01'
     &/'   2.94518E-01 2.97211E-01 2.98491E-01 3.01974E-01 3.33733E-01'
     &/'   3.68832E-01 3.87742E-01 4.07622E-01 4.50492E-01 4.97871E-01'
     &/'   5.23397E-01 5.50232E-01 5.78444E-01 6.08101E-01 6.39279E-01'
     &/'   6.72055E-01 7.06512E-01 7.42736E-01 7.80817E-01 8.20850E-01'
     &/'   8.62936E-01 9.07180E-01 9.61640E-01 1.00259E+00 1.10803E+00'
     &/'   1.16484E+00 1.22456E+00 1.28735E+00 1.35335E+00 1.42274E+00'
     &/'   1.49569E+00 1.57237E+00 1.65299E+00 1.73774E+00 1.82684E+00'
     &/'   1.92050E+00 2.01897E+00 2.12248E+00 2.23130E+00 2.30686E+00'
     &/'   2.34570E+00 2.36525E+00 2.38521E+00 2.46597E+00 2.59240E+00'
     &/'   2.72532E+00 2.86505E+00 3.01194E+00 3.16637E+00 3.32871E+00'
     &/'   3.67879E+00 4.06570E+00 4.49329E+00 4.72367E+00 4.96585E+00'
     &/'   5.22046E+00 5.48812E+00 5.76950E+00 6.06531E+00 6.37628E+00'
     &/'   6.59238E+00 6.70320E+00 7.04688E+00 7.40818E+00 7.78801E+00'
     &/'   8.18731E+00 8.60708E+00 9.04837E+00 9.51229E+00 1.00000E+01'
     &/'   1.05127E+01 1.10517E+01 1.16183E+01 1.22140E+01 1.25232E+01'
     &/'   1.28400E+01 1.34986E+01 1.38403E+01 1.41907E+01 1.45499E+01'
     &/'   1.49183E+01 1.56831E+01 1.64872E+01 1.69046E+01 1.73325E+01'
     &/'   1.96403E+01'
     &/'     unit =    1'
     &/'     axis =  dchain') !FURUTA20200601

*-----------------------------------------------------------------------
*        title
*-----------------------------------------------------------------------

                  ittll(itnm) = titll
                  ittle(itnm) = title

               if( title(1:9) .eq. '[t-yield]' ) then
                  title = '[t-dchain]'//' '//title(10:ittll(itnm))
                  ittll(itnm+2) = titll + 2
                  ittle(itnm+2) = title

                  title = '[t-track]'//' '//title(12:ittll(itnm+2))
                  write(iodw,'(4x,"title =  ",a)')
     &            title(1:ittll(itnm)+1)
               else
                  ittll(itnm+2) = titll
                  ittle(itnm+2) = title
                  write(iodw,'(4x,"title =  ",a)')
     &            title(1:ittll(itnm+2))
               end if

*-----------------------------------------------------------------------
*        angel parameter
*-----------------------------------------------------------------------

               if( langel .gt. 0 ) then

                  itanl(itnm) = langel
                  itang(itnm) = angelp

               else

                  itanl(itnm) = 0

               end if

*-----------------------------------------------------------------------
*        x-txt parameter
*-----------------------------------------------------------------------

               if( lxtxt .gt. 0 ) then

                  itaxl(itnm) = lxtxt
                  itaxt(itnm) = cxtxt

               else

                  itaxl(itnm) = 0

               end if

*-----------------------------------------------------------------------
*        y-txt parameter
*-----------------------------------------------------------------------

               if( lytxt .gt. 0 ) then

                  itayl(itnm) = lytxt
                  itayt(itnm) = cytxt

               else

                  itayl(itnm) = 0

               end if

*-----------------------------------------------------------------------
*        z-txt parameter
*-----------------------------------------------------------------------

               if( lztxt .gt. 0 ) then

                  itazl(itnm) = lztxt
                  itazt(itnm) = cztxt

               else

                  itazl(itnm) = 0

               end if

*-----------------------------------------------------------------------
*        gshow parameter
*-----------------------------------------------------------------------

               if( lgshow .gt. 0 ) then

                  itgsh(itnm) = min(5,lgshow)  ! T.Sato 2019/02/20

               else

                  itgsh(itnm) = 0

               end if

*-----------------------------------------------------------------------
*        gslat parameter
*-----------------------------------------------------------------------
                  itglt(itnm) = igslt

*-----------------------------------------------------------------------
*        rshow parameter
*-----------------------------------------------------------------------

                  itrsh(itnm) = min(4,lrshow)

*-----------------------------------------------------------------------
*        special  ! MATSUDA 2017.03.21
*-----------------------------------------------------------------------

               if( iout .eq. 1 ) then

                  itspc(itnm) = 0
                 if( ispec .gt.0 )
     &              write(6,'(a68)') 'Warning: "special" parameter is
     & effective only when output = product'

               else

                  itspc(itnm) = ispec

               end if

*-----------------------------------------------------------------------
*        output  ! MATSUDA 2017.03.21 ! 0: product, 1: cutoff (defaukt)
*-----------------------------------------------------------------------

                  itprd(itnm) = iout

*-----------------------------------------------------------------------
*        jelas; only without elastic collision for t-yield
*-----------------------------------------------------------------------

                  itdpo(itnm) = -1

*-----------------------------------------------------------------------
*        mother
*-----------------------------------------------------------------------

                  itman(itnm) = imoth
                  itmct(itnm) = jmoth

               if( imoth .gt. 0 ) then

                  itsmn(itnm) = itsmn(itnm)
     &                        + ( imoth + mod(imoth,2) ) / 2

               end if

*-----------------------------------------------------------------------
*        nucleus
*-----------------------------------------------------------------------

                  itnun(itnm) = inucl

               if( inucl .gt. 0 ) then

                  itsmn(itnm) = itsmn(itnm)
     &                        + ( inucl + mod(inucl,2) ) / 2

               end if

*-----------------------------------------------------------------------
*        material
*-----------------------------------------------------------------------

                  itmtn(itnm) = imate
                  itmcn(itnm) = jmate

            if( imate .gt. 0 ) then

                  itsmn(itnm) = itsmn(itnm)
     &                        + ( imate + mod(imate,2) ) / 2

            end if

*-----------------------------------------------------------------------
*        mesh
*-----------------------------------------------------------------------

                  itmsh(itnm) = imesh

*           mesh for DCHAIN mode is allowed only region.
         if( imesh .eq. 1 ) then

                  itrgn(itnm) = ntrn
                  itrgm(itnm) = mtrn
                  itrnv(itnm) = nvol

            if( ntrn .gt. 0 ) then

                  call moddas_reallocate_int(
     &                    itlmax, itnm, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(itnm):itreg(itnm)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)
                  call moddas_deallocate_int(idas_itreg_temporary)
                  itsmn(itnm) = itsmn(itnm)
     &                        + ( mtrn + mod(mtrn,2) ) / 2

            end if

            if( nvol .gt. 0 ) then

                  itriv(itnm) = ivl
                  itsmn(itnm) = itsmn(itnm)
     &                        + ( nvol + mod(nvol,2) ) / 2

                  itrrv(itnm) = irvl
                  itsmn(itnm) = itsmn(itnm) + nvol

            end if

         write(iodw,'(5x,"mesh = reg")')

         do iii = 1, idcrg
             write(iodw,'(200a)') dcrgm(iii)
         end do

*---
*        mesh
*CC               call echrg(1,iechrl,0,1,iodw,mtrn,
*CC  &                       idas(ndsm),ntrn,mtrn,idas(ndsm),
*CC  &                       idas(itsmn(itnm)),idas(mtrn),nvol,
*CC  &                       idas(ivl),idas(irvl),ndsm)
*    &                       idas(itrcg(itnm)),itrgn(itnm),
*    &                       itrgm(itnm),idas(itreg(itnm)),
*    &                       das(idas1),idas(idas2),itrnv(itnm),
*    &                       idas(itriv(itnm)),das(itrrv(itnm)),
*    &                       idas3)
**---
         elseif( imesh .eq. 2 ) then

CCSE change and add for mesh=r-z parameter (2017.11.30) >>>>>
            rtrx0(itnm) = rzx0
            rtry0(itnm) = rzy0

            itrty(itnm) = irtp
            itrnm(itnm) = inr
            rtrdl(itnm) = rdel
            rtrmi(itnm) = rmin
            rtrma(itnm) = rmax
            call moddas_reallocate_dbl(
     &              itlmax, itnm, inr+1, itrrg, das_itrrg)
            das_itrrg(itrrg(itnm):itrrg(itnm)+inr)
     &               = gmsh(istrg:istrg+inr)
            itsmn(itnm) = itsmn(itnm) + inr + 1

            itzty(itnm) = iztp
            itznm(itnm) = inz
            rtzdl(itnm) = zdel
            rtzmi(itnm) = zmin
            rtzma(itnm) = zmax
            call moddas_reallocate_dbl(
     &              itlmax, itnm, inz+1, itzrg, das_itzrg)
            das_itzrg(itzrg(itnm):itzrg(itnm)+inz)
     &               = gmsh(istzg:istzg+inz)
            itsmn(itnm) = itsmn(itnm) + inz + 1

            write(iodw,'(5x,"mesh = r-z")')
            do iii = 1, idcrg
               write(iodw,'(200a)') dcrgm(iii)
            end do
CCSE change and add for mesh=r-z parameter (2017.11.30) <<<<<

         else if( imesh .eq. 3 .or. lrshow .gt. 0 ) then

CCSE change and add for mesh=xyz parameter (2017.11.30) >>>>>
            itxty(itnm) = ixtp
            itxnm(itnm) = inx
            rtxdl(itnm) = xdel
            rtxmi(itnm) = xmin
            rtxma(itnm) = xmax
            call moddas_reallocate_dbl(
     &              itlmax, itnm, inx+1, itxrg, das_itxrg)
            das_itxrg(itxrg(itnm):itxrg(itnm)+inx)
     &               = gmsh(istxg:istxg+inx)
            itsmn(itnm) = itsmn(itnm) + inx + 1

            ityty(itnm) = iytp
            itynm(itnm) = iny
            rtydl(itnm) = ydel
            rtymi(itnm) = ymin
            rtyma(itnm) = ymax
            call moddas_reallocate_dbl(
     &              itlmax, itnm, iny+1, ityrg, das_ityrg)
            das_ityrg(ityrg(itnm):ityrg(itnm)+iny)
     &               = gmsh(istyg:istyg+iny)
            itsmn(itnm) = itsmn(itnm) + iny + 1

            itzty(itnm) = iztp
            itznm(itnm) = inz
            rtzdl(itnm) = zdel
            rtzmi(itnm) = zmin
            rtzma(itnm) = zmax
            call moddas_reallocate_dbl(
     &              itlmax, itnm, inz+1, itzrg, das_itzrg)
            das_itzrg(itzrg(itnm):itzrg(itnm)+inz)
     &               = gmsh(istzg:istzg+inz)
            itsmn(itnm) = itsmn(itnm) + inz + 1

            write(iodw,'(5x,"mesh = xyz")')
            do iii = 1, idcrg
               write(iodw,'(200a)') dcrgm(iii)
            end do
CCSE change and add for mesh=xyz parameter (2017.11.30) <<<<<

         elseif( imesh .eq. 4 )then

                  itrgn(itnm) = ntrn
                  itrgm(itnm) = mtrn
                  itrnv(itnm) = 0

            if( ntrn .gt. 0 ) then

                  call moddas_reallocate_int(
     &                    itlmax, itnm, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(itnm):itreg(itnm)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)

            end if

            write(iodw,'(5x,"mesh = tet")')

            do iii = 1, idcrg
               write(iodw,'(200a)') dcrgm(iii)
            end do

         end if

*-----------------------------------------------------------------------
*        axis
*-----------------------------------------------------------------------

*           axis for DCHAIN mode is modified.
                  itaxn(itnm) = 1
                  itaxs(itnm,1) = iaxis(1)


*-----------------------------------------------------------------------
*        file
*-----------------------------------------------------------------------

*           file name is modified for DCHAIN mode.
                  itfln(itnm) = infil
                  itfln(itnm+2) = infil
               do i = 1, itfln(itnm)

                  filnm = ifile(i)(1:lfile(i))
                  ifilnm = lfile(i)
***  "RENAME" modified by norihiro MATSUDA  2012.07.26  ----------------


               do j = ifilnm , 1, -1
                  if( filnm(j:j) .eq. '.' ) goto 810
               end do
                  j = ifilnm + 1

  810             itfp = j - 1

                  filnm(1:itfp+5) = filnm(1:itfp)//'.dyld'
                  ifilnm = itfp + 6

               do j = ifilnm, 100
                  filnm(j:j) = ' '
               end do

*------------------------------------------------------- "RENAME" end --

*           for [T-Yield] tally
                  itfll(itnm,i) = ifilnm
                  ctfln(itnm,i)(1:ifilnm) = filnm(1:ifilnm)

                  do k = ifilnm + 1, 100
                     ctfln(itnm,i)(k:k) = ' '
                  end do

*           for [T-Track] tally
                  itfll(itnm+1,i) = ifilnm
                  ctfln(itnm+1,i)(1:ifilnm) =
     &                      filnm(1:ifilnm-6)//'.dtrk'

                  do k = ifilnm + 1, 100
                     ctfln(itnm+1,i)(k:k) = ' '
                  end do
                  write(iodw,'(5x,"file =  ",100a1)')
     &             ( ctfln(itnm+1,i)(j:j), j = 1, 100 )

*           for [T-dchain] tally
                  itfll(itnm+2,i) = lfile(i)
                  ctfln(itnm+2,i)(1:lfile(i)) =
     &                      ifile(i)(1:lfile(i))

                  do k = lfile(i) + 1, 100
                     ctfln(itnm+2,i)(k:k) = ' '
                  end do

               end do

*-----------------------------------------------------------------------
*       restart file
*-----------------------------------------------------------------------

               call summary_param_resfile(itnm,irfflg,lrfile,irfile)

*-----------------------------------------------------------------------
*        particles
*-----------------------------------------------------------------------

            itmxptt = inpat
            itmxgpt = 1
            do i = 1, inpat
              if ( iptyp(i).lt.0 .and. -iptyp(i).gt.itmxgpt ) then
                itmxgpt = -iptyp(i)
              end if
            end do
            call ALLOCATE_PART_TEMP

            if( inpat .gt. 0 ) then

                  itpan(itnm) = inpat

               do i = 1, inpat

                  itpat(itnm,i,1) = iptyp(i)
                  itpat(itnm,i,2) = ipnkf(i)
                  itpat(itnm,i,3) = ipsub(i)  ! kitamura22/03/31

                 if( itpat(itnm,i,1) .lt. 0 ) then

                     do k = 1, -itpat(itnm,i,1)

                        jtpat(itnm,i,k,1) = imtyp(i,k)
                        jtpat(itnm,i,k,2) = imnkf(i,k)

                     end do

                  end if

               end do

            else if( inpat .eq. 0 ) then

                  itpan(itnm) = 1
                  itpat(itnm,1,1) = 20
                  itpat(itnm,1,2) = 0
                  itpat(itnm,1,3) = 0   ! kitamura22/03/31

            end if

*-----------------------------------------------------------------------
*        ndata
*-----------------------------------------------------------------------

                  itnda(itnm) = ndata

*-----------------------------------------------------------------------
*        reg echo length
*-----------------------------------------------------------------------

                  iterl(itnm) = iechrl
                  if( iechrl .ne. 72 ) then
                     write(iodw,'(3x,"iechrl =",i5)') iechrl
                  end if

*-----------------------------------------------------------------------
*        unit
*-----------------------------------------------------------------------
*           unit for DCHAIN mode is allowed '1'.
                  itunt(itnm) = 1

*-----------------------------------------------------------------------
*        info
*-----------------------------------------------------------------------

                  itout(itnm) = info

*-----------------------------------------------------------------------
*        2D-type
*-----------------------------------------------------------------------

                  ittwo(itnm) = idtyp

*-----------------------------------------------------------------------
*        factor
*-----------------------------------------------------------------------

                  rtfac(itnm) = rfact
                  if( rfact .ne. 1.0d0 ) then
                     write(iodw,'(3x,"factor =",1p1g14.7)') rfact
                  end if
                  rtfac(itnm+2) = 1.0

*-----------------------------------------------------------------------
*        volmat
*-----------------------------------------------------------------------

                  itvm(itnm) = matvol

*-----------------------------------------------------------------------
*        epsout
*-----------------------------------------------------------------------

                  iteps(itnm) = ieps
                  if( ieps .ne. 0 ) then
                     write(iodw,'(3x,"epsout =",i5)') ieps
                  end if

*-----------------------------------------------------------------------
*        resolution
*-----------------------------------------------------------------------

                  itres(itnm) = max( 1, ireso )
                  if( itgsh(itnm) .gt. 0 .or. itrsh(itnm) .gt. 0 ) then
                     write(iodw,'(4x,"resol =",i4)') max( 1, ireso )
                  end if

*-----------------------------------------------------------------------
*        width
*-----------------------------------------------------------------------

                  rtwid(itnm) = width
                  if( itgsh(itnm) .gt. 0 .or. itrsh(itnm) .gt. 0 ) then
                     write(iodw,'(4x,"width =",1p1g14.7)') width
                  end if

*-----------------------------------------------------------------------
*        counter
*-----------------------------------------------------------------------

               do i = 1, 9

                  itcnt(i,itnm) = icount(i)

               end do

            do i = 1, 3
               if( icount(i) .ne. 0 ) then
                  write(iodw,'(1x,"ctmin(",i1,") =",i5)')
     &                  i, itcnt(i*2+2,itnm)
                  write(iodw,'(1x,"ctmax(",i1,") =",i5)')
     &                  i, itcnt(i*2+3,itnm)
               end if
            end do

*-----------------------------------------------------------------------
*        stdcut
*-----------------------------------------------------------------------
                  rtstd(itnm) = stdcut

               if( rtstd(itnm) .gt. 0.d0 ) then
                  write(iodw,'("   stdcut =",1p1g14.7)') rtstd(itnm)
               end if

*-----------------------------------------------------------------------
*        transformation
*-----------------------------------------------------------------------

                     itmtr(itnm,1) = igkst
                     itmtr(itnm,2) = ktrs
                     itmtr(itnm,3) = idtt
                     itmtr(itnm,4) = 0

               if( igkst .gt. 1 ) then

                  do k = 1, 13

                     rtmtr(itnm,k) = vtrs(k)

                  end do

                  write(iodw,'(5x,"trcl = ",i6)') idtt

               end if

*-----------------------------------------------------------------------
*        mxnuclei
*-----------------------------------------------------------------------
                  if( imxnuc.eq.-1 ) then ! frtati 2022/03/25
                    itnzn(itnm) = 3000 ! default value
                  else if( imxnuc.gt.0 ) then
                    itnzn(itnm) = imxnuc + 1
                  else
                    itnzn(itnm) = 0
                  end if

*-----------------------------------------------------------------------
*        xy, yz, zx, rz storage region
*-----------------------------------------------------------------------

                     mxyz = 1
                     mrzr = 1
                     mrzz = 1


                  do i = 1, inaxi

                     if( iaxis(i) .eq. 9 ) then

                        if( inx .gt. mxyz ) mxyz = inx
                        if( iny .gt. mxyz ) mxyz = iny

                     else if( iaxis(i) .eq. 10 ) then

                        if( iny .gt. mxyz ) mxyz = iny
                        if( inz .gt. mxyz ) mxyz = inz

                     else if( iaxis(i) .eq. 11 ) then

                        if( inx .gt. mxyz ) mxyz = inx
                        if( inz .gt. mxyz ) mxyz = inz

                     else if( iaxis(i) .eq. 12 ) then

                        mrzr = inr
                        mrzz = inz

                     end if

                  end do

                     itnfn(itnm) = mxyz
                     itnfr(itnm) = mrzr
                     itnfz(itnm) = mrzz

*-----------------------------------------------------------------------
*        iz, in storage region
*-----------------------------------------------------------------------

                        izz = 1
                        inn = 1

               do i = 1, inaxi

                  if( iaxis(i) .eq. 13 .or. iaxis(i) .eq. 8 ) then

                        izz = maxpt
                        inn = maxnt

                  else if( iaxis(i) .eq. 7 ) then

                        izz = maxpt

                  else if( iaxis(i) .eq. 1 ) then

                     if( inucl .eq. 0 ) then

                        if( inaxi .eq. 1 ) then

                           izz = 1
                           inn = maxpt + maxnt

                        else

                           izz = maxpt
                           inn = maxnt

                        end if

                     else if( inucl .gt. 0 ) then

                        izz = max( izz, inucl )
                        inn = maxnt

                     end if

                  else

                        izz = max( izz, inucl )

                  end if

               end do

                      if( itnzn(itnm).eq.0 ) then ! frtati 2022/02/18

                        itndz(itnm) = izz
                        itndn(itnm) = inn
                        itndm(itnm) = 2

                      else
                        itndz(itnm) = itnzn(itnm)
                        itndn(itnm) = 1
                        itndm(itnm) = 0
                      end if

*-----------------------------------------------------------------------
*        jz, jn -> kz, kn tranfer
*-----------------------------------------------------------------------

               do jz = 1, maxpt
               do jn = 1, maxnt

                        kz = izz + 1
                        kn = inn + 1

                  if( izz. eq. maxpt .and. inn .eq. maxnt ) then

                        kz = jz
                        kn = jn

                  else if( izz .eq. 1 .and.
     &                     inn .eq. maxpt + maxnt ) then

                        kz = 1
                        kn = jz + jn

                  else if( inucl .gt. 0 .and. inn .eq. maxnt ) then

                        kn = jn

                     do i = 1, inucl

                        iz = isnuc(nnucl-1+i) / 1000

                        if( jz .eq. iz ) kz = i

                     end do

                  else if( izz. eq. maxpt .and. inn .eq. 1 ) then

                        kz = jz
                        kn = 1

                  else if( inucl .gt. 0 .and. inn .eq. 1 ) then

                        kn = 1

                     do i = 1, inucl

                        iz = isnuc(nnucl-1+i) / 1000
                        ia = isnuc(nnucl-1+i) - iz * 1000

                        if( jz .eq. iz .and.
     &                    ( ( ia .eq. 0 ) .or.
     &                      ( ia .gt. 0 .and.
     &                        jn .eq. ia - iz ) ) ) then

                              kz = i

                        end if

                     end do

                  else if( izz .eq. 1 .and. inn .eq. 1 ) then

                        kz = 1
                        kn = 1

                  end if

                     ikzz(jz,jn) = kz
                     iknn(jz,jn) = kn

               end do
               end do

*-----------------------------------------------------------------------

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
*        initialization  2012/11/19 by matsuda
*-----------------------------------------------------------------------
               do i = 1, nientry
                  itdci(i,itnm+2) = 0
               end do
*-----------------------------------------------------------------------
*        imode
*-----------------------------------------------------------------------

                  itdci(1,itnm+2) = imode

*-----------------------------------------------------------------------
*        jmode
*-----------------------------------------------------------------------

                  itdci(2,itnm+2) = jmode

*-----------------------------------------------------------------------
*        itstep
*-----------------------------------------------------------------------

                  itdci(3,itnm+2) = idcstp

*-----------------------------------------------------------------------
*        itout
*-----------------------------------------------------------------------

                  itdci(4,itnm+2) = idcout

*-----------------------------------------------------------------------
*        idivs
*-----------------------------------------------------------------------

                  itdci(5,itnm+2) = idivs

*-----------------------------------------------------------------------
*        iregon
*-----------------------------------------------------------------------
                  itdci(6,itnm+2) = iregon

*-----------------------------------------------------------------------
*        inmtcf
*-----------------------------------------------------------------------

                  itdci(7,itnm+2) = inmtcf

*-----------------------------------------------------------------------
*        ichain
*-----------------------------------------------------------------------

                  itdci(8,itnm+2) = ichain

*-----------------------------------------------------------------------
*        itdecs
*-----------------------------------------------------------------------

                  itdci(9,itnm+2) = itdecs

*-----------------------------------------------------------------------
*        itdecn
*-----------------------------------------------------------------------

                  itdci(10,itnm+2) = itdecn

*-----------------------------------------------------------------------
*        isomtr
*-----------------------------------------------------------------------

                  itdci(11,itnm+2) = isomtr

*-----------------------------------------------------------------------
*        ifisyd
*-----------------------------------------------------------------------

                  itdci(12,itnm+2) = ifisyd

*-----------------------------------------------------------------------
*        ifisye
*-----------------------------------------------------------------------

                  itdci(13,itnm+2) = ifisye

*-----------------------------------------------------------------------
*        iyild
*-----------------------------------------------------------------------

                  itdci(14,itnm+2) = iyild

*-----------------------------------------------------------------------
*        iggrp
*-----------------------------------------------------------------------

                  itdci(15,itnm+2) = iggrp

*-----------------------------------------------------------------------
*        ibetap
*-----------------------------------------------------------------------

                  itdci(16,itnm+2) = ibetap

*-----------------------------------------------------------------------
*        acmin
*-----------------------------------------------------------------------

                  rtdci(1,itnm+2) = acmin

*-----------------------------------------------------------------------
*        istabl
*-----------------------------------------------------------------------

                  itdci(17,itnm+2) = istabl

*-----------------------------------------------------------------------
*        igsdef
*-----------------------------------------------------------------------

                  itdci(18,itnm+2) = igsdef

*-----------------------------------------------------------------------
*        iprtb1
*-----------------------------------------------------------------------

                  itdci(19,itnm+2) = iprtb1

*-----------------------------------------------------------------------
*        iprtb2
*-----------------------------------------------------------------------

                  itdci(20,itnm+2) = iprtb2

*-----------------------------------------------------------------------
*        rprtb2
*-----------------------------------------------------------------------

                  rtdci(2,itnm+2) = rprtb2

*-----------------------------------------------------------------------
*        iprtb3
*-----------------------------------------------------------------------

                  itdci(21,itnm+2) = iprtb3

*-----------------------------------------------------------------------
*        igsorg
*-----------------------------------------------------------------------

                  itdci(22,itnm+2) = igsorg

*-----------------------------------------------------------------------
*        amp
*-----------------------------------------------------------------------

                  rtdci(3,itnm+2) = amp

*-----------------------------------------------------------------------
*        ebeam
*-----------------------------------------------------------------------

                  rtdci(4,itnm+2) = ebeam

*-----------------------------------------------------------------------
*        prodnp
*-----------------------------------------------------------------------

                  rtdci(5,itnm+2) = prodnp

*-----------------------------------------------------------------------
*        Card[6]
*-----------------------------------------------------------------------

                  do i = 1, idcstp
                     tbina(itnm+2,i) = temsh(1,i)
                     tbinb(itnm+2,i) = itemsh(1,i)
                     beampw(itnm+2,i) = temsh(2,i)
                  end do

*        check
                     if( itemsh(1,idcstp) .eq. 1 ) then
                        dlastt = temsh(1,idcstp)*3600.d0*24.d0*365.d0
                       if( itemsh(2,idcout) .eq. 1 ) then
                           dlastc = temsh(3,idcout)*3600.d0*24.d0*365.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 2 ) then
                           dlastc = temsh(3,idcout)*3600.d0*24.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 3 ) then
                           dlastc = temsh(3,idcout)*3600.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 4 ) then
                           dlastc = temsh(3,idcout)*60.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 5 ) then
                           dlastc = temsh(3,idcout)
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       end if

                     else if( itemsh(1,idcstp) .eq. 2 ) then
                        dlastt = temsh(1,idcstp)*3600.d0*24.d0
                       if( itemsh(2,idcout) .eq. 1 ) then
                           dlastc = temsh(3,idcout)*3600.d0*24.d0*365.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 2 ) then
                           dlastc = temsh(3,idcout)*3600.d0*24.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 3 ) then
                           dlastc = temsh(3,idcout)*3600.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 4 ) then
                           dlastc = temsh(3,idcout)*60.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 5 ) then
                           dlastc = temsh(3,idcout)
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       end if

                     else if( itemsh(1,idcstp) .eq. 3 ) then
                        dlastt = temsh(1,idcstp)*3600.d0
                       if( itemsh(2,idcout) .eq. 1 ) then
                           dlastc = temsh(3,idcout)*3600.d0*24.d0*365.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 2 ) then
                           dlastc = temsh(3,idcout)*3600.d0*24.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 3 ) then
                           dlastc = temsh(3,idcout)*3600.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 4 ) then
                           dlastc = temsh(3,idcout)*60.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 5 ) then
                           dlastc = temsh(3,idcout)
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       end if

                     else if( itemsh(1,idcstp) .eq. 4 ) then
                        dlastt = temsh(1,idcstp)*60.d0
                       if( itemsh(2,idcout) .eq. 1 ) then
                           dlastc = temsh(3,idcout)*3600.d0*24.d0*365.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 2 ) then
                           dlastc = temsh(3,idcout)*3600.d0*24.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 3 ) then
                           dlastc = temsh(3,idcout)*3600.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 4 ) then
                           dlastc = temsh(3,idcout)*60.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 5 ) then
                           dlastc = temsh(3,idcout)
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       end if

                     else if( itemsh(1,idcstp) .eq. 5 ) then
                        dlastt = temsh(1,idcstp)
                       if( itemsh(2,idcout) .eq. 1 ) then
                           dlastc = temsh(3,idcout)*3600.d0*24.d0*365.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 2 ) then
                           dlastc = temsh(3,idcout)*3600.d0*24.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 3 ) then
                           dlastc = temsh(3,idcout)*3600.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 4 ) then
                           dlastc = temsh(3,idcout)*60.d0
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       else if( itemsh(2,idcout) .eq. 5 ) then
                           dlastc = temsh(3,idcout)
                         if ( dlastt + dlastc .lt. 0.d0 ) then
                           tbina(itnm+2,idcstp) = abs(temsh(3,idcout))
                           tbinb(itnm+2,idcstp) = itemsh(2,idcout)
                         end if
                       end if

                     end if

*-----------------------------------------------------------------------
*        Card[7]
*-----------------------------------------------------------------------

                  do i = 1, idcout
                     tmina(itnm+2,i) = temsh(3,i)
                     tminb(itnm+2,i) = itemsh(2,i)
                  end do

*-----------------------------------------------------------------------
*        Card[8]
*-----------------------------------------------------------------------

                  itdci(23,itnm+2) = iancza

                  do i = 1, iancza
                     aclst(itnm+2,i) = hancza(i)
                  end do

*-----------------------------------------------------------------------
*        Card[9]
*-----------------------------------------------------------------------

c H.Ratliff START 2019.07.30

*-----------------------------------------------------------------------
*        iertdcho
*-----------------------------------------------------------------------

                  itdci(26,itnm+2) = iertdcho

*-----------------------------------------------------------------------
*        inxslib
*-----------------------------------------------------------------------

                  itdci(27,itnm+2) = inxslib

*-----------------------------------------------------------------------
*        idcylib
*-----------------------------------------------------------------------

                  itdci(28,itnm+2) = idcylib

*-----------------------------------------------------------------------
*        infylib
*-----------------------------------------------------------------------

                  itdci(29,itnm+2) = infylib

*-----------------------------------------------------------------------
*        hnxslib
*-----------------------------------------------------------------------

                  ihxll(itnm) = hxsll
                  ihxle(itnm) = hnxslib
                  ihxll(itnm+2) = hxsll
                  ihxle(itnm+2) = hnxslib

*-----------------------------------------------------------------------
*        hdcylib
*-----------------------------------------------------------------------

                  ihdll(itnm) = hdcll
                  ihdle(itnm) = hdcylib
                  ihdll(itnm+2) = hdcll
                  ihdle(itnm+2) = hdcylib

*-----------------------------------------------------------------------
*        hnfylib
*-----------------------------------------------------------------------

                  ihfll(itnm) = hfsll
                  ihfle(itnm) = hnfylib
                  ihfll(itnm+2) = hfsll
                  ihfle(itnm+2) = hnfylib

*-----------------------------------------------------------------------
*        hsfylib
*-----------------------------------------------------------------------

                  ihsll(itnm) = hssll
                  ihsle(itnm) = hsfylib
                  ihsll(itnm+2) = hssll
                  ihsle(itnm+2) = hsfylib

*-----------------------------------------------------------------------
*        iwrtchn
*-----------------------------------------------------------------------

                  itdci(30,itnm+2) = iwrtchn

*-----------------------------------------------------------------------
*        chrlvth
*-----------------------------------------------------------------------

                  rtdci(6,itnm+2) = chrlvth

*-----------------------------------------------------------------------
*        iwrchdt
*-----------------------------------------------------------------------

                  itdci(31,itnm+2) = iwrchdt

*-----------------------------------------------------------------------
*        iwrchss
*-----------------------------------------------------------------------

                  itdci(32,itnm+2) = iwrchss

c H.Ratliff START 2020.04.17

*-----------------------------------------------------------------------
*        ixsrall
*-----------------------------------------------------------------------

                  itdci(33,itnm+2) = ixsrall

*-----------------------------------------------------------------------
*        idosecf
*-----------------------------------------------------------------------

                  itdci(34,itnm+2) = idosecf

*-----------------------------------------------------------------------
*        ipltmode
*-----------------------------------------------------------------------

                  itdci(35,itnm+2) = ipltmode

*-----------------------------------------------------------------------
*        ipltaxis
*-----------------------------------------------------------------------

                  itdci(36,itnm+2) = ipltaxis


c H.Ratliff START 2020.05.29

*-----------------------------------------------------------------------
*        foamout
*-----------------------------------------------------------------------

                  itdci(37,itnm+2) = foamout

*-----------------------------------------------------------------------
*        foamvals
*-----------------------------------------------------------------------

                  itdci(38,itnm+2) = foamvals

c H.Ratliff START 2020.06.05
*-----------------------------------------------------------------------
*        iredufmt
*-----------------------------------------------------------------------

                  itdci(39,itnm+2) = iredufmt0
                  iredufmt(itnm)   = iredufmt0
                  iredufmt(itnm+1) = iredufmt0

c H.Ratliff START 2020.06.17
*-----------------------------------------------------------------------
*        irdonce
*-----------------------------------------------------------------------

                  itdci(40,itnm+2) = irdonce

c H.Ratliff START 2020.11.16
*-----------------------------------------------------------------------
*        imtcard
*-----------------------------------------------------------------------

                  itdci(41,itnm+2) = imtcard

*-----------------------------------------------------------------------
*        imtcnum
*-----------------------------------------------------------------------

                  itdci(42,itnm+2) = imtcnum

*-----------------------------------------------------------------------
*        imtcmeta
*-----------------------------------------------------------------------

                  itdci(43,itnm+2) = imtcmeta

*-----------------------------------------------------------------------
*        thmatnd
*-----------------------------------------------------------------------

                  rtdci(7,itnm+2) = thmatnd

c H.Ratliff START 2020.12.18
*-----------------------------------------------------------------------
*        idosunit
*-----------------------------------------------------------------------

                  itdci(44,itnm+2) = idosunit

*-----------------------------------------------------------------------
*        iangbpwr
*-----------------------------------------------------------------------

                  itdci(45,itnm+2) = iangbpwr

c H.Ratliff START 2021.02.17
*-----------------------------------------------------------------------
*        iphtout
*-----------------------------------------------------------------------

                  itdci(46,itnm+2) = iphtout

c H.Ratliff START 2021.04.07
*-----------------------------------------------------------------------
*        isfylib
*-----------------------------------------------------------------------

                  itdci(47,itnm+2) = isfylib

*-----------------------------------------------------------------------
*        radlev
*-----------------------------------------------------------------------

                  rtdci(8,itnm+2) = radlev

*-----------------------------------------------------------------------
*        ergthfa
*-----------------------------------------------------------------------

                  rtdci(9,itnm+2) = ergthfa

*-----------------------------------------------------------------------
*        ergfafu
*-----------------------------------------------------------------------

                  rtdci(10,itnm+2) = ergfafu

*-----------------------------------------------------------------------
*        infyerg
*-----------------------------------------------------------------------

                  itdci(48,itnm+2) = infyerg

c H.Ratliff START 2021.04.26
*-----------------------------------------------------------------------
*        ilchain
*-----------------------------------------------------------------------

                  itdci(49,itnm+2) = ilchain

*-----------------------------------------------------------------------
*        itriton  T.Sato 2025/02/08
*-----------------------------------------------------------------------

                  itdci(50,itnm+2) = itriton

*-----------------------------------------------------------------------
*        iaonucl  T.Sato 2025/02/12
*-----------------------------------------------------------------------

                  itdci(51,itnm+2) = iaonucl

*-----------------------------------------------------------------------
*        aonucl  T.Sato 2025/02/12
*-----------------------------------------------------------------------

                  itdci(52,itnm+2) = laonucl
                  if(laonucl.ne.0) ctdci(1,itnm+2) = aonucl

*-----------------------------------------------------------------------
*        aoreg  T.Sato 2025/02/12
*-----------------------------------------------------------------------

                  itdci(53,itnm+2) = laoreg
                  if(laoreg.ne.0) ctdci(2,itnm+2) = aoreg

*-----------------------------------------------------------------------
*        Card[10]
*-----------------------------------------------------------------------
                  itdci(24,itnm+2) = klst
                  i0=itdc2reg(itdc) !FURUTA20200522
                  if(i0+klst.gt.maxtdcreg)
     &                 call reallocate_talldc2(i0+klst)
                  do i = 1, klst
                   j0=ireg2itm(i0+i)                 !FURUTA20200522
                   j1=jreg2itm(i)                    !FURUTA20200522
                   itgcel(i0+i)=itgcell(i)           !FURUTA20200522
                   rtgvlm(i0+i)=rtgvol(i)            !FURUTA20200522
                   itgnum(i0+i)=itglist(i)           !FURUTA20200522
                   itdcregitm=itdcregitm+itglist(i)  !FURUTA20200522
                   if(j0+itglist(i).gt.maxtdcregitm) !FURUTA20200522
     &                  call reallocate_talldc1(j0+itglist(i)) !FURUTA20200522
                   do j = 1, itglist(i)              !FURUTA20200522
                    tglst(j0+j)=htgnzas(j1+j)        !FURUTA20200522
                    tgnlt(j0+j)=ctgnnds(j1+j)        !FURUTA20200522
                   end do                            !FURUTA20200522
                  end do

*-----------------------------------------------------------------------
*        version
*-----------------------------------------------------------------------

                  itdci(25,itnm+2) = iversion

CCSE added for dchain parameter (2018.07.31) >>>>>
*-----------------------------------------------------------------------
*        mtscore
*-----------------------------------------------------------------------
            itscore(itnm  ) = -1
            itscore(itnm+1) = -1
            itscore(itnm+2) = mtscore
CCSE added for dchain parameter (2018.07.31) <<<<<

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  921    m_err = 'KeyWord "[T-Dchain]" is not found in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:5502/R:tdchain/F:tallsm4.f'
         goto 999

  922    m_err = 'imode (Card[2]-1) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5507/R:tdchain/F:tallsm4.f'
         goto 999
  923    m_err = 'jmode (Card[2]-2) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5511/R:tdchain/F:tallsm4.f'
         goto 999

  924    m_err = 'itstep (Card[3]-1) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5516/R:tdchain/F:tallsm4.f'
         goto 999
  925    m_err = 'itout (Card[3]-2) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5520/R:tdchain/F:tallsm4.f'
         goto 999
  926    m_err = 'idivs (Card[3]-3) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5524/R:tdchain/F:tallsm4.f'
         goto 999
  927    m_err = 'iregon (Card[3]-4) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5528/R:tdchain/F:tallsm4.f'
         goto 999
  928    m_err = 'inmtcf (Card[3]-5) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5532/R:tdchain/F:tallsm4.f'
         goto 999
  929    m_err = 'ichain (Card[3]-6) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5536/R:tdchain/F:tallsm4.f'
         goto 999
  930    m_err = 'itdecs (Card[3]-7) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5540/R:tdchain/F:tallsm4.f'
         goto 999
  931    m_err = 'itdecn (Card[3]-8) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5544/R:tdchain/F:tallsm4.f'
         goto 999
  932    m_err = 'isomtr (Card[3]-9) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5548/R:tdchain/F:tallsm4.f'
         goto 999
  933    m_err = 'ifisyd (Card[3]-10) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5552/R:tdchain/F:tallsm4.f'
         goto 999
  934    m_err = 'ifisye (Card[3]-11) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5556/R:tdchain/F:tallsm4.f'
         goto 999

  935    m_err = 'iyild (Card[4]-1) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5561/R:tdchain/F:tallsm4.f'
         goto 999
  936    m_err = 'iggrp (Card[4]-2) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5565/R:tdchain/F:tallsm4.f'
         goto 999
  937    m_err = 'ibetap (Card[4]-3) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5569/R:tdchain/F:tallsm4.f'
         goto 999
  938    m_err = 'acmin (Card[4]-4) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5573/R:tdchain/F:tallsm4.f'
         goto 999
  939    m_err = 'istabl (Card[4]-5) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5577/R:tdchain/F:tallsm4.f'
         goto 999
  940    m_err = 'igsdef (Card[4]-6) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5581/R:tdchain/F:tallsm4.f'
         goto 999
  941    m_err = 'iprtb1 (Card[4]-7) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5585/R:tdchain/F:tallsm4.f'
         goto 999
  942    m_err = 'iprtb2 (Card[4]-8) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5589/R:tdchain/F:tallsm4.f'
         goto 999
  943    m_err = 'rprtb2 (Card[4]-9) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5593/R:tdchain/F:tallsm4.f'
         goto 999
  944    m_err = 'iprtb3 (Card[4]-10) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5597/R:tdchain/F:tallsm4.f'
         goto 999
  945    m_err = 'igsorg (Card[4]-11) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5601/R:tdchain/F:tallsm4.f'
         goto 999

  946    m_err = 'amp (Card[5]-1) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5606/R:tdchain/F:tallsm4.f'
         goto 999
  947    m_err = 'ebeam (Card[5]-2) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5610/R:tdchain/F:tallsm4.f'
         goto 999
  948    m_err = 'prodnp (Card[5]-3) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5614/R:tdchain/F:tallsm4.f'
         goto 999

  949    m_err = 'ac-list (Card[8]) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5619/R:tdchain/F:tallsm4.f'
         goto 999
  920    m_err = 'tg-list (Card[10]) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5623/R:tdchain/F:tallsm4.f'
         goto 999


  950    write(dkam,'(i9)') mdas
         m_err = 'Total tally storage number exceeds mdas ='//dkam
         ErrCha = ''
         ErrID = 'L:5630/R:tdchain/F:tallsm4.f'
         goto 999
  951    m_err = 'mesh selection, allowed "1", is wrong.'
         ErrCha = ''
         ErrID = 'L:5634/R:tdchain/F:tallsm4.f'
         goto 999
  952    m_err = 'mesh description of Card[6] is wrong.'
         ErrCha = ''
         ErrID = 'L:5638/R:tdchain/F:tallsm4.f'
         goto 999
  953    m_err = 'mesh description of Card[7] is wrong.'
         ErrCha = ''
         ErrID = 'L:5642/R:tdchain/F:tallsm4.f'
         goto 999
  954    m_err = 'description or number of Card[8] is wrong.'
         ErrCha = ''
         ErrID = 'L:5646/R:tdchain/F:tallsm4.f'
         goto 999
  955    m_err = 'description or number of Card[10] is wrong.'
         ErrCha = ''
         ErrID = 'L:5650/R:tdchain/F:tallsm4.f'
         goto 999
  956    m_err = 'description of density (Card[10]) is wrong.'
         ErrCha = ''
         ErrID = 'L:5654/R:tdchain/F:tallsm4.f'
         goto 999
  957    m_err = 'description of version OLD(0)/NEW(1) is wrong.'
         ErrCha = ''
         ErrID = 'L:5658/R:tdchain/F:tallsm4.f'
         goto 999

  961    m_err = 'KeyWord "[T-Track]" is not found in tally [t-track]'
         ErrCha = ''
         ErrID = 'L:5663/R:tdchain/F:tallsm4.f'
         goto 999
  962    m_err = 'KeyWord is not found in tally [t-track]'
         ErrCha = ''
         ErrID = 'L:5667/R:tdchain/F:tallsm4.f'
         goto 999

  963    m_err = '*.dyld, *.dtrk or *.dout can not use as filename.'
         ErrCha = ''
         ErrID = 'L:5672/R:tdchain/F:tallsm4.f'
         goto 999

  975    m_err = '1H or Z=A cannot be specified for nucleus in tally '
     &            //tnamed
         ErrCha = ''
         ErrID = 'L:5678/R:tdchain/F:tallsm4.f'
         goto 999

  976    m_err = 'counter should be from -9999 to 9999'
         ErrCha = ''
         ErrID = 'L:5683/R:tdchain/F:tallsm4.f'
         goto 999

  978    m_err = 'ndata should be 0, 1, 2 or 3 in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:5688/R:tdchain/F:tallsm4.f'
         goto 999

  979    m_err = 'Description of material parameter is wrong in tally '
     &            //tnamed
         ErrCha = ''
         ErrID = 'L:5694/R:tdchain/F:tallsm4.f'
         goto 999

  980    m_err = 'unit should be 1 or 2 in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:5699/R:tdchain/F:tallsm4.f'
         goto 999

  981    m_err = 'Description of nucleus parameter is wrong in tally '
     &            //tnamed
         ErrCha = ''
         ErrID = 'L:5705/R:tdchain/F:tallsm4.f'
         goto 999

  982    m_err = 'Description of mother parameter is wrong in tally '
     &            //tnamed
         ErrCha = ''
         ErrID = 'L:5711/R:tdchain/F:tallsm4.f'
         goto 999

  983    m_err = '2D-type should be 1-7 in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:5716/R:tdchain/F:tallsm4.f'
         goto 999

  986    write(dkam,'(i5)') itlmax
         m_err = 'Total tally number exceeds itlmax = '//dkam
         ErrCha = ''
         ErrID = 'L:5722/R:tdchain/F:tallsm4.f'
         goto 999

  987    m_err = 'Unknown parameter in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:5727/R:tdchain/F:tallsm4.f'
         goto 999

  988    m_err = 'special should be positive in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:5732/R:tdchain/F:tallsm4.f'
         goto 999

  989    m_err = 'Double definition of the parameter in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:5737/R:tdchain/F:tallsm4.f'
         goto 999

  990    m_err = 'Only one axis in '//tnamed
         ErrCha = ''
         ErrID = 'L:5742/R:tdchain/F:tallsm4.f'
         goto 999

  991    m_err = 'Only one axis in '//tnamed
         ErrCha = ''
         ErrID = 'L:5747/R:tdchain/F:tallsm4.f'
         goto 999

  992    m_err = 'Unknown axis name in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:5752/R:tdchain/F:tallsm4.f'
         goto 999

  993    m_err = 'In this line, [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:5757/R:tdchain/F:tallsm4.f'
         goto 999

  994    m_err = 'Name of particle is wrong in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:5762/R:tdchain/F:tallsm4.f'
         goto 999

  995    m_err = 'Number of particles is larger than 6 in tally '//
     &   tnamed
         ErrCha = ''
         ErrID = 'L:5768/R:tdchain/F:tallsm4.f'
         goto 999

  996    m_err = 'Description of mesh is wrong in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:5773/R:tdchain/F:tallsm4.f'
         goto 999

  997    m_err = 'Description of parameter is wrong in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:5778/R:tdchain/F:tallsm4.f'
         goto 999

  998    m_err = 'Unknown mesh parameter in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:5783/R:tdchain/F:tallsm4.f'
         goto 999

  821    m_err = 'iertdcho (Card[9]-5) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5788/R:tdchain/F:tallsm4.f'
         goto 999

  822    m_err = 'inxslib (Card[3b]-1) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5793/R:tdchain/F:tallsm4.f'
         goto 999

  823    m_err = 'idcylib (Card[3b]-2) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5798/R:tdchain/F:tallsm4.f'
         goto 999

  824    m_err = 'infylib (Card[3b]-3) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5803/R:tdchain/F:tallsm4.f'
         goto 999

  825    m_err = 'hnxslib (Card[3b]-4) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5808/R:tdchain/F:tallsm4.f'
         goto 999

  826    m_err = 'hdcylib (Card[3b]-5) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5813/R:tdchain/F:tallsm4.f'
         goto 999

  827    m_err = 'hnfylib (Card[3b]-6) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5818/R:tdchain/F:tallsm4.f'
         goto 999

  828    m_err = 'iwrtchn (Card[4]-12) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5823/R:tdchain/F:tallsm4.f'
         goto 999

  829    m_err = 'chrlvth (Card[4]-13) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5828/R:tdchain/F:tallsm4.f'
         goto 999

 2830    m_err = 'iwrchdt (Card[4]-14) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5833/R:tdchain/F:tallsm4.f'
         goto 999

  831    m_err = 'iwrchss (Card[4]-15) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5838/R:tdchain/F:tallsm4.f'
         goto 999

  832    m_err = 'ixsrall (Card[3]-7) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5843/R:tdchain/F:tallsm4.f'
         goto 999

  833    m_err = 'idosecf (Card[4]-16) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5848/R:tdchain/F:tallsm4.f'
         goto 999

  834    m_err = 'ipltmode (Card[9d]-1) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5853/R:tdchain/F:tallsm4.f'
         goto 999

  835    m_err = 'ipltaxis (Card[9d]-2) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5858/R:tdchain/F:tallsm4.f'
         goto 999

  836    m_err = 'foamout (Card[9d]-1) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5863/R:tdchain/F:tallsm4.f'
         goto 999

  837    m_err = 'foamvals (Card[9d]-2) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5868/R:tdchain/F:tallsm4.f'
         goto 999

  838    m_err = 'iredufmt (Card[8]-7) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5873/R:tdchain/F:tallsm4.f'
         goto 999

  839    m_err = 'irdonce  (Card[8]-8) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5878/R:tdchain/F:tallsm4.f'
         goto 999

  840    m_err = 'imtcard  (Card[4]-17) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5883/R:tdchain/F:tallsm4.f'
         goto 999

  841    m_err = 'imtcnum  (Card[4]-18) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5888/R:tdchain/F:tallsm4.f'
         goto 999

  842    m_err = 'imtcmeta  (Card[4]-19) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5893/R:tdchain/F:tallsm4.f'
         goto 999

  843    m_err = 'thmatnd  (Card[4]-20) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5898/R:tdchain/F:tallsm4.f'
         goto 999

  844    m_err = 'idosunit  (Card[4]-16b) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5903/R:tdchain/F:tallsm4.f'
         goto 999

  845    m_err = 'iangbpwr  (Card[10]-2) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5908/R:tdchain/F:tallsm4.f'
         goto 999

  846    m_err = 'iphtout   (Card[10]-1) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5913/R:tdchain/F:tallsm4.f'
         goto 999

  847    m_err = 'isfylib   (Card[3b]-7) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5918/R:tdchain/F:tallsm4.f'
         goto 999

  848    m_err = 'radlev    (Card[10]-3) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5923/R:tdchain/F:tallsm4.f'
         goto 999

  849    m_err = 'ergthfa   (Card[3b]-8) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5928/R:tdchain/F:tallsm4.f'
         goto 999

 1850    m_err = 'ergfafu   (Card[3b]-9) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5933/R:tdchain/F:tallsm4.f'
         goto 999

  851    m_err = 'infyerg   (Card[3b]-10) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5938/R:tdchain/F:tallsm4.f'
         goto 999

  852    m_err = 'ilchain   (Card[3]-12) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5943/R:tdchain/F:tallsm4.f'
         goto 999

  853    m_err = 'itriton value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5948/R:tdchain/F:tallsm4.f'
         goto 999

  854    m_err = 'iaonucl value (<14) is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:5953/R:tdchain/F:tallsm4.f'
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
      subroutine tregdc(icc,jsn,jsi,dsin,idsi,ill,ilf,
     &                   jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                   ntrn,mtrn,ndsm,nvol,ivl,irvl,icfl,
     &                   idcrg,ndim_region,idas_region)
*                                                                      *
*       read 'reg =' sub-section of input tally section                *
*       modified by K.Niita on 2011/02/03                              *
*       modified by N.Matsuda on 2013/07/10                            *
*       original sub. is "tregion" in tallsm1.f                        *
*                                                                      *
************************************************************************
      use moddas
      use moddas_region_mtrg

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200

      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

      dimension ipar(0:20)
      dimension jpar(0:20)
      dimension lpar(0:20)

      data klnmax /1000000/

      integer, intent(in) :: ndim_region
      integer, intent(inout) :: idas_region(ndim_region)

*-----------------------------------------------------------------------

            ierr = 0
            idcrg = 0

      if( icfl .eq. 1 ) goto 150

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

               idcrg = idcrg + 1
               if(idcrg.gt.maxreg)call reallocate_talldc4(idcrg) !FURUTA20200522
               dcrgm(idcrg) = chin
               icch = inumc(chin,i1,i3,'(')

  141          if( chin(icch:icch) .eq. '(' ) then
                 icci = inumc(chin,icch,i3,')')
                 iccj = inumc(chin,icch,icci,'<')
                 icck = inumc(chin,icch,icci,'[')
*    &           icch, icci, iccj, icck
*    &           i1, i2, i3, i4

                   if( iccj .lt. icci ) then
                   else
                       goto 950
                   end if

                 icch = icci

                 if( icch .ge. i3 ) then
                 else
                   icch = inumc(chin,icch+1,i3,'(')
                   goto 141
                 end if

               end if

  150 continue

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

*-----------------------------------------------------------------------
*        error not for reg =
*-----------------------------------------------------------------------

         if( icc .eq. 1 ) then

            if( chcm(i1:i1+3) .ne. 'reg=' ) goto 998

         else if( icc .eq. 2 ) then

            if( chcm(i1:i1+8) .ne. 'reginbox=' ) goto 998

         end if

               ic = inumc(chlw,i1,i3,'=') + 1

*-----------------------------------------------------------------------
*        read region number
*-----------------------------------------------------------------------

               ic = jnumc(chlw,ic,i3)

               ntrn = 0
               mtrn = 0
               ibra = 0
               ilat = 0
               klat = 0
               ifis = 0
               kpar = 0
               iuni = 0

               nvol = 0

               do i = 0, 20
                  ipar(i) = 0
                  jpar(i) = 0
                  lpar(i) = 0
               end do

               igm = 0
               call moddas_allocate_int2(0, 20, MAX_NUM_MTRG, mtrg)
               ipmax = 0

            do i = 1, klnmax

                  ic = jnumc(chlw,ic,i3)

*    &                  chlw(ic:ic+2)

                  call tregion3(icn,chlw,i1,i3,ic,
     &                          igm,ipmax,ipar,jpar,
     &                          ibra,ilat,klat,ifis,kpar,iuni
     &                          ,MAX_NUM_MTRG,mtrg)

*    &                  i3,i3,i3,i3,i3)') ilat, klat, ifis, kpar, iuni

                     if( icn .gt. 900 ) goto 900
                     if( icn .eq. 500 ) goto 500

               if( i .lt. klnmax .and. ic .gt. i3 ) then

  147             call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) goto 500

                     if( iskip .ne. 0 ) goto 147

                     ifis = ifis + 1

                     ic = i1

*    &                            ill(jsn), chin(1:i3)

                     if( chin(i1:i1) .lt. 'A'     ! MATSUDA 2015.10.09
     &                   .or. chin(i1:i1) .eq. '['
     &                   .or. chin(i1:i1) .eq. ']'
     &                   .or. chin(i1:i1) .eq. '{'
     &                   .or. chin(i1:i1) .eq. '}' ) then
                        idcrg = idcrg + 1
               if(idcrg.gt.maxreg)call reallocate_talldc4(idcrg) !FURUTA20200522
                        dcrgm(idcrg) = chin
*    &                            ill(jsn), chin(1:i3)
                        icch = inumc(chin,i1,i3,'(')

  148                  if( chin(icch:icch) .eq. '(' ) then
                         icci = inumc(chin,icch,i3,')')
                         iccj = inumc(chin,icch,icci,'<')
                         icck = inumc(chin,icch,icci,'[')
                           if( iccj .lt. icci ) then
                           else
                               goto 950
                           end if

                         icch = icci

                         if( icch .ge. i3 ) then
                         else
                           icch = inumc(chin,icch+1,i3,'(')
                           goto 148
                         end if

                       end if

                     end if

               end if

            end do

                  goto 996

*-----------------------------------------------------------------------
*        modify the input
*-----------------------------------------------------------------------

  500    continue

               ntrn = ipar(0)
               mtrn = ipar(0)

            do k = 1, mtrn
               idas_region( ndsm + k - 1 ) = mtrg(igm+0,k)
            end do

            call tregion4(ierr,ntrn,mtrn,idas_region(ndsm),ndsm,igm
     &                    ,ndim_region,idas_region)

               if( ierr .ne. 0 ) goto 983

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*        volume
*-----------------------------------------------------------------------

            if( chcm(i1:i1+5) .eq. 'volume' .or.
     &          chcm(i1:i1+4) .eq. 'value' ) then

                  call tregvol(jsn,jsi,dsin,idsi,ill,ilf,
     &                         jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                         nvol,ivl,irvl)

            end if

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

         return

  900 continue

         if( icn .eq. 980 ) goto 980
         if( icn .eq. 984 ) goto 984
         if( icn .eq. 985 ) goto 985
         if( icn .eq. 986 ) goto 986
         if( icn .eq. 987 ) goto 987
         if( icn .eq. 988 ) goto 988
         if( icn .eq. 989 ) goto 989
         if( icn .eq. 990 ) goto 990
         if( icn .eq. 991 ) goto 991
         if( icn .eq. 992 ) goto 992
         if( icn .eq. 993 ) goto 993
         if( icn .eq. 995 ) goto 995
         if( icn .eq. 997 ) goto 997

*-----------------------------------------------------------------------
  950 continue
         m_err = 'combined cell can not use'
         ErrCha = ''
         ErrID = 'L:6243/R:tregdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  982 continue
         m_err = 'volume of value section is wrong'
         ErrCha = ''
         ErrID = 'L:6254/R:tregdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  983 continue
         m_err = 'region is too many or memory is lack'
         ErrCha = ''
         ErrID = 'L:6265/R:tregdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  980 continue

         m_err = 'n1-n2 should be used like {n1-n2}'
         ErrCha = ''
         ErrID = 'L:6277/R:tregdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  984 continue

         m_err = 'maximum lattice elements by commas is 1000.'
         ErrCha = ''
         ErrID = 'L:6289/R:tregdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  985 continue

         m_err = 'maximum level of < is 10.'
         ErrCha = ''
         ErrID = 'L:6301/R:tregdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  986 continue

         m_err = 'usage of u=# is wrong.'
         ErrCha = ''
         ErrID = 'L:6313/R:tregdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  987 continue

         m_err = '< should be inside ( ).'
         ErrCha = ''
         ErrID = 'L:6325/R:tregdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  988 continue

         m_err = 'maximun level of ( ) is 10.'
         ErrCha = ''
         ErrID = 'L:6337/R:tregdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  989 continue

         m_err = 'usage of latice [i1:i2 i3:i4 i5:i6] is wrong.'
         ErrCha = ''
         ErrID = 'L:6349/R:tregdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  990 continue

         m_err = 'usage of latice [i1 i2 i3] is wrong.'
         ErrCha = ''
         ErrID = 'L:6361/R:tregdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  991 continue

         m_err = 'all or (all) is available. (all<4), (all 3) are not.'
         ErrCha = ''
         ErrID = 'L:6373/R:tregdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         m_err = 'Usage of Parenthesis { n1 - n2 } is wrong'
         ErrCha = ''
         ErrID = 'L:6385/R:tregdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'Usage of Parenthesis ( ) is wrong'
         ErrCha = ''
         ErrID = 'L:6397/R:tregdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Description of region number '//
     &           '{n1-n2} (n1<n2) is wrong.'
         ErrCha = ''
         ErrID = 'L:6410/R:tregdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Number of region is too large. max(klnmax)=1000000'
         ErrCha = ''
         ErrID = 'L:6422/R:tregdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of region number is wrong.'
         ErrCha = ''
         ErrID = 'L:6434/R:tregdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'After { mesh = reg } line should be { reg = }.'
         ErrCha = ''
         ErrID = 'L:6446/R:tregdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end subroutine


************************************************************************
*                                                                      *
      subroutine gettev(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                  igtp,igr,temsh,itemsh)
*                                                                      *
*       read irradiation/cooling time and beam power information       *
*       modified by N.Matsuda on 2012/06/30                            *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------



      dimension temsh(4,maxitm), itemsh(2,maxitm)
      character techa*1

*-----------------------------------------------------------------------
      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

            ierr  = 0
            icmp = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 800

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of the section
*-----------------------------------------------------------------------

            if( chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 800

            end if

*-----------------------------------------------------------------------
*        read mesh points
*-----------------------------------------------------------------------

            ic = i1

            do i = 1, igr

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               temsh(2*igtp-1,i) = cvvv

                  ic = ic2
                  ic = jnumc(chlw,ic,i3)

                  techa = chlw(ic:ic)
                  if( techa .eq. 'y' ) then
                     itemsh(igtp,i) = 1
                  else if( techa .eq. 'd' ) then
                     itemsh(igtp,i) = 2
                  else if( techa .eq. 'h' ) then
                     itemsh(igtp,i) = 3
                  else if( techa .eq. 'm' ) then
                     itemsh(igtp,i) = 4
                  else if( techa .eq. 's' ) then
                     itemsh(igtp,i) = 5
                  else
                     goto 995
                  end if

                  ic = ic + 1

               if( igtp .eq. 1 ) then
                 ic = jnumc(chlw,ic,i3)
                 call snum(chlw,ic,i3,ic2,cvvv,ierr)
                    temsh(2*igtp,i) = cvvv
                    ic = ic2

               else if( igtp .eq. 2 ) then
                    temsh(2*igtp,i) = 0.0
               end if

               if( i .lt. igr .and. ic .gt. i3 ) then

  149             call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                  if( ierr .ne. 0 ) return
                  if( jpn  .eq. 3 ) goto 997

                  if( iskip .ne. 0 ) goto 149

                  ic = i1

               end if

            end do

            icmp = 1

*-----------------------------------------------------------------------
*     summary and check
*-----------------------------------------------------------------------

  800 continue

            if( igtp .eq. 1 ) then
               if( igr  .le. 0 .or.
     &             icmp .eq. 0 ) goto 998
            else if( igtp .eq. 2 ) then
               if( igr  .le. 0 .or.
     &             icmp .eq. 0 ) goto 998
            end if

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  995 continue

         m_err = 'mesh parameter (required) is wrong'
         ErrCha = ''
         ErrID = 'L:6610/R:gettev/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'In this line, [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:6622/R:gettev/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of mesh points is wrong'
         ErrCha = ''
         ErrID = 'L:6634/R:gettev/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Description of mesh parameter is wrong'
         ErrCha = ''
         ErrID = 'L:6646/R:gettev/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Mesh type should be 1'
         ErrCha = ''
         ErrID = 'L:6658/R:gettev/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end subroutine


************************************************************************
*                                                                      *
      subroutine tdchech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the dchain tally                                 *
*       last modified by T.Furuta on 2025/01/17                        *
*                                                                      *
************************************************************************

cMIURA(2016.09.30) chg and active >>>
      use sumtallymod
cMIURA(2016.09.30) chg and active <<<
      use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      include 'err.inc'
      data tdcheck_CallFlag / 0 /
*-----------------------------------------------------------------------
*     for material
      common /inggs/  iog, igcel, ioa, igsuf, iob, igtrs
      common /regdm/  idmg(kvlmax)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /kmat1a/ mxmat, mxmat0, mxnel
      common /kmat1i/ kmatc(kvlmax)

*     original
      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

*     mesh
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

*     mesh = 2 and 3
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

*     mother, nucleus,
      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)

      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax),
     &                ihxll(itlmax), ihxle(itlmax),
     &                ihdll(itlmax), ihdle(itlmax),
     &                ihfll(itlmax), ihfle(itlmax),
     &                ihsll(itlmax), ihsle(itlmax)
      character ittle*80, ihxle*80, ihdle*80, ihfle*80, ihsle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall21/ rtfac(itlmax)

*     material
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall30/ itnda(itlmax)
      common /tall31/ iterl(itlmax)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall37/ itcnt(9,itlmax)

      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)
      common /tall52/ itprd(itlmax)


      common /tall67/ itstd(itlmax), rtstd(itlmax)

*     for dchain
      common /talldc/ itdci(nientry,itlmax), rtdci(nrentry,itlmax),
     &               tbina(itlmax,maxitm), tbinb(itlmax,maxitm),
     &               beampw(itlmax,maxitm),
     &               tmina(itlmax,maxitm), tminb(itlmax,maxitm),
     &               aclst(itlmax,maxitm), bqlst(itlmax,maxitm)
      common /talldc3/ ctdci(ncentry,itlmax) ! T.Sato, character variable in [t-dchain]
      character ctdci*200

CCSE added for dchain parameter, mtscore (2018.07.31) >>>>>
      common /tall77/ itscore(itlmax)
CCSE added for dchain parameter, mtscore (2018.07.31) <<<<<
      common /tall83/ itnzn(itlmax), itndm(itlmax)

      character techa
      character aclst*8
      character chaac(1:8)*8
! T.Sato 2015/8/25 Not used?
      character bq*1

      dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

*-----------------------------------------------------------------------
       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr
!-----------------------------------------------------------------------

      character asfil*100
      character esfil*100
      character erfnm*100

      character aname(20)*7
      data aname / 'mass   ',' reg   ','   x   ','   y   ',
     &             '   z   ','   r   ','charge ','chart  ',
     &             '  xy   ','  yz   ','  xz   ','  rz   ',
     &             'dchain ',' all   ',' cla   ',' gsd   ',
     &             ' gso   ',' tba   ',' tbh   ',' yld   '/

      dimension ian(6)
      dimension jmat(20)
      character chau*8
      character chaus(7)*8
      character elmnt(104)*3

*-----------------------------------------------------------------------

      data elmnt /
     &    ' H ','He ','Li ','Be ',' B ',' C ',' N ',' O ',' F ','Ne ',
     &    'Na ','Mg ','Al ','Si ',' P ',' S ','Cl ','Ar ',' K ','Ca ',
     &    'Sc ','Ti ',' V ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ',' Y ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ',' I ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ',' W ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ',' U ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------
            write(iot,'("[ T-Dchain ]")')
*-----------------------------------------------------------------------
*           Card[1] title (free format)
                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )
*-----------------------------------------------------------------------
*           Card[2]
                  write(iot,'("    imode =",i5,10x,
     &            "  # Card[2]-1 control param.")')
     &            itdci(1,m)

                  write(iot,'("    jmode =",i5,10x,
     &            "  # Card[2]-2 control param.")')
     &            itdci(2,m)
*-----------------------------------------------------------------------
*           Card[3]
                  write(iot,'("   itstep =",i5,10x,
     &            "  # (D=1) Card[3]-1 calc. condition")')
     &            itdci(3,m)

                  write(iot,'("    itout =",i5,10x,
     &            "  # (D=1) Card[3]-2 calc. condition")')
     &            itdci(4,m)

                  write(iot,'("    idivs =",i5,10x,
     &            "  # (D=50) Card[3]-3 calc. condition")')
     &            itdci(5,m)

                  write(iot,'("   iregon =",i5,10x,
     &            "  # (D=1) Card[3]-4 calc. condition")')
     &            itdci(6,m)

                  write(iot,'("   inmtcf =",i5,10x,
     &            "  # (D=1) Card[3]-5 calc. condition")')
     &            itdci(7,m)
*---
                  write(iot,'("   ichain =",i5,10x,
     &            "  # (D=100) Card[3]-6 calc. condition")')
     &            itdci(8,m)

                  write(iot,'("   itdecs =",i5,10x,
     &            "  # (D=0) Card[3]-7 calc. condition")')
     &            itdci(9,m)

                  write(iot,'("   itdecn =",i5,10x,
     &            "  # (D=1) Card[3]-8 calc. condition")')
     &            itdci(10,m)

                  write(iot,'("   isomtr =",i5,10x,
     &            "  # (D=2) Card[3]-9 calc. condition")')
     &            itdci(11,m)

                  write(iot,'("   ifisyd =",i5,10x,
     &            "  # (D=1) Card[3]-10 calc. condition")')
     &            itdci(12,m)
*---
                  write(iot,'("   ifisye =",i5,10x,
     &            "  # (D=0) Card[3]-11 calc. condition")')
     &            itdci(13,m)

c H.Ratliff START 2020.04.17 add parameter controlling reading/loading of xs lib
                  write(iot,'("  ixsrall =",i5,10x,
     &            "  # (D=1) Card[3]-12 calc. condition")')
     &            itdci(33,m)


c H.Ratliff START 2019.07.30 add library selection parameters (fission values commented out)
*-----------------------------------------------------------------------
*           Card[3b]
                  write(iot,'("  inxslib =",i5,10x,
     &            "  # (D=2) Card[3b]-1 data lib opt.")')
     &            itdci(27,m)

                  write(iot,'("  idcylib =",i5,10x,
     &            "  # (D=5) Card[3b]-2 data lib opt.")')
     &            itdci(28,m)

                  write(iot,'("  infylib =",i5,10x,
     &            "  # (D=0) Card[3b]-3 data lib opt.")')
     &            itdci(29,m)

                  if (itdci(27,m).eq.-1)
     &            write(iot,'("  hnxslib = ",80a1,
     &            "  # (D=0) Card[3b]-4 data lib opt.")')
     &            ( ihxle(m)(i:i), i = 1, ihxll(m) )

                  if (itdci(28,m).eq.-1)
     &            write(iot,'("  hdcylib = ",80a1,
     &            "  # (D=0) Card[3b]-5 data lib opt.")')
     &            ( ihdle(m)(i:i), i = 1, ihdll(m) )

                  if (itdci(29,m).eq.-1)
     &            write(iot,'("  hnfylib = ",80a1,
     &            "  # (D=0) Card[3b]-6 data lib opt.")')
     &            ( ihfle(m)(i:i), i = 1, ihfll(m) )
*-----------------------------------------------------------------------
*           Card[4]
                  write(iot,'("    iyild =",i5,10x,
     &            "  # (D=0) Card[4]-1 output opt.")')
     &            itdci(14,m)

                  write(iot,'("    iggrp =",i5,10x,
     &            "  # (D=3) Card[4]-2 output opt.")')
     &            itdci(15,m)

                  write(iot,'("   ibetap =",i5,10x,
     &            "  # (D=1) Card[4]-3 output opt.")')
     &            itdci(16,m)

                  write(iot,'("    acmin =",1p1g14.7,1x,1x,
     &            " # (D=0.0) Card[4]-4 output opt.")')
     &            rtdci(1,m)

                  write(iot,'("   istabl =",i5,10x,
     &            "  # (D=0) Card[4]-5 output opt.")')
     &            itdci(17,m)
*---
                  write(iot,'("   igsdef =",i5,10x,
     &            "  # (D=1) Card[4]-6 output opt.")')
     &            itdci(18,m)

                  write(iot,'("C  iprtb1 =",i5,10x,
     &            "  # (D=1) Card[4]-7 output opt. -DISABLE!")')
     &            itdci(19,m)

                  write(iot,'("C  iprtb2 =",i5,10x,
     &            "  # (D=1) Card[4]-8 output opt. -DISABLE!")')
     &            itdci(20,m)

                  write(iot,'("C  rprtb2 =",1p1g14.7,1x,1x,
     &            " # (D=10.) Card[4]-9 output opt. -DISABLE!")')
     &            rtdci(2,m)

                  write(iot,'("C  iprtb3 =",i5,10x,
     &            "  # (D=0) Card[4]-10 output opt. -DISABLE!")')
     &            itdci(21,m)
*---
                  write(iot,'("   igsorg =",i5,10x,
     &            "  # (D=1) Card[4]-11 output opt.")')
     &            itdci(22,m)

c H.Ratliff START 2019.07.30 add decay chain output file options
                  write(iot,'("  iwrtchn =",i5,10x,
     &            "  # (D=1) Card[4]-12 output opt.")')
     &            itdci(30,m)

                  write(iot,'("  chrlvth =",1p1g14.7,1x,
     &            "  # (D=-1) Card[4]-13 output opt.")')
     &            rtdci(6,m)

                  write(iot,'("  iwrchdt =",i5,10x,
     &            "  # (D=0) Card[4]-14 output opt.")')
     &            itdci(31,m)

                  write(iot,'("  iwrchss =",i5,10x,
     &            "  # (D=0) Card[4]-15 output opt.")')
     &            itdci(32,m)

c H.Ratliff START 2020.04.17 add control of dose rate coefficient

                  write(iot,'("  idosecf =",i5,10x,
     &            "  # (D=1) Card[4]-16 calc. condition")')
     &            itdci(34,m)

                  write(iot,'(" ipltmode =",i5,10x,
     &            "  # (D=0) Card[9d]-1 xyz mesh plot condition")')
     &            itdci(35,m)

                  write(iot,'(" ipltaxis =",i5,10x,
     &            "  # (D=1) Card[9d]-2 xyz mesh plot condition")')
     &            itdci(36,m)

c H.Ratliff START 2020.05.29 add tet mesh extra output options

                  write(iot,'("  foamout =",i5,10x,
     &            "  # (D=0) Card[9d]-1 tet mesh OpenFoam option")')
     &            itdci(37,m)

                  write(iot,'(" foamvals =",i5,10x,
     &            "  # (D=1) Card[9d]-2 tet mesh OpenFoam option")')
     &            itdci(38,m)

c H.Ratliff START 2020.06.05

                  write(iot,'(" iredufmt =",i5,10x,
     &            "  # (D=1) Card[8]-7 dtrk/dyld formatting flag")')
     &            itdci(39,m)

c H.Ratliff START 2020.06.17

                  write(iot,'(" irdonce  =",i5,10x,
     &            "  # (D=1) Card[8]-8 dtrk/dyld parsing flag")')
     &            itdci(40,m)

c H.Ratliff START 2020.11.16 add material card output options

                  write(iot,'("  imtcard =",i5,10x,
     &            "  # (D=1) Card[4]-17 output opt.")')
     &            itdci(41,m)

                  write(iot,'("  imtcnum =",i5,10x,
     &            "  # (D=0) Card[4]-18 output opt.")')
     &            itdci(42,m)

                  write(iot,'(" imtcmeta =",i5,10x,
     &            "  # (D=2) Card[4]-19 output opt.")')
     &            itdci(43,m)

                  write(iot,'("  thmatnd =",1p1g14.7,1x,
     &            "  # (D=1E-6) Card[4]-20 output opt.")')
     &            rtdci(7,m)

c H.Ratliff START 2020.12.18 add dose unit and ANGEL beam power options
                  write(iot,'(" idosunit =",i5,10x,
     &            "  # (D=333) Card[4]-16b output opt.")')
     &            itdci(44,m)

                  write(iot,'(" iangbpwr =",i5,10x,
     &            "  # (D=1) Card[10]-2 output opt.")')
     &            itdci(45,m)

c H.Ratliff START 2021.02.17 add phits output control
                  write(iot,'(" iphtout  =",i5,10x,
     &            "  # (D=0) Card[10]-1 output opt.")')
     &            itdci(46,m)

c T.Sato 2025/02/08
                  if(itdci(50,m).ge.0) write(iot,'(" itriton  =",i5,10x,
     &            "  # (D=0) Card[1]-3 output opt.")')
     &            itdci(50,m)

c T.Sato 2025/02/12
                  if(itdci(51,m).ge.0) write(iot,'(" iaonucl  =",i5,10x,
     &            "  # (D=0) angelout_nuclides")')
     &            itdci(51,m)

                  if(itdci(52,m).ge.1) write(iot,'(" aonucl  = ",a)')
     &            ctdci(1,m)(1:itdci(52,m))

                  if(itdci(53,m).ge.1) write(iot,'(" aoreg   = ",a)')
     &            ctdci(2,m)(1:itdci(53,m))

c H.Ratliff START 2021.04.07
                  write(iot,'(" isfylib  =",i5,10x,
     &            "  # (D=17) Card[3b]-7 fission lib opt.")')
     &            itdci(47,m)

                  write(iot,'("   radlev =",1p1g14.7,1x,
     &            "  # (D=0.0) Card[10]-3 output opt.")')
     &            rtdci(8,m)

                  write(iot,'("  ergthfa =",1p1g14.7,1x,
     &            "  # (D=0.1) Card[3b]-8 fission opt.")')
     &            rtdci(9,m)

                  write(iot,'("  ergfafu =",1p1g14.7,1x,
     &            "  # (D=10.0) Card[3b]-9 fission opt.")')
     &            rtdci(10,m)

                  write(iot,'(" infyerg  =",i5,10x,
     &            "  # (D=4) Card[3b]-10 fission opt.")')
     &            itdci(48,m)

                  if (itdci(47,m).eq.-1)
     &            write(iot,'("  hsfylib = ",80a1,
     &            "  # (D=0) Card[3b]-11 data lib opt.")')
     &            ( ihsle(m)(i:i), i = 1, ihsll(m) )

c H.Ratliff START 2021.04.26 add separate control of chain length from number of chains
                  write(iot,'(" ilchain  =",i5,10x,
     &            "  # (D=1) Card[3]-12 calculation opt.")')
     &            itdci(49,m)

*-----------------------------------------------------------------------
*           Card[5]
                  write(iot,'("      amp =",1p1g14.7,1x,1x,
     &            " # (D=1.0) Source Intensity(source/sec)")')
     &            rtdci(3,m)
*-----------------------------------------------------------------------
*           'factor' operate similarly to 'amp'
*
*-----------------------------------------------------------------------
                  write(iot,'("    ebeam =",1p1g14.7,1x,1x,
     &            " # (D=3.0) Card[5]-2 beam data")')
     &            rtdci(4,m)

                  write(iot,'("   prodnp =",1p1g14.7,1x,1x,
     &            " # (D=1.0) Card[5]-3 beam data")')
     &            rtdci(5,m)
*-----------------------------------------------------------------------
*           NEW! version NEW(0)/OLD(1)
                  write(iot,'(" dversion =",i5,10x,
     &            "  # version of DCHAIN OLD(0)/NEW(1)")')
     &            itdci(25,m)
*-----------------------------------------------------------------------
*           Card[6]
                  write(iot,'("  timeevo =",i5,10x,
     &            "  # Card[6] time evolution")')
     &             itdci(3,m)

                  do i = 1, itdci(3,m)
                     if( tbinb(m,i) .eq. 1 ) then
                         techa = 'y'
                     else if( tbinb(m,i) .eq. 2 ) then
                         techa = 'd'
                     else if( tbinb(m,i) .eq. 3 ) then
                         techa = 'h'
                     else if( tbinb(m,i) .eq. 4 ) then
                         techa = 'm'
                     else if( tbinb(m,i) .eq. 5 ) then
                         techa = 's'
                     end if

                     write(iot,'("    ",e14.7,a2,2x,1p1g14.7)')
     &               tbina(m,i), techa, beampw(m,i)

                  end do
*-----------------------------------------------------------------------
*           Card[7]
                  write(iot,'("  outtime =",i5,10x,
     &            "  # Card[7] output time")')
     &             itdci(4,m)
C -- D-Chain outtime Check ----------------------------------------
                  tdcheck_CallFlag = tdcheck_CallFlag + 1
C -----------------------------------------------------------------
                  do i = 1, itdci(4,m)
                     if( tminb(m,i) .eq. 1 ) then
                         techa = 'y'
                     else if( tminb(m,i) .eq. 2 ) then
                         techa = 'd'
                     else if( tminb(m,i) .eq. 3 ) then
                         techa = 'h'
                     else if( tminb(m,i) .eq. 4 ) then
                         techa = 'm'
                     else if( tminb(m,i) .eq. 5 ) then
                         techa = 's'
                     end if

C -- D-Chain outtime Check ----------------------------------------
                     tmina_tmp = tmina(m,i)
                     if(tdcheck_CallFlag == 1 ) then
                      call dchain_outtime_check(tmina_tmp,techa,i)
                     endif
C -----------------------------------------------------------------

                     write(iot,'("    ",e14.7,a2)')
     &               tmina(m,i), techa

                  end do
*-----------------------------------------------------------------------
*           Card[8]
               write(iot,'("C  Card[8] is disabled function",
     &         " for NEW DCHAIN.")')
*-----------------------------------------------------------------------
*           Card[9]
*           FFNMTCF = extension ".dyld".
*           HHNMTCF = renamed file with ".dyld" extension.
                  write(iot,'(" iertdcho =",i5,10x,
     &            "  # (D=1) Card[9]-5 propagate err option")')
     &            itdci(26,m)

*C**-----------------------------------------------------------------------
*           Card[10]
*C**               if(itdci(24,m) .ne. 0) then
*C**                  write(iot,'("  tg-list =",i5,10x,
*C**     &            "  # Card[10a] target nuclide")')
*C**     &            itdci(24,m)
*C**
*C**                     k = 0
*C**                  do j = 1, itdci(24,m)
*C**                     write(iot,'(8x,a8,a20)')
*C**     &               tglst(m,j), tgnlt(m,j)
*C**                  end do
*C**
*C**               else
*-----------------------------
*C***     get the cell No. -> material No.
*C***              if( icgg .eq. 1 ) then
*C**
*C**                  do j = 1, igcel
*C***                           ^^^^^ <- number of all cell
*C***                    if( ndcrg(m,j) .eq. idrg(j) ) goto 100
*C**                     if( 106 .eq. idrg(j) ) goto 100
*C***                        ^^^ <- cell No.
*C**                  end do
*C**  100          continue
*C***                    related No."j" => material No.: idmg(j)
*C**
*C***              end if
*-----------------------------
*C***     get the material No. -> material list
*C***              if( mxmat .gt. 0 .and. ierrg .eq. 0 ) then
*C**
*C**                  do k = 1, mxmat0
*C**                     if( idmg(j) .eq. idmn(k) ) goto 110
*C***                                     ^^^^^^^ <- material No.
*C**                  end do
*C**  110          continue
*C***                    related No."k" => material list
*C**
*C**                  nel  = nint( das(kmatc(k)+1) )
*C**                  denh =       das(kmatc(k)+2)
*C**                  nelh = nel
*C**
*C**                  if( denh .ne. 0.0d0 ) then
*C**                      nelh = nelh + 1
*C**                  end if
*C**
*C**                  write(iot,'("  tg-list =",i5,10x,
*C**     &            "  # Card[10a] target nuclide")') nelh
*C**
*C**                  if( denh .ne. 0.0d0 ) then
*C**                      chau = ' H-1'
*C**                      write(iot,'(8x,a4,8x,1p1e15.7)')
*C**     &                chau(1:4), denh
*C**                  end if
*C**
*C**                  do kk = 1, nel
*C**                     iz   = nint( das(kmatc(k)+(kk-1)*4+14) )
*C**                     ia   = nint( das(kmatc(k)+(kk-1)*4+15) )
*C**                     den  =       das(kmatc(k)+(kk-1)*4+16)
*C**                     chau = elmnt(iz)
*C**
*C**                     if( ia .eq. 0 ) then
*C**                        write(iot,'(8x,a2,10x,1p1e15.7)')
*C**     &                  chau(1:2), den
*C**                     else if( ia .lt. 10 ) then
*C**                        write(iot,'(8x,a2,"-",i1,8x,1p1e15.7)')
*C**     &                  chau(1:2), ia, den
*C**                     else if( ia .lt. 100 ) then
*C**                        write(iot,'(8x,a2,"-",i2,7x,1p1e15.7)')
*C**     &                  chau(1:2), ia, den
*C**                     else
*C**                        write(iot,'(8x,a2,"-",i3,6x,1p1e15.7)')
*C**     &                  chau(1:2), ia, den
*C**                     end if
*C**                  end do
*C**               end if
*-----------------------------------------------------------------------
*           HOLD the original programs for additional information.
*-----------------------------------------------------------------------
*           add. info. for FILENAME
                  asfil = '  # file name of output for the above axis'
                  esfil = '  # file name of output for the error'

*-----------------------------------------------------------------------
*
*-----------------------------------------------------------------------
*
              if ( itfll(m,1) .gt. 0 ) then
                  msfile = max( 14, itfll(m,1) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,1)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

*   50                   itfp = j - 1
*
              end if
*-----------------------------------------------------------------------
*           add. info. for MESH
            if( itmsh(m-2) .eq. 1 ) then

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is region-wise")')
                     idas1 = mmmax
                     idas2 = ( idas1 + itrgn(m-2) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m-2)

                  call echrgdc(1,iterl(m-2),0,1,iot,itrcm(m-2),
     &                       idas_itrcg(itrcg(m-2)),itrgn(m-2),
     &                       itrgm(m-2),idas_itreg(itreg(m-2)),
     &                       das(idas1),idas(idas2),itrnv(m-2),
     &                       idas(itriv(m-2)),das(itrrv(m-2)),
     &                       idas3,m)

CCSE active for mesh=r-z/xyz parameter (2018.04.30) >>>>>
           else if( itmsh(m-2) .eq. 2 ) then

                  write(iot,'("     mesh =  r-z",11x,
     &            " # mesh type is r-z scoring mesh")')

                  write(iot,'("       x0 = ",1p1g14.7,1x,
     &            " # center x-position of r-z mesh")') rtrx0(m-2)

                  write(iot,'("       y0 = ",1p1g14.7,1x,
     &            " # center y-position of r-z mesh")') rtry0(m-2)

                  call echmty(0,iot,'r',itrty(m-2),itrnm(m-2),
     &                      rtrdl(m-2),rtrmi(m-2),rtrma(m-2),itrrg(m-2)
     &                      ,abs(itrnm(m)+1),das_itrrg)
                  call echmty(0,iot,'z',itzty(m-2),itznm(m-2),
     &                      rtzdl(m-2),rtzmi(m-2),rtzma(m-2),itzrg(m-2)
     &                      ,abs(itznm(m)+1),das_itzrg)

           else if( itmsh(m-2) .eq. 3 ) then

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

                  call echmty(0,iot,'x',itxty(m-2),itxnm(m-2),
     &                      rtxdl(m-2),rtxmi(m-2),rtxma(m-2),itxrg(m-2)
     &                      ,abs(itxnm(m-2)+1),das_itxrg)
                  call echmty(0,iot,'y',ityty(m-2),itynm(m-2),
     &                      rtydl(m-2),rtymi(m-2),rtyma(m-2),ityrg(m-2)
     &                      ,abs(itynm(m-2)+1),das_ityrg)
                  call echmty(0,iot,'z',itzty(m-2),itznm(m-2),
     &                      rtzdl(m-2),rtzmi(m-2),rtzma(m-2),itzrg(m-2)
     &                      ,abs(itznm(m-2)+1),das_itzrg)

CCSE active for mesh=r-z/xyz parameter (2018.07.30) >>>>>
                  if( icax .ne. 0 ) then
                     ! to echo mesh part and material part in phits.out and ***.dout
                     call echxydc(iot,itxnm(m-2),itynm(m-2),itznm(m-2),
     &                            itxrg(m-2),ityrg(m-2),itzrg(m-2) ,
     &                            itscore(m),m-2 ) !FURUTA20200615
                  endif
CCSE active for mesh=r-z/xyz parameter (2018.07.30) <<<<<
CCSE active for mesh=r-z/xyz parameter (2018.04.30) <<<<<
            elseif( itmsh(m-2) .eq. 4 ) then
                  write(iot,'("     mesh =  tet",11x,
     &            " # mesh type is tetra scoring mesh ")')

                  call echtetdc(iot,itrgn(m-2),itrgm(m-2),
     &                       idas_itreg(itreg(m-2)),
     &                       itscore(m)) !FURUTA20200601

            end if
*-----------------------------------------------------------------------
*           add. info. for TRANSFORMATION
*-----------------------------------------------------------------------
*           add. info. for MATERIAL
*-----------------------------------------------------------------------
*           add. info. for MOTHER
*-----------------------------------------------------------------------
*           add. info. for NUCLEUS
*-----------------------------------------------------------------------
*           add. info. for VOLUME
*-----------------------------------------------------------------------
*           add. info. for ECHO LENGTH
*-----------------------------------------------------------------------
               if( itnzn(m-2) .ne. 0 ) then

                  write(iot,'(" mxnuclei =",i5,10x,
     &            "  # (D=3000) maximum of product nuclei.",
     &            " 0: unlimited")') itnzn(m-2)-1

               end if

*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*        sumtally subsection
*-----------------------------------------------------------------------
            if ( icntl==13 .and. nsumtalRead(m)==16 ) then

               call sumtal_echo(iot,m)

            end if
*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------

      return
      end subroutine

************************************************************************
*                                                                      *
      subroutine pdchreg(m,mz,mn,
     &                   nr,mr,nn,kr,nt,ikzz,iknn,
     &                   tr,tm,vl,lr,nvl,ivl,rvl,
     &                   nx,ny,nz,xm,ym,zm,val,ixyz,igsh,idasa)
*                                                                      *
*       output of dchain tally in region mesh                          *
*       last modified by N.Matsuda on 2012/08/01                       *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
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

      common /tall21/ rtfac(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)

*-----------------------------------------------------------------------

      common /volreg/ dvol(kvlmax)
      common /regdc/ idrg(kvlmax), idgr(kvmmax)

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
***       common /jcomon/ nabov,nobch,nocas,nomax
*** !$OMP THREADPRIVATE(/jcomon/)
      common /talout/ itall

      character fname*100, fnume*3

*-----------------------------------------------------------------------

      dimension   kr(mr)
      dimension   nt(nn)
      dimension   ikzz(maxpt,maxnt), iknn(maxpt,maxnt)
      dimension   tm(maxnt+maxpt,2)
      dimension   vl(nr)
      dimension   lr(nr)
      dimension   tr(nr,mz,mn,2)
      dimension   ivl(nvl)
      dimension   rvl(nvl)

      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)
      dimension   val(nr)
      dimension   ixyz(1)

*-----------------------------------------------------------------------

      character hsunit(2)*15

      data hsunit / '[1/source]     ',
     &              '[1/cm^3/source]'/

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

      data ipstep / 12 /

      character erfnm*100

      character chau*8
      character cha*1
      data cha /"'"/

      character rpa*1
      data rpa /'}'/
      character yen*1
      yen  = char(92)

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do 900 iax = 1, itfln(m)

*    &     ( itaxs(m,iax) .lt. 7 .or. itrsh(m) .eq. 0 ) ) goto 900

                  fname = ctfln(m,iax)(1:itfll(m,iax))
                  ifname = itfll(m,iax)
***  "RENAME" modified by norihiro MATSUDA  2012.07.26  ----------------


                  itfp = 0
               do j = ifname, 1, -1
                  if( fname(j:j) .eq. '.' ) goto 810
               end do
                  j = ifname + 1

  810             itfp = j - 1

                  fname = fname(1:itfp)//'.dout'
                  ifname = itfp + 6

               do j = ifname, 100
                  fname(j:j) = ' '
               end do

*------------------------------------------------------- "RENAME" end --

*    &       igsh .eq. 0 ) then
         if( itall .eq. 2 .and. igsh .eq. 0 ) then
         else
         end if

            iot = 31
            open(iot, file = fname, status = 'unknown' )

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

               call tdchech(iot,m,iax,1)

*-----------------------------------------------------------------------
*        dchain axis
*-----------------------------------------------------------------------
*
*-----------------------------------------------------------------------
*        output : do not change below expresion
*-----------------------------------------------------------------------

            write(iot,'(/"#",78("-"))')

*-----------------------------------------------------------------------

            close(iot)
  900 continue

*-----------------------------------------------------------------------

      call tdctoinp(m,fname) !FURUTA20200601

*-----------------------------------------------------------------------

      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine echrgdc(icc,ierl,iva,ivs,iot,mc,kc,nr,mr,kr,
     &                 vl,lr,nvl,ivl,rvl,igm,m)
*                                                                      *
*       input echo for tally region mesh                               *
*            modified by N.MATSUDA on 2012/11/05                       *
*       last modified by H.Ratliff on 2019/07/30                       *
*                                                                      *
************************************************************************
      use moddas_material

      implicit real*8 (a-h,o-z)

      include 'param.inc'
C MATSUDA 2017.10.04 (expand natural nucleus)
      include 'param01.inc'  ! maxpt
      include 'ggsparam.inc'

      common /regdc/ idrg(kvlmax), idgr(kvmmax)

*-----------------------------------------------------------------------

      character dum1*10000
      character dum2*10000
      character chdf(500)*100
      dimension ildf(500)
      character cblan*100

      dimension kc(mc)
      dimension kr(mr)

      dimension vl(nr)
      dimension lr(nr)
      dimension ivl(nvl)
      dimension rvl(nvl)

*     for material
      common /inggs/  iog, igcel, ioa, igsuf, iob, igtrs
      common /regdm/  idmg(kvlmax)
      common /kmat1a/ mxmat, mxmat0, mxnel
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /kmat1i/ kmatc(kvlmax)
      common /talldc2/ itarget
C MATSUDA 2017.10.04 (expand natural nucleus)
      common /natnuc/ natnn(maxpt), natnm(maxpt,10), patnn(maxpt,10)

      common /talldc/ itdci(nientry,itlmax), rtdci(nrentry,itlmax),
     &               tbina(itlmax,maxitm), tbinb(itlmax,maxitm),
     &               beampw(itlmax,maxitm),
     &               tmina(itlmax,maxitm), tminb(itlmax,maxitm),
     &               aclst(itlmax,maxitm), bqlst(itlmax,maxitm)

      character aclst*8
      character chau*8
      character elmnt(104)*3

      data elmnt /
     &    ' H ','He ','Li ','Be ',' B ',' C ',' N ',' O ',' F ','Ne ',
     &    'Na ','Mg ','Al ','Si ',' P ',' S ','Cl ','Ar ',' K ','Ca ',
     &    'Sc ','Ti ',' V ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ',' Y ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ',' I ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ',' W ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ',' U ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

*-----------------------------------------------------------------------
*     read from kc(mc)
*-----------------------------------------------------------------------

            call echrg2(mc,kc,dum1,lng1,icmb,igm)

*-----------------------------------------------------------------------
*     input echo for tally region mesh
*-----------------------------------------------------------------------

                  cblan = ' '

                  ild1 = 12
                  ild0 = 72 - ild1

                  ilrm = lng1
                  isqd = 0
                  isrm = 0

  430             isqd = isqd + 1

                  if( ilrm .le. ild0 ) then

                        ildf(isqd) = ilrm
                     do k = 1, ilrm
                        chdf(isqd)(k:k) = dum1(k+isrm:k+isrm)

                     end do

                  else

                     do k = ild0, 1, -1

                        if( dum1(k+isrm:k+isrm) .eq. ' ' ) goto 420

                     end do

  420                k2 = k

                        ildf(isqd) = k2 - 1

                     do k = 1, k2 - 1

                        chdf(isqd)(k:k) = dum1(k+isrm:k+isrm)

                     end do

                     isrm = isrm + k2
                     ilrm = ilrm - k2

                     goto 430

                  end if

               if( icc .eq. 1 ) then

                  write(iot,'("      reg = ",100a1)')
     &                    (chdf(1)(j:j),j=1,ildf(1))

               else if( icc .eq. 2 ) then

                  write(iot,'(" reginbox = ",100a1)')
     &                    (chdf(1)(j:j),j=1,ildf(1))

               end if

               if( isqd .gt. 1 ) then

                  do k = 2, isqd

                     write(iot,'(100a1)') (cblan(j:j),j=1,ild1),
     &                    (chdf(k)(j:j),j=1,ildf(k))

                  end do

               end if

*-----------------------------------------------------------------------
*        echo for combined regions
*-----------------------------------------------------------------------

*    &       junf .eq. 0 .and. icmb .eq. 0 .and. nvl .eq. 0 ) ) return

         if( iva .eq. 0 ) then

*    &                  "# combined, lattice or level structure ")')
            write(iot,'("   target =",i10,7x,
     &                  "# Card[10a] target nuclide ")') igcel

            write(iot,'("   non     reg      vol     #",
     &                  " reg definition")')

         else

            write(iot,'("   value  ",18x,
     &                  "# values for each region")')
            write(iot,'("   non     reg      val     #",
     &                  " reg definition")')

         end if

               call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)

*-----------------------------------------------------------------------
*     read from kr(mr)
*-----------------------------------------------------------------------

                  j = 0

         do 500 ir = 1, nr

               call echrg3(j,mr,kr,dum2,lng1,icmb,igm)

                  ild1 = 30
                  ild0 = max(10,ierl-ild1)

                  ilrm = lng1
                  isqd = 0
                  isrm = 0

  530             isqd = isqd + 1

                  if( ilrm .le. ild0 ) then

                        ildf(isqd) = ilrm

                     do k = 1, ilrm

                        chdf(isqd)(k:k) = dum2(k+isrm:k+isrm)

                     end do

                  else

                     do k = ild0, 1, -1

                        if( dum2(k+isrm:k+isrm) .eq. ' ' ) goto 520

                     end do

  520                k2 = k

                        ildf(isqd) = k2 - 1

                     do k = 1, k2 - 1

                        chdf(isqd)(k:k) = dum2(k+isrm:k+isrm)

                     end do

                     isrm = isrm + k2
                     ilrm = ilrm - k2

                     goto 530

                  end if

*-----------------------------
*     for input check

               iintrp = 0
               do irp = 1, itdci(24,m)
                if( lr(ir) .eq. itgcel(itdc2reg(itdc)+irp)
     &               .and. itarget .ne. 0 ) then !FURUTA20200522
                   iintrp = 1
                   goto 90
             end if

               end do
                     iintrp = 0

   90          if( iintrp .eq. 1 ) then
                  write(iot,'(i5,2x,i7,1p1e13.4," # ",100a1)')
     &                 ir, lr(ir), rtgvlm(itdc2reg(itdc)+irp), !FURUTA20200522
     &                 (chdf(1)(l:l),l=1,ildf(1))

               else
                  write(iot,'(i5,2x,i7,1p1e13.4," # ",100a1)')
     &                 ir, lr(ir), vl(ir), (chdf(1)(l:l),l=1,ildf(1))
               end if

               if( isqd .gt. 1 ) then

                  do k = 2, isqd

                     write(iot,'(100a1)') (cblan(l:l),l=1,ild1-2),
     &                    '#',' ',(chdf(k)(l:l),l=1,ildf(k))

                  end do

               end if

*-----------------------------------------------------------------------
*     for tg-list(material) in [T-Dchain]
*-----------------------------------------------------------------------
            if( iintrp .eq. 1 ) then

             i0=itdc2reg(itnm2tdc(m))  !FURUTA20200522
             j0=ireg2itm(i0+irp)       !FURUTA20200522

             write(iot,'("   tg-list =",i7)') itgnum(i0+irp)
             do jrp = 1, itgnum(i0+irp)
              write(iot,'(8x,a8,6x,1p1e13.7)')
     &             tglst(j0+jrp),tgnlt(j0+jrp) !FURUTA20200522
             end do

            else
*-----------------------------
*     get the cell No. -> material No.

                  jcel = 0

               if( lr(ir) .gt. 1000000 ) then
                  ic  = jnumc(chdf(1),3,100)
                  icl = inumc(chdf(1),ic,100,' ')
                  call onum(chdf(1),ic,icl,cvvv,ierr)
                  lr(ir) = nint(cvvv)
               end if

               do jcel = 1, igcel
                  if( lr(ir) .eq. idrg(jcel) ) goto 100
               end do

  100       continue

*-----------------------------
*     get the material No. -> material list

                  kmat = 0
               do kmat = 1, mxmat0
                  if( idmg(jcel) .eq. idmn(kmat) ) goto 110
                  if( idmn(idnm(idmg(jcel))) .eq. idmn(kmat) ) goto 110
               end do

  110       continue

               nel  = nint( das_kmatc(kmatc(kmat)+1) )
               denh =       das_kmatc(kmatc(kmat)+2)
               nelh = nel

               if( denh .ne. 0.0d0 ) then
                   nelh = nelh + 1
               end if

C MATSUDA 2017.10.04 (expand natural nucleus)
                  lmat = nelh
                  denstm = 0.0d0
               do kk = 1, nel
                  iz   = nint( das_kmatc(kmatc(kmat)+(kk-1)*4+34) )
                  ia   = nint( das_kmatc(kmatc(kmat)+(kk-1)*4+35) )
                  den  =       das_kmatc(kmatc(kmat)+(kk-1)*4+36)

                  if( ia .eq. 0 .and. natnn(iz) .ne. 0 ) then
                    lmat = lmat + natnn(iz) - 1
                  end if
               end do

               write(iot,'("   tg-list =",i7)') lmat

               if( denh .ne. 0.0d0 ) then
                   chau = ' H-1'
                   write(iot,'(8x,a4,8x,1p1e15.7)')
     &             chau(1:4), denh
               end if

               do kk = 1, nel
                  iz   = nint( das_kmatc(kmatc(kmat)+(kk-1)*4+34) )
                  ia   = nint( das_kmatc(kmatc(kmat)+(kk-1)*4+35) )
                  den  =       das_kmatc(kmatc(kmat)+(kk-1)*4+36)
                  chau = elmnt(iz)

                  if( ia .eq. 0 ) then
                    do kkk = 1, natnn(iz)

                      if( iz .eq. 1 .and. natnm(iz,kkk) .eq. 1 ) then
                          write(iot,'(8x,a2,"-1"8x,1p1e15.7)')
     &                    chau(1:2), denh + den * patnn(iz,kkk) / 100.0

                      else
                          denstm = den * patnn(iz,kkk) / 100.0

                        if( natnm(iz,kkk) .lt. 10 ) then
                          write(iot,'(8x,a2,"-",i1,8x,1p1e15.7)')
     &                    chau(1:2), natnm(iz,kkk), denstm
                        else if( natnm(iz,kkk) .lt. 100 ) then
                          write(iot,'(8x,a2,"-",i2,7x,1p1e15.7)')
     &                    chau(1:2), natnm(iz,kkk), denstm
                        else
                          write(iot,'(8x,a2,"-",i3,6x,1p1e15.7)')
     &                    chau(1:2), natnm(iz,kkk), denstm
                        end if

                      end if

                    end do

                  else if( ia .lt. 10 ) then
                     write(iot,'(8x,a2,"-",i1,8x,1p1e15.7)')
     &               chau(1:2), ia, den
                  else if( ia .lt. 100 ) then
                     write(iot,'(8x,a2,"-",i2,7x,1p1e15.7)')
     &               chau(1:2), ia, den
                  else
                     write(iot,'(8x,a2,"-",i3,6x,1p1e15.7)')
     &               chau(1:2), ia, den
                  end if
               end do

            end if
*-----------------------------------------------------------------------

  500    continue

*-----------------------------------------------------------------------

      return
      end subroutine

************************************************************************
*                                                                      *
      subroutine tdctoinp(m,fname)
*     Conversion program from PHITS [T-Dchain] to DCHAIN input         *
*     designed by N.Matsuda on 2012/11/08                              *
*          modified by N.Matsuda on 2015/01/16                         *
*          modified by T.Miura on 2018/07/31                           *
*     last modified by T.Furuta on 2020/06/01                          *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'err.inc'

*-----------------------------------------------------------------------
      common /emode/  emodem, ge1, ge2, iemode
      common /talldc2/ itarget
      common /ivoll/ ivoll
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /rdfisparm/ irdfisyd
      common /clionprd/  lionprd ! T.Sato 2025/02/09, to judge itriton
*-----------------------------------------------------------------------
***   for sub_readl
      integer jsn, jpn, i1, i2, i3, i4, iskip, ierr
      integer iot, iodr, iodw, iryld, iwyld, iotrk, ionfl, ioln

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension illd(0:9), ilfd(0:9)
      dimension illy(0:9), ilfy(0:9)
      dimension illt(0:9), ilft(0:9)

      character chin*200, chlw*200, chcm*200

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------
***   for section parameters
      character chlc*200
      integer icl, i5, icsu, il, ipm, ic

      parameter ( nkwargb = 111 )
CCSE added r-z/xyz mesh, schan parameters, 69 -> 73 (2018.07.31) >>>>>
      dimension lschn(nkwargb), ischn(nkwargb)
      character schan(nkwargb)*8

      common /ndemax/ dnmax(20) ! T.Sato 2025/02/08 for deterining itriton

      data icsu / nkwargb /

      data ( schan(i), i = 1, nkwargb ) /
     &    'mesh    ','special ','mother  ','nucleus ','axis    ',
     &    'file    ','title   ','angel   ','unit    ','info    ',
     &    '2d-type ','factor  ','material','x-txt   ','y-txt   ',
     &    'z-txt   ','gshow   ','rshow   ','ndata   ','iechrl  ',
     &    'volmat  ','epsout  ','ctmin(1)','ctmax(1)','ctmin(2)',
     &    'ctmax(2)','ctmin(3)','ctmax(3)','resol   ','width   ',
     &    'trcl    ','*trcl   ','part    ','gslat   ','output  ',
     &    'resfile ','imode   ','jmode   ','itstep  ','itout   ',
     &    'idivs   ','iregon  ','inmtcf  ','ichain  ','itdecs  ',
     &    'itdecn  ','isomtr  ','ifisyd  ','ifisye  ','iyild   ',
     &    'iggrp   ','ibetap  ','acmin   ','istabl  ','igsdef  ',
     &    'iprtb1  ','iprtb2  ','rprtb2  ','iprtb3  ','igsorg  ',
     &    'amp     ','ebeam   ','prodnp  ','timeevo ','outtime ',
     &    'ac-list ','target  ','dversion','sangel  ','imat    ',
     &    'mtscore ','mt-list ','imesh   ','iertdcho','inxslib ',
     &    'idcylib ','infylib ','hnxslib ','hdcylib ','hnfylib ',
     &    'iwrtchn ','chrlvth ','iwrchdt ','iwrchss ',
     &    'ixsrall ','idosecf ','ipltmode','ipltaxis',
     &    'foamout ','foamvals','iredufmt','irdonce ',
     &    'imtcard ','imtcnum ','imtcmeta','thmatnd ','idosunit',
     &    'iangbpwr','iphtout ','isfylib ','radlev  ','ergthfa ',
     &    'ergfafu ','infyerg ','hsfylib ','ilchain ','mxnuclei',
     &    'itriton ','iaonucl ','aonucl  ','aoreg   '/

      data ( lschn(i), i = 1, nkwargb ) /
     &     4,         7,         6,         7,         4,
     &     4,         5,         5,         4,         4,
     &     7,         6,         8,         5,         5,
     &     5,         5,         5,         5,         6,
     &     6,         6,         8,         8,         8,
     &     8,         8,         8,         5,         5,
     &     4,         5,         4,         5,         6,
     &     7,         5,         5,         6,         5,
     &     5,         6,         6,         6,         6,
     &     6,         6,         6,         6,         5,
     &     5,         6,         5,         6,         6,
     &     6,         6,         6,         6,         6,
     &     3,         5,         6,         7,         7,
     &     7,         6,         8,         6,         4,
     &     7,         7,         5,         8,         7,
     &     7,         7,         7,         7,         7,
     &     7,         7,         7,         7,
     &     7,         7,         8,         8,
     &     7,         8,         8,         7,
     &     7,         7,         8,         7,         8,
     &     8,         7,         7,         6,         7,
     &     7,         7,         7,         7,         8,
     &     7,         7,         6,         5/

*-----------------------------------------------------------------------
***   for individual parameters, ipm = 1, 6, 7
*     file name
      integer infil, icf
      character ifile(6)*100
      dimension lfile(6)
CCSE added for mntc_yield & n.flux_ file (2017.12.31) >>>>>
      character hhnmtcf*100
      character hnfluxs*100
CCSE added for mntc_yield & n.flux_ file (2017.12.31) <<<<<
c H. Ratliff 2019.07.30 added yield uncertainty file henmtcf
      character henmtcf*100
      character hnxslib*80
      character hdcylib*80
      character hnfylib*80
      character hsfylib*80
*     title
      integer ict, titll
      character title*80
      character aonucl*200, aoreg*200 ! T.Sato 2025/02/12

      real*8 radlev,ergthfa,ergfafu,chrlvth

*-----------------------------------------------------------------------
***   for dchain parameters, ipm = 37 - 67
      dimension temsh(4,maxitm), itemsh(2,maxitm)

      integer itetp, iottp, iancza, itgnum2
      dimension ancza(maxitm)
      character hancza(maxitm)*8
      character hl*8
      integer izam, iz, ia, im, j1, j2
      character chaia1*1, chaia2*2, chaia3*3, chaim*1

      integer klst, ic2

*-----------------------------------------------------------------------
      character element(104)*3,cnuc*3,elmnt(104)*3

      data element /
     &    'h  ','he ','li ','be ','b  ','c  ','n  ','o  ','f  ','ne ',
     &    'na ','mg ','al ','si ','p  ','s  ','cl ','ar ','k  ','ca ',
     &    'sc ','ti ','v  ','cr ','mn ','fe ','co ','ni ','cu ','zn ',
     &    'ga ','ge ','as ','se ','br ','kr ','rb ','sr ','y  ','zr ',
     &    'nb ','mo ','tc ','ru ','rh ','pd ','ag ','cd ','in ','sn ',
     &    'sb ','te ','i  ','xe ','cs ','ba ','la ','ce ','pr ','nd ',
     &    'pm ','sm ','eu ','gd ','tb ','dy ','ho ','er ','tm ','yb ',
     &    'lu ','hf ','ta ','w  ','re ','os ','ir ','pt ','au ','hg ',
     &    'tl ','pb ','bi ','po ','at ','rn ','fr ','ra ','ac ','th ',
     &    'pa ','u  ','np ','pu ','am ','cm ','bk ','cf ','es ','fm ',
     &    'md ','no ','lr ','ku '/

      data elmnt /
     &    ' H ','He ','Li ','Be ',' B ',' C ',' N ',' O ',' F ','Ne ',
     &    'Na ','Mg ','Al ','Si ',' P ',' S ','Cl ','Ar ',' K ','Ca ',
     &    'Sc ','Ti ',' V ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ',' Y ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ',' I ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ',' W ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ',' U ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

*-----------------------------------------------------------------------
***   for summary
      character filnm*100
      integer ifilnm, itfp, kk
*     Card[9]
      logical lex
      integer kase
      character filyld*10
*     Card[10]
!      dimension eandf(1969,4) ! T.Sato 2024/11/17, not used anymore
      character filnfl*10
*     DCHAIN input file
      character chau*1

*-----------------------------------------------------------------------
***   for error
      character tnamed*10
      data      tnamed /'[t-dchain]'/

*-----------------------------------------------------------------------
      common /paran/  icfn(100), ilfn(100), chfn(100)
      character       chfn*200
      character fname*100

CCSE added r-z/xyz mesh (2018.07.31) >>>>>
      character chlst*200
      integer(8),allocatable :: ildata(:,:)
      real(8),allocatable :: rldata(:,:)
CCSE added r-z/xyz mesh (2018.07.31) <<<<<
      integer foamout, foamvals
*-----------------------------------------------------------------------
      iot   =   6
      iodr  = 802
      iodw  = 803
      iryld = 804
      iwyld = 805
      iotrk = 806
      ionfl = 807
      ioln  = 808

      klst  =   0
      lex   = .false.
      cvvv  =   0.0

*-----------------------------------------------------------------------
*     default values of [T-Dchain] parameters
*-----------------------------------------------------------------------

*     Card[1] title
            titll = 0
*     Card[2] primary control data
            imode = 2
            jmode = 2
*     Card[3] first calculation condition data
           idcstp = 1
           idcout = 1
            idivs = 50
           iregon = 1
           inmtcf = 1
           ichain = 100
           itdecs = 1
           itdecn = 1
           isomtr = 2
           ifisyd = 0
           ifisye = 0
*     Card[4] output condition data
            iyild = 2
            iggrp = 3
           ibetap = 1
            acmin = 1.d-20
           istabl = 0
           igsdef = 0   ! T.Sato 2025/02/10 change 1 to 0
           iprtb1 = 1
           iprtb2 = 1
           rprtb2 = 10
           iprtb3 = 0
           igsorg = 0   ! T.Sato 2025/02/10 change 1 to 0
*     Card[5] spallation data
c iwamoto 2012/12/10 amp=source/sec
               amp= 1.0
              ampp= 1.0
            ebeam = 3.0
           prodnp = 1.0
*     Card[6] irra. / cooling info.
            itetp = 0
*     Card[7] output time
            iottp = 0
*     Card[8] pick up nuclide
           iancza = 0
*     Card[10] neutron flux
             flux = 1.0
           ifilnm = 0
            itarget = 0
         iversion = 1

CCSE added r-z/xyz mesh (2018.07.31) >>>>>
             imat = 0
          mtscore = 0
            imesh = 0

           mtlflg = 0
CCSE added r-z/xyz mesh (2018.07.31) <<<<<
c H.Ratliff 2019.07.30 added DCHAIN uncertainty propagation toggle
*     Card[9] toggle of error propagation in DCHAIN
         iertdcho = 1
*     Card[3b] data library selection
         inxslib  = 100
         idcylib  = 5
         infylib  = 17 ! H.Ratliff 2021.04.07 re-enable fission
         write(hnxslib,"('custom-lib',70x)")
         write(hdcylib,"('custom-lib',70x)")
         write(hnfylib,"('custom-lib',70x)")
         write(hsfylib,"('custom-lib',70x)")
*     Card[4] output options for chain data file
         iwrtchn  = 1
         chrlvth  = -1.d0
         iwrchdt  = 0
         iwrchss  = 0
*     Card[3] calculation option for loading all cross sections into memory
          ixsrall = 1
*     Card[4] output options for dose rate coefficient
          idosecf = 1
*     Card[9d] output plot options for mesh geometries
        ipltmode = 0
        ipltaxis = 1
*     Card[9d] OpenFoam output options for tet mesh geometries
        foamout  = 0
        foamvals = 0
        iredufmt0 = 1
        irdonce = 1
        imtcard = 1
        imtcnum = 0
        imtcmeta = 2
        thmatnd = 1.d-6
        idosunit = 333
        iangbpwr = 1
        iphtout = 0  ! T.Sato 2025/01/31 change from 1 (gamma spectrum) to 0 (no output)
        itriton = 0  ! T.Sato 2025/02/08 tritium-production cross section (0: consider, 1:ignore when HE library avaiable, 2: always ignore)
        if(lionprd.eq.1) itriton = 1 ! tritium production should be already included in yield data
        iaonucl = 0   ! T.Sato 2025/02/12 for angelout_nuclides
        laonucl = 0 ! T.Sato 2025/02/12 for character length of angelout_nuclides data
        laoreg = 0  ! T.Sato 2025/02/12 for character length of angelout_region data

! H.Ratliff 2021.04.07 re-enable fission
        isfylib = 1
        radlev  = 1.0e-3 ! T.Sato 2025/02/08
        ergthfa = 0.1
        ergfafu = 10.0
        infyerg = 4
        ilchain = 100

*-----------------------------------------------------------------------
      open(ioln,file='dch_link.dat',form='formatted',status='unknown')
        write(ioln,'(80a)') chfn(21)(1:80)
      close(ioln)
*-----------------------------------------------------------------------
*     READ START <- [T-Dchain] output
*-----------------------------------------------------------------------
      open(iodr, file = fname, form = 'formatted', status = 'unknown')
      rewind iodr

*-----------------------------------------------------------------------
*     read input file
*     initial value for include file
*-----------------------------------------------------------------------

         do i = 0, 9

            illd(i) = 0
            ilfd(i) = 2147483647 !FURUTA20200522
            illy(i) = 0
            ilfy(i) = 2147483647 !FURUTA20200522
            illt(i) = 0
            ilft(i) = 2147483647 !FURUTA20200522

         end do
c 2012/12/17 iwamoto, modified by T.Sato from 1 to maxmreg
         do i = 1, maxreg !FURUTA20200522
            rtgvol(i) = 1.0
         end do

ccse 2023/03/29, zero set itgvoll
         itgvoll = 0

         infil = 0
           jsn = 0
           jpn = 0
         iskip = 0
          ierr = 0

          itfp = 0

*-----------------------------------------------------------------------
*     search [T-Dchain] section label
*-----------------------------------------------------------------------

  130    continue

           call readl(jsn,iodr,dsin,idsi,illd,ilfd,'#!$',
     &                jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 921

               if( iskip .ne. 0 ) goto 130

               if( chcm(i1:i1+9) .eq. '[t-dchain]' ) then
                  goto 140
               else
                  goto 130
               end if

*-----------------------------------------------------------------------
*     read one line from iodr
*-----------------------------------------------------------------------

  140    continue

           call readl(jsn,iodr,dsin,idsi,illd,ilfd,'#!$',
     &                jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 800

               if( i1.ne.0 ) then
                 if( chin(i1:i1+3) .eq. '#---' ) goto 160
               end if
               if( iskip .ne. 0 ) goto 140

  150    continue
               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 800

*-----------------------------------------------------------------------
*        end of the section
*-----------------------------------------------------------------------

           if( chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 800

           end if

*-----------------------------------------------------------------------
*        end of echo read
*-----------------------------------------------------------------------

  160    continue
           if( chin(i1:i1+3) .eq. '#---' ) then

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
*.......................................................................
*        sumtally parameters
         if ( chlw(i1:i1+13)=='sumtally start' ) then
            do
               call readl(jsn,iodr,dsin,idsi,illd,ilfd,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
               if ( ierr /= 0 ) return
               if ( jpn  == 3 ) goto 800
               if ( chlw(i1:i1+11)=='sumtally end' ) go to 140
            end do
         end if
*.......................................................................

               goto 987

*-----------------------------------------------------------------------
*        read value of parameters
*-----------------------------------------------------------------------

  100    continue

               ipm = i
C        check for dual definition
*    &       ischn( ipm ) .gt. 1 ) goto 989

           ic = inumc(chlw,il+1,i3,'=') + 1
           ic = jnumc(chlw,ic,i3)

           if( ic .gt. i3 ) goto 997

           icl = inumc(chlw,ic,i3,';') - 1

*-----------------------------------------------------------------------
*        mesh
*-----------------------------------------------------------------------

         if( ipm .eq. 1 ) then
CCSE added r-z/xyz mesh (2018.07.31) >>>>>
            if( chlw(ic:ic+2) .eq. 'reg' ) then
               meshtype = 1
            else if( chlw(ic:ic+2) .eq. 'r-z' ) then
               meshtype = 2
            else if( chlw(ic:ic+2) .eq. 'xyz' ) then
               meshtype = 3
            else if( chlw(ic:ic+2) .eq. 'tet' ) then
               meshtype = 4 !FURUTA20200505
            else
               goto 920
            end if
CCSE added r-z/xyz mesh (2018.07.31) <<<<<

  210       call readl(jsn,iodr,dsin,idsi,illd,ilfd,'#!$',
     &           jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 800

               if( chlw(i1:i1+5) .eq. 'target' ) goto 150
CCSE added r-z/xyz mesh (2018.07.31) >>>>>
               if( chlw(i1:i1+3) .eq. 'imat' ) goto 150
CCSE added r-z/xyz mesh (2018.07.31) <<<<<
               if( iskip .ne. 0 ) goto 210

               goto 210

*-----------------------------------------------------------------------
*        file name
*-----------------------------------------------------------------------

         else if( ipm .eq. 6 ) then
*           file name is modified for DCHAIN mode.

  700          infil = infil + 1

               if( infil .gt. 2 ) goto 990

               icf = min( inumc(chlw,ic,icl,' ') - 1, icl )
*    &                     ic, icl, icf, i3
               lfile(infil) = icf - ic + 1
               ifile(infil)(1:icf-ic+1) = chin(ic:icf) !FURUTA20130821

               ic = jnumc(chlw,icf+2,icl)

               if( ic .le. icl ) goto 700

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        title
*-----------------------------------------------------------------------

         else if( ipm .eq. 7 ) then

               ict = min( ic + 79, i2 )

               title = chin(ic:ict)

               titll = ict - ic + 1

*-----------------------------------------------------------------------
*        The other DCHAIN option (control) ipm = 37 to 67
*-----------------------------------------------------------------------
*        imode     in Card[2]-1
*-----------------------------------------------------------------------

         else if( ipm .eq. 37 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 922

               imode = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        jmode     in Card[2]-2
*-----------------------------------------------------------------------

         else if( ipm .eq. 38 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 923

               jmode = nint( cvvv )
               IF(iemode .ne. 0) jmode = 0  ! 2013/7/31 Ogawa

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        itstep    in Card[3]-1
*-----------------------------------------------------------------------

         else if( ipm .eq. 39 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 924

               idcstp = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        itout     in Card[3]-2
*-----------------------------------------------------------------------

         else if( ipm .eq. 40 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 925

               idcout = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        idivs     in Card[3]-3
*-----------------------------------------------------------------------

         else if( ipm .eq. 41 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 926

               idivs = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        iregon    in Card[3]-4
*-----------------------------------------------------------------------

         else if( ipm .eq. 42 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 927

               iregon = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        inmtcf    in Card[3]-5
*-----------------------------------------------------------------------

         else if( ipm .eq. 43 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 928

               inmtcf = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        ichain    in Card[3]-6
*-----------------------------------------------------------------------

         else if( ipm .eq. 44 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 929

               ichain = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        itdecs    in Card[3]-7
*-----------------------------------------------------------------------

         else if( ipm .eq. 45 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 930

               itdecs = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        itdecn    in Card[3]-8
*-----------------------------------------------------------------------

         else if( ipm .eq. 46 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 931

               itdecn = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        isomtr    in Card[3]-9
*-----------------------------------------------------------------------

         else if( ipm .eq. 47 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 932

               isomtr = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        ifisyd    in Card[3]-10
*-----------------------------------------------------------------------

         else if( ipm .eq. 48 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 933

               ifisyd = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        ifisye    in Card[3]-11
*-----------------------------------------------------------------------

         else if( ipm .eq. 49 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 934

               ifisye = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        iyild     in Card[4]-1
*-----------------------------------------------------------------------

         else if( ipm .eq. 50 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 935

               iyild = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        iggrp     in Card[4]-2
*-----------------------------------------------------------------------

         else if( ipm .eq. 51 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 936

               iggrp = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        ibetap    in Card[4]-3
*-----------------------------------------------------------------------

         else if( ipm .eq. 52 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 937

               ibetap = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        acmin     in Card[4]-4
*-----------------------------------------------------------------------

         else if( ipm .eq. 53 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 938

               acmin = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        istabl    in Card[4]-5
*-----------------------------------------------------------------------

         else if( ipm .eq. 54 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 939

               istabl = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        igsdef    in Card[4]-6
*-----------------------------------------------------------------------

         else if( ipm .eq. 55 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 940

               igsdef = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        iprtb1    in Card[4]-7
*-----------------------------------------------------------------------

         else if( ipm .eq. 56 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 941

               iprtb1 = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        iprtb2    in Card[4]-8
*-----------------------------------------------------------------------

         else if( ipm .eq. 57 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 942

               iprtb2 = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        rprtb2    in Card[4]-9
*-----------------------------------------------------------------------

         else if( ipm .eq. 58 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 943

               rprtb2 = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        iprtb3    in Card[4]-10
*-----------------------------------------------------------------------

         else if( ipm .eq. 59 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 944

               iprtb3 = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        igsorg    in Card[4]-11
*-----------------------------------------------------------------------

         else if( ipm .eq. 60 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 945

               igsorg = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        amp       in Card[5]-1
*-----------------------------------------------------------------------

         else if( ipm .eq. 61 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 946

               amp = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        ebeam     in Card[5]-2
*-----------------------------------------------------------------------

         else if( ipm .eq. 62 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 947

               ebeam = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        prodnp    in Card[5]-3
*-----------------------------------------------------------------------

         else if( ipm .eq. 63 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 948

               prodnp = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        version    DCHAIN-SP2001 or DCHAIN-SP2015
*-----------------------------------------------------------------------

         else if( ipm .eq. 68 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 957

               iversion = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        tbin and beampw    in Card[6]  irra. / cooling info.
*-----------------------------------------------------------------------

         else if( ipm .eq. 64 ) then


               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 952

               itetp = nint( cvvv )

               call gettev(jsn,iodr,dsin,idsi,illd,ilfd,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     1,itetp,temsh,itemsh)

               if( ierr .ne. 0 ) goto 952
               if( itetp .ne. idcstp ) idcstp = itetp

               if( ierr .ne. 0 ) return
               goto 140

*-----------------------------------------------------------------------
*        tmin    in Card[7]  output time
*-----------------------------------------------------------------------

         else if( ipm .eq. 65 ) then


               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 953

               iottp = nint( cvvv )

               call gettev(jsn,iodr,dsin,idsi,illd,ilfd,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     2,iottp,temsh,itemsh)

               if( ierr .ne. 0 ) goto 953
               if( iottp .ne. idcout ) idcout = iottp

               if( ierr .ne. 0 ) return
               goto 140

*-----------------------------------------------------------------------
*        ancza or hancza    in Card[8]  pick up nuclide
*-----------------------------------------------------------------------

         else if( ipm .eq. 66 ) then


               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 949

               iancza = nint( cvvv )

               if( iprtb2 .eq. 3 ) rprtb2 = dble(iancza)
                   iprtb3 = iancza

*---
*---
! T.Sato 2015/8/25, delete tgnzas & ctgnnds from parameter list
                     call getndc(jsn,iodr,dsin,idsi,illd,ilfd,
     &                           jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                           1,iancza,0,hancza,ancza)

                       if( ierr .ne. 0 ) return
                       if( jpn  .eq. 3 ) goto 954

               goto 140

*-----------------------------------------------------------------------
*        hhnmtcf    in Card[9]  [T-yield] file name
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*        hregcmm    in Card[10] and [10a]  special calc. info.
*-----------------------------------------------------------------------

         else if( ipm .eq. 67 ) then

CCSE added r-z/xyz mesh (2018.07.31) >>>>>
            if( meshtype .ne. 1 ) goto 955
CCSE added r-z/xyz mesh (2018.07.31) <<<<<

*         read data at the last minute.

            klst = 0

  613       call readl(jsn,iodr,dsin,idsi,illd,ilfd,'#!$',
     &           jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 800

               if( chlw(i1:i1+2) .eq. 'non' .or.
     &             chlw(i1:i1+2) .eq. 'reg' .or.
     &             chlw(i1:i1+2) .eq. 'vol' ) goto 613
               if( iskip .ne. 0 ) goto 613
               ic = i1

*           for number of target
  614          call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 955
               klst = klst + 1

               if(klst.gt.maxreg)call reallocate_talldc4(klst) !FURUTA20200522
               itdcreg=itdcreg+1             !FURUTA20200522
               if(itdcreg.gt.maxtdcreg)call reallocate_talldc2(itdcreg) !FURUTA20200522
               ireg2itm(itdcreg)=itdcregitm  !FURUTA20200522
               jreg2itm(klst)=iregitm        !FURUTA20200522



               ic = jnumc(chlw,ic2,i3)
               if( ic .gt. i3 ) goto 955

*           for cell No.
               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 955

                  itgcell(klst) = nint(cvvv)

               ic = jnumc(chlw,ic2,i3)
               if( ic .gt. i3 ) goto 955

*           for volume
               call onum(chlw,ic,i3,cvvv,ierr)
               if( ierr .ne. 0 ) goto 955

                  rtgvol(klst) = cvvv

                  if( rtgvol(klst) .gt. 0.0d0 ) itgvoll = 1

  615       call readl(jsn,iodr,dsin,idsi,illd,ilfd,'#!$',
     &           jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 955

               if( iskip .ne. 0 ) goto 615

                  ic = i1

               if( chlw(ic:ic+6) .eq. 'tg-list' ) then
                  ic = inumc(chlw,ic,i3,'=') + 1
                  ic = jnumc(chlw,ic,i3)

                  if( ic .gt. i3 ) goto 955

               else
                  goto 955
               end if

*           for number of material
                     call onum(chlw,ic,i3,cvvv,ierr)
                        if( ierr .ne. 0 ) goto 955

                        itglist(klst) = nint( cvvv )
                        itgnum2 = itglist(klst)

*---
*---
*           for material and ratio
! T.Sato 2015/8/25, delete tgnzas & ctgnnds from parameter list
                     call getndc(jsn,iodr,dsin,idsi,illd,ilfd,
     &                           jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                           2,itgnum2,klst,hancza,ancza)

                       if( ierr .ne. 0 ) return
                       if( jpn  .eq. 3 ) goto 955

  622    call readl(jsn,iodr,dsin,idsi,illd,ilfd,'#!$',
     &        jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

            if( ierr .ne. 0 ) return
            if( jpn  .eq. 3 ) goto 800

            if( chin(1:4) .eq. '#---' ) goto 150
            if( iskip .ne. 0 ) goto 622

               ic = i1
            if( chlw(ic:ic) .ne. '[' .and. chlw(ic:ic) .lt. 'a' ) then
               goto 614
            else
               goto 150
            end if

CCSE added r-z/xyz mesh (2018.07.31) >>>>>
*-----------------------------------------------------------------------
*        imat: number of mt-list for r-z and xyz mesh
*-----------------------------------------------------------------------

         elseif( ipm .eq. 70 ) then

            if( meshtype .lt. 2 .or. meshtype .gt. 4 ) goto 955 !FURUTA20200505

            call onum(chlw,ic,icl,cvvv,ierr)
            if( ierr .ne. 0 ) goto 955

            imat = nint( cvvv )

            icl = jnumc(chlw,icl+2,i3)
            if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        mtscore: material number for DCHAIN calculation for r-z and xyz mesh
*-----------------------------------------------------------------------

         elseif( ipm .eq. 71 ) then

            if( meshtype .lt. 2 .or. meshtype .gt. 4 ) goto 955 !FURUTA20200505

            call onum(chlw,ic,icl,cvvv,ierr)
            if( ierr .ne. 0 ) goto 955

            mtscore = nint( cvvv )

            icl = jnumc(chlw,icl+2,i3)
            if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        mt-list: list of elements for r-z and xyz mesh
*-----------------------------------------------------------------------

         elseif( ipm .eq. 72 ) then

            if( meshtype .lt. 2 .or. meshtype .gt. 4 ) goto 955 !FURUTA20200505
            if( imat .le. 0 ) goto 955
            if( mtlflg .ne. 0 ) goto 955

            mtlflg = 1

            klst = 0

            do i = 1, imat

               klst = klst + 1

               if(klst.gt.maxreg)call reallocate_talldc4(klst) !FURUTA20200522
               itdcreg=itdcreg+1             !FURUTA20200522
               if(itdcreg.gt.maxtdcreg)call reallocate_talldc2(itdcreg) !FURUTA20200522
               ireg2itm(itdcreg)=itdcregitm  !FURUTA20200522
               jreg2itm(klst)=iregitm        !FURUTA20200522

               if( klst .gt. imat ) goto 955

               call onum(chlw,ic,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 955

               itglist(klst) = nint( cvvv )
               itgnum2 = itglist(klst)

               call getndc(jsn,iodr,dsin,idsi,illd,ilfd,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     2,itgnum2,klst,hancza,ancza)
               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 955

               if( i .lt. imat ) then
  623             call readl(jsn,iodr,dsin,idsi,illd,ilfd,'#!$',
     &                       jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
                  if( ierr .ne. 0 ) return
                  if( jpn  .eq. 3 ) goto 800
                  if( iskip .ne. 0 ) goto 623

                  ic = inumc(chlw,il+1,i3,'=') + 1
               endif

            enddo

*-----------------------------------------------------------------------
*        imesh: number of mesh grid for r-z and xyz mesh
*-----------------------------------------------------------------------

         elseif( ipm .eq. 73 ) then

            if( meshtype .eq. 2 .or. meshtype .eq. 3)then !FURUTA20200505

               if( imat .le. 0 ) goto 955

               call onum(chlw,ic,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 955

               imesh = nint( cvvv )

! read non ix iy ...
 624           call readl(jsn,iodr,dsin,idsi,illd,ilfd,'#!$',
     &              jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 800
               if( iskip .ne. 0 ) goto 624
               ic = i1

               if( chlw(i1:i1+2) .ne. "non" ) goto 955

               do i = i1, i3
                  if( chlw(i:i+1) .eq. "iz" ) goto 625
               enddo
               goto 955

 625           continue
               chlst = chlw(1:i+1)//"        fluxs"//chlw(i+2:i3)
               ichlst = i3 + 13

               if( meshtype .eq. 2 ) then
                  iclmesh = 3
               elseif( meshtype .eq. 3 ) then
                  iclmesh = 4
               endif

               iclvol = imat + 1

! allocate
               allocate( ildata(imesh,iclmesh) )
               allocate( rldata(imesh,iclvol)  )

               do i = 1, imesh

 626              call readl(jsn,iodr,dsin,idsi,illd,ilfd,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
                  if( ierr .ne. 0 ) return
                  if( jpn  .eq. 3 ) goto 800
                  if( iskip .ne. 0 ) goto 626
                  ic = i1

                  do j = 1, iclmesh
                     call snum(chlw,ic,i3,ic2,cvvv,ierr)
                     if( ierr .ne. 0 ) goto 955

                     ildata(i,j) = nint(cvvv)

                     ic = jnumc(chlw,ic2,i3)
                     if( ic .gt. i3 ) goto 955
                  enddo

                  do j = 1, iclvol
                     call snum(chlw,ic,i3,ic2,cvvv,ierr)
                     if( ierr .ne. 0 ) goto 955

                     rldata(i,j) = cvvv

                     if( j .lt. iclvol ) then
                        ic = jnumc(chlw,ic2,i3)
                        if( ic .gt. i3 ) goto 955
                     endif
                  enddo

               enddo

            elseif( meshtype.eq.4 )then
               if( imat .le. 0 ) goto 955

               call onum(chlw,ic,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 955

               imesh = nint( cvvv )

! read non tetra ...
 627           call readl(jsn,iodr,dsin,idsi,illd,ilfd,'#!$',
     &              jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 800
               if( iskip .ne. 0 ) goto 627
               ic = i1

               if( chlw(i1:i1+2) .ne. "non" ) goto 955

               do i = i1, i3
                  if( chlw(i:i+1) .eq. "mt" ) goto 628
               enddo
               goto 955

 628           continue
               chlst = chlw(1:i+1)//"        fluxs"//chlw(i+2:i3)
               ichlst = i3 + 13
               iclmesh=3
               iclvol=4 !FURUTA20200529
! allocate
               allocate( ildata(imesh,iclmesh) )
               allocate( rldata(imesh,iclvol)  )

               do i = 1, imesh

 629              call readl(jsn,iodr,dsin,idsi,illd,ilfd,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
                  if( ierr .ne. 0 ) return
                  if( jpn  .eq. 3 ) goto 800
                  if( iskip .ne. 0 ) goto 629
                  ic = i1

                  do j = 1, iclmesh
                     call snum(chlw,ic,i3,ic2,cvvv,ierr)
                     if( ierr .ne. 0 ) goto 955

                     ildata(i,j) = nint(cvvv)

                     ic = jnumc(chlw,ic2,i3)
                     if( ic .gt. i3 ) goto 955
                  enddo

                  do j = 1, iclvol
                     call snum(chlw,ic,i3,ic2,cvvv,ierr)
                     if( ierr .ne. 0 ) goto 955

                     rldata(i,j) = cvvv

                     if( j .lt. iclvol ) then
                        ic = jnumc(chlw,ic2,i3)
                        if( ic .gt. i3 ) goto 955
                     endif
                  enddo

               enddo

            else
               goto 955
            endif
c H.Ratliff START 2019.07.30, add uncertainty propagation, library selection,
c                             and decay chain info file output options

*-----------------------------------------------------------------------
*        iertdcho    in Card[9]-5
*-----------------------------------------------------------------------

         else if( ipm .eq. 74 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 821

               iertdcho = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200


*-----------------------------------------------------------------------
*        inxslib    in Card[3b]-1
*-----------------------------------------------------------------------

         else if( ipm .eq. 75 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 822

               inxslib = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200


*-----------------------------------------------------------------------
*        idcylib    in Card[3b]-2
*-----------------------------------------------------------------------

         else if( ipm .eq. 76 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 823

               idcylib = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200


*-----------------------------------------------------------------------
*        infylib    in Card[3b]-3
* NOTE: This option currently does nothing; it has just been implemented
*       in case this functionality is added in the future
*-----------------------------------------------------------------------

         else if( ipm .eq. 77 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 824

               infylib = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        hnxslib    in Card[3b]-4
*-----------------------------------------------------------------------

         else if( ipm .eq. 78 ) then

               ict = min( ic + 79, i2 )

               hnxslib = chin(ic:ict)

               hxsll = ict - ic + 1


*-----------------------------------------------------------------------
*        hdcylib    in Card[3b]-5
*-----------------------------------------------------------------------

         else if( ipm .eq. 79 ) then

               ict = min( ic + 79, i2 )

               hdcylib = chin(ic:ict)

               hdcll = ict - ic + 1


*-----------------------------------------------------------------------
*        hnfylib    in Card[3b]-6
*-----------------------------------------------------------------------

         else if( ipm .eq. 80 ) then

               ict = min( ic + 79, i2 )

               hnfylib = chin(ic:ict)

               hfsll = ict - ic + 1



*-----------------------------------------------------------------------
*        iwrtchn    in Card[4]-12
*-----------------------------------------------------------------------

         else if( ipm .eq. 81 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 828

               iwrtchn = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200


*-----------------------------------------------------------------------
*        chrlvth    in Card[4]-13
*-----------------------------------------------------------------------

         else if( ipm .eq. 82 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 829

               chrlvth = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200


*-----------------------------------------------------------------------
*        iwrchdt    in Card[4]-14
*-----------------------------------------------------------------------

         else if( ipm .eq. 83 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 2830

               iwrchdt = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200


*-----------------------------------------------------------------------
*        iwrchss    in Card[4]-15
*-----------------------------------------------------------------------

         else if( ipm .eq. 84 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 831

               iwrchss = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200


c H.Ratliff START 2020.04.17

*-----------------------------------------------------------------------
*        ixsrall    in Card[3]-7
*-----------------------------------------------------------------------

         else if( ipm .eq. 85 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 832

               ixsrall = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        idosecf    in Card[4]-16
*-----------------------------------------------------------------------

         else if( ipm .eq. 86 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 833

               idosecf = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        ipltmode    in Card[9d]-1
*-----------------------------------------------------------------------

         else if( ipm .eq. 87 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 834

               ipltmode = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        ipltaxis    in Card[9d]-2
*-----------------------------------------------------------------------

         else if( ipm .eq. 88 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 835

               ipltaxis = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

c H.Ratliff START 2020.05.29

*-----------------------------------------------------------------------
*        foamout    in Card[9d]-1
*-----------------------------------------------------------------------

         else if( ipm .eq. 89 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 836

               foamout = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        foamvals    in Card[9d]-2
*-----------------------------------------------------------------------

         else if( ipm .eq. 90 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 837

               foamvals = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

c H.Ratliff START 2020.06.05
*-----------------------------------------------------------------------
*        iredufmt    in Card[8]-7
*-----------------------------------------------------------------------

         else if( ipm .eq. 91 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 838

               iredufmt0 = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

c H.Ratliff START 2020.06.17
*-----------------------------------------------------------------------
*        irdonce    in Card[8]-7
*-----------------------------------------------------------------------

         else if( ipm .eq. 92 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 839

               irdonce = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
c H.Ratliff START 2020.11.16
*-----------------------------------------------------------------------
*        imtcard    in Card[4]-17
*-----------------------------------------------------------------------

         else if( ipm .eq. 93 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 840

               imtcard = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        imtcnum    in Card[4]-18
*-----------------------------------------------------------------------

         else if( ipm .eq. 94 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 841

               imtcnum = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        imtcmeta    in Card[4]-19
*-----------------------------------------------------------------------

         else if( ipm .eq. 95 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 842

               imtcmeta = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        thmatnd    in Card[4]-20
*-----------------------------------------------------------------------

         else if( ipm .eq. 96 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 843

               thmatnd = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

c H.Ratliff START 2020.12.18
*-----------------------------------------------------------------------
*        idosunit    in Card[4]-16b
*-----------------------------------------------------------------------

         else if( ipm .eq. 97 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 844

               idosunit = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
*-----------------------------------------------------------------------
*        iangbpwr    in Card[10]-2
*-----------------------------------------------------------------------

         else if( ipm .eq. 98 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 845

               iangbpwr = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
c H.Ratliff START 2020.02.17
*-----------------------------------------------------------------------
*        iphtout    in Card[10]-1
*-----------------------------------------------------------------------

         else if( ipm .eq. 99 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 846

               iphtout = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
c H.Ratliff START 2020.04.07
*-----------------------------------------------------------------------
*        isfylib    in Card[3b]-7
*-----------------------------------------------------------------------

         else if( ipm .eq. 100 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 847

               isfylib = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
*-----------------------------------------------------------------------
*        radlev    in Card[10]-3
*-----------------------------------------------------------------------

         else if( ipm .eq. 101 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 848

               radlev = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
*-----------------------------------------------------------------------
*        ergthfa    in Card[3b]-8
*-----------------------------------------------------------------------

         else if( ipm .eq. 102 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 849

               ergthfa = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
*-----------------------------------------------------------------------
*        ergfafu    in Card[3b]-9
*-----------------------------------------------------------------------

         else if( ipm .eq. 103 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 1850

               ergfafu = cvvv

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200
*-----------------------------------------------------------------------
*        infyerg    in Card[3b]-10
*-----------------------------------------------------------------------

         else if( ipm .eq. 104 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 851

               infyerg = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        hsfylib    in Card[3b]-11
*-----------------------------------------------------------------------

         else if( ipm .eq. 105 ) then

               ict = min( ic + 79, i2 )

               hsfylib = chin(ic:ict)

               hssll = ict - ic + 1

c H.Ratliff START 2021.04.26
*-----------------------------------------------------------------------
*        ilchain    in Card[3]-12
*-----------------------------------------------------------------------

         else if( ipm .eq. 106 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 852

               ilchain = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        mxnuclei
*-----------------------------------------------------------------------
         else if( ipm .eq. 107 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 855

               imxnuc = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        itriton
*-----------------------------------------------------------------------
         else if( ipm .eq. 108 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 853

               itriton = nint( cvvv )

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        iaonucl
*-----------------------------------------------------------------------
         else if( ipm .eq. 109 ) then

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 854

               iaonucl = nint( cvvv )

               if(iaonucl.ge.14) goto 854

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

*-----------------------------------------------------------------------
*        aonucl
*-----------------------------------------------------------------------

         else if( ipm .eq. 110 ) then

               ict = min( ic + 199, i2 )

               aonucl = chin(ic:ict)

               laonucl = ict - ic + 1

*-----------------------------------------------------------------------
*        aoreg
*-----------------------------------------------------------------------

         else if( ipm .eq. 111 ) then

               ict = min( ic + 199, i2 )

               aoreg = chin(ic:ict)

               laoreg = ict - ic + 1

*-----------------------------------------------------------------------
CCSE added r-z/xyz mesh (2018.07.31) <<<<<

         end if

            goto 140

*-----------------------------------------------------------------------
*     check
*-----------------------------------------------------------------------

  800 continue

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

                  filnm(1:lfile(1)) = ifile(1)(1:lfile(1))
                  ifilnm = lfile(1)

                  do j = ifilnm , 1, -1
                   if( filnm(j:j) .eq. '.' ) goto 810
                  end do
                  j = ifilnm + 1

  810             itfp = j - 1
                  ifilnm=itfp

                  do j = ifilnm+1, 100
                   filnm(j:j) = ' '
                  end do

cFURUTA20200601 unnecessary read of dyld and dtrk files were removed

*   for DCHAIN input file
      kase = 0
*        filnm(ifilnm-5:ifilnm) = '.din '

      inquire( file = ifile(1)(1:lfile(1)), exist = lex )

         if( lex ) then
         else
         end if

      open(iodw, file = ifile(1)(1:lfile(1)), status = 'unknown' )

      if ( iversion .eq. 0 ) then
*   Card[1]
        write(iodw,'(4x,"htitle = ",80a)') title(1:titll)
*   Card[2]
        write(iodw,'(i2,1x,i2)') imode, jmode
*   Card[3]
        write(iodw,'(i6,1x,i6,1x,i6,1x,i3,1x,i2,1x,i6,1x,
     &               i2,1x,i2,1x,i2,1x,i2,1x,i2)')
*    &        idcstp, idcout, idivs, iregon, inmtcf, ichain,
     &        idcstp, idcout, idivs, klst,   inmtcf, ichain,
     &        itdecs, itdecn, isomtr, ifisyd, ifisye
*   Card[4]
        write(iodw,'(i2,1x,i2,1x,i2,1x,1pe11.4,1x,i2,1x,i2,1x,
     &               i2,1x,i2,1x,1pe11.4,1x,i6,1x,i2)')
     &        iyild, iggrp, ibetap, acmin, istabl, igsdef,
     &        iprtb1, iprtb2, rprtb2, iprtb3, igsorg
*   Card[5]
c iwamoto 2012/12/10 (source/sec) -> (mA)
              ampp = amp*1.602177d-19/1.0d-3
*    &        ebeam, prodnp
        write(iodw,'(1pe11.4,1x,1pe11.4,1x,1pe11.4)')
     &        ampp, ebeam, prodnp
*   Card[6]
        do i = 1, idcstp
           if( itemsh(1,i) .eq. 1 ) then
              chau = 'y'
           else if( itemsh(1,i) .eq. 2 ) then
              chau = 'd'
           else if( itemsh(1,i) .eq. 3 ) then
              chau = 'h'
           else if( itemsh(1,i) .eq. 4 ) then
              chau = 'm'
           else if( itemsh(1,i) .eq. 5 ) then
              chau = 's'
           end if
           write(iodw,'(5x,1pe11.4,1x,a1,1x,1pe11.4)')
     &           temsh(1,i), chau, temsh(2,i)
        end do
*   Card[7]
        do i = 1, idcout
           if( itemsh(2,i) .eq. 1 ) then
              chau = 'y'
           else if( itemsh(2,i) .eq. 2 ) then
              chau = 'd'
           else if( itemsh(2,i) .eq. 3 ) then
              chau = 'h'
           else if( itemsh(2,i) .eq. 4 ) then
              chau = 'm'
           else if( itemsh(2,i) .eq. 5 ) then
              chau = 's'
           end if
           write(iodw,'(5x,1pe11.4,1x,a1)') temsh(3,i), chau
        end do
*   Card[8]
        if( iancza .ne. 0 ) then
          do i = 1, iancza
             write(iodw,'(2x,a8)') hancza(i)
          end do
        end if
*   Card[9]
*   Card[10]
  895 continue
        kase = kase + 1

           chaia3(1:3) = '000'
        if( kase .lt. 10 ) then
           write(chaia1,'(i1)') kase
           chaia3(1:3) = '00'//chaia1
        else if( kase .lt. 100 ) then
           write(chaia2,'(i2)') kase
           chaia3(1:3) = '0'//chaia2
        else if( kase .lt. 1000 ) then
           write(chaia3,'(i3)') kase
           chaia3(1:3) = chaia3
        else
           goto 999
        end if
c amp(source/sec)
            if (ivoll .eq. 1 .or. itgvoll .eq. 1 ) then
        flux = amp  * fflux(itdc2reg(itnm2tdc(m-2))+kase) !FURUTA20200601
            end if
            if (ivoll .eq. 0 .and. itgvoll .eq. 0 ) then
        flux = amp  * fflux(itdc2reg(itnm2tdc(m-2))+kase) / rtgvol(kase) !FURUTA20200601
            end if

        write(iodw,'("DUMMY",a,1x,i8,1x,i4," 0 0",
     &               " 0 0 ",1pe11.4," n.flux_",a,1x,1pe11.4)')
     &        chaia3(1:3), itgcell(kase), itglist(kase),
     &        flux, chaia3(1:3), rtgvol(kase)

*   Card[10a]
        do i = 1, itglist(kase)

         k0=jreg2itm(i)         !FURUTA20200522
         write(iodw,'(2x,a8,2x,1pe11.4)')
     &        htgnzas(k0+kase), ctgnnds(k0+kase)
*    &                            tgnzas(kase,i), ctgnnds(kase,i)
        end do

        if( klst - 1 .ge. kase ) goto 895
        close(iodw)
      else  ! if ( iversion .eq. 1 ) ! (the modern version of DCHAIN input)

*   Card[1]
        write(iodw,'(4x,"htitle = ",80a)') title(1:titll)
*   Card[2]
        write(iodw,'(" ")')
        write(iodw,'("! --- control parameters ---")')
        write(iodw,'(5x,"imode = ",i7)') imode
        write(iodw,'(5x,"jmode = ",i7)') jmode
*   Card[3]
        write(iodw,'(" ")')
        write(iodw,'("! --- calculation parameters ---")')
        write(iodw,'(5x,"idivs = ",i7)') idivs
        if (irdfisyd.ne.1 .and. maxval(tgnzas(:)).gt.86220) ifisyd = 2  ! H.Ratliff 2021.05.19 change default ifisyd if actinides present
c H.Ratliff START 2021.04.26 increase default ichain if fission is likely
        if (ichain.eq.100 .and. ifisyd.gt.0) then ! if ichain still at normal default value
           if (maxval(tgnzas(:)).gt.86220) ichain = 500 ! new default if actinides present
        endif
        write(iodw,'(4x,"ichain = ",i7)') ichain
        write(iodw,'(3x,"ilchain = ",i7)') ilchain
c H.Ratliff END 2021.04.26
        write(iodw,'(4x,"itdecs = ",i7)') itdecs
        write(iodw,'(4x,"itdecn = ",i7)') itdecn
c H.Ratliff START 2021.01.12 adjust isomer treatment based on igamma value
        if (igamma.eq.3 .and. isomtr.eq.2) then   ! if user is producing individual isomers,
           write(iodw,'(4x,"isomtr = ",i7)') 0  ! they should be tracked individually in DCHAIN too.
        else
           write(iodw,'(4x,"isomtr = ",i7)') isomtr
        endif
        write(iodw,'(4x,"ifisyd = ",i7)') ifisyd
c H.Ratliff START 2020.04.17
        write(iodw,'(3x,"ixsrall = ",i7)') ixsrall
c H.Ratliff START 2019.07.30, add data library selection parameters
*   Card[3b]
        write(iodw,'(" ")')
        write(iodw,'("! --- data library parameters ---")')
        write(iodw,'(3x,"inxslib = ",i7)') inxslib
        write(iodw,'(3x,"idcylib = ",i7)') idcylib
        write(iodw,'(3x,"infylib = ",i7)') infylib
        if (infylib.ne.0) then
           write(iodw,'(3x,"infyerg = ",i7)') infyerg
           if (infyerg.eq.4) then
              write(iodw,'(3x,"ergthfa = ",1pe11.4)') ergthfa
              write(iodw,'(3x,"ergfafu = ",1pe11.4)') ergfafu
           endif
        endif
        write(iodw,'(3x,"isfylib = ",i7)') isfylib
        if (inxslib.eq.-1) then
           write(iodw,'(3x,"hnxslib = ",a)') hnxslib
        else
           if (inxslib.eq.0.or.inxslib.eq.40) then
              hnxslib = ' '
           elseif (inxslib.eq.1) then
              hnxslib = 'Legacy2001'
           elseif (inxslib.eq.2.or.inxslib.eq.20) then
              hnxslib = 'JENDL-AD17'
           elseif (inxslib.eq.21) then
              hnxslib = 'JENDL-4-0-'
           elseif (inxslib.eq.30) then
              hnxslib = 'ENDF-B-7-1'
           elseif (inxslib.eq.31) then
              hnxslib = 'ENDF-B-8-0'
           elseif (inxslib.eq.41) then
              hnxslib = 'JEFF-3-3--'
           elseif (inxslib.eq.50) then
              hnxslib = 'FENDL-A-30'
           elseif (inxslib.eq.51) then
              hnxslib = 'EAF-2010--'
           elseif (inxslib.eq.60) then
              hnxslib = 'BROND-3-1-'
           elseif (inxslib.eq.70) then
              hnxslib = 'CENDL-3-1-'
           elseif (inxslib.eq.90) then
              hnxslib = 'TENDL-2017'
           elseif (inxslib.eq.100) then
              hnxslib = 'h-JP-US-EU'
           elseif (inxslib.eq.101) then
              hnxslib = 'h-US-JP-EU'
           elseif (inxslib.eq.102) then
              hnxslib = 'h-EU-US-JP'
           elseif (inxslib.eq.22) then ! Add by Condar, T.Sato 2024/11/19
              hnxslib = 'JENDL-5'
           elseif (inxslib.eq.120) then
              hnxslib = 'h-JENDLAD4'
           endif
           write(iodw,'("!",3x,"hnxslib = ",a)') hnxslib
        endif
        if (idcylib.eq.-1) then
           write(iodw,'(3x,"hdcylib = ",a)') hdcylib
        else
           if (idcylib.eq.0) then
              hdcylib = ' '
           elseif (idcylib.eq.1) then
              hdcylib = 'Legacy2001'
           elseif (idcylib.eq.2) then
              hdcylib = 'JENDLDDF15'
           elseif (idcylib.eq.3) then
              hdcylib = 'ENDFBVIII0'
           elseif (idcylib.eq.4) then
              hdcylib = 'EB8-J15-ES'
           elseif (idcylib.eq.5) then
              hdcylib = 'J15-EB8-ES'
           endif
           write(iodw,'("!",3x,"hdcylib = ",a)') hdcylib
        endif
        if (infylib.eq.-1) then
           write(iodw,'(3x,"hnfylib = ",a)') hnfylib
        elseif (infylib.ne.0) then
           if (infylib.eq.10) then
              hnfylib = 'endf6'
           elseif (infylib.eq.11) then
              hnfylib = 'jndcfp2'
           elseif (infylib.eq.12) then
              hnfylib = 'endfb8'
           elseif (infylib.eq.13) then
              hnfylib = 'jendlfpy2011'
           elseif (infylib.eq.14) then
              hnfylib = 'jeff33'
           elseif (infylib.eq.16) then
              hnfylib = 'h-ENDF8-JEFF'
           elseif (infylib.eq.17) then
              hnfylib = 'h-JENDL-JEFF'
           elseif (infylib.eq.18) then
              hnfylib = 'h-JEFF-ENDF8'
           elseif (infylib.eq.19) then
              hnfylib = 'h-JEFF-JENDL'
           endif
           write(iodw,'("!",3x,"hnfylib = ",a)') hnfylib
        endif
        if (isfylib.eq.-1) then
           write(iodw,'(3x,"hsfylib = ",a)') hsfylib
        elseif (isfylib.eq.1) then
           if (infylib.eq.-1) then
              write(iodw,'(3x,"hsfylib = ",a)') hnfylib
           elseif (infylib.ne.0) then
              write(iodw,'("!",3x,"hsfylib = ",a)') hnfylib
           endif
        elseif (isfylib.ne.0) then
           if (isfylib.eq.10) then
              hsfylib = 'endf6'
           elseif (isfylib.eq.11) then
              hsfylib = 'jndcfp2'
           elseif (isfylib.eq.12) then
              hsfylib = 'endfb8'
           elseif (isfylib.eq.13) then
              hsfylib = 'jendlfpy2011'
           elseif (isfylib.eq.14) then
              hsfylib = 'jeff33'
           elseif (isfylib.eq.16) then
              hsfylib = 'h-ENDF8-JEFF'
           elseif (isfylib.eq.17) then
              hsfylib = 'h-JENDL-JEFF'
           elseif (isfylib.eq.18) then
              hsfylib = 'h-JEFF-ENDF8'
           elseif (isfylib.eq.19) then
              hsfylib = 'h-JEFF-JENDL'
           endif
           write(iodw,'("!",3x,"hsfylib = ",a)') hsfylib
        endif
*   Card[4]
        write(iodw,'(" ")')
        write(iodw,'("! --- output parameters ---")')
        write(iodw,'(5x,"iyild = ",i7)') iyild
        write(iodw,'(5x,"iggrp = ",i7)') iggrp
        write(iodw,'(4x,"ibetap = ",i7)') ibetap
        write(iodw,'(5x,"acmin = ",1pe11.4)') acmin
        write(iodw,'(4x,"istabl = ",i7)') istabl
        write(iodw,'(4x,"igsdef = ",i7)') igsdef
        write(iodw,'(4x,"igsorg = ",i7)') igsorg
c H.Ratliff START 2019.07.30, add decay chain output file control parameters
        write(iodw,'(3x,"iwrtchn = ",i7)') iwrtchn
        write(iodw,'(3x,"chrlvth = ",1pe11.4)') chrlvth
        write(iodw,'(3x,"iwrchdt = ",i7)') iwrchdt
        write(iodw,'(3x,"iwrchss = ",i7)') iwrchss
c H.Ratliff START 2020.12.18
        write(iodw,'(2x,"iwrchnuc =       0")')  ! only changable in DCHAIN due to advanced usage
c H.Ratliff START 2020.04.17
        write(iodw,'(3x,"idosecf = ",i7)') idosecf
c H.Ratliff START 2020.12.18
        write(iodw,'(2x,"idosunit = ",i7)') idosunit
c H.Ratliff START 2020.11.16
        write(iodw,'(3x,"imtcard = ",i7)') imtcard
        write(iodw,'(3x,"imtcnum = ",i7)') imtcnum
        write(iodw,'(2x,"imtcmeta = ",i7)') imtcmeta
        write(iodw,'(3x,"thmatnd = ",1pe11.4)') thmatnd
*   Card[5]
c iwamoto 2012/12/10 (source/sec) -> (mA)
            ampp = amp*1.602177d-19/1.0d-3
        write(iodw,'(" ")')
        write(iodw,'("! --- Proton beam current and",
     &               " neutron flux ---")')
        write(iodw,'(7x,"amp = ",1pe11.4)') ampp
        write(iodw,'(5x,"ebeam = ",1pe11.4)') ebeam
        write(iodw,'(4x,"prodnp = ",1pe11.4)') prodnp
*   Card[6]
        write(iodw,'(" ")')
        write(iodw,'("! --- irradiation time ---")')
        write(iodw,'(4x,"itstep = ",i7)') idcstp
        do i = 1, idcstp
           if( itemsh(1,i) .eq. 1 ) then
              chau = 'y'
           else if( itemsh(1,i) .eq. 2 ) then
              chau = 'd'
           else if( itemsh(1,i) .eq. 3 ) then
              chau = 'h'
           else if( itemsh(1,i) .eq. 4 ) then
              chau = 'm'
           else if( itemsh(1,i) .eq. 5 ) then
              chau = 's'
           end if
           write(iodw,'(5x,1pe11.4,1x,a1,1x,1pe11.4)')
     &           temsh(1,i), chau, temsh(2,i)
        end do
*   Card[7]
        write(iodw,'(" ")')
        write(iodw,'("! --- output time ---")')
        write(iodw,'(5x,"itout = ",i7)') idcout
        do i = 1, idcout
           if( itemsh(2,i) .eq. 1 ) then
              chau = 'y'
           else if( itemsh(2,i) .eq. 2 ) then
              chau = 'd'
           else if( itemsh(2,i) .eq. 3 ) then
              chau = 'h'
           else if( itemsh(2,i) .eq. 4 ) then
              chau = 'm'
           else if( itemsh(2,i) .eq. 5 ) then
              chau = 's'
           end if
           write(iodw,'(5x,1pe11.4,1x,a1)') temsh(3,i), chau
        end do
C   Card[8] is disabled function for NEW DCHAIN.
*   Card[9]
        write(iodw,'(" ")')
        write(iodw,'("! --- irradiation condition ---")')
        write(iodw,'(4x,"inmtcf =       1")')
        write(iodw,'("!  ffnmtcf =    1.0")')
CCSE change for mntc_yield file (2017.12.31) >>>>>
        hhnmtcf(1:) = filnm(1:ifilnm)//'.dyld' !FURUTA20200601
        write(iodw,'("  itdchout =       1"/
     &               "   hhnmtcf = ",a)') trim(hhnmtcf)
CCSE change for mntc_yield file (2017.12.31) <<<<<
c H.Ratliff START 2019.07.30 added DCHAIN uncertainty propagation options
        if (iertdcho.eq.1) then
           henmtcf(1:) = filnm(1:ifilnm)//'_err.dyld' !FURUTA20200601
           write(iodw,'("  iertdcho =       1"/
     &               "   henmtcf = ",a)') trim(henmtcf)
        else
           henmtcf(1:) = filnm(1:ifilnm)//'_err.dyld' !FURUTA20200601
           write(iodw,'("  iertdcho =       0"/
     &               "!   henmtcf = ",a)') trim(henmtcf)
        endif
c H.Ratliff START 2020.06.05
        write(iodw,'(2x,"iredufmt = ",i7)') iredufmt0
c H.Ratliff START 2020.06.17
        write(iodw,'(3x,"irdonce = ",i7)') irdonce
c H.Ratliff START 2020.05.29
        if( meshtype .eq. 4 ) then
          write(iodw,'("! --- special tet mesh output options ---")')
           write(iodw,'(3x,"foamout = ",i7)') foamout
           write(iodw,'(2x,"foamvals = ",i7)') foamvals
        endif
c H.Ratliff START 2020.04.17
        if( meshtype .eq. 2 .or. meshtype .eq. 3 ) then
          write(iodw,'("! --- special xyz mesh output options ---")')
           write(iodw,'(2x,"ipltmode = ",i7)') ipltmode
           write(iodw,'(2x,"ipltaxis = ",i7)') ipltaxis
        endif
        write(iodw,'(" ")')
        write(iodw,'("! --- calculation region details ---")')

CCSE add for r-z/xyz mesh(2018.07.31) >>>>>
        if ( meshtype .eq. 1 ) then
CCSE add for r-z/xyz mesh(2018.07.31) <<<<<

*   Card[10]
           write(iodw,'(" ")')
           write(iodw,'(4x,"iregon = ",i6)') klst
  897 continue
           kase = kase + 1

           chaia3(1:3) = '000'
           if( kase .lt. 10 ) then
              write(chaia1,'(i1)') kase
              chaia3(1:3) = '00'//chaia1
           else if( kase .lt. 100 ) then
              write(chaia2,'(i2)') kase
              chaia3(1:3) = '0'//chaia2
           else if( kase .lt. 1000 ) then
              write(chaia3,'(i3)') kase
              chaia3(1:3) = chaia3
           else
              goto 999
           end if

c amp(source/sec)
           if (ivoll .eq. 1 .or. itgvoll .eq. 1 ) then
              flux = amp  * fflux(itdc2reg(itnm2tdc(m-2))+kase) !FURUTA20200601
           end if
           if (ivoll .eq. 0 .and. itgvoll .eq. 0 ) then
            flux = amp  * fflux(itdc2reg(itnm2tdc(m-2))+kase) /
     &           rtgvol(kase)   !FURUTA20200601
           end if

           write(iodw,'("!1)HRGCMM 2)IREGS 3)ITGNCLS 4)FLUXS",
     &                   " 5)HNFLUXS 6)VOLUMES")')
           if(flux .eq. 0.d0) then
               flux = 1.d-30
        write(*,'("Warning: No neutron flux in region",1x,i8,
     & ". Results of T-dchain may be unreliable unless statistics is
     & improved")')
     & itgcell(kase)
           endif
CCSE change for n.flux_ file (2017.12.31) >>>>>
           hnfluxs(1:) = filnm(1:ifilnm)//'.dtrk' !FURUTA20200601
           write(iodw,'("   DUMMY",a,1x,i8,1x,i4,1x,1pe11.4,
     &                  "   ",a,1x,1pe11.4)')
     &        chaia3(1:3), itgcell(kase), itglist(kase),
     &        flux, trim(hnfluxs), rtgvol(kase)
CCSE change for n.flux_ file (2017.12.31) <<<<<
*   Card[10a]
           k0=jreg2itm(kase)       !FURUTA20200522
           do i = 1, itglist(kase)

            write(iodw,'(2x,a8,2x,1pe11.4)') !FURUTA20200522
     &       htgnzas(k0+i), ctgnnds(k0+i)    !FURUTA20200522

           end do

           if( klst - 1 .ge. kase ) goto 897

CCSE add for r-z/xyz mesh(2018.07.31) >>>>>
        ! branch for mesh type ... under construction
        elseif( meshtype .ge. 2 .and. meshtype .le. 4 ) then !FURUTA20200505
*   Card[10]
           write(iodw,'(" ")')

           write(iodw,'("      imat =",i7)') imat
           write(iodw,'("   mtscore =",i7)') mtscore

           do ik = 1, klst
              write(iodw,'("   mt-list =",i7)') itglist(ik)
              j1=jreg2itm(ik) !FURUTA20200522
              do il = 1, itglist(ik)
                 write(iodw,'(2x,a8,2x,1pe11.4)')
     &              htgnzas(j1+il), ctgnnds(j1+il) !FURUTA20200522
              enddo
           enddo

           write(iodw,'("     imesh = ",i7)') imesh
           write(iodw,'(a)') chlst(1:ichlst)

ccse 2023/03/29, add write for debug

           do i = 1, imesh

              if(meshtype.eq.4)then
                flux = amp  * fflux(itdc2reg(itnm2tdc(m-2))+i) !FURUTA20250123
              elseif (ivoll .eq. 1 .or. itgvoll .eq. 1 ) then
                flux = amp  * fflux(itdc2reg(itnm2tdc(m-2))+i) !FURUTA20200630
              elseif (ivoll .eq. 0 .and. itgvoll .eq. 0 ) then
                flux = amp  * fflux(itdc2reg(itnm2tdc(m-2))+i) /
     &               rtgvol(i)  !FURUTA20200601
              end if

              if( meshtype .eq. 2 ) then
                 write(iodw,'(3i6,1p20e13.4)')
     &            (ildata(i,j),j=1,iclmesh),flux,
     &            (rldata(i,j),j=1,iclvol)
              elseif( meshtype .eq. 3 ) then
                 write(iodw,'(4i6,1p20e13.4)')
     &            (ildata(i,j),j=1,iclmesh),flux,
     &            (rldata(i,j),j=1,iclvol)
              elseif( meshtype .eq. 4 ) then !FURUTA20200505
                 write(iodw,'(3i6,1p20e13.4)')
     &              (ildata(i,j),j=1,iclmesh),    !FURUTA20200529
     &              flux,(rldata(i,j),j=1,iclvol) !FURUTA20200529
              endif
           enddo

! deallocate
           deallocate( ildata )
           deallocate( rldata )



        endif
CCSE add for r-z/xyz mesh(2018.07.31) <<<<<

*   Card[11] NEW! (card 10 in updated DCHAIN manual)
        write(iodw,'(" ")')
        write(iodw,'("! --- output option for ANGEL ---")')
c H.Ratliff START 2020.12.18
        write(iodw,'(2x,"phitsout = ",i3)')  iphtout
        write(iodw,'(2x,"itriton  = ",i3)')  itriton
        write(iodw,'(4x,"radlev = ",1pe11.4)') radlev
c H.Ratliff START 2020.12.18
        write(iodw,'(2x,"iangbpwr = ",i2)')  iangbpwr
        write(iodw,'(2x,"angelout =  1")')
        write(iodw,'("!   angelout_act        = 0")')
        write(iodw,'("!   angelout_heat       = 0")')
        write(iodw,'("!   angelout_heat_alpha = 0")')
        write(iodw,'("!   angelout_heat_beta  = 0")')
        write(iodw,'("!   angelout_heat_gamma = 0")')
        write(iodw,'("!   angelout_dose       = 0")')
      if(laoreg.eq.0) then ! aoreg is not specified
        write(iodw,'("!   angelout_region     = all, sum")')
      else
        write(iodw,'("    angelout_region = ",a)') aoreg(1:laoreg)
      endif
c H.Ratliff START 2020.12.18
      if(laonucl.eq.0) then ! aonucl is not specified
        write(iodw,'("!   angelout_nuclides   = 0")')  ! only changable in DCHAIN due to advanced usage
      else
        write(iodw,'("    angelout_nuclides   =",i3)') iaonucl
        write(iodw,'(a)') aonucl(1:laonucl)
      endif
        write(iodw,'(" ")')
        write(iodw,'("end")')

      end if
      close(iodw)

      goto 999

*-----------------------------------------------------------------------
*     error
*-----------------------------------------------------------------------
CCSE add for r-z/xyz mesh(2018.07.31) >>>>>
  920    m_err = 'Unknown mesh parameter in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:10600/R:tdctoinp/F:tallsm4.f'
         goto 998
CCSE add for r-z/xyz mesh(2018.07.31) <<<<<

  921    m_err = 'KeyWord "[T-Dchain]" is not found in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:10606/R:tdctoinp/F:tallsm4.f'
         goto 998

  922    m_err = 'imode (Card[2]-1) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10611/R:tdctoinp/F:tallsm4.f'
         goto 998
  923    m_err = 'jmode (Card[2]-2) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10615/R:tdctoinp/F:tallsm4.f'
         goto 998

  924    m_err = 'itstep (Card[3]-1) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10620/R:tdctoinp/F:tallsm4.f'
         goto 998
  925    m_err = 'itout (Card[3]-2) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10624/R:tdctoinp/F:tallsm4.f'
         goto 998
  926    m_err = 'idivs (Card[3]-3) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10628/R:tdctoinp/F:tallsm4.f'
         goto 998
  927    m_err = 'iregon (Card[3]-4) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10632/R:tdctoinp/F:tallsm4.f'
         goto 998
  928    m_err = 'inmtcf (Card[3]-5) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10636/R:tdctoinp/F:tallsm4.f'
         goto 998
  929    m_err = 'ichain (Card[3]-6) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10640/R:tdctoinp/F:tallsm4.f'
         goto 998
  930    m_err = 'itdecs (Card[3]-7) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10644/R:tdctoinp/F:tallsm4.f'
         goto 998
  931    m_err = 'itdecn (Card[3]-8) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10648/R:tdctoinp/F:tallsm4.f'
         goto 998
  932    m_err = 'isomtr (Card[3]-9) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10652/R:tdctoinp/F:tallsm4.f'
         goto 998
  933    m_err = 'ifisyd (Card[3]-10) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10656/R:tdctoinp/F:tallsm4.f'
         goto 998
  934    m_err = 'ifisye (Card[3]-11) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10660/R:tdctoinp/F:tallsm4.f'
         goto 998

  935    m_err = 'iyild (Card[4]-1) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10665/R:tdctoinp/F:tallsm4.f'
         goto 998
  936    m_err = 'iggrp (Card[4]-2) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10669/R:tdctoinp/F:tallsm4.f'
         goto 998
  937    m_err = 'ibetap (Card[4]-3) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10673/R:tdctoinp/F:tallsm4.f'
         goto 998
  938    m_err = 'acmin (Card[4]-4) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10677/R:tdctoinp/F:tallsm4.f'
         goto 998
  939    m_err = 'istabl (Card[4]-5) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10681/R:tdctoinp/F:tallsm4.f'
         goto 998
  940    m_err = 'igsdef (Card[4]-6) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10685/R:tdctoinp/F:tallsm4.f'
         goto 998
  941    m_err = 'iprtb1 (Card[4]-7) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10689/R:tdctoinp/F:tallsm4.f'
         goto 998
  942    m_err = 'iprtb2 (Card[4]-8) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10693/R:tdctoinp/F:tallsm4.f'
         goto 998
  943    m_err = 'rprtb2 (Card[4]-9) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10697/R:tdctoinp/F:tallsm4.f'
         goto 998
  944    m_err = 'iprtb3 (Card[4]-10) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10701/R:tdctoinp/F:tallsm4.f'
         goto 998
  945    m_err = 'igsorg (Card[4]-11) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10705/R:tdctoinp/F:tallsm4.f'
         goto 998

  946    m_err = 'amp (Card[5]-1) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10710/R:tdctoinp/F:tallsm4.f'
         goto 998
  947    m_err = 'ebeam (Card[5]-2) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10714/R:tdctoinp/F:tallsm4.f'
         goto 998
  948    m_err = 'prodnp (Card[5]-3) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10718/R:tdctoinp/F:tallsm4.f'
         goto 998

  949    m_err = 'ac-list (Card[8]) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10723/R:tdctoinp/F:tallsm4.f'
         goto 998
  950    m_err = 'tg-list (Card[10]) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10727/R:tdctoinp/F:tallsm4.f'
         goto 998

  951    m_err = 'mesh selection, allowed "1", is wrong.'
         ErrCha = ''
         ErrID = 'L:10732/R:tdctoinp/F:tallsm4.f'
         goto 998
  952    m_err = 'mesh description of Card[6] is wrong.'
         ErrCha = ''
         ErrID = 'L:10736/R:tdctoinp/F:tallsm4.f'
         goto 998
  953    m_err = 'mesh description of Card[7] is wrong.'
         ErrCha = ''
         ErrID = 'L:10740/R:tdctoinp/F:tallsm4.f'
         goto 998
  954    m_err = 'description or number of Card[8] is wrong.'
         ErrCha = ''
         ErrID = 'L:10744/R:tdctoinp/F:tallsm4.f'
         goto 998
  955    m_err = 'description or number of Card[10] is wrong.'
         ErrCha = ''
         ErrID = 'L:10748/R:tdctoinp/F:tallsm4.f'
         goto 998
  956    m_err = 'description of density (Card[10]) is wrong.'
         ErrCha = ''
         ErrID = 'L:10752/R:tdctoinp/F:tallsm4.f'
         goto 998
  957    m_err = 'description of version OLD(0)/NEW(1) is wrong.'
         ErrCha = ''
         ErrID = 'L:10756/R:tdctoinp/F:tallsm4.f'
         goto 999

  959    m_err = 'KeyWord "[T-Yield]" is not found in tally [t-yield]'
         ErrCha = ''
         ErrID = 'L:10761/R:tdctoinp/F:tallsm4.f'
         goto 999
  960    m_err = 'KeyWord is not found in tally [t-yield]'
         ErrCha = ''
         ErrID = 'L:10765/R:tdctoinp/F:tallsm4.f'
         goto 998

  961    m_err = 'KeyWord "[T-Track]" is not found in tally [t-track]'
         ErrCha = ''
         ErrID = 'L:10770/R:tdctoinp/F:tallsm4.f'
         goto 999
  962    m_err = 'KeyWord is not found in tally [t-track]'
         ErrCha = ''
         ErrID = 'L:10774/R:tdctoinp/F:tallsm4.f'
         goto 998


  987    m_err = 'Unknown parameter in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:10780/R:tdctoinp/F:tallsm4.f'
         goto 998

  989    m_err = 'Double definition of the parameter in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:10785/R:tdctoinp/F:tallsm4.f'
         goto 998

  990    m_err = 'Too many file name in tally '//tnamed
         ErrCha = ''
         ErrID = 'L:10790/R:tdctoinp/F:tallsm4.f'
         goto 998

  997    m_err = 'parameter is wrong in tally '//tnamed//' '//chlw
         ErrCha = ''
         ErrID = 'L:10795/R:tdctoinp/F:tallsm4.f'
         goto 998

  821    m_err = 'iertdcho (Card[9]-5) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10800/R:tdctoinp/F:tallsm4.f'
         goto 998

  822    m_err = 'inxslib (Card[3b]-1) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10805/R:tdctoinp/F:tallsm4.f'
         goto 998

  823    m_err = 'idcylib (Card[3b]-2) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10810/R:tdctoinp/F:tallsm4.f'
         goto 998

  824    m_err = 'infylib (Card[3b]-3) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10815/R:tdctoinp/F:tallsm4.f'
         goto 998

  825    m_err = 'hnxslib (Card[3b]-4) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10820/R:tdctoinp/F:tallsm4.f'
         goto 998

  826    m_err = 'hdcylib (Card[3b]-5) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10825/R:tdctoinp/F:tallsm4.f'
         goto 998

  827    m_err = 'hnfylib (Card[3b]-6) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10830/R:tdctoinp/F:tallsm4.f'
         goto 998

  828    m_err = 'iwrtchn (Card[4]-12) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10835/R:tdctoinp/F:tallsm4.f'
         goto 998

  829    m_err = 'chrlvth (Card[4]-13) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10840/R:tdctoinp/F:tallsm4.f'
         goto 998

 2830    m_err = 'iwrchdt (Card[4]-14) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10845/R:tdctoinp/F:tallsm4.f'
         goto 998

  831    m_err = 'iwrchss (Card[4]-15) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10850/R:tdctoinp/F:tallsm4.f'
         goto 998

  832    m_err = 'ixsrall (Card[3]-7) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10855/R:tdctoinp/F:tallsm4.f'
         goto 998

  833    m_err = 'idosecf (Card[4]-16) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10860/R:tdctoinp/F:tallsm4.f'
         goto 998

  834    m_err = 'ipltmode (Card[9d]-1) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10865/R:tdctoinp/F:tallsm4.f'
         goto 998

  835    m_err = 'ipltaxis (Card[9d]-2) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10870/R:tdctoinp/F:tallsm4.f'
         goto 998

  836    m_err = 'foamout (Card[9d]-1) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10875/R:tdctoinp/F:tallsm4.f'
         goto 998

  837    m_err = 'foamvals (Card[9d]-2) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10880/R:tdctoinp/F:tallsm4.f'
         goto 998

  838    m_err = 'iredufmt (Card[8]-7) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10885/R:tdctoinp/F:tallsm4.f'
         goto 998

  839    m_err = 'irdonce  (Card[8]-8) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10890/R:tdctoinp/F:tallsm4.f'
         goto 998

  840    m_err = 'imtcard  (Card[4]-17) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10895/R:tdctoinp/F:tallsm4.f'
         goto 998

  841    m_err = 'imtcnum  (Card[4]-18) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10900/R:tdctoinp/F:tallsm4.f'
         goto 998

  842    m_err = 'imtcmeta  (Card[4]-19) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10905/R:tdctoinp/F:tallsm4.f'
         goto 998

  843    m_err = 'thmatnd  (Card[4]-20) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10910/R:tdctoinp/F:tallsm4.f'
         goto 998

  844    m_err = 'idosunit  (Card[4]-16b) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10915/R:tdctoinp/F:tallsm4.f'
         goto 998

  845    m_err = 'iangbpwr  (Card[10]-2) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10920/R:tdctoinp/F:tallsm4.f'
         goto 998

  846    m_err = 'iphtout   (Card[10]-1) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10925/R:tdctoinp/F:tallsm4.f'
         goto 998

  847    m_err = 'isfylib   (Card[3b]-7) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10930/R:tdctoinp/F:tallsm4.f'
         goto 998

  848    m_err = 'radlev    (Card[10]-3) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10935/R:tdctoinp/F:tallsm4.f'
         goto 998

  849    m_err = 'ergthfa   (Card[3b]-8) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10940/R:tdctoinp/F:tallsm4.f'
         goto 998

 1850    m_err = 'ergfafu   (Card[3b]-9) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10945/R:tdctoinp/F:tallsm4.f'
         goto 998

  851    m_err = 'infyerg   (Card[3b]-10) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10950/R:tdctoinp/F:tallsm4.f'
         goto 998

  852    m_err = 'ilchain   (Card[3]-12) value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10955/R:tdctoinp/F:tallsm4.f'
         goto 998

  853    m_err = 'itriton value is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10960/R:tdctoinp/F:tallsm4.f'
         goto 998

  854    m_err = 'iaonucl value (<14) is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10965/R:tdctoinp/F:tallsm4.f'
         goto 998

  855    m_err = 'mxnuclei is wrong or unreadable.'
         ErrCha = ''
         ErrID = 'L:10970/R:tdctoinp/F:tallsm4.f'
         goto 998


*-----------------------------------------------------------------------

  998 continue
         call ErrWrite(ErrID, ErrCha)
         write(iot,'(" ERROR = ",200a1)') ( m_err(i:i), i=1,200 )
  999 continue
         close(iodr)
         close(iodw)

         l_err = illd(jsn)
         k_err = jsn
         ierr  = 1

*-----------------------------------------------------------------------
      end subroutine

************************************************************************
*                                                                      *
      subroutine pdchreg2(m)
c      original program made by N.Matsuda on 2012/08/01                *
c      see    subroutine pdchreg(m,mz,mn,                              *
c     &                   nr,mr,nn,kr,nt,ikzz,iknn,                    *
c     &                   tr,tm,vl,lr,nvl,ivl,rvl,                     *
c     &                   nx,ny,nz,xm,ym,zm,val,ixyz,igsh,idasa)       *
*                                                                      *
*       output of dchain tally in region mesh                          *
*       last modified by Y. Iwamoto on 2012/12/05                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
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

      common /tall21/ rtfac(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)

*-----------------------------------------------------------------------

      common /volreg/ dvol(kvlmax)
      common /regdc/ idrg(kvlmax), idgr(kvmmax)

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
***       common /jcomon/ nabov,nobch,nocas,nomax
*** !$OMP THREADPRIVATE(/jcomon/)
      common /talout/ itall

      character fname*100, fnume*3

*-----------------------------------------------------------------------

      character hsunit(2)*15

      data hsunit / '[1/source]     ',
     &              '[1/cm^3/source]'/

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

      data ipstep / 12 /

      character erfnm*100

      character chau*8
      character cha*1
      data cha /"'"/

      character rpa*1
      data rpa /'}'/
      character yen*1
      yen  = char(92)

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do 900 iax = 1, itfln(m)

                  fname = ctfln(m,iax)(1:itfll(m,iax))
                  ifname = itfll(m,iax)
***  "RENAME" modified by norihiro MATSUDA  2012.07.26  ----------------

                  itfp = 0
               do j = ifname, 1, -1
                  if( fname(j:j) .eq. '.' ) goto 810
               end do
                  j = ifname + 1

  810             itfp = j - 1

                  fname = fname(1:itfp)//'.dout'
                  ifname = itfp + 6

               do j = ifname, 100
                  fname(j:j) = ' '
               end do

*------------------------------------------------------- "RENAME" end --

            fname = fname

            iot = 31

            open(iot, file = fname, status = 'unknown' )

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------
               call tdchech(iot,m,iax,1)

*-----------------------------------------------------------------------
*        dchain axis
*-----------------------------------------------------------------------
*
*-----------------------------------------------------------------------
*        output : do not change below expresion
*-----------------------------------------------------------------------

            write(iot,'(/"#",78("-"))')

*-----------------------------------------------------------------------

            close(iot)
  900 continue

*-----------------------------------------------------------------------
      call tdctoinp(m,fname) !FURUTA20200601
*-----------------------------------------------------------------------

      return
      end subroutine




************************************************************************
*                                                                      *
      subroutine echrgdc2(icc,ierl,iva,ivs,iot,mc,kc,nr,mr,kr,
     &                 vl,lr,nvl,ivl,rvl,igm,m)
*                                                                      *
*       MPI test  for tally region mesh                                *
*       2012/12/5                Y. Iwamoto                            *
*                                                                      *
************************************************************************
      use moddas_material

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'ggsparam.inc'

      common /regdc/ idrg(kvlmax), idgr(kvmmax)

*-----------------------------------------------------------------------

      character dum1*10000
      character dum2*10000
      character chdf(500)*100
      dimension ildf(500)
      character cblan*100

      dimension kc(mc)
      dimension kr(mr)

      dimension vl(nr)
      dimension lr(nr)
      dimension ivl(nvl)
      dimension rvl(nvl)

*     for material
      common /inggs/  iog, igcel, ioa, igsuf, iob, igtrs
      common /regdm/  idmg(kvlmax)
      common /kmat1a/ mxmat, mxmat0, mxnel
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /kmat1i/ kmatc(kvlmax)

      common /talldc/ itdci(nientry,itlmax), rtdci(nrentry,itlmax),
     &               tbina(itlmax,maxitm), tbinb(itlmax,maxitm),
     &               beampw(itlmax,maxitm),
     &               tmina(itlmax,maxitm), tminb(itlmax,maxitm),
     &               aclst(itlmax,maxitm), bqlst(itlmax,maxitm)

      character aclst*8
      character chau*8
      character elmnt(104)*3

      data elmnt /
     &    ' H ','He ','Li ','Be ',' B ',' C ',' N ',' O ',' F ','Ne ',
     &    'Na ','Mg ','Al ','Si ',' P ',' S ','Cl ','Ar ',' K ','Ca ',
     &    'Sc ','Ti ',' V ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ',' Y ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ',' I ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ',' W ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ',' U ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

*-----------------------------------------------------------------------
*     read from kc(mc)
*-----------------------------------------------------------------------
            write(6,*) 'echrgdc2 before echrg2'
            call echrg2(mc,kc,dum1,lng1,icmb,igm)
            write(6,*) 'echrgdc2 after echrg2'
*-----------------------------------------------------------------------
*     input echo for tally region mesh
*-----------------------------------------------------------------------

                  cblan = ' '

                  ild1 = 12
                  ild0 = 72 - ild1

                  ilrm = lng1
                  isqd = 0
                  isrm = 0

  430             isqd = isqd + 1

                  if( ilrm .le. ild0 ) then

                        ildf(isqd) = ilrm
                     do k = 1, ilrm
                        chdf(isqd)(k:k) = dum1(k+isrm:k+isrm)

                     end do
            write(6,*) 'echrgdc2 9265'
                  else

                     do k = ild0, 1, -1

                        if( dum1(k+isrm:k+isrm) .eq. ' ' ) goto 420

                     end do

  420                k2 = k

                        ildf(isqd) = k2 - 1

                     do k = 1, k2 - 1

                        chdf(isqd)(k:k) = dum1(k+isrm:k+isrm)

                     end do

                     isrm = isrm + k2
                     ilrm = ilrm - k2

                     goto 430

                  end if

               if( icc .eq. 1 ) then

                  write(iot,'("      reg = ",100a1)')
     &                    (chdf(1)(j:j),j=1,ildf(1))

               else if( icc .eq. 2 ) then

                  write(iot,'(" reginbox = ",100a1)')
     &                    (chdf(1)(j:j),j=1,ildf(1))

               end if
            write(6,*) 'echrgdc2 9302'
               if( isqd .gt. 1 ) then

                  do k = 2, isqd

                     write(iot,'(100a1)') (cblan(j:j),j=1,ild1),
     &                    (chdf(k)(j:j),j=1,ildf(k))

                  end do

               end if
            write(6,*) 'echrgdc2 9313'
*-----------------------------------------------------------------------
*        echo for combined regions
*-----------------------------------------------------------------------

*    &       junf .eq. 0 .and. icmb .eq. 0 .and. nvl .eq. 0 ) ) return

         if( iva .eq. 0 ) then

*    &                  "# combined, lattice or level structure ")')
            write(iot,'("   target =",i10,7x,
     &                  "# Card[10a] target nuclide ")') igcel

            write(iot,'("   non     reg      vol     #",
     &                  " reg definition")')

         else

            write(iot,'("   value  ",18x,
     &                  "# values for each region")')
            write(iot,'("   non     reg      val     #",
     &                  " reg definition")')

         end if

            ivl=0
            rvl=0
            write(6,*) 'echrgdc2 9339'
               call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)
            write(6,*) 'echrgdc2 9341'
*-----------------------------------------------------------------------
*     read from kr(mr)
*-----------------------------------------------------------------------

                  j = 0

         do 500 ir = 1, nr

               call echrg3(j,mr,kr,dum2,lng1,icmb,igm)

                  ild1 = 30
                  ild0 = max(10,ierl-ild1)

                  ilrm = lng1
                  isqd = 0
                  isrm = 0

  530             isqd = isqd + 1

                  if( ilrm .le. ild0 ) then

                        ildf(isqd) = ilrm

                     do k = 1, ilrm

                        chdf(isqd)(k:k) = dum2(k+isrm:k+isrm)

                     end do

                  else

                     do k = ild0, 1, -1

                        if( dum2(k+isrm:k+isrm) .eq. ' ' ) goto 520

                     end do

  520                k2 = k

                        ildf(isqd) = k2 - 1

                     do k = 1, k2 - 1

                        chdf(isqd)(k:k) = dum2(k+isrm:k+isrm)

                     end do

                     isrm = isrm + k2
                     ilrm = ilrm - k2

                     goto 530

                  end if

*-----------------------------
*     for input check

               iintrp = 0
               do irp = 1, itdci(24,m)
                  if( lr(ir) .eq. itgcel(itdc2reg(itdc)+irp) )then !FURUTA20200522
                     iintrp = 1
                     goto 90
                  end if
               end do
                     iintrp = 0

   90          if( iintrp .eq. 1 ) then
                  write(iot,'(i5,2x,i7,1p1e13.4," # ",100a1)')
     &                 ir, lr(ir), rtgvlm(itdc2reg(itdc)+irp), !FURUTA20200522
     &                 (chdf(1)(l:l),l=1,ildf(1))

               else
                  write(iot,'(i5,2x,i7,1p1e13.4," # ",100a1)')
     &                 ir, lr(ir), vl(ir), (chdf(1)(l:l),l=1,ildf(1))
               end if

               if( isqd .gt. 1 ) then

                  do k = 2, isqd

                     write(iot,'(100a1)') (cblan(l:l),l=1,ild1-2),
     &                    '#',' ',(chdf(k)(l:l),l=1,ildf(k))

                  end do

               end if

*-----------------------------------------------------------------------
*     for tg-list(material) in [T-Dchain]
*-----------------------------------------------------------------------
            if( iintrp .eq. 1 ) then

             i0=itdc2reg(itnm2tdc(m)) !FURUTA20200522
             j0=ireg2itm(i0+irp)       !FURUTA20200522

             write(iot,'("   tg-list =",i7)') itgnum(i0+irp)
             do jrp = 1, itgnum(i0+irp)
              write(iot,'(8x,a8,6x,1p1e13.7)')
     &             tglst(j0+jrp),tgnlt(j0+jrp) !FURUTA20200522
             end do

            else
*-----------------------------
*     get the cell No. -> material No.

                  jcel = 0
               do jcel = 1, igcel
                  if( lr(ir) .eq. idrg(jcel) ) goto 100
               end do

  100       continue

*-----------------------------
*     get the material No. -> material list

                  kmat = 0
               do kmat = 1, mxmat0
                  if( idmg(jcel) .eq. idmn(kmat) ) goto 110
                  if( idmn(idnm(idmg(jcel))) .eq. idmn(kmat) ) goto 110
               end do

  110       continue

               nel  = nint( das_kmatc(kmatc(kmat)+1) )
               denh =       das_kmatc(kmatc(kmat)+2)
               nelh = nel

               if( denh .ne. 0.0d0 ) then
                   nelh = nelh + 1
               end if

               write(iot,'("   tg-list =",i7)') nelh

               if( denh .ne. 0.0d0 ) then
                   chau = ' H-1'
                   write(iot,'(8x,a4,8x,1p1e15.7)')
     &             chau(1:4), denh
               end if

               do kk = 1, nel
                  iz   = nint( das_kmatc(kmatc(kmat)+(kk-1)*4+34) )
                  ia   = nint( das_kmatc(kmatc(kmat)+(kk-1)*4+35) )
                  den  =       das_kmatc(kmatc(kmat)+(kk-1)*4+36)
                  chau = elmnt(iz)

                  if( ia .eq. 0 ) then
                     write(iot,'(8x,a2,10x,1p1e15.7)')
     &               chau(1:2), den
                  else if( ia .lt. 10 ) then
                     write(iot,'(8x,a2,"-",i1,8x,1p1e15.7)')
     &               chau(1:2), ia, den
                  else if( ia .lt. 100 ) then
                     write(iot,'(8x,a2,"-",i2,7x,1p1e15.7)')
     &               chau(1:2), ia, den
                  else
                     write(iot,'(8x,a2,"-",i3,6x,1p1e15.7)')
     &               chau(1:2), ia, den
                  end if
               end do

            end if
*-----------------------------------------------------------------------

  500    continue

*-----------------------------------------------------------------------

      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine getndc(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                  icrd,inucln,klst,hancza,ancza)
*                                                                      *
*       read pick up nuclide information                               *
*       made by N.Matsuda on 2013/07/17                                *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

      character hl(maxitm)*8,   hlc*8,   hle(maxitm)*8
      character hll(maxitm)*20, hllc*20
      dimension hlle(maxitm)

      character hancza(maxitm)*8
      dimension ancza(maxitm)



      dimension zaid(maxitm)

      integer icrd, inucln, klst
      integer izam, iz, ia, im, j1, j2
      character chaia1*1, chaia2*2, chaia3*3, chaim*1

*-----------------------------------------------------------------------
      character element(104)*3,cnuc*3,elmnt(104)*3

      data element /
     &    'h  ','he ','li ','be ','b  ','c  ','n  ','o  ','f  ','ne ',
     &    'na ','mg ','al ','si ','p  ','s  ','cl ','ar ','k  ','ca ',
     &    'sc ','ti ','v  ','cr ','mn ','fe ','co ','ni ','cu ','zn ',
     &    'ga ','ge ','as ','se ','br ','kr ','rb ','sr ','y  ','zr ',
     &    'nb ','mo ','tc ','ru ','rh ','pd ','ag ','cd ','in ','sn ',
     &    'sb ','te ','i  ','xe ','cs ','ba ','la ','ce ','pr ','nd ',
     &    'pm ','sm ','eu ','gd ','tb ','dy ','ho ','er ','tm ','yb ',
     &    'lu ','hf ','ta ','w  ','re ','os ','ir ','pt ','au ','hg ',
     &    'tl ','pb ','bi ','po ','at ','rn ','fr ','ra ','ac ','th ',
     &    'pa ','u  ','np ','pu ','am ','cm ','bk ','cf ','es ','fm ',
     &    'md ','no ','lr ','ku '/

      data elmnt /
     &    ' H ','He ','Li ','Be ',' B ',' C ',' N ',' O ',' F ','Ne ',
     &    'Na ','Mg ','Al ','Si ',' P ',' S ','Cl ','Ar ',' K ','Ca ',
     &    'Sc ','Ti ',' V ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ',' Y ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ',' I ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ',' W ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ',' U ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

*-----------------------------------------------------------------------
           iskip = 0
            ierr = 0
            iot  = 6

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140    continue

           call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return

**                                    ^^^ 954 or 955
               if( icrd .eq. 1 .and. jpn .eq. 3 ) goto 954
               if( icrd .eq. 2 .and. jpn .eq. 3 ) goto 955

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of the section
*-----------------------------------------------------------------------

            if( chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 800

            end if

*-----------------------------------------------------------------------
*        read nuclide information
*-----------------------------------------------------------------------

           ic = i1

              j = 0
           do i = 1, inucln * 2

             if( ic .gt. i3 ) then

  210          call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                 if( ierr .ne. 0 ) return
**                                      ^^^ 954 or 955
                 if( jpn .eq. 3 .and. icrd .eq. 1 ) goto 954
                 if( jpn .eq. 3 .and. icrd .eq. 2 ) goto 955

                 if( iskip .ne. 0 ) goto 210

                 ic = i1

             end if


             if( mod(i,2) .eq. 1 ) then
                 j = j + 1
                 hl(j)(1:8) = '        '
                 icl = inumc(chlw,ic,i3,' ') - 1
                 hl(j)(1:icl-ic+1) = chlw(ic:icl)
                 ic = icl + 1
                 ic = jnumc(chlw,ic,i3)
                 if( icrd .eq. 1 ) then
                   if( (inucln * 2 - 1) .eq. i ) then
                     hll(j)(1:20) = '                    '
                       goto 800
                   end if
                 end if
             else
                 hll(j)(1:20) = '                    '
                 if( icrd .eq. 2 ) then
                   icl = inumc(chlw,ic,i3,' ') - 1
                   hll(j)(1:icl-ic+1) = chlw(ic:icl)
                   ic = icl + 1
                   ic = jnumc(chlw,ic,i3)
                 end if
             end if


           end do

*-----------------------------------------------------------------------
*     summary and check
*-----------------------------------------------------------------------
  800    continue
* conversion

           do k = 1, inucln
               izam = 0
               iz = 0
               ia = 0
               im = 0
               j1 = 1
             hlc = hl(k)(1:8)

             if( hlc(1:1) .ge. 'a' .and. hlc(1:1) .le. 'z' ) then
*            write(iot,'(4x,"nucl name =  ",a)') hlc(j)(1:8)
               if( hlc(2:2) .eq. '-' ) hlc(2:2) = ' '

*     search the element
               do iz = 1, 104
                  if ( hlc(1:2) .eq. element(iz)(1:2) ) goto 810
               end do
               write(iot,'(" unknown element name ",a,
     &                     " was used.")') hlc(1:2)

  810    continue
               if( hl(k)(2:2) .eq. ' ' ) then
                    ia = 0
                    im = 0
                    hle(k) = elmnt(iz)(1:2)
               else if( hl(k)(2:2) .eq. '-' ) then
                    j1 = 3
               else if( hl(k)(3:3) .eq. ' ' ) then
                    ia = 0
                    im = 0
                    hle(k) = elmnt(iz)(1:2)
               else if( hl(k)(3:3) .eq. '-' ) then
                    j1 = 4
               else
                    goto 954
               end if

               if( j1 .eq. 3 .or. j1 .eq. 4 ) then
                 do kk = j1, 8
                    if( hlc(kk:kk) .eq. ' ' ) then
                           j2 = kk - 1
                      if( hlc(j2:j2) .eq. 'g' ) then
                        call onum(hlc,j1,j2-1,cvvv,ierr)
                             ia = nint( cvvv )
                             im = 0
                        if( ia .gt. 99 ) then
                          write(chaia3,'(i3.3)') ia
                          hle(k) = elmnt(iz)(1:2)//'-'//chaia3//'g'
                        else if( ia .gt. 9 .and. ia .lt. 100 ) then
                          write(chaia2,'(i2.2)') ia
                          hle(k) = elmnt(iz)(1:2)//'-'//chaia2//'g'
                        else if( ia .lt. 10 ) then
                          write(chaia1,'(i1.1)') ia
                          hle(k) = elmnt(iz)(1:2)//'-'//chaia1//'g'
                        end if
                          goto 850
                      else if( hlc(j2:j2) .eq. 'm' ) then
                        call onum(hlc,j1,j2-1,cvvv,ierr)
                             ia = nint( cvvv )
                             im = 1
                        if( ia .gt. 99 ) then
                          write(chaia3,'(i3.3)') ia
                          hle(k) = elmnt(iz)(1:2)//'-'//chaia3//'m'
                        else if( ia .gt. 9 .and. ia .lt. 100 ) then
                          write(chaia2,'(i2.2)') ia
                          hle(k) = elmnt(iz)(1:2)//'-'//chaia2//'m'
                        else if( ia .lt. 10 ) then
                          write(chaia1,'(i1.1)') ia
                          hle(k) = elmnt(iz)(1:2)//'-'//chaia1//'m'
                        end if
                          goto 850
                      else if( hlc(j2:j2) .eq. 'n' ) then
                        call onum(hlc,j1,j2-1,cvvv,ierr)
                             ia = nint( cvvv )
                             im = 2
                        if( ia .gt. 99 ) then
                          write(chaia3,'(i3.3)') ia
                          hle(k) = elmnt(iz)(1:2)//'-'//chaia3//'n'
                        else if( ia .gt. 9 .and. ia .lt. 100 ) then
                          write(chaia2,'(i2.2)') ia
                          hle(k) = elmnt(iz)(1:2)//'-'//chaia2//'n'
                        else if( ia .lt. 10 ) then
                          write(chaia1,'(i1.1)') ia
                          hle(k) = elmnt(iz)(1:2)//'-'//chaia1//'n'
                        end if
                          goto 850
                      else
                        call onum(hlc,j1,j2,cvvv,ierr)
                             ia = nint( cvvv )
                             im = 0
                        if( ia .gt. 99 ) then
                          write(chaia3,'(i3.3)') ia
                          hle(k) = elmnt(iz)(1:2)//'-'//chaia3
                        else if( ia .gt. 9 .and. ia .lt. 100 ) then
                          write(chaia2,'(i2.2)') ia
                          hle(k) = elmnt(iz)(1:2)//'-'//chaia2
                        else if( ia .lt. 10 ) then
                          write(chaia1,'(i1.1)') ia
                          hle(k) = elmnt(iz)(1:2)//'-'//chaia1
                        end if
                          goto 850
                      end if

                    end if

                 end do

  850    continue
               end if

               zaid(k) = 1000 * iz + ia + 0.1 * im

             else
               do kk = 1, 8
                  if( hlc(kk:kk) .eq. ' ' ) goto 860
               end do

  860    continue
               call onum(hlc,1,kk-1,cvvv,ierr)
                 zaid(k) = cvvv
                 izam = zaid(k) * 10

                 iz = int( izam / 10000 )
                 ia = int( ( izam - iz * 10000 ) / 10 )
                 im = izam - ( iz * 10000 ) - ( ia * 10 )

                    chaim = ' '
                 if( im .eq. 1 ) chaim = 'm'
                 if( im .eq. 2 ) chaim = 'n'

                 if( ia .gt. 99 ) then
                   write(chaia3,'(i3.3)') ia
                   hle(k) = elmnt(iz)(1:2)//'-'//chaia3//chaim
                 else if( ia .gt. 9 .and. ia .lt. 100 ) then
                   write(chaia2,'(i2.2)') ia
                   hle(k) = elmnt(iz)(1:2)//'-'//chaia2//chaim
                 else if( ia .lt. 10 ) then
                   write(chaia1,'(i1.1)') ia
                   hle(k) = elmnt(iz)(1:2)//'-'//chaia1//chaim
                 end if

             end if

           end do

           if( icrd .eq. 2 ) then

             do k = 1, inucln

               hllc = hll(k)(1:20)

               do kk = 1, 20
                 if( hllc(kk:kk) .eq. ' ' ) then
                     j2 = kk - 1
                     goto 870
                  end if
               end do

  870    continue
               if( hllc(1:1) .eq. 'B' .or. hllc(1:1) .eq. 'b' ) then
                 call onum(hllc,2,j2,cvvv,ierr)
                 hlle(k) = cvvv

               else
                 call onum(hllc,1,j2,cvvv,ierr)
                 hlle(k) = cvvv

               end if

             end do

           end if

* summary
           k0=iregitm                                      !FURUTA20200522
           iregitm=iregitm+inucln                          !FURUTA20200522
           if(iregitm.gt.maxregitm)call reallocate_talldc3(iregitm) !FURUTA20200522
           do kkk = 1, inucln

             if( icrd .eq. 1 ) then
               hancza(kkk) =   hle(kkk)
               ancza(kkk)  =  zaid(kkk)
              else
               htgnzas(k0+kkk)=hle(kkk)  !FURUTA20200522
               tgnzas(k0+kkk)=zaid(kkk)  !FURUTA20200522
               ctgnnds(k0+kkk)=hlle(kkk) !FURUTA20200522
             end if

           end do

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  954 continue

         m_err = 'description or number of Card[8] is wrong.'
         ErrCha = ''
         ErrID = 'L:11868/R:getndc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  955 continue

         m_err = 'description or number of Card[10] is wrong.'
         ErrCha = ''
         ErrID = 'L:11878/R:getndc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end subroutine


************************************************************************
*                                                                      *
      subroutine trzmeshdc(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     rzx0,rzy0,
     &                     irtp,inr,rmin,rmax,rdel,istrg,
     &                     iztp,inz,zmin,zmax,zdel,istzg,
     &                     idcrg)
*                          ^^^^^ for [T-DCHAIN]
*                                                                      *
*       read r-z mesh sub-section of input tally section               *
*       modified by T.Miura on 2017/11/30                              *
*       original sub. is "trzmesh" in tallsm1.f                        *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

CCSE add for mesh=r-z parameter (2017.11.30) >>>>>
c-----------------------------------------------------------------------
      integer   idcrg
CCSE add for mesh=r-z parameter (2017.11.30) <<<<<
*-----------------------------------------------------------------------
*     initial values
*-----------------------------------------------------------------------

            ierr  = 0

            rzx0 = 0.0
            rzy0 = 0.0

            irtp = -1
            iztp = -1

CCSE add for mesh=r-z parameter (2017.11.30) >>>>>
            idcrg = 0
CCSE add for mesh=r-z parameter (2017.11.30) <<<<<

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

CCSE add for mesh=r-z parameter (2017.11.30) >>>>>
               idcrg = idcrg + 1
               if(idcrg.gt.maxreg)call reallocate_talldc4(idcrg) !FURUTA20200522
               dcrgm(idcrg) = chin
CCSE add for mesh=r-z parameter (2017.11.30) <<<<<

  150 continue

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

*-----------------------------------------------------------------------

         if( chcm(i1:i1+2) .eq. 'x0=' ) then

               ic = inumc(chlw,i1,i3,'=') + 1

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               rzx0 = cvvv

               goto 140

         else if( chcm(i1:i1+2) .eq. 'y0=' ) then

               ic = inumc(chlw,i1,i3,'=') + 1

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               rzy0 = cvvv

               goto 140

         else if( chcm(i1:i1+6) .eq. 'r-type=' ) then

               ic = inumc(chlw,i1,i3,'=') + 1

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               irtp = nint( cvvv )

CCSE change for mesh=r-z parameter (2017.11.30) >>>>>
               call getmshdc(jsn,jsi,dsin,idsi,ill,ilf,
     &                       jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                       'r',irtp,inr,rmin,rmax,rdel,istrg,
     &                       idcrg)
CCSE change for mesh=r-z parameter (2017.11.30) <<<<<

               if( ierr .ne. 0 ) return

               goto 150

         else if( chcm(i1:i1+6) .eq. 'z-type=' ) then

               ic = inumc(chlw,i1,i3,'=') + 1

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               iztp = nint( cvvv )

CCSE change for mesh=r-z parameter (2017.11.30) >>>>>
               call getmshdc(jsn,jsi,dsin,idsi,ill,ilf,
     &                       jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                       'z',iztp,inz,zmin,zmax,zdel,istzg,
     &                       idcrg)
CCSE change for mesh=r-z parameter (2017.11.30) <<<<<

               if( ierr .ne. 0 ) return

               goto 150

         else

               if( irtp .eq. -1 ) then

                  m_err = 'r-type is missing in r-z mesh'
                  ErrCha = ''
                  ErrID = 'L:12043/R:trzmeshdc/F:tallsm4.f'
                  goto 998

               else if( iztp .eq. -1 ) then

                  m_err = 'z-type is missing in r-z mesh'
                  ErrCha = ''
                  ErrID = 'L:12050/R:trzmeshdc/F:tallsm4.f'
                  goto 998

               end if

            return

         end if

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of r-z mesh sub-section is wrong.'
         ErrCha = ''
         ErrID = 'L:12067/R:trzmeshdc/F:tallsm4.f'

  998 continue

         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end subroutine trzmeshdc

************************************************************************
*                                                                      *
      subroutine txymeshdc(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     ixtp,inx,xmin,xmax,xdel,istxg,
     &                     iytp,iny,ymin,ymax,ydel,istyg,
     &                     iztp,inz,zmin,zmax,zdel,istzg,
     &                     idcrg)
*                          ^^^^^ for [T-DCHAIN]
*                                                                      *
*       read xyz mesh sub-section of input tally section               *
*       modified by T.Miura on 2017/11/30                              *
*       original sub. is "txymesh" in tallsm1.f                        *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

CCSE add for mesh=xyz parameter (2017.11.30) >>>>>
c-----------------------------------------------------------------------
      integer   idcrg
CCSE add for mesh=xyz parameter (2017.11.30) <<<<<
*-----------------------------------------------------------------------
*     initial values
*-----------------------------------------------------------------------

            ierr  = 0

            ixtp = -1
            iytp = -1
            iztp = -1

CCSE add for mesh=xyz parameter (2017.11.30) >>>>>
            idcrg = 0
CCSE add for mesh=xyz parameter (2017.11.30) <<<<<

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

CCSE add for mesh=xyz parameter (2017.11.30) >>>>>
               idcrg = idcrg + 1
               if(idcrg.gt.maxreg)call reallocate_talldc4(idcrg) !FURUTA20200522
               dcrgm(idcrg) = chin
CCSE add for mesh=xyz parameter (2017.11.30) <<<<<

  150 continue

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

*-----------------------------------------------------------------------

         if( chcm(i1:i1+6) .eq. 'x-type=' ) then

               ic = inumc(chlw,i1,i3,'=') + 1

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               ixtp = nint( cvvv )

CCSE change for mesh=xyz parameter (2017.11.30) >>>>>
               call getmshdc(jsn,jsi,dsin,idsi,ill,ilf,
     &                       jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                       'x',ixtp,inx,xmin,xmax,xdel,istxg,
     &                       idcrg)
CCSE change for mesh=xyz parameter (2017.11.30) <<<<<

               if( ierr .ne. 0 ) return

               goto 150

         else if( chcm(i1:i1+6) .eq. 'y-type=' ) then

               ic = inumc(chlw,i1,i3,'=') + 1

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               iytp = nint( cvvv )

CCSE change for mesh=xyz parameter (2017.11.30) >>>>>
               call getmshdc(jsn,jsi,dsin,idsi,ill,ilf,
     &                       jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                       'y',iytp,iny,ymin,ymax,ydel,istyg,
     &                       idcrg)
CCSE change for mesh=xyz parameter (2017.11.30) <<<<<

               if( ierr .ne. 0 ) return

               goto 150

         else if( chcm(i1:i1+6) .eq. 'z-type=' ) then

               ic = inumc(chlw,i1,i3,'=') + 1

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               iztp = nint( cvvv )

CCSE change for mesh=xyz parameter (2017.11.30) >>>>>
               call getmshdc(jsn,jsi,dsin,idsi,ill,ilf,
     &                       jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                       'z',iztp,inz,zmin,zmax,zdel,istzg,
     &                       idcrg)
CCSE change for mesh=xyz parameter (2017.11.30) <<<<<

               if( ierr .ne. 0 ) return

               goto 150

         else

               if( ixtp .eq. -1 ) then

                  m_err = 'x-type is missing in xyz mesh'
                  ErrCha = ''
                  ErrID = 'L:12229/R:txymeshdc/F:tallsm4.f'
                  goto 998

               else if( iytp .eq. -1 ) then

                  m_err = 'y-type is missing in xyz mesh'
                  ErrCha = ''
                  ErrID = 'L:12236/R:txymeshdc/F:tallsm4.f'
                  goto 998

               else if( iztp .eq. -1 ) then

                  m_err = 'z-type is missing in xyz mesh'
                  ErrCha = ''
                  ErrID = 'L:12243/R:txymeshdc/F:tallsm4.f'
                  goto 998

               end if

            return

         end if

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of xyz mesh sub-section is wrong.'
         ErrCha = ''
         ErrID = 'L:12260/R:txymeshdc/F:tallsm4.f'

  998 continue

         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end subroutine txymeshdc

************************************************************************
*                                                                      *
      subroutine getmshdc(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                    gdum,igtp,igr,gmin,gmax,gdel,igm,
     &                    idcrg)
*                         ^^^^^  for [T-DCHAIN]
*                                                                      *
*       get mesh information                                           *
*       modified by T.Miura on 2017/11/30                              *
*       original sub. is "getmsh" in tallsm1.f                         *
*                                                                      *
************************************************************************
      use moddas
      use moddas_mesh

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------


      character gdum*1

      dimension lschn(10), ischn(10)
      character schan(10)*8


*-----------------------------------------------------------------------

      data icsu / 4 /

      data ( schan(i), i = 1, 4 ) /
     &    'n       ',' min    ',' max    ',' del    '/

      data ( lschn(i), i = 1, 4 ) /
     &     2,        4,          4,         4/

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

CCSE add for mesh=r-z/xyz parameter (2017.11.30) >>>>>
c-----------------------------------------------------------------------
      integer   idcrg
CCSE add for mesh=r-z/xyz parameter (2017.11.30) <<<<<
*-----------------------------------------------------------------------

      ierrMSG=0 ! T.Sato 2017/07/14

         do i = 1, icsu

            ischn(i) = 0

         end do

*-----------------------------------------------------------------------
*     initial memory point
*-----------------------------------------------------------------------

         igm = iaddress_gmsh(icurrent_gmsh)

*-----------------------------------------------------------------------
*     check mesh type
*-----------------------------------------------------------------------

         if( igtp .lt. 1 .or. igtp .gt. 5 ) goto 999

*-----------------------------------------------------------------------
*     mesh type
*-----------------------------------------------------------------------

            schan(1)(2:2) = gdum

         do i = 2, icsu

            schan(i)(1:1) = gdum

         end do

*-----------------------------------------------------------------------

            ierr  = 0

            icmp = 0
            irct = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 800

               if( iskip .ne. 0 ) goto 140

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
CCSE add for mesh=r-z/xyz parameter (2017.11.30) >>>>>
            idcrg = idcrg + 1
            if(idcrg.gt.maxreg)call reallocate_talldc4(idcrg) !FURUTA20200522
            dcrgm(idcrg) = chin
CCSE add for mesh=r-z/xyz parameter (2017.11.30) <<<<<

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

               ischn( ipm ) = 1

               ic = inumc(chlw,il+1,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 998

               icl = inumc(chlw,ic,i3,';') - 1

*-----------------------------------------------------------------------
*        mesh informations
*-----------------------------------------------------------------------

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 998

         if( ipm .eq. 1 ) then

               igrr = nint( cvvv )
               igr  = iabs( igrr )

               if( igr .le. 0 ) goto 997

               call moddas_reallocate_dbl(
     &                 3, icurrent_gmsh, igr+1, iaddress_gmsh, gmsh)

               if( igtp .eq. 1 ) goto 500

         else if( ipm .eq. 2 ) then

               gmin = cvvv

         else if( ipm .eq. 3 ) then

               gmax = cvvv

         else if( ipm .eq. 4 ) then

               gdel = cvvv

         end if

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

               goto 140

*-----------------------------------------------------------------------
*        read mesh points
*-----------------------------------------------------------------------

  500    continue

                  if( chlw(icl+1:icl+1) .eq. ';' ) goto 996

  148             call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                  if( ierr .ne. 0 ) return
                  if( jpn  .eq. 3 ) goto 997

                  if( iskip .ne. 0 ) goto 148

CCSE add for mesh=r-z/xyz parameter (2017.11.30) >>>>>
                  idcrg = idcrg + 1
               if(idcrg.gt.maxreg)call reallocate_talldc4(idcrg) !FURUTA20200522
                  dcrgm(idcrg) = chin
CCSE add for mesh=r-z/xyz parameter (2017.11.30) <<<<<

                  ic = i1

            do i = 1, igr + 1

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

                  gmsh(igm-1+i) = cvvv

               ic = ic2

               if( i .le. igr .and. ic .gt. i3 ) then

  149             call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                  if( ierr .ne. 0 ) return
                  if( jpn  .eq. 3 ) goto 997

                  if( iskip .ne. 0 ) goto 149

CCSE add for mesh=r-z/xyz parameter (2017.11.30) >>>>>
                  idcrg = idcrg + 1
               if(idcrg.gt.maxreg)call reallocate_talldc4(idcrg) !FURUTA20200522
                  dcrgm(idcrg) = chin
CCSE add for mesh=r-z/xyz parameter (2017.11.30) <<<<<

                  ic = i1

               end if

            end do

                  icmp = 1

               goto 140

*-----------------------------------------------------------------------
*     summary and check
*-----------------------------------------------------------------------

  800 continue

            if( igtp .eq. 1 ) then

               if( ischn(1) .eq. 0 .or.
     &             igr  .le. 0     .or.
     &             icmp .eq. 0 ) then
                ierrMSG=1 ! T.Sato 2017/07/14
                goto 998
               endif

               do i = 1, igr

                     gmins = gmsh(igm-1+i)

                  do j = i + 1, igr + 1

                     if( gmsh(igm-1+j) .lt. gmins ) then

                        gmins         = gmsh(igm-1+j)
                        gmsh(igm-1+j) = gmsh(igm-1+i)
                        gmsh(igm-1+i) = gmins

                     end if

                  end do

               end do

            else if( igtp .eq. 2 ) then

               if( ischn(1) .eq. 0 .or.
     &             ischn(2) .eq. 0 .or.
     &             ischn(3) .eq. 0 .or.
     &             gmax .le. gmin  .or.
     &             igr  .le. 0 ) then
                  ierrMSG=2 ! T.Sato 2017/07/14
                  goto 998
               endif

                  gdel = ( gmax - gmin ) / dble( igr )

               do i = 1, igr

                  gmsh(igm-1+i) = gmin + gdel * dble( i - 1 )

               end do

                  gmsh(igm-1+igr+1) = gmax

            else if( igtp .eq. 3 ) then

               if( ischn(1) .eq. 0 .or.
     &             ischn(2) .eq. 0 .or.
     &             ischn(3) .eq. 0 .or.
     &             gmin .le. 0.0   .or.
     &             gmax .le. gmin  .or.
     &             igr  .le. 0 ) then
                ierrMSG=3 ! T.Sato 2017/07/14
                goto 998
               endif

                  gdel = log( gmax / gmin ) / dble( igr )

               do i = 1, igr

                  gmsh(igm-1+i) = gmin * exp( gdel * dble( i - 1 ) )

               end do

                  gmsh(igm-1+igr+1) = gmax

            else if( igtp .eq. 4 ) then

               if( ischn(4) .eq. 0 .or.
     &             ischn(2) .eq. 0 .or.
     &             ischn(3) .eq. 0 .or.
     &             gmax .le. gmin  .or.
     &             gdel .le. 0.0 ) then
                ierrMSG=4 ! T.Sato 2017/07/14
                goto 998
               endif

                  igr = int( ( gmax - gmin ) / gdel ) + 1

               call moddas_reallocate_dbl(
     &                 3, icurrent_gmsh, igr+1, iaddress_gmsh, gmsh)

               do i = 1, igr + 1

                  gmsh(igm-1+i) = gmin + gdel * dble( i - 1 )

               end do

                  if( gmsh(igm-1+igr) .eq. gmax ) igr = igr - 1

                  igrr = igr

            else if( igtp .eq. 5 ) then

               if( ischn(4) .eq. 0 .or.
     &             ischn(2) .eq. 0 .or.
     &             ischn(3) .eq. 0 .or.
     &             gmin .le. 0.0   .or.
     &             gmax .le. gmin  .or.
     &             gdel .le. 0.0 ) then
                ierrMSG=5 ! T.Sato 2017/07/14
                goto 998
               endif

                  igr = int( log( gmax / gmin ) / gdel ) + 1

               call moddas_reallocate_dbl(
     &                 3, icurrent_gmsh, igr+1, iaddress_gmsh, gmsh)

               do i = 1, igr + 1

                  gmsh(igm-1+i) = gmin * exp( gdel * dble( i - 1 ) )

               end do

                  if( gmsh(igm-1+igr) .eq. gmax ) igr = igr - 1

                  igrr = igr

            end if

                  iaddress_gmsh(2) = iaddress_gmsh(3)
                  icurrent_gmsh = 2

                  if( mmmax .ge. mdas ) goto 995

               igr = igrr

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  995 continue

         m_err = 'Memory error: mmmax exceeds mdas '//
     &           ': Please extend mdas in param.inc'
         ErrCha = ''
         ErrID = 'L:12690/R:getmshdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'In this line, [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:12702/R:getmshdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of mesh points is wrong'
         ErrCha = ''
         ErrID = 'L:12714/R:getmshdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue
      if(ierrMSG.le.1.and.ischn(1).eq.0) then ! T.Sato 2017/07/14
       m_err='Number of mesh points should be defined for *-type =< 3'
       ErrCha = ''
       ErrID = 'L:12726/R:getmshdc/F:tallsm4.f'
      elseif(ierrMSG.eq.1) then
       m_err='Each mesh point must be defined for *-type = 1'
       ErrCha = ''
       ErrID = 'L:12730/R:getmshdc/F:tallsm4.f'
      elseif(ierrMSG.ge.2.and.(ischn(2).eq.0.or.ischn(3).eq.0)) then
       m_err='Both *min= & *max= should be defined for *-type >= 2'
       ErrCha = ''
       ErrID = 'L:12734/R:getmshdc/F:tallsm4.f'
      elseif(ierrMSG.ge.4.and.ischn(4).eq.0) then
       m_err='*del= should be defined for *-type = 4 or 5'
       ErrCha = ''
       ErrID = 'L:12738/R:getmshdc/F:tallsm4.f'
      elseif((ierrMSG.eq.3.or.ierrMSG.eq.5).and.gmin.le.0.0) then
       m_err='Minimum value should be positive for *-type = 3 or 5'
       ErrCha = ''
       ErrID = 'L:12742/R:getmshdc/F:tallsm4.f'
      else
       m_err = 'Description of mesh parameter is wrong'
       ErrCha = ''
       ErrID = 'L:12746/R:getmshdc/F:tallsm4.f'
      endif
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Mesh type should be 1 - 5'
         ErrCha = ''
         ErrID = 'L:12759/R:getmshdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end subroutine getmshdc


CCSE add new subroutine for mesh=r-z/xyz (2018.07.31) >>>>>
************************************************************************
*                                                                      *
      subroutine echxydc(iot,nx,ny,nz,idx,idy,idz,mtscore,m)
*                                                                      *
*       input echo for tally xyz mesh                                  *
*       modified by S.Abe on 2018/07/24                                *
*                                                                      *
************************************************************************
      use moddas_material
      use moddas_mesh

      implicit real*8 (a-h,o-z)

      include 'param.inc'
C MATSUDA 2017.10.04 (expand natural nucleus)
      include 'param01.inc'  ! maxpt

*-----------------------------------------------------------------------

      common /kmat1a/ mxmat, mxmat0, mxnel
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /kmat1i/ kmatc(kvlmax)
C MATSUDA 2017.10.04 (expand natural nucleus)
      common /natnuc/ natnn(maxpt), natnm(maxpt,10), patnn(maxpt,10)

      common /redufmt/ iredufmt(itlmax) !FURUTA20200601
*-----------------------------------------------------------------------

      character chau*8
      character elmnt(104)*3

      data elmnt /
     &    ' H ','He ','Li ','Be ',' B ',' C ',' N ',' O ',' F ','Ne ',
     &    'Na ','Mg ','Al ','Si ',' P ',' S ','Cl ','Ar ',' K ','Ca ',
     &    'Sc ','Ti ',' V ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ',' Y ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ',' I ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ',' W ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ',' U ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

*-----------------------------------------------------------------------

      real(8),allocatable :: vol(:,:,:,:), sumvol(:)
      integer(8),allocatable :: matlst(:)
      character(13),allocatable :: chvol(:)

*-----------------------------------------------------------------------
*     initialize
*-----------------------------------------------------------------------

      allocate( vol(nx,ny,nz,0:mxmat),sumvol(mxmat) )
      allocate( matlst(mxmat) )
      allocate( chvol(mxmat) )

      vol(:,:,:,:) = 0.d0
      sumvol(:) = 0.d0

      nmesh = nx * ny * nz

*-----------------------------------------------------------------------
*     calculate volumes for each xyz grid and material
*     - nsm might be better to change input parameter
*-----------------------------------------------------------------------

      nl = 1
      ilt = 1
      nsm = 18
      itrns = 0

      do ix = 1, nx

       xxs0 = das_itxrg(idx+ix-1)
       xxf0 = das_itxrg(idx+ix)

       do iy = 1, ny

        yys0 = das_ityrg(idy+iy-1)
        yyf0 = das_ityrg(idy+iy)

        do iz = 1, nz

         zzs0 = das_itzrg(idz+iz-1)
         zzf0 = das_itzrg(idz+iz)

         vol(ix,iy,iz,0) = vls(1,0,ilt,0,itrns,
     &                         xxs0,xxf0,yys0,yyf0,zzs0,zzf0)

         do imat = 1, mxmat
            vol(ix,iy,iz,imat) = vlsdc(1,idmn(imat),ilt,-nsm,itrns,
     &                                 xxs0,xxf0,yys0,yyf0,zzs0,zzf0)
            sumvol(imat) = sumvol(imat) + vol(ix,iy,iz,imat)
         enddo

        enddo
       enddo
      enddo

*-----------------------------------------------------------------------
*     check whether materals exist in whole mesh
*-----------------------------------------------------------------------

      mtcnt = 0
      do imat = 1, mxmat
         if( sumvol(imat) .gt. 0.d0 ) then
            mtcnt = mtcnt + 1
            matlst(mtcnt) = imat
            write(chvol(mtcnt),'(a10,i3.3)') "    volmat",mtcnt
         endif
      enddo

*-----------------------------------------------------------------------
*     write material part
*     mtscore will be change to input parameter
*-----------------------------------------------------------------------


      write(iot,'("      imat =",i7)') mtcnt
      write(iot,'("   mtscore =",i7)') mtscore

      do imat = 1, mtcnt

         kmat = matlst(imat)
         nel  = nint( das_kmatc(kmatc(kmat)+1) )
         denh =       das_kmatc(kmatc(kmat)+2)
         nelh = nel

         if( denh .ne. 0.0d0 ) then
            nelh = nelh + 1
         end if

         lmat = nelh
         denstm = 0.0d0
         do kk = 1, nel
            iz   = nint( das_kmatc(kmatc(kmat)+(kk-1)*4+34) )
            ia   = nint( das_kmatc(kmatc(kmat)+(kk-1)*4+35) )
            den  =       das_kmatc(kmatc(kmat)+(kk-1)*4+36)

            if( ia .eq. 0 .and. natnn(iz) .ne. 0 ) then
               lmat = lmat + natnn(iz) - 1
            end if
         end do

         write(iot,'("   mt-list =",i7)') lmat

         if( denh .ne. 0.0d0 ) then
            chau = ' H-1'
            write(iot,'(8x,a4,8x,1p1e15.7)') chau(1:4), denh
         end if

         do kk = 1, nel
            iz   = nint( das_kmatc(kmatc(kmat)+(kk-1)*4+34) )
            ia   = nint( das_kmatc(kmatc(kmat)+(kk-1)*4+35) )
            den  =       das_kmatc(kmatc(kmat)+(kk-1)*4+36)
            chau = elmnt(iz)

            if( ia .eq. 0 ) then

               do kkk = 1, natnn(iz)

                  if( iz .eq. 1 .and. natnm(iz,kkk) .eq. 1 ) then

                     write(iot,'(8x,a2,"-1"8x,1p1e15.7)')
     &                chau(1:2), denh + den * patnn(iz,kkk) / 100.0

                  else

                     denstm = den * patnn(iz,kkk) / 100.0

                     if( natnm(iz,kkk) .lt. 10 ) then
                        write(iot,'(8x,a2,"-",i1,8x,1p1e15.7)')
     &                   chau(1:2), natnm(iz,kkk), denstm
                     else if( natnm(iz,kkk) .lt. 100 ) then
                        write(iot,'(8x,a2,"-",i2,7x,1p1e15.7)')
     &                   chau(1:2), natnm(iz,kkk), denstm
                     else
                        write(iot,'(8x,a2,"-",i3,6x,1p1e15.7)')
     &                   chau(1:2), natnm(iz,kkk), denstm
                     end if

                  end if

               end do

            else if( ia .lt. 10 ) then
               write(iot,'(8x,a2,"-",i1,8x,1p1e15.7)')
     &          chau(1:2), ia, den
            else if( ia .lt. 100 ) then
               write(iot,'(8x,a2,"-",i2,7x,1p1e15.7)')
     &          chau(1:2), ia, den
            else
               write(iot,'(8x,a2,"-",i3,6x,1p1e15.7)')
     &          chau(1:2), ia, den
            end if

         end do

      enddo

*-----------------------------------------------------------------------
*     write mesh part without fluxs ... fluxs will be added in tdctoinp
*-----------------------------------------------------------------------

      write(iot,'("     imesh = ",i7)') nmesh
      write(iot,'("   non    ix    iy    iz       volume",
     &            20a13)') (chvol(jmat),jmat=1,mtcnt)

      ndum = 0
      if(iredufmt(m).gt.0)then !FURUTA20200601
       do iz = 1, nz
        do iy = 1, ny
         do ix = 1, nx
          ndum = ndum + 1
          write(iot,'(4i6,1p20e13.4)')
     &         ndum, ix, iy, iz, vol(ix,iy,iz,0),
     &         (vol(ix,iy,iz,matlst(jmat)),jmat=1,mtcnt)
         enddo
        enddo
       enddo
      else
       do ix = 1, nx
        do iy = 1, ny
         do iz = 1, nz
          ndum = ndum + 1
          write(iot,'(4i6,1p20e13.4)')
     &         ndum, ix, iy, iz, vol(ix,iy,iz,0),
     &         (vol(ix,iy,iz,matlst(jmat)),jmat=1,mtcnt)
         enddo
        enddo
       enddo
      endif

*-----------------------------------------------------------------------

      deallocate( vol,sumvol )
      deallocate( matlst )
      deallocate( chvol )

*-----------------------------------------------------------------------

      return
      end subroutine echxydc
CCSE add new subroutine for mesh=r-z/xyz (2018.07.31) <<<<<


CCSE add new function for mesh=xyz (2018.07.31) >>>>>
************************************************************************
*                                                                      *
      function vlsdc(nl,lt,ilt,nsm,itrns,xxs0,xxf0,yys0,yyf0,zzs0,zzf0)
*                                                                      *
*       new function to calculate volmat_xxx
*       original code : function vls                                   *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      parameter( eps  = 1.0d-08 )
      parameter( epss = 1.0d-04 )

*-----------------------------------------------------------------------

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

*-----------------------------------------------------------------------

      dimension   lt(nl)

      dimension   xm(1000)
      dimension   ym(1000)
      dimension   zm(1000)

      dimension   xx(2)
      dimension   yy(2)
      dimension   zz(2)

      dimension   nmedn(8), isn(2), uvw(3,3)

      data isn / 1, -1 /

*-----------------------------------------------------------------------

            dx = xxf0 - xxs0
            dy = yyf0 - yys0
            dz = zzf0 - zzs0

            vls0 = dx * dy * dz

*-----------------------------------------------------------------------
*     no material selection
*-----------------------------------------------------------------------

         if( nl .eq. 0 .or. nsm .eq. 0 ) goto 1000

*-----------------------------------------------------------------------
*     transform
*-----------------------------------------------------------------------

         call trnsxv(xxs0,yys0,zzs0,xxs,yys,zzs,itrns)
         call trnsxv(xxf0,yyf0,zzf0,xxf,yyf,zzf,itrns)

         call trnsuv(1.d0,0.d0,0.d0,uvw(1,1),uvw(1,2),uvw(1,3),itrns)
         call trnsuv(0.d0,1.d0,0.d0,uvw(2,1),uvw(2,2),uvw(2,3),itrns)
         call trnsuv(0.d0,0.d0,1.d0,uvw(3,1),uvw(3,2),uvw(3,3),itrns)

*-----------------------------------------------------------------------

            nsn = abs( nsm )

            xx(1) = xxs
            xx(2) = xxf
            yy(1) = yys
            yy(2) = yyf
            zz(1) = zzs
            zz(2) = zzf

*-----------------------------------------------------------------------
*     check only of the eight corners only for nsm > 0
*-----------------------------------------------------------------------

      if( nsm .gt. 0 ) then

                  in = 0

            do ix = 1, 2
            do iy = 1, 2
            do iz = 1, 2

                  in = in + 1

                  ici   = -1
                  mark  =  1
                  markp =  0

                  xi = xx(ix) + uvw(1,1) * isn(ix) * eps
                  yi = yy(iy) + uvw(2,2) * isn(iy) * eps
                  zi = zz(iz) + uvw(3,3) * isn(iz) * eps

                  u = uvw(1,1) * isn(ix)
                  v = uvw(2,2) * isn(iy)
                  w = uvw(3,3) * isn(iz)

                  rnm = sqrt( u**2 + v**2 + w**2 )
                  u = u / rnm
                  v = v / rnm
                  w = w / rnm

                  call gomsor(xi,yi,zi,u,v,w,
     &                        nmedn(in),iblz,mark,markp,ici)

                  if( mark .lt. -1 ) goto 1000

            end do
            end do
            end do

                  ichk = 0

               do i = 1, 8

                     do k = 1, nl

                        if( ( ilt .gt. 0 .and.
     &                        idmn(nmedn(i)) .eq. lt(k) ) .or.
     &                      ( ilt .lt. 0 .and.
     &                        idmn(nmedn(i)) .ne. lt(k) ) ) then

                           ichk = ichk + 1

                           goto 2000

                        end if

                     end do

 2000                continue

               end do

               if( ichk .eq. 8 ) goto 1000

      end if

*-----------------------------------------------------------------------
*     with material selection
*-----------------------------------------------------------------------

               ddx = dx / nsn
               ddy = dy / nsn
               ddz = dz / nsn

            do i = 1, nsn

               xm(i) = xx(1) + ddx / 2.0 + dble(i-1) * ddx
               ym(i) = yy(1) + ddy / 2.0 + dble(i-1) * ddy
               zm(i) = zz(1) + ddz / 2.0 + dble(i-1) * ddz

            end do

               xs = xx(1) + uvw(1,1) * ddx * eps
               xf = xx(2) - uvw(1,1) * ddx * eps
               ys = yy(1) + uvw(2,2) * ddy * eps
               yf = yy(2) - uvw(2,2) * ddy * eps
               zs = zz(1) + uvw(3,3) * ddz * eps
               zf = zz(2) - uvw(3,3) * ddz * eps

               voly = 0.0
               volx = 0.0
               volz = 0.0

*-----------------------------------------------------------------------
*     scan z direction
*-----------------------------------------------------------------------

                  u = uvw(3,1)
                  v = uvw(3,2)
                  w = uvw(3,3)

         do 100 i = 1, nsn
         do 200 j = 1, nsn

                  xi = xm(i)
                  yi = ym(j)
                  zi = zs

                  ici   = -1
                  mark  =  1
                  markp =  0

                  call gomsor(xi,yi,zi,u,v,w,
     &                        nmed,iblz,mark,markp,ici)

               if( mark .ge. -1 ) then

                  nmed0 = nmed
                  mark  = 1
                  markp = 1

               else

                  goto 1000

               end if

*-----------------------------------------------------------------------

  300       continue

                  xc = xm(i)
                  yc = ym(j)
                  zc = zf

                  call gomprp(0,xi,yi,zi,xc,yc,zc,u,v,w,
     &                        nmed,iblz,mark,markp)

                        if( mark .lt. -1 ) goto 1000

                  if( nmed0 .gt. 0 ) then

                     do k = 1, nl

                        if( ( ilt .gt. 0 .and.
     &                        idmn(nmed0) .eq. lt(k) ) .or.
     &                      ( ilt .lt. 0 .and.
     &                        idmn(nmed0) .ne. lt(k) ) ) then

                           volz = volz + zc - zi

                           goto 400

                        end if

                     end do

  400                continue

                  end if

                  if( mark .eq. 1 .or. mark .eq. -1 ) then

                           goto 200

                  else if( mark .eq. 0 .or. mark .eq. 2 ) then

                           nmed0 = nmed

                           zi = zc

                           goto 300

                  end if

*-----------------------------------------------------------------------

  200    continue
  100    continue

               volz = volz * ddx * ddy

*-----------------------------------------------------------------------
*        no boundary
*-----------------------------------------------------------------------

         if( abs( volz - vls0 ) .lt. epss ) goto 1000

*-----------------------------------------------------------------------
*     scan x direction
*-----------------------------------------------------------------------

                  u = uvw(1,1)
                  v = uvw(1,2)
                  w = uvw(1,3)

         do 110 i = 1, nsn
         do 210 j = 1, nsn

                  xi = xs
                  yi = ym(i)
                  zi = zm(j)

                  ici   = -1
                  mark  =  1
                  markp =  0

                  call gomsor(xi,yi,zi,u,v,w,
     &                        nmed,iblz,mark,markp,ici)

               if( mark .ge. -1 ) then

                  nmed0 = nmed
                  mark  = 1
                  markp = 1

               else

                  goto 1000

               end if

*-----------------------------------------------------------------------

  310       continue

                  xc = xf
                  yc = ym(i)
                  zc = zm(j)

                  call gomprp(0,xi,yi,zi,xc,yc,zc,u,v,w,
     &                        nmed,iblz,mark,markp)

                        if( mark .lt. -1 ) goto 1000

                  if( nmed0 .gt. 0 ) then

                     do k = 1, nl

                        if( ( ilt .gt. 0 .and.
     &                        idmn(nmed0) .eq. lt(k) ) .or.
     &                      ( ilt .lt. 0 .and.
     &                        idmn(nmed0) .ne. lt(k) ) ) then

                           volx = volx + xc - xi

                           goto 410

                        end if

                     end do

  410                continue

                  end if

                  if( mark .eq. 1 .or. mark .eq. -1 ) then

                           goto 210

                  else if( mark .eq. 0 .or. mark .eq. 2 ) then

                           nmed0 = nmed

                           xi = xc

                           goto 310

                  end if

*-----------------------------------------------------------------------

  210    continue
  110    continue

               volx = volx * ddy * ddz

*-----------------------------------------------------------------------
*        no boundary
*-----------------------------------------------------------------------

         if( abs( volx - vls0 ) .lt. epss ) goto 1000

*-----------------------------------------------------------------------
*     scan y direction
*-----------------------------------------------------------------------

                  u = uvw(2,1)
                  v = uvw(2,2)
                  w = uvw(2,3)

         do 120 i = 1, nsn
         do 220 j = 1, nsn

                  xi = xm(i)
                  yi = ys
                  zi = zm(j)

                  ici   = -1
                  mark  =  1
                  markp =  0

                  call gomsor(xi,yi,zi,u,v,w,
     &                        nmed,iblz,mark,markp,ici)

               if( mark .ge. -1 ) then

                  nmed0 = nmed
                  mark  = 1
                  markp = 1

               else

                  goto 1000

               end if

*-----------------------------------------------------------------------

  320       continue

                  xc = xm(i)
                  yc = yf
                  zc = zm(j)

                  call gomprp(0,xi,yi,zi,xc,yc,zc,u,v,w,
     &                        nmed,iblz,mark,markp)

                        if( mark .lt. -1 ) goto 1000

                  if( nmed0 .gt. 0 ) then

                     do k = 1, nl

                        if( ( ilt .gt. 0 .and.
     &                        idmn(nmed0) .eq. lt(k) ) .or.
     &                      ( ilt .lt. 0 .and.
     &                        idmn(nmed0) .ne. lt(k) ) ) then

                           voly = voly + yc - yi

                           goto 420

                        end if

                     end do

  420                continue

                  end if

                  if( mark .eq. 1 .or. mark .eq. -1 ) then

                           goto 220

                  else if( mark .eq. 0 .or. mark .eq. 2 ) then

                           nmed0 = nmed

                           yi = yc

                           goto 320

                  end if

*-----------------------------------------------------------------------

  220    continue
  120    continue

               voly = voly * ddx * ddz

*-----------------------------------------------------------------------
*        no boundary
*-----------------------------------------------------------------------

         if( abs( voly - vls0 ) .lt. epss ) goto 1000

*-----------------------------------------------------------------------
*        average of three scans
*-----------------------------------------------------------------------

            vlsdc = ( volx + voly + volz ) / 3.0

*-----------------------------------------------------------------------

            return

 1000 continue

            vlsdc = vls0

*-----------------------------------------------------------------------

      return
      end function vlsdc
CCSE add new function for mesh=xyz (2018.07.31) <<<<<

************************************************************************
*                                                                      *
      subroutine echtetdc(iot,nr,mr,kr,mtscore)
*                                                                      *
*       input echo for tally tetra mesh                                *
*       Modified by T.Furuta on 2025/01/17                             *
*       orignal sub. is "echtet" in tallsm1.f                          *
*                                                                      *
************************************************************************
      use TETRAMOD, only: ielem2icl,nelem
      use moddas_material

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'param01.inc'     ! maxpt

*-----------------------------------------------------------------------

      integer,intent(in) :: iot,nr,mr,mtscore
      integer,intent(in) :: kr(mr)
      character ch0*100,ch1*100
      integer ind,kk,kkk,iz,ia
      integer ir,ielem0,icl,mtcnt
      integer imat,kmat,lmat,nel,nelh
      real(8) denh,denstm,den

*-----------------------------------------------------------------------

      common /kmat1a/ mxmat, mxmat0, mxnel
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /kmat1i/ kmatc(kvlmax)
      common /natnuc/ natnn(maxpt), natnm(maxpt,10), patnn(maxpt,10)
      common /regdm/  idmg(kvlmax)

*-----------------------------------------------------------------------
      integer,allocatable :: lr(:)
      real(8),allocatable :: vl(:)
*-----------------------------------------------------------------------
      character dum*100
      character chau*8
      character elmnt(104)*3

      data elmnt /
     &    ' H ','He ','Li ','Be ',' B ',' C ',' N ',' O ',' F ','Ne ',
     &    'Na ','Mg ','Al ','Si ',' P ',' S ','Cl ','Ar ',' K ','Ca ',
     &    'Sc ','Ti ',' V ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ',' Y ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ',' I ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ',' W ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ',' U ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /
*-----------------------------------------------------------------------

      integer,allocatable :: matlst(:),imatelem(:),imat2mtcnt(:)
      logical,allocatable :: imatflag(:)
      real(8),allocatable :: xcm(:,:) !FURUTA20200529

*-----------------------------------------------------------------------
*     allocation
*-----------------------------------------------------------------------

      allocate ( lr(nr),vl(nr) )
      allocate( matlst(mxmat),imatflag(mxmat) )
      allocate( imatelem(nr),imat2mtcnt(mxmat) )
      allocate( xcm(3,nr) )

*-----------------------------------------------------------------------
*        input echo for tally region mesh
*-----------------------------------------------------------------------

      ireg=kr(1)
      itet0=kr(2)
      nr0=kr(3)
      if(nr0.gt.0)then
       dum(1:2)='( '
       l=3
       do ir=1,nr0
        write(dum(l:l+6),'(i7)') kr(3+ir)
        l=l+7
        dum(l:l)=' '
        l=l+1
       enddo
       dum(l:l+1)='< '
       l=l+2
       write(dum(l:l+6),'(i7)') ireg
       l=l+7
       dum(l:l+1)=' )'
       l=l+1
      else
       l=1
       write(dum(l:l+6),'(i7)') ireg
       l=l+6
      endif
      write(iot,'("      reg = ",100a1)')
     &     (dum(j:j),j=1,l)

*-----------------------------------------------------------------------
*        set mesh volume
*-----------------------------------------------------------------------

      call ttetvl(mr,kr,nr,vl,lr)

*-----------------------------------------------------------------------
*        set CM coordinates of tetrahedrons !FURUTA20200529
*-----------------------------------------------------------------------

      call ttetcm(mr,kr,nr,xcm)

*-----------------------------------------------------------------------
*     check whether materals exist in whole mesh
*-----------------------------------------------------------------------

      imatflag(1:mxmat)=.false.
      do ir=1,nr
         ielem0=nelem(itet0-1)
         icl=ielem2icl(ielem0+ir)
         imat=idnm(idmg(icl))
         imatelem(ir)=imat
         imatflag(imat)=.true.
      enddo
      mtcnt=0
      do imat=1,mxmat
         if(imatflag(imat))then
          mtcnt=mtcnt+1
          matlst(mtcnt)=imat
          imat2mtcnt(imat)=mtcnt
         else
          imat2mtcnt(imat)=0
         endif
      enddo

*-----------------------------------------------------------------------
*     write material part
*     mtscore will be change to input parameter
*-----------------------------------------------------------------------

      write(iot,'("      imat =",i7)') mtcnt
      write(iot,'("   mtscore =",i7)') mtscore

      do imat = 1, mtcnt

         kmat = matlst(imat)
         nel  = nint( das_kmatc(kmatc(kmat)+1) )
         denh =       das_kmatc(kmatc(kmat)+2)
         nelh = nel

         if( denh .ne. 0.0d0 ) then
            nelh = nelh + 1
         end if

         lmat = nelh
         denstm = 0.0d0
         do kk = 1, nel
            iz   = nint( das_kmatc(kmatc(kmat)+(kk-1)*4+34) )
            ia   = nint( das_kmatc(kmatc(kmat)+(kk-1)*4+35) )
            den  =       das_kmatc(kmatc(kmat)+(kk-1)*4+36)

            if( ia .eq. 0 .and. natnn(iz) .ne. 0 ) then
               lmat = lmat + natnn(iz) - 1
            end if
         end do

         write(iot,'("   mt-list =",i7)') lmat

         if( denh .ne. 0.0d0 ) then
            chau = ' H-1'
            write(iot,'(8x,a4,8x,1p1e15.7)') chau(1:4), denh
         end if

         do kk = 1, nel
            iz   = nint( das_kmatc(kmatc(kmat)+(kk-1)*4+34) )
            ia   = nint( das_kmatc(kmatc(kmat)+(kk-1)*4+35) )
            den  =       das_kmatc(kmatc(kmat)+(kk-1)*4+36)
            chau = elmnt(iz)

            if( ia .eq. 0 ) then

               do kkk = 1, natnn(iz)

                  if( iz .eq. 1 .and. natnm(iz,kkk) .eq. 1 ) then

                     write(iot,'(8x,a2,"-1"8x,1p1e15.7)')
     &                chau(1:2), denh + den * patnn(iz,kkk) / 100.0

                  else

                     denstm = den * patnn(iz,kkk) / 100.0

                     if( natnm(iz,kkk) .lt. 10 ) then
                        write(iot,'(8x,a2,"-",i1,8x,1p1e15.7)')
     &                   chau(1:2), natnm(iz,kkk), denstm
                     else if( natnm(iz,kkk) .lt. 100 ) then
                        write(iot,'(8x,a2,"-",i2,7x,1p1e15.7)')
     &                   chau(1:2), natnm(iz,kkk), denstm
                     else
                        write(iot,'(8x,a2,"-",i3,6x,1p1e15.7)')
     &                   chau(1:2), natnm(iz,kkk), denstm
                     end if

                  end if

               end do

            else if( ia .lt. 10 ) then
               write(iot,'(8x,a2,"-",i1,8x,1p1e15.7)')
     &          chau(1:2), ia, den
            else if( ia .lt. 100 ) then
               write(iot,'(8x,a2,"-",i2,7x,1p1e15.7)')
     &          chau(1:2), ia, den
            else
               write(iot,'(8x,a2,"-",i3,6x,1p1e15.7)')
     &          chau(1:2), ia, den
            end if

         end do

      enddo

*-----------------------------------------------------------------------
*     write mesh part without fluxs ... fluxs will be added in tdctoinp
*-----------------------------------------------------------------------

      write(iot,'("     imesh = ",i7)') nr

       write(iot,'("   non tetra    mt",
     &"          xCM          yCM          zCM       volume")')
      do ir = 1, nr
       write(iot,'(3i6,1p20e13.4)')
     &       ir, lr(ir), imat2mtcnt(imatelem(ir)), xcm(1:3,ir), vl(ir)
      enddo

*-----------------------------------------------------------------------

      deallocate( lr,vl )
      deallocate( matlst,imatflag )
      deallocate( imatelem,imat2mtcnt )
      deallocate( xcm )

*-----------------------------------------------------------------------

      return
      end subroutine echtetdc

************************************************************************
*                                                                      *
      subroutine ttetmeshdc(jsn,jsi,dsin,idsi,ill,ilf,
     &                   jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                   mtrn,ntrn,
     &                   ndim_region,idas_region,
     &                   idcrg)
*                        ^^^^^  for [T-DCHAIN]
*                                                                      *
*       read 'tet =' sub-section of input tally section                *
*       modified by T.Furuta on 2025/01/17                             *
*       orignal sub. is "ttetmesh" in tallsm1.f                        *
*                                                                      *
************************************************************************
      use moddas

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character(200),intent(out) :: chin, chlw, chcm

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

      integer ipar

      data klnmax /1000000/

      integer, intent(in) :: ndim_region
      integer, intent(inout) :: idas_region(ndim_region)
      integer,allocatable :: mtetreg(:)

*-----------------------------------------------------------------------

            ierr = 0
            idcrg = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

               idcrg = idcrg + 1
               if(idcrg.gt.maxreg)call reallocate_talldc4(idcrg) !FURUTA20200522
               dcrgm(idcrg) = chin

  150 continue

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

*-----------------------------------------------------------------------
*        error not for reg =
*-----------------------------------------------------------------------

               if( chcm(i1:i1+3) .ne. 'reg=' ) goto 998

               ic = inumc(chlw,i1,i3,'=') + 1

*-----------------------------------------------------------------------
*        read region number
*-----------------------------------------------------------------------

               ic = jnumc(chlw,ic,i3)

               ifis = 0
               ipar = 0
               jpar = 0
               kpar = 0

               call moddas_allocate_int(MAX_NUM_MTRG,mtetreg)

            do i = 1, klnmax

                  ic = jnumc(chlw,ic,i3)

                  call ttetmesh3(icn,chlw,i1,i3,ic,
     &                 ipar,jpar,kpar,ifis,
     &                 max_num_MTRG,mtetreg)

                     if( icn .gt. 900 ) goto 900
                     if( icn .eq. 500 ) goto 500

               if( i .lt. klnmax .and. ic .gt. i3 ) then

  147             call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) goto 500

                     if( iskip .ne. 0 ) goto 147

                     ifis = ifis + 1

                     ic = i1

               end if

            end do

                  goto 996

*-----------------------------------------------------------------------
*        modify the input
*-----------------------------------------------------------------------

  500    continue

         if(ipar.gt.1)then
          icn=996
          goto 900
         endif

         ntrn = kpar + 1
         mtrn = 3+kpar

         ndsm=1
         do k=1,mtrn
          idas_region( ndsm + k -1) = mtetreg(k)
         enddo

         call moddas_deallocate_int(mtetreg)

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

         return

  900 continue

         if( icn .eq. 996 ) goto 996
         if( icn .eq. 997 ) goto 997

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Only 1 region is allowed with mesh = tet'
         ErrCha = ''
         ErrID = 'L:13907/R:ttetmeshdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of region number is wrong.'
         ErrCha = ''
         ErrID = 'L:13919/R:ttetmeshdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue
         m_err = 'After { mesh = tet } line should be { reg = }.'
         ErrCha = ''
         ErrID = 'L:13930/R:ttetmeshdc/F:tallsm4.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1

         return

*-----------------------------------------------------------------------

      end subroutine ttetmeshdc

************************************************************************
*                                                                      *
      subroutine check_itdc2reg(m,nr)
*                                                                      *
*       check itdc2reg index which is related to talldc2 allocation    *
*       modified by T.Furuta on 2020/06/23                             *
*                                                                      *
************************************************************************
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      integer,intent(in) :: m,nr
      integer :: itdc0,i,n0,n1,nl
      itdc0=itnm2tdc(m-1)
      if(itdc0.lt.itdc)then
       if(itdc2reg(itdc0+1)-itdc2reg(itdc0).lt.nr)then
        n0=itdc2reg(itdc0)
        n1=itdc2reg(itdc0+1)
        nl=itdcreg-n1
        itdcreg=itdcreg+nr-n1+n0
        do i=itdc0+1,itdc
         itdc2reg(i)=itdc2reg(i)+nr-n1+n0
        enddo
        if(itdcreg.gt.maxtdcreg)call reallocate_talldc2(itdcreg)
        if(nl.gt.0)then
         allocate( itgceltmp(maxtdcreg) )
         allocate( rtgvlmtmp(maxtdcreg) )
         allocate( itgnumtmp(maxtdcreg) )
         allocate( ireg2itmtmp(maxtdcreg) )
         allocate( ffluxtmp(maxtdcreg) )
         itgceltmp(1:maxtdcreg)=itgcel(1:maxtdcreg)
         rtgvlmtmp(1:maxtdcreg)=rtgvlm(1:maxtdcreg)
         itgnumtmp(1:maxtdcreg)=itgnum(1:maxtdcreg)
         ireg2itmtmp(1:maxtdcreg)=ireg2itm(1:maxtdcreg)
         ffluxtmp(1:maxtdcreg)=fflux(1:maxtdcreg)
         itgcel(n0+nr+1:n0+nr+nl)=itgceltmp(n1+1:n1+nl)
         rtgvlm(n0+nr+1:n0+nr+nl)=rtgvlmtmp(n1+1:n1+nl)
         itgnum(n0+nr+1:n0+nr+nl)=itgnumtmp(n1+1:n1+nl)
         ireg2itm(n0+nr+1:n0+nr+nl)=ireg2itmtmp(n1+1:n1+nl)
         fflux(n0+nr+1:n0+nr+nl)=ffluxtmp(n1+1:n1+nl)
         deallocate(itgceltmp,rtgvlmtmp,itgnumtmp,ireg2itmtmp,ffluxtmp)
        endif
       endif
      endif
*-----------------------------------------------------------------------
      if(nr.gt.maxreg)call reallocate_talldc4(nr)
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
      end subroutine check_itdc2reg
*-----------------------------------------------------------------------

      end module TDCHAINMOD
