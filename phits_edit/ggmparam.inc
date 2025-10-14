************************************************************************
*                                                                      *
*     'ggmparam.inc'                                                   *  *
*                                                                      *
*   definition of parameter value for GG                               *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      parameter ( mipt = 34, mtop = 89, npkey = 6, mstp = 4,
     &            mbng = 51, maxi = 34, novr = 5, mink = 200,
     &            mpng = 21, mwng = ( mtop + 1 ) / 2,
     &            mxdx = 10, mcoh = 55, minc = 21 )
      parameter ( iud  = 7, iux = 8 )

      parameter ( maxsec = 8 )
      parameter ( mixs   = 12 )

      parameter (
     & aneut=1.008664967d0,
     & avogad=6.022043446928244d23,
     & avgdn=1.d-24*avogad/aneut,
     & euler=.577215664901532861d0,
     & fscon=137.0393d0,
     & planck=4.135732d-13,
     & slite=299.7925d0,
     & third=1.0d0/3.0d0 )

*-----------------------------------------------------------------------

c S.Hashimoto changed the order for mct  (2015.4.1)
      common /gm001/ mct, kpt(mipt), nmat, nmat1, mxe, mxe1,       !FURUTA
     &               mix, mnnm, mxt, indt, mcal, naw,         !FURUTA
     &               nxsc, ntal, npert, ipert, itotnu,        !FURUTA
     &               itfxs, nsr, mgegbt(mipt),                !FURUTA
     &               its30, ntop, nwng, nstp, istern,         !FURUTA
     &               jovr(novr), ink(mink), jgm(mipt),   !FURUTA
     &               lmb, igm, img, npikmt, lxs, nfer, lfatl, !FURUTA
     &               nsa, ksdef,                              !FURUTA
     &               jtty, itty, nkxs, kc8,                   !FURUTA
     &               nww(mipt), ndx(mipt),                    !FURUTA
     &               kcy, ikz, mjss, nilw, irt,               !FURUTA
     &               nch(mipt), in_erg(mipt)                  !FURUTA

      common /gm003/
     &               lxss, lexs

      common /gm004/ emx(mipt), emcf(mipt), dnb, bnum, xnum,   !FURUTA
     &               rnok, enum, ecf(mipt+1), efac,            !FURUTA
     &               elc(mipt), dxw(mipt,3), dxx(mipt,5,mxdx), !FURUTA
     &               tco(mipt)                                 !FURUTA

      common /gm005/ iunr, ides, iphot, ibad, istrg, numb,
     &               nocoh, ispn, nee

      common /gm006/ bbrem(mtop), thgf(0:50), rkt(mtop), rktc(mtop),
     &               gpt(mipt), rka(mbng), calph(maxi), dbcn(30),
     &               vco(mcoh), vic(minc), wco(mcoh),
     &               qfiss(23), mfiss(22)

      character      hnp(mipt)*11
      character*12   htn  ! Should be changed when nuclear data is add
      character*230   xsdir, hdpath, hdpth, klin
      common /gm007/ xsdir, hdpath, hdpth, klin, hnp, htn

      common /gm008/ xunru

      common /gm010/ kpik !FURUTA

      common /gm011/ mkc,ngp,ntyn,mpan,iet,ipsc,ixre,ixcos,nter !FURUTA
!$OMP THREADPRIVATE(/gm011/)

      common /gm012/ colout(3,11),rnb(5)                        !FURUTA
!$OMP THREADPRIVATE(/gm012/)

*-----------------------------------------------------------------------
