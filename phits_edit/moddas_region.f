!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!module moddas_region
      module moddas_region

! source
      integer, allocatable :: idas_nsrc(:)
      integer, allocatable :: idas_nsrc_temporary(:)
      integer, allocatable :: idas_nsrn(:)
      integer, allocatable :: iaddress_nsrn(:)

! importance
      double precision, allocatable :: das_kfimp(:)
      integer, allocatable :: idas_inimc(:)
      integer, allocatable :: idas_inimt(:)

! weight window
      double precision, allocatable :: das_kfwwp(:)
      double precision, allocatable :: das_kgwwp(:)
      integer, allocatable :: idas_inwwc(:)
      integer, allocatable :: idas_inwwt(:)

! ww bias
      double precision, allocatable :: das_kfwbp(:)
      integer, allocatable :: idas_inwbc(:)
      integer, allocatable :: idas_inwbt(:)

! forced_collisions
      double precision, allocatable :: das_kfcls(:)
      integer, allocatable :: idas_inflc(:)
      integer, allocatable :: idas_inflt(:)

! repeated_collisions
      integer, allocatable :: idas_krcls(:)
      integer, allocatable :: idas_lrcls(:)
      integer, allocatable :: idas_inrlc(:)
      integer, allocatable :: idas_inrlt(:)

      integer :: iflag_evap

! splitting
      double precision, allocatable :: das_ksplt(:)
      integer, allocatable :: idas_ipgrc(:)
      integer, allocatable :: idas_ipgrc_temporary1(:)
      integer, allocatable :: idas_ipgrc_temporary2(:)
      integer, allocatable :: idas_ipgrt(:)

! counter
      integer, allocatable :: idas_kcont(:)
      integer, allocatable :: idas_incrc(:)
      integer, allocatable :: idas_incrt(:)

! magnetic field
      double precision, allocatable :: das_kmags(:)
      integer, allocatable :: idas_ingrc(:)
      integer, allocatable :: idas_ingrt(:)

! electro magnetic field
      double precision, allocatable :: das_kelcs(:)
      integer, allocatable :: idas_inerc(:)
      integer, allocatable :: idas_inert(:)

! super mirror
      double precision, allocatable :: das_ksmir(:)
      integer, allocatable :: idas_isgrc(:)
      integer, allocatable :: idas_isgrc_temporary1(:)
      integer, allocatable :: idas_isgrc_temporary2(:)
      integer, allocatable :: idas_isgrt(:)

! timer
      integer, allocatable :: idas_ktime(:)
      integer, allocatable :: idas_intmc(:)
      integer, allocatable :: idas_intmt(:)

! tally reg=
      integer, allocatable :: idas_itreg(:)
      integer, allocatable :: idas_itreg_temporary(:)

      integer, allocatable :: idas_itrcg(:)


! tally reginbox=
      integer, allocatable :: idas_itmeg(:)
      integer, allocatable :: idas_itmeg_temporary(:)

      integer, allocatable :: idas_itmcg(:)

! t-cross
      integer, allocatable :: idas_itrcr(:)
      integer, allocatable :: idas_itrcr_temporary1(:)
      integer, allocatable :: idas_itrcr_temporary2(:)
      double precision, allocatable :: das_itrca(:)

      integer, allocatable :: idas_itrcc(:)

      integer, allocatable :: idas_ntcr(:)
      double precision, allocatable :: das_ktcr(:)

! work array sub.echrg2, echrg3 @ tallsm3.f
      integer, allocatable :: nrst(:)

! work array sub.tdepreg0 @ result.f
      integer, allocatable :: idas_region_temporary(:)

!...
      integer :: iaddress_region(2)

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!contains
      contains

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_region_initialize()

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      common /isorsr/ nsmx(isrc), nsrn(isrc), nsrc(isrc)
      common /impmsg/ iimpn, isimp, iswct, isstr, maxip,
     &                mnimp(7,0:20),
     &                kfimp(7), inimc(7), inimt(7)
      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww
      common /wwind0/ dvals, iwmsh, kecho, kgwwp(6), idval
      common /wwbiasn/ eenwb(6,100), iwbdp, mnwbp(6,0:20), kfwbp(6),
     &                 inwbc(6), ienwb(6), inwbt(6), iswbp, maxwb
      common /fclmsg/ ifcln, isfcl, maxrg,
     &                mnfcl(6,0:20), mrfcl(kvlmax),
     &                kfcls(6), inflc(6), inflt(6)
      common /rclmsg/ ircln, isrcl, maxrr,
     &                mnrcl(6,0:20), mrrcl(kvlmax),
     &                krcls(6), lrcls(6), inrlc(6), inrlt(6),
     &                irpem(6,2), irman(6), irmct(6), jsmat(6)
      common /splreg/ isptn, npreg(6), mnspt(6,0:20),
     &                ipgrc(6), ipgrt(6), ksplt(6), isplt(6),
     &                ispct(6,9), ispem(6,2)
      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /magreg/ nmreg, ingrc, ingrt, kmags, magtin, imgusr
      common /elmgrg/ nelctf, nereg, inerc, inert, kelcs
      common /smireg/ nsreg, isgrc, isgrt, ksmir, ismir
      common /tmtreg/ ntmrg, intmc, intmt, ktime

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall13/ itrcn(itlmax), itrcr(itlmax), itrca(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)
      common /tall26/ itrcs(itlmax), itrcc(itlmax), itrss(itlmax)

      nsrc(:) = 0
      allocate(iaddress_nsrn(isrc))
      iaddress_nsrn(:) = 0

      kfimp(:) = 1
      inimc(:) = 0
      inimt(:) = 0

      kfwwp(:) = 1
      kgwwp(:) = 1
      inwwc(:) = 1
      inwwt(:) = 0

      kfwbp(:) = 1
      inwbc(:) = 1
      inwbt(:) = 0

      kfcls(:) = 1
      inflc(:) = 0
      inflt(:) = 0

      krcls(:) = 1
      lrcls(:) = 1
      inrlc(:) = 0
      inrlt(:) = 0

      ksplt(:) = 0
      ipgrc(:) = 1
      ipgrt(:) = 1

      kcont(:) = 0
      incrc(:) = 0
      incrt(:) = 0

      kmags = 0
      ingrc = 0
      ingrt = 0

      kelcs = 0
      inerc = 0
      inert = 0

      ksmir = 0
      isgrc = 1
      isgrt = 1

      ktime = 0
      intmc = 0
      intmt = 0

      itreg(:) = 1
      itrcg(:) = 1

      itrcr(:) = 1
      itrca(:) = 1
      itrcc(:) = 1

      itriv(:) = 1
      itrrv(:) = 1
      allocate( idas_itreg(1) )

      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      end module
