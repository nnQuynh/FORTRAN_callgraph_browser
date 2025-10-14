!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!module moddas_variance_reduction
      module moddas_variance_reduction

      double precision, allocatable, target :: wtxwwp(:,:,:)
      double precision, allocatable, target :: wtximp(:,:,:)
      double precision, allocatable, target :: wtxcut(:,:,:)
      double precision, allocatable, target :: wtxfcl(:,:,:)

      double precision, allocatable :: das_isstr(:)

      double precision, pointer :: wtxwwp_pointer(:)
      double precision, pointer :: wtximp_pointer(:)
      double precision, pointer :: wtxcut_pointer(:)
      double precision, pointer :: wtxfcl_pointer(:)

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!contains
      contains

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_variance_reduction_initialize()

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww
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

      iswwp = 0
      isimp = 0
      iswct = 0
      isfcl = 0
      isrcl = 0
      isstr = 0
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      end module
