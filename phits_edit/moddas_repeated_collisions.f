!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!module moddas_repeated_collisions
      module moddas_repeated_collisions

      integer, allocatable :: ismat_jsmat(:)

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!contains
      contains

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_repeated_collisions_initialize()

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      common /rclmsg/ ircln, isrcl, maxrr,
     &                mnrcl(6,0:20), mrrcl(kvlmax),
     &                krcls(6), lrcls(6), inrlc(6), inrlt(6),
     &                irpem(6,2), irman(6), irmct(6), jsmat(6)

      jsmat(:) = 0
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      end module
