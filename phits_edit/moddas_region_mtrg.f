!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!module moddas_region_mtrg
      module moddas_region_mtrg

! region data (sub.tregion, tregion2, tregion3, tregion4, tregion5, tdepreg)
      integer, allocatable :: mtrg(:,:)

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!contains
      contains

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_region_mtrg_initialize()

      implicit real*8 (a-h,o-z)

      include 'param.inc'

      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      end module
