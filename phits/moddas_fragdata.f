!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!module moddas_fragdata
      module moddas_fragdata

      integer, allocatable :: ifgs01_nei(:)
      integer, allocatable :: ifgs02_kne(:)
      integer, allocatable :: ifgs03_kxs(:)
      integer, allocatable :: ifgs04_neo(:)
      integer, allocatable :: ifgs05_kef(:)
      integer, allocatable :: ifgs06_nag(:)
      integer, allocatable :: ifgs07_kaf(:)
      integer, allocatable :: ifgs08_nfrg(:)
      integer, allocatable :: ifgs09_kim(:)
      integer, allocatable :: ifgs10_ks0(:)

      double precision, allocatable :: frgne(:)
      double precision, allocatable :: frgxs(:)
      double precision, allocatable :: frgef(:)
      double precision, allocatable :: frgaf(:)
      integer, allocatable :: ifrgm(:)

      integer, allocatable :: ifrge_ksf(:)
      integer, allocatable :: ifrge_ks1(:)
      integer, allocatable :: ifrge_ks2(:)
      integer, allocatable :: ifrge_ks3(:)
      double precision, allocatable :: frgsf(:)
      double precision, allocatable :: frgee(:)
      double precision, allocatable :: frgdd(:)
      double precision, allocatable :: frgdx(:)

      integer :: mmlm = 0
      integer :: num_frgee = 0
      integer :: num_frgdd = 0
      integer :: num_frgdx = 0

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!contains
      contains

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_fragdata_initialize()

      implicit real*8 (a-h,o-z)

      include 'param.inc'

      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      end module
