!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!module moddas_tally
      module moddas_tally

! mother (t-yield, t-dpa, t-prodct, t-star, t-dchain)
      integer, allocatable :: ismat(:)
      integer, allocatable :: ismat_temporary(:)

! nucleus (t-yield, t-dchain)
      integer, allocatable :: isnuc(:)
      integer, allocatable :: isnuc_temporary(:)

! material (t-wwg, t-sed, t-deposit, t-let, t-track, t-adjnt, t-yield,
      integer, allocatable :: ismte(:)
      integer, allocatable :: ismte_temporary(:)

! material in box (t-3dshow)
      integer, allocatable :: ismte_itmbt(:)

! dpa library information
      integer, allocatable :: mlib(:,:)
      integer :: iaddress_mlib(3) = 1
      integer :: icurrent_mlib = 1

      double precision, allocatable :: flib(:)
      integer :: iaddress_flib(3) = 1
      integer :: icurrent_flib = 1

! multiplier information
      integer, allocatable :: mltp(:,:)
      integer :: iaddress_mltp(3) = 1
      integer :: icurrent_mltp = 1

      double precision, allocatable :: slib(:)
      integer :: iaddress_slib(3) = 1
      integer :: icurrent_slib = 1

! anatally

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!contains
      contains

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_tally_initialize()

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall33/ itln(itlmax,2), itli(itlmax,2), itlr(itlmax,2),
     &                rtdm(itlmax,2)
      common /tall42/ itmbn(itlmax), itmbt(itlmax)
      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)
      common /tall78/ itism(itlmax), itist(10,itlmax), itjst(10,itlmax),
     &                itkst(10,itlmax), itstt(itlmax), itsdd(itlmax)

      itmat(:) = 1
      itnuc(:) = 1
      itmtt(:) = 1
      itmbt(:) = 1

      itli(:,:) = 0
      itlr(:,:) = 0

      itmli(:,:) = 0

      itsdd(:) = 0

      allocate(ismat(1),isnuc(1),ismte(1),ismte_itmbt(1))

      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      end module
