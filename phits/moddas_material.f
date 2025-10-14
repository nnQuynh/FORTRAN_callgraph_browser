!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!module moddas_material
      module moddas_material

      double precision, allocatable :: das_kmat(:)
      integer :: num_kmat_add
      integer :: kmat0 = 0

      double precision, allocatable :: dnel_das(:)
      double precision, allocatable :: denh_das(:)
      double precision, allocatable :: zz_das(:)
      double precision, allocatable :: a_das(:)
      double precision, allocatable :: den_das(:)

      double precision, allocatable :: das_kmatg(:)
      double precision, allocatable :: das_kmatc(:)
      double precision, allocatable :: das_kmatd(:)
      double precision, allocatable :: das_kmate(:)
      double precision, allocatable :: das_kmathd(:)
      double precision, allocatable :: das_kmathe(:)

      double precision, allocatable :: sigg(:)
      double precision, allocatable :: siggm(:)

      double precision, allocatable :: edns(:)

      double precision, allocatable :: das_intum(:)

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!contains
      contains

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_material_initialize()

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      common /kmat1g/ kmat(kvlmax)
      common /kmat1h/ kmatg(kvlmax)
      common /kmat1i/ kmatc(kvlmax)
      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax)
      common /kmat1ka/ kmathd(kvlmax), kmathe(kvlmax)
      common /xgeosm/ ksig(kvlmax)
      common /delreg/ delm(kvlmax), kdelt
      common /kmat1j/ intum

      kmat(:) = 0
      kmatg(:) = 0
      kmatc(:) = 0
      kmatd(:) = 0
      kmate(:) = 0
      kmathd(:) = 0
      kmathe(:) = 0
      ksig(:) = 0
      kdelt = 0
      intum = 0
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      end module
