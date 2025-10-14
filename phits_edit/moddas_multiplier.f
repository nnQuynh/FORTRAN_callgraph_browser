!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!module moddas_multiplier
      module moddas_multiplier

      double precision, allocatable :: gmsh_ismlt(:)

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!contains
      contains

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_multiplier_initialize()

      implicit real*8 (a-h,o-z)
      include 'param.inc'

      common /multipl/ imltp,iimlt(multmax),inmlt(multmax),
C MATSUDA 2024.11.25 (lagrange: ilmlt)
     &  ilmlt(multmax),
     &  idmlt(multmax),ismlt(multmax),impan(multmax),impat(multmax,6,3),
     &                 jmpat(multmax,6,6,2), imdfl(multmax)

      ismlt(:) = 1
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      end module
