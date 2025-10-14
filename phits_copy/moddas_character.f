!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!module moddas_character
      module moddas_character

      character(len=:), allocatable :: chrg
      character(len=:), allocatable :: chtm

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!contains
      contains

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_character_initialize()

      implicit real*8 (a-h,o-z)

      include 'param.inc'

      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      end module
