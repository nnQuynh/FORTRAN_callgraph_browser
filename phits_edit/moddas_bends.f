!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!module moddas_bends
      module moddas_bends
      integer :: iflag=0
      double precision, allocatable :: bends(:,:)
!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      contains
!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_bends_allocate(n)
      implicit none
      integer,intent(in) :: n
      iflag=1
      allocate( bends(4,n) )
      end subroutine moddas_bends_allocate
!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_bends_deallocate()
      implicit none
      if(iflag.eq.1)deallocate( bends )
      end subroutine moddas_bends_deallocate
!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      end module




