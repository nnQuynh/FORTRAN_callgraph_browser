C for noncritical_ovly12 and strict_rsouin
      module mod_rsouin

      implicit none

!$    real(8) rsouin_pre
!$    real(8), allocatable :: rsouin_arr(:)

      contains

        subroutine ALLOCATE_MODRSOUIN
          integer maxbch,maxcas
          common /cparm/  maxbch,maxcas

!$        allocate(rsouin_arr(maxcas))
!$        rsouin_pre  = 0.0d0
!$        rsouin_arr  = 0.0d0
        end subroutine

        subroutine DEALLOCATE_MODRSOUIN
!$        deallocate(rsouin_arr)
        end subroutine

      end module
