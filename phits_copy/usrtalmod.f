
      module usrtalmod

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      real(8),allocatable,save,target:: udtvar(:)
      integer,allocatable,save,target:: iudtvar(:)

*-----------------------------------------------------------------------
      contains


*-----------------------------------------------------------------------
      subroutine ALLOCATE_udtvar

      implicit real*8 (a-h,o-z)

      integer nudtvar
      common /cnudtvar/ nudtvar

      allocate(udtvar(nudtvar))
      allocate(iudtvar(nudtvar))
      udtvar(1:nudtvar) = 0.0d0
      iudtvar(1:nudtvar) = 0

      end subroutine ALLOCATE_udtvar
*-----------------------------------------------------------------------


      end module
