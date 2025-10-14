************************************************************************
      module LATDATAMOD
************************************************************************
      integer,allocatable,save :: ldata(:)
      integer,allocatable :: ltmp(:)
      integer,save :: ilatind,nlatind

      contains
!------------------------------------------------------------------------
      subroutine INIT_ldata
      implicit none
      nlatind=0
      end subroutine INIT_ldata
!------------------------------------------------------------------------
      subroutine ALLOCATE_ldata(n)
      implicit none
      integer,intent(in) :: n
      nlatind=n
      ilatind=0
      allocate(ldata(nlatind))
      end subroutine ALLOCATE_ldata
!------------------------------------------------------------------------
      subroutine DEALLOCATE_ldata
      implicit none
      deallocate(ldata)
      end subroutine DEALLOCATE_ldata
!------------------------------------------------------------------------
      subroutine RESIZE_ldata(n)
      integer,intent(in) :: n
      integer :: m
      m=nlatind
      allocate(ltmp(m))
      ltmp(1:m)=ldata(1:m)
      deallocate(ldata)
      nlatind=nlatind+n
      allocate(ldata(nlatind))
      ldata(1:m)=ltmp(1:m)
      deallocate(ltmp)
      ilatind=m
      end subroutine RESIZE_ldata
*-----------------------------------------------------------------------
      subroutine REWIND_ldata
      ilatind=0
      end subroutine REWIND_ldata
!-----------------------------------------------------------------------
      end module LATDATAMOD
