************************************************************************
      module LAFDATAMOD
************************************************************************
      integer,allocatable,save :: laf(:,:)
      integer,save :: lafmax

      contains
!------------------------------------------------------------------------
      subroutine INIT_laf
      implicit none
      if(allocated(laf))then
       laf(1:3,1:lafmax)=0
      endif
      return
      end subroutine INIT_laf
!------------------------------------------------------------------------
      subroutine EXTEND_laf(n)
      implicit none
      integer,intent(in) :: n
      integer,allocatable :: lafbuff(:,:)
      if(.not. allocated(laf))then
       allocate(laf(3,n+1)) ! +1 to avoid 0 for special function
       lafmax=n+1
      else
       allocate(lafbuff(3,lafmax))
       lafbuff(1:3,1:lafmax)=laf(1:3,1:lafmax)
       deallocate(laf)
       allocate(laf(3,lafmax+n))
       laf(1:3,1:lafmax)=lafbuff(1:3,1:lafmax)
       lafmax=lafmax+n
       deallocate(lafbuff)
      endif
      return
      end subroutine EXTEND_laf
!------------------------------------------------------------------------
      subroutine DEALLOCATE_laf
      implicit none
      if(allocated(laf))deallocate(laf)
      return
      end subroutine DEALLOCATE_laf
!------------------------------------------------------------------------
      end module LAFDATAMOD
