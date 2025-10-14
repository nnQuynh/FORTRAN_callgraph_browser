! T.Sato 2019/01/14, module for magnetic field map
      module magmod

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      real(8),allocatable,save,target:: bxlist(:,:,:),bxmesh(:)
      real(8),allocatable,save,target:: bylist(:,:,:),bymesh(:)
      real(8),allocatable,save,target:: bzlist(:,:,:),bzmesh(:)

      real(8),allocatable,save,target:: brzrlist(:,:),brzrmesh(:)
      real(8),allocatable,save,target:: brzzlist(:,:),brzzmesh(:)

      real(8),allocatable,save,target:: bxmap(:,:,:)
      real(8),allocatable,save,target:: bymap(:,:,:)
      real(8),allocatable,save,target:: bzmap(:,:,:)

      real(8),allocatable,save,target:: brzrmap(:,:),brzzmap(:,:)

*-----------------------------------------------------------------------
      contains

*-----------------------------------------------------------------------
      subroutine ALLOCATE_bxyzlist(nx,ny,nz)

      implicit real*8 (a-h,o-z)

      allocate(bxlist(nx,ny,nz))
      allocate(bylist(nx,ny,nz))
      allocate(bzlist(nx,ny,nz))
      bxlist(:,:,:)=0.0
      bylist(:,:,:)=0.0
      bzlist(:,:,:)=0.0
      allocate(bxmesh(nx))
      allocate(bymesh(ny))
      allocate(bzmesh(nz))
      bxmesh(:)=0.0
      bymesh(:)=0.0
      bzmesh(:)=0.0

      end subroutine ALLOCATE_bxyzlist
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
      subroutine ALLOCATE_brzlist(nr,nz)

      implicit real*8 (a-h,o-z)

      allocate(brzrlist(nr,nz))
      allocate(brzzlist(nr,nz))
      brzrlist(:,:)=0.0
      brzzlist(:,:)=0.0
      allocate(brzrmesh(nr))
      allocate(brzzmesh(nz))
      brzrmesh(:)=0.0
      brzzmesh(:)=0.0

      end subroutine ALLOCATE_brzlist
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
      subroutine ALLOCATE_bxyzmap(nx,ny,nz)

      implicit real*8 (a-h,o-z)

      allocate(bxmap(nx,ny,nz))
      allocate(bymap(nx,ny,nz))
      allocate(bzmap(nx,ny,nz))

      bxmap(:,:,:)=0.0
      bymap(:,:,:)=0.0
      bzmap(:,:,:)=0.0

      end subroutine ALLOCATE_bxyzmap
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
      subroutine ALLOCATE_brzmap(nr,nz)

      implicit real*8 (a-h,o-z)

      allocate(brzrmap(nr,nz))
      allocate(brzzmap(nr,nz))
      brzrmap(:,:)=0.0
      brzzmap(:,:)=0.0

      end subroutine ALLOCATE_brzmap
*-----------------------------------------------------------------------


      end module
