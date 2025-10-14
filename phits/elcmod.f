! AdvanceSoft H.Hasemi 2019/11/11, module for electric field map
      module elcmod

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      real(8),allocatable,save,target:: exlist(:,:,:),exmesh(:)
      real(8),allocatable,save,target:: eylist(:,:,:),eymesh(:)
      real(8),allocatable,save,target:: ezlist(:,:,:),ezmesh(:)

      real(8),allocatable,save,target:: erzrlist(:,:),erzrmesh(:)
      real(8),allocatable,save,target:: erzzlist(:,:),erzzmesh(:)

      real(8),allocatable,save,target:: exmap(:,:,:)
      real(8),allocatable,save,target:: eymap(:,:,:)
      real(8),allocatable,save,target:: ezmap(:,:,:)

      real(8),allocatable,save,target:: erzrmap(:,:),erzzmap(:,:)

*-----------------------------------------------------------------------
      contains

*-----------------------------------------------------------------------
      subroutine ALLOCATE_exyzlist(nx,ny,nz)

            implicit real*8 (a-h,o-z)

            allocate(exlist(nx,ny,nz))
            allocate(eylist(nx,ny,nz))
            allocate(ezlist(nx,ny,nz))
            exlist(:,:,:)=0.0
            eylist(:,:,:)=0.0
            ezlist(:,:,:)=0.0
            allocate(exmesh(nx))
            allocate(eymesh(ny))
            allocate(ezmesh(nz))
            exmesh(:)=0.0
            eymesh(:)=0.0
            ezmesh(:)=0.0

      end subroutine ALLOCATE_exyzlist
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
      subroutine ALLOCATE_erzlist(nr,nz)

            implicit real*8 (a-h,o-z)

            allocate(erzrlist(nr,nz))
            allocate(erzzlist(nr,nz))
            erzrlist(:,:)=0.0
            erzzlist(:,:)=0.0
            allocate(erzrmesh(nr))
            allocate(erzzmesh(nz))
            erzrmesh(:)=0.0
            erzzmesh(:)=0.0

      end subroutine ALLOCATE_erzlist
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
      subroutine ALLOCATE_exyzmap(nx,ny,nz)

            implicit real*8 (a-h,o-z)

            allocate(exmap(nx,ny,nz))
            allocate(eymap(nx,ny,nz))
            allocate(ezmap(nx,ny,nz))

            exmap(:,:,:)=0.0
            eymap(:,:,:)=0.0
            ezmap(:,:,:)=0.0

      end subroutine ALLOCATE_exyzmap
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
      subroutine ALLOCATE_erzmap(nr,nz)

            implicit real*8 (a-h,o-z)

            allocate(erzrmap(nr,nz))
            allocate(erzzmap(nr,nz))
            erzrmap(:,:)=0.0
            erzzmap(:,:)=0.0

      end subroutine ALLOCATE_erzmap
*-----------------------------------------------------------------------


      end module
