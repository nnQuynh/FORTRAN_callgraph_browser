************************************************************************
      module GGMBANKMOD
!     pik() should be also include in principle
!     But now the dimmension of pik() is 0
!     and therefore it is not included here
************************************************************************
      implicit none
      integer,save:: nrtc_bank,ntgp_bank,nktc_bank

      real(8),allocatable,save:: rtc(:,:),tgp(:)
      integer,allocatable,save:: ktc(:,:)
      real(8),allocatable,save:: rtcel(:) ! frtati 2021/12/17
!$OMP THREADPRIVATE(rtc,tgp,ktc,rtcel) ! frtati 2021/12/17

      contains
!------------------------------------------------------------------------
      subroutine ALLOCATE_GGMBANK

      allocate(rtc(15,nrtc_bank),tgp(ntgp_bank),ktc(2,nktc_bank))
      allocate(rtcel(nrtc_bank)) ! frtati 2021/12/17

      end subroutine ALLOCATE_GGMBANK
!------------------------------------------------------------------------
      subroutine DEALLOCATE_GGMBANK

      deallocate(rtc,tgp,ktc,rtcel) ! frtati 2021/12/17

      end subroutine DEALLOCATE_GGMBANK
!------------------------------------------------------------------------
      subroutine INIT_GGMBANK

      rtc(1:15,1:nrtc_bank)=0.0d0
      ktc(1:2,1:nktc_bank)=1
      rtcel(1:nrtc_bank)=0.d0 ! frtati 2021/12/17

      end subroutine INIT_GGMBANK
!------------------------------------------------------------------------
      end module GGMBANKMOD
