************************************************************************
      module mod_counter
!     Module for event counter
!     create by NS 2020.04
!     modified by T. FURUTA on 2021/01/19
************************************************************************
      implicit none
      integer,save:: iTYinevt,iTYiregn,iTYmxmat ! set at read00
      integer,save:: iaevt, ibevt, jaevt, jbevt ! zero set at ovly12

      real(8), dimension(40)  :: rncnt
      real(8), dimension(200) :: rnint, rnintr
      real(8), dimension(40)  :: rnpnt, rnpntr
C for  REDUCTION_COUNTER
!$OMP THREADPRIVATE(rncnt)
!$OMP THREADPRIVATE(rnint, rnintr)
!$OMP THREADPRIVATE(rnpnt, rnpntr)
      real(8), dimension(40)  :: rncnt2
      real(8), dimension(200) :: rnint2, rnintr2
      real(8), dimension(40)  :: rnpnt2, rnpntr2

      real(8),allocatable,save:: aevts(:,:),aevtr(:,:)
      real(8),allocatable,save:: bevts(:,:),bevtr(:,:)
C for  REDUCTION_COUNTER
!$OMP THREADPRIVATE(aevts, bevts)
!$OMP THREADPRIVATE(aevtr, bevtr)
      real(8),allocatable,save::aevts2(:,:),aevtr2(:,:)
      real(8),allocatable,save::bevts2(:,:),bevtr2(:,:)

      contains
!------------------------------------------------------------------------
      subroutine ALLOCATE_EVTS
      allocate(aevts(iTYinevt,iTYiregn),aevtr(iTYinevt,iTYiregn))
      allocate(bevts(iTYinevt,iTYmxmat),bevtr(iTYinevt,iTYmxmat))

      end subroutine ALLOCATE_EVTS
!------------------------------------------------------------------------
      subroutine ALLOCATE_EVTS2
      allocate(aevts2(iTYinevt,iTYiregn),aevtr2(iTYinevt,iTYiregn))
      allocate(bevts2(iTYinevt,iTYmxmat),bevtr2(iTYinevt,iTYmxmat))

      end subroutine ALLOCATE_EVTS2
!------------------------------------------------------------------------
      subroutine DEALLOCATE_EVTS

      deallocate(aevts,aevtr,bevts,bevtr)

      end subroutine DEALLOCATE_EVTS
!------------------------------------------------------------------------
      subroutine DEALLOCATE_EVTS2

      deallocate(aevts2,aevtr2,bevts2,bevtr2)

      end subroutine DEALLOCATE_EVTS2
!------------------------------------------------------------------------
      subroutine INIT_EVTS

      rncnt  = 0.0d0
      rnint  = 0.0d0
      rnintr = 0.0d0
      rnpnt  = 0.0d0
      rnpntr = 0.0d0
      aevts  = 0.0d0
      bevts  = 0.0d0
      aevtr  = 0.0d0
      bevtr  = 0.0d0
!$OMP MASTER
      iaevt  = 0
      ibevt  = 0
      jaevt  = 0
      jbevt  = 0
!$OMP END MASTER
!$OMP BARRIER

      end subroutine INIT_EVTS
!------------------------------------------------------------------------
      subroutine INIT_EVTS2

      rncnt2  = 0.0d0
      rnint2  = 0.0d0
      rnintr2 = 0.0d0
      rnpnt2  = 0.0d0
      rnpntr2 = 0.0d0
      aevts2  = 0.0d0
      bevts2  = 0.0d0
      aevtr2  = 0.0d0
      bevtr2  = 0.0d0

      end subroutine INIT_EVTS2
!------------------------------------------------------------------------
      end module mod_counter
