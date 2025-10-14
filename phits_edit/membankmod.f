************************************************************************
      module MEMBANKMOD
************************************************************************
      implicit none
      integer,save:: mxmat_bank,maxnl_bank,mxnel_bank,msumnel_bank

      real(8),allocatable,save:: arg(:),zsmcs(:),zemcs(:)
      real(8),allocatable,save:: avam(:),avaz(:),avad(:)
      real(8),allocatable,save:: bari(:),s1(:),s2(:),s3(:)
      real(8),allocatable,save:: eion(:),prtrng(:,:),prtrnglog(:,:)
      real(8),allocatable,save:: rhoi(:)
      real(8),allocatable,save:: f(:),g(:),ro(:)
      real(8),allocatable,save:: rnge(:,:),erng(:)
      real(8),allocatable,save:: sigge(:),siggn(:)

      real(8),allocatable,save:: siggc(:),siggd(:)
      integer,allocatable,save:: isigc(:),isigd(:)
      real(8),allocatable,save:: siggcn(:),siggce(:),siggpn(:)
      real(8),allocatable,save:: siggcmx(:,:)
      integer,save:: isigpn0

      integer,allocatable,save:: isigza(:)

!$OMP THREADPRIVATE(arg,zsmcs,zemcs)
!$OMP THREADPRIVATE(avam,avaz,avad)
!$OMP THREADPRIVATE(bari,s1,s2,s3)
!$OMP THREADPRIVATE(eion,prtrng,prtrnglog)
!$OMP THREADPRIVATE(rhoi)
!$OMP THREADPRIVATE(f,g,ro)
!$OMP THREADPRIVATE(rnge,erng)
!$OMP THREADPRIVATE(sigge,siggn)

!$OMP THREADPRIVATE(siggc,siggd)
!$OMP THREADPRIVATE(isigc,isigd)
!$OMP THREADPRIVATE(siggcn,siggce,siggpn,siggcmx)
!$OMP THREADPRIVATE(isigpn0)

!$OMP THREADPRIVATE(isigza)

      !-- dbspaces are shared. start ---!
      real(8),allocatable,save:: dbfspace  (:,:)
      real(8),allocatable,save:: dbfspacene(:,:)
      real(8),allocatable,save:: dbfspacens(:,:)
      integer(8),allocatable,save:: dbispace  (:,:)
      integer(8),allocatable,save:: dbispacene(:,:)
      integer(8),allocatable,save:: dbispacens(:,:)
      integer(4),allocatable,save:: dbfispace  (:,:)
      integer(4),allocatable,save:: dbfispacene(:,:)
      integer(4),allocatable,save:: dbfispacens(:,:)
      !-- dbspaces are shared. end   ---!

      integer(8),save:: mdbatima,ndbatima(3)=0
      real(8),save:: dbcutoff
      integer(8),parameter:: kspc=10000,ispc=1000

cfrtati 2022/12/20 database for pseudo cross-section
      real(8),allocatable,save :: dbpscs_signe(:), dbpscs_sigel(:)
      real(8),allocatable,save :: dbpscs_epeak(:)
      integer,allocatable,save :: dbpscs_iza(:,:)
      integer,save :: mdbpseud
      integer,save :: ndbcurrent = 0
      logical,save :: ldodbps = .false., laddbps = .false.

c----------------------------------------------------------------------
      contains
!------------------------------------------------------------------------
      subroutine ALLOCATE_MEMBANK
      implicit none
      integer nspred,nwsprd,nedisp,itstep,ndedx
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx
      integer mxmat,mxnel,maxnl,msumnel

      mxmat=mxmat_bank
      mxnel=mxnel_bank

      allocate(siggc(mxnel),siggd(mxnel))
      allocate(isigc(0:mxnel),isigd(0:mxnel))
      allocate(siggcn(mxnel),siggce(mxnel),siggpn(mxnel))
      allocate(siggcmx(4,mxnel))
      siggd = 0.d0 ! frtati 2022/03/11

      allocate(isigza(mxnel))

      allocate(sigge(mxnel),siggn(mxnel+1))
      sigge = 0.d0
      siggn = 0.d0
      allocate(rnge(2048,mxmat),erng(2048))
      if( nedisp .ne. 0 )then
        allocate(avam(mxmat),avaz(mxmat),avad(mxmat))
      endif
      if( nedisp .eq. 10 .or. nspred .eq. 10 .or.
     &     ndedx .eq. 1  .or.  ndedx .eq. 3 )then   ! S.Abe 2016/08/08
        allocate(rhoi(mxmat))
      endif
      if( ndedx .le. 3 )then   ! S.Abe 2016/08/08
        msumnel=msumnel_bank
        maxnl=maxnl_bank
        allocate(bari(mxmat),s1(mxmat),s2(mxmat),s3(mxmat))
        allocate(eion(msumnel),prtrng(341,mxmat),prtrnglog(341,mxmat))
        allocate(f(maxnl),g(maxnl),ro(maxnl+1))
      endif
      if( nspred .ne. 0 )then
        allocate(arg(mxmat),zsmcs(mxmat),zemcs(mxmat))
      endif
      if( ndedx .eq. 1 .or. ndedx .eq. 3 ) then   ! S.Abe 2016/08/08
!$OMP single
        allocate(dbfspace   (kspc,mdbatima))
        allocate(dbispace   (ispc,mdbatima))
        allocate(dbfispace  (   3,mdbatima))
        dbfispace  (:,:) = 0
        if( nedisp .eq. 10 ) then
          allocate(dbfspacene (kspc,mdbatima))
          allocate(dbispacene (ispc,mdbatima))
          allocate(dbfispacene(   3,mdbatima))
          dbfispacene(:,:) = 0
        end if
        if( nspred .eq. 10 ) then
          allocate(dbfspacens (kspc,mdbatima))
          allocate(dbispacens (ispc,mdbatima))
          allocate(dbfispacens(   3,mdbatima))
          dbfispacens(:,:) = 0
        end if
        if( mdbpseud.gt.0 ) then
          allocate(dbpscs_iza(4,mdbpseud))
          allocate(dbpscs_signe(mdbpseud),dbpscs_sigel(mdbpseud))
          allocate(dbpscs_epeak(mdbpseud))
          dbpscs_iza = 0
          dbpscs_epeak = 0.d0
          ldodbps = .true.
          laddbps = .true.
        end if
!$OMP end single
!$OMP barrier
      endif
c----------------------------------------------------------------------
      end subroutine ALLOCATE_MEMBANK
!------------------------------------------------------------------------
      subroutine DEALLOCATE_MEMBANK
      implicit none
      integer nspred,nwsprd,nedisp,itstep,ndedx
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx

      deallocate(siggc,siggd)
      deallocate(isigc,isigd)
      deallocate(siggcn,siggce,siggpn,siggcmx)

      deallocate(isigza)

      deallocate(sigge,siggn)
      deallocate(rnge,erng)
      if( nedisp .ne. 0 )then
        deallocate(avam,avaz,avad)
      endif
      if( nedisp .eq. 10 .or. nspred .eq. 10 .or.
     &     ndedx .eq. 1  .or.  ndedx .eq. 3 )then   ! S.Abe 2016/08/08
        deallocate(rhoi)
      endif
      if( ndedx .le. 3 )then   ! S.Abe 2016/08/08
        deallocate(bari,s1,s2,s3)
        deallocate(eion,prtrng,prtrnglog)
        deallocate(f,g,ro)
      endif
      if( nspred .ne. 0 )then
        deallocate(arg,zsmcs,zemcs)
      endif
      if( ndedx .eq. 1 .or. ndedx .eq. 3 )then   ! S.Abe 2016/08/08
!$OMP barrier
!$OMP single
        deallocate(dbfspace  ,dbispace  ,dbfispace  )
        if( nedisp .eq. 10 ) then
          deallocate(dbfspacene,dbispacene,dbfispacene)
        end if
        if( nspred .eq. 10 ) then
          deallocate(dbfspacens,dbispacens,dbfispacens)
        end if
        if( mdbpseud.gt.0 ) then
          deallocate(dbpscs_iza,dbpscs_signe,dbpscs_sigel,dbpscs_epeak)
        end if
!$OMP end single
      endif
c----------------------------------------------------------------------
      end subroutine DEALLOCATE_MEMBANK
!------------------------------------------------------------------------

      end module MEMBANKMOD
