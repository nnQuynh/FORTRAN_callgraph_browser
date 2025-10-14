!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!module moddas
      module moddas

      integer, parameter :: MAX_NUM_LARGE = 200000
!...
! maximum size of chrg ( >ichmx ) sub.setgg  @ ggs00.f
!                               , sub.setcg  @ marscg.f
!                               , sub.setpar @ read00.f
!                               , sub.setpag @ read00.f
!                               , sub.echoi  @ read00.f
!                               , sub.read01 @ read01.f
!                               , sub.region @ read02.f
!                               , sub.cgview @ read02.f
!                               , sub.gcell  @ read03.f
      integer, parameter :: MAX_NUM_CHRG = 200000

!...
! maximum size of slib ( >itls ) sub.tmultip @ tallsm1.f
      integer, parameter :: MAX_NUM_SLIB = 999

!...
! maximum size of mtrg ( >ipmax )
!    sub.tregion, tregion2, tregion3, tregion4 @ tallsm1.f
      integer, parameter :: MAX_NUM_MTRG = MAX_NUM_LARGE

! maximum size of idas_itreg ( >mtrn )
!    sub.settal @ tallsm1.f, sub.t* @ tallsm2.f
      integer, parameter :: MAX_NUM_ITREG = MAX_NUM_LARGE

! maximum size of idas_itmeg ( >mtbn ) sub.tdshow @ tallsm2.f
      integer, parameter :: MAX_NUM_ITMEG = MAX_NUM_LARGE

! maximum size of nrst ( >nnm ) sub.echrg2, echrg3 @ tallsm3.f
      integer, parameter :: MAX_NUM_NRST = MAX_NUM_LARGE

! maximum size of idas_itrcr ( >jdsm ) sub.settal @ tallsm1.f
      integer, parameter :: MAX_NUM_ITRCR = MAX_NUM_LARGE

! maximum size of idas_ntcr ( >jdsm ) sub.tcross @ tallsm2.f
      integer, parameter :: MAX_NUM_NTCR = MAX_NUM_LARGE

! maximum size of das_ktcr ( >ntcn ) sub.tcross @ tallsm2.f
      integer, parameter :: MAX_NUM_KTCR = MAX_NUM_LARGE

! maximum size of idas_region_temporary ( > ) sub.tdepreg0 @ result.f
      integer, parameter :: MAX_NUM_REGION_TEMPORARY = MAX_NUM_LARGE

!...
! maximum size of idas_nsrn ( >jdsm ) sub.setpag @ read00.f
      integer, parameter :: MAX_NUM_NSRN = MAX_NUM_LARGE

! maximum size of idas_inimt ( >jdsm ) sub.setpag @ read00.f
      integer, parameter :: MAX_NUM_INIMT = MAX_NUM_LARGE

! maximum size of idas_inwwt ( >jdsm ) sub.setpag @ read00.f
      integer, parameter :: MAX_NUM_INWWT = MAX_NUM_LARGE

! maximum size of idas_inflt ( >jdsm ) sub.setpag @ read00.f
      integer, parameter :: MAX_NUM_INFLT = MAX_NUM_LARGE

! maximum size of idas_inrlt ( >jdsm ) sub.setpag @ read00.f
      integer, parameter :: MAX_NUM_INRLT = MAX_NUM_LARGE

! maximum size of idas_ipgrt ( >jdsm ) sub.setpag @ read00.f
      integer, parameter :: MAX_NUM_IPGRT = MAX_NUM_LARGE

! maximum size of idas_incrt ( >jdsm ) sub.setpag @ read00.f
      integer, parameter :: MAX_NUM_INCRT = MAX_NUM_LARGE

! maximum size of idas_ingrt ( >jdsm ) sub.setpag @ read00.f
      integer, parameter :: MAX_NUM_INGRT = MAX_NUM_LARGE

! maximum size of idas_inert ( >jdsm ) sub.setpag @ read00.f
      integer, parameter :: MAX_NUM_INERT = MAX_NUM_LARGE

! maximum size of idas_isgrt ( >jdsm ) sub.setpag @ read00.f
      integer, parameter :: MAX_NUM_ISGRT = MAX_NUM_LARGE

! maximum size of idas_intmt ( >jdsm ) sub.setpag @ read00.f
      integer, parameter :: MAX_NUM_INTMT = MAX_NUM_LARGE

!...
! maximum size of idas_nsrc ( >jdsm ) sub.sours @ read02.f
      integer, parameter :: MAX_NUM_NSRC = MAX_NUM_LARGE

! maximum size of idas_inimc ( >jdsm ) sub.impot @ read02.f
      integer, parameter :: MAX_NUM_INIMC = MAX_NUM_LARGE

! maximum size of idas_inwwc ( >jdsm ) sub.wwind @ read02.f
      integer, parameter :: MAX_NUM_INWWC = MAX_NUM_LARGE

! maximum size of idas_inwbc ( >jdsm ) sub.wwbias @ read02.f
      integer, parameter :: MAX_NUM_INWBC = MAX_NUM_LARGE

! maximum size of idas_inflc ( >jdsm ) sub.foccl @ read02.f
      integer, parameter :: MAX_NUM_INFLC = MAX_NUM_LARGE

! maximum size of idas_inrlc ( >jdsm ) sub.repcl @ read02.f
      integer, parameter :: MAX_NUM_INRLC = MAX_NUM_LARGE

! maximum size of idas_ipgrc ( >jdsm ) sub.split @ read02.f
      integer, parameter :: MAX_NUM_IPGRC = MAX_NUM_LARGE

! maximum size of idas_incrc ( >jdsm ) sub.counts @ read02.f
      integer, parameter :: MAX_NUM_INCRC = MAX_NUM_LARGE

! maximum size of idas_ingrc ( >jdsm ) sub.mgnet @ read02.f
      integer, parameter :: MAX_NUM_INGRC = MAX_NUM_LARGE

! maximum size of idas_inerc ( >jdsm ) sub.elmgf @ read02.f
      integer, parameter :: MAX_NUM_INERC = MAX_NUM_LARGE

! maximum size of idas_isgrc ( >jdsm ) sub.supmir @ read02.f
      integer, parameter :: MAX_NUM_ISGRC = MAX_NUM_LARGE

! maximum size of idas_intmc ( >jdsm ) sub.timers @ read02.f
      integer, parameter :: MAX_NUM_INTMC = MAX_NUM_LARGE

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
!contains
      contains

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_initialize()

      use moddas_character
      use moddas_fragdata
      use moddas_material
      use moddas_mesh
      use moddas_multiplier
      use moddas_region
      use moddas_region_mtrg
      use moddas_repeated_collisions
      use moddas_source
      use moddas_tally
      use moddas_variance_reduction

      call moddas_character_initialize()
      call moddas_fragdata_initialize()
      call moddas_material_initialize()
      call moddas_mesh_initialize()
      call moddas_multiplier_initialize()
      call moddas_region_initialize()
      call moddas_region_mtrg_initialize()
      call moddas_repeated_collisions_initialize()
      call moddas_source_initialize()
      call moddas_tally_initialize()
      call moddas_variance_reduction_initialize()
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_reallocate_dbl(
     &           naddress, icurrent, num_data_add, iaddress, data)

      include 'err.inc'

      integer, intent(in) :: naddress
      integer, intent(in) :: icurrent
      integer, intent(in) :: num_data_add
      integer, intent(inout) :: iaddress(naddress)
      double precision, intent(inout), allocatable :: data(:)

      integer :: len_realloc, len_data, isave
      double precision, allocatable :: work(:)

      if( num_data_add > 0 ) then
! additional length of array
         len_realloc = iaddress(icurrent) + num_data_add
         if( icurrent < naddress ) then
            iaddress(icurrent+1:) = len_realloc
         end if
! save current data, deallocate array
         isave = 0
         len_data = 0
         if( allocated(data) ) then
            isave = 1
            len_data = iaddress(icurrent)
            allocate( work(len_data) )
            work(1:len_data) = data(1:len_data)
            deallocate( data )
         end if
! allocate array
         allocate( data(len_realloc) )
         data(:) = 0.0d0
! restore current data
         if( isave == 1 ) then
            data(1:len_data) = work(1:len_data)
            deallocate( work )
         end if
      endif
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_reallocate_dbl2(
     &           naddress, icurrent, nd1, num_data_add, iaddress, data)

      include 'err.inc'

      integer, intent(in) :: naddress
      integer, intent(in) :: icurrent
      integer, intent(in) :: nd1
      integer, intent(in) :: num_data_add
      integer, intent(inout) :: iaddress(naddress)
      double precision, intent(inout), allocatable :: data(:,:)

      integer :: len_realloc, len_data, isave
      double precision, allocatable :: work(:,:)

      if( num_data_add > 0 ) then
! additional length of array
         len_realloc = iaddress(icurrent) + num_data_add
         if( icurrent < naddress ) then
            iaddress(icurrent+1:) = len_realloc
         end if
! save current data, deallocate array
         isave = 0
         len_data = 0
         if( allocated(data) ) then
            isave = 1
            len_data = iaddress(icurrent)
            allocate( work(nd1,len_data) )
            work(:,1:len_data) = data(:,1:len_data)
            deallocate( data )
         end if
! allocate array
         allocate( data(nd1,len_realloc) )
         data(:,:) = 0.0d0
! restore current data
         if( isave == 1 ) then
            data(:,1:len_data) = work(:,1:len_data)
            deallocate( work )
         end if
      endif
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_reallocate_int(
     &           naddress, icurrent, num_data_add, iaddress, idata)

      include 'err.inc'

      integer, intent(in) :: naddress
      integer, intent(in) :: icurrent
      integer, intent(in) :: num_data_add
      integer, intent(inout) :: iaddress(naddress)
      integer, intent(inout), allocatable :: idata(:)

      integer :: len_realloc, len_data, isave
      integer, allocatable :: iwork(:)

      if( num_data_add > 0 ) then
! additional length of array
         len_realloc = iaddress(icurrent) + num_data_add
         if( icurrent < naddress ) then
            iaddress(icurrent+1:) = len_realloc
         end if
! save current data, deallocate array
         isave = 0
         len_data = 0
         if( allocated(idata) ) then
            isave = 1
            len_data = iaddress(icurrent)
            allocate( iwork(len_data) )
            iwork(1:len_data) = idata(1:len_data)
            deallocate( idata )
         end if
! allocate array
         allocate( idata(len_realloc) )
         idata(:) = 0
! restore current data
         if( isave == 1 ) then
            idata(1:len_data) = iwork(1:len_data)
            deallocate( iwork )
         end if
      endif
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_reallocate_int2(
     &   naddress, icurrent, lb1, ub1, num_data_add, iaddress, idata)

      include 'err.inc'

      integer, intent(in) :: naddress
      integer, intent(in) :: icurrent
      integer, intent(in) :: lb1
      integer, intent(in) :: ub1
      integer, intent(in) :: num_data_add
      integer, intent(inout) :: iaddress(naddress)
      integer, intent(inout), allocatable :: idata(:,:)

      integer :: len_realloc, len_data, isave
      integer, allocatable :: iwork(:,:)

      if( num_data_add > 0 ) then
! additional length of array
         len_realloc = iaddress(icurrent) + num_data_add
         if( icurrent < naddress ) then
            iaddress(icurrent+1:) = len_realloc
         end if
! save current data, deallocate array
         isave = 0
         len_data = 0
         if( allocated(idata) ) then
            isave = 1
            len_data = iaddress(icurrent)
            allocate( iwork(lb1:ub1,len_data) )
            iwork(:,1:len_data) = idata(:,1:len_data)
            deallocate( idata )
         end if
! allocate array
         allocate( idata(lb1:ub1,len_realloc) )
         idata(:,:) = 0
! restore current data
         if( isave == 1 ) then
            idata(:,1:len_data) = iwork(:,1:len_data)
            deallocate( iwork )
         end if
      endif
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_allocate_dbl(num_data, data)

      include 'err.inc'

      integer, intent(in) :: num_data
      double precision, intent(inout), allocatable :: data(:)

! allocate array
      if( num_data > 0 ) then
         if( allocated(data) ) then
            deallocate( data )
         end if

         allocate( data(num_data) )
         data(:) = 0.0d0
      endif
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_allocate_dbl2(nd1, lb2, ub2, data)

      include 'err.inc'

      integer, intent(in) :: nd1
      integer, intent(in) :: lb2, ub2
      double precision, intent(inout), allocatable :: data(:,:)

      integer :: num_data

! allocate array
      num_data = nd1*(ub2 - lb2 + 1)
      if( num_data > 0 ) then
         if( allocated(data) ) then
            deallocate( data )
         end if

         allocate( data(nd1, lb2:ub2) )
         data(:,:) = 0.0d0
      endif
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_allocate_dbl3(nd1, nd2, nd3, data)

      include 'err.inc'

      integer, intent(in) :: nd1
      integer, intent(in) :: nd2
      integer, intent(in) :: nd3
      double precision, intent(inout), allocatable :: data(:,:,:)

      integer :: num_data

! allocate array
cFURUTA20240509      num_data = nd1*nd2*nd3
cFURUTA20240509      if( num_data > 0 ) then
      if( nd1 > 0 .and. nd2 > 0 .and. nd3> 0 ) then
         if( allocated(data) ) then
            deallocate( data )
         end if

         allocate( data(nd1, nd2, nd3) )
         data(:,:,:) = 0.0d0
      endif
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_allocate_int(num_data, idata)

      include 'err.inc'

      integer, intent(in) :: num_data
      integer, intent(inout), allocatable :: idata(:)

! allocate array
      if( num_data > 0 ) then
         if( allocated(idata) ) then
            deallocate( idata )
         end if

         allocate( idata(num_data) )
         idata(:) = 0
      endif
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_allocate_int2(lb1, ub1, nd2, idata)

      include 'err.inc'

      integer, intent(in) :: lb1
      integer, intent(in) :: ub1
      integer, intent(in) :: nd2
      integer, intent(inout), allocatable :: idata(:,:)

      integer :: num_data

! allocate array
      num_data = (ub1 - lb1 + 1)*nd2
      if( num_data > 0 ) then
         if( allocated(idata) ) then
            deallocate( idata )
         end if

         allocate( idata(lb1:ub1, nd2) )
         idata(:,:) = 0
      endif
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_reduce_dbl(
     &           naddress, icurrent, num_data_add, iaddress, data)

      include 'err.inc'

      integer, intent(in) :: naddress
      integer, intent(in) :: icurrent
      integer, intent(in) :: num_data_add
      integer, intent(inout) :: iaddress(naddress)
      double precision, intent(inout), allocatable :: data(:)

      integer :: len_realloc
      double precision, allocatable :: work(:)

      if( num_data_add > 0 ) then
! resize length of array
         len_realloc = iaddress(icurrent) + num_data_add
         if( icurrent < naddress ) then
            iaddress(icurrent+1:) = len_realloc
         end if
! save current data, deallocate array
         allocate( work(len_realloc) )
         work(1:len_realloc) = data(1:len_realloc)
         deallocate( data )
! allocate array
         allocate( data(len_realloc) )
         data(1:len_realloc) = work(1:len_realloc)
         deallocate( work )
      endif
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_reduce_int(
     &           naddress, icurrent, num_data_add, iaddress, idata)

      include 'err.inc'

      integer, intent(in) :: naddress
      integer, intent(in) :: icurrent
      integer, intent(in) :: num_data_add
      integer, intent(inout) :: iaddress(naddress)
      integer, intent(inout), allocatable :: idata(:)

      integer :: len_realloc
      integer, allocatable :: iwork(:)

      if( num_data_add > 0 ) then
! resize length of array
         len_realloc = iaddress(icurrent) + num_data_add
         if( icurrent < naddress ) then
            iaddress(icurrent+1:) = len_realloc
         end if
! save current data, deallocate array
         allocate( iwork(len_realloc) )
         iwork(1:len_realloc) = idata(1:len_realloc)
         deallocate( idata )
! allocate array
         allocate( idata(len_realloc) )
         idata(1:len_realloc) = iwork(1:len_realloc)
         deallocate( iwork )
      endif
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_deallocate_dbl(data)

      double precision, intent(inout), allocatable :: data(:)

      if( allocated(data) ) then
         deallocate( data )
      endif
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_deallocate_int(idata)

      integer, intent(inout), allocatable :: idata(:)

      if( allocated(idata) ) then
         deallocate( idata )
      endif
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_allocate_cha(num_data, cdata)

      include 'err.inc'

      integer, intent(in) :: num_data
      character(len=:), intent(inout), allocatable :: cdata

! allocate array
      if( num_data > 0 ) then
         if( allocated(cdata) ) then
            deallocate( cdata )
         end if

         allocate( character(len=num_data)::cdata )

      endif
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      subroutine moddas_deallocate_cha(cdata)

      character(len=:), intent(inout), allocatable :: cdata

      if( allocated(cdata) ) then
         deallocate( cdata )
      endif
      end subroutine

!...+....1....+....2....+....3....+....4....+....5....+....6....+....7....+....8
      end module
