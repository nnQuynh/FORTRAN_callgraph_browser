************************************************************************
      module udm_Parameter
************************************************************************

      implicit none

      integer         , save :: iudmodel = 0 ! udm_int_num > 0 or udm_part_num > 0

      double precision, save :: udm_Kin
      integer         , save :: udm_kf_incident
      integer         , save :: udm_kf_for_11
      integer         , save :: udm_kf_for_12
      integer         , save :: udm_kf_for_13
      integer         , save :: udm_kf_for_21
      double precision, save :: udm_lifetime
      double precision, save :: udm_mass
      integer         , save :: udm_charge
      logical         , save :: udm_logical
      integer         , save :: udm_integer
      integer         , save :: udm_counter_dklos
      character(200)  , save :: udm_phits_path
      integer         , save :: udm_cell
      double precision, save :: udm_fpl
      logical         , save :: udm_has_step_backed = .false.
!$OMP THREADPRIVATE(udm_Kin)
!$OMP THREADPRIVATE(udm_kf_incident)
!$OMP THREADPRIVATE(udm_kf_for_11)
!$OMP THREADPRIVATE(udm_kf_for_12)
!$OMP THREADPRIVATE(udm_kf_for_13)
!$OMP THREADPRIVATE(udm_kf_for_21)
!$OMP THREADPRIVATE(udm_lifetime)
!$OMP THREADPRIVATE(udm_mass)
!$OMP THREADPRIVATE(udm_charge)
!$OMP THREADPRIVATE(udm_logical)
!$OMP THREADPRIVATE(udm_integer)
!$OMP THREADPRIVATE(udm_counter_dklos)
!$OMP THREADPRIVATE(udm_phits_path)
!$OMP THREADPRIVATE(udm_cell)
!$OMP THREADPRIVATE(udm_fpl)
!$OMP THREADPRIVATE(udm_has_step_backed)


      integer, parameter :: udm_module_nMax = 20
      integer, parameter :: udm_int_param_nMax = 40
      integer, parameter :: udm_part_param_nMax = 40
      integer, parameter :: nStrMax = 200
C     For [ User defined interaction ]
      integer, save      :: udm_int_num = 0
      character(len=99), allocatable,save:: udm_int_name(:)
      double precision, allocatable,save:: udm_bias(:)
      double precision      , allocatable,save:: udm_int_param    (:,:)
      character(len=nStrMax), allocatable,save:: udm_int_param_str(:,:)
C     For [ User defined particle ]
      integer, save :: udm_part_num = 0
      character(len=99), allocatable,save:: udm_part_name(:)
      integer, allocatable,save:: udm_part_kf(:)
      double precision, allocatable,save::       udm_part_param(:,:)
      character(len=nStrMax), allocatable,save:: udm_part_param_str(:,:)

C     For Material info
      integer, parameter :: num_nuclide_max = 200
      integer, save :: num_nuclide
      integer, save :: mat_Z(num_nuclide_max)
      integer, save :: mat_A(num_nuclide_max)
      double precision, save :: mat_ratio(num_nuclide_max)
!$OMP THREADPRIVATE(num_nuclide)
!$OMP THREADPRIVATE(mat_Z)
!$OMP THREADPRIVATE(mat_A)
!$OMP THREADPRIVATE(mat_ratio)

C     For "tot" variables in getflt.f
      double precision, save :: udm_sigt, totudm, totudm_mul
      double precision, save :: totudm_other
      double precision, allocatable,save ::  totudm_(:)
!$OMP THREADPRIVATE(udm_sigt, totudm, totudm_mul)
!$OMP THREADPRIVATE(totudm_other)
!$OMP THREADPRIVATE(totudm_)

C     For final states
      integer, parameter :: nFSmax = 200 ! maximal number of final states
      integer, save :: set_final_state_number
      integer, save :: set_kf(nFSmax)
      integer, save :: set_isomer_level(nFSmax) ! (0: Ground, 1,2: 1st, 2nd isomer)
      double precision, save :: set_Total_Energy_in_MeV(nFSmax)
      double precision, save :: set_Px_in_MeV(nFSmax)
      double precision, save :: set_Py_in_MeV(nFSmax)
      double precision, save :: set_Pz_in_MeV(nFSmax)
      double precision, save :: set_excitation_energy_in_MeV(nFSmax)
      logical         , save :: set_decay_success
!$OMP THREADPRIVATE(set_final_state_number)
!$OMP THREADPRIVATE(set_kf)
!$OMP THREADPRIVATE(set_isomer_level)
!$OMP THREADPRIVATE(set_Total_Energy_in_MeV)
!$OMP THREADPRIVATE(set_Px_in_MeV)
!$OMP THREADPRIVATE(set_Py_in_MeV)
!$OMP THREADPRIVATE(set_Pz_in_MeV)
!$OMP THREADPRIVATE(set_excitation_energy_in_MeV)
!$OMP THREADPRIVATE(set_decay_success)
      double precision, save :: udm_tmp_save = 0d0

      contains

************************************************************************
      subroutine udm_check(action)
************************************************************************
      integer action,i,j
      logical correct
      correct=.true.
C     ------------------------------------------------------------------
      if(action .eq. 1) then
C       [interaction]
        do i=1,udm_int_num
          do j=i+1,udm_int_num
            if(udm_int_name(i) .eq. udm_int_name(j)) then
              correct=.false.
              exit
            endif
          enddo
        enddo
        if(.not. correct) then
          print*,"In [ user defined interaction ], Name is duplicated."
          call abort
        endif
C       [particle]
        do i=1,udm_part_num
          do j=i+1,udm_part_num
            if(udm_part_name(i) .eq. udm_part_name(j)) then
              correct=.false.
              exit
            endif
          enddo
        enddo
        if(.not. correct) then
          print*,"In [ user defined particle ], Name is duplicated."
          call abort
        endif
C     ------------------------------------------------------------------
      elseif(action .eq. 2) then
        do i=1,udm_part_num
          do j=i+1,udm_part_num
            if(udm_part_kf(i) .eq. udm_part_kf(j)) then
              correct=.false.
              exit
            endif
          enddo
        enddo
        if(.not. correct) then
          print*,"In [ user defined particle ], kf is duplicated."
          call abort
        endif
C     ------------------------------------------------------------------
      endif
      end subroutine udm_check

************************************************************************
      subroutine udm_initialize
************************************************************************
      allocate( totudm_( udm_int_num ) )
      udm_counter_dklos=0
      end subroutine udm_initialize

************************************************************************
      subroutine initialize_udm_event_info
************************************************************************
      integer i
      do i=1,nFSmax
        set_final_state_number=0
        set_kf(i)=0
        set_isomer_level(i)=0
        set_Total_Energy_in_MeV(i)=0d0
        set_Px_in_MeV(i)=0d0
        set_Py_in_MeV(i)=0d0
        set_Pz_in_MeV(i)=0d0
        set_excitation_energy_in_MeV(i)=0d0
      enddo
      end subroutine initialize_udm_event_info

************************************************************************
      function get_ityp(kf)
************************************************************************
      integer get_ityp, kf
      get_ityp=11
      if(kf.eq.2212) get_ityp=1
      if(kf.eq.2112) get_ityp=2
      if(kf.eq. 211) get_ityp=3
      if(kf.eq. 111) get_ityp=4
      if(kf.eq.-211) get_ityp=5
      if(kf.eq. -13) get_ityp=6
      if(kf.eq.  13) get_ityp=7
      if(kf.eq. 321) get_ityp=8
      if(kf.eq. 311) get_ityp=9
      if(kf.eq.-321) get_ityp=10
C                    get_ityp=11
      if(kf.eq.  11) get_ityp=12
      if(kf.eq. -11) get_ityp=13
      if(kf.eq.  22) get_ityp=14
      if(kf.eq.1000002) get_ityp=15
      if(kf.eq.1000003) get_ityp=16
      if(kf.eq.2000003) get_ityp=17
      if(kf.eq.2000004) get_ityp=18
      if(kf.gt.2000004) get_ityp=19
      return
      end function

************************************************************************
      function get_proton_neutron_number(kf)
************************************************************************
      integer get_proton_neutron_number(2)
      integer kf,Z,N
      get_proton_neutron_number(1)=0
      get_proton_neutron_number(2)=0
      if(kf.eq.2212) get_proton_neutron_number(1)=1
      if(kf.eq.2112) get_proton_neutron_number(2)=1
      if(kf.ge.1000002) then
        Z = kf / 1000000
        N = kf - kf / 1000000 * 1000000 - Z
        get_proton_neutron_number(1)=Z
        get_proton_neutron_number(2)=N
      endif
      return
      end function

************************************************************************
      integer function get_cell(i)
************************************************************************
      integer i,iblz1,iblz2
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
C     iblz1 : cell id at (x,y,z)                                    *
C     iblz2 : cell id after crossing                                *
      if(i==1) then; get_cell=iblz1; return; endif
      if(i==2) then; get_cell=iblz2; return; endif
      return
      end function

************************************************************************
      function udm_interp(X,ndata,Xdata,Ydata,type)
************************************************************************
      double precision udm_interp
      double precision X,Y
      integer ndata
      double precision,intent(in) :: Xdata(:)
      double precision,intent(in) :: Ydata(:)
      character(len=*) :: type
C     ------------------------------------------------------------------
      integer i,iX1,iX2
      double precision X1,X2,Y1,Y2
C     ------------------------------------------------------------------
      if(X .lt. Xdata(1)) then
        udm_interp=0d0
        return
      endif
C     -------------------
      if(X .gt. Xdata(ndata)) then
        udm_interp=0d0
        return
      endif
C     ------------------------------------------------------------------
      iX1=0
      do i = 2, ndata
        if(X .lt. Xdata(i)) then
          iX1=i-1
          exit
        endif
      enddo
C     ------------------------------------------------------------------
      if(0 .lt. iX1 .and. iX1 .lt. ndata) then
        iX2=iX1+1
        X1=Xdata(iX1); Y1=Ydata(iX1);
        X2=Xdata(iX2); Y2=Ydata(iX2);
C       --------------------------------
        if    (type .eq. "lin-lin") then
          Y=(Y2-Y1)/(X2-X1)*(X-X1)+Y1
C       --------------------------------
        elseif(type .eq. "lin-log") then
          if(Y1<=0d0 .or. Y2<=0d0) then
            Y=0d0
          else
            Y=exp((log(Y2)-log(Y1))/(X2-X1)*(X-X1)+log(Y1))
          endif
C       --------------------------------
        elseif(type .eq. "log-lin") then
          if(X1<=0d0 .or. X2<=0d0) then
            Y=0d0
          else
            Y=(Y2-Y1)/(log(X2)-log(X1))*(log(X)-log(X1))+Y1
          endif
C       --------------------------------
        elseif(type .eq. "log-log") then
          if(X1<=0d0 .or. X2<=0d0 .or. Y1<=0d0 .or. Y2<=0d0) then
            Y=0d0
          else
            Y=exp((log(Y2)-log(Y1))/(log(X2)-log(X1))*(log(X)-log(X1))
     &       +log(Y1))
          endif
C       --------------------------------
        else
          print*,"Caution: udm_interp in type"
          Y=(Y2-Y1)/(X2-X1)*(X-X1)+Y1
C       --------------------------------
        endif
        udm_interp=Y
        return
      else
        print*,"caution: udm_interp"
        udm_interp=0
        return
      endif
C     ------------------------------------------------------------------
      udm_interp=0
      return
      end function

************************************************************************
      end module udm_Parameter
************************************************************************



