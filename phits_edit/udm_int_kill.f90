!=======================================================================
module udm_int_kill
!=======================================================================

use udm_Parameter
use udm_Utility
implicit none
private ! Functions and variables are set to private by default.
public :: caller ! The 'caller' subroutine should be public.

!-----------------------------------------------------------------------
! Default variables
!-----------------------------------------------------------------------
character(len=99), parameter :: Name = "kill" ! This 'Name' is used in input files.
! Parameters entered in the input file
double precision  , allocatable, save :: Parameters(:)     ! real number format
character(len=200), allocatable, save :: Parameters_str(:) ! string format
!-----------------------------------------------------------------------
! User variables
!-----------------------------------------------------------------------
integer          , allocatable, save :: kill_kf(:)
integer          , allocatable, save :: kill_cell(:)
double precision , allocatable, save :: kill_Emin(:)
integer, save :: n_kill ! Number of kill set
contains


!=======================================================================
subroutine initialize
! This subroutine is called only once at the beginning of the calculation.
!=======================================================================
implicit none
logical :: exists
logical :: isDirectInputMode
integer :: i,kf,cell
double precision :: Emin

isDirectInputMode = nint(Parameters(2))/=0

! ------------
! count n_kill
! ------------
if(isDirectInputMode) then
  do i=1,13
    cell=nint(Parameters(3*(i-1)+2))
    if(cell==0) then
      exit
    else
      n_kill=i
    endif
  enddo
! ------------
else
  inquire(file=trim(Parameters_str(1)), exist=exists)
  if(.not. exists) then
    print*,"Error: kill. File not found:", trim(Parameters_str(1))
    stop
  endif
  open(2201,file=trim(Parameters_str(1)),status="old",err=999)
  n_kill=0
  do
888 continue
    read(2201,*,err=888,end=777) kf,cell,Emin
    n_kill=n_kill+1
  enddo
777 continue
endif
! ------------

allocate( kill_kf  (n_kill) )
allocate( kill_cell(n_kill) )
allocate( kill_Emin(n_kill) )

! ------------
! fill parameters
! ------------
if(isDirectInputMode) then
  do i=1,n_kill
    kill_kf  (i)=nint(Parameters(3*(i-1)+1))
    kill_cell(i)=nint(Parameters(3*(i-1)+2))
    kill_Emin(i)=     Parameters(3*(i-1)+3)
  enddo
! ------------
else
  rewind(2201)
  i=0
  do
666 continue
    read(2201,*,err=666) kf,cell,Emin
    i=i+1
    kill_kf  (i)=kf
    kill_cell(i)=cell
    kill_Emin(i)=Emin
    if(i==n_kill) exit
  enddo
  close(2201)
endif
! ------------


! ------------
! check by print 
! ------------
print*,"###########################################"
print*,"### initialize: kill particle condition ###"
print*,"#        kf        cell   emin[MeV]"
do i=1,n_kill
  print*,kill_kf(i),kill_cell(i),kill_Emin(i)
enddo
print*,"(kf=0 means all particles)"
print*,"###########################################"
! ---------------

return

999 continue
print*,"Error: kill."
return




end subroutine initialize


!=======================================================================
double precision function Xsec_per_atom(Kin,Z,A)
! Integrated cross section per an atom.
! Unit: barn (10^-24 cm2)
implicit none
double precision Kin ! Kinetic energy of incident particle [MeV]
           integer Z ! Atomic number of target atom
           integer A ! Mass   number of target atom
!-----------------------------------------------------------------------
! [Variables available in this function]
! udm_kf_incident : The kf-codes (particle IDs) of the incident particles.
!=======================================================================
integer i
do i=1,n_kill
  if(udm_kf_incident == kill_kf(i) .or. kill_kf(i) == 0 ) then 
  if(get_cell(1)     == kill_cell(i)) then 
  if(Kin             <  kill_Emin(i)) then 
    Xsec_per_atom=1e+20
    return
  endif
  endif
  endif
enddo

Xsec_per_atom = 0d0
return

end function



!=======================================================================
subroutine generate_final_state
!=======================================================================
! Subroutine to determine final state information (4-momenta, etc.).
!-----------------------------------------------------------------------
! [Variables available in this subroutine]
! udm_kf_incident : The kf-codes (particle IDs) of the incident particles.
! udm_Kin         : Kinetic Energy of the incident particle [MeV]
!-----------------------------------------------------------------------

call initialize_udm_event_info
! ------------------------------------
! Set number of final states to be recorded in event history.
set_final_state_number = 1
! ------------------------------------
! Set 4-momentum of X
! set_kf                 (1) = udm_kf_incident
! set_Total_Energy_in_MeV(1) = get_mass(udm_kf_incident)
set_kf                 (1) = 987654
set_Total_Energy_in_MeV(1) = get_mass(set_kf(1))
set_Px_in_MeV          (1) = 0d0
set_Py_in_MeV          (1) = 0d0
set_Pz_in_MeV          (1) = 0d0
! ------------------------------------
! The final states are recorded.
call fill_final_state
return

end subroutine generate_final_state

























!=======================================================================
!     DO NOT CHANGE BELOW
!=======================================================================
subroutine caller(action,index)
!=======================================================================
integer action,index
if( action .eq.  2 ) call check_match_name(index)
if(udm_int_name(index) .ne. Name) return
if( action .eq.  1 ) call fill_Parameters(index)
if( action .eq.  3 ) call initialize
if( action .eq. 11 ) call calc_averaged_Xsec
if( action .eq. 21 ) call generate_final_state
end subroutine caller

!=======================================================================
subroutine check_match_name(index)
!=======================================================================
integer index
udm_logical = udm_logical .or. (udm_int_name(index) .eq. Name)
end subroutine check_match_name

!=======================================================================
subroutine fill_Parameters(index)
!=======================================================================
integer index
integer i
if(allocated(Parameters)) then
  print*,"*** Error: ",trim(Name)
  print*,"*** Calling this module multiple times is not supported."
  stop
endif
allocate( Parameters    ( udm_int_param_nMax ) )
allocate( Parameters_str( udm_int_param_nMax ) )
do i=1,udm_int_param_nMax
  Parameters    (i)=     udm_int_param    (index,i)
  Parameters_str(i)=trim(udm_int_param_str(index,i))
enddo
print*,"Set [ User defined interaction ]: ",trim(Name)
! --------------------------------------------------
! --------------------------------------------------
udm_bias(index)=1d0 ! special for this module
! --------------------------------------------------
! --------------------------------------------------
end subroutine fill_Parameters

!=======================================================================
subroutine calc_averaged_Xsec
!     Averaged cross section for mixed atoms
!=======================================================================
integer i
do i = 1, num_nuclide
  udm_sigt = udm_sigt + mat_ratio(i)*Xsec_per_atom( udm_Kin, mat_Z(i), mat_A(i) )
enddo
end subroutine calc_averaged_Xsec

!=======================================================================
function get_hit_nuclide_Z_A(Kin)
!=======================================================================
integer i, get_hit_nuclide_Z_A(2)
double precision Kin, randTMP, XsecSum, tmp
!     ------------------------------------------------------------------
XsecSum=0d0
do i = 1, num_nuclide
  XsecSum = XsecSum + mat_ratio(i)*Xsec_per_atom( Kin, mat_Z(i), mat_A(i) )
enddo
!     ------------------------------------------------------------------
randTMP=get_random_0to1()
tmp=0d0
do i = 1, num_nuclide
  tmp = tmp + ( mat_ratio(i)*Xsec_per_atom( Kin, mat_Z(i), mat_A(i) ) ) / XsecSum
  if(randTMP .lt. tmp) then ! nuclide=i is accepted.
    get_hit_nuclide_Z_A(1)=mat_Z(i) ! Z of hit nuclide
    get_hit_nuclide_Z_A(2)=mat_A(i) ! A of hit nuclide
    return
  endif
enddo
!     ------------------------------------------------------------------
print*,"[CAUTION] Something wrong. hit_nuclide was not choosed."
get_hit_nuclide_Z_A(1)=mat_Z(1)
get_hit_nuclide_Z_A(2)=mat_A(1)
return
end function

!=======================================================================
end module udm_int_kill
!=======================================================================






