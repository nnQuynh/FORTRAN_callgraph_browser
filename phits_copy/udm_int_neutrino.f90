!=======================================================================
module udm_int_neutrino
!=======================================================================

use udm_Parameter
use udm_Utility

implicit none
private ! Functions and variables are set to private by default.
public :: caller ! The 'caller' subroutine should be public.

!-----------------------------------------------------------------------
! Default variables
!-----------------------------------------------------------------------
character(len=99), parameter :: Name = "neutrino" ! This 'Name' is used in input files.
integer, parameter :: num_initial = 6 ! The number of incident particles causing this interaction.
integer, save :: kf_initial(num_initial) = (/ 12,14,16,-12,-14,-16 /) ! The kf-codes (particle IDs) of the incident particles causing this interaction.
!-----------------------------------------------------------------------
! Parameters entered in the input file
double precision  , allocatable, save :: Parameters(:)     ! real number format
character(len=200), allocatable, save :: Parameters_str(:) ! string format
!-----------------------------------------------------------------------
! User variables
!-----------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The energy unit in the file is GeV !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

logical, save :: load_failed
integer, parameter :: unitNum = 2210 ! Unit number used for reading and writing files in this module.

! *** Base
integer                      , save :: nTar    ! Number of targets. From target.dat.
integer         , allocatable, save :: Zs(:)   ! Z array (nTar). From target.dat.
integer         , allocatable, save :: As(:)   ! A array (nTar). From target.dat.
integer                      , save :: nNeu    ! Number of neutrinos. From num_initial.
integer         , allocatable, save :: Neus(:) ! Neutrino array (nNeu). From kf_initial.
integer                      , save :: nEventsPerFile
character(len=999)           , save :: folder

! *** Total cross section
integer                      , save :: nEMax1     ! Maximum number of energy elements.
integer         , allocatable, save :: nE1(:,:)   ! Number of energy elements for each target and neutrino (nTar, nNeu).
double precision, allocatable, save :: Es1(:,:,:) ! E array of the cross-section file (nTar, nNeu, nE1).
double precision, allocatable, save :: XSs(:,:,:) ! Cross-section of the cross-section file (nTar, nNeu, nE1).
! nE1 is obtained by scanning each cross-section file

! *** Final states
integer                      , save :: nEMax2         ! Maximum number of energy elements.
integer         , allocatable, save :: nE2(:,:)       ! Number of energy elements for each target and neutrino (nTar, nNeu).
double precision, allocatable, save :: Es2(:,:,:)     ! E array of the event file (nTar, nNeu, nE2).
character(len=20),allocatable, save :: Es2str(:,:,:)  ! Energy string.
integer         , allocatable, save :: nEvents(:,:,:) ! Number of events in the event file (nTar, nNeu, nE2).
! nE2 is obtained by scanning n_events_list.dat

! *** Option
integer, allocatable, save :: kfAccepted(:) ! Array of kf codes for particles to be included in the final state.
logical, save :: print_detail = .false.

contains

!=======================================================================
subroutine print_comments
!=======================================================================
if(.not. load_failed) then
  if(index(trim(folder), "genie") > 0) then
    print*,"****************************************"
    print*,"This database was created using GENIE."
    print*,"----------------------------------------"
    print*,"The GENIE Neutrino Monte Carlo Generator"
    print*,"Nucl.Instrum.Meth.A 614 (2010) 87-104."
    print*,"http://www.genie-mc.org/"
    print*,"https://github.com/GENIE-MC"
    print*,"****************************************"
  endif
  print*
  print*,"Due to the small neutrino cross section, the use of [Forced Collisions] is recommended."
  print*
endif
end subroutine print_comments

!=======================================================================
integer function nLines(path)
!=======================================================================
character(len=*) :: path
integer ios
character(len=999) :: Line

! Set the element number of the file
nLines = 0

! Open the file
open(unit=unitNum, file=trim(path), status='old', action='read', iostat=ios)

! Check if the file was opened successfully
if (ios /= 0) then
  print*, 'Error opening file:', path
endif

! Read until the end of the file and count the number of lines
do
  read(unitNum, '(A)', iostat=ios) Line
  if (ios /= 0) exit
  ! Check if the line is blank
  if (trim(adjustl(Line)) /= '') then
    nLines = nLines + 1
  endif
end do

! Close the file
close(unitNum)

return
end

!=======================================================================
character*3 function str3(num)
!=======================================================================
integer num
write(str3, '(I3.3)') num
return
end

!=======================================================================
logical function exists(path)
!=======================================================================
character(len=*), intent(in) :: path
inquire(file=path, exist=exists)
return
end

!=======================================================================
subroutine warning
!=======================================================================
print*,"*** Neutrino interaction is ignored"
print*,"*** Please write the path of the data file in the Parameters section as follows:"
print*,"------------------------------------------"
print*,"[ User Defined Interaction ]"
print*,"$ Name      Bias  Parameters"
print*,"  neutrino  1     /Users/you/GENIE_R-3_02_02_G18_02a_00_000"
print*,"------------------------------------------"
print*,"A test database may be downloadable from the following link:"
print*,"https://rcwww.kek.jp/research/shield/sakaki/phits/temp.html"
print*,"The database will likely be shared on the PHITS Wiki server."
load_failed=.true.
end subroutine warning

!=======================================================================
subroutine initialize
! This subroutine is called only once at the beginning of the calculation.
!=======================================================================
character(len=999) :: path
integer iTar,iNeu,iE,i
integer nE,nEMax1,nEMax2
logical logitmp
character(len=100) :: line
integer pos
integer nkfA

folder=trim(Parameters_str(1))

if(trim(folder)=="") then
  call warning
  return
endif

! ----------
! Check if the folder or file exists

! ----------
path=trim(folder)//"/target.dat"
if(.not. exists(path)) then
  print*,"*** Failed to load: ",trim(path)
  call warning
  return
endif

! ----------
path=trim(folder)//"/events/nEventsPerFile.dat"
if(.not. exists(path)) then
  print*,"*** Failed to load: ",trim(path)
  call warning
  return
endif

! ----------
! print*,"Set: ",trim(folder)
write(*,'(a)') "Set: "//trim(folder)

! ----------
! *** nTar
path=trim(folder)//"/target.dat"
nTar = nLines(path)
allocate(Zs(nTar))
allocate(As(nTar))

! *** The atomic numbers (Zs) and mass numbers (As) where the data is stored.
open(unitNum,file=trim(path),status="old")
do iTar=1,nTar
  read(unitNum,*) Zs(iTar), As(iTar)
enddo
close(unitNum)

! *** nNeu, Neus
nNeu=num_initial
allocate(Neus(nNeu))
do iNeu=1,nNeu
  Neus(iNeu) = kf_initial(iNeu)
enddo

allocate(nE1(nTar,nNeu))
allocate(nE2(nTar,nNeu))

! *** nEventsPerFile
path=trim(folder)//"/events/nEventsPerFile.dat"
open(unitNum,file=trim(path),status="old")
read(unitNum,*) nEventsPerFile
close(unitNum)

! --------------------------------------------------
! *** nEMax1
nEMax1=0
do iTar=1,nTar
do iNeu=1,nNeu
  path=trim(folder)//"/cross-section/"//str3(Zs(iTar))//str3(As(iTar))//"_"//trim(udm_int2char(Neus(iNeu)))//".dat"
  nE=nLines(path)
  nE1(iTar,iNeu)=nE
  if(nE>nEMax1) nEMax1=nE
enddo
enddo
allocate(Es1(nTar,nNeu,nEMax1))
allocate(XSs(nTar,nNeu,nEMax1))

! Get cross-section data
do iTar=1,nTar
write(*,'(a)') "-> /cross-section/"//str3(Zs(iTar))//str3(As(iTar))
do iNeu=1,nNeu
  path=trim(folder)//"/cross-section/"//str3(Zs(iTar))//str3(As(iTar))//"_"//trim(udm_int2char(Neus(iNeu)))//".dat"
  open(unitNum,file=trim(path),status="old")
  nE=nE1(iTar,iNeu)
  do iE=1,nE
    read(unitNum,*) Es1(iTar,iNeu,iE), XSs(iTar,iNeu,iE)
  enddo
  close(unitNum)
  if(print_detail) &
& write(*,'(I6,a,1PE12.5E2,a,1PE12.5E2)') Neus(iNeu),": E =",Es1(iTar,iNeu,1)," to ",Es1(iTar,iNeu,nE)
enddo
enddo

! --------------------------------------------------
! *** nEMax2
nEMax2=0
do iTar=1,nTar
do iNeu=1,nNeu
  path=trim(folder)//"/events/"//str3(Zs(iTar))//str3(As(iTar))//"/"//trim(udm_int2char(Neus(iNeu)))//"/n_events_list.dat"
  open(unitNum,file=trim(path),status="old")
  read(unitNum,*) nE
  nE2(iTar,iNeu)=nE
  if(nE>nEMax2) nEMax2=nE
  close(unitNum)
enddo
enddo
allocate(Es2    (nTar,nNeu,nEMax2))
allocate(Es2str (nTar,nNeu,nEMax2))
allocate(nEvents(nTar,nNeu,nEMax2))

! Get the number of events in each file
do iTar=1,nTar
write(*,'(a)') "-> /events/"//str3(Zs(iTar))//str3(As(iTar))
do iNeu=1,nNeu
  path=trim(folder)//"/events/"//str3(Zs(iTar))//str3(As(iTar))//"/"//trim(udm_int2char(Neus(iNeu)))//"/n_events_list.dat"
  open(unitNum,file=trim(path),status="old")
  read(unitNum,*) nE
  do iE=1,nE
    read(unitNum,*) Es2str(iTar,iNeu,iE), nEvents(iTar,iNeu,iE)
    read(Es2str(iTar,iNeu,iE),*) Es2(iTar,iNeu,iE)
  enddo
  close(unitNum)
  if(print_detail) &
& write(*,'(I6,a)') Neus(iNeu),": E = "//trim(Es2str(iTar,iNeu,1))//" to "//trim(Es2str(iTar,iNeu,nE))
enddo
enddo

write(*,'(a)') "For interactions with other nuclei, the database of "&
&//"the nucleus with the closest atomic number and mass number will be used."

! --------------------------------------------------
! [Optionally] Write out only the particles with the kf-code entered in Parameters to the final state.
nkfA=0
do i=2,udm_int_param_nMax
  if( Parameters(i) == 0d0 ) exit
  nkfA=nkfA+1
enddo
allocate( kfAccepted(nkfA) )
if(nkfA>0) print*,"The final states written out by neutrino interactions:"
do i=2,nkfA
  kfAccepted(i)=Parameters(i)
  print*,"-> ",kfAccepted(i)
enddo

! --------------------------------------------------
load_failed=.false. ! Success loading !!
! --------------------------------------------------


end subroutine initialize

!=======================================================================
subroutine find_index_of_target(Z,A,iTar)
!=======================================================================
implicit none
! Z, A are the actual Z and A values.
! iZ, iA are the indices of the data with the Z and A values closest to the actual Z and A.
integer Z,A,iTar
integer dZ   ,dA
integer dZmin,dAmin
integer i

dZmin=999
dAmin=999
do i=1,nTar
  dZ=abs(Z-Zs(i))
  dA=abs(A-As(i))
  if(dZ<=dZmin) then
    dZmin=dZ
    iTar=i
    if(dA<dAmin) then
      dAmin=dA
      iTar=i
    endif
  endif
enddo
end subroutine find_index_of_target

!=======================================================================
subroutine find_index_of_neutrino(iNeu)
!=======================================================================
implicit none
integer i, iNeu
do i=1,nNeu
  if(udm_kf_incident==Neus(i)) then
    iNeu=i
    return
  endif
enddo
end subroutine find_index_of_neutrino

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
integer iTar, iNeu, ndata
double precision EGeV,Emax
double precision, allocatable :: Xdata(:), Ydata(:)

Xsec_per_atom=0d0

if(load_failed) return

iTar=0
iNeu=0
call find_index_of_target(Z,A,iTar)
call find_index_of_neutrino(iNeu)

allocate(Xdata(nEMax1))
allocate(Ydata(nEMax1))

EGeV=Kin/1000d0 ! MeV to GeV
ndata=nE1(iTar,iNeu)
Emax=Es1(iTar,iNeu,ndata)

if(EGeV >= Emax) then
  ! Extrapolation
  Xsec_per_atom=XSs(iTar,iNeu,ndata)*EGeV/Emax*dble(As(iTar))/dble(A)
  return
endif

Xdata=Es1(iTar,iNeu,:)
Ydata=XSs(iTar,iNeu,:)
Xsec_per_atom=udm_interp(EGeV,ndata,Xdata,Ydata,"log-log")*dble(As(iTar))/dble(A)
return
end

!=======================================================================
subroutine find_index_of_energy(iTar,iNeu,iE)
! Find the closest energy in the event data.
!-----------------------------------------------------------------------
! udm_Kin : Kinetic Energy of the incident particle [MeV]
!=======================================================================
implicit none
integer iTar,iNeu,iE,ndata,i
double precision EGeV,Emax,dE,dEmin

! If the energy is greater than the maximum value in the data, always use the maximum value in the data
EGeV=udm_Kin/1000d0
ndata=nE2(iTar,iNeu)
Emax=Es2(iTar,iNeu,ndata)
if(EGeV>=Emax) then
  iE=ndata
  return
endif

dEmin=1e+10
do i=1,ndata
  dE=dabs(EGeV-Es2(iTar,iNeu,i))
  if(dE<dEmin) then
    dEmin=dE
    iE=i
  ! Comment out for cases where the order is not ascending
  ! else
  !   return
  endif
enddo
end subroutine find_index_of_energy

!=======================================================================
subroutine generate_final_state
!=======================================================================
! Subroutine to determine final state information (4-momenta, etc.).
!-----------------------------------------------------------------------
! [Variables available in this subroutine]
! udm_kf_incident : The kf-codes (particle IDs) of the incident particles.
! udm_Kin         : Kinetic Energy of the incident particle [MeV]
!-----------------------------------------------------------------------
integer Z_A_hit(2), Z, A
integer iTar,iNeu,iE
integer iEvent,iFile,iev
character(len=999) :: path
integer i,j,nFS,kf
integer unitNum_this
double precision rKout,Kout,the,phi,phiAdd,mass,Etot,absp
INTEGER thread_id

!-----------------------------------------------------------------------
! You may use Z and A for final state variables.
! Z and A of a target material are automatically obtained using the function 'get_hit_nuclide_Z_A'.
Z_A_hit=get_hit_nuclide_Z_A(udm_Kin)
Z=Z_A_hit(1)
A=Z_A_hit(2)

! Get the indices of the data to be used.
iTar=0
iNeu=0
iE  =0
call find_index_of_target(Z,A,iTar)
call find_index_of_neutrino(iNeu)
call find_index_of_energy(iTar,iNeu,iE)

! Get the file number (iFile) and the event number within it (iev).
iEvent=get_random_int(1,nEvents(iTar,iNeu,iE))
iFile=1+(iEvent-1)/nEventsPerFile
iev=1+mod(iEvent-1,nEventsPerFile)

unitNum_this=unitNum+udm_thread_id()
path=trim(folder)//"/events/"&
  &//str3(Zs(iTar))//str3(As(iTar))//"/"&
  &//trim(udm_int2char(Neus(iNeu)))//"/"&
  &//trim(Es2str(iTar,iNeu,iE))//"/"&
  &//trim(udm_int2char(iFile))//".dat"

open(unitNum_this,file=trim(path),status="old",err=999)
! Skip until reaching the 'iev' event.
do i=1,iev-1
  read(unitNum_this,*) nFS
  do j=1,nFS
    read(unitNum_this,'()')
  enddo
enddo

! Get the 'iev' event
call initialize_udm_event_info

! Get final state number
read(unitNum_this,*) set_final_state_number
if(set_final_state_number > nFSMax) set_final_state_number = nFSMax

! Rotate azimuthal angle
phiAdd=get_random(0.0d0,2.0*3.1415d0)

do i=1,set_final_state_number
  ! ------------------------
  read(unitNum_this,*) rKout,the,phi,kf
  ! ------------------------
  ! calculate kinematic variables
  Kout=rKout*udm_Kin
  phi=phi+phiAdd
  mass=get_mass(kf)
  Etot=Kout+mass
  absp=sqrt(dabs(Etot**2-mass**2))
  ! ------------------------
  ! Restrict the kf to be output
  if( size(kfAccepted)>0 ) then
    if ( .not. any(kf==kfAccepted) ) then
      Etot=mass
      absp=0d0
    endif
  endif
  ! ------------------------
  ! fill
  set_kf                 (i) = kf
  set_Total_Energy_in_MeV(i) = Etot
  set_Px_in_MeV          (i) = absp*sin(the)*cos(phi)
  set_Py_in_MeV          (i) = absp*sin(the)*sin(phi)
  set_Pz_in_MeV          (i) = absp*cos(the)
enddo

close(unitNum_this)

call fill_final_state
return

999 continue
print*,"Warning: Not found:",trim(path),iEvent,nEventsPerFile

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
if( action .eq.  4 ) call print_comments
if(.not. match_initial()) return
if( action .eq. 11 ) call calc_averaged_Xsec
if( action .eq. 21 ) call generate_final_state
end subroutine caller

!=======================================================================
function match_initial()
!=======================================================================
logical match_initial
integer i
match_initial = .false.
do i = 1, num_initial
  if(udm_kf_incident .eq. kf_initial(i)) then
    match_initial = .true.
    exit
  endif
enddo
return
end

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
! special for this module
if(udm_bias(index)<0d0) print_detail=.true.
if(udm_bias(index)/=1d0) then
  write(*,"(a)") trim(Name)//": Bias is set to 1"
  udm_bias(index)=1d0
endif
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
end

!=======================================================================
end module udm_int_neutrino
!=======================================================================






