!***********************************************************************
subroutine read_HDF5tetraparam(filename,melem,mpoint,ierr)
!
! Read parameters of tetra file (hdf5)
!
! Created 2024/07/09 by T. Furuta
!
!***********************************************************************
#ifdef hdf5
!------------------------------------------------------------------------
!     with -Dhdf5 option 
!------------------------------------------------------------------------
  use HDF5
  implicit none
  include 'err.inc'
  character(200),intent(in) :: filename
  integer,intent(out) :: melem,mpoint
  integer,intent(out) :: ierr
  character(200),parameter :: groupname='/geometry'
  character(200),parameter :: attr1='contiguity'
  character(200),parameter :: attr2='element_type'
  character(200),parameter :: attr3='geometry_type'
  character(200),parameter :: data1='ELEMENT_REGIONS'
  character(200),parameter :: data2='ELEMENT_VERTICES'
  character(200),parameter :: data3='NODES'
  integer(HID_T) :: ID_file,ID_group,ID_data,ID_space,ID_attr
  integer(HSIZE_T),allocatable :: ndim(:),maxdim(:)
  integer,allocatable :: idummy(:)
  integer :: ielementtype,icontig,itype,igeotype

  melem=0
  mpoint=0
!------------------------------------------------------------------------
!     hdf5 file open 
!------------------------------------------------------------------------
  call H5open_f(ierr)
  call H5Fopen_f(filename,H5F_ACC_RDONLY_F,ID_file,ierr)
  if(ierr.lt.0)then
   write(ErrCha,'(''*** TETRA ERROR: HDF5 file = '',a,'' is not found'')') trim(filename)
   ErrID = 'L:777/R:read_tetraparam1/F:tetramod.f' !E80_003_002
   call ErrWrite(ErrID,ErrCha)
   ierr=1
   return
  endif
!------------------------------------------------------------------------
!     Group select (group=geomtry)
!------------------------------------------------------------------------
  call H5Gopen_f(ID_file,groupname,ID_group,ierr)
  if(ierr.lt.0)then
   write(ErrCha,'(''*** TETRA ERROR: Incompatible HDF5 format: '',a,'' is not found'')') trim(groupname)
   ErrID = 'L:790/R:read_tetraparam1/F:tetramod.f' !E80_003_003
   call ErrWrite(ErrID,ErrCha)
   ierr=1
   return
  endif
!------------------------------------------------------------------------
!     Attribute check (attr1=contiguity,attr2=element_type)
!------------------------------------------------------------------------
  icontig=0
  ielementtype=0
  igeotype=0
  allocate( idummy(1),ndim(1) )
  call H5Aopen_f(ID_group,attr1,ID_attr,ierr)
  if(ierr.lt.0)then
   write(ErrCha,'(''*** TETRA ERROR: Incompatible HDF5 format: '',a,'' is not found'')') trim(attr1)
   ErrID = 'L:790/R:read_tetraparam1/F:tetramod.f' !E80_003_003
   call ErrWrite(ErrID,ErrCha)
   ierr=1
   return
  else
   call H5Aread_f(ID_attr,H5T_NATIVE_INTEGER,idummy,ndim,ierr)
   call H5Aclose_f(ID_attr,ierr)
   icontig=idummy(1)
  endif
  call H5Aopen_f(ID_group,attr2,ID_attr,ierr)
  if(ierr.lt.0)then
   write(ErrCha,'(''*** TETRA ERROR: Incompatible HDF5 format: '',a,'' is not found'')') trim(attr2)
   ErrID = 'L:790/R:read_tetraparam1/F:tetramod.f' !E80_003_003
   call ErrWrite(ErrID,ErrCha)
   ierr=1
   return
  else
   call H5Aread_f(ID_attr,H5T_NATIVE_INTEGER,idummy,ndim,ierr)
   call H5Aclose_f(ID_attr,ierr)
   ielementtype=idummy(1)
  endif
  call H5Aopen_f(ID_group,attr3,ID_attr,ierr)
  if(ierr.lt.0)then
   write(ErrCha,'(''*** TETRA ERROR: Incompatible HDF5 format: '',a,'' is not found'')') trim(attr2)
   ErrID = 'L:790/R:read_tetraparam1/F:tetramod.f' !E80_003_003
   call ErrWrite(ErrID,ErrCha)
   ierr=1
   return
  else
   call H5Aread_f(ID_attr,H5T_NATIVE_INTEGER,idummy,ndim,ierr)
   call H5Aclose_f(ID_attr,ierr)
   igeotype=idummy(1)
  endif
  if(ielementtype.ne.1.or.igeotype.ne.1.or.(icontig.eq.0.or.icontig.eq.5.or.icontig.eq.6))then
   write(ErrCha,'(''*** TETRA ERROR: Incompatible HDF5 format: icontig = '',i1,'' ielementtype = '',i1)') icontig,ielementtype
   ErrID = 'L:790/R:read_tetraparam1/F:tetramod.f' !E80_003_003
   call ErrWrite(ErrID,ErrCha)
   ierr=1
   return
  endif
  deallocate( ndim,idummy )
!------------------------------------------------------------------------
!     Data check data1=element_regions
!                data2=element_vertices
!                data3=nodes      
!------------------------------------------------------------------------
  call H5Dopen_f(ID_group,data1,ID_data,ierr)
  if(ierr.lt.0)then
   write(ErrCha,'(''*** TETRA ERROR: Incompatible HDF5 format: '',a,'' is not found'')') trim(data1)
   ErrID = 'L:790/R:read_tetraparam1/F:tetramod.f' !E80_003_003
   call ErrWrite(ErrID,ErrCha)
   ierr=1
   return
  else
   call H5Dget_space_f(ID_data,ID_space,ierr)
   allocate( ndim(1),maxdim(1) )
   call H5Sget_simple_extent_dims_f(ID_space,ndim,maxdim,ierr)
   melem=ndim(1)
   deallocate( ndim,maxdim )
   ierr=0
  endif
  call H5Dopen_f(ID_group,data3,ID_data,ierr)
  if(ierr.lt.0)then
   write(ErrCha,'(''*** TETRA ERROR: Incompatible HDF5 format: '',a,'' is not found'')') trim(data3)
   ErrID = 'L:790/R:read_tetraparam1/F:tetramod.f' !E80_003_003
   call ErrWrite(ErrID,ErrCha)
   ierr=1
   return
  else
   call H5Dget_space_f(ID_data,ID_space,ierr)
   allocate( ndim(2),maxdim(2) )
   call H5Sget_simple_extent_dims_f(ID_space,ndim,maxdim,ierr)
   mpoint=ndim(2)
   deallocate( ndim,maxdim )
   ierr=0
  endif
  return
#else
!------------------------------------------------------------------------
!     without -Dhdf5 option 
!------------------------------------------------------------------------
  implicit none
  include 'err.inc'
  character(200),intent(in) :: filename
  integer,intent(out) :: melem,mpoint
  integer,intent(out) :: ierr
  melem=0
  mpoint=0
  write(ErrCha,'(''*** TETRA ERROR: HDF5 file can not be use with this executable.'')')
  call ErrWrite(ErrID,ErrCha)
  write(*,*)'Recompile PHITS following instruction phits/utility/HDF5'
  ierr=1
  return
#endif
end subroutine read_HDF5tetraparam

!***********************************************************************
subroutine read_HDF5tetradata(filename,melem,mpoint,ie2u,ie2p,pxyz,ierr)
!
! Read tetra file (hdf5) and set ielem2univ, ielem2point, pointxyz 
!
! Created 2024/07/09 by T. Furuta
!
!***********************************************************************
#ifdef hdf5
!------------------------------------------------------------------------
!     with -Dhdf5 option 
!------------------------------------------------------------------------
  use HDF5
  implicit none
  include 'err.inc'
  character(200),intent(in) :: filename
  integer,intent(in) :: melem,mpoint
  integer,intent(out) :: ie2u(melem),ie2p(4,melem)
  real(8),intent(out) :: pxyz(3,mpoint)
  integer,intent(out) :: ierr
  character(200),parameter :: groupname='/geometry'
  character(200),parameter :: data1='ELEMENT_REGIONS'
  character(200),parameter :: data2='ELEMENT_VERTICES'
  character(200),parameter :: data3='NODES'
  integer(HID_T) :: ID_file,ID_group,ID_data,ID_space
  integer(HSIZE_T),allocatable :: ndim(:),maxdim(:)
  integer,allocatable :: idummy(:)
  integer :: i,itype

!------------------------------------------------------------------------
!     hdf5 file open 
!------------------------------------------------------------------------
  call H5open_f(ierr)
  call H5Fopen_f(filename,H5F_ACC_RDONLY_F,ID_file,ierr)
!------------------------------------------------------------------------
!     Group select (group=geometry)
!------------------------------------------------------------------------
  call H5Gopen_f(ID_file,groupname,ID_group,ierr)
!------------------------------------------------------------------------
!     Data read data1=element_regions
!------------------------------------------------------------------------
  call H5Dopen_f(ID_group,data1,ID_data,ierr)
  call H5Dget_space_f(ID_data,ID_space,ierr)
  allocate( ndim(1),maxdim(1) )
  call H5Sget_simple_extent_dims_f(ID_space,ndim,maxdim,ierr)
  call H5Dread_f(ID_data,H5T_NATIVE_INTEGER,ie2u,ndim,ierr)
  if(ierr.lt.0)then
   write(ErrCha,'(''*** TETRA ERROR: Incompatible HDF5 format: '',a,'' is not found'')') trim(data1)
   call ErrWrite(ErrID,ErrCha)
   ierr=1
   return
  endif
  call H5Dclose_f(ID_data,ierr)
  deallocate( ndim,maxdim )
!------------------------------------------------------------------------
!     Data read data2=element_vertices
!------------------------------------------------------------------------
  call H5Dopen_f(ID_group,data2,ID_data,ierr)
  call H5Dget_space_f(ID_data,ID_space,ierr)
  allocate( ndim(2),maxdim(2) )
  call H5Sget_simple_extent_dims_f(ID_space,ndim,maxdim,ierr)
  call H5Dread_f(ID_data,H5T_NATIVE_INTEGER,ie2p,ndim,ierr)
  if(ierr.lt.0)then
   write(ErrCha,'(''*** TETRA ERROR: Incompatible HDF5 format: '',a,'' is not found'')') trim(data2)
   call ErrWrite(ErrID,ErrCha)
   ierr=1
   return
  endif
  call H5Dclose_f(ID_data,ierr)
  deallocate( ndim,maxdim )
!------------------------------------------------------------------------
!     Data read data3=nodes      
!------------------------------------------------------------------------
  call H5Dopen_f(ID_group,data3,ID_data,ierr)
  call H5Dget_space_f(ID_data,ID_space,ierr)
  allocate( ndim(2),maxdim(2) )
  call H5Sget_simple_extent_dims_f(ID_space,ndim,maxdim,ierr)
  call H5Dread_f(ID_data,H5T_NATIVE_DOUBLE,pxyz,ndim,ierr)
  if(ierr.lt.0)then
   write(ErrCha,'(''*** TETRA ERROR: Incompatible HDF5 format: '',a,'' is not found'')') trim(data3)
   call ErrWrite(ErrID,ErrCha)
   ierr=1
   return
  endif
  call H5Dclose_f(ID_data,ierr)
  ierr=0
  deallocate( ndim,maxdim )
  return
#else
!------------------------------------------------------------------------
!     without -Dhdf5 option 
!------------------------------------------------------------------------
  implicit none
  include 'err.inc'
  character(200),intent(in) :: filename
  integer,intent(in) :: melem,mpoint
  integer,intent(out) :: ie2u(melem),ie2p(4,melem)
  real(8),intent(out) :: pxyz(3,mpoint)
  integer,intent(out) :: ierr
  integer :: i
  do i=1,melem
   ie2u(i)=0
   ie2p(1:4,i)=0
  enddo
  do i=1,mpoint
   pxyz(1:3,i)=0.0d0
  enddo
  write(ErrCha,'(''*** TETRA ERROR: HDF5 file can not be use with this executable.'')')
  call ErrWrite(ErrID,ErrCha)
  write(*,*)'Recompile PHITS following instruction phits/utility/HDF5'
  ierr=1
  return
#endif  
end subroutine read_HDF5tetradata

!***********************************************************************
subroutine read_HDF5tetramat(filename,nuniv,nmatst,iu2m,du2d,nummat,ierr)
!
! Read parameters of tetra file (hdf5)
!
! Created 2024/07/09 by T. Furuta
!
!***********************************************************************
#ifdef hdf5
!------------------------------------------------------------------------
!     with -Dhdf5 option 
!------------------------------------------------------------------------
  use HDF5
  use ISO_C_BINDING
  implicit none
  include 'err.inc'
  character(200),intent(in) :: filename
  integer,intent(in) :: nuniv
  integer,intent(in) :: nmatst
  integer,intent(out) :: iu2m(nuniv)
  real(8),intent(out) :: du2d(nuniv)
  integer,intent(out) :: nummat
  integer,intent(out) :: ierr
  character(200),parameter :: groupname='/geometry'
  character(200),parameter :: matgroupname='materials'
  character(200),parameter :: data1='REGION_MATERIALS'
  character(200),parameter :: data2='REGION_DENSITY'
  character(200),parameter :: data3='ISOTOPE_NAMES'
  character(200),parameter :: data4='ISOTOPE_FRACTIONS'
  integer(HID_T) :: ID_file,ID_group,ID_mat,ID_space
  integer(HID_T) :: ID_data,ID_data2,ID_eachmat,ID_type
  integer(HSIZE_T),allocatable :: ndim(:),maxdim(:)
  integer :: i,j,k,len
  character(80) :: cmatname(nuniv)
  character(200),allocatable :: cisonames(:)
  real(8),allocatable :: disofracs(:)
  TYPE(C_PTR), DIMENSION(:), ALLOCATABLE, TARGET :: rdata
  CHARACTER(len = 8, kind=c_char),  POINTER :: data
  TYPE(C_PTR) :: f_ptr
  
!------------------------------------------------------------------------
!     hdf5 file open 
!------------------------------------------------------------------------
  call H5open_f(ierr)
  call H5Fopen_f(filename,H5F_ACC_RDONLY_F,ID_file,ierr)
!------------------------------------------------------------------------
!     Group select (group=geomtry)
!------------------------------------------------------------------------
  call H5Gopen_f(ID_file,groupname,ID_group,ierr)
!------------------------------------------------------------------------
!     Material group select (group=materials)
!------------------------------------------------------------------------
  call H5Gopen_f(ID_group,matgroupname,ID_mat,ierr)
  if(ierr.lt.0)then
   write(ErrCha,'(''*** TETRA ERROR: Incompatible HDF5 format: '',a,'' is not found'')') trim(matgroupname)
   ErrID = 'L:790/R:read_tetraparam1/F:tetramod.f' !E80_003_003
   call ErrWrite(ErrID,ErrCha)
   ierr=1
   return
  endif
  call H5GN_members_f(ID_group,matgroupname,nummat,ierr)
  if(ierr.lt.0)then
   write(ErrCha,'(''*** TETRA ERROR: Incompatible HDF5 format: '',a,'' is not found'')') trim(matgroupname)
   ErrID = 'L:790/R:read_tetraparam1/F:tetramod.f' !E80_003_003
   call ErrWrite(ErrID,ErrCha)
   ierr=1
   return
  endif
  do i=0,nummat-1
   call H5Gget_obj_info_idx_f(ID_group,matgroupname,i,cmatname(i+1),H5G_DATASET_F,ierr)
   call H5Gopen_f(ID_mat,cmatname(i+1),ID_eachmat,ierr)
   call H5Dopen_f(ID_eachmat,data3,ID_data,ierr)
   call H5Dopen_f(ID_eachmat,data4,ID_data2,ierr)
   call H5Dget_space_f(ID_data,ID_space,ierr)
   allocate( ndim(1),maxdim(1) )
   call H5Sget_simple_extent_dims_f(ID_space,ndim,maxdim,ierr)
   allocate( cisonames(ndim(1)),disofracs(ndim(1)),rdata(ndim(1)) )
   call H5Dget_type_f(ID_data,ID_type,ierr)
   f_ptr = C_LOC(rdata(1))
   call H5Dread_f(ID_data,ID_type,f_ptr,ierr)
   do j=1,ndim(1)
    cisonames(j)(1:200)=' '
    call C_F_POINTER(rdata(j),data)
    len=0
    do
     if(data(len+1:len+1).eq.C_NULL_CHAR.or.len.ge.8) exit
     len = len + 1
    enddo
    cisonames(j)(1:len)=data(1:len)
   enddo
   call H5Dread_f(ID_data2,H5T_NATIVE_DOUBLE,disofracs,ndim,ierr)
   if(ierr.lt.0)then
    write(ErrCha,'(''*** TETRA ERROR: Incompatible HDF5 format: '',a,'' is not found'')') trim(matgroupname)
    ErrID = 'L:790/R:read_tetraparam1/F:tetramod.f' !E80_003_003
    call ErrWrite(ErrID,ErrCha)
    ierr=1
    return
   endif
   call H5Dclose_f(ID_data,ierr)
   call H5Dclose_f(ID_data2,ierr)
   call H5Gclose_f(ID_eachmat,ierr)
   call tetraHDF5setmat(i+1,nmatst,ndim(1),cisonames,disofracs,cmatname(i+1))
   deallocate( ndim,maxdim,cisonames,disofracs,rdata )
  enddo
!------------------------------------------------------------------------
!     Data read data1=MATERILA_NAMES
!------------------------------------------------------------------------
  call H5Dopen_f(ID_group,data1,ID_data,ierr)
  call H5Dget_space_f(ID_data,ID_space,ierr)
  allocate( ndim(1),maxdim(1) )
  call H5Sget_simple_extent_dims_f(ID_space,ndim,maxdim,ierr)
  allocate( rdata(ndim(1)) )
  call H5Dget_type_f(ID_data,ID_type,ierr)
  f_ptr = C_LOC(rdata(1))
  call H5Dread_f(ID_data,ID_type,f_ptr,ierr)
  do j=1,ndim(1)
   call C_F_POINTER(rdata(j),data)
   len=0
   do
    if(data(len+1:len+1).eq.C_NULL_CHAR.or.len.ge.8) exit
    len = len + 1
   enddo
   do k=1,nummat
    if(data(1:len).eq.cmatname(k)(1:len))then
     iu2m(j)=k
     exit
    endif
   enddo
  enddo
  deallocate( rdata,ndim,maxdim )
!------------------------------------------------------------------------
!     Data read data2=REGION_DENSITY
!------------------------------------------------------------------------
  call H5Dopen_f(ID_group,data2,ID_data,ierr)
  call H5Dget_space_f(ID_data,ID_space,ierr)
  allocate( ndim(1),maxdim(1) )
  call H5Sget_simple_extent_dims_f(ID_space,ndim,maxdim,ierr)
  call H5Dread_f(ID_data,H5T_NATIVE_DOUBLE,du2d,ndim,ierr)
  deallocate( ndim,maxdim )
  return
#else
!------------------------------------------------------------------------
!     without -Dhdf5 option 
!------------------------------------------------------------------------
  implicit none
  include 'err.inc'
  character(200),intent(in) :: filename
  integer,intent(in) :: nuniv
  integer,intent(in) :: nmatst
  integer,intent(out) :: iu2m(nuniv)
  real(8),intent(out) :: du2d(nuniv)
  integer,intent(out) :: nummat
  integer,intent(out) :: ierr
  iu2m=0
  du2d=0.0d0
  nummat=0
  write(ErrCha,'(''*** TETRA ERROR: HDF5 file can not be use with this executable.'')')
  call ErrWrite(ErrID,ErrCha)
  write(*,*)'Recompile PHITS following instruction phits/utility/HDF5'
  ierr=1
  return
#endif
end subroutine read_HDF5tetramat

