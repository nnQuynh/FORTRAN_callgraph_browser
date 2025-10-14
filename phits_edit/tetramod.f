************************************************************************
*                                                                      *
      module TETRAMOD
*                                                                      *
*     Modules for tetrahedron mesh                                     *
*
*     Following subroutines are called externally
*       tetrainit, tetrafin
*       tetrafnd
*       tetracald
*       tetranext
*       tetraangl
*       tetrachk
*       tetrasorsinit
*       tetrasors
*
*     Created by T. Furuta on 2015/07/14
*     Last modified 2023/03/16 by T. Furuta
*                                                                      *
************************************************************************
      implicit none
      integer :: itetra,ntetsurf,ntetelem
      logical,parameter :: iwarning=.false.,iwarning2=.false.
      integer,parameter :: ndiv=9 ! max of 4 bites
      integer,parameter :: mxind=2147483647 ! max of 4 bites
      integer,parameter :: numlarge=99999999
      real(8),parameter :: xtiny=1d-5
      integer :: itetragshow !FURUTA20221209
      integer :: npoint(0:10),nsurf(0:10),nelem(0:10)
      integer :: noutsf(0:10),noutpt(0:10) !FURUTA20221102
      integer :: noctOsec(0:10),noctEblk(0:10)
      integer :: nelemtot,ntetcl
      real(8) :: xoctOsec(3,ndiv,10),xb(6,10)
      real(8) :: xoctEblk(3,ndiv,10)
      real(8) :: xqt(10)
      integer :: lbin
      character(200) :: binfile,binfile2
      real(8),allocatable :: pointxyz(:,:)
      real(8),allocatable :: surfcoefs(:,:)
      integer,allocatable :: ielem2point(:,:)
      integer,allocatable :: ielem2surf(:,:)
      integer,allocatable :: isurf2elem(:,:)
      integer,allocatable :: ielem2univ(:)
      integer,allocatable :: ielem2icl(:)
      integer,allocatable :: ioutsfmapOsec(:),ielemmapEblk(:)
      integer,allocatable :: indmapOsec(:,:),indmapEblk(:,:)
      integer,allocatable :: ioctOsec(:),ioctEblk(:)
      integer,allocatable :: ioutsf2outpt(:,:)
      integer,allocatable :: ioutsf2surf(:),ioutpt2point(:)
      integer,allocatable :: itetlist(:),itetst(:),itetelem(:)
      real(8),allocatable :: welemtot(:),welem(:)
      integer,allocatable :: ielem2id(:) !FURUTA20190110
      real(8),allocatable :: volelems(:) !FURUTA20190110
*------------------------------------------------------------------------
*     Temporary variables used in tetrainit
*------------------------------------------------------------------------
      integer :: nOsec(0:10),nindOsec(0:10)
      integer :: nEblk(0:10),nindEblk(0:10)
      integer,allocatable :: ioutpt2outsf(:,:)
      integer,allocatable :: iside2outsf(:),isidemap(:)
      integer,allocatable :: ioutsf2side(:,:),iside2outpt(:,:)
      integer,allocatable :: isfind(:,:)
*------------------------------------------------------------------------
      contains

************************************************************************
      subroutine tetrafin
*
*     Finalize tetrahedron mesh
*
************************************************************************
      call deallocate_tetratbl
      call deallocate_blocktbl
      call deallocate_outertbl
      call deallocate_itetlist
      return
      end subroutine tetrafin

************************************************************************
      subroutine tetraread(tetsfac,nlat3,ltfile,itfform,tfilename,ierr)
*
*     Initiate tetrahedron mesh
*
*     Last modified 2017/10/20 by T. Furuta
*
************************************************************************
      real(8),intent(in) :: tetsfac(10)
      integer,intent(in) :: nlat3,ltfile(10),itfform(10)
      integer,intent(out) :: ierr
      character(200),intent(in) :: tfilename(10)
      integer :: iot1,iot2,itet
      character(200) :: nodefile,elemfile,bdffile,hdf5file

      ierr=0
      iot1=68
      iot2=69
*------------------------------------------------------------------------
      if(abs(itetra).eq.1)then
       call read_infofile(iot1,nlat3)
       call set_xquota(nlat3)
       return
      endif
*------------------------------------------------------------------------
      do itet=1,nlat3
       if(itfform(itet).eq.0)then
        write(nodefile,'(200x)')
        write(elemfile,'(200x)')
        nodefile(1:ltfile(itet)+5)
     &       =tfilename(itet)(1:ltfile(itet))//'.node'
        elemfile(1:ltfile(itet)+4)
     &       =tfilename(itet)(1:ltfile(itet))//'.ele'
        open(iot1,file=nodefile,status='old')
        open(iot2,file=elemfile,status='old')
        call read_tetraparam0(itet,iot1,iot2,ierr)
        close(iot1)
        close(iot2)
       elseif(itfform(itet).eq.1)then
        write(bdffile,'(200x)')
        bdffile(1:ltfile(itet)+4)
     &       =tfilename(itet)(1:ltfile(itet))
        open(iot1,file=bdffile,status='old')
        call read_tetraparam1(itet,iot1,ierr)
        close(iot1)
       elseif(itfform(itet).eq.2)then
        write(hdf5file,'(200x)')
        hdf5file(1:ltfile(itet))
     &       =tfilename(itet)(1:ltfile(itet))
        call read_tetraparam2(itet,hdf5file,ierr)
       endif
      enddo
      if(ierr.ne.0)return
      call allocate_tetratbl1(nlat3)
      do itet=1,nlat3
       if(itfform(itet).eq.0)then
        write(nodefile,'(200x)')
        write(elemfile,'(200x)')
        nodefile(1:ltfile(itet)+5)
     &       =tfilename(itet)(1:ltfile(itet))//'.node'
        elemfile(1:ltfile(itet)+4)
     &       =tfilename(itet)(1:ltfile(itet))//'.ele'
        open(iot1,file=nodefile,status='old')
        open(iot2,file=elemfile,status='old')
        call read_tetradata0(tetsfac(itet),itet,iot1,iot2,ierr)
        close(iot1)
        close(iot2)
        if(ierr.ne.0)return
       elseif(itfform(itet).eq.1)then
        write(bdffile,'(200x)')
        bdffile(1:ltfile(itet))=tfilename(itet)(1:ltfile(itet))
        open(iot1,file=bdffile,status='old')
        call read_tetradata1(tetsfac(itet),itet,iot1,ierr)
        close(iot1)
       elseif(itfform(itet).eq.2)then
        write(hdf5file,'(200x)')
        hdf5file(1:ltfile(itet))
     &       =tfilename(itet)(1:ltfile(itet))
        call read_tetradata2(tetsfac(itet),itet,hdf5file,ierr)
       endif
       if(ierr.ne.0)return
      enddo
      return
      end subroutine tetraread

************************************************************************
      subroutine tetrainit(io,nlat3,itgchk,coincd,ierr)
*
*     Initiate tetrahedron mesh
*
*     Last modified 2023/03/16 by T. Furuta
*
************************************************************************
      integer,intent(in) :: io,nlat3,itgchk
      real(8),intent(in) :: coincd
      integer,intent(out) :: ierr
      integer :: iot1,itet,kkk,mxsf4p,nside,nside0,np0,ns0
      real(8) :: xx0(6)

      ierr=0
      iot1=68
*------------------------------------------------------------------------
      if(abs(itetra).eq.1)then !FURUTA20190331
       nelemtot=nelem(nlat3)
       return
      endif
*------------------------------------------------------------------------
      do itet=1,nlat3
       if(ierr.ne.0)exit
       kkk=10000+itet
       call tetrabox(io,kkk,xx0,ierr)
       call check_inside(xx0,itet,ierr)
      enddo
      call set_xquota(nlat3) !FURUTA20190331
      if(ierr.ne.0)return
      call reorder_ielem2point(nlat3)
      call allocate_tmptbl1(nlat3)
      do itet=1,nlat3
       call set_nsurf(itet)
      enddo
      call allocate_tetratbl2(nlat3)
      do itet=1,nlat3
       call set_tetrasurf(itet)
       call set_isurf2elem(itet)
      enddo
      call deallocate_tmptbl1
      do itet=1,nlat3
       kkk=10000+itet
       call tetrabox(io,kkk,xx0,ierr)
       call count_block(itet,
     &      nEblk(itet),nindEblk(itet),noctEblk(itet),ierr)
       if(ierr.ne.0)return
       call set_outtbl(itet)
       call calc_mxsf(itet,mxsf4p,nside0)
       call allocate_tmptbl2(itet,mxsf4p,nside0)
       call set_tmptbl(itet,mxsf4p,nside0,nside,ierr)
       if(ierr.ne.0)return
       call count_outer(xx0,itet,nside,
     &      nOsec(itet),nindOsec(itet),noctOsec(itet),ierr)
       if(ierr.ne.0)return
       call deallocate_tmptbl2
      enddo
      if(ierr.ne.0)return
      call set_numblock(nlat3)
      call set_numouter(nlat3)
      call allocate_blocktbl(nEblk(nlat3),nindEblk(nlat3),
     &     noctEblk(nlat3))
      call allocate_outertbl(nOsec(nlat3),nindOsec(nlat3),
     &     noctOsec(nlat3))
      do itet=1,nlat3
       kkk=10000+itet
       call tetrabox(io,kkk,xx0,ierr)
       call set_block(itet,ierr)
       call calc_mxsf(itet,mxsf4p,nside0)
       call allocate_tmptbl2(itet,mxsf4p,nside0)
       call set_tmptbl(itet,mxsf4p,nside0,nside,ierr)
       call set_outer(xx0,itet,nside,ierr)
       call deallocate_tmptbl2
      enddo
      if(abs(itetra).eq.2)then
        call write_infofile(iot1,nlat3)
        write(*,*)'TETRA INFO file was successfully generated!!'
        call parastop(300)
      endif
      if(itgchk.gt.0)then
       call check_geom(nlat3,itgchk,coincd,ierr)
      endif
      nelemtot=nelem(nlat3)
      return
      end subroutine tetrainit

************************************************************************
      subroutine write_infofile(iot1,nlat3)
*
*     Output processed info to file
*
*     Last modified 2016/09/08 by T. Furuta
*
************************************************************************
      integer,intent(in) :: nlat3,iot1
      integer :: i,itet,idiv
      if(itetra.lt.0)then
        open(iot1,file='Tetra.param',status='unknown')
        write(iot1,5000)0,ntetsurf,ntetelem !FURUTA20190331
        do itet=1,nlat3
         write(iot1,5000)itet,npoint(itet)-npoint(itet-1),
     &        nelem(itet)-nelem(itet-1),
     &        nsurf(itet)-nsurf(itet-1),
     &        noutpt(itet)-npoint(itet-1), !FURUTA20221220
     &        noutsf(itet)-nsurf(itet-1), !FURUTA20221220
     &        noctOsec(itet)-noctOsec(itet-1),
     &        noctEblk(itet)-noctEblk(itet-1)
        enddo
 5000   format(10(i9,' '))
        do itet=1,nlat3
         write(iot1,5001)itet,xb(1:6,itet)
        enddo
 5001   format(i9,' ',6(1pe22.15,' '))
        do itet=1,nlat3
         do idiv=1,ndiv
          write(iot1,5002)itet,idiv,
     &         xoctOsec(1:3,idiv,itet),xoctEblk(1:3,idiv,itet)
         enddo
        enddo
 5002   format(2(i9,' '),6(1pe22.15,' '))
        write(iot1,5000)nOsec(nlat3),nindOsec(nlat3),
     &       nEblk(nlat3),nindEblk(nlat3)
        close(iot1)
        open(iot1,file='Tetra.point',status='unknown')
        do i=1,npoint(nlat3)
         write(iot1,5001)i,pointxyz(1:3,i)
        enddo
        close(iot1)
        open(iot1,file='Tetra.surf',status='unknown')
        do i=1,nsurf(nlat3)
         write(iot1,5003)i,isurf2elem(1:2,i),surfcoefs(1:4,i)
        enddo
 5003   format(3(i9,' '),4(1pe22.15,' '))
        close(iot1)
        open(iot1,file='Tetra.elem',status='unknown')
        do i=1,nelem(nlat3)
         write(iot1,5000)ielem2id(i),ielem2point(1:4,i),
     &        ielem2surf(1:4,i),ielem2univ(i) !FURUTA20191028 bugfix
        enddo
        close(iot1)
        open(iot1,file='Tetra.otsf',status='unknown')
        do i=1,noutsf(nlat3)
         write(iot1,5000)i,ioutsf2outpt(1:3,i),ioutsf2surf(i)
        enddo
        close(iot1)
        open(iot1,file='Tetra.otpt',status='unknown')
        do i=1,noutpt(nlat3)
         write(iot1,5000)i,ioutpt2point(i)
        enddo
        close(iot1)
        open(iot1,file='Tetra.octO',status='unknown')
        write(iot1,5000)(ioctOsec(i),i=1,noctOsec(nlat3))
        close(iot1)
        open(iot1,file='Tetra.octE',status='unknown')
        write(iot1,5000)(ioctEblk(i),i=1,noctEblk(nlat3))
        close(iot1)
        open(iot1,file='Tetra.indO',status='unknown')
        write(iot1,5000)(indmapOsec(i,1),i=1,nindOsec(nlat3))
        write(iot1,5000)(indmapOsec(i,2),i=1,nindOsec(nlat3))
        close(iot1)
        open(iot1,file='Tetra.indE',status='unknown')
        write(iot1,5000)(indmapEblk(i,1),i=1,nindEblk(nlat3))
        write(iot1,5000)(indmapEblk(i,2),i=1,nindEblk(nlat3))
        close(iot1)
        open(iot1,file='Tetra.mapO',status='unknown')
        write(iot1,5000)(ioutsfmapOsec(i),i=1,nOsec(nlat3))
        close(iot1)
        open(iot1,file='Tetra.mapE',status='unknown')
        write(iot1,5000)(ielemmapEblk(i),i=1,nEblk(nlat3))
        close(iot1)
      else
        open(iot1,file=binfile,form='unformatted',status='unknown')
        write(iot1)0,ntetsurf,ntetelem !FURUTA20190331
        do itet=1,nlat3
         write(iot1)npoint(itet)-npoint(itet-1),
     &        nelem(itet)-nelem(itet-1),
     &        nsurf(itet)-nsurf(itet-1),
     &        noutpt(itet)-npoint(itet-1), !FURUTA20221220
     &        noutsf(itet)-nsurf(itet-1), !FURUTA20221220
     &        noctOsec(itet)-noctOsec(itet-1),
     &        noctEblk(itet)-noctEblk(itet-1)
        enddo
        do itet=1,nlat3
         write(iot1)xb(1:6,itet)
        enddo
        do itet=1,nlat3
         do idiv=1,ndiv
          write(iot1)xoctOsec(1:3,idiv,itet),xoctEblk(1:3,idiv,itet)
         enddo
        enddo
        write(iot1)nOsec(nlat3),nindOsec(nlat3),
     &       nEblk(nlat3),nindEblk(nlat3)
        do i=1,npoint(nlat3)
         write(iot1)pointxyz(1:3,i)
        enddo
        do i=1,nsurf(nlat3)
         write(iot1)isurf2elem(1:2,i),surfcoefs(1:4,i)
        enddo
        do i=1,nelem(nlat3)
         write(iot1)ielem2id(i),ielem2point(1:4,i),
     &        ielem2surf(1:4,i),ielem2univ(i)
        enddo
        do i=1,noutsf(nlat3)
         write(iot1)ioutsf2outpt(1:3,i),ioutsf2surf(i)
        enddo
        do i=1,noutpt(nlat3)
         write(iot1)ioutpt2point(i)
        enddo
        write(iot1)(ioctOsec(i),i=1,noctOsec(nlat3))
        write(iot1)(ioctEblk(i),i=1,noctEblk(nlat3))
        write(iot1)(indmapOsec(i,1),i=1,nindOsec(nlat3))
        write(iot1)(indmapOsec(i,2),i=1,nindOsec(nlat3))
        write(iot1)(indmapEblk(i,1),i=1,nindEblk(nlat3))
        write(iot1)(indmapEblk(i,2),i=1,nindEblk(nlat3))
        write(iot1)(ioutsfmapOsec(i),i=1,nOsec(nlat3))
        write(iot1)(ielemmapEblk(i),i=1,nEblk(nlat3))
        close(iot1)
      endif
      return
      end subroutine write_infofile

************************************************************************
      subroutine write_infofile2(iot1,nlat3,nunivtot,nmats,
     &     iunivmat,materials,densuniv)
*
*     Output processed info to file for itetauto=1
*
*     Created 2019/03/31 by T. Furuta
*
************************************************************************
      integer,intent(in) :: nunivtot,nlat3,iot1
      integer,intent(in) :: nmats(0:10)
      integer,intent(in) :: iunivmat(nunivtot),materials(nunivtot)
      real(8),intent(in) :: densuniv(nunivtot)
      integer :: i,itet
      if(itetra.lt.0)then
       open(iot1,file='Tetra.univ',status='unknown')
       write(iot1,'(i9)')nunivtot
       do i=1,nlat3
        write(iot1,5001),i,nmats(i)
       enddo
       do i=1,nunivtot
        write(iot1,5000)i,iunivmat(i),densuniv(i)
       enddo
       do i=1,nmats(nlat3)
        write(iot1,5001)i,materials(i)
       enddo
 5000  format(i9,' ',i9,' ',f15.8)
 5001  format(i9,' ',i9)
       close(iot1)
      else
       open(iot1,file=binfile2,form='unformatted',status='unknown')
       write(iot1)nunivtot,nmats(1:nlat3)
       do i=1,nunivtot
        write(iot1)iunivmat(i)
       enddo
       do i=1,nunivtot
        write(iot1)densuniv(i)
       enddo
       do i=1,nmats(nlat3)
        write(iot1)materials(i)
       enddo
       close(iot1)
      endif
      return
      end subroutine write_infofile2

************************************************************************
      subroutine read_infofile(iot1,nlat3)
*
*     Read info from file
*
*     Last modified 2016/09/08 by T. Furuta
*
************************************************************************
      integer,intent(in) :: nlat3,iot1
      integer :: i,i0,itet,itet0,idiv,idiv0,n1o,n2o,n1e,n2e
      integer :: iflag !FURUTA20190331
      logical :: exex
      if(itetra.lt.0)then
        open(iot1,file='Tetra.param',status='old')
c-----------------------------------------------------------------------
        read(iot1,5000)iflag,ntetsurf,ntetelem
        if(iflag.ne.0)then
         write(*,*)'*** TETRA ERROR: ',
     &        'inconsistent file format'
         write(*,*)'Tetra files should be re-created ',
     &        'with itetra=2'
         call parastop(301)
        endif
c-----------------------------------------------------------------------
        do itet=1,nlat3
         read(iot1,5000)itet0,npoint(itet),nelem(itet),nsurf(itet),
     &        noutpt(itet),noutsf(itet),noctOsec(itet),noctEblk(itet)
        enddo
 5000   format(10(i9,' '))
        do itet=1,nlat3
         read(iot1,5001)itet0,xb(1:6,itet)
        enddo
 5001   format(i9,' ',6(1pe22.15,' '))
        do itet=1,nlat3
         do idiv=1,ndiv
          read(iot1,5002)itet0,idiv0,
     &         xoctOsec(1:3,idiv,itet),xoctEblk(1:3,idiv,itet)
         enddo
        enddo
 5002   format(2(i9,' '),6(1pe22.15,' '))
        read(iot1,5000)n1o,n2o,n1e,n2e
        close(iot1)
        call allocate_tetratbl1(nlat3)
        call allocate_tetratbl2(nlat3)
        noutpt(0)=0 !FURUTA20221220
        noutsf(0)=0 !FURUTA20221220
        noctOsec(0)=0
        noctEblk(0)=0
        do itet=1,nlat3
         noutpt(itet)=noutpt(itet)+npoint(itet-1) !FURUTA20221220
         noutsf(itet)=noutsf(itet)+nsurf(itet-1)  !FURUTA20221220
         noctOsec(itet)=noctOsec(itet)+noctOsec(itet-1)
         noctEblk(itet)=noctEblk(itet)+noctEblk(itet-1)
        enddo
        call allocate_blocktbl(n1e,n2e,noctEblk(nlat3))
        call allocate_outertbl(n1o,n2o,noctOsec(nlat3))
        open(iot1,file='Tetra.point',status='old')
        do i=1,npoint(nlat3)
         read(iot1,5001)i0,pointxyz(1:3,i)
        enddo
        close(iot1)
        open(iot1,file='Tetra.surf',status='old')
        do i=1,nsurf(nlat3)
         read(iot1,5003)i0,isurf2elem(1:2,i),surfcoefs(1:4,i)
        enddo
 5003   format(3(i9,' '),4(1pe22.15,' '))
        close(iot1)
        open(iot1,file='Tetra.elem',status='old')
        do i=1,nelem(nlat3)
         read(iot1,5000)ielem2id(i),ielem2point(1:4,i),
     &        ielem2surf(1:4,i),ielem2univ(i) !FURUTA20191028 bugfix
        enddo
        close(iot1)
        open(iot1,file='Tetra.otsf',status='old')
        do i=1,noutsf(nlat3)
         read(iot1,5000)i0,ioutsf2outpt(1:3,i),ioutsf2surf(i)
        enddo
        close(iot1)
        open(iot1,file='Tetra.otpt',status='old')
        do i=1,noutpt(nlat3)
         read(iot1,5000)i0,ioutpt2point(i)
        enddo
        close(iot1)
        open(iot1,file='Tetra.octO',status='old')
        read(iot1,5000)(ioctOsec(i),i=1,noctOsec(nlat3))
        close(iot1)
        open(iot1,file='Tetra.octE',status='old')
        read(iot1,5000)(ioctEblk(i),i=1,noctEblk(nlat3))
        close(iot1)
        open(iot1,file='Tetra.indO',status='old')
        read(iot1,5000)(indmapOsec(i,1),i=1,n2o)
        read(iot1,5000)(indmapOsec(i,2),i=1,n2o)
        close(iot1)
        open(iot1,file='Tetra.indE',status='old')
        read(iot1,5000)(indmapEblk(i,1),i=1,n2e)
        read(iot1,5000)(indmapEblk(i,2),i=1,n2e)
        close(iot1)
        open(iot1,file='Tetra.mapO',status='old')
        read(iot1,5000)(ioutsfmapOsec(i),i=1,n1o)
        close(iot1)
        open(iot1,file='Tetra.mapE',status='old')
        read(iot1,5000)(ielemmapEblk(i),i=1,n1e)
        close(iot1)
      else
        inquire(file=binfile, exist = exex)
        if(.not.exex)then
         write(*,*)'*** TETRA ERROR: ',
     &        'binary file: ',binfile(1:lbin),' NOT EXIST'
         call parastop(301)
        endif
        open(iot1,file=binfile,form='unformatted',status='old')
c-----------------------------------------------------------------------
        read(iot1)iflag,ntetsurf,ntetelem
        if(iflag.ne.0)then
         write(*,*)'*** TETRA ERROR: ',
     &        'inconsistent file format'
         write(*,*)'Binary file should be re-created ',
     &        'with itetra=2'
         call parastop(301)
        endif
c-----------------------------------------------------------------------
        do itet=1,nlat3
         read(iot1)npoint(itet),nelem(itet),nsurf(itet),
     &        noutpt(itet),noutsf(itet),noctOsec(itet),noctEblk(itet)
        enddo
        do itet=1,nlat3
         read(iot1)xb(1:6,itet)
        enddo
        do itet=1,nlat3
         do idiv=1,ndiv
          read(iot1)xoctOsec(1:3,idiv,itet),xoctEblk(1:3,idiv,itet)
         enddo
        enddo
        read(iot1)n1o,n2o,n1e,n2e
        call allocate_tetratbl1(nlat3)
        call allocate_tetratbl2(nlat3)
        noutpt(0)=0 !FURUTA20221220
        noutsf(0)=0 !FURUTA20221220
        noctOsec(0)=0
        noctEblk(0)=0
        do itet=1,nlat3
         noutpt(itet)=noutpt(itet)+npoint(itet-1) !FURUTA20221220
         noutsf(itet)=noutsf(itet)+nsurf(itet-1)  !FURUTA20221220
         noctOsec(itet)=noctOsec(itet)+noctOsec(itet-1)
         noctEblk(itet)=noctEblk(itet)+noctEblk(itet-1)
        enddo
        call allocate_blocktbl(n1e,n2e,noctEblk(nlat3))
        call allocate_outertbl(n1o,n2o,noctOsec(nlat3))
        do i=1,npoint(nlat3)
         read(iot1)pointxyz(1:3,i)
        enddo
        do i=1,nsurf(nlat3)
         read(iot1)isurf2elem(1:2,i),surfcoefs(1:4,i)
        enddo
        do i=1,nelem(nlat3)
         read(iot1)ielem2id(i),ielem2point(1:4,i),
     &        ielem2surf(1:4,i),ielem2univ(i)
        enddo
        do i=1,noutsf(nlat3)
         read(iot1)ioutsf2outpt(1:3,i),ioutsf2surf(i)
        enddo
        do i=1,noutpt(nlat3)
         read(iot1)ioutpt2point(i)
        enddo
        read(iot1)(ioctOsec(i),i=1,noctOsec(nlat3))
        read(iot1)(ioctEblk(i),i=1,noctEblk(nlat3))
        read(iot1)(indmapOsec(i,1),i=1,n2o)
        read(iot1)(indmapOsec(i,2),i=1,n2o)
        read(iot1)(indmapEblk(i,1),i=1,n2e)
        read(iot1)(indmapEblk(i,2),i=1,n2e)
        read(iot1)(ioutsfmapOsec(i),i=1,n1o)
        read(iot1)(ielemmapEblk(i),i=1,n1e)
        close(iot1)
      endif
      return
      end subroutine read_infofile

************************************************************************
      subroutine read_infofile2(iot1,nlat3,nunivtot,nmats)
*
*     Read info from file for itetauto=1
*
*     Created 2019/03/31 by T. Furuta
*
************************************************************************
      integer,intent(in) :: iot1,nlat3
      integer,intent(out) :: nunivtot,nmats(0:10)
      logical :: exex
      integer :: itet,i
      nmats(0)=0
      if(itetra.lt.0)then
       open(iot1,file='Tetra.univ',status='unknown')
       read(iot1,'(i9)')nunivtot
       do itet=1,nlat3
        read(iot1,5000)i,nmats(itet)
       enddo
 5000  format(i9,' ',i9)
      else
       inquire(file=binfile2,exist=exex)
       if(.not.exex)then
        write(*,*)'*** TETRA ERROR: ',
     &       'binary file: ',binfile2(1:lbin+1),' NOT EXIST'
        call parastop(301)
       endif
       open(iot1,file=binfile2,form='unformatted',status='old')
       read(iot1)nunivtot,nmats(1:nlat3)
      endif
      return
      end subroutine read_infofile2

************************************************************************
      subroutine read_infofile3(iot1,nlat3,nunivtot,nmats,
     &     iunivmat,materials,densuniv)
*
*     Read info from file for itetauto=1
*
*     Created 2019/03/31 by T. Furuta
*
************************************************************************
      integer,intent(in) :: nlat3,nunivtot,iot1
      integer,intent(in) :: nmats(0:10)
      integer,intent(out) :: iunivmat(nunivtot),materials(nunivtot)
      real(8),intent(out) :: densuniv(nunivtot)
      integer :: i,i0
      if(itetra.lt.0)then
       do i=1,nunivtot
        read(iot1,5000)i0,iunivmat(i),densuniv(i)
       enddo
       do i=1,nmats(nlat3)
        read(iot1,5001)i0,materials(i)
       enddo
 5000  format(i9,' ',i9,' ',1pe22.15)
 5001  format(i9,' ',i9)
       close(iot1)
      else
       do i=1,nunivtot
        read(iot1)iunivmat(i)
       enddo
       do i=1,nunivtot
        read(iot1)densuniv(i)
       enddo
       do i=1,nmats(nlat3)
        read(iot1)materials(i)
       enddo
       close(iot1)
      endif
      return
      end subroutine read_infofile3

************************************************************************
      subroutine read_tetraparam0(itet,iot1,iot2,ierr)
*
*     Read first lines of tetra files (.node and .ele)
*
************************************************************************
      include 'err.inc'

      integer,intent(in) :: itet,iot1,iot2
      integer,intent(out) :: ierr
      integer :: is,ip,ntag
      character(200) :: chin
      integer :: ios,i1,i2
*------------------------------------------------------------------------
*     Read .node parameter
*------------------------------------------------------------------------
      do
       read(iot1,'(a200)',iostat=ios)chin
       if(ios.eq.-1)then
        write(ErrCha,'(''*** TETRA ERROR: '',
     &      ''.node file is ended before reading parameter'')')
         ErrID = 'L:701/R:read_tetraparam0/F:tetramod.f' !E80_003_001
         call ErrWrite(ErrID,ErrCha)
        ierr=1
        exit
       endif
       call chlngt(chin,200,i1,i2)
       if(chin(i1:i1).ne.'#')exit
      enddo
      if(ierr.ne.0)return
      read(chin,*,iostat=ios)npoint(itet),is
      if(ios.gt.0)then

       ErrCha = ''
       ErrID = 'L:714/R:read_tetraparam0/F:tetramod.f' !E80_004_001
       call ErrWrite(ErrID,ErrCha)

       write(*,'(''*** TETRA ERROR: '',
     &      ''.node file format is not compatible'')')
       write(*,'(a200)')chin(i1:i2)

      endif
      if(is.ne.3)then

       ErrCha = ''
       ErrID = 'L:725/R:read_tetraparam0/F:tetramod.f' !E80_004_004
       call ErrWrite(ErrID,ErrCha)

       write(*,'(''*** TETRA ERROR: '',
     &      ''.node file format is not compatible'')')
       write(*,'(''npoint ='',i5,'' is ='',i3)')npoint(itet),is
       ierr=1
      endif
*------------------------------------------------------------------------
*     Read .ele parameter
*------------------------------------------------------------------------
      do
       read(iot2,'(a200)',iostat=ios)chin
       if(ios.eq.-1)then
        write(ErrCha,'(''*** TETRA ERROR: '',
     &      ''.ele file is ended before reading parameter'')')
         ErrID = 'L:741/R:read_tetraparam0/F:tetramod.f' !E80_005_001
         call ErrWrite(ErrID,ErrCha)

        ierr=1
        exit
       endif
       call chlngt(chin,200,i1,i2)
       if(chin(i1:i1).ne.'#')exit
      enddo
      if(ierr.ne.0)return
      read(chin,*,iostat=ios)nelem(itet),ip,ntag
      if(ios.gt.0)then

       ErrCha = ''
       ErrID = 'L:755/R:read_tetraparam0/F:tetramod.f' !E80_005_003
       call ErrWrite(ErrID,ErrCha)

       write(*,'(''*** TETRA ERROR: '',
     &      ''.ele file format is not compatible'')')
       write(*,'(a200)')chin(i1:i2)
      endif
      if(ip.ne.4)then

       ErrCha = ''
       ErrID = 'L:765/R:read_tetraparam0/F:tetramod.f' !E80_005_004
       call ErrWrite(ErrID,ErrCha)

       write(*,'(''*** TETRA ERROR: '',
     &      ''.ele file format is not compatible'')')
       write(*,'(''nelem ='',i5,'' ip ='',i3)')nelem(itet),ip
       ierr=1
      endif
      if(ntag.lt.1)then

       ErrCha = ''
       ErrID = 'L:776/R:read_tetraparam0/F:tetramod.f' !E80_007_001
       call ErrWrite(ErrID,ErrCha)

       write(*,'(''*** TETRA ERROR: '',
     &      ''universe info is missing'')')
       write(*,'(''ntag ='',i5)')ntag
       ierr=1
      else if(ntag.gt.1)then
       if(iwarning)then

        ErrCha = ''
        ErrID = 'L:787/R:read_tetraparam0/F:tetramod.f' !W80_001_001
        call ErrWrite(ErrID,ErrCha)

        write(*,'(''*** TETRA WARNING: '',
     &      ''extra info is not used'')')
        write(*,'(''ntag ='',i5,'' more than 1'')')ntag
       endif
      endif
      return
      end subroutine read_tetraparam0

************************************************************************
      subroutine read_tetraparam1(itet,iot,ierr)
*
*     Read parameters of tetra files (BDF)
*
*     Created 2019/02/08 by T. Furuta
*
************************************************************************
      include 'err.inc'

      integer,intent(in) :: itet,iot
      integer,intent(out) :: ierr
      integer :: it,ip
      integer :: ios
      character(80) :: dummy
*------------------------------------------------------------------------
*     Read BDF parameter
*------------------------------------------------------------------------
      ip=0
      it=0
      do
       read(iot,'(a80)',iostat=ios)dummy
       if(ios.eq.-1)then
        write(ErrCha,'(''*** TETRA ERROR: '',
     &      ''BDF file is ended before reading parameter'')')
         ErrID = 'L:823/R:read_tetraparam1/F:tetramod.f' !E80_003_002
         call ErrWrite(ErrID,ErrCha)
        ierr=1
        exit
       endif
       if(dummy(1:4).eq.'GRID')ip=ip+1
       if(dummy(1:6).eq.'CTETRA')it=it+1
       if(dummy(1:7).eq.'ENDDATA')exit
      enddo
      if(ierr.ne.0)return
      if(ip.eq.0.or.it.eq.0)then
       write(ErrCha,'(''*** TETRA ERROR: '',
     &      ''BDF file format is not compatible'')')
         ErrID = 'L:836/R:read_tetraparam1/F:tetramod.f' !E80_003_003
         call ErrWrite(ErrID,ErrCha)
       ierr=1
       return
      endif
      npoint(itet)=ip
      nelem(itet)=it
      return
      end subroutine read_tetraparam1

************************************************************************
      subroutine read_tetraparam2(itet,filename,ierr)
*
*     Interface of subroutine to read parameters of tetra file (hdf5)
*
*     Created 2024/07/09 by T. Furuta
*
************************************************************************
      implicit none
      integer,intent(in) :: itet
      character(200),intent(in) :: filename
      integer,intent(out) :: ierr
      integer :: melem,mpoint

      call read_HDF5tetraparam(filename,melem,mpoint,ierr)
      nelem(itet)=melem
      npoint(itet)=mpoint
      return
      end subroutine read_tetraparam2

************************************************************************
      subroutine read_tetradata0(tsfac,itet,iot1,iot2,ierr)
*
*     Read tetra files (.node and .ele) and set pointxyz, ielem2point
*
************************************************************************
      include 'err.inc'

      integer,intent(in) :: itet,iot1,iot2
      real(8),intent(in) :: tsfac
      integer,intent(out) :: ierr
      integer :: id,is,ip,it,i,ipp(4)
      integer,allocatable :: id2ipt(:) !FURUTA20221104
      character(200) :: chin
      integer :: ios,i1,i2
*------------------------------------------------------------------------
*     Read .node data
*------------------------------------------------------------------------
      allocate( id2ipt(0:99999999) )  !FURUTA20221107
      do i=npoint(itet-1)+1-1,npoint(itet)
       do
        read(iot1,'(a200)',iostat=ios)chin
        if(ios.eq.-1)then
         write(ErrCha,'(''*** TETRA ERROR: '',
     &      ''.node file is ended before reading parameter'')')
         ErrID = 'L:891/R:read_tetradata0/F:tetramod.f' !E80_004_003
         call ErrWrite(ErrID,ErrCha)

         ierr=1
         exit
        endif
        call chlngt(chin,200,i1,i2)
        if(chin(i1:i1).ne.'#')exit
       enddo
       if(ierr.ne.0)return
       if(i.gt.npoint(itet-1))then
        read(chin,*)id,(pointxyz(is,i),is=1,3)
        if(id.gt.99999999)then !FURUTA20221107
         write(ErrCha,'(''*** TETRA ERROR: GRID ID '',i8,
     &        '' is larger than 99999999'')') id
         ErrID = 'L:906/R:read_tetradata0/F:tetramod.f' !E80_005_002
         call ErrWrite(ErrID,ErrCha)
         ierr=1
         exit
        endif
        id2ipt(id)=i !FURUTA20221104
       endif
      enddo
      if(tsfac.ne.1.0d0)then
       do i=npoint(itet-1)+1,npoint(itet)
        pointxyz(1:3,i)=tsfac*pointxyz(1:3,i)
       enddo
      endif
*------------------------------------------------------------------------
*     Read .ele data
*------------------------------------------------------------------------
      do i=nelem(itet-1)+1-1,nelem(itet)
       do
        read(iot2,'(a200)',iostat=ios)chin
        if(ios.eq.-1)then
         write(ErrCha,'(''*** TETRA ERROR: '',
     &      ''.ele file is ended before reading parameter'')')
         ErrID = 'L:928/R:read_tetradata0/F:tetramod.f' !E80_005_002
         call ErrWrite(ErrID,ErrCha)

         ierr=1
         exit
        endif
        call chlngt(chin,200,i1,i2)
        if(chin(i1:i1).ne.'#')exit
       enddo
       if(i.gt.nelem(itet-1))then
        read(chin,*)ielem2id(i),(ipp(ip),ip=1,4),ielem2univ(i) !FURUTA20190110
        do ip=1,4                          !FURUTA20221104
         ielem2point(ip,i)=id2ipt(ipp(ip)) !FURUTA20221104
        enddo                              !FURUTA20221104
       endif
      enddo
      deallocate( id2ipt ) !FURUTA20221104
      return
      end subroutine read_tetradata0

************************************************************************
      subroutine read_tetradata1(tsfac,itet,iot,ierr)
*
*     Read tetra files (BDF) and set pointxyz, ielem2point
*
*     Last Modified 2022/11/04 by T. Furuta
*
************************************************************************
      include 'err.inc'
      integer,intent(in) :: itet,iot
      real(8),intent(in) :: tsfac
      integer,intent(out) :: ierr
      integer :: id,ip,it,i,ipp(4),is,ind,iend
      integer :: iform !FURUTA20221104
      integer :: ios
      integer,allocatable :: id2ipt(:) !FURUTA20221104
      real(8) :: factor
      character(200) :: dummy,chin
      character(1) :: ccheck
*------------------------------------------------------------------------
*     Read GRID data
*------------------------------------------------------------------------
      allocate( id2ipt(0:99999999) )
      do
       write(dummy,'(200x)')
       read(iot,'(a200)',iostat=ios)dummy
       if(ios.eq.-1)then
        write(ErrCha,'(''*** TETRA ERROR: '',
     &       ''BDF file is ended before reading GRID'')')
        ErrID = 'L:977/R:read_tetradata1/F:tetramod.f' !E80_004_004
        call ErrWrite(ErrID,ErrCha)
        ierr=1
        return
       endif
       if(dummy(1:4).eq.'GRID')exit
      enddo
      ccheck=' '
      iform=0
      do i=5,8
       if(dummy(i:i).eq.'*')then
        ccheck='*'
        iform=1
        exit
       endif
      enddo
      do i=npoint(itet-1)+1,npoint(itet)
       do
        if(dummy(1:4).eq.'GRID')exit
        write(dummy,'(200x)')
        read(iot,'(a200)',iostat=ios)dummy
        if(ios.eq.-1)then
         write(ErrCha,'(''*** TETRA ERROR: '',
     &        ''BDF file is ended during reading GRID'')')
         ErrID = 'L:1001/R:read_tetradata1/F:tetramod.f' !E80_004_005
         call ErrWrite(ErrID,ErrCha)
         ierr=1
         return
        endif
       enddo
       write(chin,'(200x)')
       ind=72
       chin(1:ind)=dummy(1:ind)
       iend=ind
       write(dummy,'(200x)')
       if(ccheck.eq.'*')then
        do
         write(dummy,'(200x)')
         read(iot,'(a200)',iostat=ios)dummy
         if(ios.eq.-1)then
          if(i.eq.npoint(itet))exit
          write(ErrCha,'(''*** TETRA ERROR: '',
     &         ''BDF file is ended during reading GRID'')')
          ErrID = 'L:1020/R:read_tetradata1/F:tetramod.f' !E80_004_005
          call ErrWrite(ErrID,ErrCha)
          ierr=1
          return
         endif
         if(dummy(1:1).eq.ccheck)then
          ind=len_trim(dummy)
          if(ind.gt.200-iend)ind=200-iend
          if(ind.gt.0)then
           chin(iend+1:iend+ind-8)=dummy(1+8:ind) ! Remove head 8 bytes
           iend=iend+ind
          endif
         else
          exit
         endif
        enddo
       endif
       if(iform.eq.0)then
        read(chin,5000,iostat=ios)id,(pointxyz(is,i),is=1,3)
       else
        read(chin,5001,iostat=ios)id,(pointxyz(is,i),is=1,3)
       endif
       if(id.gt.99999999)then   !FURUTA20221107
        write(ErrCha,'(''*** TETRA ERROR: GRID ID '',i8,
     &       '' is larger than 99999999'')') id
        ErrID = 'L:1045/R:read_tetradata1/F:tetramod.f' !E80_005_002
        call ErrWrite(ErrID,ErrCha)
        ierr=1
        exit
       endif
       id2ipt(id)=i
      enddo
      factor=tsfac*100.0d0 ! m -> cm
      do i=npoint(itet-1)+1,npoint(itet)
       pointxyz(1:3,i)=factor*pointxyz(1:3,i)
      enddo
*------------------------------------------------------------------------
*     Read CTETRA data
*------------------------------------------------------------------------
      rewind(iot)
      do
       write(dummy,'(200x)')
       read(iot,'(a200)',iostat=ios)dummy
       if(ios.eq.-1)then
        write(ErrCha,'(''*** TETRA ERROR: '',
     &       ''BDF file is ended before reading CTETRA'')')
        ErrID = 'L:1066/R:read_tetradata1/F:tetramod.f' !E80_005_005
        call ErrWrite(ErrID,ErrCha)
        ierr=1
        return
       endif
       if(dummy(1:6).eq.'CTETRA')exit
      enddo
      ccheck=' '
      iform=0
      do i=7,8
       if(dummy(i:i).eq.'*')then
        ccheck='*'
        iform=1
        exit
       endif
      enddo
      do i=nelem(itet-1)+1,nelem(itet)
       do
        if(dummy(1:6).eq.'CTETRA')exit
        write(dummy,'(200x)')
        read(iot,'(a200)',iostat=ios)dummy
        if(ios.ne.0)then
         write(ErrCha,'(''*** TETRA ERROR: '',
     &        ''BDF file is ended during reading CTETRA'')')
         ErrID = 'L:1090/R:read_tetradata1/F:tetramod.f' !E80_005_006
         call ErrWrite(ErrID,ErrCha)
         ierr=1
         return
        endif
       enddo
       write(chin,'(200x)')
       ind=72
       chin(1:ind)=dummy(1:ind)
       iend=ind
       write(dummy,'(200x)')
       if(ccheck.eq.'*')then
        do
         write(dummy,'(200x)')
         read(iot,'(a200)',iostat=ios)dummy
         if(ios.eq.-1)then
          write(ErrCha,'(''*** TETRA ERROR: '',
     &         ''BDF file is ended during reading CTETRA'')')
          ErrID = 'L:1108/R:read_tetradata1/F:tetramod.f' !E80_005_006
          call ErrWrite(ErrID,ErrCha)
          ierr=1
          return
         endif
         if(dummy(1:1).eq.ccheck)then
          ind=len_trim(dummy)
          if(ind.gt.200-iend)ind=200-iend
          if(ind.gt.0)then
           chin(iend+1:iend+ind-8)=dummy(1+8:ind) ! Remove head 8 bytes
           iend=iend+ind
          endif
         else
          exit
         endif
        enddo
       endif
       if(iform.eq.0)then
        read(chin,5002)ielem2id(i),ielem2univ(i),(ipp(ip),ip=1,4)
       else
        read(chin,5003)ielem2id(i),ielem2univ(i),(ipp(ip),ip=1,4)
       endif
       do ip=1,4
        ielem2point(ip,i)=id2ipt(ipp(ip))
       enddo
      enddo
 5000 format(8x,i8,8x,3es8.0)
 5001 format(8x, i16, 16x, 3es16.0)
 5002 format(8x, 6i8)
 5003 format(8x, 6i16)
      deallocate( id2ipt )
      return
      end subroutine read_tetradata1

************************************************************************
      subroutine read_tetradata2(tsfac,itet,filename,ierr)
*
*     Interface of subroutine to read tetra file (hdf5) data
*     ielem2univ, ielem2point, pointxyz
*
*     Created 2024/07/09 by T. Furuta
*
************************************************************************
      integer,intent(in) :: itet
      real(8),intent(in) :: tsfac
      character(200),intent(in) :: filename
      integer,intent(out) :: ierr
      integer,allocatable :: ie2u(:),ie2p(:,:)
      real(8),allocatable :: pxyz(:,:)
      integer :: i,melem,mpoint
      melem=nelem(itet)-nelem(itet-1)
      mpoint=npoint(itet)-npoint(itet-1)
      allocate( ie2u(melem),ie2p(4,melem),pxyz(3,mpoint) )
      call read_HDF5tetradata(filename,melem,mpoint,ie2u,ie2p,pxyz,ierr)
      do i=1,melem
       ielem2id(nelem(itet-1)+i)=i
       ielem2univ(nelem(itet-1)+i)=ie2u(i)
       ielem2point(1:4,nelem(itet-1)+i)=npoint(itet-1)+ie2p(1:4,i)
      enddo
      do i=1,mpoint
       pointxyz(1:3,npoint(itet-1)+i)=tsfac*pxyz(1:3,i)
      enddo
      deallocate( ie2u,ie2p,pxyz )
      return
      end subroutine read_tetradata2

************************************************************************
      subroutine read_tetratxt0(itet,iot,kvlmax,nunivtot,
     &     nuniv,iorguniv,nmats,iunivmat,materials,densuniv,ierr)
*
*     Read tetra file (.txt) and set densuniv, iunivmat
*
*     Created by T. Furuta on 2019/02/08
*
************************************************************************
      implicit real*8 (a-h,o-z)
      include 'err.inc'
      integer,intent(in) :: itet,iot,kvlmax,nunivtot
      integer,intent(in) :: nuniv(0:10),iorguniv(kvlmax)
      integer,intent(inout) :: nmats(0:10)
      integer,intent(out) :: iunivmat(nunivtot),materials(nunivtot)
      real(8),intent(out) :: densuniv(nunivtot)
      integer,intent(out) :: ierr
      integer :: i,j,k,i1,i2,ios,numuniv,nummat,iuniv,imat
      real(8) :: dens
      character(200) :: chin,cname
      logical :: iflag

*-----------------------------------------------------------------------
      numuniv=nuniv(itet)-nuniv(itet-1)
      nummat=0
      do i=1,numuniv
       do
        read(iot,'(a200)',iostat=ios)chin
        if(ios.eq.-1)then
          write(ErrCha,'(''*** TETRA ERROR: '',
     &         ''.txt file is ended before reading parameter'')')
          ErrID = 'L:1205/R:read_tetratxt0/F:tetramod.f' !E80_006_003
          call ErrWrite(ErrID,ErrCha)
          ierr=1
          return
        endif
        call chlngt(chin,200,i1,i2)
        if(chin(i1:i1).ne.'#')exit
       enddo
       read(chin,*)iuniv,imat,dens
       iflag=.false.
       do k=1,nummat
        if(materials(nmats(itet-1)+k).eq.imat)iflag=.true.
       enddo
       if(.not.iflag)then
        nummat=nummat+1
        materials(nmats(itet-1)+nummat)=imat
       endif
       iflag=.false.
       do j=nuniv(itet-1)+1,nuniv(itet)
        if(iorguniv(j).eq.iuniv)then
         iflag=.true.
         densuniv(j)=dens
         iunivmat(j)=nmats(itet-1)+nummat
        endif
       enddo
       if(.not.iflag)then
        write(ErrCha,'(''*** TETRA ERROR: univ ID '',i8,
     &       '' is inconsistent between .txt and .ele'')')
     &           iuniv
        ErrID = 'L:1234/R:read_tetratxt0/F:tetramod.f' !E80_006_004
        call ErrWrite(ErrID,ErrCha)
        ierr=1
        return
       endif
      enddo
      nmats(itet)=nmats(itet-1)+nummat
      return
      end subroutine read_tetratxt0

************************************************************************
      subroutine read_tetratxt1(itet,iot,kvlmax,nunivtot,
     &     nuniv,iorguniv,nmats,iunivmat,materials,densuniv,ierr)
*
*     Read tetra file (BDF) and set densuniv, iunivmat
*
*     Created by T. Furuta on 2019/02/08
*
************************************************************************
      include 'err.inc'
      integer,intent(in) :: itet,iot,kvlmax,nunivtot
      integer,intent(in) :: nuniv(0:10),iorguniv(kvlmax)
      integer,intent(inout) :: nmats(0:10)
      integer,intent(out) :: iunivmat(nunivtot),materials(nunivtot)
      real(8),intent(out) :: densuniv(nunivtot)
      integer,intent(out) :: ierr
      integer :: i,j,k,ios,numuniv,nummat,iuniv,imat,iend,ind

      real(8) rdummy(3)

      real(8) dens
      character(200) :: dummy,chin
      character(1) :: ccheck
      integer,allocatable :: iflag(:)
      logical :: jflag
*------------------------------------------------------------------------
*     Read PSOLID data
*------------------------------------------------------------------------
      numuniv=nuniv(itet)-nuniv(itet-1)
      nummat=0
      do i=1,numuniv
       do
        write(dummy,'(200x)')
        read(iot,'(a200)',iostat=ios)dummy
        if(ios.eq.-1)then
         write(ErrCha,'(''*** TETRA ERROR: '',
     &        ''BDF file is ended before reading PSOLID'')')
         ErrID = 'L:1281/R:read_tetratxt1/F:tetramod.f' !E80_006_005
         call ErrWrite(ErrID,ErrCha)
         ierr=1
         return
        endif
        if(dummy(1:6).eq.'PSOLID')exit
       enddo
       read(dummy,5000,iostat=ios)iuniv,imat
       jflag=.false.
       do k=1,nummat
        if(materials(nmats(itet-1)+k).eq.imat)jflag=.true.
       enddo
       if(.not.jflag)then
        nummat=nummat+1
        materials(nmats(itet-1)+nummat)=imat
       endif
       jflag=.false.
       do j=nuniv(itet-1)+1,nuniv(itet)
        if(iorguniv(j).eq.iuniv)then
         jflag=.true.
         iunivmat(j)=nmats(itet-1)+nummat
        endif
       enddo
       if(.not.jflag)then
        write(ErrCha,'(''*** TETRA ERROR: univ ID '',i8,
     &       '' is inconsistent between PSOLID and CTETRA'')')
     &           iuniv
        ErrID = 'L:1308/R:read_tetratxt1/F:tetramod.f' !E80_006_006
        call ErrWrite(ErrID,ErrCha)
        ierr=1
        return
       endif
      enddo
      nmats(itet)=nmats(itet-1)+nummat
*------------------------------------------------------------------------
*     Read density data
*------------------------------------------------------------------------
      allocate(iflag(numuniv))
      iflag(1:numuniv)=0
      do
       write(dummy,'(200x)')
       read(iot,'(a200)',iostat=ios)dummy
       if(ios.eq.-1)then
        write(ErrCha,'(''*** TETRA ERROR: '',
     &       ''BDF file is ended before reading parameter'')')
        ErrID = 'L:1326/R:read_tetratxt1/F:tetramod.f' !E80_006_007
        call ErrWrite(ErrID,ErrCha)
        ierr=1
        return
       endif
       if(dummy(1:3).eq.'MAT')exit
      enddo
      do i=1,nummat
       ccheck=dummy(5:5)
       write(chin,'(200x)')
       ind=len_trim(dummy)
       chin(1:ind)=dummy(1:ind)
       iend=ind
       write(dummy,'(200x)')
       read(iot,'(a200)',iostat=ios)dummy
       if(ios.eq.-1)then
        if(sum(iflag(1:numuniv)).ne.numuniv)then
         write(ErrCha,'(''*** TETRA ERROR: '',
     &        ''BDF file is ended during reading MAT'')')
         ErrID = 'L:1345/R:read_tetratxt1/F:tetramod.f' !E80_006_008
         call ErrWrite(ErrID,ErrCha)
         ierr=1
         return
        else
         ccheck=' '
        endif
       endif
       if(ccheck.ne.' ')then
        do
         if(dummy(1:3).eq.'MAT')then
          exit
         elseif(dummy(1:1).eq.ccheck)then
          ind=len_trim(dummy)
          if(ind.gt.200-iend)ind=200-iend
          if(ind.gt.0)then
           chin(iend+1:iend+ind-8)=dummy(1+8:ind) ! Remove head 8 bytes
           iend=iend+ind
          endif
         elseif(i.eq.nummat)then
          exit
         endif
         write(dummy,'(200x)')
         read(iot,'(a200)',iostat=ios)dummy
         if(ios.eq.-1)then
          if(sum(iflag(1:numuniv)).eq.numuniv)exit
          write(ErrCha,'(''*** TETRA ERROR: '',
     &         ''BDF file is ended during reading MAT'')')
          ErrID = 'L:1373/R:read_tetratxt1/F:tetramod.f' !E80_006_009
          call ErrWrite(ErrID,ErrCha)
          ierr=1
          return
         endif
        enddo
       endif
       read(chin,5001,iostat=ios)imat,dens
       jflag=.false.
       do k=1,nummat
        if(materials(nmats(itet-1)+k).eq.imat)then
         jflag=.true.
         exit
        endif
       enddo
       if(.not.jflag)then
        write(ErrCha,'(''*** TETRA ERROR: material ID '',i8,
     &       '' is inconsistent between MAT and PSOLID'')')
     &           imat
        ErrID = 'L:1392/R:read_tetratxt1/F:tetramod.f' !E80_006_010
        call ErrWrite(ErrID,ErrCha)
        ierr=1
        return
       endif
       do j=nuniv(itet-1)+1,nuniv(itet)
        if(iunivmat(j).eq.nmats(itet-1)+k)then
         iflag(j-nuniv(itet-1))=1
         densuniv(j)=-dens*1.0d-3 ! kg/m3 -> g/cm3 !FURUTA20190607 fix
        endif
       enddo
      enddo
      if(sum(iflag(1:numuniv)).ne.numuniv)then
       write(ErrCha,'(''*** TETRA ERROR: '',
     &      ''BDF file is ended during reading MAT'')')
       ErrID = 'L:1407/R:read_tetratxt1/F:tetramod.f' !E80_006_009
       call ErrWrite(ErrID,ErrCha)
       ierr=1
       return
      endif
      deallocate( iflag )

 5000 format(8x, 2i8)
 5001 format(16x, i8, 3(16x) es16.0)
      return
      end subroutine read_tetratxt1

************************************************************************
      subroutine read_tetratxt2(itet,filename,nunivtot,nuniv,
     &     nmats,iunivmat,materials,densuniv,ierr)
*
*     Read tetra file (HDF) and set densuniv, iunivmat
*
*     Created by T. Furuta on 2019/02/08
*
************************************************************************
      include 'err.inc'
      integer,intent(in) :: itet
      character(200),intent(in) :: filename
      integer,intent(in) :: nunivtot,nuniv(0:10)
      integer,intent(inout) :: nmats(0:10)
      integer,intent(out) :: iunivmat(nunivtot),materials(nunivtot)
      real(8),intent(out) :: densuniv(nunivtot)
      integer,intent(out) :: ierr
      integer,allocatable :: iu2m(:)
      real(8),allocatable :: du2d(:)
      integer :: i,nu,nummat

      nu=nuniv(itet)-nuniv(itet-1)
      allocate( iu2m(nu),du2d(nu))

      call read_HDF5tetramat(filename,nu,nmats(itet-1),
     &     iu2m,du2d,nummat,ierr)
      do i=1,nu
       iunivmat(nuniv(itet-1)+i)=nmats(itet-1)+iu2m(i)
       densuniv(nuniv(itet-1)+i)=du2d(i)
      enddo
      nmats(itet)=nmats(itet-1)+nummat
      do i=1,nummat
       materials(nmats(itet-1)+i)=i
      enddo

      deallocate( iu2m,du2d )

      end subroutine read_tetratxt2

************************************************************************
      subroutine set_xquota(nlat3)
*
*     Set small quota to avoid round error
*
*     Last modified 2017/10/20 by T. Furuta
*
************************************************************************
      integer,intent(in) :: nlat3
      integer :: itet,is
      real(8) :: xsum
      do itet=1,nlat3
       xsum=0.0d0
       do is=1,3
        xsum=xsum+xoctEblk(is,ndiv,itet)*xoctEblk(is,ndiv,itet)
       enddo
       xqt(itet)=sqrt(xsum)
      enddo
      return
      end subroutine set_xquota

************************************************************************
      subroutine check_inside(xx0,itet,ierr)
*
*     Check all tetrahedron points are inside of box
*
*     Last modified 2016/09/02 by T. Furuta
*
************************************************************************
      include 'err.inc'

      real(8),intent(in) :: xx0(6)
      integer,intent(in) :: itet
      integer,intent(inout) :: ierr
      integer :: i,is,isig,idiv
      real(8) :: xx(6)
      xx(1)=1.0d20
      xx(2)=-1.0d20
      xx(3)=1.0d20
      xx(4)=-1.0d20
      xx(5)=1.0d20
      xx(6)=-1.0d20
      do i=npoint(itet-1)+1,npoint(itet)
       do is=1,3
        xx(2*is-1)=min(xx(2*is-1),pointxyz(is,i))
        xx(2*is)=max(xx(2*is),pointxyz(is,i))
        if(pointxyz(is,i).lt.xx0(2*is-1).or.
     &       pointxyz(is,i).gt.xx0(2*is))then
         ierr=1
        endif
       enddo
      enddo
      if(ierr.eq.1)then

       ErrCha = ''
       ErrID = 'L:1513/R:check_inside/F:tetramod.f' !E80_008_001
       call ErrWrite(ErrID,ErrCha)

       write(*,'(''*** TETRA ERROR: '',
     &      ''tetra node is not inside box'')')
       write(*,'(''No. '',i3,
     &      '' tetra box should be larger than '',6g9.2)')
     &      itet,xx(1:6)
      endif
      do is=1,3
       xb(2*is-1,itet)=xx(2*is-1)-abs(xx(2*is-1))*xtiny
       xb(2*is,itet)=xx(2*is)+abs(xx(2*is))*xtiny
      enddo
*------------------------------------------------------------------------
*     Define octree section distance !FURUTA20190331 move to here
*------------------------------------------------------------------------
      do is=1,3
       xoctOsec(is,1,itet)=0.5d0*(xx0(is*2)-xx0(is*2-1))
      enddo
      do idiv=2,ndiv
       xoctOsec(1:3,idiv,itet)=0.5d0*xoctOsec(1:3,idiv-1,itet)
      enddo
*------------------------------------------------------------------------
*     Define octree block distance !FURUTA20190331 move to here
*------------------------------------------------------------------------
      do is=1,3
       xoctEblk(is,1,itet)=0.5d0*(xb(is*2,itet)-xb(is*2-1,itet))
      enddo
      do idiv=2,ndiv
       xoctEblk(1:3,idiv,itet)=0.5d0*xoctEblk(1:3,idiv-1,itet)
      enddo
*------------------------------------------------------------------------
      return
      end subroutine check_inside

************************************************************************
      subroutine reorder_ielem2point(nlat3)
*
*     Reorder point index of ielem2point in ascending order
*
************************************************************************
      integer,intent(in) :: nlat3
      integer :: ip(4),itmp
      integer :: ielem,i,j
      do ielem=1,nelem(nlat3)
       ip(1:4)=ielem2point(1:4,ielem)
       do j=1,3
        do i=1,4-j
         if(ip(i).gt.ip(i+1))then
          itmp=ip(i)
          ip(i)=ip(i+1)
          ip(i+1)=itmp
         endif
        enddo
       enddo
       ielem2point(1:4,ielem)=ip(1:4)
      enddo
      return
      end subroutine reorder_ielem2point

************************************************************************
      subroutine set_nsurf(itet)
*
*     Setup total number of surfaces of all elements (nsurf)
*
************************************************************************
      integer,intent(in) :: itet
      integer :: i,j,k,l,n,ip(4),isfcount,ielem
      integer :: isold(3),ne0
      l=4*nelem(itet-1)
      ne0=nelem(itet)-nelem(itet-1)
      do ielem=nelem(itet-1)+1,nelem(itet)
       ip(1:4)=ielem2point(1:4,ielem)
       n=5
       do i=1,2
        do j=i+1,3
         do k=j+1,4
          l=l+1
          n=n-1
          isfind(1,l)=ip(i)
          isfind(2,l)=ip(j)
          isfind(3,l)=ip(k)
          isfind(4,l)=-ielem
          isfind(5,l)=n
         enddo
        enddo
       enddo
      enddo
      call quicksort(isfind,4*nelem(itet-1)+1,4*nelem(itet))
      isfcount=0
      isold(1:3)=0
      do i=4*nelem(itet-1)+1,4*nelem(itet)
       if(isfind(1,i).eq.isold(1)
     &      .and.isfind(2,i).eq.isold(2)
     &      .and.isfind(3,i).eq.isold(3))then
        isfind(4,i)=-isfind(4,i)
       else
        isfcount=isfcount+1
        isold(1:3)=isfind(1:3,i)
       endif
       isfind(6,i)=isfcount
      enddo
      nsurf(itet)=isfcount
      return
      end subroutine set_nsurf

************************************************************************
      subroutine set_tetrasurf(itet)
*
*     Setup
*       surfcoefs: coefficients of surface equation defined outward
*                  and normalized
*       ielem2surf: 4 surfaces belong to an element
*
************************************************************************
      integer,intent(in) :: itet
      integer :: i,m,n,ip(4),ielem
      real(8) :: vec1(3),vec2(3),a(4),aa,t
      do i=4*nelem(itet-1)+1,4*nelem(itet)
       ip(1:3)=isfind(1:3,i)
       ielem=isfind(4,i)
       n=isfind(5,i)
       m=isfind(6,i)+nsurf(itet-1)
       ielem2surf(n,abs(ielem))=sign(m,ielem)
       if(ielem.lt.0)then
        ip(4)=ielem2point(n,abs(ielem))
        vec1(1:3)=pointxyz(1:3,ip(2))-pointxyz(1:3,ip(1))
        vec2(1:3)=pointxyz(1:3,ip(3))-pointxyz(1:3,ip(1))
        a(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
        a(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
        a(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)
        aa=a(1)*a(1)+a(2)*a(2)+a(3)*a(3)
        a(4)=a(1)*pointxyz(1,ip(1))
     &       +a(2)*pointxyz(2,ip(1))
     &       +a(3)*pointxyz(3,ip(1))
        t=a(1)*pointxyz(1,ip(4))
     &       +a(2)*pointxyz(2,ip(4))
     &       +a(3)*pointxyz(3,ip(4))-a(4)
        if(t.gt.0)then
         aa=-1.0d0/sqrt(aa)
        else
         aa=1.0d0/sqrt(aa)
        endif
        surfcoefs(1:4,m)=aa*a(1:4)
       endif
      enddo
      return
      end subroutine set_tetrasurf

************************************************************************
      subroutine set_isurf2elem(itet)
*
*     Setup a conversion table for elements sharing a surface
*     (isurf2elem)
*
************************************************************************
      integer,intent(in) :: itet
      integer :: i,is,isold
      isurf2elem(1:2,nsurf(itet-1)+1:nsurf(itet))=0
      isold=0
      do i=4*nelem(itet-1)+1,4*nelem(itet)
       is=isfind(6,i)+nsurf(itet-1)
       if(is.ne.isold)then
        isurf2elem(1,is)=abs(isfind(4,i))
        isold=is
       else
        isurf2elem(2,is)=abs(isfind(4,i))
       endif
      enddo
      return
      end subroutine set_isurf2elem

************************************************************************
      subroutine set_outtbl(itet)
*
*     Setup a conversion table for outer surfaces of elements
*     (ioutsf2surf, ioutsf2outpt, ioutpt2point)
*
************************************************************************
      integer,intent(in) :: itet
      integer :: i,j,ip(4),isf(4),icount,jcount,np0,ns0,ipoint
      integer :: ielem,ip0,isurf
      integer,allocatable :: ipcount(:),ipoint2outpt(:)
*------------------------------------------------------------------------
      ns0=nsurf(itet)-nsurf(itet-1)
      np0=npoint(itet)-npoint(itet-1)
*------------------------------------------------------------------------
      allocate( ipcount(np0),ipoint2outpt(np0) )
*------------------------------------------------------------------------
      ioutsf2surf(nsurf(itet-1)+1:nsurf(itet))=0
      ioutsf2outpt(1:3,nsurf(itet-1)+1:nsurf(itet))=0
      ioutpt2point(npoint(itet-1)+1:npoint(itet))=0
      ipcount(1:np0)=0
      icount=nsurf(itet-1)
      do i=1,ns0
       isurf=nsurf(itet-1)+i
       if(isurf2elem(2,isurf).eq.0)then
        icount=icount+1
        ioutsf2surf(icount)=isurf
        ielem=isurf2elem(1,isurf)
        isf(1:4)=abs(ielem2surf(1:4,ielem))
        ip(1:4)=ielem2point(1:4,ielem)-npoint(itet-1)
        jcount=0
        do j=1,4
         if(isf(j).ne.isurf)then
          jcount=jcount+1
          ioutsf2outpt(jcount,icount)=ip(j)+npoint(itet-1)
          ipcount(ip(j))=ipcount(ip(j))+1
         endif
        enddo
       endif
      enddo
      noutsf(itet)=icount
      icount=npoint(itet-1)
      ipoint2outpt(1:np0)=0
      do i=1,np0
       if(ipcount(i).gt.0)then
        icount=icount+1
        ip0=npoint(itet-1)+i
        ioutpt2point(icount)=ip0
        ipoint2outpt(i)=icount
       endif
      enddo
      noutpt(itet)=icount
      do i=nsurf(itet-1)+1,noutsf(itet)
       do j=1,3
        ioutsf2outpt(j,i)
     &       =ipoint2outpt(ioutsf2outpt(j,i)-npoint(itet-1))
       enddo
      enddo
*------------------------------------------------------------------------
      deallocate( ipcount,ipoint2outpt )
*------------------------------------------------------------------------
      return
      end subroutine set_outtbl

************************************************************************
      subroutine calc_mxsf(itet,mxsf4p,nside0)
*
*     Calculate
*       mxsf4p: maximum number of outer surfaces for a point
*       nside: number of sides of outer surfaces
*
*     Last modified 2016/09/02 by T. Furuta
*
************************************************************************
      integer,intent(in) :: itet
      integer,intent(out) :: mxsf4p,nside0
      integer :: i,j,k,ipmax,ip(3),isum,np0,ns0
      integer,allocatable :: ipcount(:)
*------------------------------------------------------------------------
      np0=noutpt(itet)-npoint(itet-1)
      ns0=noutsf(itet)-nsurf(itet-1)
      allocate( ipcount(np0) )
*------------------------------------------------------------------------
      ipcount(1:np0)=0
      do i=nsurf(itet-1)+1,noutsf(itet)
       ip(1:3)=ioutsf2outpt(1:3,i)-npoint(itet-1)
       do j=1,3
        ipcount(ip(j))=ipcount(ip(j))+1
       enddo
      enddo
      ipmax=0
      do i=1,np0
       ipmax=max(ipmax,ipcount(i))
      enddo
      nside0=ns0/2*3
      mxsf4p=ipmax
*------------------------------------------------------------------------
      deallocate( ipcount )
*------------------------------------------------------------------------
      return
      end subroutine calc_mxsf

************************************************************************
      subroutine set_tmptbl(itet,mxsf4p,nside0,nside,ierr)
*
*     Setup conversion tables for surfaces sharing a point
*     (ioutsf2side, ioutpt2outsf, iside2outsf, iside2outpt)
*
*     ioutsf2side(1:3,ios) = iside
*     ioutpt2outsf(1,iop) = Nsurf: # of surfaces having common point iop
*     ioutpt2outsf(1+1:1+Nsurf,iop) : List of outer surfaces
*     isidemap(iside): Address of iside2outsf map
*     iside2outsf(ind+1:ind+n) : List of surfaces having common side
*     iside2outpt(1:2,iside) = iopoint
*
*     Last modified 2016/09/02 by T. Furuta
*
************************************************************************
      include 'err.inc'

      integer,intent(in) :: itet,mxsf4p,nside0
      integer,intent(out) :: nside
      integer,intent(inout) :: ierr
      integer :: i,j,k,i1,j1,k1,ip(3),jp(3),ip0,np0,ns0,nsp
      integer :: is,is0,is1,js,js0,js1,iside,ind,iscount,iflag
      integer,allocatable :: ipcount(:),ilist(:)
*------------------------------------------------------------------------
      np0=noutpt(itet)-npoint(itet-1)
      ns0=noutsf(itet)-nsurf(itet-1)
      allocate( ipcount(np0),ilist(mxsf4p) )
*------------------------------------------------------------------------
      ioutsf2side(1:3,1:ns0)=0
      ioutpt2outsf(1:mxsf4p+1,1:np0)=0
      iside2outsf(2*nside0)=0
      iside2outpt(1:2,1:nside0)=0
      ipcount(1:np0)=0
      iside=0
      ind=0
      do i=1,ns0
       ip(1:3)=ioutsf2outpt(1:3,i+nsurf(itet-1))-npoint(itet-1)
       do is=1,3
        ipcount(ip(is))=ipcount(ip(is))+1
        ioutpt2outsf(ipcount(ip(is))+1,ip(is))=i
       enddo
      enddo
      do i=1,np0
       ioutpt2outsf(1,i)=ipcount(i)
      enddo
      do ip0=1,np0-1
       nsp=ioutpt2outsf(1,ip0)
       ilist(1:mxsf4p)=0
       k1=0
       do i1=1,nsp
        i=ioutpt2outsf(1+i1,ip0)
        ip(1:3)=ioutsf2outpt(1:3,i+nsurf(itet-1))-npoint(itet-1)
        do is=1,3
         if(ip(is).eq.ip0)is0=is
        enddo
        do is=1,3
         if(ip(is).le.ip0)cycle
         iflag=0
         do k=1,k1
          if(ilist(k).eq.ip(is))then
            iflag=k
            exit
          endif
         enddo
         if(iflag.ne.0)cycle
         k1=k1+1
         ilist(k1)=ip(is)
         iside=iside+1
         if(iside.gt.nside0)then

           ErrCha = ''
           ErrID = 'L:1859/R:set_tmptbl/F:tetramod.f' !E80_009_001
           call ErrWrite(ErrID,ErrCha)

           write(*,'(''*** TETRA ERROR: '',
     &                      ''number of side exceed 3/2 * ns0'')')
           write(*,'(''iside ='',i9,'' ns0 ='',i9)')
     &                 iside,ns0
           ierr=1
           return
         endif
         iscount=1
         is1=6-is-is0
         isidemap(iside)=ind
         iside2outsf(ind+1)=i
         ioutsf2side(is1,i)=iside
         iside2outpt(1,iside)=ip0
         iside2outpt(2,iside)=ip(is)
         do j1=i1+1,nsp
          j=ioutpt2outsf(1+j1,ip0)
          jp(1:3)=ioutsf2outpt(1:3,j+nsurf(itet-1))-npoint(itet-1)
          do js=1,3
           if(jp(js).eq.ip0)js0=js
          enddo
          do js=1,3
           if(ip(is).eq.jp(js))then
             js1=6-js-js0
             iscount=iscount+1
             iside2outsf(ind+iscount)=j
             ioutsf2side(js1,j)=iside
           endif
          enddo
         enddo
         ind=ind+iscount
        enddo
       enddo
      enddo
      isidemap(iside+1)=ind
      nside=iside
*------------------------------------------------------------------------
      deallocate( ipcount,ilist )
*------------------------------------------------------------------------
      return
      end subroutine set_tmptbl

************************************************************************
      subroutine count_block(itet,n1,n2,n3,ierr)
*
*     Count outer sections to define length of
*     ioutsfmapEblk, indmapEblk, ioctEblk
*
*     Last modified 2017/10/20 by T. Furuta
*
************************************************************************
      include 'err.inc'

      integer,intent(in) :: itet
      integer,intent(inout) :: ierr
      integer,intent(out) :: n1,n2,n3
      integer :: ne0,numelem,numelemmax
      integer :: i,is,idiv,isec,ilocate(ndiv),ind(ndiv),indmax
      integer :: jdiv,jsec,isurf,indscf,jindex
      integer :: ishift(6,8)
      data ishift/ 0,-1, 0,-1, 0,-1,
     &     1, 0, 0,-1, 0,-1,
     &     0,-1, 1, 0, 0,-1,
     &     1, 0, 1, 0, 0,-1,
     &     0,-1, 0,-1, 1, 0,
     &     1, 0, 0,-1, 1, 0,
     &     0,-1, 1, 0, 1, 0,
     &     1, 0, 1, 0, 1, 0/
      real(8) xx(6),xmin(3),xmax(3)
      integer,allocatable :: jelemlist(:,:),ielemlist(:)
      logical :: idivflag(8,ndiv)
*------------------------------------------------------------------------
       ne0=nelem(itet)-nelem(itet-1)
*------------------------------------------------------------------------
      allocate( jelemlist(0:ne0,ndiv),ielemlist(ne0) )
*------------------------------------------------------------------------
      idivflag(1:8,1:ndiv)=.false.
      idiv=1
      jelemlist(0,1)=ne0
      do i=1,ne0
       jelemlist(i,1)=i
      enddo
      jindex=0
      indscf=0
      ind(1)=1
      indmax=9
      numelemmax=0
*------------------------------------------------------------------------
*     Box (xb(6)) is subdivided using octree algorithm until each
*     block contains only limited number of elements (ntetelem)
*------------------------------------------------------------------------
      outerloop: do i=1,numlarge
       isecloop: do isec=1,8
        if(.not.idivflag(isec,idiv))then
         indmax=indmax+1
         if(indmax.gt.mxind)then

          ErrCha = ''
          ErrID = 'L:1959/R:count_block/F:tetramod.f' !E80_010_001
          call ErrWrite(ErrID,ErrCha)

          write(*,'(''*** TETRA ERROR: '',
     &         ''octree index of element block exceeds mxind'')')
          write(*,'(''indmax ='',i9,''mxind ='',i9)')
     &         indmax,mxind
          write(*,'(''Larger ntetelem is required: '',
     &          ''Current ntetelem ='',i9)')ntetelem
          ierr=1
          return
         endif
         ilocate(idiv)=isec
         xx(1:6)=xb(1:6,itet)
         jdivloop: do jdiv=1,idiv
          jsec=ilocate(jdiv)
          isurfloop: do isurf=1,6
           is=(isurf+1)/2
           xx(isurf)=xx(isurf)+ishift(isurf,jsec)*xoctEblk(is,jdiv,itet)
          enddo isurfloop
         enddo jdivloop
         do is=1,3
          xmin(is)=xx(2*is-1)
          xmax(is)=xx(2*is)
         enddo
         call count_inblock(xmin,xmax,jelemlist(0,idiv),itet,ne0,
     &        ielemlist,numelem)
         if(numelem.le.ntetelem)then
          idivflag(isec,idiv)=.true.
          jindex=jindex+1
          indscf=indscf+numelem
*------------------------------------------------------------------------
*     No element in the block
*------------------------------------------------------------------------
          if(numelem.eq.0)then
            indscf=indscf+1
          endif
*------------------------------------------------------------------------
         else
          idiv=idiv+1
          if(idiv.gt.ndiv)then
           if(numelem.gt.numelemmax)numelemmax=numelem
           idiv=idiv-1
           idivflag(isec,idiv)=.true.
           jindex=jindex+1
           indscf=indscf+numelem
          else
            jelemlist(0,idiv)=numelem
            jelemlist(1:numelem,idiv)=ielemlist(1:numelem)
            ind(idiv)=indmax
            indmax=indmax+8
            if(indmax.gt.mxind)then

              ErrCha = ''
              ErrID = 'L:2013/R:count_block/F:tetramod.f' !E80_011_001
              call ErrWrite(ErrID,ErrCha)

              write(*,'(''*** TETRA ERROR: '',
     &          ''octree index of element block exceeds mxind'')')
              write(*,'(''indmax ='',i9,'' mxind ='',i9)')
     &             indmax,mxind
              write(*,'(''Larger ntetelem is required: '',
     &          ''Current ntetelem ='',i9)')ntetelem
              ierr=1
              return
            endif
            cycle outerloop
          endif
         endif
        endif
       enddo isecloop
       ilocate(idiv)=0
       idivflag(1:8,idiv)=.false.
       idiv=idiv-1
       if(idiv.eq.0)then
        exit outerloop
       endif
       idivflag(ilocate(idiv),idiv)=.true.
      enddo outerloop

      if(numelemmax.gt.0)then

       if(iwarning)then
        ErrCha = ''
        ErrID = 'L:2043/R:count_block/F:tetramod.f' !W80_002_001
        call ErrWrite(ErrID,ErrCha)

        write(*,'(''*** TETRA WARNING: '',
     &          ''division level of element block exceeds ndiv'')')
        write(*,'(''numelem ='',i5,'' ntetelem ='',i5)')
     &          numelemmax,ntetelem
        write(*,'(''May get faster by adjusting ntetelem '',
     &          ''larger than the above numelem'')')
       endif

      endif
      if(i.ge.numlarge)then

       ErrCha = ''
       ErrID = 'L:2058/R:count_block/F:tetramod.f' !E80_012_001
       call ErrWrite(ErrID,ErrCha)

       write(*,'(''*** TETRA ERROR: '',
     &   ''fail to set up element blocks'')')
       write(*,'(''i ='',i9,'' indmax ='',i9)')i,indmax
       write(*,'(''Larger ntetelem is required: '',
     &          ''Current ntetelem ='',i9)')ntetelem
       ierr=1
      endif
      n1=indscf
      n2=jindex+1
      n3=indmax
*------------------------------------------------------------------------
      deallocate( jelemlist,ielemlist )
*------------------------------------------------------------------------
      return
      end subroutine count_block

************************************************************************
      subroutine count_inblock(xmin,xmax,jelemlist,itet,ne0,
     &     ielemlist,numelem)
*
*     Find elements for which points or center of mass are in the block
*
*     ielemlist: lists indicating elements inside the block
*     jelemlist: remember status of ielemlist in upper level
*     numelem: number of elements in the block
*
*     Last modified 2017/10/20 by T. Furuta
*
************************************************************************
      real(8),intent(in) :: xmin(3),xmax(3)
      integer,intent(in) :: jelemlist(0:ne0)
      integer,intent(in) :: itet,ne0
      integer,intent(out) :: ielemlist(ne0)
      integer,intent(out) :: numelem
      integer :: i,j,isum,ip(4)
      logical :: iflag(3),jflag
      real(8) :: xsum(3)
*------------------------------------------------------------------------
*     Find elements for which points or center of mass are in the block
*------------------------------------------------------------------------
      isum=0
      do i=1,jelemlist(0)
       ip(1:4)=ielem2point(1:4,jelemlist(i)+nelem(itet-1))
       xsum(1:3)=0.0d0
       jflag=.false.
       do j=1,4
        iflag(1:3)=pointxyz(1:3,ip(j)).ge.xmin(1:3).and.
     &       pointxyz(1:3,ip(j)).le.xmax(1:3)
        if(all(iflag))then
         jflag=.true.
         isum=isum+1
         ielemlist(isum)=jelemlist(i)
         exit
        endif
        xsum(1:3)=xsum(1:3)+pointxyz(1:3,ip(j))
       enddo
       if(.not.jflag)then
        xsum(1:3)=0.25d0*xsum(1:3)
        iflag(1:3)=xsum(1:3).ge.xmin(1:3).and.
     &       xsum(1:3).le.xmax(1:3)
        if(all(iflag))then
         isum=isum+1
         ielemlist(isum)=jelemlist(i)
        endif
       endif
      enddo
      numelem=isum
      return
      end subroutine count_inblock

************************************************************************
      subroutine set_block(itet,ierr)
*
*     Setup outer volume of universe filling tetrahedron mesh
*
*     Last modified 2017/10/20 by T. Furuta
*
************************************************************************
      integer,intent(in) :: itet
      integer,intent(inout) :: ierr
      integer :: ne0,nex,numelem
      integer :: i,is,idiv,isec,ilocate(ndiv),ind(ndiv),indmax
      integer :: jdiv,jsec,isurf,indscf,jindex
      integer :: ishift(6,8)
      data ishift/ 0,-1, 0,-1, 0,-1,
     &     1, 0, 0,-1, 0,-1,
     &     0,-1, 1, 0, 0,-1,
     &     1, 0, 1, 0, 0,-1,
     &     0,-1, 0,-1, 1, 0,
     &     1, 0, 0,-1, 1, 0,
     &     0,-1, 1, 0, 1, 0,
     &     1, 0, 1, 0, 1, 0/
      real(8) xx(6),xmin(3),xmax(3)
      integer,allocatable :: jelemlist(:,:),ielemlist(:)
      logical :: idivflag(8,ndiv)
*------------------------------------------------------------------------
      ne0=nelem(itet)-nelem(itet-1)
      allocate( jelemlist(0:ne0,ndiv),ielemlist(ne0) )
*------------------------------------------------------------------------
      idivflag(1:8,1:ndiv)=.false.
      idiv=1
      jelemlist(0,1)=ne0
      do i=1,ne0
       jelemlist(i,1)=i
      enddo
      jindex=nindEblk(itet-1)
      indscf=nEblk(itet-1)
      ind(1)=noctEblk(itet-1)+1
      ioctEblk(ind(1))=0
      indmax=ind(1)+8
*------------------------------------------------------------------------
*     Box (xb(6)) is subdivided using octree algorithm until each
*     block contains only limited number of elements (ntetelem)
*------------------------------------------------------------------------
      outerloop: do i=1,numlarge
       isecloop: do isec=1,8
        if(.not.idivflag(isec,idiv))then
         indmax=indmax+1
         ioctEblk(ind(idiv)+isec)=indmax
         ilocate(idiv)=isec
         xx(1:6)=xb(1:6,itet)
         jdivloop: do jdiv=1,idiv
          jsec=ilocate(jdiv)
          isurfloop: do isurf=1,6
           is=(isurf+1)/2
           xx(isurf)=xx(isurf)+ishift(isurf,jsec)*xoctEblk(is,jdiv,itet)
          enddo isurfloop
         enddo jdivloop
         do is=1,3
          xmin(is)=xx(2*is-1)
          xmax(is)=xx(2*is)
         enddo
         call count_inblock(xmin,xmax,jelemlist(0,idiv),itet,ne0,
     &        ielemlist,numelem)
         if(numelem.le.ntetelem)then
          idivflag(isec,idiv)=.true.
          jindex=jindex+1
          ioctEblk(indmax)=-jindex
          call set_ielemblkmap(itet,jindex,idiv,ilocate,numelem,
     &         ne0,ielemlist,indscf)
*------------------------------------------------------------------------
*     No element in the block
*------------------------------------------------------------------------
          if(numelem.eq.0)then
            ielemmapEblk(indscf+1)=-1
            indscf=indscf+1
          endif
*------------------------------------------------------------------------
         else
          idiv=idiv+1
          if(idiv.gt.ndiv)then
            idiv=idiv-1
            idivflag(isec,idiv)=.true.
            jindex=jindex+1
            ioctEblk(indmax)=-jindex
            call set_ielemblkmap(itet,jindex,idiv,ilocate,numelem,
     &           ne0,ielemlist,indscf)
          else
            jelemlist(0,idiv)=numelem
            jelemlist(1:numelem,idiv)=ielemlist(1:numelem)
            ioctEblk(indmax)=0
            ind(idiv)=indmax
            indmax=indmax+8
            cycle outerloop
          endif
         endif
        endif
       enddo isecloop
       ilocate(idiv)=0
       idivflag(1:8,idiv)=.false.
       idiv=idiv-1
       if(idiv.eq.0)then
        exit outerloop
       endif
       idivflag(ilocate(idiv),idiv)=.true.
      enddo outerloop
      indmapEblk(jindex+1,1)=0
      indmapEblk(jindex+1,2)=indscf
*------------------------------------------------------------------------
      deallocate( jelemlist,ielemlist )
*------------------------------------------------------------------------
      return
      end subroutine set_block

************************************************************************
      subroutine set_ielemblkmap(itet,jindex,idiv,ilocate,numelem,ne0,
     &     ielemlist,indscf)
*
*     Create map for outer sections
*
*     ioctEblk(ind) = 0 : Contains deeper level
*                     = -jindex : outer section ID
*     ioctEblk(ind+1:ind+8) = ind
*     indmapEblk(jindex,1) = Octree index (ex. 82 = section 2-8)
*     indmapEblk(jindex,2) = indscf: Map index for elements
*     ielemmapEblk(indscf+1:indscf+numsurf)
*          : element indexes inside the block
*
*     Last modified 2017/10/20 by T. Furuta
*
************************************************************************
      integer,intent(in) :: itet,jindex,idiv,ilocate(ndiv),numelem,ne0
      integer,intent(in) :: ielemlist(1:ne0)
      integer,intent(inout) :: indscf
      integer :: i,jdiv,jsec
      jsec=0
      do jdiv=1,idiv
       jsec=jsec+ilocate(jdiv)*10**(jdiv-1)
      enddo
      indmapEblk(jindex,1)=jsec
      indmapEblk(jindex,2)=indscf
      do i=1,numelem
       ielemmapEblk(indscf+i)=ielemlist(i)+nelem(itet-1)
      enddo
      indscf=indscf+numelem
      return
      end subroutine set_ielemblkmap

************************************************************************
      subroutine count_outer(xx0,itet,nside,n1,n2,n3,ierr)
*
*     Count outer sections to define length of
*     ioutsfmapOsec, indmapOsec, ioctOsec
*
*     Last modified 2017/10/20 by T. Furuta
*
************************************************************************
      include 'err.inc'

      real(8),intent(in) :: xx0(6)
      integer,intent(in) :: itet,nside
      integer,intent(inout) :: ierr
      integer,intent(out) :: n1,n2,n3
      integer :: np0,ns0,numsurf,numpoint,numside,numsurfmax
      integer :: i,is,idiv,isec,ilocate(ndiv),ind(ndiv),indmax
      integer :: jdiv,jsec,isurf,indscf,jindex
      integer :: ishift(6,8)
      data ishift/ 0,-1, 0,-1, 0,-1,
     &     1, 0, 0,-1, 0,-1,
     &     0,-1, 1, 0, 0,-1,
     &     1, 0, 1, 0, 0,-1,
     &     0,-1, 0,-1, 1, 0,
     &     1, 0, 0,-1, 1, 0,
     &     0,-1, 1, 0, 1, 0,
     &     1, 0, 1, 0, 1, 0/
      real(8) xx(6),x(3),u(3)
      integer :: ielem,ires,lgc(4)
      integer,allocatable :: jsuflist(:,:),jplist(:,:),jsidelist(:,:)
      integer,allocatable :: isuflist(:),iplist(:),isidelist(:)
      logical,allocatable :: isufflag(:),isideflag(:)
      logical :: idivflag(8,ndiv)

*------------------------------------------------------------------------
       np0=noutpt(itet)-npoint(itet-1)
       ns0=noutsf(itet)-nsurf(itet-1)
*------------------------------------------------------------------------
      allocate( jsuflist(0:ns0,ndiv),isuflist(ns0) )
      allocate( jplist(0:np0,ndiv),iplist(np0) )
      allocate( jsidelist(0:nside,ndiv),isidelist(nside) )
      allocate( isufflag(ns0),isideflag(nside) )
*------------------------------------------------------------------------
      idivflag(1:8,1:ndiv)=.false.
      idiv=1
      jsuflist(0,1)=ns0
      do i=1,ns0
       jsuflist(i,1)=i
      enddo
      jplist(0,1)=np0
      do i=1,np0
       jplist(i,1)=i
      enddo
      jsidelist(0,1)=nside
      do i=1,nside
       jsidelist(i,1)=i
      enddo
      jindex=0
      indscf=0
      ind(1)=1
      indmax=9
      numsurfmax=0
*------------------------------------------------------------------------
*     Outer volume is subdivided using octree algorithm until each
*     section contains only limited number of triangle surfaces of outer
*     tetrahedrons (ntetsurf)
*------------------------------------------------------------------------
      outerloop: do i=1,numlarge
       isecloop: do isec=1,8
        if(.not.idivflag(isec,idiv))then
         indmax=indmax+1
         if(indmax.gt.mxind)then

          ErrCha = ''
          ErrID = 'L:2353/R:count_outer/F:tetramod.f' !E80_013_001
          call ErrWrite(ErrID,ErrCha)

          write(*,'(''*** TETRA ERROR: '',
     &         ''octree index of outer section exceeds mxind'')')
          write(*,'(''indmax ='',i9,''mxind ='',i9)')
     &         indmax,mxind
          write(*,'(''Larger ntetsurf is required: '',
     &          ''Current ntetsurf ='',i9)')ntetsurf
          ierr=1
          return
         endif
         ilocate(idiv)=isec
         xx(1:6)=xx0(1:6)
         jdivloop: do jdiv=1,idiv
          jsec=ilocate(jdiv)
          isurfloop: do isurf=1,6
           is=(isurf+1)/2
           xx(isurf)=xx(isurf)+ishift(isurf,jsec)*xoctOsec(is,jdiv,itet)
          enddo isurfloop
         enddo jdivloop
         call set_iflags(xx,
     &        jsuflist(0,idiv),jplist(0,idiv),jsidelist(0,idiv),
     &        itet,nside,np0,ns0,isuflist,iplist,isidelist,
     &        isufflag,isideflag,
     &        numsurf,numpoint,numside)
         if(numsurf.le.ntetsurf)then
          idivflag(isec,idiv)=.true.
          jindex=jindex+1
          indscf=indscf+numsurf
*------------------------------------------------------------------------
*     No outer surface cross this section
*     => either outside of tetrahedrons or inside of a tetrahedron
*------------------------------------------------------------------------
          if(numsurf.eq.0)then
           x(1)=0.5d0*(xx(1)+xx(2))
           x(2)=0.5d0*(xx(3)+xx(4))
           x(3)=0.5d0*(xx(5)+xx(6))
           xx(1:6)=xx0(1:6)
           do jdiv=1,idiv-1 ! upto one-dimension above
            jsec=ilocate(jdiv)
            do isurf=1,6
             is=(isurf+1)/2
             xx(isurf)=xx(isurf)
     &            +ishift(isurf,jsec)*xoctOsec(is,jdiv,itet)
            enddo
           enddo
           call check_center(xx,x,itet,ns0,jsuflist(0,idiv),ires) !FURUTA20221104
           if(ires.lt.0)then
            u(1:3)=1.0d0
            do ielem=nelem(itet-1)+1,nelem(itet)
             call check_inelem(x,u,0.0d0,ielem,0,ires,lgc)
             if(ires.eq.0)exit
            enddo
           endif
           if(ires.eq.0)then    ! inside of a tetrahedron
            indscf=indscf+1
           endif
          endif
*------------------------------------------------------------------------
         else
          idiv=idiv+1
          if(idiv.gt.ndiv)then
            if(numsurf.gt.numsurfmax)numsurfmax=numsurf
            idiv=idiv-1
            idivflag(isec,idiv)=.true.
            jindex=jindex+1
            indscf=indscf+numsurf
          else
           jsuflist(0,idiv)=numsurf
           jsuflist(1:numsurf,idiv)=isuflist(1:numsurf)
           jplist(0,idiv)=numpoint
           jplist(1:numpoint,idiv)=iplist(1:numpoint)
           jsidelist(0,idiv)=numside
           jsidelist(1:numside,idiv)=isidelist(1:numside)
            ind(idiv)=indmax
            indmax=indmax+8
            if(indmax.gt.mxind)then

             ErrCha = ''
             ErrID = 'L:2433/R:count_outer/F:tetramod.f' !E80_013_002
             call ErrWrite(ErrID,ErrCha)

              write(*,'(''*** TETRA ERROR: '',
     &          ''octree index of outer section exceeds mxind'')')
              write(*,'(''indmax ='',i9,'' mxind ='',i9)')
     &             indmax,mxind
              write(*,'(''Larger ntetsurf is required: '',
     &          ''Current ntetsurf ='',i9)')ntetsurf
              ierr=1
              return
            endif
            cycle outerloop
          endif
         endif
        endif
       enddo isecloop
       ilocate(idiv)=0
       idivflag(1:8,idiv)=.false.
       idiv=idiv-1
       if(idiv.eq.0)then
        exit outerloop
       endif
       idivflag(ilocate(idiv),idiv)=.true.
      enddo outerloop
      if(numsurfmax.gt.0)then

       if(iwarning)then
        ErrCha = ''
        ErrID = 'L:2462/R:count_outer/F:tetramod.f' !W80_002_002
        call ErrWrite(ErrID,ErrCha)

        write(*,'(''*** TETRA WARNING: '',
     &          ''division level of outer section exceeds ndiv'')')
        write(*,'(''numsurf ='',i5,'' ntetsurf ='',i5)')
     &           numsurfmax,ntetsurf
        write(*,'(''May get faster by adjusting ntetsurf '',
     &          ''larger than the above numsurf'')')
       endif

      endif

      if(i.ge.numlarge)then

       ErrCha = ''
       ErrID = 'L:2478/R:count_outer/F:tetramod.f' !E80_014_001
       call ErrWrite(ErrID,ErrCha)

       write(*,'(''*** TETRA ERROR: '',
     &''fail to set up outer sections'')')
       write(*,'(''i ='',i9,'' indmax ='',i9)')i,indmax
       write(*,'(''Larger ntetsurf is required: '',
     &          ''Current ntetsurf ='',i9)')ntetsurf
       ierr=1
      endif
      n1=indscf
      n2=jindex+1
      n3=indmax
*------------------------------------------------------------------------
      deallocate( jsuflist,jplist,jsidelist,isuflist,iplist,isidelist )
      deallocate( isufflag,isideflag )
*------------------------------------------------------------------------
      return
      end subroutine count_outer

************************************************************************
      subroutine check_center(xx,x,itet,ns0,jsuflist,ires)
*
*     Check center position x is inside of tetrahedrons
*
*     Last modified 2017/10/20
*
************************************************************************
      include 'err.inc'

      real(8),intent(in) :: xx(6),x(3)
      integer,intent(in) :: itet,ns0,jsuflist(0:ns0)
      integer,intent(out) :: ires
      integer :: i,j,is,ies,jes,nes,icand,ip(3)
      integer :: icheck(2)
      integer,allocatable :: isflag(:)
      real(8) :: scf(4),t1,t4,t5,d1,dmin
      real(8) :: v(3),w(3),xp(3,3),xc(3),vvec(3,3),xcm(3),v2cm(3),xcc(3)
      logical :: iflag(3)
      integer :: nflag

*------------------------------------------------------------------------
*     Outside of element block
*------------------------------------------------------------------------
      if( x(1).lt.xb(1,itet).or.x(1).gt.xb(2,itet)
     &     .or.x(2).lt.xb(3,itet).or.x(2).gt.xb(4,itet)
     &     .or.x(3).lt.xb(5,itet).or.x(3).gt.xb(6,itet) )then
       ires=1
       return
      endif
*------------------------------------------------------------------------
      nes=jsuflist(0)
      allocate( isflag(nes) )
      icheck(1:2)=0
      do ies=1,nes
       j=jsuflist(ies)
       scf(1:4)=surfcoefs(1:4,ioutsf2surf(j+nsurf(itet-1)))
       t4=scf(1)*x(1)+scf(2)*x(2)+scf(3)*x(3)-scf(4)
       if(t4.gt.xtiny*xqt(itet))then ! outside
        isflag(ies)=-1
        icheck(1)=icheck(1)+1
       elseif(t4.lt.-xtiny*xqt(itet))then ! inside
        isflag(ies)=1
        icheck(2)=icheck(2)+1
       else
        isflag(ies)=0
       endif
      enddo
      if(icheck(2).eq.nes)then  ! all surfaces are inside
       ires=0
      elseif(icheck(1).eq.nes)then ! all surfaces are outside
       ires=1
      else
*------------------------------------------------------------------------
*     Check face of nearest surface found in direction from x to CM of
*     a surface facing outside
*------------------------------------------------------------------------
       icand=0
       do ies=1,nes+8
        j=jsuflist(ies)
        if(ies.gt.nes)then  ! Vertexes of box
         xcm(1)=xx(mod(ies-nes+1,2)+1)
         xcm(2)=xx(mod((ies-nes+1)/2+1,2)+3)
         xcm(3)=xx((ies-nes-1)/4+5)
         v2cm(1:3)=xcm(1:3)-x(1:3)
         dmin=sqrt(v2cm(1)*v2cm(1)+v2cm(2)*v2cm(2)+v2cm(3)*v2cm(3))
         v2cm(1:3)=v2cm(1:3)/dmin
        else ! CM of surfaces
         ip(1:3)=ioutsf2outpt(1:3,j+nsurf(itet-1))
         xp(1:3,1)=pointxyz(1:3,ioutpt2point(ip(1)))
         xp(1:3,2)=pointxyz(1:3,ioutpt2point(ip(2)))
         xp(1:3,3)=pointxyz(1:3,ioutpt2point(ip(3)))
         xcm(1:3)=(xp(1:3,1)+xp(1:3,2)+xp(1:3,3))/3.0d0
         xc(1:3)=xcm(1:3)
         do is=1,4
          if(xcm(1).ge.xx(1).and.xcm(1).le.xx(2)
     &         .and.xcm(2).ge.xx(3).and.xcm(2).le.xx(4)
     &         .and.xcm(3).ge.xx(5).and.xcm(3).le.xx(6))then
           v2cm(1:3)=xcm(1:3)-x(1:3)
           dmin=sqrt(v2cm(1)*v2cm(1)+v2cm(2)*v2cm(2)+v2cm(3)*v2cm(3))
           v2cm(1:3)=v2cm(1:3)/dmin
           icand=ies
           exit
          endif
          if(is.eq.4)exit
          xcm(1:3)=xp(1:3,is)*0.98d0+xc(1:3)*0.02d0 !FURUTA20180517
         enddo
         if(icand.eq.0)cycle
        endif
        do jes=1,nes
         if(jes.eq.ies)cycle
         j=jsuflist(jes)
         scf(1:4)=surfcoefs(1:4,ioutsf2surf(j+nsurf(itet-1)))
         t1=v2cm(1)*scf(1)+v2cm(2)*scf(2)+v2cm(3)*scf(3)
         if(abs(t1).le.xtiny)cycle
         d1=(scf(4)-x(1)*scf(1)-x(2)*scf(2)-x(3)*scf(3))/t1
         if(d1.lt.0)cycle
         if(d1.lt.dmin+xtiny*xqt(itet))then !FURUTA20180517
          xc(1:3)=v2cm(1:3)*d1+x(1:3)
          ip(1:3)=ioutsf2outpt(1:3,j+nsurf(itet-1))
          xp(1:3,1)=pointxyz(1:3,ioutpt2point(ip(1)))
          xp(1:3,2)=pointxyz(1:3,ioutpt2point(ip(2)))
          xp(1:3,3)=pointxyz(1:3,ioutpt2point(ip(3)))
          xcm(1:3)=(xp(1:3,1)+xp(1:3,2)+xp(1:3,3))/3.0d0 !FURUTA20180517
          xcc(1:3)=1.01d0*xc(1:3)-0.01d0*xcm(1:3)        !FURUTA20180517
          vvec(1:3,1)=xp(1:3,2)-xp(1:3,1)
          vvec(1:3,2)=xp(1:3,3)-xp(1:3,2)
          vvec(1:3,3)=xp(1:3,1)-xp(1:3,3)
          do i=1,3
           w(1:3)=xcc(1:3)-xp(1:3,i) !FURUTA20180517
           v(1)=vvec(2,i)*w(3)-vvec(3,i)*w(2)
           v(2)=vvec(3,i)*w(1)-vvec(1,i)*w(3)
           v(3)=vvec(1,i)*w(2)-vvec(2,i)*w(1)
           iflag(i)=scf(1)*v(1)+scf(2)*v(2)+scf(3)*v(3).ge.0.0d0
          enddo
          nflag=count(iflag)
          if(nflag.eq.3.or.nflag.eq.0)then ! All flags are equivalent
cFURUTA20180517------------------------------------------------------
           if(d1.gt.dmin-xtiny*xqt(itet)
     &          .and.d1.lt.dmin+xtiny*xqt(itet))then
            icand=0
            exit ! Exception for non-sharing surfaces
           elseif(d1.lt.dmin)then
            icand=jes
            dmin=d1
           endif
          else
           xcc(1:3)=0.99d0*xc(1:3)+0.01d0*xcm(1:3)
           do i=1,3
            w(1:3)=xcc(1:3)-xp(1:3,i)
            v(1)=vvec(2,i)*w(3)-vvec(3,i)*w(2)
            v(2)=vvec(3,i)*w(1)-vvec(1,i)*w(3)
            v(3)=vvec(1,i)*w(2)-vvec(2,i)*w(1)
            iflag(i)=scf(1)*v(1)+scf(2)*v(2)+scf(3)*v(3).ge.0.0d0
           enddo
           nflag=count(iflag)
           if(nflag.eq.3.or.nflag.eq.0)then
            icand=0 ! Exception for surface edge
            exit
           endif
          endif
         endif
        enddo
        if(icand.gt.0)then
         if(isflag(icand).eq.0)then
          icand=0
         else
          if(isflag(icand).lt.0)then
           ires=1
          else
           ires=0
          endif
          exit
         endif
c--------------------------------------------------------------------
        endif
       enddo
       if(icand.eq.0)then
        if(iwarning)then
         write(ErrCha,'(''*** TETRA WARNING: '',
     &       ''fail to check inside tetrahedrons for center'')')
          ErrID = 'L:2659/R:check_center/F:tetramod.f' !W80_003_001
          call ErrWrite(ErrID,ErrCha)
        endif
        ires=-1 ! Check inside tetrahedrons
       endif
      endif
      return
      end subroutine check_center

************************************************************************
      subroutine set_iflags(xx,jsuflist,jplist,jsidelist,
     &     itet,nside,np0,ns0,
     &     isuflist,iplist,isidelist,
     &     isufflag,isideflag,
     &     numsurf,numpoint,numside)
*
*     Find triangle surfaces intersects the section specified by xx
*
*     ilist (isuflist, iplist, isidelist):
*          lists indicating (surface, point, side) inside the section
*     jlist (jsuflist, jplist, jsidelist):
*          remember status of ilists in upper level
*     numsurf:
*          number of triangle surfaces intersecting the section
*     numpoint:
*          number of points inside the section
*     numside:
*          number of sides inside or intersecting the section
*
*     Last modified 2017/10/20 by T. Furuta
*
************************************************************************
      real(8),intent(in) :: xx(6)
      integer,intent(in) :: jsuflist(0:ns0)
      integer,intent(in) :: jplist(0:np0),jsidelist(0:nside)
      integer,intent(in) :: itet,nside,np0,ns0
      integer,intent(out) :: isuflist(ns0)
      integer,intent(out) :: iplist(np0),isidelist(nside)
      logical,intent(out) :: isufflag(ns0),isideflag(nside)
      integer,intent(out) :: numsurf,numpoint,numside
      integer :: i,j,k,is1,is2,is3,isum,ip(3),iside
      integer :: isurf,isurf1,isurf2,isurf3,jsurf
      real(8) :: x(3),xp(3,3),v(3),u(3),vec(3),vvec(3,3),a(4),t
      real(8) :: xx1,xx2
      logical :: iflag(3),jflag
      integer :: nflag
*------------------------------------------------------------------------
      do i=1,jsuflist(0)
       isufflag(jsuflist(i))=.false.
      enddo
      do i=1,jsidelist(0)
       isideflag(jsidelist(i))=.false.
      enddo
*------------------------------------------------------------------------
*     Find points of triangles inside the section
*
*     A triangle intersects the section if any point of the triangle is
*     inside the section
*------------------------------------------------------------------------
      isum=0
      do i=1,jplist(0)
       j=jplist(i)
       x(1:3)=pointxyz(1:3,ioutpt2point(j+npoint(itet-1)))
       if((x(1).ge.xx(1)).and.(x(1).lt.xx(2))
     &      .and.(x(2).ge.xx(3)).and.(x(2).lt.xx(4))
     &      .and.(x(3).ge.xx(5)).and.(x(3).lt.xx(6)))then
        isum=isum+1
        iplist(isum)=j
        jsurf=ioutpt2outsf(1,j)
        do isurf=1,jsurf
         isurf1=ioutpt2outsf(1+isurf,j)
         isufflag(isurf1)=.true.
         do is1=1,3
          if(ioutsf2outpt(is1,isurf1).ne.j+npoint(itet-1))then
           iside=ioutsf2side(is1,isurf1)
           isideflag(iside)=.true.
          endif
         enddo
        enddo
       endif
      enddo
      numpoint=isum
*------------------------------------------------------------------------
*     Find sides of triangles intesect any surface of the section
*
*     A triangle intersects the section if any side of the triangle
*     intesects the 6 surfaces of the section
*------------------------------------------------------------------------
      isum=0
      do i=1,jsidelist(0)
       j=jsidelist(i)
       if(isideflag(j))then
        isum=isum+1
        isidelist(isum)=j
       else
        ip(1:2)=iside2outpt(1:2,j)+npoint(itet-1)
        xp(1:3,1)=pointxyz(1:3,ioutpt2point(ip(1)))
        xp(1:3,2)=pointxyz(1:3,ioutpt2point(ip(2)))
        vec(1:3)=xp(1:3,2)-xp(1:3,1)
        do isurf=1,6
         is1=(isurf+1)/2
         if(abs(vec(is1)).gt.xtiny)then
          t=(xx(isurf)-xp(is1,1))/vec(is1)
          if(t.ge.0.0d0.and.t.le.1.0d0)then
           is2=mod(is1,3)+1
           is3=mod(is1+1,3)+1
           isurf2=2*is2-1
           isurf3=2*is3-1
           xx1=vec(is2)*t+xp(is2,1)
           xx2=vec(is3)*t+xp(is3,1)
           if((xx1.gt.xx(isurf2)-abs(xx(isurf2))*xtiny)
     &          .and.(xx1.lt.xx(isurf2+1)+abs(xx(isurf2+1))*xtiny)
     &          .and.(xx2.gt.xx(isurf3)-abs(xx(isurf3))*xtiny)
     &          .and.(xx2.lt.xx(isurf3+1)+abs(xx(isurf3+1))*xtiny))then
            isideflag(j)=.true.
            exit
           endif
          endif
         endif
        enddo
        if(isideflag(j))then
         isum=isum+1
         isidelist(isum)=j
         do is1=isidemap(j)+1,isidemap(j+1)
          isufflag(iside2outsf(is1))=.true.
         enddo
        endif
       endif
      enddo
      numside=isum
*------------------------------------------------------------------------
*     Find triangles intersect any axis of the section
*
*     A triangle intersects the section if any side of the section
*     is intersected by the triangle
*------------------------------------------------------------------------
      isum=0
      do i=1,jsuflist(0)
       j=jsuflist(i)
       if(isufflag(j))then
        isum=isum+1
        isuflist(isum)=j
       else
        ip(1:3)=ioutsf2outpt(1:3,j+nsurf(itet-1))
        xp(1:3,1)=pointxyz(1:3,ioutpt2point(ip(1)))
        xp(1:3,2)=pointxyz(1:3,ioutpt2point(ip(2)))
        xp(1:3,3)=pointxyz(1:3,ioutpt2point(ip(3)))
        vvec(1:3,1)=xp(1:3,2)-xp(1:3,1)
        vvec(1:3,2)=xp(1:3,3)-xp(1:3,2)
        vvec(1:3,3)=xp(1:3,1)-xp(1:3,3)
        a(1:4)=surfcoefs(1:4,ioutsf2surf(j+nsurf(itet-1)))
        is1loop: do is1=1,3
         isurf1loop: do isurf1=1,2
          x(is1)=xx((is1-1)*2+isurf1)
          is2loop: do is2=is1+1,3
           isurf2loop: do isurf2=1,2
            x(is2)=xx((is2-1)*2+isurf2)
            is3=6-is1-is2
            if(abs(a(is3)).gt.xtiny)then
             x(is3)=(a(4)-a(is1)*x(is1)-a(is2)*x(is2))/a(is3)
            endif
            if(x(is3).gt.xx((is3-1)*2+1)
     &           .and.x(is3).lt.xx((is3-1)*2+2))then
             kloop: do k=1,3
              u(1:3)=x(1:3)-xp(1:3,k)
              v(1)=vvec(2,k)*u(3)-vvec(3,k)*u(2)
              v(2)=vvec(3,k)*u(1)-vvec(1,k)*u(3)
              v(3)=vvec(1,k)*u(2)-vvec(2,k)*u(1)
              iflag(k)=a(1)*v(1)+a(2)*v(2)+a(3)*v(3).ge.0.0d0
             enddo kloop
             nflag=count(iflag)
             if(nflag.eq.3.or.nflag.eq.0)then
              isum=isum+1
              isuflist(isum)=j
              exit is1loop
             endif
            endif
           enddo isurf2loop
          enddo is2loop
         enddo isurf1loop
        enddo is1loop
       endif
      enddo
      numsurf=isum
      return
      end subroutine set_iflags

************************************************************************
      subroutine set_outer(xx0,itet,nside,ierr)
*
*     Setup outer volume of universe filling tetrahedron mesh
*
*     Last modified 2017/10/20 by T. Furuta
*
************************************************************************
      real(8),intent(in) :: xx0(6)
      integer,intent(in) :: itet,nside
      integer,intent(inout) :: ierr
      integer :: np0,ns0,numsurf,numpoint,numside
      integer :: i,is,idiv,isec,ilocate(ndiv),ind(ndiv),indmax
      integer :: jdiv,jsec,isurf,indscf,jindex
      integer :: ishift(6,8)
      data ishift/ 0,-1, 0,-1, 0,-1,
     &     1, 0, 0,-1, 0,-1,
     &     0,-1, 1, 0, 0,-1,
     &     1, 0, 1, 0, 0,-1,
     &     0,-1, 0,-1, 1, 0,
     &     1, 0, 0,-1, 1, 0,
     &     0,-1, 1, 0, 1, 0,
     &     1, 0, 1, 0, 1, 0/
      real(8) xx(6),x(3),u(3)
      integer :: ielem,ires,lgc(4)
      integer,allocatable :: jsuflist(:,:),jplist(:,:),jsidelist(:,:)
      integer,allocatable :: isuflist(:),iplist(:),isidelist(:)
      logical,allocatable :: isufflag(:),isideflag(:)
      logical :: idivflag(8,ndiv)
*------------------------------------------------------------------------
      np0=noutpt(itet)-npoint(itet-1)
      ns0=noutsf(itet)-nsurf(itet-1)
      allocate( jsuflist(0:ns0,ndiv),isuflist(ns0) )
      allocate( jplist(0:np0,ndiv),iplist(np0) )
      allocate( jsidelist(0:nside,ndiv),isidelist(nside) )
      allocate( isufflag(ns0),isideflag(nside) )
*------------------------------------------------------------------------
      idivflag(1:8,1:ndiv)=.false.
      idiv=1
      jsuflist(0,1)=ns0
      do i=1,ns0
       jsuflist(i,1)=i
      enddo
      jplist(0,1)=np0
      do i=1,np0
       jplist(i,1)=i
      enddo
      jsidelist(0,1)=nside
      do i=1,nside
       jsidelist(i,1)=i
      enddo
      jindex=nindOsec(itet-1)
      indscf=nOsec(itet-1)
      ind(1)=noctOsec(itet-1)+1
      ioctOsec(ind(1))=0
      indmax=ind(1)+8
*------------------------------------------------------------------------
*     Outer volume is subdivided using octree algorithm until each
*     section contains only limited number of triangle surfaces of outer
*     tetrahedrons (ntetsurf)
*------------------------------------------------------------------------
      outerloop: do i=1,numlarge
       isecloop: do isec=1,8
        if(.not.idivflag(isec,idiv))then
         indmax=indmax+1
         ioctOsec(ind(idiv)+isec)=indmax
         ilocate(idiv)=isec
         xx(1:6)=xx0(1:6)
         jdivloop: do jdiv=1,idiv
          jsec=ilocate(jdiv)
          isurfloop: do isurf=1,6
           is=(isurf+1)/2
           xx(isurf)=xx(isurf)+ishift(isurf,jsec)*xoctOsec(is,jdiv,itet)
          enddo isurfloop
         enddo jdivloop
         call set_iflags(xx,
     &        jsuflist(0,idiv),jplist(0,idiv),jsidelist(0,idiv),
     &        itet,nside,np0,ns0,isuflist,iplist,isidelist,
     &        isufflag,isideflag,
     &        numsurf,numpoint,numside)
         if(numsurf.le.ntetsurf)then
          idivflag(isec,idiv)=.true.
          jindex=jindex+1
          ioctOsec(indmax)=-jindex
          call set_ioutersecmap(itet,jindex,idiv,ilocate,numsurf,
     &         ns0,isuflist,indscf)
*------------------------------------------------------------------------
*     No outer surface cross this section
*     => either outside of tetrahedrons or inside of a tetrahedron
*------------------------------------------------------------------------
          if(numsurf.eq.0)then
           x(1)=0.5d0*(xx(1)+xx(2))
           x(2)=0.5d0*(xx(3)+xx(4))
           x(3)=0.5d0*(xx(5)+xx(6))
           xx(1:6)=xx0(1:6)
           do jdiv=1,idiv-1 ! upto one-dimension above
            jsec=ilocate(jdiv)
            do isurf=1,6
             is=(isurf+1)/2
             xx(isurf)=xx(isurf)
     &            +ishift(isurf,jsec)*xoctOsec(is,jdiv,itet)
            enddo
           enddo
           call check_center(xx,x,itet,ns0,jsuflist(0,idiv),ires) !FURUTA20221121
           if(ires.lt.0)then
            u(1:3)=1.0d0
            do ielem=nelem(itet-1)+1,nelem(itet)
             call check_inelem(x,u,0.0d0,ielem,0,ires,lgc)
             if(ires.eq.0)exit
            enddo
           endif
           if(ires.eq.0)then ! inside of a tetrahedron
            ioutsfmapOsec(indscf+1)=-1
            indscf=indscf+1
           endif
          endif
*------------------------------------------------------------------------
         else
          idiv=idiv+1
          if(idiv.gt.ndiv)then
            idiv=idiv-1
            idivflag(isec,idiv)=.true.
            jindex=jindex+1
            ioctOsec(indmax)=-jindex
            call set_ioutersecmap(itet,jindex,idiv,ilocate,numsurf,
     &           ns0,isuflist,indscf)
          else
           jsuflist(0,idiv)=numsurf
           jsuflist(1:numsurf,idiv)=isuflist(1:numsurf)
           jplist(0,idiv)=numpoint
           jplist(1:numpoint,idiv)=iplist(1:numpoint)
           jsidelist(0,idiv)=numside
           jsidelist(1:numside,idiv)=isidelist(1:numside)
            ioctOsec(indmax)=0
            ind(idiv)=indmax
            indmax=indmax+8
            cycle outerloop
          endif
         endif
        endif
       enddo isecloop
       ilocate(idiv)=0
       idivflag(1:8,idiv)=.false.
       idiv=idiv-1
       if(idiv.eq.0)then
        exit outerloop
       endif
       idivflag(ilocate(idiv),idiv)=.true.
      enddo outerloop
      indmapOsec(jindex+1,1)=0
      indmapOsec(jindex+1,2)=indscf
*------------------------------------------------------------------------
      deallocate( jsuflist,jplist,jsidelist,isuflist,iplist,isidelist )
      deallocate( isufflag,isideflag)
*------------------------------------------------------------------------
      return
      end subroutine set_outer

************************************************************************
      subroutine set_ioutersecmap(itet,jindex,idiv,ilocate,numsurf,ns0,
     &     isuflist,indscf)
*
*     Create map for outer sections
*
*     ioctOsec(ind) = 0 : Contains deeper level
*                     = -jindex : outer section ID
*     ioctOsec(ind+1:ind+8) = ind
*     indmapOsec(jindex,1) = Octree index (ex. 82 = section 2-8)
*     indmapOsec(jindex,2) = indscf: Map index for secting surfaces
*     ioutsfmapOsec(indscf+1:indscf+numsurf)
*          : surface indexes intersecting the section
*
*     Last modified  2017/10/20 by T. Furuta
*
************************************************************************
      integer,intent(in) :: itet,jindex,idiv,ilocate(ndiv),numsurf,ns0
      integer,intent(in) :: isuflist(1:ns0)
      integer,intent(inout) :: indscf
      integer :: i,ioutsf,jdiv,jsec
      jsec=0
      do jdiv=1,idiv
       jsec=jsec+ilocate(jdiv)*10**(jdiv-1)
      enddo
      indmapOsec(jindex,1)=jsec
      indmapOsec(jindex,2)=indscf
      i=0
      do i=1,numsurf
       ioutsfmapOsec(indscf+i)=isuflist(i)+nsurf(itet-1)
      enddo
      indscf=indscf+numsurf
      return
      end subroutine set_ioutersecmap

************************************************************************
      subroutine tetrafnd(xx0,xxx,yyy,zzz,uuu,vvv,www,coincd,itet,
     &     ihelem,icl,ierr)
*
*     Find the section or the element at the position x
*
*     ihelem > 0 : if x is located inside the element
*            < 0 : if x is located inside the outer section
*            = 0 : otherwise
*
************************************************************************
      include 'err.inc'

      real(8),intent(in) :: xx0(6),xxx,yyy,zzz,uuu,vvv,www,coincd
      integer,intent(in) :: itet
      integer,intent(out) :: icl,ihelem
      integer,intent(inout) :: ierr
      integer :: isec,ind,idiv,ires,jindex
      integer :: is
      real(8) :: x0(3),x1(3),x(3),u(3)

      ihelem=0
      x(1)=xxx
      x(2)=yyy
      x(3)=zzz
      u(1)=uuu
      u(2)=vvv
      u(3)=www
      if((x(1).lt.xx0(1)-abs(xx0(1))*xtiny)
     &     .or.(x(1).gt.xx0(2)+abs(xx0(2))*xtiny)
     &     .or.(x(2).lt.xx0(3)-abs(xx0(3))*xtiny)
     &     .or.(x(2).gt.xx0(4)+abs(xx0(4))*xtiny)
     &     .or.(x(3).lt.xx0(5)-abs(xx0(5))*xtiny)
     &     .or.(x(3).gt.xx0(6)+abs(xx0(6))*xtiny))then
       ierr=1
       return
      endif
      x0(1)=xx0(1)
      x0(2)=xx0(3)
      x0(3)=xx0(5)
      x1(1:3)=x(1:3)-x0(1:3)
      ind=noctOsec(itet-1)+1
      do idiv=1,ndiv
       x1(1:3)=x1(1:3)-xoctOsec(1:3,idiv,itet)
       isec=1
       do is=1,3
        if(x1(is).gt.0)then
         isec=isec+2**(is-1)
        else
         x1(is)=x1(is)+xoctOsec(is,idiv,itet)
        endif
       enddo
       ind=ioctOsec(ind+isec)
       if(ioctOsec(ind).lt.0)then
        jindex=-ioctOsec(ind)
        exit
       endif
      enddo
      call check_intetras(xx0,x,u,coincd,itet,jindex,ires)
      if(ires.eq.0)then
       call fnd_ihelem(x0,x,u,coincd,itet,ihelem,ierr)
       if(ihelem.gt.0)then
        icl=ielem2icl(ihelem)
       else
        if(iwarning)then
         write(ErrCha,'(''*** TETRA WARNING RECOVERED: '',
     &       ''particle is assumed to be in outer section'')')
          ErrID = 'L:3106/R:tetrafnd/F:tetramod.f' !W80_004_001
          call ErrWrite(ErrID,ErrCha)
        endif
        ihelem=-jindex
       endif
      elseif(ires.gt.0)then
       ihelem=-jindex
      else
       ihelem=-ires
       icl=ielem2icl(ihelem)
      endif
      return
      end subroutine tetrafnd

************************************************************************
      subroutine check_intetras(xx0,x,u,coincd,itet,jindex,ires)
*
*     Check the position x is inside of tetrahedrons
*
*     Last modified 2017/10/20 by T. Furuta
*
************************************************************************
      include 'err.inc'

      real(8),intent(in) :: xx0(6),x(3),u(3),coincd
      integer,intent(in) :: itet,jindex
      integer,intent(out) :: ires
      integer :: i,is,ies,jes,nes,ioutsf,icand,ip(3),m
      integer :: ilocate(ndiv),isec,jsec,idiv,jdiv
      integer :: icheck(3)
      integer,allocatable :: isflag(:)
      integer :: ishift(6,8)
      data ishift/ 0,-1, 0,-1, 0,-1,
     &     1, 0, 0,-1, 0,-1,
     &     0,-1, 1, 0, 0,-1,
     &     1, 0, 1, 0, 0,-1,
     &     0,-1, 0,-1, 1, 0,
     &     1, 0, 0,-1, 1, 0,
     &     0,-1, 1, 0, 1, 0,
     &     1, 0, 1, 0, 1, 0/
      real(8) :: scf(4),t1,t4,t5,d1,dmin
      real(8) :: v(3),w(3),xp(3,3),xc(3),vvec(3,3),xcm(3),v2cm(3)
      real(8) :: xx(6),xcc(3)
      logical :: iflag(3)
      integer :: nflag
      integer :: lgc(4)
*------------------------------------------------------------------------
      nes=indmapOsec(jindex+1,2)-indmapOsec(jindex,2)
      allocate( isflag(nes) )
      ires=0
*------------------------------------------------------------------------
*     No crossing surface to this section and section is outside of
*     tetrahedrons
*------------------------------------------------------------------------
      if(nes.eq.0)then
       ires=2
       return
      endif
*------------------------------------------------------------------------
*     Outside of element block
*------------------------------------------------------------------------
      if( x(1).lt.xb(1,itet).or.x(1).gt.xb(2,itet)
     &     .or.x(2).lt.xb(3,itet).or.x(2).gt.xb(4,itet)
     &     .or.x(3).lt.xb(5,itet).or.x(3).gt.xb(6,itet) )then
       ires=1
       return
      endif
*------------------------------------------------------------------------
*     Check face of surfaces: outside or inside
*------------------------------------------------------------------------
      icheck(1:3)=0
      do ies=1,nes
       ioutsf=ioutsfmapOsec(indmapOsec(jindex,2)+ies)
       if(ioutsf.lt.0)then
        icheck(3)=1 ! Section is inside of a tetrahedron
        exit
       endif
       scf(1:4)=surfcoefs(1:4,ioutsf2surf(ioutsf))
       t4=scf(1)*x(1)+scf(2)*x(2)+scf(3)*x(3)-scf(4)
       t5=scf(1)*u(1)+scf(2)*u(2)+scf(3)*u(3)
       if(abs(t4).lt.coincd*xqt(itet)*abs(t5))then !FURUTA20200622
        d1=-t4/t5
        xc(1:3)=u(1:3)*d1+x(1:3)
        ip(1:3)=ioutsf2outpt(1:3,ioutsf)
        xp(1:3,1)=pointxyz(1:3,ioutpt2point(ip(1)))
        xp(1:3,2)=pointxyz(1:3,ioutpt2point(ip(2)))
        xp(1:3,3)=pointxyz(1:3,ioutpt2point(ip(3)))
        vvec(1:3,1)=xp(1:3,2)-xp(1:3,1)
        vvec(1:3,2)=xp(1:3,3)-xp(1:3,2)
        vvec(1:3,3)=xp(1:3,1)-xp(1:3,3)
        do i=1,3
         w(1:3)=xc(1:3)-xp(1:3,i)
         v(1)=vvec(2,i)*w(3)-vvec(3,i)*w(2)
         v(2)=vvec(3,i)*w(1)-vvec(1,i)*w(3)
         v(3)=vvec(1,i)*w(2)-vvec(2,i)*w(1)
         iflag(i)=scf(1)*v(1)+scf(2)*v(2)+scf(3)*v(3).ge.0.0d0
        enddo
        nflag=count(iflag)
        if(nflag.eq.3.or.nflag.eq.0)then ! All flags are equivalent
         if(t5.ge.0)then ! Judged by velocity vector & a proximity surface
          ires=2
cFURUTA20180517 for non-surface-sharing tetras         exit
         else
          icand=isurf2elem(1,ioutsf2surf(ioutsf))
          call check_inelem(x,u,coincd*xqt(itet),icand,0,ires,lgc)
          if(ires.eq.0)then
           ires=-icand
           exit
          else
           ires=0
          endif
         endif
        endif
        isflag(ies)=0
       endif
       if(t4.gt.coincd*xqt(itet))then ! outside
        isflag(ies)=-2
        icheck(1)=icheck(1)+1
       elseif(t4.gt.0)then ! outside
        isflag(ies)=-1
        icheck(2)=icheck(2)+1
       elseif(t4.le.0)then ! inside
        isflag(ies)=1
        icheck(3)=icheck(3)+1
       endif
      enddo
      if(ires.ne.0)then
       if(iwarning.and.iwarning2)then
        if(ires.lt.0)then

         ErrCha = ''
         ErrID = 'L:3237/R:check_intetras/F:tetramod.f' !W80_005_001
         call ErrWrite(ErrID,ErrCha)

         write(*,'(''*** TETRA WARNING: '',
     &       ''judged INSIDE of an element by a surface '')')
         write(*,'(''ielem ='',i9,'' isurf ='',i9)')
     &           -ires,ioutsf2surf(ioutsf)

        elseif(ires.eq.2)then

         ErrCha = ''
         ErrID = 'L:3248/R:check_intetras/F:tetramod.f' !W80_006_001
         call ErrWrite(ErrID,ErrCha)

         write(*,'(''*** TETRA WARNING: '',
     &       ''judged OUTSIDE by a surface '')')
         write(*,'('' isurf ='',i9)')
     &           ioutsf2surf(ioutsf)
        endif
       endif
       if(ires.gt.0.and.ires.ne.2)then

         ErrCha = ''
         ErrID = 'L:3260/R:check_intetras/F:tetramod.f' !E80_015_001
         call ErrWrite(ErrID,ErrCha)

        write(*,'(''*** TETRA ERROR: '',
     &       ''something wrong in check_intetras'')')
        write(*,'(''ires ='',i9)')
     &         ires
       endif
      elseif(icheck(3).eq.nes)then ! all surfaces are inside
       ires=0
      elseif(icheck(1)+icheck(2).eq.nes)then ! all surfaces are outside
       ires=1
      else
*------------------------------------------------------------------------
*     Check face of nearest surface found in direction from x to CM of
*     a surface facing outside
*------------------------------------------------------------------------
       jsec=indmapOsec(jindex,1)
       ilocate(1:ndiv)=0
       do idiv=1,ndiv
        ilocate(idiv)=mod(jsec,10)
        jsec=jsec/10
        if(jsec.eq.0)exit
       enddo
       xx(1:6)=xx0(1:6)
       do jdiv=1,idiv
        isec=ilocate(jdiv)
        do i=1,6
         is=(i+1)/2
         xx(i)=xx(i)+ishift(i,isec)*xoctOsec(is,jdiv,itet)
        enddo
       enddo
       do ies=1,nes+8
        icand=0
        if(ies.gt.nes)then ! Vertexes of box
         xcm(1)=xx(mod(ies-nes+1,2)+1)
         xcm(2)=xx(mod((ies-nes+1)/2+1,2)+3)
         xcm(3)=xx((ies-nes-1)/4+5)
         v2cm(1:3)=xcm(1:3)-x(1:3)
         dmin=sqrt(v2cm(1)*v2cm(1)+v2cm(2)*v2cm(2)+v2cm(3)*v2cm(3))
         v2cm(1:3)=v2cm(1:3)/dmin
        else ! CM of surfaces
         ioutsf=ioutsfmapOsec(indmapOsec(jindex,2)+ies)
         ip(1:3)=ioutsf2outpt(1:3,ioutsf)
         xp(1:3,1)=pointxyz(1:3,ioutpt2point(ip(1)))
         xp(1:3,2)=pointxyz(1:3,ioutpt2point(ip(2)))
         xp(1:3,3)=pointxyz(1:3,ioutpt2point(ip(3)))
         xcm(1:3)=(xp(1:3,1)+xp(1:3,2)+xp(1:3,3))/3.0d0
         xc(1:3)=xcm(1:3)
         do is=1,4
          if(xcm(1).ge.xx(1).and.xcm(1).le.xx(2)
     &         .and.xcm(2).ge.xx(3).and.xcm(2).le.xx(4)
     &         .and.xcm(3).ge.xx(5).and.xcm(3).le.xx(6))then
           v2cm(1:3)=xcm(1:3)-x(1:3)
           dmin=sqrt(v2cm(1)*v2cm(1)+v2cm(2)*v2cm(2)+v2cm(3)*v2cm(3))
           v2cm(1:3)=v2cm(1:3)/dmin
           icand=ies
           exit
          endif
          if(is.eq.4)exit
          xcm(1:3)=xp(1:3,is)*0.98d0+xc(1:3)*0.02d0 !FURUTA20180517
         enddo
         if(icand.eq.0)cycle
        endif
        do jes=1,nes
         if(jes.eq.ies)cycle
         ioutsf=ioutsfmapOsec(indmapOsec(jindex,2)+jes)
         scf(1:4)=surfcoefs(1:4,ioutsf2surf(ioutsf))
         t1=v2cm(1)*scf(1)+v2cm(2)*scf(2)+v2cm(3)*scf(3)
         if(abs(t1).le.xtiny)cycle
         d1=(scf(4)-x(1)*scf(1)-x(2)*scf(2)-x(3)*scf(3))/t1
         if(d1.lt.0)cycle
         if(d1.lt.dmin+xtiny*xqt(itet))then !FURUTA20180517
          xc(1:3)=v2cm(1:3)*d1+x(1:3)
          ip(1:3)=ioutsf2outpt(1:3,ioutsf)
          xp(1:3,1)=pointxyz(1:3,ioutpt2point(ip(1)))
          xp(1:3,2)=pointxyz(1:3,ioutpt2point(ip(2)))
          xp(1:3,3)=pointxyz(1:3,ioutpt2point(ip(3)))
          xcm(1:3)=(xp(1:3,1)+xp(1:3,2)+xp(1:3,3))/3.0d0 !FURUTA20180517
          xcc(1:3)=1.01d0*xc(1:3)-0.01d0*xcm(1:3) !FURUTA20180517
          vvec(1:3,1)=xp(1:3,2)-xp(1:3,1)
          vvec(1:3,2)=xp(1:3,3)-xp(1:3,2)
          vvec(1:3,3)=xp(1:3,1)-xp(1:3,3)
          do i=1,3
           w(1:3)=xcc(1:3)-xp(1:3,i) !FURUTA20180517
           v(1)=vvec(2,i)*w(3)-vvec(3,i)*w(2)
           v(2)=vvec(3,i)*w(1)-vvec(1,i)*w(3)
           v(3)=vvec(1,i)*w(2)-vvec(2,i)*w(1)
           iflag(i)=scf(1)*v(1)+scf(2)*v(2)+scf(3)*v(3).ge.0.0d0
          enddo
          nflag=count(iflag)
          if(nflag.eq.3.or.nflag.eq.0)then ! All flags are equivalent
cFURUTA20180517------------------------------------------------------
           if(d1.gt.dmin-xtiny*xqt(itet)
     &          .and.d1.lt.dmin+xtiny*xqt(itet))then
            icand=0
            exit ! Exception for non-sharing surfaces
           elseif(d1.lt.dmin)then
            icand=jes
            dmin=d1
           endif
          else
           xcc(1:3)=0.99d0*xc(1:3)+0.01d0*xcm(1:3)
           do i=1,3
            w(1:3)=xcc(1:3)-xp(1:3,i)
            v(1)=vvec(2,i)*w(3)-vvec(3,i)*w(2)
            v(2)=vvec(3,i)*w(1)-vvec(1,i)*w(3)
            v(3)=vvec(1,i)*w(2)-vvec(2,i)*w(1)
            iflag(i)=scf(1)*v(1)+scf(2)*v(2)+scf(3)*v(3).ge.0.0d0
           enddo
           nflag=count(iflag)
           if(nflag.eq.3.or.nflag.eq.0)then
            icand=0 ! Exception for surface edge
            exit
           endif
          endif
c--------------------------------------------------------------------
         endif
        enddo
        if(icand.gt.0)then
         if(isflag(icand).eq.0)then
          icand=0
         else
          if(isflag(icand).lt.0)then
           ires=1
          else
           ires=0
          endif
          exit
         endif
        endif
       enddo
       if(icand.eq.0)then
        if(iwarning)then
         write(ErrCha,'(''*** TETRA WARNING: '',
     &       ''fail to check inside tetrahedrons'')')
          ErrID = 'L:3396/R:check_intetras/F:tetramod.f' !W80_007_001
          call ErrWrite(ErrID,ErrCha)
        endif
        ires=0 ! Check inside tetrahedrons
       endif
      endif
      return
      end subroutine check_intetras

************************************************************************
      subroutine fnd_ihelem(x0,x,u,coincd,itet,ihelem,ierr)
*
*     Find the element at the position x
*
*     Last modified 2016/09/08 by T. Furuta
*
************************************************************************
      include 'err.inc'

      real(8),intent(in) :: x0(3),x(3),u(3),coincd
      integer,intent(in) :: itet
      integer,intent(out) :: ihelem
      integer,intent(inout) :: ierr
      integer :: i,ii,is,ires,icand,ielem,jelem,j,ip(4)
      integer :: mm,isf(4),ies,icand0,ne0,lgc(4)
      real(8) :: xmin,x2,dmin,dist,xsum(3),x1(3)
      integer :: iblk,ind,ind0,idiv,jdiv,isec,indmap(ndiv)
      logical :: idivflag(8,ndiv)
      logical,allocatable :: ielemflag(:)
      integer,allocatable :: icandlist(:)
*------------------------------------------------------------------------
*     Find the nearest block to the position x containing more than
*     one element
*------------------------------------------------------------------------
      ihelem=0
      x1(1)=x(1)-xb(1,itet)
      x1(2)=x(2)-xb(3,itet)
      x1(3)=x(3)-xb(5,itet)
      ind=noctEblk(itet-1)+1
      do idiv=1,ndiv
       x1(1:3)=x1(1:3)-xoctEblk(1:3,idiv,itet)
       isec=1
       do is=1,3
        if(x1(is).gt.0)then
          isec=isec+2**(is-1)
        else
          x1(is)=x1(is)+xoctEblk(is,idiv,itet)
        endif
       enddo
       ind0=ind
       ind=ioctEblk(ind+isec)
       if(ioctEblk(ind).lt.0)then
         iblk=-ioctEblk(ind)
         exit
       endif
      enddo
*------------------------------------------------------------------------
*     No element in the block at x
*     Find a block which contains elements from same level of division
*------------------------------------------------------------------------
      if(ielemmapEblk(indmapEblk(iblk,2)+1).lt.0)then
        indmap(idiv)=ind0
        idivflag(1:8,idiv)=.false.
        idivflag(isec,idiv)=.true.
        jdiv=idiv
        outerloop: do i=1,numlarge
         isecloop: do isec=1,8
          ind=indmap(jdiv)
          if(.not.idivflag(isec,jdiv))then
            ind=ioctEblk(ind+isec)
            if(ioctEblk(ind).lt.0)then
              iblk=-ioctEblk(ind)
              if(ielemmapEblk(indmapEblk(iblk,2)+1).gt.0)exit outerloop
              idivflag(isec,jdiv)=.true.
            else
              jdiv=jdiv+1
              indmap(jdiv)=ind
              idivflag(1:8,jdiv)=.false.
              cycle outerloop
            endif
          endif
         enddo isecloop
         indmap(jdiv)=0
         jdiv=jdiv-1
         if(jdiv.lt.idiv)then
          ErrCha = ''
          ErrID = 'L:3482/R:fnd_ihelem/F:tetramod.f' !E80_016_001
          call ErrWrite(ErrID,ErrCha)

          write(*,'(''*** TETRA ERROR: '',
     &       ''something wrong in fnd_ihelem'')')
          write(*,'(''jdiv ='',i5,'' idiv ='',i5)')jdiv,idiv
          ierr=1
          return
         endif
        enddo outerloop
        if(i.ge.numlarge)then

          ErrCha = ''
          ErrID = 'L:3495/R:fnd_ihelem/F:tetramod.f' !E80_016_002
          call ErrWrite(ErrID,ErrCha)

          write(*,'(''*** TETRA ERROR: '',
     &       ''something wrong in fnd_ihelem'')')
          write(*,'(''i ='',i10,'' numlarge ='',i10)')i,numlarge
          ierr=1
          return
        endif
      endif
*------------------------------------------------------------------------
*     Find the element containing the position x among elements in block
*------------------------------------------------------------------------
      xmin=1.0d20
      icand=0
      do jelem=indmapEblk(iblk,2)+1,indmapEblk(iblk+1,2)
       ielem=ielemmapEblk(jelem)
       mm=0
       call check_inelem(x,u,coincd*xqt(itet),ielem,mm,ires,lgc)
       if(ires.eq.0)then
         ihelem=ielem
         exit
       endif
       ip(1:4)=ielem2point(1:4,ielem)
       xsum(1:3)=-4*x(1:3)
       do j=1,4
        xsum(1:3)=xsum(1:3)+pointxyz(1:3,ip(j))
       enddo
       xsum(1:3)=0.25d0*xsum(1:3)
       x2=xsum(1)*xsum(1)+xsum(2)*xsum(2)+xsum(3)*xsum(3)
       if(x2.lt.xmin)then
        icand=ielem
        xmin=x2
       endif
      enddo
      if(ihelem.eq.0.and.icand.eq.0)then

       ErrCha = ''
       ErrID = 'L:3533/R:fnd_ihelem/F:tetramod.f' !E80_017_001
       call ErrWrite(ErrID,ErrCha)

       write(*,'(''*** TETRA ERROR: '',
     &      ''nearest element not found'')')
       write(*,'(''x2 ='',e15.3,'' xmin ='',e15.3)')x2,xmin
       ierr=1
       return
      endif
*------------------------------------------------------------------------
*     Find the element containing the position x
*------------------------------------------------------------------------
      if(ihelem.eq.0)then
*------------------------------------------------------------------------
        ne0=nelem(itet)-nelem(itet-1)
        allocate( ielemflag(nelem(itet-1)+1:nelem(itet)) )
        allocate( icandlist(ne0) )
*------------------------------------------------------------------------
        if(iwarning)then

         ErrCha = ''
         ErrID = 'L:3554/R:fnd_ihelem/F:tetramod.f' !W80_008_001
         call ErrWrite(ErrID,ErrCha)

          write(*,'(''*** TETRA WARNING: '',
     &        ''element could not find in box'')')
          write(*,'(''iblk ='',i9,'' Nelem contained ='',i9)')
     &         iblk,indmapEblk(iblk+1,2)-indmapEblk(iblk,2)
        endif
        ielemflag(nelem(itet-1)+1:nelem(itet))=.false.
        icandlist(1:ne0)=0
        icandlist(1)=icand
        ielemflag(icand)=.true.
        i=1
        j=1
        iiloop2: do ii=1,numlarge
         icand=icandlist(i)
         mm=0
         call check_inelem(x,u,coincd*xqt(itet),icand,mm,ires,lgc)
         if(ires.eq.0)then
           if(iwarning)then

            ErrCha = ''
            ErrID = 'L:3576/R:fnd_ihelem/F:tetramod.f' !W80_009_001
            call ErrWrite(ErrID,ErrCha)

             write(*,'(''*** TETRA WARNING RECOVERED: '',
     &           ''element found '')')
             write(*,'(''i ='',i9)')icand
           endif
           exit iiloop2
         endif
         isf(1:4)=abs(ielem2surf(1:4,icand))
         iresloop: do ies=1,4
          call fnd_nextelem(x0,x,u,coincd,itet,isf(ies),icand,
     &        icand0,0,ierr)
          if(icand0.gt.0)then
            if(.not.ielemflag(icand0).and.j.lt.ne0)then
              j=j+1
              icandlist(j)=icand0
              ielemflag(icand0)=.true.
            endif
          endif
         enddo iresloop
         i=i+1
         if(i.gt.j)then
           if(iwarning)then

           ErrCha = ''
           ErrID = 'L:3602/R:fnd_ihelem/F:tetramod.f' !W80_010_001
           call ErrWrite(ErrID,ErrCha)

             write(*,'(''*** TETRA WARNING: '',
     &        ''element at the position x not found'')')
             write(*,'(''i ='',i9,'' nelem ='',i9)')i,ne0
           endif
           do icand=nelem(itet-1)+1,nelem(itet)
            if(.not.ielemflag(icand))then
             call check_inelem(x,u,coincd*xqt(itet),icand,0,ires,lgc)
              if(ires.eq.0)then
                if(iwarning)then

                 ErrCha = ''
                 ErrID = 'L:3616/R:fnd_ihelem/F:tetramod.f' !W80_009_002
                 call ErrWrite(ErrID,ErrCha)

                 write(*,'(''*** TETRA WARNING RECOVERED: '',
     &           ''element found '')')
                 write(*,'(''i ='',i9)')icand
                endif
                exit iiloop2
              endif
            endif
           enddo
           if(iwarning)then

            ErrCha = ''
            ErrID = 'L:3630/R:fnd_ihelem/F:tetramod.f' !W80_010_002
            call ErrWrite(ErrID,ErrCha)

            write(*,'(''*** TETRA WARNING: '',
     &       ''element at the position x not found'')')
            write(*,'(''i ='',i9,'' nelem ='',i9)')i,
     &           nelem(itet)-nelem(itet-1)
           endif
           return
         endif
        enddo iiloop2
        if(ii.ge.numlarge)then
          if(iwarning)then

           ErrCha = ''
           ErrID = 'L:3645/R:fnd_ihelem/F:tetramod.f' !W80_010_003
           call ErrWrite(ErrID,ErrCha)

           write(*,'(''*** TETRA WARNING: '',
     &      ''element at the position x not found'')')
           write(*,'(''x,y,z ='',3f10.5)') x(1:3)
          endif
          return
        endif
        ihelem=icand
*------------------------------------------------------------------------
        deallocate( ielemflag,icandlist )
*------------------------------------------------------------------------
      endif
      return
      end subroutine fnd_ihelem

************************************************************************
      subroutine check_inelem(x,u,xsmall,icand,mm,ires,lgc)
*
*     Check the position x is inside the element
*
*     mm=0
*       ires = -1 : the section is inside tetrahedrons
*       ires =  0 : inside
*       ires >  0 : 1st bad surface
*     mm=2
*       ires =  0 : inside
*       ires >  0 : 1st bad surface
*       lgc(ies) = 0 : outside
*                = 1 : inside
*
************************************************************************
      real(8),intent(in) :: x(3),u(3),xsmall
      integer,intent(in) :: icand,mm
      integer,intent(out) :: ires,lgc(4)
      integer :: isf(4),ies,isurf
      real(8) :: scf(4),t4,t5
      ires=0
      isf(1:4)=ielem2surf(1:4,icand)
      do ies=1,4
       isurf=abs(isf(ies))
       scf(1:4)=surfcoefs(1:4,isurf)
       t4=scf(1)*x(1)+scf(2)*x(2)+scf(3)*x(3)-scf(4)
       t5=scf(1)*u(1)+scf(2)*u(2)+scf(3)*u(3)
       lgc(ies)=1
       if(abs(t4).le.xsmall*abs(t5))t4=t5
       if(isf(ies)*t4.lt.0)lgc(ies)=0
       if(lgc(ies).eq.mm)then
        ires=ies
        exit
       endif
      enddo
      if(mm.eq.2)then
       do ies=1,4
        ires=ires+1-lgc(ies)
       enddo
      endif
      return
      end subroutine check_inelem

************************************************************************
      subroutine check_inelem2(x,xsmall,icand,ires)
*
*     Check the position x is inside the element
*
*       ires =  0 : inside
*       ires >  0 : 1st bad surface
*
************************************************************************
      real(8),intent(in) :: x(3),xsmall
      integer,intent(in) :: icand
      integer,intent(out) :: ires
      integer :: isf(4),ies,isurf
      real(8) :: scf(4),t4
      ires=0
      isf(1:4)=ielem2surf(1:4,icand)
      do ies=1,4
       isurf=abs(isf(ies))
       scf(1:4)=surfcoefs(1:4,isurf)
       t4=scf(1)*x(1)+scf(2)*x(2)+scf(3)*x(3)-scf(4)
       if(abs(t4).le.xsmall)then
        ires=ies ! No judgment for close surface
        exit
       endif
       if(isf(ies)*t4.lt.0)then
        ires=ies
        exit
       endif
      enddo
      return
      end subroutine check_inelem2

************************************************************************
      subroutine tetracald(xx0,x,u,itet,ihelem,dls,jsu,ierr)
*
*     Calculate the distance to the next crossing surface in sections
*
*     dls : distance
*     jsu : surface index of the next crossing surface
*
************************************************************************
      include 'err.inc'

      real(8),intent(in) :: xx0(6),x(3),u(3)
      integer,intent(in) :: itet,ihelem
      real(8),intent(out) :: dls
      integer,intent(out) :: jsu
      integer,intent(inout) :: ierr
      if(ihelem.lt.0)then
       call calc_outersecdist(xx0,x,u,itet,ihelem,dls,jsu,ierr)
      elseif(ihelem.gt.0)then
       call calc_elemdist(x,u,ihelem,dls,jsu,ierr)
      else
       ierr=1
      endif
      if(ierr.gt.0)then

       ErrCha = ''
       ErrID = 'L:3764/R:tetracald/F:tetramod.f' !E80_018_001
       call ErrWrite(ErrID,ErrCha)

       write(*,'(''*** TETRA ERROR: '',
     &      ''fail to calculate distance. '',
     &      '' Neither in outersec nor element'')')
       write(*,'(''x,y,z ='',3f10.5)') x(1:3)
       write(*,'(''ihelem='',i9)') ihelem

      endif
      return
      end subroutine tetracald

************************************************************************
      subroutine calc_outersecdist(xx0,x,u,itet,ihelem,
     &     dls,jsu,ierr)
*
*     Calculate the distance to the next crossing surface in outer
*     sections
*
*     dls : distance
*     jsu : surface index of the next crossing surface
*
************************************************************************
      include 'err.inc'

      real(8),intent(in) :: xx0(6),x(3),u(3)
      integer,intent(in) :: itet,ihelem
      real(8),intent(out) :: dls
      integer,intent(out) :: jsu
      integer,intent(inout) :: ierr
      integer :: ilocate(ndiv)
      integer :: ishift(6,8)
      data ishift/ 0,-1, 0,-1, 0,-1,
     &     1, 0, 0,-1, 0,-1,
     &     0,-1, 1, 0, 0,-1,
     &     1, 0, 1, 0, 0,-1,
     &     0,-1, 0,-1, 1, 0,
     &     1, 0, 0,-1, 1, 0,
     &     0,-1, 1, 0, 1, 0,
     &     1, 0, 1, 0, 1, 0/
      integer :: ies,isf(4),icand
      integer :: ip(3),i,is,isurf,ioutsf,idiv,jdiv,isec,jsec
      integer :: isig,nes,ihsec,nflag
      logical :: iflag(3)
      real(8) :: scf(4),d1,t1,dmin
      real(8) :: xx(6),xp(3,3),xc(3),v(3),w(3),vvec(3,3)
*------------------------------------------------------------------------
      ihsec=-ihelem
*------------------------------------------------------------------------
*     Setting the border of the section xx
*------------------------------------------------------------------------
      jsec=indmapOsec(ihsec,1)
      ilocate(1:ndiv)=0
      do idiv=1,ndiv
       ilocate(idiv)=mod(jsec,10)
       jsec=jsec/10
       if(jsec.eq.0)exit
      enddo
      xx(1:6)=xx0(1:6)
      do jdiv=1,idiv
       isec=ilocate(jdiv)
       do isurf=1,6
        is=(isurf+1)/2
        xx(isurf)=xx(isurf)+ishift(isurf,isec)*xoctOsec(is,jdiv,itet)
       enddo
      enddo
*------------------------------------------------------------------------
*     Calculate the distance to the borders
*------------------------------------------------------------------------
      dmin=1.0d20
      icand=0
      do isurf=1,6
       is=(isurf+1)/2
       isig=(-1)**mod(isurf,2)
       if(isig*u(is).gt.0.0d0)then
        d1=(xx(isurf)-x(is))/u(is)
        if(d1.gt.0.0d0)then
         if(d1.lt.dmin)then
          icand=isurf
          dmin=d1
         endif
        endif
       endif
      enddo
*------------------------------------------------------------------------
*     Calculate the distance to the intersecting triangles
*------------------------------------------------------------------------
      nes=indmapOsec(ihsec+1,2)-indmapOsec(ihsec,2)
      do ies=1,nes
       ioutsf=ioutsfmapOsec(indmapOsec(ihsec,2)+ies)
       if(ioutsf.lt.0)then

        ErrCha = ''
        ErrID = 'L:3858/R:calc_outersecdist/F:tetramod.f' !E80_019_001
        call ErrWrite(ErrID,ErrCha)

        write(*,'(''*** TETRA ERROR: '',
     &       ''section is inside tetrahedron'')')

        write(*,'(''ihsec ='',i9,'' jsec ='',i9)')
     &       ihsec,indmapOsec(ihsec,1)

        ierr=1
        exit
        ! Section is inside tetrahedron
       endif
       scf(1:4)=surfcoefs(1:4,ioutsf2surf(ioutsf))
       t1=u(1)*scf(1)+u(2)*scf(2)+u(3)*scf(3)
       if(t1.lt.0.0d0)then ! normal vector is defined outward
        d1=(scf(4)-x(1)*scf(1)-x(2)*scf(2)-x(3)*scf(3))/t1
        if(d1.gt.0.0d0)then
         if(d1.le.dmin)then
          xc(1:3)=u(1:3)*d1+x(1:3)
          ip(1:3)=ioutsf2outpt(1:3,ioutsf)
          xp(1:3,1)=pointxyz(1:3,ioutpt2point(ip(1)))
          xp(1:3,2)=pointxyz(1:3,ioutpt2point(ip(2)))
          xp(1:3,3)=pointxyz(1:3,ioutpt2point(ip(3)))
          vvec(1:3,1)=xp(1:3,2)-xp(1:3,1)
          vvec(1:3,2)=xp(1:3,3)-xp(1:3,2)
          vvec(1:3,3)=xp(1:3,1)-xp(1:3,3)
          do i=1,3
           w(1:3)=xc(1:3)-xp(1:3,i)
           v(1)=vvec(2,i)*w(3)-vvec(3,i)*w(2)
           v(2)=vvec(3,i)*w(1)-vvec(1,i)*w(3)
           v(3)=vvec(1,i)*w(2)-vvec(2,i)*w(1)
           iflag(i)=scf(1)*v(1)+scf(2)*v(2)+scf(3)*v(3).ge.0.0d0
          enddo
          nflag=count(iflag)
          if(nflag.eq.3.or.nflag.eq.0)then ! All flags are equivalent
           icand=ioutsf+6
           dmin=d1
          endif
         endif
        endif
       endif
      enddo
      if(icand.eq.0)ierr=1
      dls=dmin
      jsu=icand
      return
      end subroutine calc_outersecdist

************************************************************************
      subroutine calc_elemdist(x,u,ihelem,dls,jsu,ierr)
*
*     Calculate the distance to the next crossing surface in elements
*
*     dls : distance
*     jsu : surface index of the next crossing surface
*
************************************************************************
      real(8),intent(in) :: x(3),u(3)
      integer,intent(in) :: ihelem
      real(8),intent(out) :: dls
      integer,intent(out) :: jsu
      integer,intent(inout) :: ierr
      integer :: ies,isf(4),icand,isurf,jsec
      real(8) :: scf(4),d1,t1,dmin
      isf(1:4)=ielem2surf(1:4,ihelem)
      dmin=1.0d20
      icand=0
      do ies=1,4
       isurf=abs(isf(ies))
       scf(1:4)=surfcoefs(1:4,isurf)
       t1=u(1)*scf(1)+u(2)*scf(2)+u(3)*scf(3)
       if(abs(t1).gt.xtiny)then
        d1=(scf(4)-x(1)*scf(1)-x(2)*scf(2)-x(3)*scf(3))/t1
        if(d1.gt.0.0d0.and.t1*isf(ies).lt.0.0d0)then
         if(d1.lt.dmin)then
          icand=isurf
          dmin=d1
         endif
        endif
       endif
      enddo
      if(icand.eq.0)ierr=1
      dls=dmin
      jsu=icand
      return
      end subroutine calc_elemdist

************************************************************************
      subroutine tetranext(x0,xxx,yyy,zzz,uuu,vvv,www,coincd,itet,
     &     jsu,ihelem,nxtelem,iap,ierr)
*
*     Find the next section
*
************************************************************************
      include 'err.inc'

      real(8),intent(in) :: x0(3),xxx,yyy,zzz,uuu,vvv,www,coincd
      integer,intent(in) :: itet,jsu,ihelem
      integer,intent(out) :: nxtelem
      integer,intent(inout) :: iap,ierr
      real(8) :: x(3),u(3)

      x(1)=xxx
      x(2)=yyy
      x(3)=zzz
      u(1)=uuu
      u(2)=vvv
      u(3)=www
      if(ihelem.lt.0)then
       call fnd_nextoutersec(x0,x,u,coincd,itet,jsu,ihelem,
     &     nxtelem,ierr)
      elseif(ihelem.gt.0)then
       call fnd_nextelem(x0,x,u,coincd,itet,jsu,ihelem,
     &     nxtelem,1,ierr)
      else

       ErrCha = ''
       ErrID = 'L:3976/R:tetranext/F:tetramod.f' !E80_018_002
       call ErrWrite(ErrID,ErrCha)

       write(*,'(''*** TETRA ERROR: '',
     &      ''fail to find next section. '',
     &      ''Neither in outersec nor element.'')')
       write(*,'(''x0,y0,z0 ='',3f10.5)') x0(1:3)
       write(*,'(''x ,y ,z  ='',3f10.5)') x(1:3)
       write(*,'(''ihsec,ihelem='',2i5)') ihelem
       ierr=1
      endif
      if(nxtelem.gt.0)then
       iap=ielem2icl(nxtelem)
      endif
      return
      end subroutine tetranext

************************************************************************
      subroutine fnd_nextoutersec(x0,x,u,coincd,itet,jsu,ihelem,
     &     nxtelem,ierr)
*
*     Find the next outersec
*
************************************************************************
      include 'err.inc'

      real(8),intent(in) :: x0(3),x(3),u(3),coincd
      integer,intent(in) :: itet,jsu,ihelem
      integer,intent(out) :: nxtelem
      integer,intent(inout) :: ierr
      integer :: i,isec,jsec,ies,nes,ioutsf,ind,iskip,isig
      integer :: ilocate(ndiv),imap(ndiv),idiv,jdiv,is,icount
      integer :: isf,ksec(4),isurf,ihsec,nxtsec
      integer :: ishift(6,8)
      data ishift/ 0,-1, 0,-1, 0,-1,
     &     1, 0, 0,-1, 0,-1,
     &     0,-1, 1, 0, 0,-1,
     &     1, 0, 1, 0, 0,-1,
     &     0,-1, 0,-1, 1, 0,
     &     1, 0, 0,-1, 1, 0,
     &     0,-1, 1, 0, 1, 0,
     &     1, 0, 1, 0, 1, 0/
      real(8) :: x1(3)
      ihsec=-ihelem
      nxtelem=0
      if(jsu.le.6)then
*------------------------------------------------------------------------
*     Expand the section index
*------------------------------------------------------------------------
       jsec=indmapOsec(ihsec,1)
       ilocate(1:ndiv)=0
       do idiv=1,ndiv
        ilocate(idiv)=mod(jsec,10)
        jsec=jsec/10
        if(jsec.eq.0)exit
       enddo
*------------------------------------------------------------------------
       do jdiv=idiv,1,-1 !FURUTA20181116 bugfix
        isec=ilocate(jdiv)
        if(ishift(jsu,isec).ne.0)exit ! Find level able to cross surface
       enddo
       idiv=jdiv
*------------------------------------------------------------------------
       if(idiv.eq.0)return
*------------------------------------------------------------------------
       iskip=(jsu+1)/2
       iskip=2**(iskip-1)
       if(mod(jsu,2).eq.1)then
        isig=-1
       else
        isig=1
       endif
       ind=noctOsec(itet-1)+1
       do jdiv=1,idiv-1
        ind=ioctOsec(ind+ilocate(jdiv))
       enddo
       ind=ioctOsec(ind+isec+isig*iskip)
       ilocate(idiv)=isec+isig*iskip
       if(ioctOsec(ind).eq.0)then
        x1(1:3)=x(1:3)-x0(1:3)
     &       +coincd*xoctOsec(1:3,ndiv,itet)*u(1:3)     ! Tiny shift to u vector
        do jdiv=1,idiv
         jsec=ilocate(jdiv)
         do is=1,3
          x1(is)=x1(is)-ishift(2*is-1,jsec)*xoctOsec(is,jdiv,itet)
         enddo
        enddo
        do jdiv=idiv+1,ndiv
         x1(1:3)=x1(1:3)-xoctOsec(1:3,jdiv,itet)
         isec=1
         do is=1,3
          if(x1(is).gt.0)then
           isec=isec+2**(is-1)
          else
           x1(is)=x1(is)+xoctOsec(is,jdiv,itet)
          endif
         enddo
         ind=ioctOsec(ind+isec)
         if(ioctOsec(ind).lt.0)then
          nxtsec=-ioctOsec(ind)
          jsec=indmapOsec(nxtsec,1)
          ilocate(1:ndiv)=0
          do idiv=1,ndiv
           ilocate(idiv)=mod(jsec,10)
           jsec=jsec/10
           if(jsec.eq.0)exit
          enddo
          exit
         endif
        enddo
       elseif(ioctOsec(ind).lt.0)then
        nxtsec=-ioctOsec(ind)
       else

        ErrCha = ''
        ErrID = 'L:4091/R:fnd_nextoutersec/F:tetramod.f' !E80_020_001
        call ErrWrite(ErrID,ErrCha)

        write(*,'(''*** TETRA ERROR: '',
     &       ''contradiction in octreemap'')')
        write(*,'(''jsu ='',i9,'' jindex ='',i9)')
     &       jsu,indmapOsec(ihsec,1)
        ierr=1
       endif
       nxtelem=-nxtsec
*------------------------------------------------------------------------
*     Find the next element
*------------------------------------------------------------------------
      else
       isurf=ioutsf2surf(jsu-6)
       nxtelem=isurf2elem(1,isurf)
       if(isurf2elem(2,isurf).ne.0)then

        ErrCha = ''
        ErrID = 'L:4110/R:fnd_nextoutersec/F:tetramod.f' !E80_021_001
        call ErrWrite(ErrID,ErrCha)

        write(*,'(''*** TETRA ERROR: '',
     &       ''table for isurf2elem has problem'')')
        write(*,
     &    '(''jsu ='',i9,'' ielem ='',i9,'' isurf2elem ='',2i9)')
     &       jsu,nxtelem,isurf2elem(1:2,jsu)
        ierr=1
       endif
      endif
*------------------------------------------------------------------------
      return
      end subroutine fnd_nextoutersec

************************************************************************
      subroutine fnd_nextelem(x0,x,u,coincd,itet,jsu,ihelem,
     &     nxtelem,imode,ierr)
*
*     Find the next element
*
************************************************************************
      include 'err.inc'

      real(8),intent(in) :: x0(3),x(3),u(3),coincd
      integer,intent(in) :: itet,jsu,ihelem,imode
      integer,intent(out) :: nxtelem
      integer,intent(inout) :: ierr
      integer :: ind,idiv,isec,is
      real(8) :: x1(3)
      if(isurf2elem(1,jsu).eq.ihelem)then
       nxtelem=isurf2elem(2,jsu)
      else
       nxtelem=isurf2elem(1,jsu)
       if(isurf2elem(2,jsu).ne.ihelem)then

        ErrCha = ''
        ErrID = 'L:4147/R:fnd_nextelem/F:tetramod.f' !E80_021_002
        call ErrWrite(ErrID,ErrCha)

        write(*,'(''*** TETRA ERROR: '',
     &       ''table for isurf2elem has a problem'')')
        write(*,
     &  '(''jsu ='',i9,'' ielem ='',i9,'' isurf2elem ='',2i9)')
     &       jsu,ihelem,isurf2elem(1:2,jsu)
        ierr=1
       endif
      endif
*------------------------------------------------------------------------
*     Find the next outersec
*------------------------------------------------------------------------
      if(nxtelem.eq.0.and.imode.eq.1)then
       x1(1:3)=x(1:3)-x0(1:3)
     &      +coincd*xoctOsec(1:3,ndiv,itet)*u(1:3) ! Tiny shift to u vector
       ind=noctOsec(itet-1)+1
       do idiv=1,ndiv
        x1(1:3)=x1(1:3)-xoctOsec(1:3,idiv,itet)
        isec=1
        do is=1,3
         if(x1(is).gt.0)then
          isec=isec+2**(is-1)
         else
          x1(is)=x1(is)+xoctOsec(is,idiv,itet)
         endif
        enddo
        ind=ioctOsec(ind+isec)
        if(ioctOsec(ind).lt.0)then
         nxtelem=ioctOsec(ind)
         exit
        endif
       enddo
      endif
      end subroutine fnd_nextelem

************************************************************************
      subroutine tetraangl(jsu,ihelem,uuu,vvv,www,cs)
*
*     Calculate angle to surface
*
************************************************************************
      include 'err.inc'

      integer,intent(in) :: jsu,ihelem
      real(8),intent(in) :: uuu,vvv,www
      real(8),intent(out) :: cs
      if(ihelem.gt.0)then
       call calc_anglelem(jsu,uuu,vvv,www,cs)
      elseif(ihelem.lt.0)then
       call calc_anglsec(jsu,uuu,vvv,www,cs)
      else

       ErrCha = ''
       ErrID = 'L:4202/R:tetraangl/F:tetramod.f' !E80_022_001
       call ErrWrite(ErrID,ErrCha)

       write(*,'(''*** TETRA ERROR: '',
     &      ''angle to surface is not obtained'')')
       write(*,'(''jsu ='',i9,'' ihelem ='',i9)')
     &      jsu,ihelem
       cs=1.0d0
      endif
      end subroutine tetraangl

************************************************************************
      subroutine calc_anglelem(jsu,uuu,vvv,www,cs)
*
*     Calculate angle to element surface
*
************************************************************************
      integer,intent(in) :: jsu
      real(8),intent(in) :: uuu,vvv,www
      real(8),intent(out) :: cs
      real(8) :: scf(4),ang(3),t4
      real(8),parameter :: one=1.0d0
      scf(1:4)=surfcoefs(1:4,jsu)
      ang(1:3)=scf(1:3)
      t4=sqrt(ang(1)**2+ang(2)**2+ang(3)**2)
      if(t4.ne.0)t4=1.0d0/t4
      ang(1:3)=ang(1:3)*t4
      cs=max(-one,min(one,ang(1)*uuu+ang(2)*vvv+ang(3)*www))
      return
      end subroutine calc_anglelem

************************************************************************
      subroutine calc_anglsec(jsu,uuu,vvv,www,cs)
*
*     Calculate angle to section surface or outer surface
*
************************************************************************
      integer,intent(in) :: jsu
      real(8),intent(in) :: uuu,vvv,www
      real(8),intent(out) :: cs
      integer :: isurf
      real(8) :: scf(4),ang(3),t4
      real(8),parameter :: one=1.0d0
      if(jsu.le.6)then
       select case(jsu)
       case(1)
        cs=uuu
       case(2)
        cs=uuu
       case(3)
        cs=vvv
       case(4)
        cs=vvv
       case(5)
        cs=www
       case(6)
        cs=www
       end select
      else
       isurf=ioutsf2surf(jsu-6)
       scf(1:4)=surfcoefs(1:4,isurf)
       ang(1:3)=scf(1:3)
       t4=sqrt(ang(1)**2+ang(2)**2+ang(3)**2)
       if(t4.ne.0)t4=1.0d0/t4
       ang(1:3)=ang(1:3)*t4
       cs=max(-one,min(one,ang(1)*uuu+ang(2)*vvv+ang(3)*www))
      endif
      return
      end subroutine calc_anglsec

************************************************************************
      subroutine tetrachk(xx0,x,u,coincd,itet,icl,ihelem,ires,ierr)
*
*     Check consistency of tetra and current position
*
************************************************************************
      integer,intent(in) :: itet
      integer,intent(inout) :: ihelem
      real(8),intent(in) :: xx0(6),x(3),u(3),coincd
      integer,intent(out) :: icl,ires,ierr
      integer :: ihelem0

      ihelem0=ihelem
      if(ihelem.gt.0)then
       call check_elem(xx0,x,u,coincd,itet,icl,ihelem,ierr)
      else
       call check_sec(xx0,x,u,coincd,itet,icl,ihelem,ierr)
      endif
      if(ihelem.eq.ihelem0)then
       ires=0
      else
       ires=1
       if(ihelem.gt.0)then
        icl=ielem2icl(ihelem)
       endif
       if(iwarning)call tetrawarning(ihelem0,ihelem)
      endif
      return
      end subroutine tetrachk

************************************************************************
      subroutine check_elem(xx0,x,u,coincd,itet,icl,ihelem,ierr)
*
*     Check consistency of element and current position
*
************************************************************************
      integer,intent(in) :: itet
      integer,intent(inout) :: ihelem
      real(8),intent(in) :: xx0(6),x(3),u(3),coincd
      integer,intent(out) :: icl,ierr
      integer :: ires,ies,ie,isf(1:4),icand,ihelem0
      integer :: lgc(4),lgc2(4)
      real(8) :: xxx,yyy,zzz,uuu,vvv,www
      ihelem0=ihelem
      call check_inelem(x,u,coincd*xqt(itet),ihelem,2,ires,lgc)

      if(ires.gt.0)then
*------------------------------------------------------------------------
*     Check neibouring TETRAs
*------------------------------------------------------------------------
       isf(1:4)=ielem2surf(1:4,ihelem)
       do ies=1,4
        if(lgc(ies).eq.0)then
         do ie=1,2
          if(isurf2elem(ie,abs(isf(ies))).ne.ihelem)
     &         icand=isurf2elem(ie,abs(isf(ies)))
         enddo
         if(icand.gt.0)then
          call check_inelem(x,u,coincd*xqt(itet),icand,0,ires,lgc2)
          if(ires.eq.0)then
           ihelem=icand
           exit
          endif
         endif
        endif
       enddo
*------------------------------------------------------------------------
       if(ihelem.eq.ihelem0)then
*------------------------------------------------------------------------
*     For cases not found
*------------------------------------------------------------------------
        xxx=x(1)
        yyy=x(2)
        zzz=x(3)
        uuu=u(1)
        vvv=u(2)
        www=u(3)
        call tetrafnd(xx0,xxx,yyy,zzz,uuu,vvv,www,coincd,itet,
     &       ihelem,icl,ierr)
       endif
*------------------------------------------------------------------------
      endif
      return
      end subroutine check_elem

************************************************************************
      subroutine check_sec(xx0,x,u,coincd,itet,icl,ihelem,ierr)
*
*     Check consistency of section and current position
*
************************************************************************
      integer,intent(in) :: itet
      integer,intent(inout) :: ihelem
      real(8),intent(in) :: xx0(6),x(3),u(3),coincd
      integer,intent(out) :: icl,ierr
      integer :: ihsec,isec,is,ies,isurf,jsec,idiv,jdiv,ioutsf,icand
      integer :: ires,ilocate(ndiv),nes
      real(8) :: xx(6),xxx,yyy,zzz,uuu,vvv,www,t4,t5
      logical :: outside
      integer :: ishift(6,8),lgc(4)
      data ishift/ 0,-1, 0,-1, 0,-1,
     &     1, 0, 0,-1, 0,-1,
     &     0,-1, 1, 0, 0,-1,
     &     1, 0, 1, 0, 0,-1,
     &     0,-1, 0,-1, 1, 0,
     &     1, 0, 0,-1, 1, 0,
     &     0,-1, 1, 0, 1, 0,
     &     1, 0, 1, 0, 1, 0/
*------------------------------------------------------------------------
      ihsec=-ihelem
*------------------------------------------------------------------------
*     Setting the border of the section xx
*------------------------------------------------------------------------
      jsec=indmapOsec(ihsec,1)
      ilocate(1:ndiv)=0
      do idiv=1,ndiv
       ilocate(idiv)=mod(jsec,10)
       jsec=jsec/10
       if(jsec.eq.0)exit
      enddo
      xx(1:6)=xx0(1:6)
      do jdiv=1,idiv
       isec=ilocate(jdiv)
       do isurf=1,6
        is=(isurf+1)/2
        xx(isurf)=xx(isurf)+ishift(isurf,isec)*xoctOsec(is,jdiv,itet)
       enddo
      enddo
*------------------------------------------------------------------------
cFURUTA20160303----------------------------------------------------------
      outside=.false.
      do isurf=1,6
       is=(isurf+1)/2
       t4=x(is)-xx(isurf)
       t5=u(is)
       if(mod(isurf,2).eq.1)then
        t4=-t4
        t5=-t5
       endif
       if(abs(t4).le.coincd*xqt(itet)*abs(t5))t4=t5
       if(t4.gt.0)then
        outside=.true.
        exit
       endif
      enddo
      if(.not.outside)then
cFURUTA20160303----------------------------------------------------------
       nes=indmapOsec(ihsec+1,2)-indmapOsec(ihsec,2)
       do ies=1,nes
        ioutsf=ioutsfmapOsec(indmapOsec(ihsec,2)+ies)
        if(ioutsf.lt.0)exit
        icand=isurf2elem(1,ioutsf2surf(ioutsf))
        call check_inelem(x,u,coincd*xqt(itet),icand,0,ires,lgc)
        if(ires.eq.0)then
         ihelem=icand
         icl=ielem2icl(ihelem)
         exit !FURUTA20180517
        endif
       enddo
      else
       xxx=x(1)
       yyy=x(2)
       zzz=x(3)
       uuu=u(1)
       vvv=u(2)
       www=u(3)
       call tetrafnd(xx0,xxx,yyy,zzz,uuu,vvv,www,coincd,itet,
     &      ihelem,icl,ierr)
      endif
      return
      end subroutine check_sec

************************************************************************
      subroutine tetravol(itet,kvlmax,vol,neleminv)
*
*     Calculate volume of tetrahedrons
*
*     Last Revised 2016/03/09 by T. Furuta
*
************************************************************************
      integer,intent(in) :: itet,kvlmax
      real(8),intent(out) :: vol(kvlmax)
      integer,intent(out) :: neleminv(kvlmax)
      integer :: icl,ielem
      real(8) :: volelem
      do ielem=nelem(itet-1)+1,nelem(itet)
       call calc_volelem(ielem,volelem)
       icl=ielem2icl(ielem)
       neleminv(icl)=neleminv(icl)+1
       vol(icl)=vol(icl)+volelem
      enddo
      return
      end subroutine tetravol

************************************************************************
      subroutine tetratvol(nlat3)
*
*     Calculate volume of all tetrahedrons
*
*     Last Revised 2019/01/10 by T. Furuta
*
************************************************************************
      integer,intent(in) :: nlat3
      integer :: ielem
      real(8) :: volelem
      call allocate_volelems(nelem(nlat3))
      do ielem=1,nelem(nlat3)
       call calc_volelem(ielem,volelem)
       volelems(ielem)=volelem
      enddo
      return
      end subroutine tetratvol

************************************************************************
      subroutine calc_volelem(ielem,vol)
*
*     Calculate volume of an element
*
************************************************************************
      integer,intent(in) :: ielem
      real(8),intent(out) :: vol
      integer :: i,ip(4)
      real(8) :: x(3,4),vec1(3),vec2(3),vec3(3),vec4(3)
      real(8),parameter :: fac=1.0d0/6.0d0
      ip(1:4)=ielem2point(1:4,ielem)
      do i=1,4
       x(1:3,i)=pointxyz(1:3,ip(i))
      enddo
      vec2(1:3)=x(1:3,2)-x(1:3,1)
      vec3(1:3)=x(1:3,3)-x(1:3,1)
      vec4(1:3)=x(1:3,4)-x(1:3,1)
      vec1(1)=vec2(2)*vec3(3)-vec2(3)*vec3(2)
      vec1(2)=vec2(3)*vec3(1)-vec2(1)*vec3(3)
      vec1(3)=vec2(1)*vec3(2)-vec2(2)*vec3(1)
      vol=0.0d0
      do i=1,3
       vol=vol+vec1(i)*vec4(i)
      enddo
      vol=abs(vol)*fac
      return
      end subroutine calc_volelem

************************************************************************
      subroutine set_numblock(nlat3)
************************************************************************
      integer,intent(in) :: nlat3
      integer :: itet
      nEblk(0)=0
      nindEblk(0)=0
      noctEblk(0)=0
      do itet=1,nlat3
       nEblk(itet)=nEblk(itet)+nEblk(itet-1)
       nindEblk(itet)=nindEblk(itet)+nindEblk(itet-1)
       noctEblk(itet)=noctEblk(itet)+noctEblk(itet-1)
      enddo
      return
      end subroutine set_numblock

************************************************************************
      subroutine set_numouter(nlat3)
************************************************************************
      integer,intent(in) :: nlat3
      integer :: itet
      nOsec(0)=0
      nindOsec(0)=0
      noctOsec(0)=0
      do itet=1,nlat3
       nOsec(itet)=nOsec(itet)+nOsec(itet-1)
       nindOsec(itet)=nindOsec(itet)+nindOsec(itet-1)
       noctOsec(itet)=noctOsec(itet)+noctOsec(itet-1)
      enddo
      return
      end subroutine set_numouter

************************************************************************
      subroutine tetralist(kvlmax)
************************************************************************
      integer,intent(in) :: kvlmax
      integer :: i,j,ielem,icl
      integer :: nelemcl(kvlmax)
      logical :: iflag
      integer,allocatable :: ilist(:)
      allocate( ilist(nelemtot) )
      ntetcl=0
      nelemcl(1:kvlmax)=0
      do ielem=1,nelemtot
       icl=ielem2icl(ielem)
       iflag=.true.
       do j=1,ntetcl
        if(icl.eq.ilist(j))then
         nelemcl(j)=nelemcl(j)+1
         iflag=.false.
         exit
        endif
       enddo
       if(iflag)then
        ntetcl=ntetcl+1
        nelemcl(ntetcl)=1
        ilist(ntetcl)=icl
       endif
      enddo
      call allocate_itetlist(nelemtot)
      itetlist(1:ntetcl)=ilist(1:ntetcl)
      itetst(1)=0
      do i=2,ntetcl+1
       itetst(i)=itetst(i-1)+nelemcl(i-1)
      enddo
      nelemcl(1:kvlmax)=0
      do ielem=1,nelemtot
       icl=ielem2icl(ielem)
       do j=1,ntetcl
        if(icl.eq.itetlist(j))then
         nelemcl(j)=nelemcl(j)+1
         exit
        endif
       enddo
       itetelem(itetst(j)+nelemcl(j))=ielem
      enddo
      deallocate( ilist )
      return
      end subroutine tetralist

************************************************************************
      subroutine tetrasorsinit(kvlmax)
************************************************************************
      integer,intent(in) :: kvlmax
      integer :: i,j,ielem
      real(8) :: voltot,volelem
      call allocate_welem(nelemtot)
      do j=1,ntetcl
       voltot=0.0d0
       do i=itetst(j)+1,itetst(j+1)
        ielem=itetelem(i)
        call calc_volelem(ielem,volelem)
        voltot=voltot+volelem
        welem(i)=voltot
       enddo
       welemtot(j)=voltot
      enddo
      return
      end subroutine tetrasorsinit

************************************************************************
      subroutine tetrasors(icl,x,ierr)
************************************************************************
      integer,intent(in) :: icl
      integer,intent(out) :: ierr
      real(8),intent(out) :: x(3)
      integer :: j
      logical :: iflag
      iflag=.false.
      do j=1,ntetcl
       if(icl.eq.itetlist(j))then
        iflag=.true.
        exit
       endif
      enddo
      if(iflag)then
       call fnd_tetsors(j,x,ierr)
      else
       ierr=1
      endif
      return
      end subroutine tetrasors

************************************************************************
      subroutine fnd_tetsors(j,x,ierr)
************************************************************************
      include 'err.inc'

      integer,intent(in) :: j
      integer,intent(out) :: ierr
      real(8),intent(out) :: x(3)
      integer ii,ih,il,i,ist,ielem,numelem
      real(8) :: prob
      real(8) :: dummy,unirn
      external :: unirn
      logical :: iflag
      prob=unirn(dummy)*welemtot(j)
      ist=itetst(j)
      numelem=itetst(j+1)-itetst(j)
      il=ist
      ih=ist+numelem
      iflag=.false.
      do ii=1,999999
       if(ih-il.le.1)then
        i=ih
        iflag=.true.
        exit
       endif
       i=(il+ih)/2
       if(welem(i).ge.prob)then
        ih=i
       else
        il=i
       endif
      enddo
      if(iflag)then
       ielem=itetelem(i)
       call fnd_xinelem(ielem,x,ierr)
      else

       ErrCha = ''
       ErrID = 'L:4675/R:fnd_tetsors/F:tetramod.f' !E80_023_001
       call ErrWrite(ErrID,ErrCha)

       write(*,'(''*** TETRA SOURCE ERROR: '',
     &      ''fail to converge in bisection'')')
       write(*,'(''il ='',i9,'' ih ='',i9,'' i ='',i9)')
     &      ist,ist+numelem,i
       ierr=1
      endif
      return
      end subroutine fnd_tetsors

************************************************************************
      subroutine fnd_xinelem(ielem,x,ierr)
************************************************************************
      include 'err.inc'

      integer,intent(in) :: ielem
      integer,intent(out) :: ierr
      real(8),intent(out) :: x(3)
      integer :: i,ip(4),ii,ires,lgc(4)
      real(8) :: x0(3,4),u(3),vec2(3),vec3(3),vec4(3)
      real(8) :: dummy,unirn,d2,d3,d4
      external :: unirn
      logical :: iflag

      ip(1:4)=ielem2point(1:4,ielem)
      do i=1,4
       x0(1:3,i)=pointxyz(1:3,ip(i))
      enddo
      vec2(1:3)=x0(1:3,2)-x0(1:3,1)
      vec3(1:3)=x0(1:3,3)-x0(1:3,1)
      vec4(1:3)=x0(1:3,4)-x0(1:3,1)
      iflag=.false.
      u(1:3)=0.0d0
      do ii=1,9999
       d2=unirn(dummy)
       d3=unirn(dummy)
       d4=unirn(dummy)
       x(1:3)=x0(1:3,1)+d2*vec2(1:3)+d3*vec3(1:3)+d4*vec4(1:3)
       call check_inelem(x,u,0.0d0,ielem,0,ires,lgc)
       if(ires.eq.0)then
        call tetrasorstrans(ielem,x) !FURUTA20230616
        iflag=.true.
        exit
       endif
      enddo
      if(iflag)then
       ierr=0
      else

       ErrCha = ''
       ErrID = 'L:4727/R:fnd_xinelem/F:tetramod.f' !E80_024_001
       call ErrWrite(ErrID,ErrCha)

       write(*,'(''*** TETRA SOURCE ERROR: '',
     &      ''never found a position inside of tetrahedron'')')
       write(*,'(''ielem ='',i9)')ielem
       ierr=1
      endif
      return
      end subroutine fnd_xinelem

************************************************************************
      subroutine check_geom(nlat3,itgchk,coincd,ierr)
*
*     Check tetrahedron-mesh geometry
*
*     Last modified 2024/06/26 by T. Furuta
*
************************************************************************
      integer,intent(in) :: nlat3,itgchk
      real(8),intent(in) :: coincd
      integer,intent(out) :: ierr
      integer :: itet
      integer,parameter :: iot=68

      write(*,'(''TETRAHEDRAL-MESH GEOM. check has been started '',
     &     '' by itgchk='',i2)')itgchk
      open(iot,file='tet_geoerr.inp',status='unknown')
      close(iot,status='delete')
      do itet=1,nlat3
       call check_volelem(itet,itgchk,ierr)
       if(ierr.ne.0)then
        write(*,'(''*** TETRA GEOMETRY ERROR: '',
     &       ''tetrahedral-mesh geometry itet='',i2,
     &       '' is not compatible'')')itet
        exit
       endif
      enddo
      if(ierr.ne.0)then
       call parastop(300)
      endif
      write(*,'(''--- TETRA VOLUME check has been done ---'')')
      write(*,'(''--- TETRA INTERSECTION check has been started ---'')')
      do itet=1,nlat3
       call check_intersect(itet,itgchk,coincd,ierr)
       if(ierr.ne.0)then
        write(*,'(''*** TETRA GEOMETRY ERROR: '',
     &       ''tetrahedral-mesh geometry itet='',i2,
     &       '' is not compatible'')')itet
        exit
       endif
      enddo
      if(ierr.ne.0)then
       call parastop(300)
      endif
      if(ierr.eq.0)then
       write(*,*)'No error found in TETRAHEDRAL-MESH GEOM. !!'
       write(*,*)'itgchk=0 to skip TETRAHEDRAL-MESH GEOM check'
       call parastop(300)
      endif
      return
      end subroutine check_geom

************************************************************************
      subroutine check_volelem(itet,itgchk,ierr)
*
*     Check intersections among tetrahedrons
*
*     Last Revised 2024/06/26 by T. Furuta
*
************************************************************************
      implicit none
      integer,intent(in) :: itet,itgchk
      integer,intent(out) :: ierr
      integer :: i,nerr
      real(8) :: volelem
      logical :: iflag
      nerr=0
      iflag=.true.
!$OMP PARALLEL
!$OMP& SHARED(itet,itgchk,iflag)
!$OMP& REDUCTION(+:nerr)
!$OMP& PRIVATE(i,volelem)

!$OMP MASTER
!$    write(*,'(''OpenMP PARALLEL PROCESS START'')')
!$OMP END MASTER
!$OMP DO SCHEDULE(guided)
      do i=nelem(itet-1)+1,nelem(itet)
       if(.not.iflag)cycle
       call calc_volelem(i,volelem)
       if(volelem.le.0d0)then
        if(nerr.gt.10)then
         nerr=nerr+1
         write(*,'(''*** TETRA VOLUME ERROR '',
     &         ''more than 10 times. Omitted further report'')')
         iflag=.false.
         cycle
        else
         nerr=nerr+1
         write(*,'(''*** TETRA VOLUME ERROR: '',
     &           ''zero volume tetra found in itet='',i2)')itet
         write(*,'(''ielem ='',i9)')
     &           ielem2id(i)
         if(nerr.eq.itgchk)then
!$OMP CRITICAL (geoerr2)
          call write_geoerr2(i)
!$OMP END CRITICAL (geoerr2)
         endif
        endif
       endif
      enddo
!$OMP END DO
!$OMP MASTER
!$    write(*,'(''OpenMP PARALLEL PROCESS END'')')
!$OMP END MASTER
!$OMP END PARALLEL
      if(nerr.gt.0)then
       ierr=1
      else
       ierr=0
      endif
      return
      end subroutine check_volelem

************************************************************************
      subroutine check_intersect(itet,itgchk,coincd,ierr)
*
*     Check intersections among tetrahedrons
*
*     Last Revised 2024/06/18 by T. Furuta
*
************************************************************************
      integer,intent(in) :: itet,itgchk
      real(8),intent(in) :: coincd
      integer,intent(out) :: ierr
      integer :: ires,lgc(4),ip(4),jp(4),j,ipoint
      integer :: iblk,ielem,jelem,kelem,lelem,indscf,n
      real(8) :: x(3)
      integer :: nerr
      nerr=0
!$OMP PARALLEL
!$OMP& SHARED(itet,itgchk,coincd)
!$OMP& REDUCTION(+:nerr)
!$OMP& PRIVATE(iblk,indscf,ielem,jelem,kelem,lelem)
!$OMP& PRIVATE(ip,jp,j,ipoint,ires,x)

!$OMP MASTER
!$    write(*,'(''OpenMP PARALLEL PROCESS START'')')
!$OMP END MASTER
      do iblk=nindEblk(itet-1)+1,nindEblk(itet)
       if(nerr.gt.11)cycle
       if(indmapEblk(iblk,1).eq.0)cycle
       indscf=indmapEblk(iblk,2)
       n=indmapEblk(iblk+1,2)-indmapEblk(iblk,2)
!$OMP DO SCHEDULE(guided)
       do jelem=1,n
        if(nerr.gt.11)cycle
        ielem=ielemmapEblk(indscf+jelem)
        if(ielem.lt.0)cycle
        ip(1:4)=ielem2point(1:4,ielem)
        do j=1,4
         if(nerr.gt.11)cycle
         ipoint=ip(j)
         x(1:3)=pointxyz(1:3,ipoint)
         do kelem=1,n
          if(nerr.gt.11)cycle
          if(kelem.eq.jelem)cycle
          lelem=ielemmapEblk(indscf+kelem)
          jp(1:4)=ielem2point(1:4,lelem)
          if(ipoint.eq.jp(1).or.ipoint.eq.jp(2)
     &         .or.ipoint.eq.jp(3).or.ipoint.eq.jp(4))cycle
          call check_inelem2(x,coincd*xqt(itet),lelem,ires)
          if(ires.eq.0)then
           if(nerr.gt.10)then
            nerr=nerr+1
            write(*,'(''*** TETRA INTERSECTION ERROR '',
     &           ''more than 10 times. Omitted further report'')')
            exit
           else
            nerr=nerr+1
            write(*,'(''*** TETRA INTERSECTION ERROR: '',
     &           ''intersection found in itet='',i2,'' between'')')itet
            write(*,'(''ielem ='',i9,'' jelem ='',i9)')
     &           ielem2id(ielem),ielem2id(lelem)
            if(nerr.eq.itgchk)then
             call write_geoerr(ielem,lelem)
            endif
           endif
          endif
         enddo
        enddo
       enddo
!$OMP END DO
      enddo
!$OMP MASTER
!$    write(*,'(''OpenMP PARALLEL PROCESS END'')')
!$OMP END MASTER
!$OMP END PARALLEL
      if(nerr.gt.0)then
       ierr=1
      else
       ierr=0
      endif
      return
      end subroutine check_intersect

************************************************************************
      subroutine write_geoerr(ielem,lelem)
*
*     Create input file to show intersection errors by PHIG3D
*
*     Last Revised 2023/02/27 by T. Furuta
*
************************************************************************
      implicit none
      integer,intent(in) :: ielem,lelem
      integer :: is,it,ivoid,iovoid,iorder
      integer :: ip(4),jp(4),isf(4),jsf(4)
      integer,parameter :: iot=68
      real(8) xmin(3),xmax(3),rmin(3),rmax(3)
      character(10) :: chnum,hash1,hash2

      open(iot,file='tet_overlap.inp',status='unknown')
      ip(1:4)=ielem2point(1:4,ielem)
      jp(1:4)=ielem2point(1:4,lelem)
      isf(1:4)=ielem2surf(1:4,ielem)
      jsf(1:4)=ielem2surf(1:4,lelem)
      xmin(1:3)=1.0d55
      xmax(1:3)=-1.0d55
      do is=1,4
       do it=1,3
        xmin(it)=min(xmin(it),pointxyz(it,ip(is)))
        xmin(it)=min(xmin(it),pointxyz(it,jp(is)))
        xmax(it)=max(xmax(it),pointxyz(it,ip(is)))
        xmax(it)=max(xmax(it),pointxyz(it,jp(is)))
       enddo
      enddo
      do it=1,3
       rmin(it)=xmin(it)-0.1d0*(xmax(it)-xmin(it))
       rmax(it)=xmax(it)+0.1d0*(xmax(it)-xmin(it))
      enddo
      if(ielem2id(ielem).ne.100.and.ielem2id(lelem).ne.100)then
       ivoid=100
      else if(ielem2id(ielem).ne.101.and.ielem2id(lelem).ne.101)then
       ivoid=101
      else
       ivoid=102
      endif
      if(ielem2id(ielem).ne.900.and.ielem2id(lelem).ne.900)then
       iovoid=900
      else if(ielem2id(ielem).ne.901.and.ielem2id(lelem).ne.901)then
       iovoid=901
      else
       iovoid=902
      endif
      write(iot,*)'[material]'
      write(iot,*)'mat[1] H 2 O 1 $ Dummy material 1'
      write(iot,*)'mat[2] H 2 O 1 $ Dummy material 2'
      write(iot,*)'[surface]'
      write(iot,5000)ivoid,'rpp',rmin(1),rmax(1)
      write(iot,5001)rmin(2),rmax(2)
      write(iot,5001)rmin(3),rmax(3)
 5000 format(i0,1x,a,1x,2(1x,1pe22.15))
 5001 format(15x,2(1x,1pe22.15))
      do is=1,4
       write(iot,5000)abs(isf(is)),'p',surfcoefs(1:2,abs(isf(is)))
       write(iot,5001)surfcoefs(3:4,abs(isf(is)))
      enddo
      do is=1,4
       write(iot,5000)abs(jsf(is)),'p',surfcoefs(1:2,abs(jsf(is)))
       write(iot,5001)surfcoefs(3:4,abs(jsf(is)))
      enddo
      write(iot,*)'[cell]'
      iorder=aint(log10(real(ielem2id(ielem))))+1
      write(chnum,'(i0)')ielem2id(ielem)
      hash1='#'//chnum
      iorder=aint(log10(real(ielem2id(lelem))))+1
      write(chnum,'(i0)')ielem2id(lelem)
      hash2='#'//chnum
      write(iot,5002)ivoid,0,-ivoid,trim(hash1),trim(hash2)
      write(iot,5002)iovoid,-1,ivoid
      write(iot,5003)ielem2id(ielem),1,-1.0,isf(1:4)
      write(iot,5003)ielem2id(lelem),2,-1.0,jsf(1:4)
 5002 format(i0,1x,i10,1x,i10,2(1x,a))
 5003 format(i0,1x,i10,1x,f5.2,4(1x,i10))
      close(iot)
      end subroutine write_geoerr

************************************************************************
      subroutine write_geoerr2(ielem)
*
*     Create input file to show zero volume tetra error by PHIG3D
*
*     Last Revised 2024/06/26 by T. Furuta
*
************************************************************************
      implicit none
      integer,intent(in) :: ielem
      integer :: is,it,ivoid,iovoid,iorder,ip(4)
      real(8) :: xmin(3),xmax(3),rmin(3),rmax(3)
      real(8) :: radius
      integer,parameter :: iot=68
      character(10) :: chnum,hash(4)

      open(iot,file='tet_geoerr.inp',status='unknown')
      ip(1:4)=ielem2point(1:4,ielem)
      xmin(1:3)=1.0d55
      xmax(1:3)=-1.0d55
      do is=1,4
       do it=1,3
        xmin(it)=min(xmin(it),pointxyz(it,ip(is)))
        xmax(it)=max(xmax(it),pointxyz(it,ip(is)))
       enddo
      enddo
      radius=0.0d0
      do it=1,3
       radius=radius+(xmax(it)-xmin(it))*(xmax(it)-xmin(it))
      enddo
      radius=sqrt(radius)*0.05

      do it=1,3
       rmin(it)=xmin(it)-0.1d0*(xmax(it)-xmin(it))
       rmax(it)=xmax(it)+0.1d0*(xmax(it)-xmin(it))
      enddo
      if(ip(1).ne.100.and.ip(2).ne.100
     &     .and.ip(3).ne.100.and.ip(4).ne.100)then
       ivoid=100
      elseif(ip(1).ne.101.and.ip(2).ne.101
     &      .and.ip(3).ne.101.and.ip(4).ne.101)then
       ivoid=101
      elseif(ip(1).ne.102.and.ip(2).ne.102
     &      .and.ip(3).ne.102.and.ip(4).ne.102)then
       ivoid=102
      elseif(ip(1).ne.103.and.ip(2).ne.103
     &      .and.ip(3).ne.103.and.ip(4).ne.103)then
       ivoid=103
      elseif(ip(1).ne.104.and.ip(2).ne.104
     &      .and.ip(3).ne.104.and.ip(4).ne.104)then
       ivoid=104
      else
       ivoid=105
      endif
      if(ip(1).ne.900.and.ip(2).ne.900
     &     .and.ip(3).ne.900.and.ip(4).ne.900)then
       iovoid=900
      elseif(ip(1).ne.901.and.ip(2).ne.901
     &      .and.ip(3).ne.901.and.ip(4).ne.901)then
       iovoid=901
      elseif(ip(1).ne.902.and.ip(2).ne.902
     &      .and.ip(3).ne.902.and.ip(4).ne.902)then
       iovoid=902
      elseif(ip(1).ne.903.and.ip(2).ne.903
     &      .and.ip(3).ne.903.and.ip(4).ne.903)then
       iovoid=903
      elseif(ip(1).ne.904.and.ip(2).ne.904
     &      .and.ip(3).ne.904.and.ip(4).ne.904)then
       iovoid=904
      else
       iovoid=905
      endif
      write(iot,*)'[material]'
      write(iot,*)'mat[1] H 2 O 1 $ Dummy material 1'
      write(iot,*)'[surface]'
      write(iot,5000)ivoid,'rpp',rmin(1),rmax(1)
      write(iot,5001)rmin(2),rmax(2)
      write(iot,5001)rmin(3),rmax(3)
 5000 format(i0,1x,a,1x,2(1x,1pe22.15))
 5001 format(15x,2(1x,1pe22.15))
      do is=1,4
       write(iot,5000)ip(is),'s',pointxyz(1:2,ip(is))
       write(iot,5001)pointxyz(3,ip(is)),radius
      enddo
      write(iot,*)'[cell]'
      do is=1,4
       iorder=aint(log10(real(ip(is))))+1
       write(chnum,'(i0)')ip(is)
       hash(is)='#'//chnum
      enddo
      write(iot,5002)ivoid,0,-ivoid
      write(iot,5003)trim(hash(1)),trim(hash(2)),
     &     trim(hash(3)),trim(hash(4))
      write(iot,5002)iovoid,-1,ivoid
      do is=1,4
       write(iot,5004)ip(is),1,-1.0,-ip(is)
      enddo
 5002 format(i0,1x,i10,1x,i10)
 5003 format(15x,4(1x,a))
 5004 format(i0,1x,i10,1x,f5.2,4(1x,i10))
      close(iot)
      end subroutine write_geoerr2

************************************************************************
      subroutine check_neighbor(ihelem,n,filename)
*
*     Check neighboring tetrahdrons
*
*     Last Revised 2025/01/31 by T. Furuta
*
************************************************************************
      implicit none
      integer,intent(in) :: ihelem,n
      character(200),intent(in) :: filename
      integer :: iend(n+1)
      integer :: i,j,k,l,nemax,npmax,ipoint,ne,np
      integer :: ie,ielem,isurf,is(4),iel(2),ip(4),iuniv,ind
      integer,allocatable :: ielist(:),ielist0(:),iplist(:),iplist0(:)
      integer,parameter :: nmax=1000
      real(8) :: xx0(6)

      nemax=nmax
      allocate(ielist(nemax))
      ne=1
      iend(1)=1
      ielist(ne)=ihelem
      ie=0
      ielem=ihelem
      do i=1,n
       do
        if(ie.eq.iend(i))then
         iend(i+1)=ne
         exit
        endif
        ie=ie+1
        ielem=ielist(ie)
        is(1:4)=iabs(ielem2surf(1:4,ielem))
        do j=1,4
         iel(1:2)=isurf2elem(1:2,is(j))
         do k=1,2
          if(iel(k).ne.ielem.and.iel(k).ne.0)then
           do l=1,ne
            if(ielist(l).eq.iel(k))exit
           enddo
           if(l.gt.ne)then
            ne=ne+1
            if(ne.gt.nemax)then
             allocate( ielist0(nemax) )
             ielist0(1:nemax)=ielist(1:nemax)
             deallocate( ielist )
             allocate( ielist(nemax+nmax) )
             ielist(1:nemax)=ielist0(1:nemax)
             nemax=nemax+nmax
            endif
            ielist(ne)=iel(k)
           endif
          endif
         enddo
        enddo
       enddo
      enddo

      npmax=nmax
      allocate(iplist(npmax))
      np=0
      do i=1,ne
       ielem=ielist(i)
       ip(1:4)=ielem2point(1:4,ielem)
       do j=1,4
        do k=1,np
         if(ip(j).eq.iplist(k))exit
        enddo
        if(k.gt.np)then
         np=np+1
         if(np.gt.npmax)then
          allocate( iplist0(npmax) )
          iplist0(1:npmax)=iplist(1:npmax)
          deallocate( iplist )
          allocate( iplist(npmax+nmax) )
          iplist(1:npmax)=iplist0(1:npmax)
          npmax=npmax+nmax
         endif
         iplist(np)=ip(j)
        endif
       enddo
      enddo

      ind=len_trim(filename)
      open(500,file=filename(1:ind)//'.ele',status='unknown')
      write(500,'(3i8)')ne,4,1
      do i=1,ne
       ielem=ielist(i)
       ip(1:4)=ielem2point(1:4,ielem)
       if(i.eq.1)then
        iuniv=5000
       else
        iuniv=5001
        is(1:4)=iabs(ielem2surf(1:4,ielem))
        do j=1,4
         if(isurf2elem(2,is(j)).eq.0)iuniv=5002
        enddo
       endif
       write(500,'(6i8)')ielem,ip(1:4),iuniv
      enddo
      close(500)
      open(501,file=filename(1:ind)//'.node',status='unknown')
      write(501,'(4i8)')np,3,0,0
      xx0(1)=1d20
      xx0(2)=-1d20
      xx0(3)=1d20
      xx0(4)=-1d20
      xx0(5)=1d20
      xx0(6)=-1d20
      do i=1,np
       ipoint=iplist(i)
       write(501,'(i8,3f15.8)')ipoint,pointxyz(1:3,ipoint)
       xx0(1)=min(xx0(1),pointxyz(1,ipoint))
       xx0(2)=max(xx0(2),pointxyz(1,ipoint))
       xx0(3)=min(xx0(3),pointxyz(2,ipoint))
       xx0(4)=max(xx0(4),pointxyz(2,ipoint))
       xx0(5)=min(xx0(5),pointxyz(3,ipoint))
       xx0(6)=max(xx0(6),pointxyz(3,ipoint))
      enddo
      write(501,'(a,6f15.8)')'#',xx0(1:6)
      close(501)
      return
      end subroutine check_neighbor

************************************************************************
      subroutine allocate_blocktbl(n1,n2,n3)
************************************************************************
      integer,intent(in) :: n1,n2,n3
      allocate( ielemmapEblk(n1) )
      allocate( indmapEblk(n2,2) )
      allocate( ioctEblk(n3)     )
      return
      end subroutine allocate_blocktbl

************************************************************************
      subroutine allocate_outertbl(n1,n2,n3)
************************************************************************
      integer,intent(in) :: n1,n2,n3
      allocate( ioutsfmapOsec(n1) )
      allocate( indmapOsec(n2,2)  )
      allocate( ioctOsec(n3)      )
      return
      end subroutine allocate_outertbl

************************************************************************
      subroutine deallocate_blocktbl
************************************************************************
      deallocate( ielemmapEblk,indmapEblk,ioctEblk )
      return
      end subroutine deallocate_blocktbl

************************************************************************
      subroutine deallocate_outertbl
************************************************************************
      deallocate( ioutsfmapOsec,indmapOsec,ioctOsec )
      return
      end subroutine deallocate_outertbl

************************************************************************
      subroutine allocate_tetratbl1(nlat3)
************************************************************************
      integer,intent(in) :: nlat3
      integer :: itet
      npoint(0)=0
      nelem(0)=0
      do itet=1,nlat3
       npoint(itet)=npoint(itet)+npoint(itet-1)
       nelem(itet)=nelem(itet)+nelem(itet-1)
      enddo
      allocate( pointxyz(3,npoint(nlat3)) )
      allocate( ielem2point(4,nelem(nlat3)) )
      allocate( ielem2univ(nelem(nlat3)),ielem2icl(nelem(nlat3)) )
      allocate( ielem2surf(4,nelem(nlat3)) )
      allocate( ielem2id(nelem(nlat3)) ) !FURUTA20190110
      return
      end subroutine allocate_tetratbl1

************************************************************************
      subroutine allocate_tetratbl2(nlat3)
************************************************************************
      integer,intent(in) :: nlat3
      integer :: itet
      nsurf(0)=0
      do itet=1,nlat3
       nsurf(itet)=nsurf(itet)+nsurf(itet-1)
      enddo
      noutpt(0)=0               !FURUTA20221104
      noutsf(0)=0               !FURUTA20221104
      allocate( isurf2elem(2,nsurf(nlat3)) )
      allocate( surfcoefs(4,nsurf(nlat3)) )
      allocate( ioutsf2surf(nsurf(nlat3)),ioutpt2point(npoint(nlat3)) )
      allocate( ioutsf2outpt(3,nsurf(nlat3)) )
      return
      end subroutine allocate_tetratbl2

************************************************************************
      subroutine deallocate_tetratbl
************************************************************************
      deallocate( pointxyz,ielem2point,ielem2surf,isurf2elem )
      deallocate( ielem2univ,ielem2icl )
      deallocate( surfcoefs,ioutsf2surf,ioutpt2point,ioutsf2outpt  )
      deallocate( ielem2id ) !FURUTA20190110
      return
      end subroutine deallocate_tetratbl

************************************************************************
      subroutine allocate_tmptbl1(nlat3)
*
************************************************************************
      integer,intent(in) :: nlat3
      allocate( isfind(6,4*nelem(nlat3)) )
      return
      end subroutine allocate_tmptbl1

************************************************************************
      subroutine deallocate_tmptbl1
************************************************************************
      deallocate( isfind )
      return
      end subroutine deallocate_tmptbl1

************************************************************************
      subroutine allocate_tmptbl2(itet,mxsf4p,nside0)
************************************************************************
      integer,intent(in) :: itet,mxsf4p,nside0
      integer :: np0,ns0
      np0=noutpt(itet)-npoint(itet-1)
      ns0=noutsf(itet)-nsurf(itet-1)
      allocate( ioutsf2side(3,ns0) )
      allocate( ioutpt2outsf(mxsf4p+1,np0) )
      allocate( isidemap(nside0+1) )
      allocate( iside2outsf(2*nside0),iside2outpt(2,nside0) )
      return
      end subroutine allocate_tmptbl2

************************************************************************
      subroutine deallocate_tmptbl2
************************************************************************
      deallocate( ioutsf2side,ioutpt2outsf )
      deallocate( isidemap,iside2outsf,iside2outpt)
      return
      end subroutine deallocate_tmptbl2

************************************************************************
      subroutine allocate_itetlist(n)
************************************************************************
      integer,intent(in) :: n
      allocate( itetlist(ntetcl),itetst(ntetcl+1) )
      allocate( itetelem(n))
      return
      end subroutine allocate_itetlist

************************************************************************
      subroutine deallocate_itetlist
************************************************************************
      deallocate( itetlist,itetst,itetelem )
      return
      end subroutine deallocate_itetlist

************************************************************************
      subroutine allocate_welem(n)
************************************************************************
      integer,intent(in) :: n
      allocate( welemtot(ntetcl),welem(n) )
      return
      end subroutine allocate_welem

************************************************************************
      subroutine deallocate_welem
************************************************************************
      deallocate( welemtot,welem )
      return
      end subroutine deallocate_welem

************************************************************************
      subroutine allocate_volelems(n)
************************************************************************
      integer,intent(in) :: n
      allocate( volelems(n) )
      return
      end subroutine allocate_volelems

************************************************************************
      subroutine deallocate_volelems
************************************************************************
      deallocate( volelems )
      return
      end subroutine deallocate_volelems

************************************************************************

      end module TETRAMOD

************************************************************************
      recursive subroutine quicksort(iii,ist,iend)
************************************************************************
      implicit none
      integer,intent(inout) :: iii(6,*)
      integer,intent(in) :: ist,iend
      integer i,j
      integer ix(5),it(5)
      ix(1:5)=iii(1:5,(ist+iend)/2)
      i=ist
      j=iend
      do
       do
        if(iii(1,i).gt.ix(1))then
         exit
        elseif(iii(1,i).eq.ix(1))then
         if(iii(2,i).gt.ix(2))then
          exit
         elseif(iii(2,i).eq.ix(2))then
          if(iii(3,i).ge.ix(3))then
           exit
          endif
         endif
        endif
        i=i+1
       enddo
       do
        if(iii(1,j).lt.ix(1))then
         exit
        elseif(iii(1,j).eq.ix(1))then
         if(iii(2,j).lt.ix(2))then
          exit
         elseif(iii(2,j).eq.ix(2))then
          if(iii(3,j).le.ix(3))then
           exit
          endif
         endif
        endif
        j=j-1
       enddo
       if(i.ge.j)exit
       it(1:5)=iii(1:5,i)
       iii(1:5,i)=iii(1:5,j)
       iii(1:5,j)=it(1:5)
       i=i+1
       j=j-1
      enddo
      if(ist.lt.i-1)call quicksort(iii,ist,i-1)
      if(j+1.lt.iend)call quicksort(iii,j+1,iend)
      end subroutine quicksort
************************************************************************
