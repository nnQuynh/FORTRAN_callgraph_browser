************************************************************************
      subroutine read_tetra(ierr)
*
*     Read tetra files
*
*     Created by T. Furuta on 2019/02/08
*
************************************************************************
      use TETRAMOD, only: tetraread,itetra,
     &     write_infofile2,read_infofile2,read_infofile3,
     &     lbin,binfile,binfile2
      implicit real*8 (a-h,o-z)
      include 'param.inc'
*-----------------------------------------------------------------------
      integer nlat3,itetcl(10),itetvol,itetauto,itgchk,i
      common /tetf1/ nlat3,itetcl,itetvol,itetauto,itgchk
      integer ltfile(10),itfform(10)
      real(8) tetsfac(10)
      character(200) tfilename(10)
      common /tetf2/ tetsfac,ltfile,itfform,tfilename
      common /paraj/  mstz(300), parz(300)
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
*-----------------------------------------------------------------------
      integer,intent(out) :: ierr
      integer :: nuniv(0:10),nmats(0:10)
      integer :: iorguniv(kvlmax),nunivtot,igcelst
      integer,allocatable :: iunivmat(:),materials(:)
      real(8),allocatable :: densuniv(:)
*-----------------------------------------------------------------------
      integer iot1
      iot1 = 69
*-----------------------------------------------------------------------
*     Defined here because SETPAR is called after
      itetra  = mstz(120)       !FURUTA20190331
      itetauto = mstz(133)
      lbin=ilfn(30)
      binfile  = chfn(30)
      do i=lbin,1,-1
       if(chfn(30)(i:i).eq.'.')exit
      enddo
      if(i.eq.1)i=lbin+1
      binfile2 = chfn(30)(1:i-1)//'2'//chfn(30)(i:199)
      if(mstz(158).gt.0.and.itetra.ne.0)then !itgchk>0
       write(*,'("*** TETRA WARNING: ",
     &           "itgchk>0 TETRAHEDRAL-MESH GEOM. check works ",
     &           "only with itetra=0")')
       write(*,'(" itetra=0 is selected")')
       itetra=0
      endif
*-----------------------------------------------------------------------
      ierr=0
      call tetraread(tetsfac,nlat3,ltfile,itfform,tfilename,ierr)
      if(ierr.ne.0)return
      if(itetauto.eq.1)then
*-----------------------------------------------------------------------
       if(abs(itetra).ne.1)then
        call count_tetrauniv(nlat3,nuniv,iorguniv)
        nunivtot=nuniv(nlat3)
        allocate( iunivmat(nunivtot),materials(nunivtot) )
        allocate( densuniv(nunivtot) )
        call count_tetramat(nlat3,nunivtot,ltfile,itfform,tfilename,
     &       nuniv,iorguniv,nmats,iunivmat,materials,densuniv,ierr)
        call add_tetraunivsurf
        call add_tetraunivcell(nunivtot,iunivmat,densuniv,igcelst)
        call check_tetraunivmat(nlat3,nunivtot,nmats,materials,ierr)
       else
*-----------------------------------------------------------------------
        call read_infofile2(iot1,nlat3,nunivtot,nmats)
        allocate( iunivmat(nunivtot),materials(nunivtot) )
        allocate( densuniv(nunivtot) )
        call read_infofile3(iot1,nlat3,nunivtot,nmats,
     &       iunivmat,materials,densuniv)
        call add_tetraunivsurf
        call add_tetraunivcell(nunivtot,iunivmat,densuniv,igcelst)
        call check_tetraunivmat(nlat3,nunivtot,nmats,materials,ierr)
       endif
*-----------------------------------------------------------------------
       if(abs(itetra).eq.2)then
        call write_infofile2(iot1,nlat3,nunivtot,nmats,
     &       iunivmat,materials,densuniv)
       endif
*-----------------------------------------------------------------------
       deallocate( iunivmat,materials,densuniv )
      endif
      return
      end subroutine read_tetra
************************************************************************
      subroutine count_tetrauniv(nlat3,nuniv,iorguniv)
*
*     Count universe for tetrahedrons
*
*     Created by T. Furuta on 2019/02/08
*
************************************************************************
      use TETRAMOD, only:nelem,ielem2univ
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      integer,intent(in) :: nlat3
      integer,intent(out) :: nuniv(0:10),iorguniv(kvlmax)
      integer itet,ielem,iuniv,nelemtot,numuniv
      integer :: iflag(kvlmax)

      nelemtot=nelem(nlat3)
      nuniv(0:10)=0
      iorguniv(1:kvlmax)=0
      numuniv=0
      do itet=1,nlat3
       iflag(1:kvlmax)=0
       do ielem=nelem(itet-1)+1,nelem(itet)
        iuniv=ielem2univ(ielem)
        if(iflag(iuniv).eq.0)then
         numuniv=numuniv+1
         iflag(iuniv)=numuniv
         iorguniv(numuniv)=iuniv
        endif
        ielem2univ(ielem)=5000+iflag(iuniv)
       enddo
       nuniv(itet)=numuniv
      enddo
      return
      end subroutine count_tetrauniv

************************************************************************
      subroutine count_tetramat(nlat3,nunivtot,ltfile,itfform,tfilename,
     &     nuniv,iorguniv,nmats,iunivmat,materials,densuniv,ierr)
*
*     Add huge sphere surface for tetra universe
*
*     Created by T. Furuta on 2019/02/08
*
************************************************************************
      use TETRAMOD, only: read_tetratxt0,read_tetratxt1,read_tetratxt2
      implicit real*8 (a-h,o-z)
      include 'err.inc'
      include 'param.inc'
      integer,intent(in) :: nlat3,nunivtot,ltfile(10),itfform(10)
      character(200),intent(in) :: tfilename(10)
      integer,intent(in) :: nuniv(0:10),iorguniv(kvlmax)
      integer,intent(out) :: nmats(0:10)
      integer,intent(out) :: iunivmat(nunivtot),materials(nunivtot)
      real(8),intent(out) :: densuniv(nunivtot)
      integer,intent(out) :: ierr
      integer :: iot,itet
      character(200) :: txtfile

      iot=68
*-----------------------------------------------------------------------
      nmats(0:10)=0
      do itet=1,nlat3
       if(itfform(itet).eq.0)then
        write(txtfile,'(200x)')
        txtfile(1:ltfile(itet)+4)
     &       =tfilename(itet)(1:ltfile(itet))//'.txt'
        open(iot,file=txtfile,status='old')
        call read_tetratxt0(itet,iot,kvlmax,nunivtot,
     &       nuniv,iorguniv,nmats,iunivmat,materials,densuniv,ierr)
        close(iot)
        if(ierr.ne.0)then
         write(ErrCha,'(''*** TETRA ERROR: in '',a200)')
     &        txtfile
         ErrID = 'L:162/R:count_tetramat/F:ggs04.f' !E80_006_001
         call ErrWrite(ErrID,ErrCha)
         call parastop(999)
        endif
       elseif(itfform(itet).eq.1)then
        write(txtfile,'(200x)')
        txtfile(1:ltfile(itet))=tfilename(itet)(1:ltfile(itet))
        open(iot,file=txtfile,status='old')
        call read_tetratxt1(itet,iot,kvlmax,nunivtot,
     &       nuniv,iorguniv,nmats,iunivmat,materials,densuniv,ierr)
        close(iot)
        if(ierr.ne.0)then
         write(ErrCha,'(''*** TETRA ERROR: in '',a200)')
     &        txtfile
         ErrID = 'L:176/R:count_tetramat/F:ggs04.f' !E80_006_002
         call ErrWrite(ErrID,ErrCha)
         call parastop(999)
        endif
       elseif(itfform(itet).eq.2)then
        write(txtfile,'(200x)')
        txtfile(1:ltfile(itet))=tfilename(itet)(1:ltfile(itet))
        call read_tetratxt2(itet,txtfile,nunivtot,nuniv,
     &       nmats,iunivmat,materials,densuniv,ierr)
        if(ierr.ne.0)then
         write(ErrCha,'(''*** TETRA ERROR: in '',a200)')txitfile
         ErrID = 'L:187/R:count_tetramat/F:ggs04.f' !E80_006_002
         call ErrWrite(ErrID,ErrCha)
         call parastop(999)
        endif
       endif
      enddo
      return
      end subroutine count_tetramat
************************************************************************
      subroutine add_tetraunivsurf
*
*     Add huge sphere surface for tetra universe
*
*     Created by T. Furuta on 2019/02/08
*
************************************************************************
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      common /inggs/  iog, igcel, ioa, igsuf, iob, igtrs
      common /celdb/  idsn(kvlmax), idtn(kvlmax)
*-----------------------------------------------------------------------
      igsuf=igsuf+1
      idsn(igsuf)=5000
      write(ioa) 0, 0, 5, 1, 1.0d30
      return
      end subroutine add_tetraunivsurf
************************************************************************
      subroutine add_tetraunivcell(nunivtot,iunivmat,densuniv,igcelst)
*
*     Add universes required by tetra
*
*     Created by T. Furuta on 2019/02/08
*
************************************************************************
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'
*-----------------------------------------------------------------------
      character m_err*200
      common /error/ m_err, l_err, k_err
      character dkam*6
*-----------------------------------------------------------------------
      common /regdm/  idmg(kvlmax)
      common /regda/  ichl(kvlmax), chsm(kvlmax), ichmx, iod
      character       chsm*10
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /regdu/  iuni(kvlmax)
      common /regde/  ichp(kvlmax), ilat(kvlmax), idct(kvlmax)
      common /inggs/  iog, igcel, ioa, igsuf, iob, igtrs
      common /celda/  deng(kvlmax)
      character(9),parameter :: cpre='-5000  u='
      character(4) cnum
      character chrg
      equivalence ( das, chrg )
*-----------------------------------------------------------------------
      integer,intent(in) :: nunivtot,iunivmat(nunivtot)
      real(8),intent(in) :: densuniv(nunivtot)
      integer,intent(out) :: igcelst
      integer :: iuniv
*-----------------------------------------------------------------------
      igc=5
      igcelst=igcel
      do iuniv=1,nunivtot
       igcel=igcel+1
       if(igcel.gt.kvlmax) goto 992
       idrg(igcel)=5000+iuniv
       if( idgr(idrg(igcel)).ne.0 ) goto 993
       idgr(idrg(igcel))=igcel
       idmg(igcel)=5000+iunivmat(iuniv)
       iuni(igcel)=5000+iuniv
       deng(igcel)=densuniv(iuniv)
       write(cnum,'(i4)')5000+iuniv
       ichl(igcel) = igc
       ichp(igcel) = igc+9
       ichmx = max(ichmx,ichp(igcel))
       write(iod) cpre(1:9)//cnum(1:4)//' '
      enddo

      return
*-----------------------------------------------------------------------
 992  continue

         write(dkam,'(i6)') kvlmax
         m_err = '# of cells exceeds kvlmax = '//dkam//
     &        '. Increase kvlmax.'
         ErrCha = ''
         ErrID = 'L:273/R:add_tetraunivcell/F:ggs04.f'
         l_err = 0
         k_err = 0
         ierr  = 1
         return
*-----------------------------------------------------------------------
 993  continue

         write(dkam,'(i6)') idrg(igcel)
         m_err = 'ID number of cell is duplicated. = '// dkam
         ErrCha = ''
         ErrID = 'L:284/R:add_tetraunivcell/F:ggs04.f'
         l_err = 0
         k_err = 0
         ierr  = 1
         return
*-----------------------------------------------------------------------
      end subroutine add_tetraunivcell
************************************************************************
      subroutine check_tetraunivmat(nlat3,nunivtot,nmats,materials,ierr)
*
*     Check materials required by tetra
*
*     Created by T. Furuta on 2019/02/08
*
************************************************************************
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /regdm/  idmg(kvlmax)
*-----------------------------------------------------------------------
      integer,intent(in) :: nlat3,nunivtot,nmats(0:10)
      integer,intent(in) :: materials(nunivtot)
      integer,intent(inout) :: ierr
*-----------------------------------------------------------------------
      integer :: itet,i
      do itet=1,nlat3
       do i=nmats(itet-1)+1,nmats(itet)
        if(idnm(5000+i).eq.0)then
         write(6,'(''*** ERROR : undefined material'')')
         write(6,'('' TETRA material (MID): '',i5'' for itet='',i2)')
     &        materials(i),itet
         write(6,'(''   should be defined as '',
     &             ''material ID number: '',i5)') 5000+i
         ierr = ierr + 1
        end if
       enddo
      enddo
      return
      end subroutine check_tetraunivmat

************************************************************************
      subroutine settetra(io,ierr)
*
*     Set cells composed of tetrahedrons
*
*     Created by T. Furuta on 2015/07/14
*
************************************************************************
      use TETRAMOD, only: tetrainit,tetralist
      use moddas_ggs !frtati20220905
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
*-----------------------------------------------------------------------
      integer,intent(in) :: io
      integer,intent(out) :: ierr
      integer itet
      real(8) xx0(6)
      character(200) nodefile,elemfile,volfile
      integer nlat3,itetcl(10),itetvol,itetauto,itgchk
      common /tetf1/ nlat3,itetcl,itetvol,itetauto,itgchk
      integer ltfile(10),itfform(10)
      real(8) tetsfac(10)
      character(200) tfilename(10)
      common /tetf2/ tetsfac,ltfile,itfform,tfilename
*-----------------------------------------------------------------------
      call tetrafingshow !FURUTA20221209 itetragshow=0
      call tetrainit(io,nlat3,itgchk,coincd,ierr)
      if(ierr.ne.0)return
      call tetrauniv(nlat3,ierr)
      if(ierr.ne.0)return
      call tetralist(kvlmax)
      if(itetvol.eq.1)call settvol(nlat3)
      return
      end subroutine settetra
************************************************************************
      subroutine tetrauniv(nlat3,ierr)
*
*     Set universe and cell number for each tetrahedron
*
*     Created by T. Furuta on 2015/07/14
*
************************************************************************
      use TETRAMOD, only:nelem,ielem2univ,ielem2icl
      use moddas_ggs !frtati20220905
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      include 'err.inc'
      integer,intent(out) :: ierr
      integer ic,ielem,nerr
      integer :: iflag(kvlmax)
      logical :: jflag

      jflag=.false.
      nerr=0
      iflag(1:kvlmax)=0
      ielem2icl(1:nelem(nlat3))=0
      do ielem=1,nelem(nlat3)
       do ic=1,mxa
        if(jun(ic).ne.ielem2univ(ielem))cycle
        if(ielem2icl(ielem).ne.0)then
         iflag(ic)=ielem2univ(ielem)
         jflag=.true.
        endif
        ielem2icl(ielem)=ic
        lat(1,ic)=-3
       enddo
       if(ielem2icl(ielem).eq.0)then
        do ic=mxa+1,mxa+nerr
         if(iflag(ic).ne.ielem2univ(ielem))cycle
         ielem2icl(ielem)=ic
        enddo
        if(ielem2icl(ielem).eq.0)then
         nerr=nerr+1
         iflag(mxa+nerr)=ielem2univ(ielem)
         jflag=.true.
        endif
       endif
      enddo
      if(jflag)then
       do ic=1,mxa
        if(iflag(ic).ne.0)then
         write(ErrCha,*)'*** TETRA ERROR: element univese contains',
     &        ' more than one cell'
         ErrID = 'L:409/R:tetrauniv/F:ggs04.f' !E80_001_001
         call ErrWrite(ErrID,ErrCha)

         write(*,'(''univ= '',i5)')iflag(ic)
         ierr=1
        endif
       enddo
       if(nerr.gt.0)then
        write(ErrCha,*)'*** TETRA ERROR: universe required by element',
     &       ' is missing'
        ErrID = 'L:419/R:tetrauniv/F:ggs04.f' !E80_001_002
        call ErrWrite(ErrID,ErrCha)

        do ic=mxa+1,mxa+nerr
         write(*,'(''missing univ= '',i5)')iflag(ic)
        enddo
        ierr=1
       endif
      endif
      return
      end subroutine tetrauniv

************************************************************************
      subroutine settvol(nlat3)
*
*     Set volume of each cell composed of tetrahedrons
*
*     Created by T. Furuta on 2015/07/14
*     Last Revised 2016/03/09 by T. Furuta
*
************************************************************************
      use TETRAMOD, only: tetravol
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      integer,intent(in) :: nlat3
      integer :: ic,iuniv,itet
      integer :: neleminv(kvlmax)
      real(8) :: vol(kvlmax)
*-----------------------------------------------------------------------
      common /volreg/ dvol(kvlmax)
*-----------------------------------------------------------------------
      neleminv(1:kvlmax)=0
      vol(1:kvlmax)=0.0d0
      do itet=1,nlat3
       call tetravol(itet,kvlmax,vol,neleminv)
      enddo
      do ic=1,mxa
       if(neleminv(ic).gt.0)then
        dvol(ic)=vol(ic)
       endif
      enddo
      return
      end subroutine settvol

************************************************************************
      subroutine tetrabox(io,k,xx0,ierr)
*
*     Set borders of box for tetrahedron mesh
*
*     Created by T. Furuta on 2015/07/14
*
************************************************************************
      use moddas_ggs !frtati20220905
      implicit real*8 (a-h,o-z)
*-----------------------------------------------------------------------
      include 'param.inc'
      include 'ggsparam.inc'
*-----------------------------------------------------------------------
      integer nlat3,itetcl(10),itetvol,itetauto,itgchk
      common /tetf1/ nlat3,itetcl,itetvol,itetauto,itgchk
*-----------------------------------------------------------------------
      integer,intent(in) :: io,k
      integer,intent(inout) :: ierr
      real(8),intent(out) :: xx0(6)
      integer n,j,js,ks,ls,ic

      ic = itetcl(k-10000)
      n = abs(lca(ic+1))-lca(ic)

      if(n.ne.6)then
       write(io,'(/''** ERROR : in lattice setup,''/
     &''    only RPP definition is allowed with LAT=3 for cell ='',i7)')
     &      ncl(ic)
       ierr=1
       return
      endif
      do j = 1, 6
       js = abs(lja(lca(ic)+j-1))
       ks = abs(kst(js))
       ls = lsc(js)
       select case(ks)
       case(1)
        xx0(j)=xx0(j-1)*1.00001d0 ! Tiny shift to avoid coincident surface
        xx0(j-1)=-scf(ls+4)
       case(2:4)
        xx0(j)=scf(ls+1)
       case default
        write(io,'(/''** ERROR : in lattice setup,''/
     &    ''    only RPP definition is allowed for LAT=3'')')
        ierr=1
        return
       end select
      enddo
      return
      end subroutine tetrabox

************************************************************************
      subroutine tetrasorstrans(ielem,x)
*
*     Transform source position in tetra to main space
*
*     Created by T. Furuta on 2023/06/16
*
************************************************************************
      use moddas_ggs !frtati20220905
      use TETRAMOD, only: nelem
      implicit real*8 (a-h,o-z)
*-----------------------------------------------------------------------
      include 'ggsparam.inc'
*-----------------------------------------------------------------------
      integer,intent(in) :: ielem
      real(8),intent(inout) :: x(3)
*-----------------------------------------------------------------------
      integer nlat3,itetcl(10),itetvol,itetauto,itgchk
      common /tetf1/ nlat3,itetcl,itetvol,itetauto,itgchk
*-----------------------------------------------------------------------
      integer :: itet,icl,ilev,iu,itrns
      real(8) :: x0(3)

      do itet=1,nlat3
       if(ielem.le.nelem(itet))exit
      enddo
      icl=itetcl(itet)          ! Cell ID for tetrabox
      do ilev=1,mxlv
       iu=jun(icl)
       if(iu.eq.0)exit
       do icl=1,mxa
        if(mfl(1,icl).eq.iu)exit
       enddo
       if(mfl(3,icl).gt.0)then
        x0(1:3)=x(1:3)
        itrns=mfl(3,icl)
        call trnsxv(x0(1),x0(2),x0(3),x(1),x(2),x(3),itrns)
       endif
      enddo
      return
      end subroutine tetrasorstrans

************************************************************************
      subroutine tetra0iii
*
*     Set parameter iii to be zero for gshow
*
*     Created by T. Furuta on 2016/07/01
*
************************************************************************
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      iii=0
      return
      end subroutine tetra0iii

************************************************************************
      subroutine tetranowgshow
*
*     Set parameter itetragshow to be one for gshow & 3dshow
*
*     Created by T. Furuta on 2022/12/09
*
************************************************************************
      use TETRAMOD, only: itetragshow
      implicit none
      itetragshow=1
      return
      end subroutine tetranowgshow

************************************************************************
      subroutine tetrafingshow
*
*     Set parameter itetragshow to be zero for transport
*
*     Created by T. Furuta on 2022/12/09
*
************************************************************************
      use TETRAMOD, only: itetragshow
      implicit none
      itetragshow=0
      return
      end subroutine tetrafingshow

************************************************************************
      subroutine tetragetmat(matnum)
*
*     Get mat from cell number icl for lost particle
*     without updating mat
*
*     Created by T. Furuta on 2020/06/12
*
************************************************************************
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      integer,intent(out) :: matnum
      common /regdm/  idmg(kvlmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      integer nmed
      nmed=idmg(icl)
      if(nmed.gt.0)then
       nmed=idnm(nmed)
      endif
      matnum=nmed !FURUTA20200612
      return
      end subroutine tetragetmat

************************************************************************
      subroutine tetrasetmat(matnum)
*
*     Updating mat for lost particle
*
*     Last modified by T. Furuta on 2020/06/12
*
************************************************************************
      implicit real*8 (a-h,o-z)
      integer,intent(in) :: matnum
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      mat=matnum
      return
      end subroutine tetrasetmat

************************************************************************
      subroutine tetrawarning(ihelem0,ihelem)
*
*     Show waring in region check of tetrahedrons
*
*     Created by T. Furuta on 2015/07/14
*
************************************************************************
      use MMBANKMOD
      implicit real*8 (a-h,o-z)
      integer,intent(in) :: ihelem0,ihelem
      include 'param.inc'
      include 'ggsparam.inc'
*-----------------------------------------------------------------------
      common /mpi00/ npe, me
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /cgerr/  nlost, ilost, igerr, icger, ncger, nrecover
*-----------------------------------------------------------------------

      icger=icger+1
      if( icger .le. nrecover ) then

       write(6,'(/''*** warning in tetra check *** no ='',i6)') icger

       if( npe .gt. 0 )
     &      write(6,'( '' my ip       = '',i3)') me

       write(6,'( '' nbch ncs no = '',3i10)')
     &      nobch, nocas, no
       write(6,'( '' ityp, e(no) = '',i4,1x,e17.8)')
     &      ityp, e(ibke+no,ipomp+1)
       write(6,'( '' tetra ini fin     ='',5i6)')
     &      ihelem0, ihelem
       write(6,'( '' x, y, z :'',3e17.8)')
     &      x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),z(ibkz+no,ipomp+1)
       write(6,'( '' u, v, w :'',3e17.8)')
     &      u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),w(ibkw+no,ipomp+1)

       write(6,'(''*** succeeded in recovering *** tetra check'')')

      endif

      return
      end subroutine tetrawarning

************************************************************************
      subroutine openfoam_create_filename(
     &                itfoam,fname,
     &                foamfIType, foamfIndex, numIType,
     &                fileIndex, sfname)
*
*       create OpenFoam filename. Copy from bitmap_create_filename.
*
*       last modified by T. Furuta on 2019/10/28
*
************************************************************************
      implicit none
      integer :: itfoam
      character(1) :: fname(100)
      character(1) :: foamfIType(*)
      integer :: foamfIndex(*)
      integer :: numIType
      integer :: fileIndex
      character(len=*) :: sfname
*-----------------------------------------------------------------------
      integer :: isize, ipos, maxdig, newisize, index
      integer :: i
      character(len=20) :: sgid
      character(len=10) :: format
      integer :: tid(numIType), gid(numIType)
      integer :: get_line_length
      integer :: last_char_index
      integer :: max_digits_integers
*-----------------------------------------------------------------------

      call clear_string(sfname)
      isize = get_line_length(fname, 100)
      call char2string(isize, fname, sfname)
      ipos = last_char_index(sfname, '.')
      isize = ipos-1

      maxdig = max_digits_integers(foamfIndex, numIType)
      write(format,'(''(I'',1I0,''.'',1I0,'')'')') maxdig, maxdig

      index = fileIndex
      do i = numIType, 2, -1
       tid(i) = product(foamfIndex(1:i-1))
       gid(i) = floor((index-1) / tid(i) * 1.0d0) + 1
       index = index - (gid(i)-1)*tid(i)
      end do
      gid(1) = index

      do i = numIType, 1, -1    ! from outer do-loop to inner
       if ( foamfIndex(i) .gt. 1 ) then
        write(sgid,format) gid(i)
        newisize = isize + 2 + len_trim(sgid)
        sfname(isize+1:newisize) = '_'//foamfIType(i)//trim(sgid)
        isize = newisize
       end if
      end do

      if(itfoam.eq.1)then
       sfname(isize+1:isize+5) = '.foam'
       isize = isize + 5
      elseif(itfoam.eq.2)then
       sfname(isize+1:isize+4) = '.csv'
       isize = isize + 4
      endif
      end subroutine openfoam_create_filename

************************************************************************
      subroutine tetraHDF5setmat(imat,nmatst,nelems,
     &     cisonames,disofracs,cmatname)
*
*     Material setting for HDF5 file
*
*     Last modified by T. Furuta on 2024/11/14
*
************************************************************************
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      integer,intent(in) :: imat,nmatst,nelems
      character(200),intent(in) :: cisonames(nelems)
      real(8),intent(in) :: disofracs(nelems)
      character(80),intent(in) :: cmatname
*-----------------------------------------------------------------------
      common /kmat1a/ mxmat, mxmat0, mxnel
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /kmat1o/ iom1, iom2, iom3
      common /mtnmc/  smtnc(kvlmax), dmtnc(kvlmax,2),
     &                mtncn, mtnc(kvlmax,2), nmtnc(kvlmax,2)
      character dmtnc*80
*-----------------------------------------------------------------------
      integer :: nel,libh,igas,istp,inlb,iplb,ielb,icnd,iulb,ihlb,imts
      integer :: idedx(20)
      real(8) :: denh
      character(200) :: chin,chlw
      integer :: i,j,ic,i3,icha,masi,libi,ierr,ii
      character(1) :: c1

!------------------------------------------------------------------------
!     Addition to [Material]
!------------------------------------------------------------------------
      nel=0
      igas=0
      istp=0
      inlb=0
      iplb=0
      ielb=0
      icnd=0
      iulb=0
      ihlb=0
      imts=0
      idedx(1:20)=0
      denh=0.0d0
      libh = ichar(' ') + 256 * ( ichar(' ') + 256 * ichar(' ') )
      open(iom2,form='unformatted',status='scratch')
      mxmat=mxmat+1
      idmat=5000+nmatst+imat
      idnm(idmat) = mxmat
      idmn(mxmat) = idmat
      do i=1,nelems
       ic=1
       i3=len_trim(cisonames(i))
       chin=cisonames(i)
       chlw=cisonames(i)
       do j=ic,i3
        c1=chlw(j:j)
        if(c1.ge.'A'.and.c1.le.'Z')
     &       c1=char(ichar(c1)+ichar('a')-ichar('A'))
        chlw(j:j)=c1
       enddo
       call readnc(chin,chlw,ic,i3,icha,masi,libi,ierr)
       if(icha.eq.1.and.masi.eq.1)then
        denh=denh+disofracs(i)
        libh=libi
       else
        nel = nel + 1
        write(iom2)icha,masi,disofracs(i),libi
       endif
      enddo
      write(iom1) nel
      write(iom1) denh, libh
      write(iom1) igas, istp, inlb, iplb, ielb, icnd
      write(iom1) iulb, ihlb
      write(iom1) imts
      write(iom1) (idedx(ii),ii=1,20)
      if(nel.gt.0)rewind(iom2)
      do i=1,nel
       read(iom2) icha,masi,denst,libi
       write(iom1)icha,masi,denst,libi
      enddo
      close(iom2)
!------------------------------------------------------------------------
!     Addition to [MatNameColor]
!------------------------------------------------------------------------
      if(mtncn+1.gt.kvlmax)return
      mtncn=mtncn+1
      mtnc(mtncn,1)=idmat
      mtnc(mtncn,2)=idmat
      nmtnc(mtncn,1)=0
      nmtnc(mtncn,2)=0
      smtnc(mtncn)=1.0d0
      i3=len_trim(cmatname)
      dmtnc(mtncn,1)(1:i3)=cmatname(1:i3)
      nmtnc(mtncn,1)=i3
!------------------------------------------------------------------------
      return
      end subroutine tetraHDF5setmat

