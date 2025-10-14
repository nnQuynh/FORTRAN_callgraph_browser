C =====================================================
C  modul dedx_file
C  [ material ] "dedxfile" card
C =====================================================
      module dedx_file
      implicit none

      character(len=60),allocatable :: dedx_filename(:)

      type dedx_value
       integer :: ecounts
       real(8),allocatable :: Ebin(:),dedxbin(:)
      end type

      type dedx_info
       integer :: kfcounts
       integer,allocatable :: kfcode(:)
       real(8),allocatable :: dedx_Efactor(:,:)
       type(dedx_value),allocatable :: dedxvalue(:)
      end type

      type(dedx_info),allocatable :: fdedx(:)

      integer :: idedx(20)

      integer :: unitdedx
      real(8) :: ufactor
      character(len=9) :: kfcode_char

      contains
C =====================================================
C  dedx_file_allocate : call from read00.f [setpar]
C =====================================================
      subroutine dedx_file_allocate(mxmat)
      implicit none

      integer :: mxmat

      if( mxmat /= 0 ) then

       allocate( dedx_filename(mxmat) )
       allocate( fdedx(mxmat) )

       dedx_filename(:) = ' '
       fdedx(:)%kfcounts = 0

      endif

      end subroutine dedx_file_allocate
C =====================================================
C  dedx_file_deallocate : call from main.f
C =====================================================
      subroutine dedx_file_deallocate(mxmat)
      implicit none

      integer :: mxmat

       if( mxmat /= 0 ) then

        deallocate( dedx_filename )
        deallocate( fdedx )

       endif

      end subroutine dedx_file_deallocate

C =====================================================
      subroutine get_ecounts(mat,kfno,j)
      implicit none

      integer :: mat,kfno,j

      j = fdedx(mat)%dedxvalue(kfno)%ecounts

      end subroutine get_ecounts

C =====================================================
      subroutine get_kfcounts(mat,dedx_kfcode,kfno)
      implicit none

      integer :: mat,dedx_kfcode,kfno
      integer :: i

      do i=1,fdedx(mat)%kfcounts
       if(dedx_kfcode == fdedx(mat)%kfcode(i)) then
        kfno = i
        exit
       endif
      enddo

      end subroutine get_kfcounts

C =====================================================
      subroutine get_dedx_kfcode(dedxmat,kfno,zp,mp)
C ref. utl02.f : function kfft
      implicit none

C memo.     implicit double precision (a-h,o-z)
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)

C ref. ovly15.f subroutine cpcasdump
      integer :: no,mat,ityp,ktyp,jtyp,mtyp,ctyp
      real(8) :: rtyp

      integer :: dedxmat,kfno,dedx_kfcode
      integer :: zp
      real(8) :: mp

      if(zp .eq. 0) zp = ctyp
      if(zp .eq. 0) zp = max(nint(mp/2),1)

      dedx_kfcode = 0
C proton
      if( ityp == 1 ) then
       dedx_kfcode = 2212
C pion+
      elseif( ityp == 3 ) then
       dedx_kfcode = 211
C pion-
      elseif( ityp == 5 ) then
       dedx_kfcode = -211
C muon+
      elseif( ityp == 6 ) then
       dedx_kfcode = -13
C muon-
      elseif( ityp == 7 ) then
       dedx_kfcode =  13
C kaon+
      elseif( ityp == 8 ) then
       dedx_kfcode =  321
C kaon-
      elseif( ityp == 10 ) then
       dedx_kfcode = -321
C other
      elseif( ityp  >=  15 ) then
       dedx_kfcode = zp * 1000000 + nint(mp)
      endif

      call get_kfcounts(dedxmat,dedx_kfcode,kfno)

      end subroutine get_dedx_kfcode

C =======================================================
C call from ggm01.f : setmd (read dedxfile information )
C =======================================================

      subroutine dedx_file_read(mat)
      use NGSDATAMOD, only : bindeg
      implicit none

      include 'err.inc'
      include 'param-physcnst.inc'

      integer :: i,j
      integer :: mat
      character buf*200,buf2*200
      integer :: i1,i2,i3,i4

      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn, ilfn
      character mltfl*200
      integer iotmp
      logical exex
      integer ipprt,ipneu

      character cmmt*2 ! T.Sato 2019/01/13

      data cmmt/'#$'/ ! comment out character

c ---------------------------------------------------------------

      mltfl = chfn(29)(1:ilfn(29))//'/'
     &//dedx_filename(mat)(1:len_trim(dedx_filename(mat)))

      inquire( file = mltfl, exist = exex )

       if( exex .eqv. .false. ) then

       ErrID = 'L:181/R:dedx_file_read/F:dedx_file.f' !E04_027_001
       ErrCha = ''
       call ErrWrite(ErrID,ErrCha)

        write(6,'(/'' Error : input data file for'',
     &            '' dedx does not exist.''/
     &            '' file name = '',200a1)')
     &            ( mltfl(i:i),i=1, ilfn(29)+30 )
       stop
       end if

c
C read00.f memo
C *-----------------------------------------------------------------------
C *    open temporary file io = 29
C *-----------------------------------------------------------------------



C main.f memo
*           io = 30-39   are reserved for input files

      iotmp = 30
      open(iotmp, file = mltfl, status = 'old')

      call dedx_file_check(iotmp,mat)

      do i=1,fdedx(mat)%kfcounts
       call get_ecounts(mat, i , j )
       allocate(fdedx(mat)%dedxvalue(i)%Ebin(j))
       allocate(fdedx(mat)%dedxvalue(i)%dedxbin(j))

       fdedx(mat)%dedxvalue(i)%Ebin(:) = 0.0d0
       fdedx(mat)%dedxvalue(i)%dedxbin(:) = 0.0d0
      enddo

      rewind(iotmp)

C  Energy , Stopping Power read
      j = 0
      do
       do
        read(iotmp,101,end=99) buf
        call chlngt(buf,200,i1,i2)  ! calculate character length
        call chcaps(buf,i1,i2,i3,cmmt) ! CAPITAL to lower character
        i1=1                         ! delete space even for 1st character
        call chcomp(buf,i1,i2,i4)   ! delete space
        if(buf(1:3).eq.'kf=') then
         j=j+1 ! number of kf code
         exit
        endif
       enddo
       do i=1,fdedx(mat)%dedxvalue(j)%ecounts
  10    read(iotmp,101) buf
        buf2=buf
        call chlngt(buf,200,i1,i2)  ! calculate character length
        call chcaps(buf,i1,i2,i3,cmmt) ! CAPITAL to lower character
        i1=1                         ! delete space even for 1st character
        call chcomp(buf,i1,i2,i4)   ! delete space
        if(i4.eq.0) goto 10 ! only comment, skip the line
        read(buf2,*)
     &    fdedx(mat)%dedxvalue(j)%Ebin(i),
     &    fdedx(mat)%dedxvalue(j)%dedxbin(i)

       enddo

C unitdedx 1 : MeV -> MeV/n
c unitdedx 2 : MeV/u -> MeV/n
c unitdedx 3 : MeV/n -> MeV/n as it is
       if( unitdedx .le. 2) then ! MeV or MeV/n
C proton
        if(fdedx(mat)%kfcode(j) == 2212) then
         ufactor = 1.0d0  ! number of nucleon
         if(unitdedx.eq.2) ufactor = ufactor / (rstms(1)/rstms(0))
c pion+
        elseif(fdedx(mat)%kfcode(j) == 211) then
         ufactor = 1.0d0  ! number of nucleon
         if(unitdedx.eq.2) ufactor = ufactor / (rstms(3)/rstms(0))
c pion-
        elseif(fdedx(mat)%kfcode(j) == -211) then
         ufactor = 1.0d0  ! number of nucleon
         if(unitdedx.eq.2) ufactor = ufactor / (rstms(5)/rstms(0))
c muon+
        elseif(fdedx(mat)%kfcode(j) == -13) then
         ufactor = 1.0d0  ! number of nucleon
         if(unitdedx.eq.2) ufactor = ufactor / (rstms(6)/rstms(0))
c muon-
        elseif(fdedx(mat)%kfcode(j) == 13) then
         ufactor = 1.0d0  ! number of nucleon
         if(unitdedx.eq.2) ufactor = ufactor / (rstms(7)/rstms(0))
c kaon+
        elseif(fdedx(mat)%kfcode(j) == 321) then
         ufactor = 1.0d0  ! number of nucleon
         if(unitdedx.eq.2) ufactor = ufactor / (rstms(8)/rstms(0))
c kaon-
        elseif(fdedx(mat)%kfcode(j) == -321) then
         ufactor = 1.0d0  ! number of nucleon
         if(unitdedx.eq.2) ufactor = ufactor / (rstms(10)/rstms(0))
C deutron
        elseif(fdedx(mat)%kfcode(j) == 1000002) then
         ufactor = 2.0d0  ! number of nucleon
         if(unitdedx.eq.2) ufactor = ufactor / (rstms(15)/rstms(0))
C triton
        elseif(fdedx(mat)%kfcode(j) == 1000003) then
         ufactor = 3.0d0  ! number of nucleon
         if(unitdedx.eq.2) ufactor = ufactor / (rstms(16)/rstms(0))
C 3He
        elseif(fdedx(mat)%kfcode(j) == 2000003) then
         ufactor = 3.0d0  ! number of nucleon
         if(unitdedx.eq.2) ufactor = ufactor / (rstms(17)/rstms(0))
C alpha
        elseif(fdedx(mat)%kfcode(j) == 2000004) then
         ufactor = 4.0d0  ! number of nucleon
         if(unitdedx.eq.2) ufactor = ufactor / (rstms(18)/rstms(0))
c heavy ions
        else
         ipprt=int(fdedx(mat)%kfcode(j)/1000000)
         if(ipprt.le.0) then ! not heavy ion
          ErrID = 'L:299/R:dedx_file_read/F:dedx_file.f' !E04_030_001
          write(ErrCha,*) 'invalid kf code',fdedx(mat)%kfcode(j),
     &    ' in dedxfile= ',( mltfl(i:i),i=1, ilfn(29)+30 )
          call ErrWrite(ErrID,ErrCha)
          stop
         endif
         ipneu=fdedx(mat)%kfcode(j)-ipprt*1000000 ! mass number
         ipneu=ipneu-ipprt
         ufactor = (ipprt+ipneu)*1.0d0
         if(unitdedx.eq.2) ufactor = ufactor / ((rstms(1)*ipprt +
     &   rstms(2)*ipneu - bindeg(ipprt,ipneu)/1000.0d0) / rstms(0))
        endif


        do i=1,fdedx(mat)%dedxvalue(j)%ecounts
         fdedx(mat)%dedxvalue(j)%Ebin(i) =
     &   fdedx(mat)%dedxvalue(j)%Ebin(i) / ufactor
        enddo
       endif
      enddo

   99 continue
      close(iotmp)

  101 format(a200)
      return
      end subroutine dedx_file_read

C =====================================================
C    dedx_file_check
C =====================================================

      subroutine dedx_file_check(iotmp,mat)
      implicit none

      include 'err.inc'

      integer :: j,jj
      character buf*200,buf2*200
      character cmmt*2 ! T.Sato 2019/01/13

      integer :: iotmp
      integer :: mat
      integer :: i1,i2,i3,i4
      real(8) :: Ebin_tmp(2)

C =====================================================
      data cmmt/'#$'/ ! comment out character

C  Default Unit : MeV
      unitdedx = 1

      j = 0

c  unit and number of kf code read
      do
       read(iotmp,101,end=99) buf
       call chlngt(buf,200,i1,i2)  ! calculate character length
       call chcaps(buf,i1,i2,i3,cmmt) ! CAPITAL to lower character
       i1=1                         ! delete space even for 1st character
       call chcomp(buf,i1,i2,i4)   ! delete space
       if(buf(1:5).eq.'unit=') then
        read(buf(6:200),*) unitdedx
       elseif(buf(1:3).eq.'kf=') then
        j=j+1 ! number of kf code
       endif
      enddo
   99 continue

      fdedx(mat)%kfcounts = j

      if(j /= 0) then
       allocate( fdedx(mat)%kfcode(j) )
       allocate( fdedx(mat)%dedxvalue(j) )
       allocate( fdedx(mat)%dedx_Efactor(j,2) )

       fdedx(mat)%kfcode(:) = 0
       fdedx(mat)%dedxvalue(:)%ecounts = 0
       fdedx(mat)%dedx_Efactor(:,:) = 0.0d0

      endif

      rewind(iotmp)

c  dedx file Array counts & sorting check

      Ebin_tmp = 0.0d0
      j = 0
      jj = 0
      do
       read(iotmp,101,end=999) buf
       buf2=buf ! keep original line
       call chlngt(buf,200,i1,i2)  ! calculate character length
       call chcaps(buf,i1,i2,i3,cmmt) ! CAPITAL to lower character
       i1=1                         ! delete space even for 1st character
       call chcomp(buf,i1,i2,i4)   ! delete space
       if(buf(1:5).eq.'unit=') then
        read(buf(6:200),*) unitdedx
       elseif(buf(1:3).eq.'kf=') then
        if(j.ne.0.and.jj.eq.0) then
         write(ErrCha,*) 'No dEdx data'
         ErrID = 'L:400/R:dedx_file_check/F:dedx_file.f' !E04_027_003
         call ErrWrite(ErrID,ErrCha)
         write(6,*) '  filename:', dedx_filename(mat)
         write(6,*) '  kfcode  :', fdedx(mat)%kfcode(j)
         stop
        endif
        jj=0
        j=j+1 ! number of kf code
        read(buf(4:200),*) fdedx(mat)%kfcode(j)
       elseif(i4.eq.0) then
        continue ! only comment, so skip the line
       else
        jj = jj + 1
        fdedx(mat)%dedxvalue(j)%ecounts
     &  = fdedx(mat)%dedxvalue(j)%ecounts + 1
        if(jj >= 2) Ebin_tmp(2) = Ebin_tmp(1)
        read(buf2,*) Ebin_tmp(1)
        if(jj >= 2) then
         if( Ebin_tmp(2) >= Ebin_tmp(1) ) then
          write(ErrCha,*) 'dEdXFile Energy Sorting Error :'
          ErrID = 'L:420/R:dedx_file_check/F:dedx_file.f' !E04_027_002
          call ErrWrite(ErrID,ErrCha)
          write(6,*) '  filename:', dedx_filename(mat)
          write(6,*) '  kfcode  :', fdedx(mat)%kfcode(j)
          write(6,*) '  data no.:', jj
          stop
         endif
        endif
       endif
      enddo

  999 continue

  101 format(a200)
      return
      end subroutine dedx_file_check

C ==========================================================
C *   log-log interpolation  [Ref. talls01.f]
C ==========================================================

      subroutine dedx_file_interpolation(i,mat,ene,dedx)
      implicit none

      integer :: i,mat
      real(8) :: ene,dedx
      integer :: j

      j = 0
      call dedx_file_binary_search(i,mat,ene,j)

       dedx  = fdedx(mat)%dedxvalue(i)%dedxbin(j)
     & * exp( log( fdedx(mat)%dedxvalue(i)%dedxbin(j-1)
     & / fdedx(mat)%dedxvalue(i)%dedxbin(j))
     & * log( ene / fdedx(mat)%dedxvalue(i)%Ebin(j))
     & / log( fdedx(mat)%dedxvalue(i)%Ebin(j-1)
     & / fdedx(mat)%dedxvalue(i)%Ebin(j)))

      return
      end subroutine dedx_file_interpolation

C ==========================================================
C *   binary search
C ==========================================================

      subroutine dedx_file_binary_search(i,mat,ene,midid)
      implicit none

      integer :: i,mat
      real(8) :: ene
      integer :: lowid,midid,highid
      integer :: j

      lowid  = 0
      midid  = 0
      highid = 0

      highid = fdedx(mat)%dedxvalue(i)%ecounts

      j=0
      do
       j = j + 1
       midid = int( ( lowid + highid ) / 2 )
       if(( fdedx(mat)%dedxvalue(i)%Ebin(midid) >= ene ).and.
     &    ( fdedx(mat)%dedxvalue(i)%Ebin(midid-1) <= ene )) then
         exit
       elseif( fdedx(mat)%dedxvalue(i)%Ebin(midid) < ene ) then
        lowid  = midid + 1
       else
        highid = midid - 1
       endif
      enddo

      return
      end subroutine dedx_file_binary_search

C ==========================================================
      subroutine dedx_file_select(m,kfno,energy,dedxver)
      implicit none

      include 'atimasys.inc'

      double precision m, sum
      integer i,j

      double precision dedx

      real(8) :: energy
      double precision dedxver
      integer dedx_kf,kfno

C- factor cal ---------------------------------
      if(( fdedx(dedxmat)%dedx_Efactor(kfno,1) == 0.0d0 ).and.
     &   ( fdedx(dedxmat)%dedx_Efactor(kfno,2) == 0.0d0 )) then

C- min factor ---------
      sum = 0.d0
      dedxver = 0.d0

      do i = 1, nnuc
      sum = sum+anuc(i) * mt(i) *
     & dedx(zp, zt(i), mp, mt(i), pot(i), rho, fntp,
     & fdedx(dedxmat)%dedxvalue(kfno)%Ebin(1),
     & gas) / m
      enddo

      dedxver = sum

      fdedx(dedxmat)%dedx_Efactor(kfno,1)
     & =(fdedx(dedxmat)%dedxvalue(kfno)%dedxbin(1)/1000.d0)/dedxver

C- max factor ---------
      j = 0
      call get_ecounts(dedxmat, kfno , j )

      sum = 0.d0
      dedxver = 0.d0

      do i = 1, nnuc
      sum = sum+anuc(i) * mt(i) *
     & dedx(zp, zt(i), mp, mt(i), pot(i), rho, fntp,
     & fdedx(dedxmat)%dedxvalue(kfno)%Ebin(j),
     & gas)/m
      enddo

      dedxver = sum

      fdedx(dedxmat)%dedx_Efactor(kfno,2)
     & =(fdedx(dedxmat)%dedxvalue(kfno)%dedxbin(j)/1000.d0)/dedxver


      endif

C- end factor cal ---------------------------------

c--- file out of range --------------------------
       j = 0
       call get_ecounts(dedxmat, kfno , j )
C- min
       if(fdedx(dedxmat)%dedxvalue(kfno)%Ebin(1) >= energy) then

       sum = 0.d0
       dedxver = 0.d0

       do i = 1, nnuc
       sum = sum+anuc(i) * mt(i) *
     &    dedx(zp, zt(i), mp, mt(i), pot(i), rho, fntp, energy, gas) / m
       enddo

       dedxver = sum * fdedx(dedxmat)%dedx_Efactor(kfno,1)

C- max
       elseif( fdedx(dedxmat)%dedxvalue(kfno)%Ebin(j) <= energy) then

       sum = 0.d0
       dedxver = 0.d0

       do i = 1, nnuc
       sum = sum+anuc(i) * mt(i) *
     &    dedx(zp, zt(i), mp, mt(i), pot(i), rho, fntp, energy, gas) / m
       enddo

       dedxver = sum * fdedx(dedxmat)%dedx_Efactor(kfno,2)

C- file range -
       elseif
     & ((fdedx(dedxmat)%dedxvalue(kfno)%Ebin(1)< energy).and.
     &  (fdedx(dedxmat)%dedxvalue(kfno)%Ebin(j)> energy))then

C- log - log interpolation
       dedxver = 0.0d0
       call dedx_file_interpolation(kfno,dedxmat,energy,dedxver)
       dedxver = dedxver / 1000.0d0

       endif

      return
      end subroutine dedx_file_select
C
      end module dedx_file

************************************************************************
*                                                                      *
      subroutine char_int(m,h,ix)
*                                                                      *
*       m=1 encode the char h into ix.                                 *
*       m=2 decode into h the char that is coded in ix.                *
*       Last modified by NAIS on 2019/11/01                            *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      dimension ix(20)
      character h*60

*-----------------------------------------------------------------------

         if(m.eq.2) goto 200

         do 100 i = 1, 20
           is=(i-1)*3 + 1
           ix(i) = ichar(h(is:is))+256*(ichar(h(is+1:is+1))
     &           + 256*ichar(h(is+2:is+2)))
  100    continue
         return

  200    continue
         do 220 i = 1, 20
         is=(i-1)*3 + 1
         h(is:is+2) = char(mod(ix(i),256))//char(mod(ix(i)/256,256))
     &             // char(ix(i)/65536)
  220    continue
*-----------------------------------------------------------------------

      return
      end
