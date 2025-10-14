c ------------------------------------------------------------
c Subroutine ErrorCheck
C contains
C  dchain_outtime_check
C  dmax1Check
C  MPI_dumpFileCheck
C  SourceCheck
C  ZspcReplace
C  ZspcLog

c ------------------------------------------------------------
c  D-Chain Outtime Error Check
c ------------------------------------------------------------
      subroutine dchain_outtime_check(tbina_tmp,techa,step)
      implicit none
      include 'err.inc'

      character(len=*) :: techa
      integer :: step
      integer :: DCHAIN_Eflag

      data before_time / 0.0 /
      data now_time / 0.0 /

      DCHAIN_Eflag = 0

      if(techa == 'y' ) now_time = tbina_tmp * 31557600.0d0
      if(techa == 'd' ) now_time = tbina_tmp * 86400.0d0
      if(techa == 'h' ) now_time = tbina_tmp * 3600.0d0
      if(techa == 'm' ) now_time = tbina_tmp * 60.0d0
      if(techa == 's' ) now_time = tbina_tmp

      if( step >= 2 ) then
       if( abs(now_time) <= abs(before_time) .and.
     &  now_time*before_time > 0) then ! the same sign, introduced by T.Sato 2021/03/11
        DCHAIN_Eflag = 1
       endif
      endif

      before_time = now_time

      if(DCHAIN_Eflag == 1) then
       write(ErrCha,'(''DCHAIN OUTTIME STEP Warning'',
     & '': outtime step'',i3,'' =< outtime step '',i3)') step,step-1

       ErrID = 'L:46/R:dchain_outtime_check/F:errcheck.f' !W64_001_001
       call ErrWrite(ErrID,ErrCha)

       write(ErrCha,'(''It may cause problem in DCHAIN-SP'')')
       ErrID = 'L:50/R:dchain_outtime_check/F:errcheck.f'
       call ErrWrite(ErrID,ErrCha)

      endif

      return
      end subroutine
c ------------------------------------------------------------
c dmax(1) check
c ------------------------------------------------------------
       subroutine dmax1Check(ht)
       use CHARVARMOD, only:close_rwtfile
       implicit real*8 (a-h,o-z)
       include 'err.inc'

       character(len=*) ht
       integer :: suffix
       common /ndemax/ dnmax(20)
       common /eparm/  esmax, esmin, emin(20)

       suffix = index(ht(10:10),'h')
       if( suffix /= 0 ) then
        if( dnmax(1) .gt. emin(1) ) then

         write(ErrCha,'(''Cannot find proton nuclear data library. '',
     &   ''Please revise xsdir or delete dmax(1) from your input'')')
         ErrID = 'L:76/R:dmax1Check/F:errcheck.f' !E02_003_001
         call ErrWrite(ErrID,ErrCha)
         call close_rwtfile(1)
         call parastop(999)
        endif
       endif

       return
       end subroutine

c ------------------------------------------------------------
c MPI Dump File Check
c ------------------------------------------------------------
       subroutine MPI_dumpFileCheck(ndumpmax)
       implicit real*8 (a-h,o-z)

       include 'err.inc'
       integer,intent(in) :: ndumpmax
       common /mpi00/ npe, me

       write(ErrCha,*) '** MPI-mode Error **'
       ErrID = 'L:97/R:MPI_dumpFileCheck/F:errcheck.f' !E84_005_001
       call ErrWrite(ErrID,ErrCha)

       write(ErrCha,*) 'The number of dump files is ',ndumpmax
       ErrID = 'L:101/R:MPI_dumpFileCheck/F:errcheck.f'
       call ErrWrite(ErrID,ErrCha)

       write(ErrCha,*) 'MPI Thread(npe-1) set is    ',npe-1
       ErrID = 'L:105/R:MPI_dumpFileCheck/F:errcheck.f'
       call ErrWrite(ErrID,ErrCha)

       write(6,*)'You need to prepare sufficient # of dump files.'

       write(6,*)'You may need to re-distribute dump files. ',
     &      'Ref. Sec. 8 of the PHITS manual.'


       return
       end subroutine

c ------------------------------------------------------------
c [source] section Check
c ------------------------------------------------------------
      subroutine SourceCheck(ierr)
      implicit none

      include 'err.inc'
      integer,intent(out) :: ierr
      ierr=0
      if(SOURCE_Eflag == 0) then
       write(ErrCha,*) '[source] section is not specified.'
       ErrID = 'L:128/R:SourceCheck/F:errcheck.f' !E03_110_001
       call ErrWrite(ErrID,ErrCha)
       ierr=1
      endif

      return
      end subroutine


c ------------------------------------------------------------
c Zenkaku Space Replace
c ------------------------------------------------------------
      subroutine ZspcReplace(chin)
      implicit none

      include 'err.inc'
      character(len=*) :: chin
      integer :: zspc_flag

      zspc_flag = 0
      do
       column = index(chin,'@')
       if( column /= 0 ) then
        chin(column:column+1) = '  '
        zspc_flag = 1
       else
        exit
       endif
      enddo

      if(zspc_flag == 1)  then
        zspc_count = zspc_count + 1
        ZenkakuSpaceLine(zspc_count) = chin
      endif

      return
      end subroutine
c ------------------------------------------------------------
c Zenkaku Space Replace Log
c ------------------------------------------------------------
      subroutine ZspcLog
      implicit none

      include 'err.inc'
      integer :: i

       i=1
       do
        if(ZenkakuSpaceLine(i) == '') exit
         write(ErrCha,'(''Warning:Include Zenkaku Space:'',a200)')
     & ZenkakuSpaceLine(i)
        ErrID = 'L:179/R:ZspcLog/F:errcheck.f' !W00_009_001
        call ErrWrite(ErrID,ErrCha)
        i=i+1
       enddo

      return
      end subroutine
