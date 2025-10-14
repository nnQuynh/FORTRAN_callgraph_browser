************************************************************************
*                                                                      *
      subroutine writeompfile(m,dmpd)
*                                                                      *
************************************************************************
      implicit none
      include 'param.inc'
*-----------------------------------------------------------------------
      integer itmdp,itmdf
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)
*-----------------------------------------------------------------------
      integer idmpomp
      common /idmpomp0/idmpomp
      integer ioomp
      common /ioomp0/ioomp
!$OMP THREADPRIVATE(/ioomp0/)
*-----------------------------------------------------------------------
      integer,intent(in) :: m
      real(8),intent(in) :: dmpd(30)
      integer k

      if(idmpomp.gt.0)then
       write(ioomp)m
       write(ioomp)
     &      (dmpd(itmdp(m,k)),k = 1,abs(itmdp(m,0)))
      else
       write(ioomp,'(i5)')m
       write(ioomp,'(30(1p1d24.15))')
     &      (dmpd(itmdp(m,k)),k = 1,abs(itmdp(m,0)))
      endif
      return
      end

************************************************************************
*                                                                      *
      subroutine openompfile
*                                                                      *
************************************************************************
      implicit none
*-----------------------------------------------------------------------
      integer npe,me
      integer ipomp,npomp
      common /mpi00/ npe, me
      common /ipomp0/ipomp,npomp
!$OMP THREADPRIVATE(/ipomp0/)
*-----------------------------------------------------------------------
      integer idmpomp
      common /idmpomp0/idmpomp
      integer ioomp
      common /ioomp0/ioomp
!$OMP THREADPRIVATE(/ioomp0/)
*-----------------------------------------------------------------------
      integer io2,iord,iord2
      character(100) filnm
      character(8) chipomp,chme
*-----------------------------------------------------------------------

      ioomp=800+ipomp
      iord=aint(log10(real(ipomp+1))) + 1
      write(chipomp,'(i8.8)') ipomp+1
      if(npe.le.1)then
       filnm='ompdmptmp-'//chipomp(9-iord:8)
      else
       iord2=aint(log10(real(me))) + 1
       write(chme,'(i8.8)') me
       filnm='ompdmptmp-'//chme(9-iord2:8)//'-'//chipomp(9-iord:8)
      endif
!$OMP CRITICAL (ompfile_crit) ! To avoid Intel Compiler bug 20220523
      if(idmpomp.gt.0)then
       open(ioomp,file=filnm,form='unformatted',status='unknown')
       write(ioomp)'START'
      else
       open(ioomp,file=filnm,form='formatted',status='unknown')
       write(ioomp,'(a)')'START'
      endif
!$OMP END CRITICAL (ompfile_crit)

      return
      end

************************************************************************
*                                                                      *
      subroutine cpompfile
*                                                                      *
************************************************************************
      implicit none
      include 'param.inc'
*-----------------------------------------------------------------------
      integer itmdp,itmdf
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)
*-----------------------------------------------------------------------
      integer idmpomp
      common /idmpomp0/idmpomp
      integer ioomp
      common /ioomp0/ioomp
!$OMP THREADPRIVATE(/ioomp0/)
*-----------------------------------------------------------------------
      integer i,k,m,ios
      real(8) dmpd(30)
      character(5)dummy

      rewind(ioomp)
      if(idmpomp.gt.0)then
       read(ioomp)dummy
       do
        read(ioomp,iostat=ios)m
        if(ios.lt.0)exit
        read(ioomp)(dmpd(itmdp(m,k)),k=1,abs(itmdp(m,0)))
        write(itmdf(m))(dmpd(itmdp(m,k)),k=1,abs(itmdp(m,0)))
       enddo
      else
       read(ioomp,'(a)')dummy
       do
        read(ioomp,'(i5)',iostat=ios)m
        if(ios.lt.0)exit
        read(ioomp,'(30(1p1d24.15))')
     &       (dmpd(itmdp(m,k)),k=1,abs(itmdp(m,0)))
        write(itmdf(m),'(30(1p1d24.15))')
     &       (dmpd(itmdp(m,k)),k=1,abs(itmdp(m,0)))
       enddo
      endif
      close(ioomp,status='DELETE')
      return
      end

************************************************************************
*                                                                      *
      subroutine chkdmpomp(io,m,idump,ierr)
*                                                                      *
************************************************************************
!$    use omp_lib
      use t4dtrack_mod, only: it4dtrack
      implicit none

      include 'err.inc'
*-----------------------------------------------------------------------
      integer idumpall
      common /dumpall/ idumpall
*-----------------------------------------------------------------------
      integer idmpomp
      common /idmpomp0/idmpomp
      data idmpomp /0/
*-----------------------------------------------------------------------
      integer,intent(in) :: io,m,idump
      integer,intent(out) :: ierr

      ierr=0
      if(idmpomp.eq.0)then
!$OMP PARALLEL
!$OMP SINGLE
!$    idmpomp=OMP_GET_NUM_THREADS()
!$OMP END SINGLE
!$OMP BARRIER
!$OMP END PARALLEL
       if(idmpomp.eq.1)idmpomp=0
       if(idmpomp.ne.0)then
        idmpomp=idump
        if(idumpall.ne.0)then
         write(ErrCha,'('' **** Error in '',i2,
     &              ''-th tally : tally dump option is not''
     &              '' compatible with idumpall option for OpenMP'')') m

         ErrID = 'L:163/R:chkdmpomp/F:ompdump.f' !E00_009_001
         call ErrWrite(ErrID,ErrCha)
         call ErrWriteIO(ErrID,ErrCha,io)
         ierr=1
        endif
        if(it4dtrack.gt.0)then
         write(ErrCha,'('' **** Error in '',i2,
     &              ''-th tally : tally dump option is not''
     &              '' compatible with T-4Dtrack for OpenMP'')') m

         ErrID = 'L:173/R:chkdmpomp/F:ompdump.f' !E00_009_001
         call ErrWrite(ErrID,ErrCha)
         call ErrWriteIO(ErrID,ErrCha,io)
         ierr=1
        endif
       endif
      else
       if(idmpomp*idump.lt.0)then
        write(ErrCha,'('' **** Error in '',i2,
     &                  ''-th tally : dump format should coinside''
     &                  '' with other tallies for OpenMP'')') m
         ErrID = 'L:184/R:chkdmpomp/F:ompdump.f' !E00_009_002
         call ErrWrite(ErrID,ErrCha)
         call ErrWriteIO(ErrID,ErrCha,io)

        ierr=1
       endif
      endif

      return
      end
