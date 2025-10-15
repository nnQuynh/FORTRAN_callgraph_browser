c ------------------------------------------------------------
C contains
C  ErrWrite
C  ErrWriteIO
c ------------------------------------------------------------
      subroutine ErrWrite(ErrID,ErrCha)
      implicit real*8 (a-h,o-z)
      common /paraj/ mstz(300), parz(300)
      character(len=*),intent(in) :: ErrCha
      character(len=*),intent(in) :: ErrID

      if(ErrCha /= "" ) write(6,120) trim(ErrCha)
      if(mstz(147).ne.0) then
       if(ErrID /= "") write(6,100) trim(ErrID)
      endif

100   format('  ---',a)
120   format(' ',a)

      return
      end subroutine

c ------------------------------------------------------------
      subroutine ErrWriteIO(ErrID,ErrCha,io)
      implicit real*8 (a-h,o-z)
      common /paraj/ mstz(300), parz(300)
      character(len=*),intent(in) :: ErrCha
      character(len=*),intent(in) :: ErrID
      integer :: io

      if(ErrCha /= "" ) write(io,120) trim(ErrCha)
      if(mstz(147).ne.0) then
       if(ErrID /= "") write(io,100) trim(ErrID)
      endif

100   format('  ---',a)
120   format(' ',a)

      return
      end subroutine

