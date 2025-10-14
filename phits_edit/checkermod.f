      module CHECKERMOD
c      implicit none
      integer inocas00
      save inocas00
      data inocas00/0/
!$OMP THREADPRIVATE(inocas00)

      contains

      subroutine CHECKERPARC(CHAR)
      include 'jam1.inc'
      include 'jam2.inc'
      character(5) CHAR

      write(100+inocas00,*)CHAR,':parc',parc(5)

      end subroutine CHECKERPARC

      end module CHECKERMOD
