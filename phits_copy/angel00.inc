************************************************************************
*                                                                      *
*     angel00.inc                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
!     when you increase mdas over 250000000, you have to define integer*8 instead of integer in the next line
      integer mdas,mmdas,mmmax
      real*8  das

      parameter ( mdas = 50000000 )

      common /mdasa/ das( mdas )
      common /mdasb/ mmmax


*-----------------------------------------------------------------------

