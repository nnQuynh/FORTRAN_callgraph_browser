************************************************************************
*                                                                      *
      subroutine a_angel(mmints,fname)
*                                                                      *
*        A-Angel.f                                                     *
*        ANGEL : put input file name from PHITS                        *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

      include 'param.inc'

      common /imagecom/ imageout ! T.Sato 2017/08/05
*-----------------------------------------------------------------------

      character indat(200)*1
      character fname*100

*-----------------------------------------------------------------------

      if(imageout.eq.1) return ! need not to draw eps file

      do i = 1, 200

         indat(i) = ' '

      end do

      do i = 1, 100

         indat(i) = fname(i:i)

      end do

         icang = 1

      call a_main0(icang,indat)

*-----------------------------------------------------------------------

      end


