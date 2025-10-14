************************************************************************
*                                                                      *
      subroutine bitmap_a_angel(
     &                mmints,fname,
     &                bmpWidth,bmpHeight,
     &                bmpfIType,bmpfIndex,numIType)
*                                                                      *
*        A-Angel.f                                                     *
*        ANGEL : put input file name from PHITS                        *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit none

      include 'param.inc'

      common /imagecom/ imageout ! T.Sato 2017/08/05
      integer :: imageout
*-----------------------------------------------------------------------

      integer :: mmints
      character indat(200)*1
      character fname*100
      integer :: bmpWidth, bmpHeight
      integer :: numIType
      character(1) bmpfIType(numIType)
      integer :: bmpfIndex(numIType)

*-----------------------------------------------------------------------

      integer :: icang
      integer :: i

*-----------------------------------------------------------------------

      if(imageout.eq.1) return ! need not to draw bmp file

      do i = 1, 200

         indat(i) = ' '

      end do

      do i = 1, 100

         indat(i) = fname(i:i)

      end do

         icang = 1

      call bitmap_a_main0(
     &        icang, indat, fname,
     &        bmpWidth, bmpHeight,
     &        bmpfIType, bmpfIndex, numIType)

*-----------------------------------------------------------------------

      end

