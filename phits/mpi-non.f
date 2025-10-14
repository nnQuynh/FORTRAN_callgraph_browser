************************************************************************
*                                                                      *
*     This file includes dummy subroutines and function                *
*     for non parallel version                                         *
*                                                                      *
************************************************************************
*                                                                      *
      subroutine paraint(mn)
*                                                                      *
************************************************************************

      return
      end

************************************************************************
*                                                                      *
      function paratim()
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

            paratim = 1.0

      return
      end


************************************************************************
*                                                                      *
      subroutine parafin
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      return
      end


************************************************************************
*                                                                      *
      subroutine paraiset(icc)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)




      return
      end

************************************************************************
*                                                                      *
      subroutine paratal(tr,br,nr)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      return
      end


************************************************************************
*                                                                      *
      subroutine paratal_sumover(br,m)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      return
      end


************************************************************************
*                                                                      *
      subroutine parabcsti(icc)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      return
      end


************************************************************************
*                                                                      *
      subroutine parasi(icc,nc,ip)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      return
      end

************************************************************************
*                                                                      *
      subroutine parari(icc,nc,ip)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      return
      end

************************************************************************
*                                                                      *
      subroutine parasr(rcc,nc,ip)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      return
      end

************************************************************************
*                                                                      *
      subroutine pararr(rcc,nc,ip)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      return
      end

************************************************************************
*                                                                      *
      subroutine paraiccp
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      return
      end

************************************************************************
*                                                                      *
      subroutine parastop(idnum)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)
      include 'err.inc'

! T.Sato 2015/03/09, to avoid forced stop due to parastop
! T.Ogawa 2018/3/14, Add idnum = 111
      common /errorcom/ncascerr
!$OMP THREADPRIVATE(/errorcom/)

      integer*8 :: iransb64 ! S.H. xorshift (2020.2.6)
      common /randtp/ rijk,rans,ranb, iransb64
!$OMP THREADPRIVATE(/randtp/)

      if(idnum.eq.111.or.idnum.eq.115.or.idnum.eq.116.or.idnum.eq.119
     &  .or.idnum.eq.120)then ! error from JAM
       if(ncascerr.eq.0) then ! first time called this routine
        write(ErrCha,*) 'parastop is called by idnum =',idnum
        ErrID = 'L:168/R:parastop/F:mpi-non.f' !E84_004_003
        call ErrWrite(ErrID,ErrCha)

       Select Case(idnum)
        case(111)
        write(*,*) 'Failure in jet production. Resampling is required'
        case(115)
        write(*,*) 'Unknown parton distribution. Resampling is required'
        case(116)
        write(*,*) 'Too low CM energy. Resampling is required'
        case(119)
        write(*,*) '(g,N)pi reaction skipped'
        case(120)
         write(*,*) 'ATIMA calculates stopping power of particle with
     &   charge <= 0 '
       end select

       endif
       ncascerr = -1
       return
      endif
! End of revision

      stop
      end

