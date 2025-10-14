************************************************************************
*                                                                      *
      function sect_a(stime)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to get elapse time                                      *
*              use function : secnds(stime)                            *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              stime    : starting time                                *
*              sect_a   : elapse time from starting time               *
*                                                                      *
*                                                                      *
************************************************************************

      implicit double precision(a-h, o-z)

      real*4 secnds

*-----------------------------------------------------------------------

            sect_a = dble( secnds( real( stime , 4 ) ) )

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      function cput_a(cputm)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to get cpu time                                         *
*              use function : dtime(tarray)                            *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              cputm    : starting time                                *
*              cput_a   : cpu time from starting time                  *
*                                                                      *
*                                                                      *
************************************************************************

      implicit double precision(a-h, o-z)

      real*4 secnds

*-----------------------------------------------------------------------

      data cpuall / 0.0 /
      save cpuall

*-----------------------------------------------------------------------

            cput_a = dble( secnds( real( cputm , 4 ) ) )

*-----------------------------------------------------------------------


      return
      end


************************************************************************
*                                                                      *
      subroutine date_a_time(iyer,imon,iday,ihor,imin,isec)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to get current date                                     *
*              use subroutine date_and_time                            *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              iyer     : 4-digit year                                 *
*              imon     : month                                        *
*              iday     : day                                          *
*              ihor     : hour                                         *
*              imin     : minite                                       *
*              isec     : secont                                       *
*                                                                      *
*                                                                      *
************************************************************************

      dimension ivalue(8)
      character*12 dum(3)

*-----------------------------------------------------------------------

            call date_and_time(dum(1),dum(2),dum(3),ivalue)

            iyer = ivalue(1)
            imon = ivalue(2)
            iday = ivalue(3)
            ihor = ivalue(5)
            imin = ivalue(6)
            isec = ivalue(7)

*-----------------------------------------------------------------------

      return
      end


