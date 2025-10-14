************************************************************************
*                                                                      *
      subroutine usrmgt1(b_mag,phase,ptime,o_mag)
*                                                                      *
*        sample subroutine for time dependence of magnetic field.      *
*        This sample is for Wobbler magnet                             *
*                                                                      *
*        input :                                                       *
*           b_mag : initial magnetic field                             *
*           phase : phase of magnet given by input                     *
*           ptime : present time (msec)                                *
*                                                                      *
*        output :                                                      *
*           o_mag : final magnetic field                               *
*                                                                      *
*        common :                                                      *
*           pmphs : initial phase of this event                        *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      common /srsph/  pmphs
!$OMP THREADPRIVATE(/srsph/)

*-----------------------------------------------------------------------


          o_mag = sin( pmphs + phase ) * b_mag


*-----------------------------------------------------------------------

      return
      end

