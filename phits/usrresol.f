************************************************************************
*                                                                      *
      function usrdefres(heat)
*                                                                      *
*        subroutine for arbitrary energy resolution for [t-deposit]    *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

              if(heat .lt. 0.3d0) then ! incoming radiation is less than 300 keV. Switch amplitude parameter
               aa = 1.d0
              else
               aa = 0.9371d0 - 1.4142d-1 * heat
              endif
              bb  = 3.834d-13 * (heat*1.d3)**4 - 1.462d-9 *
     &         (heat*1.d3) ** 3 + 2.087d-6 * (heat*1.d3)**2
     &          - 1.3875d-3 * (heat*1.d3) + 0.4387d0
              alp = -0.0807d0 * Exp(-0.00178d0 * heat * 1.d3)
              do i = 1, 100000
               heatdif = heat * (unirn(dummy) -0.5d0) * 1.d2 ! energy spread range

               if((exp(alp * (heatdif)**2) .gt. 2.d0 * unirn(dummy).and.
     &          heatdif .ge. 0.d0 )
     &          .or.
     &            (exp(alp * (heatdif)**2) + aa * exp(bb * heatdif) *
     &          (1.d0 - exp(0.6d0 * alp * heatdif**2)) )
     &          .gt. 2.d0 * unirn(dummy)  .and. heatdif .lt. 0.d0 ) exit
              enddo
              usrdefres = heat + heatdif*1.d-3

      return
      end
