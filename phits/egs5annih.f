************************************************************************
*                                                                      *
      subroutine egs5annih(ein)
*                                                                      *
*       positron annihilation
*       last modified by H.Iwase on 2012/3/6
*                                                                      *
************************************************************************
      implicit none
!      save
      include 'include/egs5_h.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'
      include 'include/egs5_epcont.f'

      real*8 rnnow,ein
      integer ircode

      integer ncpls

      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)
      real*8 atmrc

      if ( edep .lt. ein ) then

      if( np .eq. 0 ) np = 1

      call randomset(rnnow)
      costhe = rnnow
      call randomset(rnnow)
      if (rnnow .le. 0.5) costhe = -costhe
      sinthe = sqrt(1. - costhe**2)
      e(np) = RM
      iq(np) = 0

      u(np) = 0.                ! Make photon go along z-axis
      v(np) = 0.
      w(np) = 1.

      call uphi(2,1)            ! Set direction cosines

      uf(np) = 0.
      uf(np+1) = 0.
      vf(np) = 0.
      vf(np+1) = 0.
      wf(np) = 0.
      wf(np+1) = 0.

!     Set up second photon in opposite direction

      np = np + 1
      e(np) = RM
      iq(np) = 0
      ir(np) = ir(np-1)
      wt(np) = wt(np-1)
      u(np) = -u(np-1)
      v(np) = -v(np-1)
      w(np) = -w(np-1)

!     put second photon information to 1 for PHITS
      call npconv
!     put first photon information to 1 for PHITS
      call npconv
      atmrc(3,5) = atmrc(3,5) + 1.d0

      ircode = 2

cKN 2014/09/09, T.Sato 2015/2/26 (change 2 to 3)
           ncpls = 3

cABE 2022/02/24, save positron energy to use in egs2phits
      e(np+1) = ein

      call egs2phits(ncpls)

      endif

      np = 1
      return
      end

