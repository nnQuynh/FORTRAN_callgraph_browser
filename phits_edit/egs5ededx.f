
************************************************************************
*                                                                      *
      subroutine egs5ededx(ein,mat,ityp,dedx)
*                                                                      *
*      purpose : get EGS5 dE/dx of e- or e+
*                                                                      *
*       input                                                          *
*         ein  : initial energy (MeV)                                  *
*         ir   : region number                                         *
*         mat  : material number                                       *
*         ityp : particle id                                           *
*                                                                      *
*       output                                                         *
*         dedx : dE/dx of e-/e+
*                                                                      *
************************************************************************
      USE egs5_elecin_mod !<- 2015.08xx allocatable
      USE egs5_media_mod !<- 2015.08xx allocatable
      implicit none
!      save
      include 'include/egs5_h.f'
      include 'include/egs5_epcont.f'
      include 'include/egs5_misc.f'

      real*8 ein,step,eout,ek,dedx0,dedx,de
      integer ir,mat,ityp,lelke,irl,medium

      real*8 esmax, esmin, emin
      common /eparm/  esmax, esmin, emin(20)

      dedx = 0d0

      if(ein.eq.0) return

      ek  = ein
      if( ek .lt. emin(12) ) ek = emin(12)

      medium = mat

      if(medium .gt. 0) then

         elke  = log(ek)
         lelke = eke1(medium)*elke + eke0(medium)
         rhof  = 1.d0

c     for e- and e+

            if (ityp .eq. 12) then
               dedx0 = ededx1(lelke,medium)*elke + ededx0(lelke,medium)

            elseif (ityp .eq. 13) then
               dedx0 = pdedx1(lelke,medium)*elke + pdedx0(lelke,medium)

            endif

            dedx = rhof*dedx0

      endif
      end
