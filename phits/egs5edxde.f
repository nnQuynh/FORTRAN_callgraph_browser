
************************************************************************
*                                                                      *
      subroutine egs5edxde(eout,step,ein,mat,ityp)
*                                                                      *
*       input                                                          *
*         step : moving length (cm)                                    *
*         ein  : initial energy (MeV)                                  *

*         mat  : material number                                       *
*         ityp : particle id                                           *
*                                                                      *
*       output                                                         *
*         eout : energy at delt (MeV)                                  *
*                                                                      *
************************************************************************
      USE   egs5_elecin_mod  !<- 2015.08xx allocatable
      USE egs5_media_mod !<- 2015.08xx allocatable
      implicit none
!      save
      include 'include/egs5_h.f'
      include 'include/egs5_epcont.f'
      include 'include/egs5_misc.f'

      real*8 ein,step,eout,ek,dedx0,dedx,de
      integer mat,ityp,lelke,irl,medium

      real*8 esmax, esmin, emin
      common /eparm/  esmax, esmin, emin(20)

      ek  = ein
      if( ek .lt. emin(ityp) ) then
       ek = emin(ityp)
      endif

      eout = ek
      medium = mat

      if(medium .gt. 0) then

         elke = log(ek)
         lelke = eke1(medium)*elke + eke0(medium)

         rhof=1.d0

c     for photon
         if (ityp .eq. 14) then

! T.Sato 2016/1/26, avoid photon deposition energy
            eout = ein

c     for e- and e+
         else

            if (ityp .eq. 12) then
               dedx0 = ededx1(lelke,medium)*elke + ededx0(lelke,medium)

            elseif (ityp .eq. 13) then
               dedx0 = pdedx1(lelke,medium)*elke + pdedx0(lelke,medium)

            endif

            dedx = rhof*dedx0

            de = dedx * step

            eout = ein - de

         endif

      endif
      end
