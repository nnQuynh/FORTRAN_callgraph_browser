************************************************************************
*                                                                      *
      subroutine usrdfn3(mat,lmat,ityp,ktyp,jtyp,rtyp,
     &     dedx,dtrk,eini,efin,dhet)
*                                                                      *
*        sample subroutine for user defined let factor for heat        *
*                                                                      *
*        input :                                                       *
*           mat  : actual material number                              *
*           lmat : specfied material number                            *
*           ityp : type of particle                                    *
*           ktyp : kf code of particle                                 *
*           jtyp : charge of particle                                  *
*           rtyp : mass of particle (MeV)                              *
*           dedx : dE/dx (keV/um) at sqrt(eini*efin) of lmat           *
*           dtrk : track (cm)                                          *
*           eini : initial energy (MeV)                                *
*           efin : final energy (MeV)                                  *
*                                                                      *
*        output :                                                      *
*           dhet : dE (default) * user defined factor                  *
*                                                                      *
*                                                                      *
************************************************************************
!     Calculate Water-equivalent dose
!     Condition:
!      Water material should be specified for lmat (with positive No.)
!      Water material should be defined so that the density is 1 g/cm3
!      Unit=0 should be chosen
!     Reference W. Chan et al. TBA

      implicit real*8 (a-h,o-z)
      include 'param.inc'
      common /celdg/  rhog(kvlmax)
      real(8) sfactor
*-----------------------------------------------------------------------
*        default:  dE = eini - efin
*-----------------------------------------------------------------------
      if(jtyp.ne.0.and.mat.gt.0)then !FURUTA20190328
       rho0=rhog(mat)
       rho=rhog(lmat)
       if(efin.gt.0d0)then
        erg=sqrt(eini*efin)
        call dedxas(erg,dedx0,mat,ityp,ktyp,jtyp,rtyp) !FURUTA20190328
        dedx0=dedx0/10.0d0
        sfactor=rho0/rho*dedx/dedx0
       else
        sfactor=rho0/rho
       endif
       dhet=(eini-efin)*sfactor
      else
       dhet=0.0d0
      endif
*-----------------------------------------------------------------------

      return
      end

