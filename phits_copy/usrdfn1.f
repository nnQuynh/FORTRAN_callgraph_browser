************************************************************************
*                                                                      *
      subroutine usrdfn1(ityp,ktyp,jtyp,rtyp,
     &                    dedx,dtrk,eini,efin,dhet)
*                                                                      *
*        sample subroutine for user defined let factor for heat        *
*                                                                      *
*        input :                                                       *
*           ityp : type of particle                                    *
*           ktyp : kf code of particle                                 *
*           jtyp : charge of particle                                  *
*           rtyp : mass of particle (MeV)                              *
*           dedx : dE/dx (keV/um) at sqrt(eini*efin)                   *
*           dtrk : track (cm)                                          *
*           eini : initial energy (MeV)                                *
*           efin : final energy (MeV)                                  *
*                                                                      *
*        output :                                                      *
*           dhet : dE (default) * user defined factor                  *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)
      include 'angel01.inc'

      common /rval1/ cval(mxcval), aval(mxcval)

*-----------------------------------------------------------------------
*     IZ and IN for Nucleus
*-----------------------------------------------------------------------

      if(cval(900).ge.-299.0d0.and.cval(900).le.-200.0) then ! c900 is defined between -299 and -200
       id=nint(cval(900))   
       qlfactor=dosf(dedx,id) ! if you want to use Q(E) function, dedx should be replaced with sqrt(eini*efin)
      else ! Q(L) relationship
       if(dedx.lt.10.0) then
        qlfactor=1.0
       elseif(dedx.le.100.0) then
        qlfactor=0.32*dedx-2.2
       else
        qlfactor=300.0/sqrt(dedx)
       endif
      endif
      dhet = (eini - efin)*qlfactor

*-----------------------------------------------------------------------

      return
      end

