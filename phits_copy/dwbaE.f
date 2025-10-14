************************************************************************
*                                                                      *
      subroutine dwbaE(ityp,eein,mmas,mchg)
*                                                                      *
*                                                                      *
*       event generate using results of DWBA calculation               *
*       modified by S.Hashimoto on 2014/02/19                          *
*                                                                      *
*        called from: dwbain                                           *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      parameter (ninich=20,nfinch=17,nenep=10,nthmx=181)
      dimension DPINI(5,ninich),DPFIN(8,nfinch,ninich)
      dimension SXSDWBA(2,nenep,nfinch,ninich)
      dimension THPD(3,nenep,nfinch,ninich)
      dimension XSDWBA(nthmx,nenep,nfinch,ninich)
      common/dwbaxsdata/DPINI,DPFIN,SXSDWBA,THPD,XSDWBA

      common/dwbaevent1/iinich,ifinch,ienep
!$OMP THREADPRIVATE(/dwbaevent1/)
      common/dwbaevent2/TSXSDWBA
!$OMP THREADPRIVATE(/dwbaevent2/)

      common/dwbaoutput1/nopart,lk
!$OMP THREADPRIVATE(/dwbaoutput1/)
      common/dwbaoutput2/epo,alphao,betao,gammao
!$OMP THREADPRIVATE(/dwbaoutput2/)

*-----------------------------------------------------------------------

      DATA PAI/3.141592653589793D0/
c CODATA recommended values (1998)
      parameter (HC=197.3269601d0, AC=931.494013d0, EC=137.03599976d0)

*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*     number of data
*-----------------------------------------------------------------------

      nthetamx=idnint((THPD(2,ienep,ifinch,iinich)
     &     -THPD(1,ienep,ifinch,iinich))/THPD(3,ienep,ifinch,iinich))+1

*-----------------------------------------------------------------------
*     random number
*-----------------------------------------------------------------------

         Xrndm=rn(0)
         dweight=2d0*PAI/TSXSDWBA
      do id=1,ifinch
         if (ienep.gt.1) then
            dene=SXSDWBA(2,ienep,id,iinich)-SXSDWBA(2,ienep-1,id,iinich)
         else
            dene=1d0
         end if
         dthetaCMrad=THPD(3,ienep,id,iinich)*PAI/180d0
      do itheta=1,nthetamx
         thetaCMrad=(itheta-1)*dthetaCMrad
         thetaCMradm=thetaCMrad-dthetaCMrad/2d0
         thetaCMradp=thetaCMrad+dthetaCMrad/2d0
         if (dsin(thetaCMrad).le.1d-9) then
            dCOSTCM=1d0-dabs(dcos(thetaCMradp))
         else
            dCOSTCM=dcos(thetaCMradm)-dcos(thetaCMradp)
         end if
         if (ienep.gt.1) then
            XSINT=(XSDWBA(itheta,ienep,id,iinich)
     &           -XSDWBA(itheta,ienep-1,id,iinich))/dene
     &           *(eein-SXSDWBA(2,ienep,id,iinich))
     &           +XSDWBA(itheta,ienep,id,iinich)
         else
            XSINT=XSDWBA(itheta,ienep,id,iinich)
         end if
         weight=XSINT*dweight*dCOSTCM
         Xrndm=Xrndm-weight
         if (Xrndm .lt. 0d0) go to 1000
      end do
      end do
         nopart = -1
         go to 9000
 1000 continue

*-----------------------------------------------------------------------
*     initial channel
*-----------------------------------------------------------------------

      Aproj=DPINI(1,iinich)*AC
      Atarg=DPINI(3,iinich)*AC
      Tplab=eein
      Eplab=Tplab+Aproj
      Pplab=dsqrt(Tplab*(Tplab+2d0*Aproj))
      beta=Pplab/(Eplab+Atarg)
      gamma=1d0/dsqrt(1d0-beta**2)
      PpCM=gamma*Atarg/(Eplab+Atarg)*Pplab
      EpCM=gamma*(Aproj**2+Atarg*Eplab)/(Eplab+Atarg)
      EtCM=gamma*Atarg
      TpCM=EpCM-Aproj+EtCM-Atarg

*-----------------------------------------------------------------------
*     final channel
*-----------------------------------------------------------------------

      Aejec=DPFIN(5,id,iinich)*AC
      Aresi=DPFIN(7,id,iinich)*AC

*-----------------------------------------------------------------------

      thetaCM=(itheta-1)*THPD(3,ienep,ifinch,iinich)
      Arndm=rn(0)
      Afluc=THPD(3,ienep,ifinch,iinich)*(Arndm-0.5d0)
      if (itheta.eq.1 .and. Afluc.lt.0d0) then
         thetaCM=thetaCM+dabs(Afluc)
      else
         thetaCM=thetaCM+Afluc
      end if
      thetaCMrad=thetaCM*PAI/180d0
      COSTCM=dcos(thetaCMrad)
      SINTCM=dsin(thetaCMrad)

      Wrndm=rn(0)
C S.H. Revised to use Lorentz distribution (2013.11.15)
      Efluc=DPFIN(3,id,iinich)/2d0*dtan(PAI*(Wrndm-0.5d0))
      TfCM=TpCM+DPFIN(2,id,iinich)+Efluc
C S.H. Revised for the unphysical case (2013.11.17)
      if ( TfCM.lt.0d0 ) TfCM=0.1d0
      XTf=TfCM+Aejec+Aresi
      PfCM=dsqrt(XTf**2*(XTf**2-2d0*Aejec**2-2d0*Aresi**2)
     1     +(Aejec**2-Aresi**2)**2)/2d0/XTf
      EeCM=dsqrt(PfCM**2+Aejec**2)
      Pzelab=gamma*(PfCM*COSTCM+beta*EeCM)
      Pxelab=PfCM*SINTCM
      Pelab=dsqrt(Pzelab**2+Pxelab**2)
      Eelab=gamma*(EeCM+beta*PfCM*COSTCM)
      Telab=Eelab-Aejec

*-----------------------------------------------------------------------

      nopart = 1
      epo = Telab
      alphao = Pxelab / Pelab
      betao = 0d0
      gammao = Pzelab / Pelab

      ichg=idnint(DPFIN(6,id,iinich))
      if ( ichg .eq. 1 ) then
         lk = 1
      else if ( ichg .eq. 0) then
         lk = 2
      end if

*-----------------------------------------------------------------------

 9000 continue

      return
      end
