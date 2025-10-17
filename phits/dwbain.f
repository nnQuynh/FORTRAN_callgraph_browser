************************************************************************
*                                                                      *
      subroutine dwbain(ityp,eein,mmas,mchg)
*                                                                      *
*                                                                      *
*       event generate using results of DWBA calculation               *
*       modified by S.Hashimoto on 2014/02/19                          *
*                                                                      *
*        call subroutine : dwbaD, dwbaE                                *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       idwba   : =0 : not using this option                           *
*                 =1 : using this option                               *
*       ityp    : particle type of projectile                          *
*       eein    : energy of projectile (MeV)                           *
*       mmas    : mass of target                                       *
*       mchg    : charge of target                                     *
*                                                                      *
*     output:                                                          *
*                                                                      *
*       idwbahit: =0 : dwba event is not chosen                        *
*                 =1 : dwba event is chosen                            *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*        nclst   : total number of out going particles and nuclei      *
*                                                                      *
*        iclust(nclst)                                                 *
*                                                                      *
*                i = 0, nucleus                                        *
*                  = 1, proton                                         *
*                  = 2, neutron                                        *
*                  = 3, pion                                           *
*                  = 4, photon                                         *
*                  = 5, kaon                                           *
*                  = 6, muon                                           *
*                  = 7, others                                         *
*                                                                      *
*        jclust(i,nclst)                                               *
*                                                                      *
*                i = 0, angular momentum                               *
*                  = 1, proton number                                  *
*                  = 2, neutron number                                 *
*                  = 3, ip, see below                                  *
*                  = 4,                                                *
*                  = 5, charge                                         *
*                  = 6, baryon number                                  *
*                  = 7, kf code                                        *
*                  = 8, isomer level (0: Ground, 1,2: 1st, 2nd isomer) *
*                                                                      *
*        qclust(i,nclst)                                               *
*                                                                      *
*                i = 0, impact parameter                               *
*                  = 1, px (GeV/c)                                     *
*                  = 2, py (GeV/c)                                     *
*                  = 3, pz (GeV/c)                                     *
*                  = 4, etot = sqrt( p**2 + rm**2 ) (GeV)              *
*                  = 5, rest mass (GeV)                                *
*                  = 6, excitation energy (MeV)                        *
*                  = 7, kinetic energy (MeV)                           *
*                  = 8, weight change                                  *
*                  = 9, delay time                                     *
*                  = 10, x-displace                                    *
*                  = 11, y-displace                                    *
*                  = 12, z-displace                                    *
*                                                                      *
************************************************************************
      use NGSDATAMOD, only : bindeg

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'

*-----------------------------------------------------------------------

      common /cidwba/ idwba
      common /cdwba/  bdwba, idwbahit
!$OMP THREADPRIVATE(/cdwba/)

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)

*-----------------------------------------------------------------------

      data init_dwba/1/
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

      dimension rmspat(20)
      data rmspat/ 0.93827d0, 0.93958d0, 0.1396d0, 0.1350d0, 0.1396d0,
     &             0.0d0, 0.0d0, 0.4936d0, 0.4977d0, 0.4936d0, 4*0.0d0,
     &             1.8756d0, 2.8089d0, 2.8084d0, 3.7274d0, 0.0d0, 0.0d0/

*-----------------------------------------------------------------------

      DATA PAI/3.141592653589793D0/

*-----------------------------------------------------------------------
*        initialization
*-----------------------------------------------------------------------

         bdwba=0d0
         bimpacti=0d0

*-----------------------------------------------------------------------
*        condition check
*-----------------------------------------------------------------------

* deuteron induced reactions
      if (ityp .eq. 15) then

* for 6,7Li(d,n)7,8Be and 6,7Li(d,p)7,8Li at 10 - 50 MeV
         if (mchg .eq. 3 .and. (mmas .eq. 6 .or. mmas .eq. 7)
     1        .and. eein .ge. 10d0 .and. eein .le. 50d0) then
            go to 1000

* for 9Be(d,n)10B and 9Be(d,p)10Be at 5 - 25 MeV
         else if (mchg .eq. 4 .and. mmas .eq. 9
     1           .and. eein .ge. 5d0 .and. eein .le. 25d0) then
            go to 1000

* for 12,13C(d,n)13,14N and 12,13C(d,p)13,14C at 10 - 50 MeV
         else if (mchg .eq. 6 .and. (mmas .eq. 12 .or. mmas .eq. 13)
     1           .and. eein .ge. 10d0 .and. eein .le. 50d0) then
            go to 1000

         else
            go to 9000
         end if


* proton induced reactions
      else if (ityp .eq. 1) then

* for 7Li(p,n)7Be at 30 - 400 MeV
         if (mchg .eq. 3 .and. mmas .eq. 7
     1        .and. eein .ge. 30d0 .and. eein .le. 400d0) then
            go to 1000

* for 9Be(p,n)9B at 10 - 50 MeV
         else if (mchg .eq. 4 .and. mmas .eq. 9
     1           .and. eein .ge. 10d0 .and. eein .le. 50d0) then
            go to 1000

         else
            go to 9000
         end if

      else
         go to 9000
      end if

 1000 continue

*-----------------------------------------------------------------------
*        data file open
*-----------------------------------------------------------------------

!$OMP CRITICAL (dwbaD_crit)
      if ( init_dwba .eq. 1) then

         init_dwba=0
         call dwbaD

      end if
!$OMP END CRITICAL (dwbaD_crit)

*-----------------------------------------------------------------------
*        cross section
*-----------------------------------------------------------------------

         ienep=0
         TSXSDWBA=0d0

* deuteron induced reactions
      if (ityp .eq. 15) then

      ap=2d0
      zp=1d0
      at=dble(mmas)
      zt=dble(mchg)
      call sighi(ap,zp,eein,at,zt,signe,sigel,bmax)

* for 6Li(d,n)7Be and 6Li(d,p)7Li at 10 - 50 MeV
      if (mmas .eq. 6) then
         iinich=1
 2006    continue
         ienep=ienep+1
         if (SXSDWBA(2,ienep,1,iinich) .lt. eein) go to 2006
         ifinch=idnint(DPINI(5,iinich))
         do id=1,ifinch
         if (ienep.gt.1) then
            TSXSDWBA=TSXSDWBA
     &      +(SXSDWBA(1,ienep,id,iinich)-SXSDWBA(1,ienep-1,id,iinich))
     &      /(SXSDWBA(2,ienep,id,iinich)-SXSDWBA(2,ienep-1,id,iinich))
     &      *(eein-SXSDWBA(2,ienep-1,id,iinich))
     &      +SXSDWBA(1,ienep-1,id,iinich)
         else
            TSXSDWBA=TSXSDWBA+SXSDWBA(1,ienep,id,iinich)
         end if
         end do

* for 7Li(d,n)8Be and 7Li(d,p)8Li at 10 - 50 MeV
      else if (mmas .eq. 7) then
         iinich=2
 2007    continue
         ienep=ienep+1
         if (SXSDWBA(2,ienep,1,iinich) .lt. eein) go to 2007
         ifinch=idnint(DPINI(5,iinich))
         do id=1,ifinch
         if (ienep.gt.1) then
            TSXSDWBA=TSXSDWBA
     &      +(SXSDWBA(1,ienep,id,iinich)-SXSDWBA(1,ienep-1,id,iinich))
     &      /(SXSDWBA(2,ienep,id,iinich)-SXSDWBA(2,ienep-1,id,iinich))
     &      *(eein-SXSDWBA(2,ienep-1,id,iinich))
     &      +SXSDWBA(1,ienep-1,id,iinich)
         else
            TSXSDWBA=TSXSDWBA+SXSDWBA(1,ienep,id,iinich)
         end if
         end do

* for 9Be(d,n)10B and 9Be(d,p)10Be at 5 - 25 MeV
      else if (mmas .eq. 9) then
         iinich=3
 2009    continue
         ienep=ienep+1
         if (SXSDWBA(2,ienep,1,iinich) .lt. eein) go to 2009
         ifinch=idnint(DPINI(5,iinich))
         do id=1,ifinch
         if (ienep.gt.1) then
            TSXSDWBA=TSXSDWBA
     &      +(SXSDWBA(1,ienep,id,iinich)-SXSDWBA(1,ienep-1,id,iinich))
     &      /(SXSDWBA(2,ienep,id,iinich)-SXSDWBA(2,ienep-1,id,iinich))
     &      *(eein-SXSDWBA(2,ienep-1,id,iinich))
     &      +SXSDWBA(1,ienep-1,id,iinich)
         else
            TSXSDWBA=TSXSDWBA+SXSDWBA(1,ienep,id,iinich)
         end if
         end do

* for 12C(d,n)13N and 12C(d,p)13C at 10 - 50 MeV
      else if (mmas .eq. 12) then
         iinich=4
 2012    continue
         ienep=ienep+1
         if (SXSDWBA(2,ienep,1,iinich) .lt. eein) go to 2012
         ifinch=idnint(DPINI(5,iinich))
         do id=1,ifinch
         if (ienep.gt.1) then
            TSXSDWBA=TSXSDWBA
     &      +(SXSDWBA(1,ienep,id,iinich)-SXSDWBA(1,ienep-1,id,iinich))
     &      /(SXSDWBA(2,ienep,id,iinich)-SXSDWBA(2,ienep-1,id,iinich))
     &      *(eein-SXSDWBA(2,ienep-1,id,iinich))
     &      +SXSDWBA(1,ienep-1,id,iinich)
         else
            TSXSDWBA=TSXSDWBA+SXSDWBA(1,ienep,id,iinich)
         end if
         end do

* for 13C(d,n)14N and 13C(d,p)14C at 10 - 50 MeV
      else if (mmas .eq. 13) then
         iinich=5
 2013    continue
         ienep=ienep+1
         if (SXSDWBA(2,ienep,1,iinich) .lt. eein) go to 2013
         ifinch=idnint(DPINI(5,iinich))
         do id=1,ifinch
         if (ienep.gt.1) then
            TSXSDWBA=TSXSDWBA
     &      +(SXSDWBA(1,ienep,id,iinich)-SXSDWBA(1,ienep-1,id,iinich))
     &      /(SXSDWBA(2,ienep,id,iinich)-SXSDWBA(2,ienep-1,id,iinich))
     &      *(eein-SXSDWBA(2,ienep-1,id,iinich))
     &      +SXSDWBA(1,ienep-1,id,iinich)
         else
            TSXSDWBA=TSXSDWBA+SXSDWBA(1,ienep,id,iinich)
         end if
         end do

* for the case that there is no data.
      else
         TSXSDWBA=0d0
      end if


* proton induced reactions
      else if (ityp .eq. 1) then

      call sigrc(ityp,eein,mmas,mchg,sigt,signe,sigel)

* for 7Li(p,n)7Be at 30 - 400 MeV
      if (mmas .eq. 7) then
         iinich=6
 3007    continue
         ienep=ienep+1
         if (SXSDWBA(2,ienep,1,iinich) .lt. eein) go to 3007
         ifinch=idnint(DPINI(5,iinich))
         do id=1,ifinch
         if (ienep.gt.1) then
            TSXSDWBA=TSXSDWBA
     &      +(SXSDWBA(1,ienep,id,iinich)-SXSDWBA(1,ienep-1,id,iinich))
     &      /(SXSDWBA(2,ienep,id,iinich)-SXSDWBA(2,ienep-1,id,iinich))
     &      *(eein-SXSDWBA(2,ienep-1,id,iinich))
     &      +SXSDWBA(1,ienep-1,id,iinich)
         else
            TSXSDWBA=TSXSDWBA+SXSDWBA(1,ienep,id,iinich)
         end if
         end do
         if ( eein .gt. 50d0 ) then
            signe=signe*2.4d0   !To reproduce absolute values of peaks
         end if

* for 9Be(p,n)9B at 10 - 50 MeV
      else if (mmas .eq. 9) then
         iinich=7
 3009    continue
         ienep=ienep+1
         if (SXSDWBA(2,ienep,1,iinich) .lt. eein) go to 3009
         ifinch=idnint(DPINI(5,iinich))
         do id=1,ifinch
         if (ienep.gt.1) then
            TSXSDWBA=TSXSDWBA
     &      +(SXSDWBA(1,ienep,id,iinich)-SXSDWBA(1,ienep-1,id,iinich))
     &      /(SXSDWBA(2,ienep,id,iinich)-SXSDWBA(2,ienep-1,id,iinich))
     &      *(eein-SXSDWBA(2,ienep-1,id,iinich))
     &      +SXSDWBA(1,ienep-1,id,iinich)
         else
            TSXSDWBA=TSXSDWBA+SXSDWBA(1,ienep,id,iinich)
         end if
         end do

* for the case that there is no data.
      else
         TSXSDWBA=0d0
      end if

      end if

         bdwba=dsqrt(TSXSDWBA/PAI/10d0)
         XSratio=TSXSDWBA/(signe*1d3)
         Xrndm=rn(0)
      if (Xrndm .gt. XSratio) go to 9000
      idwbahit=1

*-----------------------------------------------------------------------

               eppin = eein / 1000.0d0
               rmsin = rmspat(ityp)
               if ( rmsin .le. 0. ) go to 340

            if( ityp .eq. 1 ) then

               masim = mmas + 1
               mchim = mchg + 1
               einad = 0.0d0

            else if( ityp .eq. 2 ) then

               masim = mmas + 1
               mchim = mchg
               einad = 0.0d0

            else if( ityp .eq. 3 ) then

               masim = mmas
               mchim = mchg + 1
               einad = 139.9d+0

            else if( ityp .eq. 4 ) then

               masim = mmas
               mchim = mchg
               einad = 135.1d+0

            else if( ityp .eq. 5 ) then

               masim = mmas
               mchim = mchg - 1
               einad = 139.9d+0

            else if( ityp .eq. 15 ) then

               masim = mmas + 2
               mchim = mchg + 1
               einad = -2.225d0

            else if( ityp .eq. 16 ) then

               masim = mmas + 3
               mchim = mchg + 1
               einad = -8.483d0

            else if( ityp .eq. 17 ) then

               masim = mmas + 3
               mchim = mchg + 2
               einad = -7.719d0

            else if( ityp .eq. 18 ) then

               masim = mmas + 4
               mchim = mchg + 2
               einad = -28.297d0

            else

               goto 340

            end if

*-----------------------------------------------------------------------
*        DWBA
*-----------------------------------------------------------------------

      call dwbaE(ityp,eein,mmas,mchg)

*-----------------------------------------------------------------------
*        booking of the result of DWBA
*-----------------------------------------------------------------------

               nbart = 0
               nchat = 0

               npipo = 0
               nping = 0
               npine = 0

               sume = 0.0d0
               tbene = 0.0d0

               poutx = 0.0d0
               pouty = 0.0d0
               poutz = 0.0d0

*-----------------------------------------------------------------------
*        pseudo collision
*-----------------------------------------------------------------------

               nclst = 0

         if( nopart .lt. 0 ) then

               nclst = -1

               return

         end if

*-----------------------------------------------------------------------
*        real collision
*-----------------------------------------------------------------------

         if( nopart .gt. 0 ) then

            do n = 1, nopart

                  nclst = nclst + 1

               if( lk .eq. 1 ) then

                  kf    = 2212
                  ibary = 1
                  ipid  = 1
                  ippad = 1
                  ipprt = 1
                  ipneu = 0
                  ipchg = 1

                  bene = 0.0d0
                  rms = 0.93827d0

               else if( lk .eq. 2 ) then

                  kf    = 2112
                  ibary = 1
                  ipid  = 2
                  ippad = 2
                  ipprt = 0
                  ipneu = 1
                  ipchg = 0

                  bene = 0.0d0
                  rms = 0.93958d0

               else if( lk .eq. 3 ) then

                  npipo = npipo + 1

                  kf    = 211
                  ibary = 0
                  ipid  = 3
                  ippad = 3
                  ipprt = 0
                  ipneu = 0
                  ipchg = 1

                  bene = 0.0d0
                  rms = 0.1396d0

               else if( lk .eq. 4 ) then

                  npine = npine + 1

                  kf    = 111
                  ibary = 0
                  ipid  = 3
                  ippad = 4
                  ipprt = 0
                  ipneu = 0
                  ipchg = 0

                  bene = 0.0d0
                  rms = 0.1350d0

               else if( lk .eq. 5 ) then

                  nping = nping + 1

                  kf    = -211
                  ibary = 0
                  ipid  = 3
                  ippad = 5
                  ipprt = 0
                  ipneu = 0
                  ipchg = -1

                  bene = 0.0d0
                  rms = 0.1396d0

               else if( lk .eq. 6 ) then

                  kf    = 1000002
                  ibary = 2
                  ipid  = 0
                  ippad = 15
                  ipprt = 1
                  ipneu = 1
                  ipchg = 1

                  bene = bindeg(ipprt,ipneu)
                  rms = 1.8756d0

               else if( lk .eq. 7 ) then

                  kf    = 1000003
                  ibary = 3
                  ipid  = 0
                  ippad = 16
                  ipprt = 1
                  ipneu = 2
                  ipchg = 1

                  bene = bindeg(ipprt,ipneu)
                  rms = 2.8089d0

               else if( lk .eq. 8 ) then

                  kf    = 2000003
                  ibary = 3
                  ipid  = 0
                  ippad = 17
                  ipprt = 2
                  ipneu = 1
                  ipchg = 2

                  bene = bindeg(ipprt,ipneu)
                  rms = 2.8084d0

               else if( lk .eq. 9 ) then

                  kf    = 2000004
                  ibary = 4
                  ipid  = 0
                  ippad = 18
                  ipprt = 2
                  ipneu = 2
                  ipchg = 2

                  rms = 3.7274d0
                  bene = bindeg(ipprt,ipneu)

               else if( lk .le. 128 ) then

                  ipid  = 0
                  ipprt = mod(lk,10)
                  ipneu = (lk-ipprt)/10 + 0.1 - ipprt
                  ippad = 19
                  ipchg = ipprt
                  ibary = ipprt + ipneu
                  kf    = 1000000 * ipprt + ibary

                  bene = bindeg(ipprt,ipneu)
                  rms = 0.93827d0 * ipprt + 0.93958d0 * ipneu
     &                - bene / 1000.0d0

               else

                  goto 340

               end if

*-----------------------------------------------------------------------

                  nbart = nbart + ibary
                  nchat = nchat + ipchg
                  tbene = tbene + bene

                  sume = sume + epo

                  epp = epo / 1000.0d0

                  pouta = sqrt( epp**2 + 2.0d0 * epp * rms )

                  pxrv  = pouta * alphao
                  pyrv  = pouta * betao
                  pzrv  = pouta * gammao

                  poutx = poutx + pxrv
                  pouty = pouty + pyrv
                  poutz = poutz + pzrv

*-----------------------------------------------------------------------
*        booking of outgoing particles
*-----------------------------------------------------------------------

                  iclust(nclst)    = ipid

                  jclust(0,nclst)  = 0
                  jclust(1,nclst)  = ipprt
                  jclust(2,nclst)  = ipneu
                  jclust(3,nclst)  = ippad
                  jclust(4,nclst)  = 0
                  jclust(5,nclst)  = ipchg
                  jclust(6,nclst)  = ibary
                  jclust(7,nclst)  = kf
                  jclust(8,nclst)  = 0

                  qclust(0,nclst)  = bimpacti
                  qclust(1,nclst)  = pxrv
                  qclust(2,nclst)  = pyrv
                  qclust(3,nclst)  = pzrv
                  qclust(4,nclst)  = epp + rms
                  qclust(5,nclst)  = rms
                  qclust(6,nclst)  = 0.0
                  qclust(7,nclst)  = epp * 1000.d0
                  qclust(8,nclst)  = 1.0
                  qclust(9,nclst)  = 0.0
                  qclust(10,nclst) = 0.0d0
                  qclust(11,nclst) = 0.0d0
                  qclust(12,nclst) = 0.0d0

            end do

         end if

*-----------------------------------------------------------------------
*           residual nucleus
*-----------------------------------------------------------------------

                  masrs = masim - nbart
                  mchrs = mchim - nchat

            if( masrs .gt. 0 .and. mchrs .ge. 0 ) then

                  nclst = nclst + 1

                  piabs = sqrt( eppin**2 + 2.0d0 * eppin * rmsin )

                  presx = - poutx
                  presy = - pouty
                  presz = - poutz + piabs

                  pabst = sqrt( presx**2 + presy**2 + presz**2 )
                  rsmas = DPFIN(7,ifinch,iinich)*931.494013d0/1000d0

                  etota = sqrt( pabst**2 + rsmas**2 )

                  erres = ( etota - rsmas ) * 1000.0d0

                  exres = eein + einad
     &                  - sume
     &                  - 139.9d+0 * ( npipo + nping )
     &                  - 135.1d+0 * npine
     &                  - erres
     &                  - bindeg(mchg,mmas-mchg)
     &                  + bindeg(mchrs,masrs-mchrs)
     &                  + tbene

                  if( exres .lt. -1.0d0 ) then

                     nclst = -1
                     return

                  end if

                  exres = max( 0.0d0, exres )

*-----------------------------------------------------------------------
*           booking of the residual nucleus
*-----------------------------------------------------------------------

                  iclust(nclst)    = 0

                  jclust(0,nclst)  = 0
                  jclust(1,nclst)  = mchrs
                  jclust(2,nclst)  = masrs - mchrs
                  jclust(3,nclst)  = 19
                  jclust(4,nclst)  = 0
                  jclust(5,nclst)  = mchrs
                  jclust(6,nclst)  = masrs
                  jclust(7,nclst)  = mchrs * 1000000 + masrs
                  jclust(8,nclst)  = 0

                  qclust(0,nclst)  = bimpacti
                  qclust(1,nclst)  = presx
                  qclust(2,nclst)  = presy
                  qclust(3,nclst)  = presz
                  qclust(4,nclst)  = etota
                  qclust(5,nclst)  = rsmas
                  qclust(6,nclst)  = exres
                  qclust(7,nclst)  = ( etota - rsmas ) * 1000.d0
                  qclust(8,nclst)  = 1.0
                  qclust(9,nclst)  = 0.0
                  qclust(10,nclst) = 0.0d0
                  qclust(11,nclst) = 0.0d0
                  qclust(12,nclst) = 0.0d0

            end if

*-----------------------------------------------------------------------
*           x-y rotation by randum
*-----------------------------------------------------------------------

                  theta = 2.0d0 * PAI * rn(0)

               do i = 1, nclst

                  px = qclust(1,i)
                  py = qclust(2,i)

                  qclust(1,i)  = px * cos( theta ) - py * sin( theta )
                  qclust(2,i)  = px * sin( theta ) + py * cos( theta )

               end do

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------
*        error in dwba
*-----------------------------------------------------------------------

  340 continue

               write(6,1002) ityp,lk
 1002          format(/' *** error message from dwba ***'
     &         /' invalid condition of itype or lk was found.'
     &         /' ityp =',i5,'  lk =',i5)
               call parastop( 844 )

*-----------------------------------------------------------------------

 9000 continue
      idwbahit=0
      return
      end
