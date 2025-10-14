************************************************************************
*                                                                      *
      subroutine sigmu(ein,ita,itz,sigvpi,sigbrm,sigppd)
*                                                                      *
*        calculate cross section of muon interaction                   *
*        last modified by S.Abe on 2015/07/03                          *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input  :                                                         *
*        ein    : incident muon energy [MeV]                           *
*        ita    : mass number of target                                *
*        itz    : atomic number of target                              *
*                                                                      *
*     output :                                                         *
*        sigvpi : cross section of muon photonuclear interaction [barn]*
*        sigbrm : cross section of muon-induced bremsstrahlung [barn]  *
*        sigppd : cross section of muon-induced pair production [barn] *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*     muint  :                                                         *
*        imuint : switch of muon photonuclear interaction              *
*        imubrm : switch of muon-induced bremsstrahlung                *
*        imuppd : switch of muon-induced pair production               *
*                                                                      *
*     emurng : energy range of calculation                             *
*        emumin : minimum transferred energy [MeV]                     *
*        emumax : maximum muon energy for CS calculation [MeV]         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      common /muint/ imuint, imubrm, imuppd, imucap
      common /emurng/ emumin, emumax

      common /nsigmu/ izncnt, itzvpi(200), itnvpi(200)
      common /sigmud/ vpidat(2,50,200)
      data izncnt/0/
      data itzvpi/200*0/
      data itnvpi/200*0/
      data vpidat/20000*0.d0/

*-----------------------------------------------------------------------
*     initialize
*-----------------------------------------------------------------------

      sigvpi = 0.d0
      sigbrm = 0.d0
      sigppd = 0.d0

*-----------------------------------------------------------------------
*     calculate cross section of virtual photonuclear interaction
*-----------------------------------------------------------------------

      if( imuint .eq. 1 .and.
     &    ein .ge. emumin .and. ein .le. emumax ) then

       emuc1 = 2.d+2
       emuc2 = 2.d+3
       imu1 = 20
       ddemu1 = dlog(emuc2/emuc1) / dble(imu1)

       if( emumin .lt. emuc1 ) then
        imu0 = 10
        ddemu0 = dlog(emuc1/emumin) / dble(imu0)
       else
        imu0 = 0
       endif

       if( emumax .gt. emuc2 ) then
        imu2 = 10
        ddemu2 = dlog(emumax/emuc2) / dble(imu2)
       else
        imu2 = 0
       endif

*-----------------------------------------------------------------------
*     check data table already exist or not for each material
*-----------------------------------------------------------------------
!$OMP CRITICAL (sigmu_crit)

       itn = ita - itz
       izn = 0

       if( izncnt .gt. 0 ) then
        do j = 1, izncnt
         if( itz.eq.itzvpi(j) .and. itn.eq.itnvpi(j) ) then
          izn = j
          exit
         endif
        enddo
       endif

*-----------------------------------------------------------------------
*     make cross section data table for each material
*-----------------------------------------------------------------------

       if( izn .eq. 0 ) then

        izncnt = izncnt + 1
        izn = izncnt
        itzvpi(izn) = itz
        itnvpi(izn) = itn

        if( emumin .lt. emuc1 ) then
         do j = 1, imu0
          emukin = emumin*dexp((dble(j-1))*ddemu0)
          call sigmuvpi(emukin,itz,itn,sigtmp)
          vpidat(1,j,izn) = emukin
          vpidat(2,j,izn) = sigtmp
         enddo
        endif

        do j = 1, imu1+1
         emukin = emuc1*dexp((dble(j-1))*ddemu1)
         call sigmuvpi(emukin,itz,itn,sigtmp)
         vpidat(1,imu0+j,izn) = emukin
         vpidat(2,imu0+j,izn) = sigtmp
        enddo

        if( emumax .gt. emuc2 ) then
         do j = 1, imu2
          emukin = emuc2*dexp((dble(j))*ddemu2)
          call sigmuvpi(emukin,itz,itn,sigtmp)
          vpidat(1,imu0+imu1+1+j,izn) = emukin
          vpidat(2,imu0+imu1+1+j,izn) = sigtmp
         enddo
        endif

       endif

*-----------------------------------------------------------------------
*     sampling cross section from data table
*-----------------------------------------------------------------------

       imut = imu0 + imu1 + imu2 + 1

       if( ein .ge. dmax1(emuc2,emumax) ) then
        sigvpi = ( ein - vpidat(1,imut,izn) )
     &          * ( vpidat(2,imut,izn) - vpidat(2,imut-1,izn) )
     &          / ( vpidat(1,imut,izn) - vpidat(1,imut-1,izn) )
     &          + vpidat(2,imut,izn)
       else
        do j = 1, imut+1
         if( ein .ge. vpidat(1,j,izn) .and.
     &       ein .lt. vpidat(1,j+1,izn) ) then
          sigvpi = ( ein - vpidat(1,j,izn) )
     &            * ( vpidat(2,j+1,izn) - vpidat(2,j,izn) )
     &            / ( vpidat(1,j+1,izn) - vpidat(1,j,izn) )
     &            + vpidat(2,j,izn)
          exit
         endif
        enddo
       endif

*-----------------------------------------------------------------------
!$OMP END CRITICAL (sigmu_crit)

      endif

*-----------------------------------------------------------------------
*     calculate cross section of muon bremsstrahlung
*-----------------------------------------------------------------------

      if( imubrm .eq. 1 ) then

       call sigmubrm(ein,itz,ita,sigbrm)

      endif

*-----------------------------------------------------------------------
*     calculate cross section of muon induced pair production
*-----------------------------------------------------------------------

      if( imuppd .eq. 1 ) then

       call sigmuppd(ein,itz,ita,sigppd)

      endif

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      subroutine sigmuvpi(emukin,itarz,itarn,xsmu)
*                                                                      *
*        calculate cross section of muon nuclear reaction              *
*        last modified by S.Abe on 2014/08/04                          *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input  :                                                         *
*        emukin : incident muon energy [MeV]                           *
*        itarz  : atomic number of target                              *
*        itarn  : neutron number of target                             *
*                                                                      *
*     output :                                                         *
*        xsmu   : muon-nucleus reaction cross section [barn]           *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*     muflag :                                                         *
*        imuinthit : flag for muon photonuclear interaction            *
*                                                                      *
*     emurng : energy range of calculation                             *
*        emumin : minimum transferred energy                           *
*        emumax : maximum muon energy for CS calculation               *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /emurng/ emumin, emumax
      common /muflag/ imuinthit,imubrmhit,imuppdhit,imucaphit,imucapflag
!$OMP THREADPRIVATE(/muflag/)

*-----------------------------------------------------------------------

      dimension inua(6)

      dimension enuloa(6)
      dimension enuhia(6)

      parameter ( rmumas = 105.6583715d0 )      ! muon mass [MeV/c**2]

*-----------------------------------------------------------------------
*     initialize
*-----------------------------------------------------------------------

      xsmu = 0.d0
      if( emukin .lt. 0.5d0 ) return
      if( emukin .lt. emumin ) return

      imuinthit = 1

      emutot = emukin + rmumas
      enu = 0.d0
      fintq2 = 0.d0

*-----------------------------------------------------------------------

      enumin = emumin
      enumax = emukin

      if( enumax .le. 150.0d0 ) then
       phxsck = phxs(itarz, itarn, enumax)
       if( phxsck .le. 0.d0 ) then
        imuinthit = 0
        goto 9999
       endif
      endif

*-----------------------------------------------------------------------
*     integration energy bin for transferred energy 
*-----------------------------------------------------------------------

      inua(1) = 50
      enuloa(1) = 0.5d0
      enuhia(1) = 60.0d0

      inua(2) = 5
      enuloa(2) = enuhia(1)
      enuhia(2) = 200.0d0

      inua(3) = 10
      enuloa(3) = enuhia(2)
      enuhia(3) = 500.0d0

      inua(4) = 20
      enuloa(4) = enuhia(3)
      enuhia(4) = 1500.0d0

      inua(5) = 5
      enuloa(5) = enuhia(4)
      enuhia(5) = 2000.0d0

      inua(6) = 5
      enuloa(6) = enuhia(5)
      enuhia(6) = enumax

      inuflg = 0

*-----------------------------------------------------------------------
*     search first bin for target material
*-----------------------------------------------------------------------

      inu = 100
      ddenu = dlog(enumax/enumin) / dble(inu)

      do i = 0, inu

       enuchk = enumin * exp((dble(i))*ddenu)
       phxsck = phxs(itarz,itarn,enuchk)

       if( phxsck .gt. 0.d0 ) then
        goto 100
       endif

      enddo

      imuinthit = 0
      goto 9999

  100 continue

*-----------------------------------------------------------------------

      do i = 1, 6

       if( enuchk.ge.enuloa(i) .and. enuchk.lt.enuhia(i) ) then

        inuflg = i
        enu = enuchk

        call finteg(emukin,itarz,itarn,enu,fintq2)

        inu = inua(inuflg)
        enulo = enuloa(inuflg)
        enuhi = enuhia(inuflg)
        ddenu = dlog(enuhi/enulo) / (dble(inu))

        do j = 0, inu

         enuchk = enulo * dexp(dble(j)*ddenu)

         if( enu .lt. enuchk ) then
          ista = j
          exit
         endif

        enddo

        exit

       endif

      enddo

      if( inuflg .eq. 0 ) then
       imuinthit = 0
       goto 9999
      endif

*-----------------------------------------------------------------------
*     log-bin for enu [MeV]
*-----------------------------------------------------------------------

  200 continue

      do i = ista, inu

       enuold = enu
       finold = fintq2

       enu = enulo * dexp(dble(i)*ddenu)
       enu = dmin1(enu, enumax)

       denu = enu - enuold

       call finteg(emukin,itarz,itarn,enu,fintq2)
       xsmu = xsmu + (finold + fintq2) * denu * 0.5d0

       if( enu .ge. enumax ) then
        imuinthit = 0
        goto 9999
       endif

      enddo

*-----------------------------------------------------------------------

      inuflg = inuflg + 1
      ista = 1

      if( inuflg .le. 6 ) then

       inu = inua(inuflg)
       enulo = enuloa(inuflg)
       enuhi = enuhia(inuflg)
       ddenu = dlog(enuhi/enulo) / (dble(inu))

       goto 200

      endif

*-----------------------------------------------------------------------
 9999 continue

      return
      end


************************************************************************
*                                                                      *
      subroutine vpimu(eein,wgti,ztar,atar,ireg,imat)
*                                                                      *
*        calculate muon photonuclear interaction                       *
*                                                                      *
*        last modified by S.Abe on 2015/06/24                          *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input  :                                                         *
*        eein   : incident muon energy [MeV]                           *
*        wgti   : weight                                               *
*        ztar   : atomic number of target                              *
*        atar   : mass number of target                                *
*        ireg   : region number, for sctpni                            *
*        imat   : material number, for sctpni                          *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*     muint  :                                                         *
*        imuint : switch of muon photonuclear interaction              *
*                                                                      *
*     muflag :                                                         *
*        imuinthit : flag for muon photonuclear interaction            *
*                                                                      *
*     tarzmn : to share with subroutine sctpni                         *
*        itarz  : proton number of target material                     *
*        itarm  : mass number of target material                       *
*        itarn  : neutron number of target material                    *
*                                                                      *
************************************************************************
      use NGSDATAMOD, only : bindeg

      implicit real*8 (a-h,o-z)

      include 'param00.inc'

*-----------------------------------------------------------------------

      common /muint/ imuint, imubrm, imuppd, imucap
      common /muflag/ imuinthit,imubrmhit,imuppdhit,imucaphit,imucapflag
!$OMP THREADPRIVATE(/muflag/)

*-----------------------------------------------------------------------

      common /icomon/ no,mmat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)

      common /tarzmn/ itarz, itarm, itarn
!$OMP THREADPRIVATE(/tarzmn/)

*-----------------------------------------------------------------------

cABE 2015/06/24 added
      dimension pirot(3)
      dimension pfrot(3)
      dimension pfin(3)

      parameter ( pi = 3.141592653589793d0 )
      parameter ( rpmass = 938.272046d0 )       ! proton mass [MeV/c**2]
      parameter ( rnmass = 939.565379d0 )       ! neutron mass [MeV/c**2]
      parameter ( rmumas = 105.6583715d0 )      ! muon mass [MeV/c**2]

*-----------------------------------------------------------------------
*     initialize
*-----------------------------------------------------------------------

      if( imuint .ne. 1 ) return

      imuinthit = 1     ! S.Abe 2015/08/14, added

      nclst = 0
      nclsts = 0

      do i = 0, 20
       numpal(i) = 0
       rumpal(i) = 0.d0
      enddo

cABE 2018/06/05
      wgt = wgti
      wg0 = wgti

*-----------------------------------------------------------------------
*     set mother nucleus
*-----------------------------------------------------------------------

      itarz = idnint(ztar)
      itarm = idnint(atar)
      itarn = itarm - itarz

      be = bindeg(itarz,itarn)
      rtmass = (itarz * rpmass + itarn * rnmass - be)

*-----------------------------------------------------------------------
*     sampling energy transfer (enu) and scattering angle (themu)
*-----------------------------------------------------------------------

      emukin = eein

      call getenu(emukin,itarz,itarn,enueve)
      call getdir(emukin,itarz,itarn,enueve,dirmu,dirta,ephoto)

      emu = emukin - enueve
      phi = 2.d0 * pi * unirn(dummy)

      umuot = dsqrt(1.d0-dirmu**2.d0) * dcos(phi)
      vmuot = dsqrt(1.d0-dirmu**2.d0) * dsin(phi)
      wmuot = dirmu

      utaot = dsqrt(1.d0-dirta**2.d0) * (-dcos(phi))
      vtaot = dsqrt(1.d0-dirta**2.d0) * (-dsin(phi))
      wtaot = dirta

*-----------------------------------------------------------------------
*     bank scattered muon
*-----------------------------------------------------------------------

      nclsts = nclsts + 1
      iclusts(nclsts) = 6

      jclusts(0,nclsts) = 0
      jclusts(1,nclsts) = 0
      jclusts(2,nclsts) = 0
      jclusts(3,nclsts) = ityp
      jclusts(4,nclsts) = 0
      jclusts(5,nclsts) = jtyp
      jclusts(6,nclsts) = 0
      jclusts(7,nclsts) = ktyp
      jclusts(8,nclsts) = 0

      rms = rmumas
      rmg = rms / 1.d+3
      pr = dsqrt(emu * (emu + 2.d0 * rms)) / 1.d+3
      ett = dsqrt(pr**2.d0 + rmg**2.d0)

      pxl = pr * umuot
      pyl = pr * vmuot
      pzl = pr * wmuot

      qclusts(0,nclsts) = 0.d0
      qclusts(1,nclsts) = pxl
      qclusts(2,nclsts) = pyl
      qclusts(3,nclsts) = pzl
      qclusts(4,nclsts) = ett
      qclusts(5,nclsts) = rmg
      qclusts(6,nclsts) = 0.d0
      qclusts(7,nclsts) = emu
      qclusts(8,nclsts) = wgt / wg0
      qclusts(9,nclsts) = 0.d0
      qclusts(10,nclsts) = 0.d0
      qclusts(11,nclsts) = 0.d0
      qclusts(12,nclsts) = 0.d0

      numpal(7) = numpal(7) + 1
      rumpal(7) = rumpal(7) + wgti

*-----------------------------------------------------------------------
*     doppler shift for virtual photon in moving frame
*-----------------------------------------------------------------------

      q2eve = 2.d0 * (enueve - ephoto) * rtmass

      ptamid = dsqrt(enueve**2.d0 + q2eve) - ephoto
      vtabs0 = ptamid**2.d0 / (rtmass**2.d0 + ptamid**2.d0)
      beta0 = dsqrt(vtabs0)
      gamma0 = 1.d0 / dsqrt(1.d0 - vtabs0)

      ephsft = ephoto / gamma0 / (1.d0 + beta0)

*-----------------------------------------------------------------------
*     calculate virtual photonuclear reaction
*-----------------------------------------------------------------------

      itypb = ityp
      ktypb = ktyp
c ... setting virtual photon data in common block /icomon/
      ityp = 14
      ktyp = 22
      jtyp = ichgf(ityp,ktyp)
      mtyp = ibryf(ityp,ktyp)
      rtyp = rmtyp(ityp,ktyp)
c ... setting virtual photon data in common block /gg005/
      uuu = utaot
      vvv = vtaot
      www = wtaot

cABE change @2014/11/11
c      call sctgam(ephsft,wgti,ireg,imat)
      call sctpni(ephsft,wgti,ireg,imat)

*-----------------------------------------------------------------------
*     transform for each particles
*-----------------------------------------------------------------------

      if( nclsts .ge. 2 ) then

cABE 2015/06/24, modified
       pxb = 0.d0
       pyb = 0.d0
       pzb = ptamid * 1.d-3 * wtaot    ! [GeV]

       ptb = dsqrt(pxb**2.d0 + pyb**2.d0 + pzb**2.d0)    ! [GeV/c]
       etb = dsqrt(ptb**2.d0+(rtmass*1.d-3)**2)          ! [GeV]

       betl = ptb / etb
       gaml = 1.d0 / dsqrt(1.d0 - betl**2.d0)

c...for rotation set to direction of target middle momentum
       costhe = wtaot
       rt2 = utaot**2 + vtaot**2
       if( rt2 .eq. 0.d0 ) then
        sinthe = 0.d0
        cosphi = 1.d0
        sinphi = 0.d0
       else
        sinthe = dsqrt(rt2)
        cosphi = utaot / sinthe
        sinphi = vtaot / sinthe
        if( cosphi.lt.-1.d0 .or. cosphi.gt.1.d0 ) then
         phidum = dasin(sinphi)
         cosphi = dcos(phidum)
        endif
        if( sinphi.lt.-1.d0 .or. sinphi.gt.1.d0 ) then
         phidum = dasin(sinphi)
         cosphi = dcos(phidum)
        endif
       endif

c...for each fragment
       do i = 2, nclsts

        pirot(1) = qclusts(1,i)
        pirot(2) = qclusts(2,i)
        pirot(3) = qclusts(3,i)
        eini = qclusts(4,i)

c...composition of momentum
        pfrot(1) = pirot(1)
        pfrot(2) = pirot(2)
        pfrot(3) = pirot(3) * gaml + eini * gaml * betl

c...rotation of momentum
        pfin(1) =  pfrot(1) * costhe * cosphi
     &            -pfrot(2) * sinphi
     &           + pfrot(3) * sinthe * cosphi
        pfin(2) =  pfrot(1) * costhe * sinphi
     &           + pfrot(2) * cosphi
     &           + pfrot(3) * sinthe * sinphi
        pfin(3) = -pfrot(1) * sinthe
     &           + pfrot(3) * costhe

c.. end of transform
        qclusts(1,i) = pfin(1)
        qclusts(2,i) = pfin(2)
        qclusts(3,i) = pfin(3)

        qclusts(4,i) = dsqrt(qclusts(1,i)**2.d0 + qclusts(2,i)**2.d0
     &                       + qclusts(3,i)**2.d0 + qclusts(5,i)**2.d0)
        qclusts(7,i) = (qclusts(4,i) - qclusts(5,i)) * 1.d3

       enddo

      endif

*-----------------------------------------------------------------------
*     restore information of icomon about transported muon
*-----------------------------------------------------------------------

      ityp = itypb
      ktyp = ktypb
      jtyp = ichgf(ityp,ktyp)
      mtyp = ibryf(ityp,ktyp)
      rtyp = rmtyp(ityp,ktyp)

      imuinthit = 0

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      subroutine jgdrqd_vpi(itarz,itarn,ephoto,jgdrqd)
*                                                                      *
*        judges virtual photonuclear reaction                          *
*        almost same as function jgdrqd                                *
*                                                                      *
*        last modified by S.Abe om 2016/01/25/                         *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input  :                                                         *
*        Z : atomic number                                             *
*        N : neutron number                                            *
*        E : incident photon energy [MeV]                              *
*                                                                      *
*     output  :                                                        *
*        jgdrqd_vpi =  2 ==> NRF                                       *
*        jgdrqd_vpi =  1 ==> GDR                                       *
*        jgdrqd_vpi =  0 ==> QD                                        *
*        jgdrqd_vpi = -1 ==> Pion production                           *
*        jgdrqd_vpi = -2 ==> String production                         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)
      integer Z, N

      common /prmui/ prmui1

cABE 2016/07/13
      common /muflag2/ betaeve,ishadowflag
!$OMP THREADPRIVATE(/muflag2/)

*-----------------------------------------------------------------------

cABE 2016/07/13
c      beta = prmui1
      beta =betaeve
      atar = dble(itarz+itarn)

      Z = itarz
      N = itarn

*-----------------------------------------------------------------------

cABE 2016/07/13
c      xsnrf    = beta * sigmaNRF(1,0,ephoto)
c     &          + (1.d0-beta) * sigmaNRF(itarz,itarn,ephoto)
c      xsgdr    = beta * getGDRxsec(1,0,ephoto)
c     &          + (1.d0-beta) * getGDRxsec(itarz,itarn,ephoto)
c      xsqd     = beta * sigma_qd(1,0,ephoto)
c     &          + (1.d0-beta) * sigma_qd(itarz,itarn,ephoto)
c      xspion   = beta * sigma_delta(1.d0,ephoto,1,4)
c     &          + (1.d0-beta) * sigma_delta(atar,ephoto,1,4)
c      xsstring = beta * sigma_delta(1.d0,ephoto,5,5)
c     &          + (1.d0-beta) * sigma_delta(atar,ephoto,5,5)
      xsnrf    = sigmaNRF(itarz,itarn,ephoto)
      xsgdr    = getGDRxsec(itarz,itarn,ephoto)
      xsqd     = sigma_qd(itarz,itarn,ephoto)
      xspion   = sigma_delta(atar,ephoto,1,4)
      ishadowflag = -1
      xsstring0 = sigma_delta(atar,ephoto,5,5)
      ishadowflag = 1
      xsstring1 = sigma_delta(atar,ephoto,5,5)
      xsstring = beta * xsstring0 + (1.d0-beta) * xsstring1

      xsptot = xsnrf + xsgdr + xsqd + xspion + xsstring

      if (xsptot .le. 0.d0) then
         jgdrqd = 10
         return
      endif

*-----------------------------------------------------------------------

      ratio0 = xsnrf  / xsptot
      ratio1 = xsgdr  / xsptot + ratio0
      ratio2 = xsqd   / xsptot + ratio1
      ratio3 = xspion / xsptot + ratio2

      random = unirn()

*-----------------------------------------------------------------------

      ! 0--<NRF>--ratio0--<GDR>--ratio1--<QD>--ratio2--<PD>--1
      if ( random .le. ratio0 .and. Z .ge. 3 .and. N .ge. 3) then  ! NRF
         jgdrqd = 2
      else if ( ratio0 .lt. random .and. random .le. ratio1 ) then ! GDR
         jgdrqd = 1
         if (Z .eq. 1 .and. N .eq. 1) jgdrqd = 0
      else if ( ratio1 .lt. random .and. random .le. ratio2 ) then ! QD
         jgdrqd = 0
      elseif (ratio2 .lt. random .and. random .le. ratio3) then    ! Reso
         jgdrqd = -1
      else                                                         ! Non-reso
         jgdrqd = -2
      end if

      return
      end


************************************************************************
*                                                                      *
      subroutine stringprod
*                                                                      *
*        calculate string production by photonuclear reaction.         *
*        last modified by S.Abe on 2014/08/01                          *
*                                                                      *
*     input  :                                                         *
*        erg         : incident energy of photon (MeV)                 *
*        uuu         : u of photon                                     *
*        vvv         : v of photon                                     *
*        www         : w of photon                                     *
*                                                                      *
*        nv          : number of nucleon in the target                 *
*                                                                      *
*        k(i,nv)                                                       *
*           i =  2, KF particle code                                   *
*           i =  5, collision counter                                  *
*                                                                      *
*        p(i,nv)                                                       *
*           i =  1, x-momentum (GeV/c)                                 *
*           i =  2, y-momentum (GeV/c)                                 *
*           i =  3, z-momentum (GeV/c)                                 *
*           i =  4, total energy (GeV)                                 *
*           i =  5, mass of particle (GeV/c2)                          *
*                                                                      *
*     output :                                                         *
*        k(i,nv)                                                       *
*           i =  1, status code KS                                     *
*           i =  2, KF particle code                                   *
*           i =  5, collision counter                                  *
*           i =  7, number of collisions sufferd so far                *
*                                                                      *
*        p(i,nv)                                                       *
*           i =  1, x-momentum (GeV/c)                                 *
*           i =  2, y-momentum (GeV/c)                                 *
*           i =  3, z-momentum (GeV/c)                                 *
*           i =  4, total energy (GeV)                                 *
*           i =  5, mass of particle (GeV/c2)                          *
*                                                                      *
*        v(i,nv)                                                       *
*           i =  5, life time of the particle (fm/c)                   *
*                                                                      *
*        kq(i,nv)                                                      *
*           i =  1, flavor code of quark                               *
*           i =  2, flavor code of diquark (anti-quark)                *
*                                                                      *
*        vq(i,nv)                                                      *
*           i =  1, px (GeV/c) of parton (quark)                       *
*           i =  2, py (GeV/c) of parton (quark)                       *
*           i =  3, pz (GeV/c) of parton (quark)                       *
*           i =  4, energy (GeV) of parton (quark)                     *
*           i =  5, mass (GeV/c2) of parton (quark)                    *
*           i =  6, px (GeV/c) of parton (diquark/anti-quark)          *
*           i =  7, py (GeV/c) of parton (diquark/anti-quark)          *
*           i =  8, pz (GeV/c) of parton (diquark/anti-quark)          *
*           i =  9, energy (GeV) of parton (diquark/anti-quark)        *
*           i = 10, mass (GeV/c2) of parton (diquark/anti-quark)       *
*                                                                      *
*        pare(2)     : c.m.energy of two-body system                   *
*                                                                      *
*---- local arrays ----------------------------------------------------*
*                                                                      *
*        pph(i)      : store information of incident photon            *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*        jam1.inc                                                      *
*        jam2.inc                                                      *
*        /pnint/   : information of photon                             *
*        /gg005/   : information of photon                             *
*                                                                      *
************************************************************************

      include 'jam1.inc'
      include 'jam2.inc'
      include 'ggsparam.inc'      ! use /gg005/

*-----------------------------------------------------------------------

      dimension pph(5)

*-----------------------------------------------------------------------
*     set photon information
*-----------------------------------------------------------------------

      eph = erg / 1.d+3
      pph(1) = uuu * eph
      pph(2) = vvv * eph
      pph(3) = www * eph
      pph(4) = eph
      pph(5) = 0.d0

*-----------------------------------------------------------------------
*     Choose one nucleon at random
*-----------------------------------------------------------------------

  110 isid = idint(unirn()*dble(nv)) + 1

      kf0 = k(2,isid)
      if( kf0.eq.2212 .or. kf0.eq.2112) goto 100

      nerrcnt = nerrcnt + 1
      if( nerrcnt .gt. 100 ) then 
         write(6,*)
     $        "Error in stringprod: cannot produce string"
         call parastop( 119 )
      endif
      goto 110

  100 continue

*-----------------------------------------------------------------------
*     one nucleon absorbs photon
*-----------------------------------------------------------------------

      p(1,isid) = p(1,isid) + pph(1)
      p(2,isid) = p(2,isid) + pph(2)
      p(3,isid) = p(3,isid) + pph(3)
      p(4,isid) = p(4,isid) + pph(4)
      p(5,isid) = dsqrt(  p(4,isid)**2.d0 - p(1,isid)**2.d0
     &                  - p(2,isid)**2.d0 - p(3,isid)**2.d0 )

*-----------------------------------------------------------------------
*     c.m.energy of two-body system
*-----------------------------------------------------------------------

      pare(2) = p(5,isid)

*-----------------------------------------------------------------------
*     change information of nucleon to string
*-----------------------------------------------------------------------

c ks code
      k(1,isid) = 3
c kf code
      k(2,isid) = 92
c collision counter
      k(5,isid) = k(5,isid) + 1
c number of collision suffered so far ... set -2 so as to cal. collision with other nucleons after string decay
      k(7,isid) = -2
c life time of the particle ... string is decay immediately
      v(5,isid) = 0.d0

*-----------------------------------------------------------------------
*     define information of quark and di-quark pair
*-----------------------------------------------------------------------

      call attflv(kf0,ifla,iflb)

      qmass1 = pymass(ifla)
      qmass2 = pymass(iflb)

      emj = p(5,isid)
      pcmq = (emj**2.d0 - (qmass1+qmass2)**2.d0)
     &      * (emj**2.d0 - (qmass1-qmass2)**2.d0)

      if( pcmq .gt. 0.d0 ) then
       pcmq = dsqrt(pcmq) / (2.d0*emj)
      else
       call jamerrm(30,0,'(:) pcmq<0')
      endif
 
      q11x = 0.d0
      q11y = 0.d0
      q11z = -pcmq
      q11e = dsqrt(qmass1**2.d0 + pcmq**2.d0)
 
      q12x = 0.d0
      q12y = 0.d0
      q12z = pcmq
      q12e = dsqrt(qmass2**2.d0 + pcmq**2.d0)

c Transform to comp. frame
      the = pyangl(p(3,isid),dsqrt(p(1,isid)**2.d0+p(2,isid)**2.d0)) 
      phi = pyangl(p(1,isid),p(2,isid))
      bex = p(1,isid) / p(4,isid)
      bey = p(2,isid) / p(4,isid)
      bez = p(3,isid) / p(4,isid)
      gg1 = p(4,isid) / emj

      call jamrobo(the,phi,bex,bey,bez,gg1,q11x,q11y,q11z,q11e)
      call jamrobo(the,phi,bex,bey,bez,gg1,q12x,q12y,q12z,q12e)

c flavor code of quark
      kq(1,isid) = ifla
c momentum, energy and mass of quark
      vq(1,isid) = q11x
      vq(2,isid) = q11y
      vq(3,isid) = q11z
      vq(4,isid) = dsqrt(qmass1**2.d0+q11x**2.d0+q11y**2.d0+q11z**2.d0)
      vq(5,isid) = qmass1

c flavor code of di-/anti-quark
      kq(2,isid) = iflb
c momentum, energy and mass of di-/anti-quark
      vq(6,isid) = q12x
      vq(7,isid) = q12y
      vq(8,isid) = q12z
      vq(9,isid) = dsqrt(qmass2**2.d0+q12x**2.d0+q12y**2.d0+q12z**2.d0)
      vq(10,isid) = qmass2

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      subroutine finteg(emukin,itarz,itarn,enueve,fintq2)
*                                                                      *
*        calculate integrated function                                 *
*                                                                      *
*        last modified by S.Abe on 2014/08/04                          *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input  :                                                         *
*        emukin : incident muon energy [MeV]                           *
*        itarz  : atomic number of target                              *
*        itarn  : neutron number of target                             *
*        enueve : energy transfer [MeV]                                *
*                                                                      *
*     output :                                                         *
*        fintq2 : cross section of photonuclear interaction [barn]     *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*     emurng : energy range of calculation                             *
*        emumin : minimum transferred energy [MeV]                     *
*        emumax : maximum muon energy for CS calculation [MeV]         *
*                                                                      *
*     iq2dat : information of q2 table                                 *
*        iq2cnt : number of q2 table data                              *
*                                                                      *
*     q2data : information of q2 table                                 *
*        q2(i)  : q2 table for using dir sampling                      *
*        q2dat(i): q2 data table for using dir sampling                *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /emurng/ emumin, emumax
      common /iq2dat/ iq2cnt
!$OMP THREADPRIVATE(/iq2dat/)
      common /q2data/ q2(0:1000), q2dat(0:1000)
!$OMP THREADPRIVATE(/q2data/)

cABE 2016/07/13
      common /muflag2/ betaeve,ishadowflag
!$OMP THREADPRIVATE(/muflag2/)

*-----------------------------------------------------------------------

      dimension q2mesh(100)
      dimension width(100)

      parameter ( rpmass = 938.272046d0 )       ! proton mass [MeV/c**2]
      parameter ( rnmass = 939.565379d0 )       ! neutron mass [MeV/c**2]
      parameter ( rmumas = 105.6583715d0 )      ! muon mass [MeV/c**2]

*-----------------------------------------------------------------------
*     initialize
*-----------------------------------------------------------------------

      fintq2 = 0.d0
      if (enueve .lt. 0.5d0) return
      if (enueve .lt. emumin) return

*-----------------------------------------------------------------------

      fl = 0.d0
      fc = 0.d0
      fr = 0.d0

      ephotl = 0.d0
      ephotc = 0.d0
      ephotr = 0.d0

      phxsl = 0.d0
      phxsc = 0.d0
      phxsr = 0.d0
      flphxs = 0.d0
      fcphxs = 0.d0
      frphxs = 0.d0

      emutot = emukin + rmumas
      enu = enueve

      iq2epf = 0

*-----------------------------------------------------------------------
*     integration photon energy bin, using 2 pattern
*-----------------------------------------------------------------------

      ephmin = 0.5d0
      ephcng = 60.0d0
      q2cng = 2.d0 * rpmass * (enu - ephcng)

      iq2 = 50
      iph = 50

      iq2cnt = 0

*-----------------------------------------------------------------------
*     make ephoto mesh and stored
*-----------------------------------------------------------------------

      ddeph = dlog(ephcng/ephmin) / (dble(iph))
      ephr = ephcng
      q2r = 2.d0 * rpmass * (enu - ephr)
      q2mesh(iph+1) = q2r

      do k = iph, 1, -1

       ephr = ephmin * dexp((dble(k-1))*ddeph)
       q2l = q2r
       q2r = 2.d0 * rpmass * (enu - ephr)
       dq2 = q2r - q2l

       q2mesh(k) = q2r
       width(k) = dq2

      enddo

*-----------------------------------------------------------------------
*     search q2min by q2 mesh
*-----------------------------------------------------------------------

      iq2chk = 10000
      q2min = 1.d-12
      q2max = 2.d0 * rpmass * (enu - ephmin)

      if( q2max .le. q2min ) goto 9999

      ddq2 = dlog(q2max/q2min) / dble(iq2chk)

      do i = 0, iq2chk

       q2chk =  q2min * dexp((dble(i))*ddq2)
       fchk = fmino(emukin,enu,q2chk,itarz,itarn)

       if( fchk .ge. 0.d0 ) then
        q2min = q2chk
        ephmax = enu - q2min * 0.5d0 / rpmass
        goto 100
       endif

      enddo

      goto 9999

*-----------------------------------------------------------------------
*     calculation of fmino based on more small mesh
*-----------------------------------------------------------------------

  100 continue

      ddq2 = dlog(q2max/q2min) / (dble(iq2))

      q2r = q2min
      ephotr = ephmax
cABE 2016/03/09
c      fr = fmino(emukin,enu,q2r,itarz,itarn)
c      phxsr = phxs(itarz, itarn, ephotr)
c      frphxs = phxsr * fr
cABE 2016/07/13
c      phxs0 = phxs(itarz, itarn, ephotr)
c      phxs1 = phxs(1, 0, ephotr)
      ishadowflag = -1
      phxs0 = phxs(itarz, itarn, ephotr)
      ishadowflag = 1
      phxs1 = phxs(itarz, itarn, ephotr)
      frphxs = fmino_v2(emukin,enu,q2r,itarz,itarn,phxs0,phxs1)

      do i = 1, iq2

       q2l = q2r
       q2r = q2min * dexp((dble(i))*ddq2)
       q2c = (q2l + q2r) * 0.5d0
       dq2 = q2r - q2l

       if( q2r .gt. q2cng ) then
        q2chk = q2l
        j = iph

        if( q2chk.lt.q2cng .and. dq2.gt.width(iph)) then
         ista = iph
         iq2epf = 1

         q2r = q2mesh(iph+1)
         q2c = (q2l + q2r) *0.5d0
         dq2 = q2r - q2l
        else
         do j = iph, 1, -1
          if( q2chk.ge.q2mesh(j+1) .and. q2chk.lt.q2mesh(j) ) then
           if( dq2 .gt. width(j) ) then
            ista = j
            q2r = q2l
            goto 200
           else
            exit
           endif
          endif
         enddo
        endif
       endif

       ephotl = ephotr
       ephotc = enu - q2c * 0.5d0 / rpmass
       ephotr = enu - q2r * 0.5d0 / rpmass

cABE 2016/03/09
c       fl = fr
c       fc = fmino(emukin,enu,q2c,itarz,itarn)
c       fr = fmino(emukin,enu,q2r,itarz,itarn)

*-----------------------------------------------------------------------
*     calculation of photonuclear reaction cross section
*-----------------------------------------------------------------------

cABE 2016/03/09
c       phxsl = phxsr
c       phxsc = phxs(itarz,itarn,ephotc)
c       phxsr = phxs(itarz,itarn,ephotr)
c
c       flphxs = frphxs
c       fcphxs = phxsc * fc
c       frphxs = phxsr * fr
       flphxs = frphxs

cABE 2016/07/13
c       phxs0 = phxs(itarz, itarn, ephotc)
c       phxs1 = phxs(1, 0, ephotc)
       ishadowflag = -1
       phxs0 = phxs(itarz, itarn, ephotc)
       ishadowflag = 1
       phxs1 = phxs(itarz, itarn, ephotc)
       fcphxs = fmino_v2(emukin,enu,q2c,itarz,itarn,phxs0,phxs1)

cABE 2016/07/13
c       phxs0 = phxs(itarz, itarn, ephotr)
c       phxs1 = phxs(1, 0, ephotr)
       ishadowflag = -1
       phxs0 = phxs(itarz, itarn, ephotr)
       ishadowflag = 1
       phxs1 = phxs(itarz, itarn, ephotr)
       frphxs = fmino_v2(emukin,enu,q2r,itarz,itarn,phxs0,phxs1)

       if( flphxs.gt.0.d0 .and. frphxs.gt.0.d0 ) then
        fintq2 = fintq2
     &          + (flphxs + 4.d0*fcphxs + frphxs) / 6.d0 * dq2
        iq2cnt = iq2cnt + 1
        q2(iq2cnt) = q2r
        q2dat(iq2cnt) = fintq2
       endif

      enddo

      if( iq2epf .ne. 1 ) goto 9999

*-----------------------------------------------------------------------
*     calculation of fmino by ephoto mesh
*-----------------------------------------------------------------------

  200 continue

      do i = ista, 1, -1

       q2l = q2r
       q2r = q2mesh(i)
       q2c = (q2l + q2r) * 0.5d0
       dq2 = q2r - q2l

       ephotl = ephotr
       ephotc = enu - q2c * 0.5d0 / rpmass
       ephotr = enu - q2r * 0.5d0 / rpmass

cABE 2016/03/09
c       fl = fr
c       fc = fmino(emukin,enu,q2c,itarz,itarn)
c       fr = fmino(emukin,enu,q2r,itarz,itarn)

*-----------------------------------------------------------------------
*     calculation of photonuclear reaction cross section
*-----------------------------------------------------------------------

cABE 2016/03/09
c       phxsl = phxsr
c       phxsc = phxs(itarz, itarn, ephotc)
c       phxsr = phxs(itarz, itarn, ephotr)
c
c       flphxs = frphxs
c       fcphxs = phxsc * fc
c       frphxs = phxsr * fr
       flphxs = frphxs

cABE 2016/07/13
c       phxs0 = phxs(itarz, itarn, ephotc)
c       phxs1 = phxs(1, 0, ephotc)
       ishadowflag = -1
       phxs0 = phxs(itarz, itarn, ephotc)
       ishadowflag = 1
       phxs1 = phxs(itarz, itarn, ephotc)
       fcphxs = fmino_v2(emukin,enu,q2c,itarz,itarn,phxs0,phxs1)

cABE 2016/07/13
c       phxs0 = phxs(itarz, itarn, ephotr)
c       phxs1 = phxs(1, 0, ephotr)
       ishadowflag = -1
       phxs0 = phxs(itarz, itarn, ephotr)
       ishadowflag = 1
       phxs1 = phxs(itarz, itarn, ephotr)
       frphxs = fmino_v2(emukin,enu,q2r,itarz,itarn,phxs0,phxs1)

       if( flphxs.gt.0.d0 .and. frphxs.gt.0.d0 ) then
        fintq2 = fintq2
     &          + (flphxs + 4.d0*fcphxs + frphxs) / 6.d0 * dq2
        iq2cnt = iq2cnt + 1
        q2(iq2cnt) = q2r
        q2dat(iq2cnt) = fintq2
       endif

      enddo

*-----------------------------------------------------------------------
 9999 continue

      return
      end


************************************************************************
*                                                                      *
      subroutine getenu(emukin,itarz,itarn,enueve)
*                                                                      *
*        sampling for energy transfer                                  *
*                                                                      *
*        last modified by S.Abe on 2014/06/02                          *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input  :                                                         *
*        emukin : incident muon energy [MeV]                           *
*        itarz  : atomic number of target                              *
*        itarn  : neutron number of target                             *
*                                                                      *
*     output :                                                         *
*        enueve : sampled energy transfer [MeV]                        *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*     emurng : energy range of calculation                             *
*        emumin : minimum transferred energy [MeV]                     *
*        emumax : maximum muon energy for CS calculation [MeV]         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /emurng/ emumin, emumax

*-----------------------------------------------------------------------

      dimension inua(6)

      dimension enuloa(6)
      dimension enuhia(6)

      dimension enu(0:1000)
      dimension enudat(0:1000)

      parameter ( rmumas = 105.6583715d0 )      ! muon mass [MeV/c**2]

*-----------------------------------------------------------------------
*     initialize
*-----------------------------------------------------------------------

      enueve = 0.d0

      if( emukin .lt. 0.5d0 ) return
      if( emukin .lt. emumin ) return

*-----------------------------------------------------------------------

      enumin = emumin
      enumax = emukin

      enunow = emumin

      icntnu = 1
      xsmu = 0.d0

*-----------------------------------------------------------------------
*     integration energy bin for transferred energy 
*-----------------------------------------------------------------------

      inua(1) = 50
      enuloa(1) = 0.5d0
      enuhia(1) = 60.0d0

      inua(2) = 5
      enuloa(2) = enuhia(1)
      enuhia(2) = 200.0d0

      inua(3) = 10
      enuloa(3) = enuhia(2)
      enuhia(3) = 500.0d0

      inua(4) = 20
      enuloa(4) = enuhia(3)
      enuhia(4) = 1500.0d0

      inua(5) = 5
      enuloa(5) = enuhia(4)
      enuhia(5) = 2000.0d0

      inua(6) = 5
      enuloa(6) = enuhia(5)
      enuhia(6) = enumax

      inuflg = 0

*-----------------------------------------------------------------------
*     Search first bin for target material
*-----------------------------------------------------------------------

      inu = 100
      ddenu = dlog(enumax/enumin) / dble(inu)

      do i = 0, inu
       enuchk = enumin * dexp((dble(i))*ddenu)
       phxsck = phxs(itarz, itarn, enuchk)
       if( phxsck .gt. 0.d0 ) exit
      enddo

*-----------------------------------------------------------------------

      do i = 1, 6

       if( enuchk.ge.enuloa(i) .and. enuchk.lt.enuhia(i) ) then

        inuflg = i

        inu = inua(inuflg)
        enulo = enuloa(inuflg)
        enuhi = enuhia(inuflg)
        ddenu = dlog(enuhi/enulo) / (dble(inu))

        do j = 0, inu
         enunow = enulo * dexp((j)*ddenu)
         if( enunow .ge. enuchk ) then
          call finteg(emukin,itarz,itarn,enunow,fintq2)
          if( fintq2 .gt. 0.d0 ) then
           enu(1) = enunow
           enudat(1) = 0.d0
           ista = j + 1
           exit
          endif
         endif
        enddo
        exit
       endif
      enddo

      if( inuflg .eq. 0 ) goto 9999

*-----------------------------------------------------------------------
*     log-bin for enu [MeV]
*-----------------------------------------------------------------------

  200 continue

      do i = ista, inu

       enuold = enunow
       finold = fintq2

       enunow = enulo * dexp(dble(i)*ddenu)
       enunow = dmin1(enunow, enumax)

       denu = enunow - enuold

       call finteg(emukin,itarz,itarn,enunow,fintq2)

       if( fintq2 .gt. 0.d0 ) then
        xsmu = xsmu + (finold + fintq2) * denu * 0.5d0
        icntnu = icntnu + 1
        enu(icntnu) = enunow
        enudat(icntnu) = xsmu
       else
        goto 400
       endif

       if( enunow .ge. enumax ) goto 400

      enddo

*-----------------------------------------------------------------------

  300 continue

      inuflg = inuflg + 1
      ista = 1

      if( inuflg .le. 6 ) then
       inu = inua(inuflg)
       enulo = enuloa(inuflg)
       enuhi = enuhia(inuflg)
       ddenu = log(enuhi/enulo) / (dble(inu))

       goto 200
      endif

*-----------------------------------------------------------------------
*     sampling of enu
*-----------------------------------------------------------------------

  400 continue

      if( xsmu .le. 0.d0 ) then
       write(6,*) "Error in getenu: xsmu is zero"
       call parastop( 119 )
      endif

      enudat(0) = unirn(dummy) * xsmu

      if( enudat(0) .le. enudat(1) ) then
       enueve = enu(1)
      elseif( enudat(0) .ge. enudat(icntnu) ) then
       enueve = enu(icntnu-1)
      else
       do i = 1, icntnu-1
        if( enudat(0) .gt. enudat(i) .and.
     &      enudat(0) .le. enudat(i+1) ) then
         enueve = enu(i) * dexp(dlog(enu(i+1)/enu(i))*unirn(dummy))
         exit
        endif
       enddo
      endif

*-----------------------------------------------------------------------
 9999 continue

      return
      end


************************************************************************
*                                                                      *
      subroutine getdir(emukin,itarz,itarn, enueve,
     &                  dirmu,dirta,ephoto)
*                                                                      *
*        sampling for scattered angle for muon and target              *
*                                                                      *
*        last modified by S.Abe on 2014/05/27                          *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input  :                                                         *
*        emukin : incident muon energy [MeV]                           *
*        itarz  : target proton number                                 *
*        itarn  : target neutron number                                *
*        enueve : energy transfer [MeV]                                *
*                                                                      *
*     output :                                                         *
*        dirmu  : sampled scattaring muon direction [cos]              *
*        dirta  : sampled recoil target direction [cos]                *
*        ephoto : sampled photon energy [MeV]                          *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*     emurng : energy range of calculation                             *
*        emumin : minimum transferred energy [MeV]                     *
*        emumax : maximum muon energy for CS calculation [MeV]         *
*                                                                      *
*     iq2dat : information of q2 table                                 *
*        iq2cnt : number of q2 table data                              *
*                                                                      *
*     q2data : information of q2 table                                 *
*        q2(i)  : q2 table for using dir sampling                      *
*        q2dat(i): q2 data table for using dir sampling                *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /emurng/ emumin, emumax
      common /iq2dat/ iq2cnt
!$OMP THREADPRIVATE(/iq2dat/)
      common /q2data/ q2(0:1000), q2dat(0:1000)
!$OMP THREADPRIVATE(/q2data/)

cABE 2016/07/13
      common /prmui/ prmui1
      common /muflag2/ betaeve,ishadowflag
!$OMP THREADPRIVATE(/muflag2/)

*-----------------------------------------------------------------------

      parameter ( rpmass = 938.272046d0 )       ! proton mass [MeV/c**2]
      parameter ( rnmass = 939.565379d0 )       ! neutron mass [MeV/c**2]
      parameter ( rmumas = 105.6583715d0 )      ! muon mass [MeV/c**2]

*-----------------------------------------------------------------------
*      initialze
*-----------------------------------------------------------------------

      dirmu = 0.d0
      dirta = 0.d0
      ephoto = 0.d0

      if( enueve .lt. 0.5d0 ) return
      if( enueve .lt. emumin ) return

      enu = enueve

      iq2cnt = 0
      fintq2 = 0.d0

*-----------------------------------------------------------------------

      call finteg(emukin,itarz,itarn,enu,fintq2)

*-----------------------------------------------------------------------
*     sampling of q2
*-----------------------------------------------------------------------

      if( fintq2 .le. 0.d0 ) then
       write(6,*) "Error in getdir: fintq2 is zero"
       call parastop( 119)
      endif

      q2dat(0) = unirn(dummy) * fintq2

      if( q2dat(0) .le. q2dat(1) )then
       q2eve = q2(1)
      elseif( q2dat(0) .ge. q2dat(iq2cnt) )then
       q2eve = q2(iq2cnt)
      else
       do i = 1, iq2cnt-1
        if (q2dat(0) .gt. q2dat(i) .and.
     &      q2dat(0) .le. q2dat(i+1)) then
         q2eve = q2(i) * dexp(dlog(q2(i+1)/q2(i))*unirn(dummy))

         dirmu = fdirmu(emukin, enu, q2eve)
         dirta = fdirta(emukin, enu, q2eve)
         ephoto = enu - q2eve * 0.5d0 / rpmass

         exit
        endif
       enddo
      endif

cABE 2016/07/13
      xxxmax = prmui1
      xxxeve = q2eve / (2.d0*rpmass*enueve)
      betaeve = (1.d0 - 0.22d0) / xxxmax * xxxeve + 0.22d0

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      function fmino(emukin,enu,q2,itarz,itarn)
*                                                                      *
*        calculate eq.(12) of [1]                                      *
*                                                                      *
*        last modified by S.Abe on 2014/05/14                          *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input :                                                          *
*        emukin: kinetic energy of muon [MeV]                          *
*        enu   : energy transfer [MeV]                                 *
*        q2    : four-momentum transfer [(MeV/c)**2]                   *
*        itarz : atomic number of target                               *
*        itarn : neutron number of target                              *
*                                                                      *
*     output:                                                          *
*        fmino : solution of eq.(12)                                   *
*                                                                      *
*---- Reference -------------------------------------------------------*
*                                                                      *
*     [1] Y.Minorikawa+,Il Nuovo Cimento C,4,471(1981)                 *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      common /prmui/ prmui1

      parameter ( alpha = 7.297353d-3 )         ! fine-structure constant
      parameter ( pi = 3.141592653589793d0 )
      parameter ( rpmass = 938.272046d0 )       ! proton mass [MeV/c**2]
      parameter ( rnmass = 939.565379d0 )       ! neutron mass [MeV/c**2]
      parameter ( rmumas = 105.6583715d0 )      ! muon mass [MeV/c**2]

*-----------------------------------------------------------------------

      enu2 = enu**2.d0

*-----------------------------------------------------------------------
*     from Table. 1
*-----------------------------------------------------------------------

      beta = prmui1
      sup = 1.0d0
      mbar2 = 0.36d+6
      r = q2 / enu2

      itara = itarz + itarn
      if (q2 .le. 0.25d+6) then
         delta = dble(itara)**(-0.11d0)
      else
         delta = 1.d0
      endif

*-----------------------------------------------------------------------
*     iktflg = 0: E is kinematic energy of incident muon
*            = 1: E is total energy of incident
*-----------------------------------------------------------------------

      iktflg = 1

      if (iktflg .eq. 0) then
         p2 = emukin**2.d0 + 2.d0 * rmumas * emukin
      elseif (iktflg .eq. 1) then
         emutot = emukin + rmumas
         p2 = emutot**2.d0 - rmumas**2.d0
      endif

*-----------------------------------------------------------------------

      c1 = alpha * (enu-(q2/(2.d0*rpmass))) / (2.d0*pi*q2*p2)
      c2 = 2.d0 * rmumas**2.d0 / q2
      if( iktflg .eq. 0 ) then
         c3 = (2.d0 * emukin * (emukin-enu) - q2*0.5d0) / (q2+enu2)
      elseif( iktflg .eq. 1 ) then
         c3 = (2.d0 * emutot * (emutot-enu) - q2*0.5d0) / (q2+enu2)
      endif

      gammat = c1 * (1.d0 - c2 + c3)
      epslon = c3 / (1.d0 - c2 + c3)
      cterm  = beta + (1.d0-beta) * delta * (mbar2/(mbar2+q2))**sup

      fmino = gammat * (1.d0 + epslon * r) * cterm

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      function fmino_v2(emukin,enu,q2,itarz,itarn,phxs0,phxs1)
*                                                                      *
*        calculate eq.(12) of [1]                                      *
*                                                                      *
*        last modified by S.Abe on 2016/01/14                          *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input :                                                          *
*        emukin: kinetic energy of muon [MeV]                          *
*        enu   : energy transfer [MeV]                                 *
*        q2    : four-momentum transfer [(MeV/c)**2]                   *
*        itarz : atomic number of target                               *
*        itarn : neutron number of target                              *
*        phxs0 : photonuclear cross section of point-like              *
*        phxs1 : photonuclear cross section of hadron-like             *
*                                                                      *
*     output:                                                          *
*        fmino : solution of eq.(12)                                   *
*                                                                      *
*---- Reference -------------------------------------------------------*
*                                                                      *
*     [1] Y.Minorikawa+,Il Nuovo Cimento C,4,471(1981)                 *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      common /prmui/ prmui1

      parameter ( alpha = 7.297353d-3 )         ! fine-structure constant
      parameter ( pi = 3.141592653589793d0 )
      parameter ( rpmass = 938.272046d0 )       ! proton mass [MeV/c**2]
      parameter ( rnmass = 939.565379d0 )       ! neutron mass [MeV/c**2]
      parameter ( rmumas = 105.6583715d0 )      ! muon mass [MeV/c**2]

*-----------------------------------------------------------------------

      enu2 = enu**2.d0

*-----------------------------------------------------------------------
*     from Table 1
*-----------------------------------------------------------------------

cABE 2016/07/13
c      beta = prmui1
      ephot = enu - q2 * 0.5d0 / rpmass
      xxx = q2 / (2.d0*rpmass*enu)
      xxxmax = prmui1
      beta = (1.d0 - 0.22d0) / xxxmax * xxx + 0.22d0
      if( beta .gt. 1.d0 ) beta = 1.d0

      sup = 1.0d0
      mbar2 = 0.36d+6   ![MeV^2]
      r = q2 / enu2

      itara = itarz + itarn

*-----------------------------------------------------------------------
*     iktflg = 0: E is kinematic energy of incident muon
*            = 1: E is total energy of incident
*-----------------------------------------------------------------------

      iktflg = 1

      if (iktflg .eq. 0) then
         p2 = emukin**2.d0 + 2.d0 * rmumas * emukin
      elseif (iktflg .eq. 1) then
         emutot = emukin + rmumas
         p2 = emutot**2.d0 - rmumas**2.d0
      endif

*-----------------------------------------------------------------------
*     eq.(12)
*-----------------------------------------------------------------------

cABE 2016/07/13
c      sigpnt = beta * dble(itara) * phxs1
c      sighad = (1.d0-beta) * phxs0 * (mbar2/(mbar2+q2))**sup
      sigpnt = beta * phxs0
      sighad = (1.d0-beta) * phxs1 * (mbar2/(mbar2+q2))**sup

      cterm  = sigpnt + sighad

*-----------------------------------------------------------------------
*     eq.(4)
*-----------------------------------------------------------------------

      c1 = alpha * (enu-(q2/(2.d0*rpmass))) / (2.d0*pi*q2*p2)
      c2 = 2.d0 * rmumas**2.d0 / q2
      if( iktflg .eq. 0 ) then
         c3 = (2.d0 * emukin * (emukin-enu) - q2*0.5d0) / (q2+enu2)
      elseif( iktflg .eq. 1 ) then
         c3 = (2.d0 * emutot * (emutot-enu) - q2*0.5d0) / (q2+enu2)
      endif

      gammat = c1 * (1.d0 - c2 + c3)
      epslon = c3 / (1.d0 - c2 + c3)

      fmino_v2 = gammat * (1.d0 + epslon * r) * cterm

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      function fdirmu(emukin,enu,q2)
*                                                                      *
*        calculation of scattered muon angle                           *
*                                                                      *
*        last modified by S.Abe on 2014/05/14                          *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input :                                                          *
*        emukin: kinetic energy of muon [MeV]                          *
*        enu   : energy transfer [MeV]                                 *
*        q2    : four-momentum transfer [(MeV/c)**2]                   *
*                                                                      *
*     output:                                                          *
*        fdirmu: angle of scatterd muon [cos]                          *
*                                                                      *
*---- Reference -------------------------------------------------------*
*                                                                      *
*     [1] Y.Minorikawa+,Il Nuovo Cimento C,4,471(1981)                 *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      parameter ( pi = 3.141592653589793d0 )
      parameter ( rmumas = 105.6583715d0 )      ! muon mass [MeV/c**2]

*-----------------------------------------------------------------------

      emutot = emukin + rmumas

*-----------------------------------------------------------------------

      c1 = emutot * (emutot-enu) - rmumas**2.d0 - q2 * 0.5d0
      c2 = dsqrt(emutot**2.d0 - rmumas**2.d0)
      c3 = dsqrt((emutot-enu)**2.d0 - rmumas**2.d0)

      cosmu = c1 / (c2 * c3)
      if( cosmu .le. -1.d0 ) cosmu = -1.d0
      if( cosmu .ge.  1.d0 ) cosmu =  1.d0

      fdirmu = cosmu

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      function fdirta(emukin,enu,q2)
*                                                                      *
*        calculation of recoil target angle                            *
*                                                                      *
*        last modified by S.Abe on 2014/05/14                          *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input :                                                          *
*        emukin: kinetic energy of muon [MeV]                          *
*        enu   : energy transfer [MeV]                                 *
*        q2    : four-momentum transfer [(MeV/c)**2]                   *
*                                                                      *
*     output:                                                          *
*        fdirta: angle of emitted target [cos]                         *
*                                                                      *
*---- Reference -------------------------------------------------------*
*                                                                      *
*     [1] Y.Minorikawa+,Il Nuovo Cimento C,4,471(1981)                 *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      parameter ( pi = 3.141592653589793d0 )
      parameter ( rmumas = 105.6583715d0 )      ! muon mass [MeV/c**2]

*-----------------------------------------------------------------------

      emutot = emukin + rmumas

*-----------------------------------------------------------------------

      c1 = emutot * enu + 0.5d0 * q2
      c2 = dsqrt(enu**2.d0 + q2)
      c3 = dsqrt(emutot**2.d0 - rmumas**2.d0)

      costa = c1 / (c2 * c3)
      if( costa .le. -1.d0 ) costa = -1.d0
      if( costa .ge.  1.d0 ) costa = 1.d0

      fdirta = costa

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      subroutine sigmubrm(eein,itz,ita,sigbt)
*                                                                      *
*        calculate cross section of muon-induced bremsstrahlung        *
*                                                                      *
*        last modified by S.Abe on 2015/06/03                          *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input :                                                          *
*        eein  : kinetic energy of muon [MeV]                          *
*        itz   : atomic number of target                               *
*        ita   : mass number of target                                 *
*                                                                      *
*     output:                                                          *
*        sigbt : cross section of muon-induced bremsstrahlung [barn]   *
*                                                                      *
*---- Reference -------------------------------------------------------*
*                                                                      *
*     [0] D.E.Groom+,Atom.Data Nucl.Data Tables,76,1(2001)             *
*     [1] S.R.Kelner+,preprint MSEPI 024-95,Moscow(1995)               *
*     [2] S.R.Kelner+,Phys.Atomic Nuclei,60,657(1997)                  *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      parameter ( remass = 0.510999d0 )
      parameter ( rmumas = 105.6583715d0 )      ! muon mass [MeV/c**2]

*-----------------------------------------------------------------------
*     initialize
*-----------------------------------------------------------------------

      sigbn = 0.d0
      sigbe = 0.d0
      sigbt = 0.d0

c ... for utlra relativistic energy
      if( eein .lt. 10.d+3 ) return

*-----------------------------------------------------------------------

      emutot = eein + rmumas

      etrfmin = 1.d-5
      etrfmax = eein / emutot
      netrf = 60
      detrf = dlog(etrfmax/etrfmin) / dble(netrf)

      etrfmaxe = 1.d0 / (1.d0 + rmumas**2.d0/(2.d0*remass*emutot))

*-----------------------------------------------------------------------
*     calculate differential cross section for transferred energy
*-----------------------------------------------------------------------

      etrf = etrfmin
      dsigbn = dmax1(0.d0,fmubrn(itz,ita,eein,etrf))
      dsigbe = dmax1(0.d0,fmubre(itz,ita,eein,etrf))
      dsigbt = dsigbn + dsigbe

      do i = 1, netrf-1

       etrf0 = etrf
       dsigbn0 = dsigbn
       dsigbe0 = dsigbe
       dsigbt0 = dsigbt

       etrf = etrfmin * dexp((dble(i))*detrf)

       dsigbn = dmax1(0.d0,fmubrn(itz,ita,eein,etrf))
       if( etrf .le. etrfmaxe ) then
        dsigbe = dmax1(0.d0,fmubre(itz,ita,eein,etrf))
       else
        dsigbe = 0.d0
       endif
       dsigbt = dsigbn + dsigbe

*-----------------------------------------------------------------------
*     integrate total cross section of bremsstrahlung
*-----------------------------------------------------------------------

       if( dsigbt .gt. 0.d0 ) then
        sigbt = sigbt + 0.5d0 * (dsigbt0+dsigbt) * (etrf-etrf0) * 1d+24 ! [cm2 to barn]
       else
        exit
       endif

      enddo

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      subroutine brmmu(eein,wgti,ztar,atar)
*                                                                      *
*        calculate bremsstrahlung for muon                             *
*                                                                      *
*        last modified by S.Abe on 2015/07/16                          *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input  :                                                         *
*        eein   : kinetic energy of muon [MeV]                         *
*        wgti   : weight for incident muon                             *
*        ztar   : atomic number of target                              *
*        atar   : mass number of target                                *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*     muint  :                                                         *
*        imubrm : switch of muon-induced bremsstrahlung                *
*                                                                      *
*     muflag :                                                         *
*        imubrmhit : flag for muon-induced bremsstrahlung              *
*                                                                      *
*     icomon :                                                         *
*        ityp   : particle type for incident particle                  *
*        ktyp   : kf-code for incident particle                        *
*        jtyp   : charge for incident particle                         *
*                                                                      *
*---- Reference -------------------------------------------------------*
*                                                                      *
*     [1] A.V.Ginneken,Nucl.Instr.Meth.A,251,21,(1986)                 *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param00.inc'

*-----------------------------------------------------------------------

      common /muint/ imuint, imubrm, imuppd, imucap
      common /muflag/ imuinthit,imubrmhit,imuppdhit,imucaphit,imucapflag
!$OMP THREADPRIVATE(/muflag/)

*-----------------------------------------------------------------------

      common /icomon/ no,mmat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)

*-----------------------------------------------------------------------
cABE 2017/02/08
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)

*-----------------------------------------------------------------------

      dimension etrf(1000)
      dimension sigbn(1000)
      dimension sigbe(1000)
      dimension sigbt(1000)

      parameter ( pi = 3.141592653589793d0 )
      parameter ( remass = 0.510999d0 )         ! electron mass [MeV/c**2]
      parameter ( rmumas = 105.6583715d0 )      ! muon mass [MeV/c**2]

*-----------------------------------------------------------------------
*     initialize
*-----------------------------------------------------------------------

      if (imubrm .ne. 1) return

      imubrmhit = 1

      do i = 1, 1000
       sigbn(i) = 0.d0
       sigbe(i) = 0.d0
       sigbt(i) = 0.d0
      enddo

      nclst = 0
      nclsts = 0

*-----------------------------------------------------------------------

      emutot = eein + rmumas

      wgt = wgti
      wg0 = wgti

      do i = 0, 20
       numpal(i) = 0
       rumpal(i) = 0.0d0
      enddo

      itz = idnint(ztar)
      ita = idnint(atar)
      itn = ita - itz

*-----------------------------------------------------------------------
*     1. fraction of transferred energy
*     1-1. set range and width for fraction of transferred energy
*-----------------------------------------------------------------------

      etrfmin = 1.d-5
      etrfmax = eein / emutot
      netrf = 60
      detrf = dlog(etrfmax/etrfmin) / dble(netrf)

      etrfmaxe = 1.d0 / (1.d0 + rmumas**2.d0/(2.d0*remass*emutot))

*-----------------------------------------------------------------------
*     1-2. calculate differential cross section for etrf
*-----------------------------------------------------------------------

      etrf(1) = etrfmin
      dsigbn = dmax1(0.d0,fmubrn(itz,ita,eein,etrf(1)))
      dsigbe = dmax1(0.d0,fmubre(itz,ita,eein,etrf(1)))
      dsigbt = dsigbn + dsigbe

cABE 2023/08/18
c      do i = 1, netrf-1
      do i = 1, netrf

       dsigbn0 = dsigbn
       dsigbe0 = dsigbe
       dsigbt0 = dsigbt

       etrf(i+1) = etrfmin * dexp((dble(i))*detrf)

       dsigbn = dmax1(0.d0,fmubrn(itz,ita,eein,etrf(i+1)))
       if( etrf(i+1) .le. etrfmaxe )then
        dsigbe = dmax1(0.d0,fmubre(itz,ita,eein,etrf(i+1)))
       else
        dsigbe = 0.d0
       endif
       dsigbt = dsigbn + dsigbe

       if( dsigbt .gt. 0.d0 ) then
        sigbn(i+1) = sigbn(i) + 0.5d0 * (dsigbn0+dsigbn)
     &                                * (etrf(i+1)-etrf(i)) * 1d+24 ! [cm2 to barn]
        sigbe(i+1) = sigbe(i) + 0.5d0 * (dsigbe0+dsigbe)
     &                                * (etrf(i+1)-etrf(i)) * 1d+24 ! [cm2 to barn]
        sigbt(i+1) = sigbt(i) + 0.5d0 * (dsigbt0+dsigbt)
     &                                * (etrf(i+1)-etrf(i)) * 1d+24 ! [cm2 to barn]
        nsig = i+1
       else
        if (sigbt(i) .le. 0.0d0) then
         write(6,*) "Error in brmmu: sigmubrm is zero"
         call parastop( 119 )
        endif
        exit
       endif

      enddo

*-----------------------------------------------------------------------
*     1-3. sampling fraction of transferred energy
*-----------------------------------------------------------------------

      sigset = unirn(dummy) * sigbt(nsig)

      do i = 1, nsig-1
       if( sigset.gt.sigbt(i) .and. sigset.le.sigbt(i+1) ) then
        etrfset = etrf(i) + (etrf(i+1)-etrf(i)) * (sigbt(i+1)-sigset)
     &                                          / (sigbt(i+1)-sigbt(i))
        exit
       endif
      enddo

      eph = eein * etrfset
      emukin = eein - eph

      prph = eph
      prmu = dsqrt(emukin * (emukin + 2.d0 * rmumas))

*-----------------------------------------------------------------------
*     2. angle sampling for scattering muon and bremsstrahlung photon
*        iselectangle = 1: according to EGS5 routine, egs5.f_4904-4931
*                     = 2: using parameterization, eq.(6) of [1]
*-----------------------------------------------------------------------

      iselectangle = 2

*-----------------------------------------------------------------------
*     2-a. sampling emission angle for bremsstrahlung photon
*          according to EGS5 routine, egs5.f_4904-4931
*-----------------------------------------------------------------------

      if( iselectangle .eq. 1 ) then

       ztarg = (ztar**(1.d0/3.d0)/111.d0)**2.d0

       tteie = eein/rmumas
       ttese = emukin/rmumas
       esedei = ttese/tteie
       y2max = (pi*tteie)**2.d0

       rjarg1 = 1.d0 + esedei**2.d0
       rjarg2 = 3.d0*rjarg1 - 2.d0*esedei
       rjarg3 = ((1.d0 - esedei)/(2.d0*tteie*esedei))**2.d0

       y2tst1 = (1.d0 + 0.d0)**2.d0
       rejmin = (4.d0 + dlog(rjarg3 + ztarg/y2tst1))*
     &          (4.d0*esedei*0.d0/y2tst1 - rjarg1) + rjarg2

       y2tst1 = (1.d0 + 1.d0)**2.d0
       rejmid = (4.d0 + dlog(rjarg3 + ztarg/y2tst1))*
     &          (4.d0*esedei*1.d0/y2tst1 - rjarg1) + rjarg2

       y2tst1 = (1.d0 + y2max)**2.d0
       rejmax = (4.d0 + dlog(rjarg3 + ztarg/y2tst1))*
     &          (4.d0*esedei*y2max/y2tst1 - rjarg1) + rjarg2

       rejtop = dmax1(rejmin,rejmid,rejmax)

    5  continue
       rnnow = unirn(dummy)
       y2tst =rnnow/(1.d0 - rnnow + 1.d0/y2max)

       y2tst1 = (1.d0 + y2tst)**2.d0
       rejtst = (4.d0 + dlog(rjarg3 + ztar/y2tst1))*
     &          (4.d0*esedei*y2tst/y2tst1 - rjarg1) + rjarg2
       rnnow = unirn(dummy)
       if( rnnow .gt. (rejtst/rejtop) ) goto 5      ! loop to define y2tst

       theph = dsqrt(y2tst)/tteie

       phiph = 2.d0 * pi * unirn(dummy)   ! S.Abe 2015/12/04, corrected

*-----------------------------------------------------------------------

       pxlph = prph * dsin(theph) * dcos(phiph)
       pylph = prph * dsin(theph) * dsin(phiph)
       pzlph = prph * dcos(theph)

cABE 2015/12/04
c       pxlmu = -pxlph
c       pylmu = -pylph
c       pzlmu = dsqrt(prmu**2.d0 - pxlmu**2.d0 - pylmu**2.d0)
       contmu = dsqrt(pxlph**2.d0 + pylph**2.d0)

       if (prmu .gt. contmu) then
        pxlmu = -pxlph
        pylmu = -pylph
        pzlmu = dsqrt(prmu**2.d0 - pxlmu**2.d0 - pylmu**2.d0)
       else
        pxlmu = -pxlph * prmu / contmu
        pylmu = -pylph * prmu / contmu
        pzlmu = 0.d0
       endif

*-----------------------------------------------------------------------
*     2-b. sampling emission angle for scattering muon
*          using parameterization, eq.(6) of [1]
*-----------------------------------------------------------------------

      elseif( iselectangle .eq. 2 ) then

       emutotg = emutot*1.d-3

       if( etrfset .le. 0.5d0 ) then

        rk1 = 0.092d0 * emutotg**(-1.d0/3.d0)
        rk2 = 0.052d0 * emutotg**(-1.d0) * ztar**(-0.25)
        rk3 = 0.220d0 * emutotg**(-0.92d0)
        rmsthemu = dmax1(dmin1(rk1*etrfset**0.5d0,rk2),rk3*etrfset)

       else

        rk4 = 0.26d0 * emutotg**(-0.91d0)
        rn4 = 0.81d0 * emutotg**0.5d0 / (emutotg**0.5d0 + 1.8d0)
        rmsthemu = rk4 * etrfset**(1.d0+rn4) * (1.d0-etrfset)**(-rn4)

cABE 2023/08/18
c        if( rmsthemu .ge. 0.2d0 ) then
c! binary search for etrfdum
c         etrfr = 0.5d0
c         etrfl = etrfset
c         drms = 1.d-6
c         etrfdum = (0.5d0+etrfset) * 0.5d0      
c
c         rmsthemu = rk4 * etrfr**rn4 * (1.d0-etrfr)**(-rn4)
c         if( rmsthemu .ge. 0.2d0 ) then
c          write(6,*) "Error in brmmu: unsuited sampling emission angle"
c          write(6,*) "eein = ", eein
c          write(6,*) "target atomic number =", ztar
c          write(6,*) "target mass number =", atar
c          call parastop( 119 )
c         endif
c
c  100    continue
c         rmsthemu = rk4 * etrfdum**rn4 * (1.d0-etrfdum)**(-rn4)
c         if( rmsthemu.le.0.2d0*(1.d0+drms) .and.
c     &       rmsthemu.ge.0.2d0*(1.d0-drms) ) goto 200
c         if( rmsthemu .gt. 0.2d0 ) then
c          etrfl = etrfdum
c          etrfdum = (etrfr+etrfdum) * 0.5d0
c         endif
c         if( rmsthemu .lt. 0.2d0 ) then
c          etrfr = etrfdum
c          etrfdum = (etrfdum+etrfl) * 0.5d0
c         endif
c         goto 100
c
c  200    continue
c
c         rk5 = 0.2d0 * (1.d0-etrfdum)**0.5d0
c         rmsthemu = rk5 * (1.d0-etrfset)**0.5d0
c
c        endif

        if( rmsthemu .ge. 0.2d0 ) then

         etrfl = 0.5d0
         rmsthemu = rk4 * etrfl**(1.d0+rn4) * (1.d0-etrfl)**(-rn4)

         if( rmsthemu .lt. 0.2d0 ) then      ! binary search for etrfdum

          etrfr = etrfset
          drms = 1.d-6
          etrfdum = (etrfl+etrfr) * 0.5d0

  100     continue
          rmsthemu = rk4 * etrfdum**(1.d0+rn4) * (1.d0-etrfdum)**(-rn4)
          if( rmsthemu .le. 0.2d0*(1.d0+drms) .and.
     &        rmsthemu .ge. 0.2d0*(1.d0-drms) ) goto 200
          if( rmsthemu .gt. 0.2d0 ) then
           etrfr = etrfdum
          elseif( rmsthemu .lt. 0.2d0 ) then
           etrfl = etrfdum
          endif
          etrfdum = (etrfl+etrfr) * 0.5d0
          goto 100
  200     continue

          rk5 = 0.2d0 * (1.d0-etrfdum)**0.5d0

         elseif( rmsthemu .ge. 0.2d0 ) then      ! k5 is connected to k1-k3

          rk1 = 0.092d0 * emutotg**(-1.d0/3.d0)
          rk2 = 0.052d0 * emutotg**(-1.d0) * ztar**(-0.25)
          rk3 = 0.220d0 * emutotg**(-0.92d0)
          rmsthemu = dmax1(dmin1(rk1*etrfl**0.5d0, rk2), rk3*etrfl)

          rk5 = rmsthemu * (1.d0-etrfl)**0.5d0

         endif

         rmsthemu = rk5 * (1.d0-etrfset)**0.5d0

        endif

       endif

*-----------------------------------------------------------------------
! sampling from gauss dist with zero mean and sig**2=0.5*rmsthemu**2

       themumin = 0.d0
       themumax = 5.d0 * dsqrt(rmsthemu**2.d0/2.d0)      ! range: 5*sig
       nthemu = 1000
       dthemu = (themumax-themumin)/dble(nthemu)

       sigset = unirn(dummy) * 0.5d0

       themu = themumin
       gthemu = 1.d0 / dsqrt(pi*rmsthemu**2.d0) * dthemu
       dgthemu = 0.d0

       do i = 1, nthemu

        themu0 = themu
        gthemu0 = gthemu
        dgthemu0 = dgthemu

        themu = dthemu * dble(i) + themumin
        gthemu = dexp(-(themu**2.d0/rmsthemu**2.d0))
     &          / dsqrt(pi*rmsthemu**2.d0) * dthemu
        dgthemu = dgthemu + (gthemu+gthemu0) * (themu-themu0) * 0.5d0

        if( dgthemu .ge. sigset ) then
         themu = themu0 + (themu-themu0) * (dgthemu-sigset)
     &                                   / (dgthemu-dgthemu0)
         exit
        endif

       enddo

*-----------------------------------------------------------------------

       phimu = 2.d0 * pi * unirn(dummy)   ! S.Abe 2015/12/04, corrected

       pxlmu = prmu * dsin(themu) * dcos(phimu)
       pylmu = prmu * dsin(themu) * dsin(phimu)
       pzlmu = prmu * dcos(themu)

cABE 2015/12/04
c       pxlph = -pxlmu
c       pylph = -pylmu
c       pzlph = dsqrt(prph**2.d0 - pxlph**2.d0 - pylph**2.d0)
       contph = dsqrt(pxlmu**2.d0 + pylmu**2.d0)

       if (prph .gt. contph) then
        pxlph = -pxlmu
        pylph = -pylmu
        pzlph = dsqrt(prph**2.d0 - pxlph**2.d0 - pylph**2.d0)
       else
        pxlph = -pxlmu * prph / contph
        pylph = -pylmu * prph / contph
        pzlph = 0.d0
       endif

      endif

*-----------------------------------------------------------------------
*     bank scattering muon
*-----------------------------------------------------------------------

      nclsts = nclsts + 1
      iclusts(nclsts) = 6

      jclusts(0,nclsts) = 0
      jclusts(1,nclsts) = 0
      jclusts(2,nclsts) = 0
      jclusts(3,nclsts) = ityp
      jclusts(4,nclsts) = 0
      jclusts(5,nclsts) = jtyp
      jclusts(6,nclsts) = 0
      jclusts(7,nclsts) = ktyp
      jclusts(8,nclsts) = 0

      ett = dsqrt(prmu**2.d0 + rmumas**2.d0)

      qclusts(0,nclsts) = 0.d0
      qclusts(1,nclsts) = pxlmu * 1.d-3
      qclusts(2,nclsts) = pylmu * 1.d-3
      qclusts(3,nclsts) = pzlmu * 1.d-3
      qclusts(4,nclsts) = ett * 1.d-3
      qclusts(5,nclsts) = rmumas * 1.d-3
      qclusts(6,nclsts) = 0.d0
      qclusts(7,nclsts) = emukin
      qclusts(8,nclsts) = wgt / wg0
      qclusts(9,nclsts) = 0.d0
      qclusts(10,nclsts) = 0.d0
      qclusts(11,nclsts) = 0.d0
      qclusts(12,nclsts) = 0.d0

      numpal(7) = numpal(7) + 1
      rumpal(7) = rumpal(7) + wgti

*-----------------------------------------------------------------------
*     bank bremsstrahlung photon
*-----------------------------------------------------------------------

      nclsts = nclsts + 1
      iclusts(nclsts) = 4

      jclusts(0,nclsts) = 0
      jclusts(1,nclsts) = 0
      jclusts(2,nclsts) = 0
      jclusts(3,nclsts) = 14
      jclusts(4,nclsts) = 0
      jclusts(5,nclsts) = 0
      jclusts(6,nclsts) = 0
      jclusts(7,nclsts) = 22
      jclusts(8,nclsts) = 0

      qclusts(0,nclsts) = 0.d0
      qclusts(1,nclsts) = pxlph * 1.d-3
      qclusts(2,nclsts) = pylph * 1.d-3
      qclusts(3,nclsts) = pzlph * 1.d-3
      qclusts(4,nclsts) = prph * 1.d-3
      qclusts(5,nclsts) = 0.d0
      qclusts(6,nclsts) = 0.d0
      qclusts(7,nclsts) = eph * 1.d-3
      qclusts(8,nclsts) = wgt / wg0
      qclusts(9,nclsts) = 0.d0
      qclusts(10,nclsts) = 0.d0
      qclusts(11,nclsts) = 0.d0
      qclusts(12,nclsts) = 0.d0

      numpal(14) = numpal(14) + 1
      rumpal(14) = rumpal(14) + wgti

cABE 2017/02/08
      if( ityp .eq. 7 ) atmrc(4,3) = atmrc(4,3) + wgt / wg0
      if( ityp .eq. 6 ) atmrc(5,3) = atmrc(5,3) + wgt / wg0

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      function fmubrn(itz,ita,eein,etrf)
*                                                                      *
*        calculate cross section of muon bremsstrahlung                *
*        on the screened nucleus                                       *
*                                                                      *
*        made by S.Abe on 2015/04/28                                   *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input :                                                          *
*        itz   : atomic number of target                               *
*        ita   : mass number of target                                 *
*        eein  : kinetic energy of muon [MeV]                          *
*        etrf  : flaction of muon's energy transferred to photon [-]   *
*                                                                      *
*     output:                                                          *
*        fmubrn: differential cross section of muon brems. [cm2]       *
*                                                                      *
*---- Reference -------------------------------------------------------*
*                                                                      *
*     [1] S.R.Kelner et al.,Preprint/MEPhI 024-95,Moscow(1995)         *
*     [2] D.E.Groom et al.,Atom.Data Nucl.Data Tables,76,1(2001)       *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      parameter ( pi = 3.141592653589793d0 )
      parameter ( eneper = 2.7182818d0 )        ! Napier's constant
      parameter ( alpha = 7.297353d-3 )         ! fine-structure constant
      parameter ( cre0 = 2.8179380d-13 )        ! classcal electron radius [cm]
      parameter ( remass = 0.510998928d0 )      ! electron mass [MeV/c**2]
      parameter ( rmumas = 105.6583715d0 )      ! muon mass [MeV/c**2]

*-----------------------------------------------------------------------
*     initialize
*-----------------------------------------------------------------------

      zt = dble(itz)
      at = dble(ita)

      emutot = eein + rmumas

*-----------------------------------------------------------------------

      dn = 1.54d0 * at**0.27d0                                          ! eq.(18) of [1]

      if( itz .eq. 1 ) then
       bnucl = 202.4d0
      else
       bnucl = 182.7d0
      endif

      zt3r = zt**(-1.d0/3.d0)
      enp2r = dsqrt(eneper)

*-----------------------------------------------------------------------
*     cross section of muon bremsstrahlung on the screened nucleus
*-----------------------------------------------------------------------

      delta = rmumas**2.d0 * etrf * 0.5d0 / (emutot*(1.d0-etrf))        ! eq.(7) of [1]

      phin = dlog(bnucl*rmumas*zt3r / (remass+delta*enp2r*bnucl*zt3r))  ! eq.(20) of [1]
     &      -dlog(dn*rmumas / (rmumas+delta*(dn*enp2r-2.d0)))           ! eq.(18) of [1]

      aaa = alpha * (2.d0 * zt * remass/rmumas * cre0)**2.d0
      bbb = 4.d0/3.d0 - 4.d0/3.d0*etrf + etrf**2.d0
      fmubrn = aaa * bbb * phin / etrf                                  ! eq.(19) of [1]

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      function fmubre(itz,ita,eein,etrf)
*                                                                      *
*        calculate cross section of muon bremsstrahlung                *
*        on atomic electrons                                           *
*                                                                      *
*        made by S.Abe on 2015/04/28                                   *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input :                                                          *
*        itz   : atomic number of target                               *
*        ita   : mass number of target                                 *
*        eein  : kinetic energy of muon [MeV]                          *
*        etrf  : flaction of muon's energy transferred to photon [-]   *
*                                                                      *
*     output:                                                          *
*        fmubre: differential cross section for muon brems. [cm2]      *
*                                                                      *
*---- Reference -------------------------------------------------------*
*                                                                      *
*     [1] S.R.Kelner+,Phys.Atomic Nuclei,60,657(1997)                  *
*     [2] D.E.Groom+,Atom.Data Nucl.Data Tables,76,1(2001)             *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      parameter ( pi = 3.141592653589793d0 )
      parameter ( eneper = 2.7182818d0 )        ! Napier's constant
      parameter ( alpha = 7.297353d-3 )         ! fine-structure constant
      parameter ( cre0 = 2.8179380d-13 )        ! classcal electron radius [cm]
      parameter ( remass = 0.510998928d0 )      ! electron mass [MeV/c**2]
      parameter ( rmumas = 105.6583715d0 )      ! muon mass [MeV/c**2]

*-----------------------------------------------------------------------
*     initialize
*-----------------------------------------------------------------------

      zt = dble(itz)
      at = dble(ita)

      emutot = eein + rmumas

*-----------------------------------------------------------------------

      if( itz .eq. 1 ) then
       belec = 446.d0
      else
       belec = 1429.d0
      endif

      zt23r = zt**(-2.d0/3.d0)
      enp2r = dsqrt(eneper)

*-----------------------------------------------------------------------
*     cross section of muon bremsstrahlung on atomic electrons
*-----------------------------------------------------------------------

      delta = rmumas**2.d0 * etrf * 0.5d0 / (emutot*(1.d0-etrf))

      phie = dlog(rmumas/delta/(rmumas*delta*remass**(-2.d0)+enp2r))    ! eq.(39) of [1]
     &      -dlog(1.d0 + remass/(delta*belec*zt23r*enp2r))

      aaa = alpha * zt * (2.d0 * remass/rmumas * cre0)**2.d0
      bbb = 4.d0/3.d0 - 4.d0/3.d0*etrf + etrf**2.d0
      fmubre = aaa * bbb * phie / etrf                                  ! eq.(39) of [1]

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      subroutine sigmuppd(eein,itz,ita,sigpt)
*                                                                      *
*        calculate cross section of muon-induced pair production       *
*                                                                      *
*        made by S.Abe on 2015/06/10                                   *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input :                                                          *
*        eein  : kinetic energy of muon [MeV]                          *
*        itz   : atomic number of target                               *
*        ita   : mass number of target                                 *
*                                                                      *
*     output:                                                          *
*        sigpt : cross section of muon-induced pair production [barn]  *
*                                                                      *
*---- Reference -------------------------------------------------------*
*                                                                      *
*     [1] R.P.Kokoulin+,Proc.12th ICRC,vol.6,p.2436(1971)              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      parameter ( eneper = 2.7182818d0 )        ! Napier's constant
      parameter ( remass = 0.510999d0 )         ! electron mass [MeV/c**2]
      parameter ( rmumas = 105.6583715d0 )      ! muon mass [MeV/c**2]

*-----------------------------------------------------------------------
*     initialize
*-----------------------------------------------------------------------

      zt = dble(itz)

      sigpt = 0.d0

      emutot = eein + rmumas

      zt3r = zt**(-1.d0/3.d0)
      enp2r = dsqrt(eneper)

*-----------------------------------------------------------------------

      rchk =  4.d0*remass + 3.d0 * enp2r * rmumas / zt3r / 4.d0
      if( emutot .le. rchk ) return

*-----------------------------------------------------------------------
*     set range and width for transferred energy
*-----------------------------------------------------------------------

      etrfmin = 4.d0 * remass / emutot
      etrfmin = dmax1(1.d-5,etrfmin)

      etrfmax = 1.d0 - 3.d0 * enp2r * rmumas / zt3r / (4.d0 * emutot)
      etrfmax = dmin1(1.d0-1.d-6,etrfmax)

      netrf = 100
      detrf = dlog(etrfmax/etrfmin) / dble(netrf)

*-----------------------------------------------------------------------

      asymmin = 1.d-6
      nasym = 100

*-----------------------------------------------------------------------
*     calculate differential cross section for transferred energy
*-----------------------------------------------------------------------

      etrf = 0.d0
      dsigpt = 0.d0

      do i = 1, netrf

       etrf0 = etrf
       dsigpt0 = dsigpt

       etrf = etrfmin * dexp((dble(i))*detrf)
       dsigpt = 0.d0

*-----------------------------------------------------------------------

       asymmax = (1.d0 - 6.d0*rmumas**2.d0/emutot**2.d0/(1.d0-etrf))
     &          * dsqrt(1.d0-4.d0*remass/emutot/etrf)
       if( asymmax .le. asymmin ) exit

       dasym = dlog(asymmax/asymmin) / dble(nasym)

       asym = 0.d0
       ddsigpt = dmax1(0.d0,fmupp(itz,eein,etrf,asym))

       do j = 0, nasym

        asym0 = asym
        ddsigpt0 = ddsigpt

        asym = asymmin * dexp((dble(j))*dasym)
        ddsigpt = dmax1(0.d0,fmupp(itz,eein,etrf,asym))

        dsigpt = dsigpt + 0.5d0 * (ddsigpt0+ddsigpt) * (asym-asym0)

       enddo

*-----------------------------------------------------------------------
*     integrate total cross section of pair production
*-----------------------------------------------------------------------

       if( dsigpt .lt. 0.d0 ) dsigpt = 0.d0

       if( i .ge. 2 ) then    ! avoid first step
        sigpt = sigpt + 0.5d0 * (dsigpt0+dsigpt) * (etrf-etrf0) * 1d+24 ! [cm2 to barn]
       endif

      enddo

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      subroutine ppdmu(eein,wgti,ztar,atar,ireg,imat)
*                                                                      *
*     sampling routine about pair production event                     *
*                                                                      *
*     made by S.Abe on 2015/07/10                                      *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input :                                                          *
*        eein  : kinetic energy of muon [MeV]                          *
*        wgti  : weight                                                *
*        ztar  : atomic number of target                               *
*        atar  : mass number of target                                 *
*        ireg  : region number of cell                                 *
*        imat  : material number of cell                               *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*     muint  :                                                         *
*        imuppd : switch of muon-induced pair production               *
*                                                                      *
*     muflag :                                                         *
*        imuppdhit : flag for muon-induced pair production             *
*                                                                      *
*     icomon :                                                         *
*        ityp   : particle type for incident particle                  *
*        ktyp   : kf-code for incident particle                        *
*        jtyp   : charge for incident particle                         *
*                                                                      *
*---- Reference -------------------------------------------------------*
*                                                                      *
*     [1] A.V.Ginneken,Nucl.Instr.Meth.A,251,21,(1986)                 *
*                                                                      *
************************************************************************
      implicit real*8 (a-h,o-z)

      include 'param00.inc'
      include 'ggsparam.inc'      ! use /gg005/ and /gg015/ for TTB calculation
      include 'ggmparam.inc'      ! use /gm001/ and /gm004/ for TTB calculation

*-----------------------------------------------------------------------
! switch and flag about muon interaction processes

      common /muint/ imuint, imubrm, imuppd, imucap
      common /muflag/ imuinthit,imubrmhit,imuppdhit,imucaphit,imucapflag
!$OMP THREADPRIVATE(/muflag/)

*-----------------------------------------------------------------------
! to store information of pair particles

      common /egsemi/ iegsemi, iegsout

      common /icomon/ no,mmat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)

*-----------------------------------------------------------------------
cABE 2017/02/08
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)

      common /paraj/  mstz(300), parz(300) ! S.Abe 2018/06/05

*-----------------------------------------------------------------------

      dimension sigpt(0:1000)
      dimension dumary(3)
      dimension usave(3)

      parameter ( eneper = 2.7182818d0 )        ! Napier's constant
      parameter ( remass = 0.510999d0 )         ! electron mass [MeV/c**2]
      parameter ( rmumas = 105.6583715d0 )      ! muon mass [MeV/c**2]

*-----------------------------------------------------------------------
*     initialize
*-----------------------------------------------------------------------

      if (imuppd .ne. 1) return

      imuppdhit = 1

      do i = 0, 1000
       sigpt(i) = 0.d0
      enddo

      nclst = 0
      nclsts = 0

*-----------------------------------------------------------------------

      emutot = eein + rmumas

      wgt = wgti
      wg0 = wgti

      uuu  = 0.d0
      vvv  = 0.d0
      www  = 1.d0

      uold(1) = uuu
      uold(2) = vvv
      uold(3) = www

      do i = 0, 20
       numpal(i) = 0
       rumpal(i) = 0.d0
      enddo

      itz = idnint(ztar)
      ita = idnint(atar)
      itn = ita - itz

      zt3r = ztar**(-1.d0/3.d0)
      enp2r = dsqrt(eneper)

*-----------------------------------------------------------------------

      rchk = 4.d0*remass + 3.d0 * enp2r * rmumas / zt3r / 4.d0
      if( emutot .le. rchk ) return

*-----------------------------------------------------------------------
*     set range and width for fraction of transferred energy
*     note: etrf in this routine is within the range of 0.d0 ot 1.d0
*-----------------------------------------------------------------------

      etrfmin = 4.d0 * remass / emutot
      etrfmin = dmax1(1.d-5,etrfmin)

      etrfmax = 1.d0 - 3.d0 * enp2r * rmumas / zt3r / (4.d0 * emutot)
      etrfmax = dmin1(1.d0-1.d-6,etrfmax)

      netrf = 100
      detrf = dlog(etrfmax/etrfmin) / dble(netrf)

*-----------------------------------------------------------------------

      asymmin = 1.d-6
      nasym = 100

*-----------------------------------------------------------------------
*     calculate differential cross section for etrf
*-----------------------------------------------------------------------

      etrf = 0.d0
      dsigpt = 0.d0

      nsig = 0

      do i = 1, netrf

       etrf0 = etrf
       dsigpt0 = dsigpt

       etrf = etrfmin * dexp((dble(i))*detrf)
       dsigpt = 0.d0

*-----------------------------------------------------------------------

       asymmax = (1.d0 - 6.d0*rmumas**2.d0/emutot**2.d0/(1.d0-etrf))
     &          * dsqrt(1.d0-4.d0*remass/emutot/etrf)
       if( asymmax .le. asymmin ) exit

       dasym = dlog(asymmax/asymmin) / dble(nasym)

       asym = 0.d0
       ddsigpt = dmax1(0.d0,fmupp(itz,eein,etrf,asym))

       do j = 0, nasym

        asym0 = asym
        ddsigpt0 = ddsigpt

        asym = asymmin * dexp((dble(j))*dasym)
        ddsigpt = dmax1(0.d0,fmupp(itz,eein,etrf,asym))

        dsigpt = dsigpt + 0.5d0 * (ddsigpt0+ddsigpt) * (asym-asym0)

       enddo

*-----------------------------------------------------------------------
*     integrate total cross section for pair production
*-----------------------------------------------------------------------

       if( dsigpt .lt. 0.d0 ) dsigpt = 0.d0

       if( i .ge. 2 ) then
        nsig = nsig + 1
        sigpt(nsig) = sigpt(nsig-1)
     &                 + 0.5d0 * (dsigpt0+dsigpt) * (etrf-etrf0)
       endif

      enddo

*-----------------------------------------------------------------------
*     sampling transferred energy
*-----------------------------------------------------------------------

      if( sigpt(nsig).le.0.d0 ) return

      sigset = unirn(dummy) * sigpt(nsig)
      do i = 1, nsig
       if( sigset.gt.sigpt(i-1) .and. sigset.le.sigpt(i) ) then
        etrf0 = etrfmin * dexp((dble(i-1))*detrf)
        etrf1 = etrfmin * dexp((dble(i))*detrf)
        etrfset = etrf0 + (etrf1-etrf0) * (sigpt(i)-sigset)
     &                                  / (sigpt(i)-sigpt(i-1))
        exit
       endif
      enddo

*-----------------------------------------------------------------------
*     calculate asymmetry parameter for etrfset
*-----------------------------------------------------------------------

! recycle sigpt(i)
      do i = 0, nsig
       sigpt(i) = 0.d0
      enddo
      nsig = 0

! to avoid error
      asymmax = (1.d0 - 6.d0*rmumas**2.d0/emutot**2.d0/(1.d0-etrfset))
     &         * dsqrt(1.d0-4.d0*remass/emutot/etrfset)

      if( asymmax .gt. asymmin ) then

       dasym = dlog(asymmax/asymmin) / dble(nasym)

       asym = 0.d0
       ddsigpt = dmax1(0.d0,fmupp(itz,eein,etrfset,asym))

       do i = 0, nasym

        asym0 = asym
        ddsigpt0 = ddsigpt

        asym = asymmin * dexp((dble(i))*dasym)
        ddsigpt = dmax1(0.d0,fmupp(itz,eein,etrfset,asym))

        nsig = nsig + 1
        sigpt(nsig) = sigpt(nsig-1)
     &               + 0.5d0*(ddsigpt0+ddsigpt)*(asym-asym0)

       enddo

*-----------------------------------------------------------------------
*     sampling asymmetry parameter
*-----------------------------------------------------------------------

       asymset = 0.d0

       if( sigpt(nsig).gt.0.d0 ) then

        sigset = unirn(dummy) * sigpt(nsig)
        if( sigset.ge.sigpt(0) .and. sigset.le.sigpt(1)) then

         asym0 = 0.d0
         asym1 = asymmin 
          asymset = asym0 + (asym1-asym0) * (sigpt(1)-sigset)
     &                                    / (sigpt(1)-sigpt(0))
        else

         do i = 2, nsig
          if( sigset.gt.sigpt(i-1) .and. sigset.le.sigpt(i) ) then
           asym0 = asymmin * dexp((dble(i-2))*dasym)
           asym1 = asymmin * dexp((dble(i-1))*dasym)
           asymset = asym0 + (asym1-asym0) * (sigpt(i)-sigset)
     &                                     / (sigpt(i)-sigpt(i-1))
           exit
          endif
         enddo

        endif

       endif

*-----------------------------------------------------------------------

      else

       asymset = 0.d0

      endif

*-----------------------------------------------------------------------
*     sampling scattering angles of muon and pair particle
*     by using eq.(26) of [1]
*-----------------------------------------------------------------------

      emutotg = emutot * 1.d-3
      remassg = remass * 1.d-3

      rmsthemu = (2.3d0+dlog(emutotg))/emutotg/(1.d0-etrfset)
     &          * (etrfset-2.d0*remassg/emutotg)**2.d0 /etrfset**2.d0
     &          * dmin1(8.9d-4*etrfset**0.25d0*(1.d0+1.5d-5*emutotg)
     &                  +0.032d0*etrfset/(etrfset+1.d0), 0.1d0)

      themumin = 0.d0
      themumax = 5.d0 * dsqrt(rmsthemu**2.d0/2.d0)
      nthemu = 1000
      dthemu = (themumax-themumin)/dble(nthemu)

      themu = themumin
      gthemu = 1.d0 / dsqrt(pi*rmsthemu**2.d0) * dthemu
      sigset = unirn(dummy) * 0.5d0
      dgthemu = 0.d0

      do i = 1, nthemu

       themu0 = themu
       gthemu0 = gthemu
       dgthemu0 = dgthemu

       themu = dthemu * dble(i) + themumin
       gthemu = dexp(-(themu**2.d0/rmsthemu**2.d0))
     &         / dsqrt(pi*rmsthemu**2.d0) * dthemu
       dgthemu = dgthemu + (gthemu+gthemu0) * (themu-themu0) * 0.5d0

       if( dgthemu .ge. sigset ) then
        themu = themu0 + (themu-themu0) * (dgthemu-sigset)
     &                                  / (dgthemu-dgthemu0)
        exit
       endif

      enddo

      costhemu = dcos(themu)
      sinthemu = dsin(themu)

      phimu = 2.d0 * pi * unirn(dummy)   ! S.Abe 2015/12/04, corrected
      cosphimu = dcos(phimu)
      sinphimu = dsin(phimu)

*-----------------------------------------------------------------------
*     momentum of virtual photon
*-----------------------------------------------------------------------

      pmu0abs = dsqrt(emutot**2.d0-rmumas**2.d0)
      emu1tot = (1.d0 - etrfset) * emutot
      pmu1abs = dsqrt(emu1tot**2.d0-rmumas**2.d0)

      pvpabs2 = pmu0abs**2.d0 + pmu1abs**2.d0 - 2.d0 * pmu0abs * pmu1abs
      pvpabs = dsqrt(pvpabs2)

      sinthevp = pmu1abs * sinthemu / pvpabs
      costhevp = dsqrt(1.d0 - sinthevp**2.d0)

      cosphivp = -cosphimu
      sinphivp = -sinphimu

      q2 = pvpabs2 - etrfset * emutot   ! not used but calculated

*-----------------------------------------------------------------------
*     momentum of pair particles
*-----------------------------------------------------------------------

      epair1 = (1.d0 + asymset) * 0.5d0 * etrfset * emutot
      epair2 = (1.d0 - asymset) * 0.5d0 * etrfset * emutot

      p1abs = dsqrt(epair1**2.d0 - remass**2.d0)
      p2abs = dsqrt(epair2**2.d0 - remass**2.d0)

      costhe1 = (pvpabs**2.d0 + p1abs**2.d0 - p2abs**2.d0) * 0.5d0
     &         / pvpabs / p1abs

      if( costhe1 .gt. 1.d0 ) costhe1 = 1.d0
      if( costhe1 .lt. -1.d0 ) costhe1 = -1.d0
      sinthe1 = dsqrt(1.d0 - costhe1**2.d0)

      costhe2 = (pvpabs - p1abs * costhe1) / p2abs
      if( costhe2 .gt. 1.d0 ) costhe2 = 1.d0
      if( costhe2 .lt. -1.d0 ) costhe2 = -1.d0
      sinthe2 = dsqrt(1.d0 - costhe2**2.d0)

      phi1 = 2.d0 * pi * unirn(dummy)   ! S.Abe 2015/12/04, corrected
      cosphi1 = dcos(phi1)
      sinphi1 = dsin(phi1)

      cosphi2 = -cosphi1
      sinphi2 = -sinphi1

*-----------------------------------------------------------------------
*     bank scattered muon and pair particles
*-----------------------------------------------------------------------

      nclsts = nclsts + 1
      iclusts(nclsts) = 6

      jclusts(0,nclsts) = 0
      jclusts(1,nclsts) = 0
      jclusts(2,nclsts) = 0
      jclusts(3,nclsts) = ityp
      jclusts(4,nclsts) = 0
      jclusts(5,nclsts) = jtyp
      jclusts(6,nclsts) = 0
      jclusts(7,nclsts) = ktyp
      jclusts(8,nclsts) = 0

      rms = rmumas
      rmg = rms * 1.d-3
      pr  = pmu1abs * 1.d-3
      ett = dsqrt(pr**2.d0 + rmg**2.d0)
      ek  = (ett - rmg) * 1.d+3
      pxl = pr * sinthemu * cosphimu
      pyl = pr * sinthemu * sinphimu
      pzl = pr * costhemu

      qclusts(0,nclsts) = 0.d0
      qclusts(1,nclsts) = pxl
      qclusts(2,nclsts) = pyl
      qclusts(3,nclsts) = pzl
      qclusts(4,nclsts) = ett
      qclusts(5,nclsts) = rmg
      qclusts(6,nclsts) = 0.d0
      qclusts(7,nclsts) = ek
      qclusts(8,nclsts) = wgt / wg0
      qclusts(9,nclsts) = 0.d0
      qclusts(10,nclsts) = 0.d0
      qclusts(11,nclsts) = 0.d0
      qclusts(12,nclsts) = 0.d0

      numpal(ityp) = numpal(ityp) + 1
      rumpal(ityp) = rumpal(ityp) + wgti

*-----------------------------------------------------------------------

cABE 2017/02/08
      if( ityp .eq. 7 ) atmrc(4,5) = atmrc(4,5) + wgt / wg0
      if( ityp .eq. 6 ) atmrc(5,5) = atmrc(5,5) + wgt / wg0

*-----------------------------------------------------------------------
*     bank pair particles
*-----------------------------------------------------------------------

      rjudge = unirn(dummy)
      if( rjudge .le. 0.5d0 ) then
       ip1 = 12
       jp1 = -1
       kp1 = 11
       ip2 = 13
       jp2 = 1
       kp2 = -11
      else
       ip1 = 13
       jp1 = 1
       kp1 = -11
       ip2 = 12
       jp2 = -1
       kp2 = 11
      endif

! roulette or split according to enum bias.
      if( enum .eq. 1.d0 ) then
       npa = 1
       wga = wgt
      else
       npa = idnint(enum + rang())
       wgt = wgt / enum
       wga = npa * wgt
      endif

*-----------------------------------------------------------------------
*     particle 1
*-----------------------------------------------------------------------

      ittbflag = 0
      if( mstz(85) .lt.1 .and. ides .eq. 0 .and. nee .ne. 0 .and.
     &    epair1 .gt. elc(2) .and. epair1 .lt. elc(3) ) then
       ittbflag = 1
      endif

*-----------------------------------------------------------------------

      if( ittbflag .eq. 0 ) then

       nclsts = nclsts + 1
       iclusts(nclsts) = 7

       jclusts(0,nclsts) = 0
       jclusts(1,nclsts) = 0
       jclusts(2,nclsts) = 0
       jclusts(3,nclsts) = ip1
       jclusts(4,nclsts) = 0
       jclusts(5,nclsts) = jp1
       jclusts(6,nclsts) = 0
       jclusts(7,nclsts) = kp1
       jclusts(8,nclsts) = 0

       rms = remass
       rmg = rms * 1.d-3
       pr  = p1abs * 1.d-3
       ett = dsqrt(pr**2.d0 + rmg**2.d0)
       ek  = (ett - rmg) * 1.d+3
       pxl0 = pr * sinthe1 * cosphi1
       pyl0 = pr * sinthe1 * sinphi1
       pzl0 = pr * costhe1
! rotation about direction of virtual photon 
       pxl =  pxl0 * costhevp * cosphivp
     &      - pyl0 * sinphivp
     &      + pzl0 * sinthevp * cosphivp
       pyl =  pxl0 * costhevp * sinphivp
     &      + pyl0 * cosphivp
     &      + pzl0 * sinthevp * sinphivp
       pzl = -pxl0 * sinthevp
     &      + pzl0 * costhevp

       qclusts(0,nclsts) = 0.d0
       qclusts(1,nclsts) = pxl
       qclusts(2,nclsts) = pyl
       qclusts(3,nclsts) = pzl
       qclusts(4,nclsts) = ett
       qclusts(5,nclsts) = rmg
       qclusts(6,nclsts) = 0.d0
       qclusts(7,nclsts) = ek
       qclusts(8,nclsts) = wgt / wg0
       qclusts(9,nclsts) = 0.d0
       qclusts(10,nclsts) = 0.d0
       qclusts(11,nclsts) = 0.d0
       qclusts(12,nclsts) = 0.d0

       numpal(ip1) = numpal(ip1) + 1
       rumpal(ip1) = rumpal(ip1) + wgt

! make ttb photons
      elseif( npa .gt. 0 ) then

       erg0 = erg
       usave(1) = uuu
       usave(2) = vvv
       usave(3) = www

       erg = epair1 - remass
       uuu =  sinthe1 * cosphi1 * costhevp * cosphivp
     &       -sinthe1 * sinphi1 *            sinphivp
     &      + costhe1           * sinthevp * cosphivp
       vvv =  sinthe1 * cosphi1 * costhevp * sinphivp
     &      + sinthe1 * sinphi1 *            cosphivp
     &      + costhe1           * sinthevp * sinphivp
       www =  sinthe1 * cosphi1 * sinthevp
     &      + costhe1           * costhevp
       if( erg .gt. elc(2) ) then
        do j = 1, npa
         dumary(1) = 0.d0
         dumary(2) = 0.d0
         dumary(3) = 0.d0
         call brmgam(0.d0,ireg,imat,dumary(3),0.d0)
        enddo
       endif

       erg = erg0
       uuu = usave(1)
       vvv = usave(2)
       www = usave(3)

      endif

*-----------------------------------------------------------------------
* particle2
*-----------------------------------------------------------------------

      ittbflag = 0
      if( ides .eq. 0 .and. nee .ne. 0 .and.
     &    epair2 .gt. elc(2) .and. epair2 .lt. elc(3) ) then
       ittbflag = 1
      endif

*-----------------------------------------------------------------------

      if( ittbflag .eq. 0 ) then

       nclsts = nclsts + 1
       iclusts(nclsts) = 7

       jclusts(0,nclsts) = 0
       jclusts(1,nclsts) = 0
       jclusts(2,nclsts) = 0
       jclusts(3,nclsts) = ip2
       jclusts(4,nclsts) = 0
       jclusts(5,nclsts) = jp2
       jclusts(6,nclsts) = 0
       jclusts(7,nclsts) = kp2
       jclusts(8,nclsts) = 0

       rms = remass
       rmg = rms * 1.d-3
       pr  = p2abs * 1.d-3
       ett = dsqrt(pr**2.d0 + rmg**2.d0)
       ek  = (ett - rmg) * 1.d+3
       pxl0 = pr * sinthe2 * cosphi2
       pyl0 = pr * sinthe2 * sinphi2
       pzl0 = pr * costhe2
! rotation about direction of virtual photon 
       pxl =  pxl0 * costhevp * cosphivp
     &       -pyl0 * sinphivp
     &      + pzl0 * sinthevp * cosphivp
       pyl =  pxl0 * costhevp * sinphivp
     &      + pyl0 * cosphivp
     &      + pzl0 * sinthevp * sinphivp
       pzl = -pxl0 * sinthevp
     &      + pzl0 * costhevp

       qclusts(0,nclsts) = 0.d0
       qclusts(1,nclsts) = pxl
       qclusts(2,nclsts) = pyl
       qclusts(3,nclsts) = pzl
       qclusts(4,nclsts) = ett
       qclusts(5,nclsts) = rmg
       qclusts(6,nclsts) = 0.d0
       qclusts(7,nclsts) = ek
       qclusts(8,nclsts) = wgt / wg0
       qclusts(9,nclsts) = 0.d0
       qclusts(10,nclsts) = 0.d0
       qclusts(11,nclsts) = 0.d0
       qclusts(12,nclsts) = 0.d0

       numpal(ip2) = numpal(ip2) + 1
       rumpal(ip2) = rumpal(ip2) + wgt

! make ttb photons
      elseif( npa .gt. 0 ) then

       erg0 = erg
       usave(1) = uuu
       usave(2) = vvv
       usave(3) = www

! particle 2
       erg = epair2 - remass
       uuu =  sinthe2 * cosphi2 * costhevp * cosphivp
     &       -sinthe2 * sinphi2 *            sinphivp
     &      + costhe2           * sinthevp * cosphivp
       vvv =  sinthe2 * cosphi2 * costhevp * sinphivp
     &      + sinthe2 * sinphi2 *            cosphivp
     &      + costhe2           * sinthevp * sinphivp
       www =  sinthe2 * cosphi2 * sinthevp
     &      + costhe2           * costhevp
       if( erg .gt. elc(2) ) then
        do j = 1, npa
         dumary(1) = 0.d0
         dumary(2) = 0.d0
         dumary(3) = 0.d0
         call brmgam(0.d0,ireg,imat,dumary(3),0.d0)
        enddo
       endif

       erg = erg0
       uuu = usave(1)
       vvv = usave(2)
       www = usave(3)

      endif

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      function fmupp(itz,eein,etrf,asym)
*                                                                      *
*     function of muon-induced pair prodiction                         *
*                                                                      *
*     made by S.Abe on 2015/07/10                                      *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input :                                                          *
*        itz   : atomic number of target                               *
*        eein  : kinetic energy of muon [MeV]                          *
*        etrf  : fraction of muon energy transferred to photon [-]     *
*        asym  : asymmetry coefficient of the energy of the pair [-]   *
*                                                                      *
*     output:                                                          *
*        fmupp : differential cross section for muon pair prod. [cm2]  *
*                                                                      *
*---- Reference -------------------------------------------------------*
*                                                                      *
*     [1] R.P.Kokoulin+,Proc.12th ICRC,6,2436(1969)                    *
*     [2] S.R.Kelner+,Phys.Atom.Nucl.,61,448(1998)                     *
*     [3] I.A.Sokalski+,Phys.Rev.D,64,074015(2001)                     *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      parameter ( pi = 3.141592653589793d0 )
      parameter ( eneper = 2.7182818d0 )        ! Napier's constant
      parameter ( alpha = 7.297353d-3 )         ! fine-structure constant
      parameter ( cre0 = 2.8179380d-13 )        ! classcal electron radius [cm]
      parameter ( remass = 0.510999d0 )         ! electron mass [MeV/c**2]
      parameter ( rmumas = 105.6583715d0 )      ! muon mass [MeV/c**2]
      parameter ( radlog0 = 183.d0 )            ! radiation logarithm
      parameter ( radlogh = 202.4d0 )           ! radiation logarithm for hydrogen

*-----------------------------------------------------------------------

      fmupp = 0.d0

*-----------------------------------------------------------------------
*     initialize
*-----------------------------------------------------------------------

      zt = dble(itz)

      if( itz .eq. 1 )then
       radlog = radlogh
      else
       radlog = radlog0
      endif

      emutot = eein + rmumas

      zt3r = zt**(-1.d0/3.d0)
      enp2r = dsqrt(eneper)
      asym2 = asym**2.d0

      beta = 0.5d0 * etrf**2.d0 / (1.d0 - etrf)
      xsi = (0.5d0 * rmumas * etrf / remass)**2.d0
     &     * (1.d0 - asym2) / (1.d0 - etrf)

      fmuppe = 0.d0
      fmuppm = 0.d0

c...function zeta from eq(48-50) of [2]
      if( emutot .gt. 35.d0*rmumas ) then

       if( itz .eq. 1 ) then
        gamma1 = 4.40d-5
        gamma2 = 4.80d-5
       else
        gamma1 = 1.95d-5
        gamma2 = 5.30d-5
       endif

       aaa = 0.073d0 * dlog(emutot/(rmumas+gamma1*emutot/(zt3r**2.d0)))
     &      - 0.26d0
       bbb = 0.058d0 * dlog(emutot/(rmumas+gamma2*emutot/zt3r)) - 0.14d0

       if(aaa .gt. 0.d0) then
        zeta = aaa / bbb
       else
        zeta = 0.d0
       endif

      else

       zeta = 0.d0

      endif

*-----------------------------------------------------------------------
*     calculate ddx for e-diagram
*-----------------------------------------------------------------------

c...function Ye
      aaa = 5.d0 - asym2 + 4.d0*beta*(1.d0+asym2)
      bbb = 2.d0*(1.d0+3.d0*beta)*dlog(3.d0+1.d0/xsi)
     &     - asym2 - 2.d0*beta*(2.d0-asym2)
      yelec = aaa / bbb

c...function Le
      rnumer = radlog * zt3r * dsqrt((1.d0+xsi)*(1.d0+yelec)) 
      denomi = 1.d0 + 2.d0 * remass * enp2r * radlog * zt3r
     &                * (1.d0+xsi) * (1.d0+yelec)
     &                / (emutot*etrf*(1.d0-asym2))
      aaa = dlog(rnumer/denomi)
      bbb = 0.5d0 * dlog(1.d0 + (1.5d0*remass/zt3r/rmumas)**2.d0
     &                         * (1.d0+xsi) * (1.d0+yelec)  )
      rlelec = aaa - bbb

c...function Be
      aaa = (2.d0+asym2)*(1.d0+beta) + xsi*(3.d0+asym2)
      bbb = (1.d0-asym2-beta) / (1.d0+xsi) - (3.d0+asym2)
      belec = aaa*dlog(1.d0+1.d0/xsi) + bbb

c...ddx for e-diagram
      aaa = 4.d0 / 3.d0 / pi * zt*(zt-zeta) * (alpha*cre0)**2.d0
     &     * (1.d0-etrf) / etrf
      fmuppe = aaa * belec * rlelec

*-----------------------------------------------------------------------
*     calculate ddx for mu-diagram
*-----------------------------------------------------------------------

c...function Ym
      aaa = 4.d0 + asym2 + 3.d0*beta*(1.d0+asym2)
      bbb = (1.d0+asym2)*(1.5d0 + 2.d0*beta)*dlog(3.d0+xsi)
     &     + 1.d0 - 1.5d0*asym2
      ymuon = aaa / bbb

c...function Lm
      rnumer = 2.d0/3.d0 * rmumas/remass * radlog * zt3r**2.d0
      denomi = 1.d0 + 2.d0 * remass * enp2r * radlog * zt3r
     &                * (1.d0+xsi) * (1.d0+ymuon)
     &                / (emutot*etrf*(1.d0-asym2))
      rlmuon = dlog(rnumer/denomi)

c...function Bm
      aaa = (1.d0+asym2)*(1.d0+1.5d0*beta)
     &     - (1.d0-asym2)*(1.d0+2.d0*beta)/xsi
      bbb = xsi*(1.d0-asym2-beta)/(1.d0+xsi)
     &     + (1.d0-asym2)*(1.d0+2.d0*beta)
      bmuon = aaa*dlog(1+xsi) + bbb

c...ddx for mu-diagram
      aaa = 4.d0 / 3.d0 / pi * zt * (zt-zeta) * (alpha*cre0)**2.d0
     &     * (1.d0-etrf) / etrf
      bbb = (remass/rmumas)**2.d0
      fmuppm = aaa * bbb * bmuon * rlmuon

*-----------------------------------------------------------------------

      fmupp = fmuppe + fmuppm

*-----------------------------------------------------------------------
      return
      end
