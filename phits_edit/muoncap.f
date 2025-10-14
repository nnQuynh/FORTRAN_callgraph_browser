************************************************************************
*                                                                      *
      subroutine capmujudge(ireg,imat,imucapflag,itz,ita)
*                                                                      *
*        judges captured by nucleon or decay in 1s orbit               *
*        last modified by S.Abe on 2015/04/24                          *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input  :                                                         *
*        ireg   : target region                                        *
*        imat   : target material                                      *
*                                                                      *
*     output :                                                         *
*        imucapflag : flag for capture or decay, =1: capture, =2:decay *
*        itz    : atomic number of target                              *
*        ita    : mass number of target                                *
*                                                                      *
************************************************************************
!<-20211126murofushi add
      use moddas_material
!--->

      implicit real*8 (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

!<-20211126murofushi delete
!--      dimension dnel(1), denh(1), zz(1), a(1), den(1)
!--      equivalence ( das, dnel, denh, zz, a, den )
!--->

      common /kmat1g/ kmat(kvlmax)

*-----------------------------------------------------------------------

      if( imat .eq. 0 ) then
       imucapflag = 0
       return
      endif

*-----------------------------------------------------------------------
*     select mother nucleus by Fermi-Teller law
*       nlap = 1: sum of total capture rate
*       nlap = 2: sampling for itz, ita and capture or decay
*-----------------------------------------------------------------------

      probcd = 0.d0

!<-20220916murofushi update
!--      lem  = idnint( dnel(kmat(imat)+1) )
!--      hydro = denh(kmat(imat)+2)    ! S.Abe 2015/08/10
      lem = idnint( dnel_das(kmat0+imat) )
      hydro = denh_das(kmat0+imat)    ! S.Abe 2015/08/10
!--->

      do nlap = 1, 2

       if( nlap .eq. 2 ) then
        rselect = unirn(dummy) * probcd
       endif

       do i = 1, lem

!<-20211126murofushi update
!--        zt = zz(kmat(imat)+(i-1)*3+3)
!--        at =  a(kmat(imat)+(i-1)*3+4)
        zt = zz_das(kmat(imat)+i)
        at = a_das(kmat(imat)+i)
!--->
        itz = idnint(zt)
        ita = idnint(at)

!<-20211126murofushi update
!--        capr = dmax1(0.d0, caprmu(itz,ita)) * den(kmat(imat)+(i-1)*3+5)
!--        dcyr = dmax1(0.d0, dcyrmu(itz)) * den(kmat(imat)+(i-1)*3+5)
        capr = dmax1(0.d0, caprmu(itz,ita))*den_das(kmat(imat)+i)
        dcyr = dmax1(0.d0, dcyrmu(itz))*den_das(kmat(imat)+i)
!--->

        if( nlap .eq. 1 ) then

         probcd = probcd + capr + dcyr

        elseif( nlap .eq. 2 ) then

         rselect = rselect - capr
         if( rselect .le. 0.d0 ) then
          imucapflag = 1   ! go to capture by nucleus
          goto 9999   ! S.Abe 2016/04/19, corrected
         endif

         rselect = rselect - dcyr
         if( rselect .le. 0.d0 ) then
          imucapflag = 2   ! go to decay
          goto 9999   ! S.Abe 2016/04/19, corrected
         endif

        endif

       enddo

       if( hydro .gt. 0.d0 ) then

        itz = 1
        ita = 1

        capr = dmax1(0.d0, caprmu(itz,ita)) * hydro
        dcyr = dmax1(0.d0, dcyrmu(itz)) * hydro

        if( nlap .eq. 1 ) then

         probcd = probcd + capr + dcyr

        elseif( nlap .eq. 2) then

         rselect = rselect - capr
         if( rselect .le. 0.d0 ) then
          imucapflag = 1   ! go to capture by nucleus
          goto 9999   ! S.Abe 2016/04/19, corrected
         endif

         rselect = rselect - dcyr
         if( rselect .le. 0.d0 ) then
          imucapflag = 2   ! go to decay
          goto 9999   ! S.Abe 2016/04/19, corrected
         endif

        endif

       endif

      enddo

*-----------------------------------------------------------------------
 9999 continue

      return
      end


************************************************************************
*                                                                      *
      subroutine capmuxray(imucapflag,wgti,itz,ita,emubind)
*                                                                      *
*        calculate a capture of a negative muon with an atom.          *
*        last modified by S.Abe on 2015/12/04                          *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input  :                                                         *
*        imucapflag : flag of muon capture                             *
*        wgti : weight                                                 *
*        itz  : atomic number of target                                *
*        ita  : mass number of target                                  *
*                                                                      *
*     output :                                                         *
*        emubind    : binding energy of muon on 1s orbit [MeV]         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'param00.inc'

      common /aama99/ popul(20),ptran(210,210),cctr(210,210),
     &                esp0,energy0(20,40)
!$OMP THREADPRIVATE(/aama99/)
      common /aamai99/ nmax0
!$OMP THREADPRIVATE(/aamai99/)

*-----------------------------------------------------------------------

      common /mueng/ emu1s, rtamas1, rtamas2, excit0
!$OMP THREADPRIVATE(/mueng/)

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)

cabe add @2015/04/09
      common /muxray/  pxsum, pysum, pzsum, nxray
!$OMP THREADPRIVATE(/muxray/)

*-----------------------------------------------------------------------
cABE 2017/02/08
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)
      common /tpdcta/ nflumu
!$OMP THREADPRIVATE(/tpdcta/)

*-----------------------------------------------------------------------

      dimension exray(4), probxray(4)

      parameter ( pi = 3.141592653589793d0 )

*-----------------------------------------------------------------------
*     initialize
*-----------------------------------------------------------------------

      nclsts = 0
      nflumu = 0    ! S.Abe 2017/02/08

      esp0 = 0.d0
      nmax0 = 0

      do i = 1, 20
       popul(i) = 0.d0
      enddo

      do i = 1, 210
       do j = 1, 210
        ptran(i,j) = 0.d0
        cctr(i,j) = 0.d0
       enddo
      enddo

      emubind = 0.d0

      pxsum = 0.d0
      pysum = 0.d0
      pzsum = 0.d0
      etsum = 0.d0

cabe add @2015/05/17
      do i = 0, 20
       numpal(i) = 0
       rumpal(i) = 0.0d0
      enddo

cABE 2017/02/08
      wgt = wgti
      wg0 = wgti

*-----------------------------------------------------------------------

      if( imucapflag .eq. 0 ) return
      call aamain(itz,ita)

*-----------------------------------------------------------------------
*     sampling l1
*-----------------------------------------------------------------------

      n1 = nmax0

      popult = 0.d0

      do i = 1, n1
       popult = popult + popul(i)
      enddo

      rpop = unirn(dummy) * popult

      do i = 1, n1
       rpop = rpop - popul(i)
       if( rpop .le. 0.d0 ) then
        l1 = i-1
        goto 1000
       endif
      enddo

*-----------------------------------------------------------------------
*     sampling n2 and l2 
*-----------------------------------------------------------------------

1000  continue

      ptrant = 0.d0
      erg = 0.d0

      do i = 1, n1
       n2 = i
       do j = 1, 7
        l2 = l1 - 4 + j
        if( l2 .lt. 0 ) cycle
        if( l2 .ge. n2 ) cycle
        if( n2.eq.n1 .and. l2.eq.l1 ) cycle
        kk1 = n1*(n1-1)/2 + l1 + 1
        kk2 = n2*(n2-1)/2 + l2 + 1
        ptrant = ptrant + ptran(kk1,kk2)
       enddo
      enddo

      rptran = unirn(dummy) * ptrant

      do i = 1, n1
       n2 = i
       do j = 1, 7
        l2 = l1 - 4 + j
        if( l2 .lt. 0 ) cycle
        if( l2 .ge. n2 ) cycle
        if( n2.eq.n1 .and. l2.eq.l1 ) cycle
        kk1 = n1*(n1-1)/2 + l1 + 1
        kk2 = n2*(n2-1)/2 + l2 + 1
        rptran = rptran - ptran(kk1,kk2)
        if( rptran .le. 0.d0) goto 2000
       enddo
      enddo

*-----------------------------------------------------------------------
*     defined transition energy
*-----------------------------------------------------------------------

2000  continue

      if( l1.eq.0 .and. l2.eq.0 ) goto 3000
      if( cctr(kk1,kk2) .lt. 1.d-3 ) goto 3000

      ldelt = iabs(l1-l2)
      if( ldelt .eq. 0) ldelt = 2

*-----------------------------------------------------------------------

      e12 = 0.d0
      e22 = 0.d0

      e11 = energy0(n1,2*l1+1)
      if( l1.ne.0 ) e12 = energy0(n1,2*l1)

      e21 = energy0(n2,2*l2+1)
      if( l2.ne.0 ) e22 = energy0(n2,2*l2)

      if( n1.eq.n2 ) then
       e11 = 0.d0
       ez = 0.5d0*(energy0(2,2)-energy0(2,3))
       e21 = esp0*1.d-6 - ez
       e22 = e21 + 2.d0*ez
      endif

*-----------------------------------------------------------------------

      if( ldelt .eq. 1 ) then

       lmax = max0(l1,l2)

       exray(1) = e21 - e11
       if(l1.eq.lmax) then
        exray(2) = e21 - e12
       elseif(l2.eq.lmax) then
        exray(2) = e22 - e11
       endif
       exray(3) = e22 - e12

       aprobx = cctr(kk1,kk2) / dble(4*lmax*lmax-1)
       probxray(1) = aprobx * dble((lmax+1)*(2*lmax-1))
       probxray(2) = aprobx
       probxray(3) = aprobx * dble((lmax-1)*(2*lmax+1))
       probxray(4) = 0.d0

*-----------------------------------------------------------------------

      elseif( ldelt .eq. 2 ) then

       if( l1 .ne. l2 ) then

        lmax = max0(l1,l2)

        exray(1) = e21 - e11
        if(l1.eq.lmax) then
         exray(2) = e21 - e12
        elseif(l2.eq.lmax) then
         exray(2) = e22 - e11
        endif
        exray(3) = e22 - e12

        aprobx = cctr(kk1,kk2) / dble(4*lmax*lmax-1)
        probxray(1) = aprobx * dble((lmax+1)*(2*lmax-1))
        probxray(2) = aprobx
        probxray(3) = aprobx * dble((lmax-1)*(2*lmax+1))
        probxray(4) = 0.d0

       else

        exray(1) = e21 - e11
        exray(2) = e22 - e11
        exray(3) = e21 - e12
        exray(4) = e22 - e12

        aprobx = cctr(kk1,kk2) / dble((2*l1+1)**2)
        probxray(1) = aprobx * dble((l1+2)*(2*l1-1))
        probxray(2) = 3.d0 * aprobx
        probxray(3) = 3.d0 * aprobx
        probxray(4) = aprobx * dble((l1-1)*(2*l1+3))

       endif

*-----------------------------------------------------------------------

      elseif( ldelt .eq. 3 ) then

       if( iabs(l1-l2) .ne. 1 ) then

        lmax = max0(l1,l2)

        exray(1) = e21 - e11
        if(l1.eq.lmax) then
         exray(2) = e21 - e12
        elseif(l2.eq.lmax) then
         exray(2) = e22 - e11
        endif
        exray(3) = e22 - e12

        aprobx = cctr(kk1,kk2) / dble((2*lmax-5)*(2*lmax+1))
        probxray(1) = aprobx * dble((2*lmax-5)*(lmax+1))
        probxray(2) = 3.d0 * aprobx
        probxray(3) = aprobx * dble((2*lmax+1)*(lmax-3))
        probxray(4) = 0.d0

       else

        lmax = max0(l1,l2)

        exray(1) = e21 - e11
        exray(2) = e22 - e11
        exray(3) = e21 - e12
        exray(4) = e22 - e12

        aprobx = cctr(kk1,kk2) / dble(4*ldelt*ldelt-1)
        probxray(1) = aprobx * dble((2*ldelt-3)*(ldelt+2))
        probxray(2) = 5.d0 * aprobx
        probxray(3) = 6.d0 * aprobx
        probxray(4) = aprobx * dble((2*ldelt+3)*(ldelt-2))

       endif

      endif

*-----------------------------------------------------------------------

      probxt = probxray(1) + probxray(2) + probxray(3) + probxray(4)

      rprobx = unirn(dummy) * probxt
      do i= 1, 4
        rprobx = rprobx - probxray(i)
        if( rprobx .le. 0.d0) then
         erg = exray(i)
         goto 3000
        endif
      enddo

*-----------------------------------------------------------------------
*     store information of x-rays in nclsts array
*-----------------------------------------------------------------------

3000  continue

      if( erg .gt. 1.d-6 ) then     ! cutoff energy

       nclsts = nclsts + 1
       iclusts(nclsts) = 4

       jclusts(0,nclsts) = 0
       jclusts(1,nclsts) = 0
       jclusts(2,nclsts) = 0
       jclusts(3,nclsts) = 14
       jclusts(4,nclsts) = 0     ! S.Abe 2015/06/24 corrected
       jclusts(5,nclsts) = 0
       jclusts(6,nclsts) = 0
       jclusts(7,nclsts) = 22
       jclusts(8,nclsts) = 0

       pr = erg * 1.d-3
       thetax = pi * unirn(dummy)   ! S.Abe 2015/12/04, corrected
       phix = 2.d0 * pi * unirn(dummy)

       pxl = pr * dsin(thetax) * dcos(phix)
       pyl = pr * dsin(thetax) * dsin(phix)
       pzl = pr * dcos(thetax)

       qclusts(0,nclsts)  = 0.d0
       qclusts(1,nclsts)  = pxl
       qclusts(2,nclsts)  = pyl
       qclusts(3,nclsts)  = pzl
       qclusts(4,nclsts)  = pr
       qclusts(5,nclsts)  = 0.d0
       qclusts(6,nclsts)  = 0.d0
       qclusts(7,nclsts)  = erg
       qclusts(8,nclsts)  = wgt / wg0
       qclusts(9,nclsts)  = 0.d0
       qclusts(10,nclsts) = 0.d0
       qclusts(11,nclsts) = 0.d0
       qclusts(12,nclsts) = 0.d0

       numpal(14) = numpal(14) + 1
       rumpal(14) = rumpal(14) + wgti

       pxsum = pxsum + pxl
       pysum = pysum + pyl
       pzsum = pzsum + pzl
       etsum = etsum + erg

      endif

*-----------------------------------------------------------------------
*     overwrite n1 and l1, and return to get next n2 and l2
*-----------------------------------------------------------------------

      if( n2 .ne. 1 ) then
       n1 = n2
       l1 = l2
       goto 1000
      endif

*-----------------------------------------------------------------------

 4000 continue

      emu1s = energy0(1,1)
      emubind = emu1s

      nxray = nclsts

cABE 2017/02/08
      if( nclsts .gt. 0 ) atmrc(4,1) = atmrc(4,1) + wgt / wg0
      nflumu = nclsts

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      subroutine capmu(wgti,itz,ita)
*                                                                      *
*        calculate capture of negative muon by nucleon.                *
*        last modified by S.Abe on 2015/06/24                          *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input  :                                                         *
*        wgti : weight                                                 *
*        itz  : atomic number of target                                *
*        ita  : mass number of target                                  *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use NGSDATAMOD, only : bindeg
*-----------------------------------------------------------------------

      implicit real*8(a-h,o-z)

      include 'param00.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

      common /qmdflg/ irqmd

      common /muflag/ imuinthit,imubrmhit,imuppdhit,imucaphit,imucapflag
!$OMP THREADPRIVATE(/muflag/)
      common /mueng/ emu1s, rtamas1, rtamas2, excit0
!$OMP THREADPRIVATE(/mueng/)

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn),  qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)

      common /muxray/  pxsum, pysum, pzsum, nxray
!$OMP THREADPRIVATE(/muxray/)

      common /capmudir/ thetacap, phicap
!$OMP THREADPRIVATE(/capmudir/)

      parameter ( pi = 3.141592653589793d0 )
      parameter ( rpmass = 938.272046d0 )       ! proton mass [MeV/c**2]
      parameter ( rnmass = 939.565379d0 )       ! neutron mass [MeV/c**2]
      parameter ( rmumas = 105.6583715d0 )      ! muon mass [MeV/c**2]

      dimension pmn(5)
      dimension pnucl(5)
cABE 2015/06/24, changed
      dimension pini(3)
      dimension pirot(3)
      dimension pfrot(3)
      dimension pfin(3)

cABE 2019/03/08
      common /cascid/ jcasc
!$OMP THREADPRIVATE(/cascid/)


*-----------------------------------------------------------------------
*     initialize
*-----------------------------------------------------------------------

      imucaphit = 1

      nclst = -1
      ipim = 0

      wgt  = wgti
      wg0  = wgti

      erg = emu1s

      itz1 = itz
      itn1 = ita - itz1
      be1 = bindeg(itz1,itn1)
      rtamas1 = ( itz1 * rpmass + itn1 * rnmass - be1 )

      itz2 = itz - 1
      itn2 = ita - itz2
      be2 = bindeg(itz2,itn2)
      rtamas2 = ( itz2 * rpmass + itn2 * rnmass - be2 )

*-----------------------------------------------------------------------

      jcasc = 4      ! S.Abe 2019/03/08, JQMD flag for therdec

  100 continue
      if(irqmd .eq. 0) then
       call jqmdin( 14, 22, erg, ita, itz1, bmax0)
      else
       call jqmdinR( 14, 22, erg, ita, itz1, bmax0)
      endif

cABE 2015/08/10
      call jqmdchk(excit0,ita,itz2)
      if( nclst .lt. 1 ) then
       ipim = ipim + 1
       if( ipim .le. 100 ) goto 100
      endif

*-----------------------------------------------------------------------
*     mu_neutrino and excited nucleus
*-----------------------------------------------------------------------

      excit = excit0       ! [MeV]
      theta = thetacap - pi
      phi = phicap      ! S.Abe 2015/06/24 corrected

      costhe = dcos(theta)
      sinthe = dsin(theta)
      cosphi = dcos(phi)
      sinphi = dsin(phi)

      ptmn = ((rmumas+rtamas1-emu1s)**2.d0-(rtamas2+excit)**2.d0)
     &       * 0.5d0 / (rmumas+rtamas1-emu1s)
      pmn(5) = 0.d0
      pmn(1) = ptmn * sinthe * cosphi
      pmn(2) = ptmn * sinthe * sinphi
      pmn(3) = ptmn * costhe
      pmn(4) = ptmn

      pnucl(5) = rtamas2 + excit
      pnucl(1) = -pmn(1)
      pnucl(2) = -pmn(2)
      pnucl(3) = -pmn(3)
      ptnucl = dsqrt(pnucl(1)**2.d0 + pnucl(2)**2.d0 + pnucl(3)**2.d0)
      pnucl(4) = dsqrt(ptnucl**2.d0 + pnucl(5)**2.d0)

      do i = 1, 5
       pmn(i) = pmn(i) * 1.d-3      ! [GeV]
       pnucl(i) = pnucl(i) * 1.d-3  ! [GeV]
      enddo

*-----------------------------------------------------------------------
*     transform particles to excited nucleus flame
*-----------------------------------------------------------------------

cABE 2015/06/24, modified
      if( nclst .ge. 1 ) then

       pxb = pnucl(1)
       pyb = pnucl(2)
       pzb = pnucl(3)

       ptb = dsqrt(pxb**2.d0 + pyb**2.d0 + pzb**2.d0)
       etb = pnucl(4)

       betl = ptb / etb
       gaml = 1.d0 / dsqrt(1.d0 - betl**2.d0)

       if( ptb .gt. 0.d0 ) then

        costhe = pzb / ptb
        if( costhe.eq.1.d0 .or. costhe.eq.-1.d0) then
         sinthe = 0.d0
         cosphi = 1.d0
         sinphi = 0.d0
        else
         sinthe = dsqrt(1.d0-costhe**2.d0)
         cosphi = pxb / ptb / sinthe
         sinphi = pyb / ptb / sinthe
         if( cosphi .gt. 1.0d0 .or. cosphi .lt. -1.0d0 ) then
          phi = dasin(sinphi)
          cosphi = dcos(phi)
         endif
         if( sinphi .gt. 1.0d0 .or. sinphi .lt. -1.0d0 ) then
          phi = dacos(cosphi)
          sinphi = dsin(phi)
         endif
        endif

        do i = 1, nclst

         pini(1) = qclust(1,i)
         pini(2) = qclust(2,i)
         pini(3) = qclust(3,i)
         eini = qclust(4,i)

         pirot(1) =  pini(1) * costhe * cosphi
     &             + pini(2) * costhe * sinphi
     &              -pini(3) * sinthe
         pirot(2) = -pini(1) * sinphi
     &             + pini(2) * cosphi
         pirot(3) =  pini(1) * sinthe * cosphi
     &             + pini(2) * sinthe * sinphi
     &             + pini(3) * costhe

         pfrot(1) = pirot(1)
         pfrot(2) = pirot(2)
         pfrot(3) = pirot(3) * gaml + eini * gaml * betl

         pfin(1) =  pfrot(1) * costhe * cosphi
     &             -pfrot(2) * sinphi
     &            + pfrot(3) * sinthe * cosphi
         pfin(2) =  pfrot(1) * costhe * sinphi
     &            + pfrot(2) * cosphi
     &            + pfrot(3) * sinthe * sinphi
         pfin(3) = -pfrot(1) * sinthe
     &            + pfrot(3) * costhe

         qclust(1,i) = pfin(1)
         qclust(2,i) = pfin(2)
         qclust(3,i) = pfin(3)
         qclust(4,i) = dsqrt(qclust(1,i)**2.d0 + qclust(2,i)**2.d0
     &                       + qclust(3,i)**2.d0 + qclust(5,i)**2.d0)
         qclust(7,i) = (qclust(4,i) - qclust(5,i)) * 1.d3

        enddo

       endif

*-----------------------------------------------------------------------
*     set information of excited nucleus when QMD calc. failed
*-----------------------------------------------------------------------

      else

       elab = (pnucl(4)-pnucl(5)) * 1.d+3

       nclst = 1

       iclust(nclst) = 0

       jclust(0,nclst) = 0
       jclust(1,nclst) = itz2
       jclust(2,nclst) = itn2
       jclust(3,nclst) = 19
       jclust(4,nclst) = 0
       jclust(5,nclst) = itz2
       jclust(6,nclst) = itz2 + itn2
       jclust(7,nclst) = itz2 * 1000000 + itz2 + itn2
       jclust(8,nclst) = 0

       qclust(0,nclst)  = 0.d0
       qclust(1,nclst)  = pnucl(1)
       qclust(2,nclst)  = pnucl(2)
       qclust(3,nclst)  = pnucl(3)
       qclust(4,nclst)  = pnucl(4)
       qclust(5,nclst)  = pnucl(5)
       qclust(6,nclst)  = excit
       qclust(7,nclst)  = elab
       qclust(8,nclst)  = wgt / wg0
       qclust(9,nclst)  = 0.d0
       qclust(10,nclst) = 0.d0
       qclust(11,nclst) = 0.d0
       qclust(12,nclst) = 0.d0

      endif

*-----------------------------------------------------------------------

      call nevap(0)

*-----------------------------------------------------------------------
*     store information of mu_neutrino
*-----------------------------------------------------------------------

      nclsts = nclsts + 1
      iclusts(nclsts) = 7     ! S.Abe 2015/06/18 corrected

      jclusts(0,nclsts) = 0
      jclusts(1,nclsts) = 0
      jclusts(2,nclsts) = 0
      jclusts(3,nclsts) = 11
      jclusts(4,nclsts) = 0
      jclusts(5,nclsts) = 0
      jclusts(6,nclsts) = 0
      jclusts(7,nclsts) = 14
      jclusts(8,nclsts) = 0

      qclusts(0,nclsts) = 0.d0
      qclusts(1,nclsts) = pmn(1)
      qclusts(2,nclsts) = pmn(2)
      qclusts(3,nclsts) = pmn(3)
      qclusts(4,nclsts) = pmn(4)
      qclusts(5,nclsts) = pmn(5)
      qclusts(6,nclsts) = 0.d0
      qclusts(7,nclsts) = pmn(4) * 1.d3
      qclusts(8,nclsts) = wgt / wg0
      qclusts(9,nclsts) = 0.d0
      qclusts(10,nclsts) = 0.d0
      qclusts(11,nclsts) = 0.d0
      qclusts(12,nclsts) = 0.d0

      numpal(11) = numpal(11) + 1
      rumpal(11) = rumpal(11) + wgt

*-----------------------------------------------------------------------
*     composition counter momentum of characteristic X-ray
*-----------------------------------------------------------------------

      pxb = -pxsum
      pyb = -pysum
      pzb = -pzsum

      ptb = dsqrt(pxb**2.d0 + pyb**2.d0 + pzb**2.d0)
      etb = rtamas1 + rmumas - emu1s

      betl = ptb / etb      ! S.Abe 2015/06/24 corrected
      gaml = 1.d0 / dsqrt(1.d0 - betl**2.d0)

      if( ptb .gt. 0.d0 ) then

       costhe = pzb / ptb
       if(costhe .eq. 1.0d0 .or. costhe .eq. -1.0d0) then
        sinthe = 0.d0
        cosphi = 1.d0
        sinphi = 0.d0
       else
        sinthe = dsqrt(1.d0-costhe**2.d0)
        cosphi = pxb / ptb / sinthe    ! S.Abe 2015/06/18 corrected
        sinphi = pyb / ptb / sinthe
        if( cosphi .gt. 1.0d0 .or. cosphi .lt. -1.0d0 ) then
         phi = dasin(sinphi)
         cosphi = dcos(phi)
        endif
        if( sinphi .gt. 1.0d0 .or. sinphi .lt. -1.0d0 ) then
         phi = dacos(cosphi)
         sinphi = dsin(phi)
        endif
       endif

cABE 2015/06/24 corrected
       do i = nxray+1, nclsts

        pini(1) = qclusts(1,i)
        pini(2) = qclusts(2,i)
        pini(3) = qclusts(3,i)
        eini = qclusts(4,i)

        pirot(1) =  pini(1) * costhe * cosphi
     &            + pini(2) * costhe * sinphi
     &             -pini(3) * sinthe
        pirot(2) = -pini(1) * sinphi
     &            + pini(2) * cosphi
        pirot(3) =  pini(1) * sinthe * cosphi
     &            + pini(2) * sinthe * sinphi
     &            + pini(3) * costhe

        pfrot(1) = pirot(1)
        pfrot(2) = pirot(2)
        pfrot(3) = pirot(3) * gaml + eini * gaml * betl

        pfin(1) =  pfrot(1) * costhe * cosphi
     &            -pfrot(2) * sinphi
     &           + pfrot(3) * sinthe * cosphi
        pfin(2) =  pfrot(1) * costhe * sinphi
     &           + pfrot(2) * cosphi
     &           + pfrot(3) * sinthe * sinphi
        pfin(3) = -pfrot(1) * sinthe
     &           + pfrot(3) * costhe

        qclusts(1,i) = pfin(1)
        qclusts(2,i) = pfin(2)
        qclusts(3,i) = pfin(3)
        qclusts(4,i) = dsqrt(qclusts(1,i)**2.d0 + qclusts(2,i)**2.d0
     &                       + qclusts(3,i)**2.d0 + qclusts(5,i)**2.d0)
        qclusts(7,i) = (qclusts(4,i) - qclusts(5,i)) * 1.d3

       enddo

      endif

*-----------------------------------------------------------------------
*     reset muon capture flag
*-----------------------------------------------------------------------

      imucaphit = 0

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      subroutine pnconvert_qmd
*                                                                      *
*      select neutron converted from proton in the target nucleus      *
*      for QMD calculation                                             *
*                                                                      *
*      last modified by S.Abe on 2015/06/24                            *
*                                                                      *
************************************************************************
      use QMD_COOD2_MOD, only : r0, p0
      implicit real*8(a-h,o-z)

      include 'param00.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /input1/ mstq1(mxpa1), parq1(mxpa1)
!$OMP THREADPRIVATE(/input1/)

      common /vriab0/ massal, massba, mmeson
!$OMP THREADPRIVATE(/vriab0/)

      common /qmdflg/ irqmd

      common /coodrp/ r(5,nnn), p(6,nnn)
!$OMP THREADPRIVATE(/coodrp/)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)
!$OMP THREADPRIVATE(/coodid/)

      common /capmudir/ thetacap, phicap
!$OMP THREADPRIVATE(/capmudir/)

*-----------------------------------------------------------------------


cABE 2015/06/24 changed
      dimension pini(3)
      dimension pirot(3)
      dimension pfrot(3)
      dimension pfin(3)

      parameter ( pi = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      r0 = 0.d0
      p0 = 0.d0
      itry1 = 0
      itry2 = 0
      itry3 = 0

*-----------------------------------------------------------------------
*     store initial information for r(i,id) and p(i,id)
*-----------------------------------------------------------------------

      do id = 1, massal
       do i = 1, 5
        r0(i,id) = r(i,id)
       enddo
       do i = 1, 6
        p0(i,id) = p(i,id)
       enddo
      enddo

*-----------------------------------------------------------------------
*     set excitation energy
*-----------------------------------------------------------------------

  100 continue                         ! retry point
      itry1 = itry1 + 1

cABE 2015/08/11 changed
c      call excit_mucap_gauss(excit)
      call excit_mucap_amado(excit)

*-----------------------------------------------------------------------
*     collision in CM frame
*-----------------------------------------------------------------------

      idnta  = mstq1(4)                ! mstq1(4) = 0
      massta = mstq1(5)
      mstapr = mstq1(6)

      prmas  = ulmass( 0 )       
      tamas  = ulmass( idnta )         ! ulmass(0) = rmass = 0.9383 GeV/c2

      n1 = 1
      n2 = massta - 1

*-----------------------------------------------------------------------
 
      elab = excit * 1.d-3
      plab = dsqrt( elab * ( 2.d0 * prmas + elab ) )

      ptot = plab * n1
      etot = elab * n1 + prmas * n1 + tamas * n2
      stot = dsqrt( etot**2.d0 - ptot**2.d0 )
      pstt = pcmsr( stot, prmas*n1, tamas*n2 )

      p1 = pstt / n1
      e1 = dsqrt( prmas**2.d0 + p1**2.d0 )
      vzpr = p1 / e1
      gampr = e1 / prmas

      if( n2 .eq. 0 ) then
       vzta = 0.d0
       gamta = 0.d0
      else
       p2 = -pstt / n2
       e2 = dsqrt( tamas**2.d0 + p2**2.d0 )
       vzta = p2 / e2
       gamta = e2 / tamas
      endif

*-----------------------------------------------------------------------
*     select proton which capture negative muon
*-----------------------------------------------------------------------

  200 continue                         ! retry point
      itry2 = itry2 + 1
      if(itry2 .gt. 10)then
       itry2 = 0
       goto 100
      endif

c...randomly select
  210 icapid = int( unirn(dummy)* massal) + 1
      if(ichg(icapid) .eq. 0) goto 210

c...select outermost proton
c      distmax = 0.d0
c      do i = 1, massal
c       dist = dsqrt(r(1,i)**2.d0 + r(2,i)**2.d0 + r(3,i)**2.d0)
c       if( ichg(i) .eq. 1 .and. dist .ge. distmax) then
c        distmax = dist
c        icapid = i
c       endif
c      enddo

*-----------------------------------------------------------------------
*     add momentum for each nucleon like rboost
*-----------------------------------------------------------------------

  300 continue                         ! retry point
      itry3 = itry3 + 1
      if(itry3 .gt. 10)then
       itry3 = 0
       goto 200
      endif

      thetacap = pi * unirn(dummy)
      phicap = 2.d0 * pi * unirn(dummy)

      costhe = dcos(thetacap)
      sinthe = dsin(thetacap)
      cosphi = dcos(phicap)
      sinphi = dsin(phicap)

*-----------------------------------------------------------------------

cABE 2015/06/24 corrected
      do i = 1, massal

       if( i .eq. icapid ) then
        betl  = vzpr
        gaml  = gampr
       else
        betl  = vzta
        gaml  = gamta
       end if

       pini(1) = p(1,i)
       pini(2) = p(2,i)
       pini(3) = p(3,i)

       if( irqmd .eq. 0 ) then
        eini = p(4,i)
       else
        eini = p(6,i)
        epotdl = (p(4,i)**2.d0 - p(5,i)**2.d0
     &            - p(1,i)**2.d0 - p(2,i)**2.d0
     &            - p(3,i)**2.d0) * 0.5d0 / p(5,i)
       endif

       pirot(1) =  pini(1) * costhe * cosphi
     &           + pini(2) * costhe * sinphi
     &            -pini(3) * sinthe
       pirot(2) = -pini(1) * sinphi
     &           + pini(2) * cosphi
       pirot(3) =  pini(1) * sinthe * cosphi
     &           + pini(2) * sinthe * sinphi
     &           + pini(3) * costhe

       pfrot(1) = pirot(1)
       pfrot(2) = pirot(2)
       pfrot(3) = pirot(3) * gaml + eini * gaml * betl

       pfin(1) =  pfrot(1) * costhe * cosphi
     &           -pfrot(2) * sinphi
     &          + pfrot(3) * sinthe * cosphi
       pfin(2) =  pfrot(1) * costhe * sinphi
     &          + pfrot(2) * cosphi
     &          + pfrot(3) * sinthe * sinphi
       pfin(3) = -pfrot(1) * sinthe
     &          + pfrot(3) * costhe

       p(1,i) = pfin(1)
       p(2,i) = pfin(2)
       p(3,i) = pfin(3)

       if(irqmd .eq. 0) then
        p(4,i) = dsqrt( p(5,i)**2 + p(1,i)**2 + p(2,i)**2 + p(3,i)**2)
       else
        p(4,i) = dsqrt(p(1,i)**2.d0 + p(2,i)**2.d0
     &                 + p(3,i)**2.d0 + p(5,i)**2.d0
     &                 + 2.d0 * p(5,i) * epotdl)
        p(6,i) = dsqrt(p(1,i)**2.d0 + p(2,i)**2.d0
     &                 + p(3,i)**2.d0 + p(5,i)**2.d0)
       endif

      end do

*-----------------------------------------------------------------------
*     convert proton to neutron and check pauli blocking
*-----------------------------------------------------------------------

      ichg(icapid) = 0

      if(irqmd .eq. 0) then
       call caldisa
       call pauli( icapid, ntag, phase )
      else
       call caldisaR
       call pauliR( icapid, ntag, phase )
       do i = 1, massal
        call epotprtR(i,epotp)
        p(4,i) = sqrt( p(1,i)**2+ p(2,i)**2 + p(3,i)**2
     &                + 2 * p(5,i) * epotp + p(5,i)**2 )
       end do
      endif

      if(ntag .eq. 1) then    ! return before convert

       do id = 1, massal
        do i = 1, 5
         r(i,id) = r0(i,id)
        enddo
        do i = 1, 6
         p(i,id) = p0(i,id)
        enddo
       enddo

       ichg(icapid) = 1
       goto 300

      endif

*-----------------------------------------------------------------------
cABE add @2015/04/10, for run "relcol" ... =1: target, =-1:projectile
      iavd(icapid) = -1

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      subroutine excit_mucap_amado(excit)
*                                                                      *
*      excitation function for muon capture                            *
*      calculated by using amado momentum distribution                 *
*                                                                      *
*      made by S.Abe on 2015/04/16                                     *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     output:                                                          *
*        excit : excitation energy of nucleus [MeV]                    *
*                                                                      *
*---- Reference -------------------------------------------------------*
*                                                                      *
*     [1] M.Lifshitz, P.Singer, Phys.Rev.C 22, 2135 (1980)             *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /mueng/ emu1s, rtamas1, rtamas2, excit0
!$OMP THREADPRIVATE(/mueng/)

      dimension exd(1001)
      dimension exfd(1001)

      save e0totb
!$OMP THREADPRIVATE(e0totb)
      save exd
!$OMP THREADPRIVATE(exd)
      save exfd
!$OMP THREADPRIVATE(exfd)

      parameter ( rnmass = 939.565379d0 )
      parameter ( rmumas = 105.6583715d0 )
      parameter ( hbarc  = 197.327053d0 )

*-----------------------------------------------------------------------
*     set constant value
*-----------------------------------------------------------------------

      rmstar = 0.68d0 * rnmass            ! M* = 0.68 * M

      gamma = 0.8d0                       ! [fm]
      gamma = gamma / hbarc               ! [1/MeV]

      e0tot = rtamas1 - rtamas2 + rmumas - emu1s

*-----------------------------------------------------------------------
*     defined integration bin of excitation energy and momentum
*-----------------------------------------------------------------------

      excmin = 0.d0
      excmax = e0tot
      nexcit = 100
      delte = (excmax - excmin) / dble(nexcit)

      np = 10000
      deltp = 10.0d-0

*-----------------------------------------------------------------------
*     skip exfd calc. when e0tot is same as before
*-----------------------------------------------------------------------

      if( e0totb .eq. e0tot ) then
       goto 1000
      else
       e0totb = e0tot
      endif

*-----------------------------------------------------------------------
*     calculation of excitation function
*-----------------------------------------------------------------------

      exd(1) = 0.d0
      totexf = 0.d0

      do i = 1, nexcit-1

       exc = delte * dble(i)
       eresid = excmax - exc
       p0 = 0.5d0 * dabs( eresid - 2.d0 * rmstar * exc / eresid )
       fmax = 0.d0
       finteg = 0.d0

*-----------------------------------------------------------------------

       pr = p0
       if( pr .le. 8.d0/gamma) then
        prfunc = 1.d0 / dcosh(gamma*pr)**2.d0
       else
        prfunc = 4.d0 * dexp(-2.d0*gamma*pr)
       endif

       qr = dsqrt(pr**2.d0 + 2.d0*rmstar*exc)
       if( qr .le. 8.d0/gamma) then
        qrfunc = 1.d0 / dcosh(gamma*qr)**2.d0
       else
        qrfunc = 4.d0 * dexp(-2.d0*gamma*qr)
       endif

       fr = pr * prfunc * (1.d0 - qrfunc)

       if( fr .le. 0.d0) exit

*-----------------------------------------------------------------------
*      integration by Simpson's rule
*-----------------------------------------------------------------------

       do j = 1, np
        pl = pr
        ql = dsqrt(pl**2.d0 + 2.d0*rmstar*exc)
        fl = fr

        pr = pl + deltp
        if( pr .le. 8.d0/gamma) then
         prfunc = 1.d0 / dcosh(gamma*pr)**2.d0
        else
         prfunc = 4.d0 * dexp(-2.d0*gamma*pr)
        endif
        qr = dsqrt(pr**2.d0 + 2.d0*rmstar*exc)
        if( qr .le. 8.d0/gamma) then
         qrfunc = 1.d0 / dcosh(gamma*qr)**2.d0
        else
         qrfunc = 4.d0 * dexp(-2.d0*gamma*qr)
        endif
        fr = pr * prfunc * (1.d0 - qrfunc)

        pc = (pl + pr) * 0.5d0
        if( pc .le. 8.d0/gamma) then
         pcfunc = 1.d0 / dcosh(gamma*pc)**2.d0
        else
         pcfunc = 4.d0 * dexp(-2.d0*gamma*pc)
        endif
        qc = dsqrt(pc**2.d0 + 2.d0*rmstar*exc)
        if( qc .le. 8.d0/gamma) then
         qcfunc = 1.d0 / dcosh(gamma*qc)**2.d0
        else
         qcfunc = 4.d0 * dexp(-2.d0*gamma*qc)
        endif
        fc = pc * pcfunc * (1.d0 - qcfunc)

        dpbin = pr - pl

        finteg = finteg +
     &           (fl + 4.0d0*fc + fr) / 6.0d0 * dpbin

        fmax = dmax1(fmax,fr)
        if( fr/fmax .lt. 1.d-3) exit

       enddo

*-----------------------------------------------------------------------

       exfunc = eresid * finteg

       if( exfunc .gt. 1.d-9 ) then
        totexf = totexf + exfunc
        exd(i+1) = exc
        exfd(i) = exfunc
       endif

      enddo

      do i = 1, nexcit
       exfd(i) = exfd(i) / totexf
      enddo

*-----------------------------------------------------------------------
*     sampling of excitation energy from excitation function
*-----------------------------------------------------------------------

 1000 continue

      exselect = unirn(dummy)

      do i = 1, nexcit

       exselect = exselect - exfd(i)

       if( exselect .le. 0.d0 ) then
        exselect = exselect + exfd(i)
        excit = (exd(i+1)-exd(i)) * exselect/exfd(i) + exd(i)
        excit0 = excit     ! S.Abe 2015/06/18 added
        exit
       endif

      enddo

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      subroutine excit_mucap_gauss(excit)
*                                                                      *
*      select excitation energy by capturering muon                    *
*                                                                      *
*      last modified by S.Abe on 2015/06/24                            *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     output:                                                          *
*        excit : excitation energy of nucleus [MeV]                    *
*                                                                      *
*---- Reference -------------------------------------------------------*
*                                                                      *
*     [1] P.Singer,Il Nuovo Cimento,23,669(1962)                       *
*     [2] M.Lifshitz+,Nucl.Phys.A,476,684(1988)                        *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      common /mueng/ emu1s, rtamas1, rtamas2, excit0
!$OMP THREADPRIVATE(/mueng/)

      dimension exd(1001)
      dimension exfd(1001)
      dimension exfdg(1001)

      save e0totb
!$OMP THREADPRIVATE(e0totb)
      save exd
!$OMP THREADPRIVATE(exd)
      save exfd
!$OMP THREADPRIVATE(exfd)
      save totexf
!$OMP THREADPRIVATE(totexf)

      parameter ( pi = 3.141592653589793d0 )
      parameter ( rpmass = 938.272046d0 )
      parameter ( rnmass = 939.565379d0 )
      parameter ( rmumas = 105.6583715d0 )

*-----------------------------------------------------------------------
*     set constant value
*-----------------------------------------------------------------------

      rmstar = 0.68d0 * rnmass           ! M* = 0.68 * M [MeV]
      awidth2 = 20.d0 * 2.d0 * rnmass    ! a2/2M = 20.0 [MeV]

      e0tot = rtamas1 - rtamas2 + rmumas - emu1s

*-----------------------------------------------------------------------
*     define integration bin of excitation energy and momentum
*-----------------------------------------------------------------------

      excmin = 0.d0
      excmax = e0tot
      nexcit = 100
      delext = (excmax - excmin) / dble(nexcit)

*-----------------------------------------------------------------------
*     skip exfd calc. when e0tot is same as before
*-----------------------------------------------------------------------

      if( e0totb .eq. e0tot ) then
       goto 1000
      else
       e0totb = e0tot
      endif

*-----------------------------------------------------------------------
*     calculation of excitation function from Ref.[1]
*-----------------------------------------------------------------------

      exd(1) = 0.d0
      totexf = 0.d0
c...eq.(42)
      aaa = 0.25d0 / awidth2 * e0tot**2.d0
      bbb = rmstar**2.d0 / awidth2
      ccc = e0tot * rmstar / awidth2

      do i = 1, nexcit-1
       qexcit = delext * dble(i)
c...eq.(41)
       func2 = aaa * (e0tot-qexcit)**2.d0 / e0tot**2.d0
     &         + bbb * qexcit**2.d0 / (e0tot-qexcit)**2.d0
       func1 = func2 - ccc * qexcit / e0tot
c...eq.(40)
       exfunc = (e0tot-qexcit) * dexp(-func1)
     &             - 0.5d0 * (e0tot-qexcit) * dexp(-2.d0*func2)

       if( exfunc .gt. 1.d-9 ) then
        totexf = totexf + exfunc
        exd(i+1) = qexcit
        exfd(i) = exfunc
       else
        exd(i+1) = qexcit
        exit
       endif
      enddo

      goto 1000

*-----------------------------------------------------------------------
*     add excitation function for meson exchange currents from Ref.[2]
*     not realistic, but just gaussian fitting
*-----------------------------------------------------------------------

      totexfg = 0.d0

      ccc = 16.085d+0
      bbb = 54.127d+0
      aaa = 1.05d+3 / dsqrt(2.d0*pi) / ccc

      do i = 1, nexcit-1
       qexcit = delext * dble(i)

       exfuncg = aaa * dexp(-1.d0*(qexcit-bbb)**2.d0 / (2.d0*ccc))

       totexfg = totexfg + exfuncg
       exfdg(i) = exfuncg
      enddo

      ratiog = 0.06d0 * totexf / totexfg

      totexf = 0.d0

      do i = 1, nexcit-1
       exfd(i) = exfd(i) + exfdg(i)*ratiog
       totexf = totexf + exfd(i)
      enddo

*-----------------------------------------------------------------------
*     sampling transferred energy into neutron
*-----------------------------------------------------------------------

 1000 continue

      exselect = unirn(dummy) * totexf

      do i = 1, nexcit
       exselect = exselect - exfd(i)
       if( exselect .le. 0.d0 ) then
        exselect = exselect + exfd(i)
        excit = (exd(i+1)-exd(i))*exselect/exfd(i) + exd(i)
        excit0 = excit
        exit
       endif
      enddo

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      function caprmu(itz,ita)
*                                                                      *
*     caprmu: negative muon capture rate [ms-1]                        *
*     if data does not exist, Goulard-Primakoff formula [1] is used    *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input :                                                          *
*        itz : atomic number of target                                 *
*        ita : mass number of target                                   *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*     mcrdat :                                                         *
*        mcrkf(15) : muon capture rate available kf code               *
*        crdat(15) : capture rate exp. data [ms-1]                     *
*        zeffd(iz) : effective atomic number                           *
*        gconst(2,4) : parameter set for equation of ca@tire rate      *
*                                                                      *
*---- Reference -------------------------------------------------------*
*                                                                      *
*     [1] B.Goulard+,Phys.Rev.C,10,2034(1974)                          *
*     [2] J.A.Wheeler,Rev.Mod.Phys.,31,133(1949)                       *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'param00.inc'

      common /mcrdat/ crdat(15), zeffd(118), gconst(2,4), mcrkf(15)

*-----------------------------------------------------------------------

      caprmu = 0.d0

      if( itz .eq. 0 .and. ita .eq. 1) then
       itkf = 2112
      elseif( itz .eq. 1 .and. ita .eq. 1) then
       itkf = 2212
      else
       itkf = itz * 1000000 + ita
      endif

      zt = dble(itz)
      at = dble(ita)

*-----------------------------------------------------------------------

      if( itz .lt. 9 ) then
       do i = 1, 15
        if( itkf .eq. mcrkf(i) ) then
         caprmu = crdat(i)
         exit
        endif
       enddo
      endif

      if ( caprmu .eq. 0.d0 ) then

       if ( zeffd(itz) .gt. 0.d0 ) then
        zeff = zeffd(itz)
       else
! semiempirical formula from Ref.[2]
        zeff = zt * (1.d0 + (zt/37.2d0)**1.54d0)**(-1.d0/1.54d0)
       endif

       icset = 2      ! use parameter set of "TRIUMF data"

       caprmu = 1.d-3 * zeff**4.d0 * gconst(icset,1) *
     &          ( 1.d0 + gconst(icset,2) * at/(2.d0*zt) -
     &            gconst(icset,3) * (at-2.d0*zt)/(2.d0*zt) -
     &            gconst(icset,4) * ((at-zt)/(2.d0*at) +
     &                               (at-2.d0*zt)/(8.d0*at*zt)) )
      endif

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      function dcyrmu(itz)
*                                                                      *
*     dcyrmu: decay rate [ms-1] calculated by Huff factor              *
*     if data does not exist, fitting function using in [1] is used.   *
*                                                                      *
*---- in argument -----------------------------------------------------*
*                                                                      *
*     input :                                                          *
*        itz : target atomic number                                    *
*        ita : target mass number                                      *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*     mdrdat :                                                         *
*        huffd(iz) : huff factor                                       *
*                                                                      *
*---- Reference -------------------------------------------------------*
*                                                                      *
*     [1] M.V.Kossov,Eur.Phys.J.A,33,7(2007)                           *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      common /mdrdat/ huffd(118)

      parameter ( dcytfm = 2.19703d-6 )

*-----------------------------------------------------------------------

      if( huffd(itz) .gt. 0.d0 ) then
       huff = huffd(itz)
      else
       zt = dble(itz)
       huff = 1.d0 - 3.94d-4 * zt**2.19d+0 /
     &               (1.d+0 + 1.218d+1 * exp(1.373d-2*zt) )
      endif

      dcyrmu = 1.d-3 / dcytfm * huff

*-----------------------------------------------------------------------
      return
      end

************************************************************************
*                                                                      *
      block data mcrate
*                                                                      *
*     experimental capture rate [ms-1] for light materials (Z<10)      *
*     taken from below reference                                       *
*                                                                      *
*      kf =    2212: [1] G.Bardin+,Nucl.Phys.A,352,365(1981)           *
*           1000002: [2] G.Bardin+,Nucl.Phys.A,453,591(1986)           *
*           2000003: [3] I.V.Falomkin+,Phys.Lett.,1,318(1961)          *
*           2000004: [4] M.M.Block,Il Nuovo Cimento,55,501(1968)       *
*           3000006: [5] T.Suzuki,Phys.Rev.C,35,2212(1987)             *
*           3000007: [5]                                               *
*           4000009: [5]                                               *
*           5000010: [5]                                               *
*           5000011: [5]                                               *
*           6000012: [5]                                               *
*           6000013: [5]                                               *
*           7000014: [5]                                               *
*           8000016: [5]                                               *
*           8000018: [5]                                               *
*           9000019: [5]                                               *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*     mcrdat :                                                         *
*        mcrkf(15) : muon capture rate available kf code               *
*        crdat(15) : capture rate exp. data [ms-1]                     *
*        zeffd(iz) : effective atomic number,                          *
*                    not available data is set -1.0                    *
*        gconst(2,4) : parameter set for equation of capture rate      *
*                                                                      *
*     mdrdat :                                                         *
*        huffd(iz) : huff factor,                                      *
*                    not available data is set -1.0                    *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /mcrdat/ crdat(15), zeffd(118), gconst(2,4), mcrkf(15)
      common /mdrdat/ huffd(118)

*-----------------------------------------------------------------------

      data (mcrkf(i),i=1,15) /
     &        2212,   1000002,   2000003,   2000004,   3000006,
     &     3000007,   4000009,   5000010,   5000011,   6000012,
     &     6000013,   7000014,   8000016,   8000018,   9000019/

      data (crdat(i),i=1,15 ) /
     &    4.200d-1,  4.700d-1,  2.140d+0,  3.640d-1,  4.180d+0,
     &    1.810d+0,  7.400d+0,  2.780d+1,  2.190d+1,  3.880d+1,
     &    3.760d+1,  6.930d+1,  1.026d+2,  8.800d+1,  2.290d+2/

*-----------------------------------------------------------------------

      data ((gconst(i,j),j=1,4),i=1,2 ) /
     &    2.520d+2, -3.800d-2, -2.400d-1,  3.230d+0,                    ! "pre-TRIUMF data"
     &    2.610d+2, -4.000d-2, -2.600d-1,  3.240d+0 /                   ! "TRIUMF data"

      data (zeffd(iz),iz=1,118 ) /
*                1          2          3          4          5
     &    1.000d+0,  1.980d+0,  2.940d+0,  3.890d+0,  4.810d+0,
     &    5.720d+0,  6.610d+0,  7.490d+0,  8.320d+0,  9.140d+0,
     1    9.950d+0,  1.069d+1,  1.148d+1,  1.222d+1,  1.290d+1,
     &    1.364d+1,  1.424d+1,  1.489d+1,  1.553d+1,  1.615d+1,
     2    1.677d+1,  1.738d+1,  1.804d+1,  1.849d+1,  1.906d+1,
     &    1.959d+1,  2.013d+1,  2.066d+1,  2.112d+1,  2.161d+1,
     3    2.202d+1,  2.243d+1,  2.284d+1,  2.324d+1,  2.365d+1,
     &   -1.000d+0,  2.447d+1,  2.485d+1,  2.523d+1,  2.561d+1,
     4    2.599d+1,  2.637d+1, -1.000d+0, -1.000d+0,  2.732d+1,
     &    2.763d+1,  2.795d+1,  2.820d+1,  2.842d+1,  2.864d+1,
     5    2.879d+1,  2.903d+1,  2.927d+1, -1.000d+0,  2.975d+1,
     &    2.999d+1,  3.022d+1,  3.036d+1,  3.053d+1,  3.069d+1,
     6   -1.000d+0,  3.101d+1, -1.000d+0,  3.134d+1,  3.148d+1,
     &    3.162d+1,  3.176d+1,  3.190d+1, -1.000d+0, -1.000d+0,
     7   -1.000d+0,  3.247d+1,  3.261d+1,  3.276d+1, -1.000d+0,
     &   -1.000d+0, -1.000d+0, -1.000d+0,  3.364d+1,  3.381d+1,
     8    3.421d+1,  3.418d+1,  3.400d+1, -1.000d+0, -1.000d+0,
     &   -1.000d+0, -1.000d+0, -1.000d+0, -1.000d+0,  3.473d+1,
     9   -1.000d+0,  3.494d+1,  3.505d+1,  3.516d+1, -1.000d+0,
     &   -1.000d+0, -1.000d+0, -1.000d+0, -1.000d+0, -1.000d+0,
     &   -1.000d+0, -1.000d+0, -1.000d+0, -1.000d+0, -1.000d+0,
     &   -1.000d+0, -1.000d+0, -1.000d+0, -1.000d+0, -1.000d+0,
     1   -1.000d+0, -1.000d+0, -1.000d+0, -1.000d+0, -1.000d+0,
     &   -1.000d+0, -1.000d+0, -1.000d+0/

      data (huffd(iz),iz=1,118 ) /
*                1          2          3          4          5
     &    1.000d+0,  1.000d+0,  1.000d+0,  1.000d+0,  1.000d+0,
     &    1.000d+0,  1.000d+0,  0.998d+0,  0.998d+0,  0.997d+0,
     1    0.996d+0,  0.995d+0,  0.993d+0,  0.992d+0,  0.991d+0,
     &    0.990d+0,  0.989d+0,  0.988d+0,  0.987d+0,  0.985d+0,
     2    0.983d+0,  0.981d+0,  0.980d+0,  0.978d+0,  0.976d+0,
     &    0.975d+0,  0.971d+0,  0.969d+0,  0.967d+0,  0.965d+0,
     3    0.962d+0,  0.960d+0,  0.958d+0,  0.955d+0,  0.952d+0,
     &   -1.000d+0,  0.948d+0,  0.945d+0,  0.942d+0,  0.940d+0,
     4    0.939d+0,  0.936d+0, -1.000d+0, -1.000d+0,  0.929d+0,
     &    0.927d+0,  0.925d+0,  0.921d+0,  0.920d+0,  0.918d+0,
     5    0.916d+0,  0.913d+0,  0.910d+0, -1.000d+0,  0.905d+0,
     &    0.902d+0,  0.901d+0,  0.899d+0,  0.897d+0,  0.895d+0,
     6   -1.000d+0,  0.890d+0, -1.000d+0,  0.885d+0,  0.882d+0,
     &    0.880d+0,  0.877d+0,  0.875d+0, -1.000d+0, -1.000d+0,
     7   -1.000d+0,  0.865d+0,  0.862d+0,  0.860d+0, -1.000d+0,
     &   -1.000d+0, -1.000d+0, -1.000d+0,  0.850d+0,  0.848d+0,
     8    0.846d+0,  0.844d+0,  0.840d+0, -1.000d+0, -1.000d+0,
     &   -1.000d+0, -1.000d+0, -1.000d+0, -1.000d+0,  0.824d+0,
     9   -1.000d+0,  0.820d+0,  0.818d+0,  0.816d+0, -1.000d+0,
     &   -1.000d+0, -1.000d+0, -1.000d+0, -1.000d+0, -1.000d+0,
     &   -1.000d+0, -1.000d+0, -1.000d+0, -1.000d+0, -1.000d+0,
     &   -1.000d+0, -1.000d+0, -1.000d+0, -1.000d+0, -1.000d+0,
     1   -1.000d+0, -1.000d+0, -1.000d+0, -1.000d+0, -1.000d+0,
     &   -1.000d+0, -1.000d+0, -1.000d+0/

      end
