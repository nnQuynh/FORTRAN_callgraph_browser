************************************************************************
*                                                                      *
      subroutine xstneu(mm,sigt,sigaa,icli,ein,tmei,mki)
*                                                                      *
*       calculate the neutron cross sections in material mk.           *
*       mm=1 for call from wtmult, otherwise zero.                     *
*       Last modified by K.Niita on 2016/08/08                         *
*                                                                      *
************************************************************************

      use MEMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      use GGMBANKMOD !FURUTAs
      use moddas_material

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      include 'err.inc'

*-----------------------------------------------------------------------

      common /elreg/  elarg(kvlmax,4), mlrgn, melrg(kvlmax,2), ielusr
      common /regdc/  idrg(kvlmax), idgr(kvmmax)

      common /eparm/  esmax, esmin, emin(20)
      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax)
      common /emode/  emodem, ge1, ge2, iemode
      common /sansd/  isans

*-----------------------------------------------------------------------
*     initialization of this material and particle
*-----------------------------------------------------------------------

         mk   = mki
         icl0 = icl

         erg  = ein
         tme  = tmei
         icl  = icli

         vel  = slite*sqrt(erg*(erg+2.*gpt(1)))/(erg+gpt(1))

         ipt  = 1
         jgp  = 0

         totm = 0.
         pfp  = 0.
         stp  = 0.
         deb  = huge
         sigt = 0.
         sigaa = 0.

*-----------------------------------------------------------------------
*     get the temperature ttn of cell icl.
*-----------------------------------------------------------------------

         if(mcal.ne.0) goto 310
         if(mm.ne.0) goto 50
         l = (icl-1)*mxt
         if(tme.ge.tth(mxt)) goto 40
         if(tme.le.tth(1)) goto 30

      do 10 n = 2, mxt
   10    if(tme.lt.tth(n)) goto 20

   20    ttn = tmp(l+n-1)+(tmp(l+n)-tmp(l+n-1))
     &       * (tme-tth(n-1))/(tth(n)-tth(n-1))
         goto 50
   30    ttn = tmp(l+1)
         goto 50
   40    ttn = tmp(l+mxt)

*-----------------------------------------------------------------------
*     calculate cross sections for each nuclide in the material.
*-----------------------------------------------------------------------

   50    totm = 0.
         siga = 0.
         lc = 0

         mii = 0

      do 270 m = jmd(1+mk), jmd(1+mk+1) - 1

*-----------------------------------------------------------------------

         mii = mii + 1
         siggc(mii) = 0.0d0
         siggcn(mii) = 0.d0
         siggce(mii) = 0.d0

         lem = m-jmd(1+mk)+1
         dsmax = das_kmate(kmate(mk)+(lem-1)*5+12)

      if( dsmax .gt. emin(2) ) then

*-----------------------------------------------------------------------

         iex = lme(1,m)

*-----------------------------------------------------------------------
*     find the energy bin in the cross-section table.
*-----------------------------------------------------------------------

         if(rtc(6,iex).eq.erg) goto 130
         n = nxs(3,iex)
         rtc(11,iex) = -1.
         if(nxs(16,iex).eq.lc) goto 100
         ic = jxs(1,iex)
         ib = ic+n-1
         if(erg.lt.rtc(6,iex)) ib = ic+ktc(1,iex)
         if(erg.gt.rtc(6,iex)) ic = ic+ktc(1,iex)-1

   60    if(ib-ic.eq.1) goto 80
         ih = (ic+ib)/2
         if(erg.lt.xss(ih)) goto 70
         ic = ih
         goto 60

   70    ib = ih
         goto 60

*-----------------------------------------------------------------------
*     get the bin index and interpolation fraction.
*-----------------------------------------------------------------------

   80    ktc(1,iex) = ic-jxs(1,iex)+1
         if(nty(iex).eq.2) goto 90
         rtc(1,iex) =
     &              max(zero,min(one,(erg-xss(ic))/(xss(ib)-xss(ic))))

*-----------------------------------------------------------------------
*     sample cross sections from unresolved range probability tables.
*-----------------------------------------------------------------------

         if(iunr.eq.0) goto 82
         l = jxs(23,iex)
         if(l.eq.0) goto 82
         if(erg.le.xss(l+6).or.erg.ge.xss(l+5+nint(xss(l)))) goto 82

         call urres(ic)

         goto 120

*-----------------------------------------------------------------------
*     get absorption and total cross sections.
*-----------------------------------------------------------------------

   82    rtc(3,iex) = xss(ic+2*n)+rtc(1,iex)
     &                   * (xss(ic+2*n+1)-xss(ic+2*n))
         rtc(4,iex) = xss(ic+n)+rtc(1,iex)
     &                   * (xss(ib+n)-xss(ic+n))
         rtcel(iex) = xss(ic+3*n)+rtc(1,iex)
     &                   * (xss(ic+3*n+1)-xss(ic+3*n))
         goto 120

   90    lc = nxs(16,iex)
         id = ktc(1,iex)
         goto 110

  100    ktc(1,iex) = id
         ic = id+jxs(1,iex)-1

  110    rtc(3,iex) = xss(ic+2*n)
         rtc(4,iex) = xss(ic+n)
         rtcel(iex) = xss(ic+3*n) ! frtati 2021/12/17
  120    rtc(6,iex) = erg
         rtc(7,iex) = -1.
  130    rtc(8,iex) = 0.

*-----------------------------------------------------------------------
*     adjust the cross sections for thermal effects if necessary.
*-----------------------------------------------------------------------

         if(rtc(11,iex).ge.0.) goto 260
      if(lmt(m).ne.0) then
         if(erg.lt.esa(lmt(m))) goto 160
      end if

         if(rtc(7,iex).eq.ttn) goto 260
         rtc(7,iex) = ttn
         rtc(5,iex) = rtc(4,iex)
         if(tbt(iex).ne.0.) goto 260
         a2 = awn(iex)*erg
         if(a2.gt.500.*ttn) goto 260
         if(a2.lt.4.*ttn) goto 140

cKN 2016/08/05 f -> ff
         ff = .5*ttn/a2
         goto 150


  140    a = sqrt(a2/ttn)
         b = 25.*a
         i = b

cKN 2016/08/05 f -> ff
         ff = (thgf(i)+(b-i)*(thgf(i+1)-thgf(i)))/a-1.

  150    i = ktc(1,iex)+jxs(1,iex)-1+3*nxs(3,iex)
         rtc(5,iex) = rtc(5,iex)
     &              + ff*(xss(i)+rtc(1,iex)*(xss(i+1)-xss(i)))
cKN 2016/08/05 f -> ff
         goto 260

*-----------------------------------------------------------------------
*     calculate the s(a,b) inelastic scattering cross section.
*-----------------------------------------------------------------------

  160    it = lmt(m)
         if(rtc(6,it).eq.erg) goto 250
         n = nint(xss(jxs(1,it)))
         ic = jxs(1,it)+1
         ib = jxs(1,it)+n
         if(erg.lt.rtc(6,it)) ib = ic+ktc(1,it)
         if(erg.gt.rtc(6,it)) ic = ic+ktc(1,it)-1
         rtc(6,it) = erg

  170    if(ib-ic.eq.1) goto 190
         ih = (ib+ic)/2
         if(erg.lt.xss(ih)) goto 180
         ic = ih
         goto 170

  180    ib = ih
         goto 170

  190    ktc(1,it) = ic-jxs(1,it)
         rtc(1,it) = (erg-xss(ic))/(xss(ib)-xss(ic))
         rtc(8,it) = xss(ic+n)+rtc(1,it)
     &                  * (xss(ic+n+1)-xss(ic+n))
         rtc(7,it) = rtc(8,it)

*-----------------------------------------------------------------------
*     calculate the s(a,b) SANS cross section.
*-----------------------------------------------------------------------

         if (isans.ne.0) then
            ITSANS = jxs(10,it)
            if (ITSANS /= 0) then
              SANS_MODEL = xss(ITSANS)
              if (SANS_MODEL .eq. 1) then
                SANS_sigma0 = xss(ITSANS+1)
                rtc(10,it) = SANS_sigma0/erg
                rtc(7,it) = rtc(7,it) + rtc(10,it) ! Add SANS to total scattering
              else
              ! --> Error, other models not implemented
                write(ErrCha,*) 'Error in xstneu: SANS model not impl.'
                ErrID = 'L:255/R:xstneu/F:ggm04.f' !E04_006_001
                call ErrWrite(ErrID,ErrCha)
                call parastop( 415 )
              endif
            end if
         endif
*-----------------------------------------------------------------------
*     calculate the s(a,b) elastic scattering cross section.
*-----------------------------------------------------------------------

         if(jxs(4,it).eq.0) goto 250 ! no elastic
         if(nxs(5,it).eq.5) goto 241 ! mixed elastic
         n = nint(xss(jxs(4,it)))
         ic = jxs(4,it)+1
         ib = jxs(4,it)+n
         if(nxs(5,it).ne.4) goto 200
         if(erg.le.xss(ic)) goto 250
         if(erg.ge.xss(ib)) goto 220

  200    if(ib-ic.eq.1) goto 230
         ih = (ib+ic)/2
         if(erg.lt.xss(ih)) goto 210
         ic = ih
         goto 200

  210    ib = ih
         goto 200

  220    ic = ib
  230    ktc(2,it) = ic-jxs(4,it)
         if(nxs(5,it).ne.4) goto 240
         rtc(7,it) = rtc(7,it)+xss(ic+n)/erg
         goto 250

  240    rtc(4,it) = (erg-xss(ic))/(xss(ib)-xss(ic))
         rtc(7,it) = rtc(7,it)+xss(ic+n)+rtc(4,it)
     &                  * (xss(ib+n)-xss(ic+n))
         goto 250
  241    n = nint(xss(jxs(4,it)))
         ic = jxs(4,it)+1
         ib = jxs(4,it)+n
         if(erg.le.xss(ic)) goto 246 ! skip coherent
         if(erg.ge.xss(ib)) goto 244 ! above table
  242    if(ib-ic.eq.1) goto 245
         ih = (ib+ic)/2
         if(erg.lt.xss(ih)) goto 243
         ic = ih
         goto 242

  243    ib = ih
         goto 242

  244    ic = ib
  245    ktc(2,it) = ic-jxs(4,it)
         rtc(7,it) = rtc(7,it)+xss(ic+n)/erg

  246    n = nint(xss(jxs(7,it)))
         ic = jxs(7,it)+1
         ib = jxs(7,it)+n

  247    if(ib-ic.eq.1) goto 249
         ih = (ib+ic)/2
         if(erg.lt.xss(ih)) goto 248
         ic = ih
         goto 247

  248    ib = ih
         goto 247

         ic = ib
  249    ktc(3,it) = ic-jxs(7,it)
         rtc(5,it) = (erg-xss(ic))/(xss(ib)-xss(ic))
         rtc(9,it) = xss(ic+n)+rtc(5,it)*(xss(ib+n)-xss(ic+n))
         rtc(7,it) = rtc(7,it)+rtc(9,it)
  250    rtc(5,iex) = rtc(7,it)+rtc(3,iex)

*-----------------------------------------------------------------------
*     accumulate the total cross section, totm, and the absorption
*     cross section, siga.
*-----------------------------------------------------------------------

  260    continue

*-----------------------------------------------------------------------

         totm = totm+rtc(5,iex)*fme(m)
         siga = siga+rtc(3,iex)*fme(m)

         siggc(mii) = rtc(5,iex)*fme(m)
         siggcn(mii)= (rtc(5,iex)-rtcel(iex))*fme(m) ! frtati 2021/12/17

         isigza(mii) = iza(m)
         if( iemode.gt.0 ) then
           itz = iza(m) / 1000
           ita = iza(m) - 1000 * itz
           if ( ita.eq.0 .and. itz.eq.6 ) ita = 12 ! S.H. 2022.3.8
           call sigrc(2,ein,ita,itz,sigt,signe,sigel)
           siggce(mii) = sigel*fme(m)
           if( iza(m).eq.1001 ) then ! frtati 2022/03/11
             siggcn(mii) = 0.d0
             siggce(mii) = 0.d0
           end if
         end if

      end if

*-----------------------------------------------------------------------

  270    continue

         isigc(0) = jmd(1+mk+1) - jmd(1+mk)


         if(mm.ne.0) goto 998

*-----------------------------------------------------------------------
*     set up for photon production biasing.
*-----------------------------------------------------------------------

         totgp1 = 0.
         if(npikmt.eq.0) goto 290
         if(gwt(icl).eq.-1e6) goto 290

         mii = 0

      do 280 m = jmd(1+mk), jmd(1+mk+1) - 1

*-----------------------------------------------------------------------

         mii = mii + 1
         siggc(mii) = 0.0d0

         lem = m-jmd(1+mk)+1
         dsmax = das_kmate(kmate(mk)+(lem-1)*5+11)

      if( dsmax .gt. emin(2) ) then

*-----------------------------------------------------------------------

         iex = lme(1,m)
         tgp(iex) = 0.
         if(nxs(15,iex).lt.0) goto 280
         j = jxs(12,iex)+ktc(1,iex)
         qu = xss(j-1)
         if(nty(iex).eq.1)
     &   qu = qu+rtc(1,iex)*(xss(j)-xss(j-1))
         tgp(iex) = qu
         totgp1 = totgp1+tgp(iex)*fme(m)

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

  280    continue

*-----------------------------------------------------------------------
*     calculate the fission cross section and nubar.
*-----------------------------------------------------------------------

  290    if(lfcl(icl).eq.0) goto 999 !return

         mii = 0

      do 300 m = jmd(1+mk), jmd(1+mk+1) - 1

*-----------------------------------------------------------------------

         mii = mii + 1
         siggc(mii) = 0.0d0

         lem = m-jmd(1+mk)+1
         dsmax = das_kmate(kmate(mk)+(lem-1)*5+11)

      if( dsmax .gt. emin(2) ) then

*-----------------------------------------------------------------------

         iex = lme(1,m)
         l = jxs(21,iex)
         if(l.eq.0) goto 300
         if(ktc(1,iex).lt.nint(xss(l))) goto 300
         j = l+2+ktc(1,iex)-nint(xss(l))
         rtc(8,iex) = xss(j)+rtc(1,iex)*(xss(j+1)-xss(j))
         if(rtc(11,iex).ge.0.) rtc(8,iex) = rtc(12,iex)
         fn = acenu(jxs(2,iex))
         rtc(10,iex) = max(fn,zero)

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

  300    continue
         goto 999

*-----------------------------------------------------------------------
*     compute multigroup cross sections.
*-----------------------------------------------------------------------

  310    totm = 0.
         pfp = 0.
         stp = 0.

         mii = 0

      do 330 m = jmd(1+mk), jmd(1+mk+1) - 1

*-----------------------------------------------------------------------

         mii = mii + 1
         siggc(mii) = 0.0d0

         lem = m-jmd(1+mk)+1
         dsmax = das_kmate(kmate(mk)+(lem-1)*5+11)

      if( dsmax .gt. emin(2) ) then

*-----------------------------------------------------------------------

         iex = lme(1,m)
         rtc(3,iex) = xss(jxs(6,iex)-1+jgp)
         rtc(5,iex) = xss(jxs(2,iex)-1+jgp)
         rtc(8,iex) = 0.
         if(nxs(11,iex).eq.0) goto 320
         if(jxs(7,iex).ne.0) stp = stp+fme(m)
     &                                * xss(jxs(7,iex)-1+jgp)
         if(jxs(8,iex).ne.0) pfp = pfp+fme(m)
     &                                * xss(jxs(8,iex)-1+jgp)

  320    if(lfcl(icl).eq.0) goto 330
         if(jxs(3,iex).eq.0) goto 330
         rtc(8,iex) = xss(jxs(3,iex)-1+jgp)
         rtc(10,iex) = xss(jxs(4,iex)-1+jgp)

         totm = totm+rtc(5,iex)*fme(m)

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

  330    continue

*-----------------------------------------------------------------------

  999    continue

*-----------------------------------------------------------------------

         sigt  = totm
         sigaa = siga

*-----------------------------------------------------------------------

  998    continue
         icl = icl0

*-----------------------------------------------------------------------
*     elastic option if exist
*-----------------------------------------------------------------------

            ioel = 0

         if( mlrgn .gt. 0 ) then

            do jj = 1, mlrgn

              do kk = melrg(jj,1), melrg(jj,2)

                 if( icl0 .eq. idgr(kk) ) then

                    ioel = jj
                    goto 500

                 end if

              end do

            end do

  500       continue

         end if

         if( ioel .eq. 0 ) return

*-----------------------------------------------------------------------

            if( ielusr .eq. 1 ) then

               call usrelst1(1,icels,ioel,ein,dum,sigbb)

            else

               call usrelst2(1,icels,ioel,ein,dum,sigbb)

            end if

               if( icels .eq. 1 ) return

               sigt = sigbb


               mii = 1
               siggc(mii) = sigt
               isigc(0)=mii ! T.Sato 2022/03/11
*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine urres(ic)
*                                                                      *
*       sample cross sections from probability tables for unresolved   *
*       resonance range and put in rtc( ,iex) array. ic = energy index *
*       l=jxs(23,iex) points to location for unres prob tables;   *
*           xss(l)  = number of incident energies                      *
*           xss(l+1)= length of a table (usually 20)                   *
*           xss(l+2)= incident energy interpolation parameter (2 or 5) *
*           xss(l+3)= flag for inelastic (-1 none, 0 by balance,>0=mt) *
*           xss(l+4)= flag for "other absorption" (-1 none, etc.)      *
*           xss(l+5)= flag for factors, if 0, regular cross sections   *
*           followed by incident energies                              *
*           followed by tables for each incident energy (6 columns)    *
*              cumulative probabilities                                *
*              total cross sections                                    *
*              elastic cross sections                                  *
*              fission cross sections                                  *
*              capture cross sections (102)                            *
*              neutron heating numbers                                 *
*       repeat for each incident energy.                               *
*                                                                      *
*       rtc( 1,iex) = interpolation fraction                           *
*       rtc( 2,iex) = interpolation fraction after free gas thermal    *
*       rtc( 3,iex) = absorption cross section                         *
*                   = rtc(14,iex) + any "other absorption"             *
*       rtc( 4,iex) = total cross section                              *
*       rtc( 5,iex) = rtc( 4,iex) = total xsec after free gas thermal  *
*       rtc( 6,iex) = neutron energy                                   *
*       rtc( 7,iex) = thermal temperature (unused for unresolved)      *
*       rtc( 8,iex) = zero or fission                                  *
*       rtc(11,iex) = elastic xsec (-1 if not in unresolved range)     *
*       rtc(12,iex) = fission xsec                                     *
*       rtc(13,iex) = neutron heating number                           *
*       rtc(14,iex) = (n,gamma) capture xsec                           *
*       rtc(15,iex) = random number used to sample xsecs               *
*                                                                      *
*       Last modified by K.Niita on 2009/10/06                         *
*                                                                      *
************************************************************************
      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      include 'err.inc'

*-----------------------------------------------------------------------

         n = nxs(3,iex)
         l = jxs(23,iex)
         n1 = nint(xss(l))
         ke = l+6
         n2 = nint(xss(l+1))
         je = ke+1
         jf = 0

*-----------------------------------------------------------------------
*     find energy table.
*-----------------------------------------------------------------------

   10    if(erg.lt.xss(je)) goto 20
         je = je+1
         goto 10

*-----------------------------------------------------------------------
*     get fractional distance between incident energy prob. tables.
*-----------------------------------------------------------------------

   20    r = (erg-xss(je-1))/(xss(je)-xss(je-1))
         jn = je-ke
         ja = ke+n1+(jn-1)*n2*6
         j = ja

*-----------------------------------------------------------------------
*     if correlated by temperature, find random number already used.
*-----------------------------------------------------------------------

         if(iunr.ge.0) goto 40

      do 30 ie = 1, mxe
         if(ie.eq.iex) goto 30
         if(nty(ie).ne.1) goto 30
         if(nxs(2,ie).ne.nxs(2,iex)) goto 30
         if(rtc(6,ie).ne.erg) goto 30

         lk = jxs(23,ie)
         if(lk.eq.0) goto 30
         if(erg.le.xss(lk+6)) goto 30
         if(erg.ge.xss(lk+5+nint(xss(lk)))) goto 30
         if(rtc(11,ie).lt.0.) goto 30

         rn = rtc(15,ie)
         goto 50

   30    continue
   40    rn = rang()

*-----------------------------------------------------------------------
*     sample from the probability tables.
*-----------------------------------------------------------------------

   50    rtc(15,iex) = rn
   60    if(xss(j).ge.rn) goto 70
         j = j+1
         goto 60

   70    if(j-ja+1.gt.n2) then
            write(ErrCha,*) 'Error in urres: prob table not normed.'
            ErrID = 'L:686/R:urres/F:ggm04.f' !E04_006_001
            call ErrWrite(ErrID,ErrCha)

            call parastop( 411 )
         end if

         if(jf.ne.0) goto 80

*-----------------------------------------------------------------------
*     set locator in probability table for lower incident energy.
*-----------------------------------------------------------------------

         j1 = j
         ja = ja+n2*6
         j = ja
         jf = 1
         goto 60

*-----------------------------------------------------------------------
*     set locator in probability table for upper incident energy.
*-----------------------------------------------------------------------

   80    j2 = j
         ki = nint(xss(l+2))

*-----------------------------------------------------------------------
*     obtain capture, elastic, fission, and heating from prob tables.
*     use lin-lin interpolation.
*-----------------------------------------------------------------------

         if(ki.eq.5.and.r.gt.0..and.r.lt.1.) goto 90
         rtc(3,iex) = xss(j1+n2*4)+r*(xss(j2+n2*4)-xss(j1+n2*4))
         rtc(11,iex) = xss(j1+n2*2)+r*(xss(j2+n2*2)-xss(j1+n2*2))
         rtc(12,iex) = xss(j1+n2*3)+r*(xss(j2+n2*3)-xss(j1+n2*3))
         rtc(13,iex) = xss(j1+n2*5)+r*(xss(j2+n2*5)-xss(j1+n2*5))
         goto 100

*-----------------------------------------------------------------------
*     use log-log interpolation.
*-----------------------------------------------------------------------

   90    g = log(erg/xss(je-1))/log(xss(je)/xss(je-1))
         i = n2*4
         rtc(3,iex) = 0.
         if(xss(j1+i).ne.0.)
     &      rtc(3,iex) = exp(log(xss(j1+i))
     &                      + g*log(xss(j2+i)/xss(j1+i)))
         i = n2*2
         rtc(11,iex) = exp(log(xss(j1+i))
     &                    + g*log(xss(j2+i)/xss(j1+i)))
         i = n2*3
         rtc(12,iex) = 0.
         if(xss(j1+i).ne.0.)
     &      rtc(12,iex) = exp(log(xss(j1+i))
     &                       + g*log(xss(j2+i)/xss(j1+i)))
         i = n2*5
         rtc(13,iex) = exp(log(xss(j1+i))
     &                    + g*log(xss(j2+i)/xss(j1+i)))

*-----------------------------------------------------------------------
*     get cross sections by factors (xss(l+5).ne.0) or get
*     inelastic cross section by balance (xss(l+3).eq.0).
*     ce and cf are the elastic and fission factors.
*-----------------------------------------------------------------------

  100    if(xss(l+5).eq.0..and.xss(l+3).ne.0.) goto 110
         ce = xss(ic+3*n)+rtc(1,iex)*(xss(ic+3*n+1)-xss(ic+3*n))
         cf = 0.
         i = jxs(21,iex)
         if(i.eq.0) goto 110
         ii = ktc(1,iex)-nint(xss(i))+1
         if(ii.lt.1) goto 110
         ii = ii+i+1
         cf = xss(ii)+rtc(1,iex)*(xss(ii+1)-xss(ii))

*-----------------------------------------------------------------------
*     get cc, other absorption factor.
*-----------------------------------------------------------------------

  110    if(xss(l+5).eq.0..and.xss(l+4).ne.0.) goto 150
         cc = 0.
         i = jxs(3,iex)-1

      do 120 ix = 1, nxs(4,iex)
  120    if(nint(xss(i+ix)).eq.102) goto 130
         goto 140

  130    iq = jxs(7,iex)+nint(xss(jxs(6,iex)+ix-1))
         jq = ktc(1,iex)+1-nint(xss(iq-1))
         if(jq.lt.1.or.jq.gt.nint(xss(iq))) goto 140
         if(jq.eq.nint(xss(iq)).and.rtc(1,iex).ne.0.) goto 140
         cc = xss(iq+jq)
         if(rtc(1,iex).ne.0.)
     &      cc = cc+rtc(1,iex)*(xss(iq+jq+1)-xss(iq+jq))

*-----------------------------------------------------------------------
*     probability tables are factors (xss(l+5).ne.0), so multiply
*     smooth xsecs by factors sampled from probability tables.
*     ch is heating number factor.
*-----------------------------------------------------------------------

  140    if(xss(l+5).eq.0.) goto 150
         ch = xss(ic+4*n)+rtc(1,iex)*(xss(ic+4*n+1)-xss(ic+4*n))
         rtc(11,iex) = rtc(11,iex)*ce
         rtc(12,iex) = rtc(12,iex)*cf
         rtc(3,iex) = rtc(3,iex)*cc
         rtc(13,iex) = rtc(13,iex)*ch

*-----------------------------------------------------------------------
*     obtain prob table total as sum of partials to avoid roundoff.
*     rtc(14,iex) is 102 capture and rtc(3,iex) is absorption.
*-----------------------------------------------------------------------

  150    rtc(4,iex) = rtc(11,iex)
     &                   + rtc(12,iex)+rtc(3,iex)
         if(jxs(2,iex).eq.0.and.rtc(12,iex).ne.0.) then
            write(ErrCha,*) 'Error in urres:'//
     &                 ' fission conflict with unresolved prob table.'
         ErrID = 'L:804/R:urres/F:ggm04.f' !E04_007_001
         call ErrWrite(ErrID,ErrCha)

            call parastop( 412 )
         end if

         rtc(14,iex) = rtc(3,iex)

*-----------------------------------------------------------------------
*     obtain inelastic, r1, and other absorption, r2, by balance
*     using smooth cross sections or from the mt numbers.
*-----------------------------------------------------------------------

         r1 = 0.
         r2 = 0.

*-----------------------------------------------------------------------
*     get ch, absorption factor, if (xss(l+3).eq.0) inelastic by
*     balance or (xss(l+4).eq.0) other absorption flag set.
*-----------------------------------------------------------------------
         if(xss(l+3).eq.0..or.xss(l+4).eq.0.) ca = xss(ic+2*n)
     &              + rtc(1,iex)*(xss(ic+2*n+1)-xss(ic+2*n))
         if(xss(l+3).gt.0.) goto 160
         if(xss(l+3).lt.0.) goto 190

*-----------------------------------------------------------------------
*     xss(l+3)=0:  get total inelastic by balance.  get smooth total,
*     ct, then subtract smooth elastic, fission, and absorption.
*-----------------------------------------------------------------------

         ct = xss(ic+n)+rtc(1,iex)*(xss(ic+n+1)-xss(ic+n))
         r1 = ct-ce-cf-ca
         goto 190

*-----------------------------------------------------------------------
*     xss(l+3)>0.  get inelastic from mt.
*     jxs(23,iex)+3 contains the mt of inelastic.
*-----------------------------------------------------------------------

  160    mt = nint(xss(l+3))
         i = jxs(3,iex)-1

      do 170 ix = 1, nxs(4,iex)
  170    if(nint(xss(i+ix)).eq.mt) goto 180

            write(ErrCha,*) 'Error in urres:'//
     &                 ' total inelastic mt not in table.'
            ErrID = 'L:851/R:urres/F:ggm04.f' !E04_008_001
            call ErrWrite(ErrID,ErrCha)

            call parastop( 413 )

  180    is = jxs(7,iex)+nint(xss(jxs(6,iex)+ix-1))
         ic = ktc(1,iex)+1-nint(xss(is-1))
         if(ic.lt.1.or.ic.gt.nint(xss(is))) goto 190
         if(ic.eq.nint(xss(is)).and.rtc(1,iex).ne.0.) goto 190
         r1 = xss(is+ic)
         if(rtc(1,iex).ne.0.) r1 = r1
     &              + rtc(1,iex)*(xss(is+ic+1)-xss(is+ic))

*-----------------------------------------------------------------------
*     check flag for other absorption to include in abs and total.
*     obtain other absorption by smooth (abs-102).
*-----------------------------------------------------------------------

  190    if(xss(l+4).lt.0.) goto 230
         if(xss(l+4).gt.0.) goto 200
         r2 = ca-cc
         goto 230

*-----------------------------------------------------------------------
*     jxs(23,iex)+4 contains the mt of other absorption.
*-----------------------------------------------------------------------

  200    mt = nint(xss(l+4))
         i = jxs(3,iex)-1

      do 210 ix = 1,nxs(4,iex)
  210    if(nint(xss(i+ix)).eq.mt) goto 220

            write(ErrCha,*) 'Error in urres:'//
     &                 ' other abs mt not in table.'
            ErrID = 'L:886/R:urres/F:ggm04.f' !E04_009_001
            call ErrWrite(ErrID,ErrCha)

            call parastop( 414 )

  220    is = jxs(7,iex)+nint(xss(jxs(6,iex)+ix-1))
         ic = ktc(1,iex)+1-nint(xss(is-1))
         if(ic.lt.1.or.ic.gt.nint(xss(is))) goto 230
         if(ic.eq.nint(xss(is)).and.rtc(1,iex).ne.0.) goto 230
         r2 = xss(is+ic)
         if(rtc(1,iex).ne.0.) r2 = r2
     &                  + rtc(1,iex)*(xss(is+ic+1)-xss(is+ic))

*-----------------------------------------------------------------------
*     increment cross sections for inelastic and other absorption.
*-----------------------------------------------------------------------

  230    rtc(3,iex) = rtc(3,iex)+r2
         rtc(4,iex) = rtc(4,iex)+r1+r2
         rtc(5,iex) = rtc(4,iex)

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      function acenu(l)
*                                                                      *
*        calculate the average number of neutrons from fission.        *
*        the data for all real nuclides are such that acenu will always*
*        be greater than 1 and almost always greater than 2.           *
*        data (for polynomial function):  lnu,nc,(c(i),i=1,nc)         *
*        l is the location of the nubar table.                         *
*                                                                      *
************************************************************************
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------
c
c        find out which kind of calculation to do.
      if(nint(xss(l)).ne.1)go to 20
c
c        calculate acenu by evaluating a polynomial in energy.
      n=nint(xss(l+1))-1
      acenu=xss(l+n+2)
      do 10 i=1,n
   10 acenu=acenu*erg+xss(l+n+2-i)
      return
c
c        calculate acenu by interpolation in a table.
   20 acenu=acefcn(l+1,erg,i)
      return
      end


************************************************************************
*                                                                      *
      function acefcn(l,eg,ln)
*                                                                      *
*       evaluate a function of eg from the table at xss(l).            *
*       also return ln, the length of the table.                       *
*       data:  nr,nbttc(i=1,nr),int(i=1,nr),nf,e(i=1,nf),f(i=1,nf)     *
*       Last modified by K.Niita on 2009/09/30                         *
*                                                                      *
************************************************************************
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------

         nr = nint(xss(l))
         ie = l+2*nr+1
         nf = nint(xss(ie))
         ln = 2*(nr+nf+1)

*-----------------------------------------------------------------------
*     use the extreme value if eg is off either end of the table.
*-----------------------------------------------------------------------

         if(eg.ge.xss(ie+nf)) goto 110
         if(eg.le.xss(ie+1)) goto 120

         if(dabs(eg - 20.d0) .le.1.d-5) eg = 20.d0 - 1.d-10 ! 2024/9/24 Ogawa. If primary is 20MeV neutron, ea is 20, eb is 200, fa and fb are zero causing X-section error.
*-----------------------------------------------------------------------
*     binary search for the location of eg in the table.
*-----------------------------------------------------------------------

         ic = ie+1
         ib = ie+nf
   10    if(ib-ic.eq.1) goto 30
         ih = (ic+ib)/2
         if(eg.lt.xss(ih)) goto 20
         ic = ih
         goto 10

   20    ib = ih
         goto 10

   30    ea = xss(ic)
         eb = xss(ib)
         fa = xss(ic+nf)
         fb = xss(ib+nf)

*-----------------------------------------------------------------------
*     find out which kind of interpolation should be used.
*-----------------------------------------------------------------------

         if(nr.eq.0) goto 70

      do 40 n = 1, nr
   40    if(ib-ie.le.nint(xss(l+n))) goto 50
         n = nr

*-----------------------------------------------------------------------
*     interpolate between table entries.
*-----------------------------------------------------------------------

   50    goto (60,70,80,90,100) nint(xss(l+nr+n))

   60    acefcn = fa
         return
   70    acefcn = fa+(fb-fa)*(eg-ea)/(eb-ea)
         return
   80    acefcn = fa+(fb-fa)*log(eg/ea)/log(eb/ea)
         return
   90    acefcn = fa*(fb/fa)**((eg-ea)/(eb-ea))
         return
  100    acefcn = fa*(fb/fa)**(log(eg/ea)/log(eb/ea))
         return

  110    acefcn = xss(ie+2*nf)
         return
  120    acefcn = xss(ie+nf+1)

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine xstgam(mm,sigt,sigaa,ein,mk)
*                                                                      *
*       calculate photon cross sections in cell icl.                   *
*       Last modified by K.Niita on 2009/12/22                         *
*                                                                      *
*       rtc( 1,iex) = incoherent (Compton) scattering                  *
*       rtc( 2,iex) = rtc( 1,iex) + coherent (Thomson) scattering      *
*       rtc( 3,iex) = rtc( 2,iex) + photo-electric effect              *
*       rtc( 4,iex) = rtc( 3,iex) + pair production                    *
*                                                                      *
************************************************************************
      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'


      common /kmat1g/ kmat(kvlmax)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)

      common /adjoint/ iadjnt

*-----------------------------------------------------------------------
*     initialization of this material and particle
*-----------------------------------------------------------------------

         erg  = ein
         imat = mk

*-----------------------------------------------------------------------

         totm = 0.
         siga = 0.
         if(mcal.ne.0) goto 80
         el = log(erg)

*-----------------------------------------------------------------------
*     calculate photonuclear cross sections.
*-----------------------------------------------------------------------

!           add totpn to totm for distance to collision calculations.

*-----------------------------------------------------------------------
*     calculate cross sections for each nuclide in the cell.
*-----------------------------------------------------------------------

      do 70 m = jmd(1+mk), jmd(1+mk+1) - 1

         ie = lme(2,m)
         if(erg.eq.rtc(6,ie)) goto 60
         nz = nxs(3,ie)

*-----------------------------------------------------------------------
*     binary search for the location of log(energy) in the table.
*-----------------------------------------------------------------------

         ic = jxs(1,ie)
         ib = ic+nz-1
   10    if(ib-ic.eq.1) goto 30
         ih = (ic+ib)/2
         if(el.lt.xss(ih)) goto 20
         ic = ih
         goto 10
   20    ib = ih
         goto 10

*-----------------------------------------------------------------------
*     generate incoherent, coherent, photo-electric, and pair
*     production cross sections and store in rtc.
*-----------------------------------------------------------------------

   30    t1 = (el-xss(ib))/(xss(ib-1)-xss(ib))
         k = ib-jxs(1,ie)+jxs(5,ie)
         rtc(5,ie) = xss(k)*(xss(k-1)/xss(k))**t1
         c = 0.
         t2 = 0.

      do 50 j = 1, 4
         ib = ib+nz
         t3 = 0.
         if(xss(ib)*xss(ib-1).eq.0.) goto 50
         t3 = exp(xss(ib)+t1*(xss(ib-1)-xss(ib)))

*-----------------------------------------------------------------------
*     if mcg treatment, set coherent cross section to zero and
*     adjust heating with ratio of total cross sections.
*-----------------------------------------------------------------------

         if(j.ne.2) goto 40
         if(erg.lt.emcf(2).and.nocoh.eq.0) goto 40
         c = t3
         t3 = 0.
   40    t2 = t2+t3
   50    rtc(j,ie) = t2
         if(erg.ge.emcf(2)) rtc(5,ie) = rtc(5,ie)
     &                       * (rtc(4,ie)+c)/rtc(4,ie)
         rtc(6,ie) = erg

*-----------------------------------------------------------------------
*     accumulate the total cross section, totm.
*-----------------------------------------------------------------------

   60    continue

      if( iadjnt .ne. 1 ) then ! T.Sato 2024/01/07, change .eq.0 to .ne.1 because iadjnt = 2 is introduced

         totm = totm + rtc(4,ie) * fme(m)
         siga = siga + (rtc(3,ie)-rtc(2,ie)) * fme(m)

*-----------------------------------------------------------------------
*     Adjoint mode
*-----------------------------------------------------------------------

      else
cAlex 2017/02/16 add missing lines of code
         ie = lme(2,m)
         nz = nxs(3,ie)

*-----------------------------------------------------------------------
*     binary search for the location of log(energy) in the table.
*-----------------------------------------------------------------------

         icc = jxs(6,ie)
cAlex 2017/02/15 extra cross section in adjoint library for XS_adjincoh peak
         ibb = icc+nz

  110    if(ibb-icc.eq.1) goto 130
         ihh = (icc+ibb)/2
         if(el.lt.xss(ihh)) goto 120
         icc = ihh
         goto 110
  120    ibb = ihh
         goto 110
  130    continue

*-----------------------------------------------------------------------

         t1 = (el-xss(ibb))/(xss(ibb-1)-xss(ibb))
cAlex 2017/02/15 Remove unnecessary code, add XS_adjincoh only
         ibb = ibb+nz+1
         t3 = 0.
         if(xss(ibb)*xss(ibb-1).eq.0.) goto 155
         t3 = exp(xss(ibb)+t1*(xss(ibb-1)-xss(ibb)))
  155    rtc(8,ie) = t3
         rtc(9,ie) = 0.0
  150    continue
  160    continue

cAlex 2017/02/15 totm is total forward cross section, siga is adjoint scat XS
         totm = totm+rtc(4,ie)*fme(m)
         siga=siga+(rtc(8,ie)+rtc(9,ie)+rtc(2,ie)-rtc(1,ie))*fme(m)

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

   70    continue

         goto 999

*-----------------------------------------------------------------------
c        compute multigroup cross sections.
*-----------------------------------------------------------------------

   80 do 90 m = jmd(1+mk), jmd(1+mk+1) - 1

         ie = lme(2,m)
         rtc(4,ie) = xss(jxs(2,ie)-1+jgp)
         totm = totm+rtc(4,ie)*fme(m)
   90    siga = siga+xss(jxs(6,ie)-1+jgp)*fme(m)

*-----------------------------------------------------------------------

  999    continue

         sigt  = totm
         sigaa = siga

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine xstpni(mm,sigt,ein,mk)
*                                                                      *
*       calculate photonuclear cross sections in cell icl.             *
*       Last modified by K.Niita on 2009/12/22                         *
*                                                                      *
************************************************************************
      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      use MEMBANKMOD ! frtati 2021/12/17
      use moddas_material
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      common /pnint/ ipnint


      common /kmat1g/ kmat(kvlmax)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)

      common /pnixs/ totmpni,pnixsm(pnlmax),ipnikind(pnlmax),egypni,
     &               dlibmax
!$OMP THREADPRIVATE(/pnixs/)

      common /dpnmaxcom/ dpnmax ! photon nuclear library maximum energy

*-----------------------------------------------------------------------
*     initialization of this material and particle
*-----------------------------------------------------------------------

      erg  = ein
      imat = mk
      totmpni = 0.0d0

*-----------------------------------------------------------------------
*     photonuclear part
*-----------------------------------------------------------------------

      if( ipnint .ne. 0 ) then


        isigpn0 = jmd(1+mk+1) - jmd(1+mk) ! frtati 2021/12/17

        do m = jmd(1+mk), jmd(1+mk+1) - 1

         ie = lme(2,m)
         izm = iza(m) / 1000
         ims = iza(m) - 1000 * izm
         icntm = m - jmd(1+mk) + 1
         siggpn(icntm) =0.d0

*-----------------------------------------------------------------------
*     expand natural nucleus and calculate average x-section
*-----------------------------------------------------------------------

         if( ims .eq. 0 ) then

          lemm  = nint( dnel_das(kmat0+imat) )

          sekt = 0.0d0
          sekc = 0.0d0

          do lem = 1, lemm

           izt = zz_das(kmat(imat)+lem)

           if( izt .eq. izm ) then

            sekt = sekt + den_das(kmat(imat)+lem)
            inm = nint( a_das(kmat(imat)+lem) ) - izm
            xsc = phxs(izm,inm,erg)

            sekc = sekc + xsc * den_das(kmat(imat)+lem)

           end if

          end do

          xsc = sekc / sekt

         else

          inm = ims - izm
          xsc = phxs(izm,inm,erg)

         end if

         rtc(7,ie) = xsc
ccABE bugfix @2014/11/11
         pnixsm(icntm) = xsc
         totmpni = totmpni + pnixsm(icntm) * fme(m)
         siggpn(icntm) = pnixsm(icntm) * fme(m)
         isigza(icntm) = iza(m) ! S.H. 2024.10.15

        end do

*-----------------------------------------------------------------------

      else ! ipnint = 0, no cross section for photo nuclear data

       do m = jmd(1+mk), jmd(1+mk+1) - 1
        rtc(7,ie) = 0.0
ccABE bugfix @2014/11/11
        icntm = m - jmd(1+mk) + 1
        pnixsm(icntm) = 0.0d0
       enddo

      end if

*-----------------------------------------------------------------------

      sigt  = totmpni

*-----------------------------------------------------------------------
      return
      end

************************************************************************
*                                                                      *
      subroutine heath(icli,ergi,heatr,mki,rhi,mtal,nm,mt)
*                                                                      *
*       calculate proton heat in cell icl.                             *
*       Last modified by T.Sato 2025/03/04 for adding mother parameter *
*                                                                      *
************************************************************************
      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      use moddas_material
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      common /eparm/  esmax, esmin, emin(20)
      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax)

! T.Sato 2025/03/04 for mother parameter
      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)
      dimension   mt(nm)

*-----------------------------------------------------------------------
*     initialization of this material and density
*-----------------------------------------------------------------------

         mk  = mki
         rh  = rhi
         heatr = 0.0d0
         ityp = 1

         if( mk .le. 0 ) return

         call sig_tot(ityp,sigt,sigaa,ergi,mk)

*-----------------------------------------------------------------------
*     proton heating from normal reaction
*-----------------------------------------------------------------------

      h = 0.
      ipt = 9

      do 270 m = jmd(1+mk), jmd(1+mk+1) - 1

*-----------------------------------------------------------------------

! T.Sato 2025/03/04 for mother parameter
         mchg = iza(m) / 1000
         mmas = iza(m) - 1000 * mchg
         if( nm .gt. 0 ) then
          do i = 1, nm
           iz = mt(i) / 1000
           ia = mt(i) - iz * 1000
           if( ( ia .gt. 0 .and.
     &      mmas .eq. ia .and. mchg .eq. iz ) .or.
     &      ( ia .eq. 0 .and. mchg .eq. iz ) ) then
            if( itmct(mtal) .gt. 0 ) goto 140 ! consider this nuclide
            if( itmct(mtal) .lt. 0 ) goto 270 ! skip this nuclide
           end if
          end do
          if( itmct(mtal) .gt. 0 ) goto 270 ! skip this nuclide
         end if


 140     lem = m-jmd(1+mk)+1
         dsmax = das_kmate(kmate(mk)+(lem-1)*5+11)

      if( dsmax .gt. emin(1) ) then

*-----------------------------------------------------------------------

         i = lme(ipt,m)
         l = jxs(1,i)+4*nxs(3,i)+ktc(1,i)

            tx = rtc(5,i)
     &         * ( xss(l-1) + rtc(1,i) * ( xss(l) - xss(l-1) ) )

            h2 = fme(m) * max(zero,tx)

            h = h + h2

      end if

  270 continue

            heatr = h * rh

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine heatn(icli,ergi,heatr,heatf,mki,rhi,mtal,nm,mt)
*                                                                      *
*       calculate neutron heat in cell icl.                            *
*       Last modified by T.Sato 2025/03/04 for adding mother parameter *
*                                                                      *
************************************************************************

      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      use moddas_material
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

      common /eparm/  esmax, esmin, emin(20)
      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax)

! T.Sato 2025/03/04 for mother parameter
      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)
      dimension   mt(nm)

*-----------------------------------------------------------------------
*     initialization of this material and density
*-----------------------------------------------------------------------

         mk  = mki
         rh  = rhi

         heatr = 0.0d0
         heatf = 0.0d0
         tme   = 0.

         if( mk .le. 0 ) return

         call xstneu(1,sigt,sigaa,icli,ergi,tme,mk)

*-----------------------------------------------------------------------
*     neutron heating from normal reaction
*-----------------------------------------------------------------------

      h = 0.
      ipt = 1

      do 270 m=jmd(1+mk),jmd(1+mk+1)-1

! T.Sato 2025/03/04 for mother parameter
         mchg = iza(m) / 1000
         mmas = iza(m) - 1000 * mchg
         if( nm .gt. 0 ) then
          do i = 1, nm
           iz = mt(i) / 1000
           ia = mt(i) - iz * 1000
           if( ( ia .gt. 0 .and.
     &      mmas .eq. ia .and. mchg .eq. iz ) .or.
     &      ( ia .eq. 0 .and. mchg .eq. iz ) ) then
            if( itmct(mtal) .gt. 0 ) goto 140 ! consider this nuclide
            if( itmct(mtal) .lt. 0 ) goto 270 ! skip this nuclide
           end if
          end do
          if( itmct(mtal) .gt. 0 ) goto 270 ! skip this nuclide
         end if

 140     lem = m-jmd(1+mk)+1
         dsmax = das_kmate(kmate(mk)+(lem-1)*5+12)

      if( dsmax .gt. emin(2) ) then

*-----------------------------------------------------------------------

         i = lme(ipt,m)
         l = jxs(1,i)+4*nxs(3,i)+ktc(1,i)

         if( rtc(11,i) .lt. 0. ) then

            tx = rtc(5,i)
     &         * ( xss(l-1) + rtc(1,i) * ( xss(l) - xss(l-1) ) )

         else if( rtc(11,i) .ge. 0. ) then

            tx = rtc(4,i) * rtc(13,i)

         end if

            h2 = fme(m) * max(zero,tx)

            h = h + h2

      end if

  270 continue

      heatr = h * rh

      return ! T.Sato 2024/12/22 to exclude fission contribution

*-----------------------------------------------------------------------
*     neutron heating from fission
*-----------------------------------------------------------------------

      h=0.

      do 350 m = jmd(1+mk), jmd(1+mk+1) - 1

*-----------------------------------------------------------------------

         lem = m-jmd(1+mk)+1
         dsmax = das_kmate(kmate(mk)+(lem-1)*5+12)

      if( dsmax .gt. emin(2) ) then

*-----------------------------------------------------------------------

         h2 = 0.
         i = lme(1,m)
         if(jxs(21,i).eq.0) goto 350
         l = jxs(21,i)
         if(ktc(1,i).lt.xss(l)) goto 350
         j = l+2+ktc(1,i)-nint(xss(l))

      do 330 k = 1, 22
  330    if(nxs(2,i).eq.mfiss(k)) goto 340

  340    if(rtc(11,i).lt.0.)
     &      h2 = fme(m)*qfiss(k)
     &         * (xss(j)+rtc(1,i)*(xss(j+1)-xss(j)))

         if(rtc(11,i).ge.0.)
     &      h2 = fme(m)*qfiss(k)*rtc(12,i)

         h = h + h2

      end if

  350 continue

         heatf = h * rh

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine heatp(icli,ergi,heatr,mki,rhi,mtal,nm,mt)
*                                                                      *
*       calculate photon heat in cell icl.                             *
*       Last modified by T.Sato 2025/03/04 for adding mother parameter *
*                                                                      *
************************************************************************
      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

! T.Sato 2025/03/04 for mother parameter
      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)
      dimension   mt(nm)

*-----------------------------------------------------------------------
*     initialization of this material and density
*-----------------------------------------------------------------------

         mk  = mki
         rh  = rhi

         heatr = 0.0d0

         if( mk .le. 0 ) return

         call xstgam(1,sigt,sigaa,ergi,mk)

*-----------------------------------------------------------------------
*     photon heating from normal reaction
*-----------------------------------------------------------------------

         h = 0.
  280 do 290 m = jmd(1+mk), jmd(1+mk+1) - 1

! T.Sato 2025/03/04 for mother parameter
         mchg = iza(m) / 1000
         mmas = iza(m) - 1000 * mchg
         if( nm .gt. 0 ) then
          do i = 1, nm
           iz = mt(i) / 1000
           ia = mt(i) - iz * 1000
           if( ( ia .gt. 0 .and.
     &      mmas .eq. ia .and. mchg .eq. iz ) .or.
     &      ( ia .eq. 0 .and. mchg .eq. iz ) ) then
            if( itmct(mtal) .gt. 0 ) goto 140 ! consider this nuclide
            if( itmct(mtal) .lt. 0 ) goto 290 ! skip this nuclide
           end if
          end do
          if( itmct(mtal) .gt. 0 ) goto 290 ! skip this nuclide
         end if

  140    i = lme(2,m)
         h2 = rtc(5,i)*rtc(4,i)*fme(m)
         h = h+h2
  290    continue

      heatr = h * rh

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine getcrs(ityp,icli,ergi,lma,mtid,gxss,melem)
*                                                                      *
*        calculate cross section from library.                         *
*        last modified by K.Niita on 2009/09/30                        *
*                                                                      *
*        ityp  : particle id                                           *
*        ergi  : energy                                                *
*        lma   : material id                                           *
*        mtid  : mt id                                                 *
*        gxss  : cross section for ergi                                *
*        melem : element ID, 0 for calculating all element             *
*                                                                      *
************************************************************************
      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------
*     particle, material and mt id
*-----------------------------------------------------------------------

            gxss = 0.

         if( ityp .eq. 1 ) then

            ipim = 9

         else if( ityp .eq. 2 ) then

            ipim = 1

         else if( ityp .eq. 14 ) then

            ipim = 2

cfrtati 2023/12/07
         else if( ityp .eq. 15 ) then
            ipim = 31
         else if( ityp .eq. 18 ) then
            ipim = 34

         else

            return

         end if

*-----------------------------------------------------------------------

         icl0 = icl

         erg = ergi
         mk  = lma
         mt  = mtid
         icl = icli

*-----------------------------------------------------------------------
*           neutron for total or elastic cross section
*-----------------------------------------------------------------------

            if( ( ipim .eq. 1 ) .and.
     &          ( mt .eq. 1 .or. mt .eq. 2 ) ) then

                  call xstneu(1,sigt,sigaa,icli,ergi,tme,mk)

            end if

            if( ipim .eq. 2 ) then

                  call xstgam(1,sigt,sigaa,ergi,mk)

            end if

cfrtati 2023/12/07 enable deuteron and alpha libraries
!            if( ipim .eq. 9 .and.
            if( (ipim.eq.9 .or. ipim.eq.31 .or. ipim.eq.34) .and.
     &          mt .ne. 444 .and. mt .ne. 445 ) then

                  call sig_tot(ityp,sigt,sigaa,ergi,mk)

            end if

*-----------------------------------------------------------------------
*     cross section from library y-type
*-----------------------------------------------------------------------

            sek = 0.

      if(melem.eq.0) then  ! conventional case
       mstart=jmd(1+mk)
       mend=jmd(1+mk+1)-1
      else  ! only one element, introduced on 2024/12/19
       mstart=melem
       mend=melem
      endif

      do m = mstart,mend
!      do 130 m = jmd(1+mk), jmd(1+mk+1) - 1

            iex = lme(ipim,m)

            sekc = 0.d0

*-----------------------------------------------------------------------
*   ---- dosimetry table ----
*-----------------------------------------------------------------------

         if( nty(iex) .eq. 3 ) then

               i = jxs(3,iex) - 1

            do ix = 1, nxs(4,iex)

               if( xss(i+ix) .eq. mt ) then

                  sekc = acefcn(jxs(7,iex)
     &                 + nint(xss(jxs(6,iex)+ix-1))-1,erg,i)

                  goto 130

               end if

            end do

*-----------------------------------------------------------------------
*   ---- other cross sections ----
*-----------------------------------------------------------------------

         else

*-----------------------------------------------------------------------
*        ---- neutron ----
*-----------------------------------------------------------------------

            if( ipim .eq. 1 ) then

*-----------------------------------------------------------------------

               if( rtc(6,iex) .ne. erg ) then

                  n = nxs(3,iex)
                  rtc(11,iex) = -1.
                  if(nxs(16,iex).eq.lc) goto 70
                  ic = jxs(1,iex)
                  ib = ic-1+n
                  if(erg.lt.rtc(6,iex)) ib = ic+ktc(1,iex)
                  if(erg.gt.rtc(6,iex)) ic = ic+ktc(1,iex)-1
   30             if(ib-ic.eq.1) goto 50
                  ih = (ic+ib)/2
                  if(erg.lt.xss(ih)) goto 40
                  ic = ih
                  goto 30
   40             ib = ih

*-----------------------------------------------------------------------
*                 get the bin index and interpolation fraction.
*-----------------------------------------------------------------------

                  goto 30
   50             ktc(1,iex) = ic-jxs(1,iex)+1
                  if(nty(iex).eq.2) goto 60
                  rtc(1,iex) = max(zero,min(one,(erg-xss(ic))
     &                            / (xss(ib)-xss(ic))))

*-----------------------------------------------------------------------
*                 sample cross sections from unresolved range
*                 probability tables.
*-----------------------------------------------------------------------

                  if(iunr.eq.0) goto 52
                  l = jxs(23,iex)
                  if(l.eq.0) goto 52
                  if(erg.le.xss(l+6).or.erg.ge.xss(l+5+nint(xss(l))))
     &            goto 52
                  call urres(ic)
                  rtc(8,iex) = rtc(12,iex)
                  if(lfcl(icl).eq.0.or.jxs(21,iex).eq.0)
     &             rtc(8,iex) = 0.
                  if(rtc(8,iex).ne.0.)
     &            rtc(10,iex) = acenu(jxs(2,iex))
                  goto 90

*-----------------------------------------------------------------------
*                 get absorption and total cross sections.
*-----------------------------------------------------------------------

   52             rtc(3,iex) = xss(ic+2*n)+rtc(1,iex)*
     &             (xss(ic+2*n+1)-xss(ic+2*n))
                  rtc(4,iex) = xss(ic+n)+rtc(1,iex)
     &                           *(xss(ib+n)-xss(ic+n))
                  goto 90
   60             lc = nxs(16,iex)
                  id = ktc(1,iex)
                  goto 80
   70             ktc(1,iex) = id
                  ic = id+jxs(1,iex)-1
   80             rtc(3,iex) = xss(ic+2*n)
                  rtc(4,iex) = xss(ic+n)
   90             rtc(6,iex) = erg
                  rtc(7,iex) = -1

               end if

*-----------------------------------------------------------------------

                     sekc = getxs(mt)

*-----------------------------------------------------------------------
*        ---- photon ----
*-----------------------------------------------------------------------

            else if( ipim .eq. 2 ) then

               j  = -mt
               xs = rtc(j,iex)

               if( j .gt. 1 ) xs = xs - rtc(j-1,iex)
               if( j .gt. 4 ) xs = rtc(j-1,iex)

               sekc = xs

*-----------------------------------------------------------------------
*        ---- proton ----
*-----------------------------------------------------------------------
cfrtati 2023/12/07 enable deuteron and alpha libraries
            else if( ipim .eq. 9 .or.
     &               ipim .eq. 31 .or. ipim .eq. 34 ) then

*-----------------------------------------------------------------------

                  xs = 0.

               if( mt .eq. 1 ) then

                  xs = rtc(5,iex)

cfrtati 2023/12/07 enable elastic and non-elastic xs
               else if( mt .eq. 2 ) then
                  xs = rtcel(iex)
               else if( mt .eq. 3 ) then
                  xs = rtc(5,iex) - rtcel(iex)

               else if( mt .eq. -4 ) then

                  i = jxs(1,iex)+4*nxs(3,iex)+ktc(1,iex)

                  xs = ( xss(i-1) + rtc(1,iex)
     &               * ( xss(i) - xss(i-1) ) )

               end if

               sekc = xs

*-----------------------------------------------------------------------

            end if

         end if

*-----------------------------------------------------------------------

  130    if(melem.eq.0) then ! conventional mode
          sek = sek + fme(m) * sekc
         else  ! one element mode
          sek = sek + sekc
         endif

         enddo  ! do loop for element

         gxss = sek

         icl = icl0

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function getxs(mt)
*                                                                      *
*        return with cross section or other quantity identified by mt  *
*        for energy erg in cross-section table iex.                    *
*        last modified by K.Niita on 2009/09/30                        *
*                                                                      *
************************************************************************
      use GGMBANKMOD !FURUTA
      use GGMARRAYMOD !2020ASTOM
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------

         getxs = 0.

*-----------------------------------------------------------------------
*     regular table.
*-----------------------------------------------------------------------

   30    if(mt.gt.1000) goto 170
cfrtati 2023/12/07 mt=3 added
!         goto (60,70,160,90,90,90,90,100,110,120,130) 3-mt
         goto (60,60,70,160,90,90,90,90,100,110,120,130) 4-mt

*-----------------------------------------------------------------------
* >>>>>  mt>2 -- a standard endf reaction cross section.
*-----------------------------------------------------------------------
         i = jxs(3,iex)-1

      do 40 ix = 1,nxs(4,iex)
   40    if(nint(xss(i+ix)).eq.mt) goto 50
         return

   50    if(rtc(11,iex).lt.0.) goto 55
         if(mt.eq.18.or.mt.eq.19) goto 53
         if(mt.ne.102) goto 55
         getxs = rtc(14,iex)
         return

   53    continue
         getxs = rtc(12,iex)
         return

   55    is = jxs(7,iex)+nint(xss(jxs(6,iex)+ix-1))
         ic = ktc(1,iex)+1-nint(xss(is-1))
         if(ic.lt.1.or.ic.gt.nint(xss(is))) return
         if(ic.eq.nint(xss(is)).and.rtc(1,iex).ne.0.) return
         getxs = xss(is+ic)+rtc(1,iex)*(xss(is+ic+1)-xss(is+ic))
         return

*-----------------------------------------------------------------------
* >>>>>  mt=2 -- elastic scattering cross section.
*-----------------------------------------------------------------------

   60    if(rtc(11,iex).lt.0.) goto 65
         getxs = rtc(11,iex)
         if( mt.eq.3 ) getxs = rtc(5,iex) - getxs ! frtati 2023/12/07
         return

   65    i = jxs(1,iex)+3*nxs(3,iex)+ktc(1,iex)
         getxs = xss(i-1)+rtc(1,iex)*(xss(i)-xss(i-1))+
     &    rtc(5,iex)-rtc(4,iex)
         if( mt.eq.3 ) getxs = rtc(5,iex) - getxs ! frtati 2023/12/07
         return

*-----------------------------------------------------------------------
* >>>>>  mt=1 -- total cross section.
*-----------------------------------------------------------------------

   70    getxs = rtc(5,iex)
         return

*-----------------------------------------------------------------------
* >>>>>  mt=-1,-2,-3,-4 -- total, absorption, elastic, heating.
*-----------------------------------------------------------------------

   90    if(rtc(11,iex).lt.0.) goto 99
         if(mt.eq.-1) getxs = rtc(4,iex)
         if(mt.eq.-2) getxs = rtc(3,iex)
         if(mt.eq.-3) getxs = rtc(11,iex)
         if(mt.eq.-4) getxs = rtc(13,iex)
         return

   99    i = jxs(1,iex)-mt*nxs(3,iex)+ktc(1,iex)
         getxs = max(zero,xss(i-1)+rtc(1,iex)*(xss(i)-xss(i-1)))
         return

*-----------------------------------------------------------------------
* >>>>>  mt=-5 -- gamma production cross section.
*-----------------------------------------------------------------------
*        the infinite dilute photon production cross sections are used
*        in tallies within the unresolved resonance energy range.
*-----------------------------------------------------------------------

  100    if(jxs(12,iex).eq.0) return
         i = jxs(12,iex)+ktc(1,iex)
         getxs = xss(i-1)+rtc(1,iex)*(xss(i)-xss(i-1))
         return

*-----------------------------------------------------------------------
* >>>>>  mt=-6 -- total fission cross section.
*-----------------------------------------------------------------------

  110    if(jxs(21,iex).eq.0) return
         if(rtc(11,iex).lt.0.) goto 115
         getxs = rtc(12,iex)
         return

  115    is = jxs(21,iex)+1
         ic = ktc(1,iex)+1-nint(xss(is-1))
         if(ic.ge.1) getxs = xss(is+ic)+rtc(1,iex)*
     &    (xss(is+ic+1)-xss(is+ic))
         return

*-----------------------------------------------------------------------
* >>>>>  mt=-7 -- fission nu.
*-----------------------------------------------------------------------

  120    if(jxs(21,iex).ne.0) getxs = acenu(jxs(2,iex))
         return

*-----------------------------------------------------------------------
* >>>>>  mt=-8 -- fission q.
*-----------------------------------------------------------------------

  130 do 140 i = 1, 22
  140    if(nxs(2,iex).eq.mfiss(i)) goto 150

  150    getxs = qfiss(i)
  160    return

*-----------------------------------------------------------------------
* >>>>>  mt>1000 -- partial photon-production cross section.
*-----------------------------------------------------------------------

  170    if(jxs(13,iex).eq.0) return

      do 180 ix=1,nxs(6,iex)
  180    if(nint(xss(jxs(13,iex)+ix-1)).eq.mt) goto 190
         return

  190    is = jxs(15,iex)+nint(xss(jxs(14,iex)+ix-1))+1
         if(nint(xss(is-2)).eq.12) goto 210
         ic = ktc(1,iex)-nint(xss(is-1))+1
         if(ic.lt.1) return
         if(rtc(1,iex).ne.0.) goto 200
         if(ic.le.nint(xss(is))) getxs = xss(is+ic)
         return

  200    if(ic.lt.nint(xss(is)))
     &    getxs = xss(is+ic)+rtc(1,iex)*(xss(is+ic+1)-xss(is+ic))
         return

  210    yd = acefcn(is,erg,ln)
         ix = nint(xss(is-1))
         if(ix.lt.0) goto 230
         is = jxs(7,iex)+nint(xss(jxs(6,iex)+ix-1))
         ic = ktc(1,iex)-nint(xss(is-1))+1
         if(ic.lt.1) return
         if(rtc(1,iex).ne.0.) goto 220
         if(ic.le.nint(xss(is))) getxs = xss(is+ic)*yd
         return

  220    if(ic.lt.nint(xss(is)))
     &    getxs = (xss(is+ic)+rtc(1,iex)
     &          * (xss(is+ic+1)-xss(is+ic)))*yd
         return

  230    if(jxs(21,iex).eq.0) return
         is = jxs(21,iex)+1
         ic = ktc(1,iex)+1-nint(xss(is-1))
         if(ic.ge.1)
     &    getxs = (xss(is+ic)+rtc(1,iex)
     &          * (xss(is+ic+1)-xss(is+ic)))*yd

*-----------------------------------------------------------------------

      return
      end

