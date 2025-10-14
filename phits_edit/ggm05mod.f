************************************************************************
*                                                                      *
      module SEGMENT
*                                                                      *
*       Main algorithm of the advanced multibody event generator mode  *
*       modified by T.Ogawa on 2014/04/16                              *
*                                                                      *
***********************************************************************

      implicit double precision(a-h, o-z)

*-----------------------------------------------------------------------

      parameter (ep     =  0.000001)
      parameter (isampl = 100    )

*     Arrays
      double precision, private :: prob(isampl) ! event likelyhood
      double precision, private :: pr(isampl) ! residue momentum
      double precision, private :: presx(isampl) ! residue momentum component
      double precision, private :: presy(isampl)
      double precision, private :: presz(isampl)
      double precision, private :: eres(isampl) ! residue excitation energy
      double precision, private :: ensec(10, isampl) ! particle kinetic energy. 10 particle emission can be handled
      double precision, private :: ensav(10, isampl) ! particle kinetic energy in the CM frame. Used for probability correction
      double precision, private :: ansec(10, isampl) ! particle angular cosine
      double precision, private :: bnsec(10, isampl) ! particle azimuthal angle
      double precision, private :: pnsec(10, isampl) ! particle momentum

      integer, save, private          :: np, lf, jj, lw, ln, lflg
      double precision, save, private :: bb, r1, rr
      double precision prfinx, prfiny, prfinz, eresfin, enfin(10),
     & anfin(10), bnfin(10)

!$OMP THREADPRIVATE(prob, pr, presx, presy, presz, eres)
!$OMP THREADPRIVATE(ensec, ensav, ansec, bnsec, pnsec)
!$OMP THREADPRIVATE(prfinx, prfiny, prfinz, eresfin)
!$OMP THREADPRIVATE(enfin, anfin, bnfin)
!$OMP THREADPRIVATE(np, lf, jj, lw, ln, lflg, bb, r1, rr)

      contains


************************************************************************
*                                                                      *
       subroutine e_a_bsmpl(i, ns, q, ia, ka, id, kd, pro, enmean,
     &  estdev, ein, ebind, enmax)
*                                                                      *
*      Sample e, a, b at random from corresponding cross-sections      *
*                                                                      *
*      input parameters :                                              *
*                                                                      *
*      i      : event candidate index                                  *
*      q      : reaction Q value                                       *
*      ia     :                                                        *
*      ka     :                                                        *
*      id     :                                                        *
*      kd     :                                                        *
*      ns     : the number of secondary particles                      *
*                                                                      *
*      output parameters :                                             *
*                                                                      *
*      pro   ; event likelihood                                        *
*      enmean; mean secondary particle energy                          *
*      estdev; secondary particle energy standard deviation            *
*      enmax ; secondary particle maximum energy                       *
*      ensec ; outgoing particle energy                                *
*      ansec ; outgoing particle direction cosine -> cos(theta)        *
*      bnsec ; outgoing particle azimuthal angle  -> 0 - 2pi           *
*      lflg  ; flag to remember if law is 9 -> exit SEGMENT because    *
*                 law 9 (evaporation spectrum ) do not need SEGMENT    *
*                                                                      *
*      Note: SEGMENT is currently applicable for neutron emission,     *
*          however, this subroutine is extendable for other particles  *
*                                                                      *
************************************************************************
      use NGSDATAMOD, only : bindeg
      use GGMARRAYMOD !2020ASTOM
      implicit double precision(a-h, o-z)

*-----------------------------------------------------------------------
      parameter ( rpmass = 938.27, rnmass = 939.58 )

      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'

*-----------------------------------------------------------------------

      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)


* Initialization
      do j = 1, ns
       ensec(j,i)  = 0.d0
       ansec(j,i)  = 0.d0
       bnsec(j,i)  = 0.d0
      enddo

* Hereafter, this subroutine is copy of subroutine "xstcas" and "xsrtbl"

*-----------------------------------------------------------------------
*     select the law.
*-----------------------------------------------------------------------

      IF(i .eq. 1) then

* Initialization
      np = 0
      lf = 0
      jj = 0
      lw = 0
      ln = 0
      r1 = 0.d0
      bb = 0.d0

         nx = 0
   10    t1 = unirn(dummy)
         n = id+kd
         goto 25

   20    n = id+nint(xss(n-1))
   25    if(nint(xss(n-1)).eq.0) goto 30
         t1 = t1-acefcn(n+2,erg,ln)
         if(t1.ge.0.d0) goto 20

   30    lw = nint(xss(n))
         iw = id-1+nint(xss(n+1))

*-----------------------------------------------------------------------
*     Law of (n,xn) reactions are
*
* >>>>>  law 4   -- continuous erg tabular distribution.
* >>>>>  law 44  -- kalbach-87 correlated formalism.
* >>>>>  law 61  -- correlated tab energy-angle dist.
*-----------------------------------------------------------------------

         IF(lw .eq. 9) then
          lflg = 1
          return
           return
         ENDIF

!  xsrtbl starts from here
         nr = nint(xss(iw))
         ie = iw+2*nr+1
         nf = nint(xss(ie))
         ln = 2*(nr+1)+nf
         r = 0.d0

*-----------------------------------------------------------------------
*     use extreme value if energy is off either end of the table.
*-----------------------------------------------------------------------

         if(erg.ge.xss(ie+nf)) then
          rr = 1.0d0
          goto 170
         elseif(erg.le.xss(ie+1)) then
          rr = 1.0d0
          goto 180
         ENDIF

*-----------------------------------------------------------------------
*     binary search for the location of the energy in the table.
*-----------------------------------------------------------------------

         jb = ie+1
         ib = ie+nf
  110    if(ib-jb.eq.1) goto 130
         ih = (jb+ib)/2
         if(erg.lt.xss(ih)) goto 120
         jb = ih
         goto 110
  120    ib = ih
         goto 110
  130    ii = jb-ie
         rr = erg/xss(jb)

*-----------------------------------------------------------------------
*     calculate interpolation fraction r unless int=1 (histogram).
*-----------------------------------------------------------------------

         if(nr.eq.0) goto 160

      do 140 n = 1, nr
  140    if(ib-ie.le.nint(xss(l+n))) goto 150
         n = nr
  150    if(nint(xss(l+nr+n)).eq.1) goto 185
  160    if(erg-xss(jb).lt.1.d-6*(xss(ib)-xss(jb))) goto 185
         r = (erg-xss(jb))/(xss(ib)-xss(jb))
         goto 185

  170    ii = nf
         goto 185
  180    ii = 1
         goto 185

*--xsrtbl ends here-----------------------------------------------------

  185    nr = nint(xss(iw))    ! ii: incident energy bin number, r: interpolation ratio,
         lb = iw-1+ln+ii
         lc = id+nint(xss(lb))
         ld = lc
         lf = lc

         np = int(xss(lc)+ep)
         mp = np
         jj = nint(xss(lc-1))
         nd = 0
         if(jj.lt.9999) nd = jj/10
         ld = id+nint(xss(lb+1))
         mp = int(xss(ld)+ep)
         t1 = xss(lc+nd+1)+r*(xss(ld+nd+1)-xss(lc+nd+1))
         t2 = xss(lc+np)+r*(xss(ld+mp)-xss(lc+np))

         ra = unirn(dummy)

         if(ra.ge.r) goto 190
         lf = ld
         jj = nint(xss(ld-1))

         if(jj.gt.10000) nd = 0
 190     if(jj.eq.9999) return
         jj = jj-nd*10
         r1 = unirn(dummy)

*-----------------------------------------------------------------------
*        check for a hit in the discrete part.
*        use histogram or corresponding-point interpolation.
*-----------------------------------------------------------------------

         if(nd.eq.0) goto 197


      do 195 ih = 1, nd
         cc = xss(lc+2*np+ih)+r*(xss(ld+2*mp+ih)-xss(lc+2*np+ih))
  195    if(cc.ge.r1) goto 196

*         if(np.eq.nd.and.mp.eq.nd) goto 252  ! strange operation?  t is undefined yet?
         goto 197

  196    continue

*-----------------------------------------------------------------------
*        handle a hit in the continuous part.
*        use histogram or unit-base interpolation inside table.
*        use scaled interpolation between tables.
*-----------------------------------------------------------------------

  197    np = int(xss(lf)+ep)
         nss = nd+1
         if(nd.ne.0) r1 =1.d0-(1.d0-r1)*(1.d0-xss(lf+2*np+nd))/(1.d0-cc)

*-----------------------------------------------------------------------
*        adjust for the energy cutoff if necessary.
*        see subroutine expung for law 4/44 distribution. 0 < wc < 1
*-----------------------------------------------------------------------

         if(jj.lt.9999) goto 200
         jj = jj/10000
         nss = nint(xss(lf-1))-jj*10000
         wc = 2.d0*(xss(lf)-np)
         r1 = 1.d0-wc*(1.d0-r1)
         wgt = wgt*wc

*-----------------------------------------------------------------------
*        find the energy bin in the continuous part.
*-----------------------------------------------------------------------


! calculation of mean secondary energy
  200     enmean = 0.d0
          estdev = 0.d0
          do iii = 1, np - 1
           enmean = enmean + xss(lf+1 + iii) *
     &       (xss(lf+1 + 2 * np + iii) - xss(lf+1 + 2 * np + iii - 1))
          enddo
          do iii = 1, np - 1
           estdev = estdev + (enmean - xss(lf+1 + iii) )**2 *
     &   ( xss(lf+1 + 2 * np + iii) - xss(lf+1 + 2 * np + iii - 1) )
          enddo
           estdev = sqrt(estdev)  ! standard deviation of secondary particle energy

           enmax = xss(lf + 1 + np - 1) ! Spectrum high energy end

        ENDIF

*--------sampling rule------------------------------------------

         pro = 1.d0 ! initialize
         do j = 1, ns

          indx = int( (np - 1) * unirn(dummy) )  ! random integer 0 ~ np - 2
          er = unirn(dummy)
          ln = lf + 1 + indx  ! pointer for energy
          ic = ln + 2 * np  ! pointer for cumm. pdf
          ib = ic + 1
          fa = xss(ln + np) ! X-section
          ea = xss(ln)      ! secondary part bin min energy

! secondary part energy

          IF(xss(lf - 1) .eq. 1 .or. xss(ln + np + 1) - xss(ln + np)
     &     .eq. 0) then
           ensec(j,i) = xss(ln) * (1.d0 - er) + xss(ln + 1) * er ! histogram
          ELSEIF(xss(lf - 1) .eq. 2) then
           ensec(j,i) = xss(ln) +
     &       (xss(ln + 1) - xss(ln))/(xss(ln + np + 1) - xss(ln + np))
     &       * ( - xss(ln + np) + sqrt( xss(ln + np)**2 + unirn(dummy)
     &       * (xss(ln + np + 1)**2 - xss(ln + np)**2) ))  ! lin-lin interpolation
          ENDIF

          ensec(j,i) = ensec(j,i) * rr ! correction for incident energy. Simple factorization. Caution: This operation is suspicious! Ask nuclear data professionals
          pro = pro * ( xss(ln + np + 1) + xss(ln + np) )
     &       * (xss(ln + 1) - xss(ln) )  ! probability. Multiplied by energy bin width

*---Energy sampling----------------------------------------------

          IF(jj .eq. 1) goto 240
           bb = (xss(ln+np+1)-fa)/(xss(ln+1)-ea)  ! gradient of the function
          IF(bb .eq. 0.d0) goto 240
           t = ea+(dsqrt(dmax1(zero,fa**2+2.d0*bb*(r1-xss(ic))))
     &            -fa)/bb

          IF(lw .ne. 44) goto 250  ! parameters to determine angle in case e-a are correlated
           fb = (ensec(j,i) -xss(ln))/(xss(ln+1)-xss(ln))
           tpd(1) = xss(ic+np)+fb*(xss(ic+1+np)-xss(ic+np)) ! R factor for Kalbach systematics
           tpd(2) = xss(ic+2*np)+fb*(xss(ic+1+2*np)-xss(ic+2*np)) ! A factor for Kalbach systematics
           goto 250
  240     lw = lw ! dummy line to commentout following equation sometimes 1/zero
          IF(lw .ne. 44) goto 250
           tpd(1) = xss(ic+np)
           tpd(2) = xss(ic+2*np)
  250     lw = lw ! dummy line to commentout following equation sometimes 1/zero

*-----------------------------------------------------------------

*---Angle sampling--------switching by angular distribution type----------------------------------
          IF(lw .eq. 44) then ! Kalbach-87 formalism. tpd(1) and tpd(2) are used as R and A of Kalbach-87 formalism
           ipsc = 14
           IF(unirn(dummy) .gt. tpd(1)) then
            t1 = (2.d0*unirn(dummy)-1.d0)*sinh(tpd(2))
            ansec(j,i) = log(t1 + sqrt(t1**2 + 1.d0)) / tpd(2)
           ELSE
            r2 = unirn(dummy)
            ansec(j,i) = log(r2* exp(tpd(2))+(1.d0-r2)* exp(-tpd(2)))
     &        / tpd(2)
           ENDIF
          ELSEIF(lw .eq. 61) then ! tabulated angular distribution
           ipsc = 16
           lb = ic+np

           if(jj.ne.1.and.(xss(ib)-r1.lt.r1-xss(ic))) lb = lb+1
           ixcos = 0
           lm = nint(xss(lb))
           if(lm.eq.0) ansec(j,i) = 2.d0 * unirn(dummy) - 1.d0
           if(lm.gt.0) ansec(j,i) = acecos(id,-lm)
          ELSEIF(lw .eq. 4) then
           ansec(j,i) = acecos(ia,ka)
          ENDIF

          bnsec(j,i) = 2.d0 * pi * unirn(dummy)

          ensav(j,i) = ensec(j,i)

*-----------------------------------------------------------------------
*        formulas below are from p. 2 of x-6:res-93-68.
*        these formulas assume two-body kinematics.
*        seamon's formulas specify atomic weight ratios (to neutron)
*        for incident particle (a), awr=gpt(ipt_incident)/gpt(1)
*        for exiting particle (b), awr=gpt(ipt)/gpt(1)
*        for target (a), awr=awn(iex)
*
*        Outgoing direction is converted from CM to Lab frame with consideration for relativistic kinematics.
*-----------------------------------------------------------------------

         ip    = 1 ! incident particle is neutron. See ggm03 "block data blkdat"
         ipt   = 1 ! outgoing particle is neutron. Same as above
         mrtar = mathz * rpmass + mathn * rnmass - bindeg(mathz,mathn)  ! target mass

         a1 = gpt(ip) / gpt(1)
         a2 = gpt(ipt)/ gpt(1)
         a3 = mrtar   / gpt(1) ! awn(iex)
         t2 = erg*a1*a2/(a3+a1)**2
         t1 = ensec(j,i)
         t3 = 2.*sqrt(a2*a1*erg*ensec(j,i))*ansec(j,i)/(a3+a1)
         t4 = t1+t2+t3
         s1 = ansec(j,i)*sqrt(ensec(j,i)/t4)
         s2 = sqrt(a1*a2*erg/t4)/(a3+a1)
         ensec(j,i) = t4
*         ansec(j,i) = s1+s2
         ansec(j,i) = min(1.d0, s1+s2)  ! sometimes angular cosine is above 1.0
*-----------------------------------------------------------------

         enddo

*--------------------------------------------------------------

      return
      end subroutine
************************************************************************



**************************************************************************
*                                                                        *
      subroutine mainsmt(ein, ns, ndsig, q, ia, ka, id, kd)
*                                                                        *
*      Sample a event out of 10^4 events to determine angle and energy of*
*      secondary neutrons                                                *
*                                                                        *
*      input parameters :                                                *
*                                                                        *
*      ein     ; Incident particle energy                                *
*      iz      ; target nucleus charge                                   *
*      ia      ; target nucleus mass                                     *
*      ns      ; outgoing particle number                                *
*      ndsig   ; sign denotes angular distribution (Neg: CMS, Pos: Lab)  *
*      ins(i)  ; i-th outgoing particle species (currently neutron only) *
*      n*bin   ; *-th particle spectrum bin = angular bin * energy bin   *
*      ismtflg ; flag to remember if SEGMENT was successful or not       *
*                                                                        *
*      output parameters :                                               *
*                                                                        *
*      enfin1, enfin2....enfin* : *-th outgoing particle energy          *
*      anfin1, anfin2....anfin* : *th particle anglular cosine in CMS frame  *
*      bnfin1, bnfin2....bnfin* : *th particle azimuthal angle in CMS frame  *
*      prfinx,y,z               : nucleus recoil momentum                *
*      eresfin                  : nucleus residual excitation energy     *
**************************************************************************
      use NGSDATAMOD, only : bindeg
      use GGMARRAYMOD !2020ASTOM
      implicit double precision(a-h, o-z)

      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /ismtfl/ ismtflg
!$OMP THREADPRIVATE(/ismtfl/)
      common /emode1/ nevhin, nlowlv, nefiss, ntwidt, mtprec,ipcnt(70)
!$OMP THREADPRIVATE(/emode1/)

      include 'ggmparam.inc'

*-----------------------------------------------------------------------

*     Basic Parameters
      resmas = rmtyp(19, mathz * 1000000 + (mathz + mathn + 1 - ns) ) ! residual mass (MeV)
!      ebind  = bindeg(mathz, mathz + mathn + 1 - ns)  ! T.Sato 2024/09/22, bug fix because bindeg requires neutron number
      ebind  = bindeg(mathz, mathn + 1 - ns)  
     & / dble(mathz + mathn + 1 - ns)
      fe = ein + q  ! available energy
      probsum = 0.d0
      ensec   = 0.d0
      ansec   = 0.d0
      bnsec   = 0.d0
      enfin   = 0.d0
      anfin   = 0.d0
      bnfin   = 0.d0
      prfinx  = 0.d0
      prfiny  = 0.d0
      prfinz  = 0.d0
      eresfin = 0.d0
      ismtflg = 0
      lflg    = 0
      rej     = fe - ebind

      do i = 1, isampl ! Calculate event probability for 10000 different events

       itry         = 0
*-----------------------------------------------------------------------------
* -------------Determination of secondary particles and residual nucleus -----
 100   if(itry .ge. 1000) then
        return ! available phase space was too small to sample
       endif
       p0           = dsqrt( ein**2   + 2.d0 * gpt(1) *  ein) ! primary momentum
       eres(i)      = fe  ! initial excitation energy
       pro          = 0.d0

        call e_a_bsmpl(i, ns, q, ia, ka, id, kd, pro, enmean, estdev,
     &   ein, ebind, enmax) ! sample particle energy and angle
*  e_a_bsmpl  samples following quantities
*        ensec(j,i)  : j-th secondary particle energy
*        ansec(j,i)  : j-th secondary particle direction cosine (-1 ~ 1)
*        bnsec(j,i)  : j-th secondary particle asimuthal angle (0 ~ 2 pi)
*        pro is likelyhood of this event
*        enmean is mean secondary particle energy

        IF(lflg .eq. 1) return ! law is 9 therefore SEGMENT is unnecessary

        IF(i .eq. 1) then
         exmean = fe - ns * enmean  ! mean excitation energy
        ENDIF

        prob(i)  = pro
        presx(i) = 0.d0
        presy(i) = 0.d0
        presz(i) = sqrt( ein**2 + 2.d0 * gpt(1) * ein)

       do j = 1, ns  ! determine recoil and excitation from secondary particle kinematics
        pnsec(j,i) = sqrt( ensec(j,i)**2 + 2.d0 * gpt(1) * ensec(j,i)) ! particle momentum
        presx(i)   = presx(i) -pnsec(j,i) * sqrt(1.d0 - ansec(j,i)**2)
     &  * sin(bnsec(j,i))  ! residual recoil momentum
        presy(i)   = presy(i) -pnsec(j,i) * sqrt(1.d0 - ansec(j,i)**2)
     &  * cos(bnsec(j,i))
        presz(i)   = presz(i) - pnsec(j,i) * ansec(j,i)
        eres (i)   = eres(i)  - ensec(j,i) ! excitation energy (i.e. remaining energy)
       enddo

        pr(i)      = sqrt( presx(i)**2 + presy(i)**2 + presz(i)**2 )

        eres(i)    = eres(i) - ( sqrt(pr(i)**2 + resmas**2 ) - resmas)

*-----------------------------------------------------------------------------
*-----Probability bias ------------------------------------------------------

* Record the event with its probability if energy and momentum conservation is satisfied
       IF(eres(i) .lt. 0.d0 .and. (mtprec .eq. 11 .or. mtprec .eq. 24
     &  .or. mtprec .eq. 41 .or. mtprec .eq. 42) ) then ! to emit charged pariticles, allow exc > ebind
        prob(i) = 0.d0
        IF(i .ne. 1) then
         itry = itry + 1
         goto 100 ! restart sampling
        ENDIF
        cycle                 ! if goto 100, basic data claculation is restarted taking extra CPU time
       ELSEIF(eres(i) .lt. 0.d0 .or. eres(i) .ge. ebind) then ! energy conservation unsatisfied -> event rejected
        prob(i) = 0.d0
        IF(i .ne. 1) then
         itry = itry + 1
         goto 100 ! restart sampling
        ENDIF
        cycle                 ! if goto 100, basic data claculation is restarted taking extra CPU time
       ENDIF

c       Bias for heavy nuclei with spread-out secondary neutron spectra
       IF(estdev .le. 2.d0 .or. exmean - ebind/2.d0 .gt. 5.d-1) then !
        do j = 1, ns
         prob(i) = frejec(ensav(j,i), ebind, fe) * prob(i)   ! Probability. 1 / rejection probability
        enddo
       ELSEIF(estdev .le. 4.d0) then ! bias is smoothly weaken with increase of estdev
        do j = 1, ns
         prob(i) = frejec(ensav(j,i), ebind, fe)**(2.d0 - estdev/2.d0)
     &    * prob(i)   ! Probability. 1 / rejection probability
        enddo
       ENDIF

! probability correction for excitation energy
        prob(i) = prob(i) * exp(5.d-2 * eres(i) * (exmean - ebind/2.d0))

! bias to increase spectrum high-energy end
       IF(ns .eq. 2 .and. exmean - ebind/2.d0 .lt. ein - 20.d0           ! Activate bias if spectra are high energy shifted
     & .and. enmax .ge. ensav(1,i) .and. enmax .ge. ensav(2,i)  ) then   ! Condision for bug prevention
       if(estdev .lt. 1.d-3) then ! 2017/3/6 Ogawa workaround to avoid infinity
        prob(i) = prob(i) * 1.d10
       else
        prob(i) = prob(i) * max(
     & (fermi(ensav(1,i), enmax, 2.d1) * fermi(ensav(2,i), enmax, 2.d1))
     & ** (-6.d0 / estdev), 1.d10)
       endif

       ENDIF

*-----------------------------------------------------------------------------
*------Sampling from calculated events----------------------------------------

       probsum = probsum + prob(i)

      enddo

      if(probsum .le. 0.d0) return ! Neutron failed for 1000 times. Something is wrong

      ransamp = probsum * unirn(dummy) ! sampling point

      do i = 1, isampl

       ransamp = ransamp - prob(i) ! search event to be sampled

       IF(ransamp .le. 0.d0 ) then ! sample this event

        prfinx  = presx(i) ! record on residue
        prfiny  = presy(i)
        prfinz  = presz(i)
        eresfin = eres(i)

        do j = 1, ns ! records on secondaries
         enfin(j)  = ensec(j,i)
         anfin(j)  = ansec(j,i)
         bnfin(j)  = bnsec(j,i)
        enddo

        call finsmt(ns)
        ismtflg = 1   ! flag to remember SEGMENT was successful

        exit

       ENDIF

      enddo

      return
      end subroutine
************************************************************************


**************************************************************************
*                                                                        *
      function frejec(e, ebind, fe)
*                                                                        *
*      Correction function for event rejection probability ( >= 1 )      *
*       (to be upgraded to consider sampling space dimension             *
*           higher than 2 (i.e. (n,2n) reaction )  )                     *
*                                                                        *
*      input parameters :                                                *
*                                                                        *
*      e       ; Sampled secondary particle energy                       *
*      ebind   ; Binding energy of the product                           *
*      fe      ; incident energy - Q value -> free available energy      *
*                                                                        *
**************************************************************************

      implicit double precision(a-h, o-z)

      IF(e .ge. fe - ebind ) then
       frejec = ebind / max( 0.d0, (fe - e)) ! if fe is greater than e, event is invalid
      ELSE
       frejec = 1.d0
      ENDIF

      return
      end function
************************************************************************

**************************************************************************
*                                                                        *
      function fermi(e, emax, fct)
*                                                                        *
*      Correction function for event rejection probability ( >= 1 )      *
*       (to be upgraded to consider sampling space dimension             *
*           higher than 2 (i.e. (n,2n) reaction )  )                     *
*                                                                        *
*      input parameters :                                                *
*                                                                        *
*      e       ; Sampled secondary particle energy                       *
*      ebind   ; Binding energy of the product                           *
*      fe      ; incident energy - Q value -> free available energy      *
*                                                                        *
**************************************************************************

      implicit double precision(a-h, o-z)

      fermi = 2.d0/( exp( (e - emax) / (emax / fct) ) + 1.d0) - 1.0d0

      return
      end function
************************************************************************


**************************************************************************
*                                                                        *
      subroutine finsmt(ns)
*                                                                        *
*      finalization of SEGMENT                                           *
*      1, Assign secondary particle quantities to qclusts array          *
*                                                                        *
*      input parameters :                                                *
*      ns    : secondary particle number                                 *
*                                                                        *
*      enfin1, enfin2....enfin* : *-th outgoing particle energy          *
*      anfin1, anfin2....anfin* : *th particle anglular cosine in CMS frame  *
*      bnfin1, bnfin2....bnfin* : *th particle azimuthal angle in CMS frame  *
*      prfinx,y,z               : nucleus recoil momentum                *
*      eresfin                  : nucleus residual excitation energy     *
*                                                                        *
*                                                                        *
*      output parameters :                                               *
*                                                                        *
*      qclusts : outgoing particle information array                     *
**************************************************************************
      use GGMARRAYMOD !2020ASTOM
      implicit double precision(a-h, o-z)

      include 'param00.inc'
      include 'ggmparam.inc'
      include 'ggsparam.inc'

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /usaves/ usave(3)
!$OMP THREADPRIVATE(/usaves/)

      resmas = rmtyp(19, mathz * 1000000 + (mathz + mathn + 1 - ns) ) ! residual mass (MeV)

* Rotation. Incident direction is (uuu,vvv,www) but sampling was done assuming (0,0,1) incident
      costhe =   usave(3)
      sinthe =   sqrt(1 - usave(3)**2)
      IF(abs(usave(3)) .eq. 1.d0) then ! T.Sato 2024/10/05, bugfix for dir = -1 source
       sinphi = 0.d0
       cosphi = usave(3)
      ELSE
       sinphi =   usave(2) / sqrt(1 - usave(3)**2)
       cosphi = - usave(1) / sqrt(1 - usave(3)**2)
      ENDIF
      rm11 =  costhe * cosphi
      rm12 = -costhe * sinphi
      rm13 =  sinthe
      rm21 =  sinphi
      rm22 =  cosphi
      rm23 =  0.d0
      rm31 = -sinthe * cosphi
      rm32 =  sinthe * sinphi
      rm33 =  costhe
*

* Record information on secondaries
      do i = 1, ns

       nclst = nclst + 1
       iclust(nclst) = 2 ! ityp of this particle is 2 (neutron)

       jclust(0,nclst) = 0
       jclust(1,nclst) = 0
       jclust(2,nclst) = 1
       jclust(3,nclst) = 2
       jclust(4,nclst) = 0
       jclust(5,nclst) = 0
       jclust(6,nclst) = 1
       jclust(7,nclst) = 2112
       jclust(8,nclst) = 0

       rms = gpt(1)  ! neutron mass in MeV
       rmg = rms / 1000.d0
       pre = dsqrt( enfin(i)**2 + 2.d0 * rms * enfin(i) ) / 1000.d0
       ett = dsqrt( pre**2 + rmg**2 )

! Momentum before roration (rotated to sample neutrons assuming (0,0,1) incident )
       pxbr = pre * dsqrt(1.d0 - anfin(i)**2 ) * dsin(bnfin(i))
       pybr = pre * dsqrt(1.d0 - anfin(i)**2 ) * dcos(bnfin(i))
       pzbr = pre * anfin(i)

       qclust(0,nclst)  = 0.d0
! Rotation to the original coordinate system
       qclust(1,nclst)  = pxbr * rm11 + pybr * rm21 + pzbr * rm31
       qclust(2,nclst)  = pxbr * rm12 + pybr * rm22 + pzbr * rm32
       qclust(3,nclst)  = pxbr * rm13 + pybr * rm23 + pzbr * rm33
       qclust(4,nclst)  = ett
       qclust(5,nclst)  = rmg
       qclust(6,nclst)  = 0.d0
       qclust(7,nclst)  = enfin(i)
       qclust(8,nclst)  = wgt / wg0
       qclust(9,nclst)  = 0.d0   ! time is zero. This should overwritten with the reaction time ?
       qclust(10,nclst) = 0.d0
       qclust(11,nclst) = 0.d0
       qclust(12,nclst) = 0.d0
      enddo

* Record of residual nucleus
       nclst = nclst + 1

       IF( mathz .gt. 2) then
        iclust(nclst) = 0 ! ityp of this particle is 0 (nucleus)
       ENDIF

       jclust(0,nclst) = 0 ! angular momentum. to be considered in the future
       jclust(1,nclst) = mathz  ! neutron emission is considered
       jclust(2,nclst) = mathn + 1 - ns
       jclust(3,nclst) = 2
       jclust(4,nclst) = 0
       jclust(5,nclst) = mathz
       jclust(6,nclst) = mathz + mathn + 1 - ns
       jclust(7,nclst) = 1000000 * mathz + (mathz + mathn + 1 - ns)
       jclust(8,nclst) = 0

       qclust(0,nclst)  = 0.d0
! Rotation to the original coordinate system
       qclust(1,nclst)  = (prfinx * rm11 + prfiny * rm21 + prfinz *
     &  rm31) / 1.d3
       qclust(2,nclst)  = (prfinx * rm12 + prfiny * rm22 + prfinz *
     & rm32) / 1.d3
       qclust(3,nclst)  = (prfinx * rm13 + prfiny * rm23 + prfinz *
     & rm33) / 1.d3
       qclust(4,nclst)  = dsqrt( resmas **2 + (prfinx**2 + prfiny**2
     &  + prfinz**2) ) /1.d3
       qclust(5,nclst)  = resmas /1.d3
       qclust(6,nclst)  = eresfin
       qclust(7,nclst)  = dsqrt( (prfinx**2 + prfiny**2 + prfinz**2)
     & + resmas **2 )  - resmas
       qclust(8,nclst)  = wgt / wg0
       qclust(9,nclst)  = 0.d0   ! time is zero. This should overwritten with the reaction time ?
       qclust(10,nclst) = 0.d0
       qclust(11,nclst) = 0.d0
       qclust(12,nclst) = 0.d0

*-----------------------------------------------------------------------

                           nevhin = 0
                           nlowlv = 0 ! related to low energy discrete level ?
                           nefiss = 0

                     call nevap(2)

                           nevhin = 0
                           nlowlv = 0
                           nefiss = 0

*-----------------------------------------------------------------------

      return
      end subroutine
************************************************************************



      end module SEGMENT
************************************************************************

