************************************************************************
*                                                                      *
      subroutine dexgam(aaa,zzz,exe,spin,isomle)
*                                                                      *
*       calculate gamma emission from deexciting nucleus               *
*       modified by K.Niita on 05/04/2000                              *
*                                                                      *
*     input:                                                           *
*                                                                      *
*        aaa :    mass number of nucleus                               *
*        zzz :    charge number of nucleus                             *
*        exe :    excitation energy of nucleus (MeV)                   *
*       spin :    spin of nucleus (hbar)                               *
*                                                                      *
*     output (module variable / array):                                *
*                                                                      *
*        nph    :   number of photon                                   *
*        eph(i) :   energy of photon                                   *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma

      include 'gamlib1.inc' !FURUTA20201007
      include 'gamlib2.inc' !FURUTA20201007

*-----------------------------------------------------------------------

      dimension maxa(1:100)
      dimension mina(1:100)

      data ncount/5*0/ !FURUTA

*    Mass of nucleus simulated by EBITEM de-excitation routine
      data mina /1,   3,   3,   5,   6,   8,  10,  12,  14,  16,
     &          18,  19,  21,  22,  24,  26,  28,  30,  32,  34,
     &          36,  38,  40,  42,  44,  45,  47,  48,  52,  54,
     &          56,  58,  60,  64,  67,  69,  71,  73,  76,  78,
     &          81,  83,  85,  87,  89,  91,  93,  95,  97,  99,
     &         103, 105, 107, 109, 112, 114, 116, 119, 121, 124,
     &         126, 128, 130, 133, 135, 138, 140, 142, 144, 148,
     &         150, 153, 155, 157, 159, 161, 164, 166, 169, 171,
     &         176, 178, 184, 186, 191, 193, 199, 201, 206, 208,
     &         212, 217, 219, 228, 230, 232, 234, 237, 239, 241/
      data maxa /7,  10,  13,  16,  21,  23,  25,  28,  31,  34,
     &          37,  40,  43,  45,  47,  49,  51,  53,  56,  58,
     &          61,  63,  66,  68,  71,  74,  76,  79,  82,  85,
     &          87,  90,  92,  95,  98, 101, 103, 107, 109, 112,
     &         115, 117, 120, 124, 126, 128, 130, 133, 135, 138,
     &         140, 143, 145, 148, 151, 153, 155, 157, 159, 161,
     &         163, 165, 167, 169, 171, 173, 175, 177, 179, 181,
     &         185, 189, 192, 194, 198, 202, 204, 206, 210, 216,
     &         218, 220, 224, 227, 229, 231, 233, 235, 237, 239,
     &         241, 243, 245, 247, 249, 252, 254, 256, 258, 260/


*    Charge of nucleus simulated by EBITEM
      data minz /3/
      data maxz /97/

* ogw end

*-----------------------------------------------------------------------
*     initialization move to setpar in read00.f
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*     initialize gamcal output option
*-----------------------------------------------------------------------

            icc = 0

            oldexe = exe
            wtex   = 1.0

*-----------------------------------------------------------------------
*     calc. gamma emission from deexciting nucleus
*-----------------------------------------------------------------------

      ierr = 1

* Gamma de-excitation by EBITEM for nuclei in particular mass, charge range
      IF(idint(zzz) .ge. minz   .and. idint(zzz) .le. maxz) then
        IF( idint(aaa) .ge. mina(idint(zzz) ) .and.
     &      idint(aaa) .le. maxa(idint(zzz) ) .and.
     &      abs(igamma) .ge. 2) then

          call gamcal2(idint(aaa), idint(zzz), exe, spin, isomle, ierr)
        endif
      ENDIF

      IF(ierr .eq. 1) then
         isomle = 0
         call gamcal(aaa,zzz,exe,icc,ecv,wtex)
      ENDIF

*-----------------------------------------------------------------------
      return
      end




************************************************************************
*                                                                      *
      subroutine gamcal2(iares, izres, eexe, totJ, isomle, ierr)
*                                                                      *
*       isomeric transition calculation from excited states to         *
*         the ground state or the metastable states                    *
*       created by T.Ogawa on 09/24/2012                               *
*                                                                      *
*       ENSDF Based Isomeric Transition/ isomEr production Model       *
*               -->  EBITEM                                            *
*                                                                      *
*        parameters (transferred as arguments) :                       *
*                                                                      *
*        Input :                                                       *
*        iares  : mass number of residual nucleus                      *
*        izres  : charge number of residual nucleus                    *
*        eexe  : excitation energy of the residual nucleus (MeV)       *
*        spin  : spin of the residual nucleus (hbar)                   *
*                                                                      *
*        Output :                                                      *
*           explicit argument transfer                                 *
*              isomle  : final energy level of the product             *
*                  0 : ground, 1 : lowest isomer, 2 : 2nd lowest isomer*
*                                                                      *
*                                                                      *
*           transfer as module variable & arrays                       *
*                                                                      *
*              nph       : the number of ejected photons               *
*              eph(i)    : energy of ejectile (MeV)   i = 1 - nph      *
*              kfejec(i) : kf code of ejectile (MeV)  i = 1 - nph      *
*              ejectile is gamma photons or conversion electron        *
*                                                                      *
************************************************************************

* call level energy, branching ratios ...etc.
      use levdat
      use NGSDATAMOD, only : bindeg
      implicit double precision(a-h, o-z)

      include 'err.inc'

*-----------------------------------------------------------------------
* Target material and reaction mode
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /nrfmem/ spis, spgr, levabs, lflgnrf, mpole
!$OMP THREADPRIVATE(/nrfmem/)

      data infflg /0/

      logical lflg

*-----------------------------------------------------------------------
*---- Initialization ---------------------------------------------------

      nsta   =-1    ! >= 0 : current level number.  -1 : undetermined.
      negm   = 0
      nph    = 0
      nparit = 0 ! parity
      eph    = 0.d0
      kfejec = 0
      isomle = 0

*-----------------------------------------------------------------------

      call levset(iares, izres, lflg)

* Check if gamcal2 works or not
      IF(.not. lflg) then
       return
      ELSEIF(ubound(elevel,1) .lt. 1) then ! 2014/8/26  check ubound. just in case nlevel !=0 but elevel array does not have 1st element
       call levUNset
       return
      ELSE
       ierr = 0
      ENDIF

*---------------------------------------------------------------

c Event generator mode by neutron
      IF(izres .eq. mathz .and. jcoll .eq. 10) then
       IF(abs(eexe - (bindeg(izres, iares - izres) - bindeg(izres, iares
     &- izres - 1))) .le. 1.d-1 .and. iares - izres .eq. mathn + 1) then ! neutron capture reaction
        negm = 1
       elseif(iares - izres .eq. mathn) then ! inelastic scattering
        negm = 2
       endif
      ENDIF

* 0-0 Special treatment for NRF reaction
      IF(lflgnrf .eq. 1 .and. levabs .ne. 0) then ! level
       nsta = abs(levabs) ! levabs is not meaningful for nuclei except Pb

* 0-1, If neutron capture and event generator mode, follow the special decay scheme
      ELSEIF( negm .eq. 1) then
       if(lev_nth_cap(izres, iares - izres) .ne. 0) then ! use RIPL data for deexcitation
        nsta = lev_nth_cap(izres, iares - izres)
       else
        call neufdc(iares, izres, eexe, nsta)
       endif
       negm = 0

* 0-2, If neutron inelatic and event generator mode, get on the level
      ELSEIF( negm .eq. 2) then
       do ii = 0, ubound(elevel, 1)
        if( abs(eexe - elevel(ii)) .le. 1.d-2 ) then
         nsta = ii
         exit
        endif
       enddo
       negm = 0

* 0-3, Goto ground state if below first excitation level. (GEM doesn't care about discrete level structure)
      ELSEIF( eexe .lt. (elevel(1) + 1.d-10)  .or. lflgnrf .eq. 1) then  ! Excited state below 1st state or that reached by NRF
       nph         = 1
       eph(1)      = eexe
       kfejec(nph) = 22 ! photon
       eexe        = 0.d0
       nsta        = 0

* 0-4, If (first gamma and excitation energy < 3000 keV) or mass < 40, levels are discrete. Follow ENSDF level based on theory
      ELSEIF( eexe .le. min(3.d0, elevel(ubound(elevel,1))) .or.
     & iares .le. 40 ) then ! elevel(1) .ge. 3.d0 .or.

      call firdec(iares, izres, eexe, nsta, totJ)

      endif


* Calculation of isomeric transition
      do iloop = 1, 1000

       if(nsta .ne. -1) then              ! Exception handling
         if(totJ .eq. -1.d0)   totJ   = spinpar(nsta, 1)
         if(nparit .eq. 0  )   nparit = parit(nsta)
         if(ndch(nsta) .eq. 0 .and. nsta .ne. 1) nsta   = -1  ! state whose deexcitation is unknown
       endif

       IF(nsta .eq. 0 .or. eexe .eq. 0.d0) then ! in case, ground is reached
        continue

* Select decay channel

* 1, If nucleus is not on the RIPL level. Decay based on theory.
       elseIF(nsta .eq. -1) then

        iflg = 0
        call therdec(iares, izres, eexe, nsta, totJ, nparit, iflg,
     & sumwid)   ! totJ is spin after decay if the destination level is unresolved

* 2, If excitation energy is below 3000 keV and on a discrete level, follow ENSDF.
       ELSE

        call decay_disc_lev(nsta,eexe,totJ,nparit,izres)

       ENDIF

c Ground state production check

       IF(nsta .eq. 0 .or. eexe .eq. 0.d0)  then
        nsta   = 0
        eexe   = 0.d0
        isomle = 0
        exit
* Isomer production check
       ELSEIF(nsta .eq. isolev1(izres, iares-izres)) then
         if(abs(igamma) .eq. 3) then
           isomle = 1
           exit
         elseif(abs(igamma) .eq. 2 .and. ndch(nsta) .eq. 0) then
           call decay_disc_lev(nsta,eexe,totJ,nparit,izres)
           !nsta = 0
           !nph  = nph + 1
           !eph(nph) = eexe - elevel( nsta ) ! - electron binding energy (to be determined)   ! internal conversion electron
           !kfejec(nph) = 22
           !eexe = elevel( nsta )
           !totJ = spinpar(nsta , 1)
         endif
       ELSEIF(nsta .eq. isolev2(izres, iares-izres)) then
         if(abs(igamma) .eq. 3) then
           isomle = 2
           exit
         elseif(abs(igamma) .eq. 2 .and. ndch(nsta) .eq. 0) then
           call decay_disc_lev(nsta,eexe,totJ,nparit,izres)
           !nsta = 0
           !nph  = nph + 1
           !eph(nph) = eexe - elevel( nsta ) ! - electron binding energy (to be determined)   ! internal conversion electron
           !kfejec(nph) = 22
           !eexe = elevel( nsta )
           !totJ = spinpar(nsta , 1)
         endif
       ELSEIF(nsta .eq. 1) then ! Non-isomer 1st excitation level
           call decay_disc_lev(nsta,eexe,totJ,nparit,izres)
           !nsta = 0
           !nph  = nph + 1
           !eph(nph) = eexe ! emit as gamma
           !kfejec(nph) = 22
           !eexe = 0.d0
           !totJ = spinpar(nsta , 1)
       ENDIF

      end do

       if(iloop .eq. 1001 .and. infflg .eq. 0) then ! 2015/8/19 ogawa. Infinite loop in gamma deexcitation simulation.
        write(ErrCha,*)'Warning: deexcitation path not found in EBITEM'
        ErrID = 'L:313/R:gamcal2/F:gamlib.f' !W00_003_001
        call ErrWrite(ErrID,ErrCha)
        infflg = 1
       endif

       call levUNset

*-----------------------------------------------------------------------
      return
      end

************************************************************************
*                                                                      *
      subroutine decay_disc_lev(nsta,eexe,totJ,nparit,izres)
*                                                                      *
*       isomeric transition and internal conversion based on RIPL-3    *
*       created by T.Ogawa on 12/19/2022                               *
*                                                                      *
*        parameters (transferred as arguments) :                       *
*                                                                      *
*        Input :                                                       *
*              izres  : charge number of residual nucleus              *
*                                                                      *
*        Input and output :                                            *
*              eexe  : excitation energy of the residual nucleus (MeV) *
*              nsta  : current state id number                         *
*              eexe  : current excitation energy (MeV)                 *
*              totJ  : spin of the residual nucleus (hbar)             *
*              nparit: parity of the residual nucleus                  *
*                                                                      *
*        Output :                                                      *
*           transfer as module variable & arrays                       *
*                                                                      *
*              nph       : the number of ejected photons               *
*              eph(i)    : energy of ejectile (MeV)   i = 1 - nph      *
*              kfejec(i) : kf code of ejectile (MeV)  i = 1 - nph      *
*              ejectile is gamma photons or conversion electron        *
*                                                                      *
************************************************************************

* call level energy, branching ratios ...etc.
      use levdat
      use ELEDATAMOD, only : lng_rel, eng_rel, ktp_rel, do_atom_relax,
     &                       pot_elem

      implicit double precision(a-h, o-z)

      probdc = unirn(dummy)
      probra = 0.d0
      nsta_init = nsta

      do jcont = 1, ndch(nsta)
       probra = bratio(nsta, jcont) + probra
       If( probdc .le. probra .or. jcont .eq. ndch(nsta) ) then
           nsta = nlevdn( nsta, jcont )
           nph = nph + 1
           if(braicc(nsta_init, jcont) .gt. unirn(dummy) ) then ! branching of gamma/IC electron 
              do i = 1, 31
                edif = eexe - elevel(nsta) - pot_elem(i,izres)*1.d-6
                if(edif .gt. 0) then  ! energy conserved
                    eph(nph)    = edif   ! internal conversion electron
                    kfejec(nph) = 11  ! electron
                    call do_atom_relax(izres,i) ! subsequent atomic deexcitation
                    do j = 1, lng_rel
                        nph         = nph + 1
                        eph(nph)    = eng_rel(j)
                        kfejec(nph) = ktp_rel(j)
                    enddo
                    eexe = elevel( nsta )
                    totJ = spinpar(nsta , 1)
                    nparit = parit(nsta)
                    return
                endif
              enddo
           endif
           eph(nph) = eexe - elevel( nsta ) ! gamma decay
           kfejec(nph) = 22 ! photon
           eexe = elevel( nsta )
           totJ = spinpar(nsta , 1)
           nparit = parit(nsta)
           return
       ENDIF
      end do

      If(nsta .eq. 1) then ! if decay is not given in RIPL
        nsta = 0
        nph  = nph + 1
        eph(nph) = eexe - elevel( nsta )
        kfejec(nph) = 22
        eexe = elevel( nsta )
        totJ = spinpar(nsta , 1)
      endif

      end subroutine

************************************************************************
*                                                                      *
      subroutine gemgamcomp(iares, izres, u, totJ, sumwid)
*                                                                      *
*       Calculation of gamma decay total width to be compared with     *
*         hadronic decay width in GEM                                  *
*       created by T.Ogawa on 10/14/2016                               *
*                                                                      *
*                                                                      *
*        parameters (transferred as arguments) :                       *
*                                                                      *
*        Input :                                                       *
*        iares  : mass number of residual nucleus                      *
*        izres  : charge number of residual nucleus                    *
*        totJ   : current angular momentum                             *
*        u      : excitation energy of the residual nucleus (MeV)      *
*                                                                      *
*        Output :                                                      *
*        sumwid  : decay width summed over all the levels              *
*                                                                      *
************************************************************************

* call level energy, branching ratios ...etc.
      use levdat

      implicit double precision(a-h, o-z)
      logical lflg

*-----------------------------------------------------------------------
*---- Initialization ---------------------------------------------------

      nsta   = 0   ! current level number

*-----------------------------------------------------------------------

      call levset(iares, izres, lflg)

* Check if gamcal2 works or not
      IF( .not. lflg) then ! 2014/8/26  check ubound. just in case nlevel !=0 but elevel array does not have 1st element
       sumwid = 0.d0
       return
      ENDIF

       nparit = 0
       iflg   = 1
       call therdec(iares,izres,u,nsta,totJ,nparit,iflg,sumwid) ! spin is specified in totJ_init subroutine

       call levUNset

      end subroutine


CTAT/ADD TO READ ORIGINAL DATA FILE ------------------------------( TO )
C
C======================================================================
C======================================================================
      subroutine conprb(e,e1x,iax,isp0x)
      implicit real*8 ( a-h, o-z )
      include 'gamlib1.inc' !FURUTA20201007
      include 'gamlib2.inc' !FURUTA20201007
C
      dimension fms(lms), arg(lms)
      dimension j0(lms), aux(lms)
C
      psjsum  = 0.0
      pnorm   = 0.0
      msp0    = mod(iax,2)-1
      u       = e1x-de
      if (u.le.0.0) return
      conx    = 0.7104*sqrt(ax*u)*cra(iax)**2
      con1    = log(4.0/conx)
      fm      = msp0
      dele    = e-e1x
      nup     = 0.5*(sqrt(conx*(con1+200.0))-fm)+1.0
      nup     = min0(lms,nup)
      do 10 i=1,nup
         fms(i) = fm+2*i
   10 continue
      do 20 i=1,nup
         arg(i) = con1-fms(i)**2/conx
   20 continue
      do 30 i=1,nup
         psj(i) = exp(arg(i))
   30 continue
      do 40 i=1,nup
         psj(i) = psj(i)*fms(i)
   40 continue
      do 50 i=1,nup
         pnorm = pnorm+psj(i)
   50 continue
      do 60 i=1,nup
         j0(i) = max0(iabs(msp0+i+i-isp0x)/2,1)
   60 continue
      if (dele.gt.0.0) then
         do 70 i=1,nup
            aux(i) = psj(i)*fj(j0(i))*(c1*cra(iax)*dele)**(2*j0(i)-2)
   70    continue
      else
         do 90 i=1,nup
            aux(i) = 0.
            if (j0(i).eq.1) aux(i) = psj(i)
   90    continue
      endif
      if (isp0x.eq.1) aux(1) = 0.
      psmax = aux(1)
      do 110 i=2,nup
         psmax = dmax1(psmax,aux(i))
         if (aux(i).lt.0.001*psmax) then
            nup = i
            goto 130
         endif
  110 continue
  130 continue
      nstat = nup
      do 140 i=1,nup
         psjsum = psjsum+aux(i)
         psj(i) = psjsum
  140 continue
      if (nup.lt.lms) then
         do 150 i=nup,lms
            psj(i) = psjsum
  150    continue
      endif
      if (pnorm.le.0.0) return
      do 170 i=1,lms
         psj(i) = psj(i)/pnorm
  170 continue
      psjsum = psj(lms)
      return
      end
C======================================================================
C======================================================================
      function fint(xl,xh)
      implicit real*8 ( a-h, o-z )
      include 'gamlib1.inc' !FURUTA20201007
      include 'gamlib2.inc' !FURUTA20201007
C
C     PERFORMS INTEGRAL OF EXP(X)/X**2 FROM XL TO XH
C     TIMIMG: FROM 12.8 TO 23.4 MICROSECONDS ON CDC7600
C     AND FROM 7.0 TO 12.3 MICROSECONDS ON CRAY
C     ACCURACY: BETTER THAN .05 PERCENT
C
      dimension tab1(66), tab2(80)
      dimension func(2)

C
      data tab1 /0.50000, 0.50847, 0.51724, 0.52631, 0.53570, 0.54542,
     &           0.55549, 0.56593, 0.57675, 0.58797, 0.59962, 0.61171,
     &           0.62427, 0.63732, 0.65088, 0.66499, 0.67966, 0.69493,
     &           0.71084, 0.72740, 0.74467, 0.76268, 0.78146, 0.80106,
     &           0.82153, 0.84290, 0.86525, 0.88861, 0.91304, 0.93861,
     &           0.96539, 0.99343, 1.02282, 1.05363, 1.08595, 1.11987,
     &           1.15547, 1.19287, 1.23218, 1.27350, 1.31696, 1.36269,
     &           1.41083, 1.46154, 1.51497, 1.57129, 1.63069, 1.69337,
     &           1.75953, 1.82939, 1.90320, 1.98121, 2.06370, 2.15096,
     &           2.24331, 2.34108, 2.44463, 2.55436, 2.67068, 2.79403,
     &           2.92489, 3.06377, 3.21123, 3.36785, 3.53427, 3.71116/
C
      data tab2 /0.04135, 0.22059, 0.50078, 0.80624, 1.08109, 1.29772,
     &           1.44991, 1.54374, 1.59063, 1.60296, 1.59178, 1.56602,
     &           1.53227, 1.49518, 1.45778, 1.42190, 1.38856, 1.35821,
     &           1.33093, 1.30662, 1.28505, 1.26593, 1.24899, 1.23394,
     &           1.22054, 1.20857, 1.19782, 1.18813, 1.17936, 1.17139,
     &           1.16412, 1.15745, 1.15132, 1.14566, 1.14043, 1.13556,
     &           1.13103, 1.12679, 1.12283, 1.11912, 1.11562, 1.11233,
     &           1.10922, 1.10629, 1.10350, 1.10087, 1.09836, 1.09598,
     &           1.09371, 1.09155, 1.08948, 1.08751, 1.08563, 1.08382,
     &           1.08209, 1.08043, 1.07883, 1.07730, 1.07583, 1.07441,
     &           1.07304, 1.07173, 1.07046, 1.06923, 1.06805, 1.06691,
     &           1.06580, 1.06473, 1.06370, 1.06270, 1.06173, 1.06079,
     &           1.05987, 1.05899, 1.05813, 1.05729, 1.05648, 1.05569,
     &           1.05493, 1.05418/
      x = xl
      do 10 k = 1,2
         if (x.lt.6.5) then
            i       = 10.0*x
            temp    = tab1(i+1)+(tab1(i+2)-tab1(i+1))*(10.0*x-dble(i))
            func(k) = x*temp
         else
            if (x.lt.40.0) then
               i    = 2.0*x
               temp = tab2(i)+(tab2(i+1)-tab2(i))*(2.0*x-dble(i))
            else
               temp = 1.0+(2.0+6.7/x)/x
            endif
            temp    = temp*exp(x)/x**2
            func(k) = temp
         endif
         x=xh
 10   continue
      fint=func(2)-func(1)-1.0/xh+1.0/xl+log(xh/xl)
      return
      end
C======================================================================
C======================================================================
      function fmax(e)
      implicit real*8 ( a-h, o-z )
      include 'gamlib1.inc' !FURUTA20201007
      include 'gamlib2.inc' !FURUTA20201007
C
      fact = 1.0
      ep   = e-(e-eminxx)*1.0e-06
      call conprb(e,e,ia,isp0)
      t1   = psjsum
      call conprb(e,ep,ia,isp0)
      t2   = psjsum
      if (t2.le.t1) fact = t1
      finv = 1.0/fact
      fmax = fact
      return
      end
C======================================================================
C======================================================================
C
CTAT/MOD BECAUSE OF SAME NAME ------------------------------------(FROM)
      subroutine gamcal(at,zt,e,icalc0,ecv,wt)
CTAT/MOD BECAUSE OF SAME NAME ------------------------------------( TO )
      use levdat
C
      implicit real*8 ( a-h, o-z )
      include 'gamlib1.inc' !FURUTA20201007
      include 'gamlib2.inc' !FURUTA20201007
C
      dimension ip(150), ij(150)
      data inn /150/,
     &     izz /100/
      icalc = icalc0
      ecv   = 0.
      nph   = 0
      iz    = zt
      ia    = at
      iz    = min(izz,iz)
      in    = ia-iz
      in    = min(inn,in)
      ia    = iz+in
      iza   = 1000*iz+ia
      in    = ia-iz
      tp    = 0.
      tcon  = 1.0e-11
cABE initialize nlo @2015/03/17
      nlo = 0
      it0   = 0
      itop  = 0
      ncont = 1
      k1    = iref(iz)
      eminxx  = 0.
      k2    = iref(iz+1)-1
      do 10 k=k1,k2
         if (iza.eq.izao(k)) then
            ik  = k
            nup = ipt(ik+1)-1
            nlo = ipt(ik)
            it0 = nlo
            if (e.le.elo(nlo)) then
               if (it0.gt.0.and.it0.le.len1) wtfin(it0) = wtfin(it0)+wt
               if (it0.le.0) wtnone = wtnone+wt
               return
            endif
            nleva = nup-nlo+1
            if (icalc.le.0) then
               if (nleva.gt.1) eminxx = elo(nlo+1)
               if (e.le.0.) then
                  if (it0.gt.0.and.it0.le.len1) wtfin(it0)=wtfin(it0)+wt
                  if (it0.le.0) wtnone = wtnone+wt
                  return
               endif
               ik = 0
            else
               eminxx = elo(nup)
               itop = nup
               do 50 i=nlo,nup
                  it0 = i
                  if (e.le.elo(i)) then
                     p1 = (e-elo(it0-1))/(elo(it0)-elo(it0-1))
                     if (unirn(dummy).gt.p1) it0 = it0-1
                     e    = elo(it0)
                     itop = it0-1
                     if (it0.ne.nup) then
                        ncont = 0
                        if (icalc.lt.3) then
                           iup = itop+1
                           do 115 ii=nlo,iup
                              j     = ii-nlo+1
                              ip(j) = isign(1,isp(ii))
                              ij(j) = iabs(isp(ii))
                              if (ij(j).eq.0) ip(j)=0
                              if (j.ne.1) then
                                 if (ip(j-1).eq.0) then
                                    ip(j) = 0
                                    ij(j) = 0
                                 endif
                              endif
  115                      continue
                           goto 120
                        endif
                     endif
                     go to 70
                  endif
   50          continue
               it0 = nup
            endif
            go to 70
         endif
   10 continue
      if (e.le.0.) then
         if (it0.gt.0.and.it0.le.len1) wtfin(it0) = wtfin(it0)+wt
         if (it0.le.0) wtnone = wtnone+wt
         return
      endif
      ik = 0
   70 continue

      is = 1

      if ((iz.ge.54  .and. iz.le.78)   .or.
     &    (iz.ge.86  .and. iz.le.122)  .or.
     &    (in.ge.86  .and. in.le.122)  .or.
     &    (in.ge.130 .and. in.le.182))  is = 2

      call lev(iz,ia,is)

      if (ncont.ne.0) eminxx = dmax1(eminxx,de)
      if (e.le.eminxx) ncont = 0
      if (icalc.ge.2) then
         if ((icalc.gt.2) .and. (ncont.ne.0)) then
            ipar0 = 1
            if (unirn(dummy).gt.0.5) ipar0 = -1
            isp0  = mspin(ia,e)
         endif
         iup = itop+1
         do 110 i=nlo,iup
            j     = i-nlo+1
            ip(j) = isign(1,isp(i))
            ij(j) = iabs(isp(i))
            if (ij(j).eq.0) ip(j)=0
            if ((icalc.lt.3) .and. (j.ne.1)) then
               if (ip(j-1).eq.0) then
                  ip(j) = 0
                  ij(j) = 0
               endif
            elseif ((icalc.ge.3) .and. (ip(j).eq.0)) then
               ip(j) = 1
               if (unirn(dummy).gt.0.5) ip(j) = -ip(j)
               ep    = elo(i)
               ij(j) = mspin(ia,ep)
            endif
  110    continue
      endif
C     RE-ENTRY POINT FOR MORE CASCADE
  120 continue
      nph = nph+1
      if (nph.ge.maxph) then
         nph = nph-1
         if (it0.gt.0.and.it0.le.len1) wtfin(it0) = wtfin(it0)+wt
         if (it0.le.0) wtnone = wtnone+wt
         return
      endif
      if (ik.gt.0) then
         if (it0.ne.0.and.it0.eq.nlo) then
            nph = nph-1
            if (it0.gt.0.and.it0.le.len1) wtfin(it0) = wtfin(it0)+wt
            if (it0.le.0) wtnone = wtnone+wt
            return
         endif
         if (ncont.eq.0) then
            tcon = tlev(it0)*1.0e+09
 140        klo  = idk(it0)
            khi  = idk(it0+1)-1
            ndk  = khi-klo+1
            if (ndk.le.0) go to 190
            test = unirn(dummy)
            do 150 kd=klo,khi
               test = test-dmod(pmode(kd),1.d+0)
               if (test.le.0.0) go to 160
  150       continue
            kd=khi
  160       continue
            mode = int(pmode(kd))
            if ((mode.eq.19) .or. (mode.eq.7 .and. icalc.le.3) .or.
     &         (mode.le.6.and.tlev(it0).le.1.0e-12)) go to 190
            if (mode.eq.7 .and. icalc.gt.3) then
               iflg = 3
               it0  = nlo+lfin(kd)
               e1   = elo(it0)
               if (it0.eq.nlo) then
                  e1           = 0.0
                  ncount(iflg) = ncount(iflg)+1
                  xmit1        = unirn(dummy)
                  tp           = tp-tcon*log(xmit1)
                  eph(nph)     = e-e1
                  kfejec(nph)  = 22 ! photon
                  e            = e1
                  if (it0.gt.0.and.it0.le.len1) wtfin(it0)=wtfin(it0)+wt
                  if (it0.le.0) wtnone = wtnone+wt
                  return
               endif
               itop         = it0-1
               ncount(iflg) = ncount(iflg)+1
                  xmit2        = unirn(dummy)
                  tp           = tp-tcon*log(xmit2)
               eph(nph)     = e-e1
               kfejec(nph)  = 22 ! photon
               e            = e1
               go to 120
            endif
            if (mode.eq.17) then
               it0       = nlo+lfin(kd)
               e1        = elo(it0)
               ecv       = ecv+(e-e1)
               ncount(4) = ncount(4)+1
               e         = e1
               if (it0.eq.nlo) then
                  nph = nph-1
                  if (it0.gt.0.and.it0.le.len1) wtfin(it0)=wtfin(it0)+wt
                  if (it0.le.0) wtnone = wtnone+wt
                  return
               endif
               itop = it0-1
               goto 140
            endif
            nph = nph-1
            if (it0.gt.0.and.it0.le.len1) wtfin(it0) = wtfin(it0)+wt
            if (it0.le.0) wtnone = wtnone+wt
            return
         endif
         if (icalc.ge.3) go to 190
         psum = 0.0
         lvls = 0
         do  40 j=nlo,itop
            lvls     = lvls+1
            psum     = psum+(e-elo(j))**3
            ps(lvls) = psum
   40    continue
      else
         psum  = e**3
         ps(1) = psum
         lvls  = 1
      endif
      if (ncont.ne.0) then
         temp = pgs(e,eminxx)
         if (icalc.ge.3) then
            temp = temp*(0.5+10.0*(c2/cra(ia))**2)
            temp = temp*fmax(e)
         endif
         if (temp.le.0.0) then
            ncont = 0
         else
            psum     = psum+temp
            lvls     = lvls+1
            ps(lvls) = psum
         endif
      endif
      goto 280
  190 continue
      if (icalc.ne.1) then
         if (ncont.ne.1) then
            jt0 = it0-nlo+1
            if (ip(jt0).ne.0) then
               ipar0 = ip(jt0)
               isp0  = ij(jt0)
            endif
         endif
         psum = 0.0
         lvls = 0
         if ( ((ncont.eq.1) .and. (icalc.le.2)) .or.
     &        ((ncont.ne.1) .and. (ip(jt0).eq.0)) ) then
            do 245 j=nlo,itop
               lvls     = lvls+1
               psum     = psum+(e-elo(j))**3
               ps(lvls) = psum
 245        continue
            if (ncont.ne.0) then
               temp = pgs(e,eminxx)
               if (icalc.ge.3) then
                  temp = temp*(0.5+10.0*(c2/cra(ia))**2)
                  temp = temp*fmax(e)
               endif
               if (temp.le.0.0) then
                  ncont = 0
               else
                  psum     =psum+temp
                  lvls     =lvls+1
                  ps(lvls) = psum
               endif
            endif
         elseif ( ((ncont.eq.1) .and. (icalc.gt.2)) .or.
     &            ((ncont.ne.1) .and. (ip(jt0).ne.0)) ) then
            do 220 j=nlo,itop
               i     = j-nlo+1
               isp1  = ij(i)
               ipar1 = ip(i)
               j0    = max0(1,iabs(isp1-isp0)/2)
               iparx = ipar0*ipar1
               if (mod(j0,2).ne.0) iparx = -iparx
               lvls  = lvls+1
               egam0  = e-elo(j)
               temp  = trns(egam0,j0,iparx,ia)
               if ((isp1+isp0).eq.2) temp = 1.0e-10*temp
               psum     = psum+temp
               ps(lvls) = psum
 220        continue
            if ((icalc.ge.3) .and. (ncont.ne.0)) then
               temp = pgs(e,eminxx)
               temp = temp*(0.5+10.0*(c2/cra(ia))**2)
               temp = temp*fmax(e)
               if (temp.le.0.0) then
                  ncont = 0
               else
                  psum     = psum+temp
                  lvls     = lvls+1
                  ps(lvls) = psum
               endif
            endif
         endif
      else
         psum = 0.0
         lvls = 0
         do 240 j=nlo,itop
            lvls     = lvls+1
            psum     = psum+(e-elo(j))**3
            ps(lvls) = psum
  240    continue
         if (ncont.ne.0) then
            temp = pgs(e,eminxx)
            if (icalc.ge.3) then
               temp = temp*(0.5+10.0*(c2/cra(ia))**2)
               temp = temp*fmax(e)
            endif
            if (temp.le.0.0) then
               ncont = 0
            else
               psum     = psum+temp
               lvls     = lvls+1
               ps(lvls) = psum
            endif
         endif
      endif
 280  continue
      nkk   = 0
      addon = 1.0
      pmax  = 0.0
      irep  = 0
      ratio = ps(max(1,lvls-1))/ps(lvls)
      if (ratio.le.1.0e-05) irep = 1
  290 continue
      ptest = ps(lvls)*unirn(dummy)
      linc = 1
      l1 = 1
      do 291 isrfgt=1,lvls
         if (ps(l1).gt.ptest) goto 292
         l1 = l1+linc
 291  continue
 292  ll = min(lvls,isrfgt)
      iflg  = 3
      if (ncont.eq.1) iflg = 2
      if (ll.eq.1) then
         it0          = nlo
         e1           = 0.0
         ncount(iflg) = ncount(iflg)+1
                 xmit3        = unirn(dummy)
                 tp           = tp-tcon*log(xmit3)
         eph(nph)     = e-e1
         kfejec(nph)  = 22 ! photon
         e            = e1
         if (it0.gt.0.and.it0.le.len1) wtfin(it0) = wtfin(it0)+wt
         if (it0.le.0) wtnone = wtnone+wt
         return
      endif
      if (ncont.ne.0.and.ll.eq.lvls) then
         iflg = 1
         e1   = skeme(e,eminxx)
         if (e1.lt.eminxx) then
            call exit (2)
         endif
         if (icalc.ge.3) then
            call rejgam(e)
            if (icc.eq.0) then
               if (mod(nkk,100).eq.0) then
               endif
               go to 290
            endif
         endif
      else
         it0   = nlo+ll-1
         e1    = elo(it0)
         itop  = it0-1
         ncont = 0
      endif
      ncount(iflg) = ncount(iflg)+1
                   xmit4        = unirn(dummy)
                   tp           = tp-tcon*log(xmit4)
      eph(nph)     = e-e1
      kfejec(nph)  = 22 ! photon
      e            = e1
      go to 120
C
  390 format (' PMAX=',1p,e13.4,' FINV=',1p,e13.4,' IZ=',i3,' IA=',i4,
     &' RATIO=',1p,e14.4)
  400 format ('0',1p,2e13.4,3i5)
  410 format (' ERROR 2')
      end
C======================================================================
C======================================================================
      function geta(izx,in,is,mode)
      implicit real*8 ( a-h, o-z )
      include 'gamlib1.inc' !FURUTA20201007
      include 'gamlib2.inc' !FURUTA20201007
C
      data b1 /25./,
     &     y1 /1.5/
C
      if (mode.lt.0) then
         call exit (440)
      elseif (mode.eq.0) then
         if (in.ge.9.and.izx.ge.9) then
            st=sz(izx)+sn(in)
            am=in+izx
            fact=(9.17e-3*st+con(is))
            amin=(1.0+y1*(dble(in-izx)/am)**2)/b1
            fact = dmax1(fact,amin)
            geta=fact*am
         else
            geta=0.125*(in+izx)

         endif
      else
         iax=in+izx
         if (iax.gt.240) then
            geta=0.125*iax
         else
            geta=amean(iax)
         endif
      endif
      return
      end
C======================================================================
C======================================================================
      subroutine lev(izt,iax,is)
      implicit real*8 ( a-h, o-z )
      include 'gamlib1.inc' !FURUTA20201007
      include 'gamlib2.inc' !FURUTA20201007
C
C     CALCULATES LEVEL DENSITY MODEL PARAMETERS
C
C         IZT - NUCLEUS CHARGE NUMBER
C         IAX  - NUCLEUS MASS NUMBER
C         IS - SHAPE CONTROL PARAMETER
C              1 - SPHERICAL    2 - DEFORMED
C                1, GILBERT AND CAMERON
C                2,SOME G+C PARAMETERS ARE MODIFIED
C                3,FERMI GAS LEVEL DENSITY ONLY, WITH
C                  PARAMETER MODIFICATION ALLOWED
C
      data inn /150/,
     &     izz /100/
C
C     CALCULATE PAIRING ENERGY AND PARAMETER 'A'
C
      izx    = min0(izt,izz)
      in    = iax-izx
      in    = min0(in,inn)
      am    = iax
      de    = pz(izx)+pn(in)
      mode  = 0
      ax    = geta(izx,in,is,mode)
C
C     CALCULATE CONSTANT IN SIGMA CORRELATION
C
      am3   = am**.33333333
      c1x   = 1.58193/am3
C
C     CALCULATE ENERGY BOUNDARY BETWEEN TWO LEVEL DENSITY FUNCTIONS
C
      ux    = 2.5+150./am
      exx   = ux+de
C
C     CALCULATE LOW ENERGY TEMPERATURE
C
      tn    = sqrt(ax/ux)-1.5/ux
      tn    = 1.0/tn
      c2ssq = 2.0*sqrt(ax*ux)
      p2ux  = c1x*ax*exp(c2ssq)/c2ssq**3
      return
C
   10 format (//5x,'Z = ',i6,' IS NOT CONTAINED WITHIN INTERNAL LEVEL
     &  density tables.'/)
   20 format (/5x,'NEUTRON NUMBER = ',i6,' IS NOT CONTAINED WITHIN
     &  internal level density tables.'//)
      end
C======================================================================
C======================================================================
      integer function mspin(iax,e)
      implicit real*8 ( a-h, o-z )
      include 'gamlib1.inc' !FURUTA20201007
      include 'gamlib2.inc' !FURUTA20201007
      dimension fms(lms), arg(lms)
      dimension aux(lms)
C
      psjsum  = 0.0
      msp0    = mod(iax,2)-1
      u       = e-de
      if (u.le.0.0) then
         psjsum = 0.0
         i      = 1
         msp    = msp0+2*i
         mspin  = msp
         return
      endif
      conx = 0.7104*sqrt(ax*u)*cra(iax)**2
      con1 = log(4.0/conx)
      fm   = msp0
      nup  = 0.5*(sqrt(conx*(con1+200.0))-fm)+1.0
      nup  = min0(lms,nup)
      do 10 i=1,nup
         fms(i) = fm+2*i
   10 continue
      do 20 i=1,nup
         arg(i) = con1-fms(i)**2/conx
   20 continue
      do 30 i=1,nup
         psj(i) = exp(arg(i))
   30 continue
      do 40 i=1,nup
         aux(i) = psj(i)*fms(i)
   40 continue
      psmax = aux(1)
      do 50 i=2,nup
         psmax =   dmax1(psmax,aux(i))
         if (aux(i).le.0.001*psmax) then
            nup = i
            goto 70
         endif
   50 continue
   70 continue
      psjsum = 0.0
      do 80 i=1,nup
         psjsum = psjsum+aux(i)
         psj(i) = psjsum
   80 continue
      if (nup.lt.lms) then
         do 90 i=nup,lms
            psj(i)=psjsum
   90    continue
      endif
      if (psjsum.le.0.0) then
         psjsum=0.0
         i=1
      else
        pp = psjsum*unirn(dummy)
        iinc = 1
        i1   = 1
        do 100 isrfgt = 1,nup
           if (psj(i1).gt.pp) goto 105
           i1 = i1 + iinc
 100    continue
 105    i = min (nup,isrfgt)
      endif
      msp=msp0+2*i
      mspin=msp
      return
      end
C======================================================================
C======================================================================
      function pgs(e,eminx)
      implicit real*8 ( a-h, o-z )
      include 'gamlib1.inc' !FURUTA20201007
      include 'gamlib2.inc' !FURUTA20201007
C
C
      f(x) = 6.+x*(6.+x*(3.+x))
C
      p2   = 0.0
      p1   = 0.0
      smax = (e-eminx)/tn
      if (smax.gt.0.0) then
         s0 = (e-  min(e,exx))/tn
         if (s0.gt.0.0) then
            u     = e-de
            x1    = 2.0*sqrt(u*ax)
            x0    = c2ssq
            x12   = x1**2
            x02   = x0**2
            t1    = 3.0*x12*(x12-x02) + x02**2 + x0*(6.*x12-4.*x02) +
     &              (12.*x02-6.*x12) - 24.*(x0-1.)
            t2    = x12*(x12+2.*x1+6.)-24.*(x1-1)
            t2    = t2*exp(x1-x0)
            t3    = x12**3*exp(-x0)*fint(x0,x1)
            p2    = t3-t2+t1
            p2    = 2.*p2*p2ux*x0*x02/(4.0*ax)**4
            temp2 = f(s0)-f(smax)*exp(s0-smax)
            p1    = p2ux*temp2*tn**4
         else
            temp2 = f(s0)-f(smax)*exp(s0-smax)
            p1    = (p2ux*temp2*tn**4)*exp((e-exx)/tn)
         endif
      endif
      pgs = p1+p2
      return
      end
C======================================================================
C======================================================================
      subroutine rejgam(e)
      implicit real*8 ( a-h, o-z )
      include 'gamlib1.inc' !FURUTA20201007
      include 'gamlib2.inc' !FURUTA20201007
   10 continue
      call conprb(e,e1,ia,isp0)
      prob  = psjsum*finv
      pmax  =  dmax1(pmax,psjsum)
      ptemp = addon*prob
      pp    = unirn(dummy)
      if (pp.lt.ptemp) then
         pp    = psjsum*unirn(dummy)
         iinc = 1
         i1   = 1
         do 40 isrfgt = 1,nstat
            if (psj(i1).gt.pp) goto 45
            i1 = i1 + iinc
 40      continue
 45      i = min (nstat,isrfgt)
         isp0 = mod(ia,2)-1+2*i
         temp = 10.0*(c2/cra(ia))**2
         if (unirn(dummy)*(1.0+temp).gt.temp) ipar0=-ipar0
         icc=1
         nacc=nacc+1
         return
      endif
      icc  = 0
      nrej = nrej+1
      nkk  = nkk+1
      if (irep.eq.0) return
      if (mod(nkk,100).eq.0) then
      endif
      e1 = skeme(e,eminxx)
      if (e1.ge.eminxx) go to 10
      call exit  (2)
      return
C
   80 format (' PMAX=',1p,e13.4,' FINV=',1p,e13.4,' IZ=',i3,' IA=',i4,
     1' RATIO=',1p,e14.4)
   90 format ('0',1p,2e13.4,3i5)
  100 format (' ERROR 2')
  110 format (' PTEMP=',1p,e13.4,' ADDON=',1pe13.4,' E=',1p,e13.4,2i5)
      end
C======================================================================
C======================================================================
      function skeme(e,eminx)
      implicit real*8 ( a-h, o-z )
      include 'gamlib1.inc' !FURUTA20201007
      include 'gamlib2.inc' !FURUTA20201007
      smax    = (e-eminx)/tn
   10 continue
      if (smax.lt.3.0) then
   20    continue
         x = unirn(dummy)
         r = exp(smax*(1-x))*x**3
         if (unirn(dummy).gt.r) go to 20
         s = smax*x
      else
   30    continue
         r1 = unirn(dummy)
         r2 = unirn(dummy)
         r3 = unirn(dummy)
         r4 = unirn(dummy)
         s  = -log(r1*r2*r3*r4)
         if (s.gt.smax) go to 30
      endif
      e1x = e-s*tn
      if (e1x.gt.exx) then
         u   = e1x-de
         rej = (ux/u)**3*exp(4.0*sqrt(ax*u)-6.0-2.0*(u+ux)/tn)
      if (rej.lt.1.0e-6) rej=1.0e-6
         r1  = unirn(dummy)
         if (r1**2.gt.rej) go to 10
      endif
      skeme = e1x
      return
      end
C======================================================================
C======================================================================
      function trns(e,j,ipar,iax)
      implicit real*8 ( a-h, o-z )
      include 'gamlib1.inc' !FURUTA20201007
      include 'gamlib2.inc' !FURUTA20201007
C
C
      temp = e**3
      x    = c1*cra(iax)*e
      if (j.gt.1) temp = fj(j)*temp*x**(j+j-2)
      fpar = 0.5*(1+ipar)
      q    = c2/cra(iax)
      temp = temp*(fpar+10.0*q**2+fj(j+1)*x**2/fj(j))
      trns = temp
      return
      end
C======================================================================
************************************************************************
*                                                                      *
      subroutine therdec(iares, izres, eexe, nsta, totJ, nparit, iflg,
     & sumwid)
*                                                                      *
*      This routine defines theory-based deexcitation                  *
*          (destination level and gamma energy) of spallation products *
*      from continuous levels to continuous/discrete levels            *
*                                                                      *
*     input  :                                                         *
*       iares  : mass number                                           *
*       izres  : atomic numbe r                                        *
*       eexe   : excitation energy  [MeV]                              *
*       totJ   : angular momentum                                      *
*       nparit : parity                                                *
*       iflg   : flag (0:called from dexgam2, 1: called from gammag)   *
*       nph    : number of emitted gamma rays so far                   *
*                                                                      *
*     output :                                                         *
*       nsta   : current level                                         *
*       eexe   : excitation energy  [MeV]                              *
*       sumwid : sum of decay width [MeV]                              *
************************************************************************

      use levdat
      use levdenmod

      implicit double precision(a-h, o-z)

      double precision, allocatable :: wsum1(:,:)
      double precision, allocatable :: wsum2(:)

      if(totJ .eq. -1.d0) totJ = totJ_init(iares, izres, totJ)

*---------------------------------------------------------------
*------De-exitation ------------------------------------

      nlevel = ubound(elevel, 1) ! number of levels
* Calculate transition probability depending on lambda and parity.

      wstot = 0.d0
      eshif = dmin1(elevel(nlevel), 3.d0) ! Boundary between continuum level and discrete level region. The highest of ENSDF or 3000keV.

      IF(elevel(1) .gt. 3.d0) eshif = elevel(1)

      IF(elevel(nlevel) .le. eshif) then ! # of levels below eshif, which should be treated as discrete levels
       levl3m = nlevel
      ELSE
       do icont = 1, nlevel   ! levels below 3000 keV
        If(elevel(icont) .gt. eshif) then
        levl3m = icont - 1
        exit
        ENDIF
       enddo
      ENDIF

      levcan1 = idnint( dint( (eexe - eshif) /1.d-2) )  ! number of level interval above 3MeV   subject to be weighted with level density
      levcan2 = levl3m                         ! number of levels below 3MeV           no weighting since each corresponds to a certain level

      allocate (wsum1(levcan1,5))
      wsum1 = 0.d0

* Transition to continuum level above 3MeV
      do icont = 1, levcan1         ! Loop for destination level

* summation over the decay channels (different lambda, E mode/ M Mode)        ! jcont = 1 (dipole), 2(quadrupole) are considered
       do jcont = max0(idnint(3.d0 - totJ),1), 5            ! Loop for desination level spin (spinow - 2) ` (spinow + 2)

       do kcont = max0(1, iabs(jcont - 3)), idnint(totJ)*2 + jcont - 3    ! kcont = lambda parameter (multipolarity)

        IF( (-1)**kcont * nparit .ge. 0) then ! electric mode

         wsum1(icont,jcont) = wsum1(icont,jcont) +
     &    4.4d0 * dble(kcont + 1)
     &    /dble(kcont) * (3.d0 / dble(3 + kcont) )**2
     &    * canfact1(eexe, icont, eshif, kcont)
     &    * canfact2(iares, kcont)
     &    * (2.d0 * (totJ + dble(jcont - 3) ) + 1.d0)
     &    * rho_levden(iares, izres, dble(icont) * 1.d-2 + eshif,                             ! from here Weighting by level density
     &     totJ + dble(jcont - 3), 2) * 1.d-2   ! 0.01 MeV bin

        ELSE                                                   !  magnetic mode

         wsum1(icont,jcont) = wsum1(icont,jcont) +
     &    1.9d0 * dble(kcont + 1)
     &    /dble(kcont) * (3.d0 / dble(3 + kcont) )**2
     &    * (2.d0 * (totJ + dble(jcont - 3) ) + 1.d0)
     &    * rho_levden(iares, izres, dble(icont) * 1.d-2 + eshif,                             ! from here Weighting by level density
     &     totJ + dble(jcont - 3), 2) * 1.d-2   ! 0.01 MeV bin
     &    * canfact1(eexe, icont, eshif, kcont)
     &    * canfact2(iares, kcont)
     &    * (1.2d0 * dble(iares**(1.d0/3.d0) ) ) ** (-2)

        ENDIF

       enddo

       wstot = wstot + wsum1(icont,jcont)

      enddo
      enddo

* Transition to discrete levels below 3 MeV

      allocate(wsum2(0:levcan2))
      wsum2 = 0.d0

      do icont = 0, levcan2         ! Loop for destination level.

       if( elevel(icont) .gt. eexe ) exit

* summation over the decay channels (different lambda, E mode/ M Mode)
       IF((ndch(icont) .eq. 0 .and. elevel(icont) .ge. 1.5d0) .or.
     &  spinpar(icont,1) .eq. -1.d0 ) cycle ! disregard this level if decay information of the destination level is not given or spin is unknown

       do kcont = max0(1, idnint(abs(spinpar(icont,1) - totJ))),
     &  idnint(spinpar(icont,1) + totJ)     ! kcont = lambda (multipolarity)

        IF( (-1)**kcont * parit(icont) * nparit .ge. 0) then !  electric mode

         wsum2(icont) = wsum2(icont) +
     &    4.4d0 * dble(kcont + 1)
     &    /dble(kcont) * (3.d0 / dble(3 + kcont) )**2
     &    * (2.d0 * (totJ + dble(kcont) ) + 1.d0)
     &    * canfact3(eexe, elevel(icont), kcont, iares)
        ELSE                                                   !  magnetic mode

         wsum2(icont) = wsum2(icont) +
     &    1.9d0 * dble(kcont + 1)
     &    /dble(kcont) * (3.d0 / dble(3 + kcont) )**2
     &    * (2.d0 * (totJ + dble(kcont) ) + 1.d0)
     &    * canfact3(eexe, elevel(icont), kcont, iares)
     &    * (1.2d0 * dble(iares**(1.d0/3.d0) ) )**(-2)
        ENDIF

       enddo
       wstot = wstot + wsum2(icont)

      enddo

      if(iflg .eq. 1) then ! gamma-neutron competition.
       sumwid = wstot
       return
      endif

* Randomly sample one of the transitions

      IF(wstot .eq. 0.d0) then ! Workaround to avoid infinite loop. Find level
         do icont = levcan2, 1, -1
           if(ndch(icont) .ne. 0. and. eexe .gt. elevel(icont)) exit
         enddo
         nsta = icont
         nph = nph + 1
         eph(nph) = eexe - elevel(nsta)
         kfejec(nph) = 22 ! photon
         eexe = elevel(nsta)
         totJ = spinpar(nsta,1)
       deallocate (wsum1)
       return
      ENDIF

      selran = unirn(dummy)
      prob = 0.d0

      do icont = 1, levcan1   ! sample a destination level from continuous or discrete levels
       do jcont = max0(idnint(3.d0 - totJ),1), 5

       prob = wsum1(icont,jcont)/wstot + prob
       IF(prob .ge. selran) then
        nph = nph + 1
        eph(nph) = eexe - ((dble(icont) - unirn(dummy)) * 1.d-2
     &   + eshif)  ! Not to see 10 keV step structure, randomly sample energy within 10 keV width
        kfejec(nph) = 22 ! photon
        eexe = eexe - eph(nph)
        totJ = totJ + (dble(jcont) -3.d0)
        goto 120
       ENDIF
       enddo
      enddo

      do icont = 0, levcan2
       prob = wsum2(icont)/wstot + prob

        IF(prob .ge. selran) then
         nsta = icont
         nph = nph + 1
         eph(nph) = eexe - elevel(icont)
         kfejec(nph) = 22 ! photon
         eexe = eexe - eph(nph)
        exit
        ENDIF
      enddo

 120   deallocate (wsum1)

*---------------------------------------------------------------

      return
      end subroutine
************************************************************************


************************************************************************
*                                                                      *
      subroutine firdec(iares,izres,eexe,nsta,totJ)
*                                                                      *
*      This routine defines first deecxitation                         *
*          (destination level and gamma energy)                        *
*      this subroutine is for spallation. n-capture products are sent  *
*      to subroutine neufdc                                            *
*                                                                      *
*     input  :                                                         *
*       iares  : mass number                                           *
*       izres  : atomic numbe r                                        *
*       eexe   : excitation energy  [MeV]                              *
*                                                                      *
*     output :                                                         *
*       nsta   : current level                                         *
*       eexe   : excitation energy  [MeV]                              *
*       totJ   : total angular momentum [hbar]                         *
*                                                                      *
************************************************************************

      use levdat

      implicit double precision(a-h, o-z)

*-----------------------------------------------------------------------

       totJ = totJ_init(iares, izres, totJ) ! determine angular momentum

*---------------------------------------------------------------
*------De-excite to the closest level --------------------------
      nlevel = ubound(elevel,1) ! Levels are 0 - ubound(elevel)

      If(eexe .gt. elevel(nlevel) + 1.d-3) then
* In case the energy is above evaluated levels and totJ changes less than 2
       if(abs(totJ - spinpar(nlevel,1)) .le. 2.d0 ) then
        nsta = nlevel
        nph = 1
        eph(1) = eexe - elevel(nlevel)
        kfejec(nph) = 22 ! photon
        eexe = elevel(nlevel)
        totJ = spinpar(nlevel,1)
       endif
      ELSE
* In case the energy is in between evaluated levels
       do icont = nlevel, 1, -1

        If( abs(eexe - elevel(icont)) .le. 1.d-3 .and.
     &      abs(totJ - spinpar(icont,1)) .le. 2.d0) then ! 2015/08/18 ogawa, In case the energy is at evaluated levels (neutron inelastic reaction)
         nph  = 0
         nsta = icont
         eexe = elevel(nsta)
         totJ = spinpar(nsta,1)
         exit
        elseIf( elevel(icont) .le. eexe  .and.
     &          abs(totJ - spinpar(icont-1,1)) .le. 2.d0) then
         if(ndch(icont-1) .gt. 0) then
          nsta   = icont - 1
          nph    = 1
          eph(1) = eexe - elevel( icont - 1 )
          kfejec(nph) = 22 ! photon
          eexe = elevel ( nsta )
          totJ = spinpar(nsta,1)
         endif
         exit
        ENDIF
       end do
      ENDIF

      If(nsta .eq. -1) then ! If no levels hit, eliminate angular momentum condition
       do icont = nlevel, 1, -1
        If( abs(eexe - elevel(icont)) .le. 1.d-3) then ! 2015/08/18 ogawa, In case the energy is at evaluated levels (neutron inelastic reaction)
         nph  = 0
         nsta = icont
         eexe = elevel(nsta)
         totJ = spinpar(nsta,1)
         exit
        elseIf( elevel(icont) .le. eexe) then
          nsta   = icont - 1
          nph    = 1
          eph(1) = eexe - elevel( icont - 1 )
          kfejec(nph) = 22 ! photon
          eexe = elevel ( nsta )
          totJ = spinpar(nsta,1)
         exit
        ENDIF
       end do
      endif


      return
      end subroutine
************************************************************************


************************************************************************
*                                                                      *
      function totJ_init(iares, izres, totJ)
*                                                                      *
*        Input :                                                       *
*         totJ   : spin of residual nucleus before                     *
*         iares  : mass number of residual nucleus                     *
*                                                                      *
*        Output :                                                      *
*         totJ_init  : angular momentum after evaporation              *
*                                                                      *
*      Estimate the angular momentum before de-excitation              *
*          the total angular momentum is predicted based on the        *
*          angular momentum of incident particle and target nucleus    *
*          as wel as orbital angular momentum                          *
*          Delete this function when evaporation model can assess J    *
************************************************************************
      use NGSDATAMOD, only : spin
      implicit double precision(a-h, o-z)

      include 'param00.inc'

* Target material
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

* Reaction information
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)

* Information on outgoing particles in dynamic phase
      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /cascid/ jcasc
!$OMP THREADPRIVATE(/cascid/)

*---------------------------------------------------------------------------------------
*------target, projectile and products spin before and after reactions -----------------

* Angular momentum of fragment (calculated based on kinematics of nucleons inside fragment)
      cfact = 0.d0

      If(jcasc .eq. 5 ) then        ! incl
       cfact = 0.6d0
      ELSEIf(jcasc .eq. 4 ) then    ! jqmd
       cfact = 0.8d0
      ELSE                          ! The other reaction models. EBITEM is not adjusted for the other models
       cfact = 1.d0
      ENDIF

      if(totJ .eq. -1.d0) then
        totJ_init = 0.d0
      else
        totJ_init = dble(nint ( totJ * cfact ))
      endif

      if(cfact .eq. 0.d0) then ! not incl or jqmd
* Orbital angular momentum of projectile-target combination @hbc in (GeV fm)
       ospp = 0.d0
* projectile intrinsic spin
       projsp = spin(ityp, ktyp)
* target spin
       spt = spin( 19, mathz * 1000000 + (mathz + mathn) )
       totJ_init = spincomb( spincomb( projsp, ospp), spt) ! spin of compound nucleus.
      endif

* ejectile spin    THIS should be handled by GEM. Correct some day.
      do icont = 2, nclsts
       IF(jclusts(1,icont) .eq. 0 .and. jclusts(2,icont) .eq. 1 ) then
        ejecsp = spincomb( dble( jclusts(0,icont) ), 0.5d0 )
       ELSEIF(jclusts(1,icont) .eq. 1 .and. jclusts(2,icont).eq. 0) then
        ejecsp = spincomb( dble( jclusts(0,icont) ), 0.5d0 )
       ELSEIF(jclusts(1,icont) .eq. 1 .and. jclusts(2,icont).eq. 1) then
        ejecsp = spincomb( dble( jclusts(0,icont) ), 1.d0 )
       ELSEIF(jclusts(1,icont) .eq. 1 .and. jclusts(2,icont).eq. 2) then
        ejecsp = spincomb( dble( jclusts(0,icont) ), 0.5d0 )
       ELSEIF(jclusts(1,icont) .eq. 2 .and. jclusts(2,icont).eq. 1) then
        ejecsp = spincomb( dble( jclusts(0,icont) ), 0.5d0 )
       ELSEIF(jclusts(1,icont) .eq. 2 .and. jclusts(2,icont).eq. 2) then
        ejecsp = dble( jclusts(0,icont) )
       ELSEIF(jclusts(1,icont) .ne. izres .and. jclusts(2,icont) .ne.
     &  iares - izres ) then
        ejecsp = spincomb( dble( jclusts(0,icont) ),
     &    spin(19, jclusts(1,icont) * 1000000
     &           + jclusts(1,icont) + jclusts(2,icont) )  )
       ELSE
        ejecsp = 0.d0
       ENDIF

* Spin before deexcitation
       totJ_init = spincomb( totJ_init, ejecsp) ! spin of residue

      enddo

* Odd/Even check
      if( mod(iares, 2) .eq. 1) then
       totJ_init = dble(2*int(totJ_init)+1)/dble(2)
      else
       totJ_init = dble(idint(totJ_init))
      endif

      return
      end function

************************************************************************
*                                                                      *
      function spincomb( sp1, sp2)
*                                                                      *
*    combine two angular momenta to determine final momentum at random *
*    considering Clebche-Gordan coefficient                            *
************************************************************************

      implicit double precision(a-h, o-z)

      smin = abs(sp1 - sp2)
      smax = sp1 + sp2
      r = unirn(dummy)

      do i = 0, int(min(sp1,sp2)*2.0)
          r = r - ((sp1 + sp2 - i) * 2.d0 +1.d0)
     &            /(sp1*2.d0+1.d0)/(sp2*2.d0+1.d0)
          if(r .le. 0.d0) exit
      enddo

      spincomb = smax - dble(i)

      return
      end function
************************************************************************

************************************************************************
*                                                                      *
      function dfact(iarg)
*                                                                      *
*      Double step factorial                                           *
************************************************************************

      implicit double precision(a-h, o-z)

      dfact = 1.d0
      do
       IF(iarg .le. 0) then
        exit
       ELSE
        dfact = dfact * dble(iarg)
        iarg = iarg - 2
       ENDIF
      enddo

      return
      end function
************************************************************************

************************************************************************
*                                                                      *
      function canfact1(eexe, icont, eshif, kcont)
*                                                                      *
*      Factorial canceling routine                                     *
************************************************************************

      implicit double precision(a-h, o-z)

      canfact1 = (eexe - dble(icont) * 1.d-2 - eshif)/ 1.97327d2
      krem = kcont
      do
       IF(kcont .le. 0) then
        exit
       ELSE
        canfact1 = canfact1
     &    * ((eexe - dble(icont) * 1.d-2 - eshif)/ 1.97327d2)**2
     &    / (2 * kcont + 1)
        kcont = kcont - 1
       ENDIF
      enddo

      kcont = krem
      return
      end function
************************************************************************

************************************************************************
*                                                                      *
      function canfact2(iares, kcont)
*                                                                      *
*      Factorial canceling routine                                     *
************************************************************************

      implicit double precision(a-h, o-z)

      canfact2 = 1.d0
      krem = kcont
      do
       IF(kcont .le. 0) then
        exit
       ELSE
        canfact2 = canfact2
     &      * (1.2d0 * dble(iares**(1.d0/3.d0) ) )**2
     &      / (2 * kcont + 1)
        kcont = kcont - 1
       ENDIF
      enddo

      kcont = krem
      return
      end function
************************************************************************

************************************************************************
*                                                                      *
      function canfact3(eexe, elev, kcont, iares)
*                                                                      *
*      Factorial canceling routine                                     *
************************************************************************

      implicit double precision(a-h, o-z)

      canfact3 = ( eexe - elev )/ 1.97327d2
      krem = kcont
      do
       IF(kcont .le. 0) then
        exit
       ELSE
        canfact3 = canfact3
     &    * (1.2d0 * dble(iares**(1.d0/3.d0) ) )**2
     &    * ( ( eexe - elev )/ 1.97327d2 )**2
     &    / (2 * kcont + 1)**2
        kcont = kcont - 1
       ENDIF
      enddo

      kcont = krem
      return
      end function
************************************************************************
