c***********************************************************************
c***********************************************************************
c                                                                      *
c                                                          Aug. 1999   *
c                                                                      *
c           This is JAM version 1.009   by Yasushi Nara                *
c                                                                      *
c               Jet Aa Microscopic transport model                     *
c                                                                      *
c   Many interesting phyiscs are jammed in heavy ion collisions.       *
c   I hope you have a real jam with JAM.                               *
c   JAM is a Monte Calro program to solve numerically                  *
c   the relativisitc Boltzmann equation for hadrons,                   *
c   their resonance states and on-shell partons                        *
c   from (mini-)jet using Cascade method.                              *
c                                                                      *
c   Resonance production cross section part for NN interaction is      *
c   together with Akira Ohnishi, Naohiko Ohtsuka, Koji Niita           *
c   and Satoshi Chiba.                                                 *
c   Box calculation  part has been written by Toshiki Maruyama.        *
c   Input part has been nicely improuved by T. Maruyama.               *
c   Angular distribution for resonance production and decay written    *
c   by Akira Ohnishi.                                                  *
c   BUU part together with Koji Niita.                                 *
c   Deuteron coalecence by Yuichi Hirata.                              *
c                                                                      *
c                                                                      *
c***********************************************************************
c                                                                      *
c        PART 1: Main steering routines                                *
c                                                                      *
c   List of subprograms in rough order of relevance with main purpose  *
c      (S = subroutine, F = function, B = block data, E = entry)       *
c                                                                      *
c  s  jamevt   to simulate one event                                   *
c  s  jaminit  to initialize the overall run and simulation procedure  *
c  s  jaminie  to initialize jam at each event                         *
c  s  jamfin   to finish up the overall run, to write final information*
c  s  jamread  to read input configuration file                        *
c  s  adjust   to get the first word from a text                       *
c  s  compress to remove blanc characters from a text                  *
c  s  caps     to convert a text into capital letters                  *
c  s  checkcfg to check errors in the config file                      *
c  s  readcfg  to open and read the config file                        *
c  s  iinput   to input integer variable from the config               *
c  s  rinput   to input real variable from the config                  *
c  s  tinput   to input text variable from the config                  *
c  s  jaminbm  to identifies the two incoming particles                *
c  s  jamlogo  to print front page and echo of input                   *
c  s  jamerrm  to inform user of errors in program execution           *
c  s  jamlist  to list event record or particle data                   *
c  s  jamname  to give the particle/parton/nucleus name                *
c  s  jamdisp  to display particles on terminal.                       *
c  s  jamboost to boost the ground state nuclei                        *
c  s  jamgrund to make ground state nucleus.                           *
c  s  jamjeti  to initialize string decay routine phyia                *
c  s  jamhijin to initialize jet cross section HIJING                  *
c  s  jampyin  to initialize PYTHIA for jet production                 *
c  f  jamk     to give charge/strangeness of particle i                *
c  f  jamchge  to gives three times the charge for a particle/parton   *
c  f  jamcomp  to compress standard kf flavour code to internal kc     *
c  s  jamsetpa to setup particle kc code                               *
c  s  jamupdat to update particle data                                 *
c                                                                      *
c                                                                      *
c***********************************************************************
c***********************************************************************
c                                                                      *
c     Verion 0.001  by Yasushi Nara       last revised 1997/10/29      *
c     Verion 0.083e by Yasushi Nara       last revised 1998/03/13      *
c     Verion 0.085a by Yasushi Nara       last revised 1998/04/07      *
c                             Box calculation by Toshiki Maruyama      *
c     Verion 0.088a                       last revised 1998/06/16      *
c                    B-W integral fitted by Chiba is implemented.      *
c     Verion 0.090a  many changes         last revised 1998/07/28      *
c     Verion 1.004b  buu part             last revised 1998/11/12      *
c     Verion 1.005   buu coulomb part     last revised 1998/11/15      *
c                                                                      *
c                                                                      *
c***********************************************************************

      subroutine jamevt(iev,ierr)

c...Purpose: to simulate complete one event.

      include 'jam1.inc'
      include 'jam2.inc'
c...For particle multiplicities.
      common/jamxml/mult(-500:500)
      save /jamxml/
!$OMP THREADPRIVATE(/jamxml/)
c...For random seed
      common/pydatr/mrpy(6),rrpy(100)
      save /pydatr/
!$OMP THREADPRIVATE(/pydatr/)
      common/rseed/iseed
      save /rseed/
!$OMP THREADPRIVATE(/rseed/)

      common /pionflag/ iqdflag, ipionflag, istrflag, iphreturn
!$OMP THREADPRIVATE(/pionflag/)
      common /paraj/  mstz(300), parz(300)

c...Save current event number.
      mstd(21)=iev

      imev=0
 3000 continue
      imev=imev+1
      if(imev.gt.30) call jamerrm(30,0,'(jamevt:) infinit loop???')

c...Initialize JAM parameters, etc for each event.

            call jaminie(mrun,ierr)

            if( ierr .ne. 0 ) return

      if(mrun.eq.1) return

c...Save current random seed.
      mstd(22)=mrpy(1)

c...Global time.
      pard(1)=0.0d0

c...Some initialization for summary.

c...Rest collision counters.
      mstd(29)=0   ! collision counter
      mstd(30)=0   ! number of dead particle
      mstd(41)=0   ! elastic
      mstd(42)=0   ! inelastic
      mstd(43)=0   ! absorption
      mstd(44)=0   ! BB collision
      mstd(45)=0   ! MB collision
      mstd(46)=0   ! MM collision
      mstd(47)=0   ! antiBB collision
      mstd(48)=0   ! parton-hadron collision
      mstd(49)=0   ! parton-parton collision
      mstd(50)=0   ! decay
      mstd(51)=0   ! Pauli block
      mstd(52)=0   ! low energy cut
      mstd(53)=0   ! decay after simul.
      mstd(54)=0   ! Econ block
      mstd(55)=0   ! Number of hard scatt.


c...Option: output phase space data.
      if(mstc(164).eq.1) call jamout1(pard(1))
      if(mstc(164).eq.2) call jamout2(pard(1))

c...Reset time dependent analysis.
      if(mstc(3).eq.1) then
c...Option: output particles on the terminal.
        npri=1
      else
        npri=max(1,int(parc(7)/parc(2)))
      endif

c...Start time evolution.
c======================================================================*
      do 10000 kdt = 1, mstc(3)
c======================================================================*

        mstd(23)=kdt

c...Cascade.
        if(mstc(6).ge.0) then

          if(mstc(3).ge.2.and.mod(kdt-1,npri).eq.0) then
            if(mstc(8).ge.1) call jamdisp(6,kdt)
          endif

            call jamcoll

cABE add @2014/07/29, to retry cal. for photonuclear reaction
            if (iphreturn .eq. 1) return

c...LPC:Collision according to mean free path.
        else if(mstc(6).le.-100) then

            call jamfpath


c...Multipule AA collision by Glauber geometry.
        else
          call jamglaub(igl)
          if(igl.eq.1) goto 3000

c...Option for final state hadron cascade.
          if(mod(abs(mstc(6))/10,10).eq.2) then
            mstc55=mstc(55)
            mstc(55)=1        ! resonance decay forbiden
            call jamfdec      ! fragment strings
            call jamcoll      ! hadron cascade
            mstc(55)=mstc55
          endif
        endif

c stop timing definition. Maximum is 200 fm/c
        if(mstc(190) .eq. 1 .and. kdt .eq. mstz(64)) exit ! 2015/10/16  QMD stop time given by nqtmax
        if(mstc(190) .eq. 0 .and. kdt .eq. mstc(3)) exit ! 2015/10/16  JAM stop time

c======================================================================*
10000 continue  !******************** End time evolution
c======================================================================*

c...1+1 simulation: only inelasitic event is allowed.
      if(mstc(17).eq.1.and.mstd(42).eq.0) then
        call jamerrm(1,0,' Elastic event in 1+1 simulation')
        goto 3000
      endif

c...Final decay of resonances delta, n*, rho ....

      if(mstc(41).eq.1) call jamfdec
      if(mstc(8).ne.0) call jamdisp(6,100)

*-----------------------------------------------------------------------

      if(mstd(30).gt.0) then

            call jamedit

         if( mstc(190) .ne. 0 ) then

            call jam2qmd
            call caldisaR

         end if

      end if

*-----------------------------------------------------------------------

c----------------------------------------------------------------------*
c.....Summary and final output for each run
c----------------------------------------------------------------------*
      pard(1)=pard(1)+parc(2)
      if(mstc(8).ge.2.and.mstc(6).ge.-100)
     $  call jamcheck('After Simul. i.e. Final Check')

c...Output analysis results.

c...Count hadron multiplicites.
      do i=1,nv
        if(k(1,i).gt.10) goto 300
        kf=k(2,i)
        kc=jamcomp(kf)
        if(kc.le.0) goto 300
        kch=jamchge(kf)
        mult(0)=mult(0)+1
        if(kch.ne.0) mult(41)=mult(41)+1
        if(kch.lt.0) mult(42)=mult(42)+1
        if(kch.gt.0) mult(43)=mult(43)+1
        mult(kc*isign(1,kf))=mult(kc*isign(1,kf))+1
 300  end do

c...Option for deuteron coalesence.
      if(mstc(45).ge.1) call jamdeut

c...Update collision counters.
      pard(71)=pard(71)+mstd(41)   ! elastic
      pard(72)=pard(72)+mstd(42)   ! inelastic
      pard(73)=pard(73)+mstd(43)   ! absorption
      pard(78)=pard(78)+mstd(50)   ! decay
      pard(79)=pard(79)+mstd(51)   ! Pauli block
      pard(80)=pard(80)+mstd(52)   ! low energy cut
      pard(81)=pard(81)+mstd(53)   ! decay after simul.
      pard(82)=pard(82)+nv
      pard(83)=pard(83)+nmeson
      pard(84)=pard(84)+mstd(54)   ! Econ block
      pard(85)=pard(85)+mstd(48)   ! h-p coll.
      pard(86)=pard(86)+mstd(49)   ! p-p coll.
      pard(87)=pard(87)+mstd(55)   ! hard scatt.
      pard(74)=pard(74)+mstd(44)   ! BB collision
      pard(75)=pard(75)+mstd(45)   ! MB collision
      pard(76)=pard(76)+mstd(46)   ! MM collision
      pard(77)=pard(77)+mstd(47)   ! antiBB collision

      if(mstc(8).ge.3) then
        ih=mstc(38)
        write(ih,*)'Event=           ',mstd(21)
        write(ih,*)'Total collisions ',(mstd(41)+mstd(42))/mstc(5)
        write(ih,*)'inelastic        ',mstd(42)/mstc(5)
        write(ih,*)'elastic          ',mstd(41)/mstc(5)
        write(ih,*)'absorption       ',mstd(43)/mstc(5)
        write(ih,*)'decays           ',mstd(50)/mstc(5)
        write(ih,*)'h-p coll.        ',mstd(48)/mstc(5)
        write(ih,*)'p-p coll.        ',mstd(49)/mstc(5)
        write(ih,*)'Pauli blocks     ',mstd(51)/mstc(5)
        write(ih,*)'low energy cuts  ',mstd(52)/mstc(5)
        write(ih,*)'Econ block       ',mstd(54)/mstc(5)
        write(ih,*)'decay after simul.',mstd(53)/mstc(5)
        write(ih,*)'Number of hard scatt.',mstd(55)/mstc(5)
        write(ih,*)'Total particle   ',nv/mstc(5)
        write(ih,*)'Total mesons     ',nmeson/mstc(5)
      endif

      if(nv.ne.nbary+nmeson) then
        write(check(1),'(''nv='',i9)')nv
        write(check(2),'(''nbary='',i9)')nbary
        write(check(3),'(''nmeson='',i9)')nmeson
        call jamerrm(11,3,'nv .ne.nbary+nmeson')
      endif

      end

c***********************************************************************

      subroutine jaminit(nev,bmin,bmax,dt,nsp,chfram,chbeam,chtarg,cwin)

c...Initialize the overall run and simulation procedure.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      common/jamxml/mult(-500:500)
!$OMP THREADPRIVATE(/jamxml/)

      character chfram*8,chbeam*8,chtarg*8,cwin*15

c...For random seed
      common/pydatr/mrpy(6),rrpy(100)
      save /pydatr/
!$OMP THREADPRIVATE(/pydatr/)
      common/rseed/iseed
      save /rseed/
!$OMP THREADPRIVATE(/rseed/)

c...Reset number of jaminit call.
      mstc(21)=mstc(21)+1

c...Reset par arrays.
      do i=1,200
       mstd(i)=0
       pard(i)=0.0d0
      end do
      do i=1,50
       mste(i)=0
       pare(i)=0.0d0
      end do
      do i=1,10
        check(i)=' '
      end do

c...Initialize particle multiplicities.
      do i=-500,500
        mult(i)=0
      end do


      if(mstc(21).eq.1) then

c...Setup offset for kc code.
        call jamsetpa

      endif

c...Read simulation parameters from prepared file.
      call jamread(nev,bmin,bmax,dt,nsp,chfram,chbeam,chtarg,cwin)

      if(mstc(21).eq.1) then

c...CPU time and initialzation of ramdom number.


c...Inititalize Multi run control
        if(mstc(39).gt.0) call jammrun(0,mrun)

c...Initialize BUU part.
        if(mstc(5).ge.2) call jambuuin(0)

      endif

c...Inititalize event study and analysis.

c...Inititaze soft and hard interaction part and string fragmentation.

            call jamjeti

c...Min. bias event.
      if(abs(parc(4)).gt.90)
     $   parc(4)=(pard(40)+pard(50))*sign(1d0,parc(4))

c...Impact parameter bin.
      pard(4)=abs(abs(parc(4))-parc(3))/mstc(10)
      if(mstc(10).gt.1.and.pard(4).eq.0.d0) then
          call jamerrm(1,0,'(jaminit:) Warning '
     $     //' You set impact parameter bin, but bmin=bmax')
      endif

c...Number of event for each impact parameter bin.
      mstd(24)=mstc(2)

c...Total number of event is changed to then number of bin times
c...input event number.
      mstc(2)=mstc(2)*mstc(10)

c...Current minimum impact parameter.
      pard(3)=parc(3)

c...Box by maru  define number of box.
      mstd(15)=1
      if(mstc(4).eq.10) mstd(15)=mstc(9)
c...end box

c...Print front page of input-echo ,print simulation conditions,etc.

      end

c***********************************************************************

      subroutine jaminie(mrun,ierr)

c...Purpose: to initialize JAM at each event.
      include 'jam1.inc'
      include 'jam2.inc'

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /pnint/ ipnint

      common /muflag/ imuinthit,imubrmhit,imuppdhit,imucaphit,imucapflag
!$OMP THREADPRIVATE(/muflag/)

c...Check Multi run control.
      mrun=0
      if(mstc(39).gt.0) then
        call jammrun(1,mrun)
        if(mrun.eq.1) return
      endif

c...Make ground state.

      if(mstc(6).gt.-100) then
*-----------------------------------------------------------------------

         if( mstc(190) .eq. 0 ) then

            call jamgrund

         else

            ierr = 0
            call jamqmdgd(ierr)

            if( ierr .ne. 0 ) return

         end if

*-----------------------------------------------------------------------
      else
c...LPC:
        call jaminil
      endif


c...Find impact parameter
c-------------------------

c...Box Doesn't consider.
      if(mstc(4).eq.10) then

c...Uniform distribution.
      else if(parc(4).gt.0) then
         pard(2)=pard(3)+pard(4)*rn(0)

c...b^2 distribution.
      else
        bmin=pard(3)
        bmax=pard(3)+pard(4)
        pard(2)=sqrt(bmin**2 + rn(0)*(bmax**2-bmin**2))
      endif

      if(pard(2).lt.0.0d0) pard(2)=0.0d0
      if(mod(mstd(21),mstd(24)).eq.0) pard(3)=pard(3)+pard(4)
      pard(37)=pard(2)/2
      pard(47)=-pard(2)/2

c...Determine initial condition of the JAM simulation.
c-----------------------------------------------------

c...Box Leave the ground state.
      if(mstc(4).eq.10) then

c...Boost the ground state nuclei.

      else if(mstc(4).lt.10) then

          if(mstc(6).gt.-100) then
            if( mstc(190) .eq. 0 ) then

                call jamboost

            end if
          end if

c...Hyper nuclear reactions.
      else if( mstc(7) .le. 30 ) then
      endif

cABE add @2014/07/22, to make string by photon
      if ((ipnint .ge. 1 .or. imuinthit .eq. 1) .and. ityp .eq. 14) then
         call stringprod
      endif

c...Save total momentum.
      pard(9)=0
      pard(10)=0
      pard(11)=0
      pard(12)=0
      do i = 1,nv
        pard(9)=pard(9)+p(1,i)
        pard(10)=pard(10)+p(2,i)
        pard(11)=pard(11)+p(3,i)
      end do
      pard(9)=pard(9)/mstc(5)
      pard(10)=pard(10)/mstc(5)
      pard(11)=pard(11)/mstc(5)

c..Total kinetic and potential energy.
      if(mstc(6).le.1) then
        ekin=0.0d0
        do i=1,nv
          if(k(1,i).ge.1.and.k(1,i).le.10.or.k(1,i).lt.0) then
               ekin=ekin+p(4+mstc(190)*2,i)
          endif
        end do
        pard(13)=ekin/mstd(11)/mstc(5)

c...Mean field potential.
      else
        call jambuuin(1)
        call jambuud
        call jambuue(ekin,epot,etot)
        pard(13)=etot
      endif

      end

c***********************************************************************

      subroutine jamfin

c...Purpose: to finish up the overall run, to write final information.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      common/jamxml/mult(-500:500)
!$OMP THREADPRIVATE(/jamxml/)

      weit=1d0/dble(mstc(2)*mstc(5))
      if(mstd(56).ne.0) then
        ftim1=pard(88)/mstd(56)
        fnum1=mstd(56)*weit
      else
       ftim1=0.0d0
       fnum1=0.0d0
      endif
      if(mstd(57).ne.0) then
        ftim2=pard(89)/mstd(57)
        fnum2=mstd(57)*weit
      else
       ftim2=0.0d0
       fnum2=0.0d0
      endif

      write(mstc(37),8093)mstc(2),mstc(5),
     1  (pard(71)+pard(72))*weit,
     2   pard(71)*weit,
     3   pard(72)*weit,
     4   pard(73)*weit,
     1   pard(74)*weit, ! BB collision
     1   pard(75)*weit, ! MB collision
     1   pard(76)*weit, ! MM collision
     1   pard(77)*weit, ! antiBB collision
     5   pard(85)*weit, ! h-p coll.
     6   pard(86)*weit, ! p-p coll.
     7   pard(78)*weit,
     8   pard(82)*weit,
     9   pard(83)*weit,
     $   pard(79)*weit,
     1   pard(84)*weit,
     2   pard(80)*weit,
     3   pard(81)*weit,
     4   ftim1,fnum1,ftim2,fnum2

 8093 format(1x,/,27x,'per simulation ',i6,' runs'
     $ ,2x,i6,' parallel runs'/,
     $ 6x,'-----------------------------------',/,
     $ 6x,' Total collisions      = ',f10.4,/,
     $ 6x,' Elastic collisions    = ',f10.4,/,
     $ 6x,' Inelastic collisions  = ',f10.4,/,
     $ 6x,' Absorptions           = ',f10.4,/,
     $ 6x,' baryon-baryon coll.   = ',f10.4,/,
     $ 6x,' meson-baryon coll.    = ',f10.4,/,
     $ 6x,' meson-meson coll.     = ',f10.4,/,
     $ 6x,' anti-B-B coll.        = ',f10.4,/,
     $ 6x,' hadron-parton coll.   = ',f10.4,/,
     $ 6x,' parton-parton coll.   = ',f10.4,/,
     $ 6x,' Resonance decays      = ',f10.4,/,
     $ 6x,' Total particles       = ',f10.4,/,
     $ 6x,' Mesons                = ',f10.4,/,
     $ 6x,' Pauli blocks          = ',f10.4,/,
     $ 6x,' Blocks due to E.con.  = ',f10.4,/,
     $ 6x,' Low energy cuts       = ',g10.4,/,
     $ 6x,' Decay after simul.    = ',f10.4,/,
     $ 6x,' Mean formation time          = ',f10.4,f10.4/,
     $ 6x,' Mean const. formation time   = ',f10.4,f10.4/)

      write(mstc(37),'(/,7x,''average number of jet='',f10.4)')
     $                                   pard(87)/dble(mstc(2))
      write(mstc(37),'(7x,''max. number of nv='',i9)')mstd(81)
      write(mstc(37),'(7x,''max. number of mentry='',i9)')mstd(82)

c....Summarize simulation.

C....Final hadron multiplicities.

      write(mstc(37),6100)
 6100 format(//1x,78('=')/10x,'FINAL HADRON MULTIPLICITIES'/1x,78('=')/)
      write(mstc(37),*)'total particle',mult(0)*weit
      write(mstc(37),*)'charged particle',mult(41)*weit
      write(mstc(37),*)'negative particle',mult(42)*weit
      write(mstc(37),*)'positive particle',mult(42)*weit
      do i=1,500
       if(i.eq.0.or.(i.ge.41.and.i.le.43)) goto 510
       if(mult(i).ge.1.or.mult(-i).ge.1) then
         kf=kchg(i,4)
         if(kchg(i,3).ne.0) then
           write(mstc(37),'(i5,1x,i9,2x,2(a16,2x,g12.3,2x))')
     $      i,kf,chaf(i,1),mult(i)*weit,chaf(i,2),mult(-i)*weit
         else
           write(mstc(37),'(i5,1x,i9,2x,a16,2x,g12.3)')
     $      i,kf,chaf(i,1),mult(i)*weit
         endif
       endif
  510 end do


c...CPU time.
      call jamcpu(1,mstc(1),isec,6)
      imin=isec/60
      ihrs=imin/60
      imin=imin-ihrs*60
      isec=isec-ihrs*3600-imin*60
      write(mstc(37),901)ihrs,imin,isec
901   format(/,6x,'* CPU time=',i4,' h ',i3,' m ',i3,' s')
      write(mstc(37),'(6x,31(''*''))')
      close(mstc(37))

      end

c***********************************************************************

      subroutine jamread(nevent,bmin,bmax,dt,nstep,
     $                             chfram,chbeam,chtarg,cwin)

c...Purpose: to read input configuration file

      include 'jam1.inc'
      include 'jam2.inc'
      common/pypars/mstp(200),parp(200),msti(200),pari(200)
!$OMP THREADPRIVATE(/pypars/)
      logical exex
      parameter(mxchan=30)
      dimension sigin(mxchan)
      character chinp*80,chinp1*80,pname*10,chalp1*26,chalp2*26
      character chfram*8,chbeam*8,chtarg*8,cwin*15
      character chline*5, cherr*90
      character varname*80
      dimension sigt(0:6,0:6,0:5)
c...Functions: momentum in two-particle cm.
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a)
      data chalp1/'abcdefghijklmnopqrstuvwxyz'/
      data chalp2/'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

c...In the case of no input file.
      if(fname(1)(1:1).eq.'0') then
        mstc(2)=nevent
        parc(3)=bmin
        parc(4)=bmax
        parc(2)=dt
        mstc(3)=nstep
        goto 10000
      endif

c...Check if input file really exist.
      inquire(file=fname(1),exist=exex)
      if(exex.eqv..false.) then
         write(6,*) ' error: input cfg file does not exist',fname(1)
         call parastop( 111 )
      end if

c...Open file for input parameter of JAM.
      call readcfg

c...Default values.
      chbeam='28Si'
      chtarg='27Al'
      cwin='14.6gevc'
      chfram='nn'       ! comp. frame  cm, nn, lab, collider,box
      parc(2)=100.d0    ! Time step size(fm/c)
      mstc(3)= 1        ! total number of time step.
      parc(3)=1.5d0     ! Min. impact parameter (fm).
      parc(4)=1.5d0     ! Max. impact parameter (fm).

c=======================================================================

c...Read One Line from the file fname(1).

c...Check pname.
      l=1

      do ivar=1,8
	write(varname,*)'FNAME(',ivar,')'
	call tinput(varname,fname(ivar))
      enddo

      do ivar=1,200
	write(varname,*)'MSTC(',ivar,')'
	call iinput(varname,mstc(ivar))
      enddo

      do ivar=1,200
	write(varname,*)'PARC(',ivar,')'
	call rinput(varname,parc(ivar))
      enddo

      do ivar=1,200
	write(varname,*)'MSTU(',ivar,')'
	call iinput(varname,mstu(ivar))
      enddo

      do ivar=1,200
	write(varname,*)'PARU(',ivar,')'
	call rinput(varname,paru(ivar))
      enddo

      do ivar=1,200
	write(varname,*)'MSTJ(',ivar,')'
	call iinput(varname,mstj(ivar))
      enddo

      do ivar=1,200
	write(varname,*)'PARJ(',ivar,')'
	call rinput(varname,parj(ivar))
      enddo

      do ivar=1,200
	write(varname,*)'MSTP(',ivar,')'
	call iinput(varname,mstp(ivar))
      enddo

      do ivar=1,200
	write(varname,*)'PARP(',ivar,')'
	call rinput(varname,parp(ivar))
      enddo

      call iinput('EVENT',mstc(2))
      call tinput('FRAME',chfram)
      call rinput('BMIN',parc(3))
      call rinput('BMAX',parc(4))
      call tinput('PROJ',chbeam)
      call tinput('TARG',chtarg)
      call rinput('DT',parc(2))
      call tinput('WIN',cwin)
      call iinput('TIMESTEP',mstc(3))


      call checkcfg




c-----------------------------------------------------------------------

200   continue

10000 continue
c...Open file for error message, only once after fisrt call of jaminit.
      if(mstc(38).ge.1.and.mstc(38).ne.6.and.mstc(21).eq.1)
     $      open(mstc(38),file=fname(3),status='unknown')

c...Open file for front page, only once after fisrt call of jaminit.
      if(mstc(37).ge.1.and.mstc(21).eq.1) then
        open(mstc(37),file=fname(2),status='unknown')
      endif

c...Leading particle cascade.
      if(mstc(6).le.-100) chfram='lab'

      call jaminbm(chfram,chbeam,chtarg,cwin)

c=======================================================================

c...Check consistency of swiches

c...Check dimension.
      if(mstc(5)*mstd(11).gt.nnn) then
        write(check(1),'(''mstc5 mstd11 nnn'',i8,1x,i8,1x,i8)')
     $    mstc(5),mstd(11),nnn
        call jamerrm(30,1,'dimension too small:nnn')
      endif

c....Check number of test particles.
      if(mstc(6).ge.2) then
        if(mstc(5).le.1)
     $     call jamerrm(30,0,'test particle too small in BUU mode')
        if(parc(2).ge.10.0d0)
     $     call jamerrm(30,0,'time step size too large in BUU mode')
        if(mstc(3).eq.1)
     $     call jamerrm(30,0,'time step equal to 1 in BUU mode')
      endif

c...Collision according to mean free path.
      if(mstc(6).le.-100) then
        if(mstd(2).ne.1) call jamerrm(30,0,
     $    'projectile should not nucleus with mstc(6)=-100')
        if(mstc(5).ne.1) call jamerrm(30,0,
     $    'test particle should be one with mstc(6)=-100')
      endif

c...In the case of 1+1 simulation
      if (mstd(2) .eq. 1 .and. mstd(5) .eq. 1) then
        if(mstc(6).ge.1) mstc(6)=0  ! forced to cascade mode
        parc(3)=0.0d0               ! bmin
        parc(4)=-0.5d0              ! bmax
        parc(2)=100.0d0             ! dt
        mstc(3)=1                   ! time step
        parc(5)=0.2d0               ! z-separation
        mstc(56)=0                  ! Pauli block off
        pr=pawt(pard(16),pard(34),pard(44))
        srt=pard(16)
        em1=pard(34)
        em2=pard(44)
        kf1=mstd(1)
        kf2=mstd(4)
        pare(3)=0.0d0
        ibar1=isign(kchg(jamcomp(kf1),6),kf1)
        ibar2=isign(kchg(jamcomp(kf2),6),kf2)
        icltyp=jamcltyp(kf1,kf2,ibar1,ibar2)
        call jamcross(1,icltyp,srt,pr,kf1,kf2,em1,em2,
     $                 sig,sigel,sigin,mchanel,mabsrb,ijet,icon)
        sigint=sig
        if(mstc(17).eq.1) sigint=sig-sigel

c...Option for only non-diffractive scattering.
        if(mstc(71).eq.0) then
          call jamxtot(kf1,kf2,srt,pr,sigt)
          sigint=sigt(0,0,5)   ! non-diffractive cross section.
        endif

           if( sigint*0.1d0/paru(1).gt.0.0)then
              parc(4)=-sqrt(sigint*0.1d0/paru(1))
           else
              parc(4)=0.0d0
           end if

      endif

c....AA collision by Glauber geometry.
      if(mstc(6).lt.0) then
        mstc(3)=1        ! time step
        parc(5)=0.0d0    ! z-separation
      endif

c...In the case fo Proj. number =0
      if(mstd(2).eq.0.and.mstc(4).ne.10) then

         if(mstc(6).le.0 .and. mstd(1).ne.22) then
          call jamerrm(30,0,'Cascade mode should not be used'
     $          //' with proj=non. Please change mstc(6)>0')
         endif
         pard(14)=0.0d0    ! incident energy
         parc(3)=0.0d0     ! bmin
         parc(4)=0.0d0     ! bmax

         pard(33)=0.0d0
         pard(31)=0.0d0
         pard(35)=0.0d0
         pard(36)=1.0d0
         pard(39)=0.0d0
         pard(37)=0.0d0

         pard(43)=0.0d0
         pard(41)=0.0d0
         pard(45)=0.0d0
         pard(46)=1.0d0
         pard(49)=0.0d0
         pard(47)=0.0d0

         pard(5)=0.0d0
         pard(6)=1.0d0
         mstc(155)=0     ! flow anal.
         mstc(156)=0     ! energy distribution of collisions
         mstc(162)=0     ! Output collision histroy
         mstc(163)=0     ! time evolution of directed transverse flow
         mstc(165)=0     ! Output time evolution of particle yield
      end if

c...In the case of Cascade mode
      if(mstc(6).le.0) then
         mstc(57)=0      ! recover e.con not necessary
      endif

      gam1=pard(36)
      gam2=pard(46)
      elab=pard(14)
      rat1=dble(mstd(5))/dble(mstd(5)+mstd(2))
      rat2=dble(mstd(2))/dble(mstd(5)+mstd(2))

      if(parc(5).lt.0.0d0) then
         if( elab .ge. 10.0d0 ) then
            rminm = 1.0d0
         else if( elab .ge. 1.0d0 ) then
             rminm = 2.0d0
         else if( elab .ge. 0.2d0 ) then
            rminm = 4.0d0
         else
            rminm = 6.0d0
         endif
         pard(39)=-(pard(40)/gam1+rminm*rat1)
         pard(49)=pard(50)/gam2 + rminm*rat2
      else
         pard(39)=-(pard(40)/gam1 + parc(5)*rat1)
         pard(49)=pard(50)/gam2 + parc(5)*rat2
      end if

c...Shift distance of z.

c----------------------------------------------------------------------c

      if(mstc(7).le.10) then
      else if(mstc(7).le.30) then
      end if

c...In case of proj.=0.
      if(mstd(2).eq.0) then
        pard(14)=0.0d0
        parc(3)=0.0d0
        parc(4)=0.0d0
      end if

c----------------------------------------------------------------------c
c      set initial condition of collisions
c----------------------------------------------------------------------c

c...Set max. c.m. distance for collisions.
      pard(51) = sqrt( 0.1d0*parc(32)/paru(1) )   ! pn
      pard(52) = sqrt( 0.1d0*parc(33)/paru(1) )   ! pp
      pard(53) = sqrt( 0.1d0*parc(34)/paru(1) )   ! BB
      pard(54) = sqrt( 0.1d0*parc(35)/paru(1) )   ! MB
      pard(55) = sqrt( 0.1d0*parc(36)/paru(1) )   ! MM
      pard(56) = sqrt( 0.1d0*parc(37)/paru(1) )   ! aBB

      return

c----------------------------------------------------------------------c
c...Error occured in reading number


      end

c***********************************************************************
c                                                                      *
c      INPUT part                                                      *
c      written by Toshiki Maruyama   30 Apr 1999.                      *
c                                                                      *
c                                                                      *
c***********************************************************************

c***********************************************************************

      subroutine adjust(t)

c...Get the first word from a text.
      implicit double precision(a-h, o-z)
      character*(*) t

      leng=min(len(t),100)
      do i1=1,leng
	if(t(i1:i1).ne.' ') goto 100
      enddo
      return
100   continue
      do i2=i1,leng
	if(t(i2:i2).eq.' ') goto 200
      enddo
200   continue
      i2=i2-1
      t=t(i1:i2)

      end

c***********************************************************************

      subroutine compress(t)

c...Compress to remove blanc characters from a text.
      implicit double precision(a-h, o-z)
      character*(*) t
      leng=min(len(t),100)
      do i1=1,leng-1
	if(t(i1:i1).eq.' ') then
	  do i2=i1,leng-1
	    if(t(i2:i2).ne.' ') then
	      t(i1:)=t(i2:)
	    endif
	  enddo
	endif
      enddo

      end

c***********************************************************************

      subroutine caps(t)

c...Convert a text into capital letters.
      implicit double precision(a-h, o-z)
      character*(*) t
      character*1 c
      leng=min(len(t),100)
      do i=1,leng
	c=t(i:i)
	if(c.ge.'a'.and.c.le.'z')
     &      c=char(ichar(c)+ichar('A')-ichar('a'))
	t(i:i)=c
      enddo

      end

c***********************************************************************

      subroutine checkcfg

c...Check errors in the config file.
      implicit double precision(a-h, o-z)
      common /config/nline,ifused,linepos,nameline,varline
!$OMP THREADPRIVATE(/config/)
      dimension ifused(1000),linepos(1000)
      character*80 nameline(1000),varline(1000)
      character*90 cherr

      nerr=0
      do i=1,nline
        if(ifused(i).ne.1) then

          nerr=nerr+1
          leng=index(nameline(i),' ')-1
          if(ifused(i).eq.0) then
            write(6,'(a,a,a,i4)') 'unknown [',
     $      nameline(i)(:leng),'] at line',linepos(i)
          else
            write(6,'(a,a,a,i4)') 'doubly defiend [',
     $      nameline(i)(:leng),'] at line',linepos(i)
          endif
        endif

      enddo

      write(6,'(i3,a,i3,a)') nline,' parameter(s),',
     &                       nerr,' error(s) in config file'
      if(nerr.ge.1) call parastop( 111 )


      end

c***********************************************************************

      subroutine readcfg

c...Open and read the config file.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      common /config/nline,ifused,linepos,nameline,varline
!$OMP THREADPRIVATE(/config/)
      dimension ifused(1000),linepos(1000)
      character*80 nameline(1000),varline(1000),templine
      character*80 FNCFG

      idcfg=mstc(36)
      FNCFG=fname(1)
      linepos0=0
      open(idcfg,file=FNCFG,status='OLD')
      leng=index(FNCFG,' ')-1
      write(6,*) 'config file=[',FNCFG(:leng),']'
1000  continue
	read(idcfg,'(a)',iostat=ios,err=999) templine
         if( ios .eq. -1 ) goto 888
	linepos0=linepos0+1

	icomment1=index(templine,'#')
	icomment2=index(templine,'!')
	if(icomment1*icomment2.ne.0) then
	  icomment=min(icomment1,icomment2)
        else
	  icomment=max(icomment1,icomment2)
	endif
	ii=index(templine,'=')

	if((icomment.eq.0.or.icomment.gt.ii+1)
     &               .and.ii.ge.2.and.ii.lt.80) then
	  nline=nline+1
	  linepos(nline)=linepos0
          nameline(nline)=templine(:ii-1)
	  call caps(nameline(nline))
	  call compress(nameline(nline))
	  call adjust(nameline(nline))
          varline(nline)=templine(ii+1:)
	  call adjust(varline(nline))
        else
	  call caps(templine)
	  call adjust(templine)
	  if(templine.eq.'END') then
	    goto 888
	  endif
        endif
      goto 1000
999   continue
      write(6,*) 'file read error'
888   continue
      write(6,*) 'cfg close'
      close(idcfg)

      end

c***********************************************************************

      subroutine iinput(name,ivar)

c...Input integer variable from the config.
      implicit double precision(a-h, o-z)
      common /config/nline,ifused,linepos,nameline,varline
!$OMP THREADPRIVATE(/config/)
      dimension ifused(1000),linepos(1000)
      character*80 nameline(1000),varline(1000),templine
      character*(*) name

      templine=name

      call compress(templine)
      leng1=index(templine,' ')-1
      if(leng1.le.0) leng1=len(templine)
      nfound=0
      do i=1,nline
	leng=index(nameline(i),' ')-1
	if(nameline(i)(:leng).eq.templine(:leng1)) then
	  read(varline(i),*,err=999) ivar
          nfound=nfound+1
	  ifused(i)=nfound
	endif
      enddo
      return

999   continue
      write(6,*) 'bad format:',nameline(i)(:leng),varline(i)
      end


c***********************************************************************

      subroutine rinput(name,rvar)

c...Input real variable from the config.

      implicit double precision(a-h, o-z)
      common /config/nline,ifused,linepos,nameline,varline
!$OMP THREADPRIVATE(/config/)
      dimension ifused(1000),linepos(1000)
      character*80 nameline(1000),varline(1000),templine
      character*(*) name

      templine=name

      call compress(templine)
      leng1=index(templine,' ')-1
      if(leng1.le.0) leng1=len(templine)
      nfound=0
      do i=1,nline
	leng=index(nameline(i),' ')-1
	if(nameline(i)(:leng).eq.templine(:leng1)) then
	read(varline(i),*,err=999) rvar
        nfound=nfound+1
	ifused(i)=nfound
	endif
      enddo
      return

999   continue
      write(6,*) 'bad format:',nameline(i)(:leng),varline(i)
      end


c***********************************************************************

      subroutine tinput(name,tvar)

c...Input text variable from the config.
      implicit double precision(a-h, o-z)
      common /config/nline,ifused,linepos,nameline,varline
!$OMP THREADPRIVATE(/config/)
      dimension ifused(1000),linepos(1000)
      character*80 nameline(1000),varline(1000),templine
      character*(*) name
      character*(*) tvar

      templine=name

      call compress(templine)
      leng1=index(templine,' ')-1
      if(leng1.le.0) leng1=len(templine)
      nfound=0
      do i=1,nline
	leng=index(nameline(i),' ')-1
	if(nameline(i)(:leng).eq.templine(:leng1)) then
          nfound=nfound+1
	  ifused(i)=nfound
	  lengt=index(varline(i),' ')-1
	  tvar=varline(i)(:lengt)
	endif
      enddo

      end

c*********************************************************************

      subroutine jaminbm(chfram,chbeam,chtarg,cwini)

c...Identifies the two incoming particles and the choice of frame.
c...Modified for box.

      include 'jam1.inc'
      include 'jam2.inc'


      common /jamnmtc/ kfjam
!$OMP THREADPRIVATE(/jamnmtc/)

      character chfram*8,chbeam*8,chtarg*8,cwini*15,cwin*30
      character chcom(3)*8,chalp(2)*26,chidnt(3)*8,chtemp*8
      parameter(mxbeam=34)
      character chcde(mxbeam)*8,chinit*76
      dimension len(3),kcde(mxbeam),pm(2),inuc(2),ma(2),lchcode(mxbeam)
      character element(108)*3,num(0:9)*1,cnuc*3
c...Functions: momentum in two-particle cm.
      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a)

      data chalp/'abcdefghijklmnopqrstuvwxyz',
     &           'ABCDEFGHIJKLMNOPQRSTUVWXYZ'/

      data chcde/'e-      ','e+      ','nu_e    ','nu_e~   ',
     &'mu-     ','mu+     ','nu_mu   ','nu_mu~  ','tau-    ',
     &'tau+    ','nu_tau  ','nu_tau~ ','pi+     ','pi-     ',
     &'n0      ','n~0     ','p+      ','p~-     ','gamma   ',
     &'lambda0 ','sigma-  ','sigma0  ','sigma+  ','xi-     ',
     &'xi0     ','omega-  ','pi0     ','reggeon ','pomeron ',
     &'k-      ','k+      ','kzero   ','other   ','non     '/
      data lchcode/2,2,4,5,3,3,5,6,4,4,6,7,3,3,2,3,2,3,5,
     $ 7,6,6,6,3,3,6,3,7,7,2,2,5,5,3/

      data kcde/11,-11,12,-12,13,-13,14,-14,15,-15,16,-16,
     & 211,-211,2112,-2112,2212,-2212,22,3122,3112,3212,3222,
     & 3312,3322,3334,111,28,29,-321,321,311,0,-9999/

      data element/
     & 'h  ','he ','li ','be ','b  ','c  ','n  ','o  ',
     & 'f  ','ne ','na ','mg ','al ','si ','p  ','s  ',
     & 'cl ','ar ','k  ','ca ','sc ','ti ','v  ','cr ',
     & 'mn ','fe ','co ','ni ','cu ','zn ','ga ','ge ',
     & 'as ','se ','br ','kr ','rb ','sr ','y  ','zr ',
     & 'nb ','mo ','tc ','ru ','rh ','pd ','ag ','cd ',
     & 'in ','sn ','sb ','te ','i  ','xe ','cs ','ba ',
     & 'la ','ce ','pr ','nd ','pm ','sm ','eu ','gd ',
     & 'tb ','dy ','ho ','er ','tm ','yb ','lu ','hf ',
     & 'ta ','w  ','re ','os ','ir ','pt ','au ','hg ',
     & 'tl ','pb ','bi ','po ','at ','rn ','fr ','ra ',
     & 'ac ','th ','pa ','u  ','np ','pu ','am ','cm ',
     & 'bk ','cf ','es ','fm ','md ','no ','lr ','ku ',
     $ '105','106','107','108'/
      data num/'0','1','2','3','4','5','6','7','8','9'/

      cwin(1:15) = cwini(1:15)

C...Convert character variables to lowercase and find their length.
      chcom(1)=chbeam
      chcom(2)=chtarg
      chcom(3)=chfram
      do 130 i=1,3

        len(i)=8
        do ll=8,1,-1
          if(len(i).eq.ll.and.chcom(i)(ll:ll).eq.' ') len(i)=ll-1
        do la=1,26
          if(chcom(i)(ll:ll).eq.chalp(2)(la:la))
     $        chcom(i)(ll:ll)=chalp(1)(la:la)
        end do
        end do
        chidnt(i)=chcom(i)
  130 continue

c...Loop over proj. and targ.
c-------------------------------------
      do i=1,2
c-------------------------------------
        inuc(i)=0
        ma(i)=0
c...Check if nucleus or hadrons
        do 1010 l=1,4
c...Box:maru
          if(chidnt(i)(l:l).eq.':') goto 1010
          do ln=0,9
            if(chidnt(i)(l:l).eq.num(ln)) inuc(i)=inuc(i)+1
          end do
 1010   continue

c....Nucleus case
        if(inuc(i).ge.1) then
          read(chidnt(i)(1:inuc(i)),*)ma(i)
        else

c...Fix up bar, underscore and charge in particle name (if needed).
          do ll=1,6
            if(chidnt(i)(ll:ll+2).eq.'bar') then
              chtemp=chidnt(i)
              chidnt(i)=chtemp(1:ll-1)//'~'//chtemp(ll+3:8)//'  '
            endif
          end do

          if(chidnt(i)(1:2).eq.'nu'.and.chidnt(i)(3:3).ne.'_') then
            chtemp=chidnt(i)
            chidnt(i)='nu_'//chtemp(3:7)
          elseif(chidnt(i)(1:2).eq.'n ') then
            chidnt(i)(1:3)='n0 '
          elseif(chidnt(i)(1:2).eq.'n~') then
            chidnt(i)(1:3)='n~0'
          elseif(chidnt(i)(1:2).eq.'p ') then
            chidnt(i)(1:3)='p+ '
          elseif(chidnt(i)(1:2).eq.'p~'.or.chidnt(i)(1:2).eq.'p-') then
            chidnt(i)(1:3)='p~-'
          elseif(chidnt(i)(1:6).eq.'lambda') then
            chidnt(i)(7:7)='0'
          elseif(chidnt(i)(1:3).eq.'reg') then
            chidnt(i)(1:7)='reggeon'
          elseif(chidnt(i)(1:3).eq.'pom') then
            chidnt(i)(1:7)='pomeron'
          elseif(chidnt(i)(1:3).eq.'non') then
            chidnt(i)(1:3)='non'
          endif

      endif
c-------------------------------------
      end do
c-------------------------------------


c...Identify incoming beam and target particles.
      do 150 i=1,2

c....Hadron beam
        if(inuc(i).eq.0) then

          mstd(i*3-2)=0
          do j=1,mxbeam
            ll=lchcode(j)
            if(chidnt(i)(1:ll).eq.chcde(j)(1:ll)) mstd(i*3-2)=kcde(j)
          end do

          if( mstd(i*3-2) .eq. -9999 ) then

            mstd(i*3-2)=0
            mstd(i*3-1)=0
            mstd(i*3)=0
            pard(i*10+24)=0.0d0

          else
            if( mstd(i*3-2) .eq. 0 ) mstd(i*3-2) = kfjam

cABE change @2014/07/22, for photon
            if(mstd(i*3-2) .eq. 22) then
              pm(i) = 0.0d0
              mstd(i*3-1) = 0
              mstd(i*3) = 0
              pard(i*10+24) = pm(i)
            else
              pm(i) = pymass(mstd(i*3-2))
              mstd(i*3-1) = 1
              mstd(i*3) = jamchge(mstd(i*3-2))/3
              pard(i*10+24) = pm(i)
            endif

          endif


c....Nucleus
        else

          cnuc=chidnt(i)(inuc(i)+1:len(i))
     $            //'                              '

c...Box maru
          if(cnuc(1:1).eq.':')then
            read(cnuc(2:),*,err=9977,iostat=ios)iz
            if( ios .eq. -1 ) goto 9977
            goto 155
9977        write(6,*) 'Z_nuc not input. I use Z of half mass.'
            iz=ma(i)/2
            goto 155
          endif

          do j=1,108
            if(cnuc(1:3).eq.element(j)(1:3)) then
              iz=j
              goto 155
            endif
          end do
          write(6,*)'Unknown nucleus code ',cnuc
          call parastop( 111 )
  155     continue
          mstd(i*3-2)=1000000*(ma(i)-iz)+1000*iz+1000000000
          mstd(i*3-1)=ma(i)
          mstd(i*3)=iz
          pm(i)=(pymass(2112)+pymass(2212))/2.d0
          pard(i*10+24)=pm(i)
        endif
  150 continue

c....Check if proj. and targ. are defined.
      if(mstd(1).eq.0.and.mstd(4).eq.0) then
        if(mstd(1).eq.0) write(6,5000) chtarg(1:len(3))
        if(mstd(4).eq.0) write(6,5100) chtarg(1:len(3))
        write(check(1),'(''mstd1'',3(i9,1x),g15.7)')
     $  mstd(1),mstd(2),mstd(3),pard(34)
        write(check(2),'(''mstd4'',3(i9,1x),g15.7)')
     $  mstd(4),mstd(5),mstd(6),pard(44)
        call jamerrm(30,2,'(jaminbm:) Input incorrect')
      endif

c...Total number of particle.
      mstd(11)=mstd(2)+mstd(5)

c...Proj.
      if(inuc(1).ge.1) then
        ib1=mstd(2)
        iz1=mstd(3)
        is1=0
      else
        if(mstd(1).ne.0) then
          kfsg=isign(1,mstd(1))
          kc=jamcomp(mstd(1))
          iz1=kchg(kc,1)*kfsg/3
          ib1=kchg(kc,6)*kfsg/3
          is1=kchg(kc,7)*kfsg
        else
          iz1=0
          ib1=0
          is1=0
        endif
      endif

c...Targ.
      if(inuc(2).ge.1) then
        ib2=mstd(5)
        iz2=mstd(6)
        is2=0
      else
        if(mstd(4).ne.0) then
          kfsg=isign(1,mstd(4))
          kc=jamcomp(mstd(4))
          iz2=kchg(kc,1)*kfsg/3
          ib2=kchg(kc,6)*kfsg/3
          is2=kchg(kc,7)*kfsg
        else
          iz2=0
          ib2=0
          is2=0
        endif
      endif

c...Total baryon number.
      mstd(12)=ib1+ib2
c...Total charge.
      mstd(13)=iz1+iz2
c...Total strangeness.
      mstd(14)=is1+is2

      if(mstd(2).eq.0.or.mstd(5).eq.0)then
        pard(14)=0.0d0
        pard(15)=0.0d0
        pard(16)=0.0d0
        pard(17)=0.0d0
        pard(18)=0.0d0

        pard(33)=0.0d0
        pard(31)=0.0d0
        pard(35)=0.0d0
        pard(36)=1.0d0
        pard(37)=0.0d0

        pard(43)=0.0d0
        pard(41)=0.0d0
        pard(45)=0.0d0
        pard(46)=1.0d0
        pard(47)=0.0d0
        goto 200
      endif

c...Set beam energy and momentum.
      do il=1,15
        do m=1,26
        if(chalp(2)(m:m).eq.cwin(il:il)) cwin(il:il)=chalp(1)(m:m)
        end do
      end do
      emproj=pard(34)
      emtarg=pard(44)
      do il=1,26
       if(cwin(il:il+3).eq.'gevc') then
         read(cwin(1:il-1),'(g20.0)')pproj
         elab=sqrt(emproj**2+pproj**2)-emproj
         goto 10
       else if(cwin(il:il+2).eq.'gev') then
         read(cwin(1:il-1),'(g20.0)')elab
         pproj=sqrt(elab*(2.d0*emproj+elab))
         goto 10
       else if(cwin(il:il+3).eq.'mevc') then
         read(cwin(1:il-1),'(g20.0)')pproj
         pproj=pproj/1000
         elab=sqrt(emproj**2+pproj**2)-emproj
         goto 10
       else if(cwin(il:il+2).eq.'mev') then
         read(cwin(1:il-1),'(g20.0)')elab
         elab=elab/1000
         pproj=sqrt(elab*(2.d0*emproj+elab))
         goto 10
       endif
      enddo
      write(6,*)'Input invalid energy ',cwin
      call parastop( 111 )
10    continue


c...Identify choice of frame and input energies.
      chinit=' '
         icell=13
         do icellx=-1,1
         do icelly=-1,1
         do icellz=-1,1
           icell=icell+1
           rcell(1,mod(icell,27)+1)=0.0d0
           rcell(2,mod(icell,27)+1)=0.0d0
           rcell(3,mod(icell,27)+1)=0.0d0
         enddo
         enddo
         enddo

C...Events defined in the CM frame.
      n1=mstd(2)
      n2=mstd(5)
      if(chcom(3)(1:2).eq.'cm') then
        mstc(4)=1
        emt1=pm(1)*n1
        emt2=pm(2)*n2
        etotaj=elab*n1+emt1+emt2
        plab=pproj*n1

        ylab=0.5d0*log((etotaj+plab)/(etotaj-plab))
        s2tot=etotaj**2-plab**2
        stot=sqrt(s2tot)
        pard(5)=plab/etotaj
        pard(6)=etotaj/stot
        if(stot.gt.emt1+emt2) then
          pst=sqrt((s2tot-(emt1+emt2)**2)*(s2tot-(emt1-emt2)**2))
     $             /(2*stot)
          p01=pst/n1
          p02=-pst/n2
        else
          p01=0.0d0
          p02=0.0d0
        endif
        e1=sqrt(pm(1)**2+p01**2)
        e2=sqrt(pm(2)**2+p02**2)
        srt=sqrt((e1+e2)**2-(p01+p02)**2)

        if(mstc(8).ge.1) then
          loffs=(31-(len(1)+len(2)))/2
          chinit(loffs+1:76)='JAM will be initialized for '//
     &    chcom(1)(1:len(1))//' on '//chcom(2)(1:len(2))//
     &    ' in C.M. frame'//' '
          write(6,5200) chinit
          write(6,5300) elab
          write(6,5500) srt
        endif

C...Events defined in fixed target frame.
      elseif(chcom(3)(1:3).eq.'lab') then

        mstc(4)=0

        srt=sqrt((elab+pm(1)+pm(2))**2-pproj**2)
        p01=pproj
        p02=0.0d0
        e1=sqrt(pm(1)**2+p01**2)
        e2=sqrt(pm(2)**2+p02**2)
        ylab=0.0d0

        pard(5)=0.0d0
        pard(6)=1.0d0

        if(mstc(8).ge.1) then
          loffs=(29-(len(1)+len(2)))/2
          chinit(loffs+1:76)='JAM will be initialized for '//
     &    chcom(1)(1:len(1))//' on '//chcom(2)(1:len(2))//
     &    ' Lab. frame'//' '
          write(6,5200) chinit
          write(6,5300) elab
          write(6,5500) srt
        endif

C...Frame defined in NN CM.
      elseif(chcom(3)(1:2).eq.'nn') then
        mstc(4)=2
        ee=elab+emproj+emtarg
        s2=ee**2-pproj**2
        srt=sqrt(s2)
        p01=sqrt((s2-(pm(1)+pm(2))**2)*(s2-(pm(1)-pm(2))**2))
     $             /(2*srt)
        p02=-p01
        e1=sqrt(pm(1)**2+p01**2)
        e2=sqrt(pm(2)**2+p02**2)

        ylab=0.5d0*log((ee+pproj)/(ee-pproj))
        pard(5)=pproj/ee
        pard(6)=ee/srt

        if(mstc(8).ge.1) then
          loffs=(29-(len(1)+len(2)))/2
          chinit(loffs+1:76)='JAM will be initialized for '//
     &    chcom(1)(1:len(1))//' on '//chcom(2)(1:len(2))//
     &    ' nn c.m. frame'//' '
          write(6,5200) chinit
          write(6,5300) elab
          write(6,5500) srt
        endif

c...Collider   elab means c.m. energy per nucl.
      elseif(chcom(3)(1:8).eq.'collider') then
         mstc(4)=3
         srt=abs(elab)
         elab=(srt*srt-4*pm(1)*pm(2))/(2*pm(2))
         s2=srt*srt
         p01=sqrt((s2-(pm(1)+pm(2))**2)*(s2-(pm(1)-pm(2))**2))
     $             /(2*srt)
         p02=-p01
         e1=sqrt(pm(1)**2+p01**2)
         e2=sqrt(pm(2)**2+p02**2)


         pproj=p01
         ee=elab+pm(1)+pm(2)
         dbeta=pproj/ee
         ylab=0.5d0*log((1.d0+dbeta)/(1.d0-dbeta))
         pard(5)=pproj/ee
         pard(6)=ee/pm(1)


        if(mstc(8).ge.1) then
          loffs=(29-(len(1)+len(2)))/2
          chinit(loffs+1:76)='JAM will be initialized for '//
     &    chcom(1)(1:len(1))//' on '//chcom(2)(1:len(2))//
     &    ' for collider'//' '
          write(6,5200) chinit
          write(6,5400) pproj
          write(6,5500) srt
        endif

c...Box maru
      elseif(chcom(3)(1:3).eq.'box') then
         mstc(4)=10
         mstc(54)=0  ! avoid first coll inside the same nucleus off
         mstc(41)=0  ! resonances is not decay after simul.
         read(chcom(3)(4:),*,err=9988,iostat=ios)dencell
         if( ios .eq. -1 ) goto 9988
         goto 9989
9988     write(6,*) 'Cell density not input. Set it by, eg. "box1.0".'
           dencell=1.0d0
9989      continue
         drcell=(mstd(12)/dencell/parc(21))**(1.d0/3)
         pard(21)=drcell
         pard(22)=elab

         pm(1)=parc(28)
         pm(2)=parc(28)

         ee=abs(elab)
         e1=ee+pm(1)
         e2=ee+pm(2)
         p01=sqrt(e1**2-pm(1)**2)
         p02=-sqrt(e2**2-pm(2)**2)
         pproj=p01
         srt=e1+e2
         ylab=0.0d0
         pard(5)=0.0d0
         pard(6)=1.0d0

         icell=13
         do icellx=-1,1
         do icelly=-1,1
         do icellz=-1,1
           icell=icell+1
           rcell(1,mod(icell,27)+1)=drcell*icellx
           rcell(2,mod(icell,27)+1)=drcell*icelly
           rcell(3,mod(icell,27)+1)=drcell*icellz
         enddo
         enddo
         enddo

        if(mstc(8).ge.1) then
          loffs=(29-(len(1)+len(2)))/2
          chinit(loffs+1:76)='JAM will be initialized for '//
     &    chcom(1)(1:len(1))//' on '//chcom(2)(1:len(2))//
     &    ' for box'//' '
          write(6,5200) chinit
          write(6,5300) elab
        endif


C...Unknown frame. Error for too low CM energy.
      else
        write(6,5800) chfram(1:len(3))
        call parastop( 111 )
      endif


      if(mstd(2).eq.0) then
        gam1=1.0d0
        beta1=0.0d0
      else
        gam1=e1/pm(1)
        beta1=p01/e1
      endif
      if(mstd(5).eq.0) then
        gam2=1.0d0
        beta2=0.0d0
      else
        gam2=e2/pm(2)
        beta2=p02/e2
      endif

      pard(14)=elab
      pard(15)=pproj
      pard(16)=srt
      pard(17)=ylab
      pard(18)=pawt(srt,pm(1),pm(2))

      pard(33)=p01
      pard(31)=0.0d0
      pard(35)=beta1
      pard(36)=gam1
      pard(34)=emproj
      pard(37)=0.0d0

      pard(43)=p02
      pard(41)=0.0d0
      pard(45)=beta2
      pard(46)=gam2
      pard(44)=emtarg
      pard(47)=0.0d0

 200  continue
      if(mstd(2).gt.2) then
      pard(40)=1.19d0*dble(mstd(2))**(1.d0/3.d0)
     $        - 1.61d0*dble(mstd(2))**(-1.d0/3.d0)
      else
       pard(40)=0.0d0
      endif

      if(mstd(5).gt.2) then
      pard(50)=1.19d0*dble(mstd(5))**(1.d0/3.d0)
     $        - 1.61d0*dble(mstd(5))**(-1.d0/3.d0)
      else
       pard(50)=0.0d0
      endif


C...Formats for initialization and error information.
 5000 format(1x,'Error: unrecognized beam particle ''',a,'''.'/
     &1x,'Execution stopped!')
 5100 format(1x,'Error: unrecognized target particle ''',a,'''.'/
     &1x,'Execution stopped!')

 5200 format(/1x,78('=')/1x,'I',76x,'I'/1x,'I',a76,'I')

 5300 format(1x,'I',18x,'at',1x,f10.3,1x,'GeV Lab. energy',
     &19x,'I'/1x,'I',76x,'I'/1x,78('='))

 5400 format(1x,'I',22x,'at',1x,f10.3,1x,'GeV/c lab-momentum',22x,'I')

 5500 format(1x,'I',76x,'I'/1x,'I',11x,'corresponding to',1x,f10.3,1x,
     &'GeV center-of-mass energy',12x,'I'/1x,'I',76x,'I'/1x,78('='))

 5600 format(1x,'I',76x,'I'/1x,'I',18x,'px (GeV/c)',3x,'py (GeV/c)',3x,
     &'pz (GeV/c)',6x,'E (GeV)',9x,'I')
 5700 format(1x,'I',8x,a8,4(2x,f10.3,1x),8x,'I')

 5800 format(1x,'Error: unrecognized coordinate frame ''',a,'''.'/
     &1x,'Execution stopped!')

      end

c***********************************************************************

      subroutine jamlogo

c...Purpose: to print front page and echo of input
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      character cproj*16,ctarg*16
c...Parameter for length of information block.
      parameter (irefer=17,mrefer=5)
c...Local arrays and character variables.
      integer idati(6)
      character month(12)*3, logo(48)*32, refer(2*mrefer)*36, line*79,
     &vers*1, subv*3, date*2, year*4, hour*2, minu*2, seco*2
c...Data on months, logo, titles, and references.
      data month/'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep',
     &'Oct','Nov','Dec'/
      data (logo(j),j=20,38)/
     &'                                ',
     &'                                ',
     &'                                ',
     &'                                ',
     &'                                ',
     &'                                ',
     &'                                ',
     &'                                ',
     &'                                ',
     &'                                ',
     &'                                ',
     &'                                ',
     &'                                ',
     &'                                ',
     &'                                ',
     &'                                ',
     &'                                ',
     &'                                ',
     &'                                '/
      data (logo(j),j=1,19)/
     &' Welcome to the JAM Monte Carlo!',
     &'                                ',
     &' JJJJJJJJ  AAAA      MM     MMM ',
     &'    JJ    AA  AA     MMM    MMM ',
     &'    JJ   AAAAAAAA    MM MMMM MM ',
     &' J  JJ  AA      AA   MM      MM ',
     &' JJJJJ AA        AA  MM      MM ',
     &'                                ',
     &'  This is JAM  version x.xxx    ',
     &'Last date of change: xx xxx 199x',
     &'                                ',
     &'Now is xx xxx 199x at xx:xx:xx  ',
     &'                                ',
     &'JAM is a free software which can',
     &'be distributed under the GNU Gen',
     &'eral Public License. There is no',
     &'warranty for the program.       ',
     &'                                ',
     &'Copyright (c) 1998 Yasushi Nara '/

      data (refer(j),j=1,2*mrefer)/
     &' Main author: Yasushi Nara;         ',
     &'                                    ',
     &' Advance Science Research Center, Ja',
     &'pan Atomic Energy Research Institute',
     &' Tokai-mura, Naka-gun, Ibaraki-ken  ',
     &'319-1195, Japan                     ',
     &' e-mail: ynara@hadron03.tokai.jaeri.',
     &'go.jp                               ',
     &'                                    ',
     &'                                    '/

      io=mstc(37)

c...Check that JAM linked.
      if(mstc(193)/10.ne.199) then
        write(*,'(1x,a)')
     &  'Error: JAM has not been linked.'
        write(*,'(1x,a)') 'Execution stopped!'
        call parastop( 111 )

c...Write current version number and current date+time.
      else
        write(vers,'(I1)') mstc(191)
        logo(9)(24:24)=vers
        write(subv(1:1),'(I1)') mod(mstc(192)/100,10)
        logo(9)(26:26)=subv(1:1)
        write(subv(1:1),'(I1)') mod(mstc(192)/10,10)
        logo(9)(27:27)=subv(1:1)
        write(subv(1:1),'(I1)') mod(mstc(192),10)
        logo(9)(28:28)=subv(1:1)

        write(date,'(I2)') mstc(195)
        logo(10)(22:23)=date
        logo(10)(25:27)=month(mstc(194))
        write(year,'(I4)') mstc(193)
        logo(10)(29:32)=year

        idati(1)=0
        if(idati(1).le.0) then
          logo(12)='                                '
        else
          write(date,'(I2)') idati(3)
          logo(12)(8:9)=date
          logo(12)(11:13)=month(max(1,min(12,idati(2))))
          write(year,'(I4)') idati(1)
          logo(12)(15:18)=year
          write(hour,'(I2)') idati(4)
          logo(12)(23:24)=hour
          write(minu,'(I2)') idati(5)
          logo(12)(26:27)=minu
          if(idati(5).lt.10) logo(31)(26:26)='0'
          write(seco,'(I2)') idati(6)
          logo(12)(29:30)=seco
          if(idati(6).lt.10) logo(12)(29:29)='0'
        endif
      endif

c...Loop over lines in header. Define page feed and side borders.
      do 100 ilin=1,25+mrefer+4
        line=' '
        if(ilin.eq.1) then
          line(1:1)=' '
        else
          line(2:3)='**'
          line(78:79)='**'
        endif

c...Separator lines and logos.
        if(ilin.eq.2.or.ilin.eq.3.or.ilin.gt.27+mrefer) then
          line(4:77)='***********************************************'//
     &    '***************************'
        elseif(ilin.ge.6.and.ilin.le.24) then
          line(6:37)=logo(ilin-5)
          line(44:75)=logo(ilin+14)
        elseif(ilin.ge.26.and.ilin.lt.26+mrefer) then
          line(5:40)=refer(2*ilin-51)
          line(41:76)=refer(2*ilin-50)
        endif

c...Write lines to appropriate unit.
        write(io,'(a79)') line
  100 continue


      if(mstd(2).gt.0) then

c...Pint collision type.
         write(io,860)
         write(io,870)mstd(2),mstd(3),mstd(2)-mstd(3),
     &                mstd(5),mstd(6),mstd(5)-mstd(6)

      write(io,880)
      call jamname(mstd(1),0,0,cproj)
      call jamname(mstd(4),0,0,ctarg)

      write(io,890)cproj(1:7),ctarg(1:7)
      write(io,880)

         if(pard(14).gt.1.0d0) then
           write(io,900)pard(14),pard(15),pard(16)
         else
           write(io,901)pard(14)*1000,pard(15)*1000,pard(16)
         endif

      write(io,880)
      write(io,902)
      write(io,903) pard(35),pard(45)
      write(io,904) pard(36),pard(46)
      write(io,905) pard(31),pard(41)
      write(io,906) pard(33),pard(43)
      write(io,907) pard(37),pard(47)
      write(io,908) pard(39),pard(49)


c...Impact parameter.
      write(io,880)
      if(mstc(10).gt.1) then
        write(io,910)parc(3),parc(4),pard(4),mstd(24)
      else
        if(parc(4).ge.0.0d0) then
           write(io,911) parc(3),parc(4)
        else
           write(io,912) parc(3),abs(parc(4))
        endif
      endif
      write(io,'(6x,''*'',65x,''*''/6x,67(''*''))')
      write(io,'(/)')

      else
         write(io,'(//6x,67(''*'')/6x,''*'',65x,''*'')')
         write(io,920)mstd(5),mstd(6),mstd(5)-mstd(6)

         if(mstc(4).eq.10) then
           write(io,921)
         else
           write(io,922)
         endif
      end if
c----------------------------------------------------------------------*

 600  format(/6x,'* Input output file name: ','"',a15,'"')

c...Write input file name.
      if(fname(1).ne.'0') then
        write(io,600) fname(1)
      else
        write(io,600) 'non'
      endif

c...Write comp. frame.
      if(mstc(4).eq.0) then
         write(io,605)
      else if(mstc(4).eq.1) then
         write(io,606)
      else if(mstc(4).eq.2) then
         write(io,607)
      else if(mstc(4).eq.3) then
         write(io,608)
      end if

c...Write time step,event, .....
      write(io,616)mstc(1),mstc(2),parc(2),mstc(3)
      write(io,'(//)')

c...Format
c---------------------------------------------------------------
 860  format(//6x,67('*')/6x,'*',65x,'*'/
     &             6x,'*',10x,'Reaction:',46x,'*'/
     &             6x,'*',65x,'*')
 870  format(6x,'*',14x,'mass',i4,
     &  '(',i3,',',i3,')',
     &  ' ==> mass',i4,'(',i3,',',i3,')',12x,'*')
 880  format(6x,'*',65x,'*')
 890  format(6x,'*',14x,'Collision of ',a7,' on ',a7,
     $   20x,'*')
 910  format(6x,'*',10x,
     $ 'Impact parameters are distributed: ',20x,'*'
     $  /6x,'*',65x,'*',
     $  /6x,'*',f7.3,' < b < ',f7.3,
     $ ' (fm) with bin size',f7.3,
     $ 'fm per',i5,' event',1x,'*')
 911  format(6x,'*',10x,
     $ 'Impact parametetrs are distributed uniformely:',9x,'*'
     $  /6x,'*',16x,f7.3,' < b < ',f7.3,' (fm)',23x,'*')
 912  format(6x,'*',8x,
     $ 'Impact parameters are distributed according to b^2:',6x,'*'
     $  /6x,'*',16x,f7.3,' < b < ',f7.3,' (fm)',23x,'*')
 900  format(6x,'*',65x,'*'/6x,'*',15x,
     & 'Beam energy   = ',f9.2,' A GeV',19x,'*'/
     &  6x,'*',15x,
     & 'Beam momentum = ',f9.2,' A GeV/c',17x,'*'/
     &  6x,'*',15x,
     & 'C.M. energy   = ',f9.2,' A GeV',19x,'*')
 901  format(6x,'*',65x,'*'/6x,'*',15x,
     & 'Beam energy   = ',f9.2,' A MeV',19x,'*'/
     &  6x,'*',15x,
     & 'Beam momentum = ',f9.2,' A MeV/c',17x,'*'/
     &  6x,'*',15x,
     & 'C.M. energy   = ',f9.2,' A GeV',19x,'*')
 902  format(6x,'*',31x,'Beam',8x,'Target',16x,'*')
 903  format(6x,'*',10x,'Velocity/c:',7x,f9.7,5x,f9.7,14x,'*')
 904  format(6x,'*',10x,'Gamma factor:',3x,f9.3,5x,f9.3,16x,'*')
 905  format(6x,'*',10x,'p_x(GeV/c):',5x,f9.3,5x,f9.3,16x,'*')
 906  format(6x,'*',10x,'p_z(GeV/c):',5x,f9.3,5x,f9.3,16x,'*')
 907  format(6x,'*',10x,'r_x(fm):',8x,f9.3,5x,f9.3,16x,'*')
 908  format(6x,'*',10x,'r_z(fm):',8x,f9.3,5x,f9.3,16x,'*')

 920  format(6x,'*',23x,'mass',i4,'(',i3,',',i3,')',25x,'*')
 921  format(6x,'*',65x,'*'/6x,'*',20x,
     & ' Box calculation with  ',22x,'*',/,
     & 6x,'*',65x,'*',/,6x,67('*'))
 922  format(6x,'*',65x,'*'/6x,'*',20x,
     &  'ground state properties',22x,'*',/,
     &   6x,'*',65x,'*',/,6x,67('*'))

 605  format(6x,'* Computational frame: Lab. frame')
 606  format(6x,'* Computational frame: C.M. frame')
 607  format(6x,'* Computational frame: nucl.-nucl. c.m.')
 608  format(6x,'* Computational frame: Collider')

 616  format(
     &       6x,'* Random seed     =',i10,/
     &      ,6x,'* # of event      =',i7,/
     &      ,6x,'* Time step size  =',f7.3,' fm/c'/
     &      ,6x,'* # of timestep   =',i7)

      end

c***********************************************************************

      subroutine jamerrm(merr,nerr,chmess)

c...Purpose: to inform user of errors in program execution.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      character chmess*(*)

c...Write first few "warnings", then be silent.
      if(merr.le.10) then
        mstd(25)=mstd(25)+1
        mstd(26)=merr
        if(mstc(13).ge.2.or.merr.eq.10.or.
     $         (mstc(13).gt.0.and.mstd(25).le.mstc(14))) then
          write(*,1000) merr,mstd(23),mstd(21),chmess
          do i=1,nerr
            write(*,'(5x,a)')check(i)
          end do
        endif

c...Write first few errors, then be silent or stop program.
      elseif(merr.le.20) then
        mstd(27)=mstd(27)+1
        mstd(28)=merr-10
        write(*,1100) merr-10,mstd(23),mstd(21),chmess
        do i=1,nerr
          write(*,'(5x,a)')check(i)
        end do
        if(mstc(11).eq.2) call parastop( 111 )

c...Stop program in case of irreparable error.
      else
        write(*,1400) merr-20,mstd(23),mstd(21),chmess
c...CPU time.
        call parastop( 111 )

      endif

c...Formats for output.
 1000 format('JAM warning type',i2,' given in timestep ',i4,
     &' of event ',i6,5x,a)
 1100 format(/5x,'JAM error type',i2,' has occured in timestep ',i4,
     &' of event ',i6,5x,a)
 1400 format(/5x,'irregular jet! Reaction type',i2,' has occured in time
     &step',i4,' of event ',i6,5x,a)

      end

c*********************************************************************

      subroutine jamlist(mlist)

c...Purpose: to list event record or particle data.

      include 'jam1.inc'
      include 'jam2.inc'
      character chau*16

      io=mstc(38)

      if(mlist.eq.1) then


      do i=1,nv
        if(k(1,i).le.10) then
          if(k(1,i).le.4) call pyname(k(2,i),chau)
          if(k(1,i).eq.5) call jamname(k(2,i),k(3,i),k(4,i),chau)
        else
        endif
      end do

      else if(mlist.eq.2) then

      do i=1,nv
        if(k(1,i).le.10) then
          if(k(1,i).le.4) call pyname(k(2,i),chau)
          if(k(1,i).eq.5) call jamname(k(2,i),k(3,i),k(4,i),chau)
        else
        endif
      end do

      else if(mlist.eq.3) then


      else if(mlist.eq.4) then



      else
      endif

 800  format(//,6x,'Event Listing of JAM records 1')
 801  format(//,6x,' LIST of JAM records 2')
 802  format(//,6x,' LIST of JAM records 3',i5,1x,i3,1x,i5/)
 803  format(//,6x,'Event Listing of JAM records 4')

 901  format('i k(1,) k(2,) k(10,) k(11,) r(4,) r(5,) v(4,) v(5,) p(5,)
     $ ')

 8000 format(i5,f10.3,1x,i4,1x,i10,1x,i10,1x,i9,1x,a7)
 917  format(i5,f10.3,1x,i4,1x,
     $  i7,1x,i6,1x,i9,1x,i8,1x,i9,1x,i5,1x,2i4,1x,a7)
 918  format(i5,1x,i3,i7,i6,1x,i6,1x,2(d9.4,1x,d10.4),1x,f7.4)
 919  format(i5,f7.4,1x,i3,i7,2i5,i8,1x,i8,1x,i8,1x,2i4,2i7,1x,a7)
 907  format(3d13.5,3f8.4,f7.4,i6,i6,i6,20i6)
 927  format(3d13.5,f8.4,1x,d13.5,1x,5f8.4,4f7.4,d11.5,i3,i6,20i5)
 911  format(i3,1x,i5,1x,8(i3,1x),i5,1x,3d13.5)
 912  format(1x,a9,i3,1x,i5,1x,2(i3,1x),i7,1x,3d13.5)

      end

c***********************************************************************

      subroutine jamname(kf1,kf2,kf3,chau)

c...Purpose: to give the particle/parton/nucleus name
c...as character string.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      character chau*16,chap*4,chamn*3,chaml*3,chary(7)*2
      character element(108)*4,num(0:9)*1
      dimension nl(7)
      data element/
     $ 'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ',
     $ 'F  ','Ne ','Na ','Mg ','Al ','Si ','P  ','S  ',
     $ 'Cl ','Ar ','K  ','Ca ','Sc ','Ti ','V  ','Cr ',
     $ 'Mn ','Fe ','Co ','Ni ','Cu ','Zn ','Ga ','Ge ',
     $ 'As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ',
     $ 'Nb ','Mo ','Te ','Ru ','Rh ','Pd ','Ag ','Cd ',
     $ 'In ','Sn ','Sb ','Te ',' J ','Xe ','Cs ','Ba ',
     $ 'La ','Ce ','Pr ','Nd ','Pm ','Sm ','Eu ','Gd ',
     $ 'Tb ','Dy ','Ho ','Er ','Tm ','Yb ','Lu ','Hf ',
     $ 'Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     $ 'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ',
     $ 'Ac ','Th ','Pa ','U  ','Np ','Pu ','Am ','Cm ',
     $ 'Bk ','Cf ','Es ','Fm ','Md ','No ','Lw ','Ku ',
     $ '105 ','106 ','107 ','108 '/
      data num/'0','1','2','3','4','5','6','7','8','9'/
      data chary /'L ','S-','S0','S-','X-','X0','Om'/

      chau=' '

c...Hadrons/partons.
      if(abs(kf1).lt.100000) then
        kc=jamcomp(kf1)
        if(kc.ne.0) chau=chaf(kc,(3-isign(1,kf1))/2)
        return
! y.sakaki 2023/07/19, User defined particles -->
      elseif( 900000 .le. abs(kf1) .and. abs(kf1) .le. 999999 ) then
        if(kf1 .gt. 0) write(chau,'(i6)') kf1
        if(kf1 .lt. 0) write(chau,'(i7)') kf1
        return
      endif

c...Nuclear cluster case.
      kfa1=abs(kf1)
      kfa2=abs(kf2)
      kfa3=abs(kf3)
      nn=mod(kfa1/1000000,1000)
      nz=mod(kfa1/1000,1000)
      nl(1)=mod(kfa1,1000)           ! # of lambda
      nl(2)=mod(kfa2/1000000,1000)   ! # of sigma-
      nl(3)=mod(kfa2/1000,1000)      ! # of sigma0
      nl(4)=mod(kfa2,1000)           ! # of sigma+
      nl(5)=mod(kfa3/1000000,1000)   ! # of xi-
      nl(6)=mod(kfa3/1000,1000)      ! # of xi0
      nl(7)=mod(kfa3,1000)           ! # of Omega

      ml=0
      do i=1,7
        ml=ml+nl(i)
      enddo
      m=nn+nz+ml
      if(m.le.1) then
        write(check(1),'(i10,1x,i10,1x,i10)')kf1,kf2,kf3
        call jamerrm(1,1,'(jamname:)funny m=1')
        chau='funny'
        return
      endif

      chap='  '
      if(nz.ge.1.and.nz.le.108) chap=element(nz)
      if(nz.eq.0) chap='n '

        i1=mod(m,10)
        i10=mod(m/10,10)
        i100=mod(m/100,10)
        if(m.le.9) then
          chamn=num(i1)
          chau=chamn(1:1)//chap
        else if(m.le.99) then
          chamn=num(i10)//num(i1)
          chau=chamn(1:2)//chap
        else
          chamn=num(i100)//num(i10)//num(i1)
          chau=chamn(1:3)//chap
        endif

        if(ml.eq.0) return

      do i=1,7
        len=16
        do j=1,16
          if(chau(j:j).eq.' ') then
             len=j-1
             goto 10
          endif
        end do
10      continue
        i1=mod(nl(i),10)
        i10=mod(nl(i)/10,10)
        i100=mod(nl(i)/100,10)
        if(nl(i).le.9) then
          chaml=num(i1)
          chau=chau(1:len)//chaml(1:1)//chary(i)
        else if(nl(i).le.99) then
          chaml=num(i10)//num(i1)
          chau=chau(1:len)//chaml(1:2)//chary(i)
        else
          chaml=num(i100)//num(i10)//num(i1)
          chau=chau(1:len)//chaml(1:3)//chary(i)
        endif
      end do

      end

c***********************************************************************

      subroutine jamdisp(irec,kdt)

c...Purpose: to display particles on terminal
c...Modify for box.

      include 'jam1.inc'
      include 'jam2.inc'

      parameter (iaxis=2,scale=1.0d0)
      dimension ncount(-15:15,-9:9)
c...Box:maru
      character chdisp*80
      save ifirst
      save etot0
!$OMP THREADPRIVATE(ifirst,etot0)
      data etot0/0.0d0/,ifirst/0/

      if(kdt.eq.1) then
        write(6,*)'event b=',mstd(21),pard(2),mstd(22)
      endif

      nbar=0
      nmes=0
      npart=0
      ekin=0.0d0
      do i=1,nv
        if(k(1,i).ge.1.and.k(1,i).le.10.or.k(1,i).lt.0) then
          ekin=ekin+p(4+mstc(190)*2,i)
          k9=k(9,i)
          if(abs(k9).eq.3) then
            if(r(4,i).le.pard(1)) nbar=nbar+1
          else if(k9.eq.0) then
            if(abs(k(2,i)).eq.21) then
             if(r(4,i).le.pard(1))  npart=npart+1
            else
             if(r(4,i).le.pard(1))  nmes=nmes+1
            endif
          else if(abs(k9).eq.1.or.k9.eq.2) then
            npart=npart+1
          endif
        endif
      end do
      etot=ekin/mstd(11)/mstc(5)

c...Mean field potential.
      if(mstc(6).ge.2) then
        call jambuue(ekin,epot,etot)
      endif

      if(kdt.eq.1) ifirst=0
      if(ifirst.eq.0) then
        etot0=etot
        ifirst=1
        if(mstc(8).le.1) return
      endif

      econ=abs(etot-etot0)/etot0*100
      ene0=abs(etot-etot0)*1000

      write(irec,810)pard(1),ene0,econ,
     $    nbar/mstc(5),nmes/mstc(5),npart/mstc(5)
810   format(f6.2,'(fm/c) econ(MeV/A)=',f8.4,' econ(%)=',f6.4,
     $  ' B',i3,' M',i8,' P',i3)
      if(mstd(2).eq.0.or.mstd(5).eq.0) then
         call jambuur(rms)
         write(irec,*)'B.E. rms',
     $ (ekin-parc(28))*1000,epot*1000,(etot-parc(28))*1000,rms
      endif


c...Display particles on terminal.
      if(mstc(16).ge.1) then

        do 10 i=-15,15
        do 10 j=-9,9
         ncount(i,j)=0
10      continue

        ix=mod(iaxis,3)+1
        iy=mod(ix   ,3)+1

        if(mstc(16).eq.1.or.mstc(16).eq.11) then ! all hadrons
          m1=1
          m2=nbary+nmeson
        else if(mstc(16).eq.2.or.mstc(16).eq.12)then ! only baryons
          m1=1
          m2=nbary
        else                  ! only mesons
          m1=nbary+1
          m2=nbary+nmeson
        endif

c...Box maru.
        rmax1=-1000.0d0
        rmin1=1000.0d0
        rmax2=-1000.0d0
        rmin2=1000.0d0
        drcell=pard(21)

        do 100 ii=m1,m2
          if(r(4,ii).gt.pard(1)) goto 100
          if(k(1,ii).gt.10) goto 100
          if (mstc(16).le.3) then
            dt=pard(1)-r(4,ii)
            rx=r(ix,ii)+p(ix,ii)/p(4+mstc(190)*2,ii)*dt
            ry=r(iy,ii)+p(iy,ii)/p(4+mstc(190)*2,ii)*dt
            i=nint(rx*scale/parc(6))
            j=nint(ry*scale/parc(6))
c...Box by maru
            if(mstc(4).eq.10) then
              rx=mod(mod(rx,drcell)+drcell*1.5d0,drcell)-drcell/2
              ry=mod(mod(ry,drcell)+drcell*1.5d0,drcell)-drcell/2
            endif

            i=nint(rx*scale/parc(6))
            j=nint(ry*scale/parc(6))

            if(rmax1.le.rx)rmax1=rx
            if(rmin1.ge.rx)rmin1=rx
            if(rmax2.le.ry)rmax2=ry
            if(rmin2.ge.ry)rmin2=ry
c...Box end
          else
            i=nint(p(ix,ii)*scale/parc(6))
            j=nint(p(iy,ii)*scale/parc(6))
          endif
          if (abs(i).le.15.and.abs(j).le.9) ncount(i,j)=ncount(i,j)+1
 100    continue

c...Box by maru
      icell1=-16
      icell2=16
      jcell1=-10
      jcell2=10

      if(mstc(4).eq.10) then
        icell1=-nint(pard(21)/2*scale/parc(6))
        icell2=nint(pard(21)/2*scale/parc(6))
        jcell1=-nint(pard(21)/2*scale/parc(6))
        jcell2=nint(pard(21)/2*scale/parc(6))
        write(*,*) rmin1,rmax1,rmin2,rmax2
      endif

      do j=-9,9
        write(chdisp,820)(ncount(i,j),i=-15,15)
        if(j.eq.jcell1-1 .or. j.eq.jcell2+1)
     &    write(chdisp,825) ('--',i=-15,15)
        if(icell1.ge.-15) then
          if(chdisp(2*(15+icell1)+1:2*(15+icell1)+1).eq.'-') then
            chdisp(2*(15+icell1)+1:2*(15+icell1)+1)='+'
          else
            chdisp(2*(15+icell1)+1:2*(15+icell1)+1)='|'
          endif
        endif
        if(icell2.le.15) then
          if(chdisp(2*(15+icell2)+4:2*(15+icell2)+4).eq.'-') then
            chdisp(2*(15+icell2)+4:2*(15+icell2)+4)='+'
          else
            chdisp(2*(15+icell2)+4:2*(15+icell2)+4)='|'
          endif
        endif
        write(irec,'(A)')chdisp
      enddo
c...Box end

c...Print particles.
      write(irec,830)

      endif

820   format(' ',31i2.0)
825   format(' ',31a2)
830   format(' ',/)

      end

c***********************************************************************

      subroutine jamboost

c...Purpose: to boost the ground state nuclei.
      include 'jam1.inc'
      include 'jam2.inc'

      do i=1,nv

c....Target.
        if(i.le.mstd(5)*mstc(5)) then
          rxb=pard(47)
          rzb=pard(49)
          pxb=pard(41)
          pzb=pard(43)
          gmb=pard(46)

c.....Projectail.
        else
          rxb=pard(37) !x
          rzb=pard(39) !z
          pxb=pard(31) !px
          pzb=pard(33) !pz
          gmb=pard(36) !gamma
        end if

        r(4,i)=0.0d0
        r(5,i)=0.0d0
        r(1,i)=r(1,i)+rxb
        r(3,i)=r(3,i)/gmb+rzb

        p(1,i)=p(1,i)+pxb
        p(3,i)=p(3,i)*gmb+pzb
        p(4,i)=sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
        p(6,i)=sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)

      end do

c...qmd:Calculate the relativistic distance and derivative


      end

c***********************************************************************

      subroutine jamgrund

c...Purpose: to make the goround state of target and projectile
c...Modify for box.

      include 'jam1.inc'
      include 'jam2.inc'

c...Set total particle number.
      nv=(mstd(2)+mstd(5))*mstc(5)
      nmeson=0
      nbary=nv

c...In the case of proj. mass number 1.
      if(mstd(2).eq.1) then
          kc1=jamcomp(mstd(1))
          if(kchg(kc1,6).eq.0) then
            nmeson=mstc(5)
            nbary=nv-mstc(5)
          endif
      endif

c...In the case of targ. mass number 1.
      if(mstd(5).eq.1) then
        kc2=jamcomp(mstd(4))
        if(kchg(kc2,6).eq.0) then
          nmeson=nmeson+mstc(5)
          nbary=nbary-mstc(5)
        endif
      endif

c...Min. distance between nucl.s
      dsam2=parc(13)**2
      ddif2=parc(14)**2

c....Loop over proj. and targ. in=1:targ. in=2:proj.
      do 1000 in=1,2

c...Target loop.
       if(in.eq.1) then
         if(mstc(6).le.-100) goto 1000
         nnm=mstd(5)
         nnp=mstd(6)
         if(nnm.eq.1) nnp=1
         iofset=0
c....Projectile loop.
       else if(in.eq.2) then
         nnm=mstd(2)
         nnp=mstd(3)
         if(nnm.eq.1) nnp=1
         iofset=mstd(5)*mstc(5)
         if(mstc(6).le.-100) iofset=0
       endif

c...Box by maru
       if(nnm.eq.0)goto 1000

c...Woods-Saxon R0
      if(mstc(4).ne.10) then
        rad0=1.19d0* nnm**0.333333d0-1.61d0* nnm**(-0.333333d0)
      endif

c...Loop over all particles.
      do ie=1,nnm*mstc(5)
        i=ie+iofset
        ksimul=(ie-1)/nnm+1
        iee=mod(ie,nnm)
        if(iee.eq.0) iee=mod(ie-1,nnm)+1
        itry=0
 10     continue
        itry=itry+1
        if(itry.ge.100) then
          write(6,*)'nnm nnp',nnm,nnp
          call jamerrm(30,0,'(jamgrund:)can not create g.s.')
        endif

c...Sample particle.
        if(mstc(4).ne.10) then
          rr=0.0d0
          if(nnm.gt.1) then

            if(in.eq.1) rr=hirnd(2)
            if(in.eq.2) rr=hirnd(1)


            if(mstc(157).eq.1) call jamgout(2,rr,in)
          endif
          cx=1.d0-2.d0*rn(0)
          sx=sqrt(1.d0-cx**2)
          phi=2*paru(1)*rn(0)
          r(1,i)=rr*sx*cos(phi)
          r(2,i)=rr*sx*sin(phi)
          r(3,i)=rr*cx

c...Box by maru: Distribute in a box or sphere.
        else
          r(1,i)=pard(21)*(0.5d0-rn(0))
          r(2,i)=pard(21)*(0.5d0-rn(0))
          r(3,i)=pard(21)*(0.5d0-rn(0))
        endif

        if(iee.le.nnp) then
          if(nnm.eq.1) then
           if(in.eq.2) then
             kf=mstd(1)
             kc=jamcomp(kf)
             p(5,i)=pard(34)
           else if(in.eq.1) then
             kf=mstd(4)
             kc=jamcomp(kf)
             p(5,i)=pard(44)
           endif
           k(2,i)=kf
           k(3,i)=0
           k(8,i)=ksimul
           k(9,i)=kchg(kc,6)*isign(1,kf)
           k(4,i)=0
          else
           k(2,i)=2212
           k(3,i)=0
           k(4,i)=0
           k(8,i)=ksimul
           k(9,i)=3
           p(5,i)=parc(25)
          endif
        else
         k(2,i)=2112
         k(3,i)=0
         k(4,i)=0
         k(8,i)=ksimul
         k(9,i)=3
         p(5,i)=parc(24)
        endif

c...Check if distance is large enough.
          do je=1,ie-1
            j=je+iofset
            if(k(8,i).ne.k(8,j)) goto 15
            rdis2=(r(1,i)-r(1,j))**2+(r(2,i)-r(2,j))**2
     $                                   +(r(3,i)-r(3,j))**2
            if(k(2,i).eq.k(2,j)) then
              if(rdis2.lt.dsam2) goto 10
            else
              if(rdis2.lt.ddif2) goto 10
            end if
 15       end do

        p(1,i)=0.0d0
        p(2,i)=0.0d0
        p(3,i)=0.0d0
        k(1,i)=1
        k(5,i)=-1
        k(6,i)=0
        if(in.eq.1) then
          k(7,i)=1
        else
          k(7,i)=-1
        endif
        r(4,i)= 0.0d0
        r(5,i)= 0.0d0
        k(10,i) = 0
        k(11,i) = 0
        v(1,i) = r(1,i)
        v(2,i) = r(2,i)
        v(3,i) = r(3,i)
        v(4,i) = r(4,i)
        v(5,i) =  1.d+35
        do j=1,10
        vq(j,i)=0.d0
        end do
        kq(1,i)=0
        kq(2,i)=0
c...Assign spin.

      end do

c...In the case of one particle.
      if(nnm.eq.1) then
        i=iofset+1
        p(1,i)=0.0d0
        p(2,i)=0.0d0
        p(3,i)=0.0d0
        p(4+mstc(190)*2,i)=sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
        goto 1001
      endif

c...C.M..correction
      cx=0.d0
      cy=0.d0
      cz=0.d0
      s=0.d0
      do 100 ie=1,nnm*mstc(5)
        i=iofset+ie
        cx=cx+r(1,i)*p(5,i)
        cy=cy+r(2,i)*p(5,i)
        cz=cz+r(3,i)*p(5,i)
        s=s+p(5,i)
 100  continue
      cx=-cx/s
      cy=-cy/s
      cz=-cz/s
      do 101 ie=1,nnm*mstc(5)
        i=iofset+ie
        r(1,i)=r(1,i)+cx
        r(2,i)=r(2,i)+cy
        r(3,i)=r(3,i)+cz
        v(1,i)=r(1,i)
        v(2,i)=r(2,i)
        v(3,i)=r(3,i)
 101  continue

c....Sampling (Fremi) momenta of nucleus.
      do 20 ie =1,nnm*mstc(5)
        i=iofset+ie
        itry=0
200     continue
        itry=itry+1
        if(itry.ge.3000) then
          call jamerrm(30,0,'(Infinit loop at init. momentm')
        endif

        if(mstc(4).ne.10) then
c...Calculate local Ferim momentum.
          if(mstc(43).ne.0.and.nnm.gt.2) then

            if(in.eq.1) then
              dens=wdsax2(sqrt(r(1,i)**2+r(2,i)**2+r(3,i)**2))*nnm
            else
              dens=wdsax1(sqrt(r(1,i)**2+r(2,i)**2+r(3,i)**2))*nnm
            endif
            pfermi=paru(3)*(1.5d0*paru(1)**2*dens)**(1.0d0/3.0d0)

c....Option for no Fermi motion. or Deuteron by K.NIITA
	  else
	    pfermi=0.0d0
	  endif

          pf=pfermi*rn(0)**(1.d0/3.d0)
          cth=1.d0-2.d0*rn(0)
          sth=sqrt(1.d0-cth**2)
          phi=paru(2)*rn(0)
          p(3,i)=pf*cth
          p(1,i)=pf*sth*cos(phi)
          p(2,i)=pf*sth*sin(phi)
c...Box
        else
         pf=sqrt(pard(22)**2+2*pard(22)*p(5,i))
         p(1,i)=pf*(2*rn(0)-1.0d0)
         p(2,i)=pf*(2*rn(0)-1.0d0)
         p(3,i)=pf*(2*rn(0)-1.0d0)
        endif

c....Check Pauli principle.

  20   continue

c...C.M..correction.
      cx=0.d0
      cy=0.d0
      cz=0.d0
      s=0.d0
      do 300 ie=1,nnm*mstc(5)
        i=iofset+ie
        cx=cx+p(1,i)
        cy=cy+p(2,i)
        cz=cz+p(3,i)
 300  continue
      cx=-cx/nnm
      cy=-cy/nnm
      cz=-cz/nnm
      do 301 ie=1,nnm*mstc(5)
        i=iofset+ie
        p(1,i)=p(1,i)+cx
        p(2,i)=p(2,i)+cy
        p(3,i)=p(3,i)+cz
        p(4+mstc(190)*2,i)=sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
 301  continue

1001  continue

1000  continue  !********* end loop over proj. and targ.

c...Box by maru: Adjust kinetic energy pard(22).
      if(mstc(4).eq.10) then
        itry=0
        efact=1.d0
5555    continue
        itry=itry+1
        if(itry.ge.100)
     $        call jamerrm(30,0,'(jamgrund:)box not converge')
        ekin=0.d0
        do i = 1,nv !{
          p(1,i)=p(1,i)*efact
          p(2,i)=p(2,i)*efact
          p(3,i)=p(3,i)*efact
          p(4+mstc(190)*2,i)
     $     =sqrt(p(5,i)**2+p(1,i)**2+p(2,i)**2+p(3,i)**2)
          ekin=ekin+p(4+mstc(190)*2,i)-p(5,i)
        enddo !}
        ekin=ekin/nv
        efact=pard(22)/ekin
        if(efact.gt.1.000001d0.or.efact.lt.0.999999d0) goto 5555
      endif
c...end box.


      end

c***********************************************************************

      subroutine jamjeti

c...Purpose: to initialize string decay routine jetset.

      include 'jam1.inc'
      include 'jam2.inc'
c...PYTHIA common block.
      common/pypars/mstp(200),parp(200),msti(200),pari(200)
      save /pypars/
!$OMP THREADPRIVATE(/pypars/)
c...HIJING common block.
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      common/hijdat/hidat0(10,10),hidat(10)
      save  /hiparnt/,/hijdat/
!$OMP THREADPRIVATE(/hijdat/)
!$OMP THREADPRIVATE(/hiparnt/)

c...Set up jetset
c---------------------

c...Number of lines available in the common bock LUJETS.
c...This should be changed if dimension of k(),p(),v() arrays are
c...changed.
      mstu(4)=1000

c...Titel page is not wrtten at first occasion.
      mstu(12)=0

c...Check on possible errors.
      mstu(21)=2

      if(mstc(8).eq.0) then  ! 0=no error message
        mstu(22)=0 ! (d=10) max. number of errors that are printed
        mstu(25)=0 ! (d=1) printing of warning
        mstu(26)=0 ! (d=10) max. number of warnings that are printed
      endif

c...Open a file to take the write out of the program
      mstu(11)=mstc(38)

c...Choice of mother pointers for the particles produced by
c...fragmenting parton system.

c...Forbid particle decay in pyexec.
      mstj(21)=0

C...Switch off decay of pi0, K0S, Lambda, Sigma+-, Xi0-, Omega-.
	mdcy(jamcomp(111),1)=0  ! pi0
	mdcy(jamcomp(311),1)=0  ! k0
	mdcy(jamcomp(-311),1)=0 ! ak0
	mdcy(jamcomp(411),1)=0    ! D+
	mdcy(jamcomp(-411),1)=0   ! D-
	mdcy(jamcomp(421),1)=0    ! D0
	mdcy(jamcomp(-421),1)=0   ! aD0
	mdcy(jamcomp(221),1)=0    ! eta
	mdcy(jamcomp(331),1)=0    ! eta'
	mdcy(jamcomp(441),1)=0    ! eta_c

        mdcy(jamcomp(310),1)=0
        mdcy(jamcomp(431),1)=0
        mdcy(jamcomp(-431),1)=0
        mdcy(jamcomp(511),1)=0
        mdcy(jamcomp(-511),1)=0
        mdcy(jamcomp(521),1)=0
        mdcy(jamcomp(-521),1)=0
        mdcy(jamcomp(531),1)=0
        mdcy(jamcomp(-531),1)=0
        mdcy(jamcomp(3122),1)=0
        mdcy(jamcomp(-3122),1)=0
        mdcy(jamcomp(3112),1)=0
        mdcy(jamcomp(-3112),1)=0
        mdcy(jamcomp(3212),1)=0
        mdcy(jamcomp(-3212),1)=0
        mdcy(jamcomp(3222),1)=0
        mdcy(jamcomp(-3222),1)=0
        mdcy(jamcomp(3312),1)=0
        mdcy(jamcomp(-3312),1)=0
        mdcy(jamcomp(3322),1)=0
        mdcy(jamcomp(-3322),1)=0
        mdcy(jamcomp(3334),1)=0
        mdcy(jamcomp(-3334),1)=0

c...Hadron formation point.
      mstj(10)=mstc(72)

c...Endpoint quarks obtain transverse momenta parj(21)
      mstj(13)=1

c...Width of the gaussian pt in string fragm. Lund default=0.36GeV.
      parj(21)=parc(66)  ! =hipr1(2)

c...Minimum kinetic energy in decays.
      parj(64)=parc(41)

c...Choice of baryon production model(1)
      mstj(12)=mstc(73)   ! =ihpr2(11)

      if(mstj(12).eq.3) then
c....P(qq)/P(q), the suppression of diquark-antidiquark pair prodcution.
        parj(1)=parj(1)*1.2d0
      else if(mstj(12).eq.4) then
        parj(1)=parj(1)*1.7d0
        parj(5)=parj(5)*1.2d0
      else if(mstj(12).eq.5) then
        parj(1)=parj(1)*2.0d0
        parj(18)=0.19d0
      endif

c...Prob. for production of light-mesons with spin 1.
c...Prob. for strange meson has spin 1.
c...Prob. for charm or heavier meson has spin 1.
c...Baryon resonance production parameter.

c...Reset ssbar production prob. for changed string tension.
      parj(1)=parj(1)**(1.0d0/parc(54))   ! P(qq/q)
      parj(2)=parj(2)**(1.0d0/parc(54))   ! P(s)/P(q)
      parj(3)=parj(3)**(1.0d0/parc(54))   ! (p(us)/p(ud))/(p(s)/p(d))
      parj(4)=1d0/3d0*(3d0*parj(4))**(1.0d0/parc(54)) !1/3p(ud1)/p(ud0)
      parj(21)=parj(21)*sqrt(parc(54)/1d0) ! width of Gaussian

c...a parameter of lund fragm. function.
      hipr1(3)=parj(41)

c...b parameter of lund fragm. function.
      hipr1(4)=parj(42)

c...Switch of hard scattering.
      ihpr2(8)=1   ! 1 is default value

c...Hard scattering switch off at low energy.
      if(pard(16).lt.parc(71)) then
        ihpr2(8)=0
        mstc(81)=0
      endif

c...Energy partitioning in hadron.
      mstp(92)=3
      parp(94)=1.5d0
      parp(96)=1.5d0

c...Initialize HIJING and PYTHIA.
      call jamhijin

      end

c***********************************************************************

      subroutine jamhijin

c...Purpose: to initialize HIJING for specified event type

      implicit double precision(a-h, o-z)
      include 'jam2.inc'

      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      common/hijdat/hidat0(10,10),hidat(10)
      save  /hiparnt/,/hijdat/
!$OMP THREADPRIVATE(/hijdat/)
!$OMP THREADPRIVATE(/hiparnt/)

      character frame*16,cproj*16,ctarg*16
      external fnkick,fnkickr,fnkick2,fnstru,fnstrum,fnstrus

      if(mstc(4).eq.0) then
	frame(1:3)='lab'
      else if(mstc(4).eq.1) then
	frame(1:3)='cm '
      else if(mstc(4).eq.2) then
	frame(1:3)='nn '
      else if(mstc(4).eq.3) then
	frame(1:3)='cld'
      endif

      if(mstc(8).ge.1) ihpr2(10)=1  ! print warning messages
      ihnt2(1)=mstd(2)
      ihnt2(2)=mstd(3)
      ihnt2(3)=mstd(5)
      ihnt2(4)=mstd(6)
      ihnt2(5)=2212
      ihnt2(6)=2212
      hint1(8)=parc(28) ! nucl. mass
      hint1(9)=parc(28)

      ihpr2(6)=mstc(84)  ! nuclear effect on the parton dist.

      imes1=0
      imes2=0

      if(mstd(2).eq.1) then
       kfproj=mstd(1)
       ihnt2(5)=kfproj

       kc=jamcomp(kfproj)
       if(kchg(kc,6)*isign(1,kfproj).eq.0) imes1=1

c...The rest mass of projetctile particle mass
        hint1(8)=pard(34)
        call pyname(ihnt2(5),cproj)
      else if(mstd(1).ge.1000000000) then
        call jamname(mstd(1),0,0,cproj)
      else
        cproj='non'
      endif

      if(mstd(5).eq.1) then
        ihnt2(6)=mstd(4)
        hint1(9)=pard(44)
        call pyname(ihnt2(6),ctarg)
        kc=jamcomp(mstd(4))
        if(kchg(kc,6)*isign(1,mstd(4)).eq.0) imes2=1

      else if(mstd(4).ge.1000000000) then
        call jamname(mstd(4),0,0,ctarg)
      else
        ctarg='non'
      endif


c...Set up Wood-Sax distributions.
      if(mstc(4).ne.10) then
        if(ihnt2(1).gt.1) then
          call hijwds(ihnt2(1),1,rmax)
          hipr1(34)=rmax
          pard(40)=rmax
        endif
c...Set up Wood-Sax distr for  targ.
        if(ihnt2(3).gt.1) then
          call hijwds(ihnt2(3),2,rmax)
          hipr1(35)=rmax
          pard(50)=rmax
        endif
      endif

cABE change @2014/07/24 to avoid error in soft collision
      if (mstd(1) .ne. 22) then
        if(ihnt2(1).eq.0.or.ihnt2(3).eq.0) return
      endif

      hint1(1)=pard(16) ! srt

c...Table fit
        i=0
20      i=i+1
        if(i.eq.10) go to 30
        if(hidat0(10,i).le.hint1(1)) go to 20
30      if(i.eq.1) i=2
        do 40 j=1,9
           hidat(j)=hidat0(j,i-1)+(hidat0(j,i)-hidat0(j,i-1))
     &     *(hint1(1)-hidat0(10,i-1))/(hidat0(10,i)-hidat0(10,i-1))
40      continue

      hipr1(31)=hidat(5)
      hipr1(30)=2.0d0*hidat(5)

c------JET -----------------------------------------------------------

c...Initialize pythia
c   ===================

c...Hard scattering for meson induced collision has not yet been implemented.
      if(imes1.eq.0.and.imes2.eq.0) then
       if((ihpr2(8).ne.0.or.ihpr2(3).ne.0)
     &                  .and.hint1(1).ge.parc(71)) then
        ih=mstc(37)
        call jampyin(0)

        if(mstc(8).ge.1)
     $   write(6,*)'Pythia initialization complete'

c...Calculate the jet cross sections involved with nucleon collisions.
        if(mstc(8).ge.1)write(6,*)'hijcrs calculating now'
        call hijcrs
        if(mstc(8).ge.1) then
        write(6,*)'hijcrs complete'

        write(ih,800)hint1(1)
        write(ih,805)hint1(13)
        write(ih,804)hint1(13)-hint1(12)
        write(ih,803)hint1(12)
        write(ih,801)hint1(10)
        write(ih,807)hipr1(31)
        write(ih,808)hint1(59)
        write(ih,810)
        write(ih,811)hint1(11),hint1(14),hint1(15),hint1(16),hint1(17)
        write(ih,806)(hint1(13)-hint1(12))/hint1(13)
        write(ih,812)hipr1(33)
        end if
       endif
      endif
c---------------------------------------------------------------------

c...Booking for some functions
      if(mstc(21).eq.1) then
        if(mstc(75).ne.0) then
          call hifun(3,0.0d0,36.0d0,fnkick) ! pt**2 for pt kick
        endif
        call hifun(8,0.0d0,36.0d0,fnkickr)  ! pt**2 for resonance.
      endif

c...Booking for x distribution of valence quarks

      if(mstc(8).ge.1) then
        write(mstc(37),100)frame(1:3),hint1(1)
     $   ,cproj(1:7),ihnt2(1),ihnt2(2) ,ctarg(1:7),ihnt2(3),ihnt2(4)
      endif

 800  format(/,5x,'Cross sections',
     $   ' in HIJING at c.m. energy(GeV) of',f8.3)
 801  format(5x,'h-h jet x-section    ',f8.3,'mb')
 803  format(5x,'h-h inel x-section   ',7x,f8.3,'mb')
 804  format(5x,'h-h elastic xsection ',5x,f8.3,'mb')
 805  format(5x,'h-h total x-section  ',6x,f8.3,'mb')
 806  format(5x,'elastic to total xsection',5x,f8.3)
 807  format(5x,'sigma_0              ',f8.3,'mb')
 808  format(5x,'h-h triggered jet x-section (mb)',f8.3)
 810  format(/,5x,'Jet xsection: <no shadowing part>',
     $  2x,'<proj. part>',2x,'<targ. part>',2x,'<cross term>')
 811  format(7x,2(f8.3,4x),7x,3(f8.3,6x))
 812  format(/,6x,'diffractive to inel ratio',f8.3)
 100  format(//10x,50('*')/
     &  10x,'*',48x,'*'/
     &  10x,'*',9x,'Jet xc has been initialized at',9x,'*'/
     &  10x,'*',13x,a3,'= ',f10.2,' GeV/A',14x,'*'/
     &  10x,'*',48x,'*'/
     &  10x,'*',5x,'for ',
     &  a7,'(',i3,',',i3,')',' + ',a7,'(',i3,',',i3,')',4x,'*'/
     &  10x,50('*'))

      end

c***********************************************************************

      subroutine jampyin(i_type)

c...Initialize PYTHIA for jet production.
c...i_type=0: all cuts are given through this subroutine
c...i_type>0: cuts and section of channels must be given before
c...the call to this subroutine, as each sepcified values

      implicit double precision(a-h, o-z)
      include 'jam2.inc'
c...HIJING commonblocks.
      common/hiparnt/hipr1(100),ihpr2(50),hint1(100),ihnt2(50)
      common/hipyint/mint4,mint5,atco(200,20),atxs(0:200)
      save  /hiparnt/,/hipyint/
!$OMP THREADPRIVATE(/hiparnt/)
!$OMP THREADPRIVATE(/hipyint/)
c...PYTHIA commonblocks.
      common/pyjets/njet,npad,kjet(1000,5),pjet(1000,6),vjet(1000,5)
      save /pyjets/
!$OMP THREADPRIVATE(/pyjets/)
      common/pysubs/msel,mselpd,msub(500),kfin(2,-40:40),ckin(200)
      common/pypars/mstp(200),parp(200),msti(200),pari(200)
      common/pyint1/mint(400),vint(400)
      common/pyint2/iset(500),kfpr(500,2),coef(500,20),icol(40,4,2)
      common/pyint5/ngenpd,ngen(0:500,3),xsec(0:500,3)
      save  /pysubs/,/pypars/,/pyint1/,/pyint2/,/pyint5/
!$OMP THREADPRIVATE(/pyint1/)
!$OMP THREADPRIVATE(/pyint2/)
!$OMP THREADPRIVATE(/pyint5/)
!$OMP THREADPRIVATE(/pypars/)
!$OMP THREADPRIVATE(/pysubs/)
c...Local arraies.
      dimension xsec0(0:5,0:500),coef0(0:5,500,20),ini(0:5),
     &          mint44(0:5),mint45(0:5)
      save coef0
!$OMP THREADPRIVATE(coef0)

      data ini/6*0/
      save ini !FURUTA
!$OMP THREADPRIVATE(ini)
      if(i_type.ne.0) go to 360

c...Lowest energy to do hard scattering.
      parp(2)=parc(71)
c...Second order running alpha_strong
      mstp(2)=2
c...Inclusion of K factor(d=0)
      mstp(33)=1
      parp(31)=hipr1(17)

c...Choice of nucleon structure function
      mstp(51)=mstc(82)
      ihpr2(7)=mstc(82)

c...Initial state radiation
      mstp(61)=1
      if(ihpr2(2).eq.0.or.ihpr2(2).eq.2) mstp(61)=0

c...Final state radiation
      mstp(71)=1
      if(ihpr2(2).eq.0.or.ihpr2(2).eq.1) mstp(71)=0

c...No multiple interaction
        mstp(81)=1

c...Structure of mutliple interaction
      mstp(82)=1

c...Fragmentation off (have to be done by local call)
      mstp(111)=0

c...No printout of initialization information
      if(ihpr2(10).eq.0) mstp(122)=0

c...Effective minimum transverse momentum p_Tmin for multiple
c...interactions with mstp(82)=1  (D=1.40GeV/c)
      parp(81)=hipr1(8)

      ckin(5)=hipr1(8)
      ckin(3)=hipr1(8)  ! minimum p_t transfer in (semi)hard scatt.
      ckin(4)=hipr1(9)
      if(hipr1(9).le.hipr1(8)) ckin(4)=-1.0d0

      ckin(9)=-10.0d0
      ckin(10)=10.0d0

c**** Switch on and off scattering channels

c...QCD subprocesses
      msel=0
      do isub=1,500
        msub(isub)=0
      end do
      msub(11)=1 ! f_i f_i(~) -> f_i f_i(~)
      msub(12)=1 ! f_i f_i~ -> f_k f_k~
      msub(13)=1 ! f_i f_i~ -> g g
      msub(28)=1 ! f_i g  -> f_i g
      msub(53)=1 ! g g -> f_k f_k~
      msub(68)=1 ! g g -> g g
      msub(81)=1 ! f_i f_i~ -> Q_i Q_i~
      msub(82)=1 ! g g -> Q_i Q_i~

c...Gluon decays are not allowed
      do j=1,min(8,mdcy(21,3))
        mdme(mdcy(21,2)+j-1,1)=0
      end do

      isel=4
c...Option to switch on B-quark production.
      if(hint1(1).ge.20.0d0 .and. ihpr2(18).eq.1) isel=5
      mdme(mdcy(21,2)+isel-1,1)=1

c...Direct photon production allowed
      msub(14)=1   ! q + qbar -> g     + gamma
      msub(18)=1   ! q + qbar -> gamma + gamma
      msub(29)=1   ! q + g    -> q     + gamma


      if(ini(0).ne.0) go to 800
      ihnt2(11)=1
      ihnt2(12)=1
      go to 400
360   if(ini(i_type).ne.0) go to 800
400   ini(i_type)=1
      if(ihpr2(10).eq.0) mstp(122)=0

        mstp(171)=1  ! switch for variable energies
        pfrm1=0.268d0
        pfrm2=0.268d0
        if(mstd(2).eq.1.or.mstc(4).eq.10)pfrm1=0.0d0
        if(mstd(5).eq.1.or.mstc(4).eq.10)pfrm2=0.0d0
        pjet(1,1)=pfrm1
        pjet(1,2)=pfrm1
        pjet(2,1)=-pfrm2
        pjet(2,2)=-pfrm2
        pjet(1,3)=pfrm1+pard(33)
        pjet(2,3)=-pfrm1+pard(43)
        pjet(1,4)=sqrt(pard(34)**2+pjet(1,3)**2)
        pjet(2,4)=sqrt(pard(44)**2+pjet(2,3)**2)
        pjet(1,5)=pard(34)
        pjet(2,5)=pard(44)
        call pyinit('five',ihnt2(5),ihnt2(6),hint1(1),icon)

c...Store the initialization information for late use
      mint4=mint(44)
      mint5=mint(45)
      mint44(i_type)=mint(44)
      mint45(i_type)=mint(45)
      atxs(0)=xsec(0,1)
      xsec0(i_type,0)=xsec(0,1)
      do 500 i=1,200
        atxs(i)=xsec(i,1)
        xsec0(i_type,i)=xsec(i,1)
        do 500 j=1,20
          atco(i,j)=coef(i,j)
          coef0(i_type,i,j)=coef(i,j)
500   continue

      return
800   mint(44)=mint44(i_type)
      mint(45)=mint45(i_type)
      mint4=mint(44)
      mint5=mint(45)
      xsec(0,1)=xsec0(i_type,0)
      atxs(0)=xsec(0,1)
      do 900 i=1,200
        xsec(i,1)=xsec0(i_type,i)
        atxs(i)=xsec(i,1)
      do 900 j=1,20
        coef(i,j)=coef0(i_type,i,j)
        atco(i,j)=coef(i,j)
900   continue

      end

c***********************************************************************

      function jamk(mc,i)

c...Give charge/strangeness of particle i.
      include 'jam1.inc'
      include 'jam2.inc'

      kf=k(2,i)

c...Charge.
      if(mc.eq.1) then
        if(kf.eq.92) then
          jamk=jamchge(kq(1,i))+jamchge(kq(2,i))
        else
          jamk=jamchge(kf)
        endif

c...Strangeness.
      else
        if(kf.eq.92) then
          jamk=kfprop(kq(1,i),3)+kfprop(kq(2,i),3)
        else if(k(1,i).eq.4) then
          jamk=kfprop(k(2,i),3)
        else
          kc=jamcomp(kf)
          jamk=kchg(kc,7)*isign(1,kf)
        endif
      endif

      end

c***********************************************************************

      function jamchge(kf)

c...Gives three times the charge for a particle/parton.

c...Double precision and integer declarations.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'

c...Read out charge and change sign for antiparticle.
      jamchge=0
      kc=jamcomp(kf)
      if(kc.ne.0) jamchge=kchg(kc,1)*isign(1,kf)

      end

c***********************************************************************
      function jamcomp(kf)

      implicit double precision(a-h, o-z)
      common/pydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
!$OMP THREADPRIVATE(/pydat2/)
      common/pydat4/chaf(500,2)
!$OMP THREADPRIVATE(/pydat4/)
      character chaf*16
      parameter(kfmx=35553)
      dimension kcode(kfmx)
      logical first
      save kcode
      save first
!$OMP THREADPRIVATE(kcode,first)
      data first/.true./

      if(first) then
        do i=1,kfmx
          kcode(i)=0
        end do
        first=.false.
        do kc=1,500
          kfa=kchg(kc,4)
          if(kfa.ge.1.and.kfa.le.kfmx) kcode(kfa)=kc
        end do
      endif

      if(kf.eq.0) then
        jamcomp=0
        return
      endif

      kfa=abs(kf)
      if(kfa.ge.1.and.kfa.le.kfmx) then
        if(mod(kfa/10,10).eq.0.and.kfa.lt.100000
     &     .and.mod(kfa/1000,10).gt.0) kfa=mod(kfa,10000)
        jamcomp=kcode(kfa)
      else
        jamcomp=jamcomp0(kf)
      endif


      end
c***********************************************************************

      function jamcomp0(kf)

c...Compress the standard KF codes for use in mass and decay arrays;
c...also checks whether a given code actually is defined.

      implicit double precision(a-h, o-z)
      common/pydat1/mstu(200),paru(200),mstj(200),parj(200)
      common/pydat2/kchg(500,7),pmas(500,4),parf(2000),vckm(4,4)
      save /pydat1/,/pydat2/
!$OMP THREADPRIVATE(/pydat1/)
!$OMP THREADPRIVATE(/pydat2/)

c...Local arrays and saved data.
      dimension kford(100:500),kcord(101:500)
      save kford,kcord,nford,kflast,kclast
!$OMP THREADPRIVATE(kford,kcord,nford,kflast,kclast)
c...Whenever necessary reorder codes for faster search.
      if(mstu(20).eq.0) then
        nford=100
        kford(100)=0
        do i=101,500
          kford(i)=0
          kcord(i)=0
        end do

        do 120 i=101,500
          kfa=kchg(i,4)
          if(kfa.le.100) goto 120
          nford=nford+1
          do 100 i1=nford-1,0,-1
            if(kfa.ge.kford(i1)) goto 110
            kford(i1+1)=kford(i1)
            kcord(i1+1)=kcord(i1)
  100     continue
  110     kford(i1+1)=kfa
          kcord(i1+1)=i
  120   continue
        mstu(20)=1
        kflast=0
        kclast=0
      endif

c...Fast action if same code as in latest call.
      if(kf.eq.kflast) then
        jamcomp0=kclast
        return
      endif

c...Starting values. Remove internal diquark flags.
      jamcomp0=0
      kfa=iabs(kf)
      if(mod(kfa/10,10).eq.0.and.kfa.lt.100000
     &     .and.mod(kfa/1000,10).gt.0) kfa=mod(kfa,10000)

c...Simple cases: direct translation.
      if(kfa.gt.kford(nford)) then
      elseif(kfa.le.100) then
        jamcomp0=kfa

c...Else binary search.
      else
        imin=100
        imax=nford+1
  130   iavg=(imin+imax)/2
        if(kford(iavg).gt.kfa) then
          imax=iavg
          if(imax.gt.imin+1) goto 130
        elseif(kford(iavg).lt.kfa) then
          imin=iavg
          if(imax.gt.imin+1) goto 130
        else
          jamcomp0=kcord(iavg)
        endif
      endif

c...Check if antiparticle allowed.
      if(jamcomp0.ne.0.and.kf.lt.0) then
        if(kchg(jamcomp0,3).eq.0) jamcomp0=0
      endif

c...Save codes for possible future fast action.
      kflast=kf
      kclast=jamcomp0

      end

c*********************************************************************

      subroutine jamsetpa

c...Setup particle kc code.
      implicit double precision(a-h, o-z)
      include 'jam2.inc'
      dimension ipar1(6),ipar2(6)

c...Loop over kc code.
      do i=1,6
       ipar1(i)=0
       ipar2(i)=0
      end do
      do kc=101,mstu(6)
       ibar=kchg(kc,6)
       id=kchg(kc,5)
       if(id.eq.0) goto 100
       if(id.eq.id_quark) goto 100
       if(id.eq.id_diq) goto 100
       if(id.eq.id_lept) goto 100
       if(id.eq.id_exc) goto 100
       if(id.eq.id_boson) goto 100
       if(id.eq.id_tec) goto 100
       if(id.eq.id_susy) goto 100
       if(id.eq.id_special) goto 100
       if(id.eq.id_mdiff) goto 100
       if(id.eq.id_bdiff) goto 100
       if(ibar.eq.0) then
          if(ipar1(1).eq.0)ipar1(1)=kc
          ipar2(1)=kc
       else if(ibar.eq.3) then
         if(id.eq.id_nucls) then
           if(ipar1(2).eq.0)ipar1(2)=kc
           ipar2(2)=kc
         else if(id.eq.id_delts) then
           if(ipar1(3).eq.0)ipar1(3)=kc
           ipar2(3)=kc
         else if(id.eq.id_lambs) then
           if(ipar1(4).eq.0)ipar1(4)=kc
           ipar2(4)=kc
         else if(id.eq.id_sigms) then
           if(ipar1(5).eq.0)ipar1(5)=kc
           ipar2(5)=kc
         else if(id.eq.id_xis) then
           if(ipar1(6).eq.0)ipar1(6)=kc
           ipar2(6)=kc
         endif
       endif
 100  end do



      mstc(22)=ipar1(2)
      mstc(23)=ipar2(2)-1
      mstc(24)=ipar1(3)
      mstc(25)=ipar2(3)-3
      mstc(26)=ipar1(4)
      mstc(27)=ipar2(4)
      mstc(28)=ipar1(5)
      mstc(29)=ipar2(5)-2
      mstc(30)=ipar1(6)
      mstc(31)=ipar2(6)-1
      mstc(32)=ipar1(1)
      mstc(33)=ipar2(1)

      mstc(20)=1

cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      end

c***********************************************************************

      subroutine jamupdat(mupda,lfn)

c...Facilitates the updating of particle and decay data
c...by allowing it to be done in an external file.
c...Double precision declarations.
      implicit double precision(a-h, o-z)
c...Commonblocks.
      include 'jam2.inc'
      common/pyint4/mwid(500),wids(500,5)
      save /pyint4/
!$OMP THREADPRIVATE(/pyint4/)

c...Local arrays, character variables and data.
      character chinl*130,chkf*9,chvar(26)*9,chlin*72,
     &chblk(20)*72,chold*16,chtmp*16,chnew*16,chcom*24
      dimension pmd(5),kmd(5)

      data chvar/ 'kchg(i,1)','kchg(i,2)','kchg(i,3)','kchg(i,4)',
     $'kchg(i,5)','kchg(i,6)','kchg(i,7)',
     &'pmas(i,1)','pmas(i,2)','pmas(i,3)','pmas(i,4)','mdcy(i,1)',
     &'mdcy(i,2)','mdcy(i,3)','mdme(i,1)','mdme(i,2)','mdme(i,3)',
     $'brat(i)  ',
     &'kfdp(i,1)','kfdp(i,2)','kfdp(i,3)','kfdp(i,4)','kfdp(i,5)',
     &'chaf(i,1)','chaf(i,2)','mwid(i)  '/

c...Write header if not yet done.
      if(mstu(12).ge.1) call pylist(0)

c...Write information on file for editing.
      if(mupda.eq.1) then
        do 110 kc=1,500
          write(lfn,9000) kchg(kc,4),(chaf(kc,j1),j1=1,2)
     &    ,(kchg(kc,j2),j2=1,3),(kchg(kc,j2),j2=5,7)
     $    ,(pmas(kc,j3),j3=1,4)
     &    ,mwid(kc),mdcy(kc,1)
          do 100 idc=mdcy(kc,2),mdcy(kc,2)+mdcy(kc,3)-1
            write(lfn,9100) (mdme(idc,j1),j1=1,3),brat(idc),
     $      (kfdp(idc,j3),j3=1,5)
  100     continue
  110   continue

c...Read complete set of information from edited file or
c...read partial set of new or updated information from edited file.
      elseif(mupda.eq.2.or.mupda.eq.3) then

c...Reset counters.
        kcc=100
        ndc=0
        chkf='         '
        if(mupda.eq.2) then
          do 120 i=1,mstu(6)
          do j=1,7
           kchg(i,j)=0
          end do
          do j=1,4
           pmas(i,j)=0.0d0
          end do
          chaf(i,1)=' '
          chaf(i,2)=' '
          mwid(i)=0
          mdcy(i,1)=0
          mdcy(i,2)=0
          mdcy(i,3)=0
  120     continue
        else
          do 130 kc=1,mstu(6)
            if(kc.gt.100.and.kchg(kc,4).gt.100) kcc=kc
            ndc=max(ndc,mdcy(kc,2)+mdcy(kc,3)-1)
  130     continue
        endif

c...Begin of loop: read new line; unknown whether particle or
c...decay data.
  140   read(lfn,5200,iostat=ios) chinl
        if( ios .eq. -1 ) goto 190

c...Identify particle code and whether already defined  (for MUPDA=3).
        if(chinl(2:10).ne.'         ') then
          chkf=chinl(2:10)
          read(chkf,5300) kf

          if(mupda.eq.2) then
            if(kf.le.100) then
              kc=kf
            else
              kcc=kcc+1
              kc=kcc
            endif
          else
            kcrep=0
            if(kf.le.100) then
              kcrep=kf
            else
              do 150 kcr=101,kcc
                if(kchg(kcr,4).eq.kf) kcrep=kcr
  150         continue
            endif
c...Remove duplicate old decay data.
            if(kcrep.ne.0) then
              idcrep=mdcy(kcrep,2)
              ndcrep=mdcy(kcrep,3)
              do 160 i=1,kcc
                if(mdcy(i,2).gt.idcrep) mdcy(i,2)=mdcy(i,2)-ndcrep
  160         continue
              do 180 i=idcrep,ndc-ndcrep
                mdme(i,1)=mdme(i+ndcrep,1)
                mdme(i,2)=mdme(i+ndcrep,2)
                mdme(i,3)=mdme(i+ndcrep,3)
                brat(i)=brat(i+ndcrep)
                do 170 j=1,5
                  kfdp(i,j)=kfdp(i+ndcrep,j)
  170           continue
  180         continue
              ndc=ndc-ndcrep
              kc=kcrep
            else
              kcc=kcc+1
              kc=kcc
            endif
          endif

c...Study line with particle data.
          if(kc.gt.mstu(6)) call pyerrm(27,
     &    '(JAMUPDAT:) Particle arrays full by KF ='//chkf)

          read(chinl,9000) kchg(kc,4),(chaf(kc,j1),j1=1,2)
     &    ,(kchg(kc,j2),j2=1,3),(kchg(kc,j3),j3=5,7)
     $    ,(pmas(kc,j4),j4=1,4)
     &    ,mwid(kc),mdcy(kc,1)

          mdcy(kc,2)=0
          mdcy(kc,3)=0

c...Study line with decay data.
        else
          ndc=ndc+1
          if(ndc.gt.mstu(7)) call pyerrm(27,
     &    '(JAMUPDA:) Decay data arrays full by KF ='//chkf)

          if(mdcy(kc,2).eq.0) mdcy(kc,2)=ndc
          mdcy(kc,3)=mdcy(kc,3)+1

          read(chinl,9100) (mdme(ndc,j1),j1=1,3),brat(ndc),
     $      (kfdp(ndc,j2),j2=1,5)

        endif

c...End of loop; ensure that JAMCOMP tables are updated.
        goto 140
  190   continue
        mstu(20)=0

c...Perform possible tests that new information is consistent.
        mstj24=mstj(24)
        mstj(24)=0
        mxdch=0  ! max. decay channel.
        kcmxd=0
        do 220 kc=1,mstu(6)
          kf=kchg(kc,4)
          if(kf.eq.0) goto 220
          write(chkf,5300) kf

          if(min(pmas(kc,1),pmas(kc,2),pmas(kc,3),pmas(kc,1)-pmas(kc,3),
     &    pmas(kc,4)).lt.0d0.or.mdcy(kc,3).lt.0) call pyerrm(17,
     &    '(PYUPDA:) Mass/width/life/(# channels) wrong for KF ='//chkf)

          if(kc.gt.100.and.kchg(kc,5).ne.id_susy
     $                     .and.mdcy(kc,3).gt.mxdch) then
             mxdch=mdcy(kc,3)
             kcmxd=kc
          endif
          brsum=0d0
          jdc=0
          do 210 idc=mdcy(kc,2),mdcy(kc,2)+mdcy(kc,3)-1
            if(mdme(idc,2).gt.80) goto 210
            kq=kchg(kc,1)
            pms=pmas(kc,1)-pmas(kc,3)-parj(64)
            merr=0
            jdc=jdc+1
            kdc=0
            do 200 j=1,5
              kp=kfdp(idc,j)
              if(kp.eq.0.or.kp.eq.81.or.iabs(kp).eq.82) then
                if(kp.eq.81) kq=0
              elseif(jamcomp(kp).eq.0) then
                write(6,*)'jamcomp=0',kp,jdc
                merr=3
              else
                kdc=kdc+1
                kq=kq-jamchge(kp)
                kpc=jamcomp(kp)
                pms=pms-pmas(kpc,1)+pmas(kpc,3)
                pmd(kdc)=pmas(kpc,1)-pmas(kpc,3)
                kmd(kdc)=kp
                if(mstj(24).gt.0) pms=pms+0.5d0*min(pmas(kpc,2),
     &          pmas(kpc,3))
              endif
  200       continue

            if(kq.ne.0) merr=max(2,merr)
            if(mwid(kc).eq.0.and.kf.ne.311.and.pms.lt.0d0)
     &      merr=max(1,merr)

            if(merr.eq.3) call pyerrm(17,
     &      '(JAMUPDAT:) Unknown particle code in decay of KF ='//chkf)
            if(merr.eq.2) call pyerrm(17,
     &      '(JAMUPDAT:) Charge not conserved in decay of KF ='//chkf)
            if(merr.eq.1) call pyerrm(7,
     &      '(JAMUPDAT:) Kinematically unallowed decay of KF ='//chkf)

            if(merr.eq.1) then
             io=mstc(38)
             write(io,*)'Kinema kc kf',kc,kf,pms,jdc,' ',chaf(kc,1)
             write(io,*)'p1 p2 parj64',pmas(kc,1),pmas(kc,3),parj(64)
             write(io,*)'pms0=',pmas(kc,1)-pmas(kc,3)-parj(64)
             do i=1,kdc
             write(io,*)'kfd',kmd(i),pmd(i),' ',chaf(jamcomp(kmd(i)),1)
             end do
             endif

            brsum=brsum+brat(idc)
  210     continue

          write(chtmp,5500) brsum
          if(abs(brsum).gt.0.0005d0.and.abs(brsum-1d0).gt.0.0005d0)
     &    call pyerrm(7,'(PYUPDA:) Sum of branching ratios is '//
     &    chtmp(9:16)//' for KF ='//chkf)
  220   continue
        mstj(24)=mstj24
        print *,'max. dec ch',mxdch,kcmxd,kchg(kcmxd,4),chaf(kcmxd,1)

c...Write DATA statements for inclusion in program.
      elseif(mupda.eq.4) then

c...Find out how many codes and decay channels are actually used.
        kcc=0
        ndc=0
        do 230 i=1,mstu(6)
          if(kchg(i,4).ne.0) then
            kcc=i
            ndc=max(ndc,mdcy(i,2)+mdcy(i,3)-1)
          endif
  230   continue

c...Initialize writing of DATA statements for inclusion in program.
        do 300 ivar=1,26

c...Write comment.
          if(ivar.eq.1) write(lfn,8000)
          if(ivar.eq.12) write(lfn,8100)
          if(ivar.eq.24) write(lfn,8200)
          if(ivar.eq.26) write(lfn,8300)

          ndim=mstu(6)
          if(ivar.ge.15.and.ivar.le.23) ndim=mstu(7)
          nlin=1
          chlin=' '
          chlin(7:35)='data ('//chvar(ivar)//',i=   1,    )/'
          llin=35
          chold='START'

c...Loop through variables for conversion to characters.
          do 280 idim=1,ndim
            if(ivar.eq.1) write(chtmp,5400) kchg(idim,1)
            if(ivar.eq.2) write(chtmp,5400) kchg(idim,2)
            if(ivar.eq.3) write(chtmp,5400) kchg(idim,3)
            if(ivar.eq.4) write(chtmp,5400) kchg(idim,4)
            if(ivar.eq.5) write(chtmp,5400) kchg(idim,5)
            if(ivar.eq.6) write(chtmp,5400) kchg(idim,6)
            if(ivar.eq.7) write(chtmp,5400) kchg(idim,7)

            if(ivar.eq.8) write(chtmp,5500) pmas(idim,1)
            if(ivar.eq.9) write(chtmp,5500) pmas(idim,2)
            if(ivar.eq.10) write(chtmp,5500) pmas(idim,3)
            if(ivar.eq.11) write(chtmp,5500) pmas(idim,4)

            if(ivar.eq.12) write(chtmp,5400) mdcy(idim,1)
            if(ivar.eq.13) write(chtmp,5400) mdcy(idim,2)
            if(ivar.eq.14) write(chtmp,5400) mdcy(idim,3)

            if(ivar.eq.15) write(chtmp,5400) mdme(idim,1)
            if(ivar.eq.16) write(chtmp,5400) mdme(idim,2)
            if(ivar.eq.17) write(chtmp,5400) mdme(idim,3)

            if(ivar.eq.18) write(chtmp,5600) brat(idim)

            if(ivar.eq.19) write(chtmp,5400) kfdp(idim,1)
            if(ivar.eq.20) write(chtmp,5400) kfdp(idim,2)
            if(ivar.eq.21) write(chtmp,5400) kfdp(idim,3)
            if(ivar.eq.22) write(chtmp,5400) kfdp(idim,4)
            if(ivar.eq.23) write(chtmp,5400) kfdp(idim,5)

            if(ivar.eq.24) chtmp=chaf(idim,1)
            if(ivar.eq.25) chtmp=chaf(idim,2)
            if(ivar.eq.26) write(chtmp,5400) mwid(idim)

c...Replace variables beyond what is properly defined.
            if(ivar.le.7) then
              if(idim.gt.kcc) chtmp='               0'
            elseif(ivar.le.11) then
              if(idim.gt.kcc) chtmp='             0.0'
            elseif(ivar.le.14) then
              if(idim.gt.kcc) chtmp='               0'
            elseif(ivar.le.17) then
              if(idim.gt.ndc) chtmp='               0'
            elseif(ivar.le.18) then
              if(idim.gt.ndc) chtmp='             0.0'
            elseif(ivar.le.23) then
              if(idim.gt.ndc) chtmp='               0'
            elseif(ivar.le.25) then
              if(idim.gt.kcc) chtmp='                '
            else
              if(idim.gt.kcc) chtmp='               0'
            endif

c...Length of variable, trailing decimal zeros, quotation marks.
            llow=1
            lhig=1
            do 240 ll=1,16
              if(chtmp(17-ll:17-ll).ne.' ') llow=17-ll
              if(chtmp(ll:ll).ne.' ') lhig=ll
  240       continue
            chnew=chtmp(llow:lhig)//' '
            lnew=1+lhig-llow
c......Real values.
            if((ivar.ge.8.and.ivar.le.11).or.ivar.eq.18) then
              lnew=lnew+1
  250         lnew=lnew-1
              if(lnew.ge.2.and.chnew(lnew:lnew).eq.'0') goto 250
              if(chnew(lnew:lnew).eq.'.') lnew=lnew-1
              if(lnew.eq.0) then
                chnew(1:3)='0d0'
                lnew=3
              else
                chnew(lnew+1:lnew+2)='d0'
                lnew=lnew+2
              endif
c....Characters.
            elseif(ivar.eq.24.or.ivar.eq.25) then
              do 260 ll=lnew,1,-1
                if(chnew(ll:ll).eq.'''') then
                  chtmp=chnew
                  chnew=chtmp(1:ll)//''''//chtmp(ll+1:11)
                  lnew=lnew+1
                endif
  260         continue
              lnew=min(14,lnew)
              chtmp=chnew
              chnew(1:lnew+2)=''''//chtmp(1:lnew)//''''
              lnew=lnew+2
            endif

c...Form composite character string, often including repetition counter.
            if(chnew.ne.chold) then
              nrpt=1
              chold=chnew
              chcom=chnew
              lcom=lnew
            else
              lrpt=lnew+1
              if(nrpt.ge.2) lrpt=lnew+3
              if(nrpt.ge.10) lrpt=lnew+4
              if(nrpt.ge.100) lrpt=lnew+5
              if(nrpt.ge.1000) lrpt=lnew+6
              llin=llin-lrpt
              nrpt=nrpt+1
              write(chtmp,5400) nrpt
              lrpt=1
              if(nrpt.ge.10) lrpt=2
              if(nrpt.ge.100) lrpt=3
              if(nrpt.ge.1000) lrpt=4
              chcom(1:lrpt+1+lnew)=chtmp(17-lrpt:16)//'*'//chnew(1:lnew)
              lcom=lrpt+1+lnew
            endif

c...Add characters to end of line, to new line (after storing old line),
c...or to new block of lines (after writing old block).
            if(llin+lcom.le.70) then
              chlin(llin+1:llin+lcom+1)=chcom(1:lcom)//','
              llin=llin+lcom+1
            elseif(nlin.le.19) then
              chlin(llin+1:72)=' '
              chblk(nlin)=chlin
              nlin=nlin+1
              chlin(6:6+lcom+1)='&'//chcom(1:lcom)//','
              llin=6+lcom+1
            else
              chlin(llin:72)='/'//' '
              chblk(nlin)=chlin
              write(chtmp,5400) idim-nrpt
              chblk(1)(30:33)=chtmp(13:16)
              do 270 ilin=1,nlin
                write(lfn,5700) chblk(ilin)
  270         continue
              nlin=1
              chlin=' '
              chlin(7:35+lcom+1)='data ('//chvar(ivar)//
     &        ',i=    ,    )/'//chcom(1:lcom)//','
              write(chtmp,5400) idim-nrpt+1
              chlin(25:28)=chtmp(13:16)
              llin=35+lcom+1
            endif
  280     continue

c...Write final block of lines.
          chlin(llin:72)='/'//' '
          chblk(nlin)=chlin
          write(chtmp,5400) ndim
          chblk(1)(30:33)=chtmp(13:16)
          do 290 ilin=1,nlin
            write(lfn,5700) chblk(ilin)
  290     continue
  300   continue
      endif

c...Formats for reading and writing particle data.
 9000 format(1x,i9,2x,a16,2x,a16,3i3,1x,i3,1x,2i3,3f12.5,1p,e13.5,2i3)
 9100 format(10x,3i5,f12.6,5i10)
 5000 format(1x,i9,2x,a16,2x,a16,3i3,3f12.5,1p,e13.5,2i3)
 5010 format(1x,i3,1x,i9,2x,a16,2x,a16,3i3,3f12.5,1p,e13.5,2i3)
 5100 format(10x,2i5,f12.6,5i10)
 5200 format(a130)
 5300 format(i9)
 5400 format(i16)
 5500 format(f16.5)
 5600 format(f16.6)
 5700 format(a72)
 8000 format(/'C...PYDAT2,',
     $ ' with particle data and flavour treatment parameters.')
 8100 format(/'c...PYDAT3, with particle decay parameters and data.')
 8200 format(/'c...PYDAT4, with particle names (character strings).')
 8300 format(/'c...Treatment of resonances.')

      end
