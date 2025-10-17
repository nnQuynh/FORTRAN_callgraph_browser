************************************************************************
*                                                                      *
      subroutine gemset
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              initialize GEM for PHITS                                *
*                                                                      *
*        parameters :                                                  *
*                                                                      *
*        alev    = 0.0  ; The GCCI level density parameter is used.    *
*                = 1.0  ; The level density parameter is given by a=A/8*
*                > 1.0  ; The level density parameter.                 *
*                         is given by a = A / alev.                    *
*                                                                      *
*        rcal    = 0.0  ; The Dostrovsky and Matsuse parameter set     *
*                         are used for the inverse reaction            *
*                         cross sections.                              *
*                = 10.0 ; The simple parameter set is used for         *
*                         the inverse reaction cross section           *
*                         r0 is set to 1.5.                            *
*                < 10.0 ; The simple parameter set is used for         *
*                         the inverse reaction cross section           *
*                         r0 is set to rcal.                           *
*                                                                      *
*        ifis    ne 0   ; Original parameter in the Atchison model     *
*                         is used.                                     *
*                = 0    ; New parameter set is used.                   *
*                                                                      *
*        nimax   : maximum number of ejectiles                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision(a-h, o-z)

*-----------------------------------------------------------------------

      common /options/alev, rcal, ifis
!$OMP THREADPRIVATE(/options/)
      common /exiejn/ nimax
!$OMP THREADPRIVATE(/exiejn/)

*-----------------------------------------------------------------------
*     read parameters
*-----------------------------------------------------------------------

            alev = 0.0d0
            rcal = 0.0d0
            ifis = 0

*-----------------------------------------------------------------------

            nimax = 70 !original

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine gemexec(iz,in,ex,jj,px,py,pz,pt,et,rm,wt,ierr,ipos)
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to execute statistical particle decay / fission.        *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*        iz, in     : proton and neutron number of mother              *
*        ex         : excitation energy of mother (MeV)                *
*        jj         : spin of mother (hbar)                            *
*        px,py,pz   : momentum vector of mother (GeV)                  *
*        pt         : absolute value of momentum of mother (GeV)       *
*        et         : energy of mother (sqrt(p**2+m**2) GeV)           *
*        rm         : rest mass of mother (GeV)                        *
*        wt         : weight change                                    *
*        ierr       : error flag                                       *
*        ipos  : = 0, call from ovly12 and 13 ,=1 from tally           *
*                = 2, call from sctneut                                *
*                                                                      *
*                                                                      *
************************************************************************
      use mod_ompparallel !--- NS 2020.04 del THREADPRIVATE

      use smmmod
      use fission_mod, only : zfis,afis,ufis,er,betf,fisinh

      implicit doubleprecision(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      parameter ( amu = 931.494d0 )

*-----------------------------------------------------------------------

      common /clusts/ nclust, kclust(3,nnn)
!$OMP THREADPRIVATE(/clusts/)
      common /clustu/ lclust(0:8,nnn), sclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustu/)
      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)
      common /emode1/ nevhin, nlowlv, nefiss, ntwidt, mtprec,ipcnt(70)
!$OMP THREADPRIVATE(/emode1/)

      common /smmflg/ ismm, ifbm

*-----------------------------------------------------------------------

      common /mom/ erec,bet0(3),ekin,bet1(3)
!$OMP THREADPRIVATE(/mom/)
      dimension ir(2), ir1(2)

      common /sim/nocas
!$OMP THREADPRIVATE(/sim/)



      integer, save :: initgm = 0
!$OMP THREADPRIVATE(initgm)
*-----------------------------------------------------------------------
*     initialization
*-----------------------------------------------------------------------

         if( initgm .eq. 0 ) then

            initgm = initgm + 1

            call gemset

         end if

*-----------------------------------------------------------------------
*     initialization for gamlib move to setpar in read00.f
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*        SMM loop
*-----------------------------------------------------------------------

      IF(ismm .eq. 1 .and. ex / dble(iz + in) .gt. 2.d0 .and.
     &      in + iz .gt. 90) then

       call smmexec(iz + in, iz, px, py, pz, ex, multp) ! run smmexec (SMM mode)

      ELSEIF(ifbm .eq. 1 .and. iz + in .le. 18 .and. iz + in .ge. 5
     & .and. ex .ge. 1.d0) then

       call fbmexec(iz + in, iz, px, py, pz, ex, multp) ! run fbmexec (Fermi Breakup mode)

      ELSE

       allocate(pxmfrg(1), pymfrg(1), pzmfrg(1),
     &      ibarf(1), ichaf(1), exmfrg(1), frgrm(1))

       multp = 1
       ibarf(1)  = iz + in
       ichaf(1)  = iz
       exmfrg(1) = ex
       frgrm(1)  = rm * 1.d3
       pxmfrg(1) = px * 1.d3
       pymfrg(1) = py * 1.d3
       pzmfrg(1) = pz * 1.d3
      ENDIF

      do 100 igmrp = 1, multp
       iz = ichaf(igmrp)
       in = ibarf(igmrp) - ichaf(igmrp)
       ex = exmfrg(igmrp)
       rm = frgrm(igmrp) * 1.d-3
       px = pxmfrg(igmrp) * 1.d-3
       py = pymfrg(igmrp) * 1.d-3
       pz = pzmfrg(igmrp) * 1.d-3
       pt = dsqrt(px**2.d0 + py**2.d0 + pz**2.d0)
       et = dsqrt( pt ** 2.d0 + rm ** 2.d0 )

*-----------------------------------------------------------------------
*        initial nucleus and values
*-----------------------------------------------------------------------


               ierr = 0

               erec = ( et - rm ) * 1000.0

            if( pt .gt. 0.0d0 ) then

               bet0(1) = px / pt
               bet0(2) = py / pt
               bet0(3) = pz / pt

            else

               bet0(1) = 0.0
               bet0(2) = 0.0
               bet0(3) = 1.0

            end if

               iz0 = iz
               in0 = in
               ia0 = in + iz

               e1  = ex
               a0  = dble(ia0)
               z0  = dble(iz0)

               iflag = 0
               nocas = 0

               ipcnt = 0 ! 2024/9/30 remember ejectile multiplicity

*-----------------------------------------------------------------------
*        initialize fission fragment
*-----------------------------------------------------------------------

          fisinh = .false.

          do k = 1, 2

            afis(k) = 0.0d0
            zfis(k) = 0.0d0
            ufis(k) = 0.0d0
            er(k)   = 0.0d0

            do kk = 1, 3
              betf(k,kk) = 0.d0
            end do

          end do

*-----------------------------------------------------------------------
*     Start decay.
*        zero set for booking
*-----------------------------------------------------------------------

          IF(igmrp .eq. 1) then
               nclust = 0
               i1     = 0
          ENDIF

               ifssev = 0

            do i = 1, 4

               kdecay(i) = 0

            end do

*-----------------------------------------------------------------------

            imstp = 0

 2000    continue

            imstp = imstp + 1

            if( i1 + 2 .gt. nnn ) goto 5000

*-----------------------------------------------------------------------
*        for event generator mode
*-----------------------------------------------------------------------

               ntwidt = nevhin

            if( nevhin .ge. 2 .and. imstp .ge. nevhin ) then

               ntwidt = 1

            end if

*-----------------------------------------------------------------------
*        start evaporation calculation
*
*        iflag=1 : no more emission
*        iflag=2 : fission
*        iflag=3 : emission occured
*-----------------------------------------------------------------------


            call stdcay(a0,z0,e1,sp,iaf,izf,iflag,ipos) ! T.Sato 2018/03/06

            if( .not. fisinh ) goto 40

               ifssev = ifssev + 1

*-----------------------------------------------------------------------
*     fission fragment and ejectile from fission fragment
*-----------------------------------------------------------------------

         do 510 k = 1, 2

               a0   = afis(k)
               z0   = zfis(k)
               e1   = ufis(k)
               erec = er(k)

               do kk = 1, 3
                  bet0(kk) = betf(k,kk)
               end do

 3000    continue

         if( i1 + 2 .gt. nnn ) goto 5000

               e1old = e1
               aold  = a0
               zold  = z0
               sp = dble(jj)    ! 2020/3/26 Residue inherit mothers spin. Consider angular momentum shift during evaporation someday

               call stdcay(a0,z0,e1,sp,iaf,izf,iflag,ipos)  ! T.Sato 2018/03/06

            if( iflag .eq. 3 ) then

               i1 = i1 + 1
               erex = 0.0
               sp   = 0.d0
               call gemout(i1,izf,iaf,ekin,sp,bet1,erex,wt)


            else if( iflag .eq. 1 )then

               izr = nint(z0)
               iar = nint(a0)

               i1 = i1 + 1
               sp = 0.d0
               call gemout(i1,izr,iar,erec,sp,bet0,e1,wt)

               goto 510

            end if

               goto 3000

 510    continue

               goto 6000

*-----------------------------------------------------------------------
*     ejectile and residual nucleus without fission
*-----------------------------------------------------------------------

 40      if( iflag .eq. 3 ) then

               i1 = i1 + 1
               erex = 0.d0
               sp   = 0.d0
               call gemout(i1,izf,iaf,ekin,sp,bet1,erex,wt)

         else if( iflag .eq. 1 )then

               izr = nint(z0)
               iar = nint(a0)

               i1 = i1 + 1
               sp = dble(jj)    ! 2020/3/26 Residue inherit mothers spin. Consider angular momentum shift during evaporation someday
               call gemout(i1,izr,iar,erec,sp,bet0,e1,wt)

               goto 6000

         end if

*-----------------------------------------------------------------------

      goto 2000

*-----------------------------------------------------------------------
*     end of evaporation
*-----------------------------------------------------------------------

 6000    continue

            if( ifssev .gt. 0 )  kdecay(4) = 1

  100    continue

        deallocate(exmfrg, frgrm, ibarf, ichaf, pxmfrg, pymfrg, pzmfrg)

        return

*-----------------------------------------------------------------------
*           check dimension
*-----------------------------------------------------------------------

 5000    continue

               ErrCha = ''
               ErrID = 'L:397/R:gemexec/F:gem.f' !E00_001_001
               call ErrWrite(ErrID,ErrCha)

               write(*,*) ' **** Error at gemexec, too many products'
               write(*,*) ' ========================================='
               write(*,*) ' nnn  = ', nnn

               ierr = 1
               call parastop( 999 )

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine gemout(i1,iz,ia,ekin0,spin,bet,ex,wt)
*                                                                      *
*                                                                      *
*       booking of out going particles from GEM
*                                                                      *
*     input:                                                           *
*                                                                      *
*        i1         : number of ejectile                               *
*        iz, ia     : proton and mass number of ejectile               *
*        ekin0      : kinetic energy of ejectile (MeV)                 *
*        spin       : angular momentum of ejectile (hbar)              *
*        bet        : momentum unit vector                             *
*        ex         : excitation energy of ejectile (MeV)              *
*        wt         : weight change                                    *
*                                                                      *
*     output:                                                          *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*        nclust   : total number of out going particles and nuclei     *
*                                                                      *
*        kclust(3,nclust)                                              *
*                                                                      *
*                   kclust(1,i) = 101 : final output                   *
*                   kclust(2,i) = 0                                    *
*                   kclust(3,i) = 0                                    *
*                                                                      *
*        lclust(i,nclust)                                              *
*                                                                      *
*                i = 0, angular momentum                               *
*                  = 1, proton number                                  *
*                  = 2, neutron number                                 *
*                  = 3,                                                *
*                  = 4,                                                *
*                  = 5, charge                                         *
*                  = 6,                                                *
*                  = 7,                                                *
*                  = 8, isomer level (0: Ground, 1,2: 1st, 2nd isomer) *
*                                                                      *
*        sclust(i,nclust)                                              *
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
*        kdecay(4) = 0 : no fission                                    *
*                  = 1 : with fission                                  *
*                                                                      *
*                                                                      *
************************************************************************
      use NGSDATAMOD, only : bindeg

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      parameter ( rpmass = 938.27, rnmass = 939.58 )

*-----------------------------------------------------------------------
*     common for output of evaporation and fission
*-----------------------------------------------------------------------

      include 'param00.inc'

      common /clusts/ nclust, kclust(3,nnn)
!$OMP THREADPRIVATE(/clusts/)
      common /clustu/ lclust(0:8,nnn), sclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustu/)
      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)

      dimension bet(3)

*-----------------------------------------------------------------------
*     number of ejectiles
*-----------------------------------------------------------------------

            nclust = i1

*-----------------------------------------------------------------------

            ekin = ekin0 / 1000.0
            rmsp = (dble(iz) * rpmass + dble(ia-iz) * rnmass
     &       - bindeg(iz,ia-iz) + ex) / 1000.d0
            pabs = sqrt( ekin**2 + 2.d0 * rmsp * ekin )

*-----------------------------------------------------------------------

               kclust(1,nclust)  = 101
               kclust(2,nclust)  = 0
               kclust(3,nclust)  = 0

               lclust(0,nclust)  = spin
               lclust(1,nclust)  = iz
               lclust(2,nclust)  = ia - iz
               lclust(3,nclust)  = 0
               lclust(4,nclust)  = 0
               lclust(5,nclust)  = iz
               lclust(6,nclust)  = 0
               lclust(7,nclust)  = 0
               lclust(8,nclust)  = 0

               sclust(0,nclust)  = 0.0
               sclust(1,nclust)  = pabs * bet(1)
               sclust(2,nclust)  = pabs * bet(2)
               sclust(3,nclust)  = pabs * bet(3)
               sclust(4,nclust)  = ekin + rmsp
               sclust(5,nclust)  = rmsp
               sclust(6,nclust)  = ex
               sclust(7,nclust)  = ekin * 1000.
               sclust(8,nclust)  = wt
               sclust(9,nclust)  = 0.0
               sclust(10,nclust) = 0.0d0
               sclust(11,nclust) = 0.0d0
               sclust(12,nclust) = 0.0d0

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine stdcay(a,z,u,sp,iaf,izf,iflag,ipos)

C/////////////////////////////////////////////////////////////////////
C <STDCAY>
C  Calculate evaporation process and fission process
C    - determine if emission occurs
C    - determine if fission occurs
C        -if fission , call FIS
C        -if no fission, determine kinetic energy, etc. for emittor
C                        and caluculate recoil energy, etc.
C=====================================================================
C <Subroutine>
C   gamma    :  Decay width calculation
C    fis     :  Fission calculation
C selectE    :  Select kinetic energy in the CM system
C--------------------------------------------------------------------
C <Function>
C  fprob     : Calculate fission probability
C=====================================================================
C <Variables>
C     a    : mass of parent nuclei -> residual nuclei   (IN and OUT)
C     z    : charge of parent nuclei    -> residual     (IN and OUT)
C     u    : excited energy of parent nuclei-> residual (IN and OUT)
C   iflag  : 1 = no more emission                       (OUT)
C          : 2 = fission
C          : 3 = emission occur
C  fisinh  : true=fission occur, false=no fission       (OUT)
c   erec   : recoil energy in the lab system            (IN and OUT)
c   ekin   : kinetic energy in the lab system           (OUT)
c   bet0   : unit vector of recoil momentum             (IN and OUT)
c   bet1   : unit vector of momentum of emittor         (IN and OUT)
C/////////////////////////////////////////////////////////////////////
      use NGSDATAMOD, only : energm, bindeg

      use fission_mod
      implicit doubleprecision(a-h,o-z)

      parameter (amu=931.494d0)
      parameter (pi=3.1415926535898d0)

      common /exiejn/ nimax
!$OMP THREADPRIVATE(/exiejn/)

      common /std1/ r(70),s(70),sigma,rr(70)
!$OMP THREADPRIVATE(/std1/)
      common /ejectl/ omega(70),ifa(70),ifz(70)
      common /emitr/ gj(70),q(70),V(70),delta(70),smalla(70)
!$OMP THREADPRIVATE(/emitr/)

*-----------------------------------------------------------------------
*     common for gamlib
*-----------------------------------------------------------------------

      common /qparm/ ielas,icasc,iqstep,lvlopt,igamma
      common /emode1/ nevhin, nlowlv, nefiss, ntwidt, mtprec,ipcnt(70)
!$OMP THREADPRIVATE(/emode1/)

      parameter (len1 =41600)
      parameter (len2 = 1600)

      common/ big1   /izao(len2),ipt(len2),iref(105)
      common/ big2   /elo(len1),isp(len1)

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /pnint/  ipnint
      integer ipnint

      common /muflag/ imuinthit,imubrmhit,imuppdhit,imucaphit,imucapflag
!$OMP THREADPRIVATE(/muflag/)

      common /qmdscm/ qmdscm_h0, qmdscm_d, qmdscm_rcls, iqmdscm

      common /gem/    ngem

      double precision rdef(70)
      integer counter
*-----------------------------------------------------------------------

C  Initialization
      iaf=0
      izf=0
      ia=nint(a)
      iz=nint(z)

      if(    ia .gt. 1 .and. iz .eq. ia) then
       write(ErrCha,*)
     & 'Caution: GEM may fail. Z = A at beginning of GEM.'
       ErrID = 'L:637/R:stdcay/F:gem.f' !W00_004_001
       call ErrWrite(ErrID,ErrCha)
       elseif(ia .gt. 1 .and. iz .eq.  0) then
       write(ErrCha,*)
     & 'Caution: GEM may fail. Z = 0 at beginning of GEM.'
       ErrID = 'L:642/R:stdcay/F:gem.f' !W00_004_002
       call ErrWrite(ErrID,ErrCha)
      endif

C   Evaporation calculation starts
      uran=rn(0)

C  Calculate decay width

      call gammag(a,z,u,sp,uran,ipos)  ! T.Sato 2018/03/06

C     Modify the decay widths for photonuclear reaction, Noda 2013/8/19
      if ( (ipnint .ge. 1 .or. imuinthit == 1) .and. ityp == 14 ) then

C     Copy default decay widths
cABE add @2014/07/29, to avoid break
         if (sigma .gt. 0.0d0) then
            do counter=1,70
               rdef(counter) = r(counter)
            end do

            inn = ia-iz
            call suppressalpha(iz,inn,u)
         end if
      end if

      if (sigma.le.0.0d0) then
       iflag=1  ! no more emission
       return
      endif

C... Select ejectile
      uran=uran*sigma
      sum=0.d0

      If( z .gt. 60.d0 .and. .not. fisinh) then
       if(mtprec .eq. 18) then ! Allow nothing but fission

          r           = 0.d0
          r(nimax-1)  = sigma

       elseif( nefiss .gt. 0 .and. fprob(z,a,u) .gt. 0.d0) then

          r(nimax-1) = r(1)
          r(1)       = 0.d0

       elseif( ifiss .eq. 1 .or. z .lt. 89 ) then

          r(nimax-1) = r(1) * fprob(z,a,u)
          r(1)       = r(1) - r(nimax-1)

       else  ! Actinoid && ver 2


****** temporary treatment until fission model revision is done *****
          ifiss      = 1
          r(nimax-1) = r(1) * fprob(z,a,u) ! use previous version
          r(1)       = r(1) - r(nimax-1)
          ifiss      = 2
*********************************************************************

       endif
      endif

      do 20 j = 1, nimax

        k=j
        sum=sum+r(j)
        if (sum.ge.uran) go to 30

 20   continue

       iflag=1
       u = 0.0d0
       return

 30   continue
      jemiss=k

c   2002.4.15 by Y. Watanabe (Kyushu Univ.)
c   r(jemiss) = 0 : emission of a nucleus with jemiss does
c   If allowed, an error (divided by zero) occures in sele
c
      if( iqmdscm .eq. 1 .and. r(jemiss) .eq. 0.d0 ) then
         iflag = 1  ! no more emission
         return
      endif

      if(jemiss .eq. nimax) then
        iflag=1  ! go to gamma emission
        return
      elseif(jemiss .eq. nimax-1) then
        fisinh=.true.
        call fis(a,z,u)
        iflag=2                 !: fission
        return
      endif

      if(mtprec .eq. 18) then
       if(r(1) .gt. 0.d0 .or. .not.fisinh) then ! 2017/8/25 sometimes r(1) of (fissile + n compound) is 0
        jemiss = 1 ! 2016/03/01 Ogawa. e-mode=2 + fission -> n emission
       else
        iflag=1  ! no more emission
        return
       endif
      endif

C     Reset decay width
      if ( (ipnint .ge. 1 .or. imuinthit == 1) .and. ityp == 14 ) then
         do counter=1,70
            r(counter) = rdef(counter)
         end do
      end if

*-----------------------------------------------------------------------
C...set A & Z of the ejectile

      iaf=ifa(jemiss)
      izf=ifz(jemiss)
      iflag=3   ! Emission occured

*-----------------------------------------------------------------------
C...select kinetic energy in the CM system
*-----------------------------------------------------------------------

      pran = rn(0)

      pekin = pran * r(jemiss) * rr(jemiss)

        call selectE(jemiss,a,z,u,pekin,ekin)

        if(ekin.le.0.d0) then
       iflag=1
       u = 0.0d0
       return
        endif

*-----------------------------------------------------------------------
C... if residue has no bound excitation states (A<=5), deexcite to ground   2015/4/10   Ogawa
*-----------------------------------------------------------------------

      if(int(a) - iaf .le. 5 .and.
     & (u -q(jemiss) - ekin .le. bindeg(iz - izf, ia - iz -(iaf - izf))
     &  - bindeg(iz - izf, ia - iz -(iaf - izf) -1) .or.
     & u -q(jemiss) - ekin .le. bindeg(iz - izf, ia - iz -(iaf - izf))
     &  - bindeg(iz - izf -1, ia - iz -(iaf - izf) )) ) then
       ekin = u -q(jemiss)
      endif

*-----------------------------------------------------------------------
C... Velocity of Parent nucleus in the LAB system: vres*bet0

      am=a*amu+energm(iz,ia,1)
      pres=sqrt(erec**2+2.d0*am*erec)
      vres=pres/(Erec+am)

*-----------------------------------------------------------------------

C... residual nucleus, A, Z
      a=a-dble(iaf)
      z=z-dble(izf)

*-----------------------------------------------------------------------
*     check level near ground state
*     and adjust excitation energy and exit energy
*-----------------------------------------------------------------------

      if( nlowlv .eq. 1 ) then

                  e = u - q(jemiss) - ekin

         if( e .gt. 0.0d0 ) then

                  ja  = nint( a )
                  jz  = nint( z )

                  jza = 1000*jz+ja

*-----------------------------------------------------------------------
*           special for 10B + n -> 11B -> alpha + 7Li
*           first excited state / ground = 14.9
*-----------------------------------------------------------------------

            if( ja+iaf .eq. 11 .and. jz+izf .eq. 5 .and.
     &          jza .ne. 3007 .and. u - q(jemiss) .lt. 4.63 ) then

                    iflag=1
                    u = 0.0d0
                    return

            else
     &      if( ja+iaf .eq. 11 .and. jz+izf .eq. 5 .and.
     &          jza .eq. 3007 .and. u - q(jemiss) .lt. 4.63 ) then

               if( 15.9 * rn(0) .gt. 1.d0 ) then

                  e = 0.47761d0

               else

                  e = 0.0d0

               end if

*-----------------------------------------------------------------------
*           special for 3He + n -> 4He -> 3H + 1H
*           all to ground
*-----------------------------------------------------------------------

            else if( ja+iaf .eq. 4 .and. jz+izf .eq. 2 .and.
     &               jza .eq. 1003 .and. u - q(jemiss) .lt. 1.00 ) then

                 e = 0.0d0

*-----------------------------------------------------------------------

            else

                  k1 = iref(jz)
                  k2 = iref(jz+1)-1

               do k = k1, k2

                  if( jza .eq. izao(k) ) then

                     ik  = k
                     nup = ipt(ik+1)-1
                     nlo = ipt(ik)

                     do i = nlo + 1, nup

                        if( e .le. elo(i) ) then

                              e = elo(i-1)

                           goto 10

                        end if

                     end do

                  end if

               end do

   10          continue

            end if

         end if

                  ekin = u - q(jemiss) - e

                  if(ekin.le.0.d0) then
                     iflag=1
                     u = 0.0d0
                     return
                  end if

      end if

*-----------------------------------------------------------------------
C... excitation energy

      u = u - q(jemiss) - ekin

      if(u.le.epsilon(1.0))u=0.0d0
*-----------------------------------------------------------------------

C... Momentum in the CM system : pcmx, pcmy, pcmz
      amf=iaf*amu+energm(izf,iaf,1)
      amr=(ia-iaf)*amu+energm(iz-izf,ia-iaf,1)
      redm=amf*amr/am

                  if(redm.le.0.d0) then
                     iflag=1
                     u = 0.0d0
                     return
                  end if

      pcm=sqrt(ekin**2+2.d0*redm*ekin)
      th  = acos( 2.0d0 * rn(0) - 1.0d0 )
      ph  = 2.0d0 * pi  * rn(0)

      pcmx = pcm* dsin(th) * dcos(ph)
      pcmy = pcm* dsin(th) * dsin(ph)
      pcmz = pcm* dcos(th)

*-----------------------------------------------------------------------

C... Boost ejectile momentum to Lab system : p1x,p1y,p1z
      gam=sqrt(1.d0-vres**2)
      gam=1.d0/gam
      pv=pcmx*bet0(1)+pcmy*bet0(2)+pcmz*bet0(3)
      pv=pv*vres

      e1cm=sqrt(pcm**2+amf**2)
      tr=gam*(e1cm+gam*pv/(gam+1))

      p1x=pcmx+vres*bet0(1)*tr
      p1y=pcmy+vres*bet0(2)*tr
      p1z=pcmz+vres*bet0(3)*tr

C... Boost residual momentum to Lab system: p2x,p2y,p2z
      e2cm=sqrt(pcm**2+amr**2)
      tr=gam*(e2cm+gam*(pv)/(gam+1))

      p2x=-pcmx+vres*bet0(1)*tr
      p2y=-pcmy+vres*bet0(2)*tr
      p2z=-pcmz+vres*bet0(3)*tr

C... Kinetic Energy and velocity of Residual in Lab system: erec, bet0
      pr2=p2x**2+p2y**2+p2z**2
      erec= sqrt(amr**2+ pr2)-amr

      bet0(1)=p2x/sqrt(pr2)
      bet0(2)=p2y/sqrt(pr2)
      bet0(3)=p2z/sqrt(pr2)

C... Kinetic Energy and velocity of Emittor in Lab system: ekin, bet1
      pe2=p1x**2+p1y**2+p1z**2
      ekin= sqrt(amf**2+ pe2)-amf

      bet1(1)=p1x/sqrt(pe2)
      bet1(2)=p1y/sqrt(pe2)
      bet1(3)=p1z/sqrt(pe2)

      return
      end


************************************************************************
*                                                                      *
      subroutine gammag(a,z,u,sp,uran,ipos)
C/////////////////////////////////////////////////////////////////////
C GAMMA
C  calculate decay width for each particle emission
C=====================================================================
C <Subroutine>
C   eye    : Calculate I0,I1,I2,I3
C--------------------------------------------------------------------
C <Function>
C  dost    : Calculate kp, cp, or k_alpha
C  paire   : Calculate pairing energy
C=====================================================================
C <variables>
C     a    :   residual mass before emission                     (IN)
C     z    :   charge number of res nuclei before emission       (IN)
C     u    :   excitation energy of nuclei before emission       (IN)
C   uran   :   random number                                     (IN)
C     r    :   decay width                                      (OUT)
C     rr   :   decay width enchancement factor by excited state
C                            particle emission                  (OUT)
C   ifa    :   mass of emittor                                  (OUT)
C   ifz    :   charge of emittor                                (OUT)
C  omega   :   spin of emittor                                  (OUT)
C   gj     :   gj in eq.(##)                                    (OUT)
C   q      :   Q-value                                          (OUT)
C  delta   :   pairing energy                                   (OUT)
C  smalla  :   level density parameter                          (OUT)
C   V      :   Coulomb barrier                                  (OUT)
C   gamn   :   Decay width for neutron emission                 (OUT)
C   an     :   Level density parameter for neutron emission     (OUT)
C beta,alp :  Inverse cross section parameters for neutron emission(OUT)
C/////////////////////////////////////////////////////////////////////
      use levdenmod, only : rho_levden, getag
      use NGSDATAMOD, only : energm
      use levdat

      implicit doubleprecision(a-h,o-z)

      include 'param-physcnst.inc'
      common /exiejn/ nimax

      common /ejectl/omega(70),ifa(70),ifz(70)
      common /emitr/ gj(70),q(70),V(70),delta(70),smalla(70)
!$OMP THREADPRIVATE(/emitr/)
      common /std1/ r(70),s(70),sigma,rr(70)
!$OMP THREADPRIVATE(/std1/)
      common /fiss/ beta,alp
!$OMP THREADPRIVATE(/fiss/)
      common /options/alev, rcal, ifis
!$OMP THREADPRIVATE(/options/)

      common /sim/nocas
!$OMP THREADPRIVATE(/sim/)
      common /emode1/ nevhin, nlowlv, nefiss, ntwidt, mtprec,ipcnt(70)
!$OMP THREADPRIVATE(/emode1/)
      common /emode/  emodem, ge1, ge2, iemode
      logical lflg

      real*8 hbarc
      parameter (hbarc=197.327053d0)
      parameter (amu=931.494d0)
      parameter (pi=3.1415926535898d0)

      dimension couk(70),couc(70)
      data couk /70*1.d0/
      data couc /70*0.d0/
      save couk,couc !FURUTA
!$OMP THREADPRIVATE(couk,couc)
      common /gem/    ngem

      common /qmdscm/ qmdscm_h0, qmdscm_d, qmdscm_rcls, iqmdscm

* About ngem
*
*     0: ver 0
*     1: ver 1 (recommended)
*     2: ver 2
*     100: furihata v

*-----------------------------------------------------------------------
*        Parent nucleus
*-----------------------------------------------------------------------

               ia = nint(a)
               iz = nint(z)
               q1 = energm(iz,ia,min(2,ngem))

*-----------------------------------------------------------------------

      if(q1.eq.1.d10) then

         sigma = -1.0
         return
      endif

*-----------------------------------------------------------------------
*        Set cj and kj for Dostrovsky's parameter set
*-----------------------------------------------------------------------

            if( rcal .eq. 0.0 ) then

               couk(2) = dostg(1,z-ifz(2))

               couk(3) = couk(2) + 0.06
               couk(4) = couk(2) + 0.12

               couk(6) = dostg(2,z-ifz(6))
               couk(5) = couk(6) - 0.06


               couc(2) = dostg(3,z-ifz(2))

               couc(3) = couc(2) / 2.d0
               couc(4) = couc(2) / 3.d0

            end if

*-----------------------------------------------------------------------
*     Start calcualting Gamma for each particle emission
*-----------------------------------------------------------------------

      do 1 j = 1, nimax

*-----------------------------------------------------------------------
*        initialization
*-----------------------------------------------------------------------

             s(j)  = 0.0d0
             r(j)  = 0.0d0
             rr(j) = 1.0d0
             gj(j) = 0.0d0

             if( ifa(j) .eq. 0 ) goto 1

*-----------------------------------------------------------------------
*        special calculation for K.NIITA
*        without neutron width
*-----------------------------------------------------------------------
         if( ntwidt .eq. 1 .and. j .eq. 1 ) goto 1
         if( ntwidt .ge. 2 .and. j .gt. 1 ) goto 1

*-----------------------------------------------------------------------
*        Reduce calculation time by KN
*-----------------------------------------------------------------------

            if( j .gt. 6 ) then

               if( ia .gt. 40 .and.
     &             uran .lt. 0.95 ) goto 1

               if( ia .gt. 30 .and. ia .le. 40 .and.
     &             uran .lt. 0.93 ) goto 1

               if( ia .gt. 20 .and. ia .le. 30 .and.
     &             uran .lt. 0.7 ) goto 1

            end if

*-----------------------------------------------------------------------
*        daughter nucleus mass and charge
*-----------------------------------------------------------------------

            iaa = ia - ifa(j)
            aa  = dble(iaa)
            izz = iz - ifz(j)
            zz  = dble(izz)
            nn  = iaa - izz

            if(izz .eq. 0 .and. nn .gt. 1 ) goto 1 ! 2014/10/17 ogawa. No Dineutron, Trineutron...
            if(izz .gt. 1 .and. nn .eq. 0 ) goto 1 ! No Diproton, Triproton...    to avoid troubles in jqmdin

*-----------------------------------------------------------------------
*        Check of the residual nuclei after the emission
*        and avoid double counting modified by KN
*-----------------------------------------------------------------------
            if( ntwidt .eq. 1 .and. iaa .eq. 1 .and. nn .eq. 1 ) goto 1

*-----------------------------------------------------------------------

            if( iaa .le. 0 .or.
     &          izz .lt. 0 .or.
     &          iaa .lt. izz ) go to 1

            if( iaa .lt. ifa(j) .or.
     &          izz .lt. ifz(j) ) then

               do k = 1, j - 1

                  if( iaa .eq. ifa(k) .and.
     &                izz .eq. ifz(k) ) goto 1

               end do

             end if

*-----------------------------------------------------------------------
*        Q-value
*-----------------------------------------------------------------------

            q2   = energm(izz,iaa,min(2,ngem))
            q(j) = q2 - q1 + energm(ifz(j),ifa(j),min(2,ngem))

*-----------------------------------------------------------------------

          if(q2.eq.1.d10) then

         sigma = -1.0
         return
          endif

*-----------------------------------------------------------------------
*        Coulomb potential
*-----------------------------------------------------------------------

            V(j) = vcoul(zz,aa,ifz(j),ifa(j),couk(j),j)

*-----------------------------------------------------------------------
*        Reduce Coulomb potential in case emode = 2 and charged particle emission    2014/8/6  ogawa
*-----------------------------------------------------------------------

            if(iemode .ge. 2 .and. ipos .eq. 2 ) then
             selectcase(mtprec) ! fe is not (ein + q) if something other than neutrons is emitted
              case(0,4,16,17,19:21,37,38,50:91,102,152,153,160,161,
     & 875:891)  ! no charged particle emission
              case default
               if(u - q(j) - V(j) .le. 0.d0)  V(j) = (u - q(j)) * 0.9d0 ! apply only when excitation is insufficient
             endselect
            endif

*-----------------------------------------------------------------------
*        dependence of excitation energy : temp by KN
*-----------------------------------------------------------------------

            if( j .gt. 1 ) then

               if(ngem .le. 1) then  ! T.Sato 2018/12/11, ignore this adjustment in GEM2.0
                  V(j) = V(j) / ( 1.d0 + 0.005d0 * u / dble(ifa(j)) )
               endif

            end if

*-----------------------------------------------------------------------
*        Set Coulomb potential to zero if Q-value < 0,
*        for light residual e.g. 8Be, 9B
*-----------------------------------------------------------------------

            if( q(j) .le. 0 .and. aa .le. 20 ) V(j) = 0.d0

*-----------------------------------------------------------------------
*        Reduce Coulomb potential if neutron width is zero
*-----------------------------------------------------------------------

            if( ntwidt .eq. 1 ) then
               if( ia .le. 24 ) v(j) = 0.5 * v(j)
               if( ia .le. 13 ) v(j) = 0.d0
            end if

*-----------------------------------------------------------------------
*        Paring energy
*-----------------------------------------------------------------------

            delta(j) = paire(izz,nn)

*-----------------------------------------------------------------------
*        Check whether the emission is enegetically possible or not
*-----------------------------------------------------------------------

           if(ngem .ge. 2) then
            if(j .gt. 6)then !for nid
                if(u-q(j)-V(j).le.0.d0) goto 1 !original
            else
                if(u-q(j).le.0.d0) goto 1 !using statistic csinv
            endif
           else
            if( u - q(j) - V(j) .le. 0.d0 ) goto 1
           endif

*-----------------------------------------------------------------------
*        alpha and beta for neutrons: the precise parameter set
*-----------------------------------------------------------------------

            if( rcal .eq. 0.0 ) then

               alp  = 0.76d0 + 1.93d0 / aa**.333333333d0
               beta = ( 1.66d0 / aa**.666666666667d0 -5.d-2 ) / alp
               beta = max(beta,0.d0)

            else

               alp  = 1.d0
               beta = 0.d0

            end if

               bett = -V(j)
               if( j .eq. 1 ) bett = beta

*-----------------------------------------------------------------------
*        level density parameter and gamma
*-----------------------------------------------------------------------

            smalla(j) = getag(u-q(j)-V(j),izz,nn,isdum)

            if(ngem .le. 1) then
             call eye(j,aa,izz,nn,u,q(j),V(j),delta(j),smalla(j),bett,
     &               r(j),s(j))
            else
             call eye10(j,aa,izz,nn,u,q(j),V(j),delta(j),smalla(j),
     &               r(j),s(j))
            endif

*-----------------------------------------------------------------------
*        (2Sj+1)mj * alpha
*-----------------------------------------------------------------------

         if( iqmdscm .ne. 1 ) then
            gj(j) = ( 2.d0 * omega(j) + 1.d0 ) * dble(ifa(j))
         else
            gj(j) = ( 2.d0 * omega(j) + 1.d0 ) * dble(ifa(j))
     &             * aa / ( aa + dble(ifa(j)) )
         endif


        if(ngem .ge. 2) then
         if(j.gt.6)then ! for nid
            gj(j) = gj(j) * ( 1.d0 + couc(j) ) !ori
         endif !ori

        else
         if( j .eq. 1 ) then

            gj(j) = gj(j) * alp

         else

            gj(j) = gj(j) * ( 1.d0 + couc(j) )

         endif
        endif

*-----------------------------------------------------------------------
*        geometric cross section... sigma R
*-----------------------------------------------------------------------

      if((j .gt. 6 .and. ngem .ge. 2) .or. ngem .le. 1)then !for nid
            rmass = rb(aa,ifa(j),j)
            gj(j) = rmass * rmass * gj(j)
      endif

*-----------------------------------------------------------------------
*        decay width
*-----------------------------------------------------------------------

            r(j) = gj(j) * r(j)
            r(j) = max(r(j),0.d0)

*-----------------------------------------------------------------------
*       dw using nem
*-----------------------------------------------------------------------
            if(j.le.6 .and. ngem .ge. 2)then
                call nid(j,ia,iz,u,r(j))
                r(j)=r(j)
                r(j)=max(r(j),0.d0)
            endif

*-----------------------------------------------------------------------
*     decay width calculation from excited states
*-----------------------------------------------------------------------

            call levset(ifa(j), ifz(j), lflg)

            if( j .le. 6 .or.
     &          r(j) .eq. 0 .or.
     &          ubound(elevel,1) .eq. 0 ) then
              if(lflg) call levUNset
              cycle
            endif

               rrq = r(j)

*-----------------------------------------------------------------------
*        paring eneryg for a parent nucleus
*-----------------------------------------------------------------------

               del    = paire(iz,ia-iz)
               aparnt = getag(u,iz,ia-iz,isdum)

*-----------------------------------------------------------------------
*        level density of mother : rho_i
*-----------------------------------------------------------------------

               sp = -1.d0 ! Current GEM disregards spin. Therefore spin integrated level density.
               nv = 1
               rhop = rho_levden(ia,iz,u,sp,nv)

*-----------------------------------------------------------------------
*        sum up for excited states
*-----------------------------------------------------------------------

         do 10 i = 1, ubound(elevel,1)

*-----------------------------------------------------------------------
*           Q-value for an excited state
*-----------------------------------------------------------------------

               qq = q(j) + elevel(i)

               if( u - qq - V(j) .le. 0.d0 ) goto 10

*-----------------------------------------------------------------------
*           level denisty for an excited state
*-----------------------------------------------------------------------

               aq = getag(u-qq-V(j),izz,nn,isdum)

*-----------------------------------------------------------------------
*           integral part of gamma for an excited state
*-----------------------------------------------------------------------

               call eye(j,aa,izz,nn,u,qq,V(j),delta(j),aq,bett,rq,sq)

*-----------------------------------------------------------------------
*           (2Sj+1)mj * alpha * sigma_R for an excited state
*-----------------------------------------------------------------------

            if( iqmdscm .ne. 1 ) then
               gq = ( 2.d0 * spinpar(i,1) + 1.d0 ) * dble(ifa(j))
            else
               gq = ( 2.d0 * spinpar(i,1) + 1.d0 ) * dble(ifa(j))
     &             * aa / ( aa + dble(ifa(j)) )
            endif

               gq = rmass * rmass * gq

*-----------------------------------------------------------------------
*           gamma for an excited state [MeV]
*           reject an excited state
*           if the decay width [MeV] < level widht [MeV]
*-----------------------------------------------------------------------

               rrqg = ( gq * rq / rhop ) * amu / hbarc**2 / pi

               if( thalf(i) .gt. 0.d0 ) then

                  width = physc(3) / physc(2) * 1.d15 * log(2.d0)
     &                   / thalf(i)

                  if( width .ge. rrqg ) goto 20

               endif

*-----------------------------------------------------------------------

               rrq = rrq + gq * rq

  10     continue

*-----------------------------------------------------------------------
*           rr(j) : gamma(1-6) / total gamma
*-----------------------------------------------------------------------

  20     continue

         call levUNset

               rr(j) = r(j) / rrq
               r(j)  = rrq

*-----------------------------------------------------------------------
*           Reduce decay width if product and residual are identical
*-----------------------------------------------------------------------

               if(ifa(j) .eq. iaa .and. ifz(j) .eq. izz) then
                rr(j) = rr(j) / 2.d0
                 r(j) =  r(j) / 2.d0
               endif

*-----------------------------------------------------------------------

    1 continue
        sepene  = q(1) + 5.d0 !q1:n.bind q2:p.bind
        egamth = sepene

      if(u. le. egamth .and. u. ne. 0.d0 .and. sum(r(1:6)) .ne. 0.d0
     & .and. ngem .ge. 2) then
       sp = -1.d0
       nv = 2

       call gemgamcomp(ia, iz, u, sp, sumwid)
       r(nimax) = sumwid * 1.d21 * 4.136d-21 / 12.d0 * pi **1.5d0 *
     &  hbarc**2 /amu



      endif

*-----------------------------------------------------------------------
*        Reaction channel specific treatment of the advanced event generator mode.   2014/7/25 ogawa
*-----------------------------------------------------------------------
      if(iemode .ge. 2 .and. mtprec .ne. 0 .and. ipos .eq. 2 ) then

       selectcase(mtprec) ! MT > 120 are not yet supported.
        case(4,16,17,20,21,37,38,50:90,102,875:891,999)  ! gamma emission only.   2024/10/3 Do nothing for 91 because intent of 91 is unclear (only gamma or allow hadronic emission).
         r(1:nimax) = 0.d0
c        case(4,50:91,16,875:891,17,20,21,37)  ! neutron emission only. This must not happen
        case(28,41,42,103,600:649)  ! single proton emission
         r(1)       = 0.d0
         r(3:nimax) = 0.d0
         mtprec     = 999
        case(44,111)  ! 2 proton emission
         r(1)       = 0.d0
         r(3:nimax) = 0.d0
         ipcnt(2)   = ipcnt(2) + 1
         if( ipcnt(2) .gt. 2 ) mtprec = 999
        case(11,32,104,650:699) ! deuteron emission
         r(1:2)     = 0.d0
         r(4:nimax) = 0.d0
         mtprec     = 999
        case(33,105,700:749)  ! triton emission
         r(1:3)     = 0.d0
         r(5:nimax) = 0.d0
         mtprec     = 999
        case(34,106,750:799)  ! He3 emission
         r(1:4)     = 0.d0
         r(6:nimax) = 0.d0
         mtprec     = 999
        case(22,24,25,107,800:849)  ! one He4 emission
         r(1:5)     = 0.d0
         r(7:nimax) = 0.d0
         mtprec     = 999
        case(29,30,108)  ! 2 He4 emission
         r(1:5)     = 0.d0
         r(7:nimax) = 0.d0
         ipcnt(6)   = ipcnt(6) + 1
         if( ipcnt(6) .gt. 2 ) mtprec = 999
        case(23,109)  ! 3 He4 emission
         r(1:5)     = 0.d0
         r(7:nimax) = 0.d0
         ipcnt(6)   = ipcnt(6) + 1
         if( ipcnt(6) .gt. 3 ) mtprec = 999
        case(45,112)  ! p + a
         r(1)       = 0.d0
         r(3:5)     = 0.d0
         r(7:nimax) = 0.d0
        case(115)  ! p + d
         r(1)       = 0.d0
         r(4:nimax) = 0.d0
        case(116)  ! p + t
         r(1)       = 0.d0
         r(3)       = 0.d0
         r(5:nimax) = 0.d0
        case(35,114,117)  ! a + d
         r(1:2)     = 0.d0
         r(4:5)     = 0.d0
         r(7:nimax) = 0.d0
        case(36,113)  ! a + t
         r(1:3)     = 0.d0
         r(5)       = 0.d0
         r(7:nimax) = 0.d0
       endselect
      endif

*-----------------------------------------------------------------------
*        total gamma
*-----------------------------------------------------------------------

  900    continue

            sigma = 0.d0

         do j = 1, nimax

            sigma = sigma + r(j)

         end do


*-----------------------------------------------------------------------
c for debug cKN
*-----------------------------------------------------------------------

      goto 990

      if( sigma .gt. 0.0d0 ) then
      sigmaa=0.d0

      do j = 1, 6
        sigmaa=sigmaa+r(j)
      end do

         prob = (sigma-sigmaa)/sigma

      if( ia .gt. 50 .and.
     &    prob .gt. 0.05 )
     &    write(6,'(''50-> '',2i3,e13.5)') ia,iz, prob

      if( ia .gt. 40 .and. ia .le. 50 .and.
     &    prob .gt. 0.05 )
     &    write(6,'(''40-> '',2i3,e13.5)') ia,iz, prob

      if( ia .gt. 30 .and. ia .le. 40 .and.
     &    prob .gt. 0.07 )
     &    write(6,'(''30-> '',2i3,e13.5)') ia,iz, prob

      if( ia .gt. 20 .and. ia .le. 30 .and.
     &    prob .gt. 0.3 )
     &    write(6,'(''20-> '',2i3,e13.5)') ia,iz, prob


      end if

*-----------------------------------------------------------------------

  990 continue

      return
      end


************************************************************************
*                                                                      *
C--------------------------------------------------------------------
      subroutine selectE(j,a,z,u,pekin,ekin)
C--------------------------------------------------------------------
      use levdenmod, only : getag
      implicit doubleprecision(a-h,o-z)

      common /ejectl/ omega(70),ifa(70),ifz(70)
      common /emitr/ gj(70),q(70),V(70),delta(70),smalla(70)
!$OMP THREADPRIVATE(/emitr/)
      common /std1/ r(70),s(70),sigma,rr(70)
!$OMP THREADPRIVATE(/std1/)
      common /fiss/ beta,alp
!$OMP THREADPRIVATE(/fiss/)
      common /gem/    ngem
      common /sim/nocas
!$OMP THREADPRIVATE(/sim/)

      ia=nint(a)
      iz=nint(z)
      aa=a-dble(ifa(j))
      izz=nint(z)-ifz(j)
      nn=nint(aa)-izz
      iaa=nint(aa)

      if(iaa.eq.1) then
       ekin=u-q(j)
       return
      endif

      bet=beta
      if(j.ne.1) bet = -V(j)

      ux = 2.5d0 + 150.d0 / aa
      ex = ux + delta(j)
      ax = getag(ex,izz,nn,isdum)

      tau = sqrt(ax / ux) - 1.5d0 / ux
      tau = 1.d0 / tau
      sx = 2.d0 * sqrt( ax * ux )

      e0 = ex - tau * (dlog(tau) - .25* dlog(ax)
     & -1.25*dlog(ux) + sx) ! switching must be positive

      ppt=0.d0

      imax=1000

c.....Calculate excat decay width

      do 30 ii=1,2
        do 10 i=0,imax
*-----------------------------------------------------------------------
*for nid
            if (j.le.6 .and. ngem .ge. 2)then !for nid
                 x  = (u-q(j))/dble(imax)*dble(i)
                 pp = pde(j,ia,iaa,iz,izz,x,q(j),V(j),u)! n~4He
            else
                 x  = V(j)+(u-q(j)-V(j))/dble(imax)*dble(i) !ori
                 pp = pe(u,q(j),delta(j),V(j),bet,tau,e0,ex,smalla(j),x)!ori
            endif
*-----------------------------------------------------------------------

          if(i.eq.0)then
           ppo=pp
           goto 10
          endif

          dppt=(ppo+pp)*(u-q(j)-V(j))/dble(imax)/2.d0*gj(j)
          if(ngem .ge. 2) dppt=(ppo+pp)*(u-q(j))/dble(imax)/2.d0*gj(j)
          ppt=ppt+dppt
          ppo=pp
          ppm=ppt

        if(ppm.gt.pekin) goto 20

 10   continue
      ran=pekin/r(j)/rr(j)
      if(ii.eq.1) then
       pekin=ppm*ran
       goto 30
      else

        ekin = -1.0
        return
      endif

C.....dubeg write

 30   continue

 20   ekin=x
      if(a .gt. 270 .and. x .lt. -q(j)) ekin = -q(j) * 1.01d0 ! Avoid energy gain by emission
      return
      end

************************************************************************
*                                                                      *
C--------------------------------------------------------------------
      function pe(e,q,delta,V,bet,tau,e0,ex,smalla,x)
C--------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

      real(8) x1,x2 !FURUTA

      if(e-q.le.0) then
       pe=0
       return
      endif

      x1=e-q-ex                                      !FURUTA
      x2=e-q                                         !FURUTA
                                                     !FURUTA
      if( abs(x1-x).le.epsilon(1.0)*abs(x1) )then    !FURUTA x1=x
        eag = 2.d0*sqrt(smalla*(e-q-delta-x))        !FURUTA
        pe = exp( eag )                              !FURUTA
        pe= pe/smalla**.25d0                         !FURUTA
        pe= pe/(e-q-delta-x)**1.25d0                 !FRUUTA
      elseif( (x.gt.x1.and.x.le.x2) .or.             !FURUTA
     &     (abs(x2-x).le.epsilon(1.0)*abs(x2)) )then !FURUTA x2=x

        eag = (e-q-x-e0)/tau
        pe = exp( eag ) / tau

      else if( (x.ge.V.and.x.le.x1) .or.             !FURUTA
     &       abs(V-x).le.epsilon(1.0)*abs(V) )then   !FURUTA V=x

        eag = 2.d0*sqrt(smalla*(e-q-delta-x))
        pe = exp( eag )

        pe= pe/smalla**.25d0
        pe= pe/(e-q-delta-x)**1.25d0

      else
        pe=0.d0
      endif

      pe=pe*(x+bet)
      return
      end

************************************************************************
*                                                                      *
C--------------------------------------------------------------------
      subroutine eye(j,aa,izz,nn,u,q,V,delta,smalla,beta,r,s)
C
C     variable       IN/OUT
C     j              I       emittor identifier
C     aa             I       residual mass after emission (=A-Aj)
C     izz            I       charge number of res nuclei after emission
C     nn             I       neutron number of res nuclei after emission
C     u              I       excitation energy of nuclei before emission
C     q              I       Q-value
C     V              I       kV in the equation
C     delta          I       pairing energy
C     smalla         I       level density parameter at U-Q-delta-kV
C     bet            I       beta in the equation
C     r              O       r in the equation
C                                      _______________
C     s              O       s = 2 \/a(U-Q-delta-kV)

C--------------------------------------------------------------------
      use levdenmod, only : getag
      implicit doubleprecision(a-h,o-z)

*-----------------------------------------------------------------------

            if( u - q - V .le. 0.0 .or. aa .le. 0.0 ) then

               r = 0.d0
               s = 0.d0
               return

            end if

*-----------------------------------------------------------------------

               ux  = 2.5d0 + 150.d0 / aa
               ex  = ux + delta

               ax  = getag(ex,izz,nn,isdum)

               tau = sqrt( ax / ux ) - 1.5d0 / ux
               tau = 1.d0 / tau
               t = ( u - q  - V ) / tau
               sx = 2.d0 * sqrt( ax * ux )

               ssx = sx - ex / tau
               ssx = min( 200.d0, ssx )

               eest =  exp( ssx ) * tau / ax**.25d0 / ux**1.25d0

         if( u - q - V .le. ex ) then

               eept = exp( t )

               eye1 = ( eept - 1.0 - t ) * tau
               eye0 =   eept - 1.0

               r = ( eye1 + ( beta + V ) * eye0 ) * eest

         else

               s  = 2.d0 * sqrt( smalla * ( u - q - delta - V ) )
               tx = ex / tau

               eeps = exp( s )
               eepx = exp( tx )
               sexp = exp( sx - s )
               eye0 = ( eepx - 1.0 ) * eest
               eye1 = (  ( eepx * ( t  - tx + 1.0 )
     &                 - ( t + 1.0 ) ) * tau ) * eest
               eye2 = ey2(s,sx,sexp)        * eeps
               eye3 = ey3(s,sx,smalla,sexp) * eeps
               r =  eye3 + eye1  + ( beta + V ) * ( eye0 + eye2 )


         end if

      return
      end


C--------------------------------------------------------------------
      function ey2(s,sx,sexp)
C--------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

      common /qmdscm/ qmdscm_h0, qmdscm_d, qmdscm_rcls, iqmdscm

      func ( x ) =  1/x**1.5d0 + 1.5d0/x**2.5d0 + 3.75d0/x**3.5d0

      temp  = func(s)
      tempx = func(sx)

      if( iqmdscm .eq. 1 .and. sx .eq. 0.d0 ) tempx = 0.d0      ! S.Abe 2020/09/29

      ey2 = 2.d0 * sqrt(2.d0) *( temp - sexp * tempx )

      return
      end

C--------------------------------------------------------------------
      function ey3(s,sx,a,sexp)
C--------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

      ssqr = s**2
      sxqr = sx**2

      temp =     325.125d0 / s**4.5d0
     &     - sexp * ( ( 324.8d0   * ssqr + 3.28d0 * sxqr ) / sx**6.5d0 )
      temp = temp + 60.0d0 / s**3.5d0
     &     - sexp * ( ( 59.0625d0 * ssqr + 0.9375 * sxqr ) / sx**5.5d0 )
      temp = temp + 13.5d0 / s**2.5d0
     &     - sexp * ( ( 12.875d0  * ssqr + 0.625  * sxqr ) / sx**4.5d0 )
      temp = temp +  4.d0 / s**1.5d0
     &     - sexp * ( (  3.75d0   * ssqr + 0.25   * sxqr ) / sx**3.5d0 )
      temp = temp +  2.d0 / s**.5d0
     &     - sexp * ( (   1.5d0   * ssqr + 0.5    * sxqr ) / sx**2.5d0 )
      temp = temp
     &     - sexp * ( ( ssqr - sxqr ) / sx**1.5d0 )

      ey3 = temp / sqrt(2.d0) / a

      return
      end

C--------------------------------------------------------------------
      function ey0(t)
C--------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

      ey0 = exp(t)- 1.d0

      return
      end

C--------------------------------------------------------------------
      function ey1(t,tx,tau)
C--------------------------------------------------------------------
      implicit doubleprecision(a-h,o-z)

      ey1 = ( exp(tx) * ( t  - tx + 1.d0 ) - ( t + 1.d0 ) ) * tau

      return
      end

************************************************************************
*                                                                      *
      function dostg(i,z)

C=====================================================================
C     This routine is originally in HETC code
C=====================================================================

      implicit doubleprecision(a-h,o-z)

      dimension t(3,4)

*-----------------------------------------------------------------------
*     Set Dostrovsky's parameter set,
*     the footnote of PR116(1959)683 (p.699)
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------

C      i=1 -> Calculate kp
C      i=2 -> Calculate k_alpha
C      i=3 -> Calculate cp

*-----------------------------------------------------------------------
*        kp
*-----------------------------------------------------------------------

      data ( t(1,i), i = 1, 4 ) /
     &      0.51d0,    !  t(1,1)
     &      0.60d0,    !  t(1,2)
     &      0.66d0,    !  t(1,3)
     &      0.68d0/    !  t(1,4)

*-----------------------------------------------------------------------
*        k_alpha
*-----------------------------------------------------------------------

      data ( t(2,i), i = 1, 4 ) /
     &      0.81d0,    !  t(2,1)
     &      0.85d0,    !  t(2,2)
     &      0.89d0,    !  t(2,3)
     &      0.93d0/    !  t(2,4)

*-----------------------------------------------------------------------
*        cp
*-----------------------------------------------------------------------

      data ( t(3,i), i = 1, 4 ) /
     &       0.0d0 ,   !  t(3,1)
     &      -0.06d0,   !  t(3,2)
     &      -0.10d0,   !  t(3,3)
     &      -0.10d0/   !  t(3,4)

*-----------------------------------------------------------------------

      if (z-50.d0) 30,10,10
   10 dostg=t(i,4)
   20 return
   30 if (z-20.d0) 40,40,50
   40 dostg=t(i,1)
      go to 20
   50 n=.1d0*z
      x=10.d0*(n+1.d0)
      x=(x-z)*.1d0
      dostg=x*t(i,n-1)+(1.d0-x)*t(i,n)
      go to 20
      end


************************************************************************
*                                                                      *
C********This routine is the same as the drein1 in LAHET code**********
      subroutine drein1 (j,s,a,eye1,eye0)
      implicit doubleprecision(a-h,o-z)
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,
     2 dp2th=dp2/dp3)
c
c     compute statistical theory emission integrals
c     for s<dph use a series expansion.
c     for s>dph the explicit relationship
c     return for neutrons (and compute eye0)
c     return 1 for all others.
c
c     coeficients for series expansions
c
c     correction of c1 by r. e. prael
      dimension c0(7), c1(7)
      data c0 /0.66666667d0,0.25d0,0.06666667d0,0.01388889d0,
     1 0.00238095d0,0.00034722d0,0.00004409d0/
      data c1 /0.53333333d0,0.16666667d0,0.03809524d0,0.00694444d0,
     1 0.0010582d0,0.00013889d0,0.0000160d0/
c
      exps=dp0
      if (s.lt.1.d+02) exps=exp(-s)
      if (s.lt.dph) go to 10
c///// explicit relation
      b=dph/a
      eye1=b*b*(dp3+s*(s-dp3)+exps*(dph*s*s-dp3))
      if (j.eq.1) eye0=b*(s-dp1+exps)
      return
c///// small s series expansion
   10 continue
c
      eye1=dp1
      b=dp1
      do 20 n=1,7
      b=b*s
      c=b*c1(n)
      if (c.lt.1.0d-7) go to 30
      eye1=eye1+c
   20 continue
   30 continue
      b=s*s/a
      eye1=eye1*exps*b*b*0.03125d0
      if (j.gt.1) return
c eye0 (neutrons only)=(.5/a)*s**2/2*(sum n=0 to 7:2*s**n/(n!*(n+2))
      eye0=dp1
      b=dp1
      do 40 n=1,7
      b=b*s
      c=b*c0(n)
      if (c.lt.1.0d-7) go to 50
      eye0=eye0+c
   40 continue
   50 continue
      eye0=exps*eye0*s*s*0.25d0/a
      return
      end


************************************************************************
*                                                                      *
C*********This routine is the same as the drein2 in LAHET code**********
      subroutine drein2 (s,a,eye2)
      implicit doubleprecision(a-h,o-z)
      parameter (dp0=0.d0, dp1=1.d0, dp2=2.d0, dp3=3.d0, dp4=4.d0, dph=.
     1 5d0, dp5=5.d0, dp10=1.d1, dpth=dp1/dp3, dppi=3.1415926535898d0,
     2 dp2th=dp2/dp3)
c
c     compute statistical theory emission integrals
c
c     compute third integral
c
c     for s<dph use a series expansion.
c     for s>dph the explicit relationship
c
c     coeficients for series expansions
c
      dimension c2(7)
      data c2 /0.45714286d0,0.125d0,0.02539683d0,0.00416667d0,
     1 0.0005772d0,0.00006944d0,0.0000074d0/
      exps=dp0
      if (s.lt.1.d+02) exps=exp(-s)
      if (s.lt.dph) go to 10
c///// explicit relation
      b=s*s
      eye2=0.25d0*(s*(15.d0-s*(6.d0-s))-15.d0+(15.d0+0.125d0*b*(b-12.d0)
     1 )*exps)/(a*a*a)
      return
   10 continue
c
c///// series expansion
      eye2=dp1
      b=dp1
      do 20 n=1,7
      b=b*s
      c=b*c2(n)
      if (c.lt.1.0d-7) go to 30
      eye2=eye2+c
   20 continue
   30 continue
      b=0.25d0*s*s/a
      eye2=eye2*b*b*b*exps*0.33333333d0
      return
      end


************************************************************************
*                                                                      *
      function efms(z,a)
C/////////////////////////////////////////////////////////////////////
C  EFMS
C  Fission barrier given by Myer & Swaiteski (PRC60,014606,1999)
C=====================================================================
C <variables>
C     a   :   the mass of a fissioning nucleus      (IN)
C     z   :   the charge of a fissioning nucleus    (IN)
C   efms  :   fission barrier  [MeV]                (OUT)
C/////////////////////////////////////////////////////////////////////

      use NGSDATAMOD, only : shellE
      implicit doubleprecision(a-h,o-z)

C... 8/15/1999
      parameter(x0=48.5428d0, x1=34.15d0)
      f1(t)=1.99749d-4*(x0-t)**3
      f2(t)=5.95553d-1-0.124136*(t-x1)

      efms=0.d0
      iz=nint(z)
      ia=nint(a)

      c=1.9+(z-80)/75
      ai=1.-2*(z/a)
      xx=1-c*ai**2
      ss=a**.66667*xx
      x=z**2/a/xx
      if(ia-iz.le.0.or.ia-iz.gt.250.or.iz.gt.150.or.iz.lt.1) then
       sh=0.d0
      else
       nv = 2
       sh=shellE(iz, ia, nv)
      endif
      if(x.ge.x1.and.x.le.x0) then
       efms=ss*f1(x)-sh
      else if(x.ge.20.and.x.lt.x0)then
       efms=ss*f2(x)-sh
      else
       efms=-1.0
      endif

      return
      end


************************************************************************
*                                                                      *
       function gaussn (xmean,sd)
C/////////////////////////////////////////////////////////////////////
C****************This routine is originally from the LAHET code********
C  GAUSSN
c      Gaussian randum number gemerator
c  compute random gaussian number for given
c  mean and s.d.
c  uses mean of sum of 12 uniform r.n"s
c
C/////////////////////////////////////////////////////////////////////
      implicit doubleprecision(a-h,o-z)

      a=0.d0
      do 10 n=1,12
        a=a+rn(0)
 10   continue
      gaussn=(a-6.d0)*sd+xmean
      return
      end


************************************************************************
*                                                                      *
      function paire(iz,in)

      implicit doubleprecision(a-h,o-z)

      logical isz, isn
      parameter (inn=150, iiz=98)
      common /cook/ sz(iiz), sn(inn), con(2), amean(240), pz(iiz),
     1 pn(inn), isz(iiz), isn(inn)
!$OMP THREADPRIVATE(/cook/)

      jz = min(iz,iiz)
      jn = min(in,inn)

      if(in.eq.0) jn=1
      if(iz.eq.0) jz=1
      paire=pz(jz)+pn(jn)
      return

      end


************************************************************************
*                                                                      *
      function rb(a,ia,j)
C/////////////////////////////////////////////////////////////////////
C  RB
C    Calculate Nuclear radius for a geometric cross section
C=====================================================================
C <variables>
C     a   :   mass of nucleus #1                    (IN)
C     z   :   charge  of nucleus #1                 (IN)
C    ia   :   mass of nucleus #2                    (IN)
C    iz   :   charge  of nucleus #2                 (IN)
C     j   :   type of the nufleus #2                (IN)
C    ck   :   transmission probability              (IN)
C   voul  :   Coulomb potential  [MeV]              (OUT)
C/////////////////////////////////////////////////////////////////////

      implicit doubleprecision(a-h,o-z)

      common /options/alev, rcal, ifis
!$OMP THREADPRIVATE(/options/)
      parameter (r0=1.5d0)

      if(rcal.gt.0.0) goto 20

C... Dostrovsky et. al.

      if (j.le.6) then


         if( j .gt. 2 ) then

            rb = a**.333333d0 + dble(ia)**.333333d0

         else

            rb = a**.333333d0

         end if

       rb = rb * r0

C... Matsuse et al. PRC26(1982)2338

      else

       r1=1.12d0*a**.333333d0-0.86d0/a**.333333d0
       r2=1.12d0*dble(ia)**.333333d0-0.86d0/dble(ia)**.333333d0
       rb=r1+r2+2.85d0

      endif

       return

C...Simple form

 20   continue

      if(rcal.eq.10.0)then
       rr=1.5d0
      else
       rr=rcal
      endif
      r1=rr*a**.333333d0
      r2=rr*dble(ia)**.333333d0
      rb=r1+r2
      return
      end


************************************************************************
*                                                                      *
      function vcoul(z,a,iz,ia,ck,j)
C/////////////////////////////////////////////////////////////////////
C  VCOUL
C    Calculate Coulomb potential
C=====================================================================
C <variables>
C     a   :   mass of nucleus #1                    (IN)
C     z   :   charge  of nucleus #1                 (IN)
C    ia   :   mass of nucleus #2                    (IN)
C    iz   :   charge  of nucleus #2                 (IN)
C     j   :   type of the nufleus #2                (IN)
C    ck   :   transmission probability              (IN)
C   voul  :   Coulomb potential  [MeV]              (OUT)
C/////////////////////////////////////////////////////////////////////

      implicit doubleprecision(a-h,o-z)

      parameter (rc=1.70d0, ee=137.0359895d0, hbarc=197.327053d0)
      common /options/alev, rcal, ifis
!$OMP THREADPRIVATE(/options/)


C  No coulomb potential for neutron emission

      if(j.eq.1) then
       vcoul=0.d0
       return
      endif


      if(rcal.gt.0.0) goto 30

      if(j.le.6) then

C Dostrovsky's parameter set

       vcoul=hbarc/rc/ee*ck
       r2 = dble(ia)**.333333d0

       if(j.le.6) r2=1.2d0/rc
       if(j.eq.2) r2=0.d0

       r1 =a**.333333d0
       if(a.le.4.and.z.le.2) r1=1.2d0/rc

       r0=r1+r2

       goto 20

      else

C Matsuse's parameter set ...PRC26(1982)2338

       vcoul=hbarc/ee

       r1=1.12d0*a**0.333333d0-0.86d0/a**0.333333d0
       r2=1.12d0*dble(ia)**.333333d0-0.86d0/dble(ia)**.333333d0

       r0=r1+r2+3.75d0

       goto 20

      endif

C...Simple parameter set
 30   vcoul=hbarc/ee
      if(rcal.eq.10.0) then
       rr=1.5d0
      else
       rr=rcal
      endif
      r1=rr*a**.333333d0
      r2=rr*dble(ia)**.333333d0
      r0=r1+r2

 20   vcoul = vcoul * dble(iz) * z / r0

      return
      end

************************************************************************
*                                                                      *
************************************************************************
*                                                                      *
C--------------------------------------------------------------------
      subroutine eye10(j,aa,izz,nn,u,q,V,delta,smalla,r,s)
C
C     Upgraded version of decay width calculation routine
C
C     variable       IN/OUT
C     j              I       emittor identifier
C     aa             I       residual mass after emission (=A-Aj)
C     izz            I       charge number of res nuclei after emission
C     nn             I       neutron number of res nuclei after emission
C     u              I       excitation energy of nuclei before emission
C     q              I       Q-value
C     V              I       kV in the equation
C     delta          I       pairing energy
C     smalla         I       level density parameter at U-Q-delta-kV
C     r              O       r in the equation
C                                      _______________
C     s              O       s = 2 \/a(U-Q-delta-kV)

C--------------------------------------------------------------------
      use levdenmod, only : getag
      implicit doubleprecision(a-h,o-z)
      include 'param02.inc'

*-----------------------------------------------------------------------
      common /ejectl/ omega(70),ifa(70),ifz(70)

            if( u - q - V .le. 0.0 .or. aa .le. 0.0 ) then

               r = 0.d0
               s = 0.d0
               return

            end if

*----- Basic parameters ------------------------------------------------

               ux  = 2.5d0 + 150.d0 / aa
               ex  = ux + delta

               ax  = getag(ex,izz,nn,isdum)

               tau = sqrt( ax / ux ) - 1.5d0 / ux
               tau = 1.d0 / tau

               e0  = ex - tau*(log(tau)-0.25d0*log(ax)-1.25d0*log(ux)
     &          +2.d0*sqrt(ax*ux))

*----- Numerical integration -------------------------------------------

      eylow = 0.d0
      eyhig = 0.d0
      ebin  = 0.5d0
      elow  = v
      ehig  = elow + ebin

      do while(ehig .le. u - q - ex)
       if(elow .eq. 0.d0) then
        sigma = 0
       elseif(j .le. 2) then
        call sigrc(3-j,elow,int(aa),izz,sigt,sigma,sigs)
       else
        call sighi(dble(ifa(j)),dble(ifz(j)),elow,aa,dble(izz),sigma,
     &   sdum,bdum)
       endif


       ex2aE = min(2.d0 * sqrt(ax*(u - q - elow -delta)),650.d0) ! Avoid overflow

       eylow = eylow + ebin * elow * sigma * pi/12.0
     &  * exp(ex2aE) / ax**0.25d0 / (u - q - elow - delta)**1.25d0

       elow  = elow + ebin ! integration energy lower bound
       ehig  = elow + ebin ! integration energy upper bound
      enddo

      do while(ehig .le. u - q)
       if(elow .eq. 0.d0) then
        sigma = 0
       elseif(j .le. 2) then
        call sigrc(3-j,elow,int(aa),izz,sigt,sigma,sigs)
       else
        call sighi(dble(ifa(j)),dble(ifz(j)),elow,aa,dble(izz),sigma,
     &   sdum,bdum)
       endif


       exE_t = min((u - q - elow) / tau ,650.d0)

       eyhig = eyhig + ebin * elow * sigma * pi/12.0 / tau
     &  * exp(exE_t)

       elow  = elow + ebin ! integration energy lower bound
       ehig  = elow + ebin ! integration energy upper bound
      enddo

      r = eylow + eyhig

*-----------------------------------------------------------------------

      return
      end
************************************************************************

************************************************************************
*
       subroutine nid(j,ia,iz,u,r)
*-----------------------------------------------------------------------
* function :vcoul,q,rhoc,rhod,pde,xsev,dffcs
*-----------------------------------------------------------------------
*   value
*     iaa : mass number of daughter nucleus
*     izz : charge number of d-nucleus
*-----------------------------------------------------------------------
*   numerical integral for decay width of each particle emittion
*                                               by S.S 2017/5/17
*     input
*      ia : mass number of mother nucleus
*      iz : charge number of m-nucleus
*       u : excitation energy
*       j : emitter identifier
*
*     output
*       r(1~6) : decay width [MeV]
*
***********************************************************************
      use NGSDATAMOD, only : energm
      implicit none
      common /ejectl/omega(70),ifa(70),ifz(70)
      double precision u,r,V,q,smalla,gj,s,rhod,ex,h,p0,pf,q1
      double precision ps,pi,q2,omega,qe,ck,cp,dostg
      double precision couk,couc,pde,vcoul,sigma,b,pff
      integer i,ia,iaa,aa,iz,izz,zz,N,ifa,ifz,j,ip
      parameter( pi=3.1415926535898d0 )
      dimension couk(70),couc(70)
      data couk/70*1.d0/
      data couc/70*0.d0/
*---------------------------------------------------------------------
*      open(12,file='mother.txt',position="append")
*      write(12,*)"out",ia,iz
*      close(12)
*---------------------------------------------------------------------
c initialization
        r       = 0.0d0
        q       = 0.0d0
        gj      = 0.0d0
        s       = 0.0d0
        ck      = 0.0d0
        cp      = 0.0d0
*---------------------------------------------------------------------
        ip      = 0 ! select partition option
*---------------------------------------------------------------------
c set
        iaa = ia - ifa(j)
        izz = iz - ifz(j)
*---------------------------------------------------------------------
* set ck,cp
      if (j .ge. 2) then
       couk(2)=dostg(1,dble(iz-ifz(2))) !ck
       couk(3)=couk(2)+0.06
       couk(4)=couk(2)+0.12
       couk(6)=dostg(2,dble(iz-ifz(6)))
       couk(5)=couk(6)-0.06

       couc(2)=dostg(3,dble(iz-ifz(2)))
       couc(3)=couc(2)/2.d0
       couc(4)=couc(2)/3.d0

       ck = couk(j)

       if(j .le. 4)then
        cp = couc(j)
       else
        cp = 0.d0
       endif

      else
          ck = 0.0d0
          cp = 0.0d0
      endif
*---------------------------------------------------------------------
*---------------------------------------------------------------------
*---------------------------------------------------------------------
*---------------------------------------------------------------------
c calculate Q-value
      q1 = energm(iz,ia,1)
      q2 = energm(izz,iaa,1)
      qe = energm(ifz(j),ifa(j),1)
      if (qe .gt. 1.0d10) qe = 0.d0
      q  = q2 - q1 + qe
*---------------------------------------------------------------------
*---------------------------------------------------------------------
c calculate gj
      gj = ( 2.d0 * omega(j) + 1.d0 ) * dble(ifa(j)) !for nid
        if(j.gt.6)then ! for dost formula
           gj = gj * ( 1.d0 + cp ) ! for dost formula
        end if
*---------------------------------------------------------------------
*---------------------------------------------------------------------
* Calculate Coulomb potential
       V = vcoul(dble(izz),dble(iaa),ifz(j),ifa(j),ck,j)

       if(u-q-V .lt. 0.d0) then
           r = 0.d0
           return
       endif

*--------------------------------------------------------------------
*---------------------------------------------------------------------
        if(j.gt.6)then !for nid
            ex = V ! original (for sigma-inv of dost formula
        else
            ex = 0.d0 !for sigma-inverse of niita or other systematics
        endif
*---------------------------------------------------------------------
* select bin width option
* 1 -> calc bin width by partition number
      if( ip .eq. 0 ) then
        N  = 400 ! partition number
       if(j.gt.2)then
            h =(u-q-V)/N ! for dost formula
        else
            h =(u-q)/N ! for systematics
        endif
      else
* 2 -> constant energy bin width to any excitation energy
        h = 0.01d0         ! energy bin width (MeV)
        N = int( u / h ) ! calcurate partition number
        b = mod( u , h )   ! the remainder of u/h
        N = N - 1       ! match up to roop number
      endif
*--------------------------------------------------------------------
*---------------------------------------------------------------------
* do integral

       p0 = pde(j,ia,iaa,iz,izz,ex,q,V,u)

       do i = 1,N-1
           ex = ex + h
           ps = pde(j,ia,iaa,iz,izz,ex,q,V,u)
*--------------------------------------------------------------------
*--------------------------------------------------------------------
           s  = s + ps
       enddo

       pf  = pde(j,ia,iaa,iz,izz,u-q,q,V,u)

      if( ip .eq. 0 ) then
       r   = gj * h * ( p0 + s * 2.d0 + pf ) / 2.d0
      else
       pff = pde(j,ia,iaa,iz,izz,u-q-b,q,V,u)
       r   = gj * ( h * ( p0 + s * 2.d0 - pff )
     &                     + ( pf + pff ) * b ) / 2.d0
      endif


      end
*--------------------------------------------------------------------

*********************************************************************
       function pde(j,ia,iaa,iz,izz,ex,q,V,u)
*     incp = 1 proton, 2 neutron for Niita sistematics
*     pde     : integrand
*     rhod    : revel dencity of daughter nucleous
*     csinvex : invers cross section * emittion energy
*     j       : emitt er identifier
*     ia      : mass number of mother
*     iaa     : /mass number of daughter
*     iz      : atomic number of mother
*     izz     : atomic number of daughter
*     ex      : energy bin
*     q       : Q-value
*     u       : excitation energy
*********************************************************************
      use levdenmod, only:rho_levden, getag

      implicit none
      double precision ex,q,cs,rho,rhod,p,del,paire,sigt,u,V,sp
      double precision pde,pi,rhoc,rhof,ux,exx,bmax
      double precision sigr,sigs,csinvex,smalla,xsev,a,z
      integer j,ia,iaa,iz,izz,incp,nn,isdum,nv
*--------------------------------------------------------------------
c initialization
      pde      = 0.0d0
      rhod     = 0.0d0
      csinvex  = 0.0d0
      sigr     = 0.0d0
      pi       = 3.1415926535897932d0
*--------------------------------------------------------------------
c calculate integrand
      nn     = dble(iaa-izz)

      if (j .eq. 1)then
          incp =  2 ! need reconfirmation
      else
          incp = 1
      endif

*--------------------------------------------------------------------
*--------------------------------------------------------------------

      sp = -1.d0 ! spin unknown
      nv = 2
      rhod   = rho_levden(iaa,izz,u-q-ex,sp,nv)

*--------------------------------------------------------------------
* XS-inv * ex
      csinvex =  xsev(j,ia,iaa,iz,izz,ex,V)

*-------------------------------------------------------------------
      pde = csinvex * rhod
*-------------------------------------------------------------------
      if (pde .lt. 0.d0) pde =0.d0
      return
      end
*--------------------------------------------------------------------

*--------------------------------------------------------------------
      function xsev(j,ia,iaa,iz,izz,ex,V)
*
*     ia      : mass number of mother
*     iaa     : mass number of daughter
*
*     inverse Cross section
*
*     ixnp      : 1 -> niita systematics
*     (n,p)       2 -> KUROTAMA
*
*     ixd       : 1 -> nasa
*     (d)         2 -> shen
*                 3 -> KUROTAMA
*                 4 -> MWO
*
*     ixhad     : 1 -> nasa
*     (t~4He)     2 -> shen
*                 3 -> KUROTAMA
*
*     !! 0 -> DOST formura (original gem)!!
*
*--------------------------------------------------------------------
      implicit none
      common /ejectl/omega(70),ifa(70),ifz(70)
      common /ckurotama/ dsck
      double precision ex,alp,beta,bett,xsev,rmass,rb,V,sigg
      double precision aa,zz,sigr,sigt,omega,dffcs,pi,sigs,bmax
      double precision ap,zp,at,zt,dsck,ea
      integer iaa,izz,j,ifa,ifz,ia,iz,incp,ixp,ixn,ixd,ixhad
*--------------------------------------------------------------------
      parameter (ixn = 1, ixp = 2, ixd = 4, ixhad = 1)
      parameter (pi = 3.1415926535897932d0)

*--------------------------------------------------------------------
      if (j .eq. 1)then
          incp =  2
      else
          incp = 1
      endif
*--------------------------------------------------------------------
* initialization for XS-inv
          ap = dble(ifa(j))
          zp = dble(ifz(j))
          at = dble(ia)
          zt = dble(iz)
          ea = ex / ap
*------------------------------------------------------------------
      if (j .eq. 1) then

         if (ixn .le. 0) then
          xsev = dffcs(j,ia,iaa,iz,izz,ex,V)
         else

          select case(ixn)

          case(1)
           if (ex .gt. 0.d0) then
            call sigrc(incp,ex,ia,iz,sigt,sigr,sigs) !niita
           else
            sigr = 0.d0
           endif
           sigr = sigr

          case(2)
           if (ex .gt. 0.d0) then
            call kurotama0(ap,zp,ea,at,zt,dsck,sigr,bmax) ! KUROTAMA
           sigr = sigr /100 ! fm^2 -> b
           else
            sigr = 0.d0
           endif

          end select
           sigr = sigr * 1.d2 / pi
           xsev = sigr * ex

         endif
*------------------------------------------------------------------
c p
      else if (j .eq. 2) then
         if (ixp .le. 0) then
          xsev = dffcs(j,ia,iaa,iz,izz,ex,V)
         else
          select case(ixp)

          case(1)
           if (ex .gt. 0.d0) then
            call sigrc(incp,ex,ia,iz,sigt,sigr,sigs) !niita
           else
            sigr = 0.d0
           endif
           sigr = sigr

          case(2)
           if (ex .gt. 0.d0) then
            call kurotama0(ap,zp,ea,at,zt,dsck,sigr,bmax) ! KUROTAMA
           sigr = sigr /100 ! fm^2 -> b
           else
            sigr = 0.d0
           endif

          end select
           sigr = sigr * 1.d2 / pi
           xsev = sigr * ex

         endif
*-------------------------------------------------------------------
c deuteron
      else if (j .eq.3)then

         if (ixd .le. 0) then
          xsev = dffcs(j,ia,iaa,iz,izz,ex,V)
         else

          select case(ixd)

          case(1)
           if (ex .gt. 0.d0) then
            call nasa(ap,zp,ex,at,zt,sigr,sigs,bmax) ! NASA
           else
            sigr = 0.d0
           endif
           sigr = sigr

          case(2)
           if (ex .gt. 0.d0) then
            call shen(ap,zp,ex,at,zt,sigr,sigs,bmax) ! SHEN
           else
            sigr = 0.d0
           endif
           sigr = sigr

          case(3)
           if (ex .gt. 0.d0) then
            call kurotama0(ap,zp,ea,at,zt,dsck,sigr,bmax) ! KUROTAMA
           else
            sigr = 0.d0
           endif
           sigr = sigr /100 ! fm^2 -> b
           sigr = sigr

          case(4)
           if (ex .gt. 0.d0) then
            call deumino(ex,at,zt,sigr,sigs,bmax) !MWO
           else
            sigr = 0.d0
           endif
           sigr = sigr

          end select

          sigr  = sigr * 1.d2 / pi
          xsev  = sigr * ex

         endif
*-------------------------------------------------------------------
c t~4He
      elseif (j.ge.4)then
         if (ixhad .le. 0) then
          xsev = dffcs(j,ia,iaa,iz,izz,ex,V)
         else

          select case(ixhad)

          case(1)
           if (ex .gt. 0.d0) then
            call nasa(ap,zp,ex,at,zt,sigr,sigs,bmax) ! NASA
           else
            sigr = 0.d0
           endif

          case(2)
           if (ex .gt. 0.d0) then
            call shen(ap,zp,ex,at,zt,sigr,sigs,bmax) ! SHEN
           else
            sigr = 0.d0
           endif


          case(3)
           if (ex .gt. 0.d0) then
            call kurotama0(ap,zp,ea,at,zt,dsck,sigr,bmax) ! KUROTAMA
           else
            sigr = 0.d0
           endif
           sigr = sigr /100 ! fm^2 -> b


          end select

          sigr  = sigr * 1.d2 / pi
          xsev  = sigr * ex

         endif
        endif
*-------------------------------------------------------------------

      if (xsev .lt. 0.d0) xsev = 0.d0
      return
      end

*--------------------------------------------------------------------
      function dffcs(j,ia,iaa,iz,izz,ex,V)
*     XS-inv * ex (Dostrovfsky et al)
*
*     ia      : mass number of mother
*     iaa     : mass number of daughter
*     ex      : kinetic energy of emitter
*
*--------------------------------------------------------------------
      implicit none
      common /ejectl/omega(70),ifa(70),ifz(70)
      double precision ex,alp,beta,bett,dffcs,rmass,rb,V,sigg
      double precision omega,aa,zz
      integer iaa,izz,j,ifa,ifz,ia,iz
*--------------------------------------------------------------------
* for neutron emittion

      if (j .eq. 1) then
          alp  = 0.76d0 + 1.93d0 / iaa**.333333333d0
          beta = ( 1.66d0 / iaa**.6666666667d0 -5.d-2 ) / alp
          beta = max(beta,0.d0)
          bett = beta
*--------------------------------------------------------------------
* for proton emittion

      else
          alp  = 1.d0
          bett = -V
      endif
*--------------------------------------------------------------------
* geometoric cross section

      rmass = rb(dble(iaa),ifa(j),j)
      sigg  = rmass * rmass
*--------------------------------------------------------------------

          dffcs = alp * sigg * ( ex + bett  ) ! XS * ex

      return
      end
*--------------------------------------------------------------------


