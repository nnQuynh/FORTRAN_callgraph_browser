************************************************************************
*                                                                      *
      subroutine jamin(iqmd,iprj,kprj,eein,mmas,mchg,bmax0)
*                                                                      *
*                                                                      *
*       control routine of JAM and JAMQMD                              *
*       call JAMJAMIN or JAMQMDIN                                      *
*       modified by K.Niita on 2004/11/16                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       iqmd    : 0; jam mode, 1; jamqmd mode                          *
*       iprj    : particle type of projectile                          *
*       kprj    : kf code of projectile                                *
*       eein    : energy of projectile (MeV)                           *
*       mmas    : mass of target                                       *
*       mchg    : charge of target                                     *
*       bmax0   : max impact parameter                                 *
*                                                                      *
************************************************************************
      use QMD_COOD2_MOD, only : coodalloc, lfirst

*-----------------------------------------------------------------------

      include 'jam1.inc'
      include 'jam2.inc'

*-----------------------------------------------------------------------

      common /pionflag/ iqdflag, ipionflag, istrflag, iphreturn
!$OMP THREADPRIVATE(/pionflag/)

      data initjm /0/
      save initjm !FURUTA
!$OMP THREADPRIVATE(initjm)
*-----------------------------------------------------------------------
*     initialization
*-----------------------------------------------------------------------

         if( initjm .eq. 0 ) then

            initjm = initjm + 1

            temp = jamcomp(2212)

         end if

*-----------------------------------------------------------------------

         if( iqmd .eq. 0 ) then

            call jamjin(iprj,kprj,eein,mmas,mchg,bmax0)

cABE add @2014/07/29, to retry cal. for photonuclear reaction
            if (iphreturn .eq. 1) return

         else

               IF ( .not. lfirst ) then
                lfirst = .true.
                call coodalloc
               ENDIF

            call jamqin(iprj,kprj,eein,mmas,mchg,bmax0)

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine jamqin(iprj,kprj,eein,mmas,mchg,bmax0)
*                                                                      *
*                                                                      *
*       control routine of JAMQMD                                      *
*       modified by K.Niita on 2011/07/15                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       iprj    : particle type of projectile                          *
*       kprj    : kf code of projectile                                *
*       eein    : energy of projectile (MeV)                           *
*       mmas    : mass of target                                       *
*       mchg    : charge of target                                     *
*       bmax0   : max impact parameter                                 *
*                                                                      *
*     output:                                                          *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      include 'jam1.inc'
      include 'jam2.inc'

      include 'param01.inc'

      parameter ( pi  = 3.1415926535898d0 )

*-----------------------------------------------------------------------

      common /jamnmtc/ kfjam
!$OMP THREADPRIVATE(/jamnmtc/)
      common /jamqmd/ dtqjam
!$OMP THREADPRIVATE(/jamqmd/)

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)

      common /clustp/ rumpat(0:20), numpat(0:20)
!$OMP THREADPRIVATE(/clustp/)

      common /jamincl/ efermi, ifermi
!$OMP THREADPRIVATE(/jamincl/)
      common /cputim/ stime(40), cputm(40)
      common /evaval/ sepc(2,0:100)
      common /coln01/ iccoll
!$OMP THREADPRIVATE(/coln01/)

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

*-----------------------------------------------------------------------

      common /swich2/ icfg, imany, icpus, idatm
!$OMP THREADPRIVATE(/swich2/)
      common /input1/ mstq1(mxpa1), parq1(mxpa1)
!$OMP THREADPRIVATE(/input1/)

      common /swich3/ ielst, jelst, kelst
!$OMP THREADPRIVATE(/swich3/)

      common /framtr/ betafr(0:2), gammfr(0:2)
!$OMP THREADPRIVATE(/framtr/)
      common /crshi/  bplus, icrhi, ijudg, imadj, iqmax
      common /const8/ betta, gamta, tamas, radta
!$OMP THREADPRIVATE(/const8/)

*-----------------------------------------------------------------------

      character chau*8

      character frame*8,proj*8,targ*8,cwin*15
      character elab*12

      character pname*16

      data rmnuc / 0.93895 /

*-----------------------------------------------------------------------
*     set initial values by Furuta
*-----------------------------------------------------------------------

            mstc(56) = 1         ! Pauli
            parc(5)  = -1.d0     ! separation

*-----------------------------------------------------------------------
*     initial values from PHITS
*-----------------------------------------------------------------------

         if( kprj .gt. 1000000 ) then

            maspr = ibryf(iprj,kprj)
            mchpr = ichgf(iprj,kprj)
            kfjam = 0

         else

            kfjam = kprj

         end if

*-----------------------------------------------------------------------
*        target and projectile
*-----------------------------------------------------------------------

            call chname(ic,mmas,mchg,chau)

            if( ic .eq. 1 ) then
               write(*,*) 'chau = ', chau
               write(*,*) ' Skip : Target Nucleus is strange in JAMQMD'

               nclst = -2
               return

            end if

            targ  = chau

         if( kfjam .eq. 0 ) then

            call chname(ic,maspr,mchpr,chau)

            if( ic .eq. 1 ) then
               write(*,*) 'chau = ', chau
               write(*,*) ' Skip : Proj. Nucleus is strange in JAMQMD'

               nclst = -2
               return

            end if

            proj  = chau

         else

            proj  = 'other'

         end if

*-----------------------------------------------------------------------
*     bmax for QMD calculation
*-----------------------------------------------------------------------


         bmax = 1.2d0 * ( dble(mmas)**(1./3.) + dble(maspr)**(1./3.) )
     &            + bplus ! 2016/11/04   Optimized and verified for JAMQMD.

*-----------------------------------------------------------------------
*     initialization for JAM
*     low energy cutoff for collision : default is 0.02 GeV <== changed
*                                                  2004/04/12
*-----------------------------------------------------------------------

         if( mmas .gt. 5 ) then

            parc(38) = 0.02

            parc(32) = 55.0
            parc(33) = 55.0

         else if( mmas .gt. 1 ) then

            parc(38) = 0.000

            parc(32) = 55.0
            parc(33) = 55.0

         else if( mmas .eq. 1 ) then

            parc(38) = 0.000

            parc(32) = 400.0
            parc(33) = 150.0

         end if

*-----------------------------------------------------------------------

            fname(1) = '0'

         if( imadj .eq. 0 ) then

            frame    = 'cm'

         else

            frame    = 'nn'

         end if

*-----------------------------------------------------------------------

            mstc(36) = 56
            mstc(37) = 0
            mstc(38) = 6 !58
            mstc(39) = 0
            mstc(40) = 0

            mstc(12) = 0
            mstc(13) = 0
            mstc(14) = 0

            mstc(43) = 1

            mstc(17) = 0
            mstc(71) = 1

*-----------------------------------------------------------------------
*        cascade mode and specify QMD by mstc(190)
*-----------------------------------------------------------------------

            mstc(6)   = 0
            mstc(190) = 1

*-----------------------------------------------------------------------
*        incident energy and initialization of JAM
*-----------------------------------------------------------------------

         if( kfjam .eq. 0 ) then

            ein = eein / 1000.0 / dble( maspr )

         else

            ein = eein / 1000.0

         end if

            write( elab, '(f12.5)') ein

            cwin = elab//'gev'

*-----------------------------------------------------------------------
*     initialization for JAM
*-----------------------------------------------------------------------

            if( ein .le. 3.5 ) then

               dt     = 150
               nstep  = 1
               dtqjam = 1.0d0

            else

               dt     = 100
               nstep  = 1
               dtqjam = 0.5d0

            end if

               iev    = 1
               mevent = 1

               bmin   = 0.0
               bmax   = - bmax

            if( mmas .le. 2 ) bmax = -2.0

*-----------------------------------------------------------------------
*     Initialize QMD
*-----------------------------------------------------------------------

               ieo  = 6
               ierr = 0

               mstq1(1) = 0
               mstq1(2) = ibryf(iprj,kprj)
               mstq1(3) = ichgf(iprj,kprj)

               msp = mstq1(3)
               msn = mstq1(2) - mstq1(3)

               mstq1(4) = 0
               mstq1(5) = mmas
               mstq1(6) = mchg

*-----------------------------------------------------------------------

               mstq1(190) = 0         ! 0:no, 1:moving frame adjust
               mstq1(191) = ijudg    ! 0:no, 1:new ielst judge

               mstq1(17)  = 1        ! ielst for judge

*-----------------------------------------------------------------------

               parq1(1)   = ein
               parq1(2)   = -1.0
               parq1(3)   = 0.0

               parq1(4)   = -bmax

               mstq1(11)  = 0

               mstq1(7)   = 1
               mstq1(8)   = 150
               parq1(5)   = 1.0

               mstq1(90)  = 1          ! 0:no, 1:ground energy adjust

               mstq1(9)   = imadj + 1  ! 0: lab,  1: cm,  2: nn
               mstq1(32)  = 1          ! 0:no, 1:rel.correction

*-----------------------------------------------------------------------
*           nucleon-nucleus collisions
*-----------------------------------------------------------------------

            if( mstq1(2) .eq. 1 .or. mstq1(5) .eq. 1 ) then

               mstq1(9)   = 0
               mstq1(17)  = 3        ! ielst for judge
               mstq1(190) = 0        ! 0:no, 1:moving frame adjust

            end if

*-----------------------------------------------------------------------

               parq1(120) = 1.0      ! sdmemin

               icfg = 4

*-----------------------------------------------------------------------
*        JAM part
*-----------------------------------------------------------------------

            call jaminit(mevent,bmin,bmax,dt,nstep,
     &                   frame,proj,targ,cwin)

            call jamevt(iev,ierr)

               if( ierr .ne. 0 ) then

                  nclst = -2
                  return

               end if

             iccoll =  mstd(41) + mstd(42)

*-----------------------------------------------------------------------
*        Final QMD part
*-----------------------------------------------------------------------

                  call cldistR

*-----------------------------------------------------------------------
*           elastic judge and frame trnasform to lab system
*-----------------------------------------------------------------------



                     ifrm = 0

               do i = 1, nclst

                     px = qclust(1,i)
                     py = qclust(2,i)
                     pz = qclust(3,i)
                     et = qclust(4,i)
                     rm = qclust(5,i)

                     gamm = gammfr(ifrm)
                     beta = betafr(ifrm)

                     pz = pz * gamm - beta * gamm * et
                     et = sqrt( px**2 + py**2 + pz**2 + rm**2 )

                     qclust(3,i)  = pz
                     qclust(4,i)  = et
                     qclust(7,i)  = ( et - rm ) * 1000.

               end do

*-----------------------------------------------------------------------
*           x-y rotation by randum
*-----------------------------------------------------------------------

                  theta = 2.0d0 * pi * rn(0)

               do i = 1, nclst

                  px = qclust(1,i)
                  py = qclust(2,i)

                  qclust(1,i)  = px * cos( theta ) - py * sin( theta )
                  qclust(2,i)  = px * sin( theta ) + py * cos( theta )

               end do

*-----------------------------------------------------------------------
*     error or nothing happened
*-----------------------------------------------------------------------

      return

  999 continue

                  nclst = -1

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine jamjin(iprj,kprj,eein,mmas,mchg,bmax0)
*                                                                      *
*                                                                      *
*       control routine of JAM                                         *
*       modified by K.Niita on 2011/07/15                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       iprj    : particle type of projectile                          *
*       kprj    : kf code of projectile                                *
*       eein    : energy of projectile (MeV)                           *
*       mmas    : mass of target                                       *
*       mchg    : charge of target                                     *
*       bmax0   : max impact parameter                                 *
*                                                                      *
*     output:                                                          *
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
*        numpat(i) : total number of out going particles or nuclei     *
*                                                                      *
*                i =  0, nuclei                                        *
*                  =  1, proton                                        *
*                  =  2, neutron                                       *
*                  =  3, pi+                                           *
*                  =  4, pi0                                           *
*                  =  5, pi-                                           *
*                  =  6, mu+                                           *
*                  =  7, mu-                                           *
*                  =  8, K+                                            *
*                  =  9, K0                                            *
*                  = 10, K-                                            *
*                                                                      *
*                  = 11, other particles                               *
*                                                                      *
*                  = 12, electron                                      *
*                  = 13, positron                                      *
*                  = 14, photon                                        *
*                                                                      *
*                  = 15, deuteron                                      *
*                  = 16, triton                                        *
*                  = 17, 3He                                           *
*                  = 18, Alpha                                         *
*                  = 19, residual nucleus                              *
*                                                                      *
*                                                                      *
************************************************************************
      use NGSDATAMOD, only : bindeg

*-----------------------------------------------------------------------

      include 'jam1.inc'
      include 'jam2.inc'

      include 'err.inc'

      parameter ( pi  = 3.1415926535898d0 )

*-----------------------------------------------------------------------

      common /jamnmtc/ kfjam
!$OMP THREADPRIVATE(/jamnmtc/)

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /clustp/ rumpat(0:20), numpat(0:20)
!$OMP THREADPRIVATE(/clustp/)

      common /jamincl/ efermi, ifermi
!$OMP THREADPRIVATE(/jamincl/)
      common /cputim/ stime(40), cputm(40)
      common /evaval/ sepc(2,0:100)

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /crshi/  bplus, icrhi, ijudg, imadj, iqmax

      common /pionflag/ iqdflag, ipionflag, istrflag, iphreturn
!$OMP THREADPRIVATE(/pionflag/)

*-----------------------------------------------------------------------

      character chau*8

      character frame*8,proj*8,targ*8,cwin*15
      character elab*12

      character pname*16

      data rmnuc / 0.93895 /

*-----------------------------------------------------------------------
*     2011/11/28 by Kasahara and Furuta
*-----------------------------------------------------------------------

            ierr = 0

*-----------------------------------------------------------------------
*     set initial values by Furuta
*-----------------------------------------------------------------------

            mstc(56) = 1         ! Pauli
            parc(5)  = -1.d0     ! separation

*-----------------------------------------------------------------------
*     initial values from PHITS
*-----------------------------------------------------------------------

            ibrym = ibryf(iprj,kprj)
            projm = rmtyp(iprj,kprj) / 1000.0

            kfjam =  kprj
            proj  = 'other'

*-----------------------------------------------------------------------

            if( iprj .lt. 3 ) then

                  pima = 0.d0

            else

               if( ibrym .eq. 0 ) then

                  pima = projm

               else if( ibrym .lt. 0 ) then

                  pima = projm + rmnuc

               else if( ibrym .gt. 0 ) then

                  pima = projm - rmnuc

               end if

            end if

*-----------------------------------------------------------------------
*     Coulomb barrier and Potential depth
*-----------------------------------------------------------------------

         if( mmas .gt. 2 ) then

            ebinn  = 0.0
            eclbp  = 1.44e-3 * dble(mchg)
     &             / ( sepc(1,2) + sepc(2,2) * dble(mmas)**(1./3.) )

            ecutin = 0.050
            efermi = 0.037
            ifermi = 1

         else

            ebinn  = 0.0
            eclbp  = 0.0
            ecutin = 0.0
            efermi = 0.037
            ifermi = 0

         end if

*-----------------------------------------------------------------------
*     initialization for JAM
*-----------------------------------------------------------------------

            dt     = 100.0
            nstep  = 1

            iev    = 1
            mevent = 1

            bmin   = 0.0
            bmax   = - bmax0

         if( mmas .le. 2 ) bmax = -2.0

*-----------------------------------------------------------------------
*     initialization for JAM
*     low energy cutoff for collision : default is 0.02 GeV <== changed
*                                                  2004/04/12
*-----------------------------------------------------------------------

         if( mmas .gt. 5 ) then

            parc(38) = 0.02

            parc(32) = 55.0
            parc(33) = 55.0

         else if( mmas .gt. 1 ) then

            parc(38) = 0.000

            parc(32) = 55.0
            parc(33) = 55.0

         else if( mmas .eq. 1 ) then

            parc(38) = 0.000

            parc(32) = 400.0
            parc(33) = 150.0

         end if

*-----------------------------------------------------------------------

            fname(1) = '0'
            frame    = 'lab'

            call chname(ic,mmas,mchg,chau)

            if( ic .eq. 1 ) then
               write(*,*) 'chau = ', chau
               write(*,*) 'mmas = ', mmas, 'mchg = ', mchg
         write(ErrCha,*) ' Error : Initial Nucleus is strange in JAM'
         ErrID = 'L:745/R:jamjin/F:jamin.f' !E00_006_001
         call ErrWrite(ErrID,ErrCha)

               nclst = -1
               return

            end if

            targ     = chau

            mstc(36) = 56
            mstc(37) = 0
            mstc(38) = 6 !58
            mstc(39) = 0
            mstc(40) = 0

            mstc(12) = 0
            mstc(13) = 0
            mstc(14) = 0

            mstc(43) = 1

            mstc(17) = 0
            mstc(71) = 1

*-----------------------------------------------------------------------
*        cascade mode and specify QMD by mstc(190)
*-----------------------------------------------------------------------

            mstc(6)   = 0
            mstc(190) = 0

*-----------------------------------------------------------------------
*        incident energy and initialization of JAM
*-----------------------------------------------------------------------

            ein = eein / 1000.0 + ecutin

            write( elab, '(f12.5)') ein

            cwin = elab//'gev'

*-----------------------------------------------------------------------
*        cutoff for low energy meson with nucleon blow 1 MeV
*-----------------------------------------------------------------------

         if( mmas .eq. 1 .and. ibrym .eq. 0 .and. ein .lt. 0.001 ) then

            nclst = -1
            return

         end if

*-----------------------------------------------------------------------
*        JAM
*-----------------------------------------------------------------------

            call jaminit(mevent,bmin,bmax,dt,nstep,
     &                   frame,proj,targ,cwin)

            call jamevt(iev,ierr)

cABE add @2014/07/29, to retry cal. for photonuclear reaction
            if (iphreturn .eq. 1) return

*-----------------------------------------------------------------------
*        booking of the result of JAM
*-----------------------------------------------------------------------

            eelse = 0.0

            poutx = 0.0
            pouty = 0.0
            poutz = 0.0

            nclst = 0

            nneut = 0
            nprot = 0

            mnucl = 0
            mneut = 0
            mprot = 0

            enucl = 0

            npipo = 0
            nping = 0
            npine = 0

            nkapo = 0
            nkang = 0
            nkane = 0
            nmupo = 0
            nmune = 0
            ngamm = 0
            nothe = 0

            epion = 0
            ekaon = 0
            egamm = 0
            emuon = 0

*-----------------------------------------------------------------------
*        without collision
*-----------------------------------------------------------------------

         if( mstd(41) + mstd(42) .eq. 0 ) then

cABE change @2014/07/24, to avoid return for gamma cal.
            if (iprj .ne. 14) then
               nclst = -1
               return
            endif

         end if

*-----------------------------------------------------------------------
*     do loop for all particles in JAM
*-----------------------------------------------------------------------

         do 100 i = 1, nv

            kf = k(2,i)

            if( kf .eq. 0 ) goto 100

            if( kf .eq. 2112 .or.
     &          kf .eq. 2212 ) then

               if( abs(k(7,i)) .ne. 1 .and. (
     &             ( kf .eq. 2112 .and.
     &               p(4,i) - p(5,i) .ge. ecutin + ebinn ) .or.
     &             ( kf .eq. 2212 .and.
     &               p(4,i) - p(5,i) .ge. ecutin + eclbp ) ) ) then

                     ebinp = ecutin

                     enucl = enucl + p(4,i) - p(5,i) - ebinp

                  if( kf .eq. 2212 ) then

                     nprot = nprot + 1

                     ipid  = 1
                     ippad = 1
                     ipprt = 1
                     ipneu = 0
                     ipchg = 1

                  else if( kf .eq. 2112 ) then

                     nneut = nneut + 1

                     ipid  = 2
                     ippad = 2
                     ipprt = 0
                     ipneu = 1
                     ipchg = 0

                  end if

               else

                     mnucl = mnucl + 1

                  if( kf .eq. 2212 ) then

                     mprot = mprot + 1

                  else if( kf .eq. 2112 ) then

                     mneut = mneut + 1

                  end if

                     goto 100

               end if

            else if( kf .eq. -211 .or.
     &               kf .eq.  111 .or.
     &               kf .eq.  211 ) then

                     ebinp = 0.0

                     epion = epion + p(4,i)

                     ipid  = 3
                     ipprt = 0
                     ipneu = 0

                  if( kf .eq. 211 ) then

                     npipo = npipo + 1

                     ippad = 3
                     ipchg = 1

                  else if( kf .eq.  111 ) then

                     npine = npine + 1

                     ippad = 4
                     ipchg = 0

                  else if( kf .eq. -211 ) then

                     nping = nping + 1

                     ippad = 5
                     ipchg = -1

                  end if

            else if( kf .eq. 22 ) then

                     ebinp = 0.0

                     egamm = egamm + p(4,i)

                     ipid  = 4
                     ipprt = 0
                     ipneu = 0
                     ippad = 14
                     ipchg = 0

                     ngamm = ngamm + 1

            else if( kf .eq. -321 .or.
     &               kf .eq.  311 .or.
     &               kf .eq.  321 ) then

                     ebinp = 0.0

                     ekaon = ekaon + p(4,i)

                     ipid  = 5
                     ipprt = 0
                     ipneu = 0

                  if( kf .eq. 321 ) then

                     nkapo = nkapo + 1

                     ippad = 8
                     ipchg = 1

                  else if( kf .eq. 311 ) then

                     nkane = nkane + 1

                     ippad = 9
                     ipchg = 0

                  else if( kf .eq. -321 ) then

                     nkang = nkang + 1

                     ippad = 10
                     ipchg = -1

                  end if

            else if( kf .eq. -13 .or.
     &               kf .eq.  13 ) then

                     ebinp = 0.0

                     emuon = emuon + p(4,i)

                     ipid  = 6
                     ipprt = 0
                     ipneu = 0

                  if( kf .eq. 13 ) then

                     nmupo = nmupo + 1

                     ippad = 6
                     ipchg = 1

                  else if( kf .eq. -13 ) then

                     nmune = nmune + 1

                     ippad = 7
                     ipchg = 0

                  end if

*-----------------------------------------------------------------------
*           the other particles
*-----------------------------------------------------------------------

            else

                     nothe = nothe + 1

                     ipid  = 7
                     ipprt = 0
                     ipneu = 0
                     ippad = 11
                     ipchg = ichgf(11,kf)

                     ebinp = 0.0

                  if( k(9,i) .eq. 3 ) then

                     eelss = p(4,i) - rmnuc

                  else if( k(9,i) .eq. -3 ) then

                     eelss = p(4,i) + rmnuc

                  else

                     eelss = p(4,i)

                  end if

                     eelse = eelse + eelss

*-----------------------------------------------------------------------
*                 debug of other kind of particles
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------

            end if

*-----------------------------------------------------------------------

                  pouta = sqrt( max( 0.0d0,
     &                        ( p(4,i) - ebinp )**2 - p(5,i)**2 ) )
                  psqrt = sqrt( p(1,i)**2 + p(2,i)**2 + p(3,i)**2 )

               if( psqrt .gt. 0.0d0 ) then

                  pxrv  = pouta * p(1,i) / psqrt
                  pyrv  = pouta * p(2,i) / psqrt
                  pzrv  = pouta * p(3,i) / psqrt

               else

                  pxrv  = 0.0
                  pyrv  = 0.0
                  pzrv  = 0.0

               end if

                  poutx = poutx + pxrv
                  pouty = pouty + pyrv
                  poutz = poutz + pzrv

*-----------------------------------------------------------------------
*        booking of outgoing particles
*-----------------------------------------------------------------------

                  nclst = nclst + 1

                  iclust(nclst)    = ipid

                  jclust(0,nclst)  = 0
                  jclust(1,nclst)  = ipprt
                  jclust(2,nclst)  = ipneu
                  jclust(3,nclst)  = ippad
                  jclust(4,nclst)  = 0
                  jclust(5,nclst)  = ipchg
                  jclust(6,nclst)  = k(9,i) / 3
                  jclust(7,nclst)  = kf
                  jclust(8,nclst)  = 0

                  qclust(0,nclst)  = pard(2)
                  qclust(1,nclst)  = pxrv
                  qclust(2,nclst)  = pyrv
                  qclust(3,nclst)  = pzrv
                  qclust(4,nclst)  = p(4,i) - ebinp
                  qclust(5,nclst)  = p(5,i)
                  qclust(6,nclst)  = 0.0
                  qclust(7,nclst)  = ( p(4,i) - ebinp - p(5,i) ) * 1000.
                  qclust(8,nclst)  = 1.0
                  qclust(9,nclst)  = 0.0
                  qclust(10,nclst) = 0.0d0
                  qclust(11,nclst) = 0.0d0
                  qclust(12,nclst) = 0.0d0

  100    continue

*-----------------------------------------------------------------------
*        booking of number, energy and momentum of out going particles
*-----------------------------------------------------------------------

                  numpat(1)  = nprot
                  numpat(2)  = nneut
                  numpat(3)  = npipo
                  numpat(4)  = npine
                  numpat(5)  = nping
                  numpat(6)  = nmupo
                  numpat(7)  = nmune
                  numpat(8)  = nkapo
                  numpat(9)  = nkane
                  numpat(10) = nkang
                  numpat(11) = nothe
                  numpat(14) = ngamm

                  rumpat(1)  = nprot
                  rumpat(2)  = nneut
                  rumpat(3)  = npipo
                  rumpat(4)  = npine
                  rumpat(5)  = nping
                  rumpat(6)  = nmupo
                  rumpat(7)  = nmune
                  rumpat(8)  = nkapo
                  rumpat(9)  = nkane
                  rumpat(10) = nkang
                  rumpat(11) = nothe
                  rumpat(14) = ngamm

*-----------------------------------------------------------------------
*        residual nucleus
*-----------------------------------------------------------------------

            if( mprot + mneut .gt. 0 .and.
     &          mprot .ge. 0 .and.   ! gt -> ge (S.H.20140812)
     &          mneut .ge. 0 ) then  ! gt -> ge (S.H.20140812)

C S.Hashimoto added IF statement for clusters composed purely of
C neutrons or protons. (2014.10.28)
             if ( ( mneut .gt. 1 .and. mprot .eq. 0 ) .or.
     &              ( mprot .gt. 1 .and. mneut .eq. 0 ) ) then

                nresidue = max( mneut, mprot )

                gein = eein / 1000.0
                piabs = dsqrt( gein**2 + 2.0 * gein * projm )
                presx = - poutx
                presy = - pouty
                presz = - poutz + piabs
                pabst = dsqrt( presx**2 + presy**2 + presz**2 )
                rsmas = 0.93827d0 * dble( mprot )
     &               + 0.93958d0 * dble( mneut )
                etota = dsqrt( pabst**2 + rsmas**2 )
                erres = etota - rsmas
                exres = gein + pima
     &                - enucl - epion - ekaon
     &                - emuon - egamm - eelse
     &                - erres
     &                - bindeg(mchg,mmas-mchg) / 1000.0

              if ( exres + gein .lt. 0d0 ) then
                 nclst = -1
                 return
              end if

                 absexres = dabs(exres)
              if ( absexres .gt. 0d0 ) then
                 niteration = 10000
                 epsiloni = 1d-5 !(10 keV)

               if ( exres .lt. 0d0 ) then
                  facti = 1d0 - 1d-5
               else
                  facti = 1d0 + 1d-5
               end if

               iteration = 0
               exresini = exres
               do while ( exresini*exres .gt. 0d0 .and.
     &              iteration .le. niteration )
                  iteration = iteration + 1
                  poutx = 0d0
                  pouty = 0d0
                  poutz = 0d0
                  enucl = 0d0
                  epion = 0d0
                  ekaon = 0d0
                  emuon = 0d0
                  egamm = 0d0
                  eelse = 0d0

                do iclst=1,nclst
                   qclust(1,iclst) = qclust(1,iclst) * facti
                   qclust(2,iclst) = qclust(2,iclst) * facti
                   qclust(3,iclst) = qclust(3,iclst) * facti
                   qclust(4,iclst) =
     &                  dsqrt( qclust(1,iclst)**2
     &                  + qclust(2,iclst)**2 + qclust(3,iclst)**2
     &                  + qclust(5,iclst)**2 )

                   qclust(7,iclst) = 1000d0
     &                  * ( qclust(4,iclst) - qclust(5,iclst) )
                 if ( qclust(7,iclst) .lt. 0d0 ) then
                    nclst = -1
                    return
                 end if

                 poutx = poutx + qclust(1,iclst)
                 pouty = pouty + qclust(2,iclst)
                 poutz = poutz + qclust(3,iclst)

                 if( jclust(7,iclst) .eq. 2112 .or.
     &                jclust(7,iclst) .eq. 2212 ) then
                    enucl = enucl
     &                   + qclust(4,iclst) - qclust(5,iclst)

                 else if( jclust(7,iclst) .eq. -211 .or.
     &                   jclust(7,iclst) .eq.  111 .or.
     &                   jclust(7,iclst) .eq.  211 ) then
                    epion = epion + qclust(4,iclst)

                 else if( jclust(7,iclst) .eq. -321 .or.
     &                   jclust(7,iclst) .eq.  311 .or.
     &                   jclust(7,iclst) .eq.  321 ) then
                    ekaon = ekaon + qclust(4,iclst)

                 else if( jclust(7,iclst) .eq. -13 .or.
     &                   jclust(7,iclst) .eq.  13 ) then
                    emuon = emuon + qclust(4,iclst)

                 else if( jclust(7,iclst) .eq. 22 ) then
                    egamm = egamm + qclust(4,iclst)

                 else
                  if( jclust(6,iclst)*3 .eq. 3 ) then
                     eelss = qclust(4,iclst) - rmnuc
                  else if( jclust(6,iclst) .eq. -3 ) then
                     eelss = qclust(4,iclst) + rmnuc
                  else
                     eelss = qclust(4,iclst)
                  end if
                  eelse = eelse + eelss

                 end if

                end do

                gein = eein / 1000.0
                piabs = dsqrt( gein**2 + 2.0 * gein * projm )
                presx = - poutx
                presy = - pouty
                presz = - poutz + piabs
                pabst = dsqrt( presx**2 + presy**2 + presz**2 )
                rsmas = 0.93827d0 * dble( mprot )
     &               + 0.93958d0 * dble( mneut )
                etota = dsqrt( pabst**2 + rsmas**2 )
                erres = etota - rsmas
                exres = gein + pima
     &               - enucl - epion - ekaon
     &               - emuon - egamm - eelse
     &               - erres
     &               - bindeg(mchg,mmas-mchg) / 1000.0

               end do

               if( dabs(exres) .le. epsiloni ) then
                  exres = 0d0

               else
                  nclst = -1
                  return

               end if

              end if


              if( exres .lt. 0.0d0 ) then

                 nclst = -1
                 return

              end if

              exres = max( 0.0d0, exres )

*-----------------------------------------------------------------------
*           booking of the residual nucleons (only neutrons or protons)
*-----------------------------------------------------------------------

              if ( mprot .eq. 0 ) then ! for only neutrons

               do iresidue = 1, nresidue

                  nclst = nclst + 1

                  iclust(nclst)    = 2

                  jclust(0,nclst)  = 0
                  jclust(1,nclst)  = 0
                  jclust(2,nclst)  = 1
                  jclust(3,nclst)  = 2
                  jclust(4,nclst)  = 0
                  jclust(5,nclst)  = 0
                  jclust(6,nclst)  = 1
                  jclust(7,nclst)  = 2112
                  jclust(8,nclst)  = 0

                  qclust(0,nclst)  = pard(2)
                  qclust(1,nclst)  = presx / dble( nresidue )
                  qclust(2,nclst)  = presy / dble( nresidue )
                  qclust(3,nclst)  = presz / dble( nresidue )
                  qclust(4,nclst)  = etota / dble( nresidue )
                  qclust(5,nclst)  = 0.93958d0
                  qclust(6,nclst)  = 0.0d0
                  fkinene = (qclust(4,nclst) - qclust(5,nclst)) *1000d0
                  qclust(7,nclst)  = fkinene
                  qclust(8,nclst)  = 1.0d0
                  qclust(9,nclst)  = 0.0d0
                  qclust(10,nclst) = 0.0d0
                  qclust(11,nclst) = 0.0d0
                  qclust(12,nclst) = 0.0d0

                  if( fkinene .lt. 0.0d0 ) then
                     nclst = -1
                     return
                  end if

               end do

                  numpat(2)  = nneut + mneut
                  rumpat(2)  = nneut + mneut

*-----------------------------------------------------------------------

              else if  ( mneut .eq. 0 ) then ! for only protons

               do iresidue = 1, nresidue

                  nclst = nclst + 1

                  iclust(nclst)    = 1

                  jclust(0,nclst)  = 0
                  jclust(1,nclst)  = 1
                  jclust(2,nclst)  = 0
                  jclust(3,nclst)  = 1
                  jclust(4,nclst)  = 0
                  jclust(5,nclst)  = 1
                  jclust(6,nclst)  = 1
                  jclust(7,nclst)  = 2212
                  jclust(8,nclst)  = 0

                  qclust(0,nclst)  = pard(2)
                  qclust(1,nclst)  = presx / dble( nresidue )
                  qclust(2,nclst)  = presy / dble( nresidue )
                  qclust(3,nclst)  = presz / dble( nresidue )
                  qclust(4,nclst)  = etota / dble( nresidue )
                  qclust(5,nclst)  = 0.93827d0
                  qclust(6,nclst)  = 0.0d0
                  fkinene = (qclust(4,nclst) - qclust(5,nclst)) *1000d0
                  qclust(7,nclst)  = fkinene
                  qclust(8,nclst)  = 1.0d0
                  qclust(9,nclst)  = 0.0d0
                  qclust(10,nclst) = 0.0d0
                  qclust(11,nclst) = 0.0d0
                  qclust(12,nclst) = 0.0d0

                  if( fkinene .lt. 0.0d0 ) then
                     nclst = -1
                     return
                  end if

               end do

                  numpat(1)  = nprot + mprot
                  rumpat(1)  = nprot + mprot

*-----------------------------------------------------------------------

              else

                 nclst = -1
                 return

              end if

*-----------------------------------------------------------------------

             else ! except for only neutrons or only protons

                  nclst = nclst + 1

                  gein = eein / 1000.0

                  piabs = sqrt( gein**2 + 2.0 * gein * projm )

                  presx = - poutx
                  presy = - pouty
                  presz = - poutz + piabs

                  pabst = sqrt( presx**2 + presy**2 + presz**2 )
                  rsmas = rmnuc * dble( mnucl )

                  etota = sqrt( pabst**2 + rsmas**2 )

                  erres = etota - rsmas

                  exres = gein + pima
     &                  - enucl - epion - ekaon
     &                  - emuon - egamm - eelse
     &                  - erres
     &                  - bindeg(mchg,mmas-mchg) / 1000.0
     &                  + bindeg(mprot,mneut) / 1000.0


C S.H. iteration to satisfy energy- and momentum-conservation
C      when a residue is a nucleon. (2014.8.13)
                  if ( exres + gein .lt. 0d0 ) then
                     nclst = -1
                     return
                  end if

                  absexres = dabs(exres)
                  if ( absexres .gt. 0d0 .and.
     &                 mprot + mneut .eq. 1 ) then
                     niteration = 10000
                     epsiloni = 1d-5 !(10 keV)
                     if ( exres .lt. 0d0 ) then
                        facti = 1d0 - 1d-5
                     else
                        facti = 1d0 + 1d-5
                     end if

                     iteration = 0
                     exresini = exres
                     do while ( exresini*exres .gt. 0d0 .and.
     &                    iteration .le. niteration )
                     iteration = iteration + 1

                     poutx = 0d0
                     pouty = 0d0
                     poutz = 0d0
                     enucl = 0d0
                     epion = 0d0
                     ekaon = 0d0
                     emuon = 0d0
                     egamm = 0d0
                     eelse = 0d0

                     do iclst=1,nclst-1
                        qclust(1,iclst) = qclust(1,iclst) * facti
                        qclust(2,iclst) = qclust(2,iclst) * facti
                        qclust(3,iclst) = qclust(3,iclst) * facti

                        qclust(4,iclst) =
     &                       dsqrt( qclust(1,iclst)**2
     &                       + qclust(2,iclst)**2 + qclust(3,iclst)**2
     &                       + qclust(5,iclst)**2 )

                        qclust(7,iclst) = 1000d0
     &                       * ( qclust(4,iclst) - qclust(5,iclst) )
                        if ( qclust(7,iclst) .lt. 0d0 ) then
                           nclst = -1
                           return
                        end if

                        poutx = poutx + qclust(1,iclst)
                        pouty = pouty + qclust(2,iclst)
                        poutz = poutz + qclust(3,iclst)

                        if( jclust(7,iclst) .eq. 2112 .or.
     &                       jclust(7,iclst) .eq. 2212 ) then
                           enucl = enucl
     &                          + qclust(4,iclst) - qclust(5,iclst)

                        else if( jclust(7,iclst) .eq. -211 .or.
     &                          jclust(7,iclst) .eq.  111 .or.
     &                          jclust(7,iclst) .eq.  211 ) then
                           epion = epion + qclust(4,iclst)

                        else if( jclust(7,iclst) .eq. -321 .or.
     &                          jclust(7,iclst) .eq.  311 .or.
     &                          jclust(7,iclst) .eq.  321 ) then
                           ekaon = ekaon + qclust(4,iclst)

                        else if( jclust(7,iclst) .eq. -13 .or.
     &                          jclust(7,iclst) .eq.  13 ) then
                           emuon = emuon + qclust(4,iclst)

                        else if( jclust(7,iclst) .eq. 22 ) then
                           egamm = egamm + qclust(4,iclst)

                        else
                           if( jclust(6,iclst)*3 .eq. 3 ) then
                              eelss = qclust(4,iclst) - rmnuc
                           else if( jclust(6,iclst) .eq. -3 ) then
                              eelss = qclust(4,iclst) + rmnuc
                           else
                              eelss = qclust(4,iclst)
                           end if
                           eelse = eelse + eelss

                        end if

                     end do

                     gein = eein / 1000.0
                     piabs = sqrt( gein**2 + 2.0 * gein * projm )
                     presx = - poutx
                     presy = - pouty
                     presz = - poutz + piabs
                     pabst = sqrt( presx**2 + presy**2 + presz**2 )
                     rsmas = 0.93827d0 * dble( mprot )
     &                    + 0.93958d0 * dble( mneut )
                     etota = sqrt( pabst**2 + rsmas**2 )
                     erres = etota - rsmas
                     exres = gein + pima
     &                    - enucl - epion - ekaon
     &                    - emuon - egamm - eelse
     &                    - erres
     &                    - bindeg(mchg,mmas-mchg) / 1000.0
     &                    + bindeg(mprot,mneut) / 1000.0

                  end do

                     if( dabs(exres) .le. epsiloni ) then
                        exres = 0d0

                     else
                        nclst = -1
                        return

                     end if

                  end if

                  if( exres .lt. 0.0d0 ) then

                     nclst = -1
                     return

                  end if

                  exres = max( 0.0d0, exres )

*-----------------------------------------------------------------------
*           booking of the residual nucleus
*-----------------------------------------------------------------------

                  numpat(0) = 1
                  rumpat(0) = 1

                  iclust(nclst)    = 0

                  jclust(0,nclst)  = 0
                  jclust(1,nclst)  = mprot
                  jclust(2,nclst)  = mneut
                  jclust(3,nclst)  = 19
                  jclust(4,nclst)  = 0
                  jclust(5,nclst)  = mprot
                  jclust(6,nclst)  = mprot + mneut
                  jclust(7,nclst)  = mprot * 1000000 + mprot + mneut
                  jclust(8,nclst)  = 0

                  qclust(0,nclst)  = pard(2)
                  qclust(1,nclst)  = presx
                  qclust(2,nclst)  = presy
                  qclust(3,nclst)  = presz
                  qclust(4,nclst)  = etota
                  qclust(5,nclst)  = rmnuc * mnucl
                  qclust(6,nclst)  = exres * 1000.0
                  qclust(7,nclst)  = ( etota - rmnuc * mnucl ) * 1000.
                  qclust(8,nclst)  = 1.0
                  qclust(9,nclst)  = 0.0
                  qclust(10,nclst) = 0.0d0
                  qclust(11,nclst) = 0.0d0
                  qclust(12,nclst) = 0.0d0

*-----------------------------------------------------------------------

             end if

*-----------------------------------------------------------------------

            else

                  numpat(0) = 0
                  rumpat(0) = 0

C S.H. iteration to satisfy energy- and momentum-conservation
C      when a residue is none. (2014.8.13)

                  gein = eein / 1000.0
                  piabs = sqrt( gein**2 + 2.0 * gein * projm )
                  presx = - poutx
                  presy = - pouty
                  presz = - poutz + piabs
                  pabst = sqrt( presx**2 + presy**2 + presz**2 )
                  exres = gein + pima
     &                 - enucl - epion - ekaon
     &                 - emuon - egamm - eelse
     &                 - pabst
     &                 - bindeg(mchg,mmas-mchg) / 1000.0

                  if ( exres + gein .lt. 0d0 ) then
                     nclst = -1
                     return
                  end if

                  absexres = dabs(exres)
                  if ( absexres .gt. 0d0 ) then
                     niteration = 10000
                     epsiloni = 1d-5 !(10 keV)
                     if ( exres .lt. 0d0 ) then
                        facti = 1d0 - 1d-5
                     else
                        facti = 1d0 + 1d-5
                     end if

                     iteration = 0
                     exresini = exres
                     do while ( exresini*exres .gt. 0d0 .and.
     &                    iteration .le. niteration )
                     iteration = iteration + 1

                     poutx = 0d0
                     pouty = 0d0
                     poutz = 0d0
                     enucl = 0d0
                     epion = 0d0
                     ekaon = 0d0
                     emuon = 0d0
                     egamm = 0d0
                     eelse = 0d0

                     do iclst=1,nclst
                        if ( iclst .ne. nclst ) then
                           qclust(1,iclst) = qclust(1,iclst) * facti
                           qclust(2,iclst) = qclust(2,iclst) * facti
                           qclust(3,iclst) = qclust(3,iclst) * facti
                           poutx = poutx + qclust(1,iclst)
                           pouty = pouty + qclust(2,iclst)
                           poutz = poutz + qclust(3,iclst)

                        else
                           qclust(1,iclst) = - poutx
                           qclust(2,iclst) = - pouty
                           qclust(3,iclst) = - poutz + piabs

                        end if

                        qclust(4,iclst) =
     &                       dsqrt( qclust(1,iclst)**2
     &                       + qclust(2,iclst)**2 + qclust(3,iclst)**2
     &                       + qclust(5,iclst)**2 )

                        qclust(7,iclst) = 1000d0
     &                       * ( qclust(4,iclst) - qclust(5,iclst) )
                        if ( qclust(7,iclst) .lt. 0d0 ) then
                           nclst = -1
                           return
                        end if

                        if( jclust(7,iclst) .eq. 2112 .or.
     &                       jclust(7,iclst) .eq. 2212 ) then
                           enucl = enucl
     &                          + qclust(4,iclst) - qclust(5,iclst)

                        else if( jclust(7,iclst) .eq. -211 .or.
     &                          jclust(7,iclst) .eq.  111 .or.
     &                          jclust(7,iclst) .eq.  211 ) then
                           epion = epion + qclust(4,iclst)

                        else if( jclust(7,iclst) .eq. -321 .or.
     &                          jclust(7,iclst) .eq.  311 .or.
     &                          jclust(7,iclst) .eq.  321 ) then
                           ekaon = ekaon + qclust(4,iclst)

                        else if( jclust(7,iclst) .eq. -13 .or.
     &                          jclust(7,iclst) .eq.  13 ) then
                           emuon = emuon + qclust(4,iclst)

                        else if( jclust(7,iclst) .eq. 22 ) then
                           egamm = egamm + qclust(4,iclst)

                        else
                           if( jclust(6,iclst)*3 .eq. 3 ) then
                              eelss = qclust(4,iclst) - rmnuc
                           else if( jclust(6,iclst) .eq. -3 ) then
                              eelss = qclust(4,iclst) + rmnuc
                           else
                              eelss = qclust(4,iclst)
                           end if
                           eelse = eelse + eelss

                        end if

                     end do

                     exres = gein + pima
     &                    - enucl - epion - ekaon
     &                    - emuon - egamm - eelse
     &                    - bindeg(mchg,mmas-mchg) / 1000.0

                  end do

                  if( dabs(exres) .le. epsiloni ) then
                     exres = 0d0

                  else
                     nclst = -1
                     return

                  end if

                  end if

            end if

*-----------------------------------------------------------------------
*           x-y rotation by randum
*-----------------------------------------------------------------------

                  theta = 2.0d0 * pi * rn(0)

               do i = 1, nclst

                  px = qclust(1,i)
                  py = qclust(2,i)

                  qclust(1,i)  = px * cos( theta ) - py * sin( theta )
                  qclust(2,i)  = px * sin( theta ) + py * cos( theta )

               end do

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sigjam(kf1,kf2,ein,sig,sigel,signo)
*                                                                      *
*        calculates total, nonelastic and elastic cross-sections       *
*        of two particles by JAM code                                  *
*                                                                      *
*     input:                                                           *
*        kf1    : first particle kf code                               *
*        kf2    : second particle kf code                              *
*        ein    : incident nucleon energy (MeV)                        *
*                                                                      *
*     output:                                                          *
*       sig     : total cross-section (b)                              *
*       sigel   : sigt-sigr=elastic scattering cross-section (b)       *
*       signo   : nonelastic cross-section (b)                         *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'jam2.inc'

      parameter(mxchan=30)
      dimension sigin(mxchan)

*-----------------------------------------------------------------------

      pawt(a,b,c)=sqrt((a**2-(b+c)**2)*(a**2-(b-c)**2))/(2.d0*a)

*-----------------------------------------------------------------------
cFURUTA  JAM initialization in case called before JAM reaction
            if(mstc(21).eq.0)then
              call jamsetpa
            endif

            sig   = 0.0
            sigel = 0.0
            signo = 0.0

            if( kf1 .eq. 0 .or. kf2 .eq. 0 ) return

            if( abs( kf1 ) .lt. 100 .or.
     &          abs( kf2 ) .lt. 100 ) return

*-----------------------------------------------------------------------

            elab = ein / 1000.0

            if( elab .lt. 0.002 ) elab = 0.002

*-----------------------------------------------------------------------

               em1=pymass(kf1)
               em2=pymass(kf2)

               ibar1=kchg(jamcomp(kf1),6)
               ibar2=kchg(jamcomp(kf2),6)
               icltyp=jamcltyp(kf1,kf2,ibar1,ibar2)
               pproj=sqrt(elab*(2.0*em1+elab))
               srt=sqrt((elab+em1+em2)**2-pproj**2)
               pr=pawt(srt,em1,em2)

            call jamcross(1,icltyp,srt,pr,kf1,kf2,em1,em2,
     &                    sig,sigel,sigin,mchanel,mabsrb,ijet,icon)

*-----------------------------------------------------------------------

               sig   = sig   / 1000.0
               sigel = sigel / 1000.0
               signo = sig - sigel

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine jamqmdgd(ierr)
*                                                                      *
*        Last Revised:     2004/09/21                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              make the goround state of target and projectile         *
*              for JAMQMD mode using the random packing method         *
*                                                                      *
*                                                                      *
************************************************************************

      include 'jam1.inc'
      include 'jam2.inc'

*-----------------------------------------------------------------------

      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)
!$OMP THREADPRIVATE(/coodid/)

      common /const1/ elab, rdist, bmin, bmax, ibch, ibin
!$OMP THREADPRIVATE(/const1/)
      common /impact/ bval(1000), ibnum(1000), bweight(1000), bdef
!$OMP THREADPRIVATE(/impact/)
      common /vriab1/ b, llnow, ntnow
!$OMP THREADPRIVATE(/vriab1/)

      data ibsf /0/
      save ibsf
!$OMP THREADPRIVATE(ibsf)
*-----------------------------------------------------------------------

      common /jamnmtc/ kfjam
!$OMP THREADPRIVATE(/jamnmtc/)

      common /rqmdaux/ ibsv
!$OMP THREADPRIVATE(/rqmdaux/)
      common /qmdsuc/ isucc
!$OMP THREADPRIVATE(/qmdsuc/)

*-----------------------------------------------------------------------
*     Initialize QMD
*-----------------------------------------------------------------------

               call qmdint

*-----------------------------------------------------------------------
*        Set total particle number.
*-----------------------------------------------------------------------

            nv = mstd(2) + mstd(5)
            nmeson = 0
            nbary  = nv

*-----------------------------------------------------------------------
*        In the case of proj. mass number 1.
*-----------------------------------------------------------------------

            if(mstd(2).eq.1) then
                kc1=jamcomp(mstd(1))
                if(kchg(kc1,6).eq.0) then
                  nmeson = 1
                  nbary  = nv - 1
                endif
            endif

*-----------------------------------------------------------------------
*        In the case of targ. mass number 1.
*-----------------------------------------------------------------------

            if(mstd(5).eq.1) then
              kc2=jamcomp(mstd(4))
              if(kchg(kc2,6).eq.0) then
                nmeson = nmeson + 1
                nbary  = nbary - 1
              endif
            endif

*-----------------------------------------------------------------------
*       QMD;  impact parameter and multi run control
*-----------------------------------------------------------------------
               if( ibsv .eq. 1 ) goto 10  ! skip b sampling and use previous b

               if( ibch .eq. 0 ) then

                     b = sqrt( max( 0.0d0,
     &                   bmin**2 + ( bmax**2 - bmin**2 ) * rn(0) ) )

               else if( ibch .eq. 1 ) then

                     ibsf = ibsf + 1

                  if( ibsf .gt. ibin ) then

                     ibsf = 1

                  end if

                     b = bval(ibsf) - bdef / 2.0 + bdef * rn(0)

               end if


*-----------------------------------------------------------------------
*     make ground stae  and  boost by QMD
*-----------------------------------------------------------------------

  10        do isucc = 0, 1
             call ground_preR(ierr)
            enddo
            ibsv = 0    ! skip
            isucc = 2
            call ground_finalR(ierr)

            if( ierr .ne. 0 ) return

*-----------------------------------------------------------------------
*    make index for JAM
*    Loop over proj. and targ. in=1:targ. in=2:proj.
*-----------------------------------------------------------------------

      do 1000 in = 1, 2

*-----------------------------------------------------------------------
*     Target loop.
*-----------------------------------------------------------------------

       if(in.eq.1) then

         nnm = mstd(5)
         nnp = mstd(6)

         iofset = 0

*-----------------------------------------------------------------------
*      Projectile loop.
*-----------------------------------------------------------------------

       else if(in.eq.2) then

         nnm = mstd(2)
         nnp = mstd(3)

         iofset = mstd(5)

       endif

*-----------------------------------------------------------------------
*     Loop over all particles.
*-----------------------------------------------------------------------

      do ie = 1, nnm

        i = ie + iofset

      if(ie.le.nnp) then

         if(nnm.eq.1) then

           if(in.eq.2) then
             kf=mstd(1)
             kc=jamcomp(kf)
           else if(in.eq.1) then
             kf=mstd(4)
             kc=jamcomp(kf)
           endif

             k(2,i)=kf
             k(3,i)=0
             k(4,i)=0
             k(8,i)=1
             k(9,i)=kchg(kc,6)*isign(1,kf)

         else

             k(2,i)=2212
             k(3,i)=0
             k(4,i)=0
             k(8,i)=1
             k(9,i)=3

         endif

      else

             k(2,i)=2112
             k(3,i)=0
             k(4,i)=0
             k(8,i)=1
             k(9,i)=3

      endif

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
             v(5,i) = 1.d+35

         do j=1,10
             vq(j,i)=0.d0
         end do

             kq(1,i)=0
             kq(2,i)=0

      end do

1000  continue

*-----------------------------------------------------------------------
*     Boost by QMD
*-----------------------------------------------------------------------

            call rboostR(ierr)

               if( ierr .ne. 0 ) return

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine qmdtrspt(tv0)
*                                                                      *
*        Last Revised:     2004/09/23                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              transport particle by QMD mode                          *
*                                                                      *
*                                                                      *
************************************************************************
      use QMD_COOD2_MOD, only : d1r, d1p, f0r, f0p

      include 'jam1.inc'
      include 'jam2.inc'

*-----------------------------------------------------------------------

      common /vriab0/ massal, massba, mmeson
!$OMP THREADPRIVATE(/vriab0/)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)
!$OMP THREADPRIVATE(/coodid/)

      common /jamqmd/ dtqjam
!$OMP THREADPRIVATE(/jamqmd/)

*-----------------------------------------------------------------------

      parameter ( c0 = 1.0 )

      parameter ( dt2 = c0 )
      parameter ( dt3 = 1.0 / 2.0 / c0 )
      parameter ( dt1 = 1.0 - c0 - dt3 )

      d1r = 0.d0
      d1p = 0.d0
      f0r = 0.d0
      f0p = 0.d0

*-----------------------------------------------------------------------
*     Time mesh, max mesh is dtqjam fm/c
*-----------------------------------------------------------------------

               tmax = 0.0d0
               tmin = 1.d+10

            do 10 i = 1, massal
               if( inds(i) .eq. 0 ) goto 10

               timi = max( 0.0d0, tv0 - r(4,i) )

               if( timi .gt. tmax ) then
                  tmax = timi
                  tmin = r(4,i)
               end if

   10       continue

            if( tmax .eq. 0.0d0 ) return

            if( tmax .gt. dtqjam ) then

               nt = int( tmax ) + 1
               dv = tmax / dble( nt )

            else

               nt = 1
               dv = tmax

            end if

*-----------------------------------------------------------------------

            call jam2qmd

      do ii = 1, nt

            tv = tmin + ii * dv

            call caldisaR

*-----------------------------------------------------------------------
*     RKG12 for JAMQMD
*-----------------------------------------------------------------------

               call graduR(d1r,d1p)

         do 20 i = 1, massal
               if( inds(i) .eq. 0 ) goto 20

               dt = tv - r(4,i)
               if( dt .le. 0.0d0 ) goto 20

            do j = 1, 3

               r(j,i) = r(j,i) + dt * dt3 * d1r(j,i)
               p(j,i) = p(j,i) + dt * dt3 * d1p(j,i)

               f0r(j,i) = d1r(j,i)
               f0p(j,i) = d1p(j,i)

            end do

               p(6,i) = sqrt( p(5,i)**2 + p(1,i)**2
     &                      + p(2,i)**2 + p(3,i)**2 )

   20    continue

               call caldisaR

         do i = 1, massal
               if( inds(i) .eq. 0 ) cycle

               call epotprtR(i,epotp)

               p(4,i) = sqrt( p(1,i)**2
     &                      + p(2,i)**2
     &                      + p(3,i)**2
     &                      + 2 * p(5,i) * epotp
     &                      + p(5,i)**2 )

         end do

*-----------------------------------------------------------------------

               call graduR(d1r,d1p)

         do 30 i = 1, massal
               if( inds(i) .eq. 0 ) goto 30

               dt = tv - r(4,i)
               if( dt .le. 0.0d0 ) goto 30

            do j = 1, 3

               r(j,i) = r(j,i)
     &                + dt * ( f0r(j,i) * dt1 + d1r(j,i) * dt2 )
               p(j,i) = p(j,i)
     &                + dt * ( f0p(j,i) * dt1 + d1p(j,i) * dt2 )

            end do

               p(6,i) = sqrt( p(5,i)**2 + p(1,i)**2
     &                      + p(2,i)**2 + p(3,i)**2 )

               r(4,i) = tv

               call jamtupda(i)

   30    continue

                call caldisaR

         do i = 1, massal
               if( inds(i) .eq. 0 ) cycle

               call epotprtR(i,epotp)

               p(4,i) = sqrt( p(1,i)**2
     &                      + p(2,i)**2
     &                      + p(3,i)**2
     &                      + 2 * p(5,i) * epotp
     &                      + p(5,i)**2 )

         end do

*-----------------------------------------------------------------------

      end do

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine jam2qmd
*                                                                      *
*        Last Revised:     2004/09/23                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              transform all inf. from JAM to QMD                      *
*                                                                      *
*                                                                      *
************************************************************************

      include 'jam1.inc'
      include 'jam2.inc'

*-----------------------------------------------------------------------

      common /vriab0/ massal, massba, mmeson
!$OMP THREADPRIVATE(/vriab0/)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)
!$OMP THREADPRIVATE(/coodid/)

*-----------------------------------------------------------------------
*     from JAM to QMD
*-----------------------------------------------------------------------

               massal = nv
               massba = nbary
               mmeson = nmeson

         do i = 1, nv

               kf = k(2,i)
               kc = jamcomp(kf)

               ibry(i) = 0
               if( abs(k(9,i)) .ge. 3 ) ibry(i) = 1

               ichg(i) = 0
               if( kc .gt. 0 )
     &         ichg(i) = kchg(kc,1) / 3 * isign(1,kf)

               inuc(i) = 0

            if( k(1,i) .ge. 11 ) then

               inds(i) = 0

            else if( kf .eq. 2212 .or. kf .eq. 2112 ) then

               inuc(i) = 1
               inds(i) = 1

            else if( kf .eq. 2224 .or. kf .eq. 2214 .or.
     &               kf .eq. 2114 .or. kf .eq. 1114 ) then

               inds(i) = 2

            else if( kf .eq. 12212 .or. kf .eq. 12112 ) then

               inds(i) = 3

            else if( kf .eq.  211 .or. kf .eq. 111 .or.
     &               kf .eq. -211 ) then

               inds(i) = 4

            else

               inds(i) = 5

            end if

         end do

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine jam2qmd2(i1,i2)
*                                                                      *
*        Last Revised:     2004/09/23                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              transform i1 and i2 inf. from JAM to QMD                *
*                                                                      *
*                                                                      *
************************************************************************

      include 'jam1.inc'
      include 'jam2.inc'

*-----------------------------------------------------------------------

      common /vriab0/ massal, massba, mmeson
!$OMP THREADPRIVATE(/vriab0/)
      common /coodid/ ichg(nnn), inuc(nnn), ibry(nnn), inds(nnn),
     &                inun(nnn), iavd(nnn), ihis(nnn)
!$OMP THREADPRIVATE(/coodid/)

*-----------------------------------------------------------------------
*     from JAM to QMD
*-----------------------------------------------------------------------

               massal = nv

               massba = nbary
               mmeson = nmeson

               idd = i2 - i1
               if( i1 .eq. i2 ) idd = 1

         do i = i1, i2, idd

               kf = k(2,i)
               kc = jamcomp(kf)

               ibry(i) = 0
               if( abs(k(9,i)) .ge. 3 ) ibry(i) = 1

               ichg(i) = 0
               if( kc .gt. 0 )
     &         ichg(i) = kchg(kc,1) / 3 * isign(1,kf)

               inuc(i) = 0

            if( k(1,i) .ge. 11 ) then

               inds(i) = 0

            else if( kf .eq. 2212 .or. kf .eq. 2112 ) then

               inuc(i) = 1
               inds(i) = 1

            else if( kf .eq. 2224 .or. kf .eq. 2214 .or.
     &               kf .eq. 2114 .or. kf .eq. 1114 ) then

               inds(i) = 2

            else if( kf .eq. 12212 .or. kf .eq. 12112 ) then

               inds(i) = 3

            else if( kf .eq.  211 .or. kf .eq. 111 .or.
     &               kf .eq. -211 ) then

               inds(i) = 4

            else

               inds(i) = 5

            end if

         end do

*-----------------------------------------------------------------------

      return
      end

